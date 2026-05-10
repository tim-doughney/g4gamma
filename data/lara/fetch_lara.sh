#!/usr/bin/env bash
# Fetch LARA decay-data files from LNHB.
#
# Usage:
#   ./fetch_lara.sh                       # fetch the default NORM-relevant set
#   ./fetch_lara.sh iso1 iso2 ...         # fetch specific isotopes (e.g. Eu-152)
#   ./fetch_lara.sh --all                 # brute-force-discover everything
#                                         # available on lnhb.fr/nuclides/
#                                         # (~5-10 min depending on network)
#   ./fetch_lara.sh --pack                # pack downloaded files into lara.tar.gz
#                                         # and remove the loose .lara.txt files
#   ./fetch_lara.sh --all --pack          # do both
#
# Env vars:
#   FORCE=1     re-fetch even if files exist locally
#   PARALLEL=N  parallel curl jobs for --all (default 8)

set -e
cd "$(dirname "$0")"

MODE="default"
PACK=0
ARGS=()
for a in "$@"; do
    case "$a" in
        --all)  MODE="all" ;;
        --pack) PACK=1 ;;
        --help|-h)
            sed -n '2,18p' "$0"; exit 0 ;;
        *)      ARGS+=("$a") ;;
    esac
done
[ ${#ARGS[@]} -gt 0 ] && MODE="explicit"

DEFAULT_ISOS=(
    Cs-137 K-40 Co-60 Na-22
    Ra-226 Pb-214 Bi-214 Pb-210 Po-210 Po-218
    Th-232 Th-228 Ra-228 Ac-228 Ra-224 Rn-220 Po-216 Pb-212 Bi-212 Tl-208 Po-212
    U-238 U-234 U-235 Th-230 Th-234 Pa-234m Pa-234
    Am-241 Eu-152 Ba-133 Mn-54 Fe-59 Zn-65 Y-88
    Ce-139 Hg-203 Sr-85 Co-57 Cd-109 Sn-113 Cs-134
    Bi-207 Eu-154 Eu-155
)

# Periodic table for brute-force discovery (atomic symbols). Matches the
# isotope filename convention "<Sym>-<A>[m].lara.txt".
ELEMENTS=(
    H He Li Be B C N O F Ne Na Mg Al Si P S Cl Ar K Ca
    Sc Ti V Cr Mn Fe Co Ni Cu Zn Ga Ge As Se Br Kr Rb Sr Y Zr
    Nb Mo Tc Ru Rh Pd Ag Cd In Sn Sb Te I Xe Cs Ba La Ce Pr Nd
    Pm Sm Eu Gd Tb Dy Ho Er Tm Yb Lu Hf Ta W Re Os Ir Pt Au Hg
    Tl Pb Bi Po At Rn Fr Ra Ac Th Pa U Np Pu Am Cm Bk Cf Es Fm
)

# Reasonable mass ranges per Z. LNHB's ~400 nuclides cluster within ~10
# units of the line of stability. The most-stable A for atomic number Z is
# approximated by A ≈ Z*(1 + Z/(2*c)) with c ≈ 200 (Weizsäcker form).
# This gives:  Z=1 → A≈1, Z=20 → A≈42, Z=50 → A≈112, Z=82 → A≈220, Z=92 → A≈250
# Close enough for setting probe ranges.
build_brute_list() {
    local i=0
    for elem in "${ELEMENTS[@]}"; do
        i=$((i+1))
        # A_stable ≈ 2Z + Z²/180 (matches measured β-stable line within ±2)
        local astable=$(( 2*i + (i*i)/180 ))
        local amin=$(( astable - 10 )); [ $amin -lt 1 ] && amin=1
        local amax=$(( astable + 14 ))
        for a in $(seq $amin $amax); do
            echo "${elem}-${a}"
            echo "${elem}-${a}m"
        done
    done
}

case "$MODE" in
    explicit) ISOS=("${ARGS[@]}") ;;
    all)
        echo "Building brute-force isotope list..."
        readarray -t ISOS < <(build_brute_list)
        echo "Will probe ${#ISOS[@]} candidate filenames against lnhb.fr/nuclides/"
        ;;
    *)        ISOS=("${DEFAULT_ISOS[@]}") ;;
esac

# Detect curl / wget. Use whichever is available.
FETCHER=""
if command -v curl >/dev/null 2>&1; then
    FETCHER="curl"
elif command -v wget >/dev/null 2>&1; then
    FETCHER="wget"
else
    echo "ERROR: need curl or wget. Install one and retry." >&2
    exit 1
fi

fetch_one() {
    local iso="$1"
    local out="${iso}.lara.txt"
    if [ "${FORCE:-0}" != "1" ] && [ -f "$out" ] && [ "$(wc -c < "$out")" -gt 200 ]; then
        echo "SKIP"; return
    fi

    # Primary URL: /nuclides/<Sym>-<A>.lara.txt  (works for ~95% of nuclides)
    _fetch_url "http://www.lnhb.fr/nuclides/${iso}.lara.txt" "$out.tmp"
    if [ $? -eq 0 ]; then echo "OK"; return; fi

    # Fallback: LaraWEB Results directory uses versioned names like Th-230_@03.lara.txt.
    # Query Result_Lara2.php to find the actual filename, then download it.
    local lara_name
    lara_name=$(curl -sS --max-time 12 -A "Mozilla/5.0" -X POST \
        "http://www.lnhb.fr/Laraweb/Result_Lara2.php" \
        -d "Nuclide[]=${iso//-/}&Show_data=Show+all+data&Form_Origine=1&Data_show=1&Emi_show=1&Langue=ENG" \
        2>/dev/null | grep -o "Results/${iso}_[^\"]*\.lara\.txt" | head -1)
    if [ -n "$lara_name" ]; then
        _fetch_url "http://www.lnhb.fr/Laraweb/${lara_name}" "$out.tmp"
        if [ $? -eq 0 ]; then echo "OK"; return; fi
    fi

    rm -f "$out.tmp"
    echo "FAIL"
}

_fetch_url() {
    local url="$1" dest="$2"
    if [ "$FETCHER" = "curl" ]; then
        curl -fLsS --max-time 8 -A "Mozilla/5.0" -o "$dest" "$url" 2>/dev/null
    else
        wget -q --timeout=8 -O "$dest" "$url"
    fi
    if [ $? -ne 0 ]; then rm -f "$dest"; return 1; fi
    local sz; sz=$(wc -c < "$dest" 2>/dev/null || echo 0)
    if [ "$sz" -gt 200 ] && ! grep -q "Page non trouvée\|404" "$dest" 2>/dev/null; then
        mv "$dest" "${dest%.tmp}"
        return 0
    fi
    rm -f "$dest"
    return 1
}

# Export for parallel sub-shells if needed.
export FETCHER FORCE
export -f fetch_one _fetch_url

OK=0
SKIP=0
FAIL=0
total=${#ISOS[@]}
echo "Fetching $total isotopes..."

if [ "$MODE" = "all" ] && command -v xargs >/dev/null 2>&1; then
    # Parallel mode for brute-force --all -- much faster.
    P="${PARALLEL:-8}"
    echo "  (using $P parallel jobs; each '.' = success, 'x' = 404, 's' = skip)"
    # Print results as they complete.
    printf '%s\n' "${ISOS[@]}" | \
        xargs -n1 -P "$P" -I {} bash -c '
            r=$(fetch_one "$1")
            case "$r" in
                OK)   printf "." ;;
                SKIP) printf "s" ;;
                *)    printf "x" ;;
            esac
        ' _ {}
    echo
    # Recount for the summary
    OK=$(ls -1 *.lara.txt 2>/dev/null | wc -l)
    echo "Total .lara.txt files now in directory: $OK"
else
    # Sequential -- for explicit lists or default mode.
    for iso in "${ISOS[@]}"; do
        r=$(fetch_one "$iso")
        case "$r" in
            OK)   OK=$((OK+1));   printf "." ;;
            SKIP) SKIP=$((SKIP+1));printf "s" ;;
            *)    FAIL=$((FAIL+1));printf "x" ;;
        esac
    done
    echo
    echo "Results: $OK fetched, $SKIP already-present, $FAIL failed"
fi


if [ $PACK -eq 1 ]; then
    echo
    echo "Packing lara.tar.gz..."
    files=( *.lara.txt )
    if [ ${#files[@]} -eq 0 ] || [ ! -e "${files[0]}" ]; then
        echo "  no .lara.txt files to pack"
    else
        tar --format=ustar -czf lara.tar.gz "${files[@]}"
        sz=$(wc -c < lara.tar.gz)
        echo "  packed ${#files[@]} files into lara.tar.gz (${sz} bytes)"
        rm -f "${files[@]}"
        echo "  removed loose .lara.txt files (provider reads from tarball directly)"
    fi
fi
