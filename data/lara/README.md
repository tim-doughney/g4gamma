# LARA / DDEP nuclear data

This directory holds per-nuclide ASCII files from the
[NUCLÉIDE-LARA](http://www.lnhb.fr/Laraweb/) database, maintained by
LNE-LNHB / DDEP (Decay Data Evaluation Project).

LARA is the gold-standard reference for radionuclide metrology — the values
are the recommended evaluations published by an international working group
of decay-data evaluators. About 220 isotopes are evaluated to DDEP standard;
the remaining ~180 nuclides on the LARAWEB list use ENSDF.

## Storage formats

The `LaraProvider` accepts either layout, auto-detected:

1. **Tarball** — `lara.tar.gz` or `lara.tar` containing all the
   `<symbol>.lara.txt` files (preferred for the full ~400 nuclide set;
   keeps disk footprint to a few hundred KB).
2. **Loose directory** — individual `<symbol>.lara.txt` files in this
   directory.

Resolution order: explicit constructor argument → `$LARA_DATA_DIR` env →
`<repo>/data/lara/lara.tar.gz` → `<repo>/data/lara/lara.tar` →
`<repo>/data/lara/` (loose files).

## Quick start

```bash
# Default NORM-relevant + calibration set (~40 isotopes)
./fetch_lara.sh

# Full DDEP-evaluated set (~220 isotopes, ~5 minutes)
./fetch_lara.sh --all

# Pack into a tarball after fetching (recommended for the full set):
./fetch_lara.sh --all --pack
# This builds lara.tar.gz and removes the loose .lara.txt files.

# Specific isotopes
./fetch_lara.sh Eu-152 Cs-134 Bi-207

# Force re-download of existing files
FORCE=1 ./fetch_lara.sh Cs-137
```

The script uses `curl -L` (or `wget`) so HTTP→HTTPS redirects from LNHB
are handled correctly.

## File format

```
Nuclide ; Cs-137
Daughter(s) ; (B-) ; Ba-137 ; 100
Half-life (s) ; 947.3E6 ; 0.7E6
...
Energy (keV) ; Ener. unc. (keV) ; Intensity (%) ; Int. unc. (%) ; Type ; Origin ; Lvl. start ; Lvl. end
661.6553 ; 0.0030 ; 85.01 ; 0.20 ; g ; Ba-137 ; 2 ; 0
...
```

`Type` codes: `g` (gamma), `g511` (511 keV annihilation), `XKa1`/`XKa2`/`XK'b1`/etc
(K X-rays), `XL` (L X-rays, often summed), `a` (alpha — skipped by g4gamma),
`b-`/`b+`/`ec` (particle emissions — skipped). Intensity is in *percent*.

**Important physical convention:** LARA's emission lists are flattened over
the cascade — each parent's file lists every gamma observed when *that
nuclide* decays, but stops at the direct daughter. So Bi-214's file
includes the 609 keV gamma (emitted by the Po-214 daughter level scheme
during de-excitation following β decay) but NOT the 351 keV gamma from
its grandchild Pb-214 — that's in `Pb-214.lara.txt`. g4gamma still walks
the chain and weights each parent's emissions by activity.

## License & citation

LARA data tables are made publicly available by LNHB on
<http://www.lnhb.fr/Laraweb/>. The site requests citation of the relevant
DDEP evaluation report when used in publications. Treat them
the same way you'd treat ENSDF or NuDat data: cite the source, redistribute
with attribution, don't claim authorship.

For thesis/publication purposes, cite:
- M.M. Bé et al., *Table de Radionucléides* (Monographie BIPM-5),
  Bureau International des Poids et Mesures (multiple volumes, 2004-present)
- Per-nuclide evaluations are linked from each
  `http://www.lnhb.fr/nuclides/<NAME>_com.pdf` evaluation comments file.
