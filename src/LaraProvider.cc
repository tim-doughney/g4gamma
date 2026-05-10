#include "g4gamma/LaraProvider.hh"
#include "g4gamma/Units.hh"
#include <filesystem>
#include <fstream>
#include <sstream>
#include <stdexcept>
#include <cstdlib>
#include <cstring>
#include <cmath>
#include <algorithm>
#include <zlib.h>

namespace fs = std::filesystem;

namespace g4gamma {

namespace {

static const std::unordered_map<std::string, int>& elementZ() {
    static const std::unordered_map<std::string, int> kTable = {
        {"H",1},{"He",2},{"Li",3},{"Be",4},{"B",5},{"C",6},{"N",7},{"O",8},
        {"F",9},{"Ne",10},{"Na",11},{"Mg",12},{"Al",13},{"Si",14},{"P",15},
        {"S",16},{"Cl",17},{"Ar",18},{"K",19},{"Ca",20},{"Sc",21},{"Ti",22},
        {"V",23},{"Cr",24},{"Mn",25},{"Fe",26},{"Co",27},{"Ni",28},{"Cu",29},
        {"Zn",30},{"Ga",31},{"Ge",32},{"As",33},{"Se",34},{"Br",35},{"Kr",36},
        {"Rb",37},{"Sr",38},{"Y",39},{"Zr",40},{"Nb",41},{"Mo",42},{"Tc",43},
        {"Ru",44},{"Rh",45},{"Pd",46},{"Ag",47},{"Cd",48},{"In",49},{"Sn",50},
        {"Sb",51},{"Te",52},{"I",53},{"Xe",54},{"Cs",55},{"Ba",56},{"La",57},
        {"Ce",58},{"Pr",59},{"Nd",60},{"Pm",61},{"Sm",62},{"Eu",63},{"Gd",64},
        {"Tb",65},{"Dy",66},{"Ho",67},{"Er",68},{"Tm",69},{"Yb",70},{"Lu",71},
        {"Hf",72},{"Ta",73},{"W",74},{"Re",75},{"Os",76},{"Ir",77},{"Pt",78},
        {"Au",79},{"Hg",80},{"Tl",81},{"Pb",82},{"Bi",83},{"Po",84},{"At",85},
        {"Rn",86},{"Fr",87},{"Ra",88},{"Ac",89},{"Th",90},{"Pa",91},{"U",92},
        {"Np",93},{"Pu",94},{"Am",95},{"Cm",96},{"Bk",97},{"Cf",98},{"Es",99},
        {"Fm",100},{"Md",101},{"No",102},{"Lr",103},{"Rf",104},{"Db",105},
        {"Sg",106},{"Bh",107},{"Hs",108},{"Mt",109},{"Ds",110},{"Rg",111},
        {"Cn",112},{"Nh",113},{"Fl",114},{"Mc",115},{"Lv",116},{"Ts",117},{"Og",118},
    };
    return kTable;
}

static std::string elementName(int Z) {
    for (const auto& [sym, z] : elementZ()) {
        if (z == Z) return sym;
    }
    return {};
}

static void trim(std::string& s) {
    auto notSpace = [](unsigned char c) { return !std::isspace(c); };
    s.erase(s.begin(), std::find_if(s.begin(), s.end(), notSpace));
    s.erase(std::find_if(s.rbegin(), s.rend(), notSpace).base(), s.end());
}

static std::vector<std::string> splitSemis(const std::string& line) {
    std::vector<std::string> out;
    size_t start = 0;
    while (start <= line.size()) {
        size_t end = line.find(';', start);
        if (end == std::string::npos) end = line.size();
        std::string f = line.substr(start, end - start);
        trim(f);
        out.push_back(std::move(f));
        if (end == line.size()) break;
        start = end + 1;
    }
    return out;
}

static double parseDouble(const std::string& s, double def = 0.0) {
    if (s.empty()) return def;
    try { return std::stod(s); } catch (...) { return def; }
}

static DecayMode parseLaraMode(const std::string& raw) {
    std::string s = raw;
    if (!s.empty() && s.front() == '(') s.erase(s.begin());
    if (!s.empty() && s.back() == ')') s.pop_back();
    trim(s);
    auto sep = s.find_first_of(",+ ");
    if (sep != std::string::npos) s = s.substr(0, sep);
    std::string lower; lower.reserve(s.size());
    for (char c : s) lower.push_back(std::tolower(static_cast<unsigned char>(c)));

    if (lower == "b-")    return DecayMode::BetaMinus;
    if (lower == "b+")    return DecayMode::BetaPlus;
    if (lower == "ec")    return DecayMode::EC;
    if (lower == "i.t."   || lower == "it") return DecayMode::IT;
    if (lower == "alpha"  || lower == "a")  return DecayMode::Alpha;
    if (lower == "sf")    return DecayMode::SpFission;
    if (lower == "p")     return DecayMode::Proton;
    if (lower == "n")     return DecayMode::Neutron;
    return DecayMode::Unknown;
}

struct EmissionParse {
    bool        keep;
    EmissionType type;
};
static EmissionParse parseEmissionType(const std::string& t) {
    if (t == "g")        return {true, EmissionType::Gamma};
    if (t == "g511")     return {true, EmissionType::AnnihilationPair};
    if (!t.empty() && t[0] == 'X') return {true, EmissionType::XRay};
    return {false, EmissionType::Gamma};
}

// ---------- tar reading ---------------------------------------------------
//
// USTAR format: 512-byte header blocks followed by 512-byte data blocks
// (rounded up). We only need to find regular files matching `*.lara.txt`
// and record their offset+length within the in-memory archive.
//
// We do NOT extract -- file contents stay in fTarBlob and we hand out
// substring views via offset/length.

static unsigned long parseOctal(const char* p, size_t n) {
    unsigned long v = 0;
    for (size_t i = 0; i < n; ++i) {
        char c = p[i];
        if (c == 0 || c == ' ') continue;
        if (c < '0' || c > '7') return v;
        v = (v << 3) | (c - '0');
    }
    return v;
}

// Extract the file name from a USTAR header. Honours both the `name` field
// (offset 0, 100 bytes) and the `prefix` field (offset 345, 155 bytes) used
// for long paths. Modern POSIX tar also supports PaxHeader extensions which
// we treat as "skip".
static std::string tarFilename(const char* hdr) {
    std::string name(hdr, strnlen(hdr, 100));
    // Check magic for ustar
    if (std::memcmp(hdr + 257, "ustar", 5) == 0) {
        std::string prefix(hdr + 345, strnlen(hdr + 345, 155));
        if (!prefix.empty()) name = prefix + "/" + name;
    }
    return name;
}

} // namespace


// ---------- Symbol/key conversions ---------------------------------------

IsotopeKey LaraProvider::symbolToKey(const std::string& sym) {
    if (sym.empty()) return IsotopeKey{-1, -1, -1};
    std::string s = sym;
    int M = 0;
    if (!s.empty() && (s.back() == 'm' || s.back() == 'n')) {
        char marker = s.back();
        s.pop_back();
        if (!s.empty() && std::isdigit(static_cast<unsigned char>(s.back()))) {
            s.push_back(marker);
        } else {
            M = 1;
        }
    }
    size_t dash = s.find('-');
    std::string elem, astr;
    if (dash != std::string::npos) {
        elem = s.substr(0, dash);
        astr = s.substr(dash + 1);
    } else {
        size_t i = 0;
        while (i < s.size() && std::isalpha(static_cast<unsigned char>(s[i]))) ++i;
        elem = s.substr(0, i);
        astr = s.substr(i);
    }
    if (elem.empty() || astr.empty()) return IsotopeKey{-1, -1, -1};
    auto it = elementZ().find(elem);
    if (it == elementZ().end()) return IsotopeKey{-1, -1, -1};
    int Z = it->second;
    int A;
    try { A = std::stoi(astr); } catch (...) { return IsotopeKey{-1, -1, -1}; }
    return IsotopeKey{Z, A, M};
}

std::string LaraProvider::keyToSymbol(const IsotopeKey& k) {
    std::string elem = elementName(k.Z);
    if (elem.empty()) return {};
    std::string s = elem + "-" + std::to_string(k.A);
    if (k.M == 1) s += "m";
    else if (k.M >= 2) { s += "m"; s += std::to_string(k.M); }
    return s;
}

// ---------- Path resolution ----------------------------------------------

std::string LaraProvider::locateData(const std::string& hint) {
    auto isTarFile = [](const std::string& p) {
        return !p.empty() && fs::exists(p) && fs::is_regular_file(p) &&
               (p.size() >= 4 &&
                (p.substr(p.size()-4) == ".tar" ||
                 (p.size() >= 7 && p.substr(p.size()-7) == ".tar.gz")));
    };
    auto isLooseDir = [](const std::string& p) {
        if (p.empty() || !fs::exists(p) || !fs::is_directory(p)) return false;
        for (auto& e : fs::directory_iterator(p)) {
            if (e.is_regular_file() &&
                e.path().filename().string().find(".lara.txt") != std::string::npos)
                return true;
        }
        return false;
    };

    // 1. Explicit hint -- could be a tarball OR a directory.
    if (isTarFile(hint))  return hint;
    if (isLooseDir(hint)) return hint;

    // 2. Env var
    if (const char* env = std::getenv("LARA_DATA_DIR")) {
        std::string p = env;
        if (isTarFile(p))  return p;
        if (isLooseDir(p)) return p;
    }

    // 3-5. Search a list of likely locations -- prefer tarballs over dirs
    static const std::vector<std::string> roots = {
        "data/lara", "../data/lara", "../../data/lara",
        "/usr/local/share/g4gamma/lara",
        "/usr/share/g4gamma/lara",
    };
    for (const auto& r : roots) {
        for (const std::string suffix : {"/lara.tar.gz", "/lara.tar"}) {
            std::string p = r + suffix;
            if (isTarFile(p)) return p;
        }
        if (isLooseDir(r)) return r;
    }
    return {};
}

// ---------- gzip decompression -------------------------------------------

std::string LaraProvider::gunzip(const std::string& path) {
    gzFile gf = gzopen(path.c_str(), "rb");
    if (!gf) throw std::runtime_error("LaraProvider: cannot open " + path);
    std::string out;
    char buf[64 * 1024];
    int n;
    while ((n = gzread(gf, buf, sizeof(buf))) > 0) {
        out.append(buf, n);
    }
    gzclose(gf);
    return out;
}

// ---------- Construction & indexing --------------------------------------

LaraProvider::LaraProvider(const std::string& hint) {
    fSource = locateData(hint);
    if (fSource.empty()) {
        throw std::runtime_error(
            "LaraProvider: could not locate LARA data. Pass an explicit path "
            "to a *.lara.txt directory or .tar/.tar.gz tarball, set "
            "$LARA_DATA_DIR, or place files in <repo>/data/lara/.");
    }
    if (fs::is_regular_file(fSource)) {
        fIsTar  = true;
        fIsGzip = (fSource.size() >= 7 && fSource.substr(fSource.size()-7) == ".tar.gz");
        scanTarball();
    } else {
        fIsTar = false;
        scanDirectory();
    }
}

void LaraProvider::scanDirectory() {
    fAvailable.clear();
    for (auto& entry : fs::directory_iterator(fSource)) {
        if (!entry.is_regular_file()) continue;
        std::string fname = entry.path().filename().string();
        const std::string suffix = ".lara.txt";
        if (fname.size() <= suffix.size()) continue;
        if (fname.compare(fname.size() - suffix.size(), suffix.size(), suffix) != 0)
            continue;
        std::string sym = fname.substr(0, fname.size() - suffix.size());
        IsotopeKey key = symbolToKey(sym);
        if (key.Z < 0) continue;
        Entry e;
        e.filename = fname;
        fAvailable[key] = e;
    }
}

void LaraProvider::scanTarball() {
    // Load whole archive into memory. Either the raw tar bytes or the
    // gunzipped tar bytes -- both end up in fTarBlob.
    if (fIsGzip) {
        fTarBlob = gunzip(fSource);
    } else {
        std::ifstream in(fSource, std::ios::binary);
        if (!in) throw std::runtime_error("LaraProvider: cannot open " + fSource);
        std::ostringstream ss;
        ss << in.rdbuf();
        fTarBlob = ss.str();
    }

    fAvailable.clear();
    const char* data = fTarBlob.data();
    size_t total     = fTarBlob.size();
    size_t pos       = 0;
    const size_t kBlock = 512;
    const std::string suffix = ".lara.txt";

    while (pos + kBlock <= total) {
        const char* hdr = data + pos;
        // End-of-archive: two zero blocks. Detect by leading nul.
        if (hdr[0] == 0) break;

        std::string name = tarFilename(hdr);
        char typeflag = hdr[156];
        unsigned long sz = parseOctal(hdr + 124, 12);

        // Advance pos past header + data (rounded up to block).
        size_t dataStart = pos + kBlock;
        size_t dataEnd   = dataStart + ((sz + kBlock - 1) / kBlock) * kBlock;

        // Only care about regular files (typeflag '0' or '\0') ending in .lara.txt
        bool isRegular = (typeflag == '0' || typeflag == 0);
        if (isRegular && name.size() > suffix.size() &&
            name.compare(name.size() - suffix.size(), suffix.size(), suffix) == 0) {
            // Strip leading "./" and any directory components -- we only key
            // off the basename.
            size_t slash = name.find_last_of('/');
            std::string base = (slash == std::string::npos) ? name : name.substr(slash+1);
            std::string sym = base.substr(0, base.size() - suffix.size());
            IsotopeKey key = symbolToKey(sym);
            if (key.Z >= 0) {
                Entry e;
                e.filename = base;
                e.offset   = dataStart;
                e.length   = sz;
                fAvailable[key] = e;
            }
        }
        pos = dataEnd;
    }
}

// ---------- Per-entry access ---------------------------------------------

std::string LaraProvider::readEntry(const Entry& e) const {
    if (fIsTar) {
        if (e.offset + e.length > fTarBlob.size()) return {};
        return fTarBlob.substr(e.offset, e.length);
    } else {
        std::ifstream in(fSource + "/" + e.filename);
        if (!in) return {};
        std::ostringstream ss;
        ss << in.rdbuf();
        return ss.str();
    }
}

const ParentDecayInfo* LaraProvider::get(const IsotopeKey& key) {
    auto it = fByKey.find(key);
    if (it != fByKey.end()) return &it->second;
    if (fAttempted[key]) return nullptr;
    fAttempted[key] = true;

    auto avail = fAvailable.find(key);
    if (avail == fAvailable.end()) return nullptr;

    std::string content = readEntry(avail->second);
    if (content.empty()) return nullptr;
    if (!parseContent(key, content)) return nullptr;

    auto it2 = fByKey.find(key);
    return (it2 == fByKey.end()) ? nullptr : &it2->second;
}

// ---------- Parser -------------------------------------------------------

bool LaraProvider::parseContent(const IsotopeKey& key, const std::string& content) {
    using namespace units;
    std::istringstream in(content);

    ParentDecayInfo info;
    info.isotope  = key;
    info.stable   = false;
    info.meanLife = 0.0;

    // Newer LARA files (NIST 2025+ format) include cascade gammas from all
    // daughter nuclides in one file, with a trailing "Parent" column that
    // identifies which nuclide in the chain each emission belongs to.  When
    // individual daughter files are also present we must skip emissions whose
    // Parent is not this file's primary nuclide, otherwise they are counted
    // twice.  Build the expected primary-nuclide symbol for comparison.
    std::string primarySym = elementName(key.Z);
    if (!primarySym.empty()) {
        primarySym += "-" + std::to_string(key.A);
        if (key.M > 0) primarySym += "m";
    }

    struct DaughterRec {
        DecayMode  mode;
        IsotopeKey key;
        std::string symbol;
        double     branchPercent;
    };
    std::vector<DaughterRec> daughters;

    std::string line;
    bool inEmissions = false;
    bool sawEmissionsHeader = false;

    while (std::getline(in, line)) {
        if (!line.empty() && line.back() == '\r') line.pop_back();
        if (line.empty()) continue;
        if (line[0] == '=' || line[0] == '-') continue;

        if (!inEmissions) {
            if (line.find("Energy (keV)") != std::string::npos &&
                line.find("Intensity") != std::string::npos) {
                sawEmissionsHeader = true;
                inEmissions = true;
                continue;
            }
            if (line.find("No emissions") != std::string::npos) {
                info.stable = true;
                fByKey.emplace(key, std::move(info));
                return true;
            }

            auto fields = splitSemis(line);
            if (fields.empty()) continue;
            const std::string& tag = fields[0];

            if (tag == "Half-life (s)") {
                if (fields.size() >= 2) {
                    double hl = parseDouble(fields[1], 0.0);
                    if (hl > 0) info.meanLife = (hl * second) / std::log(2.0);
                }
            } else if (tag == "Daughter(s)") {
                size_t i = 1;
                while (i + 2 < fields.size()) {
                    DaughterRec d;
                    d.mode          = parseLaraMode(fields[i]);
                    d.symbol        = fields[i+1];
                    d.key           = symbolToKey(fields[i+1]);
                    d.branchPercent = parseDouble(fields[i+2], 0.0);
                    if (d.key.Z >= 0 && d.branchPercent > 0.0) {
                        daughters.push_back(d);
                    }
                    i += 3;
                }
            }
        } else {
            auto fields = splitSemis(line);
            if (fields.size() < 5) continue;
            // New cascade format: 9th field is "Parent" — skip emissions that
            // belong to a daughter nuclide (they are counted in that nuclide's
            // own file).
            if (fields.size() >= 9 && !primarySym.empty()) {
                const std::string& parent = fields[8];
                if (!parent.empty() && parent != primarySym) continue;
            }
            double E_keV = parseDouble(fields[0], -1.0);
            double I_pct = parseDouble(fields[2], 0.0);
            const std::string& type = fields[4];
            std::string origin = fields.size() > 5 ? fields[5] : "";
            if (E_keV <= 0.0 || I_pct <= 0.0) continue;
            auto pe = parseEmissionType(type);
            if (!pe.keep) continue;

            Emission em;
            em.type      = pe.type;
            em.energy    = E_keV * keV;
            em.intensity = I_pct / 100.0;

            int matchIdx = -1;
            for (size_t i = 0; i < daughters.size(); ++i) {
                if (!origin.empty() && origin == daughters[i].symbol) {
                    matchIdx = static_cast<int>(i); break;
                }
            }
            if (matchIdx < 0) {
                double best = -1.0;
                for (size_t i = 0; i < daughters.size(); ++i) {
                    if (daughters[i].branchPercent > best) {
                        best = daughters[i].branchPercent;
                        matchIdx = static_cast<int>(i);
                    }
                }
            }
            if (matchIdx < 0) continue;
            while (info.branches.size() <= static_cast<size_t>(matchIdx)) {
                size_t i = info.branches.size();
                DecayBranch br;
                br.mode           = daughters[i].mode;
                br.branchingRatio = daughters[i].branchPercent / 100.0;
                br.daughter       = daughters[i].key;
                info.branches.push_back(br);
            }
            info.branches[matchIdx].emissions.push_back(em);
        }
    }

    if (!sawEmissionsHeader) info.stable = true;

    if (info.branches.empty() && !daughters.empty()) {
        for (const auto& d : daughters) {
            DecayBranch br;
            br.mode           = d.mode;
            br.branchingRatio = d.branchPercent / 100.0;
            br.daughter       = d.key;
            info.branches.push_back(br);
        }
    }
    if (info.branches.empty()) info.stable = true;

    fByKey.emplace(key, std::move(info));
    return true;
}

} // namespace g4gamma
