#include "g4gamma/DataPath.hh"
#include <cstdlib>
#include <fstream>
#include <sstream>
#include <stdexcept>
#include <filesystem>
#include <regex>
#include <vector>

namespace fs = std::filesystem;

namespace g4gamma {

std::string DataPath::getEnv(const char* name) {
    const char* v = std::getenv(name);
    return v ? std::string(v) : std::string();
}

// Parse a sh-style "export VAR=...something..." line out of a geant4.sh file.
// Supports backquote / $(cd ... ; pwd) constructs by evaluating them in a
// best-effort way (we only resolve a fixed set of patterns rather than spawning
// a real shell).
std::string DataPath::parseGeant4Sh(const std::string& path, const std::string& var) {
    std::ifstream in(path);
    if (!in) return {};

    // First pass: gather all "name=value" assignments (commented or not).
    // We only honour uncommented "export FOO=..." or "FOO=..." lines for the
    // requested var, but we *do* read commented assignments for variables we
    // need for substitution (like geant4_envbindir).
    std::string line;
    std::string scriptDir = fs::absolute(fs::path(path)).parent_path().string();
    std::vector<std::pair<std::string, std::string>> assignments;

    auto strip_inline_comment = [](std::string s) {
        // Remove trailing comment but be cautious about # inside quotes.
        bool in_single = false, in_double = false, in_back = false;
        for (size_t i = 0; i < s.size(); ++i) {
            char c = s[i];
            if (c == '\'' && !in_double && !in_back) in_single = !in_single;
            else if (c == '"' && !in_single && !in_back) in_double = !in_double;
            else if (c == '`' && !in_single && !in_double) in_back = !in_back;
            else if (c == '#' && !in_single && !in_double && !in_back) {
                return s.substr(0, i);
            }
        }
        return s;
    };

    auto trim = [](std::string s) {
        size_t a = s.find_first_not_of(" \t\r\n");
        if (a == std::string::npos) return std::string();
        size_t b = s.find_last_not_of(" \t\r\n");
        return s.substr(a, b - a + 1);
    };

    while (std::getline(in, line)) {
        std::string s = trim(line);
        if (s.empty()) continue;

        // Detect comment-only line; we still want to scan for assignments
        // (some installs comment out RadioactiveData= line). Strip a leading '#'
        // from a few likely lines.
        bool was_comment = false;
        if (s[0] == '#') {
            s = trim(s.substr(1));
            was_comment = true;
            if (s.empty()) continue;
        }

        if (s.rfind("export ", 0) == 0) s = s.substr(7);
        s = strip_inline_comment(s);
        s = trim(s);

        auto eq = s.find('=');
        if (eq == std::string::npos) continue;
        std::string name = trim(s.substr(0, eq));
        std::string value = trim(s.substr(eq + 1));
        if (name.empty()) continue;
        // Strip outer quotes
        if (value.size() >= 2 &&
            ((value.front() == '"' && value.back() == '"') ||
             (value.front() == '\'' && value.back() == '\''))) {
            value = value.substr(1, value.size() - 2);
        }
        // Only record the *first* uncommented assignment, but always record
        // for backquote pattern resolution if commented.
        (void)was_comment;
        assignments.emplace_back(name, value);
    }

    // Build a map by walking assignments left-to-right and substituting refs
    // using whatever was assigned earlier. The script's first line typically
    // defines geant4_envbindir using `dirname $BASH_SOURCE` which we cannot
    // execute -- so we give the script-derived value precedence by scanning
    // for the latest assignment of geant4_envbindir and replacing its value
    // with the actual script directory.
    for (auto& [name, value] : assignments) {
        if (name == "geant4_envbindir") {
            value = scriptDir;
        }
    }
    std::vector<std::pair<std::string,std::string>> resolved;
    auto lookup = [&](const std::string& key) -> std::string {
        // first check resolved (prefer last value)
        for (auto it = resolved.rbegin(); it != resolved.rend(); ++it) {
            if (it->first == key) return it->second;
        }
        // fallback: env
        return getEnv(key.c_str());
    };

    auto substitute = [&](std::string v) {
        // Replace $VAR and ${VAR}
        std::string out;
        for (size_t i = 0; i < v.size();) {
            if (v[i] == '$' && i + 1 < v.size()) {
                if (v[i+1] == '{') {
                    auto end = v.find('}', i + 2);
                    if (end == std::string::npos) { out += v[i++]; continue; }
                    std::string name = v.substr(i + 2, end - i - 2);
                    out += lookup(name);
                    i = end + 1;
                } else if (std::isalpha(static_cast<unsigned char>(v[i+1])) || v[i+1] == '_') {
                    size_t j = i + 1;
                    while (j < v.size() &&
                           (std::isalnum(static_cast<unsigned char>(v[j])) || v[j] == '_'))
                        ++j;
                    std::string name = v.substr(i + 1, j - i - 1);
                    out += lookup(name);
                    i = j;
                } else {
                    out += v[i++];
                }
            } else {
                out += v[i++];
            }
        }
        // Resolve `cd X > /dev/null ; pwd` and $(cd X > /dev/null ; pwd):
        // we only care about the directory expression.
        std::regex back_re(R"(`\s*cd\s+([^>]+?)\s*>\s*/dev/null\s*;\s*pwd\s*`)");
        std::regex paren_re(R"(\$\(\s*cd\s+([^>]+?)\s*>\s*/dev/null\s*;\s*pwd\s*\))");
        std::smatch m;
        while (std::regex_search(out, m, back_re)) {
            std::string dir = m[1].str();
            try {
                fs::path p = fs::absolute(fs::path(dir)).lexically_normal();
                out = m.prefix().str() + p.string() + m.suffix().str();
            } catch (...) { break; }
        }
        while (std::regex_search(out, m, paren_re)) {
            std::string dir = m[1].str();
            try {
                fs::path p = fs::absolute(fs::path(dir)).lexically_normal();
                out = m.prefix().str() + p.string() + m.suffix().str();
            } catch (...) { break; }
        }
        return out;
    };

    // Seed with geant4_envbindir (always the directory of the script)
    resolved.emplace_back("geant4_envbindir", scriptDir);

    for (auto& [name, value] : assignments) {
        resolved.emplace_back(name, substitute(value));
    }

    // Now look up the variable the caller asked for.
    for (auto it = resolved.rbegin(); it != resolved.rend(); ++it) {
        if (it->first == var) return it->second;
    }
    return {};
}

std::string DataPath::findVersionedSubdir(const std::string& parent,
                                          const std::string& prefix) {
    if (parent.empty() || !fs::exists(parent)) return {};
    std::string best;
    for (const auto& entry : fs::directory_iterator(parent)) {
        if (!entry.is_directory()) continue;
        std::string name = entry.path().filename().string();
        if (name.rfind(prefix, 0) == 0) {
            // pick the last (lexicographically largest) match -- highest version
            if (name > best) best = entry.path().string();
        }
    }
    return best;
}

static std::string resolveDataset(const std::string& envVar,
                                   const std::string& sh_hint,
                                   const std::string& dirPrefix) {
    // 1. Direct env var
    auto v = DataPath::getEnv(envVar.c_str());
    if (!v.empty() && fs::exists(v)) return v;

    // 2. Parse geant4.sh
    std::vector<std::string> shCandidates;
    if (!sh_hint.empty()) shCandidates.push_back(sh_hint);
    shCandidates.push_back("/opt/geant4/bin/geant4.sh");
    // also try GEANT4_DATA_DIR's parent / .. / bin / geant4.sh
    auto gdd = DataPath::getEnv("GEANT4_DATA_DIR");
    if (!gdd.empty()) {
        try {
            fs::path p = fs::path(gdd).parent_path().parent_path() / "bin" / "geant4.sh";
            shCandidates.push_back(p.string());
        } catch (...) {}
    }

    for (const auto& sh : shCandidates) {
        if (!fs::exists(sh)) continue;
        auto v2 = DataPath::parseGeant4Sh(sh, envVar);
        if (!v2.empty() && fs::exists(v2)) return v2;
        // Some geant4.sh comment out the per-dataset vars but set GEANT4_DATA_DIR.
        auto gddVal = DataPath::parseGeant4Sh(sh, "GEANT4_DATA_DIR");
        if (!gddVal.empty()) {
            auto found = DataPath::findVersionedSubdir(gddVal, dirPrefix);
            if (!found.empty()) return found;
        }
    }

    // 3. GEANT4_DATA_DIR fallback
    if (!gdd.empty()) {
        auto found = DataPath::findVersionedSubdir(gdd, dirPrefix);
        if (!found.empty()) return found;
    }

    return {};
}

std::string DataPath::radioactiveDecayDir(const std::string& sh) {
    auto p = resolveDataset("G4RADIOACTIVEDATA", sh, "RadioactiveDecay");
    if (p.empty()) {
        throw std::runtime_error(
            "Could not locate RadioactiveDecay data. "
            "Set G4RADIOACTIVEDATA, or pass a path to geant4.sh, or set GEANT4_DATA_DIR.");
    }
    return p;
}

std::string DataPath::photonEvaporationDir(const std::string& sh) {
    auto p = resolveDataset("G4LEVELGAMMADATA", sh, "PhotonEvaporation");
    if (p.empty()) {
        throw std::runtime_error(
            "Could not locate PhotonEvaporation data. "
            "Set G4LEVELGAMMADATA, or pass a path to geant4.sh, or set GEANT4_DATA_DIR.");
    }
    return p;
}

std::string DataPath::lowEnergyDir(const std::string& sh) {
    auto p = resolveDataset("G4LEDATA", sh, "G4EMLOW");
    if (p.empty()) {
        throw std::runtime_error(
            "Could not locate G4LEDATA. "
            "Set G4LEDATA, or pass a path to geant4.sh, or set GEANT4_DATA_DIR.");
    }
    return p;
}

std::string DataPath::ensdfStateDir(const std::string& sh) {
    // Optional dataset; return empty if not found rather than throwing.
    return resolveDataset("G4ENSDFSTATEDATA", sh, "G4ENSDFSTATE");
}

} // namespace g4gamma
