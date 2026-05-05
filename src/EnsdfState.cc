#include "g4gamma/EnsdfState.hh"
#include "g4gamma/Units.hh"
#include <fstream>
#include <sstream>
#include <iostream>
#include <algorithm>
#include <cmath>

namespace g4gamma {

EnsdfStateLoader::EnsdfStateLoader(std::string dir) {
    if (dir.empty()) return;
    fPath = dir + "/ENSDFSTATE.dat";
}

void EnsdfStateLoader::loadAll() {
    if (fLoadAttempted) return;
    fLoadAttempted = true;
    if (fPath.empty()) return;
    std::ifstream in(fPath);
    if (!in) return;

    using namespace units;
    int Z, A;
    double E_keV, life_ns, mu;
    int twoJ;
    std::string floatFlag;

    while (in >> Z >> A >> E_keV >> floatFlag >> life_ns >> twoJ >> mu) {
        EnsdfLevel lvl;
        lvl.excitation = E_keV * keV;
        // ENSDFSTATE encodes "infinite" lifetimes as a negative or zero value
        // depending on the dataset version; treat <=0 as "stable / unknown".
        lvl.meanLife = (life_ns > 0.0) ? life_ns : 0.0;
        lvl.floatFlag = floatFlag;
        int code = 1000 * A + Z;
        fCache[code].push_back(lvl);
        // suppress unused warning
        (void)twoJ; (void)mu;
    }
    fLoaded = true;
}

const std::vector<EnsdfLevel>* EnsdfStateLoader::load(int Z, int A) {
    if (!fLoadAttempted) loadAll();
    if (!fLoaded) return nullptr;
    auto it = fCache.find(1000 * A + Z);
    if (it == fCache.end()) return nullptr;
    return &it->second;
}

double EnsdfStateLoader::meanLife(const IsotopeKey& key, double thresh) {
    const auto* v = load(key.Z, key.A);
    if (!v) return 0.0;
    // Build the M index by walking levels in order of excitation. M=0 is the
    // ground state; M=1, 2, ... are successive levels above `thresh`.
    int curM = 0;
    for (const auto& lvl : *v) {
        // Skip nothing here -- we map M directly
        if (curM == key.M) return lvl.meanLife;
        // Advance M only if the level qualifies as a tracked isomer:
        if (lvl.excitation > 0.0 && lvl.meanLife > thresh) {
            ++curM;
            if (curM == key.M) return lvl.meanLife;
        } else if (curM == 0 && lvl.excitation == 0.0) {
            // first ground-state record; if the user asked for M=0, returned
            // above.
        }
    }
    return 0.0;
}

} // namespace g4gamma
