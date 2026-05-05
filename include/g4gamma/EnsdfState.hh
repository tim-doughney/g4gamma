// EnsdfState.hh -- parses ENSDFSTATE.dat from $G4ENSDFSTATEDATA, which is the
// canonical source of nuclide mean lives in Geant4. The half-life column in
// the RadioactiveDecay zZ.aA files is a placeholder/dummy and should not be
// used; the real values come from this file.
//
// File format (per G4NuclideTable::GenerateNuclide):
//   <Z> <A> <E_keV> <floatFlag> <meanLife_ns> <2J> <muMag>
// repeated per nuclide level. Records with E=0 and no float flag are ground
// states. Records with E>0 above the mean-life threshold are isomers (the
// threshold is set by the user; we adopt the same 1 ns convention as the
// PhysicsList).
#pragma once

#include "g4gamma/IsotopeKey.hh"
#include <string>
#include <vector>
#include <unordered_map>

namespace g4gamma {

struct EnsdfLevel {
    double excitation;      // internal units (MeV)
    double meanLife;        // internal units (ns)
    std::string floatFlag;  // "  " for none, otherwise "+X" / "+Y" etc.
};

class EnsdfStateLoader {
public:
    explicit EnsdfStateLoader(std::string ensdfStateDir);

    // Load all levels for (Z, A). Returns nullptr if file/data missing.
    const std::vector<EnsdfLevel>* load(int Z, int A);

    // Convenience: find the meanLife of the M-th isomer of (Z, A), where M=0
    // means ground state (excitation == 0), M=1 first isomer above threshold,
    // etc. Returns 0 if not found (callers should treat as unknown / stable).
    double meanLife(const IsotopeKey& key, double meanLifeThreshold);

    // Map an excitation energy (internal units) to an isomer index M for
    // (Z, A). M=0 = ground state. M>=1 = successive isomers above the
    // lifetime threshold. Returns -1 if no match.
    int excitationToM(int Z, int A, double excitation, double thresh,
                      double tolerance);

    bool ready() {
        if (!fLoadAttempted) loadAll();
        return fLoaded;
    }

private:
    std::string fPath;          // full path to ENSDFSTATE.dat
    bool fLoaded = false;
    bool fLoadAttempted = false;
    // key = 1000*A + Z
    std::unordered_map<int, std::vector<EnsdfLevel>> fCache;
    void loadAll();
};

} // namespace g4gamma
