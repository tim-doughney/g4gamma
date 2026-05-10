// AugerData.hh -- parses Geant4 G4LEDATA/auger/au-tr-pr-Z.dat to obtain
// per-element Auger/Coster-Kronig transition data following an inner-shell vacancy.
//
// File format (per G4AugerData::LoadData):
//   Block header: "<id> <id> <id> <id>"  (4 identical integers = vacancy EADL shell id)
//   Transition rows: "<secShell1> <secShell2> <prob> <energy_MeV>"
//   Block separator: "-1 -1 -1 -1"
//   EOF marker:      "-2 -2 -2 -2"
//
// The probabilities sum to (1 - omega) for each vacancy shell, where omega is
// the fluorescence yield from fl-tr-pr. Together, fl-tr-pr and au-tr-pr account
// for all deexcitation paths from a vacancy.
//
// Each Auger/CK transition creates two secondary vacancies (secShell1, secShell2),
// both in shells outer to the initial vacancy. We store only secShell1, secShell2,
// and prob; the Auger electron energy is not needed for photon-yield computation.
#pragma once

#include <string>
#include <vector>
#include <unordered_map>

namespace g4gamma {

struct AugerTransition {
    int    secShell1;    // first secondary vacancy EADL shell ID
    int    secShell2;    // second secondary vacancy EADL shell ID
    double prob;         // transition probability (NOT renormalised; sums to 1 - omega over block)
};

struct AugerVacancy {
    int vacancyShell;                          // EADL shell ID of the initial vacancy
    std::vector<AugerTransition> transitions;
};

class AugerDataLoader {
public:
    // augerDir should be "$G4LEDATA/auger" (the directory containing au-tr-pr-Z.dat files).
    explicit AugerDataLoader(std::string augerDir);

    // Load all vacancy blocks for element Z. Returns nullptr if no data file.
    const std::vector<AugerVacancy>* load(int Z);

    // Find a specific vacancy block by EADL shell ID (1=K, 3=L1, 5=L2, 6=L3, …).
    // Returns nullptr if Z has no data or the shell is absent.
    const AugerVacancy* findVacancy(int Z, int eadlShellId);

    const std::string& dir() const { return fDir; }

private:
    std::string fDir;
    std::unordered_map<int, std::vector<AugerVacancy>> fCache;
};

} // namespace g4gamma
