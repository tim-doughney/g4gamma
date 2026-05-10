// FluorData.hh -- parses Geant4 G4LEDATA/<fluor_dir>/fl-tr-pr-Z.dat to obtain
// per-element X-ray fluorescence transitions following an inner-shell vacancy.
//
// File format (per G4FluoData::LoadData): rows of
//     <originating_shell_id>  <transition_probability>  <transition_energy_MeV>
// grouped by "vacancy block" (which shell has the vacancy). End of vacancy
// block is marked by a row of -1; end of file is -2.
//
// We use the *probabilities* as-is (they are already fluorescence-yield
// weighted -- i.e. they sum to omega_K for K vacancy, omega_L for L, etc.).
//
// To compute the X-ray spectrum following an EC decay or an internal-conversion
// transition, we determine the shell vacancy (K or L or M...) and then sum up
// fluorescence transitions for that vacancy. Auger emission is ignored (no
// gammas produced). Cascading vacancies (e.g. K vacancy filled by L electron
// then L vacancy filled by M electron) are handled approximately by treating
// each vacancy independently using its own fluorescence yield -- this is a
// simplification of Geant4's full atomic relaxation cascade.
#pragma once

#include <string>
#include <vector>
#include <unordered_map>

namespace g4gamma {

struct FluorTransition {
    int    originShell;     // shell id of the electron that fills the vacancy
    double prob;            // transition probability (NOT renormalised to sum=1)
    double energy;          // X-ray photon energy (internal units, MeV)
};

// One per vacancy shell (K=0, L1=1, L2=2, L3=3, M1=4, ...).
struct FluorVacancy {
    int vacancyShell;
    std::vector<FluorTransition> transitions;
};

class FluorDataLoader {
public:
    // ledataDir is $G4LEDATA. fluorSubdir is one of "/fluor", "/fluor_Bearden",
    // "/fluor_ANSTO", "/fluor_XDB_EADL". Matches Geant4's default selection.
    FluorDataLoader(std::string ledataDir, std::string fluorSubdir = "/fluor");

    const std::vector<FluorVacancy>* load(int Z);

    // Find vacancy data by EADL shell ID (1=K, 3=L1, 5=L2, 6=L3, 8=M1, …).
    // Returns nullptr if Z has no data or the shell is absent.
    const FluorVacancy* findVacancy(int Z, int eadlShellId);

    // Convenience: total fluorescence yield from a given vacancy (sum of probs).
    double yield(int Z, int vacancyShellIndex);

    const std::string& dir() const { return fDir; }

private:
    std::string fDir;
    std::unordered_map<int, std::vector<FluorVacancy>> fCache;
};

} // namespace g4gamma
