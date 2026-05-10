// PhotonEvap.hh -- parses Geant4 PhotonEvaporation/zZ.aA files
// to obtain the level scheme and per-level gamma transitions, then computes
// gamma yields per de-excitation cascade starting from any level.
//
// File format (per G4LevelReader::LevelManager):
//   For each level i (i = 0, 1, 2, ...):
//     <i_int> <floatFlag>
//     <E_keV> <T_half_s> <spin> <ntrans>   (T_half_s is half-life in seconds)
//     repeated ntrans times:
//       <i2> <Etrans_keV> <fProb> <multipolarity> <ratio> <fAlpha>
//       if fAlpha > 0: 10 more numbers (K,L1,L2,L3,M1..,N..) = ICC weights
//
// Where:
//   i2 = lower level index that this transition goes to (i2 < i)
//   fProb = relative probability of selecting this transition
//   fAlpha = total internal conversion coefficient (gamma-emission prob = 1/(1+fAlpha))
//
// The selection probability for transition j (from level i) is:
//   P_select(j) = (1+alpha_j) * fProb_j / sum_k [(1+alpha_k) * fProb_k]
// And given that transition j was selected, a photon (gamma) is emitted with
// probability 1/(1+alpha_j) and an internal-conversion electron with prob
// alpha_j/(1+alpha_j). The shell distribution of the IC vacancy follows the
// 10-number weight vector.
//
// We compute, for any starting level, the *expected number of gammas* per
// de-excitation across the entire cascade, binned to the user's grid. We also
// compute the expected number of K/L/M-shell vacancies per de-excitation,
// which feeds into the X-ray code.
#pragma once

#include "g4gamma/IsotopeKey.hh"
#include <string>
#include <vector>
#include <unordered_map>
#include <array>

namespace g4gamma {

struct LevelTransition {
    int    lowerIndex;          // i2
    double gammaEnergy;         // internal units (MeV) -- explicit transition energy from file
    double selectionProb;       // P_select(j) -- already normalized so sum = 1
    double gammaEmitProb;       // = 1/(1+alpha)
    // 10 ICC weights (K, L1, L2, L3, M1, M2, M3, M4, M5, N+...) -- raw values from file.
    // Only valid if alpha > 0; sum is meaningful for distributing IC vacancies.
    std::array<float, 10> iccWeights;
    bool hasICC;                // true if alpha > 0 (ICC values were read)
};

struct NuclearLevel {
    double energy;                          // MeV
    double meanLifeTime;                    // ns (internal units); DBL_MAX if "stable" wrt EM
    std::vector<LevelTransition> transitions;
};

class PhotonEvapData {
public:
    explicit PhotonEvapData(std::string evapDir);

    // Load level scheme for (Z,A). Returns nullptr if no file. Result is a
    // sorted-by-energy vector of levels (level 0 = ground).
    const std::vector<NuclearLevel>* load(int Z, int A);

    // Find level index whose energy is closest to `excitation` (within tolerance).
    // Returns -1 if no level scheme. Returns 0 (ground state) if excitation < tol.
    int findLevel(int Z, int A, double excitation_internal,
                  double tol_internal /*=1 keV by default in Units.hh*/);

    const std::string& dir() const { return fDir; }

private:
    std::string fDir;
    std::unordered_map<int, std::vector<NuclearLevel>> fCache; // 1000*A + Z
};

} // namespace g4gamma
