// IDecayProvider.hh -- abstract source of decay data.
//
// Concrete providers wrap a particular nuclear-data library:
//   * Geant4Provider      -- reads $G4RADIOACTIVEDATA + PhotonEvaporation +
//                            ENSDFSTATE + (optional) G4EMLOW/fluor
//   * SandiaProvider      -- reads sandia.decay.xml from
//                            github.com/sandialabs/SandiaDecay (LGPL-2.1)
//   * (future) LaraProvider, EndfProvider, ...
//
// All providers expose the same high-level structure: for each isotope, the
// list of decay branches, each branch carrying a list of pre-aggregated
// photon emissions (γ + X-rays + annihilation 511) with their per-decay
// intensities. This means GammaSpectrum is provider-agnostic; the cascade
// computation (which Geant4 needs because its data is per-level) is hidden
// inside Geant4Provider.
#pragma once

#include "g4gamma/IsotopeKey.hh"
#include "g4gamma/DecayData.hh"   // for DecayMode + decayModeName
#include <string>
#include <vector>
#include <memory>

namespace g4gamma {

bool decayModeIsBetaPlus(DecayMode m);   // includes B+ family
bool decayModeIsEC(DecayMode m);          // any EC (K/L/M/N/generic)

// One emitted photon (already aggregated -- intensity is per parent decay
// going down THIS branch, not per primary).
enum class EmissionType { Gamma, XRay, AnnihilationPair };
struct Emission {
    EmissionType type;
    double       energy;     // internal units (MeV)
    double       intensity;  // mean photons per parent decay (going down this branch)
};

// One decay channel of a parent.
struct DecayBranch {
    DecayMode               mode;
    double                  branchingRatio;      // relative; sum over branches = 1
    IsotopeKey              daughter;            // (Z', A', M') — used when terminals is empty
    std::vector<Emission>   emissions;
    // If non-empty, overrides `daughter` for activity accounting in ChainBuilder.
    // Each entry is (terminal IsotopeKey, fraction); fractions sum to ~1.
    // The branchingRatio still governs emission weight; only daughter tracking splits.
    // Set by Geant4Provider when the photon-evaporation cascade stops at ≥1 isomeric
    // intermediate level rather than running straight through to the ground state.
    std::vector<std::pair<IsotopeKey,double>> terminals;
};

// Decay info for one parent (Z, A, M).
struct ParentDecayInfo {
    IsotopeKey                  isotope;
    bool                        stable = true;
    double                      meanLife = 0.0;     // internal units (ns); 0 if unknown
    std::vector<DecayBranch>    branches;
};

// Abstract interface.
class IDecayProvider {
public:
    virtual ~IDecayProvider() = default;

    // Return decay info for (Z, A, M). Pointer ownership stays with the
    // provider; safe to keep until the provider is destroyed. nullptr means
    // "no data" -- caller should treat as stable.
    virtual const ParentDecayInfo* get(const IsotopeKey& key) = 0;

    // A short identifier ("geant4", "sandia", ...) for diagnostics.
    virtual const char* name() const = 0;

    // True if this provider includes K X-rays in its emissions (so the user
    // doesn't get them twice if they also have include_xrays=true at the
    // GammaSpectrum level). Geant4 returns false (we add X-rays separately
    // when the user opts in); Sandia returns true (X-rays are baked in).
    virtual bool emissionsIncludeXrays() const = 0;

    // True if this provider's emissions already include 511 keV annihilation
    // pairs from B+ branches. Geant4 false (we add separately); Sandia
    // false too (Sandia's xml doesn't include annihilation -- it gives you
    // the positron energy in <positron> tags but not the resulting 2×511 keV
    // pair, so we still synthesize that ourselves).
    virtual bool emissionsIncludeAnnihilation() const = 0;

    // Convention for the `intensity` field in Emission:
    //   - false (default): "per parent decay going down THIS branch", so the
    //     binner must multiply by branchingRatio. Geant4 and Sandia use this.
    //   - true: "per parent decay" already (sum over branches built in).
    //     LARA's flattened tables use this convention.
    // Affects how the binner combines emission.intensity with branch.branchingRatio.
    virtual bool emissionsArePerDecay() const { return false; }
};

} // namespace g4gamma
