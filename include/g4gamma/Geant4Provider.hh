// Geant4Provider.hh -- IDecayProvider backed by the standard Geant4 data
// directories (RadioactiveDecay, PhotonEvaporation, ENSDFSTATE, optionally
// G4LEDATA/fluor).
//
// This wraps the legacy DecayDataLoader / PhotonEvapData / EnsdfStateLoader /
// FluorDataLoader stack and the cascade-firing logic that was previously
// inside GammaSpectrumBuilder. The cascade computation moves here so that
// GammaSpectrumBuilder doesn't need to know whether a provider is per-level
// (Geant4) or per-decay (Sandia, LARA).
#pragma once

#include "g4gamma/IDecayProvider.hh"
#include "g4gamma/DecayData.hh"
#include "g4gamma/PhotonEvap.hh"
#include "g4gamma/EnsdfState.hh"
#include "g4gamma/FluorData.hh"
#include "g4gamma/AugerData.hh"
#include "g4gamma/Units.hh"
#include <string>
#include <unordered_map>
#include <map>
#include <tuple>
#include <memory>

namespace g4gamma {

struct Geant4ProviderOptions {
    std::string radDir;          // empty -> auto-resolve
    std::string evapDir;         // empty -> auto-resolve
    std::string ensdfDir;        // empty -> auto-resolve (optional)
    std::string ledataDir;       // empty -> auto-resolve when xrays requested
    std::string fluorSubdir = "/fluor";
    std::string geant4Sh;        // optional path to geant4.sh
    bool        includeXrays    = false; // K-shell X-rays from IC/EC
    bool        fullXrayCascade = false; // K→L→M cascade (implies includeXrays)
    double      isomerThreshold = 1.0 * units::ns;
    double      levelTolerance  = 1.0 * units::keV;
};

class Geant4Provider : public IDecayProvider {
public:
    explicit Geant4Provider(const Geant4ProviderOptions& opts = {});

    const ParentDecayInfo* get(const IsotopeKey& key) override;
    const char* name() const override { return "geant4"; }
    bool emissionsIncludeXrays() const override {
        return fOpts.includeXrays || fOpts.fullXrayCascade;
    }
    bool emissionsIncludeAnnihilation() const override { return false; }

    // Diagnostic accessors
    const std::string& radDir()      const { return fOpts.radDir; }
    const std::string& evapDir()     const { return fOpts.evapDir; }
    const std::string& ensdfDir()    const { return fOpts.ensdfDir; }
    const std::string& ledataDir()   const { return fOpts.ledataDir; }

private:
    Geant4ProviderOptions fOpts;
    std::unique_ptr<DecayDataLoader>    fDecay;
    std::unique_ptr<PhotonEvapData>     fEvap;
    std::unique_ptr<EnsdfStateLoader>   fEnsdf;
    std::unique_ptr<FluorDataLoader>    fFluor;
    std::unique_ptr<AugerDataLoader>    fAuger;
    std::unordered_map<IsotopeKey, ParentDecayInfo> fCache;

    // Build (cache + return) ParentDecayInfo for `key`.
    const ParentDecayInfo* compute(const IsotopeKey& key);

    // Resolve a decay channel's daughter into a (Z',A',M') key, using
    // ENSDFSTATE when available. Identical to the old ChainBuilder logic
    // but now lives here.
    IsotopeKey resolveDaughter(const DecayChannel& ch, int parentZ, int parentA);

    // Find the right DecayParent record for an IsotopeKey, accounting for
    // missing ground-state P blocks in some daughter files (e.g. Ba-137).
    const DecayParent* parentFor(const IsotopeKey& key);

    // Cascade emissions: starting from the given level of (Z,A), accumulate
    // gammas, per-shell IC vacancy yields, and terminal-level probabilities.
    // terminalProb records where probability stopped: each entry (levelIdx, prob)
    // means `prob` fraction of cascades ended at that level (level 0 = ground
    // state; positive indices = long-lived isomers with τ > isomerThreshold).
    // The sum of terminalProb entries equals the fraction of probability that
    // did not fall into dead-end levels.
    struct Cascade {
        // Number of ICC shell slots: K, L1, L2, L3, M1, M2, M3, M4, M5, N+
        static constexpr int N_ICC = 10;
        std::vector<double> energies;       // internal units (MeV)
        std::vector<double> intensities;    // mean photons per cascade firing
        std::array<double, N_ICC> shellVacancy = {};  // mean IC vacancies per shell
        // (levelIndex, probability) pairs for where the cascade stopped.
        // Level 0 = ground state; positive indices = long-lived intermediate isomers.
        // Empty when startLevel=0 or the cascade runs straight to the ground state.
        std::vector<std::pair<int,double>> terminalProb;
    };
    std::map<std::tuple<int,int,int>, Cascade> fCascadeCache;
    const Cascade& getCascade(int Z, int A, int startLevel);

    // Map a PhotonEvaporation level index to the corresponding IsotopeKey,
    // using ENSDFSTATE to assign the isomer quantum number M.
    IsotopeKey levelIndexToKey(int Z, int A, int levelIdx);

    // EADL shell ID (as stored in fl-tr-pr files) → ICC index (0=K … 9=N+).
    // Returns -1 for outer shells not tracked in the 10-slot ICC array.
    static int eadlToICC(int eadlId);

    // Append X-ray emissions to `out` given per-shell vacancy counts.
    // In K-shell mode only K vacancies are converted; in full-cascade mode
    // secondary vacancies from each fluorescence transition are propagated
    // through L/M shells. The shellVac array is taken by value (modified internally).
    void appendXRays(int Zfluor,
                     std::array<double, Cascade::N_ICC> shellVac,
                     std::vector<Emission>& out);
};

} // namespace g4gamma
