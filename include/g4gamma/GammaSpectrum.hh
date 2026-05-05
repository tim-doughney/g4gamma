// GammaSpectrum.hh -- main public API.
//
// Constructs gamma spectra (binned histograms) from the Geant4 datasets,
// matching the conventions of the user's rdecay01-based simulations:
//   * Discrete gamma lines from PhotonEvaporation cascades following
//     beta/EC/alpha-induced excited daughter states, AND from IT decay of
//     metastable parents.
//   * 511 keV annihilation pairs (2 per beta+ decay) for B+ branches.
//   * Optionally, K-shell X-rays from EC decays and from internal-conversion
//     vacancies in level transitions (xrays=true).
//   * Atomic Auger electrons / radiative cascades beyond K-shell are NOT
//     included; this matches the user's `SetARM(false)` configuration plus a
//     pragmatic restriction to detectable energies.
//
// Output normalisation: histogram contents are *gammas per primary decay* of
// the input isotope, summed over all chain members weighted by their activity
// relative to the primary at the requested time. At secular equilibrium
// (time < 0) this reproduces what rdecay01 produces with
// thresholdForVeryLongDecayTime set to a large value (i.e. all chain decays
// happen and BRs stack along all paths from primary to each chain member).
#pragma once

#include "g4gamma/IsotopeKey.hh"
#include "g4gamma/DataPath.hh"
#include "g4gamma/DecayData.hh"
#include "g4gamma/PhotonEvap.hh"
#include "g4gamma/FluorData.hh"
#include "g4gamma/EnsdfState.hh"
#include "g4gamma/ChainBuilder.hh"
#include "g4gamma/Bateman.hh"
#include "g4gamma/Units.hh"
#include <vector>
#include <string>
#include <memory>
#include <map>
#include <tuple>

namespace g4gamma {

struct SpectrumOptions {
    bool   includeAnnihilation = true;   // 2x 511 keV per B+ decay
    bool   includeXrays        = false;  // K X-rays from EC + IC
    int    maxChainDepth       = 50;     // hard cap on chain size
    double isomerLifetimeThresh = 1.0 * units::ns;  // matches your physics list
    std::string geant4Sh = "";           // optional; if empty, env vars + defaults
    // Tolerance for matching daughter excitation to PhotonEvap levels, in
    // internal units. 1 keV matches Geant4's level rounding in the decay files.
    double levelMatchTolerance = 1.0 * units::keV;
    // Diagnostic: when true, prints data-path resolution, files loaded,
    // chain members discovered, and per-channel cascade firing to stderr.
    int    verbose = 0;
};

struct ChainContribution {
    IsotopeKey isotope;
    double     activity;       // relative to primary (=1 for primary at t=0)
    double     gammaYield;     // expected total gammas per decay of THIS nuclide
    double     meanLife;       // internal units
};

struct SpectrumResult {
    std::vector<double> binEdges;     // N+1 entries, internal units (MeV)
    std::vector<double> counts;       // N entries: gammas per primary decay
    // Diagnostic: per-chain-member breakdown.
    std::vector<ChainContribution> contributions;
};

class GammaSpectrumBuilder {
public:
    GammaSpectrumBuilder(const SpectrumOptions& opts = {});

    // Build a spectrum.
    //   primary       : the primary nuclide, with activity 1 (Bq) at t=0.
    //   t             : time elapsed since the activity reference. Negative
    //                   value (or +inf) means secular equilibrium.
    //                   Internal units (i.e. multiply by g4gamma::units::s,
    //                   year, etc.).
    //   binEdgesAsc   : ascending bin edges in internal units (MeV); N+1 long.
    SpectrumResult build(const IsotopeKey& primary,
                         double t,
                         const std::vector<double>& binEdgesAsc);

    // Override data directories at runtime (otherwise resolved from env/sh).
    void setRadDir(std::string dir)         { fRadDir = std::move(dir); }
    void setEvapDir(std::string dir)        { fEvapDir = std::move(dir); }
    void setLeDataDir(std::string dir)      { fLeDataDir = std::move(dir); }
    void setEnsdfStateDir(std::string dir)  { fEnsdfStateDir = std::move(dir); }

private:
    SpectrumOptions fOpts;
    std::string fRadDir, fEvapDir, fLeDataDir, fEnsdfStateDir;

    // For each starting level of (Z,A), the list of (energy, weight) photons
    // emitted on full de-excitation to ground. Cached.
    struct CascadeEmission {
        std::vector<double> energies;
        std::vector<double> weights;     // mean number of gammas at this energy per cascade
        // K-shell vacancy expected count per cascade (used when xrays=true).
        double kVacancyExpected = 0.0;
    };
    std::map<std::tuple<int,int,int>, CascadeEmission> fCascadeCache;

    const CascadeEmission& getCascade(int Z, int A, int startLevel,
                                      PhotonEvapData& evap);

    // Resolve data dirs lazily.
    void ensureDirs();
};

} // namespace g4gamma
