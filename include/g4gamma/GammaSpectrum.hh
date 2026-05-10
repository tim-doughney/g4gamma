// GammaSpectrum.hh -- main public API.
//
// Builds gamma spectra from any IDecayProvider. The provider yields per-decay
// photon emissions already aggregated; GammaSpectrumBuilder walks the chain,
// solves activities (Bateman), and bins emissions weighted by activity.
#pragma once

#include "g4gamma/IsotopeKey.hh"
#include "g4gamma/IDecayProvider.hh"
#include "g4gamma/Geant4Provider.hh"
#include "g4gamma/SandiaProvider.hh"
#include "g4gamma/LaraProvider.hh"
#include "g4gamma/Units.hh"
#include <vector>
#include <string>
#include <memory>

namespace g4gamma {

enum class DataSource {
    Geant4 = 0,
    Sandia = 1,
    Lara   = 2,
};

struct SpectrumOptions {
    DataSource source = DataSource::Geant4;

    bool   includeAnnihilation = true;
    bool   includeXrays        = false; // K-shell X-rays (Geant4 backend) / pass-through (Sandia/LARA)
    bool   fullXrayCascade     = false; // Geant4 backend: K→L→M fluorescence cascade; implies includeXrays
    int    maxChainDepth       = 50;
    double isomerLifetimeThresh = 1.0 * units::ns;
    double levelMatchTolerance  = 1.0 * units::keV;
    int    verbose             = 0;

    std::string geant4Sh;
    std::string sandiaXml;
    std::string laraDir;
};

struct ChainContribution {
    IsotopeKey isotope;
    double     activity;
    double     gammaYield;
    double     meanLife;
};

struct SpectrumResult {
    std::vector<double> binEdges;
    std::vector<double> counts;
    std::vector<ChainContribution> contributions;
    std::string sourceName;
};

class GammaSpectrumBuilder {
public:
    explicit GammaSpectrumBuilder(const SpectrumOptions& opts = {});
    GammaSpectrumBuilder(IDecayProvider& provider, const SpectrumOptions& opts);

    SpectrumResult build(const IsotopeKey& primary,
                         double t,
                         const std::vector<double>& binEdgesAsc);

    const char* sourceName() const;

private:
    SpectrumOptions fOpts;
    std::unique_ptr<IDecayProvider> fOwnedProvider;
    IDecayProvider* fProvider = nullptr;

    void initProviderFromOpts();
};

} // namespace g4gamma
