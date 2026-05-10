#include "g4gamma/GammaSpectrum.hh"
#include "g4gamma/ChainBuilder.hh"
#include "g4gamma/Bateman.hh"
#include <algorithm>
#include <cmath>
#include <iostream>
#include <stdexcept>

namespace g4gamma {

GammaSpectrumBuilder::GammaSpectrumBuilder(const SpectrumOptions& opts)
    : fOpts(opts) {
    initProviderFromOpts();
}

GammaSpectrumBuilder::GammaSpectrumBuilder(IDecayProvider& p, const SpectrumOptions& opts)
    : fOpts(opts), fProvider(&p) {}

void GammaSpectrumBuilder::initProviderFromOpts() {
    if (fOpts.source == DataSource::Geant4) {
        Geant4ProviderOptions go;
        go.geant4Sh        = fOpts.geant4Sh;
        go.includeXrays    = fOpts.includeXrays;
        go.fullXrayCascade = fOpts.fullXrayCascade;
        go.isomerThreshold = fOpts.isomerLifetimeThresh;
        go.levelTolerance  = fOpts.levelMatchTolerance;
        fOwnedProvider = std::make_unique<Geant4Provider>(go);
    } else if (fOpts.source == DataSource::Sandia) {
        fOwnedProvider = std::make_unique<SandiaProvider>(fOpts.sandiaXml);
    } else if (fOpts.source == DataSource::Lara) {
        fOwnedProvider = std::make_unique<LaraProvider>(fOpts.laraDir);
    } else {
        throw std::runtime_error("GammaSpectrumBuilder: unknown DataSource");
    }
    fProvider = fOwnedProvider.get();
    if (fOpts.verbose >= 1) {
        std::cerr << "[g4gamma] data source: " << fProvider->name() << "\n";
    }
}

const char* GammaSpectrumBuilder::sourceName() const {
    return fProvider ? fProvider->name() : "(none)";
}

static void binInto(std::vector<double>& counts,
                    const std::vector<double>& edges,
                    double E, double w) {
    if (w == 0.0) return;
    if (E < edges.front() || E >= edges.back()) return;
    int lo = 0, hi = static_cast<int>(edges.size()) - 1;
    while (lo + 1 < hi) {
        int mid = (lo + hi) / 2;
        if (edges[mid] <= E) lo = mid; else hi = mid;
    }
    counts[lo] += w;
}

SpectrumResult GammaSpectrumBuilder::build(const IsotopeKey& primary,
                                            double t,
                                            const std::vector<double>& edges) {
    if (edges.size() < 2) throw std::invalid_argument("Need at least 2 bin edges");
    for (size_t i = 1; i < edges.size(); ++i) {
        if (edges[i] <= edges[i-1])
            throw std::invalid_argument("Bin edges must be strictly ascending");
    }
    if (!fProvider) throw std::runtime_error("No data provider");

    if (fOpts.verbose >= 1) {
        std::cerr << "[g4gamma] building chain for " << primary.str()
                  << " using provider '" << fProvider->name() << "'\n";
    }

    ChainBuilder cb(*fProvider);
    cb.build(primary, fOpts.maxChainDepth);
    const auto& chain = cb.nodes();

    if (fOpts.verbose >= 1) {
        std::cerr << "[g4gamma] chain has " << chain.size() << " members\n";
        for (size_t i = 0; i < chain.size(); ++i) {
            const auto& n = chain[i];
            std::cerr << "  [" << i << "] " << n.isotope.str()
                      << " meanLife=" << n.meanLife / units::s << " s"
                      << " stable=" << n.stable
                      << " edges=" << n.edges.size() << "\n";
        }
    }

    Bateman bat(chain);
    std::vector<double> A;
    if (t < 0.0 || std::isinf(t)) A = bat.solveSecularEq();
    else                          A = bat.solve(t);

    SpectrumResult out;
    out.sourceName = fProvider->name();
    out.binEdges = edges;
    out.counts.assign(edges.size() - 1, 0.0);
    out.contributions.reserve(chain.size());

    bool providerHasXrays        = fProvider->emissionsIncludeXrays();
    bool providerHasAnnihilation = fProvider->emissionsIncludeAnnihilation();
    bool providerPerDecay        = fProvider->emissionsArePerDecay();

    for (size_t i = 0; i < chain.size(); ++i) {
        const auto& node = chain[i];
        ChainContribution cc;
        cc.isotope    = node.isotope;
        cc.activity   = A[i];
        cc.meanLife   = node.meanLife;
        cc.gammaYield = 0.0;
        if (node.stable || A[i] <= 0.0) {
            out.contributions.push_back(cc);
            continue;
        }
        const ParentDecayInfo* info = fProvider->get(node.isotope);
        if (!info) {
            out.contributions.push_back(cc);
            continue;
        }

        for (const auto& branch : info->branches) {
            for (const auto& em : branch.emissions) {
                bool keep = true;
                if (em.type == EmissionType::XRay &&
                    !fOpts.includeXrays && !fOpts.fullXrayCascade) keep = false;
                if (!keep) continue;
                // Annihilation skip if the user disabled it
                if (em.type == EmissionType::AnnihilationPair &&
                    !fOpts.includeAnnihilation) continue;
                double w = em.intensity;
                if (!providerPerDecay) w *= branch.branchingRatio;
                cc.gammaYield += w;
                binInto(out.counts, edges, em.energy, w * A[i]);
            }
            // Synthesise annihilation pairs from B+ if the provider doesn't
            // already supply them and the user wants them.
            if (fOpts.includeAnnihilation && !providerHasAnnihilation &&
                decayModeIsBetaPlus(branch.mode)) {
                double w = 2.0 * branch.branchingRatio;
                cc.gammaYield += w;
                binInto(out.counts, edges, units::electron_mass_c2, w * A[i]);
            }
            (void)providerHasXrays;
        }
        out.contributions.push_back(cc);
    }
    return out;
}

} // namespace g4gamma
