#include "g4gamma/GammaSpectrum.hh"
#include <algorithm>
#include <cmath>
#include <stdexcept>
#include <iostream>
#include <limits>

namespace g4gamma {

GammaSpectrumBuilder::GammaSpectrumBuilder(const SpectrumOptions& opts)
    : fOpts(opts) {}

void GammaSpectrumBuilder::ensureDirs() {
    if (fRadDir.empty())    fRadDir    = DataPath::radioactiveDecayDir(fOpts.geant4Sh);
    if (fEvapDir.empty())   fEvapDir   = DataPath::photonEvaporationDir(fOpts.geant4Sh);
    if (fEnsdfStateDir.empty()) {
        // Optional: don't throw if not found.
        try { fEnsdfStateDir = DataPath::ensdfStateDir(fOpts.geant4Sh); }
        catch (...) {}
    }
    if (fOpts.includeXrays && fLeDataDir.empty()) {
        try { fLeDataDir = DataPath::lowEnergyDir(fOpts.geant4Sh); }
        catch (const std::exception& e) {
            // X-rays requested but G4LEDATA not found -> keep going without
            // X-rays but warn loudly.
            std::cerr << "[g4gamma] WARNING: X-rays requested but " << e.what()
                      << " -> continuing without X-rays.\n";
            fOpts.includeXrays = false;
        }
    }
    if (fOpts.verbose >= 1) {
        std::cerr << "[g4gamma] RadioactiveDecay dir : " << fRadDir << "\n"
                  << "[g4gamma] PhotonEvaporation dir: " << fEvapDir << "\n"
                  << "[g4gamma] G4ENSDFSTATEDATA dir : "
                  << (fEnsdfStateDir.empty() ? "<not found>" : fEnsdfStateDir) << "\n";
        if (fOpts.includeXrays)
            std::cerr << "[g4gamma] G4LEDATA dir         : " << fLeDataDir << "\n";
    }
}

// Compute the cascade emission spectrum for a given starting level of (Z,A).
// Result is the *expected* gammas (and K-vacancy count) per de-excitation
// from level `startLevel` down to the ground state.
//
// Method: iterate from `startLevel` downwards. Maintain a probability array
// p[i] = probability of being at level i during the cascade (starts at p[startLevel]=1).
// For each level i (descending), distribute its probability among its
// transitions. For each transition j with selectionProb s_j and gammaEmitProb g_j:
//   * mean gammas at energy E_j contributed: p[i] * s_j * g_j
//   * mean K-vacancies contributed: p[i] * s_j * (1-g_j) * iccWeight[K]/sum(iccWeights)
//   * probability flowing to level lowerIndex_j: p[i] * s_j (independent of g_j;
//     the level will reach the lower state regardless of whether photon or IC electron)
// Continue until i = 0 (ground state) or no more probability left.
const GammaSpectrumBuilder::CascadeEmission&
GammaSpectrumBuilder::getCascade(int Z, int A, int startLevel,
                                  PhotonEvapData& evap) {
    auto key = std::make_tuple(Z, A, startLevel);
    auto it = fCascadeCache.find(key);
    if (it != fCascadeCache.end()) return it->second;

    CascadeEmission em;
    const auto* lvls = evap.load(Z, A);
    if (!lvls || lvls->empty() || startLevel <= 0 ||
        startLevel >= static_cast<int>(lvls->size())) {
        // No cascade -- either ground state or no level data.
        fCascadeCache[key] = em;
        return fCascadeCache[key];
    }

    int N = lvls->size();
    std::vector<double> p(N, 0.0);
    p[startLevel] = 1.0;

    // Process levels in descending order
    for (int i = N - 1; i > 0; --i) {
        if (p[i] <= 0.0) continue;
        const auto& lvl = (*lvls)[i];
        if (lvl.transitions.empty()) {
            // dead-end (no transitions defined) -- treat as IT-trapped / floating,
            // the probability is lost (matches Geant4's behaviour at the very top).
            continue;
        }
        for (const auto& tr : lvl.transitions) {
            if (tr.selectionProb <= 0.0) continue;
            if (tr.lowerIndex < 0 || tr.lowerIndex >= N) continue;
            double pSel = p[i] * tr.selectionProb;
            // Photon
            double gammaCount = pSel * tr.gammaEmitProb;
            if (gammaCount > 0.0) {
                em.energies.push_back(tr.gammaEnergy);
                em.weights.push_back(gammaCount);
            }
            // IC electron -> shell vacancy. K-vacancy is iccWeights[0].
            double icCount = pSel * (1.0 - tr.gammaEmitProb);
            if (icCount > 0.0 && tr.hasICC) {
                double sumW = 0.0;
                for (int k = 0; k < 10; ++k) sumW += tr.iccWeights[k];
                if (sumW > 0.0) {
                    em.kVacancyExpected += icCount * tr.iccWeights[0] / sumW;
                }
            }
            // probability flows to lower level regardless of photon vs IC
            p[tr.lowerIndex] += pSel;
        }
    }

    fCascadeCache[key] = std::move(em);
    return fCascadeCache[key];
}

// Helper: bin a list of (E, w) photons into the user histogram.
static void binInto(std::vector<double>& counts,
                    const std::vector<double>& edges,
                    double E, double w) {
    if (w == 0.0) return;
    if (E < edges.front()) return;
    if (E >= edges.back()) return;
    // Binary search for upper bound
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
    if (edges.size() < 2) {
        throw std::invalid_argument("Need at least 2 bin edges");
    }
    for (size_t i = 1; i < edges.size(); ++i) {
        if (edges[i] <= edges[i-1]) {
            throw std::invalid_argument("Bin edges must be strictly ascending");
        }
    }

    ensureDirs();

    DecayDataLoader  decay(fRadDir);
    PhotonEvapData   evap(fEvapDir);
    EnsdfStateLoader ensdf(fEnsdfStateDir);
    std::unique_ptr<FluorDataLoader> fluor;
    if (fOpts.includeXrays) {
        // Use default /fluor (matching Geant4 default; user can override later).
        fluor = std::make_unique<FluorDataLoader>(fLeDataDir, "/fluor");
    }

    // Build chain
    ChainBuilder cb(decay, evap, fOpts.isomerLifetimeThresh);
    cb.build(primary, fOpts.maxChainDepth);
    auto& chainMutable = const_cast<std::vector<ChainNode>&>(cb.nodes());

    // Patch in mean lives from ENSDFSTATE if available. The half-life column
    // in the RadioactiveDecay zZ.aA file is a placeholder/dummy and using it
    // for finite-time Bateman gives wrong answers (often A=0 across the
    // board). ENSDFSTATE.dat is the canonical source.
    if (ensdf.ready()) {
        for (auto& n : chainMutable) {
            if (n.stable) continue;
            double ml = ensdf.meanLife(n.isotope, fOpts.isomerLifetimeThresh);
            if (ml > 0.0) n.meanLife = ml;
        }
    }
    const auto& chain = cb.nodes();

    if (fOpts.verbose >= 1) {
        std::cerr << "[g4gamma] Chain rooted at " << primary.str() << ": "
                  << chain.size() << " members\n";
        for (size_t i = 0; i < chain.size(); ++i) {
            const auto& n = chain[i];
            std::cerr << "  [" << i << "] " << n.isotope.str()
                      << " meanLife=" << n.meanLife / units::s << " s"
                      << " stable=" << n.stable
                      << " channels=" << n.edges.size() << "\n";
        }
        // Warn about isotopes with no decay file
        for (const auto& n : chain) {
            const DecayParent* dp = decay.get(n.isotope);
            if (!dp && !n.stable) {
                std::cerr << "[g4gamma] WARNING: " << n.isotope.str()
                          << " not stable but no decay-data file found at "
                          << fRadDir << "/z" << n.isotope.Z << ".a"
                          << n.isotope.A << "\n";
            }
            if (dp && dp->channels.empty() && !n.stable) {
                std::cerr << "[g4gamma] WARNING: " << n.isotope.str()
                          << " has decay file but zero channels\n";
            }
        }
    }

    // Sanity check on the root: if it has no decay channels, the result will
    // be all zero and that's almost certainly not what the user wants.
    if (!chain.empty()) {
        const DecayParent* rootDp = decay.get(chain[0].isotope);
        if (!rootDp || rootDp->channels.empty()) {
            std::cerr << "[g4gamma] WARNING: primary " << primary.str()
                      << " has no decay channels -- output will be zero.\n"
                      << "  Looked in: " << fRadDir << "/z"
                      << primary.Z << ".a" << primary.A << "\n"
                      << "  File exists? "
                      << (rootDp ? "yes (but empty)" : "no") << "\n"
                      << "  Set verbose=1 in SpectrumOptions for more detail.\n";
        }
    }

    // Solve activities
    Bateman bat(chain);
    std::vector<double> A;
    if (t < 0.0 || std::isinf(t)) {
        A = bat.solveSecularEq();
    } else {
        A = bat.solve(t);
    }

    SpectrumResult out;
    out.binEdges = edges;
    out.counts.assign(edges.size() - 1, 0.0);
    out.contributions.reserve(chain.size());

    // For each unstable chain member: compute gammas per decay and accumulate
    // weighted by activity.
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
        const DecayParent* dp = decay.get(node.isotope);
        if (!dp) {
            out.contributions.push_back(cc);
            continue;
        }

        // Accumulate gammas per decay of this node
        // Format: vector of (energy, weight)
        std::vector<std::pair<double,double>> nodeGammas;
        // K-vacancies per decay of this node (used only if X-rays enabled)
        double nodeKVacancy = 0.0;

        // Iterate the channels and use the corresponding chain edge to know
        // whether the daughter was promoted to an isomer (M>=1). When the
        // daughter is a tracked isomer, the gamma cascade is deferred to the
        // isomer's own IT decay; we do NOT fire the cascade here, otherwise
        // we'd double-count the gammas.
        if (dp->channels.size() != node.edges.size()) {
            // Defensive: the channels and edges should be in 1:1 correspondence
            // (via parentChannelIdx). Build an index map.
        }
        // Build an edge lookup by parentChannelIdx
        std::vector<int> ch2edge(dp->channels.size(), -1);
        for (size_t ei = 0; ei < node.edges.size(); ++ei) {
            int ci = node.edges[ei].parentChannelIdx;
            if (ci >= 0 && ci < static_cast<int>(ch2edge.size())) {
                ch2edge[ci] = static_cast<int>(ei);
            }
        }

        for (size_t ci = 0; ci < dp->channels.size(); ++ci) {
            const auto& ch = dp->channels[ci];
            // Resolve the daughter's chain entry (which may have promoted M).
            int eIdx = ch2edge[ci];
            IsotopeKey resolvedDaughter = ch.daughter;
            if (eIdx >= 0) {
                int dIdx = node.edges[eIdx].daughterIndex;
                resolvedDaughter = chain[dIdx].isotope;
            }
            // If the daughter is an isomer that we are tracking separately,
            // defer the cascade to the isomer's own decay block.
            bool deferToIsomer = (ch.mode != DecayMode::IT) &&
                                  (resolvedDaughter.M >= 1) &&
                                  (ch.daughterExcitation > 0.0);

            // For IT decay: cascade starts from this nuclide's own excited level
            // (the parent's parentExcitation in the decay file). The "daughter"
            // is the same nuclide at ground state.
            // For other modes: cascade starts from the daughter's excited level
            // ch.daughterExcitation in PhotonEvap data of (Z',A').
            // 
            // Step 1: determine cascade
            int cZ, cA;
            double startE;
            if (ch.mode == DecayMode::IT) {
                cZ = node.isotope.Z;
                cA = node.isotope.A;
                startE = dp->parentExcitation;
            } else {
                cZ = ch.daughter.Z;
                cA = ch.daughter.A;
                startE = ch.daughterExcitation;
            }
            // Locate the starting level
            int startLevel = evap.findLevel(cZ, cA, startE, fOpts.levelMatchTolerance);
            if (startLevel < 0) startLevel = 0;
            // Only fire the cascade and emit annihilation/EC X-rays if we are
            // NOT deferring to a tracked isomer.
            if (!deferToIsomer) {
                const auto& em = getCascade(cZ, cA, startLevel, evap);
                for (size_t ki = 0; ki < em.energies.size(); ++ki) {
                    nodeGammas.emplace_back(em.energies[ki], em.weights[ki] * ch.branchingRatio);
                }
                // K-vacancies from internal conversion in the cascade
                nodeKVacancy += em.kVacancyExpected * ch.branchingRatio;
            }

            // Annihilation gammas (always fire on B+; the positron is born from
            // the parent decay regardless of whether the daughter is an isomer).
            if (fOpts.includeAnnihilation && ch.mode == DecayMode::BetaPlus) {
                // 2 photons of 511 keV per B+ decay
                nodeGammas.emplace_back(units::electron_mass_c2, 2.0 * ch.branchingRatio);
            }
            // K-vacancies from EC decay (the EC happens here regardless of
            // daughter level state).
            if (fOpts.includeXrays && ch.mode == DecayMode::KshellEC) {
                nodeKVacancy += ch.branchingRatio;  // 1 K-vacancy per K-EC decay
            }
            // Note: LshellEC, MshellEC, NshellEC create L/M/N vacancies which we
            // don't currently turn into X-rays (their energies are typically too
            // low to be relevant for NaI/LaBr3, and the cascade through K isn't
            // happening since the K shell isn't ionised in those modes).
        }

        // Convert K-vacancies to K X-ray emission
        if (fOpts.includeXrays && nodeKVacancy > 0.0 && fluor) {
            // For EC decay, the X-rays come from the daughter (Z' = Z-1).
            // For IC, they come from the parent (Z) -- the conversion happened
            // in the parent. To keep things simple we use the daughter Z for EC
            // and parent Z for IC. We mixed these into nodeKVacancy without
            // distinguishing -- this is a small approximation. Refine if needed.
            //
            // We approximate: use the parent's Z for IC vacancies and the
            // daughter's Z (= parent's Z - 1) for K-EC. Since both contribute,
            // and the spectra are similar (X-ray energy scales ~ Z^2), we use
            // the *daughter* Z if there's any EC channel, else the parent Z.
            int Zfluor = node.isotope.Z;
            for (const auto& ch : dp->channels) {
                if (ch.mode == DecayMode::KshellEC ||
                    ch.mode == DecayMode::BetaPlus) {
                    Zfluor = ch.daughter.Z;  // X-ray comes from daughter atom
                    break;
                }
            }
            const auto* vacs = fluor->load(Zfluor);
            if (vacs && !vacs->empty()) {
                // K vacancy is index 0 in EADL ordering.
                const auto& kvac = (*vacs)[0];
                for (const auto& tr : kvac.transitions) {
                    nodeGammas.emplace_back(tr.energy,
                                             nodeKVacancy * tr.prob);
                }
            }
        }

        // Accumulate yield total
        for (const auto& g : nodeGammas) cc.gammaYield += g.second;

        // Bin into histogram, weighted by activity
        for (const auto& g : nodeGammas) {
            binInto(out.counts, edges, g.first, g.second * A[i]);
        }
        out.contributions.push_back(cc);
    }

    return out;
}

} // namespace g4gamma
