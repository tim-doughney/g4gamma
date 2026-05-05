#include "g4gamma/ChainBuilder.hh"
#include "g4gamma/Units.hh"
#include <queue>
#include <cmath>
#include <algorithm>
#include <limits>

namespace g4gamma {

ChainBuilder::ChainBuilder(DecayDataLoader& d, PhotonEvapData& e,
                            EnsdfStateLoader* ensdf, double thr)
    : fDec(d), fEvap(e), fEnsdf(ensdf), fIsomerThreshold(thr) {}

int ChainBuilder::indexOf(const IsotopeKey& k) const {
    auto it = fIndex.find(k);
    return (it == fIndex.end()) ? -1 : it->second;
}

int ChainBuilder::addNode(const IsotopeKey& key) {
    auto it = fIndex.find(key);
    if (it != fIndex.end()) return it->second;
    ChainNode n;
    n.isotope = key;
    n.index   = static_cast<int>(fNodes.size());
    const DecayParent* dp = parentFor(key);
    if (dp) {
        n.meanLife = dp->meanLife;
        n.stable   = dp->isStable;
    } else {
        n.meanLife = 0.0;
        n.stable   = true;
    }
    fNodes.push_back(std::move(n));
    fIndex[key] = fNodes.back().index;
    return fNodes.back().index;
}

const DecayParent* ChainBuilder::parentFor(const IsotopeKey& key) {
    const auto* parents = fDec.load(key.Z, key.A);
    if (!parents || parents->empty()) return nullptr;

    // Determine the target excitation for this M from ENSDFSTATE if available.
    double targetExc = -1.0;
    if (fEnsdf) {
        const auto* lvls = fEnsdf->load(key.Z, key.A);
        if (lvls) {
            int M = 0;
            bool seenGround = false;
            for (const auto& lvl : *lvls) {
                int thisM;
                if (lvl.excitation == 0.0 && !seenGround) {
                    thisM = 0;
                    seenGround = true;
                } else if (lvl.excitation > 0.0 && lvl.meanLife > fIsomerThreshold) {
                    ++M;
                    thisM = M;
                } else {
                    thisM = -1;
                }
                if (thisM == key.M) {
                    targetExc = lvl.excitation;
                    break;
                }
            }
        }
    }

    if (targetExc >= 0.0) {
        // Find the P block whose excitation matches targetExc.
        const DecayParent* best = nullptr;
        double bestDiff = std::numeric_limits<double>::infinity();
        for (const auto& p : *parents) {
            double diff = std::abs(p.parentExcitation - targetExc);
            if (diff < bestDiff) { bestDiff = diff; best = &p; }
        }
        if (bestDiff > 1.0 * units::keV) return nullptr;
        return best;
    }

    // No ENSDFSTATE info -- fall back to file-order assumption.
    if (key.M >= 0 && key.M < static_cast<int>(parents->size())) {
        return &(*parents)[key.M];
    }
    return nullptr;
}

bool ChainBuilder::resolveDaughter(const DecayChannel& ch, IsotopeKey& outKey) {
    int Zp = ch.daughter.Z;
    int Ap = ch.daughter.A;
    double exc = ch.daughterExcitation;

    // Special case: IT means the "daughter" is the same nuclide in a lower
    // (typically ground) state. Geant4 stores no daughterExcitation for IT
    // (the short line carries no per-channel info). We resolve to (Z,A,M=0).
    if (ch.mode == DecayMode::IT) {
        outKey = IsotopeKey{Zp, Ap, 0};
        return true;
    }

    // Use ENSDFSTATE if available.
    int M = 0;
    if (fEnsdf) {
        int m = fEnsdf->excitationToM(Zp, Ap, exc, fIsomerThreshold,
                                       1.0 * units::keV);
        if (m >= 0) M = m;
    } else {
        // Fall back to file-order match in the daughter's decay file.
        if (const auto* parents = fDec.load(Zp, Ap)) {
            if (!parents->empty()) {
                double bestDiff = std::numeric_limits<double>::infinity();
                int    bestM    = 0;
                for (size_t i = 0; i < parents->size(); ++i) {
                    double d = std::abs((*parents)[i].parentExcitation - exc);
                    if (d < bestDiff) { bestDiff = d; bestM = static_cast<int>(i); }
                }
                if (bestM > 0 && bestDiff < 1.0 * units::keV) M = bestM;
            }
        }
    }

    outKey = IsotopeKey{Zp, Ap, M};
    return true;
}

void ChainBuilder::build(const IsotopeKey& root, int maxDepth) {
    fNodes.clear();
    fIndex.clear();

    addNode(root);

    std::queue<int> q;
    q.push(0);

    int processed = 0;
    while (!q.empty() && processed < maxDepth * 1000) {
        int idx = q.front();
        q.pop();
        ++processed;

        ChainNode& n = fNodes[idx];
        if (n.stable) continue;

        const DecayParent* dp = parentFor(n.isotope);
        if (!dp) { n.stable = true; continue; }

        for (size_t ci = 0; ci < dp->channels.size(); ++ci) {
            const auto& ch = dp->channels[ci];
            IsotopeKey dk;
            if (!resolveDaughter(ch, dk)) continue;
            // Skip if Z' or A' became negative or implausible (e.g. Geant4 can
            // produce SF entries we filtered out earlier).
            if (dk.Z < 0 || dk.A <= 0) continue;
            bool isNew = (fIndex.find(dk) == fIndex.end());
            int dIdx = addNode(dk);
            // refresh n reference: addNode may have reallocated fNodes
            ChainNode& parent = fNodes[idx];
            parent.edges.push_back(ChainNode::Edge{
                /*daughterIndex*/   dIdx,
                /*branchingRatio*/  ch.branchingRatio,
                /*mode*/            ch.mode,
                /*daughterExcit*/   ch.daughterExcitation,
                /*parentChannelIdx*/ static_cast<int>(ci),
            });
            if (isNew) q.push(dIdx);
        }
    }
}

} // namespace g4gamma
