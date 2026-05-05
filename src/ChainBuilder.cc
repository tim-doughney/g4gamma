#include "g4gamma/ChainBuilder.hh"
#include "g4gamma/Units.hh"
#include <queue>
#include <cmath>
#include <algorithm>
#include <limits>

namespace g4gamma {

ChainBuilder::ChainBuilder(DecayDataLoader& d, PhotonEvapData& e, double thr)
    : fDec(d), fEvap(e), fIsomerThreshold(thr) {}

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
    if (const auto* dp = fDec.get(key)) {
        n.meanLife = dp->meanLife;
        n.stable   = dp->isStable;
    } else {
        // No decay data file -> treat as stable.
        n.meanLife = 0.0;
        n.stable   = true;
    }
    fNodes.push_back(std::move(n));
    fIndex[key] = fNodes.back().index;
    return fNodes.back().index;
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

    // Decide which M-block of the daughter (Z',A') corresponds to this
    // excitation. Strategy: if the daughter's decay-data has multiple parent
    // blocks (M=0, M=1, ...), match by parentExcitation closest to exc.
    // If no match within tolerance, default to M=0.
    int M = 0;
    if (const auto* parents = fDec.load(Zp, Ap)) {
        if (!parents->empty()) {
            // First, find the candidate block with parentExcitation closest to exc.
            double bestDiff = std::numeric_limits<double>::infinity();
            int    bestM    = 0;
            for (size_t i = 0; i < parents->size(); ++i) {
                double d = std::abs((*parents)[i].parentExcitation - exc);
                if (d < bestDiff) { bestDiff = d; bestM = static_cast<int>(i); }
            }
            // Accept M >= 1 only if (a) it's actually an excited block, and
            // (b) the match is reasonably close. The Geant4 levelTolerance is
            // O(1 eV), but daughterExcitation in the decay file is rounded to
            // keV, so use a tolerance of 1 keV.
            if (bestM > 0 && bestDiff < 1.0 * units::keV) {
                M = bestM;
            } else {
                // The daughter excitation might land on an isomeric level even
                // if we'd otherwise go to M=0. Cross-check via PhotonEvap +
                // mean-life threshold.
                int li = fEvap.findLevel(Zp, Ap, exc, 1.0 * units::keV);
                if (li > 0) {
                    if (const auto* lvls = fEvap.load(Zp, Ap)) {
                        double t = (*lvls)[li].meanLifeTime;
                        if (t > fIsomerThreshold) {
                            // Long-lived excited state => treat as isomer if
                            // there is a matching decay-data block.
                            // Without one, leave as M=0 (Geant4 falls back to
                            // pure IT cascade in that case).
                            M = bestM; // best effort
                        }
                    }
                }
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

        const DecayParent* dp = fDec.get(n.isotope);
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
