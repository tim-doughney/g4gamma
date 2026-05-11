#include "g4gamma/ChainBuilder.hh"
#include <queue>

namespace g4gamma {

ChainBuilder::ChainBuilder(IDecayProvider& p) : fProvider(p) {}

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
    if (const auto* info = fProvider.get(key)) {
        n.meanLife = info->meanLife;
        n.stable   = info->stable;
    } else {
        n.meanLife = 0.0;
        n.stable   = true;
    }
    fNodes.push_back(std::move(n));
    fIndex[key] = fNodes.back().index;
    return fNodes.back().index;
}

void ChainBuilder::build(const IsotopeKey& root,
                         int maxDepth,
                         const std::vector<IsotopeKey>& cutoffs,
                         int depthLimit) {
    fNodes.clear();
    fIndex.clear();
    addNode(root);

    // Build O(1) cutoff lookup
    std::unordered_set<IsotopeKey> cutoffSet(cutoffs.begin(), cutoffs.end());

    // BFS depth of each node index (minimum hops from root)
    std::unordered_map<int, int> nodeDepth;
    nodeDepth[0] = 0;

    std::queue<int> q;
    q.push(0);

    int processed = 0;
    while (!q.empty() && processed < maxDepth * 1000) {
        int idx = q.front(); q.pop();
        ++processed;

        if (fNodes[idx].stable) continue;

        int d = nodeDepth.count(idx) ? nodeDepth.at(idx) : 0;

        // Check truncation: depth limit or explicit cutoff isotope.
        // The node itself is kept; we just don't expand its daughters.
        bool stopHere = (depthLimit >= 0 && d >= depthLimit) ||
                        cutoffSet.count(fNodes[idx].isotope);
        if (stopHere) {
            fNodes[idx].cutoff = true;
            continue;
        }

        const auto* info = fProvider.get(fNodes[idx].isotope);
        if (!info || info->stable) { fNodes[idx].stable = true; continue; }

        for (size_t bi = 0; bi < info->branches.size(); ++bi) {
            const auto& br = info->branches[bi];
            if (!br.terminals.empty()) {
                // Cascade terminated at ≥1 intermediate isomers: route daughter
                // activity proportionally, but keep the full branchingRatio for
                // emission accounting (gammas are already in br.emissions).
                for (const auto& [tk, frac] : br.terminals) {
                    if (tk.Z < 0 || tk.A <= 0 || frac < 1e-9) continue;
                    bool isNew = (fIndex.find(tk) == fIndex.end());
                    int dIdx = addNode(tk);
                    fNodes[idx].edges.push_back(ChainNode::Edge{
                        dIdx, br.branchingRatio * frac, br.mode, static_cast<int>(bi)
                    });
                    if (isNew) {
                        nodeDepth[dIdx] = d + 1;
                        q.push(dIdx);
                    }
                }
            } else {
                const auto& dk = br.daughter;
                if (dk.Z < 0 || dk.A <= 0) continue;
                bool isNew = (fIndex.find(dk) == fIndex.end());
                int dIdx = addNode(dk);
                fNodes[idx].edges.push_back(ChainNode::Edge{
                    dIdx, br.branchingRatio, br.mode, static_cast<int>(bi)
                });
                if (isNew) {
                    nodeDepth[dIdx] = d + 1;
                    q.push(dIdx);
                }
            }
        }
    }
}

} // namespace g4gamma
