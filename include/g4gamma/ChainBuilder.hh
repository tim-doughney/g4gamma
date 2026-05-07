// ChainBuilder.hh -- walks the decay chain from a root nuclide using any
// IDecayProvider.
#pragma once

#include "g4gamma/IDecayProvider.hh"
#include <vector>
#include <unordered_map>

namespace g4gamma {

struct ChainNode {
    IsotopeKey isotope;
    double     meanLife = 0.0;        // internal units (ns); 0 = unknown
    bool       stable   = true;
    int        index    = -1;
    struct Edge {
        int        daughterIndex;     // index into chain
        double     branchingRatio;
        DecayMode  mode;
        int        parentBranchIdx;   // index into ParentDecayInfo::branches
    };
    std::vector<Edge> edges;
};

class ChainBuilder {
public:
    explicit ChainBuilder(IDecayProvider& provider);
    void build(const IsotopeKey& root, int maxDepth = 50);
    const std::vector<ChainNode>& nodes() const { return fNodes; }
    int indexOf(const IsotopeKey& k) const;

private:
    IDecayProvider& fProvider;
    std::vector<ChainNode> fNodes;
    std::unordered_map<IsotopeKey, int> fIndex;
    int addNode(const IsotopeKey& key);
};

} // namespace g4gamma
