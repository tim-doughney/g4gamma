// ChainBuilder.hh -- walks the decay chain from a root nuclide using any
// IDecayProvider.
#pragma once

#include "g4gamma/IDecayProvider.hh"
#include <unordered_map>
#include <unordered_set>
#include <vector>

namespace g4gamma {

struct ChainNode {
    IsotopeKey isotope;
    double     meanLife = 0.0;        // internal units (ns); 0 = unknown
    bool       stable   = true;
    bool       cutoff   = false;      // true if chain was truncated at this node
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

    // Build the chain by BFS from root.
    //
    // maxDepth      — safety cap on total nodes processed (maxDepth*1000); not
    //                 a physics depth limit.
    // cutoffs       — isotopes at which to stop expanding (the cutoff node IS
    //                 included; its daughters are NOT). Empty = no cutoffs.
    // depthLimit    — max decay generations from root (-1 = unlimited; 0 =
    //                 root only; 1 = root + daughters; etc.).
    void build(const IsotopeKey& root,
               int maxDepth = 50,
               const std::vector<IsotopeKey>& cutoffs = {},
               int depthLimit = -1);

    const std::vector<ChainNode>& nodes() const { return fNodes; }
    int indexOf(const IsotopeKey& k) const;

private:
    IDecayProvider& fProvider;
    std::vector<ChainNode> fNodes;
    std::unordered_map<IsotopeKey, int> fIndex;
    int addNode(const IsotopeKey& key);
};

} // namespace g4gamma
