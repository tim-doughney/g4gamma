// ChainBuilder.hh -- recursively walks the decay chain from a starting parent,
// assigning isomer indices (M) consistently and reducing transitions into a
// transition matrix suitable for Bateman.
//
// The decay-data files give us, for each parent (Z,A,M), a set of channels
// that produce a daughter (Z',A') in some excitation level. Whether that
// daughter level is treated as a separate isomer (M' >= 1) depends on its
// half-life and the user's config (we use the same threshold as Geant4:
// 1 ns mean lifetime, mirroring G4NuclideTable::SetMeanLifeThreshold).
//
// Algorithm:
//   - Start with the input nuclide. Push it onto a queue.
//   - For each parent, iterate channels:
//       * Compute (Z', A')
//       * Look up the daughter's level scheme in PhotonEvap to find the level
//         closest to daughterExcitation. If that level has mean lifetime
//         > threshold, treat it as an isomer (assign M' >= 1 by matching to a
//         decay-data block with the same parentExcitation). Otherwise treat the
//         daughter as ground state (M' = 0).
//       * Add the (parent -> daughter, BR) edge.
//       * If daughter not yet in the chain, enqueue it.
//   - Stop when all daughters have been processed or are stable.
#pragma once

#include "g4gamma/IsotopeKey.hh"
#include "g4gamma/DecayData.hh"
#include "g4gamma/PhotonEvap.hh"
#include "g4gamma/EnsdfState.hh"
#include <vector>
#include <unordered_map>

namespace g4gamma {

// One node in the decay chain.
struct ChainNode {
    IsotopeKey isotope;
    double     meanLife = 0.0;          // internal units; 0 = stable
    bool       stable   = true;
    // Where this node fits in the chain: index into the chain vector (assigned
    // by ChainBuilder).
    int        index    = -1;
    // Outgoing edges: (channel BR, daughter's chain index, channel-info copy)
    struct Edge {
        int     daughterIndex;
        double  branchingRatio;     // sum over edges out of this node = 1 (if unstable)
        DecayMode mode;
        double  daughterExcitation; // MeV; used by GammaSpectrum to determine cascade start
        // The channel index in the parent's DecayParent::channels list, so we
        // can look up the original record for gamma-yield computation.
        int     parentChannelIdx;
    };
    std::vector<Edge> edges;
};

class ChainBuilder {
public:
    // ensdf may be nullptr; when non-null, it is used to assign isomer M
    // indices consistently with the ENSDFSTATE convention.
    ChainBuilder(DecayDataLoader& dec, PhotonEvapData& evap,
                 EnsdfStateLoader* ensdf,
                 double isomerLifetimeThreshold /* internal units, e.g. 1*ns */);

    // Build the chain rooted at `root`. Returns the ordered chain (root first).
    // After this, `nodes()` contains the chain in topological order (parents
    // before any daughter that has activity flowing into it).
    void build(const IsotopeKey& root, int maxDepth = 50);

    const std::vector<ChainNode>& nodes() const { return fNodes; }

    // Index into fNodes given an isotope key, or -1 if not in chain.
    int indexOf(const IsotopeKey& k) const;

private:
    DecayDataLoader&  fDec;
    PhotonEvapData&   fEvap;
    EnsdfStateLoader* fEnsdf;     // may be null
    double            fIsomerThreshold;
    std::vector<ChainNode> fNodes;
    std::unordered_map<IsotopeKey, int> fIndex;

    // Resolve a decay channel's daughter into a concrete IsotopeKey including
    // the M index. Returns true on success.
    bool resolveDaughter(const DecayChannel& ch, IsotopeKey& outKey);

    // Get-or-add a node for `key`, returning its index in fNodes.
    int addNode(const IsotopeKey& key);

public:
    // Find the DecayParent block for an IsotopeKey, accounting for the
    // possibility that the file has no ground-state P block (so M->index
    // mapping is not just identity). Uses ENSDFSTATE when available.
    // Returns nullptr if not found.
    const DecayParent* parentFor(const IsotopeKey& key);
};

} // namespace g4gamma
