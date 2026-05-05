// DecayData.hh -- parses the Geant4 RadioactiveDecay/zZ.aA file format.
//
// File layout (informally):
//   P  <parentExcitation_keV>  <floatFlag>  <halfLife>           # parent record
//   <mode> <dummy> <decayModeTotal>                               # short line: total BR for this mode
//   <mode> <a=daughterExcitation_keV> <floatFlag> <b=BR%> <c=Q_keV>  # long line: per-channel
//   ... (more modes / channels)
//   P  <parentExcitation>  ...                                    # next parent level entry
//
// Multiple "P ..." blocks in one file correspond to ground state and isomeric
// states (M=0,1,2,...) of the same (Z,A). The block's M index is determined by
// the order of appearance with non-trivial half-life: the ground state block
// is M=0, the first isomer with mean-life > threshold is M=1, etc.
//
// We decode:
//   - the parent's mean life (= halfLife / ln 2)
//   - the list of decay channels, each producing a daughter (Z',A',M') in a
//     specific excitation level, with a branching ratio that already includes
//     the mode total renormalization that Geant4 applies.
#pragma once

#include "g4gamma/IsotopeKey.hh"
#include <string>
#include <vector>
#include <unordered_map>
#include <memory>
#include <optional>

namespace g4gamma {

enum class DecayMode {
    Unknown = -1,
    IT      = 0,   // isomeric transition
    BetaMinus,
    BetaPlus,
    KshellEC,
    LshellEC,
    MshellEC,
    NshellEC,
    Alpha,
    Proton,
    Neutron,
    BDProton,      // beta-delayed proton
    BDNeutron,     // beta-delayed neutron
    Beta2Minus,
    Beta2Plus,
    Proton2,
    Neutron2,
    Triton,
    SpFission,
};

bool changesZA(DecayMode m);
const char* decayModeName(DecayMode m);

// One channel of one decay mode.
struct DecayChannel {
    DecayMode mode;
    double    branchingRatio;       // already mode-renormalized: sum over all
                                    // channels of all modes = 1 (for unstable parent)
    IsotopeKey daughter;            // (Z', A', daughterIsomerHint M')
                                    // M' is filled later by ChainBuilder using the
                                    // floating-level convention; for now we store
                                    // M=0 plus the explicit excitation energy below.
    double    daughterExcitation;   // in internal units (MeV)
};

struct DecayParent {
    IsotopeKey isotope;             // (Z, A, M)
    double     parentExcitation;    // MeV, internal units
    double     halfLife;            // internal units (ns); 0 means stable; +inf means user is unsure
    double     meanLife;            // = halfLife / ln(2). 0 for stable.
    bool       isStable;            // halfLife == 0 (no entries) or no decay channels
    std::vector<DecayChannel> channels;
};

class DecayDataLoader {
public:
    explicit DecayDataLoader(std::string radDir);

    // Loads (and caches) all parent blocks from zZ.aA file. Each block
    // corresponds to a different isomer level M. Returns nullptr if file
    // doesn't exist (which means: stable, or beyond data tables).
    const std::vector<DecayParent>* load(int Z, int A);

    // Convenience: get the DecayParent for a given (Z,A,M), or nullptr.
    const DecayParent* get(const IsotopeKey& key);

    const std::string& dir() const { return fDir; }

private:
    std::string fDir;
    std::unordered_map<int, std::vector<DecayParent>> fCache; // key = 1000*A + Z
};

} // namespace g4gamma
