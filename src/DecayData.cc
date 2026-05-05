#include "g4gamma/DecayData.hh"
#include "g4gamma/Units.hh"
#include <fstream>
#include <sstream>
#include <stdexcept>
#include <cmath>
#include <algorithm>
#include <filesystem>

namespace g4gamma {

const char* decayModeName(DecayMode m) {
    switch (m) {
        case DecayMode::IT: return "IT";
        case DecayMode::BetaMinus: return "B-";
        case DecayMode::BetaPlus: return "B+";
        case DecayMode::KshellEC: return "EC(K)";
        case DecayMode::LshellEC: return "EC(L)";
        case DecayMode::MshellEC: return "EC(M)";
        case DecayMode::NshellEC: return "EC(N)";
        case DecayMode::Alpha: return "alpha";
        case DecayMode::Proton: return "p";
        case DecayMode::Neutron: return "n";
        case DecayMode::SpFission: return "SF";
        case DecayMode::Triton: return "t";
        case DecayMode::BDProton: return "BDp";
        case DecayMode::BDNeutron: return "BDn";
        case DecayMode::Beta2Minus: return "2B-";
        case DecayMode::Beta2Plus: return "2B+";
        case DecayMode::Proton2: return "2p";
        case DecayMode::Neutron2: return "2n";
        default: return "?";
    }
}

bool changesZA(DecayMode m) {
    return m != DecayMode::IT;
}

// Map the integer mode code from the file into our enum.
// The codes match G4RadioactiveDecayMode in Geant4.
static DecayMode modeFromInt(int v) {
    switch (v) {
        case 0:  return DecayMode::IT;
        case 1:  return DecayMode::BetaMinus;
        case 2:  return DecayMode::BetaPlus;
        case 3:  return DecayMode::KshellEC;
        case 4:  return DecayMode::LshellEC;
        case 5:  return DecayMode::MshellEC;
        case 6:  return DecayMode::NshellEC;
        case 7:  return DecayMode::Alpha;
        case 8:  return DecayMode::Proton;
        case 9:  return DecayMode::Neutron;
        case 10: return DecayMode::SpFission;
        case 11: return DecayMode::BDProton;
        case 12: return DecayMode::BDNeutron;
        case 13: return DecayMode::Beta2Minus;
        case 14: return DecayMode::Beta2Plus;
        case 15: return DecayMode::Proton2;
        case 16: return DecayMode::Neutron2;
        case 17: return DecayMode::Triton;
        default: return DecayMode::Unknown;
    }
}

// Determine the (Z',A') of the daughter based on decay mode. Returns true
// on success. We do NOT set isomer level M here; that is done in
// ChainBuilder using the daughterExcitation.
static bool daughterZA(DecayMode mode, int Z, int A, int& Zp, int& Ap) {
    Zp = Z; Ap = A;
    switch (mode) {
        case DecayMode::IT:                                     return true;
        case DecayMode::BetaMinus: Zp = Z + 1;                  return true;
        case DecayMode::BetaPlus:
        case DecayMode::KshellEC:
        case DecayMode::LshellEC:
        case DecayMode::MshellEC:
        case DecayMode::NshellEC:  Zp = Z - 1;                  return true;
        case DecayMode::Alpha:     Zp = Z - 2; Ap = A - 4;      return true;
        case DecayMode::Proton:    Zp = Z - 1; Ap = A - 1;      return true;
        case DecayMode::Neutron:                Ap = A - 1;     return true;
        case DecayMode::Triton:    Zp = Z - 1; Ap = A - 3;      return true;
        case DecayMode::BDProton:  Zp = Z;     Ap = A - 1;      return true;
        case DecayMode::BDNeutron: Zp = Z + 1; Ap = A - 1;      return true;
        case DecayMode::Beta2Minus:Zp = Z + 2;                  return true;
        case DecayMode::Beta2Plus: Zp = Z - 2;                  return true;
        case DecayMode::Proton2:   Zp = Z - 2; Ap = A - 2;      return true;
        case DecayMode::Neutron2:               Ap = A - 2;     return true;
        case DecayMode::SpFission: return false;  // products handled separately; ignore for gammas
        default: return false;
    }
}

DecayDataLoader::DecayDataLoader(std::string radDir) : fDir(std::move(radDir)) {}

const std::vector<DecayParent>* DecayDataLoader::load(int Z, int A) {
    int key = 1000 * A + Z;
    auto it = fCache.find(key);
    if (it != fCache.end()) return &it->second;

    std::ostringstream oss;
    oss << fDir << "/z" << Z << ".a" << A;
    std::ifstream in(oss.str());
    if (!in) {
        // Cache an empty entry so we don't keep retrying.
        fCache[key] = {};
        return nullptr;
    }

    using namespace units;

    std::vector<DecayParent> parents;
    DecayParent currentParent;
    bool inParent = false;

    // For each parent block we need to: collect all per-mode totals
    // (modeTotalBR), and the per-channel sums (modeSumBR), then renormalize
    // each channel's BR by modeTotalBR / modeSumBR. This matches Geant4.
    constexpr int kNumModes = 32;
    std::vector<double> modeTotalBR(kNumModes, 0.0);
    std::vector<double> modeSumBR(kNumModes, 0.0);

    auto finalizeParent = [&]() {
        if (!inParent) return;
        for (auto& ch : currentParent.channels) {
            int idx = static_cast<int>(ch.mode);
            if (idx >= 0 && idx < kNumModes) {
                double sum = modeSumBR[idx];
                double tot = modeTotalBR[idx];
                if (ch.mode == DecayMode::IT) {
                    // IT channels are already absolute (no per-channel sum).
                } else if (sum > 0.0 && tot > 0.0) {
                    ch.branchingRatio *= tot / sum;
                }
            }
        }
        // Drop channels with zero or near-zero BR.
        currentParent.channels.erase(
            std::remove_if(currentParent.channels.begin(), currentParent.channels.end(),
                [](const DecayChannel& c){ return !(c.branchingRatio > 0.0); }),
            currentParent.channels.end());
        currentParent.isStable = currentParent.channels.empty();
        parents.push_back(std::move(currentParent));
        currentParent = DecayParent{};
        std::fill(modeTotalBR.begin(), modeTotalBR.end(), 0.0);
        std::fill(modeSumBR.begin(),   modeSumBR.end(),   0.0);
        inParent = false;
    };

    std::string line;
    while (std::getline(in, line)) {
        // strip trailing whitespace
        while (!line.empty() &&
               (line.back() == '\r' || line.back() == ' ' || line.back() == '\t'))
            line.pop_back();
        if (line.empty() || line[0] == '#') continue;

        std::istringstream ss(line);

        if (line[0] == 'P') {
            // Finish previous block first.
            finalizeParent();

            std::string recordType, floatingFlag;
            double parentExcitation_keV = 0.0;
            double halfLife = 0.0;
            ss >> recordType >> parentExcitation_keV >> floatingFlag >> halfLife;
            currentParent.parentExcitation = parentExcitation_keV * keV;
            currentParent.isotope.Z = Z;
            currentParent.isotope.A = A;
            currentParent.isotope.M = 0;  // M assigned later from order
            // The half-life field in these files is in seconds (nominal), but
            // it's NOT actually used by Geant4 for decay timing -- Geant4 uses
            // ENSDFSTATE for that. We capture it for reference. Some files have
            // 0 here, indicating the half-life is in ENSDFSTATE.
            currentParent.halfLife = halfLife * second;
            currentParent.meanLife =
                (halfLife > 0.0) ? (halfLife * second) / std::log(2.0) : 0.0;
            inParent = true;
            continue;
        }

        if (!inParent) continue;

        // Two flavours of subsequent lines:
        //   - short (< 72 chars): "<modeInt> <dummy> <decayModeTotal>"
        //   - long  (>= 72):      "<modeInt> <a=daughterExcitation_keV>
        //                          <floatFlag> <b=BR%> <c=Q_keV> [betaType]"
        if (line.length() < 72) {
            int modeInt = -1;
            double dummy = 0.0, decayModeTotal = 0.0;
            ss >> modeInt >> dummy >> decayModeTotal;
            DecayMode m = modeFromInt(modeInt);

            if (m == DecayMode::IT) {
                // IT shows up only as a short line: insert a single channel.
                DecayChannel ch{};
                ch.mode = DecayMode::IT;
                ch.branchingRatio = decayModeTotal;
                int Zp, Ap;
                daughterZA(m, Z, A, Zp, Ap);
                ch.daughter = IsotopeKey{Zp, Ap, 0};
                ch.daughterExcitation = 0.0;
                currentParent.channels.push_back(ch);
            } else if (m != DecayMode::Unknown) {
                int idx = static_cast<int>(m);
                if (idx >= 0 && idx < kNumModes) {
                    modeTotalBR[idx] = decayModeTotal;
                }
            }
        } else {
            // long line: a per-channel record (specific daughter excitation).
            int modeInt = -1;
            double a = 0.0, b = 0.0, c = 0.0;
            std::string daughterFloatFlag;
            // betaType (optional; present only if line >= 84). We don't use it.
            ss >> modeInt >> a >> daughterFloatFlag >> b >> c;

            DecayMode m = modeFromInt(modeInt);
            if (m == DecayMode::Unknown) continue;
            // a in keV; b in percent; c in keV; per Geant4: a/=1000, c/=1000, b/=100.
            double daughterEx = a * keV;       // already converted via *keV
            double br         = b / 100.0;     // percent -> fraction
            // c is Q-value, we ignore.
            (void)c;

            int Zp, Ap;
            if (!daughterZA(m, Z, A, Zp, Ap)) continue;

            int idx = static_cast<int>(m);
            if (idx >= 0 && idx < kNumModes) {
                modeSumBR[idx] += br;
            }

            DecayChannel ch{};
            ch.mode = m;
            ch.branchingRatio = br;
            ch.daughter = IsotopeKey{Zp, Ap, 0}; // M filled by ChainBuilder
            ch.daughterExcitation = daughterEx;
            currentParent.channels.push_back(ch);
        }
    }

    finalizeParent();

    // Assign M=0,1,2,... to successive parent blocks. The first block is the
    // ground state, subsequent ones are isomers ordered by appearance.
    // (Geant4 uses ENSDF data + level energy + floating-level convention; we
    // approximate by enumerating order.)
    for (size_t i = 0; i < parents.size(); ++i) {
        parents[i].isotope.M = static_cast<int>(i);
    }

    auto [iter, _] = fCache.emplace(key, std::move(parents));
    return &iter->second;
}

const DecayParent* DecayDataLoader::get(const IsotopeKey& k) {
    const auto* v = load(k.Z, k.A);
    if (!v) return nullptr;
    if (k.M >= 0 && k.M < static_cast<int>(v->size())) return &(*v)[k.M];
    return nullptr;
}

} // namespace g4gamma
