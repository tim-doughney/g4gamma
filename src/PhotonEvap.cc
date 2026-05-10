#include "g4gamma/PhotonEvap.hh"
#include "g4gamma/Units.hh"
#include <fstream>
#include <sstream>
#include <stdexcept>
#include <cmath>
#include <algorithm>
#include <limits>

namespace g4gamma {

PhotonEvapData::PhotonEvapData(std::string evapDir) : fDir(std::move(evapDir)) {}

const std::vector<NuclearLevel>* PhotonEvapData::load(int Z, int A) {
    int key = 1000 * A + Z;
    auto it = fCache.find(key);
    if (it != fCache.end()) return &it->second;

    std::ostringstream oss;
    oss << fDir << "/z" << Z << ".a" << A;
    std::ifstream in(oss.str());
    if (!in) {
        fCache[key] = {};
        return nullptr;
    }

    using namespace units;

    std::vector<NuclearLevel> levels;

    int    levelIdx;
    std::string floatTag;

    // Format-following loop. Read until EOF or read failure.
    while (in >> levelIdx >> floatTag) {
        double E_keV    = 0.0;
        double T_half_s = 0.0;  // half-life in seconds (Geant4 G4LevelReader convention)
        double spin     = 0.0;
        int    ntrans   = 0;

        if (!(in >> E_keV >> T_half_s >> spin >> ntrans)) break;

        // Convert T_half_s → mean life in ns, matching G4LevelReader:
        //   fTime *= CLHEP::second / ln(2)   where CLHEP::second = 1e9 ns
        static const double ln2 = std::log(2.0);
        NuclearLevel lvl;
        lvl.energy = E_keV * keV;
        if (T_half_s > 0.0) {
            lvl.meanLifeTime = T_half_s / ln2 * 1e9;  // ns
        } else if (T_half_s < 0.0) {
            lvl.meanLifeTime = std::numeric_limits<double>::infinity();  // stable
        } else {
            lvl.meanLifeTime = 0.0;  // prompt / no defined lifetime
        }

        // Ensure monotonically non-decreasing energy (Geant4 enforces this).
        if (!levels.empty() && lvl.energy < levels.back().energy) {
            lvl.energy = levels.back().energy + 1e-9 * keV;
        }

        // Read transitions
        std::vector<LevelTransition> raw;
        raw.reserve(ntrans > 0 ? ntrans : 0);
        bool readOk = true;

        // Will need normalisation: norm = sum_j (1+alpha_j) * fProb_j
        double norm = 0.0;

        for (int j = 0; j < ntrans; ++j) {
            int    i2     = 0;
            double Etr_keV= 0.0;
            double fProb  = 0.0;
            int    tnum   = 0;
            double ratio  = 0.0;
            double fAlpha = 0.0;
            if (!(in >> i2 >> Etr_keV >> fProb >> tnum >> ratio >> fAlpha)) {
                readOk = false;
                break;
            }
            (void)tnum; (void)ratio;
            LevelTransition t;
            t.lowerIndex   = i2;
            t.gammaEnergy  = Etr_keV * keV;
            // Sanitise alpha
            if (!(fAlpha >= 0.0)) fAlpha = 0.0;
            t.gammaEmitProb = 1.0 / (1.0 + fAlpha);
            t.selectionProb = (1.0 + fAlpha) * fProb;  // un-normalized for now
            norm += t.selectionProb;
            t.hasICC = (fAlpha > 0.0);
            for (int k = 0; k < 10; ++k) t.iccWeights[k] = 0.0f;
            if (fAlpha > 0.0) {
                for (int k = 0; k < 10; ++k) {
                    double v = 0.0;
                    if (!(in >> v)) { readOk = false; break; }
                    t.iccWeights[k] = static_cast<float>(v);
                }
                if (!readOk) break;
            }
            // Drop transitions that point upward (broken data); Geant4 sets
            // their probability to 0.
            if (i2 >= levelIdx) t.selectionProb = 0.0;
            raw.push_back(t);
        }

        if (!readOk) break;

        // Normalize selectionProb to sum=1 (skip if all zero -> orphan level).
        if (norm > 0.0) {
            for (auto& t : raw) t.selectionProb /= norm;
        }
        lvl.transitions = std::move(raw);
        levels.push_back(std::move(lvl));
    }

    auto [iter, _] = fCache.emplace(key, std::move(levels));
    return &iter->second;
}

int PhotonEvapData::findLevel(int Z, int A, double exc, double tol) {
    const auto* v = load(Z, A);
    if (!v || v->empty()) return -1;
    if (exc < tol) return 0;
    int best = -1;
    double bestDiff = std::numeric_limits<double>::infinity();
    for (size_t i = 0; i < v->size(); ++i) {
        double d = std::abs((*v)[i].energy - exc);
        if (d < bestDiff) { bestDiff = d; best = static_cast<int>(i); }
    }
    if (bestDiff <= tol) return best;
    // If no level matches within tolerance, return the closest one (best-effort).
    // Geant4 also has fallback behavior here.
    return best;
}

} // namespace g4gamma
