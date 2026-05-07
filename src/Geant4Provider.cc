#include "g4gamma/Geant4Provider.hh"
#include "g4gamma/DataPath.hh"
#include <iostream>
#include <stdexcept>
#include <cmath>
#include <limits>
#include <algorithm>

namespace g4gamma {

Geant4Provider::Geant4Provider(const Geant4ProviderOptions& opts) : fOpts(opts) {
    // Resolve directories. RadioactiveDecay and PhotonEvaporation are
    // mandatory; ENSDFSTATE is recommended; G4LEDATA is only needed if X-rays
    // are requested.
    if (fOpts.radDir.empty())   fOpts.radDir   = DataPath::radioactiveDecayDir(fOpts.geant4Sh);
    if (fOpts.evapDir.empty())  fOpts.evapDir  = DataPath::photonEvaporationDir(fOpts.geant4Sh);
    if (fOpts.ensdfDir.empty()) {
        try { fOpts.ensdfDir = DataPath::ensdfStateDir(fOpts.geant4Sh); }
        catch (...) {}
    }
    if (fOpts.includeXrays && fOpts.ledataDir.empty()) {
        try { fOpts.ledataDir = DataPath::lowEnergyDir(fOpts.geant4Sh); }
        catch (const std::exception& e) {
            std::cerr << "[g4gamma] WARNING: X-rays requested but " << e.what()
                      << " -- continuing without X-rays.\n";
            fOpts.includeXrays = false;
        }
    }
    fDecay = std::make_unique<DecayDataLoader>(fOpts.radDir);
    fEvap  = std::make_unique<PhotonEvapData>(fOpts.evapDir);
    fEnsdf = std::make_unique<EnsdfStateLoader>(fOpts.ensdfDir);
    if (fOpts.includeXrays) {
        fFluor = std::make_unique<FluorDataLoader>(fOpts.ledataDir, fOpts.fluorSubdir);
    }
}

const ParentDecayInfo* Geant4Provider::get(const IsotopeKey& key) {
    auto it = fCache.find(key);
    if (it != fCache.end()) return &it->second;
    return compute(key);
}

const DecayParent* Geant4Provider::parentFor(const IsotopeKey& key) {
    const auto* parents = fDecay->load(key.Z, key.A);
    if (!parents || parents->empty()) return nullptr;

    // Determine the target excitation for this M from ENSDFSTATE if available.
    double targetExc = -1.0;
    if (fEnsdf->ready()) {
        const auto* lvls = fEnsdf->load(key.Z, key.A);
        if (lvls) {
            int M = 0;
            bool seenGround = false;
            for (const auto& lvl : *lvls) {
                int thisM;
                if (lvl.excitation == 0.0 && !seenGround) {
                    thisM = 0;
                    seenGround = true;
                } else if (lvl.excitation > 0.0 && lvl.meanLife > fOpts.isomerThreshold) {
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

IsotopeKey Geant4Provider::resolveDaughter(const DecayChannel& ch, int /*parentZ*/, int /*parentA*/) {
    int Zp = ch.daughter.Z;
    int Ap = ch.daughter.A;
    double exc = ch.daughterExcitation;

    // IT: same nuclide, ground state.
    if (ch.mode == DecayMode::IT) {
        return IsotopeKey{Zp, Ap, 0};
    }

    // ENSDFSTATE-based mapping
    int M = 0;
    if (fEnsdf->ready()) {
        int m = fEnsdf->excitationToM(Zp, Ap, exc, fOpts.isomerThreshold,
                                       fOpts.levelTolerance);
        if (m >= 0) M = m;
    } else {
        // Fall back to file-order match in the daughter's decay file.
        if (const auto* parents = fDecay->load(Zp, Ap)) {
            if (!parents->empty()) {
                double bestDiff = std::numeric_limits<double>::infinity();
                int    bestM    = 0;
                for (size_t i = 0; i < parents->size(); ++i) {
                    double d = std::abs((*parents)[i].parentExcitation - exc);
                    if (d < bestDiff) { bestDiff = d; bestM = static_cast<int>(i); }
                }
                if (bestM > 0 && bestDiff < fOpts.levelTolerance) M = bestM;
            }
        }
    }
    return IsotopeKey{Zp, Ap, M};
}

const Geant4Provider::Cascade&
Geant4Provider::getCascade(int Z, int A, int startLevel) {
    auto key = std::make_tuple(Z, A, startLevel);
    auto it = fCascadeCache.find(key);
    if (it != fCascadeCache.end()) return it->second;

    Cascade em;
    const auto* lvls = fEvap->load(Z, A);
    if (!lvls || lvls->empty() || startLevel <= 0 ||
        startLevel >= static_cast<int>(lvls->size())) {
        fCascadeCache[key] = em;
        return fCascadeCache[key];
    }

    int N = static_cast<int>(lvls->size());
    std::vector<double> p(N, 0.0);
    p[startLevel] = 1.0;

    for (int i = N - 1; i > 0; --i) {
        if (p[i] <= 0.0) continue;
        const auto& lvl = (*lvls)[i];
        if (lvl.transitions.empty()) continue;
        for (const auto& tr : lvl.transitions) {
            if (tr.selectionProb <= 0.0) continue;
            if (tr.lowerIndex < 0 || tr.lowerIndex >= N) continue;
            double pSel = p[i] * tr.selectionProb;
            // Photon
            double gammaCount = pSel * tr.gammaEmitProb;
            if (gammaCount > 0.0) {
                em.energies.push_back(tr.gammaEnergy);
                em.intensities.push_back(gammaCount);
            }
            // IC -> shell vacancy
            double icCount = pSel * (1.0 - tr.gammaEmitProb);
            if (icCount > 0.0 && tr.hasICC) {
                double sumW = 0.0;
                for (int k = 0; k < 10; ++k) sumW += tr.iccWeights[k];
                if (sumW > 0.0) {
                    em.kVacancyExpected += icCount * tr.iccWeights[0] / sumW;
                }
            }
            p[tr.lowerIndex] += pSel;
        }
    }
    fCascadeCache[key] = std::move(em);
    return fCascadeCache[key];
}

const ParentDecayInfo* Geant4Provider::compute(const IsotopeKey& key) {
    ParentDecayInfo info;
    info.isotope = key;

    const DecayParent* dp = parentFor(key);
    if (!dp) {
        info.stable = true;
        info.meanLife = 0.0;
        auto [it, _] = fCache.emplace(key, std::move(info));
        return &it->second;
    }

    info.stable = dp->isStable;
    info.meanLife = dp->meanLife;

    // Patch in mean life from ENSDFSTATE if zero (RadioactiveDecay file
    // half-life column is a placeholder).
    if (info.meanLife <= 0.0 && !info.stable && fEnsdf->ready()) {
        double ml = fEnsdf->meanLife(key, fOpts.isomerThreshold);
        if (ml > 0.0) info.meanLife = ml;
    }

    if (info.stable) {
        auto [it, _] = fCache.emplace(key, std::move(info));
        return &it->second;
    }

    // Build branches by walking decay channels and firing cascades.
    info.branches.reserve(dp->channels.size());
    for (const auto& ch : dp->channels) {
        DecayBranch br;
        br.mode = ch.mode;
        br.branchingRatio = ch.branchingRatio;
        br.daughter = resolveDaughter(ch, key.Z, key.A);

        // Determine cascade source
        int cZ, cA;
        double startE;
        if (ch.mode == DecayMode::IT) {
            cZ = key.Z; cA = key.A;
            startE = dp->parentExcitation;
        } else {
            cZ = ch.daughter.Z; cA = ch.daughter.A;
            startE = ch.daughterExcitation;
        }

        // If the daughter resolved to a tracked isomer (M >= 1), the cascade
        // is deferred to the isomer's own IT decay.
        bool deferToIsomer = (ch.mode != DecayMode::IT) &&
                              (br.daughter.M >= 1) &&
                              (ch.daughterExcitation > 0.0);
        double kVacancy = 0.0;
        if (!deferToIsomer) {
            int startLevel = fEvap->findLevel(cZ, cA, startE, fOpts.levelTolerance);
            if (startLevel < 0) startLevel = 0;
            const auto& cas = getCascade(cZ, cA, startLevel);
            for (size_t i = 0; i < cas.energies.size(); ++i) {
                Emission em;
                em.type      = EmissionType::Gamma;
                em.energy    = cas.energies[i];
                em.intensity = cas.intensities[i];
                br.emissions.push_back(em);
            }
            kVacancy += cas.kVacancyExpected;
        }

        // K-shell EC produces a K-vacancy in the daughter atom directly.
        if (fOpts.includeXrays && ch.mode == DecayMode::KshellEC) {
            kVacancy += 1.0;
        }

        // Convert K-vacancy yield to K X-rays
        if (fOpts.includeXrays && kVacancy > 0.0 && fFluor) {
            // Z for fluorescence: daughter Z if EC channel, else parent
            int Zfluor = key.Z;
            if (ch.mode == DecayMode::KshellEC ||
                ch.mode == DecayMode::LshellEC ||
                ch.mode == DecayMode::MshellEC ||
                ch.mode == DecayMode::NshellEC ||
                ch.mode == DecayMode::EC ||
                ch.mode == DecayMode::BetaPlus) {
                Zfluor = ch.daughter.Z;
            }
            const auto* vacs = fFluor->load(Zfluor);
            if (vacs && !vacs->empty()) {
                const auto& kvac = (*vacs)[0]; // K is index 0 in EADL
                for (const auto& tr : kvac.transitions) {
                    Emission em;
                    em.type      = EmissionType::XRay;
                    em.energy    = tr.energy;
                    em.intensity = kVacancy * tr.prob;
                    br.emissions.push_back(em);
                }
            }
        }

        info.branches.push_back(std::move(br));
    }

    auto [it, _] = fCache.emplace(key, std::move(info));
    return &it->second;
}

} // namespace g4gamma
