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
    if ((fOpts.includeXrays || fOpts.fullXrayCascade) && fOpts.ledataDir.empty()) {
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
    if (fOpts.includeXrays || fOpts.fullXrayCascade) {
        fFluor = std::make_unique<FluorDataLoader>(fOpts.ledataDir, fOpts.fluorSubdir);
        if (fOpts.fullXrayCascade)
            fAuger = std::make_unique<AugerDataLoader>(fOpts.ledataDir + "/auger");
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
            // Photon: Geant4 uses fAlpha from file directly: gammaEmitProb = 1/(1+fAlpha)
            double gammaCount = pSel * tr.gammaEmitProb;
            if (gammaCount > 0.0) {
                em.energies.push_back(tr.gammaEnergy);
                em.intensities.push_back(gammaCount);
            }
            // IC -> shell vacancies (track all 10 ICC shells)
            double icCount = pSel * (1.0 - tr.gammaEmitProb);
            if (icCount > 0.0 && tr.hasICC) {
                double sumW = 0.0;
                for (int k = 0; k < Cascade::N_ICC; ++k) sumW += tr.iccWeights[k];
                if (sumW > 0.0) {
                    for (int k = 0; k < Cascade::N_ICC; ++k)
                        em.shellVacancy[k] += icCount * tr.iccWeights[k] / sumW;
                }
            }
            // Geant4 stops the cascade when the destination level is long-lived
            // (LifeTime > fMaxLifeTime = 1 ns). It marks the nucleus as LongLived
            // and G4RadioactiveDecay creates the daughter as a separate radioactive ion.
            // We replicate this: do not propagate probability into a long-lived intermediate
            // level (index > 0). The probability that "stops" here belongs to the isomer's
            // own decay chain, which ChainBuilder handles separately.
            //
            // Additional guard: only stop if the level actually has a G4RADIOACTIVEDATA
            // P-block. Some levels are marginally above the 1 ns mean-lifetime threshold
            // in ENSDFSTATE but are not in the RadioactiveDecay database (e.g. Ra-224
            // at 84.4 keV, τ=1.079 ns). Geant4 de-excites those immediately via
            // G4PhotonEvaporation rather than creating a separate ion; we must continue
            // cascading through them here.
            const auto& dstLvl = (*lvls)[tr.lowerIndex];
            bool dstIsIsomer = (tr.lowerIndex > 0) &&
                               (dstLvl.meanLifeTime > fOpts.isomerThreshold);
            if (dstIsIsomer) {
                // Verify the level is actually tracked by G4RADIOACTIVEDATA.
                bool hasDecayEntry = false;
                if (const auto* parents = fDecay->load(Z, A)) {
                    for (const auto& par : *parents) {
                        if (std::abs(par.parentExcitation - dstLvl.energy) <
                                fOpts.levelTolerance) {
                            hasDecayEntry = true;
                            break;
                        }
                    }
                }
                if (!hasDecayEntry) dstIsIsomer = false;
            }
            if (!dstIsIsomer) {
                p[tr.lowerIndex] += pSel;
            } else {
                // Track probability that stopped at this isomeric intermediate level.
                bool found = false;
                for (auto& [li, tp] : em.terminalProb) {
                    if (li == tr.lowerIndex) { tp += pSel; found = true; break; }
                }
                if (!found) em.terminalProb.emplace_back(tr.lowerIndex, pSel);
            }
        }
    }
    // Probability that cascaded all the way to the ground state.
    if (p[0] > 1e-9)
        em.terminalProb.emplace_back(0, p[0]);

    fCascadeCache[key] = std::move(em);
    return fCascadeCache[key];
}

IsotopeKey Geant4Provider::levelIndexToKey(int Z, int A, int levelIdx) {
    if (levelIdx == 0) return IsotopeKey{Z, A, 0};
    const auto* lvls = fEvap->load(Z, A);
    if (!lvls || levelIdx >= static_cast<int>(lvls->size()))
        return IsotopeKey{Z, A, 0};
    double exc = (*lvls)[levelIdx].energy;
    int M = 0;
    if (fEnsdf->ready()) {
        int m = fEnsdf->excitationToM(Z, A, exc, fOpts.isomerThreshold, fOpts.levelTolerance);
        if (m >= 0) M = m;
    }
    return IsotopeKey{Z, A, M};
}

// ICC index (0=K,1=L1,2=L2,3=L3,4=M1..8=M5,9=N+) → EADL shell ID in fl-tr-pr files.
// G4EMLOW uses the old EADL96 convention where group shells occupy the "missing" IDs:
// 9=M2+M3 group, 12=M4+M5 group, 15=N group, 17=N2+N3 group — none of these appear
// as individual vacancy shells in fl-tr-pr files.
static constexpr int kICCtoEADL[10] = {
    1, 3, 5, 6, 8, 10, 11, 13, 14, 16
};

int Geant4Provider::eadlToICC(int eadlId) {
    switch (eadlId) {
        case  1: return 0;   // K
        case  3: return 1;   // L1
        case  5: return 2;   // L2
        case  6: return 3;   // L3
        case  8: return 4;   // M1
        case 10: return 5;   // M2
        case 11: return 6;   // M3
        case 13: return 7;   // M4
        case 14: return 8;   // M5
        case 16: return 9;   // N1 → N+ bucket
        default: return -1;  // outer shells or group IDs — terminate cascade
    }
}

void Geant4Provider::appendXRays(
    int Zfluor,
    std::array<double, Cascade::N_ICC> shellVac,   // copy, modified in-place
    std::vector<Emission>& out)
{
    if (!fFluor) return;

    if (fOpts.fullXrayCascade) {
        // Forward sweep K→L→M→N: each fluorescence transition creates a
        // secondary vacancy in the donor shell, which is then processed in
        // turn.  Auger transitions are photon-dark and are dropped (secondary
        // vacancies from Auger are not tracked — this is the approximation).
        for (int icc = 0; icc < Cascade::N_ICC; ++icc) {
            if (shellVac[icc] <= 0.0) continue;
            const FluorVacancy* vac = fFluor->findVacancy(Zfluor, kICCtoEADL[icc]);
            if (!vac || vac->transitions.empty()) continue;
            double N = shellVac[icc];
            for (const auto& tr : vac->transitions) {
                double count = N * tr.prob;
                if (count <= 0.0) continue;
                Emission em;
                em.type     = EmissionType::XRay;
                em.energy   = tr.energy;
                em.intensity = count;
                out.push_back(em);
                // Propagate secondary vacancy into the donor shell
                int secICC = eadlToICC(tr.originShell);
                if (secICC > icc && secICC < Cascade::N_ICC)
                    shellVac[secICC] += count;
            }
            // Auger/Coster-Kronig: non-radiative transitions create two secondary
            // vacancies in outer shells. Propagate them so downstream shells see
            // the correct vacancy population.
            if (fAuger) {
                const AugerVacancy* av = fAuger->findVacancy(Zfluor, kICCtoEADL[icc]);
                if (av) {
                    for (const auto& atr : av->transitions) {
                        double count = N * atr.prob;
                        if (count <= 0.0) continue;
                        int sec1 = eadlToICC(atr.secShell1);
                        if (sec1 > icc && sec1 < Cascade::N_ICC)
                            shellVac[sec1] += count;
                        int sec2 = eadlToICC(atr.secShell2);
                        if (sec2 > icc && sec2 < Cascade::N_ICC)
                            shellVac[sec2] += count;
                    }
                }
            }
        }
    } else {
        // K-shell only (original behaviour).
        if (shellVac[0] <= 0.0) return;
        const FluorVacancy* kvac = fFluor->findVacancy(Zfluor, 1); // EADL 1 = K
        if (!kvac) return;
        for (const auto& tr : kvac->transitions) {
            Emission em;
            em.type      = EmissionType::XRay;
            em.energy    = tr.energy;
            em.intensity = shellVac[0] * tr.prob;
            out.push_back(em);
        }
    }
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

        // If the daughter resolved to a tracked isomer (M >= 1) AND that isomer
        // has a G4RADIOACTIVEDATA entry, defer the cascade to the isomer's own
        // IT decay. If no radioactive-decay data exists for the isomeric level,
        // Geant4 de-excites it immediately via G4PhotonEvaporation; we fire the
        // cascade here instead (same gammas, activity routed to ground state).
        bool deferToIsomer = (ch.mode != DecayMode::IT) &&
                              (br.daughter.M >= 1) &&
                              (ch.daughterExcitation > 0.0) &&
                              (parentFor(br.daughter) != nullptr);
        bool needXrays = fOpts.includeXrays || fOpts.fullXrayCascade;

        // IC vacancies from the photon-evaporation cascade.
        std::array<double, Cascade::N_ICC> shellVac = {};
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
            if (needXrays) shellVac = cas.shellVacancy;

            // Determine where the cascade actually terminated.  If probability
            // stopped at one or more long-lived intermediate levels (isomers) we
            // must route the daughter activity to those levels, not blindly to
            // M=0 as resolveDaughter() returned.  Record the terminal distribution
            // in br.terminals so ChainBuilder can split activity correctly.
            if (!cas.terminalProb.empty()) {
                for (const auto& [li, prob] : cas.terminalProb) {
                    if (prob < 1e-9) continue;
                    IsotopeKey tk = levelIndexToKey(cZ, cA, li);
                    bool found = false;
                    for (auto& [k, f] : br.terminals) {
                        if (k == tk) { f += prob; found = true; break; }
                    }
                    if (!found) br.terminals.emplace_back(tk, prob);
                }
            }
        }

        // IC electrons and EC/BetaPlus both produce vacancies in the daughter atom's
        // electron cloud (the nucleus has already transformed).
        int Zfluor = ch.daughter.Z;
        if (needXrays) {
            switch (ch.mode) {
                case DecayMode::KshellEC:                shellVac[0] += 1.0; break;
                case DecayMode::LshellEC:                shellVac[1] += 1.0; break; // L1
                case DecayMode::MshellEC:                shellVac[4] += 1.0; break; // M1
                case DecayMode::NshellEC:                shellVac[9] += 1.0; break; // N+
                case DecayMode::EC: case DecayMode::BetaPlus: shellVac[0] += 1.0; break;
                default: break;
            }
        }

        if (needXrays && fFluor) {
            appendXRays(Zfluor, shellVac, br.emissions);
        }

        info.branches.push_back(std::move(br));
    }

    auto [it, _] = fCache.emplace(key, std::move(info));
    return &it->second;
}

} // namespace g4gamma
