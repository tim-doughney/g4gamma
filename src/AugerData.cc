#include "g4gamma/AugerData.hh"
#include <fstream>
#include <sstream>

namespace g4gamma {

AugerDataLoader::AugerDataLoader(std::string augerDir)
    : fDir(std::move(augerDir)) {}

const std::vector<AugerVacancy>* AugerDataLoader::load(int Z) {
    auto it = fCache.find(Z);
    if (it != fCache.end()) return &it->second;

    std::ostringstream oss;
    oss << fDir << "/au-tr-pr-" << Z << ".dat";
    std::ifstream in(oss.str());
    if (!in) {
        fCache[Z] = {};
        return nullptr;
    }

    // Format per G4AugerData::LoadData:
    //   Block header: "<id> <id> <id> <id>"  (4 identical integers = vacancy shell EADL id)
    //   Transition rows: "<secShell1> <secShell2> <prob> <energy_MeV>"
    //   Block separator: "-1 -1 -1 -1"
    //   EOF marker:      "-2 -2 -2 -2"
    //
    // Distinguish header from transition: in a header all four values are equal
    // integers >= 1. In a transition, the third value (prob) is < 1.0.

    std::vector<AugerVacancy> vacancies;
    bool startBlock = true;
    AugerVacancy cur;

    while (in.good()) {
        double a, b, c, d;
        if (!(in >> a >> b >> c >> d)) break;
        if (a == -2) break;   // EOF marker
        if (a == -1) {
            // End of current block
            vacancies.push_back(std::move(cur));
            cur = AugerVacancy{};
            startBlock = true;
            continue;
        }
        if (startBlock) {
            cur.vacancyShell = static_cast<int>(a);
            startBlock = false;
            continue;
        }
        AugerTransition tr;
        tr.secShell1 = static_cast<int>(a);
        tr.secShell2 = static_cast<int>(b);
        tr.prob      = c;
        // d is Auger electron energy in MeV — not needed for photon yield
        cur.transitions.push_back(tr);
    }
    if (!cur.transitions.empty()) vacancies.push_back(std::move(cur));

    auto [iter, _] = fCache.emplace(Z, std::move(vacancies));
    return &iter->second;
}

const AugerVacancy* AugerDataLoader::findVacancy(int Z, int eadlShellId) {
    const auto* vacs = load(Z);
    if (!vacs) return nullptr;
    for (const auto& v : *vacs)
        if (v.vacancyShell == eadlShellId) return &v;
    return nullptr;
}

} // namespace g4gamma
