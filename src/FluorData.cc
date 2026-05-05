#include "g4gamma/FluorData.hh"
#include "g4gamma/Units.hh"
#include <fstream>
#include <sstream>

namespace g4gamma {

FluorDataLoader::FluorDataLoader(std::string ledataDir, std::string fluorSubdir)
    : fDir(std::move(ledataDir) + std::move(fluorSubdir)) {}

const std::vector<FluorVacancy>* FluorDataLoader::load(int Z) {
    auto it = fCache.find(Z);
    if (it != fCache.end()) return &it->second;

    std::ostringstream oss;
    oss << fDir << "/fl-tr-pr-" << Z << ".dat";
    std::ifstream in(oss.str());
    if (!in) {
        fCache[Z] = {};
        return nullptr;
    }

    using namespace units;
    std::vector<FluorVacancy> vacancies;

    // The format per G4FluoData::LoadData is annoying:
    //   - Each "block" begins with a header row whose THREE columns are all the
    //     vacancy shell id (e.g. "1 1 1"). Then transition rows of three columns
    //     follow: <originShell> <prob> <energy_MeV>.
    //   - Block separator: row of "-1 -1 -1".
    //   - End of file: row of "-2 -2 -2".
    //
    // Our state machine:
    //   - Read first three numbers of a block as the vacancy id (only first
    //     value matters; all three are the same).
    //   - Then read transitions until we hit "-1 -1 -1" (start of next block) or
    //     "-2 -2 -2" (EOF).

    bool startBlock = true;
    FluorVacancy cur;
    while (in.good()) {
        double a, b, c;
        if (!(in >> a >> b >> c)) break;
        if (a == -2) break;  // EOF marker
        if (a == -1) {
            // End of block
            vacancies.push_back(std::move(cur));
            cur = FluorVacancy{};
            startBlock = true;
            continue;
        }
        if (startBlock) {
            cur.vacancyShell = static_cast<int>(a);
            startBlock = false;
            // The header row in the file is "<id> <id> <id>" -- so we already
            // consumed three numbers, all equal. Move on.
            continue;
        }
        FluorTransition t;
        t.originShell = static_cast<int>(a);
        t.prob        = b;
        t.energy      = c * MeV;  // file is in MeV per Geant4
        cur.transitions.push_back(t);
    }
    if (!cur.transitions.empty()) vacancies.push_back(std::move(cur));

    auto [iter, _] = fCache.emplace(Z, std::move(vacancies));
    return &iter->second;
}

double FluorDataLoader::yield(int Z, int vacancyShellIndex) {
    const auto* v = load(Z);
    if (!v) return 0.0;
    if (vacancyShellIndex < 0 || vacancyShellIndex >= static_cast<int>(v->size())) return 0.0;
    double s = 0.0;
    for (const auto& t : (*v)[vacancyShellIndex].transitions) s += t.prob;
    return s;
}

} // namespace g4gamma
