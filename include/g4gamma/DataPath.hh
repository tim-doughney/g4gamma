// DataPath.hh -- locate Geant4 data directories.
//
// Resolution strategy, in order:
//   1. Environment variable (e.g. G4RADIOACTIVEDATA)
//   2. Parse a geant4.sh file (default /opt/geant4/bin/geant4.sh) if provided
//   3. Search GEANT4_DATA_DIR / share dir for a versioned sub-directory
//
// Note: in recent Geant4 installs (10.7+) the per-dataset env vars are no
// longer exported by default. Instead, GEANT4_DATA_DIR points at the parent
// 'data' folder and contains versioned subdirectories such as
// 'RadioactiveDecay6.1.2', 'PhotonEvaporation6.1', 'G4EMLOW8.7', etc.
#pragma once
#include <string>

namespace g4gamma {

class DataPath {
public:
    // Path to RadioactiveDecay dataset (the directory containing zZ.aA files).
    static std::string radioactiveDecayDir(const std::string& geant4_sh_hint = "");

    // Path to PhotonEvaporation dataset (containing per-isotope level files).
    static std::string photonEvaporationDir(const std::string& geant4_sh_hint = "");

    // Path to G4LEDATA (parent of /fluor and /fluor_Bearden subdirectories).
    static std::string lowEnergyDir(const std::string& geant4_sh_hint = "");

    // Path to G4ENSDFSTATEDATA (contains ENSDFSTATE.dat). Returns empty
    // string rather than throwing if not found, since this is optional.
    static std::string ensdfStateDir(const std::string& geant4_sh_hint = "");

    // Lower-level helpers, exposed for testing
    static std::string getEnv(const char* name);
    static std::string parseGeant4Sh(const std::string& path, const std::string& var);
    static std::string findVersionedSubdir(const std::string& parent,
                                           const std::string& prefix);
};

} // namespace g4gamma
