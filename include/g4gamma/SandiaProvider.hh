// SandiaProvider.hh -- IDecayProvider backed by SandiaDecay's sandia.decay.xml.
//
// Format reference: https://github.com/sandialabs/SandiaDecay
// Copyright 2018 NTESS, LGPL-2.1.
//
// We parse the XML lazily on construction (the file is ~6 MB minified and
// parses in well under a second). The parser is hand-rolled and minimal --
// rapidxml would work but introducing it as a dependency is not worth the
// gain on this single file.
//
// Resolution of the XML path:
//   1. constructor argument, if non-empty
//   2. environment variable SANDIA_DECAY_XML
//   3. <project_root>/data/sandia/sandia.decay.nocoinc.min.xml (bundled)
//   4. /usr/local/share/g4gamma/sandia.decay.nocoinc.min.xml (install path)
#pragma once

#include "g4gamma/IDecayProvider.hh"
#include <string>
#include <unordered_map>
#include <vector>

namespace g4gamma {

class SandiaProvider : public IDecayProvider {
public:
    // xmlPath: explicit path to sandia.decay.xml or any of its variants.
    // If empty, falls back to env / bundled-data search (see header comment).
    explicit SandiaProvider(const std::string& xmlPath = "");

    const ParentDecayInfo* get(const IsotopeKey& key) override;

    const char* name() const override { return "sandia"; }
    bool emissionsIncludeXrays() const override { return true; }
    bool emissionsIncludeAnnihilation() const override { return false; }

    // Diagnostic accessors
    const std::string& xmlPath() const { return fXmlPath; }
    size_t numNuclides() const { return fByKey.size(); }

    // Resolve the XML path used at construction (also exposed for testing).
    static std::string locateXml(const std::string& hint);

private:
    std::string fXmlPath;
    bool        fLoaded = false;
    std::unordered_map<IsotopeKey, ParentDecayInfo> fByKey;
    // Lookup symbol -> IsotopeKey, populated during parsing of the <nuclide>
    // table. We need this to resolve transition `child` and `parent` strings
    // back to (Z, A, M).
    std::unordered_map<std::string, IsotopeKey> fSymbolToKey;
    std::unordered_map<IsotopeKey, std::string> fKeyToSymbol;

    void load();
    void parseFile(const std::string& path);
};

} // namespace g4gamma
