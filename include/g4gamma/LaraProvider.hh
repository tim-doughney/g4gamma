// LaraProvider.hh -- IDecayProvider backed by LARA / DDEP per-nuclide ASCII
// files from LNE-LNHB.
//
// Reference: http://www.lnhb.fr/Laraweb/
//
// LARA is the gold-standard for radionuclide metrology. About 220 isotopes
// are DDEP-evaluated; the rest of LARAWEB's ~400 nuclides are sourced from
// ENSDF. For NORM thesis work this provider is excellent for cross-validation
// on headline isotopes (K-40, Cs-137, U/Th chain members) where DDEP
// evaluations are available.
//
// The provider supports two storage layouts:
//   1. A loose directory of `<symbol>.lara.txt` files
//   2. A tarball (`*.tar` or `*.tar.gz`) containing those files
// Tarball mode is preferred when shipping the full ~400 nuclide set since it
// keeps disk footprint small (a few hundred KB) and indexing is < 50 ms.
//
// Path resolution (in order):
//   1. constructor argument (file or directory)
//   2. environment variable LARA_DATA_DIR (file or directory)
//   3. <repo>/data/lara/lara.tar.gz, <repo>/data/lara/lara.tar
//   4. <repo>/data/lara/  (loose files)
//   5. /usr/local/share/g4gamma/lara{,.tar.gz,.tar}
#pragma once

#include "g4gamma/IDecayProvider.hh"
#include <string>
#include <vector>
#include <unordered_map>

namespace g4gamma {

class LaraProvider : public IDecayProvider {
public:
    // hint: explicit path to either a directory of *.lara.txt files OR a
    // tarball (.tar or .tar.gz) of them. If empty, auto-resolves.
    explicit LaraProvider(const std::string& hint = "");

    const ParentDecayInfo* get(const IsotopeKey& key) override;

    const char* name() const override { return "lara"; }
    bool emissionsIncludeXrays() const override { return true; }
    bool emissionsIncludeAnnihilation() const override { return true; }
    bool emissionsArePerDecay() const override { return true; }

    // Diagnostics
    const std::string& source() const { return fSource; }
    bool isTarball() const { return fIsTar; }
    size_t numAvailable() const { return fAvailable.size(); }

    // Path resolution helper (also used in tests).
    static std::string locateData(const std::string& hint);

    // Convert a LARA-style symbol like "Pa-234m" to an IsotopeKey.
    // Returns IsotopeKey{-1,-1,-1} on failure.
    static IsotopeKey symbolToKey(const std::string& sym);
    static std::string keyToSymbol(const IsotopeKey& k);

private:
    std::string fSource;     // directory or tarball path
    bool        fIsTar = false;
    bool        fIsGzip = false;

    // For directory mode: filename inside fSource.
    // For tar mode: byte offset + length inside the (possibly gunzipped)
    // archive contents stored in fTarBlob.
    struct Entry {
        std::string filename;     // dir mode: filename; tar mode: kept for diagnostics
        size_t      offset = 0;   // tar mode only
        size_t      length = 0;   // tar mode only
    };

    std::unordered_map<IsotopeKey, Entry> fAvailable;
    std::unordered_map<IsotopeKey, ParentDecayInfo> fByKey;
    std::unordered_map<IsotopeKey, bool> fAttempted;

    // For tarball mode: the entire (decompressed) tar contents loaded into
    // memory once at construction. Each Entry's offset/length point into
    // this blob.
    std::string fTarBlob;

    // Build the index. Dispatches on fIsTar.
    void scanDirectory();
    void scanTarball();

    // Decompress a gzip stream into memory.
    static std::string gunzip(const std::string& path);

    // Parse file content (from disk or tar blob) into a ParentDecayInfo.
    // Returns true on success.
    bool parseContent(const IsotopeKey& key, const std::string& content);

    // Read content for an entry into a string.
    std::string readEntry(const Entry& e) const;
};

} // namespace g4gamma
