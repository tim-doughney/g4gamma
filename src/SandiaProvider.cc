#include "g4gamma/SandiaProvider.hh"
#include "g4gamma/Units.hh"
#include <fstream>
#include <sstream>
#include <stdexcept>
#include <cstdlib>
#include <filesystem>
#include <iostream>
#include <cmath>
#include <zlib.h>

namespace fs = std::filesystem;

namespace g4gamma {

namespace {

// Map sandia.decay mode strings to our DecayMode enum.
DecayMode parseMode(const std::string& s) {
    if (s == "b-")    return DecayMode::BetaMinus;
    if (s == "b+")    return DecayMode::BetaPlus;
    if (s == "ec")    return DecayMode::EC;
    if (s == "it")    return DecayMode::IT;
    if (s == "a")     return DecayMode::Alpha;
    if (s == "p")     return DecayMode::Proton;
    if (s == "n")     return DecayMode::Neutron;        // not seen in file but for completeness
    if (s == "sf")    return DecayMode::SpFission;
    if (s == "2b-")   return DecayMode::Beta2Minus;
    if (s == "2b+")   return DecayMode::Beta2Plus;
    if (s == "2ec")   return DecayMode::EC;             // double EC -> still EC family
    if (s == "2p")    return DecayMode::Proton2;
    if (s == "b+p")   return DecayMode::BDProton;       // beta-delayed proton
    if (s == "b-n")   return DecayMode::BDNeutron;
    if (s == "b-2n")  return DecayMode::BDNeutron;
    if (s == "b-a")   return DecayMode::BetaMinus;      // β-α: handle as β-
    if (s == "b+a")   return DecayMode::BetaPlus;
    if (s == "b+2p")  return DecayMode::BDProton;
    if (s == "b+3p")  return DecayMode::BDProton;
    if (s == "ec2p")  return DecayMode::EC;
    if (s == "eca")   return DecayMode::EC;
    if (s == "ecp")   return DecayMode::EC;
    if (s == "14c")   return DecayMode::Unknown;        // cluster decay, not gamma-relevant
    return DecayMode::Unknown;
}

// Tiny stream parser for the SandiaDecay XML format. We don't need a real
// XML parser -- the file is regular and we just need to find tags by name
// and pull their attributes.
struct AttrParser {
    static std::string get(const std::string& attrs, const char* key) {
        // Scan for `key="..."` allowing for prior attribute values that may
        // contain the substring `key`. We require the previous character to
        // be whitespace or the start of the string.
        std::string needle = std::string(key) + "=\"";
        size_t p = 0;
        while ((p = attrs.find(needle, p)) != std::string::npos) {
            bool atBoundary = (p == 0) ||
                              attrs[p-1] == ' ' || attrs[p-1] == '\t' ||
                              attrs[p-1] == '<';
            if (atBoundary) {
                size_t q = p + needle.size();
                size_t e = attrs.find('"', q);
                if (e == std::string::npos) return {};
                return attrs.substr(q, e - q);
            }
            p += needle.size();
        }
        return {};
    }
    // Try a list of attribute keys, return the first one that resolves.
    static std::string getAny(const std::string& attrs,
                               std::initializer_list<const char*> keys) {
        for (auto* k : keys) {
            auto v = get(attrs, k);
            if (!v.empty()) return v;
        }
        return {};
    }
    static double getDouble(const std::string& attrs, const char* key, double def = 0.0) {
        return parseDouble(get(attrs, key), def);
    }
    static double getDoubleAny(const std::string& attrs,
                                std::initializer_list<const char*> keys, double def = 0.0) {
        return parseDouble(getAny(attrs, keys), def);
    }
    static int getInt(const std::string& attrs, const char* key, int def = 0) {
        return parseInt(get(attrs, key), def);
    }
    static int getIntAny(const std::string& attrs,
                          std::initializer_list<const char*> keys, int def = 0) {
        return parseInt(getAny(attrs, keys), def);
    }
private:
    static double parseDouble(const std::string& s, double def) {
        if (s.empty()) return def;
        if (s == "INF" || s == "inf") return std::numeric_limits<double>::infinity();
        try { return std::stod(s); } catch (...) { return def; }
    }
    static int parseInt(const std::string& s, int def) {
        if (s.empty()) return def;
        try { return std::stoi(s); } catch (...) { return def; }
    }
};

// Detect whether the file uses the verbose tag names (<nuclide>, <transition>,
// <gamma>, <xray>) or the minified ones (<n>, <t>, <g>, <x>). Return a small
// struct of which tag names to look for.
struct TagNames {
    std::string nuclide;
    std::string transition;
    std::string gamma;
    std::string xray;
    // Attribute name lists for getAny -- in priority order (verbose first).
    std::initializer_list<const char*> attrZ        = {"atomicNumber", "an"};
    std::initializer_list<const char*> attrA        = {"massNumber", "mn"};
    std::initializer_list<const char*> attrM        = {"isomerNumber", "iso"};
    std::initializer_list<const char*> attrSymbol   = {"symbol", "s"};
    std::initializer_list<const char*> attrHalfLife = {"halfLife", "hl"};
    std::initializer_list<const char*> attrBR       = {"branchRatio", "br"};
    std::initializer_list<const char*> attrParent   = {"parent", "p"};
    std::initializer_list<const char*> attrChild    = {"child", "c"};
    std::initializer_list<const char*> attrMode     = {"mode", "m"};
    std::initializer_list<const char*> attrEnergy   = {"energy", "e"};
    std::initializer_list<const char*> attrIntensity= {"intensity", "i"};
};

TagNames detectFormat(const std::string& xml) {
    TagNames t;
    // The verbose form uses <nuclide; the minified uses <n. Pick by which one
    // appears more often.
    auto cnt = [&](const char* s) {
        size_t pos = 0, c = 0;
        std::string needle(s);
        while ((pos = xml.find(needle, pos)) != std::string::npos) { ++c; pos += needle.size(); }
        return c;
    };
    bool verbose = cnt("<nuclide ") > cnt("<n ");
    if (verbose) {
        t.nuclide    = "<nuclide";
        t.transition = "<transition";
        t.gamma      = "<gamma";
        t.xray       = "<xray";
    } else {
        t.nuclide    = "<n ";       // trailing space avoids matching <ns or <nuclide
        t.transition = "<t ";       // matches "<t " and "<t>" via the find-then-end-find logic
        t.gamma      = "<g ";
        t.xray       = "<x ";
    }
    return t;
}

} // namespace


SandiaProvider::SandiaProvider(const std::string& hint) {
    fXmlPath = locateXml(hint);
    if (fXmlPath.empty()) {
        throw std::runtime_error(
            "SandiaProvider: could not locate sandia.decay.xml. "
            "Pass an explicit path, set $SANDIA_DECAY_XML, or place the file "
            "in <repo>/data/sandia/.");
    }
    load();
}

std::string SandiaProvider::locateXml(const std::string& hint) {
    auto exists = [](const std::string& p) {
        return !p.empty() && fs::exists(p) && fs::is_regular_file(p);
    };

    // 1. Explicit hint
    if (exists(hint)) return hint;

    // 2. Env var
    if (const char* env = std::getenv("SANDIA_DECAY_XML")) {
        std::string p = env;
        if (exists(p)) return p;
    }

    // 3. Search relative to a few likely locations. We don't know our own
    // install root from C++ at runtime, so try CWD-relative and common
    // installation paths.
    std::vector<std::string> candidates;
    // Names to try (in order from preferred to fallback). We try .xml.gz
    // variants first since the gzipped file is ~5x smaller and the C++
    // gzip read is essentially free.
    static const std::vector<std::string> filenames = {
        "sandia.decay.nocoinc.min.xml.gz",
        "sandia.decay.min.xml.gz",
        "sandia.decay.xml.gz",
        "sandia.decay.nocoinc.min.xml",
        "sandia.decay.min.xml",
        "sandia.decay.xml",
    };
    static const std::vector<std::string> roots = {
        "data/sandia",
        "../data/sandia",
        "../../data/sandia",
        "/usr/local/share/g4gamma",
        "/usr/share/g4gamma",
    };
    for (const auto& r : roots) {
        for (const auto& f : filenames) {
            std::string p = r + "/" + f;
            if (exists(p)) return p;
        }
    }
    return {};
}

void SandiaProvider::load() {
    if (fLoaded) return;
    parseFile(fXmlPath);
    fLoaded = true;
}

void SandiaProvider::parseFile(const std::string& path) {
    using namespace units;

    std::string xml;
    bool isGz = (path.size() >= 3 && path.substr(path.size()-3) == ".gz");
    if (isGz) {
        gzFile gf = gzopen(path.c_str(), "rb");
        if (!gf) throw std::runtime_error("SandiaProvider: cannot open " + path);
        char buf[64 * 1024];
        int n;
        while ((n = gzread(gf, buf, sizeof(buf))) > 0) {
            xml.append(buf, n);
        }
        gzclose(gf);
    } else {
        std::ifstream in(path);
        if (!in) throw std::runtime_error("SandiaProvider: cannot open " + path);
        std::ostringstream buf;
        buf << in.rdbuf();
        xml = buf.str();
    }

    TagNames T = detectFormat(xml);

    // Pass 1: index <nuclide ... /> entries (or <n ... /> in minified form).
    {
        size_t pos = 0;
        while (true) {
            pos = xml.find(T.nuclide, pos);
            if (pos == std::string::npos) break;
            size_t end = xml.find('>', pos);
            if (end == std::string::npos) break;
            std::string tag = xml.substr(pos, end - pos + 1);
            int Z = AttrParser::getIntAny(tag, T.attrZ, -1);
            int A = AttrParser::getIntAny(tag, T.attrA, -1);
            int M = AttrParser::getIntAny(tag, T.attrM, 0);
            std::string sym = AttrParser::getAny(tag, T.attrSymbol);
            double hl = AttrParser::getDoubleAny(tag, T.attrHalfLife, 0.0);
            if (Z >= 0 && A > 0 && !sym.empty()) {
                IsotopeKey key{Z, A, M};
                fSymbolToKey[sym] = key;
                fKeyToSymbol[key] = sym;
                ParentDecayInfo p;
                p.isotope = key;
                if (std::isinf(hl)) {
                    p.stable = true;
                    p.meanLife = 0.0;
                } else if (hl > 0.0) {
                    p.stable = false;
                    p.meanLife = (hl * second) / std::log(2.0);
                } else {
                    p.stable = true;
                    p.meanLife = 0.0;
                }
                fByKey.emplace(key, std::move(p));
            }
            pos = end + 1;
        }
    }

    // Pass 2: walk every transition block. For verbose form: <transition ...>
    // ... </transition>. For minified form: <t ...> ... </t>.
    // Closer name comes from T.transition with the angle bracket / space stripped:
    std::string closer;
    if (T.transition == "<transition") closer = "</transition>";
    else if (T.transition == "<t ")     closer = "</t>";
    else throw std::runtime_error("SandiaProvider: unrecognized format");

    {
        size_t pos = 0;
        while (true) {
            pos = xml.find(T.transition, pos);
            if (pos == std::string::npos) break;
            size_t headEnd = xml.find('>', pos);
            if (headEnd == std::string::npos) break;
            bool selfClose = (headEnd > 0 && xml[headEnd - 1] == '/');
            std::string head = xml.substr(pos, headEnd - pos + 1);

            std::string parent  = AttrParser::getAny(head, T.attrParent);
            std::string child   = AttrParser::getAny(head, T.attrChild);
            std::string modeStr = AttrParser::getAny(head, T.attrMode);
            double br           = AttrParser::getDoubleAny(head, T.attrBR, 1.0);

            DecayBranch branch;
            branch.mode           = parseMode(modeStr);
            branch.branchingRatio = br;

            size_t bodyStart = headEnd + 1;
            size_t bodyEnd;
            if (selfClose) {
                bodyEnd = bodyStart;
            } else {
                bodyEnd = xml.find(closer, bodyStart);
                if (bodyEnd == std::string::npos) break;
            }
            std::string body = xml.substr(bodyStart, bodyEnd - bodyStart);

            // Find emissions in body: gammas and xrays.
            auto findEmissions = [&](const std::string& tagOpener,
                                       EmissionType etype) {
                size_t p = 0;
                while ((p = body.find(tagOpener, p)) != std::string::npos) {
                    size_t e = body.find('>', p);
                    if (e == std::string::npos) break;
                    std::string tag = body.substr(p, e - p + 1);
                    double E = AttrParser::getDoubleAny(tag, T.attrEnergy, 0.0);
                    double I = AttrParser::getDoubleAny(tag, T.attrIntensity, 0.0);
                    if (E > 0.0 && I > 0.0) {
                        Emission em;
                        em.type      = etype;
                        em.energy    = E * keV;
                        em.intensity = I;
                        branch.emissions.push_back(em);
                    }
                    p = e + 1;
                }
            };
            findEmissions(T.gamma, EmissionType::Gamma);
            findEmissions(T.xray,  EmissionType::XRay);

            // Sandia XML labels EC/β+ mixed transitions as mode="ec" and stores
            // the β+ component as <positron> elements.  Read them and add the
            // two 511 keV annihilation photons explicitly so the binner doesn't
            // need to infer β+ from the mode string.
            {
                double betaPlusIntensity = 0.0;
                size_t pp = 0;
                const std::string positronTag = "<positron";
                while ((pp = body.find(positronTag, pp)) != std::string::npos) {
                    size_t e = body.find('>', pp);
                    if (e == std::string::npos) break;
                    std::string tag = body.substr(pp, e - pp + 1);
                    double I = AttrParser::getDouble(tag, "intensity", 0.0);
                    if (I == 0.0) I = AttrParser::getDouble(tag, "i", 0.0);
                    betaPlusIntensity += I;
                    pp = e + 1;
                }
                if (betaPlusIntensity > 0.0) {
                    Emission em;
                    em.type      = EmissionType::AnnihilationPair;
                    em.energy    = units::electron_mass_c2;
                    em.intensity = 2.0 * betaPlusIntensity;
                    branch.emissions.push_back(em);
                }
            }

            // Resolve parent/daughter symbols
            auto pIt = fSymbolToKey.find(parent);
            auto cIt = fSymbolToKey.find(child);
            if (pIt != fSymbolToKey.end()) {
                if (cIt != fSymbolToKey.end()) {
                    branch.daughter = cIt->second;
                } else {
                    branch.daughter = IsotopeKey{0, 0, 0};
                }
                auto& info = fByKey[pIt->second];
                info.stable = false;
                info.branches.push_back(std::move(branch));
            }

            pos = selfClose ? headEnd + 1 : (bodyEnd + closer.size());
        }
    }

    // Pass 3: tidy.
    for (auto& [key, info] : fByKey) {
        if (info.branches.empty()) info.stable = true;
    }
}

const ParentDecayInfo* SandiaProvider::get(const IsotopeKey& key) {
    auto it = fByKey.find(key);
    if (it == fByKey.end()) return nullptr;
    return &it->second;
}

} // namespace g4gamma
