// IsotopeKey.hh -- identifies a nuclide including its isomer state
#pragma once
#include <string>
#include <cstdint>
#include <functional>

namespace g4gamma {

// Z: atomic number; A: mass number; M: isomer level (0=ground, 1=first metastable, ...)
struct IsotopeKey {
    int Z;
    int A;
    int M;  // metastable index. 0 for ground state.

    bool operator==(const IsotopeKey& o) const noexcept {
        return Z == o.Z && A == o.A && M == o.M;
    }
    bool operator<(const IsotopeKey& o) const noexcept {
        if (Z != o.Z) return Z < o.Z;
        if (A != o.A) return A < o.A;
        return M < o.M;
    }

    std::string str() const;          // e.g. "Cs137" or "Ba137m"
    static std::string elementSymbol(int Z);
};

} // namespace g4gamma

namespace std {
template <>
struct hash<g4gamma::IsotopeKey> {
    std::size_t operator()(const g4gamma::IsotopeKey& k) const noexcept {
        // Pack (Z<=128, A<=300, M<=10) into a single 32-bit key.
        return std::hash<uint64_t>()(
            (static_cast<uint64_t>(k.Z) << 16) |
            (static_cast<uint64_t>(k.A) << 4)  |
             static_cast<uint64_t>(k.M));
    }
};
} // namespace std
