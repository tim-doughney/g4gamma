// Units.hh -- minimal subset of the Geant4/CLHEP system of units, in the
// same internal-unit convention as Geant4 (energy in MeV, time in nanoseconds).
// Anything you multiply by these constants becomes "internal" and can be
// safely combined with other internal-unit values.
//
// To convert a value from internal units to user-friendly units, divide:
//     double e_keV = energy_internal / keV;
#pragma once

namespace g4gamma {
namespace units {

// --- Energy (internal unit: MeV) ---
constexpr double MeV = 1.0;
constexpr double eV  = 1.0e-6 * MeV;
constexpr double keV = 1.0e-3 * MeV;
constexpr double GeV = 1.0e3  * MeV;

// --- Time (internal unit: nanosecond) ---
constexpr double ns           = 1.0;
constexpr double picosecond   = 1.0e-3 * ns;
constexpr double microsecond  = 1.0e3  * ns;
constexpr double millisecond  = 1.0e6  * ns;
constexpr double second       = 1.0e9  * ns;
constexpr double minute       = 60.0    * second;
constexpr double hour         = 3600.0  * second;
constexpr double day          = 86400.0 * second;
constexpr double year         = 365.25  * day;

// Aliases to match common Geant4 symbol names
constexpr double s   = second;
constexpr double ms  = millisecond;
constexpr double us  = microsecond;
constexpr double ps  = picosecond;
constexpr double min_ = minute;       // 'min' clashes with std::min
constexpr double h    = hour;
constexpr double d    = day;
constexpr double y    = year;

// Annihilation gamma energy
constexpr double electron_mass_c2 = 0.510998950 * MeV;

} // namespace units
} // namespace g4gamma
