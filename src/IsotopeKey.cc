#include "g4gamma/IsotopeKey.hh"

namespace g4gamma {

static const char* kElementSymbols[] = {
    "n",  "H",  "He", "Li", "Be", "B",  "C",  "N",  "O",  "F",
    "Ne", "Na", "Mg", "Al", "Si", "P",  "S",  "Cl", "Ar", "K",
    "Ca", "Sc", "Ti", "V",  "Cr", "Mn", "Fe", "Co", "Ni", "Cu",
    "Zn", "Ga", "Ge", "As", "Se", "Br", "Kr", "Rb", "Sr", "Y",
    "Zr", "Nb", "Mo", "Tc", "Ru", "Rh", "Pd", "Ag", "Cd", "In",
    "Sn", "Sb", "Te", "I",  "Xe", "Cs", "Ba", "La", "Ce", "Pr",
    "Nd", "Pm", "Sm", "Eu", "Gd", "Tb", "Dy", "Ho", "Er", "Tm",
    "Yb", "Lu", "Hf", "Ta", "W",  "Re", "Os", "Ir", "Pt", "Au",
    "Hg", "Tl", "Pb", "Bi", "Po", "At", "Rn", "Fr", "Ra", "Ac",
    "Th", "Pa", "U",  "Np", "Pu", "Am", "Cm", "Bk", "Cf", "Es",
    "Fm", "Md", "No", "Lr", "Rf", "Db", "Sg", "Bh", "Hs", "Mt",
    "Ds", "Rg", "Cn", "Nh", "Fl", "Mc", "Lv", "Ts", "Og"
};
static constexpr int kElementCount =
    sizeof(kElementSymbols) / sizeof(kElementSymbols[0]);

std::string IsotopeKey::elementSymbol(int Z) {
    if (Z >= 0 && Z < kElementCount) return kElementSymbols[Z];
    return "Z" + std::to_string(Z);
}

std::string IsotopeKey::str() const {
    std::string s = elementSymbol(Z) + std::to_string(A);
    if (M == 1) s += "m";
    else if (M >= 2) s += "m" + std::to_string(M);
    return s;
}

} // namespace g4gamma
