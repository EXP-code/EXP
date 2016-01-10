#include "atomic_constants.H"
#include <iomanip>

// Populate the periodic table with data.  Currently we only have
// element names and atomic numbers.
//
PeriodicTable::PeriodicTable()
{
  add("Hydrogen",        "H",      1,      1.008);
  add("Helium",          "He",     2,   4.002602);
  add("Lithium",         "Li",     3,       6.94);
  add("Beryllium",       "Be",     4,  9.0121831);
  add("Boron",           "B",      5,      10.81);
  add("Carbon",          "C",      6,     12.011);
  add("Nitrogen",        "N",      7,     14.007);
  add("Oxygen",          "O",      8,     15.999);
  add("Fluorine",        "F",      9,  18.998403);
  add("Neon",            "Ne",    10,    20.1797);
  add("Sodium",          "Na",    11,  22.989769);
  add("Magnesium",       "Mg",    12,     24.305);
  add("Aluminium",       "Al",    13,  26.981538);
  add("Silicon",         "Si",    14,     28.085);
  add("Phosphorus",      "P",     15,  30.973762);
  add("Sulfur",          "S",     16,      32.06);
  add("Chlorine",        "Cl",    17,      35.45);
  add("Argon",           "Ar",    18,     39.948);
  add("Potassium",       "K",     19,    39.0983);
  add("Calcium",         "Ca",    20,     40.078);
  add("Scandium",        "Sc",    21,  44.955908);
  add("Titanium",        "Ti",    22,     47.867);
  add("Vanadium",        "V",     23,    50.9415);
  add("Chromium",        "Cr",    24,    51.9961);
  add("Manganese",       "Mn",    25,  54.938044);
  add("Iron",            "Fe",    26,     55.845);
  add("Cobalt",          "Co",    27,  58.933194);
  add("Nickel",          "Ni",    28,    58.6934);
  add("Copper",          "Cu",    29,     63.546);
  add("Zinc",            "Zn",    30,      65.38);
  add("Gallium",         "Ga",    31,     69.723);
  add("Germanium",       "Ge",    32,      72.63);
  add("Arsenic",         "As",    33,  74.921595);
  add("Selenium",        "Se",    34,     78.971);
  add("Bromine",         "Br",    35,     79.904);
  add("Krypton",         "Kr",    36,     83.798);
  add("Rubidium",        "Rb",    37,    85.4678);
  add("Strontium",       "Sr",    38,      87.62);
  add("Yttrium",         "Y",     39,   88.90584);
  add("Zirconium",       "Zr",    40,     91.224);
  add("Niobium",         "Nb",    41,   92.90637);
  add("Molybdenum",      "Mo",    42,      95.95);
  add("Technetium",      "Tc",    43,         97);
  add("Ruthenium",       "Ru",    44,     101.07);
  add("Rhodium",         "Rh",    45,   102.9055);
  add("Palladium",       "Pd",    46,     106.42);
  add("Silver",          "Ag",    47,   107.8682);
  add("Cadmium",         "Cd",    48,    112.414);
  add("Indium",          "In",    49,    114.818);
  add("Tin",             "Sn",    50,     118.71);
  add("Antimony",        "Sb",    51,     121.76);
  add("Tellurium",       "Te",    52,      127.6);
  add("Iodine",          "I",     53,  126.90447);
  add("Xenon",           "Xe",    54,    131.293);
  add("Caesium",         "Cs",    55,  132.90545);
  add("Barium",          "Ba",    56,    137.327);
  add("Lanthanum",       "La",    57,  138.90547);
  add("Cerium",          "Ce",    58,    140.116);
  add("Praseodymium",    "Pr",    59,  140.90766);
  add("Neodymium",       "Nd",    60,    144.242);
  add("Promethium",      "Pm",    61,        145);
  add("Samarium",        "Sm",    62,     150.36);
  add("Europium",        "Eu",    63,    151.964);
  add("Gadolinium",      "Gd",    64,     157.25);
  add("Terbium",         "Tb",    65,  158.92535);
  add("Dysprosium",      "Dy",    66,      162.5);
  add("Holmium",         "Ho",    67,  164.93033);
  add("Erbium",          "Er",    68,    167.259);
  add("Thulium",         "Tm",    69,  168.93422);
  add("Ytterbium",       "Yb",    70,    173.045);
  add("Lutetium",        "Lu",    71,   174.9668);
  add("Hafnium",         "Hf",    72,     178.49);
  add("Tantalum",        "Ta",    73,  180.94788);
  add("Tungsten",        "W",     74,     183.84);
  add("Rhenium",         "Re",    75,    186.207);
  add("Osmium",          "Os",    76,     190.23);
  add("Iridium",         "Ir",    77,    192.217);
  add("Platinum",        "Pt",    78,    195.084);
  add("Gold",            "Au",    79,  196.96657);
  add("Mercury",         "Hg",    80,    200.592);
  add("Thallium",        "Tl",    81,     204.38);
  add("Lead",            "Pb",    82,      207.2);
  add("Bismuth",         "Bi",    83,   208.9804);
  add("Polonium",        "Po",    84,        209);
  add("Astatine",        "At",    85,        210);
  add("Radon",           "Rn",    86,        222);
  add("Francium",        "Fr",    87,        223);
  add("Radium",          "Ra",    88,        226);
  add("Actinium",        "Ac",    89,        227);
  add("Thorium",         "Th",    90,   232.0377);
  add("Protactinium",    "Pa",    91,  231.03588);
  add("Uranium",         "U",     92,  238.02891);
  add("Neptunium",       "Np",    93,        237);
  add("Plutonium",       "Pu",    94,        244);
  add("Americium",       "Am",    95,        243);
  add("Curium",          "Cm",    96,        247);
  add("Berkelium",       "Bk",    97,        247);
  add("Californium",     "Cf",    98,        251);
  add("Einsteinium",     "Es",    99,        252);
  add("Fermium",         "Fm",   100,        257);
  add("Mendelevium",     "Md",   101,        258);
  add("Nobelium",        "No",   102,        259);
  add("Lawrencium",      "Lr",   103,        262);
  add("Rutherfordium",   "Rf",   104,        267);
  add("Dubnium",         "Db",   105,        270);
  add("Seaborgium",      "Sg",   106,        271);
  add("Bohrium",         "Bh",   107,        270);
  add("Hassium",         "Hs",   108,        277);
  add("Meitnerium",      "Mt",   109,        276);
  add("Darmstadtium",    "Ds",   110,        281);
  add("Roentgenium",     "Rg",   111,        282);
  add("Copernicium",     "Cn",   112,        285);
  add("Ununtrium",       "Uut",  113,        285);
  add("Flerovium",       "Fl",   114,        289);
  add("Ununpentium",     "Uup",  115,        289);
  add("Livermorium",     "Lv",   116,        293);
  add("Ununseptium",     "Uus",  117,        294);
  add("Ununoctium",      "Uuo",  118,        294);
}

void PeriodicTable::print(std::ostream& out)
{
  out << std::right
      << std::setw(6 ) << "Z"
      << std::setw(20) << "Element"
      << std::setw(6 ) << "Abbrev"
      << std::setw(20) << "Weight"
      << std::setw(20) << "Scale"
      << std::endl;

  out << std::setw(6 ) << "-----"
      << std::setw(20) << "--------"
      << std::setw(6 ) << "-----"
      << std::setw(20) << "--------"
      << std::setw(20) << "--------"
      << std::endl;

  for (auto v : dataZ)
    out << std::setw(6 ) << std::get<2>(*v.second)
	<< std::setw(20) << std::get<0>(*v.second)
	<< std::setw(6 ) << std::get<1>(*v.second)
	<< std::setw(20) << std::get<3>(*v.second)
	<< std::setw(20) << std::get<4>(*v.second)
	<< std::endl;

  out << std::endl;
}
