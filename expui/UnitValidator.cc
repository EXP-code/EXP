#include "UnitValidator.H"

// Map aliases to their canonical (primary) unit types
std::unordered_map<std::string, std::string>
UnitValidator::createAllowedUnitTypes()
{
  std::unordered_map<std::string, std::string> allowed;

  // Canonical names
    allowed["length"]   = "length";
    allowed["mass"]     = "mass";
    allowed["time"]     = "time";
    allowed["velocity"] = "velocity";
    allowed["G"]        = "G";

    // Aliases
    allowed["len"]    = "length";
    allowed["l"]      = "length";
    allowed["L"]      = "length";
    allowed["m"]      = "mass";
    allowed["M"]      = "mass";
    allowed["vel"]    = "velocity";
    allowed["Vel"]    = "velocity";
    allowed["v"]      = "velocity";
    allowed["V"]      = "velocity";
    allowed["grav"]   = "G";
    allowed["gravitational_constant"] = "G";

    return allowed;
}

// Map aliases to their canonical (primary) unit names
std::unordered_map<std::string, std::string>
UnitValidator::createAllowedUnitNames()
{
  std::unordered_map<std::string, std::string> allowed;

  // Canonical length units
  //
  allowed["m"]        = "m";
  allowed["cm"]       = "cm";
  allowed["km"]       = "km";
  allowed["um"]       = "um";
  allowed["nm"]       = "nm";
  allowed["Angstrom"] = "Angstrom";
  allowed["AU"]       = "AU";
  allowed["ly"]       = "ly";
  allowed["pc"]       = "pc";
  allowed["kpc"]      = "kpc";
  allowed["Mpc"]      = "Mpc";

  // Astronomical mass units
  //
  allowed["Msun"]     = "Msun";
  allowed["Mearth"]   = "Mearth";
  allowed["g"]        = "g";
  allowed["kg"]       = "kg";
  
  // Time units
  //
  allowed["s"]        = "s";
  allowed["min"]      = "min";
  allowed["hr"]       = "hr";
  allowed["day"]      = "day";
  allowed["yr"]       = "yr";
  allowed["Myr"]      = "Myr";
  allowed["Gyr"]      = "Gyr";
	
  // Velocity units
  //
  allowed["m/s"]      = "m/s";
  allowed["km/s"]     = "km/s";
  allowed["km/hr"]    = "km/hr";
  allowed["km/min"]   = "km/min";
  allowed["c"]        = "c";
  

  // Length aliases
  //
  allowed["meter"]      = "m";
  allowed["centimeter"] = "cm";
  allowed["kilometer"]  = "km";
  allowed["nanometer"]  = "nm";
  allowed["micrometer"] = "um";
  allowed["micron"]     = "um";
  allowed["angstrom"]   = "Angstrom";
  allowed["AA"]         = "Angstrom";
  allowed["astronomical_unit"] = "AU";
  allowed["au"]         = "AU";
  allowed["light_year"] = "ly";
  allowed["lyr"]        = "ly";
  allowed["parsec"]     = "pc";
  allowed["kiloparsec"] = "kpc";
  allowed["megaparsec"] = "Mpc";


  // Mass aliases
  //
  allowed["solar_mass"] = "Msun";
  allowed["earth_mass"] = "Mearth";
  allowed["gram"]       = "g";
  allowed["kilograms"]  = "kg";

  // Time aliases
  //
  allowed["second"]     = "s";
  allowed["minute"]     = "min";
  allowed["hour"]       = "hr";
  allowed["year"]       = "yr";

  // Velocity aliases
  //
  allowed["meter_per_second"] = "m/s";
  allowed["m_per_s"]          = "m/s";
  allowed["km_per_s"]         = "km/s";
  allowed["km_per_hr"]        = "km/hr";
  allowed["km_per_min"]       = "km/min";
  allowed["speed_of_light"]   = "c";

  // Special non-units
  //
  allowed["mixed"]            = "mixed";
  allowed["none"]             = "none";

  return allowed;
}

