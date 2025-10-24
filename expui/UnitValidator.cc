#include "UnitValidator.H"

// Constructor/initializer
UnitValidator::UnitValidator()
{
  // Initialize the 'dictionaries'
  allowed_types = createAllowedUnitTypes();
  allowed_units = createAllowedUnitNames();
}

// Do the full check and return canonical strings in two steps.  First
// check the type.  Then the unit within that type.  Returns a tuple of
// (is_valid, canonical_type, canonical_unit).
std::tuple<bool, std::string, std::string>
UnitValidator::operator()(const std::string& type, const std::string& unit)
{
  // Check type first
  if (allowed_types.count(type) > 0) {

    // Get canonical type
    std::string canonical_type = allowed_types.at(type);
    
    // Now check unit in the type category
    if (allowed_units[canonical_type].count(unit) > 0) {

      // Get canonical unit name
      std::string canonical_unit = allowed_units[canonical_type].at(unit);

      // Return successful final results
      return {true, canonical_type, canonical_unit};
    }
  }

  // If we get here, we have a type or unit that is not recognized
  return {false, "unknown", "unknown"};
}


// Map aliases to their canonical (primary) unit types
std::unordered_map<std::string, std::string>
UnitValidator::createAllowedUnitTypes()
{
  std::unordered_map<std::string, std::string> allowed;

  // These are the canonical names
  //
  allowed["length"]   = "length";
  allowed["mass"]     = "mass";
  allowed["time"]     = "time";
  allowed["velocity"] = "velocity";
  allowed["G"]        = "G";
  
  // These are recognized aliases aliases for the types
  //
  allowed["Length"]   = "length";
  allowed["Len"]      = "length";
  allowed["len"]      = "length";
  allowed["l"]        = "length";
  allowed["L"]        = "length";

  allowed["Mass"]     = "mass";
  allowed["m"]        = "mass";
  allowed["M"]        = "mass";

  allowed["Time"]     = "time";
  allowed["t"]        = "time";
  allowed["T"]        = "time";

  allowed["vel"]      = "velocity";
  allowed["Vel"]      = "velocity";
  allowed["Velocity"] = "velocity";
  allowed["v"]        = "velocity";
  allowed["V"]        = "velocity";
  
  allowed["Grav"]     = "G";
  allowed["grav"]     = "G";
  allowed["grav_constant"] = "G";
  allowed["Grav_constant"] = "G";
  allowed["gravitational_constant"] = "G";
  allowed["Gravitational_constant"] = "G";
  
  return allowed;
}

// Map aliases to their canonical (primary) unit names for each type
// category
std::map<std::string, std::unordered_map<std::string, std::string>>
UnitValidator::createAllowedUnitNames()
{
  std::map<std::string, std::unordered_map<std::string, std::string>> allowed;

  // Allow 'none' as a unit for all types
  //
  allowed["length"]["none"]     = "none";
  allowed["mass"]["none"]       = "none";
  allowed["time"]["none"]       = "none";
  allowed["length"]["None"]     = "none";
  allowed["mass"]["None"]       = "none";
  allowed["velocity"]["none"]   = "none";

  // Special non-units
  //
  allowed["G"][""]              = "none";
  allowed["G"]["mixed"]         = "mixed";
  allowed["G"]["none"]          = "none";
  allowed["G"]["unitless"]      = "none";

  // Canonical length units
  //
  allowed["length"]["m"]        = "m";
  allowed["length"]["cm"]       = "cm";
  allowed["length"]["km"]       = "km";
  allowed["length"]["um"]       = "um";
  allowed["length"]["nm"]       = "nm";
  allowed["length"]["Angstrom"] = "Angstrom";
  allowed["length"]["AU"]       = "AU";
  allowed["length"]["ly"]       = "ly";
  allowed["length"]["pc"]       = "pc";
  allowed["length"]["kpc"]      = "kpc";
  allowed["length"]["Mpc"]      = "Mpc";

  // Standard astronomical mass units
  //
  allowed["mass"]["Msun"]       = "Msun";
  allowed["mass"]["Mearth"]     = "Mearth";
  allowed["mass"]["g"]          = "g";
  allowed["mass"]["kg"]         = "kg";
  
  // Time units
  //
  allowed["time"]["s"]          = "s";
  allowed["time"]["min"]        = "min";
  allowed["time"]["hr"]         = "hr";
  allowed["time"]["day"]        = "day";
  allowed["time"]["yr"]         = "yr";
  allowed["time"]["Myr"]        = "Myr";
  allowed["time"]["Gyr"]        = "Gyr";
	
  // Velocity units
  //
  allowed["velocity"]["cm/s"]     = "cmm/s";
  allowed["velocity"]["m/s"]      = "m/s";
  allowed["velocity"]["km/s"]     = "km/s";
  allowed["velocity"]["km/hr"]    = "km/hr";
  allowed["velocity"]["km/min"]   = "km/min";
  allowed["velocity"]["c"]        = "c";
  allowed["velocity"]["none"]     = "none";
  

  // Length aliases
  //
  allowed["length"]["meter"]      = "m";
  allowed["length"]["centimeter"] = "cm";
  allowed["length"]["kilometer"]  = "km";
  allowed["length"]["nanometer"]  = "nm";
  allowed["length"]["micrometer"] = "um";
  allowed["length"]["micron"]     = "um";
  allowed["length"]["angstrom"]   = "Angstrom";
  allowed["length"]["AA"]         = "Angstrom";
  allowed["length"]["astronomical_unit"] = "AU";
  allowed["length"]["au"]         = "AU";
  allowed["length"]["light_year"] = "ly";
  allowed["length"]["lyr"]        = "ly";
  allowed["length"]["parsec"]     = "pc";
  allowed["length"]["kiloparsec"] = "kpc";
  allowed["length"]["megaparsec"] = "Mpc";


  // Mass aliases
  //
  allowed["mass"]["solar_mass"] = "Msun";
  allowed["mass"]["earth_mass"] = "Mearth";
  allowed["mass"]["gram"]       = "g";
  allowed["mass"]["kilograms"]  = "kg";

  // Time aliases
  //
  allowed["time"]["second"]     = "s";
  allowed["time"]["minute"]     = "min";
  allowed["time"]["hour"]       = "hr";
  allowed["time"]["year"]       = "yr";
  allowed["time"]["None"]       = "none";

  // Velocity aliases
  //
  allowed["velocity"]["meter_per_second"]      = "m/s";
  allowed["velocity"]["centimeter_per_second"] = "cm/s";
  allowed["velocity"]["cm_per_s"]              = "cm/s";
  allowed["velocity"]["m_per_s"]               = "m/s";
  allowed["velocity"]["km_per_s"]              = "km/s";
  allowed["velocity"]["km_per_hr"]             = "km/hr";
  allowed["velocity"]["km_per_min"]            = "km/min";
  allowed["velocity"]["speed_of_light"]        = "c";


  return allowed;
}

