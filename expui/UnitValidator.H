#include <unordered_map>
#include <string>
#include <map>

/** Class to validate unit types and names The allowed types and units
    can be extended as needed by adding new names and aliases to
    UnitValidator::createAllowedUnitTypes() and
    UnitValidator::createAllowedUnitNames() in UnitValidator.cc

    Strategy:
    -------
    - Types and their aliases are stored in an std::unordered_map
    - Units and their aliases are stored in an std::map of
      std::unordered_maps
    - The operator()(type, unit, value) first checks if the provided
      type exists in the allowed types map.  If it does, it retrieves
      the canonical type string which is used to look up the allowed
      units for that type in the std::map.  Then, the unit type is
      checked against the allowed units for that type.  If both type
      and unit are valid, the canonical strings are returned.

    Usage:
    -----
    - Create an instance of UnitValidator.  Currently, this is done in
      the Coefficients class by default.
    - Call the operator() on an instance with the type and unit
      strings to be validated
    - The operator() returns a tuple with:
       + a boolean indicating if the type and unit are valid
       + the canonical type string, or "unknown" if invalid
       + the canonical unit string, or "unknown" if invalid
*/
class UnitValidator
{
private:

  //! Store allowed unit types and their aliases
  std::unordered_map<std::string, std::string> allowed_types;

  //! Store allowed unit names and their aliases
  std::map<std::string, std::unordered_map<std::string, std::string>> allowed_units;

  //! Create the allowed types map (length, mass, time, velocity, G
  //! [and aliases])
  std::unordered_map<std::string, std::string> createAllowedUnitTypes();

  //! Create the allowed units map
  std::map<std::string, std::unordered_map<std::string, std::string>> createAllowedUnitNames();

public:

  //! Constructor
  UnitValidator();

  //! Destructor
  ~UnitValidator() {}

  //! Do the full check and return canonical strings
  std::tuple<bool, std::string, std::string>
  operator()(const std::string& type, const std::string& unit);

};
