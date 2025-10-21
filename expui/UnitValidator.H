#include <unordered_map>
#include <string>
#include <map>

/** Class to validate unit types and names The allowed types and units
    can be extended as needed by adding new names and aliasess to
    UnitValidator::createAllowedUnitTypes() and
    UnitValidator::createAllowedUnitNames() in UnitValidator.cc
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
