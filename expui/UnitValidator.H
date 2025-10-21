#include <unordered_map>
#include <string>

class UnitValidator
{
private:

  //! Allowed unit types and their aliases
  std::unordered_map<std::string, std::string> allowed_types;

  //! Allowed unit names and their aliases
  std::unordered_map<std::string, std::string> allowed_units;

  //! Function to check if an input string is in the allowed list or is
  //! a valid alias
  bool isAllowed(const std::string& input,
		 const std::unordered_map<std::string, std::string>& allowed)
  {
    return allowed.count(input) > 0;
  }

  //! Function to get the canonical name for an input string
  std::string getCanonicalName
  (const std::string& input,
   const std::unordered_map<std::string, std::string>& allowed)
  {
    if (isAllowed(input, allowed)) {
      return allowed.at(input);
    }
    return "unknown";
  }

  std::unordered_map<std::string, std::string> createAllowedUnitTypes();
  std::unordered_map<std::string, std::string> createAllowedUnitNames();

public:

  //! Constructor
  UnitValidator()
  {
    allowed_types = createAllowedUnitTypes();
    allowed_units = createAllowedUnitNames();
  }

  //! Destructor
  ~UnitValidator() {}

  //! Check if type is allowed
  bool allowedType(const std::string& type)
  {
    return isAllowed(type, allowed_types);
  }

  //! Check if type is allowed
  bool allowedUnit(const std::string& unit)
  {
    return isAllowed(unit, allowed_units);
  }

  std::string canonicalType(const std::string& type)
  {
    return getCanonicalName(type, allowed_types);
  }

  std::string canonicalUnit(const std::string& unit)
  {
    return getCanonicalName(unit, allowed_units);
  }

};
