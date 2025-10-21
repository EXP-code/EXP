// This is a simple test program for the UnitValidator class

#include <iostream>
#include "UnitValidator.H"

int main()
{
  UnitValidator check;

  std::string type, unit;
  std::cout << "Enter type and unit: ";
  std::cin >> type >> unit;

  if (check.allowedType(type)) {
    auto canonical = check.canonicalType(type);
    std::cout << "'" << type << "' is allowed. The canonical type is '"
	      << canonical << "'." << std::endl;
  } else {
    std::cout << "'" << type << "' is not a recognized string or alias."
	      << std::endl;
  }
  
  if (check.allowedUnit(unit)) {
    std::string canonical = check.canonicalUnit(unit);
    std::cout << "'" << unit << "' is allowed. The canonical name is '"
	      << canonical << "'." << std::endl;
  } else {
    std::cout << "'" << unit << "' is not a recognized string or alias."
	      << std::endl;
  }
  
  return 0;
}
