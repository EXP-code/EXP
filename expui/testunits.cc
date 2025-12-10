// This is a simple test program for the UnitValidator class

#include <iostream>
#include "UnitValidator.H"

int main()
{
  UnitValidator check;

  std::string type, unit;
  std::cout << "Enter type and unit: ";
  std::cin >> type >> unit;

  bool valid;
  std::string canonical_type, canonical_unit;
  std::tie(valid, canonical_type, canonical_unit) = check(type, unit);

  if (valid) {
    std::cout << "The type '" << type << "' with unit '" << unit
	      << "' is valid." << std::endl;
    std::cout << "The canonical names are: Type='" << canonical_type
      	      << "', Unit='" << canonical_unit << "'" << std::endl;
  } else {
    std::cout << "The type '" << type << "' with unit '" << unit
	      << "' is not valid." << std::endl;
  }

  // G test
  std::tie(valid, canonical_type, canonical_unit) = check("G", "");
  if (valid) {
    std::cout << "The type 'G' with units '' is valid." << std::endl;
  } else {
    std::cout << "The type 'G' with units '' is not valid." << std::endl;
  }

  return 0;
}
