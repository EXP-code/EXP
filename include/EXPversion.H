#pragma once

#include <iostream>

namespace __EXP__
{
  // Runtime parser for "X.Y.Z" format

  /* Example usage:

     #define VERSION_STR "2.4.12"
     static constexpr EXPversion current_v = EXP_parse_version(VERSION_STR);

     int main() {
       static_assert(current_v.major == 2, "Major version mismatch");
       std::cout << "Parsed: " << current_v.major << "." << current_v.minor << "\n";
     return 0;}
  */


  struct EXPversion
  {
    int major, minor, patch;
  };

  // Compile-time parser for "X.Y.Z" format
  constexpr EXPversion EXP_parse_version(const char* str)
  {
    EXPversion v = {0, 0, 0};
    int* target = &v.major;
    int current = 0;
    
    // Scan the version string until we hit a dot, then move to the next target
    for (int i = 0; str[i] != '\0'; ++i) {
      if (str[i] == '.') {
	*target = current;
	if (target == &v.major) target = &v.minor;
	else if (target == &v.minor) target = &v.patch;
	current = 0;
      } else {
	// Add the current digit to the current version component
	current = current * 10 + (str[i] - '0');
      }
    }
    *target = current;
    return v;
  }
}

