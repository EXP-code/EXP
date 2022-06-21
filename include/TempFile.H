#ifndef _TempFile_H
#define _TempFile_H

#include <filesystem>
#include <iostream>
#include <random>

//! Generate a unique file name with a given prefix and 'length'
//! random characters.  The routine will check for existence and
//! attempt 'max_try' times to make unique name.  Otherwise, it throws
//! an exception reporting the last attempted file name.  We use
//! std::filesystem to obtain the path for temporary files.  The value
//! of max_try=1000 by default.
//!
//! Examples:
//!
//! auto path = temp_file("my_prog");      // Uses the default length=16
//! auto path = temp_file("my_prog", 8);   // User specfied length=8
//!
std::filesystem::path temp_file(const std::string& prefix, std::size_t length=16)
{
  // Allowed characters
  //
  const std::string alphanum = "ABCDEFGHIJKLMNOPQRSTUVWXYZabcdefghijklmnopqrstuvwxyz0123456789";

  // Make a random number generator
  //
  std::random_device random_device; // This should be seeded by
				    // /dev/urandom in linux
  std::mt19937 generator(random_device());
  std::uniform_int_distribution<> dist(0, alphanum.size() - 1);

  // Number of tries to find unique file
  //
  const int max_try = 10000;
  int try_num = 0;

  std::filesystem::path path;

  do {
    std::string rs = prefix + "_";
    for (std::size_t i = 0; i < length; ++i) rs += alphanum[dist(generator)];

    path = std::filesystem::temp_directory_path() / rs;
    try_num++;
  } while (std::filesystem::exists(path) or try_num>max_try);

  // Still can't find a unique file??
  //
  if (try_num>max_try) {
    std::ostringstream sout;
    sout << "EXP::TempFile: could not find unique temporary file after "
	 << max_try << " attempts. Last try was: " << path;
    throw std::runtime_error(sout.str());
  }

  // The path to a unique temporary file
  //
  return path;
}

#endif

