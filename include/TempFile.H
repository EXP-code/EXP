#ifndef _TempFile_H
#define _TempFile_H

#include <filesystem>
#include <iostream>
#include <random>

std::filesystem::path temp_file(const std::string& prefix, std::size_t length)
{
  // Allowed characters
  //
  const std::string alphanum = "ABCDEFGHIJKLMNOPQRSTUVWXYZabcdefghijklmnopqrstuvwxyz0123456789";

  // Make a random number generator
  //
  std::random_device random_device;
  std::mt19937 generator(random_device());
  std::uniform_int_distribution<> dist(0, alphanum.size() - 1);

  // Number of tries to find unique file
  //
  const int max_trials = 10000;

  std::filesystem::path path;

  for (int i=0; i<max_trials; i++) {

    std::string rs;
    for (std::size_t i = 0; i < length; ++i) {
      rs += alphanum[dist(generator)];
    }

    std::string fname = prefix + "_" + rs;
    path = std::filesystem::temp_directory_path() / fname;

    if (not std::filesystem::exists(path)) break;
  }

  return path;
}

#endif

