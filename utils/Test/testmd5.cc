// This test verifies that QuickDigest5 correctly computes the MD5
// hash of a file

#include <cstdlib>
#include <filesystem>
#include <iostream>
#include <string>
#include "quickdigest5.hpp"

#include "exputils.H" // For get_md5sum which uses the system's md5sum
		      // command for verification

int main(int argc, char* argv[])
{
  // Default to example.txt if no argument
  std::string filePath = argc > 1 ? argv[1] : "example.txt";

  // Check if the file exists before trying to hash it
  std::filesystem::path p(filePath);
  if (!std::filesystem::exists(p)) {
    std::cerr << "File <" << filePath << "> not found." << std::endl
              << "Usage: " << argv[0] << " [file_path]" << std::endl;
    return 1;
  }

  // One-line method to get the hex digest of a file
  std::string hash = QuickDigest5::fileToHash(filePath);

  if (!hash.empty()) {
    std::cout << "MD5: " << hash << std::endl;
  } else {
    std::cerr << "Error: Could not process file." << std::endl;
    return 1;
  }
  
  // System version of md5sum for comparison
  try {
    if (std::system("which md5sum > /dev/null 2>&1") != 0) {
      std::cerr << "Warning: md5sum command not found. Skipping system md5sum comparison." << std::endl;
      return 0;
    }
    std::string systemHash = get_md5sum(filePath);
    std::cout << "System md5sum: " << systemHash << std::endl;
    if (hash == systemHash) {
      std::cout << "Success: hashes match!" << std::endl;
    } else {
      std::cerr << "Error: hashes do not match!" << std::endl;
      return 1;
    }
  } catch (const std::exception& e) {
    std::cerr << "Error computing system md5sum: " << e.what() << std::endl;
    return 1;
  }

  return 0;
}
