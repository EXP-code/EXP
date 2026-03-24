// This test verifies that QuickDigest5 correctly computes the MD5
// hash of a file

#include <iostream>
#include <filesystem>

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
    std::cout << "File <" << filePath << "> not found." << std::endl
	      << "Usage: " << argv[0] << " [file_path]" << std::endl
	      << "Defaulting to example.txt if it exists." << std::endl;
  }

  // One-line method to get the hex digest of a file
  std::string hash = QuickDigest5::fileToHash(filePath);

  if (!hash.empty()) {
    std::cout << "MD5: " << hash << std::endl;
  } else {
    std::cerr << "Error: Could not process file." << std::endl;
  }
  
  // System version of md5sum for comparison
  try {
    std::string systemHash = get_md5sum(filePath);
    std::cout << "System md5sum: " << systemHash << std::endl;
    if (hash == systemHash) {
      std::cout << "Success: hashes match!" << std::endl;
    } else {
      std::cerr << "Error: hashes do not match!" << std::endl;
    }
  } catch (const std::exception& e) {
    std::cerr << "Error computing system md5sum: " << e.what() << std::endl;
  }

  return 0;
}
