#include <string>
#include <cstdio>
#include <array>
#include <stdexcept>
#include <memory>
#include <algorithm>

std::string get_md5sum(const std::string& filename)
{
  // Command to execute: md5sum <filename>
  std::string command = "md5sum " + filename;
  std::array<char, 128> buffer;
  std::string result = "";
    
  // Use popen to execute the command and read its output
  // "r" mode opens the pipe for reading
  std::unique_ptr<FILE, decltype(&pclose)>
    pipe(popen(command.c_str(), "r"), pclose);
  if (!pipe) {
    throw std::runtime_error("popen() failed!");
  }
    
  // Read the output line by line
  while (fgets(buffer.data(), buffer.size(), pipe.get()) != nullptr) {
    result += buffer.data();
  }
    
  // The stanard GNU/Linux md5sum output format is: "32-char-hash
  // filename".  We extract the 32-character hash part...
  if (result.length() >= 32) {
    return result.substr(0, 32);
  } else {
    throw std::runtime_error("Failed to parse md5sum output.");
  }
}
