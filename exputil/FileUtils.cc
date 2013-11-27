#include <cstdio>
#include <string>
#include <sys/stat.h>

//
// *** Check for file existence ***
// We could make this throwable
//
bool FileExists(const std::string& filename) 
{
  bool ret;
  int status;
  struct stat info;

  // Get the file attributes
  status = stat(filename.c_str(), &info);
  if(status == 0) {
    // Successful! The file MUST exist . . .
    ret = true;
  } else {
    // We could check for, and report, other issues here . . .
    ret = false;
  }
  
  return ret;
}

//
// We could make this throwable
//
bool FileRename(const std::string& ifile, const std::string& ofile) 
{
  bool ret = false;
  int status;
  struct stat info;

  // Get the file attributes
  status = stat(ifile.c_str(), &info);
  if(status == 0) {
    // Successful! The file MUST exist . . .
    status = rename(ifile.c_str(), ofile.c_str());
    if (status == 0) ret = true;
  }
  
  return ret;
}
