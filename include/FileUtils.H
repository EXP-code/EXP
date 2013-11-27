#ifndef _FileExists_H
#define _FileExists_H

#include <string>

extern bool FileExists(const std::string& filename);
extern bool FileRename(const std::string& file1, const std::string& file2);

#endif
