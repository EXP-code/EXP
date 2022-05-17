#ifndef _CenterFile_H
#define _CenterFile_H

#include <fstream>
#include <vector>
#include <array>

#include <yaml-cpp/yaml.h>

/**
   Read a file in @Orient format and provide the center interpolated
   in time.   The YAML node needs two key word--value pairs:

   file: file name of orient-style center info
   type: EJ or COM

 */
class CenterFile
{
private:
  std::vector<double> time;
  std::vector<std::array<double, 3>> data;

public:

  //! Read in the file and cache the data
  CenterFile(const YAML::Node& conf);

  //! Interpolate and provide the desired center
  std::array<double, 3> operator()(double T);

};

#endif
