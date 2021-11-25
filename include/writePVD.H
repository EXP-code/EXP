#ifndef _WRITEPVD_H
#define _WRITEPVD_H

#include <iostream>
#include <iomanip>
#include <fstream>
#include <vector>
#include <string>


void writePVD(const std::string& filename,
	      const std::vector<double>& times,
	      const std::vector<std::string>& files);

#endif
