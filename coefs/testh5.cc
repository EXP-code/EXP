#include <iostream>
#include <iomanip>
#include <sstream>
#include <complex>

#include "config.h"


#include <Eigen/Dense>

#include <highfive/H5File.hpp>
#include <highfive/H5DataSet.hpp>
#include <highfive/H5DataSpace.hpp>
#include <highfive/H5Attribute.hpp>

const std::string FILE_NAME("test_coefs.h5");

const int nmax = 10;

// Test: create a 2D dataset 10x3 of complex value with eigen matrix
// and write it to a file
int main(void)
{
  try {
    // Create a new hdf5 file
    //
    HighFive::File file(FILE_NAME,
			HighFive::File::ReadWrite |
			HighFive::File::Create);

    double time = 3.1417;
    
    // we create a new group
    std::ostringstream stim; stim << std::setprecision(8) << time;
    HighFive::Group group = file.createGroup(stim.str());

    HighFive::Attribute t = group.createAttribute<double>("Time", HighFive::DataSpace::From(time));

    t.write(time);

    std::string config = "This is a string describing the configuration "
      "that generated the data";

    HighFive::Attribute c = group.createAttribute<std::string>("config", HighFive::DataSpace::From(config));

    c.write(config);

    for (int j=0; j<5; j++) {

      Eigen::VectorXcd vector(nmax);

      for (int i=0; i<nmax; i++) {
	vector(i) = std::complex<double>(j + i * 100, 3.14158*i);
      }

      // Create the data set
      //
      std::ostringstream sout; sout << j;
      HighFive::DataSet dataset = group.createDataSet(sout.str(), vector);
      
      std::vector<int> key = {1, j};
      
      HighFive::Attribute v = dataset.createAttribute<int>
	("Index", HighFive::DataSpace::From(key));
      
      v.write(key);
    }

  } catch (HighFive::Exception& err) {
    std::cerr << err.what() << std::endl;
  }
  
  return 0;
}

