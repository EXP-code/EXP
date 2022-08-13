#include <iostream>
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
			HighFive::File::Create |
			HighFive::File::Truncate);

    // we create a new group
    HighFive::Group group = file.createGroup("Coefs");

    double time = 3.14159;
    
    HighFive::Attribute a = group.createAttribute<double>
      ("Time", HighFive::DataSpace::From(time));

    a.write(time);

    for (int j=0; j<5; j++) {

      Eigen::VectorXcd vector(nmax);

      for (int i=0; i<nmax; i++) {
	vector(i) = std::complex<double>(j + i * 100);
      }

      // Create the data set
      //
      std::ostringstream sout; sout << "DS" << j;
      HighFive::DataSet dataset = group.createDataSet(sout.str(), vector);
      
      std::vector<int> key = {1, 0, j};
      
      HighFive::Attribute v = dataset.createAttribute<int>
	("Index", HighFive::DataSpace::From(key));
      v.write(key);
    }

  } catch (HighFive::Exception& err) {
    std::cerr << err.what() << std::endl;
  }
  
  return 0;
}

