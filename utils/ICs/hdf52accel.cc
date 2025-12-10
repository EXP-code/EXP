// Compile: g++ -I/usr/include/hdf5/serial -O3 -o hdf52accel hdf52accel.cc -lhdf5_serial -lhdf5_cpp

#include <iostream>
#include <iomanip>
#include <fstream>
#include <vector>
#include <string>

#include "config_exp.h"
#ifdef HAVE_HDF5
#include <H5Cpp.h>
#include <H5public.h>		// For version info
#endif

#include "cxxopts.H"

int 
main(int ac, char **av)
{
  //====================
  // Begin opt parsing
  //====================

  std::string       hdf5file;
  std::string       outfile;
  
  cxxopts::Options options(av[0], "Get the acceleration field from a Gadget-stype HDF5 file");

  options.add_options()
    ("help,h", "Print this help message")
    ("hdf5", "HDF5 Gadget2 file",
     cxxopts::value<std::string>(hdf5file)->default_value("snapfile_001.hdf5"))
    ("output", "Force data from N-body evluation",
     cxxopts::value<std::string>(outfile)->default_value("force.data"))
    ;
       
  cxxopts::ParseResult vm;

  try {
    vm = options.parse(ac, av);
  } catch (cxxopts::OptionException& e) {
    std::cout << "Option error: " << e.what() << std::endl;
    exit(-1);
  }

  // Print help message and exit
  //
  if (vm.count("help")) {
    std::cout << options.help() << std::endl << std::endl;
    return 1;
  }

#ifdef HAVE_HDF5

  // The H5 input file name
  //
  const H5std_string FILE_NAME(hdf5file);

  // Try block to detect H5 exceptions
  //
  try {

    // Turn off the auto-printing when failure occurs so that we can
    // handle the errors appropriately
    //
    H5::Exception::dontPrint();
    
    // Open the hdf5 file
    //
    auto file = std::make_shared<H5::H5File>(FILE_NAME, H5F_ACC_RDWR);
    
    // Open header
    //
    auto header   = std::make_shared<H5::Group>(file->openGroup("Header"));
    auto attr     = std::make_shared<H5::Attribute>(header->openAttribute("MassTable"));
    auto type    = std::make_shared<H5::DataType>(attr->getDataType());

    std::vector<double> mass(6);
    attr->read(*type, &mass[0]);
    std::cout << "Mass table: ";
    for (auto v : mass) std::cout << v << " ";
    std::cout << std::endl;

    // Get the disk particle group
    //
    auto group = std::make_shared<H5::Group>(file->openGroup("PartType2"));
    
    const H5std_string SetName1( "Coordinates" );
    const H5std_string SetName2( "Acceleration" );
    H5::DataSet dset1 = group->openDataSet(SetName1);
    H5::DataSet dset2 = group->openDataSet(SetName2);
    
    // Get the class of the datatype that is used by the dataset.
    //
    H5T_class_t type_class = dset1.getTypeClass();
    
    // Get class of datatype and print message if it's an integer.
    //
    if ( type_class == H5T_FLOAT ) {
      std::cout << "Data set has FLOAT type" << std::endl;
	
      // Get the float datatype
      //
      H5::FloatType floattype = dset1.getFloatType();
      
      // Get order of datatype and print message if it's a little endian.
      //
      H5std_string order_string;
      H5T_order_t order = floattype.getOrder( order_string );
      std::cout << order_string << std::endl;
      
      // Get size of the data element stored in file and print it.
      //
      size_t size = floattype.getSize();
      std::cout << "Data size is " << size << std::endl;
    }
    
    // Get dataspace of the datasets
    //
    H5::DataSpace dspace1 = dset1.getSpace();
    H5::DataSpace dspace2 = dset2.getSpace();
    
    // Get the number of dimensions in the dataspace.
    //
    int rank = dspace1.getSimpleExtentNdims();
    
    // Get the dimension size of each dimension in the dataspace and
    // display them.
    
    hsize_t dims_out[2];
    int ndims = dspace1.getSimpleExtentDims( dims_out, NULL);
    std::cout << "rank " << rank << ", dimensions " <<
      (unsigned long)(dims_out[0]) << " x " <<
      (unsigned long)(dims_out[1]) << std::endl;
    
    std::vector<float> coords(dims_out[0]*dims_out[1]);
    std::vector<float> accels(dims_out[0]*dims_out[1]);

    // Define the memory dataspace.
    //
    H5::DataSpace memspace1( 2, dims_out );
    H5::DataSpace memspace2( 2, dims_out );

    // Read data from dataset
    // 
    dset1.read( &coords[0], H5::PredType::NATIVE_FLOAT, memspace1, dspace1 );
    dset2.read( &accels[0], H5::PredType::NATIVE_FLOAT, memspace2, dspace2 );

    std::ofstream out(outfile);

    float f;
    out.write((const char *)&dims_out[0], sizeof(int));
    for (int j=0; j<dims_out[0]; j++) {
      out.write((const char *)&(f=mass[2]), sizeof(float)  );
      out.write((const char *)&coords[j*3], sizeof(float)*3);
      out.write((const char *)&accels[j*3], sizeof(float)*3);
    }
    
  }
  // end of try block
  //
  // catch failure caused by the H5File operations
  catch( H5::FileIException error )
    {
#if (H5_VERS_MAJOR>1) || ((H5_VERS_MAJOR==1) && ((H5_VERS_MINOR>8)) || ((H5_VERSH5_VERS_MINOR==8)&&(H5_VERS_RELEASE>=20)))
      error.printErrorStack();
#else
      error.printError();
#endif
      return -1;
    }
  // catch failure caused by the DataSet operations
  catch( H5::DataSetIException error )
    {
#if (H5_VERS_MAJOR>1) || ((H5_VERS_MAJOR==1)&& ((H5_VERS_MINOR>8)) || ((H5_VERSH5_VERS_MINOR==8)&&(H5_VERS_RELEASE>=20)))
      error.printErrorStack();
#else
      error.printError();
#endif
      return -1;
    }
  // catch failure caused by the DataSpace operations
  catch( H5::DataSpaceIException error )
    {
#if (H5_VERS_MAJOR>1) || ((H5_VERS_MAJOR==1)&& ((H5_VERS_MINOR>8)) || ((H5_VERSH5_VERS_MINOR==8)&&(H5_VERS_RELEASE>=20)))
      error.printErrorStack();
#else
      error.printError();
#endif
      return -1;
    }
  // catch failure caused by the Attribute operations
  catch( H5::AttributeIException error )
    {
#if (H5_VERS_MAJOR>1) || ((H5_VERS_MAJOR==1)&& ((H5_VERS_MINOR>8)) || ((H5_VERSH5_VERS_MINOR==8)&&(H5_VERS_RELEASE>=20)))
      error.printErrorStack();
#else
      error.printError();
#endif
      return -1;
    }

#else
  std::cout << "No HDF5 support in your environment . . . sorry" << std::endl;
#endif

  return 0;
}
