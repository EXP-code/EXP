#include <algorithm>

#include "Covariance.H"

namespace BasisClasses
{
  
  // Main constructor
  SubsampleCovariance::SubsampleCovariance
  (ParamCallback func, const std::string& BasisID, bool flt, bool sum, bool cov)
    : paramCallback(func), BasisID(BasisID), floatType(flt),
      summed(sum), covar(cov)
  {
    if (summed) covar = true;
  }

  unsigned SubsampleCovariance::writeCovarH5
  (HighFive::Group& snaps, CovarData& elem, unsigned count, double time)
  {
    std::ostringstream stim;
    stim << std::setw(8) << std::setfill('0') << std::right << count++;
    
    // Make a new group for this time
    //
    HighFive::Group stanza = snaps.createGroup(stim.str());
    
    // Add a time attribute
    //
    time = roundTime(time);
    stanza.createAttribute<double>("Time", HighFive::DataSpace::From(time)).write(time);
    
    // Enable compression
    //
    auto dcpl1 = HighFive::DataSetCreateProps{}; // sample stats
    auto dcpl2 = HighFive::DataSetCreateProps{}; // coefficients
    auto dcpl3 = HighFive::DataSetCreateProps{}; // covariance
    auto dcpl4 = HighFive::DataSetCreateProps{}; // total covariance
    
    // Get the dimensions as a DSizes object
    const auto& dims = std::get<2>(elem).dimensions();

    // Number of samples
    //
    unsigned sampleSize   = dims[0];
    unsigned ltot         = dims[1];
    unsigned nmax         = dims[2];
    unsigned diagonalSize = nmax;
    
    // Properties for sample stats
    if (H5compress) {
      unsigned int csz = sampleSize;
      dcpl1.add(HighFive::Chunking({csz, 1}));
      if (H5shuffle) dcpl1.add(HighFive::Shuffle());
      dcpl1.add(HighFive::Deflate(H5compress));
    }
    
    // Add the sample statistics
    //
    HighFive::DataSet s1data = stanza.createDataSet("sampleCounts", std::get<0>(elem), dcpl1);
    HighFive::DataSet s2data = stanza.createDataSet("sampleMasses", std::get<1>(elem), dcpl1);
    
    // Save variance or full covariance
    if (covar) diagonalSize = nmax*(nmax + 1)/2;
    
    // Add data dimensions
    //
    stanza.createAttribute<unsigned>
      ("sampleSize", HighFive::DataSpace::From(sampleSize)).write(sampleSize);
    
    stanza.createAttribute<unsigned>
      ("angularSize", HighFive::DataSpace::From(ltot)).write(ltot);
    
    stanza.createAttribute<unsigned>
      ("rankSize", HighFive::DataSpace::From(nmax)).write(nmax);
    
    if (H5compress) {
      // Szip parameters
      const int options_mask = H5_SZIP_NN_OPTION_MASK;
      const int pixels_per_block = 8;
      
      // Properties for coefficients
      //
      unsigned int csz2 = nmax * ltot * sampleSize;
      HighFive::Chunking data_dims2{std::min<unsigned>(csz2, H5chunk), 1};
      
      dcpl2.add(data_dims2);
      if (H5shuffle) dcpl2.add(HighFive::Shuffle());
      if (H5szip) {
	dcpl2.add(HighFive::Szip(options_mask, pixels_per_block));
      } else {
	dcpl2.add(HighFive::Deflate(H5compress));
      }
      
      // Properties for covariance
      //
      unsigned int csz3 = ltot * diagonalSize * sampleSize;
      HighFive::Chunking data_dims3{std::min<unsigned>(csz3, H5chunk), 1};
      
      dcpl3.add(data_dims3);
      if (H5shuffle) dcpl3.add(HighFive::Shuffle());
      if (H5szip) {
	dcpl3.add(HighFive::Szip(options_mask, pixels_per_block));
      } else {
	dcpl3.add(HighFive::Deflate(H5compress));
      }

      // Properties for total covariance
      //
      if (summed) {
	unsigned int csz4 = ltot * diagonalSize;
	HighFive::Chunking data_dims4{std::min<unsigned>(csz4, H5chunk), 1};
      
	dcpl4.add(data_dims4);
	if (H5shuffle) dcpl4.add(HighFive::Shuffle());
	if (H5szip) {
	  dcpl4.add(HighFive::Szip(options_mask, pixels_per_block));
	} else {
	  dcpl4.add(HighFive::Deflate(H5compress));
	}
      }
    }
    
    // Pack the coefficient data
    //
    if (floatType) {
      // Create a vector of doubles for the real and imaginary parts
      Eigen::VectorXf real_part(nmax*ltot*sampleSize);
      Eigen::VectorXf imag_part(nmax*ltot*sampleSize);
      
      for (int T=0, c=0; T<sampleSize; T++) {
	for (int l=0; l<ltot; l++) {
	  for (int n=0; n<nmax; n++, c++) {
	    real_part(c) = std::real(std::get<2>(elem)(T, l, n));
	    imag_part(c) = std::imag(std::get<2>(elem)(T, l, n));
	  }
	}
      }
      
      // Create two separate, compressed datasets
      stanza.createDataSet("coefficients_real", real_part, dcpl2);
      stanza.createDataSet("coefficients_imag", imag_part, dcpl2);
      
      // Get the dimensions of the covariance data
      const auto& dimsC = std::get<3>(elem).dimensions();

      if (dimsC[0] and dimsC[1] and dimsC[2]) {
	
	// Store summed covariance only
	if (summed) {

	  real_part.resize(ltot*diagonalSize);
	  imag_part.resize(ltot*diagonalSize);

	  real_part.setZero();
	  imag_part.setZero();
	
	  for (int T=0, c=0; T<sampleSize; T++) {
	    for (int l=0; l<ltot; l++) {
	      for (int n1=0; n1<nmax; n1++) {
		// Pack the covariance data in an upper triangular format
		//
		for (int n2=n1; n2<nmax; n2++, c++) {
		  real_part(c) += std::real(std::get<3>(elem)(T, l, n1, n2));
		  imag_part(c) += std::imag(std::get<3>(elem)(T, l, n1, n2));
		}
	      }
	    }
	  }
	  // Create two separate, compressed datasets
	  stanza.createDataSet("covariance_real_total", real_part, dcpl4);
	  stanza.createDataSet("covariance_imag_total", imag_part, dcpl4);

	} else {

	  real_part.resize(ltot*diagonalSize*sampleSize);
	  imag_part.resize(ltot*diagonalSize*sampleSize);
	
	  for (int T=0, c=0; T<sampleSize; T++) {
	    for (int l=0; l<ltot; l++) {
	      for (int n1=0; n1<nmax; n1++) {
		// Pack the covariance data in an upper triangular format
		//
		if (covar) {
		  for (int n2=n1; n2<nmax; n2++, c++) {
		    real_part(c) = std::real(std::get<3>(elem)(T, l, n1, n2));
		    imag_part(c) = std::imag(std::get<3>(elem)(T, l, n1, n2));
		  }
		}
		// Pack the diagonal only
		//
		else {
		  real_part(c  ) = std::real(std::get<3>(elem)(T, l, n1, n1));
		  imag_part(c++) = std::imag(std::get<3>(elem)(T, l, n1, n1));
		}
	      }
	    }
	  }
	  // Create two separate, compressed datasets
	  stanza.createDataSet("covariance_real", real_part, dcpl3);
	  stanza.createDataSet("covariance_imag", imag_part, dcpl3);
	}
      }
      
    } else {
      Eigen::VectorXd real_part(ltot*nmax*sampleSize);
      Eigen::VectorXd imag_part(ltot*nmax*sampleSize);
      
      for (int T=0, c=0; T<sampleSize; T++) {
	for (int l=0; l<ltot; l++) {
	  for (int n=0; n<nmax; n++, c++) {
	    real_part(c) = std::real(std::get<2>(elem)(T, l, n));
	    imag_part(c) = std::imag(std::get<2>(elem)(T, l, n));
	  }
	}
      }
      
      // Create two separate, compressed datasets
      //
      stanza.createDataSet("coefficients_real", real_part, dcpl2);
      stanza.createDataSet("coefficients_imag", imag_part, dcpl2);
      
      // Get the dimensions of the covariance data
      const auto& dimsC = std::get<3>(elem).dimensions();

      if (dimsC[0] and dimsC[1] and dimsC[2]) {
	
	if (summed) {
	  
	  real_part.resize(ltot*diagonalSize);
	  imag_part.resize(ltot*diagonalSize);
	
	  real_part.setZero();
	  imag_part.setZero();
	
	  for (int T=0; T<sampleSize; T++) {
	    for (int l=0, c=0; l<ltot; l++) {
	      for (int n1=0; n1<nmax; n1++) {
		// Pack the covariance data in an upper triangular format
		//
		for (int n2=n1; n2<nmax; n2++, c++) {
		  real_part(c) += std::real(std::get<3>(elem)(T, l, n1, n2));
		  imag_part(c) += std::imag(std::get<3>(elem)(T, l, n1, n2));
		}
	      }
	    }
	  }
	
	  // Create two separate, compressed datasets
	  //
	  stanza.createDataSet("covariance_real_total", real_part, dcpl4);
	  stanza.createDataSet("covariance_imag_total", imag_part, dcpl4);
	}
	else {

	  real_part.resize(ltot*diagonalSize*sampleSize);
	  imag_part.resize(ltot*diagonalSize*sampleSize);
	
	  for (int T=0, c=0; T<sampleSize; T++) {
	    for (int l=0; l<ltot; l++) {
	      for (int n1=0; n1<nmax; n1++) {
		// Pack the covariance data in an upper triangular format
		//
		if (covar) {
		  for (int n2=n1; n2<nmax; n2++, c++) {
		    real_part(c) = std::real(std::get<3>(elem)(T, l, n1, n2));
		    imag_part(c) = std::imag(std::get<3>(elem)(T, l, n1, n2));
		  }
		}
		// Pack the diagonal only
		//
		else {
		  real_part(c  ) = std::real(std::get<3>(elem)(T, l, n1, n1));
		  imag_part(c++) = std::imag(std::get<3>(elem)(T, l, n1, n1));
		}
	      }
	    }
	  }
	
	  // Create two separate, compressed datasets
	  //
	  stanza.createDataSet("covariance_real", real_part, dcpl3);
	  stanza.createDataSet("covariance_imag", imag_part, dcpl3);
	}
      }
    }
    // END: sample loop
    
    return count;
  }
  
  void SubsampleCovariance::writeCoefCovariance
  (const std::string& compname, const std::string& runtag, CovarData& elem, double time)
  {
    // Only root process writes
    //
    if (myid) return;

    // The H5 filename
    //
    std::string fname = "coefcovar." + compname + "." + runtag + ".h5";
    
    writeCoefCovariance(fname, elem, time);
  }
  
  void SubsampleCovariance::writeCoefCovariance
  (const std::string& fname, CovarData& elem, double time)
  {
    // Only root process writes
    //
    if (myid) return;
    
    // Check that there is something to write
    //
    int totalCount = 0;
    totalCount += std::get<0>(elem).sum();
    
    if (totalCount==0) {
      std::cout << "SubsampleCovariance::writeCoefCovariance: no data" << std::endl;
      return;
    }
    
    // Round time
    //
    time = roundTime(time);
    
    // Check if file exists?
    //
    try {
      // Open the HDF5 file in read-write mode, creating if it doesn't
      // exist
      HighFive::File file(fname,
			  HighFive::File::ReadWrite |
			  HighFive::File::Create);
      
      // Check for version string
      std::string path = "CovarianceFileVersion"; 
      
      // Check for valid HDF file by attribute
      if (file.hasAttribute(path)) {
	extendCoefCovariance(fname, elem, time);
	return;
      }
      
      // Write the Version string
      //
      file.createAttribute<std::string>("CovarianceFileVersion", HighFive::DataSpace::From(CovarianceFileVersion)).write(CovarianceFileVersion);
      
      // Write the basis identifier string
      //
      file.createAttribute<std::string>("BasisID", HighFive::DataSpace::From(BasisID)).write(BasisID);
      
      // Write the data type size
      //
      int sz = 8; if (floatType) sz = 4;
      file.createAttribute<int>("FloatSize", HighFive::DataSpace::From(sz)).write(sz);
      
      // Write the specific parameters
      //
      writeCovarH5Params(file);
      
      // Group count variable
      //
      unsigned count = 0;
      HighFive::DataSet dataset = file.createDataSet("count", count);
      
      // Create a new group for coefficient snapshots
      //
      HighFive::Group group = file.createGroup("snapshots");
      
      // Write the coefficients
      //
      count = writeCovarH5(group, elem, count, time);
      
      // Update the count
      //
      dataset.write(count);
      
    } catch (const HighFive::Exception& err) {
      // Handle HighFive specific errors (e.g., file not found)
      throw std::runtime_error
	(std::string("SubsampleCovariance::writeCoefCovariance HighFive Error: ") + err.what());
    } catch (const std::exception& err) {
      // Handle other general exceptions
      throw std::runtime_error
	(std::string("SubsampleCovariance::writeCoefCovariance Error: ") + err.what());
    }
  }
  
  void SubsampleCovariance::extendCoefCovariance
  (const std::string& fname, CovarData& elem, double time)
  {
    try {
      // Open an hdf5 file
      //
      HighFive::File file(fname, HighFive::File::ReadWrite);
      
      // Get the dataset
      HighFive::DataSet dataset = file.getDataSet("count");
      
      unsigned count;
      dataset.read(count);
      
      HighFive::Group group = file.getGroup("snapshots");
      
      // Write the coefficients
      //
      count = writeCovarH5(group, elem, count, time);
      
      // Update the count
      //
      dataset.write(count);
      
    } catch (HighFive::Exception& err) {
      throw std::runtime_error
	(std::string("SubsampleCovariance::extendCoefCovariance: HighFive error: ") + err.what());
    }
  }
  
  // Read covariance data
  SubsampleCovariance::SubsampleCovariance
  (const std::string& filename, int stride)
  {
    try {
      // Open an existing hdf5 file for reading
      //
      HighFive::File file(filename, HighFive::File::ReadOnly);
      
      // Write the Version string
      //
      std::string version;
      file.getAttribute("CovarianceFileVersion").read(version);
      // Check for alpha version
      if (version == std::string("1.0")) {
	throw std::runtime_error("SubsampleCovariance: this is an early alpha test version. Please remake your files");
      }
      // Test for current version
      if (version != std::string("1.1")) {
	throw std::runtime_error(std::string("SubsampleCovariance: unsupported file version, ") + version);
      }
      
      // Read the basis identifier string
      //
      file.getAttribute("BasisID").read(BasisID);
      
      // Get the float size
      int sz = 8;
      file.getAttribute("FloatSize").read(sz);
      if (sz != 4 and sz != 8) {
	std::ostringstream sout;
	sout << "SubsampleCovariance: unsupported float size, " << sz;
	throw std::runtime_error(sout.str());
      }
      
      int lmax, nmax, ltot;
      
      // Current implemented spherical types
      const std::set<std::string> sphereType = {"Spherical", "SphereSL", "Bessel"};
      
      // Currently implemented cylindrical types
      const std::set<std::string> cylinderType = {"Cylindrical", "Cylinder"};
      
      std::cout << "SubsampleCovariance: reading basis type " << BasisID << std::endl;

      if (sphereType.find(BasisID) != sphereType.end()) {
	file.getAttribute("lmax").read(lmax);
	file.getAttribute("nmax").read(nmax);
	ltot = (lmax+1)*(lmax+2)/2;
      } else if (cylinderType.find(BasisID) != cylinderType.end()) {
	file.getAttribute("mmax").read(lmax);
	file.getAttribute("nmax").read(nmax);
	ltot = lmax + 1;
      } else if (BasisID == "Cube") {
	int nmaxx, nmaxy, nmaxz;
	file.getAttribute("nmaxx").read(nmaxx);
	file.getAttribute("nmaxy").read(nmaxy);
	file.getAttribute("nmaxz").read(nmaxz);
	ltot = (2*nmaxx + 1) * (2*nmaxy + 1) * (2*nmaxz + 1);
      } else {
	throw std::runtime_error(std::string("SubsampleCovariance: unknown or unimplemented covariance for basis type, ") + BasisID);
      }
      
      // Group count variable
      //
      unsigned count = 0;
      file.getDataSet("count").read(count);
      
      // Open the snapshot group
      //
      auto snaps = file.getGroup("snapshots");
      
      for (unsigned n=0; n<count; n+=stride) {
	
	std::ostringstream sout;
	sout << std::setw(8) << std::setfill('0') << std::right << n;
	
	auto stanza = snaps.getGroup(sout.str());
	
	double Time;
	stanza.getAttribute("Time").read(Time);
	Time = roundTime(Time);
	times.push_back(Time);
	
	// Get sample properties
	//
	Eigen::VectorXi counts;
	Eigen::VectorXd masses;
	
	stanza.getDataSet("sampleCounts").read(counts);
	stanza.getDataSet("sampleMasses").read(masses);
	
	// Get data attributes
	//
	int nT, lSize, rank, icov=1;
	stanza.getAttribute("sampleSize") .read(nT);
	stanza.getAttribute("angularSize").read(lSize);
	stanza.getAttribute("rankSize")   .read(rank);
	
	// Allocate sample vector for current time
	//
	CovarData & elem = covarData[Time];

	std::get<0>(elem) = counts;
	std::get<1>(elem) = masses;
	
	// Temporary storage for covariance elements

	// Storage
	//
	Eigen::VectorXcd data0, data1;
	
	// Get the flattened coefficient array
	//
	if (sz==4) {
	  // Get the real and imaginary parts
	  //
	  Eigen::VectorXf data_real =
	    stanza.getDataSet("coefficients_real").read<Eigen::VectorXf>();
	  
	  Eigen::VectorXf data_imag =
	    stanza.getDataSet("coefficients_imag").read<Eigen::VectorXf>();
	  
	  // Resize the complex array and assign
	  //
	  data0.resize(data_real.size());
	  data0.real() = data_real.cast<double>();
	  data0.imag() = data_imag.cast<double>();
	  
	  summed = false;

	  // Check for existence of total covariance
	  //
	  if (stanza.exist("covariance_real_total")) {

	    summed = true;
	    
	    data_real =
	      stanza.getDataSet("covariance_real_total").read<Eigen::VectorXf>();
	    
	    data_imag =
	      stanza.getDataSet("covariance_imag_total").read<Eigen::VectorXf>();
	    
	    // Resize the complex array and assign
	    data1.resize(data_real.size());
	    data1.real() = data_real.cast<double>();
	    data1.imag() = data_imag.cast<double>();
	  }
	  // Check for existence of subsample covariance
	  //
	  else if (stanza.exist("covariance_real")) {
	    
	    data_real =
	      stanza.getDataSet("covariance_real").read<Eigen::VectorXf>();
	    
	    data_imag =
	      stanza.getDataSet("covariance_imag").read<Eigen::VectorXf>();
	    
	    // Resize the complex array and assign
	    data1.resize(data_real.size());
	    data1.real() = data_real.cast<double>();
	    data1.imag() = data_imag.cast<double>();
	  }
	} else {
	  // Get the real and imaginary parts
	  Eigen::VectorXd data_real =
	    stanza.getDataSet("coefficients_real").read<Eigen::VectorXd>();
	  
	  Eigen::VectorXd data_imag =
	    stanza.getDataSet("coefficients_imag").read<Eigen::VectorXd>();
	  
	  // Resize the complex array and assign
	  data0.resize(data_real.size());
	  data0.real() = data_real;
	  data0.imag() = data_imag;
	  
	  summed = false;

	  // Check for existence of total covariance
	  //
	  if (stanza.exist("covariance_real_total")) {

	    summed = true;
	    
	    // Get the real and imaginary parts
	    data_real =
	      stanza.getDataSet("covariance_real_total").read<Eigen::VectorXd>();
	    
	    data_imag =
	      stanza.getDataSet("covariance_imag_total").read<Eigen::VectorXd>();
	    
	    // Resize the complex array and assign
	    data1.resize(data_real.size());
	    data1.real() = data_real;
	    data1.imag() = data_imag;
	  }
	  // Check for existence of subsample covariance
	  //
	  else if (stanza.exist("covariance_real")) {
	    
	    // Get the real and imaginary parts
	    data_real =
	      stanza.getDataSet("covariance_real").read<Eigen::VectorXd>();
	    
	    data_imag =
	      stanza.getDataSet("covariance_imag").read<Eigen::VectorXd>();
	    
	    // Resize the complex array and assign
	    data1.resize(data_real.size());
	    data1.real() = data_real;
	    data1.imag() = data_imag;
	  }
	}
	
	// Positions in data stanzas
	int sCof = 0, sCov = 0;
	
	std::get<2>(elem) = SubsampleCovariance::CoefType(nT, lSize, rank);
	if (data1.size())
	  std::get<3>(elem) =
	    SubsampleCovariance::CovrType(nT, lSize, rank, rank);

	// Loop through all indices and repack
	//
	for (int T=0; T<nT; T++) {
	  
	  // Pack the coefficient data
	  int c = 0;
	  for (int l=0; l<lSize; l++) {
	    for (int n=0; n<rank; n++) {
	      std::get<2>(elem)(T, l, n) = data0(sCof + c++);
	    }
	  }
	  sCof += c;
	  
	  // Pack the covariance data, if we have it
	  //
	  if (data1.size()) {
	    c = 0;		// Position counter
	    for (int l=0; l<lSize; l++) {
	      for (int n1=0; n1<rank; n1++) {
		for (int n2=n1; n2<rank; n2++) {
		  if (covar) {
		    std::get<3>(elem)(T, l, n1, n2) = data1(sCov + c++);
		    if (n1 != n2)
		      std::get<3>(elem)(T, l, n2, n1) =
			std::get<3>(elem)(T, l, n1, n2);
		  }
		  else {
		    if (n1==n2)
		      std::get<3>(elem)(T, l, n1, n2) = data1(sCov + c++);
		    else
		      std::get<3>(elem)(T, l, n1, n2) = 0.0;
		  }
		}
	      }
	    }

	    // Advance the covariance position for subsampled data only.
	    // For summed covariance, we will reuse the same data.
	    if (not summed) sCov += c;
	  }
	  // END: pack covariance data
	  
	}
	// END: sample loop
	
	// Split the summed covariance nT ways to emulate
	// subsample covariance
	if (summed) {
	  CovarData & elem = covarData[Time];
	  std::get<3>(elem) *= std::get<3>(elem).constant(1.0/nT);
	}
      }
      // END: snapshot loop

    } catch (HighFive::Exception& err) {
      std::cerr << err.what() << std::endl;
    }
  }
  
}
// END: namespace BasisClasses
