#include <iostream>
#include <cstring>
#include <cassert>
#include <cmath>

#include <yaml-cpp/yaml.h>	// YAML support


#include <coef.H>		// Old original header
#include <CoefStruct.H>		// New coefficient structure

namespace CoefClasses
{
  
  void CoefStruct::copyfields(std::shared_ptr<CoefStruct> ret)
  {
    ret->buf   = buf;
    ret->geom  = geom;
    ret->id    = id;
    ret->time  = time;
    ret->store = store;
    if (ret->ctr.size()) ret->ctr = ctr;
  }

  void CylStruct::create()
  {
    if (nmax>0)
      coefs = std::make_shared<coefType>(store.data(), mmax+1, nmax);
    else
      throw std::runtime_error("CylStruct::create: nmax must be >0");
  }

  void SphStruct::create()
  {
    if (nmax>0)
      coefs = std::make_shared<coefType>(store.data(), (lmax+1)*(lmax+2)/2, nmax);
    else
      throw std::runtime_error("SphStruct::create: nmax must be >0");
  }

  void SlabStruct::create()
  {
    if (nmaxx>=0 and nmaxy>=0 and nmaxz>0) {
      nx = 2*nmaxx + 1;
      ny = 2*nmaxy + 1;
      nz = nmaxz;
      dim = nx * ny * nz;
      coefs = std::make_shared<coefType>(store.data(), nx, ny, nz);
    } else
      throw std::runtime_error("SlabStruct::create: all dimensions must be >=0");
  }

  void CubeStruct::create()
  {
    if (nmaxx>0 and nmaxy>0 and nmaxz>0) {
      nx = 2*nmaxx + 1;
      ny = 2*nmaxy + 1;
      nz = 2*nmaxz + 1;
      dim = nx * ny * nz;
      coefs = std::make_shared<coefType>(store.data(), nx, ny, nz);
    } else
      throw std::runtime_error("CubeStruct::create: all dimensions must be >0");
  }

  void TblStruct::create()
  {
    if (cols>0)
      coefs = std::make_shared<coefType>(store.data(), cols);
    else
      throw std::runtime_error("TblStruct::create: cols must be >0");
  }

  void TrajStruct::create()
  {
    if (traj>0 and rank>0)
      coefs = std::make_shared<coefType>(store.data(), traj, rank);
    else
      throw std::runtime_error("TrajStruct::create: number of trajectories and phase-space size must be >0");
  }

  void SphFldStruct::create()
  {
    if (nfld>0 and nmax>0)
      coefs = std::make_shared<coefType>(store.data(), nfld, (lmax+1)*(lmax+2)/2, nmax);
    else
      throw std::runtime_error("CoefStruct::SphFldStruct: nfld and nmax must be >0");
  }

  void CylFldStruct::create()
  {
    if (nfld>0 and nmax>0)
      coefs = std::make_shared<coefType>(store.data(), nfld, mmax+1, nmax);
    else
      throw std::runtime_error("CoefStruct::CylFldStruct: nfld and nmax must be >0");
  }

  std::shared_ptr<CoefStruct> CylStruct::deepcopy()
  {
    auto ret = std::make_shared<CylStruct>();

    copyfields(ret);

    assert(("CylStruct::deepcopy dimension mismatch",
	    (mmax+1)*nmax == store.size()));

    ret->coefs = std::make_shared<coefType>(ret->store.data(), mmax+1, nmax);
    ret->mmax  = mmax;
    ret->nmax  = nmax;

    return ret;
  }

  std::shared_ptr<CoefStruct> SphStruct::deepcopy()
  {
    auto ret = std::make_shared<SphStruct>();

    copyfields(ret);

    assert(("SphStruct::deepcopy dimension mismatch",
	    (lmax+1)*(lmax+2)/2*nmax == store.size()));

    ret->coefs  = std::make_shared<coefType>
      (ret->store.data(), (lmax+1)*(lmax+2)/2, nmax);
    ret->lmax   = lmax;
    ret->nmax   = nmax;
    ret->scale  = scale;
    ret->normed = normed;

    return ret;
  }

  std::shared_ptr<CoefStruct> SlabStruct::deepcopy()
  {
    auto ret = std::make_shared<SlabStruct>();

    copyfields(ret);

    assert(("SlabStruct::deepcopy dimension mismatch", dim == store.size()));

    ret->coefs  = std::make_shared<coefType>(ret->store.data(), nx, ny, nz);
    ret->nmaxx  = nmaxx;
    ret->nmaxy  = nmaxy;
    ret->nmaxz  = nmaxz;
    ret->nx     = nx;
    ret->ny     = ny;
    ret->nz     = nz;
    ret->dim    = dim;

    return ret;
  }

  std::shared_ptr<CoefStruct> CubeStruct::deepcopy()
  {
    auto ret = std::make_shared<CubeStruct>();

    copyfields(ret);

    assert(("CubeStruct::deepcopy dimension mismatch", dim == store.size()));

    ret->coefs  = std::make_shared<coefType>(ret->store.data(), nx, ny, nz);
    ret->nmaxx  = nmaxx;
    ret->nmaxy  = nmaxy;
    ret->nmaxz  = nmaxz;
    ret->nx     = nx;
    ret->ny     = ny;
    ret->nz     = nz;
    ret->dim    = dim;

    return ret;
  }

  std::shared_ptr<CoefStruct> TblStruct::deepcopy()
  {
    auto ret = std::make_shared<TblStruct>();

    copyfields(ret);

    assert(("TblStruct::deepcopy dimension mismatch",
	    cols == store.size()));


    ret->coefs = std::make_shared<coefType>(ret->store.data(), cols);
    ret->cols  = cols;

    return ret;
  }

  std::shared_ptr<CoefStruct> TrajStruct::deepcopy()
  {
    auto ret = std::make_shared<TrajStruct>();

    copyfields(ret);

    assert(("TrajStruct::deepcopy dimension mismatch",
	    traj*rank == store.size()));


    ret->coefs = std::make_shared<coefType>(ret->store.data(), traj, rank);
    ret->traj  = traj;
    ret->rank  = rank;

    return ret;
  }

  std::shared_ptr<CoefStruct> SphFldStruct::deepcopy()
  {
    auto ret = std::make_shared<SphFldStruct>();

    copyfields(ret);

    assert(("SphFldStruct::deepcopy dimension mismatch",
	    nfld*(lmax+1)*(lmax+2)/2*nmax == store.size()));

    ret->coefs  = std::make_shared<coefType>
      (ret->store.data(), nfld, (lmax+1)*(lmax+2)/2, nmax);
    ret->nfld   = nfld;
    ret->lmax   = lmax;
    ret->nmax   = nmax;
    ret->scale  = scale;

    return ret;
  }

  std::shared_ptr<CoefStruct> CylFldStruct::deepcopy()
  {
    auto ret = std::make_shared<CylFldStruct>();

    copyfields(ret);

    assert(("CylFldStruct::deepcopy dimension mismatch",
	    nfld*(mmax+1)*nmax == store.size()));

    ret->coefs  = std::make_shared<coefType>
      (ret->store.data(), nfld, mmax+1, nmax);
    ret->nfld   = nfld;
    ret->mmax   = mmax;
    ret->nmax   = nmax;
    ret->scale  = scale;

    return ret;
  }


  bool CylStruct::read(std::istream& in, bool exp_type, bool verbose)
  {
    // iostream exception handling
    //
    in.exceptions ( std::istream::failbit | std::istream::badbit );
    
    // Save initial stream position
    //
    auto curpos = in.tellg();
    
    // Attempt to read coefficient magic number
    //
    const unsigned int cmagic = 0xc0a57a3;
    unsigned int tmagic;
    
    // Catch iostream exceptions on reading
    //
    try {
      
      in.read(reinterpret_cast<char*>(&tmagic), sizeof(unsigned int));
      
      
      if (tmagic == cmagic) {
	
	// YAML size
	//
	unsigned ssize;
	in.read(reinterpret_cast<char*>(&ssize), sizeof(unsigned int));
	
	// Make char buffer for YAML config
	//
	auto cbuf = std::make_unique<char[]>(ssize+1);

	// Read YAML string
	//
	in.read(cbuf.get(), ssize);
	cbuf[ssize] = '\0';
	buf = std::string(cbuf.get());
	
	YAML::Node node = YAML::Load(buf.c_str());
	
	// Get parameters
	//
	time = node["time"].as<double>();
	nmax = node["nmax"].as<int>();
	mmax = node["mmax"].as<int>();

	if (node["geom"]) geom = node["geom"].as<std::string>();
	else              geom = "cylinder";
	if (node["id"  ]) id   = node["id"  ].as<std::string>();
	
	if (verbose)
	  std::cerr << "New header: T=" << time << " nmax=" << nmax
		    << " mmax=" << mmax << std::endl;
      } else {
	
	// Rewind file
	//
	in.clear();
	in.seekg(curpos);
	
	CylCoefHeader header;
	in.read((char *)&header, sizeof(CylCoefHeader));
	
	time = header.time;
	nmax = header.nmax;
	mmax = header.mmax;
	geom = "cylinder";
	id   = "Cylinder";
	
	if (verbose)
	  std::cerr << "Old header: T=" << time << " nmax=" << nmax
		    << " mmax=" << mmax << std::endl;
      }
    }
    catch (std::istream::failure e) {
      if (not in.eof())
	std::cerr << "Exception reading coefficient file at T="
		  << time << ": " << e.what() << std::endl;
      return false;
    }
    
    store.resize((mmax+1)*nmax);
    coefs = std::make_shared<coefType>(store.data(), mmax+1, nmax);
    Eigen::VectorXd workR(nmax), workC(nmax);
    workC.setZero();
    
    try {
      for (int mm=0; mm<=mmax; mm++) {
	
	in.read((char *)workR.data(), sizeof(double)*nmax);
	if (mm) in.read((char *)workC.data(), sizeof(double)*nmax);
	
	coefs->row(mm).real() = workR;
	coefs->row(mm).imag() = workC;
      }
      
      if (verbose) {
	if (in)
	  std::cerr << "Coefficients successfully read at T=" << time << std::endl;
	else
	  std::cerr << "Coefficient read FAILED at T=" << time << std::endl;
      }
    }
    catch (std::istream::failure e) {
      if (not in.eof())
	std::cerr << "Error reading data at T="
		  << time << ": " << e.what() << std::endl;
      return false;
    }
    
    return true;
  }
  
  bool SphStruct::read(std::istream& in, bool exp_type, bool verbose)
  {
    in.exceptions ( std::istream::failbit | std::istream::badbit );
    
    // Save initial stream position
    //
    auto curpos = in.tellg();
    
    // Norm flag
    //
    bool normed = false;
    
    // Coefficient magic number
    //
    const unsigned int cmagic = 0xc0a57a2;
    
    // Try to read magic #
    //
    unsigned int tmagic;
    
    // Catch iostream exceptions on reading
    //
    try {
      in.read(reinterpret_cast<char *>(&tmagic), sizeof(unsigned int));
      
      // Found new-style coefficient file
      //
      if (cmagic == tmagic) {
	
	// Read YAML string size
	//
	unsigned int hsize;
	in.read(reinterpret_cast<char *>(&hsize), sizeof(unsigned int));
	
	if (in.eof()) return false;
	
	// Make char buffer for YAML config
	//
	auto cbuf = std::make_unique<char[]>(hsize+1);
	
	// Read YAML string
	//
	in.read(cbuf.get(), hsize);
	cbuf[hsize] = '\0';
	buf = std::string(cbuf.get());
	
	YAML::Node node = YAML::Load(buf.c_str());
	
	// Get parameters
	//
	lmax  = node["lmax"  ].as<int>();
	nmax  = node["nmax"  ].as<int>();
	time  = node["time"  ].as<double>();
	scale = node["scale" ].as<double>();

	// Optional parameters
	//
	if (node["geom"]) geom  = node["geom"  ].as<std::string>();
	else              geom  = "sphere";
	if (node["id"  ]) id    = node["id"    ].as<std::string>();
	
	// Look for norm flag
	//
	if (node["normed"]) normed = node["normed"].as<bool>();

      } else {
	
	// Rewind file
	//
	in.clear();
	in.seekg(curpos);
	
	SphCoefHeader header;
	in.read((char *)&header, sizeof(SphCoefHeader));
	
	time   = header.tnow;
	nmax   = header.nmax;
	lmax   = header.Lmax;
	geom   = "sphere";
	id     = std::string(header.id);
	normed = false;
      }
	
      int ldim = (lmax+1)*(lmax+2)/2;
      store.resize(ldim*nmax);
      coefs = std::make_shared<coefType>(store.data(), ldim, nmax);
      for (int ir=0; ir<nmax; ir++) {
	for (int l=0, L=0; l<=lmax; l++) {
	  for (int m=0; m<=l; m++, L++) {
	    double re, im=0.0;
	    if (m==0) {
	      in.read((char *)&re, sizeof(double));
	    } else {
	      in.read((char *)&re, sizeof(double));
	      in.read((char *)&im, sizeof(double));
	    }
	    (*coefs)(L, ir) = {re, im};
	  }
	}
      } 
    } catch (std::istream::failure e) {
      if (not in.eof())
	std::cerr << "Exception reading coefficient file: "
		  << e.what() << std::endl;
      return false;
    }
    
    if (in.eof()) return false;
    
    // Apply prefactors to make _true_ normed coefficients
    //
    if (exp_type and not normed) {
      int k = 0;
      for (int l=0; l<=lmax; l++) {
	for (int m=0; m<=l; m++) {
	  double fac = sqrt( (0.5*l+0.25)/M_PI * 
			     exp(lgamma(1.0+l-m) - lgamma(1.0+l+m)) );
	  
	  if (m != 0) fac *= M_SQRT2;
	  
	  // Cosine terms
	  for (int ir=0; ir<nmax; ir++) (*coefs)(k, ir) *= fac;
	  k++;
	  
	  // Sine terms
	  if (m != 0) {
	    for (int ir=0; ir<nmax; ir++) (*coefs)(k, ir) *= fac;
	    k++;
	  }
	}
      }
    }
    
    return true;
  }

  bool SlabStruct::read(std::istream& in, bool exp_type, bool verbose)
  {
    std::cout << "SlabStruct: no native coefficient format for this class" << std::endl;
    return false;
  }
  

  bool CubeStruct::read(std::istream& in, bool exp_type, bool verbose)
  {
    std::cout << "CubeStruct: no native coefficient format for this class" << std::endl;
    return false;
  }
  

  bool TblStruct::read(std::istream& in, bool exp_type, bool verbose)
  {
    std::vector<double> row;
    std::string line;

    if (std::getline(in, line)) {

      std::istringstream sin(line);

      // Read first value as time
      //
      sin >> time;

      // If okay, read the rest of the line as columns
      //
      double val;
      while (1) {
	sin >> val;
	if (sin) row.push_back(val);
	else     break;
      }

      // Number of cols
      cols = row.size();
      
      store.resize(cols);
      coefs = std::make_shared<coefType>(store.data(), cols);
      for (int i=0; i<cols; i++) (*coefs)(i) = row[i];

      return true;
    } else {
      return false;
    }
  }

  bool TrajStruct::read(std::istream& in, bool exp_type, bool verbose)
  {
    std::vector<double> row;
    std::string line;

    if (std::getline(in, line)) {

      std::istringstream sin(line);

      // Read first value as time
      //
      sin >> time;

      // If okay, read the rest of the line as columns
      //
      double val;
      while (1) {
	sin >> val;
	if (sin) row.push_back(val);
	else     break;
      }

      // Number of cols
      int cols = row.size();
      
      if (traj*rank != cols)
	throw std::runtime_error("TrajStruct::read: data != traj*rank");

      store.resize(cols);
      coefs = std::make_shared<coefType>(store.data(), traj, rank);
      for (int i=0; i<cols; i++) (*coefs)(i) = row[i];

      return true;
    } else {
      return false;
    }
  }

  bool SphFldStruct::read(std::istream& in, bool exp_type, bool verbose)
  {
    std::cout << "SphFldStruct: no native coefficient format for this class"
	      << std::endl;
    return false;
  }
  
  bool CylFldStruct::read(std::istream& in, bool exp_type, bool verbose)
  {
    std::cout << "CylFldStruct: no native coefficient format for this class"
	      << std::endl;
    return false;
  }
  
}
// END namespace CoefClasses
