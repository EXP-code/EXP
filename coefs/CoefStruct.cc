#include <iostream>
#include <cstring>
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
    ret->coefs = coefs;
    if (ret->ctr.size()) ret->ctr = ctr;
  }

  void CylStruct::create()
  {
    if (nmax>0)
      coefs.resize(mmax+1, nmax);
    else
      throw std::runtime_error("CylStruct::create: nmax must be >0");
  }

  void SphStruct::create()
  {
    if (nmax>0)
      coefs.resize((lmax+1)*(lmax+2)/2, nmax);
    else
      throw std::runtime_error("SphStruct::create: nmax must be >0");
  }

  void SlabStruct::create()
  {
    if (nmaxz>0)
      coefT.resize({2*nmaxx+1, 2*nmaxy+1, nmaxz});
    else
      throw std::runtime_error("labStruct::create: vertical order must be >0");
  }

  void BoxStruct::create()
  {
    coefT.resize({2*nmaxx+1, 2*nmaxy+1, 2*nmaxz+1});
  }

  void TblStruct::create()
  {
    if (cols>0)
      coefs.resize(1, cols);
    else
      throw std::runtime_error("TblStruct::create: cols must be >0");
  }

  std::shared_ptr<CoefStruct> CylStruct::deepcopy()
  {
    auto ret = std::make_shared<CylStruct>();

    copyfields(ret);

    ret->mmax  = mmax;
    ret->nmax  = nmax;

    return ret;
  }

  std::shared_ptr<CoefStruct> SlabStruct::deepcopy()
  {
    auto ret = std::make_shared<SlabStruct>();

    copyfields(ret);

    ret->nmaxx  = nmaxx;
    ret->nmaxy  = nmaxy;
    ret->nmaxz  = nmaxz;

    return ret;
  }

  std::shared_ptr<CoefStruct> BoxStruct::deepcopy()
  {
    auto ret = std::make_shared<BoxStruct>();

    copyfields(ret);

    ret->nmaxx  = nmaxx;
    ret->nmaxy  = nmaxy;
    ret->nmaxz  = nmaxz;

    return ret;
  }

  std::shared_ptr<CoefStruct> SphStruct::deepcopy()
  {
    auto ret = std::make_shared<SphStruct>();

    copyfields(ret);

    ret->lmax   = lmax;
    ret->nmax   = nmax;
    ret->scale  = scale;
    ret->normed = normed;

    return ret;
  }

  std::shared_ptr<CoefStruct> TblStruct::deepcopy()
  {
    auto ret = std::make_shared<TblStruct>();

    copyfields(ret);

    ret->cols = cols;

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
    
    coefs.resize(mmax+1, nmax);
    Eigen::VectorXd workR(nmax), workC(nmax);
    workC.setZero();
    
    try {
      for (int mm=0; mm<=mmax; mm++) {
	
	in.read((char *)workR.data(), sizeof(double)*nmax);
	if (mm) in.read((char *)workC.data(), sizeof(double)*nmax);
	
	coefs.row(mm).real() = workR;
	coefs.row(mm).imag() = workC;
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
	
      coefs.resize((lmax+1)*(lmax+2)/2, nmax);
	
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
	    coefs(L, ir) = {re, im};
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
	  for (int ir=0; ir<nmax; ir++) coefs(k, ir) *= fac;
	  k++;
	  
	  // Sine terms
	  if (m != 0) {
	    for (int ir=0; ir<nmax; ir++) coefs(k, ir) *= fac;
	    k++;
	  }
	}
      }
    }
    
    return true;
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
      
      coefs.resize(1, cols);
      for (int i=0; i<cols; i++) coefs(0, i) = row[i];

      return true;
    } else {
      return false;
    }
  }

  bool SlabStruct::read(std::istream& in, bool exp_type, bool verbose)
  { return false; }

  bool  BoxStruct::read(std::istream& in, bool exp_type, bool verbose)
  { return false; }


}
// END namespace CoefClasses
