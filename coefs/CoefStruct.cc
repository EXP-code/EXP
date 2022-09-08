#include <iostream>
#include <cstring>
#include <cmath>

#include <yaml-cpp/yaml.h>	// YAML support


#include <coef.H>		// Old original header
#include <CoefStruct.H>		// New coefficient structure

namespace Coefs
{
  
  void CoefStruct::copyfields(std::shared_ptr<CoefStruct> ret)
  {
    if (buf) {
      size_t buflen = std::strlen(buf.get());
      ret->buf = std::shared_ptr<char[]>(new char [buflen+1]);
      std::strncpy(ret->buf.get(), buf.get(), buflen+1);
    }

    ret->geom  = geom;
    ret->id    = id;
    ret->time  = time;
    ret->coefs = coefs;
  }

  std::shared_ptr<CylStruct> CylStruct::deepcopy()
  {
    auto ret = std::make_shared<CylStruct>();

    copyfields(ret);

    ret->mmax  = mmax;
    ret->nmax  = nmax;

    return ret;
  }

  std::shared_ptr<SphStruct> SphStruct::deepcopy()
  {
    auto ret = std::make_shared<SphStruct>();

    copyfields(ret);

    ret->lmax   = lmax;
    ret->nmax   = nmax;
    ret->scale  = scale;
    ret->normed = normed;

    return ret;
  }

  std::shared_ptr<TblStruct> TblStruct::deepcopy()
  {
    auto ret = std::make_shared<TblStruct>();

    copyfields(ret);

    ret->cols = cols;

    return ret;
  }

  bool CylStruct::read(std::istream& in, bool verbose)
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
	
	// Make and read char buffer for YAML config
	//
	buf = std::shared_ptr<char[]>(new char [ssize+1]);
	in.read(buf.get(), ssize);
	buf[ssize] = 0;		// Null terminate
	
	YAML::Node node = YAML::Load(buf.get());
	
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
  
  bool SphStruct::read(std::istream& in, bool exp_type)
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
	
	// Make and read char buffer for YAML config
	//
	buf = std::shared_ptr<char[]>(new char [hsize+1]);
	
	// Read YAML string
	//
	in.read(buf.get(), hsize);
	buf[hsize] = 0;		// Null terminate
	
	YAML::Node node = YAML::Load(buf.get());
	
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
      } else {
	throw std::runtime_error("no Spherical coeffcient file signature");
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

  bool TblStruct::read(std::istream& in, bool exp_type)
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

}
// END namespace Coefs
