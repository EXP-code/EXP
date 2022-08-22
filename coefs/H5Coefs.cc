#include <iostream>
#include <fstream>
#include <sstream>
#include <iomanip>
#include <complex>

#include "config.h"

#include <Eigen/Dense>

#include <highfive/H5File.hpp>
#include <highfive/H5DataSet.hpp>
#include <highfive/H5DataSpace.hpp>
#include <highfive/H5Attribute.hpp>

#include <H5Coefs.H>

SphH5::SphH5(HighFive::File& file, bool verbose) : Coefs("Sphere", verbose)
{
  std::string config, forceID;
  unsigned count;
  int lmax, nmax;
  double scale;

  file.getAttribute("lmax"   ).read(lmax   );
  file.getAttribute("nmax"   ).read(nmax   );
  file.getAttribute("scale"  ).read(scale  );
  file.getAttribute("config" ).read(config );
  file.getDataSet  ("count"  ).read(count  );
  file.getAttribute("forceID").read(forceID);

  // Open the snapshot group
  //
  auto snaps = file.getGroup("snapshots");

  for (unsigned n=0; n<count; n++) {

    std::ostringstream sout;
    sout << std::setw(8) << std::setfill('0') << std::right << n;

    auto stanza = snaps.getGroup(sout.str());
    
    double Time;
    stanza.getAttribute("Time").read(Time);

    Eigen::MatrixXcd in((lmax+1)*(lmax+2)/2, nmax);
    stanza.getDataSet("coefficients").read(in);

    // Pack the data into the coefficient variable
    //
    auto coef = std::make_shared<SphCoefs>();

    coef->header.Lmax  = lmax;
    coef->header.nmax  = nmax;
    coef->header.tnow  = Time;
    coef->header.scale = scale;
    std::strncpy(coef->header.id, forceID.c_str(), 64);

    // Store the data
    //
    coef->coefs.resize((lmax+1)*(lmax+1), nmax);

    for (int l=0, L=0, I=0; l<=lmax; l++) {
      for (int m=0; m<=l; m++, L++) {
	for (int n=0; n<nmax; n++) {
	  coef->coefs(I, n) = in(L, n).real();
	  if (m) coef->coefs(I+1, n) = in(L, n).imag();
	}
	if (m==0) I += 1;
	else      I += 2;
      }
    }

    coefs[roundTime(Time)] = coef;
  }
}

Eigen::MatrixXcd& SphH5::operator()(double time)
{
  auto it = coefs.find(roundTime(time));
  if (it == coefs.end()) {

    blank.resize(0, 0);

  } else {

    auto C = it->second;

    int lmax = C->header.Lmax, nmax= C->header.nmax;
    blank.resize((lmax+1)*(lmax+2)/2, nmax);

    for (int l=0, I=0, L=0; l<=C->header.Lmax; l++) {

      for (int m=0; m<=l; m++) {
      
	if (m) {
	  for (int n=0; n<C->header.nmax; n++)
	    blank(L, n) = {C->coefs(I, n), C->coefs(I+1, n)};
	  I += 2;
	} else {
	  for (int n=0; n<C->header.nmax; n++)
	    blank(L, n) = {C->coefs(I, n), 0.0};
	  I += 1;
	}
	L += 1;
	
      }
    }
  }

  return blank;
}

void SphH5::readNativeCoefs(const std::string& file)
{
  std::ifstream in(file);

  if (not in) {
    throw std::runtime_error("SphH5 ERROR (runtime) opening file <" + file + ">");
  }

  while (in) {
    try {
      SphCoefsPtr c = std::make_shared<SphCoefs>();
      if (not c->read(in, verbose)) break;
      
      coefs[roundTime(c->header.tnow)] = c;
    }
    catch(std::runtime_error& error) {
      std::cout << "SphH5 ERROR (runtime): " << error.what() << std::endl;
      break;
    }
    catch(std::logic_error& error) {
      std::cout << "SphH5 ERROR (logic): " << error.what() << std::endl;
      break;
    }
  }
}
  

void SphH5::WriteH5Params(HighFive::File& file)
{
  int lmax     = coefs.begin()->second->header.Lmax;
  int nmax     = coefs.begin()->second->header.nmax;
  double scale = coefs.begin()->second->header.scale;

  std::string forceID(coefs.begin()->second->header.id);

  file.createAttribute<int>("lmax", HighFive::DataSpace::From(lmax)).write(lmax);
  file.createAttribute<int>("nmax", HighFive::DataSpace::From(nmax)).write(nmax);
  file.createAttribute<double>("scale", HighFive::DataSpace::From(scale)).write(scale);
  file.createAttribute<std::string>("forceID", HighFive::DataSpace::From(forceID)).write(forceID);
}

unsigned SphH5::WriteH5Times(HighFive::Group& snaps, unsigned count)
{
  for (auto c : coefs) {
    auto C = c.second;

    std::ostringstream stim;
    stim << std::setw(8) << std::setfill('0') << std::right << count++;
    
    // Make a new group for this time
    //
    HighFive::Group stanza = snaps.createGroup(stim.str());

    // Add a time attribute
    //
    stanza.createAttribute<double>("Time", HighFive::DataSpace::From(C->header.tnow)).write(C->header.tnow);

    // Index counters
    //
    unsigned I = 0, L = 0;
    
    // Pack the data into an Eigen matrix
    //
    int lmax = C->header.Lmax, nmax= C->header.nmax;
    Eigen::MatrixXcd out((lmax+1)*(lmax+2)/2, nmax);

    // Store the data
    //
    for (int l=0; l<=C->header.Lmax; l++) {

      for (int m=0; m<=l; m++) {

	if (m) {
	  for (int n=0; n<C->header.nmax; n++)
	    out(L, n) = {C->coefs(I, n), C->coefs(I+1, n)};
	  I += 2;
	} else {
	  for (int n=0; n<C->header.nmax; n++)
	    out(L, n) = {C->coefs(I, n), 0.0};
	  I += 1;
	}
	L += 1;

      }
    }

    HighFive::DataSet dataset = stanza.createDataSet("coefficients", out);
  }

  return count;
}


std::string SphH5::getYAML()
{
  std::string ret;
  if (coefs.size()) ret = coefs.begin()->second->buf.get();
  std::cout << "YAML test: " << ret << std::endl;

  return ret;
}

void SphH5::dump(int lmin, int lmax, int nmin, int nmax)
{
  for (auto c : coefs) {
    unsigned I = 0;
    if (lmin>0) I += lmin*lmin;

    for (int ll=lmin; ll<=std::min<int>(lmax, c.second->header.Lmax); ll++) {
      for (int mm=0; mm<=ll; mm++) {
	int S = mm==0 ? 1 : 2;
	for (int s=0; s<S; s++) {
	  std::cout << std::setw(18) << c.first << std::setw(5) << ll << std::setw(5) << mm << std::setw(5) << s;
	  for (int nn=std::max<int>(nmin, 0); nn<std::min<int>(nmax, c.second->header.nmax); nn++) 
	    std::cout << std::setw(18) << c.second->coefs(I+s, nn);
	  std::cout << std::endl;
	  if (mm==0) break;
	}
	// Cosine/sine
	
	if (mm==0) I += 1;
	else       I += 2;
	
      } // M loop
      
    } // L loop
    
  } // T loop
}


bool SphH5::CompareStanzas(Coefs* check)
{
  bool ret = true;

  auto other = dynamic_cast<SphH5*>(check);

  // Check that every time in this one is in the other
  for (auto v : coefs) {
    if (other->coefs.find(v.first) == other->coefs.end()) {
      std::cout << "Can't find Time=" << v.first << std::endl;
      ret = false;
    }
  }

  if (ret) {
    std::cout << "Times are the same, now checking parameters at each time"
	      << std::endl;
    for (auto v : coefs) {
      auto it = other->coefs.find(v.first);
      if (v.second->header.Lmax != it->second->header.Lmax) ret = false;
      if (v.second->header.nmax != it->second->header.nmax) ret = false;
      if (v.second->header.tnow != it->second->header.tnow) ret = false;
    }
  }

  if (ret) {
    std::cout << "Parameters are the same, now checking coefficients time"
	      << std::endl;
    for (auto v : coefs) {
      auto it = other->coefs.find(v.first);
      for (int i=0; i<v.second->coefs.rows(); i++) 
	for (int j=0; j<v.second->coefs.cols(); j++) 
	  if (v.second->coefs(i, j) != it->second->coefs(i, j)) ret = false;
    }
  }

  return ret;
}

CylH5::CylH5(HighFive::File& file, bool verbose) : Coefs("Cylinder", verbose)
{
  int mmax, nmax;
  unsigned count;
  std::string config;
  
  file.getAttribute("mmax").read(mmax);
  file.getAttribute("nmax").read(nmax);
  file.getAttribute("config").read(config);
  file.getDataSet("count").read(count);

  // Open the snapshot group
  //
  auto snaps = file.getGroup("snapshots");

  for (unsigned n=0; n<count; n++) {

    std::ostringstream sout;
    sout << std::setw(8) << std::setfill('0') << std::right << n;

    auto stanza = snaps.getGroup(sout.str());
    
    double Time;
    stanza.getAttribute("Time").read(Time);

    Eigen::MatrixXcd in(mmax+1, nmax);
    stanza.getDataSet("coefficients").read(in);

    // Pack the data into the coefficient variable
    //
    auto coef = std::make_shared<CylCoefs>();

    coef->mmax  = mmax;
    coef->nmax  = nmax;
    coef->time  = Time;

    coef->cos_c.resize(mmax+1);
    coef->sin_c.resize(mmax+1);
    for (auto & v : coef->cos_c) v.resize(nmax);
    for (auto & v : coef->sin_c) v.resize(nmax);

    for (int m=0; m<=mmax; m++) {
      for (int n=0; n<nmax; n++) {
	coef->cos_c[m][n] = in(m, n).real();
	coef->sin_c[m][n] = in(m, n).imag();
      }
    }
    
    coefs[roundTime(Time)] = coef;
  }
}

Eigen::MatrixXcd& CylH5::operator()(double time)
{
  auto it = coefs.find(roundTime(time));

  if (it == coefs.end()) {

    blank.resize(0, 0);

  } else {
    
    auto C = it->second;
    blank.resize(C->mmax+1, C->nmax);

    for (int m=0; m<=C->mmax; m++) {

      if (m) {
	for (int n=0; n<C->nmax; n++)
	  blank(m, n) = {C->cos_c[m][n], C->sin_c[m][n]};
      } else {
	for (int n=0; n<C->nmax; n++)
	  blank(m, n) = {C->cos_c[m][n], 0.0};
      }

    }

  }

  return blank;
}

void CylH5::readNativeCoefs(const std::string& file)
{
  std::ifstream in(file);

  if (not in) {
    throw std::runtime_error("CylH5 ERROR (runtime) opening file <" + file + ">");
  }

  while (in) {
    CylCoefsPtr c = std::make_shared<CylCoefs>();
    if (not c->read(in, verbose)) break;

    coefs[roundTime(c->time)] = c;
  }
}
  

void CylH5::WriteH5Params(HighFive::File& file)
{
  int mmax = coefs.begin()->second->mmax;
  int nmax = coefs.begin()->second->nmax;

  file.createAttribute<int>("mmax", HighFive::DataSpace::From(mmax)).write(mmax);
  file.createAttribute<int>("nmax", HighFive::DataSpace::From(nmax)).write(nmax);
}

unsigned CylH5::WriteH5Times(HighFive::Group& snaps, unsigned count)
{
  for (auto c : coefs) {
    auto C = c.second;

    std::ostringstream stim;
    stim << std::setw(8) << std::setfill('0') << std::right << count++;
    HighFive::Group stanza = snaps.createGroup(stim.str());

    stanza.createAttribute<double>("Time", HighFive::DataSpace::From(C->time)).write(C->time);

    Eigen::MatrixXcd out(C->mmax+1, C->nmax);

    for (int m=0; m<=C->mmax; m++) {

      if (m) {
	for (int n=0; n<C->nmax; n++)
	  out(m, n) = {C->cos_c[m][n], C->sin_c[m][n]};
      } else {
	for (int n=0; n<C->nmax; n++)
	  out(m, n) = {C->cos_c[m][n], 0.0};
      }

    }

    HighFive::DataSet dataset = stanza.createDataSet("coefficients", out);
  }

  return count;
}


std::string CylH5::getYAML()
{
  std::string ret;
  if (coefs.size()) ret = coefs.begin()->second->buf.get();
  std::cout << "YAML test: " << ret << std::endl;

  return ret;
}

void CylH5::dump(int mmin, int mmax, int nmin, int nmax)
{

  for (auto c : coefs) {
    for (int mm=mmin; mm<=std::min<int>(mmax, c.second->mmax); mm++) {
      std::cout << std::setw(18) << c.first << std::setw(5) << mm;
      for (int nn=std::max<int>(nmin, 0); nn<std::min<int>(nmax, c.second->nmax); nn++) {
	if (mm==0) {
	  if (angle)
	    std::cout << std::setw(18) << 0.0;
	  else
	    std::cout << std::setw(18) << fabs(c.second->cos_c[mm][nn]);
	} else {
	  if (angle) {
	    double arg = atan2(c.second->sin_c[mm][nn], c.second->cos_c[mm][nn]);
	    std::cout << std::setw(18) << arg;
	  } else {
	    double amp =
	      c.second->cos_c[mm][nn] * c.second->cos_c[mm][nn] +
	      c.second->sin_c[mm][nn] * c.second->sin_c[mm][nn] ;
	    std::cout << std::setw(18) << sqrt(amp);
	  }
	}
      }
      std::cout << std::endl;
    }
    // M loop
  }
  // T loop

}

CoefClient::CoefClient(const std::string& file)
{
  // First attempt to read the file as H5
  //
  try {
    // Silence the HDF5 error stack
    //
    HighFive::SilenceHDF5 quiet;

    // Try opening the file as HDF5
    //
    HighFive::File h5file(file, HighFive::File::ReadOnly);
    
    // Get the coefficient file type
    //
    std::string Type;
    HighFive::Attribute type = h5file.getAttribute("type");
    type.read(Type);
    
    try {
      if (Type.compare("Sphere")==0) {
	coefs = std::make_shared<SphH5>(h5file);
      } else if (Type.compare("Cylinder")==0) {
	coefs = std::make_shared<CylH5>(h5file);
      } else {
	throw std::runtime_error("CoefClient: unknown H5 coefficient file type");
      }
    } catch (HighFive::Exception& err) {
      std::cerr << "**** Error reading H5 file ****" << std::endl;
      std::cerr << err.what() << std::endl;
      exit(-1);
    }
    
    return;

  } catch (HighFive::Exception& err) {
    std::cerr << "**** Error opening H5 file, will try other types ****" << std::endl;
    // std::cerr << err.what() << std::endl;
  }

  // Open file and read magic number
  std::ifstream in(file);

  in.exceptions ( std::istream::failbit | std::istream::badbit );

  // Attempt to read coefficient magic number
  //
  const unsigned int sph_magic = 0xc0a57a2;
  const unsigned int cyl_magic = 0xc0a57a3;
  unsigned int tmagic;

  in.read(reinterpret_cast<char*>(&tmagic), sizeof(unsigned int));
  in.close();
  
  if (tmagic==sph_magic) {
    coefs = std::make_shared<SphH5>();
  } else if (tmagic==cyl_magic) {
    coefs = std::make_shared<CylH5>();
  } else {
    throw std::runtime_error("CoefClient: unknown native coefficient file type");
  }

  coefs->readNativeCoefs(file);

}


bool CylH5::CompareStanzas(Coefs* check)
{
  bool ret = true;

  auto other = dynamic_cast<CylH5*>(check);

  // Check that every time in this one is in the other
  for (auto v : coefs) {
    if (other->coefs.find(v.first) == other->coefs.end()) {
      std::cout << "Can't find Time=" << v.first << std::endl;
      ret = false;
    }
  }

  if (ret) {
    std::cout << "Times are the same, now checking parameters at each time"
	      << std::endl;
    for (auto v : coefs) {
      auto it = other->coefs.find(v.first);
      if (v.second->mmax != it->second->mmax) ret = false;
      if (v.second->nmax != it->second->nmax) ret = false;
      if (v.second->time != it->second->time) ret = false;
    }
  }

  if (ret) {
    std::cout << "Parameters are the same, now checking coefficients time"
	      << std::endl;
    for (auto v : coefs) {
      auto it = other->coefs.find(v.first);
      for (int m=0; m<v.second->mmax; m++) {
	for (int n=0; n<v.second->nmax; n++) 
	  if (v.second->cos_c[m][n] != it->second->cos_c[m][n]) ret = false;

	if (m) {
	  for (int n=0; n<v.second->nmax; n++) 
	    if (v.second->sin_c[m][n] != it->second->sin_c[m][n]) ret = false;
	}
      }
    }
  }

  return ret;
}

void Coefs::WriteH5Coefs(const std::string& prefix)
{
  try {
    // Create a new hdf5 file
    //
    HighFive::File file(prefix + ".h5",
			HighFive::File::ReadWrite |
			HighFive::File::Create);

    // We write the coefficient file type
    //
    file.createAttribute<std::string>("type", HighFive::DataSpace::From(coefType)).write(coefType);
    
    // Stash the basis configuration (this is not yet implemented in EXP)
    //
    std::string config(getYAML());
    file.createAttribute<std::string>("config", HighFive::DataSpace::From(config)).write(config);

    // Write the specific parameters
    //
    WriteH5Params(file);

    // Group count variable
    //
    unsigned count = 0;
    HighFive::DataSet dataset = file.createDataSet("count", count);

    // Create a new group for coefficient snapshots
    //
    HighFive::Group group = file.createGroup("snapshots");

    // Write the coefficients
    //
    count = WriteH5Times(group, count);

    // Update the count
    //
    dataset.write(count);

  } catch (HighFive::Exception& err) {
    std::cerr << err.what() << std::endl;
  }
  
}

void Coefs::ExtendH5Coefs(const std::string& prefix)
{
  try {
    // Create a new hdf5 file
    //
    HighFive::File file(prefix + ".h5",	HighFive::File::ReadWrite);

    // Get the dataset
    HighFive::DataSet dataset = file.getDataSet("count");

    unsigned count;
    dataset.read(count);

    HighFive::Group group = file.getGroup("snapshots");

    // Write the coefficients
    //
    count = WriteH5Times(group, count);

    // Update the count
    //
    dataset.write(count);

  } catch (HighFive::Exception& err) {
    std::cerr << err.what() << std::endl;
  }
  
}

