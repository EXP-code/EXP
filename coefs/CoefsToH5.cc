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

#include <CoefsToH5.H>

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
      
      coefs[c->header.tnow] = c;
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
  

void SphH5::WriteH5Times(HighFive::Group& snaps)
{
  for (auto c : coefs) {
    auto C = c.second;

    std::ostringstream stim;
    stim << std::setprecision(8) << c.second->header.tnow;

    // Make a new group for this time
    //
    HighFive::Group stanza = snaps.createGroup(stim.str());

    // Add a time attribute
    //
    HighFive::Attribute t = stanza.createAttribute<double>("Time", HighFive::DataSpace::From(C->header.tnow));

    t.write(C->header.tnow);

    // Index counters
    //
    unsigned I = 0, J = 0;

    // Store the data
    //
    for (int l=0; l<=C->header.Lmax; l++) {

      for (int m=0; m<=l; m++) {

	std::vector<int> key = {l, m};
	
	Eigen::VectorXcd out(C->header.nmax);
	if (m) {
	  for (int n=0; n<C->header.nmax; n++)
	    out(n) = {C->coefs(I, n), C->coefs(I+1, n)};
	  I += 2;
	} else {
	  for (int n=0; n<C->header.nmax; n++)
	    out(n) = {C->coefs(I, n), 0.0};
	  I += 1;
	}

	std::ostringstream sout; sout << J++;
	HighFive::DataSet dataset = stanza.createDataSet(sout.str(), out);
	
	HighFive::Attribute v = dataset.createAttribute<int>
	  ("Index", HighFive::DataSpace::From(key));
      
	v.write(key);
      }
    }
  }
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


void CylH5::readNativeCoefs(const std::string& file)
{
  std::ifstream in(file);

  if (not in) {
    throw std::runtime_error("CylH5 ERROR (runtime) opening file <" + file + ">");
  }

  while (in) {
    CylCoefsPtr c = std::make_shared<CylCoefs>();
    if (not c->read(in, verbose)) break;

    coefs[c->time] = c;
  }
}
  

void CylH5::WriteH5Times(HighFive::Group& snaps)
{
  for (auto c : coefs) {
    std::ostringstream stim;
    stim << std::setprecision(8) << c.second->time;

    HighFive::Group stanza = snaps.createGroup(stim.str());

    HighFive::Attribute t = stanza.createAttribute<double>("Time", HighFive::DataSpace::From(c.second->time));

    t.write(c.second->time);

    auto C = c.second;

    for (int m=0; m<=C->mmax; m++) {

      std::vector<int> key = {m};
	
      Eigen::VectorXcd out(C->nmax);
      if (m) {
	for (int n=0; n<C->nmax; n++)
	  out(n) = {C->cos_c[m][n], C->sin_c[m][n]};
      } else {
	for (int n=0; n<C->nmax; n++)
	  out(n) = {C->cos_c[m][n], 0.0};
      }

      std::ostringstream sout; sout << m;
      HighFive::DataSet dataset = stanza.createDataSet(sout.str(), out);
      
      HighFive::Attribute v = dataset.createAttribute<int>
	("Index", HighFive::DataSpace::From(key));
      
      v.write(key);
    }
  }
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


void Coefs::WriteH5Coefs(const std::string& h5file)
{
  try {
    // Create a new hdf5 file
    //
    HighFive::File file(h5file + ".h5",
			HighFive::File::ReadWrite |
			HighFive::File::Create);

    // We write the coefficient file type
    //
    HighFive::Attribute type = file.createAttribute<std::string>("type", HighFive::DataSpace::From(coefType));

    type.write(coefType);
    
    // Stash the basis configuration
    //
    std::string config(getYAML());
    HighFive::Attribute yaml = file.createAttribute<std::string>("config", HighFive::DataSpace::From(config));

    yaml.write(config);

    double time = 3.1417;
    
    // Create a new group for coefficient snapshots
    //
    HighFive::Group group = file.createGroup("snapshots");

    // Write the coefficients
    //
    WriteH5Times(group);


  } catch (HighFive::Exception& err) {
    std::cerr << err.what() << std::endl;
  }
  
}

