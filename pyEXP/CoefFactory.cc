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

#include <CoefFactory.H>

namespace Coefs
{
  
  SphCoefs::SphCoefs(HighFive::File& file, bool verbose) :
    Coefs("Sphere", verbose)
  {
    std::string config, forceID;
    unsigned count;
    int lmax, nmax;
    double scale;
    
    file.getAttribute("name"   ).read(name   );
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
      auto coef = std::make_shared<SphStruct>();
      
      coef->lmax  = lmax;
      coef->nmax  = nmax;
      coef->time  = Time;
      coef->scale = scale;
      coef->coefs = in;
      
      coefs[roundTime(Time)] = coef;
    }
  }
  
  Eigen::MatrixXcd& SphCoefs::operator()(double time)
  {
    auto it = coefs.find(roundTime(time));
    if (it == coefs.end()) {
      
      mat.resize(0, 0);
      
    } else {
      
      mat = it->second->coefs;
      
    }
    
    return mat;
  }
  
  
  std::vector<Key> SphCoefs::makeKeys(Key k)
  {
    std::vector<Key> ret;
    if (coefs.size()==0) return ret;

    int lmax = coefs.begin()->second->lmax;
    int nmax = coefs.begin()->second->nmax;

    // Sanity
    if (k.size()) {
      k[0] = std::max<unsigned>(k[0], 0);
      k[0] = std::min<unsigned>(k[0], lmax);
    }

    if (k.size()>1) {
      k[1] = std::max<unsigned>(k[1], 0);
      k[1] = std::min<unsigned>(k[1], k[0]);
    }

    // Three options
    // 1. return all nkeys for a fixed l, m
    // 2. return all m, n for a fixed l
    // 3. return all keys

    // Option 3
    if (k.size()==0) {
      for (unsigned l=0; l<=lmax; l++)
	for (unsigned m=0; m<=l; m++) 
	  for (unsigned n=0; n<nmax; n++) ret.push_back({l, m, n});
    }
    // Option 2
    else if (k.size()==1) {
      for (unsigned m=0; m<=k[0]; m++) 
	for (unsigned n=0; n<nmax; n++) ret.push_back({k[0], m, n});
    }
    // Option 1
    else if (k.size()==2) {
      for (unsigned n=0; n<nmax; n++) ret.push_back({k[0], k[1], n});
    }
    // Bad sub key?
    else {
      throw std::runtime_error
	("SphCoefs::makeKeys: the subkey must have rank 0, 1 or 2");
    }

    return ret;
  }


  void SphCoefs::readNativeCoefs(const std::string& file)
  {
    std::ifstream in(file);
    
    if (not in) {
      throw std::runtime_error("SphCoefs ERROR (runtime) opening file <" + file + ">");
    }
    
    while (in) {
      try {
	SphStrPtr c = std::make_shared<SphStruct>();
	if (not c->read(in, verbose)) break;
	
	coefs[roundTime(c->time)] = c;
      }
      catch(std::runtime_error& error) {
	std::cout << "SphCoefs ERROR (runtime): " << error.what() << std::endl;
	break;
      }
      catch(std::logic_error& error) {
	std::cout << "SphCoefs ERROR (logic): " << error.what() << std::endl;
	break;
      }
    }
  }
  
  
  void SphCoefs::WriteH5Params(HighFive::File& file)
  {
    int lmax     = coefs.begin()->second->lmax;
    int nmax     = coefs.begin()->second->nmax;
    double scale = coefs.begin()->second->scale;
    
    std::string forceID(coefs.begin()->second->id);
    
    file.createAttribute<int>("lmax", HighFive::DataSpace::From(lmax)).write(lmax);
    file.createAttribute<int>("nmax", HighFive::DataSpace::From(nmax)).write(nmax);
    file.createAttribute<double>("scale", HighFive::DataSpace::From(scale)).write(scale);
    file.createAttribute<std::string>("forceID", HighFive::DataSpace::From(forceID)).write(forceID);
  }
  
  unsigned SphCoefs::WriteH5Times(HighFive::Group& snaps, unsigned count)
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
      stanza.createAttribute<double>("Time", HighFive::DataSpace::From(C->time)).write(C->time);
      
      // Index counters
      //
      unsigned I = 0, L = 0;
      
      // Pack the data into an Eigen matrix
      //
      int lmax = C->lmax, nmax= C->nmax;
      HighFive::DataSet dataset = stanza.createDataSet("coefficients", C->coefs);
    }
    
    return count;
  }
  
  
  std::string SphCoefs::getYAML()
  {
    std::string ret;
    if (coefs.size()) {
      if (coefs.begin()->second->buf.get())
	ret = std::string(coefs.begin()->second->buf.get());
    }
    return ret;
  }
  
  void SphCoefs::dump(int lmin, int lmax, int nmin, int nmax)
  {
    for (auto c : coefs) {
      unsigned I = 0;
      if (lmin>0) I += lmin*lmin;
      
      for (int ll=lmin; ll<=std::min<int>(lmax, c.second->lmax); ll++) {
	for (int mm=0; mm<=ll; mm++) {
	  std::cout << std::setw(18) << c.first << std::setw(5) << ll << std::setw(5) << mm << std::setw(5);
	  for (int nn=std::max<int>(nmin, 0); nn<std::min<int>(nmax, c.second->nmax); nn++) 
	    std::cout << std::setw(18) << c.second->coefs(I, nn);
	  std::cout << std::endl;
	  
	} // M loop
	
      } // L loop
      
    } // T loop
  }
  
  
  bool SphCoefs::CompareStanzas(CoefsPtr check)
  {
    bool ret = true;
    
    auto other = dynamic_cast<SphCoefs*>(check.get());
    
    // Check that every time in this one is in the other
    for (auto v : coefs) {
      if (other->coefs.find(roundTime(v.first)) == other->coefs.end()) {
	std::cout << "Can't find Time=" << v.first << std::endl;
	ret = false;
      }
    }
    
    if (not ret) {
      std::cout << "Times in other coeffcients are:";
      for (auto v : other->Times()) std::cout << " " << v;
      std::cout << std::endl;
    }

    if (ret) {
      std::cout << "Times are the same, now checking parameters at each time"
		<< std::endl;
      for (auto v : coefs) {
	auto it = other->coefs.find(v.first);
	if (v.second->lmax != it->second->lmax) ret = false;
	if (v.second->nmax != it->second->nmax) ret = false;
	if (v.second->time != it->second->time) ret = false;
      }
    }
    
    if (ret) {
      std::cout << "Parameters are the same, now checking coefficients"
		<< std::endl;
      for (auto v : coefs) {
	auto it = other->coefs.find(v.first);
	for (int i=0; i<v.second->coefs.rows(); i++) {
	  for (int j=0; j<v.second->coefs.cols(); j++) {
	    if (v.second->coefs(i, j) != it->second->coefs(i, j)) {
	      std::cout << "Coefficient (" << i << ", " << j << ")  "
			<< v.second->coefs(i, j) << " != "
			<< it->second->coefs(i, j) << std::endl;
	      ret = false;
	    }
	  }
	}
      }
    }
    
    return ret;
  }
  
  Eigen::MatrixXd& SphCoefs::Power()
  {
    if (coefs.size()) {
      
      int lmax = coefs.begin()->second->lmax;
      power.resize(coefs.size(), lmax+1);
      power.setZero();
      
      int T=0;
      for (auto v : coefs) {
	for (int l=0, L=0; l<=lmax; l++) {
	  for (int m=0; m<=l; m++, L++) {
	    auto rad = v.second->coefs.row(L);
	    power(T, l) += (rad.conjugate() * rad.transpose()).real()(0,0);
	  }
	}
	T++;
      }
    } else {
      power.resize(0, 0);
    }
    
    return power;
  }
  
  void SphCoefs::add(CoefStrPtr coef)
  {
    coefs[roundTime(coef->time)] = std::dynamic_pointer_cast<SphStruct>(coef);
  }

  CylCoefs::CylCoefs(HighFive::File& file, bool verbose) :
    Coefs("Cylinder", verbose)
  {
    int mmax, nmax;
    unsigned count;
    std::string config;
    
    file.getAttribute("name"   ).read(name  );
    file.getAttribute("mmax"   ).read(mmax  );
    file.getAttribute("nmax"   ).read(nmax  );
    file.getAttribute("config" ).read(config);
    file.getDataSet  ("count"  ).read(count );
    
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
      auto coef = std::make_shared<CylStruct>();
      
      coef->mmax  = mmax;
      coef->nmax  = nmax;
      coef->time  = Time;
      coef->coefs = in;
      
      coefs[roundTime(Time)] = coef;
    }
  }
  
  Eigen::MatrixXcd& CylCoefs::operator()(double time)
  {
    auto it = coefs.find(roundTime(time));
    
    if (it == coefs.end()) {
      
      mat.resize(0, 0);
      
    } else {
      
      mat = it->second->coefs;
      
    }
    
    return mat;
  }
  
  
  std::vector<Key> CylCoefs::makeKeys(Key k)
  {
    std::vector<Key> ret;
    if (coefs.size()==0) return ret;

    int mmax = coefs.begin()->second->mmax;
    int nmax = coefs.begin()->second->nmax;

    // Sanity check
    if (k.size()==1) {
      k[0] = std::max<unsigned>(k[0], 0);
      k[0] = std::min<unsigned>(k[0], mmax);
    }

    // Two options:
    // 1. return all keys with a fixed m
    // 2. return all keys
    //
    if (k.size()==0) {
      for (unsigned m=0; m<=mmax; m++)
	for (unsigned n=0; n<nmax; n++) ret.push_back({m, n});

    }
    else if (k.size()==1) {
      for (unsigned n=0; n<nmax; n++) ret.push_back({k[0], n});
    }
    // Bad sub key?
    else {
      throw std::runtime_error
	("SphCoefs::makeKeys: the subkey must have rank 1");
    }

    return ret;
  }

  void CylCoefs::readNativeCoefs(const std::string& file)
  {
    std::ifstream in(file);
    
    if (not in) {
      throw std::runtime_error("CylCoefs ERROR (runtime) opening file <" + file + ">");
    }
    
    while (in) {
      CylStrPtr c = std::make_shared<CylStruct>();
      if (not c->read(in, verbose)) break;
      
      coefs[roundTime(c->time)] = c;
    }
  }
  
  void CylCoefs::WriteH5Params(HighFive::File& file)
  {
    int mmax = coefs.begin()->second->mmax;
    int nmax = coefs.begin()->second->nmax;
    
    file.createAttribute<int>("mmax", HighFive::DataSpace::From(mmax)).write(mmax);
    file.createAttribute<int>("nmax", HighFive::DataSpace::From(nmax)).write(nmax);
  }
  
  unsigned CylCoefs::WriteH5Times(HighFive::Group& snaps, unsigned count)
  {
    for (auto c : coefs) {
      auto C = c.second;
      
      std::ostringstream stim;
      stim << std::setw(8) << std::setfill('0') << std::right << count++;
      HighFive::Group stanza = snaps.createGroup(stim.str());
      
      stanza.createAttribute<double>("Time", HighFive::DataSpace::From(C->time)).write(C->time);
      
      HighFive::DataSet dataset = stanza.createDataSet("coefficients", C->coefs);
    }
    
    return count;
  }
  
  
  std::string CylCoefs::getYAML()
  {
    std::string ret;
    if (coefs.size()) {
      if (coefs.begin()->second->buf.get())
	ret = std::string(coefs.begin()->second->buf.get());
    }

    return ret;
  }
  
  void CylCoefs::dump(int mmin, int mmax, int nmin, int nmax)
  {
    
    for (auto c : coefs) {
      for (int mm=mmin; mm<=std::min<int>(mmax, c.second->mmax); mm++) {
	std::cout << std::setw(18) << c.first << std::setw(5) << mm;
	for (int nn=std::max<int>(nmin, 0); nn<std::min<int>(nmax, c.second->nmax); nn++) {
	  if (angle)
	    std::cout << std::setw(18) << abs(c.second->coefs(mm, nn))
		      << std::setw(18) << arg(c.second->coefs(mm, nn));
	  else
	    std::cout << std::setw(18) << c.second->coefs(mm, nn).real()
		      << std::setw(18) << c.second->coefs(mm, nn).imag();
	}
	std::cout << std::endl;
      }
      // M loop
    }
    // T loop
    
  }
  
  Eigen::MatrixXd& CylCoefs::Power()
  {
    if (coefs.size()) {
      
      int mmax = coefs.begin()->second->mmax;
      
      power.resize(coefs.size(), mmax+1);
      power.setZero();
      
      int T=0;
      for (auto v : coefs) {
	for (int m=0; m<=mmax; m++) {
	  auto rad = v.second->coefs.row(m);
	  power(T, m) += (rad.conjugate() * rad.transpose()).real()(0, 0);
	}
	T++;
      }
    } else {
      power.resize(0, 0);
    }
    
    return power;
  }
  
  TableData::TableData(const std::vector<double>& times,
		       const std::vector<std::vector<double>>& data,
		       bool verbose) :
    times(times), data(data), Coefs("Table", verbose)
  {
    for (int i=0; i<times.size(); i++) {
      TblStrPtr c = std::make_shared<TblStruct>();
      c->cols = data[i].size();
      c->coefs.resize(1, c->cols);
      for (int j=0; j<c->cols; j++) c->coefs(0, j) = data[i][j];
      coefs[roundTime(c->time)] = c;
    }
  }

  TableData::TableData(std::string& file, bool verbose) :
    Coefs("Table", verbose)
  {
    std::ifstream in(file);
    
    if (not in) {
      throw std::runtime_error("TableData ERROR (runtime) opening file <" + file + ">");
    }
    
    while (in) {
      TblStrPtr c = std::make_shared<TblStruct>();
      if (not c->read(in, verbose)) break;
      
      coefs[roundTime(c->time)] = c;
    }

    for (auto c : coefs) times.push_back(c.first);
  }
    

  TableData::TableData(HighFive::File& file, bool verbose) :
    Coefs("Table", verbose)
  {
    int cols;
    unsigned count;
    std::string config;
    
    file.getAttribute("name"  ).read(name);
    file.getAttribute("cols"  ).read(cols);
    file.getAttribute("config").read(config);
    file.getDataSet  ("count" ).read(count);

    auto snaps = file.getGroup("snapshots");
    
    snaps.getDataSet("times").read(times);
    snaps.getDataSet("datatable").read(data);
      
    for (unsigned n=0; n<count; n++) {
      
      // Pack the data into the coefficient variable
      //
      auto coef = std::make_shared<TblStruct>();
      coef->cols  = cols;
      coef->time  = times[n];
      coef->coefs.resize(1, cols);
      for (int i=0; i<cols; i++) coef->coefs(0, i) = data[n][i];

      coefs[roundTime(times[n])] = coef;
    }
  }
  
  Eigen::MatrixXcd& TableData::operator()(double time)
  {
    auto it = coefs.find(roundTime(time));
    
    if (it == coefs.end()) {
      
      mat.resize(0, 0);
      
    } else {
      
      mat = it->second->coefs;
      
    }
    
    return mat;
  }
  
  void TableData::readNativeCoefs(const std::string& file)
  {
    std::ifstream in(file);
    
    if (not in) {
      throw std::runtime_error("TableData ERROR (runtime) opening file <" + file + ">");
    }
    
    while (in) {
      TblStrPtr c = std::make_shared<TblStruct>();
      if (not c->read(in, verbose)) break;
      
      coefs[roundTime(c->time)] = c;
    }
  }
  
  
  void TableData::WriteH5Params(HighFive::File& file)
  {
    int cols = coefs.begin()->second->cols;
    
    file.createAttribute<int>("cols", HighFive::DataSpace::From(cols)).write(cols);
  }
  
  unsigned TableData::WriteH5Times(HighFive::Group& snaps, unsigned count)
  {
    snaps.createDataSet("times", times);
    snaps.createDataSet("datatable", data);
    
    count = times.size();

    return count;
  }
  
  
  std::string TableData::getYAML()
  {
    std::string ret;
    if (coefs.size()) {
      if (coefs.begin()->second->buf.get())
	ret = std::string(coefs.begin()->second->buf.get());
    }

    return ret;
  }

  bool TableData::CompareStanzas(std::shared_ptr<Coefs> check)
  {
    bool ret = true;
    
    auto other = dynamic_cast<TableData*>(check.get());
    
    // Check that every time in this one is in the other
    for (auto v : coefs) {
      if (other->coefs.find(roundTime(v.first)) == other->coefs.end()) {
	std::cout << "Can't find Time=" << v.first << std::endl;
	ret = false;
      }
    }
    
    if (ret) {
      std::cout << "Times are the same, now checking parameters at each time"
		<< std::endl;
      for (auto v : coefs) {
	auto it = other->coefs.find(v.first);
	if (v.second->cols != it->second->cols) ret = false;
      }
    }
    
    if (ret) {
      std::cout << "Parameters are the same, now checking coefficients time"
		<< std::endl;
      for (auto v : coefs) {
	auto it = other->coefs.find(v.first);
	for (int m=0; m<v.second->cols; m++) {
	  if (v.second->coefs(m) != it->second->coefs(m)) {
	    ret = false;
	  }
	}
      }
    }
    
    return ret;
  }
  

  std::shared_ptr<Coefs> Coefs::factory(const std::string& file)
  {
    std::shared_ptr<Coefs> coefs;
    
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
	  coefs = std::make_shared<SphCoefs>(h5file);
	} else if (Type.compare("Cylinder")==0) {
	  coefs = std::make_shared<CylCoefs>(h5file);
	} else if (Type.compare("Table")==0) {
	  coefs = std::make_shared<TableData>(h5file);
	} else {
	  throw std::runtime_error("CoefFactory: unknown H5 coefficient file type");
	}
      } catch (HighFive::Exception& err) {
	std::cerr << "**** Error reading H5 file ****" << std::endl;
	std::cerr << err.what() << std::endl;
	exit(-1);
      }
      
      return coefs;
      
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
      coefs = std::make_shared<SphCoefs>();
    } else if (tmagic==cyl_magic) {
      coefs = std::make_shared<CylCoefs>();
    } else {
      coefs = std::make_shared<TableData>();
    }
    
    coefs->readNativeCoefs(file);
    
    return coefs;
  }
  
  std::shared_ptr<Coefs> Coefs::makecoefs(CoefStrPtr coef, std::string name)
  {
    std::shared_ptr<Coefs> ret;
    if (dynamic_cast<SphStruct*>(coef.get())) {
      ret = std::make_shared<SphCoefs>();
    } else if (dynamic_cast<CylStruct*>(coef.get())) {
      ret = std::make_shared<CylCoefs>();
    } else if (dynamic_cast<TblStruct*>(coef.get())) {
      ret = std::make_shared<TableData>();
    } else {
      throw std::runtime_error("CoefFactory: cannot deduce coefficient file type");
    }

    ret->setName(name);

    return ret;
  }
  

  std::shared_ptr<Coefs> Coefs::addcoef
  (std::shared_ptr<Coefs> coefs, CoefStrPtr coef)
  {
    std::shared_ptr<Coefs> ret = coefs;
    if (not coefs) ret = makecoefs(coef);

    ret->add(coef);
    
    return ret;
  }
  
  bool CylCoefs::CompareStanzas(std::shared_ptr<Coefs> check)
  {
    bool ret = true;
    
    auto other = dynamic_cast<CylCoefs*>(check.get());
    
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
	  for (int n=0; n<v.second->nmax; n++) { 
	    if (v.second->coefs(m, n) != it->second->coefs(m, n)) {
	      ret = false;
	    }
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
      
      // We write the coefficient mnemonic
      //
      file.createAttribute<std::string>("name", HighFive::DataSpace::From(name)).write(name);
      
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
  
  void CylCoefs::add(CoefStrPtr coef)
  {
    coefs[roundTime(coef->time)] = std::dynamic_pointer_cast<CylStruct>(coef);
  }

  void TableData::add(CoefStrPtr coef)
  {
    coefs[roundTime(coef->time)] = std::dynamic_pointer_cast<TblStruct>(coef);
  }


}
// END namespace Coefs
