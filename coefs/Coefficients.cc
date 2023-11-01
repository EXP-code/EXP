#include <filesystem>
#include <iterator>
#include <iostream>
#include <fstream>
#include <sstream>
#include <iomanip>
#include <complex>

#include <config_exp.h>
#include <localmpi.H>

#include <Eigen/Dense>
#include <yaml-cpp/yaml.h>

#include <highfive/H5File.hpp>
#include <highfive/H5DataSet.hpp>
#include <highfive/H5DataSpace.hpp>
#include <highfive/H5Attribute.hpp>

#include <Coefficients.H>

namespace CoefClasses
{
  
  void Coefs::copyfields(std::shared_ptr<Coefs> p)
  {
    // These variables will copy data, not pointers
    p->mat      = mat;
    p->power    = power;
    p->geometry = geometry;
    p->name     = name;
    p->verbose  = verbose;
    p->times    = times;
  }

  std::tuple<Eigen::MatrixXcd&, bool> Coefs::interpolate(double time)
  {
    bool onGrid = true;

    if (time < times.front()-deltaT or time > times.back()+deltaT) {
      
      const  int max_oab = 8;	// Allow 'slop' off grid attempts
      static int cnt_oab = 0;	// before triggering an off grid stop

      std::cerr << "Coefs::interpolate: time=" << time
	   << " is offgrid [" << times.front()
		<< ", " << times.back() << "] #" << ++cnt_oab << std::endl;
      
      if (cnt_oab > max_oab) onGrid = false;
    }

    auto it = std::lower_bound(times.begin(), times.end(), time);
    auto lo = it, hi = it;

    if (hi == std::prev(times.end()) or hi == times.end()) {
      hi = std::prev(times.end());
      lo = hi - 1;
    } else hi++;
    
    double A = (*hi - time)/(*hi - *lo);
    double B = (time - *lo)/(*hi - *lo);
    
    int iA = std::distance(times.begin(), lo);
    int iB = std::distance(times.begin(), hi);

    auto cA = getCoefStruct(times[iA]);
    auto cB = getCoefStruct(times[iB]);

    int rows = cA->coefs.rows(), cols = cA->coefs.cols();
    mat.resize(rows, cols);

    for (int c=0; c<rows; c++) {
      for (int n=0; n<cols; n++)
	mat(c, n) = A*(cA->coefs)(c, n) + B*(cB->coefs)(c, n);
    }

    return {mat, onGrid};
  }
  
  SphCoefs::SphCoefs(HighFive::File& file, int stride,
		     double Tmin, double Tmax, bool verbose) :
    Coefs("sphere", verbose)
  {
    std::string config, geometry, forceID;
    unsigned count;
    double scale;
    
    file.getAttribute("name"    ).read(name    );
    file.getAttribute("lmax"    ).read(Lmax    );
    file.getAttribute("nmax"    ).read(Nmax    );
    file.getAttribute("scale"   ).read(scale   );
    file.getAttribute("config"  ).read(config  );
    file.getDataSet  ("count"   ).read(count   );
    file.getAttribute("geometry").read(geometry);
    file.getAttribute("forceID" ).read(forceID );
    
    // Open the snapshot group
    //
    auto snaps = file.getGroup("snapshots");
    
    for (unsigned n=0; n<count; n+=stride) {
      
      std::ostringstream sout;
      sout << std::setw(8) << std::setfill('0') << std::right << n;
      
      auto stanza = snaps.getGroup(sout.str());
      
      double Time;
      stanza.getAttribute("Time").read(Time);
      
      // Check for center data
      //
      std::vector<double> ctr;
      if (stanza.hasAttribute("Center")) {
	stanza.getAttribute("Center").read(ctr);
      }

      if (Time < Tmin or Time > Tmax) continue;

      Eigen::MatrixXcd in((Lmax+1)*(Lmax+2)/2, Nmax);
      stanza.getDataSet("coefficients").read(in);
      
      // Pack the data into the coefficient variable
      //
      auto coef = std::make_shared<SphStruct>();
      
      if (ctr.size()) coef->ctr = ctr;

      coef->lmax  = Lmax;
      coef->nmax  = Nmax;
      coef->time  = Time;
      coef->scale = scale;
      coef->coefs = in;
      coef->geom  = geometry;
      coef->id    = forceID;
      
      coefs[roundTime(Time)] = coef;
    }

    times.clear();
    for (auto t : coefs) times.push_back(t.first);
  }
  
  std::shared_ptr<Coefs> SphCoefs::deepcopy()
  {
    auto ret = std::make_shared<SphCoefs>();

    // Copy the base-class fields
    copyfields(ret);

    // Copy the local structures from the map to the struct pointers
    // by copyfing fields, not the pointer
    for (auto v : coefs)
      ret->coefs[v.first] =
	std::dynamic_pointer_cast<SphStruct>(v.second->deepcopy());

    ret->Lmax = Lmax;
    ret->Nmax = Nmax;

    return ret;
  }

  std::shared_ptr<Coefs> CylCoefs::deepcopy()
  {
    auto ret = std::make_shared<CylCoefs>();

    // Copy the base-class fields
    copyfields(ret);

    // Copy the local structures from the map to the struct pointers
    // by copyfing fields, not the pointer
    for (auto v : coefs)
      ret->coefs[v.first] =
	std::dynamic_pointer_cast<CylStruct>(v.second->deepcopy());

    ret->Mmax  = Mmax;
    ret->Nmax  = Nmax;
    ret->angle = angle;

    return ret;
  }

  std::shared_ptr<Coefs> TableData::deepcopy()
  {
    auto ret = std::make_shared<TableData>();

    // Copy the base-class fields
    copyfields(ret);

    // Copy the local structures from the map to the struct pointers
    // by copyfing fields, not the pointer
    for (auto v : coefs)
      ret->coefs[v.first] =
	std::dynamic_pointer_cast<TblStruct>(v.second->deepcopy());

    ret->data  = data;
    ret->times = times;

    return ret;
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
  
  void SphCoefs::setMatrix(double time, const Eigen::MatrixXcd& mat)
  {
    auto it = coefs.find(roundTime(time));

    if (it == coefs.end()) {
      std::ostringstream str;
      str << "SphCoefs::setMatrix: requested time=" << time << " not found";
      throw std::runtime_error(str.str());
    } else {
      it->second->coefs = mat;
    }
  }
  
  Eigen::Tensor<std::complex<double>, 3> SphCoefs::getAllCoefs()
  {
    Eigen::Tensor<std::complex<double>, 3> ret;

    auto times = Times();
    int ntim = times.size();

    // Resize the tensor
    ret.resize((Lmax+1)*(Lmax+2)/2, Nmax, ntim);

    for (int t=0; t<ntim; t++) {
      auto cof = coefs[times[t]];
      for (int l=0; l<(Lmax+2)*(Lmax+1)/2; l++) {
	for (int n=0; n<Nmax; n++) {
	  ret(l, n, t) = cof->coefs(l, n);
	}
      }
    }

    return ret;
  }


  std::vector<Key> SphCoefs::makeKeys(Key k)
  {
    std::vector<Key> ret;
    if (coefs.size()==0) return ret;

    // Sanity
    if (k.size()) {
      k[0] = std::max<unsigned>(k[0], 0);
      k[0] = std::min<unsigned>(k[0], Lmax);
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
      for (unsigned l=0; l<=Lmax; l++)
	for (unsigned m=0; m<=l; m++) 
	  for (unsigned n=0; n<Nmax; n++) ret.push_back({l, m, n});
    }
    // Option 2
    else if (k.size()==1) {
      for (unsigned m=0; m<=k[0]; m++) 
	for (unsigned n=0; n<Nmax; n++) ret.push_back({k[0], m, n});
    }
    // Option 1
    else if (k.size()==2) {
      for (unsigned n=0; n<Nmax; n++) ret.push_back({k[0], k[1], n});
    }
    // Bad sub key?
    else {
      throw std::runtime_error
	("SphCoefs::makeKeys: the subkey must have rank 0, 1 or 2");
    }

    return ret;
  }


  void SphCoefs::readNativeCoefs(const std::string& file, int stride,
				 double tmin, double tmax)
  {
    std::ifstream in(file);
    
    if (not in) {
      throw std::runtime_error("SphCoefs ERROR (runtime) opening file <" + file + ">");
    }
    
    int count = 0;
    while (in) {
      try {
	SphStrPtr c = std::make_shared<SphStruct>();
	if (not c->read(in, verbose)) break;

	if (count++ % stride) continue;
	if (c->time < tmin or c->time > tmax) continue;

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

    if (coefs.size()) {
      Lmax = coefs.begin()->second->lmax;
      Nmax = coefs.begin()->second->nmax;
    }

    times.clear();
    for (auto t : coefs) times.push_back(t.first);

    if (myid==0)
      std::cerr << "---- Coefs::factory: "
		<< "read EXP native and created SphCoefs"
		<< std::endl;
  }
  
  
  void SphCoefs::WriteH5Params(HighFive::File& file)
  {
    double scale = coefs.begin()->second->scale;
    
    std::string forceID(coefs.begin()->second->id);
    
    file.createAttribute<int>("lmax", HighFive::DataSpace::From(Lmax)).write(Lmax);
    file.createAttribute<int>("nmax", HighFive::DataSpace::From(Nmax)).write(Nmax);
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
      
      // Add a center attribute
      //
      if (C->ctr.size()>0)
	stanza.createAttribute<double>("Center", HighFive::DataSpace::From(C->ctr)).write(C->ctr);
      
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
      ret = coefs.begin()->second->buf;
    }
    return ret;
  }
  
  void SphCoefs::dump(int lmin, int lmax, int nmin, int nmax)
  {
    for (auto c : coefs) {
      unsigned I = 0;
      if (lmin>0) I += lmin*lmin;
      
      for (int ll=lmin; ll<=std::min<int>(lmax, Lmax); ll++) {
	for (int mm=0; mm<=ll; mm++) {
	  std::cout << std::setw(18) << c.first << std::setw(5) << ll << std::setw(5) << mm << std::setw(5);
	  for (int nn=std::max<int>(nmin, 0); nn<std::min<int>(nmax, Nmax); nn++) 
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
  
  Eigen::MatrixXd& SphCoefs::Power(int min, int max)
  {
    if (coefs.size()) {
      
      int lmax = coefs.begin()->second->lmax;
      int nmax = coefs.begin()->second->nmax;
      power.resize(coefs.size(), lmax+1);
      power.setZero();
      
      int T=0;
      for (auto v : coefs) {
	for (int l=0, L=0; l<=lmax; l++) {
	  for (int m=0; m<=l; m++, L++) {
	    auto rad = v.second->coefs.row(L);
	    for (int n=std::max<int>(0, min); n<std::min<int>(nmax, max); n++) {
	      double val = std::abs(rad(n));
	      power(T, l) += val * val;
	    }
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
    auto p = std::dynamic_pointer_cast<SphStruct>(coef);
    Lmax = p->lmax;
    Nmax = p->nmax;
    coefs[roundTime(coef->time)] = p;
  }

  CylCoefs::CylCoefs(HighFive::File& file, int stride,
		     double Tmin, double Tmax, bool verbose) :
    Coefs("cylinder", verbose)
  {
    unsigned count;
    std::string config;
    
    file.getAttribute("name"   ).read(name  );
    file.getAttribute("mmax"   ).read(Mmax  );
    file.getAttribute("nmax"   ).read(Nmax  );
    file.getAttribute("config" ).read(config);
    file.getDataSet  ("count"  ).read(count );
    
    // Open the snapshot group
    //
    auto snaps = file.getGroup("snapshots");
    
    for (unsigned n=0; n<count; n+=stride) {
      
      std::ostringstream sout;
      sout << std::setw(8) << std::setfill('0') << std::right << n;
      
      auto stanza = snaps.getGroup(sout.str());
      
      double Time;
      stanza.getAttribute("Time").read(Time);
      
      // Check for center data
      //
      std::vector<double> ctr;
      if (stanza.hasAttribute("Center")) {
	stanza.getAttribute("Center").read(ctr);
      }

      if (Time < Tmin or Time > Tmax) continue;

      Eigen::MatrixXcd in(Mmax+1, Nmax);
      stanza.getDataSet("coefficients").read(in);
      
      // Pack the data into the coefficient variable
      //
      auto coef = std::make_shared<CylStruct>();
      
      if (ctr.size()) coef->ctr = ctr;

      coef->mmax  = Mmax;
      coef->nmax  = Nmax;
      coef->time  = Time;
      coef->coefs = in;
      
      coefs[roundTime(Time)] = coef;
    }

    times.clear();
    for (auto t : coefs) times.push_back(t.first);
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
  
  
  void CylCoefs::setMatrix(double time, const Eigen::MatrixXcd& mat)
  {
    auto it = coefs.find(roundTime(time));

    if (it == coefs.end()) {
      std::ostringstream str;
      str << "CylCoefs::setMatrix: requested time=" << time << " not found";
      throw std::runtime_error(str.str());
    } else {
      it->second->coefs = mat;
    }
  }

  Eigen::Tensor<std::complex<double>, 3> CylCoefs::getAllCoefs()
  {
    Eigen::Tensor<std::complex<double>, 3> ret;

    auto times = Times();

    int ntim = times.size();

    // Resize the tensor
    ret.resize(Mmax+1, Nmax, ntim);
    
    for (int t=0; t<ntim; t++) {
      auto cof = coefs[times[t]];
      for (int m=0; m<Mmax+1; m++) {
	for (int n=0; n<Nmax; n++) {
	  ret(m, n, t) = cof->coefs(m, n);
	}
      }
    }

    return ret;
  }

  std::vector<Key> CylCoefs::makeKeys(Key k)
  {
    std::vector<Key> ret;
    if (coefs.size()==0) return ret;

    // Sanity check
    if (k.size()==1) {
      k[0] = std::max<unsigned>(k[0], 0);
      k[0] = std::min<unsigned>(k[0], Mmax);
    }

    // Two options:
    // 1. return all keys with a fixed m
    // 2. return all keys
    //
    if (k.size()==0) {
      for (unsigned m=0; m<=Mmax; m++)
	for (unsigned n=0; n<Nmax; n++) ret.push_back({m, n});

    }
    else if (k.size()==1) {
      for (unsigned n=0; n<Nmax; n++) ret.push_back({k[0], n});
    }
    // Bad sub key?
    else {
      throw std::runtime_error
	("CylCoefs::makeKeys: the subkey must have rank 1");
    }

    return ret;
  }

  void CylCoefs::readNativeCoefs(const std::string& file, int stride,
				 double tmin, double tmax)
  {
    std::ifstream in(file);
    
    if (not in) {
      throw std::runtime_error("CylCoefs ERROR (runtime) opening file <" + file + ">");
    }
    
    int count = 0;
    while (in) {
      CylStrPtr c = std::make_shared<CylStruct>();
      if (not c->read(in, verbose)) break;
      
      if (count++ % stride) continue;
      if (c->time < tmin or c->time > tmax) continue;

      coefs[roundTime(c->time)] = c;
    }

    if (coefs.size()) {
      Mmax = coefs.begin()->second->mmax;
      Nmax = coefs.begin()->second->nmax;
    }

    times.clear();
    for (auto t : coefs) times.push_back(t.first);

    if (myid==0)
      std::cerr << "---- Coefs::factory: "
		<< "read EXP native and created CylCoefs"
		<< std::endl;
  }
  
  void CylCoefs::WriteH5Params(HighFive::File& file)
  {
    std::string forceID(coefs.begin()->second->id);

    file.createAttribute<int>("mmax", HighFive::DataSpace::From(Mmax)).write(Mmax);
    file.createAttribute<int>("nmax", HighFive::DataSpace::From(Nmax)).write(Nmax);
    file.createAttribute<std::string>("forceID", HighFive::DataSpace::From(forceID)).write(forceID);
  }
  
  unsigned CylCoefs::WriteH5Times(HighFive::Group& snaps, unsigned count)
  {
    for (auto c : coefs) {
      auto C = c.second;
      
      std::ostringstream stim;
      stim << std::setw(8) << std::setfill('0') << std::right << count++;
      HighFive::Group stanza = snaps.createGroup(stim.str());
      
      // Add time attribute
      //
      stanza.createAttribute<double>("Time", HighFive::DataSpace::From(C->time)).write(C->time);
      
      // Add center attribute
      //
      if (C->ctr.size()>0)
	stanza.createAttribute<double>("Center", HighFive::DataSpace::From(C->ctr)).write(C->ctr);

      // Add coefficient data
      //
      HighFive::DataSet dataset = stanza.createDataSet("coefficients", C->coefs);
    }
    
    return count;
  }
  
  
  std::string CylCoefs::getYAML()
  {
    std::string ret;
    if (coefs.size()) {
      ret = coefs.begin()->second->buf;
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
  
  Eigen::MatrixXd& CylCoefs::Power(int min, int max)
  {
    if (coefs.size()) {
      
      int mmax = coefs.begin()->second->mmax;
      int nmax = coefs.begin()->second->nmax;
      
      power.resize(coefs.size(), mmax+1);
      power.setZero();
      
      int T=0;
      for (auto v : coefs) {
	for (int m=0; m<=mmax; m++) {
	  auto rad = v.second->coefs.row(m);
	  for (int n=std::max<int>(0, min); n<std::min<int>(nmax, max); n++) {
	    double val = std::abs(rad(n));
	    power(T, m) += val * val;
	  }
	}
	T++;
      }
    } else {
      power.resize(0, 0);
    }
    
    return power;
  }
  
  std::tuple<Eigen::MatrixXd&, Eigen::MatrixXd&>
  CylCoefs::EvenOddPower(int nodd, int min, int max)
  {
    if (coefs.size()) {
      
      if (nodd<0) {

	YAML::Node node;

	// Create the YAML DB
	//
	try {
	  // Read the YAML from a string
	  //
	  node = YAML::Load(getYAML());
	}
	catch (const std::runtime_error& error) {
	  std::cout << "CylCoefs::EvenOddPower: found a problem while loading "
		    << "the YAML config" << std::endl;
	  throw;
	}

	if (node["ncylodd"]) nodd = node["ncylodd"].as<int>();
	else {
	  std::cout << "CylCoefs::EvenOddPower: ncylodd is not in the YAML "
		    << "config stanza.  Please specify this explicitly as "
		    << "the first argument to EvenOddPower()" << std::endl;
	  throw;
	}
      }

      int mmax = coefs.begin()->second->mmax;
      int nmax = coefs.begin()->second->nmax;
      
      powerE.resize(coefs.size(), mmax+1);
      powerE.setZero();
      
      powerO.resize(coefs.size(), mmax+1);
      powerO.setZero();
      
      int T=0;
      for (auto v : coefs) {
	for (int m=0; m<=mmax; m++) {
	  auto rad = v.second->coefs.row(m);
	  // Even
	  for (int n=std::max<int>(0, min); n<std::min<int>(nmax-nodd, max); n++) {
	    double val = std::abs(rad(n));
	    powerE(T, m) += val * val;
	  }
	  // Odd
	  for (int n=std::max<int>(nmax-nodd, min); n<std::min<int>(nmax, max); n++) {
	    double val = std::abs(rad(n));
	    powerO(T, m) += val * val;
	  }
	}
	T++;
      }
    } else {
      powerE.resize(0, 0);
      powerO.resize(0, 0);
    }
    
    return {powerE, powerO};
  }
  
  TableData::TableData(const std::vector<double>& Times,
		       const std::vector<std::vector<double>>& data,
		       bool verbose) :
    data(data), Coefs("table", verbose)
  {
    times = Times;
    for (int i=0; i<times.size(); i++) {
      TblStrPtr c = std::make_shared<TblStruct>();
      c->cols = data[i].size();
      c->coefs.resize(1, c->cols);
      for (int j=0; j<c->cols; j++) c->coefs(0, j) = data[i][j];
      coefs[roundTime(c->time)] = c;
    }
  }

  TableData::TableData(std::string& file, bool verbose) :
    Coefs("table", verbose)
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
    

  TableData::TableData(HighFive::File& file, int stride,
		       double Tmin, double Tmax, bool verbose) :
    Coefs("table", verbose)
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
      
    for (unsigned n=0; n<count; n+=stride) {
      
      if (times[n] < Tmin or times[n] > Tmax) continue;

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
  
  void TableData::setMatrix(double time, const Eigen::MatrixXcd& mat)
  {
    auto it = coefs.find(roundTime(time));

    if (it == coefs.end()) {
      std::ostringstream str;
      str << "TableData::setMatrix: requested time=" << time << " not found";
      throw std::runtime_error(str.str());
    } else {
      it->second->coefs = mat;
    }
  }

  void TableData::readNativeCoefs(const std::string& file, int stride,
				  double tmin, double tmax)
  {
    std::ifstream in(file);
    
    if (not in) {
      throw std::runtime_error("TableData ERROR (runtime) opening file <" + file + ">");
    }
    
    int count = 0;
    while (in) {
      TblStrPtr c = std::make_shared<TblStruct>();
      if (not c->read(in, verbose)) break;
      
      if (count++ % stride) continue;
      if (c->time < tmin or c->time > tmax) continue;

      coefs[roundTime(c->time)] = c;
    }

    times.clear();
    for (auto t : coefs) times.push_back(t.first);

    if (myid==0)
      std::cerr << "---- Coefs::factory: "
		<< "read ascii and created TableData"
		<< std::endl;
  }
  
  
  void TableData::WriteH5Params(HighFive::File& file)
  {
    int cols = coefs.begin()->second->cols;
    
    std::string geometry = "table";

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
      ret = coefs.begin()->second->buf;
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
  

  Eigen::MatrixXd TableData::getAllCoefs()
  {
    Eigen::MatrixXd ret;

    auto times = Times();

    int cols = coefs.begin()->second->cols;
    int ntim = times.size();

    ret.resize(cols, ntim);
    
    for (int t=0; t<ntim; t++) {
      auto cof = coefs[times[t]];
      for (int c=0; c<cols; c++) {
	ret(c, t) = cof->coefs(0, c).real();
      }
    }

    return ret;
  }


  std::shared_ptr<Coefs> Coefs::factory
  (const std::string& file, int stride, double tmin, double tmax)
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
      std::string geometry;
      HighFive::Attribute geom = h5file.getAttribute("geometry");
      geom.read(geometry);
      
      try {
	if (geometry.compare("sphere")==0) {
	  coefs = std::make_shared<SphCoefs>(h5file, stride, tmin, tmax);
	} else if (geometry.compare("cylinder")==0) {
	  coefs = std::make_shared<CylCoefs>(h5file, stride, tmin, tmax);
	} else if (geometry.compare("table")==0) {
	  coefs = std::make_shared<TableData>(h5file, stride, tmin, tmax);
	} else {
	  throw std::runtime_error("Coefs::factory: unknown H5 coefficient file geometry");
	}
      } catch (HighFive::Exception& err) {
	std::cerr << "**** Error reading HDF5 file ****" << std::endl;
	std::cerr << err.what() << std::endl;
	exit(-1);
      }
      
      return coefs;
      
    } catch (HighFive::Exception& err) {
      if (myid==0)
	std::cerr << "---- Coefs::factory: "
		  << "error opening as HDF5, trying EXP native and ascii table"
		  << std::endl;
    }
    
    // Check if file exists
    //
    if (not std::filesystem::exists(file)) {
      throw std::runtime_error("Coefs::factory: file <" + file
			       + "> does not exist");
    }

    // Open file and read magic number
    //
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
    
    coefs->readNativeCoefs(file, stride, tmin, tmax);
    
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
      throw std::runtime_error("Coefs::makecoefs: cannot deduce coefficient file type");
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
    //
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
      HighFive::File file(prefix,
			  HighFive::File::ReadWrite |
			  HighFive::File::Create);
      
      // We write the coefficient file geometry
      //
      file.createAttribute<std::string>("geometry", HighFive::DataSpace::From(geometry)).write(geometry);
      
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
      // Open an hdf5 file
      //
      HighFive::File file(prefix, HighFive::File::ReadWrite);
      
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
    auto p = std::dynamic_pointer_cast<CylStruct>(coef);
    Mmax = p->mmax;
    Nmax = p->nmax;
    coefs[roundTime(coef->time)] = p;
  }

  void TableData::add(CoefStrPtr coef)
  {
    coefs[roundTime(coef->time)] = std::dynamic_pointer_cast<TblStruct>(coef);
  }


}
// END namespace CoefClasses
