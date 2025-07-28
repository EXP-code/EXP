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

#include <highfive/highfive.hpp>
#include <highfive/eigen.hpp>

#include <Coefficients.H>

namespace CoefClasses
{
  
  void Coefs::copyfields(std::shared_ptr<Coefs> p)
  {
    // These variables will copy data, not pointers
    p->arr      = arr;
    p->power    = power;
    p->geometry = geometry;
    p->name     = name;
    p->verbose  = verbose;
    p->times    = times;
  }

  std::tuple<Eigen::VectorXcd&, bool> Coefs::interpolate(double time)
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

    int siz = cA->store.size();
    arr.resize(siz);

    for (int c=0; c<siz; c++) {
      arr(c) = A*cA->store(c) + B*cB->store(c);
    }

    return {arr, onGrid};
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
    
    // Look for Coef output version to toggle backward compatibility
    // with legacy storage order
    //
    bool H5back = true;
    if (file.hasAttribute("CoefficientOutputVersion")) H5back = false;

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

      auto in = stanza.getDataSet("coefficients").read<Eigen::MatrixXcd>();

      // If we have a legacy set of coefficients, re-order the
      // coefficients to match the new HighFive/Eigen ordering
      //
      if (H5back) {

	auto in2 = stanza.getDataSet("coefficients").read<Eigen::MatrixXcd>();
	in2.transposeInPlace();
    
	for (size_t c=0, n=0; c<in.cols(); c++) {
	  for (size_t r=0; r<in.rows(); r++) {
	    in(r, c) = in2.data()[n++];
	  }
	}
      }
      
      // Pack the data into the coefficient variable
      //
      auto coef = std::make_shared<SphStruct>();
      
      if (ctr.size()) coef->ctr = ctr;

      coef->lmax  = Lmax;
      coef->nmax  = Nmax;
      coef->time  = Time;
      coef->scale = scale;
      coef->geom  = geometry;
      coef->id    = forceID;

      coef->allocate();
      *coef->coefs = in;
      
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

  std::shared_ptr<Coefs> SlabCoefs::deepcopy()
  {
    auto ret = std::make_shared<SlabCoefs>();

    // Copy the base-class fields
    copyfields(ret);

    // Copy the local structures from the map to the struct pointers
    // by copyfing fields, not the pointer
    for (auto v : coefs)
      ret->coefs[v.first] =
	std::dynamic_pointer_cast<SlabStruct>(v.second->deepcopy());

    ret->NmaxX  = NmaxX;
    ret->NmaxY  = NmaxY;
    ret->NmaxZ  = NmaxZ;

    return ret;
  }

  std::shared_ptr<Coefs> CubeCoefs::deepcopy()
  {
    auto ret = std::make_shared<CubeCoefs>();

    // Copy the base-class fields
    copyfields(ret);

    // Copy the local structures from the map to the struct pointers
    // by copyfing fields, not the pointer
    for (auto v : coefs)
      ret->coefs[v.first] =
	std::dynamic_pointer_cast<CubeStruct>(v.second->deepcopy());

    ret->NmaxX  = NmaxX;
    ret->NmaxY  = NmaxY;
    ret->NmaxZ  = NmaxZ;

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


  std::shared_ptr<Coefs> TrajectoryData::deepcopy()
  {
    auto ret = std::make_shared<TrajectoryData>();

    // Copy the base-class fields
    copyfields(ret);

    // Copy the local structures from the map to the struct pointers
    // by copyfing fields, not the pointer
    for (auto v : coefs)
      ret->coefs[v.first] =
	std::dynamic_pointer_cast<TrajStruct>(v.second->deepcopy());

    ret->data  = data;
    ret->times = times;

    return ret;
  }


  SphFldCoefs::SphFldCoefs(HighFive::File& file, int stride,
			   double Tmin, double Tmax, bool verbose) :
    Coefs("sphere", verbose)
  {
    std::string config, geometry, fieldID;
    unsigned count;
    double scale;
    
    file.getAttribute("name"    ).read(name    );
    file.getAttribute("nfld"    ).read(Nfld    );
    file.getAttribute("lmax"    ).read(Lmax    );
    file.getAttribute("nmax"    ).read(Nmax    );
    file.getAttribute("scale"   ).read(scale   );
    file.getAttribute("config"  ).read(config  );
    file.getDataSet  ("count"   ).read(count   );
    file.getAttribute("geometry").read(geometry);
    file.getAttribute("fieldID" ).read(fieldID );
    
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

      std::array<long int, 3> shape;
      stanza.getAttribute("shape").read(shape);

      auto in = stanza.getDataSet("coefficients").read<Eigen::VectorXcd>();
      
      // Pack the data into the coefficient variable
      //
      auto coef = std::make_shared<SphFldStruct>();
      
      if (ctr.size()) coef->ctr = ctr;

      coef->nfld  = Nfld;
      coef->lmax  = Lmax;
      coef->nmax  = Nmax;
      coef->time  = Time;
      coef->scale = scale;
      coef->geom  = geometry;
      coef->id    = fieldID;

      coef->allocate();
      coef->store = in;
      
      coefs[roundTime(Time)] = coef;
    }

    times.clear();
    for (auto t : coefs) times.push_back(t.first);
  }
  
  std::shared_ptr<Coefs> SphFldCoefs::deepcopy()
  {
    auto ret = std::make_shared<SphFldCoefs>();

    // Copy the base-class fields
    copyfields(ret);

    // Copy the local structures from the map to the struct pointers
    // by copyfing fields, not the pointer
    for (auto v : coefs)
      ret->coefs[v.first] =
	std::dynamic_pointer_cast<SphFldStruct>(v.second->deepcopy());

    ret->Nfld = Nfld;
    ret->Lmax = Lmax;
    ret->Nmax = Nmax;

    return ret;
  }

  
  CylFldCoefs::CylFldCoefs(HighFive::File& file, int stride,
			       double Tmin, double Tmax, bool verbose) :
    Coefs("cylinder", verbose)
  {
    std::string config, geometry, fieldID;
    unsigned count;
    double scale;
    
    file.getAttribute("name"    ).read(name    );
    file.getAttribute("nfld"    ).read(Nfld    );
    file.getAttribute("mmax"    ).read(Mmax    );
    file.getAttribute("nmax"    ).read(Nmax    );
    file.getAttribute("scale"   ).read(scale   );
    file.getAttribute("config"  ).read(config  );
    file.getDataSet  ("count"   ).read(count   );
    file.getAttribute("geometry").read(geometry);
    file.getAttribute("fieldID" ).read(fieldID );
    
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

      std::array<long int, 3> shape;
      stanza.getAttribute("shape").read(shape);

      auto in = stanza.getDataSet("coefficients").read<Eigen::VectorXcd>();
      
      // Pack the data into the coefficient variable
      //
      auto coef = std::make_shared<CylFldStruct>();
      
      if (ctr.size()) coef->ctr = ctr;

      coef->nfld  = Nfld;
      coef->mmax  = Mmax;
      coef->nmax  = Nmax;
      coef->time  = Time;
      coef->scale = scale;
      coef->geom  = geometry;
      coef->id    = fieldID;

      coef->allocate();
      coef->store = in;
      
      coefs[roundTime(Time)] = coef;
    }

    times.clear();
    for (auto t : coefs) times.push_back(t.first);
  }
  
  std::shared_ptr<Coefs> CylFldCoefs::deepcopy()
  {
    auto ret = std::make_shared<CylFldCoefs>();

    // Copy the base-class fields
    copyfields(ret);

    // Copy the local structures from the map to the struct pointers
    // by copyfing fields, not the pointer
    for (auto v : coefs)
      ret->coefs[v.first] =
	std::dynamic_pointer_cast<CylFldStruct>(v.second->deepcopy());

    ret->Nfld = Nfld;
    ret->Mmax = Mmax;
    ret->Nmax = Nmax;

    return ret;
  }


  Eigen::VectorXcd& SphCoefs::getData(double time)
  {
    auto it = coefs.find(roundTime(time));

    if (it == coefs.end()) {
      arr.resize(0);
    } else {
      arr = it->second->store;
    }
    
    return arr;
  }
  
  Eigen::MatrixXcd& SphCoefs::getMatrix(double time)
  {
    auto it = coefs.find(roundTime(time));

    if (it == coefs.end()) {
      arr.resize(0);
    } else {
      arr = it->second->store;
      int ldim = (Lmax+1)*(Lmax+2)/2;
      mat = Eigen::Map<Eigen::MatrixXcd>(arr.data(), ldim, Nmax); 
    }
    
    return mat;
  }
  
  void SphCoefs::setData(double time, Eigen::VectorXcd& dat)
  {
    auto it = coefs.find(roundTime(time));

    if (it == coefs.end()) {
      std::ostringstream str;
      str << "SphCoefs::setData: requested time=" << time << " not found";
      throw std::runtime_error(str.str());
    } else {
      it->second->store = dat;
      it->second->coefs = std::make_shared<CoefClasses::SphStruct::coefType>
	(it->second->store.data(), (Lmax+1)*(Lmax+2)/2, Nmax);
    }
  }
  
  void SphCoefs::setMatrix(double time, Eigen::MatrixXcd& dat)
  {
    auto it = coefs.find(roundTime(time));

    if (it == coefs.end()) {
      std::ostringstream str;
      str << "SphCoefs::setMatrix: requested time=" << time << " not found";
      throw std::runtime_error(str.str());
    } else {
      it->second->allocate();
      *it->second->coefs = dat;
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
      auto & cof = *(coefs[roundTime(times[t])]->coefs);
      for (int l=0; l<(Lmax+2)*(Lmax+1)/2; l++) {
	for (int n=0; n<Nmax; n++) {
	  ret(l, n, t) = cof(l, n);
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
      Eigen::MatrixXcd out(*C->coefs);
      HighFive::DataSet dataset = stanza.createDataSet("coefficients", out);
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
      
      auto & cof = *(c.second->coefs);

      for (int ll=lmin; ll<=std::min<int>(lmax, Lmax); ll++) {
	for (int mm=0; mm<=ll; mm++) {
	  std::cout << std::setw(18) << c.first << std::setw(5) << ll << std::setw(5) << mm << std::setw(5);
	  for (int nn=std::max<int>(nmin, 0); nn<std::min<int>(nmax, Nmax); nn++) 
	    std::cout << std::setw(18) << cof(I, nn);
	  std::cout << std::endl;
	  
	} // M loop
	
      } // L loop
      
    } // T loop
  }
  
  
  bool SphCoefs::CompareStanzas(CoefsPtr check)
  {
    bool ret = true;
    
    auto other = std::dynamic_pointer_cast<SphCoefs>(check);
    
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
	auto & cv = *(v.second->coefs);
	auto & ci = *(it->second->coefs);
	for (int i=0; i<cv.rows(); i++) {
	  for (int j=0; j<cv.cols(); j++) {
	    if (cv(i, j) != ci(i, j)) {
	      std::cout << "Coefficient (" << i << ", " << j << ")  "
			<< cv(i, j) << " != "
			<< ci(i, j) << std::endl;
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
	    auto rad = v.second->coefs->row(L);
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
    if (not p) throw std::runtime_error("SphCoefs::add: Null coefficient structure, nothing added!");

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
    
    // Look for Coef output version to toggle backward compatibility
    // with legacy storage order
    //
    bool H5back = true;
    if (file.hasAttribute("CoefficientOutputVersion")) H5back = false;

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

      auto in = stanza.getDataSet("coefficients").read<Eigen::MatrixXcd>();

      // If we have a legacy set of coefficients, re-order the
      // coefficients to match the new HighFive/Eigen ordering
      //
      if (H5back) {

	auto in2 = stanza.getDataSet("coefficients").read<Eigen::MatrixXcd>();
	in2.transposeInPlace();
    
	for (size_t c=0, n=0; c<in.cols(); c++) {
	  for (size_t r=0; r<in.rows(); r++) {
	    in(r, c) = in2.data()[n++];
	  }
	}
      }
      
      // Work around for previous unitiaized data bug; enforces real data
      //
      for (int n=0; n<Nmax; n++) in(0, n) = std::real(in(0, n));

      // Pack the data into the coefficient variable
      //
      auto coef = std::make_shared<CylStruct>();
      
      if (ctr.size()) coef->ctr = ctr;

      coef->assign(in, Mmax, Nmax);
      coef->time = Time;
      
      coefs[roundTime(Time)] = coef;
    }

    times.clear();
    for (auto t : coefs) times.push_back(t.first);
  }
  
  Eigen::VectorXcd& CylCoefs::getData(double time)
  {
    auto it = coefs.find(roundTime(time));

    if (it == coefs.end()) {
      arr.resize(0);
    } else {
      arr = it->second->store;
    }
    
    return arr;
  }
  
  Eigen::MatrixXcd& CylCoefs::getMatrix(double time)
  {
    auto it = coefs.find(roundTime(time));

    if (it == coefs.end()) {
      arr.resize(0);
    } else {
      arr = it->second->store;
      mat = Eigen::Map<Eigen::MatrixXcd>(arr.data(), Mmax+1, Nmax); 
    }
    
    return mat;
  }
  
  void CylCoefs::setData(double time, Eigen::VectorXcd& dat)
  {
    auto it = coefs.find(roundTime(time));

    if (it == coefs.end()) {
      std::ostringstream str;
      str << "CylCoefs::setData: requested time=" << time << " not found";
      throw std::runtime_error(str.str());
    } else {
      it->second->store = dat;
      it->second->coefs = std::make_shared<CoefClasses::CylStruct::coefType>
	(it->second->store.data(), Mmax+1, Nmax);
    }
  }

  void CylCoefs::setMatrix(double time, Eigen::MatrixXcd& dat)
  {
    auto it = coefs.find(roundTime(time));

    if (it == coefs.end()) {
      std::ostringstream str;
      str << "CylCoefs::setMatrix: requested time=" << time << " not found";
      throw std::runtime_error(str.str());
    } else {
      it->second->allocate();
      *it->second->coefs = dat;
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
      auto & cof = *coefs[roundTime(times[t])]->coefs;
      for (int m=0; m<Mmax+1; m++) {
	for (int n=0; n<Nmax; n++) {
	  ret(m, n, t) = cof(m, n);
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
      Eigen::MatrixXcd out(*C->coefs);
      HighFive::DataSet dataset = stanza.createDataSet("coefficients", out);
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

      // The coefficient matrix
      auto & cof = *c.second->coefs;

      // Index loop
      for (int mm=mmin; mm<=std::min<int>(mmax, c.second->mmax); mm++) {
	std::cout << std::setw(18) << c.first << std::setw(5) << mm;
	for (int nn=std::max<int>(nmin, 0); nn<std::min<int>(nmax, c.second->nmax); nn++) {
	  if (angle)
	    std::cout << std::setw(18) << abs(cof(mm, nn))
		      << std::setw(18) << arg(cof(mm, nn));
	  else
	    std::cout << std::setw(18) << cof(mm, nn).real()
		      << std::setw(18) << cof(mm, nn).imag();
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
	  auto rad = (*v.second->coefs).row(m);
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
	  throw std::runtime_error
	    (
	     "CylCoefs::EvenOddPower: found a problem while loading the "
	     "YAML config"
	     );
	}

	if (node["ncylodd"]) nodd = node["ncylodd"].as<int>();
	else {
	  throw std::runtime_error
	    (
	    "CylCoefs::EvenOddPower: ncylodd is not in the YAML config "
	    "stanza.  Please specify this explicitly as the first argument "
	    "to EvenOddPower()"
	     );
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
	  auto rad = (*v.second->coefs).row(m);
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
  
  SlabCoefs::SlabCoefs(HighFive::File& file, int stride,
		       double Tmin, double Tmax, bool verbose) :
    Coefs("cube", verbose)
  {
    unsigned count;
    std::string config;
    
    file.getAttribute("name"   ).read(name  );
    file.getAttribute("nmaxx"  ).read(NmaxX );
    file.getAttribute("nmaxy"  ).read(NmaxY );
    file.getAttribute("nmaxz"  ).read(NmaxZ );
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

      auto in = stanza.getDataSet("coefficients").read<Eigen::VectorXcd>();
      
      Eigen::TensorMap<Eigen3d> dat(in.data(), 2*NmaxX+1, 2*NmaxY+1, NmaxZ);

      // Pack the data into the coefficient variable
      //
      auto coef = std::make_shared<SlabStruct>();
      
      coef->assign(dat);
      coef->time = Time;
      
      coefs[roundTime(Time)] = coef;
    }

    times.clear();
    for (auto t : coefs) times.push_back(t.first);
  }
  
  Eigen::VectorXcd& SlabCoefs::getData(double time)
  {
    auto it = coefs.find(roundTime(time));

    if (it == coefs.end()) {
      arr.resize(0);
    } else {
      arr = it->second->store;
    }
    
    return arr;
  }

  void SlabCoefs::setData(double time, Eigen::VectorXcd& dat)
  {
    auto it = coefs.find(roundTime(time));

    if (it == coefs.end()) {
      std::ostringstream str;
      str << "CylCoefs::setData: requested time=" << time << " not found";
      throw std::runtime_error(str.str());
    } else {
      it->second->store = dat;
      it->second->coefs = std::make_shared<CoefClasses::SlabStruct::coefType>
	(it->second->store.data(), 2*NmaxX+1, 2*NmaxY+1, NmaxZ);
    }
  }

  void SlabCoefs::setTensor(double time, const Eigen3d& dat)
  {
    auto it = coefs.find(roundTime(time));

    if (it == coefs.end()) {
      std::ostringstream str;
      str << "SlabCoefs::setTensor: requested time=" << time << " not found";
      throw std::runtime_error(str.str());
    } else {
      it->second->allocate();	// Assign storage for the flattened tensor
      *it->second->coefs = dat;	// Populate using the tensor map
    }
  }

  SlabCoefs::Eigen3d& SlabCoefs::getTensor(double time)
  {
    auto it = coefs.find(roundTime(time));

    if (it == coefs.end()) {
      arr.resize(0);
    } else {
      arr = it->second->store;
      dat = Eigen::TensorMap<Eigen3d>(arr.data(), 2*NmaxX+1, 2*NmaxY+1, NmaxZ);
    }
    
    return dat;
  }

  Eigen::Tensor<std::complex<double>, 4> SlabCoefs::getAllCoefs()
  {
    Eigen::Tensor<std::complex<double>, 4> ret;

    auto times = Times();

    int ntim = times.size();

    // Resize the tensor
    ret.resize(2*NmaxX+1, 2*NmaxY+1, NmaxZ, ntim);
    
    for (int t=0; t<ntim; t++) {
      auto cof = coefs[roundTime(times[t])];
      for (int ix=0; ix<=2*NmaxX; ix++) {
	for (int iy=0; iy<=2*NmaxY; iy++) {
	  for (int iz=0; iz<NmaxZ; iz++) {
	    ret(ix, iy, iz, t) = (*cof->coefs)(ix, iy, iz);
	  }
	}
      }
    }

    return ret;
  }

  std::vector<Key> SlabCoefs::makeKeys()
  {
    std::vector<Key> ret;
    if (coefs.size()==0) return ret;

    for (unsigned ix=0; ix<=2*NmaxX; ix++) {
      for (unsigned iy=0; iy<=2*NmaxY; iy++) {
	for (unsigned iz=0; iz<NmaxZ; iz++) {
	  ret.push_back({ix, iy, iz});
	}
      }
    }
    
    return ret;
  }

  void SlabCoefs::WriteH5Params(HighFive::File& file)
  {
    std::string forceID(coefs.begin()->second->id);

    file.createAttribute<int>("nmaxx", HighFive::DataSpace::From(NmaxX)).write(NmaxX);
    file.createAttribute<int>("nmaxy", HighFive::DataSpace::From(NmaxY)).write(NmaxY);
    file.createAttribute<int>("nmaxz", HighFive::DataSpace::From(NmaxZ)).write(NmaxZ);
    file.createAttribute<std::string>("forceID", HighFive::DataSpace::From(forceID)).write(forceID);
  }
  
  unsigned SlabCoefs::WriteH5Times(HighFive::Group& snaps, unsigned count)
  {
    for (auto c : coefs) {
      auto C = c.second;
      
      std::ostringstream stim;
      stim << std::setw(8) << std::setfill('0') << std::right << count++;
      HighFive::Group stanza = snaps.createGroup(stim.str());
      
      // Add time attribute
      //
      stanza.createAttribute<double>("Time", HighFive::DataSpace::From(C->time)).write(C->time);
      
      // Add coefficient data
      //
      HighFive::DataSet dataset = stanza.createDataSet("coefficients", C->store);
    }
    
    return count;
  }
  
  std::string SlabCoefs::getYAML()
  {
    std::string ret;
    if (coefs.size()) {
      ret = coefs.begin()->second->buf;
    }
    return ret;
  }
  
  void SlabCoefs::dump(int nmaxx, int nmaxy, int nmaxz)
  {
    for (auto c : coefs) {
      for (int ix=0; ix<std::min<int>(nmaxx, c.second->nmaxx); ix++) {
	std::cout << std::setw(18) << c.first << std::setw(5) << ix;
	for (int iy=0; iy<std::min<int>(nmaxy, c.second->nmaxy); iy++) {
	  for (int iz=0; iz<std::min<int>(nmaxz, c.second->nmaxz); iz++) {
	    std::cout << std::setw(18) << c.first
		      << std::setw(5) << ix
		      << std::setw(5) << iy
		      << std::setw(5) << iz
		      << std::setw(18) << (*c.second->coefs)(ix, iy, iz).real()
		      << std::setw(18) << (*c.second->coefs)(ix, iy, iz).imag()
		      << std::endl;
	  }
	  // Z loop
	}
	// Y loop
      }
      // X loop
    }
    // T loop
  }
  
  Eigen::MatrixXd& SlabCoefs::Power(char d, int min, int max)
  {
    if (coefs.size()) {
      
      int nmaxX = coefs.begin()->second->nmaxx;
      int nmaxY = coefs.begin()->second->nmaxy;
      int nmaxZ = coefs.begin()->second->nmaxz;
      
      // Internal sanity check
      assert(nmaxX == NmaxX && "nmaxX <==> NmaxX mismatch");
      assert(nmaxY == NmaxY && "nmaxY <==> NmaxY mismatch");
      assert(nmaxZ == NmaxZ && "nmaxZ <==> NmaxZ mismatch");

      int dim = 0;
      if (d == 'x')
	dim  = 2*nmaxX + 1;
      else if (d == 'y')
	dim  = 2*nmaxY + 1;
      else
	dim  = nmaxZ;

      power.resize(coefs.size(), dim);
      power.setZero();
      
      int T=0;
      for (auto v : coefs) {
	if (d=='x') {
	  for (int ix=0; ix<=2*NmaxX; ix++) {
	    double val(0.0);
	    for (int iy=0; iy<=2*NmaxY; iy++) {
	      if (abs(iy - nmaxY) < min) continue;
	      for (int iz=0; iz<NmaxZ; iz++) {
		if (abs(iz - nmaxZ) < min) continue;
		double val = std::abs((*v.second->coefs)(ix, iy, iz));
		power(T, ix) += val * val;
	      }
	    }
	  }
	} else if (d=='y') {
	  for (int iy=0; iy<=2*NmaxY; iy++) {
	    double val(0.0);
	    for (int ix=0; ix<=2*NmaxX; ix++) {
	      if (abs(ix - nmaxX) < min) continue;
	      for (int iz=0; iz<NmaxZ; iz++) {
		if (abs(iz - nmaxZ) < min) continue;
		double val = std::abs((*v.second->coefs)(ix, iy, iz));
		power(T, iy) += val * val;
	      }
	    }
	  }
	} else {
	  for (int iz=0; iz<NmaxZ; iz++) {
	    double val(0.0);
	    for (int ix=0; ix<=2*NmaxX; ix++) {
	      if (abs(ix - nmaxX) < min) continue;
	      for (int iy=0; iy<=2*NmaxY; iy++) {
		if (abs(iy - nmaxY) < min) continue;
		double val = std::abs((*v.second->coefs)(ix, iy, iz));
		power(T, iz) += val * val;
	      }
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
  
  void SlabCoefs::readNativeCoefs(const std::string& file,
				  int stride, double tmin, double tmax)
  {
    std::runtime_error("SlabCoefs: no native coefficients files");
  }

  bool SlabCoefs::CompareStanzas(CoefsPtr check)
  {
    bool ret = true;
    
    auto other = std::dynamic_pointer_cast<SlabCoefs>(check);
    
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
	if (v.second->nmaxx != it->second->nmaxx) ret = false;
	if (v.second->nmaxy != it->second->nmaxy) ret = false;
	if (v.second->nmaxz != it->second->nmaxz) ret = false;
	if (v.second->time  != it->second->time ) ret = false;
      }
    }
    
    if (ret) {
      std::cout << "Parameters are the same, now checking coefficients"
		<< std::endl;
      for (auto v : coefs) {
	auto it = other->coefs.find(v.first);
	auto & cv = *(v.second->coefs);
	auto & ci = *(it->second->coefs);
	auto  dim = cv.dimensions(); // This is an Eigen::Tensor map
	for (int i=0; i<dim[0]; i++) {
	  for (int j=0; j<dim[1]; j++) {
	    for (int k=0; k<dim[2]; k++) {
	      if (cv(i, j, k) != ci(i, j, k)) {
		std::cout << "Coefficient (" << i << ", " << j << ", " << k << ")  "
			  << cv(i, j, k) << " != "
			  << ci(i, j, k) << std::endl;
		ret = false;
	      }
	    }
	  }
	}
      }
    }
    
    return ret;
  }
  
  CubeCoefs::CubeCoefs(HighFive::File& file, int stride,
		     double Tmin, double Tmax, bool verbose) :
    Coefs("cube", verbose)
  {
    unsigned count;
    std::string config;
    
    file.getAttribute("name"   ).read(name  );
    file.getAttribute("nmaxx"  ).read(NmaxX );
    file.getAttribute("nmaxy"  ).read(NmaxY );
    file.getAttribute("nmaxz"  ).read(NmaxZ );
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

      auto in = stanza.getDataSet("coefficients").read<Eigen::VectorXcd>();
      
      Eigen::TensorMap<Eigen3d> dat(in.data(), 2*NmaxX+1, 2*NmaxY+1, 2*NmaxZ+1);

      // Pack the data into the coefficient variable
      //
      auto coef = std::make_shared<CubeStruct>();
      
      coef->assign(dat);
      coef->time = Time;
      
      coefs[roundTime(Time)] = coef;
    }

    times.clear();
    for (auto t : coefs) times.push_back(t.first);
  }
  
  Eigen::VectorXcd& CubeCoefs::getData(double time)
  {
    auto it = coefs.find(roundTime(time));

    if (it == coefs.end()) {
      arr.resize(0);
    } else {
      arr = it->second->store;
    }
    
    return arr;
  }

  void CubeCoefs::setData(double time, Eigen::VectorXcd& dat)
  {
    auto it = coefs.find(roundTime(time));

    if (it == coefs.end()) {
      std::ostringstream str;
      str << "CubeCoefs::setData: requested time=" << time << " not found";
      throw std::runtime_error(str.str());
    } else {
      it->second->store = dat;
      it->second->coefs = std::make_shared<CoefClasses::CubeStruct::coefType>
	(it->second->store.data(), 2*NmaxX+1, 2*NmaxY+1, 2*NmaxZ+1);
    }
  }

  void CubeCoefs::setTensor(double time, const Eigen3d& dat)
  {
    auto it = coefs.find(roundTime(time));

    if (it == coefs.end()) {
      std::ostringstream str;
      str << "CubeCoefs::setTensor: requested time=" << time << " not found";
      throw std::runtime_error(str.str());
    } else {
      it->second->allocate();	// Assign storage for the flattened tensor
      *it->second->coefs = dat;	// Populate using the tensor map
    }
  }

  CubeCoefs::Eigen3d& CubeCoefs::getTensor(double time)
  {
    auto it = coefs.find(roundTime(time));

    if (it == coefs.end()) {
      arr.resize(0);
    } else {
      arr = it->second->store;
      dat = Eigen::TensorMap<Eigen3d>(arr.data(), 2*NmaxX+1, 2*NmaxY+1, 2*NmaxZ+1);
    }
    
    return dat;
  }

  Eigen::Tensor<std::complex<double>, 4> CubeCoefs::getAllCoefs()
  {
    Eigen::Tensor<std::complex<double>, 4> ret;

    auto times = Times();

    int ntim = times.size();

    // Resize the tensor
    ret.resize(2*NmaxX+1, 2*NmaxY+1, 2*NmaxZ+1, ntim);
    
    for (int t=0; t<ntim; t++) {
      auto cof = coefs[roundTime(times[t])];
      for (int ix=0; ix<=2*NmaxX; ix++) {
	for (int iy=0; iy<=2*NmaxY; iy++) {
	  for (int iz=0; iz<=2*NmaxZ; iz++) {
	    ret(ix, iy, iz, t) = (*cof->coefs)(ix, iy, iz);
	  }
	}
      }
    }

    return ret;
  }

  std::vector<Key> CubeCoefs::makeKeys()
  {
    std::vector<Key> ret;
    if (coefs.size()==0) return ret;

    for (unsigned ix=0; ix<=2*NmaxX; ix++) {
      for (unsigned iy=0; iy<=2*NmaxY; iy++) {
	for (unsigned iz=0; iz<=2*NmaxZ; iz++) {
	  ret.push_back({ix, iy, iz});
	}
      }
    }
    
    return ret;
  }

  void CubeCoefs::WriteH5Params(HighFive::File& file)
  {
    std::string forceID(coefs.begin()->second->id);

    file.createAttribute<int>("nmaxx", HighFive::DataSpace::From(NmaxX)).write(NmaxX);
    file.createAttribute<int>("nmaxy", HighFive::DataSpace::From(NmaxY)).write(NmaxY);
    file.createAttribute<int>("nmaxz", HighFive::DataSpace::From(NmaxZ)).write(NmaxZ);
    file.createAttribute<std::string>("forceID", HighFive::DataSpace::From(forceID)).write(forceID);
  }
  
  unsigned CubeCoefs::WriteH5Times(HighFive::Group& snaps, unsigned count)
  {
    for (auto c : coefs) {
      auto C = c.second;
      
      std::ostringstream stim;
      stim << std::setw(8) << std::setfill('0') << std::right << count++;
      HighFive::Group stanza = snaps.createGroup(stim.str());
      
      // Add time attribute
      //
      stanza.createAttribute<double>("Time", HighFive::DataSpace::From(C->time)).write(C->time);
      
      // Add coefficient data
      //
      HighFive::DataSet dataset = stanza.createDataSet("coefficients", C->store);
    }
    
    return count;
  }
  
  std::string CubeCoefs::getYAML()
  {
    std::string ret;
    if (coefs.size()) {
      ret = coefs.begin()->second->buf;
    }
    return ret;
  }
  
  void CubeCoefs::dump(int nmaxx, int nmaxy, int nmaxz)
  {
    for (auto c : coefs) {
      for (int ix=0; ix<std::min<int>(nmaxx, c.second->nmaxx); ix++) {
	std::cout << std::setw(18) << c.first << std::setw(5) << ix;
	for (int iy=0; iy<std::min<int>(nmaxy, c.second->nmaxy); iy++) {
	  for (int iz=0; iz<std::min<int>(nmaxz, c.second->nmaxz); iz++) {
	    std::cout << std::setw(18) << c.first
		      << std::setw(5) << ix
		      << std::setw(5) << iy
		      << std::setw(5) << iz
		      << std::setw(18) << (*c.second->coefs)(ix, iy, iz).real()
		      << std::setw(18) << (*c.second->coefs)(ix, iy, iz).imag()
		      << std::endl;
	  }
	  // Z loop
	}
	// Y loop
      }
      // X loop
    }
    // T loop
  }
  
  Eigen::MatrixXd& CubeCoefs::Power(char d, int min, int max)
  {
    if (coefs.size()) {
      
      int nmaxX = coefs.begin()->second->nmaxx;
      int nmaxY = coefs.begin()->second->nmaxy;
      int nmaxZ = coefs.begin()->second->nmaxz;
      
      int dim = 0;
      if (d == 'x')
	dim  = 2*nmaxX + 1;
      else if (d == 'y')
	dim  = 2*nmaxY + 1;
      else
	dim  = 2*nmaxZ + 1;

      power.resize(coefs.size(), dim);
      power.setZero();
      
      int T=0;
      for (auto v : coefs) {
	if (d=='x') {
	  for (int ix=0; ix<=2*NmaxX; ix++) {
	    double val(0.0);
	    for (int iy=0; iy<=2*NmaxY; iy++) {
	      if (abs(iy - nmaxY) < min) continue;
	      for (int iz=0; iz<=2*NmaxZ; iz++) {
		if (abs(iz - nmaxZ) < min) continue;
		double val = std::abs((*v.second->coefs)(ix, iy, iz));
		power(T, ix) += val * val;
	      }
	    }
	  }
	} else if (d=='y') {
	  for (int iy=0; iy<=2*NmaxY; iy++) {
	    double val(0.0);
	    for (int ix=0; iy<=2*NmaxX; ix++) {
	      if (abs(ix - nmaxX) < min) continue;
	      for (int iz=0; iz<=2*NmaxZ; iz++) {
		if (abs(iz - nmaxZ) < min) continue;
		double val = std::abs((*v.second->coefs)(ix, iy, iz));
		power(T, iy) += val * val;
	      }
	    }
	  }
	} else {
	  for (int iz=0; iz<=2*NmaxZ; iz++) {
	    double val(0.0);
	    for (int ix=0; ix<=2*NmaxX; ix++) {
	      if (abs(ix - nmaxX) < min) continue;
	      for (int iy=0; iy<=2*NmaxY; iy++) {
		if (abs(iy - nmaxY) < min) continue;
		double val = std::abs((*v.second->coefs)(ix, iy, iz));
		power(T, iz) += val * val;
	      }
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
  
  void CubeCoefs::readNativeCoefs(const std::string& file,
				  int stride, double tmin, double tmax)
  {
    std::runtime_error("CubeCoefs: no native coefficients files");
  }

  bool CubeCoefs::CompareStanzas(CoefsPtr check)
  {
    bool ret = true;
    
    auto other = std::dynamic_pointer_cast<CubeCoefs>(check);
    
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
	if (v.second->nmaxx != it->second->nmaxx) ret = false;
	if (v.second->nmaxy != it->second->nmaxy) ret = false;
	if (v.second->nmaxz != it->second->nmaxz) ret = false;
	if (v.second->time  != it->second->time ) ret = false;
      }
    }
    
    if (ret) {
      std::cout << "Parameters are the same, now checking coefficients"
		<< std::endl;
      for (auto v : coefs) {
	auto it = other->coefs.find(v.first);
	auto & cv = *(v.second->coefs);
	auto & ci = *(it->second->coefs);
	auto  dim = cv.dimensions(); // This is an Eigen::Tensor map
	for (int i=0; i<dim[0]; i++) {
	  for (int j=0; j<dim[1]; j++) {
	    for (int k=0; k<dim[2]; k++) {
	      if (cv(i, j, k) != ci(i, j, k)) {
		std::cout << "Coefficient (" << i << ", " << j << ", " << k << ")  "
			  << cv(i, j, k) << " != "
			  << ci(i, j, k) << std::endl;
		ret = false;
	      }
	    }
	  }
	}
      }
    }
    
    return ret;
  }
  

  TrajectoryData::TrajectoryData(const std::vector<double>& Times,
				 const std::vector<Eigen::MatrixXd>& data,
				 bool verbose) :
    data(data), Coefs("table", verbose)
  {
    times = Times;
    int traj = data.size();
    int ntim = data[0].rows();
    int rank = data[0].cols();

    if (ntim != times.size()) {
      std::ostringstream msg;
      msg << "TrajectoryData ERROR (runtime) ntim [" << ntim
	  << "] != times.size() [" << times.size() << "]";
      throw std::runtime_error(msg.str());
    }

    for (int i=0; i<times.size(); i++) {
      TrajStrPtr c = std::make_shared<TrajStruct>();
      c->time = times[i];
      c->traj = traj;
      c->rank = rank;
      c->store.resize(c->traj*c->rank);
      c->coefs = std::make_shared<TrajStruct::coefType>(c->store.data(), c->traj, c->rank);
      for (int m=0; m<traj; m++) {
	for (int n=0; n<rank; n++) {
	  (*c->coefs)(m, n) = data[m](i, n);
	}
      }
      coefs[roundTime(c->time)] = c;
    }
  }

  TrajectoryData::TrajectoryData(std::string& file, bool verbose) :
    Coefs("trajectory", verbose)
  {
    std::ifstream in(file);
    
    if (not in) {
      throw std::runtime_error("TrajectoryData ERROR (runtime) opening file <" + file + ">");
    }
    
    while (in) {
      TrajStrPtr c = std::make_shared<TrajStruct>();
      if (not c->read(in, verbose)) break;
      
      coefs[roundTime(c->time)] = c;
    }

    for (auto c : coefs) times.push_back(c.first);
  }
    

  TrajectoryData::TrajectoryData(HighFive::File& file, int stride,
				 double Tmin, double Tmax, bool verbose) :
    Coefs("trajectory", verbose)
  {
    int traj, rank;
    unsigned count;
    std::string config;
    
    file.getAttribute("name"  ).read(name);
    file.getAttribute("traj"  ).read(traj);
    file.getAttribute("rank"  ).read(rank);
    file.getAttribute("config").read(config);
    file.getDataSet  ("count" ).read(count);

    auto snaps = file.getGroup("snapshots");
    
    snaps.getDataSet("times").read(times);
    snaps.getDataSet("datatable").read(data);
      
    for (unsigned n=0; n<count; n+=stride) {
      
      if (times[n] < Tmin or times[n] > Tmax) continue;

      // Pack the data into the coefficient variable
      //
      auto coef = std::make_shared<TrajStruct>();
      coef->traj  = traj;
      coef->rank  = rank;
      coef->time  = times[n];
      coef->store.resize(traj*rank);
      coef->coefs = std::make_shared<TrajStruct::coefType>
	(coef->store.data(), traj, rank);
      *coef->coefs = data[n];

      coefs[roundTime(times[n])] = coef;
    }
  }
  
  Eigen::VectorXcd& TrajectoryData::getData(double time)
  {
    auto it = coefs.find(roundTime(time));
    
    if (it == coefs.end()) {
      
      arr.resize(0);
      
    } else {
      
      arr = it->second->store;
      
    }
    
    return arr;
  }
  
  void TrajectoryData::setData(double time, Eigen::VectorXcd& dat)
  {
    auto it = coefs.find(roundTime(time));

    if (it == coefs.end()) {
      std::ostringstream str;
      str << "TrajectoryData::setData: requested time=" << time << " not found";
      throw std::runtime_error(str.str());
    } else {
      it->second->store = dat;
      it->second->coefs = std::make_shared<CoefClasses::TrajStruct::coefType>
	(it->second->store.data(), it->second->traj, it->second->rank);

    }
  }

  void TrajectoryData::readNativeCoefs(const std::string& file, int stride,
				       double tmin, double tmax)
  {
    std::ifstream in(file);
    
    if (not in) {
      throw std::runtime_error("TrajectoryData ERROR (runtime) opening file <" + file + ">");
    }
    
    int count = 0;
    while (in) {
      TrajStrPtr c = std::make_shared<TrajStruct>();
      if (not c->read(in, verbose)) break;
      
      if (count++ % stride) continue;
      if (c->time < tmin or c->time > tmax) continue;

      coefs[roundTime(c->time)] = c;
    }

    times.clear();
    for (auto t : coefs) times.push_back(t.first);

    if (myid==0)
      std::cerr << "---- Coefs::factory: "
		<< "read ascii and created TrajectoryData"
		<< std::endl;
  }
  
  
  void TrajectoryData::WriteH5Params(HighFive::File& file)
  {
    int traj = coefs.begin()->second->traj;
    int rank = coefs.begin()->second->rank;
    
    std::string geometry = "trajectory";

    file.createAttribute<int>("traj", HighFive::DataSpace::From(traj)).write(traj);
    file.createAttribute<int>("rank", HighFive::DataSpace::From(rank)).write(rank);
  }
  
  unsigned TrajectoryData::WriteH5Times(HighFive::Group& snaps, unsigned count)
  {
    snaps.createDataSet("times", times);
    snaps.createDataSet("datatable", data);
    
    count = times.size();

    return count;
  }
  
  
  std::string TrajectoryData::getYAML()
  {
    std::string ret;
    if (coefs.size()) {
      ret = coefs.begin()->second->buf;
    }
    return ret;
  }

  bool TrajectoryData::CompareStanzas(std::shared_ptr<Coefs> check)
  {
    bool ret = true;
    
    auto other = dynamic_cast<TrajectoryData*>(check.get());
    
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
	if (v.second->traj != it->second->traj) ret = false;
	if (v.second->rank != it->second->rank) ret = false;
      }
    }
    
    if (ret) {
      std::cout << "Parameters are the same, now checking coefficients time"
		<< std::endl;
      for (auto v : coefs) {
	auto it = other->coefs.find(v.first);
	for (int m=0; m<v.second->traj; m++) {
	  for (int n=0; n<v.second->rank; n++) {
	    if ((*v.second->coefs)(m, n) != (*it->second->coefs)(m, n)) {
	      ret = false;
	    }
	  }
	}
      }
    }
    
    return ret;
  }
  

  Eigen::Tensor<double, 3>
  TrajectoryData::getAllCoefs()
  {
    Eigen::Tensor<double, 3> ret;

    auto times = Times();

    int traj = coefs.begin()->second->traj;
    int rank = coefs.begin()->second->rank;
    int ntim = times.size();

    // Resize the tensor
    ret.resize(traj, rank, ntim);

    for (int t=0; t<ntim; t++) {
      auto cof = coefs[roundTime(times[t])];
      for (int m=0; m<traj; m++) {
	for (int n=0; n<traj; n++) {
	  ret(n, m, t) = (*cof->coefs)(m, n).real();
	}
      }
    }

    return ret;
  }

  std::vector<Key> TrajectoryData::makeKeys()
  {
    std::vector<Key> ret;
    if (coefs.size()==0) return ret;

    unsigned traj = coefs.begin()->second->traj;
    for (unsigned m=0; m<traj; m++) ret.push_back({m});
    
    return ret;
  }

  TableData::TableData(const std::vector<double>& Times,
		       const std::vector<std::vector<double>>& data,
		       bool verbose) :
    data(data), Coefs("table", verbose)
  {
    times = Times;
    for (int i=0; i<times.size(); i++) {
      TblStrPtr c = std::make_shared<TblStruct>();
      c->time = times[i];
      c->cols = data[i].size();
      c->store.resize(c->cols);
      c->coefs = std::make_shared<TblStruct::coefType>(c->store.data(), c->cols);
      for (int j=0; j<c->cols; j++) (*c->coefs)(j) = data[i][j];
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
      coef->store.resize(cols);
      coef->coefs = std::make_shared<TblStruct::coefType>(coef->store.data(), cols);
      for (int i=0; i<cols; i++) (*coef->coefs)(0, i) = data[n][i];

      coefs[roundTime(times[n])] = coef;
    }
  }
  
  Eigen::VectorXcd& TableData::getData(double time)
  {
    auto it = coefs.find(roundTime(time));
    
    if (it == coefs.end()) {
      
      arr.resize(0);
      
    } else {
      
      arr = it->second->store;
      
    }
    
    return arr;
  }
  
  void TableData::setData(double time, Eigen::VectorXcd& dat)
  {
    auto it = coefs.find(roundTime(time));

    if (it == coefs.end()) {
      std::ostringstream str;
      str << "TableData::setData: requested time=" << time << " not found";
      throw std::runtime_error(str.str());
    } else {
      it->second->store = dat;
      it->second->coefs = std::make_shared<CoefClasses::TblStruct::coefType>
	(it->second->store.data(), dat.size());

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
	  if ((*v.second->coefs)(m) != (*it->second->coefs)(m)) {
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
      auto cof = coefs[roundTime(times[t])];
      for (int c=0; c<cols; c++) {
	ret(c, t) = (*cof->coefs)(c).real();
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
	// Is the set a biorthogonal basis (has the forceID attribute)
	// or general basis (fieldID attribute)?
	//
	if (h5file.hasAttribute("forceID")) {
	  if (geometry.compare("sphere")==0) {
	    coefs = std::make_shared<SphCoefs>(h5file, stride, tmin, tmax);
	  } else if (geometry.compare("cylinder")==0) {
	    coefs = std::make_shared<CylCoefs>(h5file, stride, tmin, tmax);
	  } else if (geometry.compare("slab")==0) {
	    coefs = std::make_shared<SlabCoefs>(h5file, stride, tmin, tmax);
	  } else if (geometry.compare("cube")==0) {
	    coefs = std::make_shared<CubeCoefs>(h5file, stride, tmin, tmax);
	  } else if (geometry.compare("table")==0) {
	    coefs = std::make_shared<TableData>(h5file, stride, tmin, tmax);
	  } else {
	    throw std::runtime_error("Coefs::factory: unknown H5 coefficient file geometry: " + geometry);
	  }
	} else if (h5file.hasAttribute("fieldID")) {
	  // Use the fieldID to choose the basis type
	  std::string field;
	  HighFive::Attribute fieldID = h5file.getAttribute("fieldID");
	  fieldID.read(field);

	  if (field.compare("spherical field")>0) {
	    coefs = std::make_shared<SphFldCoefs>(h5file, stride, tmin, tmax);
	  } else if (field.compare("polar field")>0) {
	    coefs = std::make_shared<CylFldCoefs>(h5file, stride, tmin, tmax);
	  } else {
	    throw std::runtime_error("Coefs::factory: unknown H5 coefficient file fieldID: " + field);
	  }
	}
      } catch (HighFive::Exception& err) {
	std::string msg("Coefs::factory: error reading HDF5 file, ");
	throw std::runtime_error(msg + err.what());
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
    } else if (dynamic_cast<SphFldStruct*>(coef.get())) {
      ret = std::make_shared<SphFldCoefs>();
    } else if (dynamic_cast<CylFldStruct*>(coef.get())) {
      ret = std::make_shared<CylFldCoefs>();
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
	    if ((*v.second->coefs)(m, n) != (*it->second->coefs)(m, n)) {
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
    // Sanity check: throw runtime error if there are no coefficient
    // sets
    //
    if (Times().size() == 0) {
      throw std::runtime_error
	("Coefs::WriteH5Coefs: "
	 "we have NO coefficient sets...continuing without writing"
	 );
    }

    // Write coefficients
    //
    try {
      // Create a new hdf5 file
      //
      HighFive::File file(prefix,
			  HighFive::File::ReadWrite |
			  HighFive::File::Create);
      
      // Write the Version string
      //
      file.createAttribute<std::string>("CoefficientOutputVersion", HighFive::DataSpace::From(CoefficientOutputVersion)).write(CoefficientOutputVersion);

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
    if (not p) throw std::runtime_error("CylCoefs::add: Null coefficient structure, nothing added!");

    Mmax = p->mmax;
    Nmax = p->nmax;
    coefs[roundTime(coef->time)] = p;
  }

  void CubeCoefs::add(CoefStrPtr coef)
  {
    auto p = std::dynamic_pointer_cast<CubeStruct>(coef);
    if (not p) throw std::runtime_error("CubeCoefs::add: Null coefficient structure, nothing added!");

    NmaxX = p->nmaxx;
    NmaxY = p->nmaxy;
    NmaxZ = p->nmaxz;
    coefs[roundTime(coef->time)] = p;
  }

  void SlabCoefs::add(CoefStrPtr coef)
  {
    auto p = std::dynamic_pointer_cast<SlabStruct>(coef);
    if (not p) throw std::runtime_error("SlabCoefs::add: Null coefficient structure, nothing added!");

    NmaxX = p->nmaxx;
    NmaxY = p->nmaxy;
    NmaxZ = p->nmaxz;
    coefs[roundTime(coef->time)] = p;
  }

  void TableData::add(CoefStrPtr coef)
  {
    auto p = std::dynamic_pointer_cast<TblStruct>(coef);
    if (not p) throw std::runtime_error("TableData::add: Null coefficient structure, nothing added!");

    coefs[roundTime(coef->time)] = p;
  }

  void TrajectoryData::add(CoefStrPtr coef)
  {
    auto p = std::dynamic_pointer_cast<TrajStruct>(coef);
    if (not p) throw std::runtime_error("TrajectoryData::add: Null coefficient structure, nothing added!");

    coefs[roundTime(coef->time)] = p;
  }


  void SphFldCoefs::add(CoefStrPtr coef)
  {
    auto p = std::dynamic_pointer_cast<SphFldStruct>(coef);
    Nfld = p->nfld;
    Lmax = p->lmax;
    Nmax = p->nmax;
    coefs[roundTime(coef->time)] = p;
  }

  void CylFldCoefs::add(CoefStrPtr coef)
  {
    auto p = std::dynamic_pointer_cast<CylFldStruct>(coef);
    if (not p) throw std::runtime_error("CylFldCoefs::add: Null coefficient structure, nothing added!");

    Nfld = p->nfld;
    Mmax = p->mmax;
    Nmax = p->nmax;
    coefs[roundTime(coef->time)] = p;
  }

  Eigen::VectorXcd& SphFldCoefs::getData(double time)
  {
    auto it = coefs.find(roundTime(time));

    if (it == coefs.end()) {
      arr.resize(0);
    } else {
      arr = it->second->store;
    }
    
    return arr;
  }
  
  SphFldStruct::dataType SphFldCoefs::getMatrix(double time)
  {
    auto it = coefs.find(roundTime(time));

    if (it == coefs.end()) {
      arr.resize(0);
    } else {
      arr = it->second->store;
      int ldim = (Lmax+1)*(Lmax+2)/2;
      mat = std::make_shared<SphFldStruct::coefType>
	(arr.data(), Nfld, ldim, Nmax); 
    }
    
    return *mat;
  }
  
  void SphFldCoefs::setData(double time, Eigen::VectorXcd& dat)
  {
    auto it = coefs.find(roundTime(time));

    if (it == coefs.end()) {
      std::ostringstream str;
      str << "SphFldCoefs::setData: requested time=" << time << " not found";
      throw std::runtime_error(str.str());
    } else {
      it->second->store = dat;
      it->second->coefs = std::make_shared<SphFldStruct::coefType>
	(it->second->store.data(), Nfld, (Lmax+1)*(Lmax+2)/2, Nmax);
    }
  }
  
  void SphFldCoefs::setMatrix(double time, SphFldStruct::dataType& dat)
  {
    auto it = coefs.find(roundTime(time));

    if (it == coefs.end()) {
      std::ostringstream str;
      str << "SphFldCoefs::setMatrix: requested time=" << time << " not found";
      throw std::runtime_error(str.str());
    } else {
      it->second->allocate();
      *it->second->coefs = dat;
    }
  }
  
  Eigen::Tensor<std::complex<double>, 4> SphFldCoefs::getAllCoefs()
  {
    Eigen::Tensor<std::complex<double>, 4> ret;

    auto times = Times();
    int ntim = times.size();

    // Resize the tensor
    ret.resize(Nfld, (Lmax+1)*(Lmax+2)/2, Nmax, ntim);

    for (int t=0; t<ntim; t++) {
      auto & cof = *(coefs[roundTime(times[t])]->coefs);
      for (int i=0; i<Nfld; i++) {
	for (int l=0; l<(Lmax+2)*(Lmax+1)/2; l++) {
	  for (int n=0; n<Nmax; n++) {
	    ret(i, l, n, t) = cof(i, l, n);
	  }
	}
      }
    }

    return ret;
  }


  std::vector<Key> SphFldCoefs::makeKeys(Key k)
  {
    std::vector<Key> ret;
    if (coefs.size()==0) return ret;

    // Sanity
    if (k.size()) {
      k[0] = std::max<unsigned>(k[0], 0);
      k[0] = std::min<unsigned>(k[0], Nfld);
    }

    if (k.size()>1) {
      k[1] = std::max<unsigned>(k[1], 0);
      k[1] = std::min<unsigned>(k[1], Lmax);
    }

    if (k.size()>2) {
      k[2] = std::max<unsigned>(k[2], 0);
      k[2] = std::min<unsigned>(k[2], k[1]);
    }

    // Four options
    // 1. return all keys for a fixed i, l, m
    // 2. return all keys for a fixed i, l
    // 3. return all keys for a fixed i
    // 4. return all keys

    // Option 4
    if (k.size()==0) {
      for (unsigned i=0; i<Nfld; i++)
	for (unsigned l=0; l<=Lmax; l++)
	  for (unsigned m=0; m<=l; m++) 
	    for (unsigned n=0; n<Nmax; n++) ret.push_back({i, l, m, n});
    }
    // Option 3
    else if (k.size()==1) {
      for (unsigned l=0; l<=Lmax; l++)
	for (unsigned m=0; m<=l; m++) 
	  for (unsigned n=0; n<Nmax; n++) ret.push_back({k[0], l, m, n});
    }
    // Option 2
    else if (k.size()==2) {
      for (unsigned m=0; m<=k[1]; m++) 
	for (unsigned n=0; n<Nmax; n++) ret.push_back({k[0], k[1], m, n});
    }
    // Option 1
    else if (k.size()==3) {
      for (unsigned n=0; n<Nmax; n++) ret.push_back({k[0], k[1], k[2], n});
    }
    // Bad sub key?
    else {
      throw std::runtime_error
	("SphVelCoefs::makeKeys: the subkey must have rank 0, 1 or 2");
    }

    return ret;
  }

  std::string SphFldCoefs::getYAML()
  {
    std::string ret;
    if (coefs.size()) {
      ret = coefs.begin()->second->buf;
    }
    return ret;
  }
  
  void SphFldCoefs::WriteH5Params(HighFive::File& file)
  {
    // Identify myself
    //
    std::string fieldID("spherical velocity orthgonal function coefficients");
    file.createAttribute<std::string>("fieldID", HighFive::DataSpace::From(fieldID)).write(fieldID);
    
    double scale = coefs.begin()->second->scale;

    // Write the remaining parameters
    //
    file.createAttribute<int>   ("nfld",  HighFive::DataSpace::From(Nfld)  ).write(Nfld);
    file.createAttribute<int>   ("lmax",  HighFive::DataSpace::From(Lmax)  ).write(Lmax);
    file.createAttribute<int>   ("nmax",  HighFive::DataSpace::From(Nmax)  ).write(Nmax);
    file.createAttribute<double>("scale", HighFive::DataSpace::From(scale)).write(scale);
    file.createAttribute<int>   ("dof",   HighFive::DataSpace::From(dof)   ).write(dof);
  }
  
  unsigned SphFldCoefs::WriteH5Times(HighFive::Group& snaps, unsigned count)
  {
    for (auto c : coefs) {
      auto C = c.second;
      
      std::ostringstream stim;
      stim << std::setw(8) << std::setfill('0') << std::right << count++;
      HighFive::Group stanza = snaps.createGroup(stim.str());
    
      // Add time attribute
      //
      stanza.createAttribute<double>("Time", HighFive::DataSpace::From(C->time)).write(C->time);
    
      // Add a center attribute
      //
      if (C->ctr.size()>0)
	stanza.createAttribute<double>("Center", HighFive::DataSpace::From(C->ctr)).write(C->ctr);
      

      // Coefficient size (allow Eigen::Tensor to be easily recontructed from metadata)
      //
      const auto& d = C->coefs->dimensions();
      std::array<long int, 3> shape {d[0], d[1], d[2]};
      stanza.createAttribute<long int>("shape", HighFive::DataSpace::From(shape)).write(shape);
      
      // Add coefficient data from flattened tensor
      //
      HighFive::DataSet dataset = stanza.createDataSet("coefficients", C->store);
    }
    
    return count;
  }

  bool SphFldCoefs::CompareStanzas(CoefsPtr check)
  {
    bool ret = true;
    
    auto other = std::dynamic_pointer_cast<SphFldCoefs>(check);
    
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
	if (v.second->nfld != it->second->nfld) ret = false;
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
	auto & cv = *(v.second->coefs);
	auto & ci = *(it->second->coefs);
	const auto & d = cv.dimensions();
	for (int n=0; n<d[0]; n++) {
	  for (int i=0; i<d[1]; i++) {
	    for (int j=0; j<d[2]; j++) {
	      if (cv(n, i, j) != ci(n, i, j)) {
		std::cout << "Coefficient (" << n << ", " << i << ", " << j << ")  "
			  << cv(n, i, j) << " != "
			  << ci(n, i, j) << std::endl;
		ret = false;
	      }
	    }
	  }
	}
      }
    }
    
    return ret;
  }
  

  std::string CylFldCoefs::getYAML()
  {
    std::string ret;
    if (coefs.size()) {
      ret = coefs.begin()->second->buf;
    }
    return ret;
  }
  

  void CylFldCoefs::WriteH5Params(HighFive::File& file)
  {
    // Identify myself
    //
    std::string fieldID("polar velocity orthgonal function coefficients");
    file.createAttribute<std::string>("fieldID", HighFive::DataSpace::From(fieldID)).write(fieldID);
    
    double scale = coefs.begin()->second->scale;

    // Write the remaining parameters
    //
    file.createAttribute<int>   ("nfld",  HighFive::DataSpace::From(Nfld)  ).write(Nfld);
    file.createAttribute<int>   ("mmax",  HighFive::DataSpace::From(Mmax)  ).write(Mmax);
    file.createAttribute<int>   ("nmax",  HighFive::DataSpace::From(Nmax)  ).write(Nmax);
    file.createAttribute<double>("scale", HighFive::DataSpace::From(scale)).write(scale);
    file.createAttribute<int>   ("dof",   HighFive::DataSpace::From(dof)   ).write(dof);
  }
  
  unsigned CylFldCoefs::WriteH5Times(HighFive::Group& snaps, unsigned count)
  {
    for (auto c : coefs) {
      auto C = c.second;
      
      std::ostringstream stim;
      stim << std::setw(8) << std::setfill('0') << std::right << count++;
      HighFive::Group stanza = snaps.createGroup(stim.str());
    
      // Add time attribute
      //
      stanza.createAttribute<double>("Time", HighFive::DataSpace::From(C->time)).write(C->time);
    
      // Add a center attribute
      //
      if (C->ctr.size()>0)
	stanza.createAttribute<double>("Center", HighFive::DataSpace::From(C->ctr)).write(C->ctr);
      

      // Coefficient size (allow Eigen::Tensor to be easily recontructed from metadata)
      //
      const auto& d = C->coefs->dimensions();
      std::array<long int, 3> shape {d[0], d[1], d[2]};
      stanza.createAttribute<long int>("shape", HighFive::DataSpace::From(shape)).write(shape);
      
      // Add coefficient data from flattened tensor
      //
      HighFive::DataSet dataset = stanza.createDataSet("coefficients", C->store);
    }
    
    return count;
  }


  bool CylFldCoefs::CompareStanzas(CoefsPtr check)
  {
    bool ret = true;
    
    auto other = std::dynamic_pointer_cast<CylFldCoefs>(check);
    
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
	if (v.second->nfld != it->second->nfld) ret = false;
	if (v.second->mmax != it->second->mmax) ret = false;
	if (v.second->nmax != it->second->nmax) ret = false;
	if (v.second->time != it->second->time) ret = false;
      }
    }
    
    if (ret) {
      std::cout << "Parameters are the same, now checking coefficients"
		<< std::endl;
      for (auto v : coefs) {
	auto it = other->coefs.find(v.first);
	auto & cv = *(v.second->coefs);
	auto & ci = *(it->second->coefs);
	const auto & d = cv.dimensions();
	for (int n=0; n<d[0]; n++) {
	  for (int i=0; i<d[1]; i++) {
	    for (int j=0; j<d[2]; j++) {
	      if (cv(n, i, j) != ci(n, i, j)) {
		std::cout << "Coefficient (" << n << ", " << i << ", " << j << ")  "
			  << cv(n, i, j) << " != "
			  << ci(n, i, j) << std::endl;
		ret = false;
	      }
	    }
	  }
	}
      }
    }
    
    return ret;
  }
  

  Eigen::MatrixXd& SphFldCoefs::Power(int min, int max)
  {
    if (coefs.size()) {
      
      int nfld = coefs.begin()->second->nfld;
      int lmax = coefs.begin()->second->lmax;
      int nmax = coefs.begin()->second->nmax;
      power.resize(coefs.size(), lmax+1);
      power.setZero();
      
      int T=0;
      for (auto v : coefs) {
	for (int i=1; i<nfld; i++) {
	  for (int l=0, L=0; l<=lmax; l++) {
	    for (int m=0; m<=l; m++, L++) {
	      for (int n=std::max<int>(0, min); n<std::min<int>(nmax, max); n++) {
		power(T, l) += std::norm((*v.second->coefs)(i, l, n));
	      }
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
  
  Eigen::MatrixXd& CylFldCoefs::Power(int min, int max)
  {
    if (coefs.size()) {
      
      int nfld = coefs.begin()->second->nfld;
      int mmax = coefs.begin()->second->mmax;
      int nmax = coefs.begin()->second->nmax;
      power.resize(coefs.size(), mmax+1);
      power.setZero();
      
      int T=0;
      for (auto v : coefs) {
	for (int i=1; i<nfld; i++) {
	  for (int m=0; m<=mmax; m++) {
	    for (int n=std::max<int>(0, min); n<std::min<int>(nmax, max); n++) {
	      power(T, m) += std::norm((*v.second->coefs)(i, m, n));
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
  

  Eigen::VectorXcd& CylFldCoefs::getData(double time)
  {
    auto it = coefs.find(roundTime(time));

    if (it == coefs.end()) {
      arr.resize(0);
    } else {
      arr = it->second->store;
    }
    
    return arr;
  }
  
  CylFldStruct::dataType CylFldCoefs::getMatrix(double time)
  {
    auto it = coefs.find(roundTime(time));

    if (it == coefs.end()) {
      arr.resize(0);
    } else {
      arr = it->second->store;
      int mdim = Mmax + 1;
      mat = std::make_shared<CylFldStruct::coefType>(arr.data(), Nfld, mdim, Nmax); 
    }
    
    return *mat;
  }
  
  void CylFldCoefs::setData(double time, Eigen::VectorXcd& dat)
  {
    auto it = coefs.find(roundTime(time));

    if (it == coefs.end()) {
      std::ostringstream str;
      str << "CylFldCoefs::setData: requested time=" << time << " not found";
      throw std::runtime_error(str.str());
    } else {
      it->second->store = dat;
      it->second->coefs = std::make_shared<SphFldStruct::coefType>
	(it->second->store.data(), Nfld, Mmax+1, Nmax);
    }
  }
  
  void CylFldCoefs::setMatrix(double time, CylFldStruct::dataType& dat)
  {
    auto it = coefs.find(roundTime(time));

    if (it == coefs.end()) {
      std::ostringstream str;
      str << "CylFldCoefs::setMatrix: requested time=" << time << " not found";
      throw std::runtime_error(str.str());
    } else {
      it->second->allocate();
      *it->second->coefs = dat;
    }
  }
  
  Eigen::Tensor<std::complex<double>, 4> CylFldCoefs::getAllCoefs()
  {
    Eigen::Tensor<std::complex<double>, 4> ret;

    auto times = Times();
    int ntim = times.size();

    // Resize the tensor
    ret.resize(Nfld, Mmax+1, Nmax, ntim);

    for (int t=0; t<ntim; t++) {
      auto & cof = *(coefs[roundTime(times[t])]->coefs);
      for (int i=0; i<Nfld; i++) {
	for (int m=0; m<=Mmax; m++) {
	  for (int n=0; n<Nmax; n++) {
	    ret(i, m, n, t) = cof(i, m, n);
	  }
	}
      }
    }

    return ret;
  }

  std::vector<Key> CylFldCoefs::makeKeys(Key k)
  {
    std::vector<Key> ret;
    if (coefs.size()==0) return ret;

    // Sanity
    if (k.size()) {
      k[0] = std::max<unsigned>(k[0], 0);
      k[0] = std::min<unsigned>(k[0], Nfld);
    }

    if (k.size()>1) {
      k[1] = std::max<unsigned>(k[1], 0);
      k[1] = std::min<unsigned>(k[1], Mmax);
    }

    // Three options
    // 1. return all keys for a fixed i, m
    // 2. return all keys for a fixed i
    // 3. return all keys

    // Option 3
    if (k.size()==0) {
      for (unsigned i=0; i<3; i++) 
	for (unsigned m=0; m<=Mmax; m++)
	  for (unsigned n=0; n<Nmax; n++) ret.push_back({i, m, n});
    }
    // Option 2
    else if (k.size()==1) {
      for (unsigned m=0; m<=Mmax; m++) 
	for (unsigned n=0; n<Nmax; n++) ret.push_back({k[0], m, n});
    }
    // Option 1
    else if (k.size()==2) {
      for (unsigned n=0; n<Nmax; n++) ret.push_back({k[0], k[1],  n});
    }
    // Bad sub key?
    else {
      throw std::runtime_error
	("SphVelCoefs::makeKeys: the subkey must have rank 0, 1 or 2");
    }

    return ret;
  }

}
// END namespace CoefClasses
