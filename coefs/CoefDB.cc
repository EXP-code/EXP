#include <algorithm>
#include <iostream>
#include <iomanip>
#include <fstream>
#include <cctype>

#include <yaml-cpp/yaml.h>

#include "CoefDB.H"

namespace MSSA
{
  // Constructor from lists
  CoefDB::CoefDB(const std::string& name, Coefs::CoefsPtr coefs,
		 const std::vector<Key>& keys0, const std::vector<Key>& bkeys0, 
		 unsigned index, unsigned stride, double tmin, double tmax)
    :
    name(name), coefs(coefs), index(index),
    stride(stride), tmin(tmin), tmax(tmax)
  {
    keys  = keys0;
    bkeys = bkeys0;
    pack_channels();
  }

  // Deep copy
  std::shared_ptr<CoefDB> CoefDB::deepcopy()
  {
    auto ret = std::make_shared<CoefDB>();

    ret->name       = name;
    ret->index      = index;
    ret->keys       = keys;
    ret->bkeys      = bkeys;
    ret->complexKey = complexKey;

    // Don't try to copy a null instance
    if (coefs)  ret->coefs  = coefs-> deepcopy();

    ret->data       = data;
    ret->times      = times;

    return ret;
  }

  // Copy to workset coefficient set
  Coefs::CoefsPtr CoefDB::endUpdate()
  {
    // Make a new Coefs instance
    //
    if (dynamic_cast<Coefs::SphCoefs*>(coefs.get())) {
      unpack_sphere();
    }
    else if (dynamic_cast<Coefs::CylCoefs*>(coefs.get())) {
      unpack_cylinder();
    }
    else if (dynamic_cast<Coefs::TableData*>(coefs.get())) {
      unpack_table();
    }
    else {
      throw std::runtime_error("CoefDB::pack_channels(): can not reflect coefficient type");
    }

    return coefs;
  }
  
  

  void CoefDB::pack_channels()
  {
    if (dynamic_cast<Coefs::SphCoefs*>(coefs.get()))
      pack_sphere();
    else if (dynamic_cast<Coefs::CylCoefs*>(coefs.get()))
      pack_cylinder();
    else if (dynamic_cast<Coefs::TableData*>(coefs.get()))
      pack_table();
    else {
      throw std::runtime_error("CoefDB::pack_channels(): can not reflect coefficient type");
    }
  }

  void CoefDB::pack_cylinder()
  {
    auto cur = dynamic_cast<Coefs::CylCoefs*>(coefs.get());

    times = cur->Times();
    complexKey = true;

    auto cf = dynamic_cast<Coefs::CylStruct*>( cur->getCoefStruct(times[0]).get() );

    int mmax = cf->mmax;
    int nmax = cf->nmax;
    int ntimes = times.size();
    
    for (unsigned m=0; m<=mmax; m++) {
      for (unsigned n=0; n<nmax; n++) {
	Key key0 = {m, n, 0}, key1 = {m, n, 1};
	data[key0].resize(ntimes);
	if (m) data[key1].resize(ntimes);
      }
    }

    for (int t=0; t<ntimes; t++) {
      for (unsigned m=0; m<=mmax; m++) {
	for (unsigned n=0; n<nmax; n++) {
	  Key key0 = {m, n, 0}, key1 = {m, n, 1};
	  cf = dynamic_cast<Coefs::CylStruct*>( cur->getCoefStruct(times[t]).get() );
	  data[key0][t] = cf->coefs(m, n).real();
	  if (m) data[key1][t] = cf->coefs(m, n).imag();
	}
      }
    }

    // Promote desired keys into c/s pairs
    //
    auto keys0 = keys;
    keys.clear();
    for (auto v : keys0) {
      auto c = v, s = v;
      c.push_back(0);
      s.push_back(1);
      keys.push_back(c);
      if (v[0]) keys.push_back(s);
    }

    keys0 = bkeys;
    bkeys.clear();
    for (auto v : keys0) {
      auto c = v, s = v;
      c.push_back(0);
      s.push_back(1);
      bkeys.push_back(c);
      if (v[0]) keys.push_back(s);
    }

  }

  void CoefDB::unpack_cylinder()
  {
    for (int i=0; i<times.size(); i++) {

      auto cf = dynamic_cast<Coefs::CylStruct*>( coefs->getCoefStruct(times[i]).get() );
      
      int mmax = cf->mmax;
      int nmax = cf->nmax;
    
      for (unsigned m=0; m<=mmax; m++) {
	for (unsigned n=0; n<nmax; n++) {
	  Key key0 = {m, n, 0}, key1 = {m, n, 1};
	  
	  if (m==0) cf->coefs(m, n) = {data[key0][i], 0.0};
	  else cf->coefs(m, n) = {data[key0][i], data[key1][i]};
	}
      }
    }
    // END time loop
  }

  void CoefDB::pack_sphere()
  {
    auto cur = dynamic_cast<Coefs::SphCoefs*>(coefs.get());

    times = cur->Times();
    complexKey = true;

    auto cf = dynamic_cast<Coefs::SphStruct*>( cur->getCoefStruct(times[0]).get() );

    int lmax   = cf->lmax;
    int nmax   = cf->nmax;
    int ntimes = times.size();
    
    for (unsigned l=0, L=0; l<=lmax; l++) {
      for (unsigned m=0; m<=l; m++, L++) {
	for (unsigned n=0; n<nmax; n++) {
	  Key key0 = {l, m, n, 0}, key1 = {l, m, n, 1};
	  data[key0].resize(ntimes);
	  if (m) data[key1].resize(ntimes);
	}
      }
    }

    for (int t=0; t<ntimes; t++) {
      for (unsigned l=0, L=0; l<=lmax; l++) {
	for (unsigned m=0; m<=l; m++, L++) {
	  for (unsigned n=0; n<nmax; n++) {
	    Key key0 = {l, m, n, 0}, key1 = {l, m, n, 1};

	    cf = dynamic_cast<Coefs::SphStruct*>( cur->getCoefStruct(times[t]).get() );
	    data[key0][t] = cf->coefs(L, n).real();
	    if (m) data[key1][t] = cf->coefs(L, n).imag();
	  }
	}
      }
    }

    // Promote desired keys into c/s pairs
    //
    auto keys0 = keys;
    keys.clear();
    for (auto v : keys0) {
      auto c = v, s = v;
      c.push_back(0);
      s.push_back(1);
      keys.push_back(c);
      if (v[1]) keys.push_back(s);
    }

    keys0 = bkeys;
    bkeys.clear();
    for (auto v : keys0) {
      auto c = v, s = v;
      c.push_back(0);
      s.push_back(1);
      bkeys.push_back(c);
      if (v[1]) keys.push_back(s);
    }
  }

  void CoefDB::unpack_sphere()
  {
    for (int i=0; i<times.size(); i++) {

      auto cf = dynamic_cast<Coefs::SphStruct*>( coefs->getCoefStruct(times[i]).get() );
      
      int lmax   = cf->lmax;
      int nmax   = cf->nmax;
    
      for (unsigned l=0, L=0; l<=lmax; l++) {
	for (unsigned m=0; m<=l; m++, L++) {
	  for (unsigned n=0; n<nmax; n++) {
	    Key key0 = {l, m, n, 0}, key1 = {l, m, n, 1};
	  
	    if (m==0) cf->coefs(L, n) = {data[key0][i], 0.0          };
	    else      cf->coefs(L, n) = {data[key0][i], data[key1][i]};
	  }
	  // END n loop
	}
	// END m loop
      }
      // END l loop
    }
    // END time loop
  }
  
  void CoefDB::pack_table()
  {
    auto cur = dynamic_cast<Coefs::TableData*>(coefs.get());

    times = cur->Times();
    complexKey = false;

    auto cf = dynamic_cast<Coefs::TblStruct*>( cur->getCoefStruct(times[0]).get() );

    int cols    = cf->cols;
    int ntimes  = times.size();
    
    for (unsigned c=0; c<cols; c++) {
      Key key = {c};
      data[key].resize(ntimes);
    }

    for (int t=0; t<ntimes; t++) {
      for (unsigned c=0; c<cols; c++) {
	Key key = {c};

	cf = dynamic_cast<Coefs::TblStruct*>( cur->getCoefStruct(times[t]).get() );

	data[key][t] = cf->coefs(0, c).real();
      }
    }
  }

  void CoefDB::unpack_table()
  {
    for (int i=0; i<times.size(); i++) {

      auto cf = dynamic_cast<Coefs::TblStruct*>( coefs->getCoefStruct(times[i]).get() );

      int cols = cf->cols;

      for (unsigned c=0; c<cols; c++) {
	Key key = {c};
	cf->coefs(0, c) = data[key][i];
      }
      // End field loop
    }
    // END time loop

  }


  CoefContainer::CoefContainer(const CoefContainer& p)
  {
    // Protected data
    comps   = p.comps;
    keylist = p.keylist;
    namemap = p.namemap;

    // Public data
    tmin    = p.tmin;
    tmax    = p.tmax;
    times   = p.times;
    runtag  = p.runtag;
  }

  std::shared_ptr<CoefContainer> CoefContainer::deepcopy()
  {
    auto ret = std::make_shared<CoefContainer>();

    for (auto v : comps) ret->comps.push_back(v->deepcopy());

    ret->keylist = keylist;
    ret->namemap = namemap;

    ret->tmin    = tmin;
    ret->tmax    = tmax;
    ret->stride  = stride;
    ret->times   = times;
    ret->runtag  = runtag;

    return ret;
  }

  /*
  CoefContainer::CoefContainer(const std::string spec)
  {
    YAML::Node top = YAML::LoadFile(spec);

    // Defaults
    tmin   = -std::numeric_limits<double>::max();
    tmax   =  std::numeric_limits<double>::max();
    runtag = "mssa";
    stride = 1;

    // Overrides
    if (top["tmin"])   tmin   = top["tmin"].as<double>();
    if (top["tmax"])   tmax   = top["tmax"].as<double>();
    if (top["stride"]) stride = top["stride"].as<int>();
    if (top["runtag"]) runtag = top["runtag"].as<std::string>();

    if (top["components"]) {

      YAML::Node node = top["components"];

      for (YAML::const_iterator it=node.begin(); it!=node.end(); it++) {
	std::string name = it->first.as<std::string>();
	YAML::Node node = it->second;
	std::string ctype("sphere");
	if (node["geometry"]) {
	  ctype = node["geometry"].as<std::string>();
	  std:: transform(ctype.begin(), ctype.end(), ctype.begin(),
			  [](unsigned char c){ return std::tolower(c);});
	} else {
	  std::cout << "Geometry not specified.  You must specify 'cylinder' or 'sphere'." << std::endl;
	  exit(-1);
	}

	// Create the instances
	//
	if (ctype.find("sphere")==0) {
	  namemap[name] = comps.size();
	  comps.push_back(std::make_shared<SphDB>(name, node, comps.size(), stride, tmin, tmax));
	} else if (ctype.find("cylinder")==0) {
	  namemap[name] = comps.size();
	  comps.push_back(std::make_shared<CylDB>(name, node, comps.size(), stride, tmin, tmax));
	} else if (ctype.find("table")==0) {
	  namemap[name] = comps.size();
	  comps.push_back(std::make_shared<TableDB>(name, node, comps.size(), stride, tmin, tmax));
	} else {
	  std::cout << "Unknown geometry.  You must specify 'cylinder', 'sphere', or 'table'." << std::endl;
	  exit(-1);
	}
      }
    } else {
      std::cout << "CoefContainer: no components specified" << std::endl;
      exit(-1);
    }

    // Check times (all should be the same)
    //
				// First check lengths
    size_t tsize = comps[0]->times.size();
    for (size_t n=1; n<comps.size(); n++) {
      if (tsize != comps[n]->times.size()) {
	std::cout << "CoefContainer: times lengths do not agree!" << std::endl;
	exit(-4);
      }
				// Now check time values
      for (size_t t=0; t<tsize; t++) {
	if (fabs(comps[0]->times[t] - comps[n]->times[t]) > 1.0e-8) {
	  std::cout << "CoefContainer: times disagree for indices 0 and "
		    << n << std::endl;
	  exit(-5);
	}
      }
    }

    // Waste a little space for convenience
    //
    times = comps[0]->times;

    // Make key list
    //
    for (size_t n=0; n<comps.size(); n++) {
      for (size_t k=0; k<comps[n]->keys.size(); k++) {
	std::vector<unsigned> keyp = comps[n]->keys[k];
	keyp.push_back(comps[n]->index);
	keylist.push_back(keyp);
      }
    }
  }
  */

  CoefContainer::CoefContainer(const mssaConfig& config, const std::string spec)
  {
    YAML::Node top;
    if (spec.size()) top = YAML::Load(spec);

    // Defaults
    tmin   = -std::numeric_limits<double>::max();
    tmax   =  std::numeric_limits<double>::max();
    runtag = "mssa";
    stride = 1;

    // Overrides
    if (top["tmin"])   tmin   = top["tmin"].as<double>();
    if (top["tmax"])   tmax   = top["tmax"].as<double>();
    if (top["stride"]) stride = top["stride"].as<int>();
    if (top["runtag"]) runtag = top["runtag"].as<std::string>();

    int index = 0;
    for (auto v : config) {
      namemap[v.first] = index;
      comps.push_back(std::make_shared<CoefDB>(v.first, 
					       std::get<0>(v.second),
					       std::get<1>(v.second),
					       std::get<2>(v.second),
					       index, stride, tmin, tmax));
      index++;			// Update index
    }

    // Check times (all should be the same)
    //
				// First check lengths
    size_t tsize = comps[0]->times.size();
    for (size_t n=1; n<comps.size(); n++) {
      if (tsize != comps[n]->times.size()) {
	std::cout << "CoefContainer: times lengths do not agree!" << std::endl;
	exit(-4);
      }
				// Now check time values
      for (size_t t=0; t<tsize; t++) {
	if (fabs(comps[0]->times[t] - comps[n]->times[t]) > 1.0e-8) {
	  std::cout << "CoefContainer: times disagree for indices 0 and "
		    << n << std::endl;
	  exit(-5);
	}
      }
    }

    // Waste a little space for convenience
    //
    times = comps[0]->times;

    // Make key list
    //
    for (size_t n=0; n<comps.size(); n++) {
      for (size_t k=0; k<comps[n]->keys.size(); k++) {
	std::vector<unsigned> keyp = comps[n]->keys[k];
	keyp.push_back(comps[n]->index);
	keylist.push_back(keyp);
      }
    }
  }
  // END CoefContainer constructor

}
// END namespace MSSA

