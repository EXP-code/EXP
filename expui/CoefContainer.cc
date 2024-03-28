#include <algorithm>
#include <iostream>
#include <iomanip>
#include <fstream>
#include <cctype>

#include <yaml-cpp/yaml.h>

#include "CoefContainer.H"

namespace MSSA
{
  //! An index key
  using Key = std::vector<unsigned>;
  
  // Constructor from lists
  CoefDB::CoefDB(const std::string& name, CoefClasses::CoefsPtr coefs,
		    const std::vector<Key>& keys0, const std::vector<Key>& bkeys0, 
		    unsigned index, unsigned stride, double tmin, double tmax)
    :
    name(name), coefs(coefs), keys0(keys0), bkeys0(bkeys0), index(index),
    stride(stride), tmin(tmin), tmax(tmax)
  {
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
  CoefClasses::CoefsPtr CoefDB::endUpdate()
  {
    // Make a new Coefs instance
    //
    if (dynamic_cast<CoefClasses::SphCoefs*>(coefs.get())) {
      unpack_sphere();
    }
    else if (dynamic_cast<CoefClasses::CylCoefs*>(coefs.get())) {
      unpack_cylinder();
    }
    else if (dynamic_cast<CoefClasses::SlabCoefs*>(coefs.get())) {
      unpack_slab();
    }
    else if (dynamic_cast<CoefClasses::CubeCoefs*>(coefs.get())) {
      unpack_cube();
    }
    else if (dynamic_cast<CoefClasses::TableData*>(coefs.get())) {
      unpack_table();
    }
    else {
      throw std::runtime_error("CoefDB::unpack_channels(): can not reflect coefficient type");
    }

    return coefs;
  }
  
  void CoefDB::background()
  {
    if (dynamic_cast<CoefClasses::SphCoefs*>(coefs.get()))
      restore_background_sphere();
    else if (dynamic_cast<CoefClasses::CylCoefs*>(coefs.get()))
      restore_background_cylinder();
    else if (dynamic_cast<CoefClasses::TableData*>(coefs.get()))
      { } // Do nothing
    else {
      throw std::runtime_error("CoefDB::background(): can not reflect coefficient type");
    }
  }

  void CoefDB::pack_channels()
  {
    if (dynamic_cast<CoefClasses::SphCoefs*>(coefs.get()))
      pack_sphere();
    else if (dynamic_cast<CoefClasses::CylCoefs*>(coefs.get()))
      pack_cylinder();
    else if (dynamic_cast<CoefClasses::SlabCoefs*>(coefs.get()))
      pack_slab();
    else if (dynamic_cast<CoefClasses::CubeCoefs*>(coefs.get()))
      pack_cube();
    else if (dynamic_cast<CoefClasses::TableData*>(coefs.get()))
      pack_table();
    else {
      throw std::runtime_error("CoefDB::pack_channels(): can not reflect coefficient type");
    }
  }

  void CoefDB::pack_cylinder()
  {
    auto cur = dynamic_cast<CoefClasses::CylCoefs*>(coefs.get());

    times = cur->Times();
    complexKey = true;

    auto cf = dynamic_cast<CoefClasses::CylStruct*>( cur->getCoefStruct(times[0]).get() );

    int mmax = cf->mmax;
    int nmax = cf->nmax;
    int ntimes = times.size();
    
    // Promote desired keys into c/s pairs
    //
    keys.clear();
    for (auto v : keys0) {
      // Sanity check rank
      //
      if (v.size() != 2) {
	std::ostringstream sout;
	sout << "CoefDB::pack_cylinder: key vector should have rank 2; "
	     << "found rank " << v.size() << " instead";
	throw std::runtime_error(sout.str());
      }
      // Sanity check values
      //
      else if (v[0]>=0 and v[0]<=mmax and
	  v[1]>=0 and v[1]<=nmax ) {
	auto c = v, s = v;
	c.push_back(0);
	s.push_back(1);
	keys.push_back(c);
	if (v[0]) keys.push_back(s);
      } else {
	throw std::runtime_error("CoefDB::pack_cylinder: key is out of bounds");
      }
    }

    bkeys.clear();
    for (auto v : bkeys0) {
      // Sanity check values
      //
      if (v[0]>=0 and v[0]<=mmax and
	  v[1]>=0 and v[1]<=nmax ) {
	auto c = v, s = v;
	c.push_back(0);
	s.push_back(1);
	bkeys.push_back(c);
	if (v[0]) bkeys.push_back(s);
      }
    }

    // Only pack the keys in the list
    //
    for (auto k : keys) {
      data[k].resize(ntimes);
    }

    for (auto k : bkeys) {
      data[k].resize(ntimes);
    }

    for (int t=0; t<ntimes; t++) {
      cf = dynamic_cast<CoefClasses::CylStruct*>( cur->getCoefStruct(times[t]).get() );
      for (auto k : keys) {
	if (k[2]==0)
	  data[k][t] = (*cf->coefs)(k[0], k[1]).real();
	else
	  data[k][t] = (*cf->coefs)(k[0], k[1]).imag();
      }

      for (auto k : bkeys) {
	if (k[2]==0)
	  data[k][t] = (*cf->coefs)(k[0], k[1]).real();
	else
	  data[k][t] = (*cf->coefs)(k[0], k[1]).imag();
      }
    }
  }

  void CoefDB::unpack_cylinder()
  {
    for (int i=0; i<times.size(); i++) {
      auto cf = dynamic_cast<CoefClasses::CylStruct*>( coefs->getCoefStruct(times[i]).get() );
      
      for (auto k : keys0) {
	auto c = k, s = k;
	c.push_back(0); s.push_back(1);

	int m = k[0], n = k[1];

	if (m==0) (*cf->coefs)(m, n) = {data[c][i], 0.0};
	else      (*cf->coefs)(m, n) = {data[c][i], data[s][i]};
      }
      // END key loop
    }
    // END time loop
  }

  void CoefDB::pack_sphere()
  {
    auto cur = dynamic_cast<CoefClasses::SphCoefs*>(coefs.get());

    times = cur->Times();
    complexKey = true;

    auto cf = dynamic_cast<CoefClasses::SphStruct*>( cur->getCoefStruct(times[0]).get() );

    int lmax   = cf->lmax;
    int nmax   = cf->nmax;
    int ntimes = times.size();
    
    // Make extended key list
    //
    keys.clear();
    for (auto k : keys0) {
      // Sanity check rank
      //
      if (k.size() != 3) {
	std::ostringstream sout;
	sout << "CoefDB::pack_sphere: key vector should have rank 3; "
	     << "found rank " << k.size() << " instead";
	throw std::runtime_error(sout.str());
      }
      // Sanity check values
      //
      else if (k[0] <= lmax and k[0] >= 0 and
	       k[1] <= k[0] and k[1] >= 0 and
	       k[2] <= nmax and k[2] >= 0 ) {
	
	auto v = k;
	v.push_back(0);
	keys.push_back(v);
	data[v].resize(ntimes);

	if (k[1]>0) {
	  v[3] = 1;
	  keys.push_back(v);
	  data[v].resize(ntimes);
	}
      }
      else {
	throw std::runtime_error("CoefDB::pack_sphere: key is out of bounds");
      }
    }

    bkeys.clear();
    for (auto k : bkeys0) {
      // Sanity check values
      //
      if (k[0] <= lmax and k[0] >= 0 and
	  k[1] <= k[0] and k[1] >= 0 and
	  k[2] <= nmax and k[2] >= 0 ) {
	
	auto v = k;
	v.push_back(0);
	keys.push_back(v);
	data[v].resize(ntimes);

	if (k[1]>0) {
	  v[3] = 1;
	  keys.push_back(v);
	  data[v].resize(ntimes);
	}
      }
    }

    auto I = [](const Key& k) { return k[0]*(k[0]+1)/2 + k[1]; };

    for (int t=0; t<ntimes; t++) {
      cf = dynamic_cast<CoefClasses::SphStruct*>( cur->getCoefStruct(times[t]).get() );
      for (auto k : keys)  {
	auto c = (*cf->coefs)(I(k), k[2]);
	data[k][t] = c.real();
	if (k[3]) data[k][t] = c.imag();
      }

      for (auto k : bkeys)  {
	auto c = (*cf->coefs)(I(k), k[2]);
	data[k][t] = c.real();
	if (k[3]) data[k][t] = c.imag();
      }
    }
  }

  void CoefDB::unpack_sphere()
  {
    auto I = [](const Key& k) { return k[0]*(k[0]+1)/2 + k[1]; };

    for (int i=0; i<times.size(); i++) {

      auto cf = dynamic_cast<CoefClasses::SphStruct*>( coefs->getCoefStruct(times[i]).get() );
      
      for (auto k : keys0) {
	auto c = k, s = k;
	c.push_back(0);
	s.push_back(1);

	int m = k[1], n = k[2];

	if (m==0) (*cf->coefs)(I(k), n) = {data[c][i], 0.0       };
	else      (*cf->coefs)(I(k), n) = {data[c][i], data[s][i]};
      }
      // END key loop
    }
    // END time loop
  }
  
  void CoefDB::pack_slab()
  {
    auto cur = dynamic_cast<CoefClasses::SlabCoefs*>(coefs.get());

    times = cur->Times();
    complexKey = true;

    auto cf = dynamic_cast<CoefClasses::SlabStruct*>( cur->getCoefStruct(times[0]).get() );

    int nmaxx   = cf->nmaxx;
    int nmaxy   = cf->nmaxy;
    int nmaxz   = cf->nmaxz;
    int ntimes  = times.size();

    // Make extended key list
    //
    keys.clear();
    for (auto k : keys0) {
      // Sanity check rank
      //
      if (k.size() != 3) {
	std::ostringstream sout;
	sout << "CoefDB::pack_slab: key vector should have rank 3; "
	     << "found rank " << k.size() << " instead";
	throw std::runtime_error(sout.str());
      }
      // Sanity check values
      //
      else if (k[0] <= 2*nmaxx and k[0] >= 0 and
	       k[1] <= 2*nmaxy and k[1] >= 0 and
	       k[2] <  nmaxz   and k[2] >= 0 ) {
	
	auto v = k;
	v.push_back(0);
	keys.push_back(v);
	data[v].resize(ntimes);

	v[3] = 1;
	keys.push_back(v);
	data[v].resize(ntimes);
      }
      else {
	throw std::runtime_error("CoefDB::pack_slab: key is out of bounds");
      }
    }

    bkeys.clear();
    for (auto k : bkeys0) {
      // Sanity check values
      //
      if (k[0] <= 2*nmaxx and k[0] >= 0 and
	  k[1] <= 2*nmaxy and k[1] >= 0 and
	  k[2] <  nmaxz   and k[2] >= 0 ) {
	
	auto v = k;
	v.push_back(0);
	keys.push_back(v);
	data[v].resize(ntimes);

	v[3] = 1;
	keys.push_back(v);
	data[v].resize(ntimes);
      }
    }

    for (int t=0; t<ntimes; t++) {
      cf = dynamic_cast<CoefClasses::SlabStruct*>( cur->getCoefStruct(times[t]).get() );
      for (auto k : keys)  {
	auto c = (*cf->coefs)(k[0], k[1], k[2]);
	if (k[3]) data[k][t] = c.imag();
	else      data[k][t] = c.real();
      }

      for (auto k : bkeys)  {
	auto c = (*cf->coefs)(k[0], k[1], k[2]);
	if (k[3]) data[k][t] = c.imag();
	else      data[k][t] = c.real();
      }
    }
  }

  void CoefDB::pack_cube()
  {
    auto cur = dynamic_cast<CoefClasses::CubeCoefs*>(coefs.get());

    times = cur->Times();
    complexKey = true;

    auto cf = dynamic_cast<CoefClasses::CubeStruct*>( cur->getCoefStruct(times[0]).get() );

    int nmaxx   = cf->nmaxx;
    int nmaxy   = cf->nmaxy;
    int nmaxz   = cf->nmaxz;
    int ntimes  = times.size();

    // Make extended key list
    //
    keys.clear();
    for (auto k : keys0) {
      // Sanity check rank
      //
      if (k.size() != 3) {
	std::ostringstream sout;
	sout << "CoefDB::pack_cube: key vector should have rank 3; "
	     << "found rank " << k.size() << " instead";
	throw std::runtime_error(sout.str());
      }
      // Sanity check values
      //
      else if (k[0] <= 2*nmaxx and k[0] >= 0 and
	       k[1] <= 2*nmaxy and k[1] >= 0 and
	       k[2] <= 2*nmaxz and k[2] >= 0 ) {
	
	auto v = k;
	v.push_back(0);
	keys.push_back(v);
	data[v].resize(ntimes);

	v[3] = 1;
	keys.push_back(v);
	data[v].resize(ntimes);
      }
      else {
	throw std::runtime_error("CoefDB::pack_cube: key is out of bounds");
      }
    }

    bkeys.clear();
    for (auto k : bkeys0) {
      // Sanity check values
      //
      if (k[0] <= 2*nmaxx and k[0] >= 0 and
	  k[1] <= 2*nmaxy and k[1] >= 0 and
	  k[2] <= 2*nmaxz and k[2] >= 0 ) {
	
	auto v = k;
	v.push_back(0);
	keys.push_back(v);
	data[v].resize(ntimes);

	v[3] = 1;
	keys.push_back(v);
	data[v].resize(ntimes);
      }
    }

    for (int t=0; t<ntimes; t++) {
      cf = dynamic_cast<CoefClasses::CubeStruct*>( cur->getCoefStruct(times[t]).get() );
      for (auto k : keys)  {
	auto c = (*cf->coefs)(k[0], k[1], k[2]);
	if (k[3]) data[k][t] = c.imag();
	else      data[k][t] = c.real();
      }

      for (auto k : bkeys)  {
	auto c = (*cf->coefs)(k[0], k[1], k[2]);
	if (k[3]) data[k][t] = c.imag();
	else      data[k][t] = c.real();
      }
    }
  }

  void CoefDB::unpack_slab()
  {
    for (int i=0; i<times.size(); i++) {

      auto cf = dynamic_cast<CoefClasses::SlabStruct*>( coefs->getCoefStruct(times[i]).get() );
      
      for (auto k : keys0) {
	auto c = k, s = k;
	c.push_back(0);
	s.push_back(1);

	(*cf->coefs)(k[0], k[1], k[2]) = {data[c][i], data[s][i]};
      }
      // END key loop
    }
    // END time loop
  }
  
  void CoefDB::unpack_cube()
  {
    for (int i=0; i<times.size(); i++) {

      auto cf = dynamic_cast<CoefClasses::CubeStruct*>( coefs->getCoefStruct(times[i]).get() );
      
      for (auto k : keys0) {
	auto c = k, s = k;
	c.push_back(0);
	s.push_back(1);

	(*cf->coefs)(k[0], k[1], k[2]) = {data[c][i], data[s][i]};
      }
      // END key loop
    }
    // END time loop
  }
  
  void CoefDB::pack_table()
  {
    auto cur = dynamic_cast<CoefClasses::TableData*>(coefs.get());

    times = cur->Times();
    complexKey = false;

    auto cf = dynamic_cast<CoefClasses::TblStruct*>( cur->getCoefStruct(times[0]).get() );

    int cols    = cf->cols;
    int ntimes  = times.size();
    
    // Promote desired keys which are the data columns
    //
    keys.clear();
    for (auto v : keys0) {
      if (v.size() != 1) {
	std::ostringstream sout;
	sout << "CoefDB::pack_table: key vector should have rank 1; "
	     << "found rank " << v.size() << " instead";
	throw std::runtime_error(sout.str());
      }
      else if (v[0] >= cols) {
	std::ostringstream sout;
	sout << "CoefDB::pack_table: requested key=" << v[0]
	     << " exceeded the number of columns, " << cols;
	throw std::runtime_error(sout.str());
      }
      else keys.push_back(v);	// OKAY
    }

    // No bkeys for a table
    //
    bkeys.clear();

    for (unsigned c=0; c<cols; c++) {
      Key key = {c};
      data[key].resize(ntimes);
    }

    for (int t=0; t<ntimes; t++) {
      for (unsigned c=0; c<cols; c++) {
	Key key = {c};

	cf = dynamic_cast<CoefClasses::TblStruct*>( cur->getCoefStruct(times[t]).get() );
	data[key][t] = (*cf->coefs)(c).real();
      }
    }
  }

  void CoefDB::unpack_table()
  {
    for (int i=0; i<times.size(); i++) {

      auto cf = dynamic_cast<CoefClasses::TblStruct*>( coefs->getCoefStruct(times[i]).get() );

      int cols = cf->cols;

      for (unsigned c=0; c<cols; c++) {
	Key key = {c};
	(*cf->coefs)(c) = data[key][i];
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


  void CoefDB::restore_background_sphere()
  {
    auto cur = dynamic_cast<CoefClasses::SphCoefs*>(coefs.get());

    auto I = [](const Key& k) { return k[0]*(k[0]+1)/2 + k[1]; };

    for (int t=0; t<times.size(); t++) {
      // Cast the coefficient class to spherical
      auto cf = dynamic_cast<CoefClasses::SphStruct*>
	(cur->getCoefStruct(times[t]).get());

      // Get the coefficient map
      auto & ar = *(cf->coefs);

      for (auto k : bkeys)  {
	auto c = ar(I(k), k[2]);
	data[k][t] = c.real();
	if (k[3]) data[k][t] = c.imag();
      }
    }
  }

  void CoefDB::restore_background_cylinder()
  {
    auto cur = dynamic_cast<CoefClasses::CylCoefs*>(coefs.get());

    for (int t=0; t<times.size(); t++) {
      // Cast the coefficient class to cylindrical
      auto cf = dynamic_cast<CoefClasses::CylStruct*>
	(cur->getCoefStruct(times[t]).get());
      
      // Get the coefficient map
      auto & ar = *(cf->coefs);

      for (auto k : bkeys) {
	if (k[2]==0)
	  data[k][t] = ar(k[0], k[1]).real();
	else
	  data[k][t] = ar(k[0], k[1]).imag();
      }
    }
  }
  // END CoefDB::background

}
// END namespace MSSA

