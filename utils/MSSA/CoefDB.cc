#include <algorithm>
#include <iostream>
#include <iomanip>
#include <fstream>
#include <cctype>

#include "CoefDB.H"

namespace CoefDB {

  //! Constructor
  CoefDB::CoefDB(const std::string& name, YAML::Node& node, unsigned index) :
    name(name), index(index)
  {
    if (node["keys"]) {
      YAML::Node knode = node["keys"];
      for (YAML::const_iterator kt=knode.begin(); kt!=knode.end(); kt++) {
	keys.push_back(kt->as<std::vector<unsigned>>());
      }
    } else {
      std::cout << "CoefDB: no keys found" << std::endl;
      exit(-3);
    }
  }

  //! Constructor
  CylDB::CylDB(const std::string& name, YAML::Node& node, unsigned index,
	       unsigned stride, double tmin, double tmax) :
    CoefDB(name, node, index)
  {
    std::string file;
    
    if (node["filename"]) {
      file = node["filename"].as<std::string>();
    } else {
      std::cout << "CylDB: no coefficient filename found" << std::endl;
      exit(-4);
    }

    read(file, stride, tmin, tmax);
  }


  void CylDB::read(const std::string& file, unsigned stride, double tmin, double tmax)
  {
    std::map<std::string, int> ret;
    std::ifstream in(file);

    // Check for good opened file
    //
    if (!in) {
      std::cerr << "CylDB: could not open coefficient file <"
		<< file << ">" << std::endl;
      exit(-2);
    }
    
    unsigned counter = 0;

    while (in.good()) {
      CylCoefsPtr c = std::make_shared<CylCoefs>();
      if (not c->read(in)) break;
      if (c->time >= tmin and c->time <= tmax) {
	if (counter++ % stride == 0) data[c->time] = c;
      }
    }
    
    mmax   = data.begin()->second->mmax;
    nmax   = data.begin()->second->nmax;
    ntimes = data.size();
    
    for (auto v : data) {
      times.push_back(v.second->time);
    }

    for (unsigned m=0; m<=mmax; m++) {
      for (unsigned n=0; n<nmax; n++) {
	Key key0 = {m, n, 0}, key1 = {m, n, 1};
	coefs[key0].resize(ntimes);
	if (m) coefs[key1].resize(ntimes);
	for (int i=0; i<times.size(); i++) {
	  coefs[key0][i] = data[times[i]]->cos_c[m][n];
	  if (m) coefs[key1][i] = data[times[i]]->sin_c[m][n];
	}
      }
    }
  }


  //! Constructor
  SphDB::SphDB(const std::string& name, YAML::Node& node, unsigned index,
	       unsigned stride, double tmin, double tmax) :
    CoefDB(name, node, index)
  {
    std::string file;
    
    if (node["filename"]) {
      file = node["filename"].as<std::string>();
    } else {
      std::cout << "CylDB: no coefficient filename found" << std::endl;
      exit(-4);
    }

    read(file, stride, tmin, tmax);
  }


  void SphDB::read(const std::string& file, unsigned stride, double tmin, double tmax)
  {
    std::map<std::string, int> ret;
    std::ifstream in(file);
    
    // Check for good opened file
    if (!in) {
      std::cerr << "SphDB: could not open coefficient file <"
		<< file << ">" << std::endl;
    exit(-2);
    }

    unsigned counter = 0;

    while (in.good()) {
      SphCoefsPtr c = std::make_shared<SphCoefs>();
      if (not c->read(in)) break;
      if (c->header.tnow >= tmin and c->header.tnow <= tmax) {
	if (counter++ % stride == 0) data[c->header.tnow] = c;
      }
    }
    
    lmax   = data.begin()->second->header.Lmax;
    nmax   = data.begin()->second->header.nmax;
    ntimes = data.size();
    
    for (auto v : data) {
      times.push_back(v.second->header.tnow);
    }
    
    for (unsigned l=0; l<=lmax; l++) {
      for (unsigned m=0; m<=l; m++) {
	for (unsigned n=0; n<nmax; n++) {
	  LMkey lmk  = {l, m};
	  Key   key0 = {l, m, n, 0}, key1 = {l, m, n, 1};
	  coefs[key0].resize(ntimes);
	  if (m) coefs[key1].resize(ntimes);
	  for (int i=0; i<times.size(); i++) {
	    coefs[key0][i] = data[times[i]]->cos_c[lmk][n];
	    if (m) coefs[key1][i] = data[times[i]]->sin_c[lmk][n];
	  }
	}
      }
    }
  }

  void CylDB::write(int ncomp, const TrendType& type, bool useMean,
		    double totPow, double totVar,
		    std::map<Key, double>& mean,
		    std::map<Key, double>& var,
		    std::map<Key, Eigen::MatrixXd>& RC,
		    const std::vector<std::set<int>>& groups,
		    const std::vector<double>& times,
		    const std::string& runtag)
  {
    const unsigned int magic_word = 0x5ecede5;

    std::ostringstream fname;
    fname << runtag << "." << name << ".recon_bin";
    std::ofstream out (fname.str());

    if (out) {
      
      int numT = times.size();
      int NMAX = data.begin()->second->nmax;
      int ngrp = groups.size();

      std::set<int> Mset;
      for (auto k : keys) Mset.insert(k[0]);

      std::vector<int> MM;
      for (auto s : Mset) MM.push_back(s);
      int numM = MM.size();

      out.write((const char *)&magic_word, sizeof(unsigned int));
      out.write((const char *)&numT,       sizeof(int));
      out.write((const char *)&NMAX,       sizeof(int));
      out.write((const char *)&numM,       sizeof(int));
      out.write((const char *)&ngrp,       sizeof(int));
      out.write((const char *)&MM[0],      sizeof(int)*numM);
      out.write((const char *)&times[0],   sizeof(double)*numT);

      for (auto v : groups) {
	for (unsigned i=0; i<numT; i++) {
	  for (auto M : MM) {
	    for (unsigned n=0; n<NMAX; n++) {
	      unsigned mm = static_cast<unsigned>(M);
	      Key c = {mm, n, 0, index}, s = {mm, n, 1, index};
	      auto rc = RC.find(c);
	      auto rs = RC.find(s);
	      double valC = 0.0, valS = 0.0;
	      for (auto u : v) {
		if (u < ncomp) {
		  if (rc != RC.end()) {
		    valC += RC[c](i, u)*var[c] + mean[c];
		  }
		  if (rs != RC.end()) {
		    valS += RC[s](i, u)*var[s] + mean[s];
		  }
		}
	      }
	      out.write((const char *)&valC, sizeof(double));
	      out.write((const char *)&valS, sizeof(double));
	    }
	  }
	}
      }
    }
  }

  void SphDB::write(int ncomp, const TrendType& type, bool useMean,
		    double totPow, double totVar,
		    std::map<Key, double>& mean,
		    std::map<Key, double>& var,
		    std::map<Key, Eigen::MatrixXd>& RC,
		    const std::vector<std::set<int>>& groups,
		    const std::vector<double>& times,
		    const std::string& runtag)
  {
    const unsigned int magic_word = 0x5ecede4;

    std::ostringstream fname;
    fname << runtag << "." << name << ".recon_bin";
    std::ofstream out (fname.str());

    std::set<std::pair<int, int>> LMset;
    for (auto k : keys) LMset.insert({k[0], k[1]});

    if (out) {

      unsigned LMsize = LMset.size();
      int       Gsize = groups.size();
      int        numT = times.size();
      int        nmax = data.begin()->second->header.nmax;

      out.write((const char *)&magic_word, sizeof(unsigned int));
      out.write((const char *)&LMsize,     sizeof(int));
      out.write((const char *)&numT,       sizeof(int));
      out.write((const char *)&nmax,       sizeof(int));
      out.write((const char *)&Gsize,      sizeof(int));
      out.write((const char *)&times[0],   sizeof(double)*numT);
      
      for (auto lm : LMset) {
	int LL = lm.first;
	int MM = lm.second;
	out.write((const char *)&LL, sizeof(int));
	out.write((const char *)&MM, sizeof(int));
      }

      int igrp = 0;
      for (auto v : groups) {

	for (unsigned i=0; i<numT; i++) {
	  
	  for (auto lm : LMset) {
	    unsigned LL = lm.first;
	    unsigned MM = lm.second;
	    
	    for (unsigned n=0; n<nmax; n++) {
	      Key c = {LL, MM, n, 0, index}, s = {LL, MM, n, 1, index};
	      auto rc = RC.find(c);
	      auto rs = RC.find(s);
	      double valC = 0.0, valS = 0.0;
	      if (rc != RC.end()) {
		for (auto u : v) {
		  if (u >-1 and u < ncomp) valC += RC[c](i, u);
		}
	      }
	      if (rs != RC.end() and MM>0) {
		for (auto u : v) {
		  if (u >-1 and u < ncomp) valS += RC[s](i, u);
		}
	      }
	      
	      // Retrend
	      //
	      if (type == TrendType::totPow) {
		valC = valC*totPow;
		valS = valS*totPow;
		if (useMean) {
		  valC += mean[c];
		  valS += mean[s];
		}
	      } else if (type == TrendType::totVar) {
		valC = valC*totVar + mean[c];
		valS = valS*totVar + mean[s];
	      } else {
		valC = valC*var[c] + mean[c];
		valS = valS*var[s] + mean[s];
	      }
	      
	      out.write((const char *)&valC, sizeof(double));
	      out.write((const char *)&valS, sizeof(double));
	    }
	    // radial order loop
	  }
	  // L, M loop
	}
	// T loop
	igrp++;
      }
      // Group loop
    }
    // Output file okay
    else {
      std::cout << "Could not open <" << fname.str() << ">" << std::endl;
    }

    out.close();		// Close file
    fname.str("");		// Reset the strstream

    fname << runtag << "." << name << ".recon_cmpl";
    out.open(fname.str());	// Open next file

    if (out) {

      unsigned LMsize = LMset.size();
      int       Gsize = 1;
      int        numT = times.size();
      int        nmax = this->data.begin()->second->header.nmax;

      out.write((const char *)&magic_word, sizeof(unsigned int));
      out.write((const char *)&LMsize,     sizeof(unsigned int));
      out.write((const char *)&numT,       sizeof(int));
      out.write((const char *)&nmax,       sizeof(int));
      out.write((const char *)&Gsize,      sizeof(int));
      out.write((const char *)&times[0],   sizeof(double)*numT);
      
      for (auto lm : LMset) {
	int LL = lm.first;
	int MM = lm.second;
	out.write((const char *)&LL, sizeof(int));
	out.write((const char *)&MM, sizeof(int));
      }
      
      std::set<int> comple;
      for (int n=0; n<ncomp; n++) {
	bool found = false;
	for (auto s : groups) {
	  if (s.find(n) != s.end()) found = true;
	}
	if (!found) comple.insert(n);
      }
      
      for (unsigned i=0; i<numT; i++) {

	for (auto lm : LMset) {
	  unsigned LL = lm.first;
	  unsigned MM = lm.second;
	  
	  for (unsigned n=0; n<nmax; n++) {
	    Key c = {LL, MM, n, 0, index}, s = {LL, MM, n, 1, index};
	    auto rc = RC.find(c);
	    auto rs = RC.find(s);
	    double valC = 0.0, valS = 0.0;
	    if (rc != RC.end()) {
	      for (auto u : comple) {
		if (u >-1 and u < ncomp) valC += RC[c](i, u);
	      }
	    }
	    if (rs != RC.end() and MM>0) {
	      for (auto u : comple) {
		if (u >-1 and u < ncomp) valS += RC[s](i, u);
	      }
	    }
	    
	    // Retrend
	    //
	    if (type == TrendType::totPow) {
	      valC = valC*totPow;
	      valS = valS*totPow;
	      if (useMean) {
		valC += mean[c];
		valS += mean[s];
	      }
	    } else if (type == TrendType::totVar) {
	      valC = valC*totVar + mean[c];
	      valS = valS*totVar + mean[s];
	    } else {
	      valC = valC*var[c] + mean[c];
	      valS = valS*var[s] + mean[s];
	    }

	    out.write((const char *)&valC, sizeof(double));
	    out.write((const char *)&valS, sizeof(double));
	  }
	  // radial order loop
	}
	// L, M loop
      }
      // T loop
    }
    // Open file
    else {
      std::cout << "Could not open <" << fname.str() << ">" << std::endl;
    }
  }

  //! Constructor
  TableDB::TableDB(const std::string& name, YAML::Node& node, unsigned index,
		   unsigned stride, double tmin, double tmax) :
    CoefDB(name, node, index)
  {
    std::string file;
    
    if (node["filename"]) {
      file = node["filename"].as<std::string>();
    } else {
      std::cout << "TableDB: no data table filename found" << std::endl;
      exit(-4);
    }

    read(file, stride, tmin, tmax);
  }


  template<typename T>
  void pop_front(std::vector<T>& vec)
  {
    assert(!vec.empty());
    vec.erase(vec.begin());
  }

  void TableDB::read(const std::string& file, unsigned stride,
		     double tmin, double tmax)
  {
    std::map<std::string, int> ret;
    std::ifstream in(file);
    
    // Check for good opened file
    if (!in) {
      std::cerr << "TableDB: could not open data table file <"
		<< file << ">" << std::endl;
    exit(-2);
    }

    unsigned counter = 0;

    while (in.good()) {
      std::string line;
      std::getline(in, line);

      // Good line?
      //
      if (in.good()) {
	std::istringstream sin(line);
	std::vector<double> row;
	try {
	  double val;
	  while (1) { sin >> val; row.push_back(val); }
	}
	catch (std::istream::failure e) {
	  if (row.size()==0) {
	    std::cerr << "Empty row?" << std::endl;
	    break;
	  }
	}

	double now  = row[0];
	if (now >= tmin and now <= tmax) {
	  if (counter++ % stride == 0) {
	    pop_front(row);
	    data[now] = row;
	  }
	}
      }
    }

    ntimes  = data.size();
    nfields = data.begin()->second.size();

    for (auto v : data) {
      times.push_back(v.second[0]);
    }

    for (unsigned f=0; f<nfields; f++) {
      std::vector<double> D(ntimes);
      for (int n=0; n<ntimes; n++) D[n] = data[D[n]][f];
      coefs[{f}] = D;
    }

  }


  CoefContainer::CoefContainer(const std::string& spec)
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
	  comps.push_back(std::make_shared<SphDB>(name, node, comps.size(), stride, tmin, tmax));
	} else if (ctype.find("cylinder")==0) {
	  comps.push_back(std::make_shared<CylDB>(name, node, comps.size(), stride, tmin, tmax));
	} else if (ctype.find("table")==0) {
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


  void CylDB::write_ascii(int ncomp, const TrendType& type, bool useMean,
			  double totPow, double totVar,
			  std::map<Key, double>& mean,
			  std::map<Key, double>& var,
			  std::map<Key, Eigen::MatrixXd>& RC,
			  const std::vector<std::set<int>>& groups,
			  const std::vector<double>& times,
			  const std::string& runtag)
  {
    std::string filename = runtag + "." + name + "_tot.coefs";
    std::ofstream out(filename);

    int NMAX = this->data.begin()->second->nmax;

    std::set<unsigned> Mset;
    for (auto k : keys) Mset.insert(k[0]);

    if (out) {

      for (int i=0; i<times.size(); i++) {
	out << std::setw(15) << std::setprecision(6) << times[i] << std::endl;
	
	for (auto M : Mset) {

	  out << std::setw(15) << M << std::setw(15) << NMAX << std::endl;

	  for (unsigned n=0; n<NMAX; n++) {
	    Key c = {M, n, 0, index};
	    auto rc = RC.find(c);
	    if (rc == RC.end())
	      out << std::setw(15) << std::setprecision(6) << 0.0;
	    else {
	      double acc = 0.0;
	      for (int j=0; j<ncomp; j++) acc += RC[c](i, j);
	      out << std::setw(15) << std::setprecision(6)
		  << acc*var[c] + mean[c];
	    }
	  }
	  out << std::endl;

	  if (M) {
	    for (unsigned n=0; n<NMAX; n++) {
	      Key s = {M, n, 1, index};
	      auto rs = RC.find(s);
	      if (rs == RC.end())
		out << std::setw(15) << std::setprecision(6) << 0.0;
	      else {
		double acc = 0.0;
		for (int j=0; j<ncomp; j++) acc += RC[s](i, j);
		out << std::setw(15) << std::setprecision(6)
		    << acc*var[s] + mean[s];
	      }
	    }
	    out << std::endl;
	  }
	}
	// END: M loop
      }
      // END: time loop

    } else {
      std::cout << "Could not open <" << filename << ">" << std::endl;
      exit(-1);
    }

    for (int j=0; j<groups.size(); j++) {

      std::ostringstream filename;
      filename << runtag << "." << name << "_" << j << ".coefs";
      std::ofstream out(filename.str());

      if (out) {
	for (int i=0; i<times.size(); i++) {
	  for (auto M : Mset) {

	    out << std::setw(15) << std::setprecision(6) << times[i] << std::endl;
	    out << std::setw(15) << M << std::setw(15) << NMAX << std::endl;

	    for (unsigned n=0; n<NMAX; n++) {
	      Key c {static_cast<unsigned>(M), n, 0, index};
	      auto rc = RC.find(c);
	      if (rc == RC.end())
		out << std::setw(15) << std::setprecision(6) << 0.0;
	      else {
		double val = 0.0;
		for (auto v : groups[j]) val += RC[c](i, v);
		out << std::setw(15) << std::setprecision(6)
		    << val*var[c] + mean[c];
	      }
	    }
	    out << std::endl;

	    if (M) {
	      for (unsigned n=0; n<NMAX; n++) {
		Key s = {static_cast<unsigned>(M), n, 1, index};
		auto rs = RC.find(s);
		if (rs == RC.end())
		  out << std::setw(15) << std::setprecision(6) << 0.0;
		else {
		  double val = 0.0;
		  for (auto v : groups[j]) val += RC[s](i, v);
		  out << std::setw(15) << std::setprecision(6)
		      << val*var[s] + mean[s];
		}
	      }
	      out << std::endl;
	    }
	  }
	}
	out.close();
      } else {
	std::cout << "Could not open <" << filename.str() << ">" << std::endl;
	exit(-1);
      } // file open failure

    } // numW loop

  } // coefs


  void SphDB::write_ascii(int ncomp, const TrendType& type, bool useMean,
			  double totPow, double totVar,
			  std::map<Key, double>& mean,
			  std::map<Key, double>& var,
			  std::map<Key, Eigen::MatrixXd>& RC,
			  const std::vector<std::set<int>>& groups,
			  const std::vector<double>& times,
			  const std::string& runtag)
  {
    std::set<std::pair<unsigned, unsigned>> LMset;
    for (auto k : keys) LMset.insert({k[0], k[1]});

    for (auto lm : LMset) {
      unsigned   LL = lm.first;
      unsigned   MM = lm.second;
      unsigned NMAX = this->data.begin()->second->header.nmax;

      std::ostringstream filename;
      filename << runtag << "." << name << "_tot_" << LL << "_" << MM << ".coefs";
      std::ofstream out(filename.str());

      int q = 1;
      if (MM) q = 2;

      if (out) {
	out << "# L=" << LL << " M=" << MM << std::endl;

	for (unsigned i=0; i<times.size(); i++) {
	  out << std::setw(15) << std::setprecision(6) << times[i] << std::endl;
	  out << std::setw(15) << NMAX << std::endl;

	  for (unsigned n=0; n<NMAX; n++) {
	    Key c = {LL, MM, n, 0, index};
	    auto rc = RC.find(c);
	    if (rc == RC.end())
	      out << std::setw(15) << std::setprecision(6) << 0.0;
	    else {
	      double acc = 0.0;
	      for (int j=0; j<ncomp; j++) acc += RC[c](i, j);
	      out << std::setw(15) << std::setprecision(6) << acc;
	    }
	  }
	  out << std::endl;

	  if (MM) {
	    for (unsigned n=0; n<NMAX; n++) {
	      Key s = {LL, MM, n, 1, index};
	      auto rs = RC.find(s);
	      if (rs == RC.end())
		out << std::setw(15) << std::setprecision(6) << 0.0;
	      else {
		double acc = 0.0;
		for (int j=0; j<ncomp; j++) acc += RC[s](i, j);
		out << std::setw(15) << std::setprecision(6) << acc;
	      }
	    }
	    out << std::endl;
	  }
	}
	out.close();
      } else {
	std::cout << "Could not open <" << filename.str() << ">" << std::endl;
	exit(-1);
      }
    }

    for (int j=0; j<ncomp; j++) {

      for (auto lm : LMset) {

	unsigned   LL = lm.first;
	unsigned   MM = lm.second;
	unsigned NMAX = this->data.begin()->second->header.nmax;

	std::ostringstream filename;
	filename << runtag << "." << name << "_" << LL << "_" << MM << "_" << j << ".coefs";
	std::ofstream out(filename.str());

	if (out) {
	  for (unsigned i=0; i<times.size(); i++) {
	    out << std::setw(15) << std::setprecision(6) << times[i] << std::endl;
	    out << std::setw(15) << NMAX << std::endl;

	    for (unsigned n=0; n<NMAX; n++) {
	      Key c = {LL, MM, n, 0, index};
	      auto rc = RC.find(c);
	      if (rc == RC.end())
		out << std::setw(15) << std::setprecision(6) << 0.0;
	      else {
		double val = RC[c](i, j);
		out << std::setw(15) << std::setprecision(6) << val;
	      }
	    }
	    out << std::endl;

	    if (MM) {
	      for (unsigned n=0; n<NMAX; n++) {
		Key s = {LL, MM, n, 1, index};
		auto rs = RC.find(s);
		if (rs == RC.end())
		  out << std::setw(15) << std::setprecision(6) << 0.0;
		else {
		  double val = RC[s](i, j);
		  out << std::setw(15) << std::setprecision(6) << val;
		}
	      }
	      out << std::endl;
	    }
	  }
	  out.close();
	} else {
	  std::cout << "Could not open <" << filename.str() << ">" << std::endl;
	  exit(-1);
	} // file open failure

      } // LM pairs

    } // ncomp loop

  } // coefs

  void TableDB::write_ascii(int ncomp, const TrendType& type, bool useMean,
			    double totPow, double totVar,
			    std::map<Key, double>& mean,
			    std::map<Key, double>& var,
			    std::map<Key, Eigen::MatrixXd>& RC,
			    const std::vector<std::set<int>>& groups,
			    const std::vector<double>& times,
			    const std::string& runtag)
  {
    std::string filename = runtag + "." + name + ".data";
    std::ofstream out(filename);

    std::set<unsigned> Fset;
    for (auto k : keys) Fset.insert(k[0]);

    if (out) {

      for (int i=0; i<times.size(); i++) {
	out << std::setw(15) << std::setprecision(6) << times[i] << std::endl;
	
	for (auto F : Fset) {

	  for (unsigned n=0; n<nfields; n++) {
	    Key k = {F, index};
	    auto rc = RC.find(k);
	    if (rc == RC.end())
	      out << std::setw(15) << std::setprecision(6) << 0.0;
	    else {
	      double acc = 0.0;
	      for (int j=0; j<ncomp; j++) acc += RC[k](i, j);
	      out << std::setw(15) << std::setprecision(6)
		  << acc*var[k] + mean[k];
	    }
	  }
	  out << std::endl;
	}
	// END: F loop
      }
      // END: time loop

    } else {
      std::cout << "Could not open <" << filename << ">" << std::endl;
      exit(-1);
    }

    for (int j=0; j<groups.size(); j++) {

      std::ostringstream filename;
      filename << runtag << "." << name << "_" << j << ".data";
      std::ofstream out(filename.str());

      if (out) {
	for (int i=0; i<times.size(); i++) {
	  for (auto F : Fset) {

	    out << std::setw(15) << std::setprecision(6) << times[i] << std::endl;
	    out << std::setw(15) << F << std::setw(15) << nfields << std::endl;

	    for (unsigned n=0; n<nfields; n++) {
	      Key k = {static_cast<unsigned>(F), index};
	      auto rc = RC.find(k);
	      if (rc == RC.end())
		out << std::setw(15) << std::setprecision(6) << 0.0;
	      else {
		double val = 0.0;
		for (auto v : groups[j]) val += RC[k](i, v);
		out << std::setw(15) << std::setprecision(6)
		    << val*var[k] + mean[k];
	      }
	    }
	    out << std::endl;
	  }
	}
	out.close();
      } else {
	std::cout << "Could not open <" << filename.str() << ">" << std::endl;
	exit(-1);
      } // file open failure

    } // numW loop

  } // coefs

}

