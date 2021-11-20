#include <localmpi.H>

#include <yaml-cpp/yaml.h>	 // YAML support
#include <cstring>		 // For strncpy
#include <memory>		 // Shared pointers

#include "global.H"
#include "SphericalCoefs.H"

bool SphericalCoefs::Coefs::read(std::istream& in)
{
  in.exceptions ( std::istream::failbit | std::istream::badbit );

  // Save initial stream position
  //
  auto curpos = in.tellg();

  // Coefficient magic number
  //
  const unsigned int cmagic = 0xc0a57a2;

  // Try to read magic #
  //
  unsigned int tmagic;

  // Normed is false by default
  //
  bool normed = false;

  // Header structure
  //
  SphCoefHeader header;

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
    
      // Create buffer
      //
      auto buf = std::make_unique<char[]>(hsize+1);

      // Read YAML string
      //
      in.read(buf.get(), hsize);
      buf[hsize] = 0;		// Null terminate
      
      YAML::Node node = YAML::Load(buf.get());
      
      // Get parameters
      //
      header.Lmax  = node["lmax"  ].as<int>();
      header.nmax  = node["nmax"  ].as<int>();
      header.tnow  = node["time"  ].as<double>();
      header.scale = node["scale" ].as<double>();
      
      std::fill(header.id, header.id+64, 0);
      std::string ID = node["id"].as<std::string>();
      strncpy(header.id, ID.c_str(), std::min<int>(64, ID.size()));
      
      if (node["normed"]) normed = node["normed"].as<bool>();

      data.resize((header.Lmax+1)*(header.Lmax+1));
      for (auto & v : data) v.resize(header.nmax);

      for (int ir=0; ir<header.nmax; ir++) {
	for (int l=0, loffset=0; l<=header.Lmax; loffset+=(2*l+1), l++) {
	  for (int m=0, moffset=0; m<=l; m++) {
	    if (m==0) {
	      in.read((char *)&data[loffset+moffset+0][ir], sizeof(double));
	      moffset += 1;
	    } else {
	      in.read((char *)&data[loffset+moffset+0][ir], sizeof(double));
	      in.read((char *)&data[loffset+moffset+1][ir], sizeof(double));
	      moffset += 2;
	    }
	  }
	}
      }

    } else {
      
      // Rewind file
      //
      in.clear();
      in.seekg(curpos);
      
      in.read((char *)&header, sizeof(SphCoefHeader));

      if (not in) return false;

      data.resize((header.Lmax+1)*(header.Lmax+1));
      for (auto & v : data) v.resize(header.nmax);
    
      for (int ir=0; ir<header.nmax; ir++) {
	for (int l=0; l<(header.Lmax+1)*(header.Lmax+1); l++)
	  in.read((char *)&data[l][ir], sizeof(double));
      }
    }

    lmax = header.Lmax;
    nmax = header.nmax;
    time = header.tnow;
    rscl = header.scale;
  }
  catch (std::istream::failure e) {
    if (not in.eof())
      std::cerr << "Exception reading coefficient file: "
		<< e.what() << std::endl;
    return false;
  }
  
  // Remove prefactors from _true_ normed coefficients
  //
  if (normed) {
    int k = 0;
    for (int l=0; l<=header.Lmax; l++) {
      for (int m=0; m<=l; m++) {
	double fac = sqrt( (0.5*l+0.25)/M_PI * 
			   exp(lgamma(1.0+l-m) - lgamma(1.0+l+m)) );

	if (m != 0) fac *= M_SQRT2;

	// Cosine terms
	for (int ir=0; ir<header.nmax; ir++) data[k][ir] /= fac;
	k++;

	// Sine terms
	if (m != 0) {
	  for (int ir=0; ir<header.nmax; ir++) data[k][ir] /= fac;
	  k++;
	}
      }
    }
  }

  return true;
}


SphericalCoefs::SphericalCoefs(const std::string& file, unsigned stride)
{
  std::ifstream in(file);

  unsigned counter = 0;

  std::vector<CoefPtr> cvec;

  while (in.good()) {
    CoefPtr c = std::make_shared<Coefs>();
    if (not c->read(in)) break;
    if (counter++ % stride == 0) cvec.push_back(c);
  }

  ndigits = 1 - std::log10(cvec[1]->time - cvec[0]->time);
    
  for (auto v : cvec) {		// Convert to fixed point for safer
				// use of float as map key
    double T = to_ndigits(v->time);
    times.push_back(T);

    // Sanity check
    //
    if (times.size()==1) {
      lmax = v->lmax;
      nmax = v->nmax;
    } else {

      if (lmax != v->lmax) {
	std::cout << "SphericalCoefs: [" << myid << "] "
		  << "coefficient stanza rank mismatch: lmax=" << v->lmax
		  << ", expected " << lmax << std::endl;
	MPI_Finalize();
	exit(-31);
      }

      if (nmax != v->nmax) {
	std::cout << "SphericalCoefs: [" << myid << "] "
		  << "coefficient stanza rank mismatch: nmax=" << v->nmax
		  << ", expected " << nmax << std::endl;
	MPI_Finalize();
	exit(-32);
      }
    } 
    
    coefs[T] = v->data;
  }

  ntimes = times.size();

  ret = std::make_shared<Dvector>((lmax+1)*(lmax+1));
  for (auto & v : *ret) v.resize(nmax);
}


SphericalCoefs::DataPtr SphericalCoefs::interpolate(const double T)
{
  double time = to_ndigits(T);

  if (time < times.front() or time > times.back()) {
    std::cerr << "Time=" << time << " is offgrid [" << times.front()
	      << ", " << times.back() << "]" << std::endl;
    stop_signal = 1;
  }

  auto it = std::lower_bound(times.begin(), times.end(), time);
  auto lo = it, hi = it;

  if (hi == times.end()) {
    hi = times.end() - 1;
    lo = hi - 1;
  } else hi++;

  double A = (*hi - T)/(*hi - *lo);
  double B = (T - *lo)/(*hi - *lo);

  int iA = std::distance(times.begin(), lo);
  int iB = std::distance(times.begin(), hi);

  Dvector & cA = coefs[times[iA]];
  Dvector & cB = coefs[times[iB]];

  for (int lm=0; lm<ret->size(); lm++) {
    for (int n=0; n<nmax; n++) (*ret)[lm][n] = A*cA[lm][n] + B*cB[lm][n];
  }

  return ret;
}

SphericalCoefs::DataPtr SphericalCoefs::zero_ret()
{
  for (auto & v : *ret) {
    std::fill(v.begin(), v.end(), 0.0);
  }

  return ret;
}

SphericalCoefs::DataPtr SphericalCoefs::operator()(const double T)
{
  double time = to_ndigits(T);

  auto it = coefs.find(time);
  if (it == coefs.end()) return zero_ret();

  *ret = it->second;

  return ret;
}
