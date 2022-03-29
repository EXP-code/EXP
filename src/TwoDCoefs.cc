#include <localmpi.H>

#include <yaml-cpp/yaml.h>	 // YAML support
#include <cstring>		 // For strncpy
#include <memory>		 // Shared pointers

#include "global.H"
#include "TwoDCoefs.H"

bool TwoDCoefs::Coefs::read(std::istream& in)
{
  in.exceptions ( std::istream::failbit | std::istream::badbit );

  // Save initial stream position
  //
  auto curpos = in.tellg();

  // Coefficient magic number
  //
  const unsigned int cmagic = 0xc0a57a4;

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
      time  = node["time" ].as<double>();
      scale = node["scale"].as<double>();
      rmax  = node["rmax" ].as<double>();
      nmax  = node["nmax" ].as<int>();
      Lmax  = node["Lmax" ].as<int>();
      id    = node["id"   ].as<int>();
      std::string ID = node["id"].as<std::string>();
      
      // Get data
      //
      data.resize(2*(Lmax+1), nmax);
      in.read((char *)data.data(), data.size()*sizeof(double));

    } else {
      std::ostringstream sout;
      sout << "TwoDCoefs: magic number mismatch.  Got <"
	   << std::hex << tmagic << "> but expected <"
	   << std::hex << cmagic << ">";
      throw std::runtime_error(sout.str());
    }
  }
  catch (std::istream::failure e) {
    if (not in.eof())
      std::cerr << "Exception reading coefficient file: "
		<< e.what() << std::endl;
    return false;
  }
  
  return true;
}


TwoDCoefs::TwoDCoefs(const std::string& file, unsigned stride)
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
      Lmax = v->Lmax;
      nmax = v->nmax;
    } else {

      if (Lmax != v->Lmax) {
	std::cout << "TwoDCoefs: [" << myid << "] "
		  << "coefficient stanza rank mismatch: lmax=" << v->Lmax
		  << ", expected " << Lmax << std::endl;
	MPI_Finalize();
	exit(-31);
      }

      if (nmax != v->nmax) {
	std::cout << "TwoDCoefs: [" << myid << "] "
		  << "coefficient stanza rank mismatch: nmax=" << v->nmax
		  << ", expected " << nmax << std::endl;
	MPI_Finalize();
	exit(-32);
      }
    } 
    
    coefs[T] = v->data;
  }

  ntimes = times.size();

  ret = std::make_shared<Eigen::MatrixXd>(2*(Lmax+1), nmax);
}


TwoDCoefs::DataPtr TwoDCoefs::interpolate(const double T)
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

  Eigen::MatrixXd & cA = coefs[times[iA]];
  Eigen::MatrixXd & cB = coefs[times[iB]];

  for (int lm=0; lm<2*(Lmax+1); lm++) {
    for (int n=0; n<nmax; n++) (*ret)(lm, n) = A*cA(lm, n) + B*cB(lm, n);
  }

  return ret;
}

TwoDCoefs::DataPtr TwoDCoefs::zero_ret()
{
  std::fill(ret->data(), ret->data() + ret->size(), 0.0);
  return ret;
}

TwoDCoefs::DataPtr TwoDCoefs::operator()(const double T)
{
  double time = to_ndigits(T);

  auto it = coefs.find(time);
  if (it == coefs.end()) return zero_ret();

  *ret = it->second;

  return ret;
}
