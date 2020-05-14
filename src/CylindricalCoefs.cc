
#include <cmath>
#include "CylindricalCoefs.H"

bool CylindricalCoefs::Coefs::read(std::istream& in)
{
  CylCoefHeader header;
  in.read((char *)&header, sizeof(CylCoefHeader));
  if (not in) return false;

  time = header.time;
  nmax = header.nmax;
  mmax = header.mmax;

  cos_c.resize(mmax+1);
  sin_c.resize(mmax+1);
  
  for (int mm=0; mm<=mmax; mm++) {

    cos_c[mm].resize(nmax);
    in.read((char *)cos_c[mm].data(), sizeof(double)*nmax);
    if (not in) return false;

    if (mm) {
      sin_c[mm].resize(nmax);
      in.read((char *)sin_c[mm].data(), sizeof(double)*nmax);
      if (not in) return false;
    }
  }

  return true;
}


CylindricalCoefs::CylindricalCoefs(const std::string& file, int nd, unsigned stride)
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
    data[T] = v;
  }
    
  mmax   = data.begin()->second->mmax;
  nmax   = data.begin()->second->nmax;
  ntimes = data.size();

  // Create and initialize cached return buffer
  //
  ret = std::make_shared<Dpair>();

  std::get<0>(*ret).resize(mmax+1);
  std::get<1>(*ret).resize(mmax+1);

  for (auto & v : std::get<0>(*ret)) {
    v.resize(nmax);
    std::fill(v.begin(), v.end(), 0.0);
  }

  for (auto & v : std::get<1>(*ret)) {
    v.resize(nmax);
    std::fill(v.begin(), v.end(), 0.0);
  }

}

CylindricalCoefs::DataPtr CylindricalCoefs::interpolate(const double T)
{
  double time = to_ndigits(T);

  if (time < times.front() or time > times.back()) {
    std::cerr << "Time=" << time << " is offgrid [" << times.front()
	      << ", " << times.back() << "]" << std::endl;
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

  auto cA = data[times[iA]];
  auto cB = data[times[iB]];

  for (int m=0; m<mmax+1; m++) {
    for (int n=0; n<nmax; n++)
      std::get<0>(*ret)[m][n] = A*cA->cos_c[m][n] + B*cB->cos_c[m][n];
    
    if (m) {
      for (int n=0; n<nmax; n++)
	std::get<1>(*ret)[m][n] = A*cA->sin_c[m][n] + B*cB->sin_c[m][n];
    }    
  }

  return ret;
}

CylindricalCoefs::DataPtr CylindricalCoefs::operator()(const double time)
{
  auto it = data.find(time);
  if (it == data.end()) return ret;

  for (int m=0; m<mmax+1; m++) {
    std::get<0>(*ret)[m] = it->second->cos_c[m];
    if (m) std::get<1>(*ret)[m] = it->second->sin_c[m];
  }

  return ret;
}

