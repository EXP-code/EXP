
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
    in.read((char *)&cos_c[mm][0], sizeof(double)*nmax);
    if (not in) return false;

    if (mm) {
      sin_c[mm].resize(nmax);
      in.read((char *)&sin_c[mm][0], sizeof(double)*nmax);
      if (not in) return false;
    }
  }

  return true;
}


CylindricalCoefs::CylindricalCoefs(const std::string& file, unsigned stride)
{
  std::ifstream in(file);

  unsigned counter = 0;
  
  while (in.good()) {
    CoefPtr c = std::make_shared<Coefs>();
    if (not c->read(in)) break;
    if (counter++ % stride == 0) data[c->time] = c;
  }

  mmax   = data.begin()->second->mmax;
  nmax   = data.begin()->second->nmax;
  ntimes = data.size();

  for (auto v : data) {
    times.push_back(v.second->time);
  }
}

CylindricalCoefs::D2pair CylindricalCoefs::interpolate(const double time)
{
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

  double A = (*hi - time)/(*hi - *lo);
  double B = (time - *lo)/(*hi - *lo);

  int iA = std::distance(times.begin(), lo);
  int iB = std::distance(times.begin(), hi);

  auto cA = data[times[iA]];
  auto cB = data[times[iB]];

  D2pair ret;
  ret.first .resize(mmax+1);
  ret.second.resize(mmax+1);
  for (int m=0; m<mmax+1; m++) {
    ret.first [m].resize(nmax, 0.0);
    ret.second[m].resize(nmax, 0.0);
  }

  for (int m=0; m<mmax+1; m++) {
    for (int n=0; n<nmax; n++)
      ret.first[m][n] = A*cA->cos_c[m][n] + B*cB->cos_c[m][n];
    
    if (m) {
      for (int n=0; n<nmax; n++)
	ret.second[m][n] = A*cA->sin_c[m][n] + B*cB->sin_c[m][n];
    }    
  }

  return ret;
}

CylindricalCoefs::D2pair CylindricalCoefs::operator()(const double time)
{
  D2pair ret;
  ret.first .resize(mmax+1);
  ret.second.resize(mmax+1);
  for (int m= 0; m<mmax+1; m++) {
    ret.first [m].resize(nmax, 0.0);
    ret.second[m].resize(nmax, 0.0);
  }

  auto it = data.find(time);
  if (it == data.end()) return ret;

  for (int m=0; m<mmax+1; m++) {
    ret.first[m] = it->second->cos_c[m];
    if (m) ret.second[m] = it->second->sin_c[m];
  }

  return ret;
}

