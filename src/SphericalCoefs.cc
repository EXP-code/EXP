
#include "SphericalCoefs.H"

bool SphericalCoefs::Coefs::read(std::istream& in)
{
  SphCoefHeader header;
  in.read((char *)&header, sizeof(SphCoefHeader));
  if (not in) return false;
  
  time = header.tnow;
  nmax = header.nmax;
  lmax = header.Lmax;
  
  data.resize((lmax+1)*(lmax+1));
  for (auto & v : data) data.resize(nmax);

  for (int n=0; n<nmax; n++) {
    for (int ll=0; ll<(lmax+1)*(lmax+1); ll++) {
      in.read((char *)&data[ll][n], sizeof(double));
      if (not in) return false;
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

    lmax = v->lmax;
    nmax = v->nmax;

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
