// simanneal.c++ Implementation of a General Purpose Simulated
// Annealing Class
//
// Uses Cauchy training

#include <SimAnn.H>

SimAnn::SimAnn (SimAnn::Func1d f, const int d):
  func (f), ddwell (20), rrange (M_PI_2), t0 (0.0), K (1.0),
  rho (0.5), dt (0.1), tscale (0.1), maxit (400), c_jump (100.0), fsave(0)
{
  number_range = std::uniform_real_distribution<double>(-rrange, rrange);
  number_01    = std::uniform_real_distribution<double>(0.0, 1.0);

  y = ybest = std::numeric_limits<double>::max();
}

void SimAnn::set_up (Func1d f, const int seed)
{
  func = f;
  gen  = std::mt19937(seed); 
  y = ybest = std::numeric_limits<double>::max();
}

// increase the temperature until the system "melts"
double SimAnn::melt (const int iters)
{
  int n = iters;
  if (n < 1) n = maxit;

  bool ok = false;
  double cold=0, c = 0.0;
  double t = t0;

  for (int i=0; i<n; i++) {
    if (i > 0 && c > 0.0) {
      cold = c;
      ok = true;
    }

    t += dt;

    for (auto & v : x) 
      v  += rho * t * tan (number_range(gen));

    equilibrate (t, ddwell);

    // "energy"
    //
    double ynew = func(x);
    c = ynew - y;

    if (c < 0.0 && ynew < ybest) {
      xbest = x;
      ybest = ynew;
    }

    y = ynew;

    // phase transition
    //
    if (ok && c > (c_jump * cold)) break;
  }

  return t0 = t;
}

// iterate a few times at the present temperature to get to thermal
// equilibrium
//
int SimAnn::equilibrate (const double t, const int n)
{
  int equil = 0, i;

  for (i=0; i<n; i++) {

    xnew = x;
    for (auto & v : xnew) v += rho * t * tan (number_range(gen));

    //
    // "energy"
    //
    double ynew = func(xnew);
    double c    = ynew - y;

    if (c < 0.0) {		// keep xnew if energy is reduced
      auto & xtmp = x;
      x    = xnew;
      xnew = xtmp;
      y    = ynew;

      if (y < ybest) {
	xbest = x;
	ybest = y;
      }

      double delta = fabs (c);
      delta = (y != 0.0) ? delta / y : (ynew != 0.0) ? delta / ynew : delta;

      // equilibrium is defined as a 10% or smaller change in 10
      // iterations
      if (delta < 0.10) equil++;
      else              equil = 0;

    }
    else
      {
	/* keep xnew with probability, p, if ynew is increased */
	/*
	  p = exp( - (ynew - y) / (K * t) );
	  
	  if ( p > (*number_01)() )
	  {
	  xtmp = x;
	  x = xnew;
	  xnew = xtmp;
	  y = ynew;
	  equil = 0;
	  }
	  else
	*/
	
	equil++;
      }

    if (equil > 9) break;
  }

  return i + 1;
}

// cool the system with annealing
double SimAnn::anneal (const int iters)
{
  double t = 0.0;
  int n = iters;
  if (n < 1) n = maxit;

  equilibrate (t0, 10.0 * ddwell);

  for (int i=0; i<n; i++) {

    double t = t0 / (1.0 + i * tscale);

    equilibrate (t, ddwell);

    xnew = x;
    for (auto & v : xnew) v += rho * t * tan(number_range(gen));

    // "energy"
    //
    double ynew = func(xnew);

    if (ynew <= y) {		// keep xnew if energy is reduced
      auto & xtmp = x;
      x    = xnew;
      xnew = xtmp;
      y    = ynew;
      
      if (y < ybest) {
	xbest = x;
	ybest = y;
      }
      continue;
    }
    else {
      // keep xnew with probability, p, if ynew is increased
      double p = exp (-(ynew - y) / (K * t));

      if (p > number_01(gen)) {
	auto & xtmp = x;
	x    = xnew;
	xnew = xtmp;
	y    = ynew;
      }
      
      if (fsave) log_state(i);
    }
  }

  equilibrate (t, 10 * ddwell);

  t0 = t;

  return t;
}

void SimAnn::log_state(int i)
{
  std::vector<double> xp;
  if (pmap) xp = pmap(x);
  else      xp = x;

  if (fsave==2) {
    *str << std::setw( 5) << i
	 << std::setw(15) << ybest
	 << std::setw(15) << y;
    for (auto v : xp) *str << std::setw(15) << v;
    *str << std::endl;
  } else {
    std::ofstream out(fname.c_str(), std::ios::app);
    out << std::setw( 5) << i
	<< std::setw(15) << ybest
	<< std::setw(15) << y;
    for (auto v : xp) out << std::setw(15) << v;
    out << std::endl;
  }
}
