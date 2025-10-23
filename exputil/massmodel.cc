#include <cstdlib>
#include <cfloat>
#include <cmath>

#include <iostream>
#include <iomanip>
#include <fstream>
#include <sstream>
#include <string>
#include <vector>
#include <limits>

#include "massmodel.H"
#include "interp.H"

int wordSplit(string& x, vector<string>& words);

int SphericalModelTable::count    = 0;
int SphericalModelTable::even     = 0;
int SphericalModelTable::linear   = 1;
int SphericalModelTable::chebyN   = 0;

SphericalModelTable::SphericalModelTable
(string filename, int DIVERGE, double DIVERGE_RFAC, int EXTERNAL)
{
  count++;

  dim = 3;
  ModelID = "SphericalModelTable(" + filename + ")";
  dist_defined = false;

  /**************************************************************************
   *
   *	Routine to read in model file
   *
   *	1) r  2) rho  3) M(r)  4) U(r)
   *
   *	First two lines are for header/printed out
   *    and the third line is the number of entries
   *
   *  Flags:
   *  -----
   *  EXTERNAL        changes boundary conditions on spline for use with
   *                  model + external potential
   *
   **************************************************************************/

  const int MAXLINE = 256;
  char cbuf[MAXLINE];
  string x;
  double radius;

  ifstream from(filename.c_str());
  if (!from) {
    string errmsg = "cannot open input file " + filename;
    bomb(errmsg.c_str());
  }

  // Any line that starts with a "!" or a "#" is ignored

  do {
    from.getline(cbuf, MAXLINE);
    x = cbuf;
  } while (x.find('#') < x.size() || x.find('!') < x.size());

  istringstream ist(x, istringstream::in);
  ist >> num;

  density.x. resize(num);
  density.y. resize(num);
  density.y2.resize(num);
  mass.   x. resize(num);
  mass.   y. resize(num);
  mass.   y2.resize(num);
  pot.    x. resize(num);
  pot.    y. resize(num);
  pot.    y2.resize(num);

  density.num = num;
  mass.num    = num;
  pot.num     = num;
  
  for (int i=0; i<num; i++) {
    from.getline(cbuf, MAXLINE);
    x = cbuf;
    istringstream ist(x, istringstream::in);
    ist >> radius;
    ist >> density.y[i];
    ist >> mass.y[i];
    ist >> pot.y[i];
    density.x[i] = radius;
    mass.x[i] = radius;
    pot.x[i] = radius;
    if (DIVERGE)
      density.y[i] *= pow(radius+std::numeric_limits<double>::min(), DIVERGE_RFAC);
  }
  
  Spline(density.x, density.y, 0.0, 0.0, density.y2);
  Spline(mass.x, mass.y, 0.0, 0.0, mass.y2);
  if (EXTERNAL)
    Spline(pot.x, pot.y, 0.0, mass.y[mass.num-1]/(radius*radius),pot.y2);
  else
    Spline(pot.x, pot.y, -1.0e30, -1.0e30, pot.y2);
  
  num_params = 0;

  if (from.getline(cbuf, MAXLINE)) {
    x = cbuf;

    vector<string> words;
    num_params = wordSplit(x, words);
    if (num_params>0) {
      params.resize(num_params);
      for (int i=0; i<num_params; i++) params[i] = atof(words[i].c_str());
    }
  }
  from.close();

  diverge = DIVERGE;
  diverge_rfac = DIVERGE_RFAC;
  external = EXTERNAL;

}

SphericalModelTable::SphericalModelTable
(int NUM, 
 double *r, double *d, double *m, double *p,
 int DIVERGE, double DIVERGE_RFAC, int EXTERNAL,
 string ID)
{
  count++;

  dim = 3;
  ModelID = "SphericalModelTable(" + ID + ")";
  dist_defined = false;
  num = NUM;

  double radius = 0.0;
  
  density.x. resize(num);
  density.y. resize(num);
  density.y2.resize(num);

  mass.x. resize(num);
  mass.y. resize(num);
  mass.y2.resize(num);

  pot.x. resize(num);
  pot.y. resize(num);
  pot.y2.resize(num);

  density.num = num;
  mass.num    = num;
  pot.num     = num;
  
  for (int i=0; i<num; i++) {
    radius = r[i];
    density.y[i] = d[i];
    mass.y[i] = m[i];
    pot.y[i] = p[i];

    density.x[i] = radius;
    mass.x[i] = radius;
    pot.x[i] = radius;

    if (DIVERGE)
      density.y[i] *= pow(radius+std::numeric_limits<double>::min(), DIVERGE_RFAC);
  }
  
  Spline(density.x, density.y, 0.0, 0.0, density.y2);
  Spline(mass.x, mass.y, 0.0, 0.0, mass.y2);
  if (EXTERNAL)
    Spline(pot.x, pot.y, 0.0, mass.y[mass.num-1]/(radius*radius),pot.y2);
  else
    Spline(pot.x, pot.y, -1.0e30, -1.0e30, pot.y2);
  
  num_params = 0;

  diverge = DIVERGE;
  diverge_rfac = DIVERGE_RFAC;
  external = EXTERNAL;

}

SphericalModelTable::SphericalModelTable
(std::vector<double>& r, std::vector<double>& d,
 std::vector<double>& m, std::vector<double>& p,
 int DIVERGE, double DIVERGE_RFAC, int EXTERNAL,
 string ID)
{
  count++;

  dim = 3;
  ModelID = "SphericalModelTable(" + ID + ")";
  dist_defined = false;
  num = r.size();

  double radius = 0.0;
  
  density.x. resize(num);
  density.y. resize(num);
  density.y2.resize(num);

  mass.x. resize(num);
  mass.y. resize(num);
  mass.y2.resize(num);

  pot.x. resize(num);
  pot.y. resize(num);
  pot.y2.resize(num);

  density.num = num;
  mass.num    = num;
  pot.num     = num;
  
  for (int i=0; i<num; i++) {
    radius = r[i];
    density.y[i] = d[i];
    mass.y[i] = m[i];
    pot.y[i] = p[i];

    density.x[i] = radius;
    mass.x[i] = radius;
    pot.x[i] = radius;

    if (DIVERGE)
      density.y[i] *= pow(radius+std::numeric_limits<double>::min(), DIVERGE_RFAC);
  }
  
  Spline(density.x, density.y, 0.0, 0.0, density.y2);
  Spline(mass.x, mass.y, 0.0, 0.0, mass.y2);
  if (EXTERNAL)
    Spline(pot.x, pot.y, 0.0, mass.y[mass.num-1]/(radius*radius),pot.y2);
  else
    Spline(pot.x, pot.y, -1.0e30, -1.0e30, pot.y2);
  
  num_params = 0;

  diverge = DIVERGE;
  diverge_rfac = DIVERGE_RFAC;
  external = EXTERNAL;

}


SphericalModelTable::~SphericalModelTable()
{
  count--;
}


double SphericalModelTable::get_mass(double r)
{
  double ans;

  if (r<mass.x[0]) return mass.y[0];
  if (r>mass.x[mass.num-1]) return mass.y[mass.num-1];
  
  if (linear)
    ans = odd2(r, mass.x, mass.y, even);
  else
    Splint1(mass.x, mass.y, mass.y2, r, ans, even);
  return ans;
}

double SphericalModelTable::get_density(double r)
{
  double ans;

  if (diverge) {

    if (r<density.x[0])
      ans = density.y[0];
    else {
      if (linear)
	ans = odd2(r, density.x, density.y, even);
      else
	Splint1(density.x, density.y, density.y2, r, ans, even);
    }
    return ans*pow(r, -diverge_rfac);
  }

  if (r>density.x[density.num-1]) return density.y[density.num-1];
  
  if (linear)
    ans = odd2(r, density.x, density.y, even);
  else
    Splint1(density.x, density.y, density.y2, r, ans, even);

  return ans;
}

double SphericalModelTable::get_pot(const double r)
{
  double ans;

  if (r<pot.x[0]) {
    if (diverge)
      return pot.y[0] + 
	4.0*M_PI*density.y[0]/(3.0-diverge_rfac)*
	(
	 pow(density.x[0], 2.0-diverge_rfac) -
	 (pow(density.x[0], 3.0-diverge_rfac)/r - pow(r, 2.0-diverge_rfac))
	 )
	-4.0*M_PI*density.y[0]/(2.0-diverge_rfac)*
	(
	 pow(density.x[0], 2.0-diverge_rfac) -
	 pow(r, 2.0-diverge_rfac)
	 );
    else
      return pot.y[0];
  }

  if (r>pot.x[pot.num-1]) return pot.y[pot.num-1]*pot.x[pot.num-1]/r;
  
  if (linear)
    ans = odd2(r, pot.x, pot.y, even);
  else
    Splint1(pot.x, pot.y, pot.y2, r, ans, even);

  return ans;
}

double SphericalModelTable::get_dpot(const double r)
{
  double dum, ans;

  if (r<pot.x[0]) {

    if (diverge)
      ans = 4.0*M_PI*density.y[0]*pow(r, 2.0-diverge_rfac)/(3.0-diverge_rfac);
    else {
      if (linear)
	ans = drv2(pot.x[0], pot.x, pot.y, even);
      else
	Splint2(pot.x, pot.y, pot.y2, pot.x[0], dum, ans, even);
    }
  }
  else if (r>pot.x[pot.num-1]) 
    ans = -pot.y[pot.num-1]*pot.x[pot.num-1]/(r*r);
  else {
    if (linear)
      ans = drv2(r, pot.x, pot.y, even);
    else
      Splint2(pot.x, pot.y, pot.y2, r, dum, ans, even);
  }

  return ans;
}

void SphericalModelTable::get_pot_dpot(const double r, double& ur, double &dur)
{
  if (r<pot.x[0]) {
    if (diverge) {
      ur =  pot.y[0] + 
	4.0*M_PI*density.y[0]/(3.0-diverge_rfac)*
	(
	 pow(density.x[0], 2.0-diverge_rfac) -
	 (pow(density.x[0], 3.0-diverge_rfac)/r - pow(r, 2.0-diverge_rfac))
	 )
	-4.0*M_PI*density.y[0]/(2.0-diverge_rfac)*
	(
	 pow(density.x[0], 2.0-diverge_rfac) -
	 pow(r, 2.0-diverge_rfac)
	 );
      dur = 4.0*M_PI*density.y[0]*pow(r, 1.0-diverge_rfac)/(3.0-diverge_rfac);
    }
    else {
      if (linear) {
	ur = odd2(pot.x[0], pot.x, pot.y, even);
	dur = drv2(pot.x[0], pot.x, pot.y, even);
      }
      else
	Splint2(pot.x, pot.y, pot.y2, pot.x[0], ur, dur, even);
    }
  }
  else if (r>pot.x[pot.num-1]) {
    ur = pot.y[pot.num-1]*pot.x[pot.num-1]/r;
    dur = -pot.y[pot.num-1]*pot.x[pot.num-1]/(r*r);
  }
  else {
    if (linear) {
      ur = odd2(r, pot.x, pot.y, even);
      dur = drv2(r, pot.x, pot.y, even);
    }
    else
      Splint2(pot.x, pot.y, pot.y2, r, ur, dur, even);
  }

}

double SphericalModelTable::get_dpot2(const double r)
{
  double dum, ans;

  if (r<pot.x[0]) {
    if (diverge)
      ans = 4.0*M_PI*density.y[0]*pow(r, -diverge_rfac)*
	(1.0-diverge_rfac)/(3.0-diverge_rfac);
    else
      Splint2(pot.x, pot.y, pot.y2, pot.x[0], dum, ans, even);
  }
  else if (r>pot.x[pot.num-1]) 
    ans = 2.0*pot.y[pot.num-1]*pot.x[pot.num-1]/(r*r*r);
  else
    Splint3(pot.x, pot.y, pot.y2, r, dum, dum, ans, even);

  return ans;
}

void SphericalModelTable::print_df(const std::string& name)
{
  if (!dist_defined) bomb("distribution function not defined");

  ofstream out(name);
  if (!out) {
    cerr << "Couldn't open <" << name << "\n";
    return;
  }

  double g, d, h;

  if (df.num) {
    
    for (int i=0; i<df.Q.size(); i++) {

      if (linear) {
	d = odd2(df.Q[i], df.Q, df.fQ);
	g = odd2(df.Q[i], df.Q, df.ffQ);
      } else {
	Splint1(df.Q, df.fQ, df.fQ2, df.Q[i], d);
	Splint1(df.Q, df.ffQ, df.ffQ2, df.Q[i], g);
      }
      
      h = - d*exp(g - df.Q[i]);

      out << std::setw(16) << df.Q[i]
	  << std::setw(16) << d
	  << std::setw(16) << g
	  << std::setw(16) << h
	  << std::endl;
    }
  }
  else if (dfc.num) {
  
    for (int i=0; i<dfc.Q.size(); i++) {
      
      d = dfc.FF.eval(dfc.Q[i]);
      g = dfc.GG.eval(dfc.Q[i]);
      h = - d*exp(g - dfc.Q[i]);

      out << std::setw(16) << dfc.Q[i]
	  << std::setw(16) << d
	  << std::setw(16) << g
	  << std::setw(16) << h
	  << std::endl;
    }
  }

}


void SphericalModelTable::print_model(const std::string& name)
{
  ofstream out(name);
  if (!out) {
    cerr << "Couldn't open <" << name << "\n";
    return;
  }

  out.setf(ios::left);

  out << "# ModelID=" << ModelID << std::endl;
  out << std::setw(22) << "# Radius" 
      << std::setw(22) << "Density" 
      << std::setw(22) << "Mass" 
      << std::setw(22) << "Potential" 
      << std::endl;

  char c = out.fill('-');
  out << std::setw(22) << "# [1]"
      << std::setw(22) << "| [2]"
      << std::setw(22) << "| [3]"
      << std::setw(22) << "| [4]"
      << std::endl;
  out.fill(c);

  out << density.num << std::endl;
  out << setprecision(12) << scientific;
  for (int i=0; i<density.num; i++) {
    out << std::setw(22) << density.x[i]
	<< std::setw(22) << density.y[i]
	<< std::setw(22) << mass.y[i]
	<< std::setw(22) << pot.y[i]
	<< std::endl;
  }

}

void SphericalModelTable::print_model_eval(const std::string& name, int number)
{
  ofstream out(name);
  if (!out) {
    cerr << "Couldn't open <" << name << "\n";
    return;
  }

  out.setf(ios::left);

  out << "# ModelID=" << ModelID << std::endl;
  out << std::setw(22) << "# Radius" 
      << std::setw(22) << "Density" 
      << std::setw(22) << "Mass" 
      << std::setw(22) << "Potential" 
      << std::setw(22) << "Derivative" 
      << std::setw(22) << "Monopole" 
      << std::endl;

  char c = out.fill('-');
  out << std::setw(22) << "# [1]"
      << std::setw(22) << "| [2]"
      << std::setw(22) << "| [3]"
      << std::setw(22) << "| [4]"
      << std::setw(22) << "| [5]"
      << std::setw(22) << "| [6]"
      << std::endl;
  out.fill(c);

  double rmin = density.x[0];
  double rmax = density.x[density.x.size()-1];

  bool logscale = false;
  if (rmin>0.0) {
    logscale = true;
    rmin = log(rmin);
    rmax = log(rmax);
  }

  double dr   = (rmax - rmin)/(number - 1);

  out << setprecision(12) << scientific;
  for (int i=0; i<number; i++) {
    double r = rmin + dr*i;
    if (logscale) r = exp(r);
    out << std::setw(22) << r
	<< std::setw(22) << get_density(r)
	<< std::setw(22) << get_mass(r)
	<< std::setw(22) << get_pot(r)
	<< std::setw(22) << get_dpot(r)
	<< std::setw(22) << get_mass(r)/(r*r)
	<< std::endl;
  }

}


void EmbeddedDiskModel::verbose_df(void)
{
  if (!dist_defined) bomb("Embedded: <distf> not defined yet . . .");
  df->set_verbose();
}
  
double EmbeddedDiskModel::distf(double E, double L) 
{
  if (!dist_defined) bomb("Embedded: <distf> not defined yet . . .");
  return df->distf(E, L);
}

double EmbeddedDiskModel::dfde(double E, double L)
{
  if (!dist_defined) bomb("Embedded: <distf> not defined yet . . .");
  return df->dfdE(E, L);
}

double EmbeddedDiskModel::dfdl(double E, double L)
{
  if (!dist_defined) bomb("Embedded: <distf> not defined yet . . .");
  return df->dfdL(E, L);
}

double EmbeddedDiskModel::d2fde2(double E, double L)
{
  bomb("Embedded: <d2fde2> not implemented");
  return 0.0;
}

void EmbeddedDiskModel::save_df(string& file)
{
  df->write_state(file);
}


SphericalModelMulti::SphericalModelMulti(AxiSymModPtr Real,
					 AxiSymModPtr Fake) 
{
  real  = Real;
  fake  = Fake;
  noneg = true;			// Will requeue samples with negative
				// DF (mass) ratio if true

  rmin_gen = max<double>(Fake->get_min_radius(), Real->get_min_radius());
  rmin_gen = max<double>(rmin_gen, gen_rmin);
  rmax_gen = min<double>(Fake->get_max_radius(), Real->get_max_radius());

  dim = 3;
}


