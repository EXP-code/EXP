
#include <string>
#include <vector>

extern int wordSplit(string& x, vector<string>& words);

SphericalModelTable::SphericalModelTable(string filename, 
			       int DIVERGE, double DIVERGE_RFAC, int EXTERNAL)
{

  dim = 3;
  ModelID = "SphericalModelTable(" + filename + ")";
  dist_defined = FALSE;

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

  do {
    from.getline(cbuf, MAXLINE);
    x = cbuf;
  } while (x.find('#') < x.size() || x.find('!') < x.size());

  istrstream ist(x.c_str(), x.size());
  ist >> num;

  density.x.setsize(1,num);
  density.y.setsize(1,num);
  density.y2.setsize(1,num);
  mass.x.setsize(1,num);
  mass.y.setsize(1,num);
  mass.y2.setsize(1,num);
  pot.x.setsize(1,num);
  pot.y.setsize(1,num);
  pot.y2.setsize(1,num);
  density.num = num;
  mass.num = num;
  pot.num = num;
  
  for (int i=1; i<=num; i++) {
    from.getline(cbuf, MAXLINE);
    x = cbuf;
    istrstream ist((const char *)x, x.length());
    ist >> radius;
    ist >> density.y[i];
    ist >> mass.y[i];
    ist >> pot.y[i];
    density.x[i] = radius;
    mass.x[i] = radius;
    pot.x[i] = radius;
    if (DIVERGE)
      density.y[i] *= pow(radius+std::numeric_limits<float>::min(),DIVERGE_RFAC);
  }
  
  Spline(density.x, density.y, 0.0, 0.0, density.y2);
  Spline(mass.x, mass.y, 0.0, 0.0, mass.y2);
  if (EXTERNAL)
    Spline(pot.x, pot.y, 0.0, mass.y[mass.num]/(radius*radius),pot.y2);
  else
    Spline(pot.x, pot.y, -1.0e30, -1.0e30, pot.y2);
  
  num_params = 0;

  if (from.getline(cbuf, MAXLINE)) {
    x = cbuf;

    vector<string> words;
    num_params = wordSplit(x, words);
    if (num_params>0) {
      params = new double [num_params];
      for (int i=0; i<num_params; i++) params[i] = atof(words[i].c_str());
    }
  }
  from.close();

  diverge = DIVERGE;
  diverge_rfac = DIVERGE_RFAC;
  external = EXTERNAL;

}

SphericalModelTable::SphericalModelTable(int NUM, 
			 double *r, double *d, double *m, double *p,
			 int DIVERGE, double DIVERGE_RFAC, int EXTERNAL,
			 string ID)
{

  dim = 3;
  ModelID = "SphericalModelTable(" + ID + ")";
  dist_defined = FALSE;
  num = NUM;

  double radius;

  density.x.setsize(1,num);
  density.y.setsize(1,num);
  density.y2.setsize(1,num);
  mass.x.setsize(1,num);
  mass.y.setsize(1,num);
  mass.y2.setsize(1,num);
  pot.x.setsize(1,num);
  pot.y.setsize(1,num);
  pot.y2.setsize(1,num);
  density.num = num;
  mass.num = num;
  pot.num = num;
  
  for (int i=1; i<=num; i++) {
    radius = r[i];
    density.y[i] = d[i];
    mass.y[i] = m[i];
    pot.y[i] = p[i];

    density.x[i] = radius;
    mass.x[i] = radius;
    pot.x[i] = radius;

    if (DIVERGE)
      density.y[i] *= pow(radius+std::numeric_limits<float>::min(),DIVERGE_RFAC);
  }
  
  Spline(density.x, density.y, 0.0, 0.0, density.y2);
  Spline(mass.x, mass.y, 0.0, 0.0, mass.y2);
  if (EXTERNAL)
    Spline(pot.x, pot.y, 0.0, mass.y[mass.num]/(radius*radius),pot.y2);
  else
    Spline(pot.x, pot.y, -1.0e30, -1.0e30, pot.y2);
  
  num_params = 0;

  diverge = DIVERGE;
  diverge_rfac = DIVERGE_RFAC;
  external = EXTERNAL;

}


SphericalModelTable::~SphericalModelTable()
{
  if (num_params) delete [] params;
}


double SphericalModelTable::get_mass(double r)
{
  double ans;

  if (r<mass.x[1]) return mass.y[1];
  if (r>mass.x[mass.num]) return mass.y[mass.num];
  
  Splint1(mass.x, mass.y, mass.y2, r, ans, even);
  return ans;
}

double SphericalModelTable::get_density(double r)
{
  double ans;

  if (r<density.x[1]) return density.y[1];
  if (r>density.x[density.num]) return density.y[density.num];
  
  Splint1(density.x, density.y, density.y2, r, ans, even);
  return ans;
}

double SphericalModelTable::get_pot(const double r)
{
  double ans;

  if (r<pot.x[1]) {
    if (diverge)
      return -1.0/(r+1.0);
    else
      return pot.y[1];
  }

  if (r>pot.x[pot.num]) return pot.y[pot.num]*pot.x[pot.num]/r;
  
  Splint1(pot.x, pot.y, pot.y2, r, ans, even);
  return ans;
}

double SphericalModelTable::get_dpot(const double r)
{
  double dum, ans;

  if (r<pot.x[1]) {
    if (diverge)
      return 1.0/( (r+1.0)*(r+1.0) );
    else
      Splint2(pot.x, pot.y, pot.y2, pot.x[1], dum, ans, even);
  }
  else if (r>pot.x[pot.num]) 
    ans = -pot.y[pot.num]*pot.x[pot.num]/(r*r);
  else
    Splint2(pot.x, pot.y, pot.y2, r, dum, ans, even);

  return ans;
}

void SphericalModelTable::get_pot_dpot(const double r, double& ur, double &dur)
{
  if (r<pot.x[1]) {
    if (diverge) {
      ur = -1.0/(r+1.0);
      dur = 1.0/( (r+1.0)*(r+1.0) );
    }
    else
      Splint2(pot.x, pot.y, pot.y2, pot.x[1], ur, dur, even);
  }
  else if (r>pot.x[pot.num]) {
    ur = pot.y[pot.num]*pot.x[pot.num]/r;
    dur = -pot.y[pot.num]*pot.x[pot.num]/(r*r);
  }
  else {
    Splint2(pot.x, pot.y, pot.y2, r, ur, dur, even);
  }

}

double SphericalModelTable::get_dpot2(const double r)
{
  double dum, ans;

  if (r<pot.x[1]) {
    if (diverge)
      return -2.0/( (r+1.0)*(r+1.0)*(r+1.0) );
    else
      Splint2(pot.x, pot.y, pot.y2, pot.x[1], dum, ans, even);
  }
  else if (r>pot.x[pot.num]) 
    ans = -2.0*pot.y[pot.num]*pot.x[pot.num]/(r*r*r);
  else
    Splint3(pot.x, pot.y, pot.y2, r, dum, dum, ans, even);

  return ans;
}

