/* 
   odes - Integrate Lane-Emden equation using RK4

   Equations of motion:

   y_1 = \Psi
   y_2 = d\Psi/dx = dy_1/dx

   f_1 = y_2
   f_2 = exp(-y_1) - 2y_2/x


   Copyright (C) 2008 Martin Weinberg

   This program is free software; you can redistribute it and/or modify
   it under the terms of the GNU General Public License as published by
   the Free Software Foundation; either version 2, or (at your option)
   any later version.

   This program is distributed in the hope that it will be useful,
   but WITHOUT ANY WARRANTY; without even the implied warranty of
   MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
   GNU General Public License for more details.

   You should have received a copy of the GNU General Public License
   along with this program; if not, write to the Free Software Foundation,
   Inc., 59 Temple Place - Suite 330, Boston, MA 02111-1307, USA.  

*/

#define VERSION 0.1

#include <iostream>
#include <iomanip>
#include <sstream>
#include <fstream>
#include <vector>
#include <cmath>
#include <cstdlib>

#include <ACG.h>
#include <Normal.h>
#include <Uniform.h>

using namespace std;

#include <sys/types.h>
#include <getopt.h>

#define EXIT_FAILURE 1

char threading_on = 0;
pthread_mutex_t mem_lock;

static void usage (int status);

/* The name the program was run with, stripped of any leading path. */
char *program_name;

/* Option flags and variables */

string oname;			/* --output */

static struct option const long_options[] =
{
  {"output",  required_argument, 0,    'o'},
  {"xfinal",  required_argument, 0,    'x'},
  {"xstep",   required_argument, 0,    'd'},
  {"ratio",   required_argument, 0,    'r'},
  {"mass",    required_argument, 0,    'm'},
  {"temp",    required_argument, 0,    'T'},
  {"runit",   required_argument, 0,    'R'},
  {"number",  required_argument, 0,    'N'},
  {"help",    no_argument,       0,    'h'},
  {"version", no_argument,       0,    'V'},
  {NULL,      0,                 NULL,  0 }
};

static int decode_switches (int argc, char **argv);

double X = 1000.0;
double h = 0.01;
double M = 1.0;			// units of 10^10 solar masses
double ratio = 14.0;		// rho_c/rho_t
double T = 10000.0;		// degrees kelvin
double R = 300.0;		// Unit dimension in kpc
unsigned N = 0;			// Number of particles
unsigned S = 11;		// random # seed

void deriv(double x, vector<double>&y, vector<double>&a)
{
  /*
   y_1 = \Psi
   y_2 = d\Psi/dx = dy_1/dx

   f_1 = y_2
   f_2 = exp(-y_1) - 2y_2/x
  */

  a[0] = y[1];
  if (x>0.0) {
    a[0] = y[1];
    a[1] = exp(-y[0]) - 2.0*y[1]/x;
  } else {
    a[0] = 0.0;
    a[1] = 1.0;
  }
}

void rk4(double t, vector<double>& x, double h)
{
  unsigned dim = x.size();
  vector<double> x1(x), k1(dim), k2(dim), k3(dim), k4(dim);

  deriv(t, x1, k1);

  for (unsigned j=0; j<dim; j++) x1[j] = x[j] + 0.5*h*k1[j];
  deriv(t+0.5*h, x1, k2);

  for (unsigned j=0; j<dim; j++) x1[j] = x[j] + 0.5*h*k2[j];
  deriv(t+0.5*h, x1, k3);

  for (unsigned j=0; j<dim; j++) x1[j] = x[j] + h*k3[j];
  deriv(t+h, x1, k4);

  for (unsigned j=0; j<dim; j++) x[j] += h*(k1[j]+2.0*k2[j]+2.0*k3[j]+k4[j])/6.0;
}

vector< vector<double> > solution;

void accum(double& x, vector<double>& y)
{
  double rho = exp(-y[0]), xxy=x*x*y[1];
  vector<double> record(7);
  record[0] = x;
  record[1] = y[0];
  record[2] = y[1];
  record[3] = xxy;
  record[4] = rho;
  record[5] = sqrt(0.25*rho/M_PI)*xxy;
  record[6] = 1.0/record[4];

  solution.push_back(record);
}


int vlocate(double x, unsigned j, const vector< vector<double> > &xx)
{
  int imin = 0;
  int imax = xx.size() - 1;
  int lo = imin - 1;
  int hi = imax + 1;
  int ascnd = (xx[imax][j] > xx[imin][j]);
  int mid;
  while (hi-lo > 1) {
    mid =(lo+hi) >> 1;
    if ((x > xx[mid][j]) == ascnd)
      lo = mid;
    else
      hi = mid;
  }
  return lo;
}

double linear(double x, unsigned j, const vector< vector<double> >& tab)
{
  int imin = 0;
  int imax = tab.size() - 1;

  int index = vlocate(x, 0, tab);

  if (index <  imin) index = imin;
  if (index >= imax) index = imax-1;

  return ( tab[j][index+1]*(x-tab[0][index]  ) -
	   tab[j][index  ]*(x-tab[0][index+1]) )
    /( tab[0][index+1]-tab[0][index] ) ;
}

int
main (int argc, char **argv)
{
  int i;

  program_name = argv[0];

  i = decode_switches (argc, argv);

  //
  // Output header
  //

  ofstream fout;

  if (oname.size())
    fout.open(oname.c_str());
  else {
    fout.copyfmt(std::cout);
    fout.clear(std::cout.rdstate());
    fout.basic_ios<char>::rdbuf(std::cout.rdbuf());
  }

  if (!fout) {
    cerr << program_name << ": error opening output" << endl;
    exit(-1);
  }

  //
  // Initial condtions
  //

  vector<double> y(2), a(2);

  double x = 0.0;
  y[0] = 0.0;
  y[1] = 0.0;

  accum(x, y);

  while (x<=X) {
    x += h;
				// RK4 step
    rk4(x, y, h);
				// Output
    accum(x, y);
  }

  
  //
  // Get end point
  //

  int n;

  if (ratio>=solution.back()[6])
    n = solution.size()-2;
  else if (ratio<=solution.front()[6]) 
    n = 0;
  else 
    n = vlocate(ratio, 6, solution);

  double denom = solution[n+1][6] - solution[n][6];
  double A = (solution[n+1][6] - ratio)/denom;
  double B = (ratio - solution[n][6])/denom;

  double xt = A*solution[n][0] + B*solution[n+1][0];
  double mt = A*solution[n][5] + B*solution[n+1][5];

  //
  // Units
  //

  const double mp = 1.67262158e-24; // Proton mass (g)
  const double boltz = 1.3810e-16;  // Boltzmann constant (cgs)
  const double f_H = 0.76;	    // Hydrogen fraction
  static double pc = 3.086e18;	    // cm
  static double msun = 1.989e33;    // g
  static double G = 6.67428e-8;	    // cgs
  
  double mm = f_H*mp + (1.0-f_H)*4.0*mp;
  double cs2 = boltz*T/mm;

  double Pt = mt*cs2*cs2/(pow(G, 1.5)*M*1e10*msun); Pt *= Pt;
  double Rhot = Pt/cs2;
  double Rhoc = Rhot*ratio;

				// in kiloparsec
  double rfac   = sqrt(cs2/(4.0*M_PI*G*Rhoc))/(R*1.0e3*pc);
				// solar masses/kpc^3
  double rhofac = Rhoc*pow(1.0e3*pc, 3.0)/msun;
				// in units of 10^10 solar masses
  double mfac   = cs2*cs2/(sqrt(Pt)*pow(G, 1.5))/(1e10*msun);

  if (N==0) {

    for (unsigned i=0; i<=n; i++) {
      fout << setw(15) << solution[i][0]*rfac;
      fout << setw(15) << solution[i][4]*rhofac;
      fout << setw(15) << solution[i][5]*mfac;
      fout << endl;
    }
    fout << setw(15) << (A*solution[n][0] + B*solution[n+1][0])*rfac;
    fout << setw(15) << (A*solution[n][4] + B*solution[n+1][4])*rhofac;
    fout << setw(15) << (A*solution[n][5] + B*solution[n+1][5])*mfac;
    fout << endl;

  } else {

				// Circ velocity unit at edge in cgs
    double vcirc = sqrt(G*M*1.0e10*msun/(R*1e3*pc));
    double vfac = sqrt(cs2)/vcirc;

    ACG gen(S);
    Uniform unit(0.0, 1.0, &gen);
    Normal  norm(0.0, 1.0, &gen);
    vector<double> pos(3);
    double mass = M/N, x, phi, cost, sint, m;
    const int ITMAX = 1000;
    const double MMAX = 1.2;

    fout << " " << N << " 0 4" << endl; // PSP header
    fout.precision(10);

    for (int i=0; i<N; i++) {

				// Evaluate by rejection so we can
				// do unstable spheres later
      unsigned j;
      for (j=0; j<ITMAX; j++) {
	x = xt*unit();
	m = linear(x, 5, solution);
	if (unit() < m/MMAX) break;
      }
      if (j==ITMAX) { cerr << "Oops" << endl; }

      x *= rfac;

      phi  = 2.0*M_PI*unit();
      cost = 2.0*unit() - 1.0;
      sint = sqrt(1.0 - cost*cost);
      pos[0] = x*sint*cos(phi);
      pos[1] = x*sint*sin(phi);
      pos[2] = x*cost;
      
      fout << setw(18) << mass;
      for (int k=0; k<3; k++) fout << setw(18) << pos[k];
      for (int k=0; k<3; k++) fout << setw(18) << vfac*norm();
      for (int k=0; k<4; k++) fout << setw(18) << 0.0;
      fout << endl;
    }

  }

  exit (0);
}

// Set all the option flags according to the switches specified.
// Return the index of the first non-option argument.

static int
decode_switches (int argc, char **argv)
{
  int c;

  while ((c = getopt_long (argc, argv, 
			   "o:"	/* output  */
			   "x:"	/* xfinal  */
			   "d:"	/* xstep   */
			   "r:"	/* ratio   */
			   "m:"	/* mass    */
			   "T:"	/* temp    */
			   "R:"	/* runit   */
			   "N:"	/* number  */
			   "S:"	/* seed    */
			   "h"	/* help    */
			   "V",	/* version */
			   long_options, (int *) 0)) != EOF)
    {
      switch (c)
	{
	case 'o':		/* --output */
	  oname = string(optarg);
	  break;
	case 'x':		/* --xfinal */
	  X = atof(optarg);
	  break;
	case 'd':		/* --xstep */
	  h = atof(optarg);
	  break;
	case 'r':		/* --ratio */
	  ratio = atof(optarg);
	  break;
	case 'R':		/* --runit */
	  R = atof(optarg);
	  break;
	case 'm':		/* --mass */
	  M = atof(optarg);
	  break;
	case 'T':		/* --temp */
	  T = atof(optarg);
	  break;
	case 'N':		/* --number */
	  N = atoi(optarg);
	  break;
	case 'S':		/* --seed */
	  S = atoi(optarg);
	  break;
	case 'V':
	  cout << program_name << " " << VERSION << endl;
	  exit (0);

	case 'h':
	  usage (0);

	default:
	  usage (-1);
	}
    }

  return optind;
}


static void
usage (int status)
{
  cout << program_name << " - Integrate Lane-Emden equation using RK4" << endl;
  cout << "Usage: " << program_name << " [OPTION]... [FILE]..." << endl;
  cout << "Options:" << endl
       << "  -o, --output NAME          send output to NAME instead of standard output" << endl
       << "  --verbose                  print more information" << endl
       << "  -h, --help                 display this help and exit" << endl
       << "  -V, --version              output version information and exit" << endl;
  exit (status);
}
