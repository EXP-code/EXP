#ifndef _SLSphere_H
#define _SLSphere_H

static const char rcsid_SLSphere[] = "$Id$";


#include <biorth.h>
#include <SLGridMP2.h>

class SLSphere : public AxiSymBiorth
{
private:

  SLGridSph *ortho;
  int lmax, nmax, numr;
  double rmin, rmax;
  double xmin, xmax;

public:
  
  SLSphere(int Lmax, int Nmax, int Numr, double Rmin, double Rmax,
	   int cmap, double rs);
  ~SLSphere();

  double potl(const int n, const int l, const double x) 
  {
    ortho->get_pot(x, l, n, 0);
  }

  double dens(const int n, const int l, const double x)
  {
    ortho->get_dens(x, l, n, 0);
  }

  double potlR(const int n, const int l, const double r)
  {
    ortho->get_pot(r, l, n, 1);
  }

  double potlRZ(const int n, const int l, const double R, const double Z)
  {
    ortho->get_pot(sqrt(R*R+Z*Z), l, n, 1);
  }
				// Fast for expansion! 
				// (using recursion relation)
  void potl(const int nmax, const int l, const double r, Vector& vec)
  {
    ortho->get_pot(vec, r, l, 1);
  }

  void dens(const int nmax, const int l, const double r, Vector& vec)
  {
    ortho->get_dens(vec, r, l, 1);

  }

  double rb_to_r(double const rb)   {return ortho->xi_to_r(rb);}
  double r_to_rb(double const r )   {return ortho->r_to_xi(r);}
  double d_r_to_rb(double const rb) {return ortho->d_xi_to_r(rb);}

  double rb_min(void) {return xmin;}
  double rb_max(void) {return xmax;}

  int get_dof(void) { return dof; }
  double krnl(const int l, const int n) { return 1; }
  double norm(const int l, const int n) { return 1; }

};


#endif
