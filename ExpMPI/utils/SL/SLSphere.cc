#include <SLSphere.H>

SLSphere::SLSphere(int Lmax, int Nmax, int Numr, double Rmin, double Rmax,
		   int cmap, double rs)
{
  dof = 3;

  lmax = Lmax;
  nmax = Nmax;
  numr = Numr;
  rmin = Rmin;
  rmax = Rmax;

				// Generate Sturm-Liouville grid
  ortho = new SLGridSph(Lmax, nmax, numr, rmin, rmax, 
			cmap, rs);
}

SLSphere::~SLSphere()
{
  delete ortho;
}


