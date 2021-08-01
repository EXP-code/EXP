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
  ortho = boost::make_shared<SLGridSph>
    (Lmax, nmax, numr, rmin, rmax, false, cmap, rs);

  xmin = ortho->r_to_xi(rmin);
  xmax = ortho->r_to_xi(rmax);
}

SLSphere::SLSphere(int Lmax, int Nmax, int Numr, double Rmin, double Rmax,
		   int cmap, double rs,
		   boost::shared_ptr<SphericalModelTable> mod)
{
  dof = 3;

  lmax = Lmax;
  nmax = Nmax;
  numr = Numr;
  rmin = Rmin;
  rmax = Rmax;

				// Generate Sturm-Liouville grid
  if (mod)
    ortho = boost::make_shared<SLGridSph>
      (Lmax, nmax, numr, rmin, rmax, mod, false, cmap, rs);
  else
    ortho = boost::make_shared<SLGridSph>
      (Lmax, nmax, numr, rmin, rmax, false, cmap, rs);

  xmin = ortho->r_to_xi(rmin);
  xmax = ortho->r_to_xi(rmax);
}

SLSphere::~SLSphere()
{
  // NADA
}


