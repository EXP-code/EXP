#include "SLSphere.H"

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
  ortho = std::make_shared<SLGridSph>
    ("", Lmax, nmax, numr, rmin, rmax, false, cmap, rs, 0, 1.0);

  xmin = ortho->r_to_xi(rmin);
  xmax = ortho->r_to_xi(rmax);
}

SLSphere::SLSphere(int Lmax, int Nmax, int Numr, double Rmin, double Rmax,
		   int cmap, double rs,
		   std::shared_ptr<SphericalModelTable> mod)
{
  dof = 3;

  lmax = Lmax;
  nmax = Nmax;
  numr = Numr;
  rmin = Rmin;
  rmax = Rmax;

				// Generate Sturm-Liouville grid
  if (mod)
    ortho = std::make_shared<SLGridSph>
      (mod, Lmax, nmax, numr, rmin, rmax, false, cmap, rs);
  else
    ortho = std::make_shared<SLGridSph>
      ("", Lmax, nmax, numr, rmin, rmax, false, cmap, rs, 0, 1.0);

  xmin = ortho->r_to_xi(rmin);
  xmax = ortho->r_to_xi(rmax);
}

SLSphere::~SLSphere()
{
  // NADA
}


