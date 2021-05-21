#include <iostream>
#include <iomanip>
#include <sstream>
#include <fstream>
#include <cmath>

#include <Vector.h>
#include <SphereSL.H>

#ifdef HAVE_LIBPNGPP
#include <ColorGradient.H>	// For PNG images
#endif

int    SphereSL::NUMR = 800;
int    SphereSL::NEV  = 60;	// None by default
bool   SphereSL::mpi  = false;	// Initially off
double SphereSL::HEXP = 1.0;	// Hall exponent

//! Constructor
SphereSL::SphereSL(SphericalModelTable* mod, int LMAX, int NMAX,
		   int CMAP, double Scale, bool COVAR, int NPART)
{
  SphericalModelTable *model = mod;
  
  SLGridSph::mpi = mpi ? 1 : 0;

  double rmin = max<double>(mod->get_min_radius()*2.0, 
			    mod->get_max_radius()*1.0e-4);
  double rmax = mod->get_max_radius()*0.99;

  sl = new SLGridSph(LMAX, NMAX, NUMR, rmin, rmax, model, CMAP, Scale);

  lmax = LMAX;
  nmax = NMAX;
  compute_covar = COVAR;
  
  coefs_defined = false;
  rscl = 1.0;

  potd.setsize  (0, LMAX, 1, NMAX);
  dpot.setsize  (0, LMAX, 1, NMAX);
  dpt2.setsize  (0, LMAX, 1, NMAX);
  dend.setsize  (0, LMAX, 1, NMAX);
    
  legs.setsize  (0, LMAX, 0, LMAX);
  dlegs.setsize (0, LMAX, 0, LMAX);
  d2legs.setsize(0, LMAX, 0, LMAX);

  npart = NPART;
}

SphereSL::~SphereSL(void)
{
  delete sl;
}

void SphereSL::bomb(char *s)
{
  cerr << "ERROR from SphereSL: " << s << '\n';
  exit(-1);
}

void SphereSL::reset_coefs(void)
{
  if (expcoef.getnrows()>0 && expcoef.getncols()>0) expcoef.zero();
  if (compute_covar) {
    minSNR = std::numeric_limits<double>::max();
    maxSNR = 0.0;
    totalMass = 0.0;
    covar.resize(lmax+1);
    mean.resize(lmax+1);
    for (int L=0; L<=lmax; L++) {
      int size = (L+1)*nmax;
      covar[L] = Eigen::MatrixXd::Zero(size, size);
      mean [L] = Eigen::VectorXd::Zero(size);
    }
    if (npart) {
      curbin = 0;
      meanB.resize(npart);
      massB.resize(npart);
      for (int n=0; n<npart; n++) {
	meanB[n].resize(lmax+1);
	for (int L=0; L<=lmax; L++) {
	  int size = (L+1)*nmax;
	  meanB[n][L] = Eigen::VectorXd::Zero(size);
	}
	massB[n] = 0.0;
      }
    }
  }
}


void SphereSL::accumulate(double x, double y, double z, double mass)
{
  double fac, fac1, fac2, fac4;
  double norm = -4.0*M_PI;
  const double dsmall = 1.0e-20;

  if (!coefs_defined) {

    coefs_defined = true;

    expcoef.setsize(0, lmax*(lmax+2), 1, nmax);
    expcoef.zero();		// Need this?

    work.setsize(1, nmax);

    if (compute_covar) {
      minSNR = std::numeric_limits<double>::max();
      maxSNR = 0.0;
      totalMass = 0.0;
      covar.resize(lmax+1);
      mean.resize(lmax+1);
      for (int L=0; L<=lmax; L++) {
	int size = (L+1)*nmax;
	covar[L] = Eigen::MatrixXd::Zero(size, size);
	mean [L] = Eigen::VectorXd::Zero(size);
      }
      if (npart) {
	curbin = 0;
	meanB.resize(npart);
	massB.resize(npart);
	for (int n=0; n<npart; n++) {
	  meanB[n].resize(lmax+1);
	  for (int L=0; L<=lmax; L++) {
	    int size = (L+1)*nmax;
	    meanB[n][L] = Eigen::VectorXd::Zero(size);
	  }
	  massB[n] = 0.0;
	}
      }
    }

    factorial.setsize(0, lmax, 0, lmax);

    for (int l=0; l<=lmax; l++) {
      for (int m=0; m<=l; m++) {
	factorial[l][m] = sqrt( (0.5*l+0.25)/M_PI * 
				exp(lgamma(1.0+l-m) - lgamma(1.0+l+m)) );
	if (m != 0) factorial[l][m] *= M_SQRT2;
      }
    }

    used = 0;
  }

  //======================
  // Compute coefficients 
  //======================

  double r2 = (x*x + y*y + z*z);
  double r = sqrt(r2) + dsmall;
  double costh = z/r;
  double phi = atan2(y,x);
  double rs = r/rscl;
	
  used++;
  totalMass += mass;

  sl->get_pot(potd, rs);

  legendre_R(lmax, costh, legs);

  // L loop
  for (int l=0, loffset=0; l<=lmax; loffset+=(2*l+1), l++) {

    Eigen::VectorXd workE;
    int esize = (l+1)*nmax;

    if (compute_covar) workE = Eigen::VectorXd::Zero(esize);

    // M loop
    for (int m=0, moffset=0, moffE=0; m<=l; m++) {

      if (m==0) {
	fac = factorial[l][m] * legs[l][m];
	for (int n=1; n<=nmax; n++) {
	  fac4 = potd[l][n]*fac;
	  expcoef[loffset+moffset][n] += fac4 * norm *mass;
	  if (compute_covar) workE[m*nmax + n - 1] = fac4;
	}

	moffset++;
      }
      else {
	fac  = factorial[l][m] * legs[l][m];
	fac1 = fac*cos(phi*m);
	fac2 = fac*sin(phi*m);
	for (int n=1; n<=nmax; n++) {
	  fac4 = potd[l][n];
	  expcoef[loffset+moffset  ][n] += fac1 * fac4 * norm * mass;
	  expcoef[loffset+moffset+1][n] += fac2 * fac4 * norm * mass;
	  if (compute_covar) workE[m*nmax + n - 1] = fac * fac4;
	}

	moffset+=2;
      }
    }

    if (compute_covar) {
      for (int m=0; m<=l; m++) {
	for (int n1=0; n1<nmax; n1++) {
	  mean[l][m*nmax+n1] += workE[m*nmax + n1] * mass;

	  if (npart)
	    meanB[curbin % npart][l][m*nmax + n1] += workE[m*nmax + n1] * mass;
	  
	  if (npart==0) {
	    for (int n2=0; n2<nmax; n2++) {
	      covar[l](m*nmax + n1, m*nmax + n2) +=
		workE[m*nmax + n1] * workE[m*nmax + n2] * mass;
	    }
	  }
	}
      }

    }
    // END: compute_covar

  }

  if (compute_covar and npart) {
    massB[curbin % npart] += mass;
    curbin++;
  }

}

void SphereSL::make_coefs()
{
  if (mpi) {

    MPI_Allreduce(MPI_IN_PLACE, &used, 1, MPI_INT,
		  MPI_SUM, MPI_COMM_WORLD);

    if (compute_covar)
      MPI_Allreduce(MPI_IN_PLACE, &totalMass, 1, MPI_DOUBLE,
		    MPI_SUM, MPI_COMM_WORLD);

    for (int l=0; l<=lmax*(lmax+2); l++) {
      MPI_Allreduce(&expcoef[l][1], &work[1], nmax, MPI_DOUBLE,
		    MPI_SUM, MPI_COMM_WORLD);

      expcoef[l] = work;
    }

    if (compute_covar) {
      for (int l=0; l<=lmax; l++) {
	int esize = (l + 1)*nmax;

	if (npart==0)
	  MPI_Allreduce(MPI_IN_PLACE, covar[l].data(), esize*esize, MPI_DOUBLE,
			MPI_SUM, MPI_COMM_WORLD);

	MPI_Allreduce(MPI_IN_PLACE, mean[l].data(), esize, MPI_DOUBLE,
		      MPI_SUM, MPI_COMM_WORLD);
	
	if (npart) {
	  for (int n=0; n<npart; n++) {
	    MPI_Allreduce(MPI_IN_PLACE, meanB[n][l].data(), esize, MPI_DOUBLE,
			  MPI_SUM, MPI_COMM_WORLD);
	  }
	}
      }
      
      if (npart) {
	MPI_Allreduce(MPI_IN_PLACE, &curbin, 1, MPI_INT,
		      MPI_SUM, MPI_COMM_WORLD);
	  
	MPI_Allreduce(MPI_IN_PLACE, massB.data(), massB.size(), MPI_DOUBLE,
		      MPI_SUM, MPI_COMM_WORLD);
      }
    }
  }

  if (compute_covar) {

    for (int l=0; l<=lmax; l++) {
      int esize = (l + 1)*nmax;

      if (npart) {

	// Zero mean and covariance matrix
	//
	mean [l] *= 0.0;
	covar[l] *= 0.0;

	// Normalize by mass and zero covar
	//
	for (int n=0; n<npart; n++) {
	  if (massB[n]>0.0) meanB[n][l] /= massB[n];
	}

	for (int n=0; n<npart; n++) {
	  for (int i=0; i<esize; i++) mean[l][i] += meanB[n][l][i];
	}
	mean[l] /= npart;

	for (int n=0; n<npart; n++) {
	  for (int i=0; i<esize; i++) {
	    for (int j=0; j<esize; j++) {
	      covar[l](i, j) +=
		(meanB[n][l](i)  - mean[l](i)) *
		(meanB[n][l](j)  - mean[l](j)) / npart;
	    }
	  }
	}

	// Convert to 1-particle covariance using CLT
	//
	covar[l] *= static_cast<double>(used)/npart;

      } else {

	if (totalMass>0.0) {
	  mean[l]  /= totalMass;
	  covar[l] /= totalMass;

	  for (int i=0; i<esize; i++) {
	    for (int j=0; j<esize; j++) {
	      covar[l](i, j) -=	mean[l](i) * mean[l](j);
	    }
	  }
	}
	// END: totalMass>0.0
      }
    }
  }
}

void SphereSL::make_covar(bool verbose)
{
  if (compute_covar) {

    if (verbose and myid==0) std::cout << std::endl;

    svar.resize(lmax+1);
    uvec.resize(lmax+1);

    std::vector<double> totpow(lmax+1, 0.0);

    for (int l=0; l<=lmax; l++) {
      int esize = (l + 1)*nmax;
      
#ifdef LARGE
      Eigen::BDCSVD<Eigen::MatrixXd>
	svd(covar[l], Eigen::ComputeThinU | Eigen::ComputeThinV);
#else
      Eigen::JacobiSVD<Eigen::MatrixXd>
	svd(covar[l], Eigen::ComputeThinU | Eigen::ComputeThinV);
#endif

      // Get solution
      // ------------
      svar[l] = svd.singularValues();
      uvec[l] = svd.matrixU();

      // Get value in rotated space
      //
      // 
      Eigen::VectorXd R = uvec[l].transpose() * mean[l];

      // Compute SNR
      //
      for (int j=0; j<svar[l].size(); j++) {

	totpow[l] += R[j]*R[j];	// Accumulate total power

	if (verbose and myid==0) std::cout << std::setw( 4) << l
					   << std::setw( 4) << j
					   << std::setw(18) << svar[l][j]
					   << std::setw(18) << R[j]*R[j];
	if (svar[l][j]>0.0) {
	  double snr = R[j]*R[j]/svar[l][j];
	  minSNR = std::min<double>(minSNR, snr);
	  maxSNR = std::max<double>(maxSNR, snr);
	  if (verbose and myid==0) std::cout << std::setw(18) << snr;
	} else {
	  if (verbose and myid==0) std::cout << std::setw(18) << 0.0;
	}
	if (verbose and myid==0) std::cout << std::endl;
      }
    }

    if (verbose and myid==0) {
      std::cout << "Total power" << std::endl;
      for (int l=0; l<=lmax; l++)
	std::cout << std::setw(4) << l << std::setw(18) << totpow[l]
		  << std::endl;
    }
  }
}

Matrix SphereSL::get_trimmed(double snr, double mass, bool Hall)
{
  constexpr double norm = 4.0*M_PI;

  Matrix ret(0, lmax*(lmax+2), 1, nmax);
  ret.zero();

  if (compute_covar) {
    
    // L loop
    //
    for (int l=0, loffset=0; l<=lmax; loffset+=(2*l+1), l++) {

      int esize = (l+1)*nmax;
      
      // Test create new vector
      //
      Eigen::VectorXd W(esize);
      for (int m=0, moffset=0; m<=l; m++) {
	if (m==0) {
	  for (int n=1; n<=nmax; n++) {
	    W[m*nmax + n - 1] = fabs(expcoef[loffset+moffset+0][n])
	      / (norm*mass);
	  }
	  moffset++;

	} else {
	  for (int n=1; n<=nmax; n++) {
	    W[m*nmax + n - 1] =
	      sqrt(expcoef[loffset+moffset+0][n]*expcoef[loffset+moffset+0][n]
		   +
		   expcoef[loffset+moffset+1][n]*expcoef[loffset+moffset+1][n])
	      / (norm*mass);
	  }

	  moffset+=2;
	}
      }

      // Assign nullity
      //
      Eigen::VectorXd R = uvec[l].transpose() * W;
      for (int j=0; j<svar[l].size(); j++) {
	if (svar[l][j]>0.0) {
	  if (Hall)
	    R[j] *= 1.0/(pow(snr*svar[l][j]/(R[j]*R[j]), HEXP) + 1.0);
	  else if (R[j]*R[j]/svar[l][j] < snr)
	    R[j] = 0.0;
	} else {
	  R[j] = 0.0;		// Sanity
	}
      }

      // Rotate back to original space
      //
      Eigen::VectorXd Q = uvec[l] * R;

      // M loop
      //
      for (int m=0, moffset=0; m<=l; m++) {
	if (m==0) {
	  for (int n=1; n<=nmax; n++) {
	    ret[loffset+moffset][n] = Q[m*nmax + n - 1] * norm * mass;
	    if (expcoef[loffset+moffset][n] < 0.0)
	      ret[loffset+moffset][n] *= -1.0 ;
	  }
	  moffset++;

	} else {
	  for (int n=1; n<=nmax; n++) {
	    double phi = atan2(expcoef[loffset+moffset+1][n], expcoef[loffset+moffset+0][n]);
	    ret[loffset+moffset+0][n] = cos(phi) * Q[m*nmax + n - 1] * norm * mass;
	    ret[loffset+moffset+1][n] = sin(phi) * Q[m*nmax + n - 1] * norm * mass;
	  }

	  moffset+=2;
	}
      }
    }
  }



#ifdef HAVE_LIBPNGPP
  if (NEV) {

    const int minSize = 600;
    int ndupX = 1, ndupY = 1;
				// Sanity check
    NEV = std::min<int>(NEV, nmax);
    
    for (int L=0; L<=lmax; L++) {

      int esize = (L+1)*nmax;
      
      if (NEV < minSize)   ndupX = minSize/NEV   + 1;
      if (esize < minSize) ndupY = minSize/esize + 1;


      png::image< png::rgb_pixel > image(NEV*ndupX, esize*ndupY);
      ColorGradient color;
      // color.createFiveColorHeatMapGradient();
      color.createGrayGradient();

      std::ostringstream sout;
      sout << "SphereSL_EV." << L;

      double minV = std::numeric_limits<double>::max();
      double maxV = std::numeric_limits<double>::min();

      if (true) {

	for (int ev=0; ev<NEV; ev++) {
	  for (int n=0; n<esize; n++) {
	    maxV = std::max<double>(maxV, uvec[L](ev, n)*uvec[L](ev, n));
	  }
	}
	
	for (int i=0; i<esize; i++) {
	  for (int j=0; j<NEV; j++) {
	    png::rgb_pixel cval = color(uvec[L](i, j)*uvec[L](i, j)/maxV );
	    for (size_t yy = i*ndupY; yy < (i+1)*ndupY; yy++) {
	      for (size_t xx = j*ndupX; xx < (j+1)*ndupX; xx++) {
		image[yy][xx] = cval;
	      }
	    }
	  }
	}

      } else {

	for (int ev=0; ev<NEV; ev++) {
	  for (int n=0; n<esize; n++) {
	    minV = std::min<double>(minV, uvec[L](ev, n));
	    maxV = std::max<double>(maxV, uvec[L](ev, n));
	  }
	}
	
	for (int i=0; i<esize; i++) {
	  for (int j=0; j<NEV; j++) {
	    png::rgb_pixel cval = color( (uvec[L](i, j) - minV)/(maxV - minV) );
	    for (size_t yy = i*ndupY; yy < (i+1)*ndupY; yy++) {
	      for (size_t xx = j*ndupX; xx < (j+1)*ndupX; xx++) {
		image[yy][xx] = cval;
	      }
	    }
	  }
	}
      }

      image.write(sout.str() + ".png");

    }
  }
#endif

  return ret;
}

double SphereSL::get_power(double snr, double mass)
{
  constexpr double norm = 4.0*M_PI;

  double total = 0.0, trimmed = 0.0;

  if (compute_covar) {
    
    // L loop
    for (int l=0, loffset=0; l<=lmax; loffset+=(2*l+1), l++) {

      int esize = (l+1)*nmax;
      
      // Test create new vector
      //
      Eigen::VectorXd W(esize);
      for (int m=0, moffset=0; m<=l; m++) {
	if (m==0) {
	  for (int n=1; n<=nmax; n++) {
	    W[m*nmax + n - 1] = fabs(expcoef[loffset+moffset+0][n])
	      / (norm*mass);
	  }
	  moffset++;

	} else {
	  for (int n=1; n<=nmax; n++) {
	    W[m*nmax + n - 1] =
	      sqrt(expcoef[loffset+moffset+0][n]*expcoef[loffset+moffset+0][n]
		   +
		   expcoef[loffset+moffset+1][n]*expcoef[loffset+moffset+1][n])
	      / (norm*mass);
	  }

	  moffset+=2;
	}
      }

      // Assign nullity
      //
      Eigen::VectorXd R = uvec[l].transpose() * W;
      for (int j=0; j<svar[l].size(); j++) {
	total += svar[l][j];
	if (svar[l][j]>0.0) {
	  if (R[j]*R[j]/svar[l][j] < snr) trimmed += svar[l][j];
	}
      }
    }
  }

  if (total>0.0) return trimmed/total;
  else return 0.0;
}


void SphereSL::dens_pot_eval(double r, double costh, double phi,
			     double& dens0, double& dens, 
			     double& potl0, double& potl,
			     int L1, int L2)
{
  double fac1, cosm, sinm;

  fac1 = factorial[0][0];

  sl->get_dens(dend, r/rscl);
  sl->get_pot (potd, r/rscl);

  legendre_R(lmax, costh, legs, dlegs);

  dens0 = fac1 * expcoef[0]*dend[0];
  potl0 = fac1 * expcoef[0]*potd[0];

  dens = 0.0;
  potl = 0.0;

  // L loop
  for (int l=1, loffset=1; l<=lmax; loffset+=(2*l+1), l++) {
    if (l<L1 || l>L2) continue;

    // M loop
    for (int m=0, moffset=0; m<=l; m++) {
      fac1 = factorial[l][m];
      if (m==0) {
	dens += fac1*legs[l][m]* (expcoef[loffset+moffset] * dend[l]);

	potl += fac1*legs[l][m]* (expcoef[loffset+moffset] * potd[l]);

	moffset++;
      }
      else {
	cosm = cos(phi*m);
	sinm = sin(phi*m);

	dens += fac1*legs[l][m]*
	  ( expcoef[loffset+moffset]   * dend[l]*cosm + 
	    expcoef[loffset+moffset+1] * dend[l]*sinm );

	potl += fac1*legs[l][m]*
	  ( expcoef[loffset+moffset]   * potd[l]*cosm + 
	    expcoef[loffset+moffset+1] * potd[l]*sinm );

	moffset +=2;
      }
    }
  }

  double densfac = 1.0/(rscl*rscl*rscl) * 0.25/M_PI;
  double potlfac = 1.0/rscl;

  dens0 *= densfac;
  dens  *= densfac;
  potl0 *= potlfac;
  potl  *= potlfac;
}


void SphereSL::pot_force_eval(double r, double costh, double phi,
			      double& potl,
			      double& potr, double& pott, double& potp,
			      int L1, int L2)
{
  double fac1, cosm, sinm;
  double sinth = -sqrt(fabs(1.0 - costh*costh));

  fac1 = factorial[0][0];

  sl->get_pot  (potd, r/rscl);
  sl->get_force(dpot, r/rscl);

  legendre_R(lmax, costh, legs, dlegs);

  potl = fac1 * expcoef[0]*potd[0];
  potr = fac1 * expcoef[0]*dpot[0];
  pott = 0.0;
  potp = 0.0;

  // L loop
  for (int l=1, loffset=1; l<=lmax; loffset+=(2*l+1), l++) {
    if (l<L1 || l>L2) continue;

    // M loop
    for (int m=0, moffset=0; m<=l; m++) {
      fac1 = factorial[l][m];
      if (m==0) {
	potl += fac1*legs[l][m]* (expcoef[loffset+moffset] * potd[l]);
	dpot += fac1*legs[l][m]* (expcoef[loffset+moffset] * dpot[l]);
	pott += fac1*dlegs[l][m]* (expcoef[loffset+moffset] * potd[l]);

	moffset++;
      }
      else {
	cosm = cos(phi*m);
	sinm = sin(phi*m);

	potl += fac1*legs[l][m]*
	  ( expcoef[loffset+moffset]   * potd[l]*cosm + 
	    expcoef[loffset+moffset+1] * potd[l]*sinm );
	dpot += fac1*legs[l][m]*
	  ( expcoef[loffset+moffset]   * dpot[l]*cosm + 
	    expcoef[loffset+moffset+1] * dpot[l]*sinm );
	pott += fac1*dlegs[l][m]*
	  ( expcoef[loffset+moffset]   * potd[l]*cosm + 
	    expcoef[loffset+moffset+1] * potd[l]*sinm );
	potp += fac1*legs[l][m] * m *
	  (-expcoef[loffset+moffset]   * potd[l]*sinm + 
	    expcoef[loffset+moffset+1] * potd[l]*cosm );

	moffset +=2;
      }
    }
  }

  // double densfac = 1.0/(rscl*rscl*rscl) * 0.25/M_PI;
  double potlfac = 1.0/rscl;

  potl  *= potlfac;
  potr  *= potlfac/rscl;
  pott  *= potlfac*sinth;
  potp  *= potlfac;
}


void SphereSL::all_eval(double r, double costh, double phi,
			double& den0, double& den1,
			double& pot0, double& pot1,
			double& potr, double& pott, double& potp,
			int L1, int L2)
{
  double fac1, cosm, sinm;
  double sinth = -sqrt(fabs(1.0 - costh*costh));

  fac1 = factorial[0][0];

  sl->get_dens (dend, r/rscl);
  sl->get_pot  (potd, r/rscl);
  sl->get_force(dpot, r/rscl);

  legendre_R(lmax, costh, legs, dlegs);

  den0 = fac1 * expcoef[0]*dend[0];
  pot0 = fac1 * expcoef[0]*potd[0];
  potr = fac1 * expcoef[0]*dpot[0];
  den1 = 0.0;
  pot1 = 0.0;
  pott = 0.0;
  potp = 0.0;

  // L loop
  for (int l=1, loffset=1; l<=lmax; loffset+=(2*l+1), l++) {
    if (l<L1 || l>L2) continue;

    // M loop
    for (int m=0, moffset=0; m<=l; m++) {
      fac1 = factorial[l][m];
      if (m==0) {
	den1 += fac1*legs[1][m] * (expcoef[loffset+moffset] * dend[l]);
	pot1 += fac1*legs[l][m] * (expcoef[loffset+moffset] * potd[l]);
	dpot += fac1*legs[l][m] * (expcoef[loffset+moffset] * dpot[l]);
	pott += fac1*dlegs[l][m]* (expcoef[loffset+moffset] * potd[l]);

	moffset++;
      }
      else {
	cosm = cos(phi*m);
	sinm = sin(phi*m);

	den1 += fac1*legs[l][m]*
	  ( expcoef[loffset+moffset]   * dend[l]*cosm + 
	    expcoef[loffset+moffset+1] * dend[l]*sinm );
	pot1 += fac1*legs[l][m]*
	  ( expcoef[loffset+moffset]   * potd[l]*cosm + 
	    expcoef[loffset+moffset+1] * potd[l]*sinm );
	dpot += fac1*legs[l][m]*
	  ( expcoef[loffset+moffset]   * dpot[l]*cosm + 
	    expcoef[loffset+moffset+1] * dpot[l]*sinm );
	pott += fac1*dlegs[l][m]*
	  ( expcoef[loffset+moffset]   * potd[l]*cosm + 
	    expcoef[loffset+moffset+1] * potd[l]*sinm );
	potp += fac1*legs[l][m] * m *
	  (-expcoef[loffset+moffset]   * potd[l]*sinm + 
	    expcoef[loffset+moffset+1] * potd[l]*cosm );

	moffset +=2;
      }
    }
  }

  double densfac = 1.0/(rscl*rscl*rscl) * 0.25/M_PI;
  double potlfac = 1.0/rscl;

  den0  *= densfac;
  den1  *= densfac;
  pot0  *= potlfac;
  pot1  *= potlfac;
  potr  *= potlfac/rscl;
  pott  *= potlfac*sinth;
  potp  *= potlfac;
}


#define MINEPS 1.0e-10

void SphereSL::legendre_R(int lmax, double x, Matrix& p)
{
  double fact, somx2, pll, pl1, pl2;
  int m, l;
  
  p[0][0] = pll = 1.0;
  if (lmax > 0) {
    somx2 = sqrt( (1.0 - x)*(1.0 + x) );
    fact = 1.0;
    for (m=1; m<=lmax; m++) {
      pll *= -fact*somx2;
      p[m][m] = pll;
      if (std::isnan(p[m][m]))
	cerr << "legendre_R: p[" << m << "][" << m << "]: pll=" << pll << "\n";
      fact += 2.0;
    }
  }
  
  for (m=0; m<lmax; m++) {
    pl2 = p[m][m];
    p[m+1][m] = pl1 = x*(2*m+1)*pl2;
    for (l=m+2; l<=lmax; l++) {
      p[l][m] = pll = (x*(2*l-1)*pl1-(l+m-1)*pl2)/(l-m);
      if (std::isnan(p[l][m]))
	cerr << "legendre_R: p[" << l << "][" << m << "]: pll=" << pll << "\n";
      
      pl2 = pl1;
      pl1 = pll;
    }
  }
  
  if (std::isnan(x))
    cerr << "legendre_R: x\n";
  for(l=0; l<=lmax; l++)
    for (m=0; m<=l; m++)
      if (std::isnan(p[l][m]))
	cerr << "legendre_R: p[" << l << "][" << m << "] lmax=" 
	     << lmax << "\n";
  
}

void SphereSL::legendre_R(int lmax, double x, Matrix &p, Matrix &dp)
{
  double fact, somx2, pll, pl1, pl2;
  int m, l;
  
  p[0][0] = pll = 1.0;
  if (lmax > 0) {
    somx2 = sqrt( (1.0 - x)*(1.0 + x) );
    fact = 1.0;
    for (m=1; m<=lmax; m++) {
      pll *= -fact*somx2;
      p[m][m] = pll;
      fact += 2.0;
    }
  }
  
  for (m=0; m<lmax; m++) {
    pl2 = p[m][m];
    p[m+1][m] = pl1 = x*(2*m+1)*pl2;
    for (l=m+2; l<=lmax; l++) {
      p[l][m] = pll = (x*(2*l-1)*pl1-(l+m-1)*pl2)/(l-m);
      pl2 = pl1;
      pl1 = pll;
    }
  }
  
  if (1.0-fabs(x) < MINEPS) {
    if (x>0) x =   1.0 - MINEPS;
    else     x = -(1.0 - MINEPS);
  }
  
  somx2 = 1.0/(x*x - 1.0);
  dp[0][0] = 0.0;
  for (l=1; l<=lmax; l++) {
    for (m=0; m<l; m++)
      dp[l][m] = somx2*(x*l*p[l][m] - (l+m)*p[l-1][m]);
    dp[l][l] = somx2*x*l*p[l][l];
  }
}

void SphereSL::legendre_R(int lmax, double x, Matrix &p, Matrix &dp, Matrix& d2p)
{
  double fact, somx2, pll, pl1, pl2;
  int m, l;
  
  p[0][0] = pll = 1.0;
  if (lmax > 0) {
    somx2 = sqrt( (1.0 - x)*(1.0 + x) );
    fact = 1.0;
    for (m=1; m<=lmax; m++) {
      pll *= -fact*somx2;
      p[m][m] = pll;
      fact += 2.0;
    }
  }
  
  for (m=0; m<lmax; m++) {
    pl2 = p[m][m];
    p[m+1][m] = pl1 = x*(2*m+1)*pl2;
    for (l=m+2; l<=lmax; l++) {
      p[l][m] = pll = (x*(2*l-1)*pl1-(l+m-1)*pl2)/(l-m);
      pl2 = pl1;
      pl1 = pll;
    }
  }
  
  if (1.0-fabs(x) < MINEPS) {
    if (x>0) x =   1.0 - MINEPS;
    else     x = -(1.0 - MINEPS);
  }
  
  somx2 = 1.0/(x*x - 1.0);
  dp[0][0] = 0.0;
  if (lmax) {
    for (l=1; l<=lmax; l++) {
      for (m=0; m<l; m++)
	dp[l][m] = somx2*(x*l*p[l][m] - (l+m)*p[l-1][m]);
      dp[l][l] = somx2*x*l*p[l][l];
    }
  }
  
  for (l=0; l<=lmax; l++) {
    for (m=0; m<=l; m++)
      d2p[l][m] = -somx2*(2.0*x*dp[l][m] - p[l][m]*(somx2*m*m + l*(l+1)));
  }
  
}


void SphereSL::install_coefs(Matrix& newcoef)
{
  if (!coefs_defined) {

    coefs_defined = true;
    
    expcoef.setsize(0, lmax*(lmax+2), 1, nmax);
    expcoef.zero();		// Need this?

    work.setsize(1, nmax);

    if (compute_covar) {
      covar.resize(lmax+1);
      mean.resize(lmax+1);
      for (int L=0; L<=lmax; L++) {
	int size = (L+1)*nmax;
	covar[L] = Eigen::MatrixXd::Zero(size, size);
	mean[L]  = Eigen::VectorXd::Zero(size);
      }
    }

    factorial.setsize(0, lmax, 0, lmax);

    for (int l=0; l<=lmax; l++) {
      for (int m=0; m<=l; m++) {
	factorial[l][m] = sqrt( (0.5*l+0.25)/M_PI * 
				exp(lgamma(1.0+l-m) - lgamma(1.0+l+m)) );
	if (m != 0) factorial[l][m] *= M_SQRT2;
      }
    }

    used = 0;
  }

  // Sanity check
  if (newcoef.getrhigh() != expcoef.getrhigh() ||
      newcoef.getrlow()  != expcoef.getrlow()  ||
      newcoef.getchigh() != expcoef.getchigh() ||
      newcoef.getclow()  != expcoef.getclow()  )
    {
      cerr << "SphereSL: can not install coefficients, dimension mismatch\n";
      return;
    }

  // Do the assignment
  expcoef = newcoef;
}


void SphereSL::dump_coefs(double time, ostream& out)
{
  ostringstream sout;
  sout << "SphereSL";

  char buf[64];
  for (int i=0; i<64; i++) {
    if (i<sout.str().length())  buf[i] = sout.str().c_str()[i];
    else                        buf[i] = '\0';
  }

  out.write((char *)&buf,64*sizeof(char));
  out.write((char *)&time , sizeof(double));
  out.write((char *)&rscl,  sizeof(double));
  out.write((char *)&nmax,  sizeof(int));
  out.write((char *)&lmax,  sizeof(int));

  for (int ir=1; ir<=nmax; ir++) {
    for (int l=0; l<=lmax*(lmax+2); l++)
      out.write((char *)&expcoef[l][ir], sizeof(double));
  }

}
