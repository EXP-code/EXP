#include <stdlib.h>
#include <math.h>
#include <iostream.h>
#include <iomanip.h>
#include <strstream.h>
#include <fstream.h>
#include <String.h>
#include <vector.h>

#define GHQL
#undef DEBUG

#include "expand.h"

extern "C" {

  void get_dens_coefs_bes(int, double *, double *);
  void get_pot_coefs_bes(int, double *, double *, double *);
  void get_potl_bes(int lmax, int nmax, double r, double **p);
  void get_potl_dens_bes(int lmax, int nmax, double r, double **p, double **d);

  void get_dens_coefs_CB(int, double *, double *);
  void get_pot_coefs_CB(int, double *, double *, double *);
  void get_potl_CB(int lmax, int nmax, double r, double **p);
  void get_potl_dens_CB(int lmax, int nmax, double r, double **p, double **d);

  void get_dens_coefs_HERNQ(int, double *, double *);
  void get_pot_coefs_HERNQ(int, double *, double *, double *);
  void get_potl_HERNQ(int lmax, int nmax, double r, double **p);
  void get_potl_dens_HERNQ(int lmax, int nmax, double r, double **p, double **d);

  void get_dens_coefs_CBDisk(int, double *, double *);
  void get_pot_coefs_CBDisk(int, double *, double *, double *);
  void get_potl_CBDisk(int lmax, int nmax, double r, double **p);
  void get_potl_dens_CBDisk(int lmax, int nmax, double r, double **p, double **d);

  void get_potl_dens_SLsph(int l, int n, double r);
  void get_dens_coefs_SLsph(int l, Vector& coef, double *p);
  void get_pot_coefs_SLsph(int l, Vector& coef, double *p, double *dp);

  extern double **potd, **dend;

  double return_elapsed_time(double last_time=0.0);
}


const double TOLR=1.0e-3;
const int NUM=100;

static void zeroth_order(String name)
{
  static int done = 0;
  if (done) return;
  done = 1;

  String out1 = name + ".den0";
  String out2 = name + ".pot0";

  double r, dr, d, dd;

  ofstream cout1(out1);
  if (!cout1) {
    cerr << "Couldn't open <" << out1 << ">\n";
    exit(-1);
  }
  
  ofstream cout2(out2);
  if (!cout2) {
    cerr << "Couldn't open <" << out2 << ">\n";
    exit(-1);
  }
  
  cout1.precision(6);
  cout1.setf(ios::scientific, ios::floatfield);

  cout2.precision(6);
  cout2.setf(ios::scientific, ios::floatfield);


  dr = (rdiag-TOLR)/(double)NUM;
  
  Vector flag(1, nmax);
  flag.zero();

  for (int i=0; i<=NUM; i++) {
    r = TOLR + dr*i;
    cout1 << setw(14) << r;
    cout2 << setw(14) << r;
    
    if (c_brock)
      get_potl_dens_CB(lmax, nmax, r, potd, dend);
    else if (hernq)
      get_potl_dens_HERNQ(lmax, nmax, r, potd, dend);
    else if (bessel_sph)
      get_potl_dens_bes(lmax, nmax, r, potd, dend);
    else if (c_brock_disk)
      get_potl_dens_CBDisk(lmax, nmax, r, potd, dend);
    else if (sphereSL)
      get_potl_dens_SLsph(lmax, nmax, r);
    else {
      cerr << "tk: no expansion defined" << endl;
      _exit(-3);
    }


    for (int n=1; n<=nmax; n++) {
      flag[n] = 1.0;

      if (c_brock)
	get_dens_coefs_CB(L_pca, flag.array(1, nmax), &d);
      else if (hernq)
	get_dens_coefs_HERNQ(L_pca, flag.array(1, nmax), &d);
      else if (bessel_sph)
	get_dens_coefs_bes(L_pca, flag.array(1, nmax), &d);
      else if (c_brock_disk)
	get_dens_coefs_CBDisk(M_pca, flag.array(1, nmax), &d);
      else if (sphereSL)
	get_dens_coefs_SLsph(L_pca, flag, &d);
      else {
	cerr << "tk: no expansion defined" << endl;
	_exit(-3);
      }


      cout1 << setw(14) << d;

      if (c_brock)
	get_pot_coefs_CB(L_pca, flag.array(1, nmax), &d, &dd);
      else if (hernq)
	get_pot_coefs_HERNQ(L_pca, flag.array(1, nmax), &d, &dd);
      else if (bessel_sph)
	get_pot_coefs_bes(L_pca, flag.array(1, nmax), &d, &dd);
      else if (c_brock_disk)
	get_pot_coefs_CBDisk(L_pca, flag.array(1, nmax), &d, &dd);
      else if (sphereSL)
	get_pot_coefs_SLsph(L_pca, flag, &d, &dd);
      else {
	cerr << "tk: no expansion defined" << endl;
	_exit(-3);
      }

      cout2 << setw(14) << d;
      flag[n] = 0.0;
    }
    cout1 << endl;
    cout2 << endl;
  }
}

//
// The switching based of "dof" in these routines is so that one set
// of routines is valid for both 2d polar and 3d spherical polar
// orthogonal functions
//


extern "C" void out_pca(int k)
{
  if (!selector) return;	// Fail safe . . .

  ostrstream number;
  number << k << '\0';
  String curname = (String)outname + "." + number.str();

  zeroth_order(curname);

  String out1 = curname + "_pca.den";
  String out2 = curname + "_pca.pot";

				/* Print out stuff */
  ofstream cout1(out1);
  if (!cout1) {
    cerr << "Couldn't open <" << out1 << ">\n";
    exit(-1);
  }
  
  ofstream cout2(out2);
  if (!cout2) {
    cerr << "Couldn't open <" << out2 << ">\n";
    exit(-1);
  }
  
  cout1.precision(6);
  cout1.setf(ios::scientific, ios::floatfield);

  cout2.precision(6);
  cout2.setf(ios::scientific, ios::floatfield);


				// So one can test smoothing without epplying
  double smooth = tksmooth;
  if (smooth<=0.0) smooth=1.0;

  double r, dr, d, dd, fac, fac02, var;

  int indx=0;
  const double 
  start_timer();
				/* Compute variance matrix */

  int lpca;
  if (dof==2) {
    lpca = 0;
    fac02 = 1.0;
  }
  else {
    lpca = L_pca;
    fac02 = 16.0*M_PI*M_PI;
  }

  for (int l=0; l<lpca; l++) indx += 2*l+1;
  for (int m=1; m<M_pca; m++) indx += 2;

  Vector inv(1, nmax);
  Vector smth(1, nmax);
  Vector eval(1, nmax);
  Vector cuml(1, nmax);
  Matrix evec(1, nmax, 1, nmax);
  Matrix Tevec(1, nmax, 1, nmax);

  Matrix covar(1, nmax, 1, nmax);

  Vector normed(1, nmax);

  int n, nn;
  double b;

				// Both cos and sin for M>0
  int imult = M_pca>0 ? 2 : 1;
  for (int im=0; im<imult; im++) {

    for(n=1; n<=nmax; n++) {
      b = (cc[indx][n][n]*fac02 - expcoef[indx][n]*expcoef[indx][n]) /
	(expcoef[indx][n]*expcoef[indx][n]*used);
      expcoef[indx][n] *= 1.0/(1.0 + b);
    }
    
    for(n=1; n<=nmax; n++) {
      for(nn=n; nn<=nmax; nn++) {
	fac = sqrt(normM[L_pca][n]*normM[L_pca][nn]);
	covar[n][nn] = fac * expcoef[indx][n]*expcoef[indx][nn];
	
	if (n!=nn)
	  covar[nn][n] = covar[n][nn];
      }    
    }
				/* Diagonalize variance */


#ifdef GHQL
    eval = covar.Symmetric_Eigenvalues_GHQL(evec);
#else
    eval = covar.Symmetric_Eigenvalues(evec);
#endif
    Tevec = evec.Transpose();

    cout1 << "#\n#\n# Eigenvectors for density computation: (" << L_pca 
      << ", " << M_pca << ")";
    cout2 << "#\n#\n# Eigenvectors for potential computation: (" << L_pca 
      << ", " << M_pca << ")";
    if (M_pca>0) {
      if (im==0) {
	cout1 << "  Cosine series\n#\n";  
	cout2 << "  Cosine series\n#\n";  
      }
      if (im==1) {
	cout1 << "  Sine series\n#\n";  
	cout2 << "  Sine series\n#\n";  
      }
    }
    else {
      cout1 << "\n#\n";  
      cout2 << "\n#\n";  
    }

    cout1 << "#" << setw(3) << "N" << setw(15) << "Covariance" 
      << setw(15) << "Coefficient" << setw(15) << "Coef^2" 
	<< setw(15) << "Variance" << setw(15) << "Smooth" << "\n#\n";

    cout2 << "#" << setw(3) << "N" << setw(15) << "Covariance" 
      << setw(15) << "Coefficient" << setw(15) << "Coef^2" 
	<< setw(15) << "Variance" << setw(15) << "Smooth" << "\n#\n";

    for (n=1; n<=nmax; n++) {

      for (dd=0.0, nn=1; nn<=nmax; nn++) dd += Tevec[n][nn]*expcoef[indx][nn]
	* sqrt(normM[L_pca][nn]);

      var = eval[n]/used - dd*dd;

      fac = 1.0/(1.0 + smooth*var/(dd*dd));
      smth[n] = dd*fac;
      cout1 << "# " << setw(4) << n << setw(15) << eval[n] << setw(15) << dd
	    << setw(15) << dd*dd << setw(15) << var
	    << setw(15) << fac << endl;
    }

    inv = evec*smth;
    cout1 << "#" << endl;
    for (n=1; n<=nmax; n++)
      cout1 << "# " << setw(4) << n << setw(15) << expcoef[indx][n] 
	<< setw(15) << inv[n]/sqrt(normM[L_pca][n]) << endl;
    cout1 << "#" << endl;


    dr = (rdiag-TOLR)/(double)NUM;
  
    for (int i=0; i<=NUM; i++) {
      r = TOLR + dr*i;
      cout1 << setw(14) << r;
      cout2 << setw(14) << r;

      if (c_brock)
	get_potl_dens_CB(lmax, nmax, r, potd, dend);
      else if (hernq)
	get_potl_dens_HERNQ(lmax, nmax, r, potd, dend);
      else if (bessel_sph)
	get_potl_dens_bes(lmax, nmax, r, potd, dend);
      else if (c_brock_disk)
	get_potl_dens_CBDisk(lmax, nmax, r, potd, dend);
      else if (sphereSL)
	get_potl_dens_SLsph(lmax, nmax, r);
      else {
	cerr << "tk: no expansion defined" << endl;
	_exit(-3);
      }

      for (n=1; n<=nmax; n++) {

	for (nn=1; nn<=nmax; nn++) 
	  normed[nn] = Tevec[n][nn]/sqrt(normM[L_pca][nn]);

	if (c_brock)
	  get_dens_coefs_CB(L_pca, normed.array(1, nmax), &d);
	else if (hernq)
	  get_dens_coefs_HERNQ(L_pca, normed.array(1, nmax), &d);
	else if (bessel_sph)
	  get_dens_coefs_bes(L_pca, normed.array(1, nmax), &d);
	else if (c_brock_disk)
	  get_dens_coefs_bes(M_pca, normed.array(1, nmax), &d);
	else if (sphereSL)
	  get_dens_coefs_SLsph(M_pca, normed, &d);
	else {
	  cerr << "tk: no expansion defined" << endl;
	  _exit(-3);
	}

	cout1 << setw(14) << d;

	if (c_brock)
	  get_pot_coefs_CB(L_pca, normed.array(1, nmax), &d, &dd);
	else if (hernq)
	  get_pot_coefs_HERNQ(L_pca, normed.array(1, nmax), &d, &dd);
	else if (bessel_sph)
	  get_pot_coefs_bes(L_pca, normed.array(1, nmax), &d, &dd);
	else if (c_brock_disk)
	  get_pot_coefs_bes(M_pca, normed.array(1, nmax), &d, &dd);
	else if (sphereSL)
	  get_pot_coefs_SLsph(M_pca, normed, &d, &dd);
	else {
	  cerr << "tk: no expansion defined" << endl;
	  _exit(-3);
	}

	cout2 << setw(14) << d;
      }
      cout1 << endl;
      cout2 << endl;
    }

    indx++;

  }

  cout1 << "# out_pca cpu time=" << return_elapsed_time() << endl;
}

extern "C" void pca_hall(int compute)
{
#ifdef DEBUG
  start_timer();
#endif

  int l, Ldim;
  double dd, fac, fac02, var, b;

  static Vector smth;
  static Vector *weight;
  static Vector *b_Hall;
  static Vector inv;
  static Vector eval;
  static Vector cuml;
  static Matrix *evec;
  static Matrix Tevec;
  static Matrix sqnorm;

  static Matrix covar;

  static int setup = 0;
  
  if (!setup) {

    if (dof==3)
      Ldim = lmax*(lmax + 2) + 1;
    else
      Ldim = 2*lmax + 1;
      
    weight = new Vector [Ldim];
    b_Hall = new Vector [Ldim];
    evec = new Matrix [Ldim];
    
    for (l=0; l<Ldim; l++) {
      weight[l].setsize(1, nmax);
      b_Hall[l].setsize(1, nmax);
      evec[l].setsize(1, nmax, 1, nmax);
    }

    smth.setsize(1, nmax);
    inv.setsize(1, nmax);
    eval.setsize(1, nmax);
    cuml.setsize(1, nmax);
    Tevec.setsize(1, nmax, 1, nmax);
    covar.setsize(1, nmax, 1, nmax);
    sqnorm.setsize(0, lmax, 1, nmax);
      
    for (l=0; l<=lmax; l++)
      for (int n=1; n<=nmax; n++) sqnorm[l][n] = sqrt(normM[l][n]);

    setup = 1;
  }


  int m, loffset, moffset, n, nn, indx, L0, lm;

  if (dof==3) {
    L0 = 0;
    fac02 = 16.0*M_PI*M_PI;
  }
  else {
    L0 = lmax;
    fac02 = 1.0;
  }

  for (l=L0, loffset=0; l<=lmax; loffset+=(2*l+1), l++) {

    for (m=0, moffset=0; m<=l; m++) {

      if (dof==3) {
	lm = l;
	indx = loffset+moffset;
      }
      else {
	lm = m;
	indx = moffset;
      }

      if (m==0) {

	if (compute) {

	  for(n=1; n<=nmax; n++) {
	    b = (cc[indx][n][n]*fac02 - expcoef[indx][n]*expcoef[indx][n]) /
	      (expcoef[indx][n]*expcoef[indx][n]*used);
	    b_Hall[indx][n] = 1.0/(1.0 + b);
	  }
    
	  for(n=1; n<=nmax; n++) {
	    for(nn=n; nn<=nmax; nn++) {
	      fac = sqnorm[lm][n]*sqnorm[lm][nn];
	      covar[n][nn] = fac * expcoef[indx][n]*expcoef[indx][nn];
	      if (n!=nn)
		covar[nn][n] = covar[n][nn];
	    }    
	  }

				/* Diagonalize variance */

#ifdef GHQL
	  eval = covar.Symmetric_Eigenvalues_GHQL(evec[indx]);
#else
	  eval = covar.Symmetric_Eigenvalues(evec[indx]);
#endif
	  Tevec = evec[indx].Transpose();

	  if (tk_type == 2) {
	    cuml = eval;
	    for (n=2; n<=nmax; n++) cuml[n] += cuml[n-1];
	    var = cuml[nmax];
	    for (n=1; n<=nmax; n++) cuml[n] /= var;
	  }

	  for (n=1; n<=nmax; n++) {

	    for (dd=0.0, nn=1; nn<=nmax; nn++) 
	      dd += Tevec[n][nn]*expcoef[indx][nn]*sqnorm[lm][nn];

	    var = eval[n]/used - dd*dd;

	    if (tk_type == 1) {

	      if (tksmooth*var > dd*dd)
		weight[indx][n] = 0.0;
	      else
		weight[indx][n] = 1.0;

	    }
	    else if (tk_type == 2) {
	      
	      if (n==1 || cuml[n] <= tkcum)
		weight[indx][n] = 1.0;
	      else
		weight[indx][n] = 0.0;
		
	    }
	    else if (tk_type == 3) {
	      
	      weight[indx][n] = 1.0/(1.0 + var/(dd*dd + 1.0e-14));
		
	    }
	    else
		weight[indx][n] = 1.0;

	    smth[n] = dd * weight[indx][n];
	  }

	}
	else {
	  Tevec = evec[indx].Transpose();
	  for (n=1; n<=nmax; n++) {
	    for (dd=0.0, nn=1; nn<=nmax; nn++) 
	      dd += Tevec[n][nn]*expcoef[indx][nn] * sqnorm[lm][nn];
	    smth[n] = dd * weight[indx][n];
	  }
	}
	    
	inv = evec[indx]*smth;
	for (n=1; n<=nmax; n++) {
	  expcoef[indx][n] = inv[n]/sqnorm[lm][n];
	  if (tk_type == 0) expcoef[indx][n] *= b_Hall[indx][n];
	}
	
	moffset++;
      }
      else {

	if (compute) {

	  for(n=1; n<=nmax; n++) {
	    b = (cc[indx][n][n]*fac02 - expcoef[indx][n]*expcoef[indx][n]) /
	      (expcoef[indx][n]*expcoef[indx][n]*used);
	    b_Hall[indx][n] = 1.0/(1.0 + b);
	  }
    
	  for(n=1; n<=nmax; n++) {
	    for(nn=n; nn<=nmax; nn++) {
	      fac = sqnorm[lm][n] * sqnorm[lm][nn];
	      covar[n][nn] = fac * expcoef[indx][n]*expcoef[indx][nn];
	      if (n!=nn)
		covar[nn][n] = covar[n][nn];
	    }
	  }  

				/* Diagonalize variance */

#ifdef GHQL
	  eval = covar.Symmetric_Eigenvalues_GHQL(evec[indx]);
#else
	  eval = covar.Symmetric_Eigenvalues(evec[indx]);
#endif
	  Tevec = evec[indx].Transpose();

	  if (tk_type == 2) {
	    cuml = eval;
	    for (n=2; n<=nmax; n++) cuml[n] += cuml[n-1];
	    var = cuml[nmax];
	    for (n=1; n<=nmax; n++) cuml[n] /= var;
	  }

	  for (n=1; n<=nmax; n++) {

	    for (dd=0.0, nn=1; nn<=nmax; nn++) 
	      dd += Tevec[n][nn]*expcoef[indx][nn]*sqnorm[lm][nn];

	    var = eval[n]/used - dd*dd;

	    if (tk_type == 1) {

	      if (tksmooth*var > dd*dd)
		weight[indx][n] = 0.0;
	      else
		weight[indx][n] = 1.0;

	    }
	    else if (tk_type == 2) {
	      
	      if (n==1 || cuml[n] <= tkcum)
		weight[indx][n] = 1.0;
	      else
		weight[indx][n] = 0.0;
		
	    }
	    else if (tk_type == 3) {
	      
	      weight[indx][n] = 1.0/(1.0 + var/(dd*dd + 1.0e-14));
		
	    }
	    else
		weight[indx][n] = 1.0;

	    smth[n] = dd * weight[indx][n];
	  }
	}
	else {
	  Tevec = evec[indx].Transpose();
	  for (n=1; n<=nmax; n++) {
	    for (dd=0.0, nn=1; nn<=nmax; nn++) 
	      dd += Tevec[n][nn]*expcoef[indx][nn] * sqnorm[lm][nn];
	    smth[n] = dd * weight[indx][n];
	  }
	}
	    
	inv = evec[indx]*smth;
	for (n=1; n<=nmax; n++) {
	  expcoef[indx][n] = inv[n]/sqnorm[lm][n];
	  if (tk_type == 0) expcoef[indx][n] *= b_Hall[indx][n];
	}
	
	indx++;

	if (compute) {

	  for(n=1; n<=nmax; n++) {
	    b = (cc[indx][n][n]*fac02 - expcoef[indx][n]*expcoef[indx][n]) /
	      (expcoef[indx][n]*expcoef[indx][n]*used);
	    b_Hall[indx][n] = 1.0/(1.0 + b);
	  }
    
	  for(n=1; n<=nmax; n++) {
	    for(nn=n; nn<=nmax; nn++) {
	      fac = sqnorm[lm][n] * sqnorm[lm][nn];
	      covar[n][nn] = fac * expcoef[indx][n]*expcoef[indx][nn];
	      if (n!=nn)
		covar[nn][n] = covar[n][nn];
	    }    
	  }

				/* Diagonalize variance */

#ifdef GHQL
	  eval = covar.Symmetric_Eigenvalues_GHQL(evec[indx]);
#else
	  eval = covar.Symmetric_Eigenvalues(evec[indx]);
#endif
	  Tevec = evec[indx].Transpose();

	  if (tk_type == 2) {
	    cuml = eval;
	    for (n=2; n<=nmax; n++) cuml[n] += cuml[n-1];
	    var = cuml[nmax];
	    for (n=1; n<=nmax; n++) cuml[n] /= var;
	  }

	  for (n=1; n<=nmax; n++) {

	    for (dd=0.0, nn=1; nn<=nmax; nn++) 
	      dd += Tevec[n][nn]*expcoef[indx][nn]*sqnorm[lm][nn];

	    var = eval[n]/used - dd*dd;

	    if (tk_type == 1) {

	      if (tksmooth*var > dd*dd)
		weight[indx][n] = 0.0;
	      else
		weight[indx][n] = 1.0;

	    }
	    else if (tk_type == 2) {
	      
	      if (n==1 || cuml[n] <= tkcum)
		weight[indx][n] = 1.0;
	      else
		weight[indx][n] = 0.0;
		
	    }
	    else if (tk_type == 3) {
	      
	      weight[indx][n] = 1.0/(1.0 + var/(dd*dd + 1.0e-14));
		
	    }
	    else
		weight[indx][n] = 1.0;

	    smth[n] = dd * weight[indx][n];
	  }
	}
	else {
	  Tevec = evec[indx].Transpose();
	  for (n=1; n<=nmax; n++) {
	    for (dd=0.0, nn=1; nn<=nmax; nn++) 
	      dd += Tevec[n][nn]*expcoef[indx][nn] * sqnorm[lm][nn];
	    smth[n] = dd * weight[indx][n];
	  }
	}

	inv = evec[indx]*smth;
	for (n=1; n<=nmax; n++) {
	  expcoef[indx][n] = inv[n]/sqnorm[lm][n];
	  if (tk_type == 0) expcoef[indx][n] *= b_Hall[indx][n];
	}
	
	moffset += 2;
      }
    }
  }

#ifdef DEBUG
  cerr << "pca_hall cpu_time=" << return_elapsed_time() << endl;
#endif

}

void pca_hall(int compute, Matrix& expcoef, Matrix*& cc, Matrix& normM)
{
#ifdef DEBUG
  start_timer();
#endif

  int l, Ldim;
  double dd, fac, fac02, var, b;

  static Vector smth;
  static Vector *weight;
  static Vector *b_Hall;
  static Vector inv;
  static Vector eval;
  static Vector cuml;
  static Matrix *evec;
  static Matrix Tevec;
  static Matrix sqnorm;

  static Matrix covar;

  static int setup = 0;
  
  if (!setup) {

    if (dof==3)
      Ldim = lmax*(lmax + 2) + 1;
    else
      Ldim = 2*lmax + 1;
      
    weight = new Vector [Ldim];
    b_Hall = new Vector [Ldim];
    evec = new Matrix [Ldim];
    
    for (l=0; l<Ldim; l++) {
      weight[l].setsize(1, nmax);
      b_Hall[l].setsize(1, nmax);
      evec[l].setsize(1, nmax, 1, nmax);
    }

    smth.setsize(1, nmax);
    inv.setsize(1, nmax);
    eval.setsize(1, nmax);
    cuml.setsize(1, nmax);
    Tevec.setsize(1, nmax, 1, nmax);
    covar.setsize(1, nmax, 1, nmax);
    sqnorm.setsize(0, lmax, 1, nmax);
      
    for (l=0; l<=lmax; l++)
      for (int n=1; n<=nmax; n++) sqnorm[l][n] = sqrt(normM[l][n]);

    setup = 1;
  }


  int m, loffset, moffset, n, nn, indx, L0, lm;

  if (dof==3) {
    L0 = 0;
    fac02 = 16.0*M_PI*M_PI;
  }
  else {
    L0 = lmax;
    fac02 = 1.0;
  }

  for (l=L0, loffset=0; l<=lmax; loffset+=(2*l+1), l++) {

    for (m=0, moffset=0; m<=l; m++) {

      if (dof==3) {
	lm = l;
	indx = loffset+moffset;
      }
      else {
	lm = m;
	indx = moffset;
      }

      if (m==0) {

	if (compute) {

	  for(n=1; n<=nmax; n++) {
	    b = (cc[indx][n][n]*fac02 - expcoef[indx][n]*expcoef[indx][n]) /
	      (expcoef[indx][n]*expcoef[indx][n]*used);
	    b_Hall[indx][n] = 1.0/(1.0 + b);
	  }
    
	  for(n=1; n<=nmax; n++) {
	    for(nn=n; nn<=nmax; nn++) {
	      fac = sqnorm[lm][n]*sqnorm[lm][nn];
	      covar[n][nn] = fac * expcoef[indx][n]*expcoef[indx][nn];
	      if (n!=nn)
		covar[nn][n] = covar[n][nn];
	    }    
	  }

				/* Diagonalize variance */

#ifdef GHQL
	  eval = covar.Symmetric_Eigenvalues_GHQL(evec[indx]);
#else
	  eval = covar.Symmetric_Eigenvalues(evec[indx]);
#endif
	  Tevec = evec[indx].Transpose();

	  if (tk_type == 2) {
	    cuml = eval;
	    for (n=2; n<=nmax; n++) cuml[n] += cuml[n-1];
	    var = cuml[nmax];
	    for (n=1; n<=nmax; n++) cuml[n] /= var;
	  }

	  for (n=1; n<=nmax; n++) {

	    for (dd=0.0, nn=1; nn<=nmax; nn++) 
	      dd += Tevec[n][nn]*expcoef[indx][nn]*sqnorm[lm][nn];

	    var = eval[n]/used - dd*dd;

	    if (tk_type == 1) {

	      if (tksmooth*var > dd*dd)
		weight[indx][n] = 0.0;
	      else
		weight[indx][n] = 1.0;

	    }
	    else if (tk_type == 2) {
	      
	      if (n==1 || cuml[n] <= tkcum)
		weight[indx][n] = 1.0;
	      else
		weight[indx][n] = 0.0;
		
	    }
	    else if (tk_type == 3) {
	      
	      weight[indx][n] = 1.0/(1.0 + var/(dd*dd + 1.0e-14));
		
	    }
	    else
		weight[indx][n] = 1.0;

	    smth[n] = dd * weight[indx][n];
	  }

	}
	else {
	  Tevec = evec[indx].Transpose();
	  for (n=1; n<=nmax; n++) {
	    for (dd=0.0, nn=1; nn<=nmax; nn++) 
	      dd += Tevec[n][nn]*expcoef[indx][nn] * sqnorm[lm][nn];
	    smth[n] = dd * weight[indx][n];
	  }
	}
	    
	inv = evec[indx]*smth;
	for (n=1; n<=nmax; n++) {
	  expcoef[indx][n] = inv[n]/sqnorm[lm][n];
	  if (tk_type == 0) expcoef[indx][n] *= b_Hall[indx][n];
	}
	
	moffset++;
      }
      else {

	if (compute) {

	  for(n=1; n<=nmax; n++) {
	    b = (cc[indx][n][n]*fac02 - expcoef[indx][n]*expcoef[indx][n]) /
	      (expcoef[indx][n]*expcoef[indx][n]*used);
	    b_Hall[indx][n] = 1.0/(1.0 + b);
	  }
    
	  for(n=1; n<=nmax; n++) {
	    for(nn=n; nn<=nmax; nn++) {
	      fac = sqnorm[lm][n] * sqnorm[lm][nn];
	      covar[n][nn] = fac * expcoef[indx][n]*expcoef[indx][nn];
	      if (n!=nn)
		covar[nn][n] = covar[n][nn];
	    }
	  }  

				/* Diagonalize variance */

#ifdef GHQL
	  eval = covar.Symmetric_Eigenvalues_GHQL(evec[indx]);
#else
	  eval = covar.Symmetric_Eigenvalues(evec[indx]);
#endif
	  Tevec = evec[indx].Transpose();

	  if (tk_type == 2) {
	    cuml = eval;
	    for (n=2; n<=nmax; n++) cuml[n] += cuml[n-1];
	    var = cuml[nmax];
	    for (n=1; n<=nmax; n++) cuml[n] /= var;
	  }

	  for (n=1; n<=nmax; n++) {

	    for (dd=0.0, nn=1; nn<=nmax; nn++) 
	      dd += Tevec[n][nn]*expcoef[indx][nn]*sqnorm[lm][nn];

	    var = eval[n]/used - dd*dd;

	    if (tk_type == 1) {

	      if (tksmooth*var > dd*dd)
		weight[indx][n] = 0.0;
	      else
		weight[indx][n] = 1.0;

	    }
	    else if (tk_type == 2) {
	      
	      if (n==1 || cuml[n] <= tkcum)
		weight[indx][n] = 1.0;
	      else
		weight[indx][n] = 0.0;
		
	    }
	    else if (tk_type == 3) {
	      
	      weight[indx][n] = 1.0/(1.0 + var/(dd*dd + 1.0e-14));
		
	    }
	    else
		weight[indx][n] = 1.0;

	    smth[n] = dd * weight[indx][n];
	  }
	}
	else {
	  Tevec = evec[indx].Transpose();
	  for (n=1; n<=nmax; n++) {
	    for (dd=0.0, nn=1; nn<=nmax; nn++) 
	      dd += Tevec[n][nn]*expcoef[indx][nn] * sqnorm[lm][nn];
	    smth[n] = dd * weight[indx][n];
	  }
	}
	    
	inv = evec[indx]*smth;
	for (n=1; n<=nmax; n++) {
	  expcoef[indx][n] = inv[n]/sqnorm[lm][n];
	  if (tk_type == 0) expcoef[indx][n] *= b_Hall[indx][n];
	}
	
	indx++;

	if (compute) {

	  for(n=1; n<=nmax; n++) {
	    b = (cc[indx][n][n]*fac02 - expcoef[indx][n]*expcoef[indx][n]) /
	      (expcoef[indx][n]*expcoef[indx][n]*used);
	    b_Hall[indx][n] = 1.0/(1.0 + b);
	  }
    
	  for(n=1; n<=nmax; n++) {
	    for(nn=n; nn<=nmax; nn++) {
	      fac = sqnorm[lm][n] * sqnorm[lm][nn];
	      covar[n][nn] = fac * expcoef[indx][n]*expcoef[indx][nn];
	      if (n!=nn)
		covar[nn][n] = covar[n][nn];
	    }    
	  }

				/* Diagonalize variance */

#ifdef GHQL
	  eval = covar.Symmetric_Eigenvalues_GHQL(evec[indx]);
#else
	  eval = covar.Symmetric_Eigenvalues(evec[indx]);
#endif
	  Tevec = evec[indx].Transpose();

	  if (tk_type == 2) {
	    cuml = eval;
	    for (n=2; n<=nmax; n++) cuml[n] += cuml[n-1];
	    var = cuml[nmax];
	    for (n=1; n<=nmax; n++) cuml[n] /= var;
	  }

	  for (n=1; n<=nmax; n++) {

	    for (dd=0.0, nn=1; nn<=nmax; nn++) 
	      dd += Tevec[n][nn]*expcoef[indx][nn]*sqnorm[lm][nn];

	    var = eval[n]/used - dd*dd;

	    if (tk_type == 1) {

	      if (tksmooth*var > dd*dd)
		weight[indx][n] = 0.0;
	      else
		weight[indx][n] = 1.0;

	    }
	    else if (tk_type == 2) {
	      
	      if (n==1 || cuml[n] <= tkcum)
		weight[indx][n] = 1.0;
	      else
		weight[indx][n] = 0.0;
		
	    }
	    else if (tk_type == 3) {
	      
	      weight[indx][n] = 1.0/(1.0 + var/(dd*dd + 1.0e-14));
		
	    }
	    else
		weight[indx][n] = 1.0;

	    smth[n] = dd * weight[indx][n];
	  }
	}
	else {
	  Tevec = evec[indx].Transpose();
	  for (n=1; n<=nmax; n++) {
	    for (dd=0.0, nn=1; nn<=nmax; nn++) 
	      dd += Tevec[n][nn]*expcoef[indx][nn] * sqnorm[lm][nn];
	    smth[n] = dd * weight[indx][n];
	  }
	}

	inv = evec[indx]*smth;
	for (n=1; n<=nmax; n++) {
	  expcoef[indx][n] = inv[n]/sqnorm[lm][n];
	  if (tk_type == 0) expcoef[indx][n] *= b_Hall[indx][n];
	}
	
	moffset += 2;
      }
    }
  }

#ifdef DEBUG
  cerr << "pca_hall cpu_time=" << return_elapsed_time() << endl;
#endif

}

