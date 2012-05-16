#include "expand.h"
#include <AxisymmetricBasis.H>

AxisymmetricBasis:: AxisymmetricBasis(string& line) : Basis(line) 
{
  Lmax = 4;
  nmax = 10;
  dof = 3;
  npca = 500;
  pca = false;
  tksmooth = 3.0;
  tkcum = 0.95;
  tk_type = Hall;

  string val;

  if (get_value("Lmax", val)) Lmax = atoi(val.c_str());
  if (get_value("nmax", val)) nmax = atoi(val.c_str());
  if (get_value("dof", val)) dof = atoi(val.c_str());
  if (get_value("npca", val)) npca = atoi(val.c_str());
  if (get_value("pca", val)) {
    if (atoi(val.c_str())) pca = true; 
    else pca = false;
  }
  if (get_value("tksmooth", val)) tksmooth = atof(val.c_str());
  if (get_value("tkcum", val)) tkcum = atof(val.c_str());
  if (get_value("tk_type", val)) {
    switch (atoi(val.c_str())) {
    case Hall:			tk_type = Hall; break;
    case VarianceCut:		tk_type = VarianceCut; break;
    case CumulativeCut:		tk_type = CumulativeCut; break;
    case VarianceWeighted:	tk_type = VarianceWeighted; break;
    default:
      if (myid==0) {
	cout << "AxisymmetricBasis: no such TK type <" << val << ">"
	     << " using Hall type\n";
      }
    }
  }

  if (pca) {

    if (dof==3)
      Ldim = Lmax*(Lmax + 2) + 1;
    else
      Ldim = 2*Lmax + 1;
    
    weight = new Vector [Ldim];
    b_Hall = new Vector [Ldim];
    evec = new Matrix [Ldim];
    
    for (int l=0; l<Ldim; l++) {
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
    sqnorm.setsize(0, Lmax, 1, nmax);
      
    for (int l=0; l<=Lmax; l++)
      for (int n=1; n<=nmax; n++) sqnorm[l][n] = 1.0;

  }
}

AxisymmetricBasis::~AxisymmetricBasis()
{
  vector<Matrix *>::iterator it;
  for (it=expcoefN.begin(); it!=expcoefN.end(); it++) delete *it;
  for (it=expcoefL.begin(); it!=expcoefL.end(); it++) delete *it;

  if (pca) {
    delete [] weight;
    delete [] b_Hall;
    delete [] evec;
  }
}


void AxisymmetricBasis::pca_hall(bool compute)
{
#ifdef DEBUG
  // start_timer();
#endif

  if (dof==3)
    Ldim = Lmax*(Lmax + 2) + 1;
  else
    Ldim = 2*Lmax + 1;
    
  if (dof==3) {
    L0 = 0;
    fac02 = 16.0*M_PI*M_PI;
  }
  else {
    L0 = Lmax;
    fac02 = 1.0;
  }

  double dd, fac, var, b;
  int loffset, moffset, indx, lm, nn;


  for (int l=L0, loffset=0; l<=Lmax; loffset+=(2*l+1), l++) {

    for (int m=0, moffset=0; m<=l; m++) {

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

	  for(int n=1; n<=nmax; n++) {
	    b = (cc[indx][n][n]*fac02 - expcoef[indx][n]*expcoef[indx][n]) /
	      (expcoef[indx][n]*expcoef[indx][n]*used);
	    b_Hall[indx][n] = 1.0/(1.0 + b);
	  }
    
	  for(int n=1; n<=nmax; n++) {
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
	    for (int n=2; n<=nmax; n++) cuml[n] += cuml[n-1];
	    var = cuml[nmax];
	    for (int n=1; n<=nmax; n++) cuml[n] /= var;
	  }

	  for (int n=1; n<=nmax; n++) {

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
	  for (int n=1; n<=nmax; n++) {
	    for (dd=0.0, nn=1; nn<=nmax; nn++) 
	      dd += Tevec[n][nn]*expcoef[indx][nn] * sqnorm[lm][nn];
	    smth[n] = dd * weight[indx][n];
	  }
	}
	    
	inv = evec[indx]*smth;
	for (int n=1; n<=nmax; n++) {
	  expcoef[indx][n] = inv[n]/sqnorm[lm][n];
	  if (tk_type == 0) expcoef[indx][n] *= b_Hall[indx][n];
	}
	
	moffset++;
      }
      else {

	if (compute) {

	  for(int n=1; n<=nmax; n++) {
	    b = (cc[indx][n][n]*fac02 - expcoef[indx][n]*expcoef[indx][n]) /
	      (expcoef[indx][n]*expcoef[indx][n]*used);
	    b_Hall[indx][n] = 1.0/(1.0 + b);
	  }
    
	  for(int n=1; n<=nmax; n++) {
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
	    for (int n=2; n<=nmax; n++) cuml[n] += cuml[n-1];
	    var = cuml[nmax];
	    for (int n=1; n<=nmax; n++) cuml[n] /= var;
	  }

	  for (int n=1; n<=nmax; n++) {

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
	  for (int n=1; n<=nmax; n++) {
	    for (dd=0.0, nn=1; nn<=nmax; nn++) 
	      dd += Tevec[n][nn]*expcoef[indx][nn] * sqnorm[lm][nn];
	    smth[n] = dd * weight[indx][n];
	  }
	}
	    
	inv = evec[indx]*smth;
	for (int n=1; n<=nmax; n++) {
	  expcoef[indx][n] = inv[n]/sqnorm[lm][n];
	  if (tk_type == 0) expcoef[indx][n] *= b_Hall[indx][n];
	}
	
	indx++;

	if (compute) {

	  for(int n=1; n<=nmax; n++) {
	    b = (cc[indx][n][n]*fac02 - expcoef[indx][n]*expcoef[indx][n]) /
	      (expcoef[indx][n]*expcoef[indx][n]*used);
	    b_Hall[indx][n] = 1.0/(1.0 + b);
	  }
    
	  for(int n=1; n<=nmax; n++) {
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
	    for (int n=2; n<=nmax; n++) cuml[n] += cuml[n-1];
	    var = cuml[nmax];
	    for (int n=1; n<=nmax; n++) cuml[n] /= var;
	  }

	  for (int n=1; n<=nmax; n++) {

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
	  for (int n=1; n<=nmax; n++) {
	    for (dd=0.0, nn=1; nn<=nmax; nn++) 
	      dd += Tevec[n][nn]*expcoef[indx][nn] * sqnorm[lm][nn];
	    smth[n] = dd * weight[indx][n];
	  }
	}

	inv = evec[indx]*smth;
	for (int n=1; n<=nmax; n++) {
	  expcoef[indx][n] = inv[n]/sqnorm[lm][n];
	  if (tk_type == 0) expcoef[indx][n] *= b_Hall[indx][n];
	}
	
	moffset += 2;
      }
    }
  }

#ifdef DEBUG
  // cerr << "pca_hall cpu_time=" << return_elapsed_time() << endl;
#endif

}

void AxisymmetricBasis::parallel_gather_coefficients(void)
{

  if (myid == 0) {

    for (int l=L0, loffset=0; l<=Lmax; loffset+=(2*l+1), l++) {

      for (int m=0, moffset=0; m<=l; m++) {

	if (m==0) {
	  for (int n=1; n<=nmax; ++n) {
	    expcoef[loffset+moffset][n] = 0.0;

	    for (int nn=n; nn<=nmax; nn++)
	      cc[loffset+moffset][n][nn] = 0.0;
	  }
	  moffset++;
	}
	else {
	  for (int n=1; n<=nmax; ++n) {
	    expcoef[loffset+moffset][n] = 0.0;
	    expcoef[loffset+moffset+1][n] = 0.0;

	    for (int nn=n; nn<=nmax; nn++) {
	      cc[loffset+moffset][n][nn] = 0.0;
	      cc[loffset+moffset+1][n][nn] = 0.0;
	    }
	  }
	  moffset+=2;
	}
      }
    }
  }


  for (int l=L0, loffset=0; l<=Lmax; loffset+=(2*l+1), l++) {

    for (int m=0, moffset=0; m<=l; m++) {

      if (m==0) {
	MPI_Reduce(&expcoef1[loffset+moffset][1], 
		   &expcoef[loffset+moffset][1], nmax, 
		   MPI_DOUBLE, MPI_SUM, 0, MPI_COMM_WORLD);

	for (int n=1; n<=nmax; n++)
	  MPI_Reduce(&cc1[loffset+moffset][n][n],
		     &cc[loffset+moffset][n][n], nmax-n+1, 
		     MPI_DOUBLE, MPI_SUM, 0, MPI_COMM_WORLD);
	moffset++;
      }
      else {
	MPI_Reduce(&expcoef1[loffset+moffset][1], 
		   &expcoef[loffset+moffset][1], nmax, 
		   MPI_DOUBLE, MPI_SUM, 0, MPI_COMM_WORLD);

	MPI_Reduce(&expcoef1[loffset+moffset+1][1],
		   &expcoef[loffset+moffset+1][1], nmax, 
		   MPI_DOUBLE, MPI_SUM, 0, MPI_COMM_WORLD);
	
	for (int n=1; n<=nmax; n++) {
	  MPI_Reduce(&cc1[loffset+moffset][n][n],
		     &cc[loffset+moffset][n][n], nmax-n+1, 
		     MPI_DOUBLE, MPI_SUM, 0, MPI_COMM_WORLD);
	  MPI_Reduce(&cc1[loffset+moffset+1][n][n],
		     &cc[loffset+moffset+1][n][n], nmax-n+1, 
		     MPI_DOUBLE, MPI_SUM, 0, MPI_COMM_WORLD);
	}

	moffset+=2;
      }
    }
  }

}

void AxisymmetricBasis::parallel_distribute_coefficients(void)
{

  for (int l=L0, loffset=0; l<=Lmax; loffset+=(2*l+1), l++) {

      for (int m=0, moffset=0; m<=l; m++) {

	if (m==0) {
	  MPI_Bcast(&expcoef[loffset+moffset][1], nmax, MPI_DOUBLE,
		    0, MPI_COMM_WORLD);
	  moffset++;
	}
	else {
	  MPI_Bcast(&expcoef[loffset+moffset][1], nmax, MPI_DOUBLE,
		     0, MPI_COMM_WORLD);
	  MPI_Bcast(&expcoef[loffset+moffset+1][1], nmax, MPI_DOUBLE,
		    0, MPI_COMM_WORLD);
	  moffset+=2;
	}
      }
  }

}


