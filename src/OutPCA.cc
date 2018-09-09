
#include <OutPCA.H>

#define GHQL

OutPCA::OutPCA(string& line) : Output(line)
{
  
}

void OutPCA::zeroth_order(void)
{
  static int done = 0;
  if (done) return;
  done = 1;

  string out1 = name + ".den0";
  string out2 = name + ".pot0";

  double r, dr, d, dd;

  ofstream cout1(out1.c_str());
  if (!cout1) {
    throw FileOpenError(out1, "OutPCA::zeroth_order", __FILE__, __LINE__);
  }
  
  ofstream cout2(out2.c_str());
  if (!cout2) {
    throw FileOpenError(out2, "OutPCA::zeroth_order", __FILE__, __LINE__);
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
    
    get_potl_dens(lmax, nmax, r, potd, dend);

    for (int n=1; n<=nmax; n++) {
      flag[n] = 1.0;

      get_dens_coefs(L_pca, flag.array(1, nmax), &d);

      cout1 << setw(14) << d;

      get_pot_coefs(L_pca, flag.array(1, nmax), &d, &dd);

      cout2 << setw(14) << d;
      flag[n] = 0.0;
    }
    cout1 << endl;
    cout2 << endl;
  }
}

void OutPCA::Run(int k)
{
  cout << "Entering out_pca . . .\n";
  if (!selector) return;	// Fail safe . . .

  ostringstream number;
  number << k;
  string curname = (string)outname + "." + number.str();

  cout << "About to call zeroth_order . . .\n";

  zeroth_order(curname);

  string out1 = curname + "_pca.den";
  string out2 = curname + "_pca.pot";

				/* Print out stuff */
  ofstream cout1(out1.c_str());
  if (!cout1) {
    throw FileOpenError(out1, "OutPCA::Run", __FILE__, __LINE__);
  }
  
  ofstream cout2(out2.c_str());
  if (!cout2) {
    throw FileOpenError(out2, "OutPCA::Run", __FILE__, __LINE__);
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
	throw GenericError("tk: no expansion defined", __FILE__, __LINE__);
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
	  throw GenericError("tk: no expansion defined", __FILE__, __LINE__);
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
	  throw GenericError("tk: no expansion defined", __FILE__, __LINE__);
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
