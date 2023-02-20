/*****************************************************************************
 *  Description:
 *  -----------
 *
 *  This routines in this module compute response matrix for given l, m, 
 *  and nmax.   Fully general version with complex frequencies.
 *
 *  Powerlaw measure for kappa integral
 *
 *  2-d (disk) version
 *
 *  Includes both on-axis and time-asymptotic evaluation of matrix elements
 *
 *  Call sequence:
 *  -------------
 *  void make_matrix_2d(CMatrix *mr, Complex omp, 
 *		 int l, int lmax, int nmax, int nptsK, int nptsE, 
 *		 AxiSymModel &model, AxiSymBiorth &biorth,
 *		 int OLD, int CPREAD, int CHECKPT, int VERBOSE, int RAT)
 *
 *
 *  Parameters:
 *  ----------
 *
 *  mr       response matrix
 *  omp      frequency
 *  l        harmonic order
 *  lmax     expansion integer l_1 will vary between [-lmax,lmax]
 *  nptsK    kappa integrand will contain 2^nptsK+1 points
 *  nptsE    energy integrand will contail 2^nptsE+1 points
 *  model    potential/density model
 *  biorth   biorthgonal functions
 *  OLD      addup all contributions l1, l2 before integrating
 *  CPREAD   read in old check.point file
 *  CHECKPT  write check.point files
 *  VERBOSE  write timing info
 *
 *  Returns:
 *  -------
 *
 *  None
 *
 *  Notes:
 *  -----
 *
 *  Pay attention to the following parameters:
 *
 *       DELTA     offset in exponent for kappa scaling
 *
 *       TOL1      fractional offset from tidal and central potential values
 *                 to prevent truncation errors
 *
 *       TOL       offset from 0 and 1 in kappa integrand to prevent trunca-
 *                 tion errors
 *
 *       RTOL      tolerance for simultaneous root cancellation in
 *                 rational function approximation (make small ok)
 *
 *
 *  By:
 *  --
 *
 *  MDW last update 17 jan 89
 *                  09 jan 93
 *                  15 jul 94 (adapatation for Orbit class)
 ***************************************************************************/

#define SVD_INVERSE 1

#include <cstdlib>
#include <cstring>
#include <string>

#include <global.H>
#include <RespMat.H>
#include <interp.H>
#include <Eigen/SVD>

// Global definitions

void rat_integral_set_parameters(double eps, int idbg=0);
std::complex<double> compute_rat_integral(double a, double b, Eigen::VectorXcd& x_data, Eigen::VectorXcd& y_data);

     
double get_exact_potl(int n, double r);	// Exact potential evaluation

static double svd_tol = 1.0e-6;

std::complex<double> Crombe(double a, double b, Eigen::VectorXcd& f);
double   Vrombe(double a, double b, Eigen::VectorXd& f);

void write_contr(ostream& out,
		 double xmin, double xmax, double ymin, double ymax, 
		 Eigen::MatrixXd& mat);

double sign(double x) 
{
  if (x<0) return -1.0;
  else return 1.0;
}


#define TOL1 1.0e-2
#define RTOL 1.0e-20
#define TOL 1.0e-2


RespMat::RespMat()
{
  ID = "UnknownDimensional(UnInitialized)";

				// Parameter flags

  dof = 0;
  isotropic = 0;
  matrix_computed = 0;
  matrix_computed_pv = 0;
  pv = 0;
  RATINT = 0;
  DELTA = 0.0;
  RECS = defaultRECS;
  PVALUE = defaultPVALUE;
  time_in_inner_loop=0.0;
  time_in_biorth_gen=0.0;
  CPREAD = 0;
  CHECKPT = 0;
  CHK_NAME = "";
  disp_computed = 0;
}

RespMat::RespMat(std::complex<double> omp0, int ll, int mm,
		 int llmax, int nnmax, int nnptsK, int nnptsE, 
		 AxiSymModPtr MODEL, AxiSymBioPtr BIORTH,
		 int old, int cpread, int checkpt, int verbose, 
		 IntegrationType type, 
		 int principal, int ratint, int ISOTROPIC)
{
  omp = omp0;
  l = ll;
  m = mm;
  lmax = llmax;
  nmax = nnmax;
  nptsK = nnptsK;
  nptsE = nnptsE;

  model = MODEL;
  biorth = BIORTH;

				// Parameter flags
  OLD     = old;
  CPREAD  = cpread;
  CHECKPT = checkpt;
  VERBOSE = verbose;
  SITYPE  = type;

  dof = model->dof();

  switch (dof) {
  case 2:
    ID = "Two";
    CHK_NAME = "check.point2d";
    l = m;
    break;
  case 3:
    ID = "Three";
    CHK_NAME = "check.point3d";
    break;
  default:
    bomb("Initialization problem: only 2 or 3 dimensions!");
  }

  ID += "Dimensional(" + model->ModelID + ", " + biorth->BiorthID + ")";

  isotropic = ISOTROPIC;
  if (principal) pv = 1;
  else pv = 0;
				// Initialization

  matrix_computed = 0;
  matrix_computed_pv = 0;
  RATINT = ratint;
  DELTA = 0.0;
  RECS = defaultRECS;
  PVALUE = defaultPVALUE;
  time_in_inner_loop=0.0;
  time_in_biorth_gen=0.0;
  CPREAD = 0;
  CHECKPT = 0;
  CHK_NAME = "";
  disp_computed = 0;
}


RespMat::RespMat(const RespMat& p)
{
  mr = p.mr;
  mr2 = p.mr2;
  omp = p.omp;
  l = p.l;
  m = p.m;
  lmax = p.lmax;
  nmax = p.nmax;
  nptsK = p.nptsK;
  nptsE = p.nptsE;
  model = p.model;
  biorth = p.biorth;
  OLD = p.OLD;
  CPREAD = p.CPREAD;
  CHECKPT = p.CHECKPT;
  VERBOSE = p.VERBOSE;
  SITYPE = p.SITYPE;
  dof = p.dof;
  isotropic = p.isotropic;
  pv = p.pv;

  matrix_computed = p.matrix_computed;
  matrix_computed_pv = p.matrix_computed_pv;
  ID = p.ID;

  disp_computed = p.disp_computed;
  disp = p.disp;
}

/*
RespMat &RespMat::operator=(RespMat& p)
{
  mr = p.mr;
  mr2 = p.mr2;
  omp = p.omp;
  l = p.l;
  m = p.m;
  lmax = p.lmax;
  nmax = p.nmax;
  nptsK = p.nptsK;
  nptsE = p.nptsE;
  model = p.model;
  biorth = p.biorth;
  OLD = p.OLD;
  CPREAD = p.CPREAD;
  CHECKPT = p.CHECKPT;
  VERBOSE = p.VERBOSE;
  SITYPE = p.SITYPE;
  dof = p.dof;
  pv = p.pv;

  matrix_computed = p.matrix_computed;
  matrix_computed_pv = p.matrix_computed_pv;
  ID = p.ID;

  hold = NULL;

  disp_computed = p.disp_computed;
  disp = p.disp;

  return *this;
}
*/

void RespMat::set_principal_value(int pvalue)
{
  switch(pvalue) {
  case 0:
    PVALUE = 0;
    break;
  case 1:
    PVALUE = 1;
    break;
  case 2:
    PVALUE = 2;
    break;
  default:
    bomb("pvalue must be 0, 1, or 2");
  }
}


void RespMat::make_matrix(void)
{
  double dk, kap2, expok, fac, tol1, tol2, tot_time, tot_time_d;
  std::complex<double> temp;
  
  /* Set up global values */
  
  if (dof!=2 && dof!=3) {	// Sanity check
    cerr << "RespMat::make_matrix: Wrong number of degrees of freedom (help!)" << endl;
    exit(-1);
  }

  if (dof==2) l2 = l;
  num_E = 1 + (int)(pow(2.0,(double)nptsE)+1.0e-6);
  num_K = 1 + (int)(pow(2.0,(double)nptsK)+1.0e-6);
  pp_x.resize(num_E);
  pp_y.resize(num_E);
  mr.  resize(nmax, nmax);
  mr2. resize(nmax, nmax);
  
  expok = 1.0/(dof-1.0-DELTA);
  tol1 = pow(TOL, 1.0/expok);
  tol2 = TOL/expok;
  dk = (1.0 - tol1-tol2)/(double)(num_K-1);

  Eigen::VectorXcd mk(num_K);
  Eigen::VectorXcd wk(num_K);
  
  switch (SITYPE) {
  case ratint:
    rat_integral_set_parameters(RTOL);
    break;
  case jacoint:
    {
      JacoQuad jq(num_K, dof-2.0-DELTA, 0.0);
      for (int k=0; k<num_K; k++) {
	mk[k-1] = jq.knot(k);
	wk[k-1] = jq.weight(k);
      }
    }
    break;
  case rombint:
    for (int k=0; k<num_K; k++) {
      for (int k=0; k<num_K; k++) {
	mk[k] = TOL + dk*k;
	wk[k] = 1.0;
      }
    }
    break;
  }
  
  /* Set up for matrix element computation
     define and clean array . . . */
  
  mab.resize(nmax*nmax);
  for (auto & v : mab) {
    v.resize(num_K);
    v.setZero();
  }
  
  /* Set up for array of orbits in <E> at fixed <kappa> */
  
  Emodmin = model->get_pot(model->get_min_radius());
  Emodmax = model->get_pot(model->get_max_radius());
  
  /* Initialize orbit list */
  
  Orb = std::make_shared<OrbitTable>(getptr(), RECS);
  
  /* Set up grids to do angular momentum <kappa> integrals */
  
  ikap = 1;
  /* Read in old check point if it exists */
  if (CPREAD) get_chkpnt();
  for (; ikap<=num_K; ikap++) {
    if (SITYPE==jacoint)
      kappa = mk[ikap].real();
    else {
      kap2 = tol1 + dk*(ikap-1);
      kappa = pow(kap2,expok);
    }
    get_emax();
    
    timer.start();
    
    if (OLD)
      integrand_rolled();	// Original summed order
    else
      integrand_unroll();	// Loop swapped to better resolve poles
    
    if (VERBOSE) {

      tot_time = timer.getTime();
      cout << "\nTime for step #" << ikap << ":            " 
	   << tot_time << endl;
      
      tot_time_d = tot_time;
      if (tot_time_d <= 0.0) tot_time_d = 1.0;

      cout << "    Time in integration    : " 
	   << time_in_inner_loop << " (" << setprecision(3)
	   << time_in_inner_loop/tot_time_d*100 
	   << "%)" << endl;

      cout << "    Time in biorthog lookup: " 
	   << time_in_biorth_gen << " (" << setprecision(3)
	   << time_in_biorth_gen/tot_time_d*100 
	   << "%)" << endl;
      
      cout << "    Everything else        : " 
	   << tot_time - time_in_inner_loop - time_in_biorth_gen
	   << " (" << setprecision(3)
	   << (1.0 - time_in_biorth_gen/tot_time_d - 
	       time_in_inner_loop/tot_time_d)*100 
	   << "%)" << endl;
      timer.reset();
      time_in_inner_loop = 0.0;
      time_in_biorth_gen = 0.0;
    }
    
    if (CHECKPT) put_chkpnt();
  }
  
  if (SITYPE != jacoint) {

    for (int k=0; k<num_K; k++) {
      fac = pow(mk[k].real(), DELTA) * expok;
      
      for (alpha=0; alpha<nmax; alpha++) {
	for (beta=alpha; beta<nmax; beta++) {
	  mab[alpha*nmax+beta][k] *= fac;
	}
      }
    }

  }
  
  for (alpha=0; alpha<=nmax; alpha++) {
    for (beta=alpha; beta<=nmax; beta++) {
      
      switch (SITYPE) {
      case ratint:
	temp = compute_rat_integral(0.0, 1.0, mk, mab[alpha*nmax+beta]);
	break;
      case jacoint:
	temp = wk.adjoint().dot(mab[alpha*nmax+beta]);
	break;
      case rombint:
	temp = Crombe(mk[1].real(), mk[num_K].real(), mab[alpha*nmax+beta]);
	break;
      }
      
      if (dof==2) {
//
// Explanation of multiplicative factors:
// -------------------------------------
//
// *****Two-dimensional
//
//----------------------------------------------------------------------
//
//                /------To satify Possion's equation (G=1)
//                |
//                v
	temp *= - 2.0*M_PI * 4.0*M_PI*M_PI / (2.0*M_PI) * 2.0;
//                           ^                ^           ^
//                           |                |           |
//  From angle integration---/                |           |
//                                            |           |
//   Angular normalization for biortho. fcts--/           |
//                                                        |
//   Convert DF from v_r, v_t to E, L---------------------/
//
//----------------------------------------------------------------------
      } else {
// 
// *****Three-dimensional
//
//----------------------------------------------------------------------
//
//                /------To satify Possion's equation
//                |
//                v
	temp *= - 4.0*M_PI * 8.0*M_PI*M_PI*M_PI * 2.0/(2.0*l+1.0);
//                           ^                    ^
//                           |                    |
//                           |                    |
//                      From rotational matrices--/
//                           |
// From angle integration----/
//
//----------------------------------------------------------------------
//
      }
      mr(alpha, beta) = temp;
      mr(beta, alpha) = temp;
    }
  }

				// Set computed flag
  matrix_computed = 1;
  
}

void RespMat::make_matrix_pv(void)
{
  double dk, kap2, expok, fac, tol1, tol2, tot_time, tot_time_d;
  std::complex<double> temp, temp2;
  
  // Set up global values
  //
  if (dof!=2 && dof!=3) {	// Sanity check
    cerr << "RespMat::make_matrix_pv: Wrong number of degrees of freedom (help!)" << endl;
    exit(-1);
  }

  if (dof==2) l2 = l;
  num_E = 1 + (int)(pow(2.0,(double)nptsE)+1.0e-6);
  num_K = 1 + (int)(pow(2.0,(double)nptsK)+1.0e-6);
  pp_x. resize(num_E);
  pp_y. resize(num_E);
  pp_x2.resize(num_E);

  mr. resize(nmax, nmax);
  mr2.resize(nmax, nmax);
  
  expok = 1.0/(dof-1.0-DELTA);
  tol1 = pow(TOL, 1.0/expok);
  tol2 = TOL/expok;
  dk = (1.0 - tol1-tol2)/(double)(num_K-1);

  Eigen::VectorXcd mk(num_K);
  Eigen::VectorXcd wk(num_K);
  
  switch (SITYPE) {
  case ratint:
    rat_integral_set_parameters(RTOL);
    break;
  case jacoint:
    {
      JacoQuad jq(num_K, dof-2.0-DELTA, 0.0);
      for (int k=0; k<num_K; k++) {
	mk[k] = jq.knot(k);
	wk[k] = jq.weight(k);
      }
    }
    break;
  case rombint:
    for (int k=0; k<num_K; k++) {
      mk[k] = TOL + dk*k;
      wk[k] = 1.0;
    }
    break;
  }
  
  // Set up for matrix element computation
  // define and clean array . . .
  //
  mab.resize(nmax*nmax);
  for (auto & v : mab) {
    v.resize(num_K);
    v.setZero();
  }
  
  // Set up for array of orbits in <E> at fixed <kappa>
  //
  Emodmin = model->get_pot(model->get_min_radius());
  Emodmax = model->get_pot(model->get_max_radius());
  
  // Initialize orbit list
  //
  ORes = make_shared<OrbResTable>(getptr(), RECS);
  
  // Set up grids to do angular momentum <kappa> integrals
  //
  ikap = 1;
  // Read in old check point if it exists
  if (CPREAD) get_chkpnt();
  for (; ikap<=num_K; ikap++) {
    if (SITYPE==jacoint)
      kappa = mk[ikap].real();
    else {
      kap2 = tol1 + dk*(ikap-1);
      kappa = pow(kap2,expok);
    }
    get_emax();
    
    timer.start();
    
    integrand_pv();

    
    if (VERBOSE) {

      tot_time = timer.getTime();
      cout << "\nTime for step #" << ikap << ":            " 
	   << tot_time << '\n';

      tot_time_d = tot_time;
      if (tot_time_d <= 0.0) tot_time_d = 1.0;

      cout << "    Time in integration    : " 
	<< time_in_inner_loop << " (" << setprecision(3)
	  << time_in_inner_loop/tot_time_d*100 
	    << "%)\n";

      cout << "    Time in biorthog lookup: " 
	<< time_in_biorth_gen << " (" << setprecision(3)
	  << time_in_biorth_gen/tot_time_d*100 
	    << "%)\n";
      
      cout << "    Everything else        : " 
	<< tot_time - time_in_inner_loop - time_in_biorth_gen
	  << " (" << setprecision(3)
	    << (1.0 - time_in_biorth_gen/tot_time_d - 
		time_in_inner_loop/tot_time_d)*100 
		  << "%)\n";
      timer.reset();
      time_in_inner_loop = 0.0;
      time_in_biorth_gen = 0.0;
    }
    
    if (CHECKPT) put_chkpnt();
  }
  
  if (SITYPE != jacoint) {

    for (int k=0; k<num_K; k++) {
      fac = pow(mk[k].real(), DELTA) * expok;
      
      for (alpha=0; alpha<nmax; alpha++) {
	for (beta=alpha; beta<nmax; beta++) {
	  mab[alpha*nmax+beta][k] *= fac;
	}
      }
    }

  }
  
  for (alpha=0; alpha<nmax; alpha++) {
    for (beta=alpha; beta<nmax; beta++) {
      
      switch (SITYPE) {
      case ratint:
	temp = compute_rat_integral(0.0, 1.0, mk, mab[alpha*nmax+beta]);
	break;
      case jacoint:
	temp = wk.adjoint().dot(mab[alpha*nmax+beta]);
	break;
      case rombint:
	temp = Crombe(mk[1].real(), mk[num_K].real(), mab[alpha*nmax+beta]);
	break;
      }
      
      if (dof==2) {
//
// Explanation of multiplicative factors:
// -------------------------------------
//
// *****Two-dimensional
//
//----------------------------------------------------------------------
//
//                /------To satify Possion's equation (G=1)
//                |
//                v
	temp *= - 2.0*M_PI * 4.0*M_PI*M_PI / (2.0*M_PI) * 2.0;
//                           ^                ^           ^
//                           |                |           |
//  From angle integration---/                |           |
//                                            |           |
//   Angular normalization for biortho. fcts--/           |
//                                                        |
//   Convert DF from v_r, v_t to E, L---------------------/
//
//----------------------------------------------------------------------
      } else {
// 
// *****Three-dimensional
//
//----------------------------------------------------------------------
//
//                /------To satify Possion's equation
//                |
//                v
	temp *= - 4.0*M_PI * 8.0*M_PI*M_PI*M_PI * 2.0/(2.0*l+1.0);
//                           ^                    ^
//                           |                    |
//                           |                    |
//                      From rotational matrices--/
//                           |
// From angle integration----/
//
//----------------------------------------------------------------------
//
      }
      mr(alpha, beta) = temp;
      mr(beta, alpha) = temp;
    }
  }

  for (alpha=0; alpha<nmax; alpha++) {
    for (beta=0; beta<nmax; beta++) {
      
      mr2(beta, alpha) = std::complex<double>(
					      mr(alpha, beta).real(),
					      2.0*mr(alpha, beta).imag()
					      );
    }
  }

  // Set computed flag
  //
  matrix_computed = 1;
  matrix_computed_pv = 1;
}

#undef TOL


/****************************************************************************
 *
 *  Integrand routines for energy <E> integrals
 *
 ****************************************************************************/

// Integrand for real (off-resonant) elements (kappa ang.mom. integral)
//
OrbitTable::OrbitTable(std::shared_ptr<RespMat> pp, int recs)
{
  p = pp;			// Pointer to calling object

  orbits.resize(p->num_E);
  for (int i=0; i<p->num_E; i++) {
    orbits[i] = SphericalOrbit(p->model);
    orbits[i].set_numerical_params(recs);
    orbits[i].set_biorth(p->biorth, p->l, p->nmax, 1);
  }

  PotTrans.resize(p->num_E);
  EInt.resize(p->num_E);
  dfqE.resize(p->num_E);
  dfqL.resize(p->num_E);
  norm.resize(p->nmax);

  for (auto & v : PotTrans) v.resize(p->nmax);
  for (int i=0; i<p->nmax; i++)
    norm[i] = 1.0/sqrt((p->biorth)->norm(i, p->l));
}

OrbitTable::~OrbitTable()
{
  // None
}

void OrbitTable::setup()
{
  double dE, fac;

  dE = (p->Emax - p->Emodmin - 2.0*p->ESTEP)/(double)(p->num_E-1);
  for (int i=0; i<p->num_E; i++) {
    p->pp_x[i] = p->Emodmin + p->ESTEP + dE*(i-1);
    orbits[i].new_orbit(p->pp_x[i].real(), p->kappa);
    if (p->dof==2)
      fac = orbits[i].Jmax()/orbits[i].get_freq(1);
    else
      fac = orbits[i].Jmax()*orbits[i].Jmax()/orbits[i].get_freq(1);
    dfqE[i] = fac * 
      (p->model)->dfde(orbits[i].Energy(), orbits[i].AngMom());
    dfqL[i] = fac * 
      (p->model)->dfdl(orbits[i].Energy(), orbits[i].AngMom());
  }
}

// For each E (and K), compute run of ortho functions for each component
//
void OrbitTable::setup_p_array(void)
{
  double fac, facd;

  if (p->dof==3) 
    facd = p->Ylm02(p->l, p->l2);
  else
    facd = 1.0;

  for (int i=0; i<p->num_E; i++) {
    fac = orbits[i].get_freq(1)*p->l1 + orbits[i].get_freq(2)*p->l2;
    EInt[i] = (dfqE[i]*fac + dfqL[i]*p->l2) * facd/( p->omp + fac );

  }

  double t0 = 0.0;
  if (p->VERBOSE) {
    t0 = timer.getTime();
  }

  for (int i=0; i<p->num_E; i++)
    orbits[i].pot_trans(p->l1, p->l2, PotTrans[i]);
  if (p->VERBOSE) {
    p->time_in_biorth_gen += timer.getTime() - t0;
  }
  
}



/****************************************************************************
 *
 *  Integrand routines for locating resonances in <E>
 *
 ****************************************************************************/

/* Integrand for real (off-resonant) elements (kappa ang.mom. integral) */

				// Define a default
double OrbResTable::derivative_step = 1.0e-3;

OrbResTable::OrbResTable(std::shared_ptr<RespMat> pp, int recs)
{

  p = pp;			// Pointer to calling object

  orbits.resize(p->num_E);

  for (int i=0; i<p->num_E; i++) {
    orbits[i] = std::make_shared<SphericalOrbit>(p->model);
    orbits[i]->set_numerical_params(recs);
    orbits[i]->set_biorth(p->biorth, p->l, p->nmax, 1);
  }

  torbit = SphericalOrbit(p->model);
  torbit.set_numerical_params(recs);
  torbit.set_biorth(p->biorth, p->l, p->nmax, 1);
    
  PotTrans .resize(p->num_E);
  PotTrans2.resize(p->num_E);

  for (auto & v : PotTrans)  v.resize(p->nmax);
  for (auto & v : PotTrans2) v.resize(p->nmax);

  EInt.resize(p->num_E);
  ELoc.resize(p->num_E);
  ERes.resize(p->num_E);
  EJac.resize(p->num_E);
  dfqE.resize(p->num_E);
  dfqL.resize(p->num_E);
  norm.resize(p->nmax);
  
  for (int i=0; i<p->nmax; i++)
    norm[i] = 1.0/sqrt((p->biorth)->norm(i, p->l));
}

OrbResTable::~OrbResTable()
{
}

void OrbResTable::setup()
{
  double dE, fac;
  int i;

  dE = (p->Emax - p->Emodmin - 2.0*p->ESTEP)/(double)(p->num_E-1);
  for (i=1; i<=p->num_E; i++) {
    p->pp_x[i] = p->Emodmin + p->ESTEP + dE*(i-1);
    p->pp_x2[i] = p->Emodmin + p->ESTEP + dE*(i-1);
    orbits[i]->new_orbit(p->pp_x2[i], p->kappa);
    if (p->dof==2)
      fac = orbits[i]->Jmax()/orbits[i]->get_freq(1);
    else
      fac = orbits[i]->Jmax()*orbits[i]->Jmax()/orbits[i]->get_freq(1);
    dfqE[i] = fac *
      (p->model)->dfde(orbits[i]->Energy(), orbits[i]->AngMom());
    dfqL[i] = fac *
      (p->model)->dfdl(orbits[i]->Energy(), orbits[i]->AngMom());
  }
}

// For each E (and K), compute run of ortho functions for each component
//
void OrbResTable::setup_E_array(void)
{
  double Et, last, fac=0, facd, facdf, t0=0, dfdE, dfdL;
  double fac_p, fac_m;

  if (p->dof==3) 
    facd = p->Ylm02(p->l, p->l2);
  else
    facd = 1.0;

  number = 0;

  for (int i=0; i<p->num_E; i++) {
    last = fac;
    fac = p->omp.real() + 
      orbits[i]->get_freq(1)*p->l1 + orbits[i]->get_freq(2)*p->l2;
    EInt[i] = (dfqE[i]*(fac - p->omp.real()) + dfqL[i]*p->l2) * facd/fac;
      
    if (i==0) continue;

    if (last*fac <= 0.0) {
				// Interpolate to derive E res

      number++;

      Et = ELoc[number] = (fac*p->pp_x2[i-1] - last*p->pp_x2[i])/(fac - last);
      torbit.new_orbit(Et, p->kappa);
      
      if (p->dof==2)
	facdf = torbit.Jmax()/torbit.get_freq(1);
      else
	facdf = torbit.Jmax()*torbit.Jmax()/torbit.get_freq(1);

      dfdE = facdf * 
	(p->model)->dfde(torbit.Energy(), torbit.AngMom());

      dfdL = facdf * 
	(p->model)->dfdl(torbit.Energy(), torbit.AngMom());
      
      facdf = torbit.get_freq(1)*p->l1 + torbit.get_freq(2)*p->l2;

      ERes[number] = (dfdE*facdf + dfdL*p->l2) * facd;

      if (p->VERBOSE) {
	t0 = timer.getTime();
      }
      torbit.pot_trans(p->l1, p->l2, PotTrans2[number]);
      if (p->VERBOSE) {
	p->time_in_biorth_gen += timer.getTime() - t0;
      }
				// Jacobian of delta function

      torbit.new_orbit(Et*(1.0 - derivative_step), p->kappa);
      fac_p = torbit.get_freq(1)*p->l1 + torbit.get_freq(2)*p->l2;
      torbit.new_orbit(Et*(1.0 + derivative_step), p->kappa);
      fac_m = torbit.get_freq(1)*p->l1 + torbit.get_freq(2)*p->l2;

      EJac[number]  = (fac_p - fac_m)/(-2.0*Et*derivative_step);
      ERes[number] /= EJac[number];

    }

  }
  
  if (p->VERBOSE) {
    t0 = timer.getTime();
  }
  for (int i=1; i<=p->num_E; i++)
    orbits[i]->pot_trans(p->l1, p->l2, PotTrans[i]);
  if (p->VERBOSE) {
    p->time_in_biorth_gen += timer.getTime()  - t0;
  }
  
}


void RespMat::integrand_unroll(void)
{	   
  std::complex<double> temp;

  Orb->setup();

  int llm, llp;
  if (dof==3) {
    llm = -l;
    llp =  l;
  } else
    llm = llp = m;

				// Fourier loop

  for (l1=(-lmax); l1<=lmax;  l1++) {

    for (l2=llm; l2<=llp; l2+=2) {

      if ( (l1==0) && (l2==0)) continue;

      Orb->setup_p_array();

				// ALPHA & BETA loop
    
      for (alpha=0; alpha<nmax; alpha++) {

	for (beta=alpha; beta<nmax; beta++) {

	  
				//* Do <E> integral
	  
	  for (int k=0; k<num_E; k++)
	    pp_y[k] = Orb->EInt[k] * Orb->PotTrans[k][alpha] *
	      Orb->PotTrans[k][beta];

	  double t0=0;
	  if (VERBOSE)
	    t0 = timer.getTime();
	  if (RATINT)
	    temp = compute_rat_integral(Emodmin, Emax, pp_x, pp_y);
	  else
	    temp = Crombe(pp_x[1].real(), pp_x[num_E].real(), pp_y);

	  if (VERBOSE) {
	    time_in_inner_loop += timer.getTime() - t0;
	  }
	  mab[alpha*nmax+beta][ikap] += temp*Orb->norm[alpha]*Orb->norm[beta];
	}
      }
    }
  }
  
}

void RespMat::integrand_rolled(void)
{	   
  std::complex<double> temp;

  if (hold.size()==0) {
    hold.resize(nmax*nmax);
    for (auto & v : hold) v.resize(num_E);
  }

				// Clear storage
  for (auto & v : hold) v.setZero();

  int llm, llp;
  if (dof==3) {
    llm = -l;
    llp =  l;
  } else
    llm = llp = m;


  Orb->setup();

				// Fourier loop

  for (l1=(-lmax); l1<=lmax;  l1++) {

    for (l2=llm; l2<=llp; l2+=2) {

      if ( (l1==0) && (l2==0)) continue;

      Orb->setup_p_array();

				// ALPHA & BETA loop
      
      for (alpha=0; alpha<=nmax; alpha++) {

	for (beta=alpha; beta<=nmax; beta++) {

	  for (int k=0; k<num_E; k++)
	    hold[alpha*nmax+beta][k] += Orb->EInt[k] * 
	      Orb->PotTrans[k][alpha] * Orb->PotTrans[k][beta];
	  
	}
      }
    }
  }

  for (alpha=0; alpha<nmax; alpha++) {

    for (beta=alpha; beta<nmax; beta++) {

      double t0=0;
      if (VERBOSE) {
	t0 = timer.getTime();
      }
      if (RATINT)
	temp = compute_rat_integral(Emodmin, Emax, pp_x, hold[alpha*nmax+beta]);
      else
	temp = Crombe(pp_x[1].real(), pp_x[num_E].real(), hold[alpha*nmax+beta]);
      if (VERBOSE) {
	time_in_inner_loop += timer.getTime() - t0;
      }
      mab[alpha*nmax+beta][ikap] += temp*Orb->norm[alpha]*Orb->norm[beta];

    }
  }

}

void RespMat::integrand_pv(void)
{	   
  std::complex<double> temp, fac;
  int nres=0;
  double t0=0;

  if (hold.size()==0) {
    hold.resize(nmax*nmax);
    for (auto & v : hold) v.resize(num_E);
    holdR.resize(nmax, nmax);
  }

				// Clear storage
  for (auto & v : hold) v.setZero();
  holdR.setZero();


  int llm, llp;
  if (dof==3) {
    llm = -l;
    llp =  l;
  } else
    llm = llp = m;


  ORes->setup();

				// Fourier loop

  for (l1=(-lmax); l1<=lmax;  l1++) {

    for (l2=llm; l2<=llp; l2+=2) {

      if ( (l1==0) && (l2==0)) continue;

      ORes->setup_E_array();

      int num = ORes->number;

      if (num == 0) {

	nres++;

	for (alpha=0; alpha<nmax; alpha++) {
	  for (beta=alpha; beta<nmax; beta++) {

	    for (int k=0; k<num_E; k++) {

	      hold[alpha*nmax+beta][k] += ORes->EInt[k] * 
		ORes->PotTrans[k][alpha] * ORes->PotTrans[k][beta];
	  
	    }
	  }
	}

      }
      else if (num > 1) {

	if (VERBOSE) {
	  t0 = timer.getTime();
	}

	int n = 1;
	double lastE;

	// Principal value with singularity subtracted
	// Trapezoidal rule
	

	// For residues and singularities

	
	lastE = pp_x2[0];

	for (int k=0; k<num_E; k++) {
	      
	  for (alpha=0; alpha<nmax; alpha++) {

	    for (beta=alpha; beta<nmax; beta++) {

	      holdR(alpha, beta) += 0.5 * (pp_x2[k] - pp_x2[k-1]) *
		(
		 ORes->EInt[k]*
		 ORes->PotTrans[k][alpha]*ORes->PotTrans[k][beta] -
		 ORes->ERes[n]*
		 ORes->PotTrans2[n][alpha]*ORes->PotTrans2[n][beta]/
		 (pp_x2[k] - ORes->ELoc[n])
		 +
		 ORes->EInt[k-1]*
		 ORes->PotTrans[k-1][alpha]*ORes->PotTrans[k-1][beta] -
		 ORes->ERes[n]*
		 ORes->PotTrans2[n][alpha]*ORes->PotTrans2[n][beta]/
		 (pp_x2[k-1] - ORes->ELoc[n])
		 );
	    }
	  }

	  if (n<num && 
	      (pp_x2[k] > 0.5*(ORes->ELoc[n] + ORes->ELoc[n+1]) ||
	       pp_x2[k+1] > ORes->ELoc[n+1]) ) {
	    
	    // Add singularity and residue

	    // DEBUG

	    int tstflg=0;
	    double tst = (pp_x2[k]-ORes->ELoc[n])/(ORes->ELoc[n]-lastE);
	    if (tst < 0.0) {
	      cerr << "check domain: less than zero\n";
	      tstflg++;
	    }
	    if (std::isnan(tst)) {
	      cerr << "check domain: NaN\n";
	      tstflg++;
	    }
	    if (std::isinf(tst)) {
	      cerr << "check domain: Infinity\n";
	      tstflg++;
	    }
	    if (tst==0.0) {
	      cerr << "check domain: Zero\n";
	      tstflg++;
	    }
	    if (tstflg) {
	      cerr << tstflg << " errors\n";
	    }
	

	    // END DEBUG


	    fac =  ORes->ERes[n] *
	      (  
	       log((pp_x2[k]-ORes->ELoc[n])/(ORes->ELoc[n]-lastE)) 
	       + M_PI*std::complex<double>(0.0, 1.0)*sign(ORes->EJac[n])
	       );
	    
	    for (alpha=0; alpha<=nmax; alpha++) {
	      for (beta=alpha; beta<=nmax; beta++)
		holdR(alpha, beta) += fac *
		  ORes->PotTrans2[n][alpha]*ORes->PotTrans2[n][beta];
	    }

	    lastE = pp_x2[k];

	    // Advance resonance counter
	    n++;

	  }
	  
	}

	// Add *final* singularity and residue

	// DEBUG

	int tstflg=0;
	double tst = (pp_x2[num_E]-ORes->ELoc[n])/(ORes->ELoc[n]-lastE);
	if (tst < 0.0) {
	  cerr << "check domain: less than zero\n";
	  tstflg++;
	}
	if (std::isnan(tst)) {
	  cerr << "check domain: NaN\n";
	  tstflg++;
	}
	if (std::isinf(tst)) {
	  cerr << "check domain: Infinity\n";
	  tstflg++;
	}
	if (tst==0.0) {
	  cerr << "check domain: Zero\n";
	  tstflg++;
	}
	if (tstflg) {
	  cerr << tstflg << " errors\n";
	}
	

	// END DEBUG


	fac =  ORes->ERes[n] *
	  (  
	   log((pp_x2[num_E]-ORes->ELoc[n])/(ORes->ELoc[n]-lastE)) 
	   + M_PI*std::complex<double>(0.0, 1.0)*sign(ORes->EJac[n])
	   );
	    
	for (alpha=1; alpha<=nmax; alpha++) {
	  for (beta=alpha; beta<=nmax; beta++)
	    holdR(alpha, beta) += fac *
	      ORes->PotTrans2[n][alpha]*ORes->PotTrans2[n][beta];
	}
	
	if (VERBOSE) {
	  time_in_inner_loop += timer.getTime() - t0;
	}

      }
      else {

	if (VERBOSE) {
	  t0 = timer.getTime();
	}

	// Principal value with singularity subtracted

	for (int k=0; k<num_E; k++) {

	  for (alpha=0; alpha<nmax; alpha++) {
	    for (beta=alpha; beta<nmax; beta++)
	    
	      hold[alpha*nmax+beta][k] += ORes->EInt[k] *
		ORes->PotTrans[k][alpha] * ORes->PotTrans[k][beta] -
		ORes->ERes[num] *
		ORes->PotTrans2[num][alpha] * ORes->PotTrans2[num][beta]/
		(pp_x2[k] - ORes->ELoc[num]);
	  }

	}

	// Add residue and singularity

	// DEBUG

	int tstflg=0;
	double tst = (pp_x2[num_E]-ORes->ELoc[num])/(ORes->ELoc[num]-pp_x2[1]);
	if (tst < 0.0) {
	  cerr << "check domain: less than zero" << endl;
	  tstflg++;
	}
	if (std::isnan(tst)) {
	  cerr << "check domain: NaN" << endl;
	  tstflg++;
	}
	if (std::isinf(tst)) {
	  cerr << "check domain: Infinity" << endl;
	  tstflg++;
	}
	if (tst==0.0) {
	  cerr << "check domain: Zero\n";
	  tstflg++;
	}
	if (tstflg) {
	  cerr << tstflg << " errors\n";
	}
	

	// END DEBUG

	fac =  ORes->ERes[num] *
	  (  
	   log((pp_x2[num_E]-ORes->ELoc[num])/(ORes->ELoc[num]-pp_x2[1]))
	   + M_PI*std::complex<double>(0.0, 1.0)*static_cast<double>(PVALUE)*sign(ORes->EJac[num])
	   );
	
	for (alpha=0; alpha<nmax; alpha++) {
	  for (beta=alpha; beta<nmax; beta++)
	    holdR(alpha, beta) += fac *
	      ORes->PotTrans2[num][alpha] * ORes->PotTrans2[num][beta];
	}

	if (VERBOSE) {
	  time_in_inner_loop += timer.getTime() - t0;
	}

      }


    }
  }

  for (alpha=0; alpha<nmax; alpha++) {
    for (beta=alpha; beta<nmax; beta++) {

      if (VERBOSE) {
	t0 = timer.getTime();
      }

      if (RATINT)
	temp = compute_rat_integral(Emodmin, Emax, pp_x, hold[alpha*nmax+beta]);
      else
	temp = Crombe(pp_x[0].real(), pp_x[num_E-1].real(), hold[alpha*nmax+beta]);
      

      if (VERBOSE) {
	time_in_inner_loop += timer.getTime() - t0;
      }

      mab[alpha*nmax+beta][ikap] += (temp + holdR(alpha, beta)) *
	ORes->norm[alpha] * ORes->norm[beta];

    }
  }

}


				// For comparsion with old version

void RespMat::integrand_orig(void)
{	   
  std::complex<double> temp;
  double fac, facd;
  int i;
  int llm, llp;


  if (dof==3) {
    llm = -l;
    llp =  l;
  } else {
    llm = llp = m;
    facd = 1.0;
  }

  Orb->setup();

				// ALPHA & BETA loop
  
  for (alpha=0; alpha<nmax; alpha++) {
    
    for (beta=alpha; beta<nmax; beta++) {

      pp_y.setZero();
				// Fourier loop

      for (l1=(-lmax); l1<=lmax;  l1++) {

	for (l2=llm; l2<=llp; l2+=2) {
	
	  if (dof==3) facd = Ylm02(l, l2);

	  if ( (l1==0) && (l2==0)) continue;
	  
	  for (i=1; i<=num_E; i++) {
	    fac = Orb->orbits[i].get_freq(1)*l1 + Orb->orbits[i].get_freq(2)*l2;
	    pp_y[i] += 
	      (Orb->dfqE[i]*fac + Orb->dfqL[i]*l2) * facd/( omp + fac ) *
		Orb->orbits[i].pot_trans(l1, l2, alpha) *
		  Orb->orbits[i].pot_trans(l1, l2, beta);
	  }
	}
      }
    
				// Do <E> integral

      double t0 = 0;
      if (VERBOSE) {
	t0 = timer.getTime();
      }
      if (RATINT)
	temp = compute_rat_integral(Emodmin, Emax, pp_x, pp_y);
      else
	temp = Crombe(pp_x[1].real(), pp_x[num_E].real(), pp_y);
      // temp = compute_rat_integral(Emodmin, Emax, pp_x, pp_y);
      time_in_inner_loop += timer.getTime() - t0;
      mab[alpha*nmax+beta][ikap] += temp*Orb->norm[alpha]*Orb->norm[beta];
    }
  }
  
}



/****************************************************************************
 *
 *     Support and math routines
 *
 ***************************************************************************/


/* angular factors . . . (Y_lm)^2*/

double RespMat::Ylm02(int ll, int mm)
{
  mm = abs(mm);
  return (2.0*ll+1)/(4.0*M_PI*M_PI) * exp(
         2.0*mm*log(2.0) + lgamma(1.0+ll-mm) - lgamma(1.0+ll+mm) +
	 2.0*lgamma(0.5*(ll+mm+1)) - 2.0*lgamma(0.5*(ll-mm)+1.0) );
}


/*  These two routines put and get check point data sets to protect
    our hard earned CPU time */

bool RespMat::get_chkpnt(void)
{
  std::complex<double> omp01;
  double rtmp, itmp;
  int ll1,mm1,llmax1,nnmax1,nptsK1,nptsE1;

  std::ifstream fin((outdir + CHK_NAME).c_str());
  if (!fin) {
    std::cerr << "Couldn't open " << outdir + CHK_NAME << std::endl;
    std::cerr << "continuing . . ." << std::endl;
    return false;
  }

  
  fin.read((char *)&rtmp,    sizeof(double));
  fin.read((char *)&itmp,    sizeof(double));
  omp01 = std::complex<double>(rtmp, itmp);
  fin.read((char *)&ll1,     sizeof(int)   );
  fin.read((char *)&mm1,     sizeof(int)   );
  fin.read((char *)&llmax1,  sizeof(int)   );
  fin.read((char *)&nnmax1,  sizeof(int)   );
  fin.read((char *)&nptsK1,  sizeof(int)   );
  fin.read((char *)&nptsE1,  sizeof(int)   );

  /* write run parameters */

  char buff[128];
  fin.read((char *)buff, 128);
  string mod_name(buff);

  if (
      (omp.real() != omp01.real())  ||
      (omp.imag() != omp01.imag())  ||
      (ll1 != l)                    ||
      (llmax1 != lmax)              ||
      (nnmax1 != nmax)              ||
      (nptsK1 != nptsK)             ||
      (nptsE1 != nptsE)             ||
      (mod_name != model->ModelID)
      ) return false;

/* Ok, read em in . . . */

  fin.read((char *)&ikap,  sizeof(int));

  for (int i=0; i<nmax; i++) {
    for (int j=0; j<nmax; j++) {
      for (int k=0; k<num_K; k++) {
	fin.read((char *)&rtmp, sizeof(double));
	fin.read((char *)&itmp, sizeof(double));
	mab[i*nmax+j][k] = std::complex<double>(rtmp,itmp);
	}
    }
  }

  return true;
}


void RespMat::put_chkpnt(void)
{
  std::string chkname;
  
  if (CHK_NAME.size()) {
				// Make unique check.point name
    std::ifstream fin(std::string(outdir+CHK_NAME).c_str());
    std::ostringstream sout;
    int i = 0;
    while(fin) {
      fin.close();
      sout.str("");
      sout << i;
      chkname = "" + outdir + CHK_NAME + "_" + sout.str();
      fin.open(chkname.c_str());
      if (i++>99) {
	cerr << "Too many CHECK POINT files!\n";
	exit(-1);
      }
    }
    fin.close();
  }
  

  ofstream fout(chkname.c_str());
  if (!fout) {
    cerr << "Couldn't open " << chkname << '\n';
    cerr << "continuing . . . \n";
    return;
  }

  fout.write((char const *)&omp,   sizeof(std::complex<double>));
  fout.write((char const *)&l,     sizeof(int)   );
  fout.write((char const *)&m,     sizeof(int)   );
  fout.write((char const *)&lmax,  sizeof(int)   );
  fout.write((char const *)&nmax,  sizeof(int)   );
  fout.write((char const *)&nptsK, sizeof(int)   );
  fout.write((char const *)&nptsE, sizeof(int)   );

  // write run parameters
  
  char buff[128];
  strncpy(buff, model->ModelID.c_str(), 127);
  fout.write((const char *)buff, 128);

/* Ok, write em out . . . */

  fout.write((char const *)&ikap,  sizeof(int)   );

  for (int i=0; i<nmax; i++) {
    for (int j=0; j<nmax; j++) {
      for (int k=0; k<num_K; k++) {
	double op = mab[i*nmax+j][k].real();
	double ip = mab[i*nmax+j][k].imag();
	fout.write((char const *)&op, sizeof(double));
	fout.write((char const *)&ip, sizeof(double));
	}
    }
  }
}


#define tol 1.0e-30
#define ZFRAC 0.001
#define ITMAX 25
#define NUMS 100
void RespMat::get_emax(void)
{
  if (!isotropic) {
    Emax = Emodmax;
    ESTEP = TOL1*(Emax - Emodmin);
    return;
  }

  int i;

  double ff;
  double dE = (Emodmax - Emodmin - 2*TOL1)/NUMS;
//  for (int i=0; i<=NUMS; i++) {  // 3/10/95 MDW
  for (i=1; i<=NUMS; i++) {
    Emax = Emodmin + dE*i + TOL1;
    if (E_max(Emax) < 0.5*tol) break;
  }

  if (i>=NUMS) {
    Emax = Emodmax;
  }
  else {

    double E0 = Emax - dE;
    double E1 = Emax;

    for (i=0; i<ITMAX; i++) {
      Emax = 0.5*(E0 + E1);
      ff = E_max(Emax);
      if (ff < 0.5*tol) E1 = Emax;
      else              E0 = Emax;
    }

    Emax = 0.5*(E1 + E0);

  }
  
  ESTEP = TOL1*(Emax - Emodmin);
}

/* Function to iteratively locate radius of circular orbit with energy EE */

static double Emdl;
static std::shared_ptr<AxiSymModel> mdl;

static double Eqcirc(double r)
{
  double ur, dudr, dif;
  
  mdl->get_pot_dpot(r, ur, dudr);
  dif =  Emdl - 0.5*r*dudr - ur;
  return dif;
}

double RespMat::E_max(double EE)
{
  double dudr, rc;

  Emdl = EE;
  mdl = model;

  rc = zbrent(Eqcirc,
	      ZFRAC*model->get_min_radius(), model->get_max_radius(), tol);
  dudr = model->get_dpot(rc);
  
  return model->distf(EE, kappa*sqrt(rc*rc*rc*dudr));
}

#undef tol


void RespMat::read_in(string& file)
{
  ifstream fin(file.c_str());
  if (!fin) bomb("Couldn't open <" + file + ">for input");

  read_in(fin);
}

void RespMat::write_out(string &file)
{
  ofstream fout(file.c_str());
  if (!fout) bomb("Couldn't open <" + file + ">for output");

  write_out(fout);
}

void RespMat::write_out(ostream &fout)
{
  if (!matrix_computed) {
    if (pv) make_matrix_pv();
    else make_matrix();
  }

  fout.write((char const *)&nmax, sizeof(int));

  for (int i=0; i<nmax; i++) {
    for (int j=0; j<nmax; j++) {
      double op = mr(i, j).real();
      double ip = mr(i, j).imag();
      fout.write((char const *)&op, sizeof(double));
      fout.write((char const *)&ip, sizeof(double));
    }
  }

  for (int i=0; i<nmax; i++) {
    for (int j=0; j<nmax; j++) {
      double op = mr2(i, j).real();
      double ip = mr2(i, j).imag();
      fout.write((char const *)&op, sizeof(double));
      fout.write((char const *)&ip, sizeof(double));
    }
  }

  fout.write((char const *)&omp,  sizeof(std::complex<double>));
  fout.write((char const *)&l,    sizeof(int));
  fout.write((char const *)&m,    sizeof(int));
  fout.write((char const *)&lmax, sizeof(int));
  fout.write((char const *)&nptsK, sizeof(int));
  fout.write((char const *)&nptsE, sizeof(int));
  fout.write((char const *)&pv,   sizeof(int));
  fout.write((char const *)&dof,  sizeof(int));
  fout.write((char const *)&isotropic, sizeof(int));
  fout.write((char const *)&SITYPE, sizeof(int));

  char s[128];
  unsigned i;
  for (i=0; i<ID.length(); i++) s[i] = ID[i];
  if (i<128) s[i] = '\0';
  fout.write((char const *)s, 128);
  
}

void RespMat::read_in(istream &fin)
{
  fin.read((char *)&nmax, sizeof(int));

  mr .resize(nmax, nmax);
  mr2.resize(nmax, nmax);

  for (int i=0; i<nmax; i++) {
    for (int j=0; j<nmax; j++) {
      fin.read((char *)&mr(i, j), sizeof(std::complex<double>));
    }
  }

  for (int i=0; i<nmax; i++) {
    for (int j=0; j<nmax; j++) {
      fin.read((char *)&mr2(i, j), sizeof(std::complex<double>));
    }
  }

  fin.read((char *)&omp, sizeof(std::complex<double>));
  fin.read((char *)&l, sizeof(int));
  fin.read((char *)&m, sizeof(int));
  fin.read((char *)&lmax, sizeof(int));
  fin.read((char *)&nptsK, sizeof(int));
  fin.read((char *)&nptsE, sizeof(int));
  fin.read((char *)&pv, sizeof(int));
  fin.read((char *)&dof, sizeof(int));
  fin.read((char *)&isotropic, sizeof(int));
  fin.read((char *)&SITYPE, sizeof(int));
  char id[128];
  fin.read((char *)id, 128);
  ID = id;
  
  matrix_computed = 1;
  matrix_computed_pv = 1;

}

std::complex<double> RespMat::disper(void)
{
  if (disp_computed) return disp;
  disp_computed = 1;

  Eigen::MatrixXcd ident(mr.rows(), mr.cols());
  ident.setIdentity();
  ident -= mr;

  disp = ident.determinant();

  return disp;
}


Eigen::VectorXcd RespMat::get_response (Eigen::VectorXcd& ext, gravity grav)
{
  Eigen::VectorXcd resp;

  switch (grav) {

  case none:

    resp = ext;
    break;

  case self:
    {
      if (!matrix_computed) {
	if (pv) make_matrix_pv();
	else make_matrix();
      }

      Eigen::MatrixXcd ident(mr.rows(), mr.cols());
      ident.setIdentity();
      ident -= mr;

#ifdef SVD_INVERSE      
      {
	Eigen::JacobiSVD<Eigen::MatrixXcd>
	  svd(ident, Eigen::ComputeThinU | Eigen::ComputeThinV);
	auto sv = svd.singularValues();
	auto isv = sv;
	int cnt = 0;
	for (int i=0; i<sv.size(); i++) {
	  if (fabs(sv[i]) > fabs(sv[0])*svd_tol) {
	    isv[i] = 1.0/sv[i];
	    cnt++;
	  }  else isv[i] = 0.0;
	}
	ident = svd.matrixV().transpose()*isv.transpose()*svd.matrixU();
	cout << "svd: " << cnt << " singular values" << std::endl;
      }
#else
      ident = ident.inverse();
#endif // SVD_INVERSE

      if (pv) 
	resp = (ident*mr2) * ext;
      else
	resp = (ident*mr) * ext;
    }
    break;

  case noself:

  default:

    if (!matrix_computed) {
      if (pv) make_matrix_pv();
      else make_matrix();
    }
    
    if (pv)
      resp = mr2 * ext;
    else
      resp = mr  * ext;
    break;

  }

  return resp;
}

Eigen::MatrixXd RespMat::get_wake
(Eigen::VectorXcd& ext,
 double xmin, double xmax, double ymin, double ymax,
 int numx, int numy,
 gravity grav, response type,
 double phi, double theta)
{
  Eigen::MatrixXd wake(numx, numy);
  wake.setZero();

  get_wake(wake, ext, xmin, xmax, ymin, ymax, numx, numy, grav, type, 
	    phi, theta);

  return wake;
}

void RespMat::get_wake (Eigen::MatrixXd& wake,
			Eigen::VectorXcd& ext,
			double xmin, double xmax, double ymin, double ymax,
			int numx, int numy,
			gravity grav, response type,
			double phi, double theta)
{

				// Transformation
  Eigen::Matrix3d return_euler(double PHI, double THETA, double PSI, int BODY);
  double plgndr(int l, int m, double x);

  if (biorth->get_dof()==2) {
    phi = 0.0;
    theta = 0.0;
  }

  double onedeg = M_PI/180.0;
  auto trans = return_euler(phi*onedeg, theta*onedeg, 0.0, 0);

  Eigen::Vector3d x0, xt;
  x0.setZero();
  xt.setZero();
				// Get response

  auto resp = get_response (ext, grav);
  
				// Compute wake

  int M = abs(m);
  double Ylm = sqrt( (0.5*l + 0.25)/M_PI * exp(lgamma(1.0+l-M) - lgamma(1.0+l+M)) );
  if (m>0 && (M != 2*(int)(M/2))) Ylm *= -1.0;

  double dx = (xmax - xmin)/(numx-1);
  double dy = (ymax - ymin)/(numy-1);
  double r, factr1=0, factr2=0;
  constexpr std::complex<double> I(0.0, 1.0);
  Eigen::VectorXd f;

  for (int i=0; i<numx; i++) {

    x0[1] = ymin + dy*i;

    for (int j=0; j<numy; j++) {
      x0[0] = xmin + dx*(j-1);
    
      xt = trans*x0;

      r = sqrt(xt.adjoint()*xt);
      if (r <= model->get_max_radius()) {
	phi = atan2(xt[1],xt[0]);
      
	f = ( exp(I*phi*static_cast<double>(m)) * resp ).real();

	if (biorth->get_dof()==2) factr1 = 1.0/sqrt(2.0*M_PI);
	if (biorth->get_dof()==3) factr1 = Ylm * plgndr(l, M, xt[3]/r);

	switch (type) {
	case density:
	  if (biorth->get_dof()==2) factr2 = 1.0/(2.0*M_PI);
	  if (biorth->get_dof()==3) factr2 = 1.0/(4.0*M_PI);
	  wake(i, j) += biorth->get_dens(r, l, f) * factr1 * factr2;

	  break;
	case potential:
	  wake(i, j) += biorth->get_potl(r, l, f) * factr1;
	  break;
	}
      }
    }
  }
      
}


Eigen::MatrixXd RespMat::get_bkgrnd (
			    double xmin, double xmax, double ymin, double ymax,
			    int numx, int numy, response type,
			    double phi, double theta)
{

				// Transformation
  Eigen::Matrix3d return_euler(double PHI, double THETA, double PSI, int BODY);


  if (biorth->get_dof()==2) {
    phi = 0.0;
    theta = 0.0;
  }

  double onedeg = M_PI/180.0;
  Eigen::Matrix3d trans = return_euler(phi*onedeg, theta*onedeg, 0.0, 0);

  Eigen::Vector3d x0, xt;
  x0.setZero();
  xt.setZero();
				// Compute background

  Eigen::MatrixXd wake(numx, numy);

  double dx = (xmax - xmin)/(numx-1);
  double dy = (ymax - ymin)/(numy-1);
  double r;
  Eigen::VectorXd f;

  for (int i=0; i<numx; i++) {

    x0[1] = ymin + dy*i;

    for (int j=0; j<numy; j++) {
      x0[0] = xmin + dx*j;
    
      xt = trans*x0;
      r = sqrt(xt.adjoint()*xt);

      if (r > model->get_max_radius())
	wake(i, j) = 0.0;
      else {
	phi = atan2(xt[2], xt[0]);
      
	switch (type) {
	case density:
	  wake(i, j) = model->get_density(xt[0], xt[1], xt[2]);
	  break;
	case potential:
	  wake(i, j) = model->get_pot(xt[0], xt[1], xt[2]);
	  break;
	}
      }

    }
      
  }

  return wake;
}


void RespMat::print_wake_volume (ostream& out, Eigen::VectorXcd& ext,
				 double xmin, double xmax, 
				 double ymin, double ymax,
				 double zmin, double zmax,
				 int numx, int numy, int numz, 
				 bool bkgrnd, gravity grav, response type)
{
  double plgndr(int l, int m, double x);
  Eigen::Vector3d xt;

				// Get response
  Eigen::VectorXcd resp = get_response (ext, grav);

				// Compute wake

  int M = abs(m);
  double Ylm = sqrt( (0.5*l + 0.25)/M_PI * exp(lgamma(1.0+l-M) - lgamma(1.0+l+M)) );
  if (m>0 && (M != 2*(int)(M/2))) Ylm *= -1.0;

  double dx = (xmax - xmin)/(numx-1);
  double dy = (ymax - ymin)/(numy-1);
  double dz = (zmax - zmin)/(numz-1);
  double r, costh, phi;
  std::complex<double> I(0.0, 1.0);
  Eigen::VectorXd f;
  float z;

  out.write((const char *)&numx, sizeof(int));
  out.write((const char *)&numy, sizeof(int));
  out.write((const char *)&numz, sizeof(int));

  out.write((const char *)&(z=xmin), sizeof(float));
  out.write((const char *)&(z=xmax), sizeof(float));
  out.write((const char *)&(z=ymin), sizeof(float));
  out.write((const char *)&(z=ymax), sizeof(float));
  out.write((const char *)&(z=zmin), sizeof(float));
  out.write((const char *)&(z=zmax), sizeof(float));

  
  for (int k=0; k<numz; k++) {
    xt[2] = zmin + dz*k;

    for (int i=0; i<numy; i++) {
      xt[1] = ymin + dy*i;

      for (int j=0; j<numx; j++) {
	xt[0] = xmin + dx*j;
    
	r = sqrt(xt.adjoint()*xt);
	if (r > model->get_max_radius())
	  z = 0.0;
	else {
	  phi = atan2(xt[1], xt[0]);
	  costh = xt[2]/r;
      
	  f = ( exp(I*phi*static_cast<double>(m)) * resp ).real();

	  switch (type) {
	  case density:
	    z = biorth->get_dens(r, l, f);
	    if (biorth->get_dof()==2) z *= 1.0/(2.0*M_PI);
	    if (biorth->get_dof()==3) z *= 1.0/(4.0*M_PI);
	    break;
	  case potential:
	    z = biorth->get_potl(r, l, f);
	    break;
	  }

	  z *= Ylm * plgndr(l, M, costh);

	  if (bkgrnd) {
	    switch (type) {
	    case density:
	      z /= model->get_density(xt[1], xt[2], xt[3]);
	      break;
	    case potential:
	      z /= model->get_pot(xt[1], xt[2], xt[3]);
	      break;
	    }
	  }
	}

	out.write((const char *)&z, sizeof(float));

      }
    }
  }

}


void RespMat::get_wake_volume (std::vector<Eigen::MatrixXd>& mat,
			       Eigen::VectorXcd& ext,
			       double xmin, double xmax, 
			       double ymin, double ymax,
			       double zmin, double zmax,
			       int numx, int numy, int numz, 
			       bool bkgrnd, gravity grav, response type)
{
  double plgndr(int l, int m, double x);
  Eigen::Vector3d xt;

				// Get response
  Eigen::VectorXcd resp = get_response (ext, grav);

				// Compute wake

  int M = abs(m);
  double Ylm = sqrt( (0.5*l + 0.25)/M_PI * exp(lgamma(1.0+l-M) - lgamma(1.0+l+M)) );
  if (m>0 && (M != 2*(int)(M/2))) Ylm *= -1.0;

  double dx = (xmax - xmin)/(numx-1);
  double dy = (ymax - ymin)/(numy-1);
  double dz = (zmax - zmin)/(numz-1);
  double z=0, r, costh, phi;
  constexpr std::complex<double> I(0.0, 1.0);
  Eigen::VectorXd f;
  
  for (int k=0; k<numz; k++) {
    xt[2] = zmin + dz*k;

    for (int i=0; i<numy; i++) {
      xt[1] = ymin + dy*i;

      for (int j=0; j<numx; j++) {
	xt[0] = xmin + dx*j;
    
	r = sqrt(xt.adjoint()*xt);
	if (r > model->get_max_radius())
	  z = 0.0;
	else {
	  phi = atan2(xt[1],xt[0]);
	  costh = xt[2]/r;
      
	  f = ( exp(I*phi*static_cast<double>(m)) * resp ).real();

	  switch (type) {
	  case density:
	    z = biorth->get_dens(r, l, f);
	    if (biorth->get_dof()==2) z *= 1.0/(2.0*M_PI);
	    if (biorth->get_dof()==3) z *= 1.0/(4.0*M_PI);
	    break;
	  case potential:
	    z = biorth->get_potl(r, l, f);
	    break;
	  }

	  z *= Ylm * plgndr(l, M, costh);

	  if (bkgrnd) {
	    switch (type) {
	    case density:
	      z /= model->get_density(xt[1], xt[2], xt[3]);
	      break;
	    case potential:
	      z /= model->get_pot(xt[1], xt[2], xt[3]);
	      break;
	    }
	  }
	}

	mat[k](i, j) += z;

      }
    }
  }

}



void write_contr(ostream& out,
		 double xmin, double xmax, double ymin, double ymax, 
		 Eigen::MatrixXd& mat)
{
  float z;

  int numx = mat.cols();
  int numy = mat.rows();

  out.write((const char *)&numx, sizeof(int));
  out.write((const char *)&numy, sizeof(int));

  out.write((const char *)&(z=xmin), sizeof(float));
  out.write((const char *)&(z=xmax), sizeof(float));
  out.write((const char *)&(z=ymin), sizeof(float));
  out.write((const char *)&(z=ymax), sizeof(float));

  for (int i=0; i<numy; i++) 
    for (int j=0; j<numx; j++)
      out.write((const char *)&(z=mat(i, j)), sizeof(float));
}


void write_volume_hips (ostream& out,
			double xmin, double xmax, 
			double ymin, double ymax,
			double zmin, double zmax,
			int numx, int numy, int numz,
			std::vector<Eigen::MatrixXd>& mat)
{
  float z;

  out << "output by respmat3" << endl;
  out << "frames (z-axis)" << endl;
  out << numz << endl;
  out << "rows, columns" << endl;
  out << numx << endl;
  out << numy << endl;
  out << 32 << endl;
  out << 0 << endl;
  out << 3 << endl;
  out << "history" << endl;
  out << "comment" << endl;
  out << "." << endl;

  for (int k=0; k<numz; k++) {

    for (int i=0; i<numy; i++) {

      for (int j=0; j<numx; j++) {

	z = mat[k](i, j);
	out.write((const char *)&z, sizeof(float));

      }
    }
  }
}


typedef struct {
  double  x;
  double  y;
  double  z;
  double  nodeValue;
  int valid;
} node;


void write_volume_cubes (ostream& out,
			 double xmin, double xmax, 
			 double ymin, double ymax,
			 double zmin, double zmax,
			 int numx, int numy, int numz,
			 std::vector<Eigen::MatrixXd>& mat)
{
  long l;

  out.write((const char *)&(l=numx), sizeof(long));
  out.write((const char *)&(l=numy), sizeof(long));
  out.write((const char *)&(l=numz), sizeof(long));
  
  double dx = (xmax - xmin)/(numx-1);
  double dy = (ymax - ymin)/(numy-1);
  double dz = (zmax - zmin)/(numz-1);
  node n;
  n.valid = 1;
  
  for (int k=0; k<numz; k++) {
    n.z = zmin + dz*(k-1);

    for (int i=0; i<numy; i++) {
      n.y = ymin + dy*(i-1);

      for (int j=0; j<numx; j++) {
	n.x = xmin + dx*(j-1);
	n.nodeValue = mat[k](i, j);

	out.write((const char *)&n, sizeof(node));
      }
    }
  }
}



void RespMat::print_wake (ostream& out, Eigen::VectorXcd& ext,
			  double xmin, double xmax, double ymin, double ymax,
			  int numx, int numy,
			  gravity grav, response type,
			  double phi, double theta)
{
  Eigen::MatrixXd mat = get_wake (ext, xmin, xmax, ymin, ymax, numx, numy,
		  grav, type, phi, theta);

  write_contr (out, xmin, xmax, ymin, ymax, mat);
}


void RespMat::print_wake_volume_hips 
	(ostream& out, Eigen::VectorXcd& ext,
	 double xmin, double xmax, 
	 double ymin, double ymax,
	 double zmin, double zmax,
	 int numx, int numy, int numz, bool bkgrnd,
	 gravity grav, response type)
{

  std::vector<Eigen::MatrixXd> mat(numz);
  for (auto & m : mat) m.resize(numy, numx);

  get_wake_volume(mat, ext, xmin, xmax, ymin, ymax, zmin, zmax,
		  numx, numy, numz, bkgrnd, grav, type);

  write_volume_hips(out, xmin, xmax, ymin, ymax, zmin, zmax, 
		    numx, numy, numz, mat);
}

void RespMat::print_wake_volume_cubes
	(ostream& out, Eigen::VectorXcd& ext,
	 double xmin, double xmax, 
	 double ymin, double ymax,
	 double zmin, double zmax,
	 int numx, int numy, int numz, bool bkgrnd,
	 gravity grav, response type)
{

  std::vector<Eigen::MatrixXd> mat(numz);
  for (auto & m : mat) m.resize(numy, numx);

  get_wake_volume(mat, ext, xmin, xmax, ymin, ymax, zmin, zmax,
		  numx, numy, numz, bkgrnd, grav, type);

  write_volume_cubes(out, xmin, xmax, ymin, ymax, zmin, zmax, 
		    numx, numy, numz, mat);
}

void MatrixXcdSynchronize(Eigen::MatrixXcd& mat, int id);
void ComplexSynchronize(std::complex<double>& c, int id);

void RespMat::MPISynchronize(int id)
{
  MatrixXcdSynchronize(mr,  id);
  MatrixXcdSynchronize(mr2, id);
  ComplexSynchronize(omp, id);

  MPI_Bcast(&matrix_computed, 1, MPI_INT, id, MPI_COMM_WORLD);
  MPI_Bcast(&matrix_computed_pv, 1, MPI_INT, id, MPI_COMM_WORLD);
}
