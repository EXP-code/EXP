#define ID_STRING

#include <cstdlib>
#include <sys/time.h>
#include <sys/resource.h>
#include <values.h>

#include <iostream>
#include <fstream>
#include <sstream>
#include <map>
#include <string>

#include <kevin_complex.h>
#include <Vector.h>
#include <orbit.h>
#include <massmodel.h>
#include <biorth2d.h>
#include <clinalg.h>
#include <interp.h>

#include <RespMat3.h>
#include <model3d.h>
#include <sphereSL.h>
#include <UnboundCoefs.H>

#include "expand.h"

#include <UserSatWake.H>


inline bool isEven (int i)
{
  if ((i%2) == 0) return true;
  else return false;
}


inline bool isOdd (int i)
{
  if ((i%2) == 1) return true;
  else return false;
}

inline IntegrationType ITOSIT(int i)
{
  switch (i) {
  case ratint:
    return ratint;
  case jacoint:
    return jacoint;
  default:
    return rombint;
  }
}

double rot_matrix(int l, int m, int n, double beta);
double plgndr(int, int, double);


UserSatWake::UserSatWake(string &line) : ExternalForce(line)
{
  id = "SatelliteWake";

  Number = 8;
  Pmax = 40;
  LMIN = 0;
  LMAX = 5;
  M = 2;
  lmax = 4;
  nmax = 10;
  HALO_TRUNC = -1;
  nptsE = 5;
  nptsK = 4;
  CAUCHY = 0;
  RATINT = 0;
  PTGRID = 200;
  NRECS = 512;
  DIVERGE = 0;
  DIVEXPON = 1.0;
  OLD = 0;
  VERBOSE = 0;
  HALO_TYPE = 3;
  SITYPE = ITOSIT(2);
  RMODMAX = 20.0;
  RGRID = 1000;
  DELTA = 0.0;
  OMPI = 0.01;
  OMPS = -1;
  NUMDF = 100;
  RA = 1.0e20;
  INCLINE = 45.0;
  PSI = 0.0;
  PHIP = 0.0;
  TYPE = 3;
  MODEL = 0;
  NUMT = 20;
  DIVERGE = 0;
  DIVEXPON = 1.0;
  E = 0.0;
  Rperi = 0.1;
  Redge = 2.0;
  rmin = -1.0;
  rmax = -1.0;
  scale = 0.06667;
  numr = 400;
  deltaR = 0.01;
  deltaT = 0.01;
  Toffset = 0.0;
  satmass = 0.1;
  INFILE = "SLGridSph.model";
  CACHEDIR = "./";

  I = KComplex(0.0, 1.0);
  rtol = 1.0e-2;

  ctr_name = "";		// Default component for center

  cachefile = string(".halo_response.") + runtag;

  initialize();
							      
  if (ctr_name.size()>0) {
				// Look for the fiducial component for
				// centering
    bool found = false;
    list<Component*>::iterator cc;
    Component *c;
    for (cc=comp.components.begin(); cc != comp.components.end(); cc++) {
      c = *cc;
      if ( !ctr_name.compare(c->name) ) {
	c0 = c;
	found = true;
      break;
      }
    }

    if (!found) {
      cerr << "Process " << myid << ": can't find desired component <"
	   << ctr_name << ">" << endl;
      MPI_Abort(MPI_COMM_WORLD, 35);
    }

  }
  else
    c0 = NULL;


  initialize_coefficients();

  userinfo();
}

UserSatWake::~UserSatWake()
{
  
  delete halo_model;
  delete u;
}

void UserSatWake::userinfo()
{
  if (myid) return;		// Return if node master node

  print_divider();

  cout << "** User routine SATWAKE initialized, " ;
  cout << "Filename=" << INFILE << "  DIVERGE=" << DIVERGE
       << "  DIVEXPON=" << DIVEXPON;
  if (c0) 
    cout << ", center on component <" << ctr_name << ">";
  else
    cout << ", using inertial center";
  cout << setw(50) << setfill('-') << '-' << endl << setfill(' ');
  cout << "Main parameters:" << endl;
  cout << setw(9) << "" << setw(15) << "NUMT"     << " = " << NUMT     << endl;
  cout << setw(9) << "" << setw(15) << "LMIN"     << " = " << LMIN     << endl;
  cout << setw(9) << "" << setw(15) << "LMAX"     << " = " << LMAX     << endl;
  cout << setw(9) << "" << setw(15) << "M"        << " = " << M        << endl;
  cout << setw(9) << "" << setw(15) << "lmax"     << " = " << lmax     << endl;
  cout << setw(9) << "" << setw(15) << "nmax"     << " = " << nmax     << endl;
  cout << setw(9) << "" << setw(15) << "INCLINE"  << " = " << INCLINE  << endl;
  cout << setw(9) << "" << setw(15) << "PSI"      << " = " << PSI      << endl;
  cout << setw(9) << "" << setw(15) << "PHIP"     << " = " << PHIP     << endl;
  cout << setw(9) << "" << setw(15) << "E"        << " = " << E        << endl;
  cout << setw(9) << "" << setw(15) << "Rperi"    << " = " << Rperi    << endl;
  cout << setw(9) << "" << setw(15) << "Redge"    << " = " << Redge    << endl;
  cout << setw(9) << "" << setw(15) << "rmin"     << " = " << rmin     << endl;
  cout << setw(9) << "" << setw(15) << "rmax"     << " = " << rmax     << endl;
  cout << setw(9) << "" << setw(15) << "scale"    << " = " << scale    << endl;
  cout << setw(9) << "" << setw(15) << "numr"     << " = " << numr     << endl;
  cout << setw(9) << "" << setw(15) << "MASS"     << " = " << satmass  << endl;
  cout << setw(50) << setfill('-') << '-' << endl << setfill(' ');

  I = KComplex(0.0, 1.0);
  rtol = 1.0e-2;

  ctr_name = "";		// Default component for center



  cout << endl;

  print_divider();
}


void UserSatWake::initialize()
{
  string val;

  if (get_value("Number", val))		Number = atoi(val);
  if (get_value("Pmax", val))		Pmax = atoi(val);
  if (get_value("LMIN", val))		LMIN = atoi(val);
  if (get_value("LMAX", val))		LMAX = atoi(val);
  if (get_value("M", val))		M = atoi(val);
  if (get_value("lmax", val))		lmax = atoi(val);
  if (get_value("nmax", val))		nmax = atoi(val);
  if (get_value("numr", val))		numr = atoi(val);
  if (get_value("rmin", val))		rmin = atof(val);
  if (get_value("rmax", val))		rmax = atof(val);
  if (get_value("scale", val))		scale = atof(val);
  if (get_value("HALO_TRUNC", val))	HALO_TRUNC = atoi(val);
  if (get_value("nptsE", val))		nptsE = atoi(val);
  if (get_value("nptsK", val))		nptsK = atoi(val);
  if (get_value("CAUCHY", val))		CAUCHY = atoi(val);
  if (get_value("RATINT", val))		RATINT = atoi(val);
  if (get_value("HALO_TYPE", val))	HALO_TYPE = atoi(val);
  if (get_value("PTGRID", val))		PTGRID = atoi(val);
  if (get_value("NRECS", val))		NRECS = atoi(val);
  if (get_value("DIVERGE", val))	DIVERGE = atoi(val);
  if (get_value("DIVEXPON", val))	DIVEXPON = atoi(val);
  if (get_value("OLD", val))		OLD = atoi(val);
  if (get_value("VERBOSE", val))	VERBOSE = atoi(val);
  if (get_value("SITYPE", val))		SITYPE = ITOSIT(atoi(val));
  if (get_value("Rperi", val))		Rperi = atof(val);
  if (get_value("Redge", val))		Redge = atof(val);
  if (get_value("deltaR", val))		deltaR = atof(val);
  if (get_value("deltaT", val))		deltaT = atof(val);
  if (get_value("Toffset", val))	Toffset = atof(val);
  if (get_value("MASS", val))		satmass = atof(val);
  if (get_value("INFILE", val))		INFILE = val;
  if (get_value("CACHEDIR", val))	CACHEDIR = val;
  if (get_value("ctrname", val))	ctr_name = val;
}


void UserSatWake::initialize_coefficients()
{
  //==========================
  // Begin normal exectution
  //==========================
  
  const int id  = 0x11;
  KComplex       omp_halo;
  double	 ompr;
  
  
  //=================
  // Sanity checks
  //=================
  
  if (M<0 && myid==0) cerr << "Warning: M should be >0 . . . fixing\n";
  M = abs(M);
  if (M>LMIN) {
    if (myid==0) cout << "Setting LMIN=M=" << M << endl;
    LMIN = M;
  }
  
  
  // ===================================================================
  // Initilize HALO model
  // ===================================================================
  
  SphericalModelTable::logscale = 0;
  SphericalModelTable *m = 0;
  
  m = new SphericalModelTable(INFILE, DIVERGE, DIVEXPON);
  m->setup_df(NUMDF, RA);
  halo_model = m;
  Model3dNames[0] = INFILE;	// Assign filename to ID string
  
  // ===================================================================
  // Initilize biorthogonal functions
  // ===================================================================
  
  switch (HALO_TYPE) {
  case bessel:
    u = new BSSphere(RMODMAX, nmax, LMAX);
    break;
  case clutton_brock:
    u = new CBSphere;
    break;
  case hernquist:
    u = new HQSphere;
    break;
  case sturm:
    if (rmin<0.0) 
      rmin = halo_model->get_min_radius() * 1.2;
    if (rmax<0.0) 
      rmax = halo_model->get_max_radius() * 0.99;

    SLGridSph::sph_cache_name = ".slgrid_sph_cache." + runtag;
    SphereSL::mpi = 1;
    u = new SphereSL(LMAX, nmax, numr, rmin, rmax, scale, m);
    break;
  default:
    if (myid==0) 
      cerr << "Illegal spherical biorthongal series: " << HALO_TYPE << '\n';
    exit(-1);
  }
  
  
  // ===================================================================
  // Setup orbit
  // ===================================================================
  
  
  UnboundCoefs Coefs(E, Rperi, Redge, deltaR, halo_model);
  
  // Compute frequency vector
  double Tmax = Coefs.Tmax();
  double dOmega = 2.0*M_PI/(4.0*Tmax);
  
  vector< vector<RespMat> > total;

  int nfreqs = 2*Pmax+1;
  vector<double> freqs(nfreqs);
  vector<KComplex> Omega(nfreqs);
  for (int k=-Pmax; k<=Pmax; k++) freqs[k+Pmax] = dOmega*k;
  
  Times = vector<double>(NUMT+1);
  for (int nt=0; nt<=NUMT; nt++) Times[nt] = -Tmax + 2.0*Tmax*nt/NUMT;
  
  
  map< int, map<int, vector<CMatrix> > > coefs;
  CVector tcoefs;
  
  if (myid==0) cout << "Computing Laplace coefficients . . ." << endl;
  for (int L=LMIN; L<=LMAX; L++) {
    if (myid==0) cout << "L=" << L << endl;
    for (int L2=(isEven(L)?0:1);L2<=L; L2++) {
      Coefs.coefs(L, L2, nmax, Times, deltaT, u, coefs[L][L2], freqs);
      if (myid==0)
	cout << "    L2=" << L2 << " r=" << coefs[L][L2][0].getrlow()
	     << ", " << coefs[L][L2][0].getrhigh() << endl;
    }
  }
  if (myid==0) cout << "done" << endl;
  
  // Factor for point mass expansion for
  // orbit in x-y plane
  
  INCLINE *= M_PI/180.0;
  PSI     *= M_PI/180.0;
  PHIP    *= M_PI/180.0;
  
  
  // ===================================================================
  // Get combined response matrix for disk and halo for each L
  // ===================================================================
  
  
  ofstream to_save;
  ifstream from_save;
  unsigned short reading = 1;

  if (myid==0) {
    
    from_save.open(cachefile.c_str());
    if (!from_save) {
      if (myid==0) {
	cerr << "Couldn't open <" << cachefile <<
	  "> to read cached data" << endl;
	cerr << ". . .  will attempt to write a cache file." << endl;
      }
      reading = false;

    } else {
    
#ifdef ID_STRING
      // Read ID string and check version #
      int tid;
      from_save.read((char *)&tid, sizeof(int));
      if (tid != id) {
	
	cerr << "Incompatible save file!\n";
	MPI_Abort(MPI_COMM_WORLD, -34);
	
      } else {
    
	char tbuf[255];
	from_save.read(tbuf, 255);
	cout << tbuf << endl;
#endif // ID_STRING
      }
    }

    if (!reading) {

      to_save.open(cachefile.c_str());
      if (!to_save) {
	if (myid==0)
	  cerr << "Couldn't open <" << cachefile <<
	    "> to write cached data" << endl;
	MPI_Abort(MPI_COMM_WORLD, -35);
	exit(-1);
      }
    
#ifdef ID_STRING
      // Write ID string and version #
      time_t tp = time(NULL);
      string stime = string("UserSatWake") + ": " + ctime(&tp);
      cout << "ID string: " << stime<< endl;
    
      to_save.write((const char *)&id, sizeof(int));
      char tbuf[255];
      strncpy(tbuf, stime.c_str(), 255);
      to_save.write(tbuf, 255);
      
#endif // ID_STRING
    }
  }

  MPI_Bcast(&reading, 1, MPI_UNSIGNED_SHORT, 0, MPI_COMM_WORLD);
  
  if (myid && reading) {

    from_save.open(cachefile.c_str());
    if (!from_save) {
      cerr << "Process " << myid << " couldn't open <" << cachefile <<
	"> to read cached data!" << endl;
      MPI_Abort(MPI_COMM_WORLD, -32);
    } else {
    
#ifdef ID_STRING
      // Read ID string and check version #
      int tid;
      from_save.read((char *)&tid, sizeof(int));
      if (tid != id) {
	
	cerr << "Process " << myid << " claims an incompatible save file!"
	     << endl;
	MPI_Abort(MPI_COMM_WORLD, -34);
	
      } else {
    
	char tbuf[255];
	from_save.read(tbuf, 255);
#endif // ID_STRING
      }
    }
  }
  
  
  KComplex I(0.0, 1.0);
  string CACHE = CACHEDIR + "/" + runtag;
  
  for (int L=LMIN; L<=LMAX; L++) Lhalo.push_back(L);
  Nhalo = Lhalo.size();
  vector<RespMat> mab_halo (Nhalo);
  
  for (int nf=0; nf<nfreqs; nf++) {
    
    // Frequency from satellite in halo
    ompr = freqs[nf];
    omp_halo = Omega[nf] = KComplex(ompr, OMPI*OMPS) * OMPS;
    
    if (myid==0) cout << "Halo frequency: " << omp_halo << endl;
    
    
    if (reading) {
      
      for (int ihalo=0; ihalo<Nhalo; ihalo++)
	mab_halo[ihalo] = RespMat(from_save, halo_model, u);
      
      total.push_back(mab_halo);
      
    }
    else {
      
      
      // Get response matrices for each halo l-order
      //
      for (int ihalo=0; ihalo<Nhalo; ihalo++) {
	
	mab_halo[ihalo] = RespMat(omp_halo, Lhalo[ihalo], M, 
				  lmax, nmax, nptsK, nptsE,
				  halo_model, u, OLD, 0, 0,
				  VERBOSE, SITYPE, CAUCHY, RATINT);
	
	mab_halo[ihalo].set_params(DELTA, PTGRID);
	
	if (CAUCHY)
	  mab_halo[ihalo].make_matrix_pv();
	else 
	  mab_halo[ihalo].make_matrix();
	
#ifdef DEBUG
	if (myid==0) 
	  cout << "Main: finished ihalo=" << ihalo << "/" << Nhalo 
	       << endl << flush;
#endif
      }
      
      total.push_back(mab_halo);
      
      // Save to cache
      //
      if (myid==0) {
	for (int ihalo=0; ihalo<Nhalo; ihalo++)
	  mab_halo[ihalo].write_out(to_save);
      }
    }
    
  }
  
  // ===================================================================
  // For paranoia's sake
  // ===================================================================
  if (reading)
    from_save.close();
  else
    if (myid==0) to_save.close();
  
  // ===================================================================
  // Make movie!
  // ===================================================================
  
  ofstream tlog((runtag  + ".satwake.timelog").c_str());
  
  CVector tmp(1, nmax); tmp.zero();

  rcoefs = vector< vector<CVector> >(NUMT+1);
  for (int nt=0; nt<=NUMT; nt++) {
    for (int ihalo=0; ihalo<Nhalo; ihalo++)
      rcoefs[nt].push_back(tmp);
  }

  curcoefs = vector<CVector>(Nhalo);
  for (int ihalo=0; ihalo<Nhalo; ihalo++)
    curcoefs.push_back(tmp);


  for (int nt=0; nt<=NUMT; nt++) {
    
    tlog << setw(8) << nt << setw(16) << Times[nt] << endl;
    
    //
    // Frequency loop for current time
    //
    for (int nf=0; nf<nfreqs; nf++) {
      
      //
      // Get satellite perturbation coefficients
      //
      for (int ihalo=0; ihalo<Nhalo; ihalo++) {
	
	int L = Lhalo[ihalo];
	
	for (int L2=(isEven(L)?0:1);L2<=L; L2++) {
	  // 
	  // Truncate satellite coefficients
	  //
	  tcoefs = coefs[L][L2][nt][nf];
	  if (HALO_TRUNC>0) {
	    for (int n=HALO_TRUNC+1; n<=nmax; n++) tcoefs[n] = 0.0;
	  }
	  
	  
	  //
	  // Factor for norm in density component of biorthogonal pair
	  //         |
	  //         |          Account for pos and neg M values 
	  //         |          |
	  //         |          |     Factor from spherical harmonic
	  //         |          |     |
	  //         |          |     |
	  //         V          V     V
	  tcoefs *= 4.0*M_PI * 2.0 * sqrt( (0.5*L + 0.25)/M_PI *
            exp(lgamma(1.0+L-abs(L2)) - lgamma(1.0+L+abs(L2))) ) *
            plgndr(L, abs(L2), 0.0) * 
            rot_matrix(L, M, L2, INCLINE) * exp(I*L2*PSI)*exp(I*M*PHIP) ;
	  //                  ^                    ^
	  //                  |                    |
	  //  Rotate component to spherical harmonic in new orientation
	  //
	  
	  
	  //
	  // Leading sign for spherical harmonic with m>0 (following Edmonds
	  // to be consistent with definition of rotation matrices)
	  
	  if (L2>0 && isOdd(L2)) tcoefs *= -1.0;
	}

	tcoefs *= dOmega/(2.0*M_PI) * exp(I*Omega[nf]*Times[nt]);
	
	rcoefs[nt][ihalo] 
	  += total[nf][ihalo].get_response(tcoefs, RespMat::self);
      }
    }
  }
  
  

}


void UserSatWake::determine_acceleration_and_potential(void)
{
  compute_coefficients();
  exp_thread_fork(false);
}


void * UserSatWake::determine_acceleration_and_potential_thread(void * arg) 
{
  int id = *((int*)arg);
  unsigned nbodies = cC->Number();
  int nbeg = nbodies*id/nthrds;
  int nend = nbodies*(id+1)/nthrds;

  double pos[3], rr, r, x, y, z, phi, costh, pot0, pot1, dpot;
  double potr, potl, pott, potp, R2, fac, Ylm, dYlm;
  Vector f, g;
  
  map<unsigned long, Particle>::iterator it = cC->Particles().begin();

  for (int q=0   ; q<nbeg; q++) it++;
  for (int q=nbeg; q<nend; q++) {

    unsigned long indx = (it++)->first;
				// If we are multistepping, compute accel 
				// only at or below this level

    if (multistep && (cC->Part(indx)->level < mlevel)) continue;

    rr = 0.0;
    for (int k=0; k<3; k++) {
      pos[k] = cC->Pos(indx, k);	// Inertial by default
      if (c0) pos[k] -= c0->center[k];
      rr += pos[k]*pos[k];
    }
    r = sqrt(rr);

    potr = potl = pott = potp = 0.0;

    if (r <= halo_model->get_max_radius()) {
      
      x = pos[0];
      y = pos[1];
      z = pos[2];

      phi = atan2(y, x);
      costh = z/(r+DSMALL);
      R2 = x*x + y*y;
      
      for (unsigned ihalo=0; ihalo<Lhalo.size(); ihalo++) {

	int L = Lhalo[ihalo];

	fac = sqrt( (0.5*L + 0.25)/M_PI * 
		    exp(lgamma(1.0+L-M) - lgamma(1.0+L+M)) ) * satmass;
	Ylm = fac * plgndr(L, M, costh);
	dYlm = fac * dplgndr(L, M, costh);

	f = Re( exp(I*phi*M) * curcoefs[ihalo] );
	g = Re( I*M*exp(I*phi*M) * curcoefs[ihalo] );
	  
	pot0 = u->get_potl(r, L, f);
	pot1 = u->get_potl(r, L, g);
	dpot = (u->get_potl(r*(1.0+rtol), L, f) - 
		u->get_potl(r*(1.0-rtol), L, f) ) / ( 2.0*rtol*r );
	
	potl += Ylm  * pot0;
	potr += Ylm  * dpot;
	pott += dYlm * pot0;
	potp += Ylm  * pot1;
      }
      
      cC->AddAcc(indx, 0, -(potr*x/r - pott*x*z/(r*r*r)) );
      cC->AddAcc(indx, 1, -(potr*y/r - pott*y*z/(r*r*r)) );
      cC->AddAcc(indx, 2, -(potr*z/r + pott*R2 /(r*r*r))   );
      
      if (R2 > DSMALL) {
	cC->AddAcc(indx, 0,  potp*y/R2 );
	cC->AddAcc(indx, 1, -potp*x/R2 );
      }
      if (use_external)
	cC->AddPotExt(indx, potl);
      else
	cC->AddPot(indx, potl);
    }
  }

  return (NULL);
}


void UserSatWake::compute_coefficients()
{
  int indx;
  double a, b, curT = tnow - Toffset;

  // Find time index

  if (curT<Times.front()) {
    indx = 0;
    // a = 1.0;
    a = 0.0;
    b = 0.0;
  } else if (curT>Times.back()) {
    indx = Times.size()-2;
    a = 0.0;
    // b = 1.0;
    b = 0.0;
  } else {
    indx = Vlocate(curT, Times);
    a = (Times[indx+1] - tnow)/(Times[indx+1] - Times[indx]);
    b = (tnow - Times[indx+1])/(Times[indx+1] - Times[indx]);
  }

  CVector tmp(1, nmax); tmp.zero();
  rcoefs = vector< vector<CVector> >(NUMT+1);
  for (int nt=0; nt<=NUMT; nt++) {
    for (int ihalo=0; ihalo<Nhalo; ihalo++)
      rcoefs[nt].push_back(tmp);
  }

  curcoefs = vector<CVector>(Nhalo);
  for (int ihalo=0; ihalo<Nhalo; ihalo++)
    curcoefs[ihalo] = a*rcoefs[indx][ihalo] + b*rcoefs[indx+1][ihalo];
}


extern "C" {
  ExternalForce *satWake(string& line)
  {
    return new UserSatWake(line);
  }
}

class proxysatwake { 
public:
  proxysatwake()
  {
    factory["usersatwake"] = satWake;
  }
};

proxysatwake p;
