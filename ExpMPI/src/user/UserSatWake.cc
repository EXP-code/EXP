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
#include <TimeSeriesCoefs.H>

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

  LMIN 		= 0;
  LMAX 		= 4;
  MMIN 		= 0;
  MMAX 		= 4;
  lmax 		= 4;
  nmax 		= 10;
  nfreqs 	= 1200;
  HALO_TRUNC 	= -1;
  nptsE 	= 5;
  nptsK 	= 4;
  CAUCHY 	= 0;
  RATINT 	= 0;
  PTGRID 	= 200;
  NRECS 	= 512;
  DIVERGE 	= 0;
  DIVEXPON 	= 1.0;
  OLD		= 0;
  VERBOSE 	= 0;
  HALO_TYPE 	= 3;
  SITYPE 	= ITOSIT(2);
  RMODMAX 	= 20.0;
  DELTA 	= 0.0;
  OMPI 		= -0.01;
  NUMDF 	= 100;
  RA 		= 1.0e20;
  INCLINE 	= 45.0;
  PSI 		= 0.0;
  PHIP 		= 0.0;
  TYPE		= 3;
  MODEL 	= 0;
  NUMT		= 20;
  DIVERGE 	= 0;
  DIVEXPON 	= 1.0;
  E		= 0.0;
  Rperi 	= 0.1;
  Rfac 		= 0.1;
  Mfac 		= 2.0;
  rmin 		= -1.0;
  rmax 		= -1.0;
  scale 	= 0.06667;
  numr 		= 400;
  nint 		= 12;
  Tmax 		= 4.0;
  delT 		= 0.01;
  Toffset 	= 0.0;
  satmass 	= 0.1;
  INFILE 	= "SLGridSph.model";
  CACHEDIR 	= "./";
  UseCache 	= true;

  I 		= KComplex(0.0, 1.0);
  rtol 		= 1.0e-2;

  // Default component for center
  //
  ctr_name	= "";
  cachefile 	= string(".halo_response.") + runtag;

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
  cout << endl << setw(50) << setfill('-') << '-' << endl << setfill(' ');
  cout << "Main parameters:" << endl;
  cout << setw(9) << "" << setw(15) << "NUMT"     << " = " << NUMT     << endl;
  cout << setw(9) << "" << setw(15) << "LMIN"     << " = " << LMIN     << endl;
  cout << setw(9) << "" << setw(15) << "LMAX"     << " = " << LMAX     << endl;
  cout << setw(9) << "" << setw(15) << "MMIN"     << " = " << MMIN     << endl;
  cout << setw(9) << "" << setw(15) << "MMAX"     << " = " << MMAX     << endl;
  cout << setw(9) << "" << setw(15) << "lmax"     << " = " << lmax     << endl;
  cout << setw(9) << "" << setw(15) << "nmax"     << " = " << nmax     << endl;
  cout << setw(9) << "" << setw(15) << "nfreqs"   << " = " << nfreqs   << endl;
  cout << setw(9) << "" << setw(15) << "INCLINE"  << " = " << INCLINE  << endl;
  cout << setw(9) << "" << setw(15) << "PSI"      << " = " << PSI      << endl;
  cout << setw(9) << "" << setw(15) << "PHIP"     << " = " << PHIP     << endl;
  cout << setw(9) << "" << setw(15) << "E"        << " = " << E        << endl;
  cout << setw(9) << "" << setw(15) << "Rperi"    << " = " << Rperi    << endl;
  cout << setw(9) << "" << setw(15) << "Rfac"     << " = " << Rfac     << endl;
  cout << setw(9) << "" << setw(15) << "Mfac"     << " = " << Mfac     << endl;
  cout << setw(9) << "" << setw(15) << "Tmax "    << " = " << Tmax     << endl;
  cout << setw(9) << "" << setw(15) << "delT "    << " = " << delT     << endl;
  cout << setw(9) << "" << setw(15) << "rmin"     << " = " << rmin     << endl;
  cout << setw(9) << "" << setw(15) << "rmax"     << " = " << rmax     << endl;
  cout << setw(9) << "" << setw(15) << "scale"    << " = " << scale    << endl;
  cout << setw(9) << "" << setw(15) << "numr"     << " = " << numr     << endl;
  cout << setw(9) << "" << setw(15) << "nint"     << " = " << nint     << endl;
  cout << setw(9) << "" << setw(15) << "MASS"     << " = " << satmass  << endl;
  cout << setw(9) << "" << setw(15) << "UseCache" << " = " << UseCache << endl;
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

  if (get_value("LMIN", val))		LMIN = atoi(val);
  if (get_value("LMAX", val))		LMAX = atoi(val);
  if (get_value("MMIN", val))		MMIN = atoi(val);
  if (get_value("MMAX", val))		MMIN = atoi(val);
  if (get_value("lmax", val))		lmax = atoi(val);
  if (get_value("nmax", val))		nmax = atoi(val);
  if (get_value("nfreqs", val))		nfreqs = atoi(val);
  if (get_value("numr", val))		numr = atoi(val);
  if (get_value("nint", val))		nint = atoi(val);
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
  if (get_value("Rfac", val))		Rfac = atof(val);
  if (get_value("Mfac", val))		Mfac = atof(val);
  if (get_value("Tmax", val))		Tmax = atof(val);
  if (get_value("delT", val))		delT = atof(val);
  if (get_value("Toffset", val))	Toffset = atof(val);
  if (get_value("MASS", val))		satmass = atof(val);
  if (get_value("INFILE", val))		INFILE = val;
  if (get_value("CACHEDIR", val))	CACHEDIR = val;
  if (get_value("ctrname", val))	ctr_name = val;
  if (get_value("UseCache", val))	UseCache = atoi(val) ? true : false;
}


void UserSatWake::initialize_coefficients()
{
  //==========================
  // Begin normal exectution
  //==========================
  
  const int id  = 0x12;
  KComplex       omp_halo;
  
  //=================
  // Sanity checks
  //=================
  
  if (MMIN<0 && myid==0) cerr << "Warning: Mmin should be >=0 . . . fixing\n";
  MMIN = abs(MMIN);
  if (MMIN>LMIN) {
    if (myid==0) cout << "Setting LMIN=MMIN=" << MMIN << endl;
    LMIN = MMIN;
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
  
  
  double MaxOm = sqrt(halo_model->get_mass(Rfac)/(Rfac*Rfac*Rfac));
  if (myid==0) cout << "Omega(" << Rfac << ")=" << MaxOm  << endl;

  MaxOm *= Mfac*LMAX;
  
  vector< vector<RespMat> > total;

  vector<KComplex> freqs(nfreqs);
  double dOmega = 2.0*MaxOm/(nfreqs-1);
  for (int k=0; k<nfreqs; k++)
    freqs[k] = KComplex(-MaxOm + dOmega*k, OMPI);

  if (myid==0) {
    cout << "Max(Omega)=" << MaxOm  << endl
	 << "    NFreqs=" << nfreqs << endl;
  }

  Times = vector<double>(NUMT+1);
  for (int nt=0; nt<=NUMT; nt++) Times[nt] = -Tmax + 2.0*Tmax*nt/NUMT;
  
  TimeSeriesCoefs Coefs(E, Rperi, delT, Tmax, halo_model, runtag);
  
  map< int, map<int, vector<CMatrix> > > coefs;
  CVector tcoefs;
  int icnt;
  
  // ===================================================================
  // Get combined response matrix for disk and halo for each L
  // ===================================================================
  
  ofstream to_save;
  ifstream from_save;
  unsigned short reading = (UseCache ? 1 : 0);

  if (myid==0 && UseCache) {
    
    from_save.open(string(outdir + cachefile).c_str());

    if (!from_save) {
      cerr << "Couldn't open <" << cachefile <<
	"> to read cached data" << endl;
      cerr << ". . .  will attempt to write a cache file." << endl;
      
      reading = 0;

    } else {
      
      // Read ID string and check version #
      int tid;
      from_save.read((char *)&tid, sizeof(int));
      if (tid != id) {
	cerr << "Incompatible save file!\n";
	MPI_Abort(MPI_COMM_WORLD, -34);
      } 
    
      char tbuf[255];
      from_save.read(tbuf, 255);
      cout << tbuf << endl;

      int nfreq1, lmin1, lmax1, mmin1, mmax1;
      from_save.read((char *)&nfreq1, sizeof(int));
      from_save.read((char *)&lmin1,  sizeof(int));
      from_save.read((char *)&lmax1,  sizeof(int));
      from_save.read((char *)&mmin1,  sizeof(int));
      from_save.read((char *)&mmax1,  sizeof(int));

      if (nfreqs != nfreq1) reading = 0;
      if (LMIN   != lmin1 ) reading = 0;
      if (LMAX   != lmax1 ) reading = 0;
      if (MMIN   != mmin1 ) reading = 0;
      if (MMAX   != mmax1 ) reading = 0;

      if (reading) {
	KComplex tmp;
	for (int i=0; i<nfreqs; i++) {
	  from_save.read((char *)&tmp.real(), sizeof(double));
	  from_save.read((char *)&tmp.imag(), sizeof(double));
	  if (fabs(tmp - freqs[i])>1.0e-10) reading = 0;
	}
      }

      if (reading==0) {
	cerr << "Cache file <" << outdir << cachefile
	     << "> is incompatible with current input parameters!!" << endl;
	MPI_Abort(MPI_COMM_WORLD, 35);
	exit(-1);
      }
    }
  }
    
  MPI_Bcast(&reading, 1, MPI_UNSIGNED_SHORT, 0, MPI_COMM_WORLD);
  
  if (myid && reading) {

    from_save.open(string(outdir + cachefile).c_str());

    // Read ID string and check version #
    int tid;
    from_save.read((char *)&tid, sizeof(int));
    if (tid != id) {
      cerr << "Process " << myid << ": error reading save file!\n";
      MPI_Abort(MPI_COMM_WORLD, -34);
    } 

    char tbuf[255];
    from_save.read(tbuf, 255);

    int nfreq1, lmin1, lmax1, mmin1, mmax1;
    from_save.read((char *)&nfreq1, sizeof(int));
    from_save.read((char *)&lmin1,  sizeof(int));
    from_save.read((char *)&lmax1,  sizeof(int));
    from_save.read((char *)&mmin1,  sizeof(int));
    from_save.read((char *)&mmax1,  sizeof(int));

    KComplex tmp;
    for (int i=0; i<nfreqs; i++) {
      from_save.read((char *)&tmp.real(), sizeof(double));
      from_save.read((char *)&tmp.imag(), sizeof(double));
    }
  }
    
  if (!reading && myid==0 && UseCache) {
      
    to_save.open(string(outdir + cachefile).c_str());
    if (!to_save) {
      cerr << "Couldn't open <" << cachefile <<
	"> to write cached data" << endl;
      MPI_Abort(MPI_COMM_WORLD, -35);
      exit(-1);
    }
    
    // Write ID string and version #
    time_t tp = time(NULL);
    string stime = string("UserSatWake") + ": " + ctime(&tp);
    cout << "ID string: " << stime<< endl;
    
    to_save.write((const char *)&id, sizeof(int));
    char tbuf[255];
    strncpy(tbuf, stime.c_str(), 255);
    to_save.write(tbuf, 255);
    
    to_save.write((const char *)&nfreqs, sizeof(int));
    to_save.write((const char *)&LMIN,   sizeof(int));
    to_save.write((const char *)&LMAX,   sizeof(int));
    to_save.write((const char *)&MMIN,   sizeof(int));
    to_save.write((const char *)&MMAX,   sizeof(int));

    for (int i=0; i<nfreqs; i++) {
      to_save.write((const char *)&freqs[i].real(), sizeof(double));
      to_save.write((const char *)&freqs[i].imag(), sizeof(double));
    }
  }


  if (myid==0) cout << "Computing Laplace coefficients . . ." << endl;
  
  icnt = 0;
  for (int L=LMIN; L<=LMAX; L++) {
    for (int L2=-L;L2<=L; L2+=2) {
      int id = icnt++ % numprocs;
      if (id == myid) {
	Coefs.coefs(L, L2, nmax, nint, u, freqs, Times,  coefs[L][L2]);
      }
    }
  }


  icnt = 0;
  for (int L=LMIN; L<=LMAX; L++) {
    if (myid==0) cout << "L=" << L << endl;
    for (int L2=-L;L2<=L; L2+=2) {
      int id = icnt++ % numprocs;
      unsigned sz = coefs[L][L2].size();
      MPI_Bcast(&sz, 1, MPI_UNSIGNED, id, MPI_COMM_WORLD);
      if (id==myid) {
	vector<CMatrix>::iterator it;
	for (it=coefs[L][L2].begin(); it!=coefs[L][L2].end(); it++)
	  CMatrixSynchronize(*it, id);
      } else {
	CMatrix tmp;
	for (unsigned j=0; j<sz; j++) {
	  CMatrixSynchronize(tmp, id);
	  coefs[L][L2].push_back(tmp);
	}
      }
      
      if (myid==0) cout << "    L2=" << L2 << " rows=[" 
			<< coefs[L][L2][0].getrlow() << ", " 
			<< coefs[L][L2][0].getrhigh() << "]" << endl;
    }
  }

  if (myid==0) cout << "done" << endl;
  
  // Factor for point mass expansion for
  // orbit in x-y plane
  
  INCLINE *= M_PI/180.0;
  PSI     *= M_PI/180.0;
  PHIP    *= M_PI/180.0;
  
  KComplex I(0.0, 1.0);
  
  for (int L=LMIN; L<=LMAX; L++) {
    int mmin;
    if (isOdd(L)) {
      if (isOdd(MMIN)) mmin = MMIN;
      else             mmin = MMIN+1;
    } else {
      if (isOdd(MMIN)) mmin = MMIN+1;
      else             mmin = MMIN;
    }
    for (int M=mmin; M<=L; M+=2) {
      Lhalo.push_back(L);
      Mhalo.push_back(M);
    }
  }

  Nhalo = Lhalo.size();
  if (myid==0) {
    cout << "Harmonics:" << endl << "---------" << endl
	 << setw(3) << left << "L"
	 << setw(3) << left << "M" << endl
	 << setw(3) << left << "-"
	 << setw(3) << left << "-" << endl;
    for (int i=0; i<Nhalo; i++)
      cout << setw(3) << left << Lhalo[i]
	   << setw(3) << left << Mhalo[i]
	   << endl;
  }
  
  icnt = 0;
  for (int nf=0; nf<nfreqs; nf++) {
    
    if (myid==0) cout << "Halo frequency: " << freqs[nf] << endl;
    
    vector<RespMat> mab_halo(Nhalo);

    if (reading) {
      
      for (int ihalo=0; ihalo<Nhalo; ihalo++)
	mab_halo[ihalo] = RespMat(from_save, halo_model, u);
      
      total.push_back(mab_halo);
      
    } else {
      
      //
      // Get response matrices for each halo l-order
      //
      for (int ihalo=0; ihalo<Nhalo; ihalo++) {
	
	mab_halo[ihalo] = RespMat(freqs[nf], Lhalo[ihalo], Mhalo[ihalo], 
				  lmax, nmax, nptsK, nptsE,
				  halo_model, u, OLD, 0, 0,
				  VERBOSE, SITYPE, CAUCHY, RATINT);
	
	mab_halo[ihalo].set_params(DELTA, PTGRID);
	
	int id = icnt++ % numprocs;

	if (id == myid) {
	  if (CAUCHY)
	    mab_halo[ihalo].make_matrix_pv();
	  else 
	    mab_halo[ihalo].make_matrix();
	}
	
      }
      
      total.push_back(mab_halo);
    }

  }

  if (!reading) {

    icnt = 0;

    for (int nf=0; nf<nfreqs; nf++) {

      for (int ihalo=0; ihalo<Nhalo; ihalo++) {

	int id = icnt++ % numprocs;
	total[nf][ihalo].MPISynchronize(id);

	if (UseCache && myid==0)
	  total[nf][ihalo].write_out(to_save);

      }

    }

  }

  
  // ===================================================================
  // For paranoia's sake
  // ===================================================================

  if (UseCache) {
    if (reading)
      from_save.close();
    else
      if (myid==0) to_save.close();
  }  

  // ===================================================================
  // Compute the time series of response vectors
  // ===================================================================
  
  ofstream tlog;
  if (myid==0) tlog.open((outdir + runtag + ".satwake.timelog").c_str());
  
  CVector tmp(1, nmax); tmp.zero();

  //
  // All the time slices for each (L,M) pair
  //

  rcoefs = vector< vector<CVector> >(NUMT+1);
  for (int nt=0; nt<=NUMT; nt++) {
    for (int ihalo=0; ihalo<Nhalo; ihalo++)
      rcoefs[nt].push_back(tmp);
  }
  
  //
  // For interpolation in compute_coefficients()
  //

  curcoefs = vector<CVector>(Nhalo);
  for (int ihalo=0; ihalo<Nhalo; ihalo++)
    curcoefs.push_back(tmp);


  icnt = 0;

  for (int nt=0; nt<=NUMT; nt++) {
    
    if (myid==0) tlog << setw(8) << nt << setw(16) << Times[nt] << endl;
    
    //
    // Frequency loop for current time
    //
    for (int nf=0; nf<nfreqs; nf++) {
      
      //
      // Get satellite perturbation coefficients
      //
      for (int ihalo=0; ihalo<Nhalo; ihalo++) {
	
	int id = icnt++ % numprocs;

	if (id != myid) continue;

	int L = Lhalo[ihalo];
	int M = Mhalo[ihalo];
	
	for (int L2=-L; L2<=L; L2+=2) {
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
	  //                  |                    |
	  //  Rotate component to spherical harmonic in new orientation
	  //
	  
	  
	  //
	  // Leading sign for spherical harmonic with m>0 (following Edmonds
	  // to be consistent with definition of rotation matrices)
	  
	  if (L2>0 && isOdd(L2)) tcoefs *= -1.0;
	}

	tcoefs *= dOmega/2.0*M_PI;
	
	rcoefs[nt][ihalo] 
	  += total[nf][ihalo].get_response(tcoefs, RespMat::self);
      }
    }
  }


  //
  // Share with all processes
  //

  icnt = 0;

  for (int nt=0; nt<=NUMT; nt++) {
    
    for (int ihalo=0; ihalo<Nhalo; ihalo++) {
      int id = icnt++ % numprocs;
      CVectorSynchronize(rcoefs[nt][ihalo], id);
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
  
  PartMapItr it = cC->Particles().begin();

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
	int M = Mhalo[ihalo];

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

  //
  // Find time index
  //
  if (curT<Times.front()) {	// No perturbation if time is out of bounds
    indx = 0;
    a = 0.0;
    b = 0.0;
  } else if (curT>=Times.back()) {
    indx = Times.size()-2;
    a = 0.0;
    b = 0.0;
  } else {
				// Index bounds enforcement
    indx = max<int>(0, min<int>(Vlocate(curT, Times), Times.size()-2));
    a = (Times[indx+1] - curT)/(Times[indx+1] - Times[indx]);
    b = (curT - Times[indx])/(Times[indx+1] - Times[indx]);
  }

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
