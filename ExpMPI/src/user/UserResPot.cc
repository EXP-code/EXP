#include <math.h>
#include "expand.h"
#include <localmpi.h>

#include <SatelliteOrbit.h>
#include <AxisymmetricBasis.H>
#include <ExternalCollection.H>
#include <ResPot.H>
#include <biorth.h>
#include <sphereSL.h>
#include <UserResPot.H>
#include <BarForcing.H>

#include <sstream>

#include <pthread.h>  
#ifdef DEBUG
static pthread_mutex_t iolock = PTHREAD_MUTEX_INITIALIZER;
#endif

static int respot_mpi_id_var;
int respot_mpi_id()
{
  return respot_mpi_id_var;
}

UserResPot::UserResPot(string &line) : ExternalForce(line)
{
  LMAX = 2;
  NUMR = 800;
  L0 = 2;
  M0 = 2;
  L1 = -1;
  L2 =  2;
  rmin = 1.0e-3;
  rmax = 1.98;
  Klim = 1.0;
  scale = 0.067;
  drfac = 0.05;

  ton = -20.0;			// Turn on time
  toff = 1.0e20;		// Turn off time
  delta = 1.0;			// Turn on duration
  toffset = 0.0;		// Time offset for orbit
  omega = 18.9;			// Patern speed
  phase0 = 0.0;			// Initial phase

  NUMX = 400;			// Points in Ang mom grid
  NUME = 200;			// Points in Energy
  RECS = 100;			// Points in Angle grid
  ITMAX = 50;			// Number of iterations for mapping solution

  MASS = -1.0;			// Bar mass
  MFRAC = 0.05;			// Fraction of enclosed mass
  LENGTH = 0.067;		// Bar length
  AMP = 1.0;			// Mass prefactor
  COROT = 10;			// Corotation factor
  A21 = 0.2;			// Major to semi-minor ratio
  A32 = 0.05;			// Semi-minor to minor ratio

  domega = 0.0;			// Rate of forced slow down
  tom0 = 1.0;			// Midpoint of forced bar slow down
  dtom = -1.0;			// Width of forced bar slow down


  first = true;

  usetag = -1;			// Flag not used unless explicitly defined

				// Tabled spherical model
  model_file = "SLGridSph.model";
  ctr_name = "";		// Default component for com is none

				// Log file name
  filename = outdir + "ResPot." + runtag;

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
    }

  }
  else
    c0 = NULL;

				// Set up for resonance potential
  SphericalModelTable *hm = new SphericalModelTable(model_file);
  halo_model = hm;

				// Perturbation
  if (MASS < 0.0) MASS = hm->get_mass(LENGTH);

  BarForcing::L0 = L0;
  BarForcing::M0 = M0;
  BarForcing bar(NMAX, MFRAC*MASS, LENGTH, COROT);

  bar.set_model(halo_model);
  bar.compute_quad_parameters(A21, A32);
  omega = omega0 = bar.Omega();
  Iz = bar.get_Iz();

  pert = &bar;

  ResPot::NUMX = NUMX;
  ResPot::NUME = NUME;
  ResPot::RECS = RECS;
  ResPot::ITMAX = ITMAX;
  respot = new ResPot(halo_model, pert, L0, M0, L1, L2);

  btotn = vector<int>(8);
  bcount = vector< vector<int> >(nthrds);
  difLz = vector<double>(nthrds);
  for (int i=0; i<nthrds; i++)
    bcount[i] = vector<int>(8);
  
  userinfo();

  respot_mpi_id_var = myid;
}

UserResPot::~UserResPot()
{
  delete halo_model;
  delete respot;
}

void UserResPot::userinfo()
{
  if (myid) return;		// Return if node master node
  print_divider();
  cout << "** User routine RESONANCE POTENTIAL initialized";
  cout << " with Length=" << LENGTH 
       << ", Mass=" << MASS 
       << ", Mfrac=" << MFRAC 
       << ", Amp=" << AMP
       << ", Iz=" << Iz
       << ", Omega=" << omega 
       << ", Domega=" << domega 
       << ", tom0=" << tom0
       << ", dtom=" << dtom
       << ", b/a=" << A21
       << ", c/b=" << A32
       << ", Ton=" << ton
       << ", Toff=" << toff
       << ", Delta=" << delta
       << ", L=" << L0
       << ", M=" << M0
       << ", L1=" << L1
       << ", L2=" << L2
       << ", Klim=" << Klim;
  if (dtom>0) cout << ", T_om=" << tom0 << ", dT_om=" << dtom;
  cout << ", Domega=" << domega;
  if (usetag>=0)
    cout << ", with bad value tagging";
  cout << ", ITMAX=" << ITMAX << "\n";
  print_divider();
}

void UserResPot::initialize()
{
  string val;

  if (get_value("LMAX", val))     LMAX = atoi(val.c_str());
  if (get_value("NMAX", val))     NMAX = atoi(val.c_str());
  if (get_value("NUMR", val))     NUMR = atoi(val.c_str());

  if (get_value("L0", val))       L0 = atoi(val.c_str());
  if (get_value("M0", val))       M0 = atoi(val.c_str());
  if (get_value("L1", val))       L1 = atoi(val.c_str());
  if (get_value("L2", val))       L2 = atoi(val.c_str());

  if (get_value("rmin", val))     rmin = atof(val.c_str());
  if (get_value("rmax", val))     rmax = atof(val.c_str());
  if (get_value("Klim", val))     Klim = atof(val.c_str());
  if (get_value("scale", val))    scale = atof(val.c_str());
  if (get_value("drfac", val))    drfac = atof(val.c_str());

  if (get_value("ton", val))      ton = atof(val.c_str());
  if (get_value("toff", val))     toff = atof(val.c_str());
  if (get_value("delta", val))    delta = atof(val.c_str());
  if (get_value("toffset", val))  toffset = atof(val.c_str());
  if (get_value("phase0", val))   phase0 = atof(val.c_str());

  if (get_value("MASS", val))     MASS = atof(val.c_str());
  if (get_value("MFRAC", val))    MFRAC = atof(val.c_str());
  if (get_value("LENGTH", val))   LENGTH = atof(val.c_str());
  if (get_value("AMP", val))      AMP = atof(val.c_str());
  if (get_value("COROT", val))    COROT = atof(val.c_str());
  if (get_value("A21", val))      A21 = atof(val.c_str());
  if (get_value("A32", val))      A32 = atof(val.c_str());

  if (get_value("NUMX", val))     NUMX = atoi(val.c_str());
  if (get_value("NUME", val))     NUME = atoi(val.c_str());
  if (get_value("RECS", val))     RECS = atoi(val.c_str());
  if (get_value("ITMAX", val))    ITMAX = atoi(val.c_str());
  
  if (get_value("domega", val))   domega = atof(val.c_str());
  if (get_value("tom0", val))     tom0 = atof(val.c_str());
  if (get_value("dtom", val))     dtom = atof(val.c_str());


  if (get_value("model", val))    model_file = val;
  if (get_value("ctrname", val))  ctr_name = val;
  if (get_value("filename", val)) filename = val;
  if (get_value("usetag", val))   usetag = atoi(val.c_str());
}

void UserResPot::determine_acceleration_and_potential(void)
{

  if (first) {

    if (restart) {

      if (myid == 0) {
				// Backup up old file
	string backupfile = filename + ".bak";
	string command("cp ");
	command += filename + " " + backupfile;
	system(command.c_str());
	
				// Open new output stream for writing
	ofstream out(filename.c_str());
	if (!out) {
	  cout << "UserResPot: error opening new log file <" 
	       << filename << "> for writing\n";
	  MPI_Abort(MPI_COMM_WORLD, 121);
	  exit(0);
	}
	
				// Open old file for reading
	ifstream in(backupfile.c_str());
	if (!in) {
	  cout << "UserResPot: error opening original log file <" 
	       << backupfile << "> for reading\n";
	  MPI_Abort(MPI_COMM_WORLD, 122);
	  exit(0);
	}

	const int linesize = 1024;
	char line[linesize];
	
	in.getline(line, linesize); // Discard header
	in.getline(line, linesize); // Next line

	double phase1, omlast1, tlast1;
	bool firstline = true;

	while (in) {

	  if (line[0] != '#' && line[0] != '!') {

	    istringstream ins(line);

	    ins >> tlast1;
	    ins >> phase1;
	    ins >> omlast1;

	    if (tlast1 >= tpos) {
	      if (firstline) {
		cerr << "UserResPotN: can't read log file, aborting" << endl;
		cerr << "UserResPotN: line=" << line << endl;
		MPI_Abort(MPI_COMM_WORLD, 123);
	      }
	      break;
	    }

	    firstline = false;
	    tlast = tlast1;
	    phase = phase1;
	    omlast = omlast1;

	  }

	  out << line << "\n";
	  
	  in.getline(line, linesize); // Next line
	}
				// Trapezoidal rule step
	phase += (tnow - tlast)*0.5*(omega + omlast);
      }

      MPI_Bcast(&tlast, 1, MPI_DOUBLE, 0, MPI_COMM_WORLD);
      MPI_Bcast(&phase, 1, MPI_DOUBLE, 0, MPI_COMM_WORLD);
      MPI_Bcast(&omega, 1, MPI_DOUBLE, 0, MPI_COMM_WORLD);
    }

    if (!restart) {

      if (myid==0) {		// Write header
	ofstream out(filename.c_str(), ios::out | ios::app);
	out.setf(ios::left);
	out << setw(15) << "# Time"
	    << setw(15) << "Phase"
	    << setw(15) << "Omega"
	    << setw(15) << "dOmega(tot)"
	    << setw(15) << "Bounds"
	    << endl;
	for (int j=1; j<=8; j++)
	  out << setw(15) << ResPot::ReturnDesc[j];
	out << endl;
	
	char c = out.fill('-');
	int ncnt=1;
	out << "# " << setw(13) << ncnt++ 
	    << "| " << setw(13) << ncnt++
	    << "| " << setw(13) << ncnt++
	    << "| " << setw(13) << ncnt++
	    << "| " << setw(13) << ncnt++;
      	for (int j=1; j<=8; j++)
	  out << "| " << setw(13) << ncnt++;
	out << endl;
	out.fill(c);
      }
      
      phase = phase0;	// Initial phase 
    }

  } else {			// Trapezoidal rule integration
    phase += (tnow - tlast)*0.5*(omega + omlast);
  }

				// Store current state
  tlast = tnow;
  omlast = omega;

				// Clear bounds counter
  for (int n=0; n<nthrds; n++) {
    for (int j=0; j<8; j++) bcount[n][j] = 0;
  }

				// Clear difLz array
  for (int n=0; n<nthrds; n++) difLz[n] = 0.0;
  difLz0 = 0.0;

  // -----------------------------------------------------------
  // Compute the mapping
  // -----------------------------------------------------------

  exp_thread_fork(false);

  // -----------------------------------------------------------

				// Get total number out of bounds
  for (int n=1; n<nthrds; n++) {
    for (int j=0; j<8; j++) {
      bcount[0][j] += bcount[n][j];
      btotn[j] = 0;
    }
  }
  MPI_Reduce(&(bcount[0][0]), &btotn[0], 8, MPI_INT, MPI_SUM, 0, MPI_COMM_WORLD);

				// Get total change in angular momentum
  difLz0 = 0.0;
  for (int n=1; n<nthrds; n++) difLz[0] += difLz[n];
  MPI_Allreduce(&difLz[0], &difLz0, 1,
		MPI_DOUBLE, MPI_SUM, MPI_COMM_WORLD);

  if (dtom>0.0)
    omega = omega0*(1.0 + domega*0.5*(1.0 + erf( (tnow - tom0)/dtom )));
  else
    omega = omega0*(1.0 + domega*(tnow - tom0*0.5));

				// Write diagnostic log
  if (myid==0) {
    int btot=0;
    for (int j=0; j<8; j++) btot += btotn[j];
    ofstream out(filename.c_str(), ios::out | ios::app);
    out.setf(ios::left);
    out << setw(15) << tnow
	<< setw(15) << phase
	<< setw(15) << omega
	<< setw(15) << -difLz0/Iz
	<< setw(15) << btot;
    for (int j=0; j<8; j++)
      out << setw(15) << btotn[j];
    out << endl;
  }

  first = false;
}
void * UserResPot::determine_acceleration_and_potential_thread(void * arg) 
{
  double amp, R2, R;
  double posI[3], posO[3], velI[3], velO[3], Lz0, Lz1;
  
  unsigned nbodies = cC->Number();
  int id = *((int*)arg);
  int nbeg = nbodies*id/nthrds;
  int nend = nbodies*(id+1)/nthrds;

  amp = AMP *
    0.5*(1.0 + erf( (tnow - ton) /delta )) *
    0.5*(1.0 + erf( (toff - tnow)/delta )) ;
    
  
  vector<double> Phase(3);
  Phase[0] = phase;
  Phase[1] = phase + omega*0.5*dtime;
  Phase[2] = phase + omega*dtime;

				// Check for nan (can get rid of this
				// eventually)
  bool found_nan = false;
  ResPot::ReturnCode ret = ResPot::OK;
  double dpot;

  for (int i=nbeg; i<nend; i++) {

    if (usetag>=0 && cC->Part(i)->iattrib[usetag]) continue;

				// Initial conditions
    for (int k=0; k<3; k++) {
      posI[k] = cC->Pos(i, k);
      if (c0) posI[k] -= c0->com[k];
      velI[k] = cC->Vel(i, k);
    }
    
				// Initial Lz
    Lz0 = posI[0]*velI[1] - posI[1]*velI[0];
    
    
    R2 = 0.0;
    for (int k=0; k<3; k++) {
      R2 += posI[k]*posI[k];
    }
    R = sqrt(R2);

    if (R>rmin && R<rmax && ret == ResPot::OK) {

      if ((ret=respot-> 
	   Update(dtime, Phase, amp, posI, velI, posO, velO)) == ResPot::OK) {
	
				// Current ang mom
	Lz1 = posO[0]*velO[1] - posO[1]*velO[0];

				// Accumulate change in Lz for each resonance
	if (respot->K()<Klim)
	  difLz[id] += cC->Mass(i)*(Lz1 - Lz0);
	
      }
	
    } else {
      if (ret != ResPot::OK) bcount[id][ret-1]++;
      if (usetag>=0) cC->Part(i)->iattrib[usetag] = 1;
#ifdef DEBUG
      pthread_mutex_lock(&iolock);
      cout << "Process " << myid << " id=" << id << ":"
	   << " i=" << myid << " Error=" << ResPot::ReturnDesc[ret] << endl;
      pthread_mutex_unlock(&iolock);
#endif

      dpot = halo_model->get_dpot(R);

      for (int k=0; k<3; k++) {
	posO[k] = velI[k] * dtime;
	if (R>0.01*rmin) velO[k] = -dpot*posI[k]/R * dtime;
      }
    }

    for (int k=0; k<3; k++) {
      cC->Part(i)->pos[k] = posO[k];
      cC->Part(i)->vel[k] = velO[k];
      cC->Part(i)->acc[k] = (velO[k] - velI[k])/dtime;
      if (!found_nan) {
	if ( isnan(cC->Pos(i, k)) ||
	     isnan(cC->Vel(i, k)) ||
	     isnan(cC->Acc(i, k)) ) found_nan = true; 
      }
    }
    cC->Part(i)->potext = halo_model->get_pot(R);
    
    if (found_nan) {
      cout << "Process " << myid << ": found nan\n";
      for (int k=0; k<3; k++) cout << setw(15) << cC->Pos(i, k);
      for (int k=0; k<3; k++) cout << setw(15) << cC->Vel(i, k);
      for (int k=0; k<3; k++) cout << setw(15) << cC->Acc(i, k);
      cout << endl << flush;
      found_nan = false;
    }
    
  }

  return (NULL);
}


extern "C" {
  ExternalForce *makerResPot(string& line)
  {
    return new UserResPot(line);
  }
}

class proxy { 
public:
  proxy()
  {
    factory["userrespot"] = makerResPot;
  }
};

static proxy p;
