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

UserResPot::UserResPot(string &line) : ExternalForce(line)
{
  LMAX = 2;
  NMAX = 20;
  NUMR = 1000;
  L0 = 2;
  M0 = 2;
  L1 = -1;
  L2 =  2;
  rmin = 1.0e-4;
  rmax = 1.95;
  scale = 0.067;
  drfac = 0.05;

  ton = -20.0;			// Turn on time
  toff = 1.0e20;		// Turn off time
  delta = 1.0;			// Turn on duration
  toffset = 0.0;		// Time offset for orbit
  omega = 18.9;			// Patern speed

  NUME = 400;			// Points in Energy grid
  NUMK = 100;			// Points in Kappa grid
  NUMI = 2000;			// Points in Action grid

  MASS = 0.05;			// Bar mass
  LENGTH = 0.067;		// Bar length
  COROT = 10;			// Corotation factor
  A21 = 0.2;			// Major to semi-minor ratio
  A32 = 0.05;			// Semi-minor to minor ratio

  domega = 0.0;			// Frequency shift
  t0 = 0.5;			// Mid point of drift
  first = true;

				// Tabled spherical model
  model_file = "SLGridSph.model";
  ctr_name = "";		// Default component for com is none

				// Log file name
  filename = "ResPot." + runtag;

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

  SphereSL::cache = 0;
  SphereSL::mpi = 1;
  SphereSL *sl = new SphereSL(LMAX, NMAX, NUMR, rmin, rmax, scale, hm);
  halo_ortho = sl;

  ResPot::NUME = NUME;
  ResPot::NUMK = NUMK;
  ResPot::NUMI = NUMI;
  respot = new ResPot(halo_model, halo_ortho, L0, M0, L1, L2, NMAX);

  BarForcing::L0 = L0;
  BarForcing::M0 = M0;
  BarForcing bar(NMAX, MASS, LENGTH, COROT);
  bar.compute_quad_parameters(A21, A32);
  bar.compute_perturbation(halo_model, halo_ortho, bcoef, bcoefPP);
  omega = bar.Omega();

  bcount = new int [nthrds];

  userinfo();
}

UserResPot::~UserResPot()
{
  delete halo_model;
  delete halo_ortho;
  delete respot;
  delete [] bcount;
}

void UserResPot::userinfo()
{
  if (myid) return;		// Return if node master node
  print_divider();
  cout << "** User routine RESONANCE POTENTIAL initialized";
  cout << " with Length=" << LENGTH 
       << ", Mass=" << MASS 
       << ", Omega=" << omega 
       << ", Domega=" << domega 
       << ", T0=" << t0
       << ", b/a=" << A21
       << ", c/b=" << A32
       << ", Ton=" << ton
       << ", Toff=" << toff
       << ", Delta=" << delta
       << ", L=" << L0
       << ", M=" << M0
       << ", l_1=" << L1
       << ", l_2=" << L2
       << "\n";
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
  if (get_value("scale", val))     scale = atof(val.c_str());
  if (get_value("drfac", val))    drfac = atof(val.c_str());

  if (get_value("ton", val))      ton = atof(val.c_str());
  if (get_value("toff", val))     toff = atof(val.c_str());
  if (get_value("delta", val))    delta = atof(val.c_str());
  if (get_value("toffset", val))  toffset = atof(val.c_str());

  if (get_value("MASS", val))     MASS = atof(val.c_str());
  if (get_value("LENGTH", val))   LENGTH = atof(val.c_str());
  if (get_value("COROT", val))    COROT = atof(val.c_str());
  if (get_value("A21", val))      A21 = atof(val.c_str());
  if (get_value("A32", val))      A32 = atof(val.c_str());

  if (get_value("domega", val))   domega = atof(val.c_str());
  if (get_value("t0", val))       t0 = atof(val.c_str());

  if (get_value("NUME", val))     NUME = atoi(val.c_str());
  if (get_value("NUMK", val))     NUMK = atoi(val.c_str());
  if (get_value("NUMI", val))     NUMI = atoi(val.c_str());

  if (get_value("model", val))    model_file = val;
  if (get_value("ctrname", val))  ctr_name = val;
  if (get_value("filename", val)) filename = val;
}

void UserResPot::determine_acceleration_and_potential(void)
{
  Omega = omega*(1.0 + domega*(tnow - t0));
  
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
	  istringstream ins(line);

	  ins >> tlast1;
	  ins >> phase1;
	  ins >> omlast1;

	  if (tlast1 >= tpos) {
	    if (firstline) {
	      cerr << "UserResPot: can't read log file, aborting" << endl;
	      MPI_Abort(MPI_COMM_WORLD, 123);
	    }
	    break;
	  }

	  firstline = false;
	  tlast = tlast1;
	  phase = phase1;
	  omlast = omlast1;

	  out << line << "\n";

	  in.getline(line, linesize); // Next line
	}
				// Trapezoidal rule step
	phase += (tnow - tlast)*0.5*(Omega + omlast);
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
	    << setw(15) << "Bounds"
	    << endl;
      }
      
      phase = Omega*tnow;	// Initial phase 
    }

    first = false;

  } else {			// Trapezoidal rule integration
    phase += (tnow - tlast)*0.5*(Omega + omlast);
  }

				// Store current state
  tlast = tnow;
  omlast = Omega;

				// Clear bounds counter
  for (int n=0; n<nthrds; n++) bcount[n] = 0;

  // -----------------------------------------------------------
  // Compute the mapping
  // -----------------------------------------------------------

  exp_thread_fork(false);

  // -----------------------------------------------------------

				// Get total number out of bounds
  int btot=0;
  for (int n=1; n<nthrds; n++) bcount[0] += bcount[n];
  MPI_Reduce(&bcount[0], &btot, 1, MPI_INT, MPI_SUM, 0, MPI_COMM_WORLD);

				// Write diagnostic log
  if (myid==0) {
    ofstream out(filename.c_str(), ios::out | ios::app);
    out << setw(15) << tnow
	<< setw(15) << phase
	<< setw(15) << Omega
	<< setw(15) << btot
	<< endl;
  }

}
void * UserResPot::determine_acceleration_and_potential_thread(void * arg) 
{
  double amp, R2, R;
  double posI[3], posO[3], velI[3], velO[3];
  
  int nbodies = particles->size();
  int id = *((int*)arg);
  int nbeg = nbodies*id/nthrds;
  int nend = nbodies*(id+1)/nthrds;

  amp =
    0.5*(1.0 + erf( (tnow - ton) /delta )) *
    0.5*(1.0 + erf( (toff - tnow)/delta )) ;
    
				// Check for nan (can get rid of this
				// eventually)
  bool found_nan = false;

  for (int i=nbeg; i<nend; i++) {

    R2 = 0.0;
    for (int k=0; k<3; k++) {
      posI[k] = (*particles)[i].pos[k];
      if (c0) posI[k] -= c0->com[k];
      velI[k] = (*particles)[i].vel[k];
      R2 += posI[k]*posI[k];
    }
    R = sqrt(R2);

    bcount[id] += 
      respot->Update(dtime, phase, Omega, amp, bcoef, posI, velI, posO, velO);
		   
    for (int k=0; k<3; k++) {
      (*particles)[i].pos[k] = posO[k];
      (*particles)[i].vel[k] = velO[k];
      (*particles)[i].acc[k] = (velO[k] - velI[k])/dtime;
      if (!found_nan) {
	if ( isnan((*particles)[i].pos[k]) ||
	     isnan((*particles)[i].vel[k]) ||
	     isnan((*particles)[i].acc[k]) ) found_nan = true; 
      }
    }
    (*particles)[i].potext = halo_model->get_pot(R);

    if (found_nan) {
      cout << "Process " << myid << ": found nan\n";
      for (int k=0; k<3; k++) cout << setw(15) << (*particles)[i].pos[k];
      for (int k=0; k<3; k++) cout << setw(15) << (*particles)[i].vel[k];
      for (int k=0; k<3; k++) cout << setw(15) << (*particles)[i].acc[k];
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
