#include <math.h>
#include "expand.h"
#include <localmpi.h>

#include <SatelliteOrbit.h>
#include <AxisymmetricBasis.H>
#include <ExternalCollection.H>
#include <ResPot.H>
#include <biorth.h>
#include <sphereSL.h>
#include <UserRPtest.H>

#include <sstream>


UserRPtest::UserRPtest(string &line) : ExternalForce(line)
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

  NUME = 400;			// Points in Energy grid
  NUMK = 100;			// Points in Kappa grid

  with_ps = false;		// Don't print phase space for each particle
  npart = 5;			// Number of particles to trace

  first = true;

				// Tabled spherical model
  model_file = "SLGridSph.model";
  ctr_name = "";		// Default component for com is none

				// Log file name
  filename = "RPtest." + runtag;

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
  respot = new ResPot(halo_model, halo_ortho, L0, M0, L1, L2, NMAX);


  bcoef.setsize(1, NMAX);
  bcoef.zero();

  userinfo();
}

UserRPtest::~UserRPtest()
{
  delete halo_model;
  delete halo_ortho;
  delete respot;
}

void UserRPtest::userinfo()
{
  if (myid) return;		// Return if node master node
  print_divider();
  cout << "** User routine RESONANCE POTENTIAL TEST initialized";
  cout << ", L=" << L0
       << ", M=" << M0
       << ", l_1=" << L1
       << ", l_2=" << L2
       << "\n";
  print_divider();
}

void UserRPtest::initialize()
{
  string val;

  if (get_value("L0", val))       L0 = atoi(val.c_str());
  if (get_value("M0", val))       M0 = atoi(val.c_str());
  if (get_value("L1", val))       L1 = atoi(val.c_str());
  if (get_value("L2", val))       L2 = atoi(val.c_str());

  if (get_value("rmin", val))     rmin = atof(val.c_str());
  if (get_value("rmax", val))     rmax = atof(val.c_str());
  if (get_value("scale", val))     scale = atof(val.c_str());

  if (get_value("NUME", val))     NUME = atoi(val.c_str());
  if (get_value("NUMK", val))     NUMK = atoi(val.c_str());

  if (get_value("with_ps", val))  with_ps = atoi(val.c_str()) ? true : false;
  if (get_value("npart", val))    npart = atoi(val.c_str());

  if (get_value("model", val))    model_file = val;
  if (get_value("ctrname", val))  ctr_name = val;
  if (get_value("filename", val)) filename = val;
}

void UserRPtest::determine_acceleration_and_potential(void)
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
	  cout << "UserRPtest: error opening new log file <" 
	       << filename << "> for writing\n";
	  MPI_Abort(MPI_COMM_WORLD, 121);
	  exit(0);
	}
	
				// Open old file for reading
	ifstream in(backupfile.c_str());
	if (!in) {
	  cout << "UserRPtest: error opening original log file <" 
	       << backupfile << "> for reading\n";
	  MPI_Abort(MPI_COMM_WORLD, 122);
	  exit(0);
	}

	const int linesize = 1024;
	char line[linesize];
	
	in.getline(line, linesize); // Discard header
	in.getline(line, linesize); // Next line

	double tlast1;
	bool firstline = true;

	while (in) {
	  istringstream ins(line);

	  ins >> tlast1;

	  if (tlast1 >= tpos) {
	    if (firstline) {
	      cerr << "UserRPtest: can't read log file, aborting" << endl;
	      MPI_Abort(MPI_COMM_WORLD, 123);
	    }
	    break;
	  }

	  firstline = false;

	  out << line << "\n";

	  in.getline(line, linesize); // Next line
	}
      }

    }

    if (!restart) {

      if (myid == 0) {		// Write header
	ofstream out(filename.c_str(), ios::out | ios::app);
	out.setf(ios::left);
	out << setw(15) << "# Time";

	const int nlabels = 8;
	string labels[] = {"E", "K", "w1", "w2", "w3", "f", "beta", "psi"};

	const int nlabels2 = 6;
	string labels2[] = {"X", "Y", "Z", "U", "V", "W"};

	for (int i=1; i<=npart; i++) {
	  for (int j=0; j<nlabels; j++) {
	    ostringstream sout;
	    sout << labels[j] << "(" << i << ")";
	    out << "| " << setw(13) << sout.str();
	  }
	  if (with_ps) {
	    for (int j=0; j<nlabels2; j++) {
	      ostringstream sout;
	      sout << labels2[j] << "(" << i << ")";
	      out << "| " << setw(13) << sout.str();
	    }
	  }
	}
	out << endl;


	
	out << setw(15) << "# 1";
	int cntr = 2;

	for (int i=1; i<=npart; i++) {
	  for (int j=0; j<nlabels; j++) {
	    out << "| " << setw(13) << cntr++;
	  }
	  if (with_ps) {
	    for (int j=0; j<nlabels2; j++) {
	    out << "| " << setw(13) << cntr++;
	    }
	  }
	}
	cout << endl;
	
      }
      
    }

    first = false;

  }

  exp_thread_fork(false);
}

void * UserRPtest::determine_acceleration_and_potential_thread(void * arg) 
{
  double pos[3], vel[3];
  double R2, R, pot, dpot;
  double E, K, J, w1, w2, w3, f, beta, psi;
  
  int nbodies = particles->size();
  int id = *((int*)arg);
  int nbeg = nbodies*id/nthrds;
  int nend = nbodies*(id+1)/nthrds;

  ofstream out;

  if (myid==0 && id==0) {
    out.open(filename.c_str(), ios::out | ios::app);
    out.setf(ios::left);
  }


  for (int i=nbeg; i<nend; i++) {

    R2 = 0.0;
    for (int k=0; k<3; k++) {
      pos[k] = (*particles)[i].pos[k];
      if (c0) pos[k] -= c0->com[k];
      vel[k] = (*particles)[i].vel[k];
      R2 += pos[k]*pos[k];
    }
    R = sqrt(R2);

    halo_model->get_pot_dpot(R, pot, dpot);

    if (myid==0 && id==0 && i<npart) {
      if (i==0) out << setw(15) << tpos;
      respot->coord(pos, vel, E, K, J, w1, w2, w3, f, beta, psi);
      out << setw(15) << E
	  << setw(15) << K
	  << setw(15) << w1
	  << setw(15) << w2
	  << setw(15) << w3
	  << setw(15) << f
	  << setw(15) << beta
	  << setw(15) << psi;

      if (with_ps) {
	for (int k=0; k<3; k++) out << setw(15) << pos[k];
	for (int k=0; k<3; k++) out << setw(15) << vel[k];
      }
      if (i==npart-1) out << endl;
    }

    for (int k=0; k<3; k++) (*particles)[i].acc[k] += -dpot*pos[k]/R;
    
    (*particles)[i].potext += (*particles)[i].mass * pot;

  }

  if (myid==0 && id==0) {
    out.close();
  }

  return (NULL);
}


extern "C" {
  ExternalForce *makerRPtest(string& line)
  {
    return new UserRPtest(line);
  }
}

class proxyP { 
public:
  proxyP()
  {
    factory["userrptest"] = makerRPtest;
   }
};

static proxyP p;
