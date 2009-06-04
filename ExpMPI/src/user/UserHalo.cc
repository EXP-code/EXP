#include <sys/timeb.h>
#include <math.h>
#include <sstream>

#include "expand.h"

#include <UserHalo.H>


UserHalo::UserHalo(string &line) : ExternalForce(line)
{

  id = "SphericalHalo";		// Halo model file
  model_file = "SLGridSph.model";

  q1 = 1.0;			// Flattening of the Galactic Potentail to X-axis
  q2 = 1.0;			// Flattening of the Galactic Potentail to Y-axis
  q3 = 1.0;			// Flattening of the Galactic Potentail to Z-axis
  diverge = 0;			// Use analytic divergence (true/false)
  diverge_rfac = 1.0;		// Exponent for profile divergence
  ctr_name = "";		// Default component for center

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


  model = new SphericalModelTable(model_file, diverge, diverge_rfac);

  userinfo();
}

UserHalo::~UserHalo()
{
  delete model;
}

void UserHalo::userinfo()
{
  if (myid) return;		// Return if node master node

  print_divider();

  cout << "** User routine SPHERICAL HALO initialized, " ;
  cout << "Filename=" << model_file << "  diverge=" << diverge
       << "  diverge_rfac=" << diverge_rfac;
  if (c0) 
    cout << ", center on component <" << ctr_name << ">";
  else
    cout << ", using inertial center";
  cout << endl;

  cout << "Flattening of halo (q1 , q2 , q3) = (" 
       << q1 << " , " 
       << q2 << " , " 
       << q3 << " ) " << endl; 

  print_divider();
}

void UserHalo::initialize()
{
  string val;

  if (get_value("model_file", val))	model_file = val;
  if (get_value("q1", val))	        q1 = atof(val.c_str());
  if (get_value("q2", val))	        q2 = atof(val.c_str());
  if (get_value("q3", val))	        q3 = atof(val.c_str());
  if (get_value("diverge", val))	diverge = atoi(val.c_str());
  if (get_value("diverge_rfac", val))	diverge_rfac = atof(val.c_str());
  if (get_value("ctrname", val))	ctr_name = val;
}


void UserHalo::determine_acceleration_and_potential(void)
{
  exp_thread_fork(false);

  print_timings("UserHalo: accleration timings");
}


void * UserHalo::determine_acceleration_and_potential_thread(void * arg) 
{
  unsigned nbodies = cC->Number();
  int id = *((int*)arg);
  int nbeg = nbodies*id/nthrds;
  int nend = nbodies*(id+1)/nthrds;

  thread_timing_beg(id);

  double pos[3], rr, r, pot, dpot;
  double qq[3];
  
  qq[0]=q1*q1;
  qq[1]=q2*q2;
  qq[2]=q3*q3;

  PartMapItr it = cC->Particles().begin();

  /*
  if (myid==0) {
    cout << "UserHalo: component=" << cC->name << endl;
    cout << " First: " 
	 << setw(10) << it->first
	 << setw(10) << it->second.indx
	 << setw(18) << scientific << it->second.mass
	 << endl;
    it++;
    cout << "Second: " 
	 << setw(10) << it->first
	 << setw(10) << it->second.indx
	 << setw(18) << scientific << it->second.mass  
	 << endl;
    it++;
    cout << " Third: " 
	 << setw(10) << it->first
	 << setw(10) << it->second.indx
	 << setw(18) << scientific << it->second.mass  
	 << endl;

    struct timeb tp;
    ftime(&tp);

    cout << cC->name << " Number=" << cC->Particles().size() 
	 << " time=" << setw(9) << tp.time << '.' 
	 << setfill('0') << setw(3) << tp.millitm << setfill(' ')
	 << " ptr=" << hex << &(cC->Particles()) << endl << dec;

    it = cC->Particles().begin();
  }
  */

  for (int q=0   ; q<nbeg; q++) it++;
  for (int q=nbeg; q<nend; q++) {
    unsigned long i = (it++)->first;
				// If we are multistepping, compute accel 
				// only at or below this level

    if (multistep && (cC->Part(i)->level < mlevel)) continue;

    rr = 0.0;
    for (int k=0; k<3; k++) {
      pos[k] = cC->Pos(i, k);	// Inertial by default
      if (c0) pos[k] -= c0->center[k];
      rr += pos[k]*pos[k]/qq[k];
    }
    r = sqrt(rr);

    model->get_pot_dpot(r, pot, dpot);

    // DEBUG
    /*
    cout << "#, indx, tnow, mass, r, pot, dpot=" 
	 << setw(4) << myid << setw(8) << i 
	 << setw(10) << fixed << tnow
	 << setw(18) << scientific << cC->Mass(i)
	 << setw(18) << scientific << r 
	 << setw(18) << scientific << pot 
	 << setw(18) << scientific << dpot
	 << fixed << endl;
    */
    // END DEBUG

    // Add external accerlation
    for (int k=0; k<3; k++)
      cC->AddAcc(i, k, -dpot*pos[k]/(qq[k]*r) );

    // Add external potential
    cC->AddPotExt(i, pot );

  }

  thread_timing_end(id);

  return (NULL);
}


extern "C" {
  ExternalForce *makerHalo(string& line)
  {
    return new UserHalo(line);
  }
}

class proxyhalo { 
public:
  proxyhalo()
  {
    factory["userhalo"] = makerHalo;
  }
};

proxyhalo p;
