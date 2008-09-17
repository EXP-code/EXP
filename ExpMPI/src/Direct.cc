#include "expand.h"
#include <localmpi.h>

#include <Direct.H>

#ifdef DEBUG
static pthread_mutex_t iolock = PTHREAD_MUTEX_INITIALIZER;
#endif
static const int MSGTAG=103;

Direct::Direct(string& line) : PotAccel(line)
{
  soft_indx = 0;
  soft = 0.01;
  fixed_soft = true;
  ndim = 4;

  // Parameters for extended pm model
  pm_model = false;
  pmmodel_file = "SLGridSph.model";   // Table od the Extended pm model
  diverge = 0;                        // Use analytic divergence (true/false)
  diverge_rfac = 1.0;                 // Exponent for profile divergence

  initialize();

  if (pm_model) pmmodel = new SphericalModelTable(pmmodel_file, diverge, diverge_rfac);

				// Assign the ring topology
  to_proc = (myid+1) % numprocs;
  from_proc = (myid+numprocs-1) % numprocs;

				// Buffer pointers
  tmp_buffer = NULL;
  bod_buffer = NULL;
}

Direct::~Direct()
{
  delete [] tmp_buffer;
  delete [] bod_buffer;
}

void Direct::initialize(void)
{
  string val;

  if (get_value("soft_indx", val)) {
    soft_indx = atoi(val.c_str());
    fixed_soft = false;
    ndim = 5;
  }

  if (get_value("soft", val)) {
    soft = atof(val.c_str());
    fixed_soft = true;
    ndim = 4;
  }

  if (get_value("pm_model",val))          pm_model = atoi(val.c_str()) ? true : false;
  if (get_value("diverge", val))          diverge = atoi(val.c_str());
  if (get_value("diverge_rfac", val))     diverge_rfac = atof(val.c_str());
  if (get_value("pmmodel_file", val))     pmmodel_file = val;

}

void Direct::get_acceleration_and_potential(Component* C)
{
  cC = C;
  nbodies = cC->Number();

  /*======================================*/
  /* Determine potential and acceleration */
  /*======================================*/

  determine_acceleration_and_potential();
}

void Direct::determine_acceleration_and_potential(void)
{
				// Make sure softening is defined if needed
  if (!fixed_soft && component->ndattrib<soft_indx+1) {
    if (myid==0) cerr << "Direct: particle softening data missing\n";
    MPI_Abort(MPI_COMM_WORLD, 103);
    exit(0);
  }
				// Determine size of largest nbody list
  ninteract = component->Number();
  MPI_Allreduce(&ninteract, &max_bodies, 1, MPI_INT, MPI_MAX, MPI_COMM_WORLD);
  MPI_Allreduce(&ninteract, &used,       1, MPI_INT, MPI_SUM, MPI_COMM_WORLD);

  if (max_bodies==0) return;

#ifdef DEBUG
  cout << "Process " << myid 
       << ": max bodies=" << max_bodies
       << "  direct ninteract=" << ninteract
       << "  name=" << cC->name;
  if (use_external) 
    cout << " (external bodies)" << endl;
  else 
    cout << " (local bodies)" << endl;
  MPI_Barrier(MPI_COMM_WORLD);
#endif
  
				// Allocate buffers to handle largest list
  delete [] tmp_buffer;
  delete [] bod_buffer;

  int buffer_size = max_bodies*ndim;
  tmp_buffer = new double [buffer_size];
  bod_buffer = new double [buffer_size];
  
  // Load body buffer with local interactors
  double *p = bod_buffer;
  unsigned long i;
  map<unsigned long, Particle>::iterator it = component->Particles().begin();

  for (int q=0; q<ninteract; q++) {
    i = (it++)->first;
    *(p++) = component->Mass(i);
    *(p++) = component->Pos(i, 0);
    *(p++) = component->Pos(i, 1);
    *(p++) = component->Pos(i, 2);
    if (!fixed_soft) *(p++) = component->Part(i)->dattrib[soft_indx];
  }

				// Do the local interactors
  exp_thread_fork(false);

				// Do the ring . . . 
  for(int n=1; n<numprocs; n++) {
      
    MPI_Request req1, req2;
    MPI_Status stat;
    
				// Copy current to temp buffer
    memcpy(tmp_buffer, bod_buffer, buffer_size*sizeof(double));

				// Get NEW buffer from right
    MPI_Irecv(bod_buffer, buffer_size, MPI_DOUBLE, from_proc, MSGTAG, 
	      MPI_COMM_WORLD, &req1);

				// Send OLD buffer to left
    MPI_Isend(tmp_buffer, ninteract*ndim, MPI_DOUBLE, to_proc, MSGTAG, 
	      MPI_COMM_WORLD, &req2);

    MPI_Wait(&req2, &stat);
    MPI_Wait(&req1, &stat);

				// How many particles did we get?
    MPI_Get_count(&stat, MPI_DOUBLE, &ninteract);
    ninteract /= ndim;
	
				// Accumulate the interactions
    exp_thread_fork(false);

  }
				// Clear external potential flag
  use_external = false;
}

void * Direct::determine_acceleration_and_potential_thread(void * arg)
{
  double rr, rr0, rfac;
  double mass, pos[3], eps = soft;
  double *p;

  unsigned nbodies = cC->levlist[mlevel].size();
  int id = *((int*)arg);
  int nbeg = nbodies*id/nthrds;
  int nend = nbodies*(id+1)/nthrds;

  double adb = component->Adiabatic();

#ifdef DEBUG
  double tclausius[nthrds];
  for (int i=0; i<nthrds; i++) tclausius[id] = 0.0;
  unsigned ncnt=0;
#endif

  unsigned long j;		// Index of the current local particle

  for (int i=nbeg; i<nend; i++) {
    
    j = cC->levlist[mlevel][i];

				// Don't need acceleration for frozen particles
    if (cC->freeze(j)) continue;
    
				// Loop through the particle list
    p = bod_buffer;
    for (int n=0; n<ninteract; n++) {
				// Get current interaction particle
      mass = *(p++) * adb;
      for (int k=0; k<3; k++) pos[k] = *(p++);
      if (!fixed_soft) eps = *(p++);

				// Compute interparticle squared distance
      rr0 = 0.0;
      for (int k=0; k<3; k++)
	rr0 += 
	  (cC->Pos(j, k) - pos[k]) *
	  (cC->Pos(j, k) - pos[k]) ;
      
				// Compute softened distance
      rr = sqrt(rr0+eps*eps);

                                // Extended model for point masses
                                // Given model provides normalized mass distrbution
      if(pm_model && pmmodel->get_max_radius() > rr) 
	{
	  double mass_frac;
	  mass_frac = pmmodel->get_mass(rr) / pmmodel->get_mass(pmmodel->get_max_radius());
	  mass *= mass_frac;
	}
				// Acceleration
      rfac = 1.0/(rr*rr*rr);
	
      for (int k=0; k<3; k++)
	cC->AddAcc(j, k, -mass *(cC->Pos(j, k) - pos[k]) * rfac );
      
				// Potential
      if (use_external) {
	cC->AddPotExt(j, -mass/rr );
#ifdef DEBUG
	ncnt++;
	for (int k=0; k<3; k++)
	  tclausius[id] += -mass *
	    (cC->Pos(j, k) - pos[k]) * cC->Pos(j, k) * rfac;
#endif
      }
      else if (rr0 > 1.0e-16)	// Ignore "self" potential
	cC->AddPot(j, -mass/rr );
    }
  }
  
#ifdef DEBUG
  if (use_external) {
    pthread_mutex_lock(&iolock);
    cout << "Process " << myid << ", id=" << id << ": ninteract=" << ninteract
	 << "  nexterna=" << ncnt << "  VC=" << tclausius[id] << endl;
    pthread_mutex_unlock(&iolock);
  }
#endif
}

void Direct::determine_coefficients(void) {}
void * Direct::determine_coefficients_thread(void *arg) {}

