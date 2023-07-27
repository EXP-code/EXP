#include "expand.H"

#include <Direct.H>

#ifdef DEBUG
static pthread_mutex_t iolock = PTHREAD_MUTEX_INITIALIZER;
#endif
static const int MSGTAG=103;

const std::set<std::string>
Direct::valid_keys = {
  "soft_indx",
  "soft",
  "type",
  "mn_model",
  "a",
  "b",
  "pm_model",
  "diverge",
  "diverge_rfac",
  "pmmodel_file"
};

Direct::Direct(Component* c0, const YAML::Node& conf) : PotAccel(c0, conf)
{
  // Standard softening
  //
  soft_indx  = 0;
  soft       = 0.01;
  fixed_soft = true;
  ndim       = 4;

  // Miyamoto-Nagai-shape particles
  //
  mn_model = false;		// Turn on MN disk (off by default)
  a        = 0.01;		// Disk scale length
  b        = 0.002;		// Disk scale height
  
  // Parameters for extended pm model
  //
  pm_model = false;		    // This invokes a spherical model profile
  pmmodel_file = "SLGridSph.model"; // Table od the Extended pm model
  diverge      = 0;	            // Use analytic divergence (true/false)
  diverge_rfac = 1.0;               // Exponent for profile divergence

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
  // Remove matched keys
  //
  for (auto v : valid_keys) current_keys.erase(v);
  
  // Assign values from YAML
  //
  try {
    if (conf["soft_indx"]) {
      soft_indx = conf["soft_indx"].as<int>();
      fixed_soft = false;
      ndim = 5;
    }
    
    if (conf["soft"]) {
      soft = conf["soft"].as<double>();
      fixed_soft = true;
      ndim = 4;
    }

    if (conf["type"]) {
      std::string type = conf["type"].as<std::string>();
      if (type.compare("Spline") == 0) kernel = std::make_shared<SplineSoft>();
      else                             kernel = std::make_shared<PlummerSoft>();
    } else {
      kernel = std::make_shared<SplineSoft>();
      if (myid==0) std::cout << "Direct: using SplineSoft" << std::endl;
    }

    if (conf["mn_model"])         mn_model     = conf["mn_model"].as<bool>();
    if (conf["a"])                a            = conf["a"].as<double>();
    if (conf["b"])                b            = conf["b"].as<double>();

    if (conf["pm_model"])         pm_model     = conf["pm_model"].as<bool>();
    if (conf["diverge"])          diverge      = conf["diverge"].as<int>();
    if (conf["diverge_rfac"])     diverge_rfac = conf["diverge_rfac"].as<double>();
    if (conf["pmmodel_file"])     pmmodel_file = conf["pmmodel_file"].as<std::string>();
  }
  catch (YAML::Exception & error) {
    if (myid==0) std::cout << "Error parsing parameters in Direct: "
			   << error.what() << std::endl
			   << std::string(60, '-') << std::endl
			   << "Config node"        << std::endl
			   << std::string(60, '-') << std::endl
			   << conf                 << std::endl
			   << std::string(60, '-') << std::endl;
    MPI_Finalize();
    exit(-1);
  }

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
    std::string msg("Direct: particle softening data missing");
    throw GenericError(msg, __FILE__, __LINE__, 1019, false);
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
  PartMapItr it = component->Particles().begin();

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
  for (int n=1; n<numprocs; n++) {
      
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
  double pos[3], eps = soft;

  int id = *((int*)arg);

#ifdef DEBUG
  double tclausius[nthrds];
  for (int i=0; i<nthrds; i++) tclausius[id] = 0.0;
  unsigned ncnt=0;
#endif

  // If we are multistepping, compute accel only at or above <mlevel>
  //
  for (int lev=mlevel; lev<=multistep; lev++) {

    unsigned nbodies = cC->levlist[lev].size();
    int nbeg = nbodies*id/nthrds;
    int nend = nbodies*(id+1)/nthrds;

    double adb = component->Adiabatic();

    for (int i=nbeg; i<nend; i++) {
    
				// Index of the current local particle
      unsigned long j = cC->levlist[lev][i];

				// Don't need acceleration for frozen particles
      if (cC->freeze(j)) continue;
    
				// Loop through the particle list
      double * p = bod_buffer;
      for (int n=0; n<ninteract; n++) {

				// Get current interaction particle
				// mass from ring buffer
	double mass = *(p++) * adb;
				// Position
	for (int k=0; k<3; k++) pos[k] = *(p++);
	
	// Compute interparticle squared distance
	//
	double rr0 = 0.0;
	for (int k=0; k<3; k++) rr0 +=
				  (cC->Pos(j, k) - pos[k]) *
				  (cC->Pos(j, k) - pos[k]) ;
	double rr = sqrt(rr0);
	
	// Reject particle at current location
	//
	if (rr>rtol) {
	  
	  // BEG: Miyamoto-Nagai (MN) disk-shaped point mass
	  if (mn_model) {	

	    // Positions relative to point mass center
	    //
	    double xx = cC->Pos(j, 0) - pos[0];
	    double yy = cC->Pos(j, 1) - pos[1];
	    double zz = cC->Pos(j, 2) - pos[2];

	    // MN intermediate computation
	    //
	    double rr = sqrt( xx*xx + yy*yy );
	    double zb = sqrt( zz*zz + b * b );
	    double ab = a + zb;
	    double dn = sqrt( rr*rr + ab*ab );
	    
	    // MN potential and cylindrical force
	    //
	    double pot = -mass/dn;
	    double fr  = -mass*rr/(dn*dn*dn);
	    double fz  = -mass*zz*ab/(zb*dn*dn*dn);
	  
	    // Add acceleration by disk particle
	    //
	    cC->AddAcc(j, 0, fr*xx/(rr+1.0e-10) );
	    cC->AddAcc(j, 1, fr*yy/(rr+1.0e-10) );
	    cC->AddAcc(j, 2, fz );
	  
#ifdef DEBUG
	    if (use_external) {
	      ncnt++;
	      tclausius[id] += mass * ( (fr*xx + fr*yy)/(rr+1.0e-10) + fz);
	    }
#endif
	    // Particle potential
	    cC->AddPot(j, pot );

	  }
	  // END: Miyamoto-Nagai point mass
	  // BEG: Spherical point mass
	  else {
	  
	    if (!fixed_soft) eps = *(p++);

				// Extended model for point masses
                                // Given model provides normalized mass distrbution
	    double pot = 0.0;
	  
	    if (pm_model && pmmodel->get_max_radius() > rr) {
	      double mass_frac =
		pmmodel->get_mass(rr) / pmmodel->get_mass(pmmodel->get_max_radius());
	      pot = pmmodel->get_pot(rr)/pmmodel->get_mass(pmmodel->get_max_radius());
	      mass *= mass_frac;
	    } else {
	      auto y = (*kernel)(rr, eps);
	      pot = mass * y.second;
	      mass *= y.first;
	    }
	    
	    // Acceleration
	    double rfac = 1.0/(rr*rr*rr);
	    
	    for (int k=0; k<3; k++)
	      cC->AddAcc(j, k, -mass *(cC->Pos(j, k) - pos[k]) * rfac );
	    
	    // Potential
#ifdef DEBUG
	    if (use_external) {
	      ncnt++;
	      for (int k=0; k<3; k++)
		tclausius[id] += -mass *
		  (cC->Pos(j, k) - pos[k]) * cC->Pos(j, k) * rfac;
	    }
#endif
	    cC->AddPot(j, pot );
	  }
	}
	// END: spherical point mass
      }
      // END: buffer interaction loop
    }
    // END: local particle loop
  }
  // END: level loop
  
#ifdef DEBUG
  if (use_external) {
    pthread_mutex_lock(&iolock);
    cout << "Process " << myid << ", id=" << id << ": ninteract=" << ninteract
	 << "  nexternal=" << ncnt << "  VC=" << tclausius[id] << endl;
    pthread_mutex_unlock(&iolock);
  }
#endif

  return (NULL);
}

void Direct::determine_coefficients(void) {}
void * Direct::determine_coefficients_thread(void *arg) { return (NULL); }

