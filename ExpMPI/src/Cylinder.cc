#include <values.h>
#include "expand.h"

#include <gaussQ.h>
#include <WghtOrth3.h>
#include <EmpOrth7thd.h>

#include <Orient.H>
#include <Cylinder.H>

static char rcsid[] = "$Id$";


pthread_mutex_t Cylinder::used_lock;
pthread_mutex_t Cylinder::cos_coef_lock;
pthread_mutex_t Cylinder::sin_coef_lock;

Cylinder::Cylinder(string& line) : Basis(line)
{
  geometry = cylinder;

				// Default values
  rcylSL = 20.0;
  rcylEM = 20.0;
  zmax = 8.0;
  acyl = 0.5;
  nmax = 10;
  nfft = 6;
  mmax = 4;
  ncylnx = 32;
  ncylny = 32;
  hcyl = 1.0;
  EJcyl = 0;
  ncylkeep = 100;
  ncylamom = 500;
  ecyl0 = -0.5;
  ncylnzof = -1;
  ncylordz = 6;
  ncylrecomp = -1;
  self_consistent = true;
  selector = false;

  initialize();

  if (rcylSL < rcylEM) {
    cerr << "Cylacp: rcylSL < rcylEM . . . setting rcylSL=rcylEM\n";
    rcylSL = rcylEM;
  }

  CylindricalSL::ZMAX = zmax;
  CylindricalSL::RMAX = rcylSL;
  SLGridCyl::A = acyl;

  orthoS = new CylindricalSL();
  orthoS->reset(nmax, nfft, mmax);

  EmpCylSL::RMAX = rcylEM;
  EmpCylSL::NUMX = ncylnx;
  EmpCylSL::NUMY = ncylny;
  EmpCylSL::HSCALE = hcyl;
  if (ncylnzof>0)
    EmpCylSL::NZOF = ncylnzof;
  else
    EmpCylSL::NZOF = orthoS->get_maxNZ();
				// For debugging
#ifdef DENSITY
  EmpCylSL::DENS = true;
#endif

#ifdef SELECTOR
  EmpCylSL::SELECT = true;
#endif

  ortho = new EmpCylSL();
  ortho->reset(orthoS, zmax, ncylordz);

  ncompcyl = 0;

  if (EJcyl) orient = new Orient(ncylkeep, ncylamom, ecyl0);
    
  if (myid == 0) {		// Flag messages for diagnostics
    
    if (EJcyl & AXIS)
      cout << "Cylinder: AXIS orientation is *ON*\n";
    else
      cout << "Cylinder: AXIS orientation is *OFF*\n";

    if (EJcyl & CENTER)
      cout << "Cylinder: CENTER finding is *ON*\n";
    else
      cout << "Cylinder: CENTER finding is *OFF*\n";
  }
    
  pos = new Vector [nthrds];
  frc = new Vector [nthrds];
  for (int i=0; i<nthrds; i++) {
    pos[i].setsize(1, 3);
    frc[i].setsize(1, 3);
  }

  firstime = true;
}

Cylinder::~Cylinder()
{
  delete [] pos;
  delete [] frc;
}

void Cylinder::initialize()
{
  string val;

  if (get_value("rcylSL", val)) rcylSL = atof(val.c_str());
  if (get_value("rcylEM", val)) rcylEM = atof(val.c_str());
  if (get_value("zmax", val)) zmax = atof(val.c_str());
  if (get_value("acyl", val)) acyl = atof(val.c_str());
  if (get_value("hcyl", val)) hcyl = atof(val.c_str());
  if (get_value("nmax", val)) nmax = atoi(val.c_str());
  if (get_value("nfft", val)) nfft = atoi(val.c_str());
  if (get_value("mmax", val)) mmax = atoi(val.c_str());
  if (get_value("ncylnx", val)) ncylnx = atoi(val.c_str());
  if (get_value("ncylny", val)) ncylny = atoi(val.c_str());
  if (get_value("EJcyl", val)) EJcyl = atoi(val.c_str());
  if (get_value("ncylkeep", val)) ncylkeep = atoi(val.c_str());
  if (get_value("ncylamom", val)) ncylamom = atoi(val.c_str());
  if (get_value("ecyl0", val)) ecyl0 = atof(val.c_str());
  if (get_value("ncylzof", val)) ncylnzof = atoi(val.c_str());
  if (get_value("ncylordz", val)) ncylordz = atoi(val.c_str());
  if (get_value("ncylrecomp", val)) ncylrecomp = atoi(val.c_str());
  if (get_value("self_consistent", val)) {
    if (atoi(val.c_str())) self_consistent = true; 
    else self_consistent = false;
  }
  if (get_value("selector", val)) {
    if (atoi(val.c_str())) selector = true; 
    else selector = false;
  }
}

void Cylinder::get_acceleration_and_potential(vector<Particle>* Particles)
{
  particles = Particles;

				// Try to read cached tables rather
				// than recompute from distribution
  if (firstime && !(restart && ortho->read_cache())) {
    determine_coefficients();
    firstime = false;
  }

  /*======================*/
  /* Reset center of mass */
  /*======================*/

  Vector ctr;
  if (EJcyl & CENTER) {
    ctr = orient->currentCenter();
				// Tie halo center to disk center
    for (int i=0; i<3; i++) component->com[i] = ctr[i+1];
  }

  /*======================*/
  /* Compute coefficients */
  /*======================*/

  if (self_consistent) {
    determine_coefficients();

    if (ncompcyl == 0) {

      cyltime = tnow;

				// Dump basis
#ifdef DENSITY
      if (myid == 0) {
	ortho->dump_basis(outname.c_str(), this_step);
	dump_mzero(outname.c_str(), this_step);
      }
#endif
				// do it from table . . . 
      eof = 1;
      determine_coefficients();

				// In the end, it would be cleaner to call 
 				// for table recomputation separately
    }
    
  }


  /*======================================*/
  /* Determine potential and acceleration */
  /*======================================*/

  MPL_start_timer();

  determine_acceleration_and_potential();

  if (EJcyl) orient->accumulate(particles, &component->com[0]);

  MPL_stop_timer();


  /*=====================================*/
  /* Recompute PCA analysis on next step */
  /*=====================================*/

  ncompcyl++;
  if (ncompcyl == ncylrecomp) {
    ortho->compute_eof();
    ncompcyl = 0;
    eof = 0;
  }


				// Test
  /*
  MPI_Barrier(MPI_COMM_WORLD);
  if (myid==0) {
    cout << "Cyl: n=" << this_step << "  com=" << component->com[0] 
	 << ", " << component->com[1] << ", " << component->com[2] << endl;
  }
  */


  // Debug
  if (myid==1 && EJcyl) 
    {
      string toutfile = string(homedir) + "test.orientation";
      ofstream debugf(toutfile.c_str(), ios::app);
      Vector axis = orient->currentAxis();
      debugf << tnow << " "
	     << orient->currentAxis()[1] << " " 
	     << orient->currentAxis()[2] << " " 
	     << orient->currentAxis()[3] << " " 
	     << orient->currentAxisVar() << " "
	     << orient->currentCenter()[1] << " " 
	     << orient->currentCenter()[2] << " " 
	     << orient->currentCenter()[3] << " " 
	     << orient->currentCenterVar() << " "
	     << orient->currentCenterVarZ() << " "
	     << orient->currentE() << " "
	     << orient->currentUsed()
	     << '\n';
    }

  // Clear external potential flag
  use_external = false;

}

void * Cylinder::determine_coefficients_thread(void * arg)
{
  double r, r2, phi;
  double xx, yy, zz;

  int nbodies = particles->size();
  int id = *((int*)arg);
  int nbeg = nbodies*id/nthrds;
  int nend = nbodies*(id+1)/nthrds;

  use[id] = 0;
  cylmass0[id] = 0.0;

  for (int i=nbeg; i<nend; i++) {

				// frozen particles don't contribute to field
    if ((*particles)[i].freeze()) continue;

    for (int j=0; j<3; j++) 
      pos[id][j+1] = (*particles)[i].pos[j] - component->com[j];

    if (EJcyl & AXIS) pos[id] = orient->transformBody() * pos[id];

    xx = pos[id][1];
    yy = pos[id][2];
    zz = pos[id][3];

    r2 = xx*xx + yy*yy;
    r = sqrt(r2) + DSMALL;

    if (r<=rcylEM && fabs(zz)<=zmax) {
      if (eof) {
	use[id]++;
	cylmass0[id] += (*particles)[i].mass;
      }

      phi = atan2(yy, xx);
      ortho->accumulate(r, zz, phi, (*particles)[i].mass, id);
    }
    
  }

  return (NULL);
}


void Cylinder::determine_coefficients(void)
{
  static char routine[] = "determine_coefficients_Cylacp";
  int i;
  int use0, use1;
  double cylmassT;

#ifdef MPE_PROFILE
  MPE_Log_event(9, myid, "b_compute_coef");
#endif

  use0 = 0;
  use1 = 0;
  cylmassT = 0.0;
  cylmass = 0.0;

  ortho->setup_accumulation();
    
  cylmass0 = new double [nthrds];
  if (!cylmass0) {
    cerr << "Cylacp: problem allocating <cylmass0>\n";
    exit(-1);
  }


				/* Initialize locks */
  make_mutex(&used_lock, routine, "used_lock");
  make_mutex(&cos_coef_lock, routine, "cos_coef_lock");
  make_mutex(&sin_coef_lock, routine, "sin_coef_lock");

  exp_thread_fork(true);
  
  kill_mutex(&used_lock, routine, "used_lock");
  kill_mutex(&cos_coef_lock, routine, "cos_coef_lock");
  kill_mutex(&sin_coef_lock, routine, "sin_coef_lock");

  for (i=0; i<nthrds; i++) {
    use1 += use[i];
    cylmassT += cylmass0[i];
  }

  delete [] cylmass0;
				// Turn off timer so as not bias by 
				// communication barrier
  MPL_stop_timer();

  ortho->make_coefficients();

  MPI_Allreduce ( &use1, &use0,  1, MPI_INT, MPI_SUM, MPI_COMM_WORLD);
  used = use0;
  cout << "Process " << myid << ": used=" << PotAccel::used << endl;

  MPI_Allreduce ( &cylmassT, &cylmass,  1, MPI_DOUBLE, MPI_SUM, MPI_COMM_WORLD);

#ifdef MPE_PROFILE
  MPE_Log_event(10, myid, "e_compute_coef");
#endif

  MPL_start_timer();

}

void check_force_values(double phi, double p, double fr, double fz, double fp)
{
  if (
      isinf(phi) || isnan(phi) ||
      isinf(p  ) || isnan(p  ) ||
      isinf(fr ) || isnan(fr ) ||
      isinf(fz ) || isnan(fz ) ||
      isinf(fp ) || isnan(fp ) ) 
    {
      cerr << "check_force_values: Illegal value\n";
    }
}


void * Cylinder::determine_acceleration_and_potential_thread(void * arg)
{
  int i;
  double r, r2, r3, phi;
  double xx, yy, zz;
  double p, fr, fz, fp;

  int nbodies = particles->size();
  int id = *((int*)arg);
  int nbeg = nbodies*id/nthrds;
  int nend = nbodies*(id+1)/nthrds;

  for (i=nbeg; i<nend; i++) {

    if ((*particles)[i].freeze()) continue;

    for (int j=0; j<3; j++) 
      pos[id][j+1] = (*particles)[i].pos[j] - component->com[j];

    if (EJcyl & AXIS) pos[id] = orient->transformBody() * pos[id];

    xx = pos[id][1];
    yy = pos[id][2];
    zz = pos[id][3];

    r2 = xx*xx + yy*yy;
    r = sqrt(r2) + DSMALL;
    phi = atan2(yy, xx);

    if (r<=rcylEM && fabs(zz)<=zmax) {

      ortho->accumulated_eval(r, zz, phi, p, fr, fz, fp);
    
#ifdef DEBUG
      check_force_values(phi, p, fr, fz, fp);
#endif

      if (use_external)
	(*particles)[i].potext += p;
      else
	(*particles)[i].pot += p;

      frc[id][1] = fr*xx/r - fp*yy/r2;
      frc[id][2] = fr*yy/r + fp*xx/r2;
      frc[id][3] = fz;

      if (EJcyl & AXIS) frc[id] = orient->transformOrig() * frc[id];

      for (int j=0; j<3; j++) (*particles)[i].acc[j] += frc[id][j+1];
    }
    else {

      /*
      cerr.form("Process %d: i=%d r=%f zz=%f Cylmass=%f\n", 
		  myid, i, r, zz, cylmass);
      */

      r3 = r2 + zz*zz;
      p = -cylmass/sqrt(r3);	// -M/r
      fr = p/r3;		// -M/r^3

      if (use_external)
	(*particles)[i].potext += p;
      else
	(*particles)[i].pot += p;

      (*particles)[i].acc[0] += xx * fr;
      (*particles)[i].acc[1] += yy * fr;
      (*particles)[i].acc[2] += zz * fr;
    }

  }

  return (NULL);
}


void Cylinder::determine_acceleration_and_potential(void)
{
  static char routine[] = "determine_acceleration_and_potential_Cyl";
  
#ifdef MPE_PROFILE
  MPE_Log_event(11, myid, "b_compute_force");
#endif

  exp_thread_fork(false);

#ifdef MPE_PROFILE
  MPE_Log_event(12, myid, "e_compute_force");
#endif

}


void Cylinder::determine_fields_at_point_cyl(double r, double z, double phi,
					     double *tdens, double *tpotl, 
					     double *tpotr, double *tpotz, double *tpotp)
{
  ortho->accumulated_eval(r, z, phi, *tpotl, *tpotr, *tpotz, *tpotp);
#ifdef DENSITY
  *tdens = ortho->accumulated_dens_eval(r, z, phi);
#else
  *tdens = 0.0;
#endif
}

				// Dump coefficients to a file
void Cylinder::dump_coefs(ostream& out)
{
  /*
  fprintf(fout, "Time=%f\n", tnow);
  ortho->dump_coefs(fout);
  */
  /*
  ortho->dump_coefs_binary_last(out, cyltime);
  */
  ortho->dump_coefs_binary_curr(out, tnow);
}

				// Density debug
#include <fstream.h>
#include <strstream.h>

void Cylinder::dump_mzero(const string& name, int step)
{
  const double RMAX = 5.0*acyl;
  const double ZMAX = 5.0*hcyl;
  double r, dr = RMAX/(ncylnx-1);
  double z, dz = 2.0*ZMAX/(ncylny-1);

  float zz;
  char strbuf[128];
  string label[] = {".dens0.", ".pot0.", ".fr0.", ".fz0."};
  ofstream** out = new ofstream* [4];

  for (int i=0; i<4; i++) {
    ostrstream ins(strbuf, 128);
    ins << name << label[i] << step << '\0';
    out[i] = new ofstream(strbuf);

    out[i]->write((char *)&ncylnx, sizeof(int));
    out[i]->write((char *)&ncylny, sizeof(int));
    out[i]->write((char *)&(zz=  0.0), sizeof(float));
    out[i]->write((char *)&(zz= RMAX), sizeof(float));
    out[i]->write((char *)&(zz=-ZMAX), sizeof(float));
    out[i]->write((char *)&(zz= ZMAX), sizeof(float));
  }


				// Ok, write data
  double p, fr, fz, fp;

  for (int k=0; k<ncylny; k++) {

    z = -ZMAX + dz*k;
	
    for (int j=0; j<ncylnx; j++) {
	  
      r = dr*j;

      zz = ortho->accumulated_dens_eval(r, z, 0.0);
      out[0]->write((char *)&zz, sizeof(float));

      ortho->accumulated_eval(r, z, 0.0, p, fr, fz, fp);
      out[1]->write((char *)&(zz=p ), sizeof(float));
      out[2]->write((char *)&(zz=fr), sizeof(float));
      out[3]->write((char *)&(zz=fz), sizeof(float));
    }
  }

				// Close and delete streams
  for (int i=0; i<4; i++) {
    out[i]->close();
    delete out[i];
  }
  delete [] out;

}
