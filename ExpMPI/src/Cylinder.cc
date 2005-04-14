static char rcsid[] = "$Id$";

using namespace std;

#include <values.h>
#include "expand.h"

#include <sstream>

#include <gaussQ.h>
#include <EmpOrth9thd.h>

#include <Orient.H>
#include <Cylinder.H>

pthread_mutex_t Cylinder::used_lock;
pthread_mutex_t Cylinder::cos_coef_lock;
pthread_mutex_t Cylinder::sin_coef_lock;

Cylinder::Cylinder(string& line) : Basis(line)
{
  id = "Cylinder";
  geometry = cylinder;

				// Default values

  rcylmin = 0.001;		// Should only change these two in
  rcylmax = 20.0;		// extreme circumstances

  ncylnx = 128;			// These defaults should do fine in
  ncylny = 64;			// most cases, as well

  acyl = 1.0;
  nmax = 10;
  lmax = 36;
  mmax = 4;
  hcyl = 1.0;
  ncylorder = 10;
  ncylrecomp = -1;
  hallfile = "disk";
  hallfreq = 50;
  self_consistent = true;
  selector = false;
  density = false;
  coef_dump = true;

  initialize();


  EmpCylSL::RMIN = rcylmin;
  EmpCylSL::RMAX = rcylmax;
  EmpCylSL::NUMX = ncylnx;
  EmpCylSL::NUMY = ncylny;
  EmpCylSL::CMAP = true;	// Always use coordinate mapping

				// For debugging
  if (density) EmpCylSL::DENS = true;

  /*
  ortho = new EmpCylSL();
  ortho->reset(nmax, lmax, mmax, ncylorder, acyl, hcyl);
  */

  ortho = new EmpCylSL(nmax, lmax, mmax, ncylorder, acyl, hcyl);

  if (selector) {
    EmpCylSL::SELECT = true;
    ortho->setHall(hallfile, hallfreq);
  }

  cout << "Process " << myid << ": Cylinder parameters: "
       << " nmax=" << nmax
       << " lmax=" << lmax
       << " mmax=" << mmax
       << " ncylorder=" << ncylorder
       << " rcylmin=" << rcylmin
       << " rcylmax=" << rcylmax
       << " acyl=" << acyl
       << " hcyl=" << hcyl
       << " selector=" << selector
       << " hallfreq=" << hallfreq
       << " hallfile=" << hallfile
       << "\n";

  ncompcyl = 0;

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

  // These should not be user settable . . . but need them for now
  if (get_value("rcylmin", val)) rcylmin = atof(val.c_str());
  if (get_value("rcylmax", val)) rcylmax = atof(val.c_str());

  if (get_value("acyl", val)) acyl = atof(val.c_str());
  if (get_value("hcyl", val)) hcyl = atof(val.c_str());
  if (get_value("nmax", val)) nmax = atoi(val.c_str());
  if (get_value("lmax", val)) lmax = atoi(val.c_str());
  if (get_value("mmax", val)) mmax = atoi(val.c_str());
  if (get_value("ncylnx", val)) ncylnx = atoi(val.c_str());
  if (get_value("ncylny", val)) ncylny = atoi(val.c_str());
  if (get_value("ncylorder", val)) ncylorder = atoi(val.c_str());
  if (get_value("ncylrecomp", val)) ncylrecomp = atoi(val.c_str());
  if (get_value("hallfreq", val)) hallfreq = atoi(val.c_str());
  if (get_value("hallfile", val)) hallfile = val;
  if (get_value("self_consistent", val)) {
    if (atoi(val.c_str())) self_consistent = true; 
    else self_consistent = false;
  }
  if (get_value("selector", val)) {
    if (atoi(val.c_str())) selector = true; 
    else selector = false;
  }
  if (get_value("density", val)) {
    if (atoi(val.c_str())) density = true; 
    else density = false;
  }
}

void Cylinder::get_acceleration_and_potential(Component* C)
{
  cC = C;

				// External particles only
  if (use_external) {
    
    MPL_start_timer();
    determine_acceleration_and_potential();
    MPL_stop_timer();

    use_external = false;

    return;
  }

  //==================================
  // Try to read cached tables rather
  // than recompute from distribution
  //==================================
  
  if (firstime) {

    if (!(restart && ortho->read_cache()))
      determine_coefficients();	// Will compute EOF grid


    if (!self_consistent) {	// Compute mass and coefficients
      eof = 1;			// once and for all
      determine_coefficients();
    }

    firstime = false;
  }

  //======================
  // Compute coefficients 
  //======================

  if (self_consistent) {
    
    determine_coefficients();

    if (ncompcyl == 0) {

      cyltime = tnow;

				// Dump basis
#ifdef DENSITY
      if (myid == 0) {
	ortho->dump_basis(runtag.c_str(), this_step);
	dump_mzero(runtag.c_str(), this_step);
      }
#endif
				// Get coefficients from table if this
				// is an EOF recomputation
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


  // Debug
  /*
  if (myid==1 && component->EJ) 
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
  */

}

void * Cylinder::determine_coefficients_thread(void * arg)
{
  double r, r2, phi;
  double xx, yy, zz;

  unsigned nbodies = cC->Number();
  int id = *((int*)arg);
  int nbeg = nbodies*id/nthrds;
  int nend = nbodies*(id+1)/nthrds;
  double adb = component->Adiabatic();

  use[id] = 0;
  cylmass0[id] = 0.0;

  for (int i=nbeg; i<nend; i++) {

				// frozen particles don't contribute to field
    if (cC->freeze(*(cC->Part(i)))) continue;

    for (int j=0; j<3; j++) 
      pos[id][j+1] = cC->Pos(i, j, Component::Local | Component::Centered);

    if (cC->EJ & Orient::AXIS) 
      pos[id] = cC->orient->transformBody() * pos[id];

    xx = pos[id][1];
    yy = pos[id][2];
    zz = pos[id][3];

    r2 = xx*xx + yy*yy;
    r = sqrt(r2) + DSMALL;

    if (r2+zz*zz < rcylmax*rcylmax) {

      double mas = cC->Mass(i) * adb;

      if (eof) {
	use[id]++;
	cylmass0[id] += mas;
      }

      phi = atan2(yy, xx);
      ortho->accumulate(r, zz, phi, mas, id);
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

  unsigned nbodies = cC->Number();
  int id = *((int*)arg);
  int nbeg = nbodies*id/nthrds;
  int nend = nbodies*(id+1)/nthrds;

  for (i=nbeg; i<nend; i++) {

    if (use_external) {
      cC->Pos(&pos[id][1], i, Component::Inertial);
      component->ConvertPos(&pos[id][1], Component::Local | Component::Centered);
    } else
      cC->Pos(&pos[id][1], Component::Local | Component::Centered);

    if (cC->EJ & Orient::AXIS) 
      pos[id] = cC->orient->transformBody() * pos[id];

    xx = pos[id][1];
    yy = pos[id][2];
    zz = pos[id][3];

    r2 = xx*xx + yy*yy;
    r = sqrt(r2) + DSMALL;
    phi = atan2(yy, xx);

    if (r2 + zz*zz < rcylmax*rcylmax) {

      ortho->accumulated_eval(r, zz, phi, p, fr, fz, fp);
    
#ifdef DEBUG
      check_force_values(phi, p, fr, fz, fp);
#endif

      if (use_external)
	cC->AddPotExt(i, p);
      else
	cC->AddPot(i, p);

      frc[id][1] = fr*xx/r - fp*yy/r2;
      frc[id][2] = fr*yy/r + fp*xx/r2;
      frc[id][3] = fz;

      if (cC->EJ & Orient::AXIS) 
	frc[id] = cC->orient->transformOrig() * frc[id];

      for (int j=0; j<3; j++) cC->AddAcc(i, j, frc[id][j+1]);
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
	cC->AddPotExt(i, p);
      else
	cC->AddPot(i, p);

      cC->AddAcc(i, 0, xx*fr);
      cC->AddAcc(i, 1, yy*fr);
      cC->AddAcc(i, 2, zz*fr);
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

void Cylinder::
determine_fields_at_point_sph(double r, double theta, double phi,
			      double *tdens, double *tpotl, 
			      double *tpotr, double *tpott, 
			      double *tpotp)

{
  double R = r*sin(theta);
  double z = r*cos(theta);
  double tpotR, tpotZ;

  determine_fields_at_point_cyl(R, z, phi, tdens, tpotl, 
				&tpotR, &tpotZ, tpotp);
  
  *tpotr =   tpotR*sin(theta) + tpotZ*cos(theta) ;
  *tpott = (-tpotZ*sin(theta) + tpotR*cos(theta) )/(r+1.0e-10);
}



void Cylinder::determine_fields_at_point_cyl(double r, double z, double phi,
					     double *tdens, double *tpotl, 
					     double *tpotr, double *tpotz, double *tpotp)
{
  ortho->accumulated_eval(r, z, phi, *tpotl, *tpotr, *tpotz, *tpotp);
  // Accumulated eval returns forces not potential gradients
  *tpotr *= -1.0;
  *tpotz *= -1.0;
  *tpotp *= -1.0;
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
#include <fstream>

void Cylinder::dump_mzero(const string& name, int step)
{
  const double RMAX = 5.0*acyl;
  const double ZMAX = 5.0*hcyl;
  double r, dr = RMAX/(ncylnx-1);
  double z, dz = 2.0*ZMAX/(ncylny-1);

  float zz;
  string label[] = {".dens0.", ".pot0.", ".fr0.", ".fz0."};
  ofstream** out = new ofstream* [4];

  for (int i=0; i<4; i++) {
    ostringstream ins;
    ins << name << label[i] << step;
    out[i] = new ofstream(ins.str().c_str());

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
