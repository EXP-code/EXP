// This may look like C code, but it is really -*- C++ -*-

// #define DEBUG 1

/*****************************************************************************
 *  Description:
 *  -----------
 *
 *  This routine computes the potential, acceleration and density using
 *  the Cylindrical biorthogonal expansion
 *
 *
 *  Call sequence:
 *  -------------
 *
 *  Parameters:
 *  ----------
 *
 *
 *  Returns:
 *  -------
 *
 *  Value
 *
 *  Notes:
 *  -----
 *
 *	KL 5/27/92   Modified to allow freezing of particles beyond some cutoff
 *    radius. This is needed when running systems in tidal fields. 
 *    The function freeze_particle() is called for each particle
 *    to decide whether to freeze. *
 *
 *  By:
 *  --
 *
 *  MDW 11/13/91
 *      06/09/92 updated to use recursion relations rather than tables
 *
 ***************************************************************************/

#define DENSITY
// #define SELECTOR

#include <values.h>
#include <gaussQ.h>
#include <WghtOrth3.h>
#include <EmpOrth7thd.h>
#include <Orient.h>

#include "expand.h"
#include "exp_thread.h"

#ifdef RCSID
static char rcsid[] = "$Id$";
#endif

void determine_coefficients_Cyl(void);
void dump_mzero(char *, int);


extern "C" {
  void determine_acceleration_and_potential_Cyl(void);
  void determine_fields_at_point_Cyl(double r, double z, double phi,
				     double *tdens, double *tpotl, 
				     double *tpotr, double *tpotz, 
				     double *tpotp);
}

//======================================================================
//======================================================================
//======================================================================

enum OrientFlags {AXIS=1, CENTER=2};

static Orient* orient;
static EmpCylSL *ortho;
static CylindricalSL *orthoS;
static int eof;

static Vector pos, frc;

extern "C" 
void get_acceleration_and_potential_Cyl(void)
{
  static int firstime=1;

  if (firstime) {

    if (rcylSL < rcylEM) {
      cerr << "Cylacp: rcylSL < rcylEM . . . setting rcylSL=rcylEM\n";
      rcylSL = rcylEM;
    }

    CylindricalSL::ZMAX = zmax;
    CylindricalSL::RMAX = rcylSL;
    SLGridCyl::A = acyl;

    orthoS = new CylindricalSL();
    orthoS->reset(nmax2, nfft, mmax2);

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
	cout << "Cylacp3: AXIS orientation is *ON*\n";
      else
	cout << "Cylacp3: AXIS orientation is *OFF*\n";

      if (EJcyl & CENTER)
	cout << "Cylacp3: CENTER finding is *ON*\n";
      else
	cout << "Cylacp3: CENTER finding is *OFF*\n";
    }

    pos.setsize(1, 3);
    frc.setsize(1, 3);
  }

  /*======================*/
  /* Reset center of mass */
  /*======================*/

  Vector ctr;
  if (EJcyl & CENTER) {
    ctr = orient->currentCenter();
    // for (int i=0; i<3; i++) com2[i] = ctr[i+1];
				// Tie halo center to disk center
    for (int i=0; i<3; i++) com1[i] = com2[i] = ctr[i+1];
  }

  /*======================*/
  /* Compute coefficients */
  /*======================*/

  if (firstime || self_consistent) {
				// Try to read cached tables rather
				// than recompute from distribution
    if (!(restart && firstime && ortho->read_cache()))
      determine_coefficients_Cyl();

    firstime = 0;
    
    if (ncompcyl == 0) {

      cyltime = tnow;

				// Dump basis
#ifdef DENSITY
      if (myid == 0) {
	ortho->dump_basis(outname, this_step);
	dump_mzero(outname, this_step);
      }
#endif
				// do it from table . . . 
      eof = 1;
      determine_coefficients_Cyl();
				// In the end, it would be cleaner to call 
 				// for table recomputation separately
    }
    
  }


  /*======================================*/
  /* Determine potential and acceleration */
  /*======================================*/

  MPL_start_timer();

  if (myid>0) determine_acceleration_and_potential_Cyl();
  if (myid>0 && EJcyl) orient->accumulate(component, mass, pot, 
					  x, y, z, vx, vy, vz, com2, nbodies);

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

  /*
				// Test
  if (myid==0) {
    cout << "Cyl: n=" << this_step << "  com2=" << com2[0] 
	 << ", " << com2[1] << ", " << com2[2] << endl;
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

}

pthread_mutex_t used_lock, cos_coef_lock, sin_coef_lock;
static int *use;
static double *cylmass0;

void * determine_coefficients_Cyl_thread(void * arg)
{
  int i;
  double r, r2, phi;
  double xx, yy, zz;

  int id = *((int*)arg);
  int nbeg = 1+nbodies*id/nthrds;
  int nend = nbodies*(id+1)/nthrds;

  use[id] = 0;
  cylmass0[id] = 0.0;

  for (i=nbeg; i<=nend; i++) {


    if (freeze_particle(i)) continue;		/* frozen particles don't
						   contribute to field.
						   KL 5/27/92 */
    pos[1] = x[i] - com2[0];
    pos[2] = y[i] - com2[1];
    pos[3] = z[i] - com2[2];

    if (EJcyl & AXIS) pos = orient->transformBody() * pos;

    xx = pos[1];
    yy = pos[2];
    zz = pos[3];

    r2 = xx*xx + yy*yy;
    r = rr[i] = sqrt(r2) + DSMALL;

    if (component[i] != 2) continue;
    
    if (r<=rcylEM && fabs(zz)<=zmax) {
      if (eof) {
	use[id]++;
	cylmass0[id] += mass[i];
      }

      phi = atan2(yy, xx);
      ortho->accumulate(r, zz, phi, mass[i], id);
    }
    
  }

  return (NULL);
}


void determine_coefficients_Cyl(void)
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
    
  if (myid>0) {

    use = new int [nthrds];
    if (!use) {
      cerr << "Cylacp: problem allocating <use>\n";
      exit(-1);
    }

    cylmass0 = new double [nthrds];
    if (!cylmass0) {
      cerr << "Cylacp: problem allocating <cylmass0>\n";
      exit(-1);
    }


				/* Initialize locks */
    make_mutex(&used_lock, routine, "used_lock");
    make_mutex(&cos_coef_lock, routine, "cos_coef_lock");
    make_mutex(&sin_coef_lock, routine, "sin_coef_lock");

    exp_thread_fork(determine_coefficients_Cyl_thread, routine);

    kill_mutex(&used_lock, routine, "used_lock");
    kill_mutex(&cos_coef_lock, routine, "cos_coef_lock");
    kill_mutex(&sin_coef_lock, routine, "sin_coef_lock");

    for (i=0; i<nthrds; i++) {
      use1 += use[i];
      cylmassT += cylmass0[i];
    }

    delete [] use;
    delete [] cylmass0;
  }

				// Turn off timer so as not bias by 
				// communication barrier
  MPL_stop_timer();

  ortho->make_coefficients();

  MPI_Allreduce ( &use1, &use0,  1, MPI_INT, MPI_SUM, MPI_COMM_WORLD);
  if (myid==0) used += use0;

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


void * determine_acceleration_and_potential_Cyl_thread(void * arg)
{
  int i;
  double r, r2, r3, phi;
  double xx, yy, zz;
  double p, fr, fz, fp;

  int id = *((int*)arg);
  int nbeg = 1+nbodies*id/nthrds;
  int nend = nbodies*(id+1)/nthrds;


  for (i=nbeg; i<=nend; i++) {

    if (!disk_on_halo && component[i] != 2) continue;

    if (freeze_particle(i)) continue;

    pos[1] = x[i] - com2[0];
    pos[2] = y[i] - com2[1];
    pos[3] = z[i] - com2[2];

    if (EJcyl & AXIS) pos = orient->transformBody() * pos;

    xx = pos[1];
    yy = pos[2];
    zz = pos[3];

    r2 = xx*xx + yy*yy;
    r = rr[i] = sqrt(r2) + DSMALL;
    phi = atan2(yy, xx);

    if (r<=rcylEM && fabs(zz)<=zmax) {

      ortho->accumulated_eval(r, zz, phi, p, fr, fz, fp);
    
#ifdef DEBUG
      check_force_values(phi, p, fr, fz, fp);
#endif

      pot[i] += p;

      frc[1] = fr*xx/r - fp*yy/r2;
      frc[2] = fr*yy/r + fp*xx/r2;
      frc[3] = fz;

      if (EJcyl & AXIS) frc = orient->transformOrig() * frc;

      ax[i] += frc[1];
      ay[i] += frc[2];
      az[i] += frc[3];
    }
    else {

      /*
      if (component[i]==0) {
	cerr.form("Process %d: i=%d r=%f zz=%f Cylmass=%f\n", 
		  myid, i, r, zz, cylmass);
      }
      */

      r3 = r2 + zz*zz;
      p = -cylmass/sqrt(r3);	// -M/r
      fr = p/r3;		// -M/r^3

      pot[i] += p;

      ax[i] += xx * fr;
      ay[i] += yy * fr;
      az[i] += zz * fr;
    }

  }

  return (NULL);
}


void determine_acceleration_and_potential_Cyl(void)
{
  static char routine[] = "determine_acceleration_and_potential_Cyl";
  
#ifdef MPE_PROFILE
  MPE_Log_event(11, myid, "b_compute_force");
#endif

  exp_thread_fork(determine_acceleration_and_potential_Cyl_thread, routine);

#ifdef MPE_PROFILE
  MPE_Log_event(12, myid, "e_compute_force");
#endif

}


extern "C" 
void determine_fields_at_point_Cyl(double r, double z, double phi,
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

				/* Dump coefficients to a file */

extern "C" void dump_coefs_Cyl(FILE *fout)
{
  /*
  fprintf(fout, "Time=%f\n", tnow);
  ortho->dump_coefs(fout);
  */
  /*
  ortho->dump_coefs_binary_last(fout, cyltime);
  */
  ortho->dump_coefs_binary_curr(fout, tnow);
}

				/* Density debug */

#include <fstream.h>
#include <strstream.h>

void dump_mzero(char* name, int step)
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
