// This is really C++ code
/*****************************************************************************
 *  Description:
 *  -----------
 *
 *  Get V' for a spherical "halo" model
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
 *
 *  Notes:
 *  -----
 *
 *
 *  By:
 *  --
 *
 *  MDW 11/20/91
 *
 ***************************************************************************/

#include <stdlib.h>
#include <math.h>
#include <iostream.h>
#include <iomanip.h>
#include <fstream.h>
#include <string>
#include <Vector.h>

#ifdef USE_DMALLOC
#include <dmalloc.h>
#endif

#include <numerical.h>
#include <orbit.h>
#include <massmodel.h>
#include <model3d.h>
#include <isothermal.h>
#include <hernquist.h>

int parse_args(int, char **);


				// Input parameters
static int HMODEL = file;
static string INFILE = "w05";

static double MHALO=1.0;
static double RHALO=1.0;
static double RMODMIN=1.0e-3;
static double RMOD=20.0;

static double RBCORE=1.0;
static double MBULGE=1.0;
static double RBULGE=1.0;
static double RBMODMIN=1.0e-3;
static double RBMOD=20.0;
				// Global variables

static int first_time=1;
static AxiSymModel *model, *bmodel;

static void read_parm(char *file);

void global_setup_halo(void)
{
  char filename[] = "halo.params";
  read_parm(filename);

  //
  // Initialize model and orbit
  //

  switch (HMODEL) {
  case file:
    model = new SphericalModelTable(INFILE) // Halo model
    break;
  case isothermal:
    model = new IsothermalSphere(1.0, RMODMIN, RMOD); // Halo model
    break;
  case hernquist_model:
    model = new HernquistSphere(1.0, RMODMIN, RMOD); // Halo model
    break; 
  default:
    cerr << "No such HALO model type: " << (int)HMODEL << endl;
    exit(-2);
  }

  bmodel = new HernquistSphere(RBCORE, RBMODMIN, RBMOD);
}


extern "C" {

  extern double *mass,*x,*y,*z,*rr,*vx,*vy,*vz,*ax,*ay,*az,*pot,*potext,*size;
  extern int *component;
  extern double tpos,tvel,tnow;
  extern double **expcoef;
  extern int nbodies,used,ninteract,nmerge;
  extern int restart;
  extern int freeze_particle(int);

  void external_halo(void)
    {
      if (first_time) {
	global_setup_halo();
	first_time = 0;
      }

      double r, potl, dpot, potlB, dpotB;

      for (int i=1; i<=nbodies; i++) {
	if (freeze_particle(i)) continue;
	
				// Act on disk (Component #2) only
	if (component[i] != 2) continue;

	r = sqrt(x[i]*x[i] + y[i]*y[i] + z[i]*z[i]);

	model->get_pot_dpot(r/RHALO, potl, dpot);
	potl *= MHALO/RHALO;
	dpot *= MHALO/RHALO/RHALO;

	bmodel->get_pot_dpot(r/RBULGE, potlB, dpotB);
	potlB *= MBULGE/RBULGE;
	dpotB *= MBULGE/RBULGE/RBULGE;

	ax[i] += -(dpot + dpotB)*x[i]/r;
	ay[i] += -(dpot + dpotB)*y[i]/r;
	az[i] += -(dpot + dpotB)*z[i]/r;
	potext[i] += potl + potlB;
      }
    }
};


//----------------------------------------------------------------------
//----------------------------------------------------------------------
//----------------------------------------------------------------------

static void set_parm(char *, char *);

extern "C" {
  int get_key_value_from_file(char *, char ***, char ***);
}

static void read_parm(char *file)
{
  int iret, i;
  char **word, **valu;

  iret = get_key_value_from_file(file,&word,&valu);
	
			/* Set parameters */

  for (i=0; i<iret; i++)
    set_parm(word[i],valu[i]);
}

static void set_parm(char *word, char *valu)
{
  if (!strcmp("HMODEL",word))
    HMODEL = atoi(valu);

  if (!strcmp("INFILE",word))
    INFILE = valu;

  if (!strcmp("MHALO",word))
    MHALO = atof(valu);

  if (!strcmp("RHALO",word))
    RHALO = atof(valu);

  if (!strcmp("RMODMIN",word))
    RMODMIN = atof(valu);

  if (!strcmp("RMOD",word))
    RMOD = atof(valu);

  if (!strcmp("RBCORE",word))
    RBCORE = atof(valu);

  if (!strcmp("MBULGE",word))
    MBULGE = atof(valu);

  if (!strcmp("RBULGE",word))
    RBULGE = atof(valu);

  if (!strcmp("RBMODMIN",word))
    RBMODMIN = atof(valu);

  if (!strcmp("RBMOD",word))
    RBMOD = atof(valu);

}






