// This is really C++ code
/*****************************************************************************
 *  Description:
 *  -----------
 *
 *  Get V''[r(t)] for radial orbit (use for tidal forcing)
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
#include <iostream.h>
#include <iomanip.h>
#include <fstream.h>
#include <string>
#include <Vector.h>

#include <orbit.h>
#include <massmodel.h>
#include <numerical.h>

double fr(double);


void error(char *s, char *s2 = "")
{
  cerr << s << ' ' << s2 << '\n';
}


				// Input parameters
static double E=-0.5;
static double K=1.0e-4;
static double PER=0.25;
static double AMPL=1.0;
static char INFILE[128] = "w05";

				// Local variables

static int first_time=1;
static SphericalModelTable *model;
static SphericalOrbit *t;

				// Local functions
static void read_parm(char *file);

void global_setup_shock(void)
{
  char filename[] = "shock.params";
  read_parm(filename);

  //
  // Initialize model and orbit
  //

  model = new SphericalModelTable(INFILE);

  t = new SphericalOrbit(model, E, K);

}


double tidal_shock(double T)
{
  if (first_time) {
    global_setup_shock();
    first_time = 0;
  }

  return AMPL * model->get_dpot2(t->get_angle(6, T*PER));
}


extern "C" {
  double get_tidal_shock(double T)
    {
      if (first_time) {
	global_setup_shock();
	first_time = 0;
      }

      return AMPL * model->get_dpot2(t->get_angle(6, T*PER));
    }
};


//----------------------------------------------------------------------
//----------------------------------------------------------------------
//----------------------------------------------------------------------

extern "C" {
  int get_key_value_from_file(char *, char ***, char ***);
}

static void set_parm(char *word, char *valu);

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
  if (!strcmp("E",word))
    E = atof(valu);

  if (!strcmp("K",word))
    K = atof(valu);

  if (!strcmp("PER",word))
    PER = atof(valu);

  if (!strcmp("AMPL",word))
    AMPL = atof(valu);

  if (!strcmp("INFILE",word))
    strcpy(INFILE,valu);

}

