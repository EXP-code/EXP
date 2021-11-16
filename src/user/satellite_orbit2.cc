/*****************************************************************************
 *  Description:
 *  -----------
 *
 *  Follow satellite orbit
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
 *  revision MDW 01/10/93
 *  updated to use orbit classes
 *           MDW 07/15/94
 *
 ***************************************************************************/

static char rcsid_satellite_orbit[] = "$Id$";

#include <iostream>
#include <iomanip>
#include <cstdlib>
#include <cmath>
#include <string>

#ifdef USE_DMALLOC
#include <dmalloc.h>
#endif

#include <orbit.H>
#include <massmodel.H>

#include <model2d.H>
#include <model3d.H>
#include <isothermal.H>
#include <hernquist.H>


void parse_args(void);
void write_parm(void);

Matrix return_euler(double PHI, double THETA, double PSI, int BODY);

				// Input parameters
Models3d HALO_MODEL=file;
double E=-1.0;
double K=0.5;
double INCLINE=45.0;
double PSI=0.0;
double PHIP=0.0;
double VROT=1.0;
double RCORE=1.0;
double RMODMIN=1.0e-3;
double RMODMAX=20.0;
double RA=1.0e20;
int DIVERGE=0;
int NUMDF=100;
int NRECS=512;
double DIVERGE_RFAC=1.0;
string MODFILE="halo.model";

//

extern "C" void satellite_orbit(double T, double* X, double* Y, double* Z)
{
  //
  // Begin
  //

  static SphericalOrbit orb;
  static Eigen::Matrix3d rotate;
  static int firstime=1;
  static Eigen::Vector3d v0, v1;


  if (firstime) {

    parse_args();

// ===================================================================
// Initilize HALO model
// ===================================================================
    
    std::shared_ptr<SphericalModelTable> m;
    AxiSymModPtr halo_model;

    switch (HALO_MODEL) {
    case file:
      m = std::make_shared<SphericalModelTable>(MODFILE, DIVERGE, DIVERGE_RFAC);
      m->setup_df(NUMDF, RA);
      halo_model = m;
      Model3dNames[0] = MODFILE;	// Assign filename to ID string
      break;

    case sing_isothermal:
      halo_model = std::make_shared<SingIsothermalSphere>(1.0, RMODMIN, RMODMAX);
      break;

    case isothermal:
      halo_model = std::make_shared<IsothermalSphere>(RCORE, RMODMAX, VROT);
      break;

    case hernquist_model:
      halo_model = std::make_shared<HernquistSphere>(1.0, RMODMIN, RMODMAX);
      break; 

    default:
      cerr << "Illegal HALO model: " << HALO_MODEL << '\n';
      exit(-1);
    }
    

// ===================================================================
// Setup orbit
// ===================================================================

    SphericalOrbit create(halo_model, E, K);
    create.set_numerical_params(NRECS);

    orb = create;

    INCLINE *= M_PI/180.0;
    PSI *= M_PI/180.0;
    PHIP *= M_PI/180.0;

    rotate = return_euler(PHIP, INCLINE, PSI, 1);

    firstime = 0;
    
    cerr << "Satellite_orbit initialized [" << E << ", " << K << endl;
  }


  double r = orb.get_angle(6, T);
  double phi = orb.get_angle(7, T);

  v0[1] = r*cos(phi);
  v0[2] = r*sin(phi);
  v0[3] = 0.0;

  v1 = rotate*v0;

  *X = v1[1];
  *Y = v1[2];
  *Z = v1[3];

}


//----------------------------------------------------------------------
//----------------------------------------------------------------------
//----------------------------------------------------------------------

void usage(char *);
void set_parm(char *, char *);
void write_parm(void);
void print_parm(ostream &, char *);
void print_default(void);

extern "C" {
  int get_key_value(int, char **, char ***, char ***);
  int get_key_value_from_file(char *, char ***, char ***);
  void std_usage(char *prog, char *usage_head,char usage_data[][80]);
}

void parse_args(void)
{
  int i, iret;
  char **word, **valu;
  char pfile[] = "satellite.params";

  iret = get_key_value_from_file(pfile, &word, &valu);
  if (iret == 0) usage("satellite_orbit");

				/* Set parameters */

  for (i=0; i<iret; i++)
    set_parm(word[i],valu[i]);

}


void set_parm(char *word, char *valu)
{
  if (!strcmp("HALO_MODEL",word))
    switch (atoi(valu)) {
    case file:
      HALO_MODEL = file;
      break;
    case isothermal:
      HALO_MODEL = isothermal;
      break;
    case hernquist_model:
      HALO_MODEL = hernquist_model;
      break;
    default:
      cerr << "No such HALO model type: " << (int)HALO_MODEL << endl;
      exit(-2);
    }

  else if (!strcmp("E",word))
    E = atof(valu);

  else if (!strcmp("K",word))
    K = atof(valu);

  else if (!strcmp("INCLINE",word))
    INCLINE = atof(valu);

  else if (!strcmp("PSI",word))
    PSI = atof(valu);

  else if (!strcmp("PHIP",word))
    PHIP = atof(valu);

  else if (!strcmp("VROT",word))
    VROT = atof(valu);

  else if (!strcmp("RCORE",word))
    RCORE = atof(valu);

  else if (!strcmp("RMODMIN",word))
    RMODMIN = atof(valu);

  else if (!strcmp("RMODMAX",word))
    RMODMAX = atof(valu);

  else if (!strcmp("RA",word))
    RA = atof(valu);

  else if (!strcmp("DIVERGE",word))
    DIVERGE = atoi(valu);

  else if (!strcmp("NUMDF",word))
    NUMDF = atoi(valu);

  else if (!strcmp("NRECS",word))
    NRECS = atoi(valu);

  else if (!strcmp("DIVERGE_RFAC",word))
    DIVERGE_RFAC = atof(valu);

  else if (!strcmp("MODFILE",word))
    MODFILE = valu;

  else {
    cerr << "No such paramter: " << word << " . . . quitting\n";
    exit(-1);
  }

}

#include <sys/types.h>
#include <sys/stat.h>

void usage(char *prog)
{
  char usage_head[] = 
    "[-f file -d] [keyword=value [keyword=value] .. ]";

  char usage_data[][80] = {
    "     -f file",      "keyword/value parameter file",
    "     -d",           "print default parameters",
    "\nKeywords:",       " INFILE,OUTFILE",
    "\0"                 };

  (void)std_usage(prog,usage_head,usage_data);
  exit(-1);
}


