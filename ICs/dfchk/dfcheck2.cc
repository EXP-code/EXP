// This may look like C code, but it is really -*- C++ -*-

/*****************************************************************************
 *  Description:
 *  -----------
 *
 *  Check DF computation for a massmodel class
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

#define IS_MAIN

#include <values.h>

#include <iostream.h>
#include <iomanip.h>
#include <fstream.h>
#include <GetOpt.h>
#include <string>
#include <vector>

#include <String.h>

#include <numerical.h>
#include <gaussQ.h>
#include <model3d.h>
#include <isothermal.h>
#include <hernquist.h>
#include <genpoly.h>

#include <fpetrap.h>

static int parse_args(int, char **);

				// Global 

int HMODEL = file;
int NUMR = 100;
int NUMDF = 200;
int NUMMS = 200;
int NUMMODEL = 500;
int RLOG = 0;
double NN = 2.5;
double MM = 0.5;
double RA = 1.0e8;
double RMODMIN = 1.0e-2;
double RMOD = 100.0;
double EPS = 1.0e-5;
int DIVERGE=0;
int LOGSCALE=0;
double DIVERGE_RFAC=1.0;
string INFILE = "infile";
string OUTFILE = "dfcheck2.dat";
string MODFILE = "newmodel.dat";
string parmfile = "parm.file";
bool NEWMODEL = false;

				// Global variables

static ScalarType type;

int
main(int argc, char **argv)
{

  set_fpu_handler();

  //
  // Parse command line
  //

  if (!parse_args(argc,argv)) return -1;


  // 
  // Prepare output streams
  //

  ofstream out(OUTFILE.c_str());
  if (!out) {
    cerr << "Couldn't open <" << OUTFILE << "> for output\n";
    exit (-1);
  }

  out.precision(6);
  out.setf(ios::scientific, ios::floatfield);

  ofstream out2;

  if (NEWMODEL) {
    out2.open(MODFILE.c_str());
    if (!out2) {
      cerr << "Couldn't open <" << MODFILE << "> for output\n";
      exit (-1);
    }

    out2.precision(6);
    out2.setf(ios::scientific, ios::floatfield);
  }

  cout.precision(6);
  cout.setf(ios::scientific, ios::floatfield);

  //
  // Begin integration
  //

  AxiSymModel *hmodel;
  SphericalModelTable *htmodel;
  SphericalModelTable::even = 0;
  SphericalModelTable::logscale = LOGSCALE;
  
  if (HMODEL>=0) {

    switch (HMODEL) {
    case file:
      htmodel = new SphericalModelTable(String(INFILE.c_str()), 
					DIVERGE, DIVERGE_RFAC);
      htmodel->setup_df(NUMDF, RA);
      hmodel = htmodel;
      break;
    case isothermal:
      hmodel = new IsothermalSphere(1.0, RMODMIN, RMOD); // Halo model
      break;
    case hernquist_model:
      hmodel = new HernquistSphere(1.0, RMODMIN, RMOD); // Halo model
      break; 
    case gen_polytrope:
      hmodel = new GeneralizedPolytrope(NUMMODEL, NN, MM, EPS);
      break;
    default:
      cerr << "No such HALO model type: " << (int)HMODEL << endl;
      exit(-2);
    }

  }


  out  << "# "
       << setw(13) << "Radius"
       << setw(15) << "Dens (orig)"
       << setw(15) << "Dens (DF)"
       << setw(15) << "Difference"
       << endl;

  out  << "# "
       << setw(13) << "------"
       << setw(15) << "-----------"
       << setw(15) << "---------"
       << setw(15) << "----------"
       << endl;


  cout  << setw(15) << "Radius"
	<< setw(15) << "Dens (orig)"
	<< setw(15) << "Dens (DF)"
	<< setw(15) << "Difference"
	<< endl;

  cout  << setw(15) << "------"
	<< setw(15) << "-----------"
	<< setw(15) << "---------"
	<< setw(15) << "----------"
	<< endl;

  LegeQuad lq(NUMDF);

  double den, vmax, r, pot, E, J, x, y;
  double rmin, rmax, dr;

  vector<double> rv, dv, mv, pw, p0;

  rv.reserve(NUMR);
  dv.reserve(NUMR);

  if (RLOG) {
    rmin = log(hmodel->get_min_radius());
    rmax = log(hmodel->get_max_radius());
  }
  else {
    rmin = hmodel->get_min_radius();
    rmax = hmodel->get_max_radius();
  }
  dr = (rmax - rmin)/NUMR;

  for (int i=1; i<=NUMR; i++) {
    if (RLOG)
      rv[i-1] = r = exp(rmin + dr*((double)i-0.5));
    else
      rv[i-1] = r = rmin + dr*((double)i-0.5);

    pot = hmodel->get_pot(r);
    vmax = sqrt(2.0*fabs(hmodel->get_pot(rmax) - pot));
    den=0.0;
    for (int ix=1; ix<=NUMDF; ix++) {
      x = lq.knot(ix);

      for (int iy=1; iy<=NUMDF; iy++) {
	y = lq.knot(iy);

	E = pot + 0.5*vmax*vmax*(x*x + (1.0-x*x)*y*y);
	J = vmax*sqrt(1.0 - x*x)*y*r;
	den += lq.weight(ix)*lq.weight(iy) * vmax*vmax*vmax * (1.0 - x*x)*y *
	  hmodel->distf(E, J);
      }
    }

    den *= 2.0*M_PI * 2.0;

    dv[i-1] = den;

    out  << setw(15) << r 
	 << setw(15) << hmodel->get_density(r) 
	 << setw(15) << den 
	 << setw(15) << den - hmodel->get_density(r) 
	 << endl;

    cout << setw(15) << r 
	 << setw(15) << hmodel->get_density(r) 
	 << setw(15) << den 
	 << setw(15) << den - hmodel->get_density(r) 
	 << endl;
  }

  LegeQuad lq2(NUMMS);

  double energy=0.0, mass=0.0;
  dr = rmax - rmin;

  for (int i=1; i<=NUMMS; i++) {

    r = rmin + dr*lq2.knot(i);

    pot = hmodel->get_pot(r);
    den = hmodel->get_density(r);

    mass += r*r*den * lq2.weight(i);
    energy += 0.5*r*r*den*pot * lq2.weight(i);
  }

  mass *= 4.0*M_PI*dr;
  energy *= 4.0*M_PI*dr;

  cout << "Mass=" << mass << endl;
  cout << "Energy=" << energy << endl;

  if (NEWMODEL) {

    mv.reserve(NUMR);
    pw.reserve(NUMR);
    p0.reserve(NUMR);

    mv[0] = 0.0;
    pw[0] = 0.0;

    for (int i=1; i<NUMR; i++) {
      mv[i] = mv[i-1] +
	2.0*M_PI*(rv[i-1]*rv[i-1]*dv[i-1] + rv[i]*rv[i]*dv[i])*(rv[i] - rv[i-1]);
      pw[i] = pw[i-1] +
	2.0*M_PI*(rv[i-1]*dv[i-1] + rv[i]*dv[i])*(rv[i] - rv[i-1]);
    }

    for (int i=0; i<NUMR; i++) 
      p0[i] = -mv[i]/(rv[i]+MINDOUBLE) - (pw[NUMR-1] - pw[i]);

  
    out2 << NUMR << endl;
    for (int i=0; i<NUMR; i++) 
      out2 << setw(15) << rv[i]
	   << setw(15) << dv[i]
	   << setw(15) << mv[i]
	   << setw(15) << p0[i]
	   << endl;
  }

  return 0;
}

//----------------------------------------------------------------------
//----------------------------------------------------------------------
//----------------------------------------------------------------------

static void usage(char *);
static void set_parm(char *, char *);
static void write_parm(void);
static void print_parm(ostream &, char *);
static void print_default(void);

extern "C" {
  int get_key_value(int, char **, char ***, char ***);
  int get_key_value_from_file(char *, char ***, char ***);
  void std_usage(char *prog, char *usage_head,char usage_data[][80]);
}

static int parse_args(int argc, char **argv)
{
  char c,*prog=argv[0];
  int iparmf=0,iret,i;
  char file[128],**word,**valu;

  // Invokes ctor `GetOpt (int argc, char **argv,
  //                       char *optstring);'
  GetOpt getopt(argc, argv, "f:dh");

  while ((c = getopt ()) != EOF)
      switch (c) {
      case 'f': 
	iparmf=1;
	strcpy(file, getopt.optarg);
	break;
      case 'd':
	print_default();
	break;
      case '?':
      case 'h':
	usage(prog);
	break;
      }

  argc -= getopt.optind;
  if (iparmf)
    iret = get_key_value_from_file(file,&word,&valu);
  else
    {
      iret = get_key_value(argc,&argv[getopt.optind],&word,&valu);
      argc -= iret;
    }
  
  if (argc != 0)
    usage(prog);

				/* Set parameters */

  for (i=0; i<iret; i++)
    set_parm(word[i],valu[i]);

  return 1;
}

static void set_parm(char *word, char *valu)
{
  if (!strcmp("NUMR",word))
    NUMR = atoi(valu);

  if (!strcmp("NUMDF",word))
    NUMDF = atoi(valu);

  if (!strcmp("NUMMS",word))
    NUMMS = atoi(valu);

  if (!strcmp("NUMMODEL",word))
    NUMMODEL = atoi(valu);

  if (!strcmp("NN",word))
    NN = atof(valu);

  if (!strcmp("MM",word))
    MM = atof(valu);

  if (!strcmp("RA",word))
    RA = atof(valu);

  if (!strcmp("RMOD",word))
    RMOD = atof(valu);

  if (!strcmp("EPS",word))
    EPS = atof(valu);

  if (!strcmp("RMODMIN",word))
    RMODMIN = atof(valu);

  if (!strcmp("LOGSCALE",word))
    LOGSCALE = atoi(valu);

  if (!strcmp("DIVERGE",word))
    DIVERGE = atoi(valu);

  if (!strcmp("DIVERGE_RFAC",word))
    DIVERGE_RFAC = atof(valu);

  if (!strcmp("RLOG",word))
    RLOG = atoi(valu);

  if (!strcmp("NEWMODEL",word))
    if (atoi(valu)) NEWMODEL = true;

  if (!strcmp("parmfile",word))
    parmfile = valu;

  if (!strcmp("INFILE",word))
    INFILE = valu;

  if (!strcmp("OUTFILE",word))
    OUTFILE = valu;

  if (!strcmp("HMODEL",word))
    switch (atoi(valu)) {
    case file:
      HMODEL = file;
      break;
    case isothermal:
      HMODEL = isothermal;
      break;
    case hernquist_model:
      HMODEL = hernquist_model;
      break;
    case gen_polytrope:
      HMODEL = gen_polytrope;
      break;
    default:
      cerr << "No such HALO model type: " << (int)HMODEL << endl;
      exit(-2);
    }
}

static void write_parm()
{
  ofstream fout(parmfile.c_str());			
  if ( !fout ) {
    cerr.form("Couldn't open parameter file: %s . . . quitting\n",
	      parmfile.c_str());
    exit(-1);
  }

  print_parm(fout,"\0");
}

static void print_default()
{
  cerr.form("\nDefaults:\n");
  cerr.form("----------------------------\n");
  print_parm(cerr,"\0");
  exit(0);
}


static void print_parm(ostream& stream, char *comment)
{
  stream.form("%s%-15s = %d\n",comment,"NUMR", NUMR);
  stream.form("%s%-15s = %d\n",comment,"NUMDF", NUMDF);
  stream.form("%s%-15s = %d\n",comment,"NUMMS", NUMMS);
  stream.form("%s%-15s = %d\n",comment,"NUMMODEL", NUMMODEL);
  stream.form("%s%-15s = %f\n",comment,"NN", NN);
  stream.form("%s%-15s = %f\n",comment,"MM", MM);
  stream.form("%s%-15s = %e\n",comment,"RA", RA);
  stream.form("%s%-15s = %e\n",comment,"RMOD", RMOD);
  stream.form("%s%-15s = %e\n",comment,"EPS", EPS);
  stream.form("%s%-15s = %d\n",comment,"LOGSCALE", LOGSCALE);
  stream.form("%s%-15s = %d\n",comment,"DIVERGE", DIVERGE);
  stream.form("%s%-15s = %d\n",comment,"DIVERGE_RFAC", DIVERGE_RFAC);
  stream.form("%s%-15s = %e\n",comment,"RMODMIN", RMODMIN);
  if (NEWMODEL)
    stream.form("%s%-15s = %s\n",comment,"NEWMODEL", "true");
  else
    stream.form("%s%-15s = %s\n",comment,"NEWMODEL", "false");
  stream.form("%s%-15s = %d\n",comment,"RLOG", RLOG);
  stream.form("%s%-15s = %s\n",comment,"HMODEL", 
	      (const char *)Model3dNames[HMODEL]);
  stream.form("%s%-15s = %s\n",comment,"parmfile", parmfile.c_str());
  stream.form("%s%-15s = %s\n",comment,"INFILE", INFILE.c_str());
  stream.form("%s%-15s = %s\n",comment,"OUTFILE", OUTFILE.c_str());
}


static void usage(char *prog)
{
  char usage_head[] = 
    "[-f file -d] [keyword=value [keyword=value] .. ]";

  char usage_data[][80] = {
    "     -f file",      "keyword/value parameter file",
    "     -d",           "print default parameters",
    "\nKeywords:",       " OUT,INFILE,OUTFILE",
    "\0"                 };



  (void)std_usage(prog,usage_head,usage_data);
  exit(-1);
}


