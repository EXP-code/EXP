// This may look like C code, but it is really -*- C++ -*-

/*****************************************************************************
 *  Description:
 *  -----------
 *
 *  Generate a phase space realization from a distribution function 
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

#include <iostream.h>
#include <iomanip.h>
#include <fstream.h>
#include <GetOpt.h>
#include <String.h>

#include <numerical.h>
#include <model3d.h>
#include <isothermal.h>
#include <hernquist.h>
#include <genpoly.h>


static int parse_args(int, char **);

extern "C" {
  void fpeinit(int);
}

				// Global 

int HMODEL = file;
int N = 10;
int NUMDF = 800;
int NUMR = 400;
int NUMJ = 400;
int NUMG = 800;
int NREPORT = 1000;
int SEED = 11;
int ITMAX = 10000;
int NUMMODEL = 500;
int DIVERGE = 0;
double DIVERGE_RFAC = 1.5;
int LOGSCALE = 0;
double NN = 2.5;
double MM= 0.5;
double RA = 1.0e8;
double RMODMIN = 1.0e-2;
double RMOD = 100.0;
double EPS = 1.0e-5;
String INFILE = "infile";
String OUTFILE = "phase_space.new";
String parmfile = "parm.file";

				// Global variables

static ScalarType type;

int
main(int argc, char **argv)
{

#ifdef FPETRAP  
  fpeinit(0);
#endif

  //
  // Parse command line
  //

  if (!parse_args(argc,argv)) return -1;


  // 
  // Prepare output streams
  //

  ofstream out(OUTFILE);
  if (!out) {
    cerr << "Couldn't open <" << OUTFILE << "> for output\n";
    exit (-1);
  }

  out.precision(6);
  out.setf(ios::scientific, ios::floatfield);

  //
  // Begin integration
  //

  AxiSymModel *hmodel;
  SphericalModelTable *htmodel;
  SphericalModelTable::even = 0;
  SphericalModelTable::logscale = 1;
  
  
  AxiSymModel::numr = NUMR;
  AxiSymModel::numj = NUMJ;
  AxiSymModel::gen_N = NUMG;

  if (HMODEL>=0) {

				// Halo model
    switch (HMODEL) {
    case file:
      SphericalModelTable::even = 0;
      SphericalModelTable::logscale = LOGSCALE;
      htmodel = new SphericalModelTable(INFILE, DIVERGE, DIVERGE_RFAC);
      htmodel->setup_df(NUMDF, RA);
      hmodel = htmodel;
      break;
    case isothermal:
      hmodel = new IsothermalSphere(1.0, RMODMIN, RMOD);
      break;
    case hernquist_model:
      hmodel = new HernquistSphere(1.0, RMODMIN, RMOD);
      break; 
    case gen_polytrope:
      hmodel = new GeneralizedPolytrope(NUMMODEL, NN, MM, EPS);
      break;
    default:
      cerr << "No such HALO model type: " << (int)HMODEL << endl;
      exit(-2);
    }

  }

  hmodel->set_seed(SEED);
  hmodel->set_itmax(ITMAX);

  Vector ps;
  int ierr;
  double mass = hmodel->get_mass(hmodel->get_max_radius())/N;
  out << setw(8) << N << setw(15) << 0.0 << endl;
  for (int n=0; n<N; n++) {
    
    do {
      ps = hmodel->gen_point(ierr);
    } while (ierr);

    out << setw(15) << mass;
    for (int i=1; i<=6; i++) out << setw(15) << ps[i];
    out << endl;

    if (!((n+1)%NREPORT)) cout << "." << n+1 << flush;

  }
  cout << endl;

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
  if (!strcmp("NUMDF",word))
    NUMDF = atoi(valu);

  if (!strcmp("NUMR",word))
    NUMR = atoi(valu);

  if (!strcmp("NUMJ",word))
    NUMJ = atoi(valu);

  if (!strcmp("NUMG",word))
    NUMG = atoi(valu);

  if (!strcmp("RA",word))
    RA = atof(valu);

  if (!strcmp("RMOD",word))
    RMOD = atof(valu);

  if (!strcmp("RMODMIN",word))
    RMODMIN = atof(valu);

  if (!strcmp("EPS",word))
    EPS = atof(valu);

  if (!strcmp("N",word))
    N = atoi(valu);

  if (!strcmp("SEED",word))
    SEED = atoi(valu);

  if (!strcmp("ITMAX",word))
    ITMAX = atoi(valu);

  if (!strcmp("NUMMODEL",word))
    NUMMODEL = atoi(valu);

  if (!strcmp("NN",word))
    NN = atof(valu);

  if (!strcmp("MM",word))
    MM = atof(valu);

  if (!strcmp("NREPORT",word))
    NREPORT = atoi(valu);

  if (!strcmp("LOGSCALE",word))
    LOGSCALE = atoi(valu);

  if (!strcmp("DIVERGE",word))
    DIVERGE = atoi(valu);

  if (!strcmp("DIVERGE_RFAC",word))
    DIVERGE_RFAC = atof(valu);

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
  ofstream fout(parmfile);			
  if ( !fout ) {
    cerr.form("Couldn't open parameter file: %s . . . quitting\n",
	      (const char *)parmfile);
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
  stream.form("%s%-15s = %d\n",comment,"NUMDF", NUMDF);
  stream.form("%s%-15s = %d\n",comment,"NUMR", NUMR);
  stream.form("%s%-15s = %d\n",comment,"NUMJ", NUMJ);
  stream.form("%s%-15s = %d\n",comment,"NUMG", NUMG);
  stream.form("%s%-15s = %e\n",comment,"RA", RA);
  stream.form("%s%-15s = %e\n",comment,"RMOD", RMOD);
  stream.form("%s%-15s = %e\n",comment,"RMODMIN", RMODMIN);
  stream.form("%s%-15s = %e\n",comment,"EPS", EPS);
  stream.form("%s%-15s = %d\n",comment,"N", N);
  stream.form("%s%-15s = %d\n",comment,"SEED", SEED);
  stream.form("%s%-15s = %d\n",comment,"ITMAX", ITMAX);
  stream.form("%s%-15s = %d\n",comment,"NUMMODEL", NUMMODEL);
  stream.form("%s%-15s = %f\n",comment,"NN", NN);
  stream.form("%s%-15s = %f\n",comment,"MM", MM);
  stream.form("%s%-15s = %d\n",comment,"NREPORT", NREPORT);
  stream.form("%s%-15s = %d\n",comment,"LOGSCALE", LOGSCALE);
  stream.form("%s%-15s = %d\n",comment,"DIVERGE", DIVERGE);
  stream.form("%s%-15s = %f\n",comment,"DIVERGE_RFAC", DIVERGE_RFAC);
  stream.form("%s%-15s = %s\n",comment,"HMODEL", (const char *)Model3dNames[HMODEL]);
  stream.form("%s%-15s = %s\n",comment,"parmfile", (const char *)parmfile);
  stream.form("%s%-15s = %s\n",comment,"INFILE", (const char *)INFILE);
  stream.form("%s%-15s = %s\n",comment,"OUTFILE", (const char *)OUTFILE);
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


