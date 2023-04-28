/*****************************************************************************
 *  Description:
 *  -----------
 *
 *  Test orbit coordinate transform used in ResPot
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
 *  MDW 09/28/02
 *
 ***************************************************************************/

using namespace std;

#include <cstdlib>
#include <cmath>

#include <iostream>
#include <iomanip>
#include <fstream>
#include <sstream>
#include <memory>
#include <string>

#include <orbit.H>
#include <massmodel.H>
#include <biorth.H>
#include <gaussQ.H>

#include <biorth.H>
#include <sphereSL.H>
#include <model3d.H>
#include <isothermal.H>
#include <hernquist_model.H>
#include <plummer.H>

#include <BarForcing.H>
#include <ResPot.H>

#include <numerical.H>

int parse_args(int, char **);
void write_parm(void);
string respot_mpi_id() { return "testcoord.dbg"; }

				// Input parameters
int L1=-1;
int L2=2;
int L0=2;
int M0=2;
int LMAX=2;
int NUMX = 400;
int NUME = 200;
int RECS = 100;
int NUMR = 800;
int NMAX=10;
double rmin = 1.0e-3;
double rmax = 1.98;
double RCORE=1.0;
double EXPONENT=1.0;
double RSWITCH=0.1;
double RMODMIN=1.0e-3;
double RMODMAX=20.0;
double VROT=1.0;
double T0 = 1.0;
double TF = 1.0;
double TEND = 2.0;
double DELTA = 0.25;
double DOMEGA=0.0;
double AMP=1.0;
double DT = 0.0001;
double PHASE = 0.0;
Models3d HALO_MODEL=file;
int HALO_TYPE=sturm;
double MFRAC=0.1;
double LENGTH=0.067;
double COROT=1.0;
bool selfgrav = false;
bool logr = true;
string BODFILE = "bods";
string INFILE = "SLGridSph.model";
string OUTFILE = "coord.out";
string PARMFILE = "parm.file";

//

int myid = 0;
char threading_on = 0;
pthread_mutex_t mem_lock;

//

Eigen::Matrix3d return_euler(double PHI, double THETA, double PSI, int BODY);

AxiSymModPtr halo_model;

double pot(double r)
{
  return halo_model->get_pot(r);
}

double dpot(double r)
{
  return halo_model->get_dpot(r);
}

double dens(double r)
{
  return halo_model->get_density(r);
}

//

void dump_biorthSL(AxiSymBiorth* bio, double rmin, double rmax, 
		   int numr, int l, int nmax, const char* file)
{
  ofstream out(file);
  if (!out) {
    cerr << "dump_biorthSL: " << file << endl;
    return;
  }

  double r, dr = (rmax - rmin)/numr;
  for (int i=0; i<=numr; i++) {
    r = rmin + dr*i;
    out << setw(15) << r;
    for (int n=1; n<=nmax; n++) {
      out << setw(15) << bio->dens(n, l, bio->r_to_rb(r));
      out << setw(15) << bio->potl(n, l, bio->r_to_rb(r));
    }
    out << endl;
  }

}


double Ylm02(int ll, int mm)
{
  mm = abs(mm);
  return (2.0*ll+1)/(4.0*M_PI*M_PI) * exp(
         2.0*mm*log(2.0) + lgamma(1.0+ll-mm) - lgamma(1.0+ll+mm) +
         2.0*lgamma(0.5*(ll+mm+1)) - 2.0*lgamma(0.5*(ll-mm)+1.0) );
}


// Put x modulo 2*Pi
//
double wrap(double x)
{
  if (x<0.0)
    x +=  2.0*M_PI*(1.0 + int(-x/(2.0*M_PI)));
  else
    x += -2.0*M_PI*(0.0 + int( x/(2.0*M_PI)));

  return x;
}

int
main(int argc, char **argv)
{
  double	rmax;
  
// ===================================================================
// Parse command line
// ===================================================================

  if (!parse_args(argc,argv)) return -1;
  write_parm();

// ===================================================================
// Open output streams
// ===================================================================

  ofstream out(OUTFILE.c_str());
  if (!out) {
    cerr << "Couldn't open output file <" << OUTFILE << ">\n";
    exit(-1);
  }

				// Header
  out.setf(ios::left);
  out << "# " << setw(13) << "Number"	// 1
      << "| " << setw(13) << "E0"	// 2
      << "| " << setw(13) << "E1"	// 4
      << "| " << setw(13) << "K0"	// 3
      << "| " << setw(13) << "K1"	// 5
      << "| " << setw(13) << "I10"	// 6
      << "| " << setw(13) << "I1"	// 7
      << "| " << setw(13) << "I20"	// 8
      << "| " << setw(13) << "I2"	// 9
      << endl;

  out << "# " << setw(13) << 1;
  for (int i=2; i<=19; i++) out << "| " << setw(13) << i;
  out << endl;


// ===================================================================
// Initialize spherical model
// ===================================================================

  std::shared_ptr<SphericalModelTable> ml;

  switch (HALO_MODEL) {
  case file:
    ml =  std::make_shared<SphericalModelTable>(INFILE);
    ml->setup_df(400);			// isotropic
    halo_model = ml;
    Model3dNames[0] = INFILE;		// Assign filename to ID string
    rmax = halo_model->get_max_radius();
    break;

  case sing_isothermal:
    halo_model = std::shared_ptr<SingIsothermalSphere>(1.0, RMODMIN, RMODMAX);
    rmax = RMODMAX;
    break;

  case isothermal:
    halo_model = std::make_shared<IsothermalSphere>(RCORE, RMODMAX, VROT);
    rmax = RMODMAX;
    break;

  case hernquist_model:
    halo_model = std::make_shared<HernquistSphere>(1.0, RMODMIN, RMODMAX);
    rmax = RMODMAX;
    break; 

  case plummer:
    halo_model = std::make_shared<PlummerSphere>(1.0, RMODMIN, RMODMAX);
    rmax = RMODMAX;
    break; 

  default:
    cerr << "Illegal HALO model: " << HALO_MODEL << '\n';
    exit(-1);
  }


// -------------------------------------------------------------------
// Initilize biorthogonal functions (halo) 
// -------------------------------------------------------------------

  AxiSymBioPtr halo_ortho;
  switch (HALO_TYPE) {
  case bessel:
    halo_ortho = std::make_shared<BSSphere>(rmax, NMAX, (LMAX>0?LMAX:1) );
    break;
  case clutton_brock:
    halo_ortho = std::make_shared<CBSphere>();
    break;
  case hernquist:
    halo_ortho = std::make_shared<HQSphere>();
    break;
  case sturm:
    SphereSL::cache = 1;
    halo_ortho = std::make_shared<SphereSL>(LMAX, NMAX, NUMR, rmin, rmax, 1.0, ml);
    break;
  default:
    cerr << "Illegal spherical biorthongal series: " << HALO_TYPE << '\n';
    exit(-1);
  }


// -------------------------------------------------------------------
// Initialize bar model
// -------------------------------------------------------------------

  BarForcing::L0 = L0;
  BarForcing::M0 = M0;
  BarForcing::selfgrav = selfgrav;

  BarForcing bar(NMAX, MFRAC*halo_model->get_mass(LENGTH), LENGTH, COROT);
  bar.set_model(halo_model);
  bar.compute_quad_parameters();
  double OMEGA = bar.Omega();

// -------------------------------------------------------------------
// Initialize resonance potential
// -------------------------------------------------------------------

  cout << "Constructing ResPot . . . ";
  ResPot::NUMX = NUMX;
  ResPot::NUME = NUME;
  ResPot::RECS = RECS;
  ResPot respot(halo_model, &bar, L0, M0, L1, L2);
  cout << "done\n";


// -------------------------------------------------------------------
// Test against phase space
// -------------------------------------------------------------------

  ifstream bods(BODFILE.c_str());
  if (!bods) {
    cerr << "Couldn't open <" << BODFILE <<  "\n";
    exit(-2);
  }

  const int buffersize = 2048;
  char linebuf[buffersize];

  double posI[3], velI[3], posO[3], velO[3];
  double E0, K0, I10, I20, I1, I2, E1, K1, O1, O2;
  double W1, W2, W3, F, BETA, PSI, amp, Omega;
  int i, ret;

  SphericalOrbit orb(halo_model);

  if (!bods.getline(linebuf, buffersize)) {
    cerr << "Error reading first line <" << BODFILE <<  "\n";
    exit(-3);
  }

  istringstream ins(linebuf);

  int num, inum, dnum;

  ins >> num;
  ins >> inum;
  ins >> dnum;

  double mass;

  for (int n=0; n<num; n++) {
    
    if (!bods.getline(linebuf, buffersize)) {
      cerr << "Error reading ps line <" << BODFILE <<  "\n";
      exit(-4);
    }

    istringstream ins(linebuf);
    
    ins >> mass;
    for (int k=0; k<3; k++) ins >> posI[k];
    for (int k=0; k<3; k++) ins >> velI[k];
    
    respot.coord(posI, velI, E0, K0, I10, I20, O1, O2, 
		 W1, W2, W3, F, BETA, PSI);
    orb.new_orbit(E0, K0);
    I10 = orb.get_action(1);

    vector<double> phase(3);
    double Phase = PHASE;
    double T = 0.0;

    while (T<=TEND) {

      amp = AMP*0.25*(1.0 + erf((T - T0)/DELTA))*(1.0 + erf((TF - T)/DELTA));

      
      Omega = OMEGA*(1.0 + DOMEGA*(T-T0));
      phase[0] = PHASE;
      phase[1] = PHASE + Omega*DT*0.5;
      phase[2] = PHASE + Omega*DT*0.5;

      ret = respot.Update(DT, phase, amp, posI, velI, posO, velO);
      
      respot.coord(posO, velO, E1, K1, I1, I2, O1, O2, 
		   W1, W2, W3, F, BETA, PSI);
      
      for (int k=0; k<3; k++) {
	posI[k] = posO[k];
	velI[k] = velO[k];
      }
      
      T += DT;
      Phase += Omega*DT;
      i++;
    }
    
    out << setw(15) << n 
	<< setw(15) << E0
	<< setw(15) << E1
	<< setw(15) << K0
	<< setw(15) << K1
	<< setw(15) << I10
	<< setw(15) << I20
	<< setw(15) << I1
	<< setw(15) << I2
	<< endl;

    cout << setw(15) << n 
	 << setw(15) << (E1-E0)/E0
	 << setw(15) << (K1-K0)/K0
	 << setw(15) << (I10 - I1)/I10
	 << setw(15) << (I20 - I2)/I20
	 << endl;
  }

// -------------------------------------------------------------------
// Clean up
// -------------------------------------------------------------------

  cout << "Deleting halo_ortho . . . \n";
  delete halo_ortho;

  cout << "Deleting halo_model . . . \n";
  delete halo_model;

  return 0;
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

int parse_args(int argc, char **argv)
{
  char *prog=argv[0];
  int c, iparmf=0, iret, i;
  char file[128], **word, **valu;

  while (1) {

    c = getopt(argc, argv, "f:dh");
    if (c == -1) break;

    switch (c) {
    case 'f': 
      iparmf=1;
      strcpy(file, optarg);
      break;
    case 'd':
      print_default();
      break;
    case '?':
    case 'h':
      usage(prog);
      break;
    }
  }

  argc -= optind;
  if (iparmf)
    iret = get_key_value_from_file(file, &word, &valu);
  else
    {
      iret = get_key_value(argc,&argv[optind], &word, &valu);
      argc -= iret;
    }
  
    if (argc != 0)
      usage(prog);

				/* Set parameters */

  for (i=0; i<iret; i++)
    set_parm(word[i], valu[i]);

  return 1;
}

void set_parm(char *word, char *valu)
  {
  string wrd(word);

  if (wrd == "HALO_MODEL")
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
    case plummer:
      HALO_MODEL = plummer;
      break;
    default:
      cerr << "No such HALO model type: " << atoi(valu) << endl;
      exit(-2);
    }

  else if (wrd == "L1")
    L1 = atoi(valu);

  else if (wrd == "L2")
    L2 = atoi(valu);

  else if (wrd == "L0")
    L0 = atoi(valu);

  else if (wrd == "M0")
    M0 = atoi(valu);

  else if (wrd == "LMAX")
    LMAX = atoi(valu);

  else if (wrd == "NUMX")
    NUMX = atoi(valu);

  else if (wrd == "NUME")
    NUME = atoi(valu);

  else if (wrd == "RECS")
    RECS = atoi(valu);

  else if (wrd == "NUMR")
    NUMR = atoi(valu);

  else if (wrd == "NMAX")
    NMAX = atoi(valu);

  else if (wrd == "rmin")
    rmin = atof(valu);

  else if (wrd == "rmax")
    rmax = atof(valu);

  else if (wrd == "EXPONENT")
    EXPONENT = atof(valu);

  else if (wrd == "RSWITCH")
    RSWITCH = atof(valu);

  else if (wrd == "RCORE")
    RCORE = atof(valu);

  else if (wrd == "RMODMIN")
    RMODMIN = atof(valu);

  else if (wrd == "RMODMAX")
    RMODMAX = atof(valu);

  else if (wrd == "VROT")
    VROT = atof(valu);

  else if (wrd == "T0")
    T0 = atof(valu);

  else if (wrd == "TF")
    TF = atof(valu);

  else if (wrd == "TEND")
    TEND = atof(valu);

  else if (wrd == "DELTA")
    DELTA = atof(valu);

  else if (wrd == "AMP")
    AMP = atof(valu);

  else if (wrd == "DT")
    DT = atof(valu);

  else if (wrd == "PHASE")
    PHASE = atof(valu);

  else if (wrd == "MFRAC")
    MFRAC = atof(valu);

  else if (wrd == "LENGTH")
    LENGTH = atof(valu);

  else if (wrd == "COROT")
    COROT = atof(valu);

  else if (wrd == "HALO_TYPE")
    HALO_TYPE = atoi(valu);

  else if (wrd == "selfgrav")
    if (atoi(valu)) selfgrav = true;
    else selfgrav = false;

  else if (wrd == "logr")
    if (atoi(valu)) logr = true;
    else logr = false;

  else if (wrd == "BODFILE")
    BODFILE = valu;

  else if (wrd == "INFILE")
    INFILE = valu;

  else if (wrd == "OUTFILE")
    OUTFILE = valu;

  else if (wrd == "PARMFILE")
    PARMFILE = valu;

  else {
    cerr << "No such parameter: " << wrd.c_str() << endl;
    exit(-1);
  }

}

#include <sys/types.h>
#include <sys/stat.h>

void write_parm(void)
{
  ofstream fout(PARMFILE.c_str());			
  if ( !fout ) {
    cerr << "Couldn't open parameter file: " << PARMFILE 
	 << " . . . quitting\n";
    exit(-1);
  }

  print_parm(fout,"\0");
  fout.close();
}

void print_default()
{
  cerr << "\nDefaults:\n";
cerr << "----------------------------\n";
  print_parm(cerr,"\0");
  exit(0);
}


void print_parm(ostream& stream, char *comment)
{
  stream << comment << setw(15) << "L1" << " = " << L1 << endl;
  stream << comment << setw(15) << "L2" << " = " << L2 << endl;
  stream << comment << setw(15) << "L0" << " = " << L0 << endl;
  stream << comment << setw(15) << "M0" << " = " << M0 << endl;
  stream << comment << setw(15) << "LMAX" << " = " << LMAX << endl;
  stream << comment << setw(15) << "NUMX" << " = " << NUMX << endl;
  stream << comment << setw(15) << "NUME" << " = " << NUME << endl;
  stream << comment << setw(15) << "RECS" << " = " << RECS << endl;
  stream << comment << setw(15) << "NUMR" << " = " << NUMR << endl;
  stream << comment << setw(15) << "NMAX" << " = " << NMAX << endl;
  stream << comment << setw(15) << "rmin" << " = " << rmin << endl;
  stream << comment << setw(15) << "rmax" << " = " << rmax << endl;
  stream << comment << setw(15) << "EXPONENT" << " = " << EXPONENT << endl;
  stream << comment << setw(15) << "RSWITCH" << " = " << RSWITCH << endl;
  stream << comment << setw(15) << "RCORE" << " = " << RCORE << endl;
  stream << comment << setw(15) << "RMODMIN" << " = " << RMODMIN << endl;
  stream << comment << setw(15) << "RMODMAX" << " = " << RMODMAX << endl;
  stream << comment << setw(15) << "VROT" << " = " << VROT << endl;
  stream << comment << setw(15) << "T0" << " = " << T0 << endl;
  stream << comment << setw(15) << "TF" << " = " << TF << endl;
  stream << comment << setw(15) << "TEND" << " = " << TEND << endl;
  stream << comment << setw(15) << "DELTA" << " = " << DELTA << endl;
  stream << comment << setw(15) << "DT" << " = " << DT << endl;
  stream << comment << setw(15) << "PHASE" << " = " << PHASE << endl;
  stream << comment << setw(15) << "MFRAC" << " = " << MFRAC << endl;
  stream << comment << setw(15) << "LENGTH" << " = " << LENGTH << endl;
  stream << comment << setw(15) << "COROT" << " = " << COROT << endl;
  stream << comment << setw(15) << "selfgrav" << " = " << selfgrav << endl;
  stream << comment << setw(15) << "logr" << " = " << logr << endl;
  stream << comment << setw(15) << "HALO_MODEL" << " = " << Model3dNames[HALO_MODEL] << endl;
  stream << comment << setw(15) << "HALO_TYPE" << " = " << HALO_TYPE << endl;
  stream << comment << setw(15) << "BODFILE" << " = " << BODFILE << endl;
  stream << comment << setw(15) << "INFILE" << " = " << INFILE << endl;
  stream << comment << setw(15) << "OUTFILE" << " = " << OUTFILE << endl;
  stream << comment << setw(15) << "PARMFILE" << " = " << PARMFILE << endl;

}


void usage(char *prog)
{
  char usage_head[] = 
    "[-f file -d] [keyword=value [keyword=value] .. ]";

  char usage_data[][80] = {
    "     -f file",      "keyword/value parameter file",
    "     -d",           "print default parameters",
    "\nKeywords:",       " INFILE,OUTFILE",
    "\0"                 };



  // (void)std_usage(prog,usage_head,usage_data);
  exit(-1);
}


