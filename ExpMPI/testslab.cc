// This is really C++ code
/*****************************************************************************
 *  Description:
 *  -----------
 *
 *  This routine computes the potential, acceleration and density using
 *  periodic box expansion in X & Y and slab in Z
 *
 *
 *  Call sequence:
 *  -------------
 *
 *  Parameters:
 *  ----------
 *
 *  Returns:
 *  -------
 *
 *  Value
 *
 *  Notes:
 *  -----
 *  G=1
 *
 *  By:
 *  --
 *
 *  MDW 04/27/92; revision 04/29/92
 *
 *
 ***************************************************************************/
#define IS_MAIN
#define CPLUSPLUS

#include <stdlib.h>
#include <iostream.h>
#include <iomanip.h>
#include <fstream.h>
#include <math.h>
#include <ACG.h>
#include <Uniform.h>
#include <GetOpt.h>

#include <kevin_complex.h>
#include <biorth1d.h>
#include <logic.h>


int parse_args(int, char **);

				// Parameters

double HSCALE=0.01;
double ZMAX=0.1;
double KEPS=1.0e-6;
int NUMBER=10000;
int NMAXX=10;
int NMAXY=10;
int NMAXZ=10;
int NMINX=0;
int NMINY=0;
int OUTX=10;
int OUTZ=10;
Logic CHECK=FALSE;
String OUTFILE="stest";
String parmfile="stest.parm";

				// Global variables

static Complex *expccof;
static double *expreal, *expimag;
static double *expreal1, *expimag1;
static OneDTrig **trig;

int
main(int argc, char** argv)
{
  //
  // Parse command line
  //

  if (!parse_args(argc,argv)) return -1;

  //======================================================================

  static int firstime=1;
  static int imx,imy,imz,jmax,nnmax;
  int i, j, ix, iy, iz, iix, iiy, ii, jj, indx, dum;
  static double  dfac;
  static Complex kfac;
  Complex fac, startx, starty, startz, facx, facy, facz, potl, facf, facd;
  Complex accx, accy, accz, dens;
  Complex accxp, accyp, acczp, densp, potlp;
  Complex stepx, stepy, stepz;
  double kk, x, y, z;
  static Vector zfrc, zpot, zden;

  if (firstime) {
    firstime=0;

    imx = 1+2*NMAXX;
    imy = 1+2*NMAXY;
    imz = NMAXZ;
    jmax = imx*imy*imz;
    expccof = new Complex[jmax];

    expreal = new double [jmax];
    expreal1 = new double [jmax];
    expimag = new double [jmax];
    expimag1 = new double [jmax];
    
    dfac = 2.0*M_PI;
    kfac = Complex(0.0, dfac);

    nnmax = (NMAXX > NMAXY) ? NMAXX : NMAXY;

    trig = new OneDTrig* [nnmax+1];
    for (i=0; i<=nnmax; i++) {
      trig[i] = new OneDTrig [i+1];
      for (j=0; j<=i; j++) trig[i][j].reset(dfac*sqrt(i*i + j*j)+KEPS, ZMAX);
    }

    zpot.setsize(1, NMAXZ);
    zden.setsize(1, NMAXZ);
    zfrc.setsize(1, NMAXZ);
  }


  /*=====================*/
  /* Check orthogonality */
  /*=====================*/

  if (CHECK) {

    int NINT=400;
    double dz = 2.0*ZMAX/NINT;
    Matrix test(1, NMAXZ, 1, NMAXZ);

    while(1) {

      cout << "iix iiy: ";
      cin >> iix;
      cin >> iiy;


      if (iix < 0) break;
      if (iix>=nnmax) iix = nnmax;
      if (iiy>=nnmax) iiy = nnmax;

      if (iix<iiy) {
	int tmp = iix;
	iix = iiy;
	iiy = tmp;
      }

      test.zero();
      for (i=0; i<NINT; i++) {

	z = -ZMAX + (0.5+i)*dz;
	trig[iix][iiy].potl(dum, dum, z, zpot);
	trig[iix][iiy].dens(dum, dum, z, zden);
	
	for (ii=1; ii<=NMAXZ; ii++) {
	  for (jj=1; jj<=NMAXZ; jj++) test[ii][jj] += dz*zpot[ii]*zden[jj];
	}
      }

				// Test it!

      double min1=1e20, max1=-1e20, min0=1e20, max0=-1e20;
      for (ii=1; ii<=NMAXZ; ii++) {
	for (jj=1; jj<=NMAXZ; jj++) {
	  if (ii==jj) {
	    if (test[ii][jj] < min1) min1 = test[ii][jj];
	    if (test[ii][jj] > max1) max1 = test[ii][jj];
	  } else {
	    if (test[ii][jj] < min0) min0 = test[ii][jj];
	    if (test[ii][jj] > max0) max0 = test[ii][jj];
	  }
	}
      }
      
      cout.form("Off diag: min=%f  max=%f\n", min0, max0);
      cout.form("Diagonal: min=%f  max=%f\n", min1, max1);
    }

  }

  //
  // Check ok! 1/23/98
  //


  /*======================*/
  /* Compute coefficients */
  /*======================*/

				//  Coefficients are ordered as follows:
				//  n=-nmax,-nmax+1,...,0,...,nmax-1,nmax
				//  in a single array for each dimension
				//  with z dimension changing most rapidly
  /*		Clean */
  for (indx=0; indx<jmax; indx++) {
    expccof[indx] = 0.0;
    expreal[indx] = 0.0;
    expimag[indx] = 0.0;
  }

  
  /*            Compute coefficients */


  ACG gen(10, 20);
  Uniform Unit(0.0, 1.0, &gen);
  double mass = 1.0/NUMBER;
  

  for (i=0; i<NUMBER; i++) {

    double x = Unit();
    double y = Unit();
    double z = HSCALE*atanh(2.0*Unit() - 1.0);


				/* Recursion multipliers */
    stepx = exp(-kfac*x);
    stepy = exp(-kfac*y);

				/* Initial values */
    startx = exp(NMAXX*kfac*x);
    starty = exp(NMAXY*kfac*y);

    for (facx=startx, ix=0; ix<imx; ix++, facx*=stepx) {

      iix = NMAXX - ix;
      if (iix<0) iix = ix - NMAXX;

      for (facy=starty, iy=0; iy<imy; iy++, facy*=stepy) {
	  
	iiy = NMAXY - iy;
	if (iiy<0) iiy = iy - NMAXY;

	if (iix < 0 || iix > NMAXX) {
	  cerr << "Out of bounds: iix=" << iix << endl;
	}
	if (iiy < 0 || iiy > NMAXY) {
	  cerr << "Out of bounds: iiy=" << iiy << endl;
	}

	if (iix>=iiy)
	  trig[iix][iiy].potl(dum, dum, z, zpot);
	else
	  trig[iiy][iix].potl(dum, dum, z, zpot);

	for (iz=0; iz<imz; iz++) {

	  indx = imz*(iy + imy*ix) + iz;

	  expccof[indx] += 4.0*M_PI*mass*facx*facy*zpot[iz+1];
	  // ^
	  // |--- density in orthogonal series
	  //      is 4.0*M_PI rho
	}
      }
    }

  }



  /* Determine potential, acceleration, density */

  ofstream **out = new ofstream* [6];
  String name;

  name = OUTFILE + ".potl";
  out[0] = new ofstream (name);
  name = OUTFILE + ".dens";
  out[1] = new ofstream (name);
  name = OUTFILE + ".force";
  out[2] = new ofstream (name);

  name = OUTFILE + ".potl.proj";
  out[3] = new ofstream (name);
  name = OUTFILE + ".dens.proj";
  out[4] = new ofstream (name);
  name = OUTFILE + ".force.proj";
  out[5] = new ofstream (name);

  float f;

  for (i=0; i<6; i++) {
    out[i]->write(&OUTX, sizeof(int));
    out[i]->write(&OUTZ, sizeof(int));
    out[i]->write(&(f=0.0), sizeof(float));
    out[i]->write(&(f=1.0), sizeof(float));
    out[i]->write(&(f=-ZMAX), sizeof(float));
    out[i]->write(&(f= ZMAX), sizeof(float));
  }

  double dx = 1.0/(OUTX-1);
  double dz = 2.0*ZMAX/(OUTZ-1);

  y = 0.0;

  for (int zz=0; zz<OUTZ; zz++) {

    z = -ZMAX + dz*zz;

    for (int xx=0; xx<OUTX; xx++) {

      x = dx*xx;

      accz = dens = potl = 0.0;
      acczp = densp = potlp = 0.0;

				/* Recursion multipliers */
      stepx = exp(kfac*x);
      stepy = exp(kfac*y);

				/* Initial values (note sign change) */
      startx = exp(-NMAXX*kfac*x);
      starty = exp(-NMAXY*kfac*y);

      for (facx=startx, ix=0; ix<imx; ix++, facx*=stepx) {

				/* Compute wavenumber; recall that the */
				/* coefficients are stored as follows: */
				/* -nmax,-nmax+1,...,0,...,nmax-1,nmax */
	ii = ix - NMAXX;
	iix = abs(ii);

	for (facy=starty, iy=0; iy<imy; iy++, facy*=stepy) {

	  jj = iy - NMAXY;
	  iiy = abs(jj);

	  if (iix < 0 || iix > NMAXX) {
	    cerr << "Out of bounds: iix=" << iix << endl;
	  }
	  if (iiy < 0 || iiy > NMAXY) {
	    cerr << "Out of bounds: iiy=" << iiy << endl;
	  }

	  if (iix>=iiy) {
	    trig[iix][iiy].potl(dum, dum, z, zpot);
	    trig[iix][iiy].dens(dum, dum, z, zden);
	    trig[iix][iiy].force(dum, dum, z, zfrc);
	  }
	  else {
	    trig[iiy][iix].potl(dum, dum, z, zpot);
	    trig[iiy][iix].dens(dum, dum, z, zden);
	    trig[iiy][iix].force(dum, dum, z, zfrc);
	  }

	  for (iz=0; iz<imz; iz++) {

	    indx = imz*(iy + imy*ix) + iz;


	    fac  = facx*facy*zpot[iz+1]*expccof[indx];
	    facd = facx*facy*zden[iz+1]*expccof[indx];
	    facf = facx*facy*zfrc[iz+1]*expccof[indx];


				/* Don't add potential for 
				   zero wavenumber (constant term) */

	    if (ii==0 && jj==0 && iz==0) fac = 0.0;


				/* Limit to minimum wave number */

	    if (abs(ii)<NMINX || abs(jj)<NMINY) continue;

	    potl -= fac;
	    dens += facd;
	    
	    accx -= Complex(0.0,-dfac*ii)*fac;
	    accy -= Complex(0.0,-dfac*jj)*fac;
	    accz -= facf;


	    if (jj==0) {
	      potlp -= fac;
	      densp += facd;
	    
	      accxp -= Complex(0.0,-dfac*ii)*fac;
	      accyp -= Complex(0.0,-dfac*jj)*fac;
	      acczp -= facf;
	    }

	  }
	}
      }

      dens  /= 4.0*M_PI;
      densp /= 4.0*M_PI;

      out[0]->write(&(f=Re(potl)), sizeof(float));
      out[1]->write(&(f=Re(dens)), sizeof(float));
      out[2]->write(&(f=Re(accz)), sizeof(float));

      out[3]->write(&(f=Re(potlp)), sizeof(float));
      out[4]->write(&(f=Re(densp)), sizeof(float));
      out[5]->write(&(f=Re(acczp)), sizeof(float));

    }

  }

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

void set_parm(char *word, char *valu)
{
  if (!strcmp("HSCALE",word))
    HSCALE = atof(valu);

  if (!strcmp("ZMAX",word))
    ZMAX = atof(valu);

  if (!strcmp("KEPS",word))
    KEPS = atof(valu);

  if (!strcmp("NUMBER",word))
    NUMBER = atoi(valu);

  if (!strcmp("NMAXX",word))
    NMAXX = atoi(valu);

  if (!strcmp("NMAXY",word))
    NMAXY = atoi(valu);

  if (!strcmp("NMAXZ",word))
    NMAXZ = atoi(valu);

  if (!strcmp("NMINX",word))
    NMINX = atoi(valu);

  if (!strcmp("NMINY",word))
    NMINY = atoi(valu);

  if (!strcmp("OUTX",word))
    OUTX = atoi(valu);

  if (!strcmp("OUTZ",word))
    OUTZ = atoi(valu);

  if (!strcmp("CHECK",word))
    CHECK = ITOL(atoi(valu));

  if (!strcmp("OUTFILE",word))
    OUTFILE = valu;

}

void write_parm()
{
  ofstream fout(parmfile);			
  if ( !fout ) {
    cerr.form("Couldn't open parameter file: %s . . . quitting\n",
	      (const char *)parmfile);
    exit(-1);
  }

  print_parm(fout,"\0");
}

void print_default()
{
  cerr.form("\nDefaults:\n");
  cerr.form("----------------------------\n");
  print_parm(cerr,"\0");
  exit(0);
}


void print_parm(ostream& stream, char *comment)
{
  stream.form("%s%-15s = %e\n",comment,"HSCALE", HSCALE);
  stream.form("%s%-15s = %e\n",comment,"ZMAX", ZMAX);
  stream.form("%s%-15s = %e\n",comment,"KEPS", KEPS);
  stream.form("%s%-15s = %d\n",comment,"NUMBER", NUMBER);
  stream.form("%s%-15s = %d\n",comment,"NMAXX", NMAXX);
  stream.form("%s%-15s = %d\n",comment,"NMAXY", NMAXY);
  stream.form("%s%-15s = %d\n",comment,"NMAXZ", NMAXZ);
  stream.form("%s%-15s = %d\n",comment,"NMINX", NMINX);
  stream.form("%s%-15s = %d\n",comment,"NMINY", NMINY);
  stream.form("%s%-15s = %d\n",comment,"OUTX", OUTX);
  stream.form("%s%-15s = %d\n",comment,"OUTZ", OUTZ);
  stream.form("%s%-15s = %s\n",comment,"CHECK", (const char *)SLogic[CHECK]);
  stream.form("%s%-15s = %s\n",comment,"OUTFILE", (const char *)OUTFILE);
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



  (void)std_usage(prog,usage_head,usage_data);
  exit(-1);
}
