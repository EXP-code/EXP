// This may look like C code, but it is really -*- C++ -*-

/*****************************************************************************
 *  Description:
 *  -----------
 *
 *  C++ program header template
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
 *  updated 6/10/94
 *
 ***************************************************************************/


#include <iostream.h>
#include <fstream.h>
#include <iomanip.h>
#include <GetOpt.h>
#include <String.h>
#include <numerical.h>
#include <math.h>

#include <vector.h>
#include <gaussQ.h>

int parse_args(int, char **);
void get_potl_CBDisk(int lmax, int nmax, double r, Matrix& p);
void get_dens_CBDisk(int lmax, int nmax, double r, Matrix& d);
double normCBDisk(int n, int m);


int NUM=100;
int NMAX=5;
int L=2;
double RMIN=0.0;
double RMAX=10.0;
String parmfile = "parm.file";
String INFILE = "infile";
String OUTFILE = "infile.new";

Matrix work;

main(int argc, char **argv)
{

  //
  // Parse command line
  //

  if (!parse_args(argc,argv)) return -1;

  //
  // Begin integration
  //

  LegeQuad gw(NUM);

  Matrix test(1, NMAX, 1, NMAX);
  test.zero();
  Matrix p, d;
  work.setsize(0, L+1, 1, NMAX);
  double r;

  for (int i=1; i<=NUM; i++) {
    r = RMIN + gw.knot(i)*(RMAX - RMIN);

    get_potl_CBDisk(L, NMAX, r, p);
    get_dens_CBDisk(L, NMAX, r, d);

    for (int j=1; j<=NMAX; j++) {
      for (int k=1; k<=NMAX; k++) {
	test[j][k] += r*p[L][j]*d[L][k] * gw.weight(i);
      }
    }
  }
  
  double fac = 2.0*M_PI*(RMAX - RMIN);

  cout.setf(ios::scientific, ios::floatfield);
  cout.precision(4);

  for (int j=1; j<=NMAX; j++) {
    for (int k=1; k<=NMAX; k++)
      cout << setw(13) << test[j][k]*fac/normCBDisk(j-1, L);
    cout << endl;
  }

}

void get_potl_CBDisk(int lmax, int nmax, double r, Matrix& p)
{
  double r2 = r*r;
  double fac = 1.0/(1.0 + r2);
  double fac1 = (r2 - 1.0)*fac;
  double cur, curl1, curl2, cur0 = sqrt(fac), rcum = 1.0;
  int l, nn;

  p.setsize(0, lmax+1, 0, nmax);

  for (l=0; l<=lmax; l++) {
    cur = cur0;

    p[l][1] = cur*rcum;
    curl1 = 0.0;

    for (nn=1; nn<nmax; nn++) {
      curl2 = curl1;
      curl1 = cur;
      cur = (2.0 + (double)(2*l-1)/nn)*fac1*curl1 - 
	(1.0 + (double)(2*l-1)/nn)*curl2;
      p[l][nn+1] = cur*rcum;
    }
    cur0 *= fac*(2*(l+1) - 1);
    rcum *= r;
  }

}


void get_dens_CBDisk(int lmax, int nmax, double r, Matrix& d)
{

  double r2 = r*r;
  double fac = 1.0/(1.0 + r2);
  double fac1 = (r2 - 1.0)*fac;
  double cur, curl1, curl2, cur0 = sqrt(fac), rcum = 1.0;
  int l, nn;
  double fac0 = 1.0/(2.0*M_PI);

  d.setsize(0, lmax+1, 1, nmax);

  for (l=0; l<=lmax+1; l++) {
    cur = cur0;

    work[l][1] = cur;
    curl1 = 0.0;

    for (nn=1; nn<nmax; nn++) {
      curl2 = curl1;
      curl1 = cur;
      cur = (2.0 + (double)(2*l-1)/nn)*fac1*curl1 - 
	(1.0 + (double)(2*l-1)/nn)*curl2;
      work[l][nn+1] = cur;
    }
    cur0 *= fac*(2*(l+1) - 1);
  }

  for (l=0; l<=lmax; l++) {
    d[l][1] = work[l+1][1]*rcum*fac0;
    d[l][2] = work[l+1][2]*rcum*fac0;
    for (nn=2; nn<nmax; nn++)
      d[l][nn+1] = (work[l+1][nn+1] - work[l+1][nn-1])*rcum*fac0;

    rcum *= r;
  }

}

double normCBDisk(int n, int m)
{
  double ans = 1.0;
  int i;
 
  for (i=n+1; i<=n+2*m; i++)
    ans *= i;

  return pow(0.5, 2*m+1)*ans;
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
  if (!strcmp("NUM",word))
    NUM = atoi(valu);

  if (!strcmp("NMAX",word))
    NMAX = atoi(valu);

  if (!strcmp("L",word))
    L = atoi(valu);

  if (!strcmp("RMIN",word))
    RMIN = atof(valu);

  if (!strcmp("RMAX",word))
    RMAX = atof(valu);

  if (!strcmp("parmfile",word))
    parmfile = valu;

  if (!strcmp("INFILE",word))
    INFILE = valu;

  if (!strcmp("OUTFILE",word))
    OUTFILE = valu;

}

void write_parm()
{
  ofstream fout(parmfile);			
  if ( !fout ) {
    cerr.form("Couldn't open parameter file: %s . . . quitting\n",
	      (char *)parmfile);
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
  stream.form("%s%-15s = %d\n",comment,"NUM", NUM);
  stream.form("%s%-15s = %d\n",comment,"NMAX", NMAX);
  stream.form("%s%-15s = %d\n",comment,"L", L);
  stream.form("%s%-15s = %e\n",comment,"RMIN", RMIN);
  stream.form("%s%-15s = %e\n",comment,"RMAX", RMAX);
  stream.form("%s%-15s = %s\n",comment,"parmfile", (char *)parmfile);
  stream.form("%s%-15s = %s\n",comment,"INFILE", (char *)INFILE);
  stream.form("%s%-15s = %s\n",comment,"OUTFILE", (char *)OUTFILE);
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
