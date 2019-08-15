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


#include <iostream>
#include <iomanip>
#include <string>
#include <cmath>

#include <GetOpt.h>

#include <numerical.h>
#include <Vector.h>

using namespace std;

int parse_args(int, char **);


string parmfile = "parm.file";
string INFILE = "infile";
string OUTFILE = "infile.new";

int
main(int argc, char **argv)
{

  //
  // Parse command line
  //

  if (!parse_args(argc,argv)) return -1;

  //
  // Begin integration
  //

  Vector bessjz(int n, int m);

  int n, m;
  cout << "n, m: ";
  cin >> n;
  cin >> m;

  Vector ans = bessjz(n, m);

  cout.precision(8);

  for (int i=1; i<=m; i++)
    cout << setw(4) << i << ">  " << setw(14) << ans[i] 
      << setw(14) << jn(n, ans[i]) << endl;

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
  if (!strcmp("INFILE",word))
    INFILE = valu;
  
  if (!strcmp("OUTFILE",word))
    OUTFILE = valu;
}

void write_parm()
{
  ofstream fout(parmfile.c_str());			
  if ( !fout ) {
    cerr.form("Couldn't open parameter file: %s . . . quitting\n",
	      (const char *)parmfile.c_str());
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
