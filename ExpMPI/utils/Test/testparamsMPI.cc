/*****************************************************************************
 *  Description:
 *  -----------
 *
 *  Test parameter parsing
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
 *  MDW 02/05/04
 *
 ***************************************************************************/

#include <cstdlib>
#include <iostream>
#include <iomanip>
#include <fstream>
#include <sstream>
#include <string>
#include <cstring>
#include <vector>

#include <cstdlib>
#include <cmath>

#include <ParamParseMPI.H>

int myid=0;
char threading_on = 0;
pthread_mutex_t mem_lock;

using namespace std;

int parse_args(int argc, char **argv);
void print_parm(ostream &, const char *);


				// Parameters
int NICE=15;
bool DENS=true;
double RMAX=2.0;
string PARMFILE="test.param";

int main(int argc, char **argv)
{
  /*===================*/
  /* MPI preliminaries */
  /*===================*/

  int numprocs, slaves, myid, proc_namelen;
  char processor_name[MPI_MAX_PROCESSOR_NAME];

  MPI_Init(&argc, &argv);
  MPI_Comm_size(MPI_COMM_WORLD, &numprocs);
  MPI_Comm_rank(MPI_COMM_WORLD, &myid);
  MPI_Get_processor_name(processor_name, &proc_namelen);

  /*============================*/
  /* Parse command line:        */
  /*============================*/

  int iret = parse_args(argc, argv);

  cout.setf(ios::left);
  cout << setw(40) << setfill('-') << '-' << endl;
  for (int i=0; i<numprocs; i++) {
    if (myid == i) {
      cout << "Process #" << i << endl;
      print_parm(cout, " ");
      cout << setw(40) << setfill('-') << '-' << endl;
    }
  }

  return 0;
}

//----------------------------------------------------------------------
//----------------------------------------------------------------------
//----------------------------------------------------------------------

void usage(char *);
void set_parm(char *, char *);
void write_parm(void);
void print_default(void);

#define WBUFSIZE 80
#define VBUFSIZE 80

int parse_args(int argc, char **argv)
{
  char *prog=argv[0];
  int c, iparmf=0,iret,i;
  string file;
  char wbuf[WBUFSIZE];
  char vbuf[WBUFSIZE];

  ParamParse prse("=");

  MPI_Comm_rank(MPI_COMM_WORLD, &myid);

  if (myid==0) {

    while (1) {

      c = getopt(argc, argv, "f:dh");
      if (c==-1) break;

      switch (c) {
      case 'f':
	iparmf=1;
	file = optarg;
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
      iret = prse.parse_file(file);
    else {
      iret = prse.parse_argv(argc, &argv[optind]);
      argc -= iret;
    }

    if (argc != 0)
      usage(prog);
  }
  
  
  MPI_Bcast(&iret, 1, MPI_INT, 0, MPI_COMM_WORLD);

				// Set parameters
  spair ret;
  for (i=0; i<iret; i++) {
    if (myid==0) {
      prse.get_next(ret);
      strcpy(wbuf, ret.first.c_str());
      strcpy(vbuf, ret.second.c_str());
    }

				/* Send values to all processes */
    MPI_Bcast(wbuf, WBUFSIZE, MPI_CHAR, 0, MPI_COMM_WORLD);
    MPI_Bcast(vbuf, VBUFSIZE, MPI_CHAR, 0, MPI_COMM_WORLD);
    
    set_parm(wbuf, vbuf);
  }

}

void set_parm(char *word, char *valu)
{
  string wrd(word);

  if (wrd == "NICE")		NICE = atoi(valu);

  else if (wrd == "DENS")	DENS = atoi(valu) ? true : false;

  else if (wrd == "RMAX")	RMAX = atof(valu);

  else if (wrd == "PARMFILE")	PARMFILE = string(valu);

  else {
      cerr << "No such parameter: " << wrd << endl;
      exit(-1);
  }
}

void write_parm()
{
  ofstream fout(PARMFILE.c_str());			
  if ( !fout ) {
    cerr << "Couldn't open parameter file:" << PARMFILE << ". . . quitting\n" << endl;
    exit(-1);
  }

  print_parm(fout,"\0");
}

void print_default()
{
  cerr << "\nDefaults:\n";
  cerr << "----------------------------\n";
  print_parm(cerr,"\0");
  exit(0);
}


void print_parm(ostream& stream, const char *comment)
{
  stream.setf(ios::left);

  stream << comment << setw(20) <<  "NICE" << " = " << NICE << endl;
  stream << comment << setw(20) <<  "RMAX" << " = " << RMAX << endl;
  stream << comment << setw(20) <<  "DENS" << " = " << DENS << endl;
  stream << comment << setw(20) <<  "PARMFILE" << " = " << PARMFILE << endl;

  stream.unsetf(ios::left);
}


void usage(char *prog)
{
  char usage_head[] = 
    "[-f file -d] [keyword=value [keyword=value] .. ]";

  char usage_data[][80] = {
    "     -f file",      "keyword/value parameter file",
    "     -d",           "print default parameters",
    "\nKeywords:",       " DFFILE,OUTFILE",
    "\0"                 };


  cerr.setf(ios::left);

  cerr << "Usage: " << prog << " " << usage_head << endl << endl;

  int j = 0;
  while (*(usage_data[j]) != '\0') {
    cerr << setw(25) << usage_data[j] << usage_data[j+1] << endl;
    j+=2;
  }
  cerr << endl;
}
