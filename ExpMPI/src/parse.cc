/*
  Parse command line
*/

#include "expand.h"

#ifdef RCSID
static char rcsid[] = "$Id$";
#endif

#include <stdlib.h>
#include <errno.h>
#include <sys/types.h>
#include <sys/stat.h>
#include <unistd.h>
#include <getopt.h>

#include <iostream>
#include <string>

void exp_usage(char *prog)
{
  cerr << prog << " [-f file -d]\n";

  MPI_Abort(MPI_COMM_WORLD, 8);

  exit(-1);
}


void initialize(void)
{
  if (myid==0) {
    cout << "Parameter database:\n"
	 << "-------------------\n\n";
    parse->print_database(cout);
  }

  parse->find_list("global");
  string val;

  if (parse->find_item("nbodmax", val))		nbodmax = atoi(val.c_str());
  if (parse->find_item("nsteps", val))		nsteps = atoi(val.c_str());
  if (parse->find_item("nthrds", val))		nthrds = max<int>(1, atoi(val.c_str()));
  if (parse->find_item("nreport", val))		nreport = atoi(val.c_str());
  if (parse->find_item("nbalance", val))	nbalance = atoi(val.c_str());
  if (parse->find_item("dbthresh", val))	dbthresh = atof(val.c_str());

  if (parse->find_item("time", val))		tnow = atof(val.c_str());
  if (parse->find_item("dtime", val))		dtime = atof(val.c_str());
  if (parse->find_item("PFbufsz", val))         PFbufsz = atoi(val.c_str());
  if (parse->find_item("NICE", val))		NICE = atoi(val.c_str());
  if (parse->find_item("VERBOSE", val))		VERBOSE = atoi(val.c_str());
  if (parse->find_item("runtime", val))		runtime = atof(val.c_str());

  if (parse->find_item("multistep", val))	multistep = atoi(val.c_str());
  if (parse->find_item("centerlevl", val))	centerlevl = atoi(val.c_str());

  if (parse->find_item("dynfracV", val))	dynfracV = atof(val.c_str());

  if (parse->find_item("dynfracA", val))	dynfracA = atof(val.c_str());

  if (parse->find_item("DTold", val)) {
    if (atoi(val.c_str())) {
      DTold = true;
      if (myid==0)
	cout << "parse: using original (old) time-step algorithm" << endl;
    }
    else DTold = false;
  }

  if (parse->find_item("use_cwd", val)) {
    if (atoi(val.c_str())) use_cwd = true;
    else use_cwd = false;
  }

  if (parse->find_item("posnsync", val)) {
    if (atoi(val.c_str())) posnsync = true;
    else posnsync = false;
  }

  if (parse->find_item("eqmotion", val)) {
    if (atoi(val.c_str())) eqmotion = true;
    else eqmotion = false;
  }

  if (parse->find_item("global_cov", val)) {
    if (atoi(val.c_str())) global_cov = true;
    else global_cov = false;
  }

  if (parse->find_item("homedir", val))		homedir = val;
  if (parse->find_item("ldlibdir", val))	ldlibdir = val;
  if (parse->find_item("infile", val))		infile = val;
  if (parse->find_item("parmfile", val))	parmfile = val;
  if (parse->find_item("ratefile", val))	ratefile = val;
  if (parse->find_item("runtag", val))		runtag = val;
  if (parse->find_item("restart_cmd", val))     restart_cmd = val;

  if (parse->find_item("outdir", val)) {
    bool ok = true;
				// Root node with check existence and try
				// to create directory if need be . . .
    if (myid == 0) {
      struct stat sb;		// Stat buffer structure
      if (stat(val.c_str(), &sb) == -1) {
				// Error in opening the candidate directory
	cout << "parse: I can't open directory <" << val << ">" << endl;
				// Maybe we need to create it?
	cout << "parse: I will attempt to create it" << endl;
	if (mkdir(val.c_str(), 0755) == -1) {
	  cout << "parse: error creating directory <" << val 
	       << ">, aborting" << endl;
	  perror("mkdir");
	  MPI_Abort(MPI_COMM_WORLD, 10);
	  exit(EXIT_SUCCESS);
	}
      } else {
	cout << "parse: output path <" << val << "> is "; // 
	switch (sb.st_mode & S_IFMT) {
	case S_IFBLK:  cout << "a block device";     ok=false;  break;
	case S_IFCHR:  cout << "a character device"; ok=false;  break;
	case S_IFDIR:  cout << "a directory";        ok=true;   break;
	case S_IFIFO:  cout << "a FIFO/pipe";        ok=false;  break;
	case S_IFLNK:  cout << "a symlink";          ok=false;  break;
	case S_IFREG:  cout << "a regular file";     ok=false;  break;
	case S_IFSOCK: cout << "a socket";           ok=false;  break;
	default:       cout << "an unknown type";    ok=false;  break;
	}
	if (ok) cout << " . . . good"     << endl;
	else    cout << " . . . very bad" << endl;
      }
    }

    MPI_Barrier(MPI_COMM_WORLD);
				// Append "/" if necessary
    if (val[val.size()-1] != '/') val += "/";
    
    ostringstream tfile;
    if (myid==0 && ok) {	// Try to open a file in this path
      tfile << val << ".test.file." << rand();

      if (mknod(tfile.str().c_str(),  S_IFREG | 0666, 0)==-1) {
				// If the file exists, ok, otherwise NOT ok
	if (errno != EEXIST) ok = false;
      }
				// Clean up
      unlink(tfile.str().c_str()); 

				// Print results to output log
      cout << "parse: ";
      if (ok) 
	cout << "test file opened and deleted in output directory, good!";
      else    
	cout << "we can't open files in this path, attempted to create: " 
	     << tfile.str();
      cout << endl;
    }

				// Send directory status to all processes
    int iok = static_cast<int>(ok);
    MPI_Bcast(&iok, 1, MPI_INT, 0, MPI_COMM_WORLD);

    if (!ok) {
      MPI_Abort(MPI_COMM_WORLD, 11);
      exit(255);
    }

				// Ok, we have an output directory!
    outdir = val;
  }

}

void print_parm(ostream& out, const char *comment)
{
  out << comment << "[global]" << endl;

  out << comment << " " << "nbodmax"    << " = " << nbodmax     << endl;
  out << comment << " " << "nsteps"     << " = " << nsteps      << endl;
  out << comment << " " << "nthrds"     << " = " << nthrds      << endl;
  out << comment << " " << "nreport"    << " = " << nreport     << endl;
  out << comment << " " << "nbalance"   << " = " << nbalance    << endl;
  out << comment << " " << "dbthresh"   << " = " << dbthresh    << endl;

  out << comment << " " << "time"       << " = " << tnow        << endl;
  out << comment << " " << "dtime"      << " = " << dtime       << endl;
  out << comment << " " << "PFbufsz"    << " = " << PFbufsz     << endl;
  out << comment << " " << "NICE"       << " = " << NICE        << endl;
  out << comment << " " << "VERBOSE"    << " = " << VERBOSE     << endl;
  out << comment << " " << "runtime"    << " = " << runtime     << endl;

  out << comment << " " << "multistep"  << " = " << multistep   << endl;
  out << comment << " " << "centerlevl" << " = " << centerlevl  << endl;
  out << comment << " " << "posnsync"   << " = " << posnsync    << endl;
  out << comment << " " << "DTold"      << " = " << DTold       << endl;
  out << comment << " " << "dynfracV"   << " = " << dynfracV    << endl;
  out << comment << " " << "dynfracA"   << " = " << dynfracA    << endl;

  out << comment << " " << "use_cwd"    << " = " << use_cwd     << endl;
  out << comment << " " << "eqmotion"   << " = " << eqmotion    << endl;
  out << comment << " " << "global_cov" << " = " << global_cov  << endl;

  out << comment << " " << "homedir"    << " = " << homedir     << endl;
  out << comment << " " << "ldlibdir"   << " = " << ldlibdir    << endl;
  out << comment << " " << "infile"     << " = " << infile      << endl;
  out << comment << " " << "parmfile"   << " = " << parmfile    << endl;
  out << comment << " " << "ratefile"   << " = " << ratefile    << endl;
  out << comment << " " << "outdir"     << " = " << outdir      << endl;
  out << comment << " " << "runtag"     << " = " << runtag      << endl;
  out << comment << " " << "command"    << " = " << restart_cmd << endl;
}


void write_parm(void)
{
  if (myid!=0) return;
  string curparm(outdir + parmfile + "." + runtag);
  ofstream out(curparm.c_str());
  if (!out) {
    cerr << "write_parm: could not open <" << parmfile << ">\n";
    MPI_Abort(MPI_COMM_WORLD, 102);
    exit(0);
  }
  
  print_parm(out, "");

  out << endl
      << "--------------------" << endl
      << " Parameter database " << endl
      << "--------------------" << endl
      << endl;

  parse->print_database(out);
}


void print_default(void)
{
  print_parm(cout, "# ");
}


void MPL_parse_args(int argc, char** argv)
{
  extern char *optarg;
  char *prog=argv[0];
  int myid;
  char file[128];
  int c;
  string curparm(outdir + parmfile + "." + runtag);

  MPI_Comm_rank(MPI_COMM_WORLD, &myid);
  
  if (myid==0) {

    while ((c = getopt(argc, argv, "f:dvh")) != -1)
      switch (c) {
      case 'f': 
	curparm.erase(curparm.begin(), curparm.end());
	curparm = optarg;
	break;
      case 'd':
	print_default();
	exit(0);
	break;
      case 'v':
	exit(0);
	break;
      case '?':
      case 'h':
	exp_usage(prog);
	break;
      }
    
    ifstream in(curparm.c_str());
    if (!in) {
      char mbuf[512];
      cerr << "MAIN: can not open parameter file <" << parmfile << ">\n";
      cerr << "MAIN: pwd is <" << getcwd(mbuf, 512) << ">\n";
      MPI_Abort(MPI_COMM_WORLD, 9);
    }

    parse = new ParamParseMPI(&in, ":");

  } else {
    
    parse = new ParamParseMPI(NULL, ":");

  }

  initialize();

}
