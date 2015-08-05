/*
  Parse command line
*/

#include "expand.h"

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
  cerr << "Usage" << endl
       << "-----" << endl << endl
       << "[mpirun -np N ...] " << prog << " [-f file -d -h]" << endl << endl
       << "where" << endl << endl
       << "  -f file   specifies the input parameter file" << endl
       << "  -d        displays a default parameter file"  << endl
       << "  -h        shows this help"                    << endl  << endl
       << "See EXP/doc/html/index.html for extensive documentation" << endl
       << endl;
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
  if (parse->find_item("nbits", val))           nbits = atoi(val.c_str());
  if (parse->find_item("pkbits", val))          pkbits = atoi(val.c_str());
  if (parse->find_item("PFbufsz", val))         PFbufsz = atoi(val.c_str());
  if (parse->find_item("NICE", val))		NICE = atoi(val.c_str());
  if (parse->find_item("VERBOSE", val))		VERBOSE = atoi(val.c_str());
  if (parse->find_item("rlimit", val))          rlimit_val = atoi(val.c_str());
  if (parse->find_item("runtime", val))		runtime = atof(val.c_str());

  if (parse->find_item("multistep", val))	multistep = atoi(val.c_str());
  if (parse->find_item("centerlevl", val))	centerlevl = atoi(val.c_str());

  if (parse->find_item("dynfracS", val))	dynfracS = atof(val.c_str());
  if (parse->find_item("dynfracD", val))	dynfracD = atof(val.c_str());
  if (parse->find_item("dynfracV", val))	dynfracV = atof(val.c_str());
  if (parse->find_item("dynfracA", val))	dynfracA = atof(val.c_str());
  if (parse->find_item("dynfracP", val))	dynfracP = atof(val.c_str());

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

  if (parse->find_item("eqmotion", val)) {
    if (atoi(val.c_str())) eqmotion = true;
    else eqmotion = false;
  }

  if (parse->find_item("global_cov", val)) {
    if (atoi(val.c_str())) global_cov = true;
    else global_cov = false;
  }

  if (parse->find_item("barrier_check", val)) {
    if (atoi(val.c_str())) barrier_check = true;
    else barrier_check = false;
  }

  if (parse->find_item("barrier_debug", val)) {
    if (atoi(val.c_str())) barrier_debug = true;
    else barrier_debug = false;
  }

  if (parse->find_item("barrier_extra", val)) {
    if (atoi(val.c_str())) barrier_extra = true;
    else barrier_extra = false;
  }

  if (parse->find_item("barrier_label", val)) {
    if (atoi(val.c_str())) barrier_label = true;
    else barrier_label = false;
  }

  if (parse->find_item("barrier_heavy", val)) {
    if (atoi(val.c_str())) barrier_light = false;
    else barrier_light = true;
  }

  if (parse->find_item("barrier_quiet", val)) {
    if (atoi(val.c_str())) barrier_quiet = true;
    else barrier_quiet = false;
  }

  if (parse->find_item("barrier_verbose", val)) {
    if (atoi(val.c_str())) barrier_quiet = false;
    else barrier_quiet = true;
  }

  if (parse->find_item("homedir", val)) {
    // Check for and add trailing slash
    if (*val.rbegin() != '/') val += '/';
    homedir = val;
  }

  if (parse->find_item("ldlibdir", val))	ldlibdir = val;
  if (parse->find_item("infile", val))		infile = val;
  if (parse->find_item("parmfile", val))	parmfile = val;
  if (parse->find_item("ratefile", val))	ratefile = val;
  if (parse->find_item("runtag", val))		runtag = val;
  if (parse->find_item("restart_cmd", val))     restart_cmd = val;

  if (parse->find_item("outdir", val)) {
    bool ok = true;
				// Check for and add trailing slash
    if (*val.rbegin() != '/') val += '/';

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

  out << endl;

  out << comment << setw(20) << "nbodmax"    << " : " << nbodmax     << endl;
  out << comment << setw(20) << "nsteps"     << " : " << nsteps      << endl;
  out << comment << setw(20) << "nthrds"     << " : " << nthrds      << endl;
  out << comment << setw(20) << "nreport"    << " : " << nreport     << endl;
  out << comment << setw(20) << "nbalance"   << " : " << nbalance    << endl;
  out << comment << setw(20) << "dbthresh"   << " : " << dbthresh    << endl;

  out << comment << setw(20) << "time"       << " : " << tnow        << endl;
  out << comment << setw(20) << "dtime"      << " : " << dtime       << endl;
  out << comment << setw(20) << "nbits"      << " : " << nbits       << endl;
  out << comment << setw(20) << "pkbits"     << " : " << pkbits      << endl;
  out << comment << setw(20) << "PFbufsz"    << " : " << PFbufsz     << endl;
  out << comment << setw(20) << "NICE"       << " : " << NICE        << endl;
  out << comment << setw(20) << "VERBOSE"    << " : " << VERBOSE     << endl;
  out << comment << setw(20) << "rlimit"     << " : " << rlimit_val  << endl;
  out << comment << setw(20) << "runtime"    << " : " << runtime     << endl;

  out << comment << setw(20) << "multistep"  << " : " << multistep   << endl;
  out << comment << setw(20) << "centerlevl" << " : " << centerlevl  << endl;
  out << comment << setw(20) << "DTold"      << " : " << DTold       << endl;
  out << comment << setw(20) << "dynfracS"   << " : " << dynfracS    << endl;
  out << comment << setw(20) << "dynfracV"   << " : " << dynfracV    << endl;
  out << comment << setw(20) << "dynfracA"   << " : " << dynfracA    << endl;
  out << comment << setw(20) << "dynfracP"   << " : " << dynfracP    << endl;

  out << comment << setw(20) << "use_cwd"    << " : " << use_cwd     << endl;
  out << comment << setw(20) << "eqmotion"   << " : " << eqmotion    << endl;
  out << comment << setw(20) << "global_cov" << " : " << global_cov  << endl;

  out << comment << setw(20) << "homedir"    << " : " << homedir     << endl;
  out << comment << setw(20) << "ldlibdir"   << " : " << ldlibdir    << endl;
  out << comment << setw(20) << "infile"     << " : " << infile      << endl;
  out << comment << setw(20) << "parmfile"   << " : " << parmfile    << endl;
  out << comment << setw(20) << "ratefile"   << " : " << ratefile    << endl;
  out << comment << setw(20) << "outdir"     << " : " << outdir      << endl;
  out << comment << setw(20) << "runtag"     << " : " << runtag      << endl;
  out << comment << setw(20) << "command"    << " : " << restart_cmd << endl;

  out << endl;

  out << comment << "[components]" << endl;

  out << endl;

  out << comment << "[output]" << endl;

  out << endl;

  out << comment << setw(20) << "outlog"    << " : " << "nint=1"  << endl;
  out << comment << setw(20) << "outpsn"    << " : " << "nint=10" << endl;

  out << endl;

  out << comment << "[external]" << endl;

  out << endl;
}


void write_parm(void)
{
  if (myid!=0) return;
  string curparm(outdir + parmfile + "." + runtag);
  ofstream out(curparm.c_str());
  if (!out) {
    cerr << "write_parm: could not open <" << parmfile << ">\n";
    MPI_Abort(MPI_COMM_WORLD, 102);
    exit(EXIT_SUCCESS);
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
  
  int done = 0;

  if (myid==0) {

    while ((c = getopt(argc, argv, "f:dvh")) != -1)
      switch (c) {
      case 'f': 
	curparm.clear();
	curparm = optarg;
	break;
      case 'd':
	print_default();
	done = 1;
	break;
      case '?':
      case 'h':
	exp_usage(prog);
	done = 1;
      }
    
    MPI_Bcast(&done, 1, MPI_INT, 0, MPI_COMM_WORLD);

    if (done) {
      MPI_Finalize();
      exit(EXIT_SUCCESS);
    }

    ifstream in(curparm.c_str());
    if (!in) {
      char mbuf[512];
      cerr << "MAIN: can not open parameter file <" << parmfile << ">\n";
      cerr << "MAIN: pwd is <" << getcwd(mbuf, 512) << ">\n";
      done = 1;
    }

    MPI_Bcast(&done, 1, MPI_INT, 0, MPI_COMM_WORLD);

    if (done) {
      MPI_Finalize();
      exit(EXIT_SUCCESS);
    }

    parse = new ParamParseMPI(&in, ":");

  } else {
    
    MPI_Bcast(&done, 1, MPI_INT, 0, MPI_COMM_WORLD);

    if (done) {
      MPI_Finalize();
      exit(EXIT_SUCCESS);
    }

    MPI_Bcast(&done, 1, MPI_INT, 0, MPI_COMM_WORLD);

    if (done) {
      MPI_Finalize();
      exit(EXIT_SUCCESS);
    }

    parse = new ParamParseMPI(NULL, ":");

  }

  parse->parse_argv(argc-optind, &argv[optind]);

  initialize();

}
