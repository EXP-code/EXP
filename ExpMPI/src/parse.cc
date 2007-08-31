/*
  Parse command line
*/

#include "expand.h"

#ifdef RCSID
static char rcsid[] = "$Id$";
#endif
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
  if (parse->find_item("nthrds", val))		nthrds = atoi(val.c_str());
  if (parse->find_item("nbalance", val))	nbalance = atoi(val.c_str());
  if (parse->find_item("dbthresh", val))	dbthresh = atof(val.c_str());

  if (parse->find_item("time", val))		tnow = atof(val.c_str());
  if (parse->find_item("dtime", val))		dtime = atof(val.c_str());
  if (parse->find_item("NICE", val))		NICE = atoi(val.c_str());
  if (parse->find_item("VERBOSE", val))		VERBOSE = atoi(val.c_str());

  if (parse->find_item("multistep", val))	multistep = atoi(val.c_str());

  if (parse->find_item("dynfracV", val))	dynfracV = atof(val.c_str());

  if (parse->find_item("dynfracA", val))	dynfracA = atof(val.c_str());

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

  if (parse->find_item("restart", val)) {
    if (atoi(val.c_str())) restart = true;
    else restart = false;
  }

  if (parse->find_item("homedir", val))		homedir = val;
  if (parse->find_item("ldlibdir", val))	ldlibdir = val;
  if (parse->find_item("infile", val))		infile = val;
  if (parse->find_item("parmfile", val))	parmfile = val;
  if (parse->find_item("ratefile", val))	ratefile = val;
  if (parse->find_item("outdir", val))		outdir = val;
  if (parse->find_item("runtag", val))		runtag = val;

}

void print_parm(ostream& out, char *comment)
{
  out << comment << "[global]" << endl;

  out << comment << " " << "nbodmax"    << " = " << nbodmax     << endl;
  out << comment << " " << "nsteps"     << " = " << nsteps      << endl;
  out << comment << " " << "nthrds"     << " = " << nthrds      << endl;
  out << comment << " " << "nbalance"   << " = " << nbalance    << endl;
  out << comment << " " << "dbthresh"   << " = " << dbthresh    << endl;

  out << comment << " " << "time"       << " = " << tnow        << endl;
  out << comment << " " << "dtime"      << " = " << dtime       << endl;
  out << comment << " " << "NICE"       << " = " << NICE        << endl;
  out << comment << " " << "VERBOSE"    << " = " << VERBOSE     << endl;

  out << comment << " " << "multistep"  << " = " << multistep   << endl;
  out << comment << " " << "posnsync"   << " = " << posnsync    << endl;
  out << comment << " " << "dynfracV"   << " = " << dynfracV    << endl;
  out << comment << " " << "dynfracA"   << " = " << dynfracA    << endl;

  out << comment << " " << "use_cwd"    << " = " << use_cwd     << endl;
  out << comment << " " << "eqmotion"   << " = " << eqmotion    << endl;
  out << comment << " " << "global_cov" << " = " << global_cov  << endl;
  out << comment << " " << "restart"    << " = " << restart     << endl;

  out << comment << " " << "homedir"    << " = " << homedir     << endl;
  out << comment << " " << "ldlibdir"   << " = " << ldlibdir    << endl;
  out << comment << " " << "infile"     << " = " << infile      << endl;
  out << comment << " " << "parmfile"   << " = " << parmfile    << endl;
  out << comment << " " << "ratefile"   << " = " << ratefile    << endl;
  out << comment << " " << "outdir"     << " = " << outdir      << endl;
  out << comment << " " << "runtag"     << " = " << runtag      << endl;
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
