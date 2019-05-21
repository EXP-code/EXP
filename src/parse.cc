//
//  Parse command line
//

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
    cout << std::string(72, '=')  << std::endl
	 << "Parameter database:" << std::endl
	 << "-------------------" << std::endl << std::endl
	 << parse << std::endl
	 << std::string(72, '=')  << std::endl << std::endl;
  }

  YAML::Node _G;
  try {
    _G = parse["Global"];
  }
  catch (YAML::Exception & error) {
    if (myid==0) std::cout << "Error parsing 'global' stanza: "
			   << error.what() << std::endl
			   << std::string(60, '-') << std::endl
			   << "Config node"        << std::endl
			   << std::string(60, '-') << std::endl
			   << parse                << std::endl
			   << std::string(60, '-') << std::endl;
    MPI_Finalize();
    exit(-1);
  }

  if (_G) {
    
    if (_G["nsteps"])	     nsteps     = _G["nsteps"].as<int>();
    if (_G["nthrds"])	     nthrds     = std::max<int>(1, _G["nthrds"].as<int>());
    if (_G["nreport"])	     nreport    = _G["nreport"].as<int>();
    if (_G["nbalance"])      nbalance   = _G["nbalance"].as<int>();
    if (_G["dbthresh"])      dbthresh   = _G["dbthresh"].as<int>();
    
    if (_G["time"])          tnow       = _G["time"].as<double>();
    if (_G["dtime"])         dtime      = _G["dtime"].as<double>();
    if (_G["nbits"])         nbits      = _G["nbits"].as<int>();
    if (_G["pkbits"])        pkbits     = _G["pkbits"].as<int>();
    if (_G["PFbufsz"])       PFbufsz    = _G["PFbusz"].as<int>();
    if (_G["NICE"])          NICE       = _G["NICE"].as<int>();
    if (_G["VERBOSE"])       VERBOSE    = _G["VERBOSE"].as<int>();
    if (_G["rlimit"])        rlimit_val = _G["rlimit"].as<int>();
    if (_G["runtime"])       runtime    = _G["runtime"].as<double>();
    
    if (_G["multistep"])     multistep  = _G["multistep"].as<int>();
    if (_G["centerlevl"])    centerlevl = _G["centerlevl"].as<int>();

    if (_G["dynfracS"])	     dynfracS   = _G["dynfracS"].as<double>();
    if (_G["dynfracD"])	     dynfracD   = _G["dynfracD"].as<double>();
    if (_G["dynfracV"])	     dynfracV   = _G["dynfracV"].as<double>();
    if (_G["dynfracA"])	     dynfracA   = _G["dynfracA"].as<double>();
    if (_G["dynfracP"])	     dynfracP   = _G["dynfracP"].as<double>();

    if (_G["cuStreams"])     cuStreams  = _G["cuStreams"].as<int>();

    if (_G["DTold"]) {
      DTold = _G["DTold"].as<bool>();
      if (DTold and myid==0)
	cout << "parse: using original (old) time-step algorithm" << endl;
    }
    
    if (_G["use_cwd"])       use_cwd    = _G["use_cwd"].as<bool>();
    if (_G["eqmotion"])      eqmotion   = _G["eqmotion"].as<bool>();
    if (_G["global_cov"])    global_cov = _G["global_cov"].as<bool>();
    if (_G["cuda_prof"])     cuda_prof  = _G["cuda_prof"].as<bool>();
    if (_G["use_cuda"])      use_cuda   = _G["use_cuda"].as<bool>();
    if (_G["barrier_check"]) barrier_check = _G["barrier_check"].as<bool>();
    if (_G["barrier_debug"]) barrier_debug = _G["barrier_debug"].as<bool>();
    if (_G["barrier_extra"]) barrier_extra = _G["barrier_extra"].as<bool>();
    if (_G["barrier_label"]) barrier_extra = _G["barrier_label"].as<bool>();
    if (_G["barrier_heavy"]) barrier_extra = _G["barrier_heavy"].as<bool>();
    if (_G["barrier_quiet"]) barrier_quiet = _G["barrier_quiet"].as<bool>();
    if (_G["barrier_verbose"]) barrier_quiet = not _G["barrier_quiet"].as<bool>();

    if (_G["main_wait"])     main_wait  = _G["main_wait"].as<bool>();
    if (_G["debug_wait"]) {
      debug_wait = _G["debug_wait"].as<bool>();
      std::cout << "Found <debug_wait=" << std::boolalpha
		<< debug_wait << ">" << std::endl;
      if (debug_wait) {
	std::cout << "----" << std::endl;
	if (main_wait)
	  std::cout << "----  Main process will wait in a loop until gdb is attached and the loop is freed" << std::endl;
	else
	  std::cout << "---- All processes will wait in a loop until gdb is attached and the loop is freed" << std::endl;
	std::cout << "---- by setting 'debug_wait = false'" << std::endl
		  << "----" << std::endl;
      }
    }
    
    if (_G["mpi_wait"]) {
      if (_G["mpi_wait"].as<bool>()) {
	mpi_wait  = true;
      } else mpi_wait = false;
      if (myid==0) {
	std::cout << "Found <mpi_wait=" << std::boolalpha
		  << mpi_wait << ">" << std::endl;
	if (mpi_wait)
	  std::cout << "----" << std::endl
		    << "---- When the MPI error handler is called, process will spin, waiting for a gdb connection." << std::endl
		    << "---- Messages describing the affected node and pid will be written to the" << std::endl
		    << "---- standard output." << std::endl
		    << "----" << std::endl;
      }
    }

    if (_G["fpe_trap"]) {
      if (_G["fpe_trap"].as<bool>()) {
	fpe_trap  = true;
	fpe_trace = false;
	fpe_wait  = false;
      } else fpe_trap = false;
      if (myid==0) {
	std::cout << "Found <fpe_trap=" << std::boolalpha
		  << fpe_trap << ">" << std::endl;
	if (fpe_trap) std::cout << "----" << std::endl
				<< "---- Set a breakpoint in fpetrap.h:21 to catch and debug FP errors" << std::endl
				<< "----" << std::endl;
      }
    }

    if (_G["fpe_trace"]) {
      if (_G["fpe_trace"].as<bool>()) {
	fpe_trap  = false;
	fpe_trace = true;
	fpe_wait  = false;
      } else fpe_trace = false;
      if (myid==0) {
	std::cout << "Found <fpe_trace=" << std::boolalpha
		  << fpe_trace << ">" << std::endl;
	if (fpe_trace) std::cout << "----" << std::endl
				 << "---- Print a backtrace to stderr on detecting an FP error" << std::endl
				 << "----" << std::endl;
      }
    }

    if (_G["fpe_wait"]) {
      if (_G["fpe_wait"].as<bool>()) {
	fpe_trap  = false;
	fpe_trace = false;
	fpe_wait  = true;
      } else fpe_wait = false;
      if (myid==0) {
	std::cout << "Found <fpe_wait=" << std::boolalpha
		  << fpe_wait << ">" << std::endl;
	if (fpe_wait)
	  std::cout << "----" << std::endl
		    << "---- When an FPE is signalled, process will spin, waiting for a gdb connection." << std::endl
		    << "---- Messages describing the affected node and pid will be written to the" << std::endl
		    << "---- standard output." << std::endl
		    << "----" << std::endl;
      }
    }

    if (_G["homedir"]) {
      homedir = _G["homedir"].as<std::string>();
      // Check for and add trailing slash
      if (*homedir.rbegin() != '/') homedir += '/';
    }

    if (_G["ldlibdir"])	        ldlibdir = _G["ldlibdir"].as<std::string>();
    if (_G["infile"])		infile   = _G["infile"].as<std::string>();
    if (_G["parmfile"])	        parmfile = _G["parmfile"].as<std::string>();
    if (_G["ratefile"])	        ratefile = _G["ratefile"].as<std::string>();
    if (_G["runtag"])		runtag   = _G["runtag"].as<std::string>();
    if (_G["restart_cmd"])      restart_cmd = _G["restart_cmd"].as<std::string>();
    if (_G["ignore_info"])      ignore_info = _G["ignore_info"].as<bool>();
    
    bool ok = true;

    if (_G["outdir"]) {

      outdir = _G["outdir"].as<std::string>();

				// Check for and add trailing slash
      if (*outdir.rbegin() != '/') outdir += '/';

				// Root node with check existence and try
				// to create directory if need be . . .
      int nOK = 0;

      if (myid == 0) {
	struct stat sb;		// Stat buffer structure
	if (stat(outdir.c_str(), &sb) == -1) {
	  // Error in opening the candidate directory
	  cout << "parse: I can't open directory <" << outdir << ">" << endl;
	  // Maybe we need to create it?
	  cout << "parse: I will attempt to create it" << endl;
	  if (mkdir(outdir.c_str(), 0755) == -1) {
	    cout << "parse: error creating directory <" << outdir
		 << ">, aborting" << endl;
	    perror("mkdir");
	    nOK = 1;
	  }
	} else {
	  cout << "parse: output path <" << outdir << "> is "; // 
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

      MPI_Bcast(&nOK, 1, MPI_INT, 0, MPI_COMM_WORLD);
      if (nOK) {
	MPI_Finalize();
	exit(10);
      }

    }

    MPI_Barrier(MPI_COMM_WORLD);
				// Append "/" if necessary
    if (outdir[outdir.size()-1] != '/') outdir += "/";
    
    ostringstream tfile;
    if (myid==0 && ok) {	// Try to open a file in this path
      tfile << outdir << ".test.file." << rand();

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

    if (!iok) {
      MPI_Finalize();
      exit(11);
    }

  }

}

void update_parm()
{
  try {
    YAML::Node conf = parse["Global"];

    if (not conf["nsteps"])     conf["nsteps"]      = nsteps;
    if (not conf["nthrds"])     conf["nthrds"]      = nthrds;
    if (not conf["ngpus"])      conf["ngpus"]       = ngpus;
    if (not conf["nreport"])    conf["nreport"]     = nreport;
    if (not conf["nbalance"])   conf["nbalance"]    = nbalance;
    if (not conf["dbthresh"])   conf["dbthresh"]    = dbthresh;
    
    if (not conf["time"])       conf["time"]        = tnow;
    if (not conf["dtime"])      conf["dtime"]       = dtime;
    if (not conf["nbits"])      conf["nbits"]       = nbits;
    if (not conf["pkbits"])     conf["pkbits"]      = pkbits;
    if (not conf["PFbufsz"])    conf["PFbufsz"]     = PFbufsz;
    if (not conf["NICE"])       conf["NICE"]        = NICE;
    if (not conf["VERBOSE"])    conf["VERBOSE"]     = VERBOSE;
    if (not conf["rlimit"])     conf["rlimit"]      = rlimit_val;
    if (not conf["runtime"])    conf["runtime"]     = runtime;
    
    if (not conf["multistep"])  conf["multistep"]   = multistep;
    if (not conf["centerlevl"]) conf["centerlevl"]  = centerlevl;
    if (not conf["DTold"])      conf["DTold"]       = DTold;
    if (not conf["dynfracS"])   conf["dynfracS"]    = dynfracS;
    if (not conf["dynfracV"])   conf["dynfracV"]    = dynfracV;
    if (not conf["dynfracA"])   conf["dynfracA"]    = dynfracA;
    if (not conf["dynfracP"])   conf["dynfracP"]    = dynfracP;
    
    if (not conf["use_cwd"])    conf["use_cwd"]     = use_cwd;
    if (not conf["eqmotion"])   conf["eqmotion"]    = eqmotion;
    if (not conf["global_cov"]) conf["global_cov"]  = global_cov;

    if (not conf["homedir"])    conf["homedir"]     = homedir;
    if (not conf["ldlibdir"])   conf["ldlibdir"]    = ldlibdir;
    if (not conf["infile"])     conf["infile"]      = infile;
    if (not conf["parmfile"])   conf["parmfile"]    = parmfile;
    if (not conf["ratefile"])   conf["ratefile"]    = ratefile;
    if (not conf["outdir"])     conf["outdir"]      = outdir;
    if (not conf["runtag"])     conf["runtag"]      = runtag;
    if (not conf["command"])    conf["command"]    = restart_cmd;
    
    parse["Global"] = conf;
  }
  catch (YAML::Exception & error) {
    if (myid==0) std::cout << "Error updating parameters in parse.update_parm: "
			   << error.what() << std::endl;
    MPI_Finalize();
    exit(-1);
  }
}


void write_parm(void)
{
  int nOK = 0;

  if (myid==0) {
    string curparm(outdir + parmfile + "." + runtag + ".yml");
    ofstream out(curparm.c_str());
    if (!out) {
      std::cerr << "write_parm: could not open <" << parmfile << ">\n";
      nOK = 1;
    }
  }
  
  MPI_Bcast(&nOK, 1, MPI_INT, 0, MPI_COMM_WORLD);
  if (nOK) {
    MPI_Finalize();
    exit(102);
  }

  if (myid!=0) return;

  update_parm();

  out << std::endl
      << "---"                       << std::endl
      << "#------------------------" << std::endl
      << "# Parameter database     " << std::endl
      << "# EXP [" << VERSION << "]" << std::endl
      << "#------------------------" << std::endl
      << "#"                         << std::endl;

  out << parse << std::endl
      << "..." << std::endl;
}


void print_default(void)
{
  update_parm();
  
  std::cout << "# EXP [" << VERSION << "]" << std::endl
	    << parse << std::endl;
}


void YAML_parse_args(int argc, char** argv)
{
  extern char *optarg;
  char *prog = argv[0];
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

    try {
      parse = YAML::Load(in);
    }
    catch (YAML::Exception & error) {
      if (myid==0) std::cout << "Error parsing config file: "
			     << error.what() << std::endl;
      MPI_Finalize();
      exit(-1);
    }

    std::ostringstream serial;
    serial << parse << std::endl;

    int line_size = serial.str().size();
    MPI_Bcast(&line_size, 1, MPI_INT, 0, MPI_COMM_WORLD);
    MPI_Bcast(const_cast<char*>(serial.str().data()), line_size, MPI_CHAR, 0, MPI_COMM_WORLD);

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

    int line_size;
    MPI_Bcast(&line_size, 1, MPI_INT, 0, MPI_COMM_WORLD);

    std::string config;
    config.resize(line_size);

    MPI_Bcast(const_cast<char*>(config.data()), line_size, MPI_CHAR, 0, MPI_COMM_WORLD);
    
    std::istringstream sin(config);

    parse = YAML::Load(sin);

  }

  initialize();
}
