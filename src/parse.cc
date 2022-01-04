//
//  Parse command line
//

#include "expand.H"
#include <cxxopts.H>

#include <cstdlib>
#include <errno.h>
#include <sys/types.h>
#include <sys/stat.h>
#include <unistd.h>

#include <iostream>
#include <string>

void exp_version()
{
  
  const int W = 80;		// Full linewidth
  std::ostringstream sout;	// Get a std::string from the string
				// literal
  sout << "%%%%% This is " << PACKAGE_STRING << " ";
				// Print the info block
  std::cout << std::endl
	    << std::setw(W) << std::setfill('%') << '%' << std::endl
	    << std::left << setw(W) << sout.str() << std::endl
	    << std::setw(W) << std::setfill('%') << '%' << std::endl
	    << std::setfill(' ')
	    << std::setw(20) << "%%%%% Repository URL" << " | "
	    << std::setw(W-24) << PACKAGE_URL << '%' << std::endl
	    << std::setw(20) << "%%%%% Current branch" << " | "
	    << std::setw(W-24) << GIT_BRANCH << '%' << std::endl
	    << std::setw(20) << "%%%%% Current commit" << " | "
	    << std::setw(W-24) << GIT_COMMIT << '%' << std::endl
	    << std::setw(20) << "%%%%% Compile time"   << " | "
	    << std::setw(W-24) << COMPILE_TIME << '%' << std::endl
	    << std::setfill('%')
	    << std::setw(W) << '%' << std::setfill(' ') << std::endl
	    << std::endl;
}

void initialize(void)
{
  //  +---Do not dump yaml config to standard output unless VERRBOSE >= 7
  //  |
  //  V
  if (VERBOSE>6 and myid==0) {
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
    if (_G["ngpus"])	     ngpus      = _G["ngpus"].as<int>();
    if (_G["nreport"])	     nreport    = _G["nreport"].as<int>();
    if (_G["nbalance"])      nbalance   = _G["nbalance"].as<int>();
    if (_G["dbthresh"])      dbthresh   = _G["dbthresh"].as<double>();
    
    if (_G["time"])          tnow       = _G["time"].as<double>();
    if (_G["dtime"])         dtime      = _G["dtime"].as<double>();
    if (_G["nbits"])         nbits      = _G["nbits"].as<int>();
    if (_G["pkbits"])        pkbits     = _G["pkbits"].as<int>();
    if (_G["PFbufsz"])       PFbufsz    = _G["PFbufsz"].as<int>();
    if (_G["NICE"])          NICE       = _G["NICE"].as<int>();
    if (_G["VERBOSE"])       VERBOSE    = _G["VERBOSE"].as<int>();
    if (_G["rlimit"])        rlimit_val = _G["rlimit"].as<int>();
    if (_G["runtime"])       runtime    = _G["runtime"].as<double>();
    
    if (_G["multistep"])     multistep  = _G["multistep"].as<int>();
    if (_G["shiftlevl"])     shiftlevl  = _G["shiftlevl"].as<int>();
    if (_G["centerlevl"])    centerlevl = _G["centerlevl"].as<int>();

    if (_G["dynfracS"])	     dynfracS   = _G["dynfracS"].as<double>();
    if (_G["dynfracD"])	     dynfracD   = _G["dynfracD"].as<double>();
    if (_G["dynfracV"])	     dynfracV   = _G["dynfracV"].as<double>();
    if (_G["dynfracA"])	     dynfracA   = _G["dynfracA"].as<double>();
    if (_G["dynfracP"])	     dynfracP   = _G["dynfracP"].as<double>();

    if (_G["DTold"]) {
      DTold = _G["DTold"].as<bool>();
      if (DTold and myid==0)
	cout << "---- Using original (old) time-step algorithm" << endl;
    }
    
    if (_G["random_seed"])     random_seed   = _G["random_seed"].as<unsigned int>();
    if (myid==0) {
      cout << "---- Random seed for Node 0 is: " << random_seed << std::endl;
    }

    // Initialize random number generator
    random_gen.seed(random_seed+static_cast<unsigned int>(myid));

    if (_G["use_cwd"])         use_cwd       = _G["use_cwd"].as<bool>();
    if (_G["eqmotion"])        eqmotion      = _G["eqmotion"].as<bool>();
    if (_G["global_cov"])      global_cov    = _G["global_cov"].as<bool>();
    if (_G["cuda_prof"])       cuda_prof     = _G["cuda_prof"].as<bool>();
    if (_G["cuda"])            use_cuda      = _G["cuda"].as<bool>();
    if (_G["use_cuda"])        use_cuda      = _G["use_cuda"].as<bool>();
#if HAVE_LIBCUDA != 1
    use_cuda = false;
#endif
    if (_G["barrier_check"])   barrier_check = _G["barrier_check"].as<bool>();
    if (_G["barrier_debug"])   barrier_debug = _G["barrier_debug"].as<bool>();
    if (_G["barrier_extra"])   barrier_extra = _G["barrier_extra"].as<bool>();
    if (_G["barrier_label"])   barrier_extra = _G["barrier_label"].as<bool>();
    if (_G["barrier_heavy"])   barrier_extra = _G["barrier_heavy"].as<bool>();
    if (_G["barrier_quiet"])   barrier_quiet = _G["barrier_quiet"].as<bool>();
    if (_G["barrier_verbose"]) barrier_quiet = not _G["barrier_quiet"].as<bool>();

    if (_G["gdb_trace"])       gdb_trace  = _G["gdb_trace"].as<bool>();
    if (_G["main_wait"])       main_wait  = _G["main_wait"].as<bool>();
    if (_G["debug_wait"]) {
      debug_wait = _G["debug_wait"].as<bool>();
      if (myid==0) {
	std::cout << "---- Found <debug_wait=" << std::boolalpha
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
    }
    
    if (_G["mpi_wait"]) {
      if (_G["mpi_wait"].as<bool>()) {
	mpi_wait  = true;
      } else mpi_wait = false;
      if (myid==0) {
	std::cout << "---- Found <mpi_wait=" << std::boolalpha
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
	std::cout << "---- Found <fpe_trap=" << std::boolalpha
		  << fpe_trap << ">" << std::endl;
	if (fpe_trap) std::cout << "----" << std::endl
				<< "---- Set a breakpoint in fpetrap.h:25 to catch and debug FP errors" << std::endl
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
	std::cout << "---- Found <fpe_trace=" << std::boolalpha
		  << fpe_trace << ">" << std::endl
		  << "---- Found <gdb_trace=" << std::boolalpha
		  << gdb_trace << ">" << std::endl;
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
	std::cout << "---- Found <fpe_wait=" << std::boolalpha
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
	  cout << "---- I can't open directory <" << outdir << ">" << endl;
	  // Maybe we need to create it?
	  cout << "---- I will attempt to create it" << endl;
	  if (mkdir(outdir.c_str(), 0755) == -1) {
	    cout << "---- Error creating directory <" << outdir
		 << ">, aborting" << endl;
	    perror("mkdir");
	    nOK = 1;
	  }
	} else {
	  cout << "---- Output path <" << outdir << "> is "; // 
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
      cout << "---- ";
      if (ok) 
	cout << "Test file opened and deleted in output directory, good!";
      else    
	cout << "We can't open files in this path, attempted to create: " 
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
  std::ofstream out;
  int nOK = 0;

  if (myid==0) {
    string curparm(outdir + parmfile + "." + runtag + ".yml");
    out.open(curparm.c_str());
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

void add_others()
{
  // Blank nodes at this point (could add some example values...)
  //
  YAML::Node comp, outp, extr, intr;

  // Add other main map elements
  //
  parse["Components"]  = comp;
  parse["Output"]      = outp;
  parse["External"]    = extr;
  parse["Interaction"] = intr;

}

void print_default(void)
{
  update_parm();
  add_others();
  
  std::cout << "# EXP [" << VERSION << "]" << std::endl
	    << parse << std::endl;
}


void YAML_parse_args(int argc, char** argv)
{
  // Program name (should be 'exp')
  //
  char *prog = argv[0];
  
  // Default config file if none specified using -f or as a trailing option
  //
  std::string curparm("config.yml");

  // Get our MPI id
  int myid;
  MPI_Comm_rank(MPI_COMM_WORLD, &myid);
  
  int done = 0;

  if (myid==0) {

    cxxopts::Options options(prog, "Basis-function EXPansion code\n");

    // Allow unmatched options so that we can get the trailing config
    // file name, if specified
    //
    options.allow_unrecognised_options();

    options.add_options()
      ("h,help", "this help message")
      ("f,file", "the input YAML configuration file",
       cxxopts::value<std::string>(curparm))
      ("c,config", "the input YAML configuration file",
       cxxopts::value<std::string>(curparm))
      ("t,template", "provide a template YAML configuration file")
      ("v,git", "display verbose GIT version info")
      ;
    
    cxxopts::ParseResult vm;

    try {
      vm = options.parse(argc, argv);
    } catch (cxxopts::OptionException& e) {
      std::cout << "Option error: " << e.what() << std::endl;
      MPI_Finalize();
      exit(-1);
    }

    if (vm.count("help")) {
      std::cout << options.help() << std::endl
		<< "* The YAML config may be appended to the command line without flags" << std::endl
		<< "* See EXP/doc/html/index.html for extensive documentation" << std::endl
		<< std::endl << std::endl;
      done = 1;
    }
    
    if (vm.count("git")) {
      // Version is printed automatically, so we only need to exit
      done = 1;
    }
    
    if (vm.count("template")) {
      // Print a template config file with default parameter values
      print_default();
      done = 1;
    }

    // Termination check #1
    //
    MPI_Bcast(&done, 1, MPI_INT, 0, MPI_COMM_WORLD);

    if (done) {
      MPI_Finalize();
      exit(EXIT_SUCCESS);
    }

    // If config file is not specified, use first trailing, unmatched
    // argument as config file
    //
    bool trailing = false;
    if (vm.count("file")==0 and vm.count("config")==0) {
      auto trail = vm.unmatched();
      if (trail.size()) {
	curparm = trail[0];
	trailing = true;
      }
    }

    // Now attempt to read the YAML config file
    //
    std::ifstream in(curparm);
    if (!in) {
      // Failure to open YAML file for reading
      //
      if (myid==0) {
	std::cerr << argv[0] << ": can not open parameter file <" << curparm
		  << ">" << std::endl;

	if (trailing)		// Additional info for user debugging...
	  std::cerr << argv[0]
		    << ": the YAML config file name was the first "
		    << "unmatched command-line argument" << std::endl;

	char mbuf[512];
	std::cerr << argv[0] << ": pwd is <" << getcwd(mbuf, 512) << ">"
		  << std::endl;
      }
      done = 1;
    }
    
    // Termination check #2
    //
    MPI_Bcast(&done, 1, MPI_INT, 0, MPI_COMM_WORLD);
    
    if (done) {
      MPI_Finalize();
      exit(EXIT_SUCCESS);
    }
    
    try {
      parse = YAML::Load(in);
    }
    catch (YAML::Exception & error) {
      std::cout << "Error parsing config file: "
		<< error.what() << std::endl;
      done = 1;
    }

    // Termination check #3
    //
    MPI_Bcast(&done, 1, MPI_INT, 0, MPI_COMM_WORLD);

    // The current YAML structure will now be serialized and broadcast
    // to all processes
    //
    std::ostringstream serial;
    serial << parse << std::endl;

    int line_size = serial.str().size();
    MPI_Bcast(&line_size, 1, MPI_INT, 0, MPI_COMM_WORLD);
    MPI_Bcast(const_cast<char*>(serial.str().data()), line_size, MPI_CHAR, 0, MPI_COMM_WORLD);
    
  } else {
    
    // 3 termination checks
    //
    for (int i=0; i<3; i++) {
      MPI_Bcast(&done, 1, MPI_INT, 0, MPI_COMM_WORLD);

      if (done) {
	MPI_Finalize();
	exit(EXIT_SUCCESS);
      }
    }

    // Receive the YAML serialized length
    //
    int line_size;
    MPI_Bcast(&line_size, 1, MPI_INT, 0, MPI_COMM_WORLD);
    
    std::string config;
    config.resize(line_size);

    // Receive the YAML serialized string
    //
    MPI_Bcast(const_cast<char*>(config.data()), line_size, MPI_CHAR, 0, MPI_COMM_WORLD);
    
    std::istringstream sin(config);
    
    // Load the YAML structure
    //
    parse = YAML::Load(sin);
  }

  initialize();
}
