#include <cstdio>
#include <iostream>
#include <fstream>
#include <iomanip>
#include <sstream>
#include <chrono>

#include <unistd.h>		// For unlink

#include "expand.H"
#include <global.H>

#include <AxisymmetricBasis.H>
#include <OutCHKPT.H>

const std::set<std::string> OutCHKPT::valid_keys = {
  "mpio",
  "filename",
  "nint",
  "nintsub",
  "timer",
  "nagg"
};

OutCHKPT::OutCHKPT(const YAML::Node& conf) : Output(conf)
{
  initialize();
}

void OutCHKPT::initialize()
{
  // Remove matched keys
  for (auto v : valid_keys) current_keys.erase(v);
  
  // Assign values from YAML
  //
  try {
    if (Output::conf["mpio"])
      mpio = Output::conf["mpio"].as<bool>();
    else
      mpio = false;
				// Get file name
    if (Output::conf["filename"])
      filename = Output::conf["filename"].as<std::string>();
    else {
      filename.erase();
      filename = outdir + "OUT." + runtag + ".chkpt";
    }
    
    if (Output::conf["nint"])
      nint = Output::conf["nint"].as<int>();
    else
      nint = 100;
    
    if (Output::conf["nintsub"]) {
#ifdef ALLOW_NINTSUB
      nintsub = Output::conf["nintsub"].as<int>();
#else
      nintsub_warning("OutCHKPT");
      nintsub = std::numeric_limits<int>::max();
#endif
    } else
      nintsub = std::numeric_limits<int>::max();

				// Sanity check
    if (nintsub <= 0) nintsub = 1;

    if (Output::conf["timer"])
      timer = Output::conf["timer"].as<bool>();
    else
      timer = false;

    if (Output::conf["nagg"])
      nagg = Output::conf["nagg"].as<std::string>();
    else
      nagg = "1";
  }
  catch (YAML::Exception & error) {
    if (myid==0) std::cout << "Error parsing parameters in OutCHKPT: "
			   << error.what() << std::endl
			   << std::string(60, '-') << std::endl
			   << "Config node"        << std::endl
			   << std::string(60, '-') << std::endl
			   << Output::conf         << std::endl
			   << std::string(60, '-') << std::endl;
    MPI_Finalize();
    exit(-1);
  }
}


void OutCHKPT::Run(int n, int mstep, bool last)
{
  if (!dump_signal and !last) {
    if (n % nint) return;
    if (multistep>1 and mstep % nintsub !=0) return;
  }

  if (VERBOSE>5 && myid==0) {
    cout << " OutCHKPT::Run(): n=" << n << " psdump=" << psdump << endl;
  }
  
  int returnStatus = 1;

  if (n == psdump) {       

    if (myid==0) {
      string backfile = filename + ".bak";
      if (unlink(backfile.c_str())) {
	if (VERBOSE>5) perror("OutCHKPT::Run()");
	cout << "OutCHKPT::Run(): error unlinking old backup file <" 
	     << backfile << ">, it may not exist" << endl;
      } else {
	if (VERBOSE>5) {
	  cout << "OutCHKPT::Run(): successfully unlinked <"
	       << backfile << ">" << endl;
	}
      }
      if (rename(filename.c_str(), backfile.c_str())) {
	if (VERBOSE>5) perror("OutCHKPT::Run()");
	cout << "OutCHKPT: renaming backup file <" 
	     << backfile << ">, it may not exist" << endl;
      } else {
	if (VERBOSE>5) {
	  cout << "OutCHKPT::Run(): successfully renamed <"
	       << filename << "> to <" << backfile << ">" << endl;
	}
      }

      if (lastPS.size()) {
	if (symlink(lastPS.c_str(), filename.c_str())) {
	  if (VERBOSE>5) perror("OutCHKPT::Run()");
	  cout << "OutCHKPT::Run(): no file <" << lastPS
	       << "> to link, we will create a new checkpoint" << endl;
	  returnStatus = 0;
	} else {
	  if (VERBOSE>5) {
	    cout << "OutCHKPT::Run(): successfully linked <"
		 << lastPS << "> to new backup file <" 
		 << filename << ">" << endl;
	  }
	}
      } else {
	returnStatus = 0;
      }
    }
  }

  MPI_Bcast(&returnStatus, 1, MPI_INT, 0, MPI_COMM_WORLD);
  if (returnStatus==1) return;

  std::chrono::high_resolution_clock::time_point beg, end;
  if (timer) beg = std::chrono::high_resolution_clock::now();
  
  if (mpio) {
    static bool firsttime = true;

    // MPI variables
    //
    char err[MPI_MAX_ERROR_STRING];
    MPI_Offset offset = 0;
    MPI_Status status;
    MPI_Info   info;
    MPI_File   file;
    int        len;
    int        nOK = 0;

    // Return info about errors (for debugging)
    //
    MPI_Comm_set_errhandler(MPI_COMM_WORLD, MPI_ERRORS_RETURN); 

    // Set info to limit the number of aggregators
    //
    MPI_Info_create(&info);
    MPI_Info_set(info, "cb_nodes", nagg.c_str());
    
    // Open shared file
    //
    int ret =
      MPI_File_open(MPI_COMM_WORLD, filename.c_str(),
		    MPI_MODE_CREATE | MPI_MODE_WRONLY | MPI_MODE_UNIQUE_OPEN,
		    info, &file);
    
    if (ret != MPI_SUCCESS) {
      std::cerr << "OutCHKPT:run: rank [" << myid << "] can't open file <"
		<< filename << "> . . . quitting" << std::endl;
      nOK = 1;
    }
    
    int badCount = 0;
    MPI_Allreduce(&nOK, &badCount, 1, MPI_INT, MPI_SUM, MPI_COMM_WORLD);

    if (badCount) {
      MPI_Finalize();
      exit(33);
    }

    MPI_Info_free(&info);
    
    // Write master header
    //
    if (myid==0) {
      struct MasterHeader header;
      header.time  = tnow;
      header.ntot  = comp->ntot;
      header.ncomp = comp->ncomp;
      
      ret = MPI_File_write_at(file, offset, &header, sizeof(MasterHeader),
			      MPI_CHAR, &status);

      if (ret != MPI_SUCCESS) {
	MPI_Error_string(ret, err, &len);
	std::cout << "OutCHKPT::run: WRITE header " << err
		  << " at line " << __LINE__ << std::endl;
      }
    }
  
    offset += sizeof(MasterHeader);

    for (auto c : comp->components) {
      if (firsttime and myid==0 and not c->Indexing())
	std::cout << "OutCHKPT::run: component <" << c->name
		  << "> has not set 'indexing' so PSP particle sequence will be lost." << std::endl
		  << "If this is NOT what you want, set the component flag 'indexing=1'." << std::endl;

#ifdef HAVE_LIBCUDA
    if (use_cuda) {
      if (c->force->cudaAware() and not comp->fetched[c]) {
	comp->fetched[c] = true;
	c->CudaToParticles();
      }
    }
#endif
      c->write_binary_mpi(file, offset); 
    }
    
    // Try SYNC-BARRIER-SYNC semantic
    //
    ret = MPI_File_sync(file);

    if (ret != MPI_SUCCESS) {
      MPI_Error_string(ret, err, &len);
      std::cout << "OutCHKPT::run: SYNC " << err
		<< " at line " << __LINE__ << std::endl;
    }

    MPI_Barrier(MPI_COMM_WORLD);

    ret = MPI_File_sync(file);

    if (ret != MPI_SUCCESS) {
      MPI_Error_string(ret, err, &len);
      std::cout << "OutCHKPT::run: SYNC " << err
		<< " at line " << __LINE__ << std::endl;
    }

    ret = MPI_File_close(&file);

    if (ret != MPI_SUCCESS) {
      MPI_Error_string(ret, err, &len);
      std::cout << "OutCHKPT::run: CLOSE " << err
		<< " at line " << __LINE__ << std::endl;
    }

  } else {

    std::ofstream out;
    int nOK = 0;

    if (myid==0) {
				// Open file and write master header
      out.open(filename);

      if (out.fail()) {
	std::cerr << "OutCHKPT: can't open file <" << filename
	     << "> . . . quitting\n";
	nOK = 1;
      }
				// Open file and write master header
      if (!nOK) {
	struct MasterHeader header;
	header.time  = tnow;
	header.ntot  = comp->ntot;
	header.ncomp = comp->ncomp;
      
	out.write((char *)&header, sizeof(MasterHeader));
#ifdef DEBUG
	cout << "OutCHKPT: header written" << endl;
#endif
      }
    }
    
    MPI_Bcast(&nOK, 1, MPI_INT, 0, MPI_COMM_WORLD);
    if (nOK) {
      MPI_Finalize();
      exit(33);
    }

    for (auto c : comp->components) {
#ifdef DEBUG
      cout << "OutCHKPT: process " << myid << " trying to write name=" << c->name
	   << " force=" << c->id << endl;
#endif
      c->write_binary(&out);
#ifdef DEBUG
      cout << "OutCHKPT: process " << myid << " write completed on " << c->name << endl;
#endif
    }
    
    if (myid==0) {
      try {
	out.close();
      }
      catch (const ofstream::failure& e) {
	std::cout << "OutCHKPT: exception closing file <" << filename
		  << ": " << e.what() << std::endl;
      }
    }

  }

  chktimer.mark();

  dump_signal = 0;

  if (timer) {
    end = std::chrono::high_resolution_clock::now();
    std::chrono::duration<double> intvl = end - beg;
    if (myid==0)
      std::cout << "OutCHKPT [T=" << tnow << "] timing=" << intvl.count()
		<< std::endl;
  }
}

