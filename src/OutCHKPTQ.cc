#include <cstdio>
#include <iostream>
#include <fstream>
#include <iomanip>
#include <sstream>
#include <chrono>

#include <unistd.h>		// for unlink

#include "expand.H"
#include <global.H>

#include <AxisymmetricBasis.H>
#include <OutCHKPTQ.H>


const std::set<std::string> OutCHKPTQ::valid_keys = {
  "mpio",
  "filename",
  "nint",
  "nintsub",
  "timer",
};


OutCHKPTQ::OutCHKPTQ(const YAML::Node& conf) : Output(conf)
{
  initialize();
}

void OutCHKPTQ::initialize()
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
      filename = "SPL." + runtag + ".chkpt";
    }
    
    if (Output::conf["nint"])
      nint = Output::conf["nint"].as<int>();
    else
      nint = 100;

    if (Output::conf["nintsub"]) {
#ifdef ALLOW_NINTSUB
      nintsub = Output::conf["nintsub"].as<int>();
      if (nintsub <= 0) nintsub = 1;
#else
      nintsub_warning("OutCHKPTQ");
      nintsub = std::numeric_limits<int>::max();
#endif
    } else
      nintsub  = std::numeric_limits<int>::max();

    if (Output::conf["timer"])
      timer = Output::conf["timer"].as<bool>();
    else
      timer = false;
  }
  catch (YAML::Exception & error) {
    if (myid==0) std::cout << "Error parsing parameters in OutCHKPTQ: "
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

void OutCHKPTQ::Run(int n, int mstep, bool last)
{
  if (!dump_signal and !last) {
    if (n % nint) return;
    if (multistep>1 and mstep % nintsub !=0) return;
  }
  
  if (VERBOSE>5 && myid==0) {
    cout << " OutCHKPTQ::Run(): n=" << n << " psdump=" << psdump << endl;
  }

  int returnStatus = 1;
  
  if (n == psdump) {       

    if (myid==0) {
      string currfile = outdir + filename;
      string backfile = currfile + ".bak";
      if (unlink(backfile.c_str())) {
	if (VERBOSE>5) perror("OutCHKPTQ::Run()");
	cout << "OutCHKPTQ::Run(): error unlinking old backup file <" 
	     << backfile << ">, it may not exist" << endl;
      } else {
	if (VERBOSE>5) {
	  cout << "OutCHKPTQ::Run(): successfully unlinked <"
	       << backfile << ">" << endl;
	}
      }

      if (rename(currfile.c_str(), backfile.c_str())) {
	if (VERBOSE>5) perror("OutCHKPTQ::Run()");
	cout << "OutCHKPTQ: renaming backup file <" 
	     << backfile << ">, it may not exist" << endl;
      } else {
	if (VERBOSE>5) {
	  cout << "OutCHKPTQ::Run(): successfully renamed <"
	       << currfile << "> to <" << backfile << ">" << endl;
	}
      }

      if (lastPSQ.size()) {
	if (symlink(lastPSQ.c_str(), currfile.c_str())) {
	  if (VERBOSE>5) perror("OutCHKPTQ::Run()");
	  cout << "OutCHKPTQ::Run(): no file <" << lastPSQ
	       << "> to link, we will create a new checkpoint" << endl;
	  returnStatus = 0;
	} else {
	  if (VERBOSE>5) {
	    cout << "OutCHKPTQ::Run(): successfully linked <"
		 << lastPSQ << "> to new backup file <" 
		 << currfile << ">" << endl;
	  }
	}
      } else {
	returnStatus = 0;
      }

      int count = 0;
      for (auto c : comp->components) {
				// Component file
	std::ostringstream cname;
	cname << outdir << filename << "_" << count;
    
	for (int n=0; n<numprocs; n++) {

	  std::ostringstream sout1;
	  sout1 << cname.str() << "-" << n;
	  std::string fileN = sout1.str();

	  sout1 << ".bak";
	  std::string backN = sout1.str();

	  if (unlink(backN.c_str())) {
	    if (VERBOSE>15) perror("OutCHKPTQ::Run()");
	    if (VERBOSE>10)
	      cout << "OutCHKPTQ::Run(): error unlinking old backup file <" 
		   << backN << ">, it may not exist" << endl;
	  } else {
	    if (VERBOSE>15) {
	      cout << "OutCHKPTQ::Run(): successfully unlinked <"
		   << backN << ">" << endl;
	    }
	  }
	  if (rename(fileN.c_str(), backN.c_str())) {
	    if (VERBOSE>15) perror("OutCHKPTQ::Run()");
	    if (VERBOSE>10)
	      cout << "OutCHKPTQ: renaming backup file <" 
		   << backN << ">, it may not exist" << endl;
	  } else {
	    if (VERBOSE>15) {
	      cout << "OutCHKPTQ::Run(): successfully renamed <"
		   << filename << "> to <" << backN << ">" << endl;
	    }
	  }

	  if (lastPSQ.size()) {
	    std::ostringstream sout2;
	    sout2 << lastPSQ << "_" << count << "-" << n;
	    std::string lastN = sout2.str();
	  
	    if (symlink(lastN.c_str(), fileN.c_str())) {
	      if (VERBOSE>15) perror("OutCHKPTQ::Run()");
	      if (VERBOSE>10)
		cout << "OutCHKPTQ::Run(): no file <" << lastN
		     << "> to link, we will create a new checkpoint" << endl;
	      returnStatus = 0;
	    } else {
	      if (VERBOSE>15) {
		cout << "OutCHKPTQ::Run(): successfully linked <"
		     << lastN << "> to new backup file <" 
		     << fileN << ">" << endl;
	      }
	    }
	  } else {
	    returnStatus = 0;
	  }
	} // END: process loop

	count++;
      } // END: component loop

    } // myid==0

  } // n == psdump

  MPI_Bcast(&returnStatus, 1, MPI_INT, 0, MPI_COMM_WORLD);
  if (returnStatus==1) return;

  std::chrono::high_resolution_clock::time_point beg, end;
  if (timer) beg = std::chrono::high_resolution_clock::now();
  
  psdump = n;

  int nOK = 0;

  std::ofstream out;

  if (myid==0) {
				// Open file and write master header
    std::string master = outdir + filename;
    out.open(master);

    if (out.fail()) {
      std::cerr << "OutCHKPTQ: can't open file <" << master
		<< "> . . . quitting" << std::endl;
      nOK = 1;
    }
				// Open file and write master header
    if (nOK==0) {
      struct MasterHeader header;
      header.time  = tnow;
      header.ntot  = comp->ntot;
      header.ncomp = comp->ncomp;

      out.write((char *)&header, sizeof(MasterHeader));
    }
  }
  
  MPI_Bcast(&nOK, 1, MPI_INT, 0, MPI_COMM_WORLD);
  if (nOK) {
    MPI_Finalize();
    exit(33);
  }

  int count = 0;
  for (auto c : comp->components) {
#ifdef HAVE_LIBCUDA
    if (use_cuda) {
      if (c->force->cudaAware() and not comp->fetched[c]) {
	comp->fetched[c] = true;
	c->CudaToParticles();
      }
    }
#endif
				// Check for open failures
    nOK = 0;

				// Component file
    std::ostringstream cname;
    cname << filename << "_" << count++;
    
    if (myid==0) {
      c->write_binary_header(&out, false, cname.str());
    }

    cname << "-" << myid;
    
				// Open file and write component header
    std::string blobfile = outdir + cname.str();
    std::ofstream pout(blobfile);

    if (pout.fail()) {
      std::cerr << "[" << myid << "] OutCHKPTQ: can't open file <" << blobfile 
		<< "> . . . quitting" << std::endl;
      nOK = 1;
    } else {
      c->write_binary_particles(&pout, false);
      if (pout.fail()) {
	std::cout << "OutCHKPTQ: error writing binary particles to <"
		  << cname.str() << std::endl;
      }
    }

    
				// Check for errors in all file opening
    int sumOK;
    MPI_Allreduce(&nOK, &sumOK, 1, MPI_INT, MPI_SUM, MPI_COMM_WORLD);

    if (sumOK) {
      MPI_Finalize();
      exit(35);
    }
  }

  if (myid==0) {
    if (out.fail()) {
      std::cout << "OutCHKPTQ: error writing component to master <"
		<< outdir + filename << std::endl;
    }

    try {
      out.close();
    }
    catch (const ofstream::failure& e) {
      std::cout << "OutCHKPTQ: exception closing file <" << outdir + filename
		<< ": " << e.what() << std::endl;
    }
  }

  chktimer.mark();

  dump_signal = 0;

  if (timer) {
    end = std::chrono::high_resolution_clock::now();
    std::chrono::duration<double> intvl = end - beg;
    if (myid==0)
      std::cout << "OutCHKPTQ [T=" << tnow << "] timing=" << intvl.count()
		<< std::endl;
  }
}

