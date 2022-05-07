#include <iostream>
#include <fstream>
#include <iomanip>
#include <sstream>
#include <chrono>

#include "expand.H"
#include <global.H>

#include <AxisymmetricBasis.H>
#include <OutPSR.H>

OutPSR::OutPSR(const YAML::Node& conf) : Output(conf)
{
  initialize();
}

void OutPSR::initialize()
{
  try {
				// Get file name
    if (Output::conf["filename"])
      filename = Output::conf["filename"].as<std::string>();
    else{
      filename.erase();
      filename = "SPL." + runtag;
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
      nintsub_warning("OutPSR");
      nintsub = std::numeric_limits<int>::max();
#endif
    } else
      nintsub = std::numeric_limits<int>::max();

    if (Output::conf["nbeg"])
      nbeg = Output::conf["nbeg"].as<int>();
    else
      nbeg = 0;

    if (Output::conf["real4"])
      real4 = Output::conf["real4"].as<bool>();
    else
      real4 = true;

    if (Output::conf["timer"])
      timer = Output::conf["timer"].as<bool>();
    else
      timer = false;

    if (Output::conf["threads"])
      nth = Output::conf["threads"].as<int>();
    else
      nth = 1;
  }
  catch (YAML::Exception & error) {
    if (myid==0) std::cout << "Error parsing parameters in OutPSR: "
			   << error.what() << std::endl
			   << std::string(60, '-') << std::endl
			   << "Config node"        << std::endl
			   << std::string(60, '-') << std::endl
			   << conf                 << std::endl
			   << std::string(60, '-') << std::endl;
    MPI_Finalize();
    exit(-1);
  }


  // Determine last file
  // 
  if (restart && nbeg==0) {

    // Only root node looks for files
    //
    if (myid==0) {

      for (nbeg=0; nbeg<100000; nbeg++) {
	// Output name
	//
	ostringstream fname;
	fname << outdir << filename << "." << setw(5) << setfill('0') << nbeg;

	// See if we can open file
	//
	ifstream in(fname.str().c_str());

	if (!in) {
	  cout << "OutPSR: will begin with nbeg=" << nbeg << endl;
	  break;
	}
      }
    }

    // Communicate starting file to all nodes
    //
    MPI_Bcast(&nbeg, 1, MPI_INT, 0, MPI_COMM_WORLD);
  }
}


void OutPSR::Run(int n, int mstep, bool last)
{
  if (!dump_signal and !last) {
    if (n % nint            ) return;
    if (restart  && n==0    ) return;
    if (mstep % nintsub !=0 ) return;
  }

  std::chrono::high_resolution_clock::time_point beg, end;
  if (timer) beg = std::chrono::high_resolution_clock::now();
  
  std::ofstream out;
  std::ostringstream fname;

  // Output name prefix
  fname << filename << "." << setw(5) << setfill('0') << nbeg++;

  // Master file name
  std::string master = outdir + fname.str();

  psdump = n;

  int nOK = 0;

  if (myid==0) {
				// Open file and write master header
    out.open(master);

    if (out.fail()) {
      std::cerr << "OutPSR: can't open file <" << master
		<< "> . . . quitting" << std::endl;
      nOK = 1;
    }
				// Used by OutCHKPT to not duplicate a dump
    if (not real4) lastPSR = fname.str();
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

				// Write component header
    std::ostringstream cname;
    cname << fname.str() << "_" << count++;
    
    if (myid==0) {
      c->write_binary_header(&out, real4, cname.str(), nth);
    }

				// Stream and file name for each
				// thread
    std::vector<std::shared_ptr<std::ofstream>> tstreams(nth);
    std::vector<std::string > blobname(nth);
				// Create the streams
    for (int n=0; n<nth; n++) {
      std::ostringstream sname;
      sname << outdir << cname.str() << "-" << myid*nth + n;
				// Open particle file and write
      blobname[n] = sname.str();
      tstreams[n] = std::make_shared<std::ofstream>(sname.str());

      if (tstreams[n]->fail()) { // Check for good stream
	std::cerr << "[" << myid << ", " << n << "] OutPSR: can't open file <" << blobname[n]
		  << "> . . . quitting" << std::endl;
	nOK = 1;
      }
    }
    
				// Write particles into the streams
    if (nOK==0) {
      c->write_binary_particles(tstreams, real4);

      for (int n=0; n<nth; n++) {
	if (tstreams[n]->fail()) {
	  std::cout << "OutPSR: error writing binary particles to <" << blobname[n]
		    << std::endl;
	  nOK = 1;
	}

	try {
	  tstreams[n]->close();
	}
	catch (const ofstream::failure& e) {
	  std::cout << "OutPSR: exception closing file <" << blobname[n]
		    << ": " << e.what() << std::endl;
	  nOK = 1;
	}
      }
    }
    
				// Check for errors in all file opening
    int sumOK;
    MPI_Allreduce(&nOK, &sumOK, 1, MPI_INT, MPI_SUM, MPI_COMM_WORLD);

    if (sumOK) {
      MPI_Finalize();
      exit(39);
    }
  }

  if (myid==0) {
    if (out.fail()) {
      std::cout << "OutPSR: error writing component to master <" << master
		<< std::endl;
    }

    try {
      out.close();
    }
    catch (const ofstream::failure& e) {
      std::cout << "OutPSR: exception closing file <" << master
		<< ": " << e.what() << std::endl;
    }
  }

  chktimer.mark();

  dump_signal = 0;

  if (timer) {
    end = std::chrono::high_resolution_clock::now();
    std::chrono::duration<double> intvl = end - beg;
    if (myid==0)
      std::cout << "OutPSR [T=" << tnow << "] timing=" << intvl.count()
		<< std::endl;
  }
}

