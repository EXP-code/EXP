#include <iostream>
#include <fstream>
#include <iomanip>
#include <sstream>
#include <chrono>

#include "expand.H"
#include "global.H"

#include "AxisymmetricBasis.H"
#include "OutPSP.H"

const std::set<std::string>
OutPSP::valid_keys = {
  "filename",
  "nint",
  "nintsub",
  "nbeg",
  "real4",
  "timer",
  "nagg"
};


OutPSP::OutPSP(const YAML::Node& conf) : Output(conf)
{
  initialize();
}

void OutPSP::initialize()
{
  // Remove matched keys
  //
  for (auto v : valid_keys) current_keys.erase(v);
  
  // Assign values from YAML
  //
  try {
				// Get file name
    if (Output::conf["filename"])
      filename = Output::conf["filename"].as<std::string>();
    else {
      filename.erase();
      filename = outdir + "OUT." + runtag;
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
      nintsub_warning("OutPSN");
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

    if (Output::conf["nagg"])
      nagg = Output::conf["nagg"].as<std::string>();
    else
      nagg = "1";
  }
  catch (YAML::Exception & error) {
    if (myid==0) std::cout << "Error parsing parameters in OutPSP: "
			   << error.what() << std::endl
			   << std::string(60, '-') << std::endl
			   << "Config node"        << std::endl
			   << std::string(60, '-') << std::endl
			   << conf                 << std::endl
			   << std::string(60, '-') << std::endl;
    throw std::runtime_error("OutPSP::initialize: error parsing YAML");
  }

				// Determine last file
  if (restart && nbeg==0) {
    if (myid==0) {

      for (nbeg=0; nbeg<100000; nbeg++) {

				// Output name
	ostringstream fname;
	fname << filename << "." << setw(5) << setfill('0') << nbeg;

				// See if we can open file
	ifstream in(fname.str().c_str());
	
	if (!in) {
	  cout << "OutPSP: will begin with nbeg=" << nbeg << endl;
	  break;
	}
      }
				// All nodes need nbeg for MPI_File_open
      MPI_Bcast(&nbeg, 1, MPI_INT, 0, MPI_COMM_WORLD);
    } else {			// Get from root via Bcast
      MPI_Bcast(&nbeg, 1, MPI_INT, 0, MPI_COMM_WORLD);
    }
  }

}


void OutPSP::Run(int n, int mstep, bool last)
{
  if (!dump_signal and !last) {
    if (n % nint) return;
    if (restart && n==0) return;
    if (multistep>1 && mstep % nintsub !=0 ) return;
  }

  std::chrono::high_resolution_clock::time_point beg, end;
  if (timer) beg = std::chrono::high_resolution_clock::now();
  
  static bool firsttime = true;

  char err[MPI_MAX_ERROR_STRING];
  MPI_Offset offset = 0;
  MPI_Status status;
  MPI_File   file;
  MPI_Info   info;
  int        len;

  // Output name
  //
  ostringstream fname;
  fname << filename << "." << setw(5) << setfill('0') << nbeg++;

  // return info about errors (for debugging)
  MPI_Comm_set_errhandler(MPI_COMM_WORLD, MPI_ERRORS_RETURN); 

  // Set info to limit the number of aggregators
  //
  MPI_Info_create(&info);
  MPI_Info_set(info, "cb_nodes", nagg.c_str());

  // Open shared file and write master header
  //
  int ret =
    MPI_File_open(MPI_COMM_WORLD, fname.str().c_str(),
		  MPI_MODE_CREATE | MPI_MODE_WRONLY | MPI_MODE_UNIQUE_OPEN,
		  info, &file);

  if (ret != MPI_SUCCESS) {
    std::ostringstream sout;
    sout << "OutPSP: can't open file <" << fname.str() << "> . . . quitting";
    throw GenericError(sout.str(), __FILE__, __LINE__, 33, false);
  }

  MPI_Info_free(&info);

				// Used by OutCHKPT to not duplicate a dump
  if (!real4) lastPS = fname.str();
  
    
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
      std::cout << "OutPSP::run: WRITE header " << err
		<< " at line " << __LINE__ << std::endl;
    }
  }
  
  offset += sizeof(MasterHeader);

  for (auto c : comp->components) {

#ifdef HAVE_LIBCUDA
    if (use_cuda) {
      if (not comp->fetched[c]) {
	comp->fetched[c] = true;
	c->CudaToParticles();
      }
    }
#endif

    if (firsttime and myid==0 and not c->Indexing())
      std::cout << "OutPSP::run: component <" << c->name
		<< "> has not set 'indexing' so PSP particle sequence will be lost." << std::endl
		<< "If this is NOT what you want, set the component flag 'indexing=1'." << std::endl;
    c->write_binary_mpi(file, offset, real4); 
  }

  // Try SYNC-BARRIER-SYNC semantic
  //
  ret = MPI_File_sync(file);

  if (ret != MPI_SUCCESS) {
    MPI_Error_string(ret, err, &len);
    std::cout << "OutPSP::run: SYNC " << err
	      << " at line " << __LINE__ << std::endl;
  }
  
  MPI_Barrier(MPI_COMM_WORLD);

  ret = MPI_File_sync(file);

  if (ret != MPI_SUCCESS) {
    MPI_Error_string(ret, err, &len);
    std::cout << "OutPSP::run: SYNC " << err
	      << " at line " << __LINE__ << std::endl;
  }

  ret = MPI_File_close(&file);

  firsttime = false;

  if (ret != MPI_SUCCESS) {
    MPI_Error_string(ret, err, &len);
    std::cout << "OutPSP::run: CLOSE " << err
	      << " at line " << __LINE__ << std::endl;
  }

  chktimer.mark();

  dump_signal = 0;

  if (timer) {
    end = std::chrono::high_resolution_clock::now();
    std::chrono::duration<double> intvl = end - beg;
    if (myid==0)
      std::cout << "OutPSP [T=" << tnow << "] timing=" << intvl.count()
	      << std::endl;
  }
}
