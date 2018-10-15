#include <cstdio>
#include <iostream>
#include <fstream>
#include <iomanip>
#include <sstream>
#include <chrono>

#include "expand.h"
#include <global.H>

#include <AxisymmetricBasis.H>
#include <OutCHKPT.H>


OutCHKPT::OutCHKPT(string& line) : Output(line)
{
  initialize();
}

void OutCHKPT::initialize()
{
  std::string tmp;

  if (Output::get_value(string("mpio"), tmp))
    mpio = atoi(tmp.c_str()) ? true : false;
  else
    mpio = false;
				// Get file name
  if (!Output::get_value(string("filename"), filename)) {
    filename.erase();
    if (mpio)
      filename = outdir + "OUTS." + runtag + ".chkpt";
    else
      filename = outdir + "OUT." + runtag + ".chkpt";
  }

  if (Output::get_value(string("nint"), tmp))
    nint = atoi(tmp.c_str());
  else
    nint = 100;

  if (Output::get_value(string("timer"), tmp))
    timer = atoi(tmp.c_str()) ? true : false;
  else
    timer = false;

  if (Output::get_value(string("nagg"), tmp))
    nagg = tmp;
  else
    nagg = "1";
}


void OutCHKPT::Run(int n, bool last)
{
  if (n % nint && !last) return;
  if (VERBOSE>5 && myid==0) {
    cout << " OutCHKPT::Run(): n=" << n << " psdump=" << psdump << endl;
  }

  if (n == psdump) {       

    if (myid==0) {
      string backfile = filename + ".bak";
      if (unlink(backfile.c_str())) {
	perror("OutCHKPT::Run()");
	cout << "OutCHKPT::Run(): error unlinking old backup file <" 
	     << backfile << ">" << endl;
      } else {
	if (VERBOSE>5) {
	  cout << "OutCHKPT::Run(): successfully unlinked <"
	       << backfile << ">" << endl;
	}
      }
      if (rename(filename.c_str(), backfile.c_str())) {
	perror("OutCHKPT::Run()");
	cout << "OutCHKPT: error creating backup file <" 
	     << backfile << ">" << endl;
      } else {
	if (VERBOSE>5) {
	  cout << "OutCHKPT::Run(): successfully renamed <"
	       << filename << "> to <" << backfile << ">" << endl;
	}
      }
      if (symlink(lastPS.c_str(), filename.c_str())) {
	perror("OutCHKPT::Run()");
	cout << "OutCHKPT::Run(): error symlinking new backup file <" 
	     << filename << ">" << endl;
      } else {
	if (VERBOSE>5) {
	  cout << "OutCHKPT::Run(): successfully linked <"
	       << lastPS << "> to new backup file <" 
	       << filename << ">" << endl;
	}
      }
    }
    return;
  }

  std::chrono::high_resolution_clock::time_point beg, end;
  if (timer) beg = std::chrono::high_resolution_clock::now();
  
  if (myid==0) {
				// Attempt to move file to backup
    string backfile = filename + ".bak";
    if (rename(filename.c_str(), backfile.c_str())) {
      perror("OutCHKPT::Run()");
      cout << "OutCHKPT::Run(): error creating backup file <" 
	   << backfile << ">";
    }
  }

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
      cerr << "OutCHKPT:run: can't open file <" << filename << "> . . . quitting"
	   << std::endl;
      MPI_Abort(MPI_COMM_WORLD, 33);
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
	std::cout << "OutCHKPT::run: " << err
		  << " at line " << __LINE__ << std::endl;
      }
    }
  
    offset += sizeof(MasterHeader);

    for (auto c : comp->components) {
      if (firsttime and myid==0 and not c->Indexing())
	std::cout << "OutCHKPT::run: component <" << c->name
		  << "> has not set 'indexing' so PSP particle sequence will be lost." << std::endl
		  << "If this is NOT what you want, set the component flag 'indexing=1'." << std::endl;
      c->write_binary_mpi(file, offset); 
    }
    
    ret = MPI_File_close(&file);

    if (ret != MPI_SUCCESS) {
      MPI_Error_string(ret, err, &len);
      std::cout << "OutCHKPT::run: " << err
		<< " at line " << __LINE__ << std::endl;
    }

  } else {

    ofstream *out;

    if (myid==0) {
				// Open file and write master header
      out = new ofstream(filename.c_str());

      if (!*out) {
	cerr << "OutCHKPT: can't open file <" << filename.c_str() 
	     << "> . . . quitting\n";
	MPI_Abort(MPI_COMM_WORLD, 33);
      }
				// Open file and write master header
    
      struct MasterHeader header;
      header.time  = tnow;
      header.ntot  = comp->ntot;
      header.ncomp = comp->ncomp;
      
      out->write((char *)&header, sizeof(MasterHeader));
#ifdef DEBUG
      cout << "OutCHKPT: header written" << endl;
#endif
      
    }
    
    for (auto c : comp->components) {
#ifdef DEBUG
      cout << "OutCHKPT: process " << myid << " trying to write name=" << c->name
	   << " force=" << c->id << endl;
#endif
      c->write_binary(out);
#ifdef DEBUG
      cout << "OutCHKPT: process " << myid << " write completed on " << c->name << endl;
#endif
    }
    
    if (myid==0) {
      out->close();
      delete out;
    }

  }

  chktimer.mark();


  if (timer) {
    end = std::chrono::high_resolution_clock::now();
    std::chrono::duration<double> intvl = end - beg;
      if (myid==0)
	std::cout << "OutCHKPT [T=" << tnow << "] timing=" << intvl.count()
		  << std::endl;
  }
}

