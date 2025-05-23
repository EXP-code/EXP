// #define DEBUG_SLEDGE 1
// #define DEBUG_NAN 1

#define NOTIFY 1		// Tell user about *bad* SL solutions

#define SLEDGE_VERBOSE 1	// Turn on verbose output for DEBUG_SLEDGE

#define USE_TABLE 1		// Use all table values rather than analyic
				// evaluations where possible
#define XOFFSET (1.0e-8)

#include <filesystem>
#include <iostream>
#include <sstream>
#include <fstream>
#include <iomanip>
#include <cstdlib>
#include <vector>

#include <EXPException.H>
#include <SLGridMP2.H>
#include <massmodel.H>
#include <EXPmath.H>

#ifdef USE_DMALLOC
#include <dmalloc.h>
#endif

#include <config_exp.h>		// For config macros

// For reading and writing cache file
//
#include <highfive/highfive.hpp>
#include <highfive/eigen.hpp>

// For fortran call
// (This should work both 32-bit and 64-bit . . . )
//
typedef int	logical;
typedef double	doublereal;
typedef int	integer;

MPI_Status status;

//! Target model for slab SL
std::shared_ptr<SlabModel> slab;

extern "C" {
  int sledge_(logical* job, doublereal* cons, logical* endfin, 
	      integer* invec, doublereal* tol, logical* type, 
	      doublereal* ev, integer* numx, doublereal* xef, doublereal* ef, 
	      doublereal* pdef, doublereal* t, doublereal* rho, 
	      integer* iflag, doublereal* store);
}


static
std::string sledge_error(int flag)
{
  if (flag==0)
    return "reliable";
  else if (flag==-1)
    return "too many levels for eigenvalues";
  else if (flag==-2)
    return "too many levels for eigenvectors";
  else if (flag==1)
    return "eigenvalue cluster?";
  else if (flag==2)
    return "classification uncertainty";
  else if (flag>0) {
    std::ostringstream sout;
    sout << "unexpected warning: " << flag; 
    return sout.str();
  } else {
    std::ostringstream sout;
    sout << "unexpected fatal error: " << flag; 
    return sout.str();
  }
}

static double L2, M2, K2;
static int sl_dim;


//======================================================================
//======================================================================
//======================================================================


int SLGridSph::mpi = 0;		// initially off

extern "C" {
  int sledge_(logical* job, doublereal* cons, logical* endfin, 
	      integer* invec, doublereal* tol, logical* type, 
	      doublereal* ev, integer* numx, doublereal* xef, doublereal* ef, 
	      doublereal* pdef, doublereal* t, doublereal* rho, 
	      integer* iflag, doublereal* store);
}


static std::shared_ptr<AxiSymModel> model;

double sphpot(double r)
{
  return model->get_pot(r);
}

double sphdpot(double r)
{
  return model->get_dpot(r);
}

double sphdens(double r)
{
  // This 4pi from Poisson's eqn
  //        |
  //        |       /-- This begins the true density profile
  //        |       |
  //        v       v
  return 4.0*M_PI * model->get_density(r);
}


void SLGridSph::bomb(string oops)
{
  std::ostringstream sout;
  sout << "SLGridSph error [#=" << myid << "]: " << oops;
  throw std::runtime_error(sout.str());
}

				// Constructors

SLGridSph::SLGridSph(std::string modelname,
		     int LMAX, int NMAX, int NUMR,
		     double RMIN, double RMAX, 
		     bool CACHE, int CMAP, double RMAP,
		     int DIVERGE, double DFAC,
		     std::string cachename, bool VERBOSE)
{
  if (modelname.size()) model_file_name = modelname;
  else                  model_file_name = default_model;
  
  if (cachename.size()) sph_cache_name  = cachename;
  else throw std::runtime_error("SLGridSph: you must specify a cachename");
  
  model    = SphModTblPtr(new SphericalModelTable(model_file_name, DIVERGE, DFAC));
  tbdbg    = VERBOSE;
  diverge  = DIVERGE;
  dfac     = DFAC;

  initialize(LMAX, NMAX, NUMR, RMIN, RMAX, CACHE, CMAP, RMAP);
}

SLGridSph::SLGridSph(std::shared_ptr<SphericalModelTable> mod,
		     int LMAX, int NMAX, int NUMR, double RMIN, double RMAX, 
		     bool CACHE, int CMAP, double RMAP,
		     std::string cachename, bool VERBOSE)
{
  model    = mod;
  tbdbg    = VERBOSE;
  diverge  = 0;
  dfac     = 1;

  if (cachename.size()) sph_cache_name  = cachename;
  else throw std::runtime_error("SLGridSph: you must specify a cachename");

  initialize(LMAX, NMAX, NUMR, RMIN, RMAX, CACHE, CMAP, RMAP);
}


SLGridSph::SLGridSph(std::string cachename)
{
  if (cachename.size()) sph_cache_name  = cachename;
  else throw std::runtime_error("SLGridSph: you must specify a cachename");

  tbdbg = false;

  int LMAX, NMAX, NUMR, CMAP, DIVERGE=0;
  double RMIN, RMAX, RMAP, DFAC=1.0;

  try {
    
    auto node = getHeader(cachename);

    LMAX  = node["lmax"].as<int>();
    NMAX  = node["nmax"].as<int>();
    NUMR  = node["numr"].as<int>();
    CMAP  = node["cmap"].as<int>();
    RMIN  = node["rmin"].as<double>();
    RMAX  = node["rmax"].as<double>();
    RMAP  = node["rmapping"].as<double>();

    model_file_name = node["model"].as<std::string>();
    model = SphModTblPtr(new SphericalModelTable(model_file_name, diverge, dfac));
  }
  catch (YAML::Exception& error) {
    std::ostringstream sout;
    sout << "SLGridMP2: error parsing parameters from getHeader: "
	 << error.what();
    throw GenericError(sout.str(), __FILE__, __LINE__, 1039, false);
  }

  initialize(LMAX, NMAX, NUMR, RMIN, RMAX, false, CMAP, RMAP);
}


std::map<std::string, std::string>
SLGridSph::cacheInfo(const std::string& cachefile, bool verbose)
{
  std::map<std::string, std::string> ret;
  auto node = getHeader(cachefile);

  if (verbose)
    std::cout << std::string(60, '-') << std::endl
	      << "Cache parameters for SLGridSph: " << cachefile << std::endl
	      << std::string(60, '-') << std::endl;

  for (YAML::const_iterator it=node.begin(); it!=node.end(); ++it) {
    if (verbose)
      std::cout << std::left << std::setw(20) << it->first.as<std::string>()
		<< ": " << it->second.as<std::string>() << std::endl;
    ret[it->first.as<std::string>()] = it->second.as<std::string>();
  }
  if (verbose)
    std::cout << std::string(60, '-') << std::endl;

  return ret;
}


void SLGridSph::initialize(int LMAX, int NMAX, int NUMR,
			   double RMIN, double RMAX, 
			   bool CACHE, int CMAP, double RMAP)
{
  int l;

  lmax  = LMAX;
  nmax  = NMAX;
  numr  = NUMR;

  rmin  = std::max<double>(RMIN, model->get_min_radius());
  rmax  = std::min<double>(RMAX, model->get_max_radius());

  cache = CACHE;
  cmap  = CMAP;
  rmap  = RMAP;

  init_table();


  if (tbdbg) {
    if (mpi)
      std::cout << "Process " << myid << ": MPI is on!"  << std::endl;
    else
      std::cout << "Process " << myid << ": MPI is off!" << std::endl;
  }

  table = 0;
  
  if (not ReadH5Cache()) {

    // MPI loop
    //
    if (mpi) {

      table = table_ptr_1D(new TableSph [lmax+1]);
      
      mpi_setup();
      
      int totbad = 0;		// Count total number of sledge errors

      if (mpi_myid) {		// Begin workers
	compute_table_worker();
	
	//
	// <Receive completed table from root>
	//

	for (l=0; l<=lmax; l++) {

	  MPI_Bcast(&mpi_buf[0], mpi_bufsz, MPI_PACKED, 0, MPI_COMM_WORLD);
    
	  mpi_unpack_table();      
	}
	
      }
      else {			// BEGIN Root
	
	int worker = 0;
	int request_id = 1;

	l=0;

	while (l<=lmax) {

	  if (worker<mpi_numprocs-1) { // Send request to worker
	    worker++;
	    
	    
	    if (tbdbg)
	      std::cout << "Root sending orders to Worker " << worker 
			<< ": l=" << l << std::endl; 
	    
	    MPI_Send(&request_id, 1, MPI_INT, worker, 11, MPI_COMM_WORLD);
	    MPI_Send(&l, 1, MPI_INT, worker, 11, MPI_COMM_WORLD);
	    
	    if (tbdbg)
	      std::cout << "Root gave orders to Worker " << worker 
			<< ": l=" << l << std::endl;
	    
	    // Increment counters
	    l++;
	  }
	  
	  if (worker == mpi_numprocs-1 && l<=lmax) {
	  
	    //
	    // <Wait and receive>
	    //
	    
#ifdef SLEDGE_THROW
	    int bad;		// Get the sledge error count from a
				// worker
	    MPI_Recv(&bad, 1, MPI_INT, MPI_ANY_SOURCE, 10,
		     MPI_COMM_WORLD, &status);
	    totbad += bad;
	    int retid = status.MPI_SOURCE;
	    MPI_Recv(&mpi_buf[0], mpi_bufsz, MPI_PACKED, retid, 11,
		     MPI_COMM_WORLD, &status);
#else
	    MPI_Recv(&mpi_buf[0], mpi_bufsz, MPI_PACKED, MPI_ANY_SOURCE, 
		     MPI_ANY_TAG, MPI_COMM_WORLD, &status);
	    
	    int retid = status.MPI_SOURCE;
#endif
	  
	    mpi_unpack_table();      

	    //
	    // <Send new request>
	    //
	    
	    MPI_Send(&request_id, 1, MPI_INT, retid, 11, MPI_COMM_WORLD);
	    MPI_Send(&l, 1, MPI_INT, retid, 11, MPI_COMM_WORLD);
	    
	    if (tbdbg) 
	      std::cout << "Root g orders to Worker " << retid
			<< ": l=" << l << std::endl;
	    
	    // Increment counters
	    l++;
	  }
	}
	
	//
	// <Wait for all workers to return>
	//
      
	while (worker) {
#ifdef SLEDGE_THROW
	  int bad;		// Get the sledge error count from a
				// worker
	  MPI_Recv(&bad, 1, MPI_INT, MPI_ANY_SOURCE, 10,
		   MPI_COMM_WORLD, &status);
	  totbad += bad;

	  int retid = status.MPI_SOURCE;
	  MPI_Recv(&mpi_buf[0], mpi_bufsz, MPI_PACKED, retid, 11,
		   MPI_COMM_WORLD, &status);
#else
	  MPI_Recv(&mpi_buf[0], mpi_bufsz, MPI_PACKED, MPI_ANY_SOURCE, 11,
		   MPI_COMM_WORLD, &status);
	  
#endif
	  mpi_unpack_table();      
	  
	  worker--;
	}

	//
	// <Tell workers to continue>
	//
      
	request_id = -1;
	for (worker=1; worker < mpi_numprocs; worker++)
	  MPI_Send(&request_id, 1, MPI_INT, worker, 11, MPI_COMM_WORLD);


	//
	// <Send table to workers>
	//

	for (l=0; l<=lmax; l++) {
	  int position = mpi_pack_table(&table[l], l);

	  MPI_Bcast(&mpi_buf[0], position, MPI_PACKED, 0, MPI_COMM_WORLD);
	}

      }
      // END Root

#ifdef SLEDGE_THROW
      // Share the total sledge error count with all nodes
      MPI_Bcast(&totbad, 1, MPI_INT, 0, MPI_COMM_WORLD);

      // Emit runtime exception on sledge errors
      if (totbad) {
	std::ostringstream sout;

	if (myid==0) {
	  sout << std::endl
	       << "SLGridSph found " << totbad
	       << " tolerance errors in computing SL solutions." << std::endl
	       << "We suggest checking your model file for smoothness and for"
	       << std::endl
	       << "a sufficient number grid points that the relative difference"
	       << std::endl
	       << "between field quantities is <= 0.3";
	} else {
	  sout << std::endl << "sledge tolerance failure";
	}

	throw GenericError(sout.str(), __FILE__, __LINE__);
      }
#endif
    }
    // END MPI stanza, BEGIN single-process stanza
    else {

      table = table_ptr_1D(new TableSph [lmax+1]);

      for (l=0; l<=lmax; l++) {
	if (tbdbg) std::cerr << "Begin [" << l << "] . . ." << std::endl;
	compute_table(&(table[l]), l);
	if (tbdbg) std::cerr << ". . . done" << std::endl;
      }
    }
    // END single process stanza

    // Write cache
    //
    if (myid==0 and cache) WriteH5Cache();
  }
  // END: make tables

  if (tbdbg)
    std::cerr << "Process " << myid << ": exiting constructor" << std::endl;
  
}

void check_vector_values_SL(const Eigen::VectorXd& v)
{
  unsigned c_inf = 0;
  unsigned c_nan = 0;
  unsigned c_sub = 0;

  // Classify all numbers in the vector
  //
  for (int i=0; i<v.size(); i++) {

    switch(std::fpclassify(v[i])) {
    case FP_INFINITE:		// Count infinities
      c_inf++;
      break;
    case FP_NAN:		// Count non-a-numbers
      c_nan++;
      break;
    case FP_SUBNORMAL:		// Count denormalized numbers
      c_sub++;
      break;
    }
  }

  // Print any errors
  //
  if (c_inf+c_nan+c_sub>0) {
    std::cerr << "check_vector [size=" << v.size() << "]: "
	      << " NaN=" << c_nan << " Inf=" << c_inf << "Sub=" << c_sub
	      << std::endl;
  }
}

bool SLGridSph::ReadH5Cache(void)
{
  if (!cache) return false;

  // First attempt to read the file
  //
  try {
    // Silence the HDF5 error stack
    //
    HighFive::SilenceHDF5 quiet;
    
    // Try opening the file as HDF5
    //
    HighFive::File h5file(sph_cache_name, HighFive::File::ReadOnly);
    
    // Try checking the rest of the parameters before reading arrays
    //
    auto checkInt = [&h5file](int value, std::string name)
    {
      int v; HighFive::Attribute vv = h5file.getAttribute(name); vv.read(v);
      if (value == v) return true;
      if (myid==0)
	std::cout << "---- SLGridSph::ReadH5Cache: "
		  << "parameter " << name << ": wanted " << value
		  << " found " << v << std::endl;
      return false;
    };

    auto checkDbl = [&h5file](double value, std::string name)
    {
      double v; HighFive::Attribute vv = h5file.getAttribute(name); vv.read(v);
      if (fabs(value - v) < 1.0e-16) return true;
      if (myid==0)
	std::cout << "---- SLGridSph::ReadH5Cache: "
		  << "parameter " << name << ": wanted " << value
		  << " found " << v << std::endl;
      return false;
    };

    auto checkStr = [&h5file](std::string value, std::string name)
    {
      std::string v; HighFive::Attribute vv = h5file.getAttribute(name); vv.read(v);
      if (value.compare(v)==0) return true;
      if (myid==0)
	std::cout << "---- SLGridSph::ReadH5Cache: "
		  << "parameter " << name << ": wanted " << value
		  << " found " << v << std::endl;
      return false;
    };

    // For cache ID
    //
    std::string geometry("sphere"), forceID("SLGridSph");
    std::string modl(model_file_name);

    // ID check
    //
    if (not checkStr(geometry, "geometry"))  return false;
    if (not checkStr(forceID,  "forceID"))   return false;

    // Version check
    //
    if (h5file.hasAttribute("Version")) {
      if (not checkStr(Version, "Version"))  return false;
    } else {
      if (myid==0)
	std::cout << "---- SLGridSph::ReadH5Cache: "
		  << "recomputing cache for HighFive API change"
		  << std::endl;
      return false;
    }

    // Parameter check
    //
    if (not checkStr(modl,     "model"))     return false;
    if (not checkInt(lmax,     "lmax"))      return false;
    if (not checkInt(nmax,     "nmax"))      return false;
    if (not checkInt(numr,     "numr"))      return false;
    if (not checkInt(cmap,     "cmap"))      return false;
    if (not checkDbl(rmin,     "rmin"))      return false;
    if (not checkDbl(rmax,     "rmax"))      return false;
    if (not checkInt(diverge,  "diverge"))   return false;
    if (not checkDbl(dfac,     "dfac"))      return false;

    // Backward compatibility for old 'scale' key word
    //
    if (h5file.hasAttribute("scale")) {
      if (not checkDbl(rmap,   "scale"))     return false;
    } else {
      if (not checkDbl(rmap,   "rmapping"))  return false;
    }

    // Harmonic order
    //
    auto harmonic = h5file.getGroup("Harmonic");

    // Create table instances
    //
    table = table_ptr_1D(new TableSph [lmax+1]);

    for (int l=0; l<=lmax; l++) {
      std::ostringstream sout;
      sout << l;
      auto arrays = harmonic.getGroup(sout.str());
      
      // Table arrays will be allocated
      //
      arrays.getDataSet("ev").read<Eigen::VectorXd>(table[l].ev);
      arrays.getDataSet("ef").read<Eigen::MatrixXd>(table[l].ef);
    }
    
    if (myid==0)
      std::cout << "---- SLGridSph::ReadH5Cache: "
		<< "successfully read basis cache <" << sph_cache_name
		<< ">" << std::endl;

    return true;
    
  } catch (HighFive::Exception& err) {
    if (myid==0)
      std::cerr << "---- SLGridSph::ReadH5Cache: "
		<< "error reading <" << sph_cache_name << ">" << std::endl
		<< "---- SLGridSph::ReadH5Cache: HDF5 error is <" << err.what()
		<< ">" << std::endl;
  }

  return false;
}



void SLGridSph::WriteH5Cache(void)
{
  if (myid) return;

  try {

    // Check for new HDF5 file
    if (std::filesystem::exists(sph_cache_name)) {
      if (myid==0)
	std::cout << "---- SLGridSph::WriteH5Cache cache file <"
		  << sph_cache_name << "> exists" << std::endl;
      try {
	std::filesystem::rename(sph_cache_name, sph_cache_name + ".bak");
      }
      catch(std::filesystem::filesystem_error const& ex) {
	std::ostringstream sout;
        sout << "---- SLGridSph::WriteH5Cache write error: "
	     << "what():  " << ex.what()  << std::endl
	     << "path1(): " << ex.path1() << std::endl
	     << "path2(): " << ex.path2();
	throw GenericError(sout.str(), __FILE__, __LINE__, 12, true);
      }
      
      if (myid==0)
	std::cout << "---- SLGridSph::WriteH5Cache: existing file backed up to <"
		  << sph_cache_name + ".bak>" << std::endl;
    }
    
    // Create a new hdf5 file
    //
    HighFive::File file(sph_cache_name,
			HighFive::File::ReadWrite | HighFive::File::Create);
    
    // For cache ID
    //
    std::string geometry("sphere"), forceID("SLGridSph");

    file.createAttribute<std::string>("geometry",  HighFive::DataSpace::From(geometry)).write(geometry);
    file.createAttribute<std::string>("forceID",   HighFive::DataSpace::From(forceID)).write(forceID);
    file.createAttribute<std::string>("Version",   HighFive::DataSpace::From(Version)).write(Version);
      
    // Write parameters
    //
    file.createAttribute<std::string> ("model",    HighFive::DataSpace::From(model_file_name)).write(model_file_name);
    file.createAttribute<int>         ("lmax",     HighFive::DataSpace::From(lmax)).write(lmax);
    file.createAttribute<int>         ("nmax",     HighFive::DataSpace::From(nmax)).write(nmax);
    file.createAttribute<int>         ("numr",     HighFive::DataSpace::From(numr)).write(numr);
    file.createAttribute<int>         ("cmap",     HighFive::DataSpace::From(cmap)).write(cmap);
    file.createAttribute<double>      ("rmin",     HighFive::DataSpace::From(rmin)).write(rmin);
    file.createAttribute<double>      ("rmax",     HighFive::DataSpace::From(rmax)).write(rmax);
    file.createAttribute<double>      ("rmapping", HighFive::DataSpace::From(rmap)).write(rmap);
    file.createAttribute<int>         ("diverge",  HighFive::DataSpace::From(diverge)).write(diverge);
    file.createAttribute<double>      ("dfac",     HighFive::DataSpace::From(dfac)).write(dfac);
      
    // Harmonic order (for h5dump readability)
    //
    auto harmonic = file.createGroup("Harmonic");

    for (int l=0; l<=lmax; l++) {
      std::ostringstream sout;
      sout << l;
      auto arrays = harmonic.createGroup(sout.str());
      
      arrays.createDataSet("ev",   table[l].ev);
      arrays.createDataSet("ef",   table[l].ef);
    }
    
  } catch (HighFive::Exception& err) {
    std::cerr << err.what() << std::endl;
  }
    
  std::cout << "---- SLGridSph::WriteH5Cache: "
	    << "wrote <" << sph_cache_name << ">" << std::endl;
  
  return ;
}


SLGridSph::~SLGridSph()
{
  // Nothing
}

				// Members

double SLGridSph::eigenvalue(int l, int n)
{
  return table.get()[l].ev[n];
}

double SLGridSph::r_to_xi(double r)
{
  double ret;

  if (cmap==1) {
    if (r<0.0) bomb("radius < 0!");
    ret =  (r/rmap-1.0)/(r/rmap+1.0);
  } else if (cmap==2) {
    if (r<=0.0) bomb("radius <= 0!");
    ret = log(r);
  } else {
    ret = r;
  }    

  return ret;
}
    
double SLGridSph::xi_to_r(double xi)
{
  double ret;

  if (cmap==1) {
    if (xi<-1.0) bomb("xi < -1!");
    if (xi>=1.0) bomb("xi >= 1!");

    ret =(1.0+xi)/(1.0 - xi) * rmap;
  } else if (cmap==2) {
    ret = exp(xi);
  } else {
    if (xi<0.0) bomb("xi < 0!");

    ret = xi;
  }

  return ret;
}

double SLGridSph::d_xi_to_r(double xi)
{
  double ret;

  if (cmap==1) {
    if (xi<-1.0) bomb("xi < -1!");
    if (xi>=1.0) bomb("xi >= 1!");

    ret = 0.5*(1.0-xi)*(1.0-xi)/rmap;
  } else if (cmap==2) {
    ret = exp(-xi);
  } else {
    if (xi<0.0) bomb("xi < 0!");
    ret = 1.0;
  }

  return ret;
}

double SLGridSph::get_pot(double x, int l, int n, int which)
{
  if (which || !cmap)
    x = r_to_xi(x);
  else {
    if (cmap==1) {
      if (x<-1.0) x=-1.0;
      if (x>=1.0) x=1.0-XOFFSET;
    }
    if (cmap==2) {
      if (x<xmin) x=xmin;
      if (x>xmax) x=xmax;
    }
  }


  int indx = (int)( (x-xmin)/dxi );
  if (indx<0) indx = 0;
  if (indx>numr-2) indx = numr - 2;


  double x1 = (xi[indx+1] - x)/dxi;
  double x2 = (x - xi[indx])/dxi;
  

#ifdef USE_TABLE
  return (x1*table[l].ef(n, indx) + x2*table[l].ef(n, indx+1))/
    sqrt(table[l].ev[n]) * (x1*p0[indx] + x2*p0[indx+1]);
#else
  return (x1*table[l].ef(n, indx) + x2*table[l].ef(n, indx+1))/
    sqrt(table[l].ev[n]) * sphpot(xi_to_r(x));
#endif
}


double SLGridSph::get_dens(double x, int l, int n, int which)
{
  if (which || !cmap)
    x = r_to_xi(x);
  else {
    if (cmap==1) {
      if (x<-1.0) x=-1.0;
      if (x>=1.0) x=1.0-XOFFSET;
    }
    if (cmap==2) {
      if (x<xmin) x=xmin;
      if (x>xmax) x=xmax;
    }
  }

  int indx = (int)( (x-xmin)/dxi );
  if (indx<0) indx = 0;
  if (indx>numr-2) indx = numr - 2;


  double x1 = (xi[indx+1] - x)/dxi;
  double x2 = (x - xi[indx])/dxi;
  
#ifdef USE_TABLE
  return (x1*table[l].ef(n, indx) + x2*table[l].ef(n, indx+1)) *
    sqrt(table[l].ev[n]) * (x1*d0[indx] + x2*d0[indx+1]);
#else
  return (x1*table[l].ef(n, indx) + x2*table[l].ef(n, indx+1)) *
    sqrt(table[l].ev[n]) * sphdens(xi_to_r(x));
#endif

}

double SLGridSph::get_force(double x, int l, int n, int which)
{
  if (which || !cmap)
    x = r_to_xi(x);
  else {
    if (cmap==1) {
      if (x<-1.0) x=-1.0;
      if (x>=1.0) x=1.0-XOFFSET;
    }
    if (cmap==2) {
      if (x<xmin) x=xmin;
      if (x>xmax) x=xmax;
    }
  }


  int indx = (int)( (x-xmin)/dxi );
  if (indx<1) indx = 1;
  if (indx>numr-2) indx = numr - 2;


  double p = (x - xi[indx])/dxi;
  
				// Use three point formula

				// Point -1: indx-1
				// Point  0: indx
				// Point  1: indx+1

  return d_xi_to_r(x)/dxi * (
			     (p - 0.5)*table[l].ef(n, indx-1)*p0[indx-1]
			     -2.0*p*table[l].ef(n, indx)*p0[indx]
			     + (p + 0.5)*table[l].ef(n, indx+1)*p0[indx+1]
			     ) / sqrt(table[l].ev[n]);
}


void SLGridSph::get_pot(Eigen::MatrixXd& mat, double x, int which)
{
  if (which || !cmap)
    x = r_to_xi(x);
  else {
    if (cmap==1) {
      if (x<-1.0) x=-1.0;
      if (x>=1.0) x=1.0-XOFFSET;
    }
    if (cmap==2) {
      if (x<xmin) x=xmin;
      if (x>xmax) x=xmax;
    }
  }

  mat.resize(lmax+1, nmax);

  int indx = (int)( (x-xmin)/dxi );
  if (indx<0) indx = 0;
  if (indx>numr-2) indx = numr - 2;


  double x1 = (xi[indx+1] - x)/dxi;
  double x2 = (x - xi[indx])/dxi;
  

  for (int l=0; l<=lmax; l++) {
    for (int n=0; n<nmax; n++) {
#ifdef USE_TABLE
      mat(l, n) = (x1*table[l].ef(n, indx) + x2*table[l].ef(n, indx+1))/
	sqrt(table[l].ev[n]) * (x1*p0[indx] + x2*p0[indx+1]);
#else
      mat(l, n) = (x1*table[l].ef(n, indx) + x2*table[l].ef(n, indx+1))/
	sqrt(table[l].ev[n]) * sphpot(xi_to_r(x));
#endif
    }
  }

}


void SLGridSph::get_dens(Eigen::MatrixXd& mat, double x, int which)
{
  if (which || !cmap)
    x = r_to_xi(x);
  else {
    if (cmap==1) {
      if (x<-1.0) x=-1.0;
      if (x>=1.0) x=1.0-XOFFSET;
    }
    if (cmap==2) {
      if (x<xmin) x=xmin;
      if (x>xmax) x=xmax;
    }
  }

  mat.resize(lmax+1, nmax);

  int indx = (int)( (x-xmin)/dxi );
  if (indx<0) indx = 0;
  if (indx>numr-2) indx = numr - 2;


  double x1 = (xi[indx+1] - x)/dxi;
  double x2 = (x - xi[indx])/dxi;
  

  for (int l=0; l<=lmax; l++) {
    for (int n=0; n<nmax; n++) {
#ifdef USE_TABLE
      mat(l, n) = (x1*table[l].ef(n, indx) + x2*table[l].ef(n, indx+1))*
	sqrt(table[l].ev[n]) * (x1*d0[indx] + x2*d0[indx+1]);
#else
      mat(l, n) = (x1*table[l].ef(n, indx) + x2*table[l].ef(n, indx+1))*
	sqrt(table[l].ev[n]) * sphdens(xi_to_r(x));
#endif
    }
  }

}


void SLGridSph::get_force(Eigen::MatrixXd& mat, double x, int which)
{
  if (which || !cmap)
    x = r_to_xi(x);
  else {
    if (cmap==1) {
      if (x<-1.0) x=-1.0;
      if (x>=1.0) x=1.0-XOFFSET;
    }
    if (cmap==2) {
      if (x<xmin) x=xmin;
      if (x>xmax) x=xmax;
    }
  }

  mat.resize(lmax+1, nmax);

  int indx = (int)( (x-xmin)/dxi );
  if (indx<1) indx = 1;
  if (indx>numr-2) indx = numr - 2;


  double p = (x - xi[indx])/dxi;
  double fac = d_xi_to_r(x)/dxi;

  for (int l=0; l<=lmax; l++) {
    for (int n=0; n<nmax; n++) {
      mat(l, n) = fac * (
			 (p - 0.5)*table[l].ef(n, indx-1)*p0[indx-1]
			 -2.0*p*table[l].ef(n, indx)*p0[indx]
			 + (p + 0.5)*table[l].ef(n, indx+1)*p0[indx+1]
			 ) / sqrt(table[l].ev[n]);
    }
  }
  
}


void SLGridSph::get_pot(Eigen::VectorXd& vec, double x, int l, int which)
{
  if (which || !cmap)
    x = r_to_xi(x);
  else {
    if (cmap==1) {
      if (x<-1.0) x=-1.0;
      if (x>=1.0) x=1.0-XOFFSET;
    }
    if (cmap==2) {
      if (x<xmin) x=xmin;
      if (x>xmax) x=xmax;
    }
  }

  vec.resize(nmax);

  int indx = (int)( (x-xmin)/dxi );
  if (indx<0) indx = 0;
  if (indx>numr-2) indx = numr - 2;

  double x1 = (xi[indx+1] - x)/dxi;
  double x2 = (x - xi[indx])/dxi;
  

  for (int n=0; n<nmax; n++) {
#ifdef USE_TABLE
    vec[n] = (x1*table[l].ef(n, indx) + x2*table[l].ef(n, indx+1))/
      sqrt(table[l].ev[n]) * (x1*p0[indx] + x2*p0[indx+1]);
#else
    vec[n] = (x1*table[l].ef(n, indx) + x2*table[l].ef(n, indx+1))/
      sqrt(table[l].ev[n]) * sphpot(xi_to_r(x));
#endif
  }

}


void SLGridSph::get_dens(Eigen::VectorXd& vec, double x, int l, int which)
{
  if (which || !cmap)
    x = r_to_xi(x);
  else {
    if (cmap==1) {
      if (x<-1.0) x=-1.0;
      if (x>=1.0) x=1.0-XOFFSET;
    }
    if (cmap==2) {
      if (x<xmin) x=xmin;
      if (x>xmax) x=xmax;
    }
  }

  vec.resize(nmax);

  int indx = (int)( (x-xmin)/dxi );
  if (indx<0) indx = 0;
  if (indx>numr-2) indx = numr - 2;

  double x1 = (xi[indx+1] - x)/dxi;
  double x2 = (x - xi[indx])/dxi;
  

  for (int n=0; n<nmax; n++) {
#ifdef USE_TABLE
    vec[n] = (x1*table[l].ef(n, indx) + x2*table[l].ef(n, indx+1))*
      sqrt(table[l].ev[n]) * (x1*d0[indx] + x2*d0[indx+1]);
#else
    vec[n] = (x1*table[l].ef(n, indx) + x2*table[l].ef(n, indx+1))*
      sqrt(table[l].ev[n]) * sphdens(xi_to_r(x));
#endif
  }

}


void SLGridSph::get_force(Eigen::VectorXd& vec, double x, int l, int which)
{
  if (which || !cmap)
    x = r_to_xi(x);
  else {
    if (cmap==1) {
      if (x<-1.0) x=-1.0;
      if (x>=1.0) x=1.0-XOFFSET;
    }
    if (cmap==2) {
      if (x<xmin) x=xmin;
      if (x>xmax) x=xmax;
    }
  }

  vec.resize(nmax);

  int indx = (int)( (x-xmin)/dxi );
  if (indx<1) indx = 1;
  if (indx>numr-2) indx = numr - 2;


  double p = (x - xi[indx])/dxi;
  double fac = d_xi_to_r(x)/dxi;

  for (int n=0; n<nmax; n++) {
    vec[n] = fac * (
		    (p - 0.5)*table[l].ef(n, indx-1)*p0[indx-1]
		    -2.0*p*table[l].ef(n, indx)*p0[indx]
		    + (p + 0.5)*table[l].ef(n, indx+1)*p0[indx+1]
		    ) / sqrt(table[l].ev[n]);
  }

}

void SLGridSph::compute_table(struct TableSph* table, int l)
{

  double cons[8] = {0.0, 0.0, 0.0, 0.0,   0.0, 0.0,   0.0, 0.0};
  double tol [6] = {1.0e-4*rmap, 1.0e-6,  
                    1.0e-4*rmap, 1.0e-6,  
		    1.0e-4*rmap, 1.0e-6};
  int VERBOSE=0;
  integer NUM, N;
  logical type[8] = {0, 0, 1, 0, 0, 0, 1, 0};
  logical endfin[2] = {1, 1};
  
#ifdef DEBUG_SLEDGE
  if (myid==0) VERBOSE = SLEDGE_VERBOSE;
#endif

  cons[6] = rmin;
  cons[7] = rmax;
  L2 = l*(l+1);
  NUM = numr;
  N = nmax;
  
  // integer iflag[nmax], invec[nmax+3];
  integer *iflag = new integer [nmax];
  integer *invec = new integer [nmax+3];

  double *t=0, *rho=0;
  double *ev    = new double [N];
  double *store = new double [26*(NUM+16)];
  double *xef   = new double [NUM+16];
  double *ef    = new double [NUM*N];
  double *pdef  = new double [NUM*N];
  double f;

				// Inner BC
  f = sphpot(cons[6]);
  if (l==0) {
    cons[0] = sphdpot(cons[6])/f;
    cons[2] = 1.0/(cons[6]*cons[6]*f*f);
  }
  else
    cons[0] = 1.0;

				// Outer BC
  f = sphpot(cons[7]);
  cons[4] = (1.0 + l)/cons[7] + sphdpot(cons[7])/f;
  cons[5] = 1.0/(cons[7]*cons[7]*f*f);

  //
  //     Initialize the vector INVEC(*):
  //       estimates for the eigenvalues/functions specified
  //

  invec[0] = VERBOSE;		// little printing (1), no printing (0)
  invec[1] = 3;			// spectrum is ignored
  invec[2] = N;			// estimates for N eigenvalues/functions

  for (int i=0; i<N; i++) invec[3+i] = i;

  //
  //     Set the JOB(*) vector:
  //        estimate both eigenvalues and eigenvectors,
  //        don't estimate the spectral density function,
  //        classify,
  //        let SLEDGE choose the initial mesh
  //
  logical job[5] = {0,1,0,1,0};

  //
  //     Output mesh
  //
  for (int i=0; i<NUM; i++) xef[i] = r[i];

  //     
  //     Open file for output.
  //
  sl_dim = 3;

  sledge_(job, cons, endfin, invec, tol, type, ev, &NUM, xef, ef, pdef,
	   t, rho, iflag, store);

  //
  //     Check for errors
  //
#ifdef SLEDGE_THROW
  // Accumulate error count
  unsigned bad = 0;
  for (int i=0; i<N; i++) {
    if (iflag[i] != 0) bad++;
  }

  std::ostringstream sout;

  // Print info if we have errors and throw a runtime error
  if (bad>0) {

    if (myid==0) {

      std::cout.precision(6);
      std::cout.setf(ios::scientific);
      std::cout << std::left;
    
      std::cout << std::endl
		<< "Tolerance errors in Sturm-Liouville solver for l=" << l
		<<  std::endl << std::endl;

      std::cout << std::setw(15) << "order"
		<< std::setw(15) << "eigenvalue"
		<< std::setw(40) << "condition"
		<< std::endl
		<< std::setw(15) << "-----"
		<< std::setw(15) << "----------"
		<< std::setw(40) << "---------"
		<< std::endl;
      
      for (int i=0; i<N; i++) {
	std::cout << std::setw(15) << invec[3+i] 
		  << std::setw(15) << ev[i]
		  << std::setw(40) << sledge_error(iflag[i])
		  << std::endl;
      }
      std::cout << std::endl;
      
      sout << std::endl
	   << "SLGridSph found " << bad
	   << " tolerance errors in computing SL solutions." << std::endl
	   << "We suggest checking your model file for smoothness and ensure"
	   << std::endl
	   << "a sufficient number grid points that the relative difference"
	   << std::endl
	   << "between field quantities is <= 0.3";

      throw GenericError(sout.str(), __FILE__, __LINE__);
    }
    else {
      throw GenericError("sledge errors", __FILE__, __LINE__);
    }
  }
#endif

  //
  //     Print results:
  //
  if (tbdbg) {
    std::cout.precision(6);
    std::cout.setf(ios::scientific);

    for (int i=0; i<N; i++) {
      std::cout << std::setw(15) << invec[3+i] 
		<< std::setw(15) << ev[i]
		<< std::setw( 5) << iflag[i]
		<< std::endl;
      
      if (VERBOSE) {
	
	if (iflag[i] > -10) {
	  std::cout << std::setw(14) << "x"
		    << std::setw(25) << "u(x)"
		    << std::setw(25) << "(pu`)(x)"
		    << std::endl;
	  int k = NUM*i;
	  for (int j=0; j<NUM; j++) {
	    std::cout << std::setw(25) << xef[j]
		      << std::setw(25) << ef[j+k]
		      << std::setw(25) << pdef[j+k]
		      << std::endl;
	  }
	}
	
      }

    }
  }
  
				// Load table

  table->ev.resize(N);
  for (int i=0; i<N; i++) table->ev[i] = ev[i];

  table->ef.resize(N, numr);

  // Choose sign conventions for the ef table
  //
  int nfid = std::min<int>(nevsign, NUM) - 1;
  Eigen::VectorXi sgn = Eigen::VectorXi::Ones(N);
  for (int j=0; j<N; j++) {
    if (ef[j*NUM+nfid]<0.0) sgn(j) = -1;
  }
  
  for (int i=0; i<numr; i++) {
    for(int j=0; j<N; j++) 
      table->ef(j, i) = ef[j*NUM+i] * sgn(j);
  }

  table->l = l;

  delete [] iflag;
  delete [] invec;
  delete [] ev;
  delete [] store;
  delete [] xef;
  delete [] ef;
  delete [] pdef;

}


void SLGridSph::init_table(void)
{
  xi.resize(numr);
  r. resize(numr);
  p0.resize(numr);
  d0.resize(numr);

  if (cmap==1) {
    xmin = (rmin/rmap - 1.0)/(rmin/rmap + 1.0);
    xmax = (rmax/rmap - 1.0)/(rmax/rmap + 1.0);
  }
  else if (cmap==2) {
    xmin = log(rmin);
    xmax = log(rmax);
  } else {
    xmin = rmin;
    xmax = rmax;
  }
  dxi = (xmax-xmin)/(numr-1);
    
  for (int i=0; i<numr; i++) {
    xi[i] = xmin + dxi*i;
    r[i]  = xi_to_r(xi[i]);
    p0[i] = sphpot(r[i]);
    d0[i] = sphdens(r[i]);
  }

}


void SLGridSph::compute_table_worker(void)
{

  double cons[8] = {0.0, 0.0, 0.0, 0.0,   0.0, 0.0,   0.0, 0.0};
  double tol [6] = {1.0e-1*rmap, 1.0e-6,  
                    1.0e-1*rmap, 1.0e-6,  
		    1.0e-1*rmap, 1.0e-6};

  int VERBOSE=0;
  integer NUM;
  logical type[8] = {0, 0, 1, 0, 0, 0, 1, 0};
  logical endfin[2] = {1, 1};
  
  struct TableSph table;
  int L, N;

#ifdef DEBUG_SLEDGE
  if (myid==0) VERBOSE = SLEDGE_VERBOSE;
#endif

#ifdef DEBUG
  std::cout << "Worker " << mpi_myid << " begins . . ." << std::endl;
#endif

  //
  // <Wait for orders>
  //
  
  int request_id;

  while(1) {

    MPI_Recv(&request_id, 1, 
	     MPI_INT, 0, MPI_ANY_TAG, MPI_COMM_WORLD, &status);

    if (request_id < 0) break;	// Good-bye

    MPI_Recv(&L, 1, 
	     MPI_INT, 0, MPI_ANY_TAG, MPI_COMM_WORLD, &status);


    if (tbdbg)
      std::cout << "Worker " <<  mpi_myid << ": ordered to compute l = " << L << "" << std::endl;

    cons[0] = cons[1] = cons[2] = cons[3] = cons[4] = cons[5] = 0.0;
    cons[6] = rmin;
    cons[7] = rmax;
    L2 = L*(L+1);
    NUM = numr;
    N = nmax;

    // integer iflag[nmax], invec[nmax+3];
    integer *iflag = new integer [nmax];
    integer *invec = new integer [nmax+3];

    double *t=0, *rho=0;
    double *ev    = new double [N];
    double *store = new double [26*(NUM+16)];
    double *xef   = new double [NUM+16];
    double *ef    = new double [NUM*N];
    double *pdef  = new double [NUM*N];
    double f;

    f = sphpot(cons[6]);
    cons[2] = -1.0/(cons[6]*cons[6]*f*f);
    cons[4] = L/cons[7];
    f = sphpot(cons[7]);
    cons[5] = 1.0/(cons[7]*cons[7]*f*f);

    //
    //     Initialize the vector INVEC(*):
    //       estimates for the eigenvalues/functions specified
    //

    invec[0] = VERBOSE;		// little printing (1), no printing (0)
    invec[1] = 3;		// spectrum is ignored
    invec[2] = N;		// estimates for N eigenvalues/functions

    for (int i=0; i<N; i++) invec[3+i] = i;

    //
    //     Set the JOB(*) vector:
    //        estimate both eigenvalues and eigenvectors,
    //        don't estimate the spectral density function,
    //        classify,
    //        let SLEDGE choose the initial mesh
    //
    logical job[5] = {0,1,0,1,0};

    //
    //     Output mesh
    //
    for (int i=0; i<NUM; i++) xef[i] = r[i];

    //     
    //     Open file for output.
    //
    sl_dim = 3;

    sledge_(job, cons, endfin, invec, tol, type, ev, &NUM, xef, ef, pdef,
	    t, rho, iflag, store);

    //
    //     Check for errors
    //
#ifdef SLEDGE_THROW
    // Get sledge error count
    unsigned bad = 0;
    for (int i=0; i<N; i++) {
      if (iflag[i] != 0) bad++;
    }
    
    std::ostringstream sout;

    // Print info if we have errors.  Number of errors will be sent to
    // and accumulated by the root node.
    if (bad>0) {

      std::cout.precision(6);
      std::cout.setf(ios::scientific);
      std::cout << std::left;
    
      std::cout << std::endl
		<< "Tolerance errors in Sturm-Liouville solver for l=" << L
		<<  std::endl << std::endl;

      std::cout << std::setw(15) << "order"
		<< std::setw(15) << "eigenvalue"
		<< std::setw(40) << "condition"
		<< std::endl
		<< std::setw(15) << "-----"
		<< std::setw(15) << "----------"
		<< std::setw(40) << "---------"
		<< std::endl;
      
      for (int i=0; i<N; i++) {
	std::cout << std::setw(15) << invec[3+i] 
		  << std::setw(15) << ev[i]
		  << std::setw(40) << sledge_error(iflag[i])
		  << std::endl;
      }
      std::cout << std::endl;
      
    }
#endif

    //
    //     Print results:
    //
    if (tbdbg) {
      std::cout << "Worker " <<  mpi_myid << ": computed l = " << L << "" << std::endl;

      std::cout.precision(6);
      std::cout.setf(ios::scientific);

      for (int i=0; i<N; i++) {
	std::cout << std::setw(15) << invec[3+i] 
		  << std::setw(15) << ev[i]
		  << std::setw( 5) << iflag[i]
		  << std::endl;
	
	if (VERBOSE) {

	  if (iflag[i] > -10) {
	    std::cout << std::setw(14) << "x"
		      << std::setw(25) << "u(x)"
		      << std::setw(25) << "(pu`)(x)"
		      << std::endl;
	    int k = NUM*i;
	    for (int j=0; j<NUM; j++) {
	      std::cout << std::setw(25) << xef[j]
			<< std::setw(25) << ef[j+k]
			<< std::setw(25) << pdef[j+k]
			<< std::endl;
	    }
	  }
	}
      }
    }

    // Load table
    //
    table.ev.resize(N);
    for (int i=0; i<N; i++) table.ev[i] = ev[i];

    // Choose sign conventions for the ef table
    //
    int nfid = std::min<int>(nevsign, NUM) - 1;
    Eigen::VectorXi sgn = Eigen::VectorXi::Ones(N);
    for (int j=0; j<N; j++) {
      if (ef[j*NUM+nfid]<0.0) sgn(j) = -1;
    }

    table.ef.resize(N, numr);
    for (int i=0; i<numr; i++) {
      for (int j=0; j<N; j++) 
	table.ef(j, i) = ef[j*NUM+i] * sgn(j);
    }

    table.l = L;
  
#ifdef SLEDGE_THROW
    // Send sledge error count to root
    MPI_Send(&bad, 1, MPI_INT, 0, 10, MPI_COMM_WORLD);
#endif
    // Send sledge comptuation to root
    int position = mpi_pack_table(&table, L);
    MPI_Send(&mpi_buf[0], position, MPI_PACKED, 0, 11, MPI_COMM_WORLD);

    if (tbdbg)
      std::cout << "Worker " <<  mpi_myid << ": send to root l = " << L << "" << std::endl;

    delete [] iflag;
    delete [] invec;
    delete [] ev;
    delete [] store;
    delete [] xef;
    delete [] ef;
    delete [] pdef;

  }

}


void SLGridSph::mpi_setup(void)
{
				// Get MPI id

  MPI_Comm_rank(MPI_COMM_WORLD, &mpi_myid);
  MPI_Comm_size(MPI_COMM_WORLD, &mpi_numprocs);


  // Setup for pack/unpack

  int buf1, buf2;
  MPI_Pack_size( 1, MPI_INT, MPI_COMM_WORLD, &buf1);
  MPI_Pack_size( 1, MPI_DOUBLE, MPI_COMM_WORLD, &buf2);

  mpi_bufsz = buf1 +		// l
    nmax*buf2 +			// ev
    nmax*numr*buf2 ;		// ef

  mpi_buf = std::shared_ptr<char[]>(new char [mpi_bufsz]);
}


int SLGridSph::mpi_pack_table(struct TableSph* table, int l)
{
  int position = 0;

  MPI_Pack( &l, 1, MPI_INT, &mpi_buf[0], mpi_bufsz, 
	    &position, MPI_COMM_WORLD);

  for (int j=0; j<nmax; j++)
    MPI_Pack( &table->ev[j], 1, MPI_DOUBLE, &mpi_buf[0], mpi_bufsz, 
	      &position, MPI_COMM_WORLD);

  for (int j=0; j<nmax; j++)
    for (int i=0; i<numr; i++)
      MPI_Pack( &table->ef(j, i), 1, MPI_DOUBLE, &mpi_buf[0], mpi_bufsz, 
		&position, MPI_COMM_WORLD);

  return position;
}


void SLGridSph::mpi_unpack_table(void)
{
  int l, length, position = 0;

  /*
  MPI_Get_count( &status, MPI_PACKED, &length);
  */
  length = mpi_bufsz;

  
  int retid = status.MPI_SOURCE;

  MPI_Unpack( &mpi_buf[0], length, &position, &l, 1, MPI_INT,
	      MPI_COMM_WORLD);

  if (tbdbg)
    std::cout << "Process " <<  mpi_myid << ": unpacking table entry from Process " 
	      << retid << "  l = " << l << "" << std::endl;


  table[l].l = l;
  table[l].ev.resize(nmax);
  table[l].ef.resize(nmax, numr);

  for (int j=0; j<nmax; j++)
    MPI_Unpack( &mpi_buf[0], length, &position, &table[l].ev[j], 1, MPI_DOUBLE,
		MPI_COMM_WORLD);

  for (int j=0; j<nmax; j++)
    for (int i=0; i<numr; i++)
      MPI_Unpack( &mpi_buf[0], length, &position, &table[l].ef(j, i), 1, 
		  MPI_DOUBLE, MPI_COMM_WORLD);
}


YAML::Node SLGridSph::getHeader(const std::string& cachefile)
{
  std::ifstream in(cachefile);
  if (!in) {
    std::ostringstream sout;
    sout << "SLGridSph::getHeader: could not open cache file <" << cachefile << ">";
    std::runtime_error(sout.str());
  }

  YAML::Node node;

  try {
    // Silence the HDF5 error stack
    //
    HighFive::SilenceHDF5 quiet;
    
    // Open the hdf5 file
    //
    HighFive::File file(cachefile, HighFive::File::ReadOnly);
    
    auto getInt = [&file](std::string name)
    {
      int v;
      HighFive::Attribute vv = file.getAttribute(name);
      vv.read(v);
      return v;
    };

    auto getDbl = [&file](std::string name)
    {
      double v;
      HighFive::Attribute vv = file.getAttribute(name);
      vv.read(v);
      return v;
    };

    auto getStr = [&file](std::string name)
    {
      std::string v;
      HighFive::Attribute vv = file.getAttribute(name);
      vv.read(v);
      return v;
    };

    node["model"]    = getStr("model");
    node["lmax"]     = getInt("lmax");
    node["nmax"]     = getInt("nmax");
    node["numr"]     = getInt("numr");
    node["cmap"]     = getInt("cmap");
    node["rmin"]     = getDbl("rmin");
    node["rmax"]     = getDbl("rmax");
    node["rmapping"] = getDbl("rmapping");
    node["diverge"]  = getInt("diverge");
    node["dfac"]     = getDbl("dfac");
  }
  catch (YAML::Exception& error) {
    std::ostringstream sout;
    sout << "SLGridMP2::getHeader: invalid cache file <" << cachefile << ">. ";
    sout << "YAML error in getHeader: " << error.what();
    throw GenericError(sout.str(), __FILE__, __LINE__, 1038, false);
  }

  return node;
}


std::vector<Eigen::MatrixXd> SLGridSph::orthoCheck(int num)
{
  // Gauss-Legendre knots and weights
  LegeQuad lw(num);

  // Get the scaled coordinate limits
  double ximin = r_to_xi(rmin);
  double ximax = r_to_xi(rmax);

  // Initialize the return matrices
  std::vector<Eigen::MatrixXd> ret(lmax+1);
  for (auto & v : ret) v.resize(nmax, nmax);

  // Do each harmonic order
  for (int L=0; L<=lmax; L++) {

    // Unroll the loop for OpenMP parallelization
#pragma omp parallel for
    for (int nn=0; nn<nmax*nmax; nn++) {
      int n1 = nn/nmax;
      int n2 = nn - n1*nmax;
      
      // The inner product
      double ans=0.0;
      for (int i=0; i<num; i++) {
	
	double x = ximin + (ximax - ximin)*lw.knot(i);
	double r = xi_to_r(x);
	  
	ans += r*r*get_pot(x, L, n1, 0)*
	  get_dens(x, L, n2, 0) /
	  d_xi_to_r(x) * (ximax - ximin)*lw.weight(i);
	  
      }
      // END: inner product
	    
      // Assign the matrix element
      //
      ret[L](n1, n2) = - ans;
      //               ^
      //               |
      //               +--- Switch to normed scalar product rather
      //                    that normed gravitational energy
    }
    // END: unrolled loop
  }
  // END: harmonic order loop
  
  return ret;
}

//======================================================================
//======================================================================
//======================================================================


int    SLGridSlab::mpi   = 0;	// initially off
int    SLGridSlab::cache = 1;	// initially yes
double SLGridSlab::H     = 0.1;	// Scale height
double SLGridSlab::L     = 1.0;	// Periodic box size
double SLGridSlab::ZBEG  = 0.0;	// Offset on from origin
double SLGridSlab::ZEND  = 0.0;	// Offset on potential zero

static double KKZ;

static double poffset=0.0;

// Isothermal slab with G = M = 1
//
class IsothermalSlab : public SlabModel
{

public:

  IsothermalSlab() { id = "iso"; }

  double pot(double z)
  {
    return 2.0*M_PI*SLGridSlab::H*log(cosh(z/SLGridSlab::H)) - poffset;
  }

  double dpot(double z)
  {
    return 2.0*M_PI*tanh(z/SLGridSlab::H);
  }

  double dens(double z)
  {
    double tmp = 1.0/cosh(z/SLGridSlab::H);
    return 4.0*M_PI * 0.5/SLGridSlab::H * tmp*tmp;
  }
};


//! Constant density slab with G = M = 1
class ConstantSlab : public SlabModel
{

public:

  ConstantSlab()  { id = "const"; }

  double pot(double z)
  {
    return z*z/(4.0*SLGridSlab::H) - poffset;
  }

  double dpot(double z)
  {
    return z/(2.0*SLGridSlab::H);
  }

  double dens(double z)
  {
    return 4.0*M_PI / (2.0 * SLGridSlab::H);
  }
};

//! Parabolic density slab with G = M = 1
class ParabolicSlab : public SlabModel
{

public:

  ParabolicSlab() { id = "para"; }

  double pot(double z)
  {
    double z2 = z*z;
    double h  = SLGridSlab::H;
    double h2 = h*h;
    return z2*(6.0*h2 - z2)/(16.0*h*h2) - poffset;
  }

  double dpot(double z)
  {
    double z2 = z*z;
    double h  = SLGridSlab::H;
    double h2 = h*h;
    return z*(3.0*h2 - z2)/(4.0*h*h2);
  }

  double dens(double z)
  {
    double h  = SLGridSlab::H;
    double h2 = h*h;
    return 4.0*M_PI * 3.0*(1.0 - z*z/h2)/(4.0*h);
  }
};


std::shared_ptr<SlabModel> SlabModel::createModel(const std::string type)
{
  std::string data(type);
  std::transform(data.begin(), data.end(), data.begin(),
		 [](unsigned char c){ return std::tolower(c); });

  if (data.find("iso") != std::string::npos) {
    return std::make_shared<IsothermalSlab>();
  }

  if (data.find("para") != std::string::npos) {
    return std::make_shared<ParabolicSlab>();
  }

  if (data.find("const") != std::string::npos) {
    return std::make_shared<ConstantSlab>();
  }

  // Default
  return std::make_shared<IsothermalSlab>();
}


void SLGridSlab::bomb(string oops)
{
  std::ostringstream sout;
  sout << "SLGridSlab error [#=" << myid << "]: " << oops;
  throw std::runtime_error(sout.str());
}

				// Constructors

SLGridSlab::SLGridSlab(int NUMK, int NMAX, int NUMZ, double ZMAX,
		       const std::string TYPE, bool VERBOSE)
{
  int kx, ky;

  numk = NUMK;
  nmax = NMAX;
  numz = NUMZ;
  type = TYPE;

  zmax = ZMAX;

  slab  = SlabModel::createModel(type);

  poffset = slab->pot((1.0+ZEND)*zmax);

  tbdbg   = VERBOSE;

  // This could be controlled by a parameter...but at this point, this
  // is a fixed tuning.
  mM      = CoordMap::factory(CoordMapTypes::Sech, H);

  init_table();

  if (tbdbg) {
    if (mpi)
      std::cout << "Process " << myid << ": MPI is on!"  << std::endl;
    else
      std::cout << "Process " << myid << ": MPI is off!" << std::endl;
  }

  table = table_ptr_2D(new table_ptr_1D [numk+1]);
  for (kx=0; kx<=numk; kx++)
    table[kx] = table_ptr_1D(new TableSlab [kx+1]);

  if (mpi) {

    mpi_setup();

    int totbad = 0;		// Accumulate total sledge error count

    if (mpi_myid) {

      compute_table_worker();

      //
      // <Receive completed table from root>
      //

      for (kx=0; kx<=numk; kx++) {
	for (ky=0; ky<=kx; ky++) {
	  MPI_Recv(&mpi_buf[0], mpi_bufsz, MPI_PACKED, 0,
		   MPI_ANY_TAG, MPI_COMM_WORLD, &status);
    
	  mpi_unpack_table();      
	}
      }
      
    }
    else {			// BEGIN Root

      int worker = 0;
      int request_id = 1;

      if (!ReadH5Cache()) {

	kx = 0;
	ky = 0;

	while (kx<=numk) {

	  if (worker<mpi_numprocs-1) { // Send request to worker
	    worker++;
      
	    if (tbdbg)
	      std::cout << "Root sending orders to Worker " << worker 
			<< ": Kx=" << kx << ", Ky=" << ky << std::endl;

	    MPI_Send(&request_id, 1, MPI_INT, worker, 11, MPI_COMM_WORLD);
	    MPI_Send(&kx, 1, MPI_INT, worker, 11, MPI_COMM_WORLD);
	    MPI_Send(&ky, 1, MPI_INT, worker, 11, MPI_COMM_WORLD);
      
	    if (tbdbg)
	      std::cout << "Root gave orders to Worker " << worker 
			<< ": Kx=" << kx << ", Ky=" << ky << std::endl;

				// Increment counters
	    ky++;
	    if (ky>kx) {
	      kx++;
	      ky = 0;
	    }
	    
	  }

	  if (worker == mpi_numprocs-1 && kx<=numk) {
	  
	    //
	    // <Wait and receive>
	    //
#ifdef SLEDGE_THROW
	    int bad;		// Get sledge error count
	    MPI_Recv(&bad, 1, MPI_INT, MPI_ANY_TAG, 10,
		     MPI_COMM_WORLD, &status);
	    totbad += bad;
				// Get sledge computation result
	    int retid = status.MPI_SOURCE;
	    MPI_Recv(&mpi_buf[0], mpi_bufsz, MPI_PACKED, retid, 11,
		     MPI_COMM_WORLD, &status);
#else
	    MPI_Recv(&mpi_buf[0], mpi_bufsz, MPI_PACKED, MPI_ANY_SOURCE, 
		     MPI_ANY_TAG, MPI_COMM_WORLD, &status);
	    
	    int retid = status.MPI_SOURCE;
#endif

	    mpi_unpack_table();      

	    //
	    // <Send new request>
	    //

	    MPI_Send(&request_id, 1, MPI_INT, retid, 11, MPI_COMM_WORLD);
	    MPI_Send(&kx, 1, MPI_INT, retid, 11, MPI_COMM_WORLD);
	    MPI_Send(&ky, 1, MPI_INT, retid, 11, MPI_COMM_WORLD);
      
	    if (tbdbg)
	      std::cout << "Root gave orders to Worker " << retid 
			<< ": Kx=" << kx << ", Ky=" << ky << std::endl;

				// Increment counters
	    ky++;
	    if (ky>kx) {
	      kx++;
	      ky = 0;
	    }

	  }
	}
      
	//
	// <Wait for all workers to return>
	//
  
	while (worker) {
	
#ifdef SLEDGE_THROW
	  // Get sledge error count
	  int bad;
	  MPI_Recv(&bad, 1, MPI_INT, MPI_ANY_SOURCE, 10, MPI_COMM_WORLD,
		   &status);
	  totbad += bad;
	  // Get sledge computation result
	  int retid = status.MPI_SOURCE;
	  MPI_Recv(&mpi_buf[0], mpi_bufsz, MPI_PACKED, retid, 11,
		   MPI_COMM_WORLD, &status);
#else
	  MPI_Recv(&mpi_buf[0], mpi_bufsz, MPI_PACKED, MPI_ANY_SOURCE, 11,
		   MPI_COMM_WORLD, &status);
#endif
	  mpi_unpack_table();      

	  worker--;
	}

	if (cache) WriteH5Cache();

      }

      //
      // <Tell workers to continue>
      //

      request_id = -1;
      for (worker=1; worker < mpi_numprocs; worker++)
	MPI_Send(&request_id, 1, MPI_INT, worker, 11, MPI_COMM_WORLD);


      //
      // <Send table to workers>
      //

      for (kx=0; kx<=numk; kx++) {
	for (ky=0; ky<=kx; ky++) {
	  int position = mpi_pack_table(&table[kx][ky], kx, ky);
	  for (worker=1; worker < mpi_numprocs; worker++)
	    MPI_Send(&mpi_buf[0], position, MPI_PACKED, worker, 11, MPI_COMM_WORLD);
	}
      }


    } // END Root

#ifdef SLEDGE_THROW
    // Share total sledge error count with all nodes
    MPI_Bcast(&totbad, 1, MPI_INT, 0, MPI_COMM_WORLD);

    // Throw runtime error if sledge comtputation has errors
    if (totbad) {
      std::ostringstream sout;
      
      if (myid==0) {
	sout << std::endl << "SLGridSlab found " << totbad
	     << " tolerance errors in computing SL solutions." << std::endl
	     << "We suggest checking your model parameters to ensure a"
	     << std::endl
	     << "sufficient number of grid points that the relative difference"
	     << std::endl << "between field quantities is <= 0.3";
      } else {
	sout << std::endl << "sledge tolerance failure";
      }

      throw GenericError(sout.str(), __FILE__, __LINE__);
    }
#endif
  }
  else {

    if (!ReadH5Cache()) {

      for (kx=0; kx<=numk; kx++) {
	for (ky=0; ky<=kx; ky++) {
	  if (tbdbg) std::cerr << "Begin [" << kx << ", " << ky << "] . . ."
			       << std::endl;
	  compute_table(&(table[kx][ky]), kx, ky);
	  if (tbdbg) std::cerr << ". . . done" << std::endl;
	}
      }

      if (cache) WriteH5Cache();
    }
  }

  if (tbdbg)
    std::cerr << "Process " << myid << ": exiting constructor" << std::endl;
}


const string slab_cache_name = ".slgrid_slab_cache";


bool SLGridSlab::ReadH5Cache(void)
{
  if (!cache) return false;

  // First attempt to read the file
  //
  try {
    // Silence the HDF5 error stack
    //
    HighFive::SilenceHDF5 quiet;
    
    // Try opening the file as HDF5
    //
    HighFive::File h5file(slab_cache_name, HighFive::File::ReadOnly);
    
    // Try checking the rest of the parameters before reading arrays
    //
    auto checkInt = [&h5file](int value, std::string name)
    {
      int v; HighFive::Attribute vv = h5file.getAttribute(name); vv.read(v);
      if (value == v) return true;
      if (myid==0)
	std::cout << "---- SLGridSlab::ReadH5Cache: "
		  << "parameter " << name << ": wanted " << value
		  << " found " << v << std::endl;
      return false;
    };

    auto checkDbl = [&h5file](double value, std::string name)
    {
      double v; HighFive::Attribute vv = h5file.getAttribute(name); vv.read(v);
      if (fabs(value - v) < 1.0e-16) return true;
      if (myid==0)
	std::cout << "---- SLGridSlab::ReadH5Cache: "
		  << "parameter " << name << ": wanted " << value
		  << " found " << v << std::endl;
      return false;
    };

    auto checkStr = [&h5file](std::string value, std::string name)
    {
      std::string v; HighFive::Attribute vv = h5file.getAttribute(name); vv.read(v);
      if (value.compare(v)==0) return true;
      if (myid==0)
	std::cout << "---- SLGridSlab::ReadH5Cache: "
		  << "parameter " << name << ": wanted " << value
		  << " found " << v << std::endl;
      return false;
    };

    // For cache ID
    //
    std::string geometry("slab"), forceID("SLGridSlab");

    // ID check
    //
    if (not checkStr(geometry, "geometry"))  return false;
    if (not checkStr(forceID,  "forceID"))   return false;

    // Parameter check
    //
    if (not checkStr(type,     "type"))      return false;
    if (not checkInt(numk,     "numk"))      return false;
    if (not checkInt(nmax,     "nmax"))      return false;
    if (not checkInt(numz,     "numz"))      return false;
    if (not checkDbl(H,        "H"))         return false;
    if (not checkDbl(L,        "L"))         return false;
    if (not checkDbl(zmax,     "zmax"))      return false;
    if (not checkDbl(ZBEG,     "ZBEG"))      return false;
    if (not checkDbl(ZEND,     "ZEND"))      return false;

    // Harmonic order
    //
    auto harmonic = h5file.getGroup("Harmonic");

    // Create table instances
    //

    table = table_ptr_2D(new table_ptr_1D [numk+1]);
    for (int kx=0; kx<=numk; kx++)
      table[kx] = table_ptr_1D(new TableSlab [kx+1]);

    for (int kx=0; kx<=numk; kx++) {
      for (int ky=0; ky<=kx; ky++) {
	std::ostringstream sout;
	sout << kx << " " << ky;
	auto arrays = harmonic.getGroup(sout.str());
      
	arrays.getDataSet("ev").read(table[kx][ky].ev);
	arrays.getDataSet("ef").read(table[kx][ky].ef);
      }
    }
    
    if (myid==0)
      std::cout << "---- SLGridSlab::ReadH5Cache: "
		<< "successfully read basis cache <" << slab_cache_name
		<< ">" << std::endl;

    return true;
    
  } catch (HighFive::Exception& err) {
    if (myid==0)
      std::cerr << "---- SLGridSlab::ReadH5Cache: "
		<< "error reading <" << slab_cache_name << ">" << std::endl
		<< "---- SLGridSlab::ReadH5Cache: HDF5 error is <" << err.what()
		<< ">" << std::endl;
  }

  return false;
}



void SLGridSlab::WriteH5Cache(void)
{
  if (myid) return;

  try {

    // Check for new HDF5 file
    if (std::filesystem::exists(slab_cache_name)) {
      if (myid==0)
	std::cout << "---- SLGridSlab::WriteH5Cache cache file <"
		  << slab_cache_name << "> exists" << std::endl;
      try {
	std::filesystem::rename(slab_cache_name, slab_cache_name + ".bak");
      }
      catch(std::filesystem::filesystem_error const& ex) {
	std::ostringstream sout;
        sout << "---- SLGridSlab::WriteH5Cache write error: "
	     << "what():  " << ex.what()  << std::endl
	     << "path1(): " << ex.path1() << std::endl
	     << "path2(): " << ex.path2();
	throw GenericError(sout.str(), __FILE__, __LINE__, 12, true);
      }
      
      if (myid==0)
	std::cout << "---- SLGridSlab::WriteH5Cache: existing file backed up to <"
		  << slab_cache_name + ".bak>" << std::endl;
    }
    
    // Create a new hdf5 file
    //
    HighFive::File file(slab_cache_name,
			HighFive::File::ReadWrite | HighFive::File::Create);
    
    // For cache ID
    //
    std::string geometry("slab"), forceID("SLGridSlab");

    file.createAttribute<std::string>("geometry",  HighFive::DataSpace::From(geometry)).write(geometry);
    file.createAttribute<std::string>("forceID",   HighFive::DataSpace::From(forceID)).write(forceID);
      
    // Write parameters
    //
    file.createAttribute<std::string> ("type",     HighFive::DataSpace::From(type)).write(type);
    file.createAttribute<int>         ("numk",     HighFive::DataSpace::From(numk)).write(numk);
    file.createAttribute<int>         ("nmax",     HighFive::DataSpace::From(nmax)).write(nmax);
    file.createAttribute<int>         ("numz",     HighFive::DataSpace::From(numz)).write(numz);
    file.createAttribute<double>      ("H",        HighFive::DataSpace::From(H)).write(H);
    file.createAttribute<double>      ("L",        HighFive::DataSpace::From(L)).write(L);
    file.createAttribute<double>      ("zmax",     HighFive::DataSpace::From(ZBEG)).write(zmax);
    file.createAttribute<double>      ("ZBEG",     HighFive::DataSpace::From(ZBEG)).write(ZBEG);
    file.createAttribute<double>      ("ZEND",     HighFive::DataSpace::From(ZEND)).write(ZEND);
      
    // Harmonic order (for h5dump readability)
    //
    auto harmonic = file.createGroup("Harmonic");

    for (int kx=0; kx<=numk; kx++) {
      for (int ky=0; ky<=kx; ky++) {
	std::ostringstream sout;
	sout << kx << " " << ky;
	auto arrays = harmonic.createGroup(sout.str());
      
	arrays.createDataSet("ev",   table[kx][ky].ev);
	arrays.createDataSet("ef",   table[kx][ky].ef);
      }
    }
    
  } catch (HighFive::Exception& err) {
    std::cerr << err.what() << std::endl;
  }
    
  std::cout << "---- SLGridSlab::WriteH5Cache: "
	    << "wrote <" << slab_cache_name << ">" << std::endl;
  
  return ;
}


SLGridSlab::~SLGridSlab()
{
  // Nothing
}

double SLGridSlab::eigenvalue(int kx, int ky, int n) {
  return table.get()[kx].get()[ky].ev[n];
}

// Coordinate transformation member functions for tanh map
double SLGridSlab::TanhMap::z_to_xi  (double z)  { return tanh(z/H);       }
double SLGridSlab::TanhMap::xi_to_z  (double xi) { return H*atanh(xi);     }
double SLGridSlab::TanhMap::d_xi_to_z(double xi) { return (1.0 - xi*xi)/H; }

// Coordinate transformation member functions for sech map
double SLGridSlab::SechMap::z_to_xi  (double z)  { return z/sqrt(z*z + H*H); }
double SLGridSlab::SechMap::xi_to_z  (double xi) { return xi*H/sqrt(1.0 - xi*xi); }
double SLGridSlab::SechMap::d_xi_to_z(double xi) { return pow(1.0 - xi*xi, 1.5)/H; }

// Coordinate transformation member functions for linear map
double SLGridSlab::LinearMap::z_to_xi(double z)    { return z;   }
double SLGridSlab::LinearMap::xi_to_z(double xi)   { return xi;  }
double SLGridSlab::LinearMap::d_xi_to_z(double xi) { return 1.0; }

double SLGridSlab::get_pot(double x, int kx, int ky, int n, int which)
{
  int hold;

				// Flip sign for antisymmetric basis functions
  int sign=1;
  if (x<0 && 2*(n/2)!=n) sign=-1;
  x = fabs(x);

  if (which)			// Convert from z to x
    x = mM->z_to_xi(x);

  if (ky > kx) {
    hold = ky;
    ky = kx;
    kx = hold;
  }

  int indx = (int)( (x-xmin)/dxi );
  if (indx<0) indx = 0;
  if (indx>numz-2) indx = numz - 2;


  double x1 = (xi[indx+1] - x)/dxi;
  double x2 = (x - xi[indx])/dxi;
  

#ifdef USE_TABLE
  return (x1*table[kx][ky].ef(n, indx) + x2*table[kx][ky].ef(n, indx+1))/
    sqrt(table[kx][ky].ev[n]) * (x1*p0[indx] + x2*p0[indx+1]) * sign;
#else
  return (x1*table[kx][ky].ef(n, indx) + x2*table[kx][ky].ef(n, indx+1))/
    sqrt(table[kx][ky].ev[n]) * slab->pot(mM->xi_to_z(x)) * sign;
#endif
}


double SLGridSlab::get_dens(double x, int kx, int ky, int n, int which)
{
  int hold;

  int sign=1;
  if (x<0 && 2*(n/2)!=n) sign=-1;
  x = fabs(x);
  
  if (which)			// Convert from z to x
    x = mM->z_to_xi(x);

  if (ky > kx) {
    hold = ky;
    ky = kx;
    kx = hold;
  }

  int indx = (int)( (x-xmin)/dxi );
  if (indx<0) indx = 0;
  if (indx>numz-2) indx = numz - 2;


  double x1 = (xi[indx+1] - x)/dxi;
  double x2 = (x - xi[indx])/dxi;
  
#ifdef USE_TABLE
  return (x1*table[kx][ky].ef(n, indx) + x2*table[kx][ky].ef(n, indx+1)) *
    sqrt(table[kx][ky].ev[n]) * (x1*d0[indx] + x2*d0[indx+1]) * sign;
#else
  return (x1*table[kx][ky].ef(n, indx) + x2*table[kx][ky].ef(n, indx+1)) *
    sqrt(table[kx][ky].ev[n]) * slab->dens(mM->xi_to_z(x)) * sign;
#endif

}

double SLGridSlab::get_force(double x, int kx, int ky, int n, int which)
{
  int hold;

  int sign=1;
  if (x<0 && 2*(n/2)==n) sign = -1;
  x = fabs(x);

  if (which)			// Convert from z to x
    x = mM->z_to_xi(x);

  if (ky > kx) {
    hold = ky;
    ky = kx;
    kx = hold;
  }

  int indx = (int)( (x-xmin)/dxi );
  if (indx<1) indx = 1;
  if (indx>numz-2) indx = numz - 2;


  double p = (x - xi[indx])/dxi;
  
				// Use three point formula

				// Point -1: indx-1
				// Point  0: indx
				// Point  1: indx+1

  return mM->d_xi_to_z(x)/dxi * (
			     (p - 0.5)*table[kx][ky].ef(n, indx-1)*p0[indx-1]
			     -2.0*p*table[kx][ky].ef(n, indx)*p0[indx]
			     + (p + 0.5)*table[kx][ky].ef(n, indx+1)*p0[indx+1]
			     ) / sqrt(table[kx][ky].ev[n]) * sign;
}


void SLGridSlab::get_pot(Eigen::MatrixXd& mat, double x, int which)
{
  int sign=1, sign2;
  if (x<0) sign = -1;
  x = fabs(x);
  
  if (which)			// Convert from z to x
    x = mM->z_to_xi(x);

  int ktot = (numk+1)*(numk+2)/2;
  mat.resize(ktot+1, nmax);

  int indx = (int)( (x-xmin)/dxi );
  if (indx<0) indx = 0;
  if (indx>numz-2) indx = numz - 2;


  double x1 = (xi[indx+1] - x)/dxi;
  double x2 = (x - xi[indx])/dxi;
  

  int l=0;
  for (int kx=0; kx<=numk; kx++) {
    for (int ky=0; ky<=kx; ky++) {
      sign2 = 1;
      for (int n=0; n<nmax; n++) {
#ifdef USE_TABLE
	mat(l, n) = (x1*table[kx][ky].ef(n, indx) + x2*table[kx][ky].ef(n, indx+1))/
	  sqrt(table[kx][ky].ev[n]) * (x1*p0[indx] + x2*p0[indx+1]) * sign2;
#else
	mat(l, n) = (x1*table[kx][ky].ef(n, indx) + x2*table[kx][ky].ef(n, indx+1))/
	  sqrt(table[kx][ky].ev[n]) * slab->pot(mM->xi_to_z(x)) * sign2;
#endif
	sign2 *= sign;
      }
      l++;
    }    
  }

}


void SLGridSlab::get_dens(Eigen::MatrixXd& mat, double x, int which)
{
  int sign=1, sign2;
  if (x<0) sign = -1;
  x = fabs(x);
  
  if (which)			// Convert from z to x
    x = mM->z_to_xi(x);

  int ktot = (numk+1)*(numk+2)/2;
  mat.resize(ktot+1, nmax);

  int indx = (int)( (x-xmin)/dxi );
  if (indx<0) indx = 0;
  if (indx>numz-2) indx = numz - 2;

  double x1 = (xi[indx+1] - x)/dxi;
  double x2 = (x - xi[indx])/dxi;
  
  int l=0;
  for (int kx=0; kx<=numk; kx++) {
    for (int ky=0; ky<=kx; ky++) {
      sign2 = 1;
      for (int n=0; n<nmax; n++) {
#ifdef USE_TABLE
	mat(l, n) = (x1*table[kx][ky].ef(n, indx) + x2*table[kx][ky].ef(n, indx+1))*
	  sqrt(table[kx][ky].ev[n]) * (x1*d0[indx] + x2*d0[indx+1]) * sign2;
#else
	mat(l, n) = (x1*table[kx][ky].ef(n, indx) + x2*table[kx][ky].ef(n, indx+1))*
	  sqrt(table[kx][ky].ev[n]) * slab->dens(mM->xi_to_z(x)) * sign2;
#endif
	sign2 *= sign;
      }
      l++;
    }
  }

}


void SLGridSlab::get_force(Eigen::MatrixXd& mat, double x, int which)
{
  int sign=1, sign2;
  if (x<0) sign = -1;
  x = fabs(x);
  
  if (which)			// Convert from z to x
    x = mM->z_to_xi(x);

  int ktot = (numk+1)*(numk+2)/2;
  mat.resize(ktot+1, nmax);

  int indx = (int)( (x-xmin)/dxi );
  if (indx<1) indx = 1;
  if (indx>numz-2) indx = numz - 2;


  double p = (x - xi[indx])/dxi;
  double fac = mM->d_xi_to_z(x)/dxi;

  int l=0;
  for (int kx=0; kx<=numk; kx++) {
    for (int ky=0; ky<=kx; ky++) {
      sign2 = sign;
      for (int n=0; n<nmax; n++) {
	mat(l, n) = fac * (
			   (p - 0.5)*table[kx][ky].ef(n, indx-1)*p0[indx-1]
			   -2.0*p*table[kx][ky].ef(n, indx)*p0[indx]
			   + (p + 0.5)*table[kx][ky].ef(n, indx+1)*p0[indx+1]
			   ) / sqrt(table[kx][ky].ev[n]) * sign2;
	sign2 *= -sign;
      }
      l++;
    }
  }
  
}


void SLGridSlab::get_pot(Eigen::VectorXd& vec, double x, int kx, int ky, int which)
{
  int hold;

  int sign=1, sign2;
  if (x<0) sign = -1;
  x = fabs(x);
  
  if (which)			// Convert from z to x
    x = mM->z_to_xi(x);

  if (ky > kx) {
    hold = ky;
    ky = kx;
    kx = hold;
  }

  vec.resize(nmax);

  int indx = (int)( (x-xmin)/dxi );
  if (indx<0) indx = 0;
  if (indx>numz-2) indx = numz - 2;

  double x1 = (xi[indx+1] - x)/dxi;
  double x2 = (x - xi[indx])/dxi;
  

  sign2 = 1;
  for (int n=0; n<nmax; n++) {
#ifdef USE_TABLE
    vec[n] = (x1*table[kx][ky].ef(n, indx) + x2*table[kx][ky].ef(n, indx+1))/
      sqrt(table[kx][ky].ev[n]) * (x1*p0[indx] + x2*p0[indx+1]) * sign2;
#else
    vec[n] = (x1*table[kx][ky].ef(n, indx) + x2*table[kx][ky].ef(n, indx+1))/
      sqrt(table[kx][ky].ev[n]) * slab->pot(mM->xi_to_z(x)) * sign2;
#endif
    sign2 *= sign;
  }

}


void SLGridSlab::get_dens(Eigen::VectorXd& vec, double x, int kx, int ky, int which)
{
  int sign=1, sign2;
  if (x<0) sign = -1;
  x = fabs(x);
  
  if (which)			// Convert from z to x
    x = mM->z_to_xi(x);

  vec.resize(nmax);

  int indx = (int)( (x-xmin)/dxi );
  if (indx<0) indx = 0;
  if (indx>numz-2) indx = numz - 2;

  double x1 = (xi[indx+1] - x)/dxi;
  double x2 = (x - xi[indx])/dxi;
  

  sign2 = 1;
  for (int n=0; n<nmax; n++) {
#ifdef USE_TABLE
    vec[n] = (x1*table[kx][ky].ef(n, indx) + x2*table[kx][ky].ef(n, indx+1))*
      sqrt(table[kx][ky].ev[n]) * (x1*d0[indx] + x2*d0[indx+1]) * sign2;
#else
    vec[n] = (x1*table[kx][ky].ef(n, indx) + x2*table[kx][ky].ef(n, indx+1))*
      sqrt(table[kx][ky].ev[n]) * slab->dens(mM->xi_to_z(x)) * sign2;
#endif
    sign2 *= sign;
  }

}


void SLGridSlab::get_force(Eigen::VectorXd& vec, double x, int kx, int ky, int which)
{
  int hold;

  int sign=1, sign2;
  if (x<0) sign = -1;
  x = fabs(x);
  
  if (which)			// Convert from z to x
    x = mM->z_to_xi(x);

  if (ky > kx) {
    hold = ky;
    ky = kx;
    kx = hold;
  }

  vec.resize(nmax);

  int indx = (int)( (x-xmin)/dxi );
  if (indx<1) indx = 1;
  if (indx>numz-2) indx = numz - 2;


  double p = (x - xi[indx])/dxi;
  double fac = mM->d_xi_to_z(x)/dxi;

  sign2 = sign;
  for (int n=0; n<nmax; n++) {
    vec[n] = fac * (
		    (p - 0.5)*table[kx][ky].ef(n, indx-1)*p0[indx-1]
		    -2.0*p*table[kx][ky].ef(n, indx)*p0[indx]
		    + (p + 0.5)*table[kx][ky].ef(n, indx+1)*p0[indx+1]
		    ) / sqrt(table[kx][ky].ev[n]) * sign2;
    sign2 *= sign;
  }

}

void SLGridSlab::compute_table(struct TableSlab* table, int KX, int KY)
{
  double cons[8]    = {0.0, 0.0, 0.0, 0.0, 0.0, 0.0, ZBEG, zmax};
  double tol[6]     = {1.0e-6, 1.0e-7, 1.0e-6,1.0e-7, 1.0e-6,1.0e-7};
  logical type[8]   = {1, 0, 0, 0, 1, 0, 0, 0};
  logical endfin[2] = {1, 1};
  integer NUM, N;
  int VERBOSE=0;

#ifdef DEBUG_SLEDGE
  if (myid==0) VERBOSE = SLEDGE_VERBOSE;
#endif

  NUM = numz;
				// Divide total functions into symmetric
				// and antisymmetric (keeping equal number
				// of each or one more symmetric
  N = (int)( 0.5*nmax + 0.501);
  
  integer *iflag = new integer [nmax];
  integer *invec = new integer [nmax+3];

  double *t=0, *rho=0;
  double *ev    = new double [N];
  double *store = new double [26*(NUM+16)];
  double *xef   = new double [NUM+16];
  double *ef    = new double [NUM*N];
  double *pdef  = new double [NUM*N];
  double f, df;

  KKZ = 2.0*M_PI/L * sqrt((double)(KX*KX + KY*KY));

				// Even BC, inner has zero gradient
  f = slab->pot(cons[6]);
  cons[2] = -1.0/(f*f);

				// Outer
  if (KKZ>1.0e-4) {
    f = slab->pot(cons[7]);
    df = slab->dpot(cons[7]);
    cons[4] = (df + KKZ*f)*f;
  }
  cons[5] = 1.0;
  //  cons[5] = 1.0/(f*f);

  //
  //     Initialize the vector INVEC(*):
  //       estimates for the eigenvalues/functions specified
  //

  invec[0] = VERBOSE;		// little printing (1), no printing (0)
  invec[1] = 3;			// spectrum is ignored
  invec[2] = N;			// estimates for N eigenvalues/functions

  for (int i=0; i<N; i++) invec[3+i] = i;

  //
  //     Set the JOB(*) vector:
  //        estimate both eigenvalues and eigenvectors,
  //        don't estimate the spectral density function,
  //        classify,
  //        let SLEDGE choose the initial mesh
  //
  logical job[5] = {0,1,0,1,0};

  //
  //     Output mesh
  //
  for (int i=0; i<NUM; i++) xef[i] = z[i];

  //     
  //     Open file for output.
  //
  sl_dim = 1;

  sledge_(job, cons, endfin, invec, tol, type, ev, &NUM, xef, ef, pdef,
	   t, rho, iflag, store);

  //
  //     Check for errors
  //
#ifdef SLEDGE_THROW
  unsigned bad = 0;		// Number of non-zero iflag values
  for (int i=0; i<N; i++) {
    if (iflag[i] != 0) bad++;
  }

  std::ostringstream sout;	// Runtime error message

  // Print info if we have sledge errors and throw a runtime
  // exception
  if (bad>0) {

    if (myid==0) {

      std::cout.precision(6);
      std::cout.setf(ios::scientific);
      std::cout << std::left;
      
      std::cout << std::endl
		<< "Tolerance errors in Sturm-Liouville solver for Kx=" << KX
		<< " Ky=" << KY << ", even" <<  std::endl << std::endl;

      std::cout << std::setw(15) << "order"
		<< std::setw(15) << "eigenvalue"
		<< std::setw(40) << "condition"
		<< std::endl
		<< std::setw(15) << "-----"
		<< std::setw(15) << "----------"
		<< std::setw(40) << "---------"
		<< std::endl;
      
      for (int i=0; i<N; i++) {
	std::cout << std::setw(15) << invec[3+i] 
		  << std::setw(15) << ev[i]
		  << std::setw(40) << sledge_error(iflag[i])
		  << std::endl;
      }
      std::cout << std::endl;
      
      sout << std::endl
	   << "SLGridSlab found " << bad
	   << " tolerance errors in computing SL solutions." << std::endl
	   << "We suggest checking your model parameters to ensure a"
	   << std::endl
	   << "sufficient number of grid points that the relative difference"
	   << std::endl << "between field quantities is <= 0.3";

      throw GenericError(sout.str(), __FILE__, __LINE__);
    } else {
      throw GenericError("sledge errors", __FILE__, __LINE__);
    }
  }
#endif

  //
  //     Print results:
  //
  if (tbdbg) {

    std::cout.precision(6);
    std::cout.setf(ios::scientific);

    std::cout << "Even:" << std::endl;
    for (int i=0; i<N; i++) {
      std::cout << std::setw(15) << invec[3+i] 
		<< std::setw(15) << ev[i]
		<< std::setw( 5) << iflag[i]
		<< std::endl;
  
      if (VERBOSE) {

	if (iflag[i] > -10) {
	  cout << std::setw(14) << "x"
	       << std::setw(25) << "u(x)"
	       << std::setw(25) << "(pu`)(x)"
	       << std::endl;
	  int k = NUM*i;
	  for (int j=0; j<NUM; j++) {
	    std::cout << std::setw(25) << xef[j]
		      << std::setw(25) << ef[j+k]
		      << std::setw(25) << pdef[j+k]
		      << std::endl;
	  }
	}
	
      }

    }
  }


  // Allocate memory for table
  //
  table->ev.resize(nmax);
  table->ef.resize(nmax, numz);

  // Load table
  //
  for (int i=0; i<N; i++) table->ev[i*2] = ev[i];

  // Choose sign conventions for the ef table
  //
  {
    int nfid = std::min<int>(nevsign, NUM) - 1;
    Eigen::VectorXi sgn = Eigen::VectorXi::Ones(N);
    for (int j=0; j<N; j++) {
      if (ef[j*NUM+nfid]<0.0) sgn(j) = -1;
    }
  
    for (int i=0; i<numz; i++) {
      for (int j=0; j<N; j++) 
	table->ef(j*2, i) = ef[j*NUM+i] * sgn(j);
    }
  }

				// Odd BC, Inner zero value
  cons[0] = 1.0;
  cons[2] = 0.0;

				// Redo to get antisymmetric functions

  sledge_(job, cons, endfin, invec, tol, type, ev, &NUM, xef, ef, pdef,
	  t, rho, iflag, store);

  //
  //     Check for errors
  //
#ifdef SLEDGE_THROW
  bad = 0;			// Accumulate sledge error count
  for (int i=0; i<N; i++) {
    if (iflag[i] != 0) bad++;
  }

  sout.str("");

  // Print info if we have errors and throw a runtime exception
  if (bad>0) {

    if (myid==0) {

      std::cout.precision(6);
      std::cout.setf(ios::scientific);
      std::cout << std::left;
    
      std::cout << std::endl
		<< "Tolerance errors in Sturm-Liouville solver for Kx=" << KX
		<< " Ky=" << KY << ", odd" <<  std::endl;

      std::cout << std::setw(15) << "order"
		<< std::setw(15) << "eigenvalue"
		<< std::setw(40) << "condition"
		<< std::endl
		<< std::setw(15) << "-----"
		<< std::setw(15) << "----------"
		<< std::setw(40) << "---------"
		<< std::endl;
      
      for (int i=0; i<N; i++) {
	std::cout << std::setw(15) << invec[3+i] 
		  << std::setw(15) << ev[i]
		  << std::setw(40) << sledge_error(iflag[i])
		  << std::endl;
      }
      std::cout << std::endl;
      
      sout << std::endl << "SLGridSph found " << bad
	   << " tolerance errors in computing SL solutions." << std::endl
	   << "We suggest checking your model parameters to ensure a"
	   << std::endl
	   << "sufficient number of grid points that the relative difference"
	   << std::endl << "between field quantities is <= 0.3";

      throw GenericError(sout.str(), __FILE__, __LINE__);
    } else {
      throw GenericError("sledge errors", __FILE__, __LINE__);
    }
  }
#endif

  //
  //     Print results:
  //
  if (tbdbg) {
  
    std::cout.precision(6);
    std::cout.setf(ios::scientific);
    
    std::cout << "Odd:" << std::endl;
    for (int i=0; i<N; i++) {
      std::cout << std::setw(15) << invec[3+i] 
		<< std::setw(15) << ev[i]
		<< std::setw( 5) << iflag[i]
		<< std::endl;
  
      if (VERBOSE) {

	if (iflag[i] > -10) {
	  cout << std::setw(14) << "x"
	       << std::setw(25) << "u(x)"
	       << std::setw(25) << "(pu`)(x)"
	       << std::endl;
	  int k = NUM*i;
	  for (int j=0; j<NUM; j++) {
	    std::cout << std::setw(25) << xef[j]
		      << std::setw(25) << ef[j+k]
		      << std::setw(25) << pdef[j+k]
		      << std::endl;
	  }
	}
      }
    }
  }

				// Load table

  N = nmax - N;

  for (int i=0; i<N; i++) table->ev[i*2+1] = ev[i];

  // Choose sign conventions for the ef table
  //
  {
    int nfid = std::min<int>(nevsign, NUM) - 1;
    Eigen::VectorXi sgn = Eigen::VectorXi::Ones(N);
    for (int j=0; j<N; j++) {
      if (ef[j*NUM+nfid]<0.0) sgn(j) = -1;
    }
  
    for (int i=0; i<numz; i++) {
      for (int j=0; j<N; j++) 
	table->ef(j*2+1, i) = ef[j*NUM+i] * sgn(j);
    }
  }

				// Correct for symmetrizing
  table->ef *= 7.071067811865475e-01;

  table->kx = KX;
  table->ky = KY;

  delete [] iflag;
  delete [] invec;
  delete [] ev;
  delete [] store;
  delete [] xef;
  delete [] ef;
  delete [] pdef;
}


void SLGridSlab::init_table(void)
{
  xi.resize(numz);
  z. resize(numz);
  p0.resize(numz);
  d0.resize(numz);

  xmin = mM->z_to_xi( ZBEG);
  xmax = mM->z_to_xi( zmax);
  dxi = (xmax-xmin)/(numz-1);


  for (int i=0; i<numz; i++) {
    xi[i] = xmin + dxi*i;
    z[i]  = mM->xi_to_z(xi[i]);
    p0[i] = slab->pot(z[i]);
    d0[i] = slab->dens(z[i]);
  }

}


void SLGridSlab::compute_table_worker(void)
{
  double tol[6]     = {1.0e-6,1.0e-7,  1.0e-6,1.0e-7,  1.0e-6,1.0e-7};
  logical type[8]   = {1, 0, 0, 0, 1, 0, 0, 0};
  logical endfin[2] = {1, 1};
  double cons[8];
  int VERBOSE=0;
  integer NUM;
  
  struct TableSlab table;
  int KX, KY, N;

#ifdef DEBUG_SLEDGE
  if (myid==0) VERBOSE = SLEDGE_VERBOSE;
#endif

#ifdef DEBUG
  std::cout << "Worker " << mpi_myid << " begins . . ." << std::endl;
#endif

  //
  // <Wait for orders>
  //
  
  int request_id;

  while(1) {

    MPI_Recv(&request_id, 1, 
	     MPI_INT, 0, MPI_ANY_TAG, MPI_COMM_WORLD, &status);

    if (request_id < 0) break;	// Good-bye

    MPI_Recv(&KX, 1, 
	     MPI_INT, 0, MPI_ANY_TAG, MPI_COMM_WORLD, &status);

    MPI_Recv(&KY, 1, 
	     MPI_INT, 0, MPI_ANY_TAG, MPI_COMM_WORLD, &status);

    if (tbdbg)
      std::cout << "Worker " << mpi_myid << ": ordered to compute Kx, Ky = "
		<< KX << ", " << KY << "" << std::endl;

    cons[0] = cons[1] = cons[2] = cons[3] = cons[4] = cons[5] = 0.0;
    cons[6] =  ZBEG;
    cons[7] =  zmax;
    NUM = numz;
				// Divide total functions into symmetric
				// and antisymmetric (keeping equal number
				// of each or one more symmetric
    N = (int)( 0.5*nmax + 0.501);

    integer *iflag = new integer [nmax];
    integer *invec = new integer [nmax+3];

    double *t=0, *rho=0;
    double *ev    = new double [N];
    double *store = new double [26*(NUM+16)];
    double *xef   = new double [NUM+16];
    double *ef    = new double [NUM*N];
    double *pdef  = new double [NUM*N];
    double f, df;

    KKZ = 2.0*M_PI/L * sqrt((double)(KX*KX + KY*KY));

				// Even BC, inner has zero gradient
    f = slab->pot(cons[6]);
    cons[2] = -1.0/(f*f);

    if (tbdbg) {
      std::cout << "Worker " << mpi_myid << ": Kx, Ky = " 
		<< KX << ", " << KY << " [Even inputs]" << std::endl;
      for (int n=0; n<8; n++)
	std::cout << std::setw(5) << n << std::setw(15) << cons[n] << endl;
    }
				// Outer
    if (KKZ>1.0e-4) {
      f = slab->pot(cons[7]);
      df = slab->dpot(cons[7]);
      cons[4] = (df + KKZ*f)*f;
    }
    cons[5] = 1.0;

    //
    //     Initialize the vector INVEC(*):
    //       estimates for the eigenvalues/functions specified
    //

    invec[0] = VERBOSE;		// little printing (1), no printing (0)
    invec[1] = 3;		// spectrum is ignored
    invec[2] = N;		// estimates for N eigenvalues/functions

    for (int i=0; i<N; i++) invec[3+i] = i;

    //
    //     Set the JOB(*) vector:
    //        estimate both eigenvalues and eigenvectors,
    //        don't estimate the spectral density function,
    //        don't classify,
    //        we choose the initial mesh
    //
    logical job[5] = {0,1,0,1,0};

    //
    //     Output mesh
    //
    for (int i=0; i<NUM; i++) xef[i] = z[i];

    //     
    //     Open file for output.
    //
    sl_dim = 1;

    sledge_(job, cons, endfin, invec, tol, type, ev, &NUM, xef, ef, pdef,
	    t, rho, iflag, store);
    //
    //     Print results:
    //
    
    //
    //     Check for errors
    //
#ifdef SLEDGE_THROW
    unsigned bad = 0;		// Number of non-zero iflag values
    for (int i=0; i<N; i++) {
      if (iflag[i] != 0) bad++;
    }
    
    std::ostringstream sout;	// Runtime error message

    // Print info if we have errors.  Number of errors will be sent to
    // and accumulated by the root node.
    if (bad>0) {

      std::cout.precision(6);
      std::cout.setf(ios::scientific);
      std::cout << std::left;
    
      std::cout << std::endl
		<< "Tolerance errors in Sturm-Liouville solver for Kx=" << KX
		<< " Ky=" << KY << ", even" <<  std::endl << std::endl;

      std::cout << std::setw(15) << "order"
		<< std::setw(15) << "eigenvalue"
		<< std::setw(40) << "condition"
		<< std::endl
		<< std::setw(15) << "-----"
		<< std::setw(15) << "----------"
		<< std::setw(40) << "---------"
		<< std::endl;

      for (int i=0; i<N; i++) {
	std::cout << std::setw(15) << invec[3+i] 
		  << std::setw(15) << ev[i]
		  << std::setw(40) << sledge_error(iflag[i])
		  << std::endl;
      }
      std::cout << std::endl;
    }
#endif

    if (tbdbg) {
      std::cout << "Worker " << mpi_myid << ": computed Kx, Ky = " 
		<< KX << ", " << KY << " [Even]" << std::endl;

      std::cout.precision(6);
      std::cout.setf(ios::scientific);

      for (int i=0; i<N; i++) {
	std::cout << std::setw(15) << invec[3+i] 
		  << std::setw(15) << ev[i]
		  << std::setw( 5) << iflag[i]
		  << std::endl;
  
	if (VERBOSE) {

	  if (iflag[i] > -10) {
	    std::cout << std::setw(14) << "x"
		      << std::setw(25) << "u(x)"
		      << std::setw(25) << "(pu`)(x)"
		      << std::endl;
	    int k = NUM*i;
	    for (int j=0; j<NUM; j++) {
	      std::cout << std::setw(25) << xef[j]
			<< std::setw(25) << ef[j+k]
			<< std::setw(25) << pdef[j+k]
			<< std::endl;
	    }
	  }
	  
	}
	
      }
  
    }

    // Allocate memory for table (even and odd)
    //
    table.ev.resize(nmax);
    table.ef.resize(nmax, numz);
    
    // Load table (even)
    //
    for (int i=0; i<N; i++) table.ev[i*2] = ev[i];

    // Choose sign conventions for the ef table
    //
    {
      int nfid = std::min<int>(nevsign, NUM) - 1;
      Eigen::VectorXi sgn = Eigen::VectorXi::Ones(N);
      for (int j=0; j<N; j++) {
	if (ef[j*NUM+nfid]<0.0) sgn(j) = -1;
      }
  
      for (int i=0; i<numz; i++) {
	for (int j=0; j<N; j++) 
	  table.ef(j*2, i) = ef[j*NUM+i] * sgn(j);
      }
    }

				// Odd BC, inner zero value
    cons[0] = 1.0;
    cons[2] = 0.0;

				// Redo to get antisymmetric functions

    sledge_(job, cons, endfin, invec, tol, type, ev, &NUM, xef, ef, pdef,
	    t, rho, iflag, store);

    //
    //     Check for errors
    //
#ifdef SLEDGE_THROW
    // Accumulate sledge error count
    bad = 0;
    for (int i=0; i<N; i++) {
      if (iflag[i] != 0) bad++;
    }

    sout.str("");

    if (bad>0) {

      // Print sledge error info and throw runtime exception
      if (myid==0) {

	std::cout.precision(6);
	std::cout.setf(ios::scientific);
	std::cout << std::left;
    
	std::cout << std::endl
		  << "Tolerance errors in Sturm-Liouville solver for Kx=" << KX
		  << " Ky=" << KY << ", odd" <<  std::endl << std::endl;

	std::cout << std::setw(15) << "order"
		  << std::setw(15) << "eigenvalue"
		  << std::setw(40) << "condition"
		  << std::endl
		  << std::setw(15) << "-----"
		  << std::setw(15) << "----------"
		  << std::setw(40) << "---------"
		  << std::endl;
	
	for (int i=0; i<N; i++) {
	  std::cout << std::setw(15) << invec[3+i] 
		    << std::setw(15) << ev[i]
		    << std::setw(40) << sledge_error(iflag[i])
		    << std::endl;
	}
	std::cout << std::endl;
      
	sout << std::endl << "SLGridSlab found " << bad
	     << " tolerance errors in computing SL solutions" << std::endl
	     << "We suggest checking your model parameters to ensure a"
	     << std::endl
	     << "sufficient number of grid points that the relative difference"
	     << std::endl << "between field quantities is <= 0.3";

	throw GenericError(sout.str(), __FILE__, __LINE__);
      } else {
	throw GenericError("sledge errors", __FILE__, __LINE__);
      }
    }
#endif

    //     Print results:
    //
    if (tbdbg) {

      std::cout << "Worker " << mpi_myid << ": computed Kx, Ky = " 
		<< KX << ", " << KY << " [Odd]" << std::endl;

      std::cout.precision(6);
      std::cout.setf(ios::scientific);

      for (int i=0; i<N; i++) {
	std::cout << std::setw(15) << invec[3+i] 
		  << std::setw(15) << ev[i]
		  << std::setw( 5) << iflag[i]
		  << std::endl;
	
	if (VERBOSE) {

	  if (iflag[i] > -10) {
	    std::cout << std::setw(14) << "x"
		      << std::setw(25) << "u(x)"
		      << std::setw(25) << "(pu`)(x)"
		      << std::endl;
	    int k = NUM*i;
	    for (int j=0; j<NUM; j++) {
	      std::cout << std::setw(25) << xef[j]
			<< std::setw(25) << ef[j+k]
			<< std::setw(25) << pdef[j+k]
			<< std::endl;
	    }
	  }
	}
      }
    }
  
    // Load table (odd)
    //
    N = nmax - N;

    for (int i=0; i<N; i++) table.ev[i*2+1] = ev[i];

    // Choose sign conventions for the ef table
    //
    {
      int nfid = std::min<int>(nevsign, NUM) - 1;
      Eigen::VectorXi sgn = Eigen::VectorXi::Ones(N);
      for (int j=0; j<N; j++) {
	if (ef[j*NUM+nfid]<0.0) sgn(j) = -1;
      }
  
      for (int i=0; i<numz; i++) {
	for (int j=0; j<N; j++) 
	  table.ef(j*2+1, i) = ef[j*NUM+i] * sgn(j);
      }
    }

				// Correct for symmetrizing
    table.ef *= 7.071067811865475e-01;
    table.kx = KX;
    table.ky = KY;
  
#ifdef SLEDGE_THROW
    // Send sledge error count to root
    MPI_Send(&bad, 1, MPI_INT, 0, 10, MPI_COMM_WORLD);
#endif
    int position = mpi_pack_table(&table, KX, KY);
    MPI_Send(&mpi_buf[0], position, MPI_PACKED, 0, 11, MPI_COMM_WORLD);

    if (tbdbg)
      std::cout << "Worker " << mpi_myid << ": sent to root Kx, Ky = "
		<< KX << ", " <<  KY << std::endl;

    delete [] iflag;
    delete [] invec;
    delete [] ev;
    delete [] store;
    delete [] xef;
    delete [] ef;
    delete [] pdef;
  }

}


void SLGridSlab::mpi_setup(void)
{
				// Get MPI id

  MPI_Comm_rank(MPI_COMM_WORLD, &mpi_myid);
  MPI_Comm_size(MPI_COMM_WORLD, &mpi_numprocs);


  // Setup for pack/unpack

  int buf1, buf2;
  MPI_Pack_size( 1, MPI_INT, MPI_COMM_WORLD, &buf1);
  MPI_Pack_size( 1, MPI_DOUBLE, MPI_COMM_WORLD, &buf2);

  mpi_bufsz = 2*buf1 +		// kx, ky
    nmax*buf2 +			// ev
    nmax*numz*buf2 ;		// ef

  mpi_buf = std::shared_ptr<char[]>(new char [mpi_bufsz]);
}


int SLGridSlab::mpi_pack_table(struct TableSlab* table, int kx, int ky)
{
  int position = 0;

  MPI_Pack( &kx, 1, MPI_INT, &mpi_buf[0], mpi_bufsz, 
	    &position, MPI_COMM_WORLD);

  MPI_Pack( &ky, 1, MPI_INT, &mpi_buf[0], mpi_bufsz, 
	    &position, MPI_COMM_WORLD);

  for (int j=0; j<nmax; j++)
    MPI_Pack( &table->ev[j], 1, MPI_DOUBLE, &mpi_buf[0], mpi_bufsz, 
	      &position, MPI_COMM_WORLD);

  for (int j=0; j<nmax; j++)
    for (int i=0; i<numz; i++) {
      MPI_Pack( &table->ef(j, i), 1, MPI_DOUBLE, &mpi_buf[0], mpi_bufsz, 
		&position, MPI_COMM_WORLD);
    }

  return position;
}


void SLGridSlab::mpi_unpack_table(void)
{
  int length, position = 0;
  int kx, ky;
  int retid = status.MPI_SOURCE;

  MPI_Get_count( &status, MPI_PACKED, &length);


  MPI_Unpack( &mpi_buf[0], length, &position, &kx, 1, MPI_INT,
	      MPI_COMM_WORLD);

  MPI_Unpack( &mpi_buf[0], length, &position, &ky, 1, MPI_INT,
	      MPI_COMM_WORLD);

  if (tbdbg)
    std::cout << "Process " << mpi_myid << ": unpacking table entry from Process "
	      << retid << ": kx=" << kx << ", " << ky << "" << std::endl;

  table[kx][ky].kx = kx;
  table[kx][ky].ky = ky;
  table[kx][ky].ev.resize(nmax);
  table[kx][ky].ef.resize(nmax, numz);

  for (int j=0; j<nmax; j++)
    MPI_Unpack( &mpi_buf[0], length, &position, &table[kx][ky].ev[j], 1, MPI_DOUBLE,
		MPI_COMM_WORLD);

  for (int j=0; j<nmax; j++)
    for (int i=0; i<numz; i++)
      MPI_Unpack( &mpi_buf[0], length, &position, &table[kx][ky].ef(j, i), 1, 
		  MPI_DOUBLE, MPI_COMM_WORLD);
}

std::unique_ptr<SLGridSlab::CoordMap> SLGridSlab::CoordMap::factory
(CoordMapTypes type, double H)
{
  if (type == CoordMapTypes::Tanh) {
    return std::make_unique<TanhMap>(H);
  }
  else if (type == CoordMapTypes::Sech) {
    return std::make_unique<SechMap>(H);
  }
  else if (type == CoordMapTypes::Linear) {
    return std::make_unique<LinearMap>(H);
  }
  else {
    throw std::runtime_error("CoordMap::factory: invalid map type");
  }
}


std::vector<Eigen::MatrixXd> SLGridSlab::orthoCheck(int num)
{
  // Gauss-Legendre knots and weights
  LegeQuad lw(num);

  // Get the scaled coordinate limits
  double ximin = mM->z_to_xi(-zmax);
  double ximax = mM->z_to_xi( zmax);

  // Initialize the return matrices
  std::vector<Eigen::MatrixXd> ret((numk+1)*(numk+2)/2);
  for (auto & v : ret) {
    v.resize(nmax, nmax);
    v.setZero();
  }

  int nthrds = omp_get_max_threads();
  std::vector<Eigen::VectorXd> vpot(nthrds), vden(nthrds);
  for (auto & v : vpot) v.resize(nmax);
  for (auto & v : vden) v.resize(nmax);

#pragma omp parallel for
  for (int i=0; i<num; i++) {
    int tid = omp_get_thread_num();
    double x = ximin + (ximax - ximin)*lw.knot(i);
	    
    int indx = 0;
    for (int kx=0; kx<=numk; kx++) {
      for (int ky=0; ky<=kx; ky++, indx++) {
	
	get_pot (vpot[tid], x, kx, ky, 0);	
	get_dens(vden[tid], x, kx, ky, 0);
	
	for (int n1=0; n1<nmax; n1++)
	  for (int n2=0; n2<nmax; n2++)
#pragma omp critical
	    ret[indx](n1, n2) += -vpot[tid](n1)*vden[tid](n2)/mM->d_xi_to_z(x)
	      *(ximax - ximin)*lw.weight(i);
      }
    }
  }
  
  return ret;
}

//======================================================================
//======================================================================
//======================================================================


extern "C" int coeff_(doublereal* x, doublereal* px, doublereal* qx, 
	doublereal* rx)
{
  double f,rho;

  if (sl_dim==1) {		// 1-d slab

    f   = slab->pot(*x);
    rho = slab->dens(*x);

    *px = f*f;
    *qx = (KKZ*KKZ*f - rho)*f;
    *rx = -rho*f;
  }
  else {			// Spherical

    f   = sphpot(*x);
    rho = sphdens(*x);

    *px = (*x)*(*x)*f*f;
    *qx = (L2*f - rho*(*x)*(*x))*f;
    *rx = -rho*(*x)*(*x)*f;

    if (*px<=0) 
      std::cerr << "Process " << myid << ": "
		<< "px<=0: x=" << *x << " f=" << f << "" << std::endl;
    if (*rx<=0)
      std::cerr << "Process " << myid << ": "
		<< "rx<=0: x=" << *x << " f=" << f << " rho=" << rho <<  "" << std::endl;

  }

  return 0;
}

