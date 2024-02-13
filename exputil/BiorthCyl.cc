
#include <filesystem>
#include <algorithm>
#include <iostream>
#include <iomanip>
#include <sstream>
#include <memory>
#include <limits>
#include <string>

#include <yaml-cpp/yaml.h>	// YAML support

#include <Eigen/Eigenvalues>	// Eigen3

#ifdef HAVE_OMP_H
#include <omp.h>		// For multithreading basis construction
#endif

#include <BiorthCyl.H>		// Definition for this class
#include <PotRZ.H>		// Hankel computation for potential
#include <EXPException.H>	// For GenericError

#include <numerical.H>
#include <Progress.H>		// Progress bar
#include <gaussQ.H>

#include <libvars.H>
using namespace __EXP__;	// For reference to n-body globals

// Constructor
BiorthCyl::BiorthCyl(const YAML::Node& conf) : conf(conf)
{
  // Read and assign parameters
  //
  try {
    if (conf["Mmax"])        mmax = conf["Mmax"].as<int>();
    else                     mmax = 6;
    
    if (conf["Lmax"] and not conf["Mmax"])
                             mmax = conf["Lmax"].as<int>();
    
    if (conf["nmaxfid"])     nmaxfid = conf["nmaxfid"].as<int>();
    else                     nmaxfid = 40;

    if (conf["nmax"])        nmax = conf["nmax"].as<int>();    
    else                     nmax = nmaxfid;			       
    			                                           
    if (conf["numr"])        numr = conf["numr"].as<int>();	       
    else                     numr = 2000;			       
    			                                           
    if (conf["numx"])        numx = conf["numx"].as<int>();	       
    else                     numx = 512;			       
    			                                           
    if (conf["numy"])        numy = conf["numy"].as<int>();	       
    else                     numy = 256;			       
    			                                           
    if (conf["knots"])       knots = conf["knots"].as<int>();	       
    else                     knots = 1000;			       
    			                                           
    if (conf["NQDHT"])       NQDHT = conf["NQDHT"].as<int>();	       
    else                     NQDHT = 512;
    			                                           
    if (conf["rcylmin"])     rcylmin = conf["rcylmin"].as<double>();     
    else                     rcylmin = 0.0;			       
    			                                           
    if (conf["rcylmax"])     rcylmax = conf["rcylmax"].as<double>();     
    else                     rcylmax = 10.0;			       
    			                                           
    if (conf["cmapR"])       cmapR = conf["cmapR"].as<int>();      
    else                     cmapR = 1;			       
    			                                           
    if (conf["cmapZ"])       cmapZ = conf["cmapZ"].as<int>();      
    else                     cmapZ = 1;			       
    			                                           
    if (conf["scale"])       scale = conf["scale"].as<double>(); 
    else                     scale = 0.01;
    			                                           
    if (conf["cachename"])   cachename = conf["cachename"].as<std::string>();
    else                     cachename = default_cache;
    
    // Add output directory and runtag
    cachename = outdir + cachename;

    if (conf["verbose"])     verbose = true;
    else                     verbose = false;

    if (conf["mmin"])        mmin = conf["mmin"].as<int>();
    else                     mmin = 0;
    
    if (conf["mlim"])        mlim = conf["mlim"].as<int>();
    else                     mlim = mmax;
    
    if (conf["nmin"])        nmin = conf["nmin"].as<int>();
    else                     nmin = 0;
    
    if (conf["nlim"])        nlim = conf["nlim"].as<int>();
    else                     nlim = nmax;

    if (conf["EVEN_M"])      EVEN_M = conf["EVEN_M"].as<bool>();
    else                     EVEN_M = false;

    if (conf["diskconf"])    diskconf  = conf["diskconf"];
    else throw std::runtime_error("BiorthCyl: you must specify the diskconf stanza");
  }
  catch (YAML::Exception & error) {
    if (myid==0) std::cout << "Error parsing parameters in BiorthCyl: "
			   << error.what() << std::endl
			   << std::string(60, '-') << std::endl
			   << "Config node"        << std::endl
			   << std::string(60, '-') << std::endl
			   << conf
			   << std::string(60, '-') << std::endl;
    MPI_Finalize();
    exit(-1);
  }

  geometry = "cylinder";
  forceID  = "BiorthCyl";

  // Hardwired for testing.  Could be parameters.
  //
  logr   = false;
  biorth = "bess";
  
  initialize();

  if (not ReadH5Cache()) create_tables();

}

void BiorthCyl::initialize()
{
  // Create storage and mapping constants
  //
  xmin    = r_to_xi(rcylmin*scale);
  xmax    = r_to_xi(rcylmax*scale);
  dx      = (xmax - xmin)/(numx-1);

  ymin    = z_to_yi(0.0);
  ymax    = z_to_yi(rcylmax*scale);
  dy      = (ymax - ymin)/(numy-1);

  dens   .resize(mmax+1);
  pot    .resize(mmax+1);
  rforce .resize(mmax+1);
  zforce .resize(mmax+1);

  for (int m=0; m<=mmax; m++) {

    dens  [m].resize(nmax);
    pot   [m].resize(nmax);
    rforce[m].resize(nmax);
    zforce[m].resize(nmax);

    for (int n=0; n<nmax; n++) {
      dens  [m][n].resize(numx, numy);
      pot   [m][n].resize(numx, numy);
      rforce[m][n].resize(numx, numy);
      zforce[m][n].resize(numx, numy);
    }
  }
}

void BiorthCyl::create_tables()
{
  emp = EmpCyl2d(mmax, nmaxfid, nmax, knots, numr,
		 rcylmin*scale, rcylmax*scale, 
		 scale, cmapR, logr, diskconf, biorth);

  if (conf["basis"]) emp.basisTest(true);

  std::shared_ptr<progress::progress_display> progress;

  if (myid==0) {
    std::cout << "---- BiorthCyl::create_tables: creating 3d basis..."
	      << std::endl;
    if (verbose) {
      int total = (mmax+1)*nmax*numx;
      progress = std::make_shared<progress::progress_display>(total);
    }
  }


  // Pack basis grids
  //
  for (int m=0; m<=mmax; m++) {

    // Potential instance
    //
    PotRZ potrz(rcylmax*scale, NQDHT, m);

    for (int n=0; n<nmax; n++) {

      dens  [m][n].setZero();
      pot   [m][n].setZero();
      rforce[m][n].setZero();
      zforce[m][n].setZero();

      // Create the functor
      //
      auto func = [&, this](double R)
      {
	return emp.get_dens(R, m, n);
      };

      int count = 0;

      for (int ix=0; ix<numx; ix++) {
	double r = xi_to_r(xmin + dx*ix);

	for (int iy=0; iy<numy; iy++) {
	  double z = yi_to_z(ymin + dy*iy);
	  
	  if (count++ % numprocs == myid) {
	    if (fabs(z)<1.0e-6) dens[m][n](ix, iy) = -emp.get_dens(r, m, n);
	    pot   [m][n](ix, iy) = potrz(r, z, func, PotRZ::Field::potential);
	    rforce[m][n](ix, iy) = potrz(r, z, func, PotRZ::Field::rforce   );
	    zforce[m][n](ix, iy) = potrz(r, z, func, PotRZ::Field::zforce   );
	  }
	}
	// END: y loop
	if (verbose and myid==0) {
	  ++(*progress);
	}
      }
      // END: x loop
    }
    // END: n loop
  }
  // END: m loop

  if (verbose and myid==0) std::cout << std::endl;

  for (int m=0; m<=mmax; m++) {
    for (int n=0; n<nmax; n++) {
      MPI_Allreduce(MPI_IN_PLACE, dens  [m][n].data(), numx*numy,
		    MPI_DOUBLE, MPI_SUM, MPI_COMM_WORLD);
      MPI_Allreduce(MPI_IN_PLACE, pot   [m][n].data(), numx*numy,
		    MPI_DOUBLE, MPI_SUM, MPI_COMM_WORLD);
      MPI_Allreduce(MPI_IN_PLACE, rforce[m][n].data(), numx*numy,
		    MPI_DOUBLE, MPI_SUM, MPI_COMM_WORLD);
      MPI_Allreduce(MPI_IN_PLACE, zforce[m][n].data(), numx*numy,
		    MPI_DOUBLE, MPI_SUM, MPI_COMM_WORLD);
    }
  }

  if (myid==0) {
    std::cout << "---- BiorthCyl::create_tables: done!"
	      << std::endl;
  }

  WriteH5Cache();
}


double BiorthCyl::r_to_xi(double r)
{
  if (cmapR>0) {
    if (r<0.0) {
      ostringstream msg;
      msg << "radius=" << r << " < 0! [mapped]";
      throw GenericError(msg.str(), __FILE__, __LINE__, 2040, true);
    }
    return (r/scale - 1.0)/(r/scale + 1.0);
  } else {
    if (r<0.0)  {
      ostringstream msg;
      msg << "radius=" << r << " < 0!";
      throw GenericError(msg.str(), __FILE__, __LINE__, 2040, true);
    }
    return r;
  }
}
    
double BiorthCyl::xi_to_r(double xi)
{
  if (cmapR>0) {
    if (xi<-1.0) throw GenericError("xi < -1!", __FILE__, __LINE__, 2040, true);
    if (xi>=1.0) throw GenericError("xi >= 1!", __FILE__, __LINE__, 2040, true);

    return (1.0 + xi)/(1.0 - xi) * scale;
  } else {
    return xi;
  }

}

double BiorthCyl::d_xi_to_r(double xi)
{
  if (cmapR>0) {
    if (xi<-1.0) throw GenericError("xi < -1!", __FILE__, __LINE__, 2040, true);
    if (xi>=1.0) throw GenericError("xi >= 1!", __FILE__, __LINE__, 2040, true);

    return 0.5*(1.0 - xi)*(1.0 - xi) / scale;
  } else {
    return 1.0;
  }
}

// Compute non-dimensional vertical coordinate from Z
double BiorthCyl::z_to_yi(double z)
{
  if (cmapZ==1)
    return z/(fabs(z)+std::numeric_limits<double>::min())*asinh(fabs(z/scale));
  else if (cmapZ==2)
    return z/sqrt(z*z + scale*scale);
  else
    return z;
}

// Compute Z from non-dimensional vertical coordinate
double BiorthCyl::yi_to_z(double y)
{
  if (cmapZ==1)
    return scale*sinh(y);
  else if (cmapZ==2) {
    if (y<-1.0) throw GenericError("y < -1!", __FILE__, __LINE__, 2040, true);
    if (y>=1.0) throw GenericError("y >= 1!", __FILE__, __LINE__, 2040, true);
    return y * scale/sqrt(1.0 - y*y);
  }
  else
    return y;
}

// For measure transformation in vertical coordinate
double BiorthCyl::d_yi_to_z(double y)
{
  if (cmapZ==1)
    return scale*cosh(y);
  else if (cmapZ==2) {
    if (y<-1.0) throw GenericError("y < -1!", __FILE__, __LINE__, 2040, true);
    if (y>=1.0) throw GenericError("y >= 1!", __FILE__, __LINE__, 2040, true);
    return scale*pow(1.0-y*y, -1.5);
  } else
    return 1.0;
}

// Matrix interpolation on grid for n-body
void BiorthCyl::interp(double R, double Z,
		       const std::vector<std::vector<Eigen::MatrixXd>>& mat,
		       Eigen::MatrixXd& ret, bool anti_symmetric)
{
  ret.resize(mmax+1, nmax);
  ret.setZero();

  // Off grid in radial dimension
  if (R/scale>rcylmax) return;

  double z = fabs(Z);

  // Off grid in vertical dimension
  if (z/scale>rcylmax) return;

  // Remove mid-plane discontinuity
  if (anti_symmetric and fabs(z/scale) < 1.0e-6) return;

  double X = (r_to_xi(R) - xmin)/dx;
  double Y = (z_to_yi(z) - ymin)/dy;

  int ix = (int)X;
  int iy = (int)Y;
  
  if (ix < 0) {
    ix = 0;
    X  = 0.0;
  }
  if (iy < 0) {
    iy = 0;
    Y  = 0.0;
  }
  
  if (ix >= numx-1) {
    ix = numx-2;
    X  = numx-1;
  }
  if (iy >= numy-1) {
    iy = numy-2;
    Y  = numy-1;
  }
  
  double delx0 = (double)ix + 1.0 - X;
  double dely0 = (double)iy + 1.0 - Y;
  double delx1 = X - (double)ix;
  double dely1 = Y - (double)iy;
  
  double c00 = delx0*dely0;
  double c10 = delx1*dely0;
  double c01 = delx0*dely1;
  double c11 = delx1*dely1;
  
  for (int m=0; m<=mmax; m++) {
    for (int n=0; n<nmax; n++) {
      ret(m, n) =
	mat[m][n](ix  , iy  ) * c00 +
	mat[m][n](ix+1, iy  ) * c10 +
	mat[m][n](ix  , iy+1) * c01 +
	mat[m][n](ix+1, iy+1) * c11 ;
    }
  }

  if (Z<0.0 and anti_symmetric) ret *= -1.0;
}


double BiorthCyl::interp(int m, int n, double R, double Z,
			 const std::vector<std::vector<Eigen::MatrixXd>>& mat,
			 bool anti_symmetric)
{
  double ret = 0.0;

  // Off grid in radial dimension
  if (R/scale>rcylmax or n>=nmax or m>mmax) return ret;

  double z = fabs(Z);

  // Off grid in vertical dimension
  if (z/scale>rcylmax) return ret;

  // Remove mid-plane discontinuity
  if (anti_symmetric and z/scale < 1.0e-6) return ret;

  double X = (r_to_xi(R) - xmin)/dx;
  double Y = (z_to_yi(z) - ymin)/dy;

  int ix = (int)X;
  int iy = (int)Y;
  
  if (ix < 0) {
    ix = 0;
    X  = 0.0;
  }
  if (iy < 0) {
    iy = 0;
    Y  = 0.0;
  }
  
  if (ix >= numx-1) {
    ix = numx-2;
    X  = numx-1;
  }
  if (iy >= numy-1) {
    iy = numy-2;
    Y  = numy-1;
  }
  
  double delx0 = (double)ix + 1.0 - X;
  double dely0 = (double)iy + 1.0 - Y;
  double delx1 = X - (double)ix;
  double dely1 = Y - (double)iy;
  
  double c00 = delx0*dely0;
  double c10 = delx1*dely0;
  double c01 = delx0*dely1;
  double c11 = delx1*dely1;
  
  ret = 
    mat[m][n](ix  , iy  ) * c00 +
    mat[m][n](ix+1, iy  ) * c10 +
    mat[m][n](ix  , iy+1) * c01 +
    mat[m][n](ix+1, iy+1) * c11 ;

  if (Z<0.0 and anti_symmetric) ret *= -1.0;

  return ret;
}

void BiorthCyl::WriteH5Params(HighFive::File& file)
{
  file.createAttribute<int>         ("mmax",     HighFive::DataSpace::From(mmax)).write(mmax);
  file.createAttribute<int>         ("nmax",     HighFive::DataSpace::From(nmax)).write(nmax);
  file.createAttribute<int>         ("numr",     HighFive::DataSpace::From(numr)).write(numr);
  file.createAttribute<int>         ("nmaxfid",  HighFive::DataSpace::From(nmax)).write(nmaxfid);
  file.createAttribute<int>         ("numx",     HighFive::DataSpace::From(numx)).write(numx);
  file.createAttribute<int>         ("numy",     HighFive::DataSpace::From(numy)).write(numy);
  file.createAttribute<double>      ("rcylmin",  HighFive::DataSpace::From(rcylmin)).write(rcylmin);
  file.createAttribute<double>      ("rcylmax",  HighFive::DataSpace::From(rcylmax)).write(rcylmax);
  file.createAttribute<double>      ("scale",    HighFive::DataSpace::From(scale)).write(scale);
  file.createAttribute<int>         ("cmapR",    HighFive::DataSpace::From(cmapR)).write(cmapR);
  file.createAttribute<int>         ("cmapZ",    HighFive::DataSpace::From(cmapZ)).write(cmapZ);
}


void BiorthCyl::WriteH5Arrays(HighFive::Group& harmonic)
{
  for (int m=0; m<=mmax; m++) {
    std::ostringstream sout;
    sout << m;
    auto order = harmonic.createGroup(sout.str());
      
    for (int n=0; n<nmax; n++) {
      std::ostringstream sout;
      sout << n;
      auto arrays = order.createGroup(sout.str());

      HighFive::DataSet ds1 = arrays.createDataSet("density",   dens  [m][n]);
      HighFive::DataSet ds2 = arrays.createDataSet("potential", pot   [m][n]);
      HighFive::DataSet ds3 = arrays.createDataSet("rforce",    rforce[m][n]);
      HighFive::DataSet ds4 = arrays.createDataSet("zforce",    zforce[m][n]);
    }
  }
}

void BiorthCyl::ReadH5Arrays(HighFive::Group& harmonic)
{
  for (int m=0; m<=mmax; m++) {
    std::ostringstream sout;
    sout << m;
    auto order = harmonic.getGroup(sout.str());
      
    for (int n=0; n<nmax; n++) {
      std::ostringstream sout;
      sout << n;
      auto arrays = order.getGroup(sout.str());

      arrays.getDataSet("density")  .read(dens  [m][n]);
      arrays.getDataSet("potential").read(pot   [m][n]);
      arrays.getDataSet("rforce")   .read(rforce[m][n]);
      arrays.getDataSet("zforce")   .read(zforce[m][n]);
    }
  }

}

void BiorthCyl::WriteH5Cache()
{
  if (myid) return;

  try {
    // Create a new hdf5 file or overwrite an existing file
    //
    HighFive::File file(cachename, HighFive::File::Overwrite);
    
    // We write the basis geometry
    //
    file.createAttribute<std::string>("geometry", HighFive::DataSpace::From(geometry)).write(geometry);
      
    // We write the force ID
    //
    file.createAttribute<std::string>("forceID", HighFive::DataSpace::From(forceID)).write(forceID);
      
    // Stash the basis configuration (this is not yet implemented in EXP)
    //
    std::ostringstream sout; sout << conf;
    std::string config(sout.str());
    file.createAttribute<std::string>("config", HighFive::DataSpace::From(config)).write(config);
      
    // Write the specific parameters
    //
    WriteH5Params(file);
      
    // Harmonic order (for h5dump readability)
    //
    auto harmonic = file.createGroup("Harmonic");

    WriteH5Arrays(harmonic);

  } catch (HighFive::Exception& err) {
    std::cerr << err.what() << std::endl;
  }
    
  std::cout << "---- BiorthCyl::WriteH5Cache: "
	    << "wrote <" << cachename << ">" << std::endl;
}
  
bool BiorthCyl::ReadH5Cache()
{
  // First attempt to read the file
  //
  try {
    // Silence the HDF5 error stack
    //
    HighFive::SilenceHDF5 quiet;
    
    // Try opening the file as HDF5
    //
    HighFive::File h5file(cachename, HighFive::File::ReadOnly);
    
    // Try checking the rest of the parameters before reading arrays
    //
    auto checkInt = [&h5file](int value, std::string name)
    {
      int v; HighFive::Attribute vv = h5file.getAttribute(name); vv.read(v);
      if (value == v) return true;
      else {
	if (myid==0)  std::cout << "---- BiorthCyl::ReadH5Cache(): wanted "
				<< name << "=" << value
				<< " but found "
				<< name << "=" << v << std::endl;
	return false;
      }
    };

    auto checkDbl = [&h5file](double value, std::string name)
    {
      double v; HighFive::Attribute vv = h5file.getAttribute(name); vv.read(v);
      if (fabs(value - v) < 1.0e-16) return true;
      else {
	if (myid==0) std::cout << "---- BiortyCyl::ReadH5Cache: wanted "
			       << name << "=" << value
			       << " but found "
			       << name << "=" << v << std::endl;
	return false;
      }
    };

    auto checkStr = [&h5file](std::string value, std::string name)
    {
      std::string v; HighFive::Attribute vv = h5file.getAttribute(name); vv.read(v);
      if (value.compare(v)==0) return true;
      else {
	if (myid==0) std::cout << "ReadH5Cache(): wanted "
			       << name << "=" << value
			       << " but found "
			       << name << "=" << v << std::endl;
	return false;
      }
    };


    // ID check
    //
    if (not checkStr(geometry, "geometry"))  return false;
    if (not checkStr(forceID,  "forceID"))   return false;

    // Parameter check
    //
    if (not checkInt(mmax,     "mmax"))      return false;
    if (not checkInt(nmax,     "nmax"))      return false;
    if (not checkInt(numr,     "numr"))      return false;
    if (not checkInt(nmaxfid,  "nmaxfid"))   return false;
    if (not checkInt(numx,     "numx"))      return false;
    if (not checkInt(numy,     "numy"))      return false;
    if (not checkDbl(rcylmin,  "rcylmin"))   return false;
    if (not checkDbl(rcylmax,  "rcylmax"))   return false;
    if (not checkDbl(scale,    "scale"))     return false;
    if (not checkInt(cmapR,    "cmapR"))     return false;
    if (not checkInt(cmapZ,    "cmapZ"))     return false;

    // Open the harmonic group
    //
    auto harmonic = h5file.getGroup("Harmonic");

    ReadH5Arrays(harmonic);

    if (myid==0)
      std::cerr << "---- BiorthCyl::ReadH5Cache: "
		<< "read <" << cachename << ">" << std::endl;

    // Create the basis instance
    //
    emp = EmpCyl2d(mmax, nmaxfid, nmax, knots, numr,
		   rcylmin*scale, rcylmax*scale, 
		   scale, cmapR, logr, diskconf, biorth);
    
    if (conf["basis"]) emp.basisTest(true);

    return true;
    
  } catch (HighFive::Exception& err) {
    if (myid==0)
      std::cerr << "---- BiorthCyl::ReadH5Cache: "
		<< "error opening HDF5 basis cache <"
		<< cachename << ">" << std::endl;
  }

  return false;
}

// z coordinate is ignored here; assumed to be zero!
void BiorthCyl::get_pot(Eigen::MatrixXd& Vc, Eigen::MatrixXd& Vs,
			double r, double z)
{
  Vc.resize(std::max<int>(1, mmax)+1, nmax);
  Vs.resize(std::max<int>(1, mmax)+1, nmax);

  double X = (r_to_xi(r) - xmin)/dx;

  int ix = (int)X;

  if (ix < 0) {
    ix = 0;
    X = 0.0;
  }
  
  if (ix >= numx-1) {
    ix = numx-2;
    X  = numx-1;
  }

  double delx0 = (double)ix + 1.0 - X;
  double delx1 = X - (double)ix;

  double c0 = delx0;
  double c1 = delx1;

  for (int mm=0; mm<=std::min<int>(mlim, mmax); mm++) {
    
    // Suppress odd M terms?
    if (EVEN_M && (mm/2)*2 != mm) continue;

    for (int n=0; n<nmax; n++) {

      Vc(mm, n) = pot[mm][n](ix, 0)*c0 +  pot[mm][n](ix+1, 0)*c1;
      Vs(mm, n) = Vc(mm, n);
    }
  }
}

YAML::Node BiorthCyl::getHeader(const std::string& cachefile)
{
  std::ifstream in(cachefile);
  if (!in) {
    std::ostringstream sout;
    sout << "BiorthCyl::getHeader: could not open cache file <"
	 << cachefile << ">";
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

    node["mmax"]     = getInt("mmax");
    node["nmax"]     = getInt("nmax");
    node["numr"]     = getInt("numr");
    node["nmaxfid"]  = getInt("nmaxfid");
    node["numx"]     = getInt("numx");
    node["numy"]     = getInt("numy");
    node["rcylmin"]  = getDbl("rcylmin");
    node["rcylmax"]  = getDbl("rcylmax");
    node["scale"]    = getDbl("scale");
    node["cmapR"]    = getInt("cmapR");
    node["cmapZ"]    = getInt("cmapZ");
  }
  catch (YAML::Exception& error) {
    std::ostringstream sout;
    sout << "BiorthCyl::getHeader: invalid cache file <" << cachefile << ">. ";
    sout << "YAML error in getHeader: " << error.what() << std::endl;
    throw GenericError(sout.str(), __FILE__, __LINE__, 1038, false);
  }

  return node;
}


std::map<std::string, std::string>
BiorthCyl::cacheInfo(const std::string& cachefile, bool verbose)
{
  std::map<std::string, std::string> ret;
  auto node = getHeader(cachefile);

  if (verbose)
    std::cout << std::string(60, '-') << std::endl
	      << "Cache parameters for BiorthCyl: " << cachefile << std::endl
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


std::vector<Eigen::MatrixXd> BiorthCyl::orthoCheck()
{
  if (not emp()) {
    emp = EmpCyl2d(mmax, nmaxfid, nmax, knots, numr,
		   rcylmin*scale, rcylmax*scale, scale,
		   cmapR, logr, diskconf, biorth);
  }

  return emp.orthoCheck();
}
