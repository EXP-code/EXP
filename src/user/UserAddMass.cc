#include <math.h>
#include <sstream>
#include <numeric>

#include "expand.h"
#include <localmpi.h>

#include <ACG.h>
#include <Uniform.h>
#include <Normal.h>

#include <UserAddMass.H>

// For parsing parameters: Algorithm type
//
std::map<std::string, UserAddMass::Algorithm> UserAddMass::AlgMap
{ {"Disk",     UserAddMass::Algorithm::Disk},
  {"Halo",     UserAddMass::Algorithm::Halo} };

// For printing parameters: Algorithm type
//
std::map<UserAddMass::Algorithm, std::string> UserAddMass::AlgMapInv
{ {UserAddMass::Algorithm::Disk, "Disk"},
  {UserAddMass::Algorithm::Halo, "Halo"} };


// Compute cross product for two 3 std::arrays
//
std::array<double, 3> xprod(const std::array<double, 3>& a,
			    const std::array<double, 3>& b)
{
  std::array<double, 3> ret =
    {
     a[1]*b[2] - a[2]*b[1],
     a[2]*b[0] - a[0]*b[2],
     a[0]*b[1] - a[1]*b[0]
    };
  return ret;
}

// Compute cross product for two 3 C-style arrays
//
std::array<double, 3> xprod(const double* a, const double* b)
{
  std::array<double, 3> ret =
    {
     a[1]*b[2] - a[2]*b[1],
     a[2]*b[0] - a[0]*b[2],
     a[0]*b[1] - a[1]*b[0]
    };
  return ret;
}

// Compute the norm of a std::array<double, 3>
//
double xnorm(const std::array<double, 3>& a)
{
  return sqrt(a[0]*a[0] + a[1]*a[1] + a[2]*a[2]);
}

UserAddMass::UserAddMass(const YAML::Node &conf) : ExternalForce(conf)
{
  id = "AddMass";

  comp_name = "";               // Default component for com
  rmin      = 1.0e-4;           // Default inner radius
  rmax      = 1.0;              // Default outer radius
  numr      = 100;              // Default number of radial bins
  number    = 10;               // Default number of particles added
  mass      = 1.0e-10;          // Default mass of new particles
  logr      = true;             // Logaritmic binning
  seed      = 11;               // Random number seed
  dtime     = 0.012;            // Time interval between particle additions
  tnext     = dtime;            // Time to begin adding particles
  debug     = -1.0;		// Time interval between diagnostic output
  scale     =  1.0;             // Radius scale of mass adding
  vdisp     = 0.0;              // Velocity dispersion value
  interp    = true;		// Interpolate by default
  cforce    = false;		// Explicit force evalution for orbit IC
  planar    = false;		// Force disk to be in x-y plane
  accel     = true;             // Assign acceleration from basis
  alg       = AlgMap["Disk"];	// Assign the disk algorithm by default
  

  initialize();			// Read and set parameters

  dnext  = tnext;		// Time for next debug output
  tstart = tnext;

  if (comp_list.size()>0) {
				// Components for force evaluation
    for (auto name : comp_list) {
      bool found = false; 
      for (auto c : comp->components) { // Look for each specified component
	if ( !name.compare(c->name) ) {
	  comp_vec.push_back(c);
	  found = true;
	  break;
	}
      }

      if (!found) {
	cerr << "Process " << myid << ": can't find desired list component <"
	     << name << ">" << endl;
	MPI_Abort(MPI_COMM_WORLD, 34);
      }
    }
    cforce = true;
  }


  if (comp_name.size()>0) {
				// Look for the fiducial component
    bool found = false;
    for (auto c : comp->components) {
      if ( !comp_name.compare(c->name) ) {
	c0 = c;
	found = true;
      break;
      }
    }

    if (!found) {
      cerr << "Process " << myid << ": can't find desired component <"
	   << comp_name << ">" << endl;
      MPI_Abort(MPI_COMM_WORLD, 35);
    }

  } else {
    if (myid==0) {
      std:: cerr << "UserAddMass: desired component name must be specified"
	<< std::endl;
	   MPI_Abort(MPI_COMM_WORLD, 36);
    }
  }

  // Reflection to find components that provide forces from a basis
  // (which field forces at field points)
  for (auto cc : comp->components) {
    Basis* b = dynamic_cast<Basis*>(cc->force);
    if (b != NULL) {
      comp_acc.push_back(cc);
      std::cout << "Found <" << cc->name << ", " << cc->id << ">" << std::endl;
    }
  }

  // Make bins parameters
				// Sanity
  if (rmin == 0.0 and logr) logr = false;

  lrmin = rmin;
  lrmax = rmax;

  if (logr) {			// Log binning?
    lrmin = log(rmin);
    lrmax = log(rmax);
  }

  dr = (lrmax - lrmin)/(numr - 1);

  mas_bins_current.resize(nthrds);
  mas_bins_adding.resize(numr);
  vl_bins.resize(nthrds);
  v2_bins.resize(nthrds);
  L3_bins.resize(nthrds);

  // Making binning for each algorithm
  //
  for (auto & v : mas_bins_current) v.resize(numr);

  switch (alg) {
  case Algorithm::Halo:
    for (auto & v : vl_bins) v.resize(numr);
    for (auto & v : v2_bins) v.resize(numr);
    break;
  case Algorithm::Disk:
  default:
    for (auto & v : L3_bins) v.resize(numr);
    break;
  };

  // Random number generation
  //
  gen   = boost::make_shared<ACG>    (seed+myid);
  urand = boost::make_shared<Uniform>(0.0, 1.0, gen.get());
  nrand = boost::make_shared<Normal> (0.0, 1.0, gen.get());

  userinfo();
}

UserAddMass::~UserAddMass()
{
}

void UserAddMass::userinfo()
{
  if (myid) return;		// Return if node master node

  print_divider();

  cout << "** User routine ADD MASS initialized "
       << "using component <" << comp_name << "> with"
       << " rmin="      << rmin
       << " rmax="      << rmax
       << " numr="      << numr
       << " logr="      << std::boolalpha << logr
       << " number="    << number
       << " mass="      << mass
       << " scale="     << scale
       << " algorithm=" << AlgMapInv[alg]
       << " planar="    << std::boolalpha << planar
       << std::endl;

  if (comp_list.size()>0) {
    cout << " force evaluation from components [";
    for (auto s : comp_list) cout << " <" << s << ">";
    cout << " ]";
  }
  cout << std::endl;

  print_divider();
}

void UserAddMass::initialize()
{
  try {
    if (conf["complist"]) {	// Check for optional force components
      comp_list = conf["complist"].as<std::vector<std::string>>();
    }
    if (conf["compname"])         comp_name   = conf["compname"].as<std::string>();	      
    if (conf["rmin"])       	  rmin        = conf["rmin"].as<double>();		      
    if (conf["rmax"])       	  rmax        = conf["rmax"].as<double>();		      
    if (conf["mass"])       	  mass        = conf["mass"].as<double>();		      
    if (conf["numr"])       	  numr        = conf["numr"].as<unsigned>();		      
    if (conf["interp"])     	  interp      = conf["interp"].as<bool>();		      
    if (conf["planar"])     	  planar      = conf["planar"].as<bool>();		      
    if (conf["logr"])       	  logr        = conf["logr"].as<bool>();		      
    if (conf["seed"])       	  seed        = conf["seed"].as<long int>();		      
    if (conf["number"])     	  number      = conf["number"].as<unsigned>();	      
    if (conf["tstart"])     	  tnext       = conf["tstart"].as<double>();		      
    if (conf["scale"])      	  scale       = conf["scale"].as<double>();		      
    if (conf["dtime"])      	  dtime       = conf["dtime"].as<double>();		      
    if (conf["debug"])      	  debug       = conf["debug"].as<double>();		      
    if (conf["vdispersion"])      vdisp       = conf["vdispersion"].as<double>();		      
    if (conf["algorithm"])  	  alg         = AlgMap[conf["algorithm"].as<std::string>()];
  }
  catch (YAML::Exception & error) {
    if (myid==0) std::cout << "Error parsing parameters in UserAddMass: "
			   << error.what() << std::endl
			   << std::string(60, '-') << std::endl
			   << "Config node"        << std::endl
			   << std::string(60, '-') << std::endl
			   << conf                 << std::endl
			   << std::string(60, '-') << std::endl;
    MPI_Finalize();
    exit(-1);
  }

}

void UserAddMass::clear_bins()
{
  // Mass bins
  //

  for (auto & v : mas_bins_current) {
    std::fill(v.begin(), v.end(), 0.0);
  }
  

  std::array<double, 3> zero = {0.0, 0.0, 0.0};

  switch (alg) {
  case Algorithm::Halo:
    // Mean and squared velocity bins
    for (auto & v : vl_bins) std::fill(v.begin(), v.end(), zero);
    for (auto & v : v2_bins) std::fill(v.begin(), v.end(), zero);
    break;
  case Algorithm::Disk:
  default:
    // Angular momentum bins
    for (auto & v : L3_bins) std::fill(v.begin(), v.end(), zero);
    break;
  }
}

void UserAddMass::determine_acceleration_and_potential(void)
{
  if (cC != c0) return;		// Check that this component is the target

				// Only compute for top level
  if (multistep && mlevel>0) return;

				// Have we completed the interval?
  if (tnow < tnext) return;

  tnext += dtime;		// Increment the target time

  clear_bins();			// Clean for tabulation

  exp_thread_fork(false);	// Tabulate one-d radial bins

  // Thread reduce the bins
  //
  for (int n=1; n<nthrds; n++) {

    for (unsigned i=0; i<numr; i++) {

      mas_bins_current[0][i] += mas_bins_current[n][i];

      switch (alg) {
      case Algorithm::Halo:
	for (int k=0; k<3; k++) vl_bins[0][i][k] += vl_bins[n][i][k];
	for (int k=0; k<3; k++) v2_bins[0][i][k] += v2_bins[n][i][k];
	break;
      case Algorithm::Disk:
      default:
	for (int k=0; k<3; k++) L3_bins[0][i][k] += L3_bins[n][i][k];
      }
      break;
    }
  }

  // Reduce the bins across all nodes
  //
  std::vector<double> mas(numr), pack(numr*3);

  MPI_Allreduce(mas_bins_current[0].data(), mas.data(), numr, MPI_DOUBLE, MPI_SUM,
		MPI_COMM_WORLD);
  if ( fabs( tnow - tstart)<3E-4 ){
    for (unsigned i=0; i<numr; i++) mas_bins_adding[i] = mas[i];
  }
  switch (alg) {
  case Algorithm::Halo:
    for (unsigned i=0; i<numr; i++) {
      for (int k=0; k<3; k++) pack[i*3+k] = vl_bins[0][i][k];
    }
    MPI_Allreduce(MPI_IN_PLACE, pack.data(), numr*3, MPI_DOUBLE, MPI_SUM,
		  MPI_COMM_WORLD);
    for (unsigned i=0; i<numr; i++) {
      for (int k=0; k<3; k++) vl_bins[0][i][k] = pack[i*3+k];
      for (int k=0; k<3; k++) pack[i*3+k] = v2_bins[0][i][k];
    }
    MPI_Allreduce(MPI_IN_PLACE, pack.data(), numr*3, MPI_DOUBLE, MPI_SUM,
		  MPI_COMM_WORLD);
    for (unsigned i=0; i<numr; i++) {
      for (int k=0; k<3; k++) v2_bins[0][i][k] = pack[i*3+k];
    }

    break;
  case Algorithm::Disk:
  default:
    for (unsigned i=0; i<numr; i++) {
      for (int k=0; k<3; k++) pack[i*3+k] = L3_bins[0][i][k];
    }
    MPI_Allreduce(MPI_IN_PLACE, pack.data(), numr*3, MPI_DOUBLE, MPI_SUM,
		  MPI_COMM_WORLD);
    for (unsigned i=0; i<numr; i++) {
      for (int k=0; k<3; k++) L3_bins[0][i][k] = pack[i*3+k];
    }
    break;
  }

  // Compute cumulative mass
  //
  std::vector<double> wght(mas_bins_adding);
  for (unsigned i=1; i<numr; i++) wght[i] += wght[i-1];


  // BEG: diagnostic output
  //
  //  +------ Only root node writes diag output
  //  |
  //  |           +------ Write interval >0 flag
  //  |           |
  //  |           |               +------ Passed target time
  //  |           |               |
  //  v           v               v
  if (myid==0 and debug > 0.0 and tnow >= dnext) {
    // Create file name
    //
    static unsigned counter = 0;
    std::ostringstream sout;
    sout << comp_name << ".addmass." << runtag << "."
	 << std::setw(5) << std::setfill('0') << std::right << counter++;

    std::ofstream out(sout.str());

    if (out) {
      // File header
      //
      out << "# Time=" << tnow << std::endl;
      const std::vector<std::string> labels =
	{
	 "Radius", "Mass(current)", "Cum mass(adding)", "Velocity"
	};
      for (size_t i=0; i<labels.size(); i++) {
	if (i) out << std::setw(18) << std::right << "|";
	else   out << "# " << std::setw(16) << "|";
      }
      out << std::endl;
      for (size_t i=0; i<labels.size(); i++) {
	std::ostringstream lout; lout << labels[i] << " |";
	if (i) out << std::setw(18) << std::right << lout.str();
	else   out << "# " << std::setw(16) << lout.str();
      }
      out << std::endl;
      for (size_t i=0; i<labels.size(); i++) {
	std::ostringstream lout; lout << i+1 << " |";
	if (i) out << std::setw(18) << std::right << lout.str();
	else   out << "# " << std::setw(16) << lout.str();
      }
      out << std::endl << std::setfill('-');
      for (size_t i=0; i<labels.size(); i++) {
	std::ostringstream lout; lout << " |";
	if (i) out << std::setw(18) << std::right << lout.str();
	else   out << "# " << std::setw(16) << lout.str();
      }
      out << std::endl << std::setfill(' ');
      //
      // END: file header

      for (unsigned i=0; i<numr; i++) {
	double rr = lrmin + dr*i;
	if (logr) rr = exp(rr);

	out << std::setw(18) << rr	 // radius
	    << std::setw(18) << mas[i]	 // mass(current)
	    << std::setw(18) << wght[i]; // cum mas(adding)
    
	double vv = 0.0;	// velocity
	if (mas[i]>0.0) {
	  switch (alg) {
	  case Algorithm::Halo:
	    {
	      for (int k=0; k<3; k++) {
		double v1 = vl_bins[0][i][k]/mas[i];
		double v2 = v2_bins[0][i][k]/mas[i];
		vv += v2 - v1*v1;
	      }
	      vv = sqrt(fabs(vv));
	    }
	    break;
	  case Algorithm::Disk:
	  default:
	    if (rr>0.0) {	// Prevent 0/0
	      for (int k=0; k<3; k++) vv += L3_bins[0][i][k]*L3_bins[0][i][k];
	      vv = sqrt(vv)/mas[i]/rr;
	    }
	    break;
	  }
	}
      }

    } else {
      std::cerr << "UserAddMass: could not open diagnostic file <"
		<< sout.str() << ">" << std::endl;
    }
    
    // Increment the target time
    //
    dnext += debug;
  }
  // END: diagnostic output


  // Compute total angular momentum
  //
  std::array<std::array<double, 3>, 3> ee;
  for (int k=0; k<3; k++) ee[0][k] = 0.0;



  if (alg == Algorithm::Disk) {

    for (unsigned i=0; i<numr; i++) {
      for (int k=0; k<3; k++) ee[0][k] += L3_bins[0][i][k];
    }

    if (planar) {		// FOR TESTING

      if (ee[0][2]>0.0) {	// Right-handed spin
	      ee[0] = { 0.0,  0.0,  1.0}; // +z
	      ee[1] = { 1.0,  0.0,  0.0}; // +x
	      ee[2] = { 0.0,  1.0,  0.0}; // +y
      }  else {			// Left-handed spin
	      ee[0] = { 0.0,  0.0, -1.0}; // -z
	      ee[1] = { 1.0,  0.0,  0.0}; // +x
	      ee[2] = { 0.0, -1.0,  0.0}; // -y
      }

    } else {			// Use angular momentum axes

    // First orthogonal vector
    //
    double norm = xnorm(ee[0]);
    for (auto & v : ee[0]) v /= norm;

    // Second orthogonal vector
    //
    ee[1] = {ee[0][1], -ee[0][0], 0.0};
    norm = xnorm(ee[1]);
    for (auto & v : ee[1]) v /= norm;

    // Third orthogonal vector
    //
    ee[2] = xprod(ee[0], ee[1]);
    norm = xnorm(ee[2]);
    for (auto & v : ee[2]) v /= norm;

    }
  }

  if (fabs(tnow - tstart) < 3E-4) return;

  // Finally, we are ready to generate the new particles
  //
  int numgen = 1;
				// Fewer particles than processes
  if (static_cast<int>(number) < numprocs) {	
    if (myid >= static_cast<int>(number)) numgen = 0;
  } else {			// Otherwise . . .
    numgen = number/numprocs;
    if (myid==numprocs-1) numgen = number - numgen*(numprocs-1);
  }

  for (int n=0; n<numgen; n++) {
    bool notfound = true;
    while (notfound) {
      // Get radius
      //
      auto it = std::lower_bound(wght.begin(), wght.end(), wght.back()*(*urand)());
      unsigned indx = numr-1;
      if (it != wght.end()) indx = std::distance(wght.begin(), it);
      if (indx>0) indx--;

      // Decrease index to find non-zero bins(s)
      //

      std::function<bool(std::vector<double>&, int)> mfind;
      if (interp) {		// Interpolation predicate
        mfind = [](auto &mas, auto indx)
          { return (mas[indx]==0.0 or mas[indx+1]==0.0); };

        indx = std::min<unsigned>(indx, numr-2);
        } else {			// Nearest bin predicate
          mfind = [](auto &mas, auto indx)
		      { return mas[indx]==0.0; };
          indx = std::min<unsigned>(indx, numr-1);
        }
      while (mfind(mas, indx) and indx>0) indx--;
      // Skip this one
      //
      if (mfind(mas, indx)) break; 
      // Found good bin(s)
      //
      notfound = false;      
      // Radius selection
      //
      double lrr = lrmin + dr*((*urand)() + indx);
      double rr  = lrr, a=1.0, b=0.0;
      if (logr) rr = exp(lrr);	// Using log scaling for bins      
      if (interp) {
      	a = (lrmin + dr*(indx+1) - lrr)/dr;
      	b = (lrr - lrmin - dr*(indx+0))/dr;
      }
      
      Particle *P =  c0->GetNewPart();
      P->mass  = mass;
      P->level = multistep;

      switch (alg) {
      case Algorithm::Halo:
	{
	  // Positions
	  //
	  double cosx = 2.0*(*urand)() - 1.0;
	  double sinx = sqrt(fabs(1.0 - cosx*cosx));
	  double phi  = 2.0*M_PI*(*urand)();
	  double cosp = cos(phi), sinp = sin(phi);
	  
	  std::array<std::array<double, 3>, 3> ev;
	  
	  ev[0][0] =  sinx * cosp;
	  ev[0][1] =  sinx * sinp;
	  ev[0][2] =  cosx ;
	  
	  ev[1][0] =  cosx * cosp;
	  ev[1][1] =  cosx * sinp;
	  ev[1][2] = -sinx ;
	  
	  ev[2][0] = -sinp;
	  ev[2][1] =  cosp;
	  ev[2][2] =  0.0;
	  
	  // Positions
	  //
	  for (int k=0; k<3; k++) P->pos[k] = rr*ev[0][k];
	  
	  // Velocities
	  //
	  double sig = 0.0;
	  for (int k=0; k<3; k++) {
	    if (interp) {
	      double v1a = vl_bins[0][indx  ][k]/mas[indx  ];
	      double v2a = v2_bins[0][indx  ][k]/mas[indx  ];
	      double v1b = vl_bins[0][indx+1][k]/mas[indx+1];
	      double v2b = v2_bins[0][indx+1][k]/mas[indx+1];
	      sig += a*(v2a - v1a*v1a) + b*(v2b - v1b*v1b);
	    }
	    else {
	      double v1 = vl_bins[0][indx][k]/mas[indx];
	      double v2 = v2_bins[0][indx][k]/mas[indx];
	      sig += v2 - v1*v1;
	    } 
	  }
	  sig = sqrt(fabs(sig));
	  
	  std::array<double, 3> vv = {(*nrand)(), (*nrand)(), (*nrand)()};
	  double vnorm = xnorm(vv);
	  for (auto & v : vv) v /= vnorm;
	  
	  for (int k=0; k<3; k++) P->vel[k] = sig * vv[k];
	}
	
	break;
      case Algorithm::Disk:
      default:
	{
	  // Positions
	  //
	  double phi = 2.0*M_PI*(*urand)();
	  double cosp = cos(phi), sinp = sin(phi);
	  for (int k=0; k<3; k++) P->pos[k] = rr*(ee[1][k]*cosp + ee[2][k]*sinp);
	  //                                      ^               ^
	  // These are perpendicular to the       |               |
	  // mean angular momemtnum vector -------+---------------+
	  
	  if (cforce) {		// Use force evaluation to get
				// tangential veocity
	    double tdens0, dpotl0, tdens, tpotl, tpotr, tpotz, tpotp;
	    double z = P->pos[2], dPdR = 0.0;

	    for (auto c : comp_vec) {

	      reinterpret_cast<Basis*>(c->force)->determine_fields_at_point_cyl
		(rr, z, phi, &tdens0, &dpotl0, &tdens, &tpotl,
		 &tpotr, &tpotz, &tpotp);
	      
	      dPdR += tpotr;
	    }

	    double vt = sqrt(rr*fabs(dPdR));
	    if (L3_bins[0][indx][2] <= 0.0) vt *= -1.0;
	    
	    // Velocities
	    //
	    for (int k=0; k<3; k++) P->vel[k] = vt*(ee[2][k]*cosp - ee[1][k]*sinp);
	    //                                      ^               ^
	    // These are perpendicular to the       |               |
	    // mean angular momemtnum vector -------+---------------+

	    // Test for NaN
	    //
	    for (int k=0; k<3; k++) {
	      if (std::isnan(P->pos[k])) {
		std::cout << "UserAddMass [" << k << "] NaN position" << std::endl;
	      }
	      if (std::isnan(P->vel[k])) {
		std::cout << "UserAddMass [" << k << "] NaN velocity" << std::endl;
	      }
	    }

	    // TEST
	    double perp = 0.0;
	    for (int k=0; k<3; k++) perp += P->pos[k]*P->vel[k];
	    if (fabs(perp)>1.0e-8) {
	      std::cout << "Perp=" << perp << std::endl;
	    }
	    
	  } else {
	    // Use angular momentum to get tangential velocity.
	    // E.g. L X r/r^2 = (r X v) X r/r^2 = v - r/|r| * (v.r/|r|) = v_t
	    //
	    vec3 rv;
	    if (interp) {
  	      for (int k=0; k<3; k++)
		rv[k] = a*L3_bins[0][indx][k]/mas[indx+0] + b*L3_bins[0][indx+1][k]/mas[indx+1];
	    } else {
	      rv = L3_bins[0][indx];
	      for (auto & v : rv) v /= mas[indx];
	    }

	    std::array<double, 3> xyz = {P->pos[0], P->pos[1], P->pos[2]};
	    auto v3 = xprod(rv, xyz);
	    for (auto & v : v3) v /= rr*rr;
	    
	    // v3 is now the mean tangential velocity perpendicular to
	    // radius
	    
	    // Velocities
	    //
	    for (int k=0; k<3; k++) P->vel[k] = v3[k] *( (*nrand)() * 1/1.732 * vdisp + 1.0);
	  }
	}

	break;
      }
      
      // Add acceleration
      //
      if (accel) {
	for (auto cc : comp_acc) {
	  double tt, tX, tY, tZ;
	  dynamic_cast<Basis*>(cc->force)->
	    determine_fields_at_point(P->pos[0], P->pos[1], P->pos[2],
				      &tt, &tt, &tt, &tt, &tX, &tY, &tZ);
	  P->acc[0] = -tX;
	  P->acc[0] = -tY;
	  P->acc[0] = -tZ;
	}
      }
    }

  }

}


void * UserAddMass::determine_acceleration_and_potential_thread(void * arg) 
{
  int nbodies = c0->Number();
  int id = *((int*)arg);
  int nbeg = nbodies*id/nthrds;
  int nend = nbodies*(id+1)/nthrds;

  thread_timing_beg(id);

  PartMapItr it = cC->Particles().begin();

  std::array<double, 3> pos, vel, angm;

  for (int q=0   ; q<nbeg; q++) it++;

  for (int q=nbeg; q<nend; q++, it++) {

    auto P = it->second;

    for (int k=0; k<3; k++) {
      pos[k] = P->pos[k];
      vel[k] = P->vel[k];
    }
	
    // Convert to center reference frame
    //
    cC->ConvertPos(pos.data());
    cC->ConvertVel(vel.data());

    // Compute radius
    //
    double r = xnorm(pos);

    // Compute bin number
    //
    int indx = -1;		// Off grid, by default
    double s;
    if ( fabs( tnow - tstart)<3E-4 )
      s=scale;
    else
      s=1;
    if (logr) {			// Log binning:
      if (r>0.0) {		// ignore if r==0.0 in
	r = log(r*s);
	if (r >= lrmin) indx = floor( (r - lrmin)/dr );
      }
    } else {			// Linear binning
      if (r >= lrmin) indx = floor( (r*s - lrmin)/dr );
    }

    // Add value to bin
    //
    if (indx >= 0 and indx < static_cast<int>(numr)) {
      
      mas_bins_current[id][indx] += P->mass;

      switch (alg) {
      case Algorithm::Halo:
	for (int k=0; k<3; k++)	// Mean velocity
	  vl_bins[id][indx][k] += P->mass * vel[k];
	for (int k=0; k<3; k++)	// Square velocity
	  v2_bins[id][indx][k] += P->mass * vel[k]*vel[k];
	break;
      case Algorithm::Disk:	// Angular momentum
      default:
	angm = xprod(pos, vel);
	L3_bins[id][indx][0] += P->mass * angm[0];
	L3_bins[id][indx][1] += P->mass * angm[1];
	L3_bins[id][indx][2] += P->mass * angm[2];
	break;
      }
    }
  }

  thread_timing_end(id);

  return (NULL);
}


extern "C" {
  ExternalForce *makerAddMass(const YAML::Node& conf)
  {
    return new UserAddMass(conf);
  }
}

class proxyaddmass { 
public:
  proxyaddmass()
  {
    factory["useraddmass"] = makerAddMass;
  }
};

proxyaddmass p;
