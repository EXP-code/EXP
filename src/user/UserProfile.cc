#include <cmath>
#include <sstream>

#include <expand.H>
#include <localmpi.H>

#include <ExternalCollection.H>
#include <Basis.H>
#include <UserProfile.H>

#include <pthread.h>  

using namespace std;

UserProfile::UserProfile(const YAML::Node& conf) : ExternalForce(conf)
{
  first    = true;
  filename = "profile";
  NUMR     = 1000;
  RMIN     = 1.0e-4;
  RMAX     = 2.0;
  DT       = 0.01;
  LOGR     = true;
  NTHETA   = 16;
  NPHI     = 16;

  initialize();

  if (numComp==0) {
    if (myid==0) cerr << "You must specify component targets!" << endl;
    MPI_Abort(MPI_COMM_WORLD, 120);
  }

				// Search for component by name
  for (int i=0; i<numComp; i++) {
    bool found = false;
    for (auto c : comp->components) {
      if ( !C[i].compare(c->name) ) {
	// Check to see that the force can return field values
	if (dynamic_cast <Basis*> (c->force)) {
	  c0.push_back(c);
	  found = true;
	  break;
	} else {
	  cerr << "Process " << myid << ": desired component <"
	       << C[i] << "> is not a basis type!" << endl;

	  MPI_Abort(MPI_COMM_WORLD, 121);
	}
      }
    }

    if (!found) {
      cerr << "Process " << myid << ": can't find desired component <"
	   << C[i] << ">" << endl;
    }
  }

  userinfo();

  // -----------------------------------------------------------
  // Make grid spacing
  // -----------------------------------------------------------

  if (RMIN <= 0.0) {
    LOGR = false;
    RMIN = 0.0;
  }

  if (LOGR) {
    RMIN = log(RMIN);
    RMAX = log(RMAX);
  }

  dR = (RMAX - RMIN)/(NUMR-1);

  int nstep = NUMR/numprocs;

  // -----------------------------------------------------------
  // Compute the beginning, ending, and grid counts assigned to 
  // each node
  // -----------------------------------------------------------

  for (int n=0; n<numprocs; n++) {
    nbeg.push_back(nstep*n);
    nend.push_back(nstep*(n+1));
    numb.push_back(nstep);
  }
  nend[numprocs-1] = NUMR;
  numb[numprocs-1] = nend[numprocs-1] - nbeg[numprocs-1];

  // -----------------------------------------------------------
  // Assign array storage
  // -----------------------------------------------------------

  rho1  = vector<double>(numb[myid]);
  mass1 = vector<double>(numb[myid]);
  pot1  = vector<double>(numb[myid]);

  if (myid==0) {
    rho  = vector<double>(NUMR);
    mass = vector<double>(NUMR);
    pot  = vector<double>(NUMR);
  }
  
}

UserProfile::~UserProfile()
{
  // Nothing so far
}

void UserProfile::userinfo()
{
  if (myid) return;		// Return if node master node

  print_divider();

  cout << "** User routine PROFILE initialized"
       << " with Components <" << C[0];
  for (int i=1; i<numComp; i++) cout << " " << C[i];
  cout << ">";
  
  cout << ", NUMR="     << NUMR
       << ", RMIN="     << RMIN
       << ", RMAX="     << RMAX
       << ", DT="       << DT
       << ", LOGR="     << (LOGR ? "True" : "False")
       << ", NTHETA="   << NTHETA
       << ", NPHI="     << NPHI
       << ", filename=" << filename
       << endl;
  
  print_divider();
}

const std::set<std::string>
UserProfile::valid_keys = {
  "filename",
  "NUMR",
  "RMIN",
  "RMAX",
  "DT",
  "NTHETA",
  "NPHI",
  "LOGR"
};

void UserProfile::initialize()
{
  // Check for unmatched keys
  //
  auto unmatched = YamlCheck(conf, valid_keys);
  if (unmatched.size())
    throw YamlConfigError("UserProfile", "parameter", unmatched,
			  __FILE__, __LINE__);

  // Assign parameters from YAML config
  //
  try {

    for (numComp=0; numComp<1000; numComp++) {
      ostringstream count;
      count << "C(" << numComp+1 << ")";
      if (conf[count.str()])
	C.push_back(conf[count.str()].as<std::string>());
      else break;
    }
    
    if (numComp != (int)C.size()) {
      cerr << "UserProfile: error parsing component names, "
	   << "  Size(C)=" << C.size()
	   << "  numRes=" << numComp << endl;
      MPI_Abort(MPI_COMM_WORLD, 122);
    }
    
    if (conf["filename"])       filename           = conf["filename"].as<string>();
    if (conf["NUMR"])           NUMR               = conf["NUMR"].as<int>();
    if (conf["RMIN"])           RMIN               = conf["RMIN"].as<double>();
    if (conf["RMAX"])           RMAX               = conf["RMAX"].as<double>();
    if (conf["DT"])             DT                 = conf["DT"].as<double>();
    if (conf["NTHETA"])         NTHETA             = conf["NTHETA"].as<int>();
    if (conf["NPHI"])           NPHI               = conf["NPHI"].as<int>();
    
    if (conf["LOGR"]) {
      std::string val = conf["LOGR"].as<std::string>();
      if (val[0]=='T' || val[0]=='t') 
	LOGR = true;
      else if (val[0]=='F' || val[0]=='f')
	LOGR = false;
      else if (atoi(val.c_str()))
	LOGR = true;
      else
	LOGR = false;
    }
  }
  catch (YAML::Exception & error) {
    if (myid==0) std::cout << "Error parsing parameters in UserProfile: "
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


void UserProfile::determine_acceleration_and_potential(void)
{

  if (multistep and mlevel>0) return;

#if HAVE_LIBCUDA==1		// Cuda compatibility
  getParticlesCuda(cC);
#endif

  if (first) {

    count = 0;
    tnext = tnow;

    if (restart) {

      if (myid == 0) {

	for (count=0; count<10000; count++) {
	  ostringstream ostr;
	  ostr << outdir << runtag << "." << filename << "." << count;
	  
	  //
	  // Try to open stream for writing
	  //
	  ifstream in(ostr.str().c_str());
	  if (!in) break;
	}
      }

      if (myid==0) 
	cout << "UserProfile: beginning at frame=" << count << endl;

      MPI_Bcast(&count, 1, MPI_INT, 0, MPI_COMM_WORLD);
    }

    first = false;
  }
  
  // -----------------------------------------------------------
  // Let the root node make the decision to create the profile
  // for consistency
  // -----------------------------------------------------------

  unsigned char doit = 0;
  if (tnow >= tnext) doit = 1;
  MPI_Bcast(&doit, 1, MPI_UNSIGNED_CHAR, 0, MPI_COMM_WORLD);

  if (doit) {

    // -----------------------------------------------------------
    // Clean the data store before the first component in the list
    // -----------------------------------------------------------

    if (cC == c0.front()) {
      for (int i=0; i<numb[myid]; i++) rho1[i] = mass1[i] = pot1[i] = 0.0;
    }

    // -----------------------------------------------------------
    // Compute the profile
    // -----------------------------------------------------------

    double r, theta, phi, x;
    double dens0, potl0, dens, potl, potr, pott, potp, mass0;
    
    for (int i=0; i<numb[myid]; i++) {

      r = RMIN + dR*(nbeg[myid]+i);
      if (LOGR)	r = exp(r);

      mass0 = 0.0;

      for (int j=0; j<NTHETA; j++) {
	x = -1.0 + 2.0*(0.5+j)/NTHETA;
	theta = acos(x);

	for (int k=0; k<NPHI; k++) {
	  phi = 2.0*M_PI*k/NPHI;

	  ((Basis *)cC->force)->
	    determine_fields_at_point_sph(r, theta, phi,
					  &dens0, &potl0, 
					  &dens, &potl,
					  &potr, &pott, &potp);
	  mass0 += potr*r*r;
	}
      }

      rho1 [i] = -dens0;
      mass1[i] = mass0/(NTHETA*NPHI);
      pot1 [i] = potl0;
    }

    
    // -----------------------------------------------------------
    // Gather the data and print the profile after the last comp
    // -----------------------------------------------------------
    
    if (cC == c0.back()) {

      MPI_Gatherv(&rho1[0], numb[myid], MPI_DOUBLE,
		  &rho[0],  &numb[0], &nbeg[0], MPI_DOUBLE, 0, MPI_COMM_WORLD);

      MPI_Gatherv(&mass1[0], numb[myid], MPI_DOUBLE,
		  &mass[0], &numb[0], &nbeg[0], MPI_DOUBLE, 0, MPI_COMM_WORLD);

      MPI_Gatherv(&pot1[0], numb[myid], MPI_DOUBLE,
		  &pot[0],  &numb[0], &nbeg[0], MPI_DOUBLE, 0, MPI_COMM_WORLD);

      
      //
      // Print the file
      //
      if (myid==0) {		
	  
	ostringstream ostr;
	ostr << outdir << runtag << "." << filename << "." << count;
	  
	//
	// Try to open stream for writing
	//
	ofstream out(ostr.str().c_str());
	if (out) {

	  out << "# Profile at T=" << tnow << endl
	      << "#" << setw(19) << right << " 1) = r"
	      << setw(20) << "2) = rho" 
	      << setw(20) << "3) = M(r)"
	      << setw(20) << "4) U(r)" << endl
	      << setw(10) << NUMR << endl << setprecision(12);

	  for (int i=0; i<NUMR; i++) {

	    double r = RMIN + dR*i;
	    if (LOGR) r = exp(r);

	    out << setw(20) << r
		<< setw(20) << rho[i]
		<< setw(20) << mass[i]
		<< setw(20) << pot[i]
		<< endl;
	  }
	  
	} else {
	  cerr << "UserProfile: error opening <" << ostr.str() << ">" << endl;
	}
      }
      
      count++;
      tnext = tnow + DT;
    }

  }

}


extern "C" {
  ExternalForce *makerProfile(const YAML::Node& conf)
  {
    return new UserProfile(conf);
  }
}

class proxyprof { 
public:
  proxyprof()
  {
    factory["userprofile"] = makerProfile;
  }
};

static proxyprof p;
