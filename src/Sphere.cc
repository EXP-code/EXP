
#include <limits>

#include "expand.H"

#include <gaussQ.H>
#include <Sphere.H>
#include <plummer.H>
#include <interp.H>

extern void orthoTest(const std::vector<Eigen::MatrixXd>& tests,
		      const std::string& classname,
		      const std::string& indexname);

const std::set<std::string>
Sphere::valid_keys = {
  "rs",
  "numr",
  "nums",
  "cmap",
  "diverge",
  "dfac",
  "modelname",
  "cachename",
  "dtime",
  "logr",
  "plummer"
};

Sphere::Sphere(Component* c0, const YAML::Node& conf, MixtureBasis* m) :
  SphericalBasis(c0, conf, m)
{
  id = "Sphere SL";
				// Defaults
  rs = 0.067*rmax;
  numr = 2000;
  nums = 2000;
  cmap = 1;
  diverge = 0;
  dfac = 1.0;
  model_file = "SLGridSph.model";
  cache_file = "SLGridSph.cache";
  tnext = dtime = 0.0;
  recompute  = false;
  plummer    = true;
  logr       = false;

				// Get initialization info
  initialize();

				// Basis computation logic
  if (dtime>0.0) {
    recompute = true;
    tnext = tnow + dtime;
  }
				// Enable MPI code for more than one node
  if (numprocs>1) SLGridSph::mpi = 1;

  std::string modelname = model_file;
  std::string cachename = outdir  + cache_file + "." + runtag;

  id += ", model=" + modelname;

				// Generate Sturm-Liouville grid
  ortho = std::make_shared<SLGridSph>(modelname,
				      Lmax, nmax, numr, rmin, rmax, true,
				      cmap, rs, diverge, dfac, cachename);

				// Get the min and max expansion radii
  rmin  = ortho->getRmin();
  rmax  = ortho->getRmax();

  if (myid==0) {
    std::string sep("----    ");
    std::cout << "---- Sphere parameters: "
	      << std::endl << sep << "lmax="        << Lmax
	      << std::endl << sep << "nmax="        << nmax
	      << std::endl << sep << "cmap="        << cmap
	      << std::endl << sep << "rmin="        << rmin
	      << std::endl << sep << "rmax="        << rmax
	      << std::endl << sep << "logr="        << std::boolalpha << logr
	      << std::endl << sep << "NO_L0="       << std::boolalpha << NO_L0
	      << std::endl << sep << "NO_L1="       << std::boolalpha << NO_L1
	      << std::endl << sep << "EVEN_L="      << std::boolalpha << EVEN_L
	      << std::endl << sep << "EVEN_M="      << std::boolalpha << EVEN_M
	      << std::endl << sep << "M0_ONLY="     << std::boolalpha << M0_only
	      << std::endl << sep << "selfgrav="    << std::boolalpha << self_consistent
	      << std::endl;
  }


  setup();
}


void Sphere::initialize()
{
  // Remove matched keys
  //
  for (auto v : valid_keys) current_keys.erase(v);
  
  // Assign values from YAML
  //
  try {
    if (conf["rs"])        rs         = conf["rs"].as<double>();
    if (conf["numr"])      numr       = conf["numr"].as<int>();
    if (conf["nums"])      nums       = conf["nums"].as<int>();
    if (conf["cmap"])      cmap       = conf["cmap"].as<int>();
    if (conf["diverge"])   diverge    = conf["diverge"].as<int>();
    if (conf["dfac"])      dfac       = conf["dfac"].as<double>();
    if (conf["modelname"]) model_file = conf["modelname"].as<std::string>();
    if (conf["cachename"]) cache_file = conf["cachename"].as<std::string>();
    if (conf["dtime"])     dtime      = conf["dtime"].as<double>();
    if (conf["logr"])      logr       = conf["logr"].as<bool>();
    if (conf["plummer"])   plummer    = conf["plummer"].as<bool>();
  }
  catch (YAML::Exception & error) {
    if (myid==0) std::cout << "Error parsing parameters in Sphere: "
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

Sphere::~Sphere(void)
{
  // NADA
}


void Sphere::get_dpotl(int lmax, int nmax, double r, Eigen::MatrixXd& p,
		       Eigen::MatrixXd& dp, int tid)
{
  ortho->get_pot(p, r);
  ortho->get_force(dp, r);
}

void Sphere::get_potl(int lmax, int nmax, double r, Eigen::MatrixXd& p, int tid)
{
  ortho->get_pot(p, r);
}

void Sphere::get_dens(int lmax, int nmax, double r, Eigen::MatrixXd& p, int tid)
{
  ortho->get_dens(p, r);
}

void Sphere::get_potl_dens(int lmax, int nmax, double r, Eigen::MatrixXd& p,
			   Eigen::MatrixXd& d, int tid)
{
  ortho->get_pot(p, r);
  ortho->get_dens(d, r);
}

double Sphere::mapIntrp(const std::map<double, double> &data, double x)
{
  typedef std::map<double, double>::const_iterator it;

  it i=data.upper_bound(x);

  if (i==data.end()) {
    return (--i)->second;
  }

  if (i==data.begin()) {
    return i->second;
  }

  it l=i; --l;

  const double delta=(x- l->first)/(i->first - l->first);
  return delta*i->second +(1-delta)*l->second;
}

double Sphere::mapDeriv(const std::map<double, double> &data, double x)
{
  typedef std::map<double, double>::const_iterator it;

  it i=data.upper_bound(x);

  if (i==data.end()) {
    return (--i)->second;
  }

  if (i==data.begin()) {
    return i->second;
  }

  it l=i; --l;

  return (i->second - l->second)/(i->first - l->first);
}

void Sphere::make_model_bin()
{
  // Get total particle number
  //
  int N = component->Particles().size();
  std::vector<int> NN(numprocs);
  MPI_Allgather(&N, 1, MPI_INT, &NN[0], 1, MPI_INT, MPI_COMM_WORLD);

				// Allgatherv displacement vector
  std::vector<int> dspl(numprocs);
  dspl[0] = 0;
  int Ntot = NN[0];		// Total number of particles
  for (int n=1; n<numprocs; n++) {
    Ntot += NN[n];
    dspl[n] = dspl[n-1] + NN[n-1];
  }

  std::vector<double> rr, mm, RR(Ntot), MM(Ntot);

  for (auto v : component->Particles()) {
    // Get the particle index
    auto i = v.first;

    double r2 = 0.0;
      for (int k=0; k<3; k++) {
	double pos = component->Pos(i, k, Component::Local | Component::Centered);
	r2 += pos*pos;
      }
      rr.push_back(sqrt(r2));
      mm.push_back(component->Part(i)->mass);
  }
    
  // All processes get complete mass histogram
  //
  MPI_Allgatherv(&rr[0], N, MPI_DOUBLE, &RR[0], &NN[0], &dspl[0], MPI_DOUBLE, MPI_COMM_WORLD);

  MPI_Allgatherv(&mm[0], N, MPI_DOUBLE, &MM[0], &NN[0], &dspl[0], MPI_DOUBLE, MPI_COMM_WORLD);
  
  std::map<double, double> RM;
  for (int i=0; i<Ntot; i++) RM[RR[i]] = MM[i];
  
  // Cumulate mass
  //
  double totM = 0.0;
  for (auto & v : RM) {
    v.second += totM;
    totM = v.second;
  }

  // Make mass array
  //
  int    numR = std::min<int>(sqrt(Ntot), nums);
  double Rmin = std::max<double>(RM.begin()->first, rmin);
  double Rmax = RM.rbegin()->first;
  bool   logR = false;

  if (logr and Rmin>0.0) {
    Rmin = log(Rmin);
    Rmax = log(Rmax);
    logR = true;
  }

  double dR = (Rmax - Rmin)/numR;

  std::vector<double> R(numR), D(numR), M(numR), P(numR), P1(numR);

  for (int i=0; i<numR; i++) {

    double r1 = Rmin + dR*i;	// Bin edges
    double r2 = Rmin + dR*(i+1);
    R[i] = 0.5*(r1 + r2);	// Bin center
    double dr = dR;

    if (logR) {
      r1   = exp(r1);
      r2   = exp(r2);
      R[i] = exp(R[i]);
      dr   = r2 - r1;
    }
				// Mass in radial bin
    double mass = mapIntrp(RM, r2) - mapIntrp(RM, r1);

				// Bin density
    D[i] = mass / (4.0*M_PI/3.0 * (pow(r2, 3.0) - pow(r1, 3.0)));

    if (i==0) M[i] = 0.0;	// Consistent cumulative mass at edges
    else M[i] = M[i-1] + 2.0*M_PI*dr*(D[i]*r2*r2 + D[i-1]*r1*r1);
    
				// Outer potential
    if (i==0) P1[i] = 4.0*M_PI*dr*D[i]*r2;
    else P1[i] = P1[i-1] + 2.0*M_PI*dr*(D[i]*r2 + D[i-1]*r1);
  }

				// Full potential
  for (int i=numR-1; i>0; i--) 
    P[i] = M[i]/R[i] + P1[numR-1] - 0.5*(P1[i] + P1[i-1]);
  P[0] = M[0]/R[0] + P1[numR-1] - 0.5*P1[0];

  for (int i=0; i<numR; i++) P[i] *= -1.0;

  // Debug output
  //
  if (myid==0) {
    static int cnt = 0;
    std::ostringstream sout; sout << "SphereMassModel." << cnt++;
    std::ofstream dbg(sout.str());
    if (dbg.good()) {
      for (int i=0; i<numR; i++)
	dbg << std::setw(18) << R[i]
	    << std::setw(18) << D[i]
	    << std::setw(18) << M[i]
	    << std::setw(18) << P[i]
	    << std::setw(18) << mapIntrp(RM, R[i])
	    << std::setw(18) << mapDeriv(RM, R[i])/(4.0*M_PI*R[i]*R[i])
	    << std::endl;
    }
  }

  // Create a new spherical model
  //
  SphModTblPtr mod = std::make_shared<SphericalModelTable>(numR, &R[0]-1, &D[0]-1, &M[0]-1, &P[0]-1);

  // Regenerate Sturm-Liouville grid
  //
  std::string cachename = outdir  + cache_file + "." + runtag;
  ortho = std::make_shared<SLGridSph>(mod, Lmax, nmax, numR, Rmin, Rmax, false, 1, 1.0, cachename);

  // Test for basis consistency
  //
  orthoTest(ortho->orthoCheck(std::max<int>(nmax*5, 100)), "Sphere", "l");

  // Update time trigger
  //
  tnext = tnow + dtime;
}


void Sphere::make_model_plummer()
{
  // Get total particle number
  //
  int N = component->Particles().size();
  std::vector<int> NN(numprocs);
  MPI_Allgather(&N, 1, MPI_INT, &NN[0], 1, MPI_INT, MPI_COMM_WORLD);

				// Allgatherv displacement vector
  std::vector<int> dspl(numprocs);
  dspl[0] = 0;
  int Ntot = NN[0];		// Total number of particles
  for (int n=1; n<numprocs; n++) {
    Ntot += NN[n];
    dspl[n] = dspl[n-1] + NN[n-1];
  }

  std::vector<double> rr, mm, RR, MM;

  if (myid==0) {
    RR.resize(Ntot);
    MM.resize(Ntot);
  }

  for (auto v : component->Particles()) {
    // Get the particle index
    auto i = v.first;

    double r2 = 0.0;
      for (int k=0; k<3; k++) {
	double pos = component->Pos(i, k, Component::Local | Component::Centered);
	r2 += pos*pos;
      }
      rr.push_back(sqrt(r2));
      mm.push_back(component->Part(i)->mass);
  }
    
  // All processes get complete mass histogram
  //
  MPI_Gatherv(&rr[0], N, MPI_DOUBLE, &RR[0], &NN[0], &dspl[0], MPI_DOUBLE, 0, MPI_COMM_WORLD);

  MPI_Gatherv(&mm[0], N, MPI_DOUBLE, &MM[0], &NN[0], &dspl[0], MPI_DOUBLE, 0, MPI_COMM_WORLD);

  double pScl, rMax, rMin, totM=0.0;
  
  if (myid==0) {
    typedef std::pair<double, double> DP;
    std::vector<DP> MR(Ntot);
    for (int i=0; i<Ntot; i++) MR[i] = {RR[i], MM[i]};
    std::sort(MR.begin(), MR.end());
  
    // Cumulate mass
    //
    for (auto & v : MR) {
      v.second += totM;
      totM = v.second;
    }

    // Compute half mass
    //
    DP target(0.0, 0.5*totM);
    auto it = std::lower_bound(MR.begin(), MR.end(), target, 
			       [](DP lhs, DP rhs) -> bool { return lhs.second < rhs.second; });

    auto lt=it; --lt;

    const double delta=(0.5*totM - lt->second)/(it->second - lt->second);

    double hmass = delta*it->first +(1 - delta)*lt->first;

    pScl = 0.7664209365408798 * hmass;
    rMin = MR.begin()->first;
    rMax = MR.rbegin()->first;
  }

  MPI_Bcast(&pScl, 1, MPI_DOUBLE, 0, MPI_COMM_WORLD);
  MPI_Bcast(&rMin, 1, MPI_DOUBLE, 0, MPI_COMM_WORLD);
  MPI_Bcast(&rMax, 1, MPI_DOUBLE, 0, MPI_COMM_WORLD);
  MPI_Bcast(&totM, 1, MPI_DOUBLE, 0, MPI_COMM_WORLD);

  // Make mass array
  //
  double Rmin = std::max<double>(rMin, rmin);
  double Rmax = rMax;
  bool   logR = false;

  PlummerSphere ps(pScl, Rmin, Rmax, totM);

  if (logr and Rmin>0.0) {
    Rmin = log(Rmin);
    Rmax = log(Rmax);
    logR = true;
  }

  double dR = (Rmax - Rmin)/nums;

  std::vector<double> R(nums), D(nums), M(nums), P(nums), P1(nums);

  for (int i=0; i<nums; i++) {
				// Compute the radius
    R[i] = Rmin + dR*i;
    if (logR) R[i] = exp(R[i]);
				// Compute the field quantities
				// analytically
    D[i] = ps.get_density(R[i]);
    M[i] = ps.get_mass(R[i]);
    P[i] = ps.get_pot(R[i]);
  }

  // Debug output
  //
  if (myid==0) {
    static int cnt = 0;
    std::ostringstream sout; sout << "SphereMassModel." << cnt++;
    std::ofstream dbg(sout.str());
    if (dbg.good()) {
      for (int i=0; i<nums; i++)
	dbg << std::setw(18) << R[i]
	    << std::setw(18) << D[i]
	    << std::setw(18) << M[i]
	    << std::setw(18) << P[i]
	    << std::endl;
    }
  }

  // Create a new spherical model
  //
  SphModTblPtr mod = std::make_shared<SphericalModelTable>(nums, &R[0]-1, &D[0]-1, &M[0]-1, &P[0]-1);

  // Regenerate Sturm-Liouville grid
  //
  std::string cachename = outdir  + cache_file + "." + runtag;
  ortho = std::make_shared<SLGridSph>(mod, Lmax, nmax, numr, Rmin, Rmax, false, 1, 1.0, cachename);

  // Update time trigger
  //
  tnext = tnow + dtime;
}
