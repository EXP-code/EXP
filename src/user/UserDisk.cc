#include <math.h>
#include <sstream>

#include "expand.H"
#include <gaussQ.H>

#include <UserDisk.H>

extern double bessj0(double);
extern double bessj1(double);

const std::set<std::string>
UserDisk::valid_keys = {
  "ctrname",
  "a",
  "mass",
  "Ton",
  "Toff",
  "DeltaT",
  "Nscale",
  "Ngrid",
  "Nint",
  "debug",
  "dfac"
};

UserDisk::UserDisk(const YAML::Node& conf) : ExternalForce(conf)
{
  id = "ThinExponentialDiskPotential";

  a = 1.0;			// Disk scale length
  mass = 1.0;			// Total disk mass

  Ton = -20.0;			// Turn on start time
  Toff = 200.0;			// Turn off start time
  DeltaT = 1.0;			// Turn on duration

  Nscale = 25.0;		// Maximum grid radius in scale lengths
  Ngrid = 800;			// Number of points on grid
  Nint = 600;			// Number of k-integration points

  debug = false;		// Print out potential/force tables
  dfac = 1.2;			// Fraction of grid for interp test

  ctr_name = "";		// Default component for com
  
  initialize();

  if (ctr_name.size()>0) {
				// Look for the fiducial component for
				// centering
    bool found = false;
    for (auto c : comp->components) {
      if ( !ctr_name.compare(c->name) ) {
	c0 = c;
	found = true;
      break;
      }
    }

    if (!found) {
      std::ostringstream sout;
      sout << "Can't find desired component <" << ctr_name << ">";
      throw GenericError(sout.str(), __FILE__, __LINE__, 35, false);
    }

  }
  else
    c0 = NULL;


  userinfo();

  genTable();
}

UserDisk::~UserDisk()
{
  delete [] Rtable;
  delete [] Ztable;
  delete [] Ptable;
}

void UserDisk::userinfo()
{
  if (myid) return;		// Return if node master node

  print_divider();

  cout << "** User routine: thin exponential disk with a=" << a 
       << ", mass=" << mass
       << ", Nscale=" << Nscale
       << ", Ngrid=" << Ngrid;

  if (c0) 
    cout << ", center on component <" << ctr_name << ">";
  else
    cout << ", using inertial center";

  if (debug)
    cout << ", debugging is *ON*";

  cout << endl;

  print_divider();
}

void UserDisk::initialize()
{
  // Check for unmatched keys
  //
  auto unmatched = YamlCheck(conf, valid_keys);
  if (unmatched.size())
    throw YamlConfigError("UserDisk", "parameter", unmatched,
			  __FILE__, __LINE__);

  // Assign values from YAML
  //
  try {
    if (conf["ctrname"])   ctr_name    = conf["ctrname"].as<string>();
    if (conf["a"])         a           = conf["a"].as<double>();
    if (conf["mass"])      mass        = conf["mass"].as<double>();
    if (conf["Ton"])       Ton         = conf["Ton"].as<double>();
    if (conf["Toff"])      Toff        = conf["Toff"].as<double>();
    if (conf["DeltaT"])    DeltaT      = conf["DeltaT"].as<double>();

    if (conf["Nscale"])    Nscale      = conf["Nscale"].as<double>();
    if (conf["Ngrid"])     Ngrid       = conf["Ngrid"].as<int>();
    if (conf["Nint"])      Nint        = conf["Nint"].as<int>();

    if (conf["debug"])     debug       = conf["debug"].as<bool>();
    if (conf["dfac"])      dfac        = conf["dfac"].as<double>();
  }
  catch (YAML::Exception & error) {
    if (myid==0) std::cout << "Error parsing parameters in UserDisk: "
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


void UserDisk::getTable(double R, double Z,
			 double &pot, double &fr, double &fz)
{
  double RR = fabs(R);
  double ZZ = fabs(Z);

  if (RR>=Rmax || ZZ>=Zmax) {
    
    double r = sqrt(R*R + Z*Z);
    
    pot = -mass/r;
    fr = -mass * R/(r*r*r);
    fz = -mass * Z/(r*r*r);

  } else {
    int indR = min<int>(Ngrid-2, (int)(RR/dR));
    int indZ = min<int>(Ngrid-2, (int)(ZZ/dZ));

    double aR[2], aZ[2];
    
    aR[0] = (dR*(indR+1)-RR)/dR;
    aR[1] = 1.0 - aR[0];
    
    aZ[0] = (dZ*(indZ+1)-ZZ)/dZ;
    aZ[1] = 1.0 - aZ[0];
    
    pot = fr = fz = 0.0;	// Bilinear interpolation
    for (int i=0; i<2; i++) {
      for (int j=0; j<2; j++) {
	pot += Ptable[(indR+i)*Ngrid+indZ+j]*aR[i]*aZ[j];
	fr  += Rtable[(indR+i)*Ngrid+indZ+j]*aR[i]*aZ[j];
	fz  += Ztable[(indR+i)*Ngrid+indZ+j]*aR[i]*aZ[j];
      }
    }
    if (Z<0.0) fz *= -1.0;
  }

}

void UserDisk::genTable()
{
  Ptable = new double [Ngrid*Ngrid];
  Rtable = new double [Ngrid*Ngrid];
  Ztable = new double [Ngrid*Ngrid];

  Rmax = a*Nscale;
  Zmax = a*Nscale;

  dR = Rmax/(Ngrid-1);
  dZ = Zmax/(Ngrid-1);

  double R, Z, Q, K, ansP, ansR, ansZ, fac, b0, b1;

				// Compute table for upper quarter
  				// plane
  LegeQuad lq(Nint);

  for (int i=0; i<Ngrid; i++) {

    R = dR*i;

    for (int j=0; j<Ngrid; j++) {

      Z = dZ*j;

      ansP = ansR = ansZ = 0.0;
      for (int k=0; k<Nint; k++) {

				// Do integral over finite domain
				// through change of variables

				// B&T equation 2-167

	Q = lq.knot(k)/a;
	K = Q/sqrt(1.0 - Q*Q*a*a);

	fac = exp(-K*Z) * lq.weight(k)*mass/a;
	b0 = bessj0(K*R);
	b1 = bessj1(K*R);

	ansP += -b0 * fac;	// Potential
	ansR += -K * b1 * fac;	// -d(Potential)/dR
	ansZ += -K * b0 * fac;	// -d(Potential)/dZ
      }

      Ptable[i*Ngrid + j] = ansP;
      Rtable[i*Ngrid + j] = ansR;
      Ztable[i*Ngrid + j] = ansZ;
    }
  }


  if (debug) printTable();

}


// This routine is for debugging 
// (enabled by setting parameter debug=1)
//
void UserDisk::printTable()
{
  string filename;

  filename = outdir + "test_pot." + runtag;
  ofstream outP(filename.c_str());
  filename = outdir + "test_fr." + runtag;
  ofstream outR(filename.c_str());
  filename = outdir + "test_fz." + runtag;
  ofstream outZ(filename.c_str());

  double R, Z;

  for (int i=0; i<Ngrid; i++) {

    R = dR*i;

    for (int j=0; j<Ngrid; j++) {

      Z = dZ*j;

      outP << setw(18) << R
	   << setw(18) << Z
	   << setw(18) << Ptable[i*Ngrid+j]
	   << endl;

      outR << setw(18) << R
	   << setw(18) << Z
	   << setw(18) << Rtable[i*Ngrid+j]
	   << endl;

      outZ << setw(18) << R
	   << setw(18) << Z
	   << setw(18) << Ztable[i*Ngrid+j]
	   << endl;
    }

    outP << endl;
    outR << endl;
    outZ << endl;
  }


  outP.close();
  outR.close();
  outZ.close();

  filename = outdir + "test_pot1." + runtag;
  outP.open(filename.c_str());
  filename = outdir + "test_fr1." + runtag;
  outR.open(filename.c_str());
  filename = outdir + "test_fz1." + runtag;
  outZ.open(filename.c_str());

  const int num = 100;
  double dr = dfac*Rmax/(num-1);
  double dz = 2.0*dfac*Zmax/(num-1);
  double pot, fr, fz;

  for (int i=0; i<num; i++) {

    R = dr*i;

    for (int j=0; j<num; j++) {

      Z = -dfac*Zmax + dz*j;

      getTable(R, Z, pot, fr, fz);

      outP << setw(18) << R
	   << setw(18) << Z
	   << setw(18) << pot
	   << endl;

      outR << setw(18) << R
	   << setw(18) << Z
	   << setw(18) << fr
	   << endl;

      outZ << setw(18) << R
	   << setw(18) << Z
	   << setw(18) << fz
	   << endl;
    }

    outP << endl;
    outR << endl;
    outZ << endl;
  }

}


void UserDisk::determine_acceleration_and_potential(void)
{
#if HAVE_LIBCUDA==1		// Cuda compatibility
  getParticlesCuda(cC);
#endif

  exp_thread_fork(false);
}


void * UserDisk::determine_acceleration_and_potential_thread(void * arg) 
{
  unsigned nbodies = cC->Number();
  int id = *((int*)arg);
  int nbeg = nbodies*id/nthrds;
  int nend = nbodies*(id+1)/nthrds;

  double xx, yy, zz, rr;
  double pot=0.0, fr=0.0, fz=0.0;
  vector<double> pos(3);

  double amp = 
      0.5*(1.0 + erf( (tnow - Ton )/DeltaT ))
    * 0.5*(1.0 - erf( (tnow - Toff)/DeltaT )) ;


  PartMapItr it = cC->Particles().begin();
  unsigned long i;

  for (int q=0   ; q<nbeg; q++) it++;
  for (int q=nbeg; q<nend; q++) {
    i = (it++)->first;
				// If we are multistepping, compute accel 
				// only at or below this level

    if (multistep && (cC->Part(i)->level < mlevel)) continue;
    

				// Set center if component is
				// defined, otherwise use origin
    if (c0)
      for (int k=0; k<3; k++) 
	pos[k] = cC->Pos(i, k) - c0->center[k];
    else
      for (int k=0; k<3; k++) 
	pos[k] = cC->Pos(i, k);
    
    xx = pos[0];
    yy = pos[1];
    zz = pos[2];

    rr = sqrt( xx*xx + yy*yy );

				// Interpolate on table
    getTable(rr, zz, pot, fr, fz);

				// Add acceleration by disk
    cC->AddAccExt(i, 0, amp * fr*xx/(rr+1.0e-10) );
    cC->AddAccExt(i, 1, amp * fr*yy/(rr+1.0e-10) );
    cC->AddAccExt(i, 2, amp * fz );

				// Add external potential
    cC->AddPotExt(i, pot);

  }

  return (NULL);
}


extern "C" {
  ExternalForce *makerExpDisk(const YAML::Node& conf)
  {
    return new UserDisk(conf);
  }
}

class proxydisk { 
public:
  proxydisk()
  {
    factory["userdisk"] = makerExpDisk;
  }
};

static proxydisk p;
