#include <math.h>
#include <sstream>

#include "expand.H"
#include <localmpi.H>
#include <UserEBar.H>
#include <Timer.H>
static Timer timer_tot, timer_thrd;
static bool timing = false;


const std::set<std::string>
UserEBar::valid_keys = {
  "ctrname",
  "angmname",
  "tblname",
  "length",
  "bratio",
  "cratio",
  "amp",
  "angmomfac",
  "barmass",
  "Ton",
  "Toff",
  "TmonoOn",
  "TmonoOff",
  "DeltaT",
  "DeltaMonoT",
  "DOmega",
  "dtom",
  "T0",
  "Fcorot",
  "omega",
  "fixed",
  "self",
  "oscil",
  "Ofreq",
  "Oamp",
  "alpha",
  "x0",
  "y0",
  "z0",
  "u0",
  "v0",
  "w0",
  "monopole",
  "follow",
  "onoff",
  "monofrac",
  "quadfrac",
  "filename",
  "ctrname",
  "angmname",
  "tblname",
  "length",
  "bratio",
  "cratio",
  "amp",
  "angmomfac",
  "barmass",
  "Ton",
  "Toff",
  "TmonoOn",
  "TmonoOff",
  "DeltaT",
  "DeltaMonoT",
  "DOmega",
  "dtom",
  "T0",
  "Fcorot",
  "omega",
  "fixed",
  "self",
  "oscil",
  "Ofreq",
  "Oamp",
  "alpha",
  "x0",
  "y0",
  "z0",
  "u0",
  "v0",
  "w0",
  "monopole",
  "follow",
  "onoff",
  "monofrac",
  "quadfrac",
  "filename"
};

UserEBar::UserEBar(const YAML::Node& conf) : ExternalForce(conf)
{
  id = "RotatingBarWithMonopole";

  length            = 1.0;	// Bar length
  bratio            = 0.5;	// Ratio of b to a
  cratio            = 0.1;	// Ratio of c to b
  amplitude         = 0.3;	// Bar quadrupole amplitude
  angmomfac         = 1.0;	// Artifically change the total bar ang mom
  barmass           = 1.0;	// Total bar mass
  Ton               = -20.0;	// Turn on start time
  Toff              = 200.0;	// Turn off start time
  TmonoOn           = -20.0;	// Turn on start time for monopole
  TmonoOff          = 200.0;	// Turn off start time monopole
  DeltaT            = 1.0;	// Turn on duration
  DeltaMonoT        = 1.0;	// Turn on duration for monopole
  DOmega            = 0.0;	// Change in pattern speed
  dtom              = -1.0;	// Width of forced bar slow down
  T0                = 10.0;	// Center of pattern speed change
  Fcorot            = 1.0;	// Corotation factor
  fixed             = false;	// Constant pattern speed
  alpha             = 5.0;	// Variable sharpness bar potential
  table             = false;	// Not using tabled quadrupole
  monopole          = true;	// Use the monopole part of the potential
  monopole_follow   = true;	// Follow monopole center
  monopole_onoff    = false;	// To apply turn-on and turn-off to monopole
  monopole_frac     = 1.0;	// Fraction of monopole to turn off
  quadrupole_frac   = 1.0;	// Fraction of quadrupole to turn off
  mupdate           = 0;	// Report bar diagnostics for levels <= value

  oscil             = false;	// Oscillate quadrupole amplitude
  Oamp              = 0.5;	// Relative strength of amplitude oscillations
  Ofreq             = 2.0;	// Frequency of amplitude oscillations

				// Output file name
  filename = "BarRot." + runtag;

  firstime = true;
  omega0 = -1.0;

  ctr_name = "";		// Default component for com
  angm_name = "";		// Default component for angular momentum
  table_name = "";		// Default for input b1,b5 table
  
  ellip = 0;

  // Zero monopole variables
  teval = vector<double>(multistep+1, tnow);
  for (int k=0; k<3; k++) bps[k] = vel[k] = acc[k] = 0.0;

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
      sout << "Process " << myid << ": can't find desired component <"
	   << ctr_name << ">";
      throw GenericError(sout.str(), __FILE__, __LINE__, 35, false);
    }

  }
  else
    c0 = NULL;

  if (angm_name.size()>0) {
				// Look for the fiducial component
				// for angular momentum
    bool found = false;
    for (auto c : comp->components) {
      if ( !angm_name.compare(c->name) ) {
	c1 = c;
	found = true;
      break;
      }
    }

    if (!found) {
      std::ostringstream sout;
      sout << "Process " << myid << ": can't find desired component <"
	   << ctr_name << ">";
      throw GenericError(sout.str(), __FILE__, __LINE__, 35, false);
    }

  }
  else
    c1 = NULL;


  if (table_name.size()>0) {
				// Read in data
    ifstream in(string(outdir+table_name).c_str());
    if (!in) {
      std::ostringstream sout;
      sout << "Process " << myid << ": error opening quadrupole file <"
	   << outdir+table_name << ">";
      throw GenericError(sout.str(), __FILE__, __LINE__, 35, false);
    }
    
    string fline;
    double val;
    getline(in, fline, '\n');
    while (in) {
      if (fline.find("#") == string::npos && fline.find("!") == string::npos) {
	istringstream ins(fline.c_str());
	ins >> val;
	timeq.push_back(val);
	ins >> val;
	ampq.push_back(val*2.0);
	ins >> val;
	b5q.push_back(val);
      }

      getline(in, fline, '\n');
    }

    // Temporary debug
    if (myid==0) {
      cout << endl << "***Quadrupole***" << endl;
      for (unsigned i=0; i<timeq.size(); i++)
	cout << setw(5) << i
	     << setw(20) << timeq[i]
	     << setw(20) << ampq[i]
	     << setw(20) << b5q[i]
	     << endl;
    }

    qlast = timeq.size()-1;
    table = true;
  }

  // Zero monopole variables
  teval = vector<double>(multistep+1, tnow);

  // Assign working vectors for each thread
  tacc = new double* [nthrds];
  for (int n=0; n<nthrds; n++) tacc[n] = new double [3];

  userinfo();

  // Log file name

  name = outdir + filename + ".barstat";

  // Only turn on bar timing for extreme debugging levels
  if (VERBOSE>49) timing = true;
}

UserEBar::~UserEBar()
{
  for (int n=0; n<nthrds; n++) delete [] tacc[n];
  delete [] tacc;
  delete ellip;
}

void UserEBar::userinfo()
{
  if (myid) return;		// Return if node master node

  print_divider();

  cout << "** User routine ROTATING BAR with MONOPOLE initialized, " ;
  if (fixed) {
    cout << "prescribed pattern speed with"
	 << " domega=" <<  DOmega << " and t0=" << T0;
    if (dtom>0) cout << ", dT_om=" << dtom;
    cout << ", ";
  }
  else
    cout << "initial corotation fraction, ";
  cout << " amplitude=" << amplitude << ", ";

  if (fabs(angmomfac-1.0)>1.0e-10)
    cout << " ang mom factor=" << angmomfac << ", ";
    
  if (omega0<0.0)
    cout << "initial pattern speed to be computed, ";
  else
    cout << "initial pattern speed " << omega0 << ", ";

  cout << "quadrupole fraction=" << quadrupole_frac
       << ", Ton=" << Ton << ", Toff=" << Toff << ", DeltaT=" << DeltaT 
       << ", ";
  if (monopole) {
    if (monopole_onoff)
      cout << "using monopole with turn-on/off with fraction=" 
	   << monopole_frac << ", TmonoOn=" << TmonoOn
	   << ", TmonoOff=" << TmonoOff << ", DeltaMonoT=" << DeltaMonoT
	   << ", ";
    else
      cout << "using monopole, ";
    cout << "(x, y, z)=(" << bps[0] << ", " << bps[1] << ", " << bps[2] << "), ";
    cout << "(u, v, w)=(" << vel[0] << ", " << vel[1] << ", " << vel[2] << "), ";
    if (monopole_follow)
      cout << "self-consistent monopole centering, ";
    else
      cout << "monopole center fixed, ";
  }
  else
    cout << "without monopole, ";

  cout << "bar softness (alpha)=" << alpha << ", ";

  if (oscil)
    cout << "oscillation freq=" << Ofreq << " and oscillation amp=" << Oamp << ", ";

  if (mupdate<0)
    cout << "no log output, ";
  else 
    if (multistep)
      cout << "log output at levels " << mupdate << " and below, ";
  if (c0) 
    cout << "center on component <" << ctr_name << ">, ";
  else
    cout << "using inertial center, ";
  if (table)
    cout << "using user quadrupole table, ";
  if (c1) 
    cout << "angular momentum from <" << angm_name << ">" << endl;
  else
    cout << "no initial angular momentum (besides bar)" << endl;

  print_divider();
}

void UserEBar::initialize()
{
  // Check for unmatched keys
  //
  auto unmatched = YamlCheck(conf, valid_keys);
  if (unmatched.size())
    throw YamlConfigError("UserEbar", "parameter", unmatched,
			  __FILE__, __LINE__);

  // Assign values from YAML
  //
  try {
    if (conf["ctrname"])        ctr_name           = conf["ctrname"].as<string>();
    if (conf["angmname"])       angm_name          = conf["angmname"].as<string>();
    if (conf["tblname"])        table_name         = conf["tblname"].as<string>();
    if (conf["length"])         length             = conf["length"].as<double>();
    if (conf["bratio"])         bratio             = conf["bratio"].as<double>();
    if (conf["cratio"])         cratio             = conf["cratio"].as<double>();
    if (conf["amp"])            amplitude          = conf["amp"].as<double>();
    if (conf["angmomfac"])      angmomfac          = conf["angmomfac"].as<double>();
    if (conf["barmass"])        barmass            = conf["barmass"].as<double>();
    if (conf["Ton"])            Ton                = conf["Ton"].as<double>();
    if (conf["Toff"])           Toff               = conf["Toff"].as<double>();
    if (conf["TmonoOn"])        TmonoOn            = conf["TmonoOn"].as<double>();
    if (conf["TmonoOff"])       TmonoOff           = conf["TmonoOff"].as<double>();
    if (conf["DeltaT"])         DeltaT             = conf["DeltaT"].as<double>();
    if (conf["DeltaMonoT"])     DeltaMonoT         = conf["DeltaMonoT"].as<double>();
    if (conf["DOmega"])         DOmega             = conf["DOmega"].as<double>();
    if (conf["dtom"])           dtom               = conf["dtom"].as<double>();
    if (conf["T0"])             T0                 = conf["T0"].as<double>();
    if (conf["Fcorot"])         Fcorot             = conf["Fcorot"].as<double>();
    if (conf["omega"])          omega0             = conf["omega"].as<double>();
    if (conf["fixed"])          fixed              = conf["fixed"].as<bool>();
    if (conf["self"])           fixed              = conf["self"].as<bool>();
    if (conf["oscil"])          oscil              = conf["oscil"].as<bool>();
    if (conf["Ofreq"])          Ofreq              = conf["Ofreq"].as<double>();
    if (conf["Oamp"])           Oamp               = conf["Oamp"].as<double>();
    if (conf["alpha"])          alpha              = conf["alpha"].as<double>();
    if (conf["x0"])             bps[0]             = conf["x0"].as<double>();
    if (conf["y0"])             bps[1]             = conf["y0"].as<double>();
    if (conf["z0"])             bps[2]             = conf["z0"].as<double>();
    if (conf["u0"])             vel[0]             = conf["u0"].as<double>();
    if (conf["v0"])             vel[1]             = conf["v0"].as<double>();
    if (conf["w0"])             vel[2]             = conf["w0"].as<double>();
    if (conf["monopole"])       monopole           = conf["monopole"].as<bool>();
    if (conf["follow"])         monopole_follow    = conf["follow"].as<bool>();
    if (conf["onoff"])          monopole_onoff     = conf["onoff"].as<bool>();
    if (conf["monofrac"])       monopole_frac      = conf["monofrac"].as<double>();
    if (conf["quadfrac"])       quadrupole_frac    = conf["quadfrac"].as<double>();
    if (conf["filename"])       filename           = conf["filename"].as<string>();
    


    if (conf["ctrname"])        ctr_name           = conf["ctrname"].as<string>();
    if (conf["angmname"])       angm_name          = conf["angmname"].as<string>();
    if (conf["tblname"])        table_name         = conf["tblname"].as<string>();
    if (conf["length"])         length             = conf["length"].as<double>();
    if (conf["bratio"])         bratio             = conf["bratio"].as<double>();
    if (conf["cratio"])         cratio             = conf["cratio"].as<double>();
    if (conf["amp"])            amplitude          = conf["amp"].as<double>();
    if (conf["angmomfac"])      angmomfac          = conf["angmomfac"].as<double>();
    if (conf["barmass"])        barmass            = conf["barmass"].as<double>();
    if (conf["Ton"])            Ton                = conf["Ton"].as<double>();
    if (conf["Toff"])           Toff               = conf["Toff"].as<double>();
    if (conf["TmonoOn"])        TmonoOn            = conf["TmonoOn"].as<double>();
    if (conf["TmonoOff"])       TmonoOff           = conf["TmonoOff"].as<double>();
    if (conf["DeltaT"])         DeltaT             = conf["DeltaT"].as<double>();
    if (conf["DeltaMonoT"])     DeltaMonoT         = conf["DeltaMonoT"].as<double>();
    if (conf["DOmega"])         DOmega             = conf["DOmega"].as<double>();
    if (conf["dtom"])           dtom               = conf["dtom"].as<double>();
    if (conf["T0"])             T0                 = conf["T0"].as<double>();
    if (conf["Fcorot"])         Fcorot             = conf["Fcorot"].as<double>();
    if (conf["omega"])          omega0             = conf["omega"].as<double>();
    if (conf["fixed"])          fixed              = conf["fixed"].as<bool>();
    if (conf["self"])           fixed              = conf["self"].as<bool>();
    if (conf["oscil"])          oscil              = conf["oscil"].as<bool>();
    if (conf["Ofreq"])          Ofreq              = conf["Ofreq"].as<double>();
    if (conf["Oamp"])           Oamp               = conf["Oamp"].as<double>();
    if (conf["alpha"])          alpha              = conf["alpha"].as<double>();
    if (conf["x0"])             bps[0]             = conf["x0"].as<double>();
    if (conf["y0"])             bps[1]             = conf["y0"].as<double>();
    if (conf["z0"])             bps[2]             = conf["z0"].as<double>();
    if (conf["u0"])             vel[0]             = conf["u0"].as<double>();
    
    
    if (conf["v0"])             vel[1]             = conf["v0"].as<double>();
    if (conf["w0"])             vel[2]             = conf["w0"].as<double>();
    if (conf["monopole"])       monopole           = conf["monopole"].as<bool>();
    if (conf["follow"])         monopole_follow    = conf["follow"].as<bool>();
    if (conf["onoff"])          monopole_onoff     = conf["onoff"].as<bool>();
    if (conf["monofrac"])       monopole_frac      = conf["monofrac"].as<double>();
    if (conf["quadfrac"])       quadrupole_frac    = conf["quadfrac"].as<double>();
    if (conf["filename"])       filename           = conf["filename"].as<string>();
  }
  catch (YAML::Exception & error) {
    if (myid==0) std::cout << "Error parsing parameters in UserEBar: "
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


void UserEBar::determine_acceleration_and_potential(void)
{
  if (timing) timer_tot.start();
				// Write to bar state file, if true
  bool update = false;

  if (c1) c1->get_angmom();	// Tell component to compute angular momentum
  // cout << "Process " << myid << ": Lz=" << c1->angmom[2] << endl; // debug

  if (firstime) {
    
    ellip = new EllipForce(length, length*bratio, length*bratio*cratio,
			   barmass, 200, 200);

    if (omega0 < 0.0) {

      double R=length*Fcorot;
      double phi, theta=0.5*M_PI;
      double dens0, potl0, dens, potl, potr, pott, potp;
      double avg=0.0;
      
      for (int n=0; n<8; n++) {
	phi = 2.0*M_PI/8.0 * n;

	for (auto c : comp->components) {
	
	  if (c->force->geometry == PotAccel::sphere || 
	      c->force->geometry == PotAccel::cylinder) {
	  
	    ((Basis*)c->force)->
	      determine_fields_at_point_sph(R, theta, phi,
					    &dens0, &potl0,
					    &dens, &potl, &potr, &pott, &potp);
	  
	    avg += potr/8.0;
	  }
	}
      }

      omega0 = sqrt(avg/R);
    }

    if (dtom>0.0)
      omega = omega0*(1.0 + DOmega*0.5*(1.0 + erf( (tnow - T0)/dtom )));
    else
      omega = omega0*(1.0 + DOmega*(tnow - T0*0.5));

    const int N = 100;
    LegeQuad gq(N);

    double a1 = length;
    double a2 = bratio*a1;
    double a3 = cratio*a2;

    double geom = pow(a1*a2*a3, 1.0/3.0);

    double A12 = a1*a1/geom/geom;
    double A22 = a2*a2/geom/geom;
    double A32 = a3*a3/geom/geom;

    double u, d, t, denom, ans1=0.0, ans2=0.0;
    double mass = barmass * fabs(amplitude);

    for (int i=0; i<N; i++) {
      t = 0.5*M_PI*gq.knot(i);
      u = tan(t);
      d = cos(t);
      d = 1.0/(d*d);

      denom = sqrt( (A12+u)*(A22+u)*(A32+u) );
      ans1 += d*gq.weight(i) /( (A12+u)*denom );
      ans2 += d*gq.weight(i) /( (A22+u)*denom );
    }
    ans1 *= 0.5*M_PI;
    ans2 *= 0.5*M_PI;

    if (myid==0) {

      cout << "====================================================\n";
      cout << "Computed quadrupole fit to homogenous ellipsoid\n";
      cout << "with Mass=" << mass << " A_1=" << a1 << " A_2=" << a2 
	   << " A_3=" << a3 << "\n"
	   << "with an exact fit to asymptotic quadrupole, e.g.\n"
	   << "     U_{22} = b1 r**2/( 1+(r/b5)**5 ) or\n"
	   << "            = b1 r**2/( 1+ r/b5 )**5\n";
      cout << "====================================================\n";

      cout << "V_1=" << ans1 << endl;
      cout << "V_2=" << ans2 << endl;
      cout << "I_3=" << 0.2*mass*(a1*a1 + a2*a2) << endl;
      cout << "Omega(0)=" << omega0 << endl;

    }

    double rho = mass/(4.0*M_PI/3.0*a1*a2*a3);
    double b1 = M_PI*rho*sqrt(2.0*M_PI/15.0)*(ans1 - ans2);
    double b25 = 0.4*a1*a2*a3*(a2*a2 - a1*a1)/(ans1 - ans2);

    b5 = pow(b25, 0.2);
    afac = 2.0 * b1;
    // afac = b1;

    if (myid==0) {
      cout << "b1=" << b1 << endl;
      cout << "b5=" << b5 << endl;
      cout << "afac=" << afac << endl;
      cout << "====================================================\n" 
	   << flush;

      if (!restart) {
	ofstream out(name.c_str(), ios::out | ios::app);

	out << setw(15) << "# Time"
	    << setw(15) << "Phi"
	    << setw(15) << "Omega"
	    << setw(15) << "L_z(Bar)"
	    << setw(15) << "L_z(PS)"
	    << setw(15) << "Amp"
	    << setw(15) << "x"
	    << setw(15) << "y"
	    << setw(15) << "z"
	    << setw(15) << "u"
	    << setw(15) << "v"
	    << setw(15) << "w"
	    << setw(15) << "ax"
	    << setw(15) << "ay"
	    << setw(15) << "az"
	    << endl;
      }
    }

    Iz = 0.2*mass*(a1*a1 + a2*a2) * angmomfac;
    Lz = Iz * omega;

    if (c1)
      Lz0 = c1->angmom[2];
    else
      Lz0 = 0.0;

    posang = 0.0;
    lastomega = omega;
    lasttime = tnow;
    
    if (restart) {

      if (myid == 0) {
	
	// Backup up old file
	string backupfile = name + ".bak";
	string command("cp ");
	command += name + " " + backupfile;
	if (system(command.c_str()) == -1) {
	  std::cerr << "UserEBar: error in executing <"
		    << command << ">" << endl;
	}
	
	// Open new output stream for writing
	ofstream out(name.c_str());
	if (!out) {
	  throw FileCreateError(name, "UserEBar: error opening new log file",
				__FILE__, __LINE__);
	}
	
	// Open old file for reading
	ifstream in(backupfile.c_str());
	if (!in) {
	  throw FileOpenError(backupfile, "UserEBar: error opening original log file",
				__FILE__, __LINE__);
	}

	const int linesize = 1024;
	char line[linesize];
	
	in.getline(line, linesize); // Discard header
	in.getline(line, linesize); // Next line

	double Lzp,  am1;
	bool firstime1 = true;
	while (in) {
	  istringstream ins(line);

	  lastomega = omega;

	  ins >> lasttime;
	  ins >> posang;
	  ins >> omega;
	  ins >> Lz;
	  ins >> Lzp;
	  ins >> am1;
	  ins >> bps[0];
	  ins >> bps[1];
	  ins >> bps[2];
	  ins >> vel[0];
	  ins >> vel[1];
	  ins >> vel[2];
	  ins >> acc[0];
	  ins >> acc[1];
	  ins >> acc[2];

	  if (firstime1) {
	    Lz0 = Lzp;
	    firstime1 = false;
	  }

	  if (lasttime >= tnow) break;

	  out << line << "\n";

	  in.getline(line, linesize); // Next line
	}

	cout << "UserEBar: restart at T=" << lasttime 
	     << " with PosAng=" << posang
	     << ", Omega=" << omega
	     << ", Lz=" << Lz
	     << ", Lz0=" << Lz0
	     << endl;

      }

      MPI_Bcast(&lasttime, 1, MPI_DOUBLE, 0, MPI_COMM_WORLD);
      MPI_Bcast(&posang, 1, MPI_DOUBLE, 0, MPI_COMM_WORLD);
      MPI_Bcast(&omega, 1, MPI_DOUBLE, 0, MPI_COMM_WORLD);
      MPI_Bcast(&lastomega, 1, MPI_DOUBLE, 0, MPI_COMM_WORLD);
      MPI_Bcast(&Lz, 1, MPI_DOUBLE, 0, MPI_COMM_WORLD);
      MPI_Bcast(&Lz0, 1, MPI_DOUBLE, 0, MPI_COMM_WORLD);
      MPI_Bcast(&bps[0], 3, MPI_DOUBLE, 0, MPI_COMM_WORLD);
      MPI_Bcast(&vel[0], 3, MPI_DOUBLE, 0, MPI_COMM_WORLD);
      MPI_Bcast(&acc[0], 3, MPI_DOUBLE, 0, MPI_COMM_WORLD);

      // Recompute Lz from log output
      if (c1) Lz = Lz - Lz0 + c1->angmom[2];
    }

    firstime = false;
    if (mstep <= mupdate) update = true;

  } else {

    if (!fixed) {
      if (c1)
	omega = (Lz + Lz0 - c1->angmom[2])/Iz;
      else
	omega = Lz/Iz;
    }
    else {
      if (dtom>0.0)
	omega = omega0*(1.0 + DOmega*0.5*(1.0 + erf( (tnow - T0)/dtom )));
      else
	omega = omega0*(1.0 + DOmega*(tnow - T0*0.5));
    }
    
    if ( mlevel==0 and fabs(tnow-lasttime) > 2.0*std::numeric_limits<double>::min()) {
      posang += 0.5*(omega + lastomega)*(tnow - lasttime);
      lastomega = omega;
      lasttime = tnow;
      update = true;
    }
  }

				// Zero thread variables
  for (int n=0; n<nthrds; n++) {
    for (int k=0; k<3; k++) tacc[n][k] = 0.0;
  }

  if (timing) timer_thrd.start();
  exp_thread_fork(false);
  if (timing) timer_thrd.stop();

				// Get full contribution from all threads
  for (int k=0; k<3; k++) acc[k] = acc1[k] = 0.0;
  for (int n=0; n<nthrds; n++) {
    for (int k=0; k<3; k++) acc1[k] += tacc[n][k]/barmass;
  }

				// Get contribution from all processes
  MPI_Allreduce(&acc1[0], &acc[0], 3, MPI_DOUBLE, MPI_SUM, MPI_COMM_WORLD);

				// Backward Euler
  if (mlevel== 0 && monopole && monopole_follow) {
    for (int k=0; k<3; k++) {
      bps[k] += vel[k] * (tnow - teval[mlevel]);
      vel[k] += acc[k] * (tnow - teval[mlevel]);
    }
    for (unsigned m=mlevel; m<=multistep; m++) teval[m] = tnow;
  }

  if (myid==0 && update) 
    {
      ofstream out(name.c_str(), ios::out | ios::app);
      out.setf(ios::scientific);

      out << setw(15) << tnow
	  << setw(15) << posang
	  << setw(15) << omega;

      if (c1)
	out << setw(15) << Lz + Lz0 - c1->angmom[2]
	    << setw(15) << c1->angmom[2];
      else
	out << setw(15) << Lz
	    << setw(15) << 0.0;

      if (amplitude==0.0)
	out << setw(15) <<  0.0;
      else
	out << setw(15) << amplitude/fabs(amplitude) *  
	  0.5*(1.0 + erf( (tnow - Ton )/DeltaT )) *
	  0.5*(1.0 - erf( (tnow - Toff)/DeltaT ));

      for (int k=0; k<3; k++) out << setw(15) << bps[k];
      for (int k=0; k<3; k++) out << setw(15) << vel[k];
      for (int k=0; k<3; k++) out << setw(15) << acc[k];
      
      out << endl;
    }

  if (timing) {
    timer_tot.stop();
    cout << setw(20) << "Bar total: "
	 << setw(18) << timer_tot.getTime() << endl
	 << setw(20) << "Bar threads: "
	 << setw(18) << timer_thrd.getTime() << endl;
    timer_tot.reset();
    timer_thrd.reset();
  }

  print_timings("UserEBar: acceleration timings");
}


void * UserEBar::determine_acceleration_and_potential_thread(void * arg) 
{
  int id = *((int*)arg), nbodies, nbeg, nend, indx;
  double fac, ffac, amp = 0.0;
  double xx, yy, zz, rr, nn, pp, extpot, M0=0.0;
  vector<double> pos(3), acct(3); 
  double cos2p = cos(2.0*posang);
  double sin2p = sin(2.0*posang);

  thread_timing_beg(id);

  double fraction_on =   0.5*(1.0 + erf( (tnow - Ton )/DeltaT )) ;
  double fraction_off =  0.5*(1.0 - erf( (tnow - Toff)/DeltaT )) ;

  double quad_onoff = 
    fraction_on*( (1.0 - quadrupole_frac) + quadrupole_frac * fraction_off );

  double mono_fraction = 
    0.5*(1.0 + erf( (tnow - TmonoOn )/DeltaMonoT )) *
    0.5*(1.0 - erf( (tnow - TmonoOff)/DeltaMonoT )) ;

  double mono_onoff = 
    (1.0 - monopole_frac) + monopole_frac*mono_fraction;

  if (table) {
    if (tnow<timeq[0]) {
      afac = ampq[0];
      b5 = b5q[0];
    } else if (tnow>timeq[qlast]) {
      afac = ampq[qlast];
      b5 = b5q[qlast];
    } else {
      afac = odd2(tnow, timeq, ampq, 0);
      b5 = odd2(tnow, timeq, b5q, 0);
    }
  }

  if (amplitude==0.0) 
    amp = 0.0;
  else
    amp = afac * amplitude/fabs(amplitude) * quad_onoff;

  if (oscil)
    amp *= (1.0 + Oamp*sin(Ofreq*(tnow - Ton)))/(1.0 + fabs(Oamp));

#if HAVE_LIBCUDA==1		// Cuda compatibility
  getParticlesCuda(cC);
#endif

  for (unsigned lev=mlevel; lev<=multistep; lev++) {

    nbodies = cC->levlist[lev].size();
    nbeg = nbodies*(id  )/nthrds;
    nend = nbodies*(id+1)/nthrds;

    for (int i=nbeg; i<nend; i++) {

      indx = cC->levlist[lev][i];
      
      for (int k=0; k<3; k++) pos[k] = cC->Pos(indx, k);

      if (c0)
	for (int k=0; k<3; k++) pos[k] -= c0->center[k];
      else if (monopole)
	for (int k=0; k<3; k++) pos[k] -= bps[k];
    
      xx = pos[0];
      yy = pos[1];
      zz = pos[2];
      rr = sqrt( xx*xx + yy*yy + zz*zz );

				// Variable sharpness potential
      fac = 1.0 + pow(rr/b5, alpha);
      ffac = -amp*numfac/pow(fac, 5.0/alpha+1.0);
      pp = (xx*xx - yy*yy)*cos2p + 2.0*xx*yy*sin2p;
      nn = pp * pow(rr/b5, alpha-1.0) / ( b5*rr );

				// Quadrupole acceleration
      acct[0] = ffac*
	( 2.0*( xx*cos2p + yy*sin2p)*fac - 5.0*nn*xx );
    
      acct[1] = ffac*
	( 2.0*(-yy*cos2p + xx*sin2p)*fac - 5.0*nn*yy );

      acct[2] = ffac*
	( -5.0*nn*zz );
    
				// Quadrupole potential
      extpot = -ffac*pp*fac;
    

				// Monopole contribution
      if (monopole) {

	M0 = ellip->getMass(rr);
	
	if (monopole_onoff) M0 *= mono_onoff;
	
	for (int k=0; k<3; k++) {
				// Add monopole acceleration
	  acct[k] += -M0*pos[k]/(rr*rr*rr);

				// Force on bar (via Newton's 3rd law)
	  tacc[id][k] += -cC->Mass(indx) * acct[k];
	}

				// Monopole potential
	extpot += ellip->getPot(rr);
      }

				// Add bar acceleration to particle
      cC->AddAcc(indx, acct);
    
				// Add external potential
      cC->AddPotExt(indx, extpot);

    }
  }

  thread_timing_end(id);

  return (NULL);
}


extern "C" {
  ExternalForce *makerEBar(const YAML::Node& conf)
  {
    return new UserEBar(conf);
  }
}

class proxyebar { 
public:
  proxyebar()
  {
    factory["userebar"] = makerEBar;
  }
};

static proxyebar p;
