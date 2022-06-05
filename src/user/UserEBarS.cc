#include <math.h>
#include <sstream>

#include "expand.H"
#include <localmpi.H>
#include <UserEBarS.H>

UserEBarS::UserEBarS(const YAML::Node& conf) : ExternalForce(conf)
{
  id = "RotatingBarWithMonopoleTorque";

  length = 1.0;			// Bar length
  bratio = 0.5;			// Ratio of b to a
  cratio = 0.1;			// Ratio of c to b
  amplitude = 0.3;		// Bar quadrupole amplitude
  angmomfac = 1.0;		// Artifically change the total bar ang mom
  barmass = 1.0;		// Total bar mass
  Ton = -20.0;			// Turn on start time
  Toff = 200.0;			// Turn off start time
  TmonoOn = -20.0;		// Turn on start time for monopole
  TmonoOff = 200.0;		// Turn off start time monopole
  DeltaT = 1.0;			// Turn on duration
  DeltaMonoT = 1.0;		// Turn on duration for monopole
  DOmega = 0.0;			// Change in pattern speed
  tom0 = 1.0;			// Midpoint of forced bar slow down
  dtom = -1.0;			// Width of forced bar slow down
  T0 = 10.0;			// Center of pattern speed change
  Fcorot  = 1.0;		// Corotation factor
  fixed = false;		// Constant pattern speed
  soft = false;			// Use soft form of the bar potential
  table = false;		// Not using tabled quadrupole
  monopole = true;		// Use the monopole part of the potential
  monopole_onoff = false;	// To apply turn-on and turn-off to monopole
  monopole_follow = true;	// Follow monopole center
  monopole_frac = 1.0;		// Fraction of monopole to turn off
  quadrupole_frac = 1.0;	// Fraction of quadrupole to turn off

				// Output file name
  filename = outdir + "BarRot." + runtag;

  firstime = true;
  omega0 = -1.0;

  ctr_name = "";		// Default component for com
  table_name = "";		// Default for input b1,b5 table
  
  ellip = 0;

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
      cerr << "Process " << myid << ": can't find desired component <"
	   << ctr_name << " for centering>" << endl;
      MPI_Abort(MPI_COMM_WORLD, 35);
    }

  }
  else
    c0 = NULL;


  if (table_name.size()>0) {
				// Read in data
    ifstream in(string(outdir+table_name).c_str());
    if (!in) {
      cerr << "Process " << myid << ": error opening quadrupole file <"
	   << outdir + table_name << ">" << endl;
      MPI_Abort(MPI_COMM_WORLD, 35);
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
  tevalS = tnow;
  for (int k=0; k<3; k++) bps[k] = vel[k] = acc[k] = 0.0;

  // Assign working vectors for each thread
  okM  = 0;
  TzM.resize(comp->ncomp);
  TzN.resize(nthrds);
  tacc.resize(nthrds);
  for (auto & v : tacc) v.resize(3);
  
  userinfo();
}

UserEBarS::~UserEBarS()
{
  // Nothing
}

void UserEBarS::userinfo()
{
  if (myid) return;		// Return if node master node

  print_divider();

  cout << "** User routine ROTATING BAR with MONOPOLE initialized (self-consistent torque version), " ;
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
    if (monopole_follow)
      cout << "self-consistent monopole centering, ";
    else
      cout << "monopole center fixed, ";
  }
  else
    cout << "without monopole, ";
  if (soft)
    cout << "soft potential, ";
  else
    cout << "standard potential, ";
  if (c0) 
    cout << "center on component <" << ctr_name << ">";
  else
    cout << "using inertial center";

  if (table)
    cout << ", using user quadrupole table";

  cout << endl;
  print_divider();
}

void UserEBarS::initialize()
{
  try {
    if (conf["ctrname"])        ctr_name           = conf["ctrname"].as<string>();
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
    if (conf["tom0"])           tom0               = conf["tom0"].as<double>();
    if (conf["dtom"])           dtom               = conf["dtom"].as<double>();
    if (conf["T0"])             T0                 = conf["T0"].as<double>();
    if (conf["Fcorot"])         Fcorot             = conf["Fcorot"].as<double>();
    if (conf["omega"])          omega0             = conf["omega"].as<double>();
    if (conf["fixed"])          fixed              = conf["fixed"].as<bool>();
    if (conf["self"])           fixed              = conf["self"].as<bool>();
    if (conf["soft"])           soft               = conf["soft"].as<bool>();
    if (conf["monopole"])       monopole           = conf["monopole"].as<bool>();
    if (conf["follow"])         monopole_follow    = conf["follow"].as<bool>();
    if (conf["onoff"])          monopole_onoff     = conf["onoff"].as<bool>();
    if (conf["monofrac"])       monopole_frac      = conf["monofrac"].as<double>();
    if (conf["quadfrac"])       quadrupole_frac    = conf["quadfrac"].as<double>();
    if (conf["filename"])       filename           = conf["filename"].as<string>();
  }
  catch (YAML::Exception & error) {
    if (myid==0) std::cout << "Error parsing parameters in UserEBarS: "
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


void UserEBarS::determine_acceleration_and_potential(void)
{

  if (firstime) {
    
    ellip = std::make_shared<EllipForce>
      (length, length*bratio, length*bratio*cratio, barmass, 200, 200);

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
	      determine_fields_at_point_sph(R, theta, phi, &dens0, &potl0,
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
    // afac = 2.0 * b1;
    afac = b1;

    if (myid==0) {
      cout << "b1=" << b1 << endl;
      cout << "b5=" << b5 << endl;
      cout << "afac=" << afac << endl;
      cout << "====================================================\n" 
	   << flush;

      name = outdir + filename;
      name += ".barstat";

      if (!restart) {
	ofstream out(name.c_str(), ios::out | ios::app);

	out << setw(15) << "# Time"
	    << setw(15) << "Phi"
	    << setw(15) << "Omega"
	    << setw(15) << "L_z(Bar)";

	for (auto c : comp->components) {
	  ostringstream ostr;
	  ostr << "T_z(" << c->name << ")";
	  out << setw(15) << ostr.str().substr(0,15);
	}
	
	out << setw(15) << "T_z(total)"
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

    posang = 0.0;
    lastomega = omega;
    lasttime = tnow;
    
    if (restart) {

      if (myid == 0) {
	
	// Backup up old file
	string backupfile = outdir + name + ".bak";
	string command("cp ");
	command += outdir + name + " " + backupfile;
	if (system(command.c_str()) == -1) {
	  std::cerr << "UserEBarS: error in executing <"
		    << command << ">" << endl;
	}

	// Open new output stream for writing
	ofstream out(string(outdir+name).c_str());
	if (!out) {
	  throw FileCreateError(outdir+name, "UserEBarS: error opening new log file",
				__FILE__, __LINE__);
	}
	
	// Open old file for reading
	ifstream in(backupfile.c_str());
	if (!in) {
	  throw FileOpenError(backupfile, "UserEBarS: error opening original log file",
				__FILE__, __LINE__);
	}

	const int linesize = 1024;
	char line[linesize];
	
	in.getline(line, linesize); // Discard header
	in.getline(line, linesize); // Next line

	double am1;
	bool firstime1 = true;
	while (in) {
	  istringstream ins(line);

	  lastomega = omega;

	  ins >> lasttime;
	  ins >> posang;
	  ins >> omega;
	  ins >> Lz;
	  ins >> Tz;
	  for (int m=0; m<comp->ncomp; m++) ins >> TzM[m];
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
	    firstime1 = false;
	  }

	  if (lasttime >= tnow) break;

	  out << line << "\n";

	  in.getline(line, linesize); // Next line
	}

	cout << "UserEBarS: restart at T=" << lasttime 
	     << " with PosAng=" << posang
	     << ", Omega=" << omega
	     << ", Lz=" << Lz
	     << endl;

      }

      MPI_Bcast(&lasttime, 1, MPI_DOUBLE, 0, MPI_COMM_WORLD);
      MPI_Bcast(&posang, 1, MPI_DOUBLE, 0, MPI_COMM_WORLD);
      MPI_Bcast(&omega, 1, MPI_DOUBLE, 0, MPI_COMM_WORLD);
      MPI_Bcast(&lastomega, 1, MPI_DOUBLE, 0, MPI_COMM_WORLD);
      MPI_Bcast(&Lz, 1, MPI_DOUBLE, 0, MPI_COMM_WORLD);
      MPI_Bcast(&bps[0], 3, MPI_DOUBLE, 0, MPI_COMM_WORLD);
      MPI_Bcast(&vel[0], 3, MPI_DOUBLE, 0, MPI_COMM_WORLD);
      MPI_Bcast(&acc[0], 3, MPI_DOUBLE, 0, MPI_COMM_WORLD);
    }

    firstime = false;

  } else if (mlevel==0) {

    if (!fixed) {
      omega = Lz/Iz;
    }
    else {
      if (dtom>0.0)
	omega = omega0*(1.0 + DOmega*0.5*(1.0 + erf( (tnow - T0)/dtom )));
      else
	omega = omega0*(1.0 + DOmega*(tnow - T0*0.5));
    }
    
    if ( fabs(tnow-lasttime) > 2.0*DBL_EPSILON) {
      posang += 0.5*(omega + lastomega)*(tnow - lasttime);
      lastomega = omega;
      lasttime = tnow;
    }
  }

				// Zero out thread variables
  for (int n=0; n<nthrds; n++) {
    TzN[n] = 0.0;
    for (int k=0; k<3; k++) tacc[n][k] = 0.0;
  }

				// Identify list #
  cid=0;
  for (auto c : comp->components) {
    if (cC == c) break;
    cid++;
  }

  if (cid==comp->ncomp) {
    cerr << "Process " << myid << ": failed to identify component!\n";
    MPI_Abort(MPI_COMM_WORLD, 135);
  }

  okM++;

#if HAVE_LIBCUDA==1		// Cuda compatibility
  getParticlesCuda(cC);
#endif

  exp_thread_fork(false);

				// Get full contribution from all threads
  for (int k=0; k<3; k++) acc[k] = acc1[k] = 0.0;
  for (int n=0; n<nthrds; n++) {
    for (int k=0; k<3; k++) acc1[k] += tacc[n][k]/barmass;
  }
				// Get contribution from all processes
  MPI_Allreduce(&acc1[0], &acc[0], 3, MPI_DOUBLE, MPI_SUM, MPI_COMM_WORLD);

				// Get torque for each component 
				// from all threads
  double Tz1=0.0, Tz0=0.0;
  for (int n=0; n<nthrds; n++) Tz1 -= TzN[n];

				// Get contribution from all processes
  MPI_Allreduce(&Tz1, &Tz0, 1, MPI_DOUBLE, MPI_SUM, MPI_COMM_WORLD);

  TzM[cid] += Tz0;
  
				// Backward Euler
  if (monopole && monopole_follow) {
    for (int k=0; k<3; k++) {
      bps[k] += vel[k] * (tnow - teval[mlevel]);
      vel[k] += acc[k] * (tnow - teval[mlevel]);
    }
  }
  
				// All components computed
  if (okM == comp->ncomp) {

				// Total torque
    Tz0 = 0.0;
    for (int m=0; m<comp->ncomp; m++) Tz0 += TzM[m];

				// Update bar angular momentum
    Lz += Tz0 * (tnow - teval[mlevel]);

    okM   = 0;
    for (unsigned m=mlevel; m<=multistep; m++) teval[mlevel] = tnow;

    if (myid==0) 
      {
	ofstream out(string(outdir+name).c_str(), ios::out | ios::app);
	out.setf(ios::scientific);

	out << setw(15) << tnow
	    << setw(15) << posang
	    << setw(15) << omega
	    << setw(15) << Lz;

	for (int m=0; m<comp->ncomp; m++) 
	  out << setw(15) << TzM[m];
	
	out << setw(15) << Tz0;

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
  }

}


void * UserEBarS::determine_acceleration_and_potential_thread(void * arg) 
{
  int nbodies = cC->Number();
  int id = *((int*)arg);
  int nbeg = nbodies*id/nthrds;
  int nend = nbodies*(id+1)/nthrds;

  double fac, ffac, amp = 0.0;
  double xx, yy, zz, rr, nn, pp, extpot, M0=0.0;
  vector<double> pos(3), acct(3); 
  double cos2p = cos(2.0*posang);
  double sin2p = sin(2.0*posang);

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

  for (int i=nbeg; i<nend; i++) {

				// If we are multistepping, compute accel 
				// only at or below this level

    if (multistep && (cC->Part(i)->level < mlevel)) continue;

    for (int k=0; k<3; k++) pos[k] = cC->Pos(i, k);

    if (c0)
      for (int k=0; k<3; k++) pos[k] -= c0->center[k];
    else if (monopole)
      for (int k=0; k<3; k++) pos[k] -= bps[k];
    
    xx = pos[0];
    yy = pos[1];
    zz = pos[2];
    rr = sqrt( xx*xx + yy*yy + zz*zz );

    if (soft) {
      fac = 1.0 + rr/b5;

      ffac = -amp*numfac/pow(fac, 6.0);

      pp = (xx*xx - yy*yy)*cos2p + 2.0*xx*yy*sin2p;
      nn = pp /( b5*rr ) ;
    } else {
      fac = 1.0 + pow(rr/b5, 5.0);

      ffac = -amp*numfac/(fac*fac);
      
      pp = (xx*xx - yy*yy)*cos2p + 2.0*xx*yy*sin2p;
      nn = pp * pow(rr/b5, 3.0)/(b5*b5);
    }

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
	tacc[id][k] += -cC->Mass(i) * acct[k];
      }

				// Monopole potential
      extpot += ellip->getPot(rr);
    }

				// Add bar acceleration to particle
    cC->AddAcc(i, acct);
    
				// Add external potential
    cC->AddPotExt(i, extpot);

				// Accumulate torque
    
    TzN[id] += cC->Mass(i) * (pos[0]*acct[1] - pos[1]*acct[0]);
  }

  return (NULL);
}


extern "C" {
  ExternalForce *makerEBarS(const YAML::Node& conf)
  {
    return new UserEBarS(conf);
  }
}

class proxyebars { 
public:
  proxyebars()
  {
    factory["userebars"] = makerEBarS;
  }
};

static proxyebars p;
