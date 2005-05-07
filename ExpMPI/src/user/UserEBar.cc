#include <math.h>
#include <sstream>

#include "expand.h"
#include <localmpi.h>
#include <UserEBar.H>

UserEBar::UserEBar(string &line) : ExternalForce(line)
{
  id = "RotatingBarWithMonopole";

  length = 1.0;			// Bar length
  bratio = 0.5;			// Ratio of b to a
  cratio = 0.1;			// Ratio of c to b
  amplitude = 0.3;		// Bar quadrupole amplitude
  barmass = 1.0;		// Total bar mass
  Ton = -20.0;			// Turn on start time
  Toff = 200.0;			// Turn off start time
  DeltaT = 1.0;			// Turn on duration
  DOmega = 0.0;			// Change in pattern speed
  tom0 = 1.0;			// Midpoint of forced bar slow down
  dtom = -1.0;			// Width of forced bar slow down
  T0 = 10.0;			// Center of pattern speed change
  Fcorot  = 1.0;		// Corotation factor
  fixed = false;		// Constant pattern speed
  soft = false;			// Use soft form of the bar potential
  table = false;		// Not using tabled quadrupole
  monopole = true;		// Use the monopole part of the potential
				// Output file name
  filename = outdir + "BarRot." + runtag;

  firstime = true;
  omega0 = -1.0;

  ctr_name = "";		// Default component for com
  angm_name = "";		// Default component for angular momentum
  table_name = "";		// Default for input b1,b5 table
  
  initialize();

  if (ctr_name.size()>0) {
				// Look for the fiducial component for
				// centering
    bool found = false;
    list<Component*>::iterator cc;
    Component *c;
    for (cc=comp.components.begin(); cc != comp.components.end(); cc++) {
      c = *cc;
      if ( !ctr_name.compare(c->name) ) {
	c0 = c;
	found = true;
      break;
      }
    }

    if (!found) {
      cerr << "Process " << myid << ": can't find desired component <"
	   << ctr_name << ">" << endl;
      MPI_Abort(MPI_COMM_WORLD, 35);
    }

  }
  else
    c0 = NULL;

  if (angm_name.size()>0) {
				// Look for the fiducial component
				// for angular momentum
    bool found = false;
    list<Component*>::iterator cc;
    Component *c;
    for (cc=comp.components.begin(); cc != comp.components.end(); cc++) {
      c = *cc;
      if ( !angm_name.compare(c->name) ) {
	c1 = c;
	found = true;
      break;
      }
    }

    if (!found) {
      cerr << "Process " << myid << ": can't find desired component <"
	   << angm_name << ">" << endl;
      MPI_Abort(MPI_COMM_WORLD, 35);
    }

  }
  else
    c1 = NULL;


  if (table_name.size()>0) {
				// Read in data
    ifstream in(table_name.c_str());
    if (!in) {
      cerr << "Process " << myid << ": error opening quadrupole file <"
	   << table_name << ">" << endl;
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
  teval = tnow;
  for (int k=0; k<3; k++) bps[k] = vel[k] = acc[k] = 0.0;

  // Assign working vectors for each thread
  tacc = new double* [nthrds];
  for (int n=0; n<nthrds; n++) tacc[n] = new double [3];

  userinfo();
}

UserEBar::~UserEBar()
{
  for (int n=0; n<nthrds; n++) delete [] tacc[n];
  delete [] tacc;
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
  if (omega0<0.0)
    cout << "initial pattern speed to be computed, ";
  else
    cout << "initial pattern speed " << omega0 << ", ";
  if (monopole)
    cout << "using monopole, ";
  else
    cout << "without monopole, ";
  if (soft)
    cout << "soft potential, ";
  else
    cout << "standard potential, ";
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
  string val;

  if (get_value("ctrname", val))	ctr_name = val;
  if (get_value("angmname", val))	angm_name = val;
  if (get_value("tblname", val))	table_name = val;
  if (get_value("length", val))		length = atof(val.c_str());
  if (get_value("bratio", val))		bratio = atof(val.c_str());
  if (get_value("cratio", val))		cratio = atof(val.c_str());
  if (get_value("amp", val))		amplitude = atof(val.c_str());
  if (get_value("barmass", val))	barmass = atof(val.c_str());
  if (get_value("Ton", val))		Ton = atof(val.c_str());
  if (get_value("Toff", val))		Toff = atof(val.c_str());
  if (get_value("DeltaT", val))		DeltaT = atof(val.c_str());
  if (get_value("DOmega", val))		DOmega = atof(val.c_str());
  if (get_value("tom0", val))     	tom0 = atof(val.c_str());
  if (get_value("dtom", val))     	dtom = atof(val.c_str());
  if (get_value("T0", val))		T0 = atof(val.c_str());
  if (get_value("Fcorot", val))		Fcorot = atof(val.c_str());
  if (get_value("omega", val))		omega0 = atof(val.c_str());
  if (get_value("fixed", val))		fixed = atoi(val.c_str()) ? true:false;
  if (get_value("self", val))		fixed = atoi(val.c_str()) ? false:true;
  if (get_value("soft", val))		soft = atoi(val.c_str()) ? true:false;
  if (get_value("monopole", val))	monopole = atoi(val.c_str()) ? true:false;
  if (get_value("filename", val))	filename = val;

}


void UserEBar::determine_acceleration_and_potential(void)
{
				// Write to bar state file, if true
  bool update = false;

  if (c1) c1->get_angmom();	// Tell component to compute angular momentum
  // cout << "Process " << myid << ": Lz=" << c1->angmom[2] << endl; // debug

  if (firstime) {
    
    ellip = new EllipForce(length, length*bratio, length*bratio*cratio,
			   barmass, 200, 200);

    list<Component*>::iterator cc;
    Component *c;

    if (omega0 < 0.0) {

      double R=length*Fcorot;
      double phi, theta=0.5*M_PI;
      double dens, potl, potr, pott, potp;
      double avg=0.0;
      
      for (int n=0; n<8; n++) {
	phi = 2.0*M_PI/8.0 * n;

	for (cc=comp.components.begin(); cc != comp.components.end(); cc++) {
	  c = *cc;
	
	  if (c->force->geometry == PotAccel::sphere || 
	      c->force->geometry == PotAccel::cylinder) {
	  
	    ((Basis*)c->force)->
	      determine_fields_at_point_sph(R, theta, phi,
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
    for (int i=1; i<=N; i++) {
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

      name = filename;
      name += ".barstat";

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

    Iz = 0.2*mass*(a1*a1 + a2*a2);
    Lz = Iz * omega;

    if (c1)
      Lz0 = c1->angmom[2];
    else
      Lz0 = 0.0;

    posang = 0.0;
    lastomega = omega;
    lasttime = tvel;
    
    if (restart) {

      if (myid == 0) {
	
	// Backup up old file
	string backupfile = name + ".bak";
	string command("cp ");
	command += name + " " + backupfile;
	system(command.c_str());

	// Open new output stream for writing
	ofstream out(name.c_str());
	if (!out) {
	  cout << "UserEBar: error opening new log file <" 
	       << filename << "> for writing\n";
	  MPI_Abort(MPI_COMM_WORLD, 121);
	  exit(0);
	}
	
	// Open old file for reading
	ifstream in(backupfile.c_str());
	if (!in) {
	  cout << "UserEBar: error opening original log file <" 
	       << backupfile << "> for reading\n";
	  MPI_Abort(MPI_COMM_WORLD, 122);
	  exit(0);
	}

	const int linesize = 1024;
	char line[linesize];
	
	in.getline(line, linesize); // Discard header
	in.getline(line, linesize); // Next line

	double Lzp, Lz1, am1;
	bool firstime1 = true;
	while (in) {
	  istringstream ins(line);

	  lastomega = omega;

	  ins >> lasttime;
	  ins >> posang;
	  ins >> omega;
	  ins >> Lz1;
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
	    Lz = Lz1;
	    Lz0 = Lzp;
	    firstime1 = false;
	  }

	  if (lasttime >= tvel) break;

	  out << line << "\n";

	  in.getline(line, linesize); // Next line
	}

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
    }

    firstime = false;
    update = true;

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
    
    if ( fabs(tvel-lasttime) > 2.0*DBL_EPSILON) {
      posang += 0.5*(omega + lastomega)*dtime;
      lastomega = omega;
      lasttime = tvel;
      update = true;
    }
  }

				// Zero thread variables
  for (int n=0; n<nthrds; n++) {
    for (int k=0; k<3; k++) tacc[n][k] = 0.0;
  }

  exp_thread_fork(false);

				// Get full contribution from all threads
  for (int k=0; k<3; k++) acc[k] = acc1[k] = 0.0;
  for (int n=0; n<nthrds; n++) {
    for (int k=0; k<3; k++) acc1[k] += tacc[n][k]/barmass;
  }

				// Get contribution from all processes
  MPI_Allreduce(&acc1[0], &acc[0], 3, MPI_DOUBLE, MPI_SUM, MPI_COMM_WORLD);

				// Backward Euler
  if (monopole) {
    for (int k=0; k<3; k++) {
      bps[k] += vel[k] * (tnow - teval);
      vel[k] += acc[k] * (tnow - teval);
    }
    teval = tnow;
  }

  if (myid==0 && update) 
    {
      ofstream out(name.c_str(), ios::out | ios::app);
      out.setf(ios::scientific);

      out << setw(15) << tvel
	  << setw(15) << posang
	  << setw(15) << omega;

      if (c1)
	out << setw(15) << Lz + Lz0 - c1->angmom[2]
	    << setw(15) << c1->angmom[2];
      else
	out << setw(15) << Lz
	    << setw(15) << 0.0;

      out << setw(15) << amplitude/fabs(amplitude) *  
	0.5*(1.0 + erf( (tvel - Ton )/DeltaT )) *
	0.5*(1.0 - erf( (tvel - Toff)/DeltaT ));

      for (int k=0; k<3; k++) out << setw(15) << bps[k];
      for (int k=0; k<3; k++) out << setw(15) << vel[k];
      for (int k=0; k<3; k++) out << setw(15) << acc[k];
      
      out << endl;
    }

}


void * UserEBar::determine_acceleration_and_potential_thread(void * arg) 
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

  if (table) {
    if (tvel<timeq[0]) {
      afac = ampq[0];
      b5 = b5q[0];
    } else if (tvel>timeq[qlast]) {
      afac = ampq[qlast];
      b5 = b5q[qlast];
    } else {
      afac = odd2(tvel, timeq, ampq, 0);
      b5 = odd2(tvel, timeq, b5q, 0);
    }
  }

  if (fabs(amplitude)>0.0) amp = afac * amplitude/fabs(amplitude) 
			     * 0.5*(1.0 + erf( (tvel - Ton )/DeltaT ))
			     * 0.5*(1.0 - erf( (tvel - Toff)/DeltaT )) ;

  for (int i=nbeg; i<nend; i++) {

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

  }

  return (NULL);
}


extern "C" {
  ExternalForce *makerEBar(string& line)
  {
    return new UserEBar(line);
  }
}

class proxyebar { 
public:
  proxyebar()
  {
    factory["userebar"] = makerEBar;
  }
};

proxyebar p;
