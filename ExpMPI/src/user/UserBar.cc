#include <math.h>
#include <strstream>

#include "expand.h"

#include <UserBar.H>

UserBar::UserBar(string &line) : ExternalForce(line)
{
  id = "RotatingBar";

  length = 0.5;			// Bar length
  amplitude = 0.3;		// Bar amplitude
  Ton = -20.0;			// Turn on start time
  Toff = 200.0;			// Turn off start time
  DeltaT = 1.0;			// Turn on duration
  Fcorot  = 1.0;		// Corotation factor
  fixed = false;		// Constant pattern speed
  shallow = false;		// Use shallow form of the bar potential
  filename = "BarRot";		// Output file name

  firstime = true;

  com_name = "";		// Default component for com

  initialize();

  if (com_name.size()>0) {
				// Look for the fiducial component
    bool found = false;
    list<Component*>::iterator cc;
    Component *c;
    for (cc=comp.components.begin(); cc != comp.components.end(); cc++) {
      c = *cc;
      if ( !com_name.compare(c->name) ) {
	c0 = c;
	found = true;
      break;
      }
    }

    if (!found) {
      cerr << "Process " << myid << ": can't find desired component <"
	   << com_name << ">" << endl;
      MPI_Abort(MPI_COMM_WORLD, 35);
    }

  }
  else
    c0 = NULL;

  userinfo();
}

UserBar::~UserBar()
{
}

void UserBar::userinfo()
{
  if (myid) return;		// Return if node master node
  cout << "User routine initialized: rotating bar, " ;
  if (fixed)
    cout << "fixed pattern speed, ";
  else
    cout << "fixed corotation fraction, ";
  if (shallow)
    cout << "shallow potential, ";
  else
    cout << "standard potential, ";
  if (c0) 
    cout << "center on component <" << com_name << ">" << endl;
  else
    cout << "center on origin" << endl;
}

void UserBar::initialize()
{
  string val;

  if (get_value("comname", val))	com_name = val;
  if (get_value("length", val))		length = atof(val.c_str());
  if (get_value("amp", val))		amplitude = atof(val.c_str());
  if (get_value("Ton", val))		Ton = atof(val.c_str());
  if (get_value("Toff", val))		Toff = atof(val.c_str());
  if (get_value("DeltaT", val))		DeltaT = atof(val.c_str());
  if (get_value("Fcorot", val))		Fcorot = atof(val.c_str());
  if (get_value("fixed", val)) {
    if (atoi(val.c_str())) fixed = true;
    else fixed = false;
  }
  if (get_value("shallow", val)) {
    if (atoi(val.c_str())) shallow = true;
    else shallow = false;
  }
  if (get_value("filename", val))	filename = val;

}


void UserBar::determine_acceleration_and_potential(void)
{
				// Write to bar state file, if true
  bool update = false;

  c0->get_angmom();		// Tell component to compute angular momentum
  // cout << "Lz=" << c0->angmom[2] << endl; // debug

  if (firstime) {
    
    list<Component*>::iterator cc;
    Component *c;

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

    omega = sqrt(avg/R);

    const int N = 100;
    LegeQuad gq(N);

    double a1 = length;
    double a2 = 0.5*a1;
    double a3 = 0.1*a2;

    double geom = pow(a1*a2*a3, 1.0/3.0);

    double A12 = a1*a1/geom/geom;
    double A22 = a2*a2/geom/geom;
    double A32 = a3*a3/geom/geom;

    double u, d, t, denom, ans1=0.0, ans2=0.0;
    double mass = fabs(amplitude);
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

    }

    double rho = mass/(4.0*M_PI/3.0*a1*a2*a3);
    double b1 = M_PI*rho*sqrt(2.0*M_PI/15.0)*(ans1 - ans2);
    double b25 = 0.4*a1*a2*a3*(a2*a2 - a1*a1)/(ans1 - ans2);

    b5 = pow(b25, 0.2);
    afac = 2.0 * b1;

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
	    << endl;
      }
    }

    Iz = 0.2*mass*(a1*a1 + a2*a2);
    Lz = Iz * omega;

    Lz0 = c0->angmom[2];

    posang = 0.0;
    lastomega = omega;
    lasttime = tvel;
    
    if (restart) {

      if (myid == 0) {
	
	const int linesize = 1024;
	char line[linesize];
	
	ifstream in(name.c_str());
	if (!in) {
	  cerr << "User_bar_slow: can't open barstat file<" << name << ">\n";
	  MPI_Abort(MPI_COMM_WORLD, 101);
	  exit(0);
	}

	in.getline(line, linesize); // Discard header
	in.getline(line, linesize); // Next line

	double Lz1, am1;
	bool firstime1 = true;
	while (in) {
	  istrstream ins(line);

	  lastomega = omega;

	  ins >> lasttime;
	  ins >> posang;
	  ins >> omega;
	  ins >> Lz1;
	  ins >> am1;

	  if (firstime1) {
	    Lz = Lz1;
	    Lz0 = am1;
	    firstime1 = false;
	  }

	  in.getline(line, linesize); // Next line
	}
      }

      MPI_Bcast(&lasttime, 1, MPI_DOUBLE, 0, MPI_COMM_WORLD);
      MPI_Bcast(&posang, 1, MPI_DOUBLE, 0, MPI_COMM_WORLD);
      MPI_Bcast(&omega, 1, MPI_DOUBLE, 0, MPI_COMM_WORLD);
      MPI_Bcast(&lastomega, 1, MPI_DOUBLE, 0, MPI_COMM_WORLD);
      MPI_Bcast(&Lz, 1, MPI_DOUBLE, 0, MPI_COMM_WORLD);
      MPI_Bcast(&Lz0, 1, MPI_DOUBLE, 0, MPI_COMM_WORLD);
    }

    firstime = false;
    update = true;

  } else {

    if (!fixed)
      omega = (Lz + Lz0 - c0->angmom[2])/Iz;
    else
      omega = lastomega;
    
    if ( fabs(tvel-lasttime) > 2.0*DBL_EPSILON) {
      posang += 0.5*(omega + lastomega)*dtime;
      lastomega = omega;
      lasttime = tvel;
      update = true;
    }
  }

  exp_thread_fork(false);

  if (myid==0 && update) 
    {
      ofstream out(name.c_str(), ios::out | ios::app);
      out.setf(ios::scientific);

      out << setw(15) << tvel
	  << setw(15) << posang
	  << setw(15) << omega
	  << setw(15) << Lz + Lz0 - c0->angmom[2]
	  << setw(15) << c0->angmom[2]
	  << setw(15) << amplitude *  
	0.5*(1.0 + erf( (tvel - Ton )/DeltaT )) *
	0.5*(1.0 - erf( (tvel - Toff)/DeltaT ))
	  << endl;
    }

}


void * UserBar::determine_acceleration_and_potential_thread(void * arg) 
{
  int nbodies = particles->size();
  int id = *((int*)arg);
  int nbeg = 1+nbodies*id/nthrds;
  int nend = nbodies*(id+1)/nthrds;

  double fac, ffac, amp = afac * amplitude/fabs(amplitude) 
    * 0.5*(1.0 + erf( (tvel - Ton )/DeltaT ))
    * 0.5*(1.0 - erf( (tvel - Toff)/DeltaT )) ;
  double xx, yy, zz, rr, nn,pp;
  vector<double> pos(3); 
  double cos2p = cos(2.0*posang);
  double sin2p = sin(2.0*posang);

  for (int i=nbeg; i<nend; i++) {

    if ((*particles)[i].freeze()) continue;

    if (c0)
      for (int k=0; k<3; k++) pos[k] = (*particles)[i].pos[k] - c0->com[k];
    else
      for (int k=0; k<3; k++) pos[k] = (*particles)[i].pos[k];
    
    xx = pos[0];
    yy = pos[1];
    zz = pos[2];
    rr = sqrt( xx*xx + yy*yy + zz*zz );

    if (shallow) {
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

    (*particles)[i].acc[0] += ffac*
      ( 2.0*( xx*cos2p + yy*sin2p)*fac - 5.0*nn*xx );
    
    (*particles)[i].acc[1] += ffac*
      ( 2.0*(-yy*cos2p + xx*sin2p)*fac - 5.0*nn*yy );

    (*particles)[i].acc[2] += ffac*
      ( -5.0*nn*zz );
    
    (*particles)[i].potext += -ffac*pp*fac;
    
  }

  return (NULL);
}


extern "C" {
  ExternalForce *makerBar(string& line)
  {
    return new UserBar(line);
  }
}

class proxy { 
public:
  proxy()
  {
    factory["userbar"] = makerBar;
  }
};

proxy p;
