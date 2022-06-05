#include <math.h>
#include <sstream>

#include "expand.H"

#include <Particle.H>
#include <AxisymmetricBasis.H>
#include <ExternalCollection.H>

//
// Rotate the halo about given center
//

class UserRotF : public ExternalForce
{
private:
  
  void determine_acceleration_and_potential(void);
  void * determine_acceleration_and_potential_thread(void * arg);
  void initialize();

  int center, RROT_Flag;
  double RROT, Omega, Rsate;

  bool firstime;
  bool TidalON, CentrifugalON, CoriolisON;

  int Nint, Ngrid;
  double Rmax, Rmin;
  double rtidal;
  double alpha;
  double diverge_rfac;
  int diverge;
  string host_file, filename;
  string name;

  SphericalModelTable *host;

  vector<double> Mtable;
  double dgrid;

  void userinfo();

public:

  //! Constructor
  UserRotF(string &line);

  //! Destructor
  ~UserRotF();

};

UserRotF::UserRotF(string &line) : ExternalForce(line)
{

  //  id = "Halo RotF";	       

  // Set the parameters
  Omega = 0.;			// Rotation frequency
  Rsate = 1.0;                  // Staellite location
  center = 0;                   // 0 is origin, 1 is COM, and 2 is EJ
  RROT_Flag = 1;                // 0 for infinity, 1 for constant as initial RROT, and -1 for variation as rtidal
  RROT = 0.0;			// initial maximum radius in which the particle is affected by rotation 

  TidalON  =  0;                // Turn-on tidal effect by determining the membership of the particle
  Nint   = 1;                   // Timestep interval to recalculate the tidal radius
  Ngrid  = 1000;                // Grid number for satellite potential table               
  Rmax   = 1000.0;              // Maximum radius of the potential table
  Rmin   = 0.0;                 // Minimum radius of the potential table
  host_file = "SLGridSph.model";
  diverge = 0;                  // Use analytic divergence (true/false)
  diverge_rfac = 1.0;           // Exponent for profile divergence

  firstime = true;

  CentrifugalON = 0;            // Turn-on the centrifugal force
  CoriolisON = 0;               // Turn-on the coriolis force

  alpha = 0.0;                  // Come up in derivation of the tidal radius, Centrifugal OFF  alpha=0
                                // Accoding to derivation of the tidal radius


  filename = outdir + "Tidal." + runtag; // Output file name


  initialize();

  host = new SphericalModelTable(host_file, diverge, diverge_rfac);

  dgrid = (Rmax-Rmin)/double(Ngrid);

  if(!CentrifugalON) alpha = 0.0;

  userinfo();
}

UserRotF::~UserRotF()
{
  delete host;
}

void UserRotF::userinfo()
{
  if (myid) return;		// Return if node master node

  print_divider();

  cout << "** User routine RotF initialized, " ;
  cout << " Omega= " << Omega;
  cout << " , the rotation center is ";
  if (center == 0) 
    cout << " the Origin";
  if (center == 1) 
    cout << " the COM center";
  if (center == 2) 
    cout << " the EJ center";
  if ((center < 0) || (center > 2)) {
    throw GenericError("You need to choose a proper center (Origin, COM, EJ)",
		       __FILE__, __LINE__);
  } 
  cout << " and R_rot is ";
  if (RROT_Flag > 0.0)  
    cout << RROT;
  if (RROT_Flag == 0.0)
    cout << " Infinity";
  if (RROT_Flag < 0.0)
    cout << " Rtidal";
  if(TidalON){
    cout << " Host file= " << host_file << " diverge="  << diverge;
    cout << " diverge_rfac= " << diverge_rfac;
    cout << endl;
  }
  if(CentrifugalON) 
    cout << " Centrifugal force is on: alpha = " << alpha << endl;

  print_divider();
}

void UserRotF::initialize()
{
  string val;

  if (get_value("center",val))      center = atoi(val.c_str());
  if (get_value("Omega", val))      Omega = atof(val.c_str());
  if (get_value("Rsate", val))      Rsate = atof(val.c_str());
  if (get_value("RROT_Flag", val))  RROT_Flag = atoi(val.c_str());
  if (get_value("RROT", val))       RROT = atof(val.c_str());
  if (RROT <= 0.0) RROT_Flag = 0;
  
  if (get_value("TidalON",val))         TidalON = atol(val);
  if (get_value("Nint",val))            Nint = atoi(val.c_str());
  if (get_value("Ngrid",val))           Ngrid = atoi(val.c_str());
  if (get_value("alpha", val))          alpha = atof(val.c_str());
  if (get_value("Rmin", val))           Rmin = atof(val.c_str());
  if (get_value("Rmax", val))           Rmax = atof(val.c_str());

  if (get_value("host_file", val))      host_file = val;
  if (get_value("diverge", val))        diverge = atoi(val.c_str());
  if (get_value("diverge_rfac", val))   diverge_rfac = atof(val.c_str());
  if (get_value("filename", val))       filename = val;

  if (get_value("CentrifugalON",val))   CentrifugalON = atol(val);
  if (get_value("CoriolisON",val))      CoriolisON = atol(val);


}


void UserRotF::determine_acceleration_and_potential(void)
{

  MPI_Status status;

  // output file for tidal radius
  if (firstime && TidalON) {

    name = filename;
    name += ".dat";

    if (restart) {
      double lasttime, lastxe;
      if (myid==0) {

	// Backup up old file
	string curfile = outdir + name;
	string backupfile = curfile + ".bak";
	string command("cp ");
	command += curfile + " " + backupfile;
	system(command.c_str());

	// Open new output stream for writing
	ofstream out(curfile.c_str());
	if (!out) {
	  throw FileCreateError(curfile, "UserRotF: error opening new log file",
			      __FILE__, __LINE__);
	}

        // Open old file for reading
        ifstream in(backupfile.c_str());
        if (!in) {
	  throw FileCreateError(backupfile,
				"UserRotF: error opening original log file",
				__FILE__, __LINE__);
        }

        const int linesize = 1024;
        char line[linesize];

        in.getline(line, linesize); // Discard header
        in.getline(line, linesize); // Next line

        while (in) {
          istringstream ins(line);

          ins >> lasttime;
          ins >> lastxe;

          if (lasttime >= tnow) break;

          out << line << "\n";
          in.getline(line, linesize); // Next line
	}
      }

      MPI_Bcast(&lasttime, 1, MPI_DOUBLE, 0, MPI_COMM_WORLD);
      MPI_Bcast(&lastxe, 1, MPI_DOUBLE, 0, MPI_COMM_WORLD);

      cC->rtrunc = lastxe;
      firstime=false;
    } else {
      if(myid == 0) {
        ofstream out(string(outdir+name).c_str(), ios::out | ios::app);

        out << setw(15) << "# Time"
            << setw(15) << "Tidal radius"
            << endl;
      }
      firstime=false;
    }
  }

  Mtable = vector<double>(Ngrid);  
  if(this_step % Nint == 0  && TidalON) 
    for(int i=0; i<Ngrid; i++) Mtable[i] = 0.0;

#if HAVE_LIBCUDA==1		// Cuda compatibility
  getParticlesCuda(cC);
#endif

  exp_thread_fork(false);


  // Computing tidal radius
  if(this_step % Nint == 0 && mlevel==0 && TidalON) {
    vector<double> Mtotable;

    // Gethering the mass profile
    if (myid==0)  {
      Mtotable = Mtable;
      Mtable =vector<double>(Ngrid);
    }

    for (int n=1; n<numprocs; n++) {
      if (myid==n) 
	MPI_Send(&Mtable[0], Ngrid, MPI_DOUBLE, 0, 1234, MPI_COMM_WORLD);

      if (myid==0) {
	MPI_Recv(&Mtable[0], Ngrid, MPI_DOUBLE, n, 1234, MPI_COMM_WORLD, &status);
	for(int i=0; i<Ngrid; i++) Mtotable[i] += Mtable[i];
      }
    }
    MPI_Barrier(MPI_COMM_WORLD);


    // Computing the potential profile of the satellite
    if (myid==0) {
      vector<double> R(Ngrid), D(Ngrid), M(Ngrid), P(Ngrid), P1(Ngrid);
      vector<double> r1(Ngrid), r2(Ngrid);

      for (int i=0; i<Ngrid; i++) {
	r1[i] = Rmin + dgrid*i;
	r2[i] = Rmin + dgrid*(i+1);
	R[i] = Rmin + dgrid*(0.5+i);
	
	D[i] = Mtotable[i] / (4.0*M_PI/3.0 * (pow(r2[i], 3.0) - pow(r1[i], 3.0)));
      }

      for (int i=0; i<Ngrid; i++) {

	if (i==0) M[i] = 2.0*M_PI*dgrid*R[i]*R[i]*D[i];
	else M[i] = M[i-1] + 2.0*M_PI*(dgrid*R[i-1]*R[i-1]*D[i-1] +
				       dgrid*R[i  ]*R[i  ]*D[i  ] );

	if (i==0) P1[i] = 2.0*M_PI*dgrid*D[i]*R[i];
	else P1[i] = P1[i-1] + 2.0*M_PI*(dgrid*R[i-1]*D[i-1] +
					 dgrid*R[i  ]*D[i  ] );
      }

      for (int i=Ngrid-1; i>0; i--)
	P[i] = -M[i]/R[i] - (P1[Ngrid-1] - 0.5*(P1[i] + P1[i-1]));
      P[0] = -M[0]/R[0]  - (P1[Ngrid-1] - 0.5*P1[0]);

      // Computing Effective potential
      double Potmax, Pot;
      double rho_host, m_host, Pot_extern;
      
      rho_host = host->get_density(Rsate);
      m_host  = host->get_mass(Rsate);
      Pot_extern = 4.0*M_PI*rho_host - (2.0+alpha)*m_host/pow(Rsate,3);

      rtidal = 0.0;
      Potmax = -1.0e20;
      for(int i=0; i<Ngrid; i++){
	Pot=P[i] + 0.5*R[i]*R[i]*Pot_extern;
	if(Pot >= Potmax) {
	  rtidal = R[i];
	  Potmax = Pot;
	}
      }

      // Write onto Tidal output
      ofstream out(string(outdir+name).c_str(), ios::out | ios::app);
      out.setf(ios::scientific);

      out << setw(15) << tnow
	  << setw(15) << rtidal
	  << endl;
    }

    MPI_Bcast(&rtidal, 1, MPI_DOUBLE, 0, MPI_COMM_WORLD);
    cC->rtrunc=rtidal;
    if(RROT_Flag < 0) RROT = rtidal;
  }

}


void * UserRotF::determine_acceleration_and_potential_thread(void * arg) 
{

  //if (!this_step) return(NULL);

  unsigned nbodies = cC->Number();
  int id = *((int*)arg);
  int nbeg = nbodies*id/nthrds;
  int nend = nbodies*(id+1)/nthrds;

  double rr, angle1, angle2;
  double pos[3], vel[3], dvel[3];


  for (int i=nbeg; i<nend; i++) {
    for(int k=0; k<3; k++){
      if(center == 0) {
	pos[k] = cC->Pos(i, k);
	vel[k] = cC->Vel(i, k);
      }
      if(center == 1) {
	pos[k] = cC->Pos(i, k) - cC->com0[k];
	vel[k] = cC->Vel(i, k) - cC->cov0[k];
      }
      if(center ==2) {
	pos[k] = cC->Pos(i, k) - (cC->com0[k] + cC->center[k]);
	vel[k] = cC->Vel(i, k) - (cC->cov0[k] + cC->center[k]);
      }
    }

    rr=0.0;
    for(int k=0; k<3; k++) rr += pos[k]*pos[k];
    rr = sqrt(rr);

    // Tidal radius
    if(this_step % Nint == 0  && TidalON) {
      int nn;
      if(rr <= Rmin) nn = 0;
      else nn = int((rr-Rmin)/dgrid);
      if(nn < Ngrid) Mtable[nn] += cC->Mass(i);
    }

    
    // Centrifugal and Coriolis forces
    if((rr < RROT) || !RROT_Flag){
      angle1 = Omega*dtime;
      angle2 = Omega*Omega*dtime;
      for(int k=0; k<3; k++) dvel[k]=0;

      // Centrifugal force is taken account in tidal radius calculation.
      //      if(CentrifugalON) {
      //	dvel[0] += (pos[0]*angle2);
      //	dvel[1] += (pos[1]*angle2);
      //	dvel[2] += 0.0;
      //      }

      if(CoriolisON) {
	dvel[0] += (2.0*vel[1]*angle1);
	dvel[1] += (-2.0*vel[0]*angle1);
	dvel[2] += 0.0;
      }
	
      // Add additional force
      cC->AddVel(i, dvel);
    }
  }

  return (NULL);
}


extern "C" {
  ExternalForce *makerRotF(string& line)
  {
    return new UserRotF(line);
  }
}

class proxyhalo { 
public:
  proxyhalo()
  {
    factory["userrotf"] = makerRotF;
  }
};

static proxyhalo p;
