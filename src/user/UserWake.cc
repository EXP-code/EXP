#include <cmath>
#include <sstream>

#include <expand.H>
#include <localmpi.H>

#include <ExternalCollection.H>
#include <Basis.H>
#include <UserWake.H>

#include <pthread.h>  

using namespace std;

UserWake::UserWake(const YAML::Node& conf) : ExternalForce(conf)
{
  first = true;
  filename = "wake";

				// Output surface
  NUMX = 100;
  NUMY = 100;
  XMIN = -1.8;
  XMAX =  1.8;
  YMIN = -1.8;
  YMAX =  1.8;
				// Rotation of surface from equitorial plane
  PHI = 0.0;
  PSI = 0.0;
  THETA = 0.0;
				// Steps between frames
  NSTEP = 10;

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

				// Names of output file
  names.push_back("dens0");
  names.push_back("dens1");
  names.push_back("dens");
  names.push_back("densR");

  names.push_back("potl0");
  names.push_back("potl1");
  names.push_back("potl");
  names.push_back("potlR");

				// Data storage
  npix = NUMX*NUMY;
  for (unsigned i=0; i<names.size(); i++) {
    data0.push_back(vector<float>(npix));
    data1.push_back(vector<float>(npix));
  }

  nbeg = npix*myid/numprocs;
  nend = npix*(myid+1)/numprocs;

  double onedeg = M_PI/180.0;
  Eigen::Matrix3d rotate = return_euler_slater(PHI*onedeg, 
					       THETA*onedeg, 
					       PSI*onedeg, 1);
  Eigen::Vector3d P0, P1;
  P0.setZero();
    

  double dX = (XMAX - XMIN)/(NUMX-1);
  double dY = (YMAX - YMIN)/(NUMY-1);
  double R;

  for (int j=0; j<NUMY; j++) {

    P0[1] = YMIN + dY*j;

    for (int i=0; i<NUMX; i++) {

      P0[0] = XMIN + dX*i;

      P1 = rotate * P0;

      R = sqrt(P0.adjoint()*P0);

				// Location of the rotated plane
      r.push_back(R);
      theta.push_back(acos(P1[2]/R));
      phi.push_back(atan2(P1[1], P1[0]));
    }
  }
    
}

UserWake::~UserWake()
{
  // Nothing so far
}

void UserWake::userinfo()
{
  if (myid) return;		// Return if node master node

  print_divider();

  cout << "** User routine WAKE initialized"
       << " with Components <" << C[0];
  for (int i=1; i<numComp; i++) cout << " " << C[i];
  cout << ">";
  
  cout << ", NUMX=" << NUMX
       << ", NUMY=" << NUMY
       << ", XMIN=" << XMIN
       << ", XMAX=" << XMAX
       << ", YMIN=" << YMIN
       << ", YMAX=" << YMAX
       << ", PHI=" << PHI
       << ", THETA=" << THETA
       << ", PSI=" << PSI
       << ", NSTEP=" << NSTEP
       << ", filename=" << filename
       << endl;
  
  print_divider();
}

void UserWake::initialize()
{
  try {
    for (numComp=0; numComp<1000; numComp++) {
      ostringstream count;
      count << "C(" << numComp+1 << ")";
      if (conf[count.str()])
	C.push_back(conf[count.str()].as<std::string>());
      else break;
    }
    
    if (numComp != (int)C.size()) {
      cerr << "UserWake: error parsing component names, "
	   << "  Size(C)=" << C.size()
	   << "  numRes=" << numComp << endl;
      MPI_Abort(MPI_COMM_WORLD, 122);
    }
    
    if (conf["filename"])       filename           = conf["filename"].as<string>();
    if (conf["NUMX"])           NUMX               = conf["NUMX"].as<int>();
    if (conf["NUMY"])           NUMY               = conf["NUMY"].as<int>();
    if (conf["XMIN"])           XMIN               = conf["XMIN"].as<double>();
    if (conf["XMAX"])           XMAX               = conf["XMAX"].as<double>();
    if (conf["YMIN"])           YMIN               = conf["YMIN"].as<double>();
    if (conf["YMAX"])           YMAX               = conf["YMAX"].as<double>();
    if (conf["PHI"])            PHI                = conf["PHI"].as<double>();
    if (conf["THETA"])          THETA              = conf["THETA"].as<double>();
    if (conf["PSI"])            PSI                = conf["PSI"].as<double>();
    if (conf["NSTEP"])          NSTEP              = conf["NSTEP"].as<int>();
  }
  catch (YAML::Exception & error) {
    if (myid==0) std::cout << "Error parsing parameters in UserWake: "
			   << error.what()         << std::endl
			   << std::string(60, '-') << std::endl
			   << "Config node"        << std::endl
			   << std::string(60, '-') << std::endl
			   << conf                 << std::endl
			   << std::string(60, '-') << std::endl;
    MPI_Finalize();
    exit(-1);
  }
}  

void UserWake::determine_acceleration_and_potential(void)
{
#if HAVE_LIBCUDA==1		// Cuda compatibility
  getParticlesCuda(cC);
#endif
  
  if (first) {
    
    count = 0;
    nlast = this_step;
    nnext = this_step;
    
    if (restart) {
      
      if (myid == 0) {
	
	for (count=0; count<10000; count++) {
	  ostringstream ostr;
	  ostr << outdir << runtag << "." << filename << "." 
	       << names[0] << "." << count;
	  
	  // Try to open stream for writing
	  ifstream in(ostr.str().c_str());
	  if (!in) break;
	}
      }
      
      if (myid==0) 
	cout << "UserWake: beginning at frame=" << count << endl;
      
      MPI_Bcast(&count, 1, MPI_INT, 0, MPI_COMM_WORLD);
	
    }

    first = false;
  }

  if (this_step == nnext) {


    // -----------------------------------------------------------
    // Clean the data store before the first component in the list
    // -----------------------------------------------------------

    if (cC == c0.front()) {
      for (unsigned i=0; i<names.size(); i++) {
	for (int j=0; j<npix; j++) data0[i][j] = data1[i][j] = 0.0;
      }
    }

    // -----------------------------------------------------------
    // Compute the images
    // -----------------------------------------------------------

    // exp_thread_fork(false);

    double dens0, potl0, dens, potl, potr, pott, potp;

    for (int i=nbeg; i<nend; i++) {

      ((Basis *)cC->force)->
	determine_fields_at_point_sph(r[i], theta[i], phi[i], 
				      &dens0, &potl0, 
				      &dens, &potl,
				      &potr, &pott, &potp);
    
      data1[0][i] += dens0;
      data1[1][i] += dens - dens0;
      data1[2][i] += dens;
      
      data1[4][i] += potl0;
      data1[5][i] += potl - potl0;
      data1[6][i] += potl;
    }


    // -----------------------------------------------------------
    // Print the images after the last component in the list
    // -----------------------------------------------------------
    
    if (cC == c0.back()) {

      for (unsigned i=0; i<names.size(); i++) {
      
	MPI_Reduce(&data1[i][0], &data0[i][0], npix, 
		   MPI_FLOAT, MPI_SUM, 0, MPI_COMM_WORLD);
      
	if (myid==0) {		// Print the file
	  
	  float f;

	  for (unsigned i=0; i<names.size(); i++) {

	    ostringstream ostr;
	    ostr << outdir << runtag << "." << filename 
		 << "." << names[i] << "." << count;
	  
	    // Try to open stream for writing
	    ofstream out(ostr.str().c_str());
	    if (out) {
	      out.write((const char *)&NUMX, sizeof(int));
	      out.write((const char *)&NUMY, sizeof(int));
	      out.write((const char *)&(f=XMIN), sizeof(float));
	      out.write((const char *)&(f=XMAX), sizeof(float));
	      out.write((const char *)&(f=YMIN), sizeof(float));
	      out.write((const char *)&(f=YMAX), sizeof(float));

	      for (int j=0; j<npix; j++) {

		if (i==3) {	// Relative density
		  data0[i][j] = data0[1][j];
		  if (data0[0][j]>0.0) data0[i][j] /= fabs(data0[0][j]);
		}
		if (i==7) {	// Relative potential
		  data0[i][j] = data0[5][j];
		  if (data0[4][j]>0.0) data0[i][j] /= fabs(data0[4][j]);
		}
		out.write((const char *)&data0[i][j], sizeof(float));
	      }

	    } else {
	      cerr << "UserWake: error opening <" << ostr.str() << ">" << endl;
	    }
	  }
      
	}

      }

    }

    // -----------------------------------------------------------


    count++;
    nlast = this_step;
    nnext = this_step + NSTEP;
  }

}


void * UserWake::determine_acceleration_and_potential_thread(void * arg) 
{
  int id = *((int*)arg);
  int nbeg1 = nbeg + (nend - nbeg)*(id  )/nthrds;
  int nend1 = nbeg + (nend - nbeg)*(id+1)/nthrds;

  double dens0, potl0, dens, potl, potr, pott, potp;

  for (int i=nbeg1; i<nend1; i++) {

    ((Basis *)cC->force)->
      determine_fields_at_point_sph(r[i], theta[i], phi[i], 
				    &dens0, &potl0, 
				    &dens, &potl,
				    &potr, &pott, &potp);
    
    data1[0][i] += dens0;
    data1[1][i] += dens - dens0;
    data1[2][i] += dens;

    data1[4][i] += potl0;
    data1[5][i] += potl - potl0;
    data1[6][i] += potl;
  }


  return (NULL);
}


extern "C" {
  ExternalForce *makerWake(const YAML::Node& conf)
  {
    return new UserWake(conf);
  }
}

class proxywake { 
public:
  proxywake()
  {
    factory["userwake"] = makerWake;
  }
};

static proxywake p;
