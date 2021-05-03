/*****************************************************************************
 *  Description:
 *  -----------
 *
 *  Read in coefficients and compute gnuplot slices, 
 *  and compute volume for VTK rendering
 *
 *
 *  Call sequence:
 *  -------------
 *
 *  Parameters:
 *  ----------
 *
 *
 *  Returns:
 *  -------
 *
 *
 *  Notes:
 *  -----
 *
 *
 *  By:
 *  --
 *
 *  MDW 11/28/08
 *
 ***************************************************************************/

				// C++/STL headers
#include <cstdlib>
#include <iostream>
#include <iomanip>
#include <fstream>
#include <sstream>
#include <cmath>
#include <string>

using namespace std;

				// Boost stuff

#include <boost/program_options.hpp>
#include <boost/filesystem.hpp>

namespace po = boost::program_options;

                                // System libs
#include <sys/time.h>
#include <sys/resource.h>

				// MDW classes
#include <Vector.h>
#include <numerical.h>
#include "Particle.h"
#include <PSP2.H>
#include <interp.h>
#include <massmodel.h>
#include <SphereSL.H>
#include <VtkGrid.H>
#include <localmpi.h>
#include <foarray.H>

// Variables not used but needed for linking
//
int VERBOSE = 0;
int nthrds = 1;
int this_step = 0;
unsigned multistep = 0;
unsigned maxlev = 100;
int mstep = 1;
int Mstep = 1;
vector<int> stepL(1, 0), stepN(1, 1);
char threading_on = 0;
pthread_mutex_t mem_lock;
pthread_mutex_t coef_lock;
std::string outdir, runtag;
double tpos = 0.0;
double tnow = 0.0;
  
// Globals
//
std::string OUTFILE;
double RMIN, RMAX;
int OUTR, LMAX, NMAX, MMAX, L1, L2;
bool VOLUME, SURFACE, PROBE;

// Temporary center offset
// std::vector<double> c0 = {0.0, 0.0, -0.00217};

std::vector<double> c0 = {0.0, 0.0, 0.0};


class Histogram
{
public:

  std::vector<double> dataXY, dataXZ, dataYZ;
  double R, dR;
  int N;
  
  Histogram(int N, double R) : N(N), R(R)
  {
    N = std::max<int>(N, 2);
    dR = 2.0*R/(N-1);		// Want grid points to be on bin centers

    dataXY.resize(N*N);
    dataXZ.resize(N*N);
    dataYZ.resize(N*N);

    Reset();
  }

  void Reset() {
    std::fill(dataXY.begin(), dataXY.end(), 0.0);
    std::fill(dataXZ.begin(), dataXZ.end(), 0.0);
    std::fill(dataYZ.begin(), dataYZ.end(), 0.0);
  }

  void Syncr() { 
    if (myid==0) {
      MPI_Reduce(MPI_IN_PLACE, &dataXY[0], dataXY.size(), MPI_DOUBLE, MPI_SUM, 0, MPI_COMM_WORLD);
      MPI_Reduce(MPI_IN_PLACE, &dataXZ[0], dataXZ.size(), MPI_DOUBLE, MPI_SUM, 0, MPI_COMM_WORLD);
      MPI_Reduce(MPI_IN_PLACE, &dataYZ[0], dataYZ.size(), MPI_DOUBLE, MPI_SUM, 0, MPI_COMM_WORLD);
    }
    else {
      MPI_Reduce(&dataXY[0],         NULL, dataXY.size(), MPI_DOUBLE, MPI_SUM, 0, MPI_COMM_WORLD);
      MPI_Reduce(&dataXZ[0],         NULL, dataXZ.size(), MPI_DOUBLE, MPI_SUM, 0, MPI_COMM_WORLD);
      MPI_Reduce(&dataYZ[0],         NULL, dataYZ.size(), MPI_DOUBLE, MPI_SUM, 0, MPI_COMM_WORLD);
    }

  }

  void Add(double x, double y, double z, double m)
  {
    if (x < -R-0.5*dR or x >= R+0.5*dR or
	y < -R-0.5*dR or y >= R+0.5*dR or
	z < -R-0.5*dR or z >= R+0.5*dR) return;

    int indX = static_cast<int>(floor((x + R + 0.5*dR)/dR));
    int indY = static_cast<int>(floor((y + R + 0.5*dR)/dR));
    int indZ = static_cast<int>(floor((z + R + 0.5*dR)/dR));

    indX = std::max<int>(indX, 0);
    indY = std::max<int>(indY, 0);
    indZ = std::max<int>(indZ, 0);

    indX = std::min<int>(indX, N-1);
    indY = std::min<int>(indY, N-1);
    indZ = std::min<int>(indZ, N-1);

    dataXY[indX*N + indY] += m;
    dataXZ[indX*N + indZ] += m;
    dataYZ[indY*N + indZ] += m;
  }

};


void add_particles(PSPptr psp, int& nbods, vector<Particle>& p, Histogram& h)
{
  if (myid==0) {

    int nbody = nbods/numprocs;
    int nbody0 = nbods - nbody*(numprocs-1);

				// Send number of bodies to be received
				// by eacn non-root node
    MPI_Bcast(&nbody, 1, MPI_INT, 0, MPI_COMM_WORLD);

    vector<Particle>        t(nbody);
    vector<double>        val(nbody);
    vector<unsigned long> seq(nbody);

    SParticle *part = psp->GetParticle();
    Particle bod;
    
    //
    // Root's particles
    //
    for (int i=0; i<nbody0; i++) {
      if (part==0) {
	cerr << "Error reading particle [n=" << 0 << ", i=" << i << "]" << endl;
	exit(-1);
      }

      // Make a Particle
      //
      bod.mass = part->mass();
      for (int k=0; k<3; k++) bod.pos[k] = part->pos(k) - c0[k];
      for (int k=0; k<3; k++) bod.vel[k] = part->vel(k);
      bod.indx = part->indx();
      p.push_back(bod);

      part = psp->NextParticle();

      // Add to histogram
      //
      if (part) h.Add(part->pos(0) - c0[0],
		      part->pos(1) - c0[1],
		      part->pos(2) - c0[2],
		      part->mass());
    }

    //
    // Send the rest of the particles to the other nodes
    //
    for (int n=1; n<numprocs; n++) {
      
      for (int i=0; i<nbody; i++) {
	if (part==0) {
	  cerr << "Error reading particle [n=" 
	       << n << ", i=" << i << "]" << endl;
	  exit(-1);
	}
	t[i].mass = part->mass();
	for (int k=0; k<3; k++) t[i].pos[k] = part->pos(k) - c0[k];
	for (int k=0; k<3; k++) t[i].vel[k] = part->vel(k);
	part = psp->NextParticle();
      }
  
      for (int i=0; i<nbody; i++) val[i] = t[i].mass;
      MPI_Send(&val[0], nbody, MPI_DOUBLE, n, 11, MPI_COMM_WORLD);

      for (int i=0; i<nbody; i++) val[i] = t[i].pos[0];
      MPI_Send(&val[0], nbody, MPI_DOUBLE, n, 12, MPI_COMM_WORLD);

      for (int i=0; i<nbody; i++) val[i] = t[i].pos[1];
      MPI_Send(&val[0], nbody, MPI_DOUBLE, n, 13, MPI_COMM_WORLD);

      for (int i=0; i<nbody; i++) val[i] = t[i].pos[2];
      MPI_Send(&val[0], nbody, MPI_DOUBLE, n, 14, MPI_COMM_WORLD);

      for (int i=0; i<nbody; i++) val[i] = t[i].vel[0];
      MPI_Send(&val[0], nbody, MPI_DOUBLE, n, 15, MPI_COMM_WORLD);

      for (int i=0; i<nbody; i++) val[i] = t[i].vel[1];
      MPI_Send(&val[0], nbody, MPI_DOUBLE, n, 16, MPI_COMM_WORLD);

      for (int i=0; i<nbody; i++) val[i] = t[i].vel[2];
      MPI_Send(&val[0], nbody, MPI_DOUBLE, n, 17, MPI_COMM_WORLD);

      for (int i=0; i<nbody; i++) seq[i] = t[i].indx;
      MPI_Send(&seq[0], nbody, MPI_UNSIGNED_LONG, n, 18, MPI_COMM_WORLD);
    }

  } else {

    int nbody;
    MPI_Bcast(&nbody, 1, MPI_INT, 0, MPI_COMM_WORLD);
    vector<Particle>        t(nbody);
    vector<double>        val(nbody);
    vector<unsigned long> seq(nbody);
				// Get and pack

    MPI_Recv(&val[0], nbody, MPI_DOUBLE, 0, 11, 
	     MPI_COMM_WORLD, MPI_STATUS_IGNORE);
    for (int i=0; i<nbody; i++) t[i].mass = val[i];

    MPI_Recv(&val[0], nbody, MPI_DOUBLE, 0, 12, 
	     MPI_COMM_WORLD, MPI_STATUS_IGNORE);
    for (int i=0; i<nbody; i++) t[i].pos[0] = val[i];

    MPI_Recv(&val[0], nbody, MPI_DOUBLE, 0, 13, 
	     MPI_COMM_WORLD, MPI_STATUS_IGNORE);
    for (int i=0; i<nbody; i++) t[i].pos[1] = val[i];

    MPI_Recv(&val[0], nbody, MPI_DOUBLE, 0, 14, 
	     MPI_COMM_WORLD, MPI_STATUS_IGNORE);
    for (int i=0; i<nbody; i++) t[i].pos[2] = val[i];

    MPI_Recv(&val[0], nbody, MPI_DOUBLE, 0, 15, 
	     MPI_COMM_WORLD, MPI_STATUS_IGNORE);
    for (int i=0; i<nbody; i++) t[i].vel[0] = val[i];

    MPI_Recv(&val[0], nbody, MPI_DOUBLE, 0, 16, 
	     MPI_COMM_WORLD, MPI_STATUS_IGNORE);
    for (int i=0; i<nbody; i++) t[i].vel[1] = val[i];

    MPI_Recv(&val[0], nbody, MPI_DOUBLE, 0, 17, 
	     MPI_COMM_WORLD, MPI_STATUS_IGNORE);
    for (int i=0; i<nbody; i++) t[i].vel[2] = val[i];

    MPI_Recv(&seq[0], nbody, MPI_UNSIGNED_LONG, 0, 18, 
	     MPI_COMM_WORLD, MPI_STATUS_IGNORE);
    for (int i=0; i<nbody; i++) t[i].indx = seq[i];

    p.insert(p.end(), t.begin(), t.end());

    // Add to histogram
    //
    for (auto part : t)
      h.Add(part.pos[0], part.pos[1], part.pos[2], part.mass);
  }

  // Synchronize histogram
  //
  h.Syncr();
}

void partition(PSPptr psp, std::string& name, vector<Particle>& p, Histogram& h)
{
  p.erase(p.begin(), p.end());

  PSPstanza* stanza;

  int iok = 1, nbods = 0;

  if (myid==0) {
    stanza = psp->GetNamed(name);
    if (stanza==0)
      iok = 0;
    else
      nbods = stanza->comp.nbod;
  }

  MPI_Bcast(&iok, 1, MPI_INT, 0, MPI_COMM_WORLD);

  if (iok==0) {
    if (myid==0) std::cerr << "Could not find component named <" << name << ">" << std::endl;
    MPI_Finalize();
    exit(-1);
  }

  add_particles(psp, nbods, p, h);
}

void write_output(SphereSL& ortho, int icnt, double time, Histogram& histo)
{
  unsigned ncnt = 0;

  // ==================================================
  // Setup for output files
  // ==================================================
  
  ostringstream sstr;
  sstr << "." << std::setw(4) << std::setfill('0') << std::right << icnt;

  const int nout1 = 10;
  const int nout2 = 13;
  string suffix[13] = {"p0", "p1", "p", "fr", "ft", "fp", "d0", "d1", "d", "dd",
		       "histoXY", "histoXZ", "histoYZ"};

  if (VOLUME) {
      
    // ==================================================
    // Write volume density
    // ==================================================
    
    double v;
    int valid = 1;
      
    double dR = 2.0*RMAX/(OUTR-1);
    double x, y, z, r, phi, costh;
    double p0, p1, d0, d1, pl, fr, ft, fp;
    
    std::vector<double> data(nout1*OUTR*OUTR*OUTR, 0.0);
    
    for (int k=0; k<OUTR; k++) {
      
      z = -RMAX + dR*k;
      
      for (int l=0; l<OUTR; l++) {
	
	y = -RMAX + dR*l;
	
	for (int j=0; j<OUTR; j++) {
	  
	  x = -RMAX + dR*j;
	  
	  r = sqrt(x*x + y*y + z*z) + 1.0e-18;
	  costh = z/r;
	  phi = atan2(y, x);
	  
	  ortho.all_eval(r, costh, phi, d0, d1, p0, p1, fr, ft, fp, L1, L2);
	  
	  data[((0*OUTR + k)*OUTR + l)*OUTR + j] = p0;
	  data[((1*OUTR + k)*OUTR + l)*OUTR + j] = p1;
	  data[((2*OUTR + k)*OUTR + l)*OUTR + j] = fr;
	  data[((3*OUTR + k)*OUTR + l)*OUTR + j] = ft;
	  data[((4*OUTR + k)*OUTR + l)*OUTR + j] = fp;
	  data[((5*OUTR + k)*OUTR + l)*OUTR + j] = d0;
	  data[((6*OUTR + k)*OUTR + l)*OUTR + j] = d1;
	  if (d0>0.0)
	    data[((7*OUTR + k)*OUTR + l)*OUTR + j] = d1/d0;
	  else
	    data[((7*OUTR + k)*OUTR + l)*OUTR + j] = 0.0;
	}
      }
    }
    
    
    if (myid==0)
      MPI_Reduce(MPI_IN_PLACE, &data[0], nout1*OUTR*OUTR*OUTR, MPI_DOUBLE, MPI_SUM, 
		 0, MPI_COMM_WORLD);
    
    MPI_Reduce(&data[0], NULL, nout1*OUTR*OUTR*OUTR, MPI_DOUBLE, MPI_SUM, 
	       0, MPI_COMM_WORLD);
    
    if (myid==0) {

      VtkGrid vtk(OUTR, OUTR, OUTR, -RMAX, RMAX, -RMAX, RMAX, -RMAX, RMAX);

      std::vector<double> tmp(OUTR*OUTR*OUTR);

      for (int n=0; n<nout1; n++) {
	for (int k=0; k<OUTR; k++) {
	  for (int l=0; l<OUTR; l++) {
	    for (int j=0; j<OUTR; j++) {
	      tmp[(j*OUTR + l)*OUTR + k] = data[((n*OUTR + k)*OUTR + l)*OUTR + j];
	    }
	  }
	}
	vtk.Add(tmp, suffix[n]);
      }

      std::ostringstream sout;
      sout << outdir + "/" + OUTFILE + "_volume";
      vtk.Write(sout.str());
    }

  }
  
  if (SURFACE) {
    
    // ==================================================
    // Write surface profile
    //   --- in plane ---
    // ==================================================
    
    double v;
    float f;
    
    double dR = 2.0*RMAX/OUTR;
    double x, y, z=0.0, r, phi, costh;
    double p0, p1, d0, d1, fr, ft, fp;
    
    vector<double> data(nout1*OUTR*OUTR, 0.0);
    
    for (int l=0; l<OUTR; l++) {
      
      y = -RMAX + dR*(0.5+l);
      
      for (int j=0; j<OUTR; j++) {
	
	if ((ncnt++)%numprocs == myid) {
	  
	  x = -RMAX + dR*(0.5+j);
	  
	  r = sqrt(x*x + y*y + z*z) + 1.0e-18;
	  costh = z/r;
	  phi = atan2(y, x);

	  ortho.all_eval(r, costh, phi, d0, d1, p0, p1, fr, ft, fp, L1, L2);
	  
	  data[(0*OUTR+l)*OUTR+j] = p0;
	  data[(1*OUTR+l)*OUTR+j] = p1;
	  data[(2*OUTR+l)*OUTR+j] = p0 + p1;
	  data[(3*OUTR+l)*OUTR+j] = fr;
	  data[(4*OUTR+l)*OUTR+j] = ft;
	  data[(5*OUTR+l)*OUTR+j] = fp;
	  data[(6*OUTR+l)*OUTR+j] = d0;
	  data[(7*OUTR+l)*OUTR+j] = d1;
	  data[(8*OUTR+l)*OUTR+j] = d0 + d1;

	  if (d0>0.0)
	    data[(9*OUTR+l)*OUTR+j] = d1/d0;
	  else
	    data[(9*OUTR+l)*OUTR+j] = 0.0;
	}
      }
    }
    
    if (myid==0) 
      MPI_Reduce(MPI_IN_PLACE, &data[0], nout1*OUTR*OUTR,
		 MPI_DOUBLE, MPI_SUM, 0, MPI_COMM_WORLD);
    else
      MPI_Reduce(&data[0], NULL, nout1*OUTR*OUTR,
		 MPI_DOUBLE, MPI_SUM, 0, MPI_COMM_WORLD);
    
    
    if (myid==0) {
      
      std::vector<double> dataXY(OUTR*OUTR);

      VtkGrid vtkXY(OUTR, OUTR, 1, -RMAX, RMAX, -RMAX, RMAX, 0, 0);

      for (int n=0; n<nout1; n++) {
	for (int j=0; j<OUTR; j++) {
	  for (int l=0; l<OUTR; l++) {
	    dataXY[j*OUTR + l] = data[(n*OUTR+j)*OUTR+l];
	  }
	}
	vtkXY.Add(dataXY, suffix[n]);
      }

      vtkXY.Add(histo.dataXY, suffix[10]);
      vtkXY.Add(histo.dataXZ, suffix[11]);
      vtkXY.Add(histo.dataYZ, suffix[12]);

      std::ostringstream sout;
      sout << outdir + "/" + runtag + "_" + OUTFILE + "_surface" + sstr.str();
      vtkXY.Write(sout.str());

    }
  }

  if (PROBE) {
    
    // ==================================================
    // Write line profile along three axes
    // ==================================================
    
    double v;
    float f;
    bool use_log;
    double dR;

    if (RMIN>0.0) {
      use_log = true;
      dR = (log(RMAX) - log(RMIN))/(OUTR-1);
    } else {
      use_log = false;
      dR = RMAX/(OUTR-1);
    }
      
    double r, phi, costh;
    double p0, p1, d0, d1, fr, ft, fp;
    int indx;
    
    vector<double> data(3*nout1*OUTR, 0.0);
    
    for (int l=0; l<OUTR; l++) {
      
      r = dR*l;
      if (use_log) r = RMIN*exp(r);
      
      if ((ncnt++)%numprocs == myid) {
	  
	indx = 3*nout1*l;

	costh = 0.0;
	phi   = 0.0;
	ortho.all_eval(r, costh, phi, d0, d1, p0, p1, fr, ft, fp);
	  
	data[indx + 0] = p0;
	data[indx + 1] = p1;
	data[indx + 2] = p0 + p1;
	data[indx + 3] = fr;
	data[indx + 4] = ft;
	data[indx + 5] = fp;
	data[indx + 6] = d0;
	data[indx + 7] = d1;
	data[indx + 8] = d0 + d1;
	if (d0>0.0)
	  data[indx + 9] = d1/d0;
	else
	  data[indx + 9] = 0.0;

	costh = 0.0;
	phi   = 0.5*M_PI;
	ortho.all_eval(r, costh, phi, d0, d1, p0, p1, fr, ft, fp);
	  
	indx += nout1;
	data[indx + 0] = p0;
	data[indx + 1] = p1;
	data[indx + 2] = p0 + p1;
	data[indx + 3] = fr;
	data[indx + 4] = ft;
	data[indx + 5] = fp;
	data[indx + 6] = d0;
	data[indx + 7] = d1;
	data[indx + 8] = d0 + d1;
	if (d0>0.0)
	  data[indx + 9] = d1/d0;
	else
	  data[indx + 9] = 0.0;

	costh = 1.0;
	phi   = 0.0;
	ortho.all_eval(r, costh, phi, d0, d1, p0, p1, fr, ft, fp);
	  
	indx += nout1;
	data[indx + 0] = p0;
	data[indx + 1] = p1;
	data[indx + 2] = p0 + p1;
	data[indx + 3] = fr;
	data[indx + 4] = ft;
	data[indx + 5] = fp;
	data[indx + 6] = d0;
	data[indx + 7] = d1;
	if (d0>0.0)
	  data[indx + 8] = d1/d0;
	else
	  data[indx + 8] = 0.0;

      }
    }
    
    if (myid==0)
      MPI_Reduce(MPI_IN_PLACE, &data[0], 3*nout1*OUTR, MPI_DOUBLE, MPI_SUM, 
	       0, MPI_COMM_WORLD);
    else
      MPI_Reduce(&data[0], NULL, 3*nout1*OUTR, MPI_DOUBLE, MPI_SUM, 
	       0, MPI_COMM_WORLD);
    
    if (myid==0) {
      
      vector<string> names(nout1);
      for (int i=0; i<nout1; i++) {
	names[i] = outdir + "/" + OUTFILE + "." + suffix[i] + ".cut" + sstr.str();
      }

      foarray out(names, true);

      for (int l=0; l<OUTR; l++) {
	
	r = dR*l;
	if (use_log) r = RMIN*exp(r);
      
	indx = 3*nout1*l;
	
	for (int n=0; n<nout1; n++)
	  out[n] << setw(18) << time << setw(18) << r
		 << setw(18) << data[indx + 0*nout1 + n]
		 << setw(18) << data[indx + 1*nout1 + n]
		 << setw(18) << data[indx + 2*nout1 + n]
		 << endl;
      }
      
      for (int n=0; n<nout1; n++) out[n] << endl;
    }
  }
}


int
main(int argc, char **argv)
{
  
#ifdef DEBUG
  sleep(20);
#endif  
  
  double snr, rscale, Hexp;
  int NICE, LMAX, NMAX, NPART;
  int beg, end, stride, init;
  std::string MODFILE, INDEX, dir("./"), cname, coefs;

  // ==================================================
  // Parse command line or input parameter file
  // ==================================================
  
  std::ostringstream sout;
  sout << std::string(60, '-') << std::endl
       << "Compute halo potential, force and density profiles from" << std::endl
       << "PSP phase-space output files" << std::endl
       << std::string(60, '-') << std::endl << std::endl
       << "Allowed options";
  
  po::options_description desc(sout.str());
  desc.add_options()
    ("help,h",                                                                          "Print this help message")
    ("OUT",
     "assume original, single binary PSP files as input")
    ("SPL",
     "assume new split binary PSP files as input")
    ("CONLY",
     "make coefficient file only")
    ("NICE",                po::value<int>(&NICE)->default_value(0),
     "system priority")
    ("RMIN",                po::value<double>(&RMIN)->default_value(0.0),
     "minimum radius for output")
    ("RMAX",                po::value<double>(&RMAX)->default_value(0.1),
     "maximum radius for output")
    ("RSCALE",              po::value<double>(&rscale)->default_value(0.067),
     "coordinate mapping scale factor")
    ("LMAX",                po::value<int>(&LMAX)->default_value(4),
     "Maximum harmonic order for spherical expansion")
    ("NMAX",                po::value<int>(&NMAX)->default_value(12),
     "Maximum radial order for spherical expansion")
    ("MMAX",                po::value<int>(&MMAX)->default_value(4),
     "Maximum harmonic order")
    ("L1",                  po::value<int>(&L1)->default_value(0),
     "minimum l harmonic")
    ("L2",                  po::value<int>(&L2)->default_value(100),
     "maximum l harmonic")
    ("NPART",               po::value<int>(&NPART)->default_value(0),
     "Jackknife partition number for testing (0 means off, use standard eval)")
    ("Hexp",                po::value<double>(&Hexp)->default_value(1.0),           "default Hall smoothing exponent")
    ("OUTR",                po::value<int>(&OUTR)->default_value(40),
     "Number of radial points for output")
    ("PROBE",               po::value<bool>(&PROBE)->default_value(false),
     "Make traces along axes")
    ("SURFACE",             po::value<bool>(&SURFACE)->default_value(true),
     "Make equitorial and vertical slices")
    ("VOLUME",              po::value<bool>(&VOLUME)->default_value(false),
     "Make volume for VTK")
    ("OUTFILE",             po::value<string>(&OUTFILE)->default_value("haloprof"),
     "Filename prefix")
    ("runtag",              po::value<string>(&runtag)->default_value("run1"),
     "Phase space file")
    ("outdir",              po::value<string>(&outdir)->default_value("."),
     "Output directory path")
    ("MODFILE",             po::value<string>(&MODFILE)->default_value("SLGridSph.model"),
     "Halo model file")
    ("init",                po::value<int>(&init)->default_value(0),
     "fiducial PSP index")
    ("beg",                 po::value<int>(&beg)->default_value(0),
     "initial PSP index")
    ("end",                 po::value<int>(&end)->default_value(99999),
     "final PSP index")
    ("stride",              po::value<int>(&stride)->default_value(1),
     "PSP index stride")
    ("compname",            po::value<std::string>(&cname)->default_value("stars"),
     "train on Component (default=stars)")
    ("dir,d",               po::value<std::string>(&dir),
     "directory for SPL files")
    ("coefs,c",               po::value<std::string>(&coefs),
     "file of computed coefficients")
    ("snr,S",
     po::value<double>(&snr)->default_value(-1.0),
     "if not negative: do a SNR cut on the PCA basis")
    ("diff",
     "render the difference between the trimmed and untrimmed basis")
    ;
  
  
  // ==================================================
  // MPI preliminaries
  // ==================================================

  local_init_mpi(argc, argv);
  
  po::variables_map vm;

  try {
    po::store(po::parse_command_line(argc, argv, desc), vm);
    po::notify(vm);    
  } catch (po::error& e) {
    if (myid==0) std::cout << "Option error: " << e.what() << std::endl;
    exit(-1);
  }

  // ==================================================
  // Print help message and exit
  // ==================================================

  if (vm.count("help")) {
    std::cout << std::endl << desc << std::endl;
    return 0;
  }

  bool rendering = true;
  if (vm.count("coefs") and vm.count("CONLY")) {
    rendering = false;
  }

  bool SPL = false;
  if (vm.count("SPL")) SPL = true;
  if (vm.count("OUT")) SPL = false;

  bool Hall = false;
  if (vm.count("Hall")) Hall = true;


  // ==================================================
  // Nice process
  // ==================================================

  if (NICE>0)
    setpriority(PRIO_PROCESS, 0, NICE);

  // ==================================================
  // Make SL expansion
  // ==================================================

  SphericalModelTable halo(MODFILE);
  SphereSL::mpi  = true;
  SphereSL::NUMR = 4000;
  SphereSL::HEXP = Hexp;
  SphereSL ortho(&halo, LMAX, NMAX, 1, rscale, true, NPART);
  
  std::string file;

  std::ofstream outcoef;	// Coefficient file

  if (vm.count("coefs")) {
    std::string coeffile = outdir + "/" + OUTFILE + ".coefs";
				// Set exceptions to be thrown on failure
    outcoef.exceptions(std::ifstream::failbit | std::ifstream::badbit);

    try {
      outcoef.open(coeffile, std::ofstream::out | std::ofstream::app);
    } catch (std::system_error& e) {
      std::cerr << e.code().message() << std::endl;
    }
  }

  for (int indx=beg; indx<=end; indx+=stride) {

    // ==================================================
    // PSP input stream
    // ==================================================

    int iok = 1;
    std::ostringstream s0, s1;

    if (myid==0) {
      s1.str("");		// Clear stringstream
      if (SPL) s1 << "SPL.";
      else     s1 << "OUT.";
      s1 << runtag << "."<< std::setw(5) << std::setfill('0') << indx;
      
				// Check for existence of next file
      file = dir + s1.str();
      std::ifstream in(file);
      if (!in) {
	cerr << "Error opening <" << file << ">" << endl;
	iok = 0;
      }
    }
    
    MPI_Bcast(&iok, 1, MPI_INT, 0, MPI_COMM_WORLD);
    if (iok==0) break;

    // ==================================================
    // Open PSP file
    // ==================================================
    PSPptr psp;

    if (myid==0) {

      if (SPL) psp = std::make_shared<PSPspl>(s1.str(), dir, true);
      else     psp = std::make_shared<PSPout>(file, true);

      tnow = psp->CurrentTime();
      cout << "Beginning partition [time=" << tnow
	   << ", index=" << indx << "] . . . "  << flush;
    }
    
    Histogram histo(OUTR, RMAX);
    std::vector<Particle> particles;

    partition(psp, cname, particles, histo);
    if (myid==0) cout << "done" << endl;

    //------------------------------------------------------------ 

    if (myid==0) cout << "Accumulating particle positions . . . " << flush;

    ortho.reset_coefs();
    for (auto &i : particles) {
      ortho.accumulate(i.pos[0], i.pos[1], i.pos[2], i.mass);
    }
    MPI_Barrier(MPI_COMM_WORLD);
    if (myid==0) cout << "done" << endl;
  
    //------------------------------------------------------------ 

    if (myid==0) cout << "Making coefficients . . . " << flush;
    ortho.make_coefs();
    MPI_Barrier(MPI_COMM_WORLD);
    if (myid==0) cout << "done" << endl;

    //------------------------------------------------------------ 
    //
    // Coefficient trimming
    //
    if (snr>=0.0) {

      if (myid==0) {
	std::cout << "Computing SNR=" << snr;
	if (Hall) std::cout << " using Hall smoothing . . . " << flush;
	else      std::cout << " using truncation . . . " << flush;
      }
    
      ortho.make_covar();

      // Get the snr trimmed coefficients
      //
      Matrix origc = ortho.retrieve_coefs();
      Matrix coefs = ortho.get_trimmed(snr, Hall);

      if (vm.count("diff")) coefs = coefs - origc;
      ortho.install_coefs(coefs);
    }

    //------------------------------------------------------------ 

    double time = 0.0;
    if (myid==0) {
      time = psp->CurrentTime();
      if (outcoef.good()) {
	cout << "Writing coefficients . . . " << flush;
	try {
	  ortho.dump_coefs(time, outcoef);
	} catch (std::system_error& e) {
	  std::cerr << e.code().message() << std::endl;
	}
	cout << "done" << endl;
      }
    }
    if (rendering) {
      if (myid==0) cout << "Writing output . . . " << flush;
      write_output(ortho, indx, time, histo);
      if (myid==0) cout << "done" << endl;
    }
    MPI_Barrier(MPI_COMM_WORLD);

    //------------------------------------------------------------ 

  } // Dump loop

  MPI_Finalize();

  return 0;
}

