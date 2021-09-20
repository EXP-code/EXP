#ifndef _DiskHalo_H
#define _DiskHalo_H

#include <math.h>
#include <stdlib.h>
#include <string>

#include <iostream>
#include <iomanip>
#include <memory>

#include <boost/random/mersenne_twister.hpp>
#include <boost/random/uniform_real_distribution.hpp>
#include <boost/random/normal_distribution.hpp>

#include <exponential.H>

#include <gaussQ.H>
#include <localmpi.H>

#include <AddDisk.H>
#include <SphericalSL.H>
#include "EmpCylSL.H"
#include "SParticle.H"

class DiskHalo
{
 private:
  std::shared_ptr<AddDisk> newmod;
  std::shared_ptr<SphericalModelTable> halo, halo2;
  std::shared_ptr<ExponentialDisk> disk;
  double scaleheight, hmass, dmass;
  double center_pos[3], center_vel[3];

  std::shared_ptr<SphericalSL> expandh;
  std::shared_ptr<EmpCylSL>    expandd;

  std::vector<Eigen::MatrixXd> disktableP, disktableN;
  Eigen::MatrixXd epitable;
  double dP, dR, dZ, epiRmin;

  Eigen::MatrixXd halotable;
  double dr, dc;

  std::shared_ptr<LegeQuad> lwz;

  boost::mt19937 gen;
  boost::random::uniform_real_distribution<> rndU;
  boost::random::normal_distribution<> rndN;

  bool DF;
  bool com_cov;

  double disk_density(double R, double z);
  void write_record(ofstream &out, SParticle &p);
  void write_record(ofstream &out, Particle &p)
  {
    SParticle P(p);
    write_record(out, P);
  }

 public:
  static int NDP;		// Number of knots in disk table phi grid
				// Default: 16

  static int NDZ;		// Number of knots in disk table z grid
				// Default: 200

  static int NDZF;		// Number of subdivisions for each table
				// entry in Jeans equations evaluation

  static int NDR;		// Number of knots in disk table R grid
				// Default: 40

  static int NHR;		// Number of knots in halo table r grid
				// Default: 40

  static int NHT;		// Number of knots in halo table theta grid
				// Default: 40
  
  static double RHMIN;		// Smallest radius in halo grid

  static double RHMAX;		// Largest radius in halo grid

  static double RDMIN;		// Smallest radius in disk grid

  static double RDMAX;		// Largest radius in disk grid & model

  static double SHFACTOR;	// Minimum vertical integration size
				// in units of scale height

  static double DMFACTOR;	// Increase/decrease the deprojected disk
				// mass by this factor. Default: 1.0

  static double Q;		// Toomre Q

  static double R_DF;		// Change over points for DF and Jeans
  static double DR_DF;

  static int LOGSCALE;		// Log DF in SphericalModelTable

  static unsigned NBUF;		// Number of particles in MPI buffer
				// default: 8192

  static int SEED;		// Initial seed

  //! Constructor
  DiskHalo(std::shared_ptr<SphericalSL> expandh, std::shared_ptr<EmpCylSL> expandd,
	   double dz0, double dh, double dMh, double dMd, std::string& filename,
	   int DF=0, int DIVERGE=0, double DIVERGE_RFAC=1.0);

  //! Destructor
  ~DiskHalo();

  void set_halo_coordinates(vector<Particle>& phalo, int nhalo, int npart);
  void set_disk_coordinates(vector<Particle>& pdisk, int ndisk, int npart);

  void set_pos_origin(double& x, double& y, double& z) 
    { 
      center_pos[0] = x;
      center_pos[1] = y;
      center_pos[2] = z;
    }

  void set_vel_origin(double& x, double& y, double& z) 
    { 
      center_vel[0] = x;
      center_vel[1] = y;
      center_vel[2] = z;
    }

  double get_hpot(double xp, double yp, double zp);

  double get_dpot(double xp, double yp, double zp);

  double get_hdpot(double xp, double yp, double zp);

  double get_ddpot(double xp, double yp, double zp);

  double deri_pot(double xp, double yp, double zp, int n);

  double deri_pot_disk(double xp, double yp, double zp, int n);

  double deri_pot_halo(double xp, double yp, double zp, int n);

  void disk_eval(double R, double z, double phi,
		 double &p, double &fr, double &fz, double &fp);

  double epi(double xp, double yp, double zp);

  void table_disk(vector<Particle>& part);

  double get_dispdz(double xp,double yp,double zp);

  double vr_disp(double xp, double yp,double zp);

  double vphi(double xp, double yp, double zp);

  double v_circ(double xp, double yp, double zp);

  void set_vel_disk(vector<Particle>& part);

  void table_halo(vector<Particle>& part);

  double get_disp(double xp, double yp,double zp);

  void set_vel_halo(vector<Particle>& part);

  void write_file(ofstream &fou_halo, ofstream &fou_disk,
		  vector<Particle>& hpart, vector<Particle>& dpart);

  void virial_ratio(vector<Particle>& hpart, vector<Particle>& dpart);

  void virial_ratio(const char *, const char *);

  void zero_com_cov(bool val) { com_cov = val; }

};

#endif

