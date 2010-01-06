				// System
#include <values.h>
				// C++/STL
#include <iostream>
#include <fstream>
#include <sstream>
#include <algorithm>
#include <vector>
				// MDW
#include <interp.h>
#include <numerical.h>
#include <exponential.h>
#include <Vector.h>
#include <interp.h>
				// Local
#include "AddDisk.h"
#include "DiskHalo2.h"

				// Grid parameters and Toomre Q
double DiskHalo::RHMIN = 1.0e-4;
double DiskHalo::RHMAX = 50.0;
double DiskHalo::RDMIN = 1.0e-4;
double DiskHalo::RDMAX = 20.0;
double DiskHalo::Q = 1.2;
double DiskHalo::SHFACTOR = 16.0;
double DiskHalo::TOLE = 0.003;
double DiskHalo::COMPRESSION = 1.0;
int    DiskHalo::NDP = 16;
int    DiskHalo::NDZ = 40;
int    DiskHalo::NDR = 800;
int    DiskHalo::NHR = 800;
int    DiskHalo::NHT = 40;
int    DiskHalo::SEED = 11;

double DiskHalo::RA = 1.0e20;
int    DiskHalo::NUMDF = 800;
int    DiskHalo::RNUM = 4000;

double DiskHalo::R_DF = 20.0;
double DiskHalo::DR_DF = 5.0;

int    DiskHalo::LOGSCALE = 0;
bool   DiskHalo::LOGR = true;

bool   DiskHalo::JEANS_VR = true;

unsigned DiskHalo::VFLAG = 0;
unsigned DiskHalo::NBUF = 8192;

static AxiSymModel *model;
double targetmass;
				// Determine radius with given enclosed mass
double mass_func(double r)
{
  return targetmass - model->get_mass(r);
}

DiskHalo::
DiskHalo()
{
  disktableP = NULL;
  disktableN = NULL;
  gen = NULL;
  rndU = NULL;
  rndN = NULL;
  com = false;
  cov = false;
  DF = false;
  MULTI = false;
  halo = NULL;
  halo2 = NULL;
  disk = NULL;
}

DiskHalo::
DiskHalo(SphericalSL* haloexp, EmpCylSL* diskexp,
	 double H, double A, double DMass, 
	 string& filename, int DF1, int DIVERGE, double DIVERGE_RFAC)
{
  disktableP = NULL;
  disktableN = NULL;
  gen = new ACG(SEED+myid, 20);
  rndU = new Uniform(0.0, 1.0, gen);
  rndN = new Normal(0.0, 1.0, gen);
  com = false;
  cov = false;

  for (int k=0; k<3; k++) center_pos[k] = center_vel[k] = 0.0;

  if (DF1) DF = true;
  else     DF = false;

  MULTI = false;

  dmass = DMass;
  scaleheight = H;

  expandh = haloexp;
  expandd = diskexp;

  SphericalModelTable::even = 0;
  SphericalModelTable::logscale = LOGSCALE;

  halo = new SphericalModelTable(filename, DIVERGE, DIVERGE_RFAC);

  disk = new ExponentialDisk(A, RDMAX);

  if (myid==0 && VFLAG & 1) {
    cerr << "DiskHalo: DIVERGE=" << DIVERGE
	 << " A=" << A
	 << " RDMAX=" << RDMAX
	 << " filename=" << filename
	 << endl;
  }

  if (DF) {
    AddDisk::logarithmic = true;
    AddDisk::Rmin = RHMIN;
    AxiSymModel::numr = 400;
    AxiSymModel::numj = 400;
    AxiSymModel::gen_N = 800;
    AxiSymModel::gen_itmax = 400000;
    AxiSymModel::gen_rmin = RHMIN;
    newmod = new AddDisk(halo, disk, dmass*COMPRESSION); 
    halo2 = newmod->get_model();
    halo2->setup_df(NUMDF, RA);
    if (myid==0 && VFLAG & 2) {
      char debugname[] = "df.debug";
      halo2->print_df(debugname);
    }
  } else {
    halo2 = halo;
  }
}

DiskHalo::
DiskHalo(SphericalSL* haloexp, EmpCylSL* diskexp,
	 double H, double A, double DMass, 
	 string& filename1, int DIVERGE, double DIVERGE_RFAC,
	 string& filename2, int DIVERGE2, double DIVERGE_RFAC2)
{
  disktableP = NULL;
  disktableN = NULL;
  gen = new ACG(SEED+myid, 20);
  rndU = new Uniform(0.0, 1.0, gen);
  rndN = new Normal(0.0, 1.0, gen);
  com = false;
  cov = false;

  for (int k=0; k<3; k++) center_pos[k] = center_vel[k] = 0.0;

  DF = false;
  MULTI = true;

  dmass = DMass;
  scaleheight = H;

  expandh = haloexp;
  expandd = diskexp;

  SphericalModelTable::even = 0;
  SphericalModelTable::logscale = LOGSCALE;

  halo = new SphericalModelTable(filename1, DIVERGE, DIVERGE_RFAC);

  disk = new ExponentialDisk(A, RDMAX);

  if (myid==0 && VFLAG & 1) {
    cerr << "DiskHalo: DIVERGE=" << DIVERGE
	 << " DIVERGE2=" << DIVERGE2
	 << " A=" << A
	 << " RDMAX=" << RDMAX
	 << " filename1=" << filename1
	 << " filename2=" << filename2
	 << "\n";
  }

  AddDisk::logarithmic = true;
  AddDisk::Rmin = RHMIN;
  AxiSymModel::numr = 400;
  AxiSymModel::numj = 400;
  AxiSymModel::gen_N = 800;
  AxiSymModel::gen_itmax = 4000000;
  AxiSymModel::gen_rmin = RHMIN;
  newmod = new AddDisk(halo, disk, dmass*COMPRESSION); 
  halo2 = newmod->get_model();
  halo2->setup_df(NUMDF, RA);
  if (myid==0 && VFLAG & 2) {
    char debugname[] = "df.debug";
    halo2->print_df(debugname);
  }

  //
  // Generate "fake" profile
  //
  SphericalModelTable::even = 0;
  SphericalModelTable::logscale = LOGSCALE;
  SphericalModelTable::linear = 0;
  
  halo3 = new SphericalModelTable(filename2, DIVERGE2, DIVERGE_RFAC2);

  //
  // Packs fake density and mass model with target (real) potential
  // and reinitializes the model
  //
  double *r2 = new double [RNUM];
  double *d2 = new double [RNUM];
  double *m2 = new double [RNUM];
  double *p2 = new double [RNUM];
  
  double rmin = halo2->get_min_radius();
  double rmax = halo2->get_max_radius();
  double r, dr = (log(rmax) - log(rmin))/(RNUM-1);

  if (LOGR)
    dr = (log(rmax) - log(rmin))/(RNUM-1);
  else
    dr = (rmax - rmin)/(RNUM-1);
  
  for (int i=0; i<RNUM; i++) {
    if (LOGR)
      r2[i] = r = rmin*exp(dr*i);
    else
      r2[i] = r = rmin + dr*i;
    d2[i] = halo3 -> get_density(r);
    m2[i] = halo3 -> get_mass(r);
    p2[i] = halo2 -> get_pot(r);
  }
  
  delete halo3;

  halo3 = new SphericalModelTable(RNUM, r2-1, d2-1, m2-1, p2-1, 
				  DIVERGE2, DIVERGE_RFAC2);
  halo3->setup_df(NUMDF, RA);
  if (VFLAG & 2) {
    halo3->print_model("diskhalo2_model.multi");
    halo3->print_df("diskhalo2_df.multi");
  }
    
  delete [] r2;
  delete [] d2;
  delete [] m2;
  delete [] p2;
  
  //
  // Generate the multimass model
  //
  multi = new SphericalModelMulti(halo2, halo3);
  multi -> gen_tolE = TOLE;

}


DiskHalo::~DiskHalo()
{
  delete [] disktableP;
  delete [] disktableN;
  delete rndU;
  delete rndN;
  delete gen;
  delete halo;
  delete disk;
  if (DF) {
    delete newmod;
  }
  if (MULTI) {
    delete newmod;
    delete halo3;
    delete multi;
  }
}

DiskHalo::DiskHalo(const DiskHalo &p)
{
  newmod = p.newmod;
  halo   = p.halo;
  halo2  = p.halo2;
  halo3  = p.halo3;
  multi  = p.multi;
  disk   = p.disk;

  scaleheight = p.scaleheight;
  dmass = p.dmass;

  for (int k=0; k<3; k++) {
    center_pos[k] = p.center_pos[k];
    center_vel[k] = p.center_vel[k];
  }

  expandh = p.expandh;
  expandd = p.expandd;

  disktableP = p.disktableP;
  disktableN = p.disktableN;
  epitable = p.epitable;
  dv2table = p.dv2table;

  dP = p.dP;
  dR = p.dR;
  dZ = p.dZ;
  epiRmin = p.epiRmin;

  halotable = p.halotable;
  dr = p.dr;
  dc = p.dc;

  gen = new ACG(SEED+myid, 20);
  rndU = new Uniform(0.0, 1.0, gen);
  rndN = new Normal(0.0, 1.0, gen);

  DF = p.DF;
  MULTI = p.MULTI;
  com = p.com;
  cov = p.cov;
}


double DiskHalo::disk_density(double R, double z)
{
  double q = 1.0/cosh(z/scaleheight);
  return dmass*disk->get_density(R)*q*q*0.5/scaleheight;
}

void DiskHalo::set_halo(vector<Particle>& phalo, int nhalo, int npart)
{
  if (!MULTI) {
    string msg("DiskHalo::set_halo is only valid if MULTI is true");
    throw DiskHaloException(msg, __FILE__, __LINE__);
  }

  double rmin = max<double>(halo2->get_min_radius(), RHMIN);
  double rmax = halo2->get_max_radius();
  double mmin = halo2->get_mass(rmin);
  double mtot = halo2->get_mass(rmax);

  double pos[3], pos1[3], vel[3], vel1[3], massp, massp1;
				// Diagnostics
  double radmin1=1e30, radmax1=0.0, radmin, radmax, r;

  if (myid==0) cout << endl 
		    << "     *****"
		    << "  rmin=" << rmin
		    << "  rmax=" << rmax
		    << "  mmin=" << mmin
		    << "  mtot=" << mtot
		    << endl;

  for (int k=0; k<3; k++) {
    pos[k] = pos1[k] = 0.0;
    vel[k] = vel1[k] = 0.0;
  }
  massp = massp1 = 0.0;

  double meanmass = (mtot - mmin)/nhalo;

  Particle p;
  Vector ps(0, 6);
  int ierr;

  int count1=0, count=0;

  for (int i=0; i<npart; i++) {

    do {
      ps = multi->gen_point(ierr);
      if (ierr) count1++;
    } while (ierr);
    
    p.mass = meanmass * ps[0];
    
    for (int i=1; i<=3; i++) {
      p.pos[i-1] = ps[i];
      p.vel[i-1] = ps[i+3];
    }

    massp1 += p.mass;
    for (int k=0; k<3; k++) pos1[k] += p.mass*p.pos[k];
    for (int k=0; k<3; k++) vel1[k] += p.mass*p.vel[k];
    
    phalo.push_back(p);

    r = 0.0;
    for (int k=0; k<3; k++) r += p.pos[k]*p.pos[k];
    r = sqrt(r);

    radmin1 = min<double>(radmin1, r);
    radmax1 = max<double>(radmax1, r);
  }
  
  MPI_Reduce(&count1, &count, 1, MPI_INT, MPI_SUM, 0, MPI_COMM_WORLD);
  MPI_Reduce(&radmin1, &radmin, 1, MPI_DOUBLE, MPI_MIN, 0, MPI_COMM_WORLD);
  MPI_Reduce(&radmax1, &radmax, 1, MPI_DOUBLE, MPI_MAX, 0, MPI_COMM_WORLD);

  if (myid==0) cout << "     *****"
		    << "  min(r)=" << radmin 
		    << "  max(r)=" << radmax
		    << endl;

  if (myid==0 && count) cout << "DiskHalo::set_halo: " 
			     << count << " selection failures" << endl;
  
  MPI_Allreduce(&massp1, &massp, 1, MPI_DOUBLE, MPI_SUM, MPI_COMM_WORLD);
  MPI_Allreduce(pos1, pos, 3, MPI_DOUBLE, MPI_SUM, MPI_COMM_WORLD);
  MPI_Allreduce(vel1, vel, 3, MPI_DOUBLE, MPI_SUM, MPI_COMM_WORLD);

  if (massp>0.0) {
    for (int k=0; k<3; k++) pos[k] /= massp;
    for (int k=0; k<3; k++) vel[k] /= massp;
  }
  
  if (myid==0) {
    cout << "     *****";
    cout << " (x, y, z)=(" << pos[0] 
	 << ", " << pos[1]
	 << ", " << pos[2] << ")" << endl;
    cout << "     *****";
    cout << " (u, v, w)=(" << vel[0]
	 << ", " << vel[1]
	 << ", " << vel[2] << ")" << endl;
  }

  if (com || cov) {
    vector<Particle>::iterator ip;
    for (ip=phalo.begin(); ip!=phalo.end(); ip++) {
      if (com) for (int k=0; k<3; k++) ip->pos[k] -= pos[k];
      if (cov) for (int k=0; k<3; k++) ip->vel[k] -= vel[k];
    }
  }

  // Make dispersion vector
  //
  table_halo_disp();

  if (VFLAG & 1)
    cout << "Process " << myid << ": made " << phalo.size() << " particles"
	 << endl;
}      

void DiskHalo::
set_halo_coordinates(vector<Particle>& phalo, int nhalo, int npart)
{
  const double tol = 1.0e-8;
  double rmin = max<double>(halo2->get_min_radius(), RHMIN);
  double rmax = halo2->get_max_radius();
  double mmin = halo2->get_mass(rmin);
  double mtot = halo2->get_mass(rmax);

  double r, phi, costh, sinth;
  double massp, massp1, pos[3], pos1[3];
				// Diagnostics
  double radmin1=1.0e30, radmax1=0.0, radmin, radmax;

  if (myid==0 && VFLAG & 1) cout << "  rmin=" << rmin
				 << "  rmax=" << rmax
				 << "  mmin=" << mmin
				 << "  mtot=" << mtot;

  for (int k=0; k<3; k++) pos[k] = pos1[k] = 0.0;
  massp = massp1 = 0.0;

  Particle p;

  p.mass = (mtot - mmin)/nhalo;

  model = halo2;

  MPI_Barrier(MPI_COMM_WORLD);

  for (int i=0; i<npart; i++) {
    targetmass = mmin + (mtot - mmin)*(*rndU)();

    r = zbrent(mass_func, rmin, rmax, tol);
    
    phi = 2.0*M_PI*(*rndU)();
    costh = 2.0*(*rndU)() - 1.0;
    sinth = sqrt(1.0 - costh*costh);

    p.pos[0] = r*sinth*cos(phi);
    p.pos[1] = r*sinth*sin(phi);
    p.pos[2] = r*costh;

    massp1 += p.mass;
    for (int k=0; k<3; k++) pos1[k] += p.mass*p.pos[k];

    phalo.push_back(p);

    r = 0.0;
    for (int k=0; k<3; k++) r += p.pos[k] * p.pos[k];
    r = sqrt(r);

    radmin1 = min<double>(radmin1, r);
    radmax1 = max<double>(radmax1, r);
  }

  MPI_Reduce(&radmin1, &radmin, 1, MPI_DOUBLE, MPI_MIN, 0, MPI_COMM_WORLD);
  MPI_Reduce(&radmax1, &radmax, 1, MPI_DOUBLE, MPI_MAX, 0, MPI_COMM_WORLD);
  
  if (myid==0) cout << "  min(r)=" << radmin 
		    << "  max(r)=" << radmax;

  MPI_Allreduce(&massp1, &massp, 1, MPI_DOUBLE, MPI_SUM, MPI_COMM_WORLD);
  MPI_Allreduce(pos1, pos, 3, MPI_DOUBLE, MPI_SUM, MPI_COMM_WORLD);

  if (massp>0.0) {
    for (int k=0; k<3; k++) pos[k] /= massp;
  }

  if (myid==0) {
    cout << " (x, y, z)=(" << pos[0] 
	 << ", " << pos[1]
	 << ", " << pos[2] << ")" << endl;
  }
  

  if (com) {
    vector<Particle>::iterator ip;
    for(ip=phalo.begin(); ip!=phalo.end(); ip++) {
      for (int k=0; k<3; k++) ip->pos[k] -= pos[k];
    }
  }

  if (VFLAG & 1)
    cout << "Process " << myid << ": made " << phalo.size() << " particles"
	 << endl;
}      

void DiskHalo::
set_disk_coordinates(vector<Particle>& pdisk, int ndisk, int npart)
{
  const double tol = 1.0e-8;
  double rmin = max<double>(disk->get_min_radius(), RDMIN);
  double rmax = disk->get_max_radius();
  double mmin = disk->get_mass(rmin);
  double mtot = disk->get_mass(rmax);

  double R, phi;
  double pos[3], pos1[3], massp, massp1;

  // Diagnostics
  //
  double radmin1=1.0e30, radmax1=0.0, radmin, radmax, r;

  if (myid==0) cout << endl
		    << "     *****"
		    << "  rmin=" << rmin
		    << "  rmax=" << rmax
		    << "  mmin=" << mmin
		    << "  mtot=" << mtot
		    << endl;

  for (int k=0; k<3; k++) pos[k] = pos1[k] = 0.0;
  massp = massp1 = 0.0;

  Particle p;

  p.mass = dmass/ndisk;

  model = disk;

  for (int i=0; i<npart; i++) {
    targetmass = mmin + (mtot-mmin)*(*rndU)();
    R = zbrent(mass_func, rmin, rmax, tol);
    phi = 2.0*M_PI*(*rndU)();

    p.pos[0] = R*cos(phi);
    p.pos[1] = R*sin(phi);
    p.pos[2] = scaleheight*atanh(2.0*(*rndU)()-1.0);

    massp1 += p.mass;
    for (int k=0; k<3; k++) pos1[k] += p.mass*p.pos[k];
    
    pdisk.push_back(p);

    r = 0.0;
    for (int k=0; k<3; k++) r += p.pos[k] * p.pos[k];
    r = sqrt(r);

    radmin1 = min<double>(radmin1, r);
    radmax1 = max<double>(radmax1, r);
  }


  MPI_Reduce(&radmin1, &radmin, 1, MPI_DOUBLE, MPI_MIN, 0, MPI_COMM_WORLD);
  MPI_Reduce(&radmax1, &radmax, 1, MPI_DOUBLE, MPI_MAX, 0, MPI_COMM_WORLD);

  if (myid==0) cout << "     *****"
		    << "  min(r)=" << radmin 
		    << "  max(r)=" << radmax 
		    << endl;

  MPI_Allreduce(&massp1, &massp, 1, MPI_DOUBLE, MPI_SUM, MPI_COMM_WORLD);
  MPI_Allreduce(pos1, pos, 3, MPI_DOUBLE, MPI_SUM, MPI_COMM_WORLD);

  if (massp>0.0) {
    for (int k=0; k<3; k++) pos[k] /= massp;
  }

  if (myid==0) {
    cout << "     *****";
    cout << " (x, y, z)=(" << pos[0] 
	 << ", " << pos[1]
	 << ", " << pos[2] << ")" << endl;
  }
  
  if (com) {
    vector<Particle>::iterator ip;
    for(ip=pdisk.begin(); ip!=pdisk.end(); ip++) {
      for (int k=0; k<3; k++) ip->pos[k] -= pos[k];
    }
  }

}      

double DiskHalo::
get_hpot(double xp, double yp, double zp)
{
  if (!expandh) return 0.0;

  double r, phi, theta, dens, potl, potr, pott, potp;
  
  r = sqrt(xp*xp + yp*yp + zp*zp);
  phi = atan2(yp, xp);
  theta = acos(zp/(r+1.0e-14));

  expandh->determine_fields_at_point(r, theta, phi,
				    &dens, &potl, &potr, &pott, &potp);

  return potl;
}

void DiskHalo::disk_eval(double R, double z, double phi,
			 double &p, double &fr, double &fz, double &fp)
{
  //                   This is the table radius
  //                   ------------------------
  if (sqrt(R*R+z*z) <= M_SQRT1_2*expandd->get_ascale()*expandd->RMAX) {
    double p0;
    expandd->accumulated_eval(R, z, phi, p0, p, fr, fz, fp);

  } else {

    double fac, r2 = R*R + z*z;
    p = -dmass/sqrt(r2);	// -M/r
    fac = p/r2;			// -M/r^3
    fr = fac * R;
    fz = fac * z;
    fp = 0.0;
  }
}

double DiskHalo::
get_dpot(double xp, double yp, double zp)
{
  if (!expandd) return 0.0;

  double r, phi=0.0, pot=0.0, fr=0.0, fz=0.0, fp=0.0;
  
  r = sqrt(xp*xp + yp*yp);
  disk_eval(r, zp, phi, pot, fr, fz, fp);
  
  return pot;
}

double DiskHalo::
get_hdpot(double xp, double yp, double zp)
{
  if (!expandh) return 0.0;

  double r, theta, phi, dens, potl, potr, pott, potp;
  
  r = sqrt(xp*xp + yp*yp);
  theta = 0.5*M_PI;
  phi = atan2(yp, xp);
  expandh->determine_fields_at_point(r, theta, phi, 
				    &dens, &potl, &potr, &pott, &potp);
  
  return potr;
}

double DiskHalo::
get_ddpot(double xp, double yp, double zp)
{
  if (!expandd) return 0.0;

  double r = sqrt(xp*xp + yp*yp);
  double phi = atan2(yp, xp);
  double p, fr, fz, fp;

  disk_eval(r, 0.0, phi, p, fr, fz, fp);

  return -fr;
}

double DiskHalo::
deri_pot(double xp, double yp, double zp, int n)
{
  double dP=0.0;
  
  switch(n){
  case 0:
    dP = deri_pot_disk(xp,yp,zp,0) + deri_pot_halo(xp,yp,zp,0);
    break;
  case 1:
    dP = deri_pot_disk(xp,yp,zp,1) + deri_pot_halo(xp,yp,zp,1);
    break;
  default:
    cout << "deri_pot: derivative of order " << n 
	 << " can't be calculated" << endl;
    break;
  }
  return dP;
}


double DiskHalo::
deri_pot_disk(double xp, double yp, double zp, int n)
{
  double dP=0.0;

  switch(n) {
  case 0:
    dP = get_dpot(xp,yp,zp);
    break;
  case 1:
    dP = get_ddpot(xp,yp,zp);
    break;
  default:
    cout<< "deri_pot_disk: derivative of order " << n 
	<< " can't be calculated" << endl;
    break;
  }
  return dP;
}

double DiskHalo::
deri_pot_halo(double xp, double yp, double zp, int n)
{
  double dP=0.0;
  
  switch(n) {
  case 0:
    dP = get_hpot(xp,yp,zp);
    break;
  case 1:
    dP = get_hdpot(xp,yp,zp);
    break;
  default:
    cout<<"deri_pot_halo: derivative of order " << n 
	<< " can't be calculated" << endl;
    break;
  }
  return dP;
}

double DiskHalo::
epi(double xp, double yp, double zp)
{
  if (disktableP == NULL) {
    cerr << "DiskHalo::get_dispdz: must call table_disk first\n";
    MPI_Abort(MPI_COMM_WORLD, 100);
    exit(0);
  }

  double lR, phi, cp[2], cr[2], ans;
  int iphi1, iphi2, ir;

				// Azimuth
  phi = atan2(yp, xp);
  if (phi<0.0) phi = 2.0*M_PI + phi;

  iphi1 = (int)( phi/dP );
  if (iphi1==NDP-1) iphi2 = 0;	// Modulo 2Pi
  else iphi2 = iphi1+1;
  cp[1] = (phi - dP*iphi1)/dP;
  cp[0] = 1.0 - cp[1];
				// Cylindrical radius
  lR = log(max<double>(epiRmin, sqrt(xp*xp + yp*yp)));
  ir = (int)( (lR - log(RDMIN))/dR );
  ir = min<int>( ir, NDR-2 );
  ir = max<int>( ir, 0 );
  cr[1] = (lR - log(RDMIN) - dR*ir)/dR;
  cr[0] = 1.0 - cr[1];

  ans = 
    cp[0]*cr[0]* epitable[iphi1][ir  ] +
    cp[0]*cr[1]* epitable[iphi1][ir+1] +
    cp[1]*cr[0]* epitable[iphi2][ir  ] +
    cp[1]*cr[1]* epitable[iphi2][ir+1] ;
    
  if (ans>0.0) return sqrt(ans);
  else {
    cout << "Process " << myid << " epi range error" << endl
	 << "     R=" << sqrt(xp*xp + yp*yp)    << endl
	 << "   Phi=" << phi                    << endl
	 << "    cp=" << cp[0] << ", " << cp[1] << endl
	 << "    cr=" << cr[0] << ", " << cr[1] << endl
	 << "    ans=" << ans                   << endl
	 << "    ep1=" << epitable[iphi1][ir  ] << endl
	 << "    ep2=" << epitable[iphi1][ir+1] << endl
	 << "    ep3=" << epitable[iphi2][ir  ] << endl
	 << "    ep4=" << epitable[iphi2][ir+1] << endl
	 << endl;
    return 1.0e-8;
  }
}

// For the disk, setting vz_disp according to the Jeans' equation
// solution, B+T equation 4-29c, where here the potential is the 
// composite basis representation of the total particle distribution. 

void DiskHalo::
table_disk(vector<Particle>& part)
{
  if (disktableP != NULL) return;

  disktableP = new Matrix [NDP];
  disktableN = new Matrix [NDP];
  for (int i=0; i<NDP; i++) {
    disktableP[i].setsize(0, NDR-1, 0, NDZ-1);
    disktableN[i].setsize(0, NDR-1, 0, NDZ-1);
  }
    
  epitable.setsize(0, NDP-1, 0, NDR-1);
  if (JEANS_VR) dv2table.setsize(0, NDP-1, 0, NDR-1);

  dP = 2.0*M_PI/NDP;
  double maxr1 = 0.0, maxz1 = 0.0;
  double maxr =  0.0, maxz  = 0.0;
  vector<Particle>::iterator p;
  for (p=part.begin(); p != part.end(); p++) {
    maxr1 = max<double>(maxr1, sqrt(p->pos[0]*p->pos[0] + p->pos[1]*p->pos[1]));
    maxz1 = max<double>(maxz1, fabs(p->pos[2]));
  }

  maxz1 = max<double>(scaleheight*SHFACTOR, maxz1);

  MPI_Allreduce(&maxr1, &maxr, 1, MPI_DOUBLE, MPI_MAX, MPI_COMM_WORLD);
  MPI_Allreduce(&maxz1, &maxz, 1, MPI_DOUBLE, MPI_MAX, MPI_COMM_WORLD);

  dR = (log(maxr) - log(RDMIN))/(NDR-1);
  dZ = maxz/(NDZ-1);

  if (myid==0) {
    cout << endl
	 << "Table disk epi parametersZ:" << endl
	 << "  RDMIN=" << RDMIN           << endl
	 << "   maxr=" << maxr            << endl
	 << "   maxz=" << maxz            << endl
	 << "   dR=" << dR                << endl
	 << "   dZ=" << dZ                << endl
	 << endl;
  }

  double dZ1 = maxz/(NDZ-1);

  double R, r, x, y, z, phi, fr, fz, fp, theta;
  double pot, dens, potl, dpr, dpt, dpp, dpdz;

				// Add no force if no component exists
  pot = fr = fz = fp = dens = potl = dpr = dpt = dpp = 0.0;

  Vector workP(0, NDZ-1);
  Vector workN(0, NDZ-1);
  Vector workA(0, NDZ-1);
  Vector workZ(0, NDZ-1);

  Vector workR(0, NDR-1);
  Vector workE(0, NDR-1);
  Vector workE2(0, NDR-1);

  Matrix workV(0, 3, 0, NDR-1);

				// For debugging
  Matrix workD(0, 4, 0, NDR-1);
  Matrix workDZ(0, 6, 0, NDZ-1);

				// Compute this table in parallel

  vector<int> ibeg(numprocs), iend(numprocs);
  for (int i=0; i<numprocs; i++) {
    ibeg[i] = (i  )*NDP/numprocs;
    iend[i] = (i+1)*NDP/numprocs;
  }

  for (int i=ibeg[myid]; i<iend[myid]; i++) {

    phi = dP*i;

    for (int j=0; j<NDR; j++) {

      R = RDMIN*exp(dR*j);
      x = R*cos(phi);
      y = R*sin(phi);
      workR[j] = log(RDMIN) + dR*j;

				// For epicylic frequency
				// 
      disk_eval(R, 0.0, phi, pot, fr, fz, fp);
				// Use monopole part of expansion here, only
      if (expandh)		//
	expandh->determine_fields_at_point(R, 0.5*M_PI, phi,
					   &dens, &potl, &dpr, &dpt, &dpp);
      workE[j] = -fr + dpr;
      if (JEANS_VR) {
	workV[0][j] = R;
	workV[1][j] = disk_density(R, 0.0);
	workV[2][j] = workV[1][j]*workE[j];
      }
      if (i==0) {
	workD[3][j] = -fr;
	workD[4][j] = dpr;
      }

      for (int k=0; k<NDZ; k++) {

	z = workZ[k] = dZ1*k;

	r = sqrt(R*R + z*z) + MINDOUBLE;

				// Positive

	disk_eval(R, z, phi, pot, fr, fz, fp);
	theta = acos(z/(r + MINDOUBLE));

	if (expandh)
	  expandh->determine_fields_at_point(R, theta, phi,
					     &dens, &potl, &dpr, &dpt, &dpp);

	dpdz = -fz + dpr*z/r + dpt*R*R/(r*r*r);

				// Debugging
	workDZ[0][k] = disk_density(R, z);
	workDZ[1][k] = -fz;
	workDZ[2][k] = dpr*z/r;
	workDZ[3][k] = dpt*R*R/(r*r*r);

	workP[k] = disk_density(R, z) * dpdz;

				// Negative
	z *= -1.0;

	disk_eval(R, z, phi, pot, fr, fz, fp);
	theta = acos(z/(r + MINDOUBLE));

	if (expandh)
	  expandh->determine_fields_at_point(R, theta, phi,
					     &dens, &potl, &dpr, &dpt, &dpp);

	dpdz = -fz + dpr*z/r + dpt*R*R/(r*r*r);
	dpdz *= -1.0;		// No harm in changing sign

	workDZ[4][k] = -fz;
	workDZ[5][k] = dpr*z/r;
	workDZ[6][k] = dpt*R*R/(r*r*r);

	workN[k] = disk_density(R, z) * dpdz;
      }

				// Integrate positive
      // Splsum(workZ, workP, workA);
      Trapsum(workZ, workP, workA);
      for (int k=0; k<NDZ; k++)
	disktableP[i][j][k] = max<double>(workA[NDZ-1] - workA[k], MINDOUBLE);

      if (i==ibeg[myid] && j==0) {
	ostringstream ofile;
	ofile << "intgr_disk_P.d" << myid;
	ofstream dout(ofile.str().c_str());
	for (int k=0; k<NDZ; k++) {
	  dout << setw(15) << workZ[k] 
	       << setw(15) << workP[k] 
	       << setw(15) << workA[k];
	  for (int q=0; q<7; q++) dout << setw(15) << workDZ[q][k];
	  dout << endl;
	}
	dout.close();
      }
				// Integrate negative
      // Splsum(workZ, workN, workA);
      Trapsum(workZ, workN, workA);
      for (int k=0; k<NDZ; k++)
	disktableN[i][j][k] = max<double>(workA[NDZ-1] - workA[k], MINDOUBLE);

      if (i==ibeg[myid] && j==0) {
	ostringstream ofile;
	ofile << "intgr_disk_N.d" << myid;
	ofstream dout(ofile.str().c_str());
	for (int k=0; k<NDZ; k++) 
	  dout << setw(15) << workZ[k] 
	       << setw(15) << workN[k] 
	       << setw(15) << workA[k]
	       << "\n";
	dout.close();
      }

    }

				// Compute epicylic freqs
    for (int j=0; j<NDR; j++) {
      epitable[i][j] = drv2(workR[j], workR, workE);
      if (i==0) workD[0][j] = epitable[0][j];
      epitable[i][j] += 3.0*workE[j];
      if (i==0) workD[1][j] = epitable[0][j];
      epitable[i][j] /= exp(workR[j]);
      if (i==0) workD[2][j] = epitable[0][j];
    }
    
    if (JEANS_VR) {
      Trapsum(workV[0], workV[2], workV[3]);
      for (int j=0; j<NDR; j++) 
	dv2table[i][j] = (workV[3][NDR-1] - workV[3][j])/workV[1][j];
    }

  }

				// Update tables on all nodes
  for (int k=0; k<numprocs; k++) {
    for (int i=ibeg[k]; i<iend[k]; i++) {
      MPI_Bcast(&epitable[i][0], NDR, MPI_DOUBLE, k, MPI_COMM_WORLD);
      if (JEANS_VR)
	MPI_Bcast(&dv2table[i][0], NDR, MPI_DOUBLE, k, MPI_COMM_WORLD);
      for (int j=0; j<NDR; j++) {
	MPI_Bcast(&disktableP[i][j][0], NDZ, MPI_DOUBLE, k, MPI_COMM_WORLD);
	MPI_Bcast(&disktableN[i][j][0], NDZ, MPI_DOUBLE, k, MPI_COMM_WORLD);
      }
    }
  }

				// Check epitable
  double epirmin=0.0;
  for (int i=0; i<NDP; i++) {
    for (int j=0; j<NDR; j++) {
      
      if (epitable[i][j] < 0.0) {
	epirmin = max<double>(epirmin, RDMIN*exp(dR*j));

	if (myid==0) {
	  cout << "Epitable error: R=" << RDMIN*exp(dR*j)
	       << " Phi=" << dP*i << " ep2=" << epitable[i][j] << endl;
	}
      }
    }
  }

				// Check solution
  if (myid==0 && expandh) {
    ofstream out("ep_test.dat");
    out.setf(ios::scientific);
    out.precision(4);

    const double DR = 0.01;

    double r, r1, r2, deriv, rho, lhs, rhs;

    for (int j=0; j<NDR; j++) {
      r = RDMIN*exp(dR*j);
      r1 = r*(1.0 + dR*DR);
      r2 = r*(1.0 - dR*DR);

      rho = halo->get_density(r);

      deriv = (get_disp(0.0, r1, 0.0)*halo->get_density(r1) - 
	       get_disp(0.0, r2, 0.0)*halo->get_density(r2) ) /	(r1 - r2);
      
      lhs = halo->get_mass(r);
      rhs = -r*r*deriv/rho;

      out << setw(14) << r			// #1
	  << setw(14) << epitable[0][j]		// #2
	  << setw(14) << workR[j]		// #3
	  << setw(14) << workE[j]		// #4
	  << setw(14) << drv2(workR[j], workR, workE) // #5
	  << setw(14) << workD[0][j]		// #6
	  << setw(14) << workD[1][j]		// #7
	  << setw(14) << workD[2][j]		// #8
	  << setw(14) << workD[3][j]		// #9
	  << setw(14) << workD[4][j]		// #10
	  << setw(14) << rho			// #11
	  << setw(14) << deriv			// #12
	  << setw(14) << lhs			// #13
	  << setw(14) << rhs			// #14
	  << setw(14) << lhs - rhs		// #15
	  << endl;
    }

    ofstream dump("epitable.dump");
    for (int i=0; i<NDP; i++) {
      phi = dP*i;
      for (int j=0; j<NDR; j++)
	dump << setw(18) << phi 
	     << setw(18) << RDMIN*exp(dR*j)
	     << setw(18) << epitable[i][j]
	     << endl;
      dump << endl;
    }
    dump.close();

    if (JEANS_VR) {
      dump.open("dv2table.dump");
      for (int i=0; i<NDP; i++) {
	phi = dP*i;
	for (int j=0; j<NDR; j++)
	  dump << setw(18) << phi 
	       << setw(18) << RDMIN*exp(dR*j)
	       << setw(18) << dv2table[i][j]
	       << endl;
	dump << endl;
      }
      dump.close();
    }
  }
    
  MPI_Allreduce(&epirmin, &epiRmin, 1, MPI_DOUBLE, MPI_MAX, MPI_COMM_WORLD);
  if (myid==0) cout << "Epitable: Rmin=" << epiRmin << endl;

  if (myid==0)
  {
    ofstream out("ep_disk.dat");
    out.setf(ios::scientific);
    out.precision(8);

    for (int i=0; i<=NDP; i++) {
      phi = dP*i;
      out << "# i=" << " phi=" << phi << ", " << phi*180.0/M_PI << endl;
      for (int j=0; j<NDR; j++) {
	R = RDMIN*exp(dR*j);
	x = R*cos(phi);
	y = R*sin(phi);
	out << setw(18) << x
	    << setw(18) << y
	    << setw(18) << epitable[i%NDP][j] << endl;
      }
      out << endl;
    }
    out << flush;
    out.close();
  }

  if (myid==0)
    {
    ofstream out("table_disk.dat");
    out.setf(ios::scientific);
    out.precision(8);
    
    for (int j=0; j<NDR; j++) {
      for (int k=0; k<NDZ; k++) {
	out << setw(18) << RDMIN*exp(dR*j)
	    << setw(18) << dZ*k
	    << setw(18) << disktableP[0][j][k]
	    << setw(18) << disktableN[0][j][k] 
	    << setw(18) << disk_density(RDMIN*exp(dR*j), dZ*k)
	    << endl;
      }
      out << endl;
    }
  }

  if (myid==0) cout << "[table] " << flush;
}


/*
  Given assertion that $\sigma_r$ follows from Toomre Q condition, solve
  Jeans' equation to get ${\bar v_\phi}$ assume that 
  $\sigma_\phi^2 = \bar{v_\phi^2}$
*/
double DiskHalo::vphi(double xp, double yp, double zp)
{
  double R, ans, R1, R2, X1, Y1, X2, Y2, t2, t1, phi, vc;

  R = sqrt(xp*xp + yp*yp);
  phi = atan2(yp, xp);
  
  const double DR = 0.01;
  R1 = R*(1.0 - DR);
  R2 = R*(1.0 + DR);

  X1 = R1*cos(phi);
  Y1 = R1*sin(phi);

  X2 = R2*cos(phi);
  Y2 = R2*sin(phi);

  t1 = disk_density(R1, 0.0)*vr_disp(X1, Y1, 0.0);
  t2 = disk_density(R2, 0.0)*vr_disp(X2, Y2, 0.0);

  vc = v_circ(xp, yp, 0.0);

  ans = vc*vc + 
    vr_disp(xp, yp, 0.0) * (1.0 + (log(t1)-log(t2))/(log(R1)-log(R2)));
  
  ans = sqrt(max<double>(0.0, ans));

				// For debugging
  if (isnan(ans)) {
    cout << "Process " << myid << ": in vvP has bad answer, x=" << xp
	 << "  y="  << yp << "  z=" << zp
	 << "  R="  << R
	 << "  v0=" << vr_disp(xp, yp, 0.0)
	 << "  t1=" << t1
	 << "  t2=" << t2
	 << "  d1=" << disk_density(R1, 0.0)
	 << "  v1=" << vr_disp(X1, Y1, 0.0)
	 << "  d2=" << disk_density(R2, 0.0)
	 << "  v2=" << vr_disp(X2, Y2, 0.0)
	 << endl;
  }

  return ans;
}


/*
  Interpolate from Jeans' solution integral table
*/
double DiskHalo::get_dispdz(double xp,double yp, double zp)
{
  if (disktableP == NULL) {
    cerr << "DiskHalo::get_dispdz: must call table_disk first\n";
    MPI_Abort(MPI_COMM_WORLD, 100);
    exit(0);
  }

  double R, lR, phi, cp[2], cr[2], cz[2], resvd, zz;
  int iphi, ir, iz;

				// Azimuth
  phi = atan2(yp, xp) + M_PI;
  iphi = (int)( phi/dP );
  iphi = min<int>( iphi, NDP-2 );
  cp[1] = (phi - dP*iphi)/dP;
  cp[0] = 1.0 - cp[1];
				// Cylindrical radius
  R = sqrt(xp*xp + yp*yp);
  lR = log(R);
  ir = (int)( (lR - log(RDMIN))/dR );
  ir = min<int>( ir, NDR-2 );
  ir = max<int>( ir, 0 );
  cr[1] = (lR - log(RDMIN) - dR*ir)/dR;
  cr[0] = 1.0 - cr[1];

				// Zed
  zz = fabs(zp);
  iz = (int)( zz/dZ );
  iz = min<int>( iz, NDZ-2 );
  cz[1] = (zz - dZ*iz)/dZ;
  cz[0] = 1.0 - cz[1];

  if (zp > 0.0) {

    resvd = 

      cp[0]*cr[0]*cz[0] * disktableP[iphi  ][ir  ][iz  ]  +
      cp[0]*cr[0]*cz[1] * disktableP[iphi  ][ir  ][iz+1]  +
      cp[0]*cr[1]*cz[0] * disktableP[iphi  ][ir+1][iz  ]  +
      cp[0]*cr[1]*cz[1] * disktableP[iphi  ][ir+1][iz+1]  +

      cp[1]*cr[0]*cz[0] * disktableP[iphi+1][ir  ][iz  ]  +
      cp[1]*cr[0]*cz[1] * disktableP[iphi+1][ir  ][iz+1]  +
      cp[1]*cr[1]*cz[0] * disktableP[iphi+1][ir+1][iz  ]  +
      cp[1]*cr[1]*cz[1] * disktableP[iphi+1][ir+1][iz+1]  ;

  } else {

    resvd = 

      cp[0]*cr[0]*cz[0] * disktableN[iphi  ][ir  ][iz  ]  +
      cp[0]*cr[0]*cz[1] * disktableN[iphi  ][ir  ][iz+1]  +
      cp[0]*cr[1]*cz[0] * disktableN[iphi  ][ir+1][iz  ]  +
      cp[0]*cr[1]*cz[1] * disktableN[iphi  ][ir+1][iz+1]  +

      cp[1]*cr[0]*cz[0] * disktableN[iphi+1][ir  ][iz  ]  +
      cp[1]*cr[0]*cz[1] * disktableN[iphi+1][ir  ][iz+1]  +
      cp[1]*cr[1]*cz[0] * disktableN[iphi+1][ir+1][iz  ]  +
      cp[1]*cr[1]*cz[1] * disktableN[iphi+1][ir+1][iz+1]  ;

  }

  double dens = disk_density(R, zp);
  if (dens>0.0) resvd /= dens;

  return resvd;
}

/*
  Constant Toomre Q
*/
double DiskHalo::vr_disp(double xp, double yp,double zp)
{
  if (JEANS_VR) {
				// Azimuth
    double phi = atan2(yp, xp) + M_PI;
    int iphi = min<int>( (int)( phi/dP ), NDP-2 );
    double cp[2], cr[2];

    cp[1] = (phi - dP*iphi)/dP;
    cp[0] = 1.0 - cp[1];
				// Cylindrical radius
    double R = sqrt(xp*xp + yp*yp);
    double lR = log(R);
    int ir = (int)( (lR - log(RDMIN))/dR );
    ir = min<int>( ir, NDR-2 );
    ir = max<int>( ir, 0 );

    cr[1] = (lR - log(RDMIN) - dR*ir)/dR;
    cr[0] = 1.0 - cr[1];

    return
      cp[0]*cr[0] * dv2table[iphi  ][ir  ] +
      cp[0]*cr[1] * dv2table[iphi  ][ir+1] +
      cp[1]*cr[0] * dv2table[iphi+1][ir  ] +
      cp[1]*cr[1] * dv2table[iphi+1][ir+1] ;

  } else {
    double r = sqrt(xp*xp + yp*yp);
    double sigmar = 3.36*dmass*disk->get_density(r)*Q/sqrt(epi(xp,yp,zp));
    return sigmar*sigmar;
  }
}

/*
  Analytic rotation curve
*/
double DiskHalo::v_circ(double xp, double yp, double zp)
{
  double R = sqrt(xp*xp + yp*yp);
  double vcirc2 = R*deri_pot(xp, yp, 0.0, 1);

				// Sanity check
  if (vcirc2<=0.0) {
    cout << "DiskHalo::v_circ: circular velocity out of bounds, R="
	 << R << "  v_circ2=" << vcirc2 << endl;
    vcirc2 = 1.0e-20;
  }

  return sqrt(vcirc2);
}


void DiskHalo::
set_vel_disk(vector<Particle>& part)
{
  if (!expandd) {
    if (myid==0) cout << "[no disk particles] ";
    return;
  }

  double vvZ, vvR, vvP;
  double maxVR=-1.0e20, RVR=1e20;
  double maxVP=-1.0e20, RVP=1e20;
  double maxVZ=-1.0e20, RVZ=1e20;
  double vz, vr, vp, R, x, y, z;
  double vel[3], vel1[3], massp, massp1;

  for (int k=0; k<3; k++) vel[k] = vel1[k] = 0.0;
  massp = massp1 = 0.0;

  Normal rn(0.0, 1.0, gen);
				 // Better to make a 2-d table
  table_disk(part);
  
  vector<Particle>::iterator p;
  for (p=part.begin(); p!=part.end(); p++) {

				// From solution to Jeans' equations in
				// cylindrical coordinates
    x = p->pos[0];
    y = p->pos[1];
    z = p->pos[2];
    R = sqrt(x*x + y*y) + MINDOUBLE;

    vvZ = get_dispdz(x, y, z);
    vvR = vr_disp(x, y, z);
    vvP = vvR;
				 // For safety; should only be a problem
				 // on extrapolating the range
    vvZ = max<double>(vvZ, MINDOUBLE);
    vvR = max<double>(vvR, MINDOUBLE);
    vvP = max<double>(vvP, MINDOUBLE);
    
    if (maxVZ < vvZ) {
      maxVZ = vvZ;
      RVZ   = R;
    }
    if (maxVR < vvR) {
      maxVR = vvR;
      RVR   = R;
    }
    if (maxVP < vvP) {
      maxVP = vvP;
      RVP   = R;
    }

    vz   = rn()*sqrt(vvZ);
    vr   = rn()*sqrt(vvR);
    vp   = rn()*sqrt(vvP) + vphi(x, y, z);
    
    p->vel[0] = vr*x/R - vp*y/R;
    p->vel[1] = vr*y/R + vp*x/R;
    p->vel[2] = vz;

    massp1 += p->mass;
    for (int k=0; k<3; k++) vel1[k] += p->mass*p->vel[k];
  }

  MPI_Allreduce(&massp1, &massp, 1, MPI_DOUBLE, MPI_SUM, MPI_COMM_WORLD);
  MPI_Allreduce(vel1,    vel,    3, MPI_DOUBLE, MPI_SUM, MPI_COMM_WORLD);

  if (massp>0.0) {
    for (int k=0; k<3; k++) vel[k] /= massp;
  }
  
  if (myid==0) {
    double v1, v2;
    MPI_Status stt;
    for (int n=1; n<numprocs; n++) {
      MPI_Recv(&v1, 1, MPI_DOUBLE, MPI_ANY_SOURCE, 224, MPI_COMM_WORLD, &stt);
      MPI_Recv(&v2, 1, MPI_DOUBLE, stt.MPI_SOURCE, 225, MPI_COMM_WORLD, &stt);
      if (v1>maxVZ) {
	maxVZ = v1;
	RVZ   = v2;
      }
      MPI_Recv(&v1, 1, MPI_DOUBLE, MPI_ANY_SOURCE, 226, MPI_COMM_WORLD, &stt);
      MPI_Recv(&v2, 1, MPI_DOUBLE, stt.MPI_SOURCE, 227, MPI_COMM_WORLD, &stt);
      if (v1>maxVR) {
	maxVR = v1;
	RVR   = v2;
      }
      MPI_Recv(&v1, 1, MPI_DOUBLE, MPI_ANY_SOURCE, 228, MPI_COMM_WORLD, &stt);
      MPI_Recv(&v2, 1, MPI_DOUBLE, stt.MPI_SOURCE, 229, MPI_COMM_WORLD, &stt);
      if (v1>maxVP) {
	maxVP = v1;
	RVP   = v2;
      }
    }
    cout << "     *****";
    cout << " (u, v, w)=(" << vel[0] 
	 << ", " << vel[1]
	 << ", " << vel[2] << ")" << endl;
    cout <<   "maxVZ=" << maxVZ << " (" << RVZ << ")"
	 << ", maxVR=" << maxVR << " (" << RVR << ")"
	 << ", maxVP=" << maxVP << " (" << RVP << ")"
	 << endl;
  } else {
    MPI_Send(&maxVZ, 1, MPI_DOUBLE, 0, 224, MPI_COMM_WORLD);
    MPI_Send(&RVZ,   1, MPI_DOUBLE, 0, 225, MPI_COMM_WORLD);
    MPI_Send(&maxVR, 1, MPI_DOUBLE, 0, 226, MPI_COMM_WORLD);
    MPI_Send(&RVR,   1, MPI_DOUBLE, 0, 227, MPI_COMM_WORLD);
    MPI_Send(&maxVP, 1, MPI_DOUBLE, 0, 228, MPI_COMM_WORLD);
    MPI_Send(&RVP,   1, MPI_DOUBLE, 0, 229, MPI_COMM_WORLD);
  }

  if (cov) {
    vector<Particle>::iterator ip;
    for(ip=part.begin(); ip!=part.end(); ip++) {
      for (int k=0; k<3; k++) ip->vel[k] -= vel[k];
    }
  }

}

void DiskHalo::table_halo_disp()
{
  if (!MULTI) {
    string msg("DiskHalo::table_halo_disp is only valid if MULTI is true");
    throw DiskHaloException(msg, __FILE__, __LINE__);
  }

  if (halotable.getchigh() == NHR-1 && halotable.getrhigh() == 2) return;
  
  halotable.setsize(0, 2, 0, NHR-1);
  
  double rmin = halo2->get_min_radius();
  double rmax = halo2->get_max_radius();
  double Emax = halo2->get_pot(rmax);
  double dr, r, pot, E, v, v2, fac;

  if (rmin<=0.0 && LOGR) {
    string msg("DiskHalo::table_halo_disp: LOGR=1 but rmin<=0");
    throw DiskHaloException(msg, __FILE__, __LINE__);
  }

  if (LOGR)
    dr = (log(rmax) - log(rmin))/(NHR - 1);
  else
    dr = (rmax - rmin)/(NHR - 1);
  
  halotable.zero();
  
  const int nlq = 400;
  LegeQuad lq(nlq);

  int ibeg = (int)( myid*(double)NHR/numprocs );
  int iend = (int)( (myid+1)*(double)NHR/numprocs );

  if (myid==numprocs-1) iend = NHR;

  for (int i=ibeg; i<iend; i++) {
    if (LOGR)
      halotable[0][i] = r = exp(log(rmin) + dr*i);
    else
      halotable[0][i] = r = rmin + dr*i;

    pot = halo2->get_pot(r);

    for (int n=1; n<=nlq; n++) {
      E = pot + (Emax - pot)*lq.knot(n);
      v2 = 2.0*(E - pot);
      if (v2<0.0) v2 = 0.0;
      v = sqrt(v2);
      fac = (Emax - pot)*lq.weight(n)*halo2->distf(E, 0.5);
				// velocity dispersion
      halotable[1][i] += fac * v * v2;
				// density
      halotable[2][i] += fac * v;
    }

    if (halotable[2][i]>0.0) halotable[1][i] /= halotable[2][i];
    else halotable[1][i] = 0.0;
  }


  // Update tables on all nodes
  //
  MPI_Allreduce(&halotable[0][0], &halotable[0][0], NHR, MPI_DOUBLE, 
		MPI_SUM, MPI_COMM_WORLD);
  MPI_Allreduce(&halotable[1][0], &halotable[1][0], NHR, MPI_DOUBLE, 
		MPI_SUM, MPI_COMM_WORLD);
  MPI_Allreduce(&halotable[2][0], &halotable[2][0], NHR, MPI_DOUBLE, 
		MPI_SUM, MPI_COMM_WORLD);
  
  
  //
  // DEBUG output
  //
  if (myid==0 && VFLAG & 4) {
    ofstream out("disp_halo.dat");
    out.setf(ios::scientific);
    out.precision(2);
    
    for (int k=0; k<NHR; k++)
      out << setw(18) << halotable[0][k]
	  << setw(18) << halotable[1][k]
	  << setw(18) << 4.0*M_PI*halotable[2][k]
	  << setw(18) << halo->get_density(halotable[0][k]) << endl;
  }

}

void DiskHalo::
table_halo(vector<Particle>& part)
{
  if (halotable.getchigh() == NHR-1) return;
  
  halotable.setsize(0, NHT-1, 0, NHR-1);
  
  dc = 2.0/(NHT-1);
  double r2, maxr = 0.0, maxr1 = 0.0;
  vector<Particle>::iterator p;
  for (p=part.begin(); p!=part.end(); p++) {
    r2 = 0.0;
    for (int k=0; k<3; k++) r2 += p->pos[k]*p->pos[k];
    maxr1 = max<double>(maxr1, sqrt(r2));
  }
  
  MPI_Allreduce(&maxr1, &maxr, 1, MPI_DOUBLE, MPI_MAX, MPI_COMM_WORLD);

  dr = (log(max<double>(RHMAX, maxr)) - log(RHMIN))/(NHR-1);
  
  if (myid==0) {
    cout << "Table halo: RHMIN=" << RHMIN 
	 << " RHMAX=" << RHMAX
	 << " maxr=" << maxr
	 << " dr=" << dr
	 << endl;
  }

  double x, y, z, theta, r, fr, fz, fp, pot, costh, sinth;
  double dens, potl, dpr, dpt, dpp, dpdr;
  Vector work(0, NHR-1);
  Vector workR(0, NHR-1);
  Vector workA(0, NHR-1);
  
  // If no disk, add no force
  //
  pot = fr = fz = fp = 0.0;
  
  // Compute this table in parallel
  //
  vector<int> ibeg(numprocs), iend(numprocs);
  for (int i=0; i<numprocs; i++) {
    ibeg[i] = (i  )*NHT/numprocs;
    iend[i] = (i+1)*NHT/numprocs;
  }
  
  for (int i=ibeg[myid]; i<iend[myid]; i++) {
    
    costh = -1.0 + dc*i;
    sinth = sqrt(1.0 - costh*costh);
    
    for (int j=0; j<NHR; j++) {
      
      workR[j] = log(RHMIN) + dr*j;
      r = exp(workR[j]);
      x = r*sinth;
      y = 0.0;
      z = r*costh;
      
      if (expandd) disk_eval(x, z, 0.0, pot, fr, fz, fp);
      
      theta = acos(z/(r + MINDOUBLE));
      expandh->determine_fields_at_point(r, theta, 0.0,
					 &dens, &potl, &dpr, &dpt, &dpp);
      
      dpdr = dpr - fz*cos(theta) - fr*sin(theta);
      
      work[j] = halo->get_density(r) * dpdr * r;
      
    }
    
    // Splsum(workR, work, workA);
    Trapsum(workR, work, workA);
    for (int j=0; j<NHR; j++) {
      halotable[i][j] = max<double>(workA[NHR-1] - workA[j], MINDOUBLE);
      if (fabs(halotable[i][j])>1.0e8) {
	cerr << "Oops, val=" << halotable[i][j] << endl;
      }
    }
  }
  
  // Update tables on all nodes
  //
  for (int k=0; k<numprocs; k++) {
    for (int i=ibeg[k]; i<iend[k]; i++) {
      MPI_Bcast(&halotable[i][0], NHR, MPI_DOUBLE, k, MPI_COMM_WORLD);
    }
  }
  
  //
  // DEBUG output
  //
  if (myid==0 && expandh && VFLAG & 4) {
    ofstream out("table_halo.dat");
    out.setf(ios::scientific);
    out.precision(8);
    
    for (int j=0; j<NHT; j++) {
      for (int k=0; k<NHR; k++) {
	out << setw(18) << -1.0 + dc*j
	    << setw(18) << RHMIN*exp(dr*k)
	    << setw(18) << halotable[j][k] << endl;
      }
      out << endl;
    }
    
    ofstream out2("test_halo_pot.dat");
    out2.setf(ios::scientific);
    out2.precision(2);
    
    double rr, t, p, dp;
    for (int k=0; k<NHR; k++) {
      rr = RHMIN*exp(dr*k);
      expandh->determine_fields_at_point(rr, 0.5*M_PI, 0.0,
					 &t, &p, &dp, &t, &t);
      out2
	<< setw(14) << rr
	<< setw(14) << p
	<< setw(14) << dp
	<< setw(14) << get_disp(rr, 0.0, 0.0)
	<< endl;
    }
    
    
    ofstream out3("disp_diff.dat");
    double xx, zz, costh;
    
    for (int j=0; j<NHT; j++) {
      costh = -1.0 + dc*j;
      out3 << setw(14) << costh;
      for (int k=0; k<NHR; k++) {
	rr = RHMIN*exp(dr*k);
	xx = rr*sqrt(1.0 - costh*costh);
	zz = rr*costh;
	out3 
	  << setw(14) << halotable[j][k] -  
	  get_disp(xx, 0.0, zz) * halo->get_density(rr);
      }
      out3 << endl;
    }
  }
  
  if (myid==0) cout << "[table] " << flush;
}

double DiskHalo::get_disp(double xp,double yp, double zp)
{
  if (MULTI) {
    double r = sqrt(xp*xp + yp*yp + zp*zp);
    r = max<double>(r, halo2->get_min_radius());
    r = min<double>(r, halo2->get_max_radius());
    return odd2(r, halotable[0], halotable[1], 0);
  }

  if (halotable.getchigh() != NHR-1) {
    cerr << "DiskHalo::get_disp: must call table_halo first\n";
    MPI_Abort(MPI_COMM_WORLD, 100);
    exit(0);
  }
  
  int it, ir;
  double r, lr, t, ct[2], cr[2], resv;
  
  // Polar angle
  //
  r = sqrt(xp*xp + yp*yp + zp*zp);
  t = zp/(r + MINDOUBLE) + 1.0;
  
  it = (int)( t/dc );
  it = min<int>( it, NHT-2 );
  ct[1] = (t - dc*it)/dc;
  ct[0] = 1.0 - ct[1];
  // Polar radius
  //
  lr = log(r);
  ir = (int)( (lr - log(RHMIN))/dr );
  ir = min<int>( ir, NHR-2 );
  ir = max<int>( ir, 0);
  cr[1] = (lr - log(RHMIN) - dr*ir)/dr;
  cr[0] = 1.0 - cr[1];
  
  resv = 
    
    ct[0]*cr[0] * halotable[it  ][ir  ]  +
    ct[0]*cr[1] * halotable[it  ][ir+1]  +
    ct[1]*cr[0] * halotable[it+1][ir  ]  +
    ct[1]*cr[1] * halotable[it+1][ir+1]  ;
  
  double dens = halo->get_density(r);
  if (dens>0.0) resv /= dens;
  
  return resv;
}

void DiskHalo::set_vel_halo(vector<Particle>& part)
{
  if (!expandh) {
    if (myid==0) cout << "[no halo particles] ";
    return;
  }
  
  int nok;
  double v2r, vr, r;
  double vel[3], vel1[3], massp, massp1;
  
  for (int k=0; k<3; k++) vel[k] = vel1[k] = 0.0;
  massp = massp1 = 0.0;
  
  table_halo(part);
  
  vector<Particle>::iterator p;
  for (p=part.begin(); p!=part.end(); p++) {
    
    r = sqrt(p->pos[0]*p->pos[0] + 
	     p->pos[1]*p->pos[1] +
	     p->pos[2]*p->pos[2]);
    
				// Reset success flag
    nok = 1;
    
				// Use Eddington
    
    if (DF && 0.5*(1.0+erf((r-R_DF)/DR_DF)) > (*rndU)()) {
      halo2->gen_velocity(&p->pos[0], &p->vel[0], nok);
      
      if (nok) {
	cout << "gen_velocity failed: "
	     << p->pos[0] << " "
	     << p->pos[1] << " "
	     << p->pos[2] << "\n";
      }
    }
				// Use Jeans
    if (nok) {
      v2r = get_disp(p->pos[0], p->pos[1], p->pos[2]);
      vr = sqrt(max<double>(v2r, MINDOUBLE));
      for (int k=0; k<3; k++) p->vel[k] = vr*(*rndN)();
    }
    
    massp1 += p->mass;
    for (int k=0; k<3; k++) vel1[k] += p->mass*p->vel[k];
  }
  

  MPI_Allreduce(&massp1, &massp, 1, MPI_DOUBLE, MPI_SUM, MPI_COMM_WORLD);
  MPI_Allreduce(vel1, vel, 3, MPI_DOUBLE, MPI_SUM, MPI_COMM_WORLD);
    
  if (massp>0.0) {
    for (int k=0; k<3; k++) vel[k] /= massp;
  }
    
  if (myid==0) {
    cout << " (u, v, w)=(" << vel[0] 
	 << ", " << vel[1]
	 << ", " << vel[2] << ")" << endl;
  }

  if (cov) {
    vector<Particle>::iterator ip;
    for(ip=part.begin(); ip!=part.end(); ip++) {
      for (int k=0; k<3; k++) ip->vel[k] -= vel[k];
    }
  }
  
}

void DiskHalo::write_record(ostream &out, Particle &p)
{
  out << " " << setw(16) << setprecision(8) << p.mass;
  
  for (int k=0; k<3; k++)
    out << setw(24) << setprecision(15) << p.pos[k] + center_pos[k];
  
  for (int k=0; k<3; k++)
    out << setw(24) << setprecision(15) << p.vel[k] + center_vel[k];
  
  out << endl;
}


void DiskHalo::write_file(ostream &fou_halo,  ostream &fou_disk,
			  vector<Particle>& hpart, vector<Particle>& dpart)
{
  int l  = hpart.size();
  int l1 = dpart.size();
  
  vector<Particle> buf(NBUF);
  
  // Make MPI datatype
  //
  MPI_Datatype Particletype;
  MPI_Datatype type[4] = {MPI_UNSIGNED, MPI_DOUBLE, MPI_DOUBLE, MPI_DOUBLE};
  
  // Get displacements
  //
  MPI_Aint disp[4];
  MPI_Get_address(&buf[0].level,	&disp[0]);
  MPI_Get_address(&buf[0].mass,		&disp[1]);
  MPI_Get_address(&buf[0].pos[0],	&disp[2]);
  MPI_Get_address(&buf[0].vel[0],	&disp[3]);
  
  for (int i=3; i>=0; i--) disp[i] -= disp[0];
  
				// Block offsets
  int blocklen[4] = {1, 1, 3, 3};
  
  // Make and register the new type
  //
  MPI_Type_create_struct(4, blocklen, disp, type, &Particletype);
  MPI_Type_commit(&Particletype);
  
  // Get particle totals
  //
  int ndisk=0, nhalo=0;
  MPI_Reduce(&l,  &nhalo, 1, MPI_INT, MPI_SUM, 0, MPI_COMM_WORLD);
  MPI_Reduce(&l1, &ndisk, 1, MPI_INT, MPI_SUM, 0, MPI_COMM_WORLD);
  
  if (myid==0) {
    
    if (VFLAG & 1)
      cout << endl
	   << "Total number of particles: nhalo=" << nhalo
	   << " ndisk=" << ndisk << endl;
    
    fou_halo.setf(ios::scientific);
    fou_disk.setf(ios::scientific);
    
    fou_halo << nhalo << " " << 0 << " " << 0 << endl;
    fou_disk << ndisk << " " << 0 << " " << 0 << endl;
    
    if (VFLAG & 1) {
      cout << "Halo stream is ";
      if (fou_halo.good()) cout << "GOOD\n";
      else cout << "BAD\n";

      cout << "Disk stream is ";
      if (fou_disk.good()) cout << "GOOD\n";
      else cout << "BAD\n";
    }

    for (int i=0; i<l; i++)
      write_record(fou_halo, hpart[i]);
    
    for (int i=0; i<l1; i++)
      write_record(fou_disk, dpart[i]);
    
    if (VFLAG & 1) {
      cout << "Wrote " << l  << " HALO particles from Node 0" << endl;
      cout << "Wrote " << l1 << " DISK particles fron Node 0" << endl;
    }

    int imany, icur, ccnt;

    for (int n=1; n<numprocs; n++) {
      
      MPI_Recv(&imany, 1, MPI_INT, n, 10, MPI_COMM_WORLD, MPI_STATUS_IGNORE);
      ccnt=0;
      while (ccnt<imany) {
	MPI_Recv(&icur, 1, MPI_INT, n, 11, MPI_COMM_WORLD, MPI_STATUS_IGNORE);
	MPI_Recv(&buf[0], icur, Particletype, n, 12, MPI_COMM_WORLD, MPI_STATUS_IGNORE);
	for (int i=0; i<icur; i++) write_record(fou_halo, buf[i]);
	ccnt += icur;
      }
      
      if (VFLAG & 1)
	cout << "Wrote " << ccnt << " HALO particles from Node " << n << endl;

      MPI_Recv(&imany, 1, MPI_INT, n, 13, MPI_COMM_WORLD, MPI_STATUS_IGNORE);
      ccnt = 0;
      while (ccnt<imany) {
	MPI_Recv(&icur, 1, MPI_INT, n, 14, MPI_COMM_WORLD, MPI_STATUS_IGNORE);
	MPI_Recv(&buf[0], icur, Particletype, n, 15, MPI_COMM_WORLD, MPI_STATUS_IGNORE);
	for (int i=0; i<icur; i++) write_record(fou_disk, buf[i]);
	ccnt += icur;
      }
      
      if (VFLAG & 1)
	cout << "Wrote " << ccnt << " DISK particles from Node " << n << endl;

      MPI_Barrier(MPI_COMM_WORLD);
    }
    
  } else {
    
    int icur, ipack;

    for (int n=1; n<numprocs; n++) {

      if (myid==n) {

	MPI_Send(&l, 1, MPI_INT, 0, 10, MPI_COMM_WORLD);
	icur = 0;
	while (icur<l) {
	  ipack = min<int>(l-icur, NBUF);
	  for (int j=0; j<ipack; j++) buf[j] = hpart[icur+j];
	  MPI_Send(&ipack, 1, MPI_INT, 0, 11, MPI_COMM_WORLD);
	  MPI_Send(&buf[0], ipack, Particletype, 0, 12, MPI_COMM_WORLD);
	  icur += ipack;
	}

	if (VFLAG & 1)
	  cout << "Sent " << icur << " HALO particles from Node " << n << endl;

	MPI_Send(&l1, 1, MPI_INT, 0, 13, MPI_COMM_WORLD);
	icur = 0;
	while (icur<l1) {
	  ipack = min<int>(l1-icur, NBUF);
	  for (int j=0; j<ipack; j++) buf[j] = dpart[icur+j];
	  MPI_Send(&ipack, 1, MPI_INT, 0, 14, MPI_COMM_WORLD);
	  MPI_Send(&buf[0], ipack, Particletype, 0, 15, MPI_COMM_WORLD);
	  icur += ipack;
	}

	if (VFLAG & 1)
	  cout << "Sent " << icur << " DISK particles from Node " << n << endl;

      }
      MPI_Barrier(MPI_COMM_WORLD);
    }
  }
  
}


void DiskHalo::virial_ratio(vector<Particle>& hpart, vector<Particle>& dpart)
{
  double r, theta, phi, xx, yy, zz, axd, ayd, azd, axh, ayh, azh, R2, R;
  double dens, potl, potr, pott, potp, fr, fp, fz;
  
  double KE_disk1 = 0.0;
  double KE_halo1 = 0.0;
  double PE_disk_disk1 = 0.0;
  double PE_halo_disk1 = 0.0;
  double PE_disk_halo1 = 0.0;
  double PE_halo_halo1 = 0.0;
  double massd1 = 0.0;
  double massh1 = 0.0;
  
  vector<Particle>::iterator p;
  
  fr = fp = fz = 0.0;
  potr = pott = potp = 0.0;
				// -----------------
				// Halo contribution
				// -----------------
  for (p=hpart.begin(); p!=hpart.end(); p++) {
    r = 0.0;
    for (int k=0; k<3; k++) {
      r += p->pos[k]*p->pos[k];
      KE_halo1 += 0.5*p->mass*p->vel[k]*p->vel[k];
    }
    
    r = sqrt(r);
    xx = p->pos[0];
    yy = p->pos[1];
    zz = p->pos[2];
    
    theta = acos(zz/(r+MINDOUBLE));
    phi = atan2(yy, xx);
    
    if (expandh)
      expandh->determine_fields_at_point(r, theta, phi,
					 &dens, &potl, &potr, &pott, &potp);
    
    R2 = xx*xx + yy*yy + MINDOUBLE;
    R = sqrt(R2);
    
    if (expandd) 
      disk_eval(R, zz, phi, potl, fr, fz, fp);
    
    axd = fr*xx/R - fp*yy/R2;
    ayd = fr*yy/R + fp*xx/R2;
    azd = fz;
    
    axh = -(potr*xx/r - pott*xx*zz/(r*r*r)) + potp*yy/R2;
    ayh = -(potr*yy/r - pott*yy*zz/(r*r*r)) - potp*xx/R2;
    azh = -(potr*zz/r + pott*R2/(r*r*r));
    
				// Clausius
    PE_halo_disk1 += p->mass * (xx*axd + yy*ayd + zz*azd);
    PE_halo_halo1 += p->mass * (xx*axh + yy*ayh + zz*azh);
    
				// Mass
    massh1 += p->mass;
  }
  
				// -----------------
				// Disk contribution
				// -----------------
  for (p=dpart.begin(); p!=dpart.end(); p++) {
    r = 0.0;
    for (int k=0; k<3; k++) {
      r += p->pos[k]*p->pos[k];
      KE_disk1 += 0.5*p->mass*p->vel[k]*p->vel[k];
    }
    
    r = sqrt(r);
    xx = p->pos[0];
    yy = p->pos[1];
    zz = p->pos[2];
    
    theta = acos(zz/(r+MINDOUBLE));
    phi = atan2(yy, xx);
    
    if (expandh)
      expandh->determine_fields_at_point(r, theta, phi,
					 &dens, &potl, &potr, &pott, &potp);
    
    R2 = xx*xx + yy*yy;
    R = sqrt(R2);
    
    if (expandd) 
      disk_eval(R, zz, phi, potl, fr, fz, fp);
    
    axd = fr*xx/R - fp*yy/R2;
    ayd = fr*yy/R + fp*xx/R2;
    azd = fz;
    
    axh = -(potr*xx/r - pott*xx*zz/(r*r*r)) + potp*yy/R2;
    ayh = -(potr*yy/r - pott*yy*zz/(r*r*r)) - potp*xx/R2;
    azh = -(potr*zz/r + pott*R2/(r*r*r));
    
				// Clausius
    PE_disk_disk1 += p->mass * (xx*axd + yy*ayd + zz*azd);
    PE_disk_halo1 += p->mass * (xx*axh + yy*ayh + zz*azh);
    
				// Mass
    massd1 += p->mass;
    
  }
  
  
  double KE_disk = 0.0;
  double KE_halo = 0.0;
  
  double PE_disk_disk = 0.0;
  double PE_disk_halo = 0.0;
  
  double PE_halo_halo = 0.0;
  double PE_halo_disk = 0.0;
  
  double massh = 0.0;
  double massd = 0.0;
  
  MPI_Reduce(&KE_disk1, &KE_disk, 1, MPI_DOUBLE, MPI_SUM, 0, MPI_COMM_WORLD);
  MPI_Reduce(&KE_halo1, &KE_halo, 1, MPI_DOUBLE, MPI_SUM, 0, MPI_COMM_WORLD);
  
  MPI_Reduce(&PE_disk_disk1, &PE_disk_disk, 1, MPI_DOUBLE, MPI_SUM, 0, MPI_COMM_WORLD);
  MPI_Reduce(&PE_disk_halo1, &PE_disk_halo, 1, MPI_DOUBLE, MPI_SUM, 0, MPI_COMM_WORLD);
  
  MPI_Reduce(&PE_halo_halo1, &PE_halo_halo, 1, MPI_DOUBLE, MPI_SUM, 0, MPI_COMM_WORLD);
  MPI_Reduce(&PE_halo_disk1, &PE_halo_disk, 1, MPI_DOUBLE, MPI_SUM, 0, MPI_COMM_WORLD);
  
  MPI_Reduce(&massh1, &massh, 1, MPI_DOUBLE, MPI_SUM, 0, MPI_COMM_WORLD);
  MPI_Reduce(&massd1, &massd, 1, MPI_DOUBLE, MPI_SUM, 0, MPI_COMM_WORLD);
  
  if (myid==0) {
    
    double PE_halo = PE_halo_halo + PE_halo_disk;
    double PE_disk = PE_disk_disk + PE_disk_halo;
    
    cout << endl
	 << "****************************" << endl
	 <<"  KE_halo  = " << KE_halo << endl
	 <<"  KE_disk  = " << KE_disk << endl << endl
	 <<"  PE_halo(halo)  = " << PE_halo_halo << endl
	 <<"  PE_halo(disk)  = " << PE_halo_disk << endl
	 <<"  PE_disk(disk)  = " << PE_disk_disk << endl
	 <<"  PE_disk(halo)  = " << PE_disk_halo << endl << endl;
    if (PE_halo < 0.0)
      cout <<"-2T/W_halo = " << -2.0*KE_halo/PE_halo << endl;
    if (PE_disk < 0.0)
      cout <<"-2T/W_disk = " << -2.0*KE_disk/PE_disk << endl;
    cout << endl;
    
    cout << " Halo mass=" << massh << "  Disk mass=" << massd << endl << endl;
    
    double KE = KE_halo + KE_disk;
    double PE = PE_halo + PE_disk;
    
    cout << "  KE       = " << KE << endl
	 << "  PE       = " << PE << endl;
    if (PE<0.0)
      cout << " -2T/W     = " << -2.0*KE/PE << endl;
    cout << "****************************" << endl;
    
  } 
}


void DiskHalo::virial_ratio(const char *hfile, const char *dfile)
{
  ifstream in[2];
  
  in[0].open(hfile);
  in[1].open(dfile);

  if (!in[0]) {
    cout << "virial_ratio: can't open " << hfile << endl;
    return;
  }
  
  if (!in[1]) {
    cout << "virial_ratio: can't open " << dfile << endl;
    return;
  }
  
  double r, theta, phi, xx, yy, zz, ax, ay, az, R2, R;
  double dens, potl, potr, pott, potp, fr, fp, fz;
  
  double KE_disk = 0.0;
  double KE_halo = 0.0;
  double PE_disk = 0.0;
  double PE_halo = 0.0;
  double massd = 0.0;
  double massh = 0.0;
  
  fr = fp = fz = 0.0;
  potr = pott = potp = 0.0;
  
  int expected, count=0;
  double m, pos[3], vel[3];
  
  for (int c=0; c<2; c++) {

    const int linesize = 1024;
    char linebuffer[linesize];
    in[c].getline(linebuffer, linesize);
    istringstream ins(linebuffer);
    ins >> expected;
  
    int ibeg = (myid+0)*expected/numprocs;
    int iend = (myid+1)*expected/numprocs;
    if (myid == numprocs - 1) iend = expected;

    for (int i=0; i<ibeg; i++) 
      in[c].getline(linebuffer, linesize);
    
    for (int i=ibeg; i<iend; i++) {

      in[c].getline(linebuffer, linesize);
      if (!in) break;
    
      istringstream ins(linebuffer);
    
      ins >> m;
      for (int i=0; i<3; i++) ins >> pos[i];
      for (int i=0; i<3; i++) ins >> vel[i];
      if (!ins) break;
      
      count++;

      r = 0.0;
      for (int k=0; k<3; k++) {
	r += pos[k]*pos[k];
	if (c==1)
	  KE_halo += 0.5*m*vel[k]*vel[k];
	if (c==2)
	  KE_disk += 0.5*m*vel[k]*vel[k];
      }
      
      r = sqrt(r);
      xx = pos[0];
      yy = pos[1];
      zz = pos[2];

      theta = acos(zz/(r+MINDOUBLE));
      phi = atan2(yy, xx);
    
      if (expandh)
	expandh->determine_fields_at_point(r, theta, phi,
					 &dens, &potl, &potr, &pott, &potp);
      
      R2 = xx*xx + yy*yy + MINDOUBLE;
      R = sqrt(R2);

      if (expandd)
	disk_eval(R, zz, phi, potl, fr, fz, fp);

      ax = -(potr*xx/r - pott*xx*zz/(r*r*r)) + fr*xx/R - fp*yy/R2;
      ay = -(potr*yy/r - pott*yy*zz/(r*r*r)) + fr*yy/R + fp*xx/R2;
      az = -(potr*zz/r + pott*R2/(r*r*r))    + fz;
      ax +=  potp*yy/R2;
      ay += -potp*xx/R2;
      
				// Clausius and mass
      if (c==0) {
	PE_halo += m * (xx*ax + yy*ay + zz*az);
	massh += m;
      }
    
      if (c==1) {
	PE_disk += m * (xx*ax + yy*ay + zz*az);
	massd += m;
      }
    
    }
  }


  double KE_disk0;
  double KE_halo0;
  double PE_disk0;
  double PE_halo0;
  double massd0;
  double massh0;
  int count0;
  
  MPI_Reduce(&KE_disk, &KE_disk0, 1, MPI_DOUBLE, MPI_SUM, 0, MPI_COMM_WORLD);
  MPI_Reduce(&KE_halo, &KE_halo0, 1, MPI_DOUBLE, MPI_SUM, 0, MPI_COMM_WORLD);
  MPI_Reduce(&PE_disk, &PE_disk0, 1, MPI_DOUBLE, MPI_SUM, 0, MPI_COMM_WORLD);
  MPI_Reduce(&PE_halo, &PE_halo0, 1, MPI_DOUBLE, MPI_SUM, 0, MPI_COMM_WORLD);
  MPI_Reduce(&massd,   &massd0,   1, MPI_DOUBLE, MPI_SUM, 0, MPI_COMM_WORLD);
  MPI_Reduce(&massh,   &massh0,   1, MPI_DOUBLE, MPI_SUM, 0, MPI_COMM_WORLD);
  MPI_Reduce(&count,   &count0,   1, MPI_INT,    MPI_SUM, 0, MPI_COMM_WORLD);

  if (myid==0) {

    cout << endl
	 << "****************************" << endl
	 <<"  Total #  = " << count0 << endl
	 <<"  Expect # = " << expected << endl << endl
	 <<"  KE_halo  = " << KE_halo0 << endl
	 <<"  KE_disk  = " << KE_disk0 << endl << endl
	 <<"  PE_halo  = " << PE_halo0 << endl
	 <<"  PE_disk  = " << PE_disk0 << endl << endl;

    if (PE_halo0 < 0.0)
      cout <<"-2T/W_halo = " << -2.0*KE_halo0/PE_halo0 << endl;
    if (PE_disk < 0.0)
      cout <<"-2T/W_disk = " << -2.0*KE_disk0/PE_disk0 << endl;
    cout << endl;

    cout << " Halo mass=" << massh0 << "  Disk mass=" << massd0 
	 << endl << endl;

    double KE = KE_halo0 + KE_disk0;
    double PE = PE_halo0 + PE_disk0;

    cout << "  KE       = " << KE << endl
	 << "  PE       = " << PE << endl;
    if (PE<0.0)
      cout << " -2T/W     = " << -2.0*KE/PE << endl;
    cout << "****************************" << endl;
    
  }

}








