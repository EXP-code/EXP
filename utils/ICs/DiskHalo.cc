#define LOGCHEBY		// Test smoothing using log scaling
				// from Mike P

#undef ENFORCE_KAPPA		// Clamp kappa^2

				// System
#include <values.h>
				// C++/STL
#include <iostream>
#include <fstream>
#include <sstream>
#include <algorithm>
#include <memory>
#include <vector>
				// EXP classes
#include <interp.H>
#include <numerical.H>
#include <exponential.H>
#include <interp.H>
				// Local
#include <AddDisk.H>
#include <DiskHalo.H>
				// Grid parameters and Toomre Q
double DiskHalo::RHMIN       = 1.0e-4;
double DiskHalo::RHMAX       = 50.0;
double DiskHalo::RDMIN       = 1.0e-4;
double DiskHalo::RDMAX       = 20.0;
double DiskHalo::Q           = 0.0;
double DiskHalo::SIG0        = 0.1;
double DiskHalo::XI          = 1.0;
double DiskHalo::SHFACTOR    = 16.0;
double DiskHalo::TOLE        = 0.003;
double DiskHalo::COMPRESSION = 1.0;
int    DiskHalo::NDP         = 16;
int    DiskHalo::NDZ         = 40;
int    DiskHalo::NDR         = 800;
int    DiskHalo::NHR         = 800;
int    DiskHalo::NHT         = 40;
int    DiskHalo::SEED        = 11;
int    DiskHalo::ITMAX       = 1000000;

double DiskHalo::RA          = 1.0e20;
int    DiskHalo::NUMDF       = 800;
int    DiskHalo::RNUM        = 4000;

double DiskHalo::R_DF        = 20.0;
double DiskHalo::DR_DF       = 5.0;

bool   DiskHalo::LOGR        = true;

// these appear to be good settings, but may not be the best. use with caution!
int    DiskHalo::NCHEB       = 8;
bool   DiskHalo::CHEBY       = false;
bool   DiskHalo::ALLOW       = false;
bool   DiskHalo::use_mono    = true;

unsigned DiskHalo::VFLAG     = 7;
unsigned DiskHalo::NBUF      = 65568;

string DiskHalo::RUNTAG      = "debug";

// lower-case string to enum map for setting disk-velocity type
std::map<std::string, DiskHalo::DiskGenType>
DiskHalo::getDiskGenType = {
  {"jeans",      DiskHalo::Jeans     },
  {"asymmetric", DiskHalo::Asymmetric},
  {"epicyclic",  DiskHalo::Epicyclic }
};

static std::shared_ptr<AxiSymModel> model;
double targetmass;
				// Determine radius with given enclosed mass
double mass_func(double r)
{
  return targetmass - model->get_mass(r);
}

DiskHalo::
DiskHalo()
{
  com        = false;
  cov        = false;
  DF         = false;
  MULTI      = false;
  type       = Jeans;
  bufcnt     = 0;
}

DiskHalo::
DiskHalo(SphericalSLptr haloexp, EmpCylSLptr diskexp,
	 double H, double A, double DMass, std::string maptype,
	 string& filename, int DF1, int DIVERGE, double DIVERGE_RFAC,
	 DiskGenType type)
{
  gen.seed(SEED+myid);

  com        = false;
  cov        = false;
  mtype      = maptype;
  this->type = type;
  cmap       = make_cmap(mtype, A);

  for (int k=0; k<3; k++) center_pos[k] = center_vel[k] = 0.0;

  if (DF1) DF = true;
  else     DF = false;

  MULTI       = false;

  dmass       = DMass;
  scalelength = A;
  scaleheight = H;

  expandh     = haloexp;
  expandd     = diskexp;

  random_gen.seed(SEED*numprocs + myid);
  SphericalModelTable::even = 0;

  halo = std::make_shared<SphericalModelTable>(filename, DIVERGE, DIVERGE_RFAC);

  disk = std::make_shared<ExponentialDisk>(A, RDMAX);

  if (myid==0 && VFLAG & 1) {
    std::cerr << "DiskHalo: DIVERGE=" << DIVERGE
	 << " A=" << A
	 << " RDMAX=" << RDMAX
	 << " filename=" << filename
	 << std::endl;
  }

  if (DF) {
    AddDisk::logarithmic = true;
    AddDisk::Rmin = RHMIN;
    AxiSymModel::numr = 400;
    AxiSymModel::numj = 400;
    AxiSymModel::gen_N = 800;
    AxiSymModel::gen_itmax = 400000;
    AxiSymModel::gen_rmin = RHMIN;
    newmod = std::make_shared<AddDisk>(halo, disk, dmass*COMPRESSION); 
    halo2 = newmod->get_model();
    halo2->print_model("diskhalo.newmodel");
    halo2->setup_df(NUMDF, RA);
    if (myid==0 && VFLAG & 2) {
      char debugname[] = "df.debug";
      halo2->print_df(debugname);
    }
  } else {
    halo2 = halo;
  }

  double rmin = max<double>(disk->get_min_radius(), RDMIN);
  double rmax = disk->get_max_radius();
  double mmin = disk->get_mass(rmin);
  double mtot = disk->get_mass(rmax);

  double R, phi;
  double pos[3], pos1[3], massp, massp1;

  // For frequency computation
  //
  hDmin = log(max<double>(disk->get_min_radius(), RDMIN));
  hDmax = log(disk->get_max_radius());
  dRh = (hDmax - hDmin)/nh;
  nhN = vector<unsigned>(nh+1, 0  ); // Counts per bin
  nhD = vector<double>  (nh+1, 0.0); // Mass per bin
  nhM = vector<double>  (nh+1, 0.0); // Cumulative mass

  // For buffered ascii writes
  //
  bufcnt = 0;
}

DiskHalo::
DiskHalo(SphericalSLptr haloexp, EmpCylSLptr diskexp,
	 double H, double A, double DMass, std::string maptype,
	 std::string& filename1, int DIVERGE, double DIVERGE_RFAC,
	 std::string& filename2, int DIVERGE2, double DIVERGE_RFAC2,
	 DiskGenType type)
{
  gen.seed(SEED+myid);

  com        = false;
  cov        = false;
  mtype      = maptype;
  this->type = type;
  cmap       = make_cmap(mtype, A);

  for (int k=0; k<3; k++) center_pos[k] = center_vel[k] = 0.0;

  DF          = true;
  MULTI       = true;

  dmass       = DMass;
  scalelength = A;
  scaleheight = H;

  expandh     = haloexp;
  expandd     = diskexp;

  random_gen.seed(SEED*numprocs + myid);
  SphericalModelTable::even = 0;

  halo = std::make_shared<SphericalModelTable>(filename1, DIVERGE, DIVERGE_RFAC);

  disk = std::make_shared<ExponentialDisk>(A, RDMAX);

  if (myid==0 && VFLAG & 1) {
    std::cerr << "DiskHalo: DIVERGE=" << DIVERGE
	 << " DIVERGE2=" << DIVERGE2
	 << " A=" << A
	 << " RDMAX=" << RDMAX
	 << " filename1=" << filename1
	 << " filename2=" << filename2
	 << "\n";
  }

  AddDisk::logarithmic   = true;
  AddDisk::Rmin          = RHMIN;
  AxiSymModel::numr      = 400;
  AxiSymModel::numj      = 400;
  AxiSymModel::gen_N     = 800;
  AxiSymModel::gen_itmax = 4000000;
  AxiSymModel::gen_rmin  = RHMIN;

  newmod = std::make_shared<AddDisk>(halo, disk, dmass*COMPRESSION); 
  halo2 = newmod->get_model();
  halo2->print_model("diskhalo.newmodel");
  halo2->setup_df(NUMDF, RA);
  if (myid==0 && VFLAG & 2) {
    char debugname[] = "df.debug";
    halo2->print_df(debugname);
  }

  if (myid==0) std::cout << "DF MADE" << std::endl;

  MPI_Barrier(MPI_COMM_WORLD);

  //
  // Generate "fake" profile
  //
  SphericalModelTable::even     = 0;
  
  halo3 = std::make_shared<SphericalModelTable>(filename2, DIVERGE2, DIVERGE_RFAC2);

  //
  // Packs fake density and mass model with target (real) potential
  // and reinitializes the model
  //
  std::vector<double> r2(RNUM);
  std::vector<double> d2(RNUM);
  std::vector<double> m2(RNUM);
  std::vector<double> p2(RNUM);
  
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
  
  halo3 = std::make_shared<SphericalModelTable>
    (RNUM, r2.data(), d2.data(), m2.data(), p2.data(), DIVERGE2, DIVERGE_RFAC2);
  halo3->setup_df(NUMDF, RA);
  if (VFLAG & 2) {
    halo3->print_model("diskhalo_model.multi");
    halo3->print_model_eval("diskhalo_model_eval.multi", RNUM*5);
    halo3->print_df("diskhalo_df.multi");
  }
    
  // Generate the multimass model
  //
  multi = std::make_shared<SphericalModelMulti>(halo2, halo3);
  multi -> gen_tolE = TOLE;
  multi -> set_itmax(ITMAX);
  if (ALLOW) multi -> allowNegativeMass();

  // For frequency computation
  //
  hDmin = log(max<double>(disk->get_min_radius(), RDMIN));
  hDmax = log(disk->get_max_radius());
  dRh   = (hDmax - hDmin)/nh;
  nhN   = vector<unsigned>(nh+1, 0  ); // Counts per bin
  nhD   = vector<double>  (nh+1, 0.0); // Mass per bin
  nhM   = vector<double>  (nh+1, 0.0); // Cumulative mass

  // For buffered ascii writes
  //
  bufcnt = 0;
}


DiskHalo::~DiskHalo()
{
  // Nothing
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

  expandh    = p.expandh;
  expandd    = p.expandd;

  disktableP = p.disktableP;
  disktableN = p.disktableN;
  epitable   = p.epitable;
  dv2table   = p.dv2table;
  asytable   = p.asytable;

  dP = p.dP;
  dR = p.dR;
  dZ = p.dZ;

  halotable = p.halotable;
  dr = p.dr;
  dc = p.dc;

  gen.seed(SEED+myid);

  DF    = p.DF;
  MULTI = p.MULTI;
  com   = p.com;
  cov   = p.cov;
  mtype = p.mtype;
  cmap  = p.cmap;
  Xmin  = p.Xmin;
  Xmax  = p.Xmax;

  bufcnt = 0;
}


double DiskHalo::disk_density(double R, double z)
{
  double q = 1.0/cosh(z/scaleheight);
  return disk_surface_density(R)*q*q*0.5/scaleheight;
}

double DiskHalo::disk_surface_density(double R)
{
  return dmass*disk->get_density(R);
}

void DiskHalo::set_halo(vector<Particle>& phalo, int nhalo, int npart)
{
  if (!MULTI) {
    string msg("DiskHalo::set_halo is only valid if MULTI is true");
    throw DiskHaloException(msg, __FILE__, __LINE__);
  }

  double rmin = max<double>(halo->get_min_radius(), RHMIN);
  double rmax = halo->get_max_radius();
  double mmin = halo->get_mass(rmin);
  double mtot = halo->get_mass(rmax);

  double pos[3], pos1[3], vel[3], vel1[3], massp, massp1;

				// Diagnostics
  double radmin1=1e30, radmax1=0.0, radmin, radmax, r;
  vector<double>   DD(nh+1, 0.0), DD0(nh+1);
  vector<unsigned> NN(nh+1, 0),   NN0(nh+1);

  if (myid==0 and VFLAG & 1) std::cout << std::endl 
				       << "     *****"
				       << "  rmin=" << rmin
				       << "  rmax=" << rmax
				       << "  mmin=" << mmin
				       << "  mtot=" << mtot
				       << std::endl;

  for (int k=0; k<3; k++) {
    pos[k] = pos1[k] = 0.0;
    vel[k] = vel1[k] = 0.0;
  }
  massp = massp1 = 0.0;

  double meanmass = (mtot - mmin)/nhalo;

  Particle p;
  Eigen::VectorXd ps(7);
  int ierr;

  unsigned int count1=0, count=0;
  unsigned int badms1=0, badms=0;

  for (int i=0; i<npart; i++) {

    do {
      ps = multi->gen_point(ierr);
      if (ierr) count1++;
    } while (ierr);
    
    if (ps[0]<0.0) badms1++;
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

    // Mass distribution in spherial shells
    if (r >= rmin) {
      unsigned indx = 1 + floor( (log(r) - hDmin)/dRh );
      if (indx>0 && indx<=nh) {
	NN[indx]++;
	DD[indx] += p.mass;
      }
    }

    radmin1 = min<double>(radmin1, r);
    radmax1 = max<double>(radmax1, r);
  }
  
  MPI_Reduce(&count1,  &count,  1, MPI_UNSIGNED, MPI_SUM, 0, MPI_COMM_WORLD);
  MPI_Reduce(&badms1,  &badms,  1, MPI_UNSIGNED, MPI_SUM, 0, MPI_COMM_WORLD);
  MPI_Reduce(&radmin1, &radmin, 1, MPI_DOUBLE,   MPI_MIN, 0, MPI_COMM_WORLD);
  MPI_Reduce(&radmax1, &radmax, 1, MPI_DOUBLE,   MPI_MAX, 0, MPI_COMM_WORLD);

  if (myid==0) std::cout << "     *****"
			 << "  min(r)=" << radmin 
			 << "  max(r)=" << radmax
			 << std::endl;
  
  if (myid==0 && count)
    std::cout << "DiskHalo::set_halo: " 
	      << count << " selection failures" << std::endl;
  
  if (myid==0 && badms)
    std::cout << "DiskHalo::set_halo: " 
	      << badms << " NEGATIVE masses" << std::endl;
  
  MPI_Allreduce(&massp1, &massp, 1, MPI_DOUBLE, MPI_SUM, MPI_COMM_WORLD);
  MPI_Allreduce(pos1,    pos,    3, MPI_DOUBLE, MPI_SUM, MPI_COMM_WORLD);
  MPI_Allreduce(vel1,    vel,    3, MPI_DOUBLE, MPI_SUM, MPI_COMM_WORLD);

  if (massp>0.0) {
    for (int k=0; k<3; k++) pos[k] /= massp;
    for (int k=0; k<3; k++) vel[k] /= massp;
  }
  
  if (myid==0) {
    std::cout << "     *****";
    std::cout << " (x, y, z)=(" << pos[0] 
	 << ", " << pos[1]
	 << ", " << pos[2] << ")" << std::endl;
    std::cout << "     *****";
    std::cout << " (u, v, w)=(" << vel[0]
	 << ", " << vel[1]
	 << ", " << vel[2] << ")" << std::endl;
  }


  double massn1 = 0.0, massn = massp, mfac = (mtot - mmin)/massp;
  for (auto &p : phalo) {
    p.mass *= mfac;
    massn1   += p.mass;
    if (com) for (int k=0; k<3; k++) p.pos[k] -= pos[k];
    if (cov) for (int k=0; k<3; k++) p.vel[k] -= vel[k];
  }
  MPI_Reduce(&massn1, &massn, 1, MPI_DOUBLE, MPI_SUM, 0, MPI_COMM_WORLD);
  

  if (myid==0) std::cout << "     *****"
			 << "  M(comp)=" << massp
			 << "  M(modl)=" << mtot - mmin
			 << "  M(meas)=" << massn
			 << std::endl;
  
  // Make dispersion vector
  //
  table_halo_disp();

  if (VFLAG & 1)
    std::cout << "Process " << myid << ": made " << phalo.size() << " particles"
	      << std::endl;

  MPI_Allreduce(&NN[0], &NN0[0], nh+1, MPI_UNSIGNED, MPI_SUM, MPI_COMM_WORLD);
  MPI_Allreduce(&DD[0], &DD0[0], nh+1, MPI_DOUBLE,   MPI_SUM, MPI_COMM_WORLD);
  for (unsigned i=0; i<=nh; i++) {
    nhN[i] += NN0[i];
    nhD[i] += DD0[i];
  }
}      

void DiskHalo::
set_halo_table_multi(vector<Particle>& phalo)
{
  if (!MULTI) {
    string msg("DiskHalo::set_halo is only valid if MULTI is true");
    throw DiskHaloException(msg, __FILE__, __LINE__);
  }

  double rmin = max<double>(halo->get_min_radius(), RHMIN);
  double rmax = halo->get_max_radius();
  double mmin = halo->get_mass(rmin);
  double mtot = halo->get_mass(rmax);

				// Diagnostics
  double radmin1=1e30, radmax1=0.0, radmin, radmax;
  vector<double>   DD(nh+1, 0.0), DD0(nh+1);
  vector<unsigned> NN(nh+1, 0),   NN0(nh+1);

  // Particle loop for existing particles
  //
  for (auto p : phalo) {

    double r = 0.0;
    for (int k=0; k<3; k++) r += p.pos[k]*p.pos[k];
    r = sqrt(r);

    // Mass distribution in spherical shells
    if (r >= rmin) {
      unsigned indx = 1 + floor( (log(r) - hDmin)/dRh );
      if (indx>0 && indx<=nh) {
	NN[indx]++;
	DD[indx] += p.mass;
      }
    }

    radmin1 = min<double>(radmin1, r);
    radmax1 = max<double>(radmax1, r);
  }
  
  MPI_Reduce(&radmin1, &radmin, 1, MPI_DOUBLE, MPI_MIN, 0, MPI_COMM_WORLD);
  MPI_Reduce(&radmax1, &radmax, 1, MPI_DOUBLE, MPI_MAX, 0, MPI_COMM_WORLD);

  if (myid==0) std::cout << "     *****"
			 << "  min(r)=" << radmin 
			 << "  max(r)=" << radmax
			 << std::endl;
  
  // Make dispersion vector
  //
  table_halo_disp();

  if (VFLAG & 1)
    std::cout << "Process " << myid << ": made " << phalo.size() << " particles"
	 << std::endl;

  MPI_Allreduce(&NN[0], &NN0[0], nh+1, MPI_UNSIGNED, MPI_SUM, MPI_COMM_WORLD);
  MPI_Allreduce(&DD[0], &DD0[0], nh+1, MPI_DOUBLE,   MPI_SUM, MPI_COMM_WORLD);
  for (unsigned i=0; i<=nh; i++) {
    nhN[i] += NN0[i];
    nhD[i] += DD0[i];
  }
}      

void DiskHalo::
set_halo_coordinates(vector<Particle>& phalo, int nhalo, int npart)
{
  const double tol = 1.0e-12;
  double rmin = max<double>(halo->get_min_radius(), RHMIN);
  double rmax = halo->get_max_radius();
  double mmin = halo->get_mass(rmin);
  double mtot = halo->get_mass(rmax);

  double r, phi, costh, sinth;
  double massp, massp1, pos[3], pos1[3];
				// Diagnostics
  double radmin1=1.0e30, radmax1=0.0, radmin, radmax;
  vector<double>   DD(nh+1, 0.0), DD0(nh+1);
  vector<unsigned> NN(nh+1, 0),   NN0(nh+1);

  if (myid==0 and VFLAG & 1) std::cout << "  rmin=" << rmin
				       << "  rmax=" << rmax
				       << "  mmin=" << mmin
				       << "  mtot=" << mtot
				       << std::endl;

  for (int k=0; k<3; k++) pos[k] = pos1[k] = 0.0;
  massp = massp1 = 0.0;

  Particle p;

  p.mass = (mtot - mmin)/nhalo;

  model = halo;

  MPI_Barrier(MPI_COMM_WORLD);

  for (int i=0; i<npart; i++) {
    targetmass = mmin + (mtot - mmin)*rndU(gen);
    
    r = zbrent(mass_func, rmin, rmax, tol);
    
    phi = 2.0*M_PI*rndU(gen);
    costh = 2.0*rndU(gen) - 1.0;
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

    // Mass distribution in spherial shells
    if (r >= rmin) {
      unsigned indx = 1 + floor( (log(r) - hDmin)/dRh );
      if (indx>0 && indx<=nh) {
	NN[indx]++;
	DD[indx] += p.mass;
      }
    }

    radmin1 = min<double>(radmin1, r);
    radmax1 = max<double>(radmax1, r);
  }

  MPI_Reduce(&radmin1, &radmin, 1, MPI_DOUBLE, MPI_MIN, 0, MPI_COMM_WORLD);
  MPI_Reduce(&radmax1, &radmax, 1, MPI_DOUBLE, MPI_MAX, 0, MPI_COMM_WORLD);
  
  if (myid==0) std::cout << "     *****"
			 << "  min(r)=" << radmin 
			 << "  max(r)=" << radmax
			 << std::endl;

  MPI_Allreduce(&massp1, &massp, 1, MPI_DOUBLE, MPI_SUM, MPI_COMM_WORLD);
  MPI_Allreduce(pos1, pos, 3, MPI_DOUBLE, MPI_SUM, MPI_COMM_WORLD);

  if (massp>0.0) {
    for (int k=0; k<3; k++) pos[k] /= massp;
  }

  if (myid==0) {
    std::cout << " (x, y, z)=(" << pos[0] 
	      << ", " << pos[1]
	      << ", " << pos[2] << ")" << std::endl;
  }
  

  if (com) {
    for (auto &p : phalo) {
      for (int k=0; k<3; k++) p.pos[k] -= pos[k];
    }
  }

  if (VFLAG & 1)
    std::cout << "Process " << myid << ": made " << phalo.size() << " particles"
	 << std::endl;

  MPI_Allreduce(&NN[0], &NN0[0], nh+1, MPI_UNSIGNED, MPI_SUM, MPI_COMM_WORLD);
  MPI_Allreduce(&DD[0], &DD0[0], nh+1, MPI_DOUBLE,   MPI_SUM, MPI_COMM_WORLD);
  for (unsigned i=0; i<=nh; i++) {
    nhN[i] += NN0[i];
    nhD[i] += DD0[i];
  }
}

void DiskHalo::
set_halo_table_single(vector<Particle>& phalo)
{
  const double tol = 1.0e-12;
  double rmin = max<double>(halo->get_min_radius(), RHMIN);
  double rmax = halo->get_max_radius();
  double mmin = halo->get_mass(rmin);
  double mtot = halo->get_mass(rmax);

				// Diagnostics
  double radmin1=1.0e30, radmax1=0.0, radmin, radmax;
  vector<double>   DD(nh+1, 0.0), DD0(nh+1);
  vector<unsigned> NN(nh+1, 0),   NN0(nh+1);

  if (myid==0 and VFLAG & 1) std::cout << "  rmin=" << rmin
				       << "  rmax=" << rmax
				       << "  mmin=" << mmin
				       << "  mtot=" << mtot;
  
  for (auto p : phalo) {

    double r = 0.0;
    for (int k=0; k<3; k++) r += p.pos[k] * p.pos[k];
    r = sqrt(r);

    // Mass distribution in spherial shells
    if (r >= rmin) {
      unsigned indx = 1 + floor( (log(r) - hDmin)/dRh );
      if (indx>0 && indx<=nh) {
	NN[indx]++;
	DD[indx] += p.mass;
      }
    }

    radmin1 = min<double>(radmin1, r);
    radmax1 = max<double>(radmax1, r);
  }

  MPI_Reduce(&radmin1, &radmin, 1, MPI_DOUBLE, MPI_MIN, 0, MPI_COMM_WORLD);
  MPI_Reduce(&radmax1, &radmax, 1, MPI_DOUBLE, MPI_MAX, 0, MPI_COMM_WORLD);
  
  if (myid==0) std::cout << "     *****"
			 << "  min(r)=" << radmin 
			 << "  max(r)=" << radmax
			 << std::endl;

  MPI_Allreduce(&NN[0], &NN0[0], nh+1, MPI_UNSIGNED, MPI_SUM, MPI_COMM_WORLD);
  MPI_Allreduce(&DD[0], &DD0[0], nh+1, MPI_DOUBLE,   MPI_SUM, MPI_COMM_WORLD);
  for (unsigned i=0; i<=nh; i++) {
    nhN[i] += NN0[i];
    nhD[i] += DD0[i];
  }
}

void DiskHalo::
set_disk_coordinates(vector<Particle>& pdisk, int ndisk, int npart)
{
  const double tol = 1.0e-12;
  double rmin = max<double>(disk->get_min_radius(), RDMIN);
  double rmax = disk->get_max_radius();
  double mmin = disk->get_mass(rmin);
  double mtot = disk->get_mass(rmax);

  double R, phi;
  double pos[3], pos1[3], massp, massp1;

  // Diagnostics
  //
  double radmin1=1.0e30, radmax1=0.0, radmin, radmax, r;
  vector<double>   DD(nh+1, 0.0), DD0(nh+1);
  vector<unsigned> NN(nh+1, 0),   NN0(nh+1);

  if (myid==0 and VFLAG & 1) std::cout << std::endl
				       << "     *****"
				       << "  rmin=" << rmin
				       << "  rmax=" << rmax
				       << "  mmin=" << mmin
				       << "  mtot=" << mtot
				       << std::endl;
  
  for (int k=0; k<3; k++) pos[k] = pos1[k] = 0.0;
  massp = massp1 = 0.0;

  Particle p;

  p.mass = dmass/ndisk;

  model = disk;

  for (int i=0; i<npart; i++) {
    targetmass = mmin + (mtot-mmin)*rndU(gen);
    R = zbrent(mass_func, rmin, rmax, tol);
    phi = 2.0*M_PI*rndU(gen);

    p.pos[0] = R*cos(phi);
    p.pos[1] = R*sin(phi);
    p.pos[2] = scaleheight*atanh(2.0*rndU(gen)-1.0);

    massp1 += p.mass;
    for (int k=0; k<3; k++) pos1[k] += p.mass*p.pos[k];
    
    pdisk.push_back(p);

    r = 0.0;
    for (int k=0; k<3; k++) r += p.pos[k] * p.pos[k];
    r = sqrt(r);

    // Mass distribution in spherial shells
    if (r >= rmin) {
      unsigned indx = 1 + floor( (log(r) - hDmin)/dRh );
      if (indx>0 && indx<=nh) {
	NN[indx]++;
	DD[indx] += p.mass;
      }
    }

    radmin1 = min<double>(radmin1, r);
    radmax1 = max<double>(radmax1, r);
  }


  MPI_Reduce(&radmin1, &radmin, 1, MPI_DOUBLE, MPI_MIN, 0, MPI_COMM_WORLD);
  MPI_Reduce(&radmax1, &radmax, 1, MPI_DOUBLE, MPI_MAX, 0, MPI_COMM_WORLD);

  if (myid==0) std::cout << "     *****"
			 << "  min(r)=" << radmin 
			 << "  max(r)=" << radmax 
			 << std::endl;

  MPI_Allreduce(&massp1, &massp, 1, MPI_DOUBLE, MPI_SUM, MPI_COMM_WORLD);
  MPI_Allreduce(pos1, pos, 3, MPI_DOUBLE, MPI_SUM, MPI_COMM_WORLD);

  if (massp>0.0) {
    for (int k=0; k<3; k++) pos[k] /= massp;
  }

  if (myid==0) {
    std::cout << "     *****";
    std::cout << " (x, y, z)=(" << pos[0] 
	      << ", " << pos[1]
	      << ", " << pos[2] << ")" << std::endl;
  }
  
  if (com) {
    for (auto &p : pdisk) {
      for (int k=0; k<3; k++) p.pos[k] -= pos[k];
    }
  }

  MPI_Allreduce(&NN[0], &NN0[0], nh+1, MPI_UNSIGNED, MPI_SUM, MPI_COMM_WORLD);
  MPI_Allreduce(&DD[0], &DD0[0], nh+1, MPI_DOUBLE,   MPI_SUM, MPI_COMM_WORLD);
  for (unsigned i=0; i<=nh; i++) {
    nhN[i] += NN0[i];
    nhD[i] += DD0[i];
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
    std::cout << "deri_pot: derivative of order " << n 
	      << " can't be calculated" << std::endl;
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
    std::cout << "deri_pot_disk: derivative of order " << n 
	      << " can't be calculated" << std::endl;
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
    std::cout <<"deri_pot_halo: derivative of order " << n 
	      << " can't be calculated" << std::endl;
    break;
  }
  return dP;
}

double DiskHalo::
epi(double xp, double yp, double zp)
{
  if (disktableP.size()==0) {
    std::cerr << "epi: must call table_disk first\n";
    MPI_Abort(MPI_COMM_WORLD, 100);
    exit(0);
  }

  // Azimuth
  //
  double phi = atan2(yp, xp);
  if (phi<0.0) phi = 2.0*M_PI + phi;

  int iphi1 = floor( phi/dP );
  iphi1 = std::min<int>(iphi1, NDP-1);
  int iphi2 = iphi1 + 1;
  if (iphi1==NDP-1) iphi2 = 0;	// Modulo 2Pi

  double cp[2], cr[2];

  cp[1] = (phi - dP*iphi1)/dP;
  cp[0] = 1.0 - cp[1];

  // Cylindrical radius
  //
  double lR = log(max<double>(RDMIN, sqrt(xp*xp + yp*yp)));
  int ir1 = floor( (lR - log(RDMIN))/dR );
  int ir2 = ir1 + 1;

  // Make sure that the grid position is good, otherwise truncate to
  // lowest good radius
  if (ir1 < nzepi) {
    ir1 = nzepi;
    ir2 = ir1 + 1; 

    cr[1] = 0.0;
    cr[0] = 1.0;
  }  else {
    ir1 = std::min<int>( ir1, NDR-2 );
    ir1 = std::max<int>( ir1, 0 );
    ir2 = ir1 + 1;
  
    cr[1] = (lR - log(RDMIN) - dR*ir1)/dR;
    cr[0] = 1.0 - cr[1];
  }

  double ans = 
    cp[0]*cr[0]* epitable(iphi1, ir1) +
    cp[0]*cr[1]* epitable(iphi1, ir2) +
    cp[1]*cr[0]* epitable(iphi2, ir1) +
    cp[1]*cr[1]* epitable(iphi2, ir2) ;
    
  if (ans>0.0) return sqrt(ans);
  else {

    // Move forward
    int j;
    bool ok = false;
    for (j=ir2; j<NDR/4; j++) {
      if (epitable(iphi1, j)<=0.0) continue;
      if (epitable(iphi2, j)<=0.0) continue;
      ok = true;
      ans = 
	cp[0]* epitable(iphi1, j) +
	cp[1]* epitable(iphi2, j) ;
    }

    ostringstream sout;
    sout << "epi_range_error." << RUNTAG << "." << myid;
    ofstream out(sout.str().c_str(), ios::app);

    if (ok) {

      out << "Process " << myid << " epi range error [forward]" << std::endl
	  << "     R="  << sqrt(xp*xp + yp*yp)    << std::endl
	  << "   Phi="  << phi                    << std::endl
	  << "   ir1="  << ir1 << "/" << NDR      << std::endl
	  << "   ir2="  << ir2 << "/" << NDR      << std::endl
	  << "    ans=" << ans                    << std::endl
	  << "    ep1=" << epitable(iphi1, ir1)   << std::endl
	  << "    ep2=" << epitable(iphi1, ir2)   << std::endl
	  << "    ep3=" << epitable(iphi2, ir1)   << std::endl
	  << "    ep4=" << epitable(iphi2, ir2)   << std::endl
	  << "      j=" << j
	  << "    ep5=" << epitable(iphi1, j)     << std::endl
	  << "    ep6=" << epitable(iphi2, j)     << std::endl
	  << std::endl;

      return sqrt(ans);

    } else {

      out << "Process " << myid << " epi range error [unresolved]" << std::endl
	  << "     R="  << sqrt(xp*xp + yp*yp)    << std::endl
	  << "  Rmax="  << RDMIN*exp(dR*NDR)      << std::endl
	  << "   del="  << (lR - log(RDMIN))/dR   << std::endl
	  << "   ir1="  << ir1 << "/" << NDR      << std::endl
	  << "   ir2="  << ir2 << "/" << NDR      << std::endl
	  << "    cr="  << cr[0] << ", " << cr[1] << std::endl
	  << "   Phi="  << phi                    << std::endl
	  << "    cp="  << cp[0] << ", " << cp[1] << std::endl
	  << "    ans=" << ans                    << std::endl
	  << "    ep1=" << epitable(iphi1, ir1)   << std::endl
	  << "    ep2=" << epitable(iphi1, ir2)   << std::endl
	  << "    ep3=" << epitable(iphi2, ir1)   << std::endl
	  << "    ep4=" << epitable(iphi2, ir2)   << std::endl
	  << std::endl;

      ans = cp[0]*epitable(iphi1, nzepi) + cp[1]*epitable(iphi2, nzepi);
      return sqrt(ans);

    }
  }
}


// For the disk, setting vz_disp according to the Jeans' equation
// solution, B+T equation 4-29c, where here the potential is the 
// composite basis representation of the total particle distribution. 

void DiskHalo::
table_disk(vector<Particle>& part)
{
  if (disktableP.size()) return;

  disktableP.resize(NDP);
  disktableN.resize(NDP);
  for (int i=0; i<NDP; i++) {
    disktableP[i].resize(NDR, NDZ);
    disktableN[i].resize(NDR, NDZ);
  }
    
  epitable.resize(NDP, NDR);
  dv2table.resize(NDP, NDR);
  asytable.resize(NDP, NDR);

  dP = 2.0*M_PI/NDP;

  double maxr1 = 0.0, maxz1 = 0.0;
  double maxr =  0.0, maxz  = 0.0;
  for (auto &p : part) {
    maxr1 = max<double>(maxr1, sqrt(p.pos[0]*p.pos[0] + p.pos[1]*p.pos[1]));
    maxz1 = max<double>(maxz1, fabs(p.pos[2]));
  }

  maxz1 = max<double>(scaleheight*SHFACTOR, maxz1);

  MPI_Allreduce(&maxr1, &maxr, 1, MPI_DOUBLE, MPI_MAX, MPI_COMM_WORLD);
  MPI_Allreduce(&maxz1, &maxz, 1, MPI_DOUBLE, MPI_MAX, MPI_COMM_WORLD);

  Xmin = cmap->r_to_x(RDMIN);
  Xmax = cmap->r_to_x(maxr );

  dR = (Xmax - Xmin)/(NDR-1);
  dZ = maxz/(NDZ-1);

  if (myid==0) {
    std::cout << std::endl
	      << "Table disk epi parameters:"  << std::endl
	      << "  RDMIN=" << RDMIN           << std::endl
	      << "   maxr=" << maxr            << std::endl
	      << "   maxz=" << maxz            << std::endl
	      << "   Xmin=" << Xmin            << std::endl
	      << "   Xmax=" << Xmax            << std::endl
	      << "   dR="   << dR              << std::endl
	      << "   dZ="   << dZ              << std::endl
	      << std::endl;
  }

  double dZ1 = maxz/(NDZ-1);

  double R, X, r, x, y, z, phi, fr, fz, fp, theta;
  double pot, dens, potl, dpr, dpt, dpp, dpdz;

				// Add no force if no component exists
  pot = fr = fz = fp = dens = potl = dpr = dpt = dpp = 0.0;

  Eigen::VectorXd workP (NDZ);
  Eigen::VectorXd workN (NDZ);
  Eigen::VectorXd workA (NDZ);
  Eigen::VectorXd workZ (NDZ);

  Eigen::VectorXd workR (NDR);
  Eigen::VectorXd workE (NDR);
  Eigen::VectorXd workQ (NDR);
  Eigen::VectorXd workQ2(NDR);
  Eigen::VectorXd workQ3(NDR);
  Eigen::VectorXd workQ4(NDR);
#ifdef LOGCHEBY
  Eigen::VectorXd workQ2log   (NDR);
  Eigen::VectorXd workQ2smooth(NDR);
#endif

  Eigen::MatrixXd workV (5, NDR);

				// For debugging
  Eigen::MatrixXd workD (6, NDR);
  Eigen::MatrixXd workDZ(7, NDZ);

				// Sum mass grid and make radial mesh
  unsigned nzcnt=0;
  vector<double> nrD(nh+1);
  for (nzero=0; nzero<=nh; nzero++) if (nhN[nzero] >= mh) break;

  if (nzero>nh) nzero=0;	// Not enough particles . . . 
  if (myid==0) std::cout << "Nzero=" << nzero << "/" << nh << std::endl;
  
  nzero = floor( (hDmin + nzero*dRh - log(RDMIN))/dR ) + 1;
  if (myid==0) std::cout << "Nzero=" << nzero << "/" << NDR << std::endl;

				// X grid in log radius for mass
  for (int n=0; n<=nh; n++) nrD[n] = hDmin + dRh*n;
  nhM[0] = nhD[0];		// Compute cumulative mass
  for (int n=1; n<=nh; n++) nhM[n] = nhD[n] + nhM[n-1];

				// Compute this table in parallel

  std::vector<int> ibeg(numprocs), iend(numprocs);

  // Default values.  If not reset, process computation will be skipped
  //
  std::fill(ibeg.begin(), ibeg.end(), 0);
  std::fill(iend.begin(), iend.end(), 0);

  // Do at most NDP computations
  //
  if (numprocs > NDP) {

    for (int i=0; i<NDP; i++) {
      ibeg[i] = i;
      iend[i] = i + 1;
    }

  }
  // Fewer than numprocs phi points
  //
  else {

    for (int i=0; i<numprocs; i++) {
      ibeg[i] = (i  )*NDP/numprocs;
      iend[i] = (i+1)*NDP/numprocs;
    }
    
  }

  if (myid==0) {
    std::cout << std::endl << " *** Processor phi angles *** " << std::endl;
    for (int i=0; i<numprocs; i++)
      std::cout << "# " << setw(3) << i << ": " 
		<< setw(10) << ibeg[i]
		<< setw(10) << iend[i]
		<< std::endl;
  }

  std::shared_ptr<Cheby1d> cheb, cheb2;

				// Test cumulative mass evaluation
  MonotCubicInterpolator monoT(nrD, nhM); 

  if (true and myid==0) {
    std::ofstream tout("mass.debug");
    if (tout) {
      for (int i=0; i<nrD.size(); i++) {
	tout << std::setw(16) << nrD[i]
	     << std::setw(16) << nhD[i]
	     << std::setw(16) << nhM[i]
	     << std::setw(16) << monoT.evaluate(nrD[i])
	     << std::endl;
      }
    } else {
      std::cout << "DiskHalo: could not open test file <mass.debug>"
		<< std::endl;
    }
  }

  // If numprocs>NDP, some processes will skip this loop
  //
  for (int i=ibeg[myid]; i<iend[myid]; i++) {

    phi = dP*i;

    for (int j=0; j<NDR; j++) {

      R = RDMIN*exp(dR*j);
      x = R*cos(phi);
      y = R*sin(phi);
      
      workR[j] = Xmin + dR*j;

				// For epicyclic frequency
				// 
      disk_eval(R, 0.0, phi, pot, fr, fz, fp);
				// Use monopole part of expansion here, only
      if (expandh)		//
	expandh->determine_fields_at_point(R, 0.5*M_PI, phi,
					   &dens, &potl, &dpr, &dpt, &dpp);

      
      workV(0, j) = log(RDMIN) + dR*j;
				
      if (use_mono) {
	// Use monopole approximation for dPhi/dr
	//
	workE[j] = monoT.evaluate(workV.row(0)[j])/(R*R);

      } else
	// Use basis evaluation (dPhi/dr)
	//
	workE[j]    = std::max<double>(-fr + dpr, 1.0e-20);

      workV(1, j) = disk_surface_density(R);

				// Sigma(R)*dPhi/dr*R
      workV(2, j) = workV(1, j)*workE[j]*R;

				// [1/r dPhi/dr]^{1/2} = Omega
      workQ[j]    = sqrt(workE[j]/R);

				// r^2*[1/r dPhi/dr]^{1/2} = r^2 * Omega
      workQ2[j]   = workQ[j] * R*R;

      workQ3[j]   = -fr;	// For testing only
      workQ4[j]   = dpr;	// For testing only

      if (i==0) {
	workD(4, j) = -fr;
	workD(5, j) = dpr;
      }

      for (int k=0; k<NDZ; k++) {

	z = workZ[k] = dZ1*k;

	r = sqrt(R*R + z*z) + std::numeric_limits<double>::min();

				// Positive
	disk_eval(R, z, phi, pot, fr, fz, fp);
	theta = acos(z/(r + std::numeric_limits<double>::min()));

	if (expandh)
	  expandh->determine_fields_at_point(R, theta, phi,
					     &dens, &potl, &dpr, &dpt, &dpp);

	dpdz = -fz + dpr*z/r + dpt*R*R/(r*r*r);

				// Debugging
	workDZ(0, k) = disk_density(R, z);
	workDZ(1, k) = -fz;
	workDZ(2, k) = dpr*z/r;
	workDZ(3, k) = dpt*R*R/(r*r*r);

	workP[k] = disk_density(R, z) * dpdz;

				// Negative
	z *= -1.0;

	disk_eval(R, z, phi, pot, fr, fz, fp);
	theta = acos(z/(r + std::numeric_limits<double>::min()));

	if (expandh)
	  expandh->determine_fields_at_point(R, theta, phi,
					     &dens, &potl, &dpr, &dpt, &dpp);

	dpdz  = -fz + dpr*z/r + dpt*R*R/(r*r*r);
	dpdz *= -1.0;		// No harm in changing sign

	workDZ(4, k) = -fz;
	workDZ(5, k) = dpr*z/r;
	workDZ(6, k) = dpt*R*R/(r*r*r);

	workN[k] = disk_density(R, z) * dpdz;
      }

				// Integrate positive
      if (use_spline) Splsum (workZ, workP, workA);
      else            Trapsum(workZ, workP, workA);

      for (int k=0; k<NDZ; k++)
	disktableP[i](j, k) = max<double>(workA[NDZ-1] - workA[k], std::numeric_limits<double>::min());

      if (i==ibeg[myid] && j==0 && VFLAG & 4) {
	ostringstream ofile;
	ofile << "intgr_disk_P." << RUNTAG << ".d" << myid;
	ofstream dout(ofile.str().c_str());
	for (int k=0; k<NDZ; k++) {
	  dout << setw(15) << workZ[k] 
	       << setw(15) << workP[k] 
	       << setw(15) << workA[k];
	  for (int q=0; q<7; q++) dout << setw(15) << workDZ(q, k);
	  dout << std::endl;
	}
	dout.close();
      }
				// Integrate negative
      if (use_spline) Splsum (workZ, workN, workA);
      else            Trapsum(workZ, workN, workA);

      for (int k=0; k<NDZ; k++)
	disktableN[i](j, k) = max<double>(workA[NDZ-1] - workA[k], std::numeric_limits<double>::min());

      if (i==ibeg[myid] && j==0 && VFLAG & 4) {
	ostringstream ofile;
	ofile << "intgr_disk_N." << RUNTAG << ".d" << myid;
	ofstream dout(ofile.str().c_str());
	for (int k=0; k<NDZ; k++) 
	  dout << setw(15) << workZ[k] 
	       << setw(15) << workN[k] 
	       << setw(15) << workA[k]
	       << "\n";
	dout.close();
      }

    }

    if (CHEBY) {
#ifdef LOGCHEBY
      // try a log formula due to crazy numbers
      for (int j=0; j<NDR; j++) {
	workQ2log[j] = log(workQ2[j]);
      }

      cheb = std::make_shared<Cheby1d>(workR, workQ2log, NCHEB);

      for (int j=0; j<NDR; j++) {
	workQ2smooth[j] = exp(cheb->eval(workR[j]));
      }
#else
      cheb = std::make_shared<Cheby1d>(workR, workQ2, NCHEB);
#endif
    }

				// Compute epicylic freqs
    for (int j=0; j<NDR; j++) {

      if (i==0) {
	if (CHEBY)
	  workD(0, j) += cheb->eval(workR[j]);
	else
	  workD(0, j) += workE[j];
      }

      if (CHEBY) {
#ifdef LOGCHEBY
	epitable(i, j) = drv2(workR[j], Eigen::VectorXd(workV.row(0)),
			      workQ2smooth);
#else
        epitable(i, j) = cheb->deriv(workR[j]);
#endif
      } else {
	epitable(i, j) = drv2(workR[j], Eigen::VectorXd(workV.row(0)),
			      workQ2);
      }

#ifdef ENFORCE_KAPPA
      {
	double om2 = workQ[j]*workQ[j];
	epitable(i, j) = std::max<double>(epitable(i, j), om2    );
	epitable(i, j) = std::min<double>(epitable(i, j), om2*4.0);
      }
#endif

      if (i==0) workD(1, j) = epitable(0, j);
      epitable(i, j) *= 2.0*workQ[j]/exp(2.0*workR[j]);
      if (i==0) workD(2, j) = epitable(0, j);
      if (i==0) workD(3, j) = epitable(0, j);
    }

				// Cylindrical Jeans' equations
    for (int j=0; j<NDR; j++) {
      double ep   = epitable(i, j);
      double sf   = workV(1, j);
      double vr   = 3.36*workV(1, j)*Q/sqrt(epitable(i, j));
      workV(4, j) = log(workV(1, j)*vr*vr);
    }

    if (CHEBY)
      cheb2 = std::make_shared<Cheby1d>(workV.row(0), workV.row(4), NCHEB);
  

    Eigen::VectorXd Z(NDR);

    Trapsum(workV.row(0), workV.row(2), Z);
    workV.row(3) = Z;

    for (int j=0; j<NDR; j++) {
      dv2table(i, j) = (workV(3, NDR-1) - workV(3, j))/workV(1, j);
      if (CHEBY) {
	asytable(i, j) = cheb2->deriv(workV(0, j));
      } else {
	asytable(i, j) = drv2(workV(0, j),
			      Eigen::VectorXd(workV.row(0)),
			      Eigen::VectorXd(workV.row(4)), 1);
      }
    }
  }

  MPI_Barrier(MPI_COMM_WORLD);	// Inactive processes will wait here

  // Update tables on all nodes
  //
  Eigen::VectorXd Z(NDR);
  for (int k=0; k<numprocs; k++) {
    for (int i=ibeg[k]; i<iend[k]; i++) {
      if (k == myid) Z = epitable.row(i);
      MPI_Bcast(Z.data(), NDR, MPI_DOUBLE, k, MPI_COMM_WORLD);
      if (k != myid) epitable.row(i) = Z;
      if (k == myid) Z = dv2table.row(i);
      MPI_Bcast(Z.data(), NDR, MPI_DOUBLE, k, MPI_COMM_WORLD);
      if (k != myid) dv2table.row(i) = Z; 
      if (k == myid) Z = asytable.row(i);
      MPI_Bcast(Z.data(), NDR, MPI_DOUBLE, k, MPI_COMM_WORLD);
      if (k != myid) asytable.row(i) = Z; 
      MPI_Bcast(disktableP[i].data(), NDR*NDZ, MPI_DOUBLE, k, MPI_COMM_WORLD);
      MPI_Bcast(disktableN[i].data(), NDR*NDZ, MPI_DOUBLE, k, MPI_COMM_WORLD);
    }
  }

  // Compute minimum >zero index
  //
  nzepi = 0;
  for (int i=0; i<NDP; i++) {
    while (epitable(i, nzepi) <= 0.0 and nzepi<NDR) nzepi++;
  }
  if (myid==0) std::cout << "NZEPI=" << nzepi << "/" << NDR << std::endl;

  // Compute velocity dispersion scaling based on circular velocity at
  // a scale length
  //
  sigma0 = 1.0;
  if (Q <= 0.0) {

    if (SIG0 <= 0.0)
      throw std::runtime_error("DiskHalo: radial dispersion must be positive.");

    if (XI <= 0.0)
      throw std::runtime_error("DiskHalo: dispersion axis ratio must be positive.");

    sigma0 = SIG0*v_circ(scalelength, 0.0, 0.0);
  }

  // For debugging the solution
  //
  if (myid==0 && expandh && VFLAG & 4) {
    ostringstream sout;
    sout << "ep_test." << RUNTAG;
    ofstream out(sout.str().c_str());
    out.setf(ios::scientific);
    out.precision(4);

    const double DR = 0.01;

    double r, r1, r2, deriv, deriv2, vrq0, vrq1, rho, lhs, rhs;

    for (int j=0; j<NDR; j++) {
      r = RDMIN*exp(dR*j);
      r1 = r*(1.0 + dR*DR);
      r2 = r*(1.0 - dR*DR);

      rho = halo->get_density(r);

      deriv = (get_disp(0.0, r1, 0.0)*halo->get_density(r1) - 
	       get_disp(0.0, r2, 0.0)*halo->get_density(r2) ) /	(r1 - r2);
      
      lhs = halo->get_mass(r);
      rhs = -r*r*deriv/rho;

      if (CHEBY) {
#ifdef LOGCHEBY
	Eigen::VectorXd X(workV.row(0));
	deriv2 = drv2(workR[j], X, workQ2smooth);
#else
        deriv2 = cheb->deriv(workQ2[j]);
#endif
      } else {
        deriv2 = drv2(workR[j], workR, workQ2);
      }	

      vrq0 = 3.36*dmass*disk->get_density(r)*Q/epi(r, 0.0, 0.0);
      vrq1 = 3.36*dmass*disk->get_density(r)*Q/sqrt(epitable(0, j));

      out << setw(14) << r			// #1
	  << setw(14) << epitable(0, j)		// #2
	  << setw(14) << workE[j]		// #3
	  << setw(14) << workQ[j]		// #4
	  << setw(14) << workQ2[j]		// #5
	  << setw(14) << workQ3[j]		// #6
	  << setw(14) << workQ4[j]		// #7
	  << setw(14) << deriv2                 // #8
	  << setw(14) << vrq0                   // #9
	  << setw(14) << vrq1                   // #10
	  << setw(14) << v_circ(r, 0.0, 0.0)    // #11
	  << setw(14) << workD(0, j)		// #12  dV(tot)/dR
	  << setw(14) << workD(1, j)		// #13  d^2V(tot)/dlnR
	  << setw(14) << workD(2, j)		// #14  d^2V(tot)/dlnR + 3V(tot)
	  << setw(14) << workD(3, j)		// #15  kappa^2
	  << setw(14) << workD(4, j)		// #16  dV(disk)/dR
	  << setw(14) << workD(5, j)		// #17  dV(halo)/dR
	  << setw(14) << rho			// #18
	  << setw(14) << deriv			// #19
	  << setw(14) << lhs			// #20
	  << setw(14) << rhs			// #21
	  << setw(14) << lhs - rhs		// #22
	  << setw(14) << odd2(log(r), nrD, nhM, 1) // #23  Enclosed mass
	  << setw(14) << epi(r, 0.0, 0.0)	// #24  Epi routine
	  << std::endl;
    }

    ostringstream sout2;
    sout2 << "epitable." << RUNTAG;
    ofstream dump(sout2.str().c_str());
    for (int i=0; i<NDP; i++) {
      phi = dP*i;
      for (int j=0; j<NDR; j++)
	dump << setw(18) << phi 
	     << setw(18) << RDMIN*exp(dR*j)
	     << setw(18) << epitable(i, j)
	     << std::endl;
      dump << std::endl;
    }
    dump.close();

    dump.open("dv2table.dump");
    for (int i=0; i<NDP; i++) {
      phi = dP*i;
      double c = cos(phi), s = sin(phi);
      for (int j=0; j<NDR; j++) {
	double rr = RDMIN*exp(dR*j);
	double vc = v_circ(rr*c, rr*s, 0.0);
	double vr = vr_disp2(rr*c, rr*s, 0.0);
	dump << setw(18) << phi 			// 1 Phi
	     << setw(18) << rr				// 2 R
	     << setw(18) << dv2table(i, j)		// 3 
	     << setw(18) << asytable(i, j)		// 4 d log(rho*vr^2)/ dlog(R)
	     << setw(18) << asytable(i, j)*vr/(vc*vc)	// 5 ac/vc^2
	     << setw(18) << workV(4, j)			// 6 log(rho*vr^2)
	     << setw(18) << vr				// 7 vr^2
	     << setw(18) << vc*vc;			// 8 vc^2
	if (CHEBY) dump << setw(18) << cheb2->eval(workV(0, j));
	dump << std::endl;
      }
      dump << std::endl;
    }
    dump.close();
  }
    
  if (myid==0 && VFLAG & 4) {
    ostringstream sout;
    sout << "ep_disk." << RUNTAG;
    ofstream out(sout.str().c_str());
    out.setf(ios::scientific);
    out.precision(8);

    for (int i=0; i<=NDP; i++) {
      phi = dP*i;
      out << "# i=" << " phi=" << phi << ", " << phi*180.0/M_PI << std::endl;
      for (int j=0; j<NDR; j++) {
	R = RDMIN*exp(dR*j);
	x = R*cos(phi);
	y = R*sin(phi);
	out << setw(18) << x
	    << setw(18) << y
	    << setw(18) << epitable(i%NDP, j) << std::endl;
      }
      out << std::endl;
    }
    out << std::flush;
    out.close();
  }

  if (myid==0 && VFLAG & 4) {
    ostringstream sout;
    sout << "table_disk." << RUNTAG;
    ofstream out(sout.str().c_str());
    out.setf(ios::scientific);
    out.precision(8);
    
    for (int j=0; j<NDR; j++) {
      for (int k=0; k<NDZ; k++) {
	out << setw(18) << RDMIN*exp(dR*j)
	    << setw(18) << dZ*k
	    << setw(18) << disktableP[0](j, k)
	    << setw(18) << disktableN[0](j, k) 
	    << setw(18) << disk_density(RDMIN*exp(dR*j), dZ*k)
	    << std::endl;
      }
      out << std::endl;
    }
  }

  if (myid==0) std::cout << "[table] " << std::flush;
}



// Closure using the 2nd moment of the cylindrical CBE
//  
double DiskHalo::vp_disp2(double xp, double yp, double zp)
{
  double R     = sqrt(xp*xp + yp*yp) + std::numeric_limits<double>::min();
  double vc    = v_circ(xp, yp, zp);
  double omp   = vc/R;
  double kappa = epi(xp, yp, zp);
  
				// Bounds limit
  double fraction = kappa*kappa/(4.0*omp*omp);
  if (fraction > 1.0)  fraction = 1.0;
  if (fraction < 0.25) fraction = 0.25;

  return vr_disp2(xp, yp, zp) * fraction;
}


// Interpolate from Jeans' solution integral table
//
double DiskHalo::vz_disp2(double xp,double yp, double zp)
{
  if (disktableP.size()==0) {
    std::cerr << "DiskHalo::vz_disp2: must call table_disk first\n";
    MPI_Abort(MPI_COMM_WORLD, 100);
    exit(0);
  }

  double R, lR, phi, cp[2], cr[2], cz[2], resvd, zz;
  int iphi1, iphi2, ir1, ir2, iz1, iz2;

				// Azimuth
  phi = atan2(yp, xp);
  if (phi<0.0) phi = 2.0*M_PI + phi;
  iphi1 = floor( phi/dP );
  iphi1 = min<int>( iphi1, NDP-1 );
  if (iphi1==NDP-1) iphi2 = 0;
  else iphi2 = iphi1+1;
  
  cp[1] = (phi - dP*iphi1)/dP;
  cp[0] = 1.0 - cp[1];
				// Cylindrical radius
  R  = max<double>(sqrt(xp*xp + yp*yp), RDMIN);
  lR = log(R);
  ir1 = (int)( (lR - log(RDMIN))/dR );
  ir1 = min<int>( ir1, NDR-2 );
  ir1 = max<int>( ir1, 0 );
  ir2 = ir1 + 1;

  cr[1] = (lR - log(RDMIN) - dR*ir1)/dR;
  cr[0] = 1.0 - cr[1];

				// Zed
  zz = fabs(zp);
  iz1 = floor( zz/dZ );
  iz1 = min<int>( iz1, NDZ-2 );
  iz2 = iz1 + 1;

  cz[1] = (zz - dZ*iz1)/dZ;
  cz[0] = 1.0 - cz[1];

  if (zp > 0.0) {

    resvd = 

      cp[0]*cr[0]*cz[0] * disktableP[iphi1](ir1, iz1)  +
      cp[0]*cr[0]*cz[1] * disktableP[iphi1](ir1, iz2)  +
      cp[0]*cr[1]*cz[0] * disktableP[iphi1](ir2, iz1)  +
      cp[0]*cr[1]*cz[1] * disktableP[iphi1](ir2, iz2)  +

      cp[1]*cr[0]*cz[0] * disktableP[iphi2](ir1, iz1)  +
      cp[1]*cr[0]*cz[1] * disktableP[iphi2](ir1, iz2)  +
      cp[1]*cr[1]*cz[0] * disktableP[iphi2](ir2, iz1)  +
      cp[1]*cr[1]*cz[1] * disktableP[iphi2](ir2, iz2)  ;

  } else {

    resvd = 

      cp[0]*cr[0]*cz[0] * disktableN[iphi1](ir1, iz1)  +
      cp[0]*cr[0]*cz[1] * disktableN[iphi1](ir1, iz2)  +
      cp[0]*cr[1]*cz[0] * disktableN[iphi1](ir2, iz1)  +
      cp[0]*cr[1]*cz[1] * disktableN[iphi1](ir2, iz2)  +

      cp[1]*cr[0]*cz[0] * disktableN[iphi2](ir1, iz1)  +
      cp[1]*cr[0]*cz[1] * disktableN[iphi2](ir1, iz2)  +
      cp[1]*cr[1]*cz[0] * disktableN[iphi2](ir2, iz1)  +
      cp[1]*cr[1]*cz[1] * disktableN[iphi2](ir2, iz2)  ;
  }

  double dens = disk_density(R, zp);
  if (dens>0.0) resvd /= dens;

  return resvd;
}

// Constant Toomre Q
//
double DiskHalo::vr_disp2(double xp, double yp, double zp)
{
  double r = sqrt(xp*xp + yp*yp);

  if (Q <= 0.0) {
    double smth = 0.25*scaleheight;
    return sigma0*sigma0*exp(-sqrt(r*r + smth*smth)/scalelength);
  } else {
    if (r > 10.0*scalelength) return 0.0;
    double sigmar = 3.36*disk_surface_density(r)*Q/epi(xp, yp, zp);
    return sigmar*sigmar;
  }
}

// Asymmetric drift equation: returns v_a*(v_a - 2*v_c)/sigma_rr^2
//
double DiskHalo::a_drift(double xp, double yp, double zp)
{
				// sigma_r and sigma_p
  double vr2 = vr_disp2(xp, yp, zp);
  double vp2 = vp_disp2(xp, yp, zp);

				// Azimuth
  double phi = atan2(yp, xp);
  if (phi<0.0) phi = 2.0*M_PI + phi;

  int iphi1 = floor( phi/dP );
  iphi1 = std::min<int>(iphi1, NDP-1);
  int iphi2 = iphi1 + 1;
  if (iphi1==NDP-1) iphi2 = 0; // Modulo 2Pi

  double cp[2], cr[2];

  cp[1] = (phi - dP*iphi1)/dP;
  cp[0] = 1.0 - cp[1];
				// Cylindrical radius
  double R  = max<double>(sqrt(xp*xp + yp*yp), RDMIN);
  double lR = log(R);
  int ir1 = floor( (lR - log(RDMIN))/dR );
  ir1 = min<int>( ir1, NDR-2 );
  ir1 = max<int>( ir1, 0 );
  int ir2 = ir1 + 1;
  
  cr[1] = (lR - log(RDMIN) - dR*ir1)/dR;
  cr[0] = 1.0 - cr[1];
  
  // DEBUG
  if (1) {
    if (cp[1]>1.0 || cp[1] < 0.0 ||
	cp[0]>1.0 || cp[0] < 0.0 )
      std::cerr << "DiskHalo::a_drift: phi=" << phi
	   << ", cp=(" << cp[0] << ", " << cp[1] << ")" << std::endl;
    if (cr[1]>1.0 || cr[1] < 0.0 ||
	cr[0]>1.0 || cr[0] < 0.0 )
      std::cerr << "DiskHalo::a_drift: R=" << R
	   << ", cr=(" << cr[0] << ", " << cr[1] << ")" << std::endl;
  }
  // END DEBUG
  
  double ret = 0.0;
  if (vr2>0.0) ret = 1.0 - vp2/vr2;

  return ret + 
    cp[0]*cr[0] * asytable(iphi1, ir1) +
    cp[0]*cr[1] * asytable(iphi1, ir2) +
    cp[1]*cr[0] * asytable(iphi2, ir1) +
    cp[1]*cr[1] * asytable(iphi2, ir2) ;
}


// Analytic rotation curve
//
double DiskHalo::v_circ(double xp, double yp, double zp)
{
  double R = sqrt(xp*xp + yp*yp);
  double vcirc2 = R*deri_pot(xp, yp, 0.0, 1);

				// Sanity check
  if (vcirc2<=0.0) {
    if (VFLAG & 8)
      std::cout << "DiskHalo::v_circ: circular velocity out of bounds, R="
		<< R << "  v_circ2=" << vcirc2 << std::endl;
    vcirc2 = 1.0e-20;
  }

  return sqrt(vcirc2);
}


void DiskHalo::
set_vel_disk(vector<Particle>& part)
{
  if (!expandd) {
    if (myid==0) std::cout << "[no disk particles] ";
    return;
  }

  double vvZ, vvR, vvP;
  double maxVR=-1.0e20, RVR=1e20;
  double maxVP=-1.0e20, RVP=1e20;
  double maxVZ=-1.0e20, RVZ=1e20;
  double vz, vr, vp, R, x, y, z, ac, vc, va, as, ad;
  double vel[3], vel1[3], massp, massp1;
  unsigned num_oob = 0;

  for (int k=0; k<3; k++) vel[k] = vel1[k] = 0.0;
  massp = massp1 = 0.0;

				// Better to make a 2-d table
  table_disk(part);
  
				// Debugging
  ofstream out;
  if (VFLAG & 4) {
    std::ostringstream sout;
    sout << "test_vel." << RUNTAG << "." << myid;
    out.open(sout.str().c_str());
    if (type==DiskHalo::Epicyclic) 
      out << "# " << std::right
	  << std::setw(12) << "R |"    << std::setw(14) << "z |"   << std::setw(14) << "v_circ |"
	  << std::setw(14) << "v_R |" << std::setw(14) << "v_phi |" << std::setw(14) << "v_z |"
	  << std::setw(14) << "R1 |" <<  std::setw(14) << "X |" << std::setw(14) << "kappa |";
    else
      out << "# " << std::right
	  << std::setw(12) << "R |"    << std::setw(14) << "z |"   << std::setw(14) << "v_circ |"
	  << std::setw(14) << "v_T |" << std::setw(14) << "drift |" << std::setw(14) << "kappa |"
	  << std::setw(14) << "v_R |" << std::setw(14) << "v_phi |" << std::setw(14) << "v_z |"
	  << std::setw(14) << "vv_R |" << std::setw(14) << "vv_phi |" << std::setw(14) << "vv_z |"
	  << std::setw(14) << "v_x |" << std::setw(14) << "v_y |";
    out << std::endl;
  }


  for (auto &p : part) {
				// From solution to Jeans' equations in
				// cylindrical coordinates
    x = p.pos[0];
    y = p.pos[1];
    z = p.pos[2];

    R = sqrt(x*x + y*y) + std::numeric_limits<double>::min();

    vvZ = vz_disp2(x, y, z);
    vvR = vr_disp2(x, y, z);

    if (type == Jeans)
      vvP = vvR/(XI*XI);
    else
      vvP = vp_disp2(x, y, z);
				 // For safety; should only be a problem
				 // on extrapolating the range
    vvZ = std::max<double>(vvZ, std::numeric_limits<double>::min());
    vvR = std::max<double>(vvR, std::numeric_limits<double>::min());
    vvP = std::max<double>(vvP, std::numeric_limits<double>::min());
    
    if (maxVZ < vvZ) {
      maxVZ = vvZ;
      RVZ   = R;
    }
    if (maxVR < vvR) {
      maxVR = vvR;
      RVR   = R;
      if (VFLAG & 8)
	std::cout << "maxVR: vvR = " << vvR
		  << " x=" << x << " y=" << y
		  << " epi=" << epi(x, y, 0.0)
		  << " sig=" << disk_surface_density(R)
		  << std::endl;

    }
    if (maxVP < vvP) {
      maxVP = vvP;
      RVP   = R;
    }

    // Circular velocity
    vc   = v_circ(x, y, z);

    // No asymmetric drift correction by default
    ac = 0.0;

    switch (type) {
    case Asymmetric:
      // Asymmetric drift correction
      ad = a_drift(x, y, z);
      as = 1 + vvR*ad/(vc*vc);

      if (as > 0.0 and not std::isnan(as))
	ac = vc*(1.0-sqrt(as));
      else {
	if (as<0.0 or std::isnan(as)) {
	  ac = vc;
	  num_oob++;
	}
	if (VFLAG & 8) {
	  int op = std::cout.precision(3);
	  std::cout << "ac oob:"
		    << " as="   << std::setw(10) << as 
		    << ", R="   << std::setw(10) << R
		    << ", ac="  << std::setw(10) << ac
		    << ", ad="  << std::setw(10) << ad
		    << ", vc="  << std::setw(10) << vc
		    << ", vvR=" << std::setw(10) << vvR
		    << std::endl;
	  std::cout.precision(op);
	}
      }

    case Jeans:
      va = max<double>(vc - ac, std::numeric_limits<double>::min());
     
      vz   = rndN(gen)*sqrt(std::max<double>(vvZ, std::numeric_limits<double>::min()));
      vr   = rndN(gen)*sqrt(std::max<double>(vvR, std::numeric_limits<double>::min()));
      vp   = rndN(gen)*sqrt(std::max<double>(vvP, std::numeric_limits<double>::min()));
      
      {
	double omp   = vc/R;
	double kappa = epi(x, y, z);
	double vp2   = vc*vc +	// From radial cylindrical Jeans using
				// epicyclic closure
	  vvR * (1.0 - kappa*kappa/(4.0*omp*omp) - 2.0*R/scalelength);
	if (vp2 >= 0.0) {
	  vp += sqrt(vp2);
	} else {
	  num_oob++;
	}
      }

      if (out) 
	out << std::setw(14) << R   << std::setw(14) << z   << std::setw(14) << vc
	    << std::setw(14) << va  << std::setw(14) << ac  << std::setw(14) << epi(x, y, z)
	    << std::setw(14) << vr  << std::setw(14) << vp  << std::setw(14) << vz
	    << std::setw(14) << vvR << std::setw(14) << vvP << std::setw(14) << vvZ
	    << std::setw(14) << vr*x/R - vp*y/R
	    << std::setw(14) << vr*y/R + vp*x/R
	    << std::endl;
      break;
      
    case Epicyclic:
      /*
	Epicyclic theory provides the x position relative to the guiding
	center for an arbitrary amplitude X for phase alpha = kappa*t:

      		x     = X cos(alpha)
		dx/dt = -kappa X sin(alpha)

	The phase averaged radial velocity at the guiding center is then:

      		<dx/dt> = 0, <(dx/dt)^2> = kappa*kappa*X*X/2

	Choose X by equating <(dx/dt)^2> = \sigma^2_r:

      		<X>   = 0
		<X^2> = 2*\sigma^2_r/(kappa*kappa)

	Strictly speaking this is not correct because the contribution
	to \sigma^2_r comes from many guiding centers.  So this
	equivlance probably over estimates the second moment in X.
	Choose X ~ normal with 0 mean and variance <X^2>
     */
      {
				// The normal variant
	double Xampl = rndN(gen);	
				// The cylindrical polar angle
	double phi   = atan2(y, x);
				// The radial phase (kappa*t)
	double alpha = 2.0*M_PI*rndU(gen);

				// Initial guess for iteration uses
				// present positions
	double kappa = epi(x, y, z);
	double X     = sqrt(2.0*Xampl*Xampl*vvR/(kappa*kappa));

				// Iterate to get values at guiding
				// center
	double Xl, x1, y1, R1, Omg;
	int cnt = 0;
	for (int i=0; i<10; i++) {
	  Xl      = X;
				// Guiding center estimate
	  R1      = R - X*cos(alpha);
				// x,y positions w.r.t. this guiding center
	  x1      = R1*cos(phi);
	  y1      = R1*sin(phi);
				// Epicylic freq at guiding center
	  kappa   = epi(x1, y1, z);
	
				// New amplitude
	  X       = sqrt(2.0*Xampl*Xampl*vr_disp2(x1, y1, z)/(kappa*kappa));

	  if (fabs((X-Xl)/Xl)<1.0e-6) break;
	  cnt++;
	}
	if (cnt>=100) {
	  std::cerr << "OOPS" << std::endl;
	}
				// Aximuthal freq at guiding center
	Omg  = v_circ(x1, y1, z)/R1;

				// Compute the final velocities
	vz   = rndN(gen)*sqrt(std::max<double>(vvZ, std::numeric_limits<double>::min()));
	vr   = -kappa*X*sin(alpha);
	vp   = Omg*R1 - 2.0*Omg*X*cos(alpha);
    
	if (out) 
	  out << std::setw(14) << R   << std::setw(14) << z   << std::setw(14) << Omg*R1
	      << std::setw(14) << vr  << std::setw(14) << vp  << std::setw(14) << vz
	      << std::setw(14) << R1  << std::setw(14) << X   << std::setw(14) << kappa
	      << std::endl;
      }
      break;
    }

    p.vel[0] = vr*x/R - vp*y/R;
    p.vel[1] = vr*y/R + vp*x/R;
    p.vel[2] = vz;

    massp1 += p.mass;
    for (int k=0; k<3; k++) vel1[k] += p.mass*p.vel[k];
  }

  MPI_Allreduce(&massp1, &massp, 1, MPI_DOUBLE, MPI_SUM, MPI_COMM_WORLD);
  MPI_Allreduce(vel1,    vel,    3, MPI_DOUBLE, MPI_SUM, MPI_COMM_WORLD);

  if (VFLAG & 1)
    std::cout << "vel per node [" << myid << "] mass=" << massp1
	      << " vx=" << vel1[0]/massp1
	      << " vy=" << vel1[1]/massp1
	      << " vz=" << vel1[2]/massp1 << std::endl;

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

    MPI_Reduce(MPI_IN_PLACE, &num_oob, 1, MPI_UNSIGNED, MPI_SUM, 0, MPI_COMM_WORLD);

    std::cout << "     *****"
	      << " (u, v, w)=(" << vel[0] 
	      << ", " << vel[1]
	      << ", " << vel[2] << ")" << std::endl
	      << "     *****"
	      <<  " maxVZ=" << maxVZ << " (" << RVZ << ")"
	      << ", maxVR=" << maxVR << " (" << RVR << ")"
	      << ", maxVP=" << maxVP << " (" << RVP << ")"
	      << std::endl;

    if (type == Jeans)
      std::cout << "     *****"
		<< " # Jeans' overrides=" << num_oob << std::endl;
    if (type == Asymmetric)
      std::cout << "     *****"
		<< " # adrift overrides=" << num_oob << std::endl;
  } else {
    MPI_Send(&maxVZ, 1, MPI_DOUBLE, 0, 224, MPI_COMM_WORLD);
    MPI_Send(&RVZ,   1, MPI_DOUBLE, 0, 225, MPI_COMM_WORLD);
    MPI_Send(&maxVR, 1, MPI_DOUBLE, 0, 226, MPI_COMM_WORLD);
    MPI_Send(&RVR,   1, MPI_DOUBLE, 0, 227, MPI_COMM_WORLD);
    MPI_Send(&maxVP, 1, MPI_DOUBLE, 0, 228, MPI_COMM_WORLD);
    MPI_Send(&RVP,   1, MPI_DOUBLE, 0, 229, MPI_COMM_WORLD);

    MPI_Reduce(&num_oob, 0, 1, MPI_UNSIGNED, MPI_SUM, 0, MPI_COMM_WORLD);
  }

  if (cov) {
    for (auto &p : part) {
      for (int k=0; k<3; k++) p.vel[k] -= vel[k];
    }
  }

}

void DiskHalo::table_halo_disp()
{
  if (!MULTI) {
    string msg("DiskHalo::table_halo_disp is only valid if MULTI is true");
    throw DiskHaloException(msg, __FILE__, __LINE__);
  }

  if (halotable.cols() == NHR && halotable.rows() == 3) return;
  
  halotable.resize(3, NHR);
  
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
  
  halotable.setZero();
  
  const int nlq = 400;
  LegeQuad lq(nlq);

  int ibeg = (int)( myid*(double)NHR/numprocs );
  int iend = (int)( (myid+1)*(double)NHR/numprocs );

  if (myid==numprocs-1) iend = NHR;

  for (int i=ibeg; i<iend; i++) {
    if (LOGR)
      halotable(0, i) = r = exp(log(rmin) + dr*i);
    else
      halotable(0, i) = r = rmin + dr*i;

    pot = halo2->get_pot(r);

    for (int n=0; n<nlq; n++) {
      E = pot + (Emax - pot)*lq.knot(n);
      v2 = 2.0*(E - pot);
      if (v2<0.0) v2 = 0.0;
      v = sqrt(v2);
      fac = (Emax - pot)*lq.weight(n)*halo2->distf(E, 0.5);
				// velocity dispersion
      halotable(1, i) += fac * v * v2;
				// density
      halotable(2, i) += fac * v;
    }

    if (halotable(2, i)>0.0) halotable(1, i) /= halotable(2, i);
    else halotable(1, i) = 0.0;
  }


  // Update tables on all nodes
  //
  MPI_Allreduce(MPI_IN_PLACE, halotable.data(), halotable.size(), MPI_DOUBLE, 
		MPI_SUM, MPI_COMM_WORLD);

  //
  // DEBUG output
  //
  if (myid==0 && VFLAG & 4) {
    ostringstream sout;
    sout << "disp_halo." << RUNTAG;
    ofstream out(sout.str().c_str());
    out.setf(ios::scientific);
    out.precision(2);
    
    for (int k=0; k<NHR; k++)
      out << std::setw(18) << halotable(0, k)
	  << std::setw(18) << halotable(1, k)
	  << std::setw(18) << 4.0*M_PI*halotable(2, k)
	  << std::setw(18) << halo->get_density(halotable(0, k)) << std::endl;
  }

}

void DiskHalo::
table_halo(std::vector<Particle>& part)
{
  if (halotable.cols() == NHR) return;

  if (true) {
    constexpr int ngrid {2000};
    std::ofstream out("table_halo.test");
    double rmin = halo->get_min_radius();
    double rmax = halo->get_max_radius();
    double dlog = (log(rmax) - log(rmin))/(ngrid - 1);
    double theta = acos(0.0);
    double dens, potl, dpr, dpt, dpp;
      
    for (int i=0; i<ngrid; i++) {
      double r = rmin*exp(dlog*i);
      expandh->determine_fields_at_point(r, theta, 0.0,
					 &dens, &potl, &dpr, &dpt, &dpp);
      out << std::setw(16) << r
	  << std::setw(16) << dens
	  << std::setw(16) << potl
	  << std::setw(16) << dpr
	  << std::setw(16) << dpt
	  << std::setw(16) << dpt
	  << std::setw(16) << dpp;
      
      expandh->determine_fields_at_point(r, 0.0, 0.0,
					 &dens, &potl, &dpr, &dpt, &dpp);

      out << std::setw(16) << dens
	  << std::setw(16) << potl
	  << std::setw(16) << dpr
	  << std::setw(16) << dpt
	  << std::setw(16) << dpt
	  << std::setw(16) << dpp
	  << std::endl;
    }
  }
  
  halotable.resize(NHT, NHR);
  
  dc = 2.0/(NHT-1);
  double r2, maxr = 0.0, maxr1 = 0.0;

  const double rpct = 0.01;	// Percentile for average top radius
  const int nret = std::floor(rpct*part.size());

  std::vector<double> racc, rtop(nret);
  for (auto &p : part) {
    r2 = 0.0;
    for (int k=0; k<3; k++) r2 += p.pos[k]*p.pos[k];
    racc.push_back(sqrt(r2));
  }
  
  std::partial_sort_copy(
			 std::begin(racc), std::end(racc),
			 std::begin(rtop), std::end(rtop), 
			 std::greater()
			 );
 
  maxr1 = rtop.front();

  MPI_Allreduce(&maxr1, &maxr, 1, MPI_DOUBLE, MPI_SUM, MPI_COMM_WORLD);
  maxr /= numprocs;		// Mean percentile target

  dr = (log(max<double>(RHMAX, maxr)) - log(RHMIN))/(NHR-1);
  
  if (myid==0 and VFLAG & 1) {
    std::cout << std::endl
	      << "Table halo: RHMIN=" << RHMIN 
	      << " RHMAX=" << RHMAX
	      << " maxr=" << maxr
	      << " dr=" << dr
	      << std::endl;
  }

  maxr = std::min<double>(maxr, halo->get_max_radius()*0.99);

  double x, y, z, theta, r, fr, fz, fp, pot, costh, sinth;
  double dens, potl, dpr, dpt, dpp, dpdr;
  Eigen::VectorXd work(NHR);
  Eigen::VectorXd workR(NHR);
  Eigen::VectorXd workA(NHR);
  
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
      
      theta = acos(z/(r + std::numeric_limits<double>::min()));
      expandh->determine_fields_at_point(r, theta, 0.0,
					 &dens, &potl, &dpr, &dpt, &dpp);
      
      dpdr = dpr - fz*cos(theta) - fr*sin(theta);
      
      work[j] = halo->get_density(r) * dpdr * r;
      
      if (fabs(work[j]) > 1e20) {
	std::cout << "rho="  << halo->get_density(r)
		  << " dpr=" << dpr
		  << " fz="  << fz
		  << " fr="  << fr
		  << " r="   << r << std::endl;
      }

    }
    
    // Splsum(workR, work, workA);
    Trapsum(workR, work, workA);
    for (int j=0; j<NHR; j++) {
      halotable(i, j) = max<double>(workA[NHR-1] - workA[j], std::numeric_limits<double>::min());
      if (fabs(halotable(i, j))>1.0e8) {
	std::cerr << "Oops, val=" << halotable(i, j)
		  << " costh=" << costh << " r=" << r
		  << std::endl;
      }
    }
  }
  
  // Update tables on all nodes
  //
  Eigen::VectorXd Z(NHR);
  for (int k=0; k<numprocs; k++) {
    for (int i=ibeg[k]; i<iend[k]; i++) {
      if (k == myid) Z = halotable.row(i);
      MPI_Bcast(Z.data(), NHR, MPI_DOUBLE, k, MPI_COMM_WORLD);
      if (k != myid) halotable.row(i) = Z;
    }
  }
  
  //
  // DEBUG output
  //
  if (myid==0 && expandh && VFLAG & 4) {
    ostringstream sout;
    sout << "table_halo." << RUNTAG;
    ofstream out(sout.str().c_str());
    out.setf(ios::scientific);
    out.precision(8);
    
    for (int j=0; j<NHT; j++) {
      for (int k=0; k<NHR; k++) {
	out << std::setw(18) << -1.0 + dc*j
	    << std::setw(18) << RHMIN*exp(dr*k)
	    << std::setw(18) << halotable(j, k) << std::endl;
      }
      out << std::endl;
    }
    
    ostringstream sout2;
    sout2 << "test_halo_pot." << RUNTAG;
    ofstream out2(sout2.str().c_str());
    out2.setf(ios::scientific);
    out2.precision(2);
    
    double rr, t, p, dp;
    for (int k=0; k<NHR; k++) {
      rr = RHMIN*exp(dr*k);
      expandh->determine_fields_at_point(rr, 0.5*M_PI, 0.0,
					 &t, &p, &dp, &t, &t);
      out2
	<< std::setw(14) << rr
	<< std::setw(14) << p
	<< std::setw(14) << dp
	<< std::setw(14) << get_disp(rr, 0.0, 0.0)
	<< std::endl;
    }
    
    
    ostringstream sout3;
    sout3 << "disp_diff." << RUNTAG;
    ofstream out3(sout3.str().c_str());
    double xx, zz, costh;
    
    for (int j=0; j<NHT; j++) {
      costh = -1.0 + dc*j;
      out3 << std::setw(14) << costh;
      for (int k=0; k<NHR; k++) {
	rr = RHMIN*exp(dr*k);
	xx = rr*sqrt(1.0 - costh*costh);
	zz = rr*costh;
	out3 
	  << std::setw(14) << halotable(j, k) -  
	  get_disp(xx, 0.0, zz) * halo->get_density(rr);
      }
      out3 << std::endl;
    }
  }
  
  if (myid==0) std::cout << "[table] " << std::flush;
}

double DiskHalo::get_disp(double xp,double yp, double zp)
{
  if (MULTI) {
    double r = sqrt(xp*xp + yp*yp + zp*zp);
    r = max<double>(r, halo2->get_min_radius());
    r = min<double>(r, halo2->get_max_radius());

    Eigen::VectorXd X(halotable.row(0)), Y(halotable.row(1));
    return odd2(r, X, Y, 0);
  }

  if (halotable.cols() != NHR) {
    std::cerr << "DiskHalo::get_disp: must call table_halo first\n";
    MPI_Abort(MPI_COMM_WORLD, 100);
    exit(0);
  }
  
  int it, ir;
  double r, lr, t, ct[2], cr[2], resv;
  
  // Polar angle
  //
  r = sqrt(xp*xp + yp*yp + zp*zp);
  t = zp/(r + std::numeric_limits<double>::min()) + 1.0;
  
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
    
    ct[0]*cr[0] * halotable(it  , ir  )  +
    ct[0]*cr[1] * halotable(it  , ir+1)  +
    ct[1]*cr[0] * halotable(it+1, ir  )  +
    ct[1]*cr[1] * halotable(it+1, ir+1)  ;
  
  double dens = halo->get_density(r);
  if (dens>0.0) resv /= dens;
  
  return resv;
}

void DiskHalo::set_vel_halo(vector<Particle>& part)
{
  if (!expandh) {
    if (myid==0) std::cout << "[no halo particles] ";
    return;
  }
  
  int nok, ncntE=0, ncntJ=0;
  double v2r, vr, r;
  double vel[3], vel1[3], massp, massp1;
  
  for (int k=0; k<3; k++) vel[k] = vel1[k] = 0.0;
  massp = massp1 = 0.0;
  
  table_halo(part);
  
  for (auto &p : part) {
    
    r = sqrt(p.pos[0]*p.pos[0] + 
	     p.pos[1]*p.pos[1] +
	     p.pos[2]*p.pos[2]);
    
				// Reset success flag
    nok = 1;
    
				// Use Eddington
    
    if (DF && 0.5*(1.0+erf((r-R_DF)/DR_DF)) > rndU(gen)) {
      halo2->gen_velocity(&p.pos[0], &p.vel[0], nok);
      
      if (nok) {
	std::cout << "gen_velocity failed: "
		  << p.pos[0] << " "
		  << p.pos[1] << " "
		  << p.pos[2] << "\n";
      } else ncntE++;
    }
				// Use Jeans
    if (nok) {
      v2r = get_disp(p.pos[0], p.pos[1], p.pos[2]);
      vr = sqrt(max<double>(v2r, std::numeric_limits<double>::min()));
      for (int k=0; k<3; k++) p.vel[k] = vr*rndN(gen);
      ncntJ++;
    }
    
    massp1 += p.mass;
    for (int k=0; k<3; k++) vel1[k] += p.mass*p.vel[k];
  }
  

  MPI_Allreduce(&massp1, &massp, 1, MPI_DOUBLE, MPI_SUM, MPI_COMM_WORLD);
  MPI_Allreduce(vel1, vel, 3, MPI_DOUBLE, MPI_SUM, MPI_COMM_WORLD);

  MPI_Allreduce(MPI_IN_PLACE, &ncntE, 1, MPI_INT, MPI_SUM, MPI_COMM_WORLD);
  MPI_Allreduce(MPI_IN_PLACE, &ncntJ, 1, MPI_INT, MPI_SUM, MPI_COMM_WORLD);
    
  if (massp>0.0) {
    for (int k=0; k<3; k++) vel[k] /= massp;
  }
    
  if (myid==0) {
    std::cout << " (u, v, w)=(" << vel[0] 
	      << ", " << vel[1]
	      << ", " << vel[2] << ")" << std::endl;
  }

  if (cov) {
    for (auto &p : part) {
      for (int k=0; k<3; k++) p.vel[k] -= vel[k];
    }
  }
  
  if (myid==0)
    std::cout << "Velocity calls: Eddington=" << ncntE
	      << " Jeans=" << ncntJ << std::endl;
}

void DiskHalo::write_record(ostream &out, SParticle &p)
{
  bufout << " " << std::setw(16) << setprecision(8) << p.mass;
  
  for (int k=0; k<3; k++)
    bufout << std::setw(24) << setprecision(15) << p.pos[k] + center_pos[k];
  
  for (int k=0; k<3; k++)
    bufout << std::setw(24) << setprecision(15) << p.vel[k] + center_vel[k];
  
  bufout << std::endl;

  if (++bufcnt==bunchcnt) {	// Write the buffer
    out << bufout.str();
    bufout.str("");		// Clear the buffer
    bufcnt = 0;
  }
}

void DiskHalo::flush_buffer(ostream &out)
{
  // Write remaining characters in the buffer and clear
  //
  if (bufout.str().size()>0) out << bufout.str();
  bufout.str("");
  bufcnt = 0;
}

void DiskHalo::write_file(ostream &fou, vector<Particle>& part)
{
  std::vector<SParticle> buf(NBUF);
  
  // Make MPI datatype
  //
  SPtype spt;
  
  // Get particle totals
  //
  int npart=0, l = part.size();

  // Return if there are no particles to write
  //
  MPI_Allreduce(&l, &npart, 1, MPI_INT, MPI_SUM, MPI_COMM_WORLD);

  if (npart==0) return;
  
  // Begin writing loop
  //
  if (myid==0) {
    
    if (VFLAG & 1)
      std::cout << std::endl
		<< "Total number of particles: n=" << npart << std::endl;
    
    fou.setf(ios::scientific);
    
    fou << npart << " " << 0 << " " << 0 << std::endl;
    
    if (VFLAG & 1) {
      std::cout << "Particle stream is ";
      if (fou.good()) std::cout << "GOOD" << std::endl;
      else std::cout << "BAD" << std::endl;
    }

    for (int i=0; i<l; i++)
      write_record(fou, part[i]);
    
    if (VFLAG & 1) {
      std::cout << "Wrote " << l  << " particles from Node 0" << std::endl;
    }

    int imany, icur, ccnt;

    for (int n=1; n<numprocs; n++) {
      
      MPI_Recv(&imany, 1, MPI_INT, n, 10, MPI_COMM_WORLD, MPI_STATUS_IGNORE);
      ccnt=0;
      while (ccnt<imany) {
	MPI_Recv(&icur, 1, MPI_INT, n, 11, MPI_COMM_WORLD, MPI_STATUS_IGNORE);
	MPI_Recv(&buf[0], icur, spt(), n, 12, MPI_COMM_WORLD, MPI_STATUS_IGNORE);
	for (int i=0; i<icur; i++) write_record(fou, buf[i]);
	ccnt += icur;
      }
      
      flush_buffer(fou);

      if (VFLAG & 1)
	std::cout << "Wrote " << ccnt << " particles from Node " << n << std::endl;

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
	  for (int j=0; j<ipack; j++) buf[j][part[icur+j]];
	  MPI_Send(&ipack, 1, MPI_INT, 0, 11, MPI_COMM_WORLD);
	  MPI_Send(&buf[0], ipack, spt(), 0, 12, MPI_COMM_WORLD);
	  icur += ipack;
	}

	if (VFLAG & 1)
	  std::cout << "Sent " << icur << " particles from Node " << n << std::endl;
      }

      MPI_Barrier(MPI_COMM_WORLD);
    }
  }
  
}


void DiskHalo::virial_ratio(vector<Particle>& hpart, vector<Particle>& dpart)
{
  double r, theta, phi, xx, yy, zz, axd, ayd, azd, axh, ayh, azh, R2, R;
  double dens, potl, potr, pott, potp, fr, fp, fz;
  
  double KE_disk1      = 0.0;
  double KE_halo1      = 0.0;
  double PE_disk_disk1 = 0.0;
  double PE_halo_disk1 = 0.0;
  double PE_disk_halo1 = 0.0;
  double PE_halo_halo1 = 0.0;
  double massd1        = 0.0;
  double massh1        = 0.0;
  
  fr   = fp   = fz   = 0.0;
  potr = pott = potp = 0.0;
				// -----------------
				// Halo contribution
				// -----------------
  for (auto p : hpart) {
    r = 0.0;
    for (int k=0; k<3; k++) {
      r += p.pos[k]*p.pos[k];
      KE_halo1 += 0.5*p.mass*p.vel[k]*p.vel[k];
    }
    
    r = sqrt(r);
    xx = p.pos[0];
    yy = p.pos[1];
    zz = p.pos[2];
    
    theta = acos(zz/(r+std::numeric_limits<double>::min()));
    phi = atan2(yy, xx);
    
    if (expandh)
      expandh->determine_fields_at_point(r, theta, phi,
					 &dens, &potl, &potr, &pott, &potp);
    
    R2 = xx*xx + yy*yy + std::numeric_limits<double>::min();
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
    PE_halo_disk1 += p.mass * (xx*axd + yy*ayd + zz*azd);
    PE_halo_halo1 += p.mass * (xx*axh + yy*ayh + zz*azh);
    
				// Mass
    massh1 += p.mass;
  }
  
				// -----------------
				// Disk contribution
				// -----------------
  for (auto p : dpart) {
    r = 0.0;
    for (int k=0; k<3; k++) {
      r += p.pos[k]*p.pos[k];
      KE_disk1 += 0.5*p.mass*p.vel[k]*p.vel[k];
    }
    
    r = sqrt(r);
    xx = p.pos[0];
    yy = p.pos[1];
    zz = p.pos[2];
    
    theta = acos(zz/(r+std::numeric_limits<double>::min()));
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
    PE_disk_disk1 += p.mass * (xx*axd + yy*ayd + zz*azd);
    PE_disk_halo1 += p.mass * (xx*axh + yy*ayh + zz*azh);
    
				// Mass
    massd1 += p.mass;
    
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
    
    std::cout << std::endl
	      << "****************************" << std::endl
	      <<"  KE_halo  = " << KE_halo << std::endl
	      <<"  KE_disk  = " << KE_disk << std::endl << std::endl
	      <<"  PE_halo(halo)  = " << PE_halo_halo << std::endl
	      <<"  PE_halo(disk)  = " << PE_halo_disk << std::endl
	      <<"  PE_disk(disk)  = " << PE_disk_disk << std::endl
	      <<"  PE_disk(halo)  = " << PE_disk_halo << std::endl << std::endl;
    if (PE_halo < 0.0)
      std::cout <<"-2T/W_halo = " << -2.0*KE_halo/PE_halo << std::endl;
    if (PE_disk < 0.0)
      std::cout <<"-2T/W_disk = " << -2.0*KE_disk/PE_disk << std::endl;
    std::cout << std::endl;
    
    std::cout << " Halo mass=" << massh << "  Disk mass=" << massd << std::endl << std::endl;
    
    double KE = KE_halo + KE_disk;
    double PE = PE_halo + PE_disk;
    
    std::cout << "  KE       = " << KE << std::endl
	      << "  PE       = " << PE << std::endl;
    if (PE<0.0)
      std::cout << " -2T/W     = " << -2.0*KE/PE << std::endl;
    std::cout << "****************************" << std::endl;
    
  } 
}


void DiskHalo::virial_ratio(const char *hfile, const char *dfile)
{
  ifstream in[2];
  
  in[0].open(hfile);
  in[1].open(dfile);

  if (!in[0]) {
    std::cout << "virial_ratio: can't open " << hfile << std::endl;
    return;
  }
  
  if (!in[1]) {
    std::cout << "virial_ratio: can't open " << dfile << std::endl;
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
      if (!in[c]) break;
    
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

      theta = acos(zz/(r+std::numeric_limits<double>::min()));
      phi = atan2(yy, xx);
    
      if (expandh)
	expandh->determine_fields_at_point(r, theta, phi,
					 &dens, &potl, &potr, &pott, &potp);
      
      R2 = xx*xx + yy*yy + std::numeric_limits<double>::min();
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

    std::cout << std::endl
	      << "****************************" << std::endl
	      <<"  Total #  = " << count0 << std::endl
	      <<"  Expect # = " << expected << std::endl << std::endl
	      <<"  KE_halo  = " << KE_halo0 << std::endl
	      <<"  KE_disk  = " << KE_disk0 << std::endl << std::endl
	      <<"  PE_halo  = " << PE_halo0 << std::endl
	      <<"  PE_disk  = " << PE_disk0 << std::endl << std::endl;

    if (PE_halo0 < 0.0)
      std::cout <<"-2T/W_halo = " << -2.0*KE_halo0/PE_halo0 << std::endl;
    if (PE_disk < 0.0)
      std::cout <<"-2T/W_disk = " << -2.0*KE_disk0/PE_disk0 << std::endl;
    std::cout << std::endl;

    std::cout << " Halo mass=" << massh0 << "  Disk mass=" << massd0 
	 << std::endl << std::endl;

    double KE = KE_halo0 + KE_disk0;
    double PE = PE_halo0 + PE_disk0;

    std::cout << "  KE       = " << KE << std::endl
	 << "  PE       = " << PE << std::endl;
    if (PE<0.0)
      std::cout << " -2T/W     = " << -2.0*KE/PE << std::endl;
    std::cout << "****************************" << std::endl;
    
  }

}



void DiskHalo::profile(ostream &out, vector<Particle>& dpart,
		       double rmin, double rmax, int numr)
{
  if (dpart.size() == 0) return;

  bool logr = false;
  if (rmin>0.0) {
    logr = true;
    rmin = log(rmin);
    rmax = log(rmax);
  }

  double dR = (rmax - rmin)/numr;
  int indx;

  vector<double> mass1(numr, 0.0);
  vector<double> velR1(numr, 0.0);
  vector<double> velT1(numr, 0.0);
  vector<double> sigR1(numr, 0.0);
  vector<double> sigT1(numr, 0.0);

  for (auto p : dpart) {
    double xx = p.pos[0];
    double yy = p.pos[1];
    double R  = sqrt(xx*xx + yy*yy);

    if (logr) 
      indx = (log(R) - rmin)/dR;
    else
      indx = R/dR;

    if (indx < numr && indx >= 0) {

      double vr = ( xx*p.vel[0] + yy*p.vel[1])/(R+std::numeric_limits<double>::min());
      double vt = (-yy*p.vel[0] + xx*p.vel[1])/(R+std::numeric_limits<double>::min());

      mass1[indx] += p.mass;
      velR1[indx] += p.mass*vr;
      velT1[indx] += p.mass*vt;
      sigR1[indx] += p.mass*vr*vr;
      sigT1[indx] += p.mass*vt*vt;
    }
  }
  
  
  vector<double> mass (numr);
  vector<double> velR (numr);
  vector<double> velT (numr);
  vector<double> sigR (numr);
  vector<double> sigT (numr);


  MPI_Reduce(&mass1[0], &mass[0], numr, MPI_DOUBLE, MPI_SUM, 0, MPI_COMM_WORLD);
  MPI_Reduce(&velR1[0], &velR[0], numr, MPI_DOUBLE, MPI_SUM, 0, MPI_COMM_WORLD);
  MPI_Reduce(&velT1[0], &velT[0], numr, MPI_DOUBLE, MPI_SUM, 0, MPI_COMM_WORLD);
  MPI_Reduce(&sigR1[0], &sigR[0], numr, MPI_DOUBLE, MPI_SUM, 0, MPI_COMM_WORLD);
  MPI_Reduce(&sigT1[0], &sigT[0], numr, MPI_DOUBLE, MPI_SUM, 0, MPI_COMM_WORLD);

  
  if (myid==0) {
    out << "#"      << right 
	<< std::setw(14) << "Radius"	 // #1
	<< std::setw(15) << "Mass"	 // #2
	<< std::setw(15) << "S(mass)" // #3
	<< std::setw(15) << "Density" // #4
	<< std::setw(15) << "V_c"	 // #5
	<< std::setw(15) << "kappa"	 // #6
	<< std::setw(15) << "Omega"	 // #7
	<< std::setw(15) << "V_R"	 // #8
	<< std::setw(15) << "V_T"	 // #9
	<< std::setw(15) << "Sig_R"	 // #10
	<< std::setw(15) << "Sig_T"	 // #11
	<< std::endl;

    double smass = 0.0, rin, rout, ravg;
    for (int i=0; i<numr; i++) {
      if (logr) {
	rin  = exp(rmin + dR*i);
	rout = exp(rmin + dR*(i+1));
      } else {
	rin  = dR*i;
	rout = dR*(i+1);
      }
      ravg = 0.5*(rin + rout);

      out << std::setw(15) << ravg
	  << std::setw(15) << mass[i]
	  << std::setw(15) << (smass += mass[i])
	  << std::setw(15) << mass[i]/(M_PI*(rout*rout - rin*rin))
	  << std::setw(15) << v_circ(ravg, 0, 0)
	  << std::setw(15) << epi(ravg, 0, 0)
	  << std::setw(15) << v_circ(ravg, 0, 0)/ravg;

      if (mass[i]>0.0) {
	double vr = velR[i]/mass[i];
	double vt = velT[i]/mass[i];
	out << std::setw(15) << vr
	    << std::setw(15) << vt
	    << std::setw(15) << sqrt(fabs(sigR[i]/mass[i] - vr*vr))
	    << std::setw(15) << sqrt(fabs(sigT[i]/mass[i] - vt*vt))
	    << std::endl;
      } else {
	out << std::setw(15) << 0.0
	    << std::setw(15) << 0.0
	    << std::setw(15) << 0.0
	    << std::setw(15) << 0.0
	    << std::endl;
      }
    }
  } 
}

void DiskHalo::disk_model(const std::string &modfile)
{
  std::vector<double> r2(RNUM);
  std::vector<double> rr(RNUM);
  std::vector<double> d2(RNUM);
  std::vector<double> m2(RNUM);
  std::vector<double> p2(RNUM);
  std::vector<double> pp(RNUM);
  std::vector<double> dd(RNUM);
  std::vector<double> dm(RNUM);

  double rmin = halo->get_min_radius();
  double rmax = halo->get_max_radius();
  double r, dr;

  if (LOGR)
    dr = (log(rmax) - log(rmin))/(RNUM-1);
  else
    dr = (rmax - rmin)/(RNUM-1);
  
  for (int i=0; i<RNUM; i++) {
    if (LOGR) {
      r2[i] = r = rmin*exp(dr*i);
      rr[i] = log(rmin) + dr*i;
    } else {
      r2[i] = r = rmin + dr*i;
      rr[i] = r;
    }

    // dPhi/dR
    dd[i] = deri_pot(r, 0.0, 0.0, 1);

    // M(R) = R^2*dPhi/dR
    m2[i] = dd[i]*r*r;
  }
  
  auto spl = Spline1d(rr, m2);
  for (int i=0; i<RNUM; i++) {
    // dM(R)/dR or dM(R)/dlnR
    dm[i] = spl.deriv(rr[i]);
    if (LOGR) dm[i] /= r2[i];

    // rho(R) = dM(R)/dR/(4*pi*R^2)
    d2[i] = dm[i]/(4.0*M_PI*r2[i]*r2[i]);
  }

  // Monopole potential
  pp[0] = 0.0;
  for (int i=1; i<RNUM; i++) {
    pp[i] = pp[i-1] + 0.5*(dm[i]/r2[i] + dm[i-1]/r2[i-1])*(r2[i] - r2[i-1]);
  }

  for (int i=0; i<RNUM; i++) {
    p2[i] = -m2[i]/r2[i] + pp[i] - pp.back();
  }

  // Debug
  if (true) {
    for (int i=0; i<RNUM; i++) {
      std::cout << "** "
		<< std::setw(16) << r2[i]
		<< std::setw(16) << deri_pot_disk(r2[i], 0, 0, 1)
		<< std::setw(16) << deri_pot_halo(r2[i], 0, 0, 1)
		<< std::setw(16) << halo->get_density(r2[i])
		<< std::setw(16) << d2[i]
		<< std::endl;
    }
  }

  SphericalModelTable(r2, d2, m2, p2).print_model(modfile);
}
