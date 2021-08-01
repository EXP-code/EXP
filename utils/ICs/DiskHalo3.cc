				// C++/STL
#include <iostream>
#include <fstream>
#include <sstream>
#include <algorithm>
#include <vector>
#include <limits>
				// EXP classes
#include <interp.H>
#include <numerical.H>
#include <exponential.H>
#include <interp.H>
				// Local
#include <AddDisk.H>
#include <CylDisk.H>
#include "DiskWithHalo.H"
#include <DiskHalo3.H>

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

int    DiskHalo::NCHEB    = 16;
bool   DiskHalo::CHEBY    = false;

unsigned DiskHalo::VFLAG = 0;
unsigned DiskHalo::NBUF = 8192;

string DiskHalo::RUNTAG = "debug";

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
 }

DiskHalo::
DiskHalo(std::shared_ptr<SphericalSL> haloexp, std::shared_ptr<EmpCylSL> diskexp,
	 double H, double A, double DMass, 
	 string& filename, int DF1, int DIVERGE, double DIVERGE_RFAC)
{
  gen.seed(SEED+myid);

  com        = false;
  cov        = false;

  for (int k=0; k<3; k++) center_pos[k] = center_vel[k] = 0.0;

  if (DF1) DF = true;
  else     DF = false;

  MULTI = false;

  dmass = DMass;
  scalelength = A;
  scaleheight = H;

  expandh = haloexp;
  expandd = diskexp;

  SphericalModelTable::even = 0;
  SphericalModelTable::logscale = LOGSCALE;

  halo = std::make_shared<SphericalModelTable>(filename, DIVERGE, DIVERGE_RFAC);

  disk = std::make_shared<ExponentialDisk>(A, RDMAX, DMass);

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
    newmod = std::make_shared<AddDisk>(halo, disk, COMPRESSION); 
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
DiskHalo(std::shared_ptr<SphericalSL> haloexp, std::shared_ptr<EmpCylSL> diskexp,
	 double H, double A, double DMass, 
	 std::string& filename1, int DIVERGE, double DIVERGE_RFAC,
	 std::string& filename2, int DIVERGE2, double DIVERGE_RFAC2)
{
  gen.seed(SEED+myid);

  com         = false;
  cov         = false;

  for (int k=0; k<3; k++) center_pos[k] = center_vel[k] = 0.0;

  DF          = false;
  MULTI       = true;

  dmass       = DMass;
  scaleheight = H;

  expandh     = haloexp;
  expandd     = diskexp;

  SphericalModelTable::even = 0;
  SphericalModelTable::logscale = LOGSCALE;

  halo = std::make_shared<SphericalModelTable>(filename1, DIVERGE, DIVERGE_RFAC);

  disk = std::make_shared<ExponentialDisk>(A, RDMAX, DMass);

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
  newmod = std::make_shared<AddDisk>(halo, disk, COMPRESSION); 
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
    halo3->print_model("diskhalo2_model.multi");
    halo3->print_df("diskhalo2_df.multi");
  }
    
  //
  // Generate the multimass model
  //
  multi = std::make_shared<SphericalModelMulti>(halo2.get(), halo3.get());
  multi -> gen_tolE = TOLE;

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

  dP = p.dP;
  dR = p.dR;
  dZ = p.dZ;

  halotable = p.halotable;
  dr = p.dr;
  dc = p.dc;

  DF    = p.DF;
  MULTI = p.MULTI;
  com   = p.com;
  cov   = p.cov;
}


double DiskHalo::disk_density(double R, double z)
{
  double q = 1.0/cosh(z/scaleheight);
  return disk->get_density(R)*q*q*0.5/scaleheight;
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
  Eigen::VectorXd ps(7);
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
    for (auto &p : phalo) {
      if (com) for (int k=0; k<3; k++) p.pos[k] -= pos[k];
      if (cov) for (int k=0; k<3; k++) p.vel[k] -= vel[k];
    }
  }

  if (VFLAG & 1)
    cout << "Process " << myid << ": made " << phalo.size() << " particles"
	 << endl;
}      

void DiskHalo::
set_halo_coordinates(vector<Particle>& phalo, int nhalo, int npart)
{
  const double tol = 1.0e-8;
  double rmin = max<double>(halo->get_min_radius(), RHMIN);
  double rmax = halo->get_max_radius();
  double mmin = halo->get_mass(rmin);
  double mtot = halo->get_mass(rmax);

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
    for (auto &p : phalo) {
      for (int k=0; k<3; k++) p.pos[k] -= pos[k];
    }
  }

  if (VFLAG & 1)
    cout << "Process " << myid << ": made " << phalo.size() << " particles"
	 << endl;
}      


double DiskHalo::get_hpot(double xp, double yp, double zp)
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


void DiskHalo::make_disk_DF(bool diag)
{
  
  int QPverbose    = 0;		// verbosity
  int QPnumE       = 20;	// energy grid number
  int QPnumK       = 10;	// kappa  grid number
  int QPnumR       = 40;	// radial grid number
  bool QPmassE     = true;	// energy gridding using mass
  bool QPmassL     = false;     // energy gridding in mass using 
				// linear(true) or log(false) scaling

  double QPsigma   = 1.25;      // kernel width fraction of grid spacing
  double QPlambda  = 3.0;	// kappa power-law bias
  double QPalpha   = -4.0;	// kappa power-law bias exponent
  double QPbeta    = 2.0;	// mass power-lass exponent
  double QPkmin    = 0.35;	// minimum value of kappa
  double QPkmax    = 1.0;	// maximum value of kappa

  bool RLOG        = true;	// diag output logarithmic in radius
  int NUMS         = 100;       // number of surface grid points
  int NUMR         = 100;	// number of radial points for diag
  int NINT         = 200;       // number of integration knots for density computation from DF
  int NUMMS        = 200;       // number of integration knots for mass model computation from density
  int NUME         = 800;       // number of energy grid points

  string DTAG("qpdist");

  double rmin0 = 1.05*disk->get_min_radius();
  double rmax0 = 0.95*disk->get_max_radius();

  //
  // Compute QPDISTF
  //

  dmodel = DiskWithHalo(disk.get(), halo.get());

  qp = std::make_shared<QPDistF>(disk.get(), halo.get(), rmax0, rmax0, 
				 QPnumE, QPnumK, QPnumR, QPsigma, QPlambda, QPalpha, QPbeta,
				 1.0, 0.01, 0.05, 0.5,
				 QPkmin, QPkmax);

  if (diag && myid==0) qp->set_verbose();

  SphericalOrbit::ZFRAC=0.5;
  qp->MassEGrid  = QPmassE;
  qp->MassLinear = QPmassL; 
  qp->compute_distribution();
  qp->make_cdf(100, 100);

  if (myid || !diag) return;

  qp->dump_cdf(string("qpdist.cdf"));

  //
  // Begin integration check
  //

  vector<ofstream*> out(4);
  out[0] = new ofstream(string(DTAG + "_1.dat").c_str());
  out[1] = new ofstream(string(DTAG + "_2.dat").c_str());
  out[2] = new ofstream(string(DTAG + ".df"  ).c_str());
  out[3] = new ofstream(string(DTAG + ".surf").c_str());

  if (*out[0]) {

    *out[0]  << "# "
	     << setw(13) << "Radius"
	     << setw(15) << "Dens (orig)"
	     << setw(15) << "Dens (DF)"
	     << setw(15) << "Difference"
	     << setw(15) << "Vel tang"
	     << setw(15) << "Vel disp"
	     << endl;
    
    *out[0]  << "# "
	     << setw(13) << "------"
	     << setw(15) << "-----------"
	     << setw(15) << "---------"
	     << setw(15) << "----------"
	     << setw(15) << "----------"
	     << setw(15) << "--------"
	     << endl;
  }


  if (*out[1]) {

    *out[1]  << "# "
	     << setw(13) << "Radius"
	     << setw(15) << "Dens (DF)"
	     << setw(15) << "Potential"
	     << setw(15) << "Vel tang"
	     << setw(15) << "Vel disp"
	     << endl;
    
    *out[1] << "# "
	    << setw(13) << "------"
	    << setw(15) << "---------"
	    << setw(15) << "---------"
	    << setw(15) << "---------"
	    << setw(15) << "--------"
	    << endl;
  }
    
  cout  << setw(15) << "Radius"
	<< setw(15) << "Dens (orig)"
	<< setw(15) << "Dens (DF)"
	<< setw(15) << "Difference"
	<< setw(15) << "Vel mean"
	<< setw(15) << "Vel disp"
	<< endl;
    
  cout  << setw(15) << "------"
	<< setw(15) << "-----------"
	<< setw(15) << "---------"
	<< setw(15) << "----------"
	<< setw(15) << "----------"
	<< setw(15) << "--------"
	<< endl;


  LegeQuad lq(NINT);

  double den, vmax, r, pot, pmax, E, J, x, y, vmean, vdisp;
  double dr, rmin, rmax;

  vector<double> rv, dv, mv, pw, p0;

  rv.reserve(NUMR);
  dv.reserve(NUMR);

  if (RLOG) {
    rmin = log(dmodel.get_min_radius());
    rmax = log(dmodel.get_max_radius());
  }
  else {
    rmin = dmodel.get_min_radius();
    rmax = dmodel.get_max_radius();
  }
  dr = (rmax - rmin)/NUMR;

  pmax = dmodel.get_pot(dmodel.get_max_radius());

  for (int i=1; i<=NUMR; i++) {
    if (RLOG)
      rv[i-1] = r = exp(rmin + dr*((double)i-0.5));
    else
      rv[i-1] = r = rmin + dr*((double)i-0.5);
    
    pot = dmodel.get_pot(r);

    vmax = sqrt(2.0*fabs(pmax - pot));
    den = vmean = vdisp = 0.0;
    for (int ix=0; ix<NINT; ix++) {
      x = lq.knot(ix);

      for (int iy=0; iy<NINT; iy++) {
	y = lq.knot(iy);

	E = pot + 0.5*vmax*vmax*(x*x + (1.0-x*x)*y*y);
	J = vmax*sqrt(1.0 - x*x)*y*r;
	den += lq.weight(ix)*lq.weight(iy) * vmax*vmax * sqrt(1.0 - x*x) *
	  qp->distf(E, J);
	vmean += lq.weight(ix)*lq.weight(iy) * vmax*vmax * sqrt(1.0 - x*x) *
	  vmax*sqrt(1.0-x*x)*y * qp->distf(E, J);
	vdisp += lq.weight(ix)*lq.weight(iy) * vmax*vmax * sqrt(1.0 - x*x) *
	  vmax*vmax*(x*x + (1.0-x*x)*y*y) * qp->distf(E, J);
      }
    }
    
    if (den>0.0) {
      vmean = vmean/den;
      vdisp = vdisp/den - vmean*vmean;
    } else {
      vmean  = 0.0;
      vdisp  = 0.0;
    }

    den *= 4.0;

    dv[i-1] = den;

    if (*out[0])
      *out[0]  << setw(15) << r 
	       << setw(15) << dmodel.get_density(r) 
	       << setw(15) << den 
	       << setw(15) << den - dmodel.get_density(r) 
	       << setw(15) << vmean
	       << setw(15) << vdisp
	       << endl;
    
    pot = dmodel.get_pot(r);

    if (*out[1])
      *out[1]  << setw(15) << r 
	       << setw(15) << den 
	       << setw(15) << pot
	       << setw(15) << vmean
	       << setw(15) << vdisp
	       << endl;


    cout << setw(15) << r 
	 << setw(15) << dmodel.get_density(r) 
	 << setw(15) << den 
	 << setw(15) << den - dmodel.get_density(r) 
	 << setw(15) << vmean
	 << setw(15) << vdisp
	 << endl;
  }
    
  LegeQuad lq2(NUMMS);

  double energy=0.0, mass=0.0;
  dr = rmax - rmin;

  for (int i=0; i<NUMMS; i++) {

    if (RLOG)
      r = exp(rmin + dr*lq2.knot(i));
    else
      r = rmin + dr*lq2.knot(i);

    pot = dmodel.get_pot(r);
    den = dmodel.get_density(r);

    if (RLOG) {
      mass += r*r*den * lq2.weight(i);
      energy += 0.5*r*r*den*pot * lq2.weight(i);
    } else {
      mass += r*den * lq2.weight(i);
      energy += 0.5*r*den*pot * lq2.weight(i);
    }
  }

  mass *= 2.0*M_PI*dr;
  energy *= 2.0*M_PI*dr;

  cout << "Mass=" << mass << endl;
  cout << "Energy=" << energy << endl;

  double Emin = dmodel.get_pot(rmin0);
  double Emax = dmodel.get_pot(rmax0);

  double dE = (Emax - Emin)/(NUME-1);
  if (*out[2]) {
    *out[2] << setw(15) << "# Energy"
	    << setw(15) << "F(E)"
	    << setw(15) << "dF(E)/dE"
	    << endl;
    for (int i=0; i<NUME; i++) 
      *out[2] << setw(15) << Emin + dE*i
	      << setw(15) << qp->distf_EK(Emin + dE*i, 1.0/sqrt(2.0))
	      << setw(15) << qp->dfdE_EK(Emin + dE*i, 1.0/sqrt(2.0))
	      << endl;
  }
  
  if (*out[3]) {
    out[3]->precision(6);
    out[3]->setf(ios::scientific, ios::floatfield);
  
    double Emin = dmodel.get_pot(rmin0);
    double Emax = dmodel.get_pot(rmax0);

    double dE = (Emax - Emin)/NUMS;
    double dK = 1.0/NUMS;
    for (int i=0; i<NUMS; i++) {
      double E = Emin + dE*(0.5+i);
      for (int j=0; j<NUMS; j++) {
	double K = dK*(0.5+j);
	*out[3] << setw(15) << E
		<< setw(15) << K
		<< setw(15) << qp->distf_EK(E, K)
		<< setw(15) << qp->dfdE_EK(E, K)
		<< endl;
      }
      *out[3] << endl;
    }
  }

  for (int i=0; i<4; i++) delete out[i];

}


void DiskHalo::
set_disk(vector<Particle>& pdisk, int ndisk, int npart)
{
  const double tol = 1.0e-8;
  double rmin = max<double>(disk->get_min_radius(), RDMIN);
  double rmax = disk->get_max_radius();
  double mmin = disk->get_mass(rmin);
  double mtot = disk->get_mass(rmax);

  double R, phi, E, K, f, vt, vr, vz, w1, T;
  double pos[3], vel[3], pos1[3], vel1[3], massp, massp1;

  pair<double, double> pr;

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
  for (int k=0; k<3; k++) vel[k] = vel1[k] = 0.0;
  massp = massp1 = 0.0;

  Particle p;
  SphericalOrbit orb(&dmodel);

  p.mass = dmass/ndisk;

  model = disk;

  for (int i=0; i<npart; i++) {
				// Get an E and K from the QP solution
    pr = qp->gen_EK(rndU(gen), rndU(gen));

    E = pr.first;
    K = pr.second;
    orb.new_orbit(E, K);

    T   = 2.0*M_PI*rndU(gen)/orb.get_freq(1);
    w1  = orb.get_angle(1, T);
    f   = orb.get_angle(5, T);
    R   = orb.get_angle(6, T);
    phi = 2.0*M_PI*rndU(gen) + f;

    vt  = orb.Jmax()*K/R;
    vr  = 2.0*(E - dmodel.get_pot(R)) - vt*vt;

    if (vr>0.0) vr = (w1 < M_PI) ? sqrt(vr) : -sqrt(vr);
    else vr = 0.0;

    vz = rndN(gen)*sqrt(M_PI*dmodel.get_density(R)*scaleheight);

    p.pos[0] = R*cos(phi);
    p.pos[1] = R*sin(phi);
    p.pos[2] = scaleheight*atanh(2.0*rndU(gen)-1.0);

    p.vel[0] = vr*cos(phi) - vt*sin(phi);
    p.vel[1] = vr*sin(phi) + vt*cos(phi);
    p.vel[2] = vz;

    pdisk.push_back(p);

    massp1 += p.mass;
    for (int k=0; k<3; k++) pos1[k] += p.mass*p.pos[k];
    for (int k=0; k<3; k++) vel1[k] += p.mass*p.vel[k];

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
  MPI_Allreduce(vel1, vel, 3, MPI_DOUBLE, MPI_SUM, MPI_COMM_WORLD);

  if (massp>0.0) {
    for (int k=0; k<3; k++) {
      pos[k] /= massp;
      vel[k] /= massp;
    }
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
  
  if (com) {
    for (auto &p : pdisk) {
      for (int k=0; k<3; k++) p.pos[k] -= pos[k];
    }
  }

  if (cov) {
    for (auto &p : pdisk) {
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
      out << setw(18) << halotable(0, k)
	  << setw(18) << halotable(1, k)
	  << setw(18) << 4.0*M_PI*halotable(2, k)
	  << setw(18) << halo->get_density(halotable(0, k)) << endl;
  }

}

void DiskHalo::
table_halo(vector<Particle>& part)
{
  if (halotable.cols() == NHR) return;
  
  halotable.resize(NHT, NHR);
  
  dc = 2.0/(NHT-1);
  double r2, maxr = 0.0, maxr1 = 0.0;
  for (auto &p : part) {
    r2 = 0.0;
    for (int k=0; k<3; k++) r2 += p.pos[k]*p.pos[k];
    maxr1 = max<double>(maxr1, sqrt(r2));
  }
  
  MPI_Allreduce(&maxr1, &maxr, 1, MPI_DOUBLE, MPI_MAX, MPI_COMM_WORLD);

  dr = (log(max<double>(RHMAX, maxr)) - log(RHMIN))/(NHR-1);
  
  if (myid==0) {
    cout << endl
	 << "Table halo: RHMIN=" << RHMIN 
	 << " RHMAX=" << RHMAX
	 << " maxr=" << maxr
	 << " dr=" << dr
	 << endl;
  }

  double x, y, z, theta, r, fr, fz, fp, pot, costh, sinth;
  double dens, potl, dpr, dpt, dpp, dpdr;
  Eigen::VectorXd work (NHR);
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
      
    }
    
    // Splsum(workR, work, workA);
    Trapsum(workR, work, workA);
    for (int j=0; j<NHR; j++) {
      halotable(i, j) = max<double>(workA[NHR-1] - workA[j], std::numeric_limits<double>::min());
      if (fabs(halotable(i, j))>1.0e8) {
	cerr << "Oops, val=" << halotable(i, j) << endl;
      }
    }
  }
  
  // Update tables on all nodes
  //
  for (int k=0; k<numprocs; k++) {
    for (int i=ibeg[k]; i<iend[k]; i++) {
      MPI_Bcast(&halotable(i, 0), NHR, MPI_DOUBLE, k, MPI_COMM_WORLD);
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
	out << setw(18) << -1.0 + dc*j
	    << setw(18) << RHMIN*exp(dr*k)
	    << setw(18) << halotable(j, k) << endl;
      }
      out << endl;
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
	<< setw(14) << rr
	<< setw(14) << p
	<< setw(14) << dp
	<< setw(14) << get_disp(rr, 0.0, 0.0)
	<< endl;
    }
    
    
    ostringstream sout3;
    sout3 << "disp_diff." << RUNTAG;
    ofstream out3(sout3.str().c_str());
    double xx, zz, costh;
    
    for (int j=0; j<NHT; j++) {
      costh = -1.0 + dc*j;
      out3 << setw(14) << costh;
      for (int k=0; k<NHR; k++) {
	rr = RHMIN*exp(dr*k);
	xx = rr*sqrt(1.0 - costh*costh);
	zz = rr*costh;
	out3 
	  << setw(14) << halotable(j, k) -  
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
    Eigen::VectorXd X(halotable.row(0));
    Eigen::VectorXd Y(halotable.row(1));
    return odd2(r, X, Y, 0);
  }

  if (halotable.cols() != NHR) {
    cerr << "DiskHalo::get_disp: must call table_halo first\n";
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
    if (myid==0) cout << "[no halo particles] ";
    return;
  }
  
  int nok;
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
	cout << "gen_velocity failed: "
	     << p.pos[0] << " "
	     << p.pos[1] << " "
	     << p.pos[2] << "\n";
      }
    }
				// Use Jeans
    if (nok) {
      v2r = get_disp(p.pos[0], p.pos[1], p.pos[2]);
      vr = sqrt(max<double>(v2r, std::numeric_limits<double>::min()));
      for (int k=0; k<3; k++) p.vel[k] = vr*rndN(gen);
    }
    
    massp1 += p.mass;
    for (int k=0; k<3; k++) vel1[k] += p.mass*p.vel[k];
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
    for (auto &p : part) {
      for (int k=0; k<3; k++) p.vel[k] -= vel[k];
    }
  }
  
}

void DiskHalo::write_record(ostream &out, SParticle &p)
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
  SPtype spt;
  
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
	MPI_Recv(&buf[0], icur, spt(), n, 12, MPI_COMM_WORLD, MPI_STATUS_IGNORE);
	for (int i=0; i<icur; i++) write_record(fou_halo, buf[i]);
	ccnt += icur;
      }
      
      if (VFLAG & 1)
	cout << "Wrote " << ccnt << " HALO particles from Node " << n << endl;

      MPI_Recv(&imany, 1, MPI_INT, n, 13, MPI_COMM_WORLD, MPI_STATUS_IGNORE);
      ccnt = 0;
      while (ccnt<imany) {
	MPI_Recv(&icur, 1, MPI_INT, n, 14, MPI_COMM_WORLD, MPI_STATUS_IGNORE);
	MPI_Recv(&buf[0], icur, spt(), n, 15, MPI_COMM_WORLD, MPI_STATUS_IGNORE);
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
	  MPI_Send(&buf[0], ipack, spt(), 0, 12, MPI_COMM_WORLD);
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
	  MPI_Send(&buf[0], ipack, spt(), 0, 15, MPI_COMM_WORLD);
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
  
  fr = fp = fz = 0.0;
  potr = pott = potp = 0.0;
				// -----------------
				// Halo contribution
				// -----------------
  for (auto &p : hpart) {
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
  for (auto &p : dpart) {
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

