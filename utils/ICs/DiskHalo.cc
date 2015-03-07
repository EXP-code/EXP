				// System
#include <values.h>
				// C++/STL
#include <iostream>
#include <fstream>
#include <sstream>
#include <vector>
				// MDW
#include <interp.h>
#include <numerical.h>
#include <exponential.h>

				// Local
#include "AddDisk.h"
#include "DiskHalo.h"

				// Grid parameters and Toomre Q
double DiskHalo::RHMIN = 1.0e-4;
double DiskHalo::RHMAX = 50.0;
double DiskHalo::RDMIN = 1.0e-4;
double DiskHalo::RDMAX = 20.0;
double DiskHalo::Q = 1.2;
double DiskHalo::SHFACTOR = 16.0;
double DiskHalo::DMFACTOR = 1.0;
unsigned DiskHalo::NBUF = 8192;
int DiskHalo::NDP = 16;
int DiskHalo::NDZ = 40;
int DiskHalo::NDZF = 1;
int DiskHalo::NDR = 800;
int DiskHalo::NHR = 800;
int DiskHalo::NHT = 40;
int DiskHalo::SEED = 11;

double DiskHalo::R_DF = 20.0;
double DiskHalo::DR_DF = 5.0;

int DiskHalo::LOGSCALE = 0;

static AxiSymModel *model;
double targetmass;
				// Determine radius with given enclosed mass
double mass_func(double r)
{
  return targetmass - model->get_mass(r);
}

DiskHalo::
DiskHalo(SphericalSL* haloexp, EmpCylSL* diskexp,
	 double H, double A, double HMass, double DMass, 
	 string& filename, int DF1, int DIVERGE, double DIVERGE_RFAC)
{
  disktableP = NULL;
  disktableN = NULL;
  lwz = NULL;
  gen = new ACG(SEED+myid, 20);
  rndU = new Uniform(0.0, 1.0, gen);
  rndN = new Normal(0.0, 1.0, gen);
  com_cov = false;

  for (int k=0; k<3; k++) center_pos[k] = center_vel[k] = 0.0;

  if (DF1) DF = true;

  hmass = HMass;
  dmass = DMass;
  scaleheight = H;

  expandh = haloexp;
  expandd = diskexp;

  SphericalModelTable::even = 0;
  SphericalModelTable::logscale = LOGSCALE;

  halo = new SphericalModelTable(filename, DIVERGE, DIVERGE_RFAC);

  disk = new ExponentialDisk(A, RDMAX);

  if (myid==0) {
    cerr << "DEBUG: DIVERGE=" << DIVERGE
	 << " A=" << A
	 << " RDMAX=" << RDMAX
	 << " filename=" << filename
	 << "\n";
  }

  if (DF) {
    AddDisk::logarithmic = true;
    AxiSymModel::numr = 400;
    AxiSymModel::numj = 400;
    AxiSymModel::gen_N = 800;
    AxiSymModel::gen_itmax = 40000;
    AxiSymModel::gen_rmin = RHMIN;
    newmod = new AddDisk(halo, disk, DMFACTOR*dmass); 
    halo2 = newmod->get_model();
    halo2->setup_df(800, 1.0e10);
#ifdef DEBUG
    if (myid==0) {
      char debugname[] = "df.debug";
      halo2->print_df(debugname);
    }
#endif
  }
  MPI_Barrier(MPI_COMM_WORLD);

}

DiskHalo::~DiskHalo()
{
  delete [] disktableP;
  delete [] disktableN;
  delete lwz;
  delete rndU;
  delete rndN;
  delete gen;
  delete halo;
  delete disk;
  if (DF) delete newmod;
}

double DiskHalo::disk_density(double R, double z)
{
  double q = 1.0/cosh(z/scaleheight);
  return dmass*disk->get_density(R)*q*q*0.5/scaleheight;
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
  double pos[3], pos1[3], massp, massp1;

  for (int k=0; k<3; k++) pos[k] = pos1[k] = 0.0;
  massp = massp1 = 0.0;

  Particle p;

  p.mass = hmass*(mtot-mmin)/nhalo;

  model = halo;

  for (int i=0; i<npart; i++) {
    targetmass = mmin + (mtot-mmin)*(*rndU)();
    r = zbrent(mass_func, rmin, rmax, tol);
    phi = 2.0*M_PI*(*rndU)();
    costh = 2.0*(*rndU)() - 1.0;
    sinth = sqrt(1.0 - costh*costh);

    p.pos[0] = r*sinth*cos(phi);
    p.pos[1] = r*sinth*sin(phi);
    p.pos[2] = r*costh;

    if (com_cov) {
      massp1 += p.mass;
      for (int k=0; k<3; k++) pos1[k] += p.mass*p.pos[k];
    }

    phalo.push_back(p);
  }

  if (com_cov) {

    MPI_Allreduce(&massp1, &massp, 1, MPI_DOUBLE, MPI_SUM, MPI_COMM_WORLD);
    MPI_Allreduce(pos1, pos, 3, MPI_DOUBLE, MPI_SUM, MPI_COMM_WORLD);

    if (massp>0.0) {
      for (int k=0; k<3; k++) pos[k] /= massp;
    }

    for (auto &p : phalo) {
      for (int k=0; k<3; k++) p.pos[k] -= pos[k];
    }

  }

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

    if (com_cov) {
      massp1 += p.mass;
      for (int k=0; k<3; k++) pos1[k] += p.mass*p.pos[k];
    }

    pdisk.push_back(p);
  }


  if (com_cov) {

    MPI_Allreduce(&massp1, &massp, 1, MPI_DOUBLE, MPI_SUM, MPI_COMM_WORLD);
    MPI_Allreduce(pos1, pos, 3, MPI_DOUBLE, MPI_SUM, MPI_COMM_WORLD);

    if (massp>0.0) {
      for (int k=0; k<3; k++) pos[k] /= massp;
    }

    for (auto &p : pdisk) {
      for (int k=0; k<3; k++) p.pos[k] -= pos[k];
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
  double p0;

  //                   This is the table radius
  //                   ------------------------
  if (sqrt(R*R+z*z) <= M_SQRT1_2*expandd->get_ascale()*expandd->RMAX) {
    
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

  double r, phi, pot, fr, fz, fp;
  
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
  double dP;
  
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
  double dP;

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
  double dP;
  
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
    cout << "Process " << myid << " epi range error: phi=" << phi 
	 << "  cp=" << cp[0] << ", " << cp[1] 
	 << "  cr=" << cr[0] << ", " << cr[1] 
	 << "  ans=" << ans 
	 << "  ep1=" << epitable[iphi1][ir  ]
	 << "  ep2=" << epitable[iphi1][ir+1]
	 << "  ep3=" << epitable[iphi2][ir  ]
	 << "  ep4=" << epitable[iphi2][ir+1]
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

  dP = 2.0*M_PI/NDP;
  double maxr = 0.0, maxz = 0.0;
  for (auto &p : part) {
    maxr = max<double>(maxr, sqrt(p.pos[0]*p.pos[0] + p.pos[1]*p.pos[1]));
    maxz = max<double>(maxz, fabs(p.pos[2]));
  }

  maxz = max<double>(scaleheight*SHFACTOR, maxz);

  if (myid==0) {
    cout << "Table disk epi: RDMIN=" << RDMIN << " maxr=" << maxr << endl;
  }

  dR = (log(maxr) - log(RDMIN))/(NDR-1);
  dZ = maxz/(NDZ-1);

  int NDZ1 = NDZ*NDZF;
  double dZ1 = maxz/(NDZ1-1);

  double R, r, x, y, z, phi, fr, fz, fp, theta;
  double pot, dens, potl, dpr, dpt, dpp, dpdz;

  // Add no force if no component exists
  pot = fr = fz = fp = dens = potl = dpr = dpt = dpp = 0.0;

  Vector workP(0, NDZ1-1);
  Vector workN(0, NDZ1-1);
  Vector workA(0, NDZ1-1);
  Vector workZ(0, NDZ1-1);

  Vector workR(0, NDR-1);
  Vector workE(0, NDR-1);
  Vector workE2(0, NDR-1);

				// For debugging
  Matrix workD(0, 4, 0, NDR-1);

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

      disk_eval(R, 0.0, phi, pot, fr, fz, fp);
      if (expandh)
	expandh->determine_fields_at_point(R, 0.5*M_PI, phi,
					   &dens, &potl, &dpr, &dpt, &dpp);
      workE[j] = -fr + dpr;
      if (i==0) {
	workD[3][j] = -fr;
	workD[4][j] = dpr;
      }

      for (int k=0; k<NDZ1; k++) {

	z = workZ[k] = dZ1*k;

	r = sqrt(R*R + z*z) + MINDOUBLE;

				// Positive

	disk_eval(R, z, phi, pot, fr, fz, fp);
	theta = acos(z/(r + MINDOUBLE));

	if (expandh)
	  expandh->determine_fields_at_point(R, theta, phi,
					     &dens, &potl, &dpr, &dpt, &dpp);

	dpdz = -fz + dpr*z/r + dpt*R*R/(r*r*r);

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

	workN[k] = disk_density(R, z) * dpdz;
      }

				// Integrate positive
      // Splsum(workZ, workP, workA);
      Trapsum(workZ, workP, workA);
      for (int k=0; k<NDZ; k++)
	disktableP[i][j][k] = max<double>(workA[NDZ1-1] - workA[k*NDZF], 
					  MINDOUBLE);

      if (i==ibeg[myid] && j==0) {
	ostringstream ofile;
	ofile << "intgr_disk_P.d" << myid;
	ofstream dout(ofile.str().c_str());
	for (int k=0; k<NDZ1; k++) 
	  dout << setw(15) << workZ[k] 
	       << setw(15) << workP[k] 
	       << setw(15) << workA[k]
	       << "\n";
	dout.close();
      }
				// Integrate negative
      // Splsum(workZ, workN, workA);
      Trapsum(workZ, workN, workA);
      for (int k=0; k<NDZ; k++)
	disktableN[i][j][k] = max<double>(workA[NDZ1-1] - workA[k*NDZF], 
					  MINDOUBLE);

      if (i==ibeg[myid] && j==0) {
	ostringstream ofile;
	ofile << "intgr_disk_N.d" << myid;
	ofstream dout(ofile.str().c_str());
	for (int k=0; k<NDZ1; k++) 
	  dout << setw(15) << workZ[k] 
	       << setw(15) << workN[k] 
	       << setw(15) << workA[k]
	       << "\n";
	dout.close();
      }

    }

				// Compute epicylic freqs
    // Spline(workR, workE, -1.0e30, -1.0e30, workE2);
    for (int j=0; j<NDR; j++) {
      // Splint2(workR, workE, workE2, workR[j], dum, epitable[i][j]);
      epitable[i][j] = drv2(workR[j], workR, workE);
      if (i==0) workD[0][j] = epitable[0][j];
      epitable[i][j] += 3.0*workE[j];
      if (i==0) workD[1][j] = epitable[0][j];
      epitable[i][j] /= exp(workR[j]);
      if (i==0) workD[2][j] = epitable[0][j];
    }
    
  }

				// Update tables on all nodes
  for (int k=0; k<numprocs; k++) {
    for (int i=ibeg[k]; i<iend[k]; i++) {
      MPI_Bcast(&epitable[i][0], NDR, MPI_DOUBLE, k, MPI_COMM_WORLD);
      for (int j=0; j<NDR; j++) {
	MPI_Bcast(&disktableP[i][j][0], NDZ, MPI_DOUBLE, k, MPI_COMM_WORLD);
	MPI_Bcast(&disktableN[i][j][0], NDZ, MPI_DOUBLE, k, MPI_COMM_WORLD);
      }
    }
  }

				// Check solution
  if (myid==0 && expandh) {
    ofstream out("ep_test.dat");
    out.setf(ios::scientific);
    out.precision(4);

    const double DR = 0.01;

    double vr2, r, r1, r2, deriv, rho, lhs, rhs;

    for (int j=0; j<NDR; j++) {
      r = RDMIN*exp(dR*j);
      r1 = r*(1.0 + DR);
      r2 = r*(1.0 - DR);

      vr2 = get_disp(0.0, r, 0.0);
      rho = halo->get_density(r);

      deriv = (get_disp(0.0, r1, 0.0)*halo->get_density(r1) - 
	       get_disp(0.0, r2, 0.0)*halo->get_density(r2) ) /	(r1 - r2);
      
      lhs = halo->get_mass(r);
      rhs = -r*r*deriv/rho;

      out << setw(14) << r	             // 1
	  << setw(14) << epitable[0][j]      // 2
	  << setw(14) << workR[j]            // 3
	  << setw(14) << workE[j]            // 4
	  << setw(14) << workD[0][j]         // 5
	  << setw(14) << workD[1][j]         // 6
	  << setw(14) << workD[2][j]         // 7
	  << setw(14) << workD[3][j]         // 8
	  << setw(14) << workD[4][j]         // 9
	  << setw(14) << vr2                 // 10
	  << setw(14) << rho                 // 11
	  << setw(14) << deriv               // 12
	  << setw(14) << lhs                 // 13
	  << setw(14) << rhs                 // 14
	  << setw(14) << lhs - rhs           // 15
	  << endl;
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
  MPI_Allreduce(&epirmin, &epiRmin, 1, MPI_DOUBLE, MPI_MAX, MPI_COMM_WORLD);
  if (myid==0) cout << "Epitable: Rmin=" << epiRmin << endl;

  if (myid==0)
  {
    ofstream out("ep_disk.dat");
    for (int i=0; i<NDP; i++) {
      for (int j=0; j<NDR; j++)
	out << setw(10) << epitable[i][j];
      out << endl;
    }
  }

  if (myid==0)
  {
    ofstream out("table_disk.dat");
    out.setf(ios::scientific);
    out.precision(2);
    
    for (int k=0; k<NDZ; k++)
      out << setw(10) << dZ*k;
    out << endl << endl;

    for (int j=0; j<NDR; j++) {

      out << setw(10) << RDMIN*exp(dR*j);
      for (int k=0; k<NDZ; k++)
	out << setw(10) << disktableP[0][j][k];
      out << endl;

      out << setw(10) << RDMIN*exp(dR*j);
      for (int k=0; k<NDZ; k++)
	out << setw(10) << disktableN[0][j][k];
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

  ans = sqrt(max<double>(0.0, vc*vc + 
			 vr_disp(xp, yp, 0.0) * 
			 (log(t1)-log(t2))/(log(R1)-log(R2))));
  
				// For debugging
  if (std::isnan(ans)) {
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
  double r = sqrt(xp*xp + yp*yp);
  double sigmar = 3.36*dmass*disk->get_density(r)*Q/sqrt(epi(xp,yp,zp));
  return sigmar*sigmar;
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
  double vz, vr, vp, R, x, y, z;
  double vel[3], vel1[3], massp, massp1;

  for (int k=0; k<3; k++) vel[k] = vel1[k] = 0.0;
  massp = massp1 = 0.0;

  Normal rn(0.0, 1.0, gen);

				 // Better to make a 2-d table
  table_disk(part);
  
  for (auto &p : part) {

				// From solution to Jeans' equations in
				// cylindrical coordinates
    x = p.pos[0];
    y = p.pos[1];
    z = p.pos[2];
    R = sqrt(x*x + y*y) + MINDOUBLE;

    vvZ = get_dispdz(x, y, z);
    vvR = vr_disp(x, y, z);
    vvP = vvR;

				 // For safety; should only be a problem
				 // on extrapolating the range
    vvZ = max<double>(vvZ, MINDOUBLE);
    vvR = max<double>(vvR, MINDOUBLE);
    vvP = max<double>(vvP, MINDOUBLE);
    
    vz   = rn()*sqrt(vvZ);
    vr   = rn()*sqrt(vvR);
    vp   = rn()*sqrt(vvP) + vphi(x, y, z);
    
    p.vel[0] = vr*x/R - vp*y/R;
    p.vel[1] = vr*y/R + vp*x/R;
    p.vel[2] = vz;

    if (com_cov) {
      massp1 += p.mass;
      for (int k=0; k<3; k++) vel1[k] += p.mass*p.vel[k];
    }

  }

  if (com_cov) {

    MPI_Allreduce(&massp1, &massp, 1, MPI_DOUBLE, MPI_SUM, MPI_COMM_WORLD);
    MPI_Allreduce(vel1, vel, 3, MPI_DOUBLE, MPI_SUM, MPI_COMM_WORLD);

    if (massp>0.0) {
      for (int k=0; k<3; k++) vel[k] /= massp;
    }

    for (auto &p : part) {
      for (int k=0; k<3; k++) p.vel[k] -= vel[k];
    }

  }

}

void DiskHalo::
table_halo(vector<Particle>& part)
{
  if (halotable.getchigh() == NHR-1) return;
  
  halotable.setsize(0, NHT-1, 0, NHR-1);
  
  dc = 2.0/(NHT-1);
  double r2, maxr = 0.0;
  for (auto &p : part) {
    r2 = 0.0;
    for (int k=0; k<3; k++) r2 += p.pos[k]*p.pos[k];
    maxr = max<double>(maxr, sqrt(r2));
  }
  
  dr = (log(max<double>(RHMAX, maxr)) - log(RHMIN))/(NHR-1);
  
  
  double x, y, z, theta, r, fr, fz, fp, pot, costh, sinth;
  double dens, potl, dpr, dpt, dpp, dpdr;
  Vector work(0, NHR-1);
  Vector workR(0, NHR-1);
  Vector workA(0, NHR-1);
  
  // If no disk, add no force
  pot = fr = fz = fp = 0.0;
  
  // Compute this table in parallel
  
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
  for (int k=0; k<numprocs; k++) {
    for (int i=ibeg[k]; i<iend[k]; i++) {
      MPI_Bcast(&halotable[i][0], NHR, MPI_DOUBLE, k, MPI_COMM_WORLD);
    }
  }
  
  // DEBUG
  
  if (myid==0 && expandh) {
    ofstream out("table_halo.dat");
    out.setf(ios::scientific);
    out.precision(2);
    
    for (int k=0; k<NHR; k++)
      out << setw(14) << RHMIN*exp(dr*k);
    out << endl << endl;
    
    for (int j=0; j<NHT; j++) {
      out << setw(14) << -1.0 + dc*j;
      for (int k=0; k<NHR; k++)
	out << setw(14) << halotable[j][k];
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
  if (halotable.getchigh() != NHR-1) {
    cerr << "DiskHalo::get_disp: must call table_halo first\n";
    MPI_Abort(MPI_COMM_WORLD, 100);
    exit(0);
  }
  
  int it, ir;
  double r, lr, t, ct[2], cr[2], resv;
  
  // Polar angle
  r = sqrt(xp*xp + yp*yp + zp*zp);
  t = zp/(r + MINDOUBLE) + 1.0;
  
  it = (int)( t/dc );
  it = min<int>( it, NHT-2 );
  ct[1] = (t - dc*it)/dc;
  ct[0] = 1.0 - ct[1];
  // Polar radius
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
  
  for (auto &p : part) {
    
    r = sqrt(p.pos[0]*p.pos[0] + 
	     p.pos[1]*p.pos[1] +
	     p.pos[2]*p.pos[2]);
    
				// Reset success flag
    nok = 1;
    
				// Use Eddington
    
    if (DF && 0.5*(1.0+erf((r-R_DF)/DR_DF)) > (*rndU)()) {
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
      vr = sqrt(max<double>(v2r, MINDOUBLE));
      for (int k=0; k<3; k++) p.vel[k] = vr*(*rndN)();
    }
    
    if (com_cov) {
      massp1 += p.mass;
      for (int k=0; k<3; k++) vel1[k] += p.mass*p.vel[k];
    }
    
  }
  
  if (com_cov) {
    
    MPI_Allreduce(&massp1, &massp, 1, MPI_DOUBLE, MPI_SUM, MPI_COMM_WORLD);
    MPI_Allreduce(vel1, vel, 3, MPI_DOUBLE, MPI_SUM, MPI_COMM_WORLD);
    
    if (massp>0.0) {
      for (int k=0; k<3; k++) vel[k] /= massp;
    }
    
    for (auto &p : part) {
      for (int k=0; k<3; k++) p.vel[k] -= vel[k];
    }
    
  }
  
}

void DiskHalo::write_record(ofstream &out, Particle &p)
{
  out << " " << setw(16) << setprecision(8) << p.mass;
  
  for (int k=0; k<3; k++)
    out << setw(24) << setprecision(15) << p.pos[k] + center_pos[k];
  
  for (int k=0; k<3; k++)
    out << setw(24) << setprecision(15) << p.vel[k] + center_vel[k];
  
  out << endl;
}


void DiskHalo::write_file(ofstream &fou_halo,  ofstream &fou_disk,
			  vector<Particle>& hpart, vector<Particle>& dpart)
{
  int l  = hpart.size();
  int l1 = dpart.size();
  
  vector<Particle> buf(NBUF);
  
				// Make MPI datatype
  MPI_Datatype Particletype;
  MPI_Datatype type[4] = {MPI_UNSIGNED, MPI_DOUBLE, MPI_DOUBLE, MPI_DOUBLE};
  
  // Get displacements
  MPI_Aint disp[4];
  MPI_Get_address(&buf[0].level,	&disp[0]);
  MPI_Get_address(&buf[0].mass,		&disp[1]);
  MPI_Get_address(&buf[0].pos[0],	&disp[2]);
  MPI_Get_address(&buf[0].vel[0],	&disp[3]);
  
  for (int i=3; i>=0; i--) disp[i] -= disp[0];
  
				// Block offsets
  int blocklen[4] = {1, 1, 3, 3};
  
  // Make and register the new type
  MPI_Type_create_struct(4, blocklen, disp, type, &Particletype);
  MPI_Type_commit(&Particletype);
  
  // Get particle totals
  //
  int ndisk=0, nhalo=0;
  MPI_Reduce(&l,  &nhalo, 1, MPI_INT, MPI_SUM, 0, MPI_COMM_WORLD);
  MPI_Reduce(&l1, &ndisk, 1, MPI_INT, MPI_SUM, 0, MPI_COMM_WORLD);
  
  if (myid==0) {
    
    fou_halo.setf(ios::scientific);
    fou_disk.setf(ios::scientific);
    
    fou_halo << nhalo << " " << 0 << " " << 0 << endl;
    fou_disk << ndisk << " " << 0 << " " << 0 << endl;
    
    for (int i=0; i<l; i++)
      write_record(fou_halo, hpart[i]);
    
    for (int i=0; i<l1; i++)
      write_record(fou_disk, dpart[i]);
    
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
      
      MPI_Recv(&imany, 1, MPI_INT, n, 13, MPI_COMM_WORLD, MPI_STATUS_IGNORE);
      ccnt = 0;
      while (ccnt<imany) {
	MPI_Recv(&icur, 1, MPI_INT, n, 14, MPI_COMM_WORLD, MPI_STATUS_IGNORE);
	MPI_Recv(&buf[0], icur, Particletype, n, 15, MPI_COMM_WORLD, MPI_STATUS_IGNORE);
	for (int i=0; i<icur; i++) write_record(fou_disk, buf[i]);
	ccnt += icur;
      }
      
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

	MPI_Send(&l1, 1, MPI_INT, 0, 13, MPI_COMM_WORLD);
	icur = 0;
	while (icur<l1) {
	  ipack = min<int>(l1-icur, NBUF);
	  for (int j=0; j<ipack; j++) buf[j] = dpart[icur+j];
	  MPI_Send(&ipack, 1, MPI_INT, 0, 14, MPI_COMM_WORLD);
	  MPI_Send(&buf[0], ipack, Particletype, 0, 15, MPI_COMM_WORLD);
	  icur += ipack;
	}
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
  
  int c, expected, count=0;
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








