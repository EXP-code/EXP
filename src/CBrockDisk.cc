#include "expand.H"

#include <string>
#include <vector>
#include <iostream>
#include <sstream>

#include <Particle.H>
#include <MixtureBasis.H>
#include <CBrockDisk.H>

CBrockDisk::CBrockDisk(const YAML::Node& conf, MixtureBasis* m) :  AxisymmetricBasis(conf)
{
  id              = "Clutton-Brock two-dimensional disk";

  dof             = 2;		// Two degrees of freedom
  mix             = m;

  rmax            = 100.0;
  scale           = 1.0;
  Lmax            = 4;
  nmax            = 10;

  self_consistent = true;
  coef_dump       = true;

  initialize();

  expcoef .resize(2*(Lmax+1), nmax);
  expcoef1.resize(2*(Lmax+1), nmax);

  expcoef0.resize(nthrds);
  for (auto & v : expcoef0) v.resize((Lmax+1)*(Lmax+1), nmax);
  
  if (pcavar) {
    cc.resize((Lmax+1)*(Lmax+1));
    for (auto & v : cc) v.resize(nmax, nmax);
  
    cc1.resize((Lmax+1)*(Lmax+1));
    for (auto & v : cc1) v.resize(nmax, nmax);
    
    pthread_mutex_init(&cc_lock, NULL);
  }

  // Allocate and compute normalization matrix

  normM.resize(Lmax+1, nmax);
  dend .resize(Lmax+1, nmax);
  work .resize(Lmax+2, nmax);
  
  potd.resize(nthrds);
  dpot.resize(nthrds);

  for (auto & v : potd) v.resize(Lmax+1, nmax);
  for (auto & v : dpot) v.resize(Lmax+1, nmax);

				// Work vectors
  u .resize(nthrds);
  du.resize(nthrds);

  for (auto & v : u ) v.resize(nmax);
  for (auto & v : du) v.resize(nmax);

  for (int l=0; l<=Lmax; l++) {
    for (int n=0; n<nmax; n++) {
      normM (l, n) = norm(n, l);
      sqnorm(l, n) = sqrt(normM(l, n));
    }
  }
    
  // Sin, cos
  
  cosm.resize(nthrds);
  sinm.resize(nthrds);

  for (auto & v : cosm) v.resize(Lmax+1);
  for (auto & v : sinm) v.resize(Lmax+1);

  if (!self_consistent || initializing) determine_coefficients();
}

void CBrockDisk::initialize(void)
{
  try {
    if (conf["rmax"])            rmax   = conf["rmax"].as<double>();
    if (conf["scale"])           scale  = conf["scale"].as<double>();
    if (conf["Lmax"])            Lmax   = conf["Lmax"].as<int>();
    if (conf["nmax"])            nmax   = conf["nmax"].as<int>();
    if (conf["self_consistent"]) self_consistent = conf["self_consistent"].as<bool>();
  }
  catch (YAML::Exception & error) {
    if (myid==0) std::cout << "Error parsing parameters in CBrockDisk: "
			   << error.what() << std::endl
			   << std::string(60, '-') << std::endl
			   << "Config node"        << std::endl
			   << std::string(60, '-') << std::endl
			   << conf                 << std::endl
			   << std::string(60, '-') << std::endl;
    MPI_Finalize();
    exit(-1);
  }

}

CBrockDisk::~CBrockDisk(void)
{
  if (pcavar) {
    pthread_mutex_destroy(&cc_lock);
  }
}

void CBrockDisk::get_acceleration_and_potential(Component* curComp)
{
  cC = curComp;

  //========================================
  // No coefficients for external particles 
  //========================================

  if (use_external) {

    MPL_start_timer();
    determine_acceleration_and_potential();
    MPL_stop_timer();

    use_external = false;

    return;
  }

  //======================================
  // Determine potential and acceleration 
  //======================================

  MPL_start_timer();

  determine_acceleration_and_potential();

  MPL_stop_timer();

  // Clear external potential flag
  use_external = false;
}

void CBrockDisk::determine_coefficients(void)
{
  int compute;

  if (!self_consistent && !initializing) return;

  if (pcavar) compute = !(this_step%npca);

				// Clean
  expcoef .setZero();
  expcoef1.setZero();
  for (auto & v : expcoef0) v.setZero();
  if (pcavar && compute) {
    for (auto & v : cc1) v.setZero();
  }

  use0 = 0;
  use1 = 0;

  exp_thread_fork(true);

  for (int i=0; i<nthrds; i++) {
    use1 += use[i];
    for (int l=0; l<=Lmax*(Lmax+2); l++)
      for (int n=0; n<nmax; n++)
	expcoef1(l, n) += expcoef0[i](l, n);
  }

  MPI_Allreduce ( &use1, &use0,  1, MPI_INT, MPI_SUM, MPI_COMM_WORLD);
  used = use0;

  if (!pcavar) {
    MPI_Allreduce ( expcoef1.data(), expcoef.data(),
		    expcoef.size(), MPI_DOUBLE, MPI_SUM, MPI_COMM_WORLD);
  }

  if (pcavar) {

    parallel_gather_coefficients();

    if (myid == 0) pca_hall(compute);

    parallel_distribute_coefficients();

  }

}

void * CBrockDisk::determine_coefficients_thread(void * arg)
{
  double pos[3], xx, yy, zz, r, r2, phi, rs, fac1, fac2, mass;

  int nbodies = cC->Number();
  int id = *((int*)arg);
  int nbeg = nbodies*id/nthrds;
  int nend = nbodies*(id+1)/nthrds;
  double adb = component->Adiabatic();

  PartMapItr it=cC->Particles().begin();
  unsigned long j;

  vector<double> ctr(3, 0.0);
  if (mix) mix->getCenter(ctr);

  for (int i=0; i<nbeg; i++) it++;
  for (int i=nbeg; i<nend; i++) {

    j = (it++)->first;

    if (component->freeze(j)) // frozen particles do not respond
      continue;

    mass = cC->Mass(j)* adb;

    if (mix) {
      for (int k=0; k<3; k++) 
	pos[k] = cC->Pos(j, k, Component::Local) - ctr[k];
    } else {
      for (int k=0; k<3; k++) 
	pos[k] = cC->Pos(j, k, Component::Local | Component::Centered);
    }

    xx = pos[0];
    yy = pos[1];
    zz = pos[2];

    r2 = (xx*xx + yy*yy);
    r = sqrt(r2) + DSMALL;

    if (r<=rmax) {
      use1++;
      phi = atan2(yy,xx);
      rs = r/scale;
	
      
      sinecosine_R(Lmax, phi, cosm[id], sinm[id]);
      get_potl(Lmax, nmax, rs, potd[id]);
      
				// l loop 

      for (int n=0; n<nmax; n++) {
	expcoef0[id](0, n) += potd[id](0, n)*mass/normM(0, n);
	if (pcavar && compute) {
	  pthread_mutex_lock(&cc_lock);
	  for (int nn=n; nn<nmax; nn++)
	    cc1[0](n, nn) += potd[id](0, n)*potd[id](0, nn)*mass/
	      (normM(0, n)*normM(0, nn));
	  pthread_mutex_unlock(&cc_lock);
	}
      }
	
      for (int l=1; l<=Lmax; l++) {

	fac1 = cosm[id][l];
	fac2 = sinm[id][l];
	
	for (int n=0; n<nmax; n++) {
	  expcoef0[id](2*l - 1, n) +=  potd[id](l, n)*fac1*mass/normM(l, n);
	  expcoef0[id](2*l    , n) +=  potd[id](l, n)*fac2*mass/normM(l, n);
	  if (pcavar && compute) {
	    pthread_mutex_lock(&cc_lock);
	    for (int nn=n; nn<nmax; nn++) {
	      cc1[2*l - 1](n, nn) += potd[id](l, n)*potd[id](l, nn)*fac1*fac1*
		mass/(normM(l, n)*normM(l, nn));
	      cc1[2*l    ](n, nn) += potd[id](l, n)*potd[id](l, nn)*fac2*fac2*
		mass/(normM(l, n)*normM(l, nn));
	    }
	    pthread_mutex_unlock(&cc_lock);
	  }
	}
      }
    }
  }

  return (NULL);
}



void CBrockDisk::determine_acceleration_and_potential(void)
{
  exp_thread_fork(false);
}

void * CBrockDisk::determine_acceleration_and_potential_thread(void * arg)
{
  int l;
  double r,rs,fac,phi,dp;
  double potr,potl,pott,potp,p,pc,dpc,ps,dps;

  double pos[3];
  double xx, yy, zz, mfactor = 1.0;

  int nbodies = cC->Number();
  int id = *((int*)arg);
  int nbeg = nbodies*id/nthrds;
  int nend = nbodies*(id+1)/nthrds;

  PartMapItr it=cC->Particles().begin();
  unsigned long j;

  vector<double> ctr(3, 0.0);
  if (mix) mix->getCenter(ctr);

  for (int i=0; i<nbeg; i++) it++;
  for (int i=nbeg; i<nend; i++) {

    j = it->first;

    if (mix) {
      if (use_external) {
	cC->Pos(pos, i, Component::Inertial);
	component->ConvertPos(pos, Component::Local);
      } else
	cC->Pos(pos, i, Component::Local);

      mfactor = mix->Mixture(pos);
    } else {
      if (use_external) {
	cC->Pos(pos, i, Component::Inertial);
	component->ConvertPos(pos, Component::Local | Component::Centered);
      } else
	cC->Pos(pos, i, Component::Local | Component::Centered);
    }	
    
    xx = pos[0] - ctr[0];
    yy = pos[1] - ctr[1];
    zz = pos[2] - ctr[2];

    r = sqrt(xx*xx + yy*yy) + DSMALL;

    rs = r/scale;
    phi = atan2(yy,xx);

    sinecosine_R(Lmax, phi, cosm[id], sinm[id]);
    get_dpotl(Lmax, nmax, rs, potd[id], dpot[id]);
    get_pot_coefs_safe(0, expcoef.row(0), p, dp, potd[id], dpot[id]);

    potl = p;
    potr = dp;
    pott = potp = 0.0;
				// l loop
    
    for (int l=1; l<=Lmax; l++) {

      get_pot_coefs_safe(l, expcoef.row(2*l - 1), pc, dpc, potd[id], dpot[id]);
      get_pot_coefs_safe(l, expcoef.row(2*l    ), ps, dps, potd[id], dpot[id]);
      potl += pc*cosm[id][l] + ps*sinm[id][l];
      potr += dpc*cosm[id][l] + dps*sinm[id][l];
      potp += (-pc*sinm[id][l] + ps*cosm[id][l])*l;
    }

    fac = xx*xx + yy*yy;
    
    potr *= mfactor/(scale*scale);
    potl *= mfactor/scale;
    potp *= mfactor/scale;

    cC->AddAcc(j, 0, -potr*xx/r);
    cC->AddAcc(j, 1, -potr*yy/r);
    cC->AddAcc(j, 2, -potr*zz/r);
    if (fac > DSMALL) {
      cC->AddAcc(j, 0,  potp*yy/fac);
      cC->AddAcc(j, 1, -potp*xx/fac);
    }

    if (use_external)
      cC->AddPotExt(j, potl);
    else
      cC->AddPot(j, potl);

    it++;
  }

  return (NULL);
}

void 
CBrockDisk::determine_fields_at_point_sph(double r, double theta, double phi,
					  double *tdens0, double *tpotl0, 
					  double *tdens, double *tpotl, 
					  double *tpotr, double *tpott, 
					  double *tpotp)

{
  determine_fields_at_point_polar(r, phi, tdens0, tpotl0, tdens, tpotl, tpotr, tpotp);
  *tpott = 0.0;
}


void 
CBrockDisk::determine_fields_at_point_cyl(double r, double z, double phi,
					  double *tdens0, double *tpotl0, 
					  double *tdens, double *tpotl, 
					  double *tpotr, double *tpott, 
					  double *tpotp)

{
  determine_fields_at_point_polar(r, phi, tdens0, tpotl0, tdens, tpotl, tpotr, tpotp);
  *tpott = 0.0;
}

void 
CBrockDisk::determine_fields_at_point(double x, double y, double z,
				      double *tdens0, double *tpotl0, 
				      double *tdens, double *tpotl, 
				      double *tpotX, double *tpotY, 
				      double *tpotZ)

{
  double R   = sqrt(x*x + y*y);
  double phi = atan2(y, x);
  double cph = cos(phi), sph = sin(phi);
  double tpotR, tpotP;

  determine_fields_at_point_polar(R, phi, tdens0, tpotl0, tdens, tpotl,
				  &tpotR, &tpotP);

  *tpotZ = tpotR*cph - tpotP*sph ;
  *tpotY = tpotR*sph + tpotP*cph ;
  *tpotZ = 0.0;
}


void CBrockDisk::determine_fields_at_point_polar
(
 double r, double phi,
 double *tdens0, double *tpotl0,
 double *tdens, double *tpotl, double *tpotr, double *tpotp
 )
{
  int l;
  double rs,dp;
  double potr,potl,potp,p,pc,dpc,ps,dps,dens;

  rs = r/scale;

  sinecosine_R(Lmax, phi, cosm[0], sinm[0]);

  get_dens(Lmax, nmax, rs, dend);
  get_dpotl(Lmax, nmax, rs, potd[0], dpot[0]);

  get_dens_coefs(0, expcoef.row(0), dens);

  get_pot_coefs(0, expcoef.row(0), p, dp);

  potl = p;
  potr = dp;
  potp = 0.0;
  
  *tdens0 = dens;
  *tpotl0 = potl;
      
  // l loop
    
  for (l=1; l<=Lmax; l++) {
    
    get_dens_coefs(l,expcoef.row(2*l - 1), pc);
    get_dens_coefs(l,expcoef.row(2*l    ), ps);
    dens += pc*cosm[0][l] + ps*sinm[0][l];
    
    get_pot_coefs(l,expcoef.row(2*l - 1), pc, dpc);
    get_pot_coefs(l,expcoef.row(2*l    ), ps, dps);
    potl += pc*cosm[0][l] + ps*sinm[0][l];
    potr += dpc*cosm[0][l] + dps*sinm[0][l];
    potp += (-pc*sinm[0][l] + ps*cosm[0][l])*l;
  }

  *tdens0 /= scale*scale*scale;
  *tpotl0 /= scale;
      
  *tdens = dens/(scale*scale*scale);
  *tpotl = potl/scale;
  *tpotr = potr/(scale*scale);
  *tpotp = potp/scale;
}


void CBrockDisk::get_pot_coefs(int l, const Eigen::VectorXd& coef,
			       double& p, double& dp)
{
  double pp, dpp;

  pp = dpp = 0.0;

  for (int i=0; i<nmax; i++) {
    pp  += potd[0](l, i) * coef[i];
    dpp += dpot[0](l, i) * coef[i];
  }

  p = -pp;
  dp = -dpp;
}


void CBrockDisk::get_pot_coefs_safe
(int l, const Eigen::VectorXd& coef, double& p, double& dp,
 Eigen::MatrixXd& potd1, Eigen::MatrixXd& dpot1)
{
  double pp, dpp;

  pp = dpp = 0.0;

  for (int i=0; i<nmax; i++) {
    pp  += potd1(l, i) * coef[i];
    dpp += dpot1(l, i) * coef[i];
  }

  p  = -pp;
  dp = -dpp;
}


void CBrockDisk::get_dens_coefs
(int l, const Eigen::VectorXd& coef, double& p)
{
  double pp = 0.0;

  for (int i=0; i<nmax; i++)
    pp  += dend(l, i) * coef[i];

  p = pp;
}


				// Get potential functions by recursion

void CBrockDisk::get_dpotl(int lmax, int nmax, double r,
			   Eigen::MatrixXd& p, Eigen::MatrixXd& dp)
{
  double r2 = r*r;
  double fac = 1.0/(1.0 + r2);
  double fac1 = (r2 - 1.0)*fac;
  double cur, curl1, curl2, cur0 = sqrt(fac), rcum = 1.0;

  for (int l=0; l<=lmax+1; l++) {
    cur = cur0;

    p(l, 0) = cur*rcum;
    curl1 = 0.0;

    for (int nn=0; nn<nmax-1; nn++) {
      curl2 = curl1;
      curl1 = cur;
      cur = (2.0 + (double)(2*l-1)/(nn+1))*fac1*curl1 - 
	(1.0 + (double)(2*l-1)/(nn+1))*curl2;
      p(l, nn+1) = cur*rcum;
    }
    cur0 *= fac*(2*(l+1) - 1);
    rcum *= r;
  }


  for (int l=0; l<=lmax; l++) {

    dp(l, 0) = p(l, 0)*l/r - p(l+1, 0);
    if (nmax<1) break;

    dp(l, 1) = p(l, 1)*l/r - (p(l+1, 1) - 2*p(l+1, 0));
    if (nmax<2) break;

    for (int nn=1; nn<nmax-1; nn++) {
      dp(l, nn+1) = p(l, nn+1)*l/r - 
	( p(l+1, nn+1) - 2*p(l+1, nn) + p(l+1, nn-1) );
    }
  }

}

void CBrockDisk::get_potl(int lmax, int nmax, double r, Eigen::MatrixXd& p)
{
  double r2 = r*r;
  double fac = 1.0/(1.0 + r2);
  double fac1 = (r2 - 1.0)*fac;
  double cur, curl1, curl2, cur0 = sqrt(fac), rcum = 1.0;

  for (int l=0; l<=lmax+1; l++) {
    cur = cur0;

    work(l, 0)  = cur;
    p(l, 0) = cur*rcum;
    curl1 = 0.0;

    for (int nn=0; nn<nmax-1; nn++) {
      curl2 = curl1;
      curl1 = cur;
      cur = (2.0 + (double)(2*l-1)/(nn+1))*fac1*curl1 - 
	(1.0 + (double)(2*l-1)/(nn+1))*curl2;
      p(l, nn+1) = cur*rcum;
    }
    cur0 *= fac*(2*(l+1) - 1);
    rcum *= r;
  }

}


void CBrockDisk::get_dens(int lmax, int nmax, double r, Eigen::MatrixXd& d)
{

  double r2 = r*r;
  double fac = 1.0/(1.0 + r2);
  double fac1 = (r2 - 1.0)*fac;
  double cur, curl1, curl2, cur0 = sqrt(fac);
  double rcum = 0.5/M_PI;
  
  for (int l=0; l<=lmax+1; l++) {
    cur = cur0;

    work(l, 1) = cur;
    curl1 = 0.0;

    for (int nn=0; nn<nmax-1; nn++) {
      curl2 = curl1;
      curl1 = cur;
      cur = (2.0 + (double)(2*l-1)/(nn+1))*fac1*curl1 - 
	(1.0 + (double)(2*l-1)/(nn+1))*curl2;
      work(l, nn+1) = cur;
    }
    cur0 *= fac*(2*(l+1) - 1);
  }

  for (int l=0; l<=lmax; l++) {
    d(l, 0) = work(l+1, 0)*rcum;
    d(l, 1) = work(l+1, 1)*rcum;
    for (int nn=1; nn<nmax-1; nn++)
      d(l, nn+1) = (work(l+1, nn+1) - work(l+1, nn-1))*rcum;

    rcum *= r;
  }

}

void CBrockDisk::get_potl_dens(int lmax, int nmax, double r, 
			       Eigen::MatrixXd& p, Eigen::MatrixXd& d)
{
  double r2 = r*r;
  double fac = 1.0/(1.0 + r2);
  double fac1 = (r2 - 1.0)*fac;
  double cur, curl1, curl2, cur0 = sqrt(fac);
  double  rcump = 1.0, rcumd = 0.5/M_PI;

  for (int l=0; l<=lmax+1; l++) {
    cur = cur0;

    work(l, 0) = cur;
    curl1 = 0.0;

    for (int nn=0; nn<nmax-1; nn++) {
      curl2 = curl1;
      curl1 = cur;
      cur = (2.0 + (double)(2*l-1)/(nn+1))*fac1*curl1 - 
	(1.0 + (double)(2*l-1)/(nn+1))*curl2;
      p(l, nn+1) = cur*rcump;
      work(l, nn+1) = cur;
    }
    cur0 *= fac*(2*(l+1) - 1);
    rcump *= r;
  }

  for (int l=0; l<=lmax; l++) {
    d(l, 0) = work(l+1, 0)*rcumd;
    d(l, 1) = work(l+1, 1)*rcumd;
    for (int nn=1; nn<nmax-1; nn++)
      d(l, nn+1) = (work(l+1, nn+1) - work(l+1, nn-1))*rcumd;

    rcumd *= r;
  }

}

double CBrockDisk::norm(int n, int m)
{
  double ans = 1.0;
 
  for (int i=n+1; i<=n+2*m; i++)
    ans *= i;

  return pow(0.5, 2*m+1)*ans;
}

				// Dump coefficients to a file

void CBrockDisk::dump_coefs(ostream& out)
{
  ostringstream sout;
  sout << id;

  char buf[64];
  for (unsigned i=0; i<64; i++) {
    if (i<sout.str().length()) 
      buf[i] = sout.str().c_str()[i];
    else 
      buf[i] = '\0';
  }
  
  out.write((char *)&buf, 64*sizeof(char));

  out.write((char *)&tnow, sizeof(double));
  out.write((char *)&scale, sizeof(double));
  out.write((char *)&nmax, sizeof(int));
  out.write((char *)&Lmax, sizeof(int));
  out.write((char *)expcoef.data(), expcoef.size()*sizeof(double));
}

