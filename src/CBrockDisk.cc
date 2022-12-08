#include "expand.H"

#include <filesystem>
#include <iostream>
#include <sstream>
#include <vector>
#include <string>

#include <Particle.H>
#include <MixtureBasis.H>
#include <CBrockDisk.H>

const std::set<std::string>
CBrockDisk::valid_keys = {
  "rmax",
  "scale",
  "Lmax",
  "nmax",
  "self_consistent",
  "playback",
  "coefCompute",
  "coefMaster"
};

CBrockDisk::CBrockDisk(Component* c0, const YAML::Node& conf, MixtureBasis* m) :  AxisymmetricBasis(c0, conf)
{
  id              = "Clutton-Brock two-dimensional disk";
  geometry        = cylinder;

  dof             = 2;		// Two degrees of freedom
  mix             = m;

  rmax            = 100.0;
  scale           = 1.0;
  Lmax            = 4;
  nmax            = 10;

  self_consistent = true;
  coef_dump       = true;

  initialize();

  expcoef  = std::make_shared<Eigen::MatrixXd>(2*Lmax+1, nmax);
  expcoef1 = std::make_shared<Eigen::MatrixXd>(2*Lmax+1, nmax);

  expcoef0.resize(nthrds);
  for (auto & v : expcoef0) v.resize(2*Lmax+1, nmax);
  
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

  // Needed for dpot recursion--------+
  //                                  |
  //                                  v
  for (auto & v : potd) v.resize(Lmax+2, nmax);
  for (auto & v : dpot) v.resize(Lmax+1, nmax);
  //                                  ^
  //                                  |
  // Otherwise, we only need m<=Lmax--+

  // Work vectors
  //
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
  // Remove matched keys
  //
  for (auto v : valid_keys) current_keys.erase(v);
  
  try {
    if (conf["rmax"])            rmax   = conf["rmax"].as<double>();
    if (conf["scale"])           scale  = conf["scale"].as<double>();
    if (conf["Lmax"])            Lmax   = conf["Lmax"].as<int>();
    if (conf["nmax"])            nmax   = conf["nmax"].as<int>();
    if (conf["self_consistent"]) self_consistent = conf["self_consistent"].as<bool>();
    if (conf["playback"]) {
      std::string file = conf["playback"].as<std::string>();
				// Check the file exists
      {
	std::ifstream test(file);
	if (not test) {
	  std::cerr << "CBrockDisk: process " << myid << " cannot open <"
		    << file << "> for reading" << std::endl;
	  MPI_Finalize();
	  exit(-1);
	}
      }

      playback = std::make_shared<TwoDCoefs>(file);

      if (playback->nmax != nmax) {
	if (myid==0) {
	  std::cerr << "CBrockDisk: nmax for playback [" << playback->nmax
		    << "] does not match specification [" << nmax << "]"
		    << std::endl;
	}
	MPI_Finalize();
	exit(-1);
      }
      
      if (playback->Lmax != Lmax) {
	if (myid==0) {
	  std::cerr << "CBrockDisk: Lmax for playback [" << playback->Lmax
		    << "] does not match specification [" << Lmax << "]"
		    << std::endl;
	}
	MPI_Finalize();
	exit(-1);
      }

      play_back = true;

      if (conf["coefCompute"]) play_cnew = conf["coefCompute"].as<bool>();

      if (conf["coefMaster"]) coefMaster = conf["coefMaster"].as<bool>();

      if (myid==0) {
	std::cout << "---- Playback is ON for Component " << component->name
		  << " using Force " << component->id << std::endl;
	if (coefMaster)
	  std::cout << "---- Playback will use MPI master" << std::endl;
	if (play_cnew)
	  std::cout << "---- New coefficients will be computed from particles on playback" << std::endl;
      }
    }

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
  if (play_back) {
    determine_coefficients_playback();
    if (play_cnew) determine_coefficients_particles();
  } else {
    determine_coefficients_particles();
  }
}

void CBrockDisk::determine_coefficients_playback(void)
{
  // Do we need new coefficients?
  if (tnow <= lastPlayTime) return;
  lastPlayTime = tnow;

  if (coefMaster) {

    if (myid==0) expcoefP = playback->interpolate(tnow);

    MPI_Bcast(expcoefP->data(), expcoefP->size(), MPI_DOUBLE, 0, MPI_COMM_WORLD);
  } else {
    expcoefP = playback->interpolate(tnow);
  }
}

void CBrockDisk::determine_coefficients_particles(void)
{
  int compute;

  if (!self_consistent && !initializing) return;

  if (pcavar) compute = !(this_step%npca);

				// Clean
  std::fill(expcoef ->data(), expcoef ->data() + expcoef ->size(), 0.0);
  std::fill(expcoef1->data(), expcoef1->data() + expcoef1->size(), 0.0);

  for (auto & v : expcoef0) v.setZero();
  if (pcavar && compute) {
    for (auto & v : cc1) v.setZero();
  }

  use0 = 0;
  use1 = 0;

  exp_thread_fork(true);

  for (int i=0; i<nthrds; i++) {
    use1 += use[i];
    for (int l=0; l<2*Lmax+1; l++)
      for (int n=0; n<nmax; n++)
	(*expcoef1)(l, n) += expcoef0[i](l, n);
  }

  MPI_Allreduce ( &use1, &use0,  1, MPI_INT, MPI_SUM, MPI_COMM_WORLD);
  used = use0;

  if (!pcavar) {
    MPI_Allreduce ( expcoef1->data(), expcoef->data(),
		    expcoef->size(), MPI_DOUBLE, MPI_SUM, MPI_COMM_WORLD);
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
  int id      = *((int*)arg);
  int nbeg    = nbodies*id/nthrds;
  int nend    = nbodies*(id+1)/nthrds;
  double adb  = component->Adiabatic();

  vector<double> ctr(3, 0.0);
  if (mix) mix->getCenter(ctr);

  PartMapItr it = cC->Particles().begin();
  std::advance(it, nbeg);

  for (int i=nbeg; i<nend; i++) {

    unsigned long j = (it++)->first;

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
  if (play_back) std::swap(expcoefP, expcoef);

  exp_thread_fork(false);

  if (play_back) std::swap(expcoef, expcoefP);
}

void * CBrockDisk::determine_acceleration_and_potential_thread(void * arg)
{
  double p, dp;
  double pos[3];
  double mfactor = 1.0;

  int nbodies = cC->Number();
  int id = *((int*)arg);
  int nbeg = nbodies*id/nthrds;
  int nend = nbodies*(id+1)/nthrds;

  std::vector<double> ctr(3, 0.0);
  if (mix) mix->getCenter(ctr);

  PartMapItr it=cC->Particles().begin();
  std::advance(it, nbeg);

  for (int i=nbeg; i<nend; i++) {

    unsigned long j = it->first;

    if (mix) {
      if (use_external) {
	cC->Pos(pos, j, Component::Inertial);
	component->ConvertPos(pos, Component::Local);
      } else
	cC->Pos(pos, j, Component::Local);

      mfactor = mix->Mixture(pos);
    } else {
      if (use_external) {
	cC->Pos(pos, j, Component::Inertial);
	component->ConvertPos(pos, Component::Local | Component::Centered);
      } else
	cC->Pos(pos, j, Component::Local | Component::Centered);
    }	
    
    double xx  = pos[0] - ctr[0];
    double yy  = pos[1] - ctr[1];
    double zz  = pos[2] - ctr[2];

    double r   = sqrt(xx*xx + yy*yy) + DSMALL;

    double rs  = r/scale;
    double phi = atan2(yy,xx);

    sinecosine_R(Lmax, phi, cosm[id], sinm[id]);
    get_dpotl(Lmax, nmax, rs, potd[id], dpot[id]);
    get_pot_coefs_safe(0, expcoef->row(0), p, dp, potd[id], dpot[id]);

    double potl = p;
    double potr = dp;
    double pott = 0.0, potp = 0.0;

				// l loop
    for (int l=1; l<=Lmax; l++) {

      double pc, dpc, ps, dps;

      get_pot_coefs_safe(l, expcoef->row(2*l - 1), pc, dpc, potd[id], dpot[id]);
      get_pot_coefs_safe(l, expcoef->row(2*l    ), ps, dps, potd[id], dpot[id]);

      potl += pc*cosm[id][l]   + ps*sinm[id][l];
      potr += dpc*cosm[id][l]  + dps*sinm[id][l];
      potp += (-pc*sinm[id][l] + ps*cosm[id][l])*l;
    }

    double fac = xx*xx + yy*yy;
    
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
  double p, dp, dens;

  double rs = r/scale;

  sinecosine_R(Lmax, phi, cosm[0], sinm[0]);

  get_dens (Lmax, nmax, rs, dend);
  get_dpotl(Lmax, nmax, rs, potd[0], dpot[0]);

  get_dens_coefs(0, expcoef->row(0), dens);

  get_pot_coefs(0, expcoef->row(0), p, dp);

  double potl = p;
  double potr = dp;
  double potp = 0.0;
  
  *tdens0 = dens;
  *tpotl0 = potl;
      
  // l loop
    
  for (int l=1; l<=Lmax; l++) {
    
    double pc, ps, dpc, dps;

    get_dens_coefs(l,expcoef->row(2*l - 1), pc);
    get_dens_coefs(l,expcoef->row(2*l    ), ps);
    dens += pc*cosm[0][l] + ps*sinm[0][l];
    
    get_pot_coefs(l,expcoef->row(2*l - 1), pc, dpc);
    get_pot_coefs(l,expcoef->row(2*l    ), ps, dps);
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
  double r2   = r*r;
  double fac  = 1.0/(1.0 + r2);
  double fac1 = (r2 - 1.0)*fac;
  double cur0 = sqrt(fac), rcum = 1.0;

  // Phi recursion relation
  //
  // Need up to lmax+1 for dPhi/dr recursion
  //
  for (int l=0; l<=lmax+1; l++) {
    double cur = cur0;

    p(l, 0) = cur*rcum;
    double curl1 = 0.0;

    for (int nn=0; nn<nmax-1; nn++) {
      double curl2 = curl1;
      curl1 = cur;
      cur = (2.0 + (double)(2*l-1)/(nn+1))*fac1*curl1 -
	(1.0 + (double)(2*l-1)/(nn+1))*curl2;
      p(l, nn+1) = cur*rcum;
    }
    cur0 *= fac*(2*(l+1) - 1);
    rcum *= r;
  }


  // dPhi/dR recursion
  //
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
  double r2   = r*r;
  double fac  = 1.0/(1.0 + r2);
  double fac1 = (r2 - 1.0)*fac;
  double cur0 = sqrt(fac), rcum = 1.0;

  for (int l=0; l<=lmax; l++) {
    double cur = cur0;

    work(l, 0)  = cur;
    p(l, 0) = cur*rcum;
    double curl1 = 0.0;

    for (int nn=0; nn<nmax-1; nn++) {
      double curl2 = curl1;
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
  
  for (int l=0; l<=lmax; l++) {
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
  double r2   = r*r;
  double fac  = 1.0/(1.0 + r2);
  double fac1 = (r2 - 1.0)*fac;
  double cur0 = sqrt(fac);
  double rcump = 1.0, rcumd = 0.5/M_PI;

  for (int l=0; l<=lmax; l++) {
    double cur = cur0;

    work(l, 0) = cur;
    double curl1 = 0.0;

    for (int nn=0; nn<nmax-1; nn++) {
      double curl2 = curl1;
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
  // This is a node of simple {key: value} pairs.  More general
  // content can be added as needed.
  //
  YAML::Node node;
  
  node["id"    ] = id;
  node["time"  ] = tnow;
  node["scale" ] = scale;
  node["rmax"  ] = rmax;
  node["nmax"  ] = nmax;
  node["lmax"  ] = Lmax;

  // Serialize the node
  //
  YAML::Emitter y; y << node;
  
  // Get the size of the string
  //
  unsigned int hsize = strlen(y.c_str());
  
  // Write magic #
  //
  out.write(reinterpret_cast<const char *>(&cmagic),   sizeof(unsigned int));
  
  // Write YAML string size
  //
  out.write(reinterpret_cast<const char *>(&hsize),    sizeof(unsigned int));
  
  // Write YAML string
  //
  out.write(reinterpret_cast<const char *>(y.c_str()), hsize);

  // Write coefficient matrix
  //
  out.write((char *)expcoef->data(), expcoef->size()*sizeof(double));
}

void CBrockDisk::dump_coefs_h5(const std::string& file)
{
  // Add the current coefficients
  auto cur = std::make_shared<CoefClasses::CylStruct>();

  cur->time   = tnow;
  cur->geom   = geoname[geometry];
  cur->id     = id;
  cur->mmax   = Lmax;
  cur->nmax   = nmax;

  cur->coefs.resize(Lmax+1, nmax);

  Eigen::VectorXd cos1(nmax), sin1(nmax);
  
  for (int m=0; m<=Lmax; m++) {
    for (int ir=0; ir<nmax; ir++) {
      if (m==0) cur->coefs(m, ir) = {(*expcoef)(2*m, ir), 0.0};
      else cur->coefs(m, ir) = {(*expcoef)(2*m-1, ir), (*expcoef)(2*m, ir)};
    }
  }

  // Add center
  //
  cur->ctr = component->getCenter(Component::Local | Component::Centered);

  // Check if file exists
  //
  if (std::filesystem::exists(file + ".h5")) {
    cylCoefs.clear();
    cylCoefs.add(cur);
    cylCoefs.ExtendH5Coefs(file);
  } else {
    // Copy the YAML config.  We only need this on the first call.
    std::ostringstream sout; sout << conf;
    size_t hsize = sout.str().size() + 1;
    cur->buf = std::shared_ptr<char[]>(new char [hsize]);
    sout.str().copy(cur->buf.get(), hsize); // Copy to CoefStruct buffer

    // Add the name attribute.  We only need this on the first call.
    cylCoefs.setName(component->name);

    // Add the new coefficients and write the new HDF5
    cylCoefs.clear();
    cylCoefs.add(cur);
    cylCoefs.WriteH5Coefs(file);
  }
}

