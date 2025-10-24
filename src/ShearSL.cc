#include <filesystem>
#include <cmath>

#include "expand.H"

#include "Coefficients.H"
#include "biorth1d.H"
#include <ShearSL.H>

const std::set<std::string>
ShearSL::valid_keys = {
  "nmaxx",
  "nmaxy",
  "nmaxz",
  "nminx",
  "nminy",
  "hslab",
  "Lx",
  "Ly",
  "R",
  "Omega",
  "Kappa",
  "zmax"
};

ShearSL::ShearSL(Component* c0, const YAML::Node& conf) : PotAccel(c0, conf)
{
  id = "Shear (Sturm-Liouville)";
  NGRID  = 100;
  nminx  = 0;
  nminy  = 0;
  nmaxx  = 10;
  nmaxy  = 10;
  nmaxz  = 10;
  zmax   = 10.0;
  hslab  = 0.2;
  Omega  = 1.0;
  Kappa  = M_SQRT2;
  R      = 1.0;
  Lx     = 0.2;
  Ly     = 0.2;
  coef_dump = true;

  initialize();

  SLGridSlab::mpi  = 1;
  SLGridSlab::ZBEG = 0.0;
  SLGridSlab::ZEND = 0.1;
  SLGridSlab::H    = hslab;
  
  int nnmax = (nmaxx > nmaxy) ? nmaxx : nmaxy;

  grid = std::make_shared<SLGridSlab>(nnmax, nmaxz, NGRID, zmax);

  imx  = 1+2*nmaxx;
  imy  = 1+2*nmaxy;
  imz  = nmaxz;
  jmax = imx*imy*imz;

  expccof.resize(nthrds);
  for (auto & v : expccof) v.resize(jmax);
    
  dfacx = 2.0*M_PI/Lx;
  dfacy = 2.0*M_PI/Ly;

  normx = 1.0/sqrt(Lx);
  normy = 1.0/sqrt(Ly);

  kfacx = std::complex<double>(0.0, dfacx);
  kfacy = std::complex<double>(0.0, dfacy);
    
  nnmax = (nmaxx > nmaxy) ? nmaxx : nmaxy;

  zpot.resize(nthrds);
  zfrc.resize(nthrds);

  for (auto & v : zpot) v.resize(nmaxz);
  for (auto & v : zfrc) v.resize(nmaxz);
}

ShearSL::~ShearSL()
{
  // Nothing
}

void ShearSL::initialize()
{
  // Remove matched keys
  //
  for (auto v : valid_keys) current_keys.erase(v);
  
  // Assign values from YAML
  //
  try {
    if (conf["nmaxx"])          nmaxx       = conf["nmaxx"].as<int>();
    if (conf["nmaxy"])          nmaxy       = conf["nmaxy"].as<int>();
    if (conf["nmaxz"])          nmaxz       = conf["nmaxz"].as<int>();
    if (conf["nminx"])          nminx       = conf["nminx"].as<int>();
    if (conf["nminy"])          nminy       = conf["nminy"].as<int>();
    if (conf["hslab"])          hslab       = conf["hslab"].as<double>();
    if (conf["zmax" ])          zmax        = conf["zmax" ].as<double>();
    if (conf["Omega"])          Omega       = conf["Omega"].as<double>();
    if (conf["Kappa"])          Kappa       = conf["Kappa"].as<double>();
    if (conf["Lx"   ])          Lx          = conf["Lx"   ].as<double>();
    if (conf["Ly"   ])          Ly          = conf["Ly"   ].as<double>();
    if (conf["R"    ])          R           = conf["R"    ].as<double>();
  }
  catch (YAML::Exception & error) {
    if (myid==0) std::cout << "Error parsing parameters in ShearSL: "
			   << error.what() << std::endl
			   << std::string(60, '-') << std::endl
			   << "Config node"        << std::endl
			   << std::string(60, '-') << std::endl
			   << conf                 << std::endl
			   << std::string(60, '-') << std::endl;
    throw std::runtime_error("ShearSL::initialze: error parsing YAML");
    exit(-1);
  }
}

void ShearSL::determine_coefficients(void)
{
				//  Coefficients are ordered as follows:
				//  n=-nmax,-nmax+1,...,0,...,nmax-1,nmax
				//  in a single array for each dimension
				//  with z dimension changing most rapidly
  
  // Clean 

  for (int i=0; i<nthrds; i++) {
    use[i] = 0;
    expccof[i].setZero();
  }

  exp_thread_fork(true);

  int used1 = 0;
  used = 0;
  for (int i=1; i<nthrds; i++) {
    used1 += use[i];
    for (int j=0; j<jmax; j++) expccof[0][j] += expccof[i][j];
  }
  
  MPI_Allreduce ( &used1, &used,  1, MPI_INT, MPI_SUM, MPI_COMM_WORLD);

  MPI_Allreduce( MPI_IN_PLACE, expccof[0].data(), expccof[0].size(),
		 MPI_CXX_DOUBLE_COMPLEX, MPI_SUM, MPI_COMM_WORLD);

}

void * ShearSL::determine_coefficients_thread(void * arg)
{
  int ix, iy;
  std::complex<double> facx, facy;

  unsigned nbodies = cC->Number();
  int id = *((int*)arg);
  int nbeg = nbodies*id/nthrds;
  int nend = nbodies*(id+1)/nthrds;
  double adb = cC->Adiabatic();

  for (int i=nbeg; i<nend; i++) {
    
				// Increment particle counter
    use[id]++;

				// Truncate to box with sides in
				// [-Lx/2, Lx/2] x [ 0, Ly ]
    
    if (cC->Pos(i, 0) < -0.5*Lx) {
      int delta = floor((0.5*Lx - cC->Pos(i, 0)/Lx));
      cC->AddPos(i, 0, Lx*delta);
    }
    if (cC->Pos(i, 0) > 0.5*Lx) {
      int delta = floor((cC->Pos(i, 0) + 0.5*Lx)/Lx);
      cC->AddPos(i, 0, -Lx*delta);
    }
    
    if (cC->Pos(i, 1) < 0) {
      int delta = floor(-cC->Pos(i, 0)/Ly);
      cC->AddPos(i, 1, Ly*delta);
    }
    if (cC->Pos(i, 0) > Ly) {
      int delta = floor(cC->Pos(i, 1)/Ly);
      cC->AddPos(i, 0, -Ly*delta);
    }
    

				// Recursion multipliers
    auto stepx = std::exp(-kfacx*cC->Pos(i, 0));
    auto stepy = std::exp(-kfacy*cC->Pos(i, 1));
   
				// Initial values
    auto startx = exp(static_cast<double>(nmaxx)*kfacx*cC->Pos(i, 0));
    auto starty = exp(static_cast<double>(nmaxy)*kfacy*cC->Pos(i, 1));
    
    for (facx=startx, ix=0; ix<imx; ix++, facx*=stepx) {
      
      int ii = ix - nmaxx;
      int iix = abs(ii);
      
      for (facy=starty, iy=0; iy<imy; iy++, facy*=stepy) {
	
	int jj = iy - nmaxy;
	int iiy = abs(jj);
	
	if (iix > nmaxx) {
	  cerr << "Out of bounds: iix=" << ii << endl;
	}
	if (iiy > nmaxy) {
	  cerr << "Out of bounds: iiy=" << jj << endl;
	}
	
	double zz = cC->Pos(i, 2, Component::Centered);

	if (iix>=iiy)
	  grid->get_pot(zpot[id], zz, iix, iiy);
	else
	  grid->get_pot(zpot[id], zz, iiy, iix);


	for (int iz=0; iz<imz; iz++) {

	  int indx = imz*(iy + imy*ix) + iz;

                              // |--- density in orthogonal series
                              // |    is 4.0*M_PI rho
                              // v
	  expccof[id][indx] += -4.0*M_PI*cC->Mass(i)*adb*
	    facx*facy*zpot[id][iz+1] * normx * normy;
	}
      }
    }
  }
    
  return (NULL);
}

void ShearSL::get_acceleration_and_potential(Component* C)
{
  cC = C;

  MPL_start_timer();
  exp_thread_fork(false);
  MPL_stop_timer();
}

void * ShearSL::determine_acceleration_and_potential_thread(void * arg)
{
  int ix, iy, iz, iix, iiy, ii, jj, indx;
  std::complex<double> fac, startx, starty, facx, facy, potl, facf;
  std::complex<double> stepx, stepy;
  std::complex<double> accx, accy, accz;

  unsigned nbodies = cC->Number();
  int id = *((int*)arg);
  int nbeg = nbodies*id/nthrds;
  int nend = nbodies*(id+1)/nthrds;
  double zz;

  for (int i=nbeg; i<nend; i++) {
    
    accx = accy = accz = potl = 0.0;
    
				// Recursion multipliers
    auto stepx = exp(kfacx*cC->Pos(i, 0));
    auto stepy = exp(kfacy*cC->Pos(i, 1));

				// Initial values (note sign change)
    auto startx = exp(-static_cast<double>(nmaxx)*kfacx*cC->Pos(i, 0));
    auto starty = exp(-static_cast<double>(nmaxy)*kfacy*cC->Pos(i, 1));
    
    for (facx=startx, ix=0; ix<imx; ix++, facx*=stepx) {
      
				// Compute wavenumber; recall that the
				// coefficients are stored as follows:
				// -nmax,-nmax+1,...,0,...,nmax-1,nmax
      ii = ix - nmaxx;
      iix = abs(ii);
      
      for (facy=starty, iy=0; iy<imy; iy++, facy*=stepy) {
	
	jj = iy - nmaxy;
	iiy = abs(jj);
	
	if (iix > nmaxx) {
	  cerr << "Out of bounds: ii=" << ii << endl;
	}
	if (iiy > nmaxy) {
	  cerr << "Out of bounds: jj=" << jj << endl;
	}
	
	zz = cC->Pos(i, 2, Component::Centered);

	if (iix>=iiy) {
	  grid->get_pot  (zpot[id], zz, iix, iiy);
	  grid->get_force(zfrc[id], zz, iix, iiy);
	}
	else {
	  grid->get_pot  (zpot[id], zz, iiy, iix);
	  grid->get_force(zfrc[id], zz, iiy, iix);
	}

	
	for (iz=0; iz<imz; iz++) {
	  
	  indx = imz*(iy + imy*ix) + iz;
	  
	  fac  = facx*facy*zpot[id][iz+1]*expccof[0][indx] * normx * normy;
	  facf = facx*facy*zfrc[id][iz+1]*expccof[0][indx] * normx * normy;
	  
				// Limit to minimum wave number
	  
	  if (abs(ii)<nminx || abs(jj)<nminy) continue;
	  
	  potl += fac;
	  
	  accx += -dfacx*ii*std::complex<double>(0.0,1.0)*fac;
	  accy += -dfacy*jj*std::complex<double>(0.0,1.0)*fac;
	  accz += -facf;
	  
	}
      }
    }
    
    cC->AddAcc(i, 0, accx.real());
    cC->AddAcc(i, 1, accy.real());
    cC->AddAcc(i, 2, accz.real());
    cC->AddPot(i, potl.real());
  }

  return (NULL);
}

void ShearSL::dump_coefs_h5(const std::string& file)
{
  // Add the current coefficients
  auto cur = std::make_shared<CoefClasses::SlabStruct>();

  cur->time   = tnow;
  cur->geom   = geoname[geometry];
  cur->id     = id;
  cur->time   = tnow;
  cur->nmaxx  = nmaxx;
  cur->nmaxy  = nmaxy;
  cur->nmaxz  = nmaxz;


  cur->coefT.resize({2*nmaxx+1, 2*nmaxy+1, nmaxz});

  for (int nx=0; nx<2*nmaxx; nx++) {
    for (int ny=0; ny<2*nmaxy; ny++) {
      for (int nz=0; nz<nmaxz; nz++) {
	int indx = imz*(ny + imy*nx) + nz;
	cur->coefT(nx, ny, nz) = expccof[0][indx];
      }
    }
  }

  // Check if file exists
  //
  if (std::filesystem::exists(file)) {
    slabCoefs.clear();
    slabCoefs.add(cur);
    slabCoefs.ExtendH5Coefs(file);
  }
  // Otherwise, extend the existing HDF5 file
  //
  else {
    // Copy the YAML config.  We only need this on the first call.
    std::ostringstream sout; sout << conf;
    cur->buf = sout.str();	// Copy to CoefStruct buffer

    // Add the name attribute.  We only need this on the first call.
    slabCoefs.setName(component->name);

    // And the new coefficients and write the new HDF5
    slabCoefs.clear();
    slabCoefs.add(cur);
    slabCoefs.WriteH5Coefs(file);
  }
}

