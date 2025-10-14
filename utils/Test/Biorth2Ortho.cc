#include <Biorth2Ortho.H>

Biorth2Ortho::Biorth2Ortho(std::shared_ptr<AxiSymBiorth> b,
			   int Lmax, int Nmax, int Ngrid,
			   double Rmin, double Rmax, double scale, bool wght)
{
  biorth = b;

  nmax   = Nmax;
  lmax   = Lmax;
  ngrid  = Ngrid;
  rmin   = Rmin;
  rmax   = Rmax;
  scl    = scale;
  weight = wght;

  int flag;
  MPI_Initialized(&flag);
  // If flag is true (non-zero), MPI has been initialized.
  if (flag) {
    MPI_Comm_rank(MPI_COMM_WORLD, &myid);
    MPI_Comm_size(MPI_COMM_WORLD, &numprocs);
  }
}


Biorth2Ortho::~Biorth2Ortho(void)
{
}

void Biorth2Ortho::generate()
{
  // Progress bar
  //
  std::shared_ptr<progress::progress_display> progress;

  // Only root process instantiates progress bar
  //
  if (myid==0 and prog_bar) {
    std::cout << std::endl
	      << "Computing new basis functions and transformation matrices"
	      << std::endl;
    progress = std::make_shared<progress::progress_display>(lmax+1);
  }

  // Resize transformation matrices
  //
  Trho.resize(lmax+1);
  Tphi.resize(lmax+1);

  // Resize basis grids
  //
  Wtbl.resize(lmax+1);

  xmin = biorth->r_to_rb(rmin);
  xmax = biorth->r_to_rb(rmax);
  dx   = (xmax - xmin)/(ngrid-1);

  for (int l=0; l<=lmax; l++) {

    // Resize transformation matrices
    //
    Trho[l].resize(nmax, nmax); Trho[l].setZero();
    Tphi[l].resize(nmax, nmax); Tphi[l].setZero();
    Wtbl[l].resize(nmax, ngrid);

    if (l % numprocs == myid) {

      // Copy initial basis functions to grid
      //
      for (int n=0; n<nmax; n++) {
	
	for (int i=0; i<ngrid; i++) {
	  Real x = xmin + dx*i;
	  if (laguerre) {
	    double r = biorth->rb_to_r(x);
	    Wtbl[l](n, i) = std::assoc_laguerre(n, 2, r/scl)*exp(-0.5*r/scl);
	  }
	  else if (alternate) {
	    if (n % 2 == 0)
	      Wtbl[l](n, i) = biorth->dens(n/2+1, l, x);
	    else
	      Wtbl[l](n, i) = biorth->potl((n-1)/2+1, l, x);
	  }
	  else if (sum) {
	    Wtbl[l](n, i) = 0.5*(biorth->dens(n, l, x) +
				 biorth->potl(n, l, x) );
	  }
	  else {
	    if (stype == ScalarType::potential)
	      Wtbl[l](n, i) = biorth->potl(n, l, x);
	    else
	      Wtbl[l](n, i) = biorth->dens(n, l, x);
	    if (weight)
	      Wtbl[l](n, i) /= biorth->dens(0, 0, x);
	  }
	}
      }

      // Compute new basis
      //
      for (int n=0; n<nmax; n++) {
	  
	// Inner product on grid
	//
	auto ip = [&](int j, int k) -> Real {
	  Real ans = 0.0;
	  for (int i=0; i<ngrid; i++) {
	    Real x = std::min(xmin + dx*i, xmax);
	    Real r = biorth->rb_to_r(x);
	    Real J = 1.0/biorth->d_r_to_rb(x);
	    Real w = i==0 or i==ngrid-1 ? 0.5*dx : dx;
	    if (weight) w *= biorth->dens(1, 0, x);

	    ans += Wtbl[l](j, i)*Wtbl[l](k, i)*r*r*J*w;
	  }
	    
	  return ans;
	};
	  
	// Inner product between grid and biorthogonal function
	//
	auto scalar = [&](int j, int k) -> Real {
	  Real ans = 0.0;
	  for (int i=0; i<ngrid; i++) {
	    double x = std::min(xmin + dx*i, xmax);
	    double r = biorth->rb_to_r(x);
	    Real J = 1.0/biorth->d_r_to_rb(x);
	    Real w = i==0 or i==ngrid-1 ? 0.5*dx : dx;
	    if (weight) w *= biorth->dens(0, 0, x);

	    Real f;
	    if (laguerre) {
	      f = std::assoc_laguerre(k, 2, r/scl)*exp(-0.5*r/scl);
	    }
	    else if (alternate) {
	      if (k % 2 == 0)
		f = biorth->dens(k/2, l, x);
	      else
		f = biorth->potl((k-1)/2, l, x);
	    } else if (sum) {
	      f = 0.5*(biorth->dens(k, l, x) +
		       biorth->potl(k+1, l, x) );
	    } else {
	      if (stype == ScalarType::potential)
		f = biorth->potl(k, l, x);
	      else
		f = biorth->dens(k, l, x);
	    }
	    
	    ans += Wtbl[l](j, i)*f*r*r*J*w;
	  }

	  return ans;
	};
	  
	// Classical Gram-Schmidt
	//
	if (classic) {
	  for (int j=0; j<n; j++)
	    Wtbl[l].row(n) -= scalar(j, n) * Wtbl[l].row(j)/ip(j, j);
	}
	// Modified Gram-Schmidt
	//
	else {
	  Real norm2 = ip(n, n);
	  for (int j=n+1; j<nmax; j++) {
	    Wtbl[l].row(j) -= scalar(n, j) * Wtbl[l].row(n)/norm2;
	  }
	}

	// Normalize
	//
	Real norm = std::sqrt(static_cast<double>(ip(n, n)));
	if (norm > 0.0) Wtbl[l].row(n) /= norm;
      }

      // Compute transformation matrix
      //
      for (int n=0; n<nmax; n++) {

	// Inner product between grid and biorthogonal function
	//
	auto scalarProd = [&](ScalarType type, int j) -> Real {
	  Real ans = 0.0;
	  for (int i=0; i<ngrid; i++) {
	    Real x = std::min(xmin + dx*i, xmax);
	    Real r = biorth->rb_to_r(x);
	    Real J = 1.0/biorth->d_r_to_rb(x);
	    Real w = i==0 or i==ngrid-1 ? 0.5*dx : dx;
	    if (weight) w *= biorth->dens(0, 0, x);
	    Real f;
	    if (type == ScalarType::potential)
	      f = biorth->potl(n, l, x);
	    else
	      f = biorth->dens(n, l, x);
	    
	    ans += Wtbl[l](j, i)*f*r*r*J*w;
	  }

	  return ans;
	};
	  
	// Inner product on grid
	//
	for (int j=0; j<nmax; j++) {
	  Trho[l](n, j) = scalarProd(ScalarType::density,   j);
	  Tphi[l](n, j) = scalarProd(ScalarType::potential, j);
	}
      }
      // End of outer n loop
    }
    // End of parallel section

    // Update progress bar
    //
    if (progress) ++(*progress);
  }
  // End of loop over l

  // Broadcast results to all processes
  //
  for (int l=0; l<=lmax; l++) {
    int id = l % numprocs;
    MPI_Bcast(Trho[l].data(), Trho[l].size(), MPI_LONG_DOUBLE, id, MPI_COMM_WORLD);
    MPI_Bcast(Tphi[l].data(), Tphi[l].size(), MPI_LONG_DOUBLE, id, MPI_COMM_WORLD);
    MPI_Bcast(Wtbl[l].data(), Wtbl[l].size(), MPI_LONG_DOUBLE, id, MPI_COMM_WORLD);
    MPI_Barrier(MPI_COMM_WORLD);
  }
  // End of loop over l
}

void Biorth2Ortho::output(const std::string& PREFIX)
{
  if (myid) return;

  // Progress bar
  //
  std::shared_ptr<progress::progress_display> progress;

  if (prog_bar) {
    std::cout << std::endl
	      << "Making general output file with orthogonality tests"
	      << std::endl;
    progress = std::make_shared<progress::progress_display>(lmax+1);
  }

  std::ostringstream fname;
  fname << PREFIX << ".out";

  std::ofstream ofs(fname.str());

  if (ofs) {
    ofs << "# Biorthogonal expansion transformation" << std::endl;
    ofs << "# lmax: " << lmax << " nmax: " << nmax << std::endl;
    ofs << "# rmin: " << rmin << " rmax: " << rmax << std::endl;
    ofs << "# ngrid: " << ngrid << std::endl;

    MatrixXld tmp;

    // Output transformation matrices
    //
    for (int l=0; l<=lmax; l++) {
      Eigen::MatrixXd t1 = Trho[l].cast<double>();
      Eigen::MatrixXd t2 = Tphi[l].cast<double>();
      Eigen::MatrixXd i1 = t1.inverse();
      Eigen::MatrixXd i2 = t2.inverse();

      ofs << "# Trho l=" << l << " norm=" << t1.norm() << std::endl;
      ofs << t1 << std::endl;
      ofs << "# Tphi l=" << l << " norm=" << t2.norm() << std::endl;
      ofs << t2 << std::endl;
      ofs << "# Trho^{-1} l=" << l << " norm=" << i1.norm() << std::endl;
      ofs << i1 << std::endl;
      ofs << "# Tphi^{-1} l=" << l << " norm=" << i2.norm() << std::endl;
      ofs << i2 << std::endl;
      ofs << "# Trho*Tphi^{T} l=" << l << std::endl;
      ofs << t1*t2.transpose() << std::endl;
      ofs << "# Tphi*Trho^{T} l=" << l << std::endl;
      ofs << t2*t1.transpose() << std::endl;
      ofs << "# Trho^{-1}*Trho l=" << l << std::endl;
      ofs << i1*t1 << std::endl;
      ofs << "# Tphi^{-1}*Tphi l=" << l << std::endl;
      ofs << i2*t2 << std::endl;
      ofs << "# Trho^{T,-1}*Tphi^{-1} l=" << l << std::endl;
      ofs << i1.transpose()*i2 << std::endl;
      ofs << "# Tphi^{T,-1}*Trho^{-1} l=" << l << std::endl;
      ofs << i2.transpose()*i1 << std::endl;
      ofs << "# (Trho*Tphi^{T})^{-1} l=" << l << std::endl;
      ofs << (-t1*t2.transpose()).inverse() << std::endl;
      ofs << "# (Tphi*Trho^{T})^{-1} l=" << l << std::endl;
      ofs << (-t2*t1.transpose()).inverse() << std::endl;

      // Inner product on grid
      //
      auto ip = [&](int j, int k) -> Real {
	Real ans = 0.0;
	for (int i=0; i<ngrid; i++) {
	  Real x = xmin + dx*i;
	  Real r = biorth->rb_to_r(x);
	  Real J = 1.0/biorth->d_r_to_rb(x);
	  Real w = i==0 or i==ngrid-1 ? 0.5*dx : dx;
	  if (weight) w *= biorth->dens(0, 0, x);
	  ans += Wtbl[l](j, i)*Wtbl[l](k, i)*r*r*J*w;
	}
	
	return ans;
      };

      ofs << "# Inner product check l=" << l << std::endl;
      for (int j=0; j<nmax; j++) {
	for (int k=0; k<nmax; k++)
	  ofs << std::setw(16) << static_cast<double>(ip(j, k));
	ofs << std::endl;
      }

      ofs << "# Biorth inner product check l=" << l << std::endl;

      auto scalarTest = [&](int j, int k) -> Real {
	Real ans = 0.0;
	for (int i=0; i<ngrid; i++) {
	  Real x = std::min(xmin + dx*i, xmax);
	  Real r = biorth->rb_to_r(x);
	  Real J = 1.0/biorth->d_r_to_rb(x);
	  Real w = i==0 or i==ngrid-1 ? 0.5*dx : dx;
	  if (weight) w *= biorth->dens(0, 0, x);
	  Real f = biorth->dens(j, l, x) * biorth->potl(k, l, x);
	  ans += f*w*r*r*J;
	}							      
	return ans;
      };
	  
      for (int j=0; j<nmax; j++) {
	for (int k=0; k<nmax; k++)
	  ofs << std::setw(16) << static_cast<double>(scalarTest(j, k));
	ofs << std::endl;
      }

      if (progress) ++(*progress);
    }
  }
  else {
    std::ostringstream serr;
    serr << "Error opening output file <" << fname.str() << ">";
    throw std::runtime_error(serr.str());
  }

  // Progress bar
  if (prog_bar) {
    std::cout << std::endl
	      << "Writing basis functions and diagnostic differences"
	      << std::endl;
    progress = std::make_shared<progress::progress_display>(lmax+1);
  }

  // Output new basis functions
  //
  for (int l=0; l<=lmax; l++) {
    
    std::ostringstream ostr, dstr, pstr, fdel;

    ostr << PREFIX << ".basis." << l;
    dstr << PREFIX << ".den." << l;
    pstr << PREFIX << ".pot." << l;
    fdel << PREFIX << ".del." << l;
    
    std::ofstream ofsW(ostr.str());
    if (!ofsW) {
      std::ostringstream serr;
      serr << "Biorth2Ortho::output: error opening output file <"
	   << ostr.str() << ">";
      throw std::runtime_error(serr.str());
    }

    std::ofstream ofsD(dstr.str());
    if (!ofsD) {
      std::ostringstream serr;
      serr << "Biorth2Ortho::output: error opening output file <" << dstr.str() << ">";
      throw std::runtime_error(serr.str());
    }

    std::ofstream ofsP(pstr.str());
    if (!ofsP) {
      std::ostringstream serr;
      serr << "Biorth2Ortho::output: error opening output file <" << pstr.str() << ">";
      throw std::runtime_error(serr.str());
    }

    std::ofstream delF(fdel.str());
    if (!delF) {
      std::ostringstream serr;
      serr << "Biorth2Ortho::output: error opening output file <" << fdel.str() << ">";
      throw std::runtime_error(serr.str());
    }

    ofsW << "# l=" << l << std::endl;
    ofsW << "# x r fcts..." << std::endl;

    ofsD << "# l=" << l << std::endl;
    ofsD << "# x r fcts..." << std::endl;

    ofsP << "# l=" << l << std::endl;
    ofsP << "# x r fcts..." << std::endl;
    
    for (int j=0; j<ngrid; j++) {
      
      double x = xmin + dx*j;
      double r = biorth->rb_to_r(x);
      
      ofsW << std::setw(16) << x << std::setw(16) << r;
      ofsD << std::setw(16) << x << std::setw(16) << r;
      ofsP << std::setw(16) << x << std::setw(16) << r;
      
      for (int n=0; n<nmax; n++) {
	ofsW << std::setw(16) << static_cast<double>(Wtbl[l](n, j));
	Real D = 0.0, P = 0.0;	
	for (int m=0; m<nmax; m++) {
	  D += Trho[l](n, m)*Wtbl[l](m, j);
	  P += Tphi[l](n, m)*Wtbl[l](m, j);
	}
	ofsD << std::setw(16) << static_cast<double>(D) << std::setw(16) << biorth->dens(n, l, x);
	ofsP << std::setw(16) << static_cast<double>(P) << std::setw(16) << biorth->potl(n, l, x);
      }
      ofsW << std::endl;
      ofsD << std::endl;
      ofsP << std::endl;
    }
    
    auto drho = Trho[l] * Wtbl[l];
    auto dphi = Tphi[l] * Wtbl[l];

    for (int n=0; n<nmax; n++) {
      VectorXld D(ngrid), P(ngrid), W(ngrid);
      for (int j=0; j<ngrid; j++) {
	Real x = xmin + dx*j;
	Real r = biorth->rb_to_r(x);
	Real J = 1.0/biorth->d_r_to_rb(x);
	D(j) = biorth->dens(n, l, x);
	P(j) = biorth->potl(n, l, x);
	W(j) = r*r*J;
      }
      
      // Differences
      auto dataD = drho.transpose().col(n) - D;
      auto dataP = dphi.transpose().col(n) - P;

      // Weighted inner products
      Real facD = D.transpose() * D.cwiseProduct(W);
      Real facP = P.transpose() * P.cwiseProduct(W);
      
      // Vector norms
      Real difD =  dataD.transpose()*dataD.cwiseProduct(W);
      Real difP =  dataP.transpose()*dataP.cwiseProduct(W);
      
      delF  << std::setw(8) << n
	    << std::setw(16) << static_cast<double>(difD)
	    << std::setw(16) << static_cast<double>(facD)
	    << std::setw(16) << static_cast<double>(difP)
	    << std::setw(16) << static_cast<double>(facP)
	    << std::endl;
    }
    // Loop over n
    
    if (progress) ++(*progress);
  }
  // Loop over l
}


const double Biorth2Ortho::basis(int l, int n, double r)
{
  if (l<0 || l>lmax || n<0 || n>=nmax) {
    std::ostringstream serr;
    serr << "Biorth2Ortho::basis: l=" << l << " n=" << n
	 << " out of range!\n";
    throw std::runtime_error(serr.str());
  }

  // Get scaled radius
  //
  double x = biorth->r_to_rb(r);

  // Linear interpolation on grid
  //
  int i = (int)((x - xmin)/dx);

  // Off the grid
  if (i<0 || i>=ngrid-1) return 0.0;

  double f1 = Wtbl[l](n, i);
  double f2 = Wtbl[l](n, i+1);

  double x1 = xmin + dx*i;
  double x2 = xmin + dx*(i+1);

  return f1 + (f2 - f1)*(x - x1)/(x2 - x1);
}

void Biorth2Ortho::writeH5(const std::string& h5file)
{
  // Only root node makes the cachefile
  //
  if (myid==0) {

    try {
      // Create a new hdf5 file or overwrite an existing file
      //
      HighFive::File file(h5file, HighFive::File::Overwrite);
    

      // Algorithm switches as unsigned ints
      //
      unsigned int flg  = laguerre  ? 1 : 0;

      file.createAttribute<unsigned>("laguerre",  HighFive::DataSpace::From(flg)).write(flg);
      
      flg  = alternate  ? 1 : 0;
      file.createAttribute<unsigned>("alternate",  HighFive::DataSpace::From(flg)).write(flg);
      
      flg  = classic  ? 1 : 0;
      file.createAttribute<unsigned>("classic",  HighFive::DataSpace::From(flg)).write(flg);
      
      flg  = density  ? 1 : 0;
      file.createAttribute<unsigned>("density",  HighFive::DataSpace::From(flg)).write(flg);
    
      flg = weight ? 1 : 0;
      file.createAttribute<unsigned>("weight",  HighFive::DataSpace::From(flg)).write(flg);

      // Write parameters
      //
      file.createAttribute<int>    ("nmax",  HighFive::DataSpace::From(nmax)).write(nmax);
      file.createAttribute<int>    ("lmax",  HighFive::DataSpace::From(lmax)).write(lmax);
      file.createAttribute<int>    ("ngrid", HighFive::DataSpace::From(ngrid)).write(ngrid);

      file.createAttribute<double> ("rmin",  HighFive::DataSpace::From(rmin)).write(rmin);
      file.createAttribute<double> ("rmax",  HighFive::DataSpace::From(rmax)).write(rmax);
      file.createAttribute<double> ("scl",   HighFive::DataSpace::From(scl)).write(scl);

      file.createAttribute<double> ("xmin",  HighFive::DataSpace::From(xmin)).write(xmin);
      file.createAttribute<double> ("xmax",  HighFive::DataSpace::From(xmax)).write(xmax);
      file.createAttribute<double> ("dx",    HighFive::DataSpace::From(dx)).write(dx);


      // Save basis and transformation matrices
      //
      file.createDataSet("Trho", Trho);
      file.createDataSet("Tphi", Tphi);
      file.createDataSet("Wtbl", Wtbl);

      
    } catch (HighFive::Exception& err) {
      std::cerr << err.what() << std::endl;
    }
    
    std::cout << std::endl
	      << "---- RespMatFT::writeH5: "
	      << "wrote <" << h5file << ">" << std::endl;
  }

}


template <typename T>
bool checkValB2O(HighFive::File& h5file, const T& value, const std::string& name)
{
  if (h5file.hasAttribute(name)) {
    T v;
    h5file.getAttribute(name).read(v);
    if (value == v) return true;
    if (myid==0)
      std::cout << "---- readH5: "
		<< "parameter " << name << ": wanted " << value
		<< " found " << v << std::endl;
  }
  return false;
};
    


template <typename T>
T readValB2O(HighFive::File& h5file, const std::string& name)
{
  T v;
  h5file.getAttribute(name).read(v);
  return v;
};
    
void Biorth2Ortho::readH5(const std::string& file, bool check)
{
  try {
    // Silence the HDF5 error stack
    //
    HighFive::SilenceHDF5 quiet;
      
    // Try opening the file as HDF5
    //
    HighFive::File h5file(file, HighFive::File::ReadOnly);
    
    // Read algorithm switches
    //
    auto checkFlag = [&](HighFive::File& h5file, const std::string& name, bool flag) -> bool
    {
      unsigned int f = flag ? 1 : 0, v;
      h5file.getAttribute(name).read(v);
      if (f == v) return true;
      if (myid==0) {
	bool iflg = v ? true : false;
	std::cout << "---- Biorth2Ortho::readH5: "
		  << "parameter " << name << ": wanted " << flag
		  << " found " << std::boolalpha << iflg << std::endl;
      }
      return false;
    };
    
    auto readFlag = [](HighFive::File& h5file, const std::string& name)
    {
      unsigned int v;
      h5file.getAttribute(name).read(v);
      return v ? true : false;
    };
    
    // Check switches and parameters
    //
    if (check) {
      if (not checkFlag(h5file, "laguerre",  laguerre )) throw MismatchException("Biorth2Ortho::readH5 parameter mismatch");
      if (not checkFlag(h5file, "alternate", alternate)) throw MismatchException("Biorth2Ortho::readH5 parameter mismatch");
      if (not checkFlag(h5file, "sum",       sum      )) throw MismatchException("Biorth2Ortho::readH5 parameter mismatch");
      if (not checkFlag(h5file, "classic",   classic  )) throw MismatchException("Biorth2Ortho::readH5 parameter mismatch");
      if (not checkFlag(h5file, "density",   density  )) throw MismatchException("Biorth2Ortho::readH5 parameter mismatch");
      if (not checkFlag(h5file, "weight",    weight   )) throw MismatchException("Biorth2Ortho::readH5 parameter mismatch");
 
      if (not checkValB2O(h5file, nmax,     "nmax" ))    throw MismatchException("Biorth2Ortho::readH5 parameter mismatch");
      if (not checkValB2O(h5file, lmax,     "lmax" ))    throw MismatchException("Biorth2Ortho::readH5 parameter mismatch");
      if (not checkValB2O(h5file, ngrid,    "ngrid"))    throw MismatchException("Biorth2Ortho::readH5 parameter mismatch");

      if (not checkValB2O(h5file, rmin,     "rmin" ))    throw MismatchException("Biorth2Ortho::readH5 parameter mismatch");
      if (not checkValB2O(h5file, rmax,     "rmax" ))    throw MismatchException("Biorth2Ortho::readH5 parameter mismatch");
      if (not checkValB2O(h5file, scl,      "scl"  ))    throw MismatchException("Biorth2Ortho::readH5 parameter mismatch");

      if (not checkValB2O(h5file, xmin,     "xmin" ))    throw MismatchException("Biorth2Ortho::readH5 parameter mismatch");
      if (not checkValB2O(h5file, xmax,     "xmax" ))    throw MismatchException("Biorth2Ortho::readH5 parameter mismatch");
      if (not checkValB2O(h5file, dx,       "dx"   ))    throw MismatchException("Biorth2Ortho::readH5 parameter mismatch");
    }
      
    // Flag read
    //
    laguerre  = readFlag(h5file, "laguerre" );
    alternate = readFlag(h5file, "alternate");
    sum       = readFlag(h5file, "sum"      );
    classic   = readFlag(h5file, "classic"  );
    density   = readFlag(h5file, "density"  );
    weight    = readFlag(h5file, "weight"   );

    // Parameter read
    //
    nmax  = readValB2O<int>   (h5file, "nmax" );
    lmax  = readValB2O<int>   (h5file, "lmax" );
    ngrid = readValB2O<int>   (h5file, "ngrid");
    rmin  = readValB2O<double>(h5file, "rmin");
    rmax  = readValB2O<double>(h5file, "rmax");
    scl   = readValB2O<double>(h5file, "scl" );
    xmin  = readValB2O<double>(h5file, "xmin");
    xmin  = readValB2O<double>(h5file, "xmax");
    dx    = readValB2O<double>(h5file, "dx"  );

    // Basis and transformation matrices
    //
    h5file.getDataSet("Trho").read(Trho);
    h5file.getDataSet("Tphi").read(Tphi);
    h5file.getDataSet("Wtbl").read(Wtbl);

    if (myid==0)
      std::cout << "---- readH5: "
		<< "successfully read <" << file << ">" << std::endl;

  } catch (HighFive::Exception& err) {
    if (myid==0) {
      std::cerr << "---- readH5: "
		<< err.what() << std::endl
		<< "---- readH5: "
		<< "error reading previous orthogonal transforms" << std::endl;
    }
    throw HDF5Exception("---- Biorth2Ortho::readH5: read failure");
  }
  
}


