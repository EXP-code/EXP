#include <algorithm>
#include <iomanip>
#include <sstream>
#include <cctype>
#include <string>

#include <FieldGenerator.H>
#include <DataGrid.H>
#include <localmpi.H>

// Verbose output for checking stack usage
#ifdef DEBUG
#include <sys/time.h>
#include <sys/resource.h>
#endif

namespace Field
{
  
  FieldGenerator::FieldGenerator(const std::vector<double> &time,
				 const std::vector<double> &pmin,
				 const std::vector<double> &pmax,
				 const std::vector<int>    &grid
				 ) :
    times(time), pmin(pmin), pmax(pmax), grid(grid)
  {
    // Check whether MPI is initialized
    //
    int flag;
    MPI_Initialized(&flag);
    if (flag) use_mpi = true;
    else      use_mpi = false;

    
    // Fall back sanity (works for me but this needs to be fixed
    // generally)
    //
    if (use_mpi) {
      MPI_Comm_size(MPI_COMM_WORLD, &numprocs);
      MPI_Comm_rank(MPI_COMM_WORLD, &myid);
    }
  }
  
  void FieldGenerator::check_times(CoefClasses::CoefsPtr coefs)
  {
    std::vector<double> ctimes = coefs->Times();
    std::sort(ctimes.begin(), ctimes.end());
    
    for (auto t : times) {
      if (std::find(ctimes.begin(), ctimes.end(), t) == ctimes.end()) {
	std::ostringstream sout;
	sout << "FieldGenerator: requested time <" << t << "> "
	     << "not in DB" << std::endl;
	throw std::runtime_error(sout.str());
      }
    }
  }
  
  std::map<double, std::map<std::string, Eigen::VectorXf>>
  FieldGenerator::lines
  (BasisClasses::BasisPtr basis, CoefClasses::CoefsPtr coefs,
   std::vector<double> beg, std::vector<double> end, int num)
  {
    // Check
    //
    if (beg.size()!=3 or end.size()!=3) {
      throw std::runtime_error("FieldGenerator::lines: vectors beg and end must have rank 3");
    }

    if (num<1) {
      throw std::runtime_error("FieldGenerator::lines: number of evaluation points must be > 0");
    }

    check_times(coefs);

    std::map<double, std::map<std::string, Eigen::VectorXf>> ret;

    // Now get the coordinate type
    //
    auto ctype = basis->coordinates;

    // Coordinate and field labels
    //
    std::vector<std::string> coords {"x", "y", "z", "arc"};
    auto labels = basis->getFieldLabels(ctype);

    // Initialize the frame map
    //
    std::map<std::string, Eigen::VectorXf> frame;
    for (auto label : coords) frame[label] = Eigen::VectorXf::Zero(num);
    for (auto label : labels) frame[label] = Eigen::VectorXf::Zero(num);

    // Compute the probe length
    //
    std::vector<double> dd(3);
    for (int k=0; k<3; k++) dd[k] = (end[k] - beg[k])/num;
    double dlen = sqrt(dd[0]*dd[0] + dd[1]*dd[1] + dd[2]*dd[2]);

    for (int icnt=0; icnt<times.size(); icnt++) {

      if (icnt % numprocs == myid) {
      
	double T = times[icnt];

	if (not coefs->getCoefStruct(T)) {
	  std::cout << "Could not find time=" << T << ", continuing"
		    << std::endl;
	  continue;
	}

	basis->set_coefs(coefs->getCoefStruct(T));

	double r, phi, costh, R;
	double p0, p1, d0, d1, f1, f2, f3;
      
#pragma omp parallel for
	for (int ncnt=0; ncnt<num; ncnt++) {

	  double x = beg[0] + dd[0]*ncnt;
	  double y = beg[1] + dd[1]*ncnt;
	  double z = beg[2] + dd[2]*ncnt;
	  
	  std::vector<double> v;

	  if (ctype == BasisClasses::Basis::Coord::Spherical) {
	    r     = sqrt(x*x + y*y + z*z) + 1.0e-18;
	    costh = z/r;
	    phi   = atan2(y, x);
	    v = (*basis)(r, costh, phi, ctype);
	  } else if (ctype == BasisClasses::Basis::Coord::Cylindrical) {
	    R     = sqrt(x*x + y*y) + 1.0e-18;
	    phi   = atan2(y, x);
	    v = (*basis)(R, z, phi, ctype);
	  } else {		// A default
	    ctype = BasisClasses::Basis::Coord::Cartesian;
	    v = (*basis)(x, y, z, ctype);
	  }
	  
	  frame["x"      ](ncnt) = x;
	  frame["y"      ](ncnt) = y;
	  frame["z"      ](ncnt) = z;
	  frame["arc"    ](ncnt) = dlen*ncnt;

	  for (int n=0; n<labels.size(); n++) frame[labels[n]](ncnt) = v[n];
	}

	ret[T] = frame;
      }
    }
    
    if (use_mpi) {

      for (int n=0; n<numprocs; n++) {
	
	if (myid==n) {
	  // Send the number of time frames in this node
	  int num = ret.size();
	  MPI_Send(&num, 1, MPI_INT, 0, 200, MPI_COMM_WORLD);

	  // Iterate through the time frames
	  for (auto & F : ret) {
	    MPI_Send(&F.first, 1, MPI_DOUBLE, 0, 201, MPI_COMM_WORLD);

	    int nf = F.second.size(); // Number of fields in this frame
	    MPI_Send(&num, 1, MPI_INT, 0, 202, MPI_COMM_WORLD);

	    // Iterate through the fields
	    for (auto & f : F.second) {
	      // Send tag
	      MPI_Send(f.first.c_str(), f.first.length(), MPI_CHAR, 0, 203,
		       MPI_COMM_WORLD);
	      // Send field data
	      MPI_Send(f.second.data(), f.second.size(),
		       MPI_FLOAT, 0, 204, MPI_COMM_WORLD);
	    }

	  }
	}

	if (myid==0) {
	  std::vector<char> bf(9); // Char buffer for field label
	  MPI_Status status;	   // For MPI_Probe
	  int num, nf, l;
	  double T;

	  // Get the number of frames from node n
	  MPI_Recv(&num, 1, MPI_INT, n, 200, MPI_COMM_WORLD,
		   MPI_STATUS_IGNORE);

	  // Iterate through the frames
	  for (int j=0; j<num; j++) {
	    // Get the time for this frame
	    MPI_Recv(&T, 1, MPI_DOUBLE, n, 201, MPI_COMM_WORLD, &status);

	    // Get the number of fields in this frame
	    MPI_Recv(&nf, 1, MPI_INT, n, 202, MPI_COMM_WORLD, &status);

	    // Iterate through the fields
	    for (int k=0; k<nf; k++) {
	      MPI_Probe(n, 203, MPI_COMM_WORLD, &status);
	      MPI_Get_count(&status, MPI_CHAR, &l);
	      MPI_Recv(bf.data(), l, MPI_CHAR, n, 203, MPI_COMM_WORLD, &status);

	      // Sanity check on the field name
	      //
	      std::string s(bf.data(), l);
	      if (frame.find(s) == frame.end()) {
		std::cerr << "Error finding <" << s << "> in field map"
			  << std::endl;
	      }

	      // Get the data
	      //
	      MPI_Recv(frame[s].data(), frame[s].size(), MPI_FLOAT, n, 204,
		       MPI_COMM_WORLD, &status);
	    }
	    
	    ret[T] = frame;
	  }
	}
      }
    }

    return ret;
  }

  void FieldGenerator::file_lines(BasisClasses::BasisPtr basis,
				  CoefClasses::CoefsPtr  coefs,
				  std::vector<double>    beg,
				  std::vector<double>    end,
				  int                    num,
				  const std::string      prefix,
				  const std::string      outdir)
  {
    auto db = lines(basis, coefs, beg, end, num);

    if (myid==0) {

      int icnt = 0;

      for (auto & frame : db) {

	// Create the file name and open
	//
	std::ostringstream sout;
	sout << outdir << "/" << prefix << "_probe_" << icnt << ".txt";

	std::ofstream out(sout.str());

	if (out) {

	  // Write file header
	  //
	  out << "# T=" << frame.first << std::endl;

	  int cnt = 0;
	  for (auto u : frame.second) {
	    std::string label = u.first + " ";
	    if (cnt++==0) out << "#" << std::setw(15) << std::right << label;
	    else          out << std::setw(16) << std::right << label;
	  }
	  out << std::endl;
	  for (int n=0; n<cnt; n++) {
	    std::ostringstream sout; sout << "[" << n+1 << "] ";
	    if (n==0) out << "#" << std::setw(15) << std::right << sout.str();
	    else out << std::setw(16) << std::right << sout.str();
	  }
	  out << std::endl;
	  for (int n=0; n<cnt; n++) {
	    if (n==0) out << "#" << std::setw(15) << std::right << std::string(10, '-');
	    else out << std::setw(16) << std::right << std::string(10, '-');
	  }
	  out << std::endl;


	  for (int i=0; i<num; i++) {
	    for (auto & v : frame.second) {
	      out << std::setw(16) << v.second(i);
	    }
	    out << std::endl;
	  }
	} else {
	  throw std::runtime_error("FieldGenerator::file_lines: couldn't open <" + sout.str() + ">");
	}

	icnt++;
      }
    }
  }
  
  
  std::map<double, std::map<std::string, Eigen::MatrixXf>>
  FieldGenerator::slices(BasisClasses::BasisPtr basis,
			 CoefClasses::CoefsPtr coefs)
  {
    // Set midplane evaluation parameters
    //
    basis->setMidplane(midplane);
    basis->setColumnHeight(colheight);

    // Check
    //
    check_times(coefs);

    std::map<double, std::map<std::string, Eigen::MatrixXf>> ret;

    // Now get the desired coordinate type
    //
    auto ctype = basis->coordinates;

    // Field labels (force field labels added below)
    //
    auto labels = basis->getFieldLabels(ctype);

    // Find the first two non-zero indices
    //
    int i1=-1, i2=-1, i3=-1;
    std::vector<double> pos, del;
    for (size_t i=0; i<grid.size(); i++) {
      pos.push_back(pmin[i]);
      if (grid[i]>0) {
	if (i1<0) {
	  i1 = i;
	  del.push_back((pmax[i1] - pmin[i1])/std::max<int>(grid[i1]-1, 1));
	}
	else if (i2<0) {
	  i2 = i;
	  del.push_back((pmax[i2] - pmin[i2])/std::max<int>(grid[i2]-1, 1));
	}
      } else {
	del.push_back(0.0);
	i3 = i;
      }
    } 

    if (i1<0 or i2<0 or i3<0)
      throw std::runtime_error("FieldGenerator::slices: bad grid specification");

    int ncnt = 0;		// Process counter for MPI

    std::map<std::string, Eigen::MatrixXf> frame;
    for (auto label : labels) {
      frame[label].resize(grid[i1], grid[i2]);
    }	

    for (auto T : times) {

      if (ncnt++ % numprocs > 0) continue;

      if (not coefs->getCoefStruct(T)) {
	std::cout << "Could not find time=" << T << ", continuing" << std::endl;
	continue;
      }

      basis->set_coefs(coefs->getCoefStruct(T));

      int totpix = grid[i1] * grid[i2];

#pragma omp parallel for
      for (int k=0; k<totpix; k++) {

	// Create the pair of indices from the pixel number
	//
	int i = k/grid[i2];
	int j = k - i*grid[i2];

	// Compute the coordinates from the indices
	//
	std::vector<double> pp(pos);

	pp[i1] = pmin[i1] + del[i1]*i;
	pp[i2] = pmin[i2] + del[i2]*j;

	// Cartesian to spherical for all_eval
	//
	double x = pp[0];
	double y = pp[1];
	double z = pp[2];

	// Coordinate values
	double r, costh, phi, R;

	// Return values
	double p0, p1, d0, d1, f1, f2, f3;
	std::vector<double> v;

	if (ctype == BasisClasses::Basis::Coord::Spherical) {
	  r     = sqrt(x*x + y*y + z*z) + 1.0e-18;
	  costh = z/r;
	  phi   = atan2(y, x);
	  v = (*basis)(r, costh, phi, ctype);
	} else if (ctype == BasisClasses::Basis::Coord::Cylindrical) {
	  R     = sqrt(x*x + y*y) + 1.0e-18;
	  phi   = atan2(y, x);
	  v = (*basis)(R, z, phi, ctype);
	} else {
	  v = (*basis)(x, y, z, BasisClasses::Basis::Coord::Cartesian);
	}
	
	// Pack the frame structure
	//
	for (int n=0; n<labels.size(); n++)
	  frame[labels[n]](i, j) = v[n];
      }

      ret[T] = frame;
    }

    if (use_mpi) {
      
      std::vector<char> bf(9);

      for (int n=1; n<numprocs; n++) {

	if (myid==n) {
	  int sz = ret.size();
	  MPI_Send(&sz, 1, MPI_INT, 0, 102, MPI_COMM_WORLD);
	  for (auto & v : ret) {
	    MPI_Send(&v.first, 1, MPI_DOUBLE, 0, 103, MPI_COMM_WORLD);
	    int fsz = v.second.size();
	    MPI_Send(&fsz, 1, MPI_INT, 0, 104, MPI_COMM_WORLD);
	    
	    for (auto & f : v.second) {
	      MPI_Send(f.first.c_str(), f.first.length(), MPI_CHAR, 0, 105,
		       MPI_COMM_WORLD);
	      
	      MPI_Send(f.second.data(), f.second.size(),
		       MPI_FLOAT, 0, 106, MPI_COMM_WORLD);
	    }
	  }
	}

	if (myid==0) {
	  MPI_Status status;
	  int sz, fsz, l;
	  double T;

	  MPI_Recv(&sz, 1, MPI_INT, n, 102, MPI_COMM_WORLD, &status);
	  for (int i=0; i<sz; i++) {
	    MPI_Recv(&T, 1, MPI_DOUBLE, n, 103, MPI_COMM_WORLD, &status);
	    MPI_Recv(&fsz, 1, MPI_INT, n, 104, MPI_COMM_WORLD, &status);
	    
	    for (int j=0; j<fsz; j++) {
	      // Get the field name
	      //
	      MPI_Probe(n, 105, MPI_COMM_WORLD, &status);
	      MPI_Get_count(&status, MPI_CHAR, &l);
	      MPI_Recv(bf.data(), l, MPI_CHAR, n, 105, MPI_COMM_WORLD, &status);
	      std::string s(bf.data(), l);
	      
	      // Sanity check
	      //
	      if (frame.find(s) == frame.end()) {
		std::cerr << "Error finding <" << s << "> in field map"
			  << std::endl;
	      }

	      // Get the data
	      //
	      MPI_Recv(frame[s].data(), frame[s].size(), MPI_FLOAT, n, 106,
		       MPI_COMM_WORLD, &status);
	    }
	    
	    ret[T] = frame;
	  }
	}
      }
    }

    // Toggle off midplane evaluation
    basis->setMidplane(false);

    return ret;
  }
  
  void FieldGenerator::file_slices(BasisClasses::BasisPtr basis,
				   CoefClasses::CoefsPtr  coefs,
				   const std::string      prefix,
				   const std::string      outdir)
  {
    auto db = slices(basis, coefs);

    if (myid==0) {

      // Find the first two non-zero indices
      int i1=-1, i2=-1, i3=-1;
      for (size_t i=0; i<grid.size(); i++) {
	if (grid[i]>0) {
	  if (i1<0) i1 = i;
	  else if (i2<0) i2 = i;
	} else i3 = i;
      }

      int icnt = 0;

      for (auto & frame : db) {

	DataGrid datagrid(grid[i1], grid[i2], 1,
			  pmin[i1], pmax[i1], pmin[i2], pmax[i2], 0, 0);

	std::vector<double> tmp(grid[i1]*grid[i2]);

	for (auto & v : frame.second) {

	  for (int i=0; i<grid[i1]; i++) {
	    for (int j=0; j<grid[i2]; j++) {
	      tmp[j*grid[i1] + i] = v.second(i, j);
	    }
	  }

	  datagrid.Add(tmp, v.first);
	}

	std::ostringstream sout;
	sout << outdir << "/" << prefix << "_surface_" << icnt;
	datagrid.Write(sout.str());
	icnt++;
      }
    }
  }
  
  
  std::map<double, std::map<std::string, Eigen::Tensor<float, 3>>>
  FieldGenerator::volumes(BasisClasses::BasisPtr basis,
			  CoefClasses::CoefsPtr coefs)
  {
    std::map<double, std::map<std::string, Eigen::Tensor<float, 3>>> ret;

    // Now get the desired coordinate type
    //
    auto ctype = basis->coordinates;

    // Field labels (force field labels added below)
    //
    auto labels = basis->getFieldLabels(ctype);

    // Allocate frame storge
    //
    std::map<std::string, Eigen::Tensor<float, 3>> frame;
    for (auto label : labels) {
      frame[label].resize(grid[0], grid[1], grid[2]);
    }	
    
    std::vector<double> del =
      { (pmax[0] - pmin[0])/std::max<int>(grid[0]-1, 1),
	(pmax[1] - pmin[1])/std::max<int>(grid[1]-1, 1),
	(pmax[2] - pmin[2])/std::max<int>(grid[2]-1, 1) };
    
    int ncnt = 0;		// Process counter for MPI

    for (auto T : times) {

      if (ncnt++ % numprocs > 0) continue;

      basis->set_coefs(coefs->getCoefStruct(T));

      int totpix = grid[0] * grid[1] * grid[2];

#pragma omp parallel for
      for (int n=0; n<totpix; n++) {

	// Unpack the index triple by integer division
	//
	int i = n/(grid[1]*grid[2]);
	int j = (n - i*grid[1]*grid[2])/grid[2];
	int k = n - (i*grid[1] + j)*grid[2];

	// Compute the coordinates from the indices
	//
	double x = pmin[0] + del[0]*i;
	double y = pmin[1] + del[1]*j;
	double z = pmin[2] + del[2]*k;
	    
	std::vector<double> v;

	if (ctype == BasisClasses::Basis::Coord::Spherical) {
	  double r     = sqrt(x*x + y*y + z*z) + 1.0e-18;
	  double costh = z/r;
	  double phi   = atan2(y, x);
	  v = (*basis)(r, costh, phi, ctype);
	} else if (ctype == BasisClasses::Basis::Coord::Cylindrical) {
	  double R     = sqrt(x*x + y*y) + 1.0e-18;
	  double phi   = atan2(y, x);
	  v = (*basis)(R, z, phi, ctype);
	} else {
	  ctype = BasisClasses::Basis::Coord::Cartesian;
	  v = (*basis)(x, y, z, ctype);
	}

	// Pack the frame structure
	//
	for (int n=0; n<labels.size(); n++)
	  frame[labels[n]](i, j, k) = v[n];
      }

      ret[T] = frame;

#ifdef DEBUG
      if (myid==0) {
	rusage usage;
	int err = getrusage(RUSAGE_SELF, &usage);
	std::cout << "volumes: T=" << std::setw(8) << std::fixed<< T
		  << " Size=" << std::setw(8) << usage.ru_maxrss/1024/1024
		  << std::endl;
      }
#endif

    }

    if (use_mpi) {

      std::vector<char> bf(9);

      for (int n=1; n<numprocs; n++) {

	if (myid==n) {
	  int sz = ret.size();
	  MPI_Send(&sz, 1, MPI_INT, 0, 102, MPI_COMM_WORLD);
	  for (auto & v : ret) {
	    MPI_Send(&v.first, 1, MPI_DOUBLE, 0, 103, MPI_COMM_WORLD);
	    int fsz = v.second.size();
	    MPI_Send(&fsz, 1, MPI_INT, 0, 104, MPI_COMM_WORLD);
	    
	    for (auto & f : v.second) {
	      MPI_Send(f.first.c_str(), f.first.length(), MPI_CHAR, 0, 105,
		       MPI_COMM_WORLD);
	      
	      MPI_Send(f.second.data(), f.second.size(),
		       MPI_FLOAT, 0, 106, MPI_COMM_WORLD);
	    }
	  }
	}

	if (myid==0) {
	  MPI_Status status;
	  int sz, fsz, l;
	  double T;

	  MPI_Recv(&sz, 1, MPI_INT, n, 102, MPI_COMM_WORLD, &status);
	  for (int i=0; i<sz; i++) {
	    MPI_Recv(&T, 1, MPI_DOUBLE, n, 103, MPI_COMM_WORLD, &status);
	    MPI_Recv(&fsz, 1, MPI_INT, n, 104, MPI_COMM_WORLD, &status);
	    
	    for (int j=0; j<fsz; j++) {
	      // Get the field name
	      //
	      MPI_Probe(n, 105, MPI_COMM_WORLD, &status);
	      MPI_Get_count(&status, MPI_CHAR, &l);
	      MPI_Recv(bf.data(), l, MPI_CHAR, n, 105, MPI_COMM_WORLD, &status);
	      std::string s(bf.data(), l);
	      
	      // Sanity check
	      //
	      if (frame.find(s) == frame.end()) {
		std::cerr << "Error finding <" << s << "> in field map"
			  << std::endl;
	      }

	      // Get the data
	      //
	      MPI_Recv(frame[s].data(), frame[s].size(), MPI_FLOAT, n, 106,
		       MPI_COMM_WORLD, &status);
	    }
	    
	    ret[T] = frame;
#ifdef DEBUG
	    rusage usage;
	    int err = getrusage(RUSAGE_SELF, &usage);
	    std::cout << "volumes: T=" << std::setw(8) << std::fixed<< T
		      << " Size=" << std::setw(8) << usage.ru_maxrss/1024/1024
		      << std::endl;
#endif
	  }
	}
      }
    }

    return ret;
  }
  
  void FieldGenerator::file_volumes(BasisClasses::BasisPtr basis,
				    CoefClasses::CoefsPtr  coefs,
				    const std::string      prefix,
				    const std::string      outdir)
  {
    auto db = volumes(basis, coefs);

    int bunch = db.size()/numprocs;
    int first = bunch*myid;
    int last  = std::min<int>(bunch*(myid+1), db.size());

#pragma omp parallel for
    for (int icnt=first; icnt<last; icnt++) {

      auto it = db.begin();
      std::advance(it, icnt);

      DataGrid datagrid(grid[0], grid[1], grid[2],
			pmin[0], pmax[0],
			pmin[1], pmax[1],
			pmin[2], pmax[2]);

      std::vector<double> tmp(grid[0]*grid[1]*grid[2]);

      for (auto & v : it->second) {
	  
	for (int i=0; i<grid[0]; i++) {
	  for (int j=0; j<grid[1]; j++) {
	    for (int k=0; k<grid[2]; k++) {
	      tmp[(k*grid[1] + j)*grid[0] + i] = v.second(i, j, k);
	    }
	  }
	}
	
	datagrid.Add(tmp, v.first);
      }
      
      std::ostringstream sout;
      sout << outdir << "/" << prefix << "_volume_" << icnt;
      datagrid.Write(sout.str());
    }
  }
  
  std::map<std::string, Eigen::MatrixXf>
  FieldGenerator::histogram2d(PR::PRptr reader, std::vector<double> ctr)
  {

    std::map<std::string, Eigen::MatrixXf> ret;
    std::map<std::string, double> fac;

    std::vector<double> del(3, 0.0);
    for (int k=0; k<3; k++) {
      if (grid[k]>0) del[k] = (pmax[k] - pmin[k])/grid[k];
    }

    if (grid[0]>0 and grid[1]>0) {
      ret["xy"] = Eigen::MatrixXf::Zero(grid[0], grid[1]);
      fac["xy"] = 1.0/(del[0]*del[1]);
    }

    if (grid[0]>0 and grid[2]>0) {
      ret["xz"] = Eigen::MatrixXf::Zero(grid[0], grid[2]);
      fac["xz"] = 1.0/(del[0]*del[2]);
    }

    if (grid[1]>0 and grid[2]>0) {
      ret["yz"] = Eigen::MatrixXf::Zero(grid[1], grid[2]);
      fac["yz"] = 1.0/(del[1]*del[2]);
    }
    
    std::vector<double> pp(3);
    std::vector<bool> bb(3);

    for (auto p=reader->firstParticle(); p!=0; p=reader->nextParticle()) {

      for (int k=0; k<3; k++) {
	pp[k] = p->pos[k] - ctr[k];
	bb[k] = pp[k] >= pmin[k] and pp[k] < pmax[k] and del[k] > 0.0;
      }

      // x-y
      if (bb[0] and bb[1]) {
	int indx1 = floor( (pp[0] - pmin[0])/del[0] );
	int indx2 = floor( (pp[1] - pmin[1])/del[1] );

	if (indx1>=0 and indx1<grid[0] and indx2>=0 and indx2<grid[1])
	  ret["xy"](indx1, indx2) += p->mass * fac["xy"];
      }

      // x-z
      if (bb[0] and bb[2]) {
	int indx1 = floor( (pp[0] - pmin[0])/del[0] );
	int indx2 = floor( (pp[2] - pmin[2])/del[2] );

	if (indx1>=0 and indx1<grid[0] and indx2>=0 and indx2<grid[2] )
	  ret["xz"](indx1, indx2) += p->mass * fac["xz"];
      }

      // y-z
      if (bb[1] and bb[2]) {
	int indx1 = floor( (pp[1] - pmin[1])/del[1] );
	int indx2 = floor( (pp[2] - pmin[2])/del[2] );

	if (indx1>=0 and indx1<grid[1] and indx2>=0 and indx2<grid[2] )
	  ret["yz"](indx1, indx2) += p->mass * fac["yz"];
      }
      
    }

    if (use_mpi) {
      for (auto & v : ret) {
	if (myid==0) 
	  MPI_Reduce(MPI_IN_PLACE, v.second.data(), v.second.size(),
		     MPI_FLOAT, MPI_SUM, 0, MPI_COMM_WORLD);
	else
	  MPI_Reduce(v.second.data(), NULL, v.second.size(),
		     MPI_FLOAT, MPI_SUM, 0, MPI_COMM_WORLD);
      }
    }

    return ret;
  }
  // END histogram

  Eigen::VectorXf
  FieldGenerator::histogram1d(PR::PRptr reader, double rmax, int nbins,
			      std::string proj, std::vector<double> ctr)
  {
    const double pi = 3.14159265358979323846;

    Eigen::VectorXf ret = Eigen::VectorXf::Zero(nbins);
    double del = rmax/nbins;
    
    enum class Projection {xy, xz, yz, r} type;
    if      (proj == "xy") type = Projection::xy;
    else if (proj == "xz") type = Projection::xz;
    else if (proj == "yz") type = Projection::yz;
    else if (proj == "r" ) type = Projection::r;
    else {
      std::ostringstream sout;
      sout << "FieldGenerator::histogram1d: error parsing projection <" << proj
	   << ">.  Must be one of \"xy\", \"xz\", \"yz, \"r\".";
      std::runtime_error(sout.str());
    }

    // Make the histogram
    //
    for (auto p=reader->firstParticle(); p!=0; p=reader->nextParticle()) {
      double rad = 0.0;
      for (int k=0; k<3; k++) {
	double pp = p->pos[k] - ctr[k];
	if      (type == Projection::xy and (k==0 or k==1))
	  rad += pp*pp;		// Cylindrical radius x^2 + y^2	
	else if (type == Projection::xz and (k==0 or k==2))
	  rad += pp*pp;		// Cylindrical radius x^2 + z^2	
	else if (type == Projection::yz and (k==1 or k==2))
	  rad += pp*pp;		// Cylindrical radius y^2 + z^2	
	else if (type == Projection::r)
	  rad += pp*pp;		// Spherical radius
      }

      int indx = floor(sqrt(rad)/del);
      if (indx>=0 and indx<nbins) ret[indx] += p->mass;
    }
    
    // Accumulate between MPI nodes; return value to root node
    //
    if (use_mpi) {
      if (myid==0) 
	MPI_Reduce(MPI_IN_PLACE, ret.data(), ret.size(),
		   MPI_FLOAT, MPI_SUM, 0, MPI_COMM_WORLD);
      else
	MPI_Reduce(ret.data(), NULL, ret.size(),
		   MPI_FLOAT, MPI_SUM, 0, MPI_COMM_WORLD);
    }

    // Inverse area or volume for density norm
    //
    for (int i=0; i<nbins; i++) {
      if (type == Projection::r) // Spherical shells
	ret[i] /= 4.0*pi/3.0*del*del*del*(3*i*(i+1) + 1);
      else			 // Cylindrical shells
	ret[i] /= pi*del*del*(2*i + 1);
    }
    
    return ret;
  }
  // END histogram1d

}
// END namespace Field
