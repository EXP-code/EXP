#include <algorithm>
#include <iomanip>
#include <sstream>

#include <FieldGenerator.H>
#include <DataGrid.H>

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
    if (not use_mpi) {
      int argc = 0; char **argv = 0;
      MPI_Init(&argc, &argv);
    }
  }
  
  void FieldGenerator::check_times(Coefs::CoefsPtr coefs)
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
  
  std::map<double, std::map<std::string, Eigen::MatrixXf>>
  FieldGenerator::slices(Basis::BasisPtr basis, Coefs::CoefsPtr coefs)
  {
    std::map<double, std::map<std::string, Eigen::MatrixXf>> ret;

    std::vector<std::string> labels =
      {"p0", "p1", "p", "fr", "ft", "fp", "d0", "d1", "d", "dd"};

    // Find the first two non-zero indices
    int i1=-1, i2=-1, i3=-1;
    std::vector<double> pos, del;
    for (size_t i=0; i<grid.size(); i++) {
      pos.push_back(pmin[i]);
      if (grid[i]>0) {
	if (i1<0) {
	  i1 = i;
	  del.push_back((pmax[i1] - pmin[i1])/grid[i1]);
	}
	else if (i2<0) {
	  i2 = i;
	  del.push_back((pmax[i2] - pmin[i2])/grid[i2]);
	}
      } else {
	del.push_back(0.0);
	i3 = i;
      }
    } 

    if (i1<0 or i2<0 or i3<0)
      throw std::runtime_error("FieldGenerator::slices: bad grid specification");

    for (auto T : times) {

      if (not coefs->getCoefStruct(T)) {
	std::cout << "Could not find time=" << T << ", continuing" << std::endl;
	continue;
      }

      basis->set_coefs(coefs->getCoefStruct(T));

      std::map<std::string, Eigen::MatrixXf> frame;
      for (auto label : labels) {
	frame[label].resize(grid[i1]+1, grid[i2]+1);
      }	

      double r, phi, costh;
      double p0, p1, d0, d1, fr, ft, fp;
      
      int ncnt = 0;		// Process counter for MPI

      for (int i=0; i<grid[i1]; i++) {

	pos[i1] = pmin[i1] + del[i1]*i;

	for (int j=0; j<grid[i2]; j++) {

	  pos[i2] = pmin[i2] + del[i2]*j;

	
	  if ((ncnt++)%numprocs == myid) {
	  
	    double x = pos[i1];
	    double y = pos[i2];
	    double z = pos[i3];

	    r = sqrt(x*x + y*y + z*z) + 1.0e-18;
	    costh = z/r;
	    phi = atan2(y, x);

	    basis->all_eval(r, costh, phi,
			    d0, d1, p0, p1, fr, ft, fp);

	    frame["p0"](i, j) = p0;
	    frame["p1"](i, j) = p1;
	    frame["p" ](i, j) = p0 + p1;
	    frame["fr"](i, j) = fr;
	    frame["ft"](i, j) = ft;
	    frame["fp"](i, j) = fp;
	    frame["d0"](i, j) = d0;
	    frame["d1"](i, j) = d1;
	    frame["d" ](i, j) = d0 + d1;
	    
	    if (d0!=0.0)
	      frame["dd" ](i, j) = d1/d0;
	    else
	      frame["dd" ](i, j) = 0.0;
	  }
	}
      }
    
      if (use_mpi) {
	for (auto & f : frame) {
	  if (myid==0) 
	    MPI_Reduce(MPI_IN_PLACE, f.second.data(), f.second.size(),
		       MPI_FLOAT, MPI_SUM, 0, MPI_COMM_WORLD);
	  else
	    MPI_Reduce(f.second.data(), NULL, f.second.size(),
		       MPI_FLOAT, MPI_SUM, 0, MPI_COMM_WORLD);
	}
      }

      ret[T] = frame;
    }

    return ret;
  }
  
  void FieldGenerator::file_slices(Basis::BasisPtr basis, Coefs::CoefsPtr coefs,
				   const std::string prefix, const std::string outdir)
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
  FieldGenerator::volumes(Basis::BasisPtr basis, Coefs::CoefsPtr coefs)
  {
    std::map<double, std::map<std::string, Eigen::Tensor<float, 3>>> ret;

    std::vector<std::string> labels =
      {"p0", "p1", "p", "fr", "ft", "fp", "d0", "d1", "d", "dd"};

    auto times = coefs->Times();

    for (auto T : times) {

      basis->set_coefs(coefs->getCoefStruct(T));

      std::map<std::string, Eigen::Tensor<float, 3>> frame;
      for (auto label : labels) {
	frame[label].resize({grid[0]+1, grid[1]+1, grid[2]+1});
      }	

      double r, phi, costh;
      double p0, p1, d0, d1, fr, ft, fp;
      
      std::vector<double> del =
	{ (pmax[0] - pmin[0])/grid[0],
	  (pmax[1] - pmin[1])/grid[1],
	  (pmax[2] - pmin[2])/grid[2] };

      int ncnt = 0;

      for (int i=0; i<grid[0]; i++) {

	double x = pmin[0] + del[0]*i;

	for (int j=0; j<grid[1]; j++) {

	  double y = pmin[1] + del[1]*j;

	  for (int k=0; k<grid[2]; k++) {

	    double z = pmin[2] + del[2]*k;

	    if ((ncnt++)%numprocs == myid) {
	  
	      r = sqrt(x*x + y*y + z*z) + 1.0e-18;
	      costh = z/r;
	      phi = atan2(y, x);

	      basis->all_eval(r, costh, phi, d0, d1, p0, p1, fr, ft, fp);

	      frame["p0"](i, j, k) = p0;
	      frame["p1"](i, j, k) = p1;
	      frame["p" ](i, j, k) = p0 + p1;
	      frame["fr"](i, j, k) = fr;
	      frame["ft"](i, j, k) = ft;
	      frame["fp"](i, j, k) = fp;
	      frame["d0"](i, j, k) = d0;
	      frame["d1"](i, j, k) = d1;
	      frame["d" ](i, j, k) = d0 + d1;
	    
	      if (d0!=0.0)
		frame["dd" ](i, j, k) = d1/d0;
	      else
		frame["dd" ](i, j, k) = 0.0;
	    }
	  }
	}
      }
    
      if (use_mpi) {
	for (auto & f : frame) {
	  if (myid==0) 
	    MPI_Reduce(MPI_IN_PLACE, f.second.data(), f.second.size(),
		       MPI_FLOAT, MPI_SUM, 0, MPI_COMM_WORLD);
	  else
	    MPI_Reduce(f.second.data(), NULL, f.second.size(),
		       MPI_FLOAT, MPI_SUM, 0, MPI_COMM_WORLD);
	}
      }

      ret[T] = frame;
    }

    return ret;
  }
  
  void FieldGenerator::file_volumes(Basis::BasisPtr basis, Coefs::CoefsPtr coefs,
				    const std::string prefix, const std::string outdir)
  {
    auto db = volumes(basis, coefs);

    if (myid==0) {

      int icnt = 0;

      for (auto & frame : db) {

	DataGrid datagrid(grid[0], grid[1], grid[2], pmin[0], pmax[0], pmin[1], pmax[1], pmin[2], pmax[2]);

	std::vector<double> tmp(grid[0]*grid[1]*grid[2]);

	for (auto & v : frame.second) {
	  
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
	icnt++;
      }
    }
  }
  
}
// END namespace Field

  
