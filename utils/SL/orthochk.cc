#include <iostream>
#include <iomanip>
#include <algorithm>
#include <string>
#include <cmath>

#include "biorth1d.H"
#include "SLGridMP2.H"
#include "gaussQ.H"
#include "localmpi.H"
#include "cxxopts.H"


enum BioType1d {Trig, SL};

int 
main(int argc, char** argv)
{
  double KX = 0.5;
  double H = 0.1;
  double ZMAX = 1.0;
  int NMAX = 10;
  int IKX = 1;
  int IKY = 3;
  BioType1d Type = Trig;
  std::string cachename = ".slab_sl_cache";
  std::string slabID = "iso";
  bool use_mpi = false;

  cxxopts::Options options("orthochk", "Check orthogonality of 1D basis functions");

  options.add_options()
    ("m,mpi", "Use MPI")
    ("s,SL", "Use Sturm-Liouville slab basis")
    ("T,type", "Slab type (iso, parabolic, or constant)", cxxopts::value<std::string>())
    ("t,Trig", "Use trigonometric basis")
    ("x,ikx", "IKX for SLGridSlab (default: 1)", cxxopts::value<int>())
    ("y,iky", "IKY for SLGridSlab (default: 3)", cxxopts::value<int>())
    ("k,kx", "KX for OneDTrig (default: 0.5)", cxxopts::value<double>())
    ("z,zmax", "ZMAX for OneDTrig and SLGridSlab (default: 1.0)", cxxopts::value<double>())
    ("H,h", "Scale height H for SLGridSlab (default: 0.1)", cxxopts::value<double>())
    ("n,nmax", "NMAX for SLGridSlab (default: 10)", cxxopts::value<int>())
    ("c,cachename", "Cache file name for SLGridSlab (default: .slab_sl_cache)", cxxopts::value<std::string>())
    ("h,help", "Print usage");

  auto result = options.parse(argc, argv);

  if (result.count("mpi")) {
    local_init_mpi(argc, argv);
    use_mpi = true;
  }

  if (result.count("help")) {
    if (myid==0)
      std::cout << options.help() << std::endl;
    return 0;
  }


  //===================
  // Construct ortho
  //===================
  
  std::shared_ptr<OneDTrig> ortho;
  std::shared_ptr<SLGridSlab> orthoSL;

  switch (Type) {
  case Trig:
    ortho = std::make_shared<OneDTrig>(KX, ZMAX);
    break;

  case SL:
    {
      const int NUMZ=800;
      int KMAX = max<int>(IKX+1, IKY+1);
      SLGridSlab::ZBEG = 0.0;
      SLGridSlab::ZEND = 0.1;
      SLGridSlab::H = H;
      if (use_mpi) SLGridSlab::mpi = 1;

      orthoSL = std::make_shared<SLGridSlab>(KMAX, NMAX, NUMZ, ZMAX, cachename, slabID, true);
    }
    break;

  default:
    cerr << "No such one-dimensional orthogonal function type: " << Type 
	 << endl;
    exit(0);
  }


  //===================
  // Get info
  //===================

  if (!use_mpi || myid==0) {

    while (1) {
      
      bool done=false;
      int iwhich;
      
      cout << "Task:" << endl;
      cout << "1: Print out density, potential pairs" << endl;
      cout << "2: Check density" << endl;
      cout << "3: Check orthogonality" << endl;
      cout << "4: Quit" << endl;
      cout << "?? ";
      cin >> iwhich;
      
      if (iwhich < 1 || iwhich > 4) iwhich = 4;

      switch(iwhich) {
      case 1:
	{
	  string filename;
	  cout << "Filename? ";
	  cin >> filename;
	  ofstream out (filename.c_str());
	  if (!out) {
	    cout << "Can't open <" << filename << "> for output" << endl;
	    break;
	  }
	  
	  cout << "Number of points? ";
	  int num;
	  cin >> num;
	  
	  cout << "N? ";
	  int N;
	  cin >> N;
	  
	  double x, z;
	  
	  for (int i=0; i<num; i++) {
	    z = -ZMAX + 2.0*ZMAX*(0.5 + i)/num;
	    if (Type == Trig) {
	      x = ortho->r_to_rb(z);
	      out << setw(15) << z
		  << setw(15) << ortho->potl(N, i, z)
		  << setw(15) << ortho->force(N, i, z)
		  << setw(15) << ortho->dens(N, i, z)
		  << endl;
	    } else {
	      x = orthoSL->z_to_xi(z);
	      out << setw(15) << z
		  << setw(15) << orthoSL->get_pot(x, IKX, IKY, N)
		  << setw(15) << orthoSL->get_force(x, IKX, IKY, N)
		  << setw(15) << orthoSL->get_dens(x, IKX, IKY, N)
		  << endl;
	    }
	  }
	}
	
	break;
	
      case 2:
	{
	  string filename;
	  cout << "Filename? ";
	  cin >> filename;
	  ofstream out (filename.c_str());
	  if (!out) {
	    cout << "Can't open <" << filename << "> for output" << endl;
	    break;
	  }
	  
	  cout << "Number of points? ";
	  int num;
	  cin >> num;
	  
	  cout << "N? ";
	  int N;
	  cin >> N;
	  
	  cout << "dz? ";
	  double dz;
	  cin >> dz;
	  
	  double x, z, d1, d2, d3;
	  
	  for (int i=0; i<num; i++) {
	    z = -ZMAX + 2.0*ZMAX*(0.5 + i)/num;
	    if (Type == Trig) {
	      x = ortho->r_to_rb(z);

	      d1 = (
		    ortho->potl(N, i, z+dz)    -
		    ortho->potl(N, i, z)*2.0   +
		    ortho->potl(N, i, z-dz)
		    ) / (dz*dz);

	      d2 = (
		    ortho->force(N, i, z-0.5*dz) -
		    ortho->force(N, i, z+0.5*dz)
		    ) / dz;

	      d3 = -KX*KX*ortho->potl(N, i, z);


	      out << setw(15) << z
		  << setw(15) << d1+d3
		  << setw(15) << d2+d3
		  << setw(15) << -ortho->dens(N, i, z)
		  << endl;
	    } else {
	      x = orthoSL->z_to_xi(z);

	      d1 = (
		    orthoSL->get_pot(x+dz, IKX, IKY, N)   -
		    orthoSL->get_pot(x, IKX, IKY, N)*2.0  +
		    orthoSL->get_pot(x-dz, IKX, IKY, N) 
		    ) / (dz*dz);

	      d2 = (
		    orthoSL->get_force(x+0.5*dz, IKX, IKY, N) -
		    orthoSL->get_force(x-0.5*dz, IKX, IKY, N)
		    ) / dz;

	      d3 = -4.0*M_PI*M_PI*(IKX*IKX+IKY*IKY)*
		orthoSL->get_pot(x, IKX, IKY, N);

	      out << setw(15) << z
		  << setw(15) << d1+d3
		  << setw(15) << d2+d3
		  << setw(15) << orthoSL->get_dens(x, IKX, IKY, N)
		  << endl;

	    }
	  }
	}
	
	break;
	
      case 3:
	{
	  std::cout << "Number of knots? ";
	  int num;
	  std::cin >> num;
	  
	  LegeQuad lw(num);
	  
	  std::cout << "N1, N2? ";
	  int N1, N2;
	  std::cin >> N1;
	  std::cin >> N2;
	  
	  double ximin, ximax;
	  switch (Type) {
	  case Trig:
	    ximin = ortho->r_to_rb(-ZMAX);
	    ximax = ortho->r_to_rb( ZMAX);
	    break;
	  case SL:
	    ximin = orthoSL->z_to_xi(-ZMAX);
	    ximax = orthoSL->z_to_xi( ZMAX);
	    break;
	  }
	  
	  double x, r, ans=0.0;
	  for (int i=0; i<num; i++) {
	    
	    x = ximin + (ximax - ximin)*lw.knot(i);
	    
	    switch (Type) {
	    case Trig:
	      {
	      
	      double tmp1 = ortho->potl(N1, i, x);
	      double tmp2 = ortho->dens(N1, i, x);
	      double tmp3 = ortho->d_r_to_rb(x);

	      ans += ortho->potl(N1, i, x)*
		ortho->dens(N2, i, x) *
		ortho->d_r_to_rb(x) * (ximax - ximin)*lw.weight(i);
	      }
	      
	      break;
	      
	    case SL:
	      
	      ans += orthoSL->get_pot(x, IKX, IKY, N1)*
		orthoSL->get_dens(x, IKX, IKY, N2) /
		orthoSL->d_xi_to_z(x) * (ximax - ximin)*lw.weight(i);
	      
	      break;
	    }
	  }
	  
	  cout << "<" << N1 << "|" << N2 << "> = " << ans << endl;
	}
	
	break;
	
      default:
	done = true;
	break;
      }
      
      if (done) break;
    }
  }
  
  if (use_mpi) MPI_Finalize();
  
  return 0;
}

