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

std::map<std::string, BioType1d> bioTypeMap = {
  {"trig", Trig},
  {"sl",   SL}
};

std::map<BioType1d, std::string>  bioStrMap = {
  {Trig, "trig"},
  {SL,   "sl"}
};

int 
main(int argc, char** argv)
{
  double KX = 0.5;
  double H = 0.1;
  double ZMAX = 1.0;
  int NMAX = 10;
  int IKX = 1;
  int IKY = 3;
  int numz = 2000;
  std::string bioTypeStr = "trig";
  std::string cachename = ".slab_sl_cache";
  std::string cmap = "linear";
  std::string slabID = "iso";
  bool use_mpi = false;

  cxxopts::Options options("orthochk", "Check orthogonality of 1D basis functions");

  options.add_options()
    ("m,mpi",       "Use MPI")
    ("T,type",      "Slab type (trig or SL)", cxxopts::value<std::string>(bioTypeStr)->default_value("trig"))
    ("M,matrix",    "Print orthogonality matrix for a particular kx, ky choice")
    ("x,ikx",       "IKX for SLGridSlab (default: 1)", cxxopts::value<int>(IKX))
    ("y,iky",       "IKY for SLGridSlab (default: 3)", cxxopts::value<int>(IKY))
    ("k,kx",        "KX for OneDTrig (default: 0.5)", cxxopts::value<double>(KX))
    ("z,zmax",      "ZMAX for OneDTrig and SLGridSlab (default: 1.0)", cxxopts::value<double>(ZMAX))
    ("H,height",    "Scale height H for SLGridSlab (default: 0.1)", cxxopts::value<double>(H))
    ("n,nmax",      "NMAX for SLGridSlab (default: 10)", cxxopts::value<int>(NMAX))
    ("N,numz",      "Number of z points for SLGridSlab (default: 2000)", cxxopts::value<int>(numz)->default_value("2000"))
    ("C,cmap",       "Coordinate map for SLGridSlab (linear, tanh, or sech; default: linear)", cxxopts::value<std::string>(cmap)->default_value("linear"))
    ("s,slabid",    "Slab model ID for SLGridSlab (iso, parabolic, or constant; default: iso)", cxxopts::value<std::string>(slabID)->default_value("iso"))
    ("c,cachename", "Cache file name for SLGridSlab (default: .slab_sl_cache)", cxxopts::value<std::string>(cachename)->default_value(".slab_sl_cache"))
    ("h,help",      "Print usage");

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


  // Determine biorthogonal function type
  std::transform(bioTypeStr.begin(), bioTypeStr.end(), bioTypeStr.begin(),
		 [](unsigned char c){ return std::tolower(c); });

  BioType1d Type = std::find_if(bioTypeMap.begin(), bioTypeMap.end(),
				[bioTypeStr](const std::pair<std::string, BioType1d>& pair) { return pair.first == bioTypeStr; }) != bioTypeMap.end() ? bioTypeMap[bioTypeStr] : Trig;

  // Check for valid coordinate map choice
  std::transform(cmap.begin(), cmap.end(), cmap.begin(),
		 [](unsigned char c){ return std::tolower(c); });

  std::vector<std::string> valid_cmaps = { "tanh", "sech", "linear" };
  if (valid_cmaps.end() ==
      std::find_if(valid_cmaps.begin(), valid_cmaps.end(),
		   [cmap](const std::string& vmap) { return vmap == cmap; })) {
    std::cerr << "Invalid coordinate map choice: " << cmap << std::endl;
    return 1;
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
      int KMAX = max<int>(IKX+1, IKY+1);
      SLGridSlab::ZBEG = 0.0;
      SLGridSlab::ZEND = 0.1;
      SLGridSlab::H = H;
      if (use_mpi) SLGridSlab::mpi = 1;

      orthoSL = std::make_shared<SLGridSlab>(KMAX, NMAX, numz, ZMAX, cachename,
					     cmap, slabID, true);
    }
    break;

  default:
    cerr << "No such one-dimensional orthogonal function type: " << Type 
	 << endl;
    exit(0);
  }


  int ikx = max(IKX, IKY), iky = min(IKX, IKY);

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
      cout << "4: Internal ortho check" << endl;
      cout << "5: Quit" << endl;
      cout << "?? ";
      cin >> iwhich;
      
      if (iwhich < 1 || iwhich > 5) iwhich = 5;

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
		  << setw(15) << ortho->potl (N, i, z)
		  << setw(15) << ortho->force(N, i, z)
		  << setw(15) << ortho->dens (N, i, z)
		  << endl;
	    } else {
	      x = orthoSL->z_to_xi(z);
	      out << setw(15) << z
		  << setw(15) << orthoSL->get_pot  (x, ikx, iky, N)
		  << setw(15) << orthoSL->get_force(x, ikx, iky, N)
		  << setw(15) << orthoSL->get_dens (x, ikx, iky, N)
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
	  
	  cout << "eps? ";
	  double eps;
	  cin >> eps;
	  
	  double dz = 2.0*ZMAX/num;
	  double h  = dz*eps;
	  double x, z, d1, d2, d3;
	  
	  for (int i=0; i<num; i++) {
	    z = -ZMAX + dz*(0.5 + i);

	    out << setw(15) << z;

	    if (Type == Trig) {

	      for (int n=0; n<NMAX; n++) {

		d1 = (
		      ortho->potl(n, i, z+h)    -
		      ortho->potl(n, i, z)*2.0   +
		      ortho->potl(n, i, z-h)
		      ) / (h*h);
		
		d2 = (
		      ortho->force(n, i, z-0.5*h) -
		      ortho->force(n, i, z+0.5*h)
		      ) / h;

		d3 = -KX*KX*ortho->potl(n, i, z);


		out << setw(15) << d1+d3
		    << setw(15) << d2+d3
		    << setw(15) << -ortho->dens(n, i, z);
	      }

	    } else {

	      Eigen::VectorXd pot0(NMAX), potN(NMAX), potP(NMAX), den0(NMAX);
	      Eigen::VectorXd frcP(NMAX), frcN(NMAX);

	      orthoSL->get_pot  (pot0, z  ,     ikx, iky);
	      orthoSL->get_pot  (potN, z-h,     ikx, iky);
	      orthoSL->get_pot  (potP, z+h,     ikx, iky);
	      orthoSL->get_force(frcN, z-0.5*h, ikx, iky);
	      orthoSL->get_force(frcP, z+0.5*h, ikx, iky);
	      orthoSL->get_dens (den0, z,       ikx, iky);

	      for (int n=0; n<NMAX; n++) {

		d1 = (potP(n) - 2.0*pot0(n) + potN(n)) / (h*h);
		d2 = (frcP(n) - frcN(n)) / h;
		d3 = -4.0*M_PI*M_PI*(IKX*IKX+IKY*IKY) * pot0(n);

		out << setw(15) << d1+d3
		    << setw(15) << d2+d3
		    << setw(15) << den0(n);
	      }
	    }
	    out << endl;
	  }
	}
	
	break;
	
      case 3:
	{
	  std::cout << "Number of knots? ";
	  int num;
	  std::cin >> num;
	  
	  LegeQuad lw(num);
	  
	  if (result.count("matrix")) {

	    Eigen::VectorXd pot(NMAX), den(NMAX);
	    Eigen::MatrixXd ans(NMAX, NMAX);
	    ans.setZero();

	    for (int i=0; i<num; i++) {
	      double z = -ZMAX + 2.0*ZMAX*lw.knot(i);
	      double W = 2.0*ZMAX*lw.weight(i);

	      if (Type == Trig) {
		ortho->potl(i, i, z, pot);
		ortho->dens(i, i, z, den);
	      } else {
		orthoSL->get_pot (pot, z, ikx, iky);
		orthoSL->get_dens(den, z, ikx, iky);
	      }

	      for (int n1=0; n1<NMAX; n1++) {
		for (int n2=0; n2<NMAX; n2++) {
		  ans(n1, n2) += pot(n1) * den(n2) * W;
		}
	      }
	    }

	    std::cout << std::endl << ans << std::endl;
	  }
	  else {
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
	  
	    double z, r, ans=0.0;
	    for (int i=0; i<num; i++) {
	    
	      z = -ZMAX + 2.0*ZMAX*lw.knot(i);
	    
	      switch (Type) {
	      case Trig:
		{
		  
		  double tmp1 = ortho->potl(N1, i, z);
		  double tmp2 = ortho->dens(N1, i, z);
		  
		  ans += ortho->potl(N1, i, z)*
		    ortho->dens(N2, i, z) * 2.0*ZMAX*lw.weight(i);
		}
	      
		break;
	      
	      case SL:
	      
		ans += orthoSL->get_pot(z, IKX, IKY, N1)*
		  orthoSL->get_dens(z, IKX, IKY, N2) * 2.0*ZMAX*lw.weight(i);
		
		break;
	      }
	    }
	    std::cout << "<" << N1 << "|" << N2 << "> = " << ans << std::endl;
	  }
	  break;
	}
	case 4:
	{
	  switch (Type) {
	  case Trig:
	    std::cout << "No internal orthogonality check for OneDTrig"
		      << std::endl;
	    break;
	  case SL:
	    {
	      int knots;
	      std::cout << "Number of knots? ";
	      std::cin >> knots;
	      auto orthoMat = orthoSL->orthoCheck(knots);
	      std::cout << "Orthogonality check for SLGridSlab" << std::endl;
	      for (size_t i=0; i<orthoMat.size(); i++)
		std::cout << orthoMat[i] << std::endl << std::endl;
	    }
	    break;
	  }
	}
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

