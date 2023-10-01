/*****************************************************************************
 *  Description:
 *  -----------
 *
 *  Create EOF from a sequence of PSP files.  Compute per component
 *  grid generation for insight.
 *
 *
 *  Call sequence:
 *  -------------
 *
 *  Parameters:
 *  ----------
 *
 *
 *  Returns:
 *  -------
 *
 *
 *  Notes:
 *  -----
 *
 *
 *  By:
 *  --
 *
 *  MDW 12/11/20
 *
 ***************************************************************************/

				// C++/STL headers
#include <cstdlib>
#include <iostream>
#include <iomanip>
#include <fstream>
#include <sstream>
#include <memory>
#include <cmath>
#include <string>

#include <config_exp.h>		// CMake generated config

#ifdef HAVE_LIBPNGPP
#include <ColorGradient.H>	// For PNG images
#endif

#include <yaml-cpp/yaml.h>	// YAML support

#include <Eigen/Eigen>		// Eigen 3

#include <DataGrid.H>		// For VTK or ASCII grid output
#include <cxxopts.H>		// Command-line parser

                                // System libs
#include <sys/time.h>
#include <sys/resource.h>

				// MDW classes
#include <numerical.H>
#include <ParticleReader.H>
#include <interp.H>
#include <EmpCylSL.H>

#include <localmpi.H>
#include <foarray.H>
#include <libvars.H>

#ifdef DEBUG
#ifndef _REDUCED
#pragma message "NOT using reduced particle structure"
#endif
#endif

using CoefArray  = std::vector<Eigen::VectorXd>;
using DvarArray  = std::vector<Eigen::MatrixXd>;

const std::string overview = "Compute new EOF basis PSP phase-space output files";
 
void write_output(EmpCylSL& ortho, int t, int m, int nmin, int nord,
		  double RMAX, int OUTR, std::string& prefix,
		  std::vector<CoefArray>& cc, std::vector<CoefArray>& ss)
{
  unsigned ncnt = 0;
  
  // ==================================================
  // Write surface profile
  //   --- in plane ---
  // ==================================================
    
  double v;
  float f;
    
  double dR = 2.0*RMAX/(OUTR-1);
  double x, y, z=0.0, r, phi;
  double p0, d0, p, fr, fz, fp;

  std::vector<double> indat(2*nord*OUTR*OUTR, 0.0);
  std::vector<double> otdat(2*nord*OUTR*OUTR);
  
  for (int j=0; j<OUTR; j++) {
	
    x = -RMAX + dR*j;
    
    for (int l=0; l<OUTR; l++) {
      
      y = -RMAX + dR*l;
      
      if ((ncnt++)%numprocs == myid) {
	  
	r = sqrt(x*x + y*y);
	phi = atan2(y, x);

	double cosp = cos(phi*m);
	double sinp = sin(phi*m);
	  
	for (int n=0; n<nord; n++) {

	  double sumP = 0.0, sumD = 0.0;
	  
	  for (int k=0; k<nord; k++) {
	    double dC, dS;

	    ortho.getDensSC(m, n+nmin, r, z, dC, dS);
	    sumD += dC*cc[t][m][n]*cosp + dS*ss[t][m][n]*sinp;

	    ortho.getPotSC(m, n+nmin, r, z, dC, dS);
	    sumP += dC*cc[t][m][n]*cosp + dS*ss[t][m][n]*sinp;
	  }

	  indat[nord*((0*OUTR+j)*OUTR+l) + n] = sumD;
	  indat[nord*((1*OUTR+j)*OUTR+l) + n] = sumP;
	}
      }
    }
  }
    
  MPI_Reduce(&indat[0], &otdat[0], 2*nord*OUTR*OUTR,
	     MPI_DOUBLE, MPI_SUM, 0, MPI_COMM_WORLD);
    
    
  if (myid==0) {
      
    for (int n=0; n<nord; n++) {

      DataGrid grid(OUTR, OUTR, 1, -RMAX, RMAX, -RMAX, RMAX, 0, 0);

      std::vector<double> dataD(OUTR*OUTR), dataP(OUTR*OUTR);

      for (int j=0; j<OUTR; j++) {
	for (int l=0; l<OUTR; l++) {
	  dataD[j*OUTR + l] = otdat[nord*((0*OUTR+j)*OUTR+l) + n];
	  dataP[j*OUTR + l] = otdat[nord*((1*OUTR+j)*OUTR+l) + n];
	}
      }
      grid.Add(dataD, "d");
      grid.Add(dataP, "p");

      std::ostringstream sout;
      sout << prefix << "_rotated_" << runtag 
	   << "." << std::setfill('0') << std::setw(5) << m
	   << "." << std::setfill('0') << std::setw(5) << n+nmin
	   << "." << std::setfill('0') << std::setw(3) << t;
      
      grid.Write(sout.str());
    }
  }

}

int
main(int argc, char **argv)
{
  // ==================================================
  // MPI preliminaries
  // ==================================================
  local_init_mpi(argc, argv);
  

  // ==================================================
  // Parameter parsing
  // ==================================================
  int lmax=64, mmax, Nmin, Nmax, nmax, norder, numx, numy, cmapr=1, cmapz=1, nodd=-1;
  double rcylmin, rcylmax, rscale, vscale, RMAX;
  std::string CACHEFILE, COEFFILE, cname, prefix, fileType, psfile, delim;
  int beg, end, stride, mbeg, mend, OUTR;

  cxxopts::Options options(argv[0],
			   "Create EOF from a sequence of PSP files.\nCompute per component grid generation for insight.\n");

  options.add_options()
    ("h,help", "produce this help message")
    ("v,verbose", "verbose output")
    ("PNG", "PNG matrix output")
    ("gray", "use gray map for PNG matrix output")
    ("5col", "use five color heat map for PNG matrix output")
    ("F,filetype", "input file type",
    cxxopts::value<std::string>(fileType)->default_value("PSPout"))
    ("psfile", "List of phase space files for processing",
     cxxopts::value<std::string>(psfile))
    ("delimiter", "Phase-space file list delimiter for node index",
     cxxopts::value<std::string>(delim))
    ("r,rmax", "maximum output radius",
     cxxopts::value<double>(RMAX)->default_value("0.03"))
    ("n,nout", "number of points on a side for grid output",
     cxxopts::value<int>(OUTR)->default_value("50"))
    ("mbeg", "minimum azimuthal order for grid; output off if m<0",
     cxxopts::value<int>(mbeg)->default_value("2"))
    ("mend", "maximum azimuthal order for grid",
     cxxopts::value<int>(mend)->default_value("2"))
    ("nmin", "minimum order in EOF",
     cxxopts::value<int>(Nmin)->default_value("0"))
    ("nmax", "maximum order in EOF",
     cxxopts::value<int>(Nmax)->default_value(std::to_string(std::numeric_limits<int>::max())))
    ("outdir", "Output directory path",
     cxxopts::value<std::string>(outdir)->default_value("."))
    ("cachefile", "cachefile name",
     cxxopts::value<std::string>(CACHEFILE)->default_value(".eof.cache.file"))
    ("coeffile", "coefficient output file name",
     cxxopts::value<std::string>(COEFFILE))
    ("cname", "component name",
     cxxopts::value<std::string>(cname)->default_value("star"))
    ("runtag", "runtag for phase space files",
     cxxopts::value<std::string>(runtag)->default_value("run1"))
    ("prefix", "output prefix for distinguishing parameters",
     cxxopts::value<std::string>(prefix)->default_value("even"))
    ;
  
  cxxopts::ParseResult vm;

  try {
    vm = options.parse(argc, argv);
  } catch (cxxopts::OptionException& e) {
    std::cout << "Option error: " << e.what() << std::endl;
    exit(-1);
  }

  if (vm.count("help")) {
    if (myid==0) {
      std::cout << std::string(60, '-') << std::endl;
      std::cout << options.help() << std::endl;
      std::cout << std::string(60, '-') << std::endl << std::endl;
    }
    MPI_Finalize();
    return 0;
  }
 
  bool PNG = false;
  if (vm.count("PNG")) PNG = true;

#ifdef DEBUG
  sleep(20);
#endif  

  // ==================================================
  // Read basis cache
  // ==================================================

  std::ifstream in(CACHEFILE);
  if (!in) {
    if (myid==0)
      std::cerr << "Error opening cachefile named <" 
		<< CACHEFILE << "> . . . I quit!"
		<< std::endl;
    MPI_Finalize();
    exit(-1);

  } else {

    // Attempt to read magic number
    //
    unsigned int tmagic;
    in.read(reinterpret_cast<char*>(&tmagic), sizeof(unsigned int));

    //! Basis magic number
    const unsigned int hmagic = 0xc0a57a1;

    if (tmagic == hmagic) {
      // YAML size
      //
      unsigned ssize;
      in.read(reinterpret_cast<char*>(&ssize), sizeof(unsigned int));
      
      // Make and read char buffer
      //
      auto buf = std::make_unique<char[]>(ssize+1);
      in.read(buf.get(), ssize);
      buf[ssize] = 0;		// Null terminate
      
      YAML::Node node;
      
      try {
	node = YAML::Load(buf.get());
      }
      catch (YAML::Exception& error) {
	if (myid)
	  std::cerr << "YAML: error parsing <" << buf.get() << "> "
		    << "in " << __FILE__ << ":" << __LINE__ << std::endl
		    << "YAML error: " << error.what() << std::endl;
	throw error;
      }
      
      // Get parameters
      //
      mmax    = node["mmax"  ].as<int>();
      numx    = node["numx"  ].as<int>();
      numy    = node["numy"  ].as<int>();
      nmax    = node["nmax"  ].as<int>();
      norder  = node["norder"].as<int>();
      if (node["nodd"])
	nodd  = node["nodd"  ].as<int>();
      if (node["cmap"])
	cmapr = node["cmap"  ].as<int>();
      else
	cmapr = node["cmapr" ].as<int>();
      if (node["cmapz"])
	cmapz = node["cmapz"  ].as<int>();
      rcylmin = node["rmin"  ].as<double>();
      rcylmax = node["rmax"  ].as<double>();
      rscale  = node["ascl"  ].as<double>();
      vscale  = node["hscl"  ].as<double>();

    } else {
				// Rewind file
      in.clear();
      in.seekg(0);

      int tmp;
    
      in.read((char *)&mmax,    sizeof(int));
      in.read((char *)&numx,    sizeof(int));
      in.read((char *)&numy,    sizeof(int));
      in.read((char *)&nmax,    sizeof(int));
      in.read((char *)&norder,  sizeof(int));
      in.read((char *)&tmp,     sizeof(int));   
      in.read((char *)&cmapr,   sizeof(int)); 
      in.read((char *)&rcylmin, sizeof(double));
      in.read((char *)&rcylmax, sizeof(double));
      in.read((char *)&rscale,  sizeof(double));
      in.read((char *)&vscale,  sizeof(double));
    }
  }

  EmpCylSL::RMIN        = rcylmin;
  EmpCylSL::RMAX        = rcylmax;
  EmpCylSL::NUMX        = numx;
  EmpCylSL::NUMY        = numy;
  EmpCylSL::CMAPR       = cmapr;
  EmpCylSL::CMAPZ       = cmapz;
  EmpCylSL::logarithmic = true;

  // Create expansion instance
  //
  EmpCylSL ortho(nmax, lmax, mmax, norder, rscale, vscale, nodd, CACHEFILE);
    
  // ==================================================
  // Initialize and/or create basis
  // ==================================================
  
  if (ortho.read_cache()==0) {
    if (myid==0)
      std::cerr << "Error creating disk basis using <"
		<< CACHEFILE << "> . . . I quit!"
		<< std::endl;
    MPI_Finalize();
    exit(-1);
  }
  
  // ==================================================
  // Create new basis and coefficients
  // ==================================================
  
  // First create storage
  //
  Nmin = std::max<int>(Nmin, 0);
  Nmax = std::min<int>(Nmax, norder);
  int Nord = Nmax - Nmin;

  DvarArray D(mmax+1);
  for (auto & mat : D) {
    mat = Eigen::MatrixXd::Zero(Nord, Nord);
  }

  std::vector<CoefArray> coefsC, coefsS, coefsT, coefRC, coefRS;
  std::vector<double>    times;
  std::vector<int>       indices;

  unsigned indx = 0;

  for (auto batch : PR::ParticleReader::parseFileList(psfile, delim)) {
    
    // ==================================================
    // Open frame list
    // ==================================================
    
    PR::PRptr reader = PR::ParticleReader::createReader
      (fileType, batch, myid, true);

    reader->SelectType(cname);

    double tnow = reader->CurrentTime();

    if (myid==0) {
      cout << "Beginning disk partition [time=" << tnow
	   << ", index=" << indx << "] . . . "  << flush;
    }

    times.push_back(tnow);
    indices.push_back(indx);
    
    auto p = reader->firstParticle();

    CoefArray coefC(mmax + 1);	// Per snapshot storage for
    CoefArray coefS(mmax + 1);	// coefficients

    for (auto & v : coefC) v = Eigen::VectorXd::Zero(Nord);
    for (auto & v : coefS) v = Eigen::VectorXd::Zero(Nord);
    
    unsigned count = 0;
    
    EmpCylSL::ContribArray retC, retS;

    while (p) {
      if (count++ % numprocs == myid) {
	// Only need mass and position
	//
	double m = p->mass;
	double x = p->pos[0];
	double y = p->pos[1];
	double z = p->pos[2];
	double p = atan2(y, x);

	// Get coefficient contribution for this particle
	//
	ortho.getPotParticle(x, y, z, retC, retS);
	
	// Accumulate coefficients and D matrix
	//
	for (int mm=0; mm<=mmax; mm++) {
	  for (int n1=0; n1<Nord; n1++) {
	    // Coefficient contribution
	    coefC[mm](n1) += m * retC[mm](n1+Nmin);
	    if (mm) coefS[mm](n1) += m * retS[mm](n1+Nmin);

	    // Modulus for index n1
	    double mod1 = retC[mm](n1+Nmin) * retC[mm](n1+Nmin);
	    if (mm) mod1 += retS[mm](n1+Nmin) * retS[mm](n1+Nmin);

	    for (int n2=0; n2<Nord; n2++) {
	      // Modulus for index n2
	      double mod2 = retC[mm](n2+Nmin) * retC[mm](n2+Nmin);
	      if (mm) mod2 += retS[mm](n2+Nmin) * retS[mm](n2+Nmin);

	      D[mm](n1, n2) += m * sqrt(mod1 * mod2);
	    }
	  }
	}
      }
      p = reader->nextParticle();
    }

    coefsC.push_back(coefC);
    coefsS.push_back(coefS);

    if (myid==0) std::cout << "done" << std::endl;

    indx++;

  } // PSP loop

  // Full reduction
  //
  for (auto & mat : D) {
    MPI_Allreduce(MPI_IN_PLACE, mat.data(), mat.size(), MPI_DOUBLE, MPI_SUM,
		  MPI_COMM_WORLD);
  }

  for (auto & coef : coefsC) {
    for (auto & d : coef) 
      MPI_Allreduce(MPI_IN_PLACE, d.data(), d.size(), MPI_DOUBLE, MPI_SUM,
		    MPI_COMM_WORLD);
  }

  for (auto & coef : coefsS) {
    for (auto & d : coef) 
      MPI_Allreduce(MPI_IN_PLACE, d.data(), d.size(), MPI_DOUBLE, MPI_SUM,
		    MPI_COMM_WORLD);
  }

  // Okay, now compute SVD
  //
  std::vector<Eigen::VectorXd> S;
  std::vector<Eigen::MatrixXd> U;

  for (int mm=0; mm<=mmax; mm++) {

    Eigen::JacobiSVD<Eigen::MatrixXd>
      svd(D[mm], Eigen::ComputeThinU | Eigen::ComputeThinV);

    if (myid==0)
      std::cout << "Singular values for m=" << mm << std::endl
		<< svd.singularValues() << std::endl;

    S.push_back(svd.singularValues());
    auto u = svd.matrixU();
    U.push_back(u);


#ifdef HAVE_LIBPNGPP
    if (PNG and myid==0) {

      const int minSize = 600;
      int ndup = 1;

      if (Nord < minSize) ndup = minSize/Nord   + 1;

      png::image< png::rgb_pixel > image(Nord*ndup, Nord*ndup);
      ColorGradient color;
      if (vm.count("gray"))
	color.createGrayGradient();
      else if (vm.count("5col"))
	color.createFiveColorHeatMapGradient();
      else
      	color.createSevenColorHeatMapGradient();

      std::ostringstream sout;
      sout << runtag << "_" << prefix << "_disp." << mm;

      double minV = std::numeric_limits<double>::max();
      double maxV = std::numeric_limits<double>::min();

      for (int n1=0; n1<Nord; n1++) {
	for (int n2=0; n2<Nord; n2++) {
	  minV = std::min<double>(minV, D[mm](n1, n2));
	  maxV = std::max<double>(maxV, D[mm](n1, n2));
	}
      }

      for (int n1=0; n1<Nord; n1++) {
	for (int n2=0; n2<Nord; n2++) {
	  png::rgb_pixel cval = color( (D[mm](n1, n2) - minV)/(maxV - minV) );
	  for (size_t xx = n1*ndup; xx < (n1+1)*ndup; xx++) {
	    for (size_t yy = n2*ndup; yy < (n2+1)*ndup; yy++) {
	      image[xx][yy] = cval;
	    }
	  }
	}
      }
      
      image.write(sout.str() + ".png");
      
      sout.str("");
      sout << runtag << "_" << prefix << "_ev." << mm;

      minV = -1.0;
      maxV =  1.0;

      for (int n1=0; n1<Nord; n1++) {
	for (int n2=0; n2<Nord; n2++) {
	  png::rgb_pixel cval = color( (u(n1, n2) - minV)/(maxV - minV) );
	    for (size_t xx = n1*ndup; xx < (n1+1)*ndup; xx++) {
	      for (size_t yy = n2*ndup; yy < (n2+1)*ndup; yy++) {
		image[xx][yy] = cval;
	    }
	  }
	}
      }
      
      image.write(sout.str() + ".png");
    }
#endif
  }

  // Rotate coefficients
  //

  for (int t=0; t<times.size(); t++) {

    CoefArray totT(mmax + 1);	// Per snapshot storage for
    CoefArray totC(mmax + 1);	// coefficients
    CoefArray totS(mmax + 1);

    for (auto & v : totT) v.resize(Nord);
    for (auto & v : totC) v.resize(Nord);
    for (auto & v : totS) v.resize(Nord);


    for (int mm=0; mm<=mmax; mm++) {
      totC[mm] = U[mm].transpose() * coefsC[t][mm];
      totS[mm] = Eigen::VectorXd::Zero(Nord);
      if (mm) totS[mm] = U[mm].transpose() * coefsS[t][mm];

      Eigen::VectorXd coefT(Nord), coefR(Nord);
      for (int nn=0; nn<Nord; nn++) {
	totT[mm][nn] = coefsC[t][mm][nn]*coefsC[t][mm][nn];
	if (mm) totT[mm][nn] += coefsS[t][mm][nn]*coefsS[t][mm][nn];
      }
    }

    coefsT.push_back(totT);
    coefRC.push_back(totC);
    coefRS.push_back(totS);
  }


  // Make a coefficient file for rotating coefficients in the same
  // format as readcoefs
  //
  if (myid==0) {
    std::ofstream out(runtag + "_" + prefix + ".coefs");
    std::ofstream org(runtag + "_" + prefix + ".coefs_orig");
    int ntimes = times.size();
    for (int t=0; t<ntimes; t++) {
      for (int mm=0; mm<=mmax; mm++) {
	out << std::setw(18) << times[t] << std::setw(5) << mm;
	org << std::setw(18) << times[t] << std::setw(5) << mm;
	
	for (int nn=0; nn<Nord; nn++) {
	  out << std::setw(18) << sqrt( coefRC[t][mm][nn]*coefRC[t][mm][nn] +
					coefRS[t][mm][nn]*coefRS[t][mm][nn] );
	  org << std::setw(18) << sqrt( coefsT[t][mm][nn] );
	}
	out << std::endl;
	org << std::endl;
      }
    }
  }

  // Snapshot loop
  //
  if (mbeg>=0) {	     // Print out profiles for every desired m
    for (int m=mbeg; m<=std::min<int>(mend, mmax); m++) {
      for (int t=0; t<times.size(); t++) // and every snapshot:
	write_output(ortho, t, m, Nmin, Nord, RMAX, OUTR, prefix,
		     coefRC, coefRS);
    }
  }

  // Shutdown MPI
  //
  MPI_Finalize();
    
  // DONE
  //
  return 0;
}

