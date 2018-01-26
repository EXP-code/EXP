/*
  Separate a psp structure and make a 1-d histogram of electron
  energies in planes.  

  Trace species version.

  MDWeinberg 01/20/18
*/

using namespace std;

#include <cstdlib>

#include <algorithm>
#include <iostream>
#include <fstream>
#include <sstream>
#include <iomanip>
#include <vector>
#include <string>
#include <list>
#include <map>

#include <Species.H>

#include <StringTok.H>
#include <header.H>
#include <PSP.H>

#include <boost/program_options.hpp>

namespace po = boost::program_options;

//
// VTK stuff
//
#include <vtkSmartPointer.h>
#include <vtkFloatArray.h>
#include <vtkUnsignedCharArray.h>
#include <vtkStructuredPoints.h>
#include <vtkRectilinearGrid.h>
#include <vtkPoints.h>
#include <vtkPointData.h>
#include <vtkXMLRectilinearGridWriter.h>
#include <vtkLookupTable.h>
#include <vtkVersion.h>

using vtkRectilinearGridP       = vtkSmartPointer<vtkRectilinearGrid>;
using vtkRectilinearGridWriterP = vtkSmartPointer<vtkXMLRectilinearGridWriter>;
using vtkFloatArrayP            = vtkSmartPointer<vtkFloatArray>;

using Node    = std::vector<double>;
using Element = std::vector<int>;                           

typedef std::vector< std::vector<unsigned> > I2Vector;

void
writeGrid(const std::vector<double>& T,
	  const std::vector<double>& X,
	  const std::vector<double>& Y,
	  const I2Vector& gridI,
	  const I2Vector& gridE,
	  std::ostringstream& fileName)
{
  // Create a writer
  auto writer = vtkRectilinearGridWriterP::New();

  // Append the default extension to the file name
  fileName << "." << writer->GetDefaultFileExtension();
  writer->SetFileName((fileName.str()).c_str());

  // Set knots
  auto XX   = vtkFloatArray::New();
  auto YY   = vtkFloatArray::New();
  auto ZZ   = vtkFloatArray::New();

  auto tims = vtkFloatArray::New();
  auto ions = vtkFloatArray::New();
  auto elec = vtkFloatArray::New();

  XX   -> SetName("Position");
  YY   -> SetName("Energy");
  tims -> SetName("Times");
  ions -> SetName("Ion energy");
  elec -> SetName("Electron energy");

  float f;
  int k;

  k = 0;
  for (auto z : X) XX->InsertTuple(k++, &(f=z));

  k = 0;
  for (auto z : Y) YY->InsertTuple(k++, &(f=z));
  ZZ->InsertTuple(0, &(f=0));

  k = 0;
  for (auto z : T) YY->InsertTuple(k++, &(f=z));
  tims->InsertTuple(0, &(f=0));

  // Create a pointer to a VTK Unstructured Grid data set
  auto dataSet = vtkRectilinearGridP::New();

  dataSet->SetDimensions(X.size(), Y.size(), 1);
  dataSet->SetXCoordinates (XX);
  dataSet->SetYCoordinates (YY);
  dataSet->SetZCoordinates (ZZ);

  // Insert grid data
  //
  for (size_t i=0; i<X.size(); i++) {
    for (size_t j=0; j<Y.size(); j++) {
      vtkIdType n = dataSet->FindPoint(X[i], Y[j], 0.0);

      ions->InsertTuple(n, &(f=gridI[i][j]));
      elec->InsertTuple(n, &(f=gridE[i][j]));
    }
  }

  // Add fields
  dataSet->GetPointData()->AddArray(ions);
  dataSet->GetPointData()->AddArray(elec);

  // Remove unused memory
  dataSet->Squeeze();

  // Write the data
#if VTK_MAJOR_VERSION <= 5
  writer->SetInput(dataSet);
#else
  writer->SetInputData(dataSet);
#endif
  writer->SetDataModeToAscii();
  writer->Write();
}

				// Globals for exputil library
				// Unused here
int myid = 0;
char threading_on = 0;
pthread_mutex_t mem_lock;
string outdir, runtag;

int
main(int ac, char **av)
{
  char *prog = av[0];
  double time, Emin, Emax, dE, Xmin, Xmax, dX, Lunit, Tunit;
  bool verbose = false;
  std::string cname, oname;
  int comp, sindx, eindx, hindx, dim;

  // Parse command line

  po::options_description desc("Allowed options");
  desc.add_options()
    ("help,h",		"produce help message")
    ("verbose,v",       "verbose output")
    ("time,t",		po::value<double>(&time)->default_value(1.0e20),
     "find closest time slice to requested value")
    ("Lunit,L",		po::value<double>(&Lunit)->default_value(3.086e18),
     "System length in physical units (cgs)")
    ("Tunit,T",		po::value<double>(&Tunit)->default_value(3.15569e12),
     "System time in physical units (cgs)")
    ("Emin",		po::value<double>(&Emin)->default_value(0.0),
     "Mininum energy in eV")
    ("Emax",		po::value<double>(&Emax)->default_value(100.0),
     "Maximum energy in eV")
    ("deltaE",	po::value<double>(&dE)->default_value(0.5),
     "Bin size in eV")
    ("Xmin",		po::value<double>(&Xmin)->default_value(0.0),
     "Mininum position")
    ("Xmax",		po::value<double>(&Xmax)->default_value(1.0),
     "Maximum position")
    ("deltaX",		po::value<double>(&dX)->default_value(0.1),
     "Bin size in length")
    ("species,s",	po::value<int>(&sindx)->default_value(0),
     "position of species index")
    ("electrons,e",	po::value<int>(&eindx)->default_value(10),
     "position of electron index")
    ("dim,d",		po::value<int>(&dim)->default_value(0),
     "dimension of inhomogeity (x=0, y=1, z=2)")
    ("name,c",	        po::value<std::string>(&cname)->default_value("gas"),
     "component name")
    ("files,f",         po::value< std::vector<std::string> >(), 
     "input files")
    ("output,o",	po::value<std::string>(&oname)->default_value("out"), 
     "VTK output file")
    ;


  po::variables_map vm;

  try {
    po::store(po::parse_command_line(ac, av, desc), vm);
    po::notify(vm);    
  } catch (po::error& e) {
    std::cout << "Option error: " << e.what() << std::endl;
    exit(-1);
  }

  if (vm.count("help")) {
    std::cout << desc << std::endl;
    std::cout << "Example: " << std::endl;
    std::cout << "\t" << av[0]
	      << " -E 300 -n 100 -f OUT.run.00001" << std::endl;
    return 1;
  }

  if (vm.count("verbose")) {
    verbose = true;
  }

  // Sanity check
  dim = std::max<int>(0, std::min<int>(2, dim));

  // Units
  //
  const double amu = 1.66053892e-24; // atomic mass unit in g
  const double eV  = 1.60217653e-12; // erg per eV
  double Vunit     = Lunit/Tunit;
  double KEfac     = 0.5 * amu/eV * Vunit*Vunit;

  const std::vector<double> atomic_mass = {0.000549,  // 0  electron
					   1.00797,   // 1  H
					   4.00260,   // 2  He
					   6.941,     // 3  Li
					   9.01218,   // 4  Be
					   10.81,     // 5  B
					   12.011,    // 6  C
					   14.0067,   // 7  N
					   15.9994,   // 8  O
					   18.998403, // 9  F
					   20.179,    // 10 Ne
					   22.98977,  // 11 Na
					   24.305,    // 12 Mg
					   26.98154,  // 13 Al
					   28.0855,   // 14 Si
					   30.97376,  // 15 P
					   32.06,     // 16 S
					   35.453,    // 17 Cl
					   39.948,    // 18 Ar
					   39.0983,   // 19 K
					   40.08,     // 20 Ca
					   44.9559,   // 21 Sc
					   47.90,     // 22 Ti
					   50.9415,   // 23 V
					   51.996,    // 24 Cr
					   54.9380,   // 25 Mn
					   55.847,    // 26 Fe
					   58.9332,   // 27 Co
					   58.70,     // 28 Ni
					   63.546,    // 29 Cu
					   65.38 };   // 30 Zn

  // Molecular weight
  //
  double Xf = 0.76;
  double Yf = 0.24;
  double mu = 1.0/(Xf/atomic_mass[0] + Yf/atomic_mass[2]);

  // Compute grid parameters and set up structures
  //
  int nEbin = floor((Emax - Emin)/dE+1.0e-8*(Emax - Emin));
  Emax = Emin + dE*nEbin;

  std::vector<double> E(nEbin);
  for (int i=0; i<nEbin; i++) E[i] = Emin + dE*(0.5*i);

  int nLbin = floor((Xmax - Xmin)/dX+1.0e-8*(Xmax - Xmin));
  Xmax = Xmin + dE*nLbin;

  std::vector<double> L(nLbin);
  for (int i=0; i<nLbin; i++) L[i] = Xmin + dX*(0.5*i);

  std::vector<double> T;

  I2Vector Eion(nLbin), Eelc(nLbin);
  for (int n=0; n<nLbin; n++) {
    Eion[n].resize(nEbin, 0);
    Eelc[n].resize(nEbin, 0);
  }

  // Get file arguments
  //
  std::vector<std::string> files = vm["files"].as< std::vector<std::string> >();

  bool first = true;
  
  for (auto file : files ) {

    ifstream *in = new ifstream(file.c_str());
    if (!*in) {
      cerr << "Error opening file <" << file << "> for input\n";
      exit(-1);
    }

    if (verbose) cerr << "Using filename: " << file << endl;


				// Parse the PSP file
				// ------------------
    PSPDump psp(in);

    in->close();

				// Now write a summary
				// -------------------
    if (verbose) {
      
      psp.PrintSummary(in, cerr);
    
      cerr << "\nBest fit dump to <" << time << "> has time <" 
	   << psp.SetTime(time) << ">\n";
    } else 
      psp.SetTime(time);

				// Reopen file for data input
				// --------------------------
    delete in;
    in = new ifstream(file);
    
  
    vector<double> pos(3), vel(3);

    PSPstanza *stanza;
    SParticle* part;

    for (stanza=psp.GetStanza(); stanza!=0; stanza=psp.NextStanza()) {
    
      if (stanza->name != cname) continue;

      T.push_back(psp.CurrentTime());

				// Position to beginning of particles
				// ----------------------------------
      in->seekg(stanza->pspos);

      for (part=psp.GetParticle(in); part!=0; part=psp.NextParticle(in)) {

	int Pindx = floor( (part->pos(dim) - Xmin)/dX );

	if (part->pos(dim) < Xmin or part->pos(dim) > Xmax) continue;
	if (Pindx < 0 or Pindx >= nLbin) continue;

	double kEe = 0.0, kEi = 0.0;
	for (size_t i=0; i<3; i++) {
	  double ve = part->datr(eindx+i);
	  kEe += ve*ve;
	  double vi = part->vel(i);
	  kEi += vi*vi;
	}
	KeyConvert kc(part->iatr(sindx));

	kEe *= KEfac * atomic_mass[0];
	kEi *= KEfac * mu;

	if (kEe >= Emin and kEe < Emax) {
	  int Eindx = floor( (kEe - Emin)/dE );
	  if (Eindx >= 0 and Eindx < nEbin) Eelc[Pindx][Eindx]++;
	}

	if (kEi >= Emin and kEi < Emax) {
	  int Eindx = floor( (kEi - Emin)/dE );
	  if (Eindx >= 0 and Eindx < nEbin) Eion[Pindx][Eindx]++;
	}

      } // END: Particle loop

    } // END: stanza lop

  } // END: file loop

  
  // Write the VTK file
  //
  std::ostringstream fileName;

  fileName << oname;

  writeGrid(T, L, E, Eion, Eelc, fileName);
  std::cout << "Wrote file <" << fileName.str() << ">" << std::endl;

  return 0;
}
