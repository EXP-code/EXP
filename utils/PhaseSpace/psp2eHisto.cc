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
#include <locale>

#include <Species.H>

#include <StringTok.H>
#include <writePVD.H>
#include <cxxopts.H>
#include <header.H>
#include <PSP.H>

//
// VTK stuff
//
#include <vtkSmartPointer.h>
#include <vtkDoubleArray.h>
#include <vtkFloatArray.h>
#include <vtkIntArray.h>
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

typedef std::vector< std::vector<unsigned> > I2Vector;


// Add "time" value to a VTK dataset.
void
AddTimeToVTK(vtkDataSet *ds, double time)
{
  vtkDoubleArray *t = vtkDoubleArray::New();
  t->SetName("TIME");
  t->SetNumberOfTuples(1);
  t->SetTuple1(0, time);
  ds->GetFieldData()->AddArray(t);
}
 
// Add "cycle" value to a VTK dataset.
void
AddCycleToVTK(vtkDataSet *ds, int cycle)
{
  vtkIntArray *c = vtkIntArray::New();
  c->SetName("CYCLE");
  c->SetNumberOfTuples(1);
  c->SetTuple1(0, cycle);
  ds->GetFieldData()->AddArray(c);
}

void
writeGrid(const double T, const int C,
	  const std::vector<double>& X,
	  const std::vector<double>& Y,
	  const I2Vector& gridI,
	  const I2Vector& gridE,
	  std::ostringstream& fileName, bool logE)
{
  // Create a writer
  auto writer = vtkRectilinearGridWriterP::New();

  // Append the default extension to the file name
  fileName << "_" << std::setw(5) << std::setfill('0') << C
	   << "." << writer->GetDefaultFileExtension();
  writer->SetFileName((fileName.str()).c_str());

  // Set knots
  auto XX   = vtkFloatArray::New();
  auto YY   = vtkFloatArray::New();
  auto ZZ   = vtkFloatArray::New();

  auto tims = vtkFloatArray::New();
  auto ions = vtkFloatArray::New();
  auto elec = vtkFloatArray::New();

  XX   -> SetName("Position");
  if (logE) YY-> SetName("Log Energy");
  else      YY-> SetName("Energy");
  tims -> SetName("Times");
  ions -> SetName("Ion counts");
  elec -> SetName("Electron counts");

  float f;
  int k;

  k = 0;
  for (auto z : X) XX->InsertTuple(k++, &(f=z));

  k = 0;
  for (auto z : Y) YY->InsertTuple(k++, &(f=z));

  ZZ->InsertTuple(0, &(f=0));

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
      float x = X[i], y = Y[j];
      vtkIdType n = dataSet->FindPoint(x, y, 0.0);

      if (n>=0) {
	ions->InsertTuple(n, &(f=gridI[i][j]));
	elec->InsertTuple(n, &(f=gridE[i][j]));
      } else {
	std::cout << "Could not find point at (" << X[i] << ", " << Y[j] << ")"
		  << std::endl;
      }
    }
  }

  // Add fields
  dataSet->GetPointData()->AddArray(ions);
  dataSet->GetPointData()->AddArray(elec);

  AddTimeToVTK (dataSet, T);
  AddCycleToVTK(dataSet, C);

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


int
main(int ac, char **av)
{
  char *prog = av[0];
  int ibeg, iend, istride;
  double time, Emin, Emax, dE, Xmin, Xmax, dX, Lunit, Tunit;
  bool verbose = false, logE = false, PVD = false;
  std::string cname, rtag;
  int comp, sindx, eindx, hindx, dim;

  // Parse command line
  //
  cxxopts::Options options(prog, "Separate a psp structure and make a 1-d histogram of electron energies in planes.\nTrace species version.");

  options.add_options()
    ("h,help", "produce help message")
    ("v,verbose", "verbose output")
    ("logE", "bin logarithmically in energy")
    ("PVD", "create a ParaView PVD file")
    ("OUT", "assume that PSP files are in original format")
    ("SPL", "assume that PSP files are in split format")
    ("L,Lunit", "System length in physical units (cgs)",
     cxxopts::value<double>(Lunit)->default_value("3.086e18"))
    ("T,Tunit", "System time in physical units (cgs)",
     cxxopts::value<double>(Tunit)->default_value("3.15569e10"))
    ("Emin", "Mininum energy in eV",
     cxxopts::value<double>(Emin)->default_value("0.0"))
    ("Emax", "Maximum energy in eV",
     cxxopts::value<double>(Emax)->default_value("100.0"))
    ("deltaE", "Bin size in eV",
     cxxopts::value<double>(dE)->default_value("0.5"))
    ("Xmin", "Mininum position",
     cxxopts::value<double>(Xmin)->default_value("0.0"))
    ("Xmax", "Maximum position",
     cxxopts::value<double>(Xmax)->default_value("1.0"))
    ("deltaX", "Bin size in length",
     cxxopts::value<double>(dX)->default_value("0.1"))
    ("s,species", "position of species index",
     cxxopts::value<int>(sindx)->default_value("0"))
    ("e,electrons", "position of electron index",
     cxxopts::value<int>(eindx)->default_value("10"))
    ("d,dim", "dimension of inhomogeity (x=0, y=1, z=2)",
     cxxopts::value<int>(dim)->default_value("0"))
    ("c,name", "component name",
     cxxopts::value<std::string>(cname)->default_value("gas"))
    ("t,rtag", "runtag name",
     cxxopts::value<std::string>(rtag)->default_value("run"))
    ("begin", "initial sequence counter",
     cxxopts::value<int>(ibeg)->default_value("0"))
    ("final", "final sequence counter",
     cxxopts::value<int>(iend)->default_value("1000000"))
    ("stride", "sequence counter stride",
     cxxopts::value<int>(istride)->default_value("1"))
    ;

  cxxopts::ParseResult vm;

  try {
    vm = options.parse(ac, av);
  } catch (cxxopts::OptionException& e) {
    std::cout << "Option error: " << e.what() << std::endl;
    exit(-1);
  }

  if (vm.count("help")) {
    std::cout << options.help() << std::endl;
    std::cout << "Example: " << std::endl;
    std::cout << "\t" << av[0]
	      << " -E 300 -n 100 -f OUT.run.00001" << std::endl;
    return 1;
  }

  if (vm.count("verbose")) {
    verbose = true;
  }

  if (vm.count("logE")) {
    logE = true;
  }

  if (vm.count("PVD")) {
    PVD = true;
  }

  std::vector<double> times;
  std::vector<std::string> outfiles;

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
  double mu = 1.0/(Xf/atomic_mass[1] + Yf/atomic_mass[2]);

  // Compute grid parameters and set up structures
  //
  if (logE) {
    if (Emin==0.0 or Emax==0.0) {
      std::cerr << "Energy must be greater than zero for log scaling"
		<< std::endl;
      exit(-2);
    }
    Emin = log10(Emin);
    Emax = log10(Emax);
  }
  int nEbin = floor((Emax - Emin)/dE+1.0e-8*(Emax - Emin));
  Emax = Emin + dE*nEbin;

  std::vector<double> E(nEbin);
  for (int i=0; i<nEbin; i++) E[i] = Emin + dE*(0.5+i);

  int nLbin = floor((Xmax - Xmin)/dX+1.0e-8*(Xmax - Xmin));
  Xmax = Xmin + dX*nLbin;

  std::vector<double> L(nLbin);
  for (int i=0; i<nLbin; i++) L[i] = Xmin + dX*(0.5+i);

  int C=0;

  I2Vector Eion(nLbin), Eelc(nLbin);
  for (int n=0; n<nLbin; n++) {
    Eion[n].resize(nEbin);
    Eelc[n].resize(nEbin);
  }

  bool first = true;
  
  for (int n=ibeg; n<=iend; n+=istride) {

    std::ostringstream file;
    file << "OUT." << rtag << "." << std::setfill('0') << std::setw(5) << n;
    
    std::ifstream in(file.str());
    if (in) {
      cerr << "Error opening file <" << file.str() << "> for input."
	   << std::endl
	   << "Assuming end of sequence . . . continuing."
	   << std::endl;
      break;
    }
    in.close();

    if (verbose) cerr << "Using filename: " << file.str() << endl;


				// Parse the PSP file
				// ------------------
    PSPptr psp;
    if (vm.count("SPL")) psp = std::make_shared<PSPspl>(file.str());
    else                 psp = std::make_shared<PSPout>(file.str());

				// Now write a summary
				// -------------------
    if (verbose) {
      
      psp->PrintSummary(cerr);
    
      std::cerr << std::endl
		<< "Best fit dump to <" << time << "> has time <" 
		<< psp->CurrentTime() << ">" << std::endl;
    }

    vector<double> pos(3), vel(3);

    PSPstanza *stanza;
    SParticle* part;

    double T = psp->CurrentTime();
    if (verbose) {
      std::cerr << std::endl << "PSP time is <" << T << ">" << std::endl;
    }

    for (stanza=psp->GetStanza(); stanza!=0; stanza=psp->NextStanza()) {
    
      if (stanza->name != cname) continue;

      unsigned total = 0, gridded = 0, pout = 0, eEout = 0, eIout = 0;

				// Clear counters
				// --------------
      for (int n=0; n<nLbin; n++) {
	std::fill(Eion[n].begin(), Eion[n].end(), 0);
	std::fill(Eelc[n].begin(), Eelc[n].end(), 0);
      }

      for (part=psp->GetParticle(); part!=0; part=psp->NextParticle()) {

	total++;

	int Pindx = floor( (part->pos(dim) - Xmin)/dX );

	if (part->pos(dim) < Xmin or part->pos(dim) > Xmax) {
	  pout++;
	  continue;
	}

	if (Pindx < 0 or Pindx >= nLbin) {
	  pout++;
	  continue;
	}

	gridded++;

	double kEe = 0.0, kEi = 0.0;
	for (size_t i=0; i<3; i++) {
	  double ve = part->datr(eindx+i);
	  kEe += ve*ve;
	  double vi = part->vel(i);
	  kEi += vi*vi;
	}

	kEe *= KEfac * atomic_mass[0];
	kEi *= KEfac * mu;

	if (logE) {
	  kEe = log10(kEe);
	  kEi = log10(kEi);
	}

	if (kEe >= Emin and kEe < Emax) {
	  int Eindx = floor( (kEe - Emin)/dE );
	  if (Eindx >= 0 and Eindx < nEbin) Eelc[Pindx][Eindx]++;
	  else eEout++;
	}

	if (kEi >= Emin and kEi < Emax) {
	  int Eindx = floor( (kEi - Emin)/dE );
	  if (Eindx >= 0 and Eindx < nEbin) Eion[Pindx][Eindx]++;
	  else eIout++;
	}

      } // END: Particle loop

      std::cout << "File <" << file.str() << ">: "
		<< gridded << " out of " << total << " with "
		<< pout << " position oab, "
		<< eEout << " electron oab, "
		<< eIout << " ion oab" << std::endl;

    } // END: stanza loop


    // Write the VTK file
    //
    std::ostringstream fileName;

    fileName << rtag;

    writeGrid(T, C, L, E, Eion, Eelc, fileName, logE);
    std::cout << "Wrote file <" << fileName.str() << ">" << std::endl;
    
    if (PVD) {
      times.push_back(T);
      outfiles.push_back(fileName.str());
    }

    C++;

  } // END: file loop

  // Create PVD file
  //
  if (PVD) {
    writePVD(rtag+".pvd", times, outfiles);
  }

  return 0;
}
