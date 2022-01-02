/*
  Convert PSP to points

  MDWeinberg 03/13/20
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
#include <cxxopts.H>		// Option parsing
#include <libvars.H>		// EXP library globals
#include <header.H>
#include <writePVD.H>
#include <PSP.H>

//
// VTK stuff
//
#include <vtkSmartPointer.h>
#include <vtkDoubleArray.h>
#include <vtkFloatArray.h>
#include <vtkIntArray.h>
#include <vtkUnsignedCharArray.h>
#include <vtkUnstructuredGrid.h>
#include <vtkPoints.h>
#include <vtkPointData.h>
#include <vtkCellData.h>
#include <vtkCellArray.h>
#include <vtkXMLUnstructuredGridWriter.h>
#include <vtkLookupTable.h>
#include <vtkVersion.h>

using vtkFloatArrayP            = vtkSmartPointer<vtkFloatArray>;

// KD tree for density computation
//
#include <KDtree.H>

typedef std::vector< std::vector<unsigned> > I2Vector;

int
main(int ac, char **av)
{
  char *prog = av[0];
  int ibeg, iend, istride, Ndens, iskip;
  double time, Emin, Emax, dE, Xmin, Xmax, dX, Lunit, Tunit;
  bool verbose = false, logE = false, PVD = false;
  std::string cname, rtag;
  int comp, sindx, eindx, hindx, dim;

  // Parse command line
  //
  cxxopts::Options options(prog, "Compute a VTK point file with optional density computation from a PSP file\n");

  options.add_options()
   ("h,help", "produce help message")
   ("v,verbose", "verbose output")
   ("OUT", "assume that PSP files are in original format")
   ("SPL", "assume that PSP files are in split format")
   ("P,PVD", "create a ParaView PVD file")
   ("c,name", "component name",
     cxxopts::value<std::string>(cname)->default_value("gas"))
   ("t,rtag", "runtag name",
     cxxopts::value<std::string>(rtag)->default_value("run"))
   ("1,begin", "initial sequence counter",
     cxxopts::value<int>(ibeg)->default_value("0"))
   ("2,final", "final sequence counter",
     cxxopts::value<int>(iend)->default_value("1000000"))
   ("s,stride", "sequence counter stride",
     cxxopts::value<int>(istride)->default_value("1"))
   ("S,skip", "particle number stride to reduce particle counts for rendering",
     cxxopts::value<int>(iskip)->default_value("1"))
   ("N,Ndens", "KD density estimate count",
     cxxopts::value<int>(Ndens)->default_value("0"))
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

  std::vector<double> times;
  std::vector<std::string> outfiles;

  bool first = true;
  
  int C = 0;

  for (int n=ibeg; n<=iend; n+=istride) {

    std::ostringstream file;
    if (vm.count("OUT")) 
      file << "OUT." << rtag << "." << std::setfill('0') << std::setw(5) << n;
    else
      file << "SPL." << rtag << "." << std::setfill('0') << std::setw(5) << n;
    
    std::ifstream in(file.str());
    if (!in) {
      std::cerr << "Error opening file <" << file.str() << "> for input."
		<< std::endl
		<< "Assuming end of sequence . . . finalizing."
		<< std::endl;
      break;
    }
    in.close();

    if (verbose) std::cerr << "Using filename: " << file.str() << std::endl;
    else         std::cout << "Begin file " << file.str()
			   << " . . . " << std::flush;


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

    std::ostringstream fileName;

    for (stanza=psp->GetStanza(); stanza!=0; stanza=psp->NextStanza()) {
    
      if (stanza->name != cname) continue;

      /*
      std::cout << setw(15) << stanza->comp.nbod 
		<< setw(10) << stanza->comp.niatr 
		<< setw(10) << stanza->comp.ndatr 
		<< endl;
      */

      vtkSmartPointer<vtkUnstructuredGrid> uGrid =
	vtkSmartPointer<vtkUnstructuredGrid>::New();

      vtkSmartPointer<vtkPoints> pos = vtkSmartPointer<vtkPoints>::New();

      vtkFloatArrayP vel = vtkFloatArrayP::New();
      vel->SetNumberOfComponents(3);
      vel->SetName("velocities");

      vtkFloatArrayP mas = vtkFloatArrayP::New();
      mas->SetNumberOfComponents(1);
      mas->SetName("masses");

      std::vector<vtkFloatArrayP> attrib;
      if (stanza->comp.ndatr) {
	attrib.resize(stanza->comp.ndatr);
	for (int i=0; i<stanza->comp.ndatr; i++) {
	  attrib[i] = vtkFloatArrayP::New();
	  attrib[i]->SetNumberOfComponents(1);
	  std::ostringstream sout; sout << "Field " << i+1;
	  attrib[i]->SetName(sout.str().c_str());
	}
      }

      vtkSmartPointer<vtkCellArray> conn = vtkSmartPointer<vtkCellArray>::New();
    

      // Compute density based on N-pt balls
      //
      vtkFloatArrayP dens = vtkFloatArrayP::New();
      dens->SetNumberOfComponents(1);
      dens->SetName("density");

      if (Ndens) {
	typedef point <double, 3> point3;
	typedef kdtree<double, 3> tree3;

	std::vector<double> mass;
	std::vector<point3> points;

	for (part=psp->GetParticle(); part!=0; part=psp->NextParticle()) {
	  points.push_back({part->pos(0), part->pos(1), part->pos(2)});
	  mass.push_back(part->mass());
	}
	
	tree3 tree(points.begin(), points.end());

	for (int k=0; k<points.size(); k++) {
	  // Stride the density computation
	  if (k % iskip == 0) {
	    auto ret = tree.nearestN(points[k], Ndens);

	    double volume = 4.0*M_PI/3.0*std::pow(std::get<2>(ret), 3.0);
	    if (volume>0.0)
	      dens->InsertNextValue(std::get<1>(ret)/volume);
	    else
	      dens->InsertNextValue(1.0e-18);
	  }
	}
      }

      int NP = 0;
      for (part=psp->GetParticle(); part!=0; part=psp->NextParticle()) {
	
	// Stride the particle insertion to match the density evaluation
	if (NP++ % iskip == 0) {

	  float m = part->mass();
	  mas->InsertNextValue(m);
	  
	  float x = part->pos(0);
	  float y = part->pos(1);
	  float z = part->pos(2);
	  
	  pos->InsertNextPoint(x, y, z);

	  float u = part->vel(0);
	  float v = part->vel(1);
	  float w = part->vel(2);
	  
	  vel->InsertNextTuple3(u, v, w);
	}

      }

      // Create topology of point cloud. This is necessary!
      //
      int cnt = 0;
      for (int i=0; i<stanza->comp.nbod; i++) {
	if (i % iskip==0) {
	  vtkIdType c = cnt++;
	  conn->InsertNextCell(1, &c);
	}
      }
    
      // Add the point locations
      //
      uGrid->SetPoints(pos);

      // Add the particle masses
      //
      uGrid->GetPointData()->AddArray(mas);

      // Add density
      //
      if (dens->GetNumberOfTuples()) 
	uGrid->GetPointData()->AddArray(dens);

      // Add the atrribute fields
      //
      for (auto f : attrib)
	uGrid->GetPointData()->AddArray(f);

      // Add the velocity vectors
      //
      uGrid->SetCells(VTK_VERTEX, conn);
      uGrid->GetCellData()->SetVectors(vel);

      // Write the data to file
      //
      vtkSmartPointer<vtkXMLUnstructuredGridWriter> writer =
	vtkSmartPointer<vtkXMLUnstructuredGridWriter>::New();

      // Append the default extension to the file name
      //
      fileName << rtag << "_" << std::setw(5) << std::setfill('0') << n
	       << "." << writer->GetDefaultFileExtension();
      writer->SetFileName((fileName.str()).c_str());
      writer->SetInputData(uGrid);
      writer->Write();

    } // END: PSP

    if (PVD) {
      times.push_back(T);
      outfiles.push_back(fileName.str());
    }

    if (not verbose) std::cout << "done" << std::endl;

    C++;

  } // END: file loop

  // Create PVD file
  //
  if (PVD) {
    writePVD(rtag+".pvd", times, outfiles);
  }

  return 0;
}
