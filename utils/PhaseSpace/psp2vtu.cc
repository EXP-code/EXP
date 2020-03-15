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
#include <header.H>
#include <PSP.H>

#include <boost/program_options.hpp>
#include <boost/property_tree/ptree.hpp>
#include <boost/property_tree/xml_parser.hpp> 

namespace po = boost::program_options;
namespace pt = boost::property_tree;

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

void writePVD(const std::string& filename,
	      const std::vector<double>& times,
	      const std::vector<std::string>& files)
{
  // Sanity check
  //
  if (times.size() != files.size()) {
    std::cerr << "Mismatch in file and time arrays" << std::endl;
    exit(-3);
  }

  // Make file collection elements
  //
  pt::ptree ptC;

  for (size_t i=0; i<times.size(); i++) {
    boost::property_tree::ptree x;
    x.put("<xmlattr>.timestep", times[i]);
    x.put("<xmlattr>.part", 0);
    x.put("<xmlattr>.file", files[i]);

    ptC.add_child("DataSet", x);
  }

  // Add VTKFile attributes
  //
  pt::ptree ptP;
  
  ptP.put("<xmlattr>.type", "Collection");
  ptP.put("<xmlattr>.version", "0.1");
  ptP.put("<xmlattr>.byte_order", "LittleEndian");
  ptP.put("<xmlattr>.compressor", "vtkZLibDataCompressor");
  ptP.add_child("Collection", ptC);
  
  // Make the top-level property tree
  //
  pt::ptree PT;

  PT.add_child("VTKFile", ptP);

  // Write the property tree to the XML file.
  //
  pt::xml_parser::write_xml(filename.c_str(), PT, std::locale(), pt::xml_writer_make_settings<std::string>(' ', 4));

  std::cout << "Wrote PVD file <" << filename.c_str() << "> "
	    << " with " << times.size() << " data sets." << std::endl;
}

				// Globals for exputil library
				// Unused here
int myid = 0;
char threading_on = 0;
pthread_mutex_t mem_lock;
string outdir, runtag;
int multistep = 1;

int
main(int ac, char **av)
{
  char *prog = av[0];
  int ibeg, iend, istride, Ndens;
  double time, Emin, Emax, dE, Xmin, Xmax, dX, Lunit, Tunit;
  bool verbose = false, logE = false, PVD = false;
  std::string cname, rtag;
  int comp, sindx, eindx, hindx, dim;

  // Parse command line

  po::options_description desc("Allowed options");
  desc.add_options()
    ("help,h",		"produce help message")
    ("verbose,v",       "verbose output")
    ("PVD,P",		"create a ParaView PVD file")
    ("name,c",	        po::value<std::string>(&cname)->default_value("gas"),
     "component name")
    ("rtag,t",		po::value<std::string>(&rtag)->default_value("run"), 
     "runtag name")
    ("begin,1",		po::value<int>(&ibeg)->default_value(0),
     "initial sequence counter")
    ("final,2",		po::value<int>(&iend)->default_value(1000000),
     "final sequence counter")
    ("stride,s",	po::value<int>(&istride)->default_value(1),
     "sequence counter stride")
    ("Ndens,N",		po::value<int>(&Ndens)->default_value(0),
     "KD density estimate count")
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

  std::vector<double> times;
  std::vector<std::string> outfiles;

  bool first = true;
  
  int C = 0;

  for (int n=ibeg; n<=iend; n+=istride) {

    std::ostringstream file;
    file << "OUT." << rtag << "." << std::setfill('0') << std::setw(5) << n;
    
    ifstream *in = new ifstream(file.str());
    if (!*in) {
      cerr << "Error opening file <" << file.str() << "> for input."
	   << std::endl
	   << "Assuming end of sequence . . . continuing."
	   << std::endl;
      break;
    }

    if (verbose) cerr << "Using filename: " << file.str() << endl;


				// Parse the PSP file
				// ------------------
    PSPDump psp(in);

    in->close();

				// Now write a summary
				// -------------------
    if (verbose) {
      
      psp.PrintSummary(in, cerr);
    
      std::cerr << std::endl
		<< "Best fit dump to <" << time << "> has time <" 
		<< psp.SetTime(time) << ">" << std::endl;
    } else 
      psp.SetTime(time);

				// Reopen file for data input
				// --------------------------
    delete in;
    in = new ifstream(file.str());
    
  
    vector<double> pos(3), vel(3);

    PSPstanza *stanza;
    SParticle* part;

    double T = psp.CurrentTime();
    if (verbose) {
      std::cerr << std::endl << "PSP time is <" << T << ">" << std::endl;
    }

    std::ostringstream fileName;

    for (stanza=psp.GetStanza(); stanza!=0; stanza=psp.NextStanza()) {
    
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
      mas->SetNumberOfComponents(1);
      mas->SetName("density");

      if (Ndens) {
	in->seekg(stanza->pspos); // Move to beginning of particles

	typedef point <double, 3> point3;
	typedef kdtree<double, 3> tree3;

	std::vector<double> mass;
	std::vector<point3> points;

	for (part=psp.GetParticle(in); part!=0; part=psp.NextParticle(in)) {
	  points.push_back({part->pos(0), part->pos(1), part->pos(2)});
	  mass.push_back(part->mass());
	}
	  
	tree3 tree(points.begin(), points.end());

	for (int k=0; k<points.size(); k++) {
	  auto ret = tree.nearestN(points[k], Ndens);
	    
	  double volume = 4.0*M_PI/3.0*std::pow(std::get<2>(ret), 3.0);
	  if (volume>0.0)
	    dens->InsertNextValue(mass[k]*Ndens/volume);
	  else
	    dens->InsertNextValue(1.0e-18);
	}
      }
				// Position to beginning of particles
				// ----------------------------------
      in->seekg(stanza->pspos);

      for (part=psp.GetParticle(in); part!=0; part=psp.NextParticle(in)) {

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

      // Create topology of point cloud. This is necessary!
      //
      for (int i=0; i<stanza->comp.nbod; i++) {
	vtkIdType c = i;
	conn->InsertNextCell(1, &c);
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
      fileName << rtag << "_" << std::setw(5) << std::setfill('0') << C
	       << "." << writer->GetDefaultFileExtension();
      writer->SetFileName((fileName.str()).c_str());
      writer->SetInputData(uGrid);
      writer->Write();

    } // END: PSP

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
