#include <iostream>
#include <fstream>

#include <cstdlib>
#include <cmath>
#include <values.h>

//
// STL stuff
//
#include <vector>
#include <string>
#include <list>
#include <algorithm>

//
// BOOST stuff
//
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

//
// PSP stuff
//
#include <StringTok.H>
#include <header.H>
#include <PSP.H>

				// Globals for exputil library
				// Unused here
int myid = 0;
int VERBOSE = 4;
int nthrds = 1;
int this_step = 0;
unsigned multistep = 0;
unsigned maxlev = 100;
int mstep = 1;
int Mstep = 1;
vector<int> stepL(1, 0), stepN(1, 1);
char threading_on = 0;
pthread_mutex_t mem_lock;
pthread_mutex_t coef_lock;
string outdir, runtag;
double tpos = 0.0;
double tnow = 0.0;


int main(int argc, char**argv)
{
  int numx=20, numy=20, numz=20;
  double xmin=-1.0, xmax=1.0;
  double ymin=-1.0, ymax=1.0;
  double zmin=-1.0, zmax=1.0;
  double vscale = 1.0;
  string infile("OUT.bin");
  string outfile("OUT");
  string cname, dname, sname;
  double time = 1.0;
  unsigned long initial_dark = 0, final_dark = MAXLONG;
  unsigned long initial_star = 0, final_star = MAXLONG;
  unsigned long initial_gas  = 0, final_gas  = MAXLONG;
  bool mask = false;
  bool verbose = false;
				// Boltzmann constant (cgs)
  const double boltz = 1.3810e-16;
				// Hydrogen fraction
  const double f_H = 0.76;
  				// Proton mass (g)
  const double mp = 1.67262158e-24;
				// Mean mass
  const double mm = f_H*mp + (1.0-f_H)*4.0*mp;
				// Cp/Cv; isentropic expansion factor
  const double gamma = 5.0/3.0;

  double Vconv = 120*1e5;

  po::options_description desc("Allowed options");
  desc.add_options()
    ("help,h", "produce this help message")
    ("verbose,v", "verbose output")
    ("mask,b", "blank empty cells")
    ("debug", "turn on debugging output")
    ("numx,1", po::value<int>(&numx)->default_value(20), 
     "number of bins in x direction")
    ("numy,2", po::value<int>(&numy)->default_value(20), 
     "number of bins in x direction")
    ("numz,3", po::value<int>(&numz)->default_value(20), 
     "number of bins in x direction")
    ("xmin,x", po::value<double>(&xmin)->default_value(-1.0), 
     "minimum x value")
    ("xmax,X", po::value<double>(&xmax)->default_value(1.0), 
     "maximum x value")
    ("ymin,y", po::value<double>(&ymin)->default_value(-1.0), 
     "minimum y value")
    ("ymax,Y", po::value<double>(&ymax)->default_value(1.0), 
     "maximum y value")
    ("zmin,z", po::value<double>(&zmin)->default_value(-1.0), 
     "minimum z value")
    ("zmax,Z", po::value<double>(&zmax)->default_value(1.0), 
     "maximum z value")
    ("vscale,V", po::value<double>(&vscale)->default_value(1.0), 
     "vertical scale factor")
    ("time,t", po::value<double>(&time)->default_value(0.0), 
     "desired PSP time")
    ("dark-name,d", po::value<string>(&dname),
     "PSP dark component name")
    ("star-name,s", po::value<string>(&sname),
     "PSP star component name")
    ("gas-name,g", po::value<string>(&cname),
     "PSP gas component name")
    ("input,i", po::value<string>(&infile)->default_value("OUT.bin"),
     "input file name")
    ("output,o", po::value<string>(&outfile)->default_value("OUT"),
     "output file ename")
    ("initial-gas", po::value<unsigned long>(&initial_gas)->default_value(0), 
     "initial gas particle index")
    ("final-gas", po::value<unsigned long>(&final_gas)->default_value(MAXLONG), 
     "initial gas particle index")
    ("initial-star", po::value<unsigned long>(&initial_star)->default_value(0), 
     "initial star particle index")
    ("final-star", po::value<unsigned long>(&final_star)->default_value(MAXLONG), 
     "initial star particle index")
    ("initial-dark", po::value<unsigned long>(&initial_dark)->default_value(0), 
     "initial dark particle index")
    ("final-dark", po::value<unsigned long>(&final_dark)->default_value(MAXLONG), 
     "initial dark particle index")
    ;
  
  po::variables_map vm;

  try {
    po::store(po::parse_command_line(argc, argv, desc), vm);
    po::notify(vm);    
  } catch (po::error& e) {
    std::cout << "Option error: " << e.what() << std::endl;
    exit(-1);
  }

  if (vm.count("help")) {
    cout << desc << "\n";
    return 1;
  }
 
  if (vm.count("verbose")) verbose = true;

  if (vm.count("mask")) mask = true;

  ifstream *in = new ifstream(infile.c_str());
  if (!*in) {
    cerr << "Error opening file <" << infile << "> for input\n";
    exit(-1);
  }

  if (verbose) cerr << "Using filename: " << infile << endl;

				// Parse the PSP file
				// ------------------
  PSPDump psp(in);
  
				// Now write a summary
				// -------------------
  if (verbose) {

    psp.PrintSummary(in, cerr);
    
    cerr << "\nBest fit dump to <" << time << "> has time <" 
	 << psp.SetTime(time) << ">\n";
  } else 
    psp.SetTime(time);

  in->close();
  delete in;
  
				// Make the arrays
				// -----------------------------

  double dx = (xmax - xmin)/numx;
  double dy = (ymax - ymin)/numy;
  double dz = (zmax - zmin)/numz;

  cout << "Grid bounds: "
       << "[" << xmin << ", " << xmax << "] "
       << "[" << ymin << ", " << ymax << "] "
       << "[" << zmin << ", " << zmax << "] "
       << endl;
  cout << "Grid spacing: [" << dx << ", " << dy << ", " << dz << "]"
       << endl;

  vector<float> xyz(3), uvw(3), dattr;
  vector< vector< vector<float> > > mass(numx), gdens(numx), gtemp(numx);
  vector< vector< vector<float> > > gknud(numx), gstrl(numx), gmach(numx);
  vector< vector< vector<float> > > sdens(numx), ddens(numx);
  vector< vector< vector< vector<double> > > > pos(numx);
  vector< vector< vector< vector<float > > > > vel(numx);

  for (int i=0; i<numx; i++) {
    
    mass [i] = vector< vector<float> >(numy);
    gtemp[i] = vector< vector<float> >(numy);
    gdens[i] = vector< vector<float> >(numy);
    gknud[i] = vector< vector<float> >(numy);
    gstrl[i] = vector< vector<float> >(numy);
    gmach[i] = vector< vector<float> >(numy);
    sdens[i] = vector< vector<float> >(numy);
    ddens[i] = vector< vector<float> >(numy);

    pos[i]  = vector< vector< vector<double> > >(numy);
    vel[i]  = vector< vector< vector<float> > >(numy);
    
    for (int j=0; j<numy; j++) {
      
      mass[i][j]  = vector<float>(numz, 0.0);
      gtemp[i][j] = vector<float>(numz, 0.0);
      gdens[i][j] = vector<float>(numz, 0.0);
      gknud[i][j] = vector<float>(numz, 0.0);
      gstrl[i][j] = vector<float>(numz, 0.0);
      gmach[i][j] = vector<float>(numz, 0.0);
      sdens[i][j] = vector<float>(numz, 0.0);
      ddens[i][j] = vector<float>(numz, 0.0);

      pos[i][j]   = vector< vector<double> >(numz);
      vel[i][j]   = vector< vector<float > >(numz);

      for (int k=0; k<numz; k++) {
	pos[i][j][k] = vector<double>(3);
	vel[i][j][k] = vector<float> (3, 0.0);

	pos[i][j][k][0] = xmin + dx*(0.5 + i);
	pos[i][j][k][1] = ymin + dy*(0.5 + j);
	pos[i][j][k][2] = zmin + dz*(0.5 + k);
	pos[i][j][k][2] *= vscale;
      }
    }
  }
				// Reopen file to get component
				// -----------------------------
  in = new ifstream(infile.c_str());

  list<PSPstanza>::iterator its;
  double rtmp, ms, ps[3], vs[3];
  size_t indx;
  int itmp;

  bool found_gas  = false;
  bool found_star = false;
  bool found_dark = false;

  bool btemp = false, bdens = false, bknud = false, bstrl = false;

  vtkSmartPointer<vtkPoints>    part  = vtkPoints    ::New();
  vtkSmartPointer<vtkDataArray> dens  = vtkFloatArray::New();
  vtkSmartPointer<vtkDataArray> temp  = vtkFloatArray::New();
  vtkSmartPointer<vtkDataArray> knud  = vtkFloatArray::New();
  vtkSmartPointer<vtkDataArray> strl  = vtkFloatArray::New();
  vtkSmartPointer<vtkDataArray> mach  = vtkFloatArray::New();
  vtkSmartPointer<vtkDataArray> velo  = vtkFloatArray::New();
  
  dens->SetName("density");
  temp->SetName("temperature");
  knud->SetName("Knudsen");
  strl->SetName("Strouhal");
  mach->SetName("Mach");

  velo->SetName("velocity");
  velo->SetNumberOfComponents(3);

  int offset = 0;
  
  for (its  = psp.CurrentDump()->stanzas.begin(); 
       its != psp.CurrentDump()->stanzas.end();  its++) {

    if (dname.compare(its->name) == 0) {

      found_dark = true;

      //
      // Position to beginning of particles
      //
      in->seekg(its->pspos);
      
      for (int j=0; j<its->nbod; j++) {
	if (its->index_size) in->read((char *)&indx, its->index_size);
	else                 indx = j;

	in->read((char *)&ms, sizeof(double));
	for (int i=0; i<3; i++) {
	  in->read((char *)&ps[i], sizeof(double));
	  xyz[i] = ps[i];
	}
	for (int i=0; i<3; i++) {
	  in->read((char *)&vs[i], sizeof(double));
	  uvw[i] = vs[i];
	}
	in->read((char *)&rtmp, sizeof(double));
	for (int i=0; i<its->niatr; i++) {
	  in->read((char *)&itmp, sizeof(int));
	}
	dattr.clear();
	for (int i=0; i<its->ndatr; i++) {
	  in->read((char *)&rtmp, sizeof(double));
	  dattr.push_back(rtmp);
	}
				// Accumulate
				// 
	if (indx > initial_dark && indx <= final_dark &&
	    ps[0] >= xmin && ps[0] < xmax       &&
	    ps[1] >= ymin && ps[1] < ymax       &&
	    ps[2] >= zmin && ps[2] < zmax       ) {
	  
	  int ii = (ps[0] - xmin)/dx;
	  int jj = (ps[1] - ymin)/dy;
	  int kk = (ps[2] - zmin)/dz;
	  
	  ddens[ii][jj][kk] += ms;
	}
      }

    } else if (sname.compare(its->name) == 0) {

      found_star = true;

      //
      // Position to beginning of particles
      //
      in->seekg(its->pspos);
      
      for (int j=0; j<its->nbod; j++) {
	if (its->index_size) in->read((char *)&indx, its->index_size);
	else                 indx = j;

	in->read((char *)&ms, sizeof(double));
	for (int i=0; i<3; i++) {
	  in->read((char *)&ps[i], sizeof(double));
	  xyz[i] = ps[i];
	}
	for (int i=0; i<3; i++) {
	  in->read((char *)&vs[i], sizeof(double));
	  uvw[i] = vs[i];
	}
	in->read((char *)&rtmp, sizeof(double));
	for (int i=0; i<its->niatr; i++) {
	  in->read((char *)&itmp, sizeof(int));
	}
	dattr.clear();
	for (int i=0; i<its->ndatr; i++) {
	  in->read((char *)&rtmp, sizeof(double));
	  dattr.push_back(rtmp);
	}
				// Accumulate
				// 
	if (indx > initial_star && indx <= final_star &&
	    ps[0] >= xmin && ps[0] < xmax       &&
	    ps[1] >= ymin && ps[1] < ymax       &&
	    ps[2] >= zmin && ps[2] < zmax       ) {
	  
	  int ii = (ps[0] - xmin)/dx;
	  int jj = (ps[1] - ymin)/dy;
	  int kk = (ps[2] - zmin)/dz;
	  
	  sdens[ii][jj][kk] += ms;
	}
      }

    } else if (cname.compare(its->name) == 0) {

      found_gas = true;

      //
      // Position to beginning of particles
      //
      in->seekg(its->pspos);
      
      for (int j=0; j<its->nbod; j++) {
	if (its->index_size) in->read((char *)&indx, its->index_size);
	else                 indx = j;

	in->read((char *)&ms, sizeof(double));
	for (int i=0; i<3; i++) {
	  in->read((char *)&ps[i], sizeof(double));
	  xyz[i] = ps[i];
	}
	for (int i=0; i<3; i++) {
	  in->read((char *)&vs[i], sizeof(double));
	  uvw[i] = vs[i];
	}
	in->read((char *)&rtmp, sizeof(double));
	for (int i=0; i<its->niatr; i++) {
	  in->read((char *)&itmp, sizeof(int));
	}
	dattr.clear();
	for (int i=0; i<its->ndatr; i++) {
	  in->read((char *)&rtmp, sizeof(double));
	  dattr.push_back(rtmp);
	}
      
				// Accumulate
				// 
	if (indx > initial_gas && indx <= final_gas &&
	    ps[0] >= xmin && ps[0] < xmax       &&
	    ps[1] >= ymin && ps[1] < ymax       &&
	    ps[2] >= zmin && ps[2] < zmax       ) {
	  
	  int ii = (ps[0] - xmin)/dx;
	  int jj = (ps[1] - ymin)/dy;
	  int kk = (ps[2] - zmin)/dz;
	  
	  mass[ii][jj][kk] += ms;
	  if (its->ndatr>0) { gtemp[ii][jj][kk] += ms*dattr[0]; btemp=true; }
	  if (its->ndatr>1) { gdens[ii][jj][kk] += ms*dattr[1]; bdens=true; }
	  if (its->ndatr>4) { gknud[ii][jj][kk] += ms*dattr[4]; bknud=true; }
	  if (its->ndatr>5) { gstrl[ii][jj][kk] += ms*dattr[5]; bstrl=true; }
	  for (int ll=0; ll<3; ll++) vel[ii][jj][kk][ll] += ms*vs[ll];
	  
	  // Pack gas arrays
	  //
	  part->InsertPoint(offset, &xyz[0]);
	  if (btemp) temp->InsertTuple(offset, &dattr[0]);
	  if (bdens) dens->InsertTuple(offset, &dattr[1]);
	  if (bknud) knud->InsertTuple(offset, &dattr[4]);
	  if (bstrl) strl->InsertTuple(offset, &dattr[5]);
	  velo->InsertTuple(offset, &uvw[0]);
	  offset++;
	}
      }
    }
  }

  if (!found_dark && dname.size() > 0) {
    cerr << "Could not find dark component named <" << dname << ">\n";
    exit(-1);
  }

  if (!found_dark && sname.size() > 0) {
    cerr << "Could not find star component named <" << sname << ">\n";
    exit(-1);
  }

  if (!found_gas && cname.size() > 0) {
    cerr << "Could not find gas component named <" << cname << ">\n";
    exit(-1);
  }

  for (int i=0; i<numx; i++) {
    for (int j=0; j<numy; j++) {
      for (int k=0; k<numz; k++) {
	if (found_gas && mass[i][j][k]>0.0) {
	  for (int l=0; l<3; l++) vel[i][j][k][l] /= mass[i][j][k];
	  if (btemp) gtemp[i][j][k] /= mass[i][j][k];
	  if (bdens) gdens[i][j][k] /= mass[i][j][k];
	  if (bknud) gknud[i][j][k] /= mass[i][j][k];
	  if (bstrl) gstrl[i][j][k] /= mass[i][j][k];
	  if (btemp) {
	    double vt = 0.0;
	    for (int l=0; l<3; l++) vt += vel[i][j][k][l]*vel[i][j][k][l];
	    gmach[i][j][k] = sqrt(vt*Vconv*Vconv /
				  (gamma*boltz/mm*gtemp[i][j][k]));
	  }
	  mass [i][j][k] /= dx*dy*dz;
	}
	if (found_dark) ddens[i][j][k] /= dx*dy*dz;
	if (found_star) sdens[i][j][k] /= dx*dy*dz;
      }
    }
  }
  
  vtkSmartPointer<vtkFloatArray> XX = vtkFloatArray::New();
  vtkSmartPointer<vtkFloatArray> YY = vtkFloatArray::New();
  vtkSmartPointer<vtkFloatArray> ZZ = vtkFloatArray::New();
  float f;

  for (int i=0; i<numx; i++)  XX->InsertTuple(i, &(f=xmin + dx*(0.5 + i)));
  for (int j=0; j<numy; j++)  YY->InsertTuple(j, &(f=ymin + dy*(0.5 + j)));
  for (int k=0; k<numz; k++)  ZZ->InsertTuple(k, &(f=zmin + dz*(0.5 + k)));
    
  vtkSmartPointer<vtkRectilinearGrid> dataSet = vtkRectilinearGrid::New();
  dataSet->SetDimensions(numx, numy, numz);
  dataSet->SetXCoordinates (XX);
  dataSet->SetYCoordinates (YY);
  dataSet->SetZCoordinates (ZZ);

  vtkSmartPointer<vtkDataArray> Temp;
  vtkSmartPointer<vtkDataArray> Dens;
  vtkSmartPointer<vtkDataArray> Knud;
  vtkSmartPointer<vtkDataArray> Strl;
  vtkSmartPointer<vtkDataArray> Mach;
  vtkSmartPointer<vtkDataArray> density;
  vtkSmartPointer<vtkDataArray> velocity;
  vtkSmartPointer<vtkPoints>    points;
  vtkSmartPointer<vtkDataArray> dRho;
  vtkSmartPointer<vtkDataArray> sRho;


  if (found_gas) {
    Temp       = vtkFloatArray::New();
    Dens       = vtkFloatArray::New();
    Knud       = vtkFloatArray::New();
    Strl       = vtkFloatArray::New();
    Mach       = vtkFloatArray::New();
    density    = vtkFloatArray::New();
    velocity   = vtkFloatArray::New();
    points     = vtkPoints::New();

    Temp      -> SetName("gas_temp");
    Dens      -> SetName("gas_dens");
    Knud      -> SetName("Knudsen");
    Strl      -> SetName("Strouhal");
    Mach      -> SetName("Mach");
    density   -> SetName("density");
    velocity  -> SetName("velocity");
    velocity  -> SetNumberOfComponents(3);
  }

  if (found_dark) {
    dRho = vtkFloatArray::New();
    dRho->SetName("Dark density");
  }

  if (found_star) {
    sRho = vtkFloatArray::New();
    sRho->SetName("Star density");
  }

  vtkSmartPointer<vtkUnsignedCharArray> visible = vtkUnsignedCharArray::New();

  unsigned activ = 0;
  unsigned blank = 0;

  for (int k=0; k<numz; k++) {

    for (int j=0; j<numy; j++) {

      for (int i=0; i<numx; i++) {

	double x0 = min<double>(
				max<double>(pos[i][j][k][0], xmin+0.501*dx),
				xmax - 0.501*dx
				);

	double y0 = min<double>(
				max<double>(pos[i][j][k][1], ymin+0.501*dy),
				ymax - 0.501*dy
				);

	double z0 = min<double>(
				max<double>(pos[i][j][k][2], zmin+0.501*dz),
				zmax - 0.501*dz
				);

	vtkIdType n = dataSet->FindPoint(x0, y0, z0);

	if (found_gas) {

	  density->InsertTuple(n, &mass[i][j][k]);
	
	  velocity->InsertTuple(n, &vel[i][j][k][0]);
	
	  if (btemp) Temp->InsertTuple(n, &gtemp[i][j][k]);
	  if (bdens) Dens->InsertTuple(n, &gdens[i][j][k]);
	  if (bknud) Knud->InsertTuple(n, &gknud[i][j][k]);
	  if (bstrl) Strl->InsertTuple(n, &gstrl[i][j][k]);
	  if (btemp) Mach->InsertTuple(n, &gmach[i][j][k]);

	  if (mass[i][j][k]>0.0 || !mask) {
	    visible->InsertValue(n, 1);
	    activ++;
	  }
	  else {
	    visible->InsertValue(n, 0);
	    blank++;
	  }
	}

	if (found_dark)
	  dRho->InsertTuple(n, &ddens[i][j][k]);

	if (found_star)
	  sRho->InsertTuple(n, &sdens[i][j][k]);


	offset++;
      }
    }
  }

  // Add all the fields

  if (found_gas) {
    if (btemp) dataSet->GetPointData()->AddArray(Temp);
    if (bdens) dataSet->GetPointData()->AddArray(Dens);
    if (bknud) dataSet->GetPointData()->AddArray(Knud);
    if (bstrl) dataSet->GetPointData()->AddArray(Strl);
    if (btemp) dataSet->GetPointData()->AddArray(Mach);
    dataSet->GetPointData()->AddArray(density);
    dataSet->GetPointData()->SetVectors(velocity);
    // dataSet->SetPointVisibilityArray(visible);
  }
    
  if (found_dark)
    dataSet->GetPointData()->AddArray(dRho);
  
  if (found_star)
    dataSet->GetPointData()->AddArray(sRho);
  
  dataSet->SetDimensions(numx, numy, numz);

  // Print out the VTK file (in XML)

  string gname = outfile + ".vtr";
  vtkSmartPointer<vtkXMLRectilinearGridWriter> writer = 
    vtkXMLRectilinearGridWriter::New();
     writer->SetInput(dataSet);
     writer->SetFileName(gname.c_str());
     writer->Write();

  if (mask)
    cout << blank << " blank voxels and " << activ << " active ones" << endl;

  return (0);
}


