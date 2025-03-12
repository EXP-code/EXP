#include <iomanip>
#include <VtkGrid.H>

#ifdef HAVE_VTK

VtkGrid::VtkGrid(int nx, int ny, int nz,
		 double xmin, double xmax,
		 double ymin, double ymax,
		 double zmin, double zmax) :
  ThreeDGrid(nx, ny, nz, xmin, xmax, ymin, ymax, zmin, zmax)
{
  // Set knots
  auto XX   = vtkSmartPointer<vtkFloatArray>::New();
  auto YY   = vtkSmartPointer<vtkFloatArray>::New();
  auto ZZ   = vtkSmartPointer<vtkFloatArray>::New();

  XX -> SetName("X");
  YY -> SetName("Y");
  if (nz>1) ZZ -> SetName("Z");

  float f;
  int k = 0;

  for (int i=1; i<=nx; i++, k++) {
    float f = xmin + (xmax - xmin)*(i-1)/(nx-1);
    XX->InsertTuple(k, &f);
  }

  k = 0;
  for (int i=1; i<=ny; i++, k++) {
    float f = ymin + (ymax - ymin)*(i-1)/(ny-1);
    YY->InsertTuple(k, &f);
  }

  if (nz<2) {
    ZZ->InsertTuple(0, &(f=0));
  } else {
    k = 0;
    for (int i=1; i<=nz; i++, k++) {
      float f = zmin + (zmax - zmin)*(i-1)/(nz-1);
      ZZ->InsertTuple(k, &f);
    }
  }    

  // Create a pointer to a VTK Rectilinear Grid data set
  dataSet = vtkRectilinearGridP::New();

  dataSet->SetDimensions(nx, ny, nz);
  dataSet->SetXCoordinates (XX);
  dataSet->SetYCoordinates (YY);
  dataSet->SetZCoordinates (ZZ);
}

void VtkGrid::Add(const std::vector<double>& data, const std::string& name)
{
  vtkFloatArrayP T = vtkFloatArrayP::New();

  // Insert grid data
  //

  if (nz<2) {

    for (int i=0; i<nx; i++) {
      double x = xmin + (xmax - xmin)*i/(nx-1);
      for (int j=0; j<ny; j++) {
	double y = ymin + (ymax - ymin)*j/(ny-1);
	vtkIdType n = dataSet->FindPoint(x, y, 0);

	if (n>=0) {
	  float f = static_cast<float>(data[j*nx + i]);
	  T->InsertTuple(n, &f);
	} else {
	  std::cout << "Could not find point at (" << x << ", " << y << ")"
		    << std::endl;
	}
      }
    }

  } else {

    for (int i=0; i<nx; i++) {
      double x = xmin + (xmax - xmin)*i/(nx-1);
      for (int j=0; j<ny; j++) {
	double y = ymin + (ymax - ymin)*j/(ny-1);
	for (int k=0; k<nz; k++) {
	  double z = zmin + (zmax - zmin)*k/(nz-1);

	  vtkIdType n = dataSet->FindPoint(x, y, z);

	  if (n>=0) {
	    float f = static_cast<float>(data[(k*ny + j)*nx + i]);
	    T->InsertTuple(n, &f);
	  } else {
	    std::cout << "Could not find point at (" << x << ", " << y << ", "<< z << ")"
		      << std::endl;
	  }
	}
      }
    }
    
  }

  T -> SetName(name.c_str());
  dataSet->GetPointData()->AddArray(T);
}


void VtkGrid::Write(const std::string& name)
{
  // Create a writer
  auto writer = vtkRectilinearGridWriterP::New();

  // Create the filename with the correct extension for Paraview
  std::ostringstream filename;
  filename << name << "." << writer->GetDefaultFileExtension();

  writer->SetFileName(filename.str().c_str());

  // Remove unused memory
  dataSet->Squeeze();

  // Write the data
#if VTK_MAJOR_VERSION <= 5
  writer->SetInput(dataSet);
#else
  writer->SetInputData(dataSet);
#endif
  writer->SetDataModeToAscii();
  // writer->SetDataModeToBinary();
  writer->Write();
}

#else

VtkGrid::VtkGrid(int nx, int ny, int nz,
		 double xmin, double xmax,
		 double ymin, double ymax,
		 double zmin, double zmax) :
  ThreeDGrid(nx, ny, nz, xmin, xmax, ymin, ymax, zmin, zmax)
{
  // Set knots
  coord["X"].resize(nx);
  coord["Y"].resize(ny);
  coord["Z"].resize(nz);

  for (int i=0; i<nx; i++) coord["X"][i] = xmin + (xmax - xmin)*i/(nx-1);
  for (int i=0; i<ny; i++) coord["Y"][i] = ymin + (ymax - ymin)*i/(ny-1);

  if (nz>1) {
    for (int i=0; i<nz; i++) coord["Z"][i] = zmin + (zmax - zmin)*i/(nz-1);
  } else {
    coord["Z"][0] = 0.0;
  } 
}

void VtkGrid::replaceAll(std::string& str,
			 const std::string& from,
			 const std::string& to)
{
  // Sanity check: nothing to replace
  if (from.empty()) return;

  // Look for strings to replace
  size_t start_pos = 0;
  while((start_pos = str.find(from, start_pos)) != std::string::npos) {
    str.replace(start_pos, from.length(), to);
    start_pos += to.length();	// In case 'to' contains 'from'
  }
}


void VtkGrid::Add(const std::vector<double>& data, const std::string& name)
{
  // Remove XML sensitive characters; paraview bug?
  std::string newName(name);
  replaceAll(newName, ">", ".gt.");
  replaceAll(newName, "<", ".lt.");

  dataSet[newName].resize(nx*ny*nz);

  // Insert grid data
  //
  int I = 0;

  for (int k=0; k<nz; k++) {
    for (int j=0; j<ny; j++) {
      for (int i=0; i<nx; i++) {
	dataSet[newName][I++] = static_cast<float>(data[(k*ny + j)*nx + i]);
      }
    }
  }

}


void VtkGrid::writeBeg(std::ofstream & fout)
{
  int zero = 0;
 fout << "<?xml version=\"1.0\"?>" << std::endl;
  fout << "<VTKFile type=\"RectilinearGrid\" version=\"0.1\" byte_order=\"LittleEndian\" header_type=\"UInt32\">" << std::endl;
  fout << "  <RectilinearGrid WholeExtent=\""
       << zero << " " << nx-1 << " "
       << zero << " " << ny-1 << " "
       << zero << " " << nz-1
       << "\">" << std::endl;
  fout << "  <Piece Extent=\""
       << zero << " " << nx-1 << " "
       << zero << " " << ny-1 << " "
       << zero << " " << nz-1
       << "\">" << std::endl;
}

void VtkGrid::writeFields(std::ofstream & fout)
{
 fout << "    <PointData>" << std::endl;
  for (auto v : dataSet) {
    // Get ranges
    auto vmin = std::min_element(v.second.begin(), v.second.end());
    auto vmax = std::max_element(v.second.begin(), v.second.end());
    fout << "      <DataArray type=\"Float32\" Name=\"" << v.first << "\" "
	 << " format=\"ascii\" RangeMin=\"" << *vmin << "\" RangeMax=\"" << *vmax
	 << "\">" << std::endl;
    int cnt = 0;
    for (auto & d : v.second) {
      if (cnt % 6 == 0) fout << "         ";
      fout << std::scientific << std::setprecision(8) << d << " ";
      if (++cnt % 6 == 0) fout << std::endl;
    }
    if (cnt % 6) fout << std::endl;
    fout << "      </DataArray>" << std::endl;
  }
  fout << "    </PointData>" << std::endl;

  // No cell data here so this stanza is empty, but it's part of the spec
  fout << "    <CellData>" << std::endl;
  fout << "    </CellData>" << std::endl;
}

void VtkGrid::writeCoords(std::ofstream & fout)
{
  fout << "    <Coordinates>" << std::endl;
  for (auto & c : coord) {
    fout << "      <DataArray type=\"Float32\" Name=\""
	 << c.first << "\" format=\"ascii\" RangeMin=\""
	 << *std::min_element(c.second.begin(), c.second.end())
       << "\" RangeMax=\""
	 << *std::max_element(c.second.begin(), c.second.end()) << "\">"
	 << std::endl;

    int cnt = 0;
    for (auto & v : c.second) {
      if (cnt % 6 == 0) fout << "        ";
      fout << std::scientific << std::setprecision(8) << v << " ";
      if (++cnt % 6 == 0) fout << std::endl;
    }
    if (cnt % 6) fout << std::endl;
    fout << "      </DataArray>" << std::endl;
  }
  fout << "    </Coordinates>" << std::endl;
}

void VtkGrid::writeEnd(std::ofstream & fout)
{
  fout << "  </Piece>" << std::endl;
  fout << "  </RectilinearGrid>" << std::endl;
  fout << "</VTKFile>" << std::endl;
}

void VtkGrid::Write(const std::string& name)
{
  // Create the filename with the correct extension for Paraview
  std::ostringstream filename;
  filename << name << ".vtr";

  std::ofstream fout(filename.str());
  if (fout) {

  } else {
    throw std::runtime_error
      (std::string("VtkGrid::Write: could not open file <") + filename.str() + ">");
  }

  writeBeg(fout);		// VTK XML preamble; open stanzas
  writeFields(fout);		// Write each data field in the dataset map
  writeCoords(fout); 		// Write the coordinate grids
  writeEnd(fout);		// VTK XML finish; close open stanzas
}

#endif
