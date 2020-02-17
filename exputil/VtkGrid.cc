#include <VtkGrid.H>

VtkGrid::VtkGrid(int nx, int ny, int nz,
		 double xmin, double xmax,
		 double ymin, double ymax,
		 double zmin, double zmax) :
  nx(nx), ny(ny), nz(nz),
  xmin(xmin), xmax(xmax),
  ymin(ymin), ymax(ymax),
  zmin(zmin), zmax(zmax)
{
  // Set knots
  auto XX   = vtkFloatArray::New();
  auto YY   = vtkFloatArray::New();
  auto ZZ   = vtkFloatArray::New();

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

  // Create a pointer to a VTK Unstructured Grid data set
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

    for (int j=0; j<ny; j++) {
      float y = ymin + (ymax - ymin)*j/(ny-1);

      for (int i=0; i<nx; i++) {
	float x = xmin + (xmax - xmin)*i/(nx-1);

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

    for (int k=0; k<nz; k++) {
      float z = zmin + (zmax - zmin)*k/(nz-1);

      for (int j=0; j<ny; j++) {
	float y = ymin + (ymax - ymin)*j/(ny-1);

	for (int i=0; i<nx; i++) {
	  float x = xmin + (xmax - xmin)*i/(nx-1);
	  vtkIdType n = dataSet->FindPoint(x, y, z);
	  
	  if (n>=0) {
	    float f = static_cast<float>(data[(k*ny+j)*nx + i]);
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
