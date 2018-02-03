#include <VtkPCA.H>

VtkPCA::VtkPCA(int N, bool smooth) : nmax(N), smooth(smooth)
{
  // Set knots
  auto XX   = vtkFloatArray::New();
  auto YY   = vtkFloatArray::New();
  auto ZZ   = vtkFloatArray::New();

  XX -> SetName("Coef index");
  YY -> SetName("PCA index");

  float f;
  int k = 0;

  for (int i=1; i<=nmax; i++, k++) {
    XX->InsertTuple(k, &(f=i));
    YY->InsertTuple(k, &(f=i));
  }

  ZZ->InsertTuple(0, &(f=0));

  // Create a pointer to a VTK Unstructured Grid data set
  dataSet = vtkRectilinearGridP::New();

  dataSet->SetDimensions(nmax, nmax, 1);
  dataSet->SetXCoordinates (XX);
  dataSet->SetYCoordinates (YY);
  dataSet->SetZCoordinates (ZZ);
}

void VtkPCA::Add(const Vector& eval, const Matrix& evec, int m)
{
  vtkFloatArrayP V = vtkFloatArrayP::New();
  vtkFloatArrayP T = vtkFloatArrayP::New();

  // Insert grid data
  //
  float f;			// Temp float storage
  for (int i=0; i<nmax; i++) {
    for (int j=0; j<nmax; j++) {
      float x = i+1, y = j+1;
      vtkIdType n = dataSet->FindPoint(x, y, 0);

      if (n>=0) {
	f = evec[i+1][j+1];
	if (smooth) f *= eval[i+1];
	if (std::isnan(f)) f = 0.0;
	T->InsertTuple(n, &f);
      } else {
	std::cout << "Could not find point at (" << x << ", " << y << ")"
		  << std::endl;
      }
    }
  }

  // Add eigenvalues
  for (int i=0; i<nmax; i++) {
    f = eval[i+1];
    if (std::isnan(f)) f = 0.0;
    V->InsertTuple(i, &f);
  }

  // Add arrays
  vecs.push_back(T);
  vals.push_back(V);

  // Add label
  std::ostringstream lab;
  lab << m;
  elab.push_back(lab.str());
}

void VtkPCA::Add(const Vector& eval, const Matrix& evec, int l, int m, char tag)
{
  vtkFloatArrayP V = vtkFloatArrayP::New();
  vtkFloatArrayP T = vtkFloatArrayP::New();

  // Insert grid data
  //
  float f;			// Temp float storage
  for (int i=0; i<nmax; i++) {
    for (int j=0; j<nmax; j++) {
      float x = i+1, y = j+1;
      vtkIdType n = dataSet->FindPoint(x, y, 0);

      if (n>=0) {
	f = evec[i+1][j+1];
	if (smooth) f *= eval[i+1];
	if (std::isnan(f)) f = 0.0;
	T->InsertTuple(n, &f);
      } else {
	std::cout << "Could not find point at (" << x << ", " << y << ")"
		  << std::endl;
      }
    }
  }

  // Add eigenvalues
  //
  for (int i=0; i<nmax; i++) {
    f = eval[i+1];
    if (std::isnan(f)) f = 0.0;
    V->InsertTuple(i, &f);
  }

  // Add arrays
  vecs.push_back(T);
  vals.push_back(V);

  // Add label
  std::ostringstream lab;
  lab << l << "_" << m << "_" << tag;
  elab.push_back(lab.str());
}

void VtkPCA::Write(const std::string& name)
{
  // Create a writer
  auto writer = vtkRectilinearGridWriterP::New();

  // Create the filename with the correct extension for Paraview
  std::ostringstream filename;
  filename << name << "." << writer->GetDefaultFileExtension();

  writer->SetFileName(filename.str().c_str());

  // Add everything to the data set
  for (size_t k=0; k<elab.size(); k++) {
    std::string lab1 = "Value " + elab[k];
    vals[k] -> SetName(lab1.c_str());
    std::string lab2 = "Vector " + elab[k];
    vecs[k] -> SetName(lab2.c_str());

    // Add fields
    dataSet->GetFieldData()->AddArray(vals[k]);
    dataSet->GetPointData()->AddArray(vecs[k]);
  }

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
