#include <VtkPCA.H>

VtkPCA::VtkPCA(int N) : nmax(N)
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

void VtkPCA::AddMatrix(const Vector& eval, const Matrix& evec, int indx)
{
  vtkFloatArrayP V = vtkFloatArrayP::New();
  vtkFloatArrayP T = vtkFloatArrayP::New();

  // Insert grid data
  //
  float f;			// Temp float storage
  for (int i=0; i<nmax; i++) {
    for (int j=0; j<nmax; j++) {
      float x = i+1, y = i+1;
      vtkIdType n = dataSet->FindPoint(x, y, 0);

      if (n>=0) {
	T->InsertTuple(n, &(f=evec[i+1][j+1]));
      } else {
	std::cout << "Could not find point at (" << x << ", " << y << ")"
		  << std::endl;
      }
    }
  }

  // Add eigenvalues
  for (int i=0; i<nmax; i++)
    V->InsertTuple(i, &(f=eval[i+1]));

  // Add arrays
  vals.push_back(T);
  vecs.push_back(V);

  // Add label
  std::ostringstream lab;
  lab << "i_" << indx;
  elab.push_back(lab.str());
}

void VtkPCA::Write(const std::string& name)
{
  // Create a writer
  auto writer = vtkRectilinearGridWriterP::New();

  writer->SetFileName(name.c_str());

  // Add everything to the data set
  for (size_t k=0; k<elab.size(); k++) {
    std::string lab1 = "Eigenvalues_" + elab[k];
    vals[k] -> SetName(lab1.c_str());
    std::string lab2 = "Eigenvectors_" + elab[k];
    vecs[k] -> SetName(lab2.c_str());

    // Add fields
    dataSet->GetPointData()->AddArray(vals[k]);
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
  writer->Write();
}
