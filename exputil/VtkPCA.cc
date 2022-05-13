#include <VtkPCA.H>

#include <map>

VtkPCA::VtkPCA(int N, bool reorder, bool smooth) :
  nmax(N), reorder(reorder), smooth(smooth)
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

void VtkPCA::Add(const Eigen::VectorXd& Coef, 
		 const Eigen::VectorXd& Hall,
		 const Eigen::VectorXd& SnrV,
		 const Eigen::VectorXd& Eval,
		 const Eigen::MatrixXd& Evec,
		 const Eigen::MatrixXd& Covr,
		 int m)
{
  vtkFloatArrayP C = vtkFloatArrayP::New();
  vtkFloatArrayP H = vtkFloatArrayP::New();
  vtkFloatArrayP S = vtkFloatArrayP::New();
  vtkFloatArrayP V = vtkFloatArrayP::New();
  vtkFloatArrayP T = vtkFloatArrayP::New();
  vtkFloatArrayP Y = vtkFloatArrayP::New();
  vtkFloatArrayP U = vtkFloatArrayP::New();

  // Make reorder map
  //
  std::vector<int> R;

  if (reorder) {
    std::multimap<double, int> reord;
    for (int i=0; i<=SnrV.size(); i++) reord.insert(std::make_pair(SnrV[i], i));;
    for (auto i=reord.rbegin(); i!=reord.rend(); i++) R.push_back(i->second);
  } else {
    for (int i=0; i<=SnrV.size(); i++) R.push_back(i);
  }
    
  // Insert grid data
  //
  float f;			// Temp float storage
  for (int i=0; i<nmax; i++) {
    for (int j=0; j<nmax; j++) {
      float x = i+1, y = j+1;
      vtkIdType n = dataSet->FindPoint(x, y, 0);

      if (n>=0) {
	f = Evec(R[i], j);	// Evec
	if (smooth) f *= Hall[R[i]];
	if (std::isnan(f)) f = 0.0;
	T->InsertTuple(n, &f);
	f *= f;
	Y->InsertTuple(n, &f);  // Evec squared
	f = Covr(i, j);		// Covariance
	U->InsertTuple(n, &f);
      } else {
	std::cout << "Could not find point at (" << x << ", " << y << ")"
		  << std::endl;
      }
    }
  }

  // Add coefficients (should not be reordered)
  for (int i=0; i<nmax; i++) {
    f = Coef[i+1];
    if (std::isnan(f)) f = 0.0;
    C->InsertTuple(i, &f);
  }

  // Add Hall smoothing
  for (int i=0; i<nmax; i++) {
    f = Hall[R[i]];
    if (std::isnan(f)) f = 0.0;
    H->InsertTuple(i, &f);
  }

  // Add S/N
  for (int i=0; i<nmax; i++) {
    f = SnrV[R[i]];
    if (std::isnan(f)) f = 0.0;
    S->InsertTuple(i, &f);
  }

  // Add eigenvalues
  for (int i=0; i<nmax; i++) {
    f = Eval[R[i]];
    if (std::isnan(f)) f = 0.0;
    V->InsertTuple(i, &f);
  }

  // Add arrays
  coef.push_back(C);
  hall.push_back(H);
  snrv.push_back(S);
  eval.push_back(V);
  vecs.push_back(T);
  vec2.push_back(Y);
  covr.push_back(U);

  // Add label
  std::ostringstream lab;
  lab << m;
  elab.push_back(lab.str());
}

void VtkPCA::Add(const Eigen::VectorXd& Coef,
		 const Eigen::VectorXd& Hall,
		 const Eigen::VectorXd& SnrV,
		 const Eigen::VectorXd& Eval,
		 const Eigen::MatrixXd& Evec,
		 const Eigen::MatrixXd& Covr,
		 int l, int m)
{
  vtkFloatArrayP C = vtkFloatArrayP::New();
  vtkFloatArrayP H = vtkFloatArrayP::New();
  vtkFloatArrayP S = vtkFloatArrayP::New();
  vtkFloatArrayP V = vtkFloatArrayP::New();
  vtkFloatArrayP T = vtkFloatArrayP::New();
  vtkFloatArrayP Y = vtkFloatArrayP::New();
  vtkFloatArrayP U = vtkFloatArrayP::New();

  // Make reorder map
  //
  std::vector<int> R;

  if (reorder) {
    std::multimap<double, int> reord;
    for (int i=0; i<=SnrV.size(); i++) reord.insert(std::make_pair(SnrV[i], i));
    for (auto i=reord.rbegin(); i!=reord.rend(); i++) R.push_back(i->second);
  } else {
    for (int i=0; i<=SnrV.size(); i++) R.push_back(i);
  }

  // Insert grid data
  //
  float f;			// Temp float storage
  for (int i=0; i<nmax; i++) {
    for (int j=0; j<nmax; j++) {
      float x = i+1, y = j+1;
      vtkIdType n = dataSet->FindPoint(x, y, 0);

      if (n>=0) {
	f = Evec(R[i], j);
	if (smooth) f *= Hall[R[i]];
	if (std::isnan(f)) f = 0.0;
	T->InsertTuple(n, &f);
	f *= f;
	Y->InsertTuple(n, &f);
	
	f = Covr(i, j);		// Covariance
	U->InsertTuple(n, &f);
      } else {
	std::cout << "Could not find point at (" << x << ", " << y << ")"
		  << std::endl;
      }
    }
  }

  // Add coefficients (should not be reordered)
  //
  for (int i=0; i<nmax; i++) {
    f = Coef[i+1];
    if (std::isnan(f)) f = 0.0;
    C->InsertTuple(i, &f);
  }

  // Add Hall smoothing
  //
  for (int i=0; i<nmax; i++) {
    f = Hall[R[i]];
    if (std::isnan(f)) f = 0.0;
    H->InsertTuple(i, &f);
  }

  // Add S/N
  //
  for (int i=0; i<nmax; i++) {
    f = SnrV[R[i]];
    if (std::isnan(f)) f = 0.0;
    S->InsertTuple(i, &f);
  }

  // Add eigenvalues
  //
  for (int i=0; i<nmax; i++) {
    f = Eval[R[i]];
    if (std::isnan(f)) f = 0.0;
    V->InsertTuple(i, &f);
  }

  // Add arrays
  coef.push_back(C);
  hall.push_back(H);
  snrv.push_back(S);
  eval.push_back(V);
  vecs.push_back(T);
  vec2.push_back(Y);
  covr.push_back(U);

  // Add label
  std::ostringstream lab;
  lab << l << "_" << m;
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
    std::string lab0 = "Coefs " + elab[k];
    coef[k] -> SetName(lab0.c_str());
    std::string lab1 = "HallC " + elab[k];
    hall[k] -> SetName(lab1.c_str());
    std::string lab2 = "S/N " + elab[k];
    snrv[k] -> SetName(lab2.c_str());
    std::string lab3 = "Evals " + elab[k];
    eval[k] -> SetName(lab3.c_str());
    std::string lab4 = "Evecs " + elab[k];
    vecs[k] -> SetName(lab4.c_str());
    std::string lab5 = "Covar " + elab[k];
    covr[k] -> SetName(lab5.c_str());
    std::string lab6 = "Evec2 " + elab[k];
    vec2[k] -> SetName(lab6.c_str());

    // Add fields
    dataSet->GetFieldData()->AddArray(coef[k]);
    dataSet->GetFieldData()->AddArray(hall[k]);
    dataSet->GetFieldData()->AddArray(snrv[k]);
    dataSet->GetFieldData()->AddArray(eval[k]);
    dataSet->GetPointData()->AddArray(vecs[k]);
    dataSet->GetPointData()->AddArray(covr[k]);
    dataSet->GetPointData()->AddArray(vec2[k]);
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
