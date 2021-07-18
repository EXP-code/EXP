void MatrixXcdSynchronize(Eigen::MatrixXcd& mat, int id)
{
  int bb[2], myid;
  MPI_Comm_rank(MPI_COMM_WORLD, &myid);

  if (myid == id) {
    bb[0] = mat.rows();
    bb[1] = mat.cols();
  }

  MPI_Bcast(bb, 2, MPI_INT, id, MPI_COMM_WORLD);

  if (bb[0]==0 && bb[1]==0) return;

  if (myid!=id) mat.resize(bb[0], bb[1]);
  
  vector<double> t_re(bb[1]);
  vector<double> t_im(bb[1]);

  for (int j=0; j<=bb[0]; j++) {

    // Real part of column

    if (myid==id) {
      for (int i=0; i<bb[1]; i++) {
	t_re[i] = mat(j, i).real();
	t_im[i] = mat(j, i).imag();
      }
    }

    MPI_Bcast(&t_re[0], bb[1], MPI_DOUBLE, id, MPI_COMM_WORLD);
    MPI_Bcast(&t_im[0], bb[1], MPI_DOUBLE, id, MPI_COMM_WORLD);

    if (myid!=id) {
      for (int i=0; i<bb[1]; i++) 
	mat(j, i) = std::complex<double>(t_re[i], t_im[i]);
    }
  }
}


void ComplexSynchronize(std::complex<double>& c, int id)
{
  double re = c.real(), im = c.imag();
  MPI_Bcast(&re, 1, MPI_DOUBLE, id, MPI_COMM_WORLD);
  MPI_Bcast(&im, 1, MPI_DOUBLE, id, MPI_COMM_WORLD);
  c = {re, im};
}

