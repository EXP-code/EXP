
int lu_decomp(Matrix& a, int* indx, double& d);
int linear_solve(Matrix& a, Vector& b, Vector& x);
void lu_backsub(Matrix& a, int* indx, Vector& b);
int inverse(Matrix& a, Matrix& b);
void improve(Matrix& a, Matrix& lud, int* indx, Vector& b, Vector& x);
double determinant(Matrix& a);
double lu_determinant(Matrix& a, double& d);
Matrix sub_matrix(Matrix& in, 
		   int ibeg, int iend, int jbeg, int jend, 
		   int ioff=0, int joff=0);

void embed_matrix(Matrix& to, Matrix& from, int rbeg, int cbeg);
void embed_matrix(Matrix& to, Matrix& from, int rbeg, int cbeg);
void embed_matrix(Matrix& to, Matrix& from, int rbeg, int cbeg);
