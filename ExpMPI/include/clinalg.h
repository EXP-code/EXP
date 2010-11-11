
int lu_decomp(CMatrix& a, int* indx, double& d);
int linear_solve(CMatrix& a, CVector& b, CVector& x);
void lu_backsub(CMatrix& a, int* indx, CVector& b);
int inverse(CMatrix& a, CMatrix& b);
KComplex determinant(CMatrix& a);
KComplex lu_determinant(CMatrix& a, double& d);
CMatrix sub_matrix(CMatrix& in, 
		   int ibeg, int iend, int jbeg, int jend, 
		   int ioff=0, int joff=0);

void embed_matrix(CMatrix& to, CMatrix& from, int rbeg, int cbeg);
void embed_matrix(Matrix& to, Matrix& from, int rbeg, int cbeg);
void embed_matrix(CMatrix& to, Matrix& from, int rbeg, int cbeg);
