#include <cmath>
#include <numerical.h>
#include <Vector.h>

#include <Eigen/Eigen>

extern "C" {
  void dsyevd_(char* jobz, char* uplo, int* n, double* a, int* lda, double *w,
	       double* work, int* lwork, int *iwork, int* liwork, int* info);
}

/*

  call dsyevd(jobz, uplo, n, a, lda, w, work, lwork, iwork, liwork, info)

  Description

  The routine computes all the eigenvalues, and optionally all the
  eigenvectors, of a real symmetric matrix A. In other words, it can
  compute the spectral factorization of A as: A = Z*Λ*ZT.

  Here Λ is a diagonal matrix whose diagonal elements are the
  eigenvalues λi, and Z is the orthogonal matrix whose columns are the
  eigenvectors zi. Thus,

  A*zi = λi*zi for i = 1, 2, ..., n.

  If the eigenvectors are requested, then this routine uses a divide
  and conquer algorithm to compute eigenvalues and
  eigenvectors. However, if only eigenvalues are required, then it
  uses the Pal-Walker-Kahan variant of the QL or QR algorithm.

  Note that for most cases of real symmetric eigenvalue problems the
  default choice should be ?syevr function as its underlying algorithm
  is faster and uses less workspace. ?syevd requires more workspace
  but is faster in some cases, especially for large matrices.

  Input Parameters

  jobz

    CHARACTER*1. Must be 'N' or 'V'.

    If jobz = 'N', then only eigenvalues are computed.

    If jobz = 'V', then eigenvalues and eigenvectors are computed.

  uplo

    CHARACTER*1. Must be 'U' or 'L'.

    If uplo = 'U', a stores the upper triangular part of A.

    If uplo = 'L', a stores the lower triangular part of A.

  n

    INTEGER. The order of the matrix A (n ≥ 0).
  a

    REAL for ssyevd

    DOUBLE PRECISION for dsyevd

    Array, DIMENSION (lda, *).

    a(lda,*) is an array containing either upper or lower triangular part of the symmetric matrix A, as specified by uplo.

    The second dimension of a must be at least max(1, n).

  lda

    INTEGER. The first dimension of the array a.

    Must be at least max(1, n).

  work

    REAL for ssyevd

    DOUBLE PRECISION for dsyevd.

    Workspace array, DIMENSION at least lwork.

  lwork

    INTEGER.

    The dimension of the array work.

    Constraints:

    if n ≤ 1, then lwork ≥ 1;

    if jobz = 'N' and n > 1, then lwork ≥ 2*n + 1;

    if jobz = 'V' and n > 1, then lwork ≥ 2*n^2+ 6*n + 1.

    If lwork = -1, then a workspace query is assumed; the routine only
    calculates the required sizes of the work and iwork arrays,
    returns these values as the first entries of the work and iwork
    arrays, and no error message related to lwork or liwork is issued
    by xerbla. See Application Notes for details.

  iwork

    INTEGER.

    Workspace array, its dimension max(1, liwork).

  liwork

    INTEGER.

    The dimension of the array iwork.

    Constraints:

    if n ≤ 1, then liwork ≥ 1;

    if jobz = 'N' and n > 1, then liwork ≥ 1;

    if jobz = 'V' and n > 1, then liwork ≥ 5*n + 3.

    If liwork = -1, then a workspace query is assumed; the routine
    only calculates the required sizes of the work and iwork arrays,
    returns these values as the first entries of the work and iwork
    arrays, and no error message related to lwork or liwork is issued
    by xerbla. See Application Notes for details.

Output Parameters

  w

    REAL for ssyevd

    DOUBLE PRECISION for dsyevd

    Array, DIMENSION at least max(1, n).

    If info = 0, contains the eigenvalues of the matrix A in ascending
    order. See also info.  

  a
    If jobz = 'V', then on exit this array is overwritten by the
    orthogonal matrix Z which contains the eigenvectors of A.  

    work(1)

    On exit, if lwork > 0, then work(1) returns the required minimal
    size of lwork.  

    iwork(1)

    On exit, if liwork > 0, then iwork(1) returns the required minimal
    size of liwork.

  info

    INTEGER.

    If info = 0, the execution is successful.

    If info = i, then the algorithm failed to converge; i indicates
    the number of off-diagonal elements of an intermediate tridiagonal
    form which did not converge to zero.

    If info = i, and jobz = 'V', then the algorithm failed to compute
    an eigenvalue while working on the submatrix lying in rows and
    columns info/(n+1) through mod(info,n+1).

    If info = -i, the i-th parameter had an illegal value.
*/


Vector Symmetric_Eigenvalues_SYEVD(Matrix& a, Matrix& ef, int M)
{
  int lo = a.getrlow();
  int hi = a.getrhigh();
  int N = hi - lo + 1;

  double *A = new double [N*N];
  for (int i=lo; i<=hi; i++)
    for (int j=i; j<=hi; j++) A[N*(j-lo) + i-lo] = a[i][j];

  double *W = new double [N];

  int lwork = 2*N*N + 6*N + 1;
  double *work = new double [lwork];
  int liwork = 5*N+3;
  int *iwork = new int [liwork];
  int info;

  char jobz = 'V';
  char uplo = 'U';

  dsyevd_(&jobz, &uplo, &N, A, &N, W, work, &lwork, iwork, &liwork, &info);

  M = min<double>(N, M);

  Vector ev(1, M);
  ef.setsize(1, M, lo, hi);

  for (int i=0; i<M; i++) {
    ev[i+1] = W[N-i-1];
    for (int j=0; j<N; j++) ef[i+1][j+lo] = A[N*(N-i-1) + j];
  }

  if (info != 0) {
    std::cout << "Symmetric_Eigenvalues_SYEVD: failed to converge with i="
	      << info << std::endl;
  }

  delete [] A;
  delete [] W;
  delete [] work;
  delete [] iwork;

  return ev;
}

Vector Symmetric_Eigenvalues_SVD(Matrix& a, Matrix& ef, int M, bool Large)
{
  int lo = a.getrlow();
  int hi = a.getrhigh();
  int n = hi - lo + 1;

  // Eigen comparison
  Eigen::MatrixXd mm(n, n);
  for (int i=0; i<n; i++) {
    for (int j=i; j<n; j++) {
      mm(i, j) = a[i+lo][j+lo];
      if (i!=j) mm(j, i) = mm(i, j);
    }
  }
  
  M = std::min<double>(n, M);

  Vector ev(1, M);
  ef.setsize(1, M, lo, hi);

  Eigen::MatrixXd V;
  Eigen::VectorXd S;

  if (Large) {
    Eigen::BDCSVD<Eigen::MatrixXd> svd(mm, Eigen::ComputeThinU | Eigen::ComputeThinV);

    V = svd.matrixV();
    S = svd.singularValues();

  } else {

    Eigen::JacobiSVD<Eigen::MatrixXd> svd(mm, Eigen::ComputeThinU | Eigen::ComputeThinV);

    V = svd.matrixV();
    S = svd.singularValues();
  }

  for (int i=0; i<M; i++) {
    ev[i+1] = S[i];
    for (int j=0; j<n; j++) ef[i+1][j+lo] = V(j, i);
  }

  return ev;
}

