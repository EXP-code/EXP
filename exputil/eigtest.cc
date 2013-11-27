#include <math.h>
#include <unistd.h>
#include <stdlib.h>
#include <numerical.h>
#include <Vector.h>



int main()
{
	Matrix M(1, 3, 1, 3);
	Vector e(1, 3);

	M[1][1] = 1.0;
	M[1][2] = 0.5;
	M[1][3] = 0.1;

	M[2][2] = 1.0;
	M[2][3] = 0.2;

	M[3][3] = 1.0;

	M[2][1] = M[1][2];
	M[3][1] = M[1][3];
	M[3][2] = M[2][3];

	e = M.Symmetric_Eigenvalues();
	e.print(stdout);
}
	



