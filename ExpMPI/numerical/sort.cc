
#include <stdlib.h>
#include <Vector.h>


/* sort a vector */


void Vector::Sort(void)
{
	void heapsort(int, double *);
	double *a;
	int i;

	a = array(1, high-low+1);
	heapsort(high-low+1, a);

	for (i=low; i<=high; i++) elements[i] = a[1 + i-low];
}
		


void heapsort(int n, double *ra)
{
	int l,j,ir,i;
	double rra;

	l=(n >> 1)+1;
	ir=n;
	for (;;) {
		if (l > 1)
			rra=ra[--l];
		else {
			rra=ra[ir];
			ra[ir]=ra[1];
			if (--ir == 1) {
				ra[1]=rra;
				return;
			}
		}
		i=l;
		j=l << 1;
		while (j <= ir) {
			if (j < ir && ra[j] < ra[j+1]) ++j;
			if (rra < ra[j]) {
				ra[i]=ra[j];
				j += (i=j);
			}
			else j=ir+1;
		}
		ra[i]=rra;
	}
}
	




