#include <stdio.h>
#include <math.h>
#include <numerical.h>
#include <Vector.h>

/*
4th order sympletic integration algorithm from
Candy and Rozmus, J, Comp. Phys. 92, 230 (1991)
*/


#define A1 0.6756035959798288170238458065678177873714
#define A4 0.6756035959798288170238458065678177873714
#define A2 -0.1756035959798288170238458065678177873714
#define A3 -0.1756035959798288170238458065678177873714
#define B1 0.0
#define B2 1.3512071919596576340477441177810572122634
#define B4 1.3512071919596576340477441177810572122634
#define B3 -1.7024143839193152680951301399452727444042


static double  *ql=NULL, *qc, *pl, *pc, *f;
// static Vector A(1, 4);
// static Vector B(1, 4);
static double A[5];
static double B[5];
static int sia4_first=1;

void sia4_init(void)
{

  sia4_first=0;

  ql = nr_vector(1,3);
  qc = nr_vector(1,3);
  pl = nr_vector(1,3);
  pc = nr_vector(1,3);
  f  = nr_vector(1,3);

  A[1] = A1;
  A[2] = A2;
  A[3] = A3;
  A[4] = A4;

  B[1] = B1;
  B[2] = B2;
  B[3] = B3;
  B[4] = B4;

}

void sia4(
	  double *x,
	  double *v,
	  double *x1,
	  double *v1,
	  double t,
	  double h,
	  symp_derivs derivs)
{
  int i, j;

  if (sia4_first) sia4_init();

				// zeroth order step

  for (i=1; i<=3; i++) {
    ql[i] = x[i];
    pl[i] = v[i];
  }

  (*derivs)(t, ql, pl,f);
	
				// j = 1, ... , 4

  for (j=1; j<=4; j++) {


    for (i=1; i<=3; i++)
      pc[i] = pl[i] + B[j]*h*f[i];
    
    for (i=1; i<=3; i++)
      qc[i] = ql[i] + A[j]*h*pc[i];
    
    t += A[j]*h;

    (*derivs)(t, qc, pc, f);
    
    for (i=1; i<=3; i++) {
      ql[i] = qc[i];
      pl[i] = pc[i];
    }
  }

  for (i=1; i<=3; i++) {
    x1[i] = ql[i];
    v1[i] = pl[i];
  }
  
}


