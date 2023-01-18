#include <stdio.h>
#include <math.h>
#include <numerical.H>

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


static Eigen::VectorXd ql, qc, pl, pc, f;
static Eigen::Vector4d A, B;
static int sia4_first=1;

void sia4_init(void)
{

  sia4_first=0;

  A[0] = A1;
  A[1] = A2;
  A[2] = A3;
  A[3] = A4;

  B[0] = B1;
  B[1] = B2;
  B[2] = B3;
  B[3] = B4;

  ql.resize(3);
  qc.resize(3);
  pl.resize(3);
  pc.resize(3);
  f .resize(3);
}

void sia4(Eigen::VectorXd &x,
	  Eigen::VectorXd &v,
	  Eigen::VectorXd &x1,
	  Eigen::VectorXd &v1,
	  double t,
	  double h,
	  symp_derivs derivs)
{
  if (sia4_first) sia4_init();
				// zeroth order step

  for (int i=0; i<3; i++) {
    ql[i] = x[i];
    pl[i] = v[i];
  }

  derivs(t, ql, pl, f);

  for (int j=0; j<4; j++) {


    for (int i=0; i<3; i++)
      pc[i] = pl[i] + B[j]*h*f[i];
    
    for (int i=0; i<3; i++)
      qc[i] = ql[i] + A[j]*h*pc[i];
    
    t += A[j]*h;

    derivs(t, qc, pc, f);
    
    for (int i=0; i<3; i++) {
      ql[i] = qc[i];
      pl[i] = pc[i];
    }
  }

  for (int i=0; i<3; i++) {
    x1[i] = ql[i];
    v1[i] = pl[i];
  }
  
}

