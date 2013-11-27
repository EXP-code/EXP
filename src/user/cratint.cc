// This is really C++ code
/*****************************************************************************
 *  Description:
 *  -----------
 *
 *  This routine tests the inverse difference routine
 *
 *
 *  Call sequence:
 *  -------------
 *
 *  Parameters:
 *  ----------
 *
 *
 *  Returns:
 *  -------
 *
 *
 *  Notes:
 *  -----
 *
 *
 *  By:
 *  --
 *
 *  MDW 11/20/91
 *
 ***************************************************************************/


#include <stdlib.h>
#include <stdio.h>
#include <Vector.h>
#include <kevin_complex.h>
#include "cpoly.h"

extern "C" {
  void fpeinit(int);
}

int parse_args(int, char **);
CVector Cget_horner(KComplex z, CPoly &p);
KComplex compute_rat_integral(double a, double b,CVector& x_data, CVector& y_data);
KComplex integral_pole_contribution(double a, double b, KComplex root, int mult,
				    CPoly &n, CPoly &d);
void Csyn_poly_div(CPoly &p1, CPoly &p2, CPoly &Quotient, CPoly &Remainder);
void make_tableau(CVector& x, CVector& y, CMatrix &m);
void print_tableau(CMatrix& m);
KComplex get_rat(KComplex x, CVector& numer, CVector& denom);
void gen_coefs(CPoly& numer, CPoly& denom, CVector& x_data, CMatrix& m);
void renorm_coef(CVector& x1, CVector& y1, CVector& x2, CVector& y2);
void strike_entry1(CVector& rrn, int i, int& irn);
CVector get_multiplicity(CVector &rrn, int &dim, int *nmult);
void zroots(CVector& a, CVector& roots, int polish);

static double DEPS=0.01;
const double SMALL=1.0e-12;
static int IDBG=0;
static CPoly monad(1);

void rat_integral_set_parameters(double eps, int idbg)
{
  DEPS = eps;
  IDBG = idbg;
}


KComplex compute_rat_integral(double a, double b, CVector& x_data, CVector& y_data)
{
  int i;
  
  int min = x_data.getlow();
  int max = x_data.gethigh();
  
  // Check for constant (including zero) integrand
  
  i = min;
  KComplex val = y_data[i++];
  while(i<=max) {
    if (fabs(val-y_data[i++]) > SMALL) break;
  }
  if (i>y_data.gethigh()) {
    return val * (b-a);
  }
  
  CMatrix m(min-1, max-1, min-1, max-1);
  make_tableau(x_data, y_data, m);
  
  // print_tableau(m);
  
  CPoly Numer, Denom;
  gen_coefs(Numer, Denom, x_data, m);
  
  CVector rd(1, Denom.gethigh());
  zroots(Denom, rd, 0);
  zroots(Denom, rd, 1);
  
  if (IDBG) {
    printf("\nNumerator coefficients:\n");
    for (i=Numer.getlow(); i<=Numer.gethigh(); i++)
      printf("%d  %e  %e\n", i, Numer[i].real(), Numer[i].imag());
    
    printf("\nDenomator coefficients:\n");
    for (i=Denom.getlow(); i<=Denom.gethigh(); i++)
      printf("%d  %e  %e\n", i, Denom[i].real(), Denom[i].imag());
    
    printf("\nRoots:\n");
    for (i=1; i<=rd.gethigh(); i++)
      printf("%d  %e  %e\n", i, rd[i].real(), rd[i].imag());
    
  }
  
  int *dmult = new int[rd.gethigh()] - 1;
  int rdmax = rd.gethigh();
  CVector rrdu = get_multiplicity(rd, rdmax, dmult);
  
  // For principal value evaluation
  
  CPoly numr;
  CPoly ply;
  Csyn_poly_div(Numer, Denom, ply, numr);
  
  // Evaluate integral
  
  KComplex ans = 0.0;
  for (i=1; i<=rrdu.gethigh(); i++) {
    ans += integral_pole_contribution(a, b, rrdu[i], dmult[i], numr, Denom);
    if (IDBG) {
      printf("%2d** ",i);  
      ans.print();
    }
  }
  delete [] (dmult+1);
  
  monad[0] = 0.0;
  monad[1] = 1.0;
  ply &= monad;
  for (i=1; i<=ply.getorder(); i++) ply[i] /= i;
  ans += ply.eval(b) - ply.eval(a);
  
  return ans;
}



KComplex integral_pole_contribution(double a, double b, KComplex root, int mult,
				    CPoly &n, CPoly &d)
{
  int i;
  KComplex t1, t2, ans=0.0, c1, c2;
  CPoly denom = d;
  
  monad[0] = -root;
  monad[1] = 1.0;
  
  for (i=1; i<=mult; i++) {
    denom %= monad;
    denom.Pop(1);
  }
  
  CPoly q = CPoly(Cget_horner(root, n)) % CPoly(Cget_horner(root, denom));
  
  t1 = a-root;
  t2 = b-root;
  
  if (mult > 1) {
    c1 = pow(t1, -mult);
    c2 = pow(t2, -mult);
    
    for (i=mult; i>1; i--) {
      c1 *= t1;
      c2 *= t2;
      ans += q[mult-i]*(c2 - c1)/(1-i);
    }
  }
  
  ans += q[0]*log(t2/t1);
  
  return ans;
}


CVector get_multiplicity(CVector &rrn, int &dim, int *nmult)
{
  int i, j, icnt;
  
  CVector tmp = rrn;
  int irn = dim;
  
  for (i=1; i<=irn; i++) {
    
    icnt = 1;
    if (i!=irn) {
      for (j=i+1; j<=irn; j++) {
	if ( fabs(tmp[i] - tmp[j]) < DEPS*(fabs(tmp[i])+fabs(tmp[j])) ) {
	  strike_entry1(tmp, j, irn);
	  icnt++;
	}
      }
    }
    nmult[i] = icnt;
  }
  
  CVector ans(1,irn);
  for (i=1; i<=irn; i++) ans[i] = tmp[i];
  
  return ans;
}


void strike_entry1(CVector& rrn, int i, int& irn)
{
  int k;
  
  for (k=i; k<irn; k++) rrn[k] = rrn[k+1];
  irn--;
}

// Divided differences

void make_tableau(CVector& x, CVector& y, CMatrix &m)
{
  int i,j;
  int min = x.getlow()-1;
  int max = x.gethigh()-1;
  
  m.zero();
  
  for (j=min; j<=max; j++)
    m[min][j] = y[j+1];
  
  for (i=min+1; i<=max; i++) {
    for (j=i; j<=max; j++)
      m[i][j] = (x[i] - x[j+1])/(m[i-1][i-1] - m[i-1][j]);
  }
  
}

void print_tableau(CMatrix& m)
{
  int i,j;
  int min = m.getrlow();
  int max = m.getrhigh();
  
  for (j=min; j<=max; j++) {
    for (i=min; i<=j; i++)
      printf("  (%15.6e, %15.6e)",m[i][j].real(), m[i][j].imag());
    printf("\n");
  }
}

KComplex gen_fct(KComplex x, CVector& x_data, CMatrix& m)
{
  int min = x_data.getlow();
  int max = x_data.gethigh();
  
  // Recursion
  
  KComplex q, ql, qll, p, pl, pll;
  
  int cur = min;		// Initialization
  
  pll = m[cur][cur];
  qll = 1.0;
  cur++;
  
  pl = pll*m[cur][cur] + x - x_data[cur];
  ql = m[cur][cur];
  cur++;
  
  
  while (cur<=max) {
    // Get next Q
    q = m[cur][cur]*ql + (x - x_data[cur])*qll;
    
    // Get next P
    p = m[cur][cur]*pl + (x - x_data[cur])*pll;
    
    qll = ql;
    ql = q;
    pll = pl;
    pl = p;
    cur++;
  }
  
  return p/q;
}

void gen_coefs(CPoly& numer, CPoly& denom, CVector& x_data, CMatrix& m)
{
  
  // Initialize coefficient vectors
  int max = m.getrhigh();
  int ndim = (max+1)/2;
  int ddim = max/2;
  numer = CPoly(ndim);
  denom = CPoly(ddim);
  
  CVector numer_l = numer;
  CVector numer_ll = numer;
  CVector denom_l = denom;
  CVector denom_ll = denom;
  
  // Begin recursion
  
  int i;
  int cur = 0;
  
  numer_ll[0] = m[cur][cur];	
  denom_ll[0] = 1.0;
  cur++;
  
  numer_l[0] = numer_ll[0]*m[cur][cur] - x_data[cur];
  numer_l[1] = 1.0;
  denom_l[0] = m[cur][cur];
  cur++;
  
  while (cur<=max) {
    // Get next Q
    
    denom[0] = m[cur][cur]*denom_l[0] - x_data[cur]*denom_ll[0];
    for (i=1; i<=cur/2; i++)
      denom[i] = m[cur][cur]*denom_l[i] + denom_ll[i-1] - 
	x_data[cur]*denom_ll[i];
    
    // Get next P
    
    numer[0] = m[cur][cur]*numer_l[0] - x_data[cur]*numer_ll[0];
    for (i=1; i<=(cur+1)/2; i++)
      numer[i] = m[cur][cur]*numer_l[i] + numer_ll[i-1] - 
	x_data[cur]*numer_ll[i];
    
    denom_ll = denom_l;
    denom_l  = denom;
    numer_ll = numer_l;
    numer_l  = numer;
    
    cur++;
    renorm_coef(numer_l, denom_l, numer_ll, denom_ll);
  }
  
}

void renorm_coef(CVector& x1, CVector& y1, CVector& x2, CVector& y2)
{
  int minx = x1.getlow();
  int maxx = x1.gethigh();
  int miny = y1.getlow();
  int maxy = y1.gethigh();
  
  int i;
  
  double ymax = fabs(y1[miny]);
  double ytmp;
  for (i=miny+1; i<=maxy; i++) {
    ytmp = fabs(y1[i]);
    if (ytmp <= 1.0e-10) continue;
    if (ytmp > ymax) ymax = ytmp;
  }
  
  //  printf("Ymax=%e\n",ymax);
  
  for (i=minx; i<=maxx; i++) {
    x1[i] /= ymax;
    x2[i] /= ymax;
  }
  for (i=miny; i<=maxy; i++) {
    y1[i] /= ymax;
    y2[i] /= ymax;
  }
}

