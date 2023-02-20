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


#include <cstdlib>
#include <iostream>
#include <iomanip>

#include "cpoly.H"

extern "C" {
  void fpeinit(int);
}

int parse_args(int, char **);
Eigen::VectorXcd Cget_horner(std::complex<double> z, CPoly &p);
std::complex<double> compute_rat_integral(double a, double b,Eigen::VectorXcd& x_data, Eigen::VectorXcd& y_data);
std::complex<double> integral_pole_contribution(double a, double b, std::complex<double> root, int mult,
				    CPoly &n, CPoly &d);
void Csyn_poly_div(CPoly &p1, CPoly &p2, CPoly &Quotient, CPoly &Remainder);
void make_tableau(Eigen::VectorXcd& x, Eigen::VectorXcd& y, Eigen::MatrixXcd &m);
void print_tableau(Eigen::MatrixXcd& m);
std::complex<double> get_rat(std::complex<double> x, Eigen::VectorXcd& numer, Eigen::VectorXcd& denom);
void gen_coefs(CPoly& numer, CPoly& denom, Eigen::VectorXcd& x_data, Eigen::MatrixXcd& m);
void renorm_coef(Eigen::VectorXcd& x1, Eigen::VectorXcd& y1, Eigen::VectorXcd& x2, Eigen::VectorXcd& y2);
void strike_entry1(Eigen::VectorXcd& rrn, int i, int& irn);
Eigen::VectorXcd get_multiplicity(Eigen::VectorXcd &rrn, int &dim, std::vector<int>& nmult);
void zroots(Eigen::VectorXcd& a, Eigen::VectorXcd& roots, int polish);

static double DEPS=0.01;
const double SMALL=1.0e-12;
static int IDBG=0;
static CPoly monad(1);

void rat_integral_set_parameters(double eps, int idbg)
{
  DEPS = eps;
  IDBG = idbg;
}


std::complex<double> compute_rat_integral(double a, double b,
					  Eigen::VectorXcd& x_data,
					  Eigen::VectorXcd& y_data)
{
  int sz = x_data.size();
  
  // Check for constant (including zero) integrand
  
  int i = 0;
  std::complex<double> val = y_data[i++];
  while(i<sz) {
    if (std::abs(val-y_data[i++]) > SMALL) break;
  }
  if (i>y_data.size()) {
    return val * (b-a);
  }
  
  Eigen::MatrixXcd m(sz, sz);
  make_tableau(x_data, y_data, m);
  
  // print_tableau(m);
  
  CPoly Numer, Denom;
  gen_coefs(Numer, Denom, x_data, m);
  
  Eigen::VectorXcd rd(Denom.size());
  zroots(Denom, rd, 0);
  zroots(Denom, rd, 1);
  
  if (IDBG) {
    std::cout << "\nNumerator coefficients:" << std::endl;
    for (int i=0; i<Numer.size(); i++) 
      std::cout << std::setw(10) << i
		<< std::setw(16) << Numer[i].real()
		<< std::setw(16) << Numer[i].imag()
		<< std::endl;
    
    std::cout << std::endl << "Denomator coefficients:" << std::endl;
    for (int i=0; i<Denom.size(); i++)
      std::cout << std::setw(10) << i
		<< std::setw(16) << Denom[i].real()
		<< std::setw(16) << Denom[i].imag()
		<< std::endl;
    
    std::cout << std::endl << "Roots:" << std::endl;
    for (int i=0; i<rd.size(); i++)
      std::cout << std::setw(10) << i
		<< std::setw(16) << rd[i].real()
		<< std::setw(16) << rd[i].imag()
		<< std::endl;
  }
  
  int rdsize = rd.size();
  std::vector<int> dmult(rdsize);
  Eigen::VectorXcd rrdu = get_multiplicity(rd, rdsize, dmult);
  
  // For principal value evaluation
  
  CPoly numr;
  CPoly ply;
  Csyn_poly_div(Numer, Denom, ply, numr);
  
  // Evaluate integral
  
  std::complex<double> ans = 0.0;
  for (int i=0; i<rrdu.size(); i++) {
    ans += integral_pole_contribution(a, b, rrdu[i], dmult[i], numr, Denom);
    if (IDBG) {
      std::cout << std::setw(4) << i
		<< std::setw(16) << ans.real()
		<< std::setw(16) << ans.imag()
		<< std::endl;
    }
  }
  
  monad[0] = 0.0;
  monad[1] = 1.0;
  ply &= monad;
  for (i=1; i<=ply.getorder(); i++) ply[i] /= i;
  ans += ply.eval(b) - ply.eval(a);
  
  return ans;
}



std::complex<double> integral_pole_contribution(double a, double b, std::complex<double> root, int mult,
				    CPoly &n, CPoly &d)
{
  int i;
  std::complex<double> t1, t2, ans=0.0, c1, c2;
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
    
    for (int i=mult-1; i>0; i--) {
      c1 *= t1;
      c2 *= t2;
      ans += q[mult-i]*(c2 - c1)/static_cast<double>(-i);
    }
  }
  
  ans += q[0]*log(t2/t1);
  
  return ans;
}


Eigen::VectorXcd get_multiplicity(Eigen::VectorXcd &rrn, int &dim,
				  std::vector<int> &nmult)
{
  int icnt;
  
  Eigen::VectorXcd tmp = rrn;
  int irn = dim;
  
  for (int i=0; i<irn; i++) {
    
    icnt = 1;
    if (i!=irn) {
      for (int j=i+1; j<irn; j++) {
	if ( std::abs(tmp[i] - tmp[j]) < DEPS*(std::abs(tmp[i])+std::abs(tmp[j])) ) {
	  strike_entry1(tmp, j, irn);
	  icnt++;
	}
      }
    }
    nmult[i] = icnt;
  }
  
  Eigen::VectorXcd ans(irn);
  for (int i=0; i<irn; i++) ans[i] = tmp[i];
  
  return ans;
}


void strike_entry1(Eigen::VectorXcd& rrn, int i, int& irn)
{
  for (int k=i; k<irn; k++) rrn[k] = rrn[k+1];
  irn--;
}

// Divided differences

void make_tableau(Eigen::VectorXcd& x, Eigen::VectorXcd& y, Eigen::MatrixXcd &m)
{
  int sz = x.size()-1;
  
  m.setZero();
  
  for (int j=0; j<sz; j++)
    m(0, j) = y[j+1];
  
  for (int i=1; i<sz; i++) {
    for (int j=i; j<sz; j++)
      m(i, j) = (x[i] - x[j+1])/(m(i-1, i-1) - m(i-1, j));
  }
  
}

void print_tableau(Eigen::MatrixXcd& m)
{
  int sz = m.rows();
  
  for (int j=0; j<sz; j++) {
    for (int i=0; i<=j; i++)
      std::cout << "(" << std::setw(16) << m(i, j).real()
		<< "," << std::setw(16) << m(i, j).imag()
		<< ")";
    std::cout << std::endl;
  }
}

std::complex<double> gen_fct(std::complex<double> x, Eigen::VectorXcd& x_data, Eigen::MatrixXcd& m)
{
  int sz = x_data.size();
  
  // Recursion
  
  std::complex<double> q, ql, qll, p, pl, pll;
  
  int cur = 0;		// Initialization
  
  pll = m(cur, cur);
  qll = 1.0;
  cur++;
  
  pl = pll*m(cur, cur) + x - x_data(cur);
  ql = m(cur, cur);
  cur++;
  
  
  while (cur<sz) {
    // Get next Q
    q = m(cur, cur)*ql + (x - x_data[cur])*qll;
    
    // Get next P
    p = m(cur, cur)*pl + (x - x_data[cur])*pll;
    
    qll = ql;
    ql = q;
    pll = pl;
    pl = p;
    cur++;
  }
  
  return p/q;
}

void gen_coefs(CPoly& numer, CPoly& denom,
	       Eigen::VectorXcd& x_data, Eigen::MatrixXcd& m)
{
  
  // Initialize coefficient vectors
  int max = m.size();
  int ndim = (max+1)/2;
  int ddim = max/2;
  numer = CPoly(ndim);
  denom = CPoly(ddim);
  
  Eigen::VectorXcd numer_l = numer;
  Eigen::VectorXcd numer_ll = numer;
  Eigen::VectorXcd denom_l = denom;
  Eigen::VectorXcd denom_ll = denom;
  
  // Begin recursion
  
  int i;
  int cur = 0;
  
  numer_ll[0] = m(cur, cur);	
  denom_ll[0] = 1.0;
  cur++;
  
  numer_l[0] = numer_ll[0]*m(cur, cur) - x_data[cur];
  numer_l[1] = 1.0;
  denom_l[0] = m(cur, cur);
  cur++;
  
  while (cur<=max) {
    // Get next Q
    
    denom[0] = m(cur, cur)*denom_l[0] - x_data[cur]*denom_ll[0];
    for (i=1; i<=cur/2; i++)
      denom[i] = m(cur, cur)*denom_l[i] + denom_ll[i-1] - 
	x_data[cur]*denom_ll[i];
    
    // Get next P
    
    numer[0] = m(cur, cur)*numer_l[0] - x_data[cur]*numer_ll[0];
    for (i=1; i<=(cur+1)/2; i++)
      numer[i] = m(cur, cur)*numer_l[i] + numer_ll[i-1] - 
	x_data[cur]*numer_ll[i];
    
    denom_ll = denom_l;
    denom_l  = denom;
    numer_ll = numer_l;
    numer_l  = numer;
    
    cur++;
    renorm_coef(numer_l, denom_l, numer_ll, denom_ll);
  }
  
}

void renorm_coef(Eigen::VectorXcd& x1, Eigen::VectorXcd& y1, Eigen::VectorXcd& x2, Eigen::VectorXcd& y2)
{
  int sizx = x1.size();
  int sizy = y1.size();
  
  double ymax = std::abs(y1[0]);
  double ytmp;
  for (int i=1; i<sizy; i++) {
    ytmp = std::abs(y1[i]);
    if (ytmp <= 1.0e-10) continue;
    if (ytmp > ymax) ymax = ytmp;
  }
  
  for (int i=0; i<sizx; i++) {
    x1[i] /= ymax;
    x2[i] /= ymax;
  }
  for (int i=0; i<sizy; i++) {
    y1[i] /= ymax;
    y2[i] /= ymax;
  }
}

