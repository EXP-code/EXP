// This may look like C code, but it is really -*- C++ -*-

#include <Vector.h>
#include <logic.h>
#include <fft.h>

class WghtOrtho
{
public:

				// Error function

  void bomb(const char *s) {
    cerr << "WghtOrtho ERROR: " << s << '\n';
    exit(-1);
  }

};


class CylindricalCB : public WghtOrtho
{
private:
  double A, H, ZMAX, dk, dZ;
  int MAXR, NFFT;
  int M, n, n2;

  double *MPIin, *MPIout;
  Logic MPIset;

  Matrix *ortho, *U, *UT;
  Matrix output;
  Matrix Z;

  Matrix accum_cos;
  Matrix accum_sin;
  Matrix coscoef;
  Matrix sincoef;
  Logic coefs_made;

  Logic initialized;
  void initialize(void);
  int NINT;

  FFT *fft;

  const double defZMAX=100.0;
  const int defNINT=200;

  class PotTable {

  private:

    int N, N2, NR;
    double RMAX;
    double dR;
    CylindricalCB *ref;

    Vector** fr;
    Vector** fz;

  public:

    void make_table(int, double, CylindricalCB*);
    
    PotTable(void) : N(0) {};
    PotTable(int n, double rmax, CylindricalCB* p) {make_table(n, rmax, p);}
    ~PotTable(void);

    void feval(double r, double z, double phi,
	  double& p, double& fr, double& fz, 
	  double &fp);
  };

  Logic use_tables;
  int NPT;
  double NPTRMAX;
  PotTable pt;


public:

  CylindricalCB(void);
  CylindricalCB(int mu, int nfft, int m, double a, double h);
  ~CylindricalCB(void);

  void reset(int mu, int nfft, int m, double a, double h);
  void set_num(int nint) { NINT=nint; }
  void set_zmax(double zmax) { ZMAX=zmax; }
  void set_table(Logic onoff, int npt=0, double nptrmax=0.0) { 
    use_tables=onoff;
    NPT = npt;
    NPTRMAX = nptrmax;
  }

  double get_pot(int mu, int nu, double r, double z);
  double get_dens(int mu, int nu, double r, double z);

  Matrix inner_product_dens(double (*fct)(double, double));
  Matrix inner_product_pot(double (*fct)(double, double));

  double pot_eval(Matrix& coef, double r, double z);
  double dens_eval(Matrix& coef, double r, double z);
  double r_force_eval(Matrix& coef, double r, double z);
  double z_force_eval(Matrix& coef, double r, double z);
  void force_eval(Matrix& coef, double r, double z, double& fr, double& fz);

  double norm(int mu, int nu);

  void setup_accumulation(void);
  void accumulate(double r, double z, double phi, double mass);
  Matrix accumulated_coefficients(int sincos);
  void make_coefficients(void);
  void accumulated_eval(double r, double z, double phi,
			double& p, double& fr, double& fz, double& fp);
  double accumulated_dens_eval(double r, double z, double phi);
  
  const Matrix& get_ortho(int iz) { return ortho[iz]; }

};

