// This may look like C code, but it is really -*- C++ -*-

#ifndef _EmpOrth_h

#define _EmpOrth_h 1

#include <WghtOrth3.h>
#include <math.h>


class EmpCylSL : public WghtOrtho
{
private:
  int MMAX;
  int NORDER;
  int NKEEP;
  double ZMAX;

  double YMIN, YMAX;
  double dX, dY, dk, facZ;
  int M, NUMR, cylused, cylused1;
  Vector zn;

  LegeQuad *lwR, *lwZ;
  Matrix *facC, *facS;
  Vector** accum_cos0;
  Vector** accum_sin0;
  
  int rank1, rank2, rank3;

  double *MPIin, *MPIout;
  double *MPIin_eof, *MPIout_eof;
  double *mpi_double_buf1, *mpi_double_buf2, *mpi_double_buf3;
  int MPIbufsz, MPItable;
  Logic MPIset, MPIset_eof;
  MPI_Status status;

  double**** SC;
  double**** SS;

  Matrix var;

  Vector ev;
  Matrix ef;

  CylindricalSL *ortho;

  Matrix** potC;
  Matrix** densC;
  Matrix** rforceC;
  Matrix** zforceC;

  Matrix** potS;
  Matrix** densS;
  Matrix** rforceS;
  Matrix** zforceS;

  Matrix** table;

  Matrix* tpot;
  Matrix* tdens;
  Matrix* trforce;
  Matrix* tzforce;

  Vector* accum_cos;
  Vector* accum_sin;
  Vector* accum_cos2;
  Vector* accum_sin2;
  Matrix *vc, *vs;

  Matrix tabp, tabf, tabd;

  Logic coefs_made;
  Logic eof_recompute;
  Logic eof_made;

  void make_grid();
  void send_eof_grid();
  void receive_eof(int request_id, int m, int ir);
  void compute_eof_grid(int request_id, int m, int ir);
  void setup_eof_grid(void);
  void setup_eof(void);
  void accumulate_eof(double r, double z, double phi, double mass, int id);
  void make_eof(void);
				// 1=write, 0=read
				// return: 0=failure
  int cache_grid(int);		
  double integral(int, int, int, int);
  void get_pot(Matrix&, Matrix&, double, double);
  void pca_hall(void);

public:

				// Global parameters


  static Logic DENS;
  static Logic SELECT;
  static int NZOF;
  static int NINTR;
  static int NINTZ;
  static int NUMX;
  static int NUMY;
  static int NOUT;
  static double RMAX;
  static double HSCALE;
  static String CACHEFILE;
  static String TABLEFILE;


				// Constructors
  EmpCylSL(void);
  EmpCylSL(CylindricalSL *ortho, double zmax, int nord);
  ~EmpCylSL(void);

  void reset(CylindricalSL *ortho, double zmax, int nord);
  int read_cache(void);

				// Parameter access
  int get_order(void) {return NORDER;}
  CylindricalSL* SL(void) {return ortho;}

				// Z coordinate transformation

  /*
  inline double z_to_y(double z) { return asinh(z/HSCALE); }
  inline double y_to_z(double y) { return HSCALE*sinh(y); }
  */
  double z_to_y(double z) { return z/(fabs(z)+DBL_MIN)*asinh(fabs(z)/HSCALE); }
  double y_to_z(double y) { return HSCALE*sinh(y); }

				// Main member functions

  void compute_eof(void) { eof_recompute = TRUE; }

  void get_all(int m, int ir, int n, double r, double z, double phi,
	       double& p, double& d, double& fr, double& fz, double& fp);

  void setup_accumulation(void);
  void make_coefficients(void);

  void accumulate(double r, double z, double phi, double mass, int id);
  void accumulated_eval(double r, double z, double phi,
			double& p, double& fr, double& fz, double& fp);
  double accumulated_dens_eval(double r, double z, double phi);

  
  void dump_coefs(FILE *fout);
  void dump_coefs_binary_last(FILE *fout, double time);
  void dump_coefs_binary_curr(FILE *fout, double time);
  void dump_basis(char* name, int step);

};

#endif
