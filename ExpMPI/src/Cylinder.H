#ifndef _Cylinder_H
#define _Cylinder_H

#include <Vector.h>

#include <EmpOrth9thd.h>

#include <Orient.H>
#include <Basis.H>

//! This routine computes the potential, acceleration and density using the Cylindrical biorthogonal expansion
class Cylinder : public Basis
{
private:

  Orient* orient;
  EmpCylSL *ortho;
  int eof, ncompcyl;
  double cyltime, cylmass;

  // Three_Vector *pos, *frc;
  Vector *pos, *frc;

  double *cylmass0;
  int *offgrid;

  void initialize(void);

				// Parameters
  double rcylmin, rcylmax, zmax, acyl;
  int nmax, lmax, mmax;
  int ncylnx, ncylny;
  double hcyl;
  int ncylorder, ncylrecomp, hallfreq;
  string hallfile;
  bool self_consistent, logarithmic, density, selector, firstime;

  // These should be ok for all derived classes, hence declared private

  void determine_coefficients();

  void determine_acceleration_and_potential();

  void * determine_coefficients_thread(void * arg);

  void * determine_acceleration_and_potential_thread(void * arg);

public:

  //! Mutexes for multithreading
  static pthread_mutex_t used_lock, cos_coef_lock, sin_coef_lock;

  //! Constructor
  /** \param line contains the parameters passed to the method */
  Cylinder(string& line);

  //! Destructor
  virtual ~Cylinder();
  
  //! The main force call
  void get_acceleration_and_potential(Component*);

  //! Return the value for the fields at a spherical polar coordinate
  void 
  determine_fields_at_point_sph(double r, double theta, double phi,
				double *tdens0, double *tpotl0, 
				double *tdens, double *tpotl, 
				double *tpotr, double *tpott, 
				double *tpotp);

  //! Return the value for the fields at a cylindrical coordinate
  void 
  determine_fields_at_point_cyl(double r, double z, double phi,
				double *tdens0, double *tpotl0, 
				double *tdens, double *tpotl, 
				double *tpotr, double *tpotz, 
				double *tpotp);

  //! Save coefficients to file (need type marker to id dump, component id?)
  void dump_coefs(ostream& out);

  //! Sanity check on grid: dumps SM-style images of initial field
  void dump_mzero(const string& name, int step);
};


#endif
