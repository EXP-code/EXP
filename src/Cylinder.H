#ifndef _Cylinder_H
#define _Cylinder_H

#include <Vector.h>
#include <Orient.H>
#include <Basis.H>
#include <EmpOrth9thd.h>

#include "config.h"

#if HAVE_LIBCUDA==1
#include <cudaParticle.cuH>
#include <cudaMappingConstants.cuH>
#endif
#include <NVTX.H>

class MixtureBasis;

/** This routine computes the potential, acceleration and density 
    using the Cylindrical biorthogonal expansion 
*/
class Cylinder : public Basis
{
private:

  bool expcond;
  int rnum, pnum, tnum;
  double ashift;
  unsigned int vflag;

  EmpCylSL *ortho;
  int eof, ncompcyl;
  double cylmass;

  // Three_Vector *pos, *frc;
  Vector *pos, *frc;

  double *cylmass0;
  int *offgrid;

  void initialize(void);

				// Parameters
  double rcylmin, rcylmax, zmax, acyl;
  int nmax, lmax, mmax;
  int ncylnx, ncylny, ncylr;
  double hcyl;
  int ncylorder, ncylrecomp, hallfreq, vtkfreq;
  string eof_file;
  bool self_consistent, logarithmic, density, pca, pcainit, pcavtk;
  bool cmap, try_cache, firstime, dump_basis;

  // These should be ok for all derived classes, hence declared private

  void determine_coefficients();

  void determine_coefficients_eof();

  void determine_acceleration_and_potential();

  void * determine_coefficients_thread(void * arg);

  void * determine_acceleration_and_potential_thread(void * arg);

  /** Extrapolate and sum coefficents per multistep level to get
      a complete set of coefficients for force evaluation at an
      intermediate time step
  */
  void compute_multistep_coefficients() 
  {
    nvTracerPtr tPtr;
    if (cuda_prof)
      tPtr = nvTracerPtr(new nvTracer("Cylinder::compute_multistep_coefficients"));
    ortho->compute_multistep_coefficients(mlevel);
  }

  //! Reset used particle counter and mass registers
  virtual void multistep_reset();

  MixtureBasis *mix;


protected:

  //! CUDA method for coefficient accumulation
#if HAVE_LIBCUDA==1
  virtual void determine_coefficients_cuda();
  virtual void determine_acceleration_cuda();

  thrust::host_vector<double>   host_coefs;
  thrust::device_vector<double> dev_coefs;

  virtual void HtoD_coefs();
  virtual void DtoH_coefs(int mlevel);

  std::vector<cudaArray_t> cuInterpArray;
  thrust::host_vector<cudaTextureObject_t> tex;

  virtual void initialize_cuda()
  {
    ortho->initialize_cuda(cuInterpArray, tex);
  }

  virtual cudaMappingConstants getCudaMappingConstants()
  {
    return ortho->getCudaMappingConstants();
  }


  virtual void destroy_cuda();

  virtual void initialize_mapping_constants();

  //! Compare CPU and GPU force evaluation for debugging
  void host_dev_force_compare();
#endif

public:

  //! Mutexes for multithreading
  //@{
  static pthread_mutex_t used_lock, cos_coef_lock, sin_coef_lock;
  //@}

  //! Constructor
  //! \param line contains the parameters passed to the method
  //! \param m is the MixtureBasis for a multicenter expansion
  Cylinder(string& line, MixtureBasis *m=0);

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

  /** Update the multi time step coefficient table when moving particle 
      <code>i</code> from level <code>cur</code> to level 
      <code>next</code>.
  */
  //@{
  virtual void multistep_update_begin() { 
    if (self_consistent) ortho->multistep_update_begin(); }

  virtual void multistep_update(int cur, int next, Component* c, int i, int id);
  virtual void multistep_update_finish() { 
    if (self_consistent) ortho->multistep_update_finish(); }
  //@}

  //! Print debug info
  virtual void multistep_debug();

  //! Save coefficients to file (need type marker to id dump, component id?)
  void dump_coefs(ostream& out);

  //! Sanity check on grid: dumps SM-style images of initial field
  void dump_mzero(const string& name, int step);

};


#endif
