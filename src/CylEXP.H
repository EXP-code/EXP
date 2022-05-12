#ifndef _CylEXP_H
#define _CylEXP_H

#include <global.H>
#include <expand.H>
#include <EmpCylSL.H>

#if HAVE_LIBCUDA==1
#include <cudaParticle.cuH>
#include <cudaMappingConstants.cuH>
#endif


/** Encapsulatates a SLGridSph (Sturm-Liouville basis) for use as
    force method.

    This provides the extra routines needed for n-body integration.
    See EmpCylSL in the exputils directory.
 */
class CylEXP : public EmpCylSL
{
public:

  CylEXP() : EmpCylSL() {}

  /** Constructor with parameters

      \par Parameters:

      @param multistep is the number of coefficient levels for n-body

      @param numr is the spherical radial order of the input basis

      @param lmax is the spherical angular order of the input basis

      @param mmax is the output aximuthal order for the EOF basis

      @param nord is the output radial order for the EOF basis

      @param ascale is the target disk scale LENGTH

      @param hscale is the target disk scale HEIGHT

      @param nodd is the number of vertically odd parity basis
      functions.  If unspecified, you get eigenvalue order.
      
   */
  CylEXP(int numr, int lmax, int mmax, int nord,
	 double ascale, double hscale, int Nodd=-1,
	 std::string cachename="") :
    EmpCylSL(numr, lmax, mmax, nord, ascale, hscale, Nodd, cachename)
  {
    differC1.resize(nthrds);
    differS1.resize(nthrds);
    for (auto & v : differC1) v.resize(multistep+1);
    for (auto & v : differS1) v.resize(multistep+1);

    unsigned sz = (multistep+1)*(MMAX+1)*NORDER;
    workC1.resize(sz);
    workC .resize(sz);
    workS1.resize(sz);
    workS .resize(sz);
  }

  /** Extrapolate and sum coefficents per multistep level to get
      a complete set of coefficients for force evaluation at an
      intermediate time step
  */
  void compute_multistep_coefficients(unsigned mlevel);

  /** Update the multi time step coefficient table when moving particle 
      <code>i</code> from level <code>cur</code> to level 
      <code>next</code>
  */
  //@{
  virtual void multistep_update_begin();
  virtual void multistep_update(int from, int to, double r, double z, double phi, double mass, int id);
  virtual void multistep_update_finish();
  virtual void multistep_reset();
  //@}

  //! Print debug info
  void multistep_debug();

};

#endif

// -*- C++ -*-

