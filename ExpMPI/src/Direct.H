
#ifndef _Direct_H
#define _Direct_H

/**
   Computes the potential, acceleration and density using the direct summation

   Positions are passed from node to node on a ring.  Accelerations
   computed on requested particles array due to particles in component
   associated with the Direct force.  

   Settable parameters:
   @param soft_indx is the index in the particle double array for the 
   softening value.

   @param soft is a fixed softening value for all particles.  This can be 
   used even if a value is specified in the particle attribute array.

   @param ton is the midpoint of the adiabatic turn on for this component

   @param twid is the width for the adiabatic turn on.  Unless <code>ton</code>
   or <code>toff</code> is specified, the component begins and remains fully 
   on.
*/
class Direct : public PotAccel
{
private:

  int soft_indx;

  int to_proc, from_proc;
  int max_bodies;
  int ninteract;
  int ndim;

  double *tmp_buffer, *bod_buffer;

  double soft;
  bool fixed_soft;

  bool adiabatic;
  double ton, twid;

  void initialize();

  void determine_coefficients(void);
  void determine_acceleration_and_potential(void);

  void * determine_coefficients_thread(void * arg);
  void * determine_acceleration_and_potential_thread(void * arg);

public:

  //! The constructor
  /** <param> line passes in any explicit parameters */
  Direct(string& line);

  //! The constructor
  ~Direct();

  //! The main force call
  void get_acceleration_and_potential(vector<Particle>*);

};

#endif

