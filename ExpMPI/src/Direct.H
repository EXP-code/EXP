
#ifndef _Direct_H
#define _Direct_H

//! computes the potential, acceleration and density using the direct summation
class Direct : public PotAccel
{
private:

  int soft_indx;

  void initialize();

  void determine_coefficients(void);
  void determine_acceleration_and_potential(void);

  void * determine_coefficients_thread(void * arg);
  void * determine_acceleration_and_potential_thread(void * arg);

public:

  //! The constructor
  /** <param> line passes in any explicit parameters */
  Direct(string& line);

  //! The main force call
  void get_acceleration_and_potential(vector<Particle>*);

};

#endif

