
#ifndef _NoForce_H
#define _NoForce_H

#include <string>

#include <PotAccel.H>

/**
   A *NULL* force method

   Applies no acceleration to the particles whatsoever.  This allows the
   user to supply an alternative as an external force.
*/
class NoForce : public PotAccel
{
private:

  void initialize() {}

  void determine_coefficients(void) {}
  void determine_acceleration_and_potential(void) {}

  void * determine_coefficients_thread(void * arg) {}
  void * determine_acceleration_and_potential_thread(void * arg) {}

public:

  //! The constructor
  /*! \param line is ignored here */
  NoForce(std::string& line);

  //! The constructor
  ~NoForce();

  //! The main force call
  void get_acceleration_and_potential(Component*) {}

};

#endif

