#ifndef _UserResPoT_H
#define _UserResPoT_H

#include <massmodel.h>
#include <biorth.h>
#include <ResPot.H>

class UserResPot : public ExternalForce
{
private:
  
  string com_name;
  Component *c0;

  CVector bcoef, bcoefPP;
  AxiSymModel  *halo_model;
  AxiSymBiorth *halo_ortho;
  ResPot *respot;

  void * determine_acceleration_and_potential_thread(void * arg);
  void initialize();

  double ton, toff, delta, toffset, omega;
  
  int LMAX, NMAX, NUMR, L0, M0, L1, L2;
  double rmin, rmax, drfac;
  
  double MASS, LENGTH, COROT, A21, A32;

  int NUME, NUMK;

  string ctr_name;
  string model_file;

  void userinfo();

public:

  UserResPot(string &line);
  ~UserResPot();

};


#endif
