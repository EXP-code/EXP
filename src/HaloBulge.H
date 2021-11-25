#ifndef _HaloBulge_H
#define _HaloBulge_H

#include <ExternalForce.H>

//! Compute halo + bulge external force
class HaloBulge : public ExternalForce
{
private:
  
  void initialize();
  
  //! Halo model type (default: file)
  int HMODEL;
  //! Halo model file (default: w05)
  string INFILE;
  //! Halo mass scale factor (default: 1.0)
  double MHALO;
  //! Halo radius scale factor (default: 1.0)
  double RHALO;
  //! Halo minimum model radius (default: 1.0e-3)
  double RMODMIN;
  //! Halo maximum model radius (default: 20.0)
  double RMOD;
  //! Bulge core/scale size (default: 1.0)
  double RBCORE;
  //! Bulge mass scale factor (default: 1.0)
  double MBULGE;
  //! Bulge radius scale factor (default: 1.0)
  double RBULGE;
  //! Bulge miniimum model radius (default: 1.0e-3)
  double RBMODMIN;
  //! Bulge maximum model radius (default: 20.0)
  double RBMOD;
  
  AxiSymModPtr model, bmodel;

  void * determine_acceleration_and_potential_thread(void * arg);

public:

  //! Constructor
  HaloBulge(const YAML::Node& conf);

};


#endif
