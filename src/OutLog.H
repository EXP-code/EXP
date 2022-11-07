#ifndef _OutLog_H
#define _OutLog_H

#include <OutLog.H>

//! Log vital statistics every n steps
/**
   New version of log file puts all entries on a single line.  This makes
   for difficult human reading but great for machines.

   Centers of mass, Centers of velocity, Virial of Clausius and 2T/VC
   will be computed for each component
   
   NB: 
   the "Virial of Clausius" is a better indicatator of the true potential 
   energy for multiple components or external forces

   Current ordering (will be labeled in the file):
    - Global statistics
    	-# Current time
	-# Total mass
	-# Number of bodies
	-# COM (x, y, z)
	-# COV (x, y, z)
	-# Ang mom (x, y, z)
	-# Total KE
	-# Total PE
	-# Total E
	-# Virial of Clausius (VC)
	-# 2T/VC
	-# Time per step
	-# Total particles used
    - For each component in succession
    	-# Total mass
	-# Number of bodies
	-# COM (x, y, z)
	-# COV (x, y, z)
	-# Ang mom (x, y, z)
	-# Total KE
	-# Total PE
	-# Total E
	-# Virial of Clausius (VC)
	-# 2T/VC
	-# Particles used
*/
class OutLog : public Output
{

  typedef vector<double> dvector;

private:

  string filename;
  int nfreq;
  int laststep;

  double curwtime, lastwtime;
  double ektotxy;

  void initialize(void);
  bool firstime;

  vector<int>     nbodies, nbodies1;
  vector<int>     used, used1;
  vector<double>  mtot, mtot1;
  vector<dvector> com, com1;
  vector<dvector> cov, cov1;
  vector<dvector> angm, angm1;
  vector<dvector> ctr;
  vector<double>  ektot, eptot, eptotx, clausius;
  vector<double>  ektot1, eptot1, eptotx1, clausius1;

  vector<double>  com0, cov0, angmG, angm0, pos0, vel0;
  vector<double>  comG, covG, comL, covL, posL, velL;

  //! Valid keys for YAML configurations
  static const std::set<std::string> valid_keys;

public:

  //! Number of global entries
  static const int num_global = 19;
  //! Number of per-component entries
  static const int num_component = 20;

  //! Global column labels
  static char lab_global[][19];
  //! Per-component column labels
  static char lab_component[][20];

  /*! \brief Constructor
      \param line token list of specific parameters
  */
  OutLog(const YAML::Node& conf);

  /*! \brief Generate output
      \param nstep 
             the current step number (used e.g. to decide whether
	     or not to write output
      \param mstep is the current multistep level to decide whether or not to dump multisteps
      \param last 
      	     if true, output is generated independent of the step number
  */
  void Run(int nstep, int mstep, bool last);

};

#endif
