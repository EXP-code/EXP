#ifndef _OutHDF5_H
#define _OutHDF5_H

/** Write phase-space dumps at regular intervals from each node in
    component pieces.  This writer is indended to be Gadget-4 HDF5
    compliant.

    @param filename is the name of the output file
    @param nint is the number of steps between dumps
    @param nbeg is suffix of the first phase space %dump
    @param timer set to true turns on wall-clock timer for PS output
    @param gadget set to true for Gadget-4-type HDF5
    @param noids set to true turns off particle id writing
    @param checkpt sets checkpoint mode: write a real8 snapshot with a unique filename
    @param preserve prevents overwrite of snapshot files
    @param expconfig set to true stashes the YAML configuration file in the HDF5 file
    @param allmeta set to true writes all meta data to each HDF5 file
    @param directory set to true stores all HDF5 files in a named directory
*/
class OutHDF5 : public Output
{

private:

  std::string filename;
  bool real4=true, real8=false, ids=true, gadget4=false, chkpt=false;
  bool overwrite=false, expconfig=true, allmeta=false, directory=true, timer;
  unsigned write_flags;
  int nbeg;
  void initialize(void);
  std::vector<double> masses;
  std::vector<bool> multim;
  std::string hdf_dir, chkpt_dir;

  //! Valid keys for YAML configurations
  static const std::set<std::string> valid_keys;

  //! Check for single particle mass on the first invocation
  void checkParticleMasses();

  //! Gadget-4 format
  void RunGadget4(const std::string& path);
  
  //! PSP format
  void RunPSP(const std::string& path);
  

public:

  //! Constructor
  OutHDF5(const YAML::Node & conf);

  //! Provided by derived class to generate some output
  /*!
    \param nstep is the current time step used to decide whether or not
    to %dump
    \param mstep is the current multistep level to decide whether or not to dump multisteps
    \param last should be true on final step to force phase space %dump
    indepentently of whether or not the frequency criterion is met
    \param timer set to true turns on wall-clock timer for PS output
    \param noids set to true turns off particle id writing
    \param gadget set to true for Gadget-4-type HDF5
    \param checkpt set to true for checkpoint mode
    \param preserve set to true prevents overwrite of snapshot files
    \param expconfig set to true stashes the YAML configuration file in the HDF5 file
    \param allmeta set to true writes all meta data to each HDF5 file
    \param directory set to true stores all HDF5 files in a named directory
  */
  void Run(int nstep, int mstep, bool last);

};

#endif
