#ifndef _TreeDSMC_H
#define _TreeDSMC_H

#include <set>

#include <Particle.H>
#include <AxisymmetricBasis.H>
#include <ExternalForce.H>
#include <RunningTime.H>
#include <pHOT.H>

class Collide;
class CollideLTE;
class CollideIon;

/** DSMC particle routine

    <h2> Introduction </h2>

    The DSMC technique models a gas by following the trajectories of
    gas molecules and computes the effects of collisions between
    particles through probabilistically generated encounters. One of
    the fundamental assumptions of the DSMC method is that particle
    movement can be decoupled from collisional behavior.  Particles
    then move along paths determined by the mean force fields
    (gravitational and electromagnetic) with discrete changes in
    velocity and energy caused by collisions. In order to accomplish
    this, the simulation is broken up into time steps the magnitude of
    which are small compared to the average time between
    collisions. In this small interval of time, these two processes
    can realistically be treated as independent.

    The basic structure of the DSMC algorithm is as follows. In each
    computational cell, particles are chosen to collide in such a way
    as to reproduce the required macroscopic collision rate. The
    relative location of particles in the cell is not considered when
    selecting collision pairs. This fact restricts the size of a
    computational cell; it must be smaller than the local average mean
    free path of the flow field. Each of the colliding pairs can then
    chemically react or exchange internal energy. All of the colliding
    pairs then exchange energy elastically. Following the collision
    routines, each particle is moved throughout the domain for the
    duration of a single time step.  Particles are introduced at the
    inflow boundaries, removed at outflow boundaries, and reflected
    from surfaces. If a steady state has been reached in the
    simulation, as indicated by constancy of particle number and total
    energy in the domain, the macroscopic flow variables are sampled
    in each cell. At this point the algorithm is repeated until a
    sufficient number of samples has been collected.

    A complete DSMC simulation is a series of discrete time steps or
    iterations. During each iteration particles are moved and
    collisions are computed. The algorithm for a single iteration can
    be broken down into the following operations:
    <ol>
    <li> Move particles
    <li> Sort particles into cells
    <li> Select collision pairs 
    <li> Calculate collision mechanics 
    <li> Sample macroscopic properties 
    </ol>

    <h2> Structure </h2>

    Step (1) from the list above is performed by EXP internals,
    typical of an n-body code.  Steps (2) and (3) are initialed by the
    determine_acceleration_and_potential() member function at each
    step.  Specifically, this member calls the appropriate members of
    pHOT to repartition the particles into cells and rebuild the tree,
    if necessary.  determine_acceleration_and_potential() also
    accumulates diagnostics and time step information for future
    repartitioning and time-step level adjustment.

    Then, determine_acceleration_and_potential() requests an estimate
    for the cross sections for the particles in the cell from the
    Collide instance  Collide is an abstract base class that provides
    all of the information about the particle interactions.  Multiple
    collision schemes may be derived form Collide and selected using
    the parameter <code>ctype</code>.

    Step (4).  The cross sections are used to estimate the collision
    rates.  If the collision rates are higher the the EPSM limit, and
    EPSM is enabled (see the parameter <code>epsm</code>), the Collide
    instance will be asked for a mean heating/cooling rate and and an
    EPSM step will be performed.  Otherwise, the DSMC algorithm is
    invoked to propose collisions among each cells particles at the
    computed rate.

    Finally, Step(5), final per-cell diagnostics are performed and
    accumulated.  These are reported in various log files and used for
    adjusting time-step levels for each cell.

    <h2> See Collide for deriving a new collsion routine </h2> 

    See Collide for deriving a new collsion routine.  Then, add your
    constructor to the list following the comment <code>Create the
    collision instance from the allowed list</code>.

    <h2> Parameters </h2>

    @param Lunit is the physical scale in cm per system units (default: 300 kpc)
    @param Tunit is the physical time scale in years per system unit (default: gravitational)
    @param Munit is the physical mass in units of mp per system units (default: 1e12 Msun)
    @param boxsize is the size of a side in system units (default: 1, 1, 1)
    @param boffset is offset in the side for the origin  (default: 0, 0, 0)
    @param ncell is the desired number of particles per cell (default: 7)
    @param Ncell is the desired number of particles per sample cell (default: 64)
    @param cnum non-zero adjusts cross section to yield cnum collisions per cell
    @param epsm >0 enables EPSM computation for cells with  MFP/size ratio < epsm
    @param hsdiam is a multiplicative factor for the Bohr radius for computing the hard-sphere cross-section (converted to system units) 
    @param crossfac is a multiplicative factor for geometric scaling of the scattering cross section
    @param spc is a weight for the collision number computation describing the converstion of mass in protons to number
    @param frontier enables printing of active cell list (for debugging)
    @param tsdiag enables printing gas time step diagnostics (for debugging)
    @param voldiag enables printing volume-cell diagnostics (for debugging)
    @param mfpstat enables printing of mean free path and time step statistics (for debugging)
    @param dryrun enabled turns off collsions (primarily for debugging)
    @param nsteps is the number of steps between diagnostic output (default: -1, off)
    @param msteps is the maximum level for diagnostic output (default: -1, off)
    @param use_temp is the float particle attribute position for temperature
    @param use_dens is the float particle attribute position for density
    @param use_delt include energy loss rate in setting time step (if multistep is enabled) 
    @param use_multi attempts to set the time step (if multistep is enabled)
    @param use_pullin set to true uses the Pullin algorithm for selecting velocities (default: false)
    @param use_repair repairs the tree rather than remakes the tree
    @param use_exes track energy loss per particle for exact energy conservation
    @param use_Kn track Knudsen number in float attribute position n (default: -1)
    @param use_St track Strouhal number in float attribute position n (default: -1)
    @param use_vol track cell volume number in float attribute position n (default: -1)
    @param esol solves the LTE cooling rate equation for each cell
    @param coolfrac is the fraction of the cooling time for setting the time step (default: 0.1)
    @param cba use the Collisional Boltzmann Approximation from Alexander et al.
    @param madj levels and below perform a full domain decomposition and above madj the tree is repaired (default: 0)
    @param wght sets the cell time partitioning (default: 1 = on)
    @param sub_sample toggles the key partition algorithm that uses only a subsample of the key list to decompose the domaine (default: true)
    @param effort sets the use of DSMC effort to partition particles (default: 1=true)
    @param treechk turns on particle and cell checking in the tree (may be very time consuming) (default: 0 = off)
    @param mpichk turns on barriers for timing and sync checking (default: 0 = off)
    @param remap sets the processor remapping frequency in full steps (default: 0)
    @param enhance multiplies the standard cooling/heating rate for testing (default: 1)
    @param compname is the name of the component that is evolved using DSCM
    @param tspow is the power of two spacing between histogram bins in the cooling diagnostic (default: 4)
    @param nocool set to true turns off cooling for testing purposes
    @param use_effort set to true (default) uses the effort expended in computing the force per particle from DSMC to load balance the simulation
    @param mfpts set to true (default) uses the mean free path to estimate a time step
    @param ntc set to true (default) uses the no time counter algorithm to select collistions
    @param ntcfac is an upscaling factor for the NTC database in mixed units (default: 2.0)
    @param ctype is a string referring to the name of the desired collision type (default: LTE)
    @param rrtype is a string referring to the name of the desired radiative recombination cross section type.  This is only used by the ctype=Ion method (default: Mewe)
    @param equiptn set to true (default: true) enforces equipartition between electron and ion pairs. This is only used by the ctype=Ion method.
    @param seed is the default random number seed.  If seed=0, /dev/urandom is used.  Otherwise, Collide class default is used.

    <h2>Note on multiple species</h2>

    There are two versions of the multispecies interactions.  Each method has two different variants.
    <ol>
    <li> Each particle acts as an individual ion (i.e. atom in a
    particular ionization state which may be neutral.  This version,
    the default, is chosen by specifying the <code>keypos</code>
    attribute for the Component.  The variants are:
    <ol> 
    <li> The <code>Direct</code> method.  This is closest to pure
    DSMC.  Each particle represents the same number of atoms or ions.
    <li> The <code>Weight</code> method.  All particle are assumed to
    have the same mass but represent a different number of true
    particles according their physical chemical abundance.
    </ol>
    <li> Each particle represents a number of elements in various
    charge states.  Use of multiple species is initiated by setting
    the <code>species</code> flag to an integer \f$\ge0\f$.  The
    argument indicates the position of the species index in the
    per-particle integer parameter vector.  The species indices to not
    have to be sequential, only unique.  They are stored in the code
    using a STL <code>map</code> template.  This scheme as two variants:
    <ol>

    <li> In the <code>Hybrid</code> variant, each particle has the same
    mass but represents a different number of atoms or ions, as in the
    <code>Weight</code> above.  However, the particle is divided into
    pieces represented by the ionic level and the ionic level
    fractions are represented using the <code>species</code> map as
    described above.

    <li> In the <code>Trace</code> variant, each particle has the same
    mass divided into elements according to a predetemined chemical
    abundance <i>and</i> ionic substates.
    </ol> 
    </ol>
*/
class TreeDSMC : public ExternalForce
{
  static const std::string version;

private:
  
  int         ncell, Ncell;
  int         cnum, wght;
  int         nsteps, msteps, remap;
  unsigned    madj, seed;
  int         use_key, use_temp, use_dens, use_delt, use_exes, use_Kn;
  int         use_St, use_vol, use_Eint, tspow;
  double      hsdiam;
  double      crossfac;
  double      collfrac;
  double      boxsize;
  double      boffset;
  double      epsm;
  double      coolfrac;
  double      enhance;
  double      ntcfac;
  std::string ctype;
  std::string rrtype;

  bool frontier, tsdiag, voldiag, mfpstat, use_multi, use_pullin, use_repair;
  bool cba, cbadiag, dryrun, nocool, ntc, esol;
  bool sub_sample, treechk, mpichk, mfpts, hybrid, dumpHOT;

				// Box size and offsets
  std::vector<double> bSiz, bOff;

  unsigned cntHOT;		// For frontier dump files

				// For NTC db trimming
  unsigned ageout;

				// For MFP and time step statistics
  std::vector<double> quant, mfp_, ts_, coll_, nsel_, cool_, rate_;

  unsigned stepnum;
  double   curtime;

  Collide *collide;
  double ElostTotCollide, ElostTotEPSM;

  std::string comp_name, spec_map;
  Component *c0;

  void determine_acceleration_and_potential(void);
  void * determine_acceleration_and_potential_thread(void * arg) 
  { return (NULL); }
  void assignTempDensVol();
  void initialize();


  void outHeader0(ostream& out)
  {
    out << left << scientific << "   "
	<< setw(12) << "" << " " << setw(7) << "  min" << "  ";
    for (unsigned n=0; n<pHOT::ntile; n++) {
      ostringstream sout;
      sout << "  " << setw(2) << pHOT::qtile[n] << "%";
      out << setw(7) << sout.str() << "  ";
    }
    out << setw(7) << "  max" << "   " << right;
    out << setw(4) << "min" << "  ";
    for (unsigned n=0; n<pHOT::ntile; n++) {
      ostringstream sout;
      sout << setw(2) << pHOT::qtile[n] << "%";
      out << setw(4) << sout.str() << "  ";
    }
    out << setw(4) << "max" << endl << left;
  }
    
  void outHeader1(ostream& out)
  {
    out << "                    " 
	<< "    min     ";
    for (unsigned n=0; n<pHOT::ntile; n++)
      out << "    " << setw(2) << pHOT::qtile[n] << "%     ";
    out << "    max     " << endl;
  }


  void outHelper0(ostream& out, const char* lab, int c,
		  vector< vector<double> >& v, vector<double>& t)
  {
    unsigned nt = t.size();
    out << left << scientific << "   "
	<< setw(12) << lab << "[" << setprecision(1);

    for (unsigned k=0; k<nt-1; k++) 
      out << setw(7) << v[k][c] << ", ";

    out << setw(7) << v[nt-1][c] << "] [" << right << fixed;

    for (unsigned k=0; k<nt-1; k++) 
      out << setw(4) << v[k][c]/(t[k]+1.0e-10)*100.0 << ", ";

    out << setw(4) << v[nt-1][c]/(t[nt-1]+1.0e-10)*100.0 << "]%" 
	<< endl << left;
  }


  template<typename T>
  void outHelper1(ostream& out, const char* lab, vector<T>& v)
  {
    unsigned nt = v.size();
    out << left << setprecision(6) << right
	<< "      *** " << lab << ": [";
    for (unsigned k=0; k<nt-1; k++) 
      out << setw(10) << v[k] << ", ";
    out << setw(10) << v[nt-1] << "]" << endl << left;
  }

				// The tree
  double volume;

				// Timers
  Timer partnTime, tree1Time, tradjTime, tcellTime, tstepTime, llistTime;
  Timer clldeTime, clldeWait, partnWait, tree1Wait, tree2Wait, timerDiag;
  RunningTime overhead, fullstep;
  double tM0;

				// Elapsed time counters
  double partnSoFar, tree1SoFar, tradjSoFar, tcellSoFar, tstepSoFar;
  double waitcSoFar, waitpSoFar, wait1SoFar, wait2SoFar, timerSoFar;
  double collideSoFar;
  
				// Full time step counters
  vector<double> tt;

  void userinfo();


				// Debugging
  void triggered_cell_body_dump(double time, double radios);
  void TempHisto();


				// Make species map
  void makeSpeciesMap();

#ifdef TIMER
  Timer *timer;
#endif

protected:

  //@{
  //! Species lists
  std::map<speciesKey,  unsigned long> spec, spec1;
  std::map<speciesKey,  double> collFrac;
  std::set<speciesKey>  spec_list;
  double                TempTot;
  //@}

  //! Valid collisions types
  static std::set<std::string> colltypes;

  //! Initialize collision types
  static void initialize_colltypes()
  {
    if (colltypes.size()) return;
    colltypes.insert("LTE");
    colltypes.insert("Ion");
  }

  //! Check for valid collision types (return true if good)
  bool check_ctype(std::string type)
  {
    if (colltypes.find(type) == colltypes.end())  return false;
    else                                          return true;
  }


public:

				// System to physical scaling factors
  static double Lunit;
  static double Tunit;
  static double Vunit;
  static double Munit;
  static double Eunit;

  static std::vector<double> atomic_weights;

				// Use timing info to compute effort
  static bool use_effort;

  //! Constructor
  TreeDSMC(const YAML::Node& conf);

  //! Destructor
  ~TreeDSMC();

  //! Finish
  void finish();
};

#endif
