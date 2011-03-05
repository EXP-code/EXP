#ifndef _UserTreeDSMC_H
#define _UserTreeDSMC_H

#include <Particle.H>
#include <AxisymmetricBasis.H>
#include <ExternalForce.H>
#include <RunningTime.H>
#include <pHOT.H>

class CollideLTE;

/** DSMC particle routine
    @param Lunit is the physical scale in cm per system units (default: 300 kpc)
    @param Tunit is the physical time scale in years per system unit (default: gravitational)
    @param Munit is the physical mass in units of mp per system units (default: 1e12 Msun)
    @param boxsize half the size of a side in system units (default: 1)
    @param boxratio is the ratio of the z box height to the x or y box height (default: 1)
    @param jitter is the fraction of the box size added to the tree origin (default: 0)
    @param ncell is the desired number of particles per cell (default: 7)
    @param Ncell is the desired number of particles per sample cell (default: 64)
    @param cnum non-zero adjusts cross section for cnum collisions per cell
    @param epsm >0 enables EPSM computation for cells with  MFP/size ratio < epsm
    @param diamfac is a multiplicative factor for the Bohr radius (converted to system units
    @param collfrac is a multiplicative fraction for collsion number
    @param frontier enables printing of active cell list (for debugging)
    @param tsdiag enables printing gas time step diagnostics (for debugging)
    @param voldiag enables printing volume-cell diagnostics (for debugging)
    @param mfpstat enables printing of mean free path and time step statistics (for debugging)
    @param cbadiag enables printing of length-scale statistics for the CBA algorithm (for debugging)
    @param dryrun enabled turns off collsions (primarily for debugging)
    @param nsteps is the number of steps between diagnostic output (default: -1)
    @param msteps is the maximum level for diagnostic output (default: -1)
    @param use_temp is the particle attribute position for temperature
    @param use_dens is the particle attribute position for density
    @param use_multi attempts to set the time step (if multistep is enabled)
    @param use_repair repairs the tree rather than remakes the tree
    @param use_delt include energy loss rate in setting time step (if multistep is enabled) 
    @param use_exes track energy loss per particle for exact energy conservation
    @param use_Kn track Knudsen number in attribute position n (default: -1)
    @param use_St track Strouhal number in attribute position n (default: -1)
    @param use_vol track cell volume number in attribute position n (default: -1)
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
*/
class UserTreeDSMC : public ExternalForce
{
private:
  
  int ncell, Ncell;
  int cnum, wght;
  int nsteps, msteps, remap;
  unsigned madj;
  int use_temp, use_dens, use_delt, use_exes, use_Kn, use_St, use_vol, tspow;
  double diamfac;
  double diam;
  double collfrac;
  double boxsize;
  double boxratio;
  double jitter;
  double epsm;
  double coolfrac;
  double enhance;

  bool frontier, tsdiag, voldiag, mfpstat, use_multi, use_pullin, use_repair;
  bool cba, cbadiag, dryrun, nocool, ntc, tube, slab, esol;
  bool sub_sample, treechk, mpichk;

  string sdim;

				// For MFP and time step statistics
  vector<double> quant, mfp_, ts_, coll_, nsel_, cool_, rate_;

  unsigned stepnum;
  double curtime;

  CollideLTE *collide;
  double ElostTotCollide, ElostTotEPSM;

  string comp_name;
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
  RunningTime overhead;

  void userinfo();


				// Debugging
  void triggered_cell_body_dump(double time, double radios);
  void TempHisto();

#ifdef TIMER
  Timer *timer;
#endif

public:

				// System to physical scaling factors
  static double Lunit;
  static double Tunit;
  static double Vunit;
  static double Munit;
  static double Eunit;
				// Use timing info to compute effort
  static bool use_effort;

  //! Constructor
  UserTreeDSMC(string &line);

  //! Destructor
  ~UserTreeDSMC();

};

#endif
