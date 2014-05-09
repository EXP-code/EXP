#ifndef _CollideIon_H
#define _CollideIon_H

#include "Collide.H"
#include "interactSelect.H"
#include "Elastic.H"
#include "Ion.H"

class CollideIon : public Collide
{
private:
  //! Global variables for the calculations with the individual ions
  std::map <int, std::map<int, Ion> > IonList;

  //! Not sure what this does . . . looks like cruft.  Remove??
  std::vector<int> Ni;

  chdata ch;
  Elastic elastic;
  InteractSelect IS;
  unsigned int NUM;

  std::map<unsigned short, double> atomic_weights;
  std::map<unsigned short, std::pair<unsigned short, double> > interInfo;

  std::pair<double, double> ff_d;
  std::pair<double, double> CE_d;
  std::pair<double, double> CI_d;
  std::pair<double, double> RR_d;
  std::pair<double, double> dv;

  double eV_av;
  double eV_N;
  double eV_max;
  double eV_min;
  double eV_10;
  
  ACG *gen;
  Uniform *unit;
  
  //! Debugging file names
  std::string coll_file_debug, energy_file_debug;

  //! Initialize debugging file names
  void printCollInitialize();

protected:

  //! List of cross-ections per thread
  std::vector<sKey2Dmap> csections;
  
protected:
  
  double deltaE;
  
  
  void initialize_cell_dsmc (pHOT* tree, pCell* cell, 
			     double vrelmax, double dt, 
			     map<speciesKey, map<speciesKey, unsigned> >& nsel, int id)
  {
    initialize_cell(tree, cell, vrelmax, dt, nsel, id);
  }
  
  void initialize_cell_epsm (pHOT* tree, pCell* cell, 
			     double vrelmax, double dt, 
			     map<speciesKey, map<speciesKey, unsigned> >& nsel, int id)
  {
    initialize_cell(tree, cell, vrelmax, dt, nsel, id);
  }
  
  
  void initialize_cell(pHOT* tree, pCell* cell,
		       double rvmax, double tau,
		       map<speciesKey, map<speciesKey, unsigned> >& nsel,
		       int id);
  
  void finalize_cell(pHOT* tree, pCell* cell, double kedsp, int id);

  
  map<speciesKey, map<speciesKey, double> >& 
  totalCrossSections(double crm, pCell* c, int id);
  
  
  int inelastic(pHOT *tree, Particle* p1, Particle* p2, double *crm, int id=0);
  bool hasHeatCool() { return true; }
  double getCoolingRate(int id) { return 0; }
  
  
  //@{
  //! Diagnostic output
  double totalSoFar, massSoFar;
  vector<double> lostSoFar;
  
  void list_sizes_proc(ostream*);
  //@}
  
public:
  static double Nmin;
  static double Nmax;
  static double Tmin;
  static double Tmax;
  static double TolV;
  static unsigned Nnum;
  static unsigned Tnum;
  static string cache;
  static bool frost_warning;

  //! Constructor
  CollideIon(ExternalForce *force, double hsDiam, double diamfac, int nth=1);

  //! Destructor
  ~CollideIon();
  
  virtual void set_timestep(int DTpos) { use_delt=DTpos; }
  
  virtual double Etotal();
  virtual double Mtotal();

  void printCollSummary();
  void printCollGather();
  void resetColls();

  virtual void *timestep_thread(void *arg);
  
  virtual void Elost(double* collide, double* epsm);
  
  void Debug(double t);
};

#endif
