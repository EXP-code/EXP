#ifndef _CollideLTE_H
#define _CollideLTE_H

#include "HeatCool.H"
#include "Collide.H"

class CollideLTE : public Collide
{
private:
  vector<Precord> prec;
  HeatCool *hc;

protected:
  //! List of cross-sections per thread
  std::vector<sKey2Dmap> csections;

  //@{
  //! For debugging
  std::vector<unsigned> cellcnt;
  std::vector<double> minT, maxT, avgT, dispT;
  std::vector< std::vector<double> > tlist;

  unsigned numt;
  double tmin, tmax, dtmp;
  std::vector<double> thisto1, thisto2;

  unsigned numn;
  double nmin, nmax, ntmp;
  std::vector<double> nhisto1, nhisto2;

  static unsigned trhocnt;
  std::vector< std::vector<double> > trho;

  bool debug_enabled;
  //@}

protected:

  vector<double> deltaE;

  void initialize_cell(pHOT* tree, pCell* c, double vrelmax, int id);

  void initialize_cell_dsmc
  (pHOT* tree, pCell* c, sKey2Umap& nsel, double vrelmax, double tau, int id);

  void initialize_cell_epsm
  (pHOT* tree, pCell* c, sKey2Umap& nsel, double vrelmax, double tau, int id) 
  {
    initialize_cell_dsmc(tree, c, nsel, vrelmax, tau, id);  
  }

  void finalize_cell(pHOT* tree, pCell* c, double kedsp, int id);

  sKey2Dmap& totalCrossSections(int id)
  {    
    return csections[id];
  }

  sKey2Dmap& totalScatteringCrossSections(double crm, pCell* c, int id);

  bool hasHeatCool() { return true; }

  double crossSection(pHOT* tree, Particle* p1, Particle* p2, 
		      double crm, int id=0);


  int inelastic(pHOT *tree, Particle* p1, Particle* p2, double *crm, int id=0);

  double getCoolingRate(int id) { return coolheat[id]; }

  //@{
  //! Diagnostic output
  double totalSoFar, massSoFar;
  vector<double> lostSoFar;
  vector<double> coolheat;

  virtual void list_sizes();
  void list_sizes_proc(ostream*);
  //@}

public:
				// Cooling table
  static double Nmin;
  static double Nmax;
  static double Tmin;
  static double Tmax;
  static double TolV;
  static unsigned Nnum;
  static unsigned Tnum;
  static string cache;

  CollideLTE(ExternalForce *force, Component *comp,
	     double hsDiam, double diamfac, int nth=1);
  ~CollideLTE();

  virtual void set_timestep(int DTpos) { use_delt=DTpos; }

  virtual double Etotal();
  virtual double Mtotal();

  virtual void Elost(double* collide, double* epsm);

  virtual sKey2Umap generateSelection(pCell* c, sKeyDmap* Fn, double crm, double tau, int id,
				      double& meanLambda, double& meanCollP, 
				      double& totalNsel);

  void Debug(double t);
};

#endif
