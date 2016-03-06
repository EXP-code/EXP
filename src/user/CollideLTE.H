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
  std::vector<sKey2Amap> csections;

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

  virtual void 
  initialize_cell(pHOT* const tree, pCell* const c, double vrelmax, int id);

  virtual void 
  initialize_cell_dsmc
  (pHOT* const tree, pCell* const c, sKey2Amap& nsel, 
   double vrelmax, double tau, int id);

  virtual void
  initialize_cell_epsm
  (pHOT* const tree, pCell* const c, sKey2Amap& nsel, 
   double vrelmax, double tau, int id) 
  {
    initialize_cell_dsmc(tree, c, nsel, vrelmax, tau, id);  
  }

  virtual void 
  finalize_cell(pHOT* tree, pCell* c, sKeyDmap* const Fn, 
		double kedsp, double tau, int id);

  virtual sKey2Amap& totalCrossSections(int id)
  {    
    return csections[id];
  }

  sKey2Dmap& totalScatteringCrossSections(double crm, pCell* c, int id);

  bool hasHeatCool() { return true; }

  virtual double crossSection(int id, pCell* const tree, 
			      Particle* const p1, Particle* const p2, double crm, 
			      const Interact::T& ityp=Interact::singleton);

  virtual int inelastic(int id, pCell* const cell, 
			Particle* const p1, Particle* const p2, double *crm, 
			const Interact::T& ityp=Interact::singleton,
			double weight=-1);

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
	     double hsDiam, double crossfac, int nth=1);
  ~CollideLTE();

  virtual void set_timestep(int DTpos) { use_delt=DTpos; }

  virtual double Etotal();
  virtual double Mtotal();

  virtual void Elost(double* collide, double* epsm);

  virtual sKey2Amap generateSelection(pCell* c, sKeyDmap* Fn, double crm, double tau, int id,
				      double& meanLambda, double& meanCollP, 
				      double& totalNsel);

  void Debug(double t);
};

#endif
