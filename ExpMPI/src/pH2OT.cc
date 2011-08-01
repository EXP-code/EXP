#include <cstdlib>
#include <cstdio>
#include <ctime>

#include <iostream>
#include <iomanip>
#include <fstream>
#include <sstream>
#include <vector>
#include <list>
#include <map>
#include <algorithm>
#include <cmath>

#ifdef USE_GPTL
#include <gptl.h>
#endif

#include <execinfo.h>

using namespace std;

#include "global.H"
#include "pH2OT.H"

double   pH2OT::sides[]  = {2.0, 2.0, 2.0};
double   pH2OT::box[]    = {2.0, 2.0, 2.0};
double   pH2OT::offset[] = {1.0, 1.0, 1.0};
unsigned pH2OT::qtile[]  = {10,  50,  90 };

int tCell::live = 0;		// Track number of instances
int tTree::live = 0;		// Track number of instances
int pTree::live = 0;		// Track number of instances

pTree::pTree(pH2OT *p, tCell *q, unsigned key)
{
  live++;

  ph      = p;
  pc      = q;
  tkey    = key;
  root    = new pCell(this);
  C       = q->getCorner();
  clevels = vector< set<pCell*> >(multistep+1);
}

pTree::~pTree()
{
  live--;
  delete root;
}

//
// Test internal paritioning lists for consistency
// [set to false for production]
//
bool pH2OT::list_check = false;

//
// Turn on/off subsampling the key list for partitioning
//
bool pH2OT::sub_sample = true;

// For formatting
unsigned pH2OT::klen = 3*nbits/4+6;


// For diagnostic MPI communication
MPI_Datatype pH2OT::CellDiagType;


template<class U, class V>
struct pair_compare
{
  bool operator()(const pair<U, V>& a, const pair<U, V>& b)
  { return a.first<=b.first; }
};

void pH2OT::bomb(const string& membername, const string& msg)
{
  cerr << "pH2OT::" << membername << "(): " << msg << endl;
  MPI_Abort(MPI_COMM_WORLD, 497);
}

//
//  Constructor: initialize domain
//
pH2OT::pH2OT(Component *C)
{
  cc = C;			// Register the calling component

				// Sanity checks

  if (pkbits*tCell::ndim >= sizeof(unsigned)*8) {
    unsigned maxb = sizeof(unsigned)*8/tCell::ndim;
    if (maxb*tCell::ndim >= sizeof(unsigned)*8) maxb--;
    ostringstream mesg;
    mesg << "pkbits=" << pkbits << " but must be less than " << maxb;
    bomb("pH2OT::pH2OT", mesg.str());
  }

  if (nbits*3 >= sizeof(key_type)*8) {
    unsigned maxb = sizeof(key_type)*8/3;
    if (maxb*3 >= sizeof(key_type)*8) maxb--;
    ostringstream mesg;
    mesg << "nbits=" << nbits << " but must be less than " << maxb;
    bomb("pH2OT::pH2OT", mesg.str());
  }
  klen = 3*nbits/4+6;

				// pkfactor is the fractional length
				// scale for each HOT side
  pkfactor = 1.0/static_cast<double>(1u<<pkbits);
  for (unsigned k=0; k<tCell::ndim; k++) {
    unsigned j = tCell::idim[k];
    box[j] = sides[j]*pkfactor;
  }
  for (unsigned j=0; j<3; j++) {
    bool unclaimed = true;
    for (unsigned k=0; k<tCell::ndim; k++) 
      if (tCell::idim[k]==j) unclaimed = false;
    if (unclaimed) box[j] = sides[j];
  }

  pkmask = 0u;			// Mask per dimension
  for (unsigned k=0; k<pkbits; k++) {
    pkmask <<= 1;
    pkmask |= 1u;
  }

  tkmask = 0u;			// Total mask
  for (unsigned k=0; k<tCell::ndim; k++) {
    tkmask <<= pkbits;
    tkmask |= pkmask;
  }

				// Total # of regions in primary partition
  tottrees = (1u << pkbits*tCell::ndim);

				// Total volume
  total_volume = sides[0]*sides[1]*sides[2];

				// Total volume of oct-tree region
  volume = box[0]*box[1]*box[2];

  key_min = key_type(1u) << (nbits*3  );
  key_max = key_type(1u) << (nbits*3+1);

  // Compute tree structure for pTrees
  trees.reset(this, pkbits);

  numkeys = 0;
				// Timers for diagnostic output
  timer_cstatus.Microseconds();
  timer_xchange.Microseconds();
  timer_convert.Microseconds();
  timer_tadjust.Microseconds();
  timer_prepare.Microseconds();
  timer_cupdate.Microseconds();
  timer_scatter.Microseconds();
  timer_repartn.Microseconds();
  timer_keybods.Microseconds();
  timer_schecks.Microseconds();
  timer_waiton0.Microseconds();
  timer_waiton1.Microseconds();
  timer_waiton2.Microseconds();
  timer_bodlist.Microseconds();
  timer_celladj.Microseconds();
  timer_getsta1.Microseconds();
  timer_getsta2.Microseconds();
  timer_getsta3.Microseconds();

  use_weight = true;
 
  // Initialize timing structures
  //
  cstat3 = exchg3 = prepr3 = cnvrt3 = tadjt3 = updat3 = scatr3 
    = reprt3 = keybd3 = schks3 = wait03 = wait13 = wait23 
    = bodls3 = cladj3 = cells1 = cells2 = cells3 = barri3 
    = vector<float>(numprocs);

  numk3    = vector<unsigned>(numprocs);
  numfront = vector<int>(numprocs);
  displace = vector<int>(numprocs);

  //
  // Generate the initial particle partition; this inializes the map
  // for the getProc() function
  //
  use_number = true;
  Remap();


  // Change "true" to "false" to disable the barrier checking
  //
  barrier = new BarrierWrapper(MPI_COMM_WORLD, true);

  // Make the tree diagnostic MPI structure
  //
  const int nf = 9;
  MPI_Aint disp[nf];
  int blocklen [nf] = {1, 1, 1, 1, 1, 1, 1, 1, 1};
  MPI_Datatype type[nf] = {MPI_UNSIGNED, MPI_UNSIGNED, MPI_UNSIGNED, 
			   MPI_UNSIGNED, MPI_UNSIGNED,
			   MPI_DOUBLE, MPI_DOUBLE, MPI_DOUBLE, MPI_DOUBLE};
  CellDiag buf;
  MPI_Get_address(&buf.ntre,	&disp[0]);
  MPI_Get_address(&buf.bcel,	&disp[1]);
  MPI_Get_address(&buf.btot,	&disp[2]);
  MPI_Get_address(&buf.mind,	&disp[3]);
  MPI_Get_address(&buf.maxd,	&disp[4]);
  MPI_Get_address(&buf.navg,	&disp[5]);
  MPI_Get_address(&buf.nvar,	&disp[6]);
  MPI_Get_address(&buf.bavg,	&disp[7]);
  MPI_Get_address(&buf.bvar,	&disp[8]);

  for (int i=nf-1; i>=0; i--) disp[i] -= disp[0];
  
  MPI_Type_create_struct(nf, blocklen, disp, type, &CellDiagType);
  MPI_Type_commit(&CellDiagType);
}


pH2OT::~pH2OT()
{
  map<unsigned, tCell*>::iterator it;
  for (it=trees.frontier.begin(); it!=trees.frontier.end(); it++)
    delete it->second->ptree;
  delete barrier;
}

key_type pH2OT::getKey(double *p, vector<double>& c)
{
  const double tol = 1.0e-12;

  timer_keybods.start();

#ifdef USE_GPTL
  GPTLstart("pH2OT::getKey");
#endif

  // Bad values
  //
  for (unsigned k=0; k<3; k++) {
    if (isnan(p[k]) || isinf(p[k])) {
      d_val++;
      timer_keybods.stop();
      return key_type(0u);
    }
  }

  // Out of bounds?
  //
  double z[3];
  for (unsigned k=0; k<3; k++) { 
    z[k] = (p[k] - c[k])/box[k];
				// Deal with some roundoff/truncation error
    if (z[k] < 0.0 && z[k] > -tol)  
      z[k] = 0.0;
    else if (z[k] >= 1.0 && z[k] < 1.0+tol)
      z[k] = 1.0 - tol;
    else if (z[k] < 0.0 || z[k] >=1.0) {
#ifdef USE_GPTL
      GPTLstop("pH2OT::getKey");
#endif
      timer_keybods.stop();
      return key_type(0u);
    }
  }

#ifdef I128
  const double factor = (key_type(1u)<<nbits).toDouble();
#else
  const double factor = key_type(1u)<<nbits;
#endif
  const key_type mask = 0x1u;

  vector<key_type> bins(3, 0u);

				// Reverse the order
  for (unsigned k=0; k<3; k++)
    bins[2-k] = key_type( floor(z[k]*factor) );
  
  key_type place = 1u;
  key_type  _key = 0u;

  for (unsigned i=0; i<nbits; i++) {
    for (unsigned k=0; k<3; k++) {
      _key |= (bins[k] & mask)*place;
      place = place << 1;
      bins[k] = bins[k] >> 1;
    }
  }

  _key += place;		// Leading placeholder for cell masking

#ifdef USE_GPTL
  GPTLstop("pH2OT::getKey");
#endif

  timer_keybods.stop();
  return _key;
}

string pH2OT::printKey(key_type p)
{
  ostringstream sout, sret;

  unsigned short cnt = 0;
  for (unsigned k=0; k<nbits; k++) {
#ifdef I128
    sout << ( (p & 1u).toUint() ? '1' : '0' );
#else
    sout << ( (p & 1u) ? '1' : '0' );
#endif
    if (++cnt==3) {sout << '.'; cnt = 0;}
    p = p>>1;
  }

  string s = sout.str();	// Reverse the string
  for (unsigned k=0; k<s.size(); k++) sret << s[s.size()-1-k];

  return sret.str();
}


void pH2OT::Remap()
{
  if (use_number)
    Remap_by_number();
  else
    Remap_by_effort();
}

void pH2OT::Remap_by_number()
{
  //
  // Make a map containing the keys from the full frontier
  //
  map<unsigned, unsigned> counts;
  //      ^         ^
  // key--/         |
  //                |
  // number---------/
  //

  //
  // Make a key-counts map from the particle list
  //
  {
    PartMapItr t;
    for (t=cc->Particles().begin(); t!=cc->Particles().end(); t++) {
      t->second.tree = tCell::getKey(&(t->second.pos[0]));
      if (t->second.tree) {
	if (counts.find(t->second.tree) == counts.end())
	  counts[t->second.tree]  = 1;
	else
	  counts[t->second.tree] += 1;
      }
    }
  }

  //
  // Sum up the counts for each key on the frontier
  //
  unsigned csize1 = counts.size(), csize0 = 0;
  vector<unsigned> csize(numprocs);
  MPI_Allgather(&csize1, 1, MPI_UNSIGNED, &csize[0], 1, MPI_UNSIGNED,
		MPI_COMM_WORLD);

  for (int n=0; n<numprocs; n++) csize0 = max<unsigned>(csize0, csize[n]);


  //
  // Order the key list by cell counts
  //
  multimap<unsigned, unsigned> m;
  //          ^         ^
  // counts---/         |
  //                    |
  // key----------------/
  //

  vector<unsigned> one(csize0), two(csize0);
  for (int n=0; n<numprocs; n++) {
    if (n==myid) {
      map<unsigned, unsigned>::iterator j=counts.begin();
      for (unsigned i=0; i<csize[n]; i++) {
	one[i] = j->first;
	two[i] = j->second;
	j++;
      }
    }
  
    MPI_Bcast(&one[0], csize[n], MPI_UNSIGNED, n, MPI_COMM_WORLD);
    MPI_Bcast(&two[0], csize[n], MPI_UNSIGNED, n, MPI_COMM_WORLD);

    for (unsigned i=0; i<csize[n]; i++)
      m.insert(pair<unsigned, unsigned>(two[i], one[i]));
  }

  //
  // Assign processors round-robin in cell count order
  //
  int cnt = 0;
  keyproc.clear();
  for (multimap<unsigned, unsigned>::iterator k=m.begin(); k!=m.end(); k++)
    keyproc[k->second] = (cnt++ % numprocs);

  //
  // Diagnostic
  //
  if (myid==0) {
    ostringstream sout;
    sout << outdir << runtag << ".counts." << cc->id;
    ofstream out(sout.str().c_str());
    multimap<unsigned, unsigned>::iterator it;
    for (it=m.begin(); it!=m.end(); it++)
      out << setw(9) << dec << it->first 
	  << setw(9) << hex << it->second 
	  << setw(9) << dec << keyproc[it->second]
	  << endl;
  }
}


void pH2OT::Remap_by_effort()
{
  //
  // Make a map containing the keys from the full frontier
  //
  map<unsigned, double> counts;
  //      ^         ^
  // key--/         |
  //                |
  // weight---------/
  //

  //
  // Make a key-counts map from the particle list
  //
  {
    PartMapItr t;
    for (t=cc->Particles().begin(); t!=cc->Particles().end(); t++) {
      t->second.tree = tCell::getKey(&(t->second.pos[0]));
      if (t->second.tree) {
	if (counts.find(t->second.tree) == counts.end())
	  counts[t->second.tree]  = t->second.effort;
	else
	  counts[t->second.tree] += t->second.effort;
      }
    }
  }

  //
  // Sum up the counts for each key on the frontier
  //
  unsigned csize1 = counts.size(), csize0 = 0;
  vector<unsigned> csize(numprocs);
  MPI_Allgather(&csize1, 1, MPI_UNSIGNED, &csize[0], 1, MPI_UNSIGNED,
		MPI_COMM_WORLD);

  for (int n=0; n<numprocs; n++) csize0 = max<unsigned>(csize0, csize[n]);


  //
  // Order the key list by cell weight
  //
  multimap<double, unsigned> m;
  //          ^         ^
  // counts---/         |
  //                    |
  // key----------------/
  //

  vector<unsigned> one(csize0);
  vector<double>   two(csize0);
  for (int n=0; n<numprocs; n++) {
    if (n==myid) {
      map<unsigned, double>::iterator j=counts.begin();
      for (unsigned i=0; i<csize[n]; i++) {
	one[i] = j->first;
	two[i] = j->second;
	j++;
      }
    }
  
    MPI_Bcast(&one[0], csize[n], MPI_UNSIGNED, n, MPI_COMM_WORLD);
    MPI_Bcast(&two[0], csize[n], MPI_DOUBLE,   n, MPI_COMM_WORLD);

    for (unsigned i=0; i<csize[n]; i++)
      m.insert(pair<double, unsigned>(two[i], one[i]));
  }

  //
  // Heuristic balanced partition solution
  // -------------------------------------
  // (1) Put on object in each bin (from biggest to smallest)
  //
  // (2) In reverse, put as many objects in each bin as possible to
  //     bring measure upto the first (largest) bin, but not
  //     exceeding.
  //
  // (3) Repeat (1) and (2) until list is exhausted
  //
  keyproc.clear();
  vector<double> buckets(numprocs, 0.0);
  vector<unsigned> ncell(numprocs, 0);
  multimap<double, unsigned>::reverse_iterator k=m.rbegin(); 
  while (k!=m.rend()) {
    // Step 1
    for (int i=0; i<numprocs; i++) {
      buckets[i] += k->first;
      keyproc[k->second] = i;
      ncell[i]++;
      if (++k == m.rend()) break;
    }
    // Step 2
    for (int i=numprocs-1; i>0 && k!=m.rend(); i--) {
      while (buckets[i]+k->first <= buckets[0]) {
	buckets[i] += k->first;
	keyproc[k->second] = i;
	ncell[i]++;
	if (++k == m.rend()) break;
      };
    }
  }

  //
  // Diagnostic
  //
  if (myid==0) {
    ostringstream sout;
    sout << outdir << runtag << ".counts." << cc->id;
    ofstream out(sout.str().c_str(), ios::app);
    // Compute mean and std. dev.
    vector<double> sum(2, 0), sum2(2, 0);
    for (int i=0; i<numprocs; i++) {
      sum [0] += ncell[i];
      sum2[0] += ncell[i]*ncell[i];
      sum [1] += buckets[i];
      sum2[1] += buckets[i]*buckets[i];
    }
    for (int j=0; j<2; j++) {
      sum[j] /= numprocs;
      if (numprocs>1)
	sum2[j] = sqrt((sum2[j] - sum[j]*sum[j]*numprocs)/(numprocs-1));
      else
	sum2[j] = 0.0;
    }
    out << "# Time=" << tnow << endl
	<< "#  Ncells(mean)=" << sum[0] << "  Ncells(stdev)=" << sum2[0]
	<< endl
	<< "#  Effort(mean)=" << sum[1] << "  Effort(stdev)=" << sum2[1]
	<< endl;
    // Data
    for (int i=0; i<numprocs; i++) {
      out << setw(9)  << i
	  << setw(9)  << ncell[i]
	  << setw(12) << buckets[i]
	  << endl;
    }
  }

}


int pH2OT::getProc(unsigned key)
{
  map<unsigned, int>::iterator q=keyproc.find(key);
  if (q == keyproc.end()) {
#ifdef DEBUG
    cout << "Process " << myid << ": unassigned cell key" << endl;
#endif
    return (key % numprocs);
  }
  return q->second;
}


void pH2OT::makeTree()
{
  unsigned min1=MAXINT, max1=0, my_cells=0, d_val0=0, d_pval0=0;
  d_pval = d_val = 0;
  
  //
  // This runs through every tree on each node
  //
  for (map<unsigned, tCell*>::iterator 
	 t=trees.frontier.begin(); t!=trees.frontier.end(); t++) 
    {
      pTree *p = t->second->ptree;
      p->makeTree();
      min1 = min<unsigned>(p->minc, min1);
      max1 = max<unsigned>(p->maxc, max1);
      my_cells += p->frontier.size();
    }
  

  //
  // Compute cell density, mass, momentum, velocity dispersion
  //

  computeCellStates();

  //
  // Gather up statistics
  //

  MPI_Allreduce(&min1, &min_cell, 1, MPI_UNSIGNED, MPI_MIN, MPI_COMM_WORLD);
  MPI_Allreduce(&max1, &max_cell, 1, MPI_UNSIGNED, MPI_MAX, MPI_COMM_WORLD);

  vector<unsigned> chist1(max_cell-min_cell+1, 0);

  for (unsigned n=min_cell; n<=max_cell; n++) {
    for (map<unsigned, tCell*>::iterator 
	   t=trees.frontier.begin(); t!=trees.frontier.end(); t++) {
      chist1[n-min_cell] += t->second->ptree->getCHist(n);
    }
  }

  chist = vector<unsigned>(max_cell-min_cell+1, 0);
  MPI_Allreduce(&chist1[0], &chist[0], max_cell-min_cell+1, 
		MPI_UNSIGNED, MPI_SUM, MPI_COMM_WORLD);
  
  for (unsigned i=1; i<max_cell-min_cell+1; i++) chist[i] += chist[i-1];

  // Accumulate the total number of cells in all trees
  //
  MPI_Allreduce(&my_cells, &total_cells, 1, MPI_UNSIGNED, MPI_SUM, MPI_COMM_WORLD);

  // Debug bad values
  //
  MPI_Reduce(&d_pval, &d_pval0, 1, MPI_UNSIGNED, MPI_SUM, 0, MPI_COMM_WORLD);
  MPI_Reduce(&d_val,  &d_val0,  1, MPI_UNSIGNED, MPI_SUM, 0, MPI_COMM_WORLD);
  if (myid==0) {
    if (d_pval0 || d_val0)
      cout << "pH2OT::makeTree: bad values in [pkeys, keys]=[" 
	   << d_pval0 << ", "
	   << d_val0  << "] " << endl;
  }


  if (!checkLevelList()) {
    cout << "Process " << myid << ": "
	 << "makeTree: ERROR check level list FAILED on exit!" << endl;
  }
#ifdef DEBUG
  if (myid==0) cout << "Before checkBodycell() in makeTree" 
		    << endl << flush;
#endif
  if (!checkBodycell()) {
    cout << "Process " << myid << ": "
	 << "makeTree: ERROR body cell check FAILED on exit!" << endl;
  }    
#ifdef DEBUG
  if (myid==0) cout << "After checkBodycell() in makeTree" << endl;
#endif
}

void pTree::makeTree()
{
#ifdef USE_GPTL
  GPTLstart("pTree::makeTree");
#endif
  //
  // Clean up
  // 
  frontier.clear();
  bodycell.clear();

				// Remove the current tree
  delete root;

#ifdef DEBUG
  //-----------------------------------------------------------------      
  if (0) {
    ostringstream ostr;
    ostr << outdir << runtag << "_" << hex << tkey << ".pH2OT_storage";
    for (int n=0; n<numprocs; n++) {
      if (myid==n) {
	ofstream out(ostr.str().c_str(), ios::app);
	if (out) {
	  out << setw(18) << tnow
	      << setw(6)  << myid
	      << setw(6)  << tkey
	      << setw(12) << keybods.size()
	      << setw(12) << frontier.size()
	      << setw(12) << bodycell.size()
	      << endl;
	  if (myid==numprocs-1) out << endl;
	  out.close();
	}
      }
      (*ph->barrier)("pTree: pH2OT_storage");
    }
  }
  //-----------------------------------------------------------------      
#endif

  //
  // Make a new root cell
  //
  root = new pCell(this);

  //
  // Add the data
  //
  key_indx::iterator it;
  pCell* p = root;
  for (it=keybods.begin(); it!=keybods.end(); it++)  {

    if (it->first < ph->key_min || it->first >= ph->key_max) {
      cout << "Process " << myid << ": in makeTree, key=" 
	   << hex << it->first << " expected range ["
	   << ph->key_min << ", " << ph->key_max << "]" << dec
	   << endl;
    }
				// Do the work
    p = p->Add(*it);		
  }

#ifdef USE_GPTL
  GPTLstart("pTree::makeTree::getFrontier");
#endif

  // Find min and max cell occupation
  //
  unsigned nt;
  minc = MAXINT;
  maxc = 0;
  for (key_cell::iterator it=frontier.begin(); it != frontier.end(); it++) {
    nt = it->second->bods.size();
    minc = min<unsigned>(minc, nt);
    maxc = max<unsigned>(maxc, nt);
  }

  chist = vector<unsigned>(maxc-minc+1, 0);
  for (key_cell::iterator it=frontier.begin(); it != frontier.end(); it++)
    chist[it->second->bods.size()-minc]++;
  
#ifdef USE_GPTL
  GPTLstop("pTree::makeTree::getFrontier");
  GPTLstart("pTree::makeTree::makeCellLevelList");
#endif

  makeCellLevelList();

#ifdef USE_GPTL
  GPTLstop("pTree::makeTree::makeCellLevelList");
  GPTLstop("pTree::makeTree");
#endif

}

void pTree::adjustCellLevelList(unsigned mlevel)
{
  if (multistep==0) return;	// No need to bother if multistepping is off
				// Otherwise . . . 
#ifdef USE_GPTL
  GPTLstart("pTree::adjustCellLevelList");
#endif

  for (unsigned M=mlevel; M<=multistep; M++) {

    if ( clevels[M].size() > 0 ) {

      set<pCell*>::iterator it = clevels[M].begin(), nit;

      while (it != clevels[M].end()) {
	nit = it++;
	if ( (*nit)->bods.size() ) { 
				// Newly computed level for this cell:
				// we may move cell down but not up . . .
	  unsigned m = max<unsigned>(mlevel, (*nit)->remake_plev());
	  if (M!=m) {
	    clevels[m].insert(*nit);
	    clevlst[*nit] = m;
	    clevels[M].erase(nit);
	  }
	}
	
	if (clevels[M].empty()) break;
      }
    }
  }

#ifdef USE_GPTL
  GPTLstop("pTree::adjustCellLevelList");
#endif
}


void pH2OT::geomUpdate()
{
  // Reset subbox corners
  for (map<unsigned, tCell*>::iterator 
	 t=trees.frontier.begin(); t!=trees.frontier.end(); t++) {

    pTree *tr = t->second->ptree;
    tr->C = tr->pc->getCorner();
  }

  // Recompute the particle partition
  Remap();
}


void pH2OT::computeCellStates()
{
  pTree *tr;
  key_cell::iterator it;
  map<unsigned, tCell*>::iterator t;
  
  timer_getsta1.start();

  // Iterate through the list of trees
  //
  for (t=trees.frontier.begin(); t!=trees.frontier.end(); t++) {

    tr = t->second->ptree;

    // Each node start at root and walk down the tree to zero the
    // counts, and then, walk up the tree to accumulate the state
    //
    // tr->root->computeState();

    // Each node start at root and walk down the tree to zero the
    // counts
    //
    tr->root->zeroState();

    // March through the frontier and accumulate the counts
    //
    for (it=tr->frontier.begin(); it!=tr->frontier.end(); it++) 
      it->second->accumState();
    
  }

  timer_getsta1.stop();
  timer_getsta2.start();

  //
  // Clean up
  //
  trees.ZeroAccum();

  //
  // Each node now broadcasts its root cell states
  //
  const int stateSize = 10;
  int nf = trees.frontier.size()*stateSize, sz=0;
  MPI_Allgather(&nf, 1, MPI_INT, &numfront[0], 1, MPI_INT,
		MPI_COMM_WORLD);

  //
  // Set (1) to (0) for MPI_Allgatherv instead of blocked broadcasts
  //
  if (0) {

    for (int n=0; n<numprocs; n++) sz = max<int>(sz, numfront[n]);
    vector<double> sbuf(sz);

    for (int n=0; n<numprocs; n++) {
      
      if (n==myid) {
	map<unsigned, tCell*>::iterator it;
	unsigned icnt = 0;
	for (it=trees.frontier.begin(); it!=trees.frontier.end(); it++) {
	  for (int k=0; k<stateSize; k++) {
	    if (icnt >= sz) {
	      cerr << "Bad ERROR!!!  About to overwrite in statebuf!" << endl;
	    } else {
	      sbuf[icnt++] = it->second->ptree->root->state[k];
	    }
	  }
	}
      }
      
      MPI_Bcast(&sbuf[0], numfront[n], MPI_DOUBLE, n, MPI_COMM_WORLD);
      
      for (int j=0; j<numfront[n]/stateSize; j++) trees.AddState(&sbuf[j*stateSize]);
    }

  } else {

    vector<double> sbuf(nf);

    map<unsigned, tCell*>::iterator it;
    unsigned icnt = 0;
    for (it=trees.frontier.begin(); it!=trees.frontier.end(); it++) {
      for (unsigned k=0; k<stateSize; k++) {
	sbuf[icnt++] = it->second->ptree->root->state[k];
      }
    }

    unsigned s = 0;
    displace[0] = 0;
    for (int n=0; n<numprocs; n++) {
      if (n>0) displace[n] = displace[n-1] + numfront[n-1];
      s += numfront[n];
    }

    vector<double> sbuf0(s);


    MPI_Allgatherv(&sbuf[0], nf, MPI_DOUBLE, &sbuf0[0], &numfront[0],
		   &displace[0], MPI_DOUBLE, MPI_COMM_WORLD);

    for (int j=0; j<s/stateSize; j++) trees.AddState(&sbuf0[j*stateSize]);
  }

  timer_getsta2.stop();
  timer_getsta3.start();

  //
  // This runs through every tree on each node
  //
  for (t=trees.frontier.begin(); t!=trees.frontier.end(); t++) {
    pTree *p = t->second->ptree;

    // March through the frontier to find the sample cells
    //
    for (key_cell::iterator it=p->frontier.begin(); it != p->frontier.end(); it++) 
      it->second->findSampleCell();
  }
  
  timer_getsta3.stop();
}


void pH2OT::adjustTree(unsigned mlevel)
{
  vector< vector<unsigned long> > exchange(numprocs);

  timer_tadjust.start();

#ifdef DEBUG
  //-----------------------------------------------------------------      
  cout << "Process " 
       << myid << ": adjustTree: before checkLevelList on entrance" 
       << endl << flush;
  if (!checkLevelList()) {
    cout << "Process " << myid << ": adjustTree:"
	 << " ERROR mlevel=" << mlevel
	 << " check level list FAILED on entrance!" << endl;
  } else {
    cout << "Process " << myid << ": adjustTree: mlevel=" << mlevel
	 << " level list OK on entrance!" << endl;
  }

  cout << "Process " << myid 
       << ": adjustTree: before checkBodycell on entrance" 
       << endl << flush;
  if (!checkBodycell()) {
    cout << "Process " << myid << ": adjustTree:"
	 << " ERROR mlevel=" << mlevel 
	 << " body cell check FAILED on entrance!" << endl;
  }    
  //-----------------------------------------------------------------      
#endif

  timer_waiton0.start();   
  (*barrier)("pH2OT: begin adjust timer", __FILE__, __LINE__);
  timer_waiton0.stop();


  //
  // Determine the particles whose cells, trees, and nodes have
  // changed
  //
  for (map<unsigned, tCell*>::iterator 
	 t=trees.frontier.begin(); t!=trees.frontier.end(); t++) 
    {
      t->second->ptree->computeStatus(mlevel, exchange);
    }

  timer_waiton1.start();   
  (*barrier)("pH2OT: compute status timer", __FILE__, __LINE__);
  timer_waiton1.stop();

  //
  // Use the particle exchange list to shift bodies to new procs as
  // necessary
  //
  particleExchange(exchange);
  
  timer_waiton2.start();   
  (*barrier)("pH2OT: particle exchange timer", __FILE__, __LINE__);
  timer_waiton2.stop();

  //
  // Create, remove and delete changed cells
  //
  map<unsigned, tCell*>::iterator t;
  for (t=trees.frontier.begin(); t!=trees.frontier.end(); t++) 
    {
      t->second->ptree->cellUpdate(mlevel);
    }

  //
  // Adjust the cell levels for each tree
  //
  for (t=trees.frontier.begin(); t!=trees.frontier.end(); t++) 
    {
      t->second->ptree->adjustCellLevelList(mlevel);
    }


  timer_tadjust.stop();

  //
  // Compute the physical states in each cell for each tree
  //
  computeCellStates();

#ifdef DEBUG
  //-----------------------------------------------------------------      
  cout << "Process " << myid 
       << ": adjustTree: before checkBodycell on exit" 
       << endl << flush;
  if (!checkBodycell()) {
    cout << "Process " << myid << ": "
	 << "adjustTree: ERROR mlevel=" << mlevel 
	 << " body cell check FAILED on exit!" << endl;
  }    
  if (!checkParticles(cout)) {
    cout << "Process " << myid << ": "
	 << "adjustTree: ERROR mlevel=" << mlevel 
	 << " initial particle check FAILED on exit!" << endl;
  }    
  if (!checkFrontier(cout)) {
    cout << "Process " << myid << ": "
	 << "adjustTree: ERROR mlevel=" << mlevel
	 << ": frontier check FAILED on exit!" << endl;
  }
  
  if (!checkLevelList()) {
    cout << "Process " << myid << ": "
	 << "adjustTree: ERROR mlevel=" << mlevel 
	 << " check level list FAILED on exit!" << endl;
  }
  //-----------------------------------------------------------------      
#endif

}


void pTree::computeStatus(unsigned mlevel, vector< vector<unsigned long> >& exchange)
{
#ifdef USE_GPTL
  GPTLstart("pTree::computeStatus");
#endif
				// Barrier to make sure that the timer
				// gives a sensible measurement of key time
  ph->timer_cstatus.start();

  key_type oldkey;
  unsigned oldpkey;
  list<unsigned long>::iterator ip;
  list<unsigned long> curp;
  pCell* c=0;

  //
  // Exchange list
  //
#ifdef USE_GPTL
  GPTLstart("pTree::keyBods");
#endif

  // 
  // Make body list from frontier cells for this level
  //
  set<pCell*>::iterator it;
  // set<unsigned long>::iterator ib;
  vector<unsigned long>::iterator ib;
  
  ph->timer_bodlist.start();

  for (unsigned M=mlevel; M<=multistep; M++) {
    for (it=clevels[M].begin(); it!=clevels[M].end(); it++) {
      for (ib=(*it)->bods.begin(); ib!=(*it)->bods.end(); ib++) {
	curp.push_back(*ib);
      }
    }
  }
  
  ph->timer_bodlist.stop();

#ifdef USE_GPTL
  GPTLstop ("pTree::keyBods");
  GPTLstart("pTree::keyCells");
#endif

  //
  // Update body by body using the list
  //
  for (ip=curp.begin(); ip!=curp.end(); ip++) {

    ph->timer_bodlist.start();

    Particle *p = ph->cc->Part(*ip);
    if (p==0) {			// Sanity check
      cout << "Process " << myid 
	   << ": pTree::computeStatus: ERROR crazy particle index!" << endl;
    }
    ph->numkeys++;

    //
    // Get and recompute keys
    //
				// Current key & tree
    oldkey  = p->key;
    oldpkey = p->tree;
    key_pair oldpair(oldkey, p->indx);
				// New tree
    p->tree = tCell::getKey(&(p->pos[0]));

    //
    // Get this particle's cell
    //

    key_key ::iterator ij = bodycell.find(oldkey);

    //
    // Missing key?
    //
    if (ij == bodycell.end()) {	//
      cout << "Process " << myid 
	   << ": pTree::computeStatus: ERROR can not find cell for particle"
	   << ", key="		<< hex << oldkey << dec
	   << ", index="	<< p->indx
	   << ", pnumber="	<< ph->cc->Number() 
	   << ", bodycell="	<< bodycell.size() 
	   << ", cur cell="	<< oldpkey 
	   << ", proc="		<< ph->getProc(oldpkey)
	   << endl;
      continue;
    }

    //
    // Missing cell?
    //
    if ( frontier.find(ij->second) == frontier.end() ) {
      cout << "Process " << myid 
	   << ": pTree::computeStatus: ERROR can not find expected cell"
           << " on frontier, oldbody=" << hex << oldkey
	   << " cell="    << bodycell.find(oldkey)->second << dec
	   << " index="   << p->indx 
	   << endl;
      continue;
    } else {
      c = frontier[ij->second];
    }
    
    ph->timer_bodlist.stop();
    ph->timer_celladj.start();

				// Is particle newly out of bounds?
    if (p->tree == 0u) {
				// Remove from current cell
      c->Remove(oldpair, &change);
      p->key = 0u;
#ifdef ADJUST_INFO
      //-----------------------------------------------------------------      
      cout << "Process " << myid << ": REMOVING index=" 
           << p->indx << " cell=" << hex << c->mykey << dec << endl << flush;
      //-----------------------------------------------------------------      
#endif
				// Is tree the same?
    } else if (p->tree == oldpkey) {
				// Get the new body key
      p->key = ph->getKey(&p->pos[0], C);
      key_pair newpair(p->key, p->indx);
      
				// Is the cell the same?
      if (c->isMine(p->key)) {
				// Update key list and body-cell index
	  c->UpdateKeys(oldpair, newpair);
#ifdef ADJUST_INFO
      //-----------------------------------------------------------------      
	  cout << "Process " << myid << ": SAME CELL index=" 
	  << p->indx << " tree=" << tkey << " cell=" << hex << c->mykey
	  << " old=" << oldkey << " new=" << p->key << dec << endl << flush;
      //-----------------------------------------------------------------      
#endif
				// Otherwise, body needs to be relocated in
				// the current tree
      } else {
				// Remove from current cell
	c->Remove(oldpair, &change);
				// Add the particle to its new cell
	c->Add(newpair, &change);
#ifdef ADJUST_INFO
      //-----------------------------------------------------------------      
	cout << "Process " << myid << ": SAME TREE, NEW CELL index=" 
	     << p->indx << " tree=" << tkey << " cell=" << hex << c->mykey
	     << " old=" << oldkey << " new=" << p->key << dec << endl << flush;
      //-----------------------------------------------------------------      
#endif
      }
				// Body needs to be exchanged
    } else {
				// Remove from current cell
      c->Remove(oldpair, &change);
				// Add to exchange list
      exchange[ph->getProc(p->tree)].push_back(p->indx);
#ifdef ADJUST_INFO
      //-----------------------------------------------------------------      
      cout << "Process " << myid << ": EXCHANGING index=" 
	   << p->indx << ", to " << ph->getProc(p->tree) 
	   << endl << flush;
      //-----------------------------------------------------------------      
#endif
    }

    ph->timer_celladj.stop();

  }


#ifdef USE_GPTL
  GPTLstop("pTree::keyCells");
  GPTLstop("pTree::computeStatus");
#endif

  ph->timer_cstatus.stop();
}


bool pH2OT::checkLevelList()
{
  if (!list_check) return true;

  timer_schecks.start();
  bool ok = true;
  map<unsigned, tCell*>::iterator t;
  for (t=trees.frontier.begin(); t!=trees.frontier.end(); t++) 
    {
      if (!t->second->ptree->checkLevelList()) ok = false;
    }
  
  timer_schecks.stop();
  return ok;
}


bool pTree::checkLevelList()
{
  ph->timer_schecks.start();

#ifdef USE_GPTL
  GPTLstart("pTree::checkLevelList");
#endif

  bool ok = true;
  unsigned cnt=0;
  set<pCell*>::iterator it;
  // set<unsigned long>::iterator ib;
  vector<unsigned long>::iterator ib;
  list<unsigned long>::iterator ip;

  //
  // Make cell list from frontier
  //
  unsigned found=0, notfront=0, notlevel=0, badlevel=0;
  set<pCell*> clist;
  for (key_cell::iterator it=frontier.begin(); 
       it != frontier.end(); it++) clist.insert(it->second);

  // 
  // Make body list from frontier cells
  //
  list<unsigned long> curp;
  for (unsigned M=0; M<=multistep; M++) {
    for (it=clevels[M].begin(); it!=clevels[M].end(); it++) {
				// Look for this cell on the frontier
      if (clist.find(*it) == clist.end()) {
	notfront++;
	continue;
      }
				// Look for this cell on the level list
      if (clevlst.find(*it) == clevlst.end()) {
	notlevel++;
	continue;
      } else {
	if (M != clevlst[*it]) {
	  badlevel++;
	}
	found++;
      }
      for (ib=(*it)->bods.begin(); ib!=(*it)->bods.end(); ib++) {
	curp.push_back(*ib);
      }
    }
  }
  
  //
  // Check for all bodies
  //
  for (ip=curp.begin(); ip!=curp.end(); ip++) {

    Particle *p = ph->cc->Part(*ip);

    //
    // Get this particle's cell
    //
    key_key::iterator ij = bodycell.find(p->key);
    
    //
    // Bad key?
    //
    if (ij == bodycell.end()) {
      ok = false;
      cnt++;
    }
  }

#ifdef DEBUG
  //-----------------------------------------------------------------      
  if (cnt) {
    cout << "Process " << myid << ": checkLevelList: " 
	 << " tree id=" << tkey << ", " << cnt 
	 << " unmatched particles" << endl;
  }
  //-----------------------------------------------------------------      
#endif
  
  //-----------------------------------------------------------------      
  if (notfront || notlevel || badlevel) {
    cout << "ERRORS: summary for tkey=" << tkey << ":" << endl
	 << "   Total found = " << found << endl
	 << "  Not on front = " << notfront << endl
	 << "  Not on level = " << notlevel << endl
	 << "     Bad level = " << badlevel << endl;
  }
  //-----------------------------------------------------------------      

  ph->timer_schecks.stop();

#ifdef USE_GPTL
  GPTLstop("pTree::checkLevelList");
#endif

  return ok;
}


void pH2OT::particleExchange(vector< vector<unsigned long> >& exchange)
{
#ifdef USE_GPTL
  GPTLstart("pH2OT::particleExchange");
#endif

  Particle part;
  unsigned Tcnt=0, Fcnt=0, sum, err0=0, err1=0;
  vector<int> sdispls(numprocs), rdispls(numprocs);
  vector<int> sendcounts(numprocs, 0), recvcounts(numprocs, 0);

  timer_xchange.start();

  timer_prepare.start();
  for (int k=0; k<numprocs; k++) {
    if (k==myid) continue;
    sendcounts[k] = exchange[k].size();
    Tcnt += sendcounts[k];
  }

  MPI_Allreduce(&Tcnt, &sum, 1, MPI_UNSIGNED, MPI_SUM, MPI_COMM_WORLD);
  timer_prepare.stop();

  if (sum) {
    
    timer_scatter.start();
    
    MPI_Alltoall(&sendcounts[0], 1, MPI_INT, 
		 &recvcounts[0], 1, MPI_INT,
		 MPI_COMM_WORLD);
    
    for (int k=0; k<numprocs; k++) {
      if (k==0) {
	sdispls[0] = rdispls[0] = 0;
	Fcnt = recvcounts[0];
      } else {
	sdispls[k] = sdispls[k-1] + sendcounts[k-1];
	rdispls[k] = rdispls[k-1] + recvcounts[k-1];
	Fcnt += recvcounts[k];
      }
    }
    
    timer_scatter.stop();
    
    //
    // Exchange particles between processes
    //
    int ps;
    vector<Partstruct> psend(Tcnt), precv(Fcnt);
    
    timer_convert.start();
    for (int toID=0; toID<numprocs; toID++) {
      ps = sdispls[toID];
      for (int i=0; i<sendcounts[toID]; i++) {
	pf.Particle_to_part(psend[ps+i], cc->Particles()[exchange[toID][i]]);
	cc->Particles().erase(exchange[toID][i]);
      }
    }
    timer_convert.stop();
    
    timer_xchange.start();
    
    MPI_Alltoallv(&psend[0], &sendcounts[0], &sdispls[0], pf.Particletype, 
		  &precv[0], &recvcounts[0], &rdispls[0], pf.Particletype, 
		  MPI_COMM_WORLD);
    
    timer_xchange.stop();
    
    timer_convert.start();
    
    for (unsigned i=0; i<Fcnt; i++) {
      pf.part_to_Particle(precv[i], part);
      if (part.mass<=0.0 || isnan(part.mass)) {
	cout << "[adjustTree: crazy mass=" << part.mass 
	     << ", indx=" << part.indx 
	     << ", key=" << hex << part.key << dec << "]";
      }
      
      if (part.tree != 0u) {

	pTree *tr = trees.RetrievePtree(part.tree);
	part.key = getKey(&part.pos[0], tr->C);
	
	key_pair newpair(part.key, part.indx);
	tr->keybods.insert(newpair);
	tr->root->Add(newpair, &tr->change);
      } else part.key = 0u;

      cc->Particles()[part.indx] = part;
    }
    cc->nbodies = cc->particles.size();
    
    
    timer_convert.stop();
  }
  
  
  //
  // Add particles to new trees in the same process
  //
  for (vector<unsigned long>::iterator 
	 it=exchange[myid].begin(); it != exchange[myid].end(); it++) {

    Particle *p = cc->Part(*it);

    if (p->tree != 0u) {

      pTree *tr = trees.RetrievePtree(p->tree);
      p->key = getKey(&p->pos[0], tr->C);
      if (p->key != 0u) {
	key_pair newpair(p->key, p->indx);
	tr->keybods.insert(newpair);
	tr->root->Add(newpair, &tr->change);
      } else {
	err1++;
      }
    }

  }

#ifdef USE_GPTL
  GPTLstop("pH2OT::particleExchange");
#endif

  //
  // DEBUG
  //
  MPI_Reduce(&err1, &err0, 1, MPI_UNSIGNED, MPI_SUM, 0, MPI_COMM_WORLD);
  if (myid==0) {
    if (err0) cout << "pH2OT::particleExchange: " 
		   << err0 << " mistakes"
		   << endl;
  }

  timer_xchange.stop();
}

void pTree::cellUpdate(unsigned mlevel)
{
#ifdef USE_GPTL
  GPTLstart("pTree::cellUpdate");
#endif

  ph->timer_cupdate.start();

  // Create, remove and delete changed cells
  // ---------------------------------------
  // Need to do the creates first in case that cells enter and exit
  // the frontier during the Add() step
  //
  for (list<cell_indx>::iterator 
	 it=change.begin(); it!=change.end(); it++) {
      
    pCell *c = it->first;	// This is the cell pointer

    switch(it->second) {	// Add this one to the active lists
    case CREATE:
      {
	unsigned m = max<unsigned>(c->maxplev, mlevel);
	clevlst[c] = m;
	clevels[m].insert(c);
      }
      break;
      
    case REMOVE:		// Remove this cell from the active list
      {
#ifdef DEBUG
	//-----------------------------------------------------------------
	if (clevlst.find(c) == clevlst.end()) {
	  cout << "pTree::adjustTree: cell=" << hex << c
	       << dec << " not in level list";
	  if (frontier.find(c->mykey) == frontier.end())
	    cout << " and gone from frontier";
	  else
	    cout << " but is in frontier";
	  cout << endl;
	  continue;
	}
	//-----------------------------------------------------------------
#endif
	unsigned m = clevlst[c];
	clevlst.erase(c);
#ifdef DEBUG
	//-----------------------------------------------------------------
	if (clevels[m].find(c) == clevels[m].end()) {
	  cout << "pTree::adjustTree: cell=" << hex << c
	       << dec << " not in level " << m << endl;
	}
	//-----------------------------------------------------------------
#endif
	clevels[m].erase(c);
      }
      break;

    case KILL:			// Delete the cell
      delete c;
      break;

    case RECOMP:		// Find the sample cell
      c->findSampleCell();
      break;

    default:
      cout << "Process " << myid << ": unknown action in pH2OT::adjustTree()!!"
	   << endl;
    }
  }

  change.clear();		// Reset the change list for next time

  ph->timer_cupdate.stop();

#ifdef USE_GPTL
  GPTLstop("pTree::cellUpdate");
#endif
}

bool pTree::checkBodycell()
{
  ph->timer_schecks.start();

  bool ok = true;
  unsigned cnt=0, matched=0;
  PartMapItr n;
  for (n=ph->cc->Particles().begin(); n!=ph->cc->Particles().end(); n++) {
				// Ignore OOB particle
    if (n->second.key==0u)      continue;
				// Is this my particle?
    if (n->second.tree != tkey) continue;
				// Look for the bodycell
    key_key::iterator it = bodycell.find(n->second.key);
    if (it==bodycell.end()) {
      ok = false;
      cnt++;
#ifdef DEBUG
      //-----------------------------------------------------------------
      cout << "Process " << myid << ": checkBodycell:" 
	   << " unmatched particle: tree=" << n->second.tree
	   << " key=" << hex << n->second.key << dec
	   << " index=" << n->second.indx << endl;
      //-----------------------------------------------------------------
#endif
    } else matched++;
  }

  // #ifdef DEBUG
  //-----------------------------------------------------------------
  if (cnt) {
    cout << "Process " << setw(4) << myid << ": checkBodycell: " 
	 << " tree id=" << setw(9) << tkey << ", [" << setw(6) << cnt 
	 << " unmatched][" << setw(6) << matched << " matched] particles" 
	 << endl;
  }
  //-----------------------------------------------------------------
  // #endif

  ph->timer_schecks.stop();
  return ok;
}


bool pH2OT::checkBodycell()
{
  if (!list_check) return true;

  timer_schecks.start();

  bool ok = true;
  map<unsigned, tCell*>::iterator t;
  for (t=trees.frontier.begin(); t!=trees.frontier.end(); t++) 
    {
      if (!t->second->ptree->checkBodycell()) ok = false;
    }

  timer_schecks.stop();
  return ok;
}


unsigned pH2OT::CellCount(double pctl)
{
  pctl = min<double>( max<double>(0.0, pctl), 1.0 );
  unsigned target = static_cast<unsigned>(floor( pctl*chist.back() ));
  unsigned li=0, ui=max_cell-min_cell, mi;
  while (ui-li>1) {
    mi = (li + ui)/2;
    if (chist[mi] <= target)
      li = mi;
    else
      ui = mi;
  }

  if ( chist[ui]-target > target-chist[li] )
    return min_cell+li;
  else
    return min_cell+ui;
}


vector<unsigned> cntlev;
vector<unsigned> kidlev;
vector<double>   maslev;
vector<double>   vollev;


void pTree::densEmit(unsigned lev, pCell *p)
{
  if (p->level == lev) {
    cntlev[lev]++;
    if (p->parent) kidlev[lev] += p->parent->children.size();
    maslev[lev] += p->state[0];
#ifdef I128
    vollev[lev] += ph->volume/(key_type(1u) << (3*p->level)).toDouble();
#else
    vollev[lev] += ph->volume/(key_type(1u) << (3*p->level));
#endif
  } else {
    map<unsigned, pCell*>::iterator it;
    for (it=p->children.begin(); it!=p->children.end(); it++) 
      densEmit(lev, it->second);
  }
}

void pH2OT::densCheck()
{
  timer_schecks.start();

  unsigned MaxLev = 6;
  cntlev = vector<unsigned> (MaxLev+1, 0);
  kidlev = vector<unsigned> (MaxLev+1, 0);
  maslev = vector<double>   (MaxLev+1, 0);
  vollev = vector<double>   (MaxLev+1, 0);

  map<unsigned, tCell*>::iterator t;
  
  for (unsigned lev=0; lev<=MaxLev; lev++) {
    for (int n=0; n<numprocs; n++) {
      if (myid==n) {
				// Walk tree
	for (t=trees.frontier.begin(); t!=trees.frontier.end(); t++) 
	  t->second->ptree->densEmit(lev, t->second->ptree->root);
      }
    }
  }
  
  if (myid==0) {
    cout << endl << "Density check " << setw(60) << setfill('=') << '='
	 << setfill(' ') << endl;
    cout << endl << "Node #" << myid << endl << endl;
    cout << setw(8)  << "Level"
	 << setw(8)  << "Count"
	 << setw(8)  << "Child"
	 << setw(18) << "Mass"
	 << setw(18) << "Volume"
	 << setw(18) << "Vol/cell"
	 << endl
	 << setw(8)  << "------"
	 << setw(8)  << "------"
	 << setw(8)  << "------"
	 << setw(18) << "------"
	 << setw(18) << "------"
	 << setw(18) << "------"
	 << endl;
    
    for (unsigned k=0; k<=MaxLev; k++) {
      if (cntlev[k])
	cout << setw(8)  << k 
	     << setw(8)  << cntlev[k]
	     << setw(8)  << (int)floor((double)kidlev[k]/cntlev[k]+0.5)
	     << setw(18) << maslev[k]
	     << setw(18) << vollev[k]
	     << setw(18) << vollev[k]/cntlev[k]
	     << endl;
    }
    cout << endl;
  }
    
  for (int n=1; n<numprocs; n++) {

    if (myid==0) {
      vector<unsigned> cntlevN(MaxLev+1, 0);
      vector<unsigned> kidlevN(MaxLev+1, 0);
      vector<double>   maslevN(MaxLev+1, 0);
      vector<double>   vollevN(MaxLev+1, 0);

      MPI_Recv(&cntlevN[0], MaxLev+1, MPI_UNSIGNED, n, 144, MPI_COMM_WORLD,
	       MPI_STATUS_IGNORE);
      MPI_Recv(&kidlevN[0], MaxLev+1, MPI_UNSIGNED, n, 145, MPI_COMM_WORLD,
	       MPI_STATUS_IGNORE);
      MPI_Recv(&maslevN[0], MaxLev+1, MPI_DOUBLE,   n, 146, MPI_COMM_WORLD,
	       MPI_STATUS_IGNORE);
      MPI_Recv(&vollevN[0], MaxLev+1, MPI_DOUBLE,   n, 147, MPI_COMM_WORLD,
	       MPI_STATUS_IGNORE);

      cout << endl << "Node #" << n << endl;
      for (unsigned k=0; k<=MaxLev; k++) {
	if (cntlevN[k])
	  cout << setw(8)  << k 
	       << setw(8)  << cntlevN[k]
	       << setw(8)  << (int)floor((double)kidlevN[k]/cntlevN[k]+0.5)
	       << setw(18) << maslevN[k]
	       << setw(18) << vollevN[k]
	       << setw(18) << vollevN[k]/cntlevN[k]
	       << endl;
      }
      cout << endl;
      
    } else if (n==myid) {
      MPI_Send(&cntlev[0], MaxLev+1, MPI_UNSIGNED, 0, 144, MPI_COMM_WORLD);
      MPI_Send(&kidlev[0], MaxLev+1, MPI_UNSIGNED, 0, 145, MPI_COMM_WORLD);
      MPI_Send(&maslev[0], MaxLev+1, MPI_DOUBLE, 0,   146, MPI_COMM_WORLD);
      MPI_Send(&vollev[0], MaxLev+1, MPI_DOUBLE, 0,   147, MPI_COMM_WORLD);
    }
    (*barrier)("pH2OT: density check report", __FILE__, __LINE__);
  }

  vector<unsigned> cntlev0(MaxLev+1, 0);
  vector<unsigned> kidlev0(MaxLev+1, 0);
  vector<double>   maslev0(MaxLev+1, 0);
  vector<double>   vollev0(MaxLev+1, 0);

  MPI_Reduce(&cntlev[0], &cntlev0[0], MaxLev+1, MPI_UNSIGNED, MPI_SUM, 0,
	     MPI_COMM_WORLD);
  MPI_Reduce(&kidlev[0], &kidlev0[0], MaxLev+1, MPI_UNSIGNED, MPI_SUM, 0,
	     MPI_COMM_WORLD);
  MPI_Reduce(&maslev[0], &maslev0[0], MaxLev+1, MPI_DOUBLE,   MPI_SUM, 0,
	     MPI_COMM_WORLD);
  MPI_Reduce(&vollev[0], &vollev0[0], MaxLev+1, MPI_DOUBLE,   MPI_SUM, 0,
	     MPI_COMM_WORLD);

  if (myid==0) {
    cout << endl << "Total" << endl;
    cout << setw(8)  << "Level"
	 << setw(8)  << "Count"
	 << setw(8)  << "Child"
	 << setw(18) << "Mass"
	 << setw(18) << "Volume"
	 << setw(18) << "Vol/cell"
	 << endl
	 << setw(8)  << "------"
	 << setw(8)  << "------"
	 << setw(8)  << "------"
	 << setw(18) << "------"
	 << setw(18) << "------"
	 << setw(18) << "------"
	 << endl;

    for (unsigned n=0; n<=MaxLev; n++) {
      if (cntlev[n])
	cout << setw(8)  << n 
	     << setw(8)  << cntlev0[n]
	     << setw(8)  << (int)floor((double)kidlev0[n]/cntlev0[n]+0.5)
	     << setw(18) << maslev0[n]
	     << setw(18) << vollev0[n]
	     << setw(18) << vollev0[n]/cntlev0[n]
	     << endl;
    }
    cout << endl << setw(74) << setfill('=') << '='
	 << setfill(' ') << endl << endl;
  }

  timer_schecks.stop();
}


void pH2OT::testFrontier(string& filename)
{
  timer_schecks.start();

  key_cell::iterator c;
  pCell *p;

  const unsigned fields = 15;
  vector<unsigned> prec(fields, 18);
  for (unsigned n=0; n<5; n++) prec[n] = 10;
  prec[0] = 14;
  prec[2] = 14;

  vector<ios_base::fmtflags> fmt(fields, ios::dec);
  fmt[1] = ios::hex;

  vector<ios_base::fmtflags> typ(fields, ios::fixed);
  typ[0] = typ[5] = typ[6] = ios::scientific;

  if (myid==0) {

    char labels[][18] = {
      "Time |",
      "Proc |",
      "Key |",
      "Level |",
      "Number |",
      "Mass |",
      "Volume |",
      "Density |",
      "Temp |",
      "Mean X |",
      "Mean Y |",
      "Mean Z |",
      "Mean U |",
      "Mean V |",
      "Mean W |"};
    
    ifstream in(filename.c_str());
    in.close();

    if (in.fail()) {

      ofstream out(filename.c_str());

      out << right;
    
      out << "#" << endl << "#";
      for (unsigned n=0; n<fields; n++)
	out << setw(prec[n]) << setiosflags(fmt[n]) << labels[n];
      out << endl;
      out << "#";
      for (unsigned n=0; n<fields; n++) {
	ostringstream sout;
	sout << n+1 << " |";
	out << setw(prec[n]) << setiosflags(fmt[n]) << sout.str();
      }
      out << endl << "#" << endl;
    }
  }


  map<unsigned, tCell*>::iterator t;

  for (int n=0; n<numprocs; n++) {

    if (n == myid) {

      ofstream out(filename.c_str(), ios::app);

      for (t=trees.frontier.begin(); t!=trees.frontier.end(); t++) {

	for (c =  t->second->ptree->frontier.begin(); 
	     c != t->second->ptree->frontier.end(); c++) {
	  
	  double mass=0, temp=0, pos[]={0,0,0}, vel[]={0,0,0};
	  p = c->second;
	  // set<unsigned long>::iterator ib = p->bods.begin();
	  vector<unsigned long>::iterator ib = p->bods.begin();
	  while (ib != p->bods.end()) {
	    mass += cc->particles[*ib].mass;
	    for (int k=0; k<3; k++) {
	      pos[k] += 
		cc->particles[*ib].mass * 
		cc->particles[*ib].pos[k];
	      vel[k] += 
		cc->particles[*ib].mass * 
		cc->particles[*ib].vel[k];
	      temp += 
		cc->particles[*ib].mass * 
		cc->particles[*ib].vel[k] *
		cc->particles[*ib].vel[k];
	    }
	    ib++;
	  }
	  
	  double v2=0.0;
	  for (int k=0; k<3; k++) {
	    pos[k] /= mass;
	    vel[k] /= mass;
	    v2 += vel[k]*vel[k];
	  }

	  temp = 0.333333333333*(temp/mass - v2);

	  unsigned n=0;
#ifdef I128
	  double vol = volume/( key_type(1u) << (3*p->level)).toDouble(); 
#else
	  double vol = volume/( key_type(1u) << (3*p->level));
#endif
	  
	  out << " ";
	  out << setw(prec[n]) << setiosflags(fmt[n]|typ[n]) << tnow;
	  n++;
	  out << setw(prec[n]) << setiosflags(fmt[n]|typ[n]) << myid;
	  n++;
	  out << setw(prec[n]) << setiosflags(fmt[n]|typ[n]) << hex << c->first << dec;
	  n++;
	  out << setw(prec[n]) << setiosflags(fmt[n]|typ[n]) << p->level;
	  n++;
	  out << setw(prec[n]) << setiosflags(fmt[n]|typ[n]) << p->bods.size();
	  n++;
	  out << setw(prec[n]) << setiosflags(fmt[n]|typ[n]) << mass;
	  n++;
	  out << setw(prec[n]) << setiosflags(fmt[n]|typ[n]) << vol;
	  n++;
	  out << setw(prec[n]) << setiosflags(fmt[n]|typ[n]) << mass/vol;
	  n++;
	  out << setw(prec[n]) << setiosflags(fmt[n]|typ[n]) << temp;
	  n++;
	  for (int k=0; k<3; k++) {
	    out << setw(prec[n]) << setiosflags(fmt[n]|typ[n]) << pos[k];
	    n++;
	  }
	  for (int k=0; k<3; k++) {
	    out << setw(prec[n]) << setiosflags(fmt[n]|typ[n]) << vel[k];
	    n++;
	  }
	  out << endl;
	}
      }
    }

    (*barrier)("HOT: test frontier", __FILE__, __LINE__);
  }

  timer_schecks.stop();
}

void pH2OT::countFrontier(vector<unsigned>& ncells, vector<unsigned>& bodies)
{
  pCell *p;
  key_cell::iterator it;
  map<unsigned, pair<unsigned, unsigned> > data;
  map<unsigned, pair<unsigned, unsigned> >::iterator d;
  unsigned maxLev1=0, maxLev=0;
  
  for (map<unsigned, tCell*>::iterator 
	 t=trees.frontier.begin(); t!=trees.frontier.end(); t++) {
    
    for (it=t->second->ptree->frontier.begin(); 
	 it != t->second->ptree->frontier.end(); it++) {
      p = it->second;
      maxLev1 = max<unsigned>(maxLev1, p->level);
      if ((d=data.find(p->level)) == data.end()) {
	data[p->level] = pair<unsigned, unsigned>(1, p->bods.size());
      } else {
	d->second.first++;
	d->second.second += p->bods.size();
      }
    }
  }

  MPI_Allreduce(&maxLev1, &maxLev, 1, MPI_UNSIGNED, MPI_MAX, 
		MPI_COMM_WORLD);

  vector<unsigned> tcellcnts(maxLev+1, 0), tbodscnts(maxLev+1, 0);
  if (myid==0) {
    ncells = vector<unsigned>(maxLev+1);
    bodies = vector<unsigned>(maxLev+1);
  }

  for (d=data.begin(); d!=data.end(); d++) {
    tcellcnts[d->first] = d->second.first;
    tbodscnts[d->first] = d->second.second;
  }

  MPI_Reduce(&tcellcnts[0], &ncells[0], maxLev+1, MPI_UNSIGNED, MPI_SUM, 
	     0, MPI_COMM_WORLD);

  MPI_Reduce(&tbodscnts[0], &bodies[0], maxLev+1, MPI_UNSIGNED, MPI_SUM,
	     0, MPI_COMM_WORLD);
}

double pH2OT::minVol()
{
  unsigned MaxLev = 0;
  key_cell::iterator it;
  map<unsigned, tCell*>::iterator t;

  for (t=trees.frontier.begin(); t!=trees.frontier.end(); t++) {
    for (it=t->second->ptree->frontier.begin(); 
	 it!=t->second->ptree->frontier.end(); it++)
      MaxLev = max<unsigned>(MaxLev, it->second->level);
  }

  double vol1, vol;
#ifdef I128
  vol1 = volume/(key_type(1u) << (3*MaxLev)).toDouble();
#else
  vol1 = volume/(key_type(1u) << (3*MaxLev));
#endif
  MPI_Allreduce(&vol1, &vol, 1, MPI_DOUBLE, MPI_MIN, MPI_COMM_WORLD);

  return vol;
}

double pH2OT::maxVol()
{
  unsigned MinLev = MAXINT;
  key_cell::iterator it;
  map<unsigned, tCell*>::iterator t;

  for (t=trees.frontier.begin(); t!=trees.frontier.end(); t++) {
    for (it=t->second->ptree->frontier.begin(); 
	 it!=t->second->ptree->frontier.end(); it++)
      MinLev = min<unsigned>(MinLev, it->second->level);
  }

  double vol1, vol;
#ifdef I128
  vol1 = volume/(key_type(1u) << (3*MinLev)).toDouble();
#else
  vol1 = volume/(key_type(1u) << (3*MinLev));
#endif
  MPI_Allreduce(&vol1, &vol, 1, MPI_DOUBLE, MPI_MAX, MPI_COMM_WORLD);

  return vol;
}

double pH2OT::medianVol()
{
  unsigned mlev, num;
  vector<unsigned> lev;
  key_cell::iterator it;
  map<unsigned, tCell*>::iterator t;

  for (t=trees.frontier.begin(); t!=trees.frontier.end(); t++) {
    for (it=t->second->ptree->frontier.begin(); 
	 it!=t->second->ptree->frontier.end(); it++) 
      lev.push_back(it->second->level);
  }

  if (myid==0) {

    for (int n=1; n<numprocs; n++) {
      MPI_Recv(&num, 1, MPI_UNSIGNED, n, 61, MPI_COMM_WORLD, 
	       MPI_STATUS_IGNORE);
      vector<unsigned> lev1(num);
      MPI_Recv(&lev1[0], num, MPI_UNSIGNED, n, 62, MPI_COMM_WORLD,
	       MPI_STATUS_IGNORE);
      for (unsigned j=0; j<num; j++) lev.push_back(lev1[j]);
    }

    sort(lev.begin(), lev.end());
    mlev = lev[(unsigned)floor(0.5*(1+lev.size()))];
  } else {
    num = lev.size();
    MPI_Send(&num, 1, MPI_UNSIGNED, 0, 61, MPI_COMM_WORLD);
    MPI_Send(&lev[0], num, MPI_UNSIGNED, 0, 62, MPI_COMM_WORLD);
  }
  
  MPI_Bcast(&mlev, 1, MPI_UNSIGNED, 0, MPI_COMM_WORLD);
  
#ifdef I128
  return volume/(key_type(1u) << (3*mlev)).toDouble();
#else
  return volume/(key_type(1u) << (3*mlev));
#endif
}

void pH2OT::Repartition(unsigned mlevel)
{
  static bool firstime = true;

#ifdef USE_GPTL
  GPTLstart("pH2OT::Repartition");
  GPTLstart("pH2OT::Repartition::entrance_waiting");
  (*barrier)("pH2OT: repartition entrance wait", __FILE__, __LINE__);
  GPTLstop ("pH2OT::Repartition::entrance_waiting");
#endif

  PartMapItr it;
  
				// Total volume of oct-tree region
  volume = box[0]*box[1]*box[2]; 

				// No need to repartition 
				// if there are no bodies
  if (cc->nbodies_tot==0) {
    if (myid==0) 
      cout << "pH2OT::Repartition with ZERO bodies, continuing" << endl;
#ifdef USE_GPTL
    GPTLstop("pH2OT::Repartition");
#endif
    return;
  }

#ifdef DEBUG
  vector<unsigned long> erased;
#endif

  timer_repartn.start();

  //
  // Recompute keys and compute new processor
  //
#ifdef USE_GPTL
  GPTLstart("pH2OT::Repartition::compute_pkeys");
#endif

  oob.clear();

  int cproc;
  unsigned int pkey;
  vector< vector<unsigned long> > bodylist(numprocs);

  for (it=cc->Particles().begin(); it!=cc->Particles().end(); it++) {

    it->second.tree = pkey = tCell::getKey(&(it->second.pos[0]));

    if (pkey == 0u) {
      oob.insert(it->first);
      it->second.key = 0u;
    } else {
      if (it->second.tree < tottrees) { 
	cout << "Error!! [0] pkey=" << it->second.tree
	     << " < " << tottrees << endl;
      }
      if (it->second.tree >= 2*tottrees) { 
	cout << "Error!! [0] pkey=" << it->second.tree
	     << " > " << 2*tottrees << endl;
      }

      if ( (cproc=getProc(pkey)) != myid) bodylist[cproc].push_back(it->first);
    }
  }

#ifdef DEBUG
  //-----------------------------------------------------------------
  cout << "Process " << left << setw(4) << myid 
       << ": part #=" << setw(9) << cc->Particles().size()
       << "  oob size=" << oob.size() << endl;
  //-----------------------------------------------------------------
#endif

#ifdef USE_GPTL
  GPTLstop ("pH2OT::Repartition::compute_pkeys");
  GPTLstart("pH2OT::Repartition::compute_pkeys_waiting");
  (*barrier)("pH2OT: repartition pkey wait", __FILE__, __LINE__);
  GPTLstop ("pH2OT::Repartition::compute_pkeys_waiting");
  GPTLstart("pH2OT::spreadOOB");
#endif

  spreadOOB();

#ifdef USE_GPTL
  GPTLstop("pH2OT::spreadOOB");
  GPTLstart("pH2OT::bodyList");
#endif

  timer_prepare.start();

  //
  // Nodes compute send list
  //
  unsigned Tcnt=0, Fcnt=0, sum;
  vector<int> sendcounts(numprocs, 0), recvcounts(numprocs, 0);
  vector<int> sdispls(numprocs), rdispls(numprocs);
  
  for (int k=0; k<numprocs; k++) {
    sendcounts[k] = bodylist[k].size();
    Tcnt += sendcounts[k];
  }
  MPI_Reduce(&Tcnt, &sum, 1, MPI_UNSIGNED, MPI_SUM, 0, MPI_COMM_WORLD);

#ifdef USE_GPTL
  GPTLstop ("pH2OT::bodyList");
  GPTLstart("pH2OT::scatter");
#endif

  timer_scatter.start();
  MPI_Alltoall(&sendcounts[0], 1, MPI_INT, 
	       &recvcounts[0], 1, MPI_INT,
	       MPI_COMM_WORLD);
  timer_scatter.stop();

  for (int k=0; k<numprocs; k++) {
    if (k==0) {
      sdispls[0] = rdispls[0] = 0;
      Fcnt = recvcounts[0];
    } else {
      sdispls[k] = sdispls[k-1] + sendcounts[k-1];
      rdispls[k] = rdispls[k-1] + recvcounts[k-1];
      Fcnt += recvcounts[k];
    }
  }

				// DEBUG OUTPUT
  if (false) {			// Set to "true" to enable

    if (myid==0){
      cout <<"--------------------------------------------------------"<< endl
	   <<"---- Send and receive counts for each process "<< endl
	   <<"--------------------------------------------------------"<< endl;
    }

    for (int k=0; k<numprocs; k++) {
      if (myid==k) {
	cout << "Process " << k << endl;
	for (int m=0; m<numprocs; m++) {
	  cout << setw(4) << m << setw(15) << sendcounts[m]
	       << setw(15) << recvcounts[m] << endl;
	}
      }
      (*barrier)("pH2OT: repartition send/receive report", __FILE__, __LINE__);
    }
  }
				// END DEBUG OUTPUT
      
#ifdef USE_GPTL
  GPTLstop ("pH2OT::scatter");
  GPTLstart("pH2OT::exchange");
#endif

				// DEBUG OUTPUT
  if (false) {	 // If true, write send and receive list for each node
    for (int n=0; n<numprocs; n++) {
      if (myid==n) {
	cout<<"--------------------------------------------------------"<<endl
	    <<"---- Repartition: send and receive counts"<<endl
	    <<"--------------------------------------------------------" << endl
	    << "Process " << myid << ": Tcnt=" << Tcnt 
	    << " Fcnt=" << Fcnt << endl;
	for (int m=0; m<numprocs; m++)
	  cout << setw(5) << m 
	       << setw(8) << sendcounts[m]
	       << setw(8) << sdispls[m]
	       << setw(8) << recvcounts[m]
	       << setw(8) << rdispls[m]
	       << endl;
	cout << "--------------------------------------------------------" << endl;
      }
      (*barrier)("pH2OT: repartition send/receive counts", __FILE__, __LINE__);
    }
  }


				// END DEBUG OUTPUT

  //
  // Exchange particles between processes
  //

  int ps;
  vector<Partstruct> psend(Tcnt), precv(Fcnt);
  unsigned tcel;

  timer_convert.start();
  for (int toID=0; toID<numprocs; toID++) {
    ps = sdispls[toID];
    for (int i=0; i<sendcounts[toID]; i++) {
      // Check
      if (cc->Particles().find(bodylist[toID][i]) == cc->Particles().end()) {
	cout << "Not in body list!  Indx=" << bodylist[toID][i] 
	     << "  #" <<  myid;
	if (firstime) cout << " [first time]";
	cout << endl;
      }
      if ((tcel=cc->Particles()[bodylist[toID][i]].tree)>0u) {
	if (tcel > (1u<<(tCell::ndim*pkbits+1)) || 
	    tcel < (1u<<(tCell::ndim*pkbits+0)) ) {
	  cout << "Tree out of bounds  Indx=" << bodylist[toID][i]
	       << "  #" <<  myid << "  Tree=" << tcel;
	  if (firstime) cout << " [first time]";
	  cout << endl;
	}
      }
      // End Check
      pf.Particle_to_part(psend[ps+i], cc->Particles()[bodylist[toID][i]]);
      cc->Particles().erase(bodylist[toID][i]);
    }
  }
  timer_convert.stop();
  timer_xchange.start();

  MPI_Alltoallv(&psend[0], &sendcounts[0], &sdispls[0], pf.Particletype, 
		&precv[0], &recvcounts[0], &rdispls[0], pf.Particletype, 
		MPI_COMM_WORLD);

  timer_xchange.stop();
  timer_convert.start();

  if (Fcnt) {
    Particle part;
    for (unsigned i=0; i<Fcnt; i++) {
      pf.part_to_Particle(precv[i], part);
      if (part.mass<=0.0 || isnan(part.mass)) {
	cout << "[Repartition: crazy mass=" << part.mass
	     << ", indx=" << part.indx 
	     << ", key=" << hex << part.key << dec << "]";
      }
      cc->Particles()[part.indx] = part;
    }
    cc->nbodies = cc->particles.size();
  }
  timer_convert.stop();
  timer_prepare.stop();

  //
  // Remake key body index
  //

  for (map<unsigned, tCell*>::iterator 
	 t=trees.frontier.begin(); t!=trees.frontier.end(); t++)
    t->second->ptree->keybods.clear();

  pTree* tree;
  unsigned oob1_cnt=0, oob_cnt=0;
  for (PartMapItr n=cc->Particles().begin(); n!=cc->Particles().end(); n++) {

    if (n->second.tree==0u) {
      oob1_cnt++;
      continue;
    }

    if (n->second.indx==0u) {
      cout << "pH2OT::Repartition bad particle indx=0!" << endl;
    }

    if (n->second.tree < tottrees) { 
      cout << "Error!! [1] pkey=" << n->second.tree
	   << " < " << tottrees << endl;
    }

    tree = trees.RetrievePtree(n->second.tree);

    n->second.key = getKey(n->second.pos, tree->C);
    
    if (n->second.key==0u) {
      vector<double> xyz, XYZ;
      tree->getRange(xyz, XYZ);
      cout << "pH2OT::Repartition: particle key=0!  Tree=" 
	   << n->second.tree << "  (x, y, z)"
	   << "=(" << (n->second.pos[0] - tree->C[0])/box[0] 
	   << ", " << (n->second.pos[1] - tree->C[1])/box[1] 
	   << ", " << (n->second.pos[2] - tree->C[2])/box[2] << ")" << endl;
      cout << "(" << setw(8) << n->second.pos[0]  << ", "
	     << setw(8) << n->second.pos[1] << ", " << n->second.pos[2] << ")"
	   << endl;
      for (int k=0; k<3; k++)
	cout << "[" << setw(8) << xyz[k] << ", " << setw(8) << XYZ[k] << "] ";
      cout << endl;
    } else {
      tree->keybods.insert(key_pair(n->second.key, n->second.indx));
    }
  }

  MPI_Reduce(&oob1_cnt, &oob_cnt, 1, MPI_UNSIGNED, MPI_SUM, 0, MPI_COMM_WORLD);
  unsigned oob_tot = oobNumber();
  if (myid==0 && oob_cnt != oob_tot)
    cout << endl << "pH2OT::Repartition: " << oob_cnt << " out of bounds," 
	 << " expected " << oob_tot << endl;

  if (false) {			// Sanity checks for bad particle indices
    bool ok = true;		// and bad particle counts
    for (PartMapItr ip=cc->particles.begin(); ip!=cc->particles.end(); ip++) {
      if (ip->second.indx==0) {
	cout << "pH2OT::Repartition BAD particle in proc=" << myid
	     << " key=" << hex << ip->second.key << dec << endl;
	ok = false;
      }
    }
    if (ok) cout << "pH2OT::Repartition: leaving with good indx" << endl;

    if (cc->particles.size() != cc->nbodies) 
      cout << "pH2OT::Repartition: leaving with # mismatch!" << endl;

    unsigned nbodies1 = cc->nbodies, nbodies0=0;
    MPI_Reduce(&nbodies1, &nbodies0, 1, MPI_UNSIGNED, MPI_SUM, 0, MPI_COMM_WORLD);
    if (myid==0) {
      if (nbodies0 != cc->nbodies_tot)
	cout << "pH2OT::Repartition: leaving with total # mismatch!" << endl;
    }

  }

  timer_repartn.stop();

#ifdef USE_GPTL
  GPTLstop("pH2OT::exchange");
  GPTLstop("pH2OT::Repartition");
#endif

  firstime = false;
}

void pTree::makeCellLevelList()
{
				// Make new lists
  clevlst.clear();
  clevels = vector< set<pCell*> >(multistep+1);

  for (key_cell::iterator it=frontier.begin(); it != frontier.end(); it++) {

    // Ignore empty cells . . .
    if (it->second->bods.size() > 0) {
      it->second->remake_plev();
      clevlst[it->second] = it->second->maxplev;
      clevels[it->second->maxplev].insert(it->second);
    }
  }

}

void pH2OT::printCellLevelList(ostream &out)
{
  vector<unsigned> pcnt(multistep+1, 0);
  vector<unsigned> clev(multistep+1, 0);
  unsigned ncells = 0;
  map<unsigned, tCell*>::iterator t;
  map<pCell*, unsigned>::iterator pit;
  
  for (t=trees.frontier.begin(); t!=trees.frontier.end(); t++) {
    for (pit=t->second->ptree->clevlst.begin(); 
	 pit!=t->second->ptree->clevlst.end(); pit++) 
      {
	pcnt[pit->second]++;
	ncells++;
	for (unsigned M=0; M<=multistep; M++) 
	  {
	    clev[M] += t->second->ptree->clevels[M].size();
	  }
      }
  }
  
  out << left << setw(60) << setfill('-') << "-" << endl << setfill(' ')
      << "*** T=" << tnow << "  N=" << cc->Number()
      << "  M=" << ncells << endl
      << setw(10) << "M" << setw(10) << "number" 
      << setw(10) << "counts" << endl;
  for (unsigned M=0; M<=multistep; M++)
    out << setw(10) << M << setw(10) << clev[M]
	<< setw(10) << pcnt[M] << endl;
  out << left << setw(60) << setfill('-') << "-" << endl << setfill(' ');
}


bool pH2OT::checkParticles(ostream& out, bool pc)
{
  if (!list_check) return true;

  unsigned cnt, badb=0, badc=0;
  const string membername = "checkParticles";

  timer_schecks.start();

  // FOR CELL OUTPUT
  map<pCell*, pair<unsigned, unsigned> > outdat;

  for (map<unsigned, tCell*>::iterator 
	 t=trees.frontier.begin(); t!=trees.frontier.end(); t++) {

    for (unsigned M=0; M<=multistep; M++) {
      if (t->second->ptree->clevels[M].size()) {
	cnt = 0;
	for (set<pCell*>::iterator it=t->second->ptree->clevels[M].begin(); 
	     it!=t->second->ptree->clevels[M].end(); it++) {
	  cnt++;
	  // FOR CELL OUTPUT
	  outdat[*it] = pair<unsigned, unsigned>(M, (*it)->bods.size());

	  if ((*it)->bods.size()) {
	    if (pc) {
	      // for (set<unsigned long>::iterator ib=(*it)->bods.begin();
	      for (vector<unsigned long>::iterator ib=(*it)->bods.begin();
		   ib!=(*it)->bods.end(); ib++) {
		if (!cc->Part(*ib)) {
		  out << "pH2OT::checkParticles:: M=" << M << ", bad body at "
		      << cnt << "/" << t->second->ptree->clevels[M].size() 
		      << " cell=" << hex << (*it) << dec 
		      << " ID=" << t->second->ptree->tkey << endl;
		  badb++;
		}
	      }
	    }
	  } else {
	    out << "pH2OT::checkParticles:: M=" << M << ", zero bods at "
		<< cnt << "/" << t->second->ptree->clevels[M].size() 
		<< " cell=" << hex << (*it) << dec
		<< " ID=" << t->second->ptree->tkey << endl;
	    badc++;
	  }
	}
      }
    }
  }
  
  timer_schecks.stop();

  if (badb || badc) {
    out << "Process " << myid << ": pH2OT::checkParticles, bad cell=" 
	<< badc;
    if (pc) out << " bad bod=" << badb;
    out << endl;
    return false;
  }
  else return true;
}


bool pH2OT::checkFrontier(ostream& out)
{
  if (!list_check) return true;

  unsigned bad=0;
  bool good=true;
  map<unsigned, tCell*>::iterator t;

  timer_schecks.start();

  for (t=trees.frontier.begin(); t!=trees.frontier.end(); t++) {
    for (unsigned M=0; M<=multistep; M++) {
      if (t->second->ptree->clevels[M].size()) {
	for (set<pCell*>::iterator it=t->second->ptree->clevels[M].begin(); 
	     it!=t->second->ptree->clevels[M].end(); it++) {
	  if (t->second->ptree->frontier.find((*it)->mykey) == 
	      t->second->ptree->frontier.end()) {
	    out << "pH2OT::checkFrontier error on M=" << M
		<< ", cell=" << hex << *it << dec << endl;
	    bad++;
	    good = false;
	  }
	}
      }
    }
  }
  
  if (bad) {
    out << "Process " << myid << ": pH2OT::checkFrontier, bad cell=" 
	 << bad << endl;
  }

  timer_schecks.stop();

  return good;
}

unsigned pH2OT::oobNumber()
{
  unsigned number=0, number1=oob.size();
  MPI_Reduce(&number1, &number, 1, MPI_UNSIGNED, MPI_SUM, 0, MPI_COMM_WORLD);
  return number;
}

void pH2OT::checkOOB(vector<unsigned>& sendlist)
{
  bool aok = true;
  unsigned bcnt=0;

  timer_schecks.start();

  if (myid==0) {
    vector<unsigned> recvlist(numprocs*numprocs);
    for (int n=1; n<numprocs; n++) {
      MPI_Recv(&recvlist[0], numprocs*numprocs, MPI_UNSIGNED, n, 321, 
	       MPI_COMM_WORLD, MPI_STATUS_IGNORE);
      bool ok = true;
      for (int j=0; j<numprocs*numprocs; j++)
	if (sendlist[j] != recvlist[j]) ok = false;
      if (!ok) {
	cout << "checkOOB: sendlist from #" << n << " is in error" << endl;
	aok = false;
	bcnt++;
      } else {
	cout << "checkOOB: sendlist from #" << n << " is in OK!" << endl;
      }
    }

  } else {
    MPI_Send(&sendlist[0], numprocs*numprocs, MPI_UNSIGNED, 0, 321, 
	     MPI_COMM_WORLD);
  }

  if (myid==0) {
    if (aok) {
      cout << "checkOOB: all sendlists match!" << endl;
    } else {
      cout << "checkOOB: " << bcnt << " sendlists match!" << endl;
    }
  }

  timer_schecks.stop();
}


void pH2OT::spreadOOB()
{
#ifdef USE_GPTL
  GPTLstart("pH2OT::spreadOOB::in_reduce");
#endif
				// 3% tolerance
  const unsigned long tol = 33;

  vector<long> list0(numprocs, 0), delta(numprocs);

  long list1 = oob.size();
  MPI_Allgather(&list1, 1, MPI_LONG, &list0[0], 1, MPI_LONG, MPI_COMM_WORLD);

#ifdef USE_GPTL
  GPTLstop ("pH2OT::spreadOOB::in_reduce");
  GPTLstart("pH2OT::spreadOOB::spread_comp");
#endif

  double tot = 0.5;		// Round off
  for (int n=0; n<numprocs; n++) tot += list0[n];
  long avg = static_cast<long>(floor(tot/numprocs));

  unsigned long maxdif=0;
  map<unsigned, unsigned> nsend, nrecv;
  for (int n=0; n<numprocs; n++) {
    delta[n] = avg - list0[n];
    //
    // Positive delta ===> Receive some
    // Negative delta ===> Send some
    // Zero delta     ===> Just right
    //
    maxdif = max<long>(maxdif, abs(delta[n]));
    if (delta[n]>0) nrecv[n] =  delta[n];
    if (delta[n]<0) nsend[n] = -delta[n];
  }

				// Debug output
  if (false) {
    ostringstream sout;
    sout << "In spreadOOB, maxdif=" << maxdif << " #=" << cc->nbodies_tot
	 << " ns=" << nsend.size() << " rs=" << nrecv.size();
  }

#ifdef USE_GPTL
  GPTLstop("pH2OT::spreadOOB::spread_comp");
#endif

				// Don't bother if changes are small
  if (maxdif < static_cast<unsigned long>(cc->nbodies_tot)/tol) return;

				// Nothing to send or receive
  if (nsend.size()==0 || nrecv.size()==0) return;

#ifdef USE_GPTL
  GPTLstart("pH2OT::spreadOOB::make_list");
#endif

  map<unsigned, unsigned>::iterator isnd = nsend.begin();
  map<unsigned, unsigned>::iterator ircv = nrecv.begin();
  vector<unsigned> sendlist(numprocs*numprocs, 0);

  while (1) {
    sendlist[numprocs*isnd->first + ircv->first]++;
    if (--(isnd->second) == 0) isnd++;
    if (--(ircv->second) == 0) ircv++;
    if (isnd == nsend.end() || ircv == nrecv.end()) break;
  }

  // DEBUG
  checkOOB(sendlist);
  // END DEBUG

  unsigned Tcnt=0, Fcnt=0;
  for (int i=0; i<numprocs; i++) {
    Tcnt += sendlist[numprocs*myid + i];
    Fcnt += sendlist[numprocs*i + myid];
  }

#ifdef USE_GPTL
  GPTLstop("pH2OT::spreadOOB::make_list");
#endif

  // DEBUG (set to "true" to enable send/recv list diagnostic)
  //
  if (true) {
    if (myid==0) {
      cout <<"--------------------------------------------------------"<< endl
	   <<"---- spreadOOB" << endl
	   <<"--------------------------------------------------------"<< endl
	   << setw(70) << setfill('-') << "-" << endl << setfill(' ')
	   << left << setw(5) << "P" << setw(10) << "Orig"
	   << setw(10) << "To" << setw(10) << "From" 
	   << setw(10) << "Target"
	   << endl
	   << setw(70) << setfill('-') << "-" << endl << setfill(' ');    
      int tto=0, tfr=0;
      for (int j=0; j<numprocs; j++) {
	cout << setw(5) << j << setw(10) << list0[j];
	int to=0, fr=0;
	for (int i=0; i<numprocs; i++) {
	  to += sendlist[numprocs*j + i];
	  fr += sendlist[numprocs*i + j];
	}
	cout << setw(10) << to << setw(10) << fr 
	     << setw(10) << list0[j] - to + fr << endl;
	tto += to;
	tfr += fr;
      }
      cout << setw(70) << setfill('-') << "-" << endl << setfill(' ')
	   << left << setw(5) << "T" 
	   << setw(10) << tto << setw(10) << tfr << setw(10) << tto-tfr << endl
	   << setw(70) << setfill('-') << "-" << endl << setfill(' ');    
    }
  }
  // END DEBUG

#ifdef USE_GPTL
  GPTLstart("pH2OT::spreadOOB::exchange_particles");
#endif

  unsigned ps=0, pr=0;
  set<indx_type>::iterator ioob;
  vector<Partstruct> psend(Tcnt), precv(Fcnt);
  vector<MPI_Request> rql;
  MPI_Request r;
  int ierr;

  //
  // Exchange particles between processes
  //
  for (int frID=0; frID<numprocs; frID++) {
    for (int toID=0; toID<numprocs; toID++) {
				// 
				// Current process sends particles
      if (myid==frID) {		// 
	unsigned To = sendlist[numprocs*frID+toID];
	if (To) {
	  for (unsigned i=0; i<To; i++) {
	    ioob = oob.begin();
	    pf.Particle_to_part(psend[ps+i], cc->Particles()[*ioob]);
	    cc->Particles().erase(*ioob);
	    if (oob.find(*ioob) == oob.end())
	      cerr << "Process " << myid << ": serious error, oob="
		   << *ioob << endl;
	    else oob.erase(ioob);
	  }
	  rql.push_back(r);
	  if ( (ierr=MPI_Isend(&psend[ps], To, pf.Particletype, 
			       toID, 49, MPI_COMM_WORLD, &rql.back()))
	       != MPI_SUCCESS) {
	    cout << "Process " << myid << ": error in spreadOOP sending "
		 << To << " particles to #" << toID 
		 << " ierr=" << ierr << endl;
	  }
	  ps += To;
	}
      }
				// 
				// Current process receives particles (blocking)
      if (myid==toID) {		// 
	unsigned From = sendlist[numprocs*frID+toID];
	if (From) {
	  rql.push_back(r);
	  if ( (ierr=MPI_Irecv(&precv[pr], From, pf.Particletype, 
			       frID, 49, MPI_COMM_WORLD, &rql.back())) 
	       != MPI_SUCCESS)
	    {
	      cout << "Process " << myid << ": error in spreadOOP receiving "
		   << From << " particles from #" << frID 
		   << " ierr=" << ierr << endl;
	    }
	  pr += From;
	}
      }

    } // Receipt loop
  }

  //
  // Wait for completion of sends and receives
  //

  if ( (ierr=MPI_Waitall(rql.size(), &rql[0], MPI_STATUSES_IGNORE)) != MPI_SUCCESS ) 
    {
      cout << "Process " << myid << ": error in spreadOOB Waitall"
	   << ", ierr=" << ierr << endl;
    }

#ifdef USE_GPTL
  GPTLstop ("pH2OT::spreadOOB::exchange_particles");
  GPTLstart("pH2OT::spreadOOB::add_to_particles");
#endif

  //
  // Add particles
  //

  if (Fcnt) {
    Particle part;
    for (unsigned i=0; i<Fcnt; i++) {
      pf.part_to_Particle(precv[i], part);
      if (part.mass<=0.0 || isnan(part.mass)) {
	cout << "[spreadOOB crazy mass indx=" << part.indx 
	     << ", key=" << hex << part.key << dec << "]";
      }
      cc->Particles()[part.indx] = part;
      oob.insert(part.indx);
    }
    cc->nbodies = cc->particles.size();
  }

#ifdef USE_GPTL
  GPTLstop("pH2OT::spreadOOB::add_to_particles");
#endif

}

void pH2OT::checkCellTree()
{
  CellDiag d(trees.frontier.size()), *buf=0;
  const unsigned nh = 8;
  const unsigned cnts[nh] = {3, 6, 10, 20, 30, 100, 300, 3000};
  vector<unsigned> nums(nh+1, 0), hist(nh+1);
  map<unsigned, unsigned> levS;
  unsigned nt, nb;
  key_cell::iterator it, itb, ite;

  timer_schecks.start();

  for (map<unsigned, tCell*>::iterator 
	 t=trees.frontier.begin(); t!=trees.frontier.end(); t++) {

				// Mean and variance for the number of
				// cells per tree
    nt = t->second->ptree->frontier.size();
    d.navg += nt;
    d.nvar += nt*nt;

				// Mean, variance, min, max for the
				// number of bodies per cell in this
				// tree
    itb = t->second->ptree->frontier.begin();
    ite = t->second->ptree->frontier.end();
    for (it=itb; it!=ite; it++) {
      nb = it->second->bods.size();
      if (nb>0) {
				// Tally the level
	if (levS.find(it->second->level) == levS.end())
	  levS[it->second->level] = 1;
	else
	  levS[it->second->level]++;
				// Min and max
	d.mind = min<unsigned>(d.mind, it->second->level);
	d.maxd = max<unsigned>(d.maxd, it->second->level);
	d.bcel++;
	d.btot += nb;
	d.bavg += nb;
	d.bvar += nb*nb;
	for (unsigned i=0; i<nh; i++) {
	  if (nb<=cnts[i]) {
	    nums[i]++;
	    break;
	  }
	}
	if (nb > cnts[nh-1]) nums[nh]++;

      }
      
    }
    
  }

				// Mean and var of cells per tree
  if (d.ntre) {
    d.navg /= d.ntre;
    if (d.ntre>1) d.nvar = sqrt( (d.nvar - d.ntre*d.navg*d.navg)/(d.ntre-1) );
  }

  if (d.bcel) {
    d.bavg /= d.bcel;
    if (d.bcel>1) d.bvar = sqrt( (d.bvar - d.bcel*d.bavg*d.bavg)/(d.bcel-1) );
  }

				// Compute the maximum level for all processes
				//
  unsigned levN=1, levT;
  if (levS.size()) levN = levS.rbegin()->first + 1;
  MPI_Allreduce(&levN, &levT, 1, MPI_UNSIGNED, MPI_MAX,  MPI_COMM_WORLD);

				// Reduce the cell level tally to 
				// the root process
  vector<unsigned> lev1(levT, 0), lev0(levT);
  for (map<unsigned, unsigned>::iterator is=levS.begin(); is!=levS.end(); is++)
    lev1[is->first] = is->second;

  MPI_Reduce(&lev1[0], &lev0[0], levT, MPI_UNSIGNED, MPI_SUM, 0, MPI_COMM_WORLD);
  
				// Reduce the cell body count histogram to 
				// the root process
  MPI_Reduce(&nums[0], &hist[0], nh+1, MPI_UNSIGNED, MPI_SUM, 0, MPI_COMM_WORLD);

  if (myid) {
    MPI_Gather(&d, 1, CellDiagType, buf, 1,   CellDiagType, 0, MPI_COMM_WORLD);
  } else {
    vector<CellDiag> r(numprocs);

    MPI_Gather(&d, 1, CellDiagType, &r[0], 1, CellDiagType, 0, MPI_COMM_WORLD);


    string   sout = outdir + runtag + ".tree_stat";
    ofstream out(sout.c_str(), ios::app | ios::out);

    if (out) {

      const int nf = 10*5 + 15*4;

      out << setw(nf) << setfill('-') << "-" << endl << setfill(' ')
	  <<"---- Tree statistics, T=" << tnow << endl
	  << setw(nf) << setfill('-') << "-" << endl << setfill(' ')
	  << left << setw(5) << "ID"   << setw(10) << "N(tree)"
	  << setw(10) << "N(cell)"     << setw(10) << "N(body)"
	  << setw(10) << "Min(lev)"    << setw(10) << "Max(lev)"
	  << setw(15) << "Avg(N cell)" << setw(15) << "Std(N cell)"
	  << setw(15) << "Avg(N body)" << setw(15) << "Std(N body)"
	  << endl
	  << setw(nf) << setfill('-') << "-" << endl << setfill(' ');
      for (int j=0; j<numprocs; j++) {
	if (r[j].bcel==0) r[j].mind = r[j].maxd = 0;
	out << setw(5) << j 
	    << setw(10) << r[j].ntre 
	    << setw(10) << r[j].bcel
	    << setw(10) << r[j].btot
	    << setw(10) << r[j].mind
	    << setw(10) << r[j].maxd
	    << setw(15) << r[j].navg
	    << setw(15) << r[j].nvar
	    << setw(15) << r[j].bavg
	    << setw(15) << r[j].bvar
	    << endl;
      }
      out << setw(nf) << setfill('-') << "-" << endl << setfill(' ');
      out << "Cell body counts:" << endl
	  << "---------------- " << endl;
      for (unsigned i=0; i<nh; i++) 
	out << "**  " << setw(8) << cnts[i] << setw(10) << hist[i] << endl;
      out << "** " << ">" << setw(8) << cnts[nh-1]
	  << setw(10) << hist[nh] << endl;
      out << setw(nf) << setfill('-') << "-" << endl << setfill(' ');
      out << "Cell level tally:" << endl
	  << "---------------- " << endl;
      for (unsigned i=0; i<levT; i++) 
	out << "**  " << setw(8) << i << setw(10) << lev0[i] << endl;
      out << setw(nf) << setfill('-') << "-" << endl << setfill(' ');
    } else {
      cerr << "Error opening <" << sout << ">" << endl;
    }
  }

  timer_schecks.stop();

}

unsigned pH2OT::checkNumber()
{
  unsigned nbods1=0, nbods=0;
  
  timer_schecks.start();

  key_cell::iterator it, itb, ite;

  for (map<unsigned, tCell*>::iterator 
	 t=trees.frontier.begin(); t!=trees.frontier.end(); t++) {

    itb = t->second->ptree->frontier.begin();
    ite = t->second->ptree->frontier.end();
    for (it=itb; it!=ite; it++)
      nbods1 += it->second->bods.size();
  }
    
  MPI_Reduce(&nbods1, &nbods, 1, MPI_UNSIGNED, MPI_SUM, 0, MPI_COMM_WORLD);

  timer_schecks.stop();

  return nbods;
}


void pH2OT::CollectTiming()
{
  float fval;

  fval  = timer_cstatus.getTime()();
  MPI_Gather(&fval, 1, MPI_FLOAT, &cstat3[0], 1, MPI_FLOAT, 0, MPI_COMM_WORLD);

  fval = timer_xchange.getTime()();
  MPI_Gather(&fval, 1, MPI_FLOAT, &exchg3[0], 1, MPI_FLOAT, 0, MPI_COMM_WORLD);

  fval  = timer_prepare.getTime()();
  MPI_Gather(&fval, 1, MPI_FLOAT, &prepr3[0], 1, MPI_FLOAT, 0, MPI_COMM_WORLD);

  fval = timer_convert.getTime()();
  MPI_Gather(&fval, 1, MPI_FLOAT, &cnvrt3[0], 1, MPI_FLOAT, 0, MPI_COMM_WORLD);

  fval = timer_tadjust.getTime()();
  MPI_Gather(&fval, 1, MPI_FLOAT, &tadjt3[0], 1, MPI_FLOAT, 0, MPI_COMM_WORLD);

  fval = timer_cupdate.getTime()();
  MPI_Gather(&fval, 1, MPI_FLOAT, &updat3[0], 1, MPI_FLOAT, 0, MPI_COMM_WORLD);

  fval = timer_scatter.getTime()();
  MPI_Gather(&fval, 1, MPI_FLOAT, &scatr3[0], 1, MPI_FLOAT, 0, MPI_COMM_WORLD);

  fval = timer_repartn.getTime()();
  MPI_Gather(&fval, 1, MPI_FLOAT, &reprt3[0], 1, MPI_FLOAT, 0, MPI_COMM_WORLD);

  fval = timer_keybods.getTime()();
  MPI_Gather(&fval, 1, MPI_FLOAT, &keybd3[0], 1, MPI_FLOAT, 0, MPI_COMM_WORLD);

  fval = timer_schecks.getTime()();
  MPI_Gather(&fval, 1, MPI_FLOAT, &schks3[0], 1, MPI_FLOAT, 0, MPI_COMM_WORLD);

  fval = timer_waiton0.getTime()();
  MPI_Gather(&fval, 1, MPI_FLOAT, &wait03[0], 1, MPI_FLOAT, 0, MPI_COMM_WORLD);

  fval = timer_waiton1.getTime()();
  MPI_Gather(&fval, 1, MPI_FLOAT, &wait13[0], 1, MPI_FLOAT, 0, MPI_COMM_WORLD);

  fval = timer_waiton2.getTime()();
  MPI_Gather(&fval, 1, MPI_FLOAT, &wait23[0], 1, MPI_FLOAT, 0, MPI_COMM_WORLD);

  fval = timer_bodlist.getTime()();
  MPI_Gather(&fval, 1, MPI_FLOAT, &bodls3[0], 1, MPI_FLOAT, 0, MPI_COMM_WORLD);

  fval = timer_celladj.getTime()();
  MPI_Gather(&fval, 1, MPI_FLOAT, &cladj3[0], 1, MPI_FLOAT, 0, MPI_COMM_WORLD);

  fval = timer_getsta1.getTime()();
  MPI_Gather(&fval, 1, MPI_FLOAT, &cells1[0], 1, MPI_FLOAT, 0, MPI_COMM_WORLD);

  fval = timer_getsta2.getTime()();
  MPI_Gather(&fval, 1, MPI_FLOAT, &cells2[0], 1, MPI_FLOAT, 0, MPI_COMM_WORLD);

  fval = timer_getsta3.getTime()();
  MPI_Gather(&fval, 1, MPI_FLOAT, &cells3[0], 1, MPI_FLOAT, 0, MPI_COMM_WORLD);

  fval = barrier->getTime();
  MPI_Gather(&fval, 1, MPI_FLOAT, &barri3[0], 1, MPI_FLOAT, 0, MPI_COMM_WORLD);

  MPI_Gather(&numkeys, 1, MPI_UNSIGNED, &numk3[0], 1, MPI_UNSIGNED, 0, MPI_COMM_WORLD);

  timer_cstatus.reset();
  timer_xchange.reset();
  timer_prepare.reset();
  timer_convert.reset();
  timer_cupdate.reset();
  timer_tadjust.reset();
  timer_scatter.reset();
  timer_repartn.reset();
  timer_keybods.reset();
  timer_schecks.reset();
  timer_waiton0.reset();
  timer_waiton1.reset();
  timer_waiton2.reset();
  timer_bodlist.reset();
  timer_celladj.reset();
  timer_getsta1.reset();
  timer_getsta2.reset();
  timer_getsta3.reset();

  numkeys = 0;
}


template <typename T> 
void pH2OT::getQuant(vector<T>& in, vector<T>& out)
{
  sort(in.begin(), in.end());
  out = vector<T>(3);
  for (int k=0; k<3; k++)
    out[k] = in[static_cast<int>(floor(in.size()*0.01*qtile[k]))];
}

void pH2OT::Timing(vector<float>    &cstatus, vector<float>    &exchange, 
		   vector<float>    &prepare, vector<float>    &convert, 
		   vector<float>    &tadjust, vector<float>    &update,
		   vector<float>    &scatter, vector<float>    &repartn,
		   vector<float>    &keybods, vector<float>    &schecks, 
		   vector<float>    &waiton0, vector<float>    &waiton1,
		   vector<float>    &waiton2, vector<float>    &bodlist, 
		   vector<float>    &celladj, vector<float>    &getsta1,
		   vector<float>    &getsta2, vector<float>    &getsta3,
		   vector<float>    &treebar, vector<unsigned> &numk)
{
  getQuant<float   >(cstat3, cstatus);
  getQuant<float   >(exchg3, exchange);
  getQuant<float   >(prepr3, prepare);
  getQuant<float   >(cnvrt3, convert);
  getQuant<float   >(tadjt3, tadjust);
  getQuant<float   >(updat3, update);
  getQuant<float   >(scatr3, scatter);
  getQuant<float   >(reprt3, repartn);
  getQuant<float   >(keybd3, keybods);
  getQuant<float   >(schks3, schecks);
  getQuant<float   >(wait03, waiton0);
  getQuant<float   >(wait13, waiton1);
  getQuant<float   >(wait23, waiton2);
  getQuant<float   >(bodls3, bodlist);
  getQuant<float   >(cladj3, celladj);
  getQuant<float   >(cells1, getsta1);
  getQuant<float   >(cells2, getsta2);
  getQuant<float   >(cells3, getsta3);
  getQuant<float   >(barri3, treebar);
  getQuant<unsigned>(numk3,  numk);
}

double pH2OT::totalKE(double& KEtot, double& KEdsp)
{
  vector<double> state(10);

  // Test
  //
  vector<double> state1(10, 0.0);

  map<unsigned, tCell*>::iterator t;
  key_cell::iterator it, itb, ite;
  
  for (t=trees.frontier.begin(); t!=trees.frontier.end(); t++) {

    itb = t->second->ptree->frontier.begin(); 
    ite = t->second->ptree->frontier.end();
    for (it = itb; it != ite; it++) 
      {
	for (unsigned k=0; k<10; k++) 
	  state1[k] += it->second->state[k];
      }
  }

  MPI_Reduce(&state1[0], &state[0], 10, MPI_DOUBLE, MPI_SUM,
	     0, MPI_COMM_WORLD);
  //
  // End test

  KEtot = KEdsp = 0.0;

  double mass  = state[0];
  double *vel2 = &state[1];
  double *vel1 = &state[4];

  if (mass>0.0) {
    for (int k=0; k<3; k++) {
      KEtot += 0.5*vel2[k];
      KEdsp += 0.5*(vel2[k] - vel1[k]*vel1[k]/mass);
    }

    KEtot /= mass;
    KEdsp /= mass;
  }

  return mass;
}

void pH2OT::totalMass(unsigned& Count, double& Mass)
{
  unsigned count1 = 0,   count;
  double   mass1  = 0.0, mass;

  map<unsigned, tCell*>::iterator t;
  for (t=trees.frontier.begin(); t!=trees.frontier.end(); t++)
    {
				// Sanity check for debugging
      if ( t->first<tottrees || t->first>=2*tottrees )  {
	cout << "Process " << myid << ": STRANGE tree id=" << t->first << endl;
	if (1) {
	  const unsigned SIZE = 400;
	  void *buffer[SIZE];

	  int nptrs = backtrace(buffer, SIZE);
	  char **strings = backtrace_symbols(buffer, nptrs);
	  
	  if (strings == NULL) {
	    cout << "Process " << myid << ": error getting backtrace" << endl;
	  } else {
	    cout << "Process " << myid << ": backtrace" << endl;
	    for (int j=0; j<nptrs; j++)
	      cout << "*** <" << setw(3) << myid << "> " << strings[j] << endl;
	    cout << flush;
	    free(strings);
	    
	    cout << "Process " << myid << ": trying to continue" << endl;
	  }
	}
      } else {
	if (t->second != 0) {
	  count1  += t->second->ptree->root->count;
	  mass1   += t->second->ptree->root->state[0];
	} else {
	  cout << "Process " << myid << ": NO root cell in tree=" 
	       << t->first << endl;
	}
      }
    }

  MPI_Reduce(&count1, &count, 1, MPI_UNSIGNED, MPI_SUM, 0, MPI_COMM_WORLD);
  MPI_Reduce(&mass1,  &mass,  1, MPI_DOUBLE,   MPI_SUM, 0, MPI_COMM_WORLD);

  Count = count;
  Mass  = mass;
}


unsigned pH2OT_iterator::nextCell() 
{ 
				// Return if no more trees
  if (tr==pr->trees.frontier.end()) return 0;

				// Set to first tree's frontier to
				// begin, otherwise increment the
				// current frontier iterator
  if (first) {
    fit = tr->second->ptree->frontier.begin();
    first = false;
  } else {
    fit++;
  }
				// If we are and the end of the
				// current tree's frontier, increment
				// the tree iterator and look for a
				// non-empty frontier
  if (fit==tr->second->ptree->frontier.end()) {
    do {
      tr++;
      if (tr==pr->trees.frontier.end()) return 0;
    } while (tr->second->ptree->frontier.size()==0);
				// Found at least one cell . . .
    fit = tr->second->ptree->frontier.begin();
  }

  return fit->second->bods.size(); 
}


unsigned tCell::pbits   = 0;
tTree*   tCell::tr      = 0;
unsigned tCell::ndim    = 3;
unsigned tCell::idim[3] = {0, 1, 2};


tCell::tCell(tTree* tree, unsigned nb)
{
  live++;

  tr      = tree;
  pbits   = nb;
  owner   = 0;
				// I am the root node
  parent  = 0;
  mykey   = 1u;
  level   = 0;
  count   = 0;
  ptree   = 0;
				// My body mask
  mask    = mykey << ndim*(pbits - level);
				// Initialize state
  state   = vector<double>(10, 0.0);
}

tCell::tCell(tCell* mom, unsigned id) : parent(mom)
{
  live++;

  owner   = 0;
				// My map key
  mykey   = (parent->mykey << ndim) + id;
				// My level
  level   = parent->level + 1;
				// My body mask
  mask    = mykey << ndim*(pbits - level);
				// Initialize state
  state   = vector<double>(10, 0.0);

  count   = 0;
  ptree   = 0;
}

tCell::~tCell()
{
  live--;

  // Recursively kill all the cells
  for (map<unsigned, tCell*>::iterator it=children.begin(); 
       it!=children.end(); it++) delete it->second;
}


// Change the dimensions
void tCell::setDimensions(const string& dims)
{
  ndim = dims.length();
  for (unsigned k=0; k<ndim; k++) {
    switch (dims[k]) {
    case 'x':
    case 'X':
      idim[k] = 0;
      break;
    case 'y':
    case 'Y':
      idim[k] = 1;
      break;
    case 'z':
    case 'Z':
      idim[k] = 2;
      break;
    default:
      throw string("Dimension indicator must be in [xyzXYZ]");
    };
  }
}

unsigned tCell::getKey(const double *p)
{
  const double tol = 1.0e-12;

  // Out of bounds?
  //
  double z[ndim];
  unsigned j;
  for (unsigned k=0; k<ndim; k++) {
    j = idim[k];
    z[k] = (p[j] + pH2OT::offset[j])/pH2OT::sides[j];
				// Deal with some roundoff/truncation error
    if (z[j] < 0.0 && z[j] > -tol)  
      z[j] = 0.0;
    else if (z[j] >= 1.0 && z[j] < 1.0+tol)
      z[j] = 1.0 - tol;
    else if (z[j] < 0.0 || z[j] >=1.0) {
      return 0u;
    }
  }

  const double factor = static_cast<double>(1u<<pbits);
  const unsigned mask = 0x1u;

  vector<unsigned> bins(ndim, 0u);

				// Reverse the order
  for (unsigned k=0; k<ndim; k++)
    bins[ndim-1-k] = unsigned( floor(z[idim[k]]*factor) );
  
  unsigned place = 1u;
  unsigned  _key = 0u;

  for (unsigned i=0; i<pbits; i++) {
    for (unsigned k=0; k<ndim; k++) {
      _key |= (bins[k] & mask)*place;
      place = place << 1;
      bins[k] = bins[k] >> 1;
    }
  }

  _key += place;		// Leading placeholder for cell masking

  return _key;
}


tCell* tCell::findNode(const double *p)
{
  unsigned key = getKey(p);
  if (key) return findNode(key);
  else return NULL;
}


tCell* tCell::findNode(const unsigned& key)
{
				// Check that this key belongs to this branch
  unsigned sig = (unsigned)(key - mask) >> ndim*(pbits-level);
  
  if (sig) {
    
    if (parent == 0) {
      cout << "tCell::findNode: impossible condition, process " 
	   << myid << ": level=" << level 
	   << hex << " key=" << key << endl
	   << hex << " sig=" << sig << endl << dec;
    }

    return parent->findNode(key);
  }

				// You found me!
  if (children.size() == 0) return this;
				// Which child
  unsigned key2 = childId(key);
				// Not in my tree?
  if (children.find(key2) == children.end()) return 0;

				// Look for node amongst children
  return children[key2]->findNode(key);
}


tCell* tCell::findAddNode(const double *p)
{
  unsigned key = getKey(p);
  if (key) return findAddNode(key);
  else return NULL;
}


tCell* tCell::findAddNode(const unsigned& key)
{
				// Check that this key belongs to this branch
  unsigned sig = (unsigned)(key - mask) >> ndim*(pbits-level);
  
  if (sig) {
    
    if (parent == 0) {
      cout << "tCell::findAddNode: impossible condition, process " 
	   << myid << ": level=" << level 
	   << hex << " key=" << key << endl
	   << hex << " sig=" << sig << endl << dec;
    }

    return parent->findAddNode(key);
  }

				// You found me!
  if (key == mykey) return this;

  if (level >= tr->maxlevel) {
    cout << "Crazy level=" << level << endl;
  }

				// Which child
  unsigned key2 = childId(key);
				// Create the node if is does not
				// exist
  if (children.find(key2) == children.end()) {
    children[key2] = new tCell(this, key2);
  }

				// Look for node amongst children
  return children[key2]->findAddNode(key);
}


void tCell::getRange(vector<double>& xyz, vector<double>& XYZ)
{
  xyz = vector<double>(3, 0);
  XYZ = vector<double>(3, 0);

  unsigned _place = 1u;
  unsigned _key   = mykey;

  for (unsigned i=0; i<level; i++) {
    for (unsigned j=0; j<ndim; j++) {
      unsigned k=idim[j];
      xyz[2-k] += (_key & 1u)*_place;
      _key >>= 1;
    }
    _place <<= 1;
  }
    
  unsigned _range = 1u << pbits;
  for (unsigned j=0; j<ndim; j++) {
    unsigned k=idim[j];
    XYZ[k] = xyz[k] + static_cast<double>(_range)/(1u<<level);
    xyz[k] = xyz[k]/_range * pH2OT::sides[k] - pH2OT::offset[k];
    XYZ[k] = XYZ[k]/_range * pH2OT::sides[k] - pH2OT::offset[k];
  }

  for (unsigned k=0; k<3; k++) {
    bool nope = true;
    for (unsigned j=0; j<ndim; j++)
      if (idim[j] == k) nope = false;
    if (nope) {
      xyz[k] = -pH2OT::offset[k];
      XYZ[k] = pH2OT::sides[k] - pH2OT::offset[k];
    }
  }

}

vector<double> tCell::getCorner()
{
  vector<double> xyz(3, 0);

  unsigned _place = 1u;
  unsigned _key   = mykey;

  for (unsigned i=0; i<level; i++) {
    for (unsigned j=0; j<ndim; j++) {
      unsigned k=idim[j];
      xyz[2-k] += (_key & 1u)*_place;
      _key >>= 1;
    }
    _place <<= 1;
  }

  unsigned _range = 1u << pbits;
  for (unsigned j=0; j<ndim; j++) {
    unsigned k=idim[j];
    xyz[k] = xyz[k]/_range * pH2OT::sides[k] - pH2OT::offset[k];
  }
  for (unsigned k=0; k<3; k++) {
    bool nope = true;
    for (unsigned j=0; j<ndim; j++)
      if (idim[j] == k) nope = false;
    if (nope) {
      xyz[k] = -pH2OT::offset[k];
    }
  }

  return xyz;
}



sCell* tCell::findSampleCell(unsigned Bucket)
{
  tCell *cur = this;		// Begin with this cell
  while(cur->count < Bucket) {
				// We are at the root
    if (cur->parent == 0) break;
				// Keep walking up the tree . . 
    cur = cur->parent;
  }
				// The answer.
  return cur;
}


double tCell::Volume()
{
  return tr->ph->TotalVolume()/(1u << ndim*level);
}

double tCell::Scale()
{
  return 1.0/(1u << level);
}


void tCell::Accum(const double* s)
{
  for (int i=0; i<10; i++) state[i] += s[i];
  if (parent) parent->Accum(s);
}

tTree::tTree(pH2OT *p, unsigned nb)
{
  live++;

  ph       = p;
  maxlevel = nb;

  // Create the new root
  root = new tCell(this, nb);
}

void tTree::reset(pH2OT *p, unsigned nb)
{
  delete root;

  ph       = p;
  maxlevel = nb;

  // Create the new root
  root = new tCell(this, nb);
}

void tTree::ZeroAccum(tCell *t)
{
  for (int i=0; i<10; i++) t->state[i] = 0.0;
  map<unsigned, tCell*>::iterator c;
  for (c=t->children.begin(); c!=t->children.end(); c++) 
    ZeroAccum(c->second);
}

tCell* tTree::MakeCell(unsigned key)
{
  return TreeWalk(root, key);
}

tCell* tTree::TreeWalk(tCell *t, unsigned key)
{
  unsigned id = t->childId(key);

  if (t->children.find(id) == t->children.end())
    t->children[id] = new tCell(t, id);

  if (t->level==maxlevel-1) {
    // Sanity check
    if (t->children[id]->mykey != key) {
      cout << "TreeWalk: cell mismatch key=" << key
	   << " not equal mykey=" << t->children[id]->mykey << endl;
    }
    if (t->children[id]->ptree == 0)
      t->children[id]->ptree = new pTree(ph, t->children[id], key);
    return t->children[id];
  } else {
    return TreeWalk(t->children[id], key);
  }
}

void tTree::AddState(const double* s)
{
  vector<double> p(3);
  if (s[0]>0.0) {
    for (int k=0; k<3; k++) p[k] = s[7+k]/s[0];
  
    tCell *c = root->findAddNode(&p[0]);
    if (c==0) {
      oob++;
    } else {
      c->Accum(s);
    }
  }
}

void tTree::printCell(const double *p)
{
  tCell *c = root->findNode(p);
  if (c==0) {
    cerr << "Point out of bounds" << endl;
  } else {
    vector<double> xyz, XYZ;
    c->getRange(xyz, XYZ);
    for (int k=0; k<3; k++)
      cout << "[" << setw(8) << xyz[k] << ", "
	   << setw(8) << XYZ[k] << "] ";
    cout << endl;
  }
}


void tTree::printTree(tCell *c)
{
  vector<double> r, R;

  c->getRange(r, R);

  cout << right
       << "  0x" << setw(8) << hex << setfill('0') << c->mykey 
       << "  0x" << setw(8) << hex << setfill('0') << c->mask 
       << setw(3) << dec << setfill(' ') << c->level;
  for (int k=0; k<3; k++) 
    cout << setw(12) << r[k] << setw(12) << R[k];
  cout << endl;

  map<unsigned, tCell*>::iterator it;
  for (it=c->children.begin(); it!=c->children.end(); it++) 
    printTree(it->second);
}
