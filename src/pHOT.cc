
#include <config_exp.h>

#ifdef HAVE_OMP_H
#include <omp.h>
#endif

#include <cstdlib>
#include <cstdio>
#include <ctime>
#include <bitset>
#include <limits>

#include <iostream>
#include <iomanip>
#include <fstream>
#include <sstream>
#include <vector>
#include <list>
#include <map>
#include <set>
#include <algorithm>
#include <cmath>

#include <memory>

// Hardwired DEBUGGING statements for pHOT
// [set all to false for production]

// Chatty debugging statements with minimal CPU overhead
static bool DEBUG_NOISY   = false;

// Debugging that checks internal lists
static bool DEBUG_CHECK   = false;

// Extra-verbose output with internal checking
static bool DEBUG_EXTRA   = false;

// Extra debugging for the adjustTree algorithm
static bool DEBUG_ADJUST  = false;

// Clean up bad entries in particle list
static bool DEBUG_CLEAN   = false;

// Check for malformed particles in cells
static bool DEBUG_SANITY  = false;

// Debug keys across all nodes (uses MPI calls)
static bool DEBUG_KEYS    = false;

// Warn when root is on frontier
static bool DEBUG_ROOT    = false;

// Execute VERY SLOW body count report
static bool DEBUG_BCOUNT  = false;

// Execute body count summary report
static bool BC_SUMMARY    = true;

// Execute body count per process report
static bool BC_PROCESS    = true;

// Debug key list
static bool DEBUG_KEYLIST = false;

// Debug OOB diagnostic
static bool DEBUG_OOB     = true;

// Test the parallel merge algorithm against round robin
static bool DEBUG_MERGE   = false;

// Use round-robin distribution and scalar sort rather than parallel sort
static bool USE_RROBIN    = false;

#ifdef USE_GPTL
#include <gptl.h>
#endif

using namespace std;

#include "global.H"
#include "pHOT.H"

// Default side lengths for prism to be partitioned
std::vector<double> pHOT::sides(3, 2.0);

// Location of origin
std::vector<double> pHOT::offst(3, 1.0);

// Default quantiles for diagnostic
unsigned pHOT::ntile;
std::vector<unsigned> pHOT::qtile;

// Use computed effort per particle in domain decomposition (default: true)
bool pHOT::use_weight  = false;

double   pHOT::hystrs  = 0.5;

void pHOT::qtile_initialize()
{
  if (qtile.size()) return;

  qtile.push_back(10);
  qtile.push_back(50);
  qtile.push_back(90);

  ntile = qtile.size();
}

static bool wghtKEY(const pair<key_type, double>& a, 
		    const pair<key_type, double>& b)
{ 
  return (a.first < b.first); 
}

static bool wghtDBL(const pair<key_type, double>& a, 
		    const pair<key_type, double>& b)
{ 
  return (a.second < b.second); 
}

//
// Check sample cell sanity (for debugging)
// [set to false for production]
//
bool pHOT::samp_debug = true;

//
// Turn on/off subsampling the key list for partitioning
//
bool pHOT::sub_sample = false;

// For formatting
unsigned pHOT::klen = 3*nbits/4+6;

// For diagnostic MPI communication
MPI_Datatype pHOT::CellDiagType;

template<class U, class V>
struct pair_compare
{
  bool operator()(const pair<U, V>& a, const pair<U, V>& b)
  { return (a.first<=b.first); }
};

// Error messages
//
void pHOT::bomb(const string& membername, const string& msg)
{
  std::ostringstream sout;
  sout << "pHOT::" << membername << "(): " << msg << endl;
  throw GenericError(sout.str(), __FILE__, __LINE__, 1037, true);
}

/*
  Constructor: initialize domain
*/
pHOT::pHOT(Component *C, sKeySet spec_list)
{
  qtile_initialize();		// Quantile set up

#ifdef HAVE_OMP_H
  omp_set_num_threads(nthrds);	// OpenMP set up
#endif

  cc = C;			// Register the calling component

				// Register multiple species
  this->spec_list = spec_list;
  if (spec_list.size() == 0) spec_list.insert(Particle::defaultKey);
			
  partType = Hilbert;		// Partition type

				// Sanity check
  if (nbits*3 >= sizeof(key_type)*8) {
    unsigned maxb = sizeof(key_type)*8/3;
    if (maxb*3 >= sizeof(key_type)*8) maxb--;
    ostringstream mesg;
    mesg << "nbits=" << nbits << " but must be less than " << maxb;
    bomb("pHOT::pHOT", mesg.str());
  }
  klen = 3*nbits/4+6;

  volume = sides[0]*sides[1]*sides[2];	// Total volume of oct-tree region
  root = 0;

  offset = vector<double>(3);
  for (unsigned k=0; k<3; k++) offset[k] = offst[k];

  kbeg = vector<key_type>(numprocs);
  kfin = vector<key_type>(numprocs);

  key_min = key_type(1u) << (nbits*3);
  key_max = key_type(1u) << (nbits*3+1);

  m_xchange = n_xchange = 0;

  sumstep = sumzero = 0;

  cntr_total = cntr_new_key = cntr_mine = cntr_not_mine = cntr_ship = 0;

  numkeys = 0;
 
  // Initialize timing structures
  //
  keymk3 = exchg3 = cnvrt3 = tovlp3 = prepr3 = updat3 =
    scatr3 = reprt3 = tadjt3 = celcl3 = keycm3 = keybd3 =
    keyst3 = keygn3 = wait03 = wait13 = wait23 = keync3 =
    keyoc3 = barri3 = diagd3 = vector<float>(numprocs);
  
  numk3    = vector<unsigned>(numprocs);
  numfront = vector<int>(numprocs);
  displace = vector<int>(numprocs);

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

				// Initialize particle ferry
  pf = ParticleFerryPtr(new ParticleFerry(cc->niattrib, cc->ndattrib));

  				// Filename for debugging info
  debugf = outdir + runtag + ".pHOT_debug";
}

pHOT::~pHOT()
{
  delete root;
}

#undef NEWKEY

uint64_t split3( unsigned int a )
{
  // we only use the first 21 bits
  uint64_t x = a & 0x1fffff;
  x = (x | x << 32) & 0x001f00000000ffff; // shift left 32 bits, OR with self, and 00011111000000000000000000000000000000001111111111111111
  x = (x | x << 16) & 0x001f0000ff0000ff; // shift left 32 bits, OR with self, and 00011111000000000000000011111111000000000000000011111111
  x = (x | x << 8 ) & 0x100f00f00f00f00f; // shift left 32 bits, OR with self, and 0001000000001111000000001111000000001111000000001111000000000000
  x = (x | x << 4 ) & 0x10c30c30c30c30c3; // shift left 32 bits, OR with self, and 0001000011000011000011000011000011000011000011000011000100000000
  x = (x | x << 2 ) & 0x1249249249249249;
  return x;
}

uint64_t mortonEncode_mask( unsigned int* u )
{
  uint64_t answer = 0;
  answer |= split3(u[0]) | split3(u[1]) << 1 | split3(u[2]) << 2;
  return answer;
}

uint64_t newKey(double *d)
{
  const unsigned int maxI = 0x1fffff;
  unsigned int u[3];

  for (int k=0; k<3; k++) u[k] = d[k] * maxI;
  
  return mortonEncode_mask(&u[0]);
}


key_type pHOT::getKey(double *p)
{
#ifdef USE_GPTL
  GPTLstart("pHOT::getKey");
#endif

  // Bad values
  //
  for (unsigned k=0; k<3; k++) {
    if (std::isnan(p[k]) || std::isinf(p[k])) {
      d_val++;
      timer_keybods.stop();
      return key_type(0u);
    }
  }

  // Out of bounds?
  //
  for (unsigned k=0; k<3; k++) { 
    if (fabs((p[k]+offset[k])/sides[k])> 1.0) {	
      if (DEBUG_NOISY) {
	cout << "Coordinate out of pbounds in pHOT::key: ";
	for (int l=0; l<3; l++) cout << setw(18) << p[l] ;
	cout << endl;
      }
#ifdef USE_GPTL
      GPTLstop("pHOT::getKey");
#endif
      return key_type(0u);
    }
  }

#ifdef NEWKEY

  double dd[3];
  const uint64_t lead = (1ull << nbits*3);
  for (unsigned k=0; k<3; k++) dd[2-k] = (p[k]+offset[k])/sides[k];
  key_type _key = newKey(dd) | lead;

#else

  const double factor = static_cast<double>(key_type(1u)<<nbits);

  // const unsigned int factor = (key_type(1u)<<nbits)-1;

  const key_type mask = 0x1u;

  vector<key_type> bins(3, 0u);

				// Reverse the order
  for (unsigned k=0; k<3; k++)
    bins[2-k] = key_type( floor( ((p[k]+offset[k])/sides[k])*factor ) );
  
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

#endif
  
#ifdef USE_GPTL
  GPTLstop("pHOT::getKey");
#endif

  return _key;
}

string pHOT::printKey(key_type p)
{
  ostringstream sout, sret;

  unsigned short cnt = 0;
  for (unsigned k=0; k<nbits; k++) {
    sout << ( (p & 1u) ? '1' : '0' );
    if (++cnt==3) {sout << '.'; cnt = 0;}
    p = p>>1;
  }

  string s = sout.str();	// Reverse the string
  for (unsigned k=0; k<s.size(); k++) sret << s[s.size()-1-k];

  return sret.str();
}


void pHOT::computeCellStates()
{

  // Each node start at root and walk down the tree to zero the counts
  //
  root->zeroState();

  // Make temporary to get around OpenMP limitations.  This can be
  // changed eventually when the GnuC suite supports STL iterators.
  //
  size_t iTmp = 0, cTmp = frontier.size();
  std::vector<pCell*> tmp(cTmp);
  for (auto it : frontier) tmp[iTmp++] = it.second;

  // March through the frontier and accumulate the counts
  //
#pragma omp parallel for default(none) private(iTmp) shared(tmp, cTmp)
  for (iTmp = 0; iTmp < cTmp; iTmp++) {
    tmp[iTmp]->accumState();
  }
  

  // March through the frontier to find the sample cells
  //
#pragma omp parallel for default(none) private(iTmp) shared(tmp, cTmp)
  for (iTmp = 0; iTmp < cTmp; iTmp++) {
    tmp[iTmp]->findSampleCell("Frontier scan");
  }

  // Sanity check
  //
  if (samp_debug) {
    timer_diagdbg.start();

    ostringstream sout;
    sout << "pHOT::computeCellStates, myid=" << std::setw(5) << myid
	 << ", time=" << std::setw(14) << tnow << " [this is a BAD ERROR]";
    checkSampleCells(sout.str());

    static unsigned msgcnt=0, maxcnt=10;
    if (msgcnt < maxcnt) {
      // Look for root key on the frontier
      key_cell::iterator it=frontier.find(root->mykey);
      if (DEBUG_ROOT && it != frontier.end()) {
	cout << "computeCellStates, root on frontier"
	     << ", T=" << tnow
	     << ", owner=" << it->second->owner
	     << ", level=" << it->second->level
	     << ", count=" << it->second->ctotal
	     << ", mass="  << it->second->stotal[0]
	     << ", root key="    << hex << it->second->mykey
	     << ", sample key="  << hex << it->second->mykey
	     << ", sample cell=" << hex << it->second->sample
	     << dec << endl;
	if (++msgcnt==maxcnt)
	  cout << "computeCellStates, suppressing non-fatal "
	       << "\"root on frontier\" messages after " << maxcnt
	       << " on Proc# " << myid << endl;
      }
    }

    timer_diagdbg.stop();
  }
}

void pHOT::logFrontierStats()
{
  timer_diagdbg.start();

  vector<unsigned> fstat1(nbits, 0), fstat(nbits);
  for (auto it : frontier) fstat1[it.second->level]++;
  
  MPI_Reduce(&fstat1[0], &fstat[0], nbits, MPI_UNSIGNED, MPI_SUM, 0,
	     MPI_COMM_WORLD);
  
  if (myid==0) {
    string filename = outdir + "frontier_debug." + runtag;
    ofstream out(filename.c_str(), ios::app);
    out << "# Time=" << tnow << endl;
    out << setw(6) << left << "Proc";
    for (unsigned j=0; j<nbits; j++)
      out << left << setw(8) << j+1;
    out << endl;
    out << setw(6) << left << "----";
    for (unsigned j=0; j<nbits; j++)
      out << left << setw(8) << "----";
    out << endl;
    out << setw(6) << left << "All";
    for (unsigned j=0; j<nbits; j++) 
      out << left << setw(8) << fstat[j];
    out << endl;
  }
  
  for (int n=0; n<numprocs; n++) {
    if (myid==n) {
      string filename = outdir + "frontier_debug." + runtag;
      ofstream out(filename.c_str(), ios::app);
      out << endl << left << setw(6) << myid;
      for (unsigned j=0; j<nbits; j++) 
	out << left << setw(8) << fstat1[j];
      out << endl;
    }
    MPI_Barrier(MPI_COMM_WORLD);
  }

  timer_diagdbg.stop();
}

key_type pHOT::getHeadKey()
{
  key_type headKey = 0ul;

  if (keybods.size()) {

    if (DEBUG_CHECK) {
      timer_diagdbg.start();
				// check validity of key
      if (bodycell.find(keybods.begin()->first) == bodycell.end()) {
	cout << "pHOT::getHeadKey, process " << myid << ": bad key=" 
	     << hex << keybods.begin()->first << dec
	     << " #cells=" << bodycell.size() << endl;
      }
      timer_diagdbg.stop();
    }
    
    headKey = bodycell.find(keybods.begin()->first)->second.first;
				// Number of bodies in my head cell
    if (DEBUG_CHECK) {
      timer_diagdbg.start();
      // Debug: check for key in frontier
      if (frontier.find(headKey) == frontier.end()) {
	cout << "pHOT::getHeadKey, process " << myid << ": headKey=" 
	     << headKey << dec << " is NOT on frontier with frontier size="
	       << frontier.size() << " [1]" << endl;
	std::cout << std::string(45, '-')    << std::endl
		  << "--- Bodycell list ---" << std::endl
		  << std::string(45, '-')    << std::endl;
	for (auto k : bodycell) {
	  std::cout << std::setw(15) << k.first
		    << std::setw(15) << k.second.first
		    << std::setw(15) << k.second.second
		    << std::endl;
	}
      }
      timer_diagdbg.stop();
    }
  }

  return headKey;
}

key_type pHOT::getTailKey()
{
  key_type tailKey = 0ul;
				// Compute the tailkey
  if (keybods.size()) {
    if (DEBUG_CHECK) {
      timer_diagdbg.start();
				// check validity of key
      if (bodycell.find(keybods.rbegin()->first) == bodycell.end()) {
	cout << "pHOT::getTailKey, process " << myid << ": bad tail key=" 
	     << hex << keybods.rbegin()->first << dec
	     << " #cells=" << bodycell.size() << endl;
      }
      timer_diagdbg.stop();
    }
    
    tailKey = bodycell.find(keybods.rbegin()->first)->second.first;
    
    if (DEBUG_CHECK) {
      if (tailKey == 1ul) {
	if (frontier.find(tailKey) == frontier.end()) {
	  cout << "pHOT::getTailKey, process " << myid << ": tailKey=" 
	       << tailKey << dec << " not on frontier! [3]" << endl;
	} else {
	  // Verbose info
	  if (false) {
	    cout << "pHOT::getTailKey, process " << myid << ": tailKey=" 
		 << tailKey << dec << " IS on frontier with frontier size="
		 << frontier.size() << " [3]" << endl;
	    cout << std::string(30, '-')    << std::endl << std::setfill('-')
		 << "--- Frontier list ---" << std::endl << std::setfill(' ')
		 << std::string(30, '-')    << std::endl << std::hex;
	    for (key_cell::iterator 
		   kit=frontier.begin(); kit!=frontier.end(); kit++) 
	      {
		std::cout << std::setw(15) << kit->first
			  << std::endl;
	      }
	    std::cout << std::string(30, '-')    << std::endl << std::dec;
	  }
	}
      }
    }

    if (DEBUG_CHECK) {
      timer_diagdbg.start();
				// Debug: check for tail key in frontier
      if (frontier.find(tailKey) == frontier.end()) {
	cout << "pHOT::getTailKey, process " << myid << ": tailKey=" 
	     << tailKey << dec << " is NOT on frontier with frontier size="
	     << frontier.size() << " [1]" << endl;
	std::cout << std::string(45, '-')    << std::endl
		  << "--- Bodycell list ---" << std::endl
		  << std::string(45, '-')    << std::endl;
	for (auto k : bodycell) {
	  std::cout << std::setw(15) << k.first
		    << std::setw(15) << k.second.first
		    << std::setw(15) << k.second.second
		    << std::endl;
	}
      }
      timer_diagdbg.stop();
    }
  }

  return tailKey;
}


void pHOT::makeTree()
{
  (*barrier)("pHOT::entering makeTree", __FILE__, __LINE__);

#ifdef USE_GPTL
  GPTLstart("pHOT::makeTree");
#endif
  //
  // Clean up
  // 
  frontier.clear();
  bodycell.clear();
  adjcnt = 0;

  delete root;

  if (DEBUG_EXTRA) {
    timer_diagdbg.start();

    string sname =  runtag + ".pHOT_storage";
    for (int n=0; n<numprocs; n++) {
      if (myid==n) {
	ofstream out(sname.c_str(), ios::app);
	if (out) {
	  out << setw(18) << tnow
	      << setw(6)  << myid
	      << setw(12) << keybods.size()
	      << setw(12) << frontier.size()
	      << setw(12) << bodycell.size()
	      << endl;
	  if (myid==numprocs-1) out << endl;
	  out.close();
	}
      }
      (*barrier)("pHOT::makeTree in debug pHOT storage", __FILE__, __LINE__);
      MPI_Barrier(MPI_COMM_WORLD);
    }
    timer_diagdbg.stop();
  }

  //
  // Make the root
  //
  root = new pCell(this);

  //
  // Assign the origin offset in the rectangular prism
  //
  for (unsigned k=0; k<3; k++) offset[k] = offst[k];
  MPI_Bcast(&offset[0], 3, MPI_DOUBLE, 0, MPI_COMM_WORLD);

  //
  // Add the data
  //
  pCell* p = root;
  for (auto it : keybods) {

    if (DEBUG_KEYS and (it.first < key_min || it.first >= key_max) ) {
      double* pos = cc->Particles()[it.second]->pos;
      key_type tkey = getKey(pos);
      cout << "Process " << myid << ": in makeTree, key=" 
	   << hex << it.first << ", tkey=" << tkey
	   << "[" << key_min << ", " << key_max << "]"
	   << endl << dec
	   << " pos = [" << pos[0] << ", " << pos[1] << ", "
	   << pos[2] << "]" << endl;
    }

    p = p->Add(it);		// Do the work
  }

				// Sanity checks and debugging
  if (DEBUG_BCOUNT) {
				// Report on particle sizes for each node
    if (BC_SUMMARY) {
      unsigned long bdcel1=cc->Particles().size(), bdcelmin, bdcelmax, bdcelsum, bdcelsm2;

      (*barrier)("pHOT::makeTree initial body report", __FILE__, __LINE__);

      MPI_Reduce(&bdcel1, &bdcelmin, 1, MPI_UNSIGNED_LONG, MPI_MIN, 0, MPI_COMM_WORLD);
      MPI_Reduce(&bdcel1, &bdcelmax, 1, MPI_UNSIGNED_LONG, MPI_MAX, 0, MPI_COMM_WORLD);
      MPI_Reduce(&bdcel1, &bdcelsum, 1, MPI_UNSIGNED_LONG, MPI_SUM, 0, MPI_COMM_WORLD);
      bdcel1 = bdcel1*bdcel1;
      MPI_Reduce(&bdcel1, &bdcelsm2, 1, MPI_UNSIGNED_LONG, MPI_SUM, 0, MPI_COMM_WORLD);

      if (myid==0)
	cout << endl << "In makeTree, particles"
	     << " min="  << bdcelmin 
	     << " max="  << bdcelmax
	     << " mean=" << bdcelsum/numprocs
	     << " stdv=" << sqrt( (bdcelsm2 - bdcelsum*bdcelsum/numprocs)/(numprocs-1) )
	     << " time=" << tnow
	     << endl;
    }
				// Report on body counts for each node
    if (BC_PROCESS) {
      const int w1 = 5, w2 = 10, w3 = 14;
      for (int i=0; i<numprocs; i++) {
	if (myid==0 && i==0) {
	  cout << right << setfill('-') 
	       << setw(w1) << '+'
	       << setw(w2) << '+' << setw(w2) << '+' 
	       << setw(w3) << '+' << setw(w3) << '+' 
	       << setw(w3) << '+' << setfill(' ') << endl
	       << setw(w1) << "pid"
	       << setw(w2) << "# bodies"
	       << setw(w2) << "# cells"
	       << setw(w3) << "x"
	       << setw(w3) << "y"
	       << setw(w3) << "z"
	       << endl << setfill('-') 
	       << setw(w1) << '+'
	       << setw(w2) << '+' << setw(w2) << '+' 
	       << setw(w3) << '+' << setw(w3) << '+' 
	       << setw(w3) << '+' << setfill(' ') << endl;
	}
	if (myid==i)  {
	  std::vector<double> pos(3, 0);
	  unsigned cnt = 0;
	  for (PartMapItr 
		 it=cc->Particles().begin(); it!=cc->Particles().end(); it++) 
	    {
	      for (int k=0; k<3; k++) pos[k] += it->second->pos[k];
	      cnt++;
	    }
	  cout << right << setw(w1)  << i 
	       << setw(w2) << cnt
	       << setw(w2) << frontier.size();
	  if (cnt)
	    cout << setw(w3) << pos[0]/cnt
		 << setw(w3) << pos[1]/cnt
		 << setw(w3) << pos[2]/cnt;
	  cout << endl;
	}
	(*barrier)("pHOT::makeTree body count report", __FILE__, __LINE__);
      }
      if (myid==0) {
	cout << right << setfill('-') 
	     << setw(w1) << '+'
	     << setw(w2) << '+' << setw(w2) << '+' 
	     << setw(w3) << '+' << setw(w3) << '+' 
	     << setw(w3) << '+' << setfill(' ') << endl;
      }
    }

    if (false && bodycell.size()==0) {
      cout << "pHOT::makeTree, process " << myid 
	   << ": unusual condition #bodycell=0"
	   << " with #keybods=" << keybods.size() 
	   << " and #bodies=" << cc->Particles().size()
	   << endl;
    }
  }
  
  //
  // Adjust boundaries bodies to prevent cell duplication on the boundary
  //

#ifdef USE_GPTL
  GPTLstart("pHOT::makeTree::adjustBoundaries");
#endif

				// Exchange boundary keys
  key_type headKey=0u, tailKey=0u, prevKey=0u, nextKey=0u;
  unsigned head_num=0, tail_num=0, next_num=0, prev_num=0;

				// Do the boundaries sequentially to prevent
				// inconstencies

  for (int n=1; n<numprocs; n++) {
				// Send the next node my tail value
				// to compare with its head
    if (myid==n-1) {

      tailKey = getTailKey();
      if (tailKey>0ul) tail_num = frontier[tailKey]->bods.size();

      MPI_Send(&tailKey,  1, MPI_EXP_KEYTYPE, n, 1000, MPI_COMM_WORLD);
      MPI_Send(&tail_num, 1, MPI_UNSIGNED,    n, 1001, MPI_COMM_WORLD);

      MPI_Recv(&nextKey, 1, MPI_EXP_KEYTYPE,  n, 1002, MPI_COMM_WORLD,
	       MPI_STATUS_IGNORE);
      MPI_Recv(&next_num, 1, MPI_UNSIGNED,    n, 1003, MPI_COMM_WORLD,
	       MPI_STATUS_IGNORE);

      if (tailKey != 0u && tailKey == nextKey) {
	if (tail_num <= next_num) {
	  if (tail_num) {
	    if (frontier.find(tailKey) == frontier.end()) {
	      std::cout << "pHOT::makeTree, process " << myid << ": tailKey=" 
			<< hex << tailKey
			<< dec << " is NOT on frontier with frontier size="
			<< frontier.size() << " [2]" << endl;
	      std::cout << std::string(30, '-')    << std::endl << std::setfill('-')
			<< "--- Frontier list ---" << std::endl << std::setfill(' ')
			<< std::string(30, '-')    << std::endl << std::hex;
	      for (auto k : frontier) {
		std::cout << std::setw(15) << k.first
			  << std::endl;
	      }
	      std::cout << std::string(45, '-')    << std::endl
			<< "--- Bodycell list ---" << std::endl
			<< std::string(45, '-')    << std::endl;
	      for (auto k : bodycell) {
		std::cout << std::setw(15) << k.first
			  << std::setw(15) << k.second.first
			  << std::setw(15) << k.second.second
			  << std::endl;
	      }
	      std::cout << std::string(45, '-')    << std::endl << std::dec;
	    }
	    sendCell(tailKey, n, tail_num);
	  } else
	    cout << "pHOT::makeTree, process " << myid 
		 << ": not sending cell with zero particles" << endl;
	} else {
	  if (next_num)
	    recvCell(n, next_num);
	  else
	    cout << "Process " << myid << ": not receiving cell with zero particles" << endl;
	}
      }

    }
				// Send the previous node my head value
				// to compare with its tail
    if (myid==n) {

      headKey = getHeadKey();
      if (headKey>0ul) head_num = frontier[headKey]->bods.size();

      MPI_Send(&headKey,  1, MPI_EXP_KEYTYPE, n-1, 1002, MPI_COMM_WORLD);
      MPI_Send(&head_num, 1, MPI_UNSIGNED,    n-1, 1003, MPI_COMM_WORLD);

      MPI_Recv(&prevKey,  1, MPI_EXP_KEYTYPE, n-1, 1000, MPI_COMM_WORLD,
	       MPI_STATUS_IGNORE);
      MPI_Recv(&prev_num, 1, MPI_UNSIGNED,    n-1, 1001, MPI_COMM_WORLD,
	       MPI_STATUS_IGNORE);

      if (headKey != 0u && headKey == prevKey) {
	if (head_num < prev_num) {
	  if (head_num) {
	    if (frontier.find(headKey) == frontier.end())
	      cout << "pHOT::makeTree, process " << myid << ": headKey=" 
		   << headKey << dec << " not on frontier! [2]" << endl;
	    sendCell(headKey, n-1, head_num);
	  } else
	    cout << "pHOT::makeTree, process " << myid 
		 << ": not sending cell with zero particles" << endl;
	} else {
	  if (prev_num)
	    recvCell(n-1, prev_num);
	  else
	    cout << "pHOT::makeTree, process " << myid 
		 << ": not receiving cell with zero particles" << endl;
	}
      }

    }    

  }

  (*barrier)("pHOT::makeTree(): boundaries adjusted", __FILE__, __LINE__);

#ifdef USE_GPTL
  GPTLstop ("pHOT::makeTree::adjustBoundaries");
  GPTLstart("pHOT::makeTree::getFrontier");
#endif

  // Compute the physical states in each cell for the entire tree and
  // find the sample cells
  //
  computeCellStates();

  // Get the true partition
  key_type kbeg1 = 0xffffffffffffffff, kfin1 = 0ul;

  for (auto i : keybods) {
    kbeg1 = min<key_type>(kbeg1, i.first);
    kfin1 = max<key_type>(kfin1, i.first);
  }

  MPI_Allgather(&kbeg1, 1, MPI_EXP_KEYTYPE, &kbeg[0], 1, MPI_EXP_KEYTYPE,
		MPI_COMM_WORLD);

  MPI_Allgather(&kfin1, 1, MPI_EXP_KEYTYPE, &kfin[0], 1, MPI_EXP_KEYTYPE,
		MPI_COMM_WORLD);

  unsigned isiz = loclist.size();
  loclist = kbeg;
  loclist.push_back(kfin[numprocs-1]);	// End point for binary search

  // Find min and max cell occupation
  //
  unsigned min1=std::numeric_limits<int>::max(), max1=0, nt;
  for (auto i : frontier) {
    nt = i.second->bods.size();
    min1 = min<unsigned>(min1, nt);
    max1 = max<unsigned>(max1, nt);
  }

  MPI_Allreduce(&min1, &min_cell, 1, MPI_UNSIGNED, MPI_MIN, MPI_COMM_WORLD);
  MPI_Allreduce(&max1, &max_cell, 1, MPI_UNSIGNED, MPI_MAX, MPI_COMM_WORLD);

  vector<unsigned>  chist1(max_cell-min_cell+1, 0);
  chist = vector<unsigned>(max_cell-min_cell+1, 0);
  for (auto i : frontier)
    chist1[i.second->bods.size()-min_cell]++;
  
  MPI_Allreduce(&chist1[0], &chist[0], max_cell-min_cell+1, 
		MPI_UNSIGNED, MPI_SUM, MPI_COMM_WORLD);

  for (unsigned i=1; i<max_cell-min_cell+1; i++) chist[i] += chist[i-1];

  // Accumulate the total number of cells in the tree
  //
  unsigned my_cells = frontier.size();
  MPI_Allreduce(&my_cells, &total_cells, 1, MPI_UNSIGNED, MPI_SUM, MPI_COMM_WORLD);

  // Check the load balance based on the total effort
  checkEffort(-1);

  if (DEBUG_CHECK) {
    checkIndices();
    std::ostringstream sout;
    sout << "pHOT::makeTree after cell creation, Node " 
	 << myid << ",  at T=" << tnow;
    checkKeybods(sout.str());
  }

#ifdef USE_GPTL
  GPTLstop("pHOT::makeTree::getFrontier");
  GPTLstop("pHOT::makeTree");
#endif

  //
  // Chatty frontier list output
  //
  if (DEBUG_NOISY) {
				// Find max and min key values
    key_type k_max =  0ul;
    key_type k_min = ~0ul;
    for (auto i : frontier) {
      k_min = std::min<key_type>(k_min, i.first);
      k_max = std::max<key_type>(k_max, i.first);
    }

				// Field sizes and headers
    const int nn = 6, nf = 20;
    const int tt = nn + 3*nf;
    if (myid==0) {
      std::cout << std::string(tt, '-') << std::endl 
		<< std::setfill('-')    << std::left 
		<< std::setw(tt) << "--- Frontier keys " << std::endl
		<< std::string(tt, '-') << std::endl
		<< std::right    << std::setfill(' ')
		<< std::setw(nn)  << "Node"
		<< std::setw(nf) << "Min key"
		<< std::setw(nf) << "Max key"
		<< std::setw(nf) << "Number"
		<< std::endl     << std::setfill('-')
		<< std::setw(nn) << '+'
		<< std::setw(nf) << '+'
		<< std::setw(nf) << '+'
		<< std::setw(nf) << '+'
		<< std::endl     << std::setfill(' ');
    }

    for (int id=0; id<numprocs; id++) {
      if (id==myid) {
	std::cout << std::right
		  << std::setw(nn) << id
		  << std::setw(nf) << std::hex << k_min
		  << std::setw(nf) << std::hex << k_max
		  << std::setw(nf) << std::dec << frontier.size()
		  << std::endl;
      }
      (*barrier)("pHOT::makeTree(): frontier info", __FILE__, __LINE__);
    }
  
    if (myid==0)
      std::cout << std::setfill('-') << std::right
		<< std::setw(nn) << '+'
		<< std::setw(nf) << '+'
		<< std::setw(nf) << '+'
		<< std::setw(nf) << '+'
		<< std::endl     << std::setfill(' ');
  }

  if (DEBUG_CHECK) {
    std::ostringstream sout;
    sout << "pHOT::makeTree at END, myid=" << std::setw(5) << myid 
	 << " at T=" << tnow;
    checkKeybodsFrontier(sout.str());
  }
}

unsigned pHOT::CellCount(double pctl)
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

void pHOT::densEmit(unsigned lev, pCell *p)
{
  if (p->level == lev) {
    cntlev[lev]++;
    if (p->parent) kidlev[lev] += p->parent->children.size();
    maslev[lev] += p->stotal[0];
    vollev[lev] += volume/static_cast<double>(key_type(1u) << (3*p->level));
  } else {
    for (auto i : p->children)
      densEmit(lev, i.second);
  }
}

void pHOT::densCheck()
{
  timer_diagdbg.start();

  makeState();

  unsigned MaxLev = 6;
  cntlev = vector<unsigned> (MaxLev+1, 0);
  kidlev = vector<unsigned> (MaxLev+1, 0);
  maslev = vector<double>   (MaxLev+1, 0);
  vollev = vector<double>   (MaxLev+1, 0);

  for (unsigned lev=0; lev<=MaxLev; lev++) {
    for (int n=0; n<numprocs; n++) {
      if (myid==n) densEmit(lev, root);	// Walk tree
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
    (*barrier)("pHOT: density check report", __FILE__, __LINE__);
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

  timer_diagdbg.stop();
}

void pHOT::dumpFrontier(std::ostream& out)
{
  unsigned sum = 0, cnt=0;
  double mean=0.0, disp=0.0;
  double totmass=0.0, totvol=0.0, tmp;

  timer_diagdbg.start();

  if (myid==0) {		// Write output header
    out << "#" << std::endl
	<< "# Frontier info" << std::endl
	<< "# id  "
	<< std::setw(12) << "key"
	<< std::setw( 8) << "level"
	<< std::setw(18) << "num"
	<< std::setw(18) << "mass"
	<< std::setw(18) << "density"
	<< std::setw(10) << "pos(x)"
	<< std::setw(10) << "var(x)"
	<< std::setw(10) << "min(x)"
	<< std::setw(10) << "max(x)"
	<< std::setw(10) << "pos(y)"
	<< std::setw(10) << "var(y)"
	<< std::setw(10) << "min(y)"
	<< std::setw(10) << "max(y)"
	<< std::setw(10) << "pos(z)"
	<< std::setw(10) << "var(z)"
	<< std::setw(10) << "min(z)"
	<< std::setw(10) << "max(z)"
	<< std::endl
	<< "# [1] "
	<< std::setw(12) << "[2]"
	<< std::setw( 8) << "[3]"
	<< std::setw(18) << "[4]"
	<< std::setw(18) << "[5]"
	<< std::setw(18) << "[6]"
	<< std::setw(10) << "[7]"
	<< std::setw(10) << "[8]"
	<< std::setw(10) << "[9]"
	<< std::setw(10) << "[10]"
	<< std::setw(10) << "[11]"
	<< std::setw(10) << "[12]"
	<< std::setw(10) << "[13]"
	<< std::setw(10) << "[14]"
	<< std::setw(10) << "[15]"
	<< std::setw(10) << "[16]"
	<< std::setw(10) << "[17]"
	<< std::setw(10) << "[18]"
	<< std::endl;
  }

  // I suppose I could do this with MPI_IO, but that seems like
  // overkill for gathering debug info . . .
  //
  for (int n=0; n<numprocs; n++) {

    // One line for each cell in the frontier
    std::vector<std::string> output;

    if (n==myid) {
      
      for (auto i : frontier) {
	// The output line for this cell
	std::ostringstream line;

	std::vector<double> pmin(3, DBL_MAX), pmax(3, -DBL_MAX);
	std::vector<double> mpos(3,0.0), vpos(3,0.0);
	unsigned num = i.second->bods.size();
	double  mass = 0.0;

	for (auto j : i.second->bods) {
	  mass += cc->particles[j]->mass;
	  for (unsigned k=0; k<3; k++) {
	    tmp = cc->particles[j]->pos[k];
	    mpos[k] += tmp;
	    vpos[k] += tmp*tmp;
	    pmin[k]  = std::min<double>(pmin[k], tmp);
	    pmax[k]  = std::max<double>(pmax[k], tmp);
	  }
	}
	
	totmass += mass;
	totvol  += volume/static_cast<double>(key_type(1u)<<(3*i.second->level));

	line << setw( 6) << myid
	     << setw(12) << hex << i.first << dec
	     << setw( 8) << i.second->level
	     << setw(18) << num
	     << setw(18) << mass
	     << setw(18) << mass/(volume/static_cast<double>(key_type(1u)<<(3*i.second->level)));
	
	for (unsigned k=0; k<3; k++) {
	  mpos[k] /= num;
	  if (num>1)
	    vpos[k] = sqrt( (vpos[k] - mpos[k]*mpos[k]*num)/(num-1) );
	  else
	    vpos[k] = 0.0;

	  line << setprecision(4) << setw(10) << mpos[k] 
	       << setprecision(4) << setw(10) << vpos[k]
	       << setprecision(4) << setw(10) << pmin[k]
	       << setprecision(4) << setw(10) << pmax[k];
	}
	mean += num;
	disp += num*num;
	sum  += num;
	cnt++;

	output.push_back(line.str());
      }

      if (n != 0) {
	unsigned osize = output.size();
	MPI_Send(&osize, 1, MPI_UNSIGNED, 0, 233, MPI_COMM_WORLD);
	for (auto s : output)
	  MPI_Send(s.c_str(), s.size(), MPI_CHAR, 0, 234, MPI_COMM_WORLD);
      }

    } // END: myid==n
    
    if (myid==0) {

      if (n != 0) {
	MPI_Status status;
	unsigned osize;
	int len;

	MPI_Recv(&osize, 1, MPI_UNSIGNED, n, 233, MPI_COMM_WORLD, MPI_STATUS_IGNORE);
	
	for (unsigned k=0; k<osize; k++) {
				// Get size of string
	  MPI_Probe(n, 234, MPI_COMM_WORLD, &status);
	  MPI_Get_count(&status, MPI_CHAR, &len);
				// Receive the string
	  std::shared_ptr<char> buf(new char[len]);
	  MPI_Recv(buf.get(), len, MPI_CHAR, n, 234, MPI_COMM_WORLD, &status);
				// Add to output list
	  output.push_back(std::string(buf.get(), len));
	}
      }

      for (auto s : output) out << s << std::endl;
				
    }
    (*barrier)("pHOT: dump frontier", __FILE__, __LINE__);
  }

  unsigned sum0=0, cnt0=0;
  double mean0=0.0, disp0=0.0, totmass0=0.0, totvol0=0.0;

  MPI_Reduce(&sum,     &sum0,     1, MPI_UNSIGNED, MPI_SUM, 0, MPI_COMM_WORLD);
  MPI_Reduce(&cnt,     &cnt0,     1, MPI_UNSIGNED, MPI_SUM, 0, MPI_COMM_WORLD);
  MPI_Reduce(&mean,    &mean0,    1, MPI_DOUBLE,   MPI_SUM, 0, MPI_COMM_WORLD);
  MPI_Reduce(&disp,    &disp0,    1, MPI_DOUBLE,   MPI_SUM, 0, MPI_COMM_WORLD);
  MPI_Reduce(&totmass, &totmass0, 1, MPI_DOUBLE,   MPI_SUM, 0, MPI_COMM_WORLD);
  MPI_Reduce(&totvol,  &totvol0,  1, MPI_DOUBLE,   MPI_SUM, 0, MPI_COMM_WORLD);

  if (myid==0) {
    out << "#" << endl
	<< "#" << setw(12) << "Total"  << setw(18) << sum0       << endl
	<< "#" << setw(12) << "Mass"   << setw(18) << totmass0   << endl
	<< "#" << setw(12) << "Volume" << setw(18) << totvol0    << endl
	<< "#" << setw(12) << "Mean"   << setw(18) << mean0/cnt0 << endl
	<< "#" << setw(12) << "Sigma"  << setw(18) 
	<< sqrt((disp0 - mean0*mean0/cnt0)/cnt0) << endl << "#" << endl;
  }

  timer_diagdbg.stop();
}

void pHOT::statFrontier()
{
  timer_diagdbg.start();

  unsigned sum1=0, cnt1=0, sum=0, cnt=0, num;
  double mean1=0.0, disp1=0.0, mean=0.0, disp=0.0;
  vector<unsigned> freq1(pCell::bucket+1, 0), freq(pCell::bucket+1, 0);
  
  for (auto i : frontier) {
    num = i.second->bods.size();
    mean1 += num;
    disp1 += num*num;
    sum1  += num;
    freq1[num]++;
    cnt1++;
  }

  MPI_Reduce(&sum1, &sum, 1, MPI_UNSIGNED, MPI_SUM, 0, MPI_COMM_WORLD );
  MPI_Reduce(&cnt1, &cnt, 1, MPI_UNSIGNED, MPI_SUM, 0, MPI_COMM_WORLD );
  MPI_Reduce(&mean1, &mean, 1, MPI_DOUBLE, MPI_SUM, 0, MPI_COMM_WORLD );
  MPI_Reduce(&disp1, &disp, 1, MPI_DOUBLE, MPI_SUM, 0, MPI_COMM_WORLD );
  MPI_Reduce(&freq1[0], &freq[0], pCell::bucket+1, 
	     MPI_UNSIGNED, MPI_SUM, 0, MPI_COMM_WORLD );

  if (myid==0) {
    cout << endl;
    cout << setw(12) << "Size" << setw(18) << frontier.size() << endl;
    cout << setw(12) << "Total" << setw(18) << sum << endl;
    cout << setw(12) << "Mean" << setw(18) << mean/cnt << endl;
    cout << setw(12) << "Sigma" << setw(18) 
	 << sqrt((disp - mean*mean/cnt)/cnt) << endl;
    
    cout << endl;
    for (unsigned n=0; n<pCell::bucket+1; n++)
      cout << setw(5) << n << setw(5) << freq[n] << endl;
    cout << endl;
  }

  timer_diagdbg.stop();
}


void pHOT::testFrontier(string& filename)
{
  pCell *p;

  timer_diagdbg.start();

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


  for (int n=0; n<numprocs; n++) {

    if (n == myid) {

      ofstream out(filename.c_str(), ios::app);

      for (auto c : frontier) {
	double mass=0, temp=0, pos[]={0,0,0}, vel[]={0,0,0};
	p = c.second;
	vector<unsigned long>::iterator ib = p->bods.begin();
	while (ib != p->bods.end()) {
	  mass += cc->particles[*ib]->mass;
	  for (int k=0; k<3; k++) {
	    pos[k] += 
	      cc->particles[*ib]->mass * 
	      cc->particles[*ib]->pos[k];
	    vel[k] += 
	      cc->particles[*ib]->mass * 
	      cc->particles[*ib]->vel[k];
	    temp += 
	      cc->particles[*ib]->mass * 
	      cc->particles[*ib]->vel[k] *
	      cc->particles[*ib]->vel[k];
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

	double vol = volume/static_cast<double>( key_type(1u) << (3*p->level));

	out << " ";
	out << setw(prec[n]) << setiosflags(fmt[n]|typ[n]) << tnow;
	n++;
	out << setw(prec[n]) << setiosflags(fmt[n]|typ[n]) << myid;
	n++;
	out << setw(prec[n]) << setiosflags(fmt[n]|typ[n]) << hex << c.first << dec;
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

    (*barrier)("HOT: test frontier", __FILE__, __LINE__);
  }

  timer_diagdbg.stop();

}


void pHOT::countFrontier(vector<unsigned>& ncells, vector<unsigned>& bodies)
{
  pCell *p;
  map<unsigned, pair<unsigned, unsigned> > data;
  map<unsigned, pair<unsigned, unsigned> >::iterator d;
  unsigned maxLev1=0, maxLev=0;
  
  timer_diagdbg.start();

  for (auto i : frontier) {
    p = i.second;
    maxLev1 = max<unsigned>(maxLev1, p->level);
    if ((d=data.find(p->level)) == data.end()) {
      data[p->level] = pair<unsigned, unsigned>(1, p->bods.size());
    } else {
      d->second.first++;
      d->second.second += p->bods.size();
    }
  }

  MPI_Allreduce(&maxLev1, &maxLev, 1, MPI_UNSIGNED, MPI_MAX, 
		MPI_COMM_WORLD);

  vector<unsigned> tcellcnts(maxLev+1, 0), tbodscnts(maxLev+1, 0);
  if (myid==0) {
    ncells = vector<unsigned>(maxLev+1);
    bodies = vector<unsigned>(maxLev+1);
  }

  for (auto d : data) {
    tcellcnts[d.first] = d.second.first;
    tbodscnts[d.first] = d.second.second;
  }

  MPI_Reduce(&tcellcnts[0], &ncells[0], maxLev+1, MPI_UNSIGNED, MPI_SUM, 
	     0, MPI_COMM_WORLD);

  MPI_Reduce(&tbodscnts[0], &bodies[0], maxLev+1, MPI_UNSIGNED, MPI_SUM,
	     0, MPI_COMM_WORLD);

  timer_diagdbg.stop();
}


void pHOT::sendCell(key_type key, int to, unsigned num)
{
#ifdef USE_GPTL
  GPTLstart("pHOT::sendCell");
#endif

  key_cell::iterator kit = frontier.find(key);
  pCell *p = 0;

  if (kit == frontier.end()) {
    std::cout << "pHOT::sendCell: myid=" << myid
	      << ", key=" << key << " is NOT on frontier "
	      << " with frontier size=" << frontier.size() << std::endl;
    std::cout << std::string(30, '-')    << std::endl << std::setfill('-')
	      << "--- Frontier list ---" << std::endl << std::setfill(' ')
	      << std::string(30, '-')    << std::endl;
    for (auto k : frontier)
      std::cout << std::setw(10) << k.first
		<< std::endl;
    std::cout << std::string(30, '-')    << std::endl;
    num = 0;
  } else {
    p = frontier.find(key)->second;
  }
  
  if (DEBUG_NOISY) {
    std::cout << "Process " << std::left << std::setw(4) << myid 
	      << std::setw(12) << ": sending " << std::setw(4) << num 
	      << "  to  " << std::setw(4) << to << endl;
  }

  vector<unsigned long> erased;
  vector<double> buffer1(3*num);
  vector<unsigned> buffer2(num);
  vector<key_type> buffer3(num);

  if (DEBUG_SANITY) {
    unsigned crazy = 0;
    vector<unsigned long>::iterator ib = p->bods.begin();
    for (unsigned j=0; j<num; j++) {
      if (cc->particles[*(ib++)]->indx == 0) crazy++;
    }
    if (crazy) std::cout << "[sendCell node " << myid << " has " << crazy
			 << " crazy bodies out of " << num << "]" << std::endl;
  }

  pf->ShipParticles(to, myid, num);
  
  key_pair tpair;
  vector<unsigned long>::iterator ib = p->bods.begin();
  for (unsigned j=0; j<num; j++) {

    pf->SendParticle(cc->particles[*ib]);
    
				// Find the record and delete it
    tpair.first  = cc->particles[*ib]->key;
    tpair.second = cc->particles[*ib]->indx;
    key_indx::iterator it = keybods.find(tpair);

    if (it != keybods.end()) {
				// Remove the key from the cell list
      {
	key2Range ij = bodycell.equal_range(it->first);

	if (ij.first != ij.second) {
	  key_key::iterator ijk=ij.first, rmv;
	  while (ijk!=ij.second) {
	    if ((rmv=ijk++)->second.second == tpair.second) bodycell.erase(rmv);
	  }
	} else {
	  std::cout << "pHOT::sendCell, myid=" << setw(5) << myid
		    << ", body not in bodycell list" << std::endl;
	}
      }

      cc->particles.erase(*ib);
      if (DEBUG_CHECK) erased.push_back(*ib);

      {
	key_indx::iterator ij = keybods.find(*it);

	if (ij != keybods.end()) {
	  keybods.erase(ij);
	} else {
	  std::cout << "pHOT::sendCell, myid=" << setw(5) << myid
		    << ", body not in keybods list" << std::endl;
	}
      }


    } else {
      cerr << "pHOT::sendCell, myid=" << myid << ": error, "
	   << "removing body from keybods list that does not exist! " 
	   << endl;
    }
    ib++;
  }

  // If this cell is not the root
  //
  if (p->parent) {
    
    // Delete this cell from the parent
    p->parent->children.erase( (p->mykey & 0x7u) );
	
    if (DEBUG_CHECK) {
      timer_diagdbg.start();
      if (frontier.find(p->mykey)==frontier.end()) {
	cout << "Process " << myid << ": in pHOT:sendCell: "
	     << " key not on frontier as expected" << endl;
      }
      timer_diagdbg.stop();
    }

    // Delete this cell from the frontier
    frontier.erase(p->mykey);

    // Delete the cell altogether
    delete p;

  } else {			// Special treatment for root
    p->keys.clear();
    p->bods.clear();
  }

  if (DEBUG_CHECK) {
    timer_diagdbg.start();
    for (auto i : erased) {
      if (cc->particles.find(i) != cc->particles.end())
	cout << "pHOT::sendCell proc=" << myid
	     << " found erased index=" << i << endl;
    }
    timer_diagdbg.stop();
  }

  // Refresh size of local particle list
  cc->nbodies = cc->particles.size();

#ifdef USE_GPTL
  GPTLstop("pHOT::sendCell");
#endif
}


void pHOT::recvCell(int from, unsigned num)
{
#ifdef USE_GPTL
  GPTLstart("pHOT::recvCell");
#endif

  if (DEBUG_NOISY) 
    std::cout << "Process " << std::left << std::setw(4) << myid 
	      << std::setw(12) << ": receiving " << std::setw(4) << num 
	      << " from " << std::setw(4) << from << endl;

  pCell *p = root;

  pf->ShipParticles(myid, from, num);

  for (unsigned j=0; j<num; j++) {
    PartPtr part = pf->RecvParticle();
    if (part->indx==0 || part->mass<=0.0 || std::isnan(part->mass)) {
      cout << "[recvCell, myid=" << myid 
	   << ", will ignore crazy body with indx=" << part->indx 
	   << ", j=" << j << ", num=" << num << ", mass=" << part->mass
	   << ", key=" << hex << part->key << dec << "]"
	   << " from Node " << from << std::endl;
    } else {
      cc->particles[part->indx] = part;
      if (part->key == 0u) continue;
      if (part->key < key_min || part->key >= key_max) {
	cout << "Process " << myid << ": in recvCell, key=" 
	     << hex << part->key << dec << "]"
	  ;
      }
      key_pair tpair(part->key, part->indx);
      keybods.insert(tpair);
      p = p->Add(tpair);
    }
  }

  // Refresh size of local particle list
  cc->nbodies = cc->particles.size();

#ifdef USE_GPTL
  GPTLstop("pHOT::recvCell");
#endif
}

void pHOT::makeState()
{
  // Currently unused
}


void pHOT::State(double *x, double& dens, double& temp,
		 double& velx, double& vely, double& velz)
{
  key_type key = getKey(x);

  dens = temp = velx = vely = velz = 0.0;

  // Walk tree to get count
  //
  unsigned count, level;
  vector<double> state(7);
  root->Find(key, count, level, state);

  vector<double>   stt1(numprocs*7, 0);
  vector<double>   stt0(numprocs*7, 0);
  vector<unsigned> cnt1(numprocs, 0), lev1(numprocs, 0);
  vector<unsigned> cnt0(numprocs, 0), lev0(numprocs, 0);

  cnt1[myid] = count;
  lev1[myid] = level;
  for (int k=0; k<7; k++) stt1[myid*7+k] = state[k];


  MPI_Reduce(&cnt1[0], &cnt0[0], numprocs,   MPI_UNSIGNED, MPI_SUM, 
	     0, MPI_COMM_WORLD);

  MPI_Reduce(&lev1[0], &lev0[0], numprocs,   MPI_UNSIGNED, MPI_SUM, 
	     0, MPI_COMM_WORLD);

  MPI_Reduce(&lev1[0], &lev0[0], numprocs,   MPI_UNSIGNED, MPI_SUM, 
	     0, MPI_COMM_WORLD);

  MPI_Reduce(&stt1[0], &stt0[0], numprocs*5, MPI_DOUBLE,   MPI_SUM, 
	     0, MPI_COMM_WORLD);


  // Compute the state variables for the "deepest" cell(s)
  //
  if (myid==0) {
    vector< pair<unsigned, unsigned> > dlist;
    for (int n=0; n<numprocs; n++) 
      dlist.push_back(pair<unsigned, int>(lev0[n], n));

    sort(dlist.begin(), dlist.end(), greater< pair<unsigned, int> >());

				// This is the deepest level
    unsigned clv = dlist[0].first;
    unsigned cnt = 1;
    for (int k=0; k<7; k++) state[k] = stt0[7*dlist[0].second + k];

				// Add others at the same level
    for (int n=1; n<numprocs; n++) {
      if (clv != dlist[n].first) break;
      for (int k=0; k<7; k++) state[k] += stt0[7*dlist[n].second + k];
      cnt++;
    }

				// Mass average
    if (state[0]>0.0) {
      double disp = 0.0;
      for (int k=0; k<3; k++)
	disp += (state[1+k] - state[4+k]*state[4+k]/state[0])/state[0];
      
      dens = state[0] * static_cast<double>(key_type(1u) << (3*clv))/(volume*cnt);
      temp = 0.333333333333*disp;
      velx = state[4]/state[0];
      vely = state[5]/state[0];
      velz = state[6]/state[0];
    }
    else {
      dens = temp = velx = vely = velz = 0.0;
    }

  }

}


void pHOT::Slice(int nx, int ny, int nz, string cut, string prefix)
{
  vector<double> pt(3);
  double dens, temp, vx, vy, vz;

  double dx = sides[0]/nx;
  double dy = sides[1]/ny;
  double dz = sides[2]/nz;

  ofstream out;

  if (myid==0) out.open(string(prefix + "." + cut).c_str());
  
  makeState();

  if (cut.compare("XY")==0) {

				// X-Y slice
    pt[2] = 0.5*sides[2] - offset[2];
    for (int i=0; i<nx; i++) {
      pt[0] = (0.5+i)*dx;
      for (int j=0; j<ny; j++) {
	pt[1] = (0.5+j)*dy;
	
	State(&pt[0], dens, temp, vx, vy, vz);

	if (myid==0)
	  out << setw(15) << pt[0]
	      << setw(15) << pt[1]
	      << setw(15) << dens
	      << setw(15) << temp
	      << setw(15) << vx
	      << setw(15) << vy
	      << setw(15) << vz
	      << endl;
      }
      if (myid==0) out << endl;
    }
  }

  if (cut.compare("XZ")==0) {

				// X-Z slice
    pt[1] = 0.5*sides[1] - offset[1];
    for (int i=0; i<nx; i++) {
      pt[0] = (0.5+i)*dx;
      for (int j=0; j<nz; j++) {
	pt[2] = (0.5+j)*dz;
	
	State(&pt[0], dens, temp, vx, vy, vz);

	if (myid==0)
	  out << setw(15) << pt[0]
	      << setw(15) << pt[2]
	      << setw(15) << dens
	      << setw(15) << temp
	      << setw(15) << vx
	      << setw(15) << vy
	      << setw(15) << vz
	      << endl;
      }
      if (myid==0) out << endl;
    }
  }

  if (cut.compare("YZ")==0) {

				// Y-Z slice
    pt[0] = 0.5*sides[0] - offset[0];
    for (int i=0; i<ny; i++) {
      pt[1] = (0.5+i)*dy;
      for (int j=0; j<nz; j++) {
	pt[2] = (0.5+j)*dz;
	
	State(&pt[0], dens, temp, vx, vy, vz);
      
	if (myid==0)
	  out << setw(15) << pt[1]
	      << setw(15) << pt[2]
	      << setw(15) << dens
	      << setw(15) << temp
	      << setw(15) << vx
	      << setw(15) << vy
	      << setw(15) << vz
	      << endl;
      }
      if (myid==0) out << endl;
    }
  }

}

void pHOT::Slab
(vector<int>& n, vector<double>& pmin, vector<double>& pmax, string cut,
 vector<double>&    x, vector<double>& dens, vector<double>& temp, 
 vector<double>& velx, vector<double>& vely, vector<double>& velz)
{
  vector<double> pt(3);
  double Dens, Temp, Velx, Vely, Velz;

  double dx = (pmax[0] - pmin[0])/n[0];
  double dy = (pmax[1] - pmin[1])/n[1];
  double dz = (pmax[2] - pmin[2])/n[2];

  ofstream out;

  makeState();

  if (cut.compare("X")==0) {

    x    = vector<double>(n[0], 0.0);
    dens = vector<double>(n[0], 0.0);
    temp = vector<double>(n[0], 0.0);
    velx = vector<double>(n[0], 0.0);
    vely = vector<double>(n[0], 0.0);
    velz = vector<double>(n[0], 0.0);

				// X cut
    for (int i=0; i<n[0]; i++) {
      pt[0] = x[i] = pmin[0] + (0.5+i)*dx;

      for (int j=0; j<n[1]; j++) {
	pt[1] = pmin[1] + (0.5+j)*dy;
	
	for (int k=0; k<n[2]; k++) {
	  pt[2] = pmin[2] + (0.5+k)*dz;

	  State(&pt[0], Dens, Temp, Velx, Vely, Velz);

	  dens[i] += Dens;
	  temp[i] += Temp;
	  velx[i] += Velx;
	  vely[i] += Vely;
	  velz[i] += Velz;
	}
      }
      
      dens[i] /= n[1]*n[2];
      temp[i] /= n[1]*n[2];
      velx[i] /= n[1]*n[2];
      vely[i] /= n[1]*n[2];
      velz[i] /= n[1]*n[2];
    }
  }

  if (cut.compare("Y")==0) {

    x    = vector<double>(n[1], 0.0);
    dens = vector<double>(n[1], 0.0);
    temp = vector<double>(n[1], 0.0);
    velx = vector<double>(n[1], 0.0);
    vely = vector<double>(n[1], 0.0);
    velz = vector<double>(n[1], 0.0);

				// Y cut
    for (int j=0; j<n[1]; j++) {
      pt[1] = x[j] = pmin[1] + (0.5+j)*dy;
	
      for (int i=0; i<n[0]; i++) {
	pt[0] = pmin[0] + (0.5+i)*dx;

	for (int k=0; k<n[2]; k++) {
	  pt[2] = pmin[2] + (0.5+k)*dz;

	  State(&pt[0], Dens, Temp, Velx, Vely, Velz);

	  dens[j] += Dens;
	  temp[j] += Temp;
	  velx[j] += Velx;
	  vely[j] += Vely;
	  velz[j] += Velz;
	}
      }

      dens[j] /= n[0]*n[2];
      temp[j] /= n[0]*n[2];
      velx[j] /= n[0]*n[2];
      vely[j] /= n[0]*n[2];
      velz[j] /= n[0]*n[2];
    }
  }


  if (cut.compare("Z")==0) {

    x    = vector<double>(n[2], 0.0);
    dens = vector<double>(n[2], 0.0);
    temp = vector<double>(n[2], 0.0);
    velx = vector<double>(n[2], 0.0);
    vely = vector<double>(n[2], 0.0);
    velz = vector<double>(n[2], 0.0);

				// Z cut
    for (int k=0; k<n[2]; k++) {
      pt[2] = x[k] = pmin[2] + (0.5+k)*dz;

      for (int i=0; i<n[0]; i++) {
	pt[0] = pmin[0] + (0.5+i)*dx;

	for (int j=0; j<n[1]; j++) {
	  pt[1] = pmin[1] + (0.5+j)*dy;
	
	  State(&pt[0], Dens, Temp, Velx, Vely, Velz);

	  dens[k] += Dens;
	  temp[k] += Temp;
	  velx[k] += Velx;
	  vely[k] += Vely;
	  velz[k] += Velz;
	}
      }

      dens[k] /= n[1]*n[2];
      temp[k] /= n[1]*n[2];
      velx[k] /= n[1]*n[2];
      vely[k] /= n[1]*n[2];
      velz[k] /= n[1]*n[2];
    }
  }

}

double pHOT::minVol()
{
  unsigned MaxLev = 0;
  for (auto i : frontier)
    MaxLev = max<unsigned>(MaxLev, i.second->level);

  double vol1, vol;
  vol1 = volume/static_cast<double>(key_type(1u) << (3*MaxLev));
  MPI_Allreduce(&vol1, &vol, 1, MPI_DOUBLE, MPI_MIN, MPI_COMM_WORLD);

  return vol;
}

double pHOT::maxVol()
{
  unsigned MinLev = std::numeric_limits<int>::max();
  for (auto i : frontier)
    MinLev = min<unsigned>(MinLev, i.second->level);

  double vol1, vol;
  vol1 = volume/static_cast<double>(key_type(1u) << (3*MinLev));
  MPI_Allreduce(&vol1, &vol, 1, MPI_DOUBLE, MPI_MAX, MPI_COMM_WORLD);

  return vol;
}

double pHOT::medianVol()
{
  unsigned mlev, num;
  vector<unsigned> lev;

  for (auto i : frontier) 
    lev.push_back(i.second->level);

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

  return volume/static_cast<double>(key_type(1u) << (3*mlev));
}

void pHOT::Repartition(unsigned mlevel)
{
#ifdef USE_GPTL
  GPTLstart("pHOT::Repartition");
  GPTLstart("pHOT::Repartition::entrance_waiting");
  (*barrier)("pHOT: repartition entrance wait", __FILE__, __LINE__);
  GPTLstop ("pHOT::Repartition::entrance_waiting");
#endif

  PartMapItr it;
  
  volume = sides[0]*sides[1]*sides[2]; // Total volume of oct-tree region


				// No need to repartition 
				// if there are no bodies
  if (cc->nbodies_tot==0) {
    if (myid==0) 
      cout << "pHOT::Repartition with ZERO bodies, continuing" << endl;
#ifdef USE_GPTL
    GPTLstop("pHOT::Repartition");
#endif
    return;
  }

  // For debugging
  vector<unsigned long> erased;

  timer_repartn.start();

  //
  // Recompute keys and compute new partition
  //
#ifdef USE_GPTL
  GPTLstart("pHOT::Repartition::compute_keys");
#endif

  std::vector<key_wght> keys;

  static unsigned long d1a=0, d1b=0;

  bool have_cuda = false;
#if HAVE_LIBCUDA==1
  have_cuda = true;
#endif
  
  if (use_cuda and have_cuda) {

#if HAVE_LIBCUDA==1
    keyProcessCuda(keys);
#endif

  } // END: use_cuda
  else {

    timer_keygenr.start();

    oob.clear();
    for (it=cc->Particles().begin(); it!=cc->Particles().end(); it++) {
      
      it->second->key = getKey(&(it->second->pos[0]));
      if (it->second->key == 0u) {
	oob.insert(it->first);
      } else {
	if (use_weight) {
	  // Floor effort flag to prevent divide-by-zero
	  it->second->effort =
	    std::max<double>(Particle::effort_default, it->second->effort);
	  
	  // Push onto vector
	  keys.push_back(key_wght(it->second->key, it->second->effort));
	  
	  // Reset effort value with some hysteresis
	  it->second->effort = hystrs*(1.0 - hystrs)*it->second->effort;
	  
	} else {
	  keys.push_back(key_wght(it->second->key, 1.0));
	}
      }
    }
    
    
    if (checkDupes1("pHOT::Repartition: after new keys")) {
      cout << "Process " << myid << " at T=" << tnow 
	   << ", L=" << mlevel
	   << ": duplicate check failed after new keys, cnt="
	   << d1a << endl;
    }
    d1a++;

    if (DEBUG_NOISY) 
      std::cout << "Process " << std::left << std::setw(4) << myid 
		<< ": part #="   << std::setw(10) << cc->Particles().size()
		<< "  key size=" << std::setw(10) << keys.size()
		<< "  oob size=" << std::setw(10) << oob.size() << endl;
    
#ifdef USE_GPTL
    GPTLstop ("pHOT::Repartition::compute_keys");
    GPTLstart("pHOT::Repartition::compute_keys_waiting");
    (*barrier)("pHOT: repartition key wait", __FILE__, __LINE__);
    GPTLstop ("pHOT::Repartition::compute_keys_waiting");
    GPTLstart("pHOT::Repartition::spreadOOB");
#endif

    timer_keygenr.stop();
    timer_keysort.start();

    spreadOOB();

#ifdef USE_GPTL
    GPTLstop ("pHOT::Repartition::spreadOOB");
    GPTLstart("pHOT::Repartition::partitionKeys");
#endif

    partitionKeys(keys, kbeg, kfin);
    
#ifdef USE_GPTL
    GPTLstop ("pHOT::Repartition::partitionKeys");
    GPTLstart("pHOT::bodyList");
#endif

    timer_keysort.stop();

    if (DEBUG_KEYS) {		// Deep debug output
      static unsigned count = 0;
      std::ostringstream sout;
      sout << debugf << "." << myid << "." << count++;
      ofstream out(sout.str());

      for (unsigned i=0; i<keys.size(); i++) {
	out << std::setw( 5) << i
	    << std::setw(18) << hex << keys[i].first
	    << std::setw(18) << dec << keys[i].second
	    << std::endl;
	if (i>0 and keys[i].first < keys[i-1].first)
	  out << "####" << std::endl;
      }
    } // END: DEBUG_KEYS

  } // END: not have_cuda
  
  timer_prepare.start();

  //
  // Nodes compute send list
  //
  loclist = kbeg;
  loclist.push_back(kfin[numprocs-1]);	// End point for binary search

  unsigned Tcnt=0, Fcnt, sum;
  vector<int> sendcounts(numprocs, 0), recvcounts(numprocs, 0);
  vector<int> sdispls(numprocs), rdispls(numprocs);
  
  vector< vector<unsigned> > bodylist(numprocs);
  unsigned t;
  for (it=cc->Particles().begin(); it!=cc->Particles().end(); it++) {
				// Skip an OOB particle
    if (it->second->key == 0u) continue;
				// Look for key in this node's list
    t = find_proc(loclist, it->second->key);
    if (t == numprocs) {
      cerr << "Process " << myid << ": loclist found last entry, "
	   << " key=" << hex << it->second->key << dec
	;
      
      cerr << ", end pt="
	   << hex << loclist.back() << dec
	   << ", index=" << t << endl;
    }
    if (t == myid) continue;
    bodylist[t].push_back(it->first);
    sendcounts[t]++;
  }

  for (unsigned k=0; k<numprocs; k++) Tcnt += sendcounts[k];
  MPI_Reduce(&Tcnt, &sum, 1, MPI_UNSIGNED, MPI_SUM, 0, MPI_COMM_WORLD);

#ifdef USE_GPTL
  GPTLstop ("pHOT::bodyList");
  GPTLstart("pHOT::scatter");
#endif

  timer_scatter.start();
  MPI_Alltoall(&sendcounts[0], 1, MPI_INT, 
	       &recvcounts[0], 1, MPI_INT,
	       MPI_COMM_WORLD);
  timer_scatter.stop();

  for (unsigned k=0; k<numprocs; k++) {
    if (k==0) {
      sdispls[0] = rdispls[0] = 0;
      Fcnt = recvcounts[0];
    } else {
      sdispls[k] = sdispls[k-1] + sendcounts[k-1];
      rdispls[k] = rdispls[k-1] + recvcounts[k-1];
      Fcnt += recvcounts[k];
    }
  }

  if (DEBUG_EXTRA) {			// Set to "true" at top of
					// file to enable

    timer_diagdbg.start();

    if (myid==0){
      cout << string(60, '-') << endl
	   << setw(60) << setfill('-')
	   << "---- Send and receive counts for each process "
	   << endl << setfill(' ') << string(60, '-') << endl;
    }

    for (unsigned k=0; k<numprocs; k++) {
      if (myid==k) {
	cout << "Process " << k << endl;
	for (unsigned m=0; m<numprocs; m++) {
	  cout << setw(4) << m << setw(15) << sendcounts[m]
	       << setw(15) << recvcounts[m] << endl;
	}
      }
      (*barrier)("pHOT: repartition send/receive report", __FILE__, __LINE__);
    }

    timer_diagdbg.stop();
  }
      
#ifdef USE_GPTL
  GPTLstop ("pHOT::scatter");
  GPTLstart("pHOT::exchange");
#endif

  //
  // Debug counter
  //
  if (myid==0) {
    n_xchange += sum;
    m_xchange++;
  }
  if (DEBUG_EXTRA) { // If true, write send and receive list for each
		     // node
    timer_diagdbg.start();
    for (int n=0; n<numprocs; n++) {
      if (myid==n) {
	std::string ul(4, '-'), hdr(60, '-');
	cout << hdr << endl
	     << "---- Repartition: send and receive counts" << endl
	     << hdr << endl
	     << "Process " << myid << ": Tcnt=" << Tcnt 
	     << " Fcnt=" << Fcnt << endl
	     << setw(5) << "proc" << setw(8) << "send"
	     << setw(8) << "delt" << setw(8) << "recv"
	     << setw(8) << "delt" << endl
	     << setw(5) << ul     << setw(8) << ul
	     << setw(8) << ul     << setw(8) << ul
	     << setw(8) << ul     << endl;
	for (int m=0; m<numprocs; m++)
	  cout << setw(5) << m 
	       << setw(8) << sendcounts[m]
	       << setw(8) << sdispls[m]
	       << setw(8) << recvcounts[m]
	       << setw(8) << rdispls[m]
	       << endl;
	cout << hdr << endl;
      }
      (*barrier)("pHOT: repartition send/receive counts", __FILE__, __LINE__);
    }
    timer_diagdbg.stop();
  }


  //
  // Exchange particles between processes
  //

  int ps;
  size_t bufsiz = pf->getBufsize();

  // Allocate send and receive buffers (bytes)
  std::vector<char> psend(Tcnt*bufsiz), precv(Fcnt*bufsiz);

  timer_convert.start();
  for (int toID=0; toID<numprocs; toID++) {
    ps = sdispls[toID];
    for (unsigned i=0; i<sendcounts[toID]; i++) {
      pf->particlePack(cc->Particles()[bodylist[toID][i]], &psend[(ps+i)*bufsiz]);
      cc->Particles().erase(bodylist[toID][i]);
    }
  }
  timer_convert.stop();
  timer_xchange.start();

  // Multiply counts and displacements by particle buffer size
  for (auto & v : sendcounts) v *= bufsiz;
  for (auto & v : recvcounts) v *= bufsiz;
  for (auto & v : sdispls   ) v *= bufsiz;
  for (auto & v : rdispls   ) v *= bufsiz;

  MPI_Alltoallv(&psend[0], &sendcounts[0], &sdispls[0], MPI_CHAR,
		&precv[0], &recvcounts[0], &rdispls[0], MPI_CHAR,
		MPI_COMM_WORLD);

  timer_xchange.stop();
  timer_convert.start();

  if (Fcnt) {
    for (unsigned i=0; i<Fcnt; i++) {
      PartPtr part = std::make_shared<Particle>();
      pf->particleUnpack(part, &precv[i*bufsiz]);
      if (part->mass<=0.0 || std::isnan(part->mass)) {
	cout << "[Repartition, myid=" << myid 
	     << ": crazy body with indx=" << part->indx 
	     << ", mass=" << part->mass  << ", key="
	     << hex << part->key << dec
	     << ", i=" << i << " out of " << Fcnt << "]" << endl;
      }
      cc->Particles()[part->indx] = part;
    }

    // Refresh size of local particle list
    cc->nbodies = cc->particles.size();
  }
  timer_convert.stop();
  timer_prepare.stop();

  //
  // Remake key body index
  //
  keybods.clear();
  unsigned oob1_cnt=0, oob_cnt=0;
  for (PartMapItr n=cc->Particles().begin(); n!=cc->Particles().end(); n++) {
    if (n->second->key==0u) {
      oob1_cnt++;
      continue;
    }
    if (n->second->indx==0) {
      cout << "pHOT::Repartition bad particle indx=0!" << endl;
      oob1_cnt++;
    } else {
      keybods.insert(key_pair(n->second->key, n->second->indx));
    }
  }

  // checkBounds(2.0, "AFTER repartition");

  MPI_Reduce(&oob1_cnt, &oob_cnt, 1, MPI_UNSIGNED, MPI_SUM, 0, MPI_COMM_WORLD);
  unsigned oob_tot = oobNumber();
  if (myid==0 && oob_cnt != oob_tot)
    cout << endl << "pHOT::Repartition: " << oob_cnt << " out of bounds," 
	 << " expected " << oob_tot << endl;

  if (DEBUG_CLEAN) {
    //
    // Sanity checks for bad particle indices and bad particle counts
    //
    std::list<PartMapItr> badP;
    for (PartMapItr ip=cc->particles.begin(); ip!=cc->particles.end(); ip++) {
      if (ip->second->indx==0) {
	cout << "pHOT::Repartition BAD particle in proc=" << myid
	     << ", mass=" << ip->second->mass << ", key="
	     << hex << ip->second->key << dec
	     << endl;
	badP.push_back(ip);
      }
    }

    if (badP.size()) {
      cout << "pHOT::Repartition: removing " << badP.size() << " bad entries " 
	   << "from particle list" << std::endl;
      for (auto i : badP) cc->particles.erase(i);
    }

    // Refresh size of local particle list
    cc->nbodies = cc->particles.size();

    // Count total number of particles as sanity check
    if (DEBUG_SANITY) {
      int nbodies1 = cc->nbodies, nbodies0=0;
      MPI_Reduce(&nbodies1, &nbodies0, 1, MPI_INT, MPI_SUM, 0, MPI_COMM_WORLD);
      if (myid==0) {
	if (nbodies0 != cc->nbodies_tot)
	  std::cout << "pHOT::Repartition: leaving with total # mismatch" 
		    << ", total number="    << nbodies0
		    << ", expected number=" << cc->nbodies_tot
		    << endl;
      }
    }
  }

  if (DEBUG_KEYS) {		// Summary/diaganostic  output
				//
    unsigned int nsiz = cc->Particles().size();
    unsigned int ksiz = keys.size();
    vector<unsigned> nsize(numprocs), ksize(numprocs);

    timer_diagdbg.start();

    MPI_Gather(&nsiz, 1, MPI_UNSIGNED, &nsize[0], 1, MPI_UNSIGNED, 
	       0, MPI_COMM_WORLD);

    MPI_Gather(&ksiz, 1, MPI_UNSIGNED, &ksize[0], 1, MPI_UNSIGNED, 
	       0, MPI_COMM_WORLD);
    
    if (myid==0) {
      ofstream out(debugf.c_str(), ios::app);
      unsigned nhead = 35+2*klen;

      out << setfill('-') << setw(nhead) << '-' << endl  
	  << "---- Post-scatter summary, T = " << tnow << endl
	  << setw(nhead) << '-' << endl << setfill(' ') << left
	  << setw(5)    << "#" << right
	  << setw(klen) << "kbeg" 
	  << setw(klen) << "kfin" 
	  << setw(15)   << "nkeys"
	  << setw(15)   << "bodies" << endl << "#" << endl;
      for (int i=0; i<numprocs; i++) {
	out << left << setw(5) << i << right;
	out << setw(klen) << hex << kbeg[i];
	out << setw(klen) << hex << kfin[i];
	out << dec << setw(15) << ksize[i] << setw(15) << nsize[i] << endl;
      }
      out << setfill('-') << setw(nhead) << '-' << endl << setfill(' ') << left
	  << endl << endl;
    }

    timer_diagdbg.stop();
  }

  if (not use_cuda or not have_cuda) {
    if (checkDupes1("pHOT::Repartition: after exchange")) {
      cout << "Process " << myid << " at T=" << tnow 
	   << ", L=" << mlevel
	   << ": duplicate check failed after exchange, cnt=" 
	   << d1b << endl;
    }
    d1b++;
  }

  timer_repartn.stop();

#ifdef USE_GPTL
  GPTLstop("pHOT::exchange");
  GPTLstop("pHOT::Repartition");
#endif
}

void pHOT::checkLevelLists(const std::string& msg)
{
  unsigned cnt1=0, cnt2=0;

  timer_diagdbg.start();

  for (auto i : clevlst) {
    if (clevels[i.second].find(i.first) == clevels[i.second].end()) cnt1++;
  }

  for (auto j : clevels ) {
    for (auto i : j) {
      if (clevlst.find(i) == clevlst.end()) cnt2++;
    }

  }

  if (cnt2)
    std::cout << msg << ", " << cnt2 
	      << " entries in clevels not in clevlst" << std::endl;

  timer_diagdbg.stop();
}

void pHOT::checkCellLevelList(const std::string& msg)
{
  unsigned missing_frontier_cell = 0;
  unsigned missing_clevlst_cell  = 0;

  timer_diagdbg.start();

  for (auto i : frontier) {

    if (i.second->mykey==1u && i.second->ctotal==0u) {
      cout << "Process " << myid 
	   << ": empty root node in checkCellLevelList" << endl;
      continue;
    }
    if (clevlst.find(i.second) == clevlst.end()) 
      missing_frontier_cell++;
  }

  for (auto i : clevlst) {

      if (frontier.find(i.first->mykey) == frontier.end())
	missing_clevlst_cell++;
    }

  if (missing_frontier_cell)
    cout << msg << ", "	 << missing_frontier_cell
	 << " frontier cells not in level list" << endl;

  if (missing_clevlst_cell)
    cout << msg << ", " << missing_clevlst_cell
	 << " level list cells not on frontier" << endl;

  timer_diagdbg.stop();
}

void pHOT::checkSampleCells(const std::string& msg)
{
  timer_diagdbg.start();

  // Check for missing sample cells

  unsigned cnt=0;
  map<unsigned, unsigned> missing;
  for (auto i : frontier) {
    if (i.second->sample == 0x0) {
      cnt++;
      if (missing.find(i.second->level) == missing.end())
	missing[i.second->level] = 1;
      else
	missing[i.second->level]++;
    }
  }
  
  if (cnt) {
    cout << msg << ", " << cnt << " missing sample cells" << endl << left
	 << setw(6) << "Level" << setw(6) << "Count" << endl
	 << setw(6) << "-----" << setw(6) << "-----" << endl;
    for (auto k : missing)
      cout << left << setw(6) << k.first << setw(6) << k.second << endl;
  }

  // Check for ophaned cells: cells not in sample cell child list

  unsigned bad1=0, bad2=0;
  set<pCell*> orphan1, orphan2;
  for (auto i : frontier) {
    if (i.second->sample) {
      if (!i.second->sampleTest()) {
	bad1++;
	orphan1.insert(i.second);
      }
    }
    if (i.second->parent) {
      bool found = false;
      for (auto v : i.second->parent->children) {
	if (v.second == i.second) found = true;
      }
      if (!found) {
	bad2++;
	orphan2.insert(i.second);
      }
    }
  }
  
  if (bad1) {
    cout << msg << ", " << bad1 << " orphaned frontier cells [sample]" << endl;
    cout << setw(10) << "Cell"   << setw(10) << "Sample" 
	 << setw(10) << "Check"  << endl
	 << setw(10) << "------" << setw(10) << "------" 
	 << setw(10) << "------" << endl;
    for (auto k : orphan1) {
      cout << left << setw(10) << k 
	   << setw(10) << k->sample
	   << setw(10) << k->findSampleCell()
	   << endl;
    }
  }

  if (bad2) {
    cout << msg << ", " << bad2 << " orphaned frontier cells [parent]" << endl;
    for (auto k : orphan2)
      cout << left << setw(10) << k << endl;
  }

  timer_diagdbg.stop();
}


void pHOT::makeCellLevelList()
{
  ostringstream sout;
  ofstream out;

  if (DEBUG_EXTRA) {
    sout << "pHOT_cells." << runtag << "." << myid;
    out.open(sout.str().c_str(), ios::out | ios::app);
  }

				// Make new lists
  clevlst.clear();
  clevels = vector< set<pCell*> >(multistep+1);

  if (DEBUG_EXTRA) {
    out << "Process " << myid << " in makeCellLevelList()" 
	<< ", frontier size=" << frontier.size() << endl;
  }

  unsigned ng=0, nt=0;
  for (auto i : frontier) {
				// Check for empty root node
    if (i.second->mykey==1u && i.second->ctotal==0u) {
      if (DEBUG_EXTRA)
	out << "Process " << myid << " in makeCellLevelList()" 
	    << ", empty root node" << endl;
      
      continue;
    }

    nt++;			// Otherwise, count this one
    i.second->remake_plev();
    clevlst[i.second] = i.second->maxplev;
    clevels[i.second->maxplev].insert(i.second);
    if (i.second->bods.size() == 0) {
      cerr << "Process " << myid 
	   << ": makeCellLevelList has a broken frontier!\n";
    } else {
      ng++;
    }
  }

  if (nt!=ng && DEBUG_EXTRA) {
    out << "Process " << myid << ": made level list with " << ng
	<< " good cells out of " << nt << " expected" << endl;

    std::ostringstream sout;
    sout << "pHOT::makeLevelList, myid=" << std::setw(5) << myid;
    printCellLevelList (out, sout.str());
    checkParticles     (out, sout.str());
  }
}

void pHOT::gatherCellLevelList()
{
  timer_diagdbg.start();

  //
  // Working vectors per node
  //
  vector<unsigned> pcnt(multistep+1, 0);
  vector<unsigned> plev(multistep+1, 0);
  unsigned nlev = cc->particles.size();

  //
  // Enter data
  //
  for (auto p : clevlst) pcnt[p.second]++;

  for (unsigned M=0; M<=multistep; M++) plev[M] = CLevels(M).size();
  
  //
  // Recv vectors; make sure space has been allocated.  std::vector is
  // smart about this.
  //
  Pcnt.resize(multistep+1);
  Plev.resize(multistep+1);

  MPI_Reduce(&pcnt[0], &Pcnt[0], multistep+1, MPI_UNSIGNED, MPI_SUM, 0,
	     MPI_COMM_WORLD);

  MPI_Reduce(&plev[0], &Plev[0], multistep+1, MPI_UNSIGNED, MPI_SUM, 0,
	     MPI_COMM_WORLD);

  MPI_Reduce(&nlev,    &Nlev,    1,           MPI_UNSIGNED, MPI_SUM, 0,
	     MPI_COMM_WORLD);

  timer_diagdbg.stop();
}

void pHOT::printCellLevelList(ostream& out, const std::string& msg)
{
  // Sanity check

  if (Pcnt.size() != Plev.size() || Pcnt.size() != multistep+1) return;

  // OK

  timer_diagdbg.start();

  out << msg << endl;
  out << left << setw(60) << setfill('-') << "-" << endl << setfill(' ')
      << "*** T=" << tnow << "  N=" << Nlev << endl
      << setw(10) << "M" << setw(10) << "number" 
      << setw(10) << "counts" << endl;
  for (unsigned M=0; M<=multistep; M++)
    out << setw(10) << M << setw(10) << Plev[M]
	<< setw(10) << Pcnt[M] << endl;
  out << left << setw(60) << setfill('-') << "-" << endl << setfill(' ');

  timer_diagdbg.stop();
}


void pHOT::adjustCellLevelList(unsigned mlevel)
{
  if (multistep==0) return;	// No need to bother if multistepping is off
				// Otherwise . . . 

  ostringstream sout;
  ofstream out;

  if (DEBUG_EXTRA) {
    sout << "pHOT_cells." << runtag << "." << myid;
    out.open(sout.str().c_str(), ios::out | ios::app);
  }

#ifdef USE_GPTL
  GPTLstart("pHOT::adjustCellLevelList");
#endif

  unsigned ng=0, nt=0, ns=0, m, cnt;
  for (unsigned M=mlevel; M<=multistep; M++) {
    nt += CLevels(M).size();
    cnt = 0;
    if (CLevels(M).size()>0) {
      set<pCell*>::iterator it = CLevels(M).begin(), nit;
      while (it != CLevels(M).end()) {
				// Skip the root cell if it's empty
				// (lazy kludge)
	if ( (*it)->mykey==1u && (*it)->ctotal==0 ) { 
	  nt--; 		// Reduce the node count by one
	  it++;			// Go the next set in the set . . .
	  continue; 
	}

	cnt++;			// Count the (presumably good) cells
	
				// For diagnostic info only
	if ((*it)->bods.size()) ng++;
	else {			// This shouldn't happen
	  cout << "Process " << myid << ": pHOT::adjustCellLevelList: "
	       << cnt << "/" << CLevels(M).size()
	       << " zero!" << endl;
	}
				// Newly computed level:
				// we may move cell down but not up . . .
	m = max<unsigned>(mlevel, (*it)->remake_plev());
	nit = it++;
	if (M!=m) {
	  clevels[m].insert(*nit);
	  clevlst[*nit] = m;
	  clevels[M].erase(nit);
	  ns++;
	}
	
	if (CLevels(M).empty()) break;

      }
    }
  }

				// Diagnostic output . . .
  if (nt!=ng)
    cout << "Process " << myid << ": adjusted level list with " << ng
	 << " good cells out of " << nt << " expected, " << ns
	 << " cells moved" << endl;

  if (DEBUG_EXTRA) {
    std::ostringstream sout;
    sout << "pHOT::adjustLevelList, myid=" << std::setw(5) << myid;

    printCellLevelList (out, sout.str());
    checkParticles     (out, sout.str());
  }

#ifdef USE_GPTL
  GPTLstop("pHOT::adjustCellLevelList");
#endif
}

void pHOT::adjustTree(unsigned mlevel)
{
#ifdef USE_GPTL
  GPTLstart("pHOT::adjustTree");
  GPTLstart("pHOT::keyBods");
#endif
				// Barrier to make sure that the timer
				// gives a sensible measurement of key time
  timer_waiton0.start();   
  (*barrier)("pHOT: repartition key timer [0]", __FILE__, __LINE__);
  timer_waiton0.stop();

  timer_tadjust.start();

  if (clevels.size()==0) makeCellLevelList();

  if (DEBUG_ADJUST) {
    std::ostringstream sout;
    sout << "pHOT::adjustTree, Node " << myid 
	 << ": ERROR bodycell BEFORE adjustTree(), T=" << tnow 
	 << " mlevel=" << mlevel << endl;

    checkBodycell         (sout.str());
    checkPartKeybods      (sout.str(), mlevel);
    checkKeybods          (sout.str());
    checkKeybodsFrontier  (sout.str());
    checkCellFrontier     (sout.str());
    checkCellClevel       (sout.str(), mlevel);
    checkCellClevelSanity (sout.str(), mlevel);
    checkSampleCells      (sout.str());
    checkCellLevelList    (sout.str());
    checkLevelLists       (sout.str());
  }

  adjcnt++;			// For debug labeling only . . .
  
  timer_keymake.start();

  pCell* c;
  key_type newkey, oldkey;
  list<unsigned long> oldp;

  //
  // Exchange list
  //
  timer_keybods.start();
  vector< vector<unsigned long> > exchange(numprocs);

  // 
  // Make body list from frontier cells for this level
  // OOB particles must wait until a full tree build
  //
  for (unsigned M=mlevel; M<=multistep; M++) {
    set<pCell*>::iterator it    = CLevels(M).begin();
    set<pCell*>::iterator itend = CLevels(M).end();
    for (;it!=itend; it++) {
      oldp.insert(oldp.end(), (*it)->bods.begin(), (*it)->bods.end());
    }
  }
  timer_keybods.stop();
  
#ifdef USE_GPTL
  GPTLstop ("pHOT::keyBods");
  GPTLstart("pHOT::keyCells");
#endif

				// Barrier to make sure that the timer
				// gives a sensible measurement of key time
  timer_waiton1.start();   
  (*barrier)("pHOT: repartition key timer [1]", __FILE__, __LINE__);
  timer_waiton1.stop();

  //
  // Update body by body using the oldp list
  //
  unsigned newproc;
  for (auto ip : oldp) {

    timer_keybods.start();
    Particle *p = cc->Part(ip);
    if (p==0) {			// Sanity check
      cout << "Process " << myid 
	   << ": pHOT::adjustTree: ERROR, requested particle index "
	   << "does not exist!" << endl;
    }
    numkeys++;
    timer_keybods.stop();

    //
    // Get and recompute keys
    //
    oldkey = p->key;
    newkey = getKey(&(p->pos[0]));

    //
    // Get this particle's cell
    //
    timer_keybods.start();
    key_key::iterator ij = bodycell.find(oldkey);

				// Bad key sanity check (should NEVER happen)
    if (ij == bodycell.end()) {	//
      cout << "Process " << myid 
	   << ": pHOT::adjustTree: ERROR could not find cell for particle"
	   << " key=" << hex << oldkey << dec << ", index=" << p->indx
	   << " pnumber=" << cc->Number() << " bodycell=" << bodycell.size() 
	   << endl;
      timer_keybods.stop();
      continue;
    }

    key_cell::iterator cit = frontier.find(ij->second.first);

				// Bad cell (should NEVER happen)
    if (cit == frontier.end() ) {	//
      cout << "Process " << myid 
	   << ": pHOT::adjustTree: ERROR could not find expected cell"
	   << " on frontier, count=" << adjcnt
	   << " oldbody=" << hex << oldkey << dec
	   << " newbody=" << hex << newkey << dec
	   << " cell="    << hex << bodycell.find(oldkey)->second.first << dec
	   << " index="   << p->indx 
	   << endl;
      timer_keybods.stop();
      continue;
    }

    //
    // This is the cell for this body
    //
    c = cit->second;

    timer_keybods.stop();
      
    cntr_total++;

    //
    // Are the new and old keys the same?  I.e. same cell?
    // 
    timer_keycomp.start();

    if (newkey != oldkey) {

      cntr_new_key++;
				// Key pairs
				//
      key_pair newpair(newkey, p->indx);
      key_pair oldpair(oldkey, p->indx);

				// Put the particle in a new cell?
				// 
      if ( !(c->isMine(newkey)) ) {
	
	timer_keynewc.start();
	cntr_not_mine++;

	c->Remove(oldpair, &change);

				// Same processor and in bounds?
	if (newkey != 0u) {
	  newproc = find_proc(loclist, newkey);
	  if (newproc != myid) {
	    cntr_ship++;
				// Ship this particle elsewhere	  
	    exchange[newproc].push_back(ip);

	  } else {
				// Add the new pair to the tree 
	    keybods.insert(newpair);
	    c->Add(newpair, &change);
	  }

	} else {		// Out of bounds
	  oob.insert(p->indx);
	}
	
	p->key = newkey;	// Assign the new key to the particle

	timer_keynewc.stop();

      } else {			// Same cell: update body cell index 
				// for the new key
	timer_keyoldc.start();
	cntr_mine++;
				// Update key list
	c->UpdateKeys(oldpair, newpair);
				// Assign the new key to the particle
	p->key = newkey;

	timer_keyoldc.stop();
      }
    }
    timer_keycomp.stop();
    
  }
  timer_keymake.stop();
				// Barrier to make sure that the timer
				// gives a sensible measurement of key time
  timer_waiton2.start();   
  (*barrier)("pHOT: repartition key timer [2]", __FILE__, __LINE__);
  timer_waiton2.stop();

  timer_cupdate.start();

#ifdef USE_GPTL
  GPTLstop ("pHOT::keyCells");
  GPTLstart("pHOT::adjExchange");
#endif


  if (DEBUG_CHECK) {
    std::ostringstream sout;
    sout << "pHOT::adjustTree, Node " << myid 
	 << ", BEFORE particle exchange, T=" << tnow;
    
    checkCellFrontier     (sout.str());
    checkKeybodsFrontier  (sout.str());
    checkCellClevel       (sout.str(), mlevel);
    checkLevelLists       (sout.str());
  }

  //
  // Exchange particles
  //

  Particle part;
  unsigned Tcnt=0, Fcnt, sum;
  std::vector<int> sdispls(numprocs), rdispls(numprocs);
  std::vector<int> sendcounts(numprocs, 0), recvcounts(numprocs, 0);

  timer_prepare.start();
  for (unsigned k=0; k<numprocs; k++) {
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
    
    for (unsigned k=0; k<numprocs; k++) {
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
    size_t bufsiz = pf->getBufsize();

    // Allocate send and receive buffers
    std::vector<char> psend(Tcnt*bufsiz), precv(Fcnt*bufsiz);
    
    timer_convert.start();
    for (int toID=0; toID<numprocs; toID++) {
      ps = sdispls[toID];
      for (unsigned i=0; i<sendcounts[toID]; i++) {
	pf->particlePack(cc->Particles()[exchange[toID][i]], &psend[(ps+i)*bufsiz]);
	cc->Particles().erase(exchange[toID][i]);
      }
    }
    timer_convert.stop();
    
    timer_xchange.start();
    
    // Mulitiply counts and displacements by particle buffer size
    for (auto & v : sendcounts) v *= bufsiz;
    for (auto & v : recvcounts) v *= bufsiz;
    for (auto & v : sdispls   ) v *= bufsiz;
    for (auto & v : rdispls   ) v *= bufsiz;

    MPI_Alltoallv(&psend[0], &sendcounts[0], &sdispls[0], MPI_CHAR,
		  &precv[0], &recvcounts[0], &rdispls[0], MPI_CHAR,
		  MPI_COMM_WORLD);
    
    timer_xchange.stop();
    
    timer_convert.start();
    
    for (unsigned i=0; i<Fcnt; i++) {
      PartPtr part = std::make_shared<Particle>();
      pf->particleUnpack(part, &precv[i*bufsiz]);
      if (part->mass<=0.0 || std::isnan(part->mass)) {
	cout << "[adjustTree, myid=" << myid
	     << ": crazy body indx=" << part->indx 
	     << ", mass=" << part->mass << ", key="
	     << hex << part->key<< dec
	     << ", i=" << i << " out of " << Fcnt << "]" << endl;
      }
      
      cc->Particles()[part->indx] = part;
      
      if (part->key != 0u) {
	key_pair newpair(part->key, part->indx);
	keybods.insert(newpair);
	root->Add(newpair, &change);
      }
    }

    // Refresh size of local particle list
    cc->nbodies = cc->particles.size();
    
    timer_convert.stop();
  }
  
  //
  // Debug counter
  //
  if (myid==0) {
    n_xchange += sum;
    m_xchange++;
    sumstep++;
    if (sum==0) sumzero++;
  }


#ifdef USE_GPTL
  GPTLstop ("pHOT::adjExchange");
  GPTLstart("pHOT::overlap");
#endif


  if (DEBUG_CHECK) {
    std::ostringstream sout;
    sout << "pHOT::adjustTree, Node " << myid 
	 << ", BEFORE cell overlap, T=" << tnow;

    checkCellFrontier    (sout.str());
    checkKeybodsFrontier (sout.str());
    checkCellClevel      (sout.str(), mlevel);
    checkLevelLists      (sout.str());
  }


  //
  // Cell overlap?
  //
  
  key_type headKey=0u, tailKey=0u;
  unsigned head_num=0, tail_num=0;
  
  timer_overlap.start();
  
  size_t bufsiz = pf->getBufsize();

  for (int n=0; n<numprocs-1; n++) {
    
    if (n==myid) {
      if (keybods.size()) {
	key_indx::reverse_iterator it = keybods.rbegin();
	tailKey  = bodycell.find(it->first)->second.first;
	tail_num = frontier[tailKey]->bods.size();
      } else {
	tailKey  = 0u;
	tail_num = 0;
      }
      
      // Send my tail cell info to next node
      //
      MPI_Send(&tailKey, 1, MPI_EXP_KEYTYPE, n+1, 131, MPI_COMM_WORLD);
      MPI_Send(&tail_num, 1, MPI_UNSIGNED, n+1, 132, MPI_COMM_WORLD);
      
      // Get next node's head cell info
      //
      MPI_Recv(&headKey, 1, MPI_EXP_KEYTYPE, n+1, 133, 
	       MPI_COMM_WORLD, MPI_STATUS_IGNORE);
      MPI_Recv(&head_num, 1, MPI_UNSIGNED, n+1, 134,
	       MPI_COMM_WORLD, MPI_STATUS_IGNORE);
      
      if (tailKey==headKey && tailKey!=0u) {
	
	c = frontier[tailKey];	// Cell in question: my last cell
	
	if (tail_num>head_num) { 
	  //
	  // Receive particles
	  //
	  std::vector<char> Precv(head_num*bufsiz);

	  MPI_Recv(&Precv[0], head_num*bufsiz, MPI_CHAR,
		   n+1, 136, MPI_COMM_WORLD, MPI_STATUS_IGNORE);

	  for (int i=0; i<head_num; i++) {
	    PartPtr part = std::make_shared<Particle>();
	    pf->particleUnpack(part, &Precv[i*bufsiz]);
	    cc->Particles()[part->indx] = part;
	    key_pair newpair(part->key, part->indx);
	    keybods.insert(newpair);
	    c->Add(newpair, &change);
	  }
	  
	} else {
	  //
	  // Send particles
	  //
	  unsigned k=0;
	  std::vector<char> Psend(tail_num*bufsiz);
	  vector<unsigned long>::iterator ib;
	  for (auto b : c->bods) {
	    pf->particlePack(cc->Particles()[b], &Psend[k++*bufsiz]);
	    cc->Particles().erase(b);
	  }
	  
	  c->RemoveAll();
			
	  if (c->mykey!=1u) {	// Don't remove the root node
	    
	    // queue for removal from level lists
	    change.push_back(cell_indx(c, REMOVE));
	  
	    // queue for deletion
	    change.push_back(cell_indx(c, DELETE));
	  }

	  MPI_Send(&Psend[0], tail_num*bufsiz, MPI_CHAR, 
		   n+1, 135, MPI_COMM_WORLD);
	}
      }
    }
    
    if (n+1==myid) {

      if (keybods.size()) {
	key_indx::iterator it = keybods.begin(); // Sanity check:
	if (bodycell.find(it->first) == bodycell.end()) {
	  cerr << "In adjustTree: No cell for body=" 
	       << hex << it->first << dec
	       << " bodycell size=" << bodycell.size() << endl;
	  headKey  = 0u;
	  head_num = 0;
	} else {
	  headKey  = bodycell.find(it->first)->second.first;
	  head_num = frontier[headKey]->bods.size();
	}
      } else {
	headKey  = 0u;
	head_num = 0;
      }
      
      // Get previous nodes tail cell info
      MPI_Recv(&tailKey, 1, MPI_EXP_KEYTYPE, n, 131, 
	       MPI_COMM_WORLD, MPI_STATUS_IGNORE);
      MPI_Recv(&tail_num, 1, MPI_UNSIGNED, n, 132,
	       MPI_COMM_WORLD, MPI_STATUS_IGNORE);
      
      // Send head node into to previous node
      MPI_Send(&headKey, 1, MPI_EXP_KEYTYPE, n, 133, MPI_COMM_WORLD);
      MPI_Send(&head_num, 1, MPI_UNSIGNED, n, 134, MPI_COMM_WORLD);
      
      if (tailKey==headKey && tailKey!=0u) {
	
	c = frontier[headKey];	// Cell in question
	
	if (tail_num>head_num) { 
	  //
	  // Send particles
	  //
	  unsigned k=0;
	  std::vector<char> Psend(head_num*bufsiz);
	  for (auto b : c->bods) {
	    pf->particlePack(cc->Particles()[b], &Psend[k++*bufsiz]);
	    cc->Particles().erase(b);
	  }
	  
	  c->RemoveAll();
	  
	  if (c->mykey!=1u) { // Dont remove the root node!

	    // queue for removal from level lists
	    change.push_back(cell_indx(c, REMOVE));
	  
	    // queue for deletion
	    change.push_back(cell_indx(c, DELETE));
	  }

	  MPI_Send(&Psend[0], head_num*bufsiz, MPI_CHAR, 
		   n, 136, MPI_COMM_WORLD);
	  
	} else {		
	  //
	  // Receive particles
	  //
	  std::vector<char> Precv(tail_num*bufsiz);
	  MPI_Recv(&Precv[0], tail_num*bufsiz, MPI_CHAR,
		   n, 135, MPI_COMM_WORLD, MPI_STATUS_IGNORE);

	  for (int i=0; i<tail_num; i++) {
	    PartPtr part = std::make_shared<Particle>();
	    pf->particleUnpack(part, &Precv[i*bufsiz]);
	    cc->Particles()[part->indx] = part;
	    key_pair newpair(part->key, part->indx);
	    keybods.insert(newpair);
	    c->Add(newpair, &change);
	  }
	}
      }
    }
    
  }


  timer_overlap.stop();

#ifdef USE_GPTL
  GPTLstop ("pHOT::overlap");
  GPTLstart("pHOT::cUpdate");
#endif

  if (DEBUG_CHECK) {
    std::ostringstream sout;
    sout << "pHOT::adjustTree, Node " << myid 
	 << ", BEFORE cell list reconstruction, T=" << tnow;

    checkCellFrontier    (sout.str());
    checkKeybodsFrontier (sout.str());
    checkCellClevel      (sout.str(), mlevel);
    checkLevelLists      (sout.str());
  }

  //
  // Compute the physical states in each cell for the entire tree and
  // find the sample cells, if anything changed
  //
  timer_cellcul.start();

  computeCellStates();

  //
  // Work around for OpenMP (maybe this is better anyhow)
  //
  std::vector<pCell*> createL, removeL, deleteL, recompL;
  for (auto i : change) {
    if (i.second == CREATE) createL.push_back(i.first);
    if (i.second == REMOVE) removeL.push_back(i.first);
    if (i.second == DELETE) deleteL.push_back(i.first);
    if (i.second == RECOMP) recompL.push_back(i.first);
  }
  change.clear();		// Reset the change list for next time
  
  if (DEBUG_CHECK) {

    // Deep check (verbose)
    //
    if (false) {
      unsigned preLevlst = 0, preFrontier = 0;
      for (int i=0; i<removeL.size(); i++) {
	pCell *c = removeL[i];
	if (clevlst.find(c) == clevlst.end()) {
	  preLevlst++;
	  if (frontier.find(c->mykey) == frontier.end())
	    preFrontier++;
	}
      }
      
      if (preLevlst)
	std::cout << "pHOT::adjustTree, Node " << myid 
		  << ", before reconstruction: "
		  << preLevlst << "/" << clevlst.size()
		  << " cells missing from level list and " << preFrontier
		  << "/" << frontier.size() << " cells missing from frontier" 
		  << std::endl;
    }
    
    std::ostringstream sout;
    sout  << "pHOT::adjustTree, Node " << myid 
	  << ", before reconstruction, "
	  << "ERROR mlevel=" << mlevel << ", T=" << tnow;

    checkBodycell(sout.str());
    
    // Check the remove list for duplicates . . . 
    sort(removeL.begin(), removeL.end());
    unsigned imult = 0;
    for (size_t i=1; i<removeL.size(); i++) {
      if (removeL[i-1] == removeL[i]) imult++;
    }
    if (imult) {
      std::cout << "pHOT::adjustTree: remove list has " << imult 
		<< " multiple entries" << std::endl;
    }
  }

  //
  // Create the new leaves
  //
#pragma omp parallel for default(shared)
  for (int i=0; i<createL.size(); i++) {
    pCell *c = createL[i];

    // Only add cells with bodies.  Newly added cells may be branches
    // or scheduled for later deletion.
    if (c->bods.size()) {
#pragma omp critical
      {
	unsigned m = max<unsigned>(c->maxplev, mlevel);
	clevlst[c] = m;
	clevels[m].insert(c);
      }

      // Locate the new sample cell
      c->findSampleCell("adjustTree<create leaves>");

    }
  }

  //
  // Remove the former leaves from the lists
  //
#pragma omp parallel for default(shared)
  for (int i=0; i<removeL.size(); i++) {
    pCell *c = removeL[i];
    if (DEBUG_CHECK) {
      if (clevlst.find(c) == clevlst.end()) {
	std::cout << "pHOT::adjustTree: cell=" << std::hex << c
		  << std::dec << " not in level list";
	if (frontier.find(c->mykey) == frontier.end())
	  std::cout << " and gone from frontier";
	else
	  std::cout << " but is in frontier";
	std::cout << std::endl;
	// continue;
      }
    }
    unsigned m = clevlst[c];
#pragma omp critical
    {
      clevlst.erase(c);
      if (DEBUG_CHECK) {
	if (clevels[m].find(c) == clevels[m].end()) {
#ifdef HAVE_OMP_H
	  std::cout << "pHOT::adjustTree(REMOVE) [" << omp_get_thread_num()
		    << "]: cell=" << std::hex << c
#else
	  std::cout << "pHOT::adjustTree(REMOVE) [" << 1
		    << "]: cell=" << std::hex << c
#endif
		    << std::dec << " not in level " << m << std::endl;
	}
      }
      clevels[m].erase(c);
    }
  }
  
  //
  // Delete empty leaves
  //
#pragma omp parallel for default(shared)
  for (int i=0; i<deleteL.size(); i++) {
    delete deleteL[i];
  }
  
  //
  // Recompute sample cells for changed leaves
  //
  std::sort(deleteL.begin(), deleteL.end());
#pragma omp parallel for default(shared)
  for (int i=0; i<recompL.size(); i++) {
    if (!std::binary_search(deleteL.begin(), deleteL.end(), recompL[i])) {
      // Am I still a leaf?  If previous cell reaches bucket limit, it
      // will be split and RECOMP calls may still be in the list
      if (recompL[i]->children.size() == 0)
	recompL[i]->findSampleCell("adjustTree<recompute>");
    }
  }

  
  if (DEBUG_ADJUST) {
    std::ostringstream sout;
    sout << "pHOT::adjustTree at END, Node " << myid
	 << ", ERROR mlevel=" << mlevel;

    checkBodycell         (sout.str());
    checkKeybodsFrontier  (sout.str());
    checkParticles        (cout, sout.str());
    checkFrontier         (cout, sout.str());
    checkCellClevel       (sout.str(), mlevel);
    checkCellClevelSanity (sout.str(), mlevel);
    
    checkDupes2();
    checkIndices();
    checkSampleCells      (sout.str());
    checkCellLevelList    (sout.str());
    checkKeybods          (sout.str());
    checkPartKeybods      (sout.str(), mlevel);
  }
  
  // Check the load balance based on the total effort
  checkEffort(mlevel);

  timer_cellcul.stop();

  timer_cupdate.stop();

  timer_tadjust.stop();

  if (DEBUG_KEYS) {		// Summary/diaganostic  output
				//
    unsigned int nsiz = cc->Particles().size();
    vector<unsigned> nsize(numprocs);

    MPI_Gather(&nsiz, 1, MPI_UNSIGNED, &nsize[0], 1, MPI_UNSIGNED, 
	       0, MPI_COMM_WORLD);

    if (myid==0) {
      ofstream out(debugf.c_str(), ios::app);
      unsigned nhead = 20 + 2*klen;

      out << left << setfill('-') << setw(nhead) << '-' << endl
	  << "---- Post-adjustTree summary [" << mlevel << "], T = " 
	  << tnow << endl << setw(nhead) << '-' << endl << setfill(' ');
      out << left << setw(5) << "#" 
	  << setw(klen) << right << "kbeg" << setw(klen) << "kfin" 
	  << setw(15) << "bodies" << endl << "#" << endl;
      for (int i=0; i<numprocs; i++) {
	out << left << setw(5) << i << right
	    << hex << setw(klen) << kbeg[i]
	    << setw(klen) << kfin[i] << dec
	    << setw(15) << nsize[i]
	    << endl;
      }
      out << left << setfill('-') << setw(nhead) << '-' << endl
	  << setfill(' ') << endl;
    }
  }

#ifdef USE_GPTL
  GPTLstop("pHOT::cUpdate");
  GPTLstop("pHOT::adjustTree");
#endif

}
  

bool pHOT::checkParticles(ostream& out, const std::string& msg, bool pc)
{
  timer_diagdbg.start();

  unsigned cnt, badb=0, badc=0;
  const string membername = "checkParticles";

  // FOR CELL OUTPUT
  map<pCell*, pair<unsigned, unsigned> > outdat;

  for (unsigned M=0; M<=multistep; M++) {
    if (CLevels(M).size()) {
      cnt = 0;
      for (auto i : CLevels(M)) {
	cnt++;
	// FOR CELL OUTPUT
	outdat[i] = pair<unsigned, unsigned>(M, i->bods.size());

	if (i->bods.size()) {
	  if (pc) {
	    for (auto b : i->bods) {
	      if (!cc->Part(b)) {
		out << "pHOT::checkParticles:: M=" << M << ", bad body at "
		    << cnt << "/" << CLevels(M).size() 
		    << " cell=" << hex << i << dec << endl;
		badb++;
	      }
	    }
	  }
	} else {
	  out << "pHOT::checkParticles:: M=" << M << ", zero bods at "
	       << cnt << "/" << CLevels(M).size() 
	       << " cell=" << hex << i << dec << endl;
	  out << "pHOT::checkParticles:: More info for " 
	      << hex << i << dec << ": marked as ";
	  if (i->isLeaf) out << "LEAF, ";
	  else               out << "BRANCH, ";
	  if (i->parent) {
	    if (i->parent->isLeaf) out << "parent IS a leaf. ";
	    else                   out << "parent is NOT a leaf. ";
	  } else {
	    out << "parent does NOT exist. ";
	  }
	  out << "Node has " << i->children.size() << " children and ";
	  if (frontier.find(i->mykey) == frontier.end())
	    out << "is gone from frontier." << std::endl;
	  else
	    out << "is in frontier" << std::endl;

	  badc++;
	}
      }
    }
  }
  
  // OUTPUT CELL LIST
  ostringstream origfile1, backfile1;

  origfile1 << "chkcell." << myid;
  backfile1 << "chkcell." << myid << ".bak";
  if (rename(origfile1.str().c_str(), backfile1.str().c_str())) {
    perror("pHOT");
    ostringstream message;
    message << "error creating backup file <" << backfile1.str() << ">";
  }
  
  ostringstream origfile2, backfile2;

  origfile2 << "chklist." << myid;
  backfile2 << "chklist." << myid << ".bak";
  if (rename(origfile2.str().c_str(), backfile2.str().c_str())) {
    perror("pHOT");
    ostringstream message;
    message << "error creating backup file <" << backfile2.str() << ">";
  }
  
  ofstream out1(origfile1.str().c_str());
  ofstream out2(origfile2.str().c_str());

  if (!out1) {
    ostringstream message;
    message << "error opening output file <" << origfile1.str() << ">";
    bomb(membername, message.str());
  }

  if (!out2) {
    ostringstream message;
    message << "error opening output file <" << origfile2.str() << ">";
    bomb(membername, message.str());
  }

  for (auto l : outdat) {
    out1 << setw(15) << hex << l.first << dec 
	 << setw(8) << l.second.first
	 << setw(8) << l.second.second
	 << endl;
  }

  for (auto l : clevlst) {
    out2 << setw(15) << hex << l.first << dec 
	 << setw(8) << l.second
	 << endl;
  }
  // END OUTPUT CELL LIST

  timer_diagdbg.stop();

  if (msg.size() && (badb || badc) ) {
    out << msg << ": pHOT::checkParticles, bad cell=" << badc;
    if (pc) out << " bad bod=" << badb;
    out << endl;
    return false;
  }
  else return true;

}


bool pHOT::checkFrontier(ostream& out, const std::string& msg)
{
  unsigned bad=0;
  bool good=true;

  timer_diagdbg.start();

  for (unsigned M=0; M<=multistep; M++) {
    if (CLevels(M).size()) {
      for (auto i : CLevels(M)) {
	if (frontier.find(i->mykey) == frontier.end()) {
	  out << "pHOT::checkFrontier error on M=" << M
	      << ", cell=" << hex << i << dec << endl;
	  bad++;
	  good = false;
	}
      }
    }
  }
  
  if (bad && msg.size()) {
    out << msg << ": pHOT::checkFrontier, bad cell=" 
	<< bad << endl;
  }

  timer_diagdbg.stop();

  return good;
}

//
// This big nasty mess should only be used for debugging!
//
bool pHOT::checkDupes1(const std::string& msg)
{
  bool ret = false;
  if (!DEBUG_KEYS) return ret;
  
  timer_diagdbg.start();

  // Check for duplicate keys, duplicate sequence numbers
  multimap<indx_type, key_type> plist;
  set<key_type, less<key_type> > keydup;

#ifdef USE_GPTL
  GPTLstart("pHOT::checkDupes1::duplicate_check");
#endif

  // Make multimap of all index/key pairs
  for (PartMapItr it=cc->Particles().begin(); it!=cc->Particles().end(); it++)
    plist.insert(pair<indx_type, key_type>(it->first, it->second->key));
  
  // Hold unique duplicated indices
  multiset<indx_type> dups;
  map<indx_type, list<key_type> > kdup;
  typedef pair<indx_type, list<key_type> > ptype ;
  multimap<indx_type, key_type>::iterator k=plist.begin(), kl=plist.begin();

  if (k != plist.end()) {

    for (k++; k!=plist.end(); k++) {

      while (k->first == kl->first && k!=plist.end()) {
	// Make a list for this index, if it's a new dup
	if (kdup.find(kl->first)==kdup.end()) {
	  kdup.insert(ptype(kl->first, list<key_type>()));
	  kdup[kl->first].push_back(kl->second);
	}
	
	// Add the key to the list for this index
	kdup[k->first].push_back(k->second);
	
	// Add this index to the dup index multiset
	dups.insert(k->first);
	
	// Add the key to the unique key map
	keydup.insert(k->second);
      }

      kl = k;
    }

  }
  
  if (dups.size()) {
    ret = true;
    ostringstream sout;
    sout << outdir << runtag << ".pHOT_crazy." << myid;
    ofstream out(sout.str().c_str(), ios::app);
    out << endl
	<< "#" << setfill('-') << setw(60) << '-' << setfill(' ') << endl
	<< "#---- " << msg << endl
	<< "#---- Time=" << tnow 
	<< ", N=" << cc->Number() 
	<< ", Duplicates=" << dups.size() 
	<< ", Unique keys=" << keydup.size() 
	<< endl
	<< "#" << setfill('-') << setw(60) << '-' << setfill(' ') << endl;
    for (auto b : kdup) {
      Particle *p = cc->Part(b.first);
      out << setw(10) << b.first << setw(18) << p->mass;
      for (int k=0; k<3; k++) out << setw(18) << p->pos[k];
      out << setw(10) << p->indx
	  << "    "   << hex << p->key << dec << endl
	;
      for (auto k : b.second) {
	out << left << setw(10) << "---" << hex << k << dec << endl;
      }
    }
    out << "#" << setfill('-') << setw(60) << '-' << setfill(' ') << endl;
  }
#ifdef USE_GPTL
  GPTLstop("pHOT::checkDupes1::duplicate_check");
#endif

  timer_diagdbg.stop();

  return ret;
}


bool pHOT::checkDupes2()
{
  timer_diagdbg.start();

  vector<unsigned> indices;
  for (auto i : frontier)
    indices.insert(indices.end(), i.second->bods.begin(), i.second->bods.end());

  sort(indices.begin(), indices.end());
  
  unsigned dup = 0;
  for (unsigned n=1; n<indices.size(); n++) {
    if (indices[n-1] == indices[n]) dup++;
  }

  timer_diagdbg.stop();

  if (dup) {
    cout << "Process " << myid << ": pHOT::checkDupes, dup=" << dup << endl;
    return false;
  }

  return true;
}


void pHOT::checkIndices()
{
  timer_diagdbg.start();

  // All processes make an index list
  //
  vector<unsigned> indices;
  for (auto i : frontier)
    indices.insert(indices.end(), i.second->bods.begin(), i.second->bods.end());
  
  // Add oob particles
  //
  for (auto i : oob) indices.push_back(i);

  // All processes send the index list to the master node
  //
  for (int n=1; n<numprocs; n++) {

    if (myid==0) {
      unsigned icnt;
      MPI_Recv(&icnt, 1, MPI_UNSIGNED, n, 392, MPI_COMM_WORLD,
	       MPI_STATUS_IGNORE);
      vector<unsigned> get(icnt);
      MPI_Recv(&get[0], icnt, MPI_UNSIGNED, n, 393, MPI_COMM_WORLD,
	       MPI_STATUS_IGNORE);

      indices.insert(indices.end(), get.begin(), get.end());

    } else if (myid==n) {

      unsigned icnt = indices.size();
      MPI_Send(&icnt, 1, MPI_UNSIGNED, 0, 392, MPI_COMM_WORLD);
      MPI_Send(&indices[0], icnt, MPI_UNSIGNED, 0, 393, MPI_COMM_WORLD);
    }
    (*barrier)("pHOT: check indicies", __FILE__, __LINE__);
  }

  if (myid==0) {
    //
    // Sort the index list
    //
    sort(indices.begin(), indices.end());
    //
    // Check each entry in the list
    //
    for (unsigned n=1; n<=indices.size(); n++) {
      if (n!=indices[n-1]) {
	cout << "Process " << myid << ": "
	     << "pHOT::checkIndices ERROR: found=" << indices[n-1] 
	     << " expected " << n
	     << " total=" << indices.size() 
	     << " time=" << tnow
	     << endl;
	break;
      }
    }
    //
    // Check the total body count
    //
    if (indices.size() != cc->nbodies_tot)
      cout << "Process " << myid << ": "
	   << "pHOT::checkIndices ERROR: time=" << tnow
	   << " index count=" << indices.size() 
	   << " body count=" << cc->nbodies_tot << endl;
  }

  timer_diagdbg.stop();
}

bool pHOT::checkKeybods(const std::string& msg)
{
  timer_diagdbg.start();

  bool ok = true;
  unsigned cnt=0;
  for (PartMapItr n = cc->Particles().begin(); n!=cc->Particles().end(); n++) {

				// Skip oob particles
    if (n->second->key) {
      key_pair tpair(n->second->key, n->second->indx);
      key_indx::iterator it = keybods.find(tpair);

      if (it==keybods.end()) {
	if (DEBUG_NOISY) {
	  cout << msg << ": checkKeybods: " 
	       << cnt << " unmatched particle, (x, y, z)=("
	       << n->second->pos[0] << ", " << n->second->pos[1] 
	       << ", " << n->second->pos[2] << ")" << endl;
	}
	ok = false;
	cnt++;
      }
    }
  }
  
  if (msg.size() && cnt) {
    cout << msg << ": checkKeybods: " 
	 << cnt << " unmatched particles" << endl;
  }

  timer_diagdbg.stop();

  return ok;
}


bool pHOT::checkKeybodsFrontier(const std::string& msg)
{
  timer_diagdbg.start();

  bool ok = true;
  unsigned pcnt=0, bcnt=0;
  std::set<key_type> check;

  for (PartMapItr n = cc->Particles().begin(); n!=cc->Particles().end(); n++) {
    // Skip oob particles
    //
    if (n->second->key) {
				// a keybods (std::set) key
      key_pair tpair(n->second->key, n->second->indx);
      key_indx::iterator it = keybods.find(tpair);

      if (it==keybods.end()) {
	if (DEBUG_NOISY) {
	  cout << "pHOT::checkKeybodsFrontier, " << msg
	       << ": count=" << std::setw(5) << pcnt 
	       << " unmatched particle" << std::endl 
	       << "(x, y, z)=("
	       << std::setw(16) << n->second->pos[0] << ", " 
	       << std::setw(16) << n->second->pos[1] << ", " 
	       << std::setw(16) << n->second->pos[2] << ")" << std::endl
	       << "(u, v, w)=("
	       << std::setw(16) << n->second->vel[0] << ", " 
	       << std::setw(16) << n->second->vel[1] << ", " 
	       << std::setw(16) << n->second->vel[2] << ")" << std::endl;
	}
	ok = false;
	pcnt++;

      } else {
	//
	// Look for the particle key in the bodycell multimap
	//
	key2Range ij   = bodycell.equal_range(n->second->key);
	key_type  ckey = 0l;

	if (ij.first != ij.second) {
	  key2Itr ijk = ij.first;
	  while (ijk != ij.second) {
	    if (ijk->second.second == tpair.second) {
	      ckey = ijk->second.first;
	      break;
	    }
	    ijk++;
	  }
	}

	if (ckey == 0l) {
	  if (DEBUG_NOISY) {
				// Check for index on in bodycell with
				// a sequential scan
	    key_type fCell = 0l;
	    for (key2Itr kt=bodycell.begin(); kt!=bodycell.end(); kt++) {
	      if (kt->second.second == n->second->indx) {
		fCell = kt->second.first;
		break;
	      }
	    }

	    cout << "pHOT::checkKeybodsFrontier, " << msg
		 << ": body count=" << std::setw(5) << bcnt 
		 << ", key="  << std::hex << n->second->key << std::dec
		 << ", indx="  << n->second->indx
		 << ", found="  << std::hex << fCell << std::dec
		 << ", particle with no bodycell" << std::endl
		 << "(x, y, z)=("
		 << std::setw(16) << n->second->pos[0] << ", " 
		 << std::setw(16) << n->second->pos[1] << ", " 
		 << std::setw(16) << n->second->pos[2] << ")" << std::endl
		 << "(u, v, w)=("
		 << std::setw(16) << n->second->vel[0] << ", " 
		 << std::setw(16) << n->second->vel[1] << ", " 
		 << std::setw(16) << n->second->vel[2] << ")" << std::endl;
	  }
	  ok = false;
	  bcnt++;

	} else {
	  key_cell::iterator kt = frontier.find(ckey);
	  
	  if (kt == frontier.end()) {
	    if (DEBUG_NOISY) {
	      cout << "pHOT::checkKeybodsFrontier, " << msg
		   << ": frontier count=" << std::setw(5) << check.size()
		   << " unmatched particle" << std::endl
		   << "(x, y, z)=("
		   << std::setw(16) << n->second->pos[0] << ", " 
		   << std::setw(16) << n->second->pos[1] << ", " 
		   << std::setw(16) << n->second->pos[2] << ")" << std::endl
		   << "(u, v, w)=("
		   << std::setw(16) << n->second->vel[0] << ", " 
		   << std::setw(16) << n->second->vel[1] << ", " 
		   << std::setw(16) << n->second->vel[2] << ")" << std::endl;
	    }
	    ok = false;		// Count the missing cells
	    if (check.find(ckey) == check.end())  check.insert(ckey);
	  }
	}
      }
    }
  }
  
  
  unsigned fcnt = check.size();

  if (msg.size() && (pcnt || bcnt || fcnt)) {
    cout << msg  << ": checkKeybods: " 
	 << pcnt << " unmatched particles, " 
	 << bcnt << " body/cell entries, "
	 << fcnt << " cells not on frontier" << std::endl;
  }
  
  timer_diagdbg.stop();

  return ok;
}


bool pHOT::checkPartKeybods(const std::string& msg, unsigned mlevel)
{
  unsigned bcelbod = 0;
  unsigned bfrontr = 0;
  bool ok = true;

  timer_diagdbg.start();

  // 
  // Make body list from frontier cells for this level
  //
  for (unsigned M=mlevel; M<=multistep; M++) {
    for (auto i : CLevels(M)) {
      for (auto b : i->bods) {
	key_type key  = cc->Particles()[b]->key;
	unsigned indx = cc->Particles()[b]->indx;

	// Look for cell for this body
	//
	if (bodycell.find(key) == bodycell.end()) {
	  bcelbod++;
	  ok = false;
	}
	// Look for cell in frontier . . .
	//
	if (frontier.find(bodycell.find(key)->second.first) == frontier.end()) {
	  bfrontr++;
	  ok = false;
	}
      }
    }
  }
  //
  // Report
  //
  if (!ok && msg.size()) {
    cout << msg
	 << ": pHOT::checkPartKeybods: ERROR bad bodycell=" << bcelbod
	 << " bad frontier=" << bfrontr << endl;

  }

  timer_diagdbg.stop();

  return ok;
}


bool pHOT::checkBodycell(const std::string& msg)
{
  timer_diagdbg.start();

  bool ok = true;
  unsigned cnt=0;
  for (PartMapItr n=cc->Particles().begin(); n!=cc->Particles().end(); n++) {
    // Ignore OOB particle
    if (n->second->key==0u) continue;

    // Look for the bodycell
    key_key::iterator it = bodycell.find(n->second->key);
    if (it==bodycell.end()) {
      ok = false;
      cnt++;
      if (DEBUG_NOISY) {
	cout << "Process " << myid << ": checkBodycell: " 
	     << cnt << " unmatched particle: key=" << hex
	     << n->second->key << dec << " index=" 
	     << n->second->indx  << ", (x, y, z)=("
	     << n->second->pos[0] << ", " << n->second->pos[1] 
	     << ", " << n->second->pos[2] << ")" << endl;
      }
    }
  }

  if (msg.size() && cnt) {
    cout << msg << ": checkBodycell: " 
	 << cnt << " unmatched particles" << endl;
  }

  timer_diagdbg.stop();

  return ok;
}


bool pHOT::checkCellClevel(const std::string& msg, unsigned mlevel)
{
  bool ok = true;
  timer_diagdbg.start();

  unsigned bad = 0, total = 0;
  for (unsigned M=mlevel; M<=multistep; M++) {
    for (auto i : CLevels(M)) {
      if (clevlst.find(i) == clevlst.end()) bad++;
    }
  }

  if (msg.size() && bad) {
    cout << msg << ": "
	 << bad << "/" << total << " unmatched cells" << endl;
  }

  timer_diagdbg.stop();

  if (bad) return false;
  else return true;
}
  
bool pHOT::checkCellClevelSanity(const std::string& msg, unsigned mlevel)
{
  bool ok = true;
  timer_diagdbg.start();

  unsigned error = 0;
  std::set<pCell*> expect, found;

  for (unsigned M=mlevel; M<=multistep; M++) {
    for (auto i : CLevels(M)) {
      expect.insert(i);
      for (auto b : i->bods) {
	// Look for the cell
	key_key::iterator kk = bodycell.find(Body(b)->key);
	// Is the cell on the frontier?
	  if (kk != bodycell.end()) {
	    key_cell::iterator kc = frontier.find(kk->second.first);
	    if (kc==frontier.end()) {
	      error++;
	    } else {
	      found.insert(kc->second);
	    }
	  }
	}
    }
  }

  // Check lists
  bool size_mismatch = (expect.size() == found.size() ? false : true);
  unsigned missing = 0;
  for (auto i : found) {
    if (expect.find(i) == expect.end()) missing++;
  }

  if (size_mismatch || missing || error) {
    ok = false;

    if (msg.size()) {
      std::cout << msg << ", checkCellClevelSanity: ";
      if (size_mismatch)
	std::cout << "expected " << expect.size() << " and found "
		  << found.size();
      if (error)
	std::cout << ": " << error << " cells not on in expected list";
      if (missing)
	std::cout << ": " << missing << " cells not on frontier";
    }
  }

  return ok;
}
  


bool pHOT::checkCellFrontier(const std::string& msg)
{
  bool ok = true;
  timer_diagdbg.start();

  std::set<key_type> cells;
  for (auto n : bodycell) {
    if (n.first != 0u) cells.insert(n.second.first);
  }

  size_t found=0, empty=0, branch=0, missed=0, total=cells.size(), cnt=0;
  
  for (auto n : cells) {
    cnt++;
    key_cell::iterator it = frontier.find(n);
    if (it != frontier.end()) {
      found++;
      pCell* c = it->second;
      if (c->isLeaf) {
	if (c->bods.size() == 0) {
	  empty++;
	  if (DEBUG_NOISY) {
	    std::cout << "Process " << myid << ": checkCellFrontier: " 
		      << cnt << "/" << total << "---"
		      << "cell key=" << std::hex << c->mykey << std::dec 
		      << " is a leaf with no bodies" << std::endl;
	  }
	}
      } else {
	branch++;
	if (DEBUG_NOISY) {
	  std::cout << "Process " << myid << ": checkCellFrontier: " 
		    << cnt << "/" << total << "---"
		    << "cell key=" << std::hex << c->mykey << std::dec 
		    << " is a not a leaf bu is on the frontier" 
		    << std::endl;
	}
      }
    } else {
      missed++;
      if (DEBUG_NOISY) {
	std::cout << "Process " << myid << ": checkCellFrontier: " 
		  << cnt << "/" << total << "---"
		  << " unmatched cell key=" 
		  << std::hex << n << std::dec << std::endl;
      }
    }
  }

  if (empty || branch || missed || total != frontier.size()) {
    ok = false;
    if (msg.size()) {
      std::cout << msg <<  ", checkCellFrontier: " 
		<< cnt << "/" << total << " cells counted, " << found
		<< " on frontier, " << empty << " empty leaves, "
		<< branch << " branches, and " << missed << " frontier misses"
		<< std::endl;
    }
  }

  timer_diagdbg.stop();

  return ok;
}


void pHOT::checkBounds(double rmax, const char *msg)
{
  timer_diagdbg.start();

  int bad = 0;
  for (PartMapItr n = cc->Particles().begin(); n!=cc->Particles().end(); n++) {
    for (int k=0; k<3; k++) if (fabs(n->second->pos[k])>rmax) bad++;
  }
  if (bad) {
    cout << "Process " << myid << ": has " << bad << " out of bounds";
    if (msg) cout << ", " << msg << endl;
    else cout << endl;
  }

  timer_diagdbg.stop();
}

unsigned pHOT::oobNumber()
{
  unsigned number=0, number1=oob.size();
  MPI_Reduce(&number1, &number, 1, MPI_UNSIGNED, MPI_SUM, 0, MPI_COMM_WORLD);
  return number;
}

void pHOT::checkOOB(vector<unsigned>& sendlist)
{
  timer_diagdbg.start();

  bool aok = true;
  unsigned bcnt=0;
  if (myid==0) {
    vector<unsigned> recvlist(numprocs*numprocs);
    for (unsigned n=1; n<numprocs; n++) {
      MPI_Recv(&recvlist[0], numprocs*numprocs, MPI_UNSIGNED, n, 321, 
	       MPI_COMM_WORLD, MPI_STATUS_IGNORE);
      bool ok = true;
      for (unsigned j=0; j<numprocs*numprocs; j++)
	if (sendlist[j] != recvlist[j]) ok = false;
      if (!ok) {
	cout << "checkOOB: sendlist from #" << n << " is in error" << endl;
	aok = false;
	bcnt++;
      } else {
	cout << "checkOOB: sendlist from #" << n << " is OK!" << endl;
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

  timer_diagdbg.stop();
}


void pHOT::spreadOOB()
{
#ifdef USE_GPTL
  GPTLstart("pHOT::spreadOOB::in_reduce");
#endif
				// 3% tolerance
  const unsigned long tol = 33;

  vector<long> list0(numprocs, 0), delta(numprocs);

  long list1 = oob.size();
  MPI_Allgather(&list1, 1, MPI_LONG, &list0[0], 1, MPI_LONG, MPI_COMM_WORLD);

#ifdef USE_GPTL
  GPTLstop ("pHOT::spreadOOB::in_reduce");
  GPTLstart("pHOT::spreadOOB::spread_comp");
#endif

  double tot = 0.5;		// Round off
  for (unsigned n=0; n<numprocs; n++) tot += list0[n];
  long avg = static_cast<long>(floor(tot/numprocs));

  long maxdif=0;
  map<unsigned, unsigned> nsend, nrecv;
  for (unsigned n=0; n<numprocs; n++) {
    delta[n] = avg - list0[n];
    /* 
       Positive delta ===> Receive some
       Negative delta ===> Send some
       Zero delta     ===> Just right
    */
    maxdif = max<long>(maxdif, abs(delta[n]));
    if (delta[n]>0) nrecv[n] =  delta[n];
    if (delta[n]<0) nsend[n] = -delta[n];
  }

				// Debug output
  if (DEBUG_OOB) {
    ostringstream sout;
    sout << "In spreadOOB, proc=" << std::setw(4) << myid
	 << " maxdif=" << std::setw(12) << maxdif 
	 << " #="  << std::setw(10) << cc->nbodies_tot
	 << " ns=" << std::setw( 8) << nsend.size() 
	 << " rs=" << std::setw( 8) << nrecv.size();
  }

#ifdef USE_GPTL
  GPTLstop("pHOT::spreadOOB::spread_comp");
#endif

				// Don't bother if changes are small
  if (maxdif < cc->nbodies_tot/tol) return;

				// Nothing to send or receive
  if (nsend.size()==0 || nrecv.size()==0) return;

#ifdef USE_GPTL
  GPTLstart("pHOT::spreadOOB::make_list");
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

  if (DEBUG_CHECK) checkOOB(sendlist);

  unsigned Tcnt=0, Fcnt=0;
  for (int i=0; i<numprocs; i++) {
    Tcnt += sendlist[numprocs*myid + i];
    Fcnt += sendlist[numprocs*i + myid];
  }

#ifdef USE_GPTL
  GPTLstop("pHOT::spreadOOB::make_list");
#endif

  // DEBUG (set to "true" to enable send/recv list diagnostic)
  //
  if (DEBUG_EXTRA) {
    if (myid==0) {
      std::string hdr(60, '-');
      cout << hdr << endl
	   <<"---- spreadOOB" << endl
	   << hdr << endl
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

#ifdef USE_GPTL
  GPTLstart("pHOT::spreadOOB::exchange_particles");
#endif

  unsigned ps=0, pr=0;
  std::set<indx_type>::iterator ioob;
  size_t bufsiz = pf->getBufsize();
  char *psend=0, *precv=0;
  std::vector<MPI_Request> rql;
  MPI_Request r;
  int ierr;

  if (Tcnt) psend = new char [Tcnt*bufsiz];
  if (Fcnt) precv = new char [Fcnt*bufsiz];

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
	    pf->particlePack(cc->Particles()[*ioob], &psend[(ps+i)*bufsiz]);
	    cc->Particles().erase(*ioob);
	    if (oob.find(*ioob) == oob.end())
	      cerr << "Process " << myid << ": serious error, oob="
		   << *ioob << endl;
	    else oob.erase(ioob);
	  }
	  rql.push_back(r);
	  if ( (ierr=MPI_Isend(&psend[ps*bufsiz], To*bufsiz, MPI_CHAR, 
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
	  if ( (ierr=MPI_Irecv(&precv[pr*bufsiz], From*bufsiz, MPI_CHAR, 
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
  GPTLstop ("pHOT::spreadOOB::exchange_particles");
  GPTLstart("pHOT::spreadOOB::add_to_particles");
#endif

  //
  // Add particles
  //

  if (Fcnt) {
    Particle part;
    for (unsigned i=0; i<Fcnt; i++) {
      PartPtr part = std::make_shared<Particle>();
      pf->particleUnpack(part, &precv[i*bufsiz]);
      if (part->mass<=0.0 || std::isnan(part->mass)) {
	cout << "[spreadOOB, myid=" << myid 
	     << ", crazy body with indx=" << part->indx 
	     << ", mass=" << part->mass << ", key="
	     << hex << part->key << dec
	     << ", i=" << i << " out of " << Fcnt << "]" << endl;
      }
      cc->Particles()[part->indx] = part;
      oob.insert(part->indx);
    }

    // Refresh size of local particle list
    cc->nbodies = cc->particles.size();
  }

#ifdef USE_GPTL
  GPTLstop("pHOT::spreadOOB::add_to_particles");
#endif

}


void pHOT::partitionKeysHilbert(vector<key_wght>& keys,
				vector<key_type>& kbeg, vector<key_type>& kfin)
{
#ifdef USE_GPTL
  GPTLstart("pHOT::partitionKeysHilbert");
#endif

				// For diagnostics
  Timer *timer_debug;
  if (DEBUG_KEYS && myid==0) {
    timer_debug = new Timer();
    timer_debug->start();
  }
				// Sort the keys
  sort(keys.begin(), keys.end(), wghtKEY);

  vector<key_wght> keylist, keylist1;

  double srate = 1;		// Only used for subsampling

				// Set to <true> to enable subsampling
  if (sub_sample) {
				// Default rate (a bit more than
				// one per DSMC cell)
    const double subsample = 0.1; 

    srate = subsample;		// Actual value, may be reset for
				// consistency

				// Desired number of samples
				//
				// Want at least 32, if possible
				// since this is a good number for
				// a DSMC cell
    unsigned nsamp = 
      max<unsigned>(32, static_cast<unsigned>(floor(srate*keys.size())));
    
				// Too many for particle count
    if (nsamp > keys.size()) nsamp = keys.size();

				// Consistent sampling rate
    if (nsamp) srate = static_cast<double>(nsamp)/keys.size();

    // Subsample the key list, with weight accumulation
    //
    if (keys.size()) {
      double twght = 0.0;
      unsigned j = 1;
      for (unsigned i=0; i<keys.size(); i++) {
	twght += keys[i].second;
	if (static_cast<unsigned>(floor(srate*(i+1)-std::numeric_limits<double>::min())) == j) {
	  keylist1.push_back(key_wght(keys[i].first, twght));
	  twght = 0.0;
	  j++;
	}
      }
    }
  } else {
    keylist1 = keys;
  }

  // DEBUG 
  // (set to DEBUG_KEYS "true" at top of file to enable key range diagnostic)
  //
  if (DEBUG_KEYS) {

    timer_diagdbg.start();

    if (myid==0) {
	ofstream out(debugf.c_str(), ios::app);
	std::string hdr(60, '-');
	out << endl
	    << hdr << endl << setfill('-')
	    << setw(60) << "--- Sampling stats " << setfill(' ') << endl
	    << hdr << endl << left
	    << setw(5)  << "Id"
	    << setw(10) << "#keys" 
	    << setw(15) << "rate" 
	    << setw(10) << "samples"
	    << endl
	    << setw(5)  << "--"
	    << setw(10) << "-----" 
	    << setw(15) << "----" 
	    << setw(10) << "-------"
	    << endl;
    }
    for (int n=0; n<numprocs; n++) {
      if (n==myid) {
	ofstream out(debugf.c_str(), ios::app);
	out << left
	    << setw(5)  << myid
	    << setw(10) << keys.size()
	    << setw(15) << srate
	    << setw(10) << keylist1.size()
	    << endl;
      }

      (*barrier)("pHOT: partitionKeys debug 1", __FILE__, __LINE__);
    }
    if (myid==0) {
      ofstream out(debugf.c_str(), ios::app);
      out << std::string(60, '-') << endl << endl;
    }

    unsigned nhead = 15 + 5*klen;

    if (myid==0) {
      ofstream out(debugf.c_str(), ios::app);
      out << left << setfill('-') << setw(nhead) << "-" << endl
	  << setw(nhead) << "------ Sampled keys " << endl
	  << setw(nhead) << "-" << setfill(' ') << endl
	  << left << setw(5) << "#"
	  << setw(10) << "#keys" << right
	  << setw(klen) << "First"
	  << setw(klen) << "Q1"
	  << setw(klen) << "Median"
	  << setw(klen) << "Q3"
	  << setw(klen) << "Last"
	  << endl
	  << left << setw(5)  << "-"
	  << setw(10) << "-----"
	  << right << setw(klen) << "-----"
	  << setw(klen) << "--"
	  << setw(klen) << "------"
	  << setw(klen) << "--"
	  << setw(klen) << "----"
	  << endl;
    }
    for (int n=0; n<numprocs; n++) {
      if (myid==n) {
	ofstream out(debugf.c_str(), ios::app);
	out << left << setw(5) << n << setw(10) << keys.size();
	if (keys.size()>0) {
	  out << right
	      << hex
	      << setw(klen) << keys[0].first
	      << setw(klen) << keys[keys.size()/4].first
	      << setw(klen) << keys[keys.size()/2].first
	      << setw(klen) << keys[keys.size()*3/4].first
	      << setw(klen) << keys[keys.size()-1].first << dec
	      << endl;
	}
      }
      (*barrier)("pHOT: partitionKeys debug 2", __FILE__, __LINE__);
    }

    if (myid==0) {
      ofstream out(debugf.c_str(), ios::app);
      out << left << setfill('-') << setw(nhead) << "-" << endl 
	  << setfill('-') << endl;
    }

    //
    // set to "true" to enable key list diagnostic
    //
    if (DEBUG_KEYLIST) {
      const unsigned cols = 3;	// # of columns in output
      const unsigned cwid = 35;	// column width

      if (myid==0) {
	ofstream out(debugf.c_str(), ios::app);
	std::string hdr(60, '-');
	out << hdr << endl
	    << std::setw(60) << setfill('-') 
	    << "----- Partition keys " << endl << setfill(' ')
	    << hdr << endl;

	for (unsigned q=0; q<cols; q++) {
	  ostringstream sout;
	  sout << left << setw(4) << "Id" << setw(8) 
	       << "Size" << setw(8) << "Tot" << setw(8) << "Order";
	  out << left << setw(cwid) << sout.str();
	}
	out << endl;
	
	for (unsigned q=0; q<cols; q++) {
	  ostringstream sout;
	  sout << setfill('-') << setw(cwid-5) << '-';
	  out << left << setw(cwid) << sout.str();
	}
	out << endl << endl;
      }
      
      for (unsigned n=0; n<numprocs; n++) {
	if (myid==n) {
	  ofstream out(debugf.c_str(), ios::app);
	  bool ok = true;
	  if (keylist1.size()>1) {
	    for (unsigned j=1; j<keylist1.size(); j++)
	      if (keylist1[j-1]>keylist1[j]) ok = false;
	  }
	  ostringstream sout;
	  sout << left 
	       << setw(4) << myid 
	       << setw(8) << keylist1.size() 
	       << setw(8) << keys.size();
	  if (ok) sout << left << setw(8) << "GOOD";
	  else sout << left << setw(8) << "BAD";
	  
	  out << left << setw(cwid) << sout.str() << flush;
	  if (n % cols == cols-1 || n == numprocs-1) out << endl;
	}
	(*barrier)("pHOT: partitionKeys debug 3", __FILE__, __LINE__);
      }
    }

    timer_diagdbg.stop();
  }
  //
  // END DEBUG
  //

  // Tree aggregation (merge sort) of the entire key list
  //
  // <keylist1> is (and must be) sorted to start for parallelMerge
  //
  // Round-robin exchange (rrMerge) followed by std::sort can be
  // faster than the parallel merge (parallelMerge) for small particle
  // numbers, but was intended for testing only
  //
  (*barrier)("pHOT: partitionKeys before pMerge", __FILE__, __LINE__);
  if (USE_RROBIN)
    rrMerge(keylist1, keylist);
  else
    parallelMerge(keylist1, keylist);

  (*barrier)("pHOT: partitionKeys after pMerge", __FILE__, __LINE__);

  MPI_Barrier(MPI_COMM_WORLD);

  if (myid==0) {

				// Accumulate weight list
    for (unsigned i=1; i<keylist.size(); i++) 
      keylist[i].second += keylist[i-1].second;

				// Normalize weight list
    double ktop = keylist.back().second;
    for (unsigned i=0; i<keylist.size(); i++)  keylist[i].second /= ktop;

    vector<double> frate(numprocs);

				// Use an even rate
    frate[0] = 1.0;
    for (unsigned i=1; i<numprocs; i++) 
      frate[i] = frate[i-1] + 1.0;
    //                        ^
    //                        |
    //                        |
    // Replace with the computation rate for the node <i>
    // for load balancing
    //
    

				// The overhead for computing these is small so
				// no matter if they are not used below

    vector<double>   wbeg(numprocs), wfin(numprocs); // Weights for debugging
    vector<unsigned> pbeg(numprocs), pfin(numprocs); // Counts  for debugging
    
				// Compute the key boundaries in the partition
				//
    for (unsigned i=0; i<numprocs-1; i++) {
      if (keylist.size()) {
	key_wght k(0u, frate[i]/frate[numprocs-1]);
	vector<key_wght>::iterator
	  ret = lower_bound(keylist.begin(), keylist.end(), k, wghtDBL);
	kfin[i] = ret->first;
	wfin[i] = ret->second;
	pfin[i] = ret - keylist.begin();
      }
      else {
	kfin[i] = key_min;
	wfin[i] = 0.0;
	pfin[i] = 0;
      }
    }

    kfin[numprocs-1] = key_max;
    wfin[numprocs-1] = 1.0;
    pfin[numprocs-1] = keylist.size();

    kbeg[0] = key_min;
    wbeg[0] = 0.0;
    pbeg[0] = 0;
    for (unsigned i=1; i<numprocs; i++) {
      kbeg[i] = kfin[i-1];
      wbeg[i] = wfin[i-1];
      pbeg[i] = pfin[i-1];
    }
      
    if (DEBUG_KEYS) {	 // If true, print key ranges for each process
      ofstream out(debugf.c_str(), ios::app);
      unsigned nhead = 5 + 3*15 + 3*10 + 3*klen;
      out << setw(nhead) << setfill('-') << '-' << endl
	  << "---- partitionKeys: keys in list="<< keylist.size() << endl
	  << setw(nhead) << setfill('-') << '-' << endl << setfill(' ')
	  << left << setw(5) << "proc" << right << setw(klen) << "kbeg"
	  << setw(klen) << "kfin" << setw(klen) << "# keys" 
	  << setw(15) << "wbeg" << setw(15) << "wend" << setw(15) << "wdif" 
	  << setw(10) << "pbeg" << setw(10) << "pend" << setw(10) << "pdif" 
	  << endl
	  << left << setw(5) << "----" << right << setw(klen) << "----"
	  << setw(klen) << "----" << setw(klen) << "------" 
	  << setw(15) << "----" << setw(15) << "----" << setw(15) << "----"  
	  << setw(10) << "----" << setw(10) << "----" << setw(10) << "----"  
	  << endl;
      double mdif=0.0, mdif2=0.0;
      for (int i=0; i<numprocs; i++) {
	out << left << setw(5) << i << right
	    << hex << setw(klen) << kbeg[i]
	    << setw(klen) << kfin[i] << dec
	    << setw(klen) << (kfin[i] - kbeg[i])
	    << setw(15) << wbeg[i]
	    << setw(15) << wfin[i]
	    << setw(15) << wfin[i] - wbeg[i]
	    << setw(10) << pbeg[i]
	    << setw(10) << pfin[i]
	    << setw(10) << pfin[i] - pbeg[i]
	    << endl;
	mdif  += wfin[i] - wbeg[i];
	mdif2 += (wfin[i] - wbeg[i])*(wfin[i] - wbeg[i]);
      }
      out << setw(nhead) << setfill('-') << '-' << endl << setfill(' ') 
	  << endl;
      if (numprocs>1) {
	mdif /= numprocs;
	mdif2 = sqrt((mdif2 - mdif*mdif*numprocs)/(numprocs-1));
	out << "----  mean wght = " << setw(10) << keylist.size() << endl
	    << "----  std. dev. = " << setw(10) << keylist.size() << endl;
      }
    }
  }

  MPI_Bcast(&kbeg[0], numprocs, MPI_EXP_KEYTYPE, 0, MPI_COMM_WORLD);
  MPI_Bcast(&kfin[0], numprocs, MPI_EXP_KEYTYPE, 0, MPI_COMM_WORLD);


  if (DEBUG_KEYS) {		// If true, print key totals
    unsigned oobn = oobNumber();
    unsigned tkey1 = keys.size(), tkey0 = 0;
    MPI_Reduce(&tkey1, &tkey0, 1, MPI_UNSIGNED, MPI_SUM, 0, MPI_COMM_WORLD);
    if (myid==0) {
      ofstream out(debugf.c_str(), ios::app);
      out << endl
	  << setfill('-') << setw(60) << '-' << setfill(' ') << endl
	  << "---- partitionKeys" << endl
	  << setfill('-') << setw(60) << '-' << setfill(' ') << endl
	  << "----  list size = " << setw(10) << keylist.size() << endl
	  << "---- total keys = " << setw(10) << tkey0 << endl
	  << "----  total oob = " << setw(10) << oobn << endl
	  << "----      TOTAL = " << setw(10) << tkey0 + oobn << endl
	  << setfill('-') << setw(60) << '-' << setfill(' ') << endl
	  << "----   Time (s) = " << setw(10) << timer_debug->stop()
	  << endl
	  << setfill('-') << setw(60) << '-' << setfill(' ') << endl
	  << endl;
      delete timer_debug;
    }
  }

#ifdef USE_GPTL
  GPTLstop("pHOT::partitionKeysHilbert");
#endif
}


//
// Binary search
//
unsigned pHOT::find_proc(vector<key_type>& keys, key_type key)
{
  unsigned num = keys.size();

  if (num==0) {
    cerr << "pHOT::find_proc: crazy input, no keys!" << endl;
    return 0;
  }

  if (key<=keys[0])     return 0;
  if (key>=keys[num-1]) return num-2;

  unsigned beg=0, end=num-1, cur;

  while (end-beg>1) {
    cur = (beg + end)/2;
    if (key < keys[cur])
      end = cur;
    else
      beg = cur;
  }

  return beg;
}


//
// This routine combines two sorted vectors into one
// larger sorted vector
//
void pHOT::sortCombine(vector<key_wght>& one, vector<key_wght>& two,
		       vector<key_wght>& comb)
{
  int i=0, j=0;
  int n = one.size()-1;
  int m = two.size()-1;
  
  comb = vector<key_wght>(one.size()+two.size());

  for (int k=0; k<n+m+2; k++) {
    if (i > n)
      comb[k] = two[j++];
    else if(j > m)
      comb[k] = one[i++];
    else {
      if(one[i].first < two[j].first)
	comb[k] = one[i++];
      else
	comb[k] = two[j++];
    }
  }
}

//
// This routine combines the initial input vector on
// each node (sorted to start) by a binary merge algorithm
//
void pHOT::parallelMerge(vector<key_wght>& initl, vector<key_wght>& final)
{
  MPI_Status status;
  vector<key_wght> work;
  unsigned n;

#ifdef USE_GPTL
  GPTLstart("pHOT::parallelMerge");
#endif

  // Find the largest power of two smaller than
  // the number of processors
  // 
  int M2 = 1;
  while (M2*2 < numprocs) M2 = M2*2;

  // Combine the particles of the high nodes
  // with those of the lower nodes so that
  // all particles are within M2 nodes
  //
  // NB: if M2 == numprocs, no particles
  // will be sent or received
  //
  if (myid >= M2) {
    n = initl.size();
    if (false && barrier_debug) {
      std::cout << "pHOT::parallelMerge: myid=" << setw(5) << myid 
		<< setw(10) << " sending " << setw(6) << n 
		<< setw(10) << " to id="   << setw(5) << myid-M2 
		<< std::endl;
    }
    MPI_Send(&n, 1, MPI_UNSIGNED, myid-M2, 11, MPI_COMM_WORLD);
    if (n) {
      vector<key_type> one(n);
      vector<double>   two(n);

      for (unsigned k=0; k<n; k++) {
	one[k] = initl[k].first;
	two[k] = initl[k].second;
      }
	
      MPI_Send(&one[0], n, MPI_EXP_KEYTYPE, myid-M2, 12,
	       MPI_COMM_WORLD);

      MPI_Send(&two[0], n, MPI_DOUBLE,      myid-M2, 13,
	       MPI_COMM_WORLD);
    }
#ifdef USE_GPTL
    GPTLstop("pHOT::parallelMerge");
#endif
  }

  vector<key_wght> data = initl;

  //
  // Retrieve the excess particles
  //
  if (myid + M2 < numprocs) {
    MPI_Recv(&n, 1, MPI_UNSIGNED, myid+M2, 11, MPI_COMM_WORLD, &status);
    if (false && barrier_debug) {
      std::cout << "pHOT::parallelMerge: myid=" << setw(5) << myid 
		<< setw(10) << " received " << setw(6) << n 
		<< setw(10) << " from id="  << setw(5) << status.MPI_SOURCE
		<< ", expected id=" << myid+M2 
		<< std::endl;
    }
    if (n) {
      vector<key_type> recv1(n);
      vector<double>   recv2(n);

      MPI_Recv(&recv1[0], n, MPI_EXP_KEYTYPE, status.MPI_SOURCE, 12,
	       MPI_COMM_WORLD, MPI_STATUS_IGNORE);

      MPI_Recv(&recv2[0], n, MPI_DOUBLE,      status.MPI_SOURCE, 13,
	       MPI_COMM_WORLD, MPI_STATUS_IGNORE);
      
      vector<key_wght> recv(n);
      for (unsigned k=0; k<n; k++) {
	recv[k].first  = recv1[k];
	recv[k].second = recv2[k];
      }
				// data=data+new_data
      sortCombine(initl, recv, data);

      if (0) {
	std::cout << "pHOT::parallelMerge: myid=" << myid 
		  << " completed sortCombine" << std::endl;
      }
    }
  }

  (*barrier)("pHOT: parallelMerge, before interative merge", 
	     __FILE__, __LINE__);

  if (myid < M2) {

    //
    // Now do the iterative binary merge
    //
    while (M2 > 1) {

      M2 = M2/2;

      // When M2 = 1, we are on the the last iteration.
      // The final node left will be the root with the entire sorted array.
      
      //
      // The upper half of the nodes send to the lower half and is done
      //
      if (myid >= M2) {
	n = data.size();
	if (false && barrier_debug) {
	  std::cout << "pHOT::parallelMerge: myid=" << setw(5) << myid;
	  std::ostringstream sout;
	  sout << " [" << n << "] ";
	  std::cout << left << setw(10) << " sending" << setw(12) 
		    << sout.str() << setw(6) << " to"
		    << " node=" << setw(5) << myid-M2 << std::endl << right;
	}
	MPI_Send(&n, 1, MPI_UNSIGNED, myid-M2, 11, MPI_COMM_WORLD);
	if (n) {
	  vector<key_type> one(n);
	  vector<double>   two(n);
	  
	  for (unsigned k=0; k<n; k++) {
	    one[k] = data[k].first;
	    two[k] = data[k].second;
	  }
	  
	  MPI_Send(&one[0], n, MPI_EXP_KEYTYPE, myid-M2, 12, 
		   MPI_COMM_WORLD);
	  MPI_Send(&two[0], n, MPI_DOUBLE,      myid-M2, 13, 
		   MPI_COMM_WORLD);
	}
#ifdef USE_GPTL
	GPTLstop("pHOT::parallelMerge");
#endif
	// return;
	break;

      } else {
	MPI_Recv(&n, 1, MPI_UNSIGNED, MPI_ANY_SOURCE, 11, MPI_COMM_WORLD, 
		 &status);
	if (false && barrier_debug) {
	  std::cout << "pHOT::parallelMerge: myid=" << setw(5) << myid;
	  std::ostringstream sout;
	  sout << " [" << n << "] ";
	  std::cout << left << setw(10) << " received" << setw(12)
		    << sout.str() << setw(6) << " from "
		    << " node=" << setw(5) << status.MPI_SOURCE 
		    << std::endl << right;
	}
	if (n) {
	  
	  vector<key_type> recv1(n);
	  vector<double>   recv2(n);

	  
	  MPI_Recv(&recv1[0], n, MPI_EXP_KEYTYPE, status.MPI_SOURCE, 12, 
		   MPI_COMM_WORLD, MPI_STATUS_IGNORE);

	  MPI_Recv(&recv2[0], n, MPI_DOUBLE,      status.MPI_SOURCE, 13, 
		   MPI_COMM_WORLD, MPI_STATUS_IGNORE);
	  
	  vector<key_wght> recv(n);
	  for (unsigned k=0; k<n; k++) {
	    recv[k].first  = recv1[k];
	    recv[k].second = recv2[k];
	  }
	  
	  //
	  // The lower half sorts and loop again
	  //
	  sortCombine(data, recv, work);
	  data = work;
	}
      }
    }
  }

  (*barrier)("pHOT: parallelMerge, after interative merge",
	     __FILE__, __LINE__);

  //
  // We are done, return the result
  //

  final = data;

  /*
  if (myid == 0) 
    std::cout << "pHOT::parallelMerge: data size=" << data.size() << std::endl;
  */

#ifdef USE_GPTL
  GPTLstop("pHOT::parallelMerge");
#endif

  return;
}


//
// Trivial round robin merge for testing parallelMerge
//
void pHOT::rrMerge(vector<key_wght>& initl, vector<key_wght>& final)
{
  Timer timer;			// For debugging

  if (myid==0 && DEBUG_MERGE)  {
    timer.start();   		// Time the trivial version
  }


  final = initl;		// Put this node's keys on the list

  unsigned n = initl.size(), total = 0;
  unsigned p = n;
				// Extend final vector to total size
  MPI_Allreduce(&n, &total, 1, MPI_UNSIGNED, MPI_SUM, MPI_COMM_WORLD);
  final.resize(total);

  for (int id=0; id<numprocs; id++) {
				// Current node's size
    if (id == myid) n = initl.size();
    MPI_Bcast(&n, 1, MPI_UNSIGNED, id, MPI_COMM_WORLD);

    vector<key_type> one(n);
    vector<double>   two(n);
				// Load the temporary vectors
    if (id==myid) {
      for (unsigned k=0; k<n; k++) {
	one[k] = initl[k].first;
	two[k] = initl[k].second;
      }
    }
				// Send the temporary vectors
    MPI_Bcast(&one[0], n, MPI_EXP_KEYTYPE, id, MPI_COMM_WORLD);
    MPI_Bcast(&two[0], n, MPI_DOUBLE,      id, MPI_COMM_WORLD);
    
				// Load the final vector
    if (id!=myid) {
      for (unsigned k=0; k<n; k++) {
	final[p+k].first  = one[k];
	final[p+k].second = two[k];
      }
      p += n;
    }
  }

				// Do the global sort
  std::sort(final.begin(), final.end());

  if (DEBUG_MERGE) {

    if (myid==0) {
      std::cout << "Trivial sort in " << timer.stop()
		<< " seconds" << std::endl;

				// Time the parallel version
      timer.reset();
      timer.start();
    }

    std::vector<key_wght> final0;

    parallelMerge(initl, final0);

    if (myid==0) {
      std::cout << "Parallel sort in " << timer.stop()
		<< " seconds" << std::endl;

				// Check sizes
      if (final0.size() != final.size()) {
	std::cout << "pHOT::rrMerge: sizes differ, " << final0.size()
		  << " [parallel] vs " << final.size() << " [round-robin]" 
		  << std::endl
		  << "pHOT::rrMerge: input size is " << total << std::endl;
      } else {
	unsigned cnt = 0;
	for (int i=0; i<final.size(); i++) {
	  if (final0[i].first != final[i].first) cnt++;
	}
	if (cnt) {
	  std::cout << "pHOT::rrMerge: " << cnt << " differing key values"
		    << std::endl;
	} else {
	  std::cout << "pHOT::rrMerge: identical array" << std::endl;
	}
      }
    }
  }

  return;
}


unsigned pHOT::checkNumber()
{
  timer_diagdbg.start();

  unsigned nbods1=0, nbods=0;
  for (auto i : frontier) 
    nbods1 += i.second->bods.size();
  
  MPI_Reduce(&nbods1, &nbods, 1, MPI_UNSIGNED, MPI_SUM, 0, MPI_COMM_WORLD);

  timer_diagdbg.stop();

  return nbods;
}

void pHOT::CollectTiming()
{
  float fval;

  fval = timer_keymake.getTime();
  MPI_Gather(&fval, 1, MPI_FLOAT, &keymk3[0], 1, MPI_FLOAT, 0, MPI_COMM_WORLD);

  fval = timer_xchange.getTime();
  MPI_Gather(&fval, 1, MPI_FLOAT, &exchg3[0], 1, MPI_FLOAT, 0, MPI_COMM_WORLD);

  fval = timer_convert.getTime();
  MPI_Gather(&fval, 1, MPI_FLOAT, &cnvrt3[0], 1, MPI_FLOAT, 0, MPI_COMM_WORLD);

  fval  = timer_overlap.getTime();
  MPI_Gather(&fval, 1, MPI_FLOAT, &tovlp3[0], 1, MPI_FLOAT, 0, MPI_COMM_WORLD);

  fval  = timer_prepare.getTime();
  MPI_Gather(&fval, 1, MPI_FLOAT, &prepr3[0], 1, MPI_FLOAT, 0, MPI_COMM_WORLD);

  fval = timer_cupdate.getTime();
  MPI_Gather(&fval, 1, MPI_FLOAT, &updat3[0], 1, MPI_FLOAT, 0, MPI_COMM_WORLD);

  fval = timer_scatter.getTime();
  MPI_Gather(&fval, 1, MPI_FLOAT, &scatr3[0], 1, MPI_FLOAT, 0, MPI_COMM_WORLD);

  fval = timer_repartn.getTime();
  MPI_Gather(&fval, 1, MPI_FLOAT, &reprt3[0], 1, MPI_FLOAT, 0, MPI_COMM_WORLD);

  fval = timer_tadjust.getTime();
  MPI_Gather(&fval, 1, MPI_FLOAT, &tadjt3[0], 1, MPI_FLOAT, 0, MPI_COMM_WORLD);

  fval = timer_cellcul.getTime();
  MPI_Gather(&fval, 1, MPI_FLOAT, &celcl3[0], 1, MPI_FLOAT, 0, MPI_COMM_WORLD);

  fval = timer_keycomp.getTime();
  MPI_Gather(&fval, 1, MPI_FLOAT, &keycm3[0], 1, MPI_FLOAT, 0, MPI_COMM_WORLD);

  fval = timer_keybods.getTime();
  MPI_Gather(&fval, 1, MPI_FLOAT, &keybd3[0], 1, MPI_FLOAT, 0, MPI_COMM_WORLD);

  fval = timer_keysort.getTime();
  MPI_Gather(&fval, 1, MPI_FLOAT, &keyst3[0], 1, MPI_FLOAT, 0, MPI_COMM_WORLD);

  fval = timer_keygenr.getTime();
  MPI_Gather(&fval, 1, MPI_FLOAT, &keygn3[0], 1, MPI_FLOAT, 0, MPI_COMM_WORLD);

  fval = timer_waiton0.getTime();
  MPI_Gather(&fval, 1, MPI_FLOAT, &wait03[0], 1, MPI_FLOAT, 0, MPI_COMM_WORLD);

  fval = timer_waiton1.getTime();
  MPI_Gather(&fval, 1, MPI_FLOAT, &wait13[0], 1, MPI_FLOAT, 0, MPI_COMM_WORLD);

  fval = timer_waiton2.getTime();
  MPI_Gather(&fval, 1, MPI_FLOAT, &wait23[0], 1, MPI_FLOAT, 0, MPI_COMM_WORLD);

  fval = timer_keynewc.getTime();
  MPI_Gather(&fval, 1, MPI_FLOAT, &keync3[0], 1, MPI_FLOAT, 0, MPI_COMM_WORLD);

  fval = timer_keyoldc.getTime();
  MPI_Gather(&fval, 1, MPI_FLOAT, &keyoc3[0], 1, MPI_FLOAT, 0, MPI_COMM_WORLD);

  fval = barrier->getTime();
  MPI_Gather(&fval, 1, MPI_FLOAT, &barri3[0], 1, MPI_FLOAT, 0, MPI_COMM_WORLD);

  fval = timer_diagdbg.getTime();
  MPI_Gather(&fval, 1, MPI_FLOAT, &diagd3[0], 1, MPI_FLOAT, 0, MPI_COMM_WORLD);

  MPI_Gather(&numkeys, 1, MPI_UNSIGNED, &numk3[0], 1, MPI_UNSIGNED, 0, MPI_COMM_WORLD);


  timer_keymake.reset();
  timer_xchange.reset();
  timer_convert.reset();
  timer_overlap.reset();
  timer_prepare.reset();
  timer_cupdate.reset();
  timer_scatter.reset();
  timer_repartn.reset();
  timer_tadjust.reset();
  timer_cellcul.reset();
  timer_keycomp.reset();
  timer_keybods.reset();
  timer_keysort.reset();
  timer_keygenr.reset();
  timer_waiton1.reset();
  timer_waiton2.reset();
  timer_keynewc.reset();
  timer_keyoldc.reset();
  timer_waiton0.reset();
  timer_diagdbg.reset();

  numk = numkeys;
  numkeys = 0;
}


template <typename T> 
void pHOT::getQuant(vector<T>& in, vector<T>& out)
{
  sort(in.begin(), in.end());
  out = vector<T>(ntile+2);
  out[0]       = in.front();
  for (int k=0; k<ntile; k++)
    out[k+1]   = in[static_cast<int>(floor(in.size()*0.01*qtile[k]))];
  out[ntile+1] = in.back();
}

void pHOT::Timing(vector<float>    &keymake, vector<float>    &exchange, 
		  vector<float>    &convert, vector<float>    &overlap, 
		  vector<float>    &prepare, vector<float>    &update,
		  vector<float>    &scatter, vector<float>    &repartn,
		  vector<float>    &tadjust, vector<float>    &cellcul,
		  vector<float>    &keycomp, vector<float>    &keybods,
		  vector<float>    &keysort, vector<float>    &keygenr,
		  vector<float>    &waiton0, vector<float>    &waiton1,
		  vector<float>    &waiton2, vector<float>    &keynewc,
		  vector<float>    &keyoldc, vector<float>    &treebar,
		  vector<float>    &diagdbg, vector<unsigned> &numk)
{
  getQuant<float   >(keymk3, keymake);
  getQuant<float   >(exchg3, exchange);
  getQuant<float   >(cnvrt3, convert);
  getQuant<float   >(tovlp3, overlap);
  getQuant<float   >(prepr3, prepare);
  getQuant<float   >(updat3, update);
  getQuant<float   >(scatr3, scatter);
  getQuant<float   >(reprt3, repartn);
  getQuant<float   >(tadjt3, tadjust);
  getQuant<float   >(celcl3, cellcul);
  getQuant<float   >(keycm3, keycomp);
  getQuant<float   >(keybd3, keybods);
  getQuant<float   >(keyst3, keysort);
  getQuant<float   >(keygn3, keygenr);
  getQuant<float   >(wait03, waiton0);
  getQuant<float   >(wait13, waiton1);
  getQuant<float   >(wait23, waiton2);
  getQuant<float   >(keync3, keynewc);
  getQuant<float   >(keyoc3, keyoldc);
  getQuant<float   >(barri3, treebar);
  getQuant<float   >(diagd3, diagdbg);
  getQuant<unsigned>(numk3,  numk   );
}


double pHOT::totalKE(double& KEtot, double& KEdsp)
{
  vector<double> state(10);

  // Test
  //
  vector<double> state1(10, 0.0);

  for (auto i : frontier) {
    for (unsigned k=0; k<10; k++) 
      state1[k] += i.second->stotal[k];
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

void pHOT::totalMass(unsigned& Counts, double& Mass)
{
  double mass1    = root->stotal[0];
  unsigned count1 = root->ctotal;

  MPI_Reduce(&mass1,  &Mass,   1, MPI_DOUBLE,   MPI_SUM, 0, MPI_COMM_WORLD);
  MPI_Reduce(&count1, &Counts, 1, MPI_UNSIGNED, MPI_SUM, 0, MPI_COMM_WORLD);
}


void pHOT::adjustCounts(ostream& out)
{
  out << setfill('-') << setw(60) << '-' << setfill(' ') << endl;
  out << "  " << left << setw(20) << "Total"       << " : " << cntr_total 
      << endl;
  out << "  " << left << setw(20) << "New key"     << " : " << cntr_new_key 
      << endl;
  out << "  " << left << setw(20) << "Same cell"   << " : " << cntr_mine 
      << endl;
  out << "  " << left << setw(20) << "New cell"    << " : " << cntr_not_mine 
      << endl;
  out << "  " << left << setw(20) << "Shipped"     << " : " << cntr_ship
      << endl;
  out << setfill('-') << setw(60) << '-' << setfill(' ') << endl;

  cntr_total = cntr_new_key = cntr_mine = cntr_not_mine = cntr_ship = 0;
}

void pHOT::checkEffort(unsigned mlevel)
{
  if (!DEBUG_KEYS) return;

  timer_diagdbg.start();

  double eff=0.0;
  for (auto i : frontier) {
    vector<unsigned long>::iterator ib;
    for (auto b : i.second->bods)
      eff += cc->Particles()[b]->effort;
  }

  vector<double> eq(numprocs);
  MPI_Gather(&eff, 1, MPI_DOUBLE, &eq[0], 1, MPI_DOUBLE, 0, MPI_COMM_WORLD);

  if (myid==0) {
    ofstream out(string(outdir + runtag + ".pHOT_effort").c_str(), ios::app);
    if (out) {
      double mean=0.0, mean2=0.0, emin=1e20, emax=0.0;
      for (int n=0; n<numprocs; n++) {
	mean  += eq[n];
	mean2 += eq[n]*eq[n];
	emin = min<double>(emin, eq[n]);
	emax = max<double>(emax, eq[n]);
      }
      mean /= numprocs;
      mean2 = mean2 - mean*mean*numprocs;
      if (numprocs>1) mean2 = sqrt(mean2/(numprocs-1));
      out << setw(40) << setfill('-') << '-' << setfill(' ') << endl
	  << "---- Time  = " << tnow << endl;
      if (mlevel == -1)
	out << "---- Level = " << "full tree" << endl;
      else
	out << "---- Level = " << mlevel  << endl;
      out << "---- Mean  = " << mean      << endl
	  << "---- Stdev = " << mean2     << endl
	  << "---- E_min = " << emin      << endl
	  << "---- E_max = " << emax      << endl
	  << "---- Ratio = " << emax/emin << endl
	  << setw(8) << left << "Proc #" << setw(12) << "Effort" << endl;
      for (int n=0; n<numprocs; n++)
	out << left << setw(8) << n << setw(12) << eq[n] << endl;
    }
  }

  timer_diagdbg.stop();
}

