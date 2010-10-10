#include <values.h>

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

using namespace std;

#include "global.H"
#include "pHOT.H"

double   pHOT::sides[] = {2.0, 2.0, 2.0};
double   pHOT::offst[] = {1.0, 1.0, 1.0};
unsigned pHOT::qtile[] = {10,  50,  90 };


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
// Write verbose output to a file if true (for debugging)
// [set to false for production]
//
bool pHOT::keys_debug = true;

//
// Turn on/off subsampling the key list for partitioning
//
bool pHOT::sub_sample = true;

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

void pHOT::bomb(const string& membername, const string& msg)
{
  cerr << "pHOT::" << membername << "(): " << msg << endl;
  MPI_Abort(MPI_COMM_WORLD, 497);
}

/*
  Constructor: initialize domain
*/
pHOT::pHOT(Component *C)
{
  cc = C;			// Register the calling component
			
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

  offset = new double [3];
  for (unsigned k=0; k<3; k++) offset[k] = offst[k];

  kbeg = vector<key_type>(numprocs);
  kfin = vector<key_type>(numprocs);

  key_min = key_type(1u) << (nbits*3);
  key_max = key_type(1u) << (nbits*3+1);

  m_xchange = n_xchange = 0;

  sumstep = sumzero = 0;

  cntr_total = cntr_new_key = cntr_mine = cntr_not_mine = cntr_ship = 0;

  numkeys = 0;
  timer_keymake.Microseconds();
  timer_xchange.Microseconds();
  timer_convert.Microseconds();
  timer_overlap.Microseconds();
  timer_prepare.Microseconds();
  timer_cupdate.Microseconds();
  timer_scatter.Microseconds();
  timer_repartn.Microseconds();
  timer_tadjust.Microseconds();
  timer_keycall.Microseconds();
  timer_keycomp.Microseconds();
  timer_keybods.Microseconds();
  timer_waiton0.Microseconds();
  timer_waiton1.Microseconds();
  timer_waiton2.Microseconds();
  timer_keynewc.Microseconds();
  timer_keyoldc.Microseconds();

  use_weight = true;
 
  // Initialize timing structures
  //
  keymk3 = exchg3 = cnvrt3 = tovlp3 = prepr3 = updat3 =
    scatr3 = reprt3 = tadjt3 = keycl3 = keycm3 = keybd3 =
    wait03 = wait13 = wait23 = keync3 = keyoc3 = barri3 =
    vector<float>(numprocs);
  
  numk3    = vector<unsigned>(numprocs);
  numfront = vector<int>(numprocs);
  displace = vector<int>(numprocs);

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

pHOT::~pHOT()
{
  delete barrier;
  delete root;
  delete [] offset;
}


key_type pHOT::getKey(double *p)
{
#ifdef USE_GPTL
  GPTLstart("pHOT::getKey");
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
  for (unsigned k=0; k<3; k++) { 
    if (fabs((p[k]+offset[k])/sides[k])> 1.0) {	
#ifdef DEBUG
      cout << "Coordinate out of pbounds in pHOT::key: ";
      for (int l=0; l<3; l++) cout << setw(18) << p[l];
      cout << endl;
#endif
#ifdef USE_GPTL
      GPTLstop("pHOT::getKey");
#endif
      return key_type(0u);
    }
  }

#ifdef INT128
  const double factor = (key_type(1u)<<nbits).toDouble();
#else
  const double factor = static_cast<double>(key_type(1u)<<nbits);
#endif
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
#ifdef INT128
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


void pHOT::makeTree()
{
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

#ifdef DEBUG
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
    (*barrier)("pHOT: pHOT_storage");
  }
#endif  

  //
  // Make the root
  //
  root = new pCell(this);

  //
  // Make new offset
  //
  for (unsigned k=0; k<3; k++) offset[k] = offst[k];
  MPI_Bcast(&offset[0], 3, MPI_DOUBLE, 0, MPI_COMM_WORLD);

  //
  // Add the data
  //
  key_indx::iterator it;
  pCell* p = root;
  for (it=keybods.begin(); it!=keybods.end(); it++)  {

    bool a = (it->first < key_min);
    bool b = (it->first >= key_max);

    if (it->first < key_min || it->first >= key_max) {
#ifdef INT128
      cout << "Process " << myid << ": in makeTree, key=" 
	   << it->first.toHex() << endl << dec;
#else
      cout << "Process " << myid << ": in makeTree, key=" 
	   << hex << it->first << endl << dec;
#endif
    }
    p = p->Add(*it);		// Do the work
  }

				// Sanity checks and debugging
  if (false) {
    if (true) {			// Report on particle sizes for each node
      unsigned long bdcel1=cc->Particles().size(), bdcelmin, bdcelmax, bdcelsum, bdcelsm2;
      MPI_Reduce(&bdcel1, &bdcelmin, 1, MPI_UNSIGNED_LONG, MPI_MIN, 0, MPI_COMM_WORLD);
      MPI_Reduce(&bdcel1, &bdcelmax, 1, MPI_UNSIGNED_LONG, MPI_MAX, 0, MPI_COMM_WORLD);
      MPI_Reduce(&bdcel1, &bdcelsum, 1, MPI_UNSIGNED_LONG, MPI_SUM, 0, MPI_COMM_WORLD);
      bdcel1 = bdcel1*bdcel1;
      MPI_Reduce(&bdcel1, &bdcelsm2, 1, MPI_UNSIGNED_LONG, MPI_SUM, 0, MPI_COMM_WORLD);
      if (myid==0)
	cout << "In makeTree, particles min=" << bdcelmin << " max=" << bdcelmax
	     << " mean=" << bdcelsum/numprocs
	     << " stdv=" << sqrt( (bdcelsm2 - bdcelsum*bdcelsum/numprocs)/(numprocs-1) )
	     << endl;
    }
    if (true) {		// Report on body counts for each node
      if (myid==0) {
	cout << left << setfill('-')
	     << setw(4) << '+' << setw(15) << '+' << endl << setfill(' ')
	     << setw(4) << "pid" << setw(15) << "# bodies" << endl << setfill('-')
	     << setw(4) << '+' << setw(15) << '+' << endl << setfill(' ') << right;
      }
      for (int i=0; i<numprocs; i++) {
	if (myid==i) cout << setw(4) << i << setw(15) << cc->Particles().size() << endl;
	(*barrier)("pHOT: body count report");
      }
      if (myid==0) 
	cout << setfill('-') << setw(4) << '+' << setw(15) << '+' << endl << setfill(' ');
    }
    if (bodycell.size()==0) {
      cout << "Process " << myid << ": in makeTree, unusual condition #bodycell=0"
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

      MPI_Send(&tailKey,  1, MPI_EXP_KEYTYPE, n, 1000, MPI_COMM_WORLD);
      MPI_Send(&tail_num, 1, MPI_UNSIGNED,    n, 1001, MPI_COMM_WORLD);

      MPI_Recv(&nextKey, 1, MPI_EXP_KEYTYPE,  n, 1000, MPI_COMM_WORLD,
	       MPI_STATUS_IGNORE);
      MPI_Recv(&next_num, 1, MPI_UNSIGNED,    n, 1001, MPI_COMM_WORLD,
	       MPI_STATUS_IGNORE);

      if (tailKey != 0u && tailKey == nextKey) {
	if (tail_num <= next_num) {
	  if (tail_num)
	    sendCell(tailKey, n, tail_num);
	  else
	    cout << "Process " << myid << ": not sending cell with zero particles" << endl;
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

      if (keybods.size()) {

#ifdef DEBUG
				// check validity of key
	if (bodycell.find(keybods.begin()->first) == bodycell.end()) {
	  cout << "Process " << myid << ": bad key=" 
	       << hex << keybods.begin()->first << dec
	       << " #cells=" << bodycell.size() << endl;
	}
#endif
	
	headKey = bodycell.find(keybods.begin()->first)->second;
				// Number of bodies in my head cell
#ifdef DEBUG
	// Debug: check for key in frontier
	if (frontier.find(headKey) == frontier.end()) {
	  cout << "Process " << myid << ": headKey=" 
	       << headKey << dec << " not in frontier!" << endl;
	}
#endif
	//
	head_num = frontier[headKey]->bods.size();
      } else {
	headKey  = 0u;
	head_num = 0u;
      }

      MPI_Send(&headKey,  1, MPI_EXP_KEYTYPE, n-1, 1000, MPI_COMM_WORLD);
      MPI_Send(&head_num, 1, MPI_UNSIGNED,    n-1, 1001, MPI_COMM_WORLD);

      MPI_Recv(&prevKey,  1, MPI_EXP_KEYTYPE, n-1, 1000, MPI_COMM_WORLD,
	       MPI_STATUS_IGNORE);
      MPI_Recv(&prev_num, 1, MPI_UNSIGNED,    n-1, 1001, MPI_COMM_WORLD,
	       MPI_STATUS_IGNORE);

      if (headKey != 0u && headKey == prevKey) {
	if (head_num < prev_num) {
	  if (head_num)
	    sendCell(headKey, n-1, head_num);
	  else
	    cout << "Process " << myid << ": not sending cell with zero particles" << endl;
	} else {
	  if (prev_num)
	    recvCell(n-1, prev_num);
	  else
	    cout << "Process " << myid << ": not receiving cell with zero particles" << endl;
	}
      }

    }    

  }


#ifdef USE_GPTL
  GPTLstop ("pHOT::makeTree::adjustBoundaries");
  GPTLstart("pHOT::makeTree::getFrontier");
#endif

  // Each node start at root and walk down the tree to zero the counts
  //
  root->zeroState();

  // March through the frontier and accumulate the counts
  //
  for (key_cell::iterator it=frontier.begin(); 
       it != frontier.end(); it++) it->second->accumState();

  // March through the frontier to find the sample cells
  //
  for (key_cell::iterator it=frontier.begin(); it != frontier.end(); it++) 
    it->second->findSampleCell();

  // Find min and max cell occupation
  //
  unsigned min1=MAXINT, max1=0, nt;
  for (key_cell::iterator it=frontier.begin(); it != frontier.end(); it++) {
    nt = it->second->bods.size();
    min1 = min<unsigned>(min1, nt);
    max1 = max<unsigned>(max1, nt);
  }

  MPI_Allreduce(&min1, &min_cell, 1, MPI_UNSIGNED, MPI_MIN, MPI_COMM_WORLD);
  MPI_Allreduce(&max1, &max_cell, 1, MPI_UNSIGNED, MPI_MAX, MPI_COMM_WORLD);

  vector<unsigned>  chist1(max_cell-min_cell+1, 0);
  chist = vector<unsigned>(max_cell-min_cell+1, 0);
  for (key_cell::iterator it=frontier.begin(); it != frontier.end(); it++)
    chist1[it->second->bods.size()-min_cell]++;
  
  MPI_Allreduce(&chist1[0], &chist[0], max_cell-min_cell+1, 
		MPI_UNSIGNED, MPI_SUM, MPI_COMM_WORLD);

  for (unsigned i=1; i<max_cell-min_cell+1; i++) chist[i] += chist[i-1];

  // Accumulate the total number of cells in the tree
  //
  unsigned my_cells = frontier.size();
  MPI_Allreduce(&my_cells, &total_cells, 1, MPI_UNSIGNED, MPI_SUM, MPI_COMM_WORLD);

#ifdef USE_GPTL
  GPTLstop("pHOT::makeTree::getFrontier");
  GPTLstop("pHOT::makeTree");
#endif

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
    maslev[lev] += p->state[0];
#ifdef INT128
    vollev[lev] += volume/(key_type(1u) << (3*p->level)).toDouble();
#else
    vollev[lev] += volume/static_cast<double>(key_type(1u) << (3*p->level));
#endif
  } else {
    map<unsigned, pCell*>::iterator it;
    for (it=p->children.begin(); it!=p->children.end(); it++) 
      densEmit(lev, it->second);
  }
}

void pHOT::densCheck()
{
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
    (*barrier)("pHOT: density check report");
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

}

void pHOT::dumpFrontier()
{
  key_cell::iterator it;
  unsigned sum = 0, cnt=0;
  double mean=0.0, disp=0.0;
  double totmass=0.0, totvol=0.0, tmp;

  if (myid==0) cout << endl << "Frontier info: " << endl;

  for (int n=0; n<numprocs; n++) {
    if (n==myid) {
      for (it=frontier.begin(); it!=frontier.end(); it++) {
	vector<double> mpos(3,0.0), vpos(3,0.0);
	vector<unsigned long>::iterator ib = it->second->bods.begin();
	double mass = 0.0;
	unsigned num = 0;
	while (ib != it->second->bods.end()) {
	  mass += cc->particles[*ib].mass;
	  for (unsigned k=0; k<3; k++) {
	    tmp = cc->particles[*ib].pos[k];
	    mpos[k] += tmp;
	    vpos[k] += tmp*tmp;
	  }
	  ib++;
	  num++;
	}
	
	totmass += mass;
#ifdef INT128
	totvol  += volume/(key_type(1u)<<(3*it->second->level)).toDouble();
#else
	totvol  += volume/static_cast<double>(key_type(1u)<<(3*it->second->level));
#endif

	cout << setw(4)  << myid
#ifdef INT128
	     << setw(12) << it->first.toHex()
#else
	     << setw(12) << hex << it->first << dec
#endif
	     << setw(8)  << it->second->level
	     << setw(18) << num
	     << setw(18) << mass
#ifdef INT128
	     << setw(18) << mass/(volume/(key_type(1u)<<(3*it->second->level))).toDouble();
#else
	     << setw(18) << mass/(volume/static_cast<double>(key_type(1u)<<(3*it->second->level)));
#endif
	
	for (unsigned k=0; k<3; k++) {
	  mpos[k] /= num;
	  if (num>1)
	    vpos[k] = sqrt( (vpos[k] - mpos[k]*mpos[k]*num)/(num-1) );
	  else
	    vpos[k] = 0.0;

	  cout << setprecision(4) << setw(10) << mpos[k] 
	       << setprecision(4) << setw(10) << vpos[k];
	}
	cout << endl;
	mean += num;
	disp += num*num;
	sum +=  num;
	cnt++;
      }
    }
    (*barrier)("pHOT: dump frontier");
  }

  unsigned sum0=0, cnt0=0;
  double mean0=0.0, disp0=0.0, totmass0=0.0, totvol0=0.0;

  MPI_Reduce(&sum, &sum0, 1, MPI_UNSIGNED, MPI_SUM, 0, MPI_COMM_WORLD);
  MPI_Reduce(&cnt, &cnt0, 1, MPI_UNSIGNED, MPI_SUM, 0, MPI_COMM_WORLD);
  MPI_Reduce(&mean, &mean0, 1, MPI_DOUBLE, MPI_SUM, 0, MPI_COMM_WORLD);
  MPI_Reduce(&disp, &disp0, 1, MPI_DOUBLE, MPI_SUM, 0, MPI_COMM_WORLD);
  MPI_Reduce(&totmass, &totmass0, 1, MPI_DOUBLE, MPI_SUM, 0, MPI_COMM_WORLD);
  MPI_Reduce(&totvol, &totvol0, 1, MPI_DOUBLE, MPI_SUM, 0, MPI_COMM_WORLD);

  if (myid==0) {
    cout << endl << setw(12) << "Total" << setw(18) << sum0 << endl;
    cout << endl << setw(12) << "Mass" << setw(18) << totmass0 << endl;
    cout << endl << setw(12) << "Volume" << setw(18) << totvol0 << endl;
    cout << setw(12) << "Mean" << setw(18) << mean0/cnt0 << endl;
    cout << setw(12) << "Sigma" << setw(18) 
	 << sqrt((disp0 - mean0*mean0/cnt0)/cnt0) << endl << endl;
  }

}

void pHOT::statFrontier()
{
  key_cell::iterator it;
  unsigned sum1=0, cnt1=0, sum=0, cnt=0, num;
  double mean1=0.0, disp1=0.0, mean=0.0, disp=0.0;
  vector<unsigned> freq1(pCell::bucket+1, 0), freq(pCell::bucket+1, 0);
  
  for (it=frontier.begin(); it!=frontier.end(); it++) {
    num = it->second->bods.size();
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

}


void pHOT::testFrontier(string& filename)
{
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


  for (int n=0; n<numprocs; n++) {

    if (n == myid) {

      ofstream out(filename.c_str(), ios::app);

      for (c=frontier.begin(); c!=frontier.end(); c++) {
	double mass=0, temp=0, pos[]={0,0,0}, vel[]={0,0,0};
	p = c->second;
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
#ifdef INT128
	double vol = volume/( key_type(1u) << (3*p->level)).toDouble(); 
#else
	double vol = volume/static_cast<double>( key_type(1u) << (3*p->level));
#endif

	out << " ";
	out << setw(prec[n]) << setiosflags(fmt[n]|typ[n]) << tnow;
	n++;
	out << setw(prec[n]) << setiosflags(fmt[n]|typ[n]) << myid;
	n++;
#ifdef INT128
	out << setw(prec[n]) << setiosflags(fmt[n]|typ[n]) << c->first.toHex();
#else
	out << setw(prec[n]) << setiosflags(fmt[n]|typ[n]) << hex << c->first << dec;
#endif
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

    (*barrier)("HOT: test frontier");
  }
}


void pHOT::countFrontier(vector<unsigned>& ncells, vector<unsigned>& bodies)
{
  pCell *p;
  key_cell::iterator it;
  map<unsigned, pair<unsigned, unsigned> > data;
  map<unsigned, pair<unsigned, unsigned> >::iterator d;
  unsigned maxLev1=0, maxLev=0;
  
  for (it=frontier.begin(); it != frontier.end(); it++) {
    p = it->second;
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

  for (d=data.begin(); d!=data.end(); d++) {
    tcellcnts[d->first] = d->second.first;
    tbodscnts[d->first] = d->second.second;
  }

  MPI_Reduce(&tcellcnts[0], &ncells[0], maxLev+1, MPI_UNSIGNED, MPI_SUM, 
	     0, MPI_COMM_WORLD);

  MPI_Reduce(&tbodscnts[0], &bodies[0], maxLev+1, MPI_UNSIGNED, MPI_SUM,
	     0, MPI_COMM_WORLD);
}


void pHOT::sendCell(key_type key, int to, unsigned num)
{
#ifdef USE_GPTL
  GPTLstart("pHOT::sendCell");
#endif

  pCell *p = frontier.find(key)->second;
  
#ifdef DEBUG
    cout << "Process " << myid << ": sending " << num 
	 << " to " << to << endl;

  vector<unsigned long> erased;
#endif

  vector<double> buffer1(3*num);
  vector<unsigned> buffer2(num);
  vector<key_type> buffer3(num);

  pf.ShipParticles(to, myid, num);

  key_pair tpair;
  vector<unsigned long>::iterator ib = p->bods.begin();
  for (unsigned j=0; j<num; j++) {

    pf.SendParticle(cc->particles[*ib]);
    
				// Find the record and delete it
    tpair.first  = cc->particles[*ib].key;
    tpair.second = cc->particles[*ib].indx;
    key_indx::iterator it = keybods.find(tpair);

    if (it != keybods.end()) {
				// Remove the key from the cell list
      key_key::iterator ij = bodycell.find(it->first);
      if (ij != bodycell.end()) bodycell.erase(ij);
      cc->particles.erase(*ib);
      keybods.erase(it);
#ifdef DEBUG
      erased.push_back(*ib);
#endif
    } else {
      cerr << "Process " << myid << ": error! " << endl;
    }
    ib++;
  }
  
  // If this cell is not the root
  //
  if (p->parent) {
    // Delete this cell from the parent
#ifdef INT128
    p->parent->children.erase( (p->mykey & 0x7u).toUint());
#else
    p->parent->children.erase( (p->mykey & 0x7u) );
#endif
	
#ifdef DEBUG
    if (frontier.find(p->mykey)==frontier.end()) {
      cout << "Process " << myid << ": in pHOT:sendCell: "
	   << " key not on frontier as expected" << endl;
    }
#endif

    // Delete this cell from the frontier
    frontier.erase(p->mykey);

    // Delete the cell altogether
    delete p;

  } else {			// Special treatment for root
    p->keys.clear();
    p->bods.clear();
  }

#ifdef DEBUG
  vector<unsigned long>::iterator iq;
  for (iq=erased.begin(); iq!=erased.end(); iq++) {
    if (cc->particles.find(*iq) != cc->particles.end())
      cout << "pHOT::sendCell proc=" << myid
	   << " found erased index=" << *iq << endl;
  }
#endif

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

#ifdef DEBUG
    cout << "Process " << myid << ": receiving " << num 
	 << " from " << from << endl;
#endif
  
  Particle part;

  pCell *p = root;

  pf.ShipParticles(myid, from, num);

  for (unsigned j=0; j<num; j++) {
    pf.RecvParticle(part);
      if (part.mass<=0.0 || isnan(part.mass)) {
	cout << "[recvCell crazy mass indx=" << part.indx 
#ifdef INT128
	     << ", key=" << part.key.toHex() << "]"
#else
	     << ", key=" << hex << part.key << dec << "]"
#endif
	  ;
      }
    cc->particles[part.indx] = part;
    if (part.key == 0u) continue;
    if (part.key < key_min || part.key >= key_max) {
      cout << "Process " << myid << ": in recvCell, key=" 
#ifdef INT128
	   << part.key.toHex() << "]"
#else
	   << hex << part.key << dec << "]"
#endif
	;
    }
    if (part.indx==0) cout << "pHOT::recvCell bad particle indx=0!" << endl;
    key_pair tpair(part.key, part.indx);
    keybods.insert(tpair);
    p = p->Add(tpair);
  }

  cc->nbodies = cc->particles.size();
#ifdef USE_GPTL
  GPTLstop("pHOT::recvCell");
#endif
}

void pHOT::makeState()
{
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
      
#ifdef INT128
      dens = state[0] * (key_type(1u) << (3*clv)).toDouble()/(volume*cnt);
#else
      dens = state[0] * static_cast<double>(key_type(1u) << (3*clv))/(volume*cnt);
#endif
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
  key_cell::iterator it;
  for (it=frontier.begin(); it!=frontier.end(); it++)
    MaxLev = max<unsigned>(MaxLev, it->second->level);

  double vol1, vol;
#ifdef INT128
  vol1 = volume/(key_type(1u) << (3*MaxLev)).toDouble();
#else
  vol1 = volume/static_cast<double>(key_type(1u) << (3*MaxLev));
#endif
  MPI_Allreduce(&vol1, &vol, 1, MPI_DOUBLE, MPI_MIN, MPI_COMM_WORLD);

  return vol;
}

double pHOT::maxVol()
{
  unsigned MinLev = MAXINT;
  key_cell::iterator it;
  for (it=frontier.begin(); it!=frontier.end(); it++)
    MinLev = min<unsigned>(MinLev, it->second->level);

  double vol1, vol;
#ifdef INT128
  vol1 = volume/(key_type(1u) << (3*MinLev)).toDouble();
#else
  vol1 = volume/static_cast<double>(key_type(1u) << (3*MinLev));
#endif
  MPI_Allreduce(&vol1, &vol, 1, MPI_DOUBLE, MPI_MAX, MPI_COMM_WORLD);

  return vol;
}

double pHOT::medianVol()
{
  unsigned mlev, num;
  vector<unsigned> lev;

  key_cell::iterator it;
  for (it=frontier.begin(); it!=frontier.end(); it++) 
    lev.push_back(it->second->level);

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

#ifdef INT128
  return volume/(key_type(1u) << (3*mlev)).toDouble();
#else
  return volume/static_cast<double>(key_type(1u) << (3*mlev));
#endif
}

void pHOT::Repartition(unsigned mlevel)
{
#ifdef USE_GPTL
  GPTLstart("pHOT::Repartition");
  GPTLstart("pHOT::Repartition::entrance_waiting");
  (*barrier)("pHOT: repartition entrance wait");
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

  oob.clear();
  vector<key_wght> keys;
  for (it=cc->Particles().begin(); it!=cc->Particles().end(); it++) {

    it->second.key = getKey(&(it->second.pos[0]));
    if (it->second.key == 0u) {
      oob.insert(it->first);
    } else {
      if (use_weight) {
	keys.push_back(key_wght(it->second.key, it->second.effort));
	// Reset effort value
	it->second.effort = Particle::effort_default;
      } else {
	keys.push_back(key_wght(it->second.key, 1.0));
      }
    }
  }

  // Check for duplicate keys, duplicate sequence numbers
  if (keys_debug) {
    key_indx keylist;
    set<key_type, less<key_type> > keydup;
    set<indx_type> indxdup;

#ifdef USE_GPTL
    GPTLstart("pHOT::Repartition::duplicate_check");
#endif

    for (it=cc->Particles().begin(); it!=cc->Particles().end(); it++)
      keylist.insert(pair<key_type, indx_type>(it->second.key, it->first));

    vector<indx_type> dups;
    key_indx::iterator k, kl;
    k = kl = keylist.begin();

    for (k++; k!=keylist.end(); k++) {
      if (k->first == kl->first) {
	dups.push_back(k->second);
	indxdup.insert(k->second);
	keydup. insert(k->first);
      }
      kl = k;
    }

    if (dups.size()) {
      ostringstream sout;
      sout << runtag << ".pHOT_crazy." << myid;
      ofstream out(sout.str().c_str(), ios::app);
      out << endl
	  << "#" << setfill('-') << setw(60) << '-' << setfill(' ') << endl
	  << "#---- Time=" << tnow 
	  << ", N=" << cc->Number() 
	  << ",  Duplicates=" << dups.size() 
	  << ", Unique keys=" << keydup.size() 
	  << ", Unique indx=" << indxdup.size() 
	  << endl
	  << "#" << setfill('-') << setw(60) << '-' << setfill(' ') << endl;
      for (vector<indx_type>::iterator 
	     ib=dups.begin(); ib!=dups.end(); ib++) {
	out << setw(10) << *ib << setw(18) << cc->Mass(*ib);
	for (int k=0; k<3; k++) out << setw(18) << cc->Pos(*ib, k);
	out << setw(10) << cc->Part(*ib)->indx
#ifdef INT128
	    << "    "   << cc->Part(*ib)->key.toHex() << endl;
#else
	<< "    "   << hex << cc->Part(*ib)->key << dec << endl;
#endif
      }
      out << "#" << setfill('-') << setw(60) << '-' << setfill(' ') << endl;
    }
#ifdef USE_GPTL
    GPTLstop("pHOT::Repartition::duplicate_check");
#endif
  }


#ifdef DEBUG
    cout << "Process " << myid 
	 << ": part #=" << cc->Particles().size()
	 << "  key size=" << keys.size()
	 << "  oob size=" << oob.size() << endl;
#endif
  
#ifdef USE_GPTL
  GPTLstop ("pHOT::Repartition::compute_keys");
  GPTLstart("pHOT::Repartition::compute_keys_waiting");
  (*barrier)("pHOT: repartition key wait");
  GPTLstop ("pHOT::Repartition::compute_keys_waiting");
  GPTLstart("pHOT::Repartition::spreadOOB");
#endif

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
    if (it->second.key == 0u) continue;
				// Look for key in this node's list
    t = find_proc(loclist, it->second.key);
    if (t == numprocs) {
      cerr << "Process " << myid << ": loclist found last entry, "
#ifdef INT128
	   << " key=" << it->second.key.toHex()
#else
	   << " key=" << hex << it->second.key << dec
#endif
	;
      
      cerr << ", end pt="
#ifdef INT128
	   << loclist.back().toHex()
#else
	   << hex << loclist.back() << dec
#endif
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

				// DEBUG OUTPUT
  if (false) {			// Set to "true" to enable

    if (myid==0){
      cout <<"--------------------------------------------------------"<< endl
	   <<"---- Send and receive counts for each process "<< endl
	   <<"--------------------------------------------------------"<< endl;
    }

    for (unsigned k=0; k<numprocs; k++) {
      if (myid==k) {
	cout << "Process " << k << endl;
	for (unsigned m=0; m<numprocs; m++) {
	  cout << setw(4) << m << setw(15) << sendcounts[m]
	       << setw(15) << recvcounts[m] << endl;
	}
      }
      (*barrier)("pHOT: repartition send/receive report");
    }
  }
				// END DEBUG OUTPUT
      
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
      (*barrier)("pHOT: repartition send/receive counts");
    }
  }


				// END DEBUG OUTPUT

  //
  // Exchange particles between processes
  //

  int ps;
  vector<Partstruct> psend(Tcnt), precv(Fcnt);

  timer_convert.start();
  for (int toID=0; toID<numprocs; toID++) {
    ps = sdispls[toID];
    for (unsigned i=0; i<sendcounts[toID]; i++) {
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
	cout << "[Repartition: crazy mass indx=" << part.indx 
	     << ", key="
#ifdef INT128
	     << part.key.toHex()
#else
	     << hex << part.key << dec
#endif
	<< "]";
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
  keybods.clear();
  unsigned oob1_cnt=0, oob_cnt=0;
  for (PartMapItr n=cc->Particles().begin(); n!=cc->Particles().end(); n++) {
    if (n->second.key==0u) {
      oob1_cnt++;
      continue;
    }
    if (n->second.indx==0)
      cout << "pHOT::Repartition bad particle indx=0!" << endl;

    keybods.insert(key_pair(n->second.key, n->second.indx));
  }

  // checkBounds(2.0, "AFTER repartition");

  MPI_Reduce(&oob1_cnt, &oob_cnt, 1, MPI_UNSIGNED, MPI_SUM, 0, MPI_COMM_WORLD);
  unsigned oob_tot = oobNumber();
  if (myid==0 && oob_cnt != oob_tot)
    cout << endl << "pHOT::Repartition: " << oob_cnt << " out of bounds," 
	 << " expected " << oob_tot << endl;

  if (false) {			// Sanity checks for bad particle indices
    bool ok = true;		// and bad particle counts
    for (PartMapItr ip=cc->particles.begin(); ip!=cc->particles.end(); ip++) {
      if (ip->second.indx==0) {
	cout << "pHOT::Repartition BAD particle in proc=" << myid
	     << " key="
#ifdef INT128
	     << ip->second.key.toHex()
#else
	     << hex << ip->second.key << dec
#endif
	     << endl;
	ok = false;
      }
    }
    if (ok) cout << "pHOT::Repartition: leaving with good indx" << endl;

    if (cc->particles.size() != cc->nbodies) 
      cout << "pHOT::Repartition: leaving with # mismatch!" << endl;

    int nbodies1 = cc->nbodies, nbodies0=0;
    MPI_Reduce(&nbodies1, &nbodies0, 1, MPI_INT, MPI_SUM, 0, MPI_COMM_WORLD);
    if (myid==0) {
      if (nbodies0 != cc->nbodies_tot)
	cout << "pHOT::Repartition: leaving with total # mismatch!" << endl;
    }

  }

  if (false) {
    checkIndices();
    if (!checkKeybods())
      cout << "Process " << myid 
	   << ": particle key not in keybods list at T=" << tnow << endl;
  }

  if (keys_debug) {		// Summary/diaganostic  output
				//
    unsigned int nsiz = cc->Particles().size();
    unsigned int ksiz = keys.size();
    vector<unsigned> nsize(numprocs), ksize(numprocs);

    MPI_Gather(&nsiz, 1, MPI_UNSIGNED, &nsize[0], 1, MPI_UNSIGNED, 
	       0, MPI_COMM_WORLD);

    MPI_Gather(&ksiz, 1, MPI_UNSIGNED, &ksize[0], 1, MPI_UNSIGNED, 
	       0, MPI_COMM_WORLD);
    
    if (myid==0) {
      ofstream out(string(runtag + ".pHOT_debug").c_str(), ios::app);
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
#ifdef INT128
	out << setw(klen) << kbeg[i].toHex();
	out << setw(klen) << kfin[i].toHex();
#else
	out << setw(klen) << hex << kbeg[i];
	out << setw(klen) << hex << kfin[i];
#endif
	out << dec << setw(15) << ksize[i] << setw(15) << nsize[i] << endl;
      }
      out << setfill('-') << setw(nhead) << '-' << endl << setfill(' ') << left
	  << endl << endl;
    }
  }

  timer_repartn.stop();

#ifdef USE_GPTL
  GPTLstop("pHOT::exchange");
  GPTLstop("pHOT::Repartition");
#endif
}

void pHOT::checkCellLevelList(const char *msg)
{
  unsigned missing_frontier_cell = 0;
  unsigned missing_clevlst_cell  = 0;

  for (key_cell::iterator 
	 it=frontier.begin(); it != frontier.end(); it++) 
    {
      if (it->second->mykey==1u && it->second->count==0u) {
	cout << "Process " << myid 
	     << ": empty root node in checkCellLevelList" << endl;
	continue;
      }
      if (clevlst.find(it->second) == clevlst.end()) 
	missing_frontier_cell++;
    }

  for (map<pCell*, unsigned>::iterator 
	 it=clevlst.begin(); it != clevlst.end(); it++) 
    {
      if (frontier.find(it->first->mykey) == frontier.end())
	missing_clevlst_cell++;
    }

  if (missing_frontier_cell)
    cout << "Process " << myid << ": " << msg << ", "
	 << missing_frontier_cell
	 << " frontier cells not in level list" << endl;

  if (missing_clevlst_cell)
    cout << "Process " << myid << ": " << msg << ", "
	 << missing_clevlst_cell
	 << " level list cells not on frontier" << endl;

}


void pHOT::makeCellLevelList()
{
  ostringstream sout;
  sout << "pHOT_cells." << runtag << "." << myid;
  ofstream out(sout.str().c_str(), ios::out | ios::app);

				// Make new lists
  clevlst.clear();
  clevels = vector< set<pCell*> >(multistep+1);

  out << "Process " << myid << " in makeCellLevelList()" 
      << ", frontier size=" << frontier.size() << endl;

  unsigned ng=0, nt=0;
  for (key_cell::iterator it=frontier.begin(); it != frontier.end(); it++) {
				// Check for empty root node
    if (it->second->mykey==1u && it->second->count==0u) {
      out << "Process " << myid << " in makeCellLevelList()" 
	  << ", empty root node" << endl;

      continue;
    }

    nt++;			// Otherwise, count this one
    it->second->remake_plev();
    clevlst[it->second] = it->second->maxplev;
    clevels[it->second->maxplev].insert(it->second);
    if (it->second->bods.size() == 0) {
      cerr << "Process " << myid 
	   << ": makeCellLevelList has a broken frontier!\n";
    } else {
      ng++;
    }
  }
  if (nt!=ng)
    out << "Process " << myid << ": made level list with " << ng
	<< " good cells out of " << nt << " expected" << endl;

#ifdef DEBUG
  printCellLevelList(out);
  checkParticles(out);
#endif
}

void pHOT::printCellLevelList(ostream& out)
{
  vector<unsigned> pcnt(multistep+1, 0);
  for (map<pCell*, unsigned>::iterator
	 pit=clevlst.begin(); pit!=clevlst.end(); pit++) pcnt[pit->second]++;

  out << left << setw(60) << setfill('-') << "-" << endl << setfill(' ')
      << "*** T=" << tnow << "  N=" << cc->particles.size() << endl
      << setw(10) << "M" << setw(10) << "number" 
      << setw(10) << "counts" << endl;
  for (unsigned M=0; M<=multistep; M++)
    out << setw(10) << M << setw(10) << clevels[M].size() 
	<< setw(10) << pcnt[M] << endl;
  out << left << setw(60) << setfill('-') << "-" << endl << setfill(' ');
}

void pHOT::adjustCellLevelList(unsigned mlevel)
{
  if (multistep==0) return;	// No need to bother if multistepping is off
				// Otherwise . . . 

  ostringstream sout;
  sout << "pHOT_cells." << runtag << "." << myid;
  ofstream out(sout.str().c_str(), ios::out | ios::app);

#ifdef USE_GPTL
  GPTLstart("pHOT::adjustCellLevelList");
#endif

  unsigned ng=0, nt=0, ns=0, m, cnt;
  for (unsigned M=mlevel; M<=multistep; M++) {
    nt += clevels[M].size();
    cnt = 0;
    if (clevels[M].size()>0) {
      set<pCell*>::iterator it = clevels[M].begin(), nit;
      while (it != clevels[M].end()) {
				// Skip the root cell if it's empty
				// (lazy kludge)
	if ( (*it)->mykey==1u && (*it)->count==0 ) { 
	  nt--; 		// Reduce the node count by one
	  it++;			// Go the next set in the set . . .
	  continue; 
	}

	cnt++;			// Count the (presumably good) cells
	
				// For diagnostic info only
	if ((*it)->bods.size()) ng++;
	else {			// This shouldn't happen
	  cout << "Process " << myid << ": pHOT::adjustCellLevelList: "
	       << cnt << "/" << clevels[M].size()
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
	
	if (clevels[M].empty()) break;

      }
    }
  }

				// Diagnostic output . . .
  if (nt!=ng)
    cout << "Process " << myid << ": adjusted level list with " << ng
	 << " good cells out of " << nt << " expected, " << ns
	 << " cells moved" << endl;
#ifdef DEBUG
  printCellLevelList(out);
  checkParticles(out);
#endif

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
  (*barrier)("pHOT: repartition key timer [0]");
  timer_waiton0.stop();

  timer_tadjust.start();

#ifdef DEBUG_ADJUST
  if (!checkBodycell()) {
    cout << "Process " << myid << ": "
	 << "adjustTree: ERROR bodycell BEFORE adjustTree(), T="
	 << tnow << " mlevel=" << mlevel << endl;
  }    
  if (!checkPartKeybods(mlevel)) {
    cout << "Process " << myid 
	 << ": adjustTree: ERROR particle/keybods BEFORE adjustTree(), T=" 
	 << tnow << " mlevel=" << mlevel << endl;
  }
  
  if (!checkKeybods()) {
    cout << "Process " << myid 
	 << ": adjustTree: ERROR particle key not in keybods"
	 << " BEFORE adjustTree(), T=" << tnow << endl;
  }

  checkCellLevelList("BEFORE adjustTree()");
#endif

  adjcnt++;			// For debug labeling only . . .
  
  timer_keymake.start();

  pCell* c;
  key_type newkey, oldkey;
  list<unsigned long>::iterator ip;
  list<unsigned long> oldp;

  //
  // Exchange list
  //
  timer_keybods.start();
  vector< vector<unsigned long> > exchange(numprocs);

  // 
  // Make body list from frontier cells for this level
  //
  for (unsigned M=mlevel; M<=multistep; M++) {
				// For speed efficiency (I hope)
    set<pCell*>::iterator it    =clevels[M].begin();
    set<pCell*>::iterator itend =clevels[M].end();
    for (;it!=itend; it++) {
      vector<unsigned long>::iterator ib    =(*it)->bods.begin(); 
      vector<unsigned long>::iterator ibend =(*it)->bods.end(); 
      for (; ib!=ibend; ib++) {
	oldp.push_back(*ib);
      }
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
  (*barrier)("pHOT: repartition key timer [1]");
  timer_waiton1.stop();

  //
  // Update body by body using the list without regard to level
  //
  unsigned newproc;
  for (ip=oldp.begin(); ip!=oldp.end(); ip++) {

    timer_keybods.start();
    Particle *p = cc->Part(*ip);
    if (p==0) {			// Sanity check
      cout << "Process " << myid 
	   << ": pHOT::adjustTree: ERROR crazy particle index!" << endl;
    }
    numkeys++;
    timer_keybods.stop();

    //
    // Get and recompute keys
    //
    timer_keycall.start();
    oldkey = p->key;
    newkey = getKey(&(p->pos[0]));
    timer_keycall.stop();

    //
    // Get this particle's cell
    //

    timer_keybods.start();
    key_key ::iterator ij = bodycell.find(oldkey);

				// Bad key?
    if (ij == bodycell.end()) {	//
      cout << "Process " << myid 
	   << ": pHOT::adjustTree: ERROR could not find cell for particle"
#ifdef INT128
	   << " key=" << oldkey.toHex() << ", index=" << p->indx
#else
	   << " key=" << hex << oldkey << dec << ", index=" << p->indx
#endif
	   << " pnumber=" << cc->Number() << " bodycell=" << bodycell.size() 
	   << endl;
    } else {
				// Bad cell?
      if (0) {			//
	if ( frontier.find(ij->second) == frontier.end() ) {
	  cout << "Process " << myid 
	       << ": pHOT::adjustTree: ERROR could not find expected cell"
	       << " on frontier, count=" << adjcnt;
#ifdef INT128
	  cout << " oldbody=" << oldkey.toHex();
	  cout << " newbody=" << newkey.toHex();
	  cout << " cell="    << bodycell.find(oldkey)->second.toHex()
#else
	  cout << " oldbody=" << hex << oldkey << dec;
	  cout << " newbody=" << hex << newkey << dec;
	  cout << " cell="    << hex << bodycell.find(oldkey)->second<< dec
#endif
	       << " index="   << p->indx 
	       << endl;
	  continue;
	}
      }
      
      c = frontier[ij->second];
    }
    timer_keybods.stop();
      
    cntr_total++;

    //
    // Is the key the same?
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
	    exchange[newproc].push_back(*ip);

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

#ifdef DEBUG
	if (ij == bodycell.end()) {
	  cout << "Process " << myid 
	       << ": pHOT::adjustTree: ERROR could not find cell for"
	       << " key=" << hex << oldkey << ", index=" << dec << p->indx;
	  if (!c->isMine(oldkey))
	    cout << ", this IS NOT my key!";
	  cout << endl;
	}
#endif
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
  (*barrier)("pHOT: repartition key timer [2]");
  timer_waiton2.stop();

  timer_cupdate.start();

#ifdef USE_GPTL
  GPTLstop ("pHOT::keyCells");
  GPTLstart("pHOT::adjExchange");
#endif


  //
  // Exchange particles
  //

  Particle part;
  unsigned Tcnt=0, Fcnt, sum;
  vector<int> sdispls(numprocs), rdispls(numprocs);
  vector<int> sendcounts(numprocs, 0), recvcounts(numprocs, 0);

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
    vector<Partstruct> psend(Tcnt), precv(Fcnt);
    
    timer_convert.start();
    for (int toID=0; toID<numprocs; toID++) {
      ps = sdispls[toID];
      for (unsigned i=0; i<sendcounts[toID]; i++) {
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
	cout << "[adjustTree: crazy mass indx=" << part.indx 
	     << ", key="
#ifdef INT128
	     << part.key.toHex()
#else
	     << hex << part.key<< dec
#endif
	     << "]";
      }
      
      cc->Particles()[part.indx] = part;
      
      if (part.key != 0u) {
	key_pair newpair(part.key, part.indx);
	keybods.insert(newpair);
	root->Add(newpair, &change);
      }
    }
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


  //
  // Cell overlap?
  //
  
  key_type headKey=0u, tailKey=0u;
  unsigned head_num=0, tail_num=0;
  
  timer_overlap.start();
  
  for (int n=0; n<numprocs-1; n++) {
    
    if (n==myid) {
      if (keybods.size()) {
	key_indx::reverse_iterator it = keybods.rbegin();
	tailKey  = bodycell.find(it->first)->second;
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
	  vector<Partstruct> Precv(head_num);
	  MPI_Recv(&Precv[0], head_num, pf.Particletype, 
		   n+1, 136, MPI_COMM_WORLD, MPI_STATUS_IGNORE);
	  
	  for (int i=0; i<head_num; i++) {
	    pf.part_to_Particle(Precv[i], part);
	    cc->Particles()[part.indx] = part;
	    key_pair newpair(part.key, part.indx);
	    keybods.insert(newpair);
	    c->Add(newpair, &change);
	  }
	  
	} else {
	  //
	  // Send particles
	  //
	  unsigned k=0;
	  vector<Partstruct> Psend(tail_num);
	  vector<unsigned long>::iterator ib;
	  for (ib=c->bods.begin(); ib!=c->bods.end(); ib++) {
	    pf.Particle_to_part(Psend[k++], cc->Particles()[*ib]);
	    cc->Particles().erase(*ib);
	  }
	  
	  c->RemoveAll();
			
	  if (c->mykey!=1u) {	// Don't remove the root node
	    
	    // queue for removal from level lists
	    change.push_back(cell_indx(c, REMOVE));
	  
	    // queue for deletion
	    change.push_back(cell_indx(c, KILL));
	  }

	  MPI_Send(&Psend[0], tail_num, pf.Particletype, 
		   n+1, 135, MPI_COMM_WORLD);
	}
      }
    }
    
    if (n+1==myid) {

      if (keybods.size()) {
	key_indx::iterator it = keybods.begin(); // Sanity check:
	if (bodycell.find(it->first) == bodycell.end()) {
	  cerr << "In adjustTree: No cell for body=" 
#ifdef INT128
	       << it->first.toHex()
#else
	       << hex << it->first << dec
#endif
	       << " bodycell size=" << bodycell.size() << endl;
	  headKey  = 0u;
	  head_num = 0;
	} else {
	  headKey  = bodycell.find(it->first)->second;
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
	  vector<Partstruct> Psend(head_num);
	  vector<unsigned long>::iterator ib;
	  for (ib=c->bods.begin(); ib!=c->bods.end(); ib++) {
	    pf.Particle_to_part(Psend[k++], cc->Particles()[*ib]);
	    cc->Particles().erase(*ib);
	  }
	  
	  c->RemoveAll();
	  
	  if (c->mykey!=1u) { // Dont remove the root node!

	    // queue for removal from level lists
	    change.push_back(cell_indx(c, REMOVE));
	  
	    // queue for deletion
	    change.push_back(cell_indx(c, KILL));
	  }

	  MPI_Send(&Psend[0], head_num, pf.Particletype, 
		   n, 136, MPI_COMM_WORLD);
	  
	} else {		
	  //
	  // Receive particles
	  //
	  vector<Partstruct> Precv(tail_num);
	  MPI_Recv(&Precv[0], tail_num, pf.Particletype, 
		   n, 135, MPI_COMM_WORLD, MPI_STATUS_IGNORE);

	  for (int i=0; i<tail_num; i++) {
	    pf.part_to_Particle(Precv[i], part);
	    cc->Particles()[part.indx] = part;
	    key_pair newpair(part.key, part.indx);
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

  // timer_cupdate.start();

  // Create, remove and delete changed cells
  // ---------------------------------------
  // Need to do the creates first in case that cells enter and exit
  // the frontier during the Add() step
  //
  for (list<cell_indx>::iterator 
	 it=change.begin(); it!=change.end(); it++) {
      
    pCell *c = it->first;	// This is the cell pointer

    switch(it->second) {
    case CREATE:		// Add this one to the active lists
      {
	unsigned m = max<unsigned>(c->maxplev, mlevel);
	clevlst[c] = m;
	clevels[m].insert(c);
				// And locate the sample cell
	c->findSampleCell();
      }
      break;

    case REMOVE:		// Remove this cell from the active list
      {
#ifdef DEBUG
	if (clevlst.find(c) == clevlst.end()) {
	  cout << "pHOT::adjustTree: cell=" << hex << c
	       << dec << " not in level list";
	  if (frontier.find(c->mykey) == frontier.end())
	    cout << " and gone from frontier";
	  else
	    cout << " but is in frontier";
	  cout << endl;
	  continue;
	}
#endif
	unsigned m = clevlst[c];
	clevlst.erase(c);
#ifdef DEBUG
	if (clevels[m].find(c) == clevels[m].end()) {
	  cout << "pHOT::adjustTree: cell=" << hex << c
	       << dec << " not in level " << m << endl;
	}
#endif
	clevels[m].erase(c);
      }
      break;

    case KILL:			// Delete the cell
      delete c;
      break;

    case RECOMP:		// Relocate the sample cell (changes)
      c->findSampleCell();
      break;

    default:
      cout << "Process " << myid << ": unknown action in pHOT::adjustTree()!!"
	   << endl;
    }
  }

  change.clear();		// Reset the change list for next time
  
#ifdef DEBUG_ADJUST
  if (!checkBodycell()) {
    cout << "Process " << myid << ": "
	 << "adjustTree: ERROR mlevel=" << mlevel 
	 << " completed: body cell check FAILED!" << endl;
  }    
  if (!checkParticles(cout)) {
    cout << "Process " << myid << ": "
	 << "adjustTree: ERROR mlevel=" << mlevel 
	 << " completed: initial particle check FAILED!" << endl;
  }    
  if (!checkFrontier(cout)) {
    cout << "Process " << myid << ": "
	 << "adjustTree: ERROR mlevel=" << mlevel
	 << " completed: frontier check FAILED!" << endl;
  }
  
  checkDupes();
  checkIndices();
  checkCellLevelList("AFTER adjustTree()");
  
  if (!checkKeybods()) {
    cout << "Process " << myid 
	 << ": adjustTree: ERROR particle key not in keybods AFTER adjustTree(), T=" 
	 << tnow << endl;
  }
  
  if (!checkPartKeybods(mlevel)) {
    cout << "Process " << myid 
	 << ": adjustTree: ERROR particle/keybods AFTER adjustTree(), T=" 
	 << tnow << " mlevel=" << mlevel << endl;
  }
#endif
  
  timer_cupdate.stop();

  timer_tadjust.stop();

  if (keys_debug) {		// Summary/diaganostic  output
				//
    unsigned int nsiz = cc->Particles().size();
    vector<unsigned> nsize(numprocs);

    MPI_Gather(&nsiz, 1, MPI_UNSIGNED, &nsize[0], 1, MPI_UNSIGNED, 
	       0, MPI_COMM_WORLD);

    if (myid==0) {
      ofstream out(string(runtag + ".pHOT_debug").c_str(), ios::app);
      unsigned nhead = 20 + 2*klen;

      out << left << setfill('-') << setw(nhead) << '-' << endl
	  << "---- Post-adjustTree summary [" << mlevel << "], T = " 
	  << tnow << endl << setw(nhead) << '-' << endl << setfill(' ');
      out << left << setw(5) << "#" 
	  << setw(klen) << right << "kbeg" << setw(klen) << "kfin" 
	  << setw(15) << "bodies" << endl << "#" << endl;
      for (int i=0; i<numprocs; i++) {
	out << left << setw(5) << i << right
#ifdef INT128
	    << setw(klen) << kbeg[i].toHex()
	    << setw(klen) << kfin[i].toHex()
#else
	    << hex << setw(klen) << kbeg[i]
	    << setw(klen) << kfin[i] << dec
#endif
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
  

bool pHOT::checkParticles(ostream& out, bool pc)
{
  unsigned cnt, badb=0, badc=0;
  const string membername = "checkParticles";

  // FOR CELL OUTPUT
  map<pCell*, pair<unsigned, unsigned> > outdat;

  for (unsigned M=0; M<=multistep; M++) {
    if (clevels[M].size()) {
      cnt = 0;
      for (set<pCell*>::iterator 
	     it=clevels[M].begin(); it!=clevels[M].end(); it++) {
	cnt++;
	// FOR CELL OUTPUT
	outdat[*it] = pair<unsigned, unsigned>(M, (*it)->bods.size());

	if ((*it)->bods.size()) {
	  if (pc) {
	    for (vector<unsigned long>::iterator ib=(*it)->bods.begin();
		 ib!=(*it)->bods.end(); ib++) {
	      if (!cc->Part(*ib)) {
		out << "pHOT::checkParticles:: M=" << M << ", bad body at "
		     << cnt << "/" << clevels[M].size() 
		     << " cell=" << hex << (*it) << dec << endl;
		badb++;
	      }
	    }
	  }
	} else {
	  out << "pHOT::checkParticles:: M=" << M << ", zero bods at "
	       << cnt << "/" << clevels[M].size() 
	       << " cell=" << hex << (*it) << dec << endl;
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
    message << "error creating backup file <" << backfile1 << ">";
    // bomb(membername, message.str());
  }
  
  ostringstream origfile2, backfile2;

  origfile2 << "chklist." << myid;
  backfile2 << "chklist." << myid << ".bak";
  if (rename(origfile2.str().c_str(), backfile2.str().c_str())) {
    perror("pHOT");
    ostringstream message;
    message << "error creating backup file <" << backfile2 << ">";
    // bomb(membername, message.str());
  }
  
  ofstream out1(origfile1.str().c_str());
  ofstream out2(origfile2.str().c_str());

  if (!out1) {
    ostringstream message;
    message << "error opening output file <" << origfile1 << ">";
    bomb(membername, message.str());
  }

  if (!out2) {
    ostringstream message;
    message << "error opening output file <" << origfile2 << ">";
    bomb(membername, message.str());
  }

  for (map<pCell*, pair<unsigned, unsigned> >::iterator
	 lit=outdat.begin(); lit!=outdat.end(); lit++) {
    out1 << setw(15) << hex << lit->first << dec 
	 << setw(8) << lit->second.first
	 << setw(8) << lit->second.second
	 << endl;
  }

  for (map<pCell*, unsigned>::iterator
	 lit=clevlst.begin(); lit!=clevlst.end(); lit++) {
    out2 << setw(15) << hex << lit->first << dec 
	 << setw(8) << lit->second
	 << endl;
  }
  // END OUTPUT CELL LIST

  if (badb || badc) {
    out << "Process " << myid << ": pHOT::checkParticles, bad cell=" 
	 << badc ;
    if (pc) out << " bad bod=" << badb;
    out << endl;
    return false;
  }
  else return true;
}


bool pHOT::checkFrontier(ostream& out)
{
  unsigned bad=0;
  bool good=true;

  for (unsigned M=0; M<=multistep; M++) {
    if (clevels[M].size()) {
      for (set<pCell*>::iterator 
	     it=clevels[M].begin(); it!=clevels[M].end(); it++) {
	if (frontier.find((*it)->mykey) == frontier.end()) {
	  out << "pHOT::checkFrontier error on M=" << M
	      << ", cell=" << hex << *it << dec << endl;
	  bad++;
	  good = false;
	}
      }
    }
  }
  
  if (bad) {
    out << "Process " << myid << ": pHOT::checkFrontier, bad cell=" 
	 << bad << endl;
  }

  return good;
}


bool pHOT::checkDupes()
{
  vector<unsigned> indices;
  for (key_cell::iterator it=frontier.begin(); it!=frontier.end(); it++)
    indices.insert(indices.end(), it->second->bods.begin(), it->second->bods.end());

  sort(indices.begin(), indices.end());
  
  unsigned dup = 0;
  for (unsigned n=1; n<indices.size(); n++) {
    if (indices[n-1] == indices[n]) dup++;
  }

  if (dup) {
    cout << "Process " << myid << ": pHOT::checkDupes, dup=" << dup << endl;
    return false;
  } else true;
}


void pHOT::checkIndices()
{

  // All prcesses make an index list
  //
  vector<unsigned> indices;
  for (key_cell::iterator it=frontier.begin(); it!=frontier.end(); it++)
    indices.insert(indices.end(), it->second->bods.begin(), it->second->bods.end());

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
    (*barrier)("pHOT: check indicies");
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
}

bool pHOT::checkKeybods()
{
  bool ok = true;
  unsigned cnt=0;
  for (PartMapItr n = cc->Particles().begin(); n!=cc->Particles().end(); n++) {

				// Skip oob particles
    if (n->second.key) {
      key_pair tpair(n->second.key, n->second.indx);
      key_indx::iterator it = keybods.find(tpair);

      if (it==keybods.end()) {
#ifdef DEBUG
	cout << "Process " << myid << ": checkKeybods: " 
	     << cnt << " unmatched particle, (x, y, z)=("
	     << n->second.pos[0] << ", " << n->second.pos[1] 
	     << ", " << n->second.pos[2] << ")" << endl;
#endif
	ok = false;
	cnt++;
      }
    }
  }
  
#ifdef DEBUG
  if (cnt) {
    cout << "Process " << myid << ": checkKeybods: " 
	 << cnt << " unmatched particles" << endl;
  }
#endif

  return ok;
}


bool pHOT::checkPartKeybods(unsigned mlevel)
{
  unsigned bcelbod = 0;
  unsigned bfrontr = 0;
  bool ok = true;

  // 
  // Make body list from frontier cells for this level
  //
  for (unsigned M=mlevel; M<=multistep; M++) {
    for (set<pCell*>::iterator it=clevels[M].begin(); it!=clevels[M].end(); it++) {
      for (vector<unsigned long>::iterator ib=(*it)->bods.begin(); ib!=(*it)->bods.end(); ib++) {
	key_type key  = cc->Particles()[*ib].key;
	unsigned indx = cc->Particles()[*ib].indx;

	// Look for cell for this body
	//
	if (bodycell.find(key) == bodycell.end()) {
	  bcelbod++;
	  ok = false;
	}
	// Look for cell in frontier . . .
	//
	if (frontier.find(bodycell.find(key)->second) == frontier.end()) {
	  bfrontr++;
	  ok = false;
	}
      }
    }
  }
  //
  // Report
  //
  if (!ok) {
    cout << "Process " << myid 
	 << ": pHOT::checkPartKeybods: ERROR bad bodycell=" << bcelbod
	 << " bad frontier=" << bfrontr << endl;

  }

  return ok;
}


bool pHOT::checkBodycell()
{
  bool ok = true;
  unsigned cnt=0;
  for (PartMapItr n=cc->Particles().begin(); n!=cc->Particles().end(); n++) {
    // Ignore OOB particle
    if (n->second.key==0u) continue;

    // Look for the bodycell
    key_key::iterator it = bodycell.find(n->second.key);
    if (it==bodycell.end()) {
      ok = false;
      cnt++;
#ifdef DEBUG
      cout << "Process " << myid << ": checkBodycell: " 
	   << cnt << " unmatched particle: key=" << hex
	   << n->second.key << dec << " index=" 
	   << n->second.indx  << ", (x, y, z)=("
	   << n->second.pos[0] << ", " << n->second.pos[1] 
	   << ", " << n->second.pos[2] << ")" << endl;
#endif
    }
  }

#ifdef DEBUG
  if (cnt) {
    cout << "Process " << myid << ": checkBodycell: " 
	 << cnt << " unmatched particles" << endl;
  }
#endif

  return ok;
}


void pHOT::checkBounds(double rmax, const char *msg)
{
  int bad = 0;
  for (PartMapItr n = cc->Particles().begin(); n!=cc->Particles().end(); n++) {
    for (int k=0; k<3; k++) if (fabs(n->second.pos[k])>rmax) bad++;
  }
  if (bad) {
    cout << "Process " << myid << ": has " << bad << " out of bounds";
    if (msg) cout << ", " << msg << endl;
    else cout << endl;
  }
}

unsigned pHOT::oobNumber()
{
  unsigned number=0, number1=oob.size();
  MPI_Reduce(&number1, &number, 1, MPI_UNSIGNED, MPI_SUM, 0, MPI_COMM_WORLD);
  return number;
}

void pHOT::checkOOB(vector<unsigned>& sendlist)
{
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
  if (false) {
    ostringstream sout;
    sout << "In spreadOOB, maxdif=" << maxdif << " #=" << cc->nbodies_tot
	 << " ns=" << nsend.size() << " rs=" << nrecv.size();
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

#ifdef DEBUG
  checkOOB(sendlist);
#endif

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
  if (false) {
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
  GPTLstart("pHOT::spreadOOB::exchange_particles");
#endif

  unsigned ps=0, pr=0;
  set<indx_type>::iterator ioob;
  Partstruct *psend=0, *precv=0;
  vector<MPI_Request> rql;
  MPI_Request r;
  int ierr;

  if (Tcnt) psend = new Partstruct [Tcnt];
  if (Fcnt) precv = new Partstruct [Fcnt];

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
  GPTLstop ("pHOT::spreadOOB::exchange_particles");
  GPTLstart("pHOT::spreadOOB::add_to_particles");
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
	     << ", key="
#ifdef INT128
	     << part.key.toHex()
#else
	     << hex << part.key << dec
#endif
	     << "]";
      }
      cc->Particles()[part.indx] = part;
      oob.insert(part.indx);
    }
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
  if (keys_debug && myid==0) {
    timer_debug = new Timer(true);
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
	if (static_cast<unsigned>(floor(srate*(i+1)-DBL_MIN)) == j) {
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
  // (set to keys_debug "true" at top of file to enable key range diagnostic)
  //
  if (keys_debug) {

    if (myid==0) {
	ofstream out(string(runtag + ".pHOT_debug").c_str(), ios::app);
	out << endl
	    << "-------------------------------------------" << endl
	    << "--- Sampling stats ------------------------" << endl
	    << "-------------------------------------------" << endl << left
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
	ofstream out(string(runtag + ".pHOT_debug").c_str(), ios::app);
	out << left
	    << setw(5)  << myid
	    << setw(10) << keys.size()
	    << setw(15) << srate
	    << setw(10) << keylist1.size()
	    << endl;
      }

      (*barrier)("pHOT: partitionKeys debug 1");
    }
    if (myid==0) {
      ofstream out(string(runtag + ".pHOT_debug").c_str(), ios::app);
      out << "-------------------------------------------" << endl << endl;
    }

    unsigned nhead = 15 + 5*klen;

    if (myid==0) {
      ofstream out(string(runtag + ".pHOT_debug").c_str(), ios::app);
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
	ofstream out(string(runtag + ".pHOT_debug").c_str(), ios::app);
	out << left << setw(5) << n << setw(10) << keys.size();
	if (keys.size()>0) {
	  out << right
#ifdef INT128
	      << setw(klen) << keys[0].first.toHex()
	      << setw(klen) << keys[keys.size()/4].first.toHex()
	      << setw(klen) << keys[keys.size()/2].first.toHex()
	      << setw(klen) << keys[keys.size()*3/4].first.toHex()
	      << setw(klen) << keys[keys.size()-1].first.toHex()
#else
	      << hex
	      << setw(klen) << keys[0].first
	      << setw(klen) << keys[keys.size()/4].first
	      << setw(klen) << keys[keys.size()/2].first
	      << setw(klen) << keys[keys.size()*3/4].first
	      << setw(klen) << keys[keys.size()-1].first << dec
#endif
	      << endl;
	}
      }
      (*barrier)("pHOT: partitionKeys debug 2");
    }

    if (myid==0) {
      ofstream out(string(runtag + ".pHOT_debug").c_str(), ios::app);
      out << left << setfill('-') << setw(nhead) << "-" << endl 
	  << setfill('-') << endl;
    }

    //
    // set to "true" to enable key list diagnostic
    //
    if (false) {
      const unsigned cols = 3;	// # of columns in output
      const unsigned cwid = 35;	// column width

      if (myid==0) {
	ofstream out(string(runtag + ".pHOT_debug").c_str(), ios::app);
	out << "--------------------------" << endl;
	out << "----- Partition keys -----" << endl;
	out << "--------------------------" << endl;

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
	  ofstream out(string(runtag + ".pHOT_debug").c_str(), ios::app);
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
	(*barrier)("pHOT: partitionKeys debug 3");
      }
    }
  }
  //
  // END DEBUG
  //
				// Tree aggregation (merge sort) of the
				// entire key list

				// <keylist1> is (and must be) sorted to start
  (*barrier)("pHOT: partitionKeys before pMerge");
  parallelMerge(keylist1, keylist);
  (*barrier)("pHOT: partitionKeys after pMerge");

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
      
    if (keys_debug) {	 // If true, print key ranges for each process
      ofstream out(string(runtag + ".pHOT_debug").c_str(), ios::app);
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
      for (int i=0; i<numprocs; i++) {
	out << left << setw(5) << i << right
#ifdef INT128
	    << setw(klen) << kbeg[i].toHex()
	    << setw(klen) << kfin[i].toHex()
	    << setw(klen) << (kfin[i] - kbeg[i]).toDec()
#else
	    << hex << setw(klen) << kbeg[i]
	    << setw(klen) << kfin[i] << dec
	    << setw(klen) << (kfin[i] - kbeg[i])
#endif
	    << setw(15) << wbeg[i]
	    << setw(15) << wfin[i]
	    << setw(15) << wfin[i] - wbeg[i]
	    << setw(10) << pbeg[i]
	    << setw(10) << pfin[i]
	    << setw(10) << pfin[i] - pbeg[i]
	    << endl;
      }
      out << setw(nhead) << setfill('-') << '-' << endl << setfill(' ') 
	  << endl;
    }
  }

  MPI_Bcast(&kbeg[0], numprocs, MPI_EXP_KEYTYPE, 0, MPI_COMM_WORLD);
  MPI_Bcast(&kfin[0], numprocs, MPI_EXP_KEYTYPE, 0, MPI_COMM_WORLD);


  if (keys_debug) {		// If true, print key totals
    unsigned oobn = oobNumber();
    unsigned tkey1 = keys.size(), tkey0 = 0;
    MPI_Reduce(&tkey1, &tkey0, 1, MPI_UNSIGNED, MPI_SUM, 0, MPI_COMM_WORLD);
    if (myid==0) {
      ofstream out(string(runtag + ".pHOT_debug").c_str(), ios::app);
      out << endl
	  << setfill('-') << setw(60) << '-' << setfill(' ') << endl
	  << "---- partitionKeys" << endl
	  << setfill('-') << setw(60) << '-' << setfill(' ') << endl
	  << "----  list size=" << setw(10) << keylist.size() << endl
	  << "---- total keys=" << setw(10) << tkey0 << endl
	  << "----  total oob=" << setw(10) << oobn << endl
	  << "----      TOTAL=" << setw(10) << tkey0 + oobn << endl
	  << setfill('-') << setw(60) << '-' << setfill(' ') << endl
	  << "----   Time (s)=" << setw(10) << 1.0e-6*timer_debug->stop().getTotalTime()
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
  unsigned beg=0, end=num-1, cur;

  if (num==0) {
    cerr << "pHOT::find_proc: crazy input, no keys!" << endl;
    return 0;
  }

  if (key<=keys[beg])   return beg;
  if (key>=keys[num-1]) return num-1;

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
  while (2*M2 < numprocs) M2 = M2*2;

  // Combine the particles of the high nodes
  // with those of the lower nodes so that
  // all particles are within M2 nodes
  //
  // NB: if M2 == numprocs, no particles
  // will be sent or received
  //
  if (myid >= M2) {
    n = initl.size();
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
    return;
  }

  vector<key_wght> data = initl;

  //
  // Retrieve the excess particles
  //
  if (myid + M2 < numprocs) {
    MPI_Recv(&n, 1, MPI_UNSIGNED, myid+M2, 11, MPI_COMM_WORLD, &status);
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
    }
  }
    
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
      return;

    } else {
      MPI_Recv(&n, 1, MPI_UNSIGNED, MPI_ANY_SOURCE, 11, MPI_COMM_WORLD, 
	       &status);
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

  //
  // We are done, return the result
  //

  final = data;

#ifdef USE_GPTL
  GPTLstop("pHOT::parallelMerge");
#endif

  return;
}


unsigned pHOT::checkNumber()
{
  unsigned nbods1=0, nbods=0;
  for (key_cell::iterator it=frontier.begin(); it!=frontier.end(); it++) 
    nbods1 += it->second->bods.size();
  
  MPI_Reduce(&nbods1, &nbods, 1, MPI_UNSIGNED, MPI_SUM, 0, MPI_COMM_WORLD);

  return nbods;
}

void pHOT::CollectTiming()
{
  float fval;

  fval  = timer_keymake.getTime().getRealTime()*1.0e-6;
  MPI_Gather(&fval, 1, MPI_FLOAT, &keymk3[0], 1, MPI_FLOAT, 0, MPI_COMM_WORLD);

  fval = timer_xchange.getTime().getRealTime()*1.0e-6;
  MPI_Gather(&fval, 1, MPI_FLOAT, &exchg3[0], 1, MPI_FLOAT, 0, MPI_COMM_WORLD);

  fval = timer_convert.getTime().getRealTime()*1.0e-6;
  MPI_Gather(&fval, 1, MPI_FLOAT, &cnvrt3[0], 1, MPI_FLOAT, 0, MPI_COMM_WORLD);

  fval  = timer_overlap.getTime().getRealTime()*1.0e-6;
  MPI_Gather(&fval, 1, MPI_FLOAT, &tovlp3[0], 1, MPI_FLOAT, 0, MPI_COMM_WORLD);

  fval  = timer_prepare.getTime().getRealTime()*1.0e-6;
  MPI_Gather(&fval, 1, MPI_FLOAT, &prepr3[0], 1, MPI_FLOAT, 0, MPI_COMM_WORLD);

  fval = timer_cupdate.getTime().getRealTime()*1.0e-6;
  MPI_Gather(&fval, 1, MPI_FLOAT, &updat3[0], 1, MPI_FLOAT, 0, MPI_COMM_WORLD);

  fval = timer_scatter.getTime().getRealTime()*1.0e-6;
  MPI_Gather(&fval, 1, MPI_FLOAT, &scatr3[0], 1, MPI_FLOAT, 0, MPI_COMM_WORLD);

  fval = timer_repartn.getTime().getRealTime()*1.0e-6;
  MPI_Gather(&fval, 1, MPI_FLOAT, &reprt3[0], 1, MPI_FLOAT, 0, MPI_COMM_WORLD);

  fval = timer_tadjust.getTime().getRealTime()*1.0e-6;
  MPI_Gather(&fval, 1, MPI_FLOAT, &tadjt3[0], 1, MPI_FLOAT, 0, MPI_COMM_WORLD);

  fval = timer_keycall.getTime().getRealTime()*1.0e-6;
  MPI_Gather(&fval, 1, MPI_FLOAT, &keycl3[0], 1, MPI_FLOAT, 0, MPI_COMM_WORLD);

  fval = timer_keycomp.getTime().getRealTime()*1.0e-6;
  MPI_Gather(&fval, 1, MPI_FLOAT, &keycm3[0], 1, MPI_FLOAT, 0, MPI_COMM_WORLD);

  fval = timer_keybods.getTime().getRealTime()*1.0e-6;
  MPI_Gather(&fval, 1, MPI_FLOAT, &keybd3[0], 1, MPI_FLOAT, 0, MPI_COMM_WORLD);

  fval = timer_waiton0.getTime().getRealTime()*1.0e-6;
  MPI_Gather(&fval, 1, MPI_FLOAT, &wait03[0], 1, MPI_FLOAT, 0, MPI_COMM_WORLD);

  fval = timer_waiton1.getTime().getRealTime()*1.0e-6;
  MPI_Gather(&fval, 1, MPI_FLOAT, &wait13[0], 1, MPI_FLOAT, 0, MPI_COMM_WORLD);

  fval = timer_waiton2.getTime().getRealTime()*1.0e-6;
  MPI_Gather(&fval, 1, MPI_FLOAT, &wait23[0], 1, MPI_FLOAT, 0, MPI_COMM_WORLD);

  fval = timer_keynewc.getTime().getRealTime()*1.0e-6;
  MPI_Gather(&fval, 1, MPI_FLOAT, &keync3[0], 1, MPI_FLOAT, 0, MPI_COMM_WORLD);

  fval = timer_keyoldc.getTime().getRealTime()*1.0e-6;
  MPI_Gather(&fval, 1, MPI_FLOAT, &keyoc3[0], 1, MPI_FLOAT, 0, MPI_COMM_WORLD);

  fval = barrier->getTime();
  MPI_Gather(&fval, 1, MPI_FLOAT, &barri3[0], 1, MPI_FLOAT, 0, MPI_COMM_WORLD);

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
  timer_keycall.reset();
  timer_keycomp.reset();
  timer_keybods.reset();
  timer_waiton1.reset();
  timer_waiton2.reset();
  timer_keynewc.reset();
  timer_keyoldc.reset();
  timer_waiton0.reset();

  numk = numkeys;
  numkeys = 0;
}


template <typename T> 
void pHOT::getQuant(vector<T>& in, vector<T>& out)
{
  sort(in.begin(), in.end());
  out = vector<T>(3);
  for (int k=0; k<3; k++)
    out[k] = in[static_cast<int>(floor(in.size()*0.01*qtile[k]))];
}

void pHOT::Timing(vector<float>    &keymake, vector<float>    &exchange, 
		  vector<float>    &convert, vector<float>    &overlap, 
		  vector<float>    &prepare, vector<float>    &update,
		  vector<float>    &scatter, vector<float>    &repartn,
		  vector<float>    &tadjust, vector<float>    &keycall,
		  vector<float>    &keycomp, vector<float>    &keybods,
		  vector<float>    &waiton0, vector<float>    &waiton1,
		  vector<float>    &waiton2, vector<float>    &keynewc,
		  vector<float>    &keyoldc, vector<float>    &treebar,
		  vector<unsigned> &numk)
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
  getQuant<float   >(keycl3, keycall);
  getQuant<float   >(keycm3, keycomp);
  getQuant<float   >(keybd3, keybods);
  getQuant<float   >(wait03, waiton0);
  getQuant<float   >(wait13, waiton1);
  getQuant<float   >(wait23, waiton2);
  getQuant<float   >(keync3, keynewc);
  getQuant<float   >(keyoc3, keyoldc);
  getQuant<float   >(barri3, treebar);
  getQuant<unsigned>(numk3,  numk);
}


double pHOT::totalKE(double& KEtot, double& KEdsp)
{
  vector<double> state(10);

  // Test
  //
  vector<double> state1(10, 0.0);

  for (key_cell::iterator it=frontier.begin(); it != frontier.end(); it++) 
    {
      for (unsigned k=0; k<10; k++) 
	state1[k] += it->second->state[k];
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
  double mass1 = root->state[0];
  unsigned count1 = root->count;

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
