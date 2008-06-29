// Uncomment to use non-blocking sends and receives in adjustTree particle
// exchange

#define NON_BLOCK

#include <values.h>

# include <cstdlib>
# include <ctime>

#include <iostream>
#include <iomanip>
#include <fstream>
#include <sstream>
#include <vector>
#include <list>
#include <algorithm>
#include <cmath>

using namespace std;

#include "global.H"
#include "pHOT.H"

double pHOT::sides[] = {2.0, 2.0, 2.0};
double pHOT::offst[] = {1.0, 1.0, 1.0};
double pHOT::jittr[] = {0.0, 0.0, 0.0};
unsigned pHOT::neg_half = 0;

template<class U, class V>
struct pair_compare
{
  bool operator()(const pair<U, V>& a, const pair<U, V>& b)
  { return a.first<=b.first; }
};

/*
  Constructor: initialize domain
*/
pHOT::pHOT(Component *C)
{
  cc = C;			// Register the calling component

  volume = sides[0]*sides[1]*sides[2];	// Total volume of oct-tree region
  root = 0;

  gen = new ACG(11+myid);
  unit = new Uniform(-1.0, 1.0, gen);
  offset = new double [3];

  kbeg = vector<key_type>(numprocs);
  kfin = vector<key_type>(numprocs);

  key_min = (key_type)1 << 48;
  key_max = (key_type)1 << 49;

  timer_keymake.Microseconds();
  timer_xchange.Microseconds();
  timer_overlap.Microseconds();
} 


pHOT::~pHOT()
{
  delete root;
  delete unit;
  delete gen;
  delete [] offset;
}


key_type pHOT::getKey(double *p)
{
  // Out of bounds?
  //
  for (unsigned k=0; k<3; k++) { 
    if (fabs((p[k]+offset[k])/sides[k])> 1.0) {
#ifdef DEBUG
      cout << "Coordinate out of pbounds in pHOT::key: ";
      for (int l=0; l<3; l++) cout << setw(18) << p[l];
      cout << endl;
#endif
      return 0;
    }
  }

  const unsigned nbits = 16;
  // const double factor = 1<<(nbits-1);
  const double factor = 1<<nbits;
  const unsigned mask = 0x1;
  vector<unsigned long> bins(3, 0);

				// Reverse the order
  for (unsigned k=0; k<3; k++)
    bins[2-k] = (unsigned)floor( ((p[k]+offset[k])/sides[k]+neg_half)*factor );
  
  key_type place = 1;
  key_type _key = 0;
  for (unsigned i=0; i<nbits; i++) {
    for (unsigned k=0; k<3; k++) {
      _key |= (bins[k] & mask)*place;
      place = place << 1;
      bins[k] = bins[k] >> 1;
    }
  }

  _key += place;		// Leading placeholder for cell masking

  return _key;
}

string pHOT::printKey(key_type p)
{
  ostringstream sout, sret;

  unsigned short cnt = 0;
  unsigned nbits = sizeof(p)*8;
  for (unsigned k=0; k<nbits; k++) {
    sout << ( (p & 0x1) ? '1' : '0' );
    if (++cnt==3) {sout << '.'; cnt = 0;}
    p = p>>1;
  }

  string s = sout.str();	// Reverse the string
  for (unsigned k=0; k<s.size(); k++) sret << s[s.size()-1-k];

  return sret.str();
}


void pHOT::makeTree()
{
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
    MPI_Barrier(MPI_COMM_WORLD);
  }
#endif

  //
  // Make the root
  //
  root = new pCell(this);

  //
  // Make new offset (only 0 node's values matter, obviously)
  //
  for (unsigned k=0; k<3; k++) offset[k] = offst[k] + (*unit)()*jittr[k];
  MPI_Bcast(&offset[0], 3, MPI_DOUBLE, 0, MPI_COMM_WORLD);

  //
  // Add the data
  //
  key_indx::iterator it;
  pCell* p = root;
  for (it=keybods.begin(); it!=keybods.end(); it++)  {

    if (it->first < key_min || it->first >= key_max) {
      cout << "Process " << myid << ": in makeTree, key=" 
	   << hex << it->first << endl << dec;
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
	MPI_Barrier(MPI_COMM_WORLD);
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

				// Exchange boundary keys
  key_type headKey=0, tailKey=0, prevKey=0, nextKey=0;
  unsigned head_num=0, tail_num=0, next_num=0, prev_num=0;

				// Do the boundaries sequentially to prevent
				// inconstencies

  for (int n=1; n<numprocs; n++) {
				// Send the next node my tail value
				// to compare with its head
    if (myid==n-1) {

      MPI_Send(&tailKey, 1, MPI_UNSIGNED_LONG_LONG, n, 1000, MPI_COMM_WORLD);
      MPI_Send(&tail_num, 1, MPI_UNSIGNED, n, 1001, MPI_COMM_WORLD);

      MPI_Recv(&nextKey, 1, MPI_UNSIGNED_LONG_LONG, n, 1000, MPI_COMM_WORLD,
	       MPI_STATUS_IGNORE);
      MPI_Recv(&next_num, 1, MPI_UNSIGNED, n, 1001, MPI_COMM_WORLD,
	       MPI_STATUS_IGNORE);

      if (tailKey && (tailKey == nextKey)) {
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
	{
				// check validity of key
	  if (bodycell.find(keybods.begin()->first) == bodycell.end()) {
	    cout << "Process " << myid << ": bad key=" 
		 << hex << keybods.begin()->first << dec
		 << " #cells=" << bodycell.size() << endl;
	  }
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
	headKey = 0;
	head_num = 0;
      }

      MPI_Send(&headKey, 1, MPI_UNSIGNED_LONG_LONG, n-1, 1000, MPI_COMM_WORLD);
      MPI_Send(&head_num, 1, MPI_UNSIGNED, n-1, 1001, MPI_COMM_WORLD);

      MPI_Recv(&prevKey, 1, MPI_UNSIGNED_LONG_LONG, n-1, 1000, MPI_COMM_WORLD,
	       MPI_STATUS_IGNORE);
      MPI_Recv(&prev_num, 1, MPI_UNSIGNED, n-1, 1001, MPI_COMM_WORLD,
	       MPI_STATUS_IGNORE);

      if (headKey && headKey == prevKey) {
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
    vollev[lev] += volume/((key_type)1 << (3*p->level));
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
    MPI_Barrier(MPI_COMM_WORLD);
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
	set<unsigned>::iterator ib = it->second->bods.begin();
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
	totvol  += volume/((key_type)1<<(3*it->second->level));

	cout << setw(4)  << myid
	     << setw(12) << hex << it->first << dec
	     << setw(8)  << it->second->level
	     << setw(18) << num
	     << setw(18) << mass
	     << setw(18) << mass/(volume/((key_type)1<<(3*it->second->level)));
	
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
    MPI_Barrier(MPI_COMM_WORLD);
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

  const unsigned fields = 14;
  vector<unsigned> prec(fields, 18);
  for (unsigned n=0; n<4; n++) prec[n] = 10;
  prec[1] = 14;

  vector<ios_base::fmtflags> fmt(fields, ios::dec);
  fmt[1] = ios::hex;

  vector<ios_base::fmtflags> typ(fields, ios::fixed);
  typ[4] = typ[5] = ios::scientific;

  if (myid==0) {

    char labels[][18] = {
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
  

  for (int n=0; n<numprocs; n++) {

    if (n == myid) {

      ofstream out(filename.c_str(), ios::app);

      for (c=frontier.begin(); c!=frontier.end(); c++) {
	double mass=0, temp=0, pos[]={0,0,0}, vel[]={0,0,0};
	p = c->second;
	set<unsigned>::iterator ib = p->bods.begin();
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
	double vol = volume/((key_type)1 << (3*p->level)); 

	out << " ";
	out << setw(prec[n]) << setiosflags(fmt[n]|typ[n]) << myid;
	n++;
	out << setw(prec[n]) << setiosflags(fmt[n]|typ[n]) << c->first;
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

    MPI_Barrier(MPI_COMM_WORLD);
  }
}


void pHOT::sendCell(key_type key, int to, unsigned num)
{
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
  set<unsigned>::iterator ib = p->bods.begin();
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
    p->parent->children.erase(p->mykey & 0x7);
	
    // DEBUG
    if (frontier.find(p->mykey)==frontier.end()) {
        cout << "Process " << myid << ": in pHOT:sendCell: "
	     << " key not on frontier as expected" << endl;
    }
    // END

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
}


void pHOT::recvCell(int from, unsigned num)
{
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
	     << ", key=" << hex << part.key << dec << "]";
      }
    cc->particles[part.indx] = part;
    if (part.key == 0) continue;
    if (part.key < key_min || part.key >= key_max) {
      cout << "Process " << myid << ": in recvCell, key=" 
	   << hex << part.key << endl << dec;
    }
    if (part.indx==0) cout << "pHOT::recvCell bad particle indx=0!" << endl;
    key_pair tpair(part.key, part.indx);
    keybods.insert(tpair);
    p = p->Add(tpair);
  }

  cc->nbodies = cc->particles.size();
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
      
      dens = state[0] * ((key_type)1 << (3*clv))/(volume*cnt);
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

				// X-Y slice
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
  vol1 = volume/((key_type)1 << (3*MaxLev));
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
  vol1 = volume/((key_type)1 << (3*MinLev));
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

  return volume/((key_type)1 << (3*mlev));
}

void pHOT::Repartition()
{
  // debug_ts("Entering Repartition");

  map<unsigned long, Particle>::iterator it;
  
  volume = sides[0]*sides[1]*sides[2]; // Total volume of oct-tree region


				// No need to repartition 
				// if there are no bodies
  if (cc->nbodies_tot==0) {
    if (myid==0) 
      cout << "pHOT::Repartition with ZERO bodies, continuing" << endl;
    return;
  }

#ifdef DEBUG
  vector<unsigned long> erased;
#endif

  //
  // Recompute keys and compute new partition
  //
  oob.clear();
  vector<key_type> keys;
  for (it=cc->Particles().begin(); it!=cc->Particles().end(); it++) {
    it->second.key = getKey(&(it->second.pos[0]));
    if (it->second.key == 0) {
      oob.insert(it->first);
    } else {
      keys.push_back(it->second.key);
    }
  }
#ifdef DEBUG
  cout << "Process " << myid 
       << ": part #=" << cc->Particles().size()
       << "  key size=" << keys.size()
       << "  oob size=" << oob.size() << endl;
#endif
  // debug_ts("Before spreadOOB");
  spreadOOB();
  // debug_ts("Before partitionKeys");
  partitionKeys(keys, kbeg, kfin);
  // debug_ts("After partitionKeys");

  //
  // Nodes compute send list
  //
  loclist = kbeg;
  loclist.push_back(kfin[numprocs-1]);	// End point for binary search

  vector<unsigned> sndlist1(numprocs*numprocs, 0);
  vector<unsigned> sendlist(numprocs*numprocs);
  vector< vector<unsigned> > bodylist(numprocs);
  unsigned t;
  for (it=cc->Particles().begin(); it!=cc->Particles().end(); it++) {
				// Skip an OOB particle
    if (it->second.key == 0) continue;
				// Look for key in this node's list
    t = find_proc(loclist, it->second.key);
    if (t == numprocs) {
      cerr << "Process " << myid << ": loclist found last entry, "
	   << " key=" << hex << it->second.key 
	   << ", end pt=" << loclist.back() << dec
	   << ", index=" << t << endl;
    }
    if (t == myid) continue;
    bodylist[t].push_back(it->first);
    sndlist1[numprocs*myid + t]++;
  }

  int ntot = numprocs*numprocs;
  MPI_Allreduce(&sndlist1[0], &sendlist[0], ntot, MPI_UNSIGNED, MPI_SUM, MPI_COMM_WORLD);
  
  unsigned Tcnt=0, Fcnt=0;
  for (int i=0; i<numprocs; i++) {
    Tcnt += sendlist[numprocs*myid + i];
    Fcnt += sendlist[numprocs*i + myid];
  }

				// DEBUG OUTPUT
  if (false) {			// If true, write send and receive list for each node
    for (int n=0; n<numprocs; n++) {
      if (myid==n) {
	cout<<"--------------------------------------------------------"<<endl
	    <<"---- Repartition"<<endl
	    <<"--------------------------------------------------------"<<endl
	    << "Process " << myid << ": Tcnt=" << Tcnt 
	    << " Fcnt=" << Fcnt << endl;
	for (int m=0; m<numprocs; m++)
	  cout << setw(5) << m 
	       << setw(8) << sendlist[numprocs*n+m]
	       << setw(8) << sendlist[numprocs*m+n]
	       << endl;
	cout << "--------------------------------------------------------" << endl;
      }
      MPI_Barrier(MPI_COMM_WORLD);
    }
  }
				// END DEBUG OUTPUT

  static unsigned debug_ctr=0;
  unsigned ps=0, pr=0;
  Partstruct *psend=0, *precv=0;
  vector<MPI_Request> rql;
  MPI_Request req;
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
	    pf.Particle_to_part(psend[ps+i], cc->Particles()[bodylist[toID][i]]);
	    cc->Particles().erase(bodylist[toID][i]);
#ifdef DEBUG
	    erased.push_back(bodylist[toID][i]);
#endif
	  }
	  for (unsigned i=0; i<To; i++) {
	    if (psend[ps+i].indx == 0) {
	      cerr << "pHOT::Repartition[" << debug_ctr
		   << "]: SEND from=" << frID << " to=" << toID << " #=" << To
		   << " ps=" << ps+i << " mass=" << scientific << psend[ps+i].mass << endl;
	    }
	  }
#ifdef NON_BLOCK
	  rql.push_back(req);
#endif
	  if ( (ierr=
#ifdef NON_BLOCK
		MPI_Isend(&psend[ps], To, ParticleFerry::Particletype, toID, 49, MPI_COMM_WORLD, &rql.back())
#else
		MPI_Send(&psend[ps], To, ParticleFerry::Particletype, toID, 49, MPI_COMM_WORLD)
#endif
		) != MPI_SUCCESS  ) {
	    cout << "Process " << myid << ": error in Reparition sending "
		 << To << " particles to #" << toID << " ierr=" << ierr << endl;
	  }
	  ps += To;
	}
      }
				// 
				// Current process receives particles (non blocking)
	if (myid==toID) {		// 
	unsigned From = sendlist[numprocs*frID+toID];
	if (From) {
#ifdef NON_BLOCK
	  rql.push_back(req);
#endif
	  if ( (ierr=
#ifdef NON_BLOCK
		MPI_Irecv(&precv[pr], From, ParticleFerry::Particletype, frID, 49, MPI_COMM_WORLD, &rql.back())
#else
		MPI_Recv(&precv[pr], From, ParticleFerry::Particletype, frID, 49, MPI_COMM_WORLD, MPI_STATUS_IGNORE)
#endif
		) != MPI_SUCCESS ) {
	    cout << "Process " << myid << ": error in Reparition receiving "
		 << From << " particles from #" << frID << " ierr=" << ierr << endl;
	  }
	  pr += From;
	}
      }
      
    } // Receipt loop
    
    MPI_Barrier(MPI_COMM_WORLD); // Ok, everybody move on to next sending node
  }

  //
  // Wait for completion of sends and receives
  //

#ifdef NON_BLOCK
  if ( (ierr=MPI_Waitall(rql.size(), &rql[0], MPI_STATUSES_IGNORE)) != MPI_SUCCESS ) 
    {
      cout << "Process " << myid << ": error in Reparition Waitall"
	   << ", ierr=" << ierr << endl;
    }
#endif

  //
  // Buffer size sanity check
  //
  if (true) {			// Enabled if true

    if (Tcnt-ps) {
      cout << "Process " << myid 
	   << ": Repartition [Tcnt] found <" << ps << "> but expected <" 
	   << Tcnt << ">" << endl;
    }
    
    if (Fcnt-pr) {
      cout << "Process " << myid 
	   << ": Repartition [Fcnt] found <" << pr << "> but expected <" 
	   << Fcnt << ">" << endl;
    }

  }

  //
  // DEBUG
  //
  pr = 0;
  for (int id=0; id<numprocs; id++) {
    unsigned From = sendlist[numprocs*id+myid];
    if (From) {
      for (int i=0; i<From; i++) {
	if (precv[pr+i].indx == 0) {
	  cerr << "pHOT::Repartition[" << debug_ctr 
	       << "]: RECV to=" << myid << " from=" << id << " #=" << From
	       << " pr=" << pr+i << " mass=" << scientific << precv[pr+i].mass 
	       << endl;
	}
      }
      pr += From;
    }
  }

  debug_ctr++;

  // END DEBUG

  if (Fcnt) {
    Particle part;
    for (unsigned i=0; i<Fcnt; i++) {
      pf.part_to_Particle(precv[i], part);
      if (part.mass<=0.0 || isnan(part.mass)) {
	cout << "[adjustTree crazy mass indx=" << part.indx 
	     << ", key=" << hex << part.key << dec << "]";
      }
      cc->Particles()[part.indx] = part;
    }
    cc->nbodies = cc->particles.size();
  }
      
  //
  // Clean up temporary body storage
  //
  if (Tcnt) delete [] psend;
  if (Fcnt) delete [] precv;

  //
  // Remake key body index
  //
  keybods.clear();
  unsigned oob1_cnt=0, oob_cnt=0;
  map<unsigned long, Particle>::iterator n;
  for (n=cc->Particles().begin(); n!=cc->Particles().end(); n++) {
    if (n->second.key==0) {
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
    map<unsigned long, Particle>::iterator ip;
    for (ip=cc->particles.begin(); ip!=cc->particles.end(); ip++) {
      if (ip->second.indx==0) {
	cout << "pHOT::Repartition BAD particle in proc=" << myid
	     << " key=" << ip->second.key << endl;
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
}


void pHOT::makeCellLevelList()
{
				// Make new lists
  clevlst.clear();
  clevels = vector< set<pCell*> >(multistep+1);

  unsigned ng=0, nt=0;
  for (key_cell::iterator it=frontier.begin(); it != frontier.end(); it++) {
    nt++;
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
    cout << "Process " << myid << ": made level list with " << ng
	 << " good cells out of " << nt << " expected" << endl;

#ifdef DEBUG
  vector<unsigned> pcnt(multistep+1, 0);
  for (map<pCell*, unsigned>::iterator
	 pit=clevlst.begin(); pit!=clevlst.end(); pit++)
    pcnt[pit->second]++;

  cout << left << setw(60) << setfill('-') << "-" << endl << setfill(' ')
       << setw(10) << "M" << setw(10) << "number" 
       << setw(10) << "counts" << endl;
  for (unsigned M=0; M<=multistep; M++)
    cout << setw(10) << M << setw(10) << clevels[M].size() 
	 << setw(10) << pcnt[M] << endl;
  cout << left << setw(60) << setfill('-') << "-" << endl << setfill(' ');
  
  checkParticles();
#endif // END DEBUG
}

void pHOT::adjustCellLevelList(unsigned mlevel)
{
  if (multistep==0) return;	// No need to bother if multistepping is off
				// Otherwise . . . 
  unsigned ng=0, nt=0, ns=0, m, cnt;
  for (unsigned M=mlevel; M<=multistep; M++) {
    nt += clevels[M].size();
    cnt = 0;
    if (clevels[M].size()>0) {
      set<pCell*>::iterator it = clevels[M].begin(), nit;
      while (it != clevels[M].end()) {
	cnt++;			// Count the cells

				// Debugging: This shouldn't happen
	if ((*it)->bods.size()) ng++;
	else {
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

				// Debugging . . .
  if (nt!=ng)
    cout << "Process " << myid << ": adjusted level list with " << ng
	 << " good cells out of " << nt << " expected, " << ns
	 << " cells moved" << endl;
#ifdef DEBUG
  checkParticles();
#endif
}

void pHOT::adjustTree(unsigned mlevel)
{

#ifdef DEBUG
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
  vector< vector<unsigned long> > exchange(numprocs);

  // 
  // Make body list from frontier cells for this level
  //
  for (unsigned M=mlevel; M<=multistep; M++) {
    for (set<pCell*>::iterator it=clevels[M].begin(); 
	 it!=clevels[M].end(); it++) {
      for (set<unsigned>::iterator ib=(*it)->bods.begin(); 
	   ib!=(*it)->bods.end(); ib++) {
	oldp.push_back(*ib);
      }
    }
  }
  
  //
  // Update body by body using the list without regard to level
  //
  unsigned newproc;
  for (ip=oldp.begin(); ip!=oldp.end(); ip++) {

    Particle *p = cc->Part(*ip);
    if (p==0) {			// Sanity check
      cout << "Process " << myid 
	   << ": pHOT::adjustTree: ERROR crazy particle index!" << endl;
    }

    //
    // Get and recompute keys
    //
    oldkey = p->key;
    newkey = getKey(&(p->pos[0]));

    //
    // Get this particle's cell
    //

    //
    // Look for keys . . .
    //
    if (bodycell.find(oldkey) == bodycell.end()) {
      cout << "Process " << myid 
	   << ": pHOT::adjustTree: ERROR could not find cell for particle"
	   << " key=" << hex << oldkey << ", index=" << dec << p->indx
	   << " pnumber=" << cc->Number() << " bodycell=" << bodycell.size() 
	   << endl;
    }
    //
    // Look for cell in frontier . . .
    //
    if (frontier.find(bodycell.find(oldkey)->second) == frontier.end()) {
      cout << "Process " << myid 
	   << ": pHOT::adjustTree: ERROR could not find expected cell"
	   << " on frontier, count=" << adjcnt << hex
	   << " oldbody=" << oldkey 
	   << " newbody=" << newkey 
	   << " cell="    << bodycell.find(oldkey)->second << dec
	   << " index="   << p->indx 
	   << endl;
      continue;
    }
      
    //
    // Find this particle's previous cell assignment
    //
    c = frontier[bodycell.find(oldkey)->second];
    
    //
    // Is the key the same?
    // 
    if (newkey != oldkey) {
				// Key pairs
      key_pair newpair(newkey, p->indx);
      key_pair oldpair(oldkey, p->indx);

				// Put the particle in a new cell?
				// 
      if ( !(c->isMine(newkey)) ) {

	if (c->Remove(oldpair, &change)) {
				// Remove the old pair from the current cell
				// (only transactions added are sample cells)
				// queue for removal from level lists
	  change.push_back(cell_indx(c, REMOVE));
				// queue for deletion
	  change.push_back(cell_indx(c, KILL));
	}

				// Same processor and in bounds?
	if (newkey) {
	  newproc = find_proc(loclist, newkey);
	  if (newproc != myid) {
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

      } else {			// Same cell: update body cell index 
				// for the new key
	key_key::iterator ij = bodycell.find(oldkey);
	// DEBUG
	if (ij == bodycell.end()) {
	  cout << "Process " << myid 
	       << ": pHOT::adjustTree: ERROR could not find cell for"
	       << " key=" << hex << oldkey << ", index=" << dec << p->indx
	       << endl;
	}
	// END DEBUG
	if (ij != bodycell.end()) bodycell.erase(ij);
	bodycell.insert(key_item(newkey, c->mykey));
				// Update key list
	c->UpdateKeys(oldpair, newpair);
				// Update key body index for the new key
	key_indx::iterator ik = keybods.find(oldpair);
	if (ik != keybods.end()) keybods.erase(ik);
	// DEBUG
	else {
	  cout << "Process " << myid 
	       << ": pHOT::adjustTree: ERROR mlevel=" << mlevel 
	       << ": could not find keybods entry" << endl;
	}
	// END DEBUG
	
	p->key = newkey;	// Assign the new key to the particle
	keybods.insert(newpair);
      }
    }
    
  }
  
  timer_keymake.stop();
  timer_xchange.start();


  //
  // Exchange particles
  //

  vector<unsigned long> templist(numprocs*numprocs, 0);
  vector<unsigned long> sendlist(numprocs*numprocs, 0);

  for (unsigned k=0; k<numprocs; k++)
    templist[numprocs*myid + k] = exchange[k].size();

  MPI_Allreduce(&templist[0], &sendlist[0], numprocs*numprocs,
		MPI_UNSIGNED_LONG, MPI_SUM, MPI_COMM_WORLD);


  unsigned Tcnt=0, Fcnt=0;
  for (unsigned k=0; k<numprocs; k++) {
    Tcnt += sendlist[numprocs*myid + k   ];
    Fcnt += sendlist[numprocs*k    + myid];
  }


  static unsigned debug_ctr=0;
  unsigned ps=0, pr=0;
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
				// (non-blocking)
      if (myid==frID) {		// 
	unsigned To = sendlist[numprocs*frID+toID];
	if (To) {
	  for (unsigned i=0; i<To; i++) {
	    pf.Particle_to_part(psend[ps+i], 
				cc->Particles()[exchange[toID][i]]);
	    cc->Particles().erase(exchange[toID][i]);
	  }
	  for (unsigned i=0; i<To; i++) {
	    if (psend[ps+i].indx == 0) {
	      cerr << "pHOT::adjustTree[" << debug_ctr
		   << "]: SEND from=" << frID << " to=" << toID << " #=" << To
		   << " ps=" << ps+i 
		   << " mass=" << scientific << psend[ps+i].mass << endl;
	    }
	  }
	  rql.push_back(r);
	  if ( (ierr=MPI_Isend(&psend[ps], To, ParticleFerry::Particletype, toID, 49, MPI_COMM_WORLD, &rql.back())) != MPI_SUCCESS ) {
	    cout << "Process " << myid << ": error in adjustTree sending"
		 << To << " particles to #" << toID 
		 << " ierr=" << ierr << endl;
	  }
	  ps += To;
	}
      }	// Send block

				// Current process receives particles 
				// (non-blocking)
      if (myid==toID) {		// 
	unsigned From = sendlist[numprocs*frID+toID];
	if (From) {
	  rql.push_back(r);
	  if ( (ierr=MPI_Irecv(&precv[pr], From, ParticleFerry::Particletype, frID, 49, MPI_COMM_WORLD, &rql.back())) != MPI_SUCCESS ) {
	    cout << "Process " << myid << ": error in adjustTree receiving "
		 << From << " particles from #" << frID 
		 << " ierr=" << ierr << endl;
	  }
	  pr += From;
	}
      }

    } // Receipt block

  }

  //
  // Wait for completion of sends and receives
  //

  if ( (ierr=MPI_Waitall(rql.size(), &rql[0], MPI_STATUSES_IGNORE)) != MPI_SUCCESS ) 
    {
      cout << "Process " << myid << ": error in adjustTree Waitall"
	   << ", ierr=" << ierr << endl;
    }

  //
  // Buffer size sanity check
  //
  if (true) {			// Enabled if true

    if (Tcnt-ps) {
      cout << "Process " << myid 
	   << ": adjustTree [Tcnt] found <" << ps << "> but expected <" 
	   << Tcnt << ">" << endl;
    }
    
    if (Fcnt-pr) {
      cout << "Process " << myid 
	   << ": adjustTree [Fcnt] found <" << pr << "> but expected <" 
	   << Fcnt << ">" << endl;
    }

  }

  //
  // DEBUG
  //
  if (false) {			// Enabled if true

    pr = 0;
    for (int id=0; id<numprocs; id++) {
      unsigned From = sendlist[numprocs*id+myid];
      if (From) {
	for (int i=0; i<From; i++) {
	  if (precv[pr+i].indx == 0) {
	    cerr << "pHOT::adjustTree[" << debug_ctr 
		 << "]: RECV to=" << myid << " from=" << id << " #=" << From
		 << " pr=" << pr+i << " mass=" << scientific 
		 << precv[pr+i].mass << endl;
	  }
	}
	pr += From;
      }
    }
  }

  debug_ctr++;

  // END DEBUG

  Particle part;
  for (unsigned i=0; i<Fcnt; i++) {
    pf.part_to_Particle(precv[i], part);
      if (part.mass<=0.0 || isnan(part.mass)) {
	cout << "[spreadOOB crazy mass indx=" << part.indx 
	     << ", key=" << hex << part.key << dec << "]";
      }

    cc->Particles()[part.indx] = part;
    
    if (part.key) {
      key_pair newpair(part.key, part.indx);
      keybods.insert(newpair);
      root->Add(newpair, &change);
    }
  }
  cc->nbodies = cc->particles.size();


  timer_xchange.stop();
  timer_overlap.start();

  //
  // Cell overlap?
  //

  key_type headKey=0, tailKey=0;
  unsigned head_num=0, tail_num=0;

  for (int n=0; n<numprocs-1; n++) {

    if (n==myid) {
      if (keybods.size()) {
	key_indx::reverse_iterator it = keybods.rbegin();
	tailKey = bodycell.find(it->first)->second;
	tail_num = frontier[tailKey]->bods.size();
      } else {
	tailKey = 0;
	tail_num = 0;
      }

      // Send my tail cell info to next node
      MPI_Send(&tailKey, 1, MPI_UNSIGNED_LONG_LONG, n+1, 131, MPI_COMM_WORLD);
      MPI_Send(&tail_num, 1, MPI_UNSIGNED, n+1, 132, MPI_COMM_WORLD);
      
      // Get next node's head cell info
      MPI_Recv(&headKey, 1, MPI_UNSIGNED_LONG_LONG, n+1, 133, 
	       MPI_COMM_WORLD, MPI_STATUS_IGNORE);
      MPI_Recv(&head_num, 1, MPI_UNSIGNED, n+1, 134,
	       MPI_COMM_WORLD, MPI_STATUS_IGNORE);

      if (tailKey==headKey && tailKey!=0) {

	c = frontier[tailKey];	// Cell in question: my last cell

	if (tail_num>head_num) { 
	  //
	  // Receive particles
	  //
	  vector<Partstruct> Precv(head_num);
	  MPI_Recv(&Precv[0], head_num, ParticleFerry::Particletype, 
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
	  set<unsigned>::iterator ib;
	  for (ib=c->bods.begin(); ib!=c->bods.end(); ib++) {
	    pf.Particle_to_part(Psend[k++], cc->Particles()[*ib]);
	    cc->Particles().erase(*ib);
	  }

	  c->RemoveAll();

	  // queue for removal from level lists
	  change.push_back(cell_indx(c, REMOVE));

	  // queue for deletion
	  change.push_back(cell_indx(c, KILL));

	  MPI_Send(&Psend[0], tail_num, ParticleFerry::Particletype, 
		   n+1, 135, MPI_COMM_WORLD);
	}
      }
    }
    
    if (n+1==myid) {
      if (keybods.size()) {
	key_indx::iterator it = keybods.begin();
	headKey = bodycell.find(it->first)->second;
	head_num = frontier[headKey]->bods.size();
      } else {
	headKey = 0;
	head_num = 0;
      }

      // Get previous nodes tail cell info
      MPI_Recv(&tailKey, 1, MPI_UNSIGNED_LONG_LONG, n, 131, 
	       MPI_COMM_WORLD, MPI_STATUS_IGNORE);
      MPI_Recv(&tail_num, 1, MPI_UNSIGNED, n, 132,
	       MPI_COMM_WORLD, MPI_STATUS_IGNORE);

      // Send head node into to previous node
      MPI_Send(&headKey, 1, MPI_UNSIGNED_LONG_LONG, n, 133, MPI_COMM_WORLD);
      MPI_Send(&head_num, 1, MPI_UNSIGNED, n, 134, MPI_COMM_WORLD);

      if (tailKey==headKey && tailKey!=0) {

	c = frontier[headKey];	// Cell in question

	if (tail_num>head_num) { 
	  //
	  // Send particles
	  //
	  unsigned k=0;
	  vector<Partstruct> Psend(head_num);
	  set<unsigned>::iterator ib;
	  for (ib=c->bods.begin(); ib!=c->bods.end(); ib++) {
	    pf.Particle_to_part(Psend[k++], cc->Particles()[*ib]);
	    cc->Particles().erase(*ib);
	  }

	  c->RemoveAll();

	  // queue for removal from level lists
	  change.push_back(cell_indx(c, REMOVE));

	  // queue for deletion
	  change.push_back(cell_indx(c, KILL));
	  
	  MPI_Send(&Psend[0], head_num, ParticleFerry::Particletype, 
		   n, 136, MPI_COMM_WORLD);
	  
	} else {		
	  //
	  // Receive particles
	  //
	  vector<Partstruct> Precv(tail_num);
	  MPI_Recv(&Precv[0], tail_num, ParticleFerry::Particletype, 
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
    
    MPI_Barrier(MPI_COMM_WORLD);
  }


  timer_overlap.stop();

  // #define DEBUG
  
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

    case RECOMP:		// Find the sample cell
      c->findSampleCell();
      break;

    default:
      cout << "Process " << myid << ": unknown action in pHOT::adjustTree()!!"
	   << endl;
    }
  }

  // #undef DEBUG

  change.clear();		// Reset the change list for next time
  
#ifdef DEBUG
  if (!checkBodycell()) {
    cout << "Process " << myid << ": "
	 << "adjustTree: ERROR mlevel=" << mlevel 
	 << " completed: body cell check FAILED!" << endl;
  }    
  if (!checkParticles()) {
    cout << "Process " << myid << ": "
	 << "adjustTree: ERROR mlevel=" << mlevel 
	 << " completed: initial particle check FAILED!" << endl;
  }    
  if (!checkFrontier()) {
    cout << "Process " << myid << ": "
	 << "adjustTree: ERROR mlevel=" << mlevel
	 << ": frontier check FAILED!" << endl;
  }

  checkDupes();
  checkIndices();

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
  
}
  

bool pHOT::checkParticles(bool pc)
{
  unsigned cnt, badb=0, badc=0;

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
	    for (set<unsigned>::iterator ib=(*it)->bods.begin();
		 ib!=(*it)->bods.end(); ib++) {
	      if (!cc->Part(*ib)) {
		cout << "pHOT::checkParticles:: M=" << M << ", bad body at "
		     << cnt << "/" << clevels[M].size() 
		     << " cell=" << hex << (*it) << dec << endl;
		badb++;
	      }
	    }
	  }
	} else {
	  cout << "pHOT::checkParticles:: M=" << M << ", zero bods at "
	       << cnt << "/" << clevels[M].size() 
	       << " cell=" << hex << (*it) << dec << endl;
	  badc++;
	}
      }
    }
  }
  
  // OUTPUT CELL LIST
  ostringstream cmd1, cmd2;

  cmd1 << "mv chkcell." << myid << " chkcell." << myid << ".bak";
  system(cmd1.str().c_str());

  cmd2 << "mv chklist." << myid << " chklist." << myid << ".bak";
  system(cmd2.str().c_str());

  ostringstream newf1, newf2;
  newf1 << "chkcell." << myid;
  newf2 << "chklist." << myid;
  ofstream out1(newf1.str().c_str());
  ofstream out2(newf2.str().c_str());

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
    cout << "Process " << myid << ": pHOT::checkParticles, bad cell=" 
	 << badc ;
    if (pc) cout << " bad bod=" << badb;
    cout << endl;
    return false;
  }
  else return true;
}


bool pHOT::checkFrontier()
{
  unsigned bad=0;
  bool good=true;

  for (unsigned M=0; M<=multistep; M++) {
    if (clevels[M].size()) {
      for (set<pCell*>::iterator 
	     it=clevels[M].begin(); it!=clevels[M].end(); it++) {
	if (frontier.find((*it)->mykey) == frontier.end()) {
	  cout << "pHOT::checkFrontier error on M=" << M
	       << ", cell=" << hex << *it << dec << endl;
	  bad++;
	  good = false;
	}
      }
    }
  }
  
  if (bad) {
    cout << "Process " << myid << ": pHOT::checkFrontier, bad cell=" 
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

  /*
  if (myid==0)
    cout << "Process " << myid << ": # indices=" << indices.size() << endl;
  */

  // All processes send the index list to the master node
  //
  for (int n=1; n<numprocs; n++) {
    /*
    if (myid==n)
      cout << "Process " << myid << ": # indices=" << indices.size() << endl;
    */

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
    MPI_Barrier(MPI_COMM_WORLD);
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
  map<unsigned long, Particle>::iterator n = cc->Particles().begin();
  for (; n!=cc->Particles().end(); n++) {
    key_pair tpair(n->second.key, n->second.indx);
    key_indx::iterator it = keybods.find(tpair);

    if (it==keybods.end()) {
      ok = false;
      cnt++;
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
      for (set<unsigned>::iterator ib=(*it)->bods.begin(); ib!=(*it)->bods.end(); ib++) {
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
  map<unsigned long, Particle>::iterator n;
  for (n=cc->Particles().begin(); n!=cc->Particles().end(); n++) {
    // Ignore OOB particle
    if (n->second.key==0) continue;
    // Look for the bodycell
    key_key::iterator it = bodycell.find(n->second.key);
    if (it==bodycell.end()) {
      ok = false;
      cnt++;
#ifdef DEBUG
      cout << "Process " << myid << ": checkBodycell: " 
	   << cnt << " unmatched particle: key=" << hex
	   << n->second.key << dec << " index=" 
	   << n->second.indx << endl;
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
  map<unsigned long, Particle>::iterator n = cc->Particles().begin();
  for (; n!=cc->Particles().end(); n++) {
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
				// 3% tolerance
  const unsigned long tol = 33;

  vector<long> list1(numprocs, 0), list0(numprocs), delta(numprocs);

  list1[myid] = oob.size();
  MPI_Allreduce(&list1[0], &list0[0], numprocs, 
		MPI_LONG, MPI_SUM, MPI_COMM_WORLD);

  // debug_ts("In spreadOOB, after reduce");


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
    debug_ts(sout.str().c_str());
  }

				// Don't bother if changes are small
  if (maxdif < cc->nbodies_tot/tol) return;

				// Nothing to send or receive
  if (nsend.size()==0 || nrecv.size()==0) return;

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
  // debug_ts("Before checkOOB");
  checkOOB(sendlist);
  // debug_ts("After checkOOB");
  // END DEBUG

  unsigned Tcnt=0, Fcnt=0;
  for (int i=0; i<numprocs; i++) {
    Tcnt += sendlist[numprocs*myid + i];
    Fcnt += sendlist[numprocs*i + myid];
  }

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

  unsigned ps=0, pr=0;
  set<indx_type>::iterator ioob;
  Partstruct *psend=0, *precv=0;
  vector<MPI_Request> rql;
  MPI_Request r;
  int ierr;

  if (Tcnt) psend = new Partstruct [Tcnt];
  if (Fcnt) precv = new Partstruct [Fcnt];

  // debug_ts("In spreadOOB, before exchange");

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
	  if ( (ierr=MPI_Isend(&psend[ps], To, ParticleFerry::Particletype, 
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
	  if ( (ierr=MPI_Irecv(&precv[pr], From, ParticleFerry::Particletype, 
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

  // debug_ts("In spreadOOB, after exchange, before completion");

  //
  // Wait for completion of sends and receives
  //

  if ( (ierr=MPI_Waitall(rql.size(), &rql[0], MPI_STATUSES_IGNORE)) != MPI_SUCCESS ) 
    {
      cout << "Process " << myid << ": error in spreadOOB Waitall"
	   << ", ierr=" << ierr << endl;
    }

  // debug_ts("In spreadOOB, after completion");

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

  if (Tcnt) delete [] psend;
  if (Fcnt) delete [] precv;

}


void pHOT::partitionKeys(vector<key_type>& keys, 
			 vector<key_type>& kbeg, vector<key_type>& kfin)
{
				// Sort the keys
  sort(keys.begin(), keys.end());

				// Sample keys at desired rate
  unsigned srate = 1000;	// e.g. every srate^th key is sampled

				// Number of samples
  unsigned nsamp = keys.size()/srate;
				
  if (nsamp < 10000/numprocs)	// Too few samples
    nsamp = 10000/numprocs; 
				
  if (nsamp > keys.size())	// Too many for particle count
    nsamp = keys.size();

  if (nsamp)			// Consistent sampling rate
    srate = keys.size()/nsamp;

  // Require: srate*(2*nsamp-1)/2 < keys.size()
  //
  if (keys.size()) {
    unsigned decr=0;
    while (srate*(2*nsamp-1)/2 >= keys.size()) {
      nsamp--;
      // DEBUG
      decr++;
      // END DEBUG
    }
    // DEBUG
    if (decr) 
      cout << "partitionKeys: process " << myid 
	   << ": decreased nsamp by " << decr  << " to " << nsamp 
	   << endl;
    // END DEBUG
  }
  
  vector<key_type> keylist1, keylist;
  if (keys.size()) {
    for (unsigned i=0; i<nsamp; i++)
      keylist1.push_back(keys[srate*(2*i+1)/2]);
  }

  //
  // DEBUG (set to "true" to enable key list diagnostic)
  //
  if (false) {
    const unsigned cols = 3;	// # of columns in output
    const unsigned cwid = 35;	// column width

    if (myid==0) {
      cout << "--------------------------" << endl;
      cout << "----- Partition keys -----" << endl;
      cout << "--------------------------" << endl;

      for (unsigned q=0; q<cols; q++) {
	ostringstream sout;
	sout << left << setw(4) << "Id" << setw(8) 
	     << "Size" << setw(8) << "Tot" << setw(8) << "Order";
	cout << left << setw(cwid) << sout.str();
      }
      cout << endl;
	
      for (unsigned q=0; q<cols; q++) {
	ostringstream sout;
	sout << setfill('-') << setw(cwid-5) << '-';
	cout << left << setw(cwid) << sout.str();
      }
      cout << endl;
    }

    for (unsigned n=0; n<numprocs; n++) {
      if (myid==n) {
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

	cout << left << setw(cwid) << sout.str() << flush;
	if (n % cols == cols-1 || n == numprocs-1) cout << endl;
      }
      MPI_Barrier(MPI_COMM_WORLD);
    }
  }
  //
  // END DEBUG
  //

				// Tree aggregation (merge sort) of the
				// entire key list

				// <keylist1> is (and must be) sorted to start
  MPI_Barrier(MPI_COMM_WORLD);
  parallelMerge(keylist1, keylist);
  MPI_Barrier(MPI_COMM_WORLD);

  if (myid==0) {

    for (unsigned i=0; i<numprocs-1; i++) {
      if (keylist.size())
	kfin[i] = keylist[static_cast<unsigned>(keylist.size()*(i+1)/numprocs)];
      else
	kfin[i] = key_min;
    }

    kfin[numprocs-1] = key_max;

    kbeg[0] = key_min;
    for (unsigned i=1; i<numprocs; i++)
      kbeg[i] = kfin[i-1];
      
    if (false) {	 // If true, print key ranges for each process
      cout << "--------------------------------------------------" << endl
	   << "---- partitionKeys: keys in list=" << keylist.size() << endl
	   << "--------------------------------------------------" << endl;
      for (int i=0; i<numprocs; i++)
	cout << setw(5) << i << hex << setw(15) << kbeg[i] 
	     << setw(15) << kfin[i] << dec << endl;
      cout << "--------------------------------------------------" << endl;
    }
  }

  if (false) {			// If true, print key totals
    unsigned oobn = oobNumber();
    unsigned tkey1 = keys.size(), tkey0 = 0;
    MPI_Reduce(&tkey1, &tkey0, 1, MPI_UNSIGNED, MPI_SUM, 0, MPI_COMM_WORLD);
    if (myid==0) 
      cout << endl
	   << setfill('-') << setw(60) << '-' << setfill(' ') << endl
	   << "---- partitionKeys" << endl
	   << setfill('-') << setw(60) << '-' << setfill(' ') << endl
	   << "----  list size=" << setw(10) << keylist.size() << endl
	   << "---- total keys=" << setw(10) << tkey0 << endl
	   << "----  total oob=" << setw(10) << oobn << endl
	   << "----      TOTAL=" << setw(10) << tkey0 + oobn << endl
	   << setfill('-') << setw(60) << '-' << setfill(' ') << endl;
  }

  MPI_Bcast(&kbeg[0], numprocs, MPI_UNSIGNED_LONG_LONG, 0, MPI_COMM_WORLD);
  MPI_Bcast(&kfin[0], numprocs, MPI_UNSIGNED_LONG_LONG, 0, MPI_COMM_WORLD);
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
void pHOT::sortCombine(vector<key_type>& one, vector<key_type>& two,
		       vector<key_type>& comb)
{
  int i=0, j=0;
  int n = one.size()-1;
  int m = two.size()-1;
  
  comb = vector<key_type>(one.size()+two.size());

  for(int k=0; k<n+m+2; k++) {
    if (i > n)
      comb[k] = two[j++];
    else if(j > m)
      comb[k] = one[i++];
    else {
      if(one[i] < two[j])
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
void pHOT::parallelMerge(vector<key_type>& initl, vector<key_type>& final)
{
  MPI_Status status;
  vector<key_type> work;
  unsigned n;

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
      MPI_Send(&initl[0], n, MPI_UNSIGNED_LONG_LONG, myid-M2, 12, 
	       MPI_COMM_WORLD);
    }
    return;
  }

  vector<key_type> data = initl;

  //
  // Retrieve the excess particles
  //
  if (myid + M2 < numprocs) {
    MPI_Recv(&n, 1, MPI_UNSIGNED, myid+M2, 11, MPI_COMM_WORLD, &status);
    if (n) {
      vector<key_type> recv(n);
      MPI_Recv(&recv[0], n, MPI_UNSIGNED_LONG_LONG, status.MPI_SOURCE, 12, 
	       MPI_COMM_WORLD, MPI_STATUS_IGNORE);
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
	MPI_Send(&data[0], n, MPI_UNSIGNED_LONG_LONG, myid-M2, 12, 
		 MPI_COMM_WORLD);
      }
      return;
    } else {
      MPI_Recv(&n, 1, MPI_UNSIGNED, MPI_ANY_SOURCE, 11, MPI_COMM_WORLD, 
	       &status);
      if (n) {
	vector<key_type> recv(n);
	MPI_Recv(&recv[0], n, MPI_UNSIGNED_LONG_LONG, status.MPI_SOURCE, 12, 
		 MPI_COMM_WORLD, MPI_STATUS_IGNORE);
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


void pHOT::adjustTiming(double &keymake, double &exchange, double &overlap)
{
  keymake  = timer_keymake.getTime().getRealTime()*1.0e-6;
  exchange = timer_xchange.getTime().getRealTime()*1.0e-6;
  overlap  = timer_overlap.getTime().getRealTime()*1.0e-6;

  timer_keymake.reset();
  timer_xchange.reset();
  timer_overlap.reset();
}


void pHOT::debug_ts(const char* msg)
{
  struct timeval tv;
  struct timezone tz;

  for (int n=0; n<numprocs; n++) {
    if (n==myid) {
      
      gettimeofday(&tv, &tz);

      cout << "Process " << myid << ": " << msg 
	   << ", " << tv.tv_sec << "." 
	   << setw(6) << setfill('0') << tv.tv_usec 
	   << setfill(' ') << endl;
    }
    MPI_Barrier(MPI_COMM_WORLD);
  }
}

