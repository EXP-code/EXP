#include <values.h>

#include <iostream>
#include <iomanip>
#include <fstream>
#include <sstream>
#include <vector>
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

  key_min = (key_type)1 << 48;
  key_max = (key_type)1 << 49;
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
  if (false) {			// Set to true for check
    for (unsigned k=0; k<3; k++) { 
      if (fabs((p[k]+offset[k])/sides[k])> 1.0) {
	cout << "Coordinate out of pbounds in pHOT::key: ";
	for (int l=0; l<3; l++) cout << setw(18) << p[l];
	cout << endl;
	return 0;
      }
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

  delete root;

  //
  // DEBUG
  //
  for (int n=0; n<numprocs; n++) {
    if (myid==n) {
      ofstream out("pHOT_storage.size", ios::app);
      if (out) {
	out << setw(18) << tnow
	    << setw(6)  << myid
	    << setw(12) << keybods.size()
	    << setw(12) << frontier.size()
	    << setw(12) << bodycell.size()
	    << endl;
	if (myid==numprocs-1) out << endl;
      }
    }
    MPI_Barrier(MPI_COMM_WORLD);
  }

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

  MPI_Status status;
				// Exchange boundary keys
  key_type headKey=0, tailKey=0, prevKey=0, nextKey=0;
  unsigned head_num=0, tail_num=0, next_num=0, prev_num=0;

				// Do the boundaries sequentially to prevent
				// inconstencies

  for (int n=1; n<numprocs; n++) {
				// Send the next node my tail value
				// to compare with its head
    if (myid==n-1) {

      if (keybods.size()) {

	key_indx::reverse_iterator it = keybods.rbegin();
	{
	  // Debug: check validity of key
	  if (bodycell.find(it->first) == bodycell.end()) {
	    cout << "Process " << myid << ": bad key=" 
		 << hex << it->first << dec
		 << " #cells=" << bodycell.size() << endl;
	  }
	}
	tailKey = bodycell[it->first];
				// Number of bodies in my tail cell
	tail_num = frontier[tailKey]->bods.size();
      } else {
	tailKey = 0;
	tail_num = 0;
      }

      MPI_Send(&tailKey, 1, MPI_UNSIGNED_LONG_LONG, n, 1000, MPI_COMM_WORLD);
      MPI_Send(&tail_num, 1, MPI_UNSIGNED, n, 1001, MPI_COMM_WORLD);

      MPI_Recv(&nextKey, 1, MPI_UNSIGNED_LONG_LONG, n, 1000, MPI_COMM_WORLD, &status);
      MPI_Recv(&next_num, 1, MPI_UNSIGNED, n, 1001, MPI_COMM_WORLD, &status);

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

      sort(keybods.begin(), keybods.end());
    }

				// Send the previous node my head value
				// to compare with its tail
    if (myid==n) {

      if (keybods.size()) {

	{
	  // Debug: check validity of key
	  if (bodycell.find(keybods.begin()->first) == bodycell.end()) {
	    cout << "Process " << myid << ": bad key=" 
		 << hex << keybods.begin()->first << dec
		 << " #cells=" << bodycell.size() << endl;
	  }
	}
	headKey = bodycell[keybods.begin()->first];
				// Number of bodies in my head cell
	head_num = frontier[headKey]->bods.size();
      } else {
	headKey = 0;
	head_num = 0;
      }

      MPI_Send(&headKey, 1, MPI_UNSIGNED_LONG_LONG, n-1, 1000, MPI_COMM_WORLD);
      MPI_Send(&head_num, 1, MPI_UNSIGNED, n-1, 1001, MPI_COMM_WORLD);

      MPI_Recv(&prevKey, 1, MPI_UNSIGNED_LONG_LONG, n-1, 1000, MPI_COMM_WORLD, &status);
      MPI_Recv(&prev_num, 1, MPI_UNSIGNED, n-1, 1001, MPI_COMM_WORLD, &status);

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

      sort(keybods.begin(), keybods.end());
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

  MPI_Barrier(MPI_COMM_WORLD);

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

      MPI_Status s;
      MPI_Recv(&cntlevN[0], MaxLev+1, MPI_UNSIGNED, n, 144, MPI_COMM_WORLD, &s);
      MPI_Recv(&kidlevN[0], MaxLev+1, MPI_UNSIGNED, n, 145, MPI_COMM_WORLD, &s);
      MPI_Recv(&maslevN[0], MaxLev+1, MPI_DOUBLE,   n, 146, MPI_COMM_WORLD, &s);
      MPI_Recv(&vollevN[0], MaxLev+1, MPI_DOUBLE,   n, 147, MPI_COMM_WORLD, &s);

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
	unsigned num = it->second->bods.size();
	double mass = 0.0;
	for (unsigned n=0; n<num; n++) {
	  mass += cc->particles[it->second->bods[n]].mass;
	  for (unsigned k=0; k<3; k++) {
	    tmp = cc->particles[it->second->bods[n]].pos[k];
	    mpos[k] += tmp;
	    vpos[k] += tmp*tmp;
	  }
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
	for (unsigned n=0; n<p->bods.size(); n++) {
	  mass += cc->particles[p->bods[n]].mass;
	  for (int k=0; k<3; k++) {
	    pos[k] += 
	      cc->particles[p->bods[n]].mass * 
	      cc->particles[p->bods[n]].pos[k];
	    vel[k] += 
	      cc->particles[p->bods[n]].mass * 
	      cc->particles[p->bods[n]].vel[k];
	    temp += 
	      cc->particles[p->bods[n]].mass * 
	      cc->particles[p->bods[n]].vel[k] *
	      cc->particles[p->bods[n]].vel[k];
	  }
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
#endif

  // DEBUG
  vector<unsigned long> erased;
  // END DEBUG


  vector<double> buffer1(3*num);
  vector<unsigned> buffer2(num);
  vector<key_type> buffer3(num);

  pf.ShipParticles(to, myid, num);

  key_type tkey;
  for (unsigned j=0; j<num; j++) {

    pf.SendParticle(cc->particles[p->bods[j]]);
    
				// Find the record and delete it
    tkey = cc->particles[p->bods[j]].key;
    key_indx::iterator it = 
      lower_bound(keybods.begin(), keybods.end(), tkey, ltULL());

    bool testme = false;
    while (it->first == tkey) {
      if (it->second == p->bods[j]) {
	keybods.erase(it);
	cc->particles.erase(p->bods[j]);
	// DEBUG
	erased.push_back(p->bods[j]);
	// END DEBUG
	testme = true;
	break;
      }
      it++;
    }

    if (!testme) {
      cerr << "Process " << myid << ": error! " << endl;
    }
  }
  
  // If this cell is not the root
  //
  if (p->parent) {
    // Delete this cell from the parent
    p->parent->children.erase(p->mykey & 0x7);
	
    // Delete this cell from the frontier
    frontier.erase(p->mykey);

    // Delete the cell altogether
    delete p;

  } else {			// Special treatment for root
    p->keys.clear();
    p->bods.clear();
  }

  // BEGIN DEBUG
  vector<unsigned long>::iterator iq;
  for (iq=erased.begin(); iq!=erased.end(); iq++) {
    if (cc->particles.find(*iq) != cc->particles.end())
      cout << "pHOT::sendCell proc=" << myid
	   << " found erased index=" << *iq << endl;
  }
  // END DEBUG

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
    cc->particles[part.indx] = part;
    if (part.key == 0) continue;
    if (part.key < key_min || part.key >= key_max) {
      cout << "Process " << myid << ": in recvCell, key=" 
	   << hex << part.key << endl << dec;
    }
    if (part.indx==0) cout << "pHOT::recvCell bad particle indx=0!" << endl;
    keybods.push_back(pair<key_type, unsigned>(part.key, part.indx));
    p = p->Add(pair<key_type, unsigned>(part.key, part.indx));
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

  MPI_Barrier(MPI_COMM_WORLD);

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

  MPI_Barrier(MPI_COMM_WORLD);

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
  unsigned maxlev = 0;
  key_cell::iterator it;
  for (it=frontier.begin(); it!=frontier.end(); it++)
    maxlev = max<unsigned>(maxlev, it->second->level);

  double vol1, vol;
  vol1 = volume/((key_type)1 << (3*maxlev));
  MPI_Allreduce(&vol1, &vol, 1, MPI_DOUBLE, MPI_MIN, MPI_COMM_WORLD);

  return vol;
}

double pHOT::maxVol()
{
  unsigned minlev = MAXINT;
  key_cell::iterator it;
  for (it=frontier.begin(); it!=frontier.end(); it++)
    minlev = min<unsigned>(minlev, it->second->level);

  double vol1, vol;
  vol1 = volume/((key_type)1 << (3*minlev));
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

    MPI_Status s;

    for (int n=1; n<numprocs; n++) {
      MPI_Recv(&num, 1, MPI_UNSIGNED, n, 61, MPI_COMM_WORLD, &s);
      vector<unsigned> lev1(num);
      MPI_Recv(&lev1[0], num, MPI_UNSIGNED, n, 62, MPI_COMM_WORLD, &s);
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
  MPI_Status s;
  map<unsigned long, Particle>::iterator it;
  
  MPI_Barrier(MPI_COMM_WORLD);

  volume = sides[0]*sides[1]*sides[2]; // Total volume of oct-tree region


  // BEGIN DEBUG
  vector<unsigned long> erased;
  // END DEBUG

  // checkBounds(2.0, "BEFORE repartition");

  //
  // Recompute keys and compute new partition
  //
  vector<key_type> keys, kbeg(numprocs), kfin(numprocs);
  for (it=cc->Particles().begin(); it!=cc->Particles().end(); it++) {
    it->second.key = getKey(&(it->second.pos[0]));
    keys.push_back(it->second.key);
    if (it->second.key == 0) {
      cout << "pHOT::Repartition [" << myid << "]: out of bounds,"
	   << " indx=" << setw(12) << it->second.indx
	   << " mass=" << setw(18) << it->second.mass
	   << " pos=";
      for (int k=0; k<3; k++)
	cout << setw(18) << it->second.pos[k];
      for (int k=0; k<3; k++)
	cout << setw(18) << it->second.vel[k];
      cout << endl;
    }
  }
  partitionKeys(keys, kbeg, kfin);

  //
  // Nodes compute send list
  //
  vector<key_type> loclist(numprocs+1);
  for (unsigned i=0; i<numprocs; i++) 
    loclist[i] = kbeg[i];
  loclist[numprocs] = kfin[numprocs-1];	// End point for binary search

  MPI_Barrier(MPI_COMM_WORLD);

  vector<unsigned> sndlist1(numprocs*numprocs, 0);
  vector<unsigned> sendlist(numprocs*numprocs);
  vector< vector<unsigned> > bodylist(numprocs);
  unsigned t;
  for (it=cc->Particles().begin(); it!=cc->Particles().end(); it++) {
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
	cout << "--------------------------------------------------------" << endl;
	cout << "Process " << myid << ": Tcnt=" << Tcnt << " Fcnt=" << Fcnt << endl;
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
  Partstruct *psend, *precv = 0;
  int ierr;

  if (Tcnt) psend = new Partstruct [Tcnt];
  if (Fcnt) precv = new Partstruct [Fcnt];

  //
  // Exchange particles between processes
  //
  for (int frID=0; frID<numprocs; frID++) {
    for (int toID=0; toID<numprocs; toID++) {
				// 
				// Current process sends particles (non-blocking)
      if (myid==frID) {		// 
	unsigned To = sendlist[numprocs*frID+toID];
	if (To) {
	  for (unsigned i=0; i<To; i++) {
	    pf.Particle_to_part(psend[ps+i], cc->Particles()[bodylist[toID][i]]);
	    cc->Particles().erase(bodylist[toID][i]);
	    erased.push_back(bodylist[toID][i]);
	  }
	  for (unsigned i=0; i<To; i++) {
	    if (psend[ps+i].indx == 0) {
	      cerr << "pHOT::Repartition[" << debug_ctr
		   << "]: SEND from=" << frID << " to=" << toID << " #=" << To
		   << " ps=" << ps+i << " mass=" << scientific << psend[ps+i].mass << endl;
	    }
	  }
	  if ( (ierr=MPI_Send(&psend[ps], To, ParticleFerry::Particletype, toID, 49, 
			      MPI_COMM_WORLD)) != MPI_SUCCESS)
	    {
	      cout << "Process " << myid << ": error in Reparition sending "
		   << To << " particles to #" << toID << " ierr=" << ierr << endl;
	    }
	  ps += To;
	}
      }
				// 
				// Current process receives particles (blocking)
      if (myid==toID) {		// 
	unsigned From = sendlist[numprocs*frID+toID];
	if (From) {
	  if ( (ierr=MPI_Recv(&precv[pr], From, ParticleFerry::Particletype, frID, 49, 
			      MPI_COMM_WORLD, &s)) != MPI_SUCCESS)
	    {
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
  // Buffer size sanity check
  //
  if (true) {			// Enabled if true

    if (Tcnt-ps) {
      cout << "Process " << myid 
	   << ": [Tcnt] found <" << ps << "> but expected <" 
	   << Tcnt << ">" << endl;
    }
    
    if (Fcnt-pr) {
      cout << "Process " << myid 
	   << ": [Fcnt] found <" << pr << "> but expected <" 
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
	       << " pr=" << pr+i << " mass=" << scientific << precv[pr+i].mass << endl;
	}
      }
      pr += From;
    }
  }

  debug_ctr++;

  // END DEBUG

  Particle part;
  for (unsigned i=0; i<Fcnt; i++) {
    pf.part_to_Particle(precv[i], part);
    cc->Particles()[part.indx] = part;
  }
  cc->nbodies = cc->particles.size();
      
  //
  // Clean up temporary body storage
  //
  if (Tcnt) delete [] psend;
  if (Fcnt) delete [] precv;

  //
  // Remake key body index
  //
  keybods.clear();
  unsigned oab1=0, oab=0;
  map<unsigned long, Particle>::iterator n;
  for (n=cc->Particles().begin(); n!=cc->Particles().end(); n++) {
    if (n->second.key==0) {
      oab1++;
      continue;
    }
    if (n->second.indx==0) cout << "pHOT::Repartition bad particle indx=0!" << endl;
    keybods.push_back(pair<key_type, unsigned>(n->second.key, n->second.indx));
  }
  sort(keybods.begin(), keybods.end());

  // checkBounds(2.0, "AFTER repartition");

  MPI_Reduce(&oab1, &oab, 1, MPI_UNSIGNED, MPI_SUM, 0, MPI_COMM_WORLD);
  if (myid==0 && oab)
    cout << endl << "pHOT::Repartition: " << oab << " out of bounds" << endl;

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

  MPI_Barrier(MPI_COMM_WORLD);

  if (true) {
    checkIndices();
    if (!checkKeybods())
      cout << "Process " << myid 
	   << ": particle key not in keybods list at T=" << tnow << endl;
  }
}

void pHOT::Rectify()
{
  MPI_Status s;
  map<unsigned long, Particle>::iterator it;
  
  MPI_Barrier(MPI_COMM_WORLD);

  volume = sides[0]*sides[1]*sides[2]; // Total volume of oct-tree region


  // BEGIN DEBUG
  vector<unsigned long> erased;
  // END DEBUG

  // checkBounds(2.0, "BEFORE Rectify");

  //
  // Recompute keys
  //
  vector<unsigned long> changelist;
  key_type newkey, celkey;
  unsigned inown=0;
  for (it=cc->Particles().begin(); it!=cc->Particles().end(); it++) {
    newkey = getKey(&(it->second.pos[0]));
    if (newkey != it->second.key) {
      // Check that this key belongs to the current cell
      // If so, do nothing; otherwise, delete from current cell, and 
      // delete from bodycell sturcture, put it on the changelist
      if (bodycell.find(it->second.key) != bodycell.end())
	celkey = bodycell[it->second.key];
      else
	cout << "Process " << myid << ": pHOT::Rectify logic error in bodycell"
	     << endl;
      pCell *c;
      if (frontier.find(celkey) != frontier.end())
	c = frontier[celkey];
      else
	cout << "Process " << myid << ": pHOT::Rectify logic error in frontier"
	     << endl;
      if ( !c->isMine(newkey) ) {
				// Look for the right cell in our tree
	pCell *c2 = c->findNode(newkey);
	if (c2) {		// It IS in our tree, add the body
	  c2->bods.push_back(it->first);
	  bodycell[newkey] = c2->mykey;
	  keybods.push_back(pair<key_type, unsigned>(newkey, it->first));
	  inown++;
	} else {		// It is NOT in our tree
	  changelist.push_back(it->first);
	} 
				// Remove the key from the index list
	keybods.erase(lower_bound(keybods.begin(), keybods.end(),
				  it->second.key, ltULL()));
				// Remove the key from the cell list
	bodycell.erase(it->first);
				// Assign the new key
	it->second.key = newkey;
      }
    }
    
    if (it->second.key == 0) {
      cout << "pHOT::Rectify [" << myid << "]: out of bounds,"
	   << " indx=" << setw(12) << it->second.indx
	   << " mass=" << setw(18) << it->second.mass
	   << " pos=";
      for (int k=0; k<3; k++)
	cout << setw(18) << it->second.pos[k];
      for (int k=0; k<3; k++)
	cout << setw(18) << it->second.vel[k];
      cout << endl;
    }
  }
  
  MPI_Barrier(MPI_COMM_WORLD);

  //
  // Accumulate global list by broadcast
  //
  vector<key_type> change_key;
  vector<unsigned long> change_indx;
  vector<unsigned short> old_own;
  vector<unsigned short> new_own, new_own0;
  vector<unsigned short> checkme, checkme0;
  unsigned cnt;
  
  for (int n=0; n<numprocs; n++) {
    if (n==myid) cnt=changelist.size();
    MPI_Bcast(&cnt, 1, MPI_UNSIGNED, n, MPI_COMM_WORLD);
    if (cnt) {
      vector<unsigned long> indx(cnt);
      if (n==myid) indx = changelist;
      MPI_Bcast(&indx[0], cnt, MPI_UNSIGNED_LONG, n,  MPI_COMM_WORLD);
      change_indx.insert(change_indx.end(), indx.begin(), indx.end());

      vector<unsigned long long> key(cnt);
      if (n==myid) {
	for (unsigned k=0; k<cnt; k++)
	  key[k] = cc->Particles()[changelist[k]].key;
      }
      MPI_Bcast(&key[0], cnt, MPI_UNSIGNED_LONG_LONG, n, MPI_COMM_WORLD);
      change_key.insert(change_key.end(), key.begin(), key.end());

      for (unsigned k=0; k<cnt; k++) {
	old_own.push_back(n+1);
	new_own.push_back(0);
	checkme.push_back(0);
      }
    }
  }

  //
  // Each node looks for owner
  //
  key_cell keycache;
  unsigned csize = change_indx.size();
  for (unsigned n=0; n<csize; n++) {
    if (old_own[n]-1==myid) continue;
    pCell *c = root->findNode(change_key[n]);
    if (c) {
      /*
      cout << "Process " << myid << ": found match key=" 
	   << hex << change_key[n] << dec << " whose owner was "
	   << old_own[n]-1 << endl;
      */
      new_own[n] = myid+1;
      checkme[n] = 1;
      keycache[change_key[n]] = c;
    }
  }

  //
  // Distribute new list to all processes
  //
  new_own0 = new_own;
  MPI_Allreduce(&new_own[0], &new_own0[0], csize, MPI_UNSIGNED_SHORT, MPI_SUM, 
		MPI_COMM_WORLD);

  checkme0 = checkme;
  MPI_Allreduce(&checkme[0], &checkme0[0], csize, MPI_UNSIGNED_SHORT, MPI_SUM, 
		MPI_COMM_WORLD);

  unsigned sum=0;
  for (unsigned k=0; k<csize; k++) sum += checkme0[k];

  //
  // Do the sanity check
  //
  if (true) {
    if (myid==0) {
      bool check_ok = true;
      unsigned bad_val, unclaimed=0, overclaimed=0;
      for (unsigned n=0; n<csize; n++) {
	if (checkme0[n] != 1) {
	  if (check_ok) bad_val = checkme0[n];
	  else          bad_val = max<unsigned>(checkme0[n], bad_val);
	  check_ok = false;
	  if (checkme0[n]==0) unclaimed++;
	  if (checkme0[n] >1) overclaimed++;
	}
      }
      if (!check_ok) {
	cout << endl 
	     << "pHOT::Rectify: error in check, largest val=" << bad_val
	     << ", " << unclaimed << " unclaimed, " 
	     << overclaimed << " overclaimed" << endl;
      }
    }
  }

  //
  // Make the shipping lists
  //
  vector<unsigned> sendlist(numprocs*numprocs, 0);
  vector< vector<unsigned> > bodylist(numprocs);
  for (unsigned k=0; k<csize; k++) {
    if (old_own[k]-1 == myid)
      bodylist[new_own[k]-1].push_back(change_indx[k]);
    sendlist[numprocs*(old_own[k]-1) + new_own[k]-1]++;
  }
  
  unsigned Tcnt=0, Fcnt=0;
  for (int i=0; i<numprocs; i++) {
    Tcnt += sendlist[numprocs*myid + i   ];
    Fcnt += sendlist[numprocs*i    + myid];
  }

				// DEBUG OUTPUT
  if (false) {			// If true, write send and 
				// receive list for each node
    for (int n=0; n<numprocs; n++) {
      if (myid==n) {
	cout << "--------------------------------------------------------" << endl;
	cout << "Process " << myid << ": Tcnt=" << Tcnt << " Fcnt=" << Fcnt << endl;
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
  Partstruct *psend, *precv = 0;
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
	    pf.Particle_to_part(psend[ps+i], cc->Particles()[bodylist[toID][i]]);
	    cc->Particles().erase(bodylist[toID][i]);
	    erased.push_back(bodylist[toID][i]);
	  }
	  for (unsigned i=0; i<To; i++) {
	    if (psend[ps+i].indx == 0) {
	      cerr << "pHOT::Rectify[" << debug_ctr
		   << "]: SEND from=" << frID << " to=" << toID << " #=" << To
		   << " ps=" << ps+i << " mass=" << scientific << psend[ps+i].mass << endl;
	    }
	  }
	  if ( (ierr=MPI_Send(&psend[ps], To, ParticleFerry::Particletype, toID, 49, 
			      MPI_COMM_WORLD)) != MPI_SUCCESS)
	    {
	      cout << "Process " << myid << ": error in Reparition sending "
		   << To << " particles to #" << toID << " ierr=" << ierr << endl;
	    }
	  ps += To;
	}
      }
				// 
				// Current process receives particles 
				// (blocking)
      if (myid==toID) {		// 
	unsigned From = sendlist[numprocs*frID+toID];
	if (From) {
	  if ( (ierr=MPI_Recv(&precv[pr], From, ParticleFerry::Particletype, frID, 49, 
			      MPI_COMM_WORLD, &s)) != MPI_SUCCESS)
	    {
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
  // Buffer size sanity check
  //
  if (true) {			// Enabled if true

    if (Tcnt-ps) {
      cout << "Process " << myid 
	   << ": [Tcnt] found <" << ps << "> but expected <" 
	   << Tcnt << ">" << endl;
    }
    
    if (Fcnt-pr) {
      cout << "Process " << myid 
	   << ": [Fcnt] found <" << pr << "> but expected <" 
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
	  cerr << "pHOT::Rectify[" << debug_ctr 
	       << "]: RECV to=" << myid << " from=" << id << " #=" << From
	       << " pr=" << pr+i << " mass=" << scientific << precv[pr+i].mass << endl;
	}
      }
      pr += From;
    }
  }

  debug_ctr++;

  // END DEBUG

  Particle part;
  for (unsigned i=0; i<Fcnt; i++) {
    pf.part_to_Particle(precv[i], part);
    cc->Particles()[part.indx] = part;
    keybods.push_back(pair<key_type, unsigned>(part.key, part.indx));
    keycache[part.key]->bods.push_back(part.indx);
    bodycell[part.key] = keycache[part.key]->mykey;
  }
  sort(keybods.begin(), keybods.end());
  cc->nbodies = cc->particles.size();
      
  //
  // Clean up temporary body storage
  //
  if (Tcnt) delete [] psend;
  if (Fcnt) delete [] precv;

  // checkBounds(2.0, "AFTER Rectify");

  if (true) {			// Sanity checks for bad particle indices
    bool ok = true;		// and bad particle counts
    map<unsigned long, Particle>::iterator ip;
    for (ip=cc->particles.begin(); ip!=cc->particles.end(); ip++) {
      if (ip->second.indx==0) {
	cout << "pHOT::Rectify BAD particle in proc=" << myid
	     << " key=" << ip->second.key << endl;
	ok = false;
      }
    }
    if (ok) cout << "pHOT::Rectify: leaving with good indx" << endl;

    if (cc->particles.size() != cc->nbodies) 
      cout << "pHOT::Rectify: leaving with # mismatch!" << endl;

    int nbodies1 = cc->nbodies, nbodies0=0;
    MPI_Reduce(&nbodies1, &nbodies0, 1, MPI_INT, MPI_SUM, 0, MPI_COMM_WORLD);
    if (myid==0) {
      if (nbodies0 != cc->nbodies_tot)
	cout << "pHOT::Rectify: leaving with total # mismatch!" << endl;
    }

  }

  MPI_Barrier(MPI_COMM_WORLD);

  // #ifdef DEBUG
  checkDupes();
  checkIndices();
  // #endif
}

void pHOT::checkDupes()
{
  MPI_Status s;

  vector<unsigned> indices;
  for (key_cell::iterator it=frontier.begin(); it!=frontier.end(); it++)
    indices.insert(indices.end(), it->second->bods.begin(), it->second->bods.end());

  sort(indices.begin(), indices.end());
  
  unsigned dup = 0;
  for (unsigned n=1; n<indices.size(); n++) {
    if (indices[n-1] == indices[n]) dup++;
  }

  if (dup)
    cout << "Process " << myid << ": pHOT::checkDupes, dup=" << dup << endl;
}


void pHOT::checkIndices()
{
  MPI_Status s;

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
      MPI_Recv(&icnt, 1, MPI_UNSIGNED, n, 392, MPI_COMM_WORLD, &s);
      vector<unsigned> get(icnt);
      MPI_Recv(&get[0], icnt, MPI_UNSIGNED, n, 393, MPI_COMM_WORLD, &s);

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
	cout << "pHOT::checkIndices ERROR: found=" << indices[n-1] 
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
      cout << "pHOT::checkIndices ERROR: time=" << tnow
	   << " index count=" << indices.size() 
	   << " body count=" << cc->nbodies_tot << endl;
  }
}

bool pHOT::checkKeybods()
{
  bool ok = true;
  map<unsigned long, Particle>::iterator n = cc->Particles().begin();
  for (; n!=cc->Particles().end(); n++) {
    key_indx::iterator it = 
      lower_bound(keybods.begin(), keybods.end(), n->second.key, ltULL());
    if (it==keybods.end()) {
      ok = false;
    }
  }

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

void pHOT::partitionKeys(vector<key_type>& keys, 
			 vector<key_type>& kbeg, vector<key_type>& kfin)
{
				// Sort the keys
  sort(keys.begin(), keys.end());

				// Sample keys at desired rate
  unsigned srate = 1000;	// e.g. every srate^th key is sampled

				// Number of samples
  unsigned nsamp = keys.size()/srate;
				// Sanity checks
  if (nsamp < 10000/numprocs) nsamp = 10000/numprocs; // Too few samples
  if (nsamp > keys.size())    nsamp = keys.size();    // Too many for particle count
  if (nsamp) srate = keys.size()/nsamp;	              // Consistent sampling rate

  vector<key_type> keylist1, keylist;
  for (unsigned i=0; i<nsamp; i++)
    keylist1.push_back(keys[srate*(2*i+1)/2]);

				// Tree aggregation (merge sort) of the
				// entire key list
  parallelMerge(keylist1, keylist);


  if (myid==0) {

    for (unsigned i=0; i<numprocs-1; i++)
      kfin[i] = keylist[static_cast<unsigned>(keylist.size()*(i+1)/numprocs)];
    kfin[numprocs-1] = key_max;

    kbeg[0] = key_min;
    for (unsigned i=1; i<numprocs; i++)
      kbeg[i] = kfin[i-1];

    if (false) {	 // If true, print key ranges for each process
      cout << "--------------------------------------------------" << endl
	   << "partitionKeys: keys in list=" << keylist.size() << endl
	   << "--------------------------------------------------" << endl;
      for (int i=0; i<numprocs; i++)
	cout << setw(5) << i << hex << setw(15) << kbeg[i] 
	     << setw(15) << kfin[i] << dec << endl;
      cout << "--------------------------------------------------" << endl;
    }
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
    MPI_Send(&initl[0], n, MPI_UNSIGNED_LONG_LONG, myid-M2, 12, 
	     MPI_COMM_WORLD);
    return;
  }

  vector<key_type> data = initl;

  //
  // Retrieve the excess particles
  //
  if (myid + M2 < numprocs) {
    MPI_Recv(&n, 1, MPI_UNSIGNED, myid+M2, 11, MPI_COMM_WORLD, &status);
    vector<key_type> recv(n);
    MPI_Recv(&recv[0], n, MPI_UNSIGNED_LONG_LONG, myid+M2, 12, 
	     MPI_COMM_WORLD, &status);
				// data=data+new_data
    sortCombine(initl, recv, data);
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
      MPI_Send(&data[0], n, MPI_UNSIGNED_LONG_LONG, myid-M2, 12, 
	       MPI_COMM_WORLD);
      return;
    } else {
      MPI_Recv(&n, 1, MPI_UNSIGNED, myid+M2, 11, MPI_COMM_WORLD, &status);
      vector<key_type> recv(n);
      MPI_Recv(&recv[0], n, MPI_UNSIGNED_LONG_LONG, myid+M2, 12, 
	       MPI_COMM_WORLD, &status);
      //
      // The lower half sorts and loop again
      //
      sortCombine(data, recv, work);
      data = work;
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

