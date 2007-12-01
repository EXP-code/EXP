#include <iostream>
#include <iomanip>
#include <fstream>
#include <sstream>
#include <vector>
#include <algorithm>
#include <cmath>

using namespace std;

#include "pHOT.H"

double pHOT::sides[] = {2.0, 2.0, 2.0};
double pHOT::offst[] = {1.0, 1.0, 1.0};
unsigned pHOT::neg_half = 0;

/*
  Constructor: initialize domain
*/
pHOT::pHOT(Component *C)
{
  cc = C;			// Register the calling component

  volume = sides[0]*sides[1]*sides[2];	// Total volume of oct-tree region
  root = 0;
				// For debugging
  key_min = (key_type)1 << 48;
  key_max = (key_type)1 << 49;
} 


pHOT::~pHOT()
{
  delete root;
}


unsigned long long pHOT::getKey(double *p)
{
  // Out of bounds?
  //
  for (unsigned k=0; k<3; k++) { 
    if (fabs((p[k]+offst[k])/sides[k])> 1.0) {
      cout << "Coordinate out of pbounds in pHOT::key: ";
      for (int l=0; l<3; l++) cout << setw(18) << p[l];
      cout << endl;
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
    bins[2-k] = (unsigned)floor( ((p[k]+offst[k])/sides[k]+neg_half)*factor );
  
  unsigned long long place = 1;
  unsigned long long _key = 0;
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

string pHOT::printKey(unsigned long long p)
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
  // Make the root
  //
  root = new pCell(this);

  //
  // Add the data
  //
  key_indx::iterator it;
  pCell* p = root;
  for (it=keybods.begin(); it!=keybods.end(); it++)  {
				// Ignore out of bounds particles
    if (it->first == 0) continue;
				// Sanity check . . .
    if (it->first < key_min || it->first >= key_max) {
      cout << "Process " << myid << ": in makeTree, key=" 
	   << hex << it->first << endl << dec;
    }
    p = p->Add(*it);		// Do the work
  }

  //
  // Adjust boundaries bodies to prevent cell duplication on the boundary
  //

  MPI_Status status;
				// Exchange boundary keys
  unsigned long long headKey, tailKey, prevKey, nextKey;
  unsigned head_num=0, tail_num=0, next_num=0, prev_num=0;

				// Do the boundaries sequentially to prevent
				// inconstencies

  for (int n=1; n<numprocs; n++) {
				// Send the next node my tail value
				// to compare with its head
    if (myid==n-1) {

      key_indx::iterator it = keybods.end(); it--;
      tailKey = bodycell[it->first];
      tail_num = frontier[tailKey]->bods.size(); // Number of bodies in my tail cell

      MPI_Send(&tailKey, 1, MPI_UNSIGNED_LONG_LONG, n, 1000, MPI_COMM_WORLD);
      MPI_Send(&tail_num, 1, MPI_UNSIGNED, n, 1001, MPI_COMM_WORLD);

      MPI_Recv(&nextKey, 1, MPI_UNSIGNED_LONG_LONG, n, 1000, MPI_COMM_WORLD, &status);
      MPI_Recv(&next_num, 1, MPI_UNSIGNED, n, 1001, MPI_COMM_WORLD, &status);

      if (tailKey == nextKey) {
	if (tail_num <= next_num) {
	  sendCell(tailKey, n, tail_num);
	} else {
	  recvCell(n, next_num);
	}
      }

      sort(keybods.begin(), keybods.end());
    }
				// Send the previous node my head value
				// to compare with its tail
    if (myid==n) {
      headKey = bodycell[keybods.begin()->first];
      head_num = frontier[headKey]->bods.size(); // Number of bodies in my head cell

      MPI_Send(&headKey, 1, MPI_UNSIGNED_LONG_LONG, n-1, 1000, MPI_COMM_WORLD);
      MPI_Send(&head_num, 1, MPI_UNSIGNED, n-1, 1001, MPI_COMM_WORLD);

      MPI_Recv(&prevKey, 1, MPI_UNSIGNED_LONG_LONG, n-1, 1000, MPI_COMM_WORLD, &status);
      MPI_Recv(&prev_num, 1, MPI_UNSIGNED, n-1, 1001, MPI_COMM_WORLD, &status);

      if (headKey == prevKey) {
	if (head_num < prev_num) {
	  sendCell(headKey, n-1, head_num);
	} else {
	  recvCell(n-1, prev_num);
	}
      }

      sort(keybods.begin(), keybods.end());
    }    
  }

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
    vollev[lev] += volume/(1<<(3*p->level));
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
	totvol  += volume/(1<<(3*it->second->level));

	cout << setw(4)  << myid
	     << setw(12) << hex << it->first << dec
	     << setw(8)  << it->second->level
	     << setw(18) << num
	     << setw(18) << mass
	     << setw(18) << mass/(volume/(1<<(3*it->second->level)));
	
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


void pHOT::testDump()
{
  key_indx::iterator it;

  for (int n=0; n<numprocs; n++) {

    if (myid==n) {
      for (it=keybods.begin(); it!=keybods.end(); it++) {
	cout << setw(12) << hex << it->first;
	for (unsigned k=0; k<3; k++) 
	  cout << setw(18) << dec << (cc->particles)[it->second].pos[k];
	cout << endl;
      }
    }
    MPI_Barrier(MPI_COMM_WORLD);
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
	double vol = volume/((unsigned long long)1 << (3*p->level)); 

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


void pHOT::sendCell(unsigned long long key, int to, unsigned num)
{
  pCell *p = frontier.find(key)->second;
  
#ifdef DEBUG
  cout << "Process " << myid << ": sending " << num 
       << " to " << to << endl;
#endif


  vector<double> buffer1(3*num);
  vector<unsigned> buffer2(num);
  vector<unsigned long long> buffer3(num);

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
    keybods.push_back(pair<key_type, unsigned>(part.key, part.indx));
    if (part.key == 0) continue;
    if (part.key < key_min || part.key >= key_max) {
      cout << "Process " << myid << ": in recvCell, key=" 
	   << hex << part.key << endl << dec;
    }
    p = p->Add(pair<unsigned long long, unsigned>(part.key, part.indx));
  }

  cc->nbodies = cc->particles.size();
}

void pHOT::makeState()
{
  // Each node start at root and walk down the three to zero the counts
  //
  root->zeroState();

  // March through the frontier and accumulate the counts
  //
  for (key_cell::iterator it=frontier.begin(); 
       it != frontier.end(); it++) it->second->accumState();
}


void pHOT::State(double *x, double& dens, double& temp,
		 double& velx, double& vely, double& velz)
{
  unsigned long long key = getKey(x);

  dens = temp = velx = vely = velz = 0.0;

  // Walk tree to get count
  //
  unsigned count, level;
  vector<double> state(5);
  root->Find(key, count, level, state);

  vector<double>   stt1(numprocs*5, 0);
  vector<double>   stt0(numprocs*5, 0);
  vector<unsigned> cnt1(numprocs, 0), lev1(numprocs, 0);
  vector<unsigned> cnt0(numprocs, 0), lev0(numprocs, 0);

  cnt1[myid] = count;
  lev1[myid] = level;
  for (int k=0; k<5; k++) stt1[myid*5+k] = state[k];


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
    for (int k=0; k<5; k++) state[k] = stt0[5*dlist[0].second + k];

				// Add others at the same level
    for (int n=1; n<numprocs; n++) {
      if (clv != dlist[n].first) break;
      for (int k=0; k<5; k++) state[k] += stt0[5*dlist[n].second + k];
      cnt++;
    }

				// Mass average
    if (state[0]>0.0) {
      double v2 = 0.0;
      for (int k=0; k<3; k++) {
	state[2+k] /= state[0];
	v2 += state[2+k] * state[2+k] * state[0];
      }
      

      dens = state[0] * (1 << (3*clv))/(volume*cnt);
      temp = 0.333333333333*(state[1] - v2)/state[0];
      velx = state[2];
      vely = state[3];
      velz = state[4];
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
    pt[2] = 0.5*sides[2] - offst[2];
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
    pt[1] = 0.5*sides[1] - offst[1];
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
    pt[0] = 0.5*sides[0] - offst[0];
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

  return volume/(1 << (3*maxlev));
}

double pHOT::medianVol()
{
  unsigned mlev, num;
  vector<unsigned> lev;

  key_cell::iterator it;
  for (it=frontier.begin(); it!=frontier.end(); it++) 
    lev.push_back(it->second->level);

  // MPI_Barrier(MPI_COMM_WORLD);

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

  return volume/((unsigned long long)1 << (3*mlev));
}

void pHOT::Repartition()
{
  MPI_Status s;
  unsigned nbod, Tcnt, Fcnt;
  vector<unsigned> From(numprocs, 0), To(numprocs, 0), Tlst;
  
  MPI_Barrier(MPI_COMM_WORLD);

  volume = sides[0]*sides[1]*sides[2]; // Total volume of oct-tree region


  //
  // Recompute keys
  //
  map<unsigned long, Particle>::iterator it = cc->Particles().begin();
  for (; it!=cc->Particles().end(); it++) {
    it->second.key = getKey(&(it->second.pos[0]));
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

  //
  // Compute the new partition
  //
  if (myid==0) {
    
    typedef pair<key_type, pair<unsigned, int> > NewPair;
    vector<NewPair> newmap;
				// 
				// Root's particles first
				// 
    map<unsigned long, Particle>::iterator n = cc->Particles().begin();
    for (; n!=cc->Particles().end(); n++)
      newmap.push_back(NewPair(n->second.key, 
			       pair<unsigned, int>(n->first, 0)));
    
				// 
				// Get remote particles
				// 
    for (int id=1; id<numprocs; id++) {
      MPI_Recv(&nbod, 1, MPI_UNSIGNED, id, 41, MPI_COMM_WORLD, &s);

      vector<unsigned long long> nkeys(nbod);
      vector<unsigned> nindx(nbod);

      MPI_Recv(&nkeys[0], nbod, MPI_UNSIGNED_LONG_LONG, id, 42, MPI_COMM_WORLD, &s);
      MPI_Recv(&nindx[0], nbod, MPI_UNSIGNED, id, 43, MPI_COMM_WORLD, &s);

      for (unsigned i=0; i<nbod; i++) {
	newmap.push_back(NewPair(nkeys[i], 
				 pair<unsigned, int>(nindx[i], id)));
      }
    }
    sort(newmap.begin(), newmap.end());

    unsigned newnum = newmap.size();
    unsigned num = newnum/numprocs;
    unsigned nnb = newnum - num*(numprocs-1);

    vector<NewPair>::iterator it = newmap.begin();

    vector< vector< vector<unsigned> > > tlst(numprocs);
    vector<unsigned> tcnt(numprocs, 0), fcnt(numprocs, 0);
    vector< vector<unsigned> > to(numprocs), from(numprocs);

    for (int i=0; i<numprocs; i++) {
      tlst[i] = vector< vector<unsigned> >(numprocs);
      to[i]   = vector<unsigned>(numprocs, 0);
      from[i] = vector<unsigned>(numprocs, 0);
    }

    for (unsigned i=0; i<nnb; i++) {
      if (it->second.second != 0) {
	tlst[it->second.second][0].push_back(it->second.first);
	tcnt[it->second.second]++;
	to[it->second.second][0]++;
	from[0][it->second.second]++;
	fcnt[0]++;
      }
      it++;
    }

    for (int id=1; id<numprocs; id++) {
      for (unsigned i=0; i<num; i++) {
	if (it->second.second != id) {
	  tlst[it->second.second][id].push_back(it->second.first);
	  tcnt[it->second.second]++;
	  to[it->second.second][id]++;
	  from[id][it->second.second]++;
	  fcnt[id]++;
	}
	it++;
      }
    }

    if (0) {

      cout << setw(60) << setfill('-') << '-' << endl << setfill(' ');
    
      for (int id=0; id<numprocs; id++) {
	if (tcnt[id]) {
	  for (int id2=0; id2<numprocs; id2++) {
	    if (tlst[id][id2].size()) {
	      cout << "#" << setw(3) << id 
		   << " sending " << setw(3) << tlst[id][id2].size()
		   << " to   " << setw(3) << id2 << endl;
	    }
	  }
	}
      }

      for (int id=0; id<numprocs; id++) {
	if (fcnt[id]) {
	  for (int id2=0; id2<numprocs; id2++) {
	    if (from[id][id2]) {
	      cout << "#" << setw(3) << id
		   << " getting " << setw(3) << from[id][id2]
		   << " from " << setw(3) << id2 << endl;
	    }
	  }
	}
      }
      
      cout << setw(60) << setfill('-') << '-' << endl << setfill(' ');
    }

    //
    // Ship lists
    //
    for (int id=1; id<numprocs; id++) {
      MPI_Send(&tcnt[id], 1, MPI_UNSIGNED, id, 44, MPI_COMM_WORLD);
      MPI_Send(&fcnt[id], 1, MPI_UNSIGNED, id, 45, MPI_COMM_WORLD);
    }

    vector<unsigned> tolist;	// Work vector

    for (int id=1; id<numprocs; id++) {
      if (tcnt[id]) {
	MPI_Send(&to[id][0], numprocs, MPI_UNSIGNED, id, 46, MPI_COMM_WORLD);
	tolist.clear();
	for (int id2=0; id2<numprocs; id2++) 
	  tolist.insert(tolist.end(), tlst[id][id2].begin(), tlst[id][id2].end());
	MPI_Send(&tolist[0], tcnt[id], MPI_UNSIGNED, id, 47, MPI_COMM_WORLD);
      }
      if (fcnt[id])
	MPI_Send(&from[id][0], numprocs, MPI_UNSIGNED, id, 48, MPI_COMM_WORLD);
    }
    
				// Root node's lists
    Tcnt = tcnt[0];
    Fcnt = fcnt[0];
    for (int id=0; id<numprocs; id++) 
      Tlst.insert(Tlst.end(), tlst[0][id].begin(), tlst[0][id].end());
    To = to[0];
    From = from[0];

  } else {

    nbod = cc->Number();
    MPI_Send(&nbod, 1, MPI_UNSIGNED, 0, 41, MPI_COMM_WORLD);

    vector<unsigned long long> Nkeys;
    vector<unsigned> Nindx;

    map<unsigned long, Particle>::iterator n = cc->Particles().begin();
    for (; n!=cc->Particles().end(); n++) {
      Nkeys.push_back(n->second.key);
      Nindx.push_back(n->second.indx);
    }
    MPI_Send(&Nkeys[0], nbod, MPI_UNSIGNED_LONG_LONG, 0, 42, MPI_COMM_WORLD);
    MPI_Send(&Nindx[0], nbod, MPI_UNSIGNED, 0, 43, MPI_COMM_WORLD);

    //
    // Receive lists
    //
    MPI_Recv(&Tcnt, 1, MPI_UNSIGNED, 0, 44, MPI_COMM_WORLD, &s);
    MPI_Recv(&Fcnt, 1, MPI_UNSIGNED, 0, 45, MPI_COMM_WORLD, &s);

    if (Tcnt) {
      MPI_Recv(&To[0], numprocs, MPI_UNSIGNED, 0, 46, MPI_COMM_WORLD, &s);
      Tlst = vector<unsigned>(Tcnt);
      MPI_Recv(&Tlst[0], Tcnt, MPI_UNSIGNED, 0, 47, MPI_COMM_WORLD, &s);
    }

    if (Fcnt)
      MPI_Recv(&From[0], numprocs, MPI_UNSIGNED, 0, 48, MPI_COMM_WORLD, &s);
  }

  MPI_Barrier(MPI_COMM_WORLD);

  unsigned ps=0, pr=0;
  vector<MPI_Request> send_requests, recv_requests;
  vector<Partstruct> psend(Tcnt), precv(Fcnt);
  MPI_Request r;

  //
  // Send particles
  //

  for (int id=0; id<numprocs; id++) {
    if (To[id]) {
      for (unsigned i=0; i<To[id]; i++) {
	pf.Particle_to_part(psend[ps+i], cc->Particles()[Tlst[ps+i]]);
	cc->Particles().erase(Tlst[ps+i]);
      }
      send_requests.push_back(r);
      MPI_Isend(&psend[ps], To[id], ParticleFerry::Particletype, id, 49, 
		MPI_COMM_WORLD, &send_requests.back());
      ps += To[id];
    }
  }

  //
  // Receive particles
  //
  
  for (int id=0; id<numprocs; id++) {
    if (From[id]) {
      recv_requests.push_back(r);
      MPI_Irecv(&precv[pr], From[id], ParticleFerry::Particletype, id, 49, 
		MPI_COMM_WORLD, &recv_requests.back());
      pr += From[id];
    }
  }


  //
  // Buffer size sanity check
  //

  if (1) {

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

  for (unsigned i=0; i<send_requests.size(); i++) {
    MPI_Wait(&send_requests[i], &s);
  }

  for (unsigned i=0; i<recv_requests.size(); i++) {
    MPI_Wait(&recv_requests[i], &s);
  }
  
  Particle part;
  for (unsigned i=0; i<Fcnt; i++) {
    pf.part_to_Particle(precv[i], part);
    cc->Particles()[part.indx] = part;
  }
      
  //
  // Remake key body index
  //
  keybods.clear();
  map<unsigned long, Particle>::iterator n = cc->Particles().begin();
  for (; n!=cc->Particles().end(); n++)
    keybods.push_back(pair<key_type, unsigned>(n->second.key, n->second.indx));
  sort(keybods.begin(), keybods.end());

  MPI_Barrier(MPI_COMM_WORLD);
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
