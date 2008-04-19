#include <iostream>
#include <fstream>
#include <iomanip>
#include <sstream>
#include <string>

#include "localmpi.h"
#include "ParticleFerry.H"
#include "pCell.H"
#include "pHOT.H"

unsigned pCell::bucket = 7;	// Target microscopic (collision) bucket size
unsigned pCell::Bucket = 64;	// Target macroscopic bucket size
unsigned pCell::nbits = 16;	// Number of bits per dimension

string printKey(key_type p)
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


pCell::pCell(pHOT* tr) : tree(tr), isLeaf(true)
{
  owner = myid;
				// I am the root node
  parent = 0;
  mykey = 1;
  level = 0;
  maxplev = 0;
				// My body mask
  mask = mykey << 3*(nbits - level);
				// Initialize state
  state = vector<double>(10, 0.0);

  tree->frontier[mykey] = this;	// Root is born on the frontier
}

pCell::pCell(pCell* mom, unsigned id) : 
  tree(mom->tree), parent(mom), isLeaf(true)
{
  owner = myid;
				// My map key
  mykey = (parent->mykey << 3) + id;
				// My level
  level = parent->level + 1;
				// Maximum particle level
  maxplev = 0;
				// My body mask
  mask = mykey << 3*(nbits - level);
				// Initialize state
  state = vector<double>(10, 0.0);

  tree->frontier[mykey] = this;	// All nodes born on the frontier
}

pCell::~pCell()
{
  // Recursively kill all the cells
  map<unsigned, pCell*>::iterator it;
  for (it=children.begin(); it!=children.end(); it++) delete it->second;
}

unsigned pCell::childId(key_type key)
{
  key_type id = key - mask;
  id = id >> 3*(nbits - 1 - level);
  return id;
}


pCell* pCell::Add(const key_pair& keypair)
{
  key_type key=keypair.first, key2;

				// Check that this key belongs to this branch
  key_type sig = (key_type)(key - mask) >> 3*(nbits-level);

  if (sig!=0) {

    if (parent == 0) {
      cout << "Process " << myid << ": level=" << level 
	   << " key=" << hex << key << endl
	   << " sig=" << hex << sig << endl << dec;
				// Get the particle info
      key_indx::iterator p =
	lower_bound(tree->keybods.begin(), tree->keybods.end(), keypair, ltPAIR());

      while (p->first==keypair.first && p->second==keypair.second ) {
	cout << "pos=";
	for (int k=0; k<3; k++) 
	  cout << setw(18) << tree->cc->Particles()[p->second].pos[k];
	cout << endl;
	p++;
      }
    }

    return parent->Add(keypair);
  }

  
  if (isLeaf && keys.find(keypair)==keys.end()) {

				// I am still a leaf . . .
    if (bods.size() < bucket || level+1==nbits) {
      keys.insert(keypair);
      tree->bodycell[key] = mykey;
      bods.push_back(keypair.second);
      maxplev = max<int>(maxplev, tree->cc->Particles()[keypair.second].level);
      
      return this;
    }
    
				// I need to make leaves and become a branch
    key_set::iterator n;
    for (n=keys.begin(); n!=keys.end(); n++) {
      key2 = childId(n->first);
      if (children.find(key2) == children.end())
	children[key2] = new pCell(this, key2);
      
      children[key2]->Add(*n);
    }
				// Erase my list
    keys.clear();
    bods.clear();
				// Erase my key from the frontier
    tree->frontier.erase(mykey);
				// I'm a branch now . . .
    isLeaf = false;
  }

				// Now add the *new* key
  key2 = childId(key);
  if (children.find(key2) == children.end())
    children[key2] = new pCell(this, key2);

  return children[key2]->Add(keypair);
}

bool pCell::Remove(const key_pair& keypair)
{
  bool ret = false;
				// Sanity check
  if (isMine(keypair.first)) {
				// Remove keypair from cell list
    keys.erase(keypair);
				// Remove from body-cell list
    tree->bodycell.erase(keypair.first);
				// Remove the key-body entry
    key_indx::iterator p = 
      find(tree->keybods.begin(), tree->keybods.end(), keypair);
    if (p==tree->keybods.end()) {
      cout << "pCell::error: mising keypair entry in keybods" << endl;
    } else {
      tree->keybods.erase(p);
    }
				// Remove the index from the cell body list
    sort(bods.begin(), bods.end());
    vector<unsigned>::iterator 
      ib = find(bods.begin(), bods.end(), keypair.second);
    if (ib!=bods.end()) bods.erase(ib);

    // Remove this cell if it is now empty
    if (bods.size()==0) {
      // Find the parent delete the cell from the parent list
      map<unsigned, pCell*>::iterator ic;
      for (ic=parent->children.begin(); ic!=parent->children.end(); ic++) {
	if (ic->second == this) {
	  parent->children.erase(ic);
	  break;
	}
      }
      // Remove it from the frontier
      tree->frontier.erase(mykey);
      ret = true;
    }
  } else {
    cout << "Process " << myid 
	 << ": pCell::Remove: body not in my cell!" << endl;
  }

  return ret;
}

pCell* pCell::findNode(const key_type& key)
{
				// Check that this key belongs to this branch
  key_type sig = (key_type)(key - mask) >> 3*(nbits-level);
  
  if (sig!=0) {
    
    if (parent == 0) {
      cout << "pHOT::findNode: impossible condition, process " 
	   << myid << ": level=" << level 
	   << " key=" << hex << key << endl
	   << " sig=" << hex << sig << endl << dec;
    }

    return parent->findNode(key);
  }

				// You found me!
  if (isLeaf) return this;
				// Which child
  key_type key2 = childId(key);
				// Not in my tree?
  if (children.find(key2)==children.end()) return 0;

				// Look for node amongst children
  return children[key2]->findNode(key);
}
 

void pCell::zeroState()
{
  count = 0;
  for (int k=0; k<10; k++) state[k] = 0.0;
  map<unsigned, pCell*>::iterator it = children.begin();
  for (; it != children.end(); it++) it->second->zeroState();
}

void pCell::accumState()
{
  unsigned indx, count = bods.size();
  for (unsigned j=0; j<count; j++) {
    indx = bods[j];
    state[0] += tree->cc->Particles()[indx].mass;
    for (int k=0; k<3; k++) {
      state[1+k] += tree->cc->Particles()[indx].mass * 
	tree->cc->Particles()[indx].vel[k]*tree->cc->Particles()[indx].vel[k];
      state[4+k] += tree->cc->Particles()[indx].mass * tree->cc->Particles()[indx].vel[k];
      state[7+k] += tree->cc->Particles()[indx].mass * tree->cc->Particles()[indx].pos[k];
    }
  }
  if (parent) parent->accumState(count, state);
}


void pCell::accumState(unsigned _count, vector<double>& _state)
{
  count += _count;
  for (int k=0; k<10; k++) state[k] += _state[k];
  if (parent) parent->accumState(_count, _state);
}


void pCell::Find(key_type key, unsigned& curcnt, unsigned& lev,
		 vector<double>& st)
{
  // Check to see if this key belongs to one of the children
  //
  key_type cid = key - mask;
  cid = cid >> 3*(nbits - 1 - level);

  map<unsigned, pCell*>::iterator it = children.begin();
  for (; it != children.end(); it++) {
    if (cid == it->first) {
      it->second->Find(key, curcnt, lev, st);
      return;
    }
  }

  // Return the values from this cell
  //
  curcnt = count;
  lev    = level;
  st     = state;

  return;
}

double pCell::Mass()
{
  return state[0]; 
}

void pCell::MeanPos(double &x, double &y, double& z)
{
  if (state[0]<=0.0) {
    x = y = z = 0.0;
    return;
  }
  x = state[7]/state[0];
  y = state[8]/state[0];
  z = state[9]/state[0];
}

void pCell::MeanPos(vector<double> &p)
{
  p = vector<double>(3, 0);
  if (state[0]<=0.0) return;
  for (int k=0; k<3; k++) p[k] = state[7+k]/state[0];
}

void pCell::MeanVel(double &u, double &v, double& w)
{
  if (state[0]<=0.0) {
    u = v = w = 0.0;
    return;
  }
  u = state[4]/state[0];
  v = state[5]/state[0];
  w = state[6]/state[0];
}

void pCell::MeanVel(vector<double> &p)
{
  p = vector<double>(3, 0);
  if (state[0]<=0.0) return;
  for (int k=0; k<3; k++) p[k] = state[4+k]/state[0];
}

void pCell::KE(double &total, double &dispr)
{
  total = 0.0;
  dispr = 0.0;

  if (state[0]>0.0) {
    for (int k=0; k<3; k++) {
      total += 0.5*state[1+k];
      dispr += 0.5*(state[1+k] - state[4+k]*state[4+k]/state[0]);
    }

    if (count<2) dispr=0.0;

    // DEBUG
    //
    static int cnt = 0;
    if (dispr<0.0) {
      ostringstream sout;
      sout << "pCell_tst." << myid << "." << cnt++;
      ofstream out(sout.str().c_str());
      out << "# number=" << count << endl;
      for (unsigned i=0; i<10; i++) 
	out << setw(3) << i << setw(15) << state[i] << endl;
    }
    
    // Return energy per unit mass
    //
    total /= state[0];
    dispr /= state[0];
  }

}

void pCell::Vel(double &mass, vector<double>& v1, vector<double>& v2)
{
  mass = 0.0;
  v1 = vector<double>(3, 0.0);
  v2 = vector<double>(3, 0.0);

  if (isLeaf) {
    unsigned number = bods.size();
    for (unsigned i=0; i<number; i++) {
      unsigned indx = bods[i];
      for (int k=0; k<3; k++) {
	v1[k] += tree->cc->Particles()[indx].mass * 
	  tree->cc->Particles()[indx].vel[k];

	v2[k] += tree->cc->Particles()[indx].mass * 
	  tree->cc->Particles()[indx].vel[k] * tree->cc->Particles()[indx].vel[k];
      }
      mass += tree->cc->Particles()[indx].mass;
    }
  }

}

double pCell::Volume()
{
  return tree->volume/((unsigned long long)1 << 3*level);
}

double pCell::Scale()
{
  return 1.0/((unsigned long long)1 << level);
}

pCell* pCell::findSampleCell()
{
  pCell *cur = this;		// Begin with this cell
  while(cur->count < Bucket) {
				// We are at the root
    if (cur->parent == 0) break;
				// Keep walking up the tree . . 
    cur = cur->parent;
  }
				// The answer.
  sample = cur;
}


Particle* pCell::Body(unsigned k)
{ 
  if (k<0 || k>=bods.size()) return 0;
  return &(tree->cc->Particles()[bods[k]]); 
}

unsigned pCell::remake_plev()
{
  maxplev = 0;
  vector<unsigned>::iterator it;
  for (it=bods.begin(); it!=bods.end(); it++) {
    maxplev = max<unsigned>(maxplev, tree->cc->Particles()[*it].level);
  }
  return maxplev;
}
