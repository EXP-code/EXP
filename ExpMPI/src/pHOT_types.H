#ifndef _pHOT_types_H
#define _pHOT_types_H

#include <map>
#include <set>
#include <ext/hash_map>
#include <ext/hash_set>

using namespace std;
using namespace __gnu_cxx;

struct eqULL {
  bool operator()(const unsigned long long __i1, 
		  const unsigned long long __i2) const
  { return (__i1==__i2); }
};

struct ltPAIR {
  bool operator()(const pair<unsigned long long, unsigned> __p1,
		  const pair<unsigned long long, unsigned> __p2) const
  { 
    if (__p1.first==__p2.first) return (__p1.second<__p2.second); 
    else return (__p1.first<__p2.first);
  }
};

struct hashULL {
  size_t operator()(unsigned long long __i) const
  { return static_cast<size_t>(__i); }
};

class pCell;
struct Partstruct;

typedef unsigned long long key_type;
typedef pair<unsigned long long, unsigned> key_pair;
typedef pair<unsigned long long, unsigned long long> key_item;
typedef pair<pCell*, unsigned> cell_indx;

				// For cell transaction list
typedef list<cell_indx> change_list;

				// A list of key/index pairs for particles
typedef set< pair<unsigned long long, unsigned>, ltPAIR> key_indx;
				// Cell key points to cell pointer
typedef hash_map<unsigned long long, pCell*, hashULL, eqULL> key_cell;
				// Particle key points to cell key
				// multimap permits multiple identical pairs
typedef multimap<unsigned long long, unsigned long long> key_key;
				// Particle index points to particle key
typedef hash_map<unsigned, unsigned long long, hash<unsigned> > indx_key;


#endif
