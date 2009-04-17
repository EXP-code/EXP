#ifndef _pHOT_types_H
#define _pHOT_types_H

#include "../config.h"

#include <map>
#include <set>

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

typedef unsigned long indx_type;
typedef unsigned long long key_type;
typedef pair<key_type, unsigned> key_pair;
typedef pair<key_type, key_type> key_item;
typedef pair<pCell*, indx_type> cell_indx;

				// For cell transaction list
typedef list<cell_indx> change_list;

				// A list of key/index pairs for particles
typedef set< pair<key_type, indx_type>, ltPAIR> key_indx;
				// Particle key points to cell key
				// multimap permits multiple identical pairs
typedef multimap<key_type, key_type> key_key;


/*
  Now, deal with the transition to ANSI undordered_maps from pre-ANSI
  hash_maps using the config variables
*/

#ifdef HAVE_CXX0X
//============================================================================
#include <unordered_map>
				// Cell key points to cell pointer
typedef std::unordered_map<key_type, pCell*, hashULL, eqULL> key_cell;
				// Particle index points to particle key
typedef std::unordered_map<unsigned, key_type, std::hash<indx_type> > indx_key;
//============================================================================

#else

#ifdef HAVE_TR1

//============================================================================
#include <tr1/unordered_map>
				// Cell key points to cell pointer
typedef std::tr1::unordered_map<key_type, pCell*, hashULL, eqULL> key_cell;
				// Particle index points to particle key
typedef std::tr1::unordered_map<unsigned, key_type, std::tr1::hash<indx_type> > indx_key;
//============================================================================

#else

#ifdef HAVE_HASH

//============================================================================
#include <ext/hash_map>
				// Cell key points to cell pointer
typedef hash_map<key_type, pCell*, hashULL, eqULL> key_cell;
				// Particle index points to particle key
typedef hash_map<unsigned, key_type, hash<indx_type> > indx_key;
//============================================================================

#else

//============================================================================
				// Cell key points to cell pointer
typedef std::map<key_type, pCell*, eqULL> key_cell;
				// Particle index points to particle key
typedef std::map<unsigned, key_type> indx_key;
//============================================================================

#endif // HAVE_HASH

#endif // HAVE_TR1

#endif // HAVE_CXX0X

#endif
