#ifndef _pHOT_types_H
#define _pHOT_types_H

#include <config_exp.h>

#include <map>
#include <set>

using namespace std;

class pCell;
struct Partstruct;

typedef unsigned long indx_type;
typedef unsigned long key_type;
typedef pair<key_type, indx_type> key_pair;
typedef pair<key_type, key_pair>  key_item;
typedef pair<pCell*, indx_type>   cell_indx;

struct eqULL {
  bool operator()(const key_type __i1, 
		  const key_type __i2) const
  { return (__i1==__i2); }
};

struct ltPAIR {
  bool operator()(const pair<key_type, indx_type> __p1,
		  const pair<key_type, indx_type> __p2) const
  { 
    if (__p1.first==__p2.first) return (__p1.second<__p2.second); 
    else return (__p1.first<__p2.first);
  }
};

struct hashULL {
  size_t operator()(unsigned long __i) const
  { return static_cast<size_t>(__i); }
};

				// For cell transaction list
typedef vector<cell_indx> change_list;

				// A list of key/index pairs for
				// particles These pairs must be
				// unique because the particle indices
				// are unique
typedef set<key_pair, ltPAIR> key_indx;
				// Particle key points to cell key
				// multimap permits multiple identical pairs
typedef multimap<key_type, key_pair> key_key;
				// For searching the key_key multimap
typedef key_key::iterator key2Itr;
typedef std::pair <key2Itr, key2Itr> key2Range;


				// For load balancing and partitioning
typedef pair<key_type, double> key_wght;


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
typedef std::unordered_map<indx_type, key_type, std::hash<indx_type> > indx_key;
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
