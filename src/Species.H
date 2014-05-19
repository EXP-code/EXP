#ifndef _Species_H
#define _Species_H

#include <map>
#include <set>
#include <vector>

//! Typedef for the species key for the different types
typedef std::pair<unsigned short, unsigned short> speciesKey;

typedef std::set<speciesKey>           sKeySet;
typedef std::map<speciesKey, double>   sKeyDmap;
typedef std::map<speciesKey, sKeyDmap> sKey2Dmap;
typedef std::map<speciesKey, unsigned> sKeyUmap;
typedef std::map<speciesKey, sKeyUmap> sKey2Umap;
typedef std::map<speciesKey, std::vector<double> > sKeyvDmap;

/**
   Wrap a union for converting a speciesKey to an int
*/
class KeyConvert 
{

private:
  union {
    unsigned short us[2];
    int i;
  } c;

public:

  //! Constructor taking a speciesKey
  KeyConvert(const speciesKey& k) { c.us[0] = k.first; c.us[1] = k.second; }

  //! Constructor taking a converted int
  KeyConvert(const int i)         { c.i = i; }

  //! Copy constructor
  KeyConvert(const KeyConvert& p) { c = p.c; }

  //! Get the converted int
  int        getInt() { return c.i; }
  
  //! Get the converter speciesKey
  speciesKey getKey() { return speciesKey(c.us[0], c.us[1]); }
};

#endif
