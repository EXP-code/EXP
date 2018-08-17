#ifndef _Species_H
#define _Species_H

#include <map>
#include <set>
#include <vector>

//@{
//! Typedefs for the species key
typedef std::pair<unsigned short, unsigned short>  speciesKey;
typedef std::set<speciesKey>                       sKeySet;
typedef std::map<speciesKey, double>               sKeyDmap;
typedef std::map<speciesKey, sKeyDmap>             sKey2Dmap;
typedef std::map<speciesKey, unsigned>             sKeyUmap;
typedef std::map<speciesKey, sKeyUmap>             sKey2Umap;
typedef std::map<speciesKey, std::vector<double> > sKeyvDmap;
typedef std::pair<speciesKey, speciesKey>          sKeyPair;
//@}

/**
   Wrap a union for converting a speciesKey to an int
*/
class KeyConvert 
{

private:
  // The union datum
  union {
    unsigned short us[2];
    int i;
  } c;

public:

  //! Null constructor
  KeyConvert() { c.i = -1; }

  //! Constructor taking a speciesKey (enter values into the union)
  KeyConvert(const speciesKey& k) { c.us[0] = k.first; c.us[1] = k.second; }

  //! Constructor taking a converted int (enter value into the union)
  KeyConvert(const int i)         { c.i = i; }

  //! Copy constructor
  KeyConvert(const KeyConvert& p) { c = p.c; }

  //! Get the converted int
  int        getInt() { return c.i; }
  
  //! Get the converter speciesKey
  speciesKey getKey() { return speciesKey(c.us[0], c.us[1]); }

  //! Get Z
  unsigned short Z() { return c.us[0]; }

  //! Get C
  unsigned short C() { return c.us[1]; }

  //! Update C and return int
  int updateC(short unsigned C) { c.us[1] = C; return c.i; }
};

#endif
