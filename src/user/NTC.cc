#include <iostream>
#include <iomanip>
#include <sstream>

#include "Particle.H"
#include "pCell.H"
#include "NTC.H"

// For verbose debugging
static const bool DEBUG_V = true;

// Maximum number of cached relative velocities
size_t   NTCitem::VelCrsSZ  = 32; 

// Minimum CrossSection x Velocity value
double   NTCitem::VelCrsMin = 1.0e-24;

void NTCitem::VelCrsTest()
{
  if (VelCrsSum.size()==0)
    std::cout << "Empty caller: " << caller << " : " 
	      << " p: "  << &VelCrsSum
	      << std::endl;

  if (VelCrsList.size()>0 || VelCrsSum.size()>0) {

    std::cout << std::string(70, '-') << std::endl
	      << "VelCrs structure: "  << this << " caller: " << caller 
	      << std::endl << std::string(70, '-') << std::endl;

    for (auto k : VelCrsSum) {
      sKeyPair p = k.first;
      std::ostringstream sout;
      sout << "(" << p.first.first << "," << p.first.second
	   << "|" << p.first.first << "," << p.first.second << ")";
      std::cout << std::setw(14) << sout.str() 
		<< std::setw(10) << VelCrsNum[p]
		<< "  " << k.second << std::endl;
    }
    std::cout << std::string(70, '-') << std::endl;
  } else {
    std::cout << std::string(70, '-') << std::endl
	      << "VelCrs no points: " << this << std::endl
	      << std::string(70, '-') << std::endl;
  }
}

std::map<sKeyPair, NTCitem::vcTup> NTCitem::VelCrsAvg()
{
  std::map<sKeyPair, NTCitem::vcTup> ret;
  
  for (auto k : VelCrsSum) ret[k.first] = k.second / VelCrsNum[k.first];

  if (DEBUG_V) VelCrsTest();

  return ret;
}

NTCitem::vcTup NTCitem::VelCrsAvg(sKeyPair indx)
{
  vcMap::iterator  it = VelCrsSum.find(indx);
  numMap::iterator jt = VelCrsNum.find(indx);

  if (it == VelCrsSum.end()) return vcTup(VelCrsMin, 0, 0);
  if (jt == VelCrsNum.end()) return vcTup(VelCrsMin, 0, 0);

  return it->second / jt->second;
}

void NTCitem::VelCrsAdd(vcMap& vals)
{
  for (auto it : vals) {
    if (std::get<0>(it.second)>VelCrsMin) VelCrsAdd(it.first, it.second);
  }
}

void NTCitem::VelCrsAdd(sKeyPair indx, const vcTup& val)
{
  // value below threshold
  //
  if (std::get<0>(val)<=VelCrsMin) return;

  // Initialize running average
  //
  if (VelCrsList.find(indx) == VelCrsList.end()) {
    VelCrsSum[indx] = vcTup(0, 0, 0);
    VelCrsNum[indx] = 0;
  }

  // Add new element
  //
  VelCrsList[indx].push_back(val);
  VelCrsSum [indx] += val;
  VelCrsNum [indx] += 1;

  // Push out the oldest element if deque is full
  //
  if (VelCrsList[indx].size()>VelCrsSZ) {
    VelCrsSum [indx] -= VelCrsList[indx].front();
    VelCrsNum [indx] -= 1;
    VelCrsList[indx].pop_front();
  }

}


NTCptr NTCdb::operator[](const key_type& k)
{
  NTCdata::iterator it = data.find(k);
  NTCptr ret;

  if (it == data.end()) {
    
    // Look for an initialization parent
    //
    key_type check = k >> 3;
    it = data.find(check);
    while (it == data.end() and check != 0) {
      check = check >> 3;
      it = data.find(check);
    }

    if (DEBUG_V && it == data.end()) {
      std::cout << "No parent found for " << k << std::endl;
    }

    // If none, create an empty item.  Otherwise, initialize from the
    // parent
    //
    if (it == data.end()) {
      data[k] = ret = NTCptr(new NTCitem);
    } else {
      data[k] = ret = NTCptr(new NTCitem(it->second));
    }

  } else {
    ret = it->second;
  }

  if (DEBUG_V) ret->setKey(k);

  return ret;
}
