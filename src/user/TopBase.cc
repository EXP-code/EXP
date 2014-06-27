#include <iostream>
#include <iomanip>
#include <fstream>
#include <sstream>
#include <string>
#include <cmath>

#include <TopBase.H>

void dbg(std::istringstream& s)
{
  std::cout << std::string(30, '-')   << std::endl
	    << "  Good: " << s.good() << std::endl
	    << "   EOF: " << s.eof()  << std::endl
	    << "  Fail: " << s.fail() << std::endl
	    << "   Bad: " << s.bad()  << std::endl
	    << std::string(30, '-')   << std::endl;
}

void TopBase::readData()
{
  char * val;
  if ( (val = getenv("TOPBASE_DATA")) == 0x0) {
    // if (myid==0)
      std::cout << "Could not find TOPBASE_DATA environment variable"
		<< " . . . exiting" << std::endl;
      // MPI_Abort(MPI_COMM_WORLD, 48);
  }

  std::string fileName(val);
  fileName.append("/topbase.");

  std::string csfile(fileName), swfile(fileName);

  csfile.append("cross");
  swfile.append("weight");

  std::ifstream in(csfile.c_str());
  while (in.good()) {
    std::string line;
    getline(in, line);

    // Skip "readable" info
    //
    if (line.find("="  ) != std::string::npos || 
	line.find("RYD" )!= std::string::npos  ) continue;

    // This should be line header for a (I, NZ, NE) triple
    //
    std::istringstream sin(line.c_str());

    TBptr data = TBptr(new TBline);

    // Parse the record
    //
    sin >> data->I;
    sin >> data->NZ;
    sin >> data->NE;
    sin >> data->ISLP;
    sin >> data->ILV;
    sin >> data->Eph;
    sin >> data->NP;

    // If header has been parsed, read the line data
    //
    if ( !sin.fail() && !sin.bad() && data->NP>0) {
      data->E.resize(data->NP);
      data->S.resize(data->NP);
      for (int i=0; i<data->NP; i++) {
	getline(in, line);
	std::istringstream sin(line.c_str());
	sin >> data->E[i];
	sin >> data->S[i];
      }

      // The ion key
      //
      iKey key(data->NZ, data->NE);

      // Add new datum
      //
      ions[key][data->ISLP][data->ILV] = data;

    }
  }

  in.close();
  in.open(swfile.c_str());

  while (in.good()) {
    std::string line;
    getline(in, line);

    // Skip "readable" info
    //
    if (line.find("="  ) != std::string::npos || 
	line.find("RYD" )!= std::string::npos  ) continue;

    // This should be line header for a (I, NZ, NE) triple
    //
    std::istringstream sin(line.c_str());

    TBptr data = TBptr(new TBline);

    // Parse the record
    //
    iKey key;
    unsigned short I;
    int ISLP, ILV;
    double Eph, Eln, SW;

    sin >> I;
    sin >> key.first;
    sin >> key.second;
    sin >> ISLP;
    sin >> ILV;
    sin >> Eph;
    sin >> Eln;
    sin >> SW;

    // If header has been parsed, read the line data
    //
    if ( !sin.fail() && !sin.bad()) {
      if (SWlow.find(key) == SWlow.end()) SWlow[key] = SW;
      TBmapItr it = ions.find(key);
      if (it != ions.end()) {
	TBslp::iterator jt = it->second.find(ISLP);
	if (jt != it->second.end()) {
	  TBcfg::iterator kt = jt->second.find(ILV);
	  if (kt != jt->second.end()) {
	    kt->second->wght = SW;
	  }
	}
      }
    }
  }

}


void TopBase::printInfo()
{
  for (TBmapItr I=ions.begin(); I!=ions.end(); I++) {
    std::cout << std::string(60, '-') << std::endl
	      << " *** NZ = " << std::setw(3) << I->first.first
	      << ", NE = "    << std::setw(3) << I->first.second 
	      << std::endl;
    TBslp::iterator S    = I->second.begin();
    TBslp::iterator Send = I->second.end();
    for (; S != Send; S++) {
      TBcfg::iterator L    = S->second.begin();
      TBcfg::iterator Lend = S->second.end();
      
      for (; L != Lend; L++) {
	std::cout << "    "
		  << "  iSLP = " << std::setw(3)  << S->first 
		  << "  levl = " << std::setw(3)  << L->first
		  << "  g_n = "  << std::setw(10) << L->second->wght
		  << "  [" << std::setw(6) <<  L->second->NP << "] "
		  << "  [" << std::setw(16) << L->second->E.front()
		  << ", "  << std::setw(16) << L->second->E.back() << "]"
		  << std::endl;
      }
    }
  }
}


void TopBase::printLine(unsigned short NZ, unsigned short NE, 
			int I, int L, const std::string& file)
{
  TBmapItr ion = ions.find(iKey(NZ, NE));

  if (ion != ions.end()) {

    std::map<int, TBcfg>::iterator recl = ion->second.find(I);

    if (recl != ion->second.end()) {

      std::map<int, TBptr>::iterator line = recl->second.find(L);

      if (line != recl->second.end()) {
	std::ofstream out(file.c_str());
	TBptr l = line->second;
	if (out) {
	  for (int i=0; i<l->NP; i++) {
	    out << std::setw(16) << l->E[i]
		<< std::setw(16) << l->S[i]
		<< std::endl;
	  }

	} else throw FileOpen(file);

      } else throw NoLine(NZ, NE, I, L);
      
    } else throw NoSLP(NZ, NE, I);

  } else throw NoIon(NZ, NE);

}

double TopBase::sigmaFB(const iKey& key, double E)
{
  // Electron rest mass in eV
  //
  const double mec2 = 510.998896 * 1.0e3;

  // Return value
  //
  double cross = 0.0;

  if (key.second <= 1) return cross;

  // Neutral ion key
  //
  iKey low(key.first, key.second-1);

  // Iterate through all values for the ion
  //
  TBmapItr i = ions.find(low);

  if (i != ions.end()) {
    
    double mult0 = (key.second > key.first ? 1 : SWlow[key]);

    TBslp::iterator j = i->second.begin(), jend(i->second.end());

    for (; j != jend; j++) {

      TBcfg::iterator k = j->second.begin(), kend(j->second.end());

      for (; k != kend; k++) {

	TBptr l = k->second;

	double hnu  = E - l->Eph;
	double Erat = (hnu*hnu)/(2.0*mec2*E);
	double crs  = 0.0;

	if (E >= l->E.front() && E < l->E.back()) {

	  std::vector<double>::iterator lb = 
	    std::lower_bound(l->E.begin(), l->E.end(), E);
	  std::vector<double>::iterator ub = lb--;

	  size_t ii = lb - l->E.begin();
	  size_t jj = ub - l->E.begin();
	  
	  if (*ub > *lb)
	    crs = ( (E - *lb) * l->S[ii] + (*ub - E) * l->S[jj] ) /
	      (*ub - *lb) ;
	}

	double crossi = l->wght/mult0 * Erat * crs;

	cross += crossi;
      }
    }
  }

  // In megabarnes: (1.0e-11 m)^2
  // Convert to nm^2 (1.0e-9 m)^2

  return cross * 1.0e+04;
}
