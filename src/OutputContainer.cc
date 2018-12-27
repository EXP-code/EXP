
#include "expand.h"

#include <OutputContainer.H>

#include <OutLog.H>
#include <OrbTrace.H>
#include <OutDiag.H>
#include <OutPS.H>
#include <OutPSN.H>
#include <OutPSP.H>
#include <OutAscii.H>
#include <OutCHKPT.H>
#include <OutCoef.H>
#include <OutFrac.H>
#include <OutCalbr.H>
#include <OutMulti.H>

OutputContainer::OutputContainer() {}

void OutputContainer::initialize(void)
{
  YAML::Node outs = parse["output"];

  for (YAML::const_iterator it=outs.begin(); it!=outs.end(); ++it) {

    std::string name = it->first.as<std::string>();

    if ( !name.compare("outlog") ) {
      out.push_back(new OutLog(it->second));
    }

    else if ( !name.compare("orbtrace") ) {
      out.push_back(new OrbTrace(it->second));
    }

    else if ( !name.compare("outdiag") ) {
      out.push_back(new OutDiag(it->second));
    }

    else if ( !name.compare("outps") ) {
      out.push_back(new OutPS(it->second));
    }

    else if ( !name.compare("outpsn") ) {
      out.push_back(new OutPSN(it->second));
    }
    
    else if ( !name.compare("outpsp") ) {
      out.push_back(new OutPSP(it->second));
    }
    
    else if ( !name.compare("outascii") ) {
      out.push_back(new OutAscii(it->second));
    }
    
    else if ( !name.compare("outchkpt") ) {
      out.push_back(new OutCHKPT(it->second));
    }

    else if ( !name.compare("outcoef") ) {
      out.push_back(new OutCoef(it->second));
    }

    else if ( !name.compare("outfrac") ) {
      out.push_back(new OutFrac(it->second));
    }

    else if ( !name.compare("outmulti") ) {
      out.push_back(new OutMulti(it->second));
    }

    else if ( !name.compare("outcalbr") ) {
      out.push_back(new OutCalbr(it->second));
    }

    else {
      string msg("I don't know about the output type: ");
      msg += name;
      throw GenericError(msg, __FILE__, __LINE__);
    }
  }
}

  
OutputContainer::~OutputContainer()
{
  for (auto it : out) delete it;
}

void OutputContainer::Run(int n, bool final)
{
  for (auto it : out) it->Run(n, final);
  if (myid==0) {
#ifdef DEBUG
    cout << setw(60) << setfill('=') << "=" << endl
	 << "====== Step " << n << endl
	 << setw(60) << setfill('=') << "=" << endl
	 << setfill(' ');
#else
    cout << "." << n << flush;
#endif
    if (final) cout << "\n";
  }
}




