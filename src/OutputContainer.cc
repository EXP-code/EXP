
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
  spair data;
  
  parse->find_list("output");

  while ( parse->get_next(data) ) {

    if ( !data.first.compare("outlog") ) {
      out.push_back(new OutLog(data.second));
    }

    else if ( !data.first.compare("orbtrace") ) {
      out.push_back(new OrbTrace(data.second));
    }

    else if ( !data.first.compare("outdiag") ) {
      out.push_back(new OutDiag(data.second));
    }

    else if ( !data.first.compare("outps") ) {
      out.push_back(new OutPS(data.second));
    }

    else if ( !data.first.compare("outpsn") ) {
      out.push_back(new OutPSN(data.second));
    }
    
    else if ( !data.first.compare("outpsp") ) {
      out.push_back(new OutPSP(data.second));
    }
    
    else if ( !data.first.compare("outascii") ) {
      out.push_back(new OutAscii(data.second));
    }
    
    else if ( !data.first.compare("outchkpt") ) {
      out.push_back(new OutCHKPT(data.second));
    }

    else if ( !data.first.compare("outcoef") ) {
      out.push_back(new OutCoef(data.second));
    }

    else if ( !data.first.compare("outfrac") ) {
      out.push_back(new OutFrac(data.second));
    }

    else if ( !data.first.compare("outmulti") ) {
      out.push_back(new OutMulti(data.second));
    }

    else if ( !data.first.compare("outcalbr") ) {
      out.push_back(new OutCalbr(data.second));
    }

    else {
      string msg("I don't know about the output type: ");
      msg += data.first;
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




