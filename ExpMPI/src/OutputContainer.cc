#ifdef RCSID
static char rcsid[] = "$Id$";
#endif

#include "expand.h"

#include <OutputContainer.H>

#include <OutLog.H>
#include <OrbTrace.H>
#include <OutDiag.H>
#include <OutPS.H>
#include <OutPSN.H>
#include <OutAscii.H>
#include <OutCHKPT.H>
#include <OutCoef.H>
#include <OutFrac.H>
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

    else {
      string msg("I don't know about the output type: ");
      msg += data.first;
      if (myid) cerr << msg << endl;
      MPI_Abort(MPI_COMM_WORLD, 8);
      exit(0);
    }
  }
}

  
OutputContainer::~OutputContainer()
{
  list<Output*>::iterator it;
  for (it=out.begin(); it != out.end(); it++) delete *it;
}

void OutputContainer::Run(int n, bool final)
{
  list<Output*>::iterator it;
  for (it=out.begin(); it != out.end(); it++) (*it)->Run(n, final);
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




