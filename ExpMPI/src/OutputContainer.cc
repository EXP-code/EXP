#ifdef RCSID
static char rcsid[] = "$Id$";
#endif

#include "expand.h"

#include <OutputContainer.H>

#include <OutLog.H>
#include <OrbTrace.H>
#include <OutDiag.H>
#include <OutPS.H>

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
    cout << "." << n << flush;
    if (final) cout << "\n";
  }
}




