
#include "expand.H"

#include "OutputContainer.H"

#include "OutLog.H"
#include "OrbTrace.H"
#include "OutDiag.H"
#include "OutPS.H"
#include "OutPSN.H"
#include "OutPSP.H"
#include "OutPSQ.H"
#include "OutPSR.H"
#include "OutVel.H"
#include "OutHDF5.H"
#include "OutAscii.H"
#include "OutCHKPT.H"
#include "OutCHKPTQ.H"
#include "OutCoef.H"
#include "OutFrac.H"
#include "OutCalbr.H"
#include "OutMulti.H"
#include "OutSample.H"

OutputContainer::OutputContainer()
{
  // Mark time ahead of current time on restart
  //
  if (restart) last = tnow + 0.6*dtime;
  //
  // and behind current time on restart
  //
  else         last = tnow - 0.6*dtime;
}

void OutputContainer::initialize(void)
{
  YAML::Node outs = parse["Output"];

  if (outs.IsSequence()) {

    int nout = 0;

    while (outs[nout]) {
      std::string name = outs[nout]["id"].as<std::string>();
      const YAML::Node& node = outs[nout]["parameters"];
    
      if ( !name.compare("outlog") ) {
	out.push_back(new OutLog(node));
      }
      
      else if ( !name.compare("orbtrace") ) {
	out.push_back(new OrbTrace(node));
      }
      
      else if ( !name.compare("outdiag") ) {
	out.push_back(new OutDiag(node));
      }

      else if ( !name.compare("outps") ) {
	out.push_back(new OutPS(node));
      }

      else if ( !name.compare("outpsn") ) {
	out.push_back(new OutPSN(node));
      }
    
      else if ( !name.compare("outpsp") ) {
	out.push_back(new OutPSP(node));
      }
    
      else if ( !name.compare("outpsq") ) {
	out.push_back(new OutPSQ (node));
      }
    
      else if ( !name.compare("outhdf5") ) {
	out.push_back(new OutHDF5 (node));
      }
    
      else if ( !name.compare("outpsr") ) {
	out.push_back(new OutPSR (node));
      }
    
      else if ( !name.compare("outvel") ) {
	out.push_back(new OutVel (node));
      }
    
      else if ( !name.compare("outascii") ) {
	out.push_back(new OutAscii(node));
      }
    
      else if ( !name.compare("outchkpt") ) {
	out.push_back(new OutCHKPT(node));
      }

      else if ( !name.compare("outchkptq") ) {
	out.push_back(new OutCHKPTQ(node));
      }

      else if ( !name.compare("outcoef") ) {
	out.push_back(new OutCoef(node));
      }

      else if ( !name.compare("outfrac") ) {
	out.push_back(new OutFrac(node));
      }

      else if ( !name.compare("outmulti") ) {
	out.push_back(new OutMulti(node));
      }

      else if ( !name.compare("outcalbr") ) {
	out.push_back(new OutCalbr(node));
      }
      
      else if ( !name.compare("outsamp") ) {
	out.push_back(new OutSample(node));
      }

      else {
	string msg("I don't know about the output type: ");
	msg += name;
	throw GenericError(msg, __FILE__, __LINE__, 1025, false);
      }
      nout++;

      // Check YAML configuration
      auto unmatched = out.back()->unmatched();
      if (unmatched.size()) {
	throw YamlConfigError("OutputContainer", name, unmatched, __FILE__, __LINE__, 1026, false);
      }

    }
  } else {
    if (myid==0)
      std::cout << std::string(72, '-') << std::endl
		<< "No output entries" << std::endl
		<< std::string(72, '-') << std::endl;
  }

}

  
OutputContainer::~OutputContainer()
{
  // Delete all Output instances
  //
  for (auto it : out) delete it;
}

void OutputContainer::Run(int nstep, int mstep, bool final)
{
  // Don't rerun a step unless EXP is quitting . . . but allow for
  // multisteps to be run
  //
  if (not stop_signal and fabs(tnow - last) < 0.5*dtime/Mstep) return;

#ifdef HAVE_LIBCUDA
  // List of components for cuda fetching
  //
  if (use_cuda) {
    for (auto c : comp->components) comp->fetched[c] = false;

    // Wait check that all previous threads are finished
    //
    for (auto v : cproc) v->join();
  }
#endif
  
  cproc.clear();		// Delete the threads

  // Loop through all instances
  //
  for (auto it : out) it->Run(nstep, mstep, final);
  
  // Root node output
  //
  if (myid==0 and mstep==0) {
#ifdef DEBUG
    cout << setw(60) << setfill('=') << "=" << endl
	 << "====== Step " << nstep << "/" << mstep << endl
	 << setw(60) << setfill('=') << "=" << endl
	 << setfill(' ');
#else
    cout << "." << nstep << flush;
#endif
    if (final) cout << "\n";
  }

  // Mark: step ran at this time
  //
  last = tnow;

  // Wait check that all previous threads are finished
  //
  if (not use_cuda) {
    for (auto v : cproc) v->join();
  }

}
