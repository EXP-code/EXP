
#include <memory>
#include <BarrierWrapper.H>
#include <unistd.h>

#include <cstring>
#include <sstream>
#include <list>

//
// The next four are all for debugging backtraces
//
				// true while control is in the wrapper
bool   BarrierWrapper::inOper  = false;	
				// the label supplied to this process
string BarrierWrapper::lbOper;
				// the source file name with the call
string BarrierWrapper::flOper;
				// the source file line number
int    BarrierWrapper::lnOper;

//
// BarrierWrapper parameters
//
				// Set to true for the original version
bool   BarrierWrapper::light         = false;
				// Debugging info
bool   BarrierWrapper::verbose       = true;
				// Extra checking
bool   BarrierWrapper::extra_verbose = false;
				// Debugging info
bool   BarrierWrapper::debugging     = false;
				// Expiration timers
int    BWData::dt1                   = 5;
int    BWData::dt2                   = 20;
				// Loop wait in microseconds
int    BarrierWrapper::loop_delay    = 100;
				// Buffer size for checking working
				// labels
int    BarrierWrapper::cbufsz        = 128;

BarrierWrapper::BarrierWrapper(MPI_Comm communicator, bool label)
{
  comm        = communicator;
  check_label = label;

  MPI_Comm_size(comm, &commsize);
  MPI_Comm_rank(comm, &localid);

  buffer      = new char [cbufsz];
  bufferT     = 0;
  if (localid==0) 
    bufferT   = new char [cbufsz*commsize];
  onoff       = true;
  nrecv       = 0;
  queued      = 0;
}

BarrierWrapper::BarrierWrapper(const BarrierWrapper &p)
{
  comm        = p.comm;
  check_label = p.check_label;
  commsize    = p.commsize;
  localid     = p.localid;
  onoff       = p.onoff;

  buffer      = new char [cbufsz];
  bufferT     = 0;
  if (localid==0) 
    bufferT   = new char [cbufsz*commsize];
  nrecv       = 0;
  queued      = 0;
}

BarrierWrapper::~BarrierWrapper()
{
  delete [] buffer;
  delete [] bufferT;
}

void BarrierWrapper::light_operator(const string& label, 
				    const char* file, const int line)
{
  if (!onoff) return;

  // Set the global debug info
  //
  inOper = true;
  lbOper = label;
  flOper = string(file);
  lnOper = line;

  timer.start();

  if (check_label) {
				// Copy the label to the send buffer
    strncpy(buffer, label.c_str(), cbufsz);

				// Gather to the root node (0)
    MPI_Gather(buffer, cbufsz, MPI_CHAR, bufferT, cbufsz, MPI_CHAR, 0, comm); 

				// Compare adjacent strings in the list
    if (localid==0) {
      char tmp[cbufsz];		// Working buffer
      bool firstime   = true;
      string one, two = strncpy(tmp, &bufferT[0], cbufsz);
      for (int n=1; n<commsize; n++) {
	one = two;
	two = strncpy(tmp, &bufferT[cbufsz*n], cbufsz);
	if (one.compare(two) != 0) {
	  if (firstime) {
	    cout << std::string(60,'-') << std::endl << left
		 << std::setw(60) << "----- Barrier Error " <<  endl
		 << std::string(60,'-') << std::endl << right;
	    firstime = false;
	  }
	  cout << "Process " << setw(4) << n-1 << " has <" << one
	       << "> while Process " << setw(4) << n << " has <" 
	       << two << ">" << endl;
	}
      }

      if (!firstime)
	cout << std::string(60,'-') << std::endl;
    }

    MPI_Barrier(comm);
  }

  timer.stop();

  // Reset the global debug info
  //
  inOper = false;

}


// This loop defeats the purpose of this entire class and should ONLY
// be used for debugging
//
void BarrierWrapper::syncTest(const std::string& mesg, const std::string& label)
{
  for (int n=0; n<commsize; n++) {
    if (n == localid) {
      std::cout << std::left << std::setw(30) << mesg
		<< " {" << std::setw(5) << localid
		<< " name=" << label
		<< "} " << std::endl;
    }
    MPI_Barrier(comm);
  }
}


void BarrierWrapper::updateMap(InfoPtr p)
{
  // Enter our own data into the pending map
  //
  std::map<std::string, BWPtr>::iterator ipend = pending.find(p->s);
  if (ipend == pending.end()) { 
				// No, ADD it
				// 
    pending[p->s] = BWPtr(new BWData(p, commsize));

    nrecv  += commsize - 1;	// Expect a receive from each of the
    queued += commsize - 1;	// other nodes

    if (debugging) {
      std::cout << "Node " << std::setw(4) << localid << ", in updateMap, has "
		<< nrecv << " queued receives and buffer count of "
		<< queued << " at <" << p->s << ">"
		<< std::endl;
    }

  } else {
				// Update local info and increment
				// counter
    ipend->second->Add(p);

  }

  if (nrecv) {			// More left to go . . .queue another
				// receive
    CharPtr cptr(new char [Info::bufsz]);
    ReqPtr  rptr(new MPI_Request);

    req.push_back(RecvElem(cptr, rptr));
    MPI_Irecv(cptr.get(), Info::bufsz, MPI_CHAR, MPI_ANY_SOURCE, Info::tag, 
	      comm, rptr.get());
    nrecv--;
  }

}


void BarrierWrapper::heavy_operator(const string& label, 
				    const char* file, const int line)
{
  if (!onoff || commsize==1) return;

  //-----------------------------------------------------------------
  // Set the global debug info; reported by the EXP global exception
  // handler
  //-----------------------------------------------------------------
  inOper = true;
  lbOper = label;
  flOper = string(file);
  lnOper = line;

  timer.start();

  //-------------------
  // Create barrier id
  //-------------------
  std::ostringstream sid;
  sid << label << " [" << file << ":" << line << "]";
  
  //---------------------
  // Initial diagnostics
  //---------------------
  if (pending.size() > 1) {
    std::cout << "#" << localid << " entered BW with " << pending.size()
	      << " unresolved tickets" << std::endl;
    for(auto i : pending)
      std::cout << " ** #" << localid << " {" << std::left << i.first << "} "
		<< std::endl;
  }

  //---------------------
  // For deeper checking
  //---------------------
  std::map<int, size_t> multiplicity;

  //------------------------
  // Loop control variables
  //------------------------
  time_t entry = time(0);
  bool notdone = true;

  std::string sfinal;
  MPI_Status status;
  int good = 0;
  
  //----------------------------------------------------------------
  // Create and send info buffer in a non-blocking send to all
  // processes
  // ----------------------------------------------------------------

  InfoPtr p = InfoPtr(new Info(localid, entry, sid.str()));

  p->sendInfo(localid, comm, commsize);

  //----------------------------------------------------------------
  // Enter our ticket into the registry
  //----------------------------------------------------------------

  updateMap(p);

  //----------------------------------------------------------------
  // Enter the loop: read and process info buffers until one tag is
  // full
  //----------------------------------------------------------------

  unsigned cnt_recv = 0;

  // Are all posted messages processed?
  //
  while (notdone) {
    
    // Are more buffers expected?
    //
    if (queued) {

      if (extra_verbose) {
	static time_t next = 0, cur = time(0);
	if (next==0) next = entry + 10;
	if (cur>next) {
	  int total = 0;
	  for (auto i : pending) total += i.second->count;
	  std::cout << "EXTRA: dT=" 
		    << std::setw(4) << std::left << cur - entry
		    << ", Node " << std::setw(4) << localid << " scanning "
		    << req.size() << " pending receives, expecting "  
		    << commsize - total << std::endl;
	  next += 10;		// wait 10 more seconds
	}
      }

      // Check for buffers pending reception
      //
      MPI_Test(req.back().second.get(), &good, &status);

    } else {			// THIS IS A SANITY CHECK
				// No buffers expected . . . why are
				// we still in this loop?
      if (debugging) {
	std::cout << "No request for Node " << std::setw(3) << localid 
		  << ", done is " << (notdone ? 0 : 1) << ", received "
		  << cnt_recv << ", ";
	if (pending.size()) {
	  std::cout << std::setw(4) << pending.size() << " barrier(s)"
		    << " named ";
	  for (auto i : pending)
	    std::cout << "<" << i.first << ">, #=" << i.second->count;
	}
      } else {
	std::cout << " no pending barrriers";
      }
      std::cout << std::endl;
    }
      
				// Did we satisfy a request?
    if (good) {
    
      CharPtr c = req.back().first;
      req.pop_back();		// Pop the completed request
      InfoPtr p = InfoPtr(new Info(c));

      updateMap(p);		// Will generate another receive if
				// nrecv>0
    
      queued--;			// Decrement the expected buffers
				// count
      
      //---------------------------------------
      // Check all pending tags for completion
      //---------------------------------------
      for (auto it=pending.begin(); it!=pending.end(); it++) {

	if (it->second->count == commsize) {
	    
	  if (debugging) {
	    std::cout << "Node " << std::setw(4) << localid 
		      << " has count " << it->second->count 
		      << " with " << queued << " queued for <" 
		      << it->second->info.back()->s << ">, size=" 
		      << it->second->info.size() << std::endl;
	  }

	  // All nodes have now hit the barrier
	  sfinal = it->first;
	  pending.erase(it);
	  notdone = false;
	  break;
	}
      }
	
      //----------------------------------------------------------
      // Check for bad synchronization: more than one active tag!
      //----------------------------------------------------------

      if (pending.size() > 1 && extra_verbose) {
				// Only want to report multiple tags
				// on the first encounter
	std::map<int, size_t>::iterator im = multiplicity.find(localid);
	bool report = true;
	if (im != multiplicity.end()) {
	  if (pending.size() == im->second) report = false;
	}
	multiplicity[localid] = pending.size();
	
	if (report) {		// Only report a CHANGE in multiplicity
	  
	  std::cout << "#" << localid << ": pending list has " 
		    << pending.size() << " entries: ";
	  for (auto i : pending)
	    std::cout << " {" << std::left << i.first << "} ";
	  std::cout<< std::endl;
	}
      }
    }
    
    if (notdone) {
      //---------------
      // Update alarms
      //---------------
      time_t curtime = time(0);
      for (auto it=pending.begin(); it != pending.end(); it++) {
	if (it->second->Owner() == localid && 
	    curtime > it->second->expire) {
	  listReport("Expire", it);
	  it->second->expire += BWData::dt2 * it->second->nexpiry++;
	}
      }
  
      //--------------------
      // Print pending list
      //--------------------

      if (extra_verbose) {
	
	static time_t next = 0;
	
	if (next==0) next = entry + 10;

	if (curtime>next) {
	  
	  for (auto i : pending) {

	    if (curtime > i.second->first + BWData::dt2*4) {
	      std::vector<int> ret = getMissing(i.second);

	      std::cout << "EXTRA: dT=" << std::setw(4)  << std::left
			<<  curtime - entry << ", Node " 
			<< std::setw(4) << localid << ", barrier=" << i.first
			<< ", " << req.size() << " MPI request(s) remaining, "
			<< ret.size() << " wait(s)";
	      
	      if (ret.size()) {
		std::cout << ", nodes:";
		for (std::vector<int>::iterator 
		       n=ret.begin(); n!=ret.end(); n++)
		  std::cout << " " << *n;
	      }
	      std::cout << std::endl;

	    } // time check

	  } // pending barrier loop
	  
	  next += 10;

	} // 10 more seconds

      }	// extra_verbose
	
      //----------------------------------
      // Wait loop_delay (default:100) us
      //----------------------------------
      usleep(loop_delay);
    }
  }

  //---------------------
  // Final diagnostics
  //---------------------
  if (pending.size() > 1) {
    std::cout << "#" << localid << " leaving BW with " << pending.size()
	      << " unresolved tickets" << std::endl;;
    std::map<std::string, BWPtr>::iterator it;
    for (auto i : pending)
      std::cout << " ** #" << localid << " {" << std::left << i.first << "} "
		<< std::endl;
  }


  finalReport(sfinal);

  timer.stop();

  //--------------------------
  // Here is the REAL barrier
  //--------------------------
  MPI_Barrier(MPI_COMM_WORLD);

  //-----------------------------
  // Reset the global debug info
  //-----------------------------
  inOper = false;
}

void BarrierWrapper::listReport(const char* title, 
				std::map<std::string, BWPtr>::iterator it)
{
  std::cout << title << " [#" << localid << "," << it->second->nexpiry
	    << "]: " << it->first << " ** "
	    << time(0) - it->second->first << " secs, " 
	    << it->second->count << "/" << commsize << " waiting, ";
  for (int i=0; i<commsize; i++) std::cout << it->second->nd[i];
  std::cout << std::endl;
}

void Info::pack()
{
  if (blob == 0) {
    charsz = bufsz - (sizeof(int) + sizeof(time_t) + sizeof(unsigned));
    blob = CharPtr(new char [bufsz]);
  }
  
  unsigned cursor = 0;

  *reinterpret_cast<int*>(blob.get()+cursor) = own;
  cursor += sizeof(int);

  *reinterpret_cast<time_t*>(blob.get()+cursor) = ctm;
  cursor += sizeof(time_t);

  *reinterpret_cast<unsigned*>(blob.get()+cursor) = siz;
  cursor += sizeof(unsigned);
  
  strncpy(blob.get()+cursor, s.c_str(), bufsz-cursor);
}

Info::Info(CharPtr p)
{
  if (blob == 0) {
    charsz = bufsz - (2*sizeof(int) + sizeof(time_t) + sizeof(unsigned));
    blob = CharPtr(new char [bufsz]);
  }

  memcpy(blob.get(), p.get(), bufsz);
  
  unsigned cursor = 0;

  own = *reinterpret_cast<int*>(blob.get()+cursor);
  cursor += sizeof(int);

  ctm = *reinterpret_cast<time_t*>(blob.get()+cursor);
  cursor += sizeof(time_t);

  siz = *reinterpret_cast<unsigned*>(blob.get()+cursor);
  cursor += sizeof(unsigned);
  
  c = CharPtr(new char [siz]);

  memcpy(c.get(), blob.get()+cursor, std::min<size_t>(siz, bufsz-cursor));

  s = std::string(c.get());
}


void Info::sendInfo(int source, MPI_Comm comm, int commsize)
{
  pack();

  for (int i=0; i<commsize; i++) {
    if (i!=source) {
      req.push_back(ReqPtr(new MPI_Request));
      MPI_Isend(blob.get(), bufsz, MPI_CHAR, i, tag, comm, req.back().get());
    }
  }
}

void BarrierWrapper::finalReport(std::string& s)
{
  std::map<std::string, std::vector<bool> >::iterator is;
  std::map<std::string, BWPtr>::iterator it;

  if (localid == 0) {
				// Do the current node
    std::map<string, std::vector<bool> > table;
    for (auto i : pending) {
      table[i.first] = vector<bool>(commsize, false);
      table[i.first][localid] = true;
    }
				// Loop to receive info from all nodes
    for (int i=1; i<commsize; i++) {
      int num, siz;
      MPI_Recv(&num, 1, MPI_INT, i, 113549, comm, MPI_STATUS_IGNORE);
      for (int j=0; j<num; j++) {
	MPI_Recv(&siz, 1, MPI_INT, i, 113550, comm, MPI_STATUS_IGNORE);
	CharPtr c = CharPtr(new char [siz]);
	MPI_Recv(c.get(), siz, MPI_CHAR, i, 113551, comm, MPI_STATUS_IGNORE);
	std::string s(c.get());
	is = table.find(s);
	if (is != table.end()) {
	  table[s][i] = true;
	} else {
	  table[s] = vector<bool>(commsize, false);
	  table[s][i] = true;
	}
      }
    }

    if (table.size() == 0) {
      if (verbose) {
	std::cout << "BarrierWrapper: state is good for " << s 
		  << ", req size=" << req.size() << std::endl;
      }
    } else {
      std::cout << "BarrierWrapper: suspect sync problem for " << s 
		<< std::endl;
    }

    for (auto t : table)
      {
	std::cout << std::setw(32) << t.first << " : ";
	for (int i=0; i<commsize; i++) std::cout << t.second[i];
	std::cout << std::endl;
      }

  } else {
    
    int num = pending.size();
    MPI_Send(&num, 1, MPI_INT, 0, 113549, comm);
    for (auto i : pending) {
      int siz = i.first.size() + 1;
      MPI_Send(&siz, 1, MPI_INT, 0, 113550, comm);
      CharPtr c = CharPtr(new char [siz]);
      strncpy(c.get(), i.first.c_str(), siz);
      MPI_Send(c.get(), siz, MPI_CHAR, 0, 113551, comm);
    }
  }

}


std::vector<int> BarrierWrapper::getMissing(BWPtr p) 
{
  std::vector<int> ret;
  for (int n=0; n<commsize; n++) {
    if (!p->nd[n]) ret.push_back(n);
  }
  return ret;
}

Info::Info(int own, time_t ctm, const std::string& s)
{
  this->own = own;
  this->ctm = ctm;
  this->s   = s;
  this->siz = s.size() + 1;
  this->c   = CharPtr(new char [siz]);
  strncpy(c.get(), s.c_str(), siz);
}



BWData::BWData(InfoPtr& p, int commsize) 
{
  first       = p->ctm;
  expire      = first + dt1;
  nexpiry     = 1;
  count       = 1;
  nd          = std::vector<bool>(commsize, false);
  nd[p->own]  = true;
  
  info.push_back(p);
}

void BWData::Add(InfoPtr& p) 
{ 
  info.push_back(p); 
  count++;
  nd[p->own] = true;
}

int BWData::Owner()
{
  int i, sz = nd.size();
  for (i=0; i<sz; i++) { if (nd[i]) break; }
  return i;
}

