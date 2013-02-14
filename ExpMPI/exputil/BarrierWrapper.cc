
#include <boost/shared_ptr.hpp>
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
bool   BarrierWrapper::extra_verbose = false;
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
  timer.Microseconds();
  onoff       = true;
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
  timer.Microseconds();
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

  //----------------------------------------------
  // Temporary character buffers
  // [Use shared pointers as a garbage collector]
  //----------------------------------------------
  std::list<InfoPtr> info;

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
	      << " unresolved tickets: ";
    std::map<std::string, BWData>::iterator it;
    for(it=pending.begin(); it!=pending.end(); it++) 
      std::cout << " {" << std::left << it->first << "} ";
    std::cout<< std::endl;
  }

  //----------------
  // Loop variables
  //----------------
  unsigned siz = sid.str().size()+1;
  std::map<int, size_t> multiplicity;
  time_t entry = time(0);
  bool notdone = true;
  std::string sfinal;
  MPI_Status status;
  MPI_Request req;
  int good;

  //----------------------------------------------------------------
  // Enter and send the caller; sends info buffer in a non-blocking
  // send to all processes
  //----------------------------------------------------------------
  InfoPtr p = InfoPtr(new Info(localid, entry, sid.str()));
  pending[sid.str()] = BWData(p->own, commsize, p->ctm);
  info.push_back(p);
  p->sendInfo(localid, comm, commsize);
  
  //-------------------------------------------------------------
  // Begin receiving commands; initiate receipt of incoming info
  // buffer
  //-------------------------------------------------------------
  CharPtr buf(new char [Info::bufsz]);
  MPI_Irecv(buf.get(), Info::bufsz, MPI_CHAR, MPI_ANY_SOURCE, Info::tag, 
	    comm, &req);

  //----------------------------------------------------------------
  // Enter the loop: read and process info buffers until one tag is
  // full
  //----------------------------------------------------------------
  while (notdone) {
				// Check for command
    MPI_Test(&req, &good, &status);
				// For debugging this class!
    if (0) {
      std::cout << "#" << std::setw(4) << std::left << localid;
      for (std::map<std::string, BWData>::iterator 
	     it = pending.begin(); it != pending.end(); it++) 
	{
	  std::cout << " {" << std::left << it->first 
		    << ", time=" << it->second.owner.begin()->first 
		    << " owner=" << it->second.owner.begin()->second
		    << "} ";
	}
      std::cout<< std::endl;
    }

    if (good) {
      
      InfoPtr p = InfoPtr(new Info(buf));
      std::map<std::string, BWData>::iterator ipend = pending.find(p->s);
	
      //----------------------------
      // Do we know about this one?
      //----------------------------

      if (ipend == pending.end()) { 
				// No, ADD it
	pending[p->s] = BWData(p->own, commsize, p->ctm);
	
      } else {
				// Update local info
	ipend->second.Add(p->own, p->ctm);
	ipend->second.count++;
	ipend->second.nd[p->own] = true;
      }

      //---------------------------------------
      // Check all pending tags for completion
      //---------------------------------------
      for (std::map<std::string, BWData>::iterator 
	     it = pending.begin(); it != pending.end(); it++) 
	{
	  if (it->second.count == commsize) {
	    sfinal = it->first;
	    pending.erase(it);
	    notdone = false;
	    break;
	  }
	}
      
      //---------------------------------- 
      // Start receiving a new Info buffer
      //----------------------------------
      if (notdone)
	MPI_Irecv(buf.get(), Info::bufsz, MPI_CHAR, MPI_ANY_SOURCE, Info::tag, 
		  comm, &req);

      //----------------------------------------------------------
      // Check for synchronize problems: more than one active tag
      //----------------------------------------------------------
      if (pending.size() > 1 && extra_verbose) {
	std::map<int, size_t>::iterator im = multiplicity.find(localid);
	bool report = true;
	if (im != multiplicity.end()) {
	  if (pending.size() == im->second) report = false;
	}
	multiplicity[localid] = pending.size();

	if (report) {		// Report a change in multiplicity
	  std::cout << "#" << localid << ": pending list has " 
		    << pending.size() << " entries: ";
	  for(std::map<std::string, BWData>::iterator it = pending.begin();
	      it != pending.end(); it++) 
	    std::cout << " {" << std::left << it->first << "} ";
	  std::cout<< std::endl;
	}
      }
    }
    
    if (notdone) {
      //---------------
      // Update alarms
      //---------------
      time_t curtime = time(0);
      for (std::map<std::string, BWData>::iterator 
	     it = pending.begin(); it != pending.end(); it++) {
	if (it->second.Owner() == localid && 
	    curtime > it->second.expire) {
	  listReport("Expire", it);
	  it->second.expire += BWData::dt2;
	}
      }
  
      //----------------------------------
      // Wait loop_delay (default:100) us
      //----------------------------------
      usleep(loop_delay);
    }
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
				std::map<std::string, BWData>::iterator it)
{
  std::cout << title << " [#" << localid << "]: " << it->first << " ** "
	    << time(0) - it->second.CTime() << " secs, " 
	    << it->second.count << " waiting, ";
  for (int i=0; i<commsize; i++) std::cout << it->second.nd[i];
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

void Info::unpackInfo(CharPtr p)
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

  MPI_Request req;
  for (int i=0; i<commsize; i++) {
    if (i!=source) {
      MPI_Isend(blob.get(), bufsz, MPI_CHAR, i, tag, comm, &req);
      MPI_Request_free(&req);
    }
  }
}

void BarrierWrapper::finalReport(std::string& s)
{
  std::map<std::string, std::vector<bool> >::iterator is;
  std::map<std::string, BWData>::iterator it;

  if (localid == 0) {
				// Do the current node
    std::map<string, std::vector<bool> > table;
    for (it = pending.begin(); it != pending.end(); it++) {
      table[it->first] = vector<bool>(commsize, false);
      table[it->first][localid] = true;
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
		  << std::endl;
      }
    } else {
      std::cout << "BarrierWrapper: suspect sync problem for " << s 
		<< std::endl;
    }

    for (is=table.begin(); is!=table.end(); is++) 
      {
	std::cout << std::setw(32) << is->first << " : ";
	for (int i=0; i<commsize; i++) std::cout << is->second[i];
	std::cout << std::endl;
      }

  } else {
    
    int num = pending.size();
    MPI_Send(&num, 1, MPI_INT, 0, 113549, comm);
    for (it = pending.begin(); it != pending.end(); it++) {
      int siz = it->first.size() + 1;
      MPI_Send(&siz, 1, MPI_INT, 0, 113550, comm);
      CharPtr c = CharPtr(new char [siz]);
      strncpy(c.get(), it->first.c_str(), siz);
      MPI_Send(c.get(), siz, MPI_CHAR, 0, 113551, comm);
    }
  }

}
