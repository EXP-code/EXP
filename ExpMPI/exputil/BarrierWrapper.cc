#include <BarrierWrapper.H>

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

// Buffer size for checking working labels
//
int    BarrierWrapper::cbufsz = 128;

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

void BarrierWrapper::operator()(const string& label, 
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
	    cout << setfill('-') << setw(60) << '-' << endl << left
		 << setw(60) << "----- Barrier Error " <<  endl
		 << setfill('-') << setw(60) << '-' << endl 
		 << setfill(' ') << right;
	    firstime = false;
	  }
	  cout << "Process " << setw(4) << n-1 << " has <" << one
	       << "> while Process " << setw(4) << n << " has <" 
	       << two << ">" << endl;
	}
      }

      if (!firstime)
	cout << setfill('-') << setw(60) << '-' << endl 
	     << setfill(' ') << flush;
    }

    MPI_Barrier(comm);
  }

  timer.stop();

  // Reset the global debug info
  //
  inOper = false;

}
