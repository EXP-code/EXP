#include <BarrierWrapper.H>

BarrierWrapper::BarrierWrapper(MPI_Comm communicator, bool label)
{
  comm = communicator;
  check_label = label;

  MPI_Comm_size(comm, &commsize);
  MPI_Comm_rank(comm, &localid);

  buffer = new char [cbufsz];
  bufferT = 0;
  if (localid==0) bufferT = new char [cbufsz*commsize];
}

BarrierWrapper::~BarrierWrapper()
{
  delete [] buffer;
  if (localid==0) delete [] bufferT;
}

void BarrierWrapper::operator()(const string& label)
{
  MPI_Barrier(comm);

  if (check_label) {
    strncpy(buffer, label.c_str(), cbufsz);
    MPI_Gather(buffer, cbufsz, MPI_CHAR, bufferT, cbufsz, MPI_CHAR, 0, comm); 

				// Compare strings
    if (localid==0) {
      char tmp[cbufsz];		// Working buffer
      bool firstime = true;
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
	cout << setfill('-') << setw(60) << '-' << endl << setfill(' ');
    }
  }

}
