#include <iostream>
#include <iomanip>

#include <ParamParseMPI.H>

ParamParseMPI::ParamParseMPI(istream* in, string Delim) : 
  ParamParse(in, Delim)
{
  const int lbufsize = 512;
  char lbuf[lbufsize];

  Stanza current;
  spair cpair;
  int signal;

  if (myid == 0) {

    list<Stanza>::iterator it;
    list<spair>::iterator it1;

    for (it=database.begin(); it!=database.end(); it++) {
      
      signal = 1;
      MPI_Bcast(&signal, 1, MPI_INT, 0, MPI_COMM_WORLD);

      strncpy(lbuf, it->name.c_str(), lbufsize);
#ifdef DEBUG
      cout << "Sending: " << lbuf << endl;
#endif
      MPI_Bcast(&lbuf, lbufsize, MPI_CHAR, 0, MPI_COMM_WORLD);

      for (it1=it->elist.begin(); it1!=it->elist.end(); it1++) {
	signal = 2;
	MPI_Bcast(&signal, 1, MPI_INT, 0, MPI_COMM_WORLD);
	
	strncpy(lbuf, it1->first.c_str(), lbufsize);
	MPI_Bcast(&lbuf, lbufsize, MPI_CHAR, 0, MPI_COMM_WORLD);
	
	strncpy(lbuf, it1->second.c_str(), lbufsize);
	MPI_Bcast(&lbuf, lbufsize, MPI_CHAR, 0, MPI_COMM_WORLD);

      }
    }

    signal = 0;
    MPI_Bcast(&signal, 1, MPI_INT, 0, MPI_COMM_WORLD);
    
  } else {

    MPI_Bcast(&signal, 1, MPI_INT, 0, MPI_COMM_WORLD);

    while (signal) 
      {
	if (signal == 1) {
	  if (!current.name.empty()) database.push_back(current);
	  MPI_Bcast(&lbuf, lbufsize, MPI_CHAR, 0, MPI_COMM_WORLD);
	  current.name = lbuf;
	  current.elist.erase(current.elist.begin(), current.elist.end());
	}

	if (signal == 2) {
	  MPI_Bcast(&lbuf, lbufsize, MPI_CHAR, 0, MPI_COMM_WORLD);
	  cpair.first = lbuf;
	
	  MPI_Bcast(&lbuf, lbufsize, MPI_CHAR, 0, MPI_COMM_WORLD);
	  cpair.second = lbuf;
	  current.elist.push_back(cpair);
	}

	MPI_Bcast(&signal, 1, MPI_INT, 0, MPI_COMM_WORLD);
      }

    database.push_back(current);
  }

}




