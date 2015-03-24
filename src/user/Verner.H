
/** 
   Reads the Verner & Yakovlev (A&AS 109, 125, 1995) photoionization
   cross-section data
*/
class VernerData
{
public:

  struct VernerRec 
  {
    int pql, l;
    double eth, e0, sig0, ya, p, yw;

    void sync(int myid=0) 
    {
      MPI_Bcast(&pql,  1, MPI_INT,    myid, MPI_COMM_WORLD);
      MPI_Bcast(&l,    1, MPI_INT,    myid, MPI_COMM_WORLD);
      MPI_Bcast(&eth,  1, MPI_DOUBLE, myid, MPI_COMM_WORLD);
      MPI_Bcast(&e0,   1, MPI_DOUBLE, myid, MPI_COMM_WORLD);
      MPI_Bcast(&sig0, 1, MPI_DOUBLE, myid, MPI_COMM_WORLD);
      MPI_Bcast(&ya,   1, MPI_DOUBLE, myid, MPI_COMM_WORLD);
      MPI_Bcast(&p,    1, MPI_DOUBLE, myid, MPI_COMM_WORLD);
      MPI_Bcast(&yw,   1, MPI_DOUBLE, myid, MPI_COMM_WORLD);
    }
  };
  
  typedef boost::shared_ptr<VernerRec> vrPtr;

  std::map<std::pair<int, int>, vrPtr> data;

  VernerData() {}

  void initialize()
  {
    int nVern  = 465;
    int maxZ   = 30 + 1; 
    int maxNel = 31 + 1; 
    int nOK    = 0;

    if (myid==0) {

      char * val;
      if ( (val = getenv("CHIANTI_DATA")) == 0x0) {
	std::cout << "Could not find CHIANTI_DATA environment variable"
		  << " . . . exiting" << std::endl;
	nOK = 1;
      }

      if (nOK == 0) {

	std::string fileName(val);

	fileName.append("/continuum/verner_short.txt");
      
	std::string inLine;
	std::ifstream vdFile(fileName.c_str());
      
	if (vdFile.is_open()) {

	  while (vdFile.good()) {

	    std::vector<std::string> v;
	    std::getline(vdFile, inLine);
	    std::istringstream iss(inLine);
	    std::copy(std::istream_iterator<std::string>(iss), 
		      std::istream_iterator<std::string>(), 
		      std::back_inserter<std::vector<std::string> >(v));
	    
	    if (atoi(v[0].c_str()) == -1) break;
	  
	    int z   = atoi(v[0].c_str());
	    int nel = atoi(v[1].c_str());
	    int stg = z - nel + 1;
	    
	    std::pair<int, int> key(z, stg);

	    vrPtr dat(new VernerRec);

	    dat->pql  = atoi(v[2].c_str());
	    dat->l    = atoi(v[3].c_str());
	    dat->eth  = atof(v[4].c_str());
	    dat->e0   = atof(v[5].c_str());
	    dat->sig0 = atof(v[6].c_str());
	    dat->ya   = atof(v[7].c_str());
	    dat->p    = atof(v[8].c_str());
	    dat->yw   = atof(v[9].c_str());
	    
	    data[key] = dat;
	  }
	}
	vdFile.close();
      }
    }

    MPI_Bcast(&nOK, 1, MPI_UNSIGNED_CHAR, 0, MPI_COMM_WORLD);
    if (nOK) MPI_Abort(MPI_COMM_WORLD, 59);

    if (myid==0) {

      std::cout << "Verner data size=" << data.size() 
		<< ", expected: " << nVern << std::endl

      for (auto v : data) {
	std::pair<int, int> k = v.first;
	MPI_Bcast(&k.first,  1, MPI_INT, 0, MPI_COMM_WORLD);
	MPI_Bcast(&k.second, 1, MPI_INT, 0, MPI_COMM_WORLD);
	v.second->sync();
      }
  
      int done = 0;
      MPI_Bcast(&done,  1, MPI_INT, 0, MPI_COMM_WORLD);

    } else {
      std::pair<int, int> key;

      while (1) {
	MPI_Bcast(&key.first,  1, MPI_INT, 0, MPI_COMM_WORLD);
	if (key.first=0) break;
	MPI_Bcast(&key.second, 1, MPI_INT, 0, MPI_COMM_WORLD);

	vrPtr dat(new VernerRec);
	dat->sync();
	data[key] = dat;
      }
    }
  }

};

