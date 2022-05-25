#include <algorithm>
#include <iostream>
#include <iomanip>
#include <fstream>
#include <sstream>
#include <numeric>
#include <string>
#include <vector>


#include <yaml-cpp/yaml.h>	// YAML support

#include <mpi.h>		// MPI support

#include <H5Cpp.h>		// HDF5 C++ support

#include <ParticleReader.H>
#include <gadget.H>
#include <Sutils.H>		// For string trimming

namespace PR {

  std::vector<std::string> GadgetNative::Ptypes
  {"Gas", "Halo", "Disk", "Bulge", "Stars", "Bndry"};
  
  std::unordered_map<std::string, int> GadgetNative::findP
  { {"Gas", 0}, {"Halo", 1}, {"Disk", 2}, {"Bulge", 3}, {"Stars", 4}, {"Bndry", 5}};
  
  
  GadgetNative::GadgetNative(const std::vector<std::string>& files, bool verbose)
  {
    _files   = files;		// Copy file list (bunch)
    _verbose = verbose;
    
    ptype = 1;			// Default is halo particles

    getNumbers();		// Get the number of particles in all
				// files

    curfile = _files.begin();	// Set file to first one

    if (not nextFile()) {	// Try opening
      std::cerr << "GadgetNative: no files found" << std::endl;
    }
  }


  void GadgetNative::getNumbers()
  {
    std::set<std::string> pfound;
    std::fill(nptot, nptot+6, 0);

    for (auto f : _files) {

      // Attempt to open file
      //
      std::ifstream file(f, std::ios::binary | std::ios::in);
      if (!file.is_open())
      {
	std::cerr << "Error opening file: " << f << std::endl;
	int flag;
	MPI_Initialized(&flag);
	if (flag) MPI_Finalize();
	exit(1);
      }
      
      // read in file data
      //
      if (myid==0 and _verbose)
	std::cout << "GadgetNative: reading " << f << " header...";
    
      file.seekg(sizeof(int), std::ios::cur); // block count
    
      file.read((char*)&header, sizeof(gadget_header)); 
    
      if (myid==0 and _verbose) std::cout << "done" << std::endl;
      
      for (int n=0; n<6; n++) {
	nptot[n] += npart[n];
	if (header.npart[n] > 0) pfound.insert(Ptypes[n]);
      }
    }
    
    Pfound.clear();
    for (auto p : pfound) Pfound.push_back(p);
  }
   

  bool GadgetNative::nextFile()
  {
    if (curfile==_files.end()) return false;
    read_and_load();
    curfile++;
    return true;
  }

  void GadgetNative::read_and_load()
  {
    // attempt to open file
    //
    std::ifstream file(*curfile, std::ios::binary | std::ios::in);
    if (!file.is_open()) {
      std::ostringstream ost;
      ost << "Error opening file: " << *curfile;
      throw std::runtime_error(ost.str());
    }
    
    if (myid==0 and _verbose)
      std::cout << "GadgetNative: opened <" << *curfile << ">" << std::endl;
    
    // read in file data
    //
    file.seekg(sizeof(int), std::ios::cur); // block count
    
    file.read((char*)&header, sizeof(gadget_header)); 
    
    file.seekg(sizeof(int), std::ios::cur); // block count
    
    // Total number of particles in this stanza.  Note: this needs to be
    // generalized for reading multiple Gadget native files per snap
    //
    totalCount = header.npart[ptype];
    
    particles.clear();		// Should be empty, but enforce that
    
    Particle P;			// Temporary for packing array
    
    // Read positions
    //
    file.seekg(sizeof(int), std::ios::cur); // block count
    
    float temp[3];		// Pos/vel temporary
    
    for (int k=0; k<6; k++) {
      if ( k == ptype ) {
	for (int n=0; n<header.npart[k]; n++) {
	  file.read((char*)temp, 3*sizeof(float));
	  if (n % numprocs == myid) {
	    P.pos[0] = temp[0];
	    P.pos[1] = temp[1];
	    P.pos[2] = temp[2];
	    P.level  = 0;		// Assign level 0 to all particles
	    particles.push_back(P);
	  }
	}
      }
      else {
	file.seekg(header.npart[k]*3*sizeof(float), std::ios::cur);
      }
    }
    
    file.seekg(sizeof(int), std::ios::cur); // block count
    file.seekg(sizeof(int), std::ios::cur); // block count
    
    // Read velocities
    //
    for (int k=0; k<6; k++) {
      if ( k == ptype ) {
	int pc = 0;
	for (int n=0; n<header.npart[k]; n++) {
	  file.read((char*)temp, 3*sizeof(float));
	  if (n % numprocs == myid) {
	    particles[pc].vel[0] = temp[0];
	    particles[pc].vel[1] = temp[1];
	    particles[pc].vel[2] = temp[2];
	    pc++;
	  }
	}
      }
      else {
	file.seekg(header.npart[k]*3*sizeof(float), std::ios::cur);
      }
      
    }
    
    file.seekg(sizeof(int), std::ios::cur); // block count
    file.seekg(sizeof(int), std::ios::cur); // block count
    
    // Read particle ID
    //
    for (int k=0; k<6; k++) {
      if ( k == ptype ) {
	int pc = 0;
	for (int n=0; n<header.npart[k]; n++) {
	  int temp;
	  file.read((char*)&temp, sizeof(int));
	  if (n % numprocs == myid) {
	    particles[pc].indx = temp;
	    pc++;
	  }
	}
      }
      else {
	file.seekg(header.npart[k]*sizeof(int), std::ios::cur);
      }
    }
    
    file.seekg(sizeof(int), std::ios::cur); // block count
    
    // Do we need to read mass stanza?
    bool with_mass = std::accumulate(header.mass, header.mass+6, 0.0)>0.0;
    if (with_mass) file.seekg(sizeof(int), std::ios::cur); // block count
    
    for (int k=0; k<6; k++) {
      if ( k == ptype) {
	int pc = 0;
	for (int n=0; n<header.npart[k]; n++) {
	  if (n % numprocs == myid) {
	    if (header.mass[k]==0) {
	      file.read((char*)temp, sizeof(float));
	      particles[pc].mass = temp[0];
	    }
	    else
	      particles[pc].mass = header.mass[k];
	    pc++;
	  }
	}
      }
      else {
	if (header.mass[k]==0)
	  file.seekg(header.npart[k]*sizeof(float), std::ios::cur);
      }
    }
    
    if (with_mass) file.seekg(sizeof(int), std::ios::cur); // block count
    
    // Add other fields, as necessary. Acceleration?
    
    file.close();
    if (myid==0 and _verbose) std::cout << "done." << std::endl;
  }
  
  const Particle* GadgetNative::firstParticle()
  {
    pcount = 0;
    
    return &particles[pcount++];
  }
  
  const Particle* GadgetNative::nextParticle()
  {
    if (pcount < particles.size()) {
      return &particles[pcount++];
    } else {
      if (nextFile()) return firstParticle();
      return 0;
    }
  }
  
  
  std::vector<std::string> GadgetHDF5::Ptypes
  {"Gas", "Halo", "Disk", "Bulge", "Stars", "Bndry"};
  
  std::unordered_map<std::string, int> GadgetHDF5::findP
  { {"Gas", 0}, {"Halo", 1}, {"Disk", 2}, {"Bulge", 3}, {"Stars", 4}, {"Bndry", 5}};
  
  
  GadgetHDF5::GadgetHDF5(const std::vector<std::string>& files, bool verbose)
  {
    _files   = files;
    _verbose = verbose;
    
    ptype = 1;			// Default is halo particles

    getNumbers();
    curfile = _files.begin();

    if (not nextFile()) {
      std::cerr << "GadgetNative: no files found" << std::endl;
    }
  }

  void GadgetHDF5::getNumbers()
  {
    std::set<std::string> pfound;
    std::fill(nptot, nptot+6, 0);

    for (auto file : _files) {

      // Try to catch and HDF5 and parsing errors
      //
      try {
	// Turn off the auto-printing when failure occurs so that we can
	// handle the errors appropriately
	H5::Exception::dontPrint();
	
	const H5std_string FILE_NAME (file);
	const H5std_string GROUP_NAME_what ("/Header");
      
	H5::H5File    file( FILE_NAME, H5F_ACC_RDONLY );
	H5::Group     what(file.openGroup( GROUP_NAME_what ));
      
	// Get time
	{
	  H5::Attribute attr(what.openAttribute("Time"));
	  H5::DataType  type(attr.getDataType());
	
	  attr.read(type, &time);
	}
	
	// Get Mass table
	{
	  H5::Attribute attr(what.openAttribute("MassTable"));
	  H5::DataType  type(attr.getDataType());
	  
	  attr.read(type, mass);
	}
	
	// Get particle counts
	{
	  H5::Attribute attr(what.openAttribute("NumPart_ThisFile"));
	  H5::DataType  type(attr.getDataType());
	
	  attr.read(type, npart);
	}
	
	for (int n=0; n<6; n++) {
	  nptot[n] += npart[n];
	  if (npart[n] > 0) pfound.insert(Ptypes[n]);
	}
	
      }
      // end of try block
      
      // catch failure caused by the H5File operations
      catch(H5::FileIException error)
	{
	  error.printErrorStack();
	}
      
      // catch failure caused by the DataSet operations
      catch(H5::DataSetIException error)
	{
	  error.printErrorStack();
	}
      
      // catch failure caused by the DataSpace operations
      catch(H5::DataSpaceIException error)
	{
	  error.printErrorStack();
	}
    }
      
    Pfound.clear();
    for (auto p : pfound) Pfound.push_back(p);
  }


  bool GadgetHDF5::nextFile()
  {
    if (curfile==_files.end()) return false;
    read_and_load();
    curfile++;
    return true;
  }


  void GadgetHDF5::read_and_load()
  {
    // Try to catch and HDF5 and parsing errors
    //
    try {
      // Turn off the auto-printing when failure occurs so that we can
      // handle the errors appropriately
      H5::Exception::dontPrint();
      
      const H5std_string FILE_NAME (*curfile);
      const H5std_string GROUP_NAME_what ("/Header");
      
      H5::H5File    file( FILE_NAME, H5F_ACC_RDONLY );
      H5::Group     what(file.openGroup( GROUP_NAME_what ));
      
      // Get time
      {
	H5::Attribute attr(what.openAttribute("Time"));
	H5::DataType  type(attr.getDataType());
	
	attr.read(type, &time);
      }
      
      // Get Mass table
      {
	H5::Attribute attr(what.openAttribute("MassTable"));
	H5::DataType  type(attr.getDataType());
	
	attr.read(type, mass);
      }
      
      // Get particle counts
      {
	H5::Attribute attr(what.openAttribute("NumPart_ThisFile"));
	H5::DataType  type(attr.getDataType());
	
	attr.read(type, npart);
      }
      
      if (npart[ptype]>0) {
	std::ostringstream sout;
	sout << "PartType" << ptype;
	
	std::string grpnam = "/" + sout.str();
	H5::Group grp(file.openGroup(grpnam));
	H5::DataSet dataset = grp.openDataSet("Coordinates");
	H5::DataSpace dataspace = dataset.getSpace();
	
	// Get the number of dimensions in the dataspace.
	//
	int rank = dataspace.getSimpleExtentNdims();
	
	// Get the dimension size of each dimension in the dataspace and
	// display them.
	//
	hsize_t dims[2];
	int ndims = dataspace.getSimpleExtentDims( dims, NULL);
	if (myid==0 and _verbose)
	  std::cout << "GadgetHDF5: coordinate rank " << rank
		    << ", dimensions " <<
	    (unsigned long)(dims[0]) << " x " <<
	    (unsigned long)(dims[1]) << std::endl;
	
	// Define the memory space to read dataset.
	//
	H5::DataSpace mspace(rank, dims);
	
	std::vector<float> buf(dims[0]*dims[1]);
	dataset.read(&buf[0], H5::PredType::NATIVE_FLOAT, mspace, dataspace );
	
	if (myid==0 and _verbose)
	  std::cout << "GadgetHDF5: coordinate storage size="
		    << dataset.getStorageSize() << std::endl;
	
	// Set the particle vector
	//
	Particle P;
	
	for (int n=0; n<dims[0]; n++) {
	  if (n % numprocs ==  myid) {
	    P.mass  = mass[ptype];
	    P.level = 0;
	    for (int k=0; k<3; k++) P.pos[k] = buf[n*3+k];
	    particles.push_back(P);
	  }
	}
	
	// Get velocities
	dataspace.close();
	dataset.close();
	dataset = grp.openDataSet("Velocities");
	dataspace = dataset.getSpace();
	
	std::vector<float> vel(dims[0]*dims[1]);
	dataset.read(&vel[0], H5::PredType::NATIVE_FLOAT, mspace, dataspace );
	
	if (myid==0 and _verbose)
	  std::cout << "GadgetHDF5: velocity storage size="
		    << dataset.getStorageSize() << std::endl;
	
	auto it = particles.begin();
	for (int n=0; n<dims[0]; n++) {
	  if (n % numprocs ==  myid) {
	    for (int k=0; k<3; k++) it->vel[k] = vel[n*3+k];
	    it++;
	  }
	}
	
	dataspace.close();
	
	// Try to get Masses.  This will override the assignment from
	// the header if the data exists.
	//
	try {
	  dataset = grp.openDataSet("Masses");
	
	  if (myid==0 and _verbose)
	    std::cout << "GadgetHDF5: mass storage size="
		      << dataset.getStorageSize() << std::endl;
	  
	  if (dataset.getStorageSize()) {
	    
	    dataspace = dataset.getSpace();
	    
	    int rank = dataspace.getSimpleExtentNdims();
	    
	    hsize_t dims[rank];
	    
	    int ndims = dataspace.getSimpleExtentDims(dims, NULL);
	    
	    H5::DataSpace mspace(rank, dims);
	    
	    std::vector<float> masses(dims[0]);
	    dataset.read(&masses[0], H5::PredType::NATIVE_FLOAT, mspace, dataspace );
	  
	    auto it = particles.begin();
	    for (int n=0; n<dims[0]; n++) {
	      if (n % numprocs == myid) (it++)->mass = masses[n];
	    }
	  }
	}
	catch(H5::GroupIException error)
	  {
	    error.printErrorStack();
	  }

	
	dataspace.close();
	dataset.close();
	
	// Try to get particle ids
	//
	dataset = grp.openDataSet("ParticleIDs");
	
	
	if (myid==0 and _verbose)
	  std::cout << "GadgetHDF5: particle ID storage size="
		    << dataset.getStorageSize() << std::endl;
	
	if (dataset.getStorageSize()) {
	  
	  dataspace = dataset.getSpace();
	  
	  int rank = dataspace.getSimpleExtentNdims();
	  
	  hsize_t dims[rank];
	  
	  int ndims = dataspace.getSimpleExtentDims( dims, NULL);
	  
	  H5::DataSpace mspace(rank, dims);
	  
	  std::vector<unsigned> seq(dims[0]);
	  dataset.read(&seq[0], H5::PredType::NATIVE_UINT32, mspace, dataspace );
	  
	  auto it = particles.begin();
	  for (int n=0; n<dims[0]; n++) {
	    if (n % numprocs == myid) (it++)->indx = seq[n];
	  }
	} else {
	  for (int n=0; n<dims[0]; n++) particles[n].indx = n + 1;
	}
      } else {
	std::cerr << "GadgetHDF5:: zero pass particles for type <"
		  << Ptypes[ptype] << ">" << std::endl;
      }
      // end particle type loop
      
    }
    // end of try block
    
    // catch failure caused by the H5File operations
    catch(H5::FileIException error)
      {
	std::ostringstream ost;
	ost << "Error opening HDF5 file: " << *curfile;
	throw std::runtime_error(ost.str());
      }
    
    // catch failure caused by the DataSet operations
    catch(H5::DataSetIException error)
      {
	error.printErrorStack();
      }
    
    // catch failure caused by the DataSpace operations
    catch(H5::DataSpaceIException error)
      {
	error.printErrorStack();
      }
  }
  
  const Particle* GadgetHDF5::firstParticle()
  {
    pcount = 0;
    
    return & particles[pcount++];
  }
  
  const Particle* GadgetHDF5::nextParticle()
  {
    if (pcount < particles.size()) {
      return & particles[pcount++];
    } else {
      if (nextFile()) return firstParticle();
      else return 0;
    }
  }
  
  
  bool badstatus(istream& in)
  {
    ios::iostate i = in.rdstate();
    
    if (i & ios::eofbit) {
      std::cout << "EOF encountered" << std::endl;
      return true;
    }
    else if(i & ios::failbit) {
      std::cout << "Non-Fatal I/O error" << std::endl;;
      return true;
    }  
    else if(i & ios::badbit) {
      std::cout << "Fatal I/O error" << std::endl;
      return true;
    }
    else
      return false;
  }
  
  
  PSPout::PSPout(const std::vector<std::string>& infile, bool verbose) : PSP(verbose)
  {
    
    // Open the file
    // -------------
    try {
      in.open(infile[0]);
    } catch (...) {
      std::ostringstream sout;
      sout << "Could not open PSP file <" << infile[0] << ">";
      throw std::runtime_error(sout.str());
    }
    
    pos = in.tellg();
    
    // Read the header, quit on failure
    // --------------------------------
    try {
      in.read((char *)&header, sizeof(MasterHeader));
    } catch (...) {
      std::ostringstream sout;
      sout << "Could not read master header for <" << infile[0] << ">";
      throw std::runtime_error(sout.str());
    }
    
    for (int i=0; i<header.ncomp; i++) {
      
      PSPstanza stanza;
      unsigned long rsize;
      
      try {
	unsigned long ret;
	in.read((char *)&ret, sizeof(unsigned long));
	
	rsize = sizeof(double);
	if ( (ret & nmask) == magic ) {
	  rsize = ret & mmask;
	}
      } catch (...) {
	std::ostringstream sout;
	sout << "Error reading magic for <" << infile[0] << ">";
	throw std::runtime_error(sout.str());
      }
      
      try {
	stanza.comp.read(&in);
      } catch (...) {
	std::ostringstream sout;
	sout << "Error reading component header for <" << infile[0] << ">";
	throw std::runtime_error(sout.str());
      }
      
      stanza.pspos = in.tellg();
      
      // Parse the info string
      // ---------------------
      std::istringstream sin(stanza.comp.info.get());
      YAML::Node conf, cconf, fconf;
      bool yaml_ok = true;
      
      std::ostringstream errors;
      
      try {
	conf = YAML::Load(sin);
      }
      catch (YAML::Exception & error) {
	errors << "Error parsing component config.  Trying old-style PSP"
	       << std::endl
	       << error.what() << std::endl;
	yaml_ok = false;
      }
      
      if (yaml_ok) {
	
	cconf  = conf["parameters"];
	fconf  = conf["force"];
	
	// Output map in flow style
	//
	cconf.SetStyle(YAML::EmitterStyle::Flow);
	if (fconf.IsMap())
	  fconf["parameters"].SetStyle(YAML::EmitterStyle::Flow);
	
	// Write node to sstream
	//
	std::ostringstream csout, fsout;
	csout << cconf;
	if (fconf.IsMap())
	  fsout << fconf["parameters"];
	else
	  fsout << "<undefined>";
	
	stanza.name       = conf["name"].as<std::string>(); 
	if (fconf.IsMap())
	  stanza.id       = fconf["id"].as<std::string>();
	else
	  stanza.id       = "<undefined>";
	stanza.cparam     = csout.str();
	stanza.fparam     = fsout.str();
	stanza.index_size = 0;
	stanza.r_size     = rsize;
	
	// Check for indexing
	// -------------------
	size_t pos1 = stanza.cparam.find("indexing");
	if (cconf["indexing"]) {
	  if (cconf["indexing"].as<bool>()) 
	    stanza.index_size = sizeof(unsigned long);
	}
	
      } // END: new PSP
      else {
	
	// Parse the info string
	// ---------------------
	StringTok<string> tokens(stanza.comp.info.get());
	stanza.name       = trim_copy(tokens(":"));
	stanza.id         = trim_copy(tokens(":"));
	stanza.cparam     = trim_copy(tokens(":"));
	stanza.fparam     = trim_copy(tokens(":"));
	stanza.index_size = 0;
	stanza.r_size     = rsize;
	
	// Check for indexing
	// -------------------
	size_t pos1 = stanza.cparam.find("indexing");
	if (pos1 != string::npos) {
	  // Look for equals sign
	  size_t pos2 = stanza.cparam.find("=", pos1);
	  
	  // No equals sign?!!
	  if (pos2 == string::npos) {
	    cerr << "Bad syntax in component parameter string" << std::endl;
	    exit(-1);
	  }
	  
	  // Look for field delimiter
	  size_t pos3 = stanza.cparam.find(",", pos2);
	  if (pos3 != string::npos) pos3 -= pos2+1;
	  
	  if (atoi(stanza.cparam.substr(pos2+1, pos3).c_str()))
	    stanza.index_size = sizeof(unsigned long);
	}
	
      }
      // END: old PSP
      
      
      // Skip forward to next header
      // ---------------------------
      try {
	in.seekg(stanza.comp.nbod*(stanza.index_size                +
				   8*stanza.r_size                  + 
				   stanza.comp.niatr*sizeof(int)    +
				   stanza.comp.ndatr*stanza.r_size
				   ), ios::cur);
      } 
      catch(...) {
	std::cout << "IO error: can't find next header for time="
		  << header.time << " . . . quit reading <" << infile[0] << ">";
	break;
      }
      
      stanzas.push_back(stanza);
      
      if (verbose) {
	std::cout << errors.str();
      }
    }
    
    spos = stanzas.begin();
  }
  
  
  PSPspl::PSPspl(const std::vector<std::string>& master, bool verbose) : PSP(verbose)
  {
    
    // Open the file
    // -------------
    try {
      in.open(master[0]);
    } catch (...) {
      std::ostringstream sout;
      sout << "Could not open the master SPL file <" << master[0] << ">";
      throw std::runtime_error(sout.str());
    }
    
    if (!in.good()) {
      std::ostringstream sout;
      sout << "Error opening master SPL file <" << master[0] << ">";
      throw std::runtime_error(sout.str());
    }
    
    // Read the header, quit on failure
    // --------------------------------
    try {
      in.read((char *)&header, sizeof(MasterHeader));
    } catch (...) {
      std::ostringstream sout;
      sout << "Could not read master header for <" << master[0] << ">";
      throw std::runtime_error(sout.str());
    }
    
    for (int i=0; i<header.ncomp; i++) {
      
      unsigned long   rsize;
      unsigned long   cmagic;
      int             number;
      PSPstanza       stanza;
      
      try {
	in.read((char*)&cmagic,   sizeof(unsigned long));
	in.read((char*)&number, sizeof(int));
      } catch (...) {
	std::ostringstream sout;
	sout << "Error reading magic info for Comp #" << i << " from <"
	     << master[0] << ">";
	throw std::runtime_error(sout.str());
      }
      
      rsize = sizeof(double);
      if ( (cmagic & nmask) == magic ) {
	rsize = cmagic & mmask;
      }
      
      try {
	stanza.comp.read(&in);
      } catch (...) {
	std::ostringstream sout;
	sout << "Error reading component header Comp #" << i << " from <"
	     << master[0] << ">";
	throw std::runtime_error(sout.str());
      }
      
      // Parse the info string
      // ---------------------
      std::istringstream sin(stanza.comp.info.get());
      YAML::Node conf, cconf, fconf;
      
      try {
	conf = YAML::Load(sin);
      }
      catch (YAML::Exception & error) {
	std::ostringstream sout;
	sout << "Error parsing component config in Comp #" << i
	     << " from <" << master[0] << ">";
	throw std::runtime_error(sout.str());
      }
      
      cconf  = conf["parameters"];
      fconf  = conf["force"];
      
      // Output map in flow style
      //
      cconf.SetStyle(YAML::EmitterStyle::Flow);
      if (fconf.IsMap())
	fconf["parameters"].SetStyle(YAML::EmitterStyle::Flow);
      
      // Write node to sstream
      //
      std::ostringstream csout, fsout;
      csout << cconf;
      if (fconf.IsMap())
	fsout << fconf["parameters"];
      else
	fsout << "<undefined>";
      
      stanza.name       = conf["name"].as<std::string>(); 
      if (fconf.IsMap())
	stanza.id       = fconf["id"].as<std::string>();
      else
	stanza.id       = "<undefined>";
      stanza.cparam     = csout.str();
      stanza.fparam     = fsout.str();
      stanza.index_size = 0;
      stanza.r_size     = rsize;
      
      // Check for indexing
      // -------------------
      size_t pos1 = stanza.cparam.find("indexing");
      if (cconf["indexing"]) {
	if (cconf["indexing"].as<bool>()) 
	  stanza.index_size = sizeof(unsigned long);
      }
      
      // Get file names for parts
      // ------------------------
      std::vector<std::string> parts(number);
      
      const size_t PBUF_SIZ = 1024;
      char buf [PBUF_SIZ];
      
      for (int n=0; n<number; n++) {
	in.read((char *)buf, PBUF_SIZ);
	stanza.nparts.push_back(buf);
      }
      
      stanzas.push_back(stanza);
    }
    // END: component loop
    
    // Close current file
    //
    if (in.is_open()) in.close();
    
    spos = stanzas.begin();
  }
  
  
  void PSP::PrintSummary(ostream &out, bool stats, bool timeonly)
  {
    out << "Time=" << header.time << std::endl;
    if (!timeonly) {
      out << "   Total particle number: " << header.ntot  << std::endl;
      out << "   Number of components:  " << header.ncomp << std::endl;
      
      int cnt=1;
      
      for (auto s : stanzas) {
	
	// Print the info for this stanza
	// ------------------------------
	out << std::setw(60) << std::setfill('-') << "-" << std::endl << std::setfill(' ');
	out << "--- Component #" << std::setw(2) << cnt++          << std::endl;
	out << std::setw(20) << " name :: "      << s.name         << std::endl
	    << std::setw(20) << " id :: "        << s.id           << std::endl
	    << std::setw(20) << " cparam :: "    << s.cparam       << std::endl
	    << std::setw(20) << " fparam :: "    << s.fparam       << std::endl
	    << std::setw(20) << " nbod :: "      << s.comp.nbod    << std::endl
	    << std::setw(20) << " niatr :: "     << s.comp.niatr   << std::endl
	    << std::setw(20) << " ndatr :: "     << s.comp.ndatr   << std::endl
	    << std::setw(20) << " rsize :: "     << s.r_size       << std::endl;
	out << std::setw(60) << std::setfill('-')     << "-" << std::endl << std::setfill(' ');
	if (stats) {
	  ComputeStats();
	  out << std::endl << std::setw(20) << "*** Position" 
	      << std::setw(15) << "X" << std::setw(15) << "Y" << std::setw(15) << "Z"
	      << std::endl;
	  out << std::setw(20) << "Min :: ";
	  for (unsigned k=0; k<3; k++) out << std::setw(15) << pmin[k];
	  out << std::endl;
	  out << std::setw(20) << "Med :: ";
	  for (unsigned k=0; k<3; k++) out << std::setw(15) << pmed[k];
	  out << std::endl;
	  out << std::setw(20) << "Max :: ";
	  for (unsigned k=0; k<3; k++) out << std::setw(15) << pmax[k];
	  out << std::endl;
	  out << std::endl << std::setw(20) << "*** Velocity"
	      << std::setw(15) << "U" << std::setw(15) << "Vn" << std::setw(15) << "W"
	      << std::endl;
	  out << std::setw(20) << "Min :: ";
	  for (unsigned k=0; k<3; k++) out << std::setw(15) << vmin[k];
	  out << std::endl;
	  out << std::setw(20) << "Med :: ";
	  for (unsigned k=0; k<3; k++) out << std::setw(15) << vmed[k];
	  out << std::endl;
	  out << std::setw(20) << "Max :: ";
	  for (unsigned k=0; k<3; k++) out << std::setw(15) << vmax[k];
	  out << std::endl;
	}      
      }
    }
  }
  
  PSPstanza* PSP::GetStanza()
  {
    spos = stanzas.begin();
    cur  = &(*spos);
    if (spos != stanzas.end()) 
      return cur;
    else 
      return 0;
  }
  
  PSPstanza* PSP::NextStanza()
  {
    spos++;
    cur = &(*spos);
    if (spos != stanzas.end()) 
      return cur;
    else 
      return 0;
  }
  
  const Particle* PSPout::firstParticle()
  {
    pcount = 0;
    
    in.seekg(cur->pspos);
    
    return nextParticle();
  }
  
  const Particle *PSPout::nextParticle()
  {
    badstatus(in);		// DEBUG
    
    // Stagger on first read
    // ---------------------
    
    if (pcount==0) {
      for (int n=0; n<myid; n++) {
	if (pcount < spos->comp.nbod)
	  if (spos->r_size == 4)
	    fpart.skip(in, pcount++, spos);
	  else
	    dpart.skip(in, pcount++, spos);
      }
    }
    
    // Read partcle
    // ------------
    if (pcount < spos->comp.nbod) {
      
      if (spos->r_size == 4)
	fpart.read(in, pcount++, spos);
      else
	dpart.read(in, pcount++, spos);
      
      // Stride by numprocs-1
      // --------------------
      for (int n=0; n<numprocs-1; n++)  {
	if (pcount < spos->comp.nbod) {
	  if (spos->r_size == 4)
	    fpart.skip(in, pcount++, spos);
	  else
	    dpart.skip(in, pcount++, spos);
	}
      }
      
      if (spos->r_size == 4)
	return static_cast<Particle*>(&fpart);
      else
	return static_cast<Particle*>(&dpart);
      
    } else
      return 0;
  }
  
  const Particle* PSPspl::firstParticle()
  {
    pcount = 0;
    
    // Set iterator to beginning of vector
    fit = spos->nparts.begin();
    
    // Open next file in sequence
    openNextBlob();
    
    return nextParticle();
  }
  
  void PSPspl::openNextBlob()
  {
    // Close current file
    //
    if (in.is_open()) in.close();
    
    std::string curfile(*fit);
    
    try {
      in.open(curfile);
    } catch (...) {
      std::ostringstream sout;
      sout << "Could not open SPL blob <" << curfile << ">";
      throw std::runtime_error(sout.str());
    }
    
    if (not in.good()) {
      std::ostringstream sout;
      sout << "Could not open SPL blob <" << curfile << ">";
      throw std::runtime_error(sout.str());
    }
    
    try {
      in.read((char*)&N, sizeof(unsigned int));
    } catch (...) {
      std::ostringstream sout;
      sout << "Could not get particle count from <" << curfile << ">";
      throw std::runtime_error(sout.str());
    }
    
    fcount = 0;
    
    // Advance filename iterator
    fit++;
  }
  
  const Particle* PSPspl::nextParticle()
  {
    badstatus(in);		// DEBUG
    
    // Stagger on first read
    // ---------------------
    if (pcount==0) {
      for (int n=0; n<myid; n++) {
	if (pcount < spos->comp.nbod) {
	  if (fcount==N) openNextBlob();
	  if (spos->r_size == 4)
	    fpart.skip(in, pcount++, spos);
	  else
	    dpart.skip(in, pcount++, spos);
	  fcount++;
	}
      }
    }
    
    // Read partcle
    // ------------
    if (pcount < spos->comp.nbod) {
      
      if (fcount==N) openNextBlob();
      
      // Read blob
      if (spos->r_size == 4)
	fpart.read(in, pcount++, spos);
      else
	dpart.read(in, pcount++, spos);
      fcount++;
      
      // Stride by numprocs-1
      // --------------------
      for (int n=0; n<numprocs-1; n++)  {
	if (pcount < spos->comp.nbod) {
	  if (fcount==N) openNextBlob();
	  if (spos->r_size == 4)
	    fpart.skip(in, pcount++, spos);
	  else
	    dpart.skip(in, pcount++, spos);
	  fcount++;
	}
      }
      
      if (spos->r_size == 4)
	return static_cast<Particle*>(&fpart);
      else
	return static_cast<Particle*>(&dpart);
      
    } else
      return 0;
  }
  
  
  void PSP::ComputeStats()
  {
    cur = &(*spos);
    
    // Initialize lists
    std::vector< std::vector<double> > plist(3, std::vector<double>(spos->comp.nbod) );
    std::vector< std::vector<double> > vlist(3, std::vector<double>(spos->comp.nbod) );
    mtot = 0.0;
    
    auto P = firstParticle();
    unsigned n=0;
    while (P) {
      mtot += P->mass;
      for (unsigned k=0; k<3; k++) {
	plist[k][n] = P->pos[k];
	vlist[k][n] = P->vel[k];
      }
      P = nextParticle();
      n++;
    }
    
    pmin = vector<double>(3);
    pmed = vector<double>(3);
    pmax = vector<double>(3);
    
    vmin = vector<double>(3);
    vmed = vector<double>(3);
    vmax = vector<double>(3);
    
    for (unsigned k=0; k<3; k++) {
      std::sort(plist[k].begin(), plist[k].end());
      pmin[k] = plist[k].front();
      pmed[k] = plist[k][floor(0.5*spos->comp.nbod+0.5)];
      pmax[k] = plist[k].back();
      std::sort(vlist[k].begin(), vlist[k].end());
      vmin[k] = vlist[k].front();
      vmed[k] = vlist[k][floor(0.5*spos->comp.nbod+0.5)];
      vmax[k] = vlist[k].back();
    }
  }
  
  
  void PSP::writePSP(std::ostream& out,  bool real4)
  {
    // Write master header
    // --------------------
    out.write((char *)&header, sizeof(MasterHeader));
    
    // Write each component
    // --------------------
    for (auto its=stanzas.begin(); its!=stanzas.end(); its++) 
      write_binary(out, its, real4);
  }
  
  void PSP::write_binary(std::ostream& out,
			 list<PSPstanza>::iterator its, bool real4)
  {
    spos = its;
    
    // Write magic #
    unsigned long cmagic = magic;
    if (real4) cmagic += sizeof(float);
    else       cmagic += sizeof(double);
    out.write((const char *)&cmagic, sizeof(unsigned long));
    
    
    // Write component header
    its->comp.write(&out);
    
    // Position the stream
    if (its->pspos) in.seekg(its->pspos, ios::beg);
    
    
    // Write each particle
    unsigned count = 0;
    bool indexing = false;
    if (its->index_size>0) indexing = true;
    
    for (auto part=firstParticle(); part!=0; part=nextParticle()) 
      {
	part->writeBinary(sizeof(float), indexing, &out);
	count++;
      }
    
    if (VERBOSE)
      std::cerr << std::string(72, '-') << std::endl
		<< "Wrote " << count << " particles "<< std::endl
		<< std::string(72, '-') << std::endl
		<< spos->comp << std::endl
		<< std::string(72, '-') << std::endl;
  }
  
  
  std::vector<std::string> ParticleReader::readerTypes
  {"PSPout", "PSPspl", "GadgetNative", "GadgetHDF5", "TipsyNative", "TipsyXDR", "Bonsai"};
  
  
  std::vector<std::vector<std::string>>
  ParticleReader::parseFileList
  (const std::string& file, const std::string& delimit)
  {
    std::vector<std::vector<std::string>> batches;

    std::ifstream in(file);
    if (in) {
      std::vector<std::string> files;
      
      std::string name;
      while(in >> name) files.push_back(name);
      
      std::sort(files.begin(), files.end());
      
      std::vector<std::string> batch;
      std::string templ;
      
      for (auto f : files) {
	std::size_t found = f.find_last_of(delimit);

	// No delimiter?
	if (found == std::string::npos) {
	  batch.push_back(f);
	  batches.push_back(batch);
	  batch.clear();
	}
	// Found a delimiter
	else {
	  auto trimmed = f.substr(0, found);

	  if (batch.size()==0) {
	    templ = trimmed;
	    batch.push_back(f);
	  }
	  else if (trimmed == templ) {
	    batch.push_back(f);
	  }
	  else {		// Mismatch: new batch
	    if (batch.size()) {
	      batches.push_back(batch);
	      batch.clear();
	    }
	    templ = trimmed;
	    batch.push_back(f);
	  }
	}
      }
    } else {
      std::cerr << "Error opening file <" << file << ">" << std::endl;
    }

    return batches;
  }

  std::shared_ptr<ParticleReader>
  ParticleReader::createReader(const std::string& reader,
			       const std::vector<std::string>& file,
			       int myid, bool verbose)
  {
    if (reader.find("PSPout") == 0)
      return std::make_shared<PSPout>(file, verbose);
    else if (reader.find("PSPspl") == 0)
      return std::make_shared<PSPspl>(file, verbose);
    else if (reader.find("GadgetNative") == 0)
      return std::make_shared<GadgetNative>(file, verbose);
    else if (reader.find("GadgetHDF5") == 0)
      return std::make_shared<GadgetHDF5>(file, verbose);
    else if (reader.find("TipsyNative") == 0)
      return std::make_shared<Tipsy>(file, Tipsy::TipsyType::native, verbose);
    else if (reader.find("TipsyXDR") == 0)
#ifdef HAVE_XDR
      return std::make_shared<Tipsy>(file, Tipsy::TipsyType::xdr, verbose);
#else
    {
      if (myid==0) {
	std::cout << "ParticleReader: your build does not have RPC/XDR " 
		  << "support so Tipsy standard reading is not available."
		  << std::endl
		  << "Try installing the libtirpc package for your "
		  << "distribution or from Sourceforge directly or "
		  << "use Tipsy native format." << std::endl;
      }
      exit(1);
    }
#endif

    else if (reader.find("Bonsai") == 0)
      return std::make_shared<Tipsy>(file, Tipsy::TipsyType::bonsai, verbose);
    else {
      if (myid==0) {
	std::cout << "ParticleReader: I don't know about reader <" << reader
		  << ">" << std::endl
		  << "Available readers are:";
	for (auto s : readerTypes) std::cout << " " << s;
	std::cout << std::endl;
      }
      exit(1);
    }
  }

  std::vector<std::string> Tipsy::Ptypes
  {"Gas", "Dark", "Star"};
  
  std::unordered_map<std::string, int> Tipsy::findP
  { {"Gas", 0}, {"Dark", 1}, {"Star", 2} };
  
  
  void Tipsy::getNumbers()
  {
    Ngas = Ndark = Nstar = 0;
    std::set<std::string> types;

    for (auto file : files) {

      // Make a tipsy native reader
      if (ttype == TipsyType::native)
	ps = std::make_shared<TipsyReader::TipsyNative>(*curfile, false);
      // Native tipsy with ID conversion
      else if (ttype == TipsyType::bonsai)
	ps = std::make_shared<TipsyReader::TipsyNative>(*curfile, true);
      // Make a tipsy xdr reader
      else {
#ifdef HAVE_XDR
	ps = std::make_shared<TipsyReader::TipsyXDR>(*curfile);
#else
	ps = std::make_shared<TipsyReader::TipsyNative>(*curfile, false);
#endif
      }

      if (ps->gas_particles.size() ) {
	types.insert("Gas");
	Ngas += ps->gas_particles.size();
      }

      if (ps->dark_particles.size()) {
	types.insert("Dark");
	Ndark += ps->dark_particles.size();
      }

      if (ps->star_particles.size()) {
	types.insert("Star");
	Nstar += ps->star_particles.size();
      }

    }

    curTypes.clear();
    for (auto s : types) curTypes.push_back(s);
  }

  bool Tipsy::nextFile()
  {
    if (curfile==files.end()) return false;
    
    if (ttype == TipsyType::native)
      ps = std::make_shared<TipsyReader::TipsyNative>(*curfile, false);
    else if (ttype == TipsyType::bonsai)
      ps = std::make_shared<TipsyReader::TipsyNative>(*curfile, true);
    else {
#ifdef HAVE_XDR
      ps = std::make_shared<TipsyReader::TipsyXDR>(*curfile);
#else
      ps = std::make_shared<TipsyReader::TipsyNative>(*curfile), false;
#endif
    }

    ps->readParticles();
    curfile++;
    return true;
  }

  Tipsy::Tipsy(const std::string& file, TipsyType Type,
	       bool verbose)
  {
    ttype = Type;
    files.push_back(file);
    getNumbers();
    curfile = files.begin();
    if (not nextFile()) {
      std::cerr << "Tipsy: no files found" << std::endl;
    }
  }
  
  Tipsy::Tipsy(const std::vector<std::string>& filelist, TipsyType Type,
	       bool verbose)
  {
    ttype = Type;
    files = filelist;
    getNumbers();
    curfile = files.begin();
    if (not nextFile()) {
      std::cerr << "Tipsy: no files found" << std::endl;
    }
  }
  
  void Tipsy::SelectType(const std::string& name)
  {
    if (std::find(curTypes.begin(), curTypes.end(), name) == curTypes.end()) {
      std::cout << "Tipsy error: no particle type <" << name << ">"
		<< std::endl;
    } else {
      curName = name;
      ps->readParticles();
    }

    curfile = files.begin();	// Set to first file and open
    nextFile();
  }

  void Tipsy::packParticle()
  {
    if (curName=="Gas") {
      P.mass  = ps->gas_particles[pcount].mass;
      P.level = 0;
      for (int k=0; k<3; k++) {
	P.pos[k] = ps->gas_particles[pcount].pos[k];
	P.vel[k] = ps->gas_particles[pcount].vel[k];
      }
      if (ttype == TipsyType::bonsai) P.indx = ps->gas_particles[pcount].ID();
      pcount++;
      return;
    }
      
    if (curName=="Dark") {
      P.mass  = ps->dark_particles[pcount].mass;
      P.level = 0;
      for (int k=0; k<3; k++) {
	P.pos[k] = ps->dark_particles[pcount].pos[k];
	P.vel[k] = ps->dark_particles[pcount].vel[k];
      }
      if (ttype == TipsyType::bonsai) P.indx = ps->dark_particles[pcount].ID();
      pcount++;
      return;
    }
      
    if (curName=="Star") {
      P.mass  = ps->star_particles[pcount].mass;
      P.level = 0;
      for (int k=0; k<3; k++) {
	P.pos[k] = ps->star_particles[pcount].pos[k];
	P.vel[k] = ps->star_particles[pcount].vel[k];
      }
      if (ttype == TipsyType::bonsai) P.indx = ps->star_particles[pcount].ID();
      pcount++;
      return;
    }
    
    std::cerr << "Tipsy: logic error" << std::endl;
    exit(-1);
  }
  
  unsigned long Tipsy::CurrentNumber()
  {
    if (curName=="Gas") {
      return ps->gas_particles.size();
    } else if (curName=="Dark") {
      return ps->dark_particles.size();
    } else if (curName=="Star") {
      return ps->star_particles.size();
    } else {
      return 0;
    }
  }
  
  const Particle* Tipsy::firstParticle()
  {
    pcount = 0;
    packParticle();
    return &P;
  }
    
  const Particle* Tipsy::nextParticle()
  {
    if (curName=="Gas"  and pcount==ps->gas_particles.size()) {
      if (nextFile()) return firstParticle();
      else return NULL;
    }
    
    if (curName=="Dark" and pcount==ps->dark_particles.size()) {
      if (nextFile()) return firstParticle();
      else return NULL;
    }

    if (curName=="Star" and pcount==ps->star_particles.size()) {
      if (nextFile()) return firstParticle();
      else return NULL;
    }

    packParticle();
    return &P;
  }

}
// END: PR namespace
