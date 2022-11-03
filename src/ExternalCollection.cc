#include "expand.H"

#include <dirent.h>
#include <dlfcn.h> 
#include <stdlib.h> 
#include <unistd.h> 

#include <iostream> 
#include <fstream> 
#include <map> 
#include <list> 
#include <vector>
#include <string> 

#include <ExternalCollection.H>

#include <tidalField.H>
#include <externalShock.H>
#include <generateRelaxation.H>
#include <ScatterMFP.H>
#include <HaloBulge.H>

#if DSMC_ENABLED>0
#include <TreeDSMC.H>
#endif

ExternalCollection::ExternalCollection(void)
{
				// Do nothing
}

void ExternalCollection::initialize()
{
  (*barrier)("ExternalCollection::initialize: BEGIN", __FILE__, __LINE__);

  dynamicload();
  
  YAML::Node ext;

  try {
    ext = parse["External"];
  }
  catch (YAML::Exception & error) {
    if (myid==0) std::cout << "Error parsing 'external' stanza: "
			   << error.what() << std::endl
			   << std::string(60, '-') << std::endl
			   << "Config node"        << std::endl
			   << std::string(60, '-') << std::endl
			   << parse                << std::endl
			   << std::string(60, '-') << std::endl;
    MPI_Finalize();
    exit(-1);
  }

  if (ext.IsSequence()) {

    int next = 0;

    while (ext[next]) {

      std::string name = ext[next]["id"].as<std::string>();
      const YAML::Node& node = ext[next]["parameters"];
      
      // Instantiate the force
      if ( !name.compare("tidalField") )
	
	force_list.insert(force_list.end(), new tidalField(node));
      
      else if ( !name.compare("externalShock") )
	
	force_list.insert(force_list.end(), new externalShock(node));
      
      else if ( !name.compare("generateRelaxation") )
	
	force_list.insert(force_list.end(), new generateRelaxation(node));
      
      else if ( !name.compare("ScatterMFP") )

	force_list.insert(force_list.end(), new ScatterMFP(node));
      
#if DSMC_ENABLED>0
      else if ( !name.compare("TreeDSMC") )
	
	force_list.insert(force_list.end(), new TreeDSMC(node));
#endif
      else if ( !name.compare("HaloBulge") )
	
	force_list.insert(force_list.end(), new HaloBulge(node));
      
      else {			// First try to find it in dynamic library
    
	bool found = false;

	for (fitr=factory.begin(); fitr!=factory.end(); fitr++) {

	  if (!name.compare(fitr->first)) {
	    force_list.insert(force_list.end(), factory[name](node));
	    found = true;
	  }
	  
	}
	

	if (!found) {		// Complain to user

	  string msg("I don't know about the external force named: ");
	  msg += name;
	  throw GenericError(msg, __FILE__, __LINE__, 1021, false);
	}
      }
      
      next++;
    }
  } else {
    if (myid==0)
      std::cout << std::string(72, '-') << std::endl
		<< "No external force entries" << std::endl
		<< std::string(72, '-') << std::endl;
  }

  (*barrier)("ExternalCollection::initialize: FINISH", __FILE__, __LINE__);
}

ExternalCollection::~ExternalCollection(void)
{

				// destroy any forces we created
  int i = 0;
  for(sitr=force_list.begin(); sitr!=force_list.end(); sitr++) {
#ifdef DEBUG
    cout << "Process " << myid << ": deleting <" << ++i << "> . . .";
#endif
    delete *sitr;
#ifdef DEBUG
    cout << " done" << endl;
#endif
  }

				// close all the dynamic libs we opened
  i = 0;
  for(itr=dl_list.begin(); itr!=dl_list.end(); itr++) {
    void *dlib = *itr;
    if (dlib) {
      dlclose(dlib);
    }
  }

}


vector<string> ExternalCollection::getlibs(void)
{
  vector<string> ret;
  struct dirent **namelist;

  int n = scandir(ldlibdir.c_str(), &namelist, 0, alphasort);
  if (n < 0)
    perror("scandir");
  else {
    for (int i=0; i<n; i++) {
      char *s = namelist[i]->d_name;
      unsigned sz = strlen(s);
      if (sz>3) {
	if (strncmp(s, "lib", 3) ==0 && 
	    strncmp(&s[sz-3], ".so", 3)==0  ) 
	  ret.push_back(namelist[i]->d_name);
      }
      free(namelist[i]);
    }
    free(namelist);
  }
  
  return ret;
}

void ExternalCollection::dynamicload(void)
{
#ifdef DEBUG
  ostringstream ostr;
  ostr << "extcoll.log." << myid;
  ofstream tout(ostr.str().c_str());
#endif

  if (myid==0) cout << endl << "ExternalCollection:" << endl
		    << setw(71) << setfill('-') << "-" << endl
		    <<  "Loaded user libraries <";
  void *dlib; 
  bool first = true
;
  for (int i=0; i<numprocs; i++) {
    //
    // Put dlopen in a loop with an barrier to prevent swamping 
    // slow NFS server problems that seem to plague our cluster
    //
    if (i==myid) {		
      //
      // get the names of all the dynamic libs (.so files)
      //
      vector<string> liblist = getlibs();
      vector<string>::iterator s = liblist.begin();
      while (s != liblist.end()) {
	//
	// preappend ldlibdir to the front of the lib name
	//
	ostringstream name;
	name << ldlibdir << "/" << *s;
#ifdef DEBUG
	tout << "Process " << myid << ": call dlopen on <" 
	     << name.str() << ">" << endl;
#endif
	dlib = dlopen(name.str().c_str(), RTLD_NOW | RTLD_GLOBAL);
	if (dlib == NULL) {
	  std::cout << std::endl
		    << std::string(80, '-') << std::endl
		    << "Error opening <" << name.str() << ">: dlerror="
		    << dlerror() << std::endl
		    << std::string(80, '-') << std::endl;
	  //
	  // open the memory map? (Only set to true for super deep
	  // library debug which is why I'm providing no flag)
	  //
	  if (false) {
	    std::ifstream map("/proc/self/maps");
	    if (map.is_open()) {
	      std::string line;
	      while (getline(map, line)) {
		std::cout << line << std::endl;
	      }
	      std::cout << std::flush;
	    } else {
	      std::cout << "ExternalCollection: "
			<< "could not open /proc/self/maps" << std::endl;
	    }
	  }
	  std::string emsg = "Error opening <" + name.str() + ">";
	  throw GenericError(emsg.c_str(), __FILE__, __LINE__, 1022);
	}
	if (myid==0) {
	  if (first) {
	    cout << *s;
	    first = false;
	  } else 
	    cout << " " << *s;
	}
				// add the handle to our list
	dl_list.insert(dl_list.end(), dlib);
	s++;
      }
    }
    MPI_Barrier(MPI_COMM_WORLD);
  }

#ifdef DEBUG
  tout << "Process " << myid << ": done with dlopen" << endl;
  tout.close();
#endif

  if (myid==0) {
    cout << ">" << endl 
	 << setw(71) << setfill('-') << "-" << endl
	 << "Available user routines are <";
    first = true;
    for (fitr=factory.begin(); fitr!=factory.end(); fitr++) {
      if (first) {
	cout << fitr->first;
	first = false;
      } else
	cout << " " << fitr->first;
    }

    cout << ">" << endl 
	 << setw(71) << setfill('-') << "-" << setfill(' ') << endl << endl;
  }
}

void ExternalCollection::finish()
{ 
  for (auto p : force_list) p->finish(); 
}

