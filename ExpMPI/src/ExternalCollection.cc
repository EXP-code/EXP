#include "expand.h"

#include <dlfcn.h> 
#include <stdio.h> 
#include <unistd.h> 

#include <iostream> 
#include <fstream> 
#include <map> 
#include <list> 
#include <vector>
#include <string> 

#include <ExternalCollection.H>

ExternalCollection::ExternalCollection(void)
{
				// Do nothing
}

void ExternalCollection::initialize()
{
  spair data;

  dynamicload();
  
  parse->find_list("external");

  while (parse->get_next(data)) {
      
    string name = data.first;
    string rest = data.second;
    
				// Instantiate the force
    if ( !name.compare("tidalField") )

      force_list.insert(force_list.end(), new tidalField(rest));

    else if ( !name.compare("externalShock") )
      
      force_list.insert(force_list.end(), new externalShock(rest));

    else if ( !name.compare("generateRelaxation") )

      force_list.insert(force_list.end(), new generateRelaxation(rest));

    else if ( !name.compare("ScatterMFP") )

      force_list.insert(force_list.end(), new ScatterMFP(rest));
    
    else if ( !name.compare("HaloBulge") )

      force_list.insert(force_list.end(), new HaloBulge(rest));
    
    else {			// First try to find it in dynamic library
    
      bool found = false;

      for (fitr=factory.begin(); fitr!=factory.end(); fitr++) {

	if (!name.compare(fitr->first)) {
	  force_list.insert(force_list.end(), factory[name](rest));
	  found = true;
	}

      }


      if (!found) {		// Complain to user
	string msg("I don't know about the external force named: ");
	msg += name;
	bomb(msg);
      }
    }

  }

}

ExternalCollection::~ExternalCollection(void)
{

				// close all the dynamic libs we opened
  for(itr=dl_list.begin(); itr!=dl_list.end(); itr++) {
    void *dlib = *itr;
    if (dlib) dlclose(dlib);
  }
				// destroy any forces we created
  for(sitr=force_list.begin(); sitr!=force_list.end(); sitr++) {
    delete *sitr;
  }
}

void ExternalCollection::dynamicload(void)
{
  // size of buffer for reading in directory entries 
  const unsigned int BUF_SIZE = 1024;

				// command
  string command;
  command = "cd " + ldlibdir + "; ls *.so";

				// string to get dynamic lib names
  char in_buf[BUF_SIZE];	// input buffer for lib
  

				// get the names of all the dynamic libs (.so 
				// files) in the current dir 
  FILE *dl = popen(command.c_str(), "r");
  if(!dl) {
    perror("popen");
    exit(-1);
  }
  
  if (myid==0) cout << "ExternalCollection: loaded user libraries <";

  void *dlib; 
  char name[1024];
  while(fgets(in_buf, BUF_SIZE, dl)) {
				// trim off the whitespace 
    char *ws = strpbrk(in_buf, " \t\n"); 
    if (ws) *ws = '\0';
				// preappend ldlibdir to the front of 
				// the lib name
    sprintf(name, "%s/%s", ldlibdir.c_str(), in_buf); 
    dlib = dlopen(name, RTLD_NOW | RTLD_GLOBAL);
    if(dlib == NULL) {
      cerr << dlerror() << endl; 
      exit(-1);
    }

    if (myid==0) cout << " " << in_buf;
    
				// add the handle to our list
    dl_list.insert(dl_list.end(), dlib);
  }

  if (myid==0) cout << " >\n\n";
}

void ExternalCollection::bomb(string& msg)
{
  cerr << "ExternalCollection: " << msg << endl;
  exit(-1);
}
