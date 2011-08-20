/*
  Compute accelerations, potential, and density.
*/

#include <expand.h>

#include <vector>

#include <ComponentContainer.H>
#include <ExternalCollection.H>
#include <StringTok.H>

#ifdef USE_GPTL
#include <gptl.h>
#endif


long ComponentContainer::tinterval = 300;	// Seconds between timer dumps

ComponentContainer::ComponentContainer(void)
{
  gottapot      = false;
  gcom1         = new double [3];
  gcov1         = new double [3];

  timing        = false;
  thread_timing = false;
  state         = NONE;

  // Fine resolution for these timers (default resolution is 1 sec)
  //
  timer_posn.	Microseconds();
  timer_gcom.	Microseconds();
  timer_angmom.	Microseconds();
  timer_zero.	Microseconds();
  timer_accel.	Microseconds();
  timer_thr_acc.Microseconds();
  timer_thr_int.Microseconds();
  timer_thr_ext.Microseconds();
  timer_inter.	Microseconds();
  timer_force.	Microseconds();
  timer_expand.	Microseconds();
  timer_fixp.	Microseconds();
  timer_extrn.	Microseconds();
  timer_wait.	Microseconds();
}

void ComponentContainer::initialize(void)
{
  Component *c, *c1;

				// Set centerlevl variable

  if (centerlevl < 0) centerlevl = multistep/2;
  centerlevl = min<int>(centerlevl, multistep);


  read_rates();			// Read initial processor rates


				// Look for a restart file
  unsigned char ir = 0;
  if (myid==0) {
    string resfile = outdir + infile;
    ifstream in(resfile.c_str());
    if (in) {
      cerr << "ComponentContainer::initialize: successfully opened <"
	   << resfile << ">, assuming a restart" << endl;
      ir = 1;
    } else {
      cerr << "ComponentContainer::initialize: could not open <"
	   << resfile << ">, assuming a new run" << endl;
      ir = 0;
    }
  }

  MPI_Bcast(&ir, 1, MPI_UNSIGNED_CHAR, 0, MPI_COMM_WORLD);
  restart = ir ? true : false;

  spair data;

  if (restart) {

    struct MasterHeader master;
    ifstream *in = NULL;

				// Open file
    if (myid==0) {

      string resfile = outdir + infile;
      in = new ifstream(resfile.c_str());
      if (!*in) {
	cerr << "ComponentContainer::initialize: could not open <"
	     << resfile << ">\n";
	MPI_Abort(MPI_COMM_WORLD, 5);
	exit(0);

      }

      in->read((char *)&master, sizeof(MasterHeader));
      if (!*in) {
	cerr << "ComponentContainer::initialize: "
	     << "could not read master header from <"
	     << resfile << ">\n";
	MPI_Abort(MPI_COMM_WORLD, 6);
	exit(0);
      }

      cout << "Recovering from <" << resfile << ">:"
	   << "  Tnow="  << master.time
	   << "  Ntot="  << master.ntot
	   << "  Ncomp=" << master.ncomp << endl;

      tnow  = master.time;
      ntot  = master.ntot;
      ncomp = master.ncomp;
    }

    MPI_Bcast(&tnow,  1, MPI_DOUBLE, 0, MPI_COMM_WORLD);
    
    MPI_Bcast(&ntot,  1, MPI_INT,    0, MPI_COMM_WORLD);
      
    MPI_Bcast(&ncomp, 1, MPI_INT,    0, MPI_COMM_WORLD);
      
    for (int i=0; i<ncomp; i++) {
      c = new Component(in);
      components.push_back(c);
    }
      
    delete in;
  }
  else {
    
    parse->find_list("components");

    ncomp = 0;

    while (parse->get_next(data)) {
      
      string name = trimLeft(trimRight(data.first));
      string id, pfile, cparam, fparam;
      const int linesize = 2048;
      char line[linesize];

      if (myid==0) {
	ifstream desc(data.second.c_str());
	if (!desc) {
	  cerr << "ComponentContainer::initialize: could not open ps description file <"
	       << data.second << ">\n";
	  MPI_Abort(MPI_COMM_WORLD, 6);
	  exit(0);
	}
	
	desc.get(line, linesize, '\0');

				// Replace delimeters with spaces
	for (int i=0; i<linesize; i++) {
	  if ( line[i] == '\0' ) break;
	  if ( line[i] == '\n' ) line[i] = ' ';
	}

      }
	
      MPI_Bcast(line, linesize, MPI_CHAR, 0, MPI_COMM_WORLD);

      string sline(line);
      StringTok<string> tokens(sline);
      id     = tokens(":");
      cparam = tokens(":");
      pfile  = tokens(":");
      fparam = tokens(":");

      id     = trimLeft(trimRight(id));
      cparam = trimLeft(trimRight(cparam));
      pfile  = trimLeft(trimRight(pfile));
      fparam = trimLeft(trimRight(fparam));

      c = new Component(name, id, cparam, pfile, fparam);
      components.push_back(c);
      
      ncomp++;
    }

  }

				// Initialize components
  list<Component*>::iterator cc, cc1;

  // Do this in the component constructors
  /*
  for (cc=components.begin(); cc != components.end(); cc++) {
    c = *cc;
    c->initialize();
  }
  */

				// Initialize interactions between components
  string value;
  ntot = 0;
  
  for (cc=components.begin(); cc != components.end(); cc++) {
    c = *cc;
				// Use this loop, BTW, to sum up all bodies
    ntot += c->nbodies_tot;
    
				// A new interaction list for THIS component
    Interaction *curr = new Interaction;
    curr->c = &(*c);
				// Loop through looking for pairs, it's n^2
				// but there will not be that many . . .
    parse->find_list("interaction");
    while (parse->get_next(data)) {

				// Are we talking about THIS component?
      if (c->name.compare(data.first) == 0) {
	
	for (cc1=comp.components.begin(); cc1 != comp.components.end(); cc1++) {
	  c1 = *cc1;

				// If the second in the pair matches, use it
	  if (c1->name.compare(data.second) == 0) {
	    curr->l.push_back(&(*c1));
	  }

	}
      }
    }

    if (!curr->l.empty()) 
      interaction.push_back(curr);
    else
      delete curr;

  }

  if (myid==0 && !interaction.empty()) {
    cout << "\nUsing the following component interation list:\n";
    cout << setiosflags(ios::left)
	 << setw(30) << setfill('-') << "-"
	 << "-----------" 
	 << resetiosflags(ios::left)
	 << setw(30) << setfill('-') << "-"
	 << "\n" << setfill(' ');
    
    list<Interaction*>::iterator i;
    list<Component*>::iterator j;
    for (i=interaction.begin(); i != interaction.end(); i++) {
      Interaction *inter = *i;
      for (j=inter->l.begin(); j != inter->l.end(); j++) {
	Component *comp = *j;
	cout << setiosflags(ios::left)
	     << setw(30) << inter->c->name 
	     << "acts on" 
	     << resetiosflags(ios::left)
	     << setw(30) << comp->name
	     << "\n";
      }
      cout << setiosflags(ios::left)
	   << setw(30) << setfill('-') << "-"
	   << "-----------" 
	   << resetiosflags(ios::left)
	   << setw(30) << setfill('-') << "-"
	   << "\n" << setfill(' ');
    }
    cout << "\n";
  }

  
}


ComponentContainer::~ComponentContainer(void)
{
  Component   *p1;
  Interaction *p2;

  list<Component*>::iterator c;
  for (c=comp.components.begin(); c != comp.components.end(); c++) {
    p1 = *c;
#ifdef DEBUG
    cout << "Process " << myid 
	 << " deleting component <" << p1->name << ">" << endl;
#endif
    delete p1;
  }

  list<Interaction*>::iterator i;
  for (i=interaction.begin(); i != interaction.end(); i++) {
    p2 = *i;
#ifdef DEBUG
    cout << "Process " << myid 
	 << " deleting interaction <" << p2->c->name << ">" << endl;
#endif
    delete p2;
  }

  delete [] gcom1;
  delete [] gcov1;
}

void ComponentContainer::compute_potential(unsigned mlevel)
{
  list<Component*>::iterator cc;
  Component *c;
  
#ifdef DEBUG
  cout << "Process " << myid << ": entered <compute_potential>\n";
#endif

#ifdef USE_GPTL
  GPTLstart("ComponentContainer::compute_potential");
#endif

  // Turn on step timers or VERBOSE level 4 or greater
  //
  if (VERBOSE>3) timing        = true;
  if (VERBOSE>4) thread_timing = true;

  if (timing) {
    timer_clock.start();
    timer_force.start();
    if (levcnt.size()==0) levcnt = vector<unsigned>(multistep+1, 0);
    levcnt[mlevel]++;
  }

  // Potential/force clock
  //
  for (cc=comp.components.begin(); cc != comp.components.end(); cc++)
    (*cc)->time_so_far.reset();

  //
  // Compute accel for each component
  //
  int nbeg, nend, indx;
  unsigned ntot;

  state = SELF;

  if (timing) timer_wait.start();
#ifdef USE_GPTL
#ifdef GPTL_WAIT
  GPTLstart("ComponentContainer::waiting_acceleration");
  MPI_Barrier(MPI_COMM_WORLD);
  GPTLstop ("ComponentContainer::waiting_acceleration");
#endif
  GPTLstart("ComponentContainer::acceleration");
#endif

  for (cc=comp.components.begin(); cc != comp.components.end(); cc++) {
    c = *cc;

    if (timing) {
      timer_wait.stop();
      timer_zero.start();
    }

    for (int lev=mlevel; lev<=multistep; lev++) {
      
      ntot = c->levlist[lev].size();
      
      for (unsigned n=0; n<ntot; n++) {
				// Particle index
	indx = c->levlist[lev][n];
				// Zero-out external potential
	c->Part(indx)->potext = 0.0;
				// Zero-out potential and acceleration
	c->Part(indx)->pot = 0.0;
	for (int k=0; k<c->dim; k++) c->Part(indx)->acc[k] = 0.0;
      }
    }
    if (timing) {
      timer_zero.stop();
      timer_wait.start();
    }

				// Compute new accelerations and potential
#ifdef DEBUG
    cout << "Process " << myid << ": about to call force <"
	 << c->id << "> for mlevel=" << mlevel << endl;
#endif
    if (timing) {
      timer_wait.stop();
      timer_accel.start();
    }
    c->time_so_far.start();
    c->force->set_multistep_level(mlevel);
    c->force->get_acceleration_and_potential(c);
    c->time_so_far.stop();
    if (timing) {
      timer_accel.stop();
      timer_wait.start();
    }
#ifdef DEBUG
    cout << "Process " << myid << ": force <"
	 << c->id << "> for mlevel=" << mlevel << "done" << endl;
#endif
  }

#ifdef USE_GPTL
  GPTLstop ("ComponentContainer::acceleration");
#ifdef GPTL_WAIT
  GPTLstart("ComponentContainer::waiting_interactions");
  MPI_Barrier(MPI_COMM_WORLD);
  GPTLstop ("ComponentContainer::waiting_interactions");
#endif
  GPTLstart("ComponentContainer::interactions");
#endif


  //
  // Do the component interactions
  //
  list<Interaction*>::iterator inter;
  list<Component*>::iterator other;
  vector< pair<string, Timer> >::iterator itmr;
  
  state = INTERACTION;

  if (timing) {			// Initialize interaction timers?
				// [One for each pair in the list]
				//
    unsigned npairs = 0;	// Count the pairs
    for (inter=interaction.begin(); inter != interaction.end(); inter++) {
      for (other=(*inter)->l.begin(); other != (*inter)->l.end(); other++) {
	npairs++;
      }
    }

				// Remake the timer list?
    if (npairs != timer_sntr.size()) {
      timer_sntr.clear();	// Clear the list and make a new one
      for (inter=interaction.begin(); inter != interaction.end(); inter++) {
	for (other=(*inter)->l.begin(); other != (*inter)->l.end(); other++) {
	  ostringstream sout;
	  sout << (*inter)->c->name << " <=> " << (*other)->name;
	  timer_sntr.push_back( pair<string, Timer>(sout.str(), Timer(true)) );
	}
      }
    }
    
    timer_inter.start();
    itmr = timer_sntr.begin();
  }

  for (inter=interaction.begin(); inter != interaction.end(); inter++) {
    for (other=(*inter)->l.begin(); other != (*inter)->l.end(); other++) {

#ifdef USE_GPTL
      ostringstream sout;
      sout <<"ComponentContainer::interation run<"
	   << (*inter)->c->name << "-->" << (*other)->name << ">";
      GPTLstart(sout.str().c_str());
#endif

      if (timing) {
	timer_accel.start();
	itmr->second.start();
      }
      (*other)->time_so_far.start();
      (*inter)->c->force->SetExternal();
      (*inter)->c->force->set_multistep_level(mlevel);
      (*inter)->c->force->get_acceleration_and_potential(*other);
      (*inter)->c->force->ClearExternal();
      (*other)->time_so_far.stop();
      if (timing) {
	timer_accel.stop();
	itmr->second.stop();
	itmr++;
      }

#ifdef USE_GPTL
      GPTLstop (sout.str().c_str());
#ifdef GPTL_WAIT
      sout.str("");
      sout <<"ComponentContainer::interation wait<"
	   << (*inter)->c->name << "-->" << (*other)->name << ">";
      GPTLstart(sout.str().c_str());
      MPI_Barrier(MPI_COMM_WORLD);
      GPTLstop (sout.str().c_str());
#endif
#endif      
    }
  }

  if (timing) timer_inter.stop();
      
#ifdef USE_GPTL
  GPTLstop ("ComponentContainer::interactions");
#ifdef GPTL_WAIT
  GPTLstart("ComponentContainer::waiting_external");
  MPI_Barrier(MPI_COMM_WORLD);
  GPTLstop ("ComponentContainer::waiting_external");
#endif
  GPTLstart("ComponentContainer::external");
#endif

  //
  // Do the external forces (if there are any . . .)
  //

  state = EXTERNAL;

  if (timing) {
    timer_extrn.start();
				// Initialize external force timers?
				// [One for each in external force list]
    if (external.force_list.size() != timer_sext.size()) {
      timer_sext.clear();	// Clear the list
      list<ExternalForce *>::iterator ext;
      for (ext =  external.force_list.begin(); 
	   ext != external.force_list.end(); ext++) {
	timer_sext.push_back( pair<string, Timer>((*ext)->id, Timer(true)) );
      }
    }
  }
  if (!external.force_list.empty()) {
    
    unsigned cnt=0;
    list<ExternalForce*>::iterator ext;

    for (cc=comp.components.begin(); cc != comp.components.end(); cc++) {
      c = *cc;
      c->time_so_far.start();
      if (timing) itmr = timer_sext.begin();
      for (ext=external.force_list.begin(); 
	   ext != external.force_list.end(); ext++) {
	if (timing) itmr->second.start();
	(*ext)->set_multistep_level(mlevel);
	(*ext)->get_acceleration_and_potential(c);
	if (timing) (itmr++)->second.stop();
      }
      c->time_so_far.stop();
    }

  }

#ifdef USE_GPTL
  GPTLstop ("ComponentContainer::external");
#ifdef GPTL_WAIT
  GPTLstart("ComponentContainer::waiting_centering");
  MPI_Barrier(MPI_COMM_WORLD);
  GPTLstop ("ComponentContainer::waiting_centering");
#endif
  GPTLstart("ComponentContainer::centering");
#endif



  if (timing) timer_extrn.stop();

  if (timing) timer_force.stop();
  

  state = NONE;

  //
  // Compute new center(s)
  //
  if (mactive[mstep][centerlevl]) {

    if (timing) timer_posn.start();
    fix_positions();
    if (timing) timer_posn.stop();

#ifdef DEBUG
    cout << "Process " << myid << ": returned from <fix_positions>\n";
#endif

    //
    // Recompute global com
    //
    if (timing) timer_gcom.start();
    for (int k=0; k<3; k++) gcom[k] = 0.0;
    for (cc=comp.components.begin(); cc != comp.components.end(); cc++) {
      c = *cc;
      for (int k=0; k<3; k++) gcom[k] += c->com[k];
    }
    if (timing) timer_gcom.stop();
    
#ifdef DEBUG
    cout << "Process " << myid << ": gcom computed\n";
#endif

    //
    // Compute angular momentum for each component
    //
    if (timing) timer_angmom.start();
    for (cc=comp.components.begin(); cc != comp.components.end(); cc++) {
      (*cc)->get_angmom();
    }
    if (timing) timer_angmom.stop();
    
#ifdef DEBUG
    cout << "Process " << myid << ": angmom computed\n";
#endif
    
    
    //
    // Update center of mass system coordinates
    //
    if (timing) timer_gcom.start();
    for (cc=comp.components.begin(); cc != comp.components.end(); cc++) {
      c = *cc;
      if (c->com_system) c->update_accel();
    }
    if (timing) timer_gcom.stop();
  }

#ifdef USE_GPTL
  GPTLstop ("ComponentContainer::centering");
#ifdef GPTL_WAIT
  GPTLstart("ComponentContainer::waiting_timing");
  MPI_Barrier(MPI_COMM_WORLD);
  GPTLstop ("ComponentContainer::waiting_timing");
#endif
  GPTLstart("ComponentContainer::timing");
#endif

  if (timing && timer_clock.getTime().getRealTime()>tinterval) {
    if (myid==0) {
      vector< pair<string, Timer> >::iterator itmr;
      ostringstream sout;
      sout << "--- Timer info in comp, mlevel=" << mlevel;
      cout << endl
	   << setw(70) << setfill('-') << '-' << endl
	   << setw(70) << left << sout.str().c_str() << endl
	   << setw(70) << setfill('-') << '-' << endl << setfill(' ') << right;
      
      if (multistep) {
	cout << setw(20) << "COM: "
	     << setw(18) << timer_gcom.getTime()() << endl
	     << setw(20) << "Position: "
	     << setw(18) << timer_posn.getTime()() << endl
	     << setw(20) << "" << setw(50) << setfill('-') << '-' << endl 
	     << setfill(' ') << right
	     << setw(20) << "*** " << setw(30) << left << "fix pos" << ": " 
	     << setw(18) << timer_fixp.getTime()() << endl
	     << setw(20) << "" << setw(50) << setfill('-') << '-' << endl 
	     << setfill(' ') << right
	     << setw(20) << "Ang mom: "
	     << setw(18) << timer_angmom.getTime()() << endl
	     << setw(20) << "Zero: "
	     << setw(18) << timer_zero.getTime()() << endl
	     << setw(20) << "Accel: "
	     << setw(18) << timer_accel.getTime()() << endl;

	if (thread_timing)
	  cout << setw(20) << "" << setw(50) << setfill('-') << '-' << endl 
	       << setfill(' ') << right
	       << setw(20) << "*** " << setw(30) << left << "threaded" << ": " 
	       << right << setw(18) 
	       << timer_thr_acc.getTime()() << endl
	       << setw(20) << "" << setw(50) << setfill('-') << '-' << endl 
	       << setfill(' ') << right;

	cout << setw(20) << "Interaction: "
	     << setw(18) << timer_inter.getTime()() << endl;

	if (timer_sntr.size()) {
	  cout << setw(20) << "" << setw(50) << setfill('-') << '-' << endl 
	       << setfill(' ') << right;
	  for (itmr=timer_sntr.begin(); itmr != timer_sntr.end(); itmr++) {
	    cout << setw(20) << "*** " << setw(30) << left << itmr->first 
		 << ": " << right
		 << setw(18) << itmr->second.getTime()()
		 << endl;
	  }
	  cout << setw(20) << "" << setw(50) << setfill('-') << '-' << endl 
	       << setfill(' ') << right;
	}

	if (thread_timing)
	  cout << setw(20) << "" << setw(50) << setfill('-') << '-' << endl 
	       << setfill(' ') << right
	       << setw(20) << "*** " << setw(30) << left << "threaded" << ": "
	       << right << setw(18) 
	       << timer_thr_int.getTime()() << endl
	       << setw(20) << "" << setw(50) << setfill('-') << '-' << endl 
	       << setfill(' ') << right;

	cout << setw(20) << "External: "
	     << setw(18) << timer_extrn.getTime()() << endl;

	if (thread_timing)
	  cout << setw(20) << "" << setw(50) << setfill('-') << '-' << endl 
	       << setfill(' ') << right
	       << setw(20) << "*** " << setw(30) << left << "threaded" << ": " 
	       << right << setw(18) 
	       << timer_thr_ext.getTime()() << endl
	       << setw(20) << "" << setw(50) << setfill('-') << '-' << endl 
	       << setfill(' ') << right;

	
	if (timer_sext.size()) {
	  cout << setw(20) << "" << setw(50) << setfill('-') << '-' << endl 
	       << setfill(' ') << right;
	  for (itmr = timer_sext.begin(); itmr != timer_sext.end(); itmr++) {
	    cout << setw(20) << "*** " << setw(30) << left << itmr->first 
		 << ": " << right
		 << setw(18) << itmr->second.getTime()()
		 << endl;
	  }
	  cout << setw(20) << "" << setw(50) << setfill('-') << '-' << endl 
	       << setfill(' ') << right;
	}
	  
	cout << setw(20) << "Expand: "
	     << setw(18) << timer_expand.getTime()() << endl;

	cout << setw(20) << "Force: "
	     << setw(18) << timer_force.getTime()() << endl;
      }

      cout << setw(70) << setfill('-') << '-' << endl << setfill(' ');
      cout << endl << "mstep/Mstep=" << mstep << "/" << Mstep << endl;
      unsigned n = 0;
      while (n<=multistep) {
	for (int i=0; i<5; i++) {
	  if (n<=multistep) {
	    cout << left << setw(3) << n << "|" 
		 << setw(8) << levcnt[n];
	    levcnt[n++] = 0;
	  }
	}
	cout << endl;
      }
      cout << setw(70) << setfill('-') << '-' << endl << setfill(' ');
    }

    timer_gcom.reset();
    timer_posn.reset();
    timer_fixp.reset();
    timer_angmom.reset();
    timer_zero.reset();

    timer_accel.reset();

    timer_thr_acc.reset();
    timer_thr_int.reset();
    timer_thr_ext.reset();

    timer_inter.reset();
    timer_extrn.reset();
    timer_force.reset();
    timer_expand.reset();

    timer_clock.reset();

    vector< pair<string, Timer> >::iterator itmr;

    for (itmr=timer_sntr.begin(); itmr != timer_sntr.end(); itmr++) 
      itmr->second.reset();

    for (itmr=timer_sext.begin(); itmr != timer_sext.end(); itmr++) 
      itmr->second.reset();

  }

#ifdef USE_GPTL
  GPTLstop("ComponentContainer::timing");
  GPTLstop("ComponentContainer::compute_potential");
#endif

  gottapot = true;
}


void ComponentContainer::compute_expansion(unsigned mlevel)
{
  list<Component*>::iterator cc;
  Component *c;
  
#ifdef USE_GPTL
  GPTLstart("ComponentContainer::compute_expansion");
#endif

  if (timing) timer_expand.start();

#ifdef DEBUG
  cout << "Process " << myid << ": entered <compute_expansion>\n";
#endif

  //
  // Compute expansion for each component
  //
  for (cc=comp.components.begin(); cc != comp.components.end(); cc++) {
    c = *cc;
    
#ifdef DEBUG
    cout << "Process " << myid << ": about to compute coefficients <"
	 << c->id << "> for mlevel=" << mlevel << endl;
#endif
				// Compute coefficients
    c->force->set_multistep_level(mlevel);
    c->force->determine_coefficients(c);
#ifdef DEBUG
    cout << "Process " << myid << ": coefficients <"
	 << c->id << "> for mlevel=" << mlevel << " done" << endl;
#endif
  }

#ifdef USE_GPTL
  GPTLstop("ComponentContainer::compute_expansion");
#endif

  if (timing) timer_expand.stop();
}


void ComponentContainer::multistep_reset()
{
  //
  // Do reset for each component
  //
  list<Component*>::iterator cc;
  for (cc=comp.components.begin(); cc != comp.components.end(); cc++) {
    (*cc)->force->multistep_reset();
  }
}


void ComponentContainer::print_level_list_header()
{
  ostringstream ofil;
  ofil << runtag << ".levels";

  ifstream in(ofil.str().c_str());
  if (!in) {
    in.close();
    ofstream out(ofil.str().c_str());
    out << setw(80) << setfill('-') << '-' << endl << setfill(' ')
	<< "--- Column explanations" << endl
	<< setw(80) << setfill('-') << '-' << endl << setfill(' ');
    out << left
	<< setw(15) << "L"       << ": level" << endl
	<< setw(15) << "Number"  << ": number of particles on L" << endl
	<< setw(15) << "dN/dL"   << ": fractional occupation on L" << endl
	<< setw(15) << "N(<=L)"  << ": cumulative occupation on L" << endl
	<< setw(15) << "s"       << ": per particle scale" << endl
	<< setw(15) << "v"       << ": per particle velocity" << endl
	<< setw(15) << "a"       << ": per particle acceleration" << endl
	<< setw(15) << "int"     << ": internal time step (e.g. cooling)" 
	<< endl;
    
    if (DTold)
      out << left << setw(15) << "r" 
	  << ": coordinate radius" << endl
	  << setw(15) << "f(r/v)"
	  << ": fraction with dt=|r|/|v|" << endl
	  << setw(15) << "f(s/v)" 
	  << ": fraction with dt=s/|v|" << endl
	  << setw(15) << "f(v/a)" 
	  << ": fraction with dt=|v|/|a|" << endl
	  << setw(15) << "f(r/a)" 
	  << ": fraction with dt=sqrt(|r|/|a|)" << endl
	  << setw(15) << "f(ext)"
	  << ": fraction with dt=dt(internal)" << endl;
    else
      out << left  << setw(15) << "r" 
	  << ": grav. potential scale length, |phi|/|d(phi)/dx|"  << endl
	  << setw(15) << "f(r/v)"
	  << ": fraction with dt=|phi|/|d(phi)/dx * v|" << endl
	  << setw(15) << "f(s/v)" 
	  << ": fraction with dt=s/|v|" << endl
	  << setw(15) << "f(v/a)" 
	  << ": fraction with dt=|v|/|a|" << endl
	  << setw(15) << "f(r/a)" 
	  << ": fraction with dt=sqrt(|phi|/|a*a|)" << endl
	  << setw(15) << "f(ext)"
	  << ": fraction with dt=dt(internal)" << endl;
    
    out << endl
	<< "NB: simple particles, such as stars or dark matter, will have not" 
	<< endl << "have internal length scales or time steps" << endl << endl;
    
  }

}


void ComponentContainer::print_level_lists(double T)
{
  static bool firstime = true;
  if (firstime) {
    print_level_list_header();
    firstime = false;
  }

  //
  // Do reset for each component
  //
  list<Component*>::iterator cc;
  for (cc=comp.components.begin(); cc != comp.components.end(); cc++) {
    (*cc)->print_level_lists(T);
  }
}


void ComponentContainer::multistep_debug()
{
  list<Component*>::iterator cc;
  for (cc=comp.components.begin(); cc != comp.components.end(); cc++) {
    (*cc)->force->multistep_debug();
  }
}


void ComponentContainer::fix_acceleration(void)
{
  double axcm, aycm, azcm, mtot;
  double axcm1, aycm1, azcm1, mtot1;

  axcm = aycm = azcm = mtot = 0.0;
  axcm1 = aycm1 = azcm1 = mtot1 = 0.0;

  list<Component*>::iterator cc;
  Component *c;
  PartMapItr p, pend;

  for (cc=comp.components.begin(); cc != comp.components.end(); cc++) {
    c = *cc;

    pend = c->particles.end();
    for (p=c->particles.begin(); p != pend; p++) {
    
      if (c->freeze(p->first)) continue;

      mtot1 += p->second.mass;
      axcm1 += p->second.mass*p->second.acc[0];
      aycm1 += p->second.mass*p->second.acc[1];
      azcm1 += p->second.mass*p->second.acc[2];
    }
  }

  MPI_Allreduce(&mtot1, &mtot, 1, MPI_DOUBLE, MPI_SUM, MPI_COMM_WORLD);
  MPI_Allreduce(&axcm1, &axcm, 1, MPI_DOUBLE, MPI_SUM, MPI_COMM_WORLD);
  MPI_Allreduce(&aycm1, &aycm, 1, MPI_DOUBLE, MPI_SUM, MPI_COMM_WORLD);
  MPI_Allreduce(&azcm1, &azcm, 1, MPI_DOUBLE, MPI_SUM, MPI_COMM_WORLD);

  if (mtot>0.0) {
    axcm = axcm/mtot;
    aycm = aycm/mtot;
    azcm = azcm/mtot;
  }

  for (cc=comp.components.begin(); cc != comp.components.end(); cc++) {
    c = *cc;

    pend = c->particles.end();
    for (p=c->particles.begin(); p != pend; p++) {

      if (c->freeze(p->first)) continue;
      p->second.acc[0] -= axcm;
      p->second.acc[1] -= aycm;
      p->second.acc[2] -= azcm;

    }

  }

}



void ComponentContainer::fix_positions()
{
  double mtot1, mtot0;
  MPI_Status status;

  mtot = mtot1 = 0.0;
  for (int k=0; k<3; k++) 
    gcom[k] = gcom1[k] = gcov[k] = gcov1[k] = 0.0;

  list<Component*>::iterator cc;
  Component *c;
  PartMapItr p, pend;

  for (cc=comp.components.begin(); cc != comp.components.end(); cc++) {
    c = *cc;

    if (timing) timer_fixp.start();
    c->fix_positions();
    if (timing) timer_fixp.stop();
    
    mtot1 += c->mtot;
    for (int k=0; k<3; k++) gcom1[k] += c->com[k];
    for (int k=0; k<3; k++) gcov1[k] += c->cov[k];

    if (c->EJ && (gottapot || restart)) {
      c->orient->accumulate(tnow, c);
      c->orient->logEntry  (tnow, c);
    }
    
  }

  MPI_Barrier(MPI_COMM_WORLD);
  mtot0 = 0.0;
  MPI_Allreduce(&mtot1, &mtot0, 1, MPI_DOUBLE, MPI_SUM, MPI_COMM_WORLD);
  mtot = mtot0;
  MPI_Allreduce(gcom1, gcom, 3, MPI_DOUBLE, MPI_SUM, MPI_COMM_WORLD);
  MPI_Allreduce(gcov1, gcov, 3, MPI_DOUBLE, MPI_SUM, MPI_COMM_WORLD);

  if (global_cov) {

    for (cc=comp.components.begin(); cc != comp.components.end(); cc++) {
      c = *cc;

      pend = c->particles.end();
      for (p=c->particles.begin(); p != pend; p++) {
    
	if (c->freeze(p->first)) continue;

	for (int k=0; k<3; k++) p->second.vel[k] -= gcov[k];
      }
    }
  }

}


void ComponentContainer::read_rates(void)
{

  rates =  vector<double>(numprocs);

  if (myid == 0) {
    ifstream in(ratefile.c_str());
    
    double norm = 0.0;
    
    if (in) {			// We are reading from a file
      for (int n=0; n<numprocs; n++) {
	in >> rates[n];
	if (!in) {
	  cerr << "setup: error reading <" << ratefile << ">\n";
	  MPI_Abort(MPI_COMM_WORLD, 33);
	  exit(0);
	}
	norm += rates[n];
      }

    } else {			// Assign equal rates to all nodes
      cerr << "setup: can not find <" << ratefile << "> . . . will assume homogeneous cluster\n";
      for (int n=0; n<numprocs; n++) {
	rates[n] = 1.0;
	norm += rates[n];
      }
    }

    for (int n=0; n<numprocs; n++) rates[n] /= norm;

  }

  MPI_Bcast(&rates[0], numprocs, MPI_DOUBLE, 0, MPI_COMM_WORLD);

}


void ComponentContainer::report_numbers(void)
{
  if (!nreport || this_step % nreport)  return;

  for (int num=0; num<numprocs; num++) {

    if (myid==num) {
      string fout = outdir + runtag + ".number";
      ofstream out(fout.c_str(), ios::out | ios::app);
      list<Component*>::iterator cc;
      if (out) {
	if (myid==0) {
	  out << "# Step: " << this_step << " Time: " << tnow << endl 
	      << right << "# " << setw(5)  << "Proc";
	  for (cc=comp.components.begin(); 
	       cc != comp.components.end(); cc++) {
	    out << setw(20) << (*cc)->name << setw(20) << "Effort";
	  }
	  out << endl << "# " << setw(5) << "-----";
	  for (cc=comp.components.begin(); 
	       cc != comp.components.end(); cc++) {
	    out << setw(20) << "----------" << setw(20) << "----------";
	  }
	  out << endl;
	}
	out << setw(7) << num;
	for (cc=comp.components.begin(); cc != comp.components.end(); cc++) {
	  out << setw(20) << (*cc)->Number();
	  double toteff = 0.0;
	  PartMap::iterator tp;
	  for (tp=(*cc)->particles.begin(); tp!=(*cc)->particles.end(); tp++)
	    toteff += tp->second.effort;
	  out << setw(20) << toteff;
	}
	out << endl;
      }
    }
    MPI_Barrier(MPI_COMM_WORLD);
  }
}


void ComponentContainer::load_balance(void)
{
  if (!nbalance || this_step % nbalance)  return;

				// Query timers
  vector<double> rates1(numprocs, 0.0), trates(numprocs, 0.0);
  rates1[myid] = MPL_read_timer(1);
  MPI_Allreduce(&rates1[0], &trates[0], numprocs, 
		MPI_DOUBLE, MPI_SUM, MPI_COMM_WORLD);

				// Compute normalized rate vector
  double norm = 0.0;
  for (int i=0; i<numprocs; i++) {
    rates1[i] = 1.0/trates[i];
    norm += rates1[i];
  }
  for (int i=0; i<numprocs; i++) rates1[i] /= norm;

				// For debugging
#ifdef RANDOMTIME
  {
    if (myid==0) cout << "*** WARNING: using random time intervals for load balance testing ***\n";
    double norm = 0.0;
    for (int i=0; i<numprocs; i++) {
      rates1[i] = rand();
      norm += rates1[i];
    }
    for (int i=0; i<numprocs; i++) rates1[i] /= norm;
  }
#endif

				// Compare relative difference with threshold
  bool toobig = false;
  double curdif;
  for (int n=0; n<numprocs; n++) {
    if (rates[n]>0.0) {
      curdif = fabs(rates[n]-rates1[n])/rates[n];
      if ( curdif > dbthresh) toobig = true;
    }
  }
  

				// Print out info
  if (myid==0) {
    
    string outrates = outdir + "current.processor.rates.test." + runtag;

    ofstream out(outrates.c_str(), ios::out | ios::app);
    if (out) {
      out << "# Step: " << this_step << endl;
      out << "# "
	  << setw(5)  << "Proc"
	  << setw(15) << "Step time"
	  << setw(15) << "Norm rate"
	  << setw(15) << "Rate frac"
	  << endl
	  << "# "
	  << setw(5)  << "-----"
	  << setw(15) << "----------"
	  << setw(15) << "----------"
	  << setw(15) << "----------"
	  << endl;
      
      for (int n=0; n<numprocs; n++) {
	out << "  "
	    << setw(5) << n
	    << setw(15) << trates[n]
	    << setw(15) << rates1[n];

	if (rates[n]>0.0)
	  out << setw(15) << fabs(rates[n]-rates1[n])/rates[n] << endl;
	else
	  out << setw(15) << " ***" << endl;
      
      }
    }
  }


  if (toobig) {

				// Use new rates
    rates = rates1;

				// Initiate load balancing for each component
    list<Component*>::iterator cc;
  
    for (cc=comp.components.begin(); cc != comp.components.end(); cc++) {
      (*cc)->load_balance();
    }

  }

}

bool ComponentContainer::bad_values()
{
  bool bad = false;
  PartMapItr it;
  list<Component*>::iterator cc;
  for (cc=comp.components.begin(); cc != comp.components.end(); cc++) {
    bool badval = false;
    for (it=(*cc)->Particles().begin(); it!=(*cc)->Particles().end(); it++) {
      if (isnan(it->second.mass)) badval=true;
      for (int k=0; k<3; k++) {
	if (isnan(it->second.pos[k]))  badval=true;
	if (isnan(it->second.vel[k]))  badval=true;
	if (isnan(it->second.acc[k]))  badval=true;
      }
      if (badval) {
	cout << "Bad value in <" << (*cc)->name << ">: ";
	cout << setw(12) << it->second.indx
	     << setw(16) << hex << it->second.key << dec
	     << setw(18) << it->second.mass;
	for (int k=0; k<3; k++)
	  cout << setw(18) << it->second.pos[k];
	for (int k=0; k<3; k++)
	  cout << setw(18) << it->second.vel[k];
	for (int k=0; k<3; k++)
	  cout << setw(18) << it->second.acc[k];
	cout << endl;
	bad = true;
	break;
      }
    }
  }
  return bad;
}
