/*
  Compute accelerations, potential, and density.
*/

#ifdef RCSID
static char rcsid[] = "$Id$";
#endif

#include "expand.h"

#include <vector>

#include <ComponentContainer.H>
#include <StringTok.H>

ComponentContainer::ComponentContainer(void)
{
				// Do nothing
}

void ComponentContainer::initialize(void)
{
  spair data;

				// Initialize individual components
  if (restart) {

    struct MasterHeader master;
    ifstream *in = NULL;

				// Open file
    if (myid==0) {

      in = new ifstream(infile.c_str());
      if (!*in) {
	cerr << "ComponentContainer::initialize: could not open <"
	     << infile << ">\n";
	MPI_Abort(MPI_COMM_WORLD, 5);
	exit(0);

      }

      in->read(&master, sizeof(MasterHeader));
      if (!*in) {
	cerr << "ComponentContainer::initialize: could not read master header from <"
	     << infile << ">\n";
	MPI_Abort(MPI_COMM_WORLD, 6);
	exit(0);
      }

      cout << "Recovering from <" << infile << ">:"
	   << "  Tnow=" << master.time
	   << "  Ntot=" << master.ntot
	   << "  Ncomp=" << master.ncomp << endl;

      tnow = master.time;
      ntot = master.ntot;
      ncomp = master.ncomp;
    }

    MPI_Bcast(&tnow, 1, MPI_DOUBLE, 0, MPI_COMM_WORLD);
    tpos = tvel = tnow;
    
    MPI_Bcast(&ntot, 1, MPI_INT, 0, MPI_COMM_WORLD);
      
    MPI_Bcast(&ncomp, 1, MPI_INT, 0, MPI_COMM_WORLD);
      
    for (int i=0; i<ncomp; i++) components.push_back(new Component(in));

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
      id = tokens(":");
      cparam = tokens(":");
      pfile = tokens(":");
      fparam = tokens(":");

      id = trimLeft(trimRight(id));
      cparam = trimLeft(trimRight(cparam));
      pfile = trimLeft(trimRight(pfile));
      fparam = trimLeft(trimRight(fparam));

      components.push_back(new Component(name, id, cparam, pfile, fparam));
      
      ncomp++;
    }

  }

				// Initialize components
  list<Component*>::iterator cc, cc1;
  Component *c, *c1;

  for (cc=components.begin(); cc != components.end(); cc++) {
    c = *cc;
    c->initialize();
  }

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

    if (!curr->l.empty()) interaction.push_back(curr);

  }

  if (myid==0 && !interaction.empty()) {
    cout << "\nUsing the following component interation list:\n";
    cout << setw(-30) << setfill('-') << "-"
	 << "------" 
	 << setw(30) << setfill('-') << "-"
	 << "\n" << setfill(' ');
    
    list<Interaction*>::iterator i;
    list<Component*>::iterator j;
    for (i=interaction.begin(); i != interaction.end(); i++) {
      Interaction *inter = *i;
      for (j=inter->l.begin(); j != inter->l.end(); j++) {
	Component *comp = *j;
	cout << setw(-30) << inter->c->name 
	     << " <==> " 
	     << setw(30) << comp->name
	     << "\n";
      }
      cout << setw(-30) << setfill('-') << "-"
	   << "------" 
	   << setw(30) << setfill('-') << "-"
	   << "\n" << setfill(' ');
    }
    cout << "\n";
  }
}


ComponentContainer::~ComponentContainer(void)
{
  list<Component*>::iterator c;
  for (c=comp.components.begin(); c != comp.components.end(); c++)
    delete *c;

  list<Interaction*>::iterator i;
  for (i=interaction.begin(); i != interaction.end(); i++)
    delete *i;
}

void ComponentContainer::compute_potential(void)
{
  list<Component*>::iterator cc;
  Component *c;
  vector<Particle>::iterator p, pend;

  for (cc=comp.components.begin(); cc != comp.components.end(); cc++) {
    c = *cc;

    pend = c->particles.end();
    for (p=c->particles.begin(); p != pend; p++) {

				// Zero-out external potential
      p->potext = 0.0;
				// Zero-out potential and acceleration
      p->pot = 0.0;
      for (int i=0; i<c->dim; i++) p->acc[i] = 0.0;
    }

				// Compute new accelerations and potential

    c->force->get_acceleration_and_potential(&(c->particles));
  
  }
      

  // Do the component interactions

  list<Interaction*>::iterator inter;
  list<Component*>::iterator other;
  
  for (inter=interaction.begin(); inter != interaction.end(); inter++) {
    for (other=(*inter)->l.begin(); other != (*inter)->l.end(); other++) {

      (*inter)->c->force->SetExternal();
      (*inter)->c->force->get_acceleration_and_potential(&((*other)->particles));
      (*inter)->c->force->ClearExternal();

    }
  }
      
  // Do the external forces (if there are any . . .)

  if (!external.force_list.empty()) {

    list<ExternalForce*>::iterator ext;

    for (cc=comp.components.begin(); cc != comp.components.end(); cc++) {
      c = *cc;
      for (ext=external.force_list.begin(); 
	   ext != external.force_list.end(); ext++) {
	(*ext)->get_acceleration_and_potential(&(c->particles));
      }
    }

  }
  
  // Recompute global com
  for (int k=0; k<3; k++) gcom[k] = 0.0;
  for (cc=comp.components.begin(); cc != comp.components.end(); cc++) {
    c = *cc;
    for (int k=0; k<3; k++) gcom[k] += c->com[k];
  }

  fix_positions();

  if (fixacc) fix_acceleration();

}


void ComponentContainer::fix_acceleration(void)
{
  double axcm, aycm, azcm, mtot;
  double axcm1, aycm1, azcm1, mtot1;

  axcm = aycm = azcm = mtot = 0.0;
  axcm1 = aycm1 = azcm1 = mtot1 = 0.0;

  list<Component*>::iterator cc;
  Component *c;
  vector<Particle>::iterator p, pend;

  for (cc=comp.components.begin(); cc != comp.components.end(); cc++) {
    c = *cc;

    pend = c->particles.end();
    for (p=c->particles.begin(); p != pend; p++) {
    
      if (p->freeze()) continue;

      mtot1 += p->mass;
      axcm1 += p->mass*p->acc[0];
      aycm1 += p->mass*p->acc[1];
      azcm1 += p->mass*p->acc[2];
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

      if (p->freeze()) continue;
      p->acc[0] -= axcm;
      p->acc[1] -= aycm;
      p->acc[2] -= azcm;

    }

  }

}



void ComponentContainer::fix_positions(void)
{
  double mtot1, mtot0;
  double *gcom1 = new double [3];
  double *gcov1 = new double [3];
  MPI_Status status;

  mtot = mtot1 = 0.0;
  for (int k=0; k<3; k++) gcom[k] = gcom1[k] = gcov[k] = gcov1[k] = 0.0;

  list<Component*>::iterator cc;
  Component *c;
  vector<Particle>::iterator p, pend;

  for (cc=comp.components.begin(); cc != comp.components.end(); cc++) {
    c = *cc;

    c->fix_positions();
    
    mtot1 += c->mtot;
    for (int k=0; k<3; k++) gcom1[k] += c->com[k];
    for (int k=0; k<3; k++) gcov1[k] += c->cov[k];

    // cout << "Process " << myid << ": about to accumulate center . . . ";
    // MPI_Barrier(MPI_COMM_WORLD); // I'm not sure why I need this here . . .
    if (c->EJ) c->orient->accumulate(&c->particles, &c->center[0]);
    // cout << "Process " << myid << ": done\n";
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
    
	if (p->freeze()) continue;

	for (int k=0; k<3; k++) p->vel[k] -= gcov[k];
      }
    }
  }

  delete [] gcom1;
  delete [] gcov1;

}


