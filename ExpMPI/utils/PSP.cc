#include <PSP.H>

#include <StringTok.H>
extern string trimLeft(const string);
extern string trimRight(const string);


PSPDump::PSPDump(ifstream *in, bool tipsy)
{
  TIPSY = tipsy;

  while (1) {

    Dump dump;

    dump.pos = in->tellg();
				// Read the header, quit on failure
				// --------------------------------
    if(!in->read((char *)&dump.header, sizeof(MasterHeader))) break;


    bool ok = true;

    for (int i=0; i<dump.header.ncomp; i++) {

      Stanza stanza;
      stanza.pos = in->tellg();
      
      ComponentHeader headerC;
      if (!headerC.read(in)) {
	cerr << "Error reading component header\n";
	ok = false;
	break;
      }

      stanza.pspos = in->tellg();

				// Parse the info string
				// ---------------------
      StringTok<string> tokens(headerC.info);
      stanza.name = trimLeft(trimRight(tokens(":")));
      stanza.id = trimLeft(trimRight(tokens(":")));
      stanza.param = trimLeft(trimRight(tokens(":")));
      
				// Strip of the tipsy type
      StringTok<string> tipsytype(stanza.name);
      stanza.ttype = trimLeft(trimRight(tipsytype(" ")));
      stanza.nbod = headerC.nbod;
      stanza.niatr = headerC.niatr;
      stanza.ndatr = headerC.ndatr;

      
				// Skip forward to next header
				// ---------------------------
      in->seekg(headerC.nbod*(8*sizeof(double)             + 
			      headerC.niatr*sizeof(int)    +
			      headerC.ndatr*sizeof(double)
			      ), ios::cur);

      dump.stanzas.push_back(stanza);

      if (TIPSY) {

				// Count up Tipsy types and make
				// linked  lists
				// -----------------------------
	if (!stanza.ttype.compare("gas")) {
	  dump.ngas += stanza.nbod;
	  dump.ntot += stanza.nbod;
	  dump.gas.push_back(stanza);
	}
	if (!stanza.ttype.compare("dark")) {
	  dump.ndark += stanza.nbod;
	  dump.ntot += stanza.nbod;
	  dump.dark.push_back(stanza);
	}
	if (!stanza.ttype.compare("star")) {
	  dump.nstar += stanza.nbod;
	  dump.ntot += stanza.nbod;
	  dump.star.push_back(stanza);
	}
      }

    }
    
    if (!ok) break;

    dumps.push_back(dump);
  }

  fid = &(*dumps.begin());

}

double PSPDump::SetTime(double time) 
{
  list<Dump>::iterator it;
  double tdif = 1.0e30;

  for (it=dumps.begin(); it!=dumps.end(); it++) {

    if (fabs(time - it->header.time) < tdif) {
      fid = &(*it);
      tdif = fabs(time-it->header.time);
    }
  }
  
  return fid->header.time;
}


void PSPDump::PrintSummary(ostream &out)
{
  list<Dump>::iterator itd;
  list<Stanza>::iterator its;

  for (itd = dumps.begin(); itd != dumps.end(); itd++) {

    out << "Time=" << itd->header.time << "   [" << itd->pos << "]" << endl;
    out << "   Total particle number: " << itd->header.ntot  << endl;
    out << "   Number of components:  " << itd->header.ncomp << endl;
    if (TIPSY) {
      out << "          Gas particles:  " << itd->ngas << endl;
      out << "         Dark particles:  " << itd->ndark << endl;
      out << "         Star particles:  " << itd->nstar << endl;
    }

    int cnt=1;

    for (its = itd->stanzas.begin(); its != itd->stanzas.end(); its++) {
	
				// Print the info for this stanza
				// ------------------------------
      out << setw(60) << setfill('-') << "-" << endl << setfill(' ');
      out << "--- Component #" << setw(2) << cnt++ << endl;
      out << setw(20) << " name :: "  << its->name   << endl
	  << setw(20) << " id :: "    << its->id     << endl
	  << setw(20) << " param :: " << its->param  << endl;
      if (TIPSY) out << setw(20) << " tipsy :: " << its->ttype  << endl;
      out << setw(20) << " nbod :: "  << its->nbod  << endl
	  << setw(20) << " niatr :: " << its->niatr << endl
	  << setw(20) << " ndatr :: " << its->ndatr << endl;
      out << setw(60) << setfill('-') << "-" << endl << setfill(' ');
      
    }
  }
}

void PSPDump::PrintSummaryCurrent(ostream &out)
{
  list<Stanza>::iterator its;

  out << "Time=" << fid->header.time << "   [" << fid->pos << "]" << endl;
  out << "   Total particle number: " << fid->header.ntot  << endl;
  out << "   Number of components:  " << fid->header.ncomp << endl;
  if (TIPSY) {
    out << "          Gas particles:  " << fid->ngas << endl;
    out << "         Dark particles:  " << fid->ndark << endl;
    out << "         Star particles:  " << fid->nstar << endl;
  }

  int cnt=1;

  for (its = fid->stanzas.begin(); its != fid->stanzas.end(); its++) {
	
    // Print the info for this stanza
    // ------------------------------
    out << setw(60) << setfill('-') << "-" << endl << setfill(' ');
    out << "--- Component #" << setw(2) << cnt++ << endl;
    out << setw(20) << " name :: "  << its->name   << endl
	<< setw(20) << " id :: "    << its->id     << endl
	<< setw(20) << " param :: " << its->param  << endl;
    if (TIPSY) out << setw(20) << " tipsy :: " << its->ttype  << endl;
    out << setw(20) << " nbod :: "  << its->nbod  << endl
	<< setw(20) << " niatr :: " << its->niatr << endl
	<< setw(20) << " ndatr :: " << its->ndatr << endl;
    out << setw(60) << setfill('-') << "-" << endl << setfill(' ');
    
  }
}

Dump* PSPDump::GetDump()
{
  sdump = dumps.begin();
  fid = &(*sdump);
  if (sdump != dumps.end()) 
    return fid;
  else 
    return 0;
}  

Dump* PSPDump::NextDump()
{
  sdump++;
  fid = &(*sdump);
  if (sdump != dumps.end()) 
    return fid;
  else 
    return 0;
}


Stanza* PSPDump::GetStanza()
{
  spos = fid->stanzas.begin();
  cur = &(*spos);
  if (spos != fid->stanzas.end()) 
    return cur;
  else 
    return 0;
}

Stanza* PSPDump::NextStanza()
{
  spos++;
  cur = &(*spos);
  if (spos != fid->stanzas.end()) 
    return cur;
  else 
    return 0;
}

SParticle *PSPDump::GetParticle(istream* in)
{
				// Position to beginning of particles
  in->seekg(spos->pspos);
  pcount = 0;

				// Clear particle
  part.iatr.erase(part.iatr.begin(), part.iatr.end()); 
  if (spos->niatr) part.iatr = vector<int>(spos->niatr);

  part.datr.erase(part.datr.begin(), part.datr.end()); 
  if (spos->ndatr) part.datr = vector<double>(spos->ndatr);

  return NextParticle(in);
}

SParticle *PSPDump::NextParticle(istream* in)
{
				// Read partcle
  if (pcount < spos->nbod) {
    in->read(&part.mass, sizeof(double));
    for (int i=0; i<3; i++) in->read(&part.pos[i], sizeof(double));
    for (int i=0; i<3; i++) in->read(&part.vel[i], sizeof(double));
    in->read(&part.phi, sizeof(double));
    for (int i=0; i<spos->niatr; i++) in->read(&part.iatr[i], sizeof(int));
    for (int i=0; i<spos->ndatr; i++) in->read(&part.datr[i], sizeof(double));

    pcount++;

    return &part;

  } else
    return 0;
  
}

Stanza* PSPDump::GetGas()
{
  spos = fid->gas.begin();
  cur = &(*spos);
  if (spos != fid->gas.end()) 
    return cur;
  else 
    return 0;
}

Stanza* PSPDump::NextGas()
{
  spos++;
  cur = &(*spos);
  if (spos != fid->gas.end()) 
    return cur;
  else 
    return 0;
}

Stanza* PSPDump::GetDark()
{
  spos = fid->dark.begin();
  cur = &(*spos);
  if (spos != fid->dark.end()) 
    return cur;
  else 
    return 0;
}

Stanza* PSPDump::NextDark()
{
  spos++;
  cur = &(*spos);
  if (spos != fid->dark.end()) 
    return cur;
  else 
    return 0;
}

Stanza* PSPDump::GetStar()
{
  spos = fid->star.begin();
  cur = &(*spos);
  if (spos != fid->star.end()) 
    return cur;
  else 
    return 0;
}

Stanza* PSPDump::NextStar()
{
  spos++;
  cur = &(*spos);
  if (spos != fid->star.end()) 
    return cur;
  else 
    return 0;
}


