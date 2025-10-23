#include <ios>
#include <yaml-cpp/yaml.h>	      // YAML support
#include "Sutils.H"		      // For trim-copy

#include "PSP.H"
#include "libvars.H"		// Library support

bool badstatus(std::istream& in)
{
  std::ios::iostate i = in.rdstate();
  
  if (i & std::ios::eofbit) {
    std::cout << "EOF encountered" << std::endl;
    return true;
  }
  else if(i & std::ios::failbit) {
    std::cout << "Non-Fatal I/O error" << std::endl;;
    return true;
  }  
  else if(i & std::ios::badbit) {
    std::cout << "Fatal I/O error" << std::endl;
    return true;
  }
  else
    return false;
}


PSPout::PSPout(const std::string& infile, bool verbose) : PSP(verbose, "")
{
  // Open the file
  // -------------
  try {
    in.open(infile);
  }
  catch(const std::runtime_error& err) {
    std::ostringstream sout;
    sout << "Could not open PSP file <" << infile << ">"
	 << " Error is: " << err.what();
    throw std::runtime_error(sout.str());
  }
  catch (...) {
    std::ostringstream sout;
    sout << "Could not open PSP file <" << infile << "> Unknown error";
    throw std::runtime_error(sout.str());
  }

  pos = in.tellg();

  // Read the header, quit on failure
  // --------------------------------
  try {
    in.read((char *)&header, sizeof(MasterHeader));
  }
  catch (const std::runtime_error& err) {
    std::ostringstream sout;
    sout << "Could not open master header for <" << infile << ">"
	 << " Error is: " << err.what();
    throw std::runtime_error(sout.str());
  }
  catch (...) {
    std::ostringstream sout;
    sout << "Could not read master header for <" << infile << ">"
	 << " Unknown error";
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
    }
    catch (const std::runtime_error& err) {
      std::ostringstream sout;
      sout << "Error reading magic for <" << infile << ">"
	   << " Error is: " << err.what();
      throw std::runtime_error(sout.str());
    }
    catch (...) {
      std::ostringstream sout;
      sout << "Error reading magic for <" << infile << ">";
      throw std::runtime_error(sout.str());
    }
      
    try {
      stanza.comp.read(&in);
    } catch (...) {
      std::ostringstream sout;
      sout << "Error reading component header for <" << infile << ">";
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
		<< header.time << " . . . quit reading <" << infile << ">";
      break;
    }
      
    stanzas.push_back(stanza);

    if (verbose) {
      std::cout << errors.str();
    }
  }

  spos = stanzas.begin();
}


PSPspl::PSPspl(const std::string& master, const std::string dir, bool verbose) : PSP(verbose, dir)
{
  // Open the file
  // -------------
  try {
    in.open(new_dir+master);
  }
  catch (const std::runtime_error& err) {
    std::ostringstream sout;
    sout << "Could not open the master SPL file <" << new_dir + master << ">"
	 << " Error is: " << err.what();
    throw std::runtime_error(sout.str());
  }
  catch (...) {
    std::ostringstream sout;
    sout << "Could not open the master SPL file <" << new_dir + master << ">"
	 << " Unknown error";
    throw std::runtime_error(sout.str());
  }

  if (!in.good()) {
    std::ostringstream sout;
    sout << "Error opening master SPL file <" << new_dir + master << ">";
    throw std::runtime_error(sout.str());
  }

  // Read the header, quit on failure
  // --------------------------------
  try {
    in.read((char *)&header, sizeof(MasterHeader));
  } catch (...) {
    std::ostringstream sout;
    sout << "Could not read master header for <" << master << ">";
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
	   << master << ">";
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
	   << master << ">";
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
	   << " from <" << master << ">";
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


void PSP::check_dirname()
{
  if (new_dir.size()>0) {
    if (new_dir.back() != '/') new_dir += '/';
  }
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

SParticle* PSPout::GetParticle()
{
  pcount = 0;
  
  in.seekg(cur->pspos);

  return NextParticle();
}

SParticle *PSPout::NextParticle()
{
  badstatus(in);		// DEBUG
  
  // Read partcle
  // ------------
  if (pcount < spos->comp.nbod) {
    
    part.read(in, spos->r_size, pcount++, spos);
    
    return &part;
    
  } else
    return 0;
}

SParticle* PSPspl::GetParticle()
{
  pcount = 0;

  // Set iterator to beginning of vector
  fit = spos->nparts.begin();

  // Open next file in sequence
  openNextBlob();
  
  return NextParticle();
}

void PSPspl::openNextBlob()
{
  // Close current file
  //
  if (in.is_open()) in.close();

  std::string curfile(*fit);

  try {
    if (new_dir.size()) {
      auto pos = curfile.find_last_of("/");
      if (pos != std::string::npos) // Rewrite leading directory
	curfile = new_dir + curfile.substr(pos);
      else
	curfile = new_dir + curfile;
    }
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

  // Check for empty blob
  if (N==0) openNextBlob();

}

SParticle* PSPspl::NextParticle()
{
  badstatus(in);		// DEBUG
  
  // Read partcle
  // ------------
  if (pcount < spos->comp.nbod) {
    
    if (fcount==N) openNextBlob();

    // Read blob
    part.read(in, spos->r_size, pcount++, spos);
    fcount++;
    
    return &part;
    
  } else
    return 0;
}


void PSP::ComputeStats()
{
  cur = &(*spos);
  
  // Initialize lists
  vector< vector<float> > plist(3, vector<float>(spos->comp.nbod) );
  vector< vector<float> > vlist(3, vector<float>(spos->comp.nbod) );
  mtot = 0.0;
  
  SParticle *P = GetParticle();
  unsigned n=0;
  while (P) {
    if (spos->r_size == sizeof(float)) {
      mtot += P->f->mass;
      for (unsigned k=0; k<3; k++) {
	plist[k][n] = P->f->pos[k];
	vlist[k][n] = P->f->vel[k];
      }
    } else {
      mtot += P->d->mass;
      for (unsigned k=0; k<3; k++) {
	plist[k][n] = P->d->pos[k];
	vlist[k][n] = P->d->vel[k];
      }
    }
    P = NextParticle();
    n++;
  }
  
  pmin = vector<float>(3);
  pmed = vector<float>(3);
  pmax = vector<float>(3);
  
  vmin = vector<float>(3);
  vmed = vector<float>(3);
  vmax = vector<float>(3);
  
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
  for (SParticle* part=GetParticle(); part!=0; part=NextParticle()) 
    {
      part->write(out, real4, its->index_size);
      count++;
    }
  
  if (VERBOSE)
    std::cerr << std::string(72, '-') << std::endl
	      << "Wrote " << count << " particles "<< std::endl
	      << std::string(72, '-') << std::endl
	      << spos->comp << std::endl
	      << std::string(72, '-') << std::endl;
}

void SParticle::write(std::ostream& out, bool real4, size_t isiz)
{
  unsigned long k;
  double d;
  float f;
  int j;
  
  if (isiz) 
    out.write((const char *)&(k=indx()), sizeof(unsigned long));
  
  if (real4) {
    out.write((const char *)&(f=mass()), sizeof(float));
    for (int i=0; i<3; i++) 
      out.write((const char *)&(f=pos(i)), sizeof(float));
    for (int i=0; i<3; i++) 
      out.write((const char *)&(f=vel(i)), sizeof(float));
    out.write((const char *)&(f=phi()), sizeof(float));
    for (int i=0; i<niatr(); i++)
      out.write((const char *)&(j=iatr(i)), sizeof(int));
    for (int i=0; i<ndatr(); i++)
      out.write((const char *)&(f=datr(i)), sizeof(float));
  } else {
    out.write((const char *)&(d=mass()), sizeof(double));
    for (int i=0; i<3; i++) 
      out.write((const char *)&(d=pos(i)), sizeof(double));
    for (int i=0; i<3; i++) 
      out.write((const char *)&(d=vel(i)), sizeof(double));
    out.write((const char *)&(d=phi()), sizeof(double));
    for (int i=0; i<niatr(); i++) 
      out.write((const char *)&(j=iatr(i)), sizeof(int));
    for (int i=0; i<ndatr(); i++)
      out.write((const char *)&(d=datr(i)), sizeof(double));
  }
}

// PSP factory: choose type based on file name
std::shared_ptr<PSP> PSP::getPSP(const std::string& file, const std::string dir, bool verbose)
{
  if (file.find("SPL") != std::string::npos)
    return std::make_shared<PSPspl>(file, dir, verbose);
  else
    return std::make_shared<PSPout>(file, verbose);
}
