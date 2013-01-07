/*
  Info on each stanza of a phase space dump

  MDWeinberg 12/23/09
*/

#include <cstdlib>
#include <iostream>
#include <sstream>
#include <fstream>
#include <iomanip>
#include <vector>
#include <string>
#include <list>

#include <StringTok.H>
#include <header.H>
#include <PSP.H>

extern string trimLeft(const string);
extern string trimRight(const string);

				// Globals for exputil library
				// Unused here
int myid = 0;
char threading_on = 0;
pthread_mutex_t mem_lock;
string outdir, runtag;

				// Useful globals
int  use_pos = 0;
bool use_int = false;
int comp = 2;
const char *compnames[] = {"dark", "star", "gas"};
string suffix = "array";

//-------------
// Help message
//-------------

void Usage(char* prog) {
  cerr << prog << ": [-t time -v -h] filename\n\n";
  cerr << "    -t time         use dump closest to <time>\n";
  cerr << "    -a pos          attribute position\n";
  cerr << "    -c comp         which component (0=dark, 1=star, 2=gas [default])\n";
  cerr << "    -i              use integer rather than float component";
  cerr << "    -o name         suffix (default: array)\n";
  cerr << "    -h              print this help message\n";
  cerr << "    -v              verbose output\n\n";
  exit(0);
}

void write_array(ifstream *in, PSPDump &psp)
{
  string filename = string(compnames[comp]) + "." + suffix;
  ofstream out(filename.c_str());
  if (!out) {
    string msg = "could not open <" + filename + ">";
    throw msg.c_str();
  }

  
  PSPstanza *stanza;
  SParticle *part;

  if (comp==2) {
    
    for (stanza=psp.GetGas(); stanza!=0; stanza=psp.NextGas()) {
      
      out << setw(20) << stanza->comp.nbod << endl;

      for (part=psp.GetParticle(in); part!=0; part=psp.NextParticle(in)) {

	if (part->f.get()) {
	  if (use_int) {
	    if (use_pos<static_cast<int>(part->f->iatr.size()))
	      out << setw(20) << part->f->iatr[use_pos] << endl;
	    else {
	      ostringstream msg;
	      msg << "attribute size is " << part->f->iatr.size()
		  << " and you can not use pos="  << use_pos;
	      throw msg.str().c_str();
	    }
	  } else {
	    if (use_pos<static_cast<int>(part->f->datr.size()))
	      out << setw(20) << part->f->datr[use_pos] << endl;
	    else {
	      ostringstream msg;
	      msg << "attribute size is " << part->f->datr.size()
		  << " and you can not use pos="  << use_pos;
	      throw msg.str().c_str();
	    }
	  }
	} else {
	  if (use_int) {
	    if (use_pos<static_cast<int>(part->d->iatr.size()))
	      out << setw(20) << part->d->iatr[use_pos] << endl;
	    else {
	      ostringstream msg;
	      msg << "attribute size is " << part->d->iatr.size()
		  << " and you can not use pos="  << use_pos;
	      throw msg.str().c_str();
	    }
	  } else {
	    if (use_pos<static_cast<int>(part->d->datr.size()))
	      out << setw(20) << part->d->datr[use_pos] << endl;
	    else {
	      ostringstream msg;
	      msg << "attribute size is " << part->d->datr.size()
		  << " and you can not use pos="  << use_pos;
	      throw msg.str().c_str();
	    }
	  }
	}
      }

    }
    
  }

  if (comp==0) {

    for (stanza=psp.GetDark(); stanza!=0; stanza=psp.NextDark()) {
      
      out << setw(20) << stanza->comp.nbod << endl;

      for (part=psp.GetParticle(in); part!=0; part=psp.NextParticle(in)) {

	if (part->f.get()) {
	  if (use_int) {
	    if (use_pos<static_cast<int>(part->f->iatr.size()))
	      out << setw(20) << part->f->iatr[use_pos] << endl;
	    else {
	      ostringstream msg;
	      msg << "attribute size is " << part->f->datr.size()
		  << " and you can not use pos="  << use_pos;
	      throw msg.str().c_str();
	    }
	  } else {
	    if (use_pos<static_cast<int>(part->f->datr.size()))
	      out << setw(20) << part->f->datr[use_pos] << endl;
	    else {
	      ostringstream msg;
	      msg << "attribute size is " << part->f->datr.size()
		  << " and you can not use pos="  << use_pos;
	      throw msg.str().c_str();
	    }
	  }

	} else {

	  if (use_int) {
	    if (use_pos<static_cast<int>(part->d->iatr.size()))
	      out << setw(20) << part->d->iatr[use_pos] << endl;
	    else {
	      ostringstream msg;
	      msg << "attribute size is " << part->d->datr.size()
		  << " and you can not use pos="  << use_pos;
	      throw msg.str().c_str();
	    }
	  } else {
	    if (use_pos<static_cast<int>(part->d->datr.size()))
	      out << setw(20) << part->d->datr[use_pos] << endl;
	    else {
	      ostringstream msg;
	      msg << "attribute size is " << part->d->datr.size()
		  << " and you can not use pos="  << use_pos;
	      throw msg.str().c_str();
	    }
	  }
	}
	
      }

    }

  }
  

  if (comp==1) {

    for (stanza=psp.GetStar(); stanza!=0; stanza=psp.NextStar()) {
      
      out << setw(20) << stanza->comp.nbod << endl;

      for (part=psp.GetParticle(in); part!=0; part=psp.NextParticle(in)) {
	if (part->f.get()) {
	  if (use_int) {
	    if (use_pos<static_cast<int>(part->f->iatr.size()))
	      out << setw(20) << part->f->iatr[use_pos];
	    else {
	      ostringstream msg;
	      msg << "attribute size is " << part->f->datr.size()
		  << " and you can not use pos="  << use_pos;
	      throw msg.str().c_str();
	    }
	  } else {
	    if (use_pos<static_cast<int>(part->f->datr.size()))
	      out << setw(20) << part->f->datr[use_pos];
	    else {
	      ostringstream msg;
	      msg << "attribute size is " << part->f->datr.size()
		  << " and you can not use pos="  << use_pos;
	      throw msg.str().c_str();
	    }
	  }
	} else {
	  if (use_int) {
	    if (use_pos<static_cast<int>(part->d->iatr.size()))
	      out << setw(20) << part->d->iatr[use_pos];
	    else {
	      ostringstream msg;
	      msg << "attribute size is " << part->d->datr.size()
		  << " and you can not use pos="  << use_pos;
	      throw msg.str().c_str();
	    }
	  } else {
	    if (use_pos<static_cast<int>(part->d->datr.size()))
	      out << setw(20) << part->d->datr[use_pos];
	    else {
	      ostringstream msg;
	      msg << "attribute size is " << part->d->datr.size()
		  << " and you can not use pos="  << use_pos;
	      throw msg.str().c_str();
	    }
	  }
	}
	
      }
    }
  }  

}

int
main(int argc, char **argv)
{
  char *prog = argv[0];
  double time=1e20;
  bool verbose = false;

  // Parse command line

  while (1) {

    int c = getopt(argc, argv, "t:a:c:o:ivh");

    if (c == -1) break;

    switch (c) {

    case 't':
      time = atof(optarg);
      break;

    case 'a':
      use_pos = atoi(optarg);
      break;

    case 'c':
      comp = atoi(optarg);
      break;

    case 'o':
      suffix = optarg;
      break;

    case 'i':
      use_int = true;
      break;

    case 'v':
      verbose = true;
      break;

    case '?':
    case 'h':
    default:
      Usage(prog);
    }

  }

				// Sanity check
  if (comp<0 || comp>2) {
    cerr << "Comp=" << comp << " but must be in [0,2]" << endl;
    exit(-1);
  }

  ifstream *in;

  if (optind < argc) {

    ifstream *in2 = new ifstream(argv[optind]);
    if (!*in2) {
      cerr << "Error opening file <" << argv[optind] << "> for input\n";
      exit(-1);
    }

    cerr << "Using filename: " << argv[optind] << endl;

				// Assign file stream to input stream
    in = in2;

  }

				// Read the phase space file
				// -------------------------

  PSPDump psp(in, true);

  in->close();
  delete in;
				// Reopen file
				// -----------
  in = new ifstream(argv[optind]);

  cerr << endl << "Best fit dump to <" << time << "> has time <" 
       << psp.SetTime(time) << ">" << endl;
  
				// Write a summary
				// ---------------
  if (verbose) psp.PrintSummaryCurrent(in, cerr);
  
  try {
    write_array(in, psp);
  }
  catch (const char *error) {
    cout << "*** Error: " << error << endl;
    exit(1);
  }

  cerr << "Done" << endl;

  return 0;
}
  
