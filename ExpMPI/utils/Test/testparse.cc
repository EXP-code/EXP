#include <iostream>
#include <iomanip>
#include <fstream>
#include <string>

#include <unistd.h>
#include <getopt.h>

#include <ParamParse.H>

int myid=0;
char threading_on = 0;
pthread_mutex_t mem_lock;


//===========================================================================

void usage(char *prog)
{
  cout << "Usage:\n\n"
       << prog << " [options]\n\n"
       << setw(15) << "Option" << setw(10) << "Argument" << setw(10) << " " 
       << setiosflags(ios::left)
       << setw(40) << "Description" << endl << endl
       << resetiosflags(ios::left)
       << setw(15) << "-f or --file" << setw(10) << "string" << setw(10) << " " 
       << setiosflags(ios::left)
       << setw(40) << "File name containin the paramters to parse" << endl
       << resetiosflags(ios::left)
       << "" << endl;

  exit(0);
}

int 
main(int argc, char** argv)
{
  string parmfile("in.file");

  int c;
  while (1) {
    int this_option_optind = optind ? optind : 1;
    int option_index = 0;
    static struct option long_options[] = {
      {"file", 1, 0, 0},
      {0, 0, 0, 0}
    };

    c = getopt_long (argc, argv, "f:",
		     long_options, &option_index);

    if (c == -1) break;

    switch (c) {
    case 0:
      {
	string optname(long_options[option_index].name);

	if (!optname.compare("file")) {
	  parmfile = string(optarg);
	} else {
	  cout << "Option " << long_options[option_index].name;
	  if (optarg) cout << " with arg " << optarg;
	  cout << " is not defined " << endl;
	  exit(0);
	}
      }
      break;

    case 'f':
      parmfile = string(optarg);
      break;

    case 'h':
    default:
      usage(argv[0]);
    }

  }

  ifstream in(parmfile.c_str());
  if (!in) {
    char mbuf[512];
    cerr << argv[0] << ": can not open parameter file <" << parmfile << ">\n";
    cerr << argv[0] << ": pwd is <" << getcwd(mbuf, 512) << ">\n";
    exit(-1);
  }


  ParamParse *parse = new ParamParse(&in, ":");

  cout << "Parameter database:\n"
       << "-------------------\n\n";
  parse->print_database(cout);
  

  delete parse;

  return 0;
}

