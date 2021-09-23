/*
  Info on each stanza of a phase space dump

  MDWeinberg 03/17/02
*/

#include <stdlib.h>
#include <unistd.h>

#include <iostream>
#include <fstream>
#include <iomanip>
#include <vector>
#include <string>
#include <memory>

#include <boost/random/mersenne_twister.hpp>

#include <PSP.H>

//-------------
// Help message
//-------------

void Usage(char* prog) {
  cerr << prog << ": [[-v] [-S] [-s] [-T] [-d data_dir]] filename\n";
  cerr << "        -v      verbose output\n";
  cerr << "        -S      assume split PSP files\n";
  cerr << "        -s      particle and velocity statistics\n";
  cerr << "        -T      print time info only\n";
  cerr << "        -d dir  data directory\n";
  exit(-1);
}

int
main(int argc, char *argv[])
{
  bool spl = false;
  bool stats = false;
  bool timeonly = false;
  bool verbose = false;
  std::string new_dir("./");
  int c;
  
  while (1) {
    c = getopt(argc, argv, "vsSTd:h");
    if (c == -1) break;

    switch (c) {
      
    case 'v':
      verbose = true;
      break;

    case 'S':
      spl = true;
      break;

    case 's':
      stats = true;
      break;

    case 'T':
      timeonly = true;
      break;

    case 'd':
      new_dir.erase();
      new_dir = string(optarg);
      break;

    case 'h':
    case '?':
    default:
      Usage(argv[0]); 
    
    }
  }

  if (optind >= argc) Usage(argv[0]);

  std::string file(argv[optind]);

  if (not spl and file.find("SPL")!=std::string::npos) spl = true;

  std::shared_ptr<PSP> psp;
  try {
    if (spl) psp = std::make_shared<PSPspl>(argv[optind], new_dir, verbose);
    else     psp = std::make_shared<PSPout>(argv[optind], verbose);
  }
  catch (const std::exception& e)  {
    std::cout << "pspinfo runtime error: " << e.what() << std::endl;
    exit(-1);
  }
  catch (...) {
    std::cout << "pspinfo unknown error" << std::endl;
    exit(-2);
  }

  cerr << "Filename: " << file << endl;

  psp->PrintSummary(cout, stats, timeonly);

  return 0;
}
  
