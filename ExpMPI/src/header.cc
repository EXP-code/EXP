#include <iostream>
#include <iomanip>
#include <header.H>

int ComponentHeader::defaultInfoSize = 1024;

bool ComponentHeader::write(ostream *out)
{
  out->write(&nbod,  sizeof(int));
  out->write(&niatr, sizeof(int));
  out->write(&ndatr, sizeof(int));
  out->write(&ninfochar, sizeof(int));
  out->write(info, ninfochar*sizeof(char));

  if (*out)
    return true;
  else
    return false;
}

bool ComponentHeader::read(istream *in)
{
  int ninfo;

  in->read(&nbod,  sizeof(int));		if (!*in) return false;
  in->read(&niatr, sizeof(int));		if (!*in) return false;
  in->read(&ndatr, sizeof(int));		if (!*in) return false;
  in->read(&ninfo, sizeof(int));		if (!*in) return false;
  cout << "Ninfo=" << ninfo << "  Ninfochar=" << ninfochar << endl;
  if (ninfo != ninfochar) {
    delete [] info;
    ninfochar = ninfo;
    info = new char [ninfochar];
  }
  in->read(info, ninfochar*sizeof(char));	if (!*in) return false;

  return true;
}

