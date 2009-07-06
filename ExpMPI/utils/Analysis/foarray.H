#include <vector>
#include <string>
#include <algorithm>

using namespace std;

//
// Maintain an array of pointers to ofstreams
//
class foarray
{
private:
  vector<ofstream*> fo;
  unsigned fnum;

public:

  foarray(vector<string>& names, bool append=false)
  {
    fnum = names.size();
    fo = vector<ofstream*>(fnum);
    for (unsigned n=0; n<fnum; n++) {
      if (append)
	fo[n] = new ofstream(names[n].c_str(), ios::app);
      else
	fo[n] = new ofstream(names[n].c_str());
      if (!*fo[n]) {
	cerr << "Error opening <" << names[n] << "> for output\n";
	exit(-1);
      }
    }
  }

  ~foarray()
  {
    for (unsigned n=0; n<fnum; n++) delete fo[n];
  }

  ofstream& operator[](int n) 
  { 
    n = min<int>(max<int>(n, 0), n);
    return *(fo[n]); 
  }
};


