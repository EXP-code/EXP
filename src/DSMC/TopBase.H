#ifndef _TopBase_H
#define _TopBase_H

#include <memory>
#include <sstream>
#include <vector>
#include <map>

//! Access to the TOPbase data
class TopBase
{

public:

  typedef std::pair<unsigned short, unsigned short> iKey;

  // One line in the TOPbase data base for photoionization
  struct TBline 
  {
    unsigned short I, NZ, NE;
    int ISLP, ILV, NP;
    double Eph, wght;
    std::vector<double> E, S;
  };

  typedef std::shared_ptr<TBline> TBptr;

  //! A collection of lines for a single electronic config indexed by
  //! line level
  typedef std::map<int, TBptr>      TBcfg;

  //! A collection of electronic config indexed by iSLP value
  typedef std::map<int, TBcfg>      TBslp;

public:

  class Error : public std::exception 
  {
  public:
    std::string message;

    Error() {}
    Error(const std::string s) : message(s) {}
    virtual ~Error() throw() {}
    virtual const char *what() const throw()
    {
      std::ostringstream sout;
      sout << "TopBase: " << message;
      return sout.str().c_str();
    }
  };

  class FileOpen : public Error
  {
  public:
    FileOpen(const std::string s) {
      std::ostringstream sout;
      sout << "Failure to open file <" << s << ">";
      message = sout.str();
    }
  };

  class NoIon : public Error
  {
  public:
    NoIon(const unsigned short NZ, const unsigned short NE) {
      std::ostringstream sout;
      sout << "Ion with NZ=" << NZ << ", NE=" << NE << " not found";
      message = sout.str();
    }
  };

  class NoSLP : public Error
  {
  public:
    NoSLP(const unsigned short NZ, const unsigned short NE, 
	  const int I) 
    {
      std::ostringstream sout;
      sout << "Config=" << I << "not found in Ion(" 
	   << NZ << ", " << NE << ")";
      message = sout.str();
    }
  };

  class NoLine : public Error
  {
  public:
    NoLine(const unsigned short NZ, const unsigned short NE, 
	   const int I, const int L) 
    {
      std::ostringstream sout;
      sout << "Line level=" << L << "not found in Ion(" 
	   << NZ << ", " << NE << ") with config=" << I;
      message = sout.str();
    }
  };

private:

  void readData();
  typedef std::map<iKey, TBslp> TBmap;
  typedef TBmap::iterator       TBmapItr;
  TBmap ions;

  std::map<iKey, double>        SWlow;

public:

  //! Constructor
  TopBase() 
  { 
    readData(); 
  }

  //! Access to an ion
  const TBslp* Ion(unsigned short NZ, unsigned short NE)
  {
    TBmapItr ion = ions.find(iKey(NZ, NE));
    if (ion != ions.end()) return &ion->second;
    else                   return 0;
  }

  //! Compute Milne fb cross section
  double sigmaFB(const iKey& key, double E);

  //! Print list of lines
  void printInfo();

  //! Print line
  void printLine(unsigned short NZ, unsigned short NE, int I, int L,
		 const std::string& file);

};

#endif
