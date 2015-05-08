#ifndef _AsciiHisto_H
#define _AsciiHisto_H

#include <iostream>
#include <iomanip>
#include <vector>

#include <boost/shared_ptr.hpp>

template <typename T = double>
class AsciiHisto
{
private:

  std::vector<T> bins;
  T beg, end, dx;
  bool logH;

public:

  AsciiHisto(const std::vector<T>& data, int nbins=0, bool logV=false) : 
    logH(logV)
  {
    if (nbins==0) nbins = floor(sqrt(data.size()));
    
    beg = *std::min_element(data.begin(), data.end());
    end = *std::max_element(data.begin(), data.end());
  
    if (beg<=0.0) 
      logH = false;
    else if (logH) {
      beg = log(beg);
      end = log(end);
    }
    
    dx  = (end - beg)/nbins;
    bins.resize(nbins, 0);

    for (auto v : data) {
      if (logH) v = log(v);
      bins[floor( (v - beg)/dx )]++;
    }
  }

  void operator()(std::ostream &out, unsigned hwid=60)
  {
    unsigned wid = 10;
    if (logH) wid = 16;
    unsigned maxb = *std::max_element(bins.begin(), bins.end());
    double   hscl = static_cast<double>(hwid)/maxb;

    std::string sep(" | ");
    std::string header(2*wid + 4 + 2*sep.size() + hwid + 10, '-');

    out << header << std::endl << std::left;
    
    for (size_t j=0; j<bins.size(); j++) {
      std::ostringstream sout;
      T lower = beg + dx*(j+0);
      T upper = beg + dx*(j+1);
      sout << "[" << std::setw(wid) << lower << ", "
	   << std::setw(wid) << upper << "]" << sep;
      
      unsigned cnt = floor(hscl*bins[j]);
      out << sout.str() << std::setw(hwid) << std::string(cnt, '=')
	  << sep << std::setw(10) << bins[j] << std::endl;
    }
    
    out << header << std::endl;
  }
};

typedef boost::shared_ptr< AsciiHisto<unsigned> > chistoUPtr;
typedef boost::shared_ptr< AsciiHisto<double  > > chistoDPtr;
typedef boost::shared_ptr< AsciiHisto<int     > > chistoIPtr;

#endif
