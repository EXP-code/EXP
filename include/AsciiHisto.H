#ifndef _AsciiHisto_H
#define _AsciiHisto_H

#include <iostream>
#include <iomanip>
#include <vector>

#include <memory>

//! Process a std::vector into an ASCII histogram
template <typename T = double>
class AsciiHisto
{
private:

  std::vector<unsigned> bins;
  T beg, end, dx;
  bool logH;

  void makeBins(const std::vector<T>& d, int nbins)
  {
    if (d.size()) {
      beg = *std::min_element(d.begin(), d.end());
      end = *std::max_element(d.begin(), d.end());
  
      if (beg<=0.0) 
	logH = false;
      else if (logH) {
	beg = log(beg);
	end = log(end);
      }
      
      if (nbins==0) nbins = floor(sqrt(d.size()));
      
      dx  = (end - beg)/nbins;
      if (dx==0) {
	nbins = end - beg + 1;
	dx = 1;
      }
      
      bins.resize(nbins, 0);

      for (auto v : d) {
	if (not std::isnan(v) and not std::isinf(v)) {
	  if (logH) {
	    if (v<=0.0) continue;
	    v = log(v);
	  }
	  int indx = floor( (v - beg)/dx );
	  //! Sanity check
	  if (indx>=0 and indx < nbins) bins[floor( (v - beg)/dx )]++;
	}
      }
    }
  }

public:

  //! Constructor with bin number and trim
  AsciiHisto(const std::vector<T>& data, int nbins, double trim,
	     bool logV=false) : logH(logV)
  {
    if (trim>0.0) {

      std::vector<T> d(data);
      std::sort(d.begin(), d.end());

      size_t lower = floor(trim*d.size());
      size_t upper = floor((1.0-trim)*d.size());

      d.erase(d.begin()+upper, d.end());
      d.erase(d.begin(), d.begin()+lower);

      makeBins(d, nbins);
    } else {
      makeBins(data, nbins);
    }
  }


  //! Simple constructor
  AsciiHisto(const std::vector<T>& d, int nbins=0, bool logV=false) : logH(logV)
  {
    makeBins(d, nbins);
  }

  //! Constructor for pre-existing bins
  AsciiHisto(const std::vector<unsigned>& inbins,
	     const T& inBeg, const T& inEnd, bool logV=false) : logH(logV)
  {
    bins = inbins;
    beg  = inBeg;
    end  = inEnd;
    dx   = (end - beg)/bins.size();
  }

  //! Print the distribution
  void operator()(std::ostream &out, unsigned hwid=60)
  {
    if (bins.size() == 0) return;

    unsigned wid  = logH ? 16 : 12;
    unsigned maxb = *std::max_element(bins.begin(), bins.end());
    double   hscl = static_cast<double>(hwid)/maxb;

    std::string sep(" | ");
    std::string header(2*wid + 4 + 2*sep.size() + hwid + 10, '-');

    out << "[" 
	<< std::setw(wid) << "bin low" << ", "
	<< std::setw(wid) << "bin high"
	<< "]" 
	<< std::setw(2*sep.size() + hwid) << "" << std::left
	<< std::setw(10) << "Counts" << std::endl;

    out << header << std::endl;
    
    for (size_t j=0; j<bins.size(); j++) {
      std::ostringstream sout;
      T lower = beg + dx*(j+0);
      T upper = beg + dx*(j+1);

      if (logH) {		// Print out decimal bin edges for log
	lower = exp(lower);	// scaled bins
	upper = exp(upper);
      }

      sout << "[" 
	   << std::setw(wid) << lower << ", "
	   << std::setw(wid) << upper << "]" << sep;
      
      unsigned cnt = floor(hscl*bins[j]);
      out << sout.str() << std::setw(hwid) << std::string(cnt, '=')
	  << sep << std::setw(10) << bins[j] << std::endl;
    }
    
    out << header << std::endl;
  }
};

typedef std::shared_ptr< AsciiHisto<unsigned> > ahistoUPtr;
typedef std::shared_ptr< AsciiHisto<double  > > ahistoDPtr;
typedef std::shared_ptr< AsciiHisto<int     > > ahistoIPtr;

#endif
