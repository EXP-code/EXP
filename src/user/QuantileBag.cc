#include <exception>
#include <algorithm>
#include <iostream>
#include <iomanip>
#include <sstream>
#include <string>
#include <cmath>

#include "QuantileBag.H"

using namespace NTC;

#include "QuantileBag.H"

// Count instances for debugging
unsigned int QuantileBag::instance = 0;

// Main constructor
//
QuantileBag::QuantileBag()
{
  // No data to start
  M = 0;

  // Initialize quantiles
  for (auto v : qs) quant[v] = Quantile(v);
  hist = Quantile(Nequal);

  // For debugging
  instance++;
}

// Copy constructor
//
QuantileBag::QuantileBag(const QuantileBag& p)
{
  quant = p.quant;
  hist  = p.hist;
  M     = p.M;

  // For debugging
  instance++;
}

// Copy operator
//
QuantileBag &QuantileBag::operator=(const QuantileBag& p)
{
  quant = p.quant;
  hist  = p.hist;
  M     = p.M;

  // For debugging
  instance++;

  return *this;
}

// Add a value to all quantiles
//
void QuantileBag::add(double x)
{
  for (auto q : quant) q.second.add(x);
  hist.add(x);
  M++;
}

// Retrieve the value with desired quantile
//
double QuantileBag::operator()(double p)
{
  // Find the closest quantile value to p
  //
  std::map<double, Quantile>::iterator it = quant.find(p);
  if (it != quant.end()) return it->second();
  return hist(p);
}

// Node sends its internal data to root
//
void QuantileBag::send()
{
  for (auto q : quant) q.second.send();
  hist.send();
}
    
// Root intializes itself from node's data
//
void QuantileBag::recv(int id)
{
  for (auto q : quant) q.second.recv(id);
  hist.recv(id);
}

// Only used to decrement live count
//
QuantileBag::~QuantileBag()
{
  instance--;			// Count instances for debugging only
}

