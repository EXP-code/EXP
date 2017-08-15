#include <exception>
#include <algorithm>
#include <iostream>
#include <iomanip>
#include <sstream>
#include <string>
#include <cmath>

#include "Quantile.H"

using namespace NTC;

// Count instances for debugging
unsigned QuantileBag::instance = 0;

#include "QuantileBag.H"

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

QuantileBag(const QuantileBag& p)
{
  quant = p.quant;
  hist  = p.hist;
  M     = p.M;

  // For debugging
  instance++;
}

QuantileBag & QuantileBag::operator(const QuantileBag& p)
{
  quant = p.quant;
  hist  = p.hist;
  M     = p.M;

  // For debugging
  instance++;
}

void QuantileBag::add(double x)
{
  for (auto q : quant) q.second.add(x);
  hist.add(x);
  N++;
}

double QuantileBag::operator()(double p)
{
  // Find the closest quantile value to p
  //
  std::map<double, Quantile>::iterator it = quant.find(p);
  if (it != quant.end()) return (*it)();
  return hist(p);
}

// Node sends its internal data to root
void QuantileBag::send()
{
  for (auto q : quant) q.send();
  hist.send();
}
    
// Root intializes itself from node's data
void QuantileBag::recv(int id)
{
  for (auto q : quant) q.recv(id);
  hist.recv(id);
}

QuantileBag::~QuantileBag()
{
  instance--;			// Count instances for debugging only
}

