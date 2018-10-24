// -*- C++ -*-

#ifndef Cspline_H
#define Cspline_H

#include <iostream>
#include <vector>

/** A cubic spline evaluator templated on data types X, Y. 
    Types X and Y must have operator +, -, *,
    Y must have defined a constructor that takes a scalar.
*/
template <typename X, typename Y>
class Cspline 
{
public:

  //! Null constructor
  Cspline() {}
  
  //! Constructor for spline with x and y vectors
  Cspline(const std::vector<X>& x, const std::vector<Y>& y) 
  {
    if (x.size() != y.size()) {
      std::cerr << "X[" << x.size() << "] and Y[" << y.size()
		<< "] must be the same size " << std::endl;
      return;
    }
        
    if (x.size() < 3) {
      std::cerr << "Must have at least three points for interpolation, " 
		<< "current size is " << x.size() << std::endl;
      return;
    }
        
    typedef typename std::vector<X>::difference_type size_type;
        
    size_type n = y.size() - 1;
        
    std::vector<Y> b(n), d(n), a(n), c(n+1), l(n+1), u(n+1), z(n+1);
    std::vector<X> h(n+1);
    
    l[0] = Y(1);
    u[0] = Y(0);
    z[0] = Y(0);
    h[0] = x[1] - x[0];
            
    for (size_type i = 1; i < n; i++) {
      h[i] = x[i+1] - x[i];
      l[i] = Y(2 * (x[i+1] - x[i-1])) - Y(h[i-1]) * u[i-1];
      u[i] = Y(h[i]) / l[i];
      a[i] = (Y(3) / Y(h[i])) * (y[i+1] - y[i]) - (Y(3) / Y(h[i-1])) * (y[i] - y[i-1]);
      z[i] = (a[i] - Y(h[i-1]) * z[i-1]) / l[i];
    }
    
    l[n] = Y(1);
    z[n] = c[n] = Y(0);
        
    for (size_type j = n-1; j >= 0; j--) {
      c[j] = z[j] - u[j] * c[j+1];
      b[j] = (y[j+1] - y[j]) / Y(h[j]) - (Y(h[j]) * (c[j+1] + Y(2) * c[j])) / Y(3);
      d[j] = (c[j+1] - c[j]) / Y(3 * h[j]);
    }
        
    for (size_type i = 0; i < n; i++) {
      elems.push_back(Element(x[i], y[i], b[i], c[i], d[i]));
    }        
  }

  //! Destructor
  virtual ~Cspline() {}
    
  //! Evaluation opterator
  Y operator()(const X& x) const 
  {
    return interpolate(x);
  }
    
  //! Evaluation function
  Y interpolate(const X& x) const 
  {
    if (elems.size() == 0) return Y();
    
    typename std::vector<element_type>::const_iterator it;
    it = std::lower_bound(elems.begin(), elems.end(), element_type(x));
    if (it != elems.begin()) {
      it--;
    }   
    
    return it->eval(x);
  }
    
  //! Evaluation operator for vector input, assuming xx is sorted ascending
  std::vector<Y> operator()(const std::vector<X>& xx) const 
  {
    return interpolate(xx);
  }
  
  //! Evaluation function for vector input, assuming xx is sorted ascending
  std::vector<Y> interpolate(const std::vector<X>& xx) const 
  {
    if (elems.size() == 0) return std::vector<Y>(xx.size());
        
    typename std::vector<X>::const_iterator it;
    typename std::vector<element_type>::const_iterator it2;
    it2 = elems.begin();
    std::vector<Y> ys;
    for (it = xx.begin(); it != xx.end(); it++) {
      it2 = std::lower_bound(it2, elems.end(), element_type(*it));
      if (it2 != elems.begin()) {
	it2--;
      }
      
      ys.push_back(it2->eval(*it));
    }
    
    return ys;
  }

protected:
  
  //! A cubic spline fit for an interval
  class Element 
  {
  public:
    //! The ordinate
    X x;

    //! Cubic coefficients
    Y a, b, c, d;

    //! Define the ordinate only
    Element(X _x) 
      : x(_x), a(0), b(0), c(0), d(0) {}
    
    //! The ordinate and the cubic coefficients
    Element(X _x, Y _a, Y _b, Y _c, Y _d)
      : x(_x), a(_a), b(_b), c(_c), d(_d) {}
    
    //! Evaluate the fit
    Y eval(const X& xx) const {
      X xix(xx - x);
      return a + b * xix + c * (xix * xix) + d * (xix * xix * xix);
    }
    
    //@{
    //! Sorting operators
    bool operator<(const Element& e) const {
      return x < e.x;
    }
    bool operator<(const X& xx) const {
      return x < xx;
    }
    //@}
    
  };
  
  typedef Element element_type;
  std::vector<element_type> elems;
};


#endif
