#ifndef UINT128_H
#define UINT128_H

#include <ostream>
#include <string>
#include <cstdlib>
#include <cstring>
#include <stdint.h>

/**
 Class defining an unsigned 128-bit integer
 */
class uint128 
{
  friend std::ostream &operator<<( std::ostream&, const uint128& );

private:
  // The unsigned 128bit integer is represented in two unsigned 64 bit ints
  uint64_t    lo;
  uint64_t    hi;
        
protected:

  friend bool operator <  (const uint128 &, const uint128 &) throw ();
  friend bool operator == (const uint128 &, const uint128 &) throw ();
  friend bool operator || (const uint128 &, const uint128 &) throw ();
  friend bool operator && (const uint128 &, const uint128 &) throw ();
  
public:
  // Constructors
  inline uint128 () throw () {};
  inline uint128 (const uint128  & a) throw () : lo (a.lo), hi (a.hi) {};
  inline uint128 (const uint32_t & a) throw () : lo (a), hi (0ull) {};
  inline uint128 (const uint64_t & a) throw () : lo (a), hi (0ull) {};
  
  uint128 (const float a) throw ();
  uint128 (const double & a) throw ();
  uint128 (const long double & a) throw ();
  uint128 (const char * sz) throw ();
  
  void decompose(unsigned long & a, unsigned long & b)
  {
    a = lo;
    b = hi;
  }
  
  void recompose(const unsigned long a, const unsigned long b)
  {
    lo = a;
    hi = b;
  }
  
private:

  // For constructing an instance from communicated pieces
  uint128 (const uint64_t & a, const uint64_t & b) throw ()
    : lo (a), hi (b) {};


public:
  // Operators
  bool operator ! () const throw ();
  
  uint128 operator - () const throw ();
  uint128 operator ~ () const throw ();
  
  uint128 & operator ++ ();
  uint128 & operator -- ();
  uint128 operator ++ (int);
  uint128 operator -- (int);
  
  uint128 & operator += (const uint128 & b) throw ();
  uint128 & operator -= (const uint128 & b) throw ();
  uint128 & operator *= (const uint128 & b) throw ();
  
  uint128 & operator >>= (unsigned int n) throw ();
  uint128 & operator <<= (unsigned int n) throw ();
  
  uint128 & operator |= (const uint128 & b) throw ();
  uint128 & operator &= (const uint128 & b) throw ();
  uint128 & operator ^= (const uint128 & b) throw ();

  uint128 & operator =  (const uint128 & b) throw ();

  // Inline operators
  inline const uint128 & operator + () const throw () { return *this; };

  inline uint128 & operator /= (const uint128 & b) throw () {
    uint128 dummy;
    *this = this->div (b, dummy);
    return *this;
  };
  inline uint128 & operator %= (const uint128 & b) throw () {
    this->div (b, *this);
    return *this;
  };
        
  // Output and type conversion methods

  unsigned int toUint () const throw () {
    return (unsigned int) this->lo; };

  uint64_t toUint64 () const throw () { 
    return (uint64_t) this->lo; };

  std::string toString (unsigned int radix, const char* prefix = "") const throw ();

  std::string toHex() const throw () {
    return toString(16, "0x");
  }

  std::string toOct() const throw () {
    return toString(8, "0");
  }

  std::string toBin() const throw () {
    return toString(2);
  }

  std::string toDec() const throw () {
    return toString(10);
  }

  float toFloat () const throw ();
  double toDouble () const throw ();
  long double toLongDouble () const throw ();
                
  // Arithmetic methods
  uint128  div (const uint128 &, uint128 &) const throw ();
  
  // Bit operations
  bool    bit (unsigned int n) const throw ();
  void    bit (unsigned int n, bool val) throw ();

}
#ifdef __GNUC__
  __attribute__ ((__aligned__ (16), __packed__))
#endif
  ;


// Global operators

bool operator <  (const uint128 & a, const uint128 & b) throw ();
bool operator == (const uint128 & a, const uint128 & b) throw ();
bool operator || (const uint128 & a, const uint128 & b) throw ();
bool operator && (const uint128 & a, const uint128 & b) throw ();

inline bool operator <= (const uint128 & a, const uint128 & b) throw () {
    return !(b < a); };
inline bool operator >= (const uint128 & a, const uint128 & b) throw () {
    return !(a < b); };

// Global inlined operators

inline uint128 operator + (const uint128 & a, const uint128 & b) throw () {
  return uint128 (a) += b; };
inline uint128 operator - (const uint128 & a, const uint128 & b) throw () {
  return uint128 (a) -= b; };
inline uint128 operator * (const uint128 & a, const uint128 & b) throw () {
  return uint128 (a) *= b; };
inline uint128 operator / (const uint128 & a, const uint128 & b) throw () {
  return uint128 (a) /= b; };
inline uint128 operator % (const uint128 & a, const uint128 & b) throw () {
  return uint128 (a) %= b; };

inline uint128 operator >> (const uint128 & a, unsigned int n) throw () {
  return uint128 (a) >>= n; };
inline uint128 operator << (const uint128 & a, unsigned int n) throw () {
  return uint128 (a) <<= n; };

inline uint128 operator & (const uint128 & a, const uint128 & b) throw () {
  return uint128 (a) &= b; };
inline uint128 operator | (const uint128 & a, const uint128 & b) throw () {
  return uint128 (a) |= b; };
inline uint128 operator ^ (const uint128 & a, const uint128 & b) throw () {
  return uint128 (a) ^= b; };

inline bool operator >  (const uint128 & a, const uint128 & b) throw () {
  return   b < a; };
inline bool operator != (const uint128 & a, const uint128 & b) throw () {
  return !(a == b); };


// I'm not sure I need this . . .

typedef uint128 uint128_t;

#endif
