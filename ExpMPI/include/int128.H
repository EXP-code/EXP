#ifndef INT128_H
#define INT128_H

#include <ostream>
#include <string>
#include <cstdlib>
#include <cstring>
#include <stdint.h>

/**
 Class defining an unsigned 128-bit integer
 */
class int128
{
  friend std::ostream &operator<<( std::ostream&, const int128& );

private:
  // The signed 128bit integer is represented in one unsigned and one
  // signed 64 bit ints
  uint64_t    lo;
  int64_t     hi;
        
protected:

  friend bool operator <  (const int128 &, const int128 &) throw ();
  friend bool operator == (const int128 &, const int128 &) throw ();
  friend bool operator || (const int128 &, const int128 &) throw ();
  friend bool operator && (const int128 &, const int128 &) throw ();
  
public:
  // Constructors
  inline int128 () throw () {};
  inline int128 (const int128 & a) throw () : lo (a.lo), hi (a.hi) {};
  
  inline int128 (const unsigned int & a) throw () : lo (a), hi (0ll) {};
  inline int128 (const signed int & a) throw () : lo (a), hi (0ll) {
    if (a < 0) this->hi = -1ll;
  };
  
  inline int128 (const uint64_t & a) throw () : lo (a), hi (0ll) {};
  inline int128 (const int64_t  & a) throw () : lo (a), hi (0ll) {
    if (a < 0) this->hi = -1ll;
  };
  
  int128 (const float a) throw ();
  int128 (const double & a) throw ();
  int128 (const long double & a) throw ();
  int128 (const char * sz) throw ();
  
        
private:
  // Used internally for contructing return instances
  int128 (const uint64_t & a, const int64_t & b) throw ()
    : lo (a), hi (b) {};
        
public:
  // Operators
  bool operator ! () const throw ();
  
  int128 operator - () const throw ();
  int128 operator ~ () const throw ();
  
  int128 & operator ++ ();
  int128 & operator -- ();
  int128 operator ++ (int);
  int128 operator -- (int);
  
  int128 & operator += (const int128 & b) throw ();
  int128 & operator *= (const int128 & b) throw ();
  
  int128 & operator >>= (unsigned int n) throw ();
  int128 & operator <<= (unsigned int n) throw ();
  
  int128 & operator |= (const int128 & b) throw ();
  int128 & operator &= (const int128 & b) throw ();
  int128 & operator ^= (const int128 & b) throw ();
  
  int128 & operator =  (const int128 & b) throw ();

  // Inline operators
  inline const int128 & operator + () const throw () { return *this; };
  
  inline int128 & operator -= (const int128 & b) throw () {
    return *this += (-b);
  };
  inline int128 & operator /= (const int128 & b) throw () {
    int128 dummy;
    *this = this->div (b, dummy);
    return *this;
  };
  inline int128 & operator %= (const int128 & b) throw () {
    this->div (b, *this);
    return *this;
  };
        
  // Output and type conversion methods
  int toInt () const throw () {  return (int) this->lo; };
  int64_t toInt64 () const throw () {  return (int64_t) this->lo; };
  
  std::string toString (unsigned int radix, const char* prefix = "") const throw ();

  std::string toHex() throw () {
    return toString(16, "0x");
  }

  std::string toOct() throw () {
    return toString(8, "0");
  }

  std::string toBin() throw () {
    return toString(2);
  }

  std::string toDec() throw () {
    return toString(10);
  }

  float toFloat () const throw ();
  double toDouble () const throw ();
  long double toLongDouble () const throw ();
  
  // Arithmetic methods
  int128  div (const int128 &, int128 &) const throw ();
  
  // Bit operations
  bool    bit (unsigned int n) const throw ();
  void    bit (unsigned int n, bool val) throw ();
}
#ifdef __GNUC__
  __attribute__ ((__aligned__ (16), __packed__))
#endif
  ;


// GLOBAL OPERATORS

bool operator <  (const int128 & a, const int128 & b) throw ();
bool operator == (const int128 & a, const int128 & b) throw ();
bool operator || (const int128 & a, const int128 & b) throw ();
bool operator && (const int128 & a, const int128 & b) throw ();

inline bool operator <= (const int128 & a, const int128 & b) throw () {
    return !(b < a); };
inline bool operator >= (const int128 & a, const int128 & b) throw () {
    return !(a < b); };

// GLOBAL OPERATOR INLINES

inline int128 operator + (const int128 & a, const int128 & b) throw () {
  return int128 (a) += b; };
inline int128 operator - (const int128 & a, const int128 & b) throw () {
  return int128 (a) -= b; };
inline int128 operator * (const int128 & a, const int128 & b) throw () {
  return int128 (a) *= b; };
inline int128 operator / (const int128 & a, const int128 & b) throw () {
  return int128 (a) /= b; };
inline int128 operator % (const int128 & a, const int128 & b) throw () {
  return int128 (a) %= b; };

inline int128 operator >> (const int128 & a, unsigned int n) throw () {
  return int128 (a) >>= n; };
inline int128 operator << (const int128 & a, unsigned int n) throw () {
  return int128 (a) <<= n; };

inline int128 operator & (const int128 & a, const int128 & b) throw () {
  return int128 (a) &= b; };
inline int128 operator | (const int128 & a, const int128 & b) throw () {
  return int128 (a) |= b; };
inline int128 operator ^ (const int128 & a, const int128 & b) throw () {
  return int128 (a) ^= b; };

inline bool operator >  (const int128 & a, const int128 & b) throw () {
  return   b < a; };
inline bool operator != (const int128 & a, const int128 & b) throw () {
  return !(a == b); };


// I'm not sure I need this . . .

typedef int128 int128_t;

#endif
