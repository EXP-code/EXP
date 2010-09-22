#include "uint128.H"

#include <utility>
#include <memory>
#include <cmath>

std::string 
uint128::toString (unsigned int radix, const char* prefix) const throw () 
{
  if (radix < 2 || radix > 37) return "(invalid radix)";
  
  uint128 r;
  uint128 ii = *this;
  std::string sz;
  
  if (!*this) {
    sz = sz + '0';
  } else {
    while (!!ii) {
      ii = ii.div (radix, r);
      sz = 
	static_cast<char>(r.toUint () + ((r.toUint () > 9) ? 'a' - 10 : '0'))
	+ sz;
    };
  }
  
  sz = prefix + sz;

  return sz;
};


uint128::uint128 (const char * sz) throw ()
  : lo (0u), hi (0u) {
  
  if (!sz) return;
  if (!sz [0]) return;
  
  unsigned int radix = 10;
  unsigned int i = 0;
  bool minus = false;
  
  if (sz [i] == '-') {
    ++i;
    minus = true;
  };
  
  if (sz [i] == '0') {
    radix = 8;
    ++i;
    if (sz [i] == 'x') {
      radix = 16;
      ++i;
    };
  };
  
  for (; i < strlen (sz); ++i) {
    unsigned int n = 0;
    if (sz [i] >= '0' && sz [i] <= '9')
      n = sz [i] - '0';
    else if (sz [i] >= 'a' && sz [i] <= 'a' + (int) radix - 10)
      n = sz [i] - 'a' + 10;
    else if (sz [i] >= 'A' && sz [i] <= 'A' + (int) radix - 10)
      n = sz [i] - 'A' + 10;
    else
      break;
    
    (*this) *= radix;
    (*this) += n;
  };
  
  if (minus)
    *this = 0u - *this;
  
  return;
};

uint128::uint128 (const float a) throw ()
  : lo ((uint64_t) fmodf (a, 18446744073709551616.0f)),
    hi ((uint64_t) (a / 18446744073709551616.0f)) {};

uint128::uint128 (const double & a) throw ()
  : lo ((uint64_t) fmod (a, 18446744073709551616.0)),
    hi ((uint64_t) (a / 18446744073709551616.0)) {};

uint128::uint128 (const long double & a) throw ()
  : lo ((uint64_t) fmodl (a, 18446744073709551616.0l)),
    hi ((uint64_t) (a / 18446744073709551616.0l)) {};

float uint128::toFloat () const throw () {
  return (float) this->hi * 18446744073709551616.0f
    + (float) this->lo;
};

double uint128::toDouble () const throw () {
  return (double) this->hi * 18446744073709551616.0
    + (double) this->lo;
};

long double uint128::toLongDouble () const throw () {
  return (long double) this->hi * 18446744073709551616.0l
    + (long double) this->lo;
};

uint128 uint128::operator - () const throw () { 
  return uint128(0u) - *this;
}
  

uint128 uint128::operator ~ () const throw () {
  return uint128 (~this->lo, ~this->hi);
};

uint128 & uint128::operator ++ () {
  ++this->lo;
  if (!this->lo)
    ++this->hi;
  
  return *this;
};

uint128 & uint128::operator -- () {
  if (!this->lo)
    --this->hi;
  --this->lo;
  
  return *this;
};

uint128 uint128::operator ++ (int) {
  uint128 b = *this;
  ++ *this;
  
  return b;
};

uint128 uint128::operator -- (int) {
  uint128 b = *this;
  -- *this;
  
  return b;
};

uint128 & uint128::operator += (const uint128 & b) throw () {
  uint64_t old_lo = this->lo;
  
  this->lo += b.lo;
  this->hi += b.hi + (this->lo < old_lo);
  
  return *this;
};

uint128 & uint128::operator *= (const uint128 & b) throw () {
  if (!b)
    return *this = 0u;
  if (b == 1u)
    return *this;
  
  uint128 a = *this;
  uint128 t = b;
  
  this->lo = 0ull;
  this->hi = 0ull;
  
  for (unsigned int i = 0; i < 128; ++i) {
    if (t.lo & 1)
      *this += a << i;
    
    t >>= 1;
  };
  
  return *this;
};


uint128 & uint128::operator = (const uint128 & b) throw () 
{
  this->lo = b.lo;
  this->hi = b.hi;
  
  return *this;
};


uint128 uint128::div (const uint128 & ds, uint128 & remainder) const throw () {
  if (!ds)
    return 1u / (unsigned int) ds.lo;
  
  uint128 dd = *this;
  
  // only remainder
  if (ds > dd) {
    remainder = *this;
    return 0ul;
  };
  
  uint128 r = 0ul;
  uint128 q = 0ul;
  //    while (dd >= ds) { dd -= ds; q += 1; }; // extreme slow version
  
  unsigned int b = 127;
  while (r < ds) {
    r <<= 1;
    if (dd.bit (b--))
      r.lo |= 1;
  };
  ++b;
  
  while (true)
    if (r < ds) {
      if (!(b--)) break;
      
      r <<= 1;
      if (dd.bit (b))
	r.lo |= 1;
      
    } else {
      r -= ds;
      q.bit (b, true);
    };
  
  remainder = r;
  return q;
};

bool uint128::bit (unsigned int n) const throw () {
  n &= 0x7F;
  
  if (n < 64)
    return this->lo & (1ull << n);
  else
    return this->hi & (1ull << (n - 64));
};

void uint128::bit (unsigned int n, bool val) throw () {
  n &= 0x7F;
  
  if (val) {
    if (n < 64) this->lo |= (1ull << n);
    else this->hi |= (1ull << (n - 64));
  } else {
    if (n < 64) this->lo &= ~(1ull << n);
    else this->hi &= ~(1ull << (n - 64));
  };
};


uint128 & uint128::operator >>= (unsigned int n) throw () {
  n &= 0x7F;
  
  if (n > 63) {
    n -= 64;
    this->lo = this->hi;
    this->hi = 0ul;
  };
  
  if (n) {
    // shift low qword
    this->lo >>= n;
    
    // get lower N bits of high qword
    uint64_t mask = 0ul;
    for (unsigned int i = 0; i < n; ++i) mask |= (1 << i);
    
    // and add them to low qword
    this->lo |= (this->hi & mask) << (64 - n);
    
    // and finally shift also high qword
    this->hi >>= n;
  };
  
  return *this;
};

uint128 & uint128::operator <<= (unsigned int n) throw () {
  n &= 0x7F;
  
  if (n > 63) {
    n -= 64;
    this->hi = this->lo;
    this->lo = 0ul;
  };
  
  if (n) {
    // shift high qword
    this->hi <<= n;
    
    // get higher N bits of low qword
    uint64_t mask = 0ul;
    for (unsigned int i = 0; i < n; ++i) mask |= (1 << (63 - i));
    
    // and add them to high qword
    this->hi |= (this->lo & mask) >> (64 - n);
    
    // and finally shift also low qword
    this->lo <<= n;
  };
  
  return *this;
};

bool uint128::operator ! () const throw () {
  return !(this->hi || this->lo);
};

uint128 & uint128::operator -= (const uint128 & b) throw () {
  if (this->lo < b.lo) this->hi--;
  this->hi -= b.hi;
  this->lo -= b.lo;
  
  return *this;
};

uint128 & uint128::operator |= (const uint128 & b) throw () {
  this->hi |= b.hi;
  this->lo |= b.lo;
  
  return *this;
};

uint128 & uint128::operator &= (const uint128 & b) throw () {
  this->hi &= b.hi;
  this->lo &= b.lo;
  
  return *this;
};

uint128 & uint128::operator ^= (const uint128 & b) throw () {
  this->hi ^= b.hi;
  this->lo ^= b.lo;
  
  return *this;
};

bool operator <  (const uint128 & a, const uint128 & b) throw () {
  return (a.hi == b.hi) ? (a.lo < b.lo) : (a.hi < b.hi);
};

bool operator == (const uint128 & a, const uint128 & b) throw () {
  return (a.hi == b.hi) && (a.lo == b.lo);
};
bool operator && (const uint128 & a, const uint128 & b) throw () {
  return (a.hi || a.lo) && (b.hi || b.lo);
};
bool operator || (const uint128 & a, const uint128 & b) throw () {
  return (a.hi || a.lo) || (b.hi || b.lo);
};


std::ostream &operator<<(std::ostream& out, const uint128& u)
{
  out << u.toString(16, "0x");
}
