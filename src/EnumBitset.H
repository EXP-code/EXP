#ifndef _EnumBitset_H
#define _EnumBitset_H

#include <bitset>
#include <type_traits>

//! EnumTraits must define a max data member to set the bitset size
//! and wraps an enum class
template<typename T>
struct EnumTraits;

//! Wrap the bitset according to the provided EnumTraits class.  The
//! class allows fields in a bit set to be accessed by members of the
//! enum.
template<typename T>
class EnumBitset
{
private:

  //! The underlying bitset
  std::bitset<static_cast<typename std::underlying_type<T>::type>(EnumTraits<T>::max)> s;

  typename std::underlying_type<T>::type get_value(T v) const
  {
    return static_cast<typename std::underlying_type<T>::type>(v);
  }

public:

  //! Constructor
  EnumBitset() : s() {}

  //! Test the flag
  bool test(T pos) const
  {
    return s.test(get_value(pos));
  }

  //! operator version of the test
  bool operator[](T pos) const
  {
    return s.test(get_value(pos));
  }

  //! operator version of the test
  bool operator[](size_t i) const
  {
    if (i>=s.size()) return false;
    return s.test(i);
  }

  //! Set the flag
  EnumBitset& set(T pos)
  {
    s.set(get_value(pos));
    return *this;
  }

  //! Reset the flag
  EnumBitset& reset(T pos)
  {
    s.reset(get_value(pos));
    return *this;
  }

  //! Flip the flag
  EnumBitset& flip(T pos)
  {
    s.flip(get_value(pos));
    return *this;
  }

  //! Reset all flags
  void reset()
  {
    s.reset();
  }

  //! Return size of bitset
  size_t size()
  {
    return s.size();
  }


};

#endif
