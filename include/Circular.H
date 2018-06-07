#ifndef _CIRCULAR_H
#define _CIRCULAR_H

#include <vector>

template <typename T, typename Container = std::vector<T>,
	  typename Iterator = typename Container::iterator>
class RingIterator :
  public std::iterator <std::bidirectional_iterator_tag, T, ptrdiff_t>
{
  Container& data;
  
  Iterator   cursor;
  Iterator   begin;
  Iterator   end;

public:
  
  //! Create from v and position at the beginning
  RingIterator (Container& v) :
    data(v), cursor(v.begin()), begin(v.begin()), end(v.end()) {}
  
  //! Create from v and begin at iterator i
  RingIterator (Container& v, const Iterator& i) :
    data(v), cursor(i), begin(v.begin()), end(v.end()) {}

  //! Create from v using positions i through j
  RingIterator (Container& v, const Iterator& i, const Iterator& j) :
    data(v), cursor(i), begin(i), end(j) {}

  //! Create from v and begin at position i
  RingIterator (Container& v, size_t i) :
    data(v), cursor(v.begin() + i % v.size()), begin(v.begin()), end(v.end()) {}

  //! Test for equality
  bool operator == (const RingIterator& x) const 
  { 
    return cursor == x.cursor; 
  }
  
  //! Test for inequality
  bool operator != (const RingIterator& x) const 
  {
    return ! (*this == x); 
  }
  
  //! Return current data element
  T operator*() const 
  {
    return *cursor; 
  }
  
  //! Return current iterator
  Iterator operator->() const 
  {
    return cursor; 
  }
  
  //! Pre-increment operator
  RingIterator& operator++() 
  {
    ++cursor;
    if (cursor == end) cursor = begin;
    return *this;
  }
  
  //! Post-increment operator
  RingIterator operator++(int) 
  {
    RingIterator ring = *this;
    ++*this;
    return ring;
  }

  //! Pre-decrement operator
  RingIterator& operator--() 
  {
    if (cursor == begin) cursor = end;
    --cursor;
    return *this;
  }

  //! Post-decrement operator
  RingIterator operator--(int) 
  {
    RingIterator ring = *this;
    --*this; 
    return ring;
  }

  //! Insert a new element at the current position
  RingIterator insert (const T& x)
  {
    return RingIterator (data, data.insert (cursor, x));
  }

  //! Erase the current element
  RingIterator erase() 
  {
    return RingIterator (data, data.erase (cursor));
  }
};

template <typename T, typename Iterator>
class CyclicIterator :
  public std::iterator <std::bidirectional_iterator_tag, T, ptrdiff_t>
{
  Iterator   cursor;
  Iterator   begin;
  Iterator   end;

public:

  //! Create using positions i through j
  CyclicIterator (const Iterator& i, const Iterator& j) :
    cursor(i), begin(i), end(j) {}
  
  //! Test for equality
  bool operator == (const CyclicIterator& x) const 
  { 
    return cursor == x.cursor; 
  }

  //! Test for inequality
  bool operator != (const CyclicIterator& x) const 
  {
    return ! (*this == x); 
  }
  
  //! Return current data
  T operator*() const 
  {
    return *cursor; 
  }

  //! Return current iterator
  Iterator operator->() const 
  {
    return cursor; 
  }

  //! Pre-increment operator
  CyclicIterator& operator++() 
  {
    ++cursor;
    if (cursor == end) cursor = begin;
    return *this;
  }

  //! Post-increment operator
  CyclicIterator operator++(int) 
  {
    CyclicIterator ring = *this;
    ++*this;
    return ring;
  }
  
  //! Pre-decrement operator
  CyclicIterator& operator--() 
  {
    if (cursor == begin) cursor = end;
    --cursor;
    return *this;
  }
  
  //! Post-decrement operator
  CyclicIterator operator--(int) 
  {
    CyclicIterator ring = *this;
    --*this; 
    return ring;
  }
};

#endif
