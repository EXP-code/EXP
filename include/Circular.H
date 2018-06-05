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
  
  RingIterator (Container& v) : data(v), cursor(v.begin()), begin(v.begin()), end(v.end()) {}
  
  RingIterator (Container& v, const Iterator& i) : data(v), cursor(i), begin(v.begin()), end(v.end()) {}

  RingIterator (Container& v, const Iterator& i, const Iterator& j) : data(v), cursor(i), begin(i), end(j) {}

  RingIterator (Container& v, size_t i) : data(v), cursor(v.begin() + i % v.size()), begin(v.begin()), end(v.end()) {}

  bool operator == (const RingIterator& x) const 
  { 
    return cursor == x.cursor; 
  }
  
  bool operator != (const RingIterator& x) const 
  {
    return ! (*this == x); 
  }
  
  T operator*() const 
  {
    return *cursor; 
  }
  
  Iterator operator->() const 
  {
    return cursor; 
  }
  
  RingIterator& operator++() 
  {
    ++cursor;
    if (cursor == end)
      cursor = begin;
    return *this;
  }
  
  RingIterator operator++(int) 
  {
    RingIterator ring = *this;
    ++*this;
    return ring;
  }

  RingIterator& operator--() 
  {
    if (cursor == begin)
      cursor = end;
    --cursor;
    return *this;
  }

  RingIterator operator--(int) 
  {
    RingIterator ring = *this;
    --*this; 
    return ring;
  }

  RingIterator insert (const T& x)
  {
    return RingIterator (data, data.insert (cursor, x));
  }

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

  CyclicIterator (const Iterator& i, const Iterator& j) : cursor(i), begin(i), end(j) {}
  
  bool operator == (const CyclicIterator& x) const 
  { 
    return cursor == x.cursor; 
  }

  bool operator != (const CyclicIterator& x) const 
  {
    return ! (*this == x); 
  }
  
  T operator*() const 
  {
    return *cursor; 
  }

  Iterator operator->() const 
  {
    return cursor; 
  }

  CyclicIterator& operator++() 
  {
    ++cursor;
    if (cursor == end)
      cursor = begin;
    return *this;
  }

  CyclicIterator operator++(int) 
  {
    CyclicIterator ring = *this;
    ++*this;
    return ring;
  }
  
  CyclicIterator& operator--() 
  {
    if (cursor == begin)
      cursor = end;
    --cursor;
    return *this;
  }
  
  CyclicIterator operator--(int) 
  {
    CyclicIterator ring = *this;
    --*this; 
    return ring;
  }
};

#endif
