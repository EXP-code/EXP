#ifndef Iterable_H
#define Iterable_H

#include <iterator> // For std::forward_iterator_tag
#include <cstddef>  // For std::ptrdiff_t

template <typename T>
struct Iterator 
{
    using iterator_category = std::forward_iterator_tag;
    using difference_type   = std::ptrdiff_t;
    using value_type        = T;
    using pointer           = T*;  // or also value_type*
    using reference         = T&;  // or also value_type&
};


template <typename T>
class Iterable: public Iterator<T>
{
private:
  T* ptr;
public:
  explicit Iterable(T* _ptr) : ptr(_ptr) {}
  Iterable& operator++() {ptr++; return *this;}
  Iterable operator++(int) {Iterable retval = *this; ++(*this); return retval;}
  bool operator==(Iterable other) const {return ptr == other.ptr;}
  bool operator!=(Iterable other) const {return !(*this == other);}
  T* operator*() const {return ptr;}
};


#endif
