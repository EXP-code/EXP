#ifndef _InitContainer_H
#define _InitContainer_H

#include <map>

//@{
/**
   Wrapper template for STL containers for initializing value for new
   key.

   GetWithDef returns a value and MapWithDef returns the container.

   The Icont class derives a helper class using the STL container as a
   base class.  This version seems the cleanest but may have slightly
   more overhead.

   This version is updated for C++11 variadic arguments.  Should work
   for generic associative containers including the optional
   comparator and allocator parameters.
*/

//! Non-invasive, constant rvalue return
template < template<class,class,class...>
	   class C, typename K, typename V, typename... Args >
V GetWithDef(const C<K,V,Args...>& m, K const& key, const V & defval)
{
  typename C<K,V,Args...>::const_iterator it = m.find( key );
  if (it == m.end()) return defval;
  return it->second;
}

//! Non-invasive, assignable, lvalue return
template < template<class,class,class...>
	   class C, typename K, typename V, typename... Args >
V& GetWithDef(C<K,V,Args...>& m, K const& key, const V & defval)
{
  typename C<K,V,Args...>::iterator it = m.find( key );
  if (it == m.end()) return (m[key] = defval);
  return it->second;
}

//! This does the same thing as the non-invasive templates above using
//! inheritance to wrap the container
template < template<class,class,class...>
	   class C, typename K, typename V, typename... Args >
class Icont : public C<K, V, Args...>
{
public:
  //! Default value is zero for simple arithmetic types, but may be
  //! specialized as necessary.  If it fails to compile, it probably
  //! needs to be specialized.
  virtual V Default() { return 0; }

  //! Initializer
  V& operator[](const K& key)
  {
    std::map<K, V, Args...>* p = static_cast<std::map<K,V, Args...>*>(this);
    typename std::map<K, V>::iterator it = p->find(key);
    if (it == p->end()) {
      typename std::pair<const K, V> value(key, Default());
      it = p->insert(value).first; // Insert and return new iterator
    }
    return it->second;
  }
};

//@}

#endif
