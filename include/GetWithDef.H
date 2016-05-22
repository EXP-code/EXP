#ifndef _GetWithDef_H
#define _GetWithDef_H

#include <map>

/**
   Wrapper template for STL containers for initializing value for new
   key.

   GetWithDef returns a value and MapWithDef returns the container.

   This version is updated for C++11 variadic arguments.  Should work
   for generic associative containers including the optional
   comparator and allocator parameters.
*/
template < template<class,class,class...>
	   class C, typename K, typename V, typename... Args >
V GetWithDef(const C<K,V,Args...>& m, K const& key, const V & defval)
{
  typename C<K,V,Args...>::const_iterator it = m.find( key );
  if (it == m.end()) return defval;
  return it->second;
}

template < template<class,class,class...>
	   class C, typename K, typename V, typename... Args >
V& GetWithDef(C<K,V,Args...>& m, K const& key, const V & defval)
{
  typename C<K,V,Args...>::iterator it = m.find( key );
  if (it == m.end()) return (m[key] = defval);
  return it->second;
}

template <typename K, typename V>
class Imap : public std::map<K, V>
{
public:
  //! Default value is zero for simple arithmetic types, but may be
  //! specialized as necessary.  If it fails to compile, it probably
  //! needs to be specialized.
  virtual V Default() { return 0; }

  //! Initializer
  V& operator[](const K& key)
  {
    std::map<K, V>* p = static_cast<std::map<K,V>*>(this);
    typename std::map<K, V>::iterator it = p->find(key);
    if (it == p->end()) {
      typename std::pair<const K, V> value(key, Default());
      it = p->insert(value).first; // Insert and return new iterator
    }
    return it->second;
  }
};

#endif
