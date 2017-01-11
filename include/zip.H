#include <iterator>
#include <boost/iterator/zip_iterator.hpp>

template <typename C1, typename C2>
class zip_container
{
  C1* c1;
  C2* c2;

  typedef boost::tuple<decltype(std::begin(*c1)),
		       decltype(std::begin(*c2))> Tuple;

public:
  zip_container(C1& c1, C2& c2) : c1(&c1), c2(&c2) {}

  typedef boost::zip_iterator<Tuple> iterator;

  iterator begin()// const
  {
    return iterator(Tuple(std::begin(*c1), std::begin(*c2)));
  }

  iterator end()// const
  {
    return iterator(Tuple(std::end(*c1), std::end(*c2)));
  }
};

template <typename C1, typename C2>
zip_container<C1, C2> zip(C1& c1, C2& c2)
{
  return zip_container<C1, C2>(c1, c2);
}
