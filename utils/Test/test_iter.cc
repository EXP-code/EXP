#include <iostream>
#include <iomanip>
#include <vector>

#include "Iterable.H"

template <typename T>
class Numbers
{
private:
  std::vector<T> array;

  T *start_; 
  T *end_;
public:
  Numbers(std::vector<T>& data)
  {
    start_ = &data[0];
    end_   = &data[data.size()];
  }
  Iterable<T> begin() { return Iterable<T>(start_); }
  Iterable<T> end()   { return Iterable<T>(end_); }
};


int main()
{
  std::vector<int> data = {3, 5, 7, 9, 11};

  for (auto n : Numbers<int>(data)) std::cout << *n << ",";
  std::cout << std::endl;

  Numbers<int> nums(data);
  for (auto v : nums) {
    std::cout << *v << ",";
  }
  std::cout << std::endl;
  
  return(0);
}
