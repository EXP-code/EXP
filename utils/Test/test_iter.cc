#include <iostream>
#include <iomanip>
#include <vector>

#include "Iterable.H"

class Numbers
{
private:
  const int start_;
  const int end_;
public:
  Numbers(int start, int end) : start_(start) , end_(end) {}
  myit begin() { return myit(start_); }
  myit end()   { return myit(end_); }
};


int main()
{
  // Numbers will be 3, 4
  for (auto n : Numbers(3, 5)) std::cout << n << ",";
  std::cout << std::endl;

  // Fills vec with 7, 8, 9
  Numbers nums(7, 10);
  std::vector<int> vec{std::begin(nums), std::end(nums)};
  for (auto v : vec) {
    std::cout << v << ",";
  }
  std::cout << std::endl;
  
  return(0);
}
