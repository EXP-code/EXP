#include <iostream>

#include <Elastic.H>

int main(void)
{
  Elastic elastic;
  Geometric geometric;

  while (1) {
    unsigned short Z;
    double E;
    std::cout << "Z, E? ";
    std::cin >> Z;
    std::cin >> E;
    std::cout << "Elastic   = " << elastic(Z, E) << std::endl;
    std::cout << "Geometric = " << geometric(Z)  << std::endl;
  }

  return 0;
}
