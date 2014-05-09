#include <iostream>
#include <Elastic.H>

int main(void)
{
  Elastic elastic;

  while (1) {
    unsigned short Z;
    double E;
    std::cout << "Z, E? ";
    std::cin >> Z;
    std::cin >> E;
    std::cout << elastic(Z, E) << std::endl;
  }

  return 0;
}
