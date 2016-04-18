// Compile string: g++ -std=c++11 -o testBS testBitset.cc

#include <iostream>
#include <iomanip>
#include "EnumBitset.H"

enum BitFlags
{
  False,
  True, 
  FileNotFound,
  Write,
  Read,
  MaxVal
};

template<>
struct EnumTraits<BitFlags>
{
  static const BitFlags max = BitFlags::MaxVal;
};

int main()
{
  EnumBitset<BitFlags> f;

  f.flip(BitFlags::True);
  f.flip(BitFlags::FileNotFound);

  // f.flip(2); // fails to compile, as expected
  
  std::cout << "Size of bitset: " << f.size() << std::endl << std::endl;
  std::cout << "Is False? "        << f.test(BitFlags::False) << std::endl;
  std::cout << "Is True? "         << f.test(BitFlags::True) << std::endl;
  std::cout << "Is FileNotFound? " << f.test(BitFlags::FileNotFound) << std::endl;
  std::cout << "Is Write? "        << f.test(BitFlags::Write) << std::endl;
  std::cout << "Is Read? "         << f.test(BitFlags::Read) << std::endl;
  std::cout << std::endl << "Operator version" << std::endl << std::endl;
  std::cout << "Is False? "        << f[BitFlags::False] << std::endl;
  std::cout << "Is True? "         << f[BitFlags::True] << std::endl;
  std::cout << "Is FileNotFound? " << f[BitFlags::FileNotFound] << std::endl;
  std::cout << "Is Write? "        << f[BitFlags::Write] << std::endl;
  std::cout << "Is Read? "         << f[BitFlags::Read] << std::endl;

  f.reset();

  std::cout << std::endl << "After full reset" << std::endl << std::endl;
  std::cout << "Is False? "        << f[BitFlags::False] << std::endl;
  std::cout << "Is True? "         << f[BitFlags::True] << std::endl;
  std::cout << "Is FileNotFound? " << f[BitFlags::FileNotFound] << std::endl;
  std::cout << "Is Write? "        << f[BitFlags::Write] << std::endl;
  std::cout << "Is Read? "         << f[BitFlags::Read] << std::endl;
}
