
using namespace std;

#include <iostream>

void testme(ostream& out)
{
  out << "This is a test!\n";
  out << "This is a test!\n";
  out << "This is a test!\n";
}

int
main()
{
  testme(cout);
}

