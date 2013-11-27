
using namespace std;

#include <iostream>
#include <iomanip>
#include <sstream>
#include <vector>


void set_element(int *i, int j)
{
  *i = j;
}

int
main()
{
  vector<int> t(3, 3);

  cout << t[2] << endl;
  int *i = &t[2];
  set_element(i, 18);
  cout << t[2] << endl;
}

