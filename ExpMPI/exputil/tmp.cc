
using namespace std;

#include <iostream>
#include <iomanip>
#include <sstream>

void testme(istream* ins)
{
  const int nbuf = 1024;
  char line[nbuf];
  int icnt=0;
  while (!ins->eof()) {
    ins->getline(line, nbuf);
    cout << setw(5) << ++icnt << ": " << line << endl;
  }
}

int
main()
{
  string data("This is the first test!\n");
  data += "This is the second test!\n";
  data += "This is the third test!";

  istringstream ins(data);
  testme(&ins);
}

