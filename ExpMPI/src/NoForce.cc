#include <iostream>
#include <NoForce.H>

NoForce::NoForce(std::string& line) : PotAccel(line)
{
  cout << "NoForce created, line=" << line << "\n";
}

NoForce::~NoForce() 
{}
