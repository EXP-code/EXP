#include <iostream>
#include <NoForce.H>

NoForce::NoForce(Component* c0, const YAML::Node& conf) : PotAccel(c0, conf)
{}

NoForce::~NoForce() 
{}
