
#include <math.h>
#include <iostream.h>

#include <numerical.h>
#include <Vector.h>
#include <interp.h>
#include <isothermal.h>

main()
{
  cout << "This is a test!" << endl;
  ACG gen (11, 20);
  Uniform u(0.0, 1.0, &gen);

  for (int i=0; i<10; i++) cout << u() << endl;
}

