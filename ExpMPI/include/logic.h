
#ifndef _Logic_h
#define _Logic_h

#include <string>

using namespace std;

#ifdef IS_MAIN
string StrLogic[] = {"False", "True"};
#else
extern string *StrLogic;
#endif

#define ITOL(A) ( (A) == 0 ? FALSE : TRUE )

#endif
