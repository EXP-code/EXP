
#ifndef _Logic_h
#ifdef __GNUG__
#pragma interface
#endif
#define _Logic_h 1

#include <string>

#ifdef IS_MAIN
String StrLogic[] = {"False", "True"};
#else
extern string *StrLogic;
#endif

#define ITOL(A) ( (A) == 0 ? FALSE : TRUE )

#endif
