
/* 
macros to help with portability between ANSI/traditional C, and C/C++.
*/


#ifdef __STDC__
#	define __(bob) bob
#else
#	define __(bob) ()
#endif




#ifdef __cplusplus
#	define LINK "C"
#else
#	define LINK 
#endif








