
#ifndef _gknots_h

#define _gknots_h 1

struct GAUSS {
  int n;
  double *x,*w;
};

/*
#ifdef __GNUC__
void get_gknots(struct GAUSS *, int);
void free_gknots(struct GAUSS *, int);
#else
void get_gknots();
void free_gknots();
#endif
*/

void get_gknots(struct GAUSS *, int);
void free_gknots(struct GAUSS *, int);

#endif

