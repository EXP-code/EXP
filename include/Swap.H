#ifndef _Swap_H
#define _Swap_H

template <class T>
class Swap
{
private:

  int done;

  union {
    T n;
    unsigned char b[sizeof(T)];
  } u;

public:
  
  Swap() : done(0) {}
  Swap(T i) : done(0) {u.n=i;}

  T& in() {return u.n;}
  T& rin() {done=0; return u.n;}

  T out() {

    if (!done) {
      unsigned char t;
      short n = sizeof(T);
      for (int j=0; j<n/2; j++) {
	t = u.b[j];
	u.b[j] = u.b[n-1-j];
	u.b[n-1-j] = t;
      }
      done = 1;
    }

    return u.n;
  }

};

#endif
