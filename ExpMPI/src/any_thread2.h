                   
#include <pthread.h>

class any_thread;

struct thrd_pass {
  any_thread *t;
  int id;
};

class any_thread
{
 private:
  int nthreads;

 public:
  thrd_pass *td;
  pthread_t *t;

  any_thread(void) : t(NULL), td(NULL) {}
  ~any_thread(void);
  void start_threads(int);
				// this function will be used for thread
				// execution
  virtual void thread_call(void*) = 0;
};


