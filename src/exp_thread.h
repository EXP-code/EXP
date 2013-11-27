
#if defined(__cplusplus)
extern "C" {
#endif

#include <pthread.h>

void exp_thread_fork(void *(*func) (void *), char * caller);
void make_mutex(pthread_mutex_t *m, char *caller, char *name);
void kill_mutex(pthread_mutex_t *m, char *caller, char *name);

#if defined(__cplusplus)
}
#endif

