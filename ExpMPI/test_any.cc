extern "C"
void *
call_any_threads_thread_call(void *atp)
{
  any_thread *p = (any_thread *)atp;
  p -> thread_call();
  return p;
}
                   
any_thread::start_thread()
{
  pthread_create(&t, 0, call_any_threads_thread_call, this);
}

void
thread_call()
{
  /* this code represents the thread */
}
