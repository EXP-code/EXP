
#include <iostream>
#include <any_thread.h>

class my_thread : public any_thread
{
private:
  void thread_call(void* pid);
public:
};

void my_thread::thread_call(void* pid)
{
  int id = *((int*)pid);
  cout << "It works! [" << id << "]\n";
}

int main()
{
  my_thread p;
  
  p.start_threads(20);
  p.join_threads();
}
