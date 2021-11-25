// Library variables

char threading_on = 0;
pthread_mutex_t mem_lock;
std::string outdir, runtag;
int nthrds = 1;
int this_step = 0;
unsigned multistep = 0;
unsigned maxlev = 100;
int mstep = 1;
int Mstep = 1;
int myid = 0;
std::mt19937 random_gen;
