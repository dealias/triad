#include <cstdlib>
#include <iostream>
#include <new>

// #include <mem_test_user>

using std::cerr;
using std::endl;
using std::set_new_handler;

static size_t dynamic_memory=0;

void my_new_handler()
{
  cerr << endl << "Memory limits exceeded" << endl;
  exit(1);
}

void (*old_new_handler)()=set_new_handler(&my_new_handler);

#ifdef __DECCXX_LIBCXX_RH70
void *operator new(size_t size, const std::nothrow_t&) _RWSTD_THROW_SPEC_NULL
#else
#define _RWSTD_THROW_SPEC_NULL
void *operator new(size_t size) throw (std::bad_alloc)
#endif
{
  void *mem=malloc(size);
  if(size && !mem) my_new_handler();
  dynamic_memory += size;
  return mem;
}

void *operator new(size_t size, int len)
{
  void *mem=calloc(len,size);
  if(len && !mem) (my_new_handler)();
  dynamic_memory += len*size;
  return mem;
}

void operator delete(void *ptr) throw ()//_RWSTD_THROW_SPEC_NULL
{
  free(ptr);
}

// provide a C++ interface to vector-resize via realloc 
// Warning: this does not initialize virtual member functions.
void *operator new(size_t size, void *ptr, int new_len)
{
  size_t new_size=new_len*size;
  dynamic_memory += new_size;
  void *mem=realloc(ptr, new_size);
  if(new_size && !mem) (my_new_handler)();
  return mem;
}

// Total amount of dynamically allocated (including subsequently freed) memory.
size_t memory()
{
  return dynamic_memory;
}
