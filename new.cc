#include <cstdlib>
#include <iostream>
#include <new>
//#include <fftw3.h>

// #include <mem_test_user>

using std::cout;
using std::endl;
using std::set_new_handler;

#if defined(__linux__) && defined (__i386__)
#include <fpu_control.h>
// Setup FPU for exceptions on overflow, zero divide and NaN.
fpu_control_t __fpu_control=
_FPU_EXTENDED | _FPU_RC_NEAREST | _FPU_MASK_DM | _FPU_MASK_UM | _FPU_MASK_PM;
#endif

static size_t dynamic_memory=0;

void my_new_handler()
{
  cout << endl << "Memory limits exceeded" << endl;
  exit(1);
}

void (*old_new_handler)()=set_new_handler(&my_new_handler);

// For Compaq cxx compiler:
#ifdef _RWSTD_THROW_SPEC_NULL
void *operator new(size_t size, const std::nothrow_t&) _RWSTD_THROW_SPEC_NULL
#else
void *operator new(size_t size)
#endif
{
  void *mem=malloc(size);
//  void *mem=fftw_malloc(size);
  if(size && !mem) (my_new_handler)();
  dynamic_memory += size;
  return mem;
}

void *operator new(size_t size, int len)
{
  void *mem=calloc(len,size);
//  void *mem=fftw_malloc(len*size);
//  bzero(mem,len*size);
  if(len && !mem) (my_new_handler)();
  dynamic_memory += len*size;
  return mem;
}

#ifdef _RWSTD_THROW_SPEC_NULL
void operator delete(void *ptr) _RWSTD_THROW_SPEC_NULL
#else
void operator delete(void *ptr)
#endif
{
  free(ptr);
//  fftw_free(ptr);
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
