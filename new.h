#ifndef __new_h__
#define __new_h__ 1

#include <new.h>

void *operator new(size_t size, void *ptr, int new_len);

#endif

