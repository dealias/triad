#include <stdlib.h>
#include <iostream.h>
#include <new.h>

static size_t dynamic_memory=0;

void my_new_handler()
{
	cout << "Memory limits exceeded" << endl;
	exit(1);
}

// Total amount of dynamically allocated (including subsequently freed) memory.
size_t memory()
{
	return dynamic_memory;
}

void (*old_new_handler)()=set_new_handler(&my_new_handler);

void *operator new(size_t size)
{
	void *mem=malloc(size);
	if(size && !mem) (my_new_handler)();
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

void operator delete(void *ptr)
{
	free(ptr);
}

// provide a C++ interface to vector-resize via realloc 
void *operator new(size_t size, void *ptr, int new_len)
{
	size_t new_size=new_len*size;
	dynamic_memory += new_size;
	void *mem=realloc(ptr, new_size);
	if(new_size && !mem) (my_new_handler)();
	return mem;
}

