#include "kernel.h"

int ApproximationCompare(const void *a, const void *b)
{
	return Problem->ApproximationTable->DefaultCompare(a,b);
}

int ApproximationKeyCompare(const void *key, const void *p, const size_t n)
{
	return Problem->ApproximationTable->DefaultKeyCompare(key,p,n);
}
