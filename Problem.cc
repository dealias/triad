#include "kernel.h"

int ProblemCompare(const void *a, const void *b)
{
	return Vocabulary->ProblemTable->DefaultCompare(a,b);
}

int ProblemKeyCompare(const void *key, const void *p, const size_t n)
{
	return Vocabulary->ProblemTable->DefaultKeyCompare(key,p,n);
}
