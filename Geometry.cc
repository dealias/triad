#include "NWave.h"

NWave *GeometryProblem;
GeometryBase *Geometry;

DynVector<Weight> *weight_ptr;

int GeometryCompare(const void *a, const void *b)
{
	return GeometryProblem->GeometryTable->DefaultCompare(a,b);
}

int GeometryKeyCompare(const void *key, const void *p, const size_t n)
{
	return GeometryProblem->GeometryTable->DefaultKeyCompare(key,p,n);
}
