#include "options.h"
#include "kernel.h"
#include "Matrix.h"

class MatExp : public ProblemBase {
private:	
	const Array2<Var> A;
public:
	MatExp(const Array2<Var> &a) : A(a) {
		ny=a.Nx();
		assert(ny == A.Ny());
	}
 	virtual ~MatExp() {}
	void InitialConditions() {y=new Var[ny];} 
	void Source(Var *src, Var *Y, double) {
		for(unsigned int i=0; i < ny; i++) {
			Array1<Var> Ai=A[i];
			Var sum=0.0;
			for(unsigned int j=0; j < ny; j++) sum += Ai[j]*Y[j];
			src[i]=sum;
		}
	}
};

// Compute the matrix exponential of a square matrix
template<class T>
const Array2<T>& exp(const Array2<T>& A)
{
//	double tolmax=1e-7;
//	double tolmin=9e-8;
	
	double tolmax=1e-4;
	double tolmin=9e-5;
	MatExp MatExpProblem(A);

	MatExpProblem.InitialConditions();
	Var *y=MatExpProblem.Vector();
	unsigned int n=MatExpProblem.Size();
	
	unsigned int n2=n;
	static DynVector<T> temp(n2);
	if(n2 > temp.Alloc()) temp.Resize(n2);
	static Array2<T> B(n,n,temp);
	
	static RK5 Int;
	Int.Allocate(n);
	Int.SetParam(tolmax,tolmin,2.0,2.0,0.0,DBL_STD_MAX,INT_MAX,100,0,1);
		
	for(unsigned int j=0; j < n; j++) {
		for(unsigned int i=0; i < j; i++) y[i]=0.0;
		y[j]=1.0;
		for(unsigned int i=j+1; i < n; i++) y[i]=0.0;
		
		Real dt=1.0e-5;
		Real t=0.0;
		int iteration=0;
		Int.Integrate(MatExpProblem,y,t,1.0,dt,-1.0,iteration);
		for(unsigned int i=0; i < n; i++) B(i,j)=y[i];
	}
	
	A.Purge();
	return B;
}	
