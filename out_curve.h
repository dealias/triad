// Cfront can't seem to handle a template here.
inline void out_function(ostream& os, Real (*f)(int), char *text, int n,
						 int nperline)
{
	int i;
	os << "# " << text << endl;
	if(nperline == 4) { // This case is optimized.
		int ncont=((n-1)/nperline)*nperline;
		for(i=0; i < ncont; i += nperline)
			os << (*f)(i) << "\t" << (*f)(i+1) << "\t" << (*f)(i+2) << "\t"
			   << (*f)(i+3) << " \\" << endl;
		for(i=ncont; i < n-1; i++) os << (*f)(i) << "\t";
	} else {
		for(i=0; i < n-1;) {
			os << (*f)(i);
			if(++i % nperline) os << "\t"; else os << " \\" << endl;
		}
	}
	os << (*f)(n-1) << endl;
}

template<class T>
void out_curve(ostream& os, T *f, char *text, int n, int nperline)
{
	int i;
	os << "# " << text << endl;
	if(nperline == 4) { // This case is optimized.
		int ncont=((n-1)/nperline)*nperline;
		for(i=0; i < ncont; i += nperline)
			os << f[i] << "\t" << f[i+1] << "\t" << f[i+2] << "\t"
			   << f[i+3] << " \\" << endl;
		for(i=ncont; i < n-1; i++) os << f[i] << "\t";
	} else {
		for(i=0; i < n-1;) {
			os << f[i];
			if(++i % nperline) os << "\t"; else os << " \\" << endl;
		}
	}
	os << f[n-1] << endl;
}

