#ifndef __xstream_h__
#define __xstream_h__ 1

#include <stdio.h>
#include <iostream.h>
#include <rpc/rpc.h>

class xios : virtual public ios {
	iostate state;
public:
	int good() const { return state == 0; }
	int eof() const { return state & eofbit; }
	int fail() const { return state & (badbit|failbit); }
	int bad() const { return state & badbit; }
	void clear(iostate _state = 0) {state=_state;}
	void set(iostate flag) {state |= flag;}
    operator void*() const { return fail() ? (void*)0 : (void*)(-1); }
    int operator!() const { return fail(); }
};

class xstream : virtual public xios {
protected:
	FILE *buf;
	XDR xdrs;
public:
	void xopen(const char *filename, char *mode, xdr_op xop) {
		clear();
		buf=fopen(filename,mode);
		if(buf) xdrstdio_create(&xdrs, buf, xop);
		else set(badbit);
	}
	void close() {if(buf) fclose(buf); buf=NULL;}
};

#define IXSTREAM(T,N) ixstream& operator >> (T& x) \
{if(!xdr_##N(&xdrs, &x)) set(eofbit); return *this;}

#define OXSTREAM(T,N) oxstream& operator << (T x) \
{if(!xdr_##N(&xdrs, &x)) set(badbit); return *this;}

class ixstream : public xstream {
public:
    void open(const char *filename, openmode=in) {
		xopen(filename,"r",XDR_DECODE);
	}
	
	ixstream() {}
	ixstream(const char *filename) {open(filename);}
	ixstream(const char *filename, openmode mode) {open(filename,mode);}
	~ixstream() {close();}
	
	typedef ixstream& (*imanip)(ixstream&);
    ixstream& operator << (imanip func) { return (*func)(*this); }
	
	IXSTREAM(int,int)
	IXSTREAM(unsigned int,u_int)
	IXSTREAM(long,long)
	IXSTREAM(unsigned long,u_long)
	IXSTREAM(short,short)
	IXSTREAM(unsigned short,u_short)
	IXSTREAM(char,char)
	IXSTREAM(unsigned char,u_char)
	IXSTREAM(float,float)
	IXSTREAM(double,double)
};

class oxstream : public xstream {
public:
    void open(const char *filename, openmode mode=trunc) {
		char *smode;
		switch(mode) {
		case app:
			smode="a"; break;
		case trunc:
			smode="w"; break;
		}
		xopen(filename,smode,XDR_ENCODE);
	}
	
	oxstream() {}
	oxstream(const char *filename) {open(filename);}
	oxstream(const char *filename, openmode mode) {open(filename,mode);}
	~oxstream() {close();}

	oxstream& flush() {if(buf) fflush(buf); return *this;}
	
	typedef oxstream& (*omanip)(oxstream&);
    oxstream& operator << (omanip func) { return (*func)(*this); }
	
	OXSTREAM(int,int)
	OXSTREAM(unsigned int,u_int)
	OXSTREAM(long,long)
	OXSTREAM(unsigned long,u_long)
	OXSTREAM(short,short)
	OXSTREAM(unsigned short,u_short)
	OXSTREAM(char,char)
	OXSTREAM(unsigned char,u_char)
	OXSTREAM(float,float)
	OXSTREAM(double,double)
};

inline oxstream& endl(oxstream& s) {return s;}
inline oxstream& flush(oxstream& s) {return s; s.flush();}

#endif
