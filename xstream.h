#ifndef __xstream_h__
#define __xstream_h__ 1

#define _ALL_SOURCE 1
#ifdef _POSIX_SOURCE
#undef _POSIX_SOURCE
#include <rpc/rpc.h>
#define _POSIX_SOURCE
#else
#include <rpc/rpc.h>
#endif

#include <stdio.h>

class xios {
public:
    enum io_state {goodbit=0, eofbit=1, failbit=2, badbit=4};
    enum open_mode {in=1, out=2, app=8, trunc=16};
private:	
	int _state;
public:	
	int good() const { return _state == 0; }
	int eof() const { return _state & eofbit; }
	int fail() const { return _state & (badbit|failbit); }
	int bad() const { return _state & badbit; }
	void clear(int state = 0) {_state=state;}
	void set(int flag) {_state |= flag;}
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
	void close() {
		if(buf) {
#if !(_CRAY || (__i386__ && !__ELF__))
			xdr_destroy(&xdrs);
#endif			
			fclose(buf);
			buf=NULL;
		}
	}
	void precision(int) {}
};

#define IXSTREAM(T,N) ixstream& operator >> (T& x) \
{if(!xdr_##N(&xdrs, &x)) set(eofbit); return *this;}

#define OXSTREAM(T,N) oxstream& operator << (T x) \
{if(!xdr_##N(&xdrs, &x)) set(badbit); return *this;}

class ixstream : public xstream {
public:
    void open(const char *filename, open_mode=in) {
		xopen(filename,"r",XDR_DECODE);
	}
	
	ixstream() {}
	ixstream(const char *filename) {open(filename);}
	ixstream(const char *filename, open_mode mode) {open(filename,mode);}
	~ixstream() {close();}
	
	typedef ixstream& (*imanip)(ixstream&);
    ixstream& operator << (imanip func) { return (*func)(*this); }
	
	IXSTREAM(int,int);
	IXSTREAM(unsigned int,u_int);
	IXSTREAM(long,long);
	IXSTREAM(unsigned long,u_long);
	IXSTREAM(short,short);
	IXSTREAM(unsigned short,u_short);
	IXSTREAM(char,char);
#ifndef _CRAY		
	IXSTREAM(unsigned char,u_char);
#endif		
	IXSTREAM(float,float);
	IXSTREAM(double,double);
};

class oxstream : public xstream {
public:
    void open(const char *filename, open_mode mode=trunc) {
		xopen(filename,(mode & app) ? "a" : "w",XDR_ENCODE);
	}
	
	oxstream() {}
	oxstream(const char *filename) {open(filename);}
	oxstream(const char *filename, open_mode mode) {open(filename,mode);}
	~oxstream() {close();}

	oxstream& flush() {if(buf) fflush(buf); return *this;}
	
	typedef oxstream& (*omanip)(oxstream&);
    oxstream& operator << (omanip func) { return (*func)(*this); }
	
	OXSTREAM(int,int);
	OXSTREAM(unsigned int,u_int);
	OXSTREAM(long,long);
	OXSTREAM(unsigned long,u_long);
	OXSTREAM(short,short);
	OXSTREAM(unsigned short,u_short);
	OXSTREAM(char,char);
#ifndef _CRAY		
	OXSTREAM(unsigned char,u_char);
#endif		
	OXSTREAM(float,float);
	OXSTREAM(double,double);
};

inline oxstream& endl(oxstream& s) {s.flush(); return s;}
inline oxstream& flush(oxstream& s) {s.flush(); return s;}

#undef IXSTREAM
#undef OXSTREAM

#endif
