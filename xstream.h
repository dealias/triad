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
		else {errno=0; set(badbit);}
	}
	void close() {if(buf) fclose(buf);}
};

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
	
    ixstream& operator >> (char& x) {
		if(!xdr_char(&xdrs, &x)) set(badbit);
		return *this;
	}
    ixstream& operator >> (int& x) {
		if(!xdr_int(&xdrs, &x)) set(badbit);
		return *this;
	}
    ixstream& operator >> (float& x) {
		if(!xdr_float(&xdrs, &x)) set(badbit);
		return *this;
	}
    ixstream& operator >> (double& x) {
		if(!xdr_double(&xdrs, &x)) set(badbit);
		return *this;
	}
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

	typedef oxstream& (*omanip)(oxstream&);
    oxstream& operator << (omanip func) { return (*func)(*this); }
	
    oxstream& operator << (char x) {
		if(!xdr_char(&xdrs, &x)) set(badbit);
		return *this;
	}
    oxstream& operator << (int x) {
		if(!xdr_int(&xdrs, &x)) set(badbit);
		return *this;
	}
    oxstream& operator << (float x) {
		if(!xdr_float(&xdrs, &x)) set(badbit);
		return *this;
	}
    oxstream& operator << (double x) {
		if(!xdr_double(&xdrs, &x)) set(badbit);
		return *this;
	}
};

inline oxstream& endl(oxstream& s) {return s;}

#endif
