/* 
   Copyright (C) 1988 Free Software Foundation
   written by Doug Lea (dl@rocky.oswego.edu)

   This file is part of the GNU C++ Library.  This library is free
   software; you can redistribute it and/or modify it under the terms of
   the GNU Library General Public License as published by the Free
   Software Foundation; either version 2 of the License, or (at your
   option) any later version.	This library is distributed in the hope
   that it will be useful, but WITHOUT ANY WARRANTY; without even the
   implied warranty of MERCHANTABILITY or FITNESS FOR A PARTICULAR
   PURPOSE.  See the GNU Library General Public License for more details.
   You should have received a copy of the GNU Library General Public
   License along with this library; if not, write to the Free Software
   Foundation, 675 Mass Ave, Cambridge, MA 02139, USA.
*/

#ifndef __Complex_h__
#ifdef __GNUG__
#pragma interface
#endif
#define __Complex_h__ 1

#define __ATT_complex__

#include <iostream.h>
#include <math.h>
#include "precision.h"

class Complex
{
#ifdef __ATT_complex__
public:
#else
protected:
#endif

	Real		   re;
	Real		   im;

public:

	Complex() {}
	Complex(Real r, Real i=0) :re(r), im(i) {}
	Complex(const Complex& y) :re(y.re), im(y.im) {}
	
	~Complex() {}

	Real real(const Complex &) const {return re;}
	Real imag(const Complex &) const {return im;}

	Complex& operator = (const Complex& y);
	
	Complex& operator += (const Complex& y);
	Complex& operator += (Real y);
	Complex& operator -= (const Complex& y);
	Complex& operator -= (Real y);
	Complex& operator *= (const Complex& y);
	Complex& operator *= (Real y);
	Complex& operator /= (const Complex& y); 
	Complex& operator /= (Real y); 
	
	void error(char* msg) const;
};


// non-inline functions

Complex cos(const Complex& x);
Complex sin(const Complex& x);

Complex cosh(const Complex& x);
Complex sinh(const Complex& x);

Complex log(const Complex& x);

Complex pow(const Complex& x, int p);
Complex pow(const Complex& x, const Complex& p);
Complex pow(const Complex& x, Real y);
Complex sqrt(const Complex& x);
   
// other functions defined as inlines

inline int operator == (const Complex& x, const Complex& y);
inline int operator == (const Complex& x, Real y);
inline int operator != (const Complex& x, const Complex& y);
inline int operator != (const Complex& x, Real y);

inline Complex operator - (const Complex& x);
inline Complex operator + (const Complex& x, const Complex& y);
inline Complex operator + (const Complex& x, Real y);
inline Complex operator + (Real x, const Complex& y);
inline Complex operator - (const Complex& x, const Complex& y);
inline Complex operator - (const Complex& x, Real y);
inline Complex operator - (Real x, const Complex& y);
inline Complex operator * (const Complex& x, const Complex& y);
inline Complex operator * (const Complex& x, Real y);
inline Complex operator * (Real x, const Complex& y);
inline Complex operator / (const Complex& x, const Complex& y);
inline Complex operator / (const Complex& x, Real y);
inline Complex operator / (Real x, const Complex& y);

inline Complex conj(const Complex& x);
inline Real real(const Complex& x);
inline Real imag(const Complex& x);
inline Real abs(const Complex& x);
inline Real norm(const Complex& x);
inline Real arg(const Complex& x);

inline Complex polar(Real r, Real t = 0.0);


// inline members

inline Complex& Complex::operator = (const Complex& y) 
{ 
  re = y.re; im = y.im; return *this; 
} 

inline Complex& Complex::operator += (const Complex& y)
{ 
	re += y.re;  im += y.im; return *this; 
}

inline Complex& Complex::operator += (Real y)
{ 
	re += y; return *this; 
}

inline Complex& Complex::operator -= (const Complex& y)
{ 
	re -= y.re;  im -= y.im; return *this; 
}

inline Complex& Complex::operator -= (Real y)
{ 
	re -= y; return *this; 
}

inline Complex& Complex::operator *= (const Complex& y)
{  
	Real r = re * y.re - im * y.im;
	im = re * y.im + im * y.re; 
	re = r; 
	return *this; 
}

inline Complex& Complex::operator *= (Real y)
{  
	re *= y; im *= y; return *this; 
}

inline Complex& Complex::operator /= (const Complex& y)
{
	register double t1,t2,t3,t4;
	t3=y.re; t4=y.im; t1=t2=1.0/(t3*t3+t4*t4);
	t1*=t3; t2*=t4; t3=re; t4=im;
	re*=t1;	re+=t4*t2;
	im*=t1;	im-=t3*t2;
	return *this;
}

inline Complex& Complex::operator /= (Real y)
{
	re/=y;
	im/=y;
	return *this;
}

//	functions

inline int	operator == (const Complex& x, const Complex& y)
{
	return x.re == y.re && x.im == y.im;
}

inline int	operator == (const Complex& x, Real y)
{
	return x.im == 0.0 && x.re == y;
}

inline int	operator != (const Complex& x, const Complex& y)
{
	return x.re != y.re || x.im != y.im;
}

inline int	operator != (const Complex& x, Real y)
{
	return x.im != 0.0 || x.re != y;
}

inline Complex operator - (const Complex& x)
{
	return Complex(-x.re, -x.im);
}

inline Complex conj(const Complex& x)
{
	return Complex(x.re, -x.im);
}

inline Complex operator + (const Complex& x, const Complex& y)
{
	return Complex(x.re + y.re, x.im + y.im);
}

inline Complex operator + (const Complex& x, Real y)
{
	return Complex(x.re + y, x.im);
}

inline Complex operator + (Real x, const Complex& y)
{
	return Complex(x + y.re, y.im);
}

inline Complex operator - (const Complex& x, const Complex& y)
{
	return Complex(x.re - y.re, x.im - y.im);
}

inline Complex operator - (const Complex& x, Real y)
{
	return Complex(x.re - y, x.im);
}

inline Complex operator - (Real x, const Complex& y)
{
	return Complex(x - y.re, -y.im);
}

inline Complex operator * (const Complex& x, const Complex& y)
{
	return Complex(x.re * y.re - x.im * y.im, 
				   x.re * y.im + x.im * y.re);
}

inline Complex operator * (const Complex& x, Real y)
{
	return Complex(x.re * y, x.im * y);
}

inline Complex operator * (Real x, const Complex& y)
{
	return Complex(x * y.re, x * y.im);
}

inline Complex operator / (const Complex& x, const Complex& y)
{
	register double t1,t2,t3,t4;
	t3=y.re; t4=y.im; t1=t2=1.0/(t3*t3+t4*t4);
	t1*=t3; t2*=t4; t3=x.re; t4=x.im;
	return Complex(t4*t2+t3*t1,t4*t1-t3*t2);
}

inline Complex operator / (const Complex& x, Real y)
{
	return Complex(x.re/y,x.im/y);
}

inline Complex operator / (Real x, const Complex& y)
{
	register double factor;
	factor=1.0/(y.re*y.re+y.im*y.im);
	return Complex(x*y.re*factor,-x*y.im*factor);
}

inline Real real(const Complex& x)
{
	return x.re;
}

inline Real imag(const Complex& x)
{
	return x.im;
}

inline Real norm(const Complex& x)
{
	return (x.re * x.re + x.im * x.im);
}

inline Real abs(const Complex& x)
{
	return sqrt(norm(x));
}

inline Real arg(const Complex& x)
{
	return atan2(x.im, x.re);
}

inline Complex polar(Real r, Real t)
{
	return Complex(r * cos(t), r * sin(t));
}

#ifdef __GNUC__
istream&  operator >> (istream& s, Complex& x);
ostream&  operator << (ostream& s, const Complex& x);
#else
inline istream& operator >> (istream& s, Complex& y)
{
	Real re,im;
	s >> "(" >> re >> "," >> im >> ")";
	y=Complex(re,im);
	return s;
}

inline ostream& operator << (ostream& s, const Complex& y)
{
	s << "(" << y.re << ", " << y.im << ")";
	return s;
}
#endif

#endif
