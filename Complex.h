// This may look like C code, but it is really -*- C++ -*-
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

  Real		   real() const;
  Real		   imag() const;

				   Complex();
				   Complex(const Complex& y);
				   Complex(const Real r, Real i=0);

				  ~Complex();

  Complex&		   operator =  (const Complex& y);

  Complex&		   operator += (const Complex& y);
  Complex&		   operator += (const Real y);
  Complex&		   operator -= (const Complex& y);
  Complex&		   operator -= (const Real y);
  Complex&		   operator *= (const Complex& y);
  Complex&		   operator *= (const Real y);

  Complex&		   operator /= (const Complex& y); 
  Complex&		   operator /= (const Real y); 

  void			   error(const char* msg) const;
};


inline Complex operator /  (const Complex& x, const Complex& y);
inline Complex operator /  (const Complex& x, const Real y);
inline Complex operator /  (const Real x, const Complex& y);
inline Complex exp(const Complex& x);

// non-inline functions

Complex	  cos(const Complex& x);
Complex	  sin(const Complex& x);

Complex	  cosh(const Complex& x);
Complex	  sinh(const Complex& x);

Complex	  log(const Complex& x);

Complex	  pow(const Complex& x, int p);
Complex	  pow(const Complex& x, const Complex& p);
Complex	  pow(const Complex& x, const Real y);
Complex	  sqrt(const Complex& x);
   
#ifdef __GNUC__
#define INLINE
#else
#define INLINE inline
#endif

INLINE istream&  operator >> (istream& s, Complex& x);
INLINE ostream&  operator << (ostream& s, const Complex& x);

// other functions defined as inlines

inline int operator == (const Complex& x, const Complex& y);
inline int operator == (const Complex& x, const Real y);
inline int operator != (const Complex& x, const Complex& y);
inline int operator != (const Complex& x, const Real y);

inline Complex operator - (const Complex& x);
inline Complex conj(const Complex& x);
inline Complex operator + (const Complex& x, const Complex& y);
inline Complex operator + (const Complex& x, const Real y);
inline Complex operator + (const Real x, const Complex& y);
inline Complex operator - (const Complex& x, const Complex& y);
inline Complex operator - (const Complex& x, const Real y);
inline Complex operator - (const Real x, const Complex& y);
inline Complex operator * (const Complex& x, const Complex& y);
inline Complex operator * (const Complex& x, const Real y);
inline Complex operator * (const Real x, const Complex& y);

inline Real real(const Complex& x);
inline Real imag(const Complex& x);
inline Real abs(const Complex& x);
inline Real norm(const Complex& x);
inline Real arg(const Complex& x);

inline Complex polar(const Real r, const Real t = 0.0);


// inline members

inline Real Complex::real() const { return re; }
inline Real Complex::imag() const { return im; }

inline Complex::Complex() {}
inline Complex::Complex(const Complex& y) :re(y.real()), im(y.imag()) {}
inline Complex::Complex(const Real r, const Real i) :re(r), im(i) {}

inline Complex::~Complex() {}

inline Complex&	 Complex::operator =  (const Complex& y) 
{ 
  re = y.real(); im = y.imag(); return *this; 
} 

inline Complex&	 Complex::operator += (const Complex& y)
{ 
  re += y.real();  im += y.imag(); return *this; 
}

inline Complex&	 Complex::operator += (const Real y)
{ 
  re += y; return *this; 
}

inline Complex&	 Complex::operator -= (const Complex& y)
{ 
  re -= y.real();  im -= y.imag(); return *this; 
}

inline Complex&	 Complex::operator -= (const Real y)
{ 
  re -= y; return *this; 
}

inline Complex&	 Complex::operator *= (const Complex& y)
{  
  Real r = re * y.real() - im * y.imag();
  im = re * y.imag() + im * y.real(); 
  re = r; 
  return *this; 
}

inline Complex&	 Complex::operator *= (const Real y)
{  
  re *=	 y; im *=  y; return *this; 
}


//	functions

inline int	operator == (const Complex& x, const Complex& y)
{
  return x.real() == y.real() && x.imag() == y.imag();
}

inline int	operator == (const Complex& x, const Real y)
{
  return x.imag() == 0.0 && x.real() == y;
}

inline int	operator != (const Complex& x, const Complex& y)
{
  return x.real() != y.real() || x.imag() != y.imag();
}

inline int	operator != (const Complex& x, const Real y)
{
  return x.imag() != 0.0 || x.real() != y;
}

inline Complex	operator - (const Complex& x)
{
  return Complex(-x.real(), -x.imag());
}

inline Complex	conj(const Complex& x)
{
  return Complex(x.real(), -x.imag());
}

inline Complex	operator + (const Complex& x, const Complex& y)
{
  return Complex(x.real() + y.real(), x.imag() + y.imag());
}

inline Complex	operator + (const Complex& x, const Real y)
{
  return Complex(x.real() + y, x.imag());
}

inline Complex	operator + (const Real x, const Complex& y)
{
  return Complex(x + y.real(), y.imag());
}

inline Complex	operator - (const Complex& x, const Complex& y)
{
  return Complex(x.real() - y.real(), x.imag() - y.imag());
}

inline Complex	operator - (const Complex& x, const Real y)
{
  return Complex(x.real() - y, x.imag());
}

inline Complex	operator - (const Real x, const Complex& y)
{
  return Complex(x - y.real(), -y.imag());
}

inline Complex	operator * (const Complex& x, const Complex& y)
{
  return Complex(x.real() * y.real() - x.imag() * y.imag(), 
				 x.real() * y.imag() + x.imag() * y.real());
}

inline Complex	operator * (const Complex& x, const Real y)
{
  return Complex(x.real() * y, x.imag() * y);
}

inline Complex	operator * (const Real x, const Complex& y)
{
  return Complex(x * y.real(), x * y.imag());
}

inline Real  real(const Complex& x)
{
  return x.real();
}

inline Real  imag(const Complex& x)
{
  return x.imag();
}

inline Real  norm(const Complex& x)
{
  return (x.real() * x.real() + x.imag() * x.imag());
}

inline Real  abs(const Complex& x)
{
  return sqrt(norm(x));
}

inline Real  arg(const Complex& x)
{
  return atan2(x.imag(), x.real());
}

inline Complex	polar(const Real r, const Real t)
{
  return Complex(r * cos(t), r * sin(t));
}

#endif
