/* Inline C++ integer exponentiation routines 
   Version 1.01
   Copyright (C) 1999-2004 John C. Bowman <bowman@math.ualberta.ca>

This program is free software; you can redistribute it and/or modify
it under the terms of the GNU General Public License as published by
the Free Software Foundation; either version 2 of the License, or
(at your option) any later version.

This program is distributed in the hope that it will be useful,
but WITHOUT ANY WARRANTY; without even the implied warranty of
MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
GNU General Public License for more details.

You should have received a copy of the GNU General Public License
along with this program; if not, write to the Free Software
Foundation, Inc., 675 Mass Ave, Cambridge, MA 02139, USA. */

#ifndef __pow_h__
#define __pow_h__ 1

#include <cassert>

inline double pow(double x, int p)
{
  if(p == 0) return 1.0;
  if(x == 0.0 && p > 0) return 0.0;
  if(p < 0) {p=-p; x=1/x;}
	
  double r = 1.0;
  for(;;) {
    if(p & 1) r *= x;
    if((p >>= 1) == 0)	return r;
    x *= x;
  }
}

inline double pow(double x, size_t p)
{
  if(p == 0) return 1.0;
  if(x == 0.0) return 0.0;
	
  double r = 1.0;
  for(;;) {
    if(p & 1) r *= x;
    if((p >>= 1) == 0)	return r;
    x *= x;
  }
}

inline int pow(int x, int p)
{
  if(p == 0) return 1;
  if(x == 0 && p > 0) return 0;
  if(p < 0) {assert(x == 1 || x == -1); return (-p % 2) ? x : 1;}
	
  int r = 1;
  for(;;) {
    if(p & 1) r *= x;
    if((p >>= 1) == 0)	return r;
    x *= x;
  }
}

inline size_t pow(size_t x, size_t p)
{
  if(p == 0) return 1;
  if(x == 0) return 0;
	
  size_t r = 1;
  for(;;) {
    if(p & 1) r *= x;
    if((p >>= 1) == 0)	return r;
    x *= x;
  }
}

#endif
