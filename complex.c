/* * * * * * * * * * * * * * * * * * * * * * * */
/* complex.c * * * * * * * * * * * * * * * * * */
/* Created by: Jordan Bonecutter * * * * * * * */
/* 17 April 2019 * * * * * * * * * * * * * * * */
/* * * * * * * * * * * * * * * * * * * * * * * */

#include "complex.h"
#include <math.h>

#define sqr(x)		((x)*(x))

double c_mag(comp c)
{
	/* Magnitude of complex number */
	return sqrt(sqr(c.re)+sqr(c.im));
}

comp c_con(comp c)
{
	/* z' = a - jb */
	comp ret;
	ret.re = c.re;
	ret.im = -c.im;
	return ret;
}

comp c_mulc(comp c1, comp c2)
{
	/* Complex multiply */
	comp ret;
	ret.re = c1.re*c2.re - c2.im*c1.im;
	ret.im = c1.im*c2.re + c2.im*c1.re;
	return ret;
}

comp c_divc(comp c1, comp c2)
{
	/* Complex denom */
	comp ret;
	double denom = sqr(c2.re)+sqr(c2.im);
	ret.re = (c1.re*c2.re + c1.im*c2.im)/denom;
	ret.im = (c1.im*c2.re - c1.re*c2.im)/denom;
	return ret;
}

comp c_add(comp c1, comp c2)
{
	/* Complex add */
	comp ret;
	ret.re = c1.re + c2.re;
	ret.im = c1.im + c2.im;
	return ret;
}

comp c_mul(double re, comp c)
{
	/* Complex multiplpy by a real number */
	comp ret;
	ret.im = c.im*re;
	ret.re = c.re*re;
	return ret;
}

comp c_exp(double im)
{
	/* exp(j*im) */
	comp ret;
	ret.re = cos(im);
	ret.im = sin(im);
	return ret;
}

/* eof */
