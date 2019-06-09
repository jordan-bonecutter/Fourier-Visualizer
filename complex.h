/* * * * * * * * * * * * * * * * * * * * * * * */
/* complex.h * * * * * * * * * * * * * * * * * */
/* Created by: Jordan Bonecutter * * * * * * * */
/* 17 April 2019 * * * * * * * * * * * * * * * */
/* * * * * * * * * * * * * * * * * * * * * * * */

#ifndef __COMPLEX_H__
#define __COMPLEX_H__

/* complex structure */
typedef struct {double re, im; } comp;

/* functions */
comp c_exp(double im);
comp c_mul(double re, comp c);
comp c_mulc(comp c1, comp c2);
comp c_divc(comp c1, comp c2);
comp c_add(comp c1, comp c2);
comp c_con(comp c);
double c_mag(comp c);

#endif

/* eof */
