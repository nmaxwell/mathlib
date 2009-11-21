#ifndef APPLY_HDAF_H
#define APPLY_HDAF_H

#include <mathlib/math/std_math.h>

#include <mathlib/math/hdaf/hdaf.h>




// applies the (m,sigma,diff_order) differentiating hdaf, m,diff_order can be any non-negative integers.
// specify a range, (eps_min,eps_max) for tyhe truncation error to fall within; (1E-8,1E-7) will probably be good.
// h is the 
// 
void apply_hdaf_reg(
	int m,
	double sigma,
	int diff_order,
	double h,
	double eps_min,
	double eps_max,
    double * in,
	double * & out,
    int n,
    int in_stride=1,
    int out_stride=1) ;

















#endif



