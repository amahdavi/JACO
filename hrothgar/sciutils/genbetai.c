/* genbetai.c
 * 
 * Copyright (C) 2007-2012 Andisheh Mahdavi
 * The function beta_cont_frac was taken from specfunc/beta_inc.c
 * from the GNU Scientific Library v1.7
 * Copyright (C) 1996, 1997, 1998, 1999, 2000 Gerard Jungman
 * 
 * This program is free software; you can redistribute it and/or modify
 * it under the terms of the GNU General Public License as published by
 * the Free Software Foundation; either version 2 of the License, or (at
 * your option) any later version.
 * 
 * This program is distributed in the hope that it will be useful, but
 * WITHOUT ANY WARRANTY; without even the implied warranty of
 * MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the GNU
 * General Public License for more details.
 * 
 * You should have received a copy of the GNU General Public License
 * along with this program; if not, write to the Free Software
 * Foundation, Inc., 51 Franklin Street, Fifth Floor, 
 * Boston, MA 02110-1301, USA.
 */

#include <stdio.h>
#include <stdlib.h>
#include <gsl/gsl_math.h>
#include <gsl/gsl_sf_gamma.h>
#include <gsl/gsl_sf_psi.h>
#include <gsl/gsl_integration.h>
#include <math.h>
#include <sciutils.h>

#define MAXIT 100
#define EPS 1.0e-10
#define FPMIN 1.0e-30
#define EULER 0.5772156649

static
double
beta_cont_frac(
  const double a,
  const double b,
  const double x
  //,gsl_sf_result * result
  )
{
  const unsigned int max_iter = 512;        /* control iterations      */
  const double cutoff = 2.0 * GSL_DBL_MIN;  /* control the zero cutoff */
  unsigned int iter_count = 0;
  double cf;

  /* standard initialization for continued fraction */
  double num_term = 1.0;
  double den_term = 1.0 - (a+b)*x/(a+1.0);
  if (fabs(den_term) < cutoff) den_term = cutoff;
  den_term = 1.0/den_term;
  cf = den_term;

  while(iter_count < max_iter) {
    const int k  = iter_count + 1;
    double coeff = k*(b-k)*x/(((a-1.0)+2*k)*(a+2*k));
    double delta_frac;

    /* first step */
    den_term = 1.0 + coeff*den_term;
    num_term = 1.0 + coeff/num_term;
    if(fabs(den_term) < cutoff) den_term = cutoff;
    if(fabs(num_term) < cutoff) num_term = cutoff;
    den_term  = 1.0/den_term;

    delta_frac = den_term * num_term;
    cf *= delta_frac;

    coeff = -(a+k)*(a+b+k)*x/((a+2*k)*(a+2*k+1.0));

    /* second step */
    den_term = 1.0 + coeff*den_term;
    num_term = 1.0 + coeff/num_term;
    if(fabs(den_term) < cutoff) den_term = cutoff;
    if(fabs(num_term) < cutoff) num_term = cutoff;
    den_term = 1.0/den_term;

    delta_frac = den_term*num_term;
    cf *= delta_frac;

    if(fabs(delta_frac-1.0) < 2.0*GSL_DBL_EPSILON) break;

    ++iter_count;
  }

  return cf;
/*   result->val = cf; */
/*   result->err = iter_count * 4.0 * GSL_DBL_EPSILON * fabs(cf); */

/*   if(iter_count >= max_iter) */
/*     GSL_ERROR ("error", GSL_EMAXITER); */
/*   else */
/*     return GSL_SUCCESS; */
}

double genbeta(double a, double b)
{
  static double olda[3],oldb[3],oldbeta[3],answer;
  int i;

  for (i = 0; i < 3; ++i)
    if (fabs(a-olda[i]) < EPS && fabs(b-oldb[i]) < EPS) 
      return oldbeta[i];

  if (a < 0.1 && fabs(a-round(a)) < EPS) 
    printf("Integral a<=0 in genbeta: %f\n",a);
  if (b < -0.1 && fabs(b-round(b)) < EPS) 
    printf("Integral b<0 in genbeta: %f\n", b);

  if (a < 0) return genbeta(a+1.,b)*(1.+b/a);
  if (b < 0) return genbeta(b,a);

  answer = gsl_sf_beta(a,b);

  for (i = 2; i > 0; --i) {
    oldbeta[i] = oldbeta[i-1];
    oldb[i] = oldb[i-1];
    olda[i] = olda[i-1];
  }
  oldbeta[0] = answer;
  oldb[0] = b;
  olda[0] = a;
   
  return answer;
}


// Generalized COMPLEMENT of beta function. 
// where genbeta(a,b,0) = beta(a,b) and genbeta(a,b,1) = 0.
// Here a or b can be negative (unlike in the traditional version)
// But a cannot be zero.
// Other inputs: betaab gives beta(a,b) (this should check for poles)
//               omxbxa gives (1-x)^b x^a
//               i.e. omxbxa = exp(a*log(x)+b*log(1.-x));

double genbetaq(double a, double b, double x, 
		double betaab, double omxbxa, double psiofa)
{

  double fac;

  if (x < 0. || x > 1.) {
    printf("Error---invalid X in genbetaq\n");
    exit(-1);
  }
  
  if (a < 0.) {
    fac = 1.+b/a;
    return genbetaq(a+1.,b,x,betaab/fac,omxbxa*x,psiofa)*fac+omxbxa/a;
  }
  if (b < 0.) 
    // omxbxa is unaltered by switch in arguments!
	  return betaab-genbetaq(b,a,1.-x,betaab,omxbxa,psiofa);
  
  if (fabs(x) < 1.E-7) return 0.;
  if (fabs(x-1.0) < 1.E-7) return betaab;
  
  if (x < (a+1.)/(a+b+2.) || fabs(b) < EPS)
    return omxbxa*beta_cont_frac(a,b,x)/a;
  else 
    return betaab-omxbxa*beta_cont_frac(b,a,1.0-x)/b;
/*   else { */
/*     omx = 1.-x; */
/*     return (-log(omx)-psiofa+ */
/*       (a-1.)*omx*(2.*a*a*omx*omx+a*omx*(10*x-19)+ */
/* 		  6.*(11.+x*(2.*x-7.)))/36.); */
/*   } */
}


double genbetai(double a, double b, double x)
{
  double betaab=1.E30,omxbxa;
  static double psiofa=0.;
  
  if (fabs(b) > EPS) 
    betaab = genbeta(a,b);
  else {
    if (a < 1) {
      printf("Warning---inaccurate Beta(z,a,0) regime\n");
      return x;
    }
/*     if (fabs(a-olda) > EPS) { */
/*       psiofa = gsl_sf_psi(a); */
/*       olda = a;  */
/*     }  */
  }

  omxbxa = pow(1.-x,b)*pow(x,a);

  return genbetaq(a,b,x,betaab,omxbxa,psiofa);
}



/* int main(int argc, char *argv[]) */
/* { */
/*   double a,b,x; */

/*   a = atof(argv[1]); */
/*   b = atof(argv[2]); */
/*   x = atof(argv[3]); */
/*   printf("%.10f\n",genbetai(a,b,x)); */

/*   return 0.; */
/* } */


