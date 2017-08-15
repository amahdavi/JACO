/* ztompc.c
 * 
 * Copyright (C) 2007-2012 Andisheh Mahdavi
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
#include <string.h>
#include <math.h>
#include <gsl/gsl_integration.h>

struct ztompc_info {
  
   double omega_lambda,omega_curve,omega_matter;
};

static double oneovereofz(double z, void *p)
{
  double opz;
  struct ztompc_info *par = (struct ztompc_info *)p;

  opz = 1.+z;

  return 1./sqrt((par->omega_lambda + 
		  opz*opz*(par->omega_curve + 
			   opz*par->omega_matter)));
}

static double comoving_distance(double oM, double oL, double *oC, double z)
{
  gsl_function oozfunc;
  double answer,acc;
  struct ztompc_info par;
  static gsl_integration_workspace *ztompc_workspace = NULL;

  if (ztompc_workspace == NULL)
    ztompc_workspace = gsl_integration_workspace_alloc(1000);

  oozfunc.function = oneovereofz;
  oozfunc.params = &par;

  par.omega_matter = oM;
  if (oL < 0.) {
    par.omega_curve = 0.;
    par.omega_lambda = 1.-par.omega_matter;
  }
  else {
    par.omega_lambda = oL;
    par.omega_curve = 1.-(par.omega_matter+par.omega_lambda);
  }

  gsl_integration_qag(&oozfunc,0,z,0,1e-6,1000,GSL_INTEG_GAUSS21,
		      ztompc_workspace,&answer,&acc);
  *oC = par.omega_curve;
  return answer;
}

// Returns the transverse comoving distance ( = "angular size distance" )
// DIFFERENT FROM  the "angular diameter distance" (= TCD/(1+z))
double ztoMpc(double h0, double oM, double oL, double z)
{
  
  double comoving, hubble, sqrtk, answer,oC;

  comoving = comoving_distance(oM,oL,&oC,z);
  hubble = 299792.458/h0;

  if (oL < 0. || fabs(oC) < 1.e-6)
    answer = hubble*comoving;

  else if (oC > 0.) {
    sqrtk = sqrt(oC);
    answer = hubble*sinh(sqrtk*comoving)/sqrtk;
  } else {
    sqrtk = sqrt(-oC);
    answer = hubble*sin(sqrtk*comoving)/sqrtk;
  }

  return answer;
  
}

/* int main(int argc, char *argv[]) */
/* { */
/*   printf("%f\n",ztoMpc(atof(argv[1]),atof(argv[2]),atof(argv[3]), */
/* 		       atof(argv[4]))); */
/*   return 0.; */
/* } */
