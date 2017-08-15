/* veldisp.c
 * 
 * Copyright (C) 2009 Alison Mansheim, Andisheh Mahdavi
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
#include <math.h>
#include <gsl/gsl_integration.h>
#include <gsl/gsl_math.h>
#include <jaco.h>

//all int main should be one call velocity dispersion
//for cvs backup: cp dispersion_int7.c ../../velocitydispersion/
//gets r_nu, r_o and n value(s) from command line
//loops to create output file for each n, for a given beta function
//I removed all analytical functions in this version 
//compile with gcc -Wall -lm -lgsl -lgslcblas -o dispersion_int7 dispersion_int7.c
//full general integrator for line-of-sight velocity dispersion
//run with ./dispersion_int7 r_nu r_o n1 n2 n3 n4... 

// Profiles for testing:

//next step to define all functions within a structure, pass to JACO
/* struct f_params{  */
  
/*   double (*beta)(double r, void *p);  */
/*   double (*mass)(double r, void *p); */
/*   double (*nu)(double r, void *p); */
/*   double G_Mpc; double M_o; double r_nu; double nu_o; double r_o; double R; double n; */

/* }; */


/* //define mass distribution function */
/* double M_profile(double r, void * p) */
/* { */
/*   //cast pointer p as variable type jaco_state */
/*   //insert incoming p values into   */
/*   struct jaco_state * params = (struct jaco_state *)p; */

/*   double r_nu = (params->r_nu); */
/*   double M_o = (params->M_o); */
/*   double n = (params->n); */
/*   double x = pow(r, 3-n); //get rid of pow to increase efficiency? */
/*   double y = (r_nu + r); */
/*   double z = pow(y, 3-n); */

/*   return ( M_o * x ) / z; */
/* } */


/* //anisotropy parameter as a function of radius */
/* double beta_profile(double r, void *p) */
/* { */
/*  struct jaco_state * params = (struct jaco_state *)p; */

/*  // Kludge for CDM profile */
/*  return (r/params->r_o); */
/* } */


/* //define stellar density profile */
/* double nu_profile(double r, void * p) */
/* { */
/*   struct jaco_state * params = (struct jaco_state *)p; */
/*   double r_nu = (params->r_nu); */
/*   double nu_o = (params->nu_o); */
/*   double x = r / r_nu; */
/*   double y = pow(x, 0.25); */
/*   double z = exp( -1 * 7 * y); */
  
/*   return nu_o * z; */
/* } */


///////////////////////////////////////////////////////////////////
///////////////////////////////////////////////////////////////////
// EVERYTHING BELOW THIS LINE MUST BE FULLY GENERAL WITH RESPECT //
// TO ARBITRARY ANISOTROPY, MASS, and STELLAR LIGHT PROFILES //////
///////////////////////////////////////////////////////////////////
///////////////////////////////////////////////////////////////////


//define anisotropy param divided by r to assemble in final solution
//beta is a function  defined above
double anisotropy(double r, void *p)
{
  // Andi modification: this is the function we chose
 
  struct jaco_vars * jv = (struct jaco_vars *)p;
  struct jaco_state *js = jv->state;

  //return (js->aniso0+js->aniso1*js->raniso)/(1.+js->raniso);
  return js->aniso0+js->aniso1*r/(r+js->raniso);

  //return (params->beta(r, p)) / r;
}

//define surface brightness at each projected radius I(R)
//this integral is in terms of z = sqrt(r^2 - R^2) from fig2-3 B&T
//change integrand limits in corresponding SB_int function
double surface_brightness_def(double z, void *p) 
{

  struct jaco_vars * jv = (struct jaco_vars *)p;
  struct jaco_state *js = jv->state;

  double R = (jv->R1);

  double r = sqrt( R*R + z*z);
  double y = (js->stellardensity(r, js));  //replaced here
  //when alter function to put in terms of dz, r/sqrt(R*R - r*r) in 4-57 goes away

  return 2 * y;
}

//integrate beta function, INCLUDES e^() of answer!
//maybe put e^x outside of funct later 
double funct_beta_int(double r, struct jaco_vars *jv)
{
  double result, error;
  size_t neval;

  gsl_integration_workspace*w = gsl_integration_workspace_alloc (1000);
  gsl_function F;
  F.function = &anisotropy;
  F.params = jv;


  //gsl_integration_qags (&F, 1., r, 0, jv->state->precision/100., 1000, w, &result, &error);
  gsl_integration_qag (&F, 1., r, 0, jv->state->precision/100., 1000, 1, w, &result, &error);
  gsl_integration_workspace_free (w);

  return exp(2 * result);
}


//assemble collective integrand for 3D radial velocity dispersion
double funct_3d_nu_times_velocity_dispersion_integrand_def(
					double r, void *p)
{ 

  struct jaco_vars *jv = (struct jaco_vars *)p;
  struct jaco_state *js = jv->state;

  return ( (js->mass(r, js)) * (js->stellardensity(r, js)) * 
	   funct_beta_int(r, jv) ) / (r*r);  ;

}

//integrate collective integral for radial velocity dispersion
double funct_3d_nu_times_velocity_dispersion(double r, struct jaco_vars *jv)
{
  
  struct jaco_state *js = jv->state;

  double result, error;
  double r_o = (js->rshock);
  size_t neval;
  
  gsl_integration_workspace*w = gsl_integration_workspace_alloc (1000);
  gsl_function F;
  F.function = &funct_3d_nu_times_velocity_dispersion_integrand_def;
  F.params=jv;
  if (js->debug) printf("Entering 3D integral %E %E\n",r,r_o);
  //gsl_integration_qags (&F, r, r_o, 0, js->precision/10., 1000, w, &result, &error);
  gsl_integration_qag (&F, r, r_o, 0, js->precision/10., 1000, 1, w, &result, &error);
  if (js->debug) printf("Exiting 3D integral");
  gsl_integration_workspace_free (w);
  
  return result;
}

//integrates surface_brightness_df() 
//to find surface brightness at each projected radius I(R)
//here we have put the function in terms of z
//so have changed former limits (R, r_o), still iterating over R
double surface_brightness_int(double R, struct jaco_vars *jv)
{ 
  struct jaco_state * js = jv->state;
 
  double result, error; 
  size_t neval;

  double r_o = (js->rshock);
  if (r_o < R) {
    printf("Warning---data extends beyond truncation radius. Returning 0\n");
    return 0.;
  }
  double lim2 = sqrt(r_o * r_o - R * R); //new limit
  
  gsl_integration_workspace * w = gsl_integration_workspace_alloc (1000);
  gsl_function F;
  F.function = &surface_brightness_def; 
  F.params = jv;
  if (js->debug) printf("Entering SB integral %E %E\n",R,r_o);
  //gsl_integration_qags (&F, 0, lim2, 0, js->precision/100., 1000, w, &result, &error);
  gsl_integration_qag (&F, 0, lim2, 0, js->precision/100., 1000, 1, w, &result, &error);
  
  if (js->debug) printf("Exiting SB integral %E %E\n",R,r_o);
  gsl_integration_workspace_free (w);

  return result;
}

//defines velocity dispersion 
//already cancelled out nu from 4-60 with int factor
double func_2D_velocity_dispersion_def(double z, void *p)
{
  struct jaco_vars *jv = (struct jaco_vars *)p;
  struct jaco_state * js = jv->state;

  double R = (jv->R1); 
  double r = sqrt(R*R+z*z);

  if (js->debug) printf("Calling 3D integral R=%E z=%E\n",R,z);
  //new def in terms of dz
  double y = ( 1 - ( (anisotropy(r, jv)) * (R / r) * (R / r) ) );
  double k = funct_3d_nu_times_velocity_dispersion(r, jv) / 
    funct_beta_int(r, jv) ;

  return y * k;
}

//do entire integral in terms of z bc of singularity
double funct_2D_velocity_dispersion_int(double R, struct jaco_vars * jv)
{
  struct jaco_state * js = jv->state;
  size_t neval;
  
  double result, error;
  double r_o = (js->rshock);
  if (r_o < R) {
    printf("Warning---data extends beyond truncation radius. Returning 0\n");
    return 0.;
  }
  double lim2 = sqrt(r_o*r_o - R*R);

  gsl_integration_workspace*w = gsl_integration_workspace_alloc (1000);
  gsl_function F;
  F.function = &func_2D_velocity_dispersion_def; 
  F.params=jv;
  //gsl_integration_qags (&F, 0, lim2, 0, js->precision, 1000, w, &result, &error);
  gsl_integration_qag (&F, 0, lim2, 0, js->precision, 1000, 1, w, &result, &error);
  gsl_integration_workspace_free (w); 
  //printf("2D result %E\n", result);
  
  return 2.*result; 
}

//calculates total velocity dispersion
double velocity_dispersion(double R, struct jaco_vars * jv)
{
  
  jv->R1 = R;

  return sqrt( (2 * GDYN * funct_2D_velocity_dispersion_int(R, jv) ) 
	       / surface_brightness_int(R, jv) );
  
}

void dispersion_profile(double *Rvec, double *profile, int ndata, 
			  struct jaco_vars *jv)
{

  int i;

  for (i = 0; i < ndata; ++i) {

    profile[i] = velocity_dispersion(Rvec[i],jv);

  }

}

// Main test routine
/* int main(int argc, char *argv[]) */
/* { */

/*   double r_nu = atof(argv[1]);  */
/*   double r_o = atof(argv[2]); */
/*   double error = 1.0; */
/*   double velocity_disp = 0; */
/*   double G_Mpc = 430241; */
/*   double M_o = 1.0; */
/*   double nu_o = 1.0; */
/*   int i; */
/*   double R = 0.01; */
/*   double n = 0;//moved from loop */
/*   //just moved here from inside loop initialize functions as null */
/*   struct jaco_state params = { 0, 0, 0, G_Mpc, M_o, r_nu, nu_o, r_o, R, n }; */

/*   params.mass = M_profile;  */
/*   params.beta = beta_profile; */
/*   params.nu = nu_profile; */
  

/*  //loops v disp calculation over values of n from argv, makes output file for each one */
/*  //outputs n files, each with the same beta but a different n */
/*  //n is the power law parameter in the current mass distribution function */
/*  //plot each file on the same plot using qplot 1 dispersion5_beta* R disp disp.eps */
/*  //then convert to pdf: epstopdf */
/*  for (i = 3; i < argc; i++) */
/*    { */
     
/*      n = atof(argv[i]); */
/*      params.n = n; //added */
/*      printf("n = %f i = %d\n",n, i); */
/*      char format[]= "dispersion7_betavar_n%.lf.out"; */
/*      char filename[1000];  */
/*      sprintf(filename,format,n);  //returns number of chars in string */
/*      FILE *outp1= fopen(filename, "w"); */
      
/*      //iterate over R, 2D surface brightness */
/*      R = 0.01;  */
/*      for (; R <= r_o; R *= 1.2)      */
/*        { */
/*          params.R = R; */
         
/*          velocity_disp = velocity_dispersion(R, &params); */
 
/*          fprintf(outp1,"%10.3E %10.3E %10.3E\n", R, velocity_disp, error);  */
/*        }  */

/*      fclose(outp1);     */
/*    } */
 
/*  return 0; */
/* } */

