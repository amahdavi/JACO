/* hrothgar.c
 * 
 * Copyright (C) 2007-2014 Andisheh Mahdavi
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

#include <string.h>
#include <strings.h>
#include <math.h>
#include <stdio.h>
#include <time.h>
#include <sys/times.h>
#include <sys/types.h>
#include <sys/stat.h>
#include <getopt.h>
#include <unistd.h>
#include <hrothgar.h>
#include <hrothgar_help.h>
#include <hrothgar_proto.h>
#include <sciutils.h>
#include <gsl/gsl_math.h>
#include <gsl/gsl_eigen.h>
#include <gsl/gsl_linalg.h>
#include <regex.h>

// Report the percent done to stdout
void percentdone(long n, long tot)
{
  double perc;
  int np,i;

  perc = 100.*n/tot;
  np = printf("...%3.0f%% done  ",perc);
  for (i = 0; i < np; ++i)
    printf("\b");
  fflush(stdout);
}

// Calculate a chi square statistic, optionally storing
// the "unsquared" portion
double chisq_statistic(struct hrothgar_setup *setup)
{

  long i;
  double chisq=0.,diff,*model,*data,*error,penalty;

  model = setup->model;
  data = setup->y;
  error = setup->errors;
  penalty = setup->penalty;
  
  for (i = 0; i < setup->ndata; ++i) {
    diff = (model[i]-data[i])/error[i];
    if (diff < 0) diff -= penalty; else diff += penalty;

    if (setup->allchi) setup->allchi[i] = diff;
    chisq += diff*diff;
  }
 
  chisq += setup->logprior;
  if (setup->allchi)
    setup->allchi[setup->ndata] = chisq;

  return chisq;

}

int stepmatrix_from_covariance(struct hrothgar_setup *setup)
{

  int i,j;

  gsl_eigen_symmv_workspace *workspace = 
    gsl_eigen_symmv_alloc(setup->nfitparams);

  gsl_matrix *evectors = 
    gsl_matrix_alloc(setup->nfitparams,setup->nfitparams);
  gsl_matrix *inverse = 
    gsl_matrix_alloc(setup->nfitparams,setup->nfitparams);
  gsl_vector *evalues = gsl_vector_alloc(setup->nfitparams);
	
  if (setup->debug) {
    for (i = 0; i < setup->nfitparams; ++i) {
      for (j = 0; j < setup->nfitparams; ++j)
	printf("%6.0E ",gsl_matrix_get(setup->covar,i,j));
    printf("\n");
    }
  }

  int signum;
  gsl_eigen_symmv(setup->covar,evalues,evectors,workspace);

  gsl_permutation *permutation = gsl_permutation_alloc(setup->nfitparams);
  gsl_linalg_LU_decomp(evectors,permutation,&signum);
  gsl_linalg_LU_invert(evectors,permutation,inverse);
  gsl_permutation_free(permutation);
	
  for (i = 0; i < setup->nfitparams; ++i) {
    setup->evalues[i] = sqrt(gsl_vector_get(evalues,i));
    for (j = 0; j < setup->nfitparams; ++j)  
      setup->stepmatrix[i][j] = gsl_matrix_get(inverse,i,j);
  }
  
  if (setup->debug) {
    for (j = 0; j < setup->nfitparams; ++j) {
      printf("\nEigenvalue: %.5E\n   ",gsl_vector_get(evalues,j)); 
      for (i = 0; i < setup->nfitparams; ++i)
	printf("%5.1E ",gsl_matrix_get(evectors,j,i));
    }
    printf("\n\n");
	
    for (i = 0; i < setup->nfitparams; ++i) {
      for (j = 0; j < setup->nfitparams; ++j)
	printf("%6.0E ",setup->stepmatrix[i][j]);
    printf("\n");
    }
  }

  gsl_eigen_symmv_free(workspace);
  gsl_matrix_free(evectors);
  gsl_vector_free(evalues);

  return 0;
}

// This function returns TRUE if the pnew statistic is to be
// accepted given the pold statistic. 
int mcmc_test(double statnew, double statold, double *params,
	    struct hrothgar_setup *setup)
{
  
  int i, j, accepted;
  gsl_rng *generator;
  
  generator = (gsl_rng *)setup->generator;

  //if (setup->debug) printf("[%f %f]",statnew,statold);

  if (statnew < statold) 
    accepted = TRUE;
  else if (gsl_rng_uniform(generator) < 
	   exp((statold-statnew)/2.))
    accepted = TRUE;
  else 
    accepted = FALSE;

  // Adjust Gaussian step size based on number of accepted/rejected
  // points.
  int k;
  double del;
  for (i = 0; i < setup->nsteps; ++i) {
    j = setup->picked[i];
    if (accepted) 
      ++setup->naccepted[j];
    ++setup->ntested[j];
    //hrothgar_printf("%d ",setup->ntested[j]);
    if ((setup->ntested[j] % 100) == 0) {
      if (setup->covaronly) {
	if (setup->naccepted[j] < setup->acceptrate*50)
	  setup->gaussfrac /= 1.08;
	if (setup->naccepted[j] > setup->acceptrate*150)
	  setup->gaussfrac *= 1.08;
      } else {
	if (setup->debug)
	  hrothgar_printf("%s: %d %E %E --- stats %.2E %.2E\n",setup->paramnames[j],
			  setup->naccepted[j],setup->mcmcstep[j],setup->gaussfrac,
			  statnew,statold);
	del = (setup->pmax[j]-setup->pmin[j]);
	if (setup->naccepted[j] < setup->acceptrate*50 && setup->mcmcstep[j] > 0.00001*del)
	  setup->mcmcstep[j] /= 1.2;
	if (setup->naccepted[j] > setup->acceptrate*150 && setup->mcmcstep[j] < 0.5*del) 
	  setup->mcmcstep[j] *= 1+setup->naccepted[j]/100;
	if ((setup->naccepted[j] <= setup->acceptrate*150 && setup->naccepted[j] > 0) ||
	    setup->ntested[j] > 1000){
	    setup->mcmcgood[j] = 1;
	    ++setup->ngood;
	    if (setup->ngood < setup->nfitparams) 
	      printf("Node %d Calibrated %s with stepsize %E... %d parameters left to do.\n",setup->node,setup->paramnames[j],setup->mcmcstep[j],setup->nfitparams-setup->ngood);
	}
      }
      setup->naccepted[j] = 0.;
    }
  }
  
  if (setup->nvary <= setup->nfitparams)
    setup->nsteps=setup->nvary;

  if (setup->ngood >= setup->nfitparams) {
    for (i = 0; i < setup->nfitparams; ++i) {
      setup->chain[i][setup->nchain] = 
	params[setup->far[i]];
    }
    ++setup->nchain;
  }
  

  return accepted;
}

int *shuffle(struct hrothgar_setup *setup)
{
  int *ord,pick,save,i,j,noft;

  ord = sci_ivector(setup->nfitparams);
  
  j = 0;
  
  if (gsl_rng_uniform_int(setup->generator,2)) 
    for (i = 0; i < setup->nfitparams; ++i) 
      if (setup->varyoften[setup->far[i]])
	ord[j++] = setup->far[i];

  noft = j;

  for (i = 0; i < setup->nfitparams; ++i) 
    if (noft == 0 || !setup->varyoften[setup->far[i]])
      ord[j++] = setup->far[i];

  for (i = setup->nfitparams-1; i > noft; --i) {
    pick = noft+gsl_rng_uniform_int(setup->generator,i+1-noft);
    save = ord[i];
    ord[i] = ord[pick];
    ord[pick] = save;

  }

  //for (i = 0; i < setup->nfitparams; ++i) 
  //  printf("%d ",ord[i]);
  //printf("\n");

  return ord;
}

// Step into a new parameter space point from the previous
// point.
void mcmc_step(double *newparams, double *oldparams,
	       struct hrothgar_setup *setup)
{
  int i,j,*ord;
  double *newstep,sum;
  static int pos = 0;

/*   sum = 0.; */
/*   for (i = 0; i < setup->ntotparams; ++i) { */
/*     sum += log10(1+fabs(oldparams[i])); */
/*     newparams[i] = oldparams[i]; */
/*   } */
/*   hrothgar_printf(" %E ",sum); */
  //  if (setup->nvary) 
  //ord = shuffle(setup);
  ord = setup->far;

  int good,k;
  size_t *bad = sci_sizetvector(setup->nfitparams);
  gsl_matrix *covar = setup->covar;

  if (setup->nchain == 100*setup->nfitparams) {
    for (i = 0; i < setup->nfitparams; ++i) {
      gsl_matrix_set(covar,i,i,
		     gsl_stats_variance(setup->chain[i],1,setup->nchain));
      if (setup->debug)
	printf("%s covariance: %E\n",setup->paramnames[setup->far[i]],gsl_matrix_get(covar,i,i));
      for (j = i+1; j < setup->nfitparams; ++j) {
	gsl_matrix_set(covar,i,j,
		       gsl_stats_covariance(setup->chain[i],1,
					    setup->chain[j],1,setup->nchain));
	gsl_matrix_set(covar,j,i,gsl_matrix_get(covar,i,j));
      }
    }
    if (setup->debug)
      printf("Resetting covariance matrix.. %E\n",setup->gaussfrac);
    setup->covaronly = TRUE;
    setup->nchain = 0;
    stepmatrix_from_covariance(setup);
  }
  
  for (i = 0; i < setup->nfitparams; ++i)   
    newparams[setup->far[i]] = oldparams[setup->far[i]];

//  if (setup->covaronly) {
//    newstep = sci_dvector(setup->nfitparams);
//
//    do {
//      good = 1;
//      for (i = 0; i < setup->nfitparams; ++i)   {
//	newstep[i] = gsl_ran_gaussian(setup->generator,setup->gaussfrac*setup->evalues[i]);
//	newparams[setup->far[i]] = oldparams[setup->far[i]];
//	setup->picked[i] = setup->far[i];
//      }
//     
//      for (i = 0; i < setup->nfitparams; ++i)   {
//	k = setup->far[i];
//	
//	for (j = 0; j < setup->nfitparams; ++j) {
/* 	  if (strstr(setup->paramnames[k],"bnorm0")) */
/* 	    printf("%10d %10d %10E FOO\n",i,j,newstep[j]*setup->stepmatrix[j][i]); */
//	  newparams[k] += setup->stepmatrix[j][i]*newstep[j];
//	}
//	
//	if (newparams[k] < setup->pmin[k] ||
//	    newparams[k] > setup->pmax[k]) {
//	  good = 0;
//	  /*	  printf("Out of bounds %s %E\n",setup->paramnames[k],newparams[k]);*/
/* 	  ++bad[i]; */
/* 	  if ((bad[i] % 3) == 0) */
/* 	    setup->gaussfrac /= 1.2; */
//	}
//      }
//    } while (!good);
//    free(newstep);
//  } else {
//    for (j = 0; j < setup->nsteps; ++j) {
//      do {
//	setup->picked[j] = i = ord[pos];
//	++pos;
//	if (pos == setup->nfitparams) pos = 0;
//      } while (setup->mcmcgood[i] && setup->ngood < setup->nfitparams);
//	
//      do {
//	newparams[i] = oldparams[i]+
//	  gsl_ran_gaussian(setup->generator,setup->mcmcstep[i]);
//      } while (newparams[i] < setup->pmin[i] ||
//	       newparams[i] > setup->pmax[i]);
//    }
//
//  }
  /* Siegel Begin */
  if (setup->covaronly) {
	
    newstep = sci_dvector(setup->nfitparams);

	for (i = 0; i < setup->nfitparams; ++i)
	{
		newstep[i] = gsl_ran_gaussian(setup->generator,setup->gaussfrac*setup->evalues[i]);
		newparams[setup->far[i]] = oldparams[setup->far[i]];
		setup->picked[i] = setup->far[i];
	}

	for (i = 0; i < setup->nfitparams; ++i)
	{	
		k = setup->far[i];

		for (j = 0; j < setup->nfitparams; ++j) 
		{
			newparams[k] += setup->stepmatrix[j][i]*newstep[j];			
		}
		
		do {
			if (newparams[k] < setup->pmin[k])
			{
				newparams[k] = 2.0*setup->pmin[k] - newparams[k];
			}
			else if (newparams[k] > setup->pmax[k])
			{
				newparams[k] = 2.0*setup->pmax[k] - newparams[k];
			}
		} while (newparams[k] < setup->pmin[k] || newparams[k] > setup->pmax[k]);

	}
    free(newstep);
  } 
  else
  {
	for (j = 0; j < setup->nsteps; ++j) 
	{
    	do {
			setup->picked[j] = i = ord[pos];
			++pos;
			if (pos == setup->nfitparams) pos = 0;
		} while (setup->mcmcgood[i] && setup->ngood < setup->nfitparams);
		
		newparams[i] = oldparams[i]+gsl_ran_gaussian(setup->generator,setup->mcmcstep[i]);
		
		do {
			if (newparams[i] < setup->pmin[i])
			{
				newparams[i] = 2.0*setup->pmin[i] - newparams[i];
			}
			else if (newparams[i] > setup->pmax[i])
			{
				newparams[i] = 2.0*setup->pmax[i] - newparams[i];
			}
		} while (newparams[i] < setup->pmin[i] || newparams[i] > setup->pmax[i]);
    }
  }
  /* Siegel End */
  free(bad);
/*   for (j = 0; j < setup->nfitparams; ++j) { */
/*     k = setup->far[j]; */
/*     if (fabs(oldparams[k]-newparams[k]) > 1.e-4) */
/*       hrothgar_printf("%d %s %E -> %E\n",setup->nsteps,setup->paramnames[k],oldparams[k],newparams[k]); */
/*   } */
  //printf("%E %E\n",newparams[setup->far[0]]-oldparams[setup->far[0]],
  // newparams[setup->far[1]]-oldparams[setup->far[1]]);
  //if (setup->nvary) free(ord);
}


// This is the most used minimization routine, at the core of
// everything else. It returns the evaluated user function, after
// transformation from the unbounded space
int minimize_func_func(const gsl_vector *p, void *extra, gsl_vector *f)
{
  struct hrothgar_setup *setup;
  long    i;
  int    par;
  double sp;
  static double oldmodel[100][10000],chimin=1.E30,chisq;
  static int k=0;
  
  setup = (struct hrothgar_setup *)extra;

  // Transform from the unbounded sinusoidal space to the
  // "regular" parameter space
  for (i = 0; i < setup->nfitparams; ++i) {
    par = setup->far[i];
    sp = sin(gsl_vector_get(p,i));
    setup->wparams[0][par] = 0.5*(setup->pmin[par]*(1.-sp)+
				  setup->pmax[par]*(1.+sp));
  }

  setup->get_model(setup->x,setup->wparams[0],setup->model,setup->errors,
		   &setup->logprior,
		   setup->ndata,setup->data);

  if (!setup->dojacobian) 
    chisq = setup->statistic(setup);

  if (f) {
    for (i = 0; i < setup->ndata; ++i) {
      if (setup->dojacobian) 
	gsl_vector_set(f,i,setup->model[i]);
      else 
	gsl_vector_set(f,i,setup->allchi[i]);
    }
  }

  /* printf("-00000000-------%E %e\n",chisq,chimin); */
  /* if (chisq < chimin) { */
  /*   chimin = chisq; */
  /*   for (i = 0; i < setup->ndata; ++i) { */
  /*     oldmodel[k][i] = setup->model[i]; */
  /*     if (k > 1) */
  /* 	if (fabs(1.-oldmodel[k][i]/oldmodel[k-1][i]) > 1.E-2) */
  /* 	  printf("%ld %E %E %E\n",i,oldmodel[k-1][i],oldmodel[k][i],1.-oldmodel[k][i]/oldmodel[k-1][i]); */
  /*   } */
  /*   ++k; */
  /* } */

  return 0;
}

// Calculate nunits different merit functions serially on one
// CPU
void multiproc_onecpu(double **input, double **output, int nunits, 
		      struct hrothgar_setup *setup)
{
  int unit,j;
  double *fitparams = setup->fitparams;
  gsl_vector_view oview;
    

  gsl_vector_view iview = gsl_vector_view_array(fitparams,setup->nfitparams);

  for (unit = 0; unit < nunits; ++unit) {

    for (j = 0; j < setup->nfitparams; ++j) 
      fitparams[j] = input[unit][setup->far[j]];
    
    setup->dojacobian = (unit > 0);
    
    if (output != NULL) {
      oview = gsl_vector_view_array(output[unit],setup->ndata);
      minimize_func_func(&iview.vector,setup,&oview.vector);
    } else {
      minimize_func_func(&iview.vector,setup,NULL);
    }
  }

  setup->dojacobian = FALSE;

}

// This calculates the jacobian of the function as well as the
// function itself at point p.
int minimize_func_funcder(const gsl_vector *p, void *extra, 
				 gsl_vector *f, gsl_matrix *J)
{
  struct hrothgar_setup *setup;
  long    i,j;
  int    par,direction;
  static long nunits=0;
  double **bigmatrix,*h,**wparams,answer,*fitparams,*mappedpar;
  double **jacobian1,**jacobian2,**jacobian3,**jacobian4;
  double r2,r4,e4,dy,pmin,pmax,epsilon;

  setup = (struct hrothgar_setup *)extra;

  jacobian1 = setup->jacobian1;
  jacobian2 = setup->jacobian2;
  jacobian3 = setup->jacobian3;
  jacobian4 = setup->jacobian4;
  bigmatrix = setup->bigmatrix;
  wparams = setup->wparams;
  fitparams = setup->fitparams;
  h = sci_dvector(setup->nfitparams);
  mappedpar = sci_dvector(setup->nfitparams);
  
  nunits = 1+(setup->threepoint ? 2 : 4)*setup->nfitparams;

  // OpenMP allowed since we're not calling the 
  // user merit function yet yet
  #pragma omp parallel for private(i,j,par,pmin,pmax,direction)
  for (i = 0; i < setup->nfitparams; ++i) {
    par = setup->far[i];
    pmin = setup->pmin[par];
    pmax = setup->pmax[par];

    mappedpar[i] = wparams[0][par] = gsl_vector_get(p,i);

    // Direction to take the derivative in depends on the sign
    // of the mapped function
    direction = (mappedpar[i] > 0. ? -1: 1);

    // Untransformed tells us to take the Jacobian in an
    // untransformed space, rather than in the sine 
    // transformation space
    if (setup->untransformed) 
      mappedpar[i] = unmap_par(wparams[0][par],setup->pmin[par],
			       setup->pmax[par]);

    for (j = 1; j < nunits; ++j) 
      wparams[j][par] = wparams[0][par];
    
    h[i] = direction*setup->h[i];
    // Prevent zero stepsizes
    h[i] *= GSL_MAX(1.E-12,fabs(mappedpar[i]));

    // Calculate abscissae for five or three point rule four point rule
    if (setup->untransformed) {
      while (mappedpar[i]+h[i] < pmin || mappedpar[i]+h[i] > pmax)
	h[i] /= 2.;

      wparams[i+1][par] = 
	map_par(mappedpar[i]+h[i],pmin,pmax);
      wparams[i+1+setup->nfitparams][par] = 
	map_par(mappedpar[i]+h[i]/2.,pmin,pmax);
      if (!setup->threepoint) {
	wparams[i+1+2*setup->nfitparams][par] = 
	  map_par(mappedpar[i]+3.*h[i]/4.,pmin,pmax);
	wparams[i+1+3*setup->nfitparams][par] = 
	  map_par(mappedpar[i]+h[i]/4.,pmin,pmax);
      }
    } else {
      wparams[i+1][par] += h[i];
      wparams[i+1+setup->nfitparams][par] += h[i]/2.;
      if (!setup->threepoint) {
	wparams[i+1+2*setup->nfitparams][par] += 3.*h[i]/4.;
	wparams[i+1+3*setup->nfitparams][par] += h[i]/4.;
      }
    }
  }

  // Abscissae are calculated, now go get the actual merit
  // function values.
  setup->multifunc(wparams,bigmatrix,nunits,setup);

  // Calculate the function and the jacobians via the three or the
  // five point rules.
  for (i = 0; i < nunits; ++i) 
    if (i == 0) {
      if (f) 
	for (j = 0; j < setup->ndata; ++j) 
	  gsl_vector_set(f,j,bigmatrix[0][j]); 
    } else {
      for (j = 0; j < setup->ndata; ++j) {
	if (i < setup->nfitparams+1)
	  jacobian4[j][i-1] = bigmatrix[i][j];
	else if (i < 2*setup->nfitparams+1)
          jacobian2[j][i-setup->nfitparams-1] = bigmatrix[i][j];
	else if (i < 3*setup->nfitparams+1)
          jacobian3[j][i-2*setup->nfitparams-1] = bigmatrix[i][j];
	else 
          jacobian1[j][i-3*setup->nfitparams-1] = bigmatrix[i][j];
      }
    }

  for (i = 0; i < setup->ndata; ++i)
    for (j = 0; j < setup->nfitparams; ++j) {
      if (setup->threepoint) 
	answer = 2.*(jacobian4[i][j]-jacobian2[i][j])/h[j];
      else {

	// Calculate the truncation and roundoff errors
	// This code adapted from gsl/deriv/deriv.c (GPL)
	// (C) 2004 Brian Gough
	epsilon = fabs(setup->inputaccuracy);
	r2 = 2.*(jacobian4[i][j]-jacobian2[i][j]);
	r4 = (22./3.)*(jacobian4[i][j]-jacobian3[i][j])-
	  (62./3.)*(jacobian3[i][j]-jacobian2[i][j])+
	  (52./3.)*(jacobian2[i][j]-jacobian1[i][j]);
	e4 = (fabs(jacobian1[i][j])+fabs(jacobian2[i][j])+
	      fabs(jacobian3[i][j])+fabs(jacobian4[i][j]))*epsilon,
	  2.*20.67;
	dy = GSL_MAX(fabs(r2),fabs(r4))*fabs(mappedpar[j])*epsilon;
	answer = r4/h[j];
	setup->truncerr[j][i] = fabs((r4-r2)/h[j]);
	setup->rounderr[j][i] = fabs(e4/h[j])+dy;
      }
      setup->truncerr[j][i] /= setup->errors[i];
      setup->rounderr[j][i] /= setup->errors[i];
      par = setup->far[j];

      gsl_matrix_set(J,i,j,answer/setup->errors[i]);
    }

  free(mappedpar);
  free(h);
  return 0;
}


// Calculate the Jacobian as a function of a large number of
// step sizes; choose the step size with the smallest error
#define NSTEPS 50
int precision_gradient2(const gsl_vector *p, void *extra, 
		       gsl_vector *f, gsl_matrix *J)
{
  long i;
  int par,k,n,imax=0;
  char fname[1000];
  FILE *jacfile;
  size_t hind[NSTEPS],*ninacc;
  long j,ncount[NSTEPS];
  struct hrothgar_setup *setup;
  gsl_matrix *Jtrial[NSTEPS];
  double answer,***herror,*saveh[NSTEPS],v,pmin,pmax,*errvec,del,mmax,emax=0.;

  setup = (struct hrothgar_setup *)extra;

  if (setup->debug) printf("\n");
  setup->threepoint = FALSE;
  herror = sci_dtensor(setup->nfitparams,setup->ndata,NSTEPS);
  errvec = sci_dvector(setup->ndata);
  ninacc = sci_sizetvector(setup->nfitparams);
  del = pow(1000.,1./(NSTEPS-1));

  for (n = 0; n < NSTEPS; ++n) {

    saveh[n] = sci_dvector(setup->nfitparams);
    Jtrial[n] = gsl_matrix_alloc(setup->ndata,setup->nfitparams);

    for (i = 0; i < setup->nfitparams; ++i) {
      par = setup->far[i];
      pmin = setup->pmin[par]; pmax = setup->pmax[par];
      v = unmap_par(gsl_vector_get(p,i),pmin,pmax);

      if (n == NSTEPS-1 && (fabs(v/(pmax-pmin)) > 1.E-5) &&
	  ((v-pmin)/(pmax-pmin) < setup->stepsize || 
	   (pmax-v)/(pmax-pmin) < setup->stepsize))
	setup->h[i] = setup->stepsize*(pmax-pmin)/v;
      else
	setup->h[i] = setup->stepsize*pow(del,n-(NSTEPS-1)/2);

      saveh[n][i] = setup->h[i];
    }
    if (setup->debug) printf("%7.1E ",setup->h[0]);
    else
      if (setup->messages)
	percentdone(n,NSTEPS);
    fflush(stdout);
   
    minimize_func_funcder(p,extra,f,Jtrial[n]);
    for (i = 0; i < setup->nfitparams; ++i) 
      for (j = 0; j < setup->ndata; ++j) {
	herror[i][j][n] = setup->rounderr[i][j]+setup->truncerr[i][j];
	//herror[i][j][n] += fabs(0.5*gsl_matrix_get(Jtrial[n],j,i));
    //(GSL_DBL_EPSILON+fabs(gsl_matrix_get(Jtrial[n],j,i)));
      }
  }
  if (!setup->debug && setup->messages)
    percentdone(1,1);
  hrothgar_printf("\n");

  for (i = 0; i < setup->nfitparams; ++i) {
    k = 0;

    for (n = 0; n < NSTEPS; ++n)
      ncount[n] = 0;

    mmax = -1.E30;
    setup->h[i] = 0.;
    for (j = 0; j < setup->ndata; ++j) {
      for (n = 0; n < NSTEPS; ++n)
	hind[n] = n;
      gsl_sort_index(hind,herror[i][j],1,NSTEPS);
      ++ncount[hind[0]];
      answer = gsl_matrix_get(Jtrial[hind[0]],j,i);
      if (fabs(answer) > 0) errvec[k] = fabs(answer);
      if (answer > mmax) { mmax = answer; emax = herror[i][j][hind[0]]; imax = hind[0]; }
      k++;
      gsl_matrix_set(J,j,i,answer);
      if (setup->Jerr)
	gsl_matrix_set(setup->Jerr,j,i,herror[i][j][hind[0]]);
      if (herror[i][j][hind[0]] > fabs(0.5*answer))
	++ninacc[i];
      setup->h[i] += saveh[hind[0]][i];
    }
    setup->h[i] /= setup->ndata;
    if (setup->debug) {
      gsl_sort(errvec,1,k);
      printf("%9s ",setup->paramnames[setup->far[i]]);
      for (n = 0; n < NSTEPS; ++n) {
	printf("%5ld",ncount[n]);
	if (n == imax) printf("*"); else printf(" ");
      }
      printf("; %.1E--%.1E (%.1E) (%Zu)\n",
	     gsl_stats_quantile_from_sorted_data(errvec,1,k,0.1),
	     gsl_stats_quantile_from_sorted_data(errvec,1,k,0.9),
	     setup->h[i],ninacc[i]); 
    }
  }

  if (setup->debug && strcmp(setup->outfilename,"/dev/null")) {
    strncpy(fname,setup->outfilename,990);
    strcat(fname,".jac");
    jacfile = fopen(fname,"w+");
    for (i = 0; i < setup->ndata; ++i) {
      fprintf(jacfile,"%10ld",i);
      for (j = 0; j < setup->nfitparams; ++j) 
	fprintf(jacfile,"%10.1E",gsl_matrix_get(J,i,j));
      fprintf(jacfile,"\n");
    }
    fclose(jacfile);
  }


  setup->threepoint = FALSE;
  return 0;
}

// Calculate only the derivative, not the actual function. Just a wrapper
int minimize_func_der(const gsl_vector *p, void *extra, gsl_matrix *J)
{
  struct hrothgar_setup *setup;

  setup = (struct hrothgar_setup *)extra;
  minimize_func_funcder(p,extra,NULL,J);
    
  return 0;
}

// In order to provide bounded minimization, parameters are
// transformed via an arcsin function.
double unmap_par(double x, double min, double max)
{
  double sp;

  sp = sin(x);
  return 0.5*(min*(1.-sp)+max*(1.+sp));
}

double map_par(double x, double min, double max)
{
  x = (2.*x-max-min)/(max-min);

  x = GSL_MAX(-1,x);
  x = GSL_MIN(1,x);

  return asin(x);
}
  
void print_state(FILE *of, const char *prefix, double chi2,
		   gsl_vector *p, struct hrothgar_setup *setup)
{
  int i,par;
  long dof=0;
  
  if (setup->ndata)
    dof = setup->ndata-setup->nfitparams-1;
  if (prefix != NULL) fprintf(of,"%s",prefix);
  fprintf(of,"Stat = %f; ",chi2);
  if (chi2 > 0) {
    fprintf(of,"DOF = %ld; ",dof);
    if (setup->ndata) {
      fprintf(of,"Ratio = %f; ",chi2/dof);
      fprintf(of,"q = %f\n",gsl_sf_gamma_inc_Q(dof/2.,chi2/2.));
    }
  }
  if (prefix) fprintf(of,"%s",prefix);
  for (i = 0; i < setup->nfitparams; ++i) {
    par = setup->far[i];
      fprintf(of,"| %9.8s%13.3E ",setup->paramnames[par],
	      unmap_par(gsl_vector_get(p,i),setup->pmin[par],
			setup->pmax[par]));

    if (i % 3 == 2 || i+1 == setup->nfitparams) {
      fprintf(of,"|\n");
      if (prefix != NULL) fprintf(of,"%s",prefix);
    }
  }
  fprintf(of,"\n");
}

void init_minimize_workspace(struct hrothgar_setup *setup)
{
  int i,j,nlarge;

  nlarge = GSL_MAX(1+4*setup->nfitparams,setup->nchebtot+1);

  // Allocate minimization workspaces
  setup->fitparams = sci_dvector(setup->nfitparams+1);
  setup->jacobian1 = sci_dmatrix(setup->ndata,setup->nfitparams);
  setup->jacobian2 = sci_dmatrix(setup->ndata,setup->nfitparams);
  setup->jacobian3 = sci_dmatrix(setup->ndata,setup->nfitparams);
  setup->jacobian4 = sci_dmatrix(setup->ndata,setup->nfitparams);
  setup->truncerr = sci_dmatrix(setup->nfitparams,setup->ndata);
  setup->rounderr = sci_dmatrix(setup->nfitparams,setup->ndata);
  setup->bigmatrix = 
    sci_dmatrix(nlarge,setup->ndata);
  setup->wparams = sci_dmatrix(nlarge,setup->ntotparams+1);
  setup->allchi = sci_dvector(setup->ndata+setup->nfitparams);
  setup->sineparams = sci_dvector(setup->ntotparams);

  for (i = 0; i < setup->ntotparams; ++i) 
    setup->sineparams[i] = map_par(setup->params[i],setup->pmin[i],setup->pmax[i]);

  for (j = 0; j < setup->ntotparams; ++j) {

    if (!setup->frozen[j])
      setup->wparams[0][j] = setup->sineparams[j];
    else
      setup->wparams[0][j] = setup->params[j];

    for (i = 1; i < nlarge; ++i) 
      setup->wparams[i][j] = setup->wparams[0][j];
  }
}

void free_minimize_workspace(struct hrothgar_setup *setup)
{
  int nlarge;

  nlarge = GSL_MAX(1+4*setup->nfitparams,setup->nchebtot+1);

  free(setup->fitparams);
  sci_free_dmatrix(setup->jacobian1,setup->ndata);
  sci_free_dmatrix(setup->jacobian2,setup->ndata);
  sci_free_dmatrix(setup->jacobian3,setup->ndata);
  sci_free_dmatrix(setup->jacobian4,setup->ndata);
  sci_free_dmatrix(setup->bigmatrix,nlarge);
  sci_free_dmatrix(setup->wparams,nlarge);
  sci_free_dmatrix(setup->truncerr,setup->nfitparams);
  sci_free_dmatrix(setup->rounderr,setup->nfitparams);
  free(setup->allchi);
  free(setup->sineparams);
}

void confidence_2d(double cinfo[], double limits[], struct hrothgar_setup *setup) {

  double sxx,syy,sxy,determinant,dx,dy,dxdy;
  double x1,x2,y1,y2,sum=0.,xppeak,yppeak,xmin,xmax,ymin,ymax;
  long gridsize,i,j,k,ind;
  float *buffer;
  double *pdfvec,*x,*y,norm,xd,yd,*xsum,*ysum,newcinfo[NCINFO];
  size_t *pdfind,*xind,*yind;
  char outfile[1000];
  int ngrid,px,py;
  size_t xppos,yppos;
  char cfname[1000];

  // Create a subdirectory if necessary and cd to it

  limits[0] = xppeak = cinfo[0];
  limits[3] = yppeak = cinfo[1];
  
  x1 = cinfo[2];
  x2 = cinfo[3];
  y1 = cinfo[4];
  y2 = cinfo[5];
  
  // Find the intersection of the relevant terms from
  // the covariance matrix
  sxx = cinfo[6];
  syy = cinfo[7];
  sxy = cinfo[8];
  determinant = sxx*syy-sxy*sxy;

  px = (int)cinfo[9];
  py = (int)cinfo[10];


  xmin = setup->pmin[px]; xmax = setup->pmax[px];
  ymin = setup->pmin[py]; ymax = setup->pmax[py];

  ngrid = setup->ngrid;

  //if (setup->debug)
  // printf("%s x1:%f x2:%f y1:%f y2:%f\nxppeak: %f, yppeak: %f\nsxx:%f syy:%f sxy:%f\n",
  //   outfile,x1,x2,y1,y2,xppeak,yppeak,sxx,syy,sxy);

  gridsize = ngrid*ngrid;
  xsum = sci_dvector(ngrid);
  ysum = sci_dvector(ngrid);
  pdfvec = sci_dvector(gridsize);
  pdfind = sci_sizetvector(gridsize);
  xind = sci_sizetvector(gridsize);
  yind = sci_sizetvector(gridsize);
  
  x = sci_dvector(ngrid);
  y = sci_dvector(ngrid);

  // The normalization is technically not needed because we renormalize
  // the plot at the end (to account for step-function priors). However,
  // we'll keep it in for now.
  norm = 2.*M_PI*sqrt(determinant);
  buffer = (float *)malloc(gridsize*sizeof(float));

  // Initialize grid
  dx = (x2-x1)/ngrid;
  dy = (y2-y1)/ngrid;
  dxdy = dx*dy;

  for (i = 0; i < ngrid; ++i) {
    x[i] = x1+((double)i+0.5)*dx;
    y[i] = y1+((double)i+0.5)*dy;
    xsum[i] = ysum[i] = 0.;
  }

  xppos = gsl_interp_bsearch(x,xppeak,0,ngrid-1);
  yppos = gsl_interp_bsearch(y,yppeak,0,ngrid-1);
    
  sum = 0.;
  // Calculate probability image
  for (i = 0; i < ngrid; ++i) {
    for (j = 0; j < ngrid; ++j) {
	
      // This ordering of the grid is required by CFITSIO
      // (Rather than i*ngrid+j)
      ind = j*ngrid+i;
      pdfind[ind] = ind;
	
      xd = x[i]-xppeak;
      yd = y[j] - yppeak;
      pdfvec[ind] = exp((-syy*xd*xd-sxx*yd*yd+2.*sxy*xd*yd)/(2.*determinant));
      pdfvec[ind] /= norm;
      xind[ind] = i;
      yind[ind] = j;

      xsum[i] += pdfvec[ind];
      ysum[j] += pdfvec[ind];

      }
    }

  // Sort to get the cumulative probability
  gsl_sort_index(pdfind,pdfvec,1,gridsize);
  for (ind = gridsize-1; ind >= 0; --ind) {
    k = pdfind[ind];
    sum += pdfvec[k]*dx*dy;
    buffer[k] = (setup->cumulative ? sum : pdfvec[k]);
  }
  // Renormalize---this assumes we've gone far enough out in our plot
  // limits (true as long as setup->width is > 3* the Gaussian sigma)

  
  newcinfo[1] = newcinfo[4] = 1.E30;
  newcinfo[2] = newcinfo[5] = -1.E30;

  int foundlims;
  foundlims = 0;
  for (ind = gridsize-1; ind >= 0; --ind) {
    k = pdfind[ind];
    //buffer[k] /= sum;
    if (buffer[k] < (cinfo[11] > 0 ? 0.997: setup->conf)) {
      newcinfo[1] = GSL_MIN(newcinfo[1],x[GSL_MIN(xppos,xind[k])]);
      newcinfo[2] = GSL_MAX(newcinfo[2],x[GSL_MAX(xppos,xind[k])]);
      newcinfo[4] = GSL_MIN(newcinfo[4],y[GSL_MIN(yppos,yind[k])]);
      newcinfo[5] = GSL_MAX(newcinfo[5],y[GSL_MAX(yppos,yind[k])]);;
      foundlims = 1;
    }
  }

  if (!foundlims) {
    newcinfo[1] = xmin;
    newcinfo[2] = xmax;
    newcinfo[4] = ymin;
    newcinfo[5] = ymax;
  }

  if (setup->debug)
    printf("Round %.0f %E %E %E %E\n",cinfo[11],newcinfo[1],
	   newcinfo[2],newcinfo[4],newcinfo[5]);

  limits[1] = newcinfo[1];
  limits[2] = newcinfo[2];
  limits[4] = newcinfo[4];
  limits[5] = newcinfo[5];

  if (cinfo[11] > 0) {

/*   for (i = 0; i < ngrid; ++i) { */
/*     xsum[i] /= sum; */
/*     ysum[i] /= sum; */
/*   } */

/*   // Calculate 1D confidence contours */
/*   sum = 0.; */

/*   for (ind = gridsize-1; ind >= 0; --ind) { */
/*     k = pdfind[ind]; */
/*     sum += pdfvec[k]; */
/*     buffer[k] = (setup->cumulative ? sum : pdfvec[k]); */
/*   } */

/*   for (j = xppos; j >= 0 && sum < setup->conf/2.; --j) */
/*     sum += xsum[j]; */
/*   limits[1] = x[GSL_MAX(0,j)]; */
/*   for (j = xppos+1; j < ngrid && sum < setup->conf; ++j) */
/*     sum += xsum[j]; */
/*   limits[2] = x[GSL_MIN(ngrid-1,j)]; */
/*   sum = 0.; */
/*   for (j = yppos; j >= 0 && sum < setup->conf/2.; --j) */
/*     sum += ysum[j]; */
/*   limits[4] = y[GSL_MAX(0,j)]; */
/*   for (j = yppos+1; j < ngrid && sum < setup->conf; ++j) */
/*     sum += ysum[j]; */
/*   limits[5] = y[GSL_MIN(ngrid-1,j)]; */

  // Because of parameter limits, the initially chosen limits
  // might be too large. In this case shrink them.
/*   if (fabs(limits[1]-x1) > 3.*fabs(limits[1]-xppeak)) { */
/*     newcinfo[2] = xppeak-3.*(xppeak-limits[1]); */
/*     confidence_2d(newcinfo,limits,setup); */
/*   } else  if (fabs(limits[2]-x2) > 3.*fabs(limits[2]-xppeak)) { */
/*     newcinfo[3] = xppeak+3.*(limits[2]-xppeak); */
/*     confidence_2d(newcinfo,limits,setup); */
/*   } else if (fabs(limits[4]-y1) > 3.*fabs(limits[4]-yppeak)) { */
/*     newcinfo[4] = yppeak-3.*(yppeak-limits[4]); */
/*     confidence_2d(newcinfo,limits,setup); */
/*   } else if (fabs(limits[5]-y2) > 3.*fabs(limits[5]-yppeak)) { */
/*     newcinfo[5] = yppeak+3.*(limits[5]-yppeak); */
/*     confidence_2d(newcinfo,limits,setup); */
/*   } else  */
#ifdef HAVE_LIBCFITSIO
    if (setup->makecontours && strcmp(setup->outfilename,"/dev/null")) {
      strncpy(cfname,setup->outfilename,990);
      strcat(cfname,".conf");
      mkdir(cfname,0777);
      if (chdir(cfname))
	BYE("Could not change working directory to %s.",cfname);
      
      snprintf(outfile,1000,"%s--%s.fits",
	       setup->paramnames[px],setup->paramnames[py]);
      
      if (setup->overwrite) unlink(outfile);
      
      printf("%s %E %E %E %E\n",outfile,x1,x2,y1,y2);
      hrothgar_writeimage(outfile,buffer,x1,x2,y1,y2,ngrid);
      chdir("..");
    }
#endif

  }

  free(x); free(y);
  free(pdfind);
  free(pdfvec);
  free(buffer);
  free(xsum);
  free(ysum);
  free(xind);
  free(yind);
}

void initialize_confidence(double *cinfo,
			  int i, int j,
			  struct hrothgar_setup *setup)
{
  double xpeak,ypeak,xx,yy,*v;
  int p1,p2;

  p1 = setup->far[i];
  p2 = setup->far[j];
  
  v = gsl_vector_ptr(setup->s,0);

  cinfo[0] = xpeak = unmap_par(v[i],setup->pmin[p1],setup->pmax[p1]);
  cinfo[1] = ypeak = unmap_par(v[j],setup->pmin[p2],setup->pmax[p2]);
  
  xx = setup->width*sqrt(gsl_matrix_get(setup->covar,i,i));
  yy = setup->width*sqrt(gsl_matrix_get(setup->covar,j,j));
  
  cinfo[2] = GSL_MAX(setup->pmin[p1],xpeak-xx);
  cinfo[3] = GSL_MIN(setup->pmax[p1],xpeak+xx);
  cinfo[4] = GSL_MAX(setup->pmin[p2],ypeak-yy);
  cinfo[5] = GSL_MIN(setup->pmax[p2],ypeak+yy);
  
  // Find the intersection of the relevant terms from
  // the covariance matrix
  cinfo[6] = gsl_matrix_get(setup->covar,i,i);
  cinfo[7] = gsl_matrix_get(setup->covar,j,j);
  cinfo[8] = gsl_matrix_get(setup->covar,i,j);
  cinfo[9] = (double)p1;
  cinfo[10] = (double)p2;
  cinfo[11] = 1;

}

void confidence_report(int i, double info[],
		       struct hrothgar_setup *setup)
{
  int par;

  par = setup->far[i];

  printf("%12s %12.4E -%-12.4E +%-12.4E\n",setup->paramnames[par],info[0],
	 info[0]-info[1],info[2]-info[0]);
  fprintf(setup->outfile,
	  "#%12s %12.4E -%-12.4E +%-12.4E\n",setup->paramnames[par],info[0],
	  info[0]-info[1],info[2]-info[0]);
  setup->params[par] = info[0];
  setup->pmin[par] = info[1];
  setup->pmax[par] = info[2];
}

void confidence_looplimits(struct hrothgar_setup *setup)
{
  int i,j,k,nstep,par,ns=0;
  double **input,gsigma;

  input = sci_dmatrix(setup->nloop,setup->ntotparams+1);

  if (!setup->nvary) nstep = 1+setup->nfitparams/3;
  else nstep = setup->nvary;
  gsigma = setup->gaussfrac;

  ns = 0;
  for (i = 0; i < setup->nloop; ++i) {

    for (j = 0; j < setup->nfitparams; ++j)
      input[i][setup->far[j]] = setup->sineparams[setup->far[j]];

    if (i > setup->ncpu) {
      if (ns == nstep) ns = 1;
      else ++ns;

      //gsigma = setup->gaussfrac*cbrt((1+setup->nfitparams/3)/nstep);
      for (j = 0; j < ns; ++j) {
	par = gsl_rng_uniform_int(setup->generator,setup->nfitparams);
	k = setup->far[par];
      
	input[i][k] *= 1.+
	  gsl_ran_gaussian(setup->generator,gsigma);
      }
    }

  }

  setup->multifunc(input,NULL,setup->nloop,setup);
  sci_free_dmatrix(input,setup->nloop);
}

void adjust_limits(double **limits, int i, int j, double lims[], int sense)
{

  limits[i][0] = lims[0];
  limits[j][0] = lims[3];
  //if (sense == 0) {
      limits[i][1] = GSL_MAX(limits[i][1],lims[1]);
      limits[i][2] = GSL_MIN(limits[i][2],lims[2]);
      limits[j][1] = GSL_MAX(limits[j][1],lims[4]);
      limits[j][2] = GSL_MIN(limits[j][2],lims[5]);
      //}
  /* else { */
  /*     limits[i][1] = GSL_MIN(limits[i][1],lims[1]); */
  /*     limits[i][2] = GSL_MAX(limits[i][2],lims[2]); */
  /*     limits[j][1] = GSL_MIN(limits[j][1],lims[4]); */
  /*     limits[j][2] = GSL_MAX(limits[j][2],lims[5]); */
  /* } */
}

void confoutput_standalone(struct hrothgar_setup *setup)
{
  int i,j;
  double cinfo[NCINFO],**limits,lims[6];

  limits = sci_dmatrix(setup->nfitparams,3);

  for (i = 0; i < setup->nfitparams; ++i) {
    limits[i][1] = -1.E30;
    limits[i][2] = 1.E30;
  }

  setup->outfile = fopen(setup->outfilename,"a");
  if (!setup->outfile) BYE("Error appending to %s.",setup->outfilename);

  for (i = 0; i < setup->nfitparams; ++i) {
    for (j = i+1; j < setup->nfitparams; ++j) {
      initialize_confidence(cinfo,i,j,setup);
      confidence_2d(cinfo,lims,setup);
      adjust_limits(limits,i,j,lims,1);
    }
  }
  
  for (i = 0; i < setup->nfitparams; ++i) {
    for (j = i+1; j < setup->nfitparams; ++j)  {
      initialize_confidence(cinfo,i,j,setup);
      cinfo[2] = limits[i][1];
      cinfo[3] = limits[i][2];
      cinfo[4] = limits[j][1];
      cinfo[5] = limits[j][2];
      cinfo[11] = 0;
      confidence_2d(cinfo,lims,setup);
      adjust_limits(limits,i,j,lims,0);
    }
  
    confidence_report(i,limits[i],setup);
  }

  sci_free_dmatrix(limits,setup->nfitparams);
  fclose(setup->outfile);
}

void write_minimize_output(gsl_vector *p,struct hrothgar_setup *setup)
{
  int i,j=0;
  char *tmpstr,*fp;
  time_t t;

  if (!setup->outfilename) return;

  if (setup->welcome_string) {
    tmpstr = sci_strdup(setup->welcome_string);
    fp = index(tmpstr,'\n');
    if (fp) *fp = 0;
  } else
    tmpstr = setup->argv[0];

  if (strcmp(setup->outfilename,"/dev/stdout"))
      setup->outfile = fopen(setup->outfilename,"w+");
  else
    setup->outfile=stdout;

  t = time(NULL);
  fprintf(setup->outfile,"# Configuration file for %s\n",tmpstr);
  fprintf(setup->outfile,"# Automatically generated %s",ctime(&t));

  if (!setup->evalonly)
    fprintf(setup->outfile,"# Random number seed: %lu\n",setup->seed);

  fprintf(setup->outfile,"# Fittable parameters (must be real numbers)\n");

  for (i = 0; i < setup->ntotparams; ++i) {
    // Not yet initialized, so p is NULL
    if (!p) {
      fprintf(setup->outfile,"%-12s %12.5E %5d   ",
	      setup->paramnames[i],setup->initvals[i],setup->dofit[i]);
      fflush(setup->outfile);
    }
    else {
      fprintf(setup->outfile,"%17s ",setup->paramnames[i]);
      if (setup->frozen[i])
	fprintf(setup->outfile,"%12.5E ",setup->params[i]);
      else
	fprintf(setup->outfile,"%12.5E ",
		unmap_par(gsl_vector_get(p,j++),
			  setup->pmin[i],setup->pmax[i]));
      fprintf(setup->outfile,"%5d ",
	      setup->remember ? setup->origfrozen[i] : setup->frozen[i]);
    }

    fprintf(setup->outfile,"%12.5E %12.5E",setup->pmin[i],setup->pmax[i]);
    if (setup->pcomment) 
      if (setup->pcomment[i])
	fprintf(setup->outfile," %s%s",
		(setup->pcomment[i][0] == '#' ? "" : "# "),setup->pcomment[i]);
    fprintf(setup->outfile,"\n");
  }

  if (setup->nstringpars) {
    fprintf(setup->outfile,"\n# Non-fittable string parameters\n");
      
    for (i = 0; i < setup->nstringpars; ++i) {
      fprintf(setup->outfile,"%-12s %-17s",
	     setup->stringparname[i],setup->stringparval[i]);
      if (setup->stringparcomment[i])
	    fprintf(setup->outfile," %s%s",
		    setup->stringparcomment[i][0] == '#' ? "" : "# ",
		    setup->stringparcomment[i]);
      fprintf(setup->outfile,"\n");
    }
  }

  if (p) {
    fprintf(setup->outfile,"# Called with the folllwing command line:\n# ");
    for (i = 0; i < setup->argc; ++i) 
      fprintf(setup->outfile,"%s ",setup->argv[i]);
    fprintf(setup->outfile,"\n");
  }

  if (strcmp(setup->outfilename,"/dev/stdout"))
    fclose(setup->outfile);

}

double linemin_function(double linpar, void *minparams)
{
  int i;
  struct hrothgar_setup *s;
  gsl_vector_view p;

  s = (struct hrothgar_setup *)minparams;

  p = gsl_vector_view_array(s->sineparams,s->nfitparams);

  for (i = 0; i < s->nfitparams; ++i) 
    s->sineparams[i] = s->fitparams[i]+linpar*s->step[i];


  minimize_func_func(&p.vector,s,NULL);
  //printf("%d: %f, %E\n",s->node,linpar,s->allchi[s->ndata]);
  return s->allchi[s->ndata];
  
}

void line_minimize(struct hrothgar_setup *s,
		   gsl_min_fminimizer *minimizer,
		   double *linpar,
		   double *minchi)
{
  int i=0,istat;
  double a,b,f[3],p[3],thresh;
  gsl_function minfunc;


  // Start at a random spot and try to get at minimal spread
  // in the starting points.
  p[0] = p[1] = p[2] = 0.;
  //p[1] = 0.1-0.2*gsl_rng_uniform(s->generator);

  if (s->thorough) 
    thresh = 10.;
  else
    thresh = 3.;

  f[1] = linemin_function(p[1],s);
  do {
    p[0] -= 1.;
    p[2] += 1.;
    f[0] = linemin_function(p[0],s);
    f[2] = linemin_function(p[2],s);

    if (++i == 5) break;
  }
  while (fabs(f[1]-f[0]) < thresh || 
	 fabs(f[1]-f[2]) < thresh);

  // Take care of unconstrained parameters
  if (fabs(f[1]-f[0]) < s->eps*f[1]) {
    *linpar = p[0];
    *minchi = f[0];
    return;
  }
  if (fabs(f[1]-f[2]) < s->eps*f[1]) {
    *linpar = p[2];
    *minchi = f[2];
    return;
  }

  // Bracket the minimum
  i = 0;
  while (f[1] > f[0] || f[1] > f[2]) {
    if (f[0] < f[2]) {
      if (fabs(f[1]-f[0]) < s->eps*f[1]) {
	*linpar = p[0];
	*minchi = f[0];
	return;
      }
      p[1] = p[0];
      f[1] = f[0];
      p[0] -= 1.;
      f[0] = linemin_function(p[0],s);
    } else {
      if (fabs(f[1]-f[2]) < s->eps*f[1]) {
	*linpar = p[2];
	*minchi = f[2];
	return;
      }
      p[1] = p[2];
      f[1] = f[2];
      p[2] += 1.;
      f[2] = linemin_function(p[2],s);
    }
  }

  minfunc.function = linemin_function;
  minfunc.params = s;

  gsl_min_fminimizer_set(minimizer,&minfunc,p[1],p[0],p[2]);
  do {
    ++i;
    istat = gsl_min_fminimizer_iterate(minimizer);
    a = gsl_min_fminimizer_x_lower(minimizer);
    b = gsl_min_fminimizer_x_upper(minimizer);
    
    istat = gsl_min_test_interval(a,b,s->eps,s->eps/10.);
    
    if (istat == GSL_SUCCESS) {
      
      *linpar = gsl_min_fminimizer_x_minimum(minimizer);
      *minchi = gsl_min_fminimizer_f_minimum(minimizer);
      return;
    }
  }
  while (istat == GSL_CONTINUE && i < s->iterations);

  // Error condition
  *minchi = -1;
}

void lord_minimize(struct hrothgar_setup *setup)
{
  int  i,j,par,iter=0,status,nsuccess=0,nfloat=0,nparsfloat=0,nswap=0;
  int  saveiter;
  double *tparams,chimin=1.E30,oldstat=1.E30,chisq=0.,sigmam,sigmap,sigma0;
  char cfname[1000];
  FILE *cf;
  gsl_multifit_fdfsolver *solver;
  gsl_multifit_function_fdf minfunc;
  gsl_vector_view   initvector,w;
  gsl_vector *f,*savefit,*savefunc;
  gsl_matrix *Jcovar;
  struct tms timerbuf;
  clock_t cl1=0;

  init_minimize_workspace(setup);
  
  minfunc.f = &minimize_func_func;
  minfunc.df = &minimize_func_der;
  minfunc.fdf = &minimize_func_funcder;

  minfunc.n = setup->ndata;
  minfunc.p = setup->nfitparams;
  minfunc.params = setup;
  
  tparams = sci_dvector(setup->nfitparams);
  w = gsl_vector_view_array(setup->fitparams,setup->nfitparams);
  savefit = gsl_vector_alloc(setup->nfitparams);
  savefunc = gsl_vector_alloc(setup->ndata);

  for (i = 0; i < setup->nfitparams; ++i) {
    par = setup->far[i];
    tparams[i] = setup->wparams[0][par];
  }
  solver = gsl_multifit_fdfsolver_alloc(setup->solvertype,
					setup->ndata,setup->nfitparams);

  initvector = gsl_vector_view_array(tparams,setup->nfitparams);
  f = gsl_vector_alloc(setup->ndata);

  if (setup->evalonly || setup->covaronly) {
    gsl_vector_memcpy(solver->x,&initvector.vector);
    minimize_func_func(&initvector.vector,setup,solver->f);
    print_state(stdout,NULL,pow(gsl_blas_dnrm2(solver->f),2.),solver->x,setup);
    write_minimize_output(solver->x,setup);
    if (setup->outfilename) {
      setup->outfile = fopen(setup->outfilename,"a");
      print_state(setup->outfile,"#",pow(gsl_blas_dnrm2(solver->f),2.),
		  solver->x,setup);
      fclose(setup->outfile);
    }
    
    if (setup->evalonly)
      return;
  }

  if (setup->nloop) {
    confidence_looplimits(setup);
    return;
  }

  saveiter = setup->iterations;
  cf = NULL;
  gsl_vector_memcpy(solver->x,&initvector.vector);
  double chiold=1.E30;

  if (!setup->covaronly) {
    while (nfloat <= setup->floataround) {

      //printf("Initializing stepsizes.");
      hrothgar_printf("Evaluating starting point.\n");
      gsl_multifit_fdfsolver_set(solver,&minfunc,&initvector.vector);
      hrothgar_printf("Entering fitting loop.\n");
      chisq = pow(gsl_blas_dnrm2(solver->f),2.);
      if (setup->messages)
	print_state(stdout,NULL,chisq,solver->x,setup);
      iter = nsuccess = 0;
      oldstat = 1.e30;
      do {
	if (setup->timer) 
	  cl1 = times(&timerbuf);
	
	status = gsl_multifit_fdfsolver_iterate(solver);
	
	if (status != GSL_CONTINUE && status != GSL_SUCCESS) break;
	//status = gsl_multifit_test_delta(solver->dx,solver->x,1.E-10,setup->eps);
	//if (status == GSL_SUCCESS) ++nsuccess;

	chisq = pow(gsl_blas_dnrm2(solver->f),2.);

	if (chisq < chiold) {
	  if (chisq/chiold > 0.95 && chisq > chimin)
	    ++iter;	  
	  chiold = chisq;
	}
	else
	  ++iter;
	
	if (setup->messages) 
	  print_state(stdout,NULL,chisq,solver->x,setup);
	
	if (chisq < chimin) {
	  if (setup->messages) printf("Best fit yet. ");
	  saveiter=setup->iterations;
	  write_minimize_output(solver->x,setup);
	} else {
	  if (setup->messages) {
	    printf("Not better than %.1f.",chimin);	
	    printf("Trying %d more times.",saveiter-iter);
	  }
	}
	
	if (oldstat-chisq < setup->eps) {
	  ++nsuccess;
	  hrothgar_printf(" - Ending LM Cycle - ");
	} else {
	  if (nsuccess > 0)
	    hrothgar_printf(" - Resuming LM Cycle - ");
	  nsuccess = 0;
	}
	
	oldstat = chisq;
	
	if (setup->timer) 
	  hrothgar_printf("  [%.3f seconds per iteration]",
		 (double)(times(&timerbuf)-cl1)/setup->timer);
	hrothgar_printf("\n");
	fflush(stdout);

      }
      while (nsuccess < 3 && iter < saveiter);

      if (iter == saveiter)
	printf("Number of iterations exceeded.\n");

      gsl_vector_memcpy(&w.vector,solver->x);

      if ((status != GSL_SUCCESS && nsuccess == 0) || 
	  iter == saveiter) {
	printf("status = %d %s\n",status,gsl_strerror(status));
	printf("Could not converge after %d iterations. Try restarting.\n",iter);
      } 

      for (i = 0; i < setup->nfitparams; ++i) {
	j = setup->far[i];
	gsl_vector_set(solver->x,i,
		       map_par(unmap_par(gsl_vector_get(solver->x,i),
					 setup->pmin[j],setup->pmax[j]),
			       setup->pmin[j],setup->pmax[j]));
      }

      if (chisq < chimin) {
	nswap = 0;
	chimin = chisq;
	gsl_vector_memcpy(savefit,solver->x);
	gsl_vector_memcpy(savefunc,solver->f);
	nparsfloat = 2;
	saveiter=70;
      } else {
	nswap = 0;
	saveiter = 70;
	gsl_vector_memcpy(solver->x,savefit);
	gsl_vector_memcpy(solver->f,savefunc);
	if (nfloat < setup->floataround) {
	  if (nparsfloat == 2) 
	    nparsfloat = 1+setup->nfitparams/2;
	  else 
	    nparsfloat = 2;
	} else {
	  nfloat = setup->floataround;
	  printf("Successful convergence.");	
	}
      }
      
      ++nfloat;
      if (nparsfloat > 0 && nfloat <= setup->floataround)
	printf(" Floating parameters---%.0f%% done...",
	       100.*nfloat/setup->floataround);

      printf("\n");

      if (nfloat <= setup->floataround) {
	printf("Going for more...\n");
	for (i = 0; i < setup->nfitparams; ++i) {
	  tparams[i] = gsl_vector_get(solver->x,i);
	  if (gsl_fcmp(fabs(tparams[i]),PI/2,setup->gaussfrac) == 0)
	    tparams[i] += 
	      gsl_ran_gaussian(setup->generator,3.*setup->gaussfrac);
	}
	for (i = 0; i < nparsfloat; ++i) {
	  par = gsl_rng_uniform_int(setup->generator,setup->nfitparams);
	  j = setup->far[par];
	
	  /*Siegel Begin 
	  if ((nfloat > 0) && ((nfloat % 3) == 0)) {
		tparams[par] = PI*(gsl_rng_uniform(setup->generator) - 0.5); 
	  } else {
	  	tparams[par] +=	gsl_ran_gaussian(setup->generator,setup->gaussfrac);
	  }
	    Siegel End */	

	  tparams[par] +=	gsl_ran_gaussian(setup->generator,setup->gaussfrac);	

	  printf("%s %.2E->%.2E\n",
		 setup->paramnames[j],
		 gsl_vector_get(solver->x,par),tparams[par]);
	}
      }
	
	printf("gaussfrac = %.2E\n", setup->gaussfrac);
	
    }

  }
  
  chimin =  pow(gsl_blas_dnrm2(solver->f),2.);
  
  hrothgar_printf("***Best fit: ");
  setup->outfile = fopen(setup->outfilename,"a");
  if (!setup->outfile) BYE("Error appending to %s",setup->outfilename);
  print_state(setup->outfile,"#",chimin,solver->x,setup);
  fclose(setup->outfile);
  setup->chimin = chimin;
  print_state(stdout,NULL,chimin,solver->x,setup);
  
  if (!setup->resume) {
    if (!setup->covaronly || setup->ignorecovar) {
      hrothgar_printf("Calculating covariance matrix.");
      setup->untransformed = TRUE;
      Jcovar = gsl_matrix_calloc(setup->ndata,setup->nfitparams);
      setup->Jerr = gsl_matrix_alloc(setup->ndata,setup->nfitparams);

      chimin = pow(gsl_blas_dnrm2(solver->f),2.);
      precision_gradient2(solver->x,setup,savefunc,solver->J);
      //minimize_func_funcder(solver->x,setup,solver->f,solver->J);
      //gsl_matrix_memcpy(setup->Jerr,solver->J);
      //gsl_matrix_scale(setup->Jerr,0.5);
      
      if (setup->debug) {
	for (i = 0; i < setup->ndata; ++i) {
	  for (j = 0;  j < setup->nfitparams; ++j)
	    printf("%10.3E ",gsl_matrix_get(solver->J,i,j));
	  printf("\n");
	}
	printf("\n");
      }

    //precision_gradient(solver->x,setup,solver->f,solver->J);
      
      gsl_multifit_covar(solver->J,0.,setup->covar);
      gsl_matrix_memcpy(Jcovar,solver->J);
      gsl_matrix_add(Jcovar,setup->Jerr);
      gsl_multifit_covar(Jcovar,0.,setup->covarp);
      gsl_matrix_memcpy(Jcovar,solver->J);
      gsl_matrix_sub(Jcovar,setup->Jerr);
      gsl_multifit_covar(Jcovar,0.,setup->covarm);
      gsl_matrix_free(Jcovar);
      gsl_matrix_free(setup->Jerr);
      setup->Jerr = NULL;

      cfname[0] = 0;
      if (strcmp(setup->outfilename,"/dev/null")) 
	strncpy(cfname,setup->outfilename,990);
      else if (strcmp(setup->infilename,"/dev/null")) 
	strncpy(cfname,setup->infilename,990);

      if (strlen(cfname)) {
	strcat(cfname,".covar");
	cf = fopen(cfname,"w+");
	if (cf) {
	  gsl_matrix_fwrite(cf,setup->covar);
	  gsl_matrix_fwrite(cf,setup->covarm);
	  gsl_matrix_fwrite(cf,setup->covarp);
	} else
	  printf("hrothgar: warning: could not write %s\n",cfname);
      }
    }

    if (setup->debug) {
      for (i = 0; i < setup->nfitparams; ++i) {
	for (j = 0;  j < setup->nfitparams; ++j)
	  printf("%10.3E ",gsl_matrix_get(setup->covar,i,j));
	printf("\n");
      }
    }
    
    if (cf) fclose(cf);

    setup->s = gsl_vector_alloc(setup->nfitparams);
    gsl_vector_memcpy(setup->s,solver->x);

    for (i = 0; i < setup->nfitparams; ++i) {
      par = setup->far[i];
      sigma0 = sqrt(gsl_matrix_get(setup->covar,i,i));
      sigmam = sqrt(gsl_matrix_get(setup->covarm,i,i));
      sigmap = sqrt(gsl_matrix_get(setup->covarp,i,i));

      hrothgar_printf("%12s %12.4E %12.4E",setup->paramnames[par],
		      unmap_par(gsl_vector_get(solver->x,i),
				setup->pmin[par],setup->pmax[par]),sigma0);
      if (setup->debug) {
	hrothgar_printf(" (%-12.4E %12.4E)",sigmap,sigmam);
      }
      if (fabs(sigma0) < GSL_DBL_EPSILON) {
	BYE("\nZero covariance encountered for %s",setup->paramnames[par]);
      }
      if (fabs(sigmam/sigma0) > 2. || fabs(sigmam/sigma0) < 0.5 ||
	  fabs(sigmap/sigma0) > 2. || fabs(sigmap/sigma0) < 0.5) {
	hrothgar_printf(" (inaccurate)");
      }

      hrothgar_printf("\n");
    }

    printf("Best-fit values and errors:\n");
  }

  gsl_multifit_fdfsolver_free(solver);
  free(tparams);
  free_minimize_workspace(setup);

}

void hrothgar_init_welcomestring(struct hrothgar_setup *setup,
				char *welcome_string)
{
  static int doneinit = 0;

  if (!doneinit)
    setup->welcome_string = welcome_string;

  doneinit = 1;

}

void hrothgar_init_stringpars(struct hrothgar_setup *setup,
			      int nstringpars, char **stringparname, 
			      char **stringparval, char **stringparcomment)
{
  static int doneinit = 0;

  if (doneinit) return;

  setup->nstringpars = nstringpars;
  setup->stringparname = stringparname;
  setup->stringparval = stringparval;
  setup->stringparcomment = stringparcomment;
  setup->ninfo = 0;
  setup->infoname = NULL;
  setup->infoarray = NULL;

  doneinit = 1;

}

void hrothgar_init_pars(struct hrothgar_setup *setup,
		  int ntotparams, char **paramnames, 
		  double *pmin, double *pmax, double *initvals, int *dofit,
		  char **pcomment)
{
  static int doneinit = 0;
  int i;

  if (doneinit) return;

  if (ntotparams == 0)
    BYE("You must call hrothgar_init_pars with valid model information.\n");

  setup->ntotparams = ntotparams;
  setup->paramnames = paramnames;
  setup->oldpmin = pmin;
  setup->oldpmax = pmax;
  setup->initvals = initvals;
  setup->dofit = dofit;
  setup->nfitparams = 0;
  for (i = 0; i < setup->ntotparams; ++i)
    if (setup->dofit[i] == 0) ++setup->nfitparams;

  setup->pcomment = pcomment;

  doneinit = 1;
}


int hrothgar_init(struct hrothgar_setup *setup,int argc, char *argv[])
{
  int j,*frz,license=1,hardlimits=0,rundefaults=0;
  int globalfreeze=0,gotmin,gotmax,*pfound,*dofit,fp,t;
  long i, conflines,fpos=0,foldpos=0;
  unsigned long nread=0;
  char line[5000],*runprog,*setparvalstring[1000],**paramnames;
  char cfname[1000];
  FILE *sf;
  int c,ntotcpu,rank,nsetnames=0,setfrozen[1000],setstate[1000],ntotparams;
  char *pn,*setparname[1000],*equalsign,**parname,**parvalue,**parcomment;
  double setparval[1000],*tmppmin,*tmppmax,*pmin,*pmax,*initvals;

  static struct option long_options[] = {
    {"gridsize", 1, 0, 'G'},
    {"gaussstep", 1, 0, 'g'},
    {"2dcontours", 0, 0, '2'},
    {"seed", 1, 0, 'S'},
    {"minimizer", 1, 0, 'M'},
    {"ncpu", 1, 0, 'n'},
    {"evalonly,",0,0,'E'},
    {"timer", 0, 0, 't'},
    {"logscatter", 0, 0, 'l'},
    {"resume", 0, 0, 'r'},
    {"debug", 0, 0, 'v'},
    {"help",0,0,'h'},
    {"license",0,0,'L'},
    {"parameter",0,0,'p'},
    {"quiet",0,0,'q'},
    {"predictable", 0, 0, 'A'},
    {"floataround", 0, 0, 'f'},
    {"hardlimits", 0, 0, 'H'},
    {"accepted", 1, 0, 'a'},
    {"stepsize", 1, 0, 's'},
    {"mcmc", 1, 0, 'm'},
    {"eps", 1, 0, 'e'},
    {"iterations", 1, 0, 'i'},
    {"overwrite", 0, 0, 'X'},
    {"remember", 0, 0, 'R'},
    {"Powell",0,0,'P'},
    {"covaronly", 0, 0, 'c'},
    {"differential", 0, 0, 'F'},
    {"confidence", 1, 0, 'C'},
    {"defaults", 0, 0, 'D'},
    {"ignorecovar", 0, 0, 'I'},
    {"dumpconfig", 0, 0, 'd'},
    {"thorough", 0, 0, 'T'},
    {"width", 1, 0, 'w'},
    {"vary", 1, 0, 'V'},
    {0, 0, 0, 0}
  };
  static char short_options[] = 
    "AC:DEFG:HILM:NPRS:V:Xa:cde:f:g:hi:l:n:m:p:qrs:vw:t2";

  hrothgar_init_welcomestring(setup,"\nHrothgar Parallel Minimizer\n");
  hrothgar_init_stringpars(setup,0,NULL,NULL,NULL);
  hrothgar_init_pars(setup,0,NULL,NULL,NULL,NULL,NULL,NULL);

  ntotparams = setup->ntotparams;

  fp = 0;
  setup->argv = argv;
  setup->argc = argc;

  if (argv[1]) 
    if (strcmp(argv[1],"-d")*strcmp(argv[1],"--dumpconfig") == 0) {
      setup->outfilename = "/dev/stdout";
      setup->pmin = setup->oldpmin;
      setup->pmax = setup->oldpmax;
      write_minimize_output(NULL,setup);
      exit(0);
    }

  ntotcpu = 1;
  rank = 0;

#ifdef HAVE_MPI
  init_mpi(&ntotcpu,&rank);
#endif
  setup->node = rank;
  setup->nchain = 0;

  if (setup->welcome_string)
    hrothgar_printf(setup->welcome_string);

  setup->generator = gsl_rng_alloc(DEFAULT_MCMC_GENERATOR);
  setup->solvertype = gsl_multifit_fdfsolver_lmder;
  setup->predictable = setup->resume = setup->seed = setup->debug = 0;
  setup->timer = 0;
  setup->messages = setup->threepoint 
    = setup->cumulative = TRUE;
  setup->ignorecovar = setup->untransformed = 
    setup->covaronly = setup->floataround = 
    setup->evalonly = setup->overwrite = setup->dojacobian = FALSE;
  setup->params = sci_dvector(ntotparams);
  setup->bestfit = sci_dvector(ntotparams);
  setup->iterations = 500;
  setup->frozen = (int *)malloc(ntotparams*sizeof(int));
  setup->origfrozen = (int *)malloc(ntotparams*sizeof(int));
  setup->seed = setup->doui = 0;
  setup->nsim = 0;
  setup->nvary = 1;
  setup->inputaccuracy = GSL_SQRT_DBL_EPSILON;
  setup->acceptrate=0.25;
  setup->stepsize = 0.001;
  setup->nloop = setup->ncpu = 0;
  setup->nchebtot = 40;
  setup->eps = 1.e-3;
  setup->width = 4.;
  setup->ntrain = 0;
  setup->ngrid = 512;
  setup->conf = 0.683;
  setup->gaussfrac = 0.3;
  setup->Jerr = NULL;
  setup->makecontours = setup->remember = setup->thorough = FALSE;

  // Look for options
  while (TRUE) {

    c = getopt_long(argc, argv, short_options, long_options, NULL);

    if (c == -1) break;


    switch (c) {

    case '2': setup->makecontours = TRUE; break;
    case 'G': setup->ngrid = atoi(optarg); break;
    case 'g': setup->gaussfrac = atof(optarg); break;
    case 'w': setup->width = atof(optarg); break;
    case 'f': setup->floataround = atoi(optarg); break;
    case 'l': setup->nloop = atol(optarg); break;
    case 't': setup->timer = TRUE; break;
    case 'r': setup->resume = TRUE; break;
    case 'v': setup->debug = TRUE; break;
    case 'A': setup->predictable = TRUE; break;
    case 'n': setup->ncpu = atoi(optarg); break;
    case 'm': setup->nsim = atol(optarg); break;
    case 'e': setup->eps = atof(optarg); break;
    case 'C': setup->conf = atof(optarg); break;
    case 'c': setup->covaronly = TRUE; break;
    case 'V': setup->nvary = atoi(optarg); break;
    case 'i': setup->iterations = atoi(optarg); break;
    case 'I': setup->ignorecovar = TRUE; break;
    case 'F': setup->cumulative = FALSE; break;
    case 'D': rundefaults = TRUE; break;
    case 'E': setup->evalonly = TRUE; break;
    case 'H': hardlimits = TRUE; break;
    case 'R': setup->remember = TRUE; break;
    case 'T': setup->thorough = TRUE; break;
    case 'X': setup->overwrite = TRUE; break;
    case 'L': license=2; break;

    case 'h': 
      printf("%s",helpstr); 
      exit(0);
      break;

    case 'd': BYE("--dumpconfig must be specified by itself.");

    case 'a': 
      setup->acceptrate = atof(optarg); 
      if (setup->acceptrate < 0. || setup->acceptrate > 0.5)
	BYE("Invalid input accuracy.");
      break;

    case 's': 
      setup->stepsize = GSL_MAX(atof(optarg),GSL_SQRT_DBL_EPSILON); 
      break;

    case 'S': 
      hrothgar_printf("Seed specified---this implies --predictable.\n");
      setup->predictable = TRUE;
      setup->seed = atol(optarg);
      break;

    case 'M':
      if (strcmp(optarg,"lmsder") == 0)
	setup->solvertype = gsl_multifit_fdfsolver_lmsder;
      else if (strcmp(optarg,"lmder") == 0)
	  setup->solvertype = gsl_multifit_fdfsolver_lmder;
      else
	BYE("Unknown minimizing solver %s. Known types are lmsder and lmder.\n",
	    optarg);
      break;

    case 'p':
      setparname[nsetnames] = sci_strdup(optarg);
      equalsign = index(setparname[nsetnames],'=');
      setstate[nsetnames] = setfrozen[nsetnames] = 0;
      if (!equalsign || equalsign != rindex(setparname[nsetnames],'=')) {
	equalsign = index(setparname[nsetnames],'~');
	setstate[nsetnames] = 1;
	if (!equalsign || equalsign != rindex(setparname[nsetnames],'~')) {
	  equalsign = index(setparname[nsetnames],'@');
	  setfrozen[nsetnames] = 1;
	  if (!equalsign || equalsign != rindex(setparname[nsetnames],'@')) {
	    equalsign = index(setparname[nsetnames],'%');
	    globalfreeze=TRUE;
	    setfrozen[nsetnames] = 0;
	    if (!equalsign || equalsign != rindex(setparname[nsetnames],'%')) 
	      BYE("Improperly formatted parameter specification: %s",optarg);
	  }
	}
      } else
	if (equalsign[1] == 0) BYE("Must assign value when using =.");
      if (equalsign[1] == 0) setstate[nsetnames] = -1;
      else {
	setparval[nsetnames] = atof(&equalsign[1]);
	setparvalstring[nsetnames] = &equalsign[1];
      }
      equalsign[0] = 0;
      ++nsetnames;
      if (nsetnames == 1000) BYE("Too many parameters");
      break;

    case 'q': 
      if (license)
	license = FALSE; 
      else
	setup->messages = FALSE;
      break;

    default: BYE("Exiting..."); break;
    }
  }

  if ((license && rank == 0)) {
    printf(setup->welcome_string);
    printf("This program uses the Hrothgar parallel LM/MCMC minimizer.\n");
    printf("Hrothgar v%s Copyright (C) 2014 Andisheh Mahdavi.\n",VERSION);
    if (license == 1) {
      printf("This program is free software; you can redistribute\n"	\
	     "it and/or modify it under certain conditions.\n"
	     "Type %s -L for details.\n",argv[0]);
    } else {
      printf(long_gpl);
      return -1;
    }
  }

  setup->outfilename = "/dev/null";
  if (optind >= argc) {
    if (rank == 0 && setup->messages) 
      printf("\nNo input file specified. ");
    if (rundefaults) {
      if (setup->messages) printf("Using defaults.\n\n"); 
      setup->infilename = "/dev/null";
    } else {
      printf("\nType %s --help for help,\nor %s --defaults",argv[0],argv[0]);
      printf(" to run with default settings.\n");
      return -1; 
    }
  } else
    setup->infilename = argv[optind];

  if (optind+1 >= argc || !argv[optind+1]) {
    if (!setup->evalonly) 
      if (rank == 0 && setup->messages) 
	printf("No output file will be written.\n"); 
  } else
    setup->outfilename = argv[optind+1];

  if (setup->evalonly) {
    if (setup->ncpu > 1)
      hrothgar_printf("Only evaluating merit function, so setting ncpu=1\n");
    setup->ncpu = 1;
  }

  runprog = rindex(argv[0],'/');
  if (runprog) ++runprog; else runprog = argv[0];
  if (setup->ignorecovar && !setup->covaronly)
    BYE("--ignorecovar (-I) is meaningless without --covaronly (-c)");

  for (i = 0; i < ntotparams; ++i) {
    setup->origfrozen[i] = setup->frozen[i] = TRUE;
    setup->params[i] = 0.;
  }

  pmin = setup->oldpmin;
  pmax = setup->oldpmax;
  paramnames = setup->paramnames;
  initvals = setup->initvals;
  dofit = setup->dofit;

  setup->pmin = sci_dvector(ntotparams);
  setup->pmax = sci_dvector(ntotparams);
  setup->varyoften = sci_ivector(ntotparams);
  pfound = sci_ivector(ntotparams);

  if (hardlimits && (!pmin || !pmax))
    BYE("Hard limits requested but not received from the calling program.");

  // Open the fit parameter file
  sf = fopen(setup->infilename,"r");
  if (sf == NULL)
    BYE("Error: cannot open fit parameter file %s for reading",
	setup->infilename);
  fclose(sf);

  conflines = read_ascii(setup->infilename,"%s %s %d %f %f",
			 &parcomment,&parname,&parvalue,
			 &frz,&tmppmin,&tmppmax);

  for (i = 0; i < conflines; ++i) {

    if (parname[i]) {
      pn = parname[i];

      if (frz[i] > -1) {

	if ((isnan(tmppmin[i]) && !pmin) || (isnan(tmppmax[i]) && !pmax))
	  BYE("Neither %s nor the calling code specifies parameter limits.",
	      argv[optind]);

	gotmax = !isnan(tmppmax[i]);
	gotmin = !isnan(tmppmin[i]);      

	for (j = 0; j < ntotparams && 
	       (fp=strcmp(pn,paramnames[j])) != 0; ++j);

	if (fp == 0) {
	
	  pfound[j] = 1;
	  //printf("Found %s (%d %E %d %e): ",paramnames[j],gotmin,tmppmax[i],gotmax,tmppmin[i]);
	  if (pmin && !hardlimits && gotmin)
	    if (tmppmin[i] < pmin[j]) {
	      hrothgar_printf("File limit on %s (%E) < hard limit %E---resetting\n",pn,tmppmin[i],pmin[j]);
	      tmppmin[i] = pmin[j];
	    }
	  
	  if (pmax && !hardlimits && gotmax)
	    if (tmppmax[i] > pmax[j]) {
	      hrothgar_printf("File limit on %s (%E) > hard limit %E---resetting\n",
		  pn,tmppmax[i],pmax[j]);
	      tmppmax[i] = pmax[j];
	    }
	      
	  setup->pmin[j] = (!gotmin || hardlimits ? pmin[j] : tmppmin[i]);
	  setup->pmax[j] = (!gotmax || hardlimits ? pmax[j] : tmppmax[i]);
	  if (frz[i] > 1) { 
	    frz[i] = 0;
	    setup->varyoften[j] = 1;
	  }
	  //printf("%E %E\n",setup->pmin[j],setup->pmax[j]);

	  setup->bestfit[j] = setup->params[j] = atof(parvalue[i]);
	  //printf("%s: %f\n",parname,parvalue);
	  setup->origfrozen[j] = frz[i];
	  setup->frozen[j] = globalfreeze || frz[i];
	  setup->nfitparams += setup->dofit[j]-setup->frozen[j];
	} else 
	  hrothgar_printf("Warning: skipping unknown parameter %s.\n",parname[i]);
      } else {
      
	for (j = 0; j < setup->nstringpars && 
	     (fp=strcmp(pn,setup->stringparname[j])) != 0; ++j);

	if (fp == 0) {
	  if (parvalue[i])
	    setup->stringparval[j] = sci_strdup(parvalue[i]);
	  else 
	    hrothgar_printf("Warning: %s specified but no value given\n",pn);
			    
	  if (parcomment[i]) 
	    setup->stringparcomment[j] = sci_strdup(parcomment[i]);
	} else
	  hrothgar_printf("Warning: skipping unknown parameter %s.\n",parname[i]);
      }
    }
  }


  for (i = 0; i < ntotparams; ++i)
    if (!pfound[i]) {
      if (setup->messages)
	hrothgar_printf("* %s not in input file; default = %E (%s)\n",
			paramnames[i],initvals[i],dofit[i] ? "frozen" : "free");
      setup->params[i] = initvals[i];
      setup->frozen[i] = dofit[i];
      setup->pmin[i] = pmin[i];
      setup->pmax[i] = pmax[i];
    }

  if (globalfreeze && rank == 0)
    hrothgar_printf("Freezing all parameters.\n");

  regex_t preg;

  for (i = 0; i < nsetnames; ++i) {
    if (regcomp(&preg,setparname[i],REG_EXTENDED | REG_NOSUB)) {
      BYE("%s is not a valid regular expression ",setparname[i]);
    }
    fp = 0;
    for (j = 0; j < ntotparams; ++j) {
      if (regexec(&preg,paramnames[j],0,NULL,0) == 0) {
	++fp;
	if (rank == 0)
	  printf("Setting %s (from %E, %s) to ",
		 paramnames[j],setup->params[j],
		 setup->frozen[j] ? "frozen" : "free");
	if (setstate[i] > -1)
	  setup->params[j] = setparval[i];
	
	if (setstate[i] != 0) {
	  if (setup->frozen[j] != setfrozen[i])
	    setup->nfitparams += 1-2*setfrozen[i];
	  setup->frozen[j] = setfrozen[i];
	}
	if (rank == 0)
	  printf("%E, %s\n",setup->params[j],setup->frozen[j] ? "frozen" : "free");
      }
    }
    for (j = 0; j < setup->nstringpars; ++j) {
      if (regexec(&preg,setup->stringparname[j],0,NULL,0) == 0) {
	++fp;
	if (rank == 0)
	  printf("Setting string parameter %s (from %s) to %s\n",
		 setup->stringparname[j],setup->stringparval[j],
		 setparvalstring[i]);
	if (setparvalstring[i] != NULL)
	  setup->stringparval[j] = sci_strdup(setparvalstring[i]);
	else
	  BYE("Must specify value for string parameter.");
	
      }
    }
    if (fp == 0) {
      BYE("Regular expression %s didn't match any parameters.",setparname[i]);
    }
    if (fp > 1) 
      printf("Warning: %d parameters matched by %s\n",fp,setparname[i]);
  }

  if (setup->nfitparams < 2)
    BYE("Sorry, you must have at least two free parameters.");

  setup->h = sci_dvector(setup->nfitparams);
  for (i = 0; i < setup->nfitparams; ++i)
    setup->h[i] = setup->stepsize;

  for (i = 0; i < ntotparams; ++i) {
    t = gsl_fcmp(setup->params[i],setup->pmin[i],1.E-8);
    if (t < 0) {
      hrothgar_printf("%s smaller than minimum bound %E--resetting.\n",
		      setup->paramnames[i],setup->pmin[i]);
      setup->params[i] = setup->pmin[i];
    }
    

    t = gsl_fcmp(setup->params[i],setup->pmax[i],1.E-8);
    if (t > 0) {
      hrothgar_printf("%s larger than maximum bound %E--resetting.\n",
		      setup->paramnames[i],setup->pmin[i]);
      setup->params[i] = setup->pmax[i];
    }
  }

  // far[i] contains the ID of the ith parameter to be fit
  setup->far = (int *)malloc(setup->nfitparams*sizeof(int));
  j = 0;
  for (i = 0; i < ntotparams; ++i) {
    if (!setup->frozen[i]) 
      setup->far[j++] = i;
  }

  if (setup->evalonly) return rank;

  if (setup->nsim) {

    // Make NSIM number of evaluations per free parameter.
    setup->nsim *= setup->nfitparams;
    setup->nsteps = setup->nvary;
    
    setup->naccepted = sci_ivector(setup->ntotparams);
    setup->ntotaccepted = 0;
    setup->ntested = sci_ivector(setup->ntotparams);
    setup->picked = sci_ivector(setup->ntotparams);
  }

  // We have specified that we are resuming a previous session. 
  // For LM mode (not treated here), this
  // means that you intend to resume the fit later (and so will not
  // need a covariance matrix right away). Nothi

  // For MCMC mode, this means to resume the chain from a previous
  // point, and so the chain is read out here

  char suffix[10];

  if (setup->nsim) {
    setup->stepmatrix = sci_dmatrix(setup->nfitparams,setup->nfitparams);
    setup->evalues = sci_dvector(setup->nfitparams);

    if (setup->ncpu == 1)
      strcpy(suffix,".mcmc");
    else  {
      if (setup->ncpu == 2) {
	BYE("MPI MCMC needs > 2 processes.");
      } else  
	 sprintf(suffix,".%02d.mcmc",rank); 
    }
    setup->mcmcfilename = sci_strcat(setup->outfilename,suffix); 

    if (rank == 0) {
      char rmcmd[1000];
      sprintf(rmcmd,"rm -f %s.??.mcmc",setup->outfilename);
      system(rmcmd);
    }
      
  }

  if (setup->nsim > 0 && setup->resume) {

    BYE("Resuming MCMC not working right now.");

  } else {

    // Not resuming, so output file should not exist
    sf = fopen(setup->outfilename,"r");
    if (strstr(setup->outfilename,"/dev") == NULL)
      if (sf != NULL && !setup->overwrite)
	BYE("Output file %s exists & no --resume or --overwrite flag specified.",
	      setup->outfilename);
    if (sf) fclose(sf);
  }

  if (setup->ncpu > ntotcpu)
    BYE("%d CPUs are not available; maximum number is %d\n",setup->ncpu,
	ntotcpu)
  else if (setup->ncpu == 0)
    setup->ncpu = ntotcpu;

  // No random number seed has been read in or set, so get one
  if (rank == 0 && setup->seed == 0) {
      hrothgar_printf("Seeding the random number generator using time().\n");
    //sf = fopen("/dev/random","r");
    //if (sf == NULL) BYE("Cannot open %s","/dev/random");
    //fread(&setup->seed,sizeof(unsigned long),1,sf);
    setup->seed = (unsigned long)time(NULL);
    //fclose(sf);
  }

  gsl_rng_set(setup->generator,setup->seed);
  hrothgar_printf("Random number seed is %lu\n",setup->seed);

  setup->covar = gsl_matrix_alloc(setup->nfitparams,setup->nfitparams);
  setup->covarm = gsl_matrix_alloc(setup->nfitparams,setup->nfitparams);
  setup->covarp = gsl_matrix_alloc(setup->nfitparams,setup->nfitparams);


  setup->mcmcstep = sci_dvector(setup->ntotparams);
  setup->mcmcgood = sci_ivector(setup->ntotparams);
  setup->ngood = 0;
  setup->chain = sci_dmatrix(setup->nfitparams,102*setup->nfitparams);

  FILE *cf;
  if (setup->covaronly) {
    strncpy(cfname,setup->infilename,990);
    if (!setup->ignorecovar) {
      strcat(cfname,".covar");
      printf("Reading covariance matrix from %s\n",cfname);
      cf = fopen(cfname,"r");
      if (cf) {
	gsl_matrix_fread(cf,setup->covar);
	gsl_matrix_fread(cf,setup->covarm);
	gsl_matrix_fread(cf,setup->covarp);
	printf("Covariance matrix for MCMC:\n");
	for (i = 0; i < setup->nfitparams; ++i)  {
	  for (j = 0; j < setup->nfitparams; ++j)
	    printf("%.5E ",gsl_matrix_get(setup->covar,i,j));
	  printf("\n");
	}
      }
      
      // Use covariance matrix for stepping through MCMC space

	stepmatrix_from_covariance(setup);
    } else if (setup->nsim) {
      BYE("Cannot ignore covariance matrix in MCMC mode");
    }
  } else if (setup->nsim) 
    for (i = 0; i < setup->ntotparams; ++i)
      setup->mcmcstep[i] = 
	setup->gaussfrac*(fabs(setup->params[i]) > 1.e-4 ?
			  fabs(setup->params[i]) :
			  1.e-4);


  return rank;
}


void report_unmapped(double *unmapped_pars, double chisq, 
		     struct hrothgar_setup *setup)
{
  int i,k;
  double *mapparams = sci_dvector(setup->nfitparams);

  for (k = 0; k < setup->nfitparams; ++k) {
    i = setup->far[k];
    mapparams[k] = map_par(unmapped_pars[k],setup->pmin[i],setup->pmax[i]);
  }

  gsl_vector_view iview = 
    gsl_vector_view_array(mapparams,setup->nfitparams);
	  
  write_minimize_output(&iview.vector,setup);
  print_state(stdout,NULL,chisq,&iview.vector,setup);
  setup->outfile = fopen(setup->outfilename,"a");
  print_state(setup->outfile,"#",chisq,&iview.vector,setup);
  fclose(setup->outfile);
  free(mapparams);
}
      


void hrothgar_mcmc(int initial, long nsim, double *mcmcdata, 
		   struct hrothgar_setup *setup)
{
  static double *fitparams,*newparams,
    chisq,oldchisq=1.E30,minchisq=1.E30;
  int keep,i,j,k,nprint=0;
  static int firsttime=1;
  static long nkeep = 0,waypoint=0,onepercent=0,ntotdone=0;
  long ndone=0;
  gsl_rng *generator;
  static clock_t cl1=0;
  struct tms timerbuf;
  double *mapparams;
  

  generator = (gsl_rng *)setup->generator;

  // Main MCMC loop.
  int mcmccols = setup->nfitparams+1+setup->ninfo;
  static double *holdinfo = NULL;
 
  if (holdinfo == NULL) holdinfo = sci_dvector(mcmccols);

  double *mean=NULL, *m2=NULL;
  double *bestfit=NULL;
  if (mcmcdata != NULL) {
    mean = mcmcdata;
    m2 = &mcmcdata[setup->nfitparams];
    bestfit = &mcmcdata[2*setup->nfitparams];
  }

  if (initial) {

    if (setup->timer) 
      cl1 = times(&timerbuf);

    waypoint = onepercent = CHUNKSIZE_MCMC/10;
    if (setup->ncpu > 1)
      waypoint = onepercent*setup->node/(setup->ncpu-1);
    // Main MCMC loop.
    setup->mcmcfile = gzopen(setup->mcmcfilename,"w");

    newparams = sci_dvector(setup->ntotparams);
    fitparams = sci_dvector(setup->nfitparams);
    mapparams = sci_dvector(setup->nfitparams);

    for (i = 0; i < setup->ntotparams; ++i)
      newparams[i] = setup->params[i];

    for (k = 0; k < setup->nfitparams; ++k) 
      fitparams[k] = setup->params[setup->far[k]];

    gzwrite(setup->mcmcfile,"HROTHGAR",9*sizeof(char));
    gzwrite(setup->mcmcfile,&mcmccols,sizeof(int));
    
    for (k = 0; k < setup->nfitparams; ++k) {
      i = setup->far[k];
      j = strlen(setup->paramnames[i])+1;
      gzwrite(setup->mcmcfile,&j,sizeof(int));
      gzwrite(setup->mcmcfile,setup->paramnames[i],sizeof(char)*j);
    }

    for (k = 0; k < setup->ninfo; ++k) {
      j = strlen(setup->infoname[k])+1;
      gzwrite(setup->mcmcfile,&j,sizeof(int));
      gzwrite(setup->mcmcfile,setup->infoname[k],sizeof(char)*j);
    }
  } else
    setup->mcmcfile = gzopen(setup->mcmcfilename,"a");

  while (ndone < nsim) {
/*     setup->get_model(setup->x,setup->params,setup->model,setup->ndata,setup->data); */

/*     foochisq=setup->statistic(setup); */
/*     if (setup->node == 1) { */
/*       char fname[1000]; */
/*       FILE *foofile; */
/*       sprintf(fname,"chisq.%8E",foochisq); */
/*       foofile = fopen(fname,"r"); */
/*       if (foofile == NULL) { */
/* 	foofile = fopen(fname,"w"); */
/* 	long l; */
/* 	for (l = 0; l < setup->ntotparams; ++l) */
/* 	fprintf(foofile,"%20ld %20E\n",l,setup->params[l]); */
/* 	for (l = 0; l < setup->ndata; ++l) */
/* 	  fprintf(foofile,"%20E %20E %20E\n",setup->y[l],setup->errors[l],setup->model[l]); */
/* 	fclose(foofile); */
/*       } */
/*     } */

    if (setup->statonly)
      chisq = setup->get_stat(newparams,setup->ntotparams,setup->data);
    else {
      setup->get_model(setup->x,newparams,setup->model,setup->errors,
		       &setup->logprior,setup->ndata,setup->data);
      chisq = setup->statistic(setup);
    }


    //hrothgar_printf("-------------------------------------%E %E %E\n",chisq,foochisq,oldchisq);

    keep = mcmc_test(chisq,oldchisq,newparams,setup);

    if (keep) {
      for (i = 0; i < setup->ntotparams; ++i)
	setup->params[i] = newparams[i];

      for (k = 0; k < setup->nfitparams; ++k) 
	fitparams[k] = setup->params[setup->far[k]];


      if (chisq < minchisq) {
	  minchisq = chisq;
	  for (k = 0; k < setup->nfitparams; ++k) {
	    setup->bestfit[setup->far[k]] = fitparams[k];
	    if (bestfit != NULL) 
	      bestfit[k] = fitparams[k];
	  }


	  if (setup->ncpu == 1)
	    report_unmapped(fitparams,chisq,setup);
      }
	
      oldchisq = chisq;
      ++setup->ntotaccepted;

      for (i = 0; i < setup->ninfo; ++i)
	holdinfo[i] = setup->infoarray[i];
      
      if (ntotdone > 0)
	++nkeep;

      if (setup->timer) {
	if (ntotdone > waypoint) {
	  waypoint += onepercent;
	  nprint = printf("Node %2d: %6ld",setup->node,ntotdone);
	  nprint += printf(", %4.0f/sec, %4.0f%% accepted (avg %4.0f%%)" \
                           " chi^2 = %.2E ",
			   setup->timer*(double)nkeep/(times(&timerbuf)-cl1),
			   100.*setup->ntotaccepted/onepercent,
			   100.*nkeep/ntotdone,chisq);
	  if (setup->covaronly)
	    nprint += printf("%.2f",setup->gaussfrac);
	  else
	    nprint += printf("Blind");
	  
	  setup->ntotaccepted = 0;
	  for (j = 0; j < nprint; ++j) printf("\b");
	  fflush(stdout);
	}
      }
    } else {
      //if (firsttime) oldchisq = chisq;
      //else
	for (i = 0; i < setup->ninfo; ++i)
	  setup->infoarray[i] = holdinfo[i];
      //firsttime = 0;
    }

    if (setup->ngood >= setup->nfitparams) {
      gzwrite(setup->mcmcfile,&oldchisq,sizeof(double));
      gzwrite(setup->mcmcfile,fitparams,sizeof(double)*setup->nfitparams);
      if (setup->ninfo)
	gzwrite(setup->mcmcfile,setup->infoarray,sizeof(double)*setup->ninfo);
      
      ++ndone;
      ++ntotdone;
    }

    double delta;
    if (mean != NULL && ndone > 0) {
      for (k = 0; k < setup->nfitparams; ++k) {
	delta = fitparams[k]-mean[k];
	mean[k] += delta/ndone;
	if (m2 != NULL) 
	  m2[k] += delta*(fitparams[k]-mean[k]);
      }
    }
    
    mcmc_step(newparams,setup->params,setup);
  }
  
  if (mcmcdata != NULL)
    mcmcdata[3*setup->nfitparams] = minchisq;

  gzclose(setup->mcmcfile);

  if (setup->ncpu == 1) {
    setup->outfile = fopen(setup->outfilename,"a");
    fprintf(setup->outfile,"# Gaussfrac: %E\n",setup->gaussfrac);
    fclose(setup->outfile);
  }
}

void hrothgar_singlecpu(struct hrothgar_setup *setup)
{
  hrothgar_printf("Initialization complete.\n");
  setup->multifunc = multiproc_onecpu;
  if (setup->nsim)
    hrothgar_mcmc(1,setup->nsim,NULL,setup); 
  else {
    lord_minimize(setup);
    if (!setup->nloop && !setup->resume)
      confoutput_standalone(setup);
  }
}

void hrothgar(unsigned long ndata, double *x, 
	      double *bigdata, double *errors,
	      int (*get_model)(double *x, double *params, double *model, 
			       double *errors, double *logprior,
			       unsigned long ndata, void *dataptr),
	      void *dataptr,
	      struct hrothgar_setup *setup)
{
  int domcmc;
  char cmd[1000];
  
  setup->statonly = FALSE;
  setup->data = dataptr;
  setup->x = x;
  setup->y = bigdata;
  setup->errors = errors;
  setup->ndata = ndata;
  setup->statistic = chisq_statistic;
  setup->get_model = get_model;
  setup->allchi = NULL;
  setup->model = sci_dvector(ndata);

  domcmc = setup->nsim;

  if (setup->timer)
    setup->timer = sysconf(_SC_CLK_TCK);

  if (setup->evalonly) {
    if (setup->node == 0) lord_minimize(setup);
    return;
  }

  printf("Wait....\n");
  if (setup->ncpu == 1) 
    hrothgar_singlecpu(setup);
#ifdef HAVE_MPI
  else
    hrothgar_mpi(setup);
#endif

  // Copy parameter file to contours directory (if minimizing)
  if (!domcmc && setup->node == 0) {
    snprintf(cmd,999,"[ -d %s.conf ] && cp -f %s %s.conf",
	     setup->outfilename,
	     setup->outfilename,
	     setup->outfilename);
    system(cmd);
  }

}

void hrothgar_statonly(double (*get_stat)(double *params, int np, 
					  void *dataptr),
		       void *dataptr,
		       struct hrothgar_setup *setup)
{
  int domcmc;
  char cmd[1000];
  
  setup->statonly = TRUE;
  setup->get_stat = get_stat;
  setup->data = dataptr;

/*   if (setup->evalonly) { */
/*     if (setup->node == 0) lord_minimize(setup); */
/*     return; */
/*   } */

  domcmc = setup->nsim;
  if (domcmc <= 0)
    BYE("Stats-only mode requires MCMC.");

  if (setup->timer)
    setup->timer = sysconf(_SC_CLK_TCK);

  printf("Wait....\n");
  if (setup->ncpu == 1) 
    hrothgar_singlecpu(setup);
#ifdef HAVE_MPI
  else
    hrothgar_mpi(setup);
#endif

  // Copy parameter file to contours directory (if minimizing)
  if (!domcmc && setup->node == 0) {
    snprintf(cmd,999,"[ -d %s.conf ] && cp -f %s %s.conf",
	     setup->outfilename,
	     setup->outfilename,
	     setup->outfilename);
    system(cmd);
  }

}
