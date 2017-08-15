/* hrothgar_mpi.c
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

#include <string.h>
#include <math.h>
#include <stdio.h>
#include <unistd.h>
#include <sys/times.h>
#include <mpi.h>
#include <sciutils.h>
#include <hrothgar.h>
#include <hrothgar_proto.h>


void thane_minimize(struct hrothgar_setup *setup)
{
  int i,nreceived;
  MPI_Status status;
  double *params,*output,linpar,minchi;
  gsl_vector_view v,w;
  gsl_min_fminimizer *minimizer;

  minimizer = gsl_min_fminimizer_alloc(gsl_min_fminimizer_brent);

  init_minimize_workspace(setup);
  
  params = sci_dvector(2.*setup->ntotparams);
  output = sci_dvector(setup->ndata);
  v = gsl_vector_view_array(setup->fitparams,setup->nfitparams);
  w = gsl_vector_view_array(output,setup->ndata);

  while (TRUE) {
    MPI_Recv(params,2.*setup->ntotparams,
	     MPI_DOUBLE,0,MPI_ANY_TAG,MPI_COMM_WORLD,&status);
    
    for (i = 0; i < setup->nfitparams; ++i)
      setup->fitparams[i] = params[setup->far[i]];
    
    if (status.MPI_TAG == HROTHGAR_QUIT)
      return;
    if (status.MPI_TAG == HROTHGAR_RESET)
      MPI_Send(output,setup->ndata,MPI_DOUBLE,0,
	       HROTHGAR_RESET,MPI_COMM_WORLD);
    else {
      MPI_Get_count(&status,MPI_DOUBLE,&nreceived);
      
      if (nreceived == setup->ntotparams+1) {

	// The thane is calculating a function and returning it
	setup->dojacobian = (params[setup->ntotparams] < 0);
	minimize_func_func(&v.vector,setup,&w.vector);
	MPI_Send(output,setup->ndata,MPI_DOUBLE,0,
		 status.MPI_TAG,MPI_COMM_WORLD);
      } else {
	// The thane is doing line minimization.
	
	setup->step = &params[setup->ntotparams];
	
	line_minimize(setup,minimizer,&linpar,&minchi);

	// Check for error
	if (minchi < 0) 
	  MPI_Send(&minchi,1.,MPI_DOUBLE,0,HROTHGAR_RESET,
		   MPI_COMM_WORLD);
	else {

	  //if (fabs(linpar) < 1.e-8 && setup->debug) {
	  //printf("Node %d: Minimum not changed for step vector:\n",
	  //   setup->node);
	  //for (i = 0; i < setup->nfitparams; ++i)
	  //  printf("%8s %11.1E",
	  //     setup->paramnames[setup->far[i]],setup->step[i]);
	  //printf("\n");
	  //fflush(stdout);
	  //	  }

	  for (i = 0; i < setup->nfitparams; ++i)
	    setup->fitparams[i] += linpar*setup->step[i];
	  
	  setup->fitparams[setup->nfitparams] = minchi;

	  MPI_Send(setup->fitparams,setup->nfitparams+1,MPI_DOUBLE,0,
		   status.MPI_TAG,MPI_COMM_WORLD);

	}
      }
    }

  }
  free(params);
  free(output);
}


void multiproc_mpi(double **input, double **output, int nunits, 
		   struct hrothgar_setup *setup)
{
  int i,j,cpu=0;
  int nrec = 0, unit;
  MPI_Status status;

  for (i = 1; i < setup->ncpu && i <= nunits; ++i) {
    cpu = i;
    input[i-1][setup->ntotparams] = (i > 1 ? -1. : 1.);
    MPI_Send(input[i-1],setup->ntotparams+1,MPI_DOUBLE,cpu,
	     i-1,MPI_COMM_WORLD);
  }

  while (nrec < nunits) {

    MPI_Recv(setup->allchi,setup->ndata,MPI_DOUBLE,MPI_ANY_SOURCE,
	     MPI_ANY_TAG,MPI_COMM_WORLD,&status);
    
    unit = status.MPI_TAG;
    ++nrec;

    if (setup->nloop && (nrec % 100) == 0) {
      printf("%6.2f%% done.\b\b\b\b\b\b\b\b\b\b\b\b\b",
	     (float)100.*nrec/setup->nloop);
      fflush(stdout);
    }
      

    if (output)
      for (j = 0; j < setup->ndata; ++j) 
	output[unit][j] = setup->allchi[j];

    if (i <= nunits) {
      input[i-1][setup->ntotparams] = -1.;

      //if (cpu == setup->ncpu-1) cpu = 1; else ++cpu;
      //MPI_Send(input[i-1],setup->ntotparams+1,MPI_DOUBLE,
      //       cpu,i-1,MPI_COMM_WORLD);

      MPI_Send(input[i-1],setup->ntotparams+1,MPI_DOUBLE,
	       status.MPI_SOURCE,i-1,MPI_COMM_WORLD);
      ++i;
    }
  }

}

void lord_confoutput(struct hrothgar_setup *setup)
{
  int i,j,k,cpu,nsent=0,ii,jj;
  double cinfo[NCINFO],**limits,lims[6];
  MPI_Status status;

  limits = sci_dmatrix(setup->nfitparams,6);

  for (i = 0; i < setup->nfitparams; ++i) {
    limits[i][1] = -1.E30;
    limits[i][2] = 1.E30;
  }

  for (k = 1; k >= 0; --k) {
    for (i = 0; i < setup->nfitparams-1; ++i)
      for (j = i+1; j < setup->nfitparams; ++j) {
	initialize_confidence(cinfo,i,j,setup);
	if (k < 1) {
	  cinfo[2] = limits[i][1];
	  cinfo[3] = limits[i][2];
	  cinfo[4] = limits[j][1];
	  cinfo[5] = limits[j][2];
	}
	cinfo[11] = k;
	
	if (nsent >= setup->ncpu-1) {
	  MPI_Recv(lims,6,MPI_DOUBLE,MPI_ANY_SOURCE,MPI_ANY_TAG,
		   MPI_COMM_WORLD,&status);
	  
	  ii = status.MPI_TAG / setup->nfitparams;
	  jj = status.MPI_TAG - ii*setup->nfitparams;
	  
	  adjust_limits(limits,ii,jj,lims,k);
	  
	  cpu = status.MPI_SOURCE;
	  --nsent;
	} else
	  cpu = nsent+1;
	
	MPI_Send(cinfo,NCINFO,MPI_DOUBLE,cpu,i*setup->nfitparams+j,MPI_COMM_WORLD);
	++nsent;
      }
    
    while (nsent > 0) {
      MPI_Recv(lims,6,MPI_DOUBLE,MPI_ANY_SOURCE,MPI_ANY_TAG,
	       MPI_COMM_WORLD,&status);
      
      ii = status.MPI_TAG / setup->nfitparams;
      jj = status.MPI_TAG - ii*setup->nfitparams;
   
      adjust_limits(limits,ii,jj,lims,k);
      
      --nsent;
    }
  }

  setup->outfile = fopen(setup->outfilename,"a");
  if (!setup->outfile) BYE("Error appending to %s.",setup->outfilename);
  for (i = 0; i < setup->nfitparams; ++i)
    confidence_report(i,limits[i],setup);

  fclose(setup->outfile);
  
  reset_all_nodes(setup);

  sci_free_dmatrix(limits,setup->nfitparams);
  //reset_all_nodes(setup);
}

void thane_confoutput(struct hrothgar_setup *setup)
{
  double limits[6];
  double cinfo[NCINFO];
  MPI_Status status;

  while (TRUE) {
    MPI_Recv(cinfo,NCINFO,MPI_DOUBLE,0,MPI_ANY_TAG,
	     MPI_COMM_WORLD,&status);

    if (status.MPI_TAG == HROTHGAR_QUIT) break;

    confidence_2d(cinfo,limits,setup);

    MPI_Send(limits,6,MPI_DOUBLE,0,status.MPI_TAG,
	     MPI_COMM_WORLD);
  }

}


void reset_all_nodes(struct hrothgar_setup *setup)
{
  int i;

  for (i = 1; i < setup->ncpu; ++i)
    MPI_Send(&setup->params[0],1,
	     MPI_DOUBLE,i,HROTHGAR_QUIT,MPI_COMM_WORLD);

}

void mcmc_stats(double *meanlist, double *m2list, double *outmean,
		double *outm2, int nchunk)
{
  int i;
  double frac,delta;

  outmean[nchunk-1] = meanlist[nchunk-1];
  outm2[nchunk-1] = m2list[nchunk-1];

  for (i = nchunk-1; i > 0; --i) {

    delta = outmean[i]-meanlist[i-1];
    frac = (double)(nchunk-i)/(nchunk-i+1);
    outmean[i-1] = meanlist[i-1]+delta*frac;
    outm2[i-1] = m2list[i-1]+outm2[i]+delta*delta*CHUNKSIZE_MCMC*frac;

  }
}

#define MPI_SEED_TAG 1
#define MPI_DATA_TAG 2
#define MPI_CONVERGENCE_TAG 3
  
void lord_mcmc(struct hrothgar_setup *setup)
{
  long seed;
  int i,j,k,l,foundbetter,nreceived=0,nconverged[MAX_CHUNKS_MCMC];
  int burnin=0;
  double minchisq=1.E30;
  MPI_Status status;
  MPI_Request *request;

  // Send a different random number seed to each thane
  for (i = 1; i < setup->ncpu; ++i) {
    seed=0;
    while (seed < 1)
      seed = gsl_rng_get(setup->generator);
    MPI_Send(&seed,1,MPI_LONG,i,MPI_SEED_TAG,MPI_COMM_WORLD);
  }

  request = (MPI_Request *)malloc(setup->ncpu*sizeof(MPI_Request));

  double ***mean = sci_dtensor(setup->nfitparams,setup->ncpu,MAX_CHUNKS_MCMC);
  double ***m2 = sci_dtensor(setup->nfitparams,setup->ncpu,MAX_CHUNKS_MCMC);
  double **bestfit = sci_dmatrix(setup->ncpu,setup->nfitparams);
  double ***running_mean = sci_dtensor(setup->nfitparams,setup->ncpu,MAX_CHUNKS_MCMC);
  double ***running_sig = sci_dtensor(setup->nfitparams,setup->ncpu,MAX_CHUNKS_MCMC);
  double *mcmcdata = sci_dvector(3*setup->nfitparams+1);
  int **converged = sci_imatrix(MAX_CHUNKS_MCMC,setup->nfitparams);
  int **convsend = sci_imatrix(setup->ncpu,setup->nfitparams);
  int nchunk[MAX_CHUNKS_MCMC],minchunk=0,maxconverged=0;
  int *initial = sci_ivector(setup->ncpu);
  int *sentrequest = sci_ivector(setup->ncpu);
  double delta,minchi[MAX_CHUNKS_MCMC];

  // Now check for convergence of the MCMC chains
  int nexpected = 3*setup->nfitparams+1;

  for (i = 0; i < MAX_CHUNKS_MCMC; ++i) {
    nchunk[i] = 0;
    minchi[i] = 0;
  }

  while ((maxconverged < setup->nfitparams ||
	  minchunk*CHUNKSIZE_MCMC < 1.*setup->nsim/(setup->ncpu-1)) &&
	 minchunk*CHUNKSIZE_MCMC < 5.*setup->nsim/(setup->ncpu-1)){
    
    //printf("\n\nAwaiting data..\n");
    MPI_Recv(mcmcdata,nexpected,MPI_DOUBLE,
	     MPI_ANY_SOURCE,MPI_DATA_TAG,MPI_COMM_WORLD,&status);

    MPI_Get_count(&status,MPI_DOUBLE,&nreceived);
    int n = status.MPI_SOURCE;
    //printf("Received data from %d\n",n);
    fflush(stdout);
    if (nreceived != nexpected) 
      hrothgar_printf("Warning: Node %d return %d/%d items",n,nreceived,
		      nexpected);

    k = 0;
    ++nchunk[n];
    minchunk = nchunk[n];

    for (j = 0; j < setup->nfitparams; ++j) {
      mean[j][n][nchunk[n]-1] = mcmcdata[j];
      m2[j][n][nchunk[n]-1] = mcmcdata[j+setup->nfitparams];
      bestfit[n][j] = mcmcdata[j+2*setup->nfitparams];
    }
    minchi[n] = mcmcdata[3*setup->nfitparams];

    // Find the thane with the smallest chunk size.
    foundbetter=0;
    for (i = 1; i < setup->ncpu; ++i) {
      if (nchunk[i] < minchunk) minchunk = nchunk[i];
      if (minchi[n] < minchisq) {
	minchisq = minchi[n];
	foundbetter = 1;
      }
    }

    if (foundbetter)
      report_unmapped(bestfit[n],minchisq,setup);

    // Now calculate the mean and (n-1)*variance in cumulative bins,
    // starting with all bins and moving to just the final bin.
    if (nchunk[n] > MIN_CHUNKS_MCMC) 
      for (j = 0; j < setup->nfitparams; ++j) {
	mcmc_stats(mean[j][n],m2[j][n],running_mean[j][n],running_sig[j][n],
		   nchunk[n]);
	for (k = 0; k < nchunk[n]; ++k)
	  running_sig[j][n][k] = sqrt(running_sig[j][n][k]/
				      (CHUNKSIZE_MCMC*(nchunk[n]-k)));
      }

    maxconverged = 0;
    if (minchunk > MIN_CHUNKS_MCMC) { 
      
      for (l = 0; l < minchunk-MIN_CHUNKS_MCMC; ++l) {
	nconverged[l] = 0;
	for (j = 0; j < setup->nfitparams; ++j)  {
	  converged[l][j] = 1;
	  for (i = 1; i < setup->ncpu; ++i)
	    for (k = i; k < setup->ncpu; ++k) {
	      delta = fabs((running_mean[j][i][l]-running_mean[j][k][l])/
			   (running_sig[j][i][l]+running_sig[j][k][l]));
	      converged[l][j] *= (delta < 1);
	    }
	  if (converged[l][j]) ++nconverged[l];
	}
	if (nconverged[l] > maxconverged) {
	  burnin = l;
	  maxconverged = nconverged[l];
	}
      }

      if (setup->debug)
	for (j = 0; j < setup->nfitparams; ++j)  {
	  printf("%-10s %3d ",setup->paramnames[setup->far[j]],
		 converged[burnin][j]);
	  for (k = 1; k < setup->ncpu; ++k)
	    printf("%8.1E %8.1E | ",running_mean[j][k][burnin],
		   running_sig[j][k][burnin]);
	  printf("\n");
	}

      if (initial[n])
	MPI_Wait(&request[n],&status);
      initial[n] = 1;
      for (j = 0; j < setup->nfitparams; ++j)
	convsend[n][j] = converged[burnin][j];
      MPI_Isend(convsend[n],setup->nfitparams,MPI_INT,n,MPI_CONVERGENCE_TAG,
		MPI_COMM_WORLD,&request[n]);
      sentrequest[n] = 1;
    } else
      sentrequest[n] = 0;
  }
    

  for (i = 1; i < setup->ncpu; ++i) {
    if (sentrequest[i])
      MPI_Wait(&request[i],&status);
    converged[burnin][0] = -1;
    MPI_Send(converged[burnin],setup->nfitparams,MPI_INT,i,MPI_CONVERGENCE_TAG,
	     MPI_COMM_WORLD);
  }
  
  setup->ntrain = burnin;
  setup->outfile = fopen(setup->outfilename,"a");
  fprintf(setup->outfile,"# Burnin: %ld\n",(long)burnin*CHUNKSIZE_MCMC);
  fprintf(setup->outfile,"# Gaussfrac: %E\n",setup->gaussfrac);
  fclose(setup->outfile);
}

void thane_mcmc(struct hrothgar_setup *setup)
{
  double *mcmcdata;
  MPI_Request requestin, requestout;
  MPI_Status status;
  int requestisin, stoprequest=0;
  long seed;
  int i;

  int ndata_mpi = 3*setup->nfitparams+1;

  mcmcdata = sci_dvector(ndata_mpi);
  int *converged = sci_ivector(setup->nfitparams);

  MPI_Recv(&seed,1,MPI_LONG,0,MPI_SEED_TAG,MPI_COMM_WORLD,&status);

  if (seed > 1)
    gsl_rng_set(setup->generator,seed);
  else
    BYE("Nonpositive seed received from Lord.\n");

  int initial=1;
  requestisin = 1;
  while (!stoprequest) {
    
    if (requestisin)  {
      if (!initial) {
	for (i = 0; i < setup->nfitparams; ++i)
	  if (!converged[i]) 
	    setup->mcmcstep[setup->far[i]] *= 1.01;
      }
      MPI_Irecv(converged,setup->nfitparams,MPI_INT,0,
		  MPI_CONVERGENCE_TAG,MPI_COMM_WORLD,&requestin);
    }

    if (!initial) {
      if (setup->debug)
	printf("\nNode %d waiting for input.\n",setup->node);
      MPI_Wait(&requestout, &status);
      if (setup->debug)
	printf("\nNode %d received.\n",setup->node);
    }

    hrothgar_mcmc(initial,CHUNKSIZE_MCMC,mcmcdata,setup);
    
    MPI_Isend(mcmcdata,ndata_mpi,MPI_DOUBLE,0,MPI_DATA_TAG,
	      MPI_COMM_WORLD,&requestout);
    
    MPI_Test(&requestin,&requestisin,&status);
    stoprequest = (requestisin && converged[0] == -1);
    initial = 0;

  }

}

void hrothgar_mpi(struct hrothgar_setup *setup)
{
 
  int rank;
  int ncpu;

  ncpu = setup->ncpu;
  rank = setup->node;

  if (ncpu == 1) {
    hrothgar_singlecpu(setup);
    return;
  }

  setup->multifunc = multiproc_mpi;

  MPI_Barrier(MPI_COMM_WORLD);

  if (rank == 0) {
    printf("Initialization complete.\n");
    if (setup->nsim) 
      lord_mcmc(setup);
    else {
      lord_minimize(setup);
      reset_all_nodes(setup);

      MPI_Barrier(MPI_COMM_WORLD);

      if (!setup->nloop && !setup->resume)
	lord_confoutput(setup);

    }
    reset_all_nodes(setup);
    
  } else {
    if (setup->nsim) 
      thane_mcmc(setup);
    else {
      thane_minimize(setup);
      
      MPI_Barrier(MPI_COMM_WORLD);

      if (!setup->nloop && !setup->resume)
	thane_confoutput(setup);

    }
  }

  MPI_Barrier(MPI_COMM_WORLD);
  MPI_Finalize();

}

void init_mpi(int *ntotcpu, int *rank)
{
    MPI_Init(NULL,NULL);
    MPI_Comm_size(MPI_COMM_WORLD,ntotcpu);
    MPI_Comm_rank(MPI_COMM_WORLD,rank);
}
