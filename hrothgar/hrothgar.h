
/* hrothgar.h
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

#ifdef HAVE_CONFIG_H
#include <config.h>
#endif

#include <gsl/gsl_matrix.h>
#include <gsl/gsl_vector.h>
#include <zlib.h>

struct hrothgar_setup {

  // General variables
  char **argv;            // Command line past from main program
  int argc;               // Number of command line parameters
  int ncpu;               // total # of CPUs to use in calculations
  int node;               // # of this cpu.
  int messages;           // Whether to output non-error messages
  int ntotparams;         // Total number of parameters
  char **paramnames;      // The names of the parameters
  double *initvals;       // Their default initival values
  int  *dofit;            // Their default fit/non-fit state
  char **pcomment;        // Comment describing parameters
  long ndata;             // Total number of data points
  double *x;              // Data (x) variables
  double *y;              // Data (y) variables
  double *errors;         // Data errors on (y)
  double *model;          // Model calculated on (x,y,errors)
  int *frozen;            // TRUE if a parameter is not to be varied
  int *origfrozen;        // Original copy of frozen
  int overwrite;          // Whether we overwrite existing files
  int makecontours;       // Whether to generate sets of 2D contours
  int nfitparams;         // Number of free parameters (determined from frozen)
  double *bestfit;        // Set of best-fit parameters
  int *far;               // far[i] gives the ith free parameter
  double *params;         // The chain's starting point
  double *oldpmin,*pmin;  // Minimum allowed values for the parameters
  double *oldpmax,*pmax;  // Maximum allowed values for the parameters
  int *varyoften;         // Whether to vary a parameter often in MCMC mode
  long timer;             // Whether to time the performance of the run
  int doui;               // Whether to present the user with an interface
  int  debug;             // Output diagnostics
  long nloop;             // Number of loops to conduct after confidence calculations
  char *welcome_string;   // Welcome string to print at beginning
  int  nstringpars;       // Number of string-valued extra parameters
  char **stringparname;   // Their names
  char **stringparval;    // Their values
  char **stringparcomment;// Their comments
  char *infilename;          
  char *outfilename;          
  FILE *outfile;  
  gzFile mcmcfile;
  char   *mcmcfilename;

  // Variables related to MCMC
  int resume;             // Whether we are resume a previous simulation
  int predictable;        // Whether same-seed runs should be consistent
  unsigned long seed;     // Random number seed
  long nsim;              // Number of simulations 0=minimize instead of chain
  long ndone;              
  long nvary;             // Number of parameters to vary during normal phase
  long ntrain;            // Number of Burn-in simulations
  long nvary_train;       // Number of parameters to vary during burn-in
  double gaussfrac;       // Fractional gaussian displacement (1-sigma)
  double logprior;        // Natural Log of Bayesian prior
  long   nsteps;          // Total # of steps to take in each iteration
  int    *ntested;        // Number of points within cycle tested
  int    *naccepted;      // Number of points within cycle accepted
  double acceptrate; 
  int    *mcmcgood;      
  int    ngood;
  long  ntotaccepted;     // Total # of points within cycle accepted
  int    *picked;         // List of variables picked for current step
  double *mcmcstep;       // Step sizes in each direction for mcmc.
  void *generator;        // The GSL random number generator; 
                          // defined in gsl/gsl_rng.h

  // Extra parameters to be tracked along fit parameters
  int ninfo;              // # of info parameters
  char **infoname;        // Names of info parameters
  double *infoarray;      // List of info values.

  // Covariance matrix parameters
  gsl_matrix *covar;      // Covariance matrix in parameter coordinates
  gsl_matrix *covarm;     // Covariance matrix in parameter coordinates
  gsl_matrix *covarp;     // Covariance matrix in parameter coordinates
  gsl_matrix *tcovar;     // Covariance matrix in sin(p) coordinates
  double **stepmatrix;    // Rotation matrix for MCMC
  double *evalues;        // Eigenvalues of covariance matrix
  double **chain;         // Recent MCMC chain for covariance matrix
  int nchain;             // Keeps track of position within covar mtrx chain

  // Variables related to minimization
  int covaronly;          // Calculate covariance matrix only?
  int evalonly;           // Show merit function value only?
  int ignorecovar;        // Ignore covariance matrix even if it ex
  int threepoint;         // Whether to use the 3 point rule for jacobians
  int lognormal;           // Whether to use include lognormal scatter
  int thorough;           // Whether to make thorough Powell search
  int remember;           // Remember frozen/unfrozen state
  int floataround;
  int dojacobian; 
  double *step;           // Powell minimization step vector
  double conf;            // Which 1D confidence level to report
  double inputaccuracy;   // Accuracy of input model
  double stepsize;        // Stepsize for derivatives = sqrt(inputaccuracy)
  double *h;              // Array of stepsizes for calculating derivatives
  double eps;             // fractional accuracy for minimization
  double *sineparams;     // Holds sine-transformed fit parameters
  double width;           // width of confidence plot in sigma
  int ngrid;              // # of grid points for each axes in plot
  int cumulative;         // Whether to output cumulative PDF
  int untransformed;      // Whether to calculate an untransformed Jacobian
  void *s;                // LM Solution vector
  int iterations;         // Number of LM iterations before convergence
  int ncheb;
  int nchebtot;
  int whichpar;
  int whichdata;
  double chimin;          // Minimum overall chi square
  double penalty;         // Penalty to apply to overall chisq.

  void (*multifunc)(double **input, double **output, int units, 
		    struct hrothgar_setup *setup);

  int (*get_model)(double *x, double *params, double *model, double *errors,
		   double *logprior,
		   unsigned long ndata, void *dataptr);
  double (*statistic)(struct hrothgar_setup *setup);

  // Use stat-only method (i.e. function just gives -2 ln P instead
  // full info?
  int  statonly;          // Whether to employ stat-only method;
  // Stat-only function
  double (*get_stat)(double *params, int np, void *dataptr);

  // This represents the default minimizer. Right now the available options
  // are a scaled or a non-scaled Levenberg-Marquard solver.
  const void *solvertype; 

  // The following are workspaces for jacobian calculations
  double **wparams,**bigmatrix,*allchi,*fitparams;
  double **jacobian1,**jacobian2,**jacobian3,**jacobian4;
  double **rounderr,**truncerr,**hopt;
  void    *Jerr;

  // This is a pointer to external data
  void *data;             

};

struct treenode {

  int cpu;                            // The cpu that is processing this node
  int depth;                          // The depth of the tree node
  int completed;                      // Whether the node is done processing
  int accepted;                       // Whether the node is accepted in the chain
  int ntimes;                         // Number of times node appears in chain
  double *params;                     // Parameter set within the node
  double chisq;                       // Merit value
  struct treenode *left;              
  struct treenode *right;
  struct treenode *sibling;
  struct treenode *parent;

};

struct minstruct {

  int ntotparams;         // Total number of parameters
  int  *frozen;           // Wheter a parameter is frozen
  double *params;         // Placeholder for parameter values
  double *point;          // Current point being considered
  double *step;           // Offset from this point
  void   *data;           // Data to be passed on to chi^2 function


};

static const char long_gpl[] = 

"\n    This program is free software; you can redistribute it and/or modify\n    it under the terms of the GNU General Public License as published by\n    the Free Software Foundation; either version 2 of the License, or\n    (at your option) any later version.\n\n    This program is distributed in the hope that it will be useful,\n    but WITHOUT ANY WARRANTY; without even the implied warranty of\n    MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the\n    GNU General Public License for more details.\n\n    You should have received a copy of the GNU General Public License along\n    with this program; if not, write to the Free Software Foundation, Inc.,\n    51 Franklin Street, Fifth Floor, Boston, MA 02110-1301 USA.\n\n";

double get_chisq(double *params, void *data, double *outerrs,
		 double *outmodels,
		 double *chiarr, long *arrsize);

#ifndef DEFAULT_MCMC_GENERATOR
#define DEFAULT_MCMC_GENERATOR gsl_rng_mt19937
#endif

#define MAX(x,y) ((x) > (y) ? (x) : (y))
#define MIN(x,y) ((x) < (y) ? (x) : (y))

#define GOLDEN 0.38197
#define NCINFO 13
#define TRUE       1
#define FALSE      0
#define IGNORE    -1
#define HROTHGAR_WORK  (INT_MAX-2)
#define HROTHGAR_RESET (INT_MAX-1)
#define HROTHGAR_QUIT  (INT_MAX)

#define BYE(format...) { printf("Hrothgar Error:\n"); printf(format); printf("\n"); exit(-1); }

#define PI 3.14159265358979

#ifndef HROTHGAR_PROTO
void hrothgar_init_welcomestring(struct hrothgar_setup *setup,
				 char *welcome_string);
void hrothgar_init_stringpars(struct hrothgar_setup *setup,
			      int nstringpars, char **stringparname, 
			     char **stringparval, char **stringparcomment);
void hrothgar_init_pars(struct hrothgar_setup *setup,
			int ntotparams, char **paramnames, 
			double *pmin, double *pmax, 
			double *initvals, int *dofit,char **pcomment);
int hrothgar_init(struct hrothgar_setup *setup, 
		  int argc, char *argv[]);

void hrothgar_statonly(double (*get_stat)(double *params, int np, 
					  void *dataptr),
		       void *dataptr,
		       struct hrothgar_setup *setup);

void hrothgar(unsigned long ndata, double *x, double *data, double *error,
	      int (*get_model)(double *x, double *params, double *model, double *errors, double *prior,
			       unsigned long ndata, void *dataptr),
	      void *dataptr,
	      struct hrothgar_setup *setup);
#endif

#ifndef hrothgar_printf
#define hrothgar_printf(format...) if (setup->node < 2 && setup->messages) printf(format); 
#endif


#define CHUNKSIZE_MCMC 1000
#define MAX_CHUNKS_MCMC 10000
#define MIN_CHUNKS_MCMC 5
