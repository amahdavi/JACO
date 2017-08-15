#include <gsl/gsl_rng.h>
#include <gsl/gsl_randist.h>
#include <gsl/gsl_sort.h>
#include <hrothgar.h>
#include <math.h>

// HAVE_MPI is automatically defined in the build process for Hrothgar.
#ifdef HAVE_MPI
#include <mpi.h>
#endif

// This function simply returns a Gaussian distribution for 
// ndata values of x, storing them in model. par[0] is the mean,
// par[1] is the sigma, and par[2] is the normalization of the
// gaussian.
// The user can use this function as a template for building fits.

int gaussian_model(double *x, double *pars, double *model, 
		   double *error, double *logprior,
		   unsigned long ndata, void *extra){

  double diff;
  unsigned long i;

  for (i = 0; i < ndata; ++i) {
    diff = (x[i]-pars[0])/pars[1];
    model[i] = pars[2]*exp(-diff*diff/2.)/(pars[1]*2.5066283);
  }

  return 0;

}

/* double gaussian_model_stat(double *pars, int npar, void *extra) */
/* { */
/*   double diff; */
/*   diff = pow((1 - pars[0]),2.)+pow((2-pars[1]),2.)+pow((0.5-pars[2]),2.); */
/*   return diff; */
/*  } */

int main(int argc, char *argv[]) {

  unsigned long i,j;
  double tempmean;

  // This is a required declaration.
  struct hrothgar_setup setup;

  // There are 3 total fit parameters: the mean, sigma, and 
  // normalization of the Gaussian
  // Parameter names as they should appear in the default config file
  static char *paramnames[] = { "mean", "sigma", "norm" };

  // Which parameters should be frozen? None, in our case, we are 
  // fitting them all. This is a default only, and can be overridden.
  static int frozen[] = { 0, 0, 0 };

  // Optional comments describing what each parameter means: 
  static char *comments[] = { "Mean", "Sigma",
			    "Normalization" };

  // Use remote initial values for the minimization.
  static double initvalues[] = { -1, 3., 23. };
  
  // Minimum and maximum values allowed for minimization
  static double parmin[] = { -10., 0.001, 0.0001 };
  static double parmax[] = { 10., 100., 10000. };
  
  // String parameters. These are parameters not used during
  // the fit. In this case, we'll use string parameters to
  // generate our synthetic data.
  
  char *stringparname[] = { "simmean", "simsig", "nsims" };
  
  // By default, we'll simulate a centered Gaussian with unit sigma
  // Note that the parameters need to be string valued here.
  // By default we are doing 1000 simulations. On a modern system,
  // you will need 30000000 (3e7) simulations or more before you will
  // see a big slowdown on the fit.
  char *stringparval[] = { "0.", "1.", "1000" };
  char *stringparcomments[] = { "Simulated mean", 
				"Simulated sigma",
                                "Number of Simulations"};


  /******************** Program starts here **********************/

  // Initialize the string parameters. Pointers to them are stored
  // within the setup structure.
  hrothgar_init_stringpars(&setup,3,stringparname,stringparval,
			   stringparcomments);

  // Initialize the fit parameters 
  hrothgar_init_pars(&setup,3,paramnames,parmin,parmax,
		     initvalues,frozen,comments);

  // Initialize the hrothgar core.
  // Note that any command line modifications to the string 
  // parameters will be written to stringparval at this point
  if (hrothgar_init(&setup,argc,argv) < 0) { return; }

  // Get the mean and sigma of the Gaussian to simulate, as well
  // as the total number of simulations to perform.
  double simmean = atof(stringparval[0]);
  double simsig = atof(stringparval[1]);
  unsigned long nsim = atol(stringparval[2]);

  // Use 30 data points per bin
  if (nsim < 150) BYE("Too few simulations.");
  unsigned long nbins = nsim/30;
  printf("Init %d\n",setup.node);
  

  // Allocate room for the binned data.
  double *x = (double *)malloc(nbins*sizeof(double));
  double *y = (double *)malloc(nbins*sizeof(double));
  double *ye = (double *)malloc(nbins*sizeof(double));

  // Now it is time to generate the data. For simiplicity, we'll have
  // only the master MPI node generate the data. If we're not running
  // MPI, the node number will be 0 by default anyway.

  if (setup.node == 0) {

    // At this point, hrothgar_init has already conveniently 
    // initialized and seeded a random number generator for us.
    gsl_rng *rng = setup.generator;
    double *data = (double *)malloc(nsim*sizeof(double));
    
    // Generate the Gaussians
    for (i = 0; i < nsim; ++i) 
      data[i] = simmean+gsl_ran_gaussian(rng,simsig);
    
    // Sort them
    gsl_sort(data,1,nsim);
    
    // Bin them into a histogram of width 30
    j = 0;
    tempmean = 0.;

    for (i = 0; i < nsim; ++i) {
      tempmean += data[i];
      
      if (i % 30 == 29) {
	x[j] = tempmean/30.;
	tempmean = 0.;
	y[j] = 30./(nsim*(data[i]-data[i-29]));
	ye[j] = 5.5/(nsim*(data[i]-data[i-29]));
	++j;
      }
    }
    free(data);
  }
  

    // Conditional MPI section
#ifdef HAVE_MPI
    // Send all the slaves the binned data.
  printf("Trying %d\n",setup.node);
    MPI_Bcast(x,nbins,MPI_DOUBLE,0,MPI_COMM_WORLD);
    MPI_Bcast(y,nbins,MPI_DOUBLE,0,MPI_COMM_WORLD);
    MPI_Bcast(ye,nbins,MPI_DOUBLE,0,MPI_COMM_WORLD);
#endif
    printf("OK %d\n",setup.node);
  // Run the fit
  hrothgar(nbins,x,y,ye,gaussian_model,NULL,&setup);
  //hrothgar_statonly(gaussian_model_stat,NULL,&setup);

}
