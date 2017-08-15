#include <gsl/gsl_integration.h>
#include <gsl/gsl_multifit.h>
#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <math.h>
#include <time.h>
#include <sciutils.h>

#define PI       3.14159265358979
#define INTMAX   1000000

// Layout of PSF File.
// If we have N annuli, N^2 lines are output.
// The first N lines should be the light leaked from annulus 1 into
// the nth annulus.

double Rpint,Rint,Rpint1,Rpint2,r0,r0sq,b;
double Rpintsq,Rintsq,RRp,r2sq;
gsl_integration_workspace *intwork1, *intwork2, *intwork3; 

gsl_function func_psf,func_integ1,func_integ2;

double psf(double theta, void *params __attribute__ ((unused)))
{
  return pow(r0sq+Rintsq+Rpintsq - 2.*RRp*cos(theta),-b);
}

double integ1(double Rp, void *params __attribute__ ((unused)))
{
  double ans1,ans2,err;

  Rpint = Rp;
  Rpintsq = Rp*Rp;
  RRp = Rint*Rp;
  gsl_integration_qag(&func_psf,0,PI/2.,0.0,1e-5,INTMAX,GSL_INTEG_GAUSS15,
		      intwork1,&ans1,&err);
  gsl_integration_qag(&func_psf,PI/2.,PI,0.0,1e-5,INTMAX,GSL_INTEG_GAUSS15,
		      intwork1,&ans2,&err);

  return 2.*Rp*(ans1+ans2);
}


// Calculate leaked light from Rpint1--Rpint2 into R
double integ2(double R, void *params __attribute__ ((unused)))
{
  double ans,err;

  Rint = R;
  Rintsq = R*R;
  
  gsl_integration_qag(&func_integ1,Rpint1,Rpint2,0.0,1.e-4,INTMAX,GSL_INTEG_GAUSS15,
		      intwork2,&ans,&err);
  return R*ans;
}

// Calculate leaked light from Rpint1--Rpint2 into R1--R2
double getfrac(double Rp1, double Rp2, double R1, double R2)
{
  double ans, err;

  Rpint1 = Rp1;
  Rpint2 = Rp2;

  gsl_integration_qag(&func_integ2,R1,R2,0.0,1.e-3,INTMAX,GSL_INTEG_GAUSS15,
		      intwork3,&ans,&err);
  //return qromo(integ2,R1,R2,midpnt);
  return ans;
}

// Calculate normalization of leaked light from Rp1--Rp2
float gettot(float Rp1, float Rp2)
{
  double ans1,ans2, err;

  Rpint1 = Rp1;
  Rpint2 = Rp2;

  gsl_integration_qag(&func_integ2,0.0,10,0.0,1.e-3,INTMAX,GSL_INTEG_GAUSS15,
		      intwork3,&ans1,&err);
  gsl_integration_qag(&func_integ2,10,20.,0.0,1.e-3,INTMAX,GSL_INTEG_GAUSS15,
		      intwork3,&ans2,&err);
  return ans1+ans2;
    //return qromo(integ2,0.,r0,midpnt)+qromo(integ2,r0,30.,midpnt);
}

int main(int argc, char *argv[])
{
  unsigned long i,j,k,count,ne;
  double *Rp1,*Rp2;
  double totflux[10],rmean,ene,fac;
  double core, recoef, rrcoef, rrecoef;
  double alpha, aecoef, arcoef, arecoef, chisq;

  gsl_matrix *X,*cov;
  gsl_vector *y, *c;

  intwork1 = gsl_integration_workspace_alloc(INTMAX);
  intwork2 = gsl_integration_workspace_alloc(INTMAX);
  intwork3 = gsl_integration_workspace_alloc(INTMAX);

  if (argc < 5) {
    printf("PSF by A. Mahdavi, April 2004\n\n");
    printf("Usage:\n");
    printf("psf file col1 col2 m1|m2|pn|ai|as\n\n");
    printf("Given a file with a set of starting and ending annulus radii\n");
    printf("in columsn col1 and col2, and an instrument id, this program\n");
    printf("calculates the contribution of the PSF at each bin.\n");
    return 0;
  }

  Rp1 = double_readdata(argv[1],atoi(argv[2]),&count,0);
  Rp2 = double_readdata(argv[1],atoi(argv[3]),&count,0);

  if (strstr(argv[4],"pn") != NULL) {
    core = 6.636; recoef = -0.305; rrcoef = -0.175; rrecoef = -0.0067;
    alpha = 1.525; aecoef = -0.015; arcoef = -0.012; arecoef = -0.001;
  }
  else if (strstr(argv[4],"m1") != NULL) {
    core = 5.074; recoef = -0.236; rrcoef = 0.002; rrecoef = -0.018;
    alpha = 1.472; aecoef = -0.01; arcoef = -0.001; arecoef = -0.0016;
  }
  else if (strstr(argv[4],"m2") != NULL) {
    core = 4.759; recoef = -0.203; rrcoef = 0.014; rrecoef = -0.0229;
    alpha = 1.411; aecoef = -0.005; arcoef = -0.001; arecoef = -0.0002;
  }
  else if (strstr(argv[4],"ai") != NULL || strstr(argv[4],"as") != NULL) {
    core = 1.33695; recoef = -0.05347; rrcoef = 0.076318; rrecoef = 0.0340592;
    alpha = 6.1184; aecoef = -0.25421; arcoef = -0.65243; arecoef = 0.0508094;
    
  }
  else {
    printf("You must specify m1, m2, or pn.\n");
    return -1;
  }

  func_psf.function = psf;
  func_integ1.function = integ1;
  func_integ2.function = integ2;

  ne = 6;
  double energy[ne];
  energy[0] = 0.3;
  energy[ne-1] = 10.;
  fac = pow(energy[ne-1]/energy[0],1./(ne-1.));

  X = gsl_matrix_alloc(ne, 3);
  y = gsl_vector_alloc(ne);
  c = gsl_vector_alloc(3);
  cov = gsl_matrix_alloc(3,3);

  for (k = 1; k < ne-1; ++k) 
    energy[k] = fac*energy[k-1];
  for (k = 0; k < ne; ++k) {
    gsl_matrix_set(X,k,0,1.0);
    gsl_matrix_set(X,k,1,energy[k]);
    gsl_matrix_set(X,k,2,energy[k]*energy[k]);
  }
  gsl_multifit_linear_workspace *work = gsl_multifit_linear_alloc(ne,3);

  FILE *status;
  unsigned long totout=0;
  double tf;
  status = fopen("/dev/stderr","w+");

  // This loop corresponds to the FROM annulus
  for (i = 0; i < count; ++i) {
    tf = 0;

    rmean = (Rp1[i]+Rp2[i])/2.;
    // This loop corresponds to the TO annulus
    for (k = 0; k < count; ++k) {
      for (j = 0; j < ne; ++j) {

	ene = energy[j];
	r0 = core+recoef*ene+rmean*(rrcoef+rrecoef*ene);
	b = alpha+aecoef*ene+rmean*(arcoef+arecoef*ene);
	r0 /= 60;
	r0sq = r0*r0;
	if (k == 0)
	  totflux[j] = gettot(Rp1[i],Rp2[i]);

	
	fac = getfrac(Rp1[i],Rp2[i],Rp1[k],Rp2[k]);
	if (j == 0) {
	  	tf += fac;

		//printf("%ld %ld %E %E %E %E\n",i,k,fac,totflux[j],fac/totflux[j],tf/totflux[j]);
	}
	gsl_vector_set(y,j,fac/totflux[j]);
      }
      gsl_multifit_linear(X,y,c,cov,&chisq,work);

      totout = 0;
      if (chisq > 4)
	totout += 
	  fprintf(status,"Warning---accuracy %f, %E %E %E.\n",chisq,
		  gsl_vector_get(y,0),gsl_vector_get(y,1),
		  gsl_vector_get(y,2));
      printf("%12.3E %12.3E %12.3E\n",
	     gsl_vector_get(c,0),gsl_vector_get(c,1),gsl_vector_get(c,2));
    }
    totout += fprintf(status," %.1f%% done.",100.*i/count);
    fflush(status);
    fflush(stdout);
    for (k = 0; k < totout; ++k) fprintf(status,"\b");
  }
  gsl_multifit_linear_free(work);
  gsl_matrix_free(X);
  gsl_vector_free(y);
  gsl_vector_free(c);
  fprintf(status,"\n");
  fclose(status);
  return 0;
}
// 0.413807 +/-   0.000411
//Parameter 02:   0.006860 +/-   0.000266
//Parameter 03:   0.000249 +/-   0.000030
