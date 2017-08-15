#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <math.h>
#include <time.h>
#include <unistd.h>
#include <hrothgar.h>
#include <sciutils.h>
#include <gsl/gsl_integration.h>

#define NSTRINGPARAMS 5

char *stringparname[NSTRINGPARAMS] = {
    /* 000 */    "datafile",
    /* 002 */    "lognormal",
    /* 003 */    "plot",
    /* 004 */    "cols",
    /* 005 */    "invfit"
};
  
char *stringparval[NSTRINGPARAMS] = {
    /* 000 */    "fit.txt",
    /* 002 */    "1",
    /* 003 */    "NONE",
    /* 003 */    "2,3,4,5",
    /* 003 */    "0"
};
  
char *stringparcomments[NSTRINGPARAMS] = {
    /* 000 */    "Data (fmt: x xe y ye)",
    /* 001 */    "Use lognormal (1) or linear (0) scatter",
    /* 003 */    "Make PGPLOT?",
    /* 004 */    "List of columns (comma-separated)",
    /* 004 */    "Do an inverse fit?"
};

#ifdef HAVE_LIBCPGPLOT
#include <cpgplot.h>
#endif

#define WELCOMESTRING "Error fitter with scatter v" VERSION "\nCopyright (C) 2011 Andisheh Mahdavi \n"

#define NTOTPARAMS 4

// Model is y = a x^gam + b
static char *paramnames[] = {

 /* 001 */   "slope",
 /* 002 */   "intercept",
 /* 003 */   "intrinsic",
 /* 004 */   "covar",
};

static int dofit[] = {
 /* 001 */   0,
 /* 002 */   0,
 /* 003 */   0,1};

static char *paramunits[] = {
 /* 001 */   "normalization/slope",
 /* 002 */   "intercept",
 /* 003 */   "Fractional intrinsic scatter in y at constant x","Covariance of errors in x and y"};

static double initvalues[] = {
 /* 001 */   1.0,
 /* 002 */   1.0,
 /* 003 */   0.1,0.};


static double parmin[] ={
 /* 001 */   -6.45,
 /* 002 */   -1000,
 /* 003 */    0.,-1.};

static double parmax[] ={
 /* 001 */   6.45,
 /* 002 */   1000,
 /* 003 */   1.0,1.};

struct pparams {
  int lognormal;
  int invfit;
  double *y;
  double *ysq;
  double *xsq;
  double *xerr;
  double *yerr;
  double *yerrsq;
  double *xerrsq;
  double *infoarray;
  double ln10; // Equals ln10 if lognormal, 1 if not.

};

  

#define LN10 2.302585093

int linemodel(double *x, double *pars, double *model, double *errors,
	      double *prior, unsigned long ndata, void *extra)
{

  unsigned long i;
  double diff,intrinsicsq,totalerr,slope,intercept;

  struct pparams *p = (struct pparams *)extra;

  if (p->invfit) {
    slope = 1./pars[0];
    intercept = -pars[1]/pars[0];
  } else {
    slope = pars[0];
    intercept = pars[1];
  }

  //totalerr = pars[2]*sqrt(1+slope*slope);
  //if (totalerr > 1.) totalerr /= pow(totalerr,5.);
  double sectheta = sqrt(1+slope*slope);
  totalerr = pars[2]/p->ln10;
  intrinsicsq  = totalerr*totalerr;
  (*prior) = 0;
  for (i = 0; i < ndata; ++i) {
    // Simple model as above
    model[i] = slope*x[i]+intercept;

    errors[i] = sqrt((p->yerrsq[i]+slope*slope*p->xerrsq[i]+
		      intrinsicsq*(p->lognormal ? 1 : p->ysq[i])
		      -2*pars[3]*slope*p->xerr[i]*p->yerr[i]));

    (*prior) += 2*log(errors[i]/sectheta);

  }
  //p->infoarray[0] = sqrt(intrinsicsq > 0 ? intrinsicsq : 0);
/*   p->infoarray[0] = totalerr*p->ln10; */
/*   p->infoarray[1] = slope; */

  // Jo bovy: add an ndata log(cos(theta))
/*   From the way we derived Eqn. (32) [or (35)], you should not expect to */
/* get the same answer after multiplying the y-axis by three. This is */
/* because the marginalization over the real location of the data point */
/* on the line involves a prior, and the prior we chose is not invariant */
/* under multiplication of the y-axis by a constant. We put a uniform */
/* prior over the location along the line. If instead you would use a */
/* prior that is uniform on x [or y], Eqn. (32) would contain an extra */
/* term of +Ndata*log(|cos(theta)|) [or +Ndata*log(|sin(theta)|), since */
/* "location along the line = x * cos(theta) + y * sin(theta)". If you */
/* add this term you should get results that conform to slope = 3 * m */
/* after multiplying all of your y values by three. */

  return 0; 
}
// M ~ rdelta^3.5; rdelta ~ T^0.33; M~T^(0.3*3.5)

int main(int argc, char *argv[])
{

  struct hrothgar_setup setup;
  struct pparams p;
  double *x,*y,*xe,*ye;
  unsigned long i;

  hrothgar_init_welcomestring(&setup,WELCOMESTRING);
  hrothgar_init_stringpars(&setup,NSTRINGPARAMS,stringparname,stringparval,
			   stringparcomments);

  hrothgar_init_pars(&setup,NTOTPARAMS,
		     paramnames,parmin,parmax,initvalues,
		     dofit,paramunits);
  
  hrothgar_init(&setup,argc,argv);

  char **colslist;
  int ncols = parse_list(stringparval[3],',',&colslist);
  if (ncols != 4)
    BYE("Incorrect number of columns %d",ncols);

  x=y=xe=ye=NULL;
  unsigned long nbins;
  x = double_readdata(stringparval[0],atoi(colslist[0]),&nbins,0);
  xe = double_readdata(stringparval[0],atoi(colslist[1]),&nbins,0);
  y = double_readdata(stringparval[0],atoi(colslist[2]),&nbins,0);
  ye = double_readdata(stringparval[0],atoi(colslist[3]),&nbins,0);
  p.lognormal = atoi(stringparval[1]);
  p.invfit = atoi(stringparval[4]);
  
  
  
  p.xsq = sci_dvector(nbins);
  p.y = sci_dvector(nbins);
  p.xerr = sci_dvector(nbins);
  p.yerr = sci_dvector(nbins);
  p.ysq = sci_dvector(nbins);
  p.xerrsq = sci_dvector(nbins);
  p.yerrsq = sci_dvector(nbins);
  
  double xmin=1.E30,xmax=-1.E30,ymin=1.E30,ymax=-1.E30;
  float *plotx = sci_fvector(nbins);
  float *ploty = sci_fvector(nbins);
  float *plotxe = sci_fvector(nbins);
  float *plotye = sci_fvector(nbins);

  p.ln10 = (p.lognormal ? LN10 : 1);
  for (i = 0; i < nbins; ++i) {
    if (p.lognormal) {
      //x[i] = x[i]*x[i];
      //xe[i] = xe[i]*x[i]*2;
      xe[i] = xe[i]/(LN10*x[i]);
      ye[i] = ye[i]/(LN10*y[i]);
      x[i] = log10(x[i]);
      y[i] = log10(y[i]);
/*       xe[i] = xe[i]/(x[i]); */
/*       ye[i] = ye[i]/(y[i]); */
/*       x[i] = log(x[i]); */
/*       y[i] = log(y[i]); */
    }
    p.xsq[i] = x[i]*x[i];
    p.y[i] = y[i];
    p.ysq[i] = y[i]*y[i];
    p.xerr[i] = xe[i];
    p.yerr[i] = ye[i];
    p.xerrsq[i] = xe[i]*xe[i];
    p.yerrsq[i] = ye[i]*ye[i];
    if (x[i]+1.1*xe[i] > xmax) xmax = x[i]+1.1*xe[i];
    if (y[i]+1.1*ye[i] > ymax) ymax = y[i]+1.1*ye[i];
    if (x[i]-1.1*xe[i] < xmin) xmin = x[i]-1.1*xe[i];
    if (y[i]-1.1*ye[i] < ymin) ymin = y[i]-1.1*ye[i];
    plotx[i] = x[i];
    ploty[i] = y[i];
    plotxe[i] = xe[i];
    plotye[i] = ye[i];
  }
/*   setup.ninfo = 3; */
/*   setup.infoarray = sci_dvector(setup.ninfo); */
/*   setup.infoname = sci_stringmatrix(setup.ninfo,100); */
/*   strcpy(setup.infoname[0],"intrinsic"); */
/*   strcpy(setup.infoname[1],"slope"); */
/*   strcpy(setup.infoname[2],"invb"); */
/*   p.infoarray =setup.infoarray; */
  

  hrothgar(nbins,x,y,ye,linemodel,&p,&setup);
  printf("\n%E %E %E\n",setup.bestfit[0],setup.bestfit[1],setup.bestfit[2]); 

#ifdef HAVE_LIBCPGPLOT
  char *ofname = NULL;
  printf("%s\n",stringparval[2]);
  if (!strstr(stringparval[2],"NONE")) {
    if (!strstr(stringparval[2],"X11")) {
      ofname = sci_strcat(stringparval[2],"/cps");
      cpgopen(ofname);
    } else
      cpgopen("/XWINDOW");
  }  else
    exit(0);

  cpgask(0);
  cpgscf(2);
  cpgenv(xmin,xmax,ymin,ymax,1,30);
  cpgpt(nbins,plotx,ploty,-1);
  cpgmove(xmin,setup.bestfit[0]*xmin+setup.bestfit[1]);
  cpgdraw(xmax,setup.bestfit[0]*xmax+setup.bestfit[1]);
  cpgerrb(5,nbins,plotx,ploty,plotxe,1.);
  cpgerrb(6,nbins,plotx,ploty,plotye,1.);

  cpgsch(1.5);

  if (strstr(stringparval[2],"X11")) 
    sleep(10000);
  
  return 0;
#endif

}

