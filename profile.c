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
    /* 001 */    "project",
    /* 002 */    "model",
    /* 003 */    "write",
    /* 004 */    "conc"
};
  
char *stringparval[NSTRINGPARAMS] = {
    /* 000 */    "profile.txt",
    /* 001 */    "1",
    /* 002 */    "0",
    /* 003 */    "0",
    /* 004 */    "0"
};
  
char *stringparcomments[NSTRINGPARAMS] = {
    /* 000 */    "Data C1=r/R C2=SB C3=SB_Err",
    /* 001 */    "Project profile",
    /* 002 */    "0:NFW 1:Gamma 2:Einasto 3:BPower 4:BPowerSq",
    /* 003 */    "Output results?",
    /* 004 */    "Fit concentrations?"
};

#define WELCOMESTRING "2D/3D profile fitter v" VERSION "\nCopyright (C) 2012 Andisheh Mahdavi \n"

#define NTOTPARAMS 5

static char *paramnames[] = {

 /* 000 */   "norm",    
 /* 001 */   "r0",
 /* 002 */   "slope1",
 /* 003 */   "slope2",
 /* 004 */   "back",
};

static int dofit[] = {
 /* 000 */   0,
 /* 001 */   0,
 /* 002 */   0,
 /* 003 */   0, 0};

static char *paramunits[] = {
 /* 000 */   "Normalization",    
 /* 001 */   "Transition radius",
 /* 002 */   "Inner Slope",
 /* 003 */   "Outer Slope",
 /* 004 */   "Background"};

static double initvalues[] = {
 /* 000 */   1.,    
 /* 001 */   0.1,
 /* 002 */   1.0,
 /* 003 */   3.0, 0.};


static double parmin[] ={
 /* 000 */   1.e-8,    
 /* 001 */   0.0,
 /* 002 */   0.,
 /* 003 */   2.,0.};

static double parmax[] ={
 /* 000 */   1.e8,
 /* 001 */   1000.,
 /* 002 */   2.,
 /* 003 */   6.,1.e3
};

struct profile {
  int proj;
  int write;
  int fitconc;

  double (*model)(double r, void *params);
  double Rsq;
  double r0;
  double is;
  double back;
  double os;
  size_t limit;
  gsl_integration_workspace *workspace;
  double *y;
  double *err;

};

double nfw(double r, void *params)
{
  struct profile *p = (struct profile *)params;
  double x = r/p->r0;

  return pow(x,-p->is)*pow(1.+x,p->is-p->os);
}

double gammprof(double r, void *params)
{
  struct profile *p = (struct profile *)params;
  double x = r/p->r0;
  
  return pow(x,-p->is)*pow(1.+x,p->is-4.);
}

double surfint(double z, void *params)
{
  struct profile *p = (struct profile *)params;

  double r = sqrt(p->Rsq+z*z);
  return p->model(r,p);
}

double projprof(double R, struct profile *p)
{
  gsl_function f;
  double result, abserr;

  p->Rsq = R*R;

  f.function = surfint;
  f.params = p;

  gsl_integration_qagiu(&f,0.,1.e-6,1.e-6,p->limit,p->workspace,
			&result,&abserr);

  return 2.*result;
}

int profile(double *x, double *pars, double *model, double *errors, double *logprior,
	     unsigned long ndata, void *extra)
{

  struct profile *p = (struct profile *)extra;
  unsigned long i;
  double func;

  if (p->fitconc)
    p->r0 = 1./pars[1];
  else
    p->r0 = pars[1];

  p->is = pars[2];
  p->os = pars[3];
  p->back = pars[4];

  for (i = 0; i < ndata; ++i) {
    if (p->proj)
      func = projprof(x[i],p);
    else
      func = p->model(x[i],p);
    model[i] = pars[0]*func+p->back;
  }

  if (p->write) {
    char *ofname = sci_strcat(stringparval[0],".dat");
    FILE *of = fopen(ofname,"w+");

    for (i = 0; i < ndata; ++i) {
      fprintf(of,"%E %E %E %E\n",x[i],p->y[i],p->err[i],model[i]);
    }
    fclose(of);
    free(ofname);
  }

  return 0;
}

int main(int argc, char *argv[])
{

  struct hrothgar_setup setup;
  struct profile p;
  double *x,*y,*err;

  hrothgar_init_welcomestring(&setup,WELCOMESTRING);
  hrothgar_init_stringpars(&setup,NSTRINGPARAMS,stringparname,stringparval,
			   stringparcomments);

  hrothgar_init_pars(&setup,NTOTPARAMS,
		     paramnames,parmin,parmax,initvalues,
		     dofit,paramunits);
  
  hrothgar_init(&setup,argc,argv);

  x=y=err=NULL;
  printf("OK\n");
  unsigned long nbins = read_ascii(stringparval[0],"%f %f %f",NULL,&x,&y,&err);
  p.y = y;
  p.err = err;

  printf("%s %s %f %ld\n",stringparval[1],stringparval[2],x[0],nbins);
  p.proj = atoi(stringparval[1]);
  p.workspace = gsl_integration_workspace_alloc(1000);
  p.limit = 1000;
  p.write = atoi(stringparval[3]);
  p.fitconc = atoi(stringparval[4]);

  switch (atoi(stringparval[2]))
    {
    case 0: p.model=nfw; break;
    case 1: p.model=gammprof; break;
    default: exit(-1); break;
    }

  hrothgar(nbins,x,y,err,profile,&p,&setup);

  return 0;

}

