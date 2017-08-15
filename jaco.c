#include <jaco.h>
#include <stdlib.h>
#include <time.h>
#include <string.h>
#include <math.h>
#include <time.h>
#include <sys/types.h>
#include <sys/stat.h>
#include <fcntl.h>
#include <unistd.h>
#include <gsl/gsl_sf.h>
#include <gsl/gsl_math.h>
#include <gsl/gsl_errno.h>
#include <gsl/gsl_sort.h>
#include <hrothgar/hrothgar_proto.h>

// Constants for Gauss-Legendre integration
#define NGAUSS 10
const double xgauss[] = 
  {    0.013046736,    0.067468317,    0.160295216,    0.283302303,
       0.425562831,    0.574437169,    0.716697697,    0.839704784,
       0.932531683,    0.986953264};
const double wgauss[] = 
  {    0.033335672,    0.074725675,    0.109543181,    0.134633360,
       0.147762112,    0.147762112,    0.134633360,    0.109543181,
       0.074725675,    0.033335672};
const double xgausssq[] = 
  {    0.000170217,    0.004551974,    0.025694556,    0.080260195,
       0.181103723,    0.329978062,    0.513655589,    0.705104125,
       0.869615340,    0.974076746};
const double wgausssq[] = 
  {    0.000869843,    0.010083231,    0.035118496,    0.076283882,
       0.125764126,    0.169760099,    0.192982838,    0.183967867,
       0.139368118,    0.065801501};
const double wgaussinf[] = 
  {    5.110193510,    2.215133807,    1.366767943,    0.950457220,
       0.694431476,    0.514458744,    0.375704736,    0.260908794,
       0.160264098,    0.067552686};

// Hernquist-like gamma profile: r^-n (1+r)^(n-4)
 double dens_hernquistlike(double r, struct jaco_state *js)
{
  double x;

  JS(rdm1); JS(ndark); JS(nm3);

  x = r/rdm1;
  
  return -nm3*rdm1*exp(-ndark*log(r)+(ndark-4.)*log(r+rdm1))/FOURPI;
}

// Corresponding mass profile
 double mass_hernquistlike(double r, struct jaco_state *js)
{
  JS(rdm1); JS(nm3);

  return exp(nm3*log(1.+rdm1/r));
}

 double dens_powlawwithcore(double r, struct jaco_state *js)
{

  JS(rdm1); JS(ndark);

  return pow(r+rdm1,-ndark)/(2.*FOURPI);
}

 double mass_powlawwithcore(double r, struct jaco_state *js)
{
  double rrdm1, answer;

  JS(rdm1); JS(rdm1sq); JS(rdm3mn); JS(nm1); JS(nm2); JS(rdiff);

  rrdm1 = r+rdm1;

  if (fabs(nm2) > 1.E-4) {
    answer= (rdm3mn-(rdm1sq+r*nm1*(rdm1+r*nm2/2.))/pow(rrdm1,nm1))/rdiff;
  } else if (rdm1 < 1.E-5) {
    answer= r/2.;
  } else  {
    answer= r*(1.+rdm1/rrdm1)/2.+rdm1*log(rdm1/rrdm1);
  }
  return answer;
}

 double dens_nfwwithcore(double r, struct jaco_state *js)
{
  double rrdm1;
  
  JS(rdm0); JS(rdm1);

  rrdm1 = r+rdm1;
  return 1./(FOURPI*(rdm0+r)*rrdm1*rrdm1);
}

 double mass_nfwwithcore(double r,struct jaco_state *js)
{
  double sum,answer;
  
  JS(rdm1); JS(rdm0); JS(rdiff); JS(rdm0sq);
  
  sum = r+rdm1;

  if (rdm0 > 1e-6)
    answer = (-rdm1*r/sum+(-rdm0sq*log(rdm0/(rdm0+r))+
			   rdm1*(rdiff-rdm0)*log(sum/rdm1))/rdiff)/rdiff;
  else
    answer = log(sum/rdm1)-r/sum;

  return answer;
}


 double dens_isowithcore(double r, struct jaco_state *js)
{
  JS(rdm0); JS(rdm1);

  return 1./(FOURPI*(rdm0+r)*(r+rdm1));
}

 double mass_isowithcore(double r, struct jaco_state *js)
{
  double answer;

  JS(rdm0); JS(rdm1); JS(rdiff); JS(rdm0sq); JS(rdm1sq);

  if (rdm0 < 1.E-5)
    answer = r - rdm1*(log(r/rdm1+1.));
  else
    answer = r-(rdm0sq*log((r/rdm0+1.))-
		rdm1sq*log(r/rdm1+1.))/(rdiff);

  return answer;
}

 double dens_sersic_dark(double r,struct jaco_state *js)
{
  JS(ndark); JS(rdm1);

  return ndark*exp(-pow(r/rdm1,ndark))/
    (FOURPI*rdm1*rdm1*rdm1*gsl_sf_gamma(3./ndark));
}

 double mass_sersic_dark(double r,struct jaco_state *js)
{
  double x,gamm_a;
  
  JS(ndark); JS(rdm1);

  gamm_a = 3./ndark;
  x = pow(r/rdm1,ndark);

  //if (x > gamm_a+20) return 1.;

  return gsl_sf_gamma_inc_P(gamm_a,x);
}

 double dens_sersic_stellar(double r,struct jaco_state *js)
{
  JS(nstar); JS(rstar);

  return nstar*exp(-pow(r/rstar,nstar))/
    (FOURPI*rstar*rstar*rstar*gsl_sf_gamma(3./nstar));
}

 double mass_sersic_stellar(double r,struct jaco_state *js)
{
  double x,gamm_a,answer;
  
  JS(nstar); JS(rstar);

  gamm_a = 3./nstar;
  x = pow(r/rstar,nstar);
  answer = gsl_sf_gamma_inc_P(gamm_a,x);
  //printf("Sersic %E %E %E %E\n",r,nstar,rstar,answer);

  //if (x > gamm_a+20) return 1.;
  return answer;
}

 double dens_universal(double r,struct jaco_state *js)
{
  double rrdm1;

  JS(rdm1); JS(ndark);

  rrdm1 = r+rdm1;
  
  return pow(r/rrdm1,-ndark)/(rrdm1*rrdm1*rrdm1*FOURPI);
}

 double mass_universal(double r,struct jaco_state *js)
{
  double x;

  JS(rdm1); JS(ndark);

  x = r/(r+rdm1);

  return genbetai(3.-ndark,0.,x);
}

#define BETAMODEL(xsq,b) pow(1.+xsq,-3.*b/2)
#define BETAMODEL_A(xsq,a,b) (pow(xsq,-a/2)*pow(1.+xsq,(a-3.*b)/2))

#define R_DBETAMODEL(xsq,a,b) ((3.*b*xsq+a)/(1.+xsq))

double gasdensity(double r,struct jaco_state *js)
{
  double x1,x2,x3,x4,x1sq,x2sq,x3sq,x4sq;
  
  JS(rx1); JS(rx2); JS(rx3);  
  JS(b1); JS(b2); JS(b3); JS(b4);
  JS(alpha); 
  JS(rhog1); JS(rhog2); JS(rhog3); 

  x1 = r/rx1; x1sq = x1*x1;
  x2 = r/rx2; x2sq = x2*x2;
  x3 = r/rx3; x3sq = x3*x3;

  return 
    rhog1*BETAMODEL_A(x1sq,alpha,b1)+
    rhog2*BETAMODEL(x2sq,b2)+
    rhog3*BETAMODEL(x3sq,b3);

}

/* double gasdensity_model2_slope(double r,struct jaco_state *js) */
/* { */
/*   double x1,x2,x3,x4,x1sq,x2sq,x3sq,x4sq,bm1,bm2,bm3,bm4; */
  
/*   JS(rx1); JS(rx2); JS(rx3); JS(rx4);  */
/*   JS(b1); JS(b2); JS(b3); JS(b4); */
/*   JS(alpha); */
/*   JS(rhog1); JS(rhog2); JS(rhog3); JS(rhog4); */

/*   x1 = r/rx1; x1sq = x1*x1;  */
/*   x2 = r/rx2; x2sq = x2*x2; */
/*   x3 = r/rx3; x3sq = x3*x3; */

/*   if (alpha > 0.)  */
/*     bm1 = rhog1*BETAMODEL_A(x1sq,alpha,b1); */
/*   else */
/*     bm1 = rhog1*BETAMODEL(x1sq,b1); */
  
/*   bm2 = rhog2*BETAMODEL(x2sq,b2); */
/*   bm3 = rhog3*BETAMODEL(x3sq,b3); */

/*   return */
/*     (bm1*R_DBETAMODEL(x1sq,alpha,b1)+ */
/*      bm2*R_DBETAMODEL(x2sq,alpha,b2)+ */
/*      bm3*R_DBETAMODEL(x3sq,alpha,b3))/(bm1+bm2+bm3); */

/* } */


 double gasmass_beta(double rr, double b_alpha,
			   double rx, double rxsq, double b_beta)
{
  double r,x,rsq;

  //if (b_alpha > 0. && rr > b_alpha) r = b_alpha; else 
  r = rr;

  rsq = r*r;

  x = 1./(1+rxsq/rsq);

  if (fabs(b_beta-1.) < 1.E-4)
    return 2.*(asinh(r/rx)-sqrt(x));
  else 
    return genbetai((3.-b_alpha)/2.,1.5*(b_beta-1.),x);
}

 double gasmass(double r,struct jaco_state *js)
{
  double rsq;
  rsq = r*r;

  JS(rx1); JS(rx2); JS(rx3); JS(b1); JS(b2); JS(b3);  JS(b4);
  JS(rx1sq);  JS(rx2sq);  JS(rx3sq); 
  JS(Mg1); JS(Mg2); JS(Mg3); 


  return Mg1*gasmass_beta(r,js->alpha,rx1,rx1sq,b1)+
    Mg2*gasmass_beta(r,0,rx2,rx2sq,b2)+
    Mg3*gasmass_beta(r,0,rx3,rx3sq,b3);
  
}

double light(double r,struct jaco_state *js)
{
  double starterm=0.;

  JS(MsContr); JS(Ms);

  if (MsContr > 0.) 
    starterm = Ms*js->stellarmass(r,js);

  return gasmass(r,js)+starterm;
    
}

double allmass(double r,struct jaco_state *js) { 

  double starterm=0.,gasterm=0.,darkterm,answer; 
  
  JS(MsContr); JS(Ms); JS(Mg1Contr);
  JS(Md);  
  JS(darkmodel); 


  if (!js->totalmass) {
    
    if (MsContr > 0.)
      starterm = Ms*js->stellarmass(r,js);

    if (r < js->rshock && js->rshock > 0 && Mg1Contr > 0)
      gasterm = gasmass(r,js);
  } 
  

  if (abs(darkmodel) == MASS_FOLLOWSLIGHT)
    darkterm = Md*(gasterm+starterm);
  else
    darkterm = Md*js->darkmass(r,js);

  if (js->totalmass) 
    answer = darkterm;
  else
    answer = starterm+gasterm+darkterm;

  //printf("*****************%f %f %f %f %f\n",r,Ms,starterm,gasterm,darkterm);
  
  return answer; 

 } 

 double temperature_intfunc(double r,void *par)
{

  double darkterm,gasterm,starterm=0,answer;
  struct jaco_vars *jv = (struct jaco_vars *)par;
  struct jaco_state *js = jv->state;
  
  JS(MsContr); JS(Ms); 
  JS(Md); 

  darkterm = Md*js->darkmass(r,js);
  gasterm = starterm = 0;

  if (!js->totalmass) {
    if (MsContr > 0.) 
      starterm = Ms*js->stellarmass(r,js);

    gasterm = gasmass(r,js);

  }

  answer = (darkterm+gasterm+starterm)*gasdensity(r,js)/(r*r);

  //printf("%10.3E %10.3E %10.3E %10.3E %10.3E\n",r,starterm,gasterm,darkterm,answer);

  return answer;

}

 double massdensity(double r,struct jaco_state *js)
{
  double starterm=0.,gasterm=0.,darkterm=0.,answer;

  JS(MsContr); JS(Ms); 
  JS(darkmodel); JS(Md); JS(projtype);


  if (!js->totalmass) {
    if (MsContr > 0. && (projtype == PROJ_STARS || projtype == PROJ_ALL))
      starterm = Ms*js->stellardensity(r,js);

    if (r < js->rshock && js->rshock > 0)
      if (projtype == PROJ_GAS || projtype == PROJ_ALL)
	gasterm = gasdensity(r,js);
  } 
  

  if (projtype == PROJ_DARK || projtype == PROJ_ALL) {
    if (abs(darkmodel) == MASS_FOLLOWSLIGHT)
      darkterm = Md*(gasterm+starterm);
    else
      darkterm = Md*js->darkdens(r,js);
  }

  if (js->totalmass) 
    answer = darkterm;
  else
    answer = starterm+gasterm+darkterm;

  //printf("*****************%f %f %f %f %f\n",r,Ms,starterm,gasterm,darkterm);
  
  return answer;  

}

void abscissae(double R, double Rb, struct jaco_vars *jv)
{
  double Rn;
  int   i;

  JVS(rshock); 

  if (rshock > 0.) {
    if (R <= Rb) {
      Rn = Rb; 
      jv->nabs = 2*NGAUSS;
    } else {
      Rn = R;
      jv->nabs = NGAUSS;
    }

    for (i = 0; i < NGAUSS; ++i) 
      jv->xabs[i] = (rshock-Rn)*xgausssq[i]+Rn;
    
    if (R <= Rb) 
      for (i = 0; i < NGAUSS; ++i) 
	jv->xabs[i+NGAUSS] = (Rb-R)*xgausssq[i]+R;

  } else {
    if (R <= Rb) {
      Rn = Rb; 
      jv->nabs = 2*NGAUSS;
    } else {
      Rn = R;
      jv->nabs = NGAUSS;
    }
    
    for (i = 0; i < NGAUSS; ++i) 
      jv->xabs[i] = Rn/xgausssq[i];

    if (R <= Rb) 
      for (i = 0; i < NGAUSS; ++i) 
	jv->xabs[i+NGAUSS] = (Rb-R)*xgausssq[i]+R;
    
  }

    
}


// Carries out the integral of the given function from R to infinity,
// doubling the resolution around R = Rb, if necessary.
 double ninteg_nr(double (*func)(double r, void *par),double R,double Rb,
			struct jaco_vars *jv)
{
  double x,Rn,sum=0.,suma=0.,diff;
  int   i;

  if (R <= Rb) Rn = Rb; else Rn = R;

  for (i = 0; i < NGAUSS; ++i) 
    sum += wgaussinf[i]*func(Rn/xgausssq[i],jv);
  sum *= Rn;
  
  if (R <= Rb) {
    diff = Rb-R;
    for (i = 0; i < NGAUSS; ++i) {
      x = diff*xgausssq[i]+R;
      suma += wgausssq[i]*func(x,jv);
    }
    suma *= (Rb-R);
  } 

  return sum+suma;
    
}

// Carries out the integral of the given function from R to rshock,
// doubling the resolution around R = Rb, if necessary.
 double ninteg_shock(double (*func)(double r, void *par),
			   double R,double Rbreak, 
			   struct jaco_vars *jv)
{
  double x,Rn,Rb,Rinf,sum=0.,suma=0.,diff;
  int   i;

  JVS(rshock);

  Rinf = rshock;
  Rb = Rbreak;
  if (Rbreak < 0) { Rinf = Rb = -Rbreak; }

  if (R <= Rbreak) Rn = Rb; else Rn = R;

  diff = Rinf-Rn;
  for (i = 0; i < NGAUSS; ++i) {
    x = diff*xgausssq[i]+Rn;
    sum += wgausssq[i]*func(x,jv);
  }
  sum *= diff;
  
  if (R <= Rbreak) {
    diff = Rb-R;
    for (i = 0; i < NGAUSS; ++i) {
      x = diff*xgausssq[i]+R;
      suma += wgausssq[i]*func(x,jv);
    }
    suma *= diff;
  } 

  return sum+suma;
    
}

double ninteg_gsl(double (*func)(double r,void *par),double R1, double R2, 
		  struct jaco_vars *jv)
{
  int status,inttype;
  gsl_function func_gsl;
  double answer1=0., answer2;
  double err;

  func_gsl.function = func;
  func_gsl.params = jv;


  if (gsl_fcmp(R1,R2,1.e-4) == 0) return 0;
  status = -1;
  inttype = GSL_INTEG_GAUSS15;
  double precision = jv->state->precision;
  
  // Projected temperatures will have inherently less precision than 
  // the spectra themselves
  if (jv->calctproj) precision *= 100;

  if (R2 < R1)
    BYE("R2 < R1 in ninteg_gsl %E < %E",R2,R1);
  
  
  gsl_set_error_handler_off();
  status = gsl_integration_qag(&func_gsl,R1,R2,0., 
			       precision, 
			       INTMAX,
			       inttype, 
			       jv->workspace[0], 
			       &answer2,&err); 
  
  //if (status != GSL_SUCCESS && jv->state->tempwarn == 0) {
  //printf("JACO: could not achieve desired accuracy from %E to %E\n",R1,R2);
  // printf("Integrand  is %s\n",jv->integstate);
  //printf("GSL description: %s\n",gsl_strerror(status));
    //}

  return answer1+answer2;
}

// Metallicity profile
double metallicity(double r, struct jaco_vars *jv)
{
  double x,Z;
  
  JVS(z0); JVS(zinf); JVS(rz);

  x = r/rz;
  //y = r/rContr;
  //Z = (z0*(1.-y)+zinf*y*(1.+rContr/rz))/(1.+x);
  //if (Z < 0) Z = 0.001;
  //zac = ((rshock+20.*rz)*z0-(rshock+rz)*zinf)/(19.*rz);
  //zbc = (20.*(rshock+rz)*zinf-(rshock+20*rz)*z0)/(19.*rshock);
  //Z = (zac+zbc*x)/(1.+x);
  if (r < 0)
    Z = -r;
  else
    //Z = (z0+x*zinf)/(1.+x);
    Z = z0*pow(1+x*x,-3.*zinf/2.);

  // Determine mean molecular weight, hydrogen mass fraction,
  // And ion to hydrogen ratio for this abundance.

  // Using Grevesse & Sauval 1998 abundances
  jv->mu = 0.598782+0.006825*Z;
  jv->Hfrac = 0.740442-0.012406*Z;
  jv->nenh = 1.170201+0.011532*Z;

  return Z;
}

// Interpolates the rebinned cooling matrix to obtain the appropriate
// cooling function.
double cooling(double lT, double lZ, struct jaco_vars *jv)
{
  double  tt, zz, logT, logZ, answer, extrapolate_factor=0.;
  unsigned long tpos,mpos;

  JVS(tgridsize); JVS(mgridsize); 
  JVS(tempaxis); JVS(metalaxis); JVS(debug); JV(coolmatrix);

  logT = lT;
  logZ = lZ;
  tpos = gsl_interp_bsearch(tempaxis,logT,0,tgridsize-1);
  mpos = gsl_interp_bsearch(metalaxis,logZ,0,mgridsize-1);
  jv->bolocorr = jv->state->bolometric[tpos][mpos];

  // Temperature request is smaller than that tabulated. We assume that
  // at such temperatures, virtually no emission will fall into the 
  // satellite's sensitivity range.
  if (logT < tempaxis[0]) return 0.;

  // Metallicity requested is smaller than tabulated. We assume that
  // the smallest available metallicity is good enough (i.e., very
  // close to 0.
  if (logZ < metalaxis[0]) { logZ = metalaxis[0]; mpos = 0; }

  // Temperature or metallicity exceed expectations. This is most
  // likely the result of the fitting program overstepping
  // its bounds. Return zero emission.

  if (logT > tempaxis[tgridsize-1]) {
    if (!jv->state->tempwarn && debug > 0) {
      printf("Warning: temperature limit exceeded(%f>%f); ",
	     logT,tempaxis[tpos]);
      jv->state->tempwarn=1;

    }
    --tpos;
    extrapolate_factor = 0.5*(logT-tempaxis[tpos]);
    logT = tempaxis[tpos];
  }
  if (logZ > metalaxis[mgridsize-1]) {
    if (!jv->state->metalwarn && debug > 0) {
      printf("Warning: upper metallicity limit exceeded. %f > %f; ",
	     logZ,metalaxis[mpos]);
      jv->state->metalwarn=1;
    }
    --mpos;
    logZ = metalaxis[mpos];
  }


  // This is the bilinear interpolation forumla from Numerical
  // Recipies. We do it in log-log-log space.
  tt = (logT - tempaxis[tpos])/(tempaxis[tpos+1]-tempaxis[tpos]);
  zz = (logZ - metalaxis[mpos])/(metalaxis[mpos+1]-metalaxis[mpos]);


  answer = ((1.-tt)*(1.-zz)*coolmatrix[tpos][mpos]+
	    tt*(1.-zz)*coolmatrix[tpos+1][mpos]+
	    tt*zz*coolmatrix[tpos+1][mpos+1]+
	    (1.-tt)*zz*coolmatrix[tpos][mpos+1]);

  //printf("***%E %E %E %E %E***\n",coolmatrix[tpos][mpos],
  // coolmatrix[tpos+1][mpos],coolmatrix[tpos+1][mpos+1],
  // coolmatrix[tpos][mpos+1],answer);

  if (answer > 10) 
    printf("Warning---large cooling function! %f %f %f\n",logT,logZ,answer);

  return pow(10.,answer+extrapolate_factor);
}

// Returns the function to be integrated in order to obtain the
// integrated emissivity between projected radii R1 and R2. This is
//
//   rho_g^2 cooling(T(r),metallicity(r)) k(r)
//
//   where k(r) = sqrt(r^2-R1^2) if r < R2
//              = sqrt(r^2-R1^2)-sqrt(r^2-R2^2) if r > R2
// 
// This function is integrate from R1 to infinity, with a mandatory
// split of the integral at R2 (where there is a sharp discontinuity)

double emintegrand(double r, void *par)
{
  double rsq,kern,answer;
  struct jaco_vars *jv = (struct jaco_vars *)par;

  JV(R1sq); JV(R2sq); JVS(interp_coefs); JVS(rspline); JVS(emissivity);
  JVS(wgt_interp_coefs);  JVS(wgt_emissivity); JVS(unwgt_emissivity);
  JVS(nsplines);

  rsq = r*r;
  
  kern = sqrt(rsq-R1sq);
  if (rsq > R2sq) kern -= sqrt(rsq-R2sq);

  kern *= r;

  if (r < rspline[0] || r > rspline[nsplines-1]) {
    printf("Insufficient convergence of inner gas density profile.\n");
    printf("The integrator is requesting access to the emissivity around\n");
    printf("%f Mpc, which indicates a non-integrable singularity.\n",r);
    exit(-1);
  }

  if (jv->calctproj == 1)
    answer = kern*
      gsl_interp_eval(wgt_interp_coefs[jv->specbin],
		      rspline,wgt_emissivity[jv->specbin],
		      r,jv->accel);
  else if (jv->calctproj == 2)
    answer = kern*
      gsl_interp_eval(wgt_interp_coefs[jv->specbin],
		      rspline,unwgt_emissivity[jv->specbin],
		      r,jv->accel);
  else
    answer = kern*
      gsl_interp_eval(interp_coefs[jv->specbin],
		      rspline,emissivity[jv->specbin],
		      r,jv->accel);
  
  return answer;

}


// This is an integrand which returns interpolated function along the z
// direction for an x and y coordinate stored in jv->x,y. Triaxial
// Euler angles are taken into account via x_x,x_y, etc. The function has
// triaxial radius abscissae rspline, interpolating coefficients
// integrate_coefs, and ordinates integrate_ordinates.
double integrand_xy(double z, void *par)
{

  double answer; 

  struct jaco_vars *jv = (struct jaco_vars *)par;
  struct jaco_state *js = jv->state;

  JV(x); JS(x_x); JS(x_y); JS(x_z);
  JV(y); JS(y_x); JS(y_y); JS(y_z);
  JS(z_x); JS(z_y); JS(z_z);
  double xp = x_x*x+x_y*y+x_z*z;
  double yp = y_x*x+y_y*y+y_z*z;
  double zp = z_x*x+z_y*y+z_z*z;

  double rt = sqrt(xp*xp+yp*yp/js->asq+zp*zp/js->bsq);
  //printf("%E %E %E \n",x_x,y_y,z_z);

  // No gas outside of termination shock
  if (rt > js->rshock)
    return 0.;
  //printf("%E %E %E    %E %E %E    %E %E %E \n",x,y,z,xp,yp,zp,rt,js->rshock,js->rspline[js->nsplines-1]);

  // rspline[0] should be sufficiently close to the center and
  // encompass sufficiently little volume such that the following
  // line is a good approximation
  if (rt < js->rspline[0])
    rt = js->rspline[0];

  int status = gsl_interp_eval_e(jv->integrand_coefs,js->rspline,
				 jv->integrand_ordinate,rt,jv->accel, &answer);

  if (status != GSL_SUCCESS) {
    BYE("Error interpolating profile at triaxial r=%E (error %s)",rt,
	gsl_strerror(status))
  }
  
  if (jv->dolog)
    return pow(10.,answer);
  else
    return answer;

}

void emissivity_profile(struct jaco_state *js)
{
  long i,j;
  int reallocate;
  double r,tfullint,rb,logZ,logT;
  double *pinteg,*normdensity;
  double Hfrac_nenh,Hfracsq_nenh;
  struct jaco_vars jv;

  JS(rshock); JS(rspline); JS(rebinnedcooling); 
  JS(nsplines); JS(nlastbin); JS(pshock);
  JS(binningchanged); JS(tprojpower); JS(samepars);

  if (!binningchanged && samepars)
    return;

  reallocate = binningchanged;

  if (reallocate && js->interp_coefs != NULL) {
    for (i = 0; i < nlastbin; ++i) {
      gsl_interp_free(js->interp_coefs[i]);
      gsl_interp_free(js->wgt_interp_coefs[i]);
    }
    gsl_interp_free(js->szinterp_coefs);
    free(js->interp_coefs);
    free(js->wgt_interp_coefs);
    js->interp_coefs = NULL;
    if (js->calcxray) {
      sci_free_dmatrix(js->emissivity,nlastbin);
      sci_free_dmatrix(js->wgt_emissivity,nlastbin);
    }
  }
  
  if (!js->interp_coefs) {
    js->interp_coefs = (gsl_interp **)
      malloc(nlastbin*sizeof(gsl_interp *));
    js->wgt_interp_coefs = (gsl_interp **)
      malloc(nlastbin*sizeof(gsl_interp *));
    for (i = 0; i < nlastbin; ++i) {
      js->interp_coefs[i] = gsl_interp_alloc(gsl_interp_cspline,nsplines);
      js->wgt_interp_coefs[i] = gsl_interp_alloc(gsl_interp_cspline,nsplines);
    }
    js->szinterp_coefs = gsl_interp_alloc(gsl_interp_cspline,nsplines);
    if (js->calcxray) {
      js->emissivity = sci_dmatrix(nlastbin,nsplines);
      js->wgt_emissivity = sci_dmatrix(nlastbin,nsplines);
      js->unwgt_emissivity = sci_dmatrix(nlastbin,nsplines);
    }
    js->pressure = sci_dvector(nsplines);
  }
  
  pinteg = sci_dvector(nsplines);
  normdensity = sci_dvector(nsplines);

  
  gsl_set_error_handler_off();
#pragma omp parallel private(i,r,rb,jv) default(shared)
  {
    jaco_thread_init(&jv);
    jv.state = js;
    strcpy(jv.integstate,"Emissivity profile");

    jv.workspace[0] = gsl_integration_workspace_alloc(INTMAX);
    jv.workspace[1] = gsl_integration_workspace_alloc(INTMAX);

    #pragma omp for 
    for (i = nsplines-1; i >= 0; --i) {
      r = rspline[i];
      rb = (i == nsplines-1 ? rshock : rspline[i+1]);
      normdensity[i] = gasdensity(r,js);

      // Refer to the temperature T(r) formula above. TMASS gives us
      // G in the corect coordinates; rx2 is required because the
      // other factors are dimensionless.
    
      pinteg[i] = js->ninteg(temperature_intfunc,r,rb,&jv);
    }

    gsl_integration_workspace_free(jv.workspace[0]);
    gsl_integration_workspace_free(jv.workspace[1]);
    jaco_thread_shutdown(&jv);
  }
  gsl_set_error_handler(NULL);

  for (i = nsplines-2; i >= 0; --i) 
    pinteg[i] += pinteg[i+1];

  jv.state = js;

  //FILE *fdf = fopen("fdf3","w+");
  for (i = nsplines-1; i >= 0; --i) {
    js->pressure[i] = TMASS*(pshock+pinteg[i]);

    r = rspline[i];
    double ntpressure = (1.-jv.state->nonthermal*r);
    if (ntpressure < 0.2) ntpressure=0.2;
    js->pressure[i] *= ntpressure;

    logZ = log10(metallicity(r,&jv));
    tfullint = jv.mu*js->pressure[i]/normdensity[i];
    logT = log10(tfullint);
    Hfrac_nenh = jv.Hfrac*jv.nenh;
    Hfracsq_nenh = jv.Hfrac*Hfrac_nenh;
    
    if (js->calcxray) {
      for (j = 0; j < nlastbin; ++j) {
	jv.coolmatrix = rebinnedcooling[j];
	
	js->emissivity[j][i] = 
	  normdensity[i]*normdensity[i]*cooling(logT,logZ,&jv);
	js->emissivity[j][i] *= Hfracsq_nenh;
	/* if (j == 100) */
	/*   fprintf(fdf,"%E %E %E %E %E %E\n",rspline[i],logT,logZ,js->pressure[i], */
	/* 	  normdensity[i],js->emissivity[j][i]); */

	if (js->fitwrite > 2) {
	  js->wgt_emissivity[j][i] = js->emissivity[j][i]*
	    pow(tfullint,tprojpower);
	  js->unwgt_emissivity[j][i] = js->emissivity[j][i]*
	    pow(tfullint,tprojpower-1.);
	}
      }
    }
    js->pressure[i] = log10(Hfrac_nenh*js->pressure[i]);
  }
  //fclose(fdf);

  if (js->calcxray)
    for (j = 0; j < nlastbin; ++j) {
      gsl_interp_init(js->interp_coefs[j],
		      rspline,js->emissivity[j],nsplines);
      gsl_interp_init(js->wgt_interp_coefs[j],
		      rspline,js->emissivity[j],nsplines);
    }
  
  gsl_interp_init(js->szinterp_coefs,rspline,js->pressure,nsplines);

  free(normdensity);
  free(pinteg);

}

void emintegrand_fastfit(struct jaco_vars *jv)
{
  int i;
  double normdensity,rsq,pinteg=0.,rb;


  // Refer to the temperature T(r) formula above. TMASS gives us
  // G in the corect coordinates; rx2 is required because the
  // other factors are dimensionless.

  JVS(rshock); JV(nabs); JVS(calcsz); JV(xabs);
  JVS(pshock); JV(R1sq); JV(R2sq); JVS(tprojpower);
  double ntpressure;
  
  rb = rshock;
  for (i = nabs-1; i >= 0; --i) {
    if (xabs[i] > rb) { rb = rshock; pinteg = 0.; }
    normdensity = gasdensity(xabs[i],jv->state);
    pinteg += 
      jv->state->ninteg(temperature_intfunc,xabs[i],-rb,jv);
    rb = xabs[i];
    // Linear Nonthermal pressure support
    jv->szyinteg[i] = TMASS*(pshock+pinteg);
    ntpressure = (1.-jv->state->nonthermal*xabs[i]);
    if (ntpressure < 0.2) ntpressure=0.2;

    jv->szyinteg[i] *= ntpressure;
    //printf("%E %E %E\n",xabs[i],pshock,pinteg);
    jv->logZ[i] = log10(metallicity(xabs[i],jv));
    jv->tint[i] = jv->mu*jv->szyinteg[i]/normdensity;
    jv->logT[i] = log10(jv->tint[i]);
    rsq = xabs[i]*xabs[i];

    //jv->entropy[i] = jv->tint[i]/cbrt(normdensity*normdensity);
    jv->kernel[i] = sqrt(rsq-R1sq);
    if (rsq > R2sq) jv->kernel[i] -= sqrt(rsq-R2sq);
    jv->kernel[i] *= jv->Hfrac*jv->nenh*xabs[i];

    if (calcsz)
      jv->szyinteg[i] *= jv->kernel[i];

    if (jv->calctproj == 1) jv->kernel[i] *= pow(jv->tint[i],tprojpower);
    if (jv->calctproj == 2) jv->kernel[i] *= pow(jv->tint[i],tprojpower-1.);

    jv->kernel[i] *= jv->Hfrac*normdensity*normdensity;

    // Introduce penalty for declining entropy profile


  }

/*   for (i = nabs-1; i > 0; --i) { */
/*     rb = (xabs[i] - xabs[i-1])*(jv->entropy[i] - jv->entropy[i-1]); */
/*     if (rb < 0.) { */
/*       if (jv->state->debug > 5 || !jv->state->entropywarning) { */
/* 	printf("Warning---entropy bounds exceeded: (%E,%E) to (%E,%E)\n", */
/* 	       xabs[i],jv->entropy[i],xabs[i-1],jv->entropy[i-1]); */
/* 	jv->state->entropywarning=1; */
/*       } */
      
/*       jv->kernel[i] *= 1.e7*rb; */
/*     } */
/*   } */

}

 double emintegrand_fastev(double r, void *par)
{
  int i;
  double answer,cf;
  struct jaco_vars *jv = (struct jaco_vars *)par;

  JV(nabs);

  i = jv->integcounter;

  cf = cooling(jv->logT[i],jv->logZ[i],jv);
  answer = jv->kernel[i]*cf;

  ++jv->integcounter;
  if (jv->integcounter == nabs) jv->integcounter = 0;

  return answer;
}

 double szintegrand(double r, void *par)
{
  double rsq,kern,answer;
  struct jaco_vars *jv = (struct jaco_vars *)par;

  JVS(logfile);

  rsq = r*r;

  JV(R1sq); JV(R2sq); JVS(szinterp_coefs); JVS(rspline); JVS(pressure);
  JVS(debug);

  kern = sqrt(rsq-R1sq);
  if (rsq > R2sq) kern -= sqrt(rsq-R2sq);

  kern *= r;

  answer = kern*
    pow(10.,gsl_interp_eval(szinterp_coefs,rspline,pressure,r,
    			    jv->szaccel));

  if (debug > 1)
    fprintf(logfile,"SZ %E %E %E %E\n",sqrt(R1sq),sqrt(R2sq),r,answer);
  return answer;
}

 double szintegrand_fastev(double r, void *par)
{
  int i;
  double answer;
  struct jaco_vars *jv = (struct jaco_vars *)par;

  JVS(debug); JV(R1sq); JV(R2sq); JV(xabs); JV(szyinteg);
  JVS(logfile);

  i = jv->integcounter;

  answer = szyinteg[i];
  if (debug > 1)
    fprintf(logfile,"SZ %E %E %E %E\n",sqrt(R1sq),sqrt(R2sq),xabs[i],answer);

  ++jv->integcounter;
  if (jv->integcounter == jv->nabs) jv->integcounter = 0;

  return answer;
}

#define read_e(x,y,z) if(!read(x,y,z)) { printf("read error.\n"); exit(1); }

void read_cooling_function(struct jaco_state *js)
{
  long i,j,k;

  int coolfile;
  
  //printf("Taking H_0 = %d km/s/Mpc, Flat Universe, Omega_matter = %.2f\n",
  // HUBBLE,OMEGA_M);

  //printf("Initializing Mahdavi thin plasma code...\n");

  if (js->debug > 1)
    printf("Reading cooling function from %s\n",js->spectralcode);

  // Now read in the description of the cooling function.
  coolfile = open(js->spectralcode,O_RDONLY);
  if (coolfile < 0)
    BYE("Could not open cooling file %s.",js->spectralcode);

  // The first row contains the grid size in energy, temperature, and
  // metallicity. 
  read_e(coolfile,&js->egridsize,sizeof(js->egridsize));
  read_e(coolfile,&js->tgridsize,sizeof(js->tgridsize));
  read_e(coolfile,&js->mgridsize,sizeof(js->mgridsize));
  //printf("%d %d %d\n",js->egridsize,js->tgridsize,js->mgridsize);

  // Allocate memory to contain the whole cooling function
  js->mastercooling = sci_ftensor(js->tgridsize,js->mgridsize,js->egridsize);
  js->bolometric = sci_dmatrix(js->tgridsize,js->mgridsize);
  js->energyaxis = sci_dvector(js->egridsize);
  js->tempaxis   = sci_dvector(js->tgridsize);
  js->metalaxis  = sci_dvector(js->mgridsize);

  // Next we read in the centers of the bins over which the cooling
  // function was integrated. This appears only once.
  read_e(coolfile,js->energyaxis,js->egridsize*sizeof(double));
  read_e(coolfile,js->tempaxis,js->tgridsize*sizeof(double));
  read_e(coolfile,js->metalaxis,js->mgridsize*sizeof(double));

  js->ebinsize = js->energyaxis[2]-js->energyaxis[1];

  // Now we read the temperature and metallicity abscissae, followed
  // by the cooling function itself. Note that the former repeat
  // themselves redundantly throughout the coolfile; preceding each
  // array (of length egridsize) containing the cooling function,
  // there are two numbers specifying the temperature and metallicity
  // abscissae assigned to that particular function, and these must
  // self-consistently yield a uniform grid in (T,Z) space.
  double totlum,bandlum;
  double e1 = js->efit1;
  double e2 = js->efit2;

  for (i = 0; i < js->tgridsize; ++i) 
    for (j = 0; j < js->mgridsize; ++j) {
      totlum = bandlum = 0;
      read_e(coolfile,js->mastercooling[i][j],js->egridsize*sizeof(float));
      for (k = 0; k < js->egridsize; ++k) {
	if (js->energyaxis[k] >= e1 && js->energyaxis[k] <= e2)
	  bandlum += js->energyaxis[k]*js->mastercooling[i][j][k];
	if (js->energyaxis[k] >= js->lxe1 && js->energyaxis[k] <= js->lxe2)
	  totlum += js->energyaxis[k]*js->mastercooling[i][j][k];
      }
      js->bolometric[i][j] = totlum/bandlum;
    }

  close(coolfile);

  if (js->debug > 0) {
    printf("Done. Coolingfunction has %d x %d x %d elements\n",
	   js->egridsize,js->tgridsize,js->mgridsize);

    printf("Energy range: %E to %E keV\n",
	   js->energyaxis[0],js->energyaxis[js->egridsize-1]);

    printf("Temperature range: %E to %E keV\n",
	   js->tempaxis[0],js->tempaxis[js->tgridsize-1]);

    printf("Metallicity range: %E to %E keV\n",
	   js->metalaxis[0],js->metalaxis[js->mgridsize-1]);
  }

  js->coolaccel = gsl_interp_accel_alloc();
}

// Morrison & McCammon cross section in units of 10^-22 cm^2.
int photoabs(double nH, double *ener, double dener, int ne, double *abscoef)
{
  int i;
  double inve,c1,c2,c3,c4,de1,de2,dsum;
  unsigned long pos;
  
  double pe[14] = {0.03,0.1,0.284,0.4,0.532,0.707,0.867,1.303,1.840,2.471,
		  3.210,4.038,7.111,8.331};
  double p0[14] = {17.3,34.6,78.1,71.4,95.5,308.9,120.6,141.3,202.7,342.7,
		  352.2,433.9,629.0,701.2};
  double p1[14] = {608.1,267.9,18.8,66.8,145.8,-380.6,169.3,146.8,104.7,
		  18.7,18.7,-2.4,30.9,25.2};
  double p2[14] = {-2150.,-476.1,4.3,-51.4,-61.1,294.0,-47.7,-31.5,-17.0,
		  0.,0.,0.75,0.,0.};

  pos = gsl_interp_bsearch(pe,ener[0],0,13);

  for (i = 0; i < ne; ++i) {
    inve = 1./ener[i]; 
    c1 =  exp(-nH*inve*(p2[pos]+inve*(p1[pos]+inve*p0[pos]))/100.);
    //printf("%9.4f%9.4f",ener[i],c1);
    if (pos < 13 && ener[i]+dener > pe[pos+1]) {
      inve = 1./pe[pos+1];
      c3 =  exp(-nH*inve*(p2[pos]+inve*(p1[pos]+inve*p0[pos]))/100.);
      ++pos;
      c4 =  exp(-nH*inve*(p2[pos]+inve*(p1[pos]+inve*p0[pos]))/100.);
      de1 = pe[pos]-ener[i];
      de2 = dener-de1;
      dsum = 2.*dener;
    } else {
      de1 = de2 = 1.; dsum = 2.;
      c3 = c4 = 0.; 
    }
    inve = 1./(ener[i]+dener); 
    c2 =  exp(-nH*inve*(p2[pos]+inve*(p1[pos]+inve*p0[pos]))/100.);
    
    abscoef[i] = ((c1+c3)*de1+(c2+c4)*de2)/dsum;
  }
  return 0;
}
    
// Take the finely-grained cooling function in mastercooling
// and rebin it to match the detector binning given in binlo and binhi
int rebincoolingfunction(double *bincenter, int nbin, struct jaco_state *js)
{
  long      i,j,k;
  unsigned long e1pos,e2pos;
  double    mincool = 1.E30,e1fac,e2fac,e1,e2,logmincool=-1E30,countbinsize;
  int      haschanged=0,l;

  js->l1 = js->l2 = -1;
  // If this is the first time we're being called, nlastbin will be negative
  
  if (js->debug > 7) {
    if (js->nlastbin == 0)
      for (i = 0; (int)i < nbin; ++i)
	printf("%f\n",bincenter[i]);
    else for (i = 0; i < js->nlastbin; ++i)
      printf("%f %f\n",bincenter[i],js->lastbin[i]);
  }

  if (js->nlastbin > 0) {

    // The array lastbin holds the latest set of bins in which the
    // cooling function has been calculated. This array has nlastbin
    // elements. Now, given a new series of bins stored in the
    // bincenter array, we want to find whether bincenter is a
    // subset of lastbin. If so, haschanged=0, and no rebinning is
    // required, saving us a lot of CPU time. Otherwise, haschanged=1,
    // and the cooling function must be recalculated..

    // Find l1 and l2 such that lastbin[l1..l2] is the same array
    // as bincenter[1..nbin]. This is done in four steps


    // A - find l1
    for (i = 0; i < js->nlastbin && js->l1 == -1; ++i) 
      if (fabs(bincenter[0]-js->lastbin[i]) < 1.E-4) 
	js->l1 = i;
	
    // B - find l2
    for (i = js->nlastbin-1; i >= 0 && js->l2 == -1; --i) 
      if (fabs(bincenter[nbin-1]-js->lastbin[i]) < 1.E-4) 
	js->l2 = i;
	
    i = 0;
    if (js->debug > 2)
      printf("New limits: %d,%d (%E %E) (%d bins)\n",
	     js->l1,js->l2,bincenter[0],bincenter[nbin-1],nbin);

    // C - check that l1..l2 are specify the same number of
    //     elements as 1..nbin
    if (js->l2-js->l1+1 == nbin && js->l1 > -1 && js->l2 > -1)
      // D - check that each nlastbin[l1..l2] is equal to
      //     bincenter[1..nbin]
      for  (l = js->l1+1; l < js->l2 && !haschanged; ++l) 
	haschanged = (fabs(bincenter[++i]-js->lastbin[l]) > 1.E-4);

    else
      haschanged = 1;

    if (haschanged) {
      // Delete the last cooling function matrix.
      sci_free_dtensor(js->rebinnedcooling,js->nlastbin,js->tgridsize);
      if (js->lastbin) free(js->lastbin);
    } else {
      if (js->l2-js->l1+1 != nbin)
	printf("WARNING----mismatch in energy grid.\n");
    }
  } else haschanged=1;

  if (haschanged) {

    ++js->nrebin;
    // Generate a new coolingfunction matrix
    js->l1 = 0; js->l2 = nbin-1;
    js->rebinnedcooling = sci_dtensor(nbin,js->tgridsize,js->mgridsize);
    js->lastbin = sci_dvector(nbin);

    if (js->debug > 0) 
      printf("Rebinning cooling function: %d bins requested.\n",nbin);
    if (js->nrebin > 5)
      printf("...Inoptimal binning slows calculations---rebin the ARF/RMF.\n");
    for  (l = 0; l < nbin; ++l) 
      js->lastbin[l] = bincenter[l];
    js->nlastbin = nbin;
  } else 
    // We have checked that rebinnedcooling[1..nbin] must be a subset
    // of the previously calculated rebinnedcooling[1..nlastbin];
    // specifically, it equals rebinnedcooling[l1..l2]. Time to exit.

    return 0;
      
  // We assume that all the bin sizes are the SAME, and that
  // the grid is REGULAR. 
  countbinsize = (bincenter[1]-bincenter[0]);
  if (js->debug > 0)
    printf("spectral bin size is %f keV\n",countbinsize);
  if (js->debug > 1)
    printf("Bins start at %f and end at %f\n",
	   bincenter[0],bincenter[nbin-1]);
  
  for (l = 0; l < nbin; ++l) {

    // Calculate the rest-frame bin corresponding to the given
    // energy range
    e1 = js->opz*(bincenter[l]); 
    e2 = js->opz*(bincenter[l]+countbinsize);

    if (e1 < js->energyaxis[0] || e2 < js->energyaxis[0] || 
	e1 > js->energyaxis[js->egridsize-1] || 
	e2 > js->energyaxis[js->egridsize-1]) {
      printf("For bounds: %f--%f\n",e1,e2);
      printf("Energy bounds of cooling matrix were exceeded.\n");
    }

    // Find the position of the bins corresponding to this range
    e1pos = gsl_interp_accel_find(js->coolaccel,
				  js->energyaxis,js->egridsize,e1);
    e2pos = gsl_interp_accel_find(js->coolaccel,
				  js->energyaxis,js->egridsize,e2);

    if (e1pos > e2pos)
      BYE("Lower energy bound %ld > upper energy bound.",e1pos);

    if (e1pos != e2pos) {
      // Upper and lower energy limits occur in different bins
      e1fac = (js->energyaxis[e1pos+1]-e1);
      e2fac = (e2-js->energyaxis[e2pos]);
    } else {
      // Upper and lower energy limits occur in the same bin
      // (not desireable, since this means the cooling function is not
      // binned with fine enough resolution)
      e1fac = (e2-e1);
      e2fac = 0;
    }

    //printf("%E %E %E %E\n",e1,e2,e1fac,e2fac);
    // Add together the contributions from the starting and ending bins
    for (j = 0; j < js->tgridsize; ++j)
      for (k = 0; k < js->mgridsize; ++k) {
	js->rebinnedcooling[l][j][k] = 
	  js->mastercooling[j][k][e1pos]*e1fac+
	  js->mastercooling[j][k][e2pos]*e2fac;
	if (js->rebinnedcooling[l][j][k] > 0. &&
	    js->rebinnedcooling[l][j][k] < mincool)
	  mincool = js->rebinnedcooling[l][j][k];
      }
    
    
    // Add together the contributions from the intermediate bins
    for (i = (long)e1pos+1; i < (long)e2pos; ++i)
      for (j = 0; j < js->tgridsize; ++j)
	for (k = 0; k < js->mgridsize; ++k) {
	  js->rebinnedcooling[l][j][k] += 
	    js->mastercooling[j][k][i]*js->ebinsize;
	  if (js->rebinnedcooling[l][j][k] > 0. &&
	      js->rebinnedcooling[l][j][k] < mincool)
	    mincool = js->rebinnedcooling[l][j][k];
	}
    
    logmincool = log10(mincool);
    
    // Now take the log of the rebinned function, for interpolation in
    // log space. Truly null cooling rates are simply set to the
    // minimum rate seen altogether.
    for (j = 0; j < js->tgridsize; ++j)
      for (k = 0; k < js->mgridsize; ++k) {
	js->rebinnedcooling[l][j][k] = 
	  (js->rebinnedcooling[l][j][k] < mincool ? logmincool :
	   log10(js->rebinnedcooling[l][j][k]));
      }
    
  }

  if (js->xymodel == NULL)
    js->xymodel = sci_dtensor(nbin,js->xnx,js->xny);

  return 1;
  
}


#if defined(WITH_SZ) || defined(WITH_BOLOCAM)

void jaco_sz(int i, struct jaco_vars *jv)

{

  jv->R1 = jv->state->rspline[i];
  jv->R2 = jv->state->rspline[i+1];
  jv->R1sq = jv->R1*jv->R1;
  jv->R2sq = jv->R2*jv->R2;

  double szint = jv->state->ninteg(szintegrand,jv->R1,jv->state->rshock,jv);
  
  jv->state->szy[i] =
    log10(4.*SZNORM*MPCM3*szint/(jv->R2sq-jv->R1sq));

  jv->state->szrad[i] = 
    log10((jv->R1+jv->R2)/(2.*ARCMINTORAD*jv->state->angdist*60));

}

void generate_ymap(struct jaco_state *js)
{
  int i,j,obs,nannuli;
  double halfwidth,rad,xrad,yrad,Asz,bsz,xoff,yoff,raoffset,decoffset;
  gsl_interp *y_coefs;
  gsl_interp_accel *y_accel;
  FILE *outf=NULL;
  BigSZControlParams szparams;
  char on[50];
  struct jaco_vars jv;

  szparams = js->szparams;
  JS(debug); 


  sprintf(on,"sz%03d",js->mpirank);
  if (debug > 0)
    outf = fopen(on,"w+");

  if (js->dofast)
    nannuli = js->lastszannulus;
  else {
    nannuli = js->nsplines-1;
  }

  if (debug > 0) {
    for (i = 0; i < nannuli; ++i)
      fprintf(outf,"%f %f\n",js->szrad[i],js->szy[i]);
  }

  i = nannuli;

  bsz = (js->szy[i-1]-js->szy[i-2])/(js->szrad[i-1]-js->szrad[i-2]);
  Asz = pow(10.,js->szy[i-2]-bsz*js->szrad[i-2]);

  int k;
 
  //41/50 vs. 46/50

  for (obs = 0; obs < szparams.ndataset; ++obs) {

    y_coefs = gsl_interp_alloc(gsl_interp_cspline,nannuli);
    gsl_interp_init(y_coefs,js->szrad,js->szy,nannuli);

    float *buffer = sci_fvector(szparams.npix[obs]*szparams.npix[obs]);
    k = 0;

    raoffset = (js->xrayra-szparams.racent[obs])*cos(js->xraydec*PI/180.);
    decoffset = szparams.deccent[obs]-js->xraydec;

    if (szparams.npix[obs] % 2 == 0) {
      xoff = 0.5+raoffset/szparams.pixsize[obs];
      yoff = 0.5+decoffset/szparams.pixsize[obs];
    } else {
      xoff = raoffset/szparams.pixsize[obs];
      yoff = decoffset/szparams.pixsize[obs];
    }
    halfwidth = szparams.npix[obs]/2;

    gsl_set_error_handler_off();
#pragma omp parallel  default(none) private(jv,i,j,k,rad,xrad,yrad,y_accel) shared(xoff,yoff,obs,js,szparams,debug,stdout,nannuli,bsz,Asz,gsl_interp_cspline,halfwidth,raoffset,decoffset,y_coefs)

    {
      jaco_thread_init(&jv);
      jv.state = js;
      jv.integrand_ordinate = js->pressure;
      jv.integrand_coefs = js->szinterp_coefs;
      y_accel = jv.accel;
      jv.dolog = 1;
      jv.workspace[0] = gsl_integration_workspace_alloc(INTMAX);

      #pragma omp for  schedule(dynamic)
      for (i = 0; i < szparams.npix[obs]; ++i) {

	for (j = 0; j < szparams.npix[obs]; ++j) {
	  xrad = (double)(i-halfwidth)+xoff;
	  yrad = (double)(j-halfwidth)+yoff;
	  if (js->gasmodel == TRIAXIAL) {
	    jv.x = xrad*szparams.pixsize[obs]*60.*ARCMINTORAD*js->angdist;
	    jv.y = yrad*szparams.pixsize[obs]*60.*ARCMINTORAD*js->angdist;
	    szparams.ymap[obs][i][j] = js->ninteg(integrand_xy,-js->rshock,0,&jv)+
	      js->ninteg(integrand_xy,0,js->rshock,&jv);
	    szparams.ymap[obs][i][j] *= SZNORM*MPCM3;
	  } else {
	    rad = sqrt(xrad*xrad+yrad*yrad);
	    if (gsl_fcmp(rad,0.25,1.E-4) < 0)
	      rad = 0.25;
	    rad = log10(rad*szparams.pixsize[obs]);
	    //fflush(stdout);
	    if (rad > js->szrad[nannuli-1])
	      szparams.ymap[obs][i][j] = Asz*pow(10.,bsz*rad);
	    else
	      szparams.ymap[obs][i][j] = 
		pow(10.,gsl_interp_eval(y_coefs,js->szrad,js->szy,rad,y_accel));
	  }
	  
	  /* 	if (debug > 9) { */
	  /* 	  fprintf(outf," %9E\n", */
	  /* 		  //js->xrayra,szparams.racent[obs],js->xraydec, */
	  /* 		  //szparams.deccent[obs], */
/* 		  szparams.ymap[obs][i][j]); */
/* /\* 	  fflush(outf); *\/ */
/* 	} */
	}
      }
      gsl_integration_workspace_free(jv.workspace[0]);
      jaco_thread_shutdown(&jv);
    }
    gsl_set_error_handler(NULL);

    for (i = 0; i < szparams.npix[obs]; ++i)
      for (j = 0; j < szparams.npix[obs]; ++j) 
	buffer[k++] = szparams.ymap[obs][j][i];

    if (js->calcsz > 1) {
      char *imfile = sci_strcat(js->profname,".ymap.fits");
      unlink(imfile);
      hrothgar_writeimage(imfile,buffer,
			  (-halfwidth)*szparams.pixsize[obs],
			  halfwidth*szparams.pixsize[obs],
			  (-halfwidth)*szparams.pixsize[obs],
			  halfwidth*szparams.pixsize[obs],
			  szparams.npix[obs]);
    }

    free(buffer);
    gsl_interp_free(y_coefs);
  }

  if (debug > 0)
    fclose(outf);
}
#endif


// Everything must be in radians for WCS routines.
void xytowcs(double x, double y, double *ra, double *dec, struct jaco_state *js)
{
  double xi, eta, crpix1, crpix2, cdelt1, cdelt2, ra0, dec0,sinra,cosra;
  JS(cosdec0); JS(sindec0); JS(cosra0); JS(sinra0);

  crpix1 = js->wcs[0][0]; ra0 = js->wcs[1][0]; cdelt1 = js->wcs[2][0];
  crpix2 = js->wcs[0][1]; dec0 = js->wcs[1][1]; cdelt2 = js->wcs[2][1];

  //printf("*** %E %E %E    %E %E %E ***\n",crpix1,ra0,cdelt1,crpix2,dec0,cdelt2);
  xi = (x-crpix1)*cdelt1;
  eta = (y-crpix2)*cdelt2;

  *ra = ra0+atan2(xi,(cosdec0-eta*sindec0));
  *dec = asin((sindec0+eta*cosdec0)/sqrt(1+xi*xi+eta*eta));


  return;
  
}


void wcstoxy(double ra, double dec, double *x, double *y, struct jaco_state *js)
{
  double xi, eta, crpix1, crpix2, cdelt1, cdelt2, ra0, dec0, sindec, cosdec;
  double dra, sindra, cosdra;
  JS(cosdec0); JS(sindec0);

  crpix1 = js->wcs[0][0]; ra0 = js->wcs[1][0]; cdelt1 = js->wcs[2][0];
  crpix2 = js->wcs[0][1]; dec0 = js->wcs[1][1]; cdelt2 = js->wcs[2][1];

  sincos(dec,&sindec,&cosdec);
  dra = ra-ra0;
  sincos(dra,&sindra,&cosdra);

  double denom = sindec*sindec0+cosdec*cosdec0*cosdra;
  xi = cosdec*sindra/denom;
  eta = (sindec*cosdec0-cosdec*sindec0*cosdra)/denom;

  *x = crpix1+xi/cdelt1;
  *y = crpix2+eta/cdelt2;

  return;
  
}

double wcsdist(double ra1, double dec1, double ra2, double dec2) 
{
  double sinra1, sinra2, cosra1, cosra2;

  sincos(ra1,&sinra1,&cosra1);
  sincos(ra2,&sinra2,&cosra2);
  return acos(sinra1*sinra2+cosra1*cosra2*cos(dec1-dec2));
}

void generate_xmap(struct jaco_state *js,double **spectrum)
{
  int b,i,j,obs,nannuli;
  double xcen,ycen,ra,dec;
  gsl_interp_accel *x_accel;
  FILE *outf=NULL;
  char on[50];
  struct jaco_vars jv;

  JS(debug); 

  int k;
 
  //41/50 vs. 46/50
  
  k = 0;
  
  double xrayra = PIOVER180*js->xrayra;
  double xraydec = PIOVER180*js->xraydec;

  wcstoxy(xrayra,xraydec,&xcen,&ycen,js);
  //printf("------%E      %E -----------\n",xcen,ycen);

  gsl_set_error_handler_off();
  
  int pixinring[MAXANN];
  int **isinannulus = sci_imatrix(js->xnx,js->xny);
  double **dist = sci_dmatrix(js->xnx,js->xny);
  double **xval = sci_dmatrix(js->xnx,js->xny);
  double **yval = sci_dmatrix(js->xnx,js->xny);
  double pixarea = js->xpixscale*js->angdist*PIOVER180/3600.;
  pixarea *= pixarea;

  #pragma omp parallel for schedule(dynamic) default(none) private(i,j,ra,dec) shared(isinannulus,xval,yval,xrayra,xraydec,js,pixinring,dist,pixarea,xcen,ycen)
  for (i = 0; i < js->xnx; ++i) 
    for (j = 0; j < js->xny; ++j) {
      isinannulus[i][j] = -1;
      //printf("%d %d %E %E\n",i,j,ra/PIOVER180,dec/PIOVER180);
      //jv.x = wcsdist(ra,dec,xrayra,dec)*js->angdist;
      //jv.y = wcsdist(ra,dec,ra,xraydec)*js->angdist;
      xval[i][j] = (xcen-1.*i)*js->xpixscale;
      yval[i][j] = (1.*j-ycen)*js->xpixscale;
      dist[i][j] = sqrt(xval[i][j]*xval[i][j]+yval[i][j]*yval[i][j])/60.;
      xval[i][j] *= js->angdist*PIOVER180/3600.;
      yval[i][j] *= js->angdist*PIOVER180/3600.;
      js->softmodel[i][j] = 0.;
    }

#pragma omp parallel for schedule(dynamic) default(none) private(i,j,k,ra,dec,b) shared(isinannulus,xval,yval,xrayra,xraydec,js,pixinring,dist,spectrum,pixarea)
  for (k = 0; k < js->nannuli; ++k) {
    pixinring[k]=0;
    for (i = 0; i < js->xnx; ++i) 
      for (j = 0; j < js->xny; ++j) 
	if (js->softexp[i][j] > 0 && js->softimage[i][j] > SOFTIMAGE_MIN) 
	  if (dist[i][j] > js->r1list[k] && dist[i][j] < js->r2list[k]) {
	    ++pixinring[k];
	    isinannulus[i][j] = k;
	  } 
    for (b = 0;  b < js->nlastbin; ++b)
      spectrum[k][b] = 0.;
  }

#pragma omp parallel  default(none) private(jv,i,j,b,k,ra,dec) shared(ycen,xcen,js,debug,stdout,nannuli,xrayra,xraydec,spectrum,isinannulus,xval,yval,pixinring,pixarea,dist)

    {
      jaco_thread_init(&jv);
      jv.state = js;
      jv.workspace[0] = gsl_integration_workspace_alloc(INTMAX);
      jv.workspace[1] = gsl_integration_workspace_alloc(INTMAX);
      jv.dolog = 0;
      //abscissae(0,0.1,&jv);
	  
#pragma omp for schedule(dynamic)
      for (b = 0; b < js->nlastbin; ++b) {
	jv.integrand_ordinate = js->emissivity[b];
	jv.integrand_coefs = js->interp_coefs[b];
	for (i = 0; i < js->xnx; ++i) 
	  for (j = 0; j < js->xny; ++j) {
	    if (js->softexp[i][j] > 0 && js->softimage[i][j] > SOFTIMAGE_MIN) {
	      jv.x = xval[i][j];
	      jv.y = yval[i][j];
	      js->xymodel[b][i][j] = js->ninteg(integrand_xy,-js->rshock,js->rshock,&jv);
	      //js->ninteg(integrand_xy,0,js->rshock,&jv);
	    } else
	      js->xymodel[b][i][j] = 0.;
	  }
      }
      gsl_integration_workspace_free(jv.workspace[0]);
      gsl_integration_workspace_free(jv.workspace[1]);
      gsl_set_error_handler(NULL);

#pragma omp for schedule(dynamic)
      for (b = 0; b < js->nlastbin; ++b)  {
	for (i = 0; i < js->xnx; ++i) 
	  for (j = 0; j < js->xny; ++j) 
	    if (isinannulus[i][j] > -1)
	      spectrum[isinannulus[i][j]][b] += js->xymodel[b][i][j];

	for (k = 0; k < js->nannuli; ++k) {
	  //printf("%d %e %e %e %e %e \n",pixinring[k],js->r1list[k],js->r2list[k],js->xpixscale,js->ringarea[k],(pixinring[k]*js->xpixscale*js->xpixscale/3600.));

	  if (pixinring[k] > 0) {
	    spectrum[k][b] *= pixarea*js->specnorm/(4*PI);
	    spectrum[k][b] *= 
	      js->ringarea[k]/
	      (pixinring[k]*js->xpixscale*js->xpixscale/3600.);

	  }
	}
      }

#pragma omp for schedule(dynamic)
      for (i = 0; i < js->xnx; ++i) 
	for (j = 0; j < js->xny; ++j)
	  for (b = 0; b < js->nlastbin; ++b)  {
	    if (js->softexp[i][j] > 0 && js->softimage[i][j] > SOFTIMAGE_MIN)
	      if (js->lastbin[b] > js->xime1 && js->lastbin[b] < js->xime2)
		js->softmodel[i][j] += js->xymodel[b][i][j];
	  }
	  
      jaco_thread_shutdown(&jv);
    }


    if (debug > 0) {
      float *buffer = sci_fvector(js->xnx*js->xny);
      float imwcs[3][2];
      for (j = 0; j < 2; ++j) {
	imwcs[0][j] = js->wcs[0][j];
	for (i = 1; i < 3; ++i) 
	  imwcs[i][j] = js->wcs[i][j]/PIOVER180;
      }
      
      k = 0;
      for (i = 0; i < js->xnx; ++i)
	for (j = 0; j < js->xny; ++j) {
	  ++k;
	  buffer[k] += js->softmodel[j][i];
	}
      
      char *imfile = sci_strcat(js->profname,".xmap.fits");
      unlink(imfile);
      hrothgar_wcsimage(imfile,buffer,
			imwcs,
			js->xnx);
      free(buffer);
    }


    sci_free_imatrix(isinannulus,js->xnx);
    sci_free_dmatrix(dist,js->xnx);
    sci_free_dmatrix(xval,js->xnx);
    sci_free_dmatrix(yval,js->xnx);
}


double wl_total(double r, void *par)
{
  double intR1sq,intR2sq,answer,rsq;
  struct jaco_state *js = (struct jaco_state *)par;

  rsq = r*r;
  intR1sq = js->intR1sq; intR2sq = js->intR2sq;
  
  answer = sqrt(rsq-intR1sq);
  if (rsq > intR2sq) 
    answer -= sqrt(rsq-intR2sq);

  return r*answer*massdensity(r,js);
}

double wl_surf(double z, void *par)
{
  double zsq,r,Rsq;
  struct jaco_state *js = (struct jaco_state *)par;

  zsq = z*z;
  Rsq = js->intR2sq; 
  
  r = sqrt(zsq+Rsq);

  return massdensity(r,js);
}


// Returns the projected mass divided by FOURPI
double projected_mass(double R1, double R2, struct jaco_state *js)
{
  gsl_function wl_total_intfunc;
  gsl_integration_workspace *wl_integration_workspace=NULL;
  double acc,result,surfacemass;
  int status;

  JS(precision);

  wl_total_intfunc.function = &wl_total;
  wl_total_intfunc.params = js;

  js->intR1sq = R1*R1;
  js->intR2sq = R2*R2;

  wl_integration_workspace=gsl_integration_workspace_alloc(INTMAX);
  status = 
    gsl_integration_qags(&wl_total_intfunc,R1,R2,0,precision,INTMAX,
			 wl_integration_workspace,&result,&acc);

  if (status != GSL_SUCCESS)
    printf("1)Trouble converging in surface mass between %e and %e\n",
	   R1,R2);
  surfacemass = result;

/*   status = gsl_integration_qags(&wl_total_intfunc,R2,js->rshock,0,precision, */
/* 				INTMAX, */
/* 				 wl_integration_workspace,&result,&acc); */
  status = gsl_integration_qagiu(&wl_total_intfunc,R2,0,precision,
				INTMAX,
				 wl_integration_workspace,&result,&acc);
  if (status != GSL_SUCCESS)
    printf("2)Trouble converging in surface mass between %e and infinity: %E %E %s\n", js->rdm1, js->Md,
	   R2,gsl_strerror(status));

  surfacemass += result;
  gsl_integration_workspace_free(wl_integration_workspace);

  return surfacemass;
}

void tangential_shear(double *R, double *shear, int nrad, struct jaco_state *js)
{
  double surfacemass = 0,Rsq,result,acc,oldR=0.,oldRsq=0.,*wlMpcR;
  double sigmacrit,kappa;
  gsl_function wl_diff_intfunc;
  gsl_integration_workspace *wl_integration_workspace=NULL;
  int i,status;

  JS(wlbeta); JS(wlcorrfac); JS(angdist); JS(precision);

  sigmacrit = SIGMAC/(wlbeta*angdist);
  
  wlMpcR = sci_dvector(nrad);

  wl_integration_workspace=gsl_integration_workspace_alloc(INTMAX);
  wl_diff_intfunc.function = &wl_surf;
  for (i = 0; i < nrad; ++i) 
    wlMpcR[i] = R[i]*ARCMINTORAD*angdist;
    
  js->projtype = PROJ_ALL;
  gsl_set_error_handler_off();
  for (i = 0; i < nrad; ++i) {
    Rsq = wlMpcR[i]*wlMpcR[i];
    //printf("%d %f %f\n",i,sigmacrit,Rsq);
    wl_diff_intfunc.params = js;
    //printf("1) OK %f %f %d %d\n",oldR,R[i],i,nrad);

    surfacemass = projected_mass(0.,wlMpcR[i],js);

    js->intR1sq = oldRsq;
    js->intR2sq = Rsq;

    shear[i] = 4.*surfacemass/Rsq;
    //printf("3) OK %f\n",result);
    /*    status = gsl_integration_qags(&wl_diff_intfunc,0.,js->rshock,0,precision,INTMAX,
	  wl_integration_workspace,&result,&acc);*/
    status = gsl_integration_qagiu(&wl_diff_intfunc,0.,0,precision,
				   INTMAX,
				   wl_integration_workspace,&result,&acc);
    if (status != GSL_SUCCESS)
      printf("3)Trouble converging in shear calculation between %e and %e\n",
	     oldR,wlMpcR[i]);
    //printf("4) OK %f\n",result);

    kappa = 2.*result;
    shear[i] -= kappa;
    shear[i] /= sigmacrit;
    kappa /= sigmacrit;
    shear[i] /= (1.-kappa);
    shear[i] /= 1.+wlcorrfac*kappa;
    oldR = wlMpcR[i];
    oldRsq = Rsq;
  }
  gsl_set_error_handler(NULL);
  free(wlMpcR);
  gsl_integration_workspace_free(wl_integration_workspace);
}

// Calculates the spectrum obtained by integrating the observed light
// between R1 and R2. This is done over several energy bins. 

 void getspectrum(long annulus, double *spec,  
			int nbin, struct jaco_vars *jv)
{
  long i;
  clock_t clin1;
  double szint,err;
    

  JVS(precision); 
  JVS(debug); JVS(dofast); JV(breakradius); 
  JVS(rshock); JVS(lumdistancesq);
  JVS(calcsz); JV(R1deg); JV(R2deg); JV(R1); JV(R2);

  err = precision;

  if (debug > 0)
    clin1 = clock();
  if (dofast) {
    abscissae(R1,breakradius,jv);
    emintegrand_fastfit(jv);
  }

  strcpy(jv->integstate,"full spectrum");
  for (i = 0; i < nbin; ++i) {

    // The 4pi in the projection integral cancels the four pi
    // in the denominator of this expression
    jv->coolmatrix = jv->state->rebinnedcooling[i];
    if (dofast) {
      spec[i] = 
	jv->state->ninteg(emintegrand_fastev,R1,R2,jv);
    }
    else {
      jv->specbin = i;
      //if (R1 < jv->state->rspline[0]) R1 = jv->state->rspline[0];
      spec[i] =
	jv->state->ninteg(emintegrand,R1,rshock,jv);
    }
    
    if (!FINITE(spec[i])) {
      if (1) 
	printf("Error encountered: %ld %ld, %E %E %E\n",annulus,i,R1,R2,rshock);
      spec[i] = 1.E30;
    }

    if (debug > 6)
      printf("annulus %E %f %ld %ld %E %f\n",rshock,R2,
	     annulus, i,spec[i],lumdistancesq);
    
  }

  if (calcsz && dofast) {
    strcpy(jv->integstate,"sz");
    szint = jv->state->ninteg(szintegrand_fastev,R1,R2,jv);
    jv->state->szy[annulus] =
      log10(4.*SZNORM*MPCM3*szint/(R2*R2-R1*R1));
    jv->state->szrad[annulus] = 
      log10((R1deg+R2deg)/2.);
  }

}

double mdiffcalc(double r, void *par)
{  
  double mofr,contr,answer;
  struct jaco_state *js = (struct jaco_state *)par;

  JS(Md); JS(rhocrit); JS(Ms); JS(MsContr);

  contr = js->contr;
  mofr = Md*js->darkmass(r,js);

  if (!js->totalmass) {
    mofr += gasmass(r,js);

    if (MsContr > 0.) 
      mofr += Ms*js->stellarmass(r,js);
  }

  answer = 3.*mofr/(FOURPI*r*r*r*rhocrit)-contr;
  //printf("---- **%f %f %f %f %f\n",r,Md,js->darkmass(r,js),r,answer);

  return answer;
}

/* double mcdiffcalc(double c, void *par) */
/* { */
/*   struct jaco_state *js = (struct jaco_state *)par; */

/*   double M200 = js->mcm200fac*pow(c,js->mcm200slope); */
/*   double c200 = js->mcA*pow(M200/js->mcM0,js->mcslope); */

/*   return js->Mtot/M200- */

/* } */

double overdensity_radius(double od, struct jaco_state *js)
{
  gsl_function F;
  int i=0,status;
  double r1,r2;

  F.function = mdiffcalc;
  F.params = js;
  js->contr = od;

  // Typically we will be looking for r{2500} to r{100}, so this should
  // a reasonable guess for the bracketing
  if (od > 1.e4) {
    r1 = 0.001;
    r2 = 0.1;
  } else {
    r1 = 0.1;
    r2 = 5.;
  }

  // Return too large a radius if contrast is too small
  if (od < 10) return 5;
  
  // Bracket the root if it's not already
  while (mdiffcalc(r1,js)*mdiffcalc(r2,js) > 0.) {
    r1 /= 1.5;
    r2 *= 2;
    ++i;
    if (i > 5) {
      printf("Warning: Could not bracket overdensity radius r_%f.\n",od);
      return r2;
    }
  }
  
/*   printf("---- %f %f\n",rContr,mdiffcalc(1.6,&contrast)); */
/*   printf("---- %f %f\n",rContr,mdiffcalc(rContr,&contrast)); */
/*   printf("---- %f %f\n",1.7,mdiffcalc(1.7,&contrast)); */
/*   printf("---- %f %f\n",1.8,mdiffcalc(1.8,&contrast)); */
/*   printf("---- %f %f\n",1.9,mdiffcalc(1.9,&contrast)); */
/*   printf("---- %f %f\n",2.0,mdiffcalc(2.0,&contrast)); */
/*   printf("---- %f %f\n",2.1,mdiffcalc(2.1,&contrast)); */

  i = 0; 
  gsl_root_fsolver_set(js->solver,&F,r1,r2);
  do {
    ++i;
    status = gsl_root_fsolver_iterate(js->solver);
    r1 = gsl_root_fsolver_x_lower(js->solver);
    r2 = gsl_root_fsolver_x_upper(js->solver);
    //printf("Hi %E %E %E\n",r1,r2,od);
    //printf("---%f %f %f %f\n",r1,r2,mdiffcalc(r1,js),mdiffcalc(r2,js));
    //printf("Hi\n");
    status = gsl_root_test_delta(r1,r2,0.,0.0001);
  } while (status == GSL_CONTINUE && i < 500);

  if (status != GSL_SUCCESS)
    BYE("Could not find overdensity radius r_%f.\n",od);

  return r1;
    
}


// Initialize the various variables according to the params[] array.
void jaco_update_models(double *params, struct jaco_state *js)
{
  double rcubed,r200,shock_overdensity,maxrad,step,Mtot;
  int *changedpar;
  long i;
 
  changedpar = js->changedpar;

  JS(opz); JS(rhocrit); JS(dist); JS(angdist); JS(lumdistancesq);
  JSET(redshift,params[REDSHIFT]);
  JSET(contrast,params[CONTRAST]);

  if (changedpar[REDSHIFT]) {
    JRESET(opz,1.+redshift);
    JRESET(rhocrit,RHOC*js->hubble*js->hubble*
      (opz*opz*opz*js->omegam + js->omegal));
    JRESET(dist,ztoMpc(js->hubble,js->omegam,-1.,redshift));
    JRESET(angdist,dist/opz);
    JRESET(lumdistancesq,dist*dist*opz*opz);
  } 
  
#ifndef NONGRAVITY
  JSET(tempwarn,0);
  JSET(metalwarn,0);
  JSET(entropywarning,0);
  JSET(Mg1Contr,params[MG1]); 
  JSET(alpha,params[ALPHA]);
  JSET(b1,params[BETA1]);
  JSET(b2,params[BETA2]);
  JSET(b3,params[BETA3]);
  JSET(rx1,params[RX1]);
  JSET(rx2,params[RX2]);
  JSET(rx3,params[RX3]);
  JSET(rz,params[RZ]);
  JSET(Mg2Contr,params[MG2]);
  JSET(Mg3Contr,params[MG3]);
  JSET(MdContr,params[DARKMASS]);
  JSET(ndark,params[DARKSLOPE]);
  JSET(rdm0,params[RDM0]);
  JSET(rdm1,params[RDM1]);
  JSET(z0,params[Z0]);
  JSET(zinf,params[ZINF]);
  JSET(MsContr,params[STARMASS]);
  JSET(nstar,params[STARSLOPE]);
  JSET(rstar,params[RSTAR]);
  JSET(aniso0,params[ANISO0]);
  JSET(aniso1,params[ANISO1]);
  JSET(raniso,params[RANISO]);
  JSET(xrayra,params[XRAYRA]);
  JSET(xraydec,params[XRAYDEC]);
  JSET(ximnorm,params[XIMNORM]);
  JSET(ximback,params[XIMBACK]);
  JSET(nonthermal,params[NONTHERMAL]);
  shock_overdensity = params[RSHOCK];
  JSET(tshock,params[SHOCKQ]);
  JS(gasmodel);
  JS(darkmodel);
  JS(totalmass);

  JS(rContr); JS(dofast); JS(ninteg); JS(nsplines);
  JS(rx3sq);   JS(precision);
  JS(rx2sq); JS(rx1sq); 
  JS(rdm1cu); JS(darkmass); JS(darkdens); JS(rdm3mn); JS(rdm1sq);
  JS(rdiff); JS(rdm0sq); JS(stellarmass); JS(stellardensity);
  JS(Ms); JS(Md); JS(nm3); JS(rshock); JS(pshock); 
  double sinphi,cosphi,sintheta,costheta,sinpsi,cospsi;

  //triaxial parameters
  if (gasmodel == TRIAXIAL) {
    js->asq = params[AXIS1]*params[AXIS1];
    js->bsq = params[AXIS2]*params[AXIS2];
    sincos(PIOVER180*params[PHI],&sinphi,&cosphi);
    sincos(PIOVER180*params[THETA],&sintheta,&costheta);
    sincos(PIOVER180*params[PSI],&sinpsi,&cospsi);
    js->x_x = cosphi*cospsi-costheta*sinphi*sinpsi;
    js->x_y = sinphi*cospsi+costheta*cosphi*sinpsi;
    js->x_z = sintheta*sinpsi;
    js->y_x = -cosphi*sinpsi-costheta*sinphi*cospsi;
    js->y_y = -sinphi*sinpsi+costheta*cosphi*cospsi;
    js->y_z = sintheta*cospsi;
    js->z_x = sinphi*sintheta;
    js->z_y = -cosphi*sintheta;
    js->z_z = costheta;
  }

  double mtrack;

  if (totalmass) {
    Mtot = MdContr;
    if (js->fitfgas) 
      mtrack = Mg1Contr*Mtot;
    else
      mtrack = Mg1Contr;
  } else {
    if (js->fitfgas) {
      Mtot = (MsContr+MdContr)/(1.-Mg1Contr);
      mtrack = Mg1Contr*Mtot;
    } else {
      mtrack = Mg1Contr;
      Mtot = MsContr+MdContr+mtrack;
    }
  }

    //JRESET(MdContr,Mtot-mtrack-MsContr);
  
  JRESET(Mg2Contr,Mg2Contr*mtrack);
  mtrack -= Mg2Contr;
  JRESET(Mg3Contr,Mg3Contr*mtrack);
  mtrack -= Mg3Contr;
  JRESET(Mg1Contr,mtrack);
  
  if (contrast > 0.) {
    rcubed = 3.*Mtot/(FOURPI*rhocrit);
    JRESET(rContr,cbrt(rcubed/(contrast)));
  } else {
    JRESET(rContr,-contrast);
  }
  if (DIFFER(contrast,200.) && js->fitmc)
    BYE("Fixed MC relation mode only currently available for Contrast=200");
  

  double cnow;
  cnow = rdm1;
  // A negative DARKMODEL indicates that instead of the dark matter
  // scale radius, we would like to measure the concentration.

  if (darkmodel < 0) {

    // Duffy et al. (2008) MC relation
    // Currently this is only valid if contrast=200
    if (js->fitmc)
      cnow = 5.71*pow(100.*Mtot/2,-0.084)*pow(1.+js->redshift,-0.44);
      
    JRESET(rdm1,rContr/cnow);
    //rdm0 *= rContr;
    //rz *= rContr;
    //rx1 = rContr/rx1;
    //rx2 = rContr/rx2;
    //rx3 = rContr/rx3;
  }
  JSET(nm2,ndark-2.);
  JSET(nm1,ndark-1.);

  if (precision > 1.e-2 || (js->fitwrite && precision > 1.e-3)) {
    printf("Minimum precision of 1\% required; 0.1\% if writing profiles.\n");
    exit(-1);
  }

  if (js->dofast) {
    if (shock_overdensity > 0) {
      JRESET(ninteg,ninteg_shock);
    } else
      JRESET(ninteg,ninteg_nr);
    JRESET(nsplines,0);
  } else {
      JRESET(ninteg,ninteg_gsl);
      JRESET(nsplines,GSL_MIN(2./precision,300.));
  }
  
  JRESET(rx1sq,rx1*rx1);
  JRESET(rx2sq,rx2*rx2);
  JRESET(rx3sq,rx3*rx3);

  JSET(Mg1,Mg1Contr/gasmass_beta(rContr,js->alpha,rx1,rx1sq,b1));
  JSET(Mg2,Mg2Contr/gasmass_beta(rContr,0.,rx2,rx2sq,b2));
  JSET(Mg3,Mg3Contr/gasmass_beta(rContr,0.,rx3,rx3sq,b3));
  JSET(rhog1,2.*Mg1/(FOURPI*rx1*rx1sq));
  JSET(rhog2,2.*Mg2/(FOURPI*rx2*rx2sq));
  JSET(rhog3,2.*Mg3/(FOURPI*rx3*rx3sq));

  switch (abs(darkmodel)) {
    
  case MASS_HERNQUISTLIKE:
    JRESET(nm3,ndark-3.);
    JRESET(rdm1cu,rdm1*rdm1*rdm1);
    JRESET(darkmass,mass_hernquistlike);
    JRESET(darkdens,dens_hernquistlike);
    break;
    
  case MASS_POWLAWWITHCORE:
    JRESET(nm1,ndark-1.);
    JRESET(nm2,ndark-2.);
    JRESET(nm3,ndark-3.);
    JRESET(rdm1,rdm0);
    JRESET(rdm3mn,pow(rdm1,-nm3));
    JRESET(rdm1sq,rdm1*rdm1);
    JRESET(rdiff,nm3*nm2*nm1);
    JRESET(darkmass,mass_powlawwithcore);
    JRESET(darkdens,dens_powlawwithcore);
    break;
      
  case MASS_NFWWITHCORE:
    JRESET(rdm0sq,rdm0*rdm0);
    JRESET(rdiff,rdm1-rdm0);
    JRESET(darkmass,mass_nfwwithcore);
    JRESET(darkdens,dens_nfwwithcore);
    break;
    
  case MASS_ISOWITHCORE:
    JRESET(rdm0sq,rdm0*rdm0);
    JRESET(rdm1sq,rdm1*rdm1);
    JRESET(rdiff,rdm0-rdm1);
    JRESET(darkmass,mass_isowithcore);
    JRESET(darkdens,dens_isowithcore);
    break;
    
  case MASS_UNIVERSAL:
    JRESET(darkmass,mass_universal);
    JRESET(darkdens,dens_universal);
    break;
    
  case MASS_SERSIC:
    JRESET(darkmass,mass_sersic_dark);
    JRESET(darkdens,dens_sersic_dark);
    break;
    
  case MASS_FOLLOWSLIGHT:
    JRESET(darkmass,light);
    break;
    
  }

  JS(mass);
  JRESET(mass,allmass);

  if (MsContr > 0) {
    JRESET(stellarmass,mass_sersic_stellar);
    JRESET(stellardensity,dens_sersic_stellar);

    JRESET(Ms,MsContr/js->stellarmass(rContr,js));
  }

  JRESET(Md,MdContr/darkmass(rContr,js));

  if (shock_overdensity > 0.) {
    JRESET(rshock,overdensity_radius(shock_overdensity,js));
    //JRESET(pshock,tshock*(Md*darkmass(rshock,js)+gasmass(rshock,js))*
    //	   gasdensity(rshock,js)/(rshock*gasdensity_model2_slope(rshock,js)));
    
    JRESET(pshock,tshock*gasdensity(rshock,js)/(0.59*TMASS));
    //printf("%E %E %E %E\n",pshock,tshock,rshock,gasdensity(rshock,js));
    //JRESET(pshock,tshock*gasdensity(rshock,js));
    if (js->debug > 1)
      printf("truncation radius is %E\n",rshock);
  } else {
    JRESET(rshock,0.); 
    JRESET(pshock,0.);
  }
  
  // Output mass at various contrast if requested
  if (js->fitwrite > 1) {
    
    int k;
    k = INFO_SIZE*js->nannuli;
    for (i = 0; i < js->ndelta; ++i) {
      if (js->deltas[i] > 0)
	r200 = overdensity_radius(js->deltas[i],js);
      else 
	r200 = -js->deltas[i];
      double mtot = Md*darkmass(r200,js);
      double mg = gasmass(r200,js);
      double ms = (MsContr > 0 ? Ms*stellarmass(r200,js) : 0.);
      if (! js->totalmass)
	mtot+= mg+ms;
      //if (i == 0 && fabs(js->rContr-0.33) > 0.03)
      //printf("DIFFERENCE %E %E    %E %E\n",js->rContr,r200,js->MdContr,mtot);
      js->projtype = PROJ_ALL;
      js->infoarray[k++] = fabs(js->deltas[i]);
      js->infoarray[k++] = r200;
      js->infoarray[k++] = r200/rdm1;
      js->infoarray[k++] = mg;
      js->infoarray[k++] = mtot-mg-ms;
      js->infoarray[k++] = ms;
      js->infoarray[k++] = mtot;
      js->infoarray[k++] = FOURPI*projected_mass(0.,r200,js);
      /* snprintf(js->infoarray[js->nannuli+i],INFO_STRLEN, */
      /* 	       "%6.0f %E %E %E %E %E %E\n", */
      /* 	       deltas[i],r200,r200/rdm1,mg,md,ms,md+mg+ms); */
    }

  }

  // Splines for high-precision integration
  if (!dofast) {
    maxrad = (rshock > 0. ? rshock : overdensity_radius(200.,js));
    if (js->rspline) free(js->rspline);
    js->rspline = sci_dvector(nsplines);
    js->rspline[0] = 1.e-8;
    js->rspline[1] = 1.e-4;
    js->rspline[2] = 0.001;
    step = pow(maxrad/js->rspline[2],1./(nsplines-3));
    for (i = 3; i < nsplines; ++i)
      js->rspline[i] = js->rspline[2]*pow(step,i-2);
  }

  // Final normalization of output spectrum
  js->specnorm = js->opz*MPSQCM5/lumdistancesq;

#endif
  
}

// Initialize variables required for a single thread's execution
// These variables are not to be shared among threads
void jaco_thread_init(struct jaco_vars *jv) {


  jv->accel = gsl_interp_accel_alloc();
  jv->szaccel = gsl_interp_accel_alloc();
  jv->szyinteg = sci_dvector(2*NGAUSS);
  jv->xabs = sci_dvector(2*NGAUSS);
  jv->kernel = sci_dvector(2*NGAUSS);
  jv->tint = sci_dvector(2*NGAUSS);
  jv->logT = sci_dvector(2*NGAUSS);
  jv->logZ = sci_dvector(2*NGAUSS);
  jv->entropy = sci_dvector(2*NGAUSS);
  jv->calctproj = jv->integcounter = 0;
  
}

void jaco_thread_shutdown(struct jaco_vars *jv) {

  gsl_interp_accel_free(jv->accel);
  gsl_interp_accel_free(jv->szaccel);
  free(jv->szyinteg);
  free(jv->xabs);
  free(jv->kernel);
  free(jv->tint);
  free(jv->logT);
  free(jv->entropy);
  free(jv->logZ);

}

// Return X-ray photon spectrum in given annulus using hydrostatic
// equilibrium
void jaco_xray(int annulus, int nbin, double *spectrum, struct jaco_vars *jv)
{
  int   i;
  double R1,R2,Reff,coolingtime,elum,logZ=0.;
  double *w1spectrum=NULL,*w2spectrum=NULL,pressure=0.,logT;
  double tsum1,tsum2,normdensity=0.,tfullint=0.,temp,entropy;
  double shockfrac,mgas,mdark,mstar,mtot,mgascontr,mtotcontr;
  //double penalty=0.,polyindex=0.;

  struct jaco_state *js = jv->state;

  JS(angdist);  JS(r1list); JS(r2list); JS(effrad);
  JS(Md); JS(Ms); JS(MsContr);

  jv->R1deg = r1list[annulus]/60.;
  jv->R2deg = r2list[annulus]/60.;

  // Transform R1 and R2 from arcminutes to Mpc 
  // Angular diameter distance is dist/1+z 
  jv->R1 = R1 = r1list[annulus]*ARCMINTORAD*angdist;
  jv->R2 = R2 = r2list[annulus]*ARCMINTORAD*angdist;
  jv->R1sq = R1*R1;
  jv->R2sq = R2*R2;
  jv->Reff = Reff = effrad[annulus]*ARCMINTORAD*angdist;
  
  if (jv->R1 > js->rshock && js->rshock > 0.) {
    /* if (js->mpirank == 0) */
    /*   printf("Warning: annulus %d wholly outside truncation radius (%f Mpc)\n", */
    /* 	     annulus, */
    /* 	     js->rshock); */

    for (i = 0; i < nbin; ++i) 
      spectrum[i] = jv->R1-js->rshock;
    
    if (js->fitwrite > 1)
      js->infoarray[annulus*INFO_SIZE] = 1.E30;
    /* snprintf(js->infoarray[annulus],INFO_STRLEN, */
    /* 	     "%6.5f   --- Outside the truncation radius\n",jv->Reff); */
    return;

  } else if (jv->R2 > js->rshock && js->rshock > 0.) {
    R2 = jv->R2 = js->rshock;
    jv->R2sq = R2*R2;
    Reff = jv->Reff = 2.*(jv->R2sq*R2-jv->R1sq*R1)/
      (3.*(jv->R2sq-jv->R1sq));
    jv->R2deg = R2/(60.*ARCMINTORAD*angdist);
  }
  if (js->lastszannulus < annulus) js->lastszannulus = annulus;

  // Thus far, best results with spherical integration result when we
  // split the integration at R2
  jv->breakradius = R2;

  if (js->fitwrite > 2) {

    normdensity = gasdensity(jv->Reff,js);

    if (js->dofast)
      pressure = TMASS*
	(js->pshock+js->ninteg(temperature_intfunc,jv->Reff,
			       js->rshock,jv));
    else
      pressure = pow(10.,gsl_interp_eval(js->szinterp_coefs,js->rspline,js->pressure,jv->Reff,
			      jv->szaccel));
  
    tfullint = pressure/normdensity;

  // Introduce a penalty for a convectively unstable solution--inactive for now
/*     polyindex = -TMASS*mass(jv->Reff,js)/ */
/*     (jv->Reff*tfullint*gasdensity_model2_slope(jv->Reff,js)); */

/*     if (polyindex > 5/3.) { */
/*       //penalty = 3.*polyindex-5.; */
/*       if (js->debug != 0) printf("Polytopic index: %E\n",polyindex); */
/*     } else  */
/*       penalty = 0.; */

    w1spectrum = sci_dvector(nbin);
    w2spectrum = sci_dvector(nbin);
    jv->calctproj = 1;
    getspectrum(annulus,w1spectrum,nbin,jv);
    jv->calctproj = 2;
    getspectrum(annulus,w2spectrum,nbin,jv);
    jv->calctproj = 0;
  }

  if (js->fitwrite <= 1 || DIFFER(js->tprojpower,1.)) 
    getspectrum(annulus,spectrum,nbin,jv);
  else
    for (i = 0; i < nbin; ++i)
      spectrum[i] = w2spectrum[i];
    
  if (js->fitwrite > 2) {

    tsum1 = tsum2 = elum = 0.;
    logZ = log10(metallicity(Reff,jv));

    tfullint *= jv->mu;
    pressure *= jv->mu;

    logT = log10(tfullint);

    for (i = 0; (int)i < nbin; ++i) {
      tsum1 += w1spectrum[i];
      tsum2 += w2spectrum[i];
      jv->coolmatrix = jv->state->rebinnedcooling[i];
      elum += js->lastbin[i]*cooling(logT,logZ,jv);
    }

    elum *= (1+js->redshift);
    //double tshock = jv->mu*TMASS*js->pshock/normdensity;

    // Cooling time definition from  Soker, Blanton & Sarazin 2002
    shockfrac = jv->mu*TMASS*js->pshock/normdensity;
    normdensity *= jv->Hfrac*MPCM3;
    elum *= (1.e-14*jv->nenh*normdensity*normdensity);
    normdensity *= jv->nenh;
    pressure *= jv->nenh*jv->Hfrac*MPCM3;
    coolingtime = 1.5*pressure/elum;
    //printf("%E %E\n",coolingtime/GYR,3*13.7*(0.72/0.5)*sqrt(tfullint/8.6)/(normdensity/0.001));
    elum *= MPC*MPC*MPC*FOURPI*(R2*R2*R2-R1*R1*R1)/3.;
    temp = cbrt(normdensity);
    entropy = tfullint/(temp*temp);
    mgascontr = js->Mg1Contr+js->Mg2Contr+js->Mg3Contr;
    mgas = gasmass(Reff,js);

    if (MsContr > 0)
      mstar = Ms*js->stellarmass(Reff,js);
    else
      mstar = 0;

    mtot = mdark = Md*js->darkmass(Reff,js);
    if (js->totalmass)
      mdark -= mgas+mstar;
    else
      mtot = mgas+mdark+mstar;
    
    js->infoarray[annulus*INFO_SIZE] = Reff;
    js->infoarray[annulus*INFO_SIZE+1] = tsum1/tsum2;
    js->infoarray[annulus*INFO_SIZE+2] = mdark;
    js->infoarray[annulus*INFO_SIZE+3] = mgas;
    js->infoarray[annulus*INFO_SIZE+4] = mstar;
    js->infoarray[annulus*INFO_SIZE+5] = mgas/mtot;
    js->infoarray[annulus*INFO_SIZE+6] = tfullint;
    js->infoarray[annulus*INFO_SIZE+7] = normdensity;
    js->infoarray[annulus*INFO_SIZE+8] = pressure*KEV;
    js->infoarray[annulus*INFO_SIZE+9] = entropy;
    js->infoarray[annulus*INFO_SIZE+10] = coolingtime/GYR;
    js->infoarray[annulus*INFO_SIZE+11] = shockfrac;
    js->infoarray[annulus*INFO_SIZE+12] = elum*KEV/1.e42;
    js->infoarray[annulus*INFO_SIZE+13] = pow(10.,logZ);
    js->infoarray[annulus*INFO_SIZE+14] = mtot;
    
    /* snprintf(js->infoarray[annulus],INFO_STRLEN, */
    /* 	     "%6.5f %7.3f %10.3E %10.3E %10.3E %10.3E %10.3E  %10.3E "\ */
    /* 	     "%10.3E %10.3E %10.3E %10.3E %10.3E %10.3E\n", */
    /* 	     Reff,tsum1/tsum2, */
    /* 	     mdark,mgas,mstar,mgas/mtot,tfullint, */
    /* 	     normdensity,pressure*KEV,entropy, */
    /* 	     coolingtime/GYR,shockfrac, */
    /* 	     elum*KEV/1.e42,pow(10.,logZ)); */

    free(w1spectrum);
    free(w2spectrum);
  }
  
  for (i = 0; i < nbin; ++i)
    spectrum[i] *= //(1.+penalty*i*i*i*i)*
      js->specnorm;
      
  return;

}

void deallocatevectors()
{

  // Eventually, try to deallocate the vectors here.
}
