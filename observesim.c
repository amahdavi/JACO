#include <fitsio.h>
#include <sciutils.h>
#include <stdio.h>
#include <stdlib.h>
#include <time.h>
#include <string.h>
#include <math.h>
#include <unistd.h>
#include <jaco.h>
#include <gsl/gsl_sort.h>
#include <gsl/gsl_rng.h>
#include <gsl/gsl_randist.h>
#include <standalone.h>
#include <omp.h>


#define GRAVG  430157.2  // G in units of km^2 Mpc / 1e14 sunmass s^2
#define PDENSE  247.06557 // Protonmass/cm^3 in units of 1e14 sunmass/Mpc^3
#define CMPERMPC  3.08567E10 // 1e-14 * cm / Mpc

// Load a tipsy-format file into memory.
void read_sim(struct jaco_state *js)

{
  FILE *tf;
  int ndims;
  float simtime;
  long  i,j;
  int k;
  double velmean[3],nenh,Hfrac;


  JS(massunit); JS(mpc); JS(velunit); 

  tf = fopen(js->simfile,"r");

  fscanf(tf,"%ld %ld %ld",&js->simcount, &js->np[GAS], &js->np[STAR]);
  fscanf(tf,"%d",&ndims);
  if (ndims != 3) 
    BYE("Tipsy file must be 3 dimensional.\n");
  fscanf(tf,"%f",&simtime);
  js->np[DARK] = js->simcount-js->np[GAS]-js->np[STAR];
  printf("Dark, Gas, star parts: %ld %ld %ld\n",
	 js->np[DARK],js->np[GAS],js->np[STAR]);
  
  js->part = (struct particle *)
    malloc(js->simcount*sizeof(struct particle));

  js->ptype[GAS]  =  js->gas  = js->part;
  js->ptype[DARK] =  js->dark = &js->part[js->np[GAS]];
  js->ptype[STAR] =  js->star = &js->part[js->np[GAS]+js->np[DARK]];

  for (i = 0; i < js->simcount; ++i)  {
    fscanf(tf,"%f",&js->part[i].mass);
    js->part[i].mass *= massunit;
    if (i < js->np[GAS]) js->part[i].type=GAS;
    else if (i < js->np[GAS]+js->np[DARK]) js->part[i].type=DARK;
    else js->part[i].type=STAR;
  }
  
  for (i = 0; i < 3; ++i) 
    for (j = 0; j < 3; ++j) {
      js->boxmin[i][j] = 1.E30;
      js->boxmax[i][j] = -1.E30;
    }


  for (j = 0; j < 3; ++j) {
    for (i = 0; i < js->simcount; ++i)  {
      fscanf(tf,"%f",&js->part[i].pos[j]);
      js->part[i].pos[j] *= mpc;
      k = js->part[i].type;
      if (js->part[i].pos[j] < js->boxmin[k][j]) 
	js->boxmin[k][j] = js->part[i].pos[j];
      if (js->part[i].pos[j] > js->boxmax[k][j]) 
	js->boxmax[k][j] = js->part[i].pos[j];

    }
    velmean[j] = 0;
  }


  for (j = 0; j < 3; ++j) {
    for (i = 0; i < js->simcount; ++i)  {
      fscanf(tf,"%f",&js->part[i].vel[j]);
      js->part[i].vel[j] *= velunit;
      velmean[j] += js->part[i].vel[j];
    }
  }

  for (j = 0; j < 3; ++j) {
    velmean[j] /= js->simcount;
    for (i = 0; i < js->np[GAS]; ++i)  
      js->gas[i].opz[j] = 1.+(js->gas[i].vel[j]-velmean[j])/299792;
  }

  for (i = 0; i < js->np[DARK]; ++i)  fscanf(tf,"%f",&js->dark[i].eps);
  for (i = 0; i < js->np[STAR]; ++i)  fscanf(tf,"%f",&js->star[i].eps);

  for (i = 0; i < js->np[GAS]; ++i)  {
    fscanf(tf,"%f",&js->gas[i].rho);
    js->gas[i].rho *= massunit/(mpc*mpc*mpc);
  }

  for (i = 0; i < js->np[GAS]; ++i)  {
    fscanf(tf,"%f",&js->gas[i].temp);
    js->gas[i].temp /= 11604506;
    js->gas[i].logtemp = log10(js->gas[i].temp);
  }

  for (i = 0; i < js->np[GAS]; ++i) 
    fscanf(tf,"%f",&js->gas[i].hsmooth);

  for (i = 0; i < js->np[GAS]; ++i) {
    fscanf(tf,"%f",&js->gas[i].metal);
    js->gas[i].metal /= 0.0189;
    js->gas[i].metal=0.3;
    nenh = 1.170201+0.011532*js->gas[i].metal;
    Hfrac = 0.740442-0.012406*js->gas[i].metal;
    js->gas[i].logmetal = log10(js->gas[i].metal);
    js->gas[i].nenp = Hfrac*Hfrac*nenh*js->gas[i].rho*js->gas[i].rho
      *pow(js->gas[i].hsmooth*mpc,3.);

  }

  for (i = 0; i < js->np[STAR]; ++i)  fscanf(tf,"%f",&js->star[i].tform);

  for (i = 0; i < js->simcount; ++i)  fscanf(tf,"%f",&js->part[i].phi);

  printf("Done reading.\n");

  fclose(tf);
}

// Calculate the center of the densest peak on the map in O(N) time.
// Do this by dividing the whole area into a 10x10x10 grid. Find the 
// single densest cell.  Create a 10x10x10 grid centered on the
// densest cell, but with 25% the dimensions of the previous grid.
// Continue until no cell in a new grid is denser than the densest
// cell in the old grid.
void calc_centroids(int t, struct jaco_state *js)
{
  int cs = 10;
  double mingrid[3],maxgrid[3],density=1.,lastdensity=0.;
  float del[3];
  long i,nmax,gp[3],cell[cs][cs][cs]; 
  int k,l,m;
  struct particle *p;
  
  for (k = 0; k < 3; ++k) mingrid[k] = js->boxmin[t][k];
  for (k = 0; k < 3; ++k) maxgrid[k] = js->boxmax[t][k];

  p = js->ptype[t];

  while (density > lastdensity) { 
  
    for (k = 0; k < 3; ++k) 
      del[k] = (maxgrid[k]-mingrid[k])/cs;
    
    nmax = -1;
    for (k = 0; k < cs; ++k)
      for (l = 0; l < cs; ++l)
	for (m = 0; m < cs; ++m)
	  cell[k][l][m] = 0;

    for (i = 0; i < js->np[t]; ++i) {
      if (p[i].pos[0] > mingrid[0] && p[i].pos[0] < maxgrid[0] &&
	  p[i].pos[1] > mingrid[1] && p[i].pos[1] < maxgrid[1] &&
	  p[i].pos[2] > mingrid[2] && p[i].pos[2] < maxgrid[2]) {
	for (k = 0; k < 3; ++k) 
	  gp[k] = (long)floorl((p[i].pos[k]-mingrid[k])/del[k]);
	++cell[gp[0]][gp[1]][gp[2]];
	if (cell[gp[0]][gp[1]][gp[2]] > nmax) {
	  for (k = 0; k < 3; ++k) 
	    js->center[t][k] = mingrid[k]+del[k]*gp[k];
	  nmax = cell[gp[0]][gp[1]][gp[2]];
	}
      }
    }

    lastdensity = density;
    density = log(nmax)-log(del[0])-log(del[1])-log(del[2]);
    
    for (k = 0; k < 3; ++k) {
      mingrid[k] = js->center[t][k]-cs*del[k]/4.;
      maxgrid[k] = js->center[t][k]+cs*del[k]/4.;
      //printf("%f ",center[k]/boxmax[k]);
    }
    //printf(" %ld %f %f\n",nmax,lastdensity,density);

  }

  for (i = 0; i < 3; ++i)
    printf("%ld: %.2E (%.2E to %.3E)\n",
	   i,js->center[t][i],js->boxmin[t][i],js->boxmax[t][i]);
  printf("\n");

}

// Calculated projected and 3D distances of the parts to the
// center.
void calc_distances(long n, struct particle *p, double *center, 
		    int axis, double **distsq)
{
  long i;
  double diffsq[3];
  int j,k1,k2,k3;

  k1 = axis; k2 = (axis+1)%3; k3 = (axis+2)%3;

  #pragma omp parallel for private(i,j,diffsq) schedule(dynamic,500)
  for (i = 0; i < n; ++i) {
    for (j = 0; j < 3; ++j) {
      diffsq[j] = p[i].pos[j]-center[j];
      diffsq[j] *= diffsq[j];
    }
    distsq[1][i] = diffsq[k2]+diffsq[k3];
    distsq[0][i] = distsq[1][i]+diffsq[k1];
  }

}

void make_simimages(struct jaco_state *js)
{

  struct particle *p;
  double mpctopix;
  float **ximage[3];
  long xcen,ycen,i,x,y;
  int k1, k2, k3, t;
  char of[200];

  k1 = js->axis; k2 = (js->axis+1)%3; k3 = (js->axis+2)%3;

  ximage[0] = readimage(js->ximage,&js->pixscale,&js->nx,&js->ny,NULL);
  sci_free_fmatrix(ximage[0]);

  mpctopix = 60./(ARCMINTORAD*js->angdist*js->pixscale);
  xcen = js->nx/2;
  ycen = js->ny/2;

  // Make a mass map of each particle type
  # pragma omp parallel for private(t,i,p,x,y) schedule(dynamic,1)
  for (t = 0; t < 3; ++t) {

    ximage[t] = sci_fmatrix(js->nx,js->ny);
    p = js->ptype[t];

    for (i = 0; i < js->np[t]; ++i) {

      x = xcen+(p[i].pos[k2]-js->center[GAS][k2])*mpctopix;
      y = ycen+(p[i].pos[k3]-js->center[GAS][k3])*mpctopix;

      //printf("%ld %ld %E %E\n",x,y,mpctopix,js->center[GAS][k2]);
      if (x >= 0 && x < (long)js->nx && y >= 0 && y < (long)js->ny)
	ximage[t][x][y] += p[i].mass;
    
    }

  }

  for (t = 0; t < 3; ++t)  {
    snprintf(of,199,"%s.sim%s.fits",js->profname,
	     (t == DARK ? "dark" : (t == GAS ? "gas" : "star")));
    unlink(of);
    writeimage(of,ximage[t],js->nx,js->ny,js->ximage);
    sci_free_fmatrix(ximage[t]);
  }
}

void make_xrayimage(struct jaco_state *js)
{

  struct particle *p;
  double mpctopix;
  float **ximage;
  long xcen,ycen,i,x,y;
  int k1, k2, k3;
  char of[200];

  struct fitdata *data;

  data = (struct fitdata *)js->data;

  k1 = js->axis; k2 = (js->axis+1)%3; k3 = (js->axis+2)%3;

  mpctopix = 60./(ARCMINTORAD*js->angdist*js->pixscale);
  xcen = js->nx/2;
  ycen = js->ny/2;

  ximage = sci_fmatrix(js->nx,js->ny);
  p = js->ptype[GAS];

  for (i = 0; i < js->np[GAS]; ++i) {

    x = xcen+(p[i].pos[k2]-js->center[GAS][k2])*mpctopix;
    y = ycen+(p[i].pos[k3]-js->center[GAS][k3])*mpctopix;

    if (x >= 0 && x < (long)js->nx && y >= 0 && y < (long)js->ny)
      ximage[x][y] += p[i].spec*js->specnorm*
	data->phalist[0].exposure;
    
  }

  printf("EEXP %E\n",data->phalist[0].exposure);
  
  #pragma omp parallel for private(x,y)
  for (x = 0; x < (long)js->nx; ++x)
    for (y = 0; y < (long)js->ny; ++y) 
      ximage[x][y] = gsl_ran_poisson(js->rng,10000.*ximage[x][y]+2.);

  snprintf(of,199,"%s.simxray.fits",js->profname);
  unlink(of);
  writeimage(of,ximage,js->nx,js->ny,js->ximage);
  sci_free_fmatrix(ximage);
  printf("OK\n");

}

void simulate_spec(int annulus, int nbin, double *spectrum, 
		     struct jaco_vars *jv)
{

  
  int t;
  long p1, p2, i, j, k, l, peff;
  struct particle p, *ppoint;
  double unweightnorm=0.,weightnorm3D=0.,weightnormproj=0.;
  double tweight3D=0.,tweightproj=0.,mweight=0.,dist3dsq;
  double t1,t2,tmass[3];

  gsl_interp_accel *ene_accel = gsl_interp_accel_alloc();
  gsl_interp* ene_interp = gsl_interp_alloc(gsl_interp_linear,nbin);

  struct jaco_state *js = jv->state;

  struct fitdata *data;

  data = (struct fitdata *)js->data;

  JS(r1list); JS(r2list); JS(angdist); JS(axis); JS(effrad);

  double R1 = r1list[annulus]*ARCMINTORAD*angdist;
  double R2 = r2list[annulus]*ARCMINTORAD*angdist;
  double Reff = effrad[annulus]*ARCMINTORAD*angdist;

  double R1sq = R1*R1;
  double R2sq = R2*R2;

  p1 = gsl_interp_bsearch(js->sortdist2D,R1sq,0,js->simcount-1);

  if (annulus == js->nannuli-1 && js->dosim > 3)
    p2 = gsl_interp_bsearch(js->sortdist2D,36.*R2sq,0,js->simcount-1);
  else 
    p2 = gsl_interp_bsearch(js->sortdist2D,R2sq,0,js->simcount-1);
  peff = gsl_interp_bsearch(js->sortdist3D,Reff*Reff,0,js->simcount-1);

  double *ene = sci_dvector(nbin);
  double *restspec = sci_dvector(nbin);
  double *redspec = sci_dvector(nbin);

  tmass[GAS] = tmass[DARK] = tmass[STAR] = 0;

  for (j = 0; j < nbin; ++j) {
    spectrum[j] = 0.;
    ene[j] = js->opz*js->lastbin[j];
  }
  
  for (i = p1; i <= p2; ++i) {

    p = js->part[js->order2D[i]];
    ppoint = &js->part[js->order2D[i]];
    dist3dsq = js->distsq[0][js->order2D[i]];

    t = p.type;
    ppoint->spec = 0.;
    
    if (dist3dsq < R2sq) 
      tmass[p.type] += p.mass;

    if (t == GAS) {

      if (p.nenp < 1.E90) {

	for (j = 0; j < nbin; ++j) {
	  jv->coolmatrix = js->rebinnedcooling[j];
	  restspec[j] = p.nenp*cooling(p.logtemp,p.logmetal,jv);
	}

	gsl_interp_init(ene_interp,ene,restspec,nbin);

	for (j = 0; j < nbin; ++j) {
	  redspec[j] = 
	    gsl_interp_eval(ene_interp,ene,restspec,ene[j]*p.opz[axis],
			    ene_accel);
	  spectrum[j] += redspec[j];

	  if (js->fitwrite > 2) {
	    t1 = pow(p.temp,js->tprojpower)*redspec[j];
	    t2 = t1/p.temp;

	    tweightproj += t1;
	    weightnormproj += t2;

	    if (dist3dsq < R2sq) {

	      tweight3D += t1;
	      weightnorm3D += t2;

	      mweight += p.metal*redspec[j];
	      unweightnorm += redspec[j];
	    }
	  }

	  if (js->dosim > 2) 
	    for (k = 0; k < data->rmflist[annulus].ngrp[j]; ++k) 
	      for (l = 0; l < data->rmflist[annulus].nchannels[j][k]; ++l) 
		ppoint->spec += data->rmflist[annulus].response[j][k][l]*redspec[j]*(data->arflist[annulus].ehi[j]+data->arflist[annulus].elo[j])/2.;
	}
	
      }
    }
  }

  if (js->fitwrite > 2) {

    t1 = 4.1887902*(R2*R2sq-R1*R1sq);
    for (j = 0; j < 3; ++j)
      tmass[j] /= t1;

    // Fix this for observesim mode
    /* snprintf(js->infoarray[annulus],INFO_STRLEN, */
    /* 	     "%6.5f %10.3E %10.3E %10.3E %10.3E %10.3E $10.3E %10.3E "\ */
    /*           "%10.3E %10.3E %10.3E %10.3E\n", */
    /* 	     Reff,tweightproj/weightnormproj, */
    /* 	     tmass[DARK],tmass[GAS],tmass[STAR], */
    /* 	     js->simmofr[DARK][peff], */
    /* 	     js->simmofr[GAS][peff],js->simmofr[STAR][peff],  */
    /* 	     js->simmofr[GAS][peff]/(js->simmofr[DARK][peff]+  */
    /* 				    js->simmofr[GAS][peff]+ */
    /* 				    js->simmofr[STAR][peff]), */
    /* 	     tweight3D/weightnorm3D, */
    /* 	     mweight/unweightnorm); */

  }
  
  gsl_interp_accel_free(ene_accel);
  gsl_interp_free(ene_interp);

  for (j = 0; j < nbin; ++j)
    spectrum[j] *= js->specnorm;

  free(ene);
  free(restspec);
  free(redspec);

}
	      

int init_sim(struct jaco_state *js)
{
  int j;
  long k;
  double tmass[3];
  struct particle p;

  JS(axis);
  
  js->velunit = sqrt(GRAVG*js->massunit/js->mpc);
  read_sim(js);

  js->center = sci_dmatrix(3,3);

  calc_centroids(GAS,js);
  calc_centroids(DARK,js);
  calc_centroids(STAR,js);
  
  js->distsq = sci_dmatrix(2,js->simcount);

  calc_distances(js->simcount,js->part,js->center[GAS],axis,js->distsq);

  js->order2D = sci_sizetvector(js->simcount);
  js->order3D = sci_sizetvector(js->simcount);
  js->sortdist2D = sci_dvector(js->simcount);
  js->sortdist3D = sci_dvector(js->simcount);
  js->simmofr = sci_dmatrix(3,js->simcount);

  gsl_sort_index(js->order3D,js->distsq[0],1,js->simcount);
  gsl_sort_index(js->order2D,js->distsq[1],1,js->simcount);

  tmass[GAS] = tmass[DARK] = tmass[STAR] = 0;

  for (k = 0; k < js->simcount; ++k)  {
    js->sortdist2D[k] = js->distsq[1][js->order2D[k]];
    js->sortdist3D[k] = js->distsq[0][js->order3D[k]];

    // Calculate mass profile for different particle types
    p = js->part[js->order3D[k]];

    tmass[p.type] += p.mass;

    for (j = 0; j < 3; ++j)
      js->simmofr[j][k] = tmass[j];

  }

  return 0;

}

/* void rotateview (struct particle *p, float *center,  */
/* 		 long npart, float degr, int axis) */
/* { */

/*   long i; */
/*   int xaxis,yaxis,j,k; */
/*   float x,y; */
/*   float costheta = cos(degr*PI/180.); */
/*   float sintheta = sin(degr*PI/180.); */

/*   xaxis = yaxis = -1; */
/*   for (j = 0; j < 3; ++j)  */
/*     if (xaxis < 0 && j != axis) xaxis = j; */
/*     else if (yaxis < 0 && j != axis) yaxis = j; */

/*   for (i = 0; i < npart; ++i) { */
/*     x = p[i].pos[xaxis]-center[xaxis]; */
/*     y = p[i].pos[yaxis]-center[yaxis]; */
/*     p[i].pos[xaxis] = x*costheta+y*sintheta+center[xaxis]; */
/*     p[i].pos[yaxis] = y*costheta-x*sintheta+center[yaxis]; */
/*   } */
/* } */

/* int main(int argc, char *argv[]) */
/* { */
/*   char outf[400]; */
/*   double *ebins,z,edelta,e1,e2; */
/*   float *eax,*eay, effectivearea, esum, pixsize, exptime;; */
/*   unsigned long eacount; */
/*   int i,j,k; */
/*   long l; */

/*   read_units(argv[1]); */
/*   read_tipsy(argv[8]); */

/*   eax = float_readdata("earea.dat",1,&eacount); */
/*   eay = float_readdata("earea.dat",2,&eacount); */

/*   initvars("/home/amahdavi/xray/models/mekal.txt"); */
/*   //calc_centroids(np[GAS],gas,&mass[GAS],gascenter); */
/*   for (k = 0; k < 3; ++k) */
/*     gascenter[k] = 0; */
/*   for (l = 0; l < np[GAS]; ++l)  */
/*     for (k = 0; k < 3; ++k) */
/*       gascenter[k] += gas[l].pos[k]; */
/*   for (k = 0; k < 3; ++k) */
/*     gascenter[k] /= np[GAS]; */

/*   z0 = 0.3; zinf =  0.3; rz = 1.; z= atof(argv[2]); */
/*   opz = 1.+z; */
/*   ebins = dvector(1,2); */
/*   ebins[1] = 3.00; */
/*   ebins[2] = 8.00; */
/*   edelta = 0.75; */
/*   pixsize = atof(argv[3]); */
/*   exptime = atof(argv[4]); */

/*   debug = 8; */
/*   rebincoolingfunction(ebins,2,atof(argv[5])); */

/*   for (i = 1; i <= 1; ++i) { */

/*     e1 = ebins[i]-edelta; */
/*     e2 = ebins[i]+edelta; */

/*     background = 5.*part(e1,e2)*pixsize*pixsize*exptime; */

/*     effectivearea=esum=0.; */
/*     for (j = 1; j <= eacount; ++j)  */
/*       if (eax[j]/100. >= e1 && eax[j]/100. <= e2) { */
/* 	effectivearea += eay[j]/eax[j]; */
/* 	esum += 1./eax[j]; */
/*       } */
    
/*     effectivearea /= esum; */
/*     printf("For %f-%f weighted area is %f cm^2\n",e1,e2,effectivearea); */
/*     coolmatrix = rebinnedcooling[i]; */

/*     for (j = 8; j < argc; ++j) { */

/*       if (j > 8) { */
/* 	free(gas); free(dark); free(star); */
/* 	read_tipsy(argv[j]); */
/*       } */

/*       sprintf(outf,"%sx-%d-%.2f-%04.0f.fits",argv[j],atoi(argv[6]),atof(argv[7]),ebins[i]*1000); */
/*       makeobs(0,0,z,atof(argv[3]),atof(argv[4]),gascenter,effectivearea, */
/* 	      atoi(argv[6]),atof(argv[7]),outf); */

/* /\*     if (i <= 360) *\/ */
/* /\*       rotateview(gas,gascenter,np[GAS],0.5,1); *\/ */
/* /\*     else *\/ */
/* /\*       rotateview(gas,gascenter,np[GAS],0.5,2); *\/ */
/*       sprintf(outf,"%sy-%d-%.2f-%04.0f.fits",argv[j],atoi(argv[6]),atof(argv[7]),ebins[i]*1000); */
/*       makeobs(1,0,z,atof(argv[3]),atof(argv[4]),gascenter,effectivearea, */
/* 	      atoi(argv[6]),atof(argv[7]),outf); */
/*       sprintf(outf,"%sz-%d-%.2f-%04.0f.fits",argv[j],atoi(argv[6]),atof(argv[7]),ebins[i]*1000); */
/*       makeobs(2,0,z,atof(argv[3]),atof(argv[4]),gascenter,effectivearea, */
/* 	      atoi(argv[6]),atof(argv[7]),outf); */
/*     } */
/*   } */
  
/*   return 0.; */
/* } */

/* // Write various profiles: density, temperature, velocity, anisotropy etc. */
/* void write_profile(long n, struct particle *p, float *dist,  */
/* 		  long nmin, int dim, float *center, char *outf, */
/* 		  char *metalfile, char *tempfile, char *anisofile) */
/* { */
/*   FILE *of,*mf=NULL,*tf=NULL,*af=NULL; */
/*   float *ord,*velr,*veltheta,*velphi; */
/*   long i; */
/*   long nannulus,nbin; */
/*   float mbin,d1,d2,volume,dd1,dd2,meand,densq,velsq,meanvr,meanvt,sigmar,sigmat; */
/*   float diff,pos[3],vel[3],meanvp,sigmap; */
/*   double meant=0.,meanm=0.,rhoweight = 0.; */
/*   struct particle currp; */
/*   int k; */
 

/*   of = fopen(outf,"w+"); */
/*   if (metalfile != NULL) mf = fopen(metalfile,"w+"); */
/*   if (tempfile != NULL) tf = fopen(tempfile,"w+"); */
/*   if (anisofile != NULL) af = fopen(anisofile,"w+"); */
/*   ord = (float *)malloc(n*sizeof(float)); */
/*   for (i = 0; i < n; ++i) { ord[i] = (float)i;  } */

/*   sort2((unsigned long)n,dist-1,ord-1); */

/*   mbin = 0.; */
/*   nannulus = nbin = 0; i = -1; */
/*   d1 = dist[0]; */

/*   velr = vector(1,nmin); */
/*   veltheta = vector(1,nmin); */
/*   velphi = vector(1,nmin); */

/*   while (nbin < nmin && i < n) { */
/*     ++i; */
/*     ++nbin; */
/*     currp = p[(long)ord[i]]; */
/*     mbin += currp.mass; */

/*     // Emission-weight the gas. */
/*     if (metalfile != NULL || tempfile != NULL)  */
/*       if (currp.temp > 0.3 && currp.rho < 0.001) { */
/* 	densq = currp.rho; */
/* 	densq *= densq*sqrt(currp.temp); */
/* 	meant += densq*currp.temp; */
/* 	meanm += densq*currp.metal; */
/* 	rhoweight += densq; */
/*       } */

/*     if (anisofile != NULL && dist[i] > 0.) { */
/*       velr[nbin] = velsq = 0.; */
/*       meand = sqrt(dist[i]); */
/*       for (k = 0; k < 3; ++k) { */
/* 	vel[k] = currp.vel[k]; */
/* 	//currp.vel[k]=1.-2.*ran1(&idum); */
/* 	pos[k] = currp.pos[k]-center[k]; */
/* 	// Radial component of the velocity */
/* 	velr[nbin] += vel[k]*pos[k]/meand; */
/* 	velsq += vel[k]*vel[k]; */
/*       } */
/*       diff = sqrt(pos[0]*pos[0]+pos[1]*pos[1]); */
/*       // Tangential components of the velocity. */
/*       velphi[nbin] =  */
/* 	vel[0]*pos[0]*pos[2]/(meand*diff)+ */
/* 	vel[1]*pos[1]*pos[2]/(meand*diff)- */
/* 	vel[2]*diff/meand; */
/*       veltheta[nbin] = (vel[1]*pos[0]-vel[0]*pos[1])/diff; */
/*     } */

/*     for (k = 0; k < 3; ++k) */
/*       if (k != dim) */
/* 	if (center[k]-sqrt(dist[i]) < boxmin[k] || */
/* 	    center[k]+sqrt(dist[i]) > boxmax[k]) { */
/* 	  printf("----Quitting at %f %f %ld  (%f%%)\n", */
/* 	       p[(long)ord[i]].pos[k],sqrt(dist[i]),i,100.*i/n); */
/* 	return;  */
/*       } */
/*     if (nbin == nmin || i == n-1) { */
/*       ++nannulus; */
/*       //printf("%s %d %ld \n",outf,dim,nannulus); */
/*       d2 = dist[i]; */
/*       dd1 = sqrt(d1); */
/*       dd2 = sqrt(d2); */
/*       if (dim > -1) volume = PI*(d2-d1); */
/*       else volume = 4.*PI*(d2*dd2-d1*dd1)/3.; */
/*       meand = (dd1+dd2)/2.; */
/*       fprintf(of,"%10E %10E\n",meand,mbin/volume); */
/*       if (metalfile != NULL)  */
/* 	fprintf(mf,"%10E %10E\n",meand,meanm/rhoweight); */
/*       if (tempfile != NULL)  */
/* 	fprintf(tf,"%10E %10E\n",meand,meant/rhoweight); */
/*       meanvr = meanvt = meanvp = sigmar = sigmap = sigmat = 0.; */
/*       if (anisofile != NULL) { */
/* 	for (k = 1; k <= nbin; ++k) { */
/* 	  meanvr += velr[k]; */
/* 	  meanvt += veltheta[k]; */
/* 	  meanvp += velphi[k]; */
/* 	} */
/* 	meanvr /= nbin; */
/* 	meanvt /= nbin; */
/* 	meanvp /= nbin; */
/* 	//if (nannulus == 1) {  */
/* 	//  cvelr = meanvr; cvelt = meanvt; cvelp = meanvp;  */
/* 	//} */
/* 	//meanvr -= cvelr; */
/* 	//meanvt -= cvelt; */
/* 	//meanvp -= cvelp; */
/* 	for (k = 1; k <= nbin; ++k) { */
/* 	  diff = velr[k]-meanvr; */
/* 	  sigmar += diff*diff; */
/* 	  diff = veltheta[k]-meanvt; */
/* 	  sigmat += diff*diff; */
/* 	  diff = velphi[k]-meanvp; */
/* 	  sigmap += diff*diff; */
/* 	} */
/* 	fprintf(af,"%10E %10E %10E %10E %10E \n", */
/* 		meand,meanvr,meanvt,meanvp,1.-(sigmat+sigmap)/(2.*sigmar)); */
/*       } */
/*       fflush(of); */
/*       meant = meanm = rhoweight = 0.; */
/*       d1 = d2; */
/*       nbin = 0; */
/*       mbin = 0.; */
/*     } */
/*   } */

/*   free(ord); */
/*   fclose(of); */
/*   if (metalfile != NULL) fclose(mf); */
/*   if (tempfile != NULL) fclose(tf); */
/*   if (anisofile != NULL) fclose(af); */
/* } */

/* void makeimage(int axis, int instrument, float z, float pixsize, */
/* 	       long exptime, float *center, float effectivearea,  */
/* 	       int dofullobs, float tpow, */
/* 	       char *outfile) */
/* { */
/*   fitsfile *fout; */
/*   long nx, ny, i, xp, yp, zp, xx, yy, naxes[2], fpixel=1, npix, x, y; */
/*   float *buff,val,comoving,lumdist,lumfac,sphsize; */
/*   float pixpermpc,hubble=70.,omega_m=0.3,pvol,psize,psizesq,crad,cbeta; */
/*   float fieldsizex=30.,fieldsizey=30.,*emission,diff,ysize,ysizesq; */
/*   float axisdist,spectrum; */
/*   int xaxis,yaxis,j,status=0,nchar; */
/*   long percent; */


/*   naxes[0] = nx = (long)(60.*fieldsizex/pixsize); */
/*   naxes[1] = ny = (long)(60.*fieldsizey/pixsize); */
/*   npix = nx*ny; */
/*   buff = (float *)malloc(npix*sizeof(float)); */

/*   emission = sci_fvector(1,np[GAS]); */

/*   xaxis = yaxis = -1; */
/*   for (j = 0; j < 3; ++j)  */
/*     if (xaxis < 0 && j != axis) xaxis = j; */
/*     else if (yaxis < 0 && j != axis) yaxis = j; */

/*   comoving = ztoMpc(hubble,omega_m,-1.,z); */
/*   pixpermpc = 60./(pixsize*ARCMINTORAD*comoving/(1.+z)); */

/*   for (i = 0; i < np[GAS]; ++i) { */
/*     diff = gas[i].rho/gas[i].mass; */
/*     if (gas[i].temp > 0.05) */
/*       emission[i] = diff*diff*cooling(gas[i].logtemp,gas[i].logmetal); */
/*     else */
/*       emission[i] = 0; */

/*     pressure[i] = diff*gas[i].temp;     */
/*   } */

/*   percent = (long)np[GAS]/50; */
/*   for (i = 0; i < np[GAS]; ++i) { */
/*     // Only consider parts within a cubic box  */
/*     zp = (long)floor(pixpermpc*(gas[i].pos[axis]-center[axis])+nx/2); */
/*     if (gas[i].temp > 0.05 && zp >= 0 && zp < nx) { */
/*       //sphsize = cbrt(32.*3.*gas[i].mass/(FOURPI*PDENSE*gas[i].rho)); */
/*       //if (gas[i].rho > 0.005) gas[i].rho = 0.005; */
/*       xp = (long)floor(pixpermpc*(gas[i].pos[xaxis]-center[xaxis])+nx/2); */
/*       yp = (long)floor(pixpermpc*(gas[i].pos[yaxis]-center[yaxis])+ny/2); */
/*       spectrum = emission[i]; */
/*       psize = pixpermpc*sphsize; */
/*       psizesq = psize*psize; */
/*       if (i % percent == 0) { */
/* 	nchar = printf("%.0f%% complete.",100.*(float)i/(np[GAS]-1.)); */
/* 	for (j = 0; j < nchar; ++j) printf("\b"); */
/* 	fflush(stdout); */
/*       } */
/*       spectrum /= 0.312561*psize*psizesq; */

/*       spectrum =  */

/*       for (x = -(long)psize; x <= (long)psize; ++x) { */
/* 	ysizesq = psizesq-x*x; */
/* 	ysize = sqrt(ysizesq < 0 ? 0. : ysizesq); */
/* 	for (y = -(long)ysize; y <= (long)ysize; ++y) { */
/* 	  xx = xp+x; yy = yp+y; */
/* 	  if (xx >= 0 && xx < nx && yy >=0 && yy < ny) { */
/* 	    usq = (x*x+y*y)/psizesq; */

/* 	    if (usq < 1.) { */
/* 	      u = sqrt(usq); */
/* 	      omu = 1.-u; */
	      
/* 	      if (u < 0.5) */
/* 		kern = 1.818913635-u*u*omu*10.91348181; */
/* 	      else */
/* 		kern = 3.637827271*omu*omu*omu; */

/* 	      buff[yy*nx+xx] += spectrum*kern; */
/* 	    } */
/* 	  } */

/* 	} */
/*       } */
/*     } */
/*   } */

/*   free_vector(emission,1,np[GAS]); */
/*   // Add soft X-ray background; */

/*   lumdist = comoving*(1.+z); */
/*   lumfac = effectivearea*exptime*CMPERMPC/(lumdist*lumdist); */
/*   printf("%f\n",background); */
/*   if (dofullobs) { */
/*     spectrum = 1.e-12*cooling(0.25,0.); */
/*     for (x = 0; x < nx; ++x) */
/*       for (y = 0; y < ny; ++y) { */
/* 	diff = x-nx/2; */
/* 	axisdist = diff*diff; */
/* 	diff = y-ny/2; */
/* 	axisdist += diff*diff; */
/* 	axisdist = sqrt(axisdist)*pixsize; */
/* 	diff = axisdist/430; */
/* 	buff[y*nx+x] += (exp(-diff*diff/2)+0.54)*spectrum; */
/*       } */
/*     buff[i] = lumfac*buff[i]+background*(1.+0.1*gasdev(&idum)); */
/*     for (i = 0; i < npix; ++i)  */
/*       buff[i] = poidev(buff[i],&idum); */
/*   } else { */
/*     for (i = 0; i < npix; ++i)  */
/*       buff[i] = lumfac*buff[i]; */
/*   } */

/*   unlink(outfile); */
/*   if(fits_create_file(&fout,outfile,&status)) { */
/*     return; */
/*   } */

/*   //nfaint = 206*pow(0.1,-0.793)*(fieldsizex*fieldsizex/3600.); */
  
      
/*   fits_create_img(fout,-32,2,naxes,&status); */
/*   fits_write_img(fout,TFLOAT,fpixel,npix,buff,&status); */
/*   free(buff); */

/*   fits_write_key(fout,TSTRING,"WCSNAMEP","PHYSICAL","Parameter Coord. System", */
/* 		      &status); */
/*   fits_write_key(fout,TSTRING,"WCSTY1P","PHYSICAL","Parameter Coord. System", */
/* 		      &status); */
/*   fits_write_key(fout,TSTRING,"CTYPE1P","x","Parameter 1",&status); */
/*   val = 0.; */
/*   fits_write_key(fout,TFLOAT,"CRVAL1P",&val,"Parameter 1",&status); */
/*   val = nx/2.; */
/*   fits_write_key(fout,TFLOAT,"CRPIX1P",&val,"Parameter 1",&status); */
/*   val = 1./pixpermpc; */
/*   fits_write_key(fout,TFLOAT,"CDELT1P",&val,"Parameter 1",&status); */
/*   fits_write_key(fout,TSTRING,"WCSTY2P","PHYSICAL","Parameter Coord. System", */
/* 		      &status); */
/*   fits_write_key(fout,TSTRING,"CTYPE2P","y","Parameter 2",&status); */
/*   val = 0.;    */
/*   fits_write_key(fout,TFLOAT,"CRVAL2P",&val,"Parameter 2",&status); */
/*   val = ny/2.;    */
/*   fits_write_key(fout,TFLOAT,"CRPIX2P",&val,"Parameter 2",&status); */
/*   val = 1./pixpermpc;    */
/*   fits_write_key(fout,TFLOAT,"CDELT2P",&val,"Parameter 2",&status); */

/*   fits_close_file(fout,&status); */


/* } */
