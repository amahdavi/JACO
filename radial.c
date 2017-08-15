#include <fitsio.h>
#include <math.h>
#include <stdio.h>
#include <unistd.h>
#include <stdlib.h>
#include <strings.h>
#include <sciutils.h>
#include <gsl/gsl_sort.h>
#include <gsl/gsl_interp.h>

#define EPS 1.e-6

void reject(size_t *id, size_t *xxn, size_t *yyn, float **image,
	    unsigned long n1, unsigned long n2, int *rej, 
	    unsigned long *nused, float *counts, float **exposure,
	    float nsigma, float maxcounts )
{

  unsigned long p, nrej, ngood=100, i;
  float *values, mean, sigma, diff, ct;

  nrej = 1;

  values = sci_fvector(n2-n1+1);
  ct = 0.;

  while (nrej > 0 && ngood > 10) {

    ct = 0.;
    nrej = 0;
    ngood = 0;
    mean = 0.;
    unsigned long ndetect = 0;
    for (i = n1; i <= n2; ++i) {
      p = (unsigned long)id[i];
      if (exposure != NULL)
	if (!rej[p] && exposure[xxn[p]][yyn[p]] <= 0) {
	  ++nrej;
	  rej[p] = 1;
	}
      if (!rej[p] && 
	  (image[xxn[p]][yyn[p]] > maxcounts || 
	   image[xxn[p]][yyn[p]] < -20.)) {
	++nrej;
	rej[p] = 1;
      }
      if (!rej[p]) {
	//printf("************* %ld %ld %f\n",xxn[p],yyn[p],image[xxn[p]][yyn[p]]);
	values[ngood] = image[xxn[p]][yyn[p]];
	if (exposure != NULL) 
	  values[ngood] /= exposure[xxn[p]][yyn[p]];
	if (fabs(image[xxn[p]][yyn[p]]) > 1e-10) ++ndetect;
	mean += values[ngood];
	++ngood;
      }
    }

    if (ndetect > 3) {
      mean /= ndetect;
      
      sigma = 0.;
      for (i = 0; i <= ngood-1; ++i) 
	if (fabs(values[i]) > 1e-10) {
	  diff = values[i]-mean;
	  sigma += diff*diff;
	}
      
      sigma = sqrt(sigma/(ndetect-1));
      
      // sort(ngood,values);
      // median = (ngood % 2 ? values[(ngood+1)/2] : 
      //	            (values[ngood/2]+values[ngood/2-1])/2.);
      
      //printf("------%f %f %ld %ld %ld\n",mean,sigma,ndetect,nrej, ngood);
      for (i = n1; i <= n2; ++i) {
	p = id[i];
	if (!rej[p]) {
	  diff = image[xxn[p]][yyn[p]];
	  if (fabs(diff) > 1.e-4) {
	    if (exposure != NULL) diff /= exposure[xxn[p]][yyn[p]];
	    if ((fabs(diff-mean) > nsigma*sigma)) {
	      rej[p] = 1;
	      ++nrej;
	    } else
	      ct += image[xxn[p]][yyn[p]];
	  }
	}
      }
      //printf("------%f %f %ld %ld %ld\n",mean,sigma,ndetect,nrej, ngood);
    }
    *nused = ngood;
  }
  *counts = ct;

  free(values);
  
}

int main(int argc, char *argv[])
{
  float **image, pixscale, xc, yc, dist1, dist2;
  double r1, r2, pi;
  double *dist,*di;
  float **exposure, exppixscale, rmax,nsigma,maxcounts;
  float threshfac, mindist, counts, area, mincts, rate, sb = 0, errsq = 0;
  unsigned long nx,ny,npix,i,j,k,n=0, n1, n2=0, ngood, p, ntot;
  size_t *xxn, *yyn, *id;
  unsigned long expnx, expny, x, y, nt=0;
  int *rej, grow=1, *temprej, nloop, doexp = 0, doback = 0, arcmin =0;

  pi = 2.*atan2(1.,0.);
  if (argc < 6) {
    printf("radial fitsim xc yc mindist mincts rmax [maxcounts] [nsigma] [expmap] [thresh]\n");
    printf("(C) 2008 A. Mahdavi\n");
    printf("Given fitsim, an image with a valid WCS, will calculate\n");
    printf("a 1D radial profile centered on xc, yc. The minimum  \n");
    printf("radial bin size in arcsec is mindist, and the minimum \n");
    printf("counts in each bin is mincts. The profile is computed to a\n");
    printf("maximum radius rmax in arcsec, or to the edge of the image\n");
    printf("if rmax is 0. If mindist or rmax are negative, they will be\n");
    printf("taken to be |mindist| arcmin or |rmax| rmin. In either\n");
    printf("case, all output will be in units of arcminutes. If\n");
    printf("yc is negative, the output will show both the inner\n");
    printf("and the outer annulus radius; otherwise an average is shown.\n");
    printf("Only the required subset pixels will be read. Maxcounts sets\n");
    printf("the saturation threshold; all pixels with counts > maxcounts\n");
    printf("will be ignored. An nsigma (default 3) rejection of the rest\n");
    printf("of the pixels pixels is also done. Expmap is the name of an\n");
    printf("optional exposure map used for weighting the data. \n");
    printf("If mincts is negative, background subtraction of the\n");
    printf("profile will be attempted, with |mincts| counts per bin, and\n");
    printf("with thresh*|mincts| per bin in the background annulus. \n");
    printf("Normally thresh=10; you can optionally specify another value.\n");
    exit(-1);
  }
  int annulusmode=0;
  int subsetmode=1;

  xc = atof(argv[2]);
  yc = atof(argv[3]);
  if (xc < 0) { subsetmode=0; xc = -xc; }
  if (yc < 0) { annulusmode=1; yc = -yc; }
  xc -=1; yc -= 1;
  mindist = atof(argv[4]);
  if (mindist < 0) { arcmin = 1; mindist *= -60.; }
  mincts = atof(argv[5]);
  if (mincts < 0) { doback = 1; mincts = -mincts; }
  rmax = atof(argv[6]);
  if (index(argv[6],'-') != NULL) { arcmin = 1; rmax *= -60; }

  if (subsetmode) {
    image = readimage_subset(argv[1],&pixscale,&nx,&ny,xc,yc,2*rmax+1.,NULL);
  
    xc = nx/2.+0.5;
    yc = ny/2.+0.5;
  } else
    image = readimage(argv[1],&pixscale,&nx,&ny,NULL);

  npix = nx*ny;
  
  if (image == NULL) { printf ("No image!\n"); return 0; }

  if (argc > 7) 
    maxcounts = atof(argv[7]);
  else
    maxcounts = 1.E30;

  if (maxcounts < 0.) maxcounts = 1.E30;

  if (argc > 8) 
    nsigma = atof(argv[8]);
  else
    nsigma = 3;

  if (argc > 9) {
    exposure = readimage(argv[9],&exppixscale,&expnx,&expny,NULL);
    if (exposure != NULL) {
      doexp = 1;
      if (fabs(exppixscale-pixscale) > EPS || expnx != nx || expny != ny) {
	fprintf(stderr,"Error---the data and exposure images are not commensurate\n");
	return -1;
      }
    }
  } else
    exposure = NULL;

  if (argc > 10)
    threshfac = atof(argv[10]);
  else 
    threshfac = 10.;

  //printf("OK %ld %ld %ld\n",nx,ny,npix);

  dist = sci_dvector(npix);
  di = sci_dvector(npix);
  xxn  = sci_sizetvector(npix);
  yyn  = sci_sizetvector(npix);
  id  = sci_sizetvector(npix);
  rej = sci_ivector(npix);
  temprej = sci_ivector(npix);
  
  for (i = 0; i < nx; ++i)
    for (j = 0; j < ny; ++j) {
      dist1 = i-xc;
      dist2 = j-yc;
      di[n] = sqrt(dist1*dist1+dist2*dist2);
      xxn[n] = i;
      yyn[n] = j;
      id[n] = n;
      if (exposure != NULL) { temprej[n] = (exposure[i][j] < 0.); }
      else temprej[n] = 0;
      nt += temprej[n];
      ++n;
    }
  
  gsl_sort_index(id,di,1,npix);
  for (i = 0; i < npix; ++i) 
    dist[i] = di[id[i]];

  r1 = r2 = 0.; n1 = 0;
  rmax /= pixscale;

  while (r2 < rmax) {
    
    r2 = r1+mindist/pixscale;
    counts = mincts-1.;
    ngood = 0;
    //printf("Doing %e %e %e %e %e\n",r1,r2,rmax,mindist,pixscale);
    while ((counts < mincts || (mincts > 1.E-6 && ngood < 3)) && r2 <= rmax) {
    
      n2 = gsl_interp_bsearch(dist,r2,0,npix-1);
      
      for (i = n1; i <= n2; ++i) {
	p = id[i];
	rej[p] = temprej[p];
      }

      //printf("%E %E %E %E %ld %ld %ld\n",pixscale,r1,r2,counts,n1,n2,npix-1);
      reject(id,xxn,yyn,image,n1,n2,rej,&ngood,&counts,
	     exposure,nsigma,maxcounts);
      //printf("%E %E %E %ld %ld %ld %ld\n",r1,r2,counts,n1,n2,npix-1,ngood);

      if (mincts < 1.e-6) counts = 1;
      else if (counts < mincts || ngood < 3) { r2 += 1.; }
      //printf("Reject done. %ld %ld %f %f\n",n1,ngood,counts,mincts);

    }

    if (r2 > rmax) 
      rmax = r1;

    for (i = n1; i <= n2; ++i) {
      p = id[i];
      temprej[p] = rej[p];
    }

    //printf("%f %f %f\n",r1,r2,counts/(ngood*pixscale*pixscale));
    n1 = n2;
    r1 = r2;
  }
  
  for (i = 0; i <= npix-1; ++i) {
    p = (unsigned long)id[i];
    if (rej[p] && di[p] > 5)
      for (j = xxn[p]-grow; j <= xxn[p]+grow; ++j)
	for (k = yyn[p]-grow; k <= yyn[p]+grow; ++k)
	  if (j > 0 && k > 0 && j < nx && k < ny)
	    image[j][k] = -1.E30;
  }

  r1 = r2 = dist[npix];
  n2 = n1 = npix;

  unsigned long nmax;
  nmax = gsl_interp_bsearch(dist,rmax,0,npix-1);
  r2 = rmax;
  n2 = nmax;
  
  double back =0.,backerrsq = 0.,thresh;
  
  while (n1 > 1) {

    r1 = r2-mindist/pixscale;
    if (r1 <= 0) r1 = 0;
    counts = 0.;
    ngood=0;
    nloop = 0;
    if (n2 == nmax && doback) thresh = mincts*threshfac;
    else thresh = mincts;
    if (thresh < 1.e-6) thresh=1.;

    while (counts < thresh  && n1 > 1) {

      ++nloop;
      n1 = gsl_interp_bsearch(dist,r1,0,npix-1);
      if (n1 < 1) n1 = 1;
      ngood = ntot = 0;
      counts = sb = errsq = 0;

      for (i = n1; i <= n2; ++i) {
	p = id[i];
	x = xxn[p];
	y = yyn[p];
	++ntot;
	if (image[x][y] > -1.E10) {
	  ++ngood;
	  if (doexp) {
	    if (exposure[x][y] > 0) {
	      counts += image[x][y]-back*exposure[x][y];
	      rate = image[x][y]/exposure[x][y]-back;
	      sb += rate;
	      errsq += image[x][y]/(exposure[x][y]*exposure[x][y])+backerrsq;
	    } else
	      --ngood;
	  } else {
  	    counts += image[x][y]-back;
	    errsq += fabs(image[x][y])+backerrsq;
	  }
	}
      }

      for (i = n1; i <= n2; ++i) {
	p = id[i];
	x = xxn[p];
	y = yyn[p];
      }

      if (doexp)
	if (sb < sqrt(errsq)) counts = 0;

      if (mincts < 1.e-6) thresh=counts-1;
      else if (counts < thresh && r1 > 0) { r1 -= 1.; }
    }
    //printf("*************************************************");

    if (r1 < 0) r1 = 0;
    float rp1,rp2;
    char diststr[80];

    if (ngood > 3) {
      area = ngood*pixscale*pixscale/(arcmin ? 3600. : 1.);
      //printf("%d %E %E %ld %e %e %e\n",nloop,r1, r2,ngood,3.14159*(r2*r2-r1*r1)*pixscale*pixscale,area,counts);
      rp1 = r1*pixscale/(arcmin ? 60. : 1.);
      rp2 = r2*pixscale/(arcmin ? 60. : 1.);
      if (annulusmode)
	sprintf(diststr,"%f %f",rp1,rp2);
      else
	sprintf(diststr,"%f",(rp1+rp2)/2.);
      if (doback && n2 == nmax) {
	back = counts/ngood;
	if (doexp) back = sb/ngood;
	backerrsq = errsq/(ngood*ngood);
      } else {
	if (doexp)
	  printf("%s %E %E\n",
		 diststr,sb/area,sqrt(errsq)/area);
	else
	  printf("%s %E %E\n",diststr,counts/area,sqrt(errsq)/area);
      }
    }

    n2 = n1;
    r2 = r1;
  }

  for (i = 0; i < nx; ++i) 
    for (j = 0; j < ny; ++j) 
      image[i][j] = (image[i][j] > -1.E30 ? image[i][j] : 0);

  if (!subsetmode) {
    unlink("mask.fits");
    writeimage("mask.fits",image,nx,ny,argv[1]);
  }

  return 0;
}
  
