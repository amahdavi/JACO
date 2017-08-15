/* readxray.c
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

#include <fitsio.h>
#include <string.h>
#include <math.h>
#include <stdio.h>
#include <sciutils.h>
#include <unistd.h>

int readarf(char *arfname, struct arf *result)
{
  int status = 0;
  fitsfile *arffile;
  long nrows;
  int fnull;

  fits_open_file(&arffile,arfname,READONLY,&status);
  if (status) BYE("cfitsio error %d",status);
  fits_movnam_hdu(arffile,BINARY_TBL,"SPECRESP",0,&status);
  if (status) BYE("cfitsio error %d",status);

  fits_get_num_rows(arffile,&nrows,&status);
  result->nebins = (int)nrows;

  if (status) BYE("cfitsio error %d",status);
  result->elo = sci_dvector(result->nebins);
  result->ehi = sci_dvector(result->nebins);
  result->effarea = sci_dvector(result->nebins);
  
  fits_read_col_dbl(arffile,1,1,1,result->nebins,0,result->elo,&fnull,&status);
  fits_read_col_dbl(arffile,2,1,1,result->nebins,0,result->ehi,&fnull,&status);
  fits_read_col_dbl(arffile,3,1,1,result->nebins,0,result->effarea,
		    &fnull,&status);

  fits_close_file(arffile,&status);
  if (status) BYE("cfitsio error %d",status);
  return 0;
}
  
int readrmf(char *rmfname,struct rmf *result,struct ebounds *elimits)
{
  int status = 0;
  int i,j;
  fitsfile *rmffile;
  int fnull,elem;
  long nrows,ntotchannels;
  double *longvec;

  fits_open_file(&rmffile,rmfname,READONLY,&status);
  if (status) BYE("cfitsio error %d",status);
  fits_movnam_hdu(rmffile,BINARY_TBL,"MATRIX",0,&status);
  if (status) {
    status = 0;
    fits_movnam_hdu(rmffile,BINARY_TBL,"SPECRESP MATRIX",0,&status);
    if (status) BYE("Couldn't find MATRIX or SPECRESP MATR");
  }

  fits_get_num_rows(rmffile,&nrows,&status);
  if (status) BYE("cfitsio error %d",status);
  result->nebins = (int)nrows;

  result->elo = sci_dvector(result->nebins);
  result->ehi = sci_dvector(result->nebins);
  result->ngrp = sci_ivector(result->nebins);
  result->chanstart = (int **)malloc(result->nebins*sizeof(int *));
  result->nchannels = (int **)malloc(result->nebins*sizeof(int *));
  result->response = (double ***)malloc(result->nebins*sizeof(double **));

  fits_read_col_dbl(rmffile,1,1,1,result->nebins,0,result->elo,&fnull,&status);
  if (status) BYE("cfitsio error %d\n",status);
  fits_read_col_dbl(rmffile,2,1,1,result->nebins,0,result->ehi,&fnull,&status);
  if (status) BYE("cfitsio error %d\n",status);
  fits_read_col_int(rmffile,3,1,1,result->nebins,0,result->ngrp,&fnull,&status);
  if (status) BYE("cfitsio error %d\n",status);

  for (i = 0; i < result->nebins; ++i) {

    if (result->ngrp[i] > 0) {
      result->nchannels[i] = sci_ivector(result->ngrp[i]);
      result->chanstart[i] = sci_ivector(result->ngrp[i]);
      
      result->response[i] = (double **)malloc(result->ngrp[i]*sizeof(double *));
      
      fits_read_col_int(rmffile,4,i+1,1,result->ngrp[i],0,
			result->chanstart[i],&fnull,&status);
      if (status) BYE("cfitsio error %d\n",status);
      fits_read_col_int(rmffile,5,i+1,1,result->ngrp[i],0,
		      result->nchannels[i],&fnull,&status);
      if (status) BYE("cfitsio error %d\n",status);
      
      ntotchannels = 0;
      for (j = 0; j < result->ngrp[i]; ++j) 
	ntotchannels += result->nchannels[i][j];
      longvec = sci_dvector(ntotchannels);
      
      fits_read_col_dbl(rmffile,6,i+1,1,ntotchannels,0,longvec,&fnull,&status);
      
      elem = 0;
      for (j = 0; j < result->ngrp[i]; ++j) {
	result->response[i][j] = &longvec[elem];
	elem += result->nchannels[i][j];
      }
    }
  }
  fits_movnam_hdu(rmffile,BINARY_TBL,"EBOUNDS",0,&status);
  if (status) BYE("cfitsio error %d",status);

  fits_get_num_rows(rmffile,&nrows,&status);
  if (status) BYE("cfitsio error %d",status);
  elimits->nchannels = (int)nrows;

  elimits->emin = sci_dvector(nrows);
  elimits->emax = sci_dvector(nrows);
  fits_read_col_dbl(rmffile,2,1,1,nrows,0,
		    elimits->emin,&fnull,&status);
  if (status) BYE("cfitsio error %d\n",status);
  fits_read_col_dbl(rmffile,3,1,1,nrows,0,
		    elimits->emax,&fnull,&status);
  if (status) BYE("cfitsio error %d\n",status);
  
  fits_close_file(rmffile,&status);

  return 0;

}

int readpha(char *phaname, struct pha *result, int read_back)
{
  int status = 0,colnum;
  char comment[100];
  fitsfile *phafile;
  int fnull;
  long nrows;

  fits_open_file(&phafile,phaname,READONLY,&status);
  if (status) BYE("cfitsio errors %d",status);
  fits_movnam_hdu(phafile,BINARY_TBL,"SPECTRUM",0,&status);
  if (status) BYE("cfitsio errors %d",status);

  fits_get_num_rows(phafile,&nrows,&status);
  if (status) BYE("cfitsio errors %d",status);

  result->filename = phaname;
  result->nchannels = (int)nrows;

  result->channel = sci_ivector(result->nchannels);
  result->quality = sci_ivector(result->nchannels);
  result->grouping = sci_ivector(result->nchannels);
  result->counts = sci_dvector(result->nchannels);
  result->model = sci_dvector(result->nchannels);
  result->stat_err = sci_dvector(result->nchannels);

  fits_get_colnum(phafile,CASEINSEN,"CHANNEL",&colnum,&status);
  fits_read_col_int(phafile,colnum,1,1,result->nchannels,0,
		    result->channel,&fnull,&status);
  if (status) BYE("Cannot read CHANNEL column");

  fits_get_colnum(phafile,CASEINSEN,"QUALITY",&colnum,&status);
  fits_read_col_int(phafile,colnum,1,1,result->nchannels,0,
		    result->quality,&fnull,&status);

  // Ignore if there is no quality information
  if (status) status = 0;

  fits_read_key_dbl(phafile,"EXPOSURE",&result->exposure,comment,&status);
  if (status) BYE("Cannot read EXPOSURE keyword");
  fits_read_key_dbl(phafile,"BACKSCAL",&result->arcminarea,comment,&status);
  if (status) BYE("Cannot read BACKSCAL keyword");

  fits_get_colnum(phafile,CASEINSEN,"COUNTS",&colnum,&status);
  fits_read_col_dbl(phafile,colnum,1,1,result->nchannels,0,
		    result->counts,&fnull,&status);
  if (status) BYE("Cannot read COUNTS column");

  if (read_back) return;

  result->arfname= (char *)malloc(100*sizeof(char));
  fits_read_key_str(phafile,"ANCRFILE",result->arfname,comment,&status);
  if (status) BYE("Cannot read ARF file name");
  result->rmfname= (char *)malloc(100*sizeof(char));
  fits_read_key_str(phafile,"RESPFILE",result->rmfname,comment,&status);
  if (status) BYE("Cannot read RMF file name");

  fits_get_colnum(phafile,CASEINSEN,"GROUPING",&colnum,&status);
  fits_read_col_int(phafile,colnum,1,1,result->nchannels,0,
		    result->grouping,&fnull,&status);
  if (status) BYE("Cannot read GROUPING column");

  fits_get_colnum(phafile,CASEINSEN,"STAT_ERR",&colnum,&status);
  fits_read_col_dbl(phafile,colnum,1,1,result->nchannels,0,
		    result->stat_err,&fnull,&status);
  if (status) {
    status = 0;
    result->backname= (char *)malloc(100*sizeof(char));
    fits_read_key_str(phafile,"BACKFILE",result->backname,comment,&status);
    if (status) BYE("No STAT_ERR column and no BACKFILE for %s.",phaname);

    printf("%s: background=%s\n",phaname,result->backname);
    result->background = (struct pha *)malloc(sizeof(struct pha));
    readpha(result->backname,result->background,1);
  } else
    result->backname = NULL;

  
  fits_close_file(phafile,&status);

  return 0.;
}

int writepha(char *outname, struct pha *result)
{
  int status = 0,chancol,errcol;
  fitsfile *fin,*phafile;
  int i,j,chansize;
  double meanflux,meanerr;
  long nrows;

  unlink(outname);
  fits_open_file(&fin,result->filename,READONLY,&status);
  if (status) BYE("cfitsio errors %d",status);
  fits_create_file(&phafile,outname,&status);
  if (status) BYE("cfitsio errors %d",status);
  fits_copy_file(fin,phafile,1,1,1,&status);
  if (status) BYE("cfitsio errors %d",status);
  fits_close_file(fin,&status);

  fits_movnam_hdu(phafile,BINARY_TBL,"SPECTRUM",0,&status);
  if (status) BYE("cfitsio errors %d",status);
  fits_get_num_rows(phafile,&nrows,&status);
  if (status) BYE("cfitsio errors %d",status);

  if (result->nchannels != nrows) 
    BYE("Writing to non-identically formatted template in writepha.");

/*   result->channel = sci_ivector(result->nchannels); */
/*   result->quality = sci_ivector(result->nchannels); */
/*   result->grouping = sci_ivector(result->nchannels); */
/*   result->counts = sci_dvector(result->nchannels); */
/*   result->model = sci_dvector(result->nchannels); */
/*   result->stat_err = sci_dvector(result->nchannels); */


  fits_get_colnum(phafile,CASEINSEN,"COUNTS",&chancol,&status);
  fits_get_colnum(phafile,CASEINSEN,"STAT_ERR",&errcol,&status);

  for (i = 0; i < result->nbins; ++i) {
    chansize = result->chanhi[i]-result->chanlo[i]+1;
    if (chansize <= 0) BYE("Invalid channel size in writepha");

    meanflux = result->exposure*result->mcounts[i]/chansize;
    meanerr = result->exposure*sqrt(result->sigsq[i]/chansize);
    
    for (j = result->chanlo[i]; j <= result->chanhi[i]; ++j)  {
      fits_write_col_dbl(phafile,chancol,j,1,1,&meanflux,&status);
      fits_write_col_dbl(phafile,errcol,j,1,1,&meanerr,&status);
    }

  }

  fits_close_file(phafile,&status);

  return 0.;
}

void bin_pha(struct pha *spec, struct ebounds *elimits, 
	     double efit1, double efit2, double syserr)  
{
  int i,j=-1,chanstart=-1,chanstop=-1,lastchan=-1;
  double expsq;

  spec->nbins=-1;
  for (i = 0; i < spec->nchannels; ++i)
    if (elimits->emin[i]-efit1 > 1.E-6 && 
	elimits->emax[i]-efit2 < 1.E-6) {
      lastchan = i;
      if (spec->quality[i] == 0 && spec->grouping[i] != -1) {
	++spec->nbins;
	if (spec->nbins == 0) chanstart = i;
	chanstop = i-1;
      }
    }

  if (lastchan < spec->nchannels-1)
    if (spec->quality[lastchan+1] == 0 && spec->grouping[lastchan+1] != -1) {
      ++spec->nbins;
      chanstop = lastchan;
    }

  if (spec->nbins <= 0)
    BYE("Spectrum %s improperly grouped or has too few counts  (%d).\n",
	spec->filename,spec->nbins);

  struct pha *background = (struct pha*)spec->background;

  spec->sigsq = sci_dvector(spec->nbins);
  spec->sig = sci_dvector(spec->nbins);
  spec->chisq = sci_dvector(spec->nbins);
  spec->flux = sci_dvector(spec->nbins);
  spec->mcounts = sci_dvector(spec->nbins);
  spec->chanlo = sci_ivector(spec->nbins);
  spec->chanhi = sci_ivector(spec->nbins);

  double back_k=0.;
  if (spec->backname)
    back_k = spec->exposure*spec->arcminarea/
      (background->exposure*background->arcminarea);

  for (i = chanstart; i <= chanstop; ++i) 
    if (spec->quality[i] == 0) {
      if (spec->grouping[i] != -1) {	
	++j;
	spec->sigsq[j] = 0.;
	spec->flux[j] = 0.;
	spec->chanlo[j] = i;
      } 
      spec->chanhi[j] = i;
      if (j < 0 || j >= spec->nbins) 
	BYE("Error in PHA file format at channel %d\n",i);
      spec->flux[j] += spec->counts[i];
      //printf("zzzzzzzzz %f %d %d %f %f\n",spec->counts[i],j,spec->nbins,elimits->emin[i],spec->flux[j]);
      if (spec->backname) {
	spec->flux[j] -= back_k*background->counts[i];
	spec->stat_err[i] = sqrt(spec->counts[i]+back_k*back_k*
				 background->counts[i]);
      }
      spec->sigsq[j] += spec->stat_err[i]*spec->stat_err[i];
    }
  

  if (!(spec->chanhi[j] > 0.))
    BYE("Error in PHA file format: high channel is %d\n",spec->chanhi[j]);

  expsq = spec->exposure*spec->exposure;
  for (i = 0; i < spec->nbins; ++i) {
    spec->flux[i] /= spec->exposure;
    spec->sigsq[i] /= expsq;
    spec->sigsq[i] += pow(spec->flux[i]*syserr,2.0);
    spec->sig[i] = sqrt(spec->sigsq[i]);
  }

  return;
}

float **readimage_subset(char *fname, float *pixscale,
			 unsigned long *nx, unsigned long *ny,
			 long xcen, long ycen, float dist, float wcs[3][2])
{
  fitsfile *fptr;
  int status=0,nfound,anynul;
  long i, j, fpixel[2], lpixel[2], naxes[2], inc[2];
  unsigned long npixels;
  float nullval=0.;
  float **m,*buff,wcsscale[2];
  char *comment=NULL;

  if (fits_open_file(&fptr,fname,READONLY,&status)) {
    fprintf(stderr,"Error opening image %s. Error %d\n",fname,status);
    return NULL;
  }

  if (fits_read_keys_lng(fptr,"NAXIS",1,2,naxes,&nfound,&status)) {
    fprintf(stderr,"Could not find NAXIS keyword in %s.\n",fname);
    return NULL;
  }
  
  status = 1;
  for (i = 1; i <= 6 && status != 0; ++i) {
    status = 0;
    switch(i) {

    case 1:
      nfound = 0;
      if (fits_read_key_flt(fptr,"CD1_1",&wcsscale[0],comment,&status)) { 
	wcsscale[0] = 0.; 
	status = 0; 
      } else ++nfound;

      if (fits_read_key_flt(fptr,"CD1_2",&wcsscale[1],comment,&status)) { 
	wcsscale[1] = 0.; 
	status = 0; 
      } else ++nfound;

      *pixscale = 3600.*sqrt(wcsscale[0]*wcsscale[0]+
			     wcsscale[1]*wcsscale[1]);
      //printf("Boo %E %E\n",wcsscale[0],wcsscale[1]);
      status = (nfound == 0);
      break;
    case 2:
      if (!fits_read_key_flt(fptr,"CDELT1",pixscale,comment,&status)) {
	if (wcs != NULL) {
	  if (fits_read_key_flt(fptr,"CRPIX1",&wcs[0][0],comment,&status) ||
	      fits_read_key_flt(fptr,"CRPIX2",&wcs[0][1],comment,&status) ||
	      fits_read_key_flt(fptr,"CRVAL1",&wcs[1][0],comment,&status) ||
	      fits_read_key_flt(fptr,"CRVAL2",&wcs[1][1],comment,&status) ||
	      fits_read_key_flt(fptr,"CDELT2",&wcs[2][1],comment,&status)) {
	    fprintf(stderr,"Incomplete WCS header.\n");
	    return NULL;
	  } else 
	    wcs[2][0] = *pixscale;
	}
	*pixscale = fabs(*pixscale)*3600.;
      }
      break;
    case 3:
      if (!fits_read_key_flt(fptr,"PLTSCALE",&wcsscale[0],comment,&status)) {
	fits_read_key_flt(fptr,"XPIXELSZ",&wcsscale[1],comment,&status);
	*pixscale = wcsscale[0]*wcsscale[1]/1000.;
      }
      break;
    case 4:
      fits_read_key_flt(fptr,"PIXSCALE",pixscale,comment,&status);
      break;
    case 5:
      fits_read_key_flt(fptr,"SECPIX",pixscale,comment,&status);
      break;
    case 6:
      fprintf(stderr,"No pixel scale information found; assuming 1\"/pixel.\n");
      *pixscale = 1;
      break;

    }
  }
  //printf("%f %ld\n",*pixscale,npixels);
    
  if (xcen < 0 || ycen < 0) {
    *nx = naxes[0]; *ny = naxes[1];
    npixels = naxes[0]*naxes[1];
    fpixel[0] = 1; fpixel[1] = 1;
    lpixel[0] = *nx; lpixel[1] = *ny;
  } else {
    *ny = *nx = (long)ceil(dist/(*pixscale));
    if (!((*nx) % 2)) { ++(*nx); ++(*ny); }
    npixels = (*nx)*(*ny);
    //printf("%ld %ld %E %E\n",*nx, *ny, *pixscale, dist);
    fpixel[0] = xcen-(*nx-1)/2;
    fpixel[1] = ycen-(*ny-1)/2;
    lpixel[0] = xcen+(*nx-1)/2;
    lpixel[1] = ycen+(*ny-1)/2;
    
    if (fpixel[0] < 1 || fpixel[1] < 1 ||
	lpixel[0] > naxes[0] || lpixel[1] > naxes[1])
      BYE("image subset exceeds image border.");
  }

  status=0;

  m = sci_fmatrix((*nx),(*ny));
  buff = sci_fvector(npixels);
  //printf("%ld %ld %ld\n",*nx,*ny,npixels);
  inc[0] = inc[1] = 1;

  if (fits_read_subset(fptr,TFLOAT,fpixel,lpixel,inc,&nullval,buff,
  	    &anynul,&status)) {
    printf("ERROR!\n");
    return NULL;
  }

  for (j = 0; j < (long)(*ny); ++j) 
    for (i = 0; i < (long)(*nx); ++i){
      m[i][j] = buff[j*(*nx)+i];
      //printf("%ld %ld %E\n",fpixel[0]+i,fpixel[1]+j,m[i][j]);
    }

  free(buff);

  fits_close_file(fptr,&status);
  
  return m;
}

float **readimage(char *fname, float *pixscale,
		  unsigned long *nx, unsigned long *ny, float wcs[3][2])
{
  return readimage_subset(fname, pixscale, nx, ny, -1, -1, 0., wcs);
}


void writeimage(char *outfile, float **image, unsigned long nx,
		unsigned long ny, char *infile)
{
  fitsfile *fout,*fin;
  unsigned long i,j;
  int status=0;
  long fpixel = 1;
  unsigned long npix;
  float *buff;

  if(fits_create_file(&fout,outfile,&status)) {
    fprintf(stderr,"Error creating file %s. Code %d\n",outfile,status);
    return;
  }

  if (fits_open_file(&fin,infile,READONLY,&status)) {
    fprintf(stderr,"Error opening image %s. Error %d\n",infile,status);
    return;
  }

  if (fits_copy_header(fin, fout, &status)) {
    fprintf(stderr,"Error copying headers from %s to %s. Error %d\n",
     outfile,infile,status);
    return;
  }
  
  fits_close_file(fin,&status);
  npix = nx*ny;
  buff = sci_fvector(npix);

  for (j = 0; j < ny; ++j) 
    for (i = 0; i < nx; ++i)
      buff[j*nx+i] = image[i][j];

  fits_write_img(fout,TFLOAT,fpixel,npix,buff,&status);

  free(buff);

  fits_close_file(fout,&status);
  
  
}

