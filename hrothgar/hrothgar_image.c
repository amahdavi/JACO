/* hrothgar_image.c
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
#include <hrothgar.h>
#include <hrothgar_proto.h>
#include <sciutils.h>

void hrothgar_writeimage(char *outfile,
			 float *buffer, 
			 double x1, double x2, double y1, double y2,
			 long ngrid)
{
  fitsfile *fout;
  int status=0;
  long fpixel = 1,naxes[2];
  unsigned long npix;
  float val;

  
  if(fits_create_file(&fout,outfile,&status)) 
    BYE("Error creating fits file (code %d). Try --overwrite.",status);

  naxes[0] = ngrid;
  naxes[1] = ngrid;
  
  double xstep = (x2-x1)/ngrid;
  double ystep = (y2-y1)/ngrid;
  
  fits_create_img(fout,-32,2,naxes,&status);
  
  npix = ngrid*ngrid;

  fits_write_img(fout,TFLOAT,fpixel,npix,buffer,&status);

  fits_write_key(fout,TSTRING,"WCSNAMEP","PHYSICAL","Parameter Coord. System",
		      &status);
  fits_write_key(fout,TSTRING,"WCSTY1P","PHYSICAL","Parameter Coord. System",
		      &status);
  fits_write_key(fout,TSTRING,"CTYPE1P","x","Parameter 1",&status);
  val = x1+xstep/2.;
  fits_write_key(fout,TFLOAT,"CRVAL1P",&val,"Parameter 1",&status);
  val = 0.5;
  fits_write_key(fout,TFLOAT,"CRPIX1P",&val,"Parameter 1",&status);
  val = xstep;
  fits_write_key(fout,TFLOAT,"CDELT1P",&val,"Parameter 1",&status);
  fits_write_key(fout,TSTRING,"WCSTY2P","PHYSICAL","Parameter Coord. System",
		      &status);
  fits_write_key(fout,TSTRING,"CTYPE2P","y","Parameter 2",&status);
  val = y1+ystep/2.;   
  fits_write_key(fout,TFLOAT,"CRVAL2P",&val,"Parameter 2",&status);
  val = 0.5;   
  fits_write_key(fout,TFLOAT,"CRPIX2P",&val,"Parameter 2",&status);
  val = ystep;   
  fits_write_key(fout,TFLOAT,"CDELT2P",&val,"Parameter 2",&status);

  fits_close_file(fout,&status);
  
  
}

void hrothgar_writeimage_double(char *outfile,
				double **dbuffer,
				double x1, double x2, double y1, double y2,
				long ngrid)
{
  float *buffer = sci_fvector(ngrid*ngrid);
  
  long i,j;

  for (i = 0; i < ngrid; ++i)
    for (j = 0; j < ngrid; ++j)
      buffer[j+ngrid*i] = dbuffer[i][j];

  hrothgar_writeimage(outfile,buffer,x1,x2,y1,y2,ngrid);
  

  free(buffer);

}

void hrothgar_wcsimage(char *outfile,
		       float *buffer, 
		       float wcs[3][2],
		       long ngrid)
{
  fitsfile *fout;
  int status=0;
  long fpixel = 1,naxes[2];
  unsigned long npix;
  float val;

  
  if(fits_create_file(&fout,outfile,&status)) 
    BYE("Error creating fits file (code %d). Try --overwrite.",status);

  naxes[0] = ngrid;
  naxes[1] = ngrid;
  
  fits_create_img(fout,-32,2,naxes,&status);
  
  npix = ngrid*ngrid;

  fits_write_img(fout,TFLOAT,fpixel,npix,buffer,&status);

  fits_write_key(fout,TSTRING,"CTYPE1","RA--TAN","",
		      &status);
  fits_write_key(fout,TSTRING,"CUNIT1","deg","",
		      &status);
  fits_write_key(fout,TFLOAT,"CRVAL1",&wcs[1][0],"",&status);
  fits_write_key(fout,TFLOAT,"CRPIX1",&wcs[0][0],"",&status);
  fits_write_key(fout,TFLOAT,"CDELT1",&wcs[2][0],"",&status);
  fits_write_key(fout,TSTRING,"CTYPE2","DEC--TAN","",
		      &status);
  fits_write_key(fout,TSTRING,"CUNIT2","deg","",
		      &status);
  fits_write_key(fout,TFLOAT,"CRVAL2",&wcs[1][1],"",&status);
  fits_write_key(fout,TFLOAT,"CRPIX2",&wcs[0][1],"",&status);
  fits_write_key(fout,TFLOAT,"CDELT2",&wcs[2][1],"",&status);

  fits_close_file(fout,&status);
  
  
}
