/* alloc.c
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

#include <stdio.h>
#include <stdlib.h>
#include <sciutils.h>
#include <string.h>

#if HAVE_CONFIG_H
# include <config.h>
#endif
#undef malloc
    
#include <sys/types.h>

void *malloc ();

/* Allocate an N-byte block of memory from the heap.
   If N is zero, allocate a 1-byte block.  */
   
void* rpl_malloc (size_t n)
{
 if (n == 0)
    n = 1;
 return malloc (n);
}

float *sci_fvector(long nelem)
{
  long i;
  float *x=NULL;

  if (nelem > 0)
    x = (float *)malloc(nelem*sizeof(float));
  if (x == NULL) 
    BYE("Could not allocate float vector of length %ld",nelem)

  for (i = 0; i < nelem; ++i) x[i] = 0.;
  return x;
}

void sci_fvector_resize(float **fvec, long newsize)
{
  float *tempvec;

  tempvec = (float *)realloc(*fvec,newsize*sizeof(float));
  *fvec = tempvec;
}

double *sci_dvector(long nelem)
{
  long i;
  double *x=NULL;

  if (nelem > 0)
    x = (double *)malloc(nelem*sizeof(double));
  if (x == NULL) 
    BYE("Could not allocate double vector of length %ld",nelem)

  for (i = 0; i < nelem; ++i) x[i] = 0.;
  return x;
}

void sci_dvector_resize(double **dvec, long newsize)
{
  double *tempvec;

  tempvec = (double *)realloc(*dvec,newsize*sizeof(double));
  *dvec = tempvec;
}

int *sci_ivector(long nelem)
{
  int *x=NULL;
  long i;
 
  if (nelem > 0)
    x = (int *)malloc(nelem*sizeof(int));
  if (x == NULL) 
    BYE("Could not allocate int vector of length %ld",nelem)

  for (i = 0; i < nelem; ++i) x[i] = 0;
  return x;
}

size_t *sci_sizetvector(long nelem)
{
  size_t *x=NULL;
  long i;
 
  if (nelem > 0)
    x = (size_t *)malloc(nelem*sizeof(size_t));
  if (x == NULL) 
    BYE("Could not allocate size_t vector of length %ld",nelem)

  for (i = 0; i < nelem; ++i) x[i] = 0;
  return x;
}

double **sci_dmatrix(long nx, long ny)
{
  double **x=NULL;
  long i;

  if (nx > 0 && ny > 0) 
    x = (double **)malloc(nx*sizeof(double *));
  if (x == NULL) 
    BYE("Could not allocate double matrix with %ld rows",nx)

  for (i = 0; i < nx; ++i) x[i] = sci_dvector(ny);
  return x;
}

float **sci_fmatrix(long nx, long ny)
{
  float **x=NULL;
  long i;

  if (nx > 0 && ny > 0) 
    x = (float **)malloc(nx*sizeof(float *));
  if (x == NULL) 
    BYE("Could not allocate float matrix with %ld rows",nx);

  x[0] = (float *)malloc(nx*ny*sizeof(float));
  if (x[0] == NULL) 
    BYE("Could not allocate %ldx%ld float matrix",nx,ny);

  for (i = 1; i < nx; ++i) 
    x[i] = x[0]+i*ny;

  return x;
}

int **sci_imatrix(long nx, long ny)
{
  int **x=NULL;
  long i;

  if (nx > 0 && ny > 0) 
    x = (int **)malloc(nx*sizeof(int *));
  if (x == NULL) 
    BYE("Could not allocate int matrix with %ld rows",nx)

  for (i = 0; i < nx; ++i) x[i] = sci_ivector(ny);
  return x;
}

char **sci_stringmatrix(long nx, long slen)
{
  char **x=NULL;
  long i;

  if (nx > 0 && slen > 0)
    x = (char **)malloc(sizeof(char *)*nx);
  if (x == NULL)
    BYE("Could not allocate string matrix with %ld rows",nx);

  for (i = 0; i < nx; ++i) {
    x[i] = (char *)malloc(sizeof(char)*slen);
    if (x[i] == NULL)
      BYE("Could not allocate string matrix with %ld rows",nx);
  }

  return x;
}

void sci_free_stringmatrix(char **x, long nx)
{
  long i;

  for (i = 0; i < nx; ++i) 
    free(x[i]);

  free(x);

}

double ***sci_dtensor(long nx, long ny, long nz)
{
  double ***x=NULL;
  long i;

  if (nx > 0 && ny > 0 && nz > 0) 
    x = (double ***)malloc(nx*sizeof(double **));
  if (x == NULL) 
    BYE("Could not allocate double tensor  with %ld rows",nx)

  for (i = 0; i < nx; ++i) x[i] = sci_dmatrix(ny,nz);
  return x;
}

float ***sci_ftensor(long nx, long ny, long nz)
{
  float ***x=NULL;
  long i;

  if (nx > 0 && ny > 0 && nz > 0) 
    x = (float ***)malloc(nx*sizeof(float **));
  if (x == NULL) 
    BYE("Could not allocate float tensor with %ld rows",nx)

  for (i = 0; i < nx; ++i) x[i] = sci_fmatrix(ny,nz);
  return x;
}

void sci_free_fmatrix(float **x)
{
  long i;
  
  free(x[0]);
  free(x);
}

void sci_free_dmatrix(double **x, long nx)
{
  long i;
  
  for (i = 0; i < nx; ++i) free(x[i]);
  free(x);
}

void sci_free_imatrix(int **x, long nx)
{
  long i;
  
  for (i = 0; i < nx; ++i) free(x[i]);
  free(x);
}

void sci_free_dtensor(double ***x, long nx, long ny)
{
  long i;
  
  for (i = 0; i < nx; ++i) sci_free_dmatrix(x[i],ny);
  free(x);
}

void sci_free_ftensor(float ***x, long nx, long ny)
{
  long i;
  
  for (i = 0; i < nx; ++i) sci_free_fmatrix(x[i]);
  free(x);
}

char *sci_strdup(const char *s)
{
  size_t sl;
  char *outs;

  sl = strlen(s)+1;
  if (sl > 10000) return NULL;
  outs = (char *)malloc(sl*sizeof(char));
  strncpy(outs,s,sl);
  return outs;
  
}

char *sci_strcat(const char *s1, const char *s2)
{
  char *scat;

  scat = (char *)malloc(sizeof(char)*(strlen(s1)+strlen(s2)+1));
  strcpy(scat,s1);
  strcat(scat,s2);

  return scat;
}
