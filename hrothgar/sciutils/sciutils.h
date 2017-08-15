#include <stddef.h>
#ifdef HAVE_CONFIG_H
#include <config.h>
#endif

#ifndef SCIUTILS_H
#include <execinfo.h>

struct arf {

  int nebins;
  double *elo;
  double *ehi;
  double *effarea;

};

struct rmf {

  int   nebins;
  double *elo;
  double *ehi;
  int   *ngrp;
  int   **nchannels;
  int   **chanstart;
  double ***response;

};

struct ebounds {
  
  int nchannels;
  double *emin;
  double *emax;

};

struct pha {

  int nchannels;
  char *filename;
  double exposure;
  double arcminarea;
  char  *arfname;
  char  *rmfname;
  char  *backname;
  void *background;
  int  *channel;
  double *stat_err;
  double *counts;
  double *model;
  int  *quality;
  int  *grouping;
  int nbins;
  int *chanlo;
  int *chanhi;
  double *flux;
  double *mcounts;
  double *sigsq;
  double *sig;
  double *chisq;
  double totchisq;

};


#ifndef BYE
#define BYE(format...) { void **buffer = (void **)malloc(100*sizeof(void *)); backtrace(buffer,100); backtrace_symbols_fd(buffer,100,1); printf("SciUtils Error:\n"); printf(format); printf("\n"); exit(-1); }
#endif

double genbeta(double a, double b);
double genbetai(double a, double b, double x);
double genbetaq(double a, double b, double x, 
		double betaab, double omxbxa, double psiofa);
double ztoMpc(double h0, double oM, double oL, double z);

// readdata.c
double *double_readdata(char *rd_filename, 
			int rd_col, unsigned long *rd_count, int offset);
float *float_readdata(char *rd_filename, int rd_col, 
		      unsigned long *rd_count, int offset);
char *sci_strdup(const char *s);
void grow_string_array(char ***ins, long nchar, char *string);
long read_ascii(char *filenames, const char *fmt, char ***comments, ...);
int parse_list(const char *inputlist, char delimit, char ***parsedlist);
int parse_list_double(const char *inputlist, char delimit, double **parsedlist);

			     

/* readxray.c */
int readarf(char *arfname, struct arf *result);
int readrmf(char *rmfname, struct rmf *result, struct ebounds *elimits);
int readpha(char *phaname, struct pha *result, int read_back);
int writepha(char *outname, struct pha *result);
void bin_pha(struct pha *spec, struct ebounds *elimits, double efit1, 
	     double efit2, double syserrors);
float **readimage(char *fname, float *pixscale,
		  unsigned long *nx, unsigned long *ny, float wcs[3][2]);
float **readimage_subset(char *fname, float *pixscale,
			 unsigned long *nx, unsigned long *ny,
			 long xcen, long ycen, float dist, float wcs[3][2]);
void writeimage(char *outfile, float **image, unsigned long nx,
		unsigned long ny, char *infile);

/* alloc.c */
double *sci_dvector(long nelem);
void sci_dvector_resize(double **dvec, long newsize);
float *sci_fvector(long nelem);
void sci_fvector_resize(float **dvec, long newsize);
int *sci_ivector(long nelem);
size_t *sci_sizetvector(long nelem);
float **sci_fmatrix(long nx, long ny);
double **sci_dmatrix(long nx, long ny);
int **sci_imatrix(long nx, long ny);
void sci_free_fmatrix(float **x);
void sci_free_dmatrix(double **x, long nx);
void sci_free_imatrix(int **x, long nx);
char **sci_stringmatrix(long nx, long slen);
void sci_free_stringmatrix(char **x, long nx);
double ***sci_dtensor(long nx, long ny, long nz);
void sci_free_dtensor(double ***x, long nx, long ny);
float ***sci_ftensor(long nx, long ny, long nz);
void sci_free_ftensor(float ***x, long nx, long ny);
char *sci_strcat(const char *s1, const char *s2);

#define SCIUTILS_H
#endif
