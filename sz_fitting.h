#ifndef sz_fitting_h
#define sz_fitting_h

#define SZC "SZ Code Copyright (C) 2007 Jonathan Sievers, Stephen Myers.\n"

#define SZSTRINGLEN 512


#include <complex.h>
#define MAX_SZ_DATASETS 10 /*maximum number of allowed separate SZ datasets*/





/*---------------------------------------------------------------------------*/
typedef struct dump_data_s{
  double *svec;
  double *zvec;
  int *uvec;
  int *vvec;
  int maxuv;
  int gridsize;
  int npix_use;
  int n;
  int n_nonzero;
  double pixel_size;

  int nbeam;
  double du;
  complex double  ***rgrid;



} DumpData;


/*---------------------------------------------------------------------------*/
typedef struct svec_data_s{
  double *svec;
  double *zvec;
  float *uvec;
  float *vvec;
  double ra;
  double dec;
  int n;

} SvecData;

/*--------------------------------------------------------------------------------*/



struct sz_control_params_s {

  char ParamFile[SZSTRINGLEN];
  int do_defaults;
  
  double racent;
  double deccent;
  double GridderRACent;
  double GridderDecCent;
  double pixsize;
  int npix;
  int node;                // Number of MPI nodes
  int quiet;               // Suppress all non-error messages

  double datamin;
  double datamax;
  
  double freq;

  double calib_fac;  //Apply a calibration scale to the data if not 1.  Is 0.97 for CBI
  double cmb_scale_fac;  //Apply a scale factor to the CMB noise.  CBI is Main-(Lead+Trail)/2
                         //to get rid of ground, but has the effect of multiplying CMB noise by 1.5


  int nest;
  double **eigvecs;
  double *data;
  double *work;
  double *model;
  double *data_org;

  double *eigvals;  


  char configfile[SZSTRINGLEN];

  char GridderExec[SZSTRINGLEN];
  char MockCBIExec[SZSTRINGLEN];
  char GridderScriptName[SZSTRINGLEN];
  char MockCBIScriptName[SZSTRINGLEN];
  char FitsFile[SZSTRINGLEN];

  char UVFModelName[SZSTRINGLEN];
  char UVFDataName[SZSTRINGLEN];
  char GridderDataBase[SZSTRINGLEN];
  char GridderNoiseName[SZSTRINGLEN];
  char EstimatorName[SZSTRINGLEN];

  int skip_cmb;


  char DerotDataFile[SZSTRINGLEN];
  char DifferenceDataFile[SZSTRINGLEN];
  char JacoTag[SZSTRINGLEN];




  int WriteDifference;
  int FlipEndian;
  float *mapcol;

  void *szof,*szef;       // Pointers to stdout and stderr

  double **ymap;
  double y_to_t;


  char dump_name_rvec[SZSTRINGLEN];
  char dump_name_svec[SZSTRINGLEN];
  char svec_name[SZSTRINGLEN];



  DumpData *dump;
  SvecData *svec_data;

  int oversamp;

} sz_control_params;
typedef struct sz_control_params_s SZControlParams;



/*--------------------------------------------------------------------------------*/

struct big_sz_control_params_s {
  int ndataset;
  SZControlParams *szdatasets;
  /*double *datavec;*/
  /*double *errs;  */

  int nest;
  char DerotDataFile[SZSTRINGLEN];
  double *deccent;
  double *racent;
  double *pixsize;
  int *npix;
  double ***ymap;
  int node;                // Number of MPI nodes
  int quiet;               // Suppress all non-error messages

  double *data;
  double *errors;
  double *models;

  /*char **paramfiles;*/
  /*char paramfiles[MAX_SZ_DATASETS][SZ_STRING_LEN];*/
  char **paramfiles;
} big_sz_control_params;
typedef struct big_sz_control_params_s BigSZControlParams;


/*--------------------------------------------------------------------------------*/

void print_sz_params(SZControlParams *params);
int initialize_sz(SZControlParams *szparams);
int get_sz_model(double *model, SZControlParams *params);
//int get_sz_model_dump(double *model, SZControlParams *params);
//int get_sz_model_dump(double *model, SZControlParams *params);
int write_simple_fits_header_params(FILE *outfile, SZControlParams *szparams);
void flip_float_endian(float *x);
int GetNEstimator(char *filename);
int ReadSvec(double *data,char *filename);
int ReadNoise(double **noise, char *filename, int nest_in);
int ReadCMBMat(double **noise, char *filename, int nest_in);
void dflip_line(double *line, int nelem);
int WriteDiff(SZControlParams *params, double *data);
int append_tag(char *str1, char *str2, int maxlen);
int read_szmodel_file(char *filename, double **x_out, double **y_out, double **err_out);
int big_initialize_sz(BigSZControlParams *params);
int big_get_sz_model(double *model, BigSZControlParams *params);


/*void dgemv_(char *trans,int *m,&n,&alpha,a,&lda,x,&incx,&beta,y,&incy,1);*/

void  dsyev_( char *jobz, char *uplo, int *n, double *a, int *lda, double *w, double *work, int *lwork, int *info,int jobzlen, int uplolen );
void  cdsyev( char jobz, char uplo, int n, double *a, int lda, double *w, double *work, int lwork, int *info);

void lengthen_string(char **cur_string,int *len_in, int n_to_add);
char *read_all_stdin(char *filename, SZControlParams *params);
void dd2dms(double dec, int *dd, int *dm, double *ds, char *sign);

complex double **complex_dmatrix(int n, int m);
void fft_2d_map_wshifts(double **mat, complex double **matft, int n, int m);

DumpData  *read_tt_dump_svec(char *fname);
int read_tt_dump_rvec(DumpData *data, char *fname);
void clear_dump_zeros(DumpData *dump);
SvecData *ReadSvecStruct(char *filename);
int get_npix_from_oversamp(int npix_in, int oversamp);
void transpose_matrix(double **mat, int n);
void dflip_mat(double **mat, int nelem);
double yfac(SZControlParams *data);

void dump_wisdom();
void setup_wisdom();

#endif

