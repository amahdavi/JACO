struct weaklensing {

  int nwl;
  double *wlrad;
  double *modelshear;
  double *shear;
  double *shearerr;

  /* Siegel Begin */
  double *rot_modelshear;
  double *rot_shear;
  double *rot_shearerr;

  double *eigenval;
  double **eigenvec;
  double **inv_eigenvec;
  /* Siegel End */
};

struct dispersion {

  int nvel;
  double *rad;
  double *modeldisp;
  double *disp;
  double *disperr;

};

struct sz {

  int nsz;
  double *szx;
  double *szy;
  double *sze;
  double *modely;

  /* Siegel Begin */
  double *rot_szy;
  double *rot_sze;
  double *rot_modely;

  double *eigenval;
  double **eigenvec;
  double **inv_eigenvec;
  /* Siegel End */

};

struct fitdata {

  struct pha *phalist;
  struct rmf *rmflist;
  struct arf *arflist;
  struct ebounds *elimits;
  struct weaklensing *wldata;
  struct dispersion *vels;
  struct sz *szdata;
  struct jaco_state *js;
  double *spectrum;
  double syserr;
  double Rcut[NINSTR];
  int    fitannuli;
  int    sz0;

  long ndata;
  double *bigdata;
  double *bigerr;
  double *bigmodel;

};

int get_bigmodel(double *x,double *params,double *model, double *errors,
		 double *prior,
		     unsigned long ndata, void *datapointer);
int init_standalone(struct fitdata *data);

/* Siegel Begin */
int read_lensing_covariance_and_rotate_data(const char *fits_file, struct weaklensing *wldata);
int rotate_lensing_model(struct weaklensing *wldata);
int read_sz_covariance_and_rotate_data(const char *fits_file, struct sz *szdata);
int rotate_sz_model(struct sz *szdata);
/* Siegel End */