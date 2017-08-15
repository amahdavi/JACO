struct weaklensing {

  int nwl;
  double *wlrad;
  double *modelshear;
  double *shear;
  double *shearerr;

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
