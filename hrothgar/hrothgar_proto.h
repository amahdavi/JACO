#include <gsl/gsl_chebyshev.h>
#include <gsl/gsl_blas.h>
#include <gsl/gsl_min.h>
#include <gsl/gsl_rng.h>
#include <gsl/gsl_sort.h>
#include <gsl/gsl_randist.h>
#include <gsl/gsl_multifit_nlin.h>
#include <gsl/gsl_statistics.h>
#include <gsl/gsl_sf_gamma.h>
#include <gsl/gsl_interp.h>

/* hrothgar.c */
double *hrothgar_dvector(long nelem);
size_t *hrothgar_sizetvector(long nelem);
double **hrothgar_dmatrix(long nx, long ny);
void hrothgar_free_dmatrix(double **x, long nx);
double ***hrothgar_dtensor(long nx, long ny, long nz);
void hrothgar_free_dtensor(double ***x, long nx, long ny);
void percentdone(long n, long tot);
int mcmc_test(double statnew, double statold, double *params,
	      struct hrothgar_setup *setup);
void mcmc_step(double *newparams, double *oldparams, struct hrothgar_setup *setup);
int minimize_func_func(const gsl_vector *p, void *extra, gsl_vector *f);
void multiproc_onecpu(double **input, double **output, int nunits, struct hrothgar_setup *setup);
int minimize_func_funcder(const gsl_vector *p, void *extra, gsl_vector *f, gsl_matrix *J);
int precision_gradient(const gsl_vector *p, void *extra, gsl_vector *f, gsl_matrix *J);
int precision_gradient2(const gsl_vector *p, void *extra, gsl_vector *f, gsl_matrix *J);
int minimize_func_der(const gsl_vector *p, void *extra, gsl_matrix *J);
double unmap_par(double x, double min, double max);
double map_par(double x, double min, double max);
void print_state(FILE *of, const char *prefix, double chi2, gsl_vector *s, struct hrothgar_setup *setup);
void init_minimize_workspace(struct hrothgar_setup *setup);
void free_minimize_workspace(struct hrothgar_setup *setup);
void confidence_2d(double cinfo[], double lims[], struct hrothgar_setup *setup);
void confidence_looplimits(struct hrothgar_setup *setup);
void initialize_confidence(double *cinfo, int i, int j, struct hrothgar_setup *setup);
void confidence_report(int i, double info[], struct hrothgar_setup *setup);
void confoutput_standalone(struct hrothgar_setup *setup);
void write_minimize_output(gsl_vector *p, struct hrothgar_setup *setup);
double linemin_function(double linpar, void *minparams);
void line_minimize(struct hrothgar_setup *s, gsl_min_fminimizer *minimizer, double *linpar, double *minchi);
void lord_minimize(struct hrothgar_setup *setup);
void report_unmapped(double *unmapped_pars, double chisq, 
		     struct hrothgar_setup *setup);
void hrothgar_mcmc(int initial, long nsim, double *mcmcdata, struct hrothgar_setup *setup);
void hrothgar_singlecpu(struct hrothgar_setup *setup);
/* hrothgar_image.c */
void hrothgar_writeimage(char *outfile, float *buffer, double x1, double x2, double y1, double y2, long ngrid);
void hrothgar_writeimage_double(char *outfile,
				double **dbufffer,
				double x1, double x2, double y1, double y2,
				long ngrid);
void hrothgar_wcsimage(char *outfile,
		       float *buffer, 
		       float wcs[3][2],
		       long ngrid);
/* hrothgar_mpi.c */
void log_mpi_node(struct treenode *p, struct hrothgar_setup *setup);
struct treenode *allocate_node(int ntotparams);
void delete_node(struct treenode *delnode);
struct treenode *find_shallowest(struct treenode *p);
void delete_subtree(struct treenode *p, struct treenode **cpunodes);
void update_tree(struct treenode *p, int depth);
int traverse_tree(struct treenode **root, struct treenode **cpunodes, struct hrothgar_setup *setup);
struct treenode *allocate_work(int cpu, struct treenode *root, struct hrothgar_setup *setup);
void lord_mcmc(struct hrothgar_setup *setup);
void thane_mcmc(struct hrothgar_setup *setup);
void lord_powell(struct hrothgar_setup *setup, gsl_vector *p, 
		 double *chimin, double globchimin);
void thane_minimize(struct hrothgar_setup *setup);
void multiproc_mpi(double **input, double **output, int nunits, struct hrothgar_setup *setup);
void lord_confoutput(struct hrothgar_setup *setup);
void thane_confoutput(struct hrothgar_setup *setup);
void reset_all_nodes(struct hrothgar_setup *setup);
void hrothgar_mpi(struct hrothgar_setup *setup);
void init_mpi(int *ntotcpu, int *rank);
void adjust_limits(double **limits, int i, int j, double *lims, int sense);
int stepmatrix_from_covariance(struct hrothgar_setup *setup);
