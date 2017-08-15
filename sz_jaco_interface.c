#pragma warning(disable: 4996)
#include "sz_jaco_interface.h"
#include<stdlib.h>
#include<stdio.h>
#include<string.h>
#include<math.h>
#include "fitsio.h"
#include "longnam.h"
#include "fftw3.h"

/* constant variable declaration */

static const double DPI = 3.14159265358979;
static const double SIGMA_TO_FWHM = 2.354820045030949;
static const double ARCSEC_TO_DEG = 1.0/3600.0;
static const double ARCMIN_TO_DEG = 1.0/60.0;
static const double DEG_TO_RAD = 3.14159265358979/180.0;

static const double MAX_RES = 10.0/3600.0;
static const double SIGMA_POINTING = 5.0/3600.0;
static const double MIN_MODEL_SIZE = 25.0/60.0;


/* functions */

int big_get_sz_model(double *models,struct BigSZControlParams *params)
{
	int ds, npix_data, npix_model, nshift, npad, navg, nblock, N,
		ra, dec, ra_start, dec_start, ii, jj, kk, index_start, index_stop;
	double npix_data_sqrd, npix_model_sqrd, real_component, imag_component, yavg, delta_T;

	/* Fill in nest */

	int nest = params->nest;

	/* Set data, errors, and models */

	double *data = params->data;
	double *errors = params->errors;

	/*  This requires a function call of the form:
			big_get_sz_model(int *nest, double **data, double **errors, double **models, ...)
		In this approach, you do not need to copy the entire array from one
		location to another.  Instead you simply copy the address where the
		calling function can find the data and errors.  The calling function
		issues something like this:
			int nest
			double *data;
			double *errors;

			big_get_sz_model(&nest, &data, &model, ...)
	*/


	/* Loop through datasets */

	index_start = 0;

	for (ds=0; ds<params->ndataset; ds++)
	{
		npix_data = params->npix_data[ds];
		npix_model = params->npix[ds];

		npix_data_sqrd = ((double) npix_data)*((double) npix_data);
		npix_model_sqrd = 4.0*((double) npix_model)*((double) npix_model);

		N = params->nest_per_dataset[ds];

		index_stop = index_start + N;

		nshift = (npix_model-1)/2;

		nblock = params->scaling_factor[ds];
		navg = nblock*nblock;
		npad = (npix_model - nblock*npix_data)/2;

		ra_start = nshift + npix_model - npad - 1;
		dec_start = nshift + npad;

		/* Transfer params->ymap to params->fft_model */

		for (ra=0; ra<npix_model; ra++)
		{
			for (dec=0; dec<npix_model; dec++)
			{
				params->fft_model[ds][ra][dec][0] = params->ymap[ds][ra][dec];
				params->fft_model[ds][ra][dec][1] = 0.0;
			}
		}

		/* Zero-pad */

		zero_pad_sz(params->fft_model[ds], npix_model, 2*npix_model, npix_model, 2*npix_model);

		/* Execute fft_model_forward */

		fftw_execute(params->fft_model_plan_forward[ds]);

		/* Convolve the model ymap with the PSF and pointing uncertainty */

		for (ra=0; ra<2*npix_model; ra++)
		{
			for (dec=0; dec<2*npix_model; dec++)
			{
				real_component = (params->fft_model[ds][ra][dec][0])*(params->fft_beam[ds][ra][dec][0]) - (params->fft_model[ds][ra][dec][1])*(params->fft_beam[ds][ra][dec][1]);
				imag_component = (params->fft_model[ds][ra][dec][0])*(params->fft_beam[ds][ra][dec][1]) + (params->fft_model[ds][ra][dec][1])*(params->fft_beam[ds][ra][dec][0]);
				params->fft_model[ds][ra][dec][0] = real_component/npix_model_sqrd;
				params->fft_model[ds][ra][dec][1] = imag_component/npix_model_sqrd;
			}
		}

		/* Execute fft_model_backward */

		fftw_execute(params->fft_model_plan_backward[ds]);

		/* Discard outskirts of model ymap and rebin */
		/* Also correct for the shift that occured during convolution */

		for (ii=0; ii<npix_data; ii++)
		{
			for (jj=0; jj<npix_data; jj++)
			{
				yavg = 0.0;
				
				for (ra=0; ra<nblock; ra++)
				{
					for (dec=0; dec<nblock; dec++)
					{
						yavg += params->fft_model[ds][(ra_start - (ii*nblock + ra))%(2*npix_model)][(dec_start + (jj*nblock + dec))%(2*npix_model)][0];
					}
				}

				yavg /= navg;

				params->fft_data[ds][ii][jj][0] = yavg;
				params->fft_data[ds][ii][jj][1] = 0.0;

			}
		}

		/* Execute fft_data_forward */

		fftw_execute(params->fft_data_plan_forward[ds]);

		/* Convolve the map with the transfer function */

		for (ii=0; ii<npix_data; ii++)
		{
			for (jj=0; jj<npix_data; jj++)
			{
				real_component = (params->fft_data[ds][ii][jj][0])*(params->signal_transfer_function[ds][ii][jj][0]) - 
								 (params->fft_data[ds][ii][jj][1])*(params->signal_transfer_function[ds][ii][jj][1]);
				imag_component = (params->fft_data[ds][ii][jj][0])*(params->signal_transfer_function[ds][ii][jj][1]) + 
								 (params->fft_data[ds][ii][jj][1])*(params->signal_transfer_function[ds][ii][jj][0]);
				params->fft_data[ds][ii][jj][0] = real_component/npix_data_sqrd;
				params->fft_data[ds][ii][jj][1] = imag_component/npix_data_sqrd;
			}
		}

		/* Execute fft_data_backward */

		fftw_execute(params->fft_data_plan_backward[ds]);

		/* Transfer to the model vector */

		for (kk=index_start; kk<index_stop; kk++)
		{
			(models)[kk] = (*(params->coverage_mask[kk]))*(params->unit_conversion[ds]);
		}

		/* Determine least squares estimate of the DC offset */

		delta_T = 0;

		for (kk=index_start; kk<index_stop; kk++)
		{
			delta_T += ((data)[kk] - (models)[kk])/((errors)[kk]*(errors)[kk]);
		}

		delta_T /= params->sum_inverse_variance[ds];

		/* Add this offset to the model */
		
		for (kk=index_start; kk<index_stop; kk++)
		{
			(models)[kk] += delta_T;
		}
		
		params->delta_T[ds] = delta_T;

		/* Update the start index */

		index_start += N;

	}

	return(1);

}

int big_initialize_sz(struct BigSZControlParams *ptr_to_szstruct)
{
	const char *current_paramfile;
	double *data, *errors;
	double ra0, dec0, rmax, temp_err;
	struct BigSZKeywords data_keywords, model_keywords;
	int ds, noerr, ndataset, npts, index_begin, ra, dec, count, ii;
	fftw_plan beam_plan;
	fftw_complex **beam;
	
	if (ptr_to_szstruct == NULL)
	{
		fprintf(stderr,"input ptr_to_szstruct is NULL\n");
		return(0);
	}

	ptr_to_szstruct->nest = 0;

	ndataset = ptr_to_szstruct->ndataset;

	ptr_to_szstruct->nest_per_dataset = (int *) malloc(ndataset*sizeof(int));
	if (ptr_to_szstruct->nest_per_dataset == NULL)
	{
		fprintf(stderr,"ran out of memory during allocation of szstruct.nest_per_dataset\n");
		return(0);
	}

	ptr_to_szstruct->rmax = (double *) malloc(ndataset*sizeof(double));
	if (ptr_to_szstruct->rmax == NULL)
	{
		fprintf(stderr,"ran out of memory during allocation of szstruct.rmax\n");
		return(0);
	}

	/* In order to allocate vectors that store the data and errors, we need to know the total number of pixels over all of the data sets.
	This requires a preliminary loop over the datasets.  We also check that all of the fits images are in order during this step. */

	for (ds=0; ds<ndataset; ds++)
	{
		current_paramfile = (const char*) ptr_to_szstruct->paramfiles[ds];

		noerr = error_check_bolocam_fits_file(current_paramfile, &npts, &ptr_to_szstruct->rmax[ds]);
		if (!noerr)
		{
			fprintf(stderr,"failure during error_check_bolocam_fits_file of paramsfile number %i of %i\n",
							ds+1, ndataset);
			return(0);
		}

		ptr_to_szstruct->nest_per_dataset[ds] = npts;

		ptr_to_szstruct->nest += npts;
	}


	/* Allocate space for decent, racent, pixsize, and npix */
	
	ptr_to_szstruct->deccent = (double *) malloc(ndataset*sizeof(double));
	if (ptr_to_szstruct->deccent == NULL)
	{
		fprintf(stderr,"ran out of memory during allocation of szstruct.deccent\n");
		return(0);
	}

	ptr_to_szstruct->racent = (double *) malloc(ndataset*sizeof(double));
	if (ptr_to_szstruct->racent == NULL)
	{
		fprintf(stderr,"ran out of memory during allocation of szstruct.raccent\n");
		return(0);
	}

	ptr_to_szstruct->pixsize = (double *) malloc(ndataset*sizeof(double));
	if (ptr_to_szstruct->pixsize == NULL)
	{
		fprintf(stderr,"out of memory during allocation of szstruct.pixsize\n");
		return(0);
	}

	ptr_to_szstruct->npix = (int *) malloc(ndataset*sizeof(int));
	if (ptr_to_szstruct->npix == NULL)
	{
		fprintf(stderr,"out of memory during allocation of szstruct.npix\n");
		return(0);
	}

	ptr_to_szstruct->npix_data = (int *) malloc(ndataset*sizeof(int));
	if (ptr_to_szstruct->npix_data == NULL)
	{
		fprintf(stderr,"out of memory during allocation of szstruct.npix_data\n");
		return(0);
	}

	ptr_to_szstruct->scaling_factor = (int *) malloc(ndataset*sizeof(int));
	if (ptr_to_szstruct->scaling_factor == NULL)
	{
		fprintf(stderr,"out of memory during allocation of szstruct.scaling_factor\n");
		return(0);
	}

	ptr_to_szstruct->unit_conversion = (double *) malloc(ndataset*sizeof(double));
	if (ptr_to_szstruct->unit_conversion == NULL)
	{
		fprintf(stderr,"out of memory during allocation of szstruct.unit_conversion\n");
		return(0);
	}

	ptr_to_szstruct->sum_inverse_variance = (double *) malloc(ndataset*sizeof(double));
	if (ptr_to_szstruct->sum_inverse_variance == NULL)
	{
		fprintf(stderr,"out of memory during allocation of szstruct.sum_inverse_variance\n");
		return(0);
	}

	ptr_to_szstruct->delta_T = (double *) malloc(ndataset*sizeof(double));
	if (ptr_to_szstruct->delta_T == NULL)
	{
		fprintf(stderr,"out of memory during allocation of szstruct.delta_T\n");
		return(0);
	}

	/* Allocate space for the fftw plans */

	ptr_to_szstruct->fft_model_plan_forward = (fftw_plan *) malloc(ndataset*sizeof(fftw_plan));
	if (ptr_to_szstruct->fft_model_plan_forward == NULL)
	{
		fprintf(stderr,"out of memory during allocation of szstruct.fft_model_plan_forward\n");
		return(0);
	}

	ptr_to_szstruct->fft_model_plan_backward = (fftw_plan *) malloc(ndataset*sizeof(fftw_plan));
	if (ptr_to_szstruct->fft_model_plan_backward == NULL)
	{
		fprintf(stderr,"out of memory during allocation of szstruct.fft_model_plan_backward\n");
		return(0);
	}

	ptr_to_szstruct->fft_data_plan_forward = (fftw_plan *) malloc(ndataset*sizeof(fftw_plan));
	if (ptr_to_szstruct->fft_data_plan_forward == NULL)
	{
		fprintf(stderr,"out of memory during allocation of szstruct.fft_data_plan_forward\n");
		return(0);
	}

	ptr_to_szstruct->fft_data_plan_backward = (fftw_plan *) malloc(ndataset*sizeof(fftw_plan));
	if (ptr_to_szstruct->fft_data_plan_backward == NULL)
	{
		fprintf(stderr,"out of memory during allocation of szstruct.fft_data_plan_backward\n");
		return(0);
	}


	/* Allocate the vector of pointers to the arrays that will be used for performing fft's */

	ptr_to_szstruct->fft_model = (fftw_complex ***) malloc(ndataset*sizeof(fftw_complex **));
	if (ptr_to_szstruct->fft_model == NULL)
	{
		fprintf(stderr,"out of memory during allocation of szstruct.fft_model\n");
		return(0);
	}

	ptr_to_szstruct->fft_data = (fftw_complex ***) malloc(ndataset*sizeof(fftw_complex **));
	if (ptr_to_szstruct->fft_data == NULL)
	{
		fprintf(stderr,"out of memory during allocation of szstruct.fft_data\n");
		return(0);
	}

	ptr_to_szstruct->fft_beam = (fftw_complex ***) malloc(ndataset*sizeof(fftw_complex **));
	if (ptr_to_szstruct->fft_beam == NULL)
	{
		fprintf(stderr,"out of memory during allocation of szstruct.fft_beam\n");
		return(0);
	}

	ptr_to_szstruct->signal_transfer_function = (fftw_complex ***) malloc(ndataset*sizeof(fftw_complex **));
	if (ptr_to_szstruct->fft_beam == NULL)
	{
		fprintf(stderr,"out of memory during allocation of szstruct.signal_transfer_function\n");
		return(0);
	}

	/* Allocate space for data, errors, models, and coverage mask */

	ptr_to_szstruct->data = (double *) malloc(ptr_to_szstruct->nest*sizeof(double));
	if (ptr_to_szstruct->data == NULL)
	{
		fprintf(stderr,"ran out of memory during allocation of szstruct.data\n");
		return(0);
	}

	ptr_to_szstruct->errors = (double *) malloc(ptr_to_szstruct->nest*sizeof(double));
	if (ptr_to_szstruct->errors == NULL)
	{
		fprintf(stderr,"ran out of memory during allocation of szstruct.errors\n");
		return(0);
	}

	ptr_to_szstruct->models = (double *) malloc(ptr_to_szstruct->nest*sizeof(double));
	if (ptr_to_szstruct->models == NULL)
	{
		fprintf(stderr,"ran out of memory during allocation of szstruct.models\n");
		return(0);
	}

	ptr_to_szstruct->coverage_mask = (double **) malloc(ptr_to_szstruct->nest*sizeof(double *));
	if (ptr_to_szstruct->coverage_mask == NULL)
	{
		fprintf(stderr,"ran out of memory during allocation of szstruct.coverage_mask\n");
		return(0);
	}


	/* Zero model vector */

	for (ii=0; ii<ptr_to_szstruct->nest; ii++)
	{
		ptr_to_szstruct->models[ii] = 0.0;
	}


	/* Loop through input files */

	index_begin = 0;

	for (ds=0; ds<ndataset; ds++)
	{
		/* Read in keywords */

		current_paramfile = (const char*) ptr_to_szstruct->paramfiles[ds];

		noerr = read_bolocam_fits_file_keywords(current_paramfile, &data_keywords);
		if (!noerr)
		{
			fprintf(stderr,"read_bolocam_fits_file_keywords error for paramsfile number %i of %i\n",
							ds+1, ndataset);
			return(0);
		}


		data_keywords.rmax = ptr_to_szstruct->rmax[ds];

		/* Determine the desired range and resolution of the output model y-map */

		model_keywords = data_keywords;

		define_sz_model_map(&model_keywords);


		/* Save values to szstruct */

		ptr_to_szstruct->npix_data[ds] = data_keywords.npix;
		ptr_to_szstruct->racent[ds] = model_keywords.crval1;
		ptr_to_szstruct->deccent[ds] = model_keywords.crval2;
		ptr_to_szstruct->pixsize[ds] = model_keywords.pixsize;
		ptr_to_szstruct->npix[ds] = model_keywords.npix;
		ptr_to_szstruct->scaling_factor[ds] = model_keywords.scaling_factor;
		ptr_to_szstruct->unit_conversion[ds] = 1.0/model_keywords.unit_con;


		/* Allocate memory for the model, data, and beam fft */

		noerr = allocate_fftw_complex_matrix_contiguous(&(ptr_to_szstruct->fft_model[ds]), 2*model_keywords.npix, 2*model_keywords.npix);
		if (!noerr)
		{
			fprintf(stderr,"allocate_fftw_complex_matrix error for paramsfile number %i of %i\n",
							ds+1, ndataset);
			return(0);
		}

		noerr = allocate_fftw_complex_matrix_contiguous(&(ptr_to_szstruct->fft_beam[ds]), 2*model_keywords.npix, 2*model_keywords.npix);
		if (!noerr)
		{
			fprintf(stderr,"allocate_fftw_complex_matrix error for paramsfile number %i of %i\n",
							ds+1, ndataset);
			return(0);
		}

		noerr = allocate_fftw_complex_matrix_contiguous(&(ptr_to_szstruct->fft_data[ds]), data_keywords.npix, data_keywords.npix);
		if (!noerr)
		{
			fprintf(stderr,"allocate_fftw_complex_matrix error for paramsfile number %i of %i\n",
							ds+1, ndataset);
			return(0);
		}


		/* Create plans now for the ffts that will be performed in big_get_sz_model */

		ptr_to_szstruct->fft_model_plan_forward[ds] =  fftw_plan_dft_2d(2*model_keywords.npix, 2*model_keywords.npix,
																		 &(ptr_to_szstruct->fft_model[ds][0][0]), 
																		 &(ptr_to_szstruct->fft_model[ds][0][0]),
																		 FFTW_FORWARD, FFTW_MEASURE);

		ptr_to_szstruct->fft_model_plan_backward[ds] = fftw_plan_dft_2d(2*model_keywords.npix, 2*model_keywords.npix,
																		 &(ptr_to_szstruct->fft_model[ds][0][0]), 
																		 &(ptr_to_szstruct->fft_model[ds][0][0]),
																		 FFTW_BACKWARD, FFTW_MEASURE);

		ptr_to_szstruct->fft_data_plan_forward[ds] =  fftw_plan_dft_2d(data_keywords.npix, data_keywords.npix,
																		 &(ptr_to_szstruct->fft_data[ds][0][0]), 
																		 &(ptr_to_szstruct->fft_data[ds][0][0]),
																		 FFTW_FORWARD, FFTW_MEASURE);

		ptr_to_szstruct->fft_data_plan_backward[ds] = fftw_plan_dft_2d(data_keywords.npix, data_keywords.npix,
																		 &(ptr_to_szstruct->fft_data[ds][0][0]), 
																		 &(ptr_to_szstruct->fft_data[ds][0][0]),
																		 FFTW_BACKWARD, FFTW_MEASURE);


		/* Generate and fourier transform a beam kernel that can be convolved with the output model y-map */

		noerr = allocate_fftw_complex_matrix_contiguous(&beam, 2*model_keywords.npix, 2*model_keywords.npix);
		if (!noerr)
		{
			fprintf(stderr,"allocate_fftw_complex_matrix error for paramsfile number %i of %i\n",
							ds+1, ndataset);
			return(0);
		}

		beam_plan = fftw_plan_dft_2d(2*model_keywords.npix, 2*model_keywords.npix,
									 &beam[0][0], &(ptr_to_szstruct->fft_beam[ds][0][0]), 
									 FFTW_FORWARD, FFTW_MEASURE);

		noerr = generate_bolocam_beam_kernel(beam, model_keywords);

		fftw_execute(beam_plan);

		fftw_destroy_plan(beam_plan);
		free_fftw_complex_matrix_contiguous(&beam);


		/* Read in data, errors, and the transfer function */

		noerr = allocate_fftw_complex_matrix_contiguous(&(ptr_to_szstruct->signal_transfer_function[ds]), data_keywords.npix, data_keywords.npix);
		if (!noerr)
		{
			fprintf(stderr,"allocate_fftw_complex_matrix_contiguous error for paramsfile number %i of %i\n",
							ds+1, ndataset);
			return(0);
		}

		data = (double *) malloc(data_keywords.npix*data_keywords.npix*sizeof(double));
		errors = (double *) malloc(data_keywords.npix*data_keywords.npix*sizeof(double));

		noerr = read_bolocam_fits_file(current_paramfile, data, errors, ptr_to_szstruct->signal_transfer_function[ds], data_keywords.npix);


		/* Discard data points outside of aperture */

		ra0 = data_keywords.crpix1 - 1.0;
		dec0 = data_keywords.crpix2 - 1.0;

		rmax = data_keywords.rmax/data_keywords.pixsize;

		count = 0;

		ptr_to_szstruct->sum_inverse_variance[ds] = 0.0;

		if (ptr_to_szstruct->nest_per_dataset[ds] < data_keywords.npix*data_keywords.npix)
		{
			for (ra=0; ra < data_keywords.npix; ra++)
			{
				for (dec=0; dec < data_keywords.npix; dec++)
				{
					if (((ra-ra0)*(ra-ra0) + (dec-dec0)*(dec-dec0)) <= rmax*rmax)
					{
						temp_err = errors[ra*data_keywords.npix + dec];

						ptr_to_szstruct->data[index_begin + count] = data[ra*data_keywords.npix + dec];
						ptr_to_szstruct->errors[index_begin + count] = temp_err;

						ptr_to_szstruct->sum_inverse_variance[ds] += 1.0/(temp_err*temp_err);

						ptr_to_szstruct->coverage_mask[index_begin + count] = &(ptr_to_szstruct->fft_data[ds][ra][dec][0]);

						count++;
					}
				}
			}
		}
		else
		{

			memcpy(&(ptr_to_szstruct->data[index_begin]), data, data_keywords.npix*data_keywords.npix*sizeof(double));
			memcpy(&(ptr_to_szstruct->errors[index_begin]), errors, data_keywords.npix*data_keywords.npix*sizeof(double));

			for (ra=0; ra < data_keywords.npix; ra++)
			{
				for (dec=0; dec < data_keywords.npix; dec++)
				{
					temp_err = errors[ra*data_keywords.npix + dec];

					ptr_to_szstruct->sum_inverse_variance[ds] += 1.0/(temp_err*temp_err);

					ptr_to_szstruct->coverage_mask[index_begin + ra*data_keywords.npix + dec] = &(ptr_to_szstruct->fft_data[ds][ra][dec][0]);
				}
			}

		}

		free(data);
		free(errors);

		index_begin += ptr_to_szstruct->nest_per_dataset[ds];

		if (!ptr_to_szstruct->quiet)
		{
			printf("npix*npix:    %li\n",data_keywords.npix*data_keywords.npix);
			printf("nest     :	  %li\n",ptr_to_szstruct->nest_per_dataset[ds]);
		}

	}

	return(1);
}

int free_big_sz(struct BigSZControlParams *psz)
{
	int ndataset, ds;

	ndataset = psz->ndataset;

	free(psz->nest_per_dataset);
	free(psz->deccent);
	free(psz->racent);
	free(psz->pixsize);
	free(psz->npix);
	free(psz->npix_data);
	free(psz->rmax);
	free(psz->scaling_factor);
	free(psz->unit_conversion);
	free(psz->sum_inverse_variance);
	free(psz->delta_T);

	free(psz->data);
	free(psz->errors);
	free(psz->models);
	free(psz->coverage_mask);

	for (ds=0; ds<ndataset; ds++)
	{
		fftw_destroy_plan(psz->fft_model_plan_backward[ds]);
		fftw_destroy_plan(psz->fft_model_plan_forward[ds]);
		fftw_destroy_plan(psz->fft_data_plan_backward[ds]);
		fftw_destroy_plan(psz->fft_data_plan_forward[ds]);

		free_fftw_complex_matrix_contiguous(&psz->fft_model[ds]);
		free_fftw_complex_matrix_contiguous(&psz->fft_beam[ds]);
		free_fftw_complex_matrix_contiguous(&psz->fft_data[ds]);
		free_fftw_complex_matrix_contiguous(&psz->signal_transfer_function[ds]);
	}

	free(psz->fft_model_plan_backward);
	free(psz->fft_model_plan_forward);
	free(psz->fft_data_plan_backward);
	free(psz->fft_data_plan_forward);

	free(psz->fft_model);
	free(psz->fft_beam);
	free(psz->fft_data);
	free(psz->signal_transfer_function);

	return(1);
}

int error_check_bolocam_fits_file(const char *fits_file, int *ndata, double *rmax)
{
	fitsfile *fptr;
	int naxis, naxis1, naxis2, next, iext, cnaxis, cbitpix;
	long cN[2];
	double bmaj, bmin;
	char ctype1[9], ctype2[9];
	int status = 0, close = 0;

	int ninside;
	double pixsize, x0, y0, rm;
	double *errors;

	*ndata = 0;

	if (fits_file == NULL)
	{
		fprintf(stderr,"input pointer to the fits_file is NULL\n");
		return(0);
	}

	/* Open the fits file */

	fits_open_file(&fptr, fits_file, READONLY, &status);

	/* Read in relevant keywords */

	fits_read_key(fptr, TINT, "NAXIS", &naxis, NULL, &status);
	fits_read_key(fptr, TINT, "NAXIS1", &naxis1, NULL, &status);
	fits_read_key(fptr, TINT, "NAXIS2", &naxis2, NULL, &status);
	fits_read_key(fptr, TDOUBLE, "BMAJ", &bmaj, NULL, &status);
	fits_read_key(fptr, TDOUBLE, "BMIN", &bmin, NULL, &status);
	fits_read_key(fptr, TSTRING, "CTYPE1", &ctype1, NULL, &status);
	fits_read_key(fptr, TSTRING, "CTYPE2", &ctype2, NULL, &status);

	ctype1[8] = '\0';
	ctype2[8] = '\0';

	/* Error check */
	if (status)
	{
		fits_report_error(stderr, status);
		fits_close_file(fptr, &close);
		return(0);
	}
	if (naxis != 2)
	{
		fprintf(stderr,"sz image must be two-dimensional\n");
		fits_close_file(fptr, &close);
		return(0);
	}
	if (naxis1 != naxis2)
	{
		fprintf(stderr,"sz code currently only supports square images\n");
		fits_close_file(fptr, &close);
		return(0);
	}
	if (bmaj != bmin)
	{
		fprintf(stderr,"sz code currently only supports symettric PSFs\n");
		fits_close_file(fptr, &close);
		return(0);
	}
	if (strncmp(ctype1, "RA", 2) || strncmp(ctype2, "DEC", 3))
	{
		fprintf(stderr,"first dimension of sz image must be RA, second dimension must be DEC\n");
		fits_close_file(fptr, &close);
		return(0);
	}


	/* Loop through main image and extensions, check that all images are the same size */

	fits_get_num_hdus(fptr, &next, &status);

	if (status)
	{
		fits_report_error(stderr, status);
		fits_close_file(fptr, &close);
		return(0);
	}
	if (next != 4)
	{
		fprintf(stderr,"fits image must contain 4 extensions:  Map, Map RMS, Re Trans, Im Trans\n");
		fits_close_file(fptr, &close);
		return(0);
	}

	errors = (double *) malloc(naxis1*naxis2*sizeof(double));

	for (iext=1; iext<=next; iext++)
	{
		fits_movabs_hdu(fptr, iext, NULL, &status);

		fits_get_img_param(fptr, 2, &cbitpix, &cnaxis, cN, &status);

		if (status)
		{
			fits_report_error(stderr, status);
			fits_close_file(fptr, &close);
			return(0);
		}

		/* Make sure the image is two dimensions of type double */

		if (cnaxis != 2)
		{
			fprintf(stderr,"fits extension %i does not contain 2 axis\n", iext);
			fits_close_file(fptr, &close);
			return(0);
		}
		if ((cN[0] != naxis1) || (cN[1] != naxis2))
		{
			fprintf(stderr,"fits extension %i does not match dimensions specified in header keywords\n", iext);
			fits_close_file(fptr, &close);
			return(0);
		}
		if (cbitpix != -64)
		{
			fprintf(stderr,"fits extension %i is not of type double\n", iext);
			fits_close_file(fptr, &close);
			return(0);
		}

	}

	fits_movabs_hdu(fptr, 1, NULL, &status);
	
	if (status)
	{
		fits_report_error(stderr, status);
		fits_close_file(fptr, &close);
		return(0);
	}


	/* Determine if an aperture needs to be applied */

	fits_read_key(fptr, TDOUBLE, "RMAX", &rm, NULL, &status);

	if (!status)
	{
		if (rm)
		{
			fits_read_key(fptr, TDOUBLE, "PIXSIZE", &pixsize, NULL, &status);
			fits_read_key(fptr, TDOUBLE, "CRPIX1", &x0, NULL, &status);
			fits_read_key(fptr, TDOUBLE, "CRPIX2", &y0, NULL, &status);

			if (status)
			{
				fits_report_error(stderr, status);
				fits_close_file(fptr, &close);
				return(0);
			}

			x0 += - 1.0;
			y0 += - 1.0;

			rm /= pixsize;

			ninside = count_pixels_within_aperture(naxis1, naxis2, x0, y0, rm);

			if (!ninside)
			{
				fprintf(stderr,"check inputs to count_pixels_within_aperture\n");
				fits_close_file(fptr, &close);
				return(0);
			}

			*ndata = ninside;
			*rmax = rm*pixsize*ARCSEC_TO_DEG;
		}
		else
		{
			*ndata = naxis1*naxis1;				/* RMAX = 0 means that no aperture will be applied */
			*rmax = 0.0;
		}
	}
	else
	{
		if (status == 202)						/* Status indicating that the keyword RMAX does not exist */
		{
			*ndata = naxis1*naxis1;
			*rmax = 0.0;
		}
		else
		{
			fits_report_error(stderr, status);
			fits_close_file(fptr, &close);
			return(0);
		}
	} 


	/* Close fits file */

	fits_close_file(fptr, &close);

	return(1);
}

int read_bolocam_fits_file(const char *fits_file, double *data, double *errors, fftw_complex **transfer_function, int npix)
{
	fitsfile *fptr;
	double **img, **rms, **rez, **imz;
	int status=0, anynul=0, noerr, nfill, iext, ra, dec;
	long fpixel[2] = {1, 1};

	if (fits_file == NULL)
	{
		fprintf(stderr,"input pointer to the fits_file is NULL\n");
		return(0);
	}

	/* Open the fits file */

	fits_open_file(&fptr, fits_file, READONLY, &status);

	fits_get_hdu_num(fptr, &iext);
	if (iext != 1)
	{
		fits_movabs_hdu(fptr, 1, NULL, &status);
	}

	/* Read in the image */

	nfill = npix*npix;

	noerr = allocate_double_matrix_contiguous(&img, npix, npix);
	if (!noerr)
	{
		fprintf(stderr,"out of memory during allocation of temporary image matrix\n");
		return(0);
	}

	fits_read_pix(fptr, TDOUBLE, fpixel, nfill, NULL, &(img[0][0]), &anynul, &status);
	if (status)
	{
		fits_report_error(stderr, status);
		status = 0;
		fits_close_file(fptr, &status);
		return(0);
	}

	/* Place the image in the data vector, which unfortunately requires a transposition.  This is because C stores arrays as Arr[row][column]
	 -- with elements sharing the same row contiguous in memory --  so when CFITSIO reads in our image from the fits file 
	it stores it as img[dec][ra].  However, we will eventually want to compare this image with the model ymap.  The model ymap
	is stored as ymap[ra][dec].  So we either need to transpose the data now, or the model ymap at each iteration of the 
	Markov Chain. */

	for (dec=0; dec<npix; dec++)
	{
		for (ra=0; ra<npix; ra++)
		{
			data[ra*npix + dec] = img[dec][ra];
		}
	}

	free_double_matrix_contiguous(&img);


	/* Read in the errors */

	noerr = allocate_double_matrix_contiguous(&rms, npix, npix);
	if (!noerr)
	{
		fprintf(stderr,"out of memory during allocation of temporary rms matrix\n");
		return(0);
	}

	fits_movabs_hdu(fptr, 2, NULL, &status);
	fits_read_pix(fptr, TDOUBLE, fpixel, nfill, NULL, &(rms[0][0]), &anynul, &status);
	if (status)
	{
		fits_report_error(stderr, status);
		status = 0;
		fits_close_file(fptr, &status);
		return(0);
	}

	/* Store the rms matrix in the error vector */

	for (dec=0; dec<npix; dec++)
	{
		for (ra=0; ra<npix; ra++)
		{
			errors[ra*npix + dec] = rms[dec][ra];
		}
	}

	free_double_matrix_contiguous(&rms);


	/* Read in the real component of the transfer function */

	noerr = allocate_double_matrix_contiguous(&rez, npix, npix);
	if (!noerr)
	{
		fprintf(stderr,"out of memory during allocation of temporary real signal transfer function matrix\n");
		return(0);
	}

	fits_movabs_hdu(fptr, 3, NULL, &status);
	fits_read_pix(fptr, TDOUBLE, fpixel, nfill, NULL, &(rez[0][0]), &anynul, &status);
	if (status)
	{
		fits_report_error(stderr, status);
		status = 0;
		fits_close_file(fptr, &status);
		return(0);
	}

	/* Read in the imaginary component of the transfer function */

	noerr = allocate_double_matrix_contiguous(&imz, npix, npix);
	if (!noerr)
	{
		fprintf(stderr,"out of memory during allocation of temporary imaginary signal transfer function matrix\n");
		return(0);
	}

	fits_movabs_hdu(fptr, 4, NULL, &status);
	fits_read_pix(fptr, TDOUBLE, fpixel, nfill, NULL, &(imz[0][0]), &anynul, &status);
	if (status)
	{
		fits_report_error(stderr, status);
		status = 0;
		fits_close_file(fptr, &status);
		return(0);
	}


	/* Now store in the specified structure */

	for (dec=0; dec<npix; dec++)
	{
		for (ra=0; ra<npix; ra++)
		{
			transfer_function[ra][dec][0] = rez[dec][ra];
			transfer_function[ra][dec][1] = imz[dec][ra];
		}
	}

	free_double_matrix_contiguous(&rez);
	free_double_matrix_contiguous(&imz);


	/* Close the fits file */

	fits_close_file(fptr, &status);

	return(1);
}

int read_bolocam_fits_file_keywords(const char *fits_file, struct BigSZKeywords *ptr_to_keywords)
{
	fitsfile *fptr;
	int naxis, naxis1, naxis2;
	double pixsize, rmax, bmaj, bmin, unit_con, cd1, cd2, crpix1, crpix2, crval1, crval2;
	char ctype1[9], ctype2[9];
	int status = 0;

	if (ptr_to_keywords == NULL)
	{
		fprintf(stderr,"input pointer to BigSZKeywords structure is NULL\n");
		return(0);
	}

	if (fits_file == NULL)
	{
		fprintf(stderr,"input pointer to the fits_file is NULL\n");
		return(0);
	}

	/* Open the fits file */

	fits_open_file(&fptr, fits_file, READONLY, &status);

	/* Read in standard keywords */

	fits_read_key(fptr, TINT, "NAXIS", &naxis, NULL, &status);
	fits_read_key(fptr, TINT, "NAXIS1", &naxis1, NULL, &status);
	fits_read_key(fptr, TINT, "NAXIS2", &naxis2, NULL, &status);
	fits_read_key(fptr, TDOUBLE, "PIXSIZE", &pixsize, NULL, &status);
	fits_read_key(fptr, TDOUBLE, "BMAJ", &bmaj, NULL, &status);
	fits_read_key(fptr, TDOUBLE, "BMIN", &bmin, NULL, &status);
	fits_read_key(fptr, TDOUBLE, "UNIT_CON", &unit_con, NULL, &status);
	fits_read_key(fptr, TSTRING, "CTYPE1", &ctype1, NULL, &status);
	fits_read_key(fptr, TSTRING, "CTYPE2", &ctype2, NULL, &status);
	fits_read_key(fptr, TDOUBLE, "CD1_1", &cd1, NULL, &status);
	fits_read_key(fptr, TDOUBLE, "CD2_2", &cd2, NULL, &status);
	fits_read_key(fptr, TDOUBLE, "CRPIX1", &crpix1, NULL, &status);
	fits_read_key(fptr, TDOUBLE, "CRPIX2", &crpix2, NULL, &status);
	fits_read_key(fptr, TDOUBLE, "CRVAL1", &crval1, NULL, &status);
	fits_read_key(fptr, TDOUBLE, "CRVAL2", &crval2, NULL, &status);

	/* Error check */

	if (status)
	{
		fits_report_error(stderr, status);
		status = 0;
		fits_close_file(fptr, &status);
		return(0);
	}

	/* Check if RMAX keyword exists */

	fits_read_key(fptr, TDOUBLE, "RMAX", &rmax, NULL, &status);

	if (status)
	{
		if (status == 202)
		{
			rmax = 0.0;
			status = 0;
		}
		else
		{
			fits_report_error(stderr, status);
			status = 0;
			fits_close_file(fptr, &status);
			return(0);
		}
	}

	/* Save to keyword structure */

	ptr_to_keywords->scaling_factor = 1;
	ptr_to_keywords->npix = naxis1;
	ptr_to_keywords->pixsize = pixsize*ARCSEC_TO_DEG;			/* Convert from arc-seconds to degrees */
	ptr_to_keywords->rmax = rmax*ARCSEC_TO_DEG;
	ptr_to_keywords->fwhm = bmaj;
	ptr_to_keywords->sigma = bmaj/SIGMA_TO_FWHM;
	ptr_to_keywords->unit_con = unit_con;
	ptr_to_keywords->cd1 = cd1;
	ptr_to_keywords->cd2 = cd2;
	ptr_to_keywords->crpix1 = crpix1;
	ptr_to_keywords->crpix2 = crpix2;
	ptr_to_keywords->crval1 = crval1;
	ptr_to_keywords->crval2 = crval2;

	ctype1[8] = '\0';
	ctype2[8] = '\0';
	strncpy(ptr_to_keywords->ctype1, ctype1, 9);
	strncpy(ptr_to_keywords->ctype2, ctype2, 9);

	/* Close fits file */

	fits_close_file(fptr, &status);

	return(1);
}

int define_sz_model_map(struct BigSZKeywords *pkey)
{
	int scaling_factor, npix_beam, npix_model, min_npix_model, npad;

	/* Convolution of two gaussians with zero mean is also a gaussian with
	zero mean and sigma_1x2 = sqrt(sigma_1^2 + sigma_2^2).  Convolve 
	the beam PSF and the pointing uncertainty. */

	pkey->sigma = sqrt(pkey->sigma*pkey->sigma + SIGMA_POINTING*SIGMA_POINTING);
	pkey->fwhm = pkey->sigma*SIGMA_TO_FWHM;

	/* Ensure that the pixel size is not larger than MAX_RES */
	
	if (pkey->pixsize > MAX_RES)
		scaling_factor = (int) ceil(pkey->pixsize/MAX_RES);
	else
		scaling_factor = 1;

	pkey->scaling_factor = scaling_factor;
	pkey->pixsize /= scaling_factor;
	pkey->npix *= scaling_factor;
	pkey->crpix1 = (pkey->crpix1 - 0.5)*scaling_factor + 0.5;
	pkey->crpix2 = (pkey->crpix2 - 0.5)*scaling_factor + 0.5;

	/* Expand the map so that the region where data exists is not
	influenced by edge effects during convolution */

	npix_beam = (int) floor(4.0*pkey->fwhm/pkey->pixsize);

	npix_model = pkey->npix + npix_beam*2;
	min_npix_model = (int) ceil(MIN_MODEL_SIZE/pkey->pixsize);

	if (npix_model < min_npix_model)
		npad = (int) ceil((min_npix_model - pkey->npix)/2.0);
	else
		npad = (int) ceil((npix_model - pkey->npix)/2.0);

	pkey->npix += 2*npad;
	pkey->crpix1 += npad;
	pkey->crpix2 += npad;

	return(1);

}

int generate_bolocam_beam_kernel(fftw_complex **beam, struct BigSZKeywords key)
{
	double x0, y0, sigma, norm;
	int xx, yy, npix;

	npix = key.npix - !(key.npix%2);
	sigma = key.sigma/key.pixsize;

	x0 = (npix-1)/2;
	y0 = (npix-1)/2;

	norm = 2.0*DPI*sigma*sigma;

	for (xx=0; xx<npix; xx++)
	{
		for (yy=0; yy<npix; yy++)
		{
			beam[xx][yy][0] = exp(-((xx-x0)*(xx-x0) + (yy-y0)*(yy-y0))/(2.0*sigma*sigma))/norm;
			beam[xx][yy][1] = 0.0;	
		}
	}

	/* Zero-pad */

	zero_pad_sz(beam, npix, 2*key.npix, npix, 2*key.npix);

	return(1);
}

int count_pixels_within_aperture(int nx, int ny, double x0, double y0, double rmax)
{
	int xx, yy;
	int nin = 0;
	double xm, ym, rmax2;

	if ((x0 < 0) || (x0 > (nx-1)) || (y0 < 0) || (y0 > (ny-1)))
		return(0);

	rmax2 = rmax*rmax;

	xm = x0 > (nx - 1 - x0) ? x0 : (nx - 1 - x0);
	ym = y0 > (ny - 1 - y0) ? y0 : (ny - 1 - y0);

	if (rmax2 < (xm*xm + ym*ym))
	{
		for (xx=0; xx < nx; xx++)
		{
			for (yy=0; yy < ny; yy++)
			{
				nin += (((xx-x0)*(xx-x0) + (yy-y0)*(yy-y0)) <= rmax2);
			}
		}
	}
	else
	{
		nin = nx*ny;
	}

	return(nin);
}

int zero_pad_sz(fftw_complex **map, int nx1, int nx2, int ny1, int ny2)
{
	int xx, yy;

	for (xx=0; xx<nx1; xx++)
	{
		for (yy=ny1; yy<ny2; yy++)
		{
			map[xx][yy][0] = 0.0;
			map[xx][yy][1] = 0.0;
		}
	}

	for (xx=nx1; xx<nx2; xx++)
	{
		for (yy=0; yy<ny2; yy++)
		{
			map[xx][yy][0] = 0.0;
			map[xx][yy][1] = 0.0;	
		}
	}

	return(1);
}


/* functions for allocating multi-dimensional arrays */
int allocate_double_matrix_contiguous(double ***matrix, int nrows, int ncols)
{
	int rr;

	*matrix = (double **) malloc(nrows*sizeof(double *));
	if (*matrix == NULL)
	{
		fprintf(stderr,"out of memory during allocation of first matrix dimension\n");
		return(0);
	}

	(*matrix)[0] = (double *) malloc(nrows*ncols*sizeof(double));
	if ((*matrix)[0] == NULL)
	{
		fprintf(stderr,"out of memory during allocation of second matrix dimension\n");
		free(*matrix);
		return(0);
	}

	for (rr=1; rr < nrows; rr++)
		(*matrix)[rr] = (*matrix)[0] + rr*ncols;

	return(1);
}

int free_double_matrix_contiguous(double ***matrix)
{
	free((*matrix)[0]);
	free((*matrix));

	return(1);
}

int allocate_fftw_complex_matrix_contiguous(fftw_complex ***matrix, int nrows, int ncols)
{
	int rr;

	*matrix = (fftw_complex **) malloc(nrows*sizeof(fftw_complex *));
	if (*matrix == NULL)
	{
		fprintf(stderr,"out of memory during allocation of first matrix dimension\n");
		return(0);
	}

	(*matrix)[0] = (fftw_complex *) fftw_malloc(nrows*ncols*sizeof(fftw_complex));
	if ((*matrix)[0] == NULL)
	{
		fprintf(stderr,"out of memory during allocation of second matrix dimension\n");
		free(*matrix);
		return(0);
	}

	for (rr=1; rr < nrows; rr++)
		(*matrix)[rr] = (*matrix)[0] + rr*ncols;

	return(1);
}

int free_fftw_complex_matrix_contiguous(fftw_complex ***matrix)
{
	fftw_free((*matrix)[0]);
	free((*matrix));

	return(1);
}
