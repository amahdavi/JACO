#include <jaco.h>
#include <fitsio.h>
#include <string.h>
#include <math.h>
#include <stdio.h>
#include <standalone.h>
#include <gsl/gsl_randist.h>
#include <gsl/gsl_rng.h>
#include <gsl/gsl_math.h>
#include <gsl/gsl_eigen.h>
#include <gsl/gsl_linalg.h>
#include <gsl/gsl_cblas.h>
#include <sys/time.h>

unsigned long int get_random_seed()
{
	unsigned int seed;
	struct timeval tv;
	FILE *devrandom;

	if ((devrandom = fopen("/dev/random","r")) == NULL) {
	   gettimeofday(&tv,0);
	   seed = tv.tv_sec + tv.tv_usec;
	} else {
	   fread(&seed,sizeof(seed),1,devrandom);
	   fclose(devrandom);
	}

	return(seed);
}


int init_standalone(struct fitdata *data)
{
  struct pha *phalist=NULL;
  struct rmf *rmflist=NULL;
  struct arf *arflist=NULL;
  struct ebounds *elimits=NULL;
  //double *wlr1,*wlr2,ratio,x;
  long i,j,k,l;
  long ndata, wlcount=0,szcount=0,velcount=0;
  struct jaco_state *js = data->js;
  char simulation_prefix[] = "sim";

  js->data = data;

  if (jaco_read_data(js)) return -1;

  JS(xraycount); JS(fitsname); 
  JS(calcwl); JS(calcsz); JS(calcxray); JS(calcvel);
  JS(efit1); JS(efit2); JS(r1full);

  // This snippet depends on the instrument being only Chandra&XMM
  // Determine minimum radius to fit

  // Only implement Rcut for XMM data when Chandra data is also present
  
  data->bigdata = NULL;
  data->bigerr = NULL;

  ndata = 0;
  // CFITSIO not thread safe! The following OpenMP directive won't work
  //#pragma omp parallel for private(i) reduction(+:ndata)

  k = 0;
  int nback;

  if (calcxray || calcsz) {


    nback = parse_list(js->backlist,',',&js->backsyserr);
    if (nback != NINSTR-3) 
      BYE("You must provide %d background systematic values, comma-separated. You said %s.",NINSTR-3,js->backlist);

    js->havecxo = js->havexmm = 0;
    int id;
    int *checkback = sci_ivector(NINSTR);

    for (i = 0; i < xraycount; ++i) {
      id = js->instrument[i];
      if (id > PN)
	js->havecxo=1;
      if (id <= PN)
	js->havexmm=1;

    }  
    /*       if (!checkback[id]) { */
/*  	checkback[id] = 1;  */
/* 	if (j == nback) { */
/* 	  jaco_printf("No background file found for instrument %s",idcode[id]); */
/* 	  js->residual_background[id].nbins = 0; */
/* 	} else { */
/* 	  readpha(backfiles[j],&js->residual_background[id]); */
/* 	  readarf(js->residual_background[id].arfname,&js->backarf[id]); */
/* 	  jaco_printf("Background: %d from %s (%f,%d).\n", */
/* 		 js->residual_background[id].nchannels, */
/* 		 backfiles[j], */
/* 		 js->residual_background[id].arcminarea, */
/* 		 js->backarf[id].nebins); */
/* 	  for (k = 0; k < js->residual_background[id].nchannels; ++k) */
/* 	    js->residual_background[id].counts[k] /= */
/* 	      js->residual_background[id].arcminarea* */
/* 	      js->residual_background[id].exposure; */
/* 	} */
/*       } */
/*     } */
    
    if (js->havecxo)
      for (i = 0; i < NINSTR; ++i) {
	if (js->cutradius > 0. && i >= MOS1 && i <= PN) 
	  data->Rcut[i] = js->cutradius;
	else 
	  data->Rcut[i] = 0.;
      }

    data->arflist = (struct arf *)malloc(xraycount*sizeof(struct arf));
    data->rmflist = (struct rmf *)malloc(xraycount*sizeof(struct rmf));
    data->phalist = (struct pha *)malloc(xraycount*sizeof(struct pha));
    data->elimits = (struct ebounds *)malloc(xraycount*sizeof(struct ebounds));
    data->wldata = NULL;
    js->fitcomp = (char **)malloc(xraycount*sizeof(char *));
    
    phalist = data->phalist;
    arflist = data->arflist;
    rmflist = data->rmflist;
    elimits = data->elimits;
  }

  k = 0;
  if (calcxray) {
    for (i = 0; i < xraycount; ++i) {
      phalist[i].filename = js->fitsname[i];
	  phalist[i].is_sim = !strncmp(phalist[i].filename, simulation_prefix, strlen(simulation_prefix));
      readpha(fitsname[i],&phalist[i],0);
      //js->fitcomp[i] = sci_strcat(fitsname[i],".dat");

      if (strlen(phalist[i].arfname))
	readarf(phalist[i].arfname,&arflist[i]);

      if (strlen(phalist[i].rmfname))
	readrmf(phalist[i].rmfname,&rmflist[i],&elimits[i]);

      // Handle the case of missing ARF (e.g. Con-X RSP files)
      if (!strlen(phalist[i].arfname)) {
	printf("No ARF found in %s.\n",fitsname[i]);
	arflist[i].nebins = rmflist[i].nebins;
	arflist[i].elo = rmflist[i].elo;
	arflist[i].ehi = rmflist[i].ehi;
	arflist[i].effarea = sci_dvector(arflist[i].nebins);
	printf("OK\n");
	for (j = 0; j < arflist[i].nebins; ++j)
	  arflist[i].effarea[j] = 1.;
      }

      bin_pha(&phalist[i],&elimits[i],efit1,efit2,data->js->syserr);

      //if ((gsl_fcmp(phalist[i].arcminarea,0.01*js->ringarea[i],0.01)<0) &&
      //  gsl_fcmp(phalist[i].arcminarea,js->goodarea[i],0.01)) {
      //printf("Spec %d: BACKSCAL=%e datafile=%e.\n",i,phalist[i].arcminarea,
      //       js->goodarea[i]);
      //printf("------>Using datafile value instead of BACKSCAL\n");
      phalist[i].arcminarea = js->goodarea[i];
	//}

      if (data->Rcut[js->instrument[i]]-r1full[i] < 1.E-6) 
	ndata += phalist[i].nbins;
      
      //printf("%d %d %s %d\n",i,xraycount,phalist[i].rmfname,phalist[i].nchannels);
      //printf("%s %s (%d) %s (%d %d)\n",
      //fitsname[i],phalist[i].arfname,arflist[i].nebins,
      //phalist[i].rmfname,rmflist[i].nebins,rmflist[i].nchannels[0]);

    }

    if (!ndata)
      BYE("No data to model.\n");

    data->bigdata = sci_dvector(ndata);
    data->bigerr = sci_dvector(ndata);

    // Check to see if the files are all on the same energy/channels grid
    for (i = 0; i < xraycount; ++i) {
      if ((arflist[i].nebins != arflist[0].nebins) ||
	  (rmflist[i].nebins != rmflist[0].nebins) ||
	  (rmflist[i].nebins != arflist[i].nebins))
	BYE("All RMF and ARF files must have identical energy bins.%s","");
      
      if (arflist[i].elo[0] > efit1)
	BYE("your lower energy limit > calibration grid lower limit (%f).",
	    arflist[i].elo[0]);
      
      if (arflist[i].ehi[arflist[i].nebins-1] < efit2)
	BYE("your higher energy limit < calibration grid lower limit (%f).",
	    arflist[i].ehi[arflist[i].nebins-1]);
    
      if (data->Rcut[js->instrument[i]]-r1full[i] < 1.E-6) 
	for (j = 0; j < phalist[i].nbins; ++j) {
	  data->bigdata[k] = phalist[i].flux[j];
	  if (!phalist[i].is_sim)
			phalist[i].sig[j] = sqrt(pow(phalist[i].sig[j],2.) +
				   atof(js->backsyserr[js->instrument[i]])*phalist[i].arcminarea/phalist[i].exposure);
	    
	  data->bigerr[k] = phalist[i].sig[j];
	  //if (js->cutradius > r1full[i])
	  //data->bigerr[k] += 0.3*data->bigdata[k];
	  ++k;
	}
    }

    #pragma omp parallel for private(i,j,k,l) schedule(dynamic)
    for (i = 0; i < xraycount; ++i)
      for (j = 0; j < rmflist[i].nebins; ++j) {
	for (k = 0; k < rmflist[i].ngrp[j]; ++k)
	  for (l = 0; l < rmflist[i].nchannels[j][k]; ++l) 
	    rmflist[i].response[j][k][l] *= arflist[i].effarea[j];
      }
    
  } else if (calcsz) {
    
    for (i = 0; i < xraycount; ++i)  {
      arflist[i].nebins = 2;
      arflist[i].elo = sci_dvector(2);
      arflist[i].ehi = sci_dvector(2);
      arflist[i].elo[0] = 1.;
      arflist[i].elo[1] = 1.1;
      arflist[i].ehi[0] = 1.1;
      arflist[i].ehi[1] = 1.2;

    }
  }

  if (calcwl) {
    k = ndata;
    data->wldata = (struct weaklensing *)malloc(sizeof(struct weaklensing));
    wlcount = read_ascii(js->wldata,"%f %f %f",NULL,
			 &data->wldata->wlrad,
			 &data->wldata->shear,
			 &data->wldata->shearerr);
    //data->wldata->wlrad = sci_dvector(wlcount);
    data->wldata->modelshear = sci_dvector(wlcount);
	data->wldata->nwl = wlcount;
    ndata += wlcount;									// Siegel: brought this up from below
    sci_dvector_resize(&data->bigdata,ndata);
    sci_dvector_resize(&data->bigerr,ndata);

    /* Siegel Begin */
    /* If a covariance matrix was input for the lensing data,
       then we determine the eigenvectors and eigenvalues of 
       the inverse covariance matrix.  We use the eigenvectors
       to rotate wlshear to a basis where the inverse of the 
       covariance matrix is diagonal.  We set wlshearerr 
       equal to 1/sqrt(eigenvalues). */
	
	if (js->wlcov)
	{		
		if (read_lensing_covariance_and_rotate_data(js->wlcov, data->wldata))
			BYE("Error reading lensing covariance file %s.\n", js->wlcov);
	}
	else
	{
		data->wldata->rot_modelshear = data->wldata->modelshear;
		data->wldata->rot_shear = data->wldata->shear;
		data->wldata->rot_shearerr = data->wldata->shearerr;
	}	
	/* Siegel End */


    for (i = 0; i < (int)wlcount; ++i) {
      //ratio = wlr1[i]/wlr2[i];
      //x = ratio*ratio;
      //data->wldata->wlrad[i] = (2./3)*(wlr2[i]-x*wlr1[i])/(1.-x);
      
      data->wldata->wlrad[i] /= 60.;
      data->wldata->rot_shearerr[i] *= js->wllss;

      data->bigdata[k] = data->wldata->rot_shear[i];
      data->bigerr[k++] = data->wldata->rot_shearerr[i];
    }
    //free(wlr1); free(wlr2);
  } else
    data->wldata = NULL;

  if (calcvel) {
    k = ndata;
    data->vels = (struct dispersion *)malloc(sizeof(struct dispersion));
    velcount = read_ascii(js->veldata,"%f %f %f",NULL,
			 &data->vels->rad,
			 &data->vels->disp,
			 &data->vels->disperr);
    data->vels->modeldisp = sci_dvector(velcount);
    ndata += velcount;
    sci_dvector_resize(&data->bigdata,ndata);
    sci_dvector_resize(&data->bigerr,ndata);

    for (i = 0; i < (int)velcount; ++i) {
      data->bigdata[k] = data->vels->disp[i];
      data->bigerr[k++] = data->vels->disperr[i];
    }
    data->vels->nvel = velcount;
    //free(wlr1); free(wlr2);
  } else
    data->vels = NULL;

#if defined(WITH_SZ) || defined(WITH_BOLOCAM)
  if (calcsz) {
    k = ndata;
    data->szdata = (struct sz *)malloc(sizeof(struct sz));
    szcount = js->szparams.nest;

    //These are now transferred internally within big_initialize_sz
    /* szcount = read_ascii(js->szdata,"%f %f %f",NULL, */
    /* 			 &data->szdata->szx,&data->szdata->szy, */
    /* 			 &data->szdata->sze); */

    /* if (!szcount) */
    /*   BYE("Zero-length SZ data file %s.\n",js->szdata); */

    /* for (i = 0; i < szcount; ++i)  */
    /*   if (data->szdata->szx[i] >= js->sz0) ++ndata; */
    
    data->szdata->modely = sci_dvector(szcount); 
 
    data->szdata->szy = js->szparams.data;
    data->szdata->sze = js->szparams.errors;
    data->szdata->nsz = szcount;

    ndata += szcount;
    sci_dvector_resize(&data->bigdata,ndata);
    sci_dvector_resize(&data->bigerr,ndata);

    /* Siegel Begin */
    /* If a covariance matrix was input for the sz data,
       then we determine the eigenvectors and eigenvalues of 
       the inverse covariance matrix.  We use the eigenvectors
       to rotate szy to a basis where the inverse of the 
       covariance matrix is diagonal.  We set sze 
       equal to 1/sqrt(eigenvalues). */
	
	if (js->szcov)
	{		
		if (read_sz_covariance_and_rotate_data(js->szcov, data->szdata))
			BYE("Error reading lensing covariance file %s.\n", js->szcov);
	}
	else
	{
		data->szdata->rot_modely = data->szdata->modely;
		data->szdata->rot_szy = data->szdata->szy;
		data->szdata->rot_sze = data->szdata->sze;
	}	
	/* Siegel End */


    for (i = 0; i < szcount; ++i) 
      /* if (data->szdata->szx[i] >= js->sz0) { */
    {
		data->bigdata[k] = data->szdata->rot_szy[i];
		data->bigerr[k++] = data->szdata->rot_sze[i];
    }
      /* } */

  } else
#endif
    data->szdata = NULL;

  if (ndata) {
    data->ndata = ndata;
    data->bigmodel = sci_dvector(ndata);
  } else
    BYE("No suitable data found!\n");

  return 0.;
}

#define CHANGEDPAR(x) (DIFFER(params[(x)].value,oldparams[(x)].value))

int get_bigmodel(double *x,double *params,double *model, double *errors,
		 double *prior, unsigned long ndata, void *datapointer)
{
  int i,j,k,l;
  FILE *fitfile=NULL, *rotfitfile=NULL, *simfile=NULL;
  char proffilename[100];
  long npoints=0;
  struct pha *s, *phalist;
  struct arf *arflist;
  struct rmf *rmflist;
  struct ebounds *elimits;
  struct fitdata *data = (struct fitdata *)datapointer;
  struct weaklensing *wldata = data->wldata;
  struct dispersion *vels = data->vels;
  struct sz *szdata = data->szdata;
  struct jaco_state *js = data->js;
  double *spectrum,totchisq=0.,diff;
  char *fname;
  char simulation_prefix[100] = "sim.";
  char *temp_string1, *temp_string2, *ipb;            // Siegel

  /* Siegel Begin */
  /* If we are creating a simulation, then we need a random seed
     and gsl_generator that is different from js->rng.
     This is because the js->rng seed is always set to 0
 	 when running in evaluation mode, which would result in
	 the same random fluctuations being added to every
	 simulation. */
  //double meanflux=0.0;
  //struct pha *background;
  gsl_rng *simulation_generator;
  unsigned long simulation_seed;

  if (js->dosim > 0)
  {
	simulation_seed = get_random_seed();
	simulation_generator = gsl_rng_alloc(gsl_rng_mt19937);
	gsl_rng_set(simulation_generator,simulation_seed);
	
	/* Check to see whether the simulation prefix has been 
     specified in the simfile nonfittable string parameter */
  	if ((strncmp(js->simfile, simulation_prefix, 3) == 0) && (js->simfile[strlen(js->simfile)-1] == '.'))
			strcpy(simulation_prefix, js->simfile);
  }
  /* Siegel End */

  phalist = data->phalist;
  arflist = data->arflist;
  rmflist = data->rmflist;
  elimits = data->elimits;
  spectrum = data->spectrum;

  // No Bayesian prior
  (*prior) = 0.;
  
  JS(calcxray); JS(calcsz); JS(xraycount); JS(instrument); 
  JS(r1full); JS(annuli);

  params[ANNULUS] = 0;
  if (calcxray || calcsz) 
    for (i = 0; i < xraycount; ++i) {
      params[INSTRUMENT] = (double)instrument[i];
      params[ANNULUS] = (double)annuli[i];
      
      
      jaco_calc_models(params,&arflist[i].nebins,spectrum,arflist[i].elo,
		       data->js);

      if (calcxray && data->Rcut[instrument[i]]-r1full[i] < 1.E-6) {
	if (fabs(spectrum[0]) > 1.E6 && (js->debug > 0 || !js->specwarn)) {
	  js->specwarn = 1;
	  printf("ERROR: SPECTRUM: %E %d\n",spectrum[0],annuli[i]);
	}
	s = &phalist[i];
	
	for (j = 0; j < s->nchannels; ++j) 
	  s->model[j] = 0;

	for (j = 0; j < rmflist[i].nebins; ++j) {

	  // This used to be here but got moved up the data reading
	  // section for better performance
	  // spectrum[j] *= arflist[i].effarea[j];

	  for (k = 0; k < rmflist[i].ngrp[j]; ++k) 
	    for (l = 0; l < rmflist[i].nchannels[j][k]; ++l)
	      s->model[rmflist[i].chanstart[j][k]+l-phalist[i].channel[0]] += 
		rmflist[i].response[j][k][l]*spectrum[j];
	}
      
	if (js->fitwrite == 1) // && r1full[i] > js->cutradius)
	{
	  /* Siegel Begin */	
	  temp_string1 = sci_strdup(js->fitsname[i]);
	  ipb = strstr(temp_string1, ".fits");
	  if (ipb != NULL) *ipb = '\0';
	  temp_string2 = sci_strcat(temp_string1, ".xrfit.dat");
	  free(temp_string1);
	  temp_string1 = sci_strcat(js->profname, ".");
	  js->fitcomp[i] = sci_strcat(temp_string1, temp_string2);
	  free(temp_string1);  free(temp_string2);
	  /* Siegel End */
	  fitfile = fopen(js->fitcomp[i],"w+");
	}
	else
	{
	  fitfile = NULL;
	}

	for (j = 0; j < s->nbins; ++j) {
	  int id = js->instrument[i];

	  s->mcounts[j] = 0;
/* 	  // Only correct Chandra data if no XMM */
 	    s->mcounts[j] = params[PLNORM0+id]*(pow(elimits[i].emax[s->chanhi[j]],1.-params[PLSLOPE0+id])-pow(elimits[i].emin[s->chanlo[j]],1.-params[PLSLOPE0+id]))/(1.-params[PLSLOPE0+id]);  

	  for (k = s->chanlo[j]; k <= s->chanhi[j]; ++k) {
	    s->mcounts[j] += s->model[k];
	  }

	  s->mcounts[j] *= s->arcminarea;
	  
	  if (js->instrument[i] == ACISI || js->instrument[i] == ACISS)
	    //if (elimits[i].emax[s->chanhi[j]] < 2.0 && 
	    //elimits[i].emax[s->chanhi[j]] > 1.0)
	    s->mcounts[j] *= pow(elimits[i].emax[s->chanhi[j]],js->tbias);

	  if (js->fitwrite > 1) {
	    // These (statistic calculation) are now moved to hrothgar.
	    // Only do them to evaluate a profile's likelihood
	    //diff = (s->mcounts[j]-s->flux[j])/s->sig[j];
	    //totchisq += diff*diff;
	  }
	  
	  if (!data->fitannuli) 
	    model[npoints+j] = s->mcounts[j];
	  
	  if (fitfile) {
	    diff =  elimits[i].emax[s->chanhi[j]]-
	      elimits[i].emin[s->chanlo[j]];
	    fprintf(fitfile,"%15.8E %15.8E %15.8E %15.8E %15.8E\n",
		    (elimits[i].emin[s->chanlo[j]]+
		     elimits[i].emax[s->chanhi[j]])/2.,
		    s->flux[j]/diff,s->sig[j]/diff,s->mcounts[j]/diff, diff);
		/*printf("%E %E %E %E\n",
		    (elimits[i].emin[s->chanlo[j]]+
		     elimits[i].emax[s->chanhi[j]])/2.,
		    s->flux[j]/diff,s->sig[j]/diff,s->mcounts[j]/diff);*/
	  }

	  if (js->dosim > 0) 
	  {
		/* Gaussian Noise */
	    s->mcounts[j] += gsl_ran_gaussian(simulation_generator,s->sig[j]);
		
		/* Poisson Noise
		meanflux = s->exposure*s->mcounts[j];
		meanflux = gsl_ran_poisson(simulation_generator, meanflux);
		s->mcounts[j] = meanflux / s->exposure;
		
		if (s->backname)
		{
			background = (struct pha *) s->background;
			for (k = s->chanlo[j]; k <= s->chanhi[j]; ++k)
			{
				background->counts[k] = gsl_ran_poisson(simulation_generator, background->counts[k]);
			}
		}*/
	  }
	
	}
	
	/* Siegel Begin -- Took this outside of the loop over nbins */
	if (js->dosim > 0)
	{
		fname = sci_strcat(simulation_prefix,s->filename);
		printf("Outputting simulation:  %s\n", fname);
	    writepha(fname,s);
	    free(fname);
	}
	/* Siegel End */
	
	
        //s->totchisq = totchisq;
	//finalchisq += s->totchisq;
        npoints += s->nbins;
	if (fitfile) fclose(fitfile);
      
	// Fisher's transofrmation to make X^2 approximately normal
	if (data->fitannuli) {
	  model[npoints] = -sqrt(2.*s->totchisq);
	  ++npoints;
	}
      }
    }

  if (wldata) {

    if (js->fitwrite == 1) {
      fname = sci_strcat(js->profname,".wlfit.dat");
      fitfile = fopen(fname,"w+");
      free(fname);

      fname = sci_strcat(js->profname,".rotwlfit.dat");
      rotfitfile = fopen(fname,"w+");
      free(fname);
    }

    params[INSTRUMENT] = WL;
    params[ANNULUS] = calcxray+(vels!=NULL);
    jaco_calc_models(params,&wldata->nwl,wldata->modelshear,wldata->wlrad,js);

    /* Siegel Begin */
	if (js->wlcov)
	{
		rotate_lensing_model(wldata);
	}
    /* Siegel End */

    for (i = 0; i < wldata->nwl; ++i) {
				
      // Moved to hrothgar:
      diff = (wldata->rot_modelshear[i]-wldata->rot_shear[i])/wldata->rot_shearerr[i];
      //if (chiarr) chiarr[npoints] = diff;
      totchisq += diff*diff;

      model[npoints] = wldata->rot_modelshear[i];

      if (fitfile)
	  {
		fprintf(fitfile,"%15.8E %15.8E %15.8E %15.8E\n",wldata->wlrad[i],wldata->shear[i],
										wldata->shearerr[i],wldata->modelshear[i]);
	  }
	
	  if (rotfitfile)
	  {
		fprintf(rotfitfile,"%15.8E %15.8E %15.8E\n",wldata->rot_shear[i], 
										wldata->rot_shearerr[i],model[npoints]);
	  }

      ++npoints;
    }
    if (fitfile)
      fclose(fitfile);
	if (rotfitfile)
		fclose(rotfitfile);


	/* Siegel Begin */
	/* If dosim == 1, then output model simulation similair to what is done for xray.
	   Should move this to separate function if it gets more complicated (i.e. covariance) */
	if (js->dosim == 1)
	{
		fname = sci_strcat(simulation_prefix,js->wldata);
		printf("Outputting simulation:  %s\n", fname);
		simfile = fopen(fname,"w+");
		free(fname);
		
		for (i = 0; i < wldata->nwl; ++i)
		{
			fprintf(simfile, "%f\t%f\t%f\n", wldata->wlrad[i]*60., 
						wldata->modelshear[i] +  gsl_ran_gaussian(simulation_generator,wldata->shearerr[i]), 
						wldata->shearerr[i]);
		}
		
		fclose(simfile);
	}
	/* Siegel End */

  }
  //printf("WL %f %ld\n",finalchisq,npoints);

  if (vels) {

    if (js->fitwrite == 1) {
      fname = sci_strcat(js->profname,".vels.dat");
      fitfile = fopen(fname,"w+");
      free(fname);
    }

    params[INSTRUMENT] = VELS;
    params[ANNULUS] = (wldata != NULL)+calcxray;
    jaco_calc_models(params,&vels->nvel,vels->modeldisp,vels->rad,js);

    for (i = 0; i < vels->nvel; ++i) {

      // Moved to hrothgar:
      diff = (vels->modeldisp[i]-vels->disp[i])/vels->disperr[i];
      //if (chiarr) chiarr[npoints] = diff;
      totchisq += diff*diff;

      model[npoints] = vels->modeldisp[i];
      if (fitfile)
	fprintf(fitfile,"%15.8E %15.8E %15.8E %15.8E\n",vels->rad[i],vels->disp[i],
		vels->disperr[i],model[npoints]);

      ++npoints;
    }
    if (fitfile) 
      fclose(fitfile);

  }

  if (szdata) {
    params[INSTRUMENT] = SZ;
    jaco_calc_models(params,&szdata->nsz,szdata->modely,NULL,js);

    if (js->fitwrite == 1) {
      fname = sci_strcat(js->profname,".szfit.dat");
      fitfile = fopen(fname,"w+");
      free(fname);

	  /* Siegel Begin */
      fname = sci_strcat(js->profname,".rotszfit.dat");
      rotfitfile = fopen(fname,"w+");
      free(fname);
	  /* Siegel End */
    }

   /* Siegel Begin */
	if (js->szcov)
	{
		rotate_sz_model(szdata);		
	}
   /* Siegel End */
	
    for (i = 0; i < szdata->nsz; ++i) 
      /* if (szdata->szx[i] >= data->js->sz0)  */
	{
		// Moved to hrothgar:
		diff = (szdata->rot_modely[i]-szdata->rot_szy[i])/szdata->rot_sze[i];
		//if (chiarr) chiarr[npoints] = diff;
		totchisq += diff*diff;
	
		model[npoints] = szdata->rot_modely[i];

		if (fitfile)
		  fprintf(fitfile,"%d %15.8E %15.8E %15.8E\n",i,szdata->szy[i],
			  szdata->sze[i],szdata->modely[i]);
		
		/* Siegel Begin */
		if (rotfitfile)
		  fprintf(rotfitfile,"%d %15.8E %15.8E %15.8E\n",i,szdata->rot_szy[i],
			  szdata->rot_sze[i],szdata->rot_modely[i]);
		/* Siegel End */		

		++npoints;
      }

    if (fitfile)
      fclose(fitfile);
	if (rotfitfile)
		fclose(rotfitfile);		// Siegel
  }

  
  /* Siegel Begin */
  if (js->dosim > 0)
  {
	gsl_rng_free(simulation_generator);
  }
  /* Siegel End */
    
/*   int infochanged; */
/*   long nchanged=0; */
/*   static int initialprof=1; */
/*   static double minchisq=1.E30; */
  
/*   if (gsl_fcmp(totchisq,minchisq,1.E-4) < 0) { */
/*     if (!initialprof) */
/*       printf("Found better chisq %E < %E!\n",totchisq,minchisq); */
/*     initialprof=1; */
/*     minchisq=totchisq; */
/*   } */

/*   if (js->fitwrite > 1) { */
/*     infochanged = 0; */
/*     if (totchisq < minchisq+js->nfitparams) { */
/*       for (i = 0; i < js->nannuli+js->ndelta; ++i) { */
/* 	for (j = 0; j < INFO_SIZE; ++j) { */
/* 	  if (js->infoarray[i][j] < js->outprof[0][i][j] || initialprof) { */
/* 	    infochanged = 1; */
/* 	    js->outprof[0][i][j] = js->infoarray[i][j]; */
/* 	  } */
/* 	  if (js->infoarray[i][j] > js->outprof[1][i][j] || initialprof) { */
/* 	    infochanged = 1;  */
/* 	    js->outprof[1][i][j] = js->infoarray[i][j]; */
/* 	  } */
/* 	} */
/*       } */
/*     } */

/*     initialprof=0; */

/*     if (infochanged) { */
/*       nchanged = 0; */
/*       for (k = 0; k < 2; ++k) { */
/* 	snprintf(proffilename,99, */
/* 		 "%s-profiles-%03d.%s",js->profname,js->mpirank, */
/* 		 (k ? "max" : "min")); */

/* 	fitfile = fopen(proffilename,"w+"); */
	
/* 	if (js->dosim > 1) */
/* 	  fprintf(fitfile,"# Reff  Tproj rhod rhog rhostar Md Mg Mstar"	\ */
/* 		  " fgas T3D Metal\n"); */
/* 	else */
/* 	  fprintf(fitfile,"# Reff Tproj Md Mg Mstar fgas T3D n_e "	\ */
/* 		  "pressure entropy coolingtime Tshock LX metal\n"); */
	
/* 	for (i = 0; i < js->nannuli; ++i) { */
/* 	  for (j = 0; j < INFO_SIZE; ++j) */
/* 	    fprintf(fitfile,"%E ",js->outprof[k][i][j]); */
/* 	  fprintf(fitfile,"\n"); */
/* 	} */

/* 	fclose(fitfile); */
	
/* 	snprintf(proffilename,99, */
/* 		 "%s-contrasts-%03d.%s",js->profname,js->mpirank, */
/* 		 (k ? "max" : "min")); */

/* 	fitfile = fopen(proffilename,"w+"); */
/* 	fprintf(fitfile,"# Delta rDelta cdelta MgDelta MdDelta MsDelta MtDelta\n"); */
/* 	for (i = js->nannuli; i < js->nannuli+js->ndelta; ++i) { */
/* 	  for (j = 0; j < 8; ++j) */
/* 	    fprintf(fitfile,"%E ",js->outprof[k][i][j]); */
/* 	  fprintf(fitfile,"\n"); */
/* 	} */
	
/* 	fclose(fitfile); */
/*       } */
/*     } */
/*     else { */
/*       ++nchanged; */
/*       if ((nchanged % 100) == 0) */
/* 	printf("No changes for past %ld iterations\n",nchanged); */
/*     } } */
	
      
      


  //printf("SZ %f %ld\n",finalchisq,npoints);

  // if (arrsize) *arrsize = npoints;
  //printf("WL %f %ld\n",finalchisq,npoints);
  //printf("%10E %ld\n",finalchisq,npoints);
  return 0;
}


int read_lensing_covariance_and_rotate_data(const char *fits_file, struct weaklensing *wldata)
{
	fitsfile *fptr;
	double **img;
	int status=0, anynul=0, check=0, iext, signum, naxis;
	long naxis1, naxis2, nbin, nfill;
	long rr, cc;
	long fpixel[2] = {1, 1};
	gsl_vector *eigenval;
	gsl_matrix *cov, *inv_cov, *eigenvec, *inv_eigenvec;
	gsl_permutation *permutation;
	gsl_eigen_symmv_workspace *eigen_workspace;

	if (fits_file == NULL)
	{
		fprintf(stderr,"input pointer to the fits_file is NULL\n");
		return(1);
	}

	/* Open the fits file */

	fits_open_file(&fptr, fits_file, READONLY, &status);

	fits_get_hdu_num(fptr, &iext);
	if (iext != 1)
	{
		fits_movabs_hdu(fptr, 1, NULL, &status);
	}

	/* Read in relevant keywords */

	fits_read_key(fptr, TINT, "NAXIS", &naxis, NULL, &status);
	fits_read_key(fptr, TLONG, "NAXIS1", &naxis1, NULL, &status);
	fits_read_key(fptr, TLONG, "NAXIS2", &naxis2, NULL, &status);	
	
	/* Error check */
	
	if (status)
	{
		fits_report_error(stderr, status);
		fits_close_file(fptr, &close);
		return(1);
	}
	if (naxis != 2)
	{
		fprintf(stderr,"covariance must be two-dimensional\n");
		fits_close_file(fptr, &close);
		return(1);
	}
	if (naxis1 != naxis2)
	{
		fprintf(stderr,"covariance must be a square matrix\n");
		fits_close_file(fptr, &close);
		return(1);
	}
	if (naxis1 != wldata->nwl)
	{
		fprintf(stderr,"covariance is %i x %i but data is %i elements long.\n", 
						naxis1, naxis2, wldata->nwl);
		fits_close_file(fptr, &close);
		return(1);
	}
	
	nbin = naxis1;
	nfill = nbin*nbin;
	
	/* Read in data */
		
	allocate_double_matrix_contiguous(&img, nbin, nbin);

	fits_read_pix(fptr, TDOUBLE, fpixel, nfill, NULL, &(img[0][0]), &anynul, &status);
	if (status)
	{
		fits_report_error(stderr, status);
		status = 0;
		fits_close_file(fptr, &status);
		return(1);
	}
	
	/* Check that the diagonal elements are equal to the input errors squared */
	
	for (rr=0; rr<nbin; rr++)
		check += abs(gsl_fcmp(sqrt(img[rr][rr]), wldata->shearerr[rr], 0.0001));
	
	if (check != 0)
	{
		printf("diag(cov) not equal to input errors squared, check = %i.\n", check);
		return(1);
	}
	
	/* Place data into gsl_matrix */
	
	cov = gsl_matrix_alloc(nbin, nbin);
	
	for (rr=0; rr<nbin; rr++)
	{	
		for (cc=0; cc<nbin; cc++)
		{
			gsl_matrix_set(cov, rr, cc, img[rr][cc]);
		}
	}
		
	/* Invert the covariance matrix */
	
	inv_cov = gsl_matrix_alloc(nbin, nbin);
	
	permutation = gsl_permutation_alloc(nbin);
	gsl_linalg_LU_decomp(cov, permutation, &signum);
	gsl_linalg_LU_invert(cov, permutation, inv_cov);
	gsl_permutation_free(permutation);
	
	/* Calculate eigenvalues and eigenvectors of covariance matrix */
	
	eigenval = gsl_vector_alloc(nbin);	
	eigenvec = gsl_matrix_alloc(nbin, nbin);
	
	eigen_workspace = gsl_eigen_symmv_alloc(nbin);
	gsl_eigen_symmv(inv_cov, eigenval, eigenvec, eigen_workspace);
	gsl_eigen_symmv_free(eigen_workspace);
	
	/* Place the eigenvalues into wldata structure */
	
	wldata->eigenval = sci_dvector(nbin);
	
	for (rr=0; rr<nbin; rr++)
	{
		wldata->eigenval[rr] = gsl_vector_get(eigenval,rr);
	}
	
	/* Place the eigenvector matrix in the wldata structure.
 	   Transpose in the process so that eigenvectors run along rows. */

	allocate_double_matrix_contiguous(&(wldata->eigenvec), nbin, nbin);
	
	for (rr=0; rr<nbin; rr++)
	{
		for (cc=0; cc<nbin; cc++)
		{
			wldata->eigenvec[rr][cc] = gsl_matrix_get(eigenvec,cc,rr);	
		}
	}
		
	/* Calculate inverse of eigenvector matrix for rotation back to original basis */
	
	inv_eigenvec = gsl_matrix_alloc(nbin, nbin);
	
	permutation = gsl_permutation_alloc(nbin);
	gsl_linalg_LU_decomp(eigenvec, permutation, &signum);
	gsl_linalg_LU_invert(eigenvec, permutation, inv_eigenvec);
	gsl_permutation_free(permutation);
	
	/* Place the inverted eigenvector matrix in the wldata structure. 
	   Transpose in the process so that eigenvectors run along rows. */
	
	allocate_double_matrix_contiguous(&(wldata->inv_eigenvec), nbin, nbin);
	
	for (rr=0; rr<nbin; rr++)
	{
		for (cc=0; cc<nbin; cc++)
		{
			wldata->inv_eigenvec[rr][cc] = gsl_matrix_get(inv_eigenvec,cc,rr);	
		}
	}
	
	/* Allocate vectors for rotated data */
	
	wldata->rot_modelshear = sci_dvector(nbin);
	wldata->rot_shear = sci_dvector(nbin);
	wldata->rot_shearerr = sci_dvector(nbin);
	
	/* Rotate the data and errors */
	
	for (rr=0; rr<nbin; rr++)
	{
		wldata->rot_shear[rr] = 0.0;	
		for (cc=0; cc<nbin; cc++)
		{
			wldata->rot_shear[rr] += wldata->eigenvec[rr][cc]*wldata->shear[cc];	
		}
		
		wldata->rot_shearerr[rr] = 1.0/sqrt(fabs(wldata->eigenval[rr]));
	}
		
	/* Finished */
	
	free_double_matrix_contiguous(&img);
	gsl_vector_free(eigenval);
	gsl_matrix_free(cov);
	gsl_matrix_free(inv_cov);
	gsl_matrix_free(eigenvec);
	gsl_matrix_free(inv_eigenvec);
	
	fits_close_file(fptr, &close);
	return(0);
	
}

int rotate_lensing_model(struct weaklensing *wldata)
{
	long rr, cc;
		
	for (rr=0; rr<wldata->nwl; rr++)
	{
		wldata->rot_modelshear[rr] = 0.0;	
		for (cc=0; cc<wldata->nwl; cc++)
		{
			wldata->rot_modelshear[rr] += wldata->eigenvec[rr][cc]*wldata->modelshear[cc];	
		}
	}
	
	return(0);
}

int read_sz_covariance_and_rotate_data(const char *fits_file, struct sz *szdata)
{
	fitsfile *fptr;
	double **img;
	int status=0, anynul=0, check=0, iext, signum, naxis;
	long naxis1, naxis2, nbin, nfill;
	long rr, cc;
	long fpixel[2] = {1, 1};
	gsl_vector *eigenval;
	gsl_matrix *cov, *inv_cov, *eigenvec, *inv_eigenvec;
	gsl_permutation *permutation;
	gsl_eigen_symmv_workspace *eigen_workspace;

	if (fits_file == NULL)
	{
		fprintf(stderr,"input pointer to the fits_file is NULL\n");
		return(1);
	}

	/* Open the fits file */

	fits_open_file(&fptr, fits_file, READONLY, &status);

	fits_get_hdu_num(fptr, &iext);
	if (iext != 1)
	{
		fits_movabs_hdu(fptr, 1, NULL, &status);
	}

	/* Read in relevant keywords */

	fits_read_key(fptr, TINT, "NAXIS", &naxis, NULL, &status);
	fits_read_key(fptr, TLONG, "NAXIS1", &naxis1, NULL, &status);
	fits_read_key(fptr, TLONG, "NAXIS2", &naxis2, NULL, &status);	
	
	/* Error check */
	
	if (status)
	{
		fits_report_error(stderr, status);
		fits_close_file(fptr, &close);
		return(1);
	}
	if (naxis != 2)
	{
		fprintf(stderr,"covariance must be two-dimensional\n");
		fits_close_file(fptr, &close);
		return(1);
	}
	if (naxis1 != naxis2)
	{
		fprintf(stderr,"covariance must be a square matrix\n");
		fits_close_file(fptr, &close);
		return(1);
	}
	if (naxis1 != szdata->nsz)
	{
		fprintf(stderr,"covariance is %i x %i but data is %i elements long.\n", 
						naxis1, naxis2, szdata->nsz);
		fits_close_file(fptr, &close);
		return(1);
	}
	
	nbin = naxis1;
	nfill = nbin*nbin;
	
	/* Read in data */
		
	allocate_double_matrix_contiguous(&img, nbin, nbin);

	fits_read_pix(fptr, TDOUBLE, fpixel, nfill, NULL, &(img[0][0]), &anynul, &status);
	if (status)
	{
		fits_report_error(stderr, status);
		status = 0;
		fits_close_file(fptr, &status);
		return(1);
	}
	
	/* Check that the diagonal elements are equal to the input errors squared */
	
	for (rr=0; rr<nbin; rr++)
		check += abs(gsl_fcmp(sqrt(img[rr][rr]), szdata->sze[rr], 0.0001));
	
	if (check != 0)
	{
		printf("diag(cov) not equal to input errors squared, check = %i.\n", check);
		return(1);
	}
	
	/* Place data into gsl_matrix */
	
	cov = gsl_matrix_alloc(nbin, nbin);
	
	for (rr=0; rr<nbin; rr++)
	{	
		for (cc=0; cc<nbin; cc++)
		{
			gsl_matrix_set(cov, rr, cc, img[rr][cc]);
		}
	}
		
	/* Invert the covariance matrix */
	
	inv_cov = gsl_matrix_alloc(nbin, nbin);
	
	permutation = gsl_permutation_alloc(nbin);
	gsl_linalg_LU_decomp(cov, permutation, &signum);
	gsl_linalg_LU_invert(cov, permutation, inv_cov);
	gsl_permutation_free(permutation);
	
	/* Calculate eigenvalues and eigenvectors of covariance matrix */
	
	eigenval = gsl_vector_alloc(nbin);	
	eigenvec = gsl_matrix_alloc(nbin, nbin);
	
	eigen_workspace = gsl_eigen_symmv_alloc(nbin);
	gsl_eigen_symmv(inv_cov, eigenval, eigenvec, eigen_workspace);
	gsl_eigen_symmv_free(eigen_workspace);
	
	/* Place the eigenvalues into szdata structure */
	
	szdata->eigenval = sci_dvector(nbin);
	
	for (rr=0; rr<nbin; rr++)
	{
		szdata->eigenval[rr] = gsl_vector_get(eigenval,rr);
	}
	
	/* Place the eigenvector matrix in the szdata structure.
 	   Transpose in the process so that eigenvectors run along rows. */

	allocate_double_matrix_contiguous(&(szdata->eigenvec), nbin, nbin);
	
	for (rr=0; rr<nbin; rr++)
	{
		for (cc=0; cc<nbin; cc++)
		{
			szdata->eigenvec[rr][cc] = gsl_matrix_get(eigenvec,cc,rr);	
		}
	}
		
	/* Calculate inverse of eigenvector matrix for rotation back to original basis */
	
	inv_eigenvec = gsl_matrix_alloc(nbin, nbin);
	
	permutation = gsl_permutation_alloc(nbin);
	gsl_linalg_LU_decomp(eigenvec, permutation, &signum);
	gsl_linalg_LU_invert(eigenvec, permutation, inv_eigenvec);
	gsl_permutation_free(permutation);
	
	/* Place the inverted eigenvector matrix in the szdata structure. 
	   Transpose in the process so that eigenvectors run along rows. */
	
	allocate_double_matrix_contiguous(&(szdata->inv_eigenvec), nbin, nbin);
	
	for (rr=0; rr<nbin; rr++)
	{
		for (cc=0; cc<nbin; cc++)
		{
			szdata->inv_eigenvec[rr][cc] = gsl_matrix_get(inv_eigenvec,cc,rr);	
		}
	}
	
	/* Allocate vectors for rotated data */
	
	szdata->rot_szy = sci_dvector(nbin);
	szdata->rot_sze = sci_dvector(nbin);
	szdata->rot_modely = sci_dvector(nbin);
	
	/* Rotate the data and errors */
	
	for (rr=0; rr<nbin; rr++)
	{
		szdata->rot_szy[rr] = 0.0;	
		for (cc=0; cc<nbin; cc++)
		{
			szdata->rot_szy[rr] += szdata->eigenvec[rr][cc]*szdata->szy[cc];	
		}
		
		szdata->rot_sze[rr] = 1.0/sqrt(fabs(szdata->eigenval[rr]));
	}
		
	/* Finished */
	
	free_double_matrix_contiguous(&img);
	gsl_vector_free(eigenval);
	gsl_matrix_free(cov);
	gsl_matrix_free(inv_cov);
	gsl_matrix_free(eigenvec);
	gsl_matrix_free(inv_eigenvec);
	
	fits_close_file(fptr, &close);
	return(0);
	
}

int rotate_sz_model(struct sz *szdata)
{
	long rr, cc;
		
	for (rr=0; rr<szdata->nsz; rr++)
	{
		szdata->rot_modely[rr] = 0.0;	
		for (cc=0; cc<szdata->nsz; cc++)
		{
			szdata->rot_modely[rr] += szdata->eigenvec[rr][cc]*szdata->modely[cc];	
		}
	}
	
	return(0);
}
