#include <jaco.h>
#include <fitsio.h>
#include <string.h>
#include <math.h>
#include <stdio.h>
#include <standalone.h>
#include <gsl/gsl_randist.h>

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
      readpha(fitsname[i],&phalist[i],0);
      js->fitcomp[i] = sci_strcat(fitsname[i],".dat");

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
	  phalist[i].sig[j] = sqrt(pow(phalist[i].sig[j],2.)+
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
   
    if (strcmp(js->ximage,"NONE")) {
      k = ndata;
      ndata += js->nsoftdata;
      sci_dvector_resize(&data->bigdata,ndata);
      sci_dvector_resize(&data->bigerr,ndata);
      for (i = 0 ; i < js->xnx; ++i)
	for (j = 0; j < js->xny; ++j)  {
	  if (js->softexp[i][j] > 0 && js->softimage[i][j] > SOFTIMAGE_MIN) {
	    data->bigdata[k] = js->softimage[i][j];
	    // Chi-square Gehrels error
	    data->bigerr[k] = 0.5*(1+sqrt(js->softimage[i][j]+0.75));
	    k++;
	  }
	}
    }
 
  } else if (calcsz && js->dofast) {
    
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
    ndata += wlcount;
    sci_dvector_resize(&data->bigdata,ndata);
    sci_dvector_resize(&data->bigerr,ndata);

    for (i = 0; i < (int)wlcount; ++i) {
      //ratio = wlr1[i]/wlr2[i];
      //x = ratio*ratio;
      //data->wldata->wlrad[i] = (2./3)*(wlr2[i]-x*wlr1[i])/(1.-x);
      
      data->wldata->wlrad[i] /= 60.;
      data->wldata->shearerr[i] *= js->wllss;

      data->bigdata[k] = data->wldata->shear[i];
      data->bigerr[k++] = data->wldata->shearerr[i];
    }
    data->wldata->nwl = wlcount;
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


    for (i = 0; i < szcount; ++i) 
      /* if (data->szdata->szx[i] >= js->sz0) { */
      {
	data->bigdata[k] = data->szdata->szy[i];
	data->bigerr[k++] = data->szdata->sze[i];
      }
      /* } */

  } else
#endif
    data->szdata = NULL;

  if (ndata) {
    data->ndata = ndata;
    data->bigmodel = sci_dvector(ndata);
  } else
    BYE("Nothing to do! (all 'mode' flags were turned off)\n");

  return 0.;
}

#define CHANGEDPAR(x) (DIFFER(params[(x)].value,oldparams[(x)].value))

int get_bigmodel(double *x,double *params,double *model, double *errors,
		 double *prior, unsigned long ndata, void *datapointer)
{
  int i,j,k,l;
  FILE *fitfile=NULL;
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

  phalist = data->phalist;
  arflist = data->arflist;
  rmflist = data->rmflist;
  elimits = data->elimits;
  spectrum = data->spectrum;

  // No Bayesian prior
  (*prior) = 0.;
  
  JS(calcxray); JS(calcsz); JS(xraycount); JS(instrument); 
  JS(r1full); JS(annuli);

  int ndatasets = 0;
  if (calcxray || calcsz) {
    for (i = 0; i < xraycount; ++i) {
      
      jaco_calc_models(instrument[i],annuli[i],
		       params,&arflist[i].nebins,spectrum,arflist[i].elo,
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
	  fitfile = fopen(js->fitcomp[i],"w+");
	else
	  fitfile = NULL;

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
	  }
	  /* diff = (s->mcounts[j]-s->flux[j])/s->sig[j]; */
	  /* totchisq += diff*diff; */
	  /* printf("%7E %7E\n",s->mcounts[j],totchisq); */
	  
	  if (!data->fitannuli) 
	    model[npoints+j] = s->mcounts[j];
	  
	  if (fitfile) {
	    diff =  elimits[i].emax[s->chanhi[j]]-
	      elimits[i].emin[s->chanlo[j]];
	    fprintf(fitfile,"%E %E %E %E\n",
		    (elimits[i].emin[s->chanlo[j]]+
		     elimits[i].emax[s->chanhi[j]])/2.,
		    s->flux[j]/diff,s->sig[j]/diff,s->mcounts[j]/diff);
	  }

	  if (js->dosim > 0) {
	    fname = sci_strcat("sim.",s->filename);
	    s->mcounts[j] += gsl_ran_gaussian(js->rng,s->sig[j]);
	    writepha(fname,s);
	    free(fname);
	  }
	}
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
      ++ndatasets;
    }

    if (js->softimage != NULL) {
      for (i = 0; i < js->xnx; ++i)
	for (j = 0; j < js->xny; ++j)
	  if (js->softexp[i][j] > 0 && js->softimage[i][j] > SOFTIMAGE_MIN) {
	    model[npoints] = js->softmodel[i][j]*js->ximnorm+js->ximback;
	    ++npoints;
	  }
    }
  }

  if (wldata) {

    if (js->fitwrite == 1) {
      fname = sci_strcat(js->profname,".wlfit.dat");
      fitfile = fopen(fname,"w+");
      free(fname);
    }

    jaco_calc_models(WL,ndatasets,params,&wldata->nwl,wldata->modelshear,wldata->wlrad,js);

    for (i = 0; i < wldata->nwl; ++i) {

      // Moved to hrothgar:
      diff = (wldata->modelshear[i]-wldata->shear[i])/wldata->shearerr[i];
      //if (chiarr) chiarr[npoints] = diff;
      totchisq += diff*diff;

      model[npoints] = wldata->modelshear[i];
      if (fitfile)
	fprintf(fitfile,"%E %E %E %E\n",wldata->wlrad[i],wldata->shear[i],
		wldata->shearerr[i],model[npoints]);

      ++npoints;
    }
    if (fitfile)
      fclose(fitfile);

    ++ndatasets;
  }
  //printf("WL %f %ld\n",finalchisq,npoints);

  if (vels) {

    if (js->fitwrite == 1) {
      fname = sci_strcat(js->profname,".vels.dat");
      fitfile = fopen(fname,"w+");
      free(fname);
    }

    jaco_calc_models(VELS,ndatasets,params,&vels->nvel,vels->modeldisp,vels->rad,js);

    for (i = 0; i < vels->nvel; ++i) {

      // Moved to hrothgar:
      diff = (vels->modeldisp[i]-vels->disp[i])/vels->disperr[i];
      //if (chiarr) chiarr[npoints] = diff;
      totchisq += diff*diff;

      model[npoints] = vels->modeldisp[i];
      if (fitfile)
	fprintf(fitfile,"%E %E %E %E\n",vels->rad[i],vels->disp[i],
		vels->disperr[i],model[npoints]);

      ++npoints;
    }
    if (fitfile) 
      fclose(fitfile);

    ++ndatasets;
  }

  if (szdata) {
    jaco_calc_models(SZ,ndatasets,params,&szdata->nsz,szdata->modely,NULL,js);

    if (js->fitwrite == 1) {
      fname = sci_strcat(js->profname,".szfit.dat");
      fitfile = fopen(fname,"w+");
      free(fname);
    }

    for (i = 0; i < szdata->nsz; ++i) 
      /* if (szdata->szx[i] >= data->js->sz0)  */
	{

	// Moved to hrothgar:
	diff = (szdata->modely[i]-szdata->szy[i])/szdata->sze[i];
	//if (chiarr) chiarr[npoints] = diff;
	totchisq += diff*diff;
	if (fitfile)
	  fprintf(fitfile,"%d %E %E %E\n",i,szdata->szy[i],
		  szdata->sze[i],szdata->modely[i]);

	model[npoints] = szdata->modely[i];
	++npoints;
      }

    if (fitfile)
      fclose(fitfile);
    ++ndatasets;
  }

    
  return 0;
}
