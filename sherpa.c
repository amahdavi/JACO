#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <math.h>
#include <time.h>
#include <unistd.h>
#include <gsl/gsl_rng.h>
#include <gsl/gsl_randist.h>

#define JACO_INIT
#include <jaco.h>

  char **ng_paramnames, **ng_paramunits;
  double ng_parmin[NTOTPARAMS], ng_parmax[NTOTPARAMS], 
    ng_initvalues[NTOTPARAMS];
  int ng_dofit[NTOTPARAMS];

#ifdef NONGRAVITY
   
void hrothgar_init_nongravity(struct hrothgar_setup *setup,
			      struct jaco_state *js)
{

  int i,j=0;

  ng_paramnames = (char **)malloc(NTOTPARAMS*sizeof(char *));
  ng_paramunits = (char **)malloc(NTOTPARAMS*sizeof(char *));

  // Generate parameter strings as needed.
  for (i = 0; i < 100; ++i) {

    int good=0;
    while (!good)
      switch(j) {
	
      case REDSHIFT: 
      case CONTRAST: 
      case NONTHERMAL: 
      case TEMPBIAS: 
      case EXT: 
      case EXZ:  
      case EXNORM:  
      case W0: 
      case WR:   
      case PLSLOPE0:   
      case PLSLOPE1:   
      case PLSLOPE2:   
      case PLSLOPE3:   
      case PLSLOPE4:   
      case PLNORM0: 
      case PLNORM1: 
      case PLNORM2: 
      case PLNORM3: 
      case PLNORM4: 
      case BACKT: 
      case BACKZ: 
      case BNORM0: 
      case BNORM1: 
      case BNORM2: 
      case BNORM3: 
      case BNORM4: 
	ng_paramnames[j] = paramnames[j];
	ng_paramunits[j] = paramunits[j];
	ng_parmin[j] = parmin[j];
	ng_parmax[j] = parmax[j];
	ng_initvalues[j] = initvalues[j];
	ng_dofit[j] = dofit[j];
	++j; good=0; break;
	
      default: good=1; break;
	
      }

    js->Tpar[i] = j;
    ng_paramnames[j] = (char *)malloc(10*sizeof(char));
    ng_paramunits[j] = (char *)malloc(10*sizeof(char));
    sprintf(ng_paramnames[j],"T%03d",i);
    sprintf(ng_paramunits[j],"keV");
    if (i >= 0) {
      ng_parmin[j] = 0.3;
      ng_parmax[j] = 18;
      ng_initvalues[j] = 2.;
    } else {
      ng_parmin[j] = 0.5;
      ng_parmax[j] = 2;
      ng_initvalues[j] = 1.;
    }
    ng_dofit[j] = (i != 0);
    ++j;
    good=0;

  }

  for (i = 0; i < 100; ++i) {
    js->Zpar[i] = j;
    ng_paramnames[j] = (char *)malloc(10*sizeof(char));
    ng_paramunits[j] = (char *)malloc(10*sizeof(char));
    sprintf(ng_paramnames[j],"Z%03d",i);
    sprintf(ng_paramunits[j],"Zsun");
    ng_parmin[j] = 0.0001;
    ng_parmax[j] = 2.5;
    ng_initvalues[j] = 0.3;
    ng_dofit[j] = (i != 0);
    ++j;
  }

  for (i = 0; i < 100; ++i) {
    js->nepar[i] = j;
    ng_paramnames[j] = (char *)malloc(10*sizeof(char));
    ng_paramunits[j] = (char *)malloc(10*sizeof(char));
    sprintf(ng_paramnames[j],"n%03d",i);
    sprintf(ng_paramunits[j],"cm^(-3)");
    ng_parmin[j] = -5;
    ng_parmax[j] = 0.;
    ng_initvalues[j] = -3;
    ng_dofit[j] = (i != 0);
    ++j;
  }

  hrothgar_init_pars(setup,NTOTPARAMS,ng_paramnames,ng_parmin,
		     ng_parmax,ng_initvalues,ng_dofit,ng_paramunits);

}

#endif

// Set state variables to default values, including nulling pointers
// which should be NULL at startup.
// Read the configuration file as appropriate
void jaco_init(struct jaco_state *js, 
	       struct hrothgar_setup *setup, int argc, char **argv) 
{
  int  nwllist;
  char **wllist,*endptr;

  js->setup = setup;

  js->interp_coefs = NULL; js->rspline = NULL; 
  js->nannuli = js->lastszannulus = js->tempwarn = 
    js->metalwarn = js->specwarn=0;
  js->nlastbin = js->standalone = js->mpirank = 0;
  js->calcsz  = js->calcwl = 0;
  js->samepars = 1; js->nrebin = 0;
  js->background = NULL;

  hrothgar_init_welcomestring(setup,WELCOMESTRING);
  hrothgar_init_stringpars(setup,NSTRINGPARAMS,stringparname,stringparval,
			   stringparcomments);

  #ifdef NONGRAVITY
  
  hrothgar_init_nongravity(setup,js);

  #else

  hrothgar_init_pars(setup,NTOTPARAMS,
		     paramnames,parmin,parmax,initvalues,
		     dofit,paramunits);
  #endif

  js->mpirank = hrothgar_init(setup,argc,argv);

  if (js->mpirank < 0) exit(0);
  js->nfitparams = setup->nfitparams;

  js->efit1 = atof(stringparval[E1]);
  js->efit2 = atof(stringparval[E2]);
  js->syserr = atof(stringparval[SYSERR]);
  js->calcxray = atoi(stringparval[XRAYMODE]);
  js->calcwl = atoi(stringparval[WLMODE]);
  js->calcsz = atoi(stringparval[SZMODE]);
  js->calcvel = atoi(stringparval[VELMODE]);
  js->backlist = stringparval[BACKLIST];
  js->fitmc = atoi(stringparval[FITMC]);
  js->fitfgas = atoi(stringparval[FITFGAS]);


  js->cutradius = atof(stringparval[RCUT]);

  js->jacofilename = stringparval[XRAYDATA];
  js->spectralcode = stringparval[SPECTRALCODE];

  js->psflist = stringparval[XRAYPSF];
  js->wldata=stringparval[WLDATA];
  js->veldata=stringparval[VELDATA];
  js->fitwrite = atoi(stringparval[WRITE]);
  js->dosim = atoi(stringparval[DOSIM]);
  js->axis = atoi(stringparval[AXIS]);
  js->simfile = stringparval[SIMFILE];
  js->massunit = atof(stringparval[SIMMASS]);
  js->mpc = atof(stringparval[SIMLENGTH]);
  js->debug = atoi(stringparval[DEBUG]);
  js->ximage = stringparval[XIMAGE];
  js->xime1 = atof(stringparval[XIMAGEE1]);
  js->xime2 = atof(stringparval[XIMAGEE2]);
  js->lxe1 = atof(stringparval[LXE1]);
  js->lxe2 = atof(stringparval[LXE2]);
  js->precision = atof(stringparval[PRECISION]);
  js->dofast = !atoi(stringparval[FINEGRID]);
  js->gasmodel = (int)floor(atof(stringparval[GASMODEL]));
  js->darkmodel = (int)floor(atoi(stringparval[DARKMODEL]));
  if (js->gasmodel < 0) {
    js->totalmass = 1;
    js->gasmodel = -js->gasmodel;
  } else
    js->totalmass = 0;

  if (js->gasmodel == TRIAXIAL) {
    if (js->dofast)
      BYE("Triaxial mode requires finegrid=1.");
    if (strcmp(js->ximage,"NONE") == 0 && js->calcxray)
      BYE("X-ray triaxial mode requires input images and exposure maps.");
  }

  float xpixscale,wcs[3][2];
  unsigned long xnx, xny;
  int i, j;
  char **imglist;
  if (js->dosim > 1)
    init_sim(js);
  else
    if (strcmp(js->ximage,"NONE")) {
      parse_list(js->ximage,',',&imglist);
      js->softimage = readimage(imglist[0],&xpixscale,&xnx,&xny,wcs);
      js->softmodel = sci_dmatrix(xnx,xny);
      if (js->softimage == NULL) 
	BYE("Could not open soft X-ray photon image");
      js->softexp  = readimage(imglist[1],&js->xpixscale,&js->xnx,&js->xny,js->wcs);
      if (js->softexp == NULL) 
	BYE("Could not open soft X-ray exposure mask");
      if (js->xnx != xnx|| js->xny != xny)
	BYE("Exposure map size does not match photon image scaling: %ld %ld %ld %ld",
	    js->xnx,xnx,js->xny,xny);
      for (i = 1; i < 3; ++i)
	for (j = 0; j < 2; ++j) {
	  js->wcs[i][j] *= PIOVER180;
	  printf("%d %d %E %E\n",i,j,js->wcs[i][j],wcs[i][j]);
	}
      js->nsoftdata = 0;
      double ra,dec;
      for (i = 0; i < js->xnx; ++i)
	for (j = 0; j < js->xny; ++j)
	  if (js->softexp[i][j]> 0 && js->softimage[i][j] > SOFTIMAGE_MIN) 
	    ++js->nsoftdata;
      

      sincos(js->wcs[1][0],&js->sinra0,&js->cosra0);
      sincos(js->wcs[1][1],&js->sindec0,&js->cosdec0);
      js->xymodel = NULL;
    } else 
      js->softimage = NULL;

  js->ndelta = parse_list_double(stringparval[DELTAS],',',&js->deltas);
  nwllist = parse_list(stringparval[WLCONFIG],',',&wllist);

  if (nwllist < 0 || nwllist > 3)
    BYE("Error parsing weak lensing parameters %s",stringparval[WLCONFIG]);
  
  js->wlbeta = strtod(wllist[0],&endptr);

  if (endptr[0]) 
    BYE("Error parsing weak lensing parameters %s",stringparval[WLCONFIG]);

  if (nwllist < 2) 
    js->wlcorrfac = 0;
  else 
    js->wlcorrfac = strtod(wllist[1],&endptr)/(js->wlbeta*js->wlbeta)-1.;


  if (endptr[0]) 
    BYE("Error parsing weak lensing parameters %s",stringparval[WLCONFIG]);

  if (nwllist < 3) 
    js->wllss = 1.;
  else 
    js->wllss = sqrt(strtod(wllist[2],&endptr));

  if (endptr[0]) 
    BYE("Error parsing weak lensing parameters %s",stringparval[WLCONFIG]);


  js->szdata=stringparval[SZDATA];
  js->szconfig=stringparval[SZCONFIGS];
  js->tprojpower=atof(stringparval[TPROJPOWER]);
  js->hubble = atof(stringparval[HUBBLE]);
  js->omegam = atof(stringparval[OMEGA_M]);
  js->omegal = atof(stringparval[OMEGA_L]);

  if (js->xrayra < 0. && js->calcsz)
    jaco_printf("X-ray center not specified. Will assume SZ center = X-ray center");

  if (js->fitwrite > 2) 
    system("rm -f jaco-profiles*.txt");

  js->solver = gsl_root_fsolver_alloc(gsl_root_fsolver_falsepos);

}
	   

// Read X-ray input files (including PSF files) and calculate which
// annuli are unique
int jaco_read_data(struct jaco_state *js)
{
  int m,annulus,newannulus,instr,npsf;
  int speccount[NINSTR];
  char **psffiles;
  long i,j=0,k,ii,jj;
  unsigned long count,psfcount[NINSTR];
  double r1sq,r2sq;

  JS(standalone); JS(jacofilename);

  if (!standalone) {
    printf("Joint Analysis of Cluster Observations v%s\n",VERSION);
    printf("Copyright (C) 2006 Andisheh Mahdavi\n");
  }
				  
  // Initialize weak lensing and SZ components, if necessary
#if defined(WITH_SZ) || defined(WITH_BOLOCAM)
  if (!standalone)
    printf(SZC);

  js->szparams.node = js->mpirank;
  js->szparams.quiet = (js->debug == 0);
  #ifdef WITH_SZ
  sprintf(js->szparams.DerotDataFile,js->szdata);
  #endif
  
  js->szparams.ndataset = parse_list(js->szconfig,',',
				     &js->szparams.paramfiles);

  if (js->calcsz) {
    #ifdef WITH_SZ
    if (big_initialize_sz(&js->szparams) != 0) js->calcsz = 0;
    #else
    if (big_initialize_sz(&js->szparams) == 0) js->calcsz = 0;
    js->szparams.ymap = sci_dtensor(js->szparams.ndataset,
				js->szparams.npix[0],
				js->szparams.npix[0]);
    #endif
  }
  if (js->calcsz && js->mpirank == 0) printf("JACO: SZ Fitting enabled.\n");

/* #else */
/*   js->szparams.ndataset = 1; */
/*   js->szparams.racent = sci_dvector(1); */
/*   js->szparams.deccent = sci_dvector(1); */
/*   js->szparams.pixsize = sci_dvector(1); */
/*   js->szparams.npix = sci_ivector(1); */
/*   js->szparams.racent[0] = js->xrayra; */
/*   js->szparams.deccent[0] = js->xraydec; */
/*   js->szparams.npix[0] = 256; */
/*   js->szparams.pixsize[0] = 0.1/60.; */
/*   js->szparams.ymap = sci_dtensor(2,js->szparams.npix[0],js->szparams.npix[0]); */
#endif



  if (js->calcwl && 
      js->mpirank == 0) printf("JACO: Weak Lensing Fitting enabled.\n");

  if (!standalone) {
    printf("JACO comes with ABSOLUTELY NO WARRANTY. This is free\n");
    printf("software; you are welcome to redistribute it under certain\n");
    printf("conditions. See the documentation for details.\n");
  }

  // Read in the X-ray cooling function and initialize interpolators for it
  if (js->calcxray || (js->calcsz && js->dofast)) {
    read_cooling_function(js);
  
    // Read the inner and outer radii of the spectra
    count = read_ascii(js->jacofilename,"%s %f %f %f",NULL,
		       &js->fitsname,&js->r1full,&js->r2full,
		       &js->goodarea);
    js->xraycount = (int)count;
    
/* >>>>>>> 1.119 */

    JS(r1full); JS(r2full);

    js->annuluscount = sci_ivector(NINSTR);

    for (i = 0; i < NINSTR; ++i) speccount[i] = js->annuluscount[i] = 0;
    
    // Allocate space for the various annulus data
    js->r1list = sci_dvector(count+1);
    js->r2list = sci_dvector(count+1);
    js->effrad = sci_dvector(count);
    js->ringarea = sci_dvector(count);
    js->annuli = sci_ivector(count);
    js->annulusid = sci_imatrix(NINSTR,count);
    js->psfrow = sci_imatrix(NINSTR,count);
    js->instrument = sci_ivector(count);
    
    // Calculate the number of unique annuli for the whole dataset
    // And the number of unique annuli per instrument
    js->nannuli = 0;
    for (i = 0; i < js->xraycount; ++i) {
      
      // Use the spectrum file name to figure out which instrument it comes from
      
      if (strstr(js->fitsname[i],"/") != NULL) 
	BYE("Spectrum file %s contains / characters---not allowed.",
	    jacofilename);
      
      instr = -1;
      for (j = 0; j < NINSTR; ++j)
	if (strstr(js->fitsname[i],idcode[j])) instr = j;
      if (instr < 0)
	BYE("Unrecognized instrument in %s.",js->fitsname[i]);
      
      js->instrument[i] = instr;
      
      // See if this annulus has been seen before.
      newannulus = 1;
      annulus = -1;
      for (j = 0; j < i && newannulus; ++j) {
	newannulus = DIFFER(r1full[i],r1full[j])+DIFFER(r2full[i],r2full[j]);
	if (!newannulus) annulus = js->annuli[j];
      }
      
      js->annuli[i] = annulus;
      // If this is a new annulus; calculate/store its area and effective radius
      if (newannulus) {
	annulus = js->nannuli;
	js->annuli[i] = annulus;
	js->r1list[annulus] = r1full[i];
	js->r2list[annulus] = r2full[i];
	r1sq = r1full[i]*r1full[i];
	r2sq = r2full[i]*r2full[i];
	js->ringarea[annulus] = (r2sq-r1sq);
	js->effrad[annulus] = 2.*(r2sq*r2full[i]-r1sq*r1full[i])/
	  (3.*(js->ringarea[annulus]));
	js->ringarea[annulus] *= PI;
	++js->nannuli;
      }
      
      // See if this annulus has been seen before for this particular instrument
      newannulus = 1;
      for (j = 0; j < js->annuluscount[instr] && newannulus; ++j) {
	k = js->annulusid[instr][j];
      newannulus = DIFFER(r1full[i],js->r1list[k])+
	DIFFER(r2full[i],js->r2list[k]);
      }
      // If it has not been seen, then add it to the list.
      if (newannulus) {
	js->annulusid[instr][js->annuluscount[instr]] = annulus;
	js->psfrow[instr][annulus] = js->annuluscount[instr];
	++js->annuluscount[instr];
	//printf("New annulus %d for %s (%d)\n",annulus,idcode[instr],
	//     js->psfrow[instr][annulus]);
      }
    }
    js->r1list[js->nannuli] = js->r2list[js->nannuli-1];
    js->r2list[js->nannuli] = 3.*js->r1list[js->nannuli];
    
    // Geometric factors for spherical projection
    js->geomfactors = sci_dmatrix(js->nannuli+1,js->nannuli+1);
    for (i = 0; i < js->nannuli+1; ++i) {
      for (j = 0; j < js->nannuli+1; ++j) 
	if (j > i-2) {
	  //printf("%2d %2d ",i,j);
	  js->geomfactors[i][j] = 
	    pow(js->r2list[j]*js->r2list[j]-
		js->r1list[i]*js->r1list[i],1.5)/3.;
	  //printf("%.1E   ",js->geomfactors[i][j]);
	}
      //printf("\n");
    }
    
    // Read in the relevant psf files
    k = 0;
    js->psf = (double **)malloc(3*NINSTR*sizeof(double *));
    npsf = parse_list(js->psflist,',',&psffiles);
    for (m = 0; m < NINSTR; ++m) 
      if (js->annuluscount[m] > 0) {
	long npsfrows = js->annuluscount[m]*js->annuluscount[m];
	
	for (j = 0; j < npsf && !strstr(psffiles[j],idcode[m]); ++j);
	
	if (j == npsf) 
	  jaco_printf("No PSF file for instrument %s... assuming "	\
		      "a perfect response\n",idcode[m]);
	
	for (i = 1; i <= 3; ++i) 
	  if (j == npsf) {
	    js->psf[k] = sci_dvector(npsfrows);
	    for (ii = 0; ii < js->annuluscount[m]; ++ii)
	      for (jj = 0; jj < js->annuluscount[m]; ++jj)
		js->psf[k][js->annuluscount[m]*ii+jj] = 1.*(ii == jj)*(i==1);
	    k++;
	  } else {
	    js->psf[k++] = double_readdata(psffiles[j],i,&psfcount[m],0);
	    if (psfcount[m] != npsfrows) {
	      BYE("   Mismatch of PSF (%ld) and data files(%d).\n",
		  psfcount[m],js->annuluscount[m]);
	    }
	  }
      }
      else
	for (i = 1; i <= 3; ++i) 
	  js->psf[k++] = NULL;
    
  } else  
    js->xraycount = 0;
  //else
  //printf("   --%ld verified data points observed by %s (%ld annuli) \n",
  //     speccount[m],detector[m],annuluscount[m]);
  //printf("Done.\n");

#ifndef NONGRAVITY
  if (js->fitwrite > 1) {
    const char *profparams[INFO_SIZE] = 
      {"Reff","Tproj","Md","Mg","Mstar",
       "fgas","T3D","n_e","pressure","entropy",
       "coolingtime","Tshock","LX","metal","Mtot"};
    
    const char *contrparams[DELTA_SIZE] = 
      {"Delta","rDelta","cdelta","MgDelta","MdDelta",
       "MsDelta","MtDelta","MprojDelta"};
    
    js->setup->ninfo = js->ninfo = INFO_SIZE*js->nannuli+DELTA_SIZE*js->ndelta;
    js->setup->infoname = js->infoname = sci_stringmatrix(js->ninfo,100);
    js->setup->infoarray = js->infoarray = 
      sci_dvector(js->ninfo);
    
    k = 0;
    for(i = 0; i < js->nannuli; ++i) 
      for (j = 0; j < INFO_SIZE; ++j) 
	sprintf(js->infoname[k++],"prof%s%03ld",profparams[j],i);
    
    for(i = 0; i < js->ndelta; ++i)
      for (j = 0; j < DELTA_SIZE; ++j) {
	sprintf(js->infoname[k++],"contr%s%03ld",contrparams[j],i);
      }
  }
#else
  js->setup->ninfo = js->ninfo = 0;  
#endif

  return 0;

}


// Calculate various physical models and return spectra/SZ decrements/shear
int jaco_calc_models(int instr, 
		     int  annulus,
		     double *params, int *ndata, double *outValues,
		     double *x0, struct jaco_state *js)
{

  double c0,c1,c2,nclock,logT,logZ,normspec;
  static unsigned long int elem,nrun=0;
  int i,j,k,l,samepars=1,sameeverything=1,coresame=0,newmodel;
  int recalcback = 0,id;
  static double oldp[NTOTPARAMS];
  static clock_t cl0=0,cl1=0,runtime=0,timenow;
  static int donexray=0,nolddata=0,firstinstr=-1;
  static double *finalspectrum=NULL,*backspec=NULL,
    *corespectrum=NULL,**absfac=NULL,*plbackspec=NULL,
    *softbackground=NULL,**spectrum=NULL,*background=NULL,
    *savespectrum=NULL, **singleTspec=NULL, 
    *currT=NULL,*currZ=NULL,*currn=NULL;
  static double dener=0.,volcubed;
  int ebins=-1;
  struct jaco_vars jv;

  static double *savesz = NULL;

  // If the basic physical parameters have not changed since the last call,
  // pass that information on to the main program

  JS(debug);
  JS(gasmodel);
  js->tbias = params[TEMPBIAS];

#ifdef NONGRAVITY

  for (i = 0; i < NTOTPARAMS; ++i) {
    
    // Do the checking only if neither parameter is close to zero
    coresame = sameeverything = 0;
    js->changedpar[i] = (fabs(params[i]-oldp[i]) > 0.);
    samepars *= 1.-js->changedpar[i];
    
    //printf("%20.10E",params[i]);
    // Save copies of the current parameters for later checking
    oldp[i] = params[i];
  }
  //samepars = 0;
  //printf("%d %.0f\n",samepars,params[ANNULUS]);

#else
  for (i = 1; i < NTOTPARAMS; ++i) {
    
    // Do the checking only if neither parameter is close to zero
    js->changedpar[i] = (fabs(params[i]-oldp[i]) > 0.);
    
    //printf("%20.10E",params[i]);
    // Save copies of the current parameters for later checking
    oldp[i] = params[i];
  }
  coresame = js->changedpar[EXNORM]*js->changedpar[EXT]*js->changedpar[EXZ];

#endif

  samepars = 0;

  newmodel = !annulus;
  if (newmodel) donexray = 0;
  
  // Initialize the various model parameters using the given inputs.
  if (newmodel) 
    jaco_update_models(params,js);

  jv.state = js;

#if defined(WITH_SZ) || defined(WITH_BOLOCAM)

  // The SZ model is calculated during the X-ray integration process
  // and saved in the savesz array.  There is currently no way to have
  // the SZ model on its own.
  if (instr == SZ) {
    if (js->szparams.nest != *ndata) {
      BYE("ERROR -- mismatch in length of SZ data. %d vs. %d\n",
	     js->szparams.nest,*ndata);
    }
    if (js->dofast) {
      if (savesz == NULL || newmodel || !donexray) {
	BYE("ERROR --- in fast mode you need an X-ray model before the SZ model.\n");
      } else 
	for (i = 0; i < js->szparams.nest; ++i)
	  outValues[i] = savesz[i];
      return 0;
    }
  }
#endif  

  if (instr == WL) {
    tangential_shear(x0,outValues,*ndata,js);
    return 0;
  }

  if (instr == VELS) {
    printf("Calculating vdp\n");
    #pragma omp parallel for private(i) firstprivate(jv) schedule(dynamic)
    for (i = 0; i < (*ndata); ++i)
       outValues[i] = velocity_dispersion(x0[i],&jv);
    return 0;
  }
  
  donexray = 1;
  if (debug > 0)
    js->logfile=fopen("testlog","a+");


  js->samepars = samepars;


  // This code only gets called if an X-ray instrument is here.
  if (instr < SZ) {
    // Determine whether the energy binning of the cooling function needs
    // to be changed, and adjust parameters accordingly.
    js->binningchanged = rebincoolingfunction(x0,*ndata,js);
    ebins = js->nlastbin;
    //printf("%d %d\n",ebins,*ndata);

    //printf("OK\n %d %d\n",js->nannuli,*ndata);
    if (js->binningchanged) {
      if (spectrum) sci_free_dmatrix(spectrum,js->nannuli);
      spectrum = sci_dmatrix(js->nannuli,ebins);
    }
    finalspectrum = sci_dvector(*ndata);

    if (js->binningchanged) {
      dener = x0[1]-x0[0];
      if (corespectrum != NULL) free(corespectrum);
      if (absfac != NULL) sci_free_dmatrix(absfac,js->nannuli);
      if (plbackspec != NULL) free(plbackspec);
      if (backspec != NULL) free(backspec);
      if (softbackground != NULL) free(softbackground);
      if (background != NULL) free(background);\
      if (singleTspec != NULL) sci_free_dmatrix(singleTspec,js->nannuli+1);
      if (js->background != NULL) free(js->background);
      currn = sci_dvector(js->nannuli+1);
      currT = sci_dvector(js->nannuli+1);
      currZ = sci_dvector(js->nannuli+1);
      corespectrum = sci_dvector(ebins);
      absfac = sci_dmatrix(js->nannuli,ebins);
      singleTspec = sci_dmatrix(js->nannuli+1,ebins);
      plbackspec = sci_dvector(ebins);
      backspec = sci_dvector(ebins);
      softbackground = sci_dvector(ebins);
      background = sci_dvector(ebins);
      js->background = sci_dvector(ebins);
    }
  

    // Recalculate absorption coefficients, if necessary
    if (recalcback += (js->changedpar[W0] || js->changedpar[WR] || 
		       js->binningchanged)) {
      if (DIFFER(params[WR],0.)) {
#pragma omp parallel for private(i)  schedule(dynamic)
	for (i = 0; i < js->nannuli; ++i)  {
	  photoabs(params[W0]+params[WR]*js->effrad[i],
		   js->lastbin,dener,ebins,absfac[i]);
	}
      } else {
	if (annulus == 0 || js->binningchanged) {
	  photoabs(params[W0],js->lastbin,dener,ebins,absfac[0]);
	  for (i = 1; i < js->nannuli; ++i) {
	    for (j = 0; j < ebins; ++j) {
	      absfac[i][j] = absfac[0][j];
	      //printf("ABS %E %E\n",js->lastbin[j],absfac[0][j]);
	    }
	  }
	
	}
      }
    }
    
    // Recalculate powerlaw extragalactic background, if necessary
    /*   if (recalcback += (js->changedpar[PLSLOPE] || */
    /* 		     js->binningchanged)) { */
    /*     for (i = 0; i < ebins; ++i) */
    /*       plbackspec[i] = */
    /* 	(pow(js->lastbin[i],1.-params[PLSLOPE])- */
    /* 	 pow(js->lastbin[i]+dener,1.-params[PLSLOPE]))/(params[PLSLOPE]-1.); */
    /*   } */

    // Recalculate unnnormalized soft X-ray background spectrum,
    // if necessary
    if (recalcback += (js->changedpar[BACKT] || js->changedpar[BACKZ] || 
		       js->binningchanged)) {
      logT = log10(params[BACKT]);
      logZ = log10(params[BACKZ]);
      for (i = 0; i < ebins; ++i) {
	jv.coolmatrix = js->rebinnedcooling[i];
	backspec[i] = cooling(logT,logZ,&jv);
      }
    }
      
    jv.state = js;

    if (params[EXNORM] > 0. && (js->binningchanged || !annulus)) {
      logT = log10(params[EXT]);
      logZ = log10(params[EXZ]);
      for (j = 0; j < ebins; ++j) {
	jv.coolmatrix = js->rebinnedcooling[j];
	//corespectrum[j] = params[EXNORM]*cooling(logT,logZ,&jv);
	corespectrum[j] = params[EXNORM]*((pow(js->lastbin[j],1.-params[EXT])- pow(js->lastbin[j]+dener,1.-params[EXT])))/(params[EXT]-1); 
      
      }
    }
  }

  jv.state = js;


  if (!annulus || js->binningchanged || instr == SZ) {
    if (js->dosim > 2)
      make_simimages(js);
    

#ifdef NONGRAVITY
    double norm,esum,logn=0.,penalty;

    if (js->changedpar[REDSHIFT])
      volcubed = pow(ARCMINTORAD*js->angdist,3.);

    logZ = logT = logn = 0;
      
    for (i = 0; i < js->nannuli+1; ++i) {
      
      if ((((i-1) % js->calcxray) == 0) || i == 0 ) {
	
	logT = log10(params[js->Tpar[i]]);
	logZ = log10(metallicity(-params[js->Zpar[i]],&jv));
      }
      logn = params[js->nepar[i]];
      
      esum = 0.;
      for (j = 0; j < ebins; ++j) {
	jv.coolmatrix=js->rebinnedcooling[j];
	singleTspec[i][j] = cooling(logT,logZ,&jv);
	if (gasmodel > 0) {
	  esum += js->lastbin[j]*singleTspec[i][j];
	}
      }
      //fclose(zoom);
      // }
      
      // nepar Gives us the total luminosity in units of 1e47 Lsun
      if (gasmodel > 0) {
	esum *= jv.bolocorr;
	if (js->debug > 2) 
	  printf("Bolometric correction %E\n",jv.bolocorr);
	norm = pow(10.,logn)*1.e47/
	  (KEV*esum*FOURPI*MPC*MPC*js->lumdistancesq/js->opz);
      }
      // nepar Gives us the gas density
      else {
	norm = 1.e-14*pow(10.,2.*logn)/jv.nenh;
	norm *= js->opz*volcubed*MPC/js->lumdistancesq;
      }
      
      //printf("%d %E %E %E foo\n",i,params[js->nepar[i]],logT,logZ);
      
      for (j = 0; j < ebins; ++j)  
	singleTspec[i][j] *= norm;
    }
	

/*       penalty=0.; */
/*       double crit; */
      
/*       for (k = 1; k < 4; k += 2) { */
/* 	crit = k*0.1; */

/* 	for (i = k; i < js->nannuli-k; ++i) { */
/* 	  if ((i % js->calcxray) == 0) { */
/* 	    j = i/js->calcxray; */
	    
/* 	    norm = fabs(currT[j]-(currT[j-k]+currT[j+k])/2.); */
/* 	    if (norm > crit) penalty += norm-crit; */
	    
/* 	    norm = fabs(currZ[j]-(currZ[j-k]+currZ[j+k])/2.); */
/* 	    if (norm > 2*crit) penalty += norm-crit; */
	    
/* 	    norm = fabs(currn[j]-(currn[j-k]+currn[j+k])/2.); */
/* 	    if (norm > crit) penalty += norm-crit; */
/* 	  } */
	  
/* 	} */
/*       } */

/*       js->setup->penalty = 0*penalty/10.; */
/*       if (debug > 0) */
/* 	printf("Penalty-------------------%E\n",penalty); */

/* /\*       if (js->tprojpower > 0) *\/ */
/* /\*       for (i = 0; i < js->nannuli+1; ++i) *\/ */
/* /\*       	for (j = 1; j < ebins; ++j)  { *\/ */
/* /\* 	  if (j % 6 == 0) *\/ */
/* /\* 	    singleTspec[i][j] -= js->tprojpower*singleTspec[i][0]*penalty; *\/ */
/* /\* 	  if (j % 6 == 3) *\/ */
/* /\* 	    singleTspec[i][j] += js->tprojpower*singleTspec[i][0]*penalty; *\/ */
/* /\* 	} *\/ */
      
/*       //#pragma omp parallel for private(i,j,k,norm) firstprivate (jv)  */
      
/*       //double ** foo = sci_dmatrix(js->nannuli,js->nannuli); */
    
    double ownfrac;
    
    if (gasmodel > 0) {
      if (js->nannuli > 1) {
	BYE("Error- luminosity mode supported only for 1 annulus so far.");
      }
      for (j = 0; j < ebins; ++j) 
	spectrum[0][j] = singleTspec[0][j];
    } else {      
      for (i = 0; i < js->nannuli; ++i) {
	/* 	  printf("%3d------------------%6.0E\n",i,js->geomfactors[i][i]/ */
	/* 	      ((1./3.)*(pow(js->r2list[i],3.)- */
	/* 			  pow(js->r1list[i],3.)))); */
	for (j = 0; j < ebins; ++j) {
	  spectrum[i][j] = singleTspec[i][j]*js->geomfactors[i][i];
	  
	  ownfrac = spectrum[i][j];
	  
	  for (k = i+1; k < js->nannuli+1; ++k) {
	    
	    norm = (js->geomfactors[i][k]+js->geomfactors[i+1][k-1]- 
		    js->geomfactors[i][k-1]-js->geomfactors[i+1][k]); 
	    /* 	      if (k < js->nannuli && j ==0) */
	    /* 		printf("Proj: %.1E-%.1E; 3D: %.1E-%.1E:  %6.0E\n", */
	    /* 		       js->r1list[i],js->r2list[i], */
	    /* 		       js->r1list[k],js->r2list[k] */
/* 		       ,norm/((1./3.)*(pow(js->r2list[k],3.)- */
/* 				       pow(js->r1list[k],3.)))); */
	    
	    norm *= singleTspec[k][j];
	    
	    spectrum[i][j] += norm;
	    
	  }
	}
	
      }
	//for (i = 0; i < js->nannuli; ++i)  {
	//printf("%d ",i);
	//for (j=i; j < js->nannuli; ++j)
	//printf("%.1E ",foo[i][j]*singleTspec[j][50]/spectrum[i][50]);
	//printf("\n");
	// }
    }
  

#else

    if (!js->dofast) 
      emissivity_profile(js);
 
    gsl_set_error_handler_off();

    if (js->gasmodel == TRIAXIAL && js->calcxray && instr < SZ)
      generate_xmap(js,spectrum);
    else {
      
     #pragma omp parallel private(i,j) firstprivate(jv) shared(js,spectrum,ebins,instr) default(none)
      {


      jaco_thread_init(&jv);
      jv.workspace[0] = gsl_integration_workspace_alloc(INTMAX);
      jv.workspace[1] = gsl_integration_workspace_alloc(INTMAX);
    

      if (js->dosim > 1) {

	//#pragma omp for nowait schedule(dynamic)
        for (i = 0; i < js->nannuli; ++i) 
	  simulate_spec(i,ebins,spectrum[i],&jv);

      } else {

	if (instr < SZ) {
          #pragma omp for  schedule(dynamic) 
	    for (i = 0; i < js->nannuli; ++i) 
	      jaco_xray(i,ebins,spectrum[i],&jv);
	}
      
#if defined(WITH_SZ) | defined(WITH_BOLOCAM)
	if (instr == SZ){
          #pragma omp for schedule(dynamic) 
	  for (i = 0; i < js->nsplines-1; ++i)
	    jaco_sz(i,&jv);
	}

#endif
      
      }

      gsl_integration_workspace_free(jv.workspace[0]);
      gsl_integration_workspace_free(jv.workspace[1]);
      jaco_thread_shutdown(&jv);
    }
  }
    
#endif

     gsl_set_error_handler(NULL);
  
#if defined(WITH_SZ) | defined(WITH_BOLOCAM)
    if (js->calcsz) {
      if ((js->dofast && newmodel) || instr == SZ) {
	generate_ymap(js);
	if (savesz == NULL) savesz = sci_dvector(js->szparams.nest);
	big_get_sz_model(savesz,&js->szparams);

	if (instr == SZ)
	  for (i = 0; i < js->szparams.nest; ++i)
	    outValues[i] = savesz[i];
      }
    }
#endif

    // No model has any business beyond this point other than the X-ray model.
    if (instr >= SZ)
      return 0;


    if (params[EXNORM] > 0.) 
      for (j = 0; j < ebins; ++j) {
	for (k = 0; k < js->nannuli; ++k) {
	  //if (fabs(params[GASMODEL]) > 2) spectrum[0][j] = 0;
	  //spectrum[0][j] += corespectrum[j]*absfac[0][j];
	  if (js->r1list[k] < js->cutradius || 
	      k == 0) 
	    spectrum[k][j] += corespectrum[j]*js->ringarea[k];
	}
      }

    if (js->dosim > 2)
      make_xrayimage(js);
    
  }

  if (js->dosim > 3)
    exit(1);

  for (j = 0; j < *ndata; ++j) 
    finalspectrum[j] = 0;

  for (j = 0; j < ebins; ++j) {
    background[j] = plbackspec[j];
    softbackground[j] = backspec[j];
    js->background[j] = plbackspec[j];
  }

  // Go through all annuli
    double ff;
    ff=0;
  for (i = 0; i < js->annuluscount[instr]; ++i) {

    // There are N blocks of N rows. Go to the top of the ith block,
    // and move down to the row corresponding *this* annulus
    elem = i*js->annuluscount[instr]+js->psfrow[instr][annulus];
    //printf("--%ld--",elem);

    // These coefficients represent that amount of light scattered by
    // the ith annulus into *this* annulus.
    c0 = js->psf[3*instr][elem];
    c1 = js->psf[3*instr+1][elem];
    c2 = js->psf[3*instr+2][elem];

    id = js->annulusid[instr][i];

    //printf("%d %E %E\n",i,spectrum[id][0+js->l1],absfac[id][js->l1]);
    for (j = 0; j < *ndata; ++j) {
      normspec = absfac[id][j+js->l1]*
	spectrum[id][j+js->l1]+(softbackground[j+js->l1]*params[BNORM0+instr])*js->ringarea[id];
	//+background[j+js->l1])**
      //	params[BNORM0+instr];
      //}
      if (!FINITE(absfac[id][j+js->l1])) printf("%d %d %E spectrum error\n",
						  j,js->l1,normspec);

      finalspectrum[j] += normspec*(c0+x0[j]*(c1+x0[j]*c2));
      ff += finalspectrum[j];
      
    }
    //printf("\n");
    //j = 0; if (id == 0)
    //printf("%ld %ld %10.3E %10.3E %10.3E %10.3E\n",id,j,absfac[j+js->l1],background[j+js->l1],js->ringarea[annulus],spectrum[id][j]);
  }
  /* printf("Annulus %2d Total: %.3E\n",annulus,ff/js->ringarea[annulus]); */
      /* 	printf("Instr %d Annulus %d: getting ",instr,annulus); */
      /* 	printf("%E from %02d\n", (c0+x0[j]*(c1+x0[j]*c2)),i); */
      /* } */
  //printf("\n");
  //printf("%20.10E%20.10E",finalspectrum[10],spectrum[id][10]);

  
  // 40 50 47 52
  if (debug > 0) {

    timenow = clock();

    
    runtime += timenow-cl1;
    nclock = (double)runtime/CLOCKS_PER_SEC;

    if (nclock > 5) {
      printf("Speed: %.1f unique runs per second (%ld/%.0f); "\
             "total = %.0fs\n",
	     nrun/nclock,nrun,nclock,(double)(timenow-cl0)/CLOCKS_PER_SEC);

      nrun = 0;
      nclock = 0;
      runtime = 0;
      cl0 = timenow;
    }
  }

  for (j = 0; j < *ndata; ++j) {  
    //finalspectrum[j] +=  background[j+js->l1]*params[PLNORM0+instr];
    finalspectrum[j] /= js->ringarea[annulus]; 
    outValues[j] = finalspectrum[j];
    //printf("%f %f\n",x0[j],outValues[j]);
  } 

  if (debug > 2)
    for (j = 0; j < *ndata; ++j) 
      fprintf(js->logfile,"%f %f\n",x0[j],outValues[j]);
  
  //printf("%d %d %d %f\n",samepars,annulus,instr,ringarea[instr][annulus]);

  if (debug > 0)
    fclose(js->logfile);

  nolddata = *ndata;
  free(finalspectrum);

  return 0.;
}

