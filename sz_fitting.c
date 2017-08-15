#include <stdio.h>
#include <stdlib.h>
#include "sz_fitting.h"
#include <zlib.h>
#include <nrutil.h>
#include <tweaked_recipes.h>
#include <zlib.h>
#include <math.h>
#include <string.h>
#include <assert.h>
#include <complex.h>
#include <fftw3.h>
#define SZ_NEW_GRIDDER

#define USE_RGRID
#define USE_CBLAS

#define HACK_FAC 1.0//-5.31375 //appears to be an extra T_CMB*1.95 that I want to get rid of in the dump format.
//#define DEBUG_SZ

#ifdef USE_CBLAS
#include <clapack.h>
#include <cblas.h>
#endif
/*---------------------------------------------------------------------------*/
/*void cdgemv(char trans, int m, int n, double alpha, double *a, int lda, double *x, int incx, double beta, double *y, int incy)
{
  dgemv_(&trans,&m,&n,&alpha,a,&lda,x,&incx,&beta,y,&incy,1);
}
*/
/*---------------------------------------------------------------------------*/

#if 1

void  cdsyev( char jobz, char uplo, int n, double *a, int lda, double *w, double *work, int lwork, int *info)
{
#ifdef USE_CBLAS
  dsyev_(&jobz,&uplo,&n,a,&lda,w,work,&lwork,info,1,1);
#endif
}

/*---------------------------------------------------------------------------*/



void print_sz_params(SZControlParams *params)
{
  fprintf(params->szof,"\n======================\n");
  fprintf(params->szof,"Printing SZ Parameters\n");
  fprintf(params->szof,"======================\n\n");
  fprintf(params->szof,"mapcent is %12.6f %12.6f\n",params->racent,params->deccent);
  fprintf(params->szof,"griddercent is %12.6f %12.6f\n",params->GridderRACent,params->GridderDecCent);
  fprintf(params->szef,"Expecting %d pixels in map of %8.3f arcminutes.\n",params->npix,60.0*params->pixsize);
  fprintf(params->szef,"datamin/max are %14.6e %14.6e\n",params->datamin,params->datamax);
  fprintf(params->szef,"freq is %14.4e\n",params->freq);
  fprintf(params->szef,"have %d estimators\n",params->nest);
  fprintf(params->szef,"configfile is %s\n",params->configfile);
  fprintf(params->szef,"GridderExec is %s\n",params->GridderExec);
  fprintf(params->szef,"MockCBIExec is %s\n",params->MockCBIExec);
  fprintf(params->szef,"GridderScriptName is %s\n",params->GridderScriptName);
  fprintf(params->szef,"MockCBIScriptName is %s\n",params->MockCBIScriptName);
  fprintf(params->szef,"FitsFile is %s\n",params->FitsFile);
  fprintf(params->szef,"UVFModelName is %s\n",params->UVFModelName);
  fprintf(params->szef,"UvfDataName is %s\n",params->UVFDataName);
  fprintf(params->szef,"GridderDataBase is %s\n",params->GridderDataBase);
  fprintf(params->szef,"EstiamtorName is %s\n",params->EstimatorName);
  fprintf(params->szef,"DerotDataFile is %s\n",params->DerotDataFile);
  fprintf(params->szef,"FlipEndian is %d\n",params->FlipEndian);


  fprintf(params->szef,"UVF File mocks will be based on is %s\n",params->UVFDataName);
  fprintf(params->szef,"SZ Estimator base lives in %s\n",params->GridderDataBase);

  fprintf(params->szef,"MockCBI executable=%s\n",params->MockCBIExec);
  fprintf(params->szef,"Gridder executable=%s\n",params->GridderExec);

  fprintf(params->szef,"Temporary FITS Image=%s\n",params->FitsFile);
  fprintf(params->szef,"Temporary UVF Model=%s\n",params->UVFModelName);
  fprintf(params->szef,"Temporary Gridded Model=%s\n",params->EstimatorName);
}

/*---------------------------------------------------------------------------*/
void parse_sz_default_params(SZControlParams *params, char *line_in)
{
  char **argv;
  int argc,*found_list,was_odd,i,have_tag;

  argv=create_argv(line_in,&argc," \n");
  found_list=ivector(0,argc);
  for (i=1;i<=argc;i++)
    found_list[i]=0;
  was_odd=0;

  fprintf(params->szof,"\n");

  if ((strlen(params->JacoTag)>0)&&(strlen(params->JacoTag)<SZSTRINGLEN-1))
    have_tag=1;
  else
    have_tag=0;

  if (params->do_defaults)
    sprintf(params->MockCBIExec,"/nfs/rab0/home/sievers/mockcbi");
  if (get_command_line_string(argc,argv,"-mockcbi_exec",params->MockCBIExec,found_list)==0)
    fprintf(params->szof,"MockCBI executable is %s\n",params->MockCBIExec);

  if (params->do_defaults)
    sprintf(params->GridderExec,"/nfs/rab0/home/sievers/gridder_mpi/mpigridr53x");
  if (get_command_line_string(argc,argv,"-gridder_exec",params->GridderExec,found_list)==0)
    fprintf(params->szof,"Gridder  executable is %s\n",params->GridderExec);

  if (params->do_defaults)
    sprintf(params->UVFDataName,"/nfs/rab0/home/sievers/szsrc/a478sub_all-sub_uvsub.uvf");
  if (get_command_line_string(argc,argv,"-uvf_data_name",params->UVFDataName,found_list)==0)
    fprintf(params->szof,"UVF Data Name is %s\n",params->UVFDataName);

  if (params->do_defaults)
    sprintf(params->GridderDataBase,"/nfs/rab0/home/sievers/szsrc/a478_shape_1bin");
  if (get_command_line_string(argc,argv,"-gridder_data_base",params->GridderDataBase,found_list)==0)
    fprintf(params->szof,"Gridder data base is %s\n",params->GridderDataBase);
  
  if (params->do_defaults)
    sprintf(params->GridderNoiseName,"%s_NOIS.gz",params->GridderDataBase);
  if (get_command_line_string(argc,argv,"-gridder_noise_name",params->GridderNoiseName,found_list)==0)
    fprintf(params->szof,"Gridder noise name is %s\n",params->GridderNoiseName);

  if (params->do_defaults)
    params->skip_cmb=0;
  if (exists_in_command_line(argc,argv,"-skip_cmb",found_list)==0) {
    params->skip_cmb=1;
    printf("Going to just use the noise matrix for the noise.\n");
  }
  if (params->do_defaults)
    params->freq=30.0;
  if (get_command_line_double(argc,argv,"-freq",&(params->freq),found_list)==0)
    printf("setting frequency to %8.3f GHz.\n",params->freq);
  

  if (params->do_defaults)
    sprintf(params->GridderScriptName,"./run_cbigridder_sz.scr");
  if (get_command_line_string(argc,argv,"-gridder_script_name",params->GridderScriptName,found_list)==0)
    fprintf(params->szof,"Gridder script base name is %s\n",params->GridderScriptName);
#ifdef SZ_NEW_GRIDDER
  fprintf(stderr,"Hello from new gridder!\n");
  sprintf(params->GridderScriptName,"%s",params->GridderScriptName);
#else
  if (have_tag)
    append_tag(params->GridderScriptName,params->JacoTag,SZSTRINGLEN-1);
  sprintf(params->GridderScriptName,"%s.%03d",params->GridderScriptName,params->node);
#endif



  if (params->do_defaults)
    sprintf(params->MockCBIScriptName,"./run_cbimockcbi_sz.scr");
  if (get_command_line_string(argc,argv,"-mockcbi_script_name",params->MockCBIScriptName,found_list)==0)
    fprintf(params->szof,"Mockcbi script base name is %s\n",params->MockCBIScriptName);
  if (have_tag)
    append_tag(params->MockCBIScriptName,params->JacoTag,SZSTRINGLEN-1);
  sprintf(params->MockCBIScriptName,"%s.%03d",params->MockCBIScriptName,params->node);

  if (params->do_defaults)
    sprintf(params->EstimatorName,"./cbigridder_forestimators");
  if (get_command_line_string(argc,argv,"-estimator_name",params->EstimatorName,found_list)==0)
    fprintf(params->szof,"Estimator name is %s.%03d\n",params->EstimatorName,params->node);
  if (have_tag)
    append_tag(params->EstimatorName,params->JacoTag,SZSTRINGLEN-1);
  sprintf(params->EstimatorName,"%s.%03d",params->EstimatorName,params->node);

  if (params->do_defaults)
    sprintf(params->FitsFile,"./temporary_szfile.fits");
  if (get_command_line_string(argc,argv,"-fits_file",params->FitsFile,found_list)==0)
    fprintf(params->szof,"Fits file name is %s\n",params->FitsFile);
  if (have_tag)
    append_tag(params->FitsFile,params->JacoTag,SZSTRINGLEN-1);
  sprintf(params->FitsFile,"%s.%03d",params->FitsFile,params->node);

  if (params->do_defaults)
    sprintf(params->UVFModelName,"./temporary_szfile.uvf");
  if (get_command_line_string(argc,argv,"-uvf_model_name",params->UVFModelName,found_list)==0)
    fprintf(params->szof,"UVF model name is %s\n",params->UVFModelName);
  if (have_tag)
    append_tag(params->UVFModelName,params->JacoTag,SZSTRINGLEN-1);
  sprintf(params->UVFModelName,"%s.%03d",params->UVFModelName,params->node);

  if (params->do_defaults)
    sprintf(params->configfile,"./cbiconfig.txt");
  if (get_command_line_string(argc,argv,"-configfile",params->configfile,found_list)==0)
    fprintf(params->szof,"Configfile name is %s\n",params->configfile);

  if (params->do_defaults)
    sprintf(params->DerotDataFile,"./cbi_szdata.dat");
  if (get_command_line_string(argc,argv,"-derot_data_file",params->DerotDataFile,found_list)==0)
    fprintf(params->szof,"Derot data file name is %s\n",params->DerotDataFile);
  if (have_tag)
    append_tag(params->DerotDataFile,params->JacoTag,SZSTRINGLEN-1);
  sprintf(params->DerotDataFile,"%s.%03d",params->DerotDataFile,params->node);

  if (params->do_defaults)
    {
      sprintf(params->DifferenceDataFile,"./data_diff_SVEC");
      params->WriteDifference=0;
    }
  if (get_command_line_string(argc,argv,"-data_diff_file",params->DifferenceDataFile,found_list)==0)
    {
      fprintf(params->szof,"Difference data file name is %s\n",params->DifferenceDataFile);
      params->WriteDifference=1;
      if (have_tag)
	append_tag(params->DifferenceDataFile,params->JacoTag,SZSTRINGLEN-1);
      sprintf(params->DifferenceDataFile,"%s.%03d",params->DifferenceDataFile,
	      params->node);
    }




  if (params->do_defaults)
    params->FlipEndian=1;
  if (get_command_line_int(argc,argv,"-flip_endian",&params->FlipEndian,found_list)==0)
    {
      if (params->FlipEndian)
	fprintf(params->szof,"Going to flip endian\n");
      else
	fprintf(params->szof,"Not going to flip endian.\n");
    }

  if (params->do_defaults)
    params->npix=512;
  if (get_command_line_int(argc,argv,"-npix",&params->npix,found_list)==0)
    fprintf(params->szof,"Map will have %d pixels\n",params->npix);

  if (params->do_defaults)
    params->pixsize=0.5/60;
  /*do pixel size in arcminutes*/
  if (get_command_line_double(argc,argv,"-pixsize",&params->pixsize,found_list)==0)
    {
      fprintf(params->szof,"Map will have %8.3f arcminute pixels\n",params->pixsize);
      params->pixsize=params->pixsize/60.0;
    }

  if (params->do_defaults)
    params->calib_fac=1.0;
  if (get_command_line_double(argc,argv,"-calib_factor",&params->calib_fac,found_list)==0)
    fprintf(params->szof,"SZ data will be scaled by a factor of %6.3f.\n",params->calib_fac);


  if (params->do_defaults)
    params->cmb_scale_fac=1.0;
  if (get_command_line_double(argc,argv,"-cmb_scale_factor",&params->cmb_scale_fac,found_list)==0)
    fprintf(params->szof,"CMB covariance matrix will be scaled by a factor of %6.3f.\n",params->cmb_scale_fac);



  if (params->do_defaults)
    {
      /*don't default overwrite if we already have something here*/
      if (params->racent==0)
	params->racent=63.33525; /*why not...  If we didn't specify anything, probably wanted A478 anyways*/
    }
  if (get_command_line_double(argc,argv,"-racent",&params->racent,found_list)==0)
    fprintf(params->szef,"Setting map RA\n");

  if (params->do_defaults)
    {
      /*don't default overwrite if we already have something here*/
      if (params->deccent==0)
	params->deccent=10.465; /*why not...  If we didn't specify anything, probably wanted A478 anyways*/
    }
  if (get_command_line_double(argc,argv,"-deccent",&params->deccent,found_list)==0)
    fprintf(params->szef,"Setting map declination\n");
  
  fprintf(params->szof,"Map center is at RA=%0.6f, DEC=%0.6f\n",params->racent,params->deccent);


  if (params->do_defaults)
    {
      /*don't default overwrite if we already have something here*/
      if (params->GridderRACent==0)
	params->GridderRACent=params->racent; /*default to original map setting*/
    }
  if (get_command_line_double(argc,argv,"-gridder_racent",&params->GridderRACent,found_list)==0)
    fprintf(params->szef,"Setting gridder RA\n");

  if (params->do_defaults)
    {
      /*don't default overwrite if we already have something here*/
      if (params->GridderDecCent==0)
	params->GridderDecCent=params->deccent; /*default to original map setting*/
    }
  if (get_command_line_double(argc,argv,"-gridder_deccent",&params->GridderDecCent,found_list)==0)
    fprintf(params->szef,"Setting gridder declination\n");

  fprintf(params->szof,"Gridder center is at RA=%0.6f, DEC=%0.6f\n",params->GridderRACent,params->GridderDecCent);


  if (get_command_line_string(argc,argv,"-rdump_name",params->dump_name_rvec,found_list)==0)
    fprintf(params->szof,"Rvec dumpfile is %s\n",params->dump_name_rvec);

  if (get_command_line_string(argc,argv,"-sdump_name",params->dump_name_svec,found_list)==0)
    fprintf(params->szof,"Svec dumpfile is %s\n",params->dump_name_svec);

  sprintf(params->svec_name,"");
  if (get_command_line_string(argc,argv,"-svec_name",params->svec_name,found_list)==0)
    fprintf(params->szof,"Svec datafile is %s\n",params->svec_name);
  

  params->oversamp=0;
  if (get_command_line_int(argc,argv,"-oversamp",&params->oversamp,found_list)==0)
    fprintf(params->szof,"Map will be oversampled by a factor of %d.\n",params->oversamp);
  
  


  /*print out the stuff we didn't use that it looked like we wanted to.*/
  fprintf(params->szef,"\n");
  for (i=1;i<=argc;i++)
    if (found_list[i]==0)
      fprintf(stderr,"Unkown/unused command %s\n",argv[i]);

  
}
/*---------------------------------------------------------------------------*/
int initialize_sz(SZControlParams *params)
{
  FILE *outfile,*infile;
  char to_exec[SZSTRINGLEN],filename[SZSTRINGLEN],*line_in;
  double **noise, **cmb,**cov,*eigvals,*workvec,worksize,*data_derot;
  int lwork,info,i;

  if (params->node == 0 && !params->quiet) {
    params->szof = stdout;
    params->szef = stderr;
  } else {
    params->szof = fopen("/dev/null","w+");
    params->szef = fopen("/dev/null","w+");
  }

#if 1
  if ((strlen(params->ParamFile)>0)&&(strlen(params->ParamFile)<SZSTRINGLEN-1))
    line_in=read_all_stdin(params->ParamFile,params);
  else
    line_in=read_all_stdin("szfitting_defaults.conf",params);
  if (line_in == NULL) return -1;

  params->do_defaults=1;
  parse_sz_default_params(params, line_in);

#else

  params->FlipEndian=1;
  sprintf(params->configfile,"cbiconfig.txt");
  sprintf(params->DerotDataFile,"cbi_szdata.dat");

  params->npix=512;
  params->pixsize=0.5/60; /*pixelsize = .5 arcmin*/
  params->npix=256;
  params->pixsize=0.5/60; /*pixelsize = .5 arcmin*/

  /*a478 4:13:25.26, 10:27:54.2*/

  params->racent=63.33525;  /*going to be set externally*/
  params->deccent=10.465;  /*going to be set externally*/
  params->GridderRACent=params->racent;
  params->GridderDecCent=params->deccent;


  sprintf(params->GridderScriptName,"./run_cbigridder_sz.scr");
  sprintf(params->MockCBIScriptName,"./run_mockcbi_sz.scr");
  sprintf(params->EstimatorName,"cbigridder_forestimators");
  sprintf(params->FitsFile,"temporary_szfile.fits");
  sprintf(params->UVFModelName,"temporary_szfile.uvf");


  /*sprintf(params->GridderExec,"/cita/d/raid-sievers2/sievers/clusters/mpigridr53x_beamsky");
  sprintf(params->MockCBIExec,"/cita/d/raid-sievers2/sievers/clusters/mockcbi");
  sprintf(params->UVFDataName,"/cita/d/raid-sievers/sievers/clusters/oldcbi/a478/a478sub_all-sub_uvsub.uvf");
  sprintf(params->GridderDataBase,"/cita/d/raid-sievers2/sievers/clusters/cres/a478_shape_1bin");*/

  sprintf(params->GridderExec,"/nfs/rab0/home/sievers/gridder_mpi/mpigridr53x");
  sprintf(params->MockCBIExec,"/nfs/rab0/home/sievers/mockcbi");
  sprintf(params->UVFDataName,"/nfs/rab0/home/sievers/szsrc/a478sub_all-sub_uvsub.uvf");
  sprintf(params->GridderDataBase,"/nfs/rab0/home/sievers/szsrc/a478_shape_1bin");
  
  
#endif


#ifndef USE_RGRID

#ifdef SZ_NEW_GRIDDER
  outfile=fopen(params->GridderScriptName,"r");  
#else
  outfile=fopen(params->GridderScriptName,"w");
#endif
  if (!outfile)
    {
      fprintf(params->szef,"Error - unable to open gridder script %s in initialize_sz.\n",params->GridderScriptName);
      exit(EXIT_FAILURE);
    }
  fclose(outfile);
#ifndef SZ_NEW_GRIDDER
  sprintf(to_exec,"chmod a+x %s",params->GridderScriptName);
  system(to_exec);
#endif

  outfile=fopen(params->MockCBIScriptName,"w");
  if (!outfile)
    {
      fprintf(params->szef,"Error - unable to write to mockcbi script %s in initialize_sz.\n",params->MockCBIScriptName);
      exit(EXIT_FAILURE);
    }
  fclose(outfile);
  //#ifndef SZ_NEW_GRIDDER
  sprintf(to_exec,"chmod a+x %s",params->MockCBIScriptName);
  system(to_exec);
  //#endif
  
#endif
#if 0  
  if (!infile)
    return -1;
  return 0;
#endif

  if (strlen(params->svec_name)==0)
    sprintf(filename,"%s_SVEC",params->GridderDataBase);
  else
    strncpy(filename,params->svec_name,SZSTRINGLEN-1);
  params->nest=GetNEstimator(filename);
  fprintf(params->szof,"Expect %d estimators\n",params->nest);
  printf("Expect %d estimators\n",params->nest);
  params->data=dvector(0,params->nest-1);


  ReadSvec(params->data,filename);
  /*read in the original data and save it so we can have a
    gander at difference maps*/
  params->data_org=dvector(0,params->nest-1);
  ReadSvec(params->data_org,filename);

  //start change
  //Apply a calibration factor (if needed) to the data
  for (i=0;i<params->nest;i++)
    {
      params->data[i]*=params->calib_fac;
      params->data_org[i]*=params->calib_fac;
    }
  //stop change


#ifdef USE_RGRID
  params->dump=read_tt_dump_svec(params->dump_name_svec);
  read_tt_dump_rvec(params->dump,params->dump_name_rvec);
  clear_dump_zeros(params->dump);
  params->svec_data=ReadSvecStruct(filename);
  params->dump->npix_use=get_npix_from_oversamp(params->dump->gridsize,params->oversamp);
  params->dump->pixel_size/=(double)(1+2*params->oversamp);
  params->npix=params->dump->npix_use;
  params->pixsize=params->dump->pixel_size;
  fprintf(stderr,"Setting npix and pixel size to %d %14.6g\n",params->npix,params->pixsize);
  //make sure there's a position in the data file.
  assert(!((params->svec_data->ra==0)&&(params->svec_data->dec==0)));
  params->racent=params->svec_data->ra*180/M_PI;
  params->deccent=params->svec_data->dec*180/M_PI;

#endif


  //sprintf(filename,"%s_NOIS.gz",params->GridderDataBase);
  noise=dpmatrix(0,params->nest-1,0,params->nest-1);
  if (ReadNoise(noise,params->GridderNoiseName,params->nest)!=params->nest)
    {
      fprintf(params->szef,"Error reading Noise matrix from file %s\n",params->GridderNoiseName);
      return -1;
    }  


  cov=dpmatrix(0,params->nest-1,0,params->nest-1);
  //start change
  //multiply the noise by calib_factor^2 (since we rescaled the data, noise needs to be rescaled as well)
  //  also, if the data need a CMB scaling factor, apply that as well.
  memset(cov[0],0,sizeof(double)*params->nest*params->nest);
  /*cblas_dcopy(params->nest*params->nest,noise[0],1,cov[0],1);*/
#ifdef USE_CBLAS
  cblas_daxpy(params->nest*params->nest,params->calib_fac*params->calib_fac,noise[0],1,cov[0],1);
  //cblas_daxpy(params->nest*params->nest,params->cmb_scale_fac,cmb[0],1,cov[0],1);
#endif

  if (!params->skip_cmb) {
    cmb=dpmatrix(0,params->nest-1,0,params->nest-1);
    sprintf(filename,"%s_TT.gz",params->GridderDataBase);
    if (ReadCMBMat(cmb,filename,params->nest)!=params->nest)
      {
	fprintf(params->szef,"Error reading CMB signal matrix from file %s\n",filename);
	return -1;
      }  
#ifdef USE_CBLAS
    cblas_daxpy(params->nest*params->nest,params->cmb_scale_fac,cmb[0],1,cov[0],1);
#endif
    free(cmb[0]);
    free(cmb);

  }

  //stop change
  

  /*
  fprintf(params->szef,"noise[0][0]=%14.6g\n",noise[0][0]);
  fprintf(params->szef,"noise[19][12]=%14.6g\n",noise[19][12]);
  fprintf(params->szef,"noise[12][19]=%14.6g\n",noise[12][19]);

  fprintf(params->szef,"cmb[0][0]=%14.6g\n",cmb[0][0]);
  fprintf(params->szef,"cmb[19][12]=%14.6g\n",cmb[19][12]);
  fprintf(params->szef,"cmb[12][19]=%14.6g\n",cmb[12][19]);


  fprintf(params->szef,"cov[0][0]=%14.6g\n",cov[0][0]);
  fprintf(params->szef,"cov[19][12]=%14.6g\n",cov[19][12]);
  fprintf(params->szef,"cov[12][19]=%14.6g\n",cov[12][19]);
  */
  eigvals=dvector(0,params->nest-1);

  outfile=fopen("covmat_out.bin","w");
  fwrite(&params->nest,sizeof(int),1,outfile);
  fwrite(cov[0],sizeof(double),params->nest*params->nest,outfile);
  fclose(outfile);
  
  outfile=fopen("data_out.bin","w");
  fwrite(&params->nest,sizeof(int),1,outfile);
  fwrite(params->data,sizeof(double),params->nest,outfile);
  fclose(outfile);

  


  cdsyev('V','U',params->nest,cov[0],params->nest,eigvals,&worksize,-1,&info);
  lwork=worksize+1;
  workvec=dvector(0,lwork);

  fprintf(params->szef,"Work matrix size is %10.0f\n",worksize);  
  
  cdsyev('V','U',params->nest,cov[0],params->nest,eigvals,workvec,lwork,&info);

  fprintf(params->szef,"Back from cdsyev\n");
  data_derot=dvector(0,params->nest-1);
#ifdef USE_CBLAS
  cblas_dgemv(CblasColMajor,CblasTrans,params->nest,params->nest,1.0,cov[0],params->nest,params->data,1,0.0,data_derot,1);
#endif
  
  outfile=fopen(params->DerotDataFile,"w");
  for (i=0;i<params->nest;i++)
    fprintf(outfile,"%4d %14.6e %14.6e\n",i,data_derot[i],sqrt(eigvals[i]));
  fclose(outfile);

  for (i=0;i<params->nest;i++)
    params->data[i]=data_derot[i];
  params->eigvecs=cov;

  
  params->work=dvector(0,params->nest-1); /*workspace for the rotation down the road*/
  params->model=dvector(0,params->nest-1);
  params->eigvals=eigvals;
#ifdef USE_RGRID
  params->ymap=dpmatrix(0,params->npix-1,0,params->npix-1);
#else
  params->ymap=dpmatrix(0,params->nest-1,0,params->nest-1);
#endif
  /*You would change the -1.95 factor if you had SZ maps at different frequencies*/
  params->y_to_t=-1.95*2.725;
  

  /*
  free(data_derot);
  fprintf(params->szef,"Freed data_derot.\n");
  free(&noise[0][0]);
  fprintf(params->szef,"Freed big noise.\n");
  free(noise);
  fprintf(params->szef,"Freed small noise.\n");
  free(&cmb[0][0]);
  fprintf(params->szef,"Freed big cmb.\n");
  free(cmb);
  fprintf(params->szef,"Freed small cmb.\n");
  free(workvec);
  fprintf(params->szef,"Freed workvec.\n");
  */
  params->mapcol=vector(0,params->npix-1);
  fprintf(params->szef,"Finished SZ initialization.\n");
  return 0;
}
/*---------------------------------------------------------------------------*/
void flip_float_endian(float *x)
{
  char *cc,hold;

  cc=(char *)x;
  hold=cc[3];
  cc[3]=cc[0];
  cc[0]=hold;
  hold=cc[2];
  cc[2]=cc[1];
  cc[1]=hold;
}
/*--------------------------------------------------------------------------------*/
int get_npix_from_oversamp(int npix_in, int oversamp)
{
  //npix_in--;
  //return 1+npix_in*(1+2*oversamp);
  
  return (npix_in-1)*(1+2*oversamp)+1;
}

/*---------------------------------------------------------------------------*/
double yfac(SZControlParams *data)
{
  //at some point, put in SZ spectrum correction.  For now, assume R-J
  double xnu=data->freq/20.83674/2.725;
  double val=(xnu*(exp(xnu)+1.0)/(exp(xnu)-1.0))-4.0;
  //return -1/2.725;
  //return 0;

  return val;
  
}
/*---------------------------------------------------------------------------*/
//int get_sz_model_dump(double *model, double **map, DumpData *data)

int get_sz_model_dump(double *model, SZControlParams *params)
{

  double **map=params->ymap;
  DumpData *data=params->dump;
  

  if (0)
    {
      FILE *outfile=fopen("sz_map_dump.in","w");
      int nnn=data->npix_use;
      fwrite(&nnn,1,sizeof(int),outfile);
      fwrite(map[0],nnn*nnn,sizeof(double),outfile);
      fclose(outfile);
    }
  dflip_mat(map,data->npix_use);
  transpose_matrix(map,data->npix_use);
  memset(params->model,0,params->nest*sizeof(double));
  //printf("mapsize and pixel size are %d %14.4e\n",data->npix_use,data->pixel_size);
  complex double **mapft=complex_dmatrix(data->npix_use,data->npix_use);
  int icent=(data->npix_use-1)/2;
  int nb2=(data->nbeam-1)/2;
  fft_2d_map_wshifts(map,mapft,data->npix_use,data->npix_use);
  //printf("corner element is %14.4e %14.4e\n",creal(mapft[0][0]),cimag(mapft[0][0]));
  double normfac=data->npix_use*data->npix_use;
  double yyfac=yfac(params);
  //fprintf(stderr,"Data->n_nonzero is %d\n",data->n_nonzero);
#ifdef DEBUG_SZ
  int i1min=1000;
  int i1max=0;
  int i2min=1000;
  int i2max=0;
#endif
  //#pragma omp parallel for shared(data,params,yyfac,normfac,mapft,icent,nb2) default(none)
  for (int ii=0;ii<data->n_nonzero;ii++) {
    //this appears to be correct, but very slow.
    //please remove this comment when fixed!
    complex double val=0;
    for (int i=-nb2;i<=nb2;i++)
      for (int j=-nb2;j<=nb2;j++) {
#ifdef DEBUG_SZ
	int ind1=i+icent+data->vvec[ii];
	int ind2=j+icent+data->uvec[ii];
	if (ind1<i1min)
	  i1min=ind1;
	if (ind1>i1max)
	  i1max=ind1;
	if (ind2<i2min)
	  i2min=ind2;
	if (ind2>i2max)
	  i2max=ind2;
#endif
	val+=data->rgrid[ii][i+nb2][j+nb2]*mapft[i+icent+data->vvec[ii]][j+icent+data->uvec[ii]]; //original model
	//val+=data->rgrid[ii][i+nb2][j+nb2]*mapft[j+icent+data->uvec[ii]][i+icent+data->vvec[ii]];

      }
    params->model[2*ii]=creal(val)/data->zvec[ii]*yyfac/normfac/HACK_FAC;
    params->model[2*ii+1]=cimag(val)/data->zvec[ii]*yyfac/normfac/HACK_FAC;
  }
#ifdef DEBUG_SZ
  printf("i1min/max are %d/%d, and i2min/max are %d/%d, with nb %d\n",i1min,i1max,i2min,i2max,nb2);
#endif
  free(mapft[0]);
  free(mapft);
  if (0) {
    FILE *outfile=fopen("sz_estimators_raw.out","w");
    fwrite(&(params->nest),sizeof(int),1,outfile);
    fwrite(params->model,params->nest,sizeof(double),outfile);
    fwrite(data->uvec,params->nest/2,sizeof(int),outfile);
    fwrite(data->vvec,params->nest/2,sizeof(int),outfile);
    fclose(outfile);
  }
  if (0) {
    
    FILE *outfile=fopen("sz_rgrid_used.out","w");
    fwrite(&(params->nest),sizeof(int),1,outfile);
    fwrite(&(data->nbeam),sizeof(int),1,outfile);
    for (int i=0;i<params->nest/2;i++)
      fwrite(data->rgrid[i][0],data->nbeam*data->nbeam,sizeof(complex double),outfile);
    fclose(outfile);
    
  }
  
  
#ifdef USE_CBLAS
  cblas_dgemv(CblasColMajor,CblasTrans,params->nest,params->nest,1.0,params->eigvecs[0],params->nest,params->model,1,0.0,params->work,1);
  memcpy(params->model,params->work,params->nest*sizeof(double));
#endif
  if (0)
    {
      FILE *outfile=fopen("sz_estimators_rotated.out","w");
      fwrite(&(params->nest),sizeof(int),1,outfile);
      fwrite(params->model,params->nest,sizeof(double),outfile);
      fwrite(data->uvec,params->nest/2,sizeof(int),outfile);
      fwrite(data->vvec,params->nest/2,sizeof(int),outfile);
      fclose(outfile);
    }

  memcpy(model,params->model,sizeof(double)*params->nest);
  return 0;
}
/*---------------------------------------------------------------------------*/
int get_sz_model(double *model, SZControlParams *params)
{
  FILE *FitsFile,*MockScript,*GridderScript;
  float *mapcol;
  double datamin,datamax,chisq,dx,**map;
  int i,j,nbytes,rah,ram,dd,dm;
  double ras,ds;
  char c,to_exec[SZSTRINGLEN],filename[SZSTRINGLEN],sign;
  
  map=params->ymap;
  /*mapcol=vector(0,params->npix-1);*/
  mapcol=params->mapcol;
  
  /*first set the min/max data value for the FITS header*/
  datamin=map[0][0];
  datamax=map[0][0];
  for (i=0;i<params->npix;i++)
    for (j=0;j<params->npix;j++)
      {
	if (map[i][j]<datamin)
	  datamin=map[i][j];
	if (map[i][j]>datamax)
	  datamax=map[i][j];
      }
  params->datamin=datamin;
  params->datamax=datamax;
  
  FitsFile=fopen(params->FitsFile,"w");
  if (!FitsFile)
    {
      fprintf(params->szef,"Error - unable to open file %s for writing FITS image in get_sz_model.\n",params->FitsFile);
      exit(EXIT_FAILURE);
    }
  
  write_simple_fits_header_params(FitsFile,params);
  
  for (i=0;i<params->npix;i++)
    {
      
      /*for (j=0;j<params->npix;j++)
	mapcol[j]=map[i][j];*/
      for (j=0;j<params->npix;j++)
	mapcol[j]=map[j][i];
      if (params->FlipEndian)
	for (j=0;j<params->npix;j++)
	  flip_float_endian(mapcol+j);
      fwrite(mapcol,sizeof(float),params->npix,FitsFile);
    }
  nbytes=sizeof(float)*params->npix*params->npix;
  while (nbytes>2880)
    nbytes -= 2880;
  nbytes=2880-nbytes;
  c=' ';
  for (i=0;i<nbytes;i++)
    fwrite(&c,sizeof(char),1,FitsFile);
  fclose(FitsFile);
  /*fprintf(params->szef,"Wrote fits file.\n");*/
  
  
  MockScript=fopen(params->MockCBIScriptName,"w");
  if (!MockScript)
    {
      fprintf(params->szef,"Error - unable to open %s for writing MockCBIScript in get_sz_model.\n",params->MockCBIScriptName);
      exit(EXIT_FAILURE);
    }
  
  fprintf(MockScript,"#!/bin/csh\n");
  fprintf(MockScript,"%s << EOF >&/dev/null\n",params->MockCBIExec);
  fprintf(MockScript,"cmb_image %s\n",params->FitsFile);
  fprintf(MockScript,"beam CBI\n");
  fprintf(MockScript,"read %s\n",params->UVFDataName);
  fprintf(MockScript,"clear\n");
  fprintf(MockScript,"add_cmb +1\n");
  fprintf(MockScript,"write !%s\n",params->UVFModelName);
  fprintf(MockScript,"exit\n");
  fprintf(MockScript,"EOF\n");
  fclose(MockScript);
  
  system(params->MockCBIScriptName);
#ifdef SZ_NEW_GRIDDER
  {
    char exec_new_gridder[SZSTRINGLEN];
    sprintf(exec_new_gridder,"%s %s %s ",params->GridderScriptName,params->UVFModelName,params->EstimatorName);
    //fprintf(stderr,"exec is:  %s\n",exec_new_gridder);
    //fprintf(stderr,"Greetings from new gridder.\n");
    system(exec_new_gridder);
    
  }
#else
  fprintf(stderr,"Greetings from old gridder.\n");
  GridderScript=fopen(params->GridderScriptName,"w");
  if (!GridderScript)
    {
      fprintf(params->szef,"Error - unable to open gridder script %s for writing in get_sz_model.\n",params->GridderScriptName);
      exit(EXIT_FAILURE);
    }
  fprintf(GridderScript,"#!/bin/csh\n");
  fprintf(GridderScript,"%s << EOF >&/dev/null\n",params->GridderExec);
  fprintf(GridderScript,"-1 1 1 1 0 1 0 0 1 1 0\n\n 0 2700\n3 3\n\n\n\n1 3500\n");
  /**/
  dd2dms(params->GridderRACent/15,&rah,&ram,&ras,&sign);
  dd2dms(params->GridderDecCent,&dd,&dm,&ds,&sign);
  /*fprintf(GridderScript,"4:13:26.2 10:27:57.6\n");*/ /*need to update this at some point*/
  fprintf(GridderScript,"%02d:%02d:%05.2f %c%02d:%02d:%04.1f\n",rah,ram,ras,sign,dd,dm,ds);
  fprintf(GridderScript,"%s\n1.0 1.0 1.0\n\n\n\n -1\n\n%s\n\nEOF\n",params->UVFModelName,params->EstimatorName);
  fclose(GridderScript);
  
  
  
  system(params->GridderScriptName);
#endif

  /*free(mapcol);*/


  sprintf(filename,"%s_SVEC",params->EstimatorName);
  /*fprintf(params->szef,"zeroing workspace.\n");*/
  for (i=0;i<params->nest;i++)
    params->work[i]=0;

  ReadSvec(params->work,filename);
  for (i=0;i<params->nest;i++)
    params->work[i] *= params->y_to_t;
  printf("y_to_t is %14.4e\n",params->y_to_t);
  if (params->WriteDifference)
    {
      if (WriteDiff(params,params->work))
	fprintf(stderr,"Unable to write difference file to %s\n",params->DifferenceDataFile);
    }

#ifdef USE_CBLAS
  cblas_dgemv(CblasColMajor,CblasTrans,params->nest,params->nest,1.0,params->eigvecs[0],params->nest,params->work,1,0.0,params->model,1);
#endif

  /*fprintf(params->szef,"derotated.\n");*/
  chisq=0;
  for (i=0;i<params->nest;i++)
    {
      dx=params->data[i]-params->model[i];
      chisq += dx*dx/params->eigvals[i];
    }
  /*fprintf(params->szef,"chisq is %14.4g\n",chisq);*/
  for (i=0;i<params->nest;i++)
    model[i]=params->model[i];

  return 0;
}

/*---------------------------------------------------------------------------*/

void dd2dms(double dec, int *dd, int *dm, double *ds, char *sign)
{
  double temp;
  int doflip;
  temp=dec;

  if (temp<0)
    {
      doflip=1;
      temp*= -1;
    }
  else
    doflip=0;
  /*fprintf(params->szef,"Temp is %12.7f\n",temp);*/
  *dd=temp;
  temp=temp-(*dd);
  temp *= 60.0;

  *dm=temp;
  temp=temp-(*dm);
  temp *= 60.0;

  *ds=temp;

  if (doflip)
    *sign='-';
  else
    *sign='+'; 
  
}

/*--------------------------------------------------------------------------------*/

int write_simple_fits_header_params(FILE *outfile, SZControlParams *szparams)
     /*FILE *outfile,float racent, float deccent,int npix,float pixsize,float datamin, float datamax)*/
{
  float racent, deccent,pixsize,datamin,datamax;
  int npix;
  float pixcent;
  npix=szparams->npix;
  racent=szparams->racent;
  deccent=szparams->deccent;
  pixsize=szparams->pixsize;
  datamin=szparams->datamin;
  datamax=szparams->datamax;

  //pixcent=((float)npix/2)+1;
  pixcent=(((float)npix+1.0)/2.0);
  fprintf(outfile,"SIMPLE  =                     T                                                 BITPIX  =                   -32                                                 NAXIS   =                     4                                                 NAXIS1  =                 %5d                                                 NAXIS2  =                 %5d                                                 NAXIS3  =                     1                                                 NAXIS4  =                     1                                                 EXTEND  =                     T                                                 COMMENT =   FITS (Flexible Image Transport System) format is defined in 'AstronoCOMMENT =   and Astrophysics', volume 376, page 359; bibcode: 2001A&A...376..359OBJECT  =  'CMB     '                                                           TELESCOP=  'Simulation'                                                         INSTRUME=  'Simulation'                                                         OBSERVER=  'sievers '                                                           DATE-OBS=  '2006-03-07'                                                         BUNIT   =  'K       '                                                           RADECSYS=  'FK5     '                                                           EQUINOX =               2000.00                                                 EPOCH   =               2000.00                                                 ORIGIN  =  'CBISKY (Caltech)'                                                   DATE    =  '2006-03-07T21:09:50'                                                DATAMAX =         %13.5e                                                 DATAMIN =         %13.5e                                                 CTYPE1  =  'RA---SIN'                                                           CRVAL1  =          %12.7f                                                 CRPIX1  =            %10.4f                                                 CDELT1  =             -%08.5f                                                 CROTA1  =              0.000000                                                 CTYPE2  =  'DEC--SIN'                                                           CRVAL2  =            %10.7f                                                 CRPIX2  =            %10.4f                                                 CDELT2  =              %8.5f                                                 CROTA2  =              0.000000                                                 CTYPE3  =  'FREQ    '                                                           CRVAL3  =          3.100000E+10                                                 CRPIX3  =              1.000000                                                 CDELT3  =              0.000000                                                 CROTA3  =              0.000000                                                 CTYPE4  =  'STOKES  '                                                           CRVAL4  =              1.000000                                                 CRPIX4  =              1.000000                                                 CDELT4  =              1.000000                                                 CROTA4  =              0.000000                                                 HISTORY = Spectrum file: wmap.cmb                                               HISTORY = Column:            1                                                  HISTORY = Seed:       -10000                                                    HISTORY = Program:  CBISKY                                                      HISTORY = Version:  4.0 - 2001 Jan 29                                           HISTORY = Run by: sievers on bob (Linux 2.6.12.3-cita)                          HISTORY = Date:  7-Mar-2006 16:09:50                                            END                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                             ",npix,npix,datamax,datamin,racent,pixcent,pixsize,deccent,pixcent,pixsize);
}

#if 0
/*===========================================================================*/
int write_fits_map(float **map, int npix, float pixsize,float racent, float deccent,char *filename)
{
  FILE *outfile;
  int i,j,nbytes;
  float datamin,datamax;
  char c;
  outfile=fopen(filename,"w");
  if (!outfile)
    {
      fprintf(stderr,"Error - unable to write fits file %s\n",filename);
      return -1;
    }
  datamin=map[0][0];
  datamax=map[0][0];
  for (i=0;i<npix;i++)
    for (j=0;j<npix;j++)
      {
	if (map[i][j]<datamin)
	  datamin=map[i][j];
	if (map[i][j]>datamax)
	  datamax=map[i][j];
	flip_float_endian(&map[i][j]);
      }
  write_simple_fits_header_params(outfile,racent,deccent,npix,pixsize,datamin,datamax);
  fwrite(map[0],sizeof(float),npix*npix,outfile);
  nbytes=sizeof(float)*npix*npix;
  while (nbytes>2880)
    nbytes -= 2880;
  nbytes=2880-nbytes;
  c=' ';
  for (i=0;i<nbytes;i++)
    fwrite(&c,sizeof(char),1,outfile);
  fclose(outfile);
  for (i=0;i<npix;i++)
    for (j=0;j<npix;j++)
      flip_float_endian(&map[i][j]);
  
  return 0;
}
#endif
/*---------------------------------------------------------------------------*/
int GetNEstimator(char *filename)
{
  FILE *infile;
  int nbytes_in,nest;
  infile=fopen(filename,"r");
  if (!infile)
    {
      fprintf(stderr,"Unable to open file %s for checking size in GetNEstimator.\n",filename);
      return 0;
    }
  fread(&nbytes_in,sizeof(int),1,infile);
  if (nbytes_in!=12)
    {
      fprintf(stderr,"Error - nybtes in is messed up inside of GetNEstimator on %s.  Got %d\n",filename,nbytes_in);
      return 0;
    }
  fread(&nest,sizeof(int),1,infile);
  fclose(infile);
  return nest;
}


/*---------------------------------------------------------------------------*/

/*--------------------------------------------------------------------------------*/
double *read_fortran_dvec(FILE *infile, int *n)
{
  int nbytes;
  fread(&nbytes,1,sizeof(int),infile);
  *n=nbytes/sizeof(double);
  double *vec=dvector(0,*n-1);
  fread(vec,*n,sizeof(double),infile);
  int nb2=nbytes;
  fread(&nbytes,1,sizeof(int),infile);
  assert(nb2==nbytes);
  return vec;
}

/*--------------------------------------------------------------------------------*/
int *read_fortran_ivec(FILE *infile, int *n)
{
  int nbytes;
  fread(&nbytes,1,sizeof(int),infile);
  *n=nbytes/sizeof(int);
  int *vec=ivector(0,*n-1);
  fread(vec,*n,sizeof(int),infile);
  int nb2=nbytes;
  fread(&nbytes,1,sizeof(int),infile);
  assert(nb2==nbytes);
  return vec;
}


/*--------------------------------------------------------------------------------*/
float *read_fortran_vec(FILE *infile, int *n)
{
  int nbytes;
  fread(&nbytes,1,sizeof(int),infile);
  *n=nbytes/sizeof(float);
  float *vec=vector(0,*n-1);
  fread(vec,*n,sizeof(float),infile);
  int nb2=nbytes;
  fread(&nbytes,1,sizeof(int),infile);
  assert(nb2==nbytes);
  return vec;
}


/*---------------------------------------------------------------------------*/
DumpData  *read_tt_dump_svec(char *fname)
{
  DumpData *data=(DumpData *)malloc(sizeof(DumpData));
  int n;
  FILE *infile=fopen(fname,"r");
  assert(infile!=NULL);
  
  int *nn=read_fortran_ivec(infile,&n);
  assert(n==1);
  printf("nn is %d\n",*nn);
  data->n=*nn;
  data->svec=read_fortran_dvec(infile,&n);
  assert(n==2*data->n);
  data->uvec=read_fortran_ivec(infile,&n);
  data->vvec=read_fortran_ivec(infile,&n);
  data->zvec=read_fortran_dvec(infile,&n);
  data->n_nonzero=0;
  for (int i=0;i<data->n;i++)
    if (data->zvec[i]>0)
      data->n_nonzero++;

  free(nn);
  fclose(infile);

  data->maxuv=0;
  for (int i=0;i<data->n;i++) 
    if (data->zvec[i]>0) {
      if (fabs(data->uvec[i])>data->maxuv)
	data->maxuv=fabs(data->uvec[i]);
      if (fabs(data->vvec[i])>data->maxuv)
	data->maxuv=fabs(data->vvec[i]);
    }

  return data;

}

/*---------------------------------------------------------------------------*/
complex double **complex_dmatrix(int n, int m)
{
  complex double *vec=(complex double *)malloc(sizeof(complex double)*n*m);
  complex double **mat=(complex double **)malloc(sizeof(complex double *)*n);
  for (int i=0;i<n;i++)
    mat[i]=vec+i*m;
  memset(vec,0,n*m*sizeof(complex double));
  return mat;
}
/*---------------------------------------------------------------------------*/
int read_tt_dump_rvec(DumpData *data, char *fname)
{
  FILE *infile=fopen(fname,"r");
  assert(infile);
  int nb;
  fread(&nb,1,sizeof(int),infile);

  fread(&(data->nbeam),1,sizeof(int),infile);
  int nn;
  fread(&nn,1,sizeof(int),infile);
  fread(&(data->du),1,sizeof(double),infile);
  fread(&nb,1,sizeof(int),infile);
  printf("nbeam and du are %4d %14.4f\n",data->nbeam,data->du);
  
  complex double *tmp = (complex double *)malloc(sizeof(complex double)*data->nbeam*data->nbeam);
  data->rgrid=(complex double ***)malloc(sizeof(complex double **)*data->n_nonzero);
  
  fread(&nb,1,sizeof(int),infile);
  assert(nb==sizeof(complex double)*data->n*data->nbeam*data->nbeam);
  int icur=0;
  for (int i=0;i<data->n;i++) {
    if (data->zvec[i]>0) {
      assert(icur<data->n_nonzero);
      data->rgrid[icur]=complex_dmatrix(data->nbeam,data->nbeam);
      fread(data->rgrid[icur][0],data->nbeam*data->nbeam,sizeof(complex double),infile);
      icur++;
    }
    else
      fread(tmp,data->nbeam*data->nbeam,sizeof(complex double),infile);

  }
  assert(icur==data->n_nonzero);  

  fread(&nb,1,sizeof(int),infile);
  assert(nb==sizeof(complex double)*data->n*data->nbeam*data->nbeam);
  printf("finished reading rgrid %s\n",fname);
  fclose(infile);
  free(tmp);

  data->gridsize=2*data->maxuv+data->nbeam+2;
  //data->pixel_size= 180.0*60.0*60.0/data->du/((double)data->gridsize)/3.1415926535897;
  data->pixel_size= 180.0/data->du/((double)data->gridsize)/3.1415926535897;
  
  printf("pixel size is %14.6e, grid size is %d\n",data->pixel_size,data->gridsize);

  return 0;
  
}
/*---------------------------------------------------------------------------*/
int dcut_zeros(double **vec, double *to_check, int n) 
{
  int nn=0;
  double *v1=*vec;
  for (int i=0;i<n;i++)
    if (to_check[i]>0)
      nn++;
  double *vv=dvector(0,nn-1);
  int icur=0;
  for (int i=0;i<n;i++)
    if (to_check[i]>0) {
      vv[icur]=v1[i];
      icur++;
    }
  free(v1);
  *vec=vv;
  assert(icur==nn);
  return nn;
  
}
/*---------------------------------------------------------------------------*/
int dcut_zeros_double_len(double **vec, double *to_check, int n) 
{
  int nn=0;
  double *v1=*vec;
  for (int i=0;i<n;i++)
    if (to_check[i]>0)
      nn++;
  double *vv=dvector(0,nn*2-1);
  int icur=0;
  for (int i=0;i<n;i++)
    if (to_check[i]>0) {
      vv[icur]=v1[2*i];
      icur++;
      vv[icur]=v1[2*i+1];
      icur++;
      
    }
  free(v1);
  *vec=vv;
  assert(icur==2*nn);

  return nn;
  
}
/*---------------------------------------------------------------------------*/
int icut_zeros(int **vec, double *to_check, int n) 
{
  int nn=0;
  int *v1=*vec;
  for (int i=0;i<n;i++)
    if (to_check[i]>0)
      nn++;
  int *vv=ivector(0,nn-1);
  int icur=0;
  for (int i=0;i<n;i++)
    if (to_check[i]>0) {
      vv[icur]=v1[i];
      icur++;
    }
  free(v1);
  *vec=vv;
  assert(icur==nn);
  return nn;
  
}
/*--------------------------------------------------------------------------------*/
void clear_dump_zeros(DumpData *dump)
{
  assert(dcut_zeros_double_len(&dump->svec,dump->zvec,dump->n)==dump->n_nonzero);
  assert(icut_zeros(&dump->uvec,dump->zvec,dump->n)==dump->n_nonzero);
  assert(icut_zeros(&dump->vvec,dump->zvec,dump->n)==dump->n_nonzero);


  //this one must come last because we're checking the zvec on the other guys.
  assert(dcut_zeros(&dump->zvec,dump->zvec,dump->n)==dump->n_nonzero);
}
/*---------------------------------------------------------------------------*/
SvecData *ReadSvecStruct(char *filename)
{
  FILE *infile=fopen(filename,"r");
  assert(infile);
  SvecData *data=(SvecData *)malloc(sizeof(SvecData));

  int n;
  int *crud=read_fortran_ivec(infile,&n);
  data->n=crud[0];
  free(crud);
  data->svec=read_fortran_dvec(infile,&n);
  data->uvec=read_fortran_vec(infile,&n);
  data->vvec=read_fortran_vec(infile,&n);
  data->zvec=read_fortran_dvec(infile,&n);
  double *radec=read_fortran_dvec(infile,&n);
  assert(n==2);
  data->ra=radec[0];
  data->dec=radec[1];

  fprintf(stderr,"ra/dec from svec are %14.7f %14.7f\n",data->ra*180/M_PI/15,data->dec*180/M_PI);
  //fprintf(stderr,"ra/dec from svec are %14.7f %14.7f\n",radec[0],radec[1]);
  free(radec);
  fclose(infile);
  return data;

}
/*---------------------------------------------------------------------------*/
int ReadSvec(double *data,char *filename)
{

  printf("reading %s\n",filename);
  FILE *infile;
  int nbytes_in,nest,junk;
  infile=fopen(filename,"r");
  if (!infile)
    {
      fprintf(stderr,"Unable to open file %s for checking size in ReadSvec.\n",filename);
      return 0;
    }
  fread(&nbytes_in,sizeof(int),1,infile);
  if (nbytes_in!=12)
    {
      fprintf(stderr,"Error - nybtes in is messed up inside of ReadSvec on %s.  Got %d\n",filename,nbytes_in);
      return 0;
    }
  fread(&nest,sizeof(int),1,infile);
  fread(&junk,sizeof(int),1,infile);
  fread(&junk,sizeof(int),1,infile);
  fread(&nbytes_in,sizeof(int),1,infile);

  /*fprintf(params->szef,"going to read %d estimators\n",nest);*/

  fread(&nbytes_in,sizeof(int),1,infile);
  if (nbytes_in !=nest*sizeof(double))
    {
      fprintf(stderr,"Got wrong chunk size in ReadSvec.  Expected %ld, got %d\n",nest*sizeof(double),nbytes_in);
      return 0;
    }
  /*fprintf(params->szef,"reading svec\n");*/
  fread(data,sizeof(double),nest,infile);
  /*fprintf(params->szef,"finished reading svec\n");*/

  fclose(infile);
  return nest;
}

/*---------------------------------------------------------------------------*/
int ReadNoise(double **noise, char *filename, int nest_in)
{
  gzFile infile;
  int nbytes_in,nest,junk,n,nx,i,j,jj,nr,ni;


  fprintf(stderr,"Reading %s\n",filename);
  infile=gzopen(filename,"r");
  if (!infile)
    {
      fprintf(stderr,"Error - unable to read noise file %s in ReadNoise\n",filename);
      /*return 0;*/
    }

  gzread(infile,&nbytes_in,sizeof(int));
  /*fprintf(params->szef,"Nbytes_in is %d\n",nbytes_in);*/

  gzread(infile,&nest,sizeof(int));
  if (nest!=nest_in)
    {
      fprintf(stderr,"Error - mismatch in sizes in ReadNoise.  Expected %d, got %d\n",nest_in,nest);
      return 0;
    }
  gzread(infile,&junk,sizeof(int));
  gzread(infile,&junk,sizeof(int));
  gzread(infile,&nbytes_in,sizeof(int));
  /*fprintf(params->szef,"Nbytes_in is %d\n",nbytes_in);*/

  nx=nest/2;
  for (j=0;j<nx;j++)
    {
      gzread(infile,&nbytes_in,sizeof(int));
      gzread(infile,&jj,sizeof(int));
      gzread(infile,&nbytes_in,sizeof(int));
      if (jj!=j+1)
	{
	  fprintf(stderr,"Error - mismatch in sizes in ReadNoise at j=%d, got %d\n",j,jj);
	  return 0;
	}
      nr=2*j+1;
      ni=2*j+2;
      
      gzread(infile,&nbytes_in,sizeof(int));
      gzread(infile,noise[nr-1],nr*sizeof(double));
      gzread(infile,&nbytes_in,sizeof(int));
      dflip_line(noise[nr-1],nr);

      gzread(infile,&nbytes_in,sizeof(int));
      gzread(infile,noise[ni-1],ni*sizeof(double));
      gzread(infile,&nbytes_in,sizeof(int));
      dflip_line(noise[ni-1],ni);
    }


  gzclose(infile);

  for (i=0;i<nest;i++)
    for (j=0;j<i;j++)
      noise[j][i]=noise[i][j];
  return nest;
}



/*---------------------------------------------------------------------------*/
int ReadCMBMat(double **noise, char *filename, int nest_in)
{
  gzFile infile;
  int nbytes_in,nest,junk,n,nx,i,j,jj,nr,ni,ipol,nmat;
  float bandmin,bandmax;

  infile=gzopen(filename,"r");
  if (!infile)
    {
      fprintf(stderr,"Error - unable to read noise file %s in ReadNoise\n",filename);
      /*return 0;*/
    }

  gzread(infile,&nbytes_in,sizeof(int));
  /*fprintf(stderr,"Nbytes_in is %d\n",nbytes_in);*/

  gzread(infile,&nest,sizeof(int));
  if (nest!=nest_in)
    {
      fprintf(stderr,"Error - mismatch in sizes in ReadNoise.  Expected %d, got %d\n",nest_in,nest);
      return 0;
    }
  gzread(infile,&junk,sizeof(int));
  gzread(infile,&junk,sizeof(int));
  gzread(infile,&nbytes_in,sizeof(int));
  /*fprintf(params->szef,"Nbytes_in is %d\n",nbytes_in);*/


  gzread(infile,&nbytes_in,sizeof(int));
  gzread(infile,&ipol,sizeof(int));
  gzread(infile,&nmat,sizeof(int));
  gzread(infile,&nbytes_in,sizeof(int));
  if (nmat!=1)
    {
      fprintf(stderr,"Error - Signal file has %d matrices in it, we should only have 1.\n",nmat);
      return 0;
    }

  gzread(infile,&nbytes_in,sizeof(int));
  gzread(infile,&bandmin,sizeof(float));
  gzread(infile,&nbytes_in,sizeof(int));
  gzread(infile,&nbytes_in,sizeof(int));
  gzread(infile,&bandmax,sizeof(float));
  gzread(infile,&nbytes_in,sizeof(int));
  
  //fprintf(params->szef,"Signal matrix spans %8.2g to %8.2g\n",bandmin,bandmax);


  nx=nest/2;
  for (j=0;j<nx;j++)
    {
      gzread(infile,&nbytes_in,sizeof(int));
      gzread(infile,&jj,sizeof(int));
      gzread(infile,&nbytes_in,sizeof(int));
      if (jj!=j+1)
	{
	  fprintf(stderr,"Error - mismatch in sizes in ReadNoise at j=%d, got %d\n",j,jj);
	  return 0;
	}
      nr=2*j+1;
      ni=2*j+2;
      
      gzread(infile,&nbytes_in,sizeof(int));
      gzread(infile,noise[nr-1],nr*sizeof(double));
      gzread(infile,&nbytes_in,sizeof(int));
      dflip_line(noise[nr-1],nr);

      gzread(infile,&nbytes_in,sizeof(int));
      gzread(infile,noise[ni-1],ni*sizeof(double));
      gzread(infile,&nbytes_in,sizeof(int));
      dflip_line(noise[ni-1],ni);
    }


  gzclose(infile);

  for (i=0;i<nest;i++)
    for (j=0;j<i;j++)
      noise[j][i]=noise[i][j];
  return nest;
}


/*----------------------------------------------------------------------*/
void dflip_line(double *line, int nelem)
{
  int i;
  double d;
  for (i=0;i<nelem/2;i++)
    {
      d=line[i];
      line[i]=line[nelem-1-i];
      line[nelem-1-i]=d;
    }
}
/*----------------------------------------------------------------------*/
void dflip_mat(double **mat, int nelem) 
{
  for (int i=0;i<nelem;i++)
    dflip_line(mat[i],nelem);
}
/*----------------------------------------------------------------------*/

char *read_all_stdin(char *filename, SZControlParams *params)
{
  FILE *infile;
  char line_in[SZSTRINGLEN],*big_line,*cur_spot;
  int line_len,total_len,i,finished;

  infile=fopen(filename,"r");
  if (!infile) {
    if (!params->quiet)
      fprintf(params->szef,"unable to open parameter file %s for reading.\n",filename);
    return NULL;
  }

  total_len=1;
  big_line=scvector(0,1);
  sprintf(big_line,"");

  finished=0;
  while ((fgets(line_in,SZSTRINGLEN-1,infile))&&(finished==0))
    {
      /*fprintf(params->szef,"line_in is %s\n",line_in);*/
      if (strncmp(line_in,"#finished",strlen("#finished")-1)==0)
        finished=1;
      for (i=0;i<(int)strlen(line_in);i++)
        if ((line_in[i]=='#')||(line_in[i]=='!')||(line_in[i]=='\n'))
          line_in[i]='\0';
      /*fprintf(params->szef,"line_in is .%s.\n",line_in);*/
      lengthen_string(&big_line,&total_len,strlen(line_in)+3);
      sprintf(big_line,"%s %s",big_line,line_in);
    }
  return big_line;
}
/*--------------------------------------------------------------------------------*/
void lengthen_string(char **cur_string,int *len_in, int n_to_add)
{
  char *new_string;
  int len;

  len=*len_in;
  new_string=scvector(0,len+n_to_add);
  strncpy(new_string,*cur_string,len);

  if (len>0)
    free(*cur_string);
  *len_in = len+n_to_add;
  *cur_string=new_string;
}
/*---------------------------------------------------------------------------*/
int WriteDiff(SZControlParams *params, double *data)
{
  int i,junk,nest,nbytes_in;
  double *diff;
  FILE *outfile,*infile;
  char filename[SZSTRINGLEN];

  diff=dvector(0,params->nest-1);
  for (i=0;i<params->nest;i++)
    diff[i]=params->data_org[i]-data[i];

  /*
  sprintf(filename,"%s.txt",params->DifferenceDataFile);
  outfile=fopen(filename,"w");
  for (i=0;i<params->nest;i++)
    fprintf(outfile,"%14.6e %14.6e %14.6e\n",params->data_org[i],data[i],diff[i]);
  fclose(outfile);
  */
  
  sprintf(filename,"%s_SVEC",params->GridderDataBase);
  infile=fopen(filename,"r");
  if (!outfile)
    {
      fprintf(stderr,"unable to open %s for reading template in WriteDiff\n",filename);
      return -2;
    }

  outfile=fopen(params->DifferenceDataFile,"w");
  if (!outfile)
    {
      fprintf(stderr,"unable to open %s for writing in WriteDiff\n",params->DifferenceDataFile);
      return -1;
    }
  
  fread(&nbytes_in,sizeof(int),1,infile);
  if (nbytes_in!=12)
    {
      fprintf(stderr,"Error - nybtes in is messed up inside of ReadSvec on %s.  Got %d\n",filename,nbytes_in);
      return -3;
    }
  fwrite(&nbytes_in,sizeof(int),1,outfile);
  fread(&nest,sizeof(int),1,infile);
  fwrite(&nest,sizeof(int),1,outfile);
  fread(&junk,sizeof(int),1,infile);
  fwrite(&junk,sizeof(int),1,outfile);
  fread(&junk,sizeof(int),1,infile);
  fwrite(&junk,sizeof(int),1,outfile);
  fread(&nbytes_in,sizeof(int),1,infile);
  fwrite(&nbytes_in,sizeof(int),1,outfile);
 
  if (params->nest!=nest)
    {
      fprintf(stderr,"Error in writing difference - expected %d estimators, got %d\n",params->nest,nest);
      return -4;
    }
  
  fread(&nbytes_in,sizeof(int),1,infile);
  if (nbytes_in!=sizeof(double)*nest)
    {
      fprintf(stderr,"Size mismatch in reading svec - expected %d, got %d\n",nest*sizeof(double),nbytes_in);
      return -5;
    }
  fwrite(&nbytes_in,sizeof(int),1,outfile);
  fwrite(diff,sizeof(double),nest,outfile);
  fwrite(&nbytes_in,sizeof(int),1,outfile);

  /*now we'll use diff as workspace.  pardon the ugliness...*/
  fread(diff,sizeof(double),nest,infile);
  fread(&nbytes_in,sizeof(int),1,infile);


  /*do u*/
  fread(&nbytes_in,sizeof(int),1,infile);
  if (nbytes_in!=sizeof(float)*nest/2)
    {
      fprintf(stderr,"Size mismatch in reading uvec - expected %d, got %d\n",nest*sizeof(float)/2,nbytes_in);
      return -6;
    }
  fread(diff,sizeof(float),nest/2,infile);
  fread(&nbytes_in,sizeof(int),1,infile);

  fwrite(&nbytes_in,sizeof(int),1,outfile);
  fwrite(diff,sizeof(float),nest/2,outfile);
  fwrite(&nbytes_in,sizeof(int),1,outfile);

  /*do v*/
  fread(&nbytes_in,sizeof(int),1,infile);
  if (nbytes_in!=sizeof(float)*nest/2)
    {
      fprintf(stderr,"Size mismatch in reading vvec - expected %d, got %d\n",nest*sizeof(float)/2,nbytes_in);
      return -7;
    }

  fread(diff,sizeof(float),nest/2,infile);
  fread(&nbytes_in,sizeof(int),1,infile);

  fwrite(&nbytes_in,sizeof(int),1,outfile);
  fwrite(diff,sizeof(float),nest/2,outfile);
  fwrite(&nbytes_in,sizeof(int),1,outfile);

  /*do weights*/
  fread(&nbytes_in,sizeof(int),1,infile);
  if (nbytes_in!=sizeof(double)*nest/2)
    {
      fprintf(params->szef,"Size mismatch in reading zvec - expected %d, got %d\n",nest*sizeof(double)/2,nbytes_in);
      return -8;
    }
  fread(diff,sizeof(double),nest/2,infile);
  fread(&nbytes_in,sizeof(int),1,infile);

  fwrite(&nbytes_in,sizeof(int),1,outfile);
  fwrite(diff,sizeof(double),nest/2,outfile);
  fwrite(&nbytes_in,sizeof(int),1,outfile);

  fclose(infile);
  fclose(outfile);
  free(diff);
  return 0;
}
/*--------------------------------------------------------------------------------*/
int append_tag(char *str1, char *str2, int maxlen)
{
  int len1,len2,len;
  len1=strlen(str1);
  len2=strlen(str2);
  len=len1+len2;
  if ((len<maxlen-1)&&(len1>=0)&&(len2>=0))
    {
      sprintf(str1,"%s%s",str1,str2);
      return 0;
    }
  else
    {
      fprintf(stderr,"Incoming string funniness in append_tag.  Skipping appendage\n");
      fprintf(stderr,"Lengths are %d %d, maximum is %d\n",len1,len2,maxlen);
      return -1;

    }  
}
/*--------------------------------------------------------------------------------*/
int read_szmodel_file(char *filename, double **x_out, double **y_out, double **err_out)
{
  char line_in[SZSTRINGLEN];
  int i,n;
  FILE *infile;
  double *x,*y,*err;

  infile=fopen(filename,"r");
  if (!infile)
    {
      fprintf(stderr,"Error opening %s for reading in read_szmodel_file.\n",filename);
      exit(EXIT_FAILURE);
    }
  n=0;
  
  while (fgets(line_in,SZSTRINGLEN,infile))
    {
      if (strlen(line_in)>3)
	n++;
    }

  rewind(infile);
  x=dvector(0,n-1);
  y=dvector(0,n-1);
  err=dvector(0,n-1);
  i=0;
  while (fgets(line_in,SZSTRINGLEN,infile))
    {
      if (strlen(line_in)>3)
	{
	  sscanf(line_in,"%lf%lf%lf",&x[i],&y[i],&err[i]);
	  i++;
	}
    }
  fclose(infile);
  *x_out=x;
  *y_out=y;
  *err_out=err;

  if (i!=n)
    {
      fprintf(stderr,"Wierdness in read_sz_model_file, n changed from %d to %d\n",n,i);
      exit(EXIT_FAILURE);
    }
  return n;
  
}
/*--------------------------------------------------------------------------------*/
int big_initialize_sz(BigSZControlParams *params)
/*initialize multiple SZ datasets*/
{
  int i,j,*nfound,info;
  double **x, **y,**err;
  char tag[SZSTRINGLEN];
  FILE *outfile;
  if (params->ndataset<=0)
    {
      fprintf(stderr,"No SZ Datafiles requested.  Skipping...\n");
      return 0;
    }
  
  params->szdatasets=(SZControlParams *)malloc(params->ndataset*sizeof(SZControlParams));

  for (i=0;i<params->ndataset;i++)    
    {
      params->szdatasets[i].node = params->node;
      params->szdatasets[i].quiet = params->quiet;
      strncpy(params->szdatasets[i].ParamFile,params->paramfiles[i],SZSTRINGLEN-1);
      sprintf(params->szdatasets[i].JacoTag,"JACO_%d",i);
      
      info=initialize_sz(&(params->szdatasets[i]));
      /*if we hit an error, run away.*/
      if (info)
	{
	  fprintf(stderr,"Error in sz initialization %d with code %d from %s\n",i,info,params->szdatasets[i].ParamFile);
	  return info;
	}
    }
  /*OK - everybody is initialized now.*/

  if ((strlen(params->DerotDataFile) <=0)||(strlen(params->DerotDataFile)>=SZSTRINGLEN-1))
    sprintf(params->DerotDataFile,"szdata.dat.000");
  
  // Andi: only one of the nodes gets to write the MASTER data file.
  // Otherwise they overwrite each other
  // < 2 covers node 0 (1 CPU) or node 1 (multiple CPUs)

  if (params->node < 2) {
    outfile=fopen(params->DerotDataFile,"w");
    if (!outfile)
      {
	fprintf(stderr,"Unable to open %s for writing master derotated data file in big_initialize_sz.\n",params->DerotDataFile);
	exit(EXIT_FAILURE);
      }
  }

  params->nest=0;
  for (i=0;i<params->ndataset;i++)
    {
      params->nest+=params->szdatasets[i].nest;
      for (j=0;j<params->szdatasets[i].nest;j++)
	if (params->node < 2) 
	  fprintf(outfile,"%4d %18.10e %18.10e\n",j,params->szdatasets[i].data[j],sqrt(params->szdatasets[i].eigvals[j]));
    }
  if (params->node < 2)
    fclose(outfile);

  params->racent=dvector(0,params->ndataset-1);
  params->deccent=dvector(0,params->ndataset-1);
  params->pixsize=dvector(0,params->ndataset-1);
  params->npix=ivector(0,params->ndataset-1);
  params->ymap=(double ***)malloc(sizeof(double **)*params->ndataset);
  for (i=0;i<params->ndataset;i++)
    {
      params->ymap[i]=params->szdatasets[i].ymap;
      params->racent[i]=params->szdatasets[i].racent;
      params->deccent[i]=params->szdatasets[i].deccent;
      params->pixsize[i]=params->szdatasets[i].pixsize;      
      params->npix[i]=params->szdatasets[i].npix;

      fprintf(stderr,"In big initialize, set npix[%d] to %d %12.4g\n",i,params->npix[i],params->pixsize[i]);
      fprintf(stderr,"In big initialize, set ra/dec centers to  %14.6f %14.6f\n",params->racent[i],params->deccent[i]);
    }
  
  return 0;
  
  
}
/*--------------------------------------------------------------------------------*/
int big_get_sz_model(double *model, BigSZControlParams *params)
{
  int i,j,ncur,info=0;
  
  //44/60 vs. 48/59
  #pragma omp parallel for shared(params,stderr,model) private(info,ncur,j) default(none)  
  for (i=0;i<params->ndataset;i++)
    {      
      ncur=0;
      for (j = 0; j <i; ++j)
	ncur += params->szdatasets[j].nest;
#ifdef USE_RGRID
      info=get_sz_model_dump(&model[ncur],&(params->szdatasets[i]));      
      //info=get_sz_model_dump(&model[ncur],params->szdatasets[i].ymap, params->szdatasets[i].dump);
      //FILE *outfile=fopen("sz_model_dump_new.txt","w");
#else
      info=get_sz_model(&model[ncur],&(params->szdatasets[i]));
      //FILE *outfile=fopen("sz_model_dump_old.txt","w");

#endif
      //for (int j=0;j<params->szdatasets[i].nest;j++)
      //fprintf(outfile,"%4d %16.6e\n",j,model[j]);
      //fclose(outfile);
      //ncur += params->szdatasets[i].nest;

      
      if (info)
	  fprintf(stderr,"info is %d on %d\n",info,i);

    }
  return info;
  
}



/*--------------------------------------------------------------------------------*/
void fftshift(double **mat, int n, int m)
{
  double **mat_copy=dmatrix(0,n-1,0,m-1);

  for (int i=0;i<n;i++)
    for (int j=0;j<m;j++) {
      int ii=i+n/2;
      if (ii>=n)
	ii-=n;
      int jj=j+m/2;
      if (jj>=m)
	jj-=m;
      mat_copy[ii][jj]=mat[i][j];      
    }
  memcpy(mat[0],mat_copy[0],sizeof(double)*n*m);
  


  free(mat_copy[0]);
  free(mat_copy);
}


/*--------------------------------------------------------------------------------*/
void ifftshift(double **mat, int n, int m)
{
  double **mat_copy=dmatrix(0,n-1,0,m-1);

  for (int i=0;i<n;i++)
    for (int j=0;j<m;j++) {
      int ii=i+n/2;
      if (ii>=n)
	ii-=n;
      int jj=j+m/2;
      if (jj>=m)
	jj-=m;
#ifdef DEBUG_SZ
      assert(ii>=0);
      assert(ii<n);
      assert(jj>=0);
      assert(jj<m);
#endif
      mat_copy[i][j]=mat[ii][jj];      
    }
  memcpy(mat[0],mat_copy[0],sizeof(double)*n*m);
  


  free(mat_copy[0]);
  free(mat_copy);
}


/*--------------------------------------------------------------------------------*/
void cfftshift(complex double **mat, int n, int m)
{
  complex double **mat_copy=complex_dmatrix(n,m);

  for (int i=0;i<n;i++)
    for (int j=0;j<m;j++) {
      int ii=i+n/2;
      if (ii>=n)
	ii-=n;
      int jj=j+m/2;
      if (jj>=m)
	jj-=m;
#ifdef DEBUG_SZ
      assert(ii>=0);
      assert(ii<n);
      assert(jj>=0);
      assert(jj<m);
#endif
      mat_copy[ii][jj]=mat[i][j];      
    }
  memcpy(mat[0],mat_copy[0],sizeof(complex double)*n*m);
  


  free(mat_copy[0]);
  free(mat_copy);
}
/*--------------------------------------------------------------------------------*/
void cifftshift(complex double **mat, int n, int m)
{
  complex double **mat_copy=complex_dmatrix(n,m);

  for (int i=0;i<n;i++)
    for (int j=0;j<m;j++) {
      int ii=i+n/2;
      if (ii>=n)
	ii-=n;
      int jj=j+m/2;
      if (jj>=m)
	jj-=m;
#ifdef DEBUG_SZ
      assert(ii>=0);
      assert(ii<n);
      assert(jj>=0);
      assert(jj<m);
#endif

      mat_copy[i][j]=mat[ii][jj];      
    }
  memcpy(mat[0],mat_copy[0],sizeof(complex double)*n*m);
  


  free(mat_copy[0]);
  free(mat_copy);
}


/*--------------------------------------------------------------------------------*/
void print_mat(double **mat, int n, int m)
{
  for (int i=0;i<n;i++) {
    for (int j=0;j<m;j++)
      printf("%14.4g ",mat[i][j]);
    printf("\n");
  }
}




/*--------------------------------------------------------------------------------*/
void unpack_cmat(double complex **mat, int n, int m)
{
  double complex **cp=complex_dmatrix(n,m);
  memset(cp[0],0,sizeof(double complex)*n*m);
  int i,j;
  long ind=0;
  for (j=0;j<m/2+1;j++) {    
    cp[0][j]=mat[0][j];
    ind++;

  }
  for (j=m/2+1;j<m;j++) {
#ifdef DEBUG_SZ
    assert((m-j)>=0);
    assert((m-j)<m);
#endif
    cp[0][j]=conj(cp[0][m-j]);
  }
  
  
  for (i=1;i<n;i++) {
    for (j=0;j<m/2+1;j++) {    
      cp[i][j]=mat[0][ind];
#ifdef DEBUG_SZ
      assert(ind<n*m);
#endif
      ind++;
    }
  }     
  for (i=1;i<n;i++) 
    for (j=m/2+1;j<m;j++) {
      cp[i][j]=conj(cp[n-i][m-j]);
      //cp[i][j]=3.234;
    }
  
  
  
  memcpy(mat[0],cp[0],sizeof(double complex)*n*m);

  free(cp[0]);
  free(cp);
    
}

/*--------------------------------------------------------------------------------*/
void fft_2d_map_wshifts(double **mat, complex double **matft, int n, int m)
{

  fftw_plan p;

  #pragma omp critical
  {
    p=fftw_plan_dft_r2c_2d(n,m,mat[0],matft[0],FFTW_ESTIMATE);
  }
  ifftshift(mat,n,m);
  fftw_execute(p);
  #pragma omp critical
  {
    fftw_destroy_plan(p);
  }
  unpack_cmat(matft,n,m);
  ifftshift(mat,n,m);
  cfftshift(matft,n,m);
  
}
/*--------------------------------------------------------------------------------*/
void print_complex(char *format, double complex val)
{
  printf(format,creal(val));
  if (cimag(val)<0) {
    printf(" -");
    printf(format,-cimag(val));
  }
  else {
    printf(" +");
    printf(format,cimag(val));
  }
  printf("i");
}



/*--------------------------------------------------------------------------------*/
void print_cmat(double complex **mat, int n, int m)
{
  for (int i=0;i<n;i++) {
    for (int j=0;j<m;j++) {
      printf("  ");
      print_complex("%12.5g",mat[i][j]);
    }
    printf("\n");
  }
  printf("\n");
}


/*================================================================================*/
#endif
//build in a simple main just for testing purposes
#if 0

int main(int argc, char *argv[])
{

  printf("hello world.\n");
  DumpData *dump=read_tt_dump_svec("/home/sievers/clusters/a1689/cres/a1689_cbi1_63x2_dump_Svec.dmp");
  read_tt_dump_rvec(dump,"/home/sievers/clusters/a1689/cres/a1689_cbi1_63x2_dump_Rvec.dmp");
  int oversamp=1;
  dump->npix_use=get_npix_from_oversamp(dump->gridsize,oversamp);
  dump->pixel_size/=(1.0+2.0*(double)oversamp);
  clear_dump_zeros(dump);


  SvecData *svec=ReadSvecStruct("/home/sievers/clusters/a1689/cres/a1689_cbi1_63x2_SVEC");
  printf("Ra and dec are %14.6f %14.6f\n",svec->ra,svec->dec);


  double **map=dmatrix(0,dump->npix_use-1,0,dump->npix_use-1);
  memset(map[0],0,sizeof(double)*dump->npix_use*dump->npix_use);
  
  int ii=90*(1+2*oversamp)+1;
  int jj=90*(1+2*oversamp)+2;
  printf("ii and jj are %d %d\n",ii,jj);
  map[ii][jj]=1;


  double *vec=dvector(0,dump->n_nonzero*2-1);
  get_sz_model_dump(vec,map,dump);
  printf("data[1,2]=%14.6e %14.6e\n",vec[0],vec[1]);
  printf("data[9,10]=%14.6e %14.6e\n",vec[8],vec[9]);

  int ind=76;
  printf("vals[%d] are %14.4g %14.4g %14.4g %14.4g\n",ind,svec->svec[ind],dump->svec[ind],dump->zvec[ind/2],dump->svec[ind]/dump->zvec[ind/2]);


  int n=atoi(argv[1]);
  int m=atoi(argv[2]);
  double **mat=dmatrix(0,n-1,0,m-1);
  for (int i=0;i<n;i++)
    for (int j=0;j<m;j++)
      mat[i][j]=10*i+j;

  fftshift(mat,n,m);
  ifftshift(mat,n,m);
  print_mat(mat,n,m);
  complex double **cmat=complex_dmatrix(n,m);
  fft_2d_map_wshifts(mat,cmat,n,m);
  printf("\n\n");
  print_cmat(cmat,n,m);
  return 0;
}
#endif


/*--------------------------------------------------------------------------------*/
void transpose_matrix(double **mat, int n)
{
  for (int i=0;i<n;i++) 
    for (int j=i;j<n;j++) {
      double val=mat[i][j];
      mat[i][j]=mat[j][i];
      mat[j][i]=val;
    }
}
