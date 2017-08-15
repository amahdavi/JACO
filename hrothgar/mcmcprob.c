#include <string.h>
#include <strings.h>
#include <math.h>
#include <stdio.h>
#include <sys/times.h>
#include <sys/types.h>
#include <sys/stat.h>
#include <sciutils.h>
#include <hrothgar.h>
#include <hrothgar_proto.h>
#define _XOPEN_SOURCE 500 // Required by SUSv2 for usleep to work.
#include <unistd.h>
#include <stdlib.h>
#include <zlib.h>
#include <omp.h>
#include <regex.h>
#include <gsl/gsl_roots.h>

#ifdef HAVE_LIBREADLINE

#include <readline/readline.h>
#include <readline/history.h>

char *getstring(char *prompt) {

  static char *line_read = (char *)NULL;

  if (line_read) {
    free(line_read);
    line_read = (char *)NULL;
  }

  line_read = readline(prompt);
  if (line_read && *line_read) add_history(line_read);

  return line_read;
}

#else 

char *getstring() {

  static char f[80];

  printf("prompt");
  fgets(f,79,stdin);
  f[strlen(f)-1]=0;

  return f;
}

#endif

#ifdef HAVE_LIBCPGPLOT
#include <cpgplot.h>
#endif

#define MAXSET 5

//int usleep(__useconds_t __useconds);

int read_files(gzFile *files, int nfiles, char ***paramnames,
	       double ***data, int *totcols, long *count,
	       long *countmax, long burnin, int fromstart, int *new)
{
  int i,j,firstfile;
  char *tmpstr,line[10];
  double *row;
  
  fflush(stdout);
  gzFile mcfile;
  long lastpos;

  for (i = 0; i < nfiles; ++i) {
    
    mcfile = files[i];
    
    if (fromstart) {
      
      gzrewind(mcfile);
      gzread(mcfile,line,9*sizeof(char));
      line[8] = 0;
      if (strcmp(line,"HROTHGAR")) {
	BYE("Sorry -- input file not a Hrothgar MCMC File");
      }
      gzread(mcfile,&j,sizeof(int));
      if (*totcols < 0)
	*totcols = j;
      else
	if (*totcols != j) 
	  BYE("Mismatch in number of columns %d, %d",*totcols,j);
      
      if ((*paramnames) == NULL) {
	firstfile=1;
	(*paramnames) = sci_stringmatrix((*totcols)+1,200);
      } else
	firstfile=0;

      strcpy((*paramnames)[0],"chisq");
      strcpy((*paramnames)[(*totcols)],"ord");
      int plength;
      for (j = 1; j < (*totcols); ++j) {
	gzread(mcfile,&plength,sizeof(int));

	if (plength > 200 || plength < 1)
	  BYE("Invalid header length %d",plength);

	tmpstr = (char *)malloc(plength*sizeof(char));
	gzread(mcfile,tmpstr,plength*sizeof(char));
	if (firstfile) 
	  strcpy((*paramnames)[j],tmpstr);
	else {
	  if (strcmp((*paramnames)[j],tmpstr)) {
	    BYE("Error---mismatch in MCMC file headers (%s,%s)",
		*paramnames[j],tmpstr);
	  }
	}
	free(tmpstr);
      }
      
      if ((*data) == NULL)
	(*data) = sci_dmatrix(*totcols+1,*countmax);

      row = sci_dvector(*totcols);
      for (j = 0; j < burnin; ++j)
	gzread(mcfile,row,(*totcols)*sizeof(double));
      free(row);

    }

    row = sci_dvector(*totcols);
    
    lastpos = gztell(mcfile);
    while (gzread(mcfile,row,sizeof(double)*(*totcols))
	   == sizeof(double)*(*totcols)) {
      
      lastpos = gztell(mcfile);
      ++(*count);
      if (*count == *countmax-4) {
	*countmax += 10000;
	for (j = 0; j < *totcols; ++j) 
	  sci_dvector_resize(&(*data)[j],*countmax);
      }
      
      for (j = 0; j < *totcols; ++j) 
	(*data)[j][(*count)-1] = row[j];
      *new = 1;
    }
    free(row);
    gzrewind(mcfile);
    gzseek(mcfile,lastpos,SEEK_SET);
  }

  sci_dvector_resize(&(*data)[*totcols],*count);
  for (i = 0; i < *count; ++i)
    (*data)[*totcols][i] = (double) i;
      
  return 0;
}

// If *matches is not null, perform regexp match, stores matches in matches,
// and return number. If matches is null, perform exact match, and return
// match.
int match_column(char *col, char **paramnames, int totcols, int *matches)
{
  int i, n, nmatch;
  char colname[80],parname[80], *delimit;

  delimit = index(col,'@');
  if (delimit == NULL) n = 0; else n = strlen(delimit);
  strncpy(colname,col,79);
  colname[strlen(col)-n] = 0;

  regex_t preg;
  // Look for exact matches of portions prior to "@"
  if (matches != NULL)
    if (regcomp(&preg,colname,REG_EXTENDED | REG_NOSUB)) {
      BYE("%s is not a valid regular expression ",colname);
    }
  nmatch = 0;
  for (i = 0; i < totcols; ++i) {
    delimit = index(paramnames[i],'@');
    if (delimit == NULL) n = 0; else n = strlen(delimit);
    strncpy(parname,paramnames[i],79);
    parname[strlen(paramnames[i])-n] = 0;

    delimit = index(paramnames[i],'-');
    if (delimit == NULL) n = 0; else n = strlen(delimit);
    parname[strlen(paramnames[i])-n] = 0;

    if (matches != NULL) {
      if (regexec(&preg,parname,0,NULL,0) == 0) {
	matches[nmatch] = i;
	++nmatch;
      }
    } else
      if (strcmp(colname,parname) == 0)
	return i;
  }

  if (matches == NULL) return -1;

  return nmatch;
}

int get_stats(double *x, long count,
	      double *mean, double *sigma, double *l1, double *l2,
	      double quantile, long chimin_ind)
{

  int j;
  double q1 = 0.5-quantile/2;
  double q2 = 0.5+quantile/2;

  size_t *ord = sci_sizetvector(count);
  double *temp = sci_dvector(count);

  gsl_sort_index(ord,x,1,count);
      
  for (j = 0; j < count; ++j) {
    temp[j] = x[ord[j]];
    //if ((j % 100) == 1)
    //  printf("%f %.3E\n",1.*j/count,temp[j]);
  }

  int usevariance=(chimin_ind < -1);
  if (chimin_ind > -1) {
    j = gsl_interp_bsearch(temp,x[chimin_ind],0,count-1);
    double q0 = (double)j/count;
    q1 = q0-quantile/2.;
    q2 = q0+quantile/2.;
    if (q1 < 0 || q2 > 1) {
      q1 = (1.-quantile)/2.;
      q2 = 1.-(1.-quantile)/2.;
    }
    //printf("***** %f %f %f %f *****\n",x[chimin_ind],q0,q1,q2);
    if (q1 > q0 || q2 < q0) {
      /* printf("Warning---MCMC distribution is too skewed for proper\n"); */
      /* printf("estimation using max-likelihood '-d' method. Options are to\n"); */
      /* printf("run longer chains, or else to use variance\n"); */
      /* printf("of the distribution instead (method '-v'), which I am \n"); */
      /* printf("currently using for the parameter below only:\n"); */
      usevariance=1;
      //if (q1 < 0) { q2 = quantile; q1 = 0; }
      //else if (q2 > 1) { q1 = 1.-quantile; q2 = 1.; }
    }
    //printf("%E ---  %E %E %E\n",x[chimin_ind],q0,q1,q2);
    *mean = x[chimin_ind];
  }
  else 
    *mean = gsl_stats_mean(temp,1,count);
  
  
  *l1 = gsl_stats_quantile_from_sorted_data(temp,1,count,q1);
  *l2 = gsl_stats_quantile_from_sorted_data(temp,1,count,q2);
  //printf("--------%f (%f) %f (%f)\n",q1,*l1,q2,*l2);


  if (!usevariance) {
    *sigma = (*l2-*l1)/2.;
  }
  else
    *sigma = gsl_stats_sd(temp,1,count);

  free(ord);
  free(temp);

  return 0.;
}

int gaussian_blur(float *input,int gridsize,int width)
{

  static int oldwidth=0;
  static double *gf=NULL;
  int i,j,k;
  double diff, norm=0.;

  if (oldwidth != width) {
    if (gf != NULL) free(gf);
    gf = sci_dvector(2*width+1);
    for (i = -width; i <= width; ++i) {
      diff = 3.*i/width;
      diff = exp(-diff*diff/2.);
      gf[i+width] = diff;
    }
    oldwidth = width;
  }
  
  double **newf = sci_dmatrix(gridsize,gridsize);
  double **newg = sci_dmatrix(gridsize,gridsize);
  for (i = 0; i < gridsize; ++i)
    for (j = 0; j < gridsize; ++j) {
      for (k = -width; k <= width; ++k)
	if (j+k >= 0 && j+k < gridsize)
	  newf[i][j] += input[i+gridsize*(j+k)]*gf[width+k];
    }
  for (i = 0; i < gridsize; ++i)
    for (j = 0; j < gridsize; ++j) 
      for (k = -width; k <= width; ++k)
	if (i+k >= 0 && i+k < gridsize) {
	  diff = newf[i+k][j]*gf[width+k];
	  norm += diff;
	  newg[i][j] += diff;
	} 

  for (i = 0; i < gridsize; ++i) {
    for (j = 0; j < gridsize; ++j) { 
      input[i+gridsize*j] = newg[i][j]/norm;
    }
  }

  sci_free_dmatrix(newf,gridsize);
  sci_free_dmatrix(newg,gridsize);

  return 0;
}

struct cpars {
  double cdelta;
  double delta;
};

double c200_zero_func(double c200, void *params)
{
  struct cpars *c = (struct cpars *)params;
  
  double answer = c200/c->cdelta;
  answer = pow(answer,5.)*(200./c->delta);
  answer *= 1.-1./(1.+c->cdelta)-log(1.+c->cdelta);
  answer /= 1.-1./(1.+c200)-log(1.+c200);
  return answer-1.;

}

double get_c200_from_cdelta(double cdelta, double delta)
{
  struct cpars c;

  gsl_root_fsolver *s = gsl_root_fsolver_alloc(gsl_root_fsolver_falsepos);
  gsl_function F;
  F.function = c200_zero_func;

  c.delta = delta;
  c.cdelta = cdelta;
  F.params = &c;

  gsl_root_fsolver_set(s,&F,0.001,100.);

  int i=0, status;
  double r1, r2;

  do {
    ++i;
    status = gsl_root_fsolver_iterate(s);
    r1 = gsl_root_fsolver_x_lower(s);
    r2 = gsl_root_fsolver_x_upper(s);
    //printf("---%f %f %f %f\n",r1,r2,mdiffcalc(r1,&od),mdiffcalc(r2,&od));
    status = gsl_root_test_delta(r1,r2,0.,0.001);
  } while (status == GSL_CONTINUE && i < 500);

  if (status != GSL_SUCCESS)
    BYE("Error in c200");

  return r1;


}

int main(int argc, char *argv[])
{
  long i, j, count[MAXSET], ngrid, ngridsq, k0, smooth=6, x0, y0, x, y, k, nsig;
  long burnin=5000,chimin_ind=-2;
  int col1, col2, pos=1, doxlog=0, doylog=0, withplots=0,contour,postscript;
  double xdel, ydel, smoothsq, cumul, sum=0.;
  double quantile=0.684,chimin[MAXSET];
  double *chisq=NULL,**data[MAXSET];
  gzFile mcfile,files[MAXSET][100];
  char line[10],**paramnames[MAXSET],*mylabel;
  char *outf;
  int nset = 1;
  int nfiles[MAXSET];
  char *pfile;
  float mytextsize=1.5;
  char *delimit;

  nfiles[0] = 0;

  int datafile=0;
  int usemode=1;
  for (i = 1; i < argc; ++i) {
    if (strcmp(argv[i],"-i")==0) {
      withplots=1;
      postscript=0;
      burnin = atol(argv[i+1]);
      i = argc;
    }
    else if (strcmp(argv[i],"-p")==0) {
      withplots=1;
      postscript=1;
      pfile = argv[i+1];
      delimit = index(pfile,'@');
      mylabel=NULL;
      if (delimit) {
	mylabel = &delimit[1];
	delimit[0]=0;
      }
      
      mytextsize=1.5;
      
      burnin = atol(argv[i+2]);
      ++pos;
      i = argc;
    }
    else if (strcmp(argv[i],"-a")==0) {
      burnin = atol(argv[i+1]);
      i = argc;
    }
    else if (strcmp(argv[i],"-d")==0) {
      burnin = atol(argv[i+1]);
      datafile = 1;
      i = argc;
    }
    else if (strcmp(argv[i],"-v")==0) {
      burnin = atol(argv[i+1]);
      datafile = 1;
      usemode=0;
      i = argc;
    }
    else if (strcmp(argv[i],"-r")==0) {
      burnin = atol(argv[i+1]);
      datafile = 2;
      i = argc;
    }
    else if (strcmp(argv[i],"+")==0) {
      ++nset;
      ++pos;
      nfiles[nset-1] = 0;
      if (nset == MAXSET) BYE("Too many different sets of data.");
    }
    else {      
      //printf(" --Bee %d %d %s\n",nset,nfiles[nset-1],argv[i]);
      files[nset-1][nfiles[nset-1]] = gzopen(argv[i],"r");
      //printf(" -00000--Bee %d %s\n",nset,argv[i]);
      if (files[nset-1][nfiles[nset-1]] == NULL) 
	printf("File %s not found---skipping.\n",argv[i]);
      else 
	++nfiles[nset-1];
      ++pos;
      if (nfiles[nset-1] == 100)
	BYE("Too many files.");
    }
  }


  int totcols[MAXSET];
  int filecols  = -1;
  long countmax = 10000;

  int new, ncols, *col[MAXSET];
  double *mean[MAXSET], *sigma[MAXSET],*l1[MAXSET],*l2[MAXSET];

  size_t *chiord;
  ncols = argc-pos-2;
  char **colnames = &argv[pos+2];
  if (ncols <= 0) BYE("Invalid number of columns %d",ncols);
  int maxtotcols = -1;

  for (i = 0; i < nset; ++i) {
    if (nfiles[i] < 1) {
      BYE("No useful files found in set %ld\n",i);
    }
    count[i] = 0;
    paramnames[i] = NULL;
    data[i] = NULL;
    filecols = -1;
    read_files(files[i],nfiles[i],&paramnames[i],&data[i],&filecols,
	       &count[i],&countmax,burnin,1,&new);
  
    totcols[i] = filecols + 1;

    //printf("%s %ld\n",argv[i+1],count[i]);
    chiord = sci_sizetvector(count[i]);
  
    gsl_sort_index(chiord,data[i][0],1,count[i]);
    chimin[i] = data[i][0][chiord[0]];
    

    col[i] = sci_ivector(totcols[i]);

    mean[i] = sci_dvector(totcols[i]);
    sigma[i] = sci_dvector(totcols[i]);
    l1[i] = sci_dvector(totcols[i]);
    l2[i] = sci_dvector(totcols[i]);

    int n, nmatch,m;
    int *matches;
    if (!withplots) {
      matches = sci_ivector(count[i]);
      if (!datafile)
	printf("%ld total records for set %ld.\n",count[i],i);
      for (n = 0; n < ncols; ++n) {
	nmatch = match_column(colnames[n],paramnames[i],totcols[i],matches);
	if (nmatch == 0)
	  BYE("Column %s not found",colnames[n]);
	for (m = 0; m < nmatch; ++m) {
	  j = matches[m];
	  get_stats(data[i][j],count[i],&mean[i][0],&sigma[i][0],&l1[i][0],
		    &l2[i][0],0.684,(usemode ? chiord[0] : -2));
	  double diff;
	  
	  /* mean[i][0] = 0.; */
	  /* for (k = 0; k < count[i] && data[i][0][chiord[k]] < chimin[i]+2; ++k) { */
	  /*   //printf("%E %ld\n",data[i][j][chiord[k]],chiord[k]); */
	  /*   if (l1[i][0] > data[i][j][chiord[k]]) { */
	  /*     l1[i][0] = data[i][j][chiord[k]]; */
	  /*   } */
	  /*   if (l2[i][0] < data[i][j][chiord[k]]) { */
	  /*     l2[i][0] = data[i][j][chiord[k]]; */
	  /*   } */
	  /*   mean[i][0] += data[i][j][chiord[k]]; */
	  /* } */
	  /* mean[i][0] /= k; */
	  /* diff = mean[i][0]-l1[i][0]; */
	  /* if (diff > sigma[i][0]) sigma[i][0] = diff; */
	  /* diff = l2[i][0]-mean[i][0]; */
	  /* if (diff > sigma[i][0]) sigma[i][0] = diff; */
	  if (datafile == 1)
	    printf("%-20s %20E %20E\n",
		   paramnames[i][j],mean[i][0],sigma[i][0]);
	  else if (datafile == 2)
	    printf("%-20s %20E %20E\n",
		   paramnames[i][j],l1[i][0],l2[i][0]);
	  else {
	    printf("Axis %d Name: %s\n",n,paramnames[i][j]);
	    printf("Axis %d Mean: %E\n",n,mean[i][0]);
	    printf("Axis %d Sigma: %E\n",n,sigma[i][0]);
	    printf("Axis %d %.2f%% range: %E - %E\n",n,100.*quantile,
		   l1[i][0],l2[i][0]);
	  }
	}
      }
      free(matches);
    }
  }
  if (!withplots) return 0;

#ifndef HAVE_LIBCPGPLOT
  BYE("Withplots mode requires libpgplot.");
#endif

  char newresponse[1000],response[1000],cmd[1000],cmd1[1000],*tok;
  char curresponse[1000];
  int gs = 4;

  int waiting=1;
  int newcmd = 1;
  int nsleep = 0;

  double maxforce[1000];
  double minforce[1000];

#pragma omp parallel num_threads(3)
{
    

  #pragma omp sections 

  {
    #pragma omp section
    {
      int mynew;
       while (1) {
      
        sleep(10);
	mynew = 0;
   
	if (waiting) {
	  for (i = 0; i < nset; ++i) {
	    filecols = totcols[i]-1;
	    read_files(files[i],nfiles[i],&paramnames[i],&data[i],&filecols,
		       &count[i],&countmax,burnin,0,&mynew);
	  }
      
	  new = mynew;
	}
       }
    }

    #pragma omp section
    {
      while (1 && !postscript) {
	while (!waiting) {
	  usleep(100000);
	}
	strcpy(newresponse,getstring("HrothgarMCMC> "));
	waiting = 0;
	newcmd = 1;
      }
    }
    
    #pragma omp section 
    {
      int n,doset[MAXSET];
      char *lab[100];

      if (postscript) {
	char *ofname;
	ofname = sci_strcat(pfile,"/CPS");
	//ofname= "/TPNG";
	//ofname = "?";
	cpgopen(ofname);
      } else
	cpgopen("/XWINDOW");
      cpgscf(2);
      cpgask(0);
      strncpy(response,colnames[0],200);
      for (i = 1; i < ncols; ++i) {
	strncat(response," ",999-strlen(response));
	strncat(response,colnames[i],999-strlen(response));
	strncat(response," ",999-strlen(response));
      }
      for (n = 0; n < nset; ++n) {
	doset[n] = 1;
      }
      do {
	strcpy(curresponse,response);
	strcpy(newresponse,response);
	cpgsch(1.5);
	if(sscanf(response,"%s",cmd) == 1 && (newcmd || new)) {
	  if (strcmp(cmd,"q")==0) {
	    cpgclos();
	    exit(0);
	  } else
	    if (!postscript) cpgeras();

	  if (newcmd && cmd[0] == '%' && strlen(cmd) > 1) {
	    gs = atoi(&cmd[1]);
	    printf("Gaussian width set to %d.\n",gs);
	  }
	  if (newcmd && cmd[0] == 'f' && strlen(cmd) > 1) 
	    doset[atoi(&cmd[1])] = 2;
	  for (n = 0; n < nset; ++n) {
	    if (doset[n] != 2) doset[n] = 1;
	  }


	  if (newcmd && strcmp(cmd,"l")==0) {
	    for (n = 0; n < nset; ++n) {
	      for (i = 0; i < totcols[n]; ++i) {
		printf("%20s",paramnames[n][i]);
		if (i % 4 ==0) printf("\n");
	      }
	    }
	    printf("\n");
	  }
	  if (strcmp(cmd,"c")==0) {
	    sscanf(response,"%s %s",cmd1,cmd);
	    contour=1;
	  } else
	    contour=0;
	  if ((col[0][0]=match_column(cmd,paramnames[0],totcols[0],NULL)) >= 0) {
	    if (new) new = 0;


	    tok = strtok(response," ");
	    if (strcmp(tok,"c") == 0)
	      tok = strtok(NULL," ");

	    ncols = 0;
	    char *l;
	    do {
	      col[0][ncols] = match_column(tok,paramnames[0],totcols[0],NULL);
	      if (col[0][ncols] < 0) {
		  printf("Column not found: %s\n",tok);
		--ncols;
	      }  else {
		for (i = 1; i < nset; ++i) {
		  col[i][ncols] = 
		    match_column(tok,paramnames[i],totcols[i],NULL);
		  if (col[i][ncols] < 0) doset[i] = 0;
		}
	      }
	      l = index(tok,'@');
	      if (l == NULL) 
		l = index(paramnames[0][col[0][ncols]],'@');

	      lab[ncols] = tok;

	      if (l != NULL) {
		if (strlen(l) > 1)
		  lab[ncols] = &l[1];
	      }

	      l = index(tok,'+');
		
	      if (l != NULL) {
		maxforce[ncols] = atof(&l[1]);
		l[0] = 0;
	      } else
		maxforce[ncols] = -1.E30;

	      l = index(tok,'-');
		
	      if (l != NULL) {
		minforce[ncols] = atof(&l[1]);
		l[0] = 0;
	      } else
		minforce[ncols] = 1.E30;


	      l = strstr(lab[ncols],"solar");
	      if (l != NULL) {
		l[0] = '\\'; l[1] = 'd'; l[2] = 9; l[3] = '\\'; l[4] = 'u';
	      }
	      while(l = index(lab[ncols],'_'))
		l[0] = ' ';

	      ++ncols;
	    }
	    while((tok = strtok(NULL," ")) != NULL);
	    
	    double l1big[8],l2big[8];
	    int nv = ncols*(ncols-1)/2;
	    if (ncols < 2 || ncols > 7)
	      printf("Inappropriate number of columns %d\n",ncols);
	    else {
	      long nx = ceil(sqrt(1.*nv));
	      long ny = ceil(1.*nv/nx);
	      
	      cpgsubp(nx,ny);
	      for (i = 0; i < ncols; ++i)  {
		l1big[i] = 1.E30;
		l2big[i] = -1.E30;
		double tmp;
		for (j = 0; j < nset; ++j) {
		  if (doset[j]) {
		    get_stats(data[j][col[j][i]],count[j],
			      &mean[j][i],&sigma[j][i],
			      &l1[j][i],&l2[j][i],
			      0.98,-1);
		    if (col[j][i] == 0) {
		      l1[j][i] = chimin[j];
		      l2[j][i] = chimin[j]+20;
		    }
		    l1[j][i] -= 0.1*fabs(l1[j][i]);
		    l2[j][i] += 0.1*fabs(l2[j][i]);
		
		  }
		  if (doset[j] == 1) {
		    l1big[i] = GSL_MIN(l1big[i],l1[j][i]);
		    l2big[i] = GSL_MAX(l2big[i],l2[j][i]);
		    if (maxforce[i] > -1.E10) l2big[i] = maxforce[i];
		    if (minforce[i] < 1.E10) l1big[i] = minforce[i];
		  }
		}

	      }
	      int gridsize=64;
	      int gridsizesq=gridsize*gridsize;

	      float *pgd = sci_fvector(gridsizesq);
	      size_t *ord = sci_sizetvector(gridsizesq);
	      float delx,dely;
	      long cell,ncell;

	      for (i = 0; i < ncols; ++i)  
		for (j = i+1; j < ncols; ++j) {
		  cpgenv(l1big[i],l2big[i],l1big[j],l2big[j],0,0);
		  for (n = 0; n < nset; ++n) {
		    if (doset[n]) {
		      if (contour) {
			delx=(l2big[i]-l1big[i])/gridsize;
			dely=(l2big[j]-l1big[j])/gridsize;
			
			int xx,yy;
			for (k = 0; k < gridsize*gridsize; ++k) 
			  pgd[k] = 0.;

			ncell=0;
			for (k = 0; k < count[n]; ++k) {
			  xx = (data[n][col[n][i]][k]-l1big[i])/delx;
			  yy = (data[n][col[n][j]][k]-l1big[j])/dely;
			  cell = xx+gridsize*yy;
			  if (cell >= 0 && cell < gridsizesq) 
			    ++ncell;
			}
			for (k = 0; k < count[n]; ++k) {
			  xx = (data[n][col[n][i]][k]-l1big[i])/delx;
			  yy = (data[n][col[n][j]][k]-l1big[j])/dely;
			  cell = xx+gridsize*yy;
			  if (cell >= 0 && cell < gridsizesq) 
			    pgd[cell]+=1./ncell;
			}
			//for (k = 0; k < gridsize*gridsize-1; ++k) 
			//printf("%5E ",pgd[k]);
			gaussian_blur(pgd,gridsize,gs);
			  
			gsl_sort_float_index(ord,pgd,1,gridsize*gridsize);
			for (k = gridsize*gridsize-2; k>=0; --k) {
			  pgd[ord[k]] += pgd[ord[k+1]];
		      //printf("%E\n",pgd[ord[k]]);
			}
			float c[] = {0.68,0.95,0.999};
			float tr[] = {l1big[i],delx,0,l1big[j],0,dely};

			int linewidth[] = {1,1,2,4,3,5};
			int linethick[] = {1,15,15,8,8};
			cpgsls(linewidth[n]);
			cpgslw(linethick[n]);
			if (n)
			  cpgcont(pgd,gridsize,gridsize,1,gridsize,
				  1,gridsize,c,-2,tr);
			else {
			  cpgsci(3);
			  cpgconf(pgd,gridsize,gridsize,1,gridsize,
				  1,gridsize,c[0],c[1],tr);
			  cpgsci(4);
			  cpgconf(pgd,gridsize,gridsize,1,gridsize,
				  1,gridsize,0.,c[0],tr);
			  cpgsci(1);
			}
			cpgsls(linewidth[0]);
			cpgslw(linethick[0]);
		    
		      }
		      else {
			float **pgd2 = sci_fmatrix(2,count[n]);
			cpgsci(1+n);
			for (k = 0; k < count[n]; ++k) {
			  pgd2[0][k] = data[n][col[n][i]][k];
			  pgd2[1][k] = data[n][col[n][j]][k];
			}
			cpgpt(count[n],pgd2[0],pgd2[1],-1);
			cpgsci(6);
			//cpgpt(20,&(pgd2[0][count[n]-21]),&(pgd2[1][count[n]-21]),2);
			cpgsci(1);
			sci_free_fmatrix(pgd2);
		      }
		    }
		  
		    if (n == 0) {
		      char cwd[1000], label[1000];
		      int beg = 0;
		      getcwd(cwd,900);
		      strcpy(label,"");
		      if (mylabel) 
		         strcpy(label,mylabel); 
		      /* else { */
		      /* 	if (!postscript)  */
		      /* 	  sprintf(label,"%s: %ld links",&cwd[beg],count[n]); */
		      /* 	else { */
		      /* 	  strcpy(label,pfile); */
		      /* 	  if (strstr(label,".eps")) */
		      /* 	    *(strstr(label,".eps")) = 0; */
		      /* 	  if (strstr(label,".ps")) */
		      /* 	  *(strstr(label,".ps")) = 0; */
		      /* 	} */
		      /* } */
		      //cpgsch(3.0);
		      if (strstr(lab[i],"M\\d200") &&
			  strstr(lab[j],"c\\d200")) {
			double z;
			if (strstr(lab[i],"/")) {
			  z = atof(strstr(lab[i],"/")+1);
			  *(strstr(lab[i],"/")) = 0;
			}
			double cnorm=5.74;
			if (strstr(lab[j],"/")) {
			  cnorm = atof(strstr(lab[j],"/")+1);
			  *(strstr(lab[j],"/")) = 0;
			}
			double delta = atof(strstr(lab[i],"d")+1);
			cnorm *= pow(1.+z,-0.47);
			float **pgd2 = sci_fmatrix(2,count[n]);
			double M200,c200,Mdelta,cdelta;
			for (k = 0; k < count[n]; ++k) {
			  cdelta = pgd2[1][k] = 0.01+30.*k/count[n];
			  c200 = get_c200_from_cdelta(cdelta,delta);
			  M200 = 0.02*pow(c200/cnorm,-10.3);
			  //if (!(k % 100))
			  //printf("delta: %f cdelta: %f M200: %f c200: %f\n",
			  //   delta,cdelta,M200,c200);
			  Mdelta = M200*delta*pow(cdelta/c200,3.)/200.;
			  //pgd2[1][k] = cnorm*pow(73*pgd2[0][k]/2,-0.084);
			  pgd2[0][k] = Mdelta;
			}
			cpgpt(count[n],pgd2[0],pgd2[1],-1);
			sci_free_fmatrix(pgd2);
		      }
		      //printf("%f\n",mytextsize);
		      cpgsch(mytextsize);
		      cpglab(lab[i],lab[j],label);
		      cpgsch(1.5);
		    }
		  }
		}
	    }
	  }
	}

	if (postscript) {
	  strcpy(newresponse,"q ");
	  newcmd = 1;
	} else
	  newcmd = 0;
	nsleep = 0;
	strcpy(response,curresponse);
	while (!newcmd && nsleep < 100) {
	  waiting = 1;
	  usleep(100000);
	  ++nsleep;
	}
	waiting = 0;
	if (newcmd) {
	  newcmd = 1;
	  strcpy(response,newresponse);
	  strcpy(curresponse,newresponse);
	} 
      }
      while(response != NULL);
    }
  }
 }
}
  
    
  /* printf("\nCorrelation coefficient: %E\n", */
  /* 	 gsl_stats_covariance(xdata,1,ydata,1,count)); */


  /* for (i = 0; i < 2; ++i) {  */
  /*   //l1[i] /= pow(1.2,GSL_SIGN(l1[i])); */
  /*   l2[i] *= pow(1.4,GSL_SIGN(l2[i])); */
  /*  } */


  /* if (argc-pos > 5) { */
    
  /*   ngrid = atol(argv[pos+3]); */
  /*   ngridsq = ngrid*ngrid; */
  /*   smooth = atol(argv[pos+4]); */
  /*   smoothsq = smooth*smooth; */
  /*   nsig = 3*smooth; */

  /*   float *prob2 = sci_fvector(ngridsq*ngridsq); */
  /*   float *prob = sci_fvector(ngridsq*ngridsq); */
  /*   double **gaussian = sci_dmatrix(nsig+1,nsig+1); */

  /*   xdel = (l2[0]-l1[0])/(ngrid-1); */
  /*   ydel = (l2[1]-l1[1])/(ngrid-1); */


  /*   for (i = 0; i < nsig; ++i) */
  /*     for (j = 0; j < nsig; ++j) */
  /* 	gaussian[i][j] = exp(-(i*i+j*j)/(2.*smoothsq))/smooth; */


  /*   for (i = 0; i < count; ++i)   */
  /*     //if (chisq[i]-chimin < 5.)  */
  /*     { */
  /* 	x0 = (xdata[i]-l1[0])/xdel; */
  /* 	y0 = (ydata[i]-l1[1])/ydel; */
  /* 	//printf("%E %E %ld %ld\n",xdata[i],ydata[i],x0,y0); */
  /* 	if (x0 >= 0 && x0 < ngrid && y0 >= 0 && y0 < ngrid) {  */
  /* 	  k = y0*ngrid+x0; */
  /* 	  prob[k] += 1.; */
  /* 	} */
  /*     } */

  /*   for (x0 = 0; x0 < ngrid; ++x0) */
  /*     for (y0 = 0; y0 < ngrid; ++y0) { */
  /* 	k0 = y0*ngrid+x0; */
  /* 	for (x = x0-nsig; x <= x0+nsig; ++x) */
  /* 	  for (y = y0-nsig; y <= y0+nsig; ++y) */
  /* 	    if (x >= 0 && x < ngrid && y >= 0 && y < ngrid) { */
  /* 	      k = y*ngrid+x; */
  /* 	      prob2[k0] += prob[k]*gaussian[abs(x-x0)][abs(y-y0)]; */
  /* 	    } */
  /*     } */


  /*   for (i = 0; i < ngridsq; ++i)   */
  /*     sum += prob2[i]; */

  /*   ord = sci_sizetvector(ngridsq); */
  /*   gsl_sort_float_index(ord,prob2,1,ngridsq); */
    
  /*   cumul = 0.; */
  /*   for (k = ngridsq-1; k >= 0; --k) { */
  /*     cumul += prob2[ord[k]]; */
  /*     prob2[ord[k]] = cumul/sum; */
  /*   } */

  /*   unlink(argv[pos+5]); */
  /*   hrothgar_writeimage(argv[pos+5],prob2,l1[0],l2[0],l1[1],l2[1],ngrid); */

  /* } */
