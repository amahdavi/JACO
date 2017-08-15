/* readdata.c
 * 
 * Copyright (C) 2007-2012 Andisheh Mahdavi
 * 
 * This program is free software; you can redistribute it and/or modify
 * it under the terms of the GNU General Public License as published by
 * the Free Software Foundation; either version 2 of the License, or (at
 * your option) any later version.
 * 
 * This program is distributed in the hope that it will be useful, but
 * WITHOUT ANY WARRANTY; without even the implied warranty of
 * MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the GNU
 * General Public License for more details.
 * 
 * You should have received a copy of the GNU General Public License
 * along with this program; if not, write to the Free Software
 * Foundation, Inc., 51 Franklin Street, Fifth Floor, 
 * Boston, MA 02110-1301, USA.
 */

#include <string.h>
#include <strings.h>
#include <stdio.h>
#include <string.h>
#include <stdlib.h>
#include "sciutils.h"
#include <stdarg.h>
#include <errno.h>
#include <limits.h>
#include <math.h>

/************************************************************************
 This procedure opens a filename pointed to by rd_filename. It reads all
 the data contained in column rd_col, and stores it in 
 a double vector. Returned are the number of records read.
 ************************************************************************/
double *double_readdata(char *rd_filename, int rd_col, 
			unsigned long *rd_count, int offset)
{
  unsigned long   i, count;
  int j, comment;
  FILE  *datafile;
  double *rd_data;
  char  line[10000];

  datafile = fopen(rd_filename,"r");
  if (datafile == NULL) {
    printf("Error opening file %s.\n",rd_filename);
    exit(1);
  }

  if (rd_col < 1) {
    printf("Error -- illegal column number %d specified.\n",rd_col);
    exit(1);
  }

  count = 0;
  while(!feof(datafile)) {    
    ++count;
    fgets(line,9999,datafile);
    if (strlen(line) == 9999) {
      printf("Error---line length too long in %s.\n",rd_filename);
      exit(1);
    }
  }
  --count;

  fclose(datafile);
  datafile = fopen(rd_filename,"r");
  rd_data  = (double *)malloc((offset+count)*sizeof(double));

  long k = offset;
  for (i = offset; i < count+offset; ++i) {
    comment = 0;
    for (j = 1; j < rd_col && !comment; ++j) {
      fscanf(datafile,"%s",line); 
      if (line[0] == '#') comment = 1;
    }
    fscanf(datafile,"%lf",&rd_data[k]);
    if (!comment) ++k;
    fgets(line,9999,datafile);
  }
  
  fclose(datafile);
  *rd_count = k-offset;
  return(rd_data);
}

float *float_readdata(char *rd_filename, int rd_col, 
			unsigned long *rd_count, int offset)
{
  unsigned long   i, count;
  int j, comment;
  FILE  *datafile;
  float *rd_data;
  char  line[10000];

  datafile = fopen(rd_filename,"r");
  if (datafile == NULL) {
    printf("Error opening file %s.\n",rd_filename);
    exit(1);
  }

  if (rd_col < 1) {
    printf("Error -- illegal column number %d specified.\n",rd_col);
    exit(1);
  }

  count = 0;
  while(!feof(datafile)) {    
    ++count;
    fgets(line,9999,datafile);
    if (strlen(line) == 9999) {
      printf("Error---line length too long in %s.\n",rd_filename);
      exit(1);
    }
  }
  --count;

  fclose(datafile);
  datafile = fopen(rd_filename,"r");
  rd_data  = (float *)malloc((offset+count)*sizeof(float));

  long k = offset;
  for (i = offset; i < count+offset; ++i) {
    comment = 0;
    for (j = 1; j < rd_col && !comment; ++j) {
      fscanf(datafile,"%s",line); 
      if (line[0] == '#') comment = 1;
    }
    fscanf(datafile,"%f",&rd_data[k]);
    if (!comment) ++k;
    fgets(line,9999,datafile);
  }
  
  fclose(datafile);
  *rd_count = k-offset;
  return(rd_data);
}

// Add a string to an existing char ** array
void grow_string_array(char ***ins, long nchar, char *string) {

  char **outs;
  long i;

  outs = (char **)malloc((nchar+1)*sizeof(char *));
  for (i = 0; i < nchar; ++i) {
    outs[i] = sci_strdup((*ins)[i]);
    free((*ins)[i]);
  }
  outs[nchar] = sci_strdup(string);
  *ins = outs;
}

// Parse a comma separated list, allocating a char ** to point to
// its individual elements, returning the number of elements.
// The original list is unaltered. 
int parse_list(const char *inputlist, char delimit, char ***parsedlist)
{
  int nitems=0,noncomma=0,firstitem=-1;
  char *clist;
  size_t len,i;

  *parsedlist = NULL;
  clist = sci_strdup(inputlist);

  len = strlen(clist);
  if (len > INT_MAX || !len) return 0;

  for (i = 0; i < len; ++i) 
    if (clist[i] == delimit || clist[i] == 0) {
      clist[i] = 0;
      if (noncomma) {
	noncomma = 0;
	++nitems;
      } 
    } else {
      if (firstitem < 0) firstitem = i;
      ++noncomma;
    }

  if (noncomma) ++nitems;

  if (!nitems) return 0;
      
  *parsedlist = (char **)malloc(nitems*sizeof(char *));

  (*parsedlist)[0] = &clist[firstitem];
  nitems = 1;
  for (i = firstitem; i < len-1; ++i) {
    if (clist[i] == 0 && clist[i+1] > 0) 
      (*parsedlist)[nitems++] = &clist[i+1];
  }

  return nitems;
}

int parse_list_double(const char *inputlist, char delimit, double **parsedlist)
{
  int nitems=0,noncomma=0,firstitem=-1;
  char *clist;
  size_t len,i;

  *parsedlist = NULL;
  clist = sci_strdup(inputlist);

  len = strlen(clist);
  if (len > INT_MAX || !len) return 0;

  for (i = 0; i < len; ++i) 
    if (clist[i] == delimit || clist[i] == 0) {
      clist[i] = 0;
      if (noncomma) {
	noncomma = 0;
	++nitems;
      } 
    } else {
      if (firstitem < 0) firstitem = i;
      ++noncomma;
    }

  if (noncomma) ++nitems;

  if (!nitems) return 0;
      
  *parsedlist = (double *)malloc(nitems*sizeof(double));

  (*parsedlist)[0] = atof(&clist[firstitem]);
  nitems = 1;
  for (i = firstitem; i < len-1; ++i) {
    if (clist[i] == 0 && clist[i+1] > 0) 
      (*parsedlist)[nitems++] = atof(&clist[i+1]);
  }

  return nitems;
}

			     
// Read multiple ascii files using a scanf-like format.
long read_ascii(char *filename, 
		const char *fmt, char ***comments, ...)
{
  char buffer[5001], *buff, tmpstr[5001], ***str_arr, *maxpos, *endptr;
  double **double_arr;
  char **filenames;
  int **int_arr,result,nf,j;
  long **long_arr,i,longval;
  const char *format;
  long count;
  va_list args;
  FILE *datafile;
  
  nf = parse_list(filename,',',&filenames);

  count = 0;
  for (i = 0; i < nf; ++i) {

    datafile = fopen(filenames[i],"r");
    if (datafile == NULL) {
      printf("read_ascii: Error opening file %s.\n",filenames[i]);
      exit(1);
    }

    // Count the number of lines
    while (fgets(buffer,5000,datafile) != NULL) {
      if (strlen(buffer) == 5000) {
	printf("read_ascii: line too long in %s\n",filenames[i]);
	return -1;
      }
      result = sscanf(buffer,"%s",tmpstr);
      if ((comments || (tmpstr[0] != '#')) && result > 0)
	++count;
    }
    fclose(datafile);
  }

  if (!count) return count;

  // Allocate the various vectors
  va_start(args, comments);
  format = fmt;

  if (comments)
    *comments = (char **)malloc(count*sizeof(char *));

  while (*format) {
    switch(*format++) {
      
    case 's':
      str_arr = va_arg(args, char ***);
      *str_arr = (char **)malloc(count*sizeof(char *));
      break;

    case 'l':
      if (*format == 'd') {
	long_arr = va_arg(args, long **);
	*long_arr = (long *)malloc(count*sizeof(long));
	format++;
      }
      break;

    case 'd':
      int_arr = va_arg(args, int **);
      *int_arr = (int *)malloc(count*sizeof(int));
      break;

    case 'f':
      double_arr = va_arg(args, double **);
      *double_arr = (double *)malloc(count*sizeof(double));
      break;

    }
  }
  va_end(args);


  i = 0;
  for (j = 0; j < nf; ++j) {
    datafile = fopen(filenames[j],"r");
    
    // Reread the file
    while (fgets(buffer,5000,datafile) != NULL) {
      
      maxpos = index(buffer,'\n');
      if (maxpos) *maxpos = 0;
      
      result = sscanf(buffer,"%s",tmpstr);
      if ((comments || tmpstr[0] != '#') && result > 0) {
	maxpos = index(buffer,'#');
	if (comments) (*comments)[i] = NULL;
	if (maxpos) { 
	  if (comments) (*comments)[i] = sci_strdup(maxpos);
	  *maxpos = 0;
	}
	
	buff = buffer;
	format = fmt;
	va_start(args, comments);
      
	while (*format) {
	  errno = 0;
	
	  switch(*format++) {
	  
	  case 's':
	    str_arr = va_arg(args, char ***);
	    if (sscanf(buff,"%s",tmpstr) <= 0)
	      (*str_arr)[i] = NULL;
	    else {
	      buff = strstr(buff,tmpstr)+strlen(tmpstr);
	      (*str_arr)[i] = sci_strdup(tmpstr);
	    }
	  
	    break;
	  
	  case 'l':
	    if (*format == 'd') {
	      long_arr = va_arg(args, long **);
	      if (sscanf(buff,"%s",tmpstr) <= 0) 
		(*long_arr)[i] = -9999;
	      else {
		buff = strstr(buff,tmpstr)+strlen(tmpstr);
		(*long_arr)[i] = strtol(tmpstr,&endptr,10);
		if (errno == ERANGE || endptr[0]) {
		  printf("read_ascii: Error converting %s to long\n",tmpstr);
		  return -1;
		}
	      }
	      format++;
	    }
	    break;
	  
	  case 'f':
	    double_arr = va_arg(args, double **);
	    if (sscanf(buff,"%s",tmpstr) <= 0) 
	      (*double_arr)[i] = nan("F");
	    else {
	      buff = strstr(buff,tmpstr)+strlen(tmpstr);
	      (*double_arr)[i] = strtod(tmpstr,&endptr);
	      if (errno == ERANGE || endptr[0]) {
		printf("read_ascii: Error converting %s to double.\n",tmpstr);
		return -1;
	      }
	    }
	    break;
	
	  case 'd':
	    int_arr = va_arg(args, int **);
	    if (sscanf(buff,"%s",tmpstr) <= 0) 
	      (*int_arr)[i] = -9999;
	    else {
	      buff = strstr(buff,tmpstr)+strlen(tmpstr);
	      longval = strtol(tmpstr,&endptr,10);
	      if (errno == ERANGE || endptr[0] || 
		  longval < INT_MIN || longval > INT_MAX) {
		printf("read_ascii: Error converting %s to int.\n",tmpstr);
		return -1;
	      }
	      (*int_arr)[i] = (int)longval;
	    }
	    break;
	  }
	}
	va_end(args);
      } else i--;
      i++;
    }
    fclose(datafile);
  }
  return count;
}

/* int main(int argc, char *argv[]) */
/* { */
/*   long count,i; */
/*   char **par,**comments; */
/*   long *value; */

/*   count = read_ascii(argv[1],"%s %ld",&comments,&par,&value); */

/*   for (i = 0; i < count; ++i) { */
/*     if (par[i])  */
/*       printf("%s %ld ",par[i],value[i]); */
/*     if (comments[i]) */
/*       printf("%s",comments[i]); */
/*     printf("\n"); */
/*   } */
    
/*   return 0; */
/* } */

/* int main(int argc, char *argv[]) */
/* { */
/*   int z; */
/*   char **mylist; */

/*   z = parse_list(argv[1],'-',&mylist); */
/*   printf("*%d*\n",z); */
/*   while (z > 0) */
/*     printf("%s\n",mylist[--z]); */

/* } */
