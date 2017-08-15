#include <jaco.h>
#include <standalone.h>
#include <string.h>
#include <strings.h>
#include <mpi.h>

int main(int argc, char *argv[])
{
  struct fitdata data;
  struct jaco_state js;
  struct hrothgar_setup setup;
  char cmd[1000],*dotpos;

  data.js = &js;

  jaco_init(&js,&setup,argc,argv);
  js.standalone = 1;
  data.fitannuli=FALSE;

  // Read in ARF, RMF, and PHA files, as well as the parameter file
  init_standalone(&data);
  // printf("%f\n",phalist[0].flux[phalist[0].nbins-1]);

  if (js.xraycount > 0)
    data.spectrum = sci_dvector(data.arflist[0].nebins);

  js.profname = sci_strdup((setup.outfilename == NULL ||
			    strcmp(setup.infilename,"/dev/null") == 0) ? 
			   setup.infilename :( 
			   (strcmp(setup.outfilename,"/dev/null") ? 
			    setup.outfilename : "jaco")));

  dotpos = rindex(js.profname,'.');
  if (dotpos != NULL)
    dotpos[0] = 0;
  
  js.rng = setup.generator;


  if (js.fitwrite > 0 && js.mpirank == 0) {
    if (setup.overwrite)
      snprintf(cmd,999,"find . -name  'bin*fits.dat' -exec rm '{}' ';'");
    else
      snprintf(cmd,999,"find . -name  'bin*fits.dat' -ok rm '{}' ';'");
    system(cmd);
    if (js.fitwrite > 2) {
      if (setup.overwrite)
	snprintf(cmd,999,"find . -name '%s-[pc]*' -exec rm '{}' ';'",
		 js.profname);
      else
	snprintf(cmd,999,"find . -name '%s-[pc]*' -ok rm '{}' ';'",
		 js.profname);
      system(cmd);
    }
  }
 

  hrothgar(data.ndata,NULL,data.bigdata,data.bigerr,
	   get_bigmodel,&data,&setup);

  return 0;

}
  
