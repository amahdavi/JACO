#include <sciutils.h>
#include <stdio.h>
#include <string.h>
#include <stdlib.h>

#define PI 3.1415927

int main(int argc, char *argv[])
{
  float answer,angle,z;
  int   unit=1;

  z = atof(argv[4]);
  if (z > 100) z /= 299792.458;

  
  answer = ztoMpc(atof(argv[1]),atof(argv[2]),atof(argv[3]),z);

  if (argc == 5) {
    printf("%f %f\n",answer/(1.+z),answer*(1+z));
    return 0;
  }

  answer /= (1.+z);

  angle = atof(argv[5]);

  if (argc > 6) {
    if (strstr(argv[6],"arcsec") != NULL ||
	strstr(argv[6],"asec") != NULL)
      unit = 0;

    if (strstr(argv[6],"arcmin") != NULL ||
	strstr(argv[6],"amin") != NULL)
      unit = 1;

    if (strstr(argv[6],"deg") != NULL)
      unit = 2;

    if (strstr(argv[6],"rad") != NULL)
      unit = 3;

    if (strstr(argv[6],"Mpc") != NULL)
      unit = 4;

    if (strstr(argv[6],"kpc") != NULL)
      unit = 5;

  }

  switch (unit) {
    
  case 0: answer *= angle*PI/(180.*3600.);
    printf("Mpc: ");
    break;

  case 1: answer *= angle*PI/(180.*60.);
    printf("Mpc: ");
    break;

  case 2: answer *= angle*PI/(180.);
    printf("Mpc: ");
    break;

  case 3: answer *= angle;
    printf("Mpc: ");
    break;

  case 4: answer = (180.*60/PI)*angle/answer;
    printf("arcmin: ");
    break;

  case 5: answer = (180.*3600/PI)*(angle/1000.)/answer;
    printf("arcsec: ");
    break;

  }

  printf("%E\n",answer);

  return 0;
					  

}
