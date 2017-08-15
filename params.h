static char *detector[] __attribute__ ((unused)) = {"MOS 1","MOS 2","pn","ACIS-I","ACIS-S",
                           "SZ","WL"};
#ifdef NONGRAVITY

#define NTOTPARAMS 328

#else

#define NTOTPARAMS 60

#endif


#define XRAYRA               0 
#define XRAYDEC		     1
#define MG1                  2
#define MG2                  3
#define MG3                  4
#define ALPHA                5
#define BETA1                6
#define BETA2                7
#define BETA3                8
#define RX1                  9
#define RX2                 10
#define RX3                 11
#define AXIS1               12
#define AXIS2               13
#define PHI                 14
#define THETA               15
#define PSI                 16
#define DARKMASS            17
#define DARKSLOPE           18
#define RDM0                19
#define RDM1                20
#define Z0                  21
#define ZINF                22
#define RZ                  23
#define STARMASS            24
#define STARSLOPE           25
#define RSTAR               26
#define ANISO0              27
#define ANISO1              28
#define RANISO              29
#define RSHOCK              30
#define SHOCKQ              31
#define REDSHIFT            32
#define CONTRAST            33
#define NONTHERMAL          34
#define TEMPBIAS            35
#define EXT                 36
#define EXZ                 37
#define EXNORM              38
#define W0                  39
#define WR                  40
#define PLSLOPE0            41
#define PLSLOPE1            42
#define PLSLOPE2            43
#define PLSLOPE3            44
#define PLSLOPE4            45
#define PLNORM0             46
#define PLNORM1             47
#define PLNORM2             48
#define PLNORM3             49
#define PLNORM4             50
#define BACKT               51
#define BACKZ               52
#define BNORM0              53
#define BNORM1              54
#define BNORM2              55
#define BNORM3              56
#define BNORM4              57
#define XIMNORM             58
#define XIMBACK             59

static char *paramnames[] __attribute__ ((unused)) = {

  /* 000 */     "xrayra",
  /* 001 */     "xraydec",
 /* 002 */     "Mg1",		 
 /* 003 */     "Mg2",		 
 /* 004 */     "Mg3",		 
 /* 005 */     "alpha", 		 
 /* 006 */     "b1",    		 
 /* 007 */     "b2",    		 
 /* 008 */     "b3",    		 
 /* 009 */     "rx1",		 
 /* 010 */     "rx2",		 
 /* 011 */     "rx3",		 
 /* 012 */     "axis1",		 
 /* 013 */     "axis2",		 
 /* 014 */     "phi",		 
 /* 015 */     "theta",		 
 /* 016 */     "psi",		 
 /* 017 */     "Md",	   	 
 /* 018 */     "ndark",	   	 
 /* 019 */     "rdm0",	   	 
 /* 020 */     "rdm1",	   	 
 /* 021 */     "z0",	   	 
 /* 022 */     "zinf",	   	 
 /* 023 */     "rz",	   	 
 /* 024 */     "Mstar",	   	 
 /* 025 */     "nstar",	   	 
 /* 026 */     "rstar",	   	 
 /* 027 */     "aniso0",	   	 
 /* 028 */     "aniso1",	   	 
 /* 029 */     "raniso",	   	 
 /* 030 */     "rtrunc",	   	 
 /* 031 */     "Ttrunc",	   	 
 /* 032 */     "redshift",   	 
 /* 033 */     "Contrast",   	 
 /* 034 */     "nonthermal",  	 
 /* 035 */     "bias",  		 
 /* 036 */     "exT",	   	 
 /* 037 */     "exZ",	   	 
 /* 038 */     "exnorm",	   	 
 /* 039 */     "nH0",	   	 
 /* 040 */     "dnHdr",	   	 
 /* 041 */     "plbg0",	   	 
 /* 042 */     "plbg1",	   	 
 /* 043 */     "plbg2",	   	 
 /* 044 */     "plbg3",	   	 
 /* 045 */     "plbg4",	   	 
 /* 046 */     "plbnorm0",	 
 /* 047 */     "plbnorm1",	 
 /* 048 */     "plbnorm2",	 
 /* 049 */     "plbnorm3",	 
 /* 050 */     "plbnorm4",	 
 /* 051 */     "sxbT",	   	 
 /* 052 */     "sxbZ",	   	 
 /* 053 */     "bnorm0",	   	 
 /* 054 */     "bnorm1",	   	 
 /* 055 */     "bnorm2",	   	 
 /* 056 */     "bnorm3",	   	 
 /* 057 */     "bnorm4",
 /* 058 */     "ximnorm",          
 /* 059 */     "ximback"            

};

static int dofit[] __attribute__ ((unused)) = {
/* 000 */   1,
/* 001 */   1,
/* 002 */   0,
/* 003 */   1,
/* 004 */   1,
/* 005 */   0,
/* 006 */   0,
/* 007 */   1,
/* 008 */   1,
/* 009 */   0,
/* 010 */   1,
/* 011 */   1,
/* 012 */   1,
/* 013 */   1,
/* 014 */   1,
/* 015 */   1,
/* 016 */   1,
/* 017 */   0,
/* 018 */   1,
/* 019 */   1,
/* 020 */   0,
/* 021 */   1,
/* 022 */   1,
/* 023 */   1,
/* 024 */   1,
/* 025 */   1,
/* 026 */   1,
/* 027 */   1,
/* 028 */   1,
/* 029 */   1,
/* 030 */   1,
/* 031 */   1,
/* 032 */   1,
/* 033 */   1,
/* 034 */   1,
/* 035 */   1,
/* 036 */   1,
/* 037 */   1,
/* 038 */   1,
/* 039 */   1,
/* 040 */   1,
/* 041 */   1,
/* 042 */   1,
/* 043 */   1,
/* 044 */   1,
/* 045 */   1,
/* 046 */   1,
/* 047 */   1,
/* 048 */   1,
/* 049 */   1,
/* 050 */   1,
/* 051 */   1,
/* 052 */   1,
/* 053 */   1,
/* 054 */   1,
/* 055 */   1,
/* 056 */   1,
/* 057 */   1,
/* 058 */   1,
/* 050 */   1
};

static char *paramunits[] __attribute__ ((unused)) = {
/* 000 */   "Right ascension of X-ray center",
/* 001 */   "Declination of X-ray center",
/* 002 */   "TOTAL gas fraction/gas mass at rContrast (see fitfgas)",		      
/* 003 */   "2nd beta model fractional contribution to gas mass",		      
/* 004 */   "3rd beta model fractional contribution to gas mass",		      
/* 005 */   "Multiply the first beta model by r^{-alpha}", 			      
/* 006 */   "Beta-parameter for 1st beta model",	    			      
/* 007 */   "Beta-parameter for 2nd beta model",	    			      
/* 008 */   "Beta-parameter for 3rd beta model",           			      
/* 009 */   "1st beta core radius in Mpc",  					      
/* 010 */   "2nd beta core radius in Mpc",  					      
/* 011 */   "third beta core radius in Mpc",					      
/* 012 */   "Ratio of X' axis to Z' axis",					      
/* 013 */   "Ratio of Y' axis to Z' axis",					      
/* 014 */   "First Euler angle (plane of sky)",					      
/* 015 */   "Second Euler angle (line of sight)",				      	
/* 016 */   "Third Euler angle (Z' axis)",					      	
/* 017 */   "TOTAL (gasmodel<0) or DARK (gasmodel>0) mass @rContrast, 1e14 Msun",    	
/* 018 */   "Inner slope of Generalized NFW (|darkmodel|=5, else see docs)",	      	
/* 019 */   "Not used for Generalized NFW; see docs",				      	
/* 020 */   "Scale radius (darkmodel=5) or concentration (darkmodel=-5), see docs ",  	
/* 021 */   "Inner metalicity/Zsun: z(r) = z0*(1+r*r/rz*rz)^(-3*zinf/2)",	      	
/* 022 */   "Metallicity slope: z(r) = z0*(1+r*r/rz*rz)^(-3*zinf/2)",		      	
/* 023 */   "Metallicity core radius in z(r) = z0*(1+r*r/rz*rz)^(-3*zinf/2)",	      
/* 024 */   "Total stellar mass at rContrast, 1e14 Msun",			      	
/* 025 */   "Einasto alpha parameter for stellar mass profile",			      	
/* 026 */   "Einasto r^-2 for stellar mass profile",				      	
/* 027 */   "Central velocity dispersion anisotropy",				      	
/* 028 */   "Outer velocity dispersion anisotropy",				       	
/* 029 */   "Anisotropy transition radius, Mpc",				      	
/* 030 */   "Gas distribution trunction radius in units of overdensity",	      
/* 031 */   "Temperature at the gas truncation surface (keV)",	       		       
/* 032 */   "Redshift of the cluster",	       					      
/* 033 */   "Overdensity for Mg1, Md, and Mstar",                                     
/* 034 */   "fractional nonthermal pressure support  at 1 Mpc",	       		      
/* 035 */   "Multiply Chandra Effective area by E^bias",	       		      
/* 036 */   "Nuisance central power law source: slope",	       			      
/* 037 */   "Nuisance central power law source: metallicity",	       		      
/* 038 */   "Spectrum normalization of central nuisance source",  		      
/* 039 */   "Galactic absorbption in 10^22 cm^-2",	       			      
/* 040 */   "Radial gradient in absorption, 10^22 cm^-2/arcmin",   		      
/* 041 */   "Slope of residual power law particle background for MOS1",	  	      
/* 042 */   "Slope of residual power law particle background for MOS2 ",  	      
/* 043 */   "Slope of residual power law particle background for pn ",	  	      
/* 044 */   "Slope of residual power law particle background for ACIS-I ",	      
/* 045 */   "Slope of residual power law particle background for ACIS-S ",	      	       
/* 046 */   "Normalization of residual power law particle background for MOS1",	       	       
/* 047 */   "Normalization of residual power law particle background for MOS2 ",      	       
/* 048 */   "Normalization of residual power law particle background for pn ",	      	       
/* 049 */   "Normalization of residual power law particle background for ACIS-I ",    
/* 050 */   "Normalization of residual power law particle background for ACIS-S ",    
/* 051 */   "Temperature of soft X-ray background",	       			      
/* 052 */   "Metallicity of soft X-ray background",	       			      
/* 053 */   "Normalization of residual soft X-ray background for MOS1",	  	      
/* 054 */   "Normalization of residual soft X-ray background for MOS2 ",  	      
/* 055 */   "Noramlization of residual soft X-ray background for pn ",	  	      
/* 056 */   "Normalization of residual soft X-ray background for ACIS-I ",	      
/* 057 */   "Normalization of residual soft X-ray background for ACIS-S ",            
/* 058 */   "Normalization of soft X-ray image ",            
/* 059 */   "Background in soft X-ray image "            
};

static double initvalues[] __attribute__ ((unused)) = {
/* 000 */   0., 
/* 001 */   0.,
/* 002 */   0.1,      
/* 003 */   0.0,      
/* 004 */   0.0,      
/* 005 */   0.0,      
/* 006 */   0.6,      
/* 007 */   0.6,      
/* 008 */   0.6,      
/* 009 */   0.3,      
/* 010 */   0.3,      
/* 011 */   0.3,      
/* 012 */   1.0,      
/* 013 */   1.0,      
/* 014 */   0.,	      
/* 015 */   0.,	      
/* 016 */   0.,	      
/* 017 */   1.,	      
/* 018 */   1.,	      
/* 019 */   0.,	      
/* 020 */   0.3,      
/* 021 */   0.3,      
/* 022 */   0.65,     
/* 023 */   0.3,      
/* 024 */   0.,	      
/* 025 */   0.25,     
/* 026 */   1.e-6,    	      
/* 027 */   0.,	      	      
/* 028 */   0.3,      	      
/* 029 */   0.3,      
/* 030 */   100.,     	      
/* 031 */   3.,	      
/* 032 */   0.1,      
/* 033 */   500.,     
/* 034 */   0.,	      
/* 035 */   0.,	      
/* 036 */   4.,	      
/* 037 */   1.,	      	      
/* 038 */   0.,	      
/* 039 */   0.05,     	      
/* 040 */   0.,	      	      
/* 041 */   1.4,      	      
/* 042 */   1.4,      	      
/* 043 */   1.4,      	      
/* 044 */   1.4,      
/* 045 */   1.4,      
/* 046 */   0 ,	      
/* 047 */   0 ,	      
/* 048 */   0 ,	      
/* 049 */   0 ,	      	      
/* 050 */   0 ,	      	      
/* 051 */   0.25,     
/* 052 */   1.0,      
/* 053 */   0,	      
/* 054 */   0.,	      
/* 055 */   0.,	      
/* 056 */   0.,	      
/* 057 */   0.,	                      
/* 058 */   1.,	                      
/* 059 */   0.1,	                      
};


static double parmin[] __attribute__ ((unused)) = {
  /* 000 */   0.,
  /* 001 */   -90.,
  /* 002 */   0.0000,	   
  /* 003 */   0.,	   
  /* 004 */   0.,	   
  /* 005 */   0.,	   
  /* 006 */   0.4,	   
  /* 007 */   0.4,	   
  /* 008 */   0.4,	   
  /* 009 */   0.0005,	   
  /* 010 */   0.0005,	   
  /* 011 */   0.0005,	   
  /* 012 */   0.1,	   
  /* 013 */   0.1,	   
  /* 014 */   0.0,	   
  /* 015 */   0.0,	   
  /* 016 */   0.0,	   
  /* 017 */   0.05,	   
  /* 018 */   0.,	   
  /* 019 */   0.,	   
  /* 020 */   0.05,	   
  /* 021 */   0.1,	   
  /* 022 */   0.0,	   
  /* 023 */   0.005,	   
  /* 024 */   0.,	   
  /* 025 */   0.15,	   
  /* 026 */   1e-8,    	   
  /* 027 */   -3.,     	   
  /* 028 */   -3.,     	   
  /* 029 */   0.,      	   
  /* 030 */   0.,      	   
  /* 031 */   0.5,     	   
  /* 032 */   0.,      	   
  /* 033 */   -100.,   	   
  /* 034 */   0.,   	   
  /* 035 */   0.,     	   
  /* 036 */   0.1,    	   
  /* 037 */   0.1,    	   
  /* 038 */   0.,     	   
  /* 039 */   1.e-3,  	   
  /* 040 */   -1.,    	   
  /* 041 */   -0.5,   	   
  /* 042 */   -0.5,   	   
  /* 043 */   -0.5,   	   
  /* 044 */   -0.5,   	   
  /* 045 */   -0.5,   	   
  /* 046 */   -1.e-2,	   
  /* 047 */   -1.e-2,	   
  /* 048 */   -1.e-2,	   
  /* 049 */   -1.e-2,	   
  /* 050 */   -1.e-2,	   
  /* 051 */   0.1,    	   
  /* 052 */   0.1,    	   
  /* 053 */   -1.e-2, 	   
  /* 054 */   -1.e-2, 	   
  /* 055 */   -1.e-2, 	   
  /* 056 */   -1.e-2, 	   
  /* 057 */   -1.e-2,      
  /* 058 */   0,      
  /* 059 */   0     
};	                   
	                   
static double parmax[] __attribute__ ((unused)) = {
  /* 000 */   360.,
  /* 001 */   90.,
/* 002 */   10.,      
/* 003 */   1.0,      
/* 004 */   1.0,      
  /* 005 */   1.5,      
/* 006 */   5.0,      
/* 007 */   5.0,      
/* 008 */   5.0,      
/* 009 */   5.0,      
/* 010 */   5.0,      
/* 011 */   5.0,      
/* 012 */   1.0,      
/* 013 */   1.0,      
/* 014 */   360.0,    
/* 015 */   180.0,    
/* 016 */   360.0,    
/* 017 */   100.,     
/* 018 */   5.0,      
/* 019 */   .1,	      
/* 020 */   25.,      
/* 021 */   2.9,      
/* 022 */   2.0,      
/* 023 */   1.0,      
/* 024 */   0.1,      
/* 025 */   1.0,      
/* 026 */   0.3,      
/* 027 */   0.999,    
/* 028 */   0.999,    
/* 029 */   2.,	      
/* 030 */   5000.,    
/* 031 */   8.,       
/* 032 */   5.,       
/* 033 */   20000.,   
/* 034 */   0.7,      
/* 035 */   1.2,      
/* 036 */   15.,      
/* 037 */   2.7,      
/* 038 */   1.,       
/* 039 */   10.,      
/* 040 */   1.,       
/* 041 */   2.0,      
/* 042 */   2.0,      
/* 043 */   2.0,      
/* 044 */   2.0,      
/* 045 */   2.0,      
/* 046 */   1.e-2,    
/* 047 */   1.e-2,    
/* 048 */   1.e-2,    
/* 049 */   1.e-2,    
/* 050 */   1.e-2,    
/* 051 */   0.8,      
/* 052 */   2.7,      
/* 053 */   1.e-2,    
/* 054 */   1.e-2,    
/* 055 */   1.e-2,    
/* 056 */   1.e-2,    
/* 057 */   1.e-2,    
/* 058 */   1000.,      
/* 059 */   15.
};

#define NSTRINGPARAMS 40

#define DEBUG 0
#define SPECTRALCODE  1
#define E1 2
#define E2 3
#define SYSERR 4
#define WLMODE 5
#define WLDATA 6
#define WLCONFIG 7
#define SZMODE 8
#define SZDATA 9
#define SZCONFIGS 10
#define VELMODE 11
#define VELDATA 12
#define XRAYMODE 13
#define XRAYDATA 14
#define XRAYPSF 15
#define RCUT 16
#define GASMODEL 17
#define DARKMODEL 18
#define TPROJPOWER 19
#define HUBBLE 20
#define OMEGA_M 21
#define OMEGA_L 22
#define WRITE 23
#define DOSIM 24
#define SIMFILE 25
#define AXIS 26
#define SIMMASS 27
#define SIMLENGTH 28
#define XIMAGE 29
#define XIMAGEE1 30
#define XIMAGEE2 31
#define DELTAS 32
#define BACKLIST 33
#define FITMC 34
#define FITFGAS 35
#define LXE1 36
#define LXE2 37
#define PRECISION 38
#define FINEGRID 39
  
#ifdef JACO_INIT

char *stringparname[NSTRINGPARAMS] = {
    /* 000 */    "debug",
    /* 001 */    "spectralcode",
    /* 002 */    "e1",
    /* 003 */    "e2",
    /* 004 */    "xsyserr",
    /* 005 */    "wlmode",
    /* 006 */    "wldata",
    /* 007 */    "wlbeta",
    /* 008 */    "szmode",
    /* 009 */    "szdata",
    /* 010 */    "szconf",
    /* 011 */    "velmode",
    /* 012 */    "veldata",
    /* 013 */    "xraymode",   
    /* 014 */    "xraydata",   
    /* 015 */    "psffiles",   
    /* 016 */    "Rcut",       
    /* 017 */    "gasmodel",     
    /* 018 */    "darkmodel",    
    /* 019 */    "tprojpow",   
    /* 020 */    "hubble",     
    /* 021 */    "omegam",     
    /* 022 */    "omegal",     
    /* 023 */    "write",      
    /* 024 */    "dosim",      
    /* 025 */    "simfile",    
    /* 026 */    "axis",       
    /* 027 */    "simmass",    
    /* 028 */	 "simlength",  
    /* 029 */	 "ximage",
    /* 030 */	 "xime1",
    /* 031 */	 "xime2",
    /* 032 */    "deltas",
    /* 033 */    "backsyserr",
    /* 034 */    "fitmc",
    /* 035 */    "fitfgas",
    /* 036 */    "lxe1",
    /* 037 */    "lxe2",
    /* 038 */    "precision",
    /* 039 */    "finegrid"
};  
  
char *stringparval[NSTRINGPARAMS] = {
    /* 000 */    "0",
    /* 001 */    "spectralcode.txt",
    /* 002 */    "0.6",
    /* 003 */    "6",
    /* 004 */    "0.04",
    /* 005 */    "0",
    /* 006 */    "wl.dat",
    /* 007 */    "1",
    /* 008 */    "0",
    /* 009 */    "sz.dat",
    /* 010 */    "cbi.conf",
    /* 011 */    "0",
    /* 012 */    "veldisp.dat",
    /* 013 */    "1",
    /* 014 */    "jaco-xray.dat",		    
    /* 015 */    "m1psf,m2psf,pnpsf,aspsf,aipsf",   
    /* 016 */    "0",				    
    /* 017 */    "-2",				    
    /* 018 */    "-5",				    
    /* 019 */    "0.25",			    
    /* 020 */    "70.",				    
    /* 021 */    "0.3",				    
    /* 022 */    "0.5",				    
    /* 023 */    "0",				    
    /* 024 */    "0",				    
    /* 025 */    "simul.dat",			    
    /* 026 */    "1",				    
    /* 027 */    "1",				    
    /* 028 */	 "1",				    
    /* 029 */	 "NONE",
    /* 030 */	 "0.2",
    /* 031 */	 "0.7",
    /* 032 */    "2500.,1000.,500.,200.",
    /* 033 */    "1.e-4,1.e-4,1.e-4,7.e-6,7.e-6",
    /* 034 */    "0",
    /* 035 */    "1",
    /* 036 */    "0.01",
    /* 037 */    "100",
    /* 038 */    "1.e-4",
    /* 049 */    "0",
};                           
  
char *stringparcomments[NSTRINGPARAMS] = {
    /* 000 */    "Debug level (0-10)",
    /* 001 */    "Spectral code file",
    /* 002 */    "Minimum fitting energy",
    /* 003 */    "Maximum fitting energy",
    /* 004 */    "X-ray fractional systematic error",
    /* 005 */    "Include weak lensing? 0/1",
    /* 006 */    "Lensing data (r1 r2 shear err)",
    /* 007 */    "Lensing beta, beta squared",
    /* 008 */    "Include SZ? 0/1",
    /* 009 */    "Intermediate SZ observables file",
    /* 010 */    "SZ configuration files (comma-separated)",
    /* 011 */    "Include velocity dispersions? 0/1",
    /* 012 */    "Velocity dispersion data file",
    /* 013 */    "0=do not fit X-ray data; 1=fit X-ray data",	 
    /* 014 */    "X-ray data files (comma-separated)",		 
    /* 015 */    "List of X-ray psf files (comma-separated)",	 
    /* 016 */    "X-ray cutoff radius",				 
    /* 017 */    "+/-1=triaxial;+/-2=sphere; sets Md as total(+) or dark(-) mass"
    /* 018 */    "+/-5=GNFW profile; see docs; sign sets rdm1 behavior",	
    /* 019 */    "Temp. projection power (1=emission-weight)",	 
    /* 020 */    "Hubble constant in in km/s/Mpc",		 
    /* 021 */    "Omega Matter",				 
    /* 022 */    "Omega Lambda",				 
    /* 023 */    "1=best-fit spectra; 2=profiles",		 
    /* 024 */    "0=no; 1=model; 2=simulations; 3=with image",	 
    /* 025 */    "simulation input file",			 
    /* 026 */    "simulation axis",				 
    /* 027 */    "simulation mass unit / 1e14 Msun",		 
    /* 028 */	 "simulation length unit / Mpc",		 
    /* 029 */	 "X-ray image+optional exposure mask (comma-separated)",
    /* 030 */    "Lower energy of X-ray image",
    /* 031 */    "Upper energy of X-ray image",
    /* 032 */    "deltas for profile output",
    /* 033 */    "List of residual background systematic errors",
    /* 034 */    "1=Use Duffy et al. mass-concentration relation",
    /* 035 */    "0: Mg1=total gas mass; 1: Mg1=total gas fraction",
    /* 036 */    "Lower energy of restframe band to report Lx in",
    /* 037 */    "Upper energy of restframe band to report Lx in",
    /* 038 */    "Desired precision",
    /* 039 */    "Use superfine pressure grid (slower)? (0/1)"
};
#endif

