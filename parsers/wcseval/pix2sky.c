/* Convert between pixel and celestial coordinates */

/* Syntax: pix2sky [options] fits_file [x y | infile] */

/* e.g. pix2sky 01u56798p0301o0.fits.fz 300 400 */
/*      sky2pix 01u56798p0301o0.fits.fz 187.232 -9.412 */

// For reference see:
//   fits.gsfc.nasa.gov/registry
//   web.ipac.caltech.edu/staff/shupe/reprints/SIP_to_PV_SPIE2012.pdf

/* v1.0 140611 John Tonry */
#include <stdio.h>
#include <unistd.h>
#include <stdlib.h>
#include <string.h>
#include <math.h>

#include "wcseval.h"

int VERB;

void syntax(char *prog);

int main(int argc, char **argv)
{
   char *fitsname, *dataname, line[256];
   FILE *fp=NULL;
   int i, narg, err, forcesip, forcetpv, forcecd, inverse, cmdata, tskpix;
   double x, y, ra, dec;
   WCS wcs;

/* No arguments: print syntax */
   if(argc < 2) syntax(argv[0]);

/* First unknown arg is always the name of the FITS file */
   fitsname = NULL;

/* Verbosity level */
   VERB = 0;

/* Test SIP vs PV */
   forcesip = forcetpv = forcecd = 0;

/* Map sky to xy? */
   inverse = (strstr(argv[0], "pix2sky") == NULL) ;

/* Input file with coordinates to be converted */
   dataname = NULL;

/* Command line data? */
   cmdata = 0;

/* TSK pixels?  TSK lower left center (0.5,0.5), WCS LL ctr at (1.0,1.0) */
/* Default is to conform to WCS convention */
   tskpix = 0;

/* Parse the arguments */
   for(narg=1; narg<argc; narg++) {
      if(strcmp(argv[narg], "-o") == 0) {
	 fitsname = argv[++narg];

      } else if(strcmp(argv[narg], "-pv") == 0) {
	 forcetpv = 1;

      } else if(strcmp(argv[narg], "-sip") == 0) {
	 forcesip = 1;

      } else if(strcmp(argv[narg], "-cd") == 0) {
	 forcecd = 1;
	 forcesip = forcetpv = -1;

      } else if(strcmp(argv[narg], "-sky2pix") == 0) {
	 inverse = 1;

      } else if(strcmp(argv[narg], "-pix2sky") == 0) {
	 inverse = 0;

      } else if(strcmp(argv[narg], "-TSKpix") == 0) {
	 tskpix = 1;

      } else if(strcmp(argv[narg], "-WCSpix") == 0) {
	 tskpix = 0;

      } else if(strcmp(argv[narg], "-verb") == 0) {
	 VERB = 1;

      } else if(strcmp(argv[narg], "-VERB") == 0) {
	 VERB = 2;

      } else if(strcmp(argv[narg], "-VERBOSE") == 0) {
	 VERB = 3;

      } else {
	 fitsname = argv[narg];
	 if(argc <= narg+1) syntax(argv[0]);
	 if(argc > narg+2) {
	    i  = sscanf(argv[narg+1], "%lf", &x);
	    i += sscanf(argv[narg+2], "%lf", &y);
	    if(i == 2) cmdata = 1;
	 } else {
	    dataname = argv[narg+1];
	 }
	 break;
      }   
   }

/* Read WCS information */
   if((i=getwcs(fitsname, &wcs)) != 0) {
      fprintf(stderr, "Error: return %d from getwcs\n", i);
      exit(1);
   }

/* Open the file with the coordinate data */
   if(dataname != NULL) {
      if(strcmp(dataname, "-") == 0) {
	 fp = stdin;
      } else if((fp=fopen(dataname, "r")) == NULL) {
         fprintf(stderr, "Error: %s cannot open input file '%s'\n", 
                 argv[0], dataname);
         exit(1);
      }
   } else if(!cmdata) {
      fprintf(stderr, "Error: %s has no data\n", argv[0]);
      syntax(argv[0]);
   }

/* Convert all coordinates */
   while(1) {
      if(dataname != NULL) {
	 if(fgets(line, 256, fp) == NULL) break;
	 if(sscanf(line, "%lf %lf", &x, &y) != 2) {
	    fprintf(stderr, "Error: cannot scan coords from %s\n", line);
	    exit(1);
	 }
      }
      ra = x;
      dec = y;

      if(!inverse) {

/* If TSK (LL ctr at 0.5,0.5) add 0.5 to get to WCS convention (LL ctr 1,1) */
	 if(tskpix) {
	    x += 0.5;
	    y += 0.5;
	 }

	 if(forcesip == 0 && forcetpv == 0 && forcesip == 0) {
	    err = xy2sky(&wcs, x, y, &ra, &dec);

	 } else if(forcesip == -1 && forcetpv == -1) {
	    err = xy2sky_tan(&wcs, x, y, &ra, &dec);

	 } else if((forcesip || (wcs.sipok && !wcs.pvok)) && !forcetpv) {
	    err = xy2sky_sip(&wcs, x, y, &ra, &dec);

	 } else {
	    err = xy2sky_tpv(&wcs, x, y, &ra, &dec);
	 }

	 if(err) {
	    fprintf(stderr, "Error: cannot convert coord %.4g %.4g\n", x, y);
	    exit(1);
	 }
	 printf("%10.6f %10.6f\n", ra, dec);

      } else {

	 if(forcesip == 0 && forcetpv == 0 && forcecd == 0) {
	    err = sky2xy(&wcs, ra, dec, &x, &y);

	 } else if(forcesip == -1 && forcetpv == -1) {
	    err = sky2xy_tan(&wcs, ra, dec, &x, &y);

	 } else if((forcesip || (wcs.sipok && !wcs.pvok)) && !forcetpv) {
	    err = sky2xy_sip(&wcs, ra, dec, &x, &y);

	 } else {
	    err = sky2xy_tpv(&wcs, ra, dec, &x, &y);
	 }

/* If TSK (LL ctr at 0.5,0.5) subtract 0.5 from WCS convention (LL ctr 1,1) */
	 if(tskpix) {
	    x -= 0.5;
	    y -= 0.5;
	 }

	 if(err) {
	    fprintf(stderr, "Error: cannot convert coord %.4g %.4g\n", ra, dec);
	    exit(1);
	 }
	 printf("%7.2f %7.2f\n", x, y);
      }


      if(cmdata) break;
   }

   exit(0);
}

void syntax(char *prog)
{
   printf("%s [options] fitsfile [x y | dataname]\n", prog);
   printf("\nRead WCS terms from fitsfile and convert x,y or coordinates found in\n");
   printf("  datafile to RA,Dec[deg].  fitsfile or datafile can be '-' to indicate stdin.\n");
   printf("  Default usage (-WCSpix) has pixels in WCS convention: lower left center 1,1\n");
   printf("\n[options] include\n");
   printf("  -sky2pix    Read RA,Dec[deg] and convert to x,y\n");
   printf("  -pix2sky    Read x,y and convert to RA,Dec[deg]\n");
   printf("  -pv         Force use of PV distortion coefficients if present\n");
   printf("  -sip        Force use of SIP distortion coefficients if present\n");
   printf("  -cd         Use only gnomonic projection and CD matrix\n");
   printf("  -TSKpix     x,y are in TSK (LL @ 0.5,0.5), not WCS (LL @ 1,1)\n");
   printf("  -WCSpix     x,y are in WCS (LL @ 1,1), not TSK (LL @ 0.5,0.5)\n");
   printf("  -verb       Some verbosity\n");
   printf("  -VERB       More verbosity\n");
   printf("  -VERBOSE    Full verbosity\n");
   exit(0);
}
