/* Routines to convert between pixel and celestial coordinates */
/* v1.0 140611 John Tonry */
#include <stdio.h>
#include <unistd.h>
#include <stdlib.h>
#include <sys/file.h>
#include <string.h>
#include <math.h>

#include "wcseval.h"

#define ABS(a) (((a) > 0) ? (a) : -(a))
#define MAX(a,b) (((a) > (b)) ? (a) : (b))
#define MIN(a,b) (((a) < (b)) ? (a) : (b))

#define NFITS 2880

#define DEBUG			/* Debug output? */
// #define MAINTEST		/* Stand-alone main() for testing? */

#define TOL_DOUBLE (1e-12)	/* Double precision inversion tolerance */

#ifdef MAINTEST
static int VERB=0;		/* Verbosity level */
#else
extern int VERB;		/* Verbosity level */
#endif

/* Powers of TPV1 terms 0-39.  NOTE: swap roles of x vs y for TPV2!!! */
static int xpv_pwr[MAXPV]={0, 1,0, 0, 2,1,0, 3,2,1,0, 0, 4,3,2,1,0,
			   5,4,3,2,1,0, 0, 6,5,4,3,2,1,0, 7,6,5,4,3,2,1,0, 0};
static int ypv_pwr[MAXPV]={0, 0,1, 0, 0,1,2, 0,1,2,3, 0, 0,1,2,3,4,
			   0,1,2,3,4,5, 0, 0,1,2,3,4,5,6, 0,1,2,3,4,5,6,7, 0};
static int rpv_pwr[MAXPV]={0, 0,0, 1, 0,0,0, 0,0,0,0, 3, 0,0,0,0,0,
			   0,0,0,0,0,0, 5, 0,0,0,0,0,0,0, 0,0,0,0,0,0,0,0, 7};

/* Apply a TPV polynomial */
void apply_tpv(int npv, int nx, int ny, int nr, 
	       double *xpv, double *ypv, double *xi, double *eta);
/* Apply an inverse TPV polynomial */
void apply_tpv_inverse(int npv, int nx, int ny, int nr, 
	       double *xpv, double *ypv, double *xi, double *eta);

/* Apply a ZPN radial polynomial: sky theta -> R(theta) projection plane */
void apply_zpn(int ncoeff, double *coeff, double *x, double *y);
/* Apply an inverse ZPN polynomial */
void apply_zpn_inverse(int ncoeff, double *coeff, double *xi, double *eta);

/* Apply a ZPX distortion polynomial and ZPN radial distortion */
void apply_zpx(WCS *wcs, double *xi, double *eta);
/* Apply a ZPX distortion and ZPN radial inverse: pix->sky */
void apply_zpx_inverse(WCS *wcs, double *xi, double *eta);

/* Convert x,y to RA,Dec using gnomonic projection only */
int xy2sky_tan(WCS *wcs, double x, double y, double *ra, double *dec);
/* Convert x,y to RA,Dec using stereographic projection only */
int xy2sky_stg(WCS *wcs, double x, double y, double *ra, double *dec);
/* Convert x,y to RA,Dec using orthograhic projection only */
int xy2sky_sin(WCS *wcs, double x, double y, double *ra, double *dec);
/* Convert x,y to RA,Dec using ARC projection only */
int xy2sky_arc(WCS *wcs, double x, double y, double *ra, double *dec);
/* Convert x,y to RA,Dec using ARC projection and zenithal distortion */
int xy2sky_zpn(WCS *wcs, double x, double y, double *ra, double *dec);
/* Convert x,y to RA,Dec using ARC projection and ZPX zenithal distortion */
int xy2sky_zpx(WCS *wcs, double x, double y, double *ra, double *dec);
/* Convert x,y to RA,Dec using just the CD matrix */
int xy2sky_cd(WCS *wcs, double x, double y, double *ra, double *dec);

/* Convert RA,Dec to x,y using gnomonic projection only */
int sky2xy_tan(WCS *wcs, double ra, double dec, double *x, double *y);
/* Convert RA,Dec to x,y using stereographic projection only */
int sky2xy_stg(WCS *wcs, double ra, double dec, double *x, double *y);
/* Convert RA,Dec to x,y using orthographic projection only */
int sky2xy_sin(WCS *wcs, double ra, double dec, double *x, double *y);
/* Convert RA,Dec to x,y using arc projection only */
int sky2xy_arc(WCS *wcs, double ra, double dec, double *x, double *y);
/* Convert RA,Dec to x,y using arc projection and zenithal distortion */
int sky2xy_zpn(WCS *wcs, double ra, double dec, double *x, double *y);
/* Convert RA,Dec to x,y using arc projection and ZPX zenithal distortion */
int sky2xy_zpx(WCS *wcs, double ra, double dec, double *x, double *y);
/* Convert RA,Dec to x,y using just the CD matrix */
int sky2xy_cd(WCS *wcs, double ra, double dec, double *x, double *y);

/* Convert x,y to RA,Dec [deg] */
int xy2sky(WCS *wcs, double x, double y, double *ra, double *dec)
{
   switch(wcs->projection) {
      case PROJ_TAN:
	 return(xy2sky_tan(wcs, x, y, ra, dec));

      case PROJ_TAN_TPV:
	 return(xy2sky_tpv(wcs, x, y, ra, dec));

      case PROJ_TAN_SIP:
	 return(xy2sky_sip(wcs, x, y, ra, dec));

      case PROJ_STG:
	 return(xy2sky_stg(wcs, x, y, ra, dec));

      case PROJ_SIN:
	 return(xy2sky_sin(wcs, x, y, ra, dec));

      case PROJ_ARC:
	 return(xy2sky_arc(wcs, x, y, ra, dec));

      case PROJ_ZPN:
	 return(xy2sky_zpn(wcs, x, y, ra, dec));
	 break;

      case PROJ_ZPX:
	 return(xy2sky_zpx(wcs, x, y, ra, dec));
	 break;

      case PROJ_CD:
	 return(xy2sky_cd(wcs, x, y, ra, dec));

      case PROJ_UNKNOWN:
      default:
	 fprintf(stderr, "Unknown projection %d\n", wcs->projection);
	 break;
   }
   return(-1);
}

/* Convert RA,Dec [deg] to x,y */
int sky2xy(WCS *wcs, double ra, double dec, double *x, double *y)
{
   switch(wcs->projection) {
      case PROJ_TAN:
	 return(sky2xy_tan(wcs, ra, dec, x, y));

      case PROJ_TAN_TPV:
	 return(sky2xy_tpv(wcs, ra, dec, x, y));

      case PROJ_TAN_SIP:
	 return(sky2xy_sip(wcs, ra, dec, x, y));

      case PROJ_STG:
	 return(sky2xy_stg(wcs, ra, dec, x, y));

      case PROJ_SIN:
	 return(sky2xy_sin(wcs, ra, dec, x, y));

      case PROJ_ARC:
	 return(sky2xy_arc(wcs, ra, dec, x, y));

      case PROJ_ZPN:
	 return(sky2xy_zpn(wcs, ra, dec, x, y));
	 break;

      case PROJ_ZPX:
	 return(sky2xy_zpx(wcs, ra, dec, x, y));
	 break;

      case PROJ_CD:
	 return(sky2xy_cd(wcs, ra, dec, x, y));

      case PROJ_UNKNOWN:
      default:
	 fprintf(stderr, "Unknown projection %d\n", wcs->projection);
	 break;
   }
   return(-1);
}

/* Convert x,y to RA,Dec using SIP distortion terms */
int xy2sky_sip(WCS *wcs, double x, double y, double *ra, double *dec)
{
   int i, j, n;
   double u, v, xi, eta;
   double upow[10], vpow[10];

/* Initialize power terms */
   upow[0] = vpow[0] = 1.0;

   u = x - wcs->crpix1;
   v = y - wcs->crpix2;
   wcs->u = u;
   wcs->v = v;

/* If SIP, apply the polynomial */
   for(j=1; j<=wcs->usipmax; j++) upow[j] = upow[j-1] * u;
   for(j=1; j<=wcs->vsipmax; j++) vpow[j] = vpow[j-1] * v;

   for(i=0; i<=wcs->usipmax; i++) {
      for(j=0; i+j<=MAX(wcs->usipmax,wcs->vsipmax); j++) {
	 n = ((i+j)*((i+j)+1))/2 + j;
	 u += wcs->sipa[n] * upow[i] * vpow[j];
	 v += wcs->sipb[n] * upow[i] * vpow[j];
#ifdef DEBUG
	 if(VERB > 1) {
	    printf("%2d  %8.3f %8.3f  A_%d_%d %10.3e  B_%d_%d %10.3e\n", 
		   n, u,v, i,j, wcs->sipa[n], i,j, wcs->sipb[n]);
	 }
#endif
      }
   }

/* Apply the CD matrix */
   xi  = wcs->cd11*u + wcs->cd12*v;
   eta = wcs->cd21*u + wcs->cd22*v;

/* Do the gnomonic projection from tangent plane to sky */
   tp2sky_tan(xi, eta, wcs->crval1, wcs->crval2, ra, dec);

/* In case somebody is interested... */
   wcs->xi = xi;
   wcs->eta = eta;
   return(0);
}

/* Convert x,y to RA,Dec using PV distortion terms */
int xy2sky_tpv(WCS *wcs, double x, double y, double *ra, double *dec)
{
   double u, v, xi, eta;

/* Start the conversion with the CRPIX offset */
   u = x - wcs->crpix1;
   v = y - wcs->crpix2;

/* Apply the CD matrix */
   xi  = wcs->cd11*u + wcs->cd12*v;
   eta = wcs->cd21*u + wcs->cd22*v;

/* If PV, apply the polynomial */
   apply_tpv(wcs->npv, wcs->xpvmax, wcs->ypvmax, wcs->rpvmax, 
	     wcs->pv[0], wcs->pv[1], &xi, &eta);

/* Do the gnomonic projection from tangent plane to sky */
   tp2sky_tan(xi, eta, wcs->crval1, wcs->crval2, ra, dec);

/* In case somebody is interested... */
   wcs->u = u;
   wcs->v = v;
   wcs->xi = xi;
   wcs->eta = eta;

   return(0);
}

/* Convert x,y to RA,Dec using gnomonic projection only */
int xy2sky_tan(WCS *wcs, double x, double y, double *ra, double *dec)
{
   double u, v, xi, eta;

   u = x - wcs->crpix1;
   v = y - wcs->crpix2;
   wcs->u = u;
   wcs->v = v;

/* Apply the CD matrix */
   xi  = wcs->cd11*u + wcs->cd12*v;
   eta = wcs->cd21*u + wcs->cd22*v;

/* Do the gnomonic projection from tangent plane to sky */
   tp2sky_tan(xi, eta, wcs->crval1, wcs->crval2, ra, dec);

/* In case somebody is interested... */
   wcs->xi = xi;
   wcs->eta = eta;
   return(0);
}

/* Convert x,y to RA,Dec using just the CD matrix */
int xy2sky_cd(WCS *wcs, double x, double y, double *ra, double *dec)
{
   double u, v;

   u = x - wcs->crpix1;
   v = y - wcs->crpix2;
   wcs->u = u;
   wcs->v = v;

/* Apply the CD matrix */
   *ra  = wcs->cd11*u + wcs->cd12*v + wcs->crval1;
   *dec = wcs->cd21*u + wcs->cd22*v + wcs->crval2;

/* In case somebody is interested... */
   wcs->xi = *ra - wcs->crval1;
   wcs->eta = *dec - wcs->crval2;
   return(0);
}

/* Convert x,y to RA,Dec using stereographic projection only */
int xy2sky_stg(WCS *wcs, double x, double y, double *ra, double *dec)
{
   double u, v, xi, eta;

   u = x - wcs->crpix1;
   v = y - wcs->crpix2;
   wcs->u = u;
   wcs->v = v;

/* Apply the CD matrix */
   xi  = wcs->cd11*u + wcs->cd12*v;
   eta = wcs->cd21*u + wcs->cd22*v;

/* Do the stereographic projection from tangent plane to sky */
   tp2sky_stg(xi, eta, wcs->crval1, wcs->crval2, ra, dec);

/* In case somebody is interested... */
   wcs->xi = xi;
   wcs->eta = eta;
   return(0);
}

/* Convert x,y to RA,Dec using orthograhic projection only */
int xy2sky_sin(WCS *wcs, double x, double y, double *ra, double *dec)
{
   double u, v, xi, eta;

   u = x - wcs->crpix1;
   v = y - wcs->crpix2;
   wcs->u = u;
   wcs->v = v;

/* Apply the CD matrix */
   xi  = wcs->cd11*u + wcs->cd12*v;
   eta = wcs->cd21*u + wcs->cd22*v;

/* Do the orthographic projection from tangent plane to sky */
   tp2sky_sin(xi, eta, wcs->crval1, wcs->crval2, ra, dec);

/* In case somebody is interested... */
   wcs->xi = xi;
   wcs->eta = eta;
   return(0);
}

/* Convert x,y to RA,Dec using ARC projection only */
int xy2sky_arc(WCS *wcs, double x, double y, double *ra, double *dec)
{
   double u, v, xi, eta;

   u = x - wcs->crpix1;
   v = y - wcs->crpix2;
   wcs->u = u;
   wcs->v = v;

/* Apply the CD matrix */
   xi  = wcs->cd11*u + wcs->cd12*v;
   eta = wcs->cd21*u + wcs->cd22*v;

/* Do the arc projection from tangent plane to sky */
   tp2sky_arc(xi, eta, wcs->crval1, wcs->crval2, ra, dec);

/* In case somebody is interested... */
   wcs->xi = xi;
   wcs->eta = eta;
   return(0);
}

/* Convert x,y to RA,Dec using ARC projection and zenithal distortion */
int xy2sky_zpn(WCS *wcs, double x, double y, double *ra, double *dec)
{
   double u, v, xi, eta;

   u = x - wcs->crpix1;
   v = y - wcs->crpix2;
   wcs->u = u;
   wcs->v = v;

/* Apply the CD matrix */
   xi  = wcs->cd11*u + wcs->cd12*v;
   eta = wcs->cd21*u + wcs->cd22*v;

/* If PV, apply the inverse polynomial: pixels R(theta) to sky angle theta */
   apply_zpn_inverse(wcs->zprmax, wcs->zpr, &xi, &eta);

/* Do the arc projection from tangent plane to sky */
   tp2sky_arc(xi, eta, wcs->crval1, wcs->crval2, ra, dec);

/* In case somebody is interested... */
   wcs->xi = xi;
   wcs->eta = eta;
   return(0);
}

/* Convert x,y to RA,Dec using ARC projection and zenithal distortion */
int xy2sky_zpx(WCS *wcs, double x, double y, double *ra, double *dec)
{
   double u, v, xi, eta;

   u = x - wcs->crpix1;
   v = y - wcs->crpix2;
   wcs->u = u;
   wcs->v = v;

/* Apply the CD matrix */
   xi  = wcs->cd11*u + wcs->cd12*v;
   eta = wcs->cd21*u + wcs->cd22*v;

/* If PV, apply the polynomial */
   apply_zpx_inverse(wcs, &xi, &eta);

/* Do the arc projection from tangent plane to sky */
   tp2sky_arc(xi, eta, wcs->crval1, wcs->crval2, ra, dec);

/* In case somebody is interested... */
   wcs->xi = xi;
   wcs->eta = eta;
   return(0);
}

/* Convert RA,Dec to x,y using SIP distortion terms */
int sky2xy_sip(WCS *wcs, double ra, double dec, double *x, double *y)
{
   int i, j, n;
   double u, v, xi, eta;
   double upow[10], vpow[10];

/* Initialize power terms */
   upow[0] = vpow[0] = 1.0;

/* Do the gnomonic projection from sky to tangent plane */
   sky2tp_tan(ra, dec, wcs->crval1, wcs->crval2, &xi, &eta);

/* Apply the inverse CD matrix */
   u = wcs->dc11*xi + wcs->dc12*eta;
   v = wcs->dc21*xi + wcs->dc22*eta;

/* If SIP, apply the inverse polynomial */
   for(j=1; j<=wcs->usippmax; j++) upow[j] = upow[j-1] * u;
   for(j=1; j<=wcs->vsippmax; j++) vpow[j] = vpow[j-1] * v;

   for(i=0; i<=wcs->usippmax; i++) {
      for(j=0; i+j<=MAX(wcs->usippmax,wcs->vsippmax); j++) {
	 n = ((i+j)*((i+j)+1))/2 + j;
	 u += wcs->sipap[n] * upow[i] * vpow[j];
	 v += wcs->sipbp[n] * upow[i] * vpow[j];
#ifdef DEBUG
	 if(VERB > 1) {
	    printf("%2d  %8.3f %8.3f  AP_%d_%d %10.3e  BP_%d_%d %10.3e\n", 
		   n, u,v, i,j, wcs->sipap[n], i,j, wcs->sipbp[n]);
	 }
#endif
      }
   }

/* Finish the conversion with the CRPIX offset */
   *x = u + wcs->crpix1;
   *y = v + wcs->crpix2;

/* In case somebody is interested... */
   wcs->u = u;
   wcs->v = v;
   wcs->xi = xi;
   wcs->eta = eta;
   return(0);
}

/* Convert RA,Dec to x,y using PV distortion terms */
int sky2xy_tpv(WCS *wcs, double ra, double dec, double *x, double *y)
{
   double u, v, xi, eta;

/* Do the gnomonic projection from sky to tangent plane */
   sky2tp_tan(ra, dec, wcs->crval1, wcs->crval2, &xi, &eta);
   wcs->xi = xi;
   wcs->eta = eta;

/* Invert the polynomial */
//   apply_tpv_inverse_revert(wcs->npv, wcs->xpvmax, wcs->ypvmax, wcs->rpvmax, 
//	     wcs->pv[0], wcs->pv[1], &xi, &eta);
   apply_tpv_inverse(wcs->npv, wcs->xpvmax, wcs->ypvmax, wcs->rpvmax, 
	     wcs->pv[0], wcs->pv[1], &xi, &eta);

/* Apply the inverse CD matrix */
   u = wcs->dc11*xi + wcs->dc12*eta;
   v = wcs->dc21*xi + wcs->dc22*eta;

/* Finish the conversion with the CRPIX offset */
   *x = u + wcs->crpix1;
   *y = v + wcs->crpix2;

/* In case somebody is interested... */
   wcs->u = u;
   wcs->v = v;

   return(0);
}


/* Convert RA,Dec to x,y using gnomonic projection only */
int sky2xy_tan(WCS *wcs, double ra, double dec, double *x, double *y)
{
   double u, v, xi, eta;

/* Do the gnomonic projection from sky to tangent plane */
   sky2tp_tan(ra, dec, wcs->crval1, wcs->crval2, &xi, &eta);

/* Apply the inverse CD matrix */
   u = wcs->dc11*xi + wcs->dc12*eta;
   v = wcs->dc21*xi + wcs->dc22*eta;

/* Finish the conversion with the CRPIX offset */
   *x = u + wcs->crpix1;
   *y = v + wcs->crpix2;

/* In case somebody is interested... */
   wcs->u = u;
   wcs->v = v;
   wcs->xi = xi;
   wcs->eta = eta;
   return(0);
}

/* Convert RA,Dec to x,y using just the CD matrix */
int sky2xy_cd(WCS *wcs, double ra, double dec, double *x, double *y)
{
   double u, v;

/* Apply the inverse CD matrix */
   u = wcs->dc11*(ra-wcs->crval1) + wcs->dc12*(dec-wcs->crval2);
   v = wcs->dc21*(ra-wcs->crval1) + wcs->dc22*(dec-wcs->crval2);

/* Finish the conversion with the CRPIX offset */
   *x = u + wcs->crpix1;
   *y = v + wcs->crpix2;

/* In case somebody is interested... */
   wcs->u = u;
   wcs->v = v;
   wcs->xi = ra - wcs->crval1;
   wcs->eta = dec - wcs->crval2;
   return(0);
}

/* Convert RA,Dec to x,y using stereographic projection only */
int sky2xy_stg(WCS *wcs, double ra, double dec, double *x, double *y)
{
   double u, v, xi, eta;

/* Do the stereographic projection from sky to tangent plane */
   sky2tp_stg(ra, dec, wcs->crval1, wcs->crval2, &xi, &eta);

/* Apply the inverse CD matrix */
   u = wcs->dc11*xi + wcs->dc12*eta;
   v = wcs->dc21*xi + wcs->dc22*eta;

/* Finish the conversion with the CRPIX offset */
   *x = u + wcs->crpix1;
   *y = v + wcs->crpix2;

/* In case somebody is interested... */
   wcs->u = u;
   wcs->v = v;
   wcs->xi = xi;
   wcs->eta = eta;
   return(0);
}

/* Convert RA,Dec to x,y using orthographic projection only */
int sky2xy_sin(WCS *wcs, double ra, double dec, double *x, double *y)
{
   double u, v, xi, eta;

/* Do the orthographic projection from sky to tangent plane */
   sky2tp_sin(ra, dec, wcs->crval1, wcs->crval2, &xi, &eta);

/* Apply the inverse CD matrix */
   u = wcs->dc11*xi + wcs->dc12*eta;
   v = wcs->dc21*xi + wcs->dc22*eta;

/* Finish the conversion with the CRPIX offset */
   *x = u + wcs->crpix1;
   *y = v + wcs->crpix2;

/* In case somebody is interested... */
   wcs->u = u;
   wcs->v = v;
   wcs->xi = xi;
   wcs->eta = eta;
   return(0);
}

/* Convert RA,Dec to x,y using arc projection only */
int sky2xy_arc(WCS *wcs, double ra, double dec, double *x, double *y)
{
   double u, v, xi, eta;

/* Do the arc projection from sky to tangent plane */
   sky2tp_arc(ra, dec, wcs->crval1, wcs->crval2, &xi, &eta);

/* Apply the inverse CD matrix */
   u = wcs->dc11*xi + wcs->dc12*eta;
   v = wcs->dc21*xi + wcs->dc22*eta;

/* Finish the conversion with the CRPIX offset */
   *x = u + wcs->crpix1;
   *y = v + wcs->crpix2;

/* In case somebody is interested... */
   wcs->u = u;
   wcs->v = v;
   wcs->xi = xi;
   wcs->eta = eta;
   return(0);
}

/* Convert RA,Dec to x,y using arc projection and zenithal distortion */
int sky2xy_zpn(WCS *wcs, double ra, double dec, double *x, double *y)
{
   double u, v, xi, eta;

/* Do the arc projection from sky to tangent plane */
   sky2tp_arc(ra, dec, wcs->crval1, wcs->crval2, &xi, &eta);
   wcs->xi = xi;
   wcs->eta = eta;

/* Apply the polynomial: sky theta to pixels R(theta) */
   apply_zpn(wcs->zprmax, wcs->zpr, &xi, &eta);

/* Apply the inverse CD matrix */
   u = wcs->dc11*xi + wcs->dc12*eta;
   v = wcs->dc21*xi + wcs->dc22*eta;

/* Finish the conversion with the CRPIX offset */
   *x = u + wcs->crpix1;
   *y = v + wcs->crpix2;

/* In case somebody is interested... */
   wcs->u = u;
   wcs->v = v;
   wcs->xi = xi;
   wcs->eta = eta;
   return(0);
}

/* Convert RA,Dec to x,y using arc projection and ZPX zenithal distortion */
int sky2xy_zpx(WCS *wcs, double ra, double dec, double *x, double *y)
{
   double u, v, xi, eta;

/* Do the arc projection from sky to tangent plane */
   sky2tp_arc(ra, dec, wcs->crval1, wcs->crval2, &xi, &eta);
   wcs->xi = xi;
   wcs->eta = eta;

/* Apply the polynomial: sky theta to pixels R(theta) */
   apply_zpx(wcs, &xi, &eta);

/* Apply the inverse CD matrix */
   u = wcs->dc11*xi + wcs->dc12*eta;
   v = wcs->dc21*xi + wcs->dc22*eta;

/* Finish the conversion with the CRPIX offset */
   *x = u + wcs->crpix1;
   *y = v + wcs->crpix2;

//   printf("%12.7f %12.7f\n", *x, *y);

/* In case somebody is interested... */
   wcs->u = u;
   wcs->v = v;
   wcs->xi = xi;
   wcs->eta = eta;
   return(0);
}

/* Apply a TPV polynomial */
void apply_tpv(int npv, int nx, int ny, int nr, 
	       double *xpv, double *ypv, double *xi, double *eta)
{
   int j;
   double upow[10], vpow[10], rpow[10], r;


   upow[0] = vpow[0] = rpow[0] = 1.0;		/* Initialize power terms */
   for(j=1; j<=nx; j++) upow[j] = upow[j-1] * (*xi);
   for(j=1; j<=ny; j++) vpow[j] = vpow[j-1] * (*eta);
   if(nr > 0) {
      r = sqrt((*xi)*(*xi) + (*eta)*(*eta));
      for(j=1; j<=nr; j++) rpow[j] = rpow[j-1] * r;
   }
/* Apply the polynomial */
   *xi = *eta = 0.0;
   for(j=0; j<=npv; j++) {
      if(xpv[j] == 0 && ypv[j] == 0) continue;
      if(rpv_pwr[j] == 0) {
	 *xi  += xpv[j] * upow[xpv_pwr[j]] * vpow[ypv_pwr[j]];
	 *eta += ypv[j] * upow[ypv_pwr[j]] * vpow[xpv_pwr[j]];
      } else if(nr >= rpv_pwr[j]) {
	 *xi  += xpv[j] * rpow[rpv_pwr[j]];
	 *eta += ypv[j] * rpow[rpv_pwr[j]];
      }
#ifdef DEBUG
      if(VERB > 1) {
	 printf("%9.5f %9.5f  PV1_%d %10.3e   PV2_%d %10.3e  %d %d %d\n", 
		*xi, *eta, j, xpv[j], j, ypv[j], 
		xpv_pwr[j], ypv_pwr[j], rpv_pwr[j]);
      }
#endif
   }
}

/* Apply an inverse TPV polynomial */
void apply_tpv_inverse(int npv, int nx, int ny, int nr, 
	       double *xpv, double *ypv, double *xi, double *eta)
{
   int iter, j;
   double upow[10], vpow[10], rpow[10], r, u, v, det;
   double xp, yp, dxdu, dxdv, dydu, dydv;


   upow[0] = vpow[0] = rpow[0] = r = 1.0;	/* Initialize power terms */
   u = *xi;					/* Initialize solution */
   v = *eta;
   for(iter=0; iter<10; iter++) {
      for(j=1; j<=nx; j++) upow[j] = upow[j-1] * u;
      for(j=1; j<=ny; j++) vpow[j] = vpow[j-1] * v;
      if(nr > 0) {
	 r = sqrt(u*u + v*v);
	 for(j=1; j<=nr; j++) rpow[j] = rpow[j-1] * r;
      }
      xp = yp = dxdu = dxdv = dydu = dydv = 0.0;
/* Evaluate function and Jacobian at u,v */
      for(j=0; j<=npv; j++) {
	 if(rpv_pwr[j] == 0) {
	    xp += xpv[j] * upow[xpv_pwr[j]] * vpow[ypv_pwr[j]];
	    yp += ypv[j] * upow[ypv_pwr[j]] * vpow[xpv_pwr[j]];
	    if(xpv_pwr[j] > 0) dxdu += xpv[j] * xpv_pwr[j] * 
				  upow[xpv_pwr[j]-1] * vpow[ypv_pwr[j]];
	    if(ypv_pwr[j] > 0) dxdv += xpv[j] * ypv_pwr[j] * 
				  upow[xpv_pwr[j]] * vpow[ypv_pwr[j]-1];
	    if(ypv_pwr[j] > 0) dydu += ypv[j] * ypv_pwr[j] * 
				  upow[ypv_pwr[j]-1] * vpow[xpv_pwr[j]];
	    if(xpv_pwr[j] > 0) dydv += ypv[j] * xpv_pwr[j] * 
				  upow[ypv_pwr[j]] * vpow[xpv_pwr[j]-1];
	 } else if(nr > rpv_pwr[j]) {
	    xp += xpv[j] * rpow[rpv_pwr[j]];
	    yp += ypv[j] * rpow[rpv_pwr[j]];
	    dxdu += xpv[j] * rpv_pwr[j] * rpow[rpv_pwr[j]-1] * u/(r*r);
	    dxdv += xpv[j] * rpv_pwr[j] * rpow[rpv_pwr[j]-1] * v/(r*r);
	    dydu += ypv[j] * rpv_pwr[j] * rpow[rpv_pwr[j]-1] * u/(r*r);
	    dydv += ypv[j] * rpv_pwr[j] * rpow[rpv_pwr[j]-1] * v/(r*r);
	 }
      }
      det = dxdu*dydv - dxdv*dydu;
      u -= (dydv*(xp-(*xi)) - dydu*(yp-(*eta))) / det;
      v -= (dxdu*(yp-(*eta)) - dxdv*(xp-(*xi))) / det;
#ifdef DEBUG
      if(VERB > 1) {
	 printf("%9.7f %9.7f  %9.7f %9.7f  %9.7f %9.7f  %10.3e %10.3e\n", 
		*xi, *eta, xp, yp, u, v, *xi-xp, *eta-yp);
      }
#endif
      if(ABS(xp-(*xi))<TOL_DOUBLE && ABS(yp-(*eta))<TOL_DOUBLE) break;
   }
   *xi  = u;
   *eta = v;
}

/* ZPN definition:
 * This is a little squirrelly.  According to C&G, 
 *
 *   R(theta) = 180/pi * Sum_0^20 Pm * ((pi/180)*(90-theta))^m
 *
 * but although P0 is semi-meaningless (it pushes the ZPN off the pole
 * so you can't get to certain areas) C&G reserve a bizarre meaning
 * for it.  They say that the Pm are given by the latitude coordinate
 * (2, I think) as PV2_0, PV2_1, ... PV2_20.  Note that for C&G theta
 * is measured from the equator.
 */

/* Apply a ZPN radial polynomial: sky theta -> R(theta) projection plane */
void apply_zpn(int ncoeff, double *coeff, double *x, double *y)
{
   int j;
   double r, R, dr=atan(1.0)/45;

/* coeff[1] to coeff[ncoeff] are r^1 to r^n coefficients */
/* Deliberately avoid the constant term */
   r = dr * sqrt((*x)*(*x) + (*y)*(*y));
   for(j=ncoeff, R=0.0; j>0; j--) R = r * (R + coeff[j]);

   *x *= R/r;
   *y *= R/r;
}

/* Apply an inverse ZPN polynomial: proj plane R(theta) -> theta sky */
void apply_zpn_inverse(int ncoeff, double *coeff, double *xi, double *eta)
{
   int i, j;
   double dr=atan(1.0)/45, r, rp, R, Rp, deriv;

   R = dr * sqrt((*xi)*(*xi) + (*eta)*(*eta));

   r = rp = R / coeff[1];
   for(i=0; i<10; i++) {
      for(j=ncoeff, Rp=0.0; j>=1; j--) Rp = r * (Rp + coeff[j]);
      for(j=ncoeff, deriv=0.0; j>=1; j--) deriv = r * deriv + j*coeff[j];
      r += (R - Rp) / deriv;
      if(ABS(r-rp) < TOL_DOUBLE) break;
      rp = r;
      if(VERB > 1) {
	 printf("%d %8.5f %8.5f %8.5f %10.3e\n", i, rp, R, r, r-rp);
      }
   }
   *xi  *= r/R;
   *eta *= r/R;
}

/* Apply a ZPX distortion polynomial and ZPN radial distortion: sky->pix */
void apply_zpx(WCS *wcs, double *xi, double *eta)
{
//   printf("%10.8f %10.8f", *xi, *eta);

/* Apply the ZPN distortion using the radial terms */
   apply_zpn(wcs->zprmax, wcs->zpr, xi, eta);

//   printf(" %10.8f %10.8f", *xi, *eta);

/* ZPX is an augmentation of xi,eta, so adjust the coefficients temporarily */
   wcs->pv[0][1] += 1;
   wcs->pv[1][1] += 1;

/* Apply the inverse TPV polynomial, suppress evaluation of radial terms */
   apply_tpv_inverse(wcs->npv, wcs->xpvmax, wcs->ypvmax, 0,
	     wcs->pv[0], wcs->pv[1], xi, eta);

//   printf(" %10.8f %10.8f\n", *xi, *eta);

/* Restore the linear coefficients */
   wcs->pv[0][1] -= 1;
   wcs->pv[1][1] -= 1;
}

/* Apply a ZPX distortion and ZPN radial inverse: pix->sky */
void apply_zpx_inverse(WCS *wcs, double *xi, double *eta)
{
//   printf("start: %10.8f %10.8f\n", *xi, *eta);

/* ZPX is an augmentation of xi,eta, so adjust the coefficients temporarily */
   wcs->pv[0][1] += 1;
   wcs->pv[1][1] += 1;

/* Apply the TPV polynomial, suppress evaluation of radial terms */
   apply_tpv(wcs->npv, wcs->xpvmax, wcs->ypvmax, 0,
	     wcs->pv[0], wcs->pv[1], xi, eta);

//   printf(" +tpv %10.8f %10.8f\n", *xi, *eta);

/* Restore the linear coefficients */
   wcs->pv[0][1] -= 1;
   wcs->pv[1][1] -= 1;

/* Apply the inverse ZPN distortion using the radial terms */
   apply_zpn_inverse(wcs->zprmax, wcs->zpr, xi, eta);
//   printf("  %10.8f %10.8f\n", *xi, *eta);
}


/* Initialize a WCS to a vanilla projection */
void wcs_init(
   double a0, double d0,        // [deg] center point on sphere
   double cx, double cy,        // [pix] center point on detector
   double scale,                // [deg/pix] scale
   int proj,                    // projection type from wcseval.h
   WCS *wcs)                    // wcs structure
{
   int i, j;
   
/* [deg] Convergence tolerance for inverse */
   wcs->TOL = 1e-6;

/* Initialize all the mapping variables */
   wcs->crval1 = a0;
   wcs->crval2 = d0;
   wcs->crpix1 = cx;
   wcs->crpix2 = cy;
   wcs->cd11 = wcs->cd22 = scale;
   wcs->cd12 = wcs->cd21 = 0.0;

/* Inverse of CD matrix */
   wcs->dc11 =  wcs->cd22 / (wcs->cd11*wcs->cd22-wcs->cd12*wcs->cd21);
   wcs->dc22 =  wcs->cd11 / (wcs->cd11*wcs->cd22-wcs->cd12*wcs->cd21);
   wcs->dc12 = -wcs->cd12 / (wcs->cd11*wcs->cd22-wcs->cd12*wcs->cd21);
   wcs->dc21 = -wcs->cd21 / (wcs->cd11*wcs->cd22-wcs->cd12*wcs->cd21);

/* Defaults for PV terms */
   for(i=0; i<2; i++) {
      for(j=0; j<MAXPV; j++) wcs->pv[i][j] = 0.0;
   }
   wcs->pv[0][0] = wcs->pv[0][2] = wcs->pv[1][0] = wcs->pv[1][2] = 0.0;
   wcs->pv[0][1] = wcs->pv[1][1] = 1.0;

/* Defaults for ZPN terms */
   for(j=0; j<MAXZPN; j++) wcs->zpr[j] = 0.0;
   wcs->zprmax = 0;

/* Defaults for SIP terms */
   for(j=0; j<MAXSIP; j++) 
      wcs->sipa[j] = wcs->sipb[j] = wcs->sipap[j] = wcs->sipbp[j] = 0.0;

/* PV and SIP polynomial orders */
   wcs->usipmax = wcs->vsipmax = wcs->usippmax = wcs->vsippmax = 0;
   wcs->npv = wcs->xpvmax = wcs->ypvmax = 1;
   wcs->rpvmax = 0;

/* Projection type */
   wcs->projection = proj;

   return;
}


/* Read the WCS terms from a FITS file's header */
int getwcs(char *fitsname, WCS *wcs)
{
   int i, j, k, n, fdin, nhead, endseen, nchunk, horrible_crota=0;
   int x1=0, x2=0, y1=0, y2=0;
   int naxis, naxes[10], hdrcnt, fpack, znaxis, znaxes[10];
   int gnomonic;
   int pvdistort, sipdistort;
   char header[NFITS];
   double crota, cdelt1, cdelt2, dr=atan(1.0)/45.0;
   int aorder, border, aporder, bporder;
   char watstring[2][NFITS], fmt[80], *q1, *q2;

/* [deg] Convergence tolerance for inverse */
   wcs->TOL = 1e-6;

/* Initialize all the mapping variables */
   wcs->crval1 = wcs->crval2 = wcs->crpix1 = wcs->crpix2 = 0.0;
   wcs->cd11 = wcs->cd22 = 1.0;
   wcs->cd12 = wcs->cd21 = 0.0;
   crota = 0.0;
   cdelt1 = cdelt2 = 1.0;

/* Defaults for PV terms */
   for(i=0; i<2; i++) {
      for(j=0; j<MAXPV; j++) wcs->pv[i][j] = 0.0;
   }
   wcs->pv[0][0] = wcs->pv[0][2] = wcs->pv[1][0] = wcs->pv[1][2] = 0.0;
   wcs->pv[0][1] = wcs->pv[1][1] = 1.0;
   pvdistort = 0;

/* Defaults for ZPN terms */
   for(j=0; j<MAXZPN; j++) wcs->zpr[j] = 0.0;
   wcs->zprmax = 0;

/* Defaults for SIP terms */
   for(j=0; j<MAXSIP; j++) 
      wcs->sipa[j] = wcs->sipb[j] = wcs->sipap[j] = wcs->sipbp[j] = 0.0;
   sipdistort = 0;

/* PV and SIP polynomial orders */
   aorder = border = aporder = bporder = 0;
   wcs->usipmax = wcs->vsipmax = wcs->usippmax = wcs->vsippmax = 0;
   wcs->npv = wcs->xpvmax = wcs->ypvmax = 1;
   wcs->rpvmax = 0;

/* Defaults for WAT info */
   watstring[0][0] = '\0';
   watstring[1][0] = '\0';

/* Open the file with the WCS information */
   if(strcmp(fitsname, "-") == 0) {
      fdin = STDIN_FILENO;

/* Failed to open? */
   } else if((fdin=open(fitsname, O_RDONLY, 0)) < 0) {

/* Test for cfitsio file.fits[x1:x2,y1:y2] subarray extended file syntax */
      if(sscanf(fitsname, "%*[^[][%d:%d,%d:%d]", &x1, &x2, &y1, &y2) != 4) {
         fprintf(stderr,"wcseval: cannot read from file '%s'\n", fitsname);
         return(-1);
      }
      if(VERB > 0) {
         printf("scanf: extended fits subarray %d %d %d %d\n", x1, x2, y1, y2);
      }
      q2 = malloc(strlen(fitsname)+1);
      strcpy(q2, fitsname);
      q1 = index(q2, '[');
      *q1 = '\0';
      if((fdin=open(q2, O_RDONLY, 0)) < 0) {
         fprintf(stderr,"wcseval: cannot read from file '%s'\n", q2);
         return(-1);
      }
      free(q2);
   }

   wcs->fname = (char *)malloc(strlen(fitsname)+1);
   strcpy(wcs->fname, fitsname);

/* Read the header in FITS-sized chunks */
   fpack = znaxis = 0;
   nhead = endseen = hdrcnt = 0;
   for(nchunk=0;  ; nchunk++) {
      if( read(fdin, header, NFITS) != NFITS) {
/* Short read and we've seen an END?  Call it good and continue */
	 if(endseen) {
	    nhead = endseen;
	    break;
	 }
	 fprintf(stderr, "wcseval: short read from %s\n", fitsname);
	 exit(1);
      }
/* Conforming FITS file? */
      if(nchunk == 0 && strncmp(header, "SIMPLE  =", 8) != 0) {
	 fprintf(stderr, "wcseval: does not appear to be FITS file\n");
	 exit(1);
      }
/* Scan through the header for various things... */
      for(k=0; k<NFITS; k+=80) {

/* Look for NAXIS; naxis=0 means compressed file */
	 if(strncmp(header+k, "NAXIS", 5) == 0) {
	    if(strncmp(header+k, "NAXIS   ", 8) == 0) {
	       sscanf(header+k, "NAXIS   = %d", &naxis);
	       if(naxis == 0) fpack = 1;
	    }
	    if(sscanf(header+k, "NAXIS%d", &j) == 1) {
	       sscanf(header+k+9, "%d", &naxes[j]);
	    }
	 }

/* Look for ZNAXIS and its little pals in case of compression */
	 if(strncmp(header+k, "ZNAXIS", 6) == 0) {
	    if(strncmp(header+k, "ZNAXIS   ", 8) == 0)
	       sscanf(header+k, "ZNAXIS  = %d", &znaxis);
	    if(sscanf(header+k, "ZNAXIS%d", &j) == 1) {
	       sscanf(header+k+9, "%d", &znaxes[j]);
	    }
	 }

/* Pick up all the various values for WCS mapping that may be present */

/* Coordinate reference points */
	 if(strncmp(header+k, "CRVAL1  ", 8) == 0) {
	    sscanf(header+k+9, "%lf", &wcs->crval1);
	 } else if(strncmp(header+k, "CRVAL2  ", 8) == 0) {
	    sscanf(header+k+9, "%lf", &wcs->crval2);
	 } else if(strncmp(header+k, "CRPIX1  ", 8) == 0) {
	    sscanf(header+k+9, "%lf", &wcs->crpix1);
            if(x1 > 0) wcs->crpix1 -= x1-1;       // subarray?
	 } else if(strncmp(header+k, "CRPIX2  ", 8) == 0) {
	    sscanf(header+k+9, "%lf", &wcs->crpix2);
            if(y1 > 0) wcs->crpix2 -= y1-1;       // subarray?

/* Old style CDELT and CROTA2 */
	 } else if(strncmp(header+k, "CROTA2  ", 8) == 0) {
	    sscanf(header+k+9, "%lf", &crota);
	    horrible_crota += 1<<0;
	 } else if(strncmp(header+k, "CDELT1  ", 8) == 0) {
	    sscanf(header+k+9, "%lf", &cdelt1);
	    horrible_crota += 1<<1;
	 } else if(strncmp(header+k, "CDELT2  ", 8) == 0) {
	    sscanf(header+k+9, "%lf", &cdelt2);
	    horrible_crota += 1<<2;
/* No longer used, I trust? */
//      } else if(strncmp(header+k, "CROTA   ", 8) == 0) {
//	 sscanf(header+k+9, "%lf", &crota);
//      } else if(strncmp(header+k, "CRDELT1 ", 8) == 0) {
//	 sscanf(header+k+9, "%lf", &cdelt1);
//      } else if(strncmp(header+k, "CRDELT2 ", 8) == 0) {
//	 sscanf(header+k+9, "%lf", &cdelt2);

/* CD matrix */
	 } else if(strncmp(header+k, "CD1_1   ", 8) == 0) {
	    sscanf(header+k+9, "%lf", &wcs->cd11);
	    horrible_crota += 1<<4;
	 } else if(strncmp(header+k, "CD1_2   ", 8) == 0) {
	    sscanf(header+k+9, "%lf", &wcs->cd12);
	    horrible_crota += 1<<5;
	 } else if(strncmp(header+k, "CD2_1   ", 8) == 0) {
	    sscanf(header+k+9, "%lf", &wcs->cd21);
	    horrible_crota += 1<<6;
	 } else if(strncmp(header+k, "CD2_2   ", 8) == 0) {
	    sscanf(header+k+9, "%lf", &wcs->cd22);
	    horrible_crota += 1<<7;

/* Projection CTYPE */
	 } else if(strncmp(header+k, "CTYPE1  ", 8) == 0) {
	    q1 = index(header+k, '\'');
	    q2 = index(q1+1, '\'');
	    if(q2-q1 < 8 || q2-q1 >13 ) {
	       fprintf(stderr, "Problematic CTYPE `%60.60s` has wrong length\n", 
		       header+k);
	       if(strcmp(fitsname, "-") != 0) close(fdin);
	       return(-1);
	    }
	    strncpy(wcs->ctype1, q1+1, q2-q1-1);
	    wcs->ctype1[q2-q1-1] = '\0';

	 } else if(strncmp(header+k, "CTYPE2  ", 8) == 0) {
	    q1 = index(header+k, '\'');
	    q2 = index(q1+1, '\'');
	    if(q2-q1 < 8 || q2-q1 >13 ) {
	       fprintf(stderr, "Problematic CTYPE `%60.60s` has wrong length\n", 
		       header+k);
	       if(strcmp(fitsname, "-") != 0) close(fdin);
	       return(-1);
	    }
	    strncpy(wcs->ctype2, q1+1, q2-q1-1);
	    wcs->ctype2[q2-q1-1] = '\0';

/* SIP distortion order */
	 } else if(strncmp(header+k, "A_ORDER ", 8) == 0) {
	    sscanf(header+k+9, "%d", &aorder);
	 } else if(strncmp(header+k, "B_ORDER ", 8) == 0) {
	    sscanf(header+k+9, "%d", &border);
	 } else if(strncmp(header+k, "AP_ORDER", 8) == 0) {
	    sscanf(header+k+9, "%d", &aporder);
	 } else if(strncmp(header+k, "BP_ORDER", 8) == 0) {
	    sscanf(header+k+9, "%d", &bporder);

/* SIP coefficients */
	 } else if(sscanf(header+k, "A_%d_%d ", &i, &j) == 2) {
	    n = ((i+j)*(i+j+1))/2 + j;
	    sscanf(header+k+9, "%lf", &wcs->sipa[n]);
	    wcs->usipmax = MAX(wcs->usipmax, i);
	    wcs->vsipmax = MAX(wcs->vsipmax, j);
	 } else if(sscanf(header+k, "B_%d_%d ", &i, &j) == 2) {
	    n = ((i+j)*(i+j+1))/2 + j;
	    sscanf(header+k+9, "%lf", &wcs->sipb[n]);
	    wcs->usipmax = MAX(wcs->usipmax, i);
	    wcs->vsipmax = MAX(wcs->vsipmax, j);
	 } else if(sscanf(header+k, "AP_%d_%d ", &i, &j) == 2) {
	    n = ((i+j)*(i+j+1))/2 + j;
	    sscanf(header+k+9, "%lf", &wcs->sipap[n]);
	    wcs->usippmax = MAX(wcs->usippmax, i);
	    wcs->vsippmax = MAX(wcs->vsippmax, j);
	 } else if(sscanf(header+k, "BP_%d_%d ", &i, &j) == 2) {
	    n = ((i+j)*(i+j+1))/2 + j;
	    sscanf(header+k+9, "%lf", &wcs->sipbp[n]);
	    wcs->usippmax = MAX(wcs->usippmax, i);
	    wcs->vsippmax = MAX(wcs->vsippmax, j);

/* PV coefficients */
	 } else if(sscanf(header+k, "PV%d_%d ", &i, &j) == 2) {
	    sscanf(header+k+9, "%lf", &wcs->pv[i-1][j]);
	    pvdistort = 1;
	    wcs->npv = MAX(wcs->npv, j);
	    wcs->xpvmax = MAX(wcs->xpvmax, xpv_pwr[j]);
	    wcs->ypvmax = MAX(wcs->ypvmax, ypv_pwr[j]);
	    wcs->rpvmax = MAX(wcs->rpvmax, rpv_pwr[j]);

/* NOAO/IRAF style WAT (world attribute) cards */
/*
WAT0_001= 'system=image'       / Coordinate system
WAT1_001= 'wtype=zpx axtype=ra projp0=0. projp1=1. projp2=0. projp3=337.74 proj'
WAT1_002= 'p4=0. projp5=632052. lngcor = "3. 3. 3. 2. 0.001876397956622823 0.29'
WAT1_003= '99113930557312 0.1542460039112511 0.3032873851581314 1.9247409545894'
WAT1_004= '95E-5 -1.348328290485618E-5 1.414186703253352E-4 -1.792784764381400E'
WAT1_005= '-4 -1.276226238774833E-4 4.339217671825231E-4 "'
WAT2_001= 'wtype=zpx axtype=dec projp0=0. projp1=1. projp2=0. projp3=337.74 pro'
WAT2_002= 'jp4=0. projp5=632052. latcor = "3. 3. 3. 2. 0.001876397956622823 0.2'
WAT2_003= '999113930557312 0.1542460039112511 0.3032873851581314 9.963957331149'
WAT2_004= '402E-5 -1.378185066830135E-4 1.559892401479664E-4 -8.280442729203771'
WAT2_005= 'E-4 3.966701903249366E-4 0.001678960379199465 "'
*/
	 } else if(sscanf(header+k, "WAT%d_%d ", &i, &j) == 2) {
	    if(i < 1 || i > 2) continue;
	    q1 = index(header+k, '\'');
	    q2 = index(q1+1, '\'');
	    if(q1 == NULL || q2 == NULL || q2 <= q1) {
	       fprintf(stderr, "Malformed WAT string of length %ld `%.80s`\n", q2-q1, header+k);
	       exit(1);
	    }
	    strncat(watstring[i-1], q1+1, q2-q1-1);
	 }

/* Quit when we see the "END     " key */
	 if(strncmp(header+k, "END     ", 8) == 0) {
	    endseen = hdrcnt;
/* If naxis=0 it's either fpacked or some sort of table... keep going */
	    if(naxis > 0) {
	       nhead = endseen;
	       break;
	    }
	 }
	 hdrcnt++;
      }
      if(nhead) break;		/* nhead does *not* include END */
   }
   if(strcmp(fitsname, "-") != 0) close(fdin);

   if(fpack) {
      naxis = znaxis;
      for(i=1; i<10; i++) naxes[i] = znaxes[i];
   }
   wcs->nx = naxes[1];
   wcs->ny = naxes[2];

#ifdef DEBUG
   if(VERB > 0) {
      printf("nx= %d ny= %d\n", wcs->nx, wcs->ny);
   }
#endif

/* Patch things up if only CDELT and CROTA but no CD matrix */
   if(horrible_crota == 0x0007) {
      wcs->cd11 =                     cdelt1*cos(crota*dr);
      wcs->cd12 =  ABS(cdelt2/cdelt1)*cdelt1*sin(crota*dr);
      wcs->cd21 = -ABS(cdelt1/cdelt2)*cdelt2*sin(crota*dr);
      wcs->cd22 =                     cdelt2*cos(crota*dr);
#ifdef DEBUG
      if(VERB > 0) {
	 printf("Horrible %d\n", horrible_crota);
      }
#endif
   } else if( (horrible_crota&0x00f0) != 0x00f0) {
      fprintf(stderr, "Error: wcseval does not see all CD matrix terms %d\n",
	      horrible_crota);
      return(-1);
   }

/* Inverse of CD matrix */
   wcs->dc11 =  wcs->cd22 / (wcs->cd11*wcs->cd22-wcs->cd12*wcs->cd21);
   wcs->dc22 =  wcs->cd11 / (wcs->cd11*wcs->cd22-wcs->cd12*wcs->cd21);
   wcs->dc12 = -wcs->cd12 / (wcs->cd11*wcs->cd22-wcs->cd12*wcs->cd21);
   wcs->dc21 = -wcs->cd21 / (wcs->cd11*wcs->cd22-wcs->cd12*wcs->cd21);

/* Verify that the RA and DEC CTYPEs are the same */
   if(strcmp(wcs->ctype1+5, wcs->ctype2+5) != 0) {
      fprintf(stderr, "RA CTYPE `%s' and DEC CTYPE `%s' do not match\n",
	      wcs->ctype1, wcs->ctype2);
      return(-1);
   }

/* Is this an identifiable CTYPE? */
   wcs->projection = PROJ_UNKNOWN;
   if(strncmp(wcs->ctype1+5, "TPV", 3) == 0) {
      wcs->projection = PROJ_TAN_TPV;

   } else if(strncmp(wcs->ctype1+5, "TAN-SIP", 7) == 0) {
      wcs->projection = PROJ_TAN_SIP;

   } else if(strncmp(wcs->ctype1+5, "TAN", 3) == 0) {
      wcs->projection = PROJ_TAN;

   } else if(strncmp(wcs->ctype1+5, "STG", 3) == 0) {
      wcs->projection = PROJ_STG;

   } else if(strncmp(wcs->ctype1+5, "SIN", 3) == 0) {
      wcs->projection = PROJ_SIN;

   } else if(strncmp(wcs->ctype1+5, "ARC", 3) == 0) {
      wcs->projection = PROJ_ARC;

   } else if(strncmp(wcs->ctype1+5, "ZPN", 3) == 0) {
      wcs->projection = PROJ_ZPN;

   } else if(strncmp(wcs->ctype1+5, "ZPX", 3) == 0) {
      wcs->projection = PROJ_ZPX;

   } else if(strncmp(wcs->ctype1+5, "CAR", 3) == 0) {
      wcs->projection = PROJ_CD;

   } else if(strncmp(wcs->ctype1+5, "CD", 2) == 0) {
      wcs->projection = PROJ_CD;
   } 

/* Generic gnomonic test */
   gnomonic = wcs->projection == PROJ_TAN || 
              wcs->projection == PROJ_TAN_TPV || 
              wcs->projection == PROJ_TAN_SIP;

/* Problems if CTYPE is TAN-SIP, and there are no A_ORDER and B_ORDER */
   sipdistort = wcs->projection == PROJ_TAN_SIP && aorder > 0 && border > 0;

/* If ZPN move the r powers delivered in PV2_x terms to wcs->zpr[] */
   if(wcs->projection == PROJ_ZPN) {
      pvdistort = 0;
      for(i=0; i<=wcs->npv; i++) {
	 wcs->zpr[i] = wcs->pv[1][i];
	 if(i>1) wcs->pv[1][i] = 0.0;
      }
      wcs->zprmax = wcs->npv;
      wcs->npv = wcs->rpvmax = 0;
   }

/* If ZPX and WAT, try to figure it all out.  This is really retarded. */
   if(wcs->projection == PROJ_ZPX && 
      strlen(watstring[0]) > 0 && strlen(watstring[1]) > 0) {
      for(i=0; i<2; i++) {
	 if(strstr(watstring[i], "wtype=zpx") == NULL) {
	    fprintf(stderr, "wtype=zpx not found in ZPX WAT string\n");
	    exit(1);
	 }
      }

/* Pick up the radial terms "projpX".  Unclear how high IRAF will go,
 * example is 5 Assume the terms are the same for both axes, anything
 * else does not make sense. */
      for(j=0; j<10; j++) {
	 sprintf(fmt, "projp%d=", j);
	 if( (q1 = strstr(watstring[0], fmt)) != NULL) {
	    sprintf(fmt, "projp%d=%%lf", j);
	    sscanf(q1, fmt, &wcs->zpr[j]);
	    wcs->zprmax = MAX(wcs->zprmax, j);
	 }
//	 printf("%d %d %d %10.3f\n", i, j, nradterm, pradterm[j]);
      }

/* Get the longitude coefficient count */
      double pt, px, py, pf;
      int ptype=0, ptype2=0, npxi=0, npxi2=0, npeta=0, npeta2=0, pxterm=0, pxterm2=0;
      if( (q1 = strstr(watstring[0], "lngcor = \"")) != NULL) {
	 sscanf(q1, "lngcor = \" %lf %lf %lf %lf", &pt, &px, &py, &pf);
	 ptype  = pt+0.5;
	 npxi   = px+0.5;
	 npeta  = py+0.5;
	 pxterm = pf+0.5;
//	 printf("lngcor: %d %d %d %d\n", ptype, npxi, npeta, pxterm);
      }
/* Get the latitude coefficient count */
      if( (q1 = strstr(watstring[1], "latcor = \"")) != NULL) {
	 sscanf(q1, "latcor = \" %lf %lf %lf %lf", &pt, &px, &py, &pf);
	 ptype2  = pt+0.5;
	 npxi2   = px+0.5;
	 npeta2  = py+0.5;
	 pxterm2 = pf+0.5;
//	 printf("latcor: %d %d %d %d\n", ptype2, npxi2, npeta2, pxterm2);
      }
      if(ptype != ptype2 || npxi != npxi2 || 
	 npeta != npeta2 || pxterm != pxterm2) {
	 fprintf(stderr, "Inconsistent ZPX polys: %d %d  %d %d  %d %d  %d %d\n",
		 ptype, ptype2, npxi, npxi2, npeta, npeta2, pxterm, pxterm2);
	 exit(1);
      }
      if(ptype != 3) {
	 fprintf(stderr, "ZPX WAT calls for non-simple poly %d\n", ptype);
	 exit(1);
      }
      if(pxterm != 2) {
	 fprintf(stderr, "ZPX WAT calls for no or full cross terms %d\n", pxterm);
	 exit(1);
      }
      if(npxi != npeta) {
	 fprintf(stderr, "ZPX WAT has non-identical nxi, neta %d %d\n", npxi, npeta);
	 exit(1);
      }

      double xmin, xmax, ymin, ymax, xcoeff[16], ycoeff[16];
//      int ncoeff = (npxi+1)*(npxi+2)/2;
      int ncoeff = npxi*(npxi+1)/2;

/* Read the lng coefficients */
      q1 = strstr(watstring[0], "lngcor = \"");
      j = sscanf(q1, "lngcor = \"%lf %lf %lf %lf %lf %lf %lf %lf %lf %lf %lf %lf %lf %lf %lf %lf %lf %lf %lf %lf %lf %lf %lf %lf",
		 &pt, &px, &py, &pf, &xmin, &xmax, &ymin, &ymax, 
		 &xcoeff[0], &xcoeff[1], &xcoeff[2], &xcoeff[3], 
		 &xcoeff[4], &xcoeff[5], &xcoeff[6], &xcoeff[7], 
		 &xcoeff[8], &xcoeff[9], &xcoeff[10], &xcoeff[11], 
		 &xcoeff[12], &xcoeff[13], &xcoeff[14], &xcoeff[15]);
      if(j != 8+ncoeff) {
	 fprintf(stderr, "ZPX lng WAT has only %d instead of %d coeffs\n", j, 8+ncoeff);
	 exit(1);
      }

/* Read the lat coefficients */
      q1 = strstr(watstring[1], "latcor = \"");
      j = sscanf(q1, "latcor = \"%lf %lf %lf %lf %lf %lf %lf %lf %lf %lf %lf %lf %lf %lf %lf %lf %lf %lf %lf %lf %lf %lf %lf %lf",
		 &pt, &px, &py, &pf, &xmin, &xmax, &ymin, &ymax, 
		 &ycoeff[0], &ycoeff[1], &ycoeff[2], &ycoeff[3], 
		 &ycoeff[4], &ycoeff[5], &ycoeff[6], &ycoeff[7], 
		 &ycoeff[8], &ycoeff[9], &ycoeff[10], &ycoeff[11], 
		 &ycoeff[12], &ycoeff[13], &ycoeff[14], &ycoeff[15]);
      if(j != 8+ncoeff) {
	 fprintf(stderr, "ZPX lat WAT has only %d instead of %d coeffs\n", j, 8+ncoeff);
	 exit(1);
      }
//      printf("%8.3f %8.3f %8.3f %8.3f\n", xmin, xmax, ymin, ymax);
//      for(j=0; j<ncoeff; j++) printf(" %10.3e", xcoeff[j]);
//      printf("\n");
//      for(j=0; j<ncoeff; j++) printf(" %10.3e", ycoeff[j]);
//      printf("\n");

/* Distribute these parameters into the TAN PV coefficients */
/* Note that the min/max scaling params are ignored for a simple polynomial */

/* The x,y terms */
      for(j=n=0; j<npxi; j++) {
	 for(i=0; i<npxi-j; i++, n++) {
	    for(k=0; k<MAXPV && (xpv_pwr[k] != i || ypv_pwr[k] != j); k++);
	    wcs->pv[0][k] = xcoeff[n];
	    wcs->pv[1][k] = ycoeff[n];
	    wcs->npv = MAX(wcs->npv, k);
	    wcs->xpvmax = MAX(wcs->xpvmax, xpv_pwr[k]);
	    wcs->ypvmax = MAX(wcs->ypvmax, ypv_pwr[k]);
	 }
      }
   }

#ifdef DEBUG
   if(VERB > 0) {
      printf("CRPIX= %7.2f %7.2f   CRVAL= %9.4f %9.4f\n", 
	     wcs->crpix1, wcs->crpix2, wcs->crval1, wcs->crval2);
      printf("CD= %9.6f %9.6f %9.6f %9.6f\n", wcs->cd11, wcs->cd12, wcs->cd21, wcs->cd22);
      printf("CTYPE= %s %s  Projection %d\n", wcs->ctype1, wcs->ctype2, wcs->projection);

      printf("gnomonic= %d   sip= %d   pv= %d\n", gnomonic, sipdistort, pvdistort);
      printf("NSIP= %d %d %d %d\n", aorder, border, aporder, bporder);
      for(i=0; i<9; i++) {
	 for(j=0; i+j<9; j++) {
	    n = ((i+j)*((i+j)+1))/2 + j;
	    if(wcs->sipa[n] != 0.0 || wcs->sipb[n] != 0.0 ||
	       wcs->sipap[n] != 0.0 || wcs->sipbp[n] != 0.0) {
	       printf("%2d  A_%d_%d %12.4e  B_%d_%d %12.4e  AP_%d_%d %12.4e  BP_%d_%d %12.4e\n", 
		      n, i,j, wcs->sipa[n], i,j, wcs->sipb[n], i,j, wcs->sipap[n], i,j, wcs->sipbp[n]);
	    }
	 }
      }

      printf("npv  %d  x,y,rmax= %d %d %d\n", 
	     wcs->npv, wcs->xpvmax, wcs->ypvmax, wcs->rpvmax);
      for(j=0; j<MAXPV; j++) {
	 if(wcs->pv[0][j] != 0.0 || wcs->pv[1][j] != 0.0)
	    printf("PV1_%d\t %12.4e  x^%d y^%d r^%d       PV2_%d\t %12.4e  x^%d y^%d r^%d\n", j, wcs->pv[0][j], xpv_pwr[j], ypv_pwr[j], rpv_pwr[j], j, wcs->pv[1][j], xpv_pwr[j], ypv_pwr[j], rpv_pwr[j]);
      }

      printf("nZPN  %d\n", wcs->zprmax);
      for(j=0; j<=wcs->zprmax; j++) {
	 if(wcs->zpr[j] != 0.0)
	    printf("ZPN_%d\t %12.4e  r^%d\n", j, wcs->zpr[j], j);
      }
   }
//   if(strlen(watstring[0]) > 0) printf("WAT1= %s\n", watstring[0]);
//   if(strlen(watstring[1]) > 0) printf("WAT2= %s\n", watstring[1]);


#endif

/* Decide how to map */
   if(wcs->projection == PROJ_UNKNOWN) {
      fprintf(stderr, "Error: libwcs sorry, unknown map CTYPE %s %s, exiting...\n", wcs->ctype1, wcs->ctype2);
      exit(1);
   }

/* Definite call for PV distortion? */
   wcs->pvok = wcs->sipok = 0;
   if( wcs->projection == PROJ_TAN_TPV && pvdistort) {
      wcs->pvok = 1;

/* Definite call for SIP distortion? */
   } else if(sipdistort && aorder > 0 && border > 0) {
      wcs->sipok = 1;

/* Probable PV? */
   } else if(pvdistort) {
      wcs->pvok = 1;

/* Probable SIP? */
   } else if(aorder > 0 && border > 0 && aporder > 0 && bporder > 0) {
      wcs->sipok = 1;

   }
   
   return(0);
}

////////////////////////////////////////////////////////////////
////////////////////////////////////////////////////////////////

/* Convert sky coords to gnomonic tangent plane, all in deg */
void sky2tp_tan(double ra, double dec, double a0, double d0,
		double *xi, double *eta)
{
   double dr=atan(1.0)/45.0, Xp, Stinv;
   Xp = cos(dec*dr) * cos((ra-a0)*dr);
   Stinv = 1 / (sin(d0*dr) * sin(dec*dr) + Xp*cos(d0*dr));
   *xi  = Stinv * cos(dec*dr) * sin((ra-a0)*dr) / dr;
   *eta = Stinv * (cos(d0*dr) * sin(dec*dr) - Xp*sin(d0*dr)) / dr;
   return;
}

/* Convert gnomonic tangent plane to sky coords, all in deg */
void tp2sky_tan(double xi, double eta, double a0, double d0,
		double *ra, double *dec)
{
   double dr=atan(1.0)/45.0, theta, phi, X, Y;
   theta = atan(dr*sqrt(xi*xi+eta*eta));
   phi = atan2(xi, -eta);
   X = sin(d0*dr)*sin(theta)*cos(phi) + cos(d0*dr)*cos(theta);
   Y = sin(theta)*sin(phi);
   *ra = fmod(a0 + atan2(Y, X)/dr + 4*360, 360.0);
   *dec = asin(sin(d0*dr)*cos(theta) - cos(d0*dr)*sin(theta)*cos(phi)) / dr;
   return;
}

/* Convert sky coords to stereographic tangent plane, all in deg */
void sky2tp_stg(double ra, double dec, double a0, double d0,
	    double *xi, double *eta)
{
   double dr=atan(1.0)/45.0, Xp, Stg/* 1/(1+sin(t) */;
   Xp = cos(dec*dr) * cos((ra-a0)*dr);
   Stg = 1 / (1 + sin(d0*dr) * sin(dec*dr) + Xp*cos(d0*dr));
   *xi  = 2 * Stg * cos(dec*dr) * sin((ra-a0)*dr) / dr;
   *eta = 2 * Stg * (cos(d0*dr) * sin(dec*dr) - Xp*sin(d0*dr)) / dr;
   return;
}

/* Convert stereographic tangent plane to sky coords, all in deg */
void tp2sky_stg(double xi, double eta, double a0, double d0,
	    double *ra, double *dec)
{
   double dr=atan(1.0)/45.0, theta, phi, X, Y;
   theta = 2*atan(dr*sqrt(xi*xi+eta*eta)/2);
   phi = atan2(xi, -eta);
   X = sin(d0*dr)*sin(theta)*cos(phi) + cos(d0*dr)*cos(theta);
   Y = sin(theta)*sin(phi);
   *ra = fmod(a0 + atan2(Y, X)/dr + 4*360, 360.0);
   *dec = asin(sin(d0*dr)*cos(theta) - cos(d0*dr)*sin(theta)*cos(phi)) / dr;
   return;
}

/* Convert sky coords to orthographic tangent plane, all in deg */
void sky2tp_sin(double ra, double dec, double a0, double d0,
	    double *xi, double *eta)
{
   double dr=atan(1.0)/45.0, Xp;
   Xp = cos(dec*dr) * cos((ra-a0)*dr);
   *xi  = cos(dec*dr) * sin((ra-a0)*dr) / dr;
   *eta = (cos(d0*dr) * sin(dec*dr) - Xp*sin(d0*dr)) / dr;
   return;
}

/* Convert orthographic tangent plane to sky coords, all in deg */
void tp2sky_sin(double xi, double eta, double a0, double d0,
	    double *ra, double *dec)
{
   double dr=atan(1.0)/45.0, theta, phi, X, Y;
   theta = asin(dr*sqrt(xi*xi+eta*eta));
   phi = atan2(xi, -eta);
   X = sin(d0*dr)*sin(theta)*cos(phi) + cos(d0*dr)*cos(theta);
   Y = sin(theta)*sin(phi);
//   printf("\n%7.3f %7.3f %7.3f %7.3f\n", theta, phi, X, Y);
   *ra = fmod(a0 + atan2(Y, X)/dr + 4*360, 360.0);
   *dec = asin(sin(d0*dr)*cos(theta) - cos(d0*dr)*sin(theta)*cos(phi)) / dr;
   return;
}

/* Convert sky coords to arc, all in deg */
void sky2tp_arc(double ra, double dec, double a0, double d0,
	    double *xi, double *eta)
{
   double dr=atan(1.0)/45.0, X, Y, phi, theta;
   X = cos(dec*dr)*sin((ra-a0)*dr);
   Y = (cos(d0*dr)*sin(dec*dr) - cos(dec*dr)*sin(d0*dr)*cos((ra-a0)*dr));
   phi = atan2(-X, Y);
   theta = acos(sin(dec*dr)*sin(d0*dr) + cos(dec*dr)*cos(d0*dr)*cos((ra-a0)*dr));
   *xi = -theta/dr * sin(phi);
   *eta = theta/dr * cos(phi);
   return;
}

/* Convert arc to sky coords, all in deg */
void tp2sky_arc(double xi, double eta, double a0, double d0,
	    double *ra, double *dec)
{
   double dr=atan(1.0)/45.0, theta, phi, X, Y;
   theta = dr * sqrt(xi*xi+eta*eta);
   phi = atan2(xi, -eta);
   X = sin(d0*dr)*sin(theta)*cos(phi) + cos(d0*dr)*cos(theta);
   Y = sin(theta)*sin(phi);
   *ra = fmod(a0 + atan2(Y, X)/dr + 4*360, 360.0);
   *dec = asin(sin(d0*dr)*cos(theta) - cos(d0*dr)*sin(theta)*cos(phi)) / dr;
   if(VERB > 1) {
      printf("%8.3f %8.3f %8.3f %8.3f\n", theta, phi, X, Y);
   }

   return;
}

////////////////////////////////////////////////////////////////
/// (LEGACY) ///

/* Convert sky coords to gnomonic tangent plane, all in deg */
void sky2tp(double ra, double dec, double a0, double d0,
		double *xi, double *eta)
{
   sky2tp_tan(ra, dec, a0, d0, xi, eta);
   return;
}

/* Convert gnomonic tangent plane to sky coords, all in deg */
void tp2sky(double xi, double eta, double a0, double d0,
		double *ra, double *dec)
{
   tp2sky_tan(xi, eta, a0, d0, ra, dec);
   return;
}

////////////////////////////////////////////////////////////////



#ifdef MAINTEST
int main(int argc, char **argv)
{
   int i, j;
   double xi, eta, a0, d0, ra, dec;

   a0 = 105;
   d0 = 65;
   a0 = 205;
   d0 = -65;
//   a0 = 90;
//   d0 = 90;

   for(j=-60; j<=60; j+=10) {
      for(i=-60; i<=60; i+=10) {
	 printf("%4d %4d", i, j);

	 tp2sky_tan((double)i, (double)j, a0, d0, &ra, &dec);
	 sky2tp_tan(ra, dec, a0, d0, &xi, &eta);
	 printf(" %s: %7.3f %7.3f %7.3f %7.3f", "tan", ra, dec, xi-i, eta-j);

	 tp2sky_stg((double)i, (double)j, a0, d0, &ra, &dec);
	 sky2tp_stg(ra, dec, a0, d0, &xi, &eta);
	 printf(" %s: %7.3f %7.3f %7.3f %7.3f", "stg", ra, dec, xi-i, eta-j);

	 tp2sky_sin((double)i, (double)j, a0, d0, &ra, &dec);
	 sky2tp_sin(ra, dec, a0, d0, &xi, &eta);
	 printf(" %s: %7.3f %7.3f %7.3f %7.3f", "sin", ra, dec, xi-i, eta-j);

	 tp2sky_arc((double)i, (double)j, a0, d0, &ra, &dec);
	 sky2tp_arc(ra, dec, a0, d0, &xi, &eta);
	 printf(" %s: %7.3f %7.3f %7.3f %7.3f", "arc", ra, dec, xi-i, eta-j);

	 printf("\n");
      }
   }
   return(0);
}
#endif



#if 0
/* Apply a TPV polynomial */
void apply_tpv_orig(WCS *wcs, double *xi, double *eta)
{
   int j;
   double upow[10], vpow[10], rpow[10], r;

/* Initialize power terms */
   upow[0] = vpow[0] = rpow[0] = 1.0;
   for(j=1; j<=wcs->xpvmax; j++) upow[j] = upow[j-1] * (*xi);
   for(j=1; j<=wcs->ypvmax; j++) vpow[j] = vpow[j-1] * (*eta);
   if(wcs->rpvmax > 0) {
      r = sqrt((*xi)*(*xi) + (*eta)*(*eta));
      for(j=1; j<=wcs->rpvmax; j++) rpow[j] = rpow[j-1] * r;
   }
/* Apply the polynomial */
   *xi = *eta = 0.0;
   for(j=0; j<=wcs->npv; j++) {
      *xi  += wcs->pv[0][j] * upow[xpv_pwr[j]] * vpow[ypv_pwr[j]] * rpow[rpv_pwr[j]];
      *eta += wcs->pv[1][j] * upow[ypv_pwr[j]] * vpow[xpv_pwr[j]] * rpow[rpv_pwr[j]];
#ifdef DEBUG
      if(VERB > 1) {
	 printf("%9.5f %9.5f  PV1_%d %10.3e   PV2_%d %10.3e  %d %d %d\n", 
		*xi, *eta, j, wcs->pv[0][j], j, wcs->pv[1][j], xpv_pwr[j], ypv_pwr[j], rpv_pwr[j]);
      }
#endif
   }
}

/* Apply an inverse TPV polynomial by reversion */
void apply_tpv_inverse_revert(int npv, int nx, int ny, int nr, 
	       double *xpv, double *ypv, double *xi, double *eta)
{
   int i, j;
   double upow[10], vpow[10], rpow[10], r, pv01, pv11, xp, yp, dx, dy;

/* Initialize power terms */
   upow[0] = vpow[0] = rpow[0] = 1.0;

/* Remove a linear term and set aside the constant */
   pv01 = xpv[1];
   pv11 = ypv[1];
   xpv[1] -= 1.0;
   ypv[1] -= 1.0;
   xp = *xi;
   yp = *eta;
   for(i=0; i<10; i++) {
      dx = xp;
      dy = yp;
      for(j=1; j<=nx; j++) upow[j] = upow[j-1] * xp;
      for(j=1; j<=ny; j++) vpow[j] = vpow[j-1] * yp;
      if(nr > 0) {
	 r = sqrt(xp*xp + yp*yp);
	 for(j=1; j<=nr; j++) rpow[j] = rpow[j-1] * r;
      }
      xp = *xi;
      yp = *eta;
      for(j=1; j<=npv; j++) {
	 xp -= xpv[j] * upow[xpv_pwr[j]] * vpow[ypv_pwr[j]] * rpow[rpv_pwr[j]];
	 yp -= ypv[j] * upow[ypv_pwr[j]] * vpow[xpv_pwr[j]] * rpow[rpv_pwr[j]];
#ifdef DEBUG
	 if(VERB > 2) {
	    printf("%9.5f %9.5f %9.5f %9.5f  PV1_%d %10.3e   PV2_%d %10.3e  %d %d %d\n", 
		   *xi, *eta, xp, yp, j, xpv[j], j, ypv[j], xpv_pwr[j], ypv_pwr[j], rpv_pwr[j]);
	 }
#endif
      }
      dx -= xp;
      dy -= yp;
#ifdef DEBUG
      if(VERB > 1) {
	 printf("%9.5f %9.5f  %9.5f %9.5f %10.3e %10.3e\n", 
		*xi, *eta, xp, yp, dx, dy);
      }
#endif
      if(ABS(dx)<TOL_DOUBLE && ABS(dy)<TOL_DOUBLE) break;
   }
/* Restore linear term and apply the constant */
   xpv[1] = pv01;
   ypv[1] = pv11;
   *xi  = xp - xpv[0];
   *eta = yp - ypv[0];
}

/* Apply an inverse TPV polynomial */
void apply_tpv_inverse_orig(WCS *wcs, double *xi, double *eta)
{
   int i, j;
   double upow[10], vpow[10], rpow[10], r, pv01, pv11, xp, yp, dx, dy;

/* Initialize power terms */
   upow[0] = vpow[0] = rpow[0] = 1.0;

/* Remove a linear term and set aside the constant */
   pv01 = wcs->pv[0][1];
   pv11 = wcs->pv[1][1];
   wcs->pv[0][1] -= 1.0;
   wcs->pv[1][1] -= 1.0;
   xp = *xi;
   yp = *eta;
   for(i=0; i<10; i++) {
      dx = xp;
      dy = yp;
      for(j=1; j<=wcs->xpvmax; j++) upow[j] = upow[j-1] * xp;
      for(j=1; j<=wcs->ypvmax; j++) vpow[j] = vpow[j-1] * yp;
      if(wcs->rpvmax > 0) {
	 r = sqrt(xp*xp + yp*yp);
	 for(j=1; j<=wcs->rpvmax; j++) rpow[j] = rpow[j-1] * r;
      }
      xp = *xi;
      yp = *eta;
      for(j=1; j<=wcs->npv; j++) {
	 xp -= wcs->pv[0][j] * upow[xpv_pwr[j]] * vpow[ypv_pwr[j]] * rpow[rpv_pwr[j]];
	 yp -= wcs->pv[1][j] * upow[ypv_pwr[j]] * vpow[xpv_pwr[j]] * rpow[rpv_pwr[j]];
#ifdef DEBUG
	 if(VERB > 2) {
	    printf("%9.5f %9.5f %9.5f %9.5f  PV1_%d %10.3e   PV2_%d %10.3e  %d %d %d\n", 
		   *xi, *eta, xp, yp, j, wcs->pv[0][j], j, wcs->pv[1][j], xpv_pwr[j], ypv_pwr[j], rpv_pwr[j]);
	 }
#endif
      }
      dx -= xp;
      dy -= yp;
#ifdef DEBUG
      if(VERB > 1) {
	 printf("%9.5f %9.5f  %9.5f %9.5f %10.3e %10.3e\n", 
		*xi, *eta, xp, yp, dx, dy);
      }
#endif
      if(ABS(dx)<wcs->TOL && ABS(dy)<wcs->TOL) break;
   }
/* Restore linear term and apply the constant */
   wcs->pv[0][1] = pv01;
   wcs->pv[1][1] = pv11;
   *xi  = xp - wcs->pv[0][0];
   *eta = yp - wcs->pv[1][0];
}
#endif

#if 0
/* Apply a ZPX distortion polynomial and ZPN radial distortion: sky->pix */
void apply_zpx_orig(WCS *wcs, double *xi, double *eta)
{
   int j, nr;
   double upow[10], vpow[10], rpow[10], r, R, dr=atan(1.0)/45;
   double coeff[10];

// FIXME: this is a stupid indexing into the PV list

   printf("%10.6f %10.6f\n", *xi, *eta);

/* Initialize power terms.  Radial terms are ignored! */
   for(j=1, upow[0]=1.0; j<=wcs->xpvmax; j++) upow[j] = upow[j-1] * (*xi);
   for(j=1, vpow[0]=1.0; j<=wcs->ypvmax; j++) vpow[j] = vpow[j-1] * (*eta);
/* Add the polynomial to the coords */
   for(j=0; j<=wcs->npv; j++) {
      if(wcs->pv[0][j] == 0 && wcs->pv[1][j] == 0) continue;
      if(rpv_pwr[j] != 0) continue;
      *xi  += wcs->pv[0][j] * upow[xpv_pwr[j]] * vpow[ypv_pwr[j]];
      *eta += wcs->pv[1][j] * upow[ypv_pwr[j]] * vpow[xpv_pwr[j]];
	 printf("%9.5f %9.5f  PV1_%d %10.3e   PV2_%d %10.3e  %d %d\n", 
		*xi, *eta, j, wcs->pv[0][j], j, wcs->pv[1][j], xpv_pwr[j], ypv_pwr[j]);

#ifdef DEBUG
      if(VERB > 1) {
	 printf("%9.5f %9.5f  PV1_%d %10.3e   PV2_%d %10.3e  %d %d\n", 
		*xi, *eta, j, wcs->pv[0][j], j, wcs->pv[1][j], xpv_pwr[j], ypv_pwr[j]);
      }
#endif
   }

/* Apply the inverse ZPN distortion using the radial terms */
   for(j=0, nr=0, R=0.0; j<=wcs->npv; j++) {
      if(rpv_pwr[j] == 0) continue;
      coeff[rpv_pwr[j]] = wcs->pv[0][j];
      nr = MAX(nr, rpv_pwr[j]);
   }
   apply_zpn_inverse(nr, coeff, xi, eta);
}
#endif

