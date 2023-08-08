/* Header file to convert between pixel and celestial coordinates */

/*
 * Follow the conventions of Calabretta and Greisen
 * see http://fits.gsfc.nasa.gov/fits_wcs.html
 *
 *   Greisen, E.W., and Calabretta, M.R., A&A, 395, 1061-1075, 2002.
 *   Calabretta, M.R., and Greisen, E.W., A&A, 395, 1077-1122, 2002.
 *
 *  - pixel coords p1,p2 are defined such that the center of LL is 1.0,1.0
 *  - focal plane relative coords are u,v[pix] = p{1,2} - CRPIX{1,2}
 *  - intermediate WCS coords are x,y[deg] = [CD{1,2}_{1,2}] u,v
 *    ("x,y", is the C+G nomenclature for "xi,eta").
 *  - x,y is a *left-handed* coordinate system, with the intention
 *    that RA ~ RA0 + x, Dec ~ Dec0 + y (hold-over from AIPS?)
 *
 *  - CDELT{1,2} may also be present, CROTA is severely deprecated.
 *  - CUNIT{1,2} may be present as units of CD, CDELT, CRVAL
 *  - CNPIX{1,2} is a Vista'ism that should be ignored except by Vista
 *
 * Spherical projections are accomplished deprojecting the
 * "intermediate WCS (tangent plane) coordinates x,y" to a "native
 * coordinate system of longitude and latitude: phi,theta".  These are
 * then subjected to a spherical rotation to RA,Dec using three Euler
 * angles.  Two are provided by CRVAL{1,2} that specify the location
 * in RA,Dec of the pixel location CRPIX{1,2}.
 *
 * Be clear about a couple of peculiar conventions adopted by C+G!
 *  a) C+G use "theta" to mean latitude, angle up from equator
 *  b) C+G apparently like to define phi=0 at -Dec in analogy to 
 *     some sort of backward, southern cartographic azimuth instead
 *     of mathematical polar angle.  phi increases clockwise.
 * See Calabretta and Greisen Fig 3.
 *
 * The location of CRPIX{1,2} in this "native coordinate system" is
 * (phi0,theta0).  For zenithal (plane) projections it is the north
 * pole: (phi0,theta0) = (0,90).  In order to make things more general
 * (confusing) and able to accomodate arcane projections such as
 * cylindrical and conical, the standard also introduces the location
 * of the pole of the "native coordinate system" (phip,thetap), known
 * as {LON,LAT}POLE.  This is the coordinate of the native pole in
 * either sky or native coords.  Since there is a degeneracy between
 * phi0 and phip, the standard exploits this to implement the third
 * Euler angle rotation.  For zenithal projections the default for
 * LONPOLE is 0 for delta0>=theta0 or 180 for delta0<theta0.  Since
 * theta0 is normally 90 deg this means LONPOLE's default is 180.  The
 * default for LATPOLE is 90 deg for zenithal projects; it is normally
 * not given.
 *
 * Of course the third Euler angle can be provided by the CD matrix,
 * e.g. rotating pixel p2 to a direction that will become Dec north.
 *
 * For zenithal projections, (phi0,theta0)=(0,90), and (phip,thetap)=(180,90)
 *      
 *      x =  R(theta) sin(phi)
 *      y = -R(theta) cos(phi)
 *
 *      phi      = atan2(x,-y)
 *      R(theta) = sqrt(x*x+y*y)
 *
 * Spherical rotations:
 *
 *      a,d      = RA,Dec 
 *      p,t      = phi,theta
 *
 *      pp,tp = phi,theta of celestial pole = LONPOLE,LATPOLE [180,90]
 *      ap,dp = RA,Dec of native pole.  For zenithal ap,dp = a0,d0 = CRVAL{1,2}
 *
 *      a = a0 + atan2(-cos(t)sin(p-pp), sin(t)cos(d0)-cos(t)sin(d0)cos(p-pp))
 *      d = asin(sin(t)sin(d0) + cos(t)cos(d0)cos(p-pp))
 *
 *      p = pp + atan2(-cos(d)sin(a-a0), sin(d)cos(d0)-cos(d)sin(d0)cos(a-a0))
 *      t = asin(sin(d)sin(d0) + cos(d)cos(d0)cos(a-a0))
 *
 *      l'     r11 r12 r13   l
 *      m'  =  r21 r22 r23   m
 *      n'     r31 r32 r33   n
 *
 *      l,m,n    = cos(d)cos(a), cos(d)sin(a), sin(d)
 *      l',m',n' = cos(t)cos(p), cos(t)sin(p), sin(t)
 *
 *      r11= -sin(ap)sin(pp) - cos(ap)cos(pp)sin(dp)
 *      r12=  cos(ap)sin(pp) - sin(ap)cos(pp)sin(dp)
 *      r13=                          cos(pp)cos(dp)
 *      r21=  sin(ap)cos(pp) - cos(ap)sin(pp)sin(dp)
 *      r22= -cos(ap)cos(pp) - sin(ap)sin(pp)sin(dp)
 *      r23=                          sin(pp)cos(dp)
 *      r31=                   cos(ap)       cos(dp)
 *      r32=                   sin(ap)       cos(dp)
 *      r33=                                 sin(dp)
 *
 * For the zenithal, plane projection case of pp=LONPOLE=180; ap=a0 dp=d0
 *
 *      a = a0 + atan2(cos(t)sin(p), sin(t)cos(d0)+cos(t)sin(d0)cos(p))
 *      d = asin(sin(t)sin(d0) - cos(t)cos(d0)cos(p))
 *
 *      p = 180 + atan2(-cos(d)sin(a-a0), sin(d)cos(d0)-cos(d)sin(d0)cos(a-a0))
 *      t = asin(sin(d)sin(d0) + cos(d)cos(d0)cos(a-a0))
 *
 *      cos(d) cos(a-a0) = sin(t)cos(d0) + cos(t)sin(d0)cos(p)
 *      cos(d) sin(a-a0) = cos(t)sin(p)
 *
 *      cos(t) cos(p) = -sin(d)cos(d0) + cos(d)sin(d0)cos(a-a0)
 *      cos(t) sin(p) =  cos(d)sin(a-a0)
 *
 *      r11=  cos(a0)sin(d0)
 *      r12=  sin(a0)sin(d0)
 *      r13=        -cos(d0)
 *      r21= -sin(a0)
 *      r22=  cos(a0)
 *      r23=  0
 *      r31=  cos(a0)cos(d0)
 *      r32=  sin(a0)cos(d0)
 *      r33=         sin(d0)
 *
 *
 * Various azimuthal projections have different R(theta) functions.
 *
 * TAN = Gnomonic (project from center)
 *
 *   R(theta) = 180/pi cot(theta)
 *   theta = atan(180/pi/R(theta))
 *
 * TAN-SIP
 *
 *   x,y = CD u+f(u,v),v+g(u,v)
 *
 *     f(u,v) = Sum A_p_q u^p v^q  p+q<A_ORDER
 *     g(u,v) = Sum B_p_q u^p v^q  p+q<B_ORDER
 *
 * TAN-TPV
 *
 *   x',y' = CD u,v
 *   x = f(x',y')
 *   y = g(x',y')
 *
 *   f(x',y') = Sum PV1_m x'^i y'^j r'^k  m=0,39
 *   g(x',y') = Sum PV2_m y'^i x'^j r'^k  i(m), j(m), k(m)...
 *
 * STG = Stereographic (project from antipode)
 *
 *   R(theta) = 180/pi 2*cos(theta)/(1+sin(theta))
 *            = 180/pi 2*tan((90-theta)/2)
 *   theta    = 90 - 2*atan(pi/180*R(theta)/2)
 *
 * SIN = Orthographic (project from infinity)
 *
 *   R(theta) = 180/pi cos(theta)
 *     theta  = acos(pi/180 R(theta))
 *
 * ARC = Arc from pole
 *
 *   R(theta) = 90-theta
 *   theta = 90-R(theta)
 *
 * ZPN = Arc from pole with polynomial terms for R(theta)
 *
 *   R(theta) = 180/pi Sum_0^20 PV2_m (pi/180*(90-theta))^m
 *     theta  = non-analytic inverse for R(theta)
 */

/* Known projections */
#define PROJ_UNKNOWN   0	/* Unsupported projection */
#define PROJ_CD      666	/* Just the CD matrix */
#define PROJ_TAN     100	/* Gnomonic projection */
#define PROJ_TAN_TPV 101	/* Gnomonic projection with PV terms */
#define PROJ_TAN_SIP 102	/* Gnomonic projection with SIP terms */
#define PROJ_STG     200	/* Stereographic projection */
#define PROJ_SIN     300	/* Orthographic projection */
#define PROJ_ARC     400	/* Arc projection */
#define PROJ_ZPN     401	/* Arc projection with theta polynomial */
#define PROJ_ZPX     402	/* Reg but non-std IRAF arc with x,y poly */

/* Max number of correction terms */
#define MAXSIP 55		/* Max SIP terms 9th order */
#define MAXPV 40		/* Max TPV terms 7th order */
#define MAXZPN 20		/* Max ZPN terms 20th order */

typedef struct wcs {
   int projection;		/* Projection type */
   char ctype1[80];		/* Axis 1 type */
   char ctype2[80];		/* Axis 2 type */
   double crpix1;		/* Coord reference pixel axis 1 */
   double crpix2;		/* Coord reference pixel axis 2 */
   double crval1;		/* Coord reference sky axis 1 */
   double crval2;		/* Coord reference sky axis 2 */
   double cd11;			/* CD matrix */
   double cd12;			/* CD matrix */
   double cd21;			/* CD matrix */
   double cd22;			/* CD matrix */
   double dc11;			/* Inverse CD matrix */
   double dc12;			/* Inverse CD matrix */
   double dc21;			/* Inverse CD matrix */
   double dc22;			/* Inverse CD matrix */
   double TOL;			/* Convergence tolerance [deg] */
   int nx;			/* Image width */
   int ny;			/* Image height */
   int usipmax;			/* Max order of u seen in SIP terms */
   int vsipmax;			/* Max order of v seen in SIP terms */
   int usippmax;		/* Max order of u seen in inverse SIP terms */
   int vsippmax;		/* Max order of v seen in inverse SIP terms */
   double sipa[MAXSIP];		/* SIP A coefficients (up to 9th order) */
   double sipb[MAXSIP];		/* SIP B coefficients */
   double sipap[MAXSIP];	/* SIP A inverse coefficients */
   double sipbp[MAXSIP];	/* SIP B inverse coefficients */
   int npv;			/* Max PV terms seen */
   int xpvmax;			/* Max order of x seen in PV terms */
   int ypvmax;			/* Max order of y seen in PV terms */
   int rpvmax;			/* Max order of r seen in PV terms */
   double pv[2][MAXPV];		/* PV coefficients (up to 7th order) */
   int zprmax;			/* Max order of r seen in ZPN terms */
   double zpr[MAXZPN];		/* ZPN radial coefficients */
   int pvok;			/* Legitimate TPV appears to be present */
   int sipok;			/* Legitimate SIP appears to be present */
   double u;			/* Focal plane relative coord */
   double v;			/* Focal plane relative coord */
   double xi;			/* Tangent plane relative coord */
   double eta;			/* Tangent plane relative coord */
   char *fname;			/* File from which WCS was obtained */
} WCS;

/* Read the WCS terms from a FITS file's header */
int getwcs(char *fitsname, WCS *wcs);

/* Initialize a WCS to a vanilla projection */
void wcs_init(
   double a0, double d0,        // [deg] center point on sphere
   double cx, double cy,        // [pix] center point on detector
   double scale,                // [deg/pix] scale
   int proj,                    // projection type from wcseval.h
   WCS *wcs);                   // wcs structure

/* Convert x,y [pix] to RA,Dec [deg] */
int xy2sky(WCS *wcs, double x, double y, double *ra, double *dec);

/* Convert RA,Dec [deg] to x,y [pix] */
int sky2xy(WCS *wcs, double x, double y, double *ra, double *dec);

/* Convert RA,Dec to x,y using SIP distortion terms and gnonomic */
int sky2xy_sip(WCS *wcs, double ra, double dec, double *x, double *y);
/* Convert RA,Dec to x,y using PV distortion terms and gnonomic */
int sky2xy_tpv(WCS *wcs, double ra, double dec, double *x, double *y);
/* Convert x,y to RA,Dec using SIP distortion terms and gnonomic */
int xy2sky_sip(WCS *wcs, double x, double y, double *ra, double *dec);
/* Convert x,y to RA,Dec using PV distortion terms and gnonomic */
int xy2sky_tpv(WCS *wcs, double x, double y, double *ra, double *dec);
/* Convert x,y to RA,Dec using no distortion terms and gnonomic */
int xy2sky_tan(WCS *wcs, double x, double y, double *ra, double *dec);
/* Convert RA,Dec to x,y using no distortion terms and gnonomic */
int sky2xy_tan(WCS *wcs, double ra, double dec, double *x, double *y);


/* Convert sky coords to gnomonic tangent plane, all in deg */
void sky2tp_tan(double ra, double dec, double a0, double d0,
		double *xi, double *eta);
/* Convert gnomonic tangent plane to sky coords, all in deg */
void tp2sky_tan(double xi, double eta, double a0, double d0,
		double *ra, double *dec);

/* Convert sky coords to stereographic tangent plane, all in deg */
void sky2tp_stg(double ra, double dec, double a0, double d0,
		double *xi, double *eta);
/* Convert stereographic tangent plane to sky coords, all in deg */
void tp2sky_stg(double xi, double eta, double a0, double d0,
		double *ra, double *dec);

/* Convert sky coords to orthographic tangent plane, all in deg */
void sky2tp_sin(double ra, double dec, double a0, double d0,
		double *xi, double *eta);
/* Convert orthographic tangent plane to sky coords, all in deg */
void tp2sky_sin(double xi, double eta, double a0, double d0,
		double *ra, double *dec);

/* Convert sky coords to arc, all in deg */
void sky2tp_arc(double ra, double dec, double a0, double d0,
		double *xi, double *eta);
/* Convert arc to sky coords, all in deg */
void tp2sky_arc(double xi, double eta, double a0, double d0,
		double *ra, double *dec);

/* LEGACY */

/* Convert sky coords to gnomonic tangent plane, all in deg */
void sky2tp(double ra, double dec, double a0, double d0,
	    double *xi, double *eta);
/* Convert gnomonic tangent plane to sky coords, all in deg */
void tp2sky(double xi, double eta, double a0, double d0,
	    double *ra, double *dec);

