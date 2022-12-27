#include <stdio.h>
#include <math.h>
#define RAD_PER_DEG	0.0174532925199432955

/*
#define B2A_SQ	0.99327733 
#define B2A_SQ	0.9931177
*/
#define B2A_SQ	0.993305521
/* 
  convert geographic latitude (degree) to geocentric latitude (degree)
*/
double geog_to_geoc__(double *lat)
{
	return(atan(B2A_SQ*tan(RAD_PER_DEG*(*lat)))/RAD_PER_DEG);
}

/*
  convert geocentric latitude (degree) to geographic latitude (degree)
*/
double geoc_to_geog__(double *lat)
{
	return(atan(tan(RAD_PER_DEG*(*lat))/B2A_SQ)/RAD_PER_DEG);
}

/*
  calculate distance, azimuth between source (c0,l0) and receiver (c1,l1).
  colatitude:[0,180] (geographic), longitude [0,360].
  distance <=180.
  azimuth [0,360), measured at source location clockwise from north.
*/
void distaz_(double *c0_f,double *l0_f,double *c1_f,double *l1_f,double *dist,double *az)
{
	double colat0,lon0,colat1,lon1,del,azimuth;
	double cosdel,cosaz,sinaz,tmp;
	double c0,l0,c1,l1;
	double eps=1.0e-8;
//	float geographic_to_geocentric(float lat);

	/*convert geographic to geocentric*/
//	c0=90.-geographic_to_geocentric(90.-c0);
//	c1=90.-geographic_to_geocentric(90.-c1);

	c0=*c0_f;
	l0=*l0_f;
	c1=*c1_f;
	l1=*l1_f;
	colat0=c0*RAD_PER_DEG;
	colat1=c1*RAD_PER_DEG;
	lon0=l0*RAD_PER_DEG;
	lon1=l1*RAD_PER_DEG;


	/*calculate distance*/
	cosdel=cos(colat0)*cos(colat1)+
		sin(colat0)*sin(colat1)*cos(lon1-lon0);
	if(cosdel >1.0) cosdel -=eps;
	if(cosdel <-1.0) cosdel +=eps;
	del=acos(cosdel);
	*dist =del/RAD_PER_DEG;

	/*calculate azimuth*/
	tmp=sin(colat0)*sin(del);
	if (tmp <=eps){ 
		/*special case: source at pole or del=0 or 180*/
		if(c0 <=eps) *az=180.;
		if(c0 >=(180.-eps)) *az=0.;
		if((*dist) <=eps) *az=-999;
		if((*dist) >=(180.-eps)) *az=-999;
		return;
	}
	cosaz=(cos(colat1)-cos(colat0)*cos(del))/tmp;
	sinaz=sin(colat1)*sin(lon1-lon0)/sin(del);
	azimuth=atan2(sinaz,cosaz)/RAD_PER_DEG;
	if(azimuth <0.) azimuth +=360.;
	*az =azimuth;
}

/* find a point dis1 degree away from the first point on the 
   great cirlce of two points.  Use OA=a1*OA1+a2*OA2 
   x1,y1,z1; x2,y2,z2:  cartesian coordinates of the points. 
			They should be in geocentric coordinates.
   dis:	  distance of the two point in degree, typically calculated by distaz above
   dis1:  distance from the first point. positive means towards the second point
   lat,lon:  latitude and longitude of the point dis1 degree away from the first point.
             it is in geocentric coordinates.		
*/
void gcpt_(double *x1_f,double *y1_f,double *z1_f,double *x2_f,double *y2_f,double *z2_f,double *dis_f,double *dis1_f,double *lat,double *lon)
{
	double sin2,a1,a2;
	double x,y,z;
	double eps=1.0e-7;
	double x1,y1,z1,x2,y2,z2,dis,dis1;

	x1=*x1_f;
	y1=*y1_f;
	z1=*z1_f;
	x2=*x2_f;
	y2=*y2_f;
	z2=*z2_f;
	dis=*dis_f;
	dis1=*dis1_f;
	dis *=RAD_PER_DEG;
//	dis1 *=RAD_PER_DEG;
	sin2=sin(dis);
	sin2 *=sin2;
	a1=(cos(dis1)-cos(dis-dis1)*cos(dis))/sin2;
	a2=(cos(dis-dis1)-cos(dis)*cos(dis1))/sin2;
	x=a1*x1+a2*x2;
	y=a1*y1+a2*y2;
	z=a1*z1+a2*z2;
	z=z/6371.0;
/*
	if(z>=1.0) z-=eps; if(z<=-1.0) z+=eps;
*/
	if(z>=1.0) z=1.0-eps; if(z<=-1.0) z=-1.0+eps;
	*lat=90.-acos(z)/RAD_PER_DEG;
//	*lat=geocentric_to_geographic((*lat));
	*lon=atan2(y,x)/RAD_PER_DEG;
}

