/*
 *surfpath.c: find surface projection of raypath between sta and evt
 */

#include <stdio.h>
#include <stdlib.h>
#include <math.h>

void surfpath_(double *lat1_f,double *lat2_f,double *lon1_f,double *lon2_f,double *ddel_f,int *npt,double *ptlat,double *ptlon,double *distance)
{
	double the1,the2,lat,phe1,phe2,lon,dis1,dis;
	double x1,x2,x,y1,y2,y,z1,z2,z,px,py,pz,p,plat,plon;
	double a1,a2,a3,b1,b2,b3;
	double pi=3.14165926; double zero=1.0E-10;
	double lat1,lat2,lon1,lon2,ddel;
	int i,n;

lat1=*lat1_f;
lat2=*lat2_f;
lon1=*lon1_f;
lon2=*lon2_f;
ddel=*ddel_f;
the1=(90.0-lat1)*pi/180.0;
phe1=lon1*pi/180.0;
the2=(90.0-lat2)*pi/180.0;
phe2=lon2*pi/180.0;

x1=sin(the1)*cos(phe1);
y1=sin(the1)*sin(phe1);
z1=cos(the1);
x2=sin(the2)*cos(phe2);
y2=sin(the2)*sin(phe2);
z2=cos(the2);

 /* dist between the two points */
dis=acos(x1*x2+y1*y2+z1*z2)*180.0/pi; 
//fprintf(stderr,"dist between this two points: %f\n",dis);
*distance=dis;

 /* the pole of the great circle */
px=y1*z2-y2*z1;
py=x2*z1-x1*z2;
pz=x1*y2-x2*y1;
p=sqrt(px*px+py*py+pz*pz);

if (p <= zero){
     fprintf(stderr,"Two same points. Program exits!\n");
     exit(-1);
   }

px=px/p;
py=py/p;
pz=pz/p;
plat=90.0-acos(pz)*180.0/pi;
plon=atan2(py,px)*180.0/pi;

n=(int)(dis/ddel);

ptlon[0]=lon1;
ptlat[0]=lat1;
if (ptlon[0]<0)
	ptlon[0]=ptlon[0]+360.0;

for(i=1;i<=n;i++){
   dis1=(double)i * (double)ddel;
   dis *=pi/180.0;
   dis1*=pi/180.0;

/* point dis1 degree away from the first point.
    use cross product relation of vectors. */ 
   a1=sin(dis1)/sin(dis-dis1);
   a2=(y1+y2*a1)/(x1+x2*a1);
   a3=(z1+z2*a1)/(x1+x2*a1);
   b1=sin(dis1)/sin(dis);
   b2=a2*x1-y1;

   if(b2<=zero){
      b3=a3*x1-z1;
      x=b1*(x1*z2-x2*z1)/b3;
      }
   else{
     x=b1*(x1*y2-x2*y1)/b2;
     }
   y=a2*x;
   z=a3*x;
   lat=90.0-acos(z)*180.0/pi;
   lon=atan2(y,x)*180.0/pi;
   dis1 *=180.0/pi;
   dis *=180.0/pi;
   ptlon[i]=lon;
   ptlat[i]=lat;
   if (ptlon[i]<0.0)
	ptlon[i]=ptlon[i]+360;
} /*for*/
   *npt=n+1;
   ptlon[*npt]=lon2;
   ptlat[*npt]=lat2;
   if (ptlon[*npt]<0)
	ptlon[*npt]=ptlon[*npt]+360;
}
