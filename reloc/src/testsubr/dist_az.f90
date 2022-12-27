      program main
        real*8 dist,az
        real*8 elat,elong,edep,slat,slong,selv

        elat=21.572
        elong=105.763
        edep=14.100
        slat=21.465
        slong=105.646
        selv=-1.200

        call dist_az(elat,elong,edep,slat,slong,selv,dist,az)
        print*,"dist",dist
        print*,"az",az
      end program
!   Add subroutine dist_az to replace the distaz on wrote on C
      subroutine dist_az(ey,ex,ez,sy,sx,sz,dis,azi)
!    Subroutine calculate distace from 2 point using Euclidean formula
!        1- Calculate Earth radius at each point (used average value)
!        2- Calculate geocentric latitude at each point
!        3- Convert spherical coordinates of each point to cartesian coordinates, 
!        that is, from calculated radius, geocentric latitude and longitude to x,y,z.
!        4- Calculate distance using Euclidean distance formula
!     References: https://planetcalc.com/7725/
        implicit real(a-h,o-z)
         parameter(ro=6371.000)
         parameter(pi=3.14159265359)
         parameter(a=6378.1370)!semi-major axis of the ellipsoid along lattitude
         parameter(b=6356.7523)!semi-major axis of the ellipsoid along longitude
         real*8 ey,ex,ez,sy,sx,sz
         real*8 coelat,coslat,dist,az,coelat1,coslat1
         real*8 radperdeg,xe,ye,ze,xs,ys,zs
         real*8 geog_to_geoc, geoc_to_geog
         real*8 dis, azi
        
!     Calculate the mediate conversion rate
         radperdeg = pi/180.
!         print*,radperdeg
!     Convert latitude of geographic to geocentric
!         print*,"#",ey,ex,ez,sy,sx,sz
         coelat=atan((b/a)*tan(ey*radperdeg))/radperdeg
         coslat=atan((b/a)*tan(sy*radperdeg))/radperdeg
!         coelat1=geog_to_geoc(ey)
!         coslat1=geog_to_geoc(sy)
!         print*,coelat,coslat,coelat1,coslat1
!     Now calculate x,y,z point for station and event
         xe=ro*cos(coelat*radperdeg)*cos(ex*radperdeg)
         ye=ro*cos(coelat*radperdeg)*sin(ex*radperdeg)
         ze=ro*sin(coelat*radperdeg)
        
         xs=ro*cos(coslat*radperdeg)*cos(sx*radperdeg)
         ys=ro*cos(coslat*radperdeg)*sin(sx*radperdeg)
         zs=ro*sin(coslat*radperdeg)         
!     Now calculate the Euclidean distance
         dis=sqrt(((xs-xe)*(xs-xe))+((ys-ye)*(ys-ye))+((zs-ze)*(zs-ze)))
!     Calculate the azimuth
         azi=azimuth(coelat,ex,coslat,sx) 


      end subroutine

      real*8 function azimuth(x1,y1,x2,y2)
       implicit none
       real*8 x1,y1,x2,y2,a,b,dpi
       dpi=asin(1.0)/90.0
       a=x2-x1
       b=y2-y1
       if (a>0) then
         azimuth=90.-(atan(b/a)/dpi)
       else
         azimuth=(atan(b/(-1*a))/dpi)+270.
       endif
       return
      end function

      real*8 function geog_to_geoc(xla)
       implicit none
       real*8 xla,RAD_PER_DEG,B2A_SQ
       RAD_PER_DEG=0.0174532925199432955
       B2A_SQ=0.993305521
       geog_to_geoc = atan(B2A_SQ*tan(RAD_PER_DEG*xla)) / RAD_PER_DEG
       return
      end function geog_to_geoc

      real*8 function geoc_to_geog(xla)
       implicit none
       real*8 xla,RAD_PER_DEG,B2A_SQ
       RAD_PER_DEG=0.0174532925199432955
       B2A_SQ=0.993305521
       geoc_to_geog = atan(tan(RAD_PER_DEG*xla)/B2A_SQ) / RAD_PER_DEG
       return
      end function geoc_to_geog