real*8 function earthr(xlat)
    !  this routine establishes the short distance conversion factors
    !  given the origin of coordinates
    !  the rotation angle is converted to radians also
    !  common block variables:
  
    !  local variables:
    double precision dlt1,dxlt,drad,drlt,xlat
    data re/6378.163/, ell/298.26/
  
    drad=1.7453292d-2
    drlt=9.9330647d-1
    dxlt=dble(xlat*60.0)
  
    !  conversion factor for latitude
    dlt1=datan(drlt*dtan(dxlt*drad/60.d0))
    earthr=re*(1.0-sngl(dsin(dlt1)**2)/ell)
end function earthr


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

function velocity(r,pa,ra)
    implicit real*8 (a-h,o-z)    
    include "para.inp"
    real*8 lat,lon,dep,shiftlo   
    real*8 r,pa,ra,r2d           
    real*8 v,velocity            
    common /coord/ shiftlo       
    r2d = 90./asin(1.)           
    lat=geoc_to_geog(90.0-pa*r2d)
    lon=ra*r2d+shiftlo           
    dep=ro-r                     
    call vel_3d(lon,lat,dep,v)     
    velocity=v                   
    return                       
end function velocity


real*8 function MakeMStime(iyr,month,iday,ihr,min,second)
    integer iyr,month,iday,ihr,min    !Note that year must include century 
    real*8 second
    integer basedate    !Julian day for 01/01/1970
    integer SperDay          !seconds in a 24 hour day
    integer julianDay,DeltaDays 
    data basedate/2440588/
    data SperDay/86400/

    julianDay=JULDAY(month,iday,iyr)    ! get astronomical Julian day
    DeltaDays=julianDay-basedate
    MakeMStime=DeltaDays*SperDay+ihr*3600+min*60+second
    return
end function MakeMStime


FUNCTION JULDAY(MM,ID,IYYY)                                  
    PARAMETER (IGREG=15+31*(10+12*1582))                         
    IF (IYYY.EQ.0) STOP 'There is no Year Zero.'                
    IF (IYYY.LT.0) IYYY=IYYY+1                                   
    IF (MM.GT.2) THEN                                            
       JY=IYYY                                                    
       JM=MM+1                                                    
    ELSE                                                         
       JY=IYYY-1                                                  
       JM=MM+13                                                   
    ENDIF                                                        
    JULDAY=INT(365.25*JY)+INT(30.6001*JM)+ID+1720995             
    IF (ID+31*(MM+12*IYYY).GE.IGREG) THEN                        
       JA=INT(0.01*JY)                                            
       JULDAY=JULDAY+2-JA+INT(0.25*JA)                            
    ENDIF                                                        
    RETURN                                                       
END function JULDAY


real*8 function weighting(jj,dist,res)
    integer j,jj
    real*8 dist,res,dpi
    data xnear,xfar/80.0,800.0/                                        
    data tres/1.0/     
    dpi=asin(1.0)/90.0                                                
    j=mod(jj,5)                                                        
    wt=1.0
    !-picking factor                                                           
    if (j.ge.0.and.j.le.4) then                                          
!       wt=wt*float(4-j)/4.                                              
    else                                                               
!       wt=0.0                                                           
    endif
  
    !-distance factor
!    if (dist.gt.xnear) then                                              
!       wt=wt*(xfar-xnear)/(9.*dist+xfar-10.*xnear)                      
!    endif 
!    wt=wt*(dist/(180.*dpi*ro))
  
    !-residual factor                                                             
!    wt=wt*(tres/(tres+abs(res)))**2                                    
    weighting=wt                                                       
end function weighting                                               


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


real*8 function angle(w,pt1,pt2)
    implicit integer (i-n)
    include "para.inp"
    integer pt1,pt2
    real*8 w(3,maxtrpts),dpi
    real*8 a,b,c,dx,dy
    dpi=asin(1.0)/90.0

    call cal_delta(w(2,pt1),w(3,pt1)-0.5,w(2,pt1),w(3,pt1)+0.5,dx) ! calculate the local degree-2-km unit of event in x direction
    call cal_delta(w(2,pt1)-0.5,w(3,pt1),w(2,pt1)+0.5,w(3,pt1),dy) ! calculate the local degree-2-km unit of event in y direction
    a=sqrt(((w(3,pt1)-w(3,pt2))*dx)**2+((w(2,pt1)-w(2,pt2))*dy)**2)
    b=w(1,pt2)-w(1,pt1)
    c=sqrt(a**2+b**2)
    angle=acos(b/c)/dpi

    return
end function
  

                                                                   
