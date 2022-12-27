subroutine pn_path(evt_lon,evt_lat,evt_dep,sta_lon,sta_lat,sta_elv,w,npt,ttime,ie,is)
	implicit real*8 (a-h,o-z)
	include "para.inp"
	real*8 length
	real*8 lat_1,lon_1,lat_2,lon_2,lat,lon
        real*8 theta1,theta2,vc,vm
        real*8 ptlat(2000),ptlon(2000)
        real*8 w(3,maxtrpts+1)
	parameter(len_seg=2.0)
        parameter(r0=6371.0)
        integer ie,is !point to correct raytracing connect line
        real*8 geoc_to_geog, geog_to_geoc
	common/npoint/npt_e,npt_s,npt_pn
	dpi=asin(1.0)/90.0
        vc = 0.0
        vm = 0.0
	do k=2,ndep_c-1
	   do j=2,nlat_c-1
	      do i=2,nlon_c-1
		 vc=vc+vp(i,j,k)
	      enddo
	   enddo
	enddo
	vc=vc/((ndep_c-2)*(nlat_c-2)*(nlon_c-2))
	do j=2,nlat_m-1
	   do i=2,nlon_m-1
	      vm=vm+vpn(i,j)
	   enddo
	enddo
	vm=vm/((nlat_m-2)*(nlon_m-2))
!        print*,"vc:",vc,"vm:",vm
!        pause
! convert geographic latitude to geocentric latitude
!	print*,evt_lon,evt_lat,sta_lon,sta_lat
!	call depth(evt_lon,evt_lat,dm_evt)
!	call depth(sta_lon,sta_lat,dm_sta)
	dm_evt=dm
	dm_sta=dm
!        theta1=geog_to_geoc(evt_lat)
!        theta2=geog_to_geoc(sta_lat)
        theta1=evt_lat
        theta2=sta_lat
!	print*,"theta1,theta2",theta1,theta2
        sina=vc/vm
        alfa=asin(sina)
        rm_evt=r0-dm_evt
        raypar=rm_evt/vm
        re=r0-evt_dep
        sina1=raypar*vc/re
        alfa1=asin(sina1)
        dis1=alfa-alfa1
	rs=r0-sta_elv
	rm_sta=r0-dm_sta
	raypar=rm_sta/vm
        sina2=raypar*vc/rs
        alfa2=asin(sina2)
        dis2=alfa-alfa2
!       theta1=lat_1
!       theta2=lat_2
!       print*,"lat_1=",theta1
!       print*,"lat_2=",theta2
        theta1=90-theta1
        theta2=90-theta2
        phi1=evt_lon
        phi2=sta_lon
! calculate the distance and aizmuth
	call distaz(theta1,phi1,theta2,phi2,dis,az)
!       print*,"dis=",dis
!       print*,"az=",az
! convert to cartesian coordinates
        theta1=theta1*dpi
        theta2=theta2*dpi
        phi1=phi1*dpi
        phi2=phi2*dpi
        x1=r0*sin(theta1)*cos(phi1)
        y1=r0*sin(theta1)*sin(phi1)
        z1=r0*cos(theta1)
        x2=r0*sin(theta2)*cos(phi2)
        y2=r0*sin(theta2)*sin(phi2)
        z2=r0*cos(theta2)
        call gcpt(x1,y1,z1,x2,y2,z2,dis,dis1,lat_1,lon_1)
!       lat_1=geoc_to_geog(lat_1)
        call gcpt(x2,y2,z2,x1,y1,z1,dis,dis2,lat_2,lon_2)
!       lat_2=geoc_to_geog(lat_2)
!        print*,"lat_1,lon_1",lat_1,lon_1
!  	print*,"lat_2,lon_2",lat_2,lon_2
!       print*,"lon=",lon
! calculate the ray points in the mantle
        call surfpath(lat_1,lat_2,lon_1,lon_2,ddel,nseg,ptlat,ptlon,distance)
        npt_pn=nseg+1
!       print*,"npt_pn,distance",npt_pn,distance
!       open(1,file="raypath",status="unknown")
! calculate travel time
!	call depth(lon_1,geoc_to_geog(lat_1),dm)
	rm=r0-dm
        colat_1=(90-lat_1)*dpi
        lon_1=lon_1*dpi
        x1=re*sin(theta1)*cos(phi1)
        y1=re*sin(theta1)*sin(phi1)
        z1=re*cos(theta1)
        x0=rm*sin(colat_1)*cos(lon_1)
        y0=rm*sin(colat_1)*sin(lon_1)
        z0=rm*cos(colat_1)
        length=sqrt((x1-x0)**2+(y1-y0)**2+(z1-z0)**2)
	ni=int(length/len_seg)
	if (ni.lt.2) ni=2
	dx=(x0-x1)/ni
        dy=(y0-y1)/ni
        dz=(z0-z1)/ni
	ttime=0.0
	npt=0
!	print*,"evt_lon,evt_lat,evt_dep",evt_lon,evt_lat,evt_dep
!	pause
	call vel_3d(evt_lon,evt_lat,evt_dep,v1)
	w(1,1)=evt_dep
	w(2,1)=evt_lat
	w(3,1)=evt_lon
	do i=2,ni+1
	   x=x1+dx
	   y=y1+dy
	   z=z1+dz
	   r=sqrt(x*x+y*y+z*z)
	   acosa=z/r
	   if (acosa.lt.-1) acosa=-1
	   if (acosa.gt.1) acosa=1
	   a=acos(acosa)
	   acosa=x/r/sin(a)
	   if (acosa.lt.-1) acosa=-1
           if (acosa.gt.1) acosa=1
	   b=acos(acosa)
           if (y.lt.0.00) b=360*dpi-b
	   ca=90.0-a/dpi
	   b=b/dpi
	   d=r0-r
	   w(1,i)=d
	   w(2,i)=geoc_to_geog(ca)
	   w(3,i)=b
	   call vel_3d(w(3,i),w(2,i),w(1,i),v2)
	   dl=sqrt((x-x1)**2+(y-y1)**2+(z-z1)**2)
	   ttime=ttime+dl*((1/v1+1/v2)/2)
	   x1=x
	   y1=y
	   z1=z
	   v1=v2
	enddo
	npt_e=ni+1
	npt=npt_e
!	print*,"ttime",ttime
!	print*,"npt_e",npt_e
!	pause

!	call depth(lon_2,geoc_to_geog(lat_2),dm)
	rm=r0-dm
        colat_2=(90-lat_2)*dpi
        lon_2=lon_2*dpi
        x0=rm*sin(colat_2)*cos(lon_2)
        y0=rm*sin(colat_2)*sin(lon_2)
        z0=rm*cos(colat_2)
!	print*,x0,y0,z0
	x2=rs*sin(theta2)*cos(phi2)
        y2=rs*sin(theta2)*sin(phi2)
        z2=rs*cos(theta2)
!	print*,x2,y2,z2
        length=sqrt((x2-x0)**2+(y2-y0)**2+(z2-z0)**2)
!	print*,length
	ni=int(length/len_seg)
	if (ni.lt.2) ni=2
	dx=(x0-x2)/ni
	dy=(y0-y2)/ni
	dz=(z0-z2)/ni
	call vel_3d(sta_lon,sta_lat,sta_elv,v1)
!	print*,v1
        w(1,1+npt)=sta_elv
        w(2,1+npt)=sta_lat
        w(3,1+npt)=sta_lon
	do i=2,ni+1
	   x=x2+dx
           y=y2+dy
           z=z2+dz
           r=sqrt(x*x+y*y+z*z)
           acosa=z/r
           if (acosa.lt.-1) acosa=-1
           if (acosa.gt.1) acosa=1
           a=acos(acosa)
           acosa=x/r/sin(a)
           if (acosa.lt.-1) acosa=-1
           if (acosa.gt.1) acosa=1
           b=acos(acosa)
           if (y.lt.0.00) b=360*dpi-b
           ca=90.0-a/dpi
           b=b/dpi
           d=r0-r
           w(1,i+npt)=d
           w(2,i+npt)=geoc_to_geog(ca)
           w(3,i+npt)=b
           call vel_3d(w(3,i+npt),w(2,i+npt),w(1,i+npt),v2)
!	   print*,v2
           dl=sqrt((x-x2)**2+(y-y2)**2+(z-z2)**2)
!	   print*,dl
!	   pause
           ttime=ttime+dl*((1/v1+1/v2)/2)
           x2=x
           y2=y
	   z2=z
           v1=v2
        enddo
        npt_s=ni+1
        npt=npt+npt_s
!	print*,"ttime",ttime
!       print*,"npt_s",npt_s
!       pause
! for pn part
!	call depth(ptlon(1),geoc_to_geog(ptlat(1)),dm)
	rm=r0-dm
        x1=rm*sin((90.0-ptlat(1))*dpi)*cos(ptlon(1)*dpi)
        y1=rm*sin((90.0-ptlat(1))*dpi)*sin(ptlon(1)*dpi)
        z1=rm*cos((90.0-ptlat(1))*dpi)
	w(1,1+npt)=dm
        w(2,1+npt)=geoc_to_geog(ptlat(1))
        w(3,1+npt)=ptlon(1)
        call vel_pn(w(3,1+npt),w(2,1+npt),v1)
        do ipt=2,npt_pn
!          ptlat(ipt)=geoc_to_geog(ptlat(ipt))
!          write(1,*)ptlon(ipt),ptlat(ipt)
!	   call depth(ptlon(ipt),geoc_to_geog(ptlat(ipt)),dm)
	   rm=r0-dm
           x2=rm*sin((90.0-ptlat(ipt))*dpi)*cos(ptlon(ipt)*dpi)
           y2=rm*sin((90.0-ptlat(ipt))*dpi)*sin(ptlon(ipt)*dpi)
           z2=rm*cos((90.0-ptlat(ipt))*dpi)
           length=sqrt((x1-x2)**2+(y1-y2)**2+(z1-z2)**2)
	   w(1,npt+ipt)=dm
           w(2,npt+ipt)=geoc_to_geog(ptlat(ipt))
           w(3,npt+ipt)=ptlon(ipt)
           call vel_pn(w(3,ipt+npt),w(2,ipt+npt),v2)
!          print*,"length=",length
!          print*,"v1=",v1
!          print*,"v2=",v2
!	   pause
           ttime=ttime+length*((1/v1+1/v2)/2)
	   x1=x2
           y1=y2
           z1=z2
           v1=v2
        enddo
        npt=npt+npt_pn
        w(2,1)=geoc_to_geog(w(2,1))
        w(2,npt_e+1)=geoc_to_geog(w(2,npt_e+1))
        ie=npt_e
        is=npt_s
!	print*,"npt,ttime",npt,ttime
	end
