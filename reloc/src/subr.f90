!-Read absolute travel time of Haijiang's format. Only read the P wave travel time.
subroutine getdata1(evt_fn,sta_fn,abs_fn)
   implicit integer (i-n)
   include "para.inp"
   character*20 evt_fn,sta_fn,abs_fn
   character*90 line
   character*1 str
   character*2 phcard
   character*5 stnm
   real*8 junk1,junk2,junk3
   real*8 ttime,wei
   real*8 evt_lat,evt_lon,evt_dep,sta_lat,sta_lon,sta_elv
   real*8 evtla,stala,coevtla,costala,dist,az
   real*8 dist1(maxrds),az1(maxrds)
   real*8 geog_to_geoc,geoc_to_geog

   open(1,file=evt_fn,status="old")
   ievt=1
   do
      read(1,*,iostat=iret)evt_date(ievt),evt_time(ievt),(evt_loc(ievt,j),j=1,3),evt_mag(ievt),junk1,junk2,junk3,evt_id(ievt)
      if (iret<0) exit
      !-giving a fake time code for relocation
      if (evt_date(ievt)==0) evt_date(ievt)=20001201
      if (evt_time(ievt)==0) evt_time(ievt)=14302005
!       print*,evt_date(ievt),evt_time(ievt),(evt_loc(ievt,j),j=1,3),evt_mag(ievt),junk1,junk2,junk3,evt_id(ievt)
      ievt=ievt+1
   enddo
   close(1)
   nevt=ievt-1
   if (nevt>maxevt) stop "event array is not enough!"  ! check event dimension 

   open(2,file=sta_fn,status="old")
   ista=1
   do
      read(2,*,iostat=iret)sta_nm(ista),(sta_loc(ista,j),j=1,3),(sta_cor(ista,j),j=1,2)
      if (iret<0) exit
!      print*,sta_nm(ista),(sta_loc(ista,j),j=1,3),(sta_cor(ista,j),j=1,2)
      ista=ista+1
   enddo
   close(2)
   sta_loc(:,3)=sta_loc(:,3)/(-1000)  ! compatible with the earthquake depth in km unit (downward "+")
   nsta=ista-1
   if (nsta>maxsta) stop "station array is not enough!"  ! check station dimension

   open(3,file=abs_fn,status="old")
   ipair=0
   do
      read(3,"(a)",iostat=iret)line
      if (iret<0) exit
      if (line(1:1)=="#") then
         read(line,*)str,idevt
         !-find event index
         idxevt=0
         do i=1,nevt
            if (idevt==evt_id(i)) then
               idxevt=i
               exit
            endif
         enddo
         cycle
      else
         read(line,*)stnm,ttime,wei,phcard
         !-find station index
         idxsta=0
         do i=1,nsta
            if (stnm==sta_nm(i)) then
               idxsta=i
               exit
            endif
         enddo
      endif
      !-then pair the evt, sta, and tt data
      if ((idxevt.ne.0).and.(idxsta.ne.0)) then
         ipair=ipair+1
         evt_pack(idxevt+1)=ipair
         evt_idx(ipair)=idxevt
         sta_idx(ipair)=idxsta
         abs_tt(ipair)=ttime
         ray_wei(ipair)=int(wei)
         ! P-phase group
         if (phcard=="Pg") ray_pha(ipair)=1
         if (phcard=="Pn") ray_pha(ipair)=3
         ! S-phase group
         if (phcard=="Sg") ray_pha(ipair)=2
         if (phcard=="Sn") ray_pha(ipair)=4
!         print*,evt_id(idxevt),sta_nm(idxsta),ttime,wei,phcard,ray_pha(ipair)
      endif
   enddo
   close(3)
   nray=ipair
   if (nray>maxrds) stop "readings array is not enough!"  ! check reading dimension
!
   do i=1,nray
         evt_lat=evt_loc(evt_idx(i),1)
         evt_lon=evt_loc(evt_idx(i),2)
         evt_dep=evt_loc(evt_idx(i),3)
         sta_lat=sta_loc(sta_idx(i),1)
         sta_lon=sta_loc(sta_idx(i),2)
         sta_elv=-sta_loc(sta_idx(i),3)/1000
!		print*,evtlat,"#",stalat
!
!	  print*,evtlat,evtlon,evtdep,"#",stalat,stalon,staelv
!	  pause
!
         evtla=geog_to_geoc(evt_lat)
         stala=geog_to_geoc(sta_lat)
!		print*,evtla,"#",stala
!
         coevla=90.-evtlat
         costla=90.-stalat
!		print*,coevla1,"#",costla1
!
         call distaz(90.0-evtla,evt_lon,90.0-stala,sta_lon,dist,az)
!		print*,coevla1,evtlon,costla1,stalon,dist1,az1
!		print*,dist1,"#",az1   
!	    pause
        dist1(i)=dist
        az1(i)=az

   enddo

   open(4,file="init_data",status="unknown")
   do i=1,nray
      write(4,10)sta_nm(sta_idx(i)),evt_id(evt_idx(i)),evt_date(evt_idx(i)),evt_time(evt_idx(i)),&
                 (evt_loc(evt_idx(i),j),j=1,3),evt_mag(evt_idx(i)),(sta_loc(sta_idx(i),j),j=1,3),&
                 abs_tt(i),dist1(i)*111,az1(i),i,ray_pha(i)
   enddo
   close(4)
10  format(a,1x,3(i9,1x),10(f8.3,1x),2i6)
endsubroutine


subroutine getdata2(evt_fn,sta_fn)
    implicit integer (i-n)
    include "para.inp"
    character*90 line
    character*5 stnm
    character*20 evt_fn,sta_fn
    integer yy,mm,dd,hr,min,lat,lon,min0
    real*8 sec,xlat,xlon,dep,mag,psec0,ssec0,ptt,stt

    open(1,file=sta_fn,status="old")
    ista=1
    do
       read(1,*,iostat=iret)sta_nm(ista),(sta_loc(ista,j),j=1,3),(sta_cor(ista,j),j=1,2)
       if (iret<0) exit
       ista=ista+1
    enddo
    close(1)
    sta_loc(:,3)=sta_loc(:,3)/(-1000)  ! compatible with the earthquake depth in km unit (downward "+")
    nsta=ista-1
    if (nsta>maxsta) stop "station array is not enough!"  ! check station dimension

    open(2,file=evt_fn,status="old")
    ievt=0
    idxevt=0
    ipair=0
    do
       read(2,"(a)",iostat=iret)line
       if (iret<0) exit
       if (line(2:2)>="0".and.line(2:2)<="9") then
          read(line,10)yy,mm,dd,hr,min,sec,lat,xlat,lon,xlon,dep,mag
          ievt=ievt+1
          evt_id(ievt)=ievt
          idxevt=ievt
          evt_date(ievt)=yy*10000+mm*100+dd
          evt_time(ievt)=hr*1000000+min*10000+int(sec*100)
          evt_loc(ievt,1)=real(lat)+xlat/60.
          evt_loc(ievt,2)=real(lon)+xlon/60.
          evt_loc(ievt,3)=dep
          evt_mag(ievt)=mag
          cycle
       else
          read(line,20)stnm,min0,psec0,ipwei,ssec0,iswei
          !-find station index
          ptt=(real(min0)*60.+psec0)-(real(min)*60.+sec)
          stt=(real(min0)*60.+ssec0)-(real(min)*60.+sec)
          idxsta=0
          do i=1,nsta
             if (stnm==sta_nm(i)) then
                idxsta=i
                exit
             endif
          enddo
       endif
       !-then pair the evt, sta, and tt data
       if (ipwei>4.or.iswei>4) cycle   ! S-P records have still not been capable in this program
       if (idxevt.ne.0.and.idxsta.ne.0) then
          if (psec0.ne.0.0) then
             ipair=ipair+1
             evt_pack(idxevt+1)=ipair
             evt_idx(ipair)=idxevt
             sta_idx(ipair)=idxsta
             abs_tt(ipair)=ptt
             ray_wei(ipair)=ipwei
             ray_pha(ipair)=1
          endif
          if (ssec0.ne.0.0) then
             ipair=ipair+1
             evt_pack(idxevt+1)=ipair
             evt_idx(ipair)=idxevt
             sta_idx(ipair)=idxsta
             abs_tt(ipair)=stt
             ray_wei(ipair)=iswei
             ray_pha(ipair)=2             
          endif
       endif
    enddo
    close(2)
    nevt=ievt
    nray=ipair
    if (nevt>maxevt) stop "event array is not enough!"  ! check event dimension 

10  format(1x,i4,4i2,f6.2,i2,f5.2,i3,f5.2,f6.2,f4.2)
20  format(1x,a5,15x,i2,f6.2,6x,i1,3x,f6.2,6x,i1)
end subroutine


!-Read 3d velocity-model
subroutine loadvel(vel_fn,nmod)
    implicit integer (i-n)
    include "para.inp"
    character*20 vel_fn
    real*8 lat_ave,earthr
    integer nmod

    print*,"nmod = ",nmod
    open(1,file=vel_fn)
    read(1,*)bld1,nlon_m,nlat_m
!    read(1,*)bld2,nlon_moho,nlat_moho
    read(1,*)bld3,bld4,nlon_c,nlat_c,ndep_c
    read(1,*)(lon_m(i),i=1,nlon_m)
    read(1,*)(lat_m(i),i=1,nlat_m)
!    read(1,*)(lon_moho(i),i=1,nlon_moho)
!    read(1,*)(lat_moho(i),i=1,nlat_moho)
    read(1,*)(lon_c(i),i=1,nlon_c)
    read(1,*)(lat_c(i),i=1,nlat_c)
    read(1,*)(dep_c(i),i=1,ndep_c)

    !-read Pn velocity model
    do j=1,nlat_m
       read(1,*)(vpn(i,j),i=1,nlon_m)
    enddo
    !-read moho depth
!    do j=1,nlat_moho
!       read(1,*)(mhd(i,j),i=1,nlon_moho)
!    enddo
    !-read P velocity model
    do k=1,ndep_c
       do j=1,nlat_c
          read(1,*)(vp(i,j,k),i=1,nlon_c)
       enddo
    enddo
    if(nmod.eq.1) then
      goto 222
    elseif(nmod.eq.2) then
      goto 111
    else
      write(*,'(a)')"Error! nmod input wrong! pls chekc input file"
    endif
 111 continue 
    !-read Sn velocity model
!   do j=1,nlat_m
!     read(1,*)(vsn(i,j),i=1,nlon_m)
!  enddo
    !-read S velocity model
!  do k=1,ndep_c
!     do j=1,nlat_c
!        read(1,*)(vs(i,j,k),i=1,nlon_c)
!     enddo
!  enddo
 222 close(1)

    call bldmap

    !-obtain the average earth radius to velocity model
    lat_ave=0.0
    do i=1,nlat_c
       lat_ave=lat_ave+lat_c(i)
    enddo
    lat_ave=lat_ave/real(nlat_c)
!   ro=earthr(lat_ave)
    ro=6371.0
end subroutine


!-correct original time code
subroutine ot_corr(dt,date,time)
    implicit integer (i-n)
    include "para.inp"
    integer date,time,yy,mm,dd,hr,min
    real*8 MakeMStime,sec,secp,dt,ot,ct
 
    !-decode original date & time
    yy=int(evt_date(evtnow)/10000)
    mm=int((evt_date(evtnow)-yy*10000)/100)
    dd=evt_date(evtnow)-yy*10000-mm*100
    hr=int(evt_time(evtnow)/1000000)
    min=int((evt_time(evtnow)-hr*1000000)/10000)
    sec=real(evt_time(evtnow)-hr*1000000-min*10000)/100.

    ot=MakeMStime(yy,mm,dd,hr,min,sec)
    ot=ot+dt

    imatch=0
    do i=0,23
       do j=0,59
          do k=0,59
             ct=MakeMStime(yy,mm,dd,i,j,dble(real(k)))
             if (int(ct)==int(ot)) then
                imatch=1
                secp=ot-ct
                hr=i
                min=j
                sec=k+secp
                !-encode original time
                time=int(hr*1000000+min*10000+sec*100)
                goto 99
             endif
          enddo
       enddo
    enddo
99  continue
    if (imatch==1) then
       date=evt_date(evtnow)
    else
       print*,"original time correction error!!!"
    endif
end subroutine


!-This subroutine is to change Lon. Lat. to Km unit
subroutine cal_delta(elat,elon,slat,slon,delta)
     real*8 elat,elon,slat,slon,delta
     real*8 avlat,a,b,dlat,dlon,dx,dy
    
     avlat=0.5*(elat+slat)
     a=1.840708+avlat*(.0015269+avlat*(-.00034+avlat*(1.02337e-6)))
     b=1.843404+avlat*(-6.93799e-5+avlat*(8.79993e-6+avlat*(-6.47527e-8)))
     dlat=slat-elat
     dlon=slon-elon
     dx=a*dlon*60.
     dy=b*dlat*60.
     delta=sqrt(dx*dx+dy*dy)
endsubroutine


!=====
! find the gap angle from azimuth array
!=====
subroutine get_gap(a,n,gap)
   integer n
   integer a(n),gap,tmp(n)

   if (n==1) then
      gap=360
      return
   endif

   call qk_sort(a,n,1,n)

   do i=1,n-1
      tmp(i)=a(i+1)-a(i)
   enddo
   tmp(n)=360-(a(n)-a(1))

   gap=maxval(tmp(1:n))
   return
endsubroutine


!=====
! Quick sort method
!=====
recursive subroutine qk_sort(a,n,s,e)
    integer n    ! the dimension of array
    integer a(n) ! the array
    integer s    ! starting position
    integer e    ! ending position
    integer l,r  ! a(l)>k, a(r)<k
    integer k    ! key value
    integer temp ! temporary memory
    
    l=s
    r=e+1
    if(r<=l) return
     
    k=a(s)
    do while (.true.)
       do while (.true.)
          l=l+1
          if ((a(l)>k).or.(l>=e)) exit
       enddo
       
       do while (.true.)
          r=r-1
          if ((a(r)<k).or.(r<=s)) exit
       enddo
 
       if (r<=l) exit
       
       temp=a(l)
       a(l)=a(r)
       a(r)=temp
    enddo
     
    temp=a(s)
    a(s)=a(r)
    a(r)=temp
     
    call qk_sort(a,n,s,r-1)
    call qk_sort(a,n,r+1,e)    
     
    return
end subroutine


!=====
! researching the minimum by lottery way
!=====
subroutine lottery(lx,ly,lz)
    include "para.inp"
    real*8 rand(3),db(2),nbfactor
    real*8 lx,ly,lz,sta_min(2),sta_max(2)
    data db /5.0,25.0/
    data nbfactor/0.8/
    
    call random_number(rand)
!   write(1,*)(rand(i),i=1,3)
    sta_min(1)=sta_loc(sta_idx(ray_idx(1)),1)
    sta_max(1)=sta_loc(sta_idx(ray_idx(1)),1)
    sta_min(2)=sta_loc(sta_idx(ray_idx(1)),2)
    sta_max(2)=sta_loc(sta_idx(ray_idx(1)),2)
    do i=2,nrow
       if (sta_loc(sta_idx(ray_idx(i)),1)>sta_max(1)) sta_max(1)=sta_loc(sta_idx(ray_idx(i)),1)
       if (sta_loc(sta_idx(ray_idx(i)),1)<sta_min(1)) sta_min(1)=sta_loc(sta_idx(ray_idx(i)),1)
       if (sta_loc(sta_idx(ray_idx(i)),2)>sta_max(2)) sta_max(2)=sta_loc(sta_idx(ray_idx(i)),2)
       if (sta_loc(sta_idx(ray_idx(i)),2)<sta_min(2)) sta_min(2)=sta_loc(sta_idx(ray_idx(i)),2)
    enddo
    !-perturb to research another minimum
    ly=ly*nbfactor+(sta_min(1)*(1.-rand(1))+sta_max(1)*rand(1))*(1.-nbfactor)
    lx=lx*nbfactor+(sta_min(2)*(1.-rand(2))+sta_max(2)*rand(2))*(1.-nbfactor)
    if (lz<db(1)) lz=db(1)
    if (lz>db(2)) lz=db(2)
    lz=lz*nbfactor+(db(1)*(1.-rand(3))+db(2)*rand(3))*(1.-nbfactor)

    return
endsubroutine


!=====
! earthquake initial guess
!=====
subroutine ini_hypo(x0,y0,z0)
    implicit integer (i-n)
    include "para.inp"
    real*8 x0,y0,z0
 
    x0=0.0
    y0=0.0
    z0=0.0
    if ((evt_loc(evtnow,1).ne.0.0).and.(evt_loc(evtnow,2).ne.0.0)) then  ! use initial located hypo
       y0=evt_loc(evtnow,1)
       x0=evt_loc(evtnow,2)
       z0=evt_loc(evtnow,3)
!      y0=22.0
!      x0=104.0
!      z0=17.0
    else  ! take average point of recorded stations
       icount=0
       do i=evt_pack(evtnow)+1,evt_pack(evtnow+1)
          icount=icount+1
          y0=y0+sta_loc(sta_idx(i),1)
          x0=x0+sta_loc(sta_idx(i),2)
       enddo 
       y0=y0/icount
       x0=x0/icount
       z0=20.0
       print*,x0,y0,z0
    endif
!   print*,"initial guess:"
!   write(*,"(3(a4,g20.10))")"x0:",x0,"y0:",y0,"z0:",z0
!   print*,"---------------------------------------"                                                 
end subroutine


subroutine bldmap
    implicit integer (i-n)
    include "para.inp"
    real*8 lon_now,lat_now,dep_now

    !-for crustal velocity 3d-model
    lon1_c=bld3-lon_c(1)
    ilonmax=(1e-10)+(lon_c(nlon_c)+lon1_c)/bld3
    lat1_c=bld3-lat_c(1)
    ilatmax=(1e-10)+(lat_c(nlat_c)+lat1_c)/bld3
    dep1_c=bld4-dep_c(1)
!    print*,"dep1_c",dep1_c
    idepmax=(1e-10)+(dep_c(ndep_c)+dep1_c)/bld4
    if ((ilonmax.gt.ilondeg).or.(ilatmax.gt.ilatdeg).or.(idepmax.gt.idepkm)) then
       print*,"Error, model dimension out of range!"
       stop
    endif
    ilon=1
    do i=1,ilonmax
       ilon1=ilon+1
!      print*,ilon1
       lon_now=float(i)*bld3-lon1_c
!      print*,lon_now,lon_c(ilon1)
       if (lon_now.ge.lon_c(ilon1)) ilon=ilon1
       ilonloc_c(i)=ilon
    enddo
    do i=ilonmax+1,ilondeg
       ilonloc_c(i)=0
    enddo
    ilat=1
    do i=1,ilatmax
       ilat1=ilat+1
!      print*,ilat1
       lat_now=float(i)*bld3-lat1_c
!      print*,lat_now,lat_c(ilat1)
       if (lat_now.ge.lat_c(ilat1)) ilat=ilat1
       ilatloc_c(i)=ilat
    enddo
    do i=ilatmax+1,ilatdeg
       ilatloc_c(i)=0
    enddo
    idep=1
    do i=1,idepmax
       idep1=idep+1
!      print*,idep1
       dep_now=float(i)*bld4-dep1_c
!      print*,dep_now,dep_c(idep1)
       if (dep_now.ge.dep_c(idep1)) idep=idep1
       ideploc_c(i)=idep
    enddo
    do i=idepmax+1,idepkm
       ideploc_c(i)=0
    enddo
!    return

    !-for pn velocity 2d-map
    lon1_m=bld1-lon_m(1)
    ilonmax=(1e-10)+(lon_m(nlon_m)+lon1_m)/bld1
    lat1_m=bld1-lat_m(1)
    ilatmax=(1e-10)+(lat_m(nlat_m)+lat1_m)/bld1
!   print*,"lon1_m,lat1_m,ilonmax,ilatmax,bld1:",lon1_m,lat1_m,ilonmax,ilatmax,bld1
    if ((ilonmax.gt.ilondeg).or.(ilatmax.gt.ilatdeg)) then
       print*,"Error, model dimension out of range!"
       stop
    endif
    ilon=1
!   open(6,file="test2.out")
    do i=1,ilonmax
       ilon1=ilon+1
       lon_now=float(i)*bld1-lon1_m
       if (lon_now.ge.lon_m(ilon1)) ilon=ilon1
       ilonloc_m(i)=ilon
!      write(6,*)lon_now,ilonloc_m(i)
    enddo
    do i=ilonmax+1,ilondeg
       ilonloc_m(i)=0
    enddo
    ilat=1
    do i=1,ilatmax
       ilat1=ilat+1
       lat_now=float(i)*bld1-lat1_m
       if (lat_now.ge.lat_m(ilat1)) ilat=ilat1
       ilatloc_m(i)=ilat
!      write(6,*)lat_now,ilatloc_m(i)
    enddo
    do i=ilatmax+1,ilatdeg
       ilatloc_m(i)=0
    enddo
    return

    !-for moho map
    lon1_moho=bld2-lon_moho(1)
!   print*,nlon
!   print*,lon1,lon_mod(nlon)
!   print*,(lon_mod(i),i=1,nlon)
    ilonmax=(1e-10)+(lon_moho(nlon_moho)+lon1_moho)/bld2
!   print*,ilonm
!   stop
    lat1_moho=bld2-lat_moho(1)
    ilatmax=(1e-10)+(lat_moho(nlat_moho)+lat1_moho)/bld2
    if ((ilonmax.gt.ilondeg).or.(ilatmax.gt.ilatdeg)) then
       print*,"Error, moho depth model dimension out of range!"
       stop
    endif
    ilon=1
!   open(6,file="test2.out")
    do i=1,ilonmax
       ilon1=ilon+1
       lon_now=float(i)*bld2-lon1_moho
       if (lon_now.ge.lon_moho(ilon1)) ilon=ilon1
       ilonloc_moho(i)=ilon
!      write(6,*)lon_now,ilonloc(i)
    enddo
    do i=ilonmax+1,ilondeg
       ilonloc_moho(i)=0
    enddo
    ilat=1
    do i=1,ilatmax
       ilat1=ilat+1
       lat_now=float(i)*bld2-lat1_moho
       if (lat_now.ge.lat_moho(ilat1)) ilat=ilat1
       ilatloc_moho(i)=ilat
!      write(6,*)lat_now,ilatloc(i)
    enddo
    do i=ilatmax+1,ilatdeg
       ilatloc_moho(i)=0
    enddo
end subroutine


subroutine intmap_2d(lon,lat,dxy,ip,jp)
    implicit real*8 (a-h,o-z)
    include "para.inp"
    real*8 lon,lat,dxy
    ip=int(1e-10+(lon+lon1_m)/dxy)
    jp=int(1e-10+(lat+lat1_m)/dxy)
!    print*,lon,lat,dxy
!    print*,lon1_m,lat1_m
!    print*,ip,jp
    if ((ip.le.0).or.(jp.le.0)) then
       print*,"Error,lon,lat out of range!"
       print*,"lon=",lon,"lat=",lat
       print*,"ip,jp",ip,jp
       stop
    endif
    ip=ilonloc_m(ip)
    jp=ilatloc_m(jp)
!    print*,ip,jp
!    pause
    if ((ip.eq.0).or.(jp.eq.0)) then
       print*,"Error,lon,lat out of range!"
       print*,"lon=",lon,"lat=",lat
       print*,"ip,jp",ip,jp
       stop
    endif
    return
end subroutine


subroutine intmap_3d(lon,lat,dep,dxy,dz,ip,jp,kp)
    implicit integer (i-n)
    include "para.inp"
    real*8 lon,lat,dep,dxy,dz
    ip=int(1e-10+(lon+lon1_c)/dxy)
    jp=int(1e-10+(lat+lat1_c)/dxy)
    kp=int(1e-10+(dep+dep1_c)/dz)
    if ((ip.le.0).or.(jp.le.0).or.(kp.le.0)) then
       print*,"Error,lon,lat,dep out of range!"
       print*,"lon=",lon,"lat=",lat,"dep=",dep
       print*,"ip,jp,kp",ip,jp,kp
       print*,'intmap_3d'
       stop
    endif
    ip=ilonloc_c(ip)
    jp=ilatloc_c(jp)
    kp=ideploc_c(kp)
    if ((ip.eq.0).or.(jp.eq.0).or.(kp.eq.0)) then
       print*,"Error,crust lon,lat out of range!"
       print*,"lon=",lon,"lat=",lat,"dep=",dep
       print*,"ip,jp,kp",ip,jp,kp
       print*,'intmap_3d 1'
       stop
! Veiw only
!      else
!       print*,"ip,jp,kp",ip,jp,kp
    endif
    return
end subroutine


subroutine vel_3d(lon,lat,dep,v)
    implicit integer (i-n)
    include "para.inp"
    real*8 lon,lat,dep,v
    real*8 lonf,lonf1,latf,latf1,depf,depf1
    real*8 wv(2,2,2)
    common /weight/ wv,ip,jp,kp

    !-longitude spherical correction
    if (lon>180.0) then
       lon=lon-360.0
    elseif (lon<-180.0) then
       lon=lon+360.
    else
       continue
    endif
    !-latitude spherical correction
    if (lat>90.0.or.lat<-90.0) then
       if (lat>90.0) then
          lat=90.0-(abs(lat)-90.0)
       elseif (lat<-90.0) then
          lat=-90.0+(abs(lat)-90.0)
       else
          continue
       endif
       !-if latitude over range (-90.0, 90.0), longitude opposite
       if (lon<0.0) then
          lon=lon+180.0
       elseif (lon>0.0) then
          lon=lon-180.0
       else
          continue
       endif
    endif
    !-depth spherical correction
    if (dep>ro) dep=2*ro-dep
    call intmap_3d(lon,lat,dep,bld3,bld4,ip,jp,kp)

    ip1=ip+1
    jp1=jp+1
    kp1=kp+1
    if ((ip1.gt.nlon_c).or.(jp1.gt.nlat_c).or.(kp1.gt.ndep_c)) then
       print*,"Error, ip1,jp1,kp1 out of range!"
       print*,"ip1=",ip1,"jp1=",jp1,"kp1=",kp1
       stop
    endif
    lonf=(lon-lon_c(ip))/(lon_c(ip1)-lon_c(ip))
    latf=(lat-lat_c(jp))/(lat_c(jp1)-lat_c(jp))
    depf=(dep-dep_c(kp))/(dep_c(kp1)-dep_c(kp))
    lonf1=1.0-lonf
    latf1=1.0-latf
    depf1=1.0-depf
    wv(1,1,1)=lonf1*latf1*depf1
    wv(2,1,1)=lonf*latf1*depf1
    wv(1,2,1)=lonf1*latf*depf1
    wv(2,2,1)=lonf*latf*depf1
    wv(1,1,2)=lonf1*latf1*depf
    wv(2,1,2)=lonf*latf1*depf
    wv(1,2,2)=lonf1*latf*depf
    wv(2,2,2)=lonf*latf*depf
    if (ivmod.eq.0) then
       if (vs(ip,jp,kp)==0.0) print*,"S wave velocity model loading error!"
       v=wv(1,1,1)*vs(ip,jp,kp)+wv(2,1,1)*vs(ip1,jp,kp)&
        +wv(1,2,1)*vs(ip,jp1,kp)+wv(2,2,1)*vs(ip1,jp1,kp)&
        +wv(1,1,2)*vs(ip,jp,kp1)+wv(2,1,2)*vs(ip1,jp,kp1)&
        +wv(1,2,2)*vs(ip,jp1,kp1)+wv(2,2,2)*vs(ip1,jp1,kp1)
    endif
    if (ivmod.eq.1) then
       if (vp(ip,jp,kp)==0.0) print*,"P wave velocity model loading error!"
       v=wv(1,1,1)*vp(ip,jp,kp)+wv(2,1,1)*vp(ip1,jp,kp)&
        +wv(1,2,1)*vp(ip,jp1,kp)+wv(2,2,1)*vp(ip1,jp1,kp)&
        +wv(1,1,2)*vp(ip,jp,kp1)+wv(2,1,2)*vp(ip1,jp,kp1)&
        +wv(1,2,2)*vp(ip,jp1,kp1)+wv(2,2,2)*vp(ip1,jp1,kp1)
    endif
    return
end subroutine


subroutine vel_pn(lon,lat,v)
    implicit real*8 (a-h,o-z)
    include "para.inp"
    real*8 lon,lat,v
    real*8 lonf,lonf1,latf,latf1
    real*8 wv(2,2,2)
    common /weight/ wv,ip,jp,kp
!    v=8.035  ! use constant Pn velocity derived from 1D inversion
!    return
    call intmap_2d(lon,lat,bld1,ip,jp)
    ip1=ip+1
    jp1=jp+1
!    print*,"lon,lat",lon,lat
!    print*,ip,jp
    !-check if surrounding grid within the range
    if ((ip1.gt.nlon_m).or.(jp1.gt.nlat_m)) then
       print*,"Error, ip1,jp1 out of range!"
       print*,"ip1=",ip1,"jp1=",jp1
       stop
    endif
    lonf=(lon-lon_m(ip))/(lon_m(ip1)-lon_m(ip))
    latf=(lat-lat_m(jp))/(lat_m(jp1)-lat_m(jp))
    lonf1=1.0-lonf
    latf1=1.0-latf
    wv(1,1,1)=lonf1*latf1
    wv(2,1,1)=lonf*latf1
    wv(1,2,1)=lonf1*latf
    wv(2,2,1)=lonf*latf
    v=wv(1,1,1)*vpn(ip,jp)+wv(2,1,1)*vpn(ip1,jp)&
      +wv(1,2,1)*vpn(ip,jp1)+wv(2,2,1)*vpn(ip1,jp1)
!    print*,"wv(1,1,1)",wv(1,1,1)
!    print*,"wv(2,1,1)",wv(2,1,1)
!    print*,"wv(1,2,1)",wv(1,2,1)
!    print*,"wv(2,2,1)",wv(2,2,1)
!    print*,"v=",v
!    pause
    return
end subroutine


!subroutine depth(lon,lat,dm)
!    implicit real*8 (a-h,o-z)
!    include "para.inp"
!    real*8 lon,lat,dm
!    real*8 lonf,lonf1,latf,latf1
!    real*8 wt_moho(2,2)
!    common /wt/ wt_moho,ip_moho,jp_moho

!    call intmap_2d(lon,lat,bld2,ip,jp)
!    ip_moho=ip
!    jp_moho=jp
!    ip1_moho=ip_moho+1
!    jp1_moho=jp_moho+1
!    print*,"lon=",lon,"lat=",lat
!    print*,"ip_moho=",ip_moho,"jp_moho=",jp_moho
!    pause
    !-check if surrounding grid within the range
!    if ((ip1_moho.gt.nlon_moho).or.(jp1_moho.gt.nlat_moho)) then
!       print*,"Error, ip1_moho,jp1_moho out of range!"
!       print*,"ip1_moho=",ip1_moho,"jp1_moho=",jp1_moho
!       print*,lon,lat
!       stop
!    endif
!    lonf=(lon-lon_moho(ip_moho))/(lon_moho(ip1_moho)-lon_moho(ip_moho))
!    latf=(lat-lat_moho(jp_moho))/(lat_moho(jp1_moho)-lat_moho(jp_moho))
!    lonf1=1.0-lonf
!    latf1=1.0-latf
!    wt_moho(1,1)=lonf1*latf1
!    wt_moho(2,1)=lonf*latf1
!    wt_moho(1,2)=lonf1*latf
!    wt_moho(2,2)=lonf*latf
!    print*,"wt_moho(1,1)=",wt_moho(1,1)
!    print*,"wt_moho(2,1)=",wt_moho(2,1)
!    print*,"wt_moho(1,2)=",wt_moho(1,2)
!    print*,"wt_moho(2,2)=",wt_moho(2,2)
        
!    dm=wt_moho(1,1)*mhd(ip_moho,jp_moho)+wt_moho(2,1)*mhd(ip1_moho,jp_moho)&
!       +wt_moho(1,2)*mhd(ip_moho,jp1_moho)+wt_moho(2,2)*mhd(ip1_moho,jp1_moho)
!    print*,"dm=",dm
!    pause
!end subroutine

    subroutine output1(x0,y0,z0,dt,m)
      implicit integer (i-n)
      include "para.inp"
      real*8 x0,y0,z0,dx,dy,dt,ttnew,m(4),erh,erz
      integer date,time
      character*2 phase
  
      call ot_corr(dt,date,time)
  
      call cal_delta(y0,x0-0.5,y0,x0+0.5,dx)
      call cal_delta(y0-0.5,x0,y0+0.5,x0,dy)
      erz=abs(m(3))
      erh=sqrt((m(1)*dy)**2+(m(2)*dx)**2)  
     if(abs(dt).le.2.0) then
      write(3,10)date,time,y0,x0,z0,evt_mag(evtnow),rms,erh,erz,evt_id(evtnow),&
                 (evt_loc(evtnow,j),j=1,3) 
      write(4,"(1a,i7)")"#",evt_id(evtnow) 
      write(7,"(1a,i7)")"#",evt_id(evtnow) 
      do i=1,nrow
         !-wave phase
         if (ray_pha(ray_idx(i))==1) phase="Pg"
         if (ray_pha(ray_idx(i))==2) phase="Sg"
         if (ray_pha(ray_idx(i))==3) phase="Pn"
         if (ray_pha(ray_idx(i))==4) phase="Sn"
         !-new travel time
         ttnew=abs_tt(ray_idx(i))-dt
         if (ttnew<=0.) then
            print*,"travel time error!!"
            write(4,30)sta_nm(sta_idx(ray_idx(i))),ttnew,ray_wei(ray_idx(i)),phase,"!"
         else
            write(4,20)sta_nm(sta_idx(ray_idx(i))),ttnew,ray_wei(ray_idx(i)),phase
            !write(7,20)sta_nm(sta_idx(ray_idx(i))),ttnew,ray_wei(ray_idx(i)),phase
            write(7,22)sta_nm(sta_idx(ray_idx(i))),ttnew,abs_tt(ray_idx(i)),dt,ray_wei(ray_idx(i)),phase
         endif
      enddo
    else
       print*,"time-shift greater than 2.0 sec, remove this event"
    endif
  10  format(2(i9,1x),2(f8.3,1x),f7.2,1x,4(f4.1,1x),i6,1x,2(f8.3,1x),f7.2)
  20  format(1x,a5,f12.3,2x,i1,2x,a2)
  22  format(1x,a5,3(f12.3,2x),i1,2x,a2)
  30  format(1x,a5,f12.2,2x,i1,2x,a2,1x,a1)
  
    !-update hypo info for next correction round
    evt_loc(evtnow,2)=x0
    evt_loc(evtnow,1)=y0
    evt_loc(evtnow,3)=z0
    evt_date(evtnow)=date
    evt_time(evtnow)=time
end subroutine


!=====
! output results in phase-format
!=====
subroutine output2(x0,y0,z0,dt,m)
    implicit integer (i-n)
    include "para.inp"
    character*90 card(maxsta)
    integer date,time,yy,mm,dd,hr,min,ilat,ilon,gap
    real*8 sec,lat,lon,dep,ttnew(maxsta)
    real*8 x0,y0,z0,erz,erh,dx,dy,dt,m(4)

    call ot_corr(dt,date,time)  ! correct original time
    
    !-decode original date & time
    yy=int(date/10000)
    mm=int((date-yy*10000)/100)
    dd=date-yy*10000-mm*100
    hr=int(time/1000000)
    min=int((time-hr*1000000)/10000)
    sec=real(time-hr*1000000-min*10000)/100.
!    print*,yy,mm,dd,hr,min,sec

    !-decode hypocenter location
    ilat=int(y0)
    ilon=int(x0)
    lat=(y0-ilat)*60.
    lon=(x0-ilon)*60.
    dep=z0

    call cal_delta(y0,x0-0.5,y0,x0+0.5,dx)
    call cal_delta(y0-0.5,x0,y0+0.5,x0,dy)
    erz=abs(m(3))
    erh=sqrt((m(1)*dy)**2+(m(2)*dx)**2)

    do i=1,nrow
       ttnew(i)=abs_tt(ray_idx(i))-dt+sec
       if (ttnew(i)<=0.) stop "travel time error!!"
    enddo
    imin=0
    if (minval(ttnew)>60.0) imin=int(minval(ttnew)/60.0)

    ista=0
    do i=1,nrow
!       print*,ray_pha(ray_idx(i))
       do j=1,ista
          if (sta_nm(sta_idx(ray_idx(i)))==card(j)(2:5)) then
             if (ray_pha(ray_idx(i))==1) write(card(j)(24:39),"(f6.2,2f5.2)")ttnew(i)-real(imin)*60.,res(i),&
                                               real(ray_wei(ray_idx(i)))
             if (ray_pha(ray_idx(i))==2) write(card(j)(40:55),"(f6.2,2f5.2)")ttnew(i)-real(imin)*60.,res(i),&
                                               real(ray_wei(ray_idx(i)))
             goto 99
          endif
       enddo
       ista=ista+1
       if (ray_pha(ray_idx(i))==1) write(card(ista),20)sta_nm(sta_idx(ray_idx(i))),epid(i),nint(azi(i)),nint(tkofag(i)),min+imin,&
                                         ttnew(i)-real(imin)*60.,res(i),real(ray_wei(ray_idx(i)))
       if (ray_pha(ray_idx(i))==2) write(card(ista),30)sta_nm(sta_idx(ray_idx(i))),epid(i),nint(azi(i)),nint(tkofag(i)),min+imin,&
                                         ttnew(i)-real(imin)*60.,res(i),real(ray_wei(ray_idx(i)))
99     continue
    enddo

    call get_gap(nint(azi(1:nrow)),nrow,gap)
    if (nrow>99) nrow=99
    write(3,40)date,time,y0,x0,z0,evt_mag(evtnow),rms,erh,erz,evt_id(evtnow)
    write(4,10)yy,mm,dd,hr,min,sec,ilat,lat,ilon,lon,dep,evt_mag(evtnow),ista,minval(epid(1:nrow)),gap,rms,erh,erz,nrow
    do i=1,ista
       write(4,*)card(i)
    enddo
10  format(1x,i4,4i2,f6.2,i2,f5.2,i3,f5.2,f6.2,f4.2,i2,f5.1,i3,f4.2,2f4.1,2x,i3)
20  format(1x,a5,f5.1,2i4,1x,i3,f6.2,2f5.2)
30  format(1x,a5,f5.1,2i4,1x,i3,16x,f6.2,2f5.2)
40  format(2(i9,1x),2(f8.3,1x),f7.2,1x,4(f4.1,1x),i6)

    !-update hypo info for next correction round
    evt_loc(evtnow,2)=x0
    evt_loc(evtnow,1)=y0
    evt_loc(evtnow,3)=z0
    evt_date(evtnow)=date
    evt_time(evtnow)=time
end subroutine


!=====
! damped and weighted least square inversion for generalized inverse
! written by Hsin-Hua Huang at 11/09/2010
!=====
subroutine LS_inv(G,m,d,Wd,dp,nd,nm,GG)
    include "para.inp"
    integer nd,nm
    real*8 m(4),d(nd),G(maxsta,4),GtG(4,4),Gtd(4),Wd(maxsta,maxsta)
    real*8 dp,GG(4,4)

    !--derive matrix GT*W*G
    do i=1,nm
       do j=1,nm
          GtG(i,j)=0.0
          do k=1,nd
             GtG(i,j)=GtG(i,j)+G(k,i)*Wd(k,k)*G(k,j) !--Wd is a weighting matrix to the data norm
          enddo
       enddo
       !--damping factor/avoid the matrix being a singular matrix 
       GtG(i,i)=GtG(i,i)+dp**2
    enddo
  
    !--inverse GT*G to [GT*G]^-1
    call matrix_inv(GtG,nm)
    do i=1,nm
!      write(*,"(i5,4g20.10)")i,GtG(i,1),GtG(i,2),GtG(i,3),GtG(i,4)
    enddo

    !--derive matrix GT*d
    do i=1,nm
       Gtd(i)=0.0
       do j=1,nd
          Gtd(i)=Gtd(i)+G(j,i)*Wd(j,j)*d(j)
       enddo
    enddo
        
    !--m = ([GT*G]^-1)*(GT*d)  
    do i=1,nm
       m(i)=0.0
       do j=1,nm
          m(i)=m(i)+GtG(i,j)*Gtd(j)
       enddo
    enddo
end subroutine


!=====
! derive the inverse matrix
! written by Hsin-Hua Huang at 11/11/2010
!=====
subroutine matrix_inv(A,n)
  integer r(n-1),n
  real*8 A(n,n),X(n,n)
  real*8 dia,multi(n-1)

  !--creat a unit matrix for inverse
  X=0.0
  do i=1,n
    X(i,i)=1.0
  end do

  r=0
  do i=1,n
    do k=1,n-1
      r(k)=mod(i+k,n)
      if(r(k)==0) r(k)=n
      multi(k)=A(r(k),i)
    enddo

    !--check if any diagonal element is 0 during process
    if(A(i,i)==0.0)then  !-if yes, sum the other rows up to be non-zero
      do k=1,n-1
        do j=1,n
          !--change sign to avoid the sum of the other rows being zero
          if(A(r(k),i)<0.0) then
            A(i,j)=A(i,j)-A(r(k),j)
            X(i,j)=X(i,j)-X(r(k),j)
          else
            A(i,j)=A(i,j)+A(r(k),j)
            X(i,j)=X(i,j)+X(r(k),j)
          end if
        enddo
      enddo

      !--clean the sum value at ahead elements
      do k=1,i-1                                        
        if(A(i,i-k)<0.0)then
          do j=1,n
            A(i,j)=A(i,j)+A(i-k,j)
            X(i,j)=X(i,j)+X(i-k,j)
          enddo
        else
          do j=1,n
            A(i,j)=A(i,j)-A(i-k,j)
            X(i,j)=X(i,j)-X(i-k,j)
          enddo
        end if
      enddo
    end if

    dia=A(i,i)
    do j=1,n
      A(i,j)=A(i,j)/dia  !-diagonal normalized
      X(i,j)=X(i,j)/dia
    end do
  
    do k=1,n-1
      do j=1,n
        A(r(k),j)=A(r(k),j)-A(i,j)*multi(k)  !-make the other rows zero on the ith column
        X(r(k),j)=X(r(k),j)-X(i,j)*multi(k)
      enddo
    end do
  end do
  
  A=X  !-return inverse matrix to main program
  return
end subroutine

     subroutine distazi(coevla,evtlon,costla,stalon,dist,az)
      implicit real(a-h,o-z)
      parameter(eps=1.0e-8)
      parameter(pi=3.14159265359)
      parameter(rad_per_deg=0.0174532925199432955)
      real*8 coevla,evtlon,costla,stalon,dist,az
      real*8 colat0, colat1, lon0, lon1
      real del,cosdel, tmp, azimuth, cosaz, sinaz
!
      colat0=coevla*rad_per_deg
      colat1=costla*rad_per_deg
      lon0=evtlon*rad_per_deg
      lon1=stalon*rad_per_deg
      !rad_per_deg = pi/180.
!     Distance calculation
      cosdel=cos(colat0)*cos(colat1)+sin(colat0)*sin(colat1)*cos(lon1-lon0)
      if(cosdel.gt.1.0) cosdel=cosdel-eps
      if(cosdel.lt.-1.0) cosdel=cosdel+eps
      del=acos(cosdel)
      dist=del/rad_per_deg
!     Azimuth calculation
      tmp=sin(colat0)*sin(del)
      if(tmp.le.eps) then
         if(coevla.le.eps ) az=180
         if(coevla.ge.(180.-eps )) az=0
         if(dist.le.eps) az=-999.0
         if(dist.ge.180.0-eps) az=-999.0
      endif
      
      cosaz=(cos(colat1)-cos(colat0)*cos(del))/tmp
      sinaz=sin(colat1)*sin(lon1-lon0)/sin(del)
      azimuth=atan2(sinaz,cosaz)/rad_per_deg
      if (azimuth.lt.0.) azimuth = azimuth+360.
      az=azimuth


     end subroutine

     subroutine epidist(ey,ex,ez,sy,sx,sz,dis,azi)
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
 !              azi=azimuth(coelat,ex,coslat,sx) 
      
      
            end subroutine
