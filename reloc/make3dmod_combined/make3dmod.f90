        program make3mod
!     ==================================================================================
!     Make the 3 dimentional velocity model from input 1 dimensional model
!     Output use for 3D spherical relocation package
!     The original idea from Hsin-Hua Huang
!     2020/01/01: Transfer to fortran 90 language and change making for both P&S models
!     ================================================================================== 
      parameter (maxgrid=100)
      character*32 fnp,fns
      character*80 line
      integer nlon,nlat,ndep,nmod
      integer nlon1,nlat1,ndep1,nmod1
      real lon(maxgrid),lat(maxgrid),dep(maxgrid)
      real lon1(maxgrid),lat1(maxgrid),dep1(maxgrid)
      real lon_min,lon_max,lat_min,lat_max
      real lon_min1,lon_max1,lat_min1,lat_max1
      real vpm,dm,vsm,dm1        
      real bldz,bldxy,span,vmi(maxgrid,maxgrid)
      real bldz1,bldxy1,span1,vmi1(maxgrid,maxgrid)
      real v3d(maxgrid,maxgrid,maxgrid),v1d(maxgrid)
      real v3ds(maxgrid,maxgrid,maxgrid),v1ds(maxgrid)

      print*,"This version make the model for P, S or both P-S wave speed"
      write(*,*)"Reading the paramter file 3dmod.params"
!   Read the parameter file
      open(999,file='3dmod.params',status='old')
      l=0
      do
        read(999,'(a)',iostat=ios)line
!        print*,line
        if(ios.lt.0) exit
        if((line(1:1).eq."*").or.(line(1:1).eq."#").or.(line(1:1).eq."!")) cycle
        l=l+1
        if(l.eq.1) read(line,*)lon_min,lon_max
        if(l.eq.2) read(line,*)lat_min,lat_max
        if(l.eq.3) read(line,*)bldxy
        if(l.eq.4) read(line,*)span
        if(l.eq.5) read(line,*)nmod
        if(l.eq.6) read(line,'(a)')fnp
        if(l.eq.7) read(line,'(a)')fns
      enddo
      !   Calculate from paramter
      nlon=(lon_max-lon_min)/bldxy+3
      nlat=(lat_max-lat_min)/bldxy+3
      write(*,*)"nlon,nlat",nlon,nlat

      lon(1)=lon_min-span
      lon(nlon)=lon_max+span
      lat(1)=lat_min-span
      lat(nlat)=lat_max+span

!      
      do i=2,nlon-1
        lon(i)=lon_min+(i-2)*bldxy
      enddo

      do i=2,nlat-1
        lat(i)=lat_min+(i-2)*bldxy
      enddo

      close(999)
!      write(*,*)lon_min,lon_max,lat_min,lat_max,bldxy,span,nmod,fnp,fns
      fnp=trim(fnp)
      fns=trim(fns)
!   Check file exsiting
      if(nmod.eq.3) then
      call existf(fnp)
      call existf(fns)
      goto 1
      elseif(nmod.eq.1) then
        call existf(fnp)
        goto 1
      elseif(nmod.eq.2) then
        call existf(fns)
        goto 2
      else
        write(*,*)"nmod input wrong!!" 
        stop
      endif

!
!   - Reading the phases
!   Read P-wave velocity model
 1    write(*,*)"read 1D model for P wave in file: ",fnp
      open(1,file=fnp,status='old')     
      read(1,*)bldz,ndep
      read(1,*)(dep(i),i=1,ndep)
      read(1,*)(v1d(i),i=1,ndep)
      read(1,*)vpm ! v_moho
!       print*,vm
      read(1,*)dm ! d_moho
!   Condition to skip read S model
!   Buiding the 3D model for P wave
! make grid matrix
      do k=1,ndep
        do j=1,nlat
            do i=1,nlon
                v3d(i,j,k)=v1d(k)
            enddo
        enddo
      enddo
 ! insert moho velocity layer = nlon x nlat grid
      do j=1,nlat
        do i=1,nlon
            vmi(i,j)=vpm
        enddo
      enddo
      close(1)
!   If nmod = 1 then not read the S wave model. goto 
      if(nmod.eq.1) then
      goto 3  
      endif
!   Read S-wave velocity model
 2    write(*,*)"read 1D model for S wave in file: ",fns
      open(2,file=fns,status='old')
      read(2,*)bldz1,ndep1
      read(2,*)(dep1(i),i=1,ndep1)
      read(2,*)(v1ds(i),i=1,ndep1)
      read(2,*)vsm
      read(2,*)dm1
!   
!   Buiding the 3D model for S-wave
! make grid matrix
      do k=1,ndep1
        do j=1,nlat
            do i=1,nlon
                v3ds(i,j,k)=v1ds(k)
            enddo
        enddo
      enddo
!         print*,v3ds(1,2,3)
! insert moho velocity layer = nlon x nlat grid
      do j=1,nlat
        do i=1,nlon
            vmi1(i,j)=vsm
        enddo
      enddo  
      close(2)  
!        print*,vmi1(3,4)

! output model file
 3    continue
      print*,"3D models build up"
      print*,"Write the output model on ...MOD3D file"
      open(11,file="MOD3D")
!   Model dimension:
      if(nmod.eq.1) then
        print*,"writing 3D model for the P wave only"
        write(11,100)bldxy,nlon,nlat ! moho grid
        write(11,101)bldxy,bldz,nlon,nlat,ndep ! crust grid
        write(11,110)(lon(i),i=1,nlon)
        write(11,110)(lat(i),i=1,nlat)
        write(11,110)(lon(i),i=1,nlon)
        write(11,110)(lat(i),i=1,nlat)
        write(11,110)(dep(i),i=1,ndep) 
       
        do j=1,nlat
            write(11,120)(vmi(i,j),i=1,nlon)
        enddo

        do k=1,ndep
            do j=1,nlat
!                write(1,120)(vmi(i,j),i=1,nlon)
                write(11,120)(v3d(i,j,k),i=1,nlon)
            enddo
        enddo
        goto 9999
!   Great, it work
      elseif(nmod.eq.2) then
        print*,"writing 3D model for the S wave only"
!
        write(11,100)bldxy,nlon,nlat ! moho grid
        write(11,101)bldxy,bldz1,nlon,nlat,ndep1 ! crust grid
        write(11,110)(lon(i),i=1,nlon)
        write(11,110)(lat(i),i=1,nlat)
        write(11,110)(lon(i),i=1,nlon)
        write(11,110)(lat(i),i=1,nlat)
        write(11,110)(dep(i),i=1,ndep1)
        do jj=1,nlat
            write(11,120)(vmi1(ii,jj),ii=1,nlon)
        enddo
!
        do kk=1,ndep
            do jj=1,nlat
                write(11,120)(v3ds(ii,jj,kk),ii=1,nlon)
            enddo
        enddo
        goto 9999
      elseif(nmod.eq.3) then
!   Check the equal of 2 model
        if((ndep.ne.ndep1).or.(dm.ne.dm1)) then
            print*,"The grid of P and S model not equal, exit"
            print*,"ndepp: ",ndep," ndeps: ",ndep1
            print*,"nlong_p: ",nlon," nlong_s: ",nlon1
            print*,"nlat_s: ",nlat," nlat_s: ",nlat1
            print*,"moho depth_p: ",dm," homo depth_s: ",dm1
            stop
        else
            write(11,100)bldxy,nlon,nlat ! moho grid
            write(11,101)bldxy,bldz1,nlon,nlat,ndep1 ! crust grid
            write(11,110)(lon(i),i=1,nlon)
            write(11,110)(lat(i),i=1,nlat)
            write(11,110)(lon(i),i=1,nlon)
            write(11,110)(lat(i),i=1,nlat)
            write(11,110)(dep(i),i=1,ndep)

            print*,"writing the P wave model"
!   Pn wave part
            do j=1,nlat
                write(11,120)(vmi(i,j),i=1,nlon)
            enddo
!   Pg wave part
            do k=1,ndep
                do j=1,nlat
                    write(11,120)(v3d(i,j,k),i=1,nlon)
                enddo
            enddo
!
            print*,"writing the S wave model"
!                print*,"vmi",vmi(1,1),"vmi1",vmi1(1,1)
!   Sn wave part
            do jj=1,nlat
                write(11,120)(vmi1(ii,jj),ii=1,nlon)
            enddo
!   Sg wave part
!                print*,"v3ds",v3ds(1,2,3)
            do kk=1,ndep
                do jj=1,nlat
                    write(11,120)(v3ds(ii,jj,kk),ii=1,nlon)
                enddo
            enddo
        endif
        goto 9999
      else
        print*,"Error, number of model out of program contribution"
        stop
      endif
      
 9999 close (11)
      write(*,*)"Done!"
!   Wirte 2 models
100     format(f4.2,1x,i2,1x,i2)
101     format(f4.2,1x,f4.2,1x,i3,1x,i3,1x,i3)
110     format(100f7.2)
120     format(100f8.3)
        end program

      subroutine existf(fn)
!   Subroutine wrote by HHHuang
        character*32 fn
        logical here
        inquire(file=fn,exist=here)
        if (.not.here) then
           print*,'file ',trim(fn),' does not exist! please check the file'//&
                 ' name in 3dmod.params'
           stop
        endif 
      end subroutine

