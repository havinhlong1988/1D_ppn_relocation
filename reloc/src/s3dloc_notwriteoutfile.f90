      program spherical_3dloc
        implicit integer (i-n)
        include "para.inp"
        real*8 x0,y0,z0,dt
        real*8 m(4),dx,dy,rms_sum(maxitr+5),sta_sum(maxsta,2)
        integer sta_count(maxsta),sta_azi(maxevt,maxsta),sta_gap,nrms(maxitr+5),nmod
        integer iunit, iiunit, ifile
        character*20 sta_fn,evt_fn,abs_fn,vel_fn,line
        character*3 label
        
!      ================= Modification history =============================
!       2019/10/24: Write out station corection loop to stacor.log for test
!       2019/10/29: Station correction may not accumulated after each itteration
!     --------------------------------------------------------------------------------------------------------       
!       2019/10/30: Fix the station correction didn't accumulated 
!       (now output of s3dloc.cor was accumulated of every itteraction)
!       2019/10/30: Change print out itteration number loop and station correction accumulated to screen
!     --------------------------------------------------------------------------------------------------------
!       2019/11/25: Force every earthquake have new depth greater than
!       moho-depth will become bad event (cos this program not considered the under-moho earthquake depth)
!     --------------------------------------------------------------------------------------------------------
!       2019/12/06: Add s3dloc.tt file to compare calculated & observated travel time (output1)
!     --------------------------------------------------------------------------------------------------------
!       2020/01/02: Add ray_path.chk file, which record the ray_tracing
!       2020/03/12: Fix ray_path.chk file record; Add the location.chk file to record location changing!
!     --------------------------------------------------------------------------------------------------------
!       2020/07/27: Modification the raypath report file to record the ray USED to relocation each iteration.
!       change the ray_path.chk to ray_path_xxx file. Only apply for icor = 1. If icor > 1, ray_path_xxx
!       are represented for last icor run!
!     -------------------------------------------------------------------------------------------------------- 
        call random_seed()  
      
        print*,"read <input> file ..."
        open(1,file='input',status='old')
        l=0
      do
        read(1,'(a20)',iostat=ios)line
        if(ios.lt.0) exit
        if((line(1:1).eq.'*').or.(line(1:1).eq.'#').or.(line(1:1).eq.'!')) cycle
        l=l+1
        if(l==1)    read(line,'(a)')evt_fn
        if(l==2)    read(line,'(a)')sta_fn
        if(l==3)    read(line,'(a)')abs_fn
        if(l==4)    read(line,'(a)')vel_fn
        if(l==5)    read(line,*)dm
        if(l==6)    read(line,*)nmod
      enddo
      close(1)
        print*,"event data file"    
        print*,evt_fn
        print*,"station data file:"
        print*,sta_fn
        print*,"absolute travel time file:"
        print*,abs_fn
        print*,"input velocity model file:"
        print*,vel_fn
        print*,"the moho depth:"
        print*,dm
        print*,"relocation using",nmod," type of phases (1:P, 2: P&S)"

        print*,"Processing ..."
        call getdata1(evt_fn,sta_fn,abs_fn)  ! Separate files of evt, sta, and tt in Haijiang format
!        call getdata2(evt_fn,sta_fn)         ! Phase file in CWB format)
        call loadvel(vel_fn,nmod)   
        print*,'. . . run 3D relocation with input of: '
        print*,nevt,"events, ",nsta," stations, and ",nray," rays"
!   History record file
        open(1,file="s3dloc.chk",status="unknown") ! Running buletin file
        open(10,file="stacor.log",status="unknown") ! Station record buletin file
    !-residual check file annotations
        write(1,*)"/rms, adjustment, tshift, x0, y0, z0/"



    !-loops of station corrections
!    icor=0
    do icor=1,maxcor
!       icor=icor+1
       print*,"===icor==>",icor
        write(1,"(a28,i4)")"==STATION CORRECTION NUMBER==",icor
        write(1,"(a)")"Station name - station correction (P-S)"
!
        write(10,'(a15,i3)')"correction number",icor
        write(10,"(a)")"Station name - station correction (P-S)"
!
        do i=1,nsta
        write(1,"(a5,2(1x,f8.4))")sta_nm(i),(sta_cor(i,j),j=1,2)
        write(10,"(a5,2(1x,f8.4))")sta_nm(i),(sta_cor(i,j),j=1,2)
        enddo
!       write(1,"(a)")"==================*========================"

       open(2,file="s3dloc.err",status="unknown") ! List of event error!
       open(3,file="s3dloc.lis",status="unknown") ! List of good relocation event
       open(4,file="s3dloc.out",status="unknown") ! List of travel time out (calculated values)
       open(5,file="s3dloc.rms",status="unknown") ! Rms of each iteration (random search included)
       open(7,file="s3dloc.tt",status="unknown")   !observed and predicted travel time to check (output1 subroutine) 
!       open(8,file='ray_path.chk',status='unknown') ! ray tracing file
       open(9,file='location.chk',status='unknown') ! Compare input and output location
!       open(10,file="stacor.log",status="unknown")
!     File for raytracing check. The file unit = 20 + iteration number. EG 1st iteration -> unit file = 21
       iunit=20;
       label="xxx"
       do ilc=1,maxitr
         write(label(1:3),"(i3.3)") ilc
         iiunit=20+ilc
         open(iiunit,file="ray_path_"//label,status='unknown')
       enddo
!       station correction record

       !-travel-time calibration by station corrections
       do i=1,nray
          j=int(mod(ray_pha(i),2))
          abs_tt(i)=abs_tt(i)-sta_cor(sta_idx(i),j) ! tt calibrate by station correction
!        print*,"stacorr:",sta_cor(sta_idx(i),j) ! loop works
          write(1,"(a6,1x,f8.4)")sta_nm(sta_idx(i)),sta_cor(sta_idx(i),j)       
       enddo
!       print*,sta_nm(sta_idx(1)),abs_tt(1) 
       sta_count=0
       nrms=0.
       rms_sum=0. 
       do ievt=1,nevt
          write(*,'(a28,i10,2x,a11,i10)'),'#------> calculating event: ',ievt,"# eventid: ",evt_id(ievt)
          evtnow=ievt
          nrds=evt_pack(ievt+1)-evt_pack(ievt)
          write(1,*)"#",evt_id(evtnow),nrds
!         print*,nrds,evt_pack(ievt)+1,evt_pack(ievt+1)
          write(9,*)"# evt:",evt_id(evtnow)," nray: ",nrds
          write(9,'(a2,2x,3f10.3)')"I:",evt_loc(evtnow,2),evt_loc(evtnow,1),evt_loc(evtnow,3)
          call ini_hypo(x0,y0,z0)  
          write(*,'(a4,f7.3,2(a6,f7.3))')"x0: ",x0,"- y0: ",(y0)," - z0: ",z0
          !-loops of relocation iterations
          itr=0
          irs=0
          dt=0.0
          epid=0.0 
          ifile=iunit+0  
          do
             itr=itr+1
             ifile = iunit+itr
          print*,"--> iteration number:",itr
!          print*,"file write out unit: ",ifile
             iflag=0
             nrow=0
             call eqkloc(x0,y0,z0,dt,m,ifile)
!             write(*,'(a2,f10.4)')"dt",dt
             write(*,"(a19,2x,a4,f7.3,2(a6,f7.3))")"+ Update location: ","x: ",x0," - y: ",y0," - z: ",z0
             write(9,'(a4,i6,a7,i6)')"it: ",itr, " nray: ",nrow
             write(9,"(3f10.3)")x0,y0,z0
             rms_sum(itr)=rms_sum(itr)+rms*rms
             nrms(itr)=nrms(itr)+1
             !-condition notes
             if (iflag==1) then
                exit
             else if (iflag==2) then          
                write(1,*)"Out-of-range convergence..",x0,y0,z0
                !-perturb earthquake location to search another local minimum as being out of range
                irs=irs+1
!                print*,"search new random location, felt in luck", irs
                if(irs>maxirs) exit     ! round limit of researching
                call lottery(x0,y0,z0)  ! applied for small damping
                print*,"* researching minimum..."
                write(1,"(a24,3g15.5)")"* researching minimum:",x0,y0,z0
                itr=0  ! reset iteration count
             else if (iflag==3) then
                write(1,*)"Rays less than 3.."
                exit
             else if (iflag==4) then
                write(1,*)"Relocated, but not improved.."
                exit
             else
                continue
             endif
          enddo
          print*,'event condition flag: ',iflag
          if (iflag.ne.1) then
             !-illed-condition events, ouput they errors`
             write(2,50)evt_date(evtnow),evt_time(evtnow),evt_loc(evtnow,1),evt_loc(evtnow,2),evt_loc(evtnow,3),evt_mag(evtnow),&
                        0.0,0.0,0.0,evt_id(evtnow),y0,x0,z0,iflag
          else
             !-output the final report of well-condition event
             call output1(x0,y0,z0,dt,m)  ! Separate files of evt, sta, and tt in Haijiang format
!            call output2(x0,y0,z0,dt,m)  ! Phase file in CWB format
             !-cumulate station delay times
             do i=1,nrow
                j=int(mod(ray_pha(i),2))
                sta_cor(sta_idx(ray_idx(i)),j)=sta_cor(sta_idx(ray_idx(i)),j)+res(i)
                sta_count(sta_idx(ray_idx(i)))=sta_count(sta_idx(ray_idx(i)))+1
                sta_azi(sta_count(sta_idx(ray_idx(i))),sta_idx(ray_idx(i)))=nint(azi(i))
             enddo
          endif
       enddo
       do i=1,maxitr
         rms_sum(i)=sqrt(rms_sum(i)/float(nrms(i)))
          write(5,'(i6,1x,i6,1x,f10.3)')i,nrms(i),rms_sum(i)
       enddo
50     format(2(i9,1x),2(f8.3,1x),f7.2,1x,4(f4.1,1x),i6,1x,2(f8.3,1x),f7.2,1x,i6)
       close(2)
       close(3)
       close(4)
       close(5)
       !-obtain station corrections
       write(10,'(a)')"sta name/ stat count num/ number of res/cor value/sta gap "
       do i=1,nsta
!        if (sta_count(i)<10) cycle
        if (sta_count(i).gt.10) then
          call get_gap(sta_azi(1:sta_count(i),i),sta_count(i),sta_gap)
!         print*,sta_count(i),"gap:",sta_gap,"azi:",sta_azi(1:sta_count(i),i)
!        if (sta_gap>=90.) cycle
          if(sta_gap.lt.90)then
          j=int(mod(ray_pha(i),2))
          sta_cor(i,j)=sta_cor(i,j)/sta_count(i)
!
!       write(10,"(a)"),"======== output loop==========="
          else
          sta_cor(i,j)=0.0000 ! gap > 90 for station cor = 0
          endif
        else
          sta_cor(i,j)=0.0000 ! ray less than 10 force sta-cor = 0 
        endif
        write(*,"(a4,1x,f8.4,i3)")sta_nm(i),sta_cor(i,j),j
        write(10,"(a14,a5)"),"sta-corr ---->",sta_nm(i)
        write(10,55)sta_nm(i),i,sta_count(i),sta_cor(i,j),sta_gap
!      Try to accumulate the station correction
        sta_sum(i,j)=sta_sum(i,j)+sta_cor(i,j)
!1       write(*,"(a4,1x,f8.4,i3)")sta_nm(i),sta_sum(i,j),j
       write(10,"(a4,1x,f8.4,i3)")sta_nm(i),sta_sum(i,j),i
       enddo
       write(10,"(a10,i3)")"--endloop-->",icor

 55     format(a4,2(1x,i4),1x,f10.4,1x,i4)!,2(1x,i4))
!       From here, station correction rewind as input of next correction
       if (icor==maxcor) then  ! program terminates
          open(6,file="s3dloc.cor",status="unknown")
          do i=1,nsta
             write(6,"(a5,2x,5(f8.3,1x))")sta_nm(i),(sta_loc(i,j),j=1,3),(sta_sum(i,k),k=1,2)
          enddo
          close(6)
          exit
       endif
      write(10,'(a)')"=========== end loop ============"
    enddo
    close(1)
    close(10)
    close(11)
    close(8)
    close(9)
    do i=1,maxitr
      print*,"close file ",iunit +i
      close(iunit+i)
    enddo
    call reloc2inv1d
end program

subroutine eqkloc(ex,ey,ez,dt,m,idp)
    implicit integer (i-n)
    include "para.inp"
    real*8 ex,ey,ez,dt,sx,sy,sz,tt,tt0,dx,dy,junk1,junk2
    real*8 m(4),d(maxsta*2),G(maxsta,4),Wd(maxsta,maxsta)
    real*8 w(3,maxtrpts),dpi,colat,v0,re,p,pave
    real*8 geog_to_geoc,geoc_to_geog,weighting,angle,azimuth
    real*8 col_norm(4)
    real*8 ra(maxnorm)
    integer ja(maxnorm),na(maxrds*2)
    integer ncln,itnlim,leniw,lenrw,nout
    real*8 atol,btol,conlim
    logical wantse
    real*8 x(maxvar),v(maxvar),ww(maxvar)
    real*8 se(maxvar),dtt
    real*8 rw(1),ddt,ddt1
    integer iw(1),idp
    character*2 ph_ph
    common /matrix/ ra,ja,na
    dpi=asin(1.)/90.
!     20200727: Add "id" to denoted the variable of write out raypath file id
    !-initialize matrixes
    G=0.0
    d=0.0
    Wd=0.0
    col_norm=0.0
    adj=adjcut+1e-10
    ddt=0.0
    ddt1=0.0
    
    irow=0
!    ey=geog_to_geoc(ey)
!    print*,"evt_pack(evtnow)+1,evt_pack(evtnow+1)",evt_pack(evtnow)+1,evt_pack(evtnow+1)
    do i=evt_pack(evtnow)+1,evt_pack(evtnow+1)
!       if (ray_pha(i).ne.1) cycle      ! only use crutal P
!      if ((abs_tt(i)-dt)<=0.0) cycle  ! wrong time shifting, skip this ray

       ivmod=int(mod(ray_pha(i),2))    ! distinguish which velocity model should be used
       sy=sta_loc(sta_idx(i),1)
       sx=sta_loc(sta_idx(i),2)
       sz=sta_loc(sta_idx(i),3)
!       sy=geog_to_geoc(sy)
!      print*,"evt:",ex,ey,ez,sta_nm(sta_idx(i)),sx,sy,sz
       iskip=0
!       dm=32.321
       if (ray_pha(i)==1) ph_ph="Pg"
       if (ray_pha(i)==2) ph_ph="Sg"
       if (ray_pha(i)==3) ph_ph="Pn"
       if (ray_pha(i)==4) ph_ph="Sn"

       if (ray_pha(i)==1.or.ray_pha(i)==2) then  ! Pg & Sg
          call pbr(ey,ex,ez,sy,sx,sz,w,np,tt)
          dtt=abs(abs_tt(i)-tt)
          if (abs(dtt).ge.phcut) cycle
!          if(itr==maxitr) then
             write(idp,'(a1,2i5,4x,a2)')"#",i,np,ph_ph
!             write(*,*)"#",idp,ph_ph
           do iii=1,np
             write(idp,'(3f10.4,3x,i1)')(w(j,iii),j=2,3),(w(1,iii)),ray_pha(i)
           enddo
!          endif
!          do j=1,np
!            call depth(w(3,j),w(2,j),dm)
!             if (dm<=w(1,j)) then
!                iskip=1     ! skip the ray pass down to moho
!                exit
!             endif
!          enddo
       else if (ray_pha(i)==3.or.ray_pha(i)==4) then  ! Pn & Sn
!         call depth(ex,ey,dm)
          if (dm<=ez) then
             iskip=1        ! skip the ray as hypo below moho
             exit
          endif
          call pn_path(ex,ey,ez,sx,sy,sz,w,np,tt)
          dtt1=abs(abs_tt(i)-tt)
          if (abs(dtt1).ge.phcut) cycle
!          if(itr==maxitr) then
            write(idp,'(a1,2i5,4x,a2)')"#",i,np,ph_ph
!            write(*,*)"#",idp,ph_ph
          do ii=1,np
            write(idp,'(3f10.4,3x,i1)')(w(j,ii),j=2,3),(w(1,ii)),ray_pha(i)
          enddo
!         endif
       else
          stop "Phase out of the ray-tracing list!!"
       endif

!       write(*,270)"raynum:",i,"|pha_id:",ray_pha(i),"|obs_tt:",abs_tt(i),'|cal_tt:',tt,"|dif_tt:",(abs_tt(i)-tt)
 270   format(a6,2x,i3,2x,a8,i3,2x,a8,f10.4,2x,a8,f10.4,2x,a8,f10.4)
!       if (abs(abs_tt(i)-tt)>2.0) cycle ! skip any ray shift larger than 2.0 sec
!   add contrain for ray different between observation - calculation
       dtt=(abs_tt(i)-tt)
       if (abs(dtt).ge.phcut) cycle
       if (iskip==1) cycle
       irow=irow+1
       ray_idx(irow)=i  ! ray used in inversion
 
       !-built up d matrix 
       d(irow)=(abs_tt(i)-dt)-tt
       res(irow)=d(irow)
!      print*,irow,res(irow),ray_pha(i)
!      print*,abs_tt(i),dt,tt,d(irow)

       !-built up G matrix
       ! calculate the average ray parameter
       !pave=0.
!      do j=1,np-1
!         print*,w(3,j),w(2,j),w(1,j)
!         re=ro-w(1,j)
!         tkofag(j)=angle(w,j,j+1)
!         call vel_3d(w(3,j),w(2,j),w(1,j),v0)
!         p=(re*sin(tkofag(j)*dpi))/v0
!         pave=pave+p
!      enddo
!      pave=pave/(np-1)
!      p=pave
       colat=90.0-ey
       re=ro-ez
       tkofag(irow)=angle(w,1,2)
       call distaz(90.-ey,ex,90.-sy,sx,epid(irow),azi(irow))   ! get epicenter distance & azimuth
!       write(*,'(a23)')"w(3,1),w(2,1),w(1,1),v0"
!       write(*,'(4f10.4)')w(3,1),w(2,1),w(1,1),v0
       call vel_3d(w(3,1),w(2,1),w(1,1),v0)
       epid(irow)=epid(irow)*dpi*ro
       p=(re*sin(tkofag(irow)*dpi))/v0
!      print*,"epid:",nint(epid(irow)),"azi:",nint(azi(irow)),"tkofag:",nint(tkofag(irow)),"v0:",v0,"p:",p
       G(irow,1)=(-1.)*(p*sin(azi(irow)*dpi)*sin(colat*dpi))*dpi
       G(irow,2)=(p*cos(azi(irow)*dpi))*dpi
       G(irow,3)=(-1)*cos(tkofag(irow)*dpi)/v0
       G(irow,4)=1.0
!      print*,sta_nm(sta_idx(i)),G(irow,1:4)
        
       !-calculate each coloumn norm
       col_norm(1)=col_norm(1)+G(irow,1)
       col_norm(2)=col_norm(2)+G(irow,2)
       col_norm(3)=col_norm(3)+G(irow,3)
       col_norm(4)=col_norm(4)+G(irow,4)

       !-estimate weighting matrix Wd
       if (int(mod(ray_pha(i),2))==1) Wd(irow,irow)=weighting(ray_wei(i),epid(irow),res(irow),tkofag(irow))
       if (int(mod(ray_pha(i),2))==2) Wd(irow,irow)=weighting(ray_wei(i),epid(irow),res(irow),tkofag(irow))*Ws2p  ! weight relatively poor picking of S
!      print*,sta_nm(i),ray_wei(i),epid(irow),abs_tt(i),tt,res(irow)
    enddo
    nrow=irow
    write(*,'(a26,i10)'),'+ Rays used for location =',nrow
    if (nrow<3) then
!1       print*,"Not enough rays to relocate!!"
       iflag=3
       return
    endif

    !-derive the rms
    do i=1,nrow
       rms=rms+res(i)**2
    enddo
    rms=sqrt(rms)/nrow
    if (icor==1.and.itr==1.and.irs==0) rms_ini=rms
!    print*,"rayused:",nrow,
    write(*,'(a7,f15.10)')"+  rms:",rms 
 
    !-output check file
    write(1,"(i4,i6,4g15.5,3g15.5)")itr,nrow,rms,adj,dt,ex,geoc_to_geog(ey),ez
 
    !-stop criteria
    if ((adj.lt.adjcut).or.(rms.lt.rmscut).or.(itr.eq.maxitr)) then
       ey=geoc_to_geog(ey)
       if (rms.le.rms_ini) iflag=1  ! well relocated
       if (rms.gt.rms_ini) iflag=4  ! illed relocated
!   Set relocated EQK has new depth greater than moho depth also illed relocation
!   Date modified: 2019/11/25 19:50 - 
       if(ez.gt.dm) iflag=4
       return
    endif
!   do i=1,nrow
!      Wd(i,i)=1.0
!      print*,Wd(i,1:nrow)
!   enddo
    m=0. ! 2019/11/13 why here m=0.0 

    call LS_inv(G,m,d,Wd,damp,nrow,4,GG)

    !-normalize the coloumn of matrix G
!   col_norm(1:4)=sqrt(abs(col_norm(1:4)))
!   do i=1,nrow
!      do j=1,4
!         G(i,j)=G(i,j)/col_norm(j)
!      enddo
!   enddo

    !-forming vectors for lsqr inversion
!   nonnull=0
!   do i=1,nrow
!      na(i)=4
!      do j=1,4
!         if (abs(G(i,j)).gt.1e-10) then
!            nonnull=nonnull+1
!            ra(nonnull)=G(i,j)
!            ja(nonnull)=j
!         endif
!      enddo
!   enddo
!   print*,nonnull,ra(1:nonnull)

    ncln=4   ! number of variables
    wantse=.true.
!   damp=0   ! given in "para.inp"
    atol=1.0e-4
    btol=1.0e-4
    conlim=1.0e7
    itnlim=4*ncln
    nout=46
    leniw=1
    lenrw=1
!    open(nout,file="lsqr.log")

!   v is not velocity vector but used in lsqr subroutine, x is the output.
!   call lsqr(nrow,ncln,damp,wantse,leniw,lenrw,iw,rw,d,v,ww,x,se,atol,btol,&
!             conlim,itnlim,nout,istop,itn,anorm,acond,rnorm,arnorm,xnorm)

    !-recover the solution vector
!   do i=1,4
!      x(i)=x(i)/col_norm(i)
!   enddo
    
!   print*,m
!   print*,x(1:4)
    !-obtain the dx,dy,dz,dt and update the new hypo coordinates
!    write(*,'(a48,4f10.4)')"Input location(long lat depth) and time shift:",ex,ey,ez,dt
    ex=ex+m(1)
    ey=ey-m(2)  ! "-" because the m(2) is the colatitude
    ez=ez+m(3)
    dt=dt+m(4)
!    write(*,'(a48,4f10.4)')"updated location(long lat depth) and time shift:",ex,ey,ez,dt
   !print*,ex,ey,ez,dt,m(4)

    ey=geoc_to_geog(ey)

    !-derive the adjustments of location
    call cal_delta(ey,ex-0.5,ey,ex+0.5,dx)
    call cal_delta(ey-0.5,ex,ey+0.5,ex,dy)
    adj=adj+(x(1)*dy)**2
    adj=adj+(x(2)*dx)**2
    adj=adj+x(3)**2
    adj=sqrt(adj)

    !-check if the relocated hypo is out of range
    if ((ex-lon_c(2))*(ex-lon_c(nlon_c-1))>0..or.(ey-lat_c(2))*(ey-lat_c(nlat_c-1))>0..or.(ez-0.)*&
        (ez-dep_c(ndep_c-1))>0.) then
        iflag=2
       return
    endif

    !-view results
!   write(*,"(3(a12,g20.10))")"rms:",rms,"adjustment:",adj,"tshift:",m(4)
!     write(*,*)"Update location:"
!     write(*,"(a4,f7.3,2(a6,f7.3))")"x: ",ex," - y: ",ey," - z: ",ez
end subroutine
