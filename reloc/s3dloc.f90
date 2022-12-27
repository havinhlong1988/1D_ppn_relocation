program spherical_3dloc
    implicit integer (i-n)
    include "para.inp"
    character*20 sta_fn,evt_fn,abs_fn,vel_fn,evtoutf,outfna,outfnb
    character*3 label
    integer sta_count(maxsta),sta_azi(maxevt,maxsta),sta_gap
    integer iunit,idp,unit,nrms_sum(maxitr)
    real*8 x0,y0,z0,dt,rms_sum(maxitr),rmsout
    real*8 m(4),dx,dy,msum(4)
    


    call random_seed()  

    print*,"input event data file name:"
    read*,evt_fn
    print*,trim(evt_fn)
    print*,"input station data file name:"
    read*,sta_fn
    print*,trim(sta_fn)
    print*,"input absolute travel time file name:"
    read*,abs_fn
    print*,trim(abs_fn)
    print*,"input velocity model file name:"
    read*,vel_fn
    print*,trim(vel_fn)
    print*,"input the moho depth:"
    read*,dm
    print*,dm

    call getdata1(evt_fn,sta_fn,abs_fn)  ! Separate files of evt, sta, and tt in Haijiang format
!   call getdata2(evt_fn,sta_fn)         ! Phase file in CWB format)
    call loadvel(vel_fn)   
!
!      print*,"imod=", imod
!      pause
!
    open(1,file="s3dloc.chk",status="unknown")
    open(99,file="run.log",status="unknown")
    label="xxx"
    unit=20
    do i=1,maxitr
      write(label(1:3),'(i3.3)') i
      iunit=unit+i
      open(iunit,file="ray_path_"//label,status='unknown') ! Raypath write out
      open(iunit+10,file="residual_"//label,status='unknown') ! residual write out
!
      nrms_sum(i)=0
      rms_sum(i)=0.0
!      print*,iunit
    enddo
    !-residual check file annotations
    write(1,*)"/rms, adjustment, tshift, x0, y0, z0/"
    print*,"nevt:",nevt,"nsta:",nsta,"nray:",nray
    write(99,*)"nevt:",nevt,"nsta:",nsta,"nray:",nray
    !-loops of station corrections
    icor=0
    do 
       icor=icor+1
       write(0,"(a)")"*********************************"
       write(0,"(a50,2x,i3)")"Run the reloction with station correction number >",icor
       write(99,"(a)")"*********************************"
       write(99,'(a8,2x,i3)')"Icor ===>",icor
       write(1,"(a8,i4,a9)")"========",icor," ========"

       open(2,file="s3dloc.err",status="unknown")
       open(3,file="s3dloc.lis",status="unknown")
       open(4,file="s3dloc.out",status="unknown")
       open(8,file="00_evt_cut.out",status="unknown")
       open(10,file="s3dloc.tt",status="unknown")
       open(11,file="s3dloc.rms",status="unknown")
       open(55,file="err_check",status="unknown")
       !-travel-time calibration by station corrections
       do i=1,nray
          j=int(mod(ray_pha(i),2))
          abs_tt(i)=abs_tt(i)-sta_cor(sta_idx(i),j)
       enddo
!       print*,sta_nm(sta_idx(1)),abs_tt(1) 
       sta_count=0
       do ievt=1,nevt
          write(0,'(a)')"---------------------------------------------------"
          write(*,'(a10,i6,a10,i6)')'# event: ',ievt," #id: ",evt_id(ievt)
          write(99,'(a)')"---------------------------------------------------"
          write(99,'(a10,i6,a10,i6)')'# event: ',ievt," #id: ",evt_id(ievt)
          evtnow=ievt
          nrds=evt_pack(ievt+1)-evt_pack(ievt)
          write(1,*)"#",evt_id(evtnow),nrds
!         print*,nrds,evt_pack(ievt)+1,evt_pack(ievt+1)
          call ini_hypo(x0,y0,z0)  
          write(*,"(a15,3f7.3)")"Input location: ",x0,y0,z0
          write(99,"(a15,3f7.3)")"Input location: ",x0,y0,z0
          !-loops of relocation iterations
          itr=0
          irs=0
          dt=0.0
          epid=0.0
          idp=0
          do i=1,4
            msum(i)=0.0
          enddo
!          write(55,'(4f10.4)')msum
          write(55,'(a54)')"Model pertubation:/ long, colat, dep, origintime shift"
          write(55,'(a8,i6)')"evt ->: ",evtnow
          do
             itr=itr+1
             idp=unit+itr
             print*,"--> itr: ",itr
             write(99,*)"--> itr: ",itr
!             print*,"#-> idp: ",idp
             iflag=0
             nrow=0
             write(55,'(a6,2x,i6)')"> itr:",itr
             call eqkloc(x0,y0,z0,dt,m,idp,rmsout)
             do iii=1,4
               msum(iii)=msum(iii)+m(iii)
             enddo
             write(55,'(4f10.4)')m
             write(55,'(4f10.4)')msum
             !
             write(99,"(a52,2x,4f10.3)")"Model variation/ long, colat,dep,origin time shift: ",m 
             write(99,"(a58,2x,4f10.3)")"Total model variation/ long, colat,dep,origin time shift: ",msum
             write(*,"(a21,3f10.3)")"Update the location: ",x0,y0,z0
             write(99,"(a21,3f10.3)")"Update the location: ",x0,y0,z0             
             ! rms for all event
             rms_sum(itr)=rms_sum(itr)+rmsout*rmsout
             nrms_sum(itr)=nrms_sum(itr)+1

             !-condition notes
             if (iflag==1) then
                exit
             else if (iflag==2) then          
                write(1,*)"Out-of-range convergence..",x0,y0,z0
                !-perturb earthquake location to search another local minimum as being out of range
                irs=irs+1
                if(irs>maxirs) exit     ! round limit of researching
                call lottery(x0,y0,z0)  ! applied for small damping
                print*,"* researching minimum..."
                write(99,*)"* researching minimum..."
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
          if (iflag.ne.1) then
             !-illed-condition events, ouput they errors
             write(2,50)evt_date(evtnow),evt_time(evtnow),evt_loc(evtnow,1),evt_loc(evtnow,2),evt_loc(evtnow,3),evt_mag(evtnow),&
                        0.0,0.0,0.0,evt_id(evtnow),y0,x0,z0,iflag
          else
             !-output the final report of well-condition event
            if((iluflag.eq.0).or.(iluflag.eq.1)) then
               if(iluflag.eq.0) then
                  call output1(x0,y0,z0,dt,m)  ! Separate files of evt, sta, and tt in Haijiang format
               else
                  call output1(x0,y0,z0,dt,msum)  ! Separate files of evt, sta, and tt in Haijiang format
               endif
            else
               print*,"error!! wrong setting of the location uncertainty flag <iluflag> on para.inp"
               print*,"iluflag=",iluflag
               stop
            endif
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
!     Write out RMS for each iteration
       write(*,'(a34)')"Iteration number - nresidual - RMS"
       write(99,'(a34)')"Iteration number - nresidual - RMS"       
       do i=1,maxitr
         rms_sum(i)=sqrt(rms_sum(i)/float(nrms_sum(i)))
         write(*,'(2(i6,1x),f10.3)')i,nrms_sum(i),rms_sum(i)
         write(99,'(2(i6,1x),f10.3)')i,nrms_sum(i),rms_sum(i)
         write(11,'(2(i6,1x),f10.3)')i,nrms_sum(i),rms_sum(i)
       enddo
50     format(2(i9,1x),2(f8.3,1x),f7.2,1x,4(f4.1,1x),i6,1x,2(f8.3,1x),f7.2,1x,i6)
       close(2)
       close(3)
       close(4)
       !-obtain station corrections
       write(*,'(a)')"Station number - nray used to correction - correction value"
       write(99,'(a)')"Station number - nray used to correction - correction value"
       do i=1,nsta
          if (sta_count(i)<10) cycle
          call get_gap(sta_azi(1:sta_count(i),i),sta_count(i),sta_gap)
!         print*,sta_count(i),"gap:",sta_gap,"azi:",sta_azi(1:sta_count(i),i)
          if (sta_gap>=90.) cycle
          j=int(mod(ray_pha(i),2))
          sta_cor(i,j)=sta_cor(i,j)/sta_count(i)
          write(*,'(i8,2x,i18,f15.5)')i,sta_count(i),sta_cor(i,j)
          write(99,'(i8,2x,i18,f15.5)')i,sta_count(i),sta_cor(i,j)
       enddo
       if (icor==maxcor) then  ! program terminates
          open(5,file="s3dloc.cor",status="unknown")
          do i=1,nsta
             write(5,"(a5,2x,5(f8.3,1x))")sta_nm(i),(sta_loc(i,j),j=1,3),(sta_cor(i,k),k=1,2)
          enddo
          close(5)
          exit
       endif
    enddo
    close(1)
    close(10)
    close(11)
    close(99)
    close(8)
    close(55)
    do i=1,maxitr
      close(20+i)
    enddo
    call reloc2inv1d
    outfna=trim("00_evt_cut.out")
    outfnb=trim("output_data_cut")
!  Cut data
    call outdata(outfna,sta_fn,abs_fn,outfnb)
!  no cut data
    outfna=trim("00_evt.out")
    outfnb=trim("output_data")
    call outdata(outfna,sta_fn,abs_fn,outfnb)
!  plot compare the cut data phase and origion
    call system('./00.compare_phases.sh')
    call system('wc -l 00_evt.out')
end program

subroutine eqkloc(ex,ey,ez,dt,m,index,rmsout)
    implicit integer (i-n)
    include "para.inp"
    character*2 pha
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
    real*8 se(maxvar)
    real*8 rw(1),rmsout
    integer iw(1),index,ie,is
    common /matrix/ ra,ja,na
    dpi=asin(1.)/90.

    !-initialize matrixes
    G=0.0
    d=0.0
    Wd=0.0
    col_norm=0.0
    adj=adjcut+1e-10

!   print*,dt
    irow=0
    ey=geog_to_geoc(ey)
    do i=evt_pack(evtnow)+1,evt_pack(evtnow+1)
      if(ijoint==0) then
         if(imod==1) then
           if (ray_pha(i).ne.1) cycle      ! only use crutal P
         elseif(imod==2) then
           if ((ray_pha(i).ne.1).or.(ray_pha(i).ne.2)) cycle ! only use crutal P & S
         endif
      endif
!      if ((abs_tt(i)-dt)<=0.0) cycle  ! wrong time shifting, skip this ray
       ivmod=int(mod(ray_pha(i),2))    ! distinguish which velocity model should be used
       sy=sta_loc(sta_idx(i),1)
       sx=sta_loc(sta_idx(i),2)
       sz=sta_loc(sta_idx(i),3)
       sy=geog_to_geoc(sy)
!      print*,"evt:",ex,ey,ez,sta_nm(sta_idx(i)),sx,sy,sz
       iskip=0
!       dm=32.321
       if (ray_pha(i)==1.or.ray_pha(i)==2) then  ! Pg & Sg
          call pbr(ey,ex,ez,sy,sx,sz,w,np,tt)
          do j=1,np
!            call depth(w(3,j),w(2,j),dm)
             if (dm<=w(1,j)) then
                iskip=1     ! skip the ray pass down to moho
                exit
             endif
             if(dm<=ez) then
               iskip=1      ! skip the ray as hypo below moho
               exit
             endif
          enddo
       else if (ray_pha(i)==3.or.ray_pha(i)==4) then  ! Pn & Sn
!         call depth(ex,ey,dm)
          if (dm<=ez) then
             iskip=1        ! skip the ray as hypo below moho
             exit
          endif
          call pn_path(ex,ey,ez,sx,sy,sz,w,np,tt,ie,is)
       else
          stop "Phase out of the ray-tracing list!!"
       endif
!     Write out the calculated ray path (ray tracing)
!       print*,"index",index
       if(iskip.ne.1) then
         if((ray_pha(i)==1).or.(ray_pha(i)==2)) then
            if((ray_pha(i)==1)) then
               write(index,'(a1,2i5,4x,a)')">",i,np,"Pg"
            elseif((ray_pha(i)==2)) then
               write(index,'(a1,2i5,4x,a)')">",i,np,"Sg"
            endif
!            
            do ii=1,np
               write(index,'(3f10.4,3x,i1)')(w(j,ii),j=2,3),w(1,ii),ray_pha(i)
            enddo
         elseif((ray_pha(i)==3).or.(ray_pha(i)==4)) then
            if((ray_pha(i)==3)) then
               write(index,'(a1,2i5,4x,a)')">",i,np,"Pn"
            elseif((ray_pha(i)==4)) then
               write(index,'(a1,2i5,4x,a)')">",i,np,"Sn"
            endif
            do ii=1,ie
               write(index,'(3f10.4,3x,i1)')(w(j,ii),j=2,3),w(1,ii),ray_pha(i)
            enddo
            do ii=ie+is+1,np
               write(index,'(3f10.4,3x,i1)')(w(j,ii),j=2,3),w(1,ii),ray_pha(i)
            enddo
            do ii=ie,ie+is-1
               write(index,'(3f10.4,3x,i1)')(w(j,(2*ie+is-ii)),j=2,3),w(1,(2*ie+is-ii)),ray_pha(i)
            enddo
         endif
       endif
!
       if (iskip==1) cycle
       if((abs_tt(i)-tt).gt.cut_off) cycle
       irow=irow+1
       ray_idx(irow)=i  ! ray used in inversion
 
       !-built up d matrix 
       d(irow)=(abs_tt(i)-dt)-tt
       res(irow)=d(irow)
!      print*,irow,res(irow),ray_pha(i)
!      print*,abs_tt(i),dt,tt,d(irow)

       !-built up G matrix
       ! calculate the average ray parameter
       pave=0.
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
       if (int(mod(ray_pha(i),2))==0) Wd(irow,irow)=weighting(ray_wei(i),epid(irow),res(irow),tkofag(irow))*Ws2p  ! weight relatively poor picking of S
   !   print*,irow,sta_nm(irow),ray_wei(irow),epid(irow),abs_tt(i),tt,res(irow),tkofag(irow),Wd(irow,irow)
      !  print*,irow,ray_pha(irow),ray_wei(irow),epid(irow),res(irow),tkofag(irow),Wd(irow,irow)
    enddo
    nrow=irow
    if (nrow<3) then
       print*,"Not enough rays to relocate!!"
       iflag=3
       return
    endif

    !-derive the rms
    do i=1,nrow
       rms=rms+res(i)**2
       write(index+10,'(i6,2x,f7.3)')i,res(i)
    enddo
    rms=sqrt(rms)/nrow
    if (icor==1.and.itr==1.and.irs==0) rms_ini=rms
    write(99,'(a11,i6,a6,f10.8)')"nray used: ",nrow," rms: ",rms 
    write(*,'(a11,i6,a6,f10.8)')"nray used: ",nrow," rms: ",rms 
    rmsout=rms
    !-output check file
    write(1,"(i4,i6,4g15.5,3g15.5)")itr,nrow,rms,adj,dt,ex,geoc_to_geog(ey),ez
   
    !-stop criteria
    if ((adj.lt.adjcut).or.(rms.lt.rmscut).or.(itr.eq.maxitr)) then
       ey=geoc_to_geog(ey)
       if (rms.le.rms_ini) iflag=1  ! well relocated
       if (rms.gt.rms_ini) iflag=4  ! illed relocated
       return
    endif
!   do i=1,nrow
!      Wd(i,i)=1.0
!      print*,Wd(i,1:nrow)
!   enddo
    m=0.0

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
    open(nout,file="lsqr.log")

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
    ex=ex+m(1)
    ey=ey-m(2)  ! "-" because the m(2) is the colatitude
    ez=ez+m(3) 
    dt=dt+m(4)
!  2022/12/26 Long - tried to force the negative depth to 0
!  Of course it make artificial error in depth, but can not accepth by RMS 
    if ((iforcedepth==1).and.(ez<0)) then
        ez=0
    endif
!   print*,ex,ey,ez,dt

    ey=geoc_to_geog(ey)

    !-derive the adjustments of location
!    call cal_delta(ey,ex-0.5,ey,ex+0.5,dx)
!    call cal_delta(ey-0.5,ex,ey+0.5,ex,dy)
!    adj=adj+(x(1)*dy)**2
!    adj=adj+(x(2)*dx)**2
!    adj=adj+x(3)**2
!    adj=sqrt(adj)
!   print*,"x(1),x(2),x(3)",x(1),x(2),x(3)
    !-check if the relocated hypo is out of range
    if ((ex-lon_c(2))*(ex-lon_c(nlon_c-1))>0..or.(ey-lat_c(2))*(ey-lat_c(nlat_c-1))>0..or.(ez-0.)*&
        (ez-dep_c(ndep_c-1))>0.) then
         if(iluck.eq.1) then
            iflag=2   ! Turn on the option lottery
         elseif(iluck.eq.0) then
            iflag=4  ! Turn off the option lottery
         else
            write(*,*)"Warning! Wrong set on iluck flag! Check the para.inp file"
         endif
       return
    endif
!    if(ez.gt.dm) iflag=4

    !-view results
!   write(*,"(3(a12,g20.10))")"rms:",rms,"adjustment:",adj,"tshift:",m(4)
!   write(*,"(4(a4,g20.10))")"x:",ex,"y:",ey,"z:",ez
end subroutine
