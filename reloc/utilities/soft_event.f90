      program soft_event
!   Program reading the output file of reloctation of ppn inversion
!   The sort-out event based location uncertainty criteria (Husen and Hardebeck 2010)
!   Input: + s3dloc.lis file
!   Output: + 00_evt.out file
!           + 00_evt1.out file
        integer date,time,id
        real*8 lat,long,dep,mag,rms,erh,erz

        open(1,file='s3dloc.lis',status='old')
        open(2,file="00_evt.out",status="unknown")
        open(3,file="00_evt1.out",status="unknown")

        do
          read(1,*,iostat=ios)date,time,lat,long,dep,mag,rms,erh,erz,id
          if(ios.lt.0) exit
!   Apply the sort-out criteria
          if((erh.lt.10).and.(erz.lt.25)) then
            write(2,111)date,time,lat,long,dep,mag,rms,erh,erz,id
          endif
          if((erh.lt.20).and.(erz.lt.50)) then
            write(3,111)date,time,lat,long,dep,mag,rms,erh,erz,id
          endif
        enddo

        close(1)
        close(2)
        close(3)
 111   format(2(2x,i8.8),2(f7.3,1x),f7.2,1x,4(f5.1,1x),i6)
      end



      subroutine outdata(evt_fn1,sta_fn1,abs_fn1,outfn1)
        implicit integer (i-n)
        character*20 evt_fn1,sta_fn1,abs_fn1,outfn1
        character*90 line
        character*1 str
        character*2 phcard
        character*5 stnm
        real*8 junk1,junk2,junk3,evt_lat,evt_lon,evt_dep,sta_lat,sta_lon,sta_elv
        real*8 ttime,wei,evtla,stala,dist,az
        real*8 geog_to_geoc, geoc_to_geog,dist1(maxrds),az1(maxrds)
     
     !       outfn1=trim(outfn1)
        print*,"produce T-D data"
        open(1,file=evt_fn1,status="old")
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
     
        open(2,file=sta_fn1,status="old")
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
     
        open(3,file=abs_fn1,status="old")
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
             sta_elv=sta_loc(sta_idx(i),3)
     
             evtla=geog_to_geoc(evt_lat)
             stala=geog_to_geoc(sta_lat)
     
             call distaz(90.0-evtla,evt_lon,90.0-stala,sta_lon,dist,az) ! estimate the distance and azimuth angle
             dist1(i)=dist
             az1(i)=az
          enddo
     !   
        open(4,file=outfn1,status="unknown")
        do i=1,nray
           write(4,10)sta_nm(sta_idx(i)),evt_id(evt_idx(i)),evt_date(evt_idx(i)),evt_time(evt_idx(i)),&
                      (evt_loc(evt_idx(i),j),j=1,3),evt_mag(evt_idx(i)),(sta_loc(sta_idx(i),j),j=1,3),&
                      abs_tt(i),dist1(i)*111,az1(i),i,ray_pha(i)
        enddo
        close(4)
 10  format(a,1x,i6,1x,2(i8.8,1x),10(f8.3,1x),2i6)
     endsubroutine