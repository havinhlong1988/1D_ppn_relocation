       subroutine reloc2inv1d
        character*99 line
        character*4 sta_nm
        integer iday,itime,id
        real elat,elong,edp,emag,junk1,junk2,junk3
        real elato,elongo,edpo,sta_lat,sta_long,sta_el,cor1,cor2



!   event file after reloc and new evt for 1d inversion
            open(1,file='s3dloc.lis',status='old')
            open(2,file='00_evt.out',status='unknown')
!   station file after reloc and new station file for 1d inversion
            open(3,file='s3dloc.cor',status='old')
            open(4,file='00_sta.out',status='unknown')

!   Read event data
            print*,""
            print*,"process!! read the relocated event file: <- s3dloc.lis"
            print*,""
        do
            read(1,*,iostat=ios)iday,itime,elat,elong,edp,emag&
            ,junk1,junk2,junk3,id,elato,elongo,edpo
            if(ios.lt.0) exit
            
            write(2,101)iday,itime,elat,elong,edp,emag&
            ,junk1,junk2,junk3,id
        enddo

!   read station updated information
        print*,'process!!!,reading updated station file: s3dloc.cor'
        print*,""

        do
            read(3,*,iostat=irr)sta_nm,sta_lat,sta_long,sta_el,cor1,cor2
            if(irr.lt.0) exit
!            if((cor1.lt.2).and.(cor1.gt.-2)) then
!                if((cor2.lt.2).and.(cor2.gt.-2)) then
                 write(4,200)sta_nm,sta_lat,sta_long,&
                 -(sta_el*1000),cor1,cor2
!                endif
!            endif
        enddo
!        print*,"warning: station correction -2< and >2 have been remove"
        print*,""
        print*,"write to --> output file: -> 00_evt.out"
!
        print*,""
        print*,"write to output file: -> 00_sta.out"
        print*,""


 100    format(2(i9,1x),2(f8.3,1x),f7.2,1x,4(f4.1,1x),&
               i6,1x,2(f8.3,1x),f7.2)

 101    format(2(2x,i8),2(f7.3,1x),f7.2,1x,4(f4.1,1x),i6)
 200    format(a4,4x,2(f8.4,1x),f6.1,1x,2(f7.3,1x))



      end subroutine
