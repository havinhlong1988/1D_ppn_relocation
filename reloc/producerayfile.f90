       program producerayfile

        character*80 line, infile
        character*2 ph
        character*1 char
        integer iray,np,iph
        real*8 lat,long,dep

        read(*,*) infile

        open(2,file="pnxray.plot",status='unknown')
        open(3,file="pnyray.plot",status='unknown')
        open(4,file="pgxray.plot",status='unknown')
        open(5,file="pgyray.plot",status='unknown')
        open(6,file="pghray.plot",status='unknown')
        open(7,file="pnhray.plot",status='unknown')
!
        open(1,file=infile,status='old')
        l=0
        do
            read(1,'(a)',iostat=ios) line
            if(ios.lt.0) exit
            l=l+1
            if(line(1:1).eq.">") then
                read(line,'(a1,2i5,4x,a2)')char,iray,np,ph
                if(ph=='Pn') then
                    write(2,'(a1)') char
                    write(3,'(a1)') char
                    write(7,'(a1)') char
                elseif(ph=='Pg') then
                    write(4,'(a1)') char
                    write(5,'(a1)') char
                    write(6,'(a1)') char                    
                endif
            elseif(line(1:1).eq." ") then
                read(line,'(3f10.4,3x,i1)')lat,long,dep,iph
                if(iph==3) then
                    write(2,'(2f10.4)')long,dep
                    write(3,'(2f10.4)')dep,lat
                    write(7,'(2f10.4)')long,lat
                elseif(iph==1) then
                    write(4,'(2f10.4)')long,dep
                    write(5,'(2f10.4)')dep,lat
                    write(6,'(2f10.4)')long,lat
                endif
            else
                print*,"Wrong on raypath file"
            endif
                
        enddo

        do i=1,7
            close(i)
        enddo
            
       end program