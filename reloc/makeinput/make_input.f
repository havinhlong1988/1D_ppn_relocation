c make input file for 3d inversion from 1d inversion result
	program make_input
	character*8 st
	character*8 stnm(1000)
	integer evid(100000),id
	integer evdate,evtime
	real evla,evlo,evdp,evmg,stla,stlo,stel,ttime
	integer phase
	character*20 ray_fn

	print*,"input the ray-info file name:"
	read(*,*)ray_fn
	nevt=0
	nsta=0
	mark=0
	open(10,file=ray_fn)
	open(11,file="sta_inp")
	open(12,file="evt_inp")
	open(13,file="tt_inp")
	evdate=0
	evtime=0
	evmg=0.0

5	read(10,*,end=888)st,junk1,evla,evlo,evdp,stla,stlo,
     &		          stel,ttime,trash,phase,id
	if (id.ne.mark) then		!new observation
	   write(13,"(a1,1x,i10)")"#",id
	   mark=id
	endif
	ista=0
	do i=1,nsta
	   if (stnm(i).eq.st) then 	!old station
	      ista=i
	      exit
	   endif
	enddo 
	if (ista.eq.0) then  !new station
	   nsta=nsta+1
	   stnm(nsta)=st
	   write(11,100)st,stla,stlo,-stel*1000
	endif
	ievt=0
	do i=1,nevt
	   if (evid(i).eq.id) then	!old event
	      ievt=i
	      exit
	   endif
	enddo
	if (ievt.eq.0) then  !new event
	   nevt=nevt+1
	   evid(nevt)=id
	   write(12,200)evdate,evtime,evla,evlo,evdp,evmg,0.00,0.00,0.00,id
	endif
	if (phase.eq.1) write(13,300)st,ttime,1.0,"Pg"
	if (phase.eq.2) write(13,300)st,ttime,1.0,"Pn"
	goto 5
888	close(10)
	close(11)
	close(12)
	close(13)
100	format(a8,1x,f7.4,1x,f8.4,1x,f6.1)
200	format(i10,1x,i10,1x,f6.2,1x,f7.2,1x,f7.2,1x,4(f3.1,1x),i10)
300	format(a5,1x,f7.3,1x,f3.1,1x,a2)
	end
