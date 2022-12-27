       program tes_sincostan
        parameter(pi=3.14159265359)
        real*8 a,b,c,d,a1,b1


        a = 20.808
        b = 21.3338
        c = geog_to_geoc(a)
        d = geog_to_geoc(b)
        a1 = geoc_to_geog(c)
        b1 = geoc_to_geog(d)
        write(*,'(2f20.15)')a,b
        write(*,'(2f20.15)')c,d
        write(*,'(2f20.15)')a1,b1

       end program

       real*4 function geog_to_geoc(xla)
       implicit none
       real*8 xla,RAD_PER_DEG,B2A_SQ
       RAD_PER_DEG=0.0174532925199432955
       B2A_SQ=0.993305521
       geog_to_geoc = atan(B2A_SQ*tan(RAD_PER_DEG*xla)) / RAD_PER_DEG
       return
       end function geog_to_geoc

      real*4 function geoc_to_geog(xla)
       implicit none
       real*8 xla,RAD_PER_DEG,B2A_SQ
       RAD_PER_DEG=0.0174532925199432955
       B2A_SQ=0.993305521
       geoc_to_geog = atan(tan(RAD_PER_DEG*xla)/B2A_SQ) / RAD_PER_DEG
       return
      end function geoc_to_geog