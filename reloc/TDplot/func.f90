real*8 function geog_to_geoc(xla)
    implicit none
    real*8 xla,RAD_PER_DEG,B2A_SQ
    RAD_PER_DEG=0.0174532925199432955
    B2A_SQ=0.993305521
    geog_to_geoc = atan(B2A_SQ*tan(RAD_PER_DEG*xla)) / RAD_PER_DEG
    return
end function geog_to_geoc


real*8 function geoc_to_geog(xla)
    implicit none
    real*8 xla,RAD_PER_DEG,B2A_SQ
    RAD_PER_DEG=0.0174532925199432955
    B2A_SQ=0.993305521
    geoc_to_geog = atan(tan(RAD_PER_DEG*xla)/B2A_SQ) / RAD_PER_DEG
    return
end function geoc_to_geog



FUNCTION JULDAY(MM,ID,IYYY)                                  
    PARAMETER (IGREG=15+31*(10+12*1582))                         
    IF (IYYY.EQ.0) STOP 'There is no Year Zero.'                
    IF (IYYY.LT.0) IYYY=IYYY+1                                   
    IF (MM.GT.2) THEN                                            
       JY=IYYY                                                    
       JM=MM+1                                                    
    ELSE                                                         
       JY=IYYY-1                                                  
       JM=MM+13                                                   
    ENDIF                                                        
    JULDAY=INT(365.25*JY)+INT(30.6001*JM)+ID+1720995             
    IF (ID+31*(MM+12*IYYY).GE.IGREG) THEN                        
       JA=INT(0.01*JY)                                            
       JULDAY=JULDAY+2-JA+INT(0.25*JA)                            
    ENDIF                                                        
    RETURN                                                       
END function JULDAY                                                                 
