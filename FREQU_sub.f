      SUBROUTINE FREQU                                                  
      IMPLICIT REAL*8(A-H,O-Z)
      DOUBLE PRECISION WUV,WCD
      COMMON idummy,NOP,NOS,IB,U,UL,UU,DU,IRUN,IPUNCH                          
      COMMON /AIR/  AI(100,1200),AR(100,1200)                               
      COMMON /POL/  AIT(100),ART(100),BM(100),UA(100)                       
      DIMENSION WUV(500),WCD(500)
      DIMENSION WORD(500)
      DIMENSION UUU(500)
C     INSERITA PER PLOT
      JJJ=0
      NI=1.0+(UU-UL)/DU +1.0E-10                                        
      DO 1 J=1,NI                                                       
      ZZ=J-1                                                            
      U=UL+ZZ*DU                                                        
      UUU(J)=U
C     INSERITA PER PLOT
      DO 2 I=1,100                                                       
      IF(UA(I).EQ.0.)  GO TO 2                                          
      K=1.0+(U-UA(I))/DU+1.0E-10                                        
      AIT(I)=-AI(I,K)                                                   
      ART(I)=AR(I,K)                                                    
    2 CONTINUE                                                          
  200 FORMAT (1H+,F7.3)                                                 
      CALL INTRAC                                                       
      CALL AMAT                                                         
      IF (IRUN.NE.0)  GO TO 3                                           
      WRITE (6,201)                                                     
      WRITE (6,202)                                                     
    3 CONTINUE                                                          
  201 FORMAT(///,1X,12HOPTICAL DATA  )                                  
  202 FORMAT(/,3X,5HFREQZ,9X,5HTRACE,11X,15HAXIS 1, X, Y,Z ,17X,15HAXIS 
     12, X, Y, Z ,17X,15HAXIS 3, X, Y, Z  ,/)                           
  199 FORMAT (1X)                                                       
C     WRITE (6,199)                                                     
C     WRITE (6,200) U                                                   
      JJJ=JJJ+1
      CALL UVSHPE(JJJ,WUV)
C     INSERITA PER PLOT
      CALL RTSHPE(JJJ,WCD,WORD)
C     INSERITA PER PLOT
  198 FORMAT(1X,I5,2F15.6)
      IRUN=IRUN+1                                                       
    1 CONTINUE                                                          
      WRITE(6,203)
      WRITE(6,204)(UUU(I),WUV(I),WCD(I),WORD(I),I=1,NI)
      write(10,1000)(10000./UUU(I),WUV(I),I=1,NI)
      write(11,1000)(10000./UUU(I),WCD(I),I=1,NI)
 1000 format(f10.2,5x,f10.2)
  203 FORMAT(5X,'FREQZ',7X,'UV',10X,'CD',3X)
  204 FORMAT(3(F5.1,2X,F9.1,2X,F8.2,2X,F9.2,4X))
      UUUNI=UUU(NI)
      CALL PRNPLT(UUU,WUV,UUUNI,DU,UVMAX,UVINCR,1,1,NI)
      CALL PRNPLT(UUU,WCD,UUUNI,DU,CDMAX,CDINCR,1,1,NI)
C     INSERITA PER PLOT
      RETURN                                                            
      END                   
