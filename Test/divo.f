      IMPLICIT REAL*8(A-H,O-Z)
      COMMON idummy,NOP,NOS,IB,U,UL,UU,DU,IRUN,IPUNCH                          
      COMMON /IMU/ noint(1200,1200),nnoint
	dimension jj(1200)

      READ(5,*) nnoint,ifile
c ifile=0 legge nel formato standard di dikrat
c ifile=1 legge le coordinate atomiche in formato .mol
c ifile=2 legge le coordinate atomiche in formato .ent

      model=1

      do l=1,nnoint
      READ(5,*)n,(jj(k), k=1,n)
	DO K=1,N-1
	I=jj(K)
	DO M=k+1,n
	J=jj(M)
	noint(I,J)=1
	ENDDO
      ENDDO
      enddo

    5 CALL INPUT(ifile,model)
      IRUN=0                                                            
      CALL FREQU 
	model=model+1
      GO TO 5                                                           
      END

      SUBROUTINE INPUT(ifile,model)
      IMPLICIT REAL*8(A-H,O-Z)
      REAL*4 COSS,SENN,TETA,COSB
      REAL*4 COSA,SENA
      COMMON idummy,NOP,NOS,IB,U,UL,UU,DU,IRUN,IPUNCH                          
      COMMON /ROSC/ XR(1200),YR(1200),ZR(1200)                             
      COMMON /EEOS/ XE(1200),YE(1200),ZE(1200)                             
      COMMON /EMOS/ XM(1200),YM(1200),ZM(1200)                             
      COMMON /POL/  AIT(100),ART(100),BM(100),UA(100)                       
      COMMON /IDET/ IDA(1200),IDB(1200),IDE(1200),IW(4800)                  
      DIMENSION X(10000),Y(10000),Z(10000),iref(4,1200),f(1200)                                    
      character title*80,molfile*20,att(10000)*2,dum(80),
     &control*6 


      if(model.eq.1.or.ifile.eq.0) then 
	READ (5,100) title
	if(title.eq.'')stop
	else
	isep=0
	endif
  100 FORMAT (A78)  
   
	                                      
      ISUP=ISEP
      if(ifile.le.1.or.imodel.eq.1) WRITE (6,200)title                                       
  200 FORMAT(a78//'Coupled oscillator model to infinite order'
     &/'based on the program DIKRAT written by W. HUG'/
     &'and modified by M. Zandomeneghi and L. Di Bari'//)   

	if(ifile.eq.1) then
	write(6,1005)
 1005 format('By selecting the flag "1" in line 1
     & a ".mol" coordinate file was chosen' )
	read(5,1004)molfile
 1004 format(a20)
	write(6,1003)molfile
 1003 format(/'Coordinates taken from the file: ',a20//)
	open(unit=7,file=molfile,status='old')
      read(7,1000)title
 1000 format(a80//)
      read(7,1001)nat,nbond
 1001 format(2i3)
      do i=1,nat
      read(7,1002)x(i),y(i),z(i),att(i)
 1002 format(3f10.4,1x,a2)
      enddo
	go to 25
	endif

      if(ifile.eq.2) then
	if(model.eq.1) then
	write(6,1006)
 1006 format('By selecting the flag "2" in line 1 a ".ent" coordinate'
     &' file was chosen')
	read(5,1004)molfile
	write(6,1003)molfile
	open(unit=7,file=molfile,status='old')
	isep=0
	endif
	
	i=1
 1100	continue
      read(7,1007,ERR=1100,END=9999)control,xx,yy,zz

 1007 format(a6,24x,3f8.3)
      if(control.eq.'ATOM  '.or.control.eq.'HETATM')
     %then
	  x(i)=xx
	  y(i)=yy
	  z(i)=zz
	  i=i+1
	  endif
	if(control.eq.'END   ') go to 1101
	if(control.eq.'ENDMDL') go to 1101
	go to 1100

 1101 go to 25  
	endif


C      IF (IPUNCH.NE.0) WRITE (7,220) title                     
C 220 FORMAT(1X,A76)                                                 
C      IF (IPUNCH.NE.0) WRITE (6,221)                                    
C 221 FORMAT(//1X,5HPUNCH  )                                            
      WRITE (6,202)                                                     
  202 FORMAT(//1X,10HINPUT DATA  )                                      
      WRITE (6,203)                                                     
  203 FORMAT(//1X,16HREFERENCE POINTS  )                                
      WRITE (6,204)                                                     
  204 FORMAT(/9X,1HX,9X,1HY,9X,1HZ,/)   
      

      

	   
                                 
       DO 1 I=1,10000                                                      
         IF(ISUP.EQ.8) GO TO 2
         READ (5,101) X(I),Y(I),Z(I),ISEP                                  
  101    FORMAT (3E10.0,49X,I1)                                            
           IF (ISEP.NE.0)  GO TO 2                                           
         WRITE (6,205) I,X(I),Y(I),Z(I)                                    
    1  CONTINUE
    
                                                              
    2 CONTINUE                                                          
  205 FORMAT (1X,I3,3F10.4)                                             
      DO I=1,1200                                                     
      IDB(I)=0                                                          
      IW(I)=0                                                           
      IW(I*2)=0                                                         
      IW(I*3)=0                                                         
      IW(I*4)=0                                                         
      XM(I)=0.                                                          
      YM(I)=0.                                                          
      ZM(I)=0.
	enddo  
	
                                                             
  25  write(6,1008)model
 1008 format(////40x,'***',3x,'Model structure no. ',i2,
     1'   ***'//)                                                          
      NOP=0                                                             
      NOS=0                                                             
      NOM=0 
	                                                            
      DO 3 I=1,1200 

      if(ifile.eq.0) then
	READ(5,102)RX,RY,RZ,i1,i2,i3,ka,IW1,IW2,IW3,IW4,ISEP
  102 FORMAT (3E10.0,10X,4I3,8X,4I3,7X,I1) 

	else

	if(model.eq.1.or.ifile.eq.1)
     &READ(5,*)(iref(kk,i),kk=1,3),f(i),iref(4,i)
      
	     i1=iref(1,i)
	     i2=iref(2,i)
	     i3=iref(3,i)
           ka=iref(4,i)
	     if(i1.eq.0)isep=9
           CALL NEORED(RX,RY,RZ,I1,I2,I3,KA,F,ISEP,X,Y,Z)
           iw1=0   	     
	endif


   30 continue
      DO 20 J=1,1200
      IF (NOM.GE.NOS)  GO TO 21                                         
      NOM=NOM+1                                                         
   20 IDB(NOM)=IDA(NOS)                                                 
   21 IF(ISEP.NE.0)  GO TO 4                                            
      IF (KA.NE.0)  GO TO 5                                             
      NOP=NOP+1                                                         
      IDE(NOP)=0
	                                                                                              
      IF (f(i).NE.0.)GO TO 6                                               
      IF (I2.NE.0) GO TO 7                                              
      IF (I1.NE.0)  GO TO 8                                             
      GO TO 9                                                           
    6 RX=X(I1)+(X(I2)-X(I1))*f(i)                                         
      RY=Y(I1)+(Y(I2)-Y(I1))*f(i)                                          
      RZ=Z(I1)+(Z(I2)-Z(I1))*f(i)
      if(ifile.le.1.or.model.eq.1) write(6,400)nop,f(i),i1,i2
  400 format('Location(',i3,') between (',F4.1'%) atoms',i5,' and',i5)                                         
      GO TO 9                                                           
    7 RX=(X(I1)+X(I2))/2.                                               
      RY=(Y(I1)+Y(I2))/2.                                               
      RZ=(Z(I1)+Z(I2))/2. 
      if(ifile.le.1.or.model.eq.1) write(6,401)nop,i1,i2
  401 format('Location(',i3,') midway betweeen atoms ',i3,' and',i3)                                             
      GO TO 9                                                            
    8 RX=X(I1)                                                          
      RY=Y(I1)                                                          
      RZ=Z(I1)
      if(ifile.le.1.or.model.eq.1) write(6,402)nop,i1
  402 format('Location(',i3,') centered on atom ',i3)                                         
    
    9 XR(NOP)=RX                                                        
      YR(NOP)=RY                                                        
      ZR(NOP)=RZ 
c      write(6,103)nop,xr(nop),yr(nop),zr(nop)   
c  103 format('Location  (',i3,') [',2(f8.3,','),f8.3,']')
      GO TO 3  
	
    5 IF (I3.NE.0)  GO TO 10                                            
      IF (I2.NE.0)  GO TO 11                                            
      IF (I1.NE.0)  GO TO 12                                            
      GO TO 13                                                          
   10 DX3=X(I3)-X(I1)                                                   
      DX2=X(I2)-X(I1)                                                   
      DY3=Y(I3)-Y(I1)                                                   
      DY2=Y(I2)-Y(I1)                                                   
      DZ3=Z(I3)-Z(I1)                                                   
      DZ2=Z(I2)-Z(I1)                                                   
      RX=DZ3*DY2-DZ2*DY3                                                
      RY=DX3*DZ2-DX2*DZ3                                                
      RZ=DY3*DX2-DY2*DX3                                                
      GO TO 13                                                          
   11 RX=X(I2)-X(I1)                                                    
      RY=Y(I2)-Y(I1)                                                    
      RZ=Z(I2)-Z(I1)                                                    
      GO TO 13                                                          
   12 CONTINUE                                                          
      DX3=X(I1)-XR(NOP)                                                 
      DY3=Y(I1)-YR(NOP)                                                 
      DZ3=Z(I1)-ZR(NOP)                                                 
      DX2=DX                                                            
      DY2=DY                                                            
      DZ2=DZ                                                            
      DX1=DZ3*DY2-DZ2*DY3                                               
      DY1=DX3*DZ2-DX2*DZ3                                               
      DZ1=DY3*DX2-DY2*DX3                                               
      RX=DZ1*DY2-DZ2*DY1                                                
      RY=DX1*DZ2-DX2*DZ1                                                
      RZ=DY1*DX2-DY2*DX1
	                                                
   13 CONTINUE
      R=RX*RX+RY*RY+RZ*RZ
      R=DSQRT(R)                                                        
      DX=RX/R                                                           
      DY=RY/R                                                           
      DZ=RZ/R                                                           
      IF (IW1.NE.0)  GO TO 15  
      NOS=NOS+1                                                         
      IDE(NOP)=IDE(NOP)+1                                               
      IDA(NOS)=KA                                                       
      XE(NOS)=DX                                                        
      YE(NOS)=DY                                                        
      ZE(NOS)=DZ                                                        
c      write(6,104)nos,dx,dy,dz,ka   
c  104 format('Direction (',i3,') [',2(f8.3,','),f8.3,']; ',
c     &'Model polarizability ',i3/)	                        
      GO TO 3                                                           
   15 NOM=NOM+1                                                         
      IDB(NOM)=KA                                                       
      XM(NOS)=DX                                                        
      YM(NOS)=DY                                                        
      ZM(NOS)=DZ        
                     
    3 CONTINUE                                                          
    4 CONTINUE  
    
                                                            
      WRITE (6,206)                                                     
  206 FORMAT(//12X,20HOSCILLATOR LOCATIONS ,12X,26HELECTRIC OSC  UNIT VE
     1CTORS ,15X,26HMAGNETIC OSC  UNIT VECTORS  )              
      WRITE (6,207)                                                     
  207 FORMAT(/9X,1HX,9X,1HY,9X,1HZ,9X,3HOSC,2X,3HPOL,
     14X,1HX,9X,1HY,9X,1HZ,8X,3HOSC,2X,3HPOL,5X,1HX,9X,1HY,9X,1HZ  ) 
      K=0                                                               
      L=0                                                               
      DO 16 I=1,NOP
      WRITE (6,208)                                                     
      WRITE (6,211) I,XR(I),YR(I),ZR(I)         
  211 FORMAT (1H+,I3,3F10.4)                                     
      LINE=0                                                            
      NS=IDE(I)                                                         
      DO 17 J=1,NS                                                      
      K=K+1                                                             
      IF(LINE.EQ.1)  WRITE (6,208)                                      
      WRITE (6,209) K,IDA(K),XE(K),YE(K),ZE(K)                          
      IF (XM(K).EQ.0.AND.YM(K).EQ.0.AND.ZM(K).EQ.0) GO TO 18            
      L=L+1                                                             
      WRITE (6,210)L,IDB(K),XM(K),YM(K),ZM(K)                           
   18 LINE=1
  208 FORMAT (1X)                                                       
  209 FORMAT (1H+,35X,2I5,3F10.4)                                                                                                      
C  c'era  INSERZIONE CALCOLO   FATORE GEOMETRICO BINARIO
   17 CONTINUE                                                          
   16 CONTINUE                                                          

  210 FORMAT (1H+,75X,2I5,3F10.4)                                       
C  FINE alre INSERZIONE

      CALL ASHAPE(ISUP,ifile,model)
      RETURN
 9999 stop                                              
      END 
	
                                                                   
      subroutine neored(rx,ry,rz,i1,i2,i3,ka,theta,isep,x,y,z)
      IMPLICIT REAL*8(A-H,O-Z)
      DIMENSION X(1200),Y(1200),Z(1200)
C     I1 INDICE PRIMO PUNTO,I1I2 DIREZIONE DI RIFERIMENTO,I1I3 SECONDO
C     RIFERIMENTO,theta GRADI DI SPOSTAMENTO DA I1I2 A I1I3
      IF(KA.EQ.0) RETURN
      THETA=THETA*3.141593/180.
      COSA=DCOS(THETA)
      SENA=DSIN(THETA)
      AX=X(I2)-X(I1)
      AY=Y(I2)-Y(I1)
      AZ=Z(I2)-Z(I1)

	if(i3.ne.0) then
        BX=X(I3)-X(I1)
        BY=Y(I3)-Y(I1)
        BZ=Z(I3)-Z(I1)
        AB=AX*BX+AY*BY+AZ*BZ
        A=(AX*AX+AY*AY+AZ*AZ)**0.5
        B=(BX*BX+BY*BY+BZ*BZ)**0.5
        ZETA=DACOS(AB/(A*B))
        COSFI=DCOS(ZETA-THETA)
        AA=A*A
        BA=A*B
        B=(COSFI-COSA*AB/BA)/((1.-(AB/BA)**2)*B)
        A=COSA/A-B*AB/AA
	else
	  a=1.
	  b=0.
	endif

      RX=A*AX+B*BX
      RY=A*AY+B*BY
      RZ=A*AZ+B*BZ

      I1=0
      I2=0
      I3=0
      RETURN
      END


      SUBROUTINE ASHAPE(ISUP,ifile, model)
      IMPLICIT REAL*8(A-H,O-Z)
      COMMON idummy,NOP,NOS,IB,U,UL,UU,DU,IRUN,IPUNCH                          
      COMMON /AIR/  AI(100,1200),AR(100,1200)                               
      COMMON /POL/  AIT(100),ART(100),BM(100),UA(100)                       
      DIMENSION UE(100),UD(100),tt(10,100)  
	                                         
      DO 25 I=1,100                                                      
      UA(I)=0.                                                          
      UE(I)=0.                                                          
      UD(I)=0.                                                          
   25 BM(I)=0.                                                          
      WRITE (6,200)                                                     
  200 FORMAT(//1X,16HPOLARIZABILITIES  )  
                                
      DO 1 I=1,100                                                       
      If(ifile.le.1.or.model.eq.1)READ(5,*)(tt(m,i),m=1,7)
      Tua=tt(1,i)
      Tue=tt(2,i)
      Tud=tt(3,i)
      Tbm=tt(4,i)
      Ds=tt(5,i)
      Uk=tt(6,i)
      Gk=tt(7,i)
      IF(tua.eq.0)GO TO 2
      K=I                  
	WRITE (6,201) K,TUA,TUE,TUD,TBM                                   
  201 FORMAT( 1X,4HNBR=,I3,20X,4HNUA=,F7.3, 6H  NUE=,F7.3,6H  NUD=,F7.3,
     17H  BMAG=,F7.3  )                                                 
      UA(K)=TUA                                                         
      UE(K)=TUE                                                         
      UD(K)=TUD                                                         
      BM(K)=TBM                                                         
      ZZ=(TUE-TUA)/TUD                                                  
      IF (ZZ.LT.0.0)  ZZ=ZZ-1.0E-10                                     
      IF (ZZ.GT.0.0)  ZZ=ZZ+1.0E-10                                     
      IZ=ZZ                                                             
      IAB=IABS(IZ)+1                                                    
      IF (IAB.GT.201)  GO TO 5                                          
      IF (UK.NE.0)  GO TO 3                                             
      WRITE (6,202)                                                     
  202 FORMAT( 1X,7HREAD IN )                                            
      IF(ISEP.EQ.1) GO TO 7
      READ (5,101) (AI(K,J),J=1,IAB)                                    
  101 FORMAT(10F8.2)                                                    
      READ (5,101) (AR(K,J),J=1,IAB)                                    
      DO 1001 J= 1,IAB
      IF(AI(K,J).LT.0.001)AI(K,J)=0.001
C     AI(K,J)=AI(K,J)*0.738
C     AR(K,J)=AR(K,J)*0.738
 1001 CONTINUE
      GO TO 7                                                           
    3 IF (GK.EQ.0)  GK=0.5*DSQRT(UK)
      WRITE (6,203) DS,UK,GK                                            
  203 FORMAT(28X,9HDIPOLSTR= ,F9.5,6H  NUK=,F7.3,10H  DAMPING=,F7.3)    
      DO 6 J=1,IAB                                                      
      ZZ=J-1                                                            
      UDEL=ZZ*TUD                                                       
      IF (IZ.LT.0)  UDEL=-UDEL                                          
      U=TUA+UDEL                                                        
      UKU=UK*UK                                                         
      UU=U*U                                                            
      UMU=UKU-UU                                                        
      PRINT *, U, UK, UU, UMU
      UU=DS*UK*10.061/(UMU*UMU+UU*GK*GK)                                
      TBM=UMU*UU                                                        
      PRINT *, UU, TBM
      IF (DABS(TBM).GT.1000000.)  GO TO 5
      AR(K,J)=TBM                                                       
      TBM=U*GK*UU                                                       
      IF (TBM.GT.1000000.)  GO TO 5                                     
      AI(K,J)=TBM                                                       
    6 CONTINUE                                                          
    7 WRITE (6,204)                                                     
      IF(ISUP.EQ.8) GO TO 301
  204 FORMAT(1X,9HIMAGINARY )                                           
      WRITE (6,205) (AI(K,J),J=1,IAB)                                   
  205 FORMAT (1X,10F12.5)                                               
      WRITE (6,206)                                                     
  206 FORMAT(/1X,4HREAL  )                                              
      WRITE (6,205) (AR(K,J),J=1,IAB)                                   
  301 CONTINUE
      IF (IZ.GT.0)  GO TO 1                                             
       IF(ISEP.EQ.1.AND.IZ.LT.0) GO TO 299
      IA=IAB/2                                                          
      DO 4 J=1,IA                                                       
      JT=IAB+1-J                                                        
      TBM=AI(K,J)                                                       
      AI(K,J)=AI(K,JT)                                                  
      AI(K,JT)=TBM                                                      
      TBM=AR(K,J)                                                       
      AR(K,J)=AR(K,JT)                                                  
      AR(K,JT)=TBM                                                      
    4 CONTINUE                                                          
  299 CONTINUE
      UA(K)=TUE                                                         
      UE(K)=TUA                                                         
    1 CONTINUE                                                          
    5 WRITE (6,220)                                                     
  220 FORMAT(1X,16HPOL OUT OF RANGE )                                   
      STOP     
	                                                         
    2 UL=0.                                                             
      UU=10E5                                                           
      IB=0          
	                                                    
      DO 22 I=1,100                                                      
      IF (BM(I).NE.0.)  IB=1                                            
      UL=DMAX1(UL,UA(I))                                                
      UT=DMIN1(UU,UE(I))                                                
      IF(UT.NE.0.)  UU=UT                                               
      IF(UD(I).EQ.0.)  GO TO 22                                         
      IF (I.LT.10)  I1=I+1                                              
      DU=UD(I)                                                          
      IF(UD(I1).EQ.0.) GO TO 22                                         
      IF (UD(I).EQ.UD(I1))  GO TO 22                                    
      WRITE (6,270)                                                     
  270 FORMAT(1X,17HPOL INCOMPATIBLE  )                                  
      STOP                                                              
   22 CONTINUE                                                          
      RETURN                                                            
      END    
	
                                                                
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
	                                            
      SUBROUTINE INTRAC                                                 
      IMPLICIT REAL*8(A-H,O-Z)
      COMMON idummy,NOP,NOS,IB,U,UL,UU,DU,IRUN,IPUNCH                          
      COMMON /ROSC/ XR(1200),YR(1200),ZR(1200)                             
      COMMON /EEOS/ XE(1200),YE(1200),ZE(1200)                             
      COMMON /IDET/ IDA(1200),IDB(1200),IDE(1200),IW(4800)
      COMMON/IMU/noint(1200,1200),nnoint
      DIMENSION R(3)                                                    
      COMMON /CBL/ C(1440000)                                             
      DIMENSION G(1200,1200),v(1200,1200)                                                
      EQUIVALENCE (C(1),G(1,1))                                         
      DIMENSION EI(3),EJ(3)                                             
      DO 25 I=1,NOS                                                     
      DO 25 J=1,NOS                                                     
   25 G(I,J)=0.                                                         
      IF (IRUN.EQ.0)  WRITE (6,200)                                     
      NK=0                                                              
      DO 1 K=1,NOP                                                      
      NSK=IDE(K)                                                        
      IF(K.EQ.1)  GO TO 1                                               
      NL=0                                                              
      LMAX=K-1                                                          
      DO 2 L=1,LMAX                                                     
      NSL=IDE(L)                                                        
      R(1)=XR(L)-XR(K)                                                  
      R(2)=YR(L)-YR(K)                                                  
      R(3)=ZR(L)-ZR(K)                                                  
      DD=0.                                                             
      DO 5 I=1,3                                                        
    5 DD=DD+R(I)*R(I)                                                   
      D=DSQRT(DD)
      IF (IRUN.EQ.0)  WRITE (6,201) K,L,D                               
  200 FORMAT(///1X,9HDISTANCES  )                                       
  201 FORMAT (1X,2I8,10X,F10.4)                                         
c      DO 7 I=1,4                                                        
c      J=(L-1)*4+I                                                       
c      IF(IW(J).EQ.K)  GO TO 2                                           
c    7 CONTINUE                                                          
      R5=-3./(DD*DD*D)                                                  
      R3=1./(DD*D)                                                      
    8 DO 3 KI=1,NSK                                                     
      I=NK+KI                                                           
      EI(1)=XE(I)                                                       
      EI(2)=YE(I)                                                       
      EI(3)=ZE(I)                                                       
      DO 4 LJ=1,NSL                                                     
      J=NL+LJ                                                           
      EJ(1)=XE(J)                                                       
      EJ(2)=YE(J)                                                       
      EJ(3)=ZE(J)                                                       
      GG=0.                                                             
      DO 6 II=1,3                                                       
      DO 6 JJ=1,3                                                       
      TT=R(II)*R(JJ)*R5                                                 
      IF(II.EQ.JJ)  TT=TT+R3                                            
    6 GG=GG+EI(II)*TT*EJ(JJ)                                           
      G(I,J)=GG
      IF(nnoint.NE.0.AND.noint(I,J).EQ.1.or.noint(j,i).eq.1) G(I,J)=0.
      G(J,I)=G(I,J)                                                     
    4 CONTINUE                                                          
    3 CONTINUE                                                          
    2 NL=NL+NSL                                                         
    1 NK=NK+NSK                                                         
      IF (IRUN.NE.0)  GO TO 9                                           
      WRITE (6,202)                                                     
  202 FORMAT(//,1X,18HINTERACTION MATRIX   )                            
      CALL MATOUT (G,NOS)                                               
    9 RETURN                                                            
      END    
	
                                                                
      SUBROUTINE AMAT                                                   
      IMPLICIT REAL*8(A-H,O-Z)
      COMMON idummy,NOP,NOS,IB,U,UL,UU,DU,IRUN,IPUNCH                          
      COMMON /IDET/ IDA(1200),IDB(1200),IDE(1200),IW(4800)                  
      COMMON /POL/  AIT(100),ART(100),BM(100),UA(100)                       
      COMMON /CBL/ C(1440000)                                             
      COMPLEX*16 A(1200,1200),ARI
      DIMENSION G(1200,1200)
      EQUIVALENCE (C(1),G(1,1))                                         
      EQUIVALENCE (C(1),A(1,1))                                         
      DIMENSION S(2)                                                    
      EQUIVALENCE (S(1),ARI),(S(1),AR),(S(2),AI)  
	                      
      DO 1 K=1,NOS                                                      
      I=NOS+1-K                                                         
      DO 1 L=1,NOS                                                      
      J=NOS+1-L                                                         
    1  A(J,I)=G(J,I)                                                    
      DO 2 I=1,NOS                                                      
      K=IDA(I)                                                          
      AR=ART(K)                                                         
      AI=AIT(K)                                                         
      ARI=1./ARI                                                        
    2 A(I,I)=ARI                                                        
      CALL INVERT (A,NOS)                                               
      RETURN                                                            
      END           
	                                                    
      SUBROUTINE UVSHPE(KK,WUV)
C     INSERITA PER PLOT
      IMPLICIT REAL*8(A-H,O-Z)
      DOUBLE PRECISION WUV
      COMMON idummy,NOP,NOS,IB,U,UL,UU,DU,IRUN,IPUNCH                          
      COMMON /EEOS/ XE(1200),YE(1200),ZE(1200)                             
      DIMENSION E(3,3),TI(3,3),TR(3,3),W(3),EI(3),EJ(3)                 
      COMMON /CBL/ C(1440000)                                             
      COMPLEX*16 A(1200,1200),CI
      EQUIVALENCE (C(1),A(1,1))                                         
      DIMENSION S(2)                                                    
      DIMENSION WUV(KK)
C     INSERITA PER PLOT
      EQUIVALENCE (S(1),CI),(S(1),AR),(S(2),AI)                         
      DO 3 K=1,3                                                        
      DO 3 L=1,3                                                        
      TI(K,L)=0.                                                        
    3 TR(K,L)=0.                                                        
      DO 1 I=1,NOS                                                      
      EI(1)=XE(I)                                                       
      EI(2)=YE(I)                                                       
      EI(3)=ZE(I)                                                       
      DO 1 J=1,NOS                                                      
      CI=A(I,J)                                                         
      EJ(1)=XE(J)                                                       
      EJ(2)=YE(J)                                                       
      EJ(3)=ZE(J)                                                       
      DO 2 K=1,3                                                        
      DO 2 L=1,3                                                        
      T=EI(K)*EJ(L)                                                     
      TI(K,L)=TI(K,L)+T*AI                                              
      TR(K,L)=TR(K,L)+T*AR                                              
    2 CONTINUE                                                          
    1 CONTINUE                                                          
      V=-6.89*U                                                         
C      IF (IPUNCH.NE.0)  CALL PUNCH (TI,V,U)                             
      CALL DIAG (TI,E,W,3)                                              
  200 FORMAT(1H+,9X,'EX'    )
      CALL PRINT (W,E,V)                                                
C     WRITE (6,200)                                                     
      WW=W(1)+W(2)+W(3)
      WUV(KK)=WW
C     INSERITA PER PLOT
      V=0.838                                                           
C      IF (IPUNCH.EQ.2)  CALL PUNCH (TR,V,U)                             
      CALL DIAG (TR,E,W,3)                                              
  201 FORMAT(1H+,9X,'RF'   )
      CALL PRINT (W,E,V)                                                
C     WRITE (6,201)                                                     
      RETURN                                                            
      END               
	                                                
      SUBROUTINE RTSHPE(KWK,WCD,WORD)
C     INSERITA PER PLOT
      IMPLICIT REAL*8(A-H,O-Z)
      DOUBLE PRECISION WCD
      COMMON idummy,NOP,NOS,IB,U,UL,UU,DU,IRUN,IPUNCH                          
      COMMON /ROSC/ XR(1200),YR(1200),ZR(1200)                             
      COMMON /EEOS/ XE(1200),YE(1200),ZE(1200)                             
      COMMON /EMOS/ XM(1200),YM(1200),ZM(1200)                             
      COMMON /POL/  AIT(100),ART(100),BM(100),UA(100)                       
      COMMON /IDET/ IDA(1200),IDB(1200),IDE(1200),IW(4800)                  
      DIMENSION E(3,3),TI(3,3),TR(3,3),W(3),EI(3),EJ(3),EC(3),EM(3),    
     1D(3),TIM(3,3),TRM(3,3),TIP(3,3),TRP(3,3)                          
      COMMON /CBL/ C(1440000)                                             
      COMPLEX*16 A(1200,1200),CI
      EQUIVALENCE (C(1),A(1,1))                                         
      DIMENSION S(2)                                                    
      DIMENSION WCD(KWK)
      DIMENSION WORD(KWK)
C     INSERITA PER PLOT
      EQUIVALENCE (S(1),CI),(S(1),AR),(S(2),AI)                         
      DO 5 K=1,3                                                        
      DO 5 L=1,3                                                        
      TI(K,L)=0.                                                        
      TIP(K,L)=0.                                                       
      TIM(K,L)=0.                                                       
      TR(K,L)=0.                                                        
      TRP(K,L)=0.                                                       
    5 TRM(K,L)=0.                                                       
      NK=0                                                              
      DO 1 K=1,NOP                                                      
      NSK=IDE(K)                                                        
      NL=0                                                              
      DO 2 L=1,K                                                        
      NSL=IDE(L)                                                        
      D(1)=XR(L)-XR(K)                                                  
      D(2)=YR(L)-YR(K)                                                  
      D(3)=ZR(L)-ZR(K)                                                  
      DO 3 KI=1,NSK                                                     
      I=NK+KI                                                           
      EI(1)=XE(I)                                                       
      EI(2)=YE(I)                                                       
      EI(3)=ZE(I)                                                       
      DO 4 LJ=1,NSL                                                     
      J=NL+LJ                                                           
      CI=A(I,J)                                                         
      IF (I.NE.J)  GO TO 10                                             
      AI=AI/2.                                                          
      AR=AR/2.                                                          
   10 CONTINUE                                                          
      EJ(1)=XE(J)                                                       
      EJ(2)=YE(J)                                                       
      EJ(3)=ZE(J)                                                       
      EM(1)=XM(J)                                                       
      EM(2)=YM(J)                                                       
      EM(3)=ZM(J)                                                       
      IM=IDB(J)                                                         
      B=-BM(IM)*4.                                                      
      DO 6 II=1,3                                                       
      M=II/3                                                            
      KK=II+1-M*3                                                       
      M=(II+1)/3                                                        
      LL=II+2-M*3                                                       
      EC(II)=EI(KK)*EJ(LL)-EI(LL)*EJ(KK)                                
      DO 7 JJ=1,3                                                       
      TP=EC(II)*D(JJ)                                                   
      TM=EI(II)*EM(JJ)*B                                                
      TIP(II,JJ)=TIP(II,JJ)+TP*AI                                       
      TIM(II,JJ)=TIM(II,JJ)+TM*AI                                       
      TRP(II,JJ)=TRP(II,JJ)+TP*AR                                       
      TRM(II,JJ)=TRM(II,JJ)+TM*AR                                       
    7 CONTINUE                                                          
    6 CONTINUE                                                          
    4 CONTINUE                                                          
    3 CONTINUE                                                          
    2 NL=NL+NSL                                                         
    1 NK=NK+NSK                                                         
      DO 11 II=1,3                                                      
      DO 11 JJ=1,II                                                     
      TIP(II,JJ)=TIP(II,JJ)+TIP(JJ,II)                                  
      TIM(II,JJ)=TIM(II,JJ)+TIM(JJ,II)                                  
      TRP(II,JJ)=TRP(II,JJ)+TRP(JJ,II)                                  
      TRM(II,JJ)=TRM(II,JJ)+TRM(JJ,II)                                  
      TIP(JJ,II)=TIP(II,JJ)                                             
      TIM(JJ,II)=TIM(II,JJ)                                             
      TRP(JJ,II)=TRP(II,JJ)                                             
      TRM(JJ,II)=TRM(II,JJ)                                             
      TI(II,JJ)=TIP(II,JJ)+TIM(II,JJ)                                   
      TR(II,JJ)=TRP(II,JJ)+TRM(II,JJ)                                   
      TI(JJ,II)=TI(II,JJ)                                               
      TR(JJ,II)=TR(II,JJ)                                               
   11 CONTINUE                                                          
      IF(IB.EQ.0)  GO TO 8                                              
      WRITE(6,200)                                                      
      V=U*U*0.000433                                                    
      CALL DIAG (TIP,E,W,3)                                             
      CALL PRINT (W,E,V)                                                
      WRITE(6,201)                                                      
      CALL DIAG(TIM,E,W,3)                                              
      CALL PRINT (W,E,V)                                                
      WRITE(6,202)                                                      
      V=-V*3300.                                                        
      CALL DIAG (TRP,E,W,3)                                             
      CALL PRINT (W,E,V)                                                
      WRITE(6,203)                                                      
      CALL DIAG(TRM,E,W,3)                                              
      CALL PRINT (W,E,V)                                                
  200 FORMAT (1H+,9X,'CP')                                              
  201 FORMAT(1H+,9X,'CM')                                               
  202 FORMAT(1H+,9X,'RP')                                               
  203 FORMAT(1H+,9X,'RM')                                               
  8   CONTINUE
  206 FORMAT(1H+,9X,'CD')                                               
      V=U*U*0.000433                                                    
C      IF(IPUNCH.NE.0) CALL PUNCH (TI,V,U)                               
      CALL DIAG(TI,E,W,3)                                               
      CALL PRINT(W,E,V)                                                 
C     WRITE(6,206)
      WW=W(1)+W(2)+W(3)
      WCD(KWK)=WW
C     INSERITA PER PLOT
  205 FORMAT(1H+,9X,'RD')                                               
      V=-V*3300.                                                        
C      IF(IPUNCH.EQ.2) CALL PUNCH (TR,V,U)                               
      CALL DIAG(TR,E,W,3)                                               
      CALL PRINT(W,E,V)                                                 
C     WRITE(6,205)                                                      
      WW=W(1)+W(2)+W(3)
       WORD(KWK)=WW
      RETURN                                                            
      END           
	                                                    
      SUBROUTINE INVERT(T,N)                                            
      IMPLICIT REAL*8(A-H,O-Z)
       COMPLEX*16 T(1200,1200),R(1200),S(1200),D
      DO 100 L=1,N                                                      
      D=1./T(L,L)                                                       
      DO 30 I=1,N                                                       
      IF(L-I) 15,20,25                                                  
   15 R(I)=T(I,L)*D                                                     
      S(I)=T(I,L)                                                       
      GO TO 30                                                          
   20 R(I)=-D                                                           
      GO TO 30                                                          
   25 R(I)=T(L,I)*D                                                     
      S(I)=T(L,I)                                                       
   30 CONTINUE                                                          
      DO 100 I=1,N                                                      
      DO 100 J=1,I                                                      
      IF (I-L) 35,50,35                                                 
   35 IF (J-L) 40,45,40                                                 
   40 T(I,J)=T(I,J)-R(I)*S(J)                                           
      GO TO 100                                                         
   45 T(I,J)=R(I)                                                       
      GO TO 100                                                         
   50 T(I,J)=R(J)                                                       
  100 CONTINUE                                                          
      DO 200 I=1,N                                                      
      DO 200 J=1,I                                                      
      T(I,J)=-T(I,J)                                                    
  200 T(J,I)=T(I,J)                                                     
      RETURN                                                            
      END                   
	                                            
      SUBROUTINE DIAG (A,EIVR,E,N)                                      
      IMPLICIT REAL*8(A-H,O-Z)
      DIMENSION A(3,3),EIVR(3,3),E(3)                                   
      IEGEN=0                                                           
      IF(N-1) 2,2,1                                                     
    2 EIVR(1,1)=1.0                                                     
      E(1)=A(1,1)                                                       
      RETURN                                                            
    1 IF(IEGEN) 102,99,102                                              
   99 DO 101 J=1,N                                                      
      DO 100 I=1,N                                                      
  100 EIVR(I,J)=0.0                                                     
  101 EIVR(J,J)=1.0                                                     
      DO 888 I=1,N                                                      
      DO 888 J=1,N                                                      
      IF (DABS(A(I,J)).GT.10E-8)  GO TO 887
  888 CONTINUE                                                          
      GO TO 134                                                         
  887 CONTINUE                                                          
  102 ATOP=0.                                                           
      DO 111 I=1,N                                                      
      DO 111J=I,N                                                       
      IF(ATOP-DABS(A(I,J)))104,111,111
  104 ATOP=DABS(A(I,J))
  111 CONTINUE                                                          
      IF(ATOP)109,109,113                                               
  109 RETURN                                                            
  113 AVGF=DFLOAT(N*(N-1))*.55
      D=0.0                                                             
      DO 114 JJ=2,N                                                     
      DO 114 II=2,JJ                                                    
      S=A(II-1,JJ)/ATOP                                                 
  114 D=S*S+D                                                           
      DSTOP=(1.E-06)*D                                                  
      THRSH =DSQRT(D/AVGF)*ATOP                                         
  115 IFLAG=0                                                           
      DO 130 JCOL=2,N                                                   
      JCOL1=JCOL-1                                                      
      DO 130 IROW=1,JCOL1                                               
      AIJ=A(IROW,JCOL)                                                  
      IF(DABS(AIJ)-THRSH)130,130,117
  117 AII=A(IROW,IROW)                                                  
      AJJ=A(JCOL,JCOL)                                                  
      S=AJJ-AII                                                         
      IF(DABS(AIJ)-1.E-09*DABS(S))130,130,118
  118 IFLAG=1                                                           
      IF(1.E-10*DABS(AIJ)-DABS(S))116,119,119
  119 S=.70710678118655                                                 
      C=S                                                               
      GO TO 120                                                         
  116 T=AIJ/S                                                           
      S=0.25/DSQRT(0.25+T*T)
      C=DSQRT(0.5+S)
      S=2.*T*S/C                                                        
  120 DO 121 I=1,IROW                                                   
      T=A(I,IROW)                                                       
      U=A(I,JCOL)                                                       
      A(I,IROW)=C*T-S*U                                                 
  121 A(I,JCOL)=S*T+C*U                                                 
      I2=IROW+2                                                         
      IF(I2-JCOL)127,127,123                                            
  127 CONTINUE                                                          
      DO 122 I=I2,JCOL                                                  
      T=A(I-1,JCOL)                                                     
      U=A(IROW,I-1)                                                     
      A(I-1,JCOL)=S*U+C*T                                               
  122 A(IROW,I-1)=C*U-S*T                                               
  123 A(JCOL,JCOL)=S*AIJ+C*AJJ                                          
      A(IROW,IROW)=C*A(IROW,IROW)-S*(C*AIJ-S*AJJ)                       
      DO 124 J=JCOL,N                                                   
      T=A(IROW,J)                                                       
      U=A(JCOL,J)                                                       
      A(IROW,J)=C*T-S*U                                                 
  124 A(JCOL,J)=S*T+C*U                                                 
      IF(IEGEN) 126,131,126                                             
  131 DO 125 I=1,N                                                      
      T=EIVR(I,IROW)                                                    
      EIVR(I,IROW)=C*T-EIVR(I,JCOL)*S                                   
  125 EIVR(I,JCOL)=S*T+EIVR(I,JCOL)*C                                   
  126 CONTINUE                                                          
      S=AIJ/ATOP                                                        
      D=D-S*S                                                           
      IF(D-DSTOP)1260,129,129                                           
 1260 D=0.                                                              
      DO 128 JJ=2,N                                                     
      DO 128 II=2,JJ                                                    
      S=A(II-1,JJ)/ATOP                                                 
  128 D=S*S+D                                                           
      DSTOP=(1.E-06)*D                                                  
  129 THRSH=DSQRT(D/AVGF)*ATOP
  130 CONTINUE                                                          
      IF(IFLAG)115,134,115                                              
  134 DO 998 I=1,N                                                      
  998 E(I)=A(I,I)                                                       
      RETURN                                                            
      END                
	                                               
      SUBROUTINE MATOUT (A,N)                                           
      IMPLICIT REAL*8(A-H,O-Z)
      DIMENSION A(1200,1200)                                                
      K=0                                                               
    1 L=K+1                                                             
      K=K+15                                                            
      K=MIN0(K,N)                                                       
      WRITE (6,2) (I,I=L,K)                                             
    2 FORMAT (1H0//2X,15(6X,I2)//)                                      
      DO 3 I=1,N                                                        
    3 WRITE (6,4) I,(A(I,J),J=L,K)                                      
    4 FORMAT(1X,I4,2X,15F8.4)                                           
      IF(K.LT.N)  GO TO 1                                               
      RETURN                                                            
      END                     
	                                          
      SUBROUTINE PRINT (W,E,V)                                        
      IMPLICIT REAL*8(A-H,O-Z)
      DIMENSION W(3),E(3,3)                                             
      WW=0.                                                             
      DO 1 I=1,3                                                        
      W(I)=W(I)*V                                                       
      WW=WW+W(I)                                                        
      DO 1 J=1,3                                                        
      IF (DABS(E(I,J)).GT.1.)  E(I,J)=1.
    1 E(I,J)=DACOS(E(I,J))*57.3
C     WRITE (6,200) WW,((W(J),(E(I,J),I=1,3)),J=1,3)                    
  200 FORMAT (12X,F12.3,3X,3(F12.3,2X,3F6.1))                           
      RETURN                                                            
      END            
	                                                   
      SUBROUTINE PUNCH (T,V,U)                                          
      IMPLICIT REAL*8(A-H,O-Z)
      DIMENSION T(3,3),W(6)                                             
      JI=0                                                              
      DO 1 I=1,3                                                        
      DO 1 J=1,I                                                        
      JI=JI+1                                                           
    1 W(JI)=T(J,I)*V                                                    
      WRITE (7,100) U,(W(JI),JI=1,6)                                    
  100 FORMAT (F8.3,6F12.2)                                              
      RETURN                                                            
      END
	                                                               
      SUBROUTINE PRNPLT(X,Y,XMAX,XINCR,YMAX,YINCR,ISX,ISY,NPTS)
      DOUBLE PRECISION X,Y,XMAX,XINCR,YMAX,YINCR
C  PRINTER PLOT ROUTINE       M.S.ITZKOWITZ    MAY,1967
C
C   PLOTS THE @NPTS@ POINTS GIVEN BY @X(I),Y(I)@ ON A 51 X 101 GRID
C   USING A TOTAL OF 56 LINES ON THE PRINTER
C     IF @ISX@ OR @ISY@ ARE NON-ZERO, THE CORRESPONDING MAXIMUM AND
C      INCREMENTAL STEP SIZE ARE COMPUTED
C     IF EITHER INCREMENTAL STEP SIZE IS ZERO, THE PROGRAM EXITS
C     NEITHER OF THE INPUT ARRAYS ARE DESTROYED.  IF SCALING IS DONE
C     THE CORRESPONDING NEW VALUES OF MAXIMUM AND STEP SIZE ARE RETURNED
C
      DIMENSION X(NPTS),Y(NPTS),IGRID(105),XAXIS(11)
C
      INTEGER BLANK,DOT,STAR,IGRID,PLUS
      DATA BLANK,DOT,STAR,PLUS / 1H ,1H.,1H*,1H. /
C
901   FORMAT(14X,105A1)
902   FORMAT(1XE10.3,2X,1H+,105A1,1H+)
903   FORMAT(15X,103(1H.))
904   FORMAT(7X,11(F10.0),2H (,I4,5H PTS)  )
905   FORMAT(16X,11(1H+,9X))
9800  FORMAT(46H1SCALING ERROR IN PRNPLT, EXECUTION TERMINATED )
       WRITE(6,899)
  899 FORMAT(1H1)
C
      IF(ISX.NE.0) CALL PLSCAL(X,XMAX,XINCR,NPTS,100)
      IF(ISY.NE.0) CALL PLSCAL(Y,YMAX,YINCR,NPTS,50)
      IF(XINCR.EQ.0..OR.YINCR.EQ.0.) GO TO 800
      YAXMIN=0.01*YINCR
      XAXMIN=0.01*XINCR
      IZERO=YMAX/YINCR+1.5
      JZERO=103.5-XMAX/XINCR
      IF(JZERO.GT.103.OR.JZERO.LT.4) JZERO=2
      WRITE(6,905)
      WRITE(6,903)
      DO 10 I=1,51
      IF ( I.NE.IZERO) GO TO 16
      DO 14 J=1,105
14    IGRID(J)=PLUS
      GO TO 15
16    DO 11 J=1,105
11    IGRID(J)=BLANK
15    IGRID(JZERO)=PLUS
      IGRID(104)=DOT
      IGRID(2)=DOT
      DO 12 K=1,NPTS
      ITEST =(YMAX-Y(K))/YINCR+1.5
      IF(ITEST .NE.I) GO TO 12
      J=103.5-(XMAX-X(K))/XINCR
      IF(J.GT.103)J=105
      IF(J.LT.3) J=1
      IGRID(J)=STAR
12    CONTINUE
      IF(MOD(I,10).EQ.1) GO TO 13
      WRITE(6,901) IGRID
      GO TO 10
13    YAXIS=YMAX-(I-1)*YINCR
      IF(ABS(YAXIS).LT.YAXMIN) YAXIS=0.
      WRITE(6,902) YAXIS,(IGRID(J),J=1,105)
10    CONTINUE
      WRITE(6,903)
      WRITE(6,905)
      DO 20 M=1,11
      XAXIS(M)=XMAX-XINCR*(FLOAT(11-M))*10.0
      IF(ABS(XAXIS(M)).LT.XAXMIN)XAXIS(M)=0.
20    CONTINUE
      WRITE(6,904) XAXIS,NPTS
      RETURN
  800 WRITE(6,9800)
      STOP
      END

      SUBROUTINE PLSCAL(V,VMAX,VINCR,NPTS,NDIVIS)
      DOUBLE PRECISION V,VMAX,VINCR
C
C    SCALING PROGRAM FOR USE WITH PRNPLT   M.S.ITZKOWITZ  MAY,1967
C  THIS VERSION ADJUSTS THE FULL SCALE TO 2.5,5.0, OR 10. TIMES 10**N
C  AND ADJUSTS THE MAXIMUM POINT TO AN INTEGER MULTIPLE OF 5*VINCR
C
      DIMENSION V(NPTS)
C
   98 FORMAT(2X,10F12.6)
   99 FORMAT (1X,2I5,2F10.7)
      VMIN=V(1)
      VMAX=V(1)
      DO 10 I=1,NPTS
      IF(V(I).LT.VMIN) VMIN=V(I)
      IF(V(I).GT.VMAX) VMAX=V(I)
      QRANGE=VMAX-VMIN
10    CONTINUE
      IF(QRANGE.EQ.0.) GO TO 8000
      QRANGE=0.4342944*ALOG(QRANGE)
      IF(QRANGE)20,20,30
30    IRANGE=QRANGE
      GO TO 40
20    IRANGE=-QRANGE
      IRANGE=-IRANGE-1
40    QRANGE=QRANGE-FLOAT(IRANGE)
      RANGE=10.**QRANGE
C
C   RANGE IS BETWEEN 1.0 AND 10.0
C
43    IF(RANGE.GT.2.5) GO TO 41
      RANGE=2.5
      GO TO 50
41    IF(RANGE.GT.5.0) GO TO 42
      RANGE=5.0
      GO TO 50
42    RANGE=10.0
50    TRANGE=RANGE*(10.**IRANGE)
C   TRANGE IS NOW 2.5,5.0, OR 10.0 TIMES A POWER OF TEN
C
      VINCR=TRANGE/FLOAT(NDIVIS)
      IF(VMAX)51,51,52
52    IMAX=VMAX/(5.0*VINCR)
      XMAX=5.0*VINCR*FLOAT(IMAX+1)
      GO TO 53
51    IMAX=-VMAX/(5.0*VINCR)
      XMAX=5.0*VINCR*FLOAT(-IMAX+1)
53    IF(VMIN.GT.XMAX-TRANGE) GO TO 100
      RANGE=RANGE*2.0
      IF(RANGE-10.) 43,43,54
54    RANGE=RANGE/10.
      IRANGE=IRANGE+1
      GO TO 43
100   VMAX=XMAX
      VMIN=XMAX-TRANGE
      RETURN
 8000 WRITE(6,9800)
9800  FORMAT(45H1PLSCAL CALLED TO SCALE ARRAY WITH ZERO RANGE)
      STOP
      END
