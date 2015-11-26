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
