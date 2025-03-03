C      ######################################################################
C      #################      CAE Assistant Company          ################
C      ##############         www CAEassistant com              #############
C      ###########   Copy right by CAE Assistant Company    ###############
C      ######################################################################
C      ONLY the BUYER  of this package has permission to use its codes.
C	 Any distribution of this subroutine is illegal and will be prosecuted 
C      ######################################################################
C      ######################################################################
C      CAE Assisitant Services: 
C      Toturial Packages,Consultancy,Articles,Q&A,Video Gallery,Online Course
C      ######################################################################
C      Need help with your project? 
C      You can get initial free consultation from (Support CAEassistant com)
C      ######################################################################
	   SUBROUTINE UMAT(STRESS,STATEV,DDSDDE,SSE,SPD,SCD,
     1 RPL,DDSDDT,DRPLDE,DRPLDT,
     2 STRAN,DSTRAN,TIME,DTIME,TEMP,DTEMP,PREDEF,DPRED,CMNAME,
     3 NDI,NSHR,NTENS,NSTATV,PROPS,NPROPS,COORDS,DROT,PNEWDT,
     4 CELENT,DFGRD0,DFGRD1,NOEL,NPT,LAYER,KSPT,KSTEP,KINC)
C
      INCLUDE 'ABA_PARAM.INC'
C
      CHARACTER*80 CMNAME
      DIMENSION STRESS(NTENS),STATEV(NSTATV),
     1 DDSDDE(NTENS,NTENS),DDSDDT(NTENS),DRPLDE(NTENS),
     2 STRAN(NTENS),DSTRAN(NTENS),TIME(2),PREDEF(1),DPRED(1),
     3 PROPS(NPROPS),COORDS(3),DROT(3,3),DFGRD0(3,3),DFGRD1(3,3)

c============UMAT_MAT3==================================


	   IF (CMNAME(1:4) .EQ. 'USER') THEN
		CALL UMAT_MAT3(STRESS,STATEV,DDSDDE,SSE,SPD,SCD,
     1 RPL,DDSDDT,DRPLDE,DRPLDT,
     2 STRAN,DSTRAN,TIME,DTIME,TEMP,DTEMP,PREDEF,DPRED,CMNAME,
     3 NDI,NSHR,NTENS,NSTATV,PROPS,NPROPS,COORDS,DROT,PNEWDT,
     4 CELENT,DFGRD0,DFGRD1,NOEL,NPT,LAYER,KSPT,KSTEP,KINC)	 
         END IF
      print*,'matname=',CMNAME(1:4)
      RETURN
      END	  
c==========================================================
c=========== subroutine=============================
      SUBROUTINE UMAT_MAT3(STRESS,STATEV,DDSDDE,SSE,SPD,SCD,
     1 RPL,DDSDDT,DRPLDE,DRPLDT,
     2 STRAN,DSTRAN,TIME,DTIME,TEMP,DTEMP,PREDEF,DPRED,CMNAME,
     3 NDI,NSHR,NTENS,NSTATV,PROPS,NPROPS,COORDS,DROT,PNEWDT,
     4 CELENT,DFGRD0,DFGRD1,NOEL,NPT,LAYER,KSPT,KSTEP,KINC)
C
      INCLUDE 'ABA_PARAM.INC'
C
      CHARACTER*80 CMNAME
      DIMENSION STRESS(NTENS),STATEV(NSTATV),
     1 DDSDDE(NTENS,NTENS),DDSDDT(NTENS),DRPLDE(NTENS),
     2 STRAN(NTENS),DSTRAN(NTENS),TIME(2),PREDEF(1),DPRED(1),
     3 PROPS(NPROPS),COORDS(3),DROT(3,3),DFGRD0(3,3),DFGRD1(3,3)

      real*4 E_11,E_22,E_33,
     1 G_12,G_23,G_31, 
     1 nu_12,nu_23,!nu31,
     1 !nu21,nu32,nu13,
     1 !nu123,
     1 CC(6,6),
     1 E11_0,  E22_0,  E33_0 ,  G12_0, G23_0, G31_0,
     1 E11_f,  E22_f,  E33_f ,  G12_f, G23_f, G31_f,
     1 XT_0, XC_0, YT_0, YC_0, ZT_0, ZC_0, S12_0, S23_0, S31_0,
     1 XT, XC, XX,XXX,
     1 YT, YC, YY,YYY,
     1 ZT, ZC, ZZ,ZZZ,
     1 S12, S23, S31,
     1 eps(6), sig(6),
     1 sig_max(6),sig_min(6),
     1 sig_m(6), sig_a(6), 
     1 Fiber_Failure, Matrix_Failure,
     1 fatigue_f, fatigue_m,
     1  RR(6),
     1 sig_a_eq(6), 
     1 Nf(6),Nf_norm(6),NN,
     1 cycletime,
     1 DW_star,unifiedNf,
     1 Alfa1T,Alfa1C,Alfa2T,Alfa2C,Alfa3T,Alfa3C,Alfa4,Alfa5,Alfa6,
     1 Beta1T,Beta1C,Beta2T,Beta2C,Beta3T,Beta3C,Beta4,Beta5,Beta6,
     1 Landa1,Landa2,Landa3,Landa4,Landa5,Landa6,
     1 Gamma1,Gamma2,Gamma3,Gamma4,Gamma5,Gamma6
     
     
      integer i,j,FiberStaticFailure,MatrixStaticFailure
 
      metricpar=1.d0
c     Material elastic properties for intact material
      E11_0 = 147000.0*metricpar  
      ! Hidden Lines
      S31_0=S12_0
      
c     final strains      
      ! Hidden Lines
      eps31_f = 0.101   

c     Stiffness degradation parameters
      Landa1=14.57    ;      Gamma1=0.3024
      ! Hidden Lines
      
c     Strength degradation parameters
      Alfa1T=8.867   ;     Beta1T=0.545
      ! Hidden Lines
      
      timeperiod=1.d0                     
            
      IF (TIME(2).lt. timeperiod) then 
    
c     Material elastic properties 
      E11_0 = 147000.0*metricpar   
      ! Hidden Lines

c     Material strength properties
      XT=XT_0    
      XC=XC_0
      YT=YT_0    
      YC=YC_0
      ZT=ZT_0 
      ZC=ZC_0 
      S12=S12_0   
      S23=S23_0
      S31=S31_0
      
c     zero values    
      sig_max(:)=0
      sig_min(:)=0
      sig_m(:)=0
      sig_a(:)=0
      RR(:)=0
      sig_a_eq(:)=0
      
      Nf_norm(:)=0
      Nf(:)=0
      
      ! extra
      Fiber_Failure=0   !extra
      Matrix_Failure=0  !extra
      DW_star=0
      unifiedNf=0
      FiberStaticFailure=0
      MatrixStaticFailure=0
      
c     Life & time for the fist fatigue step
      NN=100.0  !wil be multiplied by 10 in the fatigue step
      cycletime=1.d0 
      
      
      ELSE   ! ---------------------------------------------------------
      
      E_11=STATEV(7)
      E_22=STATEV(8)
      E_33=STATEV(9)
      G_12=STATEV(10)
      G_23=STATEV(11)
      G_31=STATEV(12)
      nu_12=STATEV(13) 
      nu_23=STATEV(14)
      
      XT=STATEV(15) 
      XC=STATEV(16)
      YT=STATEV(17)    
      YC=STATEV(18) 
      ZT=STATEV(19)    
      ZC=STATEV(20)        
      S12=STATEV(21)   
      S23=STATEV(22)
      S31=STATEV(23)
      
      sig_max(1)= STATEV(24)
      sig_max(2)= STATEV(25)
      sig_max(3)= STATEV(26)
      sig_max(4)= STATEV(27)
      sig_max(5)= STATEV(28)
      sig_max(6)= STATEV(29)
      
      sig_min(1)= STATEV(30)
      sig_min(2)= STATEV(31)
      sig_min(3)= STATEV(32)
      sig_min(4)= STATEV(33)
      sig_min(5)= STATEV(34)
      sig_min(6)= STATEV(35)
     
      ! Aded in 14.for (no use)
      sig_m(1)=STATEV(36)
      sig_m(2)=STATEV(37)
      sig_m(3)=STATEV(38)
      sig_m(4)=STATEV(39)
      sig_m(5)=STATEV(40)
      sig_m(6)=STATEV(41)
      
      
      sig_a(1)=STATEV(42)
      sig_a(2)=STATEV(43)
      sig_a(3)=STATEV(44)
      sig_a(4)=STATEV(45)
      sig_a(5)=STATEV(46)
      sig_a(6)=STATEV(47)
 
      RR(1)=STATEV(48)
      RR(2)=STATEV(49)
      RR(3)=STATEV(50)
      RR(4)=STATEV(51)
      RR(5)=STATEV(52)
      RR(6)=STATEV(53)
      
      
      sig_a_eq(1)=STATEV(54)
      sig_a_eq(2)=STATEV(55)
      sig_a_eq(3)=STATEV(56)
      sig_a_eq(4)=STATEV(57)
      sig_a_eq(5)=STATEV(58)
      sig_a_eq(6)=STATEV(59)  
      
      Nf_norm(1)=STATEV(60)
      Nf_norm(2)=STATEV(61)
      Nf_norm(3)=STATEV(62)
      Nf_norm(4)=STATEV(63)
      Nf_norm(5)=STATEV(64)
      Nf_norm(6)=STATEV(65)
      
       Nf(1)=STATEV(66)
       Nf(2)=STATEV(67)
       Nf(3)=STATEV(68)
       Nf(4)=STATEV(69)
       Nf(5)=STATEV(70)
       Nf(6)=STATEV(71)
       
       Fiber_Failure=STATEV(72)  
       Matrix_Failure= STATEV(73)
       DW_star=STATEV(74)
       unifiedNf=STATEV(75)
       FiberStaticFailure= STATEV(76)
       MatrixStaticFailure=STATEV(77)
       

c     For information
      cycletime=STATEV(79)      
      NN=STATEV(80)
                 
      END IF
      
     
c     Engineering total strain 
      eps(1)=STATEV(1) + DSTRAN(1)
      eps(2)=STATEV(2) + DSTRAN(2)
      eps(3)=STATEV(3) + DSTRAN(3)
      eps(4)=STATEV(4) + DSTRAN(4)
      eps(5)=STATEV(5) + DSTRAN(5) 
      eps(6)=STATEV(6) + DSTRAN(6)
           
      CALL StiffnesMatrix(E_11,E_22,E_33, 
     1                    G_12,G_23, G_31,
     1                    nu_12,nu_23, CC )
     
         
        sig(:)=0      !Setting zero for confidence
         DO I=1,6
           do J=1,6
          sig(I) = sig(I) + 
     1                      CC(I,J)*eps(J)  ! W is not needed
           end do
         END DO
C      print*,'sig(I)**=',sig(I)   
c  ****************************************************************************             
c  *******************  failure analysis block ************************
c     Checking Failure: 

c     0- Previous failure
      if (FiberStaticFailure==1.OR. MatrixStaticFailure==1) then
      go to 20
      end if
C      print*,'FiberStaticFailure**=',FiberStaticFailure 

c     1- Fiber  failure  
      call FailureCriterionFiber(sig(1), sig(4), sig(6)    
     1              , XT, XC, S12, Fiber_Failure ) 
      
     
       if (Fiber_Failure>0) then
       call SuddenDegFiber( E_11,E_22,G_12,
     1                        nu_12,nu_23,
     1                     XT,XC,YT,YC,S12,S23 )
       FiberStaticFailure=1
       go to 20     ! Go to the end of the program
       end if     
                                                 
c     2- Transverse  failure  
      call FailureCriterionMatrix 
     1           ( sig(2), sig(3),sig(4), sig(5), sig(6),
     1           YT, YC, S12, S23, Matrix_Failure )
     
C       print*,'Matrix_Failure**=',Matrix_Failure 
       if (Matrix_Failure>0) then   
      call SuddenDegMatrix (E_22,G_12,nu_12,nu_23
     1                    ,YT,YC,S12,S23  )  
        MatrixStaticFailure=1  
        go to 20 
        end if  

c     updating Max & Min values of the stress      
        DO I=1,6
c         if ( sig(I) > sig_max(I) ) then
        ! Hidden Lines
       end do
c      print*,'sig_max**=',sig_max
c      print*,'sig_min**=',sig_min
       
       
c   ******************* failure analysis block ***********************
c   ***************************************************************************


!       if (  TIME(2) == cycletime ) then
      if ( (TIME(2)) > cycletime ) then
c   *************************************************************************************************************************************  
c   ***************** fatige failure analysis block *******************

      cycletime= cycletime + 0.01  ! do not move this to the end
      NN=2*NN
      
      
      
c      Set mean stress, amplitude stress and stress ratio (R)
c     (six individual components)      
        DO I=1,6
         sig_m(I)=0.5*( sig_max(I) + sig_min(I) )
!        ! Hidden Lines
       end do
      

c     --------------------------------------------------------
c     --------- Calculating final life (Nf)-------------------

 
      if  ( sig_m(1) > 0 ) then 
      XX=XT
      else
      XX=XC
      end if
         
      
      if  ( sig_m(2) > 0 ) then 
      YY=YT
      else
      YY=YC
      end if
      
      if  ( sig_m(3) > 0 ) then 
      ZZ=ZT
      else
      ZZ=ZC
      end if
            
       DW_star = ( sig_max(1)**2-sig_min(1)**2 )/ (XX**2) 
     1              ! Hidden Lines
     
       
      unifiedNf=( DW_star/1.088 )**! Hidden Lines
      
!       
!c     --------------------------------------------------------
!c     --------- Calculating final life (Nf)-------------------      
!    
!    
!c     -------------------------------------------------------------------------
!c     ----------- stiffness & strength degradation (5 components)--------------

c     In 11 direction:      
      if ( NN > unifiedNf ) then
          Nf_norm(1)=1.000
          call SuddenDegFiber( E_11,E_22,G_12
     1                ,nu_12,nu_23
     1                ,XT,XC,YT,YC,S12,S23 )
          
          else
          
      Nf_norm(1)= ! Hidden Lines


      E11_f =  sig_max(1)/0.0136   
      E_11= E11_f + ! Hidden Lines
     

  
         if ( sig_m(1)>=0 ) then                          
          XT= sig_max(1)+ ! Hidden Lines
         else
         XC= sig_max(1)+ ! Hidden Lines
         end if
         
      end if     
   
   
c     In 22 direction: 
      if ( NN > unifiedNf ) then
          Nf_norm(2)=1.000
          call SuddenDegMatrix (E_22,G_12,nu_12,nu_23
     1                    ,YT,YC,S12,S23  ) 

         
          
          else            
      Nf_norm(2)= ! Hidden Lines
      
      
      
      
c     stiffness degradation rule (22 direction)      
      E22_f =  sig_max(2)/0.0068   
      E_22= E22_f + (E_22 - E22_f)*                          
     1 ! Hidden Lines
       
      
c     strength degradation rule (22 direction)
      if ( sig_m(2)>=0 ) then                                 
          YT= sig_max(2)+ (YT - sig_max(2))*
     1 ! Hidden Lines
         else
         YC= sig_max(2)+ (YC - sig_max(2))*
     1 ! Hidden Lines
         end if
      
      
      
      end if   
      
      
c     In 33 direction:            
      if ( NN > unifiedNf ) then
           Nf_norm(3)=1.000 
           call SuddenDegMatrix (E_22,G_12,nu_12,nu_23
     1                    ,YT,YC,S12,S23  ) 

               
          else
      Nf_norm(3)= ! Hidden Lines
      
      
      
c     stiffness degradation rule (33 direction)
            E33_f =  sig_max(3)/0.0068   
      E_33= E33_f + ! Hidden Lines

c     strength degradation rule (33 direction)      
         if ( sig_m(3)>=0 ) then                              
          ZT= sig_max(3)+ (ZT - sig_max(3))*
     1 ! Hidden Lines
         else
         ZC= sig_max(3)+ (ZC - sig_max(3))*
     1 ! Hidden Lines
         end if
      
      
       
      end if
      
      
c     In 12 direction:     
      if ( NN > unifiedNf ) then
           Nf_norm(4)=1.000
           call SuddenDegMatrix (E_22,G_12,nu_12,nu_23
     1                    ,YT,YC,S12,S23  ) 
          
          else
      Nf_norm(4)= log(NN)/log( unifiedNf )
      
      
c     stiffness degradation rule (12 direction)      
      G12_f =  sig_max(4)/0.101   
      G_12= G12_f + ! Hidden Lines
      
      
c     strength degradation rule (12 direction)
          S12= sig_max(4)+ ! Hidden Lines

    
      end if
      
      
c     In 23 direction:      
      if ( NN > unifiedNf ) then
          Nf_norm(5)=1.000
          call SuddenDegMatrix (E_22,G_12,nu_12,nu_23
     1                    ,YT,YC,S12,S23  ) 

          
          else
      Nf_norm(5)= log(NN)/log( unifiedNf )
      
c     stiffness degradation rule (23 direction)      
c     is not required. in SUBROUTINE StiffnesMatrix: G23 = C(5,5)= ( C(2,2)-C(2,3) )/2
      
      
c     strength degradation rule (33 direction)
      S23= sig_max(5)+ (S23 - sig_max(5))*
     1 ! Hidden Lines

     
      end if
  
      
c     In 31 direction:      
      if ( NN > unifiedNf ) then
          Nf_norm(6)=1.000
          call SuddenDegMatrix (E_22,G_12,nu_12,nu_23
     1                    ,YT,YC,S12,S23  ) 

          
          else
      Nf_norm(6)= log(NN)/log( unifiedNf )
      
      
c     stiffness degradation rule (31 direction) 
      G31_f =  sig_max(6)/0.101   
      G_31= ! Hidden Lines
      
      
c     strength degradation rule (31 direction)      
      S31= ! Hidden Lines


      end if
      

c      Reset sigma-max & sigma-min for the next step        
c       sig_max(:)=0
c       sig_min(:)=0      
     
     
c     ----------- stiffness & strength degradation (5 components) --------------
c     -------------------------------------------------------------------------     
        
        
      end if  ! end of fatigue analysis     
c   ***************** fatige failure analysis block *******************
c   ******************************************************************* 
       
  20   CALL StiffnesMatrix(E_11,E_22,E_33, 
     1                    G_12,G_23, G_31,
     1                    nu_12,nu_23, CC ) 
     
     
      sig(:)=0      !Setting zero is necessary for confidence
         DO I=1,6
           do J=1,6
          sig(I) = sig(I) + 
     1                      ! Hidden Lines
           end do
         END DO

      STRESS(:)=sig(:)    !<<<      
c     setting ddsdde
         DO I=1,6
           DO J=1,6
          DDSDDE(I,J) = CC(I,J)
           END Do
         END DO     
      
c     Setting the new state variables:
c     Total strains:               
      STATEV(1) = eps(1) 
      STATEV(2) = eps(2) 
      STATEV(3) = eps(3) 
      STATEV(4) = eps(4)
      STATEV(5) = eps(5)
      STATEV(6) = eps(6)
      
c     Stiffness engineering constants:
      STATEV(7)=E_11 
      STATEV(8)=E_22
      STATEV(9)=E_33
      STATEV(10)=G_12
      STATEV(11)=G_23
      STATEV(12)=G_31 
      STATEV(13)=nu_12 
      STATEV(14)=nu_23

c     Strength parameters:
      STATEV(15) = XT
      STATEV(16) = XC 
      STATEV(17) = YT 
      STATEV(18) = YC
      STATEV(19) = ZT
      STATEV(20) = ZC 
      STATEV(21) = S12
      STATEV(22) = S23
      STATEV(23) = S31

c     Max_stress components:     
      STATEV(24) = sig_max(1)
      STATEV(25) = sig_max(2) 
      STATEV(26) = sig_max(3)
      STATEV(27) = sig_max(4) 
      STATEV(28) = sig_max(5)
      STATEV(29) = sig_max(6)
      
c     Min_stress components:      
      STATEV(30) = sig_min(1)
      STATEV(31) = sig_min(2) 
      STATEV(32) = sig_min(3)
      STATEV(33) = sig_min(4) 
      STATEV(34) = sig_min(5)
      STATEV(35) = sig_min(6) 
      
c     Mean_stress components:      
      STATEV(36) = sig_m(1)
      STATEV(37) = sig_m(2)
      STATEV(38) = sig_m(3)
      STATEV(39) = sig_m(4)
      STATEV(40) = sig_m(5)
      STATEV(41) = sig_m(6)
      
c     Stress amplitude components:      
      STATEV(42) = sig_a(1)
      STATEV(43) = sig_a(2)
      STATEV(44) = sig_a(3)
      STATEV(45) = sig_a(4)
      STATEV(46) = sig_a(5)
      STATEV(47) = sig_a(6)
 
c     Stress ratio components 
      STATEV(48) =  RR(1)
      STATEV(49) =  RR(2)
      STATEV(50) =  RR(3)
      STATEV(51) =  RR(4)
      STATEV(52) =  RR(5)
      STATEV(53) =  RR(6)
    
c     Equivalent stress components(Goodman method only)    
      STATEV(54) = sig_a_eq(1)
      STATEV(55) = sig_a_eq(2)
      STATEV(56) = sig_a_eq(3)
      STATEV(57) = sig_a_eq(4)
      STATEV(58) = sig_a_eq(5)
      STATEV(59) = sig_a_eq(6)     
      
c     Normalized life components:      
      STATEV(60) = Nf_norm(1)
      STATEV(61) = Nf_norm(2)
      STATEV(62) = Nf_norm(3)
      STATEV(63) = Nf_norm(4)
      STATEV(64) = Nf_norm(5)
      STATEV(65) = Nf_norm(6)
      
c     Final life components:      
      STATEV(66) = Nf(1)
      STATEV(67) = Nf(2)
      STATEV(68) = Nf(3)
      STATEV(69) = Nf(4)
      STATEV(70) = Nf(5)
      STATEV(71) = Nf(6)
      
c     other usefull quatities:      
      STATEV(72) = Fiber_Failure
      STATEV(73) = Matrix_Failure  
      STATEV(74) = DW_star
      STATEV(75) = unifiedNf
      STATEV(76) = FiberStaticFailure
      STATEV(77) = MatrixStaticFailure
      
c     For information      
      STATEV(79)=cycletime
      STATEV(80) = NN
     
C      print*,'sig(4)****=',sig(4)
C      print*,'sig(6)****=',sig(6)
           

      return
      end SUBROUTINE UMAT_MAT3
      
  ! ----------------------------------------------------------------------------------------------------        
       
      SUBROUTINE StiffnesMatrix(E11,E22,E33,G12,G23,G31,nu12,nu23,  C ) 
      IMPLICIT NONE 
      REAL, INTENT(IN) :: E11,E22,E33,
     1                    G12,G23,G31,nu12,nu23
      REAL, INTENT(OUT) :: C(6,6) 
      real  nu13,nu21,nu32,nu31,nu123
      
c     DEFINING STIFFNESS MATRIX C6*6 IN TERMS OF ENGINEERING CONSTANTS
      nu13=nu12       
      ! Hidden Lines
     
c     2- constructing C
      C(:,:)=0.0
      
      C(1,1)=E11*(1.0-nu23*nu32)/! Hidden Lines
      ! Hidden Lines
      C(6,6)=G31
      END SUBROUTINE StiffnesMatrix

c------------------------------------------------------
 
       SUBROUTINE SuddenDegFiber( E11,E22,G12,nu12,nu23
     1                          ,XT,XC,YT,YC,S12,S23 ) 
      IMPLICIT NONE 
      REAL, INTENT(OUT) :: E11,E22,G12,nu12,nu23
     1                          ,XT,XC,YT,YC,S12,S23 
         E11=  29 
! Hidden Lines  
     
      END SUBROUTINE SuddenDegFiber
c-------------------------------------------------------      
     
      SUBROUTINE SuddenDegMatrix(E22,G12,nu12,nu23
     1                          ,YT,YC,S12,S23  ) 
      IMPLICIT NONE 
      REAL, INTENT(OUT) :: E22,G12,nu12,nu23
     1                          ,YT,YC,S12,S23  
      
         E22=  10.7  
         ! Hidden Lines 

      END SUBROUTINE SuddenDegMatrix
c-------------------------------------------------------     
     
      SUBROUTINE FailureCriterionFiber(sig1, sig4, sig6,
     1                          XT, XC, S12, FiberFailure ) 
     
      IMPLICIT NONE
      REAL, INTENT(IN) :: sig1, sig4, sig6, XT, XC, S12
      REAL, INTENT(OUT) :: FiberFailure
      real XX
      
      if  ( sig1 > 0 ) then 
      XX=XT
      else
      XX=XC
      end if
         
      FiberFailure = (sig1**2)/(XX**2) +    
     1               ! Hidden Lines         

      END SUBROUTINE FailureCriterionFiber
      
c-------------------------------------------------------     

     
      SUBROUTINE FailureCriterionMatrix(sig2, sig3, sig4, sig5, sig6,
     1                          YT, YC, S12, S23, MatrixFailure ) 
     
      IMPLICIT NONE
      REAL, INTENT(IN) :: sig2, sig3, sig4, sig5, sig6, S12, S23, YT, YC
      REAL, INTENT(OUT) :: MatrixFailure
      real YY
      
      if  ( (sig2+sig3) > 0 ) then 
      YY=YT
      else
      YY=YC
      end if
     
      MatrixFailure =  (sig2+sig3)**2/YY**2 +
     1                 ! Hidden Lines   

      END SUBROUTINE FailureCriterionMatrix
c====================================================================	  
	  