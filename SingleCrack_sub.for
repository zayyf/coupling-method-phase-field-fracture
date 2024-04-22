! ======================================================================
! User SUBROUTINE UEL for Abaqus PHASE-Field implementation:
! All rights of reproduction or distribution in any form are reserved.
! ======================================================================
! N_ELEM : Numder of elements.
! Material properties to be given through the input file, where
! PROPS(1) = Young's Modulus
! PROPS(2) = Poisson's Ratio
! PROPS(3) = AT1 (xi = 1) and AT2 (xi = 0)
! PROPS(4) = Length Scale Parameter, Characteristic Crack Width
! PROPS(5) = Strain Energy Release Rate
! ======================================================================
! Used variables 
!
! NSVARS/3 -  solution dependent variables for the phase field element 
!             and displacement field emement
!            (1/phase, 1/engery,1/von mises stress)  
!     sdv(1)=phase
! ======================================================================
      SUBROUTINE UEL(RHS,AMATRX,SVARS,ENERGY,NDOFEL,NRHS,NSVARS,
     1     PROPS,NPROPS,COORDS,MCRD,NNODE,U,DU,V,A,JTYPE,TIME,DTIME,
     2     KSTEP,KINC,JELEM,PARAMS,NDLOAD,JDLTYP,ADLMAG,PREDEF,
     3     NPREDF,LFLAGS,MLVARX,DDLMAG,MDLOAD,PNEWDT,JPROPS,NJPROP,
     4     PERIOD)
!     ==================================================================
      INCLUDE 'ABA_PARAM.INC'
!     ==================================================================
      PARAMETER(ZERO=0.D0,ONE=1.D0,TWO=2.D0,THREE=3.D0,FOUR=4.D0,
     1 SIX=6.D0,TOLER=1.0D-6,RP25=0.25D0,HALF=0.5D0,RP1=0.1D0,
     2 RP7=0.7D0,NSTV=3,N_ELEM=12570)
!     ==================================================================
      DIMENSION RHS(MLVARX,*),AMATRX(NDOFEL,NDOFEL),
     1     SVARS(NSVARS),ENERGY(8),PROPS(*),COORDS(MCRD,NNODE),
     2     U(NDOFEL),DU(MLVARX,*),V(NDOFEL),A(NDOFEL),TIME(2),
     3     PARAMS(3),JDLTYP(MDLOAD,*),ADLMAG(MDLOAD,*),
     4     DDLMAG(MDLOAD,*),PREDEF(2,NPREDF,NNODE),LFLAGS(*),
     5     JPROPS(*)
!
       REAL*8 AINTW(8),XII(8,MCRD),XI(MCRD),dNdxi(NNODE,MCRD),  
     1 dNdx(NNODE,MCRD),VJACOB(MCRD,MCRD),VJABOBINV(MCRD,MCRD),  
     2 BB(6,NDOFEL),EPS(6),STRESS(6),CMAT(6,6),DDSDDE(6,6), 
     3 SDV(NSVARS),AN(NNODE),BP(MCRD,NNODE),DP(MCRD),VNI(MCRD,NDOFEL),
     4 ULOC(MCRD),B_B(NNODE,NNODE),EIGV(3),ALPHAI(3),VECTI(3),
     5 DCMAT(6,6),EPSC(6)
	
!
       REAL*8 DTM,EG,EG2,EG3,EBULK,EBULK3,ELAM,EMOD,ENU
      
       REAL*8 SAVG,SDIF,SDEV,SMAX,SMIN,EAVG,EDIF,EDEV,EMAX,EMIN
	 
       REAL*8 ENG0,ENG,ENGN,PHASE,PARK,CLPAR,GCPAR,PHASE0,DPHASE,
     1 HIST,HISTN,ALPHAE,THCK,VM_STRESS
	 
       REAL*8 DALPHA,DDALPHA,FT,C0,OMEGA,DOMEGA,DDOMEGA,
     1 PHASE_SOURCE,DPHASE_SOURCE,ELAMEG,ELAMEL,ENGP,ENGEG
!
       INTEGER I,J,K,L,K1,K2,K3,K4,NM,NN,NP     
      
       COMMON/KUSER/USRVAR(N_ELEM,NSTV,8)
!     History variables
        ENGKG=ZERO
        ENGEG=ZERO	 
!
!     ==================================================================
!     Material parameters
!     ==================================================================
       EMOD =PROPS(1)
       ENU  =PROPS(2)
       ATPAR=PROPS(3)
       CLPAR=PROPS(4)
       GCPAR=PROPS(5) 
       THCK =ONE
       ELAMEG=EMOD/(TWO*(ONE+ENU))
       ELAMEL=ELAMEG*TWO*ENU/(ONE-TWO*ENU) 	   
!     ==================================================================
!     Initial preparetions
!     ==================================================================
       DO K1 = 1, NDOFEL                      
        DO KRHS = 1, NRHS
         RHS(K1,KRHS) = ZERO
        END DO
        DO K2 = 1, NDOFEL
         AMATRX(K2,K1) = ZERO
        END DO
       END DO
!     ==================================================================
!     Local coordinates and weights
!     ==================================================================   	   
       IF (JTYPE.EQ.ONE) THEN
! U1: 3-node triangular elements    
        XII(1,1) = ONE/THREE
        XII(1,2) = ONE/THREE
        INNODE = 1.D0
        AINTW(1) = HALF
!
       ELSEIF (JTYPE.EQ.TWO) THEN
! U2: 4-node square elements	   
        XII(1,1) = -ONE/THREE**HALF
        XII(1,2) = -ONE/THREE**HALF
        XII(2,1) = ONE/THREE**HALF
        XII(2,2) = -ONE/THREE**HALF
        XII(3,1) = ONE/THREE**HALF
        XII(3,2) = ONE/THREE**HALF
        XII(4,1) = -ONE/THREE**HALF
        XII(4,2) = ONE/THREE**HALF
        INNODE = 4.D0
        DO I=1,INNODE
         AINTW(I) = ONE
        END DO
!		
       ELSEIF (JTYPE.EQ.THREE) THEN
! U3: 4-node linear tetrahedron elements 	
        XII(1,1) = ONE/FOUR
        XII(1,2) = ONE/FOUR
        XII(1,3) = ONE/FOUR
        INNODE=1.D0
        AINTW(1) = ONE/SIX
!
       ELSEIF (JTYPE.EQ.FOUR) THEN
! U4: 8-node linear brick elements 
        XII(1,1) = -ONE/THREE**HALF
        XII(1,2) = -ONE/THREE**HALF
        XII(1,3) = -ONE/THREE**HALF
        XII(2,1) = ONE/THREE**HALF
        XII(2,2) = -ONE/THREE**HALF
        XII(2,3) = -ONE/THREE**HALF
        XII(3,1) = ONE/THREE**HALF
        XII(3,2) = ONE/THREE**HALF
        XII(3,3) = -ONE/THREE**HALF
        XII(4,1) = -ONE/THREE**HALF
        XII(4,2) = ONE/THREE**HALF
        XII(4,3) = -ONE/THREE**HALF
        XII(5,1) = -ONE/THREE**HALF
        XII(5,2) = -ONE/THREE**HALF
        XII(5,3) = ONE/THREE**HALF
        XII(6,1) = ONE/THREE**HALF
        XII(6,2) = -ONE/THREE**HALF
        XII(6,3) = ONE/THREE**HALF
        XII(7,1) = ONE/THREE**HALF
        XII(7,2) = ONE/THREE**HALF
        XII(7,3) = ONE/THREE**HALF
        XII(8,1) = -ONE/THREE**HALF
        XII(8,2) = ONE/THREE**HALF
        XII(8,3) = ONE/THREE**HALF
        INNODE = 8.D0
        DO I=1,INNODE
         AINTW(I) = ONE
        END DO
       ENDIF	   	   
!
!     ==================================================================
!     Calculating properties at each integration point
!     ==================================================================   
       DO INPT=1,INNODE
!     Initializing solution dependent variables
        DO I=1,NSTV
           SDV(I)=ZERO
        END DO
        DO I=1,NSTV
          SDV(I)=SVARS(NSTV*(INPT-1)+I)
        END DO
!
!     Local coordinates of the integration point
        DO I=1,MCRD
         XI(I) = XII(INPT,I)
        END DO
! 
!     Shape functions and local derivatives
        IF (JTYPE.EQ.ONE) THEN
         CALL SHAPEFUNTRI(AN,dNdxi,XI)
        ELSEIF (JTYPE.EQ.TWO) THEN
         CALL SHAPEFUNQUAD(AN,dNdxi,XI)
        ELSEIF (JTYPE.EQ.THREE) THEN
         CALL SHAPEFUNTET(AN,dNdxi,XI)
        ELSEIF (JTYPE.EQ.FOUR) THEN
         CALL SHAPEFUNBRICK(AN,dNdxi,XI)
        ENDIF
! 
!     Shape functions
        IF (MCRD.EQ.2) THEN
! --------- 2D Elements ------------ 		
         IY=ZERO
         DO I = 1,NNODE
          IX=IY+1
          IY=IX+1
          VNI(1,IX)=AN(I)
          VNI(1,IY)=ZERO
          VNI(2,IX)=ZERO
          VNI(2,IY)=AN(I)
         END DO
        ELSE
! --------- 3D Elements ------------  
         IZ=ZERO
         DO I = 1,NNODE
          IX=IZ+ONE
          IY=IX+ONE
          IZ=IY+ONE
          VNI(1,IX)=AN(I)
          VNI(2,IX)=ZERO
          VNI(3,IX)=ZERO
          VNI(1,IY)=ZERO
          VNI(2,IY)=AN(I)
          VNI(3,IY)=ZERO
          VNI(1,IZ)=ZERO
          VNI(2,IZ)=ZERO
          VNI(3,IZ)=AN(I)
         END DO
        ENDIF
!		
!     Jacobian
        DO I = 1,MCRD
         DO J = 1,MCRD
          VJACOB(I,J) = ZERO
          DO K = 1,NNODE
           VJACOB(I,J) = VJACOB(I,J) + COORDS(I,K)*dNdxi(K,J)
          END DO
         END DO
        END DO
!        
        DTM = ZERO
        IF (MCRD.EQ.2) THEN
! --------- 2D Elements ------------        
         DTM = VJACOB(1,1)*VJACOB(2,2)-VJACOB(1,2)*VJACOB(2,1)
        ELSE
! --------- 3D Elements ------------        
        DTM = VJACOB(1,1)*VJACOB(2,2)*VJACOB(3,3)+VJACOB(1,2)*
     1   VJACOB(2,3)*VJACOB(3,1)+VJACOB(1,3)*VJACOB(2,1)*
     2   VJACOB(3,2)-VJACOB(3,1)*VJACOB(2,2)*VJACOB(1,3)-
     3   VJACOB(3,2)*VJACOB(2,3)*VJACOB(1,1)-VJACOB(3,3)*
     4   VJACOB(2,1)*VJACOB(1,2)        
        ENDIF 
!		
        IF (DTM.LT.ZERO) THEN
         WRITE(7,*) 'Negative Jacobian',DTM
         CALL XIT	
        ENDIF				

!     Inverse of Jacobian
        IF (MCRD.EQ.2) THEN
! --------- 2D Elements ------------        
         VJABOBINV(1,1)= VJACOB(2,2)/DTM
         VJABOBINV(1,2)=-VJACOB(1,2)/DTM
         VJABOBINV(2,1)=-VJACOB(2,1)/DTM
         VJABOBINV(2,2)= VJACOB(1,1)/DTM
        ELSE
! --------- 3D Elements ------------        
         VJABOBINV(1,1)=(VJACOB(2,2)*VJACOB(3,3)-VJACOB(2,3)*
     1    VJACOB(3,2))/DTM
         VJABOBINV(1,2)=-(VJACOB(1,2)*VJACOB(3,3)-VJACOB(3,2)*
     1    VJACOB(1,3))/DTM
         VJABOBINV(1,3)=(VJACOB(1,2)*VJACOB(2,3)-VJACOB(1,3)*
     1    VJACOB(2,2))/DTM
         VJABOBINV(2,1)=-(VJACOB(2,1)*VJACOB(3,3)-VJACOB(2,3)*
     1    VJACOB(3,1))/DTM
         VJABOBINV(2,2)=(VJACOB(1,1)*VJACOB(3,3)-VJACOB(1,3)*
     1    VJACOB(3,1))/DTM
         VJABOBINV(2,3)=-(VJACOB(1,1)*VJACOB(2,3)-VJACOB(1,3)*
     1    VJACOB(2,1))/DTM
         VJABOBINV(3,1)=(VJACOB(2,1)*VJACOB(3,2)-VJACOB(2,2)*
     1    VJACOB(3,1))/DTM
         VJABOBINV(3,2)=-(VJACOB(1,1)*VJACOB(3,2)-VJACOB(1,2)*
     1    VJACOB(3,1))/DTM
         VJABOBINV(3,3)=(VJACOB(1,1)*VJACOB(2,2)-VJACOB(1,2)*
     1    VJACOB(2,1))/DTM      
        ENDIF 
!        
!     Derivatives of shape functions respect to global ccordinates
        DO K = 1,NNODE
         DO I = 1,MCRD
          dNdx(K,I) = ZERO
          DO J = 1,MCRD
           dNdx(K,I) = dNdx(K,I) + dNdxi(K,J)*VJABOBINV(J,I)
          END DO
         END DO
        END DO
!
!     Calculating B matrix (B=LN)
       IF (MCRD.EQ.2) THEN
! --------- 2D Elements ------------ 
        IY=0
        DO INODE=1,NNODE
         IX=IY+1
         IY=IX+1
         BB(1,IX)= dNdx(INODE,1)
         BB(1,IY)= ZERO
         BB(2,IX)= ZERO
         BB(2,IY)= dNdx(INODE,2)
         BB(3,IX)= dNdx(INODE,2)
         BB(3,IY)= dNdx(INODE,1)
        END DO
!
       ELSE
! --------- 3D Elements ------------ 
        IZ=ZERO
        DO INODE=1,NNODE
         IX=IZ+ONE
         IY=IX+ONE
         IZ=IY+ONE
         BB(1,IX)= dNdx(INODE,1)
         BB(2,IX)= ZERO
         BB(3,IX)= ZERO
         BB(4,IX)= dNdx(INODE,2)
         BB(5,IX)= dNdx(INODE,3)
         BB(6,IX)= ZERO
         BB(1,IY)= ZERO
         BB(2,IY)= dNdx(INODE,2)
         BB(3,IY)= ZERO
         BB(4,IY)= dNdx(INODE,1)
         BB(5,IY)= ZERO
         BB(6,IY)= dNdx(INODE,3)
         BB(1,IZ)= ZERO
         BB(2,IZ)= ZERO
         BB(3,IZ)= dNdx(INODE,3)
         BB(4,IZ)= ZERO
         BB(5,IZ)= dNdx(INODE,1)
         BB(6,IZ)= dNdx(INODE,2)
        END DO
       ENDIF	   
!		
       DO INODE=1,NNODE
        DO K1=1,MCRD
         BP(K1,INODE)=dNdx(INODE,K1)
        END DO
       END DO	   
    
!
       DO I=1,NNODE
        DO J=1,NNODE
           B_B(I,J)=ZERO
        END DO 
       END DO	   
!      
       B_B=MATMUL(TRANSPOSE(BP),BP)
!
!     ==================================================================
!     Nodal displacements
!     ==================================================================
        DO J=1,MCRD
         ULOC(J)=ZERO
        END DO
        DO J=1,MCRD
         DO I=1,MCRD*NNODE
          ULOC(J)=ULOC(J)+VNI(J,I)*U(I)
         END DO
        END DO  
!        DO J=1,MCRD
!         SDV(J)=ULOC(J)
!        END DO
! 
!     ==================================================================
!     Nodal phase-field
!     ==================================================================		
        PHASE=ZERO
        DO I=1,NNODE
        PHASE=PHASE+AN(I)*U(I+MCRD*NNODE)
        END DO
!         
        PHASE0=SDV(1)
        IF (PHASE.LT.SDV(1)) THEN
         PHASE=SDV(1)
        ENDIF
!
        IF (PHASE.LT.ZERO) THEN
         PHASE=ZERO
        ELSEIF (PHASE.GT.ONE) THEN
         PHASE=ONE
        ENDIF
        SDV(1)=PHASE      
!      
        DPHASE=ZERO 
        DPHASE=PHASE-PHASE0 
        IF (DPHASE.LT.ZERO) THEN
         DPHASE=ZERO
        ENDIF

!       IF (PHASE.LT.RP7 .AND. DPHASE.GE.HALF) TIME(1)=RP1*DTIME

!       Gradient
        DO I=1,MCRD
         DP(I)=ZERO
        END DO
        DO I=1,MCRD
         DO J=1,NNODE
          DP(I)=DP(I)+BP(I,J)*U(J+MCRD*NNODE)
         END DO
        END DO		
!		
!     ==================================================================
!     Calculating strain
!     ================================================================== 
       IF (MCRD.EQ.2) THEN
! --------- 2D Elements ------------ 
        NP=3		
       ELSE	   
! --------- 3D Elements ------------       
        NP=6
       ENDIF

        DO J=1,NP
         EPS(J)=ZERO
        END DO
        DO I=1,NP
         DO J=1,MCRD*NNODE
           EPS(I)=EPS(I)+BB(I,J)*U(J)    
         END DO
        END DO
!		
        DO K1=1,6
         EPSC(K1)=ZERO
        END DO		
        IF (MCRD.EQ.2) THEN
! --------- 2D Elements ------------ 
         EPSC(1)=EPS(1)
         EPSC(2)=EPS(2)
         EPSC(4)=EPS(3) 		 
        ELSE	   
! --------- 3D Elements ------------       
         EPSC=EPS
        ENDIF		
!        DO J=1,BP
!         SDV(J)=EPS(J)
!        END DO
!     ==================================================================
!     Calculating principal materials stiffness matrix
!     ==================================================================
       IF (MCRD.EQ.2) THEN
! --------- 2D Elements ------------ 
        DO I=1,3
         DO J=1,3
          CMAT(I,J)=ZERO
         END DO
        END DO
        CMAT(1,1)=EMOD/((ONE+ENU)*(ONE-TWO*ENU))*(ONE-ENU)
        CMAT(2,2)=EMOD/((ONE+ENU)*(ONE-TWO*ENU))*(ONE-ENU)
        CMAT(3,3)=EMOD/((ONE+ENU)*(ONE-TWO*ENU))*(HALF-ENU)
        CMAT(1,2)=EMOD/((ONE+ENU)*(ONE-TWO*ENU))*ENU
        CMAT(2,1)=EMOD/((ONE+ENU)*(ONE-TWO*ENU))*ENU
!		
       ELSE	   
! --------- 3D Elements ------------       
        DO I=1,6
         DO J=1,6
          CMAT(I,J)=ZERO
         END DO
        END DO
        CMAT(1,1)=EMOD/((ONE+ENU)*(ONE-TWO*ENU))*(ONE-ENU)
        CMAT(2,2)=EMOD/((ONE+ENU)*(ONE-TWO*ENU))*(ONE-ENU)
        CMAT(3,3)=EMOD/((ONE+ENU)*(ONE-TWO*ENU))*(ONE-ENU)
        CMAT(1,2)=EMOD/((ONE+ENU)*(ONE-TWO*ENU))*ENU
        CMAT(2,1)=EMOD/((ONE+ENU)*(ONE-TWO*ENU))*ENU
        CMAT(1,3)=EMOD/((ONE+ENU)*(ONE-TWO*ENU))*ENU
        CMAT(3,1)=EMOD/((ONE+ENU)*(ONE-TWO*ENU))*ENU
        CMAT(2,3)=EMOD/((ONE+ENU)*(ONE-TWO*ENU))*ENU
        CMAT(3,2)=EMOD/((ONE+ENU)*(ONE-TWO*ENU))*ENU
        CMAT(4,4)=EMOD/((ONE+ENU)*(ONE-TWO*ENU))*(HALF-ENU)
        CMAT(5,5)=EMOD/((ONE+ENU)*(ONE-TWO*ENU))*(HALF-ENU)
        CMAT(6,6)=EMOD/((ONE+ENU)*(ONE-TWO*ENU))*(HALF-ENU)
       ENDIF
!
!     ==================================================================
!     Calculating stresses
!     ==================================================================

        CALL GEOMETRICFUNC(PROPS,DALPHA,DDALPHA,PHASE,FT,C0) ! GEOMETRIC FUNCTION
        CALL ENERGETICFUNC(OMEGA,DOMEGA,DDOMEGA,PHASE) ! ENERGETIC FUNCTION
	   
        DO K1=1,NP
         STRESS(K1)=ZERO
        END DO
!
        DO K1=1,NP
         DO K2=1,NP
          STRESS(K2)=STRESS(K2)+CMAT(K2,K1)*EPS(K1)
         END DO
        END DO

!      Calculating von Mises Stress		
        IF (MCRD.EQ.2) THEN
! --------- 2D Elements ------------ 
        VM_STRESS=(HALF*((STRESS(1)-STRESS(2))**TWO+STRESS(1)**TWO+ 
     1  STRESS(2)**TWO+SIX*(STRESS(3)**TWO)))**HALF	 
!		
        ELSE	   
! --------- 3D Elements ------------  
        VM_STRESS=(HALF*((STRESS(1)-STRESS(2))**TWO+
     1  (STRESS(2)-STRESS(3))**TWO+(STRESS(3)-STRESS(1))**TWO+ 
     2  SIX*(STRESS(4)**TWO+STRESS(5)**TWO+STRESS(6)**TWO)))**HALF	 
        ENDIF
		
        SDV(3)=VM_STRESS		
!		
!        DO J=1,NP
!         SDV(J+6)=STRESS(J)*OMEGA
!        END DO
! 
!     ==================================================================
!     Calculating elastic ENERGY
!     ==================================================================

        ENG0=HALF*FT**TWO/EMOD

        CALL EIGOWN(EPSC,EIGV,ALPHAE,ALPHAI,VECTI)
!       
        ENG=(ELAMEL*(ALPHAE*(EIGV(1)+EIGV(2)+EIGV(3)))**TWO)/
     1       TWO+ELAMEG*((EIGV(1)*ALPHAI(1))**TWO+(EIGV(2)*
     2       ALPHAI(2))**TWO+(EIGV(3)*ALPHAI(3))**TWO)
		

        ENGE=ZERO
        DO K2=1,NP
         ENGE=ENGE+STRESS(K2)*EPS(K2)*HALF
        END DO
		
!        WRITE(7,*) 'EIGV',EIGV	

        ENGN=max(ENG,ENG0)
		
        ENGEG=ENGEG+(ENGE*OMEGA)*DTM*THCK  

        ENERGY(2)=ENGEG
!
!     ==================================================================
!     Calculating elastic ENERGY history
!     ==================================================================         
!
        HISTN=SDV(2)
        IF (ENGN.GT.HISTN) THEN
         HIST=ENGN
        ELSE
         HIST=HISTN
        ENDIF
        SDV(2)=HIST
!
!     ==================================================================
!     Calculating element stiffness matrix(K_uu)
!     ==================================================================
!
        DO K=1,MCRD*NNODE
         DO L=1,MCRD*NNODE
          DO I=1,NP
           DO J=1,NP
            AMATRX(K,L)=AMATRX(K,L)+AINTW(INPT)*BB(I,K)*CMAT(I,J)*
     1                  BB(J,L)*DTM*OMEGA*THCK
           END DO
          END DO
         END DO
        END DO
!       
!     ==================================================================
!     Internal forces (residual vector) (r_u)
!     ==================================================================

        DO K1=1,MCRD*NNODE
         DO K4=1,NP
           RHS(K1,1)=RHS(K1,1)-AINTW(INPT)*BB(K4,K1)*STRESS(K4)*DTM*
     1               OMEGA*THCK
         END DO
        END DO
!  
!     ==================================================================
!     Calculating element stiffness matrix (K_dd)
!     ==================================================================
 
        PHASE_SOURCE=DOMEGA*HIST+GCPAR/(FOUR*C0*CLPAR)*DALPHA
        DPHASE_SOURCE=DDOMEGA*HIST+GCPAR/(FOUR*C0*CLPAR)*DDALPHA
!		
        NM=MCRD*NNODE
        DO I=1,NNODE
        NM=NM+1
        NN=MCRD*NNODE
         DO J=1,NNODE
            NN=NN+1
            AMATRX(NM,NN)=AMATRX(NM,NN)+GCPAR*CLPAR/(TWO*C0)*
     1      B_B(I,J)*AINTW(INPT)*DTM+AN(I)*AN(J)*AINTW(INPT)*
     2      DTM*(DPHASE_SOURCE)*THCK
         END DO
        END DO      
!            
!     ==================================================================
!     Internal forces (residual vector) (r_d)
!     ==================================================================

        NM=MCRD*NNODE
        DO I=1,NNODE
        NM=NM+1
         DO J=1,MCRD
           RHS(NM,1)=RHS(NM,1)-BP(J,I)*DP(J)*GCPAR*CLPAR/(TWO*C0)*
     1     AINTW(INPT)*DTM*THCK
         END DO
         RHS(NM,1)=RHS(NM,1)-AN(I)*AINTW(INPT)*DTM*(PHASE_SOURCE)*THCK  
        END DO		
!       
!     ==================================================================
!     Uploading solution dep. variables
!     ==================================================================    
        DO I=1,NSTV
          SVARS(NSTV*(INPT-1)+I)=SDV(I)
          USRVAR(JELEM,I,INPT)=SDV(I)
        END DO
       END DO
!
      RETURN
      END
!
!----------------------------------------------------------------------
!
      SUBROUTINE ENERGETICFUNC(OMEGA,DOMEGA,DDOMEGA,PHASE)
!
! ----------------------------------------------------------------------
      INCLUDE 'ABA_PARAM.INC'
!      
      REAL*8 OMEGA,DOMEGA,DDOMEGA,PHASE,PARK
      PARAMETER(ZERO=0.D0,ONE=1.D0,TWO=2.D0,THREE=3.D0,FOUR=4.D0)
!	  
      PARK=1.D-6    ! Parameter of well-Conditioned System
      
      OMEGA   =  (ONE-PHASE)**TWO+PARK        
      DOMEGA  = -TWO*(ONE-PHASE)
      DDOMEGA =  TWO
      
      RETURN
      END
      
! ----------------------------------------------------------------------
!
      SUBROUTINE GEOMETRICFUNC(PROPS,DALPHA,DDALPHA,PHASE,FT,C0)
!
! ----------------------------------------------------------------------
      INCLUDE 'ABA_PARAM.INC'
!     
      DIMENSION PROPS(*)
      REAL*8 DALPHA,DDALPHA,PHASE,FT,C0,EMOD,ATPAR,CLPAR,GCPAR
      PARAMETER(ZERO=0.D0,ONE=1.D0,TWO=2.D0,THREE=3.D0,
     1 RP75=0.375D0,HALF=0.5D0)	  
!
      EMOD=PROPS(1)   ! PROPS(1) -- YOUNG'S MODULUS
      ATPAR=PROPS(3)  ! PROPS(3) -- AT1(XI=1) AND AT2(XI=0)
      CLPAR=PROPS(4)  ! PROPS(4) -- LENGTH SCALE
      GCPAR=PROPS(5)  ! PROPS(5) -- FRACTURE ENERGY 
	  
      DALPHA =ATPAR+TWO*(ONE-ATPAR)*PHASE
      DDALPHA=TWO*(ONE-ATPAR)
	  
!     THE CRITICAL STRENGTH UPON CRACK NUCLEATION (NOT
!     NECESSARILY THE PEAK STRENGTH)	  
      IF (NINT(ATPAR)==0) THEN ! AT2
         FT = ZERO
         C0 = HALF
      ELSEIF (NINT(ATPAR)==1) THEN ! AT1
         FT = SQRT(RP75*EMOD*GCPAR/CLPAR)
         C0 = TWO/THREE
      ELSE
         WRITE (*,*) '**ERROR: MODEL LAW NO. ', ATPAR, 
     1                 'DOES NOT EXIST!'
         CALL XIT
      ENDIF           

      RETURN 
      END 	  
!
! ==============================================================
! Eigenstrains from Voigt notation
! ==============================================================
!
      SUBROUTINE EIGOWN(EPS,EIGV,ALPHAE,ALPHAI,VECTI)
      INCLUDE 'ABA_PARAM.INC'
      PARAMETER(ZERO=0.D0,ONE=1.D0,MONE=-1.D0,TWO=2.D0,
     1  TS=27.D0,THREE=3.D0,HALF=0.5D0,TOLER=1.0D-12,FOUR=4.D0,
     2  CNTN=100,TOLERE=1.0D-12)
      INTEGER I, J, K
      REAL*8 EPS(6), EIGV(3), ALPHAI(3), VECTI(3)
      REAL*8 PC, QC, ALPHAE, DISC, PI, CNT
!
       PI=FOUR*ATAN(ONE)
!   Scaling the strain vector
       VMAXE=MAXVAL(ABS(EPS))
       IF (VMAXE.GT.TOLERE) THEN
        DO K1=1,6
         EPS(K1)=EPS(K1)/VMAXE
        END DO
       ENDIF
!    
!   Calculating eigenvalues
       VECTI(1)=EPS(1)+EPS(2)+EPS(3)
       VECTI(2)=EPS(1)*EPS(2)+EPS(2)*EPS(3)+EPS(1)*EPS(3)-
     1 EPS(4)**TWO/FOUR-EPS(5)**TWO/FOUR-EPS(6)**TWO/FOUR
       VECTI(3)=EPS(1)*EPS(2)*EPS(3)+EPS(4)*EPS(5)*EPS(6)/
     1 FOUR-EPS(1)*EPS(6)**TWO/FOUR-EPS(2)*EPS(5)**TWO/FOUR-
     2 EPS(3)*EPS(4)**TWO/FOUR
!
!   Depressed coefficients    
       PC=VECTI(2)-VECTI(1)**TWO/THREE
       QC=VECTI(1)*VECTI(2)/THREE-TWO*VECTI(1)**THREE/TS-VECTI(3)
       DISC=MONE*(FOUR*PC**THREE+TS*QC**TWO)
!
       DO I=1,3
        EIGV(I)=ZERO
       END DO
       CNT=ZERO
       IF (ABS(DISC).LT.TOLER) THEN
        IF ((ABS(QC).LT.TOLER).AND.(ABS(PC).LT.TOLER)) THEN
         EIGV(1)=VECTI(1)/THREE
         EIGV(2)=VECTI(1)/THREE
         EIGV(3)=VECTI(1)/THREE
        ELSE
         EIGV(1)=-THREE*QC/TWO/PC+VECTI(1)/THREE
         EIGV(2)=-THREE*QC/TWO/PC+VECTI(1)/THREE
         EIGV(3)=THREE*QC/PC+VECTI(1)/THREE
         IF (EIGV(1).GT.EIGV(3)) THEN
          EONE=EIGV(1)
          EIGV(1)=EIGV(3)
          EIGV(3)=EONE
         ENDIF
        ENDIF
       ELSE
        DO I=1,3
         EIGV(I)=VECTI(1)/THREE+TWO*(MONE*PC/THREE)**HALF*
     1   COS(ONE/THREE*ACOS(MONE*QC/TWO*(TS/(MONE*PC**THREE))**
     2   HALF)+TWO*I*PI/THREE)
        END DO
       ENDIF
!       
       ALPHAE=ZERO
       IF ((EIGV(1)+EIGV(2)+EIGV(3)).GT.TOLER) THEN
        ALPHAE=ONE
       ENDIF
       DO K1=1,3
        ALPHAI(K1)=ZERO
        IF (EIGV(K1).GT.TOLER) THEN
         ALPHAI(K1)=ONE
        ENDIF
       END DO
!
!    Rescaling eigenvalues       
       IF (VMAXE.GT.TOLERE) THEN
        DO K1=1,6
         EPS(K1)=EPS(K1)*VMAXE
        END DO
        DO K1=1,3
         EIGV(K1)=EIGV(K1)*VMAXE
        END DO
         VECTI(1)=EIGV(1)+EIGV(2)+EIGV(3)
         VECTI(2)=EIGV(1)*EIGV(2)+EIGV(1)*EIGV(3)+EIGV(3)*EIGV(2)
         VECTI(3)=EIGV(1)*EIGV(2)*EIGV(3)
       ENDIF
!   
       RETURN
       END
!
! ======================================================================
! U1: 3-node linear elements
! U2: 4-node bilinear elements
! ======================================================================
!	 	          
! ----------------------------------------------------------------------      
! Shape functions for U1: 3-node linear elements
! ----------------------------------------------------------------------      
      SUBROUTINE SHAPEFUNTRI(AN,dNdxi,XI)
      INCLUDE 'ABA_PARAM.INC'
      Real*8 AN(4),dNdxi(3,2)
      Real*8 XI(2)
      PARAMETER(ZERO=0.D0,ONE=1.D0,MONE=-1.D0,FOUR=4.D0)

!     Values of shape functions as a function of local coord.
      AN(1) = XI(1)
      AN(2) = XI(2)
      AN(3) = ONE-XI(1)-XI(2)
!
!     Derivatives of shape functions respect to local ccordinates
      DO I=1,4
        DO J=1,2
            dNdxi(I,J) =  ZERO
        END DO
      END DO
      dNdxi(1,1) =  ONE
      dNdxi(1,2) =  ZERO
      dNdxi(2,1) =  ZERO
      dNdxi(2,2) =  ONE
      dNdxi(3,1) =  MONE
      dNdxi(3,2) =  MONE
      RETURN
      END
!
! -----------------------------------------------------------------------      
! Shape functions for U2: 4-node bilinear elements
! -----------------------------------------------------------------------      
      SUBROUTINE SHAPEFUNQUAD(AN,dNdxi,xi)
      INCLUDE 'ABA_PARAM.INC'
      Real*8 AN(4),dNdxi(4,2)
      Real*8 XI(2)
      PARAMETER(ZERO=0.D0,ONE=1.D0,MONE=-1.D0,FOUR=4.D0)
!
!     Values of shape functions as a function of local coord.
      AN(1) = ONE/FOUR*(ONE-XI(1))*(ONE-XI(2))
      AN(2) = ONE/FOUR*(ONE+XI(1))*(ONE-XI(2))
      AN(3) = ONE/FOUR*(ONE+XI(1))*(ONE+XI(2))
      AN(4) = ONE/FOUR*(ONE-XI(1))*(ONE+XI(2))
!
!     Derivatives of shape functions respect to local coordinates
      DO I=1,4
        DO J=1,2
            dNdxi(I,J) =  ZERO
        END DO
      END DO
      dNdxi(1,1) =  MONE/FOUR*(ONE-XI(2))
      dNdxi(1,2) =  MONE/FOUR*(ONE-XI(1))
      dNdxi(2,1) =  ONE/FOUR*(ONE-XI(2))
      dNdxi(2,2) =  MONE/FOUR*(ONE+XI(1))
      dNdxi(3,1) =  ONE/FOUR*(ONE+XI(2))
      dNdxi(3,2) =  ONE/FOUR*(ONE+XI(1))
      dNdxi(4,1) =  MONE/FOUR*(ONE+XI(2))
      dNdxi(4,2) =  ONE/FOUR*(ONE-XI(1))
      RETURN
      END
!
! ======================================================================
! U3: 4-node linear tetrahedron elements
! U4: 8-node linear brick elements
! ======================================================================
! 
! ----------------------------------------------------------------------      
! Shape functions for U3: 4-node linear tetrahedron elements
! ----------------------------------------------------------------------  
      SUBROUTINE SHAPEFUNTET(AN,dNdxi,XI)
      INCLUDE 'ABA_PARAM.INC'
      Real*8 AN(4),dNdxi(4,3)
      Real*8 XI(3)
      PARAMETER(ZERO=0.D0,ONE=1.D0,MONE=-1.D0,FOUR=4.D0)

!     Values of shape functions as a function of local coord.
      AN(1) = ONE-XI(1)-XI(2)-XI(3)
      AN(2) = XI(1)
      AN(3) = XI(2)
      AN(4) = XI(3)
!
!     Derivatives of shape functions respect to local coordinates
      DO I=1,4
        DO J=1,3
            dNdxi(I,J) =  ZERO
        END DO
      END DO
      dNdxi(1,1) =  MONE
      dNdxi(1,2) =  MONE
      dNdxi(1,3) =  MONE
!
      dNdxi(2,1) =  ONE
      dNdxi(2,2) =  ZERO
      dNdxi(2,3) =  ZERO
!
      dNdxi(3,1) =  ZERO
      dNdxi(3,2) =  ONE
      dNdxi(3,3) =  ZERO
!
      dNdxi(4,1) =  ZERO
      dNdxi(4,2) =  ZERO
      dNdxi(4,3) =  ONE
!
      RETURN
      END
!
! --------------------------------------------------------------------      
! Shape functions for U4: 8-node linear brick elements
! --------------------------------------------------------------------     
      SUBROUTINE SHAPEFUNBRICK(AN,dNdxi,XI)
      INCLUDE 'ABA_PARAM.INC'
      REAL*8 AN(8),dNdxi(8,3)
      REAL*8 XI(3)
      PARAMETER(ZERO=0.D0,ONE=1.D0,MONE=-1.D0,FOUR=4.D0,EIGHT=8.D0)

!     Values of shape functions as a function of local coord.
      AN(1) = ONE/EIGHT*(ONE-XI(1))*(ONE-XI(2))*(ONE-XI(3))
      AN(2) = ONE/EIGHT*(ONE+XI(1))*(ONE-XI(2))*(ONE-XI(3))
      AN(3) = ONE/EIGHT*(ONE+XI(1))*(ONE+XI(2))*(ONE-XI(3))
      AN(4) = ONE/EIGHT*(ONE-XI(1))*(ONE+XI(2))*(ONE-XI(3))
      AN(5) = ONE/EIGHT*(ONE-XI(1))*(ONE-XI(2))*(ONE+XI(3))
      AN(6) = ONE/EIGHT*(ONE+XI(1))*(ONE-XI(2))*(ONE+XI(3))
      AN(7) = ONE/EIGHT*(ONE+XI(1))*(ONE+XI(2))*(ONE+XI(3))
      AN(8) = ONE/EIGHT*(ONE-XI(1))*(ONE+XI(2))*(ONE+XI(3))
      
!     Derivatives of shape functions respect to local coordinates
      DO I=1,8
        DO J=1,3
            dNdxi(I,J) =  ZERO
        END DO
      END DO
      dNdxi(1,1) =  MONE/EIGHT*(ONE-XI(2))*(ONE-XI(3))
      dNdxi(1,2) =  MONE/EIGHT*(ONE-XI(1))*(ONE-XI(3))
      dNdxi(1,3) =  MONE/EIGHT*(ONE-XI(1))*(ONE-XI(2))
      dNdxi(2,1) =  ONE/EIGHT*(ONE-XI(2))*(ONE-XI(3))
      dNdxi(2,2) =  MONE/EIGHT*(ONE+XI(1))*(ONE-XI(3))
      dNdxi(2,3) =  MONE/EIGHT*(ONE+XI(1))*(ONE-XI(2))
      dNdxi(3,1) =  ONE/EIGHT*(ONE+XI(2))*(ONE-XI(3))
      dNdxi(3,2) =  ONE/EIGHT*(ONE+XI(1))*(ONE-XI(3))
      dNdxi(3,3) =  MONE/EIGHT*(ONE+XI(1))*(ONE+XI(2))
      dNdxi(4,1) =  MONE/EIGHT*(ONE+XI(2))*(ONE-XI(3))
      dNdxi(4,2) =  ONE/EIGHT*(ONE-XI(1))*(ONE-XI(3))
      dNdxi(4,3) =  MONE/EIGHT*(ONE-XI(1))*(ONE+XI(2))
      dNdxi(5,1) =  MONE/EIGHT*(ONE-XI(2))*(ONE+XI(3))
      dNdxi(5,2) =  MONE/EIGHT*(ONE-XI(1))*(ONE+XI(3))
      dNdxi(5,3) =  ONE/EIGHT*(ONE-XI(1))*(ONE-XI(2))
      dNdxi(6,1) =  ONE/EIGHT*(ONE-XI(2))*(ONE+XI(3))
      dNdxi(6,2) =  MONE/EIGHT*(ONE+XI(1))*(ONE+XI(3))
      dNdxi(6,3) =  ONE/EIGHT*(ONE+XI(1))*(ONE-XI(2))
      dNdxi(7,1) =  ONE/EIGHT*(ONE+XI(2))*(ONE+XI(3))
      dNdxi(7,2) =  ONE/EIGHT*(ONE+XI(1))*(ONE+XI(3))
      dNdxi(7,3) =  ONE/EIGHT*(ONE+XI(1))*(ONE+XI(2))
      dNdxi(8,1) =  MONE/EIGHT*(ONE+XI(2))*(ONE+XI(3))
      dNdxi(8,2) =  ONE/EIGHT*(ONE-XI(1))*(ONE+XI(3))
      dNdxi(8,3) =  ONE/EIGHT*(ONE-XI(1))*(ONE+XI(2))
!      
      RETURN
      END
! 
! ======================================================================
! !!! NOTE: N_ELEM has to be changed according to the UEL !!!!!
! ======================================================================
!
      
       SUBROUTINE UMAT(STRESS, STATEV, DDSDDE, SSE, SPD, SCD, RPL,
     1 DDSDDT, DRPLDE, DRPLDT, STRAN, DSTRAN, TIME, DTIME, TEMP, DTEMP,
     2 PREDEF, DPRED, CMNAME, NDI, NSHR, NTENS, NSTATV, PROPS, NPROPS,
     3 COORDS, DROT, PNEWDT, CELENT, DFGRD0, DFGRD1, NOEL, NPT, LAYER,
     4 KSPT, KSTEP, KINC)

       INCLUDE 'ABA_PARAM.INC'

       CHARACTER*80 CMNAME

       DIMENSION STRESS(NTENS), STATEV(NSTATV), DDSDDE(NTENS, NTENS),
     1 DDSDDT(NTENS), DRPLDE(NTENS), STRAN(NTENS), DSTRAN(NTENS),
     2 PREDEF(1), DPRED(1), PROPS(NPROPS), COORDS(3), DROT(3,3),
     3 DFGRD0(3, 3),DFGRD1(3,3)

       PARAMETER (ONE=1.D0,TWO=2.D0,THREE=3.D0,SIX=6.D0,HALF=0.5D0,
     1 N_ELEM=12570,NSTV=3) 
       DATA NEWTON,TOLER/40,1.D-6/ 
!       
       COMMON/KUSER/USRVAR(N_ELEM,NSTV,8)
! 
! ---------------------------------------------------------------------- 
!
!  	   Stiffness tensor
       DDSDDE=0.D0
!	   
       NELEMAN=NOEL-N_ELEM
       IF (NPT.EQ.3) THEN
        NPT=4
       ELSEIF (NPT.EQ.4) THEN
        NPT=3
       ENDIF
!	   
       IF (NPT.EQ.7) THEN
        NPT=8
       ELSEIF (NPT.EQ.8) THEN
        NPT=7
       ENDIF
!      
       DO I=1,NSTATV
        STATEV(I)=USRVAR(NELEMAN,I,NPT)
       END DO
!       
       RETURN
       END 


	   