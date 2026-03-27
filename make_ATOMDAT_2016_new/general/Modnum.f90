!==============================================================================
!             *****  Modnum contains numerical MODULES  *****
!==============================================================================
!  Modnum.f90
!
!  Modules containing   numerical routines - 03.01/3.05/8.06 ADI
!                     + special routines   - 10.04 ADI
!
! Content of Modnum (this file):
! ==============================
!  M_Special_Fpe_Cpu                , ONLY: solely ONLY  -> maintain
!  M_Simpson_Integration_weights    , ONLY: solely ONLY  -> maintain
!  M_Angular_Integration_weights    , ONLY: solely ONLY  -> maintain
!  M_Frequency_Integration          , ONLY: solely ONLY  -> maintain
!  M_Inter_Extrapolation            , ONLY: solely ONLY  -> maintain
!  M_Solve_lin_eq_systems           , ONLY: solely ONLY  -> maintain
!  M_Gaunt_factor_freefree          , ONLY: solely ONLY  -> maintain
!  M_Mathematical_functions         , ONLY: solely ONLY  -> maintain
!  M_Physical_functions             , ONLY: solely ONLY  -> maintain
!  M_3d_min_healpix                 , ONLY: solely ONLY  -> maintain
!
!==============================================================================
!                 ***** MODULE M_Special_Fpe_Cpu *****
!==============================================================================
module M_Special_Fpe_Cpu
!--- routines representing special functions
!     DigComInt
!     CPUTIME
!     BFSTOP
!--- 3.05 - ADI
!
   USE M_type_definitions
!
   IMPLICIT NONE
!
CONTAINS
!------------------------------------------------------------------------------
  subroutine DigComInt
!
!--- Special Digital-Compaq-Intel routines are called.

!--- regarding floating point exceptions:
!     bits corresponding to exceptions to be allowed are set (1),
!     bits corresponding to exceptions to be trapped are cleared (0)
!
!DEC$ IF DEFINED(_WIN32)
   USE DFLIB, ONLY: SETCONTROLFPQQ
!DEC$ ENDIF
!
      INTEGER(2), PARAMETER :: FPCW_MCW_EM     = Z'003F' ! EXCEPTION MASK
      INTEGER(2), PARAMETER :: FPCW_INVALID    = Z'0001' ! INVALID
      INTEGER(2), PARAMETER :: FPCW_DENORMAL   = Z'0002' ! DENORMAL
      INTEGER(2), PARAMETER :: FPCW_ZERODIVIDE = Z'0004' ! ZERO DIVIDE
      INTEGER(2), PARAMETER :: FPCW_OVERFLOW   = Z'0008' ! OVERFLOW
      INTEGER(2), PARAMETER :: FPCW_UNDERFLOW  = Z'0010' ! UNDERFLOW
      INTEGER(2), PARAMETER :: FPCW_INEXACT    = Z'0020' ! INEXACT (PRECISION)
!
      INTEGER(2), PARAMETER :: FPCW_MCW_PC     = Z'0300' ! PRECISION CONTROL MASK
      INTEGER(2), PARAMETER :: FPCW_64         = Z'0300' ! 64 BITS
      INTEGER(2), PARAMETER :: FPCW_53         = Z'0200' ! 53 BITS
      INTEGER(2), PARAMETER :: FPCW_24         = Z'0000' ! 24 BITS
!
      INTEGER(2), PARAMETER :: FPCW_MCW_RC     = Z'0C00' ! ROUNDING CONTROL MASK
      INTEGER(2), PARAMETER :: FPCW_CHOP       = Z'0C00' ! CHOP
      INTEGER(2), PARAMETER :: FPCW_UP         = Z'0800' ! UP
      INTEGER(2), PARAMETER :: FPCW_DOWN       = Z'0400' ! DOWN
      INTEGER(2), PARAMETER :: FPCW_NEAR       = Z'0000' ! NEAR
!
      INTEGER(2), PARAMETER :: FPCW_MCW_IC     = Z'1000' ! INFINITY CONTROL MASK
      INTEGER(2), PARAMETER :: FPCW_AFFINE     = Z'1000' ! AFFINE
      INTEGER(2), PARAMETER :: FPCW_PROJECTIVE = Z'0000' ! PROJECTIVE
!-------------------------------------------------------------------------
      CALL SETCONTROLFPQQ ( &
        IOR(FPCW_DENORMAL, &
        IOR(FPCW_UNDERFLOW, &
        IOR(FPCW_INEXACT, &
        IOR(FPCW_64, &
        IOR(FPCW_NEAR, &
            FPCW_AFFINE &
        ))))))

      RETURN
!-------------------------------------------------------------------------
  end subroutine DigComInt
!------------------------------------------------------------------------------
  subroutine CPUTIME
!
!!DIR$ INLINE ALWAYS CPUTIME
!
!--- CALLS AND PRINTS SYSTEM CPU-TIME
!
    REAL(8), SAVE       :: ELAPSED_TIME, START_TIME
    LOGICAL, SAVE       :: START=.TRUE.
!--------------------------------------------------------------------
!
!   ELAPSED_TIME = TIMEF( )
    CALL CPU_TIME ( ELAPSED_TIME )
     IF (START)THEN
      START_TIME = ELAPSED_TIME
      START = .FALSE.
     ELSE
      WRITE(*,'("                 ******************************")')
      WRITE(*,'("                 CPU-TIME: ",F9.1," SEC")') &
                                  ELAPSED_TIME-START_TIME
      WRITE(*,'("                 ******************************")')
     ENDIF
!
    RETURN
!--------------------------------------------------------------------
  end subroutine CPUTIME
!------------------------------------------------------------------------------
  subroutine BFSTOP (K,RIJ,RJI,XRAD)
!
!!DIR$ INLINE ALWAYS BFSTOP
!
!--- A core dump is provoked by causing a floating point exception
!--- CALLED BY 'RJIC' ETC.
!
!input
      INTEGER(i4b)  :: K
      REAL   (dp )  :: RIJ,RJI
      REAL   (dp ), DIMENSION (3) :: XRAD
!local
      INTEGER(i4b)  :: I
      REAL   (dp )  :: DUMP
!-----------------------------------------------------------------------
      PRINT*,'STOP ',K,': RIJ ODER RJI NEGATIV !'
      PRINT*,'RIJ=',RIJ,'   RJI=',RJI
      PRINT*,'XRAD=',XRAD
!
!-----CALL POST MORTEM DUMP --- ahemmm... doesn't work on the Cray T90
!
!     DO 10 I=10000,100000
!        XI = FLOAT(I)
!        DUMP = 10.**XI
!        PRINT*,'DUMP=',DUMP
!  10 CONTINUE

      ! Try to provoke a core dump by causing a floating point exception
      DO I=-1,1
         DUMP = 1._dp/REAL(I,dp)
         PRINT*,'DUMP=',DUMP
      END DO
!
      STOP 'RIJ ODER RJI NEGATIV !'
!
   RETURN
!-----------------------------------------------------------------------
  end subroutine BFSTOP
!------------------------------------------------------------------------------
!
end module M_Special_Fpe_Cpu
!==============================================================================
!                 ***** END MODULE M_Special_Fpe_Cpu *****
!==============================================================================
!
!==============================================================================
!             ***** MODULE M_Simpson_Integration_weights *****
!==============================================================================
Module M_Simpson_Integration_weights
!--- Simpson integration weights
!--- used by MUEINT and FORMAL
!--- resources, routines and functions
!     GWEIGHT
!
!--- 8.06 - ADI
!
   USE M_type_definitions
   USE M_global_settings               , ONLY:
!
   IMPLICIT NONE
!output
  ! 9-POINT SIMPSON-INTEGRATION WEIGHTS
    REAL   (dp ), DIMENSION (9), SAVE   :: G
!
CONTAINS
!------------------------------------------------------------------------------
  subroutine GWEIGHT
!
! CALCULATES WEIGHTS FOR 9-POINT SIMPSON-INTEGRATION
! M ADI 11.95/3.05 FILE3
!
!-------------------------------------------------------------------------
      G(1)=(1._dp-1._dp/15._dp)/3._dp
      G(2)=(1._dp+1._dp/15._dp)*4._dp/3._dp
      G(3)=1.6_dp/3._dp
      G(4)=G(2)
      G(5)=(1._dp-1._dp/15._dp)*2._dp/3._dp
      G(6)=G(2)
      G(7)=G(3)
      G(8)=G(2)
      G(9)=G(1)
!
     RETURN
!-------------------------------------------------------------------------
  end subroutine GWEIGHT
!------------------------------------------------------------------------------
!
end Module M_Simpson_Integration_weights
!==============================================================================
!           ***** END MODULE M_Simpson_Integration_weights *****
!==============================================================================
!
!==============================================================================
!             ***** MODULE M_Angular_Integration_weights *****
!==============================================================================
Module M_Angular_Integration_weights
!--- angular integration weights
!--- used by CONT_SRT
!--- resources, routines and functions
!     WEIGHT11
!     WEIGHPN
!     WEIGH1
!
!--- M 10.01 - DNS      2.05/6.06/12.14 - ADI
!
   USE M_type_definitions
   USE M_global_settings               , ONLY: NUMNT,NCORE
!
   IMPLICIT NONE
!output
  !***  0. AND 2. MOMENT -J,K.
  ! THE INTEGRATION IS PERFORMED IN THE Z VARIABLE
    REAL(dp) ,DIMENSION(:,:) ,ALLOCATABLE        :: VP,VP2
  !***  1. AND 3. MOMENT -H,N.
  ! THE INTEGRATION IS PERFORMED IN THE P VARIABLE
    REAL(dp) ,DIMENSION(:,:) ,ALLOCATABLE        :: WP,WP2
  !***  WEIGHTS FOR INTEGRAL mu*u dmu.
  ! EVALUATED ON THE REGULAR GRID POINTS - THE Z VARIABLE
    REAL(dp) ,DIMENSION(:,:) ,ALLOCATABLE        :: WX1
  !***  INTEGRATION WEIGHT FOR P = RMAX AND Z=0
    REAL(dp)                                     :: VPLAST
  !***  INTERPOLATION WEIGHTS
    REAL(dp) ,DIMENSION(NUMNT-2)                 :: Q1,Q2,Q12,Q01,   &
                                                    Q4,Q21,Q02,Q3
!
    INTEGER(i4b) ,DIMENSION(NUMNT) ,SAVE         :: JPI
    INTEGER(i4b) ,DIMENSION(:    ) ,ALLOCATABLE  :: LPI
!
    REAL(dp) ,DIMENSION(:,:) ,ALLOCATABLE        :: U1,U2
    REAL(dp) ,DIMENSION(:,:) ,ALLOCATABLE        :: TF
!
CONTAINS
!------------------------------------------------------------------------------
  subroutine Allocate_Mangintw (O_other)
!input
      LOGICAL       :: O_other
!
!-------------------------------------------------------------------------
        IF(.NOT. O_other) &
        ALLOCATE(               VP (NUMNT,   NUMNT+NCORE-1)         , &
                                VP2(NUMNT,   NUMNT+NCORE-1)         , &
                                WP (NUMNT+1, NUMNT+NCORE-1)         , &
                                WP2(NUMNT+1, NUMNT+NCORE-1)         , &
                                WX1(NUMNT,   NUMNT+NCORE)           , &
                                LPI(NUMNT+NCORE-1)                  , &
                                U1 (NUMNT+NCORE,   NUMNT+NCORE)     , &
                                U2 (NUMNT+NCORE,   NUMNT+NCORE)     , &
                                TF (NUMNT+NCORE-1, NUMNT)           )
!-------------------------------------------------------------------------
  end subroutine Allocate_Mangintw
!------------------------------------------------------------------------------
  subroutine Deallocate_Mangintw (O_other)
!input
      LOGICAL       :: O_other
!
!-------------------------------------------------------------------------
        IF(.NOT. O_other) &
        DEALLOCATE(             VP,VP2,WP,WP2,WX1,LPI,U1,U2,TF)
!-------------------------------------------------------------------------
  end subroutine Deallocate_Mangintw
!------------------------------------------------------------------------------
  subroutine WEIGH11 (ND,R,NP,P,NDIM,NPDIM,Z)
!
!-----CALLED BY 'CONT_SRT' - M 4.97/2.05/6.06/12.14 - ADI
!
!     CALCULATION OF ANGULAR INTEGRATION WEIGHTS
!     THE WEIGHTS W1  ARE CALCULATED FOR INTERMESH POINTS
!     AND INTERPOLATION WEIGHTS Q01,Q02,Q1,Q2,Q3,Q4 (Q12,Q21) FOR H,N
!
!input
    INTEGER(i4b)            :: ND,NP,NDIM,NPDIM
    REAL(dp)    , DIMENSION (ND)         :: R
    REAL(dp)    , DIMENSION (NP)         :: P
    REAL(dp)    , DIMENSION (NDIM,NPDIM) :: Z
!local
    INTEGER(i4b)            :: JP,LMAX,LZ,L
    REAL(dp)                :: RL,RL2,RRR12,A,B,C,AA,BB,CC,WW1,WW3,W1LZ,&
                               W3LZ,RZ,RZQ6,RRRR20,RD1,RD2
!-------------------------------------------------------------------------
!
        VPLAST = Z(1,NP-1)/(R(1)+R(1))
!
     DO JP=1,NP-1
        LMAX=MIN(NP+1-JP,ND)
        LZ=LMAX-1
!
!***  0. AND 2. MOMENT. THE INTEGRATION IS PERFORMED IN THE Z VARIABLE
        DO L=1,LMAX
          RL   =R(L)
          RL2  =RL+RL
          RRR12=RL*RL*RL*12._dp
   !***  FIRST STEP IF JP=1
          IF(JP.EQ.1) THEN
            B=Z(L,1)
            A=Z(L,2)
            VP(L,JP)=(B-A)/RL2
            AA=A*A
            BB=B*B
            VP2(L,JP)=(B*(3._dp*BB-AA)-A*(BB+AA))/RRR12
          ELSE
            IF (L.NE.LMAX .OR. JP.LE.(NP-ND) ) THEN
     !***  INTERMEDIATE STEP
              A=Z(L,JP+1)
              B=Z(L,JP)
              C=Z(L,JP-1)
              VP(L,JP)=(C-A)/RL2
              AA=A*A
              BB=B*B
              CC=C*C
              VP2(L,JP)=(B*(CC-AA)+C*(CC+BB)-A*(BB+AA))/RRR12
            ELSE
     !***  LAST STEP, IMPLYING Z(L,JMAX)=0
              B=Z(L,JP-1)
              VP (L,JP)=B/RL2
              VP2(L,JP)=B*B*B/RRR12
            ENDIF
          ENDIF
        END DO
!
!**** 1.MOMENT.
!***  FIRST STEP
        IF (JP.EQ.1 ) THEN
          WW1=P(2)*P(2)
          WW3=WW1*WW1
        ELSE
!***  INTERMEDIATE STEPS
          A=P(JP-1)
          B=P(JP)
          C=P(JP+1)
          WW1=(A+B+C)*(C-A)
          AA=A*A
          BB=B*B
          CC=C*C
          WW3=(B-A)*(AA*(A+2._dp*B)+BB*(3._dp*A+4._dp*B))+ &
              (C-B)*(CC*(C+2._dp*B)+BB*(3._dp*C+4._dp*B))
!***  FOR THE LAST INTERVAL (L=LZ TO THE P AXIS), THE NEXT POINT IS AN
!***  INTERMESH-POINT IN P
          C=.5_dp*(B+C)
          W1LZ=(A+B+C)*(C-A)
          CC=C*C
          W3LZ=(B-A)*(AA*(A+2._dp*B)+BB*(3._dp*A+4._dp*B))+ &
               (C-B)*(CC*(C+2._dp*B)+BB*(3._dp*C+4._dp*B))
        ENDIF
!***  NO WEIGHT FOR THE Z=0 POINT IS CALCULATED, AS V=0 THERE FOR SYMMETRY .
!***  LOOP OVER DEPTH INDEX L
        DO L=1,LZ
          RZ=.5_dp*(R(L)+R(L+1))
          RZQ6=RZ*RZ*6._dp
          RRRR20=RZ*RZ*RZ*RZ*20._dp
          IF (L.NE.LZ .OR. JP.LE.(NP-ND) ) THEN
            WP (L+1,JP)=WW1/RZQ6
            WP2(L+1,JP)=WP(L+1,JP)-WW3/RRRR20
          ELSE
            WP (L+1,JP)=W1LZ/RZQ6
            WP2(L+1,JP)=WP(L+1,JP)-W3LZ/RRRR20
          ENDIF
        END DO
!***  SPECIAL WEIGHTS AT THE OUTER BOUNDARY FOR H
        RL=R(1)
        WP (1,JP)=WW1/RL/RL/6._dp
        WP2(1,JP)=WP(1,JP)-WW3/RL/RL/RL/RL/20._dp
!***  SPECIAL WEIGHTS AT THE INNER BOUNDARY
      IF (LMAX.LT.ND) CYCLE
        IF (JP.GT.(NP-ND))THEN
  !***  CORE TANGENT RAY, LAST STEP OF INTEGRATION
          WW1=(B-A)*(2._dp*B+A)
          WW3=(B-A)*(AA*(A+2._dp*B)+BB*(3._dp*A+4._dp*B))
        ENDIF
        WP (ND+1,JP)=WW1/6._dp
        WP2(ND+1,JP)=WP(ND+1,JP)-WW3/20._dp
     END DO
!
!***  INTERPOLATION WEIGHTS
     A=-LOG10(2._dp)
      DO L=1,ND-2
       IF (ABS(R(L)-1.).LT.1.E-11) R(L)=1._dp
        RD1   =(R(L)+R(L+1))/ R(L+1)
        RD2   =(R(L)+R(L+1))/(R(L+1)+R(L+2))
        Q02(L)=(R(L+1)-R(L))/(R(L+2)-R(L))
        Q01(L)=(1._dp-Q02(L))/R(L+1)**2
        Q02(L)=       Q02(L) /R(L+1)**2
        Q4 (L)=Q02(L)*(R(L+1)+R(L+2))**2/4._dp
        Q3 (L)=Q01(L)*(R(L  )+R(L+1))**2/4._dp
        Q2 (L)=(A+LOG10(RD1))/LOG10(RD2)
        Q1 (L)=1._dp-Q2(L)
        Q21(L)=1._dp/R (L+1)**2
        Q12(L)=((R(L  )+R(L+1))/2._dp)**(2._dp*Q1(L))*  &
               ((R(L+1)+R(L+2))/2._dp)**(2._dp*Q2(L))*Q21(L)
      END DO
!
    RETURN
!-------------------------------------------------------------------------
  end subroutine WEIGH11
!------------------------------------------------------------------------------
  subroutine WEIGHPN (ND,R,NP,P,NDIM,NPDIM,Z,RRATIO)
!
!-----CALLED BY 'CONT_SRT' - M 4.97/2.05/6.06 - ADI
!
!     CALCULATION OF ANGULAR INTEGRATION WEIGHTS for planetary nebula
!     THE WEIGHTS W1 AND W3 ARE CALCULATED FOR INTERMESH POINTS
!     AND INTERPOLATION WEIGHTS Q01,Q02,Q1,Q2,Q3,Q4 (Q12,Q21) FOR H,N
!
!input
    INTEGER(i4b)            :: ND,NP,NDIM,NPDIM
    REAL(dp)    , DIMENSION (ND)         :: R
    REAL(dp)    , DIMENSION (NP)         :: P
    REAL(dp)    , DIMENSION (NDIM,NPDIM) :: Z
    REAL(dp)                :: RRATIO
!local
    INTEGER(i4b)            :: JP,LMAX,LZ,L
    REAL(dp)                :: RL,RL2,RRR12,A,B,C,AA,BB,CC,WW1,WW3,W1LZ,&
                               W3LZ,RZ,RZQ6,RRRR20,RD1,RD2
!-------------------------------------------------------------------------
      DO 100 JP=1,NP-1

      LMAX=MIN(NP+1-JP,ND)
      LZ=LMAX-1

! 0. AND 2. MOMENT. THE INTEGRATION IS PERFORMED IN THE Z VARIABLE
      DO 10 L=1,LMAX
        RL=R(L)
        RL2=RL+RL
        RRR12=RL*RL*RL*12._dp
! CENTER RAY = FIRST STEP
        IF (JP.EQ.1) THEN
          VP (L,JP)=.5_dp*(RRATIO/RL)**2
          VP2(L,JP)=.5_dp*(RRATIO/RL)**2
        ELSEIF (L.NE.LMAX .OR. JP.LE.(NP-ND) ) THEN
! INTERMEDIATE STEP
          A=Z(L,JP+1)
          B=Z(L,JP)
          C=Z(L,JP-1)
          VP(L,JP)=(C-A)/RL2
          AA=A*A
          BB=B*B
          CC=C*C
          VP2(L,JP)=(B*(CC-AA)+C*(CC+BB)-A*(BB+AA))/RRR12
        ELSE
! LAST STEP, IMPLYING Z(L,JMAX)=0
          B=Z(L,JP-1)
          VP(L,JP)=B/RL2
          VP2(L,JP)=B*B*B/RRR12
        ENDIF
10    CONTINUE

! 1. AND 3. MOMENT
! CENTER RAY = FIRST STEP
      IF (JP.EQ.1) THEN
        WP (1,JP)=.5_dp*(RRATIO/R(1))**2
        WP2(1,JP)=.5_dp*(RRATIO/R(1))**2
        DO L=1,LZ
          WP (L+1,JP)=.5_dp*(RRATIO*2/(R(L)+R(L+1)))**2
          WP2(L+1,JP)=.5_dp*(RRATIO*2/(R(L)+R(L+1)))**2
        ENDDO
        WP (ND+1,JP)=.5_dp*RRATIO**2
        WP2(ND+1,JP)=.5_dp*RRATIO**2
      ELSE
! INTERMEDIATE STEPS
        A=P(JP-1)
        B=P(JP)
        C=P(JP+1)
        WW1=(A+B+C)*(C-A)
        AA=A*A
        BB=B*B
        CC=C*C
        WW3=(B-A)*(AA*(A+2._dp*B)+BB*(3._dp*A+4._dp*B))+ &
            (C-B)*(CC*(C+2._dp*B)+BB*(3._dp*C+4._dp*B))
! FOR THE LAST INTERVAL (L=LZ TO THE P AXIS), THE NEXT POINT IS AN
! INTERMESH-POINT IN P
        C=.5_dp*(B+C)
        W1LZ=(A+B+C)*(C-A)
        CC=C*C
        W3LZ=(B-A)*(AA*(A+2._dp*B)+BB*(3._dp*A+4._dp*B))+ &
             (C-B)*(CC*(C+2._dp*B)+BB*(3._dp*C+4._dp*B))
! NO WEIGHT FOR THE Z=0 POINT IS CALCULATED, AS V=0 THERE FOR SYMMETRY.

! LOOP OVER DEPTH INDEX L
        DO 20 L=1,LZ
          RZ    =.5_dp*(R(L)+R(L+1))
          RZQ6  =RZ*RZ*6._dp
          RRRR20=RZ*RZ*RZ*RZ*20._dp
          IF (L.NE.LZ .OR. JP.LE.(NP-ND) ) THEN
            WP (L+1,JP)=WW1/RZQ6
            WP2(L+1,JP)=WP(L+1,JP)-WW3/RRRR20
          ELSE
            WP (L+1,JP)=W1LZ/RZQ6
            WP2(L+1,JP)=WP(L+1,JP)-W3LZ/RRRR20
          ENDIF
20      CONTINUE

! SPECIAL WEIGHTS AT THE OUTER BOUNDARY FOR H
        RL=R(1)
        WP (1,JP)=WW1/RL/RL/6._dp
        WP2(1,JP)=WP(1,JP)-WW3/RL/RL/RL/RL/20._dp
! SPECIAL WEIGHTS AT THE INNER BOUNDARY
        IF (LMAX.LT.ND) GOTO 100
        IF (JP.GT.(NP-ND))THEN
! CORE TANGENT RAY, LAST STEP OF INTEGRATION
          WW1=(B-A)*(2._dp*B+A)
          WW3=(B-A)*(AA*(A+2._dp*B)+BB*(3._dp*A+4._dp*B))
        ENDIF
        WP (ND+1,JP)=WW1/6._dp
        WP2(ND+1,JP)=WP(ND+1,JP)-WW3/20._dp

      ENDIF

100   CONTINUE

! INTERPOLATION WEIGHTS
      DO L=1,ND-2
       IF (ABS(R(L)-1._dp).LT.1.E-11_dp) R(L)=1._dp
       RD1   =(R(L)+R(L+1))/ R(L+1)
       RD2   =(R(L)+R(L+1))/(R(L+1)+R(L+2))
       Q02(L)=(R(L+1)-R(L))/(R(L+2)-R(L))
       Q01(L)=(1._dp-Q02(L))/R(L+1)**2
       Q02(L)=       Q02(L) /R(L+1)**2
       Q4 (L)=Q02(L)*(R(L+1)+R(L+2))**2/4._dp
       Q3 (L)=Q01(L)*(R(L  )+R(L+1))**2/4._dp
       Q2 (L)=(-LOG10(2._dp)+LOG10(RD1))/LOG10(RD2)
       Q1 (L)=1._dp-Q2(L)
       Q21(L)=1._dp/R (L+1)**2
       Q12(L)=((R(L  )+R(L+1))/2._dp)**(2._dp*Q1(L))* &
              ((R(L+1)+R(L+2))/2._dp)**(2._dp*Q2(L))*Q21(L)
      END DO
!
      RETURN
!-------------------------------------------------------------------------
   end subroutine WEIGHPN
!------------------------------------------------------------------------------
   subroutine WEIGHT1(R,P,Z,ND,NCORE)

! WEIGHTS FOR INTEGRAL mu*u dmu, EVALUATED ON THE REGULAR GRID POINTS
! CALLED BY 'CONT_SRT' - M 2.05 - ADI
!
!input
    INTEGER(i4b)            :: ND,NCORE
    REAL(dp)    , DIMENSION(ND,ND+NCORE):: Z
    REAL(dp)    , DIMENSION(ND+NCORE)   :: P
    REAL(dp)    , DIMENSION(ND)         :: R
!local
    INTEGER(i4b), DIMENSION(ND)         :: JP
    INTEGER(i4b)            :: IR,IP
    REAL(dp)                :: A,B
!-------------------------------------------------------------------------
  DO IR = 1,ND
    JP(IR) = ND + NCORE - IR + 1
  ENDDO

  DO IR = 1,ND
    DO IP = 1,JP(IR)-1
      A = Z(IR,IP+1)
      B = Z(IR,IP)
      WX1(IR,IP) = B**2/3. - B*A/6. - A**2/6.
    ENDDO
    DO IP = 2,JP(IR)
      A = Z(IR,IP)
      B = Z(IR,IP-1)
      WX1(IR,IP) =WX1(IR,IP) + B**2/6. + B*A/6. - A**2/3.
    ENDDO
    WX1(IR,:)=WX1(IR,:)/R(IR)**2
  ENDDO

! DO IR=1,ND
!   PRINT*,IR,JP(IR),SUM(W1(IR,1:JP(IR)))
! ENDDO

   RETURN
!-------------------------------------------------------------------------
  end subroutine WEIGHT1
!------------------------------------------------------------------------------
!
end Module M_Angular_Integration_weights
!==============================================================================
!           ***** END MODULE M_Angular_Integration_weights *****
!==============================================================================
!
!==============================================================================
!               ***** MODULE M_Frequency_Integration *****
!==============================================================================
module M_Frequency_Integration
!
   USE M_type_definitions
   USE M_physical_constants             , ONLY: CLIGHT

CONTAINS
!------------------------------------------------------------------------------
  function FREINT(QUAN,W,ILOW,IMAX,N,X,IOPT,CHAR)
!
!!DIR$ INLINE ALWAYS FREINT
!
!     INTEGRATES QUAN OVER FREQUENCY FROM ILOW TO IMAX
!
!     CHAR = 'NEW'..... NEW INTEGRATION WEIGHTS ARE CALCULATED IN FREWE
!     CHAR = 'NOI'..... ONLY NEW INTEGRATION WEIGHTS ARE CALCULATED IN FREWE
!     CHAR = 'OLD'..... OLD INTEGRATION WEIGHTS ARE USED
!
!     X..... INPUT FOR FREWE AS DESCRIBED THERE
!     IOPT.. INPUT FOR FREWE DESCRIBING X
!   M 2.05 - ADI
!
   IMPLICIT NONE
!
!input
    INTEGER(i4b)            :: N
    REAL(dp), DIMENSION (N) :: QUAN,X, W
!
    CHARACTER(LEN=3)        :: CHAR
    INTEGER(i4b)            :: IOPT
    INTEGER(i4b)            :: ILOW,IMAX
!output
    REAL(dp)                :: FREINT
!local
    INTEGER(i4b)            :: I
    REAL(dp)                :: SUM
!-------------------------------------------------------------------------
      IF(CHAR.EQ.'NEW' .or. CHAR.EQ.'NOI') THEN
       CALL FREWE
       IF(CHAR.EQ.'NOI')THEN
        FREINT=0._dp;RETURN
       ENDIF
      ELSE IF(CHAR.NE.'OLD') THEN
       STOP 'WRONG OPTION FOR CHAR IN FREINT'
      ENDIF
!
      SUM=0._dp
       DO I=ILOW,IMAX
        SUM=SUM+QUAN(I)*W(I)
       END DO
!
      FREINT=SUM
!
      RETURN
!-------------------------------------------------------------------------
 contains
!------------------------------------------------------------------------------
    subroutine FREWE
!
!-----CALCULATES INTEGRATION WEIGHTS: INPUT X, DIM(X) = N
!-----CALLED BY 'FREINT' - M 1.96/8.03 ADI - FILE3
!
!     IOPT=1 FREQUENCIES ( 1/S)
!     IOPT=2 WAVELENGTHS (   A)
!     IOPT=3 WAVELENGTHS (  CM)
!     IOPT=4 KAYSER      (1/CM)
!
!     X CAN BE ORDERED FROM HIGH TO LOW OR VICE VERSA
!     OUTPUT W, FREQUENCY INTEGRATION WEIGHTS SUCH THAT
!     INTEGRAL(X-NUE DNUE) = SUM OVER I(X-I*W-I)
!
!input
!    INTEGER  IOPT
!    REAL(dp), DIMENSION (N) :: X
!output
!    REAL(dp), DIMENSION (N) :: W
!local
    LOGICAL  OPTINV,OPTOUT
!
    INTEGER(i4b), DIMENSION (N) :: INDEXF
    REAL(dp)    , DIMENSION (N) :: NUE,DEL,DEL1,DEL2
    INTEGER(i4b)            :: I,K,I2,KINT,NINT,KK
    REAL(dp)                :: SUM,DELNUE
!-------------------------------------------------------------------------
      OPTOUT=.FALSE.
      OPTINV=.FALSE.
!
      IF(IOPT.EQ.1.OR.IOPT.EQ.4) THEN
       IF(X(1).LT.X(N)) THEN
       OPTINV=.TRUE.
        DO I=1,N
         NUE(I)=X(N+1-I)
        END DO
       ELSE
        DO I=1,N
         NUE(I)=X(I)
        END DO
       END IF
       IF(IOPT.EQ.4) THEN
        DO I=1,N
         NUE(I)=CLIGHT*NUE(I)
        END DO
       END IF
      ELSE IF(IOPT.EQ.2.OR.IOPT.EQ.3) THEN
       IF(X(1).GT.X(N)) THEN
        OPTINV=.TRUE.
         DO I=1,N
          NUE(I)=X(N+1-I)
         END DO
       ELSE
        DO I=1,N
         NUE(I)=X(I)
        END DO
       END IF
        IF(IOPT.EQ.2) THEN
         DO I=1,N
          NUE(I)=CLIGHT*1.E8_dp/NUE(I)
         END DO
        ELSE
         DO I=1,N
          NUE(I)=CLIGHT/NUE(I)
         END DO
        END IF
      ELSE
       STOP 'WRONG INPUT OPTION'
      END IF
!
       DO I=1,N-1
!       DEL1(I)=NUE(I)-NUE(I+1)
        DEL1(I)=MAX(NUE(I)-NUE(I+1),1.E-250_dp) ! hack
       END DO
      K=1
      I=1
      INDEXF(1)=1
!
 150  IF(I+1.EQ.N) THEN
       DEL(K)=DEL1(I)
       K=K+1
       INDEXF(K)=N
       GOTO 300
      END IF
!
       I2=I+1
       IF(ABS(1.-DEL1(I)/DEL1(I2)).GT.1.E-10) GOTO 250
      DO I2=I+2,N-1
       DEL2(I2)=ABS(1._dp-DEL1(I)/DEL1(I2))
      END DO
      DO I2=I+2,N-1
       IF(DEL2(I2).GT.1.E-10) GOTO 250
      END DO
!
      DEL(K)=DEL1(I)
      K=K+1
      INDEXF(K)=N
      GOTO 300
!
 250  DEL(K)=DEL1(I)
      K=K+1
      INDEXF(K)=I2
      I=I2
      GOTO 150
!
 300  KINT=K
       IF(OPTOUT) THEN
        PRINT*
        PRINT*,' NUMBER OF EQUIDISTANT SUBINTERVALS = ',KINT-1
        PRINT*
       END IF
!
       W(1)=0._dp
      DO KK=2,KINT
       NINT=INDEXF(KK)-INDEXF(KK-1)+1
!-----------------------------------------------------------------------
      IF(NINT.EQ.1) THEN
       STOP ' WRONG NUMBER OF SUBINTERVALS'
      ELSE IF(NINT.EQ.2) THEN
       CALL FREWE2(INDEXF(KK-1),INDEXF(KK),DEL(KK-1))
      ELSE IF(NINT.EQ.3) THEN
       CALL FREWE3(INDEXF(KK-1),INDEXF(KK),DEL(KK-1))
      ELSE IF(NINT.EQ.4) THEN
       CALL FREWE4(INDEXF(KK-1),INDEXF(KK),DEL(KK-1))
      ELSE IF(MOD(NINT,4).EQ.1) THEN
       CALL FREWE5(INDEXF(KK-1),INDEXF(KK),DEL(KK-1))
      ELSE IF(MOD(NINT,4).EQ.2) THEN
       CALL FREWE6(INDEXF(KK-1),INDEXF(KK),DEL(KK-1))
      ELSE IF(MOD(NINT,4).EQ.3) THEN
       CALL FREWE7(INDEXF(KK-1),INDEXF(KK),DEL(KK-1))
      ELSE IF(MOD(NINT,4).EQ.0) THEN
       CALL FREWE8(INDEXF(KK-1),INDEXF(KK),DEL(KK-1))
      END IF
!-----------------------------------------------------------------------
      END DO
!
      SUM=0._dp
       DO K=1,N
        SUM=SUM+W(K)
       END DO
      DELNUE=NUE(1)-NUE(N)
      IF(ABS(1.-DELNUE/SUM).GT.1.E-10) STOP 'ERROR IN FREWE!'
!
      IF(.NOT.OPTINV) RETURN
!
       DO K=1,N
        NUE(K)=W(N+1-K)
       END DO
      DO K=1,N
       W(K)=NUE(K)
      END DO
!
      RETURN
!-------------------------------------------------------------------------
    end subroutine
!-----------------------------------------------------------------------
    subroutine FREWE2(I1,I2,DEL)
!-----------------------------------------------------------------------
!!DIR$ INLINE ALWAYS FREWE2
!
!     TRAPEZOIDAL RULE
!
    REAL(dp)      :: DEL
    INTEGER(i4b)  :: I1,I2
!
      W(I1)=W(I1)+.5_dp*DEL
      W(I2)=.5_dp*DEL
!
      RETURN
    end subroutine
!-----------------------------------------------------------------------
    subroutine FREWE3(I1,I2,DEL)
!-----------------------------------------------------------------------
!!DIR$ INLINE ALWAYS FREWE3
!
!     SIMPSON'S RULE
!
    REAL(dp)      :: DEL
    INTEGER(i4b)  :: I1,I2
!
      W(I1)  =W(I1)+DEL/3._dp
      W(I1+1)=4._dp*DEL/3._dp
      W(I2)  =DEL/3._dp
!
      RETURN
    end subroutine
!-----------------------------------------------------------------------
    subroutine FREWE4(I1,I2,DEL)
!-----------------------------------------------------------------------
!!DIR$ INLINE ALWAYS FREWE4
!
!     3/8 RULE
!
    REAL(dp)      :: DEL
    INTEGER(i4b)  :: I1,I2
!
      W(I1)  =W(I1)+3._dp*DEL/8._dp
      W(I1+1)=9._dp*DEL/8._dp
      W(I1+2)=9._dp*DEL/8._dp
      W(I2)  =3._dp*DEL/8._dp
!
      RETURN
    end subroutine
!-----------------------------------------------------------------------
    subroutine FREWE5(I1,I2,DEL)
!-----------------------------------------------------------------------
!!DIR$ INLINE ALWAYS FREWE5
!
!     MODIFIED SIMPSON'S RULE
!
    REAL(dp)      :: DEL
    INTEGER(i4b)  :: I1,I2,NINT,K,NK
!
      NINT=I2-I1+1
      DO 10 K=I1,I2
      NK=K-I1+1
      IF(NK.EQ.1) THEN
      W(K)=W(K)+(1._dp-1._dp/15._dp)*DEL/3._dp
      ELSE IF(NK.EQ.NINT) THEN
      W(K)=(1._dp-1._dp/15._dp)*DEL/3._dp
      ELSE IF(MOD(NK,2).EQ.0) THEN
      W(K)=4._dp*(1._dp+1._dp/15._dp)*DEL/3._dp
      ELSE IF(MOD(NK,4).EQ.1) THEN
      W(K)=2._dp*(1._dp-1._dp/15._dp)*DEL/3._dp
      ELSE
      W(K)=1.6_dp*DEL/3._dp
      ENDIF
10    CONTINUE
!
      RETURN
    end subroutine
!-----------------------------------------------------------------------
    subroutine FREWE6(I1,I2,DEL)
!-----------------------------------------------------------------------
!!DIR$ INLINE ALWAYS FREWE6
!
!     MODIFIED SIMPSON'S RULE + TRAPEZOIDAL RULE FOR LAST INTERVAL
!
    REAL(dp)      :: DEL
    INTEGER(i4b)  :: I1,I2,NINT,K,NK
!
      NINT=I2-I1+1
      DO 10 K=I1,I2
      NK=K-I1+1
      IF(NK.EQ.1) THEN
      W(K)=W(K)+(1._dp-1._dp/15._dp)*DEL/3._dp
      ELSE IF(NK.EQ.NINT-1) THEN
      W(K)=(1._dp-1._dp/15._dp)*DEL/3._dp+DEL/2._dp
      ELSE IF(NK.EQ.NINT) THEN
      W(K)=DEL/2._dp
      ELSE IF(MOD(NK,2).EQ.0) THEN
      W(K)=4._dp*(1._dp+1._dp/15._dp)*DEL/3._dp
      ELSE IF(MOD(NK,4).EQ.1) THEN
      W(K)=2._dp*(1._dp-1._dp/15._dp)*DEL/3._dp
      ELSE
      W(K)=1.6_dp*DEL/3._dp
      ENDIF
10    CONTINUE
!
      RETURN
    end subroutine
!-----------------------------------------------------------------------
    subroutine FREWE7(I1,I2,DEL)
!-----------------------------------------------------------------------
!!DIR$ INLINE ALWAYS FREWE7
!
!     SIMPSON'S RULE
!
    REAL(dp)      :: DEL
    INTEGER(i4b)  :: I1,I2,NINT,K,NK
!
      NINT=I2-I1+1
      DO 10 K=I1,I2
      NK=K-I1+1
      IF(NK.EQ.1) THEN
      W(K)=W(K)+DEL/3._dp
      ELSE IF(NK.EQ.NINT) THEN
      W(K)=DEL/3._dp
      ELSE IF(MOD(NK,2).EQ.0) THEN
      W(K)=4._dp*DEL/3._dp
      ELSE
      W(K)=2._dp*DEL/3._dp
      ENDIF
10    CONTINUE
!
      RETURN
    end subroutine
!-----------------------------------------------------------------------
    subroutine FREWE8(I1,I2,DEL)
!-----------------------------------------------------------------------
!!DIR$ INLINE ALWAYS FREWE8
!
!     SIMPSON'S RULE + TRAPEZOIDAL RULE FOR LAST INTERVAL
!
    REAL(dp)      :: DEL
    INTEGER(i4b)  :: I1,I2,NINT,K,NK
!
      NINT=I2-I1+1
      DO 10 K=I1,I2
      NK=K-I1+1
      IF(NK.EQ.1) THEN
      W(K)=W(K)+DEL/3._dp
      ELSE IF(NK.EQ.NINT-1) THEN
      W(K)=DEL/3._dp+DEL/2._dp
      ELSE IF(NK.EQ.NINT) THEN
      W(K)=DEL/2._dp
      ELSE IF(MOD(NK,2).EQ.0) THEN
      W(K)=4._dp*DEL/3._dp
      ELSE
      W(K)=2._dp*DEL/3._dp
      ENDIF
10    CONTINUE
!
      RETURN
    end subroutine
!-----------------------------------------------------------------------
  end function FREINT
!------------------------------------------------------------------------------
!
end module M_Frequency_Integration
!==============================================================================
!             ***** END MODULE M_Frequency_Integration *****
!==============================================================================
!
!==============================================================================
!                ***** MODULE M_Inter_Extrapolation *****
!==============================================================================
Module M_Inter_Extrapolation
!--- routines used for interpolation or extrapolation
!     Monpcf
!     Spline
!     Aitken
!     Linear
!--- 2.05 - ADI

   USE M_type_definitions
!
   IMPLICIT NONE
!
CONTAINS
!------------------------------------------------------------------------------
  subroutine MONPCF(X,Y,A,B,C,N)
!
!-----2.94 M 2.05 ADI
!-----COMPUTATION OF A MONOTONIC PIECEWISE CUBIC FUNCTION
!-----METHOD DESCRIBED BY STEFFEN A&A 239,443
!-----WITH COEFFICIENTS A,B,C A DATA ARRAY (X,Y) IS APPROXIMATTED
!-----BY F(X)=Y(I)+DEL*(A(I)+DEL*(B(I)+C(I)*DEL)) IN THE INTERVAL
!-----X(I)-X-X(I+1), WHERE DEL=X-X(I).
!-----VANISHING CURVATURE AT THE BOUNDARIES IS ASSUMED:
!-----F"(X(1))=F"(X(N))=0
!
!input
      INTEGER(i4b)  :: N
      REAL(dp)    , DIMENSION (N)            :: X,Y
!output
      REAL(dp)    , DIMENSION (N)            :: A,B,C
!local
      REAL(dp)    , DIMENSION (N), AUTOMATIC :: XX,Z
      INTEGER(i4b)  :: I
!-------------------------------------------------------------------------
       DO I=1,N-1
        XX(I)= X(I+1)-X(I)
        Z (I)=(Y(I+1)-Y(I))/XX(I)
       END DO
!
       DO I=2,N-1
        A (I)=(Z(I-1)*XX(I)+Z(I)*XX(I-1)) &
                  /(XX(I-1)+XX(I))
        A (I)=(SIGN(1._dp,Z(I-1))+SIGN(1._dp,Z(I)))* &
               MIN(ABS(Z(I-1)),ABS(Z(I)),.5_dp*ABS(A(I)))
       END DO
!----------------------------------------------------------
!-      A (1)=Z(1)*(1.+XX(1)/(XX(1)+XX(2)))
!-   *       -Z(2)*XX(1)/(XX(1)+XX(2))
!-      A (1)=(SIGN(1.,A(1))+SIGN(1.,Z(1)))*
!-   *        MIN(ABS(Z(1)),.5*ABS(A(1)))
!-      A (N)=Z(N-1)*(1.+XX(N-1)/(XX(N-1)+XX(N-2)))
!-   *       -Z(N-2)*XX(N-1)/(XX(N-1)+XX(N-2))
!-      A (N)=(SIGN(1.,A(N))+SIGN(1.,Z(N-1)))*
!-   *        MIN(ABS(Z(N-1)),.5*ABS(A(N)))
!----------------------------------------------------------
        A (1)=.5_dp*(3._dp*Z(  1)-A(  2))
        A (N)=.5_dp*(3._dp*Z(N-1)-A(N-1))
!----------------------------------------------------------
!
      DO I=1,N-1
       B(I)=( 3._dp*Z(I)-2._dp*A(I)-A(I+1))/XX(I)
       C(I)=(-2._dp*Z(I)+      A(I)+A(I+1))/XX(I)**2
      END DO
!
      RETURN
!-------------------------------------------------------------------------
  end subroutine MONPCF
!------------------------------------------------------------------------------
  subroutine SPLINE(X,Y,A,B,C,N)
!
!-----M 2.94/2.05  ADI - VECTORIZED VERSION!
!
!-----COMPUTATION OF CUBIC SPLINE COEFFICIENTS A,B,C
!-----FOR APPROXIMATION OF A DATA ARRAY (X,Y) SUCH THAT
!-----F(X)=Y(I)+DEL*(A(I)+DEL*(B(I)+C(I)*DEL)) IN THE
!-----INTERVAL X(I)-X-X(I+1), WHERE DEL=X-X(I).
!-----VANISHING CURVATURE AT THE BOUNDARIES IS ASSUMED:
!-----F"(X(1))=F"(X(N))=0
!
!input
      INTEGER(i4b)  :: N
      REAL(dp)    , DIMENSION (N)            :: X,Y
!output
      REAL(dp)    , DIMENSION (N)            :: A,B,C
!local
      REAL(dp)    , DIMENSION (N), AUTOMATIC :: XX,Z
      REAL(dp)      :: X1
      INTEGER(i4b)  :: I,J
!-------------------------------------------------------------------------
       DO I=2,N
        XX(I)=X(I)-X(I-1)
       END DO
       DO I=2,N-1
        A(I)=XX(I+1)/(XX(I)+XX(I+1))
        B(I)=1._dp-A(I)
        C(I)=6._dp*((Y(I+1)-Y(I))/XX(I+1)-(Y(I)-Y(I-1))/XX(I)) &
                   /(XX(I)+XX(I+1))
       END DO
!--------------------------------------------------------------------
       DO I=1,N-1
        XX(I)=2._dp**(I-1)
       END DO
        DO I=3,N-1
         Z (I)=-B(I)*A(I-1)
        END DO
       DO I=3,N-1
        X1=XX(I-2)
         DO J=I,N-1
          XX(J)=XX(J)+Z(I)*X1*2._dp**(J-I)
         END DO
       END DO
        DO I=3,N-1
         B(I)=XX(I)/XX(I-1)
        END DO
!----------------------------------------------------
        DO I=3,N-1
         XX(I)=-(1._dp-A(I))/B(I-1)
         A(I)=C(I-1)
        END DO
       DO I=3,N-1
        Z(I)=1._dp
         DO J=3,I
          Z(J)=Z(J)*XX(I)
          C(I)=C(I)+Z(J)*A(J)
         END DO
       END DO
!----------------------------------------------------
        DO I=2,N
         XX(I)=X(I)-X(I-1)
        END DO
       DO I=2,N-1
        A(I)=XX(I+1)/(XX(I)+XX(I+1))
       END DO
!----------------------------------------------------
      C(N-1)=C(N-1)/B(N-1)
       DO I=N-2,2,-1
        C(I)= C(I)/B(I)
        B(I)=-A(I)/B(I)
        A(I)= C(I+1)
       END DO
       DO I=N-2,2,-1
        Z(I)=1._dp
         DO J=N-2,I,-1
          Z(J)=Z(J)*B(I)
          C(I)=C(I)+Z(J)*A(J)
         END DO
       END DO
      C(1)=.0_dp
      C(N)=.0_dp
!--------------------------------------------------------------------
      DO I=1,N-1
       A(I)=(Y(I+1)-Y(I))/XX(I+1)-XX(I+1)*(2._dp*C(I)+C(I+1))/6._dp
       B(I)=.5_dp*C(I)
       C(I)=(C(I+1)-C(I))/XX(I+1)/6._dp
      END DO
!
      RETURN
!-------------------------------------------------------------------------
  end subroutine SPLINE
!------------------------------------------------------------------------------
  function AITKEN(X,Y,M,UP)
!
!-----M 2.05  ADI
!
!-----AITKEN CALCULATES A NEW GUESS FOR THE RADIUS OF THE SINGULAR POINT.
!-----THE CALCULATION RUNS WITH AITKENINTERPOLATION.
!-----UP IS THE CORRECT PHOTOSPHERICAL POINTS
!-----X(M)=ARE THE GUESSED PHOTOSPHERICAL POINTS
!-----Y(M)=ARE THE GUESSED SINGULAR POINTS
!-----M HAS TO BE GREATER(EQUAL) 2.
!
!input
      INTEGER(i4b)  :: M
      REAL(dp)      :: UP
      REAL(dp)    , DIMENSION (M)            :: X,Y
!output
      REAL(dp)      :: AITKEN
!local
      INTEGER(i4b)  :: I,K,J
!-------------------------------------------------------------------------
      DO I=1,M-1
         K=I+1
       DO J=K,M
        Y(J)=((X(J)-UP)*Y(I)-(X(I)-UP)*Y(J))/(X(J)-X(I))
       END DO
      END DO
!
      AITKEN=Y(M)
!
      RETURN
!-------------------------------------------------------------------------
  end function AITKEN
!------------------------------------------------------------------------------
  subroutine LINEAR (X,Y,A,N,M)
!
!-----COMPUTATION OF LINEAR INTERPOLATION COEFFICIENTS A.
!
!input
      INTEGER(i4b)  :: N,M
      REAL   (dp ), DIMENSION(N) :: X,Y
!output
      REAL   (dp ), DIMENSION(N) :: A
!local
      INTEGER(i4b)  :: I
!-------------------------------------------------------------------------
      DO I=M,N-1
       A(I)=(Y(I+1)-Y(I))/(X(I+1)-X(I))
      END DO
!
      RETURN
!-------------------------------------------------------------------------
  end subroutine LINEAR
!------------------------------------------------------------------------------
!
end module M_Inter_Extrapolation
!==============================================================================
!              ***** END MODULE M_Inter_Extrapolation *****
!==============================================================================
!
!==============================================================================
!               ***** MODULE M_Solve_lin_eq_systems *****
!==============================================================================
module M_Solve_lin_eq_systems
!--- routines used for solving linear equations
!     FNLUD
!     FNLUDS
!     LUDCMP
!     LUBKSB
!     Simple_mat_inv
!
!--- 8.06 - ADI
!
   USE M_type_definitions
!-----
   USE M_global_settings                , ONLY:
!-----
   IMPLICIT NONE
!
CONTAINS
!------------------------------------------------------------------------------
  subroutine FNLUD (N,A,LDA,IPIV,TEMP,INFO,OPTPIV)
!
!!DIR$ INLINE ALWAYS FNLUD
!
!  =====================================================================
!     BASIC ALGORITHM FROM NAG(F07ADF) - 11.95 ADI
!                                      M  3.96/3.05/8.06 ADI
!
!  Purpose
!  =======
!
!  An LU factorization of a general n-by-n matrix A
!  using partial pivoting with row interchanges is computed.
!
!  The factorization has the form
!     A = P * L * U
!  where P is a transposition matrix, L is lower triangular with unit
!               --------------------
!  diagonal elements , and U is upper triangular .
!
!  Crout algorithm.
!
!  Arguments
!  =========
!
!  N       (input) INTEGER
!          The number of rows and columns of the matrix A.  N >= 0.
!
!  A       (input/output) REAL array, dimension (LDA,N)
!          On entry, the m by n matrix to be factored.
!          On exit, the factors L and U from the factorization
!          A = P*L*U; the unit diagonal elements of L are not stored.
!
!  LDA     (input) INTEGER
!          The leading dimension of the array A.  LDA >= max(1,N).
!
!  IPIV    (input) INTEGER array, dimension (min(M,N)) IPIV(i)=i
!          (output)                                    ---------
!          The pivot indices; for 1 <= i <= min(M,N), row i of the
!          matrix was interchanged with row IPIV(i).
!
!  INFO    (output) INTEGER
!          = 0: successful exit
!          < 0: if INFO = -k, the k-th argument had an illegal value
!          > 0: if INFO = k, U(k,k) is exactly zero. The factorization
!               has been completed, but the factor U is exactly
!               singular, and division by zero will occur if it is used
!               to solve a system of equations.
!
!  =====================================================================
!
!in-output
      LOGICAL       :: OPTPIV
      INTEGER(i4b)  :: N,LDA,INFO
      REAL    (dp), DIMENSION (:,:) :: A
      INTEGER(i4b), DIMENSION (:)   :: IPIV
      REAL    (dp), DIMENSION (LDA) :: TEMP
!local
      INTEGER(i4b)  :: IX,IJP,J,JB,JP,JX,L,M
      INTEGER(i4b), DIMENSION (LDA) :: IPIV1
!     .. Parameters ..STANDARD...NB=24
      INTEGER(i4b), PARAMETER :: NB=64
      REAL    (dp), PARAMETER :: TINY=1.0E-20_dp
!------------------------------------------------------------------------------
!
       INFO = 0
!
!---------------------------------------------------------------------
!        Use blocked code.
!
         DO 40 J = 1, N, NB
           JB = MIN(N-J+1,NB)
           L  = J-1
           M  = JB+L
!
!           Update diagonal and subdiagonal blocks.
!           Form  A(J,J):=A(J,J)-A(1,J)*A(J,1).
!
!            IF (J.NE.1) THEN
!               DO IY = J, M
!                DO IZ = 1, L
!                 DO IX = J, N
!                  A(IX,IY) = A(IX,IY) - A(IZ,IY)*A(IX,IZ)
!                 END DO
!                END DO
!               END DO
!            END IF
!--------------------------------------------------------------
!          The _GEMM routine performs the following operations:
!          C  = alpha(op)A(op)B +beta*C - BLAS3
!          (transa,transb,m,n,k,alpha,a,lda,b,ldb,beta,c,ldc)
!
           CALL SGEMM('N','N',N-L,JB,L,-1., &
                       A(J,1),LDA,A(1,J),LDA,1.,A(J,J),LDA)
!--------------------------------------------------------------
!
!---------------------------------------------------------------------
         IF (.NOT.OPTPIV) THEN
!---------------------------------------------------------------------
!           Factorize diagonal and subdiagonal blocks and test for
!           singularity.
!
          DO JX = J, M
!
!---------------------------------------------------------------
!         Update diagonal and subdiagonal elements in column J.
!         Form  Ay := Ay -Ax*A .
!
!            DO  IX = J, JX-1
!              DO  IY = JX, N
!               A(IY,JX) = A(IY,JX) - A(IX,JX)*A(IY,IX)
!              END DO
!            END DO
!--------------------------------------------------------------
!          The _GEMV routine performs the following operations:
!          y  =  alpha*Ax + beta*y - BLAS2
!          (trans, m, n, alpha, a, lda, x, incx, beta, y, incy)
!
           CALL SGEMV('N',N-JX+1,JX-J,-1.,A(JX,J),LDA,A(J,JX),1 &
                      ,1.,A(JX,JX),1)
!--------------------------------------------------------------
!
!---------------------------------------------------------------
!         Test for singularity.
!
          IF (A(JX,JX).EQ.0.) THEN
!
!          If A( JP, J ) is zero, set INFO to indicate that a
!          zero pivot has been found.
!
            IF (INFO.EQ.0) INFO = JX
            A(JX,JX) = TINY
          END IF
!
!         Compute elements J+1:N of J-th column.
!
!          DO  IY = JX+1, N
!            A(IY,JX) = A(IY,JX)/A(JX,JX)
!          END DO
!-----------------------------------------------------
!            {S,D}SCAL (n, alpha, x, incx) - BLAS1
!
           IF (JX < N) & !11.99
             CALL SSCAL(N-JX,1./A(JX,JX),A(JX+1,JX),1)
!-----------------------------------------------------
!
!---------------------------------------------------------------
!           Compute block row of U.
!           Form  Ay := Ay -A'*Ax .
!
          IF (JX.NE.J .and. JX.NE.N) THEN  !3.05 Adi
!           DO  IX = 1, M-JX
!            TEMP(IX) = 0._dp
!           END DO
!            DO  IX = 1, M-JX
!             DO  IY = 1, JX-J
!              TEMP(IX)=TEMP(IX)+A(IY+L,JX+IX)*A((IY-1)*LDA+JX,J)
!             END DO
!            END DO
!           DO IY = 1, M-JX
!            A((IY-1)*LDA+JX,JX+1)=A((IY-1)*LDA+JX,JX+1)-TEMP(IY)
!           END DO
!---------------------------------------------------------------
!           The _GEMV routine performs the following operations:
!           y  =  alpha*Ax + beta*y - BLAS2
!           (trans, m, n, alpha, a, lda, x, incx, beta, y, incy)
!
            CALL SGEMV('T',JX-J,M-JX,-1.,A(L+1,JX+1),LDA,A(JX,J) &
                       ,LDA,1.,A(JX,JX+1),LDA)
!---------------------------------------------------------------
          END IF
!
          END DO
!---------------------------------------------------------------------
         ELSE
!---------------------------------------------------------------------
!         Factorize diagonal and subdiagonal blocks and test for
!         singularity using partial pivoting with row interchanges.
!
          DO JX = J, M
!---------------------------------------------------------------
!         Update diagonal and subdiagonal elements in column J.
!         Form  Ay := Ay -Ax*A .
!
!            DO  IX = J, JX-1
!              DO  IY = JX, N
!               A(IY,JX) = A(IY,JX) - A(IX,JX)*A(IY,IX)
!              END DO
!            END DO
!--------------------------------------------------------------
           CALL SGEMV('N',N-JX+1,JX-J,-1.,A(JX,J),LDA,A(J,JX),1 &
                      ,1.,A(JX,JX),1)
!--------------------------------------------------------------
!
!---------------------------------------------------------------
!         Find pivot.
!
!          JP = JX
!          XMAX = ABS(A(JX,JX))
!            DO  IY = 1, N-JX
!               IF( XMAX.LT.ABS( A(JX+IY,JX) ) )THEN
!                  XMAX = ABS( A(JX+IY,JX) )
!                  JP = JX+IY
!               END IF
!            END DO
!---------------------------------------------
!          I{S,D}AMAX (n,x,incx) - BLAS1
!
           JP = JX-1+ISAMAX(N-JX+1,A(JX,JX),1)
!---------------------------------------------
           IPIV1(JX) = JP
           IJP       = IPIV(JX)
           IPIV(JX)  = IPIV(JP)
           IPIV(JP) = IJP
!
!           Apply interchange to columns 1:N.
!
!            DO  IY = 1, JB
!               TEMP (IY)          = A(JX+(IY-1)*LDA,J)
!               A(JX+(IY-1)*LDA,J) = A(JP+(IY-1)*LDA,J)
!               A(JP+(IY-1)*LDA,J) = TEMP(IY)
!            END DO
!----------------------------------------------------
!           {S,D}SWAP (n,x,incx,y,incy) - BLAS1
!
            IF (JP.NE.JX) &
            CALL SSWAP(JB,A(JX,J),LDA,A(JP,J),LDA)
!----------------------------------------------------
!
!---------------------------------------------------------------
!         Test for singularity.
!
          IF (A(JX,JX).EQ.0.) THEN
!
!          If A( JP, J ) is zero, set INFO to indicate that a
!          zero pivot has been found.
!
            IF (INFO.EQ.0) INFO = JX
            A(JX,JX) = TINY
          END IF
!
!         Compute elements J+1:N of J-th column.
!
!          DO  IY = JX+1, N
!            A(IY,JX) = A(IY,JX)/A(JX,JX)
!          END DO
!----------------------------------------------------
           IF (JX < N) & !11.99 Adi
            CALL SSCAL(N-JX,1./A(JX,JX),A(JX+1,JX),1)
!----------------------------------------------------
!
!---------------------------------------------------------------
!           Compute block row of U.
!           Form  Ay := Ay -A'*Ax .
!
          IF (JX.NE.J .and. JX.NE.N) THEN !3.05 Adi
!           DO  IX = 1, M-JX
!            TEMP(IX) = 0._dp
!           END DO
!            DO  IX = 1, M-JX
!             DO  IY = 1, JX-J
!              TEMP(IX)=TEMP(IX)+A(IY+L,JX+IX)*A((IY-1)*LDA+JX,J)
!             END DO
!            END DO
!           DO IY = 1, M-JX
!            A((IY-1)*LDA+JX,JX+1)=A((IY-1)*LDA+JX,JX+1)-TEMP(IY)
!           END DO
!----------------------------------------------------------------
            CALL SGEMV('T',JX-J,M-JX,-1.,A(L+1,JX+1),LDA,A(JX,J) &
                       ,LDA,1.,A(JX,JX+1),LDA)
!----------------------------------------------------------------
          END IF
!
          END DO
!---------------------------------------------------------------------
!           Apply the interchanges to the columns on either
!           side of the current block.
!
            DO IX = J, M
!             DO  IY = 1, L
!              TEMP (IY)                =A(IX+(IY-1)*LDA,1)
!              A(IX+(IY-1)*LDA,1)       =A(IPIV1(IX)+(IY-1)*LDA,1)
!              A(IPIV1(IX)+(IY-1)*LDA,1)=TEMP(IY)
!             END DO
!             DO  IY = M+1,N
!              TEMP (IY)                =A(IX+(IY-1)*LDA,1)
!              A(IX+(IY-1)*LDA,1)       =A(IPIV1(IX)+(IY-1)*LDA,1)
!              A(IPIV1(IX)+(IY-1)*LDA,1)=TEMP(IY)
!             END DO
!-----------------------------------------------------------------
!             {S,D}SWAP (n,x,incx,y,incy) - BLAS1
!
              IF (IPIV1(IX).NE.IX) THEN
               CALL SSWAP( L  ,A(IX,  1),LDA,A(IPIV1(IX),  1),LDA)
               IF (M.NE.N) &  !3.05 Adi
               CALL SSWAP( N-M,A(IX,M+1),LDA,A(IPIV1(IX),M+1),LDA)
              END IF
!-----------------------------------------------------------------
            END DO
!---------------------------------------------------------------------
         END IF
!---------------------------------------------------------------------
!         Compute block row of U.
!         Form  A(J  ,J+JB):=A(J  ,J+JB)-A(1,J+JB)*A(J  ,1).
!         and then
!         Form  A(J+1,J+JB):=A(J+1,J+JB)-A(J,J+JB)*A(J+1,J).
!
            IF (J.NE.1 .and. M.NE.N) THEN  !3.05 Adi
!               DO IY = M+1, N
!                DO IZ = 1, L
!                 DO IX = J, M
!                  A(IX,IY) = A(IX,IY) - A(IZ,IY)*A(IX,IZ)
!                 END DO
!                END DO
!               END DO
!---------------------------------------------------------------
!           The _GEMM routine performs the following operations:
!           C  = alpha(op)A(op)B +beta*C - BLAS3
!           (transa,transb,m,n,k,alpha,a,lda,b,ldb,beta,c,ldc)
!
             CALL SGEMM('N','N',JB,N-M,L,-1.,A(J,1), &
                        LDA,A(1,M+1),LDA,1.,A(J,M+1),LDA)
!---------------------------------------------------------------
            END IF
!               DO IY = M+1, N
!                DO IZ = J, M
!                 DO IX = IZ+1, M
!                  A(IX,IY) = A(IX,IY) - A(IZ,IY)*A(IX,IZ)
!                 END DO
!                END DO
!               END DO
!---------------------------------------------------------
!           The _TRSM routines solve a triangular system of
!           equations.A is a triangular matrix:
!           AX  =  alpha * B - BLAS3
!           alpha is a scalar, X and B are m by n matrices,
!           and A is a unit , lower triangular matrix.
!           (side,uplo,transa,diag,m,n,alpha,a,lda,b,ldb)
!
             IF (M.NE.N) &  !3.05 Adi
             CALL STRSM('L','L','N','U',JB, &
                        N-M,1.,A(J,J),LDA,A(J,M+1),LDA)
!---------------------------------------------------------
   40    CONTINUE
!
   RETURN
!------------------------------------------------------------------------------
  end subroutine FNLUD
!------------------------------------------------------------------------------
  subroutine FNLUDS (N,A,LDA,IPIV,B,TEMP,OPTPIV)
!
!!DIR$ INLINE ALWAYS FNLUDS
!
!  =====================================================================
!     BASIC ALGORITHM FROM NAG(F07AEF) -M 11.95/3.05/8.06 ADI
!
!  Purpose
!  =======
!
!   FNLUDS solves a system of linear equations
!     A * X = B
!  with a general n by n matrix A using the LU factorization computed
!  by FNLUD.
!
!  Arguments
!  =========
!
!  N       (input) INTEGER
!          The order of the matrix A.  N >= 0.
!
!  A       (input) REAL array, dimension (LDA,N)
!          The factors L and U from the factorization A = P*L*U
!          as computed by F07ADF.
!
!  LDA     (input) INTEGER
!          The leading dimension of the array A.  LDA >= max(1,N).
!
!  IPIV    (input) INTEGER array, dimension (N)
!          The pivot indices from F07ADF; for 1<=i<=N, row i of the
!          matrix was interchanged with row IPIV(i).
!
!  B       (input/output) REAL array, dimension (LDA)
!          On entry, the right hand side vectors B for the system of
!          linear equations.
!          On exit, the solution vectors, X.
!
!  =====================================================================
!
!inoutput
      LOGICAL       :: OPTPIV
      INTEGER(i4b)  :: N,LDA
      REAL    (dp), DIMENSION (:,:) :: A
      REAL    (dp), DIMENSION (:)   :: B,TEMP
      INTEGER(i4b), DIMENSION (:)   :: IPIV
!local
      INTEGER(i4b)  :: I
!------------------------------------------------------------------------------
!
!        Solve A * X = B.
!---------------------------------------------------------------------
!        Apply row interchanges to the right hand sides.
!
         IF (OPTPIV) THEN
          DO  I = 1, N
            TEMP( I ) = B( I )
          END DO
          DO  I = 1, N
            B( I )    = TEMP( IPIV(I) )
          END DO
         END IF
!---------------------------------------------------------------------
!
!           The _TRSV routines solve a triangular system of equations.
!           A is a triangular matrix:
!           AX  =  B - BLAS2
!           alpha is a scalar, X and B are n vectors,
!           and A is a Unit/Non-unit , Lower/Upper triangular matrix.
!           (uplo, trans, diag, n, a, lda, x, incx)
!---------------------------------------------------------------------
!        Solve L*X = B, overwriting B with X.
!
!            DO K = 1, N
!              DO I = K + 1, N
!                B(I) = B(I) - B(K)*A(I,K)
!              END DO
!            END DO
         CALL STRSV('L','N','U',N,A(1,1),LDA,B(1),1)
!---------------------------------------------------------------------
!        Solve U*X = B, overwriting B with X.
!
!            DO K = N, 1, -1
!               DO I = 1, K - 1
!                  B(I) = B(I) - B(K)*A(I,K)/A(K,K)
!               END DO
!            END DO
!
!            DO K = 1,N
!               B(K) = B(K)/A(K,K)
!            END DO
         CALL STRSV('U','N','N',N,A(1,1),LDA,B(1),1)
!---------------------------------------------------------------------
!
   RETURN
!------------------------------------------------------------------------------
  end subroutine FNLUDS
!------------------------------------------------------------------------------
  subroutine LUDCMP(A,N,NP,INDX,VV,D)
!
!  An LU factorization of a general n-by-n matrix A
!  using partial pivoting with row interchanges is computed.
!  =========================================
!  CALLED BY 'LUDFIT' - ADI 10.95/3.05 FILE3
!  =========================================
!
!input
     INTEGER(i4b)  :: N,NP
!output
     REAL    (dp)  :: D
     INTEGER(i4b), DIMENSION(NP)   :: INDX
     REAL    (dp), DIMENSION(NP)   :: VV
     REAL    (dp), DIMENSION(NP,NP):: A
!local
      INTEGER(i4b)  :: I,J,K,IMAX
      REAL    (dp)  :: SUM,DUM,AAMAX
      REAL    (dp), PARAMETER :: TINY=1.0E-20_dp
!------------------------------------------------------------------------------
!
      D=1._dp
       DO I=1,N
        VV(I)=1._dp
       END DO
!-----------------------------------------------------
!      DO I=1,N
!        AAMAX=0.
!        DO J=1,N
!          IF (ABS(A(I,J)).GT.AAMAX) AAMAX=ABS(A(I,J))
!        END DO
!        IF (AAMAX.EQ.0.) PAUSE 'SINGULAR MATRIX.'
!        VV(I)=1./AAMAX
!      END DO
!-----------------------------------------------------
      DO 11 J=1,N
          DO I=1,J-1
            SUM=A(I,J)
              DO K=1,I-1
                SUM=SUM-A(I,K)*A(K,J)
              END DO
              A(I,J)=SUM
          END DO
         AAMAX=0._dp
!
        DO I=J,N
          SUM=A(I,J)
            DO K=1,J-1
              SUM=SUM-A(I,K)*A(K,J)
            END DO
         A(I,J)=SUM
         DUM=VV(I)*ABS(SUM)
          IF (DUM.GE.AAMAX) THEN
            IMAX=I
            AAMAX=DUM
          ENDIF
        END DO
!
        IF (J.NE.IMAX)THEN
          DO K=1,N
            DUM=A(IMAX,K)
            A(IMAX,K)=A(J,K)
            A(J,K)=DUM
          END DO
          D=-D
          VV(IMAX)=VV(J)
        ENDIF
         INDX(J)=IMAX
         IF(A(J,J).EQ.0.)A(J,J)=TINY
        IF(J.NE.N)THEN
          DUM=1._dp/A(J,J)
          DO I=J+1,N
            A(I,J)=A(I,J)*DUM
          END DO
        ENDIF
   11 CONTINUE
!
   RETURN
!------------------------------------------------------------------------------
  end subroutine LUDCMP
!------------------------------------------------------------------------------
  subroutine LUBKSB(A,N,NP,INDX,B)
!
!  Solves a system of linear equations
!     A * X = B
!  with a general n by n matrix A using the LU factorization
!  ======================================================
!  CALLED BY 'LUDFIT' AND 'MPROVE' - ADI 10.95/3.05 FILE3
!  ======================================================
!
!input
     INTEGER(i4b)  :: N,NP
     REAL    (dp), DIMENSION(NP,NP):: A
     INTEGER(i4b), DIMENSION(NP)   :: INDX
!output
     REAL    (dp), DIMENSION(NP)   :: B
!local
      INTEGER(i4b)  :: I,J,II,LL
      REAL    (dp)  :: SUM
!------------------------------------------------------------------------------
!
      II=0
      DO I=1,N
        LL=INDX(I)
        SUM=B(LL)
        B(LL)=B(I)
        IF (II.NE.0)THEN
          DO J=II,I-1
            SUM=SUM-A(I,J)*B(J)
          END DO
        ELSE IF (SUM.NE.0.) THEN
          II=I
        ENDIF
        B(I)=SUM
      END DO
      DO I=N,1,-1
        SUM=B(I)
          DO J=I+1,N
            SUM=SUM-A(I,J)*B(J)
          END DO
        B(I)=SUM/A(I,I)
      END DO
!
   RETURN
!------------------------------------------------------------------------------
  end subroutine LUBKSB
!------------------------------------------------------------------------------
  subroutine Simple_mat_inv (N,NDIM,A)
! ==========================================
!     SCALED MATRIX INVERSION (column wise)
!     M 12.14 Adi
!     Originally called by CONT as INVC
! ==========================================
!
!input/output
      REAL   (dp ), DIMENSION (NDIM,NDIM) :: A
!input
      INTEGER(i4b) :: N,NDIM
!local
      REAL   (dp ), DIMENSION (NDIM)  , AUTOMATIC :: Q
      REAL   (dp ) :: QI
      INTEGER(i4b) :: I,II
!local add MATINV values
      REAL   (dp ) :: DIV,SUMM
      INTEGER(i4b) :: IM1,IP1,J,JJ,JM1,JP1,K,KO,L
! 
      DO 10 I=1,N
       QI=0.0_dp
      DO 20 II=1,N
20     QI=MAX(QI,ABS(A(II,I)))
       Q(I)=QI
      DO 10 II=1,N
10     A(II,I)=A(II,I)/QI

! =====================================================================
!     SUBROUTINE/CALL MATINV(A,N,NDIM) - inlined!
!
!-----matinv executes matrix inversion by lu decomposition
!-----inversion is accomplished in place,
!-----and the original matrix is replaced by its inverse
!-----note that n must be smaller than or equal to ndim
!
       ILOOP1: DO I = 2,N
          IM1 = I - 1
          DO J = 1,IM1
               JM1 = J - 1
               IF (A(J,J).EQ.0.0_dp) A(J,J) = 1.0E-20_dp
               DIV = A(J,J)
               SUMM = 0.0_dp
               IF (JM1.LT.1) GO TO 70

               DO L = 1,JM1
                    SUMM = SUMM + A(I,L)*A(L,J)
               END DO

        70     CONTINUE
               A(I,J) = (A(I,J)-SUMM)/DIV
          END DO

          DO J = I,N
               SUMM = 0.0_dp
               DO L = 1,IM1
                    SUMM = SUMM + A(I,L)*A(L,J)
               END DO
               A(I,J) = A(I,J) - SUMM
          END DO
       END DO ILOOP1

       IILOOP1: DO II = 2,N
          I = N + 2 - II
          IM1 = I - 1
          IF (IM1.LT.1) CYCLE
          DO JJ = 1,IM1
               J = I - JJ
               JP1 = J + 1
               SUMM = 0.0_dp
               IF (JP1.GT.IM1) GO TO 80
               DO K = JP1,IM1
                    SUMM = SUMM + A(I,K)*A(K,J)
               END DO
        80     CONTINUE
               A(I,J) = -A(I,J) - SUMM
          END DO
       END DO IILOOP1

       IILOOP2: DO II = 1,N
          I = N + 1 - II
          IF (A(I,I).EQ.0.0_dp) A(I,I) = 1.0E-20_dp
          DIV = A(I,I)
          IP1 = I + 1
          IF (IP1.GT.N) GO TO 130
          DO JJ = IP1,N
               J = N + IP1 - JJ
               SUMM = 0.0_dp
               DO K = IP1,J
                    SUMM = SUMM + A(I,K)*A(K,J)
               END DO
               A(I,J) = -SUMM/DIV
          END DO
       130 CONTINUE
          A(I,I) = 1._dp/A(I,I)
       END DO IILOOP2

       ILOOP2: DO I = 1,N
          DO J = 1,N
               KO = MAX0(I,J)
               IF (KO.EQ.J) GO TO 170
               SUMM = 0.0_dp

       150     CONTINUE

               DO K = KO,N
                    SUMM = SUMM + A(I,K)*A(K,J)
               END DO
               GO TO 180

       170     CONTINUE

               SUMM = A(I,KO)
               IF (KO.EQ.N) GO TO 180
               KO = KO + 1
               GO TO 150

       180     CONTINUE
               A(I,J) = SUMM

          END DO
       END DO ILOOP2
! =====================================================================
      DO 30 II=1,N
       QI=Q(II)
      DO 30 I=1,N
30     A(II,I)=A(II,I)/QI
!
   RETURN
!------------------------------------------------------------------------------
  end subroutine Simple_mat_inv
!------------------------------------------------------------------------------
!
end module M_Solve_lin_eq_systems
!==============================================================================
!              ***** END MODULE M_Solve_lin_eq_systems *****
!==============================================================================
!
!==============================================================================
!               ***** MODULE M_Gaunt_factor_freefree *****
!==============================================================================
Module M_Gaunt_factor_freefree
!
!--- routines representing functions for the GAUNT-FACTORS
!     GAUNTFFF
!     GFREEAP
!--- Used by MULTI_TRANSPORT and GREYTEMP
!--- 3.05 - ADI

   USE M_type_definitions
   USE M_global_settings               , ONLY: KIS,NUMNT
   USE M_physical_constants            , ONLY: CLIGHT,HPLDRYDI,XLRYDDK
!
   IMPLICIT NONE
!input
    LOGICAL, SAVE :: OPTGFINIT = .TRUE.
!output
    REAL   (dp ), DIMENSION (KIS,NUMNT) :: GFFREE
!local
    REAL   (sp ), DIMENSION (  KIS,NUMNT), PRIVATE, SAVE :: CONX = 0._sp
    REAL   (sp ), DIMENSION (8,KIS,NUMNT), PRIVATE, SAVE :: CX   = 0._sp
!
CONTAINS
!------------------------------------------------------------------------------
  subroutine GAUNTFFF (FREQCM1,TEFFG,NTZ)
!
!-----GENERATES THERMALLY AVERAGED FREE-FREE NON-RELATIVISTIC GAUNT
!-----FACTOR FOR A HYDROGENIC ION OF CHARGE Z.
!     M 3.05 - ADI
!
!input
    REAL   (dp ), DIMENSION (NUMNT)     :: TEFFG
    REAL   (dp )                        :: FREQCM1
    INTEGER(i4b)                        :: NTZ
!output
!   GFFREE
!local
    REAL   (sp ), DIMENSION (NUMNT)     :: XLRKT,CONLR
    REAL   (sp ), DIMENSION (KIS)       :: Z2,ZF
    REAL   (sp ), DIMENSION (KIS,NUMNT) :: TXG
    INTEGER(i2b)                        :: I,NT
    REAL   (sp )                        :: CLF,XLF,FNM1
!-------------------------------------------------------------------------
     ! INITIALIZATION OF THE NUE-INDEPENDENT PART OF THE f-f GAUNT-FACTOR
     ! INITIALIZE CON AND C IF TEFFG HAS CHANGED
       IF (OPTGFINIT) THEN
           DO NT=1,INT(NTZ,i2b)
          ! XLRXT IS LOG(RYD/KT)
            XLRKT(NT)= LOG10(REAL(XLRYDDK,sp)/REAL(TEFFG(NT),sp))
            CONLR(NT)= 0.72727273_sp*XLRKT(NT)+0.90909091_sp
            XLRKT(NT)= 0.66666667_sp*XLRKT(NT)
           END DO
           DO I = 1,INT(KIS,i2b)
            ZF(I) = 0.66666667_sp*2.0_sp*LOG10(REAL(I,sp))
           END DO
         DO I = 1,INT(KIS,i2b)
           DO NT=1,INT(NTZ,i2b)
            TXG(I,NT)= ZF(I)+XLRKT(NT)
           END DO
         END DO
         DO I = 1,INT(KIS,i2b)
           DO NT=1,INT(NTZ,i2b)
             IF(ABS(TXG(I,NT)) > 2.0) THEN
              CONX(I,NT)=100._sp
             ELSE
              CONX(I,NT)= CONLR(NT)
             END IF
           END DO
         END DO
         DO I = 1,INT(KIS,i2b)
           DO NT=1,INT(NTZ,i2b)
             IF(CONX(I,NT) /= 100._sp) THEN
              CALL GFFINIT(TXG(I,NT), CX(1,I,NT))
!             ============
             END IF
           END DO
         END DO
        OPTGFINIT=.FALSE.
       END IF
!
!-------------------------------------------------------------
      CLF  = REAL(FREQCM1,sp) * REAL(CLIGHT,sp)
!     HPLDRYDI=h[erg/s]/Ryd[erg]=1/(c*Ryd[cm^-1])=3.0397E-16
!     1 Ryd = 2.18E-11 [erg] ; 1 Ryd = 1.09737E+05 [cm^-1]
      XLF  = 0.72727273_sp * LOG10( CLF * REAL(HPLDRYDI,sp) )
      FNM1 = REAL(FREQCM1,sp) * 1.E-04_sp
!-------------------------------------------------------------
          DO NT=1,INT(NTZ,i2b)
           XLRKT(NT)= 9.77_sp+1.26_sp*LOG10(REAL(TEFFG(NT),sp)**1.5/CLF)
           CONLR(NT)= 5.0404E3_sp/REAL(TEFFG(NT),sp)
          END DO
          DO I = 1,INT(KIS,i2b)
           Z2(I) = REAL(I**2,sp)
           ZF(I) = 1.26_sp*LOG10(REAL(I,sp))
          END DO
         DO I = 1,INT(KIS,i2b)
           DO NT=1,INT(NTZ,i2b)
            TXG(I,NT)= XLF+CONX(I,NT)
           END DO
         END DO
        DO I = 1,INT(KIS,i2b)
!!DIR$    IVDEP
         DO NT=1,INT(NTZ,i2b)
!----------------
           IF (ABS(TXG(I,NT)) > 2.0_sp .or. CONX(I,NT) == 100._sp) THEN
            IF (CLF < 5.E12_sp) THEN
           ! LAMERS/WATERS/84/A+A136,37
             GFFREE(I,NT)= REAL(XLRKT(NT) - ZF(I),dp)
!            ==============
            ELSE
             GFFREE(I,NT)= REAL(GFREE(FNM1,CONLR(NT),Z2(I)),dp)
!            ==============
            END IF
           ELSE
             GFFREE(I,NT)= REAL(GFF(TXG(I,NT),Z2(I),CX(1,I,NT)),dp)
!            ==============
           END IF
!----------------
         END DO
        END DO
!
      RETURN
!-------------------------------------------------------------------------
!
  CONTAINS
!------------------------------------------------------------------------------
    subroutine GFFINIT (TX,C)
!
!---INITIALIZATION OF THE NUE-INDEPENDENT PART OF THE f-f GAUNT-FACTOR
!---M 4.91/3.05 ADI
!
!input
    REAL   (sp )               :: TX
!output
    REAL   (sp ), DIMENSION(8) :: C
!local
    REAL   (sp ), DIMENSION(11)         :: B
    REAL   (sp ), DIMENSION(8,11), SAVE :: D
!
    INTEGER(i2b) :: IR,I,J
!
    DATA ((D(I,J),I=1,8),J=1,6)/ &
     8.986940175E+00_sp, -4.009515855E+00_sp,  8.808871266E-01_sp, &
     2.640245111E-02_sp, -4.580645915E-02_sp, -3.568055702E-03_sp, &
     2.827798067E-03_sp,  3.365860195E-04_sp, -8.006936989E-01_sp, &
     9.466021705E-01_sp,  9.043402532E-02_sp, -9.608451450E-02_sp, &
    -1.885629865E-02_sp,  1.050313890E-02_sp,  2.800889961E-03_sp, &
    -1.078209202E-03_sp, -3.781305103E-01_sp,  1.102726332E-01_sp, &
    -1.543619180E-02_sp,  8.310561114E-03_sp,  2.179620525E-02_sp, &
     4.259726289E-03_sp, -4.181588794E-03_sp, -1.770208330E-03_sp, &
     1.877213132E-02_sp, -1.004885705E-01_sp, -5.483366378E-02_sp, &
    -4.520154409E-03_sp,  8.366530426E-03_sp,  3.700273930E-03_sp, &
     6.889320423E-04_sp,  9.460313195E-05_sp,  7.300158392E-02_sp, &
     3.576785497E-03_sp, -4.545307025E-03_sp, -1.017965604E-02_sp, &
    -9.530211924E-03_sp, -3.450186162E-03_sp,  1.040482914E-03_sp, &
     1.407073544E-03_sp, -1.744671550E-03_sp,  2.864013856E-02_sp, &
     1.903394837E-02_sp,  7.091074494E-03_sp, -9.668371391E-04_sp, &
    -2.999107465E-03_sp, -1.820642230E-03_sp, -3.874082085E-04_sp/
    DATA ((D(I,J),I=1,8),J=7,11)/ &
    -1.707268366E-02_sp, -4.694254776E-03_sp,  1.311691517E-03_sp, &
     5.316703136E-03_sp,  5.178193095E-03_sp,  2.451228935E-03_sp, &
    -2.277321615E-05_sp, -8.182359057E-04_sp,  2.567331664E-04_sp, &
    -9.155339970E-03_sp, -6.997479192E-03_sp, -3.571518641E-03_sp, &
    -2.096101038E-04_sp,  1.553822487E-03_sp,  1.509584686E-03_sp, &
     6.212627837E-04_sp,  4.098322531E-03_sp,  1.635218463E-03_sp, &
    -5.918883504E-04_sp, -2.333091048E-03_sp, -2.484138313E-03_sp, &
    -1.359996060E-03_sp, -5.371426147E-05_sp,  5.553549563E-04_sp, &
     3.837562402E-05_sp,  2.938325230E-03_sp,  2.393747064E-03_sp, &
     1.328839809E-03_sp,  9.135013312E-05_sp, -7.137252303E-04_sp, &
    -7.656848158E-04_sp, -3.504683798E-04_sp, -8.491991820E-04_sp, &
    -3.615327726E-04_sp,  3.148015257E-04_sp,  8.909207650E-04_sp, &
     9.869737522E-04_sp,  6.134671184E-04_sp,  1.068883394E-04_sp, &
    -2.046080100E-04_sp/
!-------------------------------------------------------------------------
       DO J=1,8
           IR=9
           B(11)=D(J,11)
           B(10)=TX*B(11)+D(J,10)
             DO I=1,9
                B(IR)=TX*B(IR+1)-B(IR+2)+D(J,IR)
                IR   =IR-1
             END DO
           C(J)=0.25_sp*(B(1)-B(3))
       END DO
!
      RETURN
!-------------------------------------------------------------------------
    end subroutine GFFINIT
!------------------------------------------------------------------------------
    function GFF (TX,Z_2,C)
!
!-----GENERATES THERMALLY AVERAGED FREE-FREE NON-RELATIVISTIC GAUNT
!-----FACTOR FOR A HYDROGENIC ION OF CHARGE Z, WITH A MAXIMUM RELATIVE
!-----ERROR OF .OO7, (RMS FITTING ERROR = .001) FOR TEMPERATURE AND
!-----FREQUENCY IN INTERVALS:
!                     10**-4 LE U LE 10**1.5,
!                     10**-3 LE GAMS LE 10**3,
!-----WHERE
!             GAMS = Z**2*RYD/K*T = Z**2*1.579E5/T
!                                   MAX          T=1.6E8K
!                                   MIN          T=5000.K
!
!             U = H*NU/K*T = 1.44*FR(CM-1)/T
!                            MAX  .17      T=2400.K
!                                 10.      T=140000.K
!                            MIN  3.2E7    T=1.4E6K
!                                 1.E6     T=45000.K
!                                 1.1E5    T=5000.K
!
!-----TO OBTAIN THE
!-----STATED ACCURACY, THE FULL NUMBER OF SIGNIFICANT FIGURES IN THE
!-----COEFFICIENTS MUST BE RETAINED.
!
!-----THIS SUBROUTINE USES A TWO-DIMENSIONAL CHEBYSHEV EXPANSION
!-----COMPUTED FROM EXPRESSIONS GIVEN BY KARZAS AND LATTER (AP.J.
!-----SUPPL., V.6, P.167, 1961) AUGMENTED BY VARIOUS LIMITING FORMS
!-----OF ENERGY-SPECIFIC GAUNT FACTORS.
!     D.G.HUMMER, JILA, MAY 1987.
!
!       SUBROUTINE ARGUMENTS ARE:
!         Z    =  NUCLEAR CHARGE
!         T    =  TEMPERATURE IN DEGREES KELVIN
!         XLF  =  CONTAINING INPUT VALUES OF LOG10(H*NU/RYD)
!         IFLAG=  1 IF GAMS IS OUT OF RANGE; 2 IF U IS OUT OF RANGE
!
!-----M 4.91/11.92/2.93/4.94/3.05 ADI
!     =======================================================
!
!input
      REAL   (sp )                :: TX,Z_2
      REAL   (sp ), DIMENSION (8) :: C
!output
      REAL   (sp ) :: GFF
!local
      REAL   (sp ) :: T1,T2,T3,T4,T5,T6,T7
!-------------------------------------------------------------------------
         T1=TX
         T2=T1*TX
         T3=T2*TX
         T4=T3*TX
         T5=T4*TX
         T6=T5*TX
         T7=T6*TX
!        ======
         GFF = &
              ( (T7-8._sp*T5+18._sp*T3-9._sp*T1)*C(8) &
               +(T6-6._sp*T4+ 9._sp*T2-2._sp)*C(7)    &
               +(T5-5._sp*T3+ 5._sp*T1)*C(6)          &
               +(T4-4._sp*T2+ 2._sp)*C(5)             &
               +(T3-3._sp*T1)*C(4)                    &
               +(T2-2._sp)*C(3)                       &
               +(T1)*C(2)                             &
               +C(1) ) * Z_2
!        ======
      RETURN
!-------------------------------------------------------------------------
    end function GFF
!------------------------------------------------------------------------------
    function GFREE(FN_M1,CDT,Z_2)
!
!-----FREE-FREE-GAUNT-FACTORS
!-----M 11.92/2.93/4.94/3.05 ADI
!     ==========================
!
!input
      REAL   (sp ) :: FN_M1,CDT,Z_2
!output
      REAL   (sp ) :: GFREE
!local
      REAL   (sp ) :: TH,X,C1,C2,C3,C4
!-------------------------------------------------------------------------
      TH = MAX( CDT   * Z_2 , 4.0E-2_sp)
      X  = MAX( FN_M1 / Z_2 ,    0.2_sp)
       IF(X < 1.0) THEN
!      ======
       GFREE=(1.0823_sp+2.98E-2_sp/TH)+(6.7E-3_sp+1.12E-2_sp/TH)/X
!      ======
       ELSE
        C1=(3.9999187E-3_sp-7.8622889E-5_sp/TH)/TH+1.070192_sp
        C2=(6.4628601E-2_sp-6.1953813E-4_sp/TH)/TH+2.6061249E-1_sp
        C3=(3.7542343E-2_sp+1.3983474E-5_sp/TH)/TH+5.7917786E-1_sp
        C4= 3.4169006E-1_sp+1.1852264E-2_sp/TH
!      ======
       GFREE=((C4/X-C3)/X+C2)/X+C1
!      ======
       END IF
!
      RETURN
!-------------------------------------------------------------------------
    end function GFREE
!------------------------------------------------------------------------------
  end subroutine GAUNTFFF
!------------------------------------------------------------------------------
  function GFREEAP (A)
!
!!DIR$ INLINE ALWAYS GFREEAP
!
!----APPROXIMATIVE SOLUTION OF GAUNT-FACTORS 20/07/90 M 3.05 ADI
!----CALLED BY 'ETASHM'
!
!input
      REAL   (dp ) :: A
!output
      REAL   (dp ) :: GFREEAP
!local
!-------------------------------------------------------------------------
!
      IF (A.GT.1.) THEN
      GFREEAP=1._dp
      ELSE
      GFREEAP=.551329_dp*LOG(2.35_dp/A)
      END IF
!
      RETURN
!-------------------------------------------------------------------------
  end function GFREEAP
!------------------------------------------------------------------------------
!
end module M_Gaunt_factor_freefree
!==============================================================================
!               ***** END MODULE M_Gaunt_factor_freefree *****
!==============================================================================
!
!==============================================================================
!                 ***** MODULE M_Mathematical_functions *****
!==============================================================================
Module M_Mathematical_functions
!--- routines representing mathematical functions
!     EXPINT
!     ERF7
!     VOIGT
!     VOIGTSTEP
!     SORTH
!     INDSORT
!     CINTF
!--- 3.05 - ADI
!
   USE M_type_definitions
!
   IMPLICIT NONE
!
CONTAINS
!------------------------------------------------------------------------------
  function EXPINT(X,N)
!
!!DIR$ INLINE ALWAYS EXPINT
!
!-----COMPUTATION OF AN EXPONENTIAL INTEGRAL
!-----M 11.91/3.05 ADI   CALLED BY 'RJIWC'
!
!input
     INTEGER(i4b) :: N
     REAL    (dp) :: X
!output
     REAL    (dp) :: EXPINT
!local
     INTEGER(i4b) :: K
     REAL    (dp) :: RN,XN,XNQ,SUM,XK
!-------------------------------------------------------------------------
      IF(X.LT.0.) STOP 'ARGUMENT OF EXPONENTIAL INTEGRAL LOWER ZERO!'
      IF(X.EQ.0..AND.N.LT.2) STOP 'ARGUMENT OF EXP.INT. EQUAL ZERO!'
!
      IF(X.GT.1000.) THEN
        STOP 'EXPINT CAN`T HELP YOU'
      ENDIF
!
      IF(N.EQ.0) THEN
        EXPINT=EXP(-X)/X
        RETURN
      ENDIF
!
      IF(N.EQ.1) THEN
        EXPINT=EXPINT1(X)
        RETURN
      ENDIF
!
      IF(X.EQ.0.) THEN
        EXPINT=1._dp/REAL(N-1,dp)
        RETURN
      ENDIF
!
      IF(N.GE.20.OR.(N.GE.8.AND.X.GT.5.))THEN
        RN=REAL(N,dp)
        XN=X + RN
        XNQ=XN*XN
        SUM=1._dp + RN/XNQ + RN*(        RN - 2._dp*X           )/XNQ**2 + &
                             RN*(6._dp*X**2 - 8._dp*RN*X + RN**2)/XNQ**3
        EXPINT=EXP(-X)/XN*SUM
        RETURN
      ENDIF
!
      SUM=0._dp
      XK =1._dp
        DO K=0,N-2
         IF (K.GT.0) XK=XK*(-X)
         SUM=SUM + XK*REAL(IFAKUL(N-K-2),dp)
        END DO
      SUM=SUM*EXP(-X)+(-1._dp)**(N-1)*X**(N-1)*EXPINT1(X)
      EXPINT=SUM/REAL(IFAKUL(N-1),dp)
!
      RETURN
!-------------------------------------------------------------------------
  CONTAINS
!------------------------------------------------------------------------------
    function EXPINT1(X)
!
!-----COMPUTATION OF THE FIRST EXPONENTIAL INTEGRAL
!
!input
     REAL    (dp) :: X
!output
     REAL    (dp) :: EXPINT1
!local
     REAL    (dp) :: SUM1,SUM2
     REAL    (dp), SAVE :: &
                   A0=-.57721566_dp,A1= .99999193_dp,A2=-.24991055_dp, &
                   A3= .05519968_dp,A4=-.00976004_dp,A5= .00107857_dp, &
                   B1= 8.5733287401_dp,B2=18.0590169730_dp, &
                   B3= 8.6347608925_dp,B4=  .2677737343_dp, &
                   C1= 9.5733223454_dp,C2=25.6329561486_dp, &
                   C3=21.0996530827_dp,C4= 3.9584969228_dp
!-------------------------------------------------------------------------
        IF(X.LE.1.) THEN
          EXPINT1=A0+X*(A1+X*(A2+X*(A3+X*(A4+X*A5))))-LOG(X)
          RETURN
        ENDIF
      SUM1= B4+X*(B3+X*(B2+X*(B1+X)))
      SUM2= C4+X*(C3+X*(C2+X*(C1+X)))
      EXPINT1=EXP(-X)/X*SUM1/SUM2
!
     RETURN
!-------------------------------------------------------------------------
    end function EXPINT1
!-------------------------------------------------------------------------
!-------------------------------------------------------------------------
    function IFAKUL(M)
!
!!DIR$ INLINE ALWAYS IFAKUL
!
!-----COMPUTATION OF A FACTORIAL
!
!input
     INTEGER(i4b) :: M
!output
     INTEGER(i4b) :: IFAKUL
!local
     INTEGER(i4b) :: J,I
!----------------------------------------------------------
      IF(M.LT.0) STOP 'ARGUMENT OF FAKULTAET NEGATIVE!'
        IF(M.LE.1) THEN
           IFAKUL=1
           RETURN
        ENDIF
      J=1
        DO 10 I=1,M
10      J=J*I
      IFAKUL=J
     RETURN
!----------------------------------------------------------
    end function IFAKUL
!------------------------------------------------------------------------------
  end function EXPINT
!------------------------------------------------------------------------------
  function EXPM1(X)
!
!----NUMERICALLY STABLE COMPUTATION OF EXP(X)-1.
!----FOR SMALL ABOLUTE VALUES OF X, EXP(X) IS VERY COLSE TO X.
!----AS A RESULT, WHEN COMPUTING EXP(X)-1, a "SMALL DIFFERNECE OF
!----LARGE NUMBERS" PROBLEM OCCURS, WHICH CAN, FOR EXAMPLE, LEAD
!----TO A DIVISION-BY-ZERO PROBLEM.
!----THE THRESHOLD OF 1E-3 HAS BEEN FOUND BY COMPARING BOTH THE FORTRAN
!----EXPRESSION "EXP(X)-1" AND THE TAYLOR SERIES UP the THE X**4 TERM
!----WITH THE EXPM1 FUNCTION CONTAINED IN THE C STANDARD LIBRARY  
!
!======================================================================
!  CALLED BY BNUE2 formal_sup.f90 -- M 16.5 JWEBER
!======================================================================
!
    IMPLICIT NONE
!
! input
    REAL   (dp) :: X
! output 
    REAL   (dp) :: EXPM1
    IF (ABS(X) .LT. (1E-3_dp)) THEN
       EXPM1=X*(1+(1._dp/2._dp)*X*(1+(1._dp/3._dp)*X*(1+(1._dp/4._dp)*X**2)))
    ELSE
       EXPM1=EXP(X)-1._dp
    END IF
!-------------------------------------------------------------------------
  end function EXPM1
!-------------------------------------------------------------------------
  function ERF7(X)
!
!!DIR$ INLINE ALWAYS ERF
!
! ERRORFUNCTION WITH 7-sf DIGITS
!
!input
      REAL   (dp )  :: X
!output
      REAL   (dp )  :: ERF7
!local
      REAL   (dp )  :: SX,T
!---------------------------------------------------------------------
      SX= SIGN(1._dp,X)
      T =  1._dp / ( 1._dp + .3275911_dp * ABS(X) )
      ERF7 = SX  * ( 1._dp - (((((1.061405429_dp*T -1.453152027_dp)*T  &
                           +     1.421413741_dp)*T - .284496736_dp)*T  &
                           +      .254829592_dp)*T) * EXP(-X*X) )
      RETURN
!---------------------------------------------------------------------
  end function ERF7
!------------------------------------------------------------------------------
  function VOIGT(X,Y,EPSP,EPSM,YBORD)
!
!!DIR$ INLINE ALWAYS VOIGT
!
!-----FIT-FORMULA FOR VOIGT-FUNCTION H(X,Y) --> PHI=H/DNUED/WPI
!-----Y > 1.E-3: HUI ET AL 1978, J. QUANT. SPEC. RAD. TRANSF. 19, 509
!-----Y < 1.E-3: FS
!
!input
      REAL(dp)    :: X,Y,EPSP,EPSM,YBORD
!output
      REAL(dp)    :: VOIGT
!local
      REAL(dp), PARAMETER :: &
                 A0=37.24429446739879_dp, A1=57.90331938807185_dp, &
                 A2=43.16280063072749_dp, A3=18.64649990312317_dp, &
                 A4= 4.67506018267650_dp, A5= 0.56418958297228_dp, &
                 B0=    37.2442945086_dp, B1=    99.9290005933_dp, &
                 B2=   118.6763981260_dp, B3=    80.6459493922_dp, &
                 B4=    33.5501020941_dp, B5=     8.2863279156_dp, &
                 WPIHALB   =0.8862269_dp, EPS=           1.E-6_dp
      COMPLEX(dp) :: ZC, ZD
!-------------------------------------------------------------------------
       IF (Y.GT.YBORD) THEN !1.E-2
        ZC = CMPLX(Y,-X,dp)
        VOIGT = REAL( (((((A5*ZC+A4)*ZC+A3)*ZC+A2)*ZC+A1)*ZC+A0) / &
                  ((((((ZC+B5)*ZC+B4)*ZC+B3)*ZC+B2)*ZC+B1)*ZC+B0),dp )
       ELSE
        ZC = CMPLX(Y,-X-EPSP,dp)
        ZD = CMPLX(Y,-X+EPSM,dp)
        VOIGT=(Y*(AIMAG( (((((A5*ZD+A4)*ZD+A3)*ZD+A2)*ZD+A1)*ZD+A0) /   &
                     ((((((ZD+B5)*ZD+B4)*ZD+B3)*ZD+B2)*ZD+B1)*ZD+B0) )- &
                   AIMAG( (((((A5*ZC+A4)*ZC+A3)*ZC+A2)*ZC+A1)*ZC+A0) /  &
                     ((((((ZC+B5)*ZC+B4)*ZC+B3)*ZC+B2)*ZC+B1)*ZC+B0) )) &
                     +WPIHALB*(ERF7(X+EPSP)-ERF7(X-EPSM)))/(EPSP+EPSM)
       ENDIF
!
      RETURN
!-------------------------------------------------------------------------
  end function VOIGT
!------------------------------------------------------------------------------
  function VOIGTSTEP(X,Y,DX)
!
!!DIR$ INLINE ALWAYS VOIGTSTEP
!
!-----CONVOLUTION OF VOIGT-PROFILE AND STEP-FUNCTION:
!-----VOIGTSTEP(x) = (f*g)(x)  with    f(x)=H(x,y)
!-----                         and     g(x)=theta(x)*theta(dx-x)/dx
!-----ACCURATE FOR 0 < Y < 0.1
!-----CALLED BY BLTBLOCK -- 5.95 M 3.05 ADI
!
!input
      REAL(dp)    :: X,Y,DX
!output
      REAL(dp)    :: VOIGTSTEP
!local
      REAL(dp)    :: DERF
      REAL(dp), PARAMETER :: WPIHALB   =0.8862269_dp
!---------------------------------------------------------------------
      IF (Y.LT.1.E-8) THEN
        DERF=WPIHALB*(ERF7(X)-ERF7(X-DX))
        VOIGTSTEP=DERF/DX
      ELSE
        VOIGTSTEP=VOIGT(X,Y,0._dp,DX,0.3_dp)
      ENDIF
       IF (DX.GT.0.3_dp .and. Y.GT.0.3_dp) THEN
        DERF=WPIHALB*(ERF7(X)-ERF7(X-DX))
        VOIGTSTEP=MAX( VOIGTSTEP, DERF/DX )
       ENDIF
!
      RETURN
!---------------------------------------------------------------------
  end function VOIGTSTEP
!------------------------------------------------------------------------------
  subroutine SORTH ( N, RA, ASCND, RB, IC, IMI )
!
!
!-----SORTS AN ARRAY 'RA' OF LENGTH 'N' INTO SELECTED NUMERICAL ORDER
!-----USING THE HEAPSORT ALGORITHM, WHILE MAKING THE CORRESPONDING
!-----REARRANGEMENT OF THE ARRAYS 'RB'/'IC'.'ASCND' IS INPUT AND SELECTS
!-----SORTING TO ASCENDING ORDER IF .TRUE. AND TO DESCENDING ORDER
!-----OTHERWISE.
!-----M 12.93/3.05 ADI    CALLED BY 'INDATR','CALETAS/SETSHFR'
!
!input
      LOGICAL       :: ASCND
      INTEGER(i4b)  :: N
!inoutput
      REAL   (dp ), DIMENSION (N)     :: RA
      REAL   (dp ), DIMENSION (N), OPTIONAL    :: RB
      INTEGER(i4b), DIMENSION (N), OPTIONAL    :: IC
!output
      INTEGER(i4b), DIMENSION (N), OPTIONAL, TARGET :: IMI
!local
      INTEGER(i4b)  :: L,IR,I,J,IRC
      REAL   (dp )  :: RRA
      INTEGER(i4b), DIMENSION (:), POINTER     :: IM
      INTEGER(i4b), DIMENSION (:), ALLOCATABLE :: ID
      REAL   (dp ), DIMENSION (:), ALLOCATABLE :: RD
!---------------------------------------------------------------------
      IF (PRESENT(IMI)) THEN
       IM => IMI
      ELSE
       ALLOCATE (IM(N))
      END IF
       DO I = 1,N
         IM(I) = I
       END DO
!--------------------------------
      L  = N/2 + 1
      IR = N
   10 CONTINUE
!--------------------------------
        IF ( L .GT. 1 ) THEN
          L = L - 1
          RRA = RA(L)
          IRC = IM(L)
        ELSE
          RRA = RA(IR)
          IRC = IM(IR)
          RA(IR) = RA(1)
          IM(IR) = IM(1)
          IR = IR - 1
          IF ( IR .EQ. 1 ) THEN
            RA(1) = RRA
            IM(1) = IRC
            GOTO 30
          END IF
        END IF
!--------------------------------
        I = L
        J = L + L
!-------------------------------------------------
   20   IF ( J .LE. IR ) THEN
          IF ( J .LT. IR ) THEN
            IF ( RA(J) .LT. RA(J+1) )   J = J + 1
          END IF
          IF ( RRA .LT. RA(J) ) THEN
            RA(I) = RA(J)
            IM(I) = IM(J)
            I = J
            J = J + J
          ELSE
            J = IR + 1
          END IF
        GOTO 20
        END IF
!-------------------------------------------------
        RA(I) = RRA
        IM(I) = IRC
      GOTO 10
!
   30 CONTINUE
!-----------------------------
      IF ( .NOT. ASCND ) THEN
        DO I = 1,N/2
          J = N - I + 1
          RRA = RA(I)
          IRC = IM(I)
          RA(I) = RA(J)
          IM(I) = IM(J)
          RA(J) = RRA
          IM(J) = IRC
        END DO
      END IF
!-----------------------------
      IF (PRESENT(RB)) THEN
       ALLOCATE (RD(N))
        DO I = 1,N
          RD(I) = RB(I)
        END DO
        DO I = 1,N
          RB(I) = RD(IM(I))
        END DO
       DEALLOCATE (RD)
      END IF
      IF (PRESENT(IC)) THEN
       ALLOCATE (ID(N))
        DO I = 1,N
          ID(I) = IC(I)
        END DO
        DO I = 1,N
          IC(I) = ID(IM(I))
        END DO
       DEALLOCATE (ID)
      END IF
!-----------------------------
      IF (PRESENT(IMI)) THEN
       NULLIFY (IM)
      ELSE
       DEALLOCATE (IM)
      END IF
!-----------------------------
!
      RETURN
!---------------------------------------------------------------------
  end subroutine SORTH
!------------------------------------------------------------------------------
  subroutine INDSORT(M,MM,NN,ARRIN,INDX)
!
!!DIR$ INLINE ALWAYS INDSORT
!
!-----INDEXES AN ARRAY ARRIN OF LENGTH N, I.E. OUTPUTS THE ARRAY INDX
!-----SUCHTHAT ARRIN(INDX(J)) IS IN ASCENDING ORDER FOR J=1..N.
!-----THE INPUT QUANTITIES N AND ARRIN ARE NOT CHANGED!
!-----NOT ALL VALUES OF STREN ARE SORTED, JUST THE RANGE ENTERED.
!        INPUTS:  M-BEGINNING POSITION OF ARRIN TO SORT
!                MM-ENDING    POSITION OF ARRIN TO SORT
!-----M 2.94/3.05  ADI  - CALLED BY 'SELOPA','L200'
!
!input
      INTEGER(i4b)  :: M,MM,NN
      REAL   (dp ), DIMENSION (NN)            :: ARRIN
!output
      INTEGER(i4b), DIMENSION (NN)            :: INDX
!local
      REAL   (dp )  :: Q
      INTEGER(i4b)  :: N,MN,I,J,L,IR,INDXT
!---------------------------------------------------------------------
      !PRINT*, 'M', M
      !PRINT*, 'MM', MM
      !PRINT*, 'NN', NN
      !PRINT*, 'ARRIN', ARRIN
      !PRINT*, 'INDX', INDX

      MN=M-1 ; N=MM-MN
!
       DO J=1,N
         INDX(J)=J
       END DO
!
      L=N/2+1
      IR=N
!
 10   IF (L.GT.1) THEN
         L=L-1
         INDXT=INDX(L)
         Q=ARRIN(MN+INDXT)
      ELSE
         INDXT=INDX(IR)
         !PRINT*, 'IR', IR
         !PRINT*, 'MN', MN
         !!PRINT*, 'INDX(IR)', INDX(IR)
         !PRINT*, 'INDXT', INDXT
         !PRINT*, 'MN+INDXT', MN+INDXT
         !PRINT*, ''
         Q=ARRIN(MN+INDXT)
         INDX(IR)=INDX(1)
         IR=IR-1
         IF (IR.EQ.1) THEN
            INDX(1)=INDXT
            RETURN
         ENDIF
      ENDIF
      I=L
      J=L+L
!
 20   IF (J.LE.IR) THEN
         IF (J.LT.IR) THEN
            IF (ARRIN(MN+INDX(J)) .LT. ARRIN(MN+INDX(J+1))) J=J+1
         ENDIF
         IF (Q.LT.ARRIN(MN+INDX(J))) THEN
            INDX(I)=INDX(J)
            I=J
            J=J+J
         ELSE
            J=IR+1
         END IF
         GO TO 20
      END IF
      INDX(I)=INDXT
      GO TO 10
!---------------------------------------------------------------------
  end subroutine INDSORT
!------------------------------------------------------------------------------
  function CINTF (A)
!
!!DIR$ INLINE ALWAYS CINTF
!
!----SOLUTION OF AN EXPONENTIAL INTEGRAL 20.07.90 M 3.05 ADI
!----CALLED BY 'HYDTRAN'
!
!input
      REAL   (dp )  :: A
!output
      REAL   (dp )  :: CINTF
!local
      INTEGER(i4b)  :: I,I1,N1,NM,K
      REAL   (dp )  :: X,S,A1,AN
      REAL   (dp ), DIMENSION(2) :: AF,AI
      REAL   (dp ), SAVE :: FU=3.2E7_dp,FL=1.6E6_dp,C=0.5772156649_dp
!------------------------------------------------------------------------------
      AF(1)=A*FU
      AF(2)=A*FL
!
      IF (AF(2).GT.10.) STOP 'THE SHOCK THEORY SHOULD NOT BE USED !'
!
       IF (AF(1).GT.15.) THEN
        AI(1)=-C
        I1=2
       ELSE
        I1=1
       END IF
!------------------------------------------------------------------------------
      DO I=I1,2
       N1=INT(AF(I)+1._dp)
       X =REAL(N1,dp)
        IF (N1.LT.5) THEN
         S=1._dp
        ELSE
         S=4._dp
        END IF
       NM=INT(X*(1._dp+X/S)+(7.6_dp-1.5_dp*LOG10(X)+.4343_dp*X)/LOG10(1._dp+X/S))
       A1   = -AF(I)
       AI(I)= LOG(AF(I))+A1
       AN=A1
        DO K=2,NM
         AN   = AN*A1*(K-1)/K**2
         AI(I)= AI(I)+AN
        END DO
      END DO
!------------------------------------------------------------------------------
!
      CINTF=AI(1)-AI(2)
!     ======
!-----------------------------------------------------------------------
!      PRINT*,' CINTF: NM=',NM,' A(1)=',AI(1),' A(2)=',AI(2), &
!             ' A=',AI(1)-AI(2),' AF(1)=',AF(1),' AF(2)=',AF(2)
!-----------------------------------------------------------------------
      RETURN
!------------------------------------------------------------------------------
  end function CINTF
!------------------------------------------------------------------------------
  !function
  !end function
!------------------------------------------------------------------------------
!
end module M_Mathematical_functions
!==============================================================================
!               ***** END MODULE M_Mathematical_functions *****
!==============================================================================
!
!==============================================================================
!                  ***** MODULE M_Physical_functions *****
!==============================================================================
Module M_Physical_functions
!--- routines representing physical functions
!     EDDPLA
!     TRADFUN
!     HYDREC
!--- 3.05 - ADI
!
   USE M_type_definitions
!
   IMPLICIT NONE
!
CONTAINS
!------------------------------------------------------------------------------
  function EDDPLA (AL,TEMP)
!
!!DIR$ INLINE ALWAYS EDDPLA
!
!-----EDDPLA GIVES THE PLANCK FUNCTION FROM A
!-----BLACK BODY (INDEX IS NUE, AL IN CM). - FILE2 ADI
!
!-----
   USE M_physical_constants             , ONLY: HCK,BFACT
!-----
!input
      REAL(dp)    :: AL,TEMP
!output
      REAL(dp)    :: EDDPLA
!local
      REAL(dp)    :: AUX
!-----------------------------------------------------------------------
!
!     EDDPLA = BFACT/(AL**3*(EXP(HCK/(AL*TEMP))-1._dp))
        AUX = HCK/(AL*TEMP)
      IF      (EXP(AUX).EQ.1._dp) THEN
        EDDPLA = BFACT/(AL**3*AUX)
      ELSE IF (AUX.LT.37._dp)     THEN
        EDDPLA = BFACT/(AL**3*(EXP(AUX)-1._dp))
      ELSE
        EDDPLA = BFACT/(AL**3)*(EXP(-AUX))
      ENDIF
!
      RETURN
!-----------------------------------------------------------------------
  end function EDDPLA
!------------------------------------------------------------------------------
  function TRADFUN (AIC,AJNUEW,FREQ)
!
!!DIR$ INLINE ALWAYS TRADFUN
!
!-----CALCULATES RADIATION TEMPERATURES--FREQ IN (CM)
!-----CALLED BY 'CONTO' ('CONTOUT'),'BFRCF'
!-----M 3.05 ADI
!
!-----
   USE M_physical_constants             , ONLY: HCK,BFACT
!-----
!input
      REAL(dp)    :: AIC,AJNUEW,FREQ
!output
      REAL(dp)    :: TRADFUN
!local
      REAL(dp)    :: FAK3
!---------------------------------------------------------------------
      FAK3 = BFACT / FREQ**3 / ( AIC * AJNUEW )
!
      IF (FAK3 <= 1.E-8_dp) THEN
        TRADFUN = HCK / FREQ / (FAK3 - 0.5_dp * FAK3 * FAK3)
      ELSE
        TRADFUN = HCK / FREQ / LOG (FAK3 + 1._dp)
      END IF
!
      RETURN
!---------------------------------------------------------------------
  end function TRADFUN
!------------------------------------------------------------------------------
  function HYDREC(T,N)
!
! Computes total recombination rates and total recombination cooling
! coefficients for hydrogen by interpolating the table of
! D.G.Hummer, MNRAS 268, 109-112 (1994)
! Arguments: T is temperature
!            N is 1 for alpha_A
!                 2 for beta_A
!                 3 for beta_ff
! Table entries are:
! Log T, Q*alpha_1, Q*alpha_B, Q*beta_1, Q*beta_B, Q*beta_ff, Q*beta_B_tot
! where Q=Sqrt(T) and beta_B_tot=beta_B+beta_ff
! M 3.05 ADI
!
! called by 
!
!input
      REAL   (dp )  :: T
      INTEGER(i4b)  :: N
!output
      REAL   (dp )  :: HYDREC
!local
      INTEGER(i4b)  :: I,J
      REAL   (dp )  :: Q,S,A,X
      REAL   (dp ), DIMENSION(31,7) :: C
!
      DATA ((C(I,J),J=1,7),I=1,31)/ &
 1.0_dp,1.646E-11_dp,9.283E-11_dp,1.646E-11_dp,8.287E-11_dp,1.061E-11_dp,9.348E-11_dp,&
 1.2_dp,1.646E-11_dp,8.823E-11_dp,1.646E-11_dp,7.821E-11_dp,1.068E-11_dp,8.889E-11_dp,&
 1.4_dp,1.646E-11_dp,8.361E-11_dp,1.646E-11_dp,7.356E-11_dp,1.076E-11_dp,8.432E-11_dp,&
 1.6_dp,1.646E-11_dp,7.898E-11_dp,1.646E-11_dp,6.892E-11_dp,1.085E-11_dp,7.977E-11_dp,&
 1.8_dp,1.646E-11_dp,7.435E-11_dp,1.646E-11_dp,6.430E-11_dp,1.095E-11_dp,7.525E-11_dp,&
 2.0_dp,1.646E-11_dp,6.973E-11_dp,1.645E-11_dp,5.971E-11_dp,1.106E-11_dp,7.077E-11_dp,&
 2.2_dp,1.645E-11_dp,6.512E-11_dp,1.644E-11_dp,5.515E-11_dp,1.118E-11_dp,6.633E-11_dp,&
 2.4_dp,1.645E-11_dp,6.054E-11_dp,1.643E-11_dp,5.062E-11_dp,1.132E-11_dp,6.194E-11_dp,&
 2.6_dp,1.644E-11_dp,5.599E-11_dp,1.641E-11_dp,4.614E-11_dp,1.145E-11_dp,5.758E-11_dp,&
 2.8_dp,1.642E-11_dp,5.147E-11_dp,1.638E-11_dp,4.170E-11_dp,1.161E-11_dp,5.332E-11_dp,&
 3.0_dp,1.640E-11_dp,4.700E-11_dp,1.633E-11_dp,3.734E-11_dp,1.181E-11_dp,4.915E-11_dp,&
 3.2_dp,1.636E-11_dp,4.258E-11_dp,1.625E-11_dp,3.306E-11_dp,1.202E-11_dp,4.508E-11_dp,&
 3.4_dp,1.629E-11_dp,3.823E-11_dp,1.613E-11_dp,2.888E-11_dp,1.224E-11_dp,4.112E-11_dp,&
 3.6_dp,1.620E-11_dp,3.397E-11_dp,1.594E-11_dp,2.484E-11_dp,1.248E-11_dp,3.733E-11_dp,&
 3.8_dp,1.605E-11_dp,2.983E-11_dp,1.565E-11_dp,2.098E-11_dp,1.274E-11_dp,3.373E-11_dp,&
 4.0_dp,1.582E-11_dp,2.584E-11_dp,1.522E-11_dp,1.736E-11_dp,1.303E-11_dp,3.039E-11_dp,&
 4.2_dp,1.548E-11_dp,2.204E-11_dp,1.460E-11_dp,1.402E-11_dp,1.335E-11_dp,2.737E-11_dp,&
 4.4_dp,1.499E-11_dp,1.846E-11_dp,1.374E-11_dp,1.103E-11_dp,1.369E-11_dp,2.472E-11_dp,&
 4.6_dp,1.431E-11_dp,1.520E-11_dp,1.260E-11_dp,8.442E-12_dp,1.403E-11_dp,2.247E-11_dp,&
 4.8_dp,1.341E-11_dp,1.226E-11_dp,1.119E-11_dp,6.279E-12_dp,1.434E-11_dp,2.062E-11_dp,&
 5.0_dp,1.227E-11_dp,9.696E-12_dp,9.571E-12_dp,4.539E-12_dp,1.460E-11_dp,1.914E-11_dp,&
 5.2_dp,1.093E-11_dp,7.514E-12_dp,7.844E-12_dp,3.192E-12_dp,1.478E-11_dp,1.797E-11_dp,&
 5.4_dp,9.454E-12_dp,5.710E-12_dp,6.146E-12_dp,2.185E-12_dp,1.485E-11_dp,1.704E-11_dp,&
 5.6_dp,7.920E-12_dp,4.257E-12_dp,4.601E-12_dp,1.458E-12_dp,1.482E-11_dp,1.628E-11_dp,&
 5.8_dp,6.427E-12_dp,3.117E-12_dp,3.295E-12_dp,9.484E-13_dp,1.468E-11_dp,1.563E-11_dp,&
 6.0_dp,5.058E-12_dp,2.244E-12_dp,2.262E-12_dp,6.023E-13_dp,1.444E-11_dp,1.505E-11_dp,&
 6.2_dp,3.866E-12_dp,1.590E-12_dp,1.494E-12_dp,3.738E-13_dp,1.414E-11_dp,1.451E-11_dp,&
 6.4_dp,2.877E-12_dp,1.110E-12_dp,9.520E-13_dp,2.268E-13_dp,1.380E-11_dp,1.402E-11_dp,&
 6.6_dp,2.089E-12_dp,7.642E-13_dp,5.878E-13_dp,1.348E-13_dp,1.344E-11_dp,1.358E-11_dp,&
 6.8_dp,1.485E-12_dp,5.199E-13_dp,3.528E-13_dp,7.859E-14_dp,1.311E-11_dp,1.318E-11_dp,&
 7.0_dp,1.036E-12_dp,3.498E-13_dp,2.066E-13_dp,4.499E-14_dp,1.280E-11_dp,1.285E-11_dp/
!
!--------------------------------------------------------------------------------------
! statement function for table interpolation
      X(A,I,J)=C(I,J)+A*(C(I+1,J)-C(I,J))
!-------------------------------------------------------------------------
!
      Q=SQRT (T)
      S=LOG10(T)
      IF(S.GT.7.0.OR.S.LT.1.0) STOP 'HYDREC: Temperature out of range'
       DO I=1,30
        IF(S.GE.C(I,1).AND.S.LE.C(I+1,1)) GOTO 1
       END DO
      STOP 'HYDREC: Unexpected error'
    1 CONTINUE
      A=(S-C(I,1))/(C(I+1,1)-C(I,1))
       IF    (N.EQ.1) THEN
        HYDREC=(X(A,I,2)+X(A,I,3))/Q
       ELSEIF(N.EQ.2) THEN
        HYDREC=(X(A,I,4)+X(A,I,5))/Q
       ELSEIF(N.EQ.3) THEN
        HYDREC=X(A,I,6)/Q
       ELSE
        STOP 'HYDREC: Bad parameter'
       ENDIF
!
    RETURN
!-------------------------------------------------------------------------
  end function HYDREC
!------------------------------------------------------------------------------
!
  function HELREC(T,N)
!
! Computes total recombination rates and total recombination cooling
! coefficients for hydrogen by interpolating the table of
! D.G.Hummer, MNRAS 268, 109-112 (1994)
! Arguments: T is temperature
!            N is 1 for alpha_A
!                 2 for beta_A
!                 3 for beta_ff
!                 4 for alpha_b
!                 5 for beta_B 
!
! Table entries are:
! Log T, Q*alpha_1, Q*alpha_B, Q*beta_1, Q*beta_B, Q*beta_ff, Q*beta_B_tot
! where Q=Sqrt(T) and beta_B_tot=beta_B+beta_ff
!
! MOD 03-15 WEB
! called by additional_helium_recombination (contained in MATRIXC, file1.f90)
!
    IMPLICIT NONE
!
!input
    REAL   (dp )  :: T
    INTEGER(i4b)  :: N
!output
    REAL   (dp )  :: HELREC
!local
    INTEGER(i4b)  :: I,J
    REAL   (dp )  :: Q,S,A,X
    REAL   (dp ), DIMENSION(18,7) :: C
!
    DATA ((C(I,J),J=1,7),I=1,18)/ &
!        alpha_1       alpha_B       beta_1        beta_b         beta_ff        beta_b_tot
 1.0_dp, 1.569E-11_dp, 9.284E-11_dp, 1.569e-11_dp, 8.347E-11_dp,  1.061E-11_dp,  9.408e-11_dp,&
 1.2_dp, 1.569E-11_dp, 8.847E-11_dp, 1.569E-11_dp, 7.889E-11_dp,  1.068E-11_dp,  8.957e-11_dp,&
 1.4_dp, 1.569E-11_dp, 8.403E-11_dp, 1.569E-11_dp, 7.430E-11_dp,  1.076E-11_dp,  8.506e-11_dp,&
 1.6_dp, 1.599E-11_dp, 7.952E-11_dp, 1.569E-11_dp, 6.971E-11_dp,  1.085E-11_dp,  8.056e-11_dp,&
 1.8_dp, 1.569E-11_dp, 7.499E-11_dp, 1.569E-11_dp, 6.512E-11_dp,  1.095E-11_dp,  7.607e-11_dp,&
 2.0_dp, 1.569E-11_dp, 7.044E-11_dp, 1.569E-11_dp, 6.056E-11_dp,  1.106E-11_dp,  7.162e-11_dp,&
 2.2_dp, 1.569E-11_dp, 6.589E-11_dp, 1.570E-11_dp, 5.603E-11_dp,  1.118E-11_dp,  6.721e-11_dp,&
 2.4_dp, 1.560E-11_dp, 6.136E-11_dp, 1.570E-11_dp, 5.154E-11_dp,  1.132E-11_dp,  6.286e-11_dp,&
 2.6_dp, 1.570E-11_dp, 5.685E-11_dp, 1.571E-11_dp, 4.710E-11_dp,  1.145E-11_dp,  5.855e-11_dp,&
 2.8_dp, 1.571E-11_dp, 5.238E-11_dp, 1.572E-11_dp, 4.274E-11_dp,  1.162E-11_dp,  5.436e-11_dp,&
 3.0_dp, 1.572E-11_dp, 4.797E-11_dp, 1.574E-11_dp, 3.847E-11_dp,  1.181E-11_dp,  5.028e-11_dp,&
 3.2_dp, 1.573E-11_dp, 4.364E-11_dp, 1.578E-11_dp, 3.431E-11_dp,  1.202E-11_dp,  4.633e-11_dp,&
 3.4_dp, 1.576E-11_dp, 3.940E-11_dp, 1.583E-11_dp, 3.031E-11_dp,  1.224E-11_dp,  4.255e-11_dp,&
 3.6_dp, 1.580E-11_dp, 3.528E-11_dp, 1.591E-11_dp, 2.650E-11_dp,  1.248E-11_dp,  3.898e-11_dp,&
 3.8_dp, 1.586E-11_dp, 3.132E-11_dp, 1.602E-11_dp, 2.291E-11_dp,  1.274E-11_dp,  3.565e-11_dp,&
 4.0_dp, 1.595E-11_dp, 2.755E-11_dp, 1.619E-11_dp, 1.960E-11_dp,  1.336E-11_dp,  3.296e-11_dp,&
 4.2_dp, 1.608E-11_dp, 2.401E-11_dp, 1.641E-11_dp, 1.660E-11_dp,  1.336E-11_dp,  2.996e-11_dp,&
 4.4_dp, 1.626E-11_dp, 2.073E-11_dp, 1.670E-11_dp, 1.394E-11_dp,  1.369E-11_dp,  2.763e-11_dp/
!
!--------------------------------------------------------------------------------------
! statement function for table interpolation
    X(A,I,J)=C(I,J)+A*(C(I+1,J)-C(I,J))
!-------------------------------------------------------------------------
!
    Q=SQRT (T)
    S=LOG10(T)
!      
    IF(S .LE. 4.4 .AND. S.GE.1.0) THEN
       DO I=1,17
          IF(S.GE.C(I,1).AND.S.LE.C(I+1,1)) GOTO 1
       END DO
       STOP 'HELREC: Unexpected error'
1      CONTINUE         
       A=(S-C(I,1))/(C(I+1,1)-C(I,1))         
       IF    (N.EQ.1) THEN
          HELREC=(X(A,I,2)+X(A,I,3))/Q
       ELSEIF(N.EQ.2) THEN
          HELREC=(X(A,I,4)+X(A,I,5))/Q
       ELSEIF(N.EQ.3) THEN
          HELREC=X(A,I,6)/Q
       ELSEIF(N .EQ. 4) THEN
          HELREC=X(A,I,3)/Q
       ELSEIF(N .EQ. 5) THEN
          HELREC=X(A,I,5)/Q
       ELSE
          STOP 'HYDREC: Bad parameter'
       ENDIF
!
!     The values are tabulated for tebpertures with LOG10(T) < 4.4, i.e. T < 25118 K
!     For higer temperatures, it is assumed that the ratio between the recombination
!     rate coefficients of hydrogen and helium are the same as for log(T)=4.4.
!     This should be changed when better values are available.
!
      ELSEIF (S .LE. 7 .and. S .gt. 1.0) THEN
         IF    (N.EQ.1) THEN
            HELREC=HYDREC(T,1)*((1.626E-11_dp+2.073E-11_dp)/(1.499E-11_dp+1.847E-11_dp))
         ELSEIF(N.EQ.2) THEN
            HELREC=HYDREC(T,2)*((1.670E-11_dp+1.394E-11_dp)/(1.374E-11_dp+1.103E-11_dp))
         ELSEIF(N.EQ.3) THEN
            HELREC=HYDREC(T,3)* 1.0_dp                              ! free-free ist the same
         ELSEIF(N .EQ. 4) THEN
            HELREC=HYDREC(T,4)*(2.073_dp/1.847_dp)
         ELSEIF(N .EQ. 5) THEN
            HELREC=HYDREC(T,5)*(1.369E-11_dp/1.103E-11_dp)
         ELSE
            STOP 'HELREC: Bad parameter'
         ENDIF
      ELSE
         STOP 'HYDREC: Temperature out of range'
      ENDIF
      RETURN
!------------------------------------------------------------------------------
    end function HELREC
!------------------------------------------------------------------------------
!
    function SIVREC (T)
!
!   Computes the total recombination rate of SIV to SIII accounting
!   for both dielectronic and radiative processes.
!   The formulae and coefficients have been taken from
!   Abdel-Naby et al. A&A 537, A40 (2012)
!  
!   MOD 03-15 WEB
!   CALLED BY additional_siv_recombination (contained in MATRIXC, file1.f90)
!
      IMPLICIT none
!
      REAL(dp):: SIVREC,T
      REAL(dp):: A,B,T_0,T_1
!
      SIVREC=0._dp
!   dielectronic recombination rate
!   1/(T^(3/2) * sum c_i exp (-E_i/T)
      SIVREC=SIVREC+5.817e-7*exp(-3.628e+2/T)
      SIVREC=SIVREC+1.391e-6*exp(-1.058e+3/T)
      SIVREC=SIVREC+1.123e-5*exp(-7.160e+3/T)
      SIVREC=SIVREC+1.521e-4*exp(-3.260e+4/T)
      SIVREC=SIVREC+1.875e-3*exp(-1.235e+5/T)
      SIVREC=SIVREC+2.097e-2*exp(-2.070e+5/T)
      SIVREC=SIVREC*T**(-3._dp/2._dp)
!
!    radiative recombination rate
!    B= B + C *exp(-T2/T)
      B = 0.6896+ 0.084*exp(-6.752e+4/T)
!     RR = A* sqrt(T_0/T)*(((1+sqrt(T/T_0)**(1._dp-B))
!          *((1+sqrt(T/T_1)**(1._dp+B)))
      A=2.664e-10
      T_0=2.107e+1
      T_1=2.028e+7
      sivrec=sivrec+&
             A* sqrt(T_0/T)*(((1+sqrt(T/T_0))**(1._dp-B))*&
             ((1+sqrt(T/T_1))**(1._dp+B)))**(-1)
!------------------------------------------------------------------------------
    end function SIVREC
!------------------------------------------------------------------------------
!
end module M_Physical_functions
!==============================================================================
!                ***** END MODULE M_Physical_functions *****
!==============================================================================
!
!==============================================================================
!                ***** MODULE M_3d_min_healpix *****
!==============================================================================
Module M_3d_min_healpix
!
!--- Routine taken from the Healpix project ((c) Gorski et. al).
!--- It gives the components of an unit vector according to the Healpix Scheme
!--- Selection of the directions of the rays in the 3d case:
!     WM3D_MK_PIX2XY
!     WM3D_PIX2VEC_NEST
!
!--- M 2008/12.14 Web
!
  USE M_type_definitions  
  USE M_physical_constants , ONLY : PI
!  
 IMPLICIT NONE
! 
  INTEGER, PARAMETER:: LGT=4
  ! 2^13 : largest nside available
  INTEGER(i4b), PRIVATE, PARAMETER :: NS_MAX=8192 
  ! initialise array x2pix, y2pix and pix2x, pix2y used in several routines
  INTEGER(i4b), PRIVATE, SAVE, DIMENSION(128)    :: X2PIX=0, Y2PIX=0
  INTEGER(i4b), PRIVATE, SAVE, DIMENSION(0:1023) :: PIX2X=0, PIX2Y=0
!
CONTAINS
!------------------------------------------------------------------------------
  subroutine WM3D_PIX2VEC_NEST(NSIDE, IPIX, VECTOR, VERTEX)
!
!--- renders vector (x,y,z) coordinates of the nominal pixel center
!--- for the pixel number ipix (NESTED scheme)
!--- given the map resolution parameter nside
!--- also returns the (x,y,z) position of the 4 pixel vertices (=corners)
!--- in the order N,W,S,E
!
!--- CALLED BY WM3D_INIT_GEOMETRY
!--- M Web 08, Web 12-2014
!
    INTEGER(i4b), INTENT(IN)                    :: NSIDE, IPIX
    REAL   (dp ), INTENT(OUT), DIMENSION(1:)    :: VECTOR
    REAL   (dp ), OPTIONAL   , DIMENSION(1:,1:) :: VERTEX
!
    INTEGER(i4b) :: NPIX, NPFACE, &
                    IPF, IP_LOW, IP_TRUNC, IP_MED, IP_HI, &
                    JRT, JR, NR, JPT, JP, KSHIFT, NL4
    REAL   (dp ) :: Z, FN, FACT1, FACT2, STH, PHI
!
    INTEGER(i4b) :: IX, IY, FACE_NUM
    ! coordinate of the lowest corner of each face
    INTEGER(i4b), DIMENSION(1:12) :: JRLL=(/2, 2, 2, 2, 3, 3, 3, 3, 4, 4, 4, 4/) 
    ! in unit of nside
    INTEGER(i4b), DIMENSION(1:12) :: JPLL=(/1, 3, 5, 7, 0, 2, 4, 6, 1, 3, 5, 7/) 
    ! in unit of nside/2
!
    REAL   (dp ) :: PHI_NV, PHI_WV, PHI_SV, PHI_EV, PHI_UP, PHI_DN
    REAL   (dp ) :: Z_NV, Z_SV, STH_NV, STH_SV
    REAL   (dp ) :: HDELTA_PHI
    REAL   (dp ) :: HALFPI
    INTEGER(i4b) :: IPHI_MOD, IPHI_RAT
    LOGICAL(lgt) :: DO_VERTEX
!------------------------------------------------------------------------------
!    
    HALFPI=PI/2.0_dp
!
    IF (NSIDE<1 .OR. NSIDE>NS_MAX) STOP "nside out of range"
    npix = 12 * NSIDE**2
    IF (IPIX <0 .OR. IPIX>NPIX-1) STOP "ipix out of range"
!
    ! initiates the array for the pixel number -> (x,y) mapping
    IF (PIX2X(1023) <= 0)&
         CALL wm3d_mk_pix2xy()
!
    FN = REAL(NSIDE,KIND=dp)
    FACT1 = 1.0_dp/(3.0_dp*FN*FN)
    FACT2 = 2.0_dp/(3.0_dp*FN)
    NL4   = 4*NSIDE
!
    DO_VERTEX = .FALSE.
    IF (PRESENT(VERTEX)) THEN
       IF (SIZE(VERTEX,DIM=1) >= 3 .AND. SIZE(VERTEX,DIM=2) >= 4) THEN
          DO_VERTEX = .TRUE.
       ELSE
          STOP " pix2vec_ring : vertex array has wrong size "
       ENDIF
    ENDIF
!
    ! finds the face, and the number in the face
    NPFACE = NSIDE**2
    FACE_NUM = IPIX/NPFACE          ! face number in {0,11}
    IPF = MODULO(IPIX,NPFACE)       ! pixel number in the face {0,npface-1}
!
    ! finds the x,y on the face (starting from the lowest corner)
    ! from the pixel number
    IP_LOW = MODULO(IPF,1024)       ! content of the last 10 bits
    IP_TRUNC =   IPF/1024           ! truncation of the last 10 bits
    IP_MED = MODULO(IP_TRUNC,1024)  ! content of the next 10 bits
    IP_HI  =     IP_TRUNC/1024      ! content of the high weight 10 bits
!
    IX = 1024*PIX2X(IP_HI) + 32*PIX2X(IP_MED) + PIX2X(IP_LOW)
    IY = 1024*pix2y(IP_HI) + 32*pix2y(IP_MED) + pix2y(IP_LOW)
!
    ! transforms this in (horizontal, vertical) coordinates
    JRT = IX + IY                   ! 'vertical' in {0,2*(nside-1)}
    JPT = IX - IY                   ! 'horizontal' in {-nside+1,nside-1}
!
    ! computes the z coordinate on the sphere
    JR =  JRLL(FACE_NUM+1)*NSIDE - JRT - 1   ! ring number in {1,4*nside-1}
    NR = NSIDE                      ! equatorial region (the most frequent)
    Z  = (2*NSIDE-JR)*FACT2
    KSHIFT = MODULO(JR - NSIDE, 2)
    IF (DO_VERTEX) THEN
       Z_NV = (2*NSIDE-JR+1)*FACT2
       Z_SV = (2*NSIDE-JR-1)*FACT2
       IF (JR == NSIDE) THEN        ! northern transition
          Z_NV =  1.0_dp - (NSIDE-1)**2 * FACT1
       ELSEIF (JR == 3*NSIDE) THEN  ! southern transition
          Z_SV = -1.0_dp + (NSIDE-1)**2 * FACT1
       ENDIF
    ENDIF
    IF (JR < NSIDE) THEN            ! north pole region
       NR = JR
       Z = 1.0_dp - NR*NR*FACT1
       KSHIFT = 0
       IF (DO_VERTEX) THEN
          Z_NV = 1.0_dp - (NR-1)**2*FACT1
          Z_SV = 1.0_dp - (NR+1)**2*FACT1
       ENDIF
    ELSE IF (JR > 3*NSIDE) THEN     ! south pole region
       NR = NL4 - JR
       Z = - 1.0_dp + NR*NR*FACT1
       KSHIFT = 0
       IF (DO_VERTEX) THEN
          Z_NV = - 1.0_dp + (NR+1)**2*FACT1
          z_sv = - 1.0_dp + (NR-1)**2*FACT1
       ENDIF
    ENDIF
!    
    ! computes the phi coordinate on the sphere, in [0,2Pi]
    ! 'phi' number in the ring in {1,4*nr}
    JP = (JPLL(FACE_NUM+1)*NR + JPT + 1 + KSHIFT)/2
    IF (JP > NL4) JP = JP - NL4
    IF (jp < 1)   JP = JP + NL4
    PHI = (JP - (KSHIFT+1)*0.5_dp) * (HALFPI / NR)
    STH = SQRT((1.0_dp-Z)*(1.0_dp+Z))
    VECTOR(1) = STH * COS(phi)
    VECTOR(2) = STH * SIN(phi)
    VECTOR(3) = Z
!
    RETURN
!
  end subroutine WM3D_PIX2VEC_NEST
!-----------------------------------------------------------------------------
  subroutine wm3d_mk_pix2xy()
!
!--- constructs the array giving x and y in the face from pixel number
!--- for the nested (quad-cube like) ordering of pixels
!--- the bits corresponding to x and y are interleaved in the pixel number
!--- one breaks up the pixel number by even and odd bits
!
!--- CALED BY WM3D_PIX2VEC_NEST
!--- M WEB 08
!
    INTEGER(i4b) ::  KPIX, JPIX, IX, IY, IP, ID
!-----------------------------------------------------------------------
    ! print *, 'initiate pix2xy'
    DO KPIX=0,1023             ! pixel number
       JPIX = KPIX
       IX = 0
       IY = 0
       IP = 1                  ! bit position (in x and y)
       ! do while (jpix/=0)    ! go through all the bits
       DO
          IF (JPIX == 0) EXIT  ! go through all the bits
          ID = MODULO(JPIX,2)  ! bit value (in kpix), goes in ix
          JPIX = JPIX/2
          IX = ID*IP+IX

          ID = MODULO(JPIX,2)  ! bit value (in kpix), goes in iy
          jpix = JPIX/2
          IY = ID*IP+IY

          IP = 2*IP            ! next bit (in x and y)
       ENDDO
       PIX2X(KPIX) = IX        ! in 0,31
       PIX2Y(KPIX) = IY        ! in 0,31
    ENDDO
!    
   RETURN
!
  end  subroutine wm3d_mk_pix2xy
!------------------------------------------------------------------------------
!
end module M_3d_min_healpix
!==============================================================================
!                ***** END MODULE M_3d_min_healpix *****
!==============================================================================
!
!==============================================================================
!                          ***** MODULE M_A *****
!==============================================================================
!Module
!--- routines used for
!     xx
!     xx
!     xx
!--- x.05 - ADI

   !USE M_type_definitions
!
!CONTAINS
!------------------------------------------------------------------------------
  !subroutine
  !end subroutine
!------------------------------------------------------------------------------
  !function
  !end function
!------------------------------------------------------------------------------
!
!end module
!==============================================================================
!                        ***** END MODULE M_A *****
!==============================================================================

