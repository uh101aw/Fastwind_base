!==============================================================================
! **** Datadefs contains MODULES defining the KIND numbers and Intefaces ****
!==============================================================================
!  Datadefs.f90
!
!  Modules defining the KIND numbers IMPLICIT types and Intefaces.
!  Note that IMPLICIT declarations are not importable.
!  3.98/7.98/8.03/2.05/12.14 - ADI
!  3.09 - WEB
!
! Content of Datadefs (this file):
! ================================
!  M_type_definitions               , ONLY: seldom
!  M_all_physical_constants           NO ONLY BUT ONLY USED VIA 
!                                     M_physical_constants
!                                     WHICH ACTS AS A WRAPPER
!
!==============================================================================
!                  *****  MODULE M_type_definitions  *****
!==============================================================================
module M_type_definitions
!
!   IMPLICIT REAL   (KIND=selected_real_kind(6 ,37 )) (a-h,o-z),& !Real4
!            INTEGER(KIND=selected_int_kind(9)      ) (i-n    )   !Intr4
!   IMPLICIT REAL   (KIND=selected_real_kind(15,307)) (a-h,o-z),& !Real8
!            INTEGER(KIND=selected_int_kind(9)      ) (i-n    )   !Intr4

! Indicator that the KIND= is not available for this compiler/host
    integer, parameter :: not_available = -1

! Real and Complex numbers - also used in the Interfaces below
!   Single precision
!   Real4                     = kind(0.0  )
    integer, parameter :: sp  = selected_real_kind(6 ,37 )
!   Double precision
!   Real8                     = kind(0.0d0)
    integer, parameter :: dp  = selected_real_kind(15,307)
!   Quadruple precision
    integer, parameter :: r16 = selected_real_kind( p=30 )
!   Dummy precision
    integer, parameter :: r17 = selected_real_kind( p=34 )

! Integers numbers         - also used in the Interfaces below
!   Default byte integer
    integer, parameter :: i0b   =              kind(0)
!   Single byte integer
    integer, parameter :: i1b   = selected_int_kind(2)
!   Two byte integer
    integer, parameter :: i2b   = selected_int_kind(4)
!   Four byte integer
    integer, parameter :: i4b   = selected_int_kind(9)
!   Default value
!   Eight byte integer
    integer, parameter :: i8b   = selected_int_kind(18)

! Logical values
!   Single byte logical
!   integer, parameter :: byte    = 1
!   Two byte logical
!   integer, parameter :: twobyte = 2
!   Four byte logical
!   integer, parameter :: word    = kind(.TRUE.)

! Character type
!   Normal single byte character (ASCII sequence)
    integer, parameter :: ascii   = kind('x')

! Generic names for Blas1/2/3 routines with overloaded versions
! =============================================================
!
   INTERFACE SSWAP
    SUBROUTINE SSWAP (n,x,incx,y,incy)
     REAL(kind(1.0e0)) :: x, y ! DIMENSION( : )
    END SUBROUTINE SSWAP
    SUBROUTINE DSWAP (n,x,incx,y,incy)
     REAL(kind(1.0d0)) :: x, y ! DIMENSION( : )
    END SUBROUTINE DSWAP
   END INTERFACE

   INTERFACE ISAMAX
    FUNCTION ISAMAX (n, x, incx)
     REAL(kind(1.0e0)) :: x ! DIMENSION( : )
    END FUNCTION ISAMAX
    FUNCTION IDAMAX (n, x, incx)
     REAL(kind(1.0d0)) :: x ! DIMENSION( : )
    END FUNCTION IDAMAX
   END INTERFACE

   INTERFACE SSCAL
    SUBROUTINE SSCAL (n, alpha, x, incx)
     REAL(kind(1.0e0)) :: alpha
     REAL(kind(1.0e0)) :: x ! DIMENSION( : )
    END SUBROUTINE SSCAL
    SUBROUTINE DSCAL (n, alpha, x, incx)
     REAL(kind(1.0d0)) :: alpha
     REAL(kind(1.0d0)) :: x ! DIMENSION( : )
    END SUBROUTINE DSCAL
   END INTERFACE

   INTERFACE SGEMV
    SUBROUTINE SGEMV (trans, m, n, alpha, a, lda, x, incx, beta, y, incy)
     REAL(kind(1.0e0)) :: alpha, beta
     REAL(kind(1.0e0)) :: a    ! DIMENSION(:,:)
     REAL(kind(1.0e0)) :: x, y ! DIMENSION( : )
     CHARACTER(LEN=*)  :: trans
    END SUBROUTINE SGEMV
    SUBROUTINE DGEMV (trans, m, n, alpha, a, lda, x, incx, beta, y, incy)
     REAL(kind(1.0d0)) :: alpha, beta
     REAL(kind(1.0d0)) :: a    ! DIMENSION(:,:)
     REAL(kind(1.0d0)) :: x, y ! DIMENSION( : )
     CHARACTER(LEN=*)  :: trans
    END SUBROUTINE DGEMV
   END INTERFACE

   INTERFACE STRSV
    SUBROUTINE STRSV (uplo, trans, diag, n, a, lda, x, incx)
     REAL(kind(1.0e0)) :: a ! DIMENSION(:,:)
     REAL(kind(1.0e0)) :: x ! DIMENSION( : )
     CHARACTER(LEN=*)  :: uplo, trans, diag
    END SUBROUTINE STRSV
    SUBROUTINE DTRSV (uplo, trans, diag, n, a, lda, x, incx)
     REAL(kind(1.0d0)) :: a ! DIMENSION(:,:)
     REAL(kind(1.0d0)) :: x ! DIMENSION( : )
     CHARACTER(LEN=*)  :: uplo, trans, diag
    END SUBROUTINE DTRSV
   END INTERFACE

   INTERFACE STRMV
    SUBROUTINE STRMV (uplo, trans, diag, n, a, lda, x, incx)
     REAL(kind(1.0e0)) :: a ! DIMENSION(:,:)
     REAL(kind(1.0e0)) :: x ! DIMENSION( : )
     CHARACTER(LEN=*)  :: uplo, trans, diag
    END SUBROUTINE STRMV
    SUBROUTINE DTRMV (uplo, trans, diag, n, a, lda, x, incx)
     REAL(kind(1.0d0)) :: a ! DIMENSION(:,:)
     REAL(kind(1.0d0)) :: x ! DIMENSION( : )
     CHARACTER(LEN=*)  :: uplo, trans, diag
    END SUBROUTINE DTRMV
   END INTERFACE

   INTERFACE SGEMM
    SUBROUTINE SGEMM (transa,transb,m,n,k,alpha,a,lda,b,ldb,beta,c,ldc)
     REAL(kind(1.0e0)) :: alpha, beta
     REAL(kind(1.0e0)) :: a, b, c ! DIMENSION(:,:)
     CHARACTER(LEN=*)  :: transa, transb
    END SUBROUTINE SGEMM
    SUBROUTINE DGEMM (transa,transb,m,n,k,alpha,a,lda,b,ldb,beta,c,ldc)
     REAL(kind(1.0d0)) :: alpha, beta
     REAL(kind(1.0d0)) :: a, b, c ! DIMENSION(:,:)
     CHARACTER(LEN=*)  :: transa, transb
    END SUBROUTINE DGEMM
   END INTERFACE

   INTERFACE STRSM
    SUBROUTINE STRSM (side,uplo,transa,diag,m,n,alpha,a,lda,b,ldb)
     REAL(kind(1.0e0)) :: alpha
     REAL(kind(1.0e0)) :: a, b ! DIMENSION(:,:)
     CHARACTER(LEN=*)  :: side, uplo, transa, diag
    END SUBROUTINE STRSM
    SUBROUTINE DTRSM (side,uplo,transa,diag,m,n,alpha,a,lda,b,ldb)
     REAL(kind(1.0d0)) :: alpha
     REAL(kind(1.0d0)) :: a, b ! DIMENSION(:,:)
     CHARACTER(LEN=*)  :: side, uplo, transa, diag
    END SUBROUTINE DTRSM
   END INTERFACE
!
! Interfaces for routines with optional arguments
! =============================================================================
!
   INTERFACE
   ! not necessary/allowed, since SORTH is now part of module
   ! M_Mathematical_functions - this INTERFACE is an outlook on future stuff
     SUBROUTINE SORTH_M ( N, RA, ASCND, RB, IC )
      IMPLICIT NONE
       !input
       LOGICAL       :: ASCND
       INTEGER(KIND=selected_int_kind(9)      )  :: N
       !inoutput
       REAL   (KIND=selected_real_kind(15,307)), DIMENSION (N)  :: RA
       REAL   (KIND=selected_real_kind(15,307)), DIMENSION (N), OPTIONAL :: RB
       INTEGER(KIND=selected_int_kind(9)      ), DIMENSION (N), OPTIONAL :: IC
     END SUBROUTINE
   END INTERFACE
!
end module M_type_definitions
!==============================================================================
!                  *****  MODULE M_type_definitions  *****
!==============================================================================
!
!==============================================================================
!                 ****** MODULE M_all_physical_constants *****
!==============================================================================
module M_all_physical_constants
!
  USE M_type_definitions
!
        ! Physical constants (cgs):
        !==========================
     REAL(dp)    , PARAMETER :: &
          GRAVC   =      6.67259E-08_dp,                  &
          CLIGHT  =   2.99792458E+10_dp,                  &
          XKBOL   =     1.380658E-16_dp,                  &
          HPLA    =    6.6260755E-27_dp,                  &
          HMASS   =       1.6733E-24_dp,                  &
          AUMASS  =   1.66053886E-24_dp,                  &
          PMASS   =    1.6726231E-24_dp,                  &
          EMASS   =    9.1093897E-28_dp,                  &
          ECHARGE =   4.80320680E-10_dp,                  &
          ELECVOLT=   1.60217733E-12_dp,                  &
          RYDCH   =     1.096788E+05_dp,                  &
          RYDINF  = 1.0973731569E+05_dp,                  &
          BOHRR   =   5.29177249E-09_dp
        ! BOHRR  a0 = (h/[2pi])^2/ (m e^2)
        ! RYDINF Ry = (h/[2pi])^2/(2m a0^2)
!
     REAL(dp)    , PARAMETER :: &
          SUNR    =       6.9599E+10_dp,                  &
          SUNM    =        1.989E+33_dp,                  &
          SUNL    =        3.826E+33_dp,                  &
          SUNMBOL =        4.830E+00_dp
!
     REAL(dp)    , PARAMETER :: &
          PARSEC  =   3.08567808E+18_dp,                  &
          ASTRUNIT=   1.49597870E+13_dp,                  &
          YEARINS =  3.155849984E+07_dp,                  &
          EARTHM  =        5.976E+27_dp
!
     REAL(dp)    , PARAMETER :: &
          PI      =3.141592653589793_dp
!
     ! Constants which are a combination of the constants above (cgs):
     !================================================================
        ! PI4      =          4 pi          = 12.5663706143592
        ! WPI      =          pi^1/2        = 1.77245385090552
        ! HPLAC    =           h c          = 1.986447461038579E+16
        ! HCK      =         h c / k        = 1.43876866033339
        ! HCKA     =      h c / k * 10^8    = 1.438768660333390E+08
        ! SIGBOL   =                        = 5.670508538505437E-05
        ! HCC2     =        h c^2/2         = 2.977609835163074E-06
        ! SAHAREC  = 1/2 h^3/(2pi me k)^3/2 = 2.070647865937159E-16
        ! BFACT    =      (2*h/c^2)*c^3     = 3.972894922077158E-16
        ! BFACTA3  = (2*h/c^2)*c^3 * 10^24  = 3.972894922077158E+08
        ! CROSSSC  =       pi e^2/(me c)    = 2.654009413735217E-02
        ! CROS4PHC =  pi e^2/(me c) 4pi/hc  = 1.678940246904783E+15
        ! SIGMATH  =                        = 6.652461528081342E-25
        ! SIGMATM  =        SIGMATH/mp      = 0.397726273664482
        ! RFINESC  =       hc/(2*pi e^2)    = 137.035990694620
        ! RYDNUE   =       1 RYD IN [HZ]    = 3.289841960502706E+15
        ! HIONEN   =h*c * Ry[cm^-1] =Ry[erg]= 2.179874121335895E-11
        ! HPLDRYDI =h[erg/s]/Ry[erg]=1/c/RyI= 3.039659691881352E-16
        ! EVDHPLDC =   eV / h / c [cm^-1]   = 8.065540928841530E+03
        ! AHPLCEV  =         h * c / eV [cm]= 1.239842446802427E-04
        ! AHPLCEVA =  10^8 * h * c / eV [A] = 1.239842446802427E+04
        ! XLRYDDK  =      Ry[erg] / k       = 1.578866106838837E+05
        ! HTCSSBF  = 7.91E-18 * Ry [cm^-1]  = 8.672380383492172E-13
        ! CFFOPAH  =    3.69234E-8 / c^3    = 1.370375497532096E-23
        ! COLBBCP  =                        = 8.629109466950885E-06
        ! COLBBC0  =                        = 5.465383915473385E-11
        ! COLBFCP  =                        = 1.552093168158680E+13
!        
     REAL(dp)    , PARAMETER :: &
          PI4     = 4._dp*PI,                                               &
          WPI     = PI**.5_dp,                                              &
          HPLAC   = HPLA * CLIGHT,                                          &
          HCK     = HPLA * CLIGHT / XKBOL,                                  &
          HCKA    = HPLA * CLIGHT / XKBOL * 1.E8_dp,                        &
          SIGBOL  = 2._dp*PI**5*XKBOL**4 / (15._dp*CLIGHT**2*HPLA**3),      &
          HCC2    = HPLA * CLIGHT**2 / 2._dp,                               &
          SAHAREC = .5_dp*(HPLA**3) / ((2._dp*PI*EMASS*XKBOL)**1.5),        &
          BFACT   = 2._dp*HPLA*CLIGHT,                                      &
          BFACTA3 = 2._dp*HPLA*CLIGHT * 1.E24_dp,                           &
          CROSSSC = PI*ECHARGE**2 / (EMASS*CLIGHT),                         &
          CROS4PHC= PI*ECHARGE**2 / (EMASS*CLIGHT) * 4._dp*PI/(HPLA*CLIGHT),&
          SIGMATH = 8._dp*PI*ECHARGE**4 / (3._dp*EMASS**2*CLIGHT**4),       &
          SIGMATM = 8._dp*PI*ECHARGE**4 / (3._dp*EMASS**2*CLIGHT**4)/PMASS, &
          RFINESC = HPLA * CLIGHT / (2._dp*PI*ECHARGE**2),                  &
          RYDNUE  = CLIGHT * RYDINF,                                        &
          HIONEN  = HPLA * CLIGHT * RYDINF,                                 &
          HPLDRYDI= 1._dp / (CLIGHT * RYDINF),                              &
          EVDHPLDC= ELECVOLT / (HPLA * CLIGHT),                             &
          AHPLCEV =           HPLA * CLIGHT / ELECVOLT,                     &
          AHPLCEVA= 1.E8_dp * HPLA * CLIGHT / ELECVOLT,                     &
          XLRYDDK = HPLA * CLIGHT * RYDINF / XKBOL,                         &
          HTCSSBF = 64._dp*PI**4*EMASS*ECHARGE**10/                         &
                    (3._dp*3._dp**.5*CLIGHT*HPLA**6)/                       &
                    (RYDINF**3*CLIGHT**3)   *   RYDCH ,                     &
          CFFOPAH = 4._dp*ECHARGE**6/(3._dp*CLIGHT*HPLA)*                   &
                    (2._dp*PI/(3._dp*XKBOL*EMASS**3))**.5_dp / CLIGHT**3 ,  &
          COLBBCP = (2._dp*PI/XKBOL)**.5*(HPLA/2._dp/PI)**2/EMASS**1.5,     &
          COLBBC0 = PI*BOHRR**2*(8._dp*XKBOL/PI/EMASS)**0.5,                &
          COLBFCP = 2._dp**1.5/(PI*EMASS*XKBOL)**.5 * 2._dp/(3._dp)**.5*    &
                    HPLA*CLIGHT/(2._dp*PI*ECHARGE**2) * HPLA*CLIGHT*RYDINF
!
        ! Constants which are a combination of the constants above,
        !==========================================================
        ! but are calculated elsewhere:
        !==============================
          ! EPSFACC  = COLBBC0*14.5*.5*RYDINF^2*HCK/BFACTA3/CROS4PHC
          !                                   = 1.029238437468084E-23
          ! EPSF3DTH = EPSFACC/SIGMATH * 10^24= 1.547154287361854E+25
       REAL(dp)    , SAVE      :: &
            EPSFACC  , &
            EPSF3DTH
!
         ! ATOMIC WEIGHTS OF THE ELEMENTS - 2001.
         ! Scaled to A_r(12C) = 12, where 12C is a neutral atom
         ! in its nuclear and electronic ground state - AUMASS.
         ! http://www.chem.qmw.ac.uk/iupac/AtWt/
         !
         ! List of Elements in Atomic Number Order:
         !=========================================
      REAL(dp), DIMENSION (92), PARAMETER :: ALL_AWEIGHT = &
          (/   1.00794_dp  , & !  1  H   Hydrogen
               4.002602_dp , & !  2  He  Helium
               6.941_dp    , & !  3  Li  Lithium
               9.012182_dp , & !  4  Be  Beryllium
              10.811_dp    , & !  5  B   Boron
              12.0107_dp   , & !  6  C   Carbon
              14.0067_dp   , & !  7  N   Nitrogen
              15.9994_dp   , & !  8  O   Oxygen
              18.9984032_dp, & !  9  F   Fluorine
              20.1797_dp   , & ! 10  Ne  Neon
              22.989770_dp , & ! 11  Na  Sodium
              24.3050_dp   , & ! 12  Mg  Magnesium
              26.981538_dp , & ! 13  Al  Aluminium
              28.0855_dp   , & ! 14  Si  Silicon
              30.973761_dp , & ! 15  P   Phosphorus
              32.065_dp    , & ! 16  S   Sulfur
              35.453_dp    , & ! 17  Cl  Chlorine
              39.948_dp    , & ! 18  Ar  Argon
              39.0983_dp   , & ! 19  K   Potassium
              40.078_dp    , & ! 20  Ca  Calcium
              44.955910_dp , & ! 21  Sc  Scandium
              47.867_dp    , & ! 22  Ti  Titanium
              50.9415_dp   , & ! 23  V   Vanadium
              51.9961_dp   , & ! 24  Cr  Chromium
              54.938049_dp , & ! 25  Mn  Manganese
              55.845_dp    , & ! 26  Fe  Iron
              58.933200_dp , & ! 27  Co  Cobalt
              58.6934_dp   , & ! 28  Ni  Nickel
              63.546_dp    , & ! 29  Cu  Copper
              65.409_dp    , & ! 30  Zn  Zinc
              69.723_dp    , & ! 31  Ga  Gallium
              72.64_dp     , & ! 32  Ge  Germanium
              74.92160_dp  , & ! 33  As  Arsenic
              78.96_dp     , & ! 34  Se  Selenium
              79.904_dp    , & ! 35  Br  Bromine
              83.798_dp    , & ! 36  Kr  Krypton
              85.4678_dp   , & ! 37  Rb  Rubidium
              87.62_dp     , & ! 38  Sr  Strontium
              88.90585_dp  , & ! 39  Y   Yttrium
              91.224_dp    , & ! 40  Zr  Zirconium
              92.90638_dp  , & ! 41  Nb  Niobium
              95.94_dp     , & ! 42  Mo  Molybdenum
              98.0_dp      , & ! 43  Tc  Technetium
             101.07_dp     , & ! 44  Ru  Ruthenium
             102.90550_dp  , & ! 45  Rh  Rhodium
             106.42_dp     , & ! 46  Pd  Palladium
             107.8682_dp   , & ! 47  Ag  Silver
             112.411_dp    , & ! 48  Cd  Cadmium
             114.818_dp    , & ! 49  In  Indium
             118.710_dp    , & ! 50  Sn  Tin
             121.760_dp    , & ! 51  Sb  Antimony
             127.60_dp     , & ! 52  Te  Tellurium
             126.90447_dp  , & ! 53  I   Iodine
             131.293_dp    , & ! 54  Xe  Xenon
             132.90545_dp  , & ! 55  Cs  Caesium
             137.327_dp    , & ! 56  Ba  Barium
             138.9055_dp   , & ! 57  La  Lanthanum
             140.116_dp    , & ! 58  Ce  Cerium
             140.90765_dp  , & ! 59  Pr  Praseodymium
             144.24_dp     , & ! 60  Nd  Neodymium
             145.0_dp      , & ! 61  Pm  Promethium
             150.36_dp     , & ! 62  Sm  Samarium
             151.964_dp    , & ! 63  Eu  Europium
             157.25_dp     , & ! 64  Gd  Gadolinium
             158.92534_dp  , & ! 65  Tb  Terbium
             162.500_dp    , & ! 66  Dy  Dysprosium
             164.93032_dp  , & ! 67  Ho  Holmium
             167.259_dp    , & ! 68  Er  Erbium
             168.93421_dp  , & ! 69  Tm  Thulium
             173.04_dp     , & ! 70  Yb  Ytterbium
             174.967_dp    , & ! 71  Lu  Lutetium
             178.49_dp     , & ! 72  Hf  Hafnium
             180.9479_dp   , & ! 73  Ta  Tantalum
             183.84_dp     , & ! 74  W   Tungsten
             186.207_dp    , & ! 75  Re  Rhenium
             190.23_dp     , & ! 76  Os  Osmium
             192.217_dp    , & ! 77  Ir  Iridium
             195.078_dp    , & ! 78  Pt  Platinum
             196.96655_dp  , & ! 79  Au  Gold
             200.59_dp     , & ! 80  Hg  Mercury
             204.3833_dp   , & ! 81  Tl  Thallium
             207.2_dp      , & ! 82  Pb  Lead
             208.98038_dp  , & ! 83  Bi  Bismuth
             209.0_dp      , & ! 84  Po  Polonium
             210.0_dp      , & ! 85  At  Astatine
             222.0_dp      , & ! 86  Rn  Radon
             223.0_dp      , & ! 87  Fr  Francium
             226.0_dp      , & ! 88  Ra  Radium
             227.0_dp      , & ! 89  Ac  Actinium
             232.0381_dp   , & ! 90  Th  Thorium
             231.03588_dp  , & ! 91  Pa  Protactinium
             238.02891_dp    & ! 92  U   Uranium
      /)
!
         ! SOLAR ABUNDANCES OF THE ELEMENTS - 1998.
         ! Standard solar composition from N. Grevesse and A. Sauval
         ! in: Space Science Reviews 85, 161, 1998.
         ! The values shown as comments on the right hand side
         ! are from Holweger, 1979.
         ! The values of the element abundances EL_sol are scaled to
         ! H_sol = 12 in the following way EL_sol=10^(EL_sol-H_sol).
         !
         ! List of Elements in Atomic Number Order:
         !=========================================
      REAL(dp), DIMENSION (30), PARAMETER :: ALL_SOLAR_ABUND_OLD = &
          (/   12.00_dp , & ! 12.00   1  H   Hydrogen
               11.00_dp , & ! 11.00   2  He  Helium
                1.10_dp , & !  3.28   3  Li  Lithium
                1.40_dp , & !  1.41   4  Be  Beryllium
                2.55_dp , & !  3.15   5  B   Boron
                8.52_dp , & !  8.67   6  C   Carbon
                7.92_dp , & !  7.99   7  N   Nitrogen
                8.83_dp , & !  8.92   8  O   Oxygen
                4.56_dp , & !  4.50   9  F   Fluorine
                8.08_dp , & !  7.73  10  Ne  Neon
                6.33_dp , & !  6.28  11  Na  Sodium
                7.58_dp , & !  7.53  12  Mg  Magnesium
                6.47_dp , & !  6.43  13  Al  Aluminium
                7.55_dp , & !  7.60  14  Si  Silicon
                5.45_dp , & !  5.35  15  P   Phosphorus
                7.33_dp , & !  7.20  16  S   Sulfur
                5.50_dp , & !  5.26  17  Cl  Chlorine
                6.40_dp , & !  6.83  18  Ar  Argon
                5.12_dp , & !  5.05  19  K   Potassium
                6.36_dp , & !  6.36  20  Ca  Calcium
                3.17_dp , & !  2.99  21  Sc  Scandium
                5.02_dp , & !  4.88  22  Ti  Titanium
                4.00_dp , & !  3.91  23  V   Vanadium
                5.67_dp , & !  5.61  24  Cr  Chromium
                5.39_dp , & !  5.47  25  Mn  Manganese
                7.50_dp , & !  7.60  26  Fe  Iron
                4.92_dp , & !  4.85  27  Co  Cobalt
                6.25_dp , & !  6.18  28  Ni  Nickel
                4.21_dp , & !  4.24  29  Cu  Copper
                4.60_dp  /) !  4.60  30  Zn  Zinc
!
         ! SOLAR ABUNDANCES OF THE ELEMENTS - 2009.
         ! Standard solar composition from M. Asplund, N. Grevesse, 
         ! A.J. Sauval, and  P. Scott
         ! "The Chemical Composition of the Sun"
         ! Annu. Rev. Astron. Astrophys. 2009, 47:481-522.
         ! To maintain compatibility between older WM-Basic test runs
         ! and newer ones the data from 1998 can still be used.
         ! The middle column contains the data from
         ! N. Grevesse and A. Sauval
         ! in: Space Science Reviews 85, 161, 1998,
         ! and the right-hand column the data from Hollweger 1979. 
         ! The values of the element abundances EL_sol are scaled to
         ! H_sol = 12 in the following way EL_sol=10^(EL_sol-H_sol).
         !
         ! List of Elements in Atomic Number Order:
         !=========================================
        REAL(dp), DIMENSION (30), PARAMETER :: ALL_SOLAR_ABUND = &
            (/   12.00_dp , & ! 12.00  ! 12.00   1  H   Hydrogen
                 10.93_dp , & ! 11.00  ! 11.00   2  He  Helium
                  1.05_dp , & !  1.10  !  3.28   3  Li  Lithium
                  1.38_dp , & !  1.40  !  1.41   4  Be  Beryllium
                  2.70_dp , & !  2.55  !  3.15   5  B   Boron
                  8.43_dp , & !  8.52  !  8.67   6  C   Carbon
                  7.83_dp , & !  7.92  !  7.99   7  N   Nitrogen
                  8.69_dp , & !  8.83  !  8.92   8  O   Oxygen
                  4.56_dp , & !  4.56  !  4.50   9  F   Fluorine
                  7.93_dp , & !  8.08  !  7.73  10  Ne  Neon
                  6.24_dp , & !  6.33  !  6.28  11  Na  Sodium
                  7.60_dp , & !  7.58  !  7.53  12  Mg  Magnesium
                  6.45_dp , & !  6.47  !  6.43  13  Al  Aluminium
                  7.51_dp , & !  7.55  !  7.60  14  Si  Silicon
                  5.41_dp , & !  5.45  !  5.35  15  P   Phosphorus
                  7.12_dp , & !  7.33  !  7.20  16  S   Sulfur
                  5.50_dp , & !  5.50  !  5.26  17  Cl  Chlorine
                  6.40_dp , & !  6.40  !  6.83  18  Ar  Argon
                  5.03_dp , & !  5.12  !  5.05  19  K   Potassium
                  6.34_dp , & !  6.36  !  6.36  20  Ca  Calcium
                  3.15_dp , & !  3.17  !  2.99  21  Sc  Scandium
                  4.95_dp , & !  5.02  !  4.88  22  Ti  Titanium
                  3.93_dp , & !  4.00  !  3.91  23  V   Vanadium
                  5.64_dp , & !  5.67  !  5.61  24  Cr  Chromium
                  5.43_dp , & !  5.39  !  5.47  25  Mn  Manganese
                  7.50_dp , & !  7.50  !  7.60  26  Fe  Iron
                  4.99_dp , & !  4.92  !  4.85  27  Co  Cobalt
                  6.22_dp , & !  6.25  !  6.18  28  Ni  Nickel
                  4.19_dp , & !  4.21  !  4.24  29  Cu  Copper
                  4.56_dp  /) !  4.60  !  4.60  30  Zn  Zinc
!
         ! List of Elements considered in Atomic Number Order:
         !====================================================
      CHARACTER(2), DIMENSION(30), PARAMETER :: ALL_STRING1 = &
              (/   ' H' , & !  1  H   Hydrogen
                   'He' , & !  2  He  Helium
                   'Li' , & !  3  Li  Lithium
                   'Be' , & !  4  Be  Beryllium
                   ' B' , & !  5  B   Boron
                   ' C' , & !  6  C   Carbon
                   ' N' , & !  7  N   Nitrogen
                   ' O' , & !  8  O   Oxygen
                   ' F' , & !  9  F   Fluorine
                   'Ne' , & ! 10  Ne  Neon
                   'Na' , & ! 11  Na  Sodium
                   'Mg' , & ! 12  Mg  Magnesium
                   'Al' , & ! 13  Al  Aluminium
                   'Si' , & ! 14  Si  Silicon
                   ' P' , & ! 15  P   Phosphorus
                   ' S' , & ! 16  S   Sulfur
                   'Cl' , & ! 17  Cl  Chlorine
                   'Ar' , & ! 18  Ar  Argon
                   ' K' , & ! 19  K   Potassium
                   'Ca' , & ! 20  Ca  Calcium
                   'Sc' , & ! 21  Sc  Scandium
                   'Ti' , & ! 22  Ti  Titanium
                   ' V' , & ! 23  V   Vanadium
                   'Cr' , & ! 24  Cr  Chromium
                   'Mn' , & ! 25  Mn  Manganese
                   'Fe' , & ! 26  Fe  Iron
                   'Co' , & ! 27  Co  Cobalt
                   'Ni' , & ! 28  Ni  Nickel
                   'Cu' , & ! 29  Cu  Copper
                   'Zn'  /) ! 30  Zn  Zinc
!
         ! Ionization stages considered:
         !==============================
      CHARACTER(5), DIMENSION(12) , PARAMETER :: STRING2(1:12) = &
        (/  '_I   ', '_II  ', '_III ', '_IV  ', '_V   ', '_VI  ', &
            '_VII ', '_VIII', '_IX  ', '_X   ', '_XI  ', '_XII '   /)
!
end module M_all_physical_constants
!==============================================================================
!               ****** END MODULE M_all_physical_constants *****
!==============================================================================


