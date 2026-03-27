!==============================================================================
! *****  Modgen contains MODULES defining global variables and settings  *****
!==============================================================================
! Modgen.f90
! ==========
!
!  Modules defining global variables and settings
!  6-7.98/11.99/2.02/2.05/12.14 - ADI
!  M 3.09 - WEB
!
! Content of Modgen (this file):
! ==============================
!  M_optional_settings              , USE ONLY: solely ONLY  -> maintain
!  M_global_settings                , USE ONLY: solely ONLY  -> maintain
!  M_physical_constants             , USE ONLY: solely ONLY  -> maintain
!
!==============================================================================
!                 ****** MODULE M_optional_settings *****
!==============================================================================
module M_optional_settings
!
   USE M_type_definitions
!
!   II      - Total number of ionisationstages
!   KEL     - Maximum number of calculated elements
!   KEL3D   - In original WM-Basic only dummy variable... only used in WM-Basic_3D
!   NUMEL   - number of elements for NITC=2
!   LEX     - Maximum number of levels per ionis.st. used as bf opacity sources
!   NFR     - Maximum number of ionisationstages as opacity sources
!   NFRT    - number of nue-points of the standard Kurucz grid
!   NUMNT   - Maximum number of depth points
!   NMERGE  - step size of add. Kurucz grid points
!   NADD    - number of add. nue-points due to shocks
!   NFS     - number of add. nue-points due to line blocking - Sampling M.
!   IXBL    - number of add. nue-points due to line blocking - Formal Int.
!             if OPTADD = T (see MAININP) IXBL should be equal NFS
!   IX1     - max of total number of nue-points
!   NSPEC   - Number of Kurucz-spectra in the OPTNEB_EXTSEPC case
!   NRF     - Number of radiation field calculations:
!             NRF=2: Continuum and cont.+all lines, old self shadow correction
!             NRF=3: Continuum and cont.+weak lines and cont.+all lines
!
     INTEGER(i4b), PARAMETER :: NSPEC  =  315
     INTEGER(i4b), PARAMETER :: NRF    =    2
     INTEGER(i4b), PARAMETER :: II     =  149
     INTEGER(i4b), PARAMETER :: KEL    =   30
     INTEGER(i4b), PARAMETER :: KEL3D  =   30
     INTEGER(i4b), PARAMETER :: NUMEL  =   12 !  9
     INTEGER(i4b), PARAMETER :: LEX    =   10 !  5,10,50
     INTEGER(i4b), PARAMETER :: NFR    =   48 ! 36 ! =NIONS*NUMEL  27!
     INTEGER(i4b), PARAMETER :: NFRT   =  342
     INTEGER(i4b), PARAMETER :: NUMNT  =   41
     INTEGER     , PARAMETER :: NMERGE =    1 !MOD
     INTEGER(i4b), PARAMETER :: NADD   = 2000
     INTEGER(i4b), PARAMETER :: NFS    = 8000 !2000
     INTEGER(i4b), PARAMETER :: IXBL   = 8000 !2000
     INTEGER     , PARAMETER :: IX1    = NFR*6*LEX+21+NFRT/NMERGE+2+NADD+NFS
!                                        4375 !MOD
     INTEGER(i4b), SAVE :: KILL1=20,KILL2=10,KILL3=4,KILL4=10
end module M_optional_settings
!==============================================================================
!                ****** END MODULE M_optional_settings *****
!==============================================================================
!
!!########### Use this module for 3-dimensional radiative transfer! ############
!!==============================================================================
!!                 ****** MODULE M_optional_settings *****
!!==============================================================================
!!
!module M_optional_settings
!USE M_type_definitions
!!   II      - Total number of ionisationstages
!!   KEL     - Maximum number of calculated elements
!!   KEL3D   - Set to value of original WM-Basic KEL to read in original data
!!             and subsequently remap and reduce it to ABUND with KEL=NUMEL elements
!!   NUMEL   - number of elements for NITC=2
!!   LEX     - Maximum number of levels per ionis.st. used as bf opacity sources
!!   NFR     - Maximum number of ionisationstages as opacity sources
!!   NFRT    - number of nue-points of the standard Kurucz grid
!!   NUMNT   - Maximum number of depth points
!!   NMERGE  - step size of add. Kurucz grid points
!!   NADD    - number of add. nue-points due to shocks
!!   NFS     - number of add. nue-points due to line blocking - Sampling M.
!!   IXBL    - number of add. nue-points due to line blocking - Formal Int.
!!             if OPTADD = T (see MAININP) IXBL should be equal NFS
!!   IX1     - max of total number of nue-points
!!   NSPEC   - Number of Kurucz-spectra in the OPTNEB_EXTSPEC case
!   NRF     - Number of radiation field calculations:
!             NRF=2: Continuum and cont.+all lines, old self shadow correction
!             NRF=3: Continuum and cont.+weak lines and cont.+all lines
!
!     INTEGER(i4b), PARAMETER :: NSPEC  =   3
!     INTEGER(i4b), PARAMETER :: NRF    =   2
!     INTEGER(i4b), PARAMETER :: II     =  46
!     INTEGER(i4b), PARAMETER :: KEL    =   9
!     INTEGER(i4b), PARAMETER :: KEL3D  =  30
!     INTEGER(i4b), PARAMETER :: NUMEL  =   9 ! 9
!     INTEGER(i4b), PARAMETER :: LEX    =   1 !  5
!     INTEGER(i4b), PARAMETER :: NFR    =  24 ! =NIONS*NUMEL  27!
!     INTEGER(i4b), PARAMETER :: NFRT   = 342
!     INTEGER(i4b), PARAMETER :: NUMNT  = 5*5*5
!     INTEGER     , PARAMETER :: NMERGE = 200 !MOD
!     INTEGER(i4b), PARAMETER :: NADD   =   1
!     INTEGER(i4b), PARAMETER :: NFS    =   1
!     INTEGER(i4b), PARAMETER :: IXBL   =   1 !2000
!     INTEGER     , PARAMETER :: IX1    =NFR*5*LEX+NFRT/NMERGE+2+NADD+NFS
!!                                       4375 !MOD
!     INTEGER(i4b), SAVE :: KILL1=20,KILL2=10,KILL3=0,KILL4=0
!   !  INTEGER(i4b), SAVE :: KILL1=20,KILL2=10,KILL3=4,KILL4=10   for O_3D_TIME = .FALSE.
!end module M_optional_settings
!!==============================================================================
!!                ****** END MODULE M_optional_settings *****
!!==============================================================================
!
!==============================================================================
!                 ****** MODULE M_global_settings *****
!==============================================================================
module M_global_settings
!
   USE M_type_definitions
   USE M_optional_settings
!
! Atomic data values and basic parameters
!   KIS     - Maximum of ionisationstages per element
!   LDR     - Maximum number of dielec. data sets per ionis.st.
!   LDRL    - Maximum number of levels per ionis.st. for LDR
!   IPHOTNM - Maximum number of photoiois. data sets per ionis.st.
     INTEGER(i4b), PARAMETER :: KIS    =   8
     INTEGER(i4b), PARAMETER :: LDR    = 394
     INTEGER(i4b), PARAMETER :: LDRL   =  50
     INTEGER(i4b), PARAMETER :: IPHOTNM= 106
!
!   NTZO    - Maximum number of FM intervals
!   INMAX   - Maximum number of simultaneously calculated ionis.st.
!   LPH     - Maximum number of levels per ionis.st.
!   NIONS   - Number of ionization stages per element treated in NLTE
!   LEVM    - Maximum number of levels per element
!   NCORE   - Number of core rays
     INTEGER(i4b), PARAMETER :: NTZO   =   7
     INTEGER(i4b), PARAMETER :: INMAX  =  99 ! 24*4+1+2  !75
     INTEGER(i4b), PARAMETER :: LPH    =  50
     INTEGER(i4b), SAVE      :: NIONS  =   3
     INTEGER(i4b), SAVE      :: LEVM   = 151 ! NIONS*LPH+1
     INTEGER(i4b), SAVE      :: NCORE  =   5
     INTEGER(i4b), PARAMETER :: NCCMF  =   5 ! HAS TO BE EQUAL TO NCORE FOR CMF
!
!   MODELN  - available number of Kurucz-models (EMTFLUX_NEW)
!   MODPL   - number of external  HOPF-PARAMETER-SETS
!   NFRE    - number of nue-points for exact photo-data
!   NQ      - max of main quantum number of approx. photo data
!   NFM     - number of grid-points for which the line-force is calc.
!   L3      - length of a record of the line-list
!   L3R     - Reduced L3 value used for temporary arrays
     INTEGER(i4b), PARAMETER :: MODELN=  60
     INTEGER(i4b), PARAMETER :: MODPL =  75
     INTEGER(i4b), PARAMETER :: NFRE  =  20
     INTEGER(i4b), PARAMETER :: NQ    =  10
     INTEGER(i4b), PARAMETER :: NFM   =NUMNT*3-2
     INTEGER(i4b), PARAMETER :: L3    =50000
     INTEGER(i4b), PARAMETER :: L3R   =L3
!
!  parameters for storage of photoionization values
!   LP1     - Number of levels for 2. ioiz. - a change of  LP. requires
!   LP2     - Number of levels for 3. ioiz. - a modification of the
!   LP3     - Number of levels for 4. ioiz. - structure of LMERK, IGYM
!                                           - NOTE ADI IS WATCHING YOU!*?
!   LPN     - Max. number of photoioniz. per level
!   IIR     - Integr. parts of photoioniz. rates - changed from  2 to  3
!   IIG     - I. of ph. r. for special exponents - changed from 11 to 13
!   KAUM    - Number of elements for which Auger-ioniz. rates are calc.
!   IAUM    - Number of iois.st. for which Auger-ioniz. rates are calc.
!   LMERK   - Max. number of photoioniz. per iois.st.
     INTEGER(i4b), PARAMETER :: LP1   =  34 !from 30-11.99
     INTEGER(i4b), PARAMETER :: LP2   =  22
     INTEGER(i4b), PARAMETER :: LP3   =   8
     INTEGER(i4b), PARAMETER :: LPN   =   4
     INTEGER(i4b), PARAMETER :: IIR   =   3
     INTEGER     , PARAMETER :: IIG   =  13 !MOD
     INTEGER(i4b), PARAMETER :: KAUM  =  11
     INTEGER(i4b), PARAMETER :: IAUM  =   4
     INTEGER(i4b), PARAMETER :: LMERK =LPH+LP1+LP2+LP3
!
!  only for formal.f90
!   NUMWAVE - Max. number of nue-points representing the spectral resol.
!             at the moment NUMWAVE < IX1 (reason: XHCONT,XJCONT)
!   NDOUT   - Number of depth points to print in stdout
     INTEGER(i4b), PARAMETER :: NUMWAVE= IX1 !<= IX1 from 601+IXBL 25.1.00
     INTEGER(i4b), PARAMETER :: NDOUT  =   9         !9
!
! other parameters
!   NCALC   - number of most important elements
!   NUMEL1  - number of elements for NITC=1
!   IONMAX  - max number of ionis.st. for dep. coef. c.
!   LEVNUM  - max of total number of levels
!   LINM    - max number of lines per ionis.st.
!   LINMC   - max number of col. trans. per ionis.st.
!   LINMAX  - max number of transitions
!   LININC  - LPH*INMAX
!   NUMKEL  - NUMNT*KEL
!   KELKIS  - KEL*KIS
     INTEGER(i4b), PARAMETER :: NCALC =  16
     INTEGER(i4b), PARAMETER :: NUMEL1=   2
     INTEGER(i4b), PARAMETER :: IONMAX=   4 !3 *DNS*
     INTEGER(i4b), PARAMETER :: LEVNUM=LPH*II
     INTEGER(i4b), PARAMETER :: LINM  =LPH*(LPH-1)/4
     INTEGER(i4b), PARAMETER :: LINMC =LINM*3/5
     INTEGER(i4b), PARAMETER :: LINMAX=LPH*(LPH-1)/2+LPH
     INTEGER(i4b), PARAMETER :: LININC=LPH*INMAX
     INTEGER(i4b), PARAMETER :: NUMKEL=NUMNT*KEL
     INTEGER(i4b), PARAMETER :: KELKIS=KEL*KIS
!
! the basic directory name for atomic data files is initialized
     INTEGER(i4b), PARAMETER :: LPNAM   =260
     CHARACTER (LEN=LPNAM)   :: ABDIRNAM='./'
!
end module M_global_settings
!==============================================================================
!                ****** END MODULE M_global_settings *****
!==============================================================================
!
!==============================================================================
!                 ****** MODULE M_physical_constants *****
!==============================================================================
!
module M_physical_constants
!
   USE M_type_definitions
   USE M_global_settings                  , ONLY: KEL
   USE M_all_physical_constants
!
  IMPLICIT NONE
!
   INTEGER(i4b), PRIVATE :: I
!
   INTEGER(i4b), PARAMETER, DIMENSION(KEL) :: ELEM_ARRAY_MASK=(/(I,I=1,KEL)/)
!
       ! All physical constants and constants which are a combination of those
       !======================================================================
       ! are to be found in M_all_physical_constants.
       !=============================================
!
       ! Arrays of constants to be found in M_all_physical_constants:
       !=============================================================
!         
         ! ATOMIC WEIGHTS OF THE ELEMENTS - 2001:
         !=======================================
      REAL(dp), DIMENSION (KEL), PARAMETER :: AWEIGHT(1:KEL) = &
                                          ALL_AWEIGHT(ELEM_ARRAY_MASK)
!
         ! SOLAR ABUNDANCES OF THE ELEMENTS - 1998.
         ! Standard solar composition from N. Grevesse and A. Sauval
         ! in: Space Science Reviews 85, 161, 1998:
         !=========================================
      REAL(dp), DIMENSION (KEL), PARAMETER :: SOLAR_ABUND_OLD(1:KEL) = &
                                          ALL_SOLAR_ABUND_OLD(ELEM_ARRAY_MASK)
!
         ! SOLAR ABUNDANCES OF THE ELEMENTS - 2009.
         ! Standard solar composition from M. Asplund, N. Grevesse, 
         ! A.J. Sauval, and  P. Scott
         ! "The Chemical Composition of the Sun"
         ! Annu. Rev. Astron. Astrophys. 2009, 47:481-522:
         !================================================
      REAL(dp), DIMENSION (KEL), PARAMETER :: SOLAR_ABUND(1:KEL) = &
                                          ALL_SOLAR_ABUND(ELEM_ARRAY_MASK)
!
         ! List of Elements considered in Atomic Number Order:
         !====================================================
      CHARACTER(2), DIMENSION(KEL), PARAMETER :: STRING1(1:KEL) = &
                                             ALL_STRING1(ELEM_ARRAY_MASK)
!
end module M_physical_constants
!==============================================================================
!               ****** END MODULE M_physical_constants *****
!==============================================================================
