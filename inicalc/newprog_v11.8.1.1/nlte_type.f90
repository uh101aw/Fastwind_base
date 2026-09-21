!version 1.2.1 (adapted for X-ray treatment)
!version 1.2 (for nlte.f90 vers. 9.0 and higher)
!new constants sqrt(pi) and pi/2.
!new precison qp: possible for recent intel compilers
!version 1.3: new constant eV
!version 1.4: new precision i1b, i2b
!version 1.5: additional fundamental constants
!version 1.6: polished

Module nlte_type
  Integer, Parameter :: i1b = selected_int_kind(2)
  Integer, Parameter :: i2b = selected_int_kind(4)
  Integer, Parameter :: i4b = selected_int_kind(9)
  Integer, Parameter :: sp = kind(1.0)
  Integer, Parameter :: dp = kind(1.0D0)
  Integer, Parameter :: qp = selected_real_kind(33,4931)
!  Integer, Parameter :: qp = 3 ! for NAG compiler

  Type :: rbb
    Integer (i1b) :: low
    Integer (i1b) :: lup
    Real (dp) :: wave
    Real (dp) :: gf
    Integer (i4b) :: index
    Integer (i4b) :: index_jbar
    Integer (i4b) :: index_comp
    Integer (i4b) :: quality !      0: OK;
!                                   1: differences with unpacked levels;
!                                   2: no corresponding transition in long file
    Integer (i2b) :: no_comp   !    number of corresponding unpacked components
    Integer (i4b) :: index1_ng !    index for ng-extrapolation
    Real (dp) :: sum_gf
  End Type

  Type :: cbb
    Integer (i1b) :: low
    Integer (i1b) :: lup
    Real (dp) :: omega
    Real (dp) :: slope
  End Type

  Type :: rdr
    Integer (i1b) :: low
    Integer (i2b) :: ncomp
    Integer (i4b) :: index
    Real (dp) :: ener
    Real (dp) :: flu
  End Type

  Type :: depacked_data
    Character (6) :: levlo
    Character (6) :: levup
    Integer (i4b) :: lineno
    Integer (i4b) :: nocomp
    Integer (i4b) :: index ! in subr. line_list, points to index in new line-list
!                            in subr. fgrid and later, points to index of lamcmf
    Real (dp) :: wave
    Real (dp) :: flu
  End Type

End Module

Module fund_const

  Use :: nlte_type
  Implicit None

! coconst

  Real (dp), Parameter :: pi = acos(-1.D0)
  Real (qp), Parameter :: piqp = acos(-1._qp)
!  Real(dp), Parameter :: piqp = acos(-1._dp) !for NAG compiler
  Real (dp), Parameter :: sqpi = sqrt(pi), pihalf = 0.5D0*pi

  Real (dp), Parameter :: hh = 6.6262D-27, akb = 1.38062D-16
  Real (dp), Parameter :: pm = 1.67265D-24, amh = 1.67352D-24, clight = 2.997925D10
  Real (dp), Parameter :: emass = 9.1095D-28
  Real (dp), Parameter :: gconst = 6.673D-8
  Real (dp), Parameter :: saha = 4.141584D-16, sigsb = 5.66956D-5
  Real (dp), Parameter :: sigmae = 0.3977246D0, sigmae_cross = 6.65D-25
  Real (dp), Parameter :: hc2 = 2.D0*hh*clight, hk = hh/akb
  Real (dp), Parameter :: hkl = hh*clight/akb
  Real (dp), Parameter :: e2mc2 = 2.817508719D-13
  Real (dp), Parameter :: ev = 1.602D-12 ! CONVERSION eV to erg
  Real (dp), Parameter :: log10e = log10(exp(1.D0)) ! used for new interpol. in intermepho
!
! line-parameters
! typical resonance line with f=1 at 1000 A.
  Real (sp), Parameter :: cross_class = 0.02654

  Real (dp), Parameter :: const_cross = cross_class*1.*1000.D-8
  Real (dp), Parameter :: c_vanreg = 0.5*2.055D-23
  Real (dp), Parameter :: c_forbid = 8.63D-6/(0.5*1.3707D-7)

  Real (dp), Parameter :: c_vanreg1 = 0.5*1.3707D-7 ! AVERAGE VALUE SINCE ALSO SUBORDINATE LINES
! CONVERSION TO ANGSTROM
  Real (dp), Parameter :: c_tau = 1.D-8*cross_class

! EINSTEIN COEFFICIENTS
  Real (dp), Parameter :: c_aul = 6.6702082D15, c_bul = 1.D-24/hc2
  Real (dp), Parameter :: c_blu = e2mc2/hh*4.D0*pi**2


  Real (dp), Parameter :: xmsun = 1.989D33, rsun = 6.96D10

End Module
