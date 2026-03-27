Module formalsol_var

  Use :: nlte_type
  Use :: nlte_dim
  Implicit None


  Integer (i4b), Allocatable, Dimension (:, :) :: inverted
! ----------------------------------------------------------------------
! cprepr

  Real (dp), Dimension (id_depfi) :: vzrp ! changed Sept 2015 (previous dim
! too large)
  Real (dp), Dimension (2*id_depfi) :: xcmfp
  Real (dp), Dimension (id_cores+2+4*id_nocor) :: wp, wp1
! ----------------------------------------------------------------------
! cvoigt etc

  Real (dp), Allocatable, Dimension (:) :: gammal, gammac, weig, xnue

! ----------------------------------------------------------------------
! doppler,nemin

  Real (dp) :: vmax, xnemin, vturb, vturbmin, vturbmax
  Real (dp), Allocatable, Dimension (:) :: vturbv

! ----------------------------------------------------------------------
! exi

  Integer (i4b) :: nfcmf
  Real (dp), Dimension (id_nfesc) :: xcmfe
  Real (dp), Dimension (id_depfi, id_nfesc) :: sconte
! ----------------------------------------------------------------------
! exi1

  Real (dp), Dimension (id_nfesc) :: sce, desce
  Real (dp), Dimension (2*id_depfi, id_nfesc) :: sceray
! ----------------------------------------------------------------------
! file

  Character :: file*60, line*20

! ----------------------------------------------------------------------
! starkin

  Integer (i4b), Allocatable, Dimension (:) :: nstark, nlevl, nlevu, nws, nts, &
    nes

! ----------------------------------------------------------------------
! starklo

  Logical, Allocatable, Dimension (:) :: qhalf
! ----------------------------------------------------------------------
! starkre

  Real (dp), Allocatable, Dimension (:, :) :: dws
  Real (dp), Allocatable, Dimension (:, :) :: ts
  Real (dp), Allocatable, Dimension (:, :) :: es
  Real (dp), Allocatable, Dimension (:, :) :: ps

! additional clumping parameters for optically thick clumping
  Real (dp), Dimension (id_ndept) :: fic, tcl_fac_line, tcl_fac_cont

  Logical :: optthick


! indices, continuum fluxes from FLUXCONT

  Integer (i4b) :: klow, kup
  Real (dp) :: lamtrans, fcontlow, fcontup
  Real (dp), Parameter :: uvlimit = 3000.


End Module

!***********************************************************************

!main program

!***********************************************************************
Program forsol

  Use :: nlte_type
  Use :: nlte_dim
  Use :: fund_const, Only: clight
  Use :: ffr_error
  Use :: formalsol_var, Only: nstark, gammal, gammac, weig, xnue, vmax, vturb, &
    file, line, xnemin, vturbv, klow, kup, uvlimit
  Use :: formalsol_var, Only: inverted, nlevl, nlevu, qhalf, nws, nts, nes, &
    dws, ts, es, ps, optthick, fic, tcl_fac_line, tcl_fac_cont
  Use :: preformal_var, Only: at

  Implicit None

! formal solution for singlets/doublets/multiplets with continuum

! stark broadening and doppler broadening included
! including depth dependendance through t and n_e
!
! ---- from a program by puls, modified and generalized by
! ---- e. santolaya (1993) and j.p. (1996)

! ---- present philosophy: if stark-broadening, optiov = .true.
! ----                     else: dependency on delta
!
! ---- version 1.1 april 30th  1997 (j.puls) inclusion of preformal,
!   different input

! version 1.2 july 4th 1997
!   gf-values instead of f values, broadening option 2
!   for voigt profiles

! version 1.3 dec 16th 1997
!   new input-structure(formal_input)
!   new frequency grid

! version 1.4 feb 17th 1998
!   polished version, i4 now possible
!   exi files replaced by scratch files

! version 3.0.0 jan 22nd 1999
!   changed to f90, modifications by Chr. Wiethaus incorpo-
!   rated: calc. of xmax,  wavelength in air

! version 3.0.1 march 2nd 1999
!   f90 polished

! version 3.1.0 march 12th 1999
!   more f90 changes, obsolete features removed

! version 3.1.1 march 31st 1999
!   no large vturb, if stark-option

! version 4.0 april 16th 1999
!   different treatment of vturb (now separated from temp),
!   "all" values (however constant) allowed

! version 4.1 april 20th 1999
!   calculation of nemin changed. Now evaluated at
!   tauc = taucmin (0.03)

! version 4.2 may 20th 1999
!   inclusion of NSTARK = 3 option

! version 4.3 may  4th 2000
!   "cheat" to force overlap

! version 4.4 may  30th 2000
!   treatment of inverted levels changed: SL = 0. set
!   different xmax for Stark broadening

! version 4.5 oct 18th 2000 catalog can be 30 characters long

! version 5.0 july 2003 inclusion of approx. treatment of clumping
!   (to be consistent with version nlte_8.2 and higher
!   NOTE: until further evidence, the profile functions
!   are evaluated at the given (i.e., increased) n_e's

! version 5.1 april 2004 YHe read for use in preformal

! version 5.2 febr  2005 log(xnemin) set always to 10.5, &
!   perfil marginally modified to allow to be called
!   by calcxmax when xnemin has not been calculated

! version 5.3 oct   2005 xgrid changed for NB=2, to allow for sufficient
!   resolution of 2nd component

! version 5.3.1 sept  2006 new formulation to allaw for a better (though
!   not perfect) treatment of inverted levels

! version 5.3.2 june 2009 treatment of inverted levels even improved
!   (CONTIN: interpolation of SL and OPAL, EXPUNO)

!   STILL NOT TESTED: case of doublets, where one comp is in inversion
!   but note: since opa*sl > 0 always, this is also true for
!   Stot=Sum(opai*si)/Sum(opai): If Sum(opai) <0, then also Stot!

! version 6.0 july 2009 generalization to arbitrary number of components
!   (in case, change MAXNCOMP=10)
!   for NB > 2, general approach (obsfram/obsfram1) always

! version 6.0.1 sept 2009 some constants taken from nlte_type(v1.2)

! version 6.1   april 2012: approximate treatment of transitions between
!   high lying levels (mm and submm range), for ALMA
!   coding by NSTARK = 0 or 1 and LINENO=1 (NCOMP=1).
!   used departures are controlled by OPTDEP in CONTIN
!   OPTDEP=.FALSE -> DEP = 1, LTE
!   OPTDEP=.TRUE  -> DEP_LOW = DEP(19), DEP_UP = DEP(20)
!   to obtain similar ratios.

! version 7.0 april/may 2013: improvement of incoherent electron scattering:
!   modified approach allows to consider more than one line!
!   Inclusion of Stark-profiles from Lemke (1997, A&AS)
!   So far, only Balmer lines can be used.
!   To use them, uncomment lines after 'BROAD_BALMER'
!   in main program, and provide corresponding INPUT keyword
!   Files IX.DAT and STARK.BIN are now saved in model-directory,
!   to allow for simultaneous jobs


! version 7.1.0 november 2014: inclusion of optically thick clumping
!   (programmed by Jon Sundqvist), i.e. porosity and vorosity
!   needs nlte.f90 v10.3.0/nlte_approx.f90 v7.4.0 and higher
!   (CONT_FORMAL now includes OPA_EFF_RAT = Chi_eff/<Chi>,
!   additional file OUT_CLUMPING required).
!   incoherent e-scattering adapted for optically thick clumping


! version 8.0 july 2015: emergent spectrum in specific range, from
!   occupation numbers resulting from cmf-calculations
!   line-name = 'UV'
!   function perfil renamed to perfil1 (same name as in princesa)
!
!   missing:
!   vturb(r)
!   correct treatment of inversion
!   voigt broadening
!   unpack packed lines in expl. elements
!   e-scat???
!
!   NOTE: so far, we calculate individual lines using CONT_FORMAL,
!   whilst for the 'complete' spectra we use CONT_FORMAL_CMF.
!   The difference is that CONT_FORMAL uses pseudo-continua
!   for OPAC and SCONT (somewhat inconsistently calculated with J(cmf)),
!   whilst CONT_FORMAL_CMF uses pure continuum values. This gives,
!   for unbroadened spectra, slight differences due to a different
!   normalization and background, but incorporates the effects from the
!   additional lines present in the complete treatment. We checked that
!   using CONT_FORMAL_CMF for individual lines gives profiles that are
!   identical with those from the complete approach.

! version 8.1 july 2016: inclusion of depth-dependent micro-turbulence
!   for standard path (from v7.1 in path standard_xrays_v10.2.1)
!   NOTE: preformal (Stark etc. profiles) not affected,
!   since photospheric, i.e., calculated with Min(vturb).
!   IN THIS FIRST VERSION, ONLY ESCAT = 0 can be treated with
!   vturb(r). For ESCAT = 1, use constant vturb.
!   UPDATE of corresponding 'ranges' still required

!   THIS VERSION SHOULD BE CONSISTENT with
!   standard_v10.4.1/formalsol v7.2

! version 8.2  jan 11 2018: use SCON and OPAC from CONT_FORMAL_CMF
!   when calculating individual lines, and complete CMF model
!   has been calculated, since these values are somewhat different
!   from those provided by the approx. approach (saved in CONT_FORMAL)

! version 8.3  march 2018: depth-dependent vturb now also allowed
!   in "range"-branch, in the same spirit as for individual line complexes
!   NOTE: for individual lines, VTURBV is defined on the micro-grid,
!   whilst in the "range"-branch it is defined on the macro-grid.

! version 8.4  oct. 2020: compatible with gfortran

! version 8.5   nov 2021: triggered by problems in OV1371, a new
!   approach to find suitable cont. points (for opacities
!   and normalization) is used (when calculating individual
!   profiles). Only continuum point that behave smoothly
!   and are not contaminated by strong features in the
!   pseudo background are used.
!   Moreover, for lines below UVLIMIT (set to 3000 A),
!   the interpolation of cont. opacities and source-functions
!   is performed in the observer's frame (contrasted to the
!   CMF for optical/IR ranges), to avoid spurious results.
!   Finally, the consistency with the FLUXCONT output
!   is checked (differences arise because of differential
!   vs. integral methods). If an unreasonably large
!   inconsistency is found, the corresponding output
!   profile-file contains only zeros for cont. fluxes
!   and profiles.

!   NOTE: not checked for individual profiles based on
!   full cmf-treatment. In that case, however, there should
!   be no problems, because continuum is differently treated
!   anyway.

! version 8.5.1 nov 2022: inclusion of Lemke's Stark broadening for
!   potentially all hydrogen levels
!   default: Balmer-lines + important HeII lines: Butler
!   Paschen & Brackett lines until
!   upper level = 10: Lemke (correct)
!   other hydrogen and HeII lines: Griem

! version 8.5.2  feb 2023 -- affects only optically thick clumping:
!   bug found in contin (identified by Miguel).
!   Whenever a large number of components, opac was no
!   longer correctly corrected, and continuum opacity
!   became too small. Bug cured now, but there is still
!   a certain inconsistency between the pure continuum
!   case and the coupled transfer line+continuum. In
!   the latter case, the opac-correction would need to
!   be done frequency dependent, which is currently not
!   possible.

! version 8.5.3   small changes regarding gfortran

! version 8.5.4   polished

! !! WARNING !!!!
! for future calculations, remember that opac_nolines is mean, not
! effective opacity
! !! WARNING !!!!

! ----------------------------------------------------------------------

! .. parameters ..


  Integer (i4b), Parameter :: maxno = 500
  Integer (i4b), Parameter :: maxncomp = 310
  Integer (i4b), Parameter :: nd1 = id_ndept
  Integer (i4b), Parameter :: ndm = id_depfi, ltom = 2*ndm, nc = id_cores
  Integer (i4b), Parameter :: np = nc + 2 + 4*id_nocor
  Integer (i4b), Parameter :: nfobs = id_nfobs
! ..
! .. local scalars ..
  Real (dp) :: delta, deltax, realo, rmax, sr, srvmax, relxmax, tem, tgrad, &
    vsini, vturb1, xcblue, xcred, yhe, dum, dum1, lam_blu, lam_red, resol, &
    time_cpu, vtmi, vtma
  Integer (i4b) :: i, i1, i2, i3, iescat, integ, is, iturb, iu, nb, nch, nco, &
    nli, nlin, nsum, j, nd, iblu, ired

  Logical :: escat, optiov, exi, escat_uv, first

  Logical :: balmer_lemke, pb_lemke
! if true, Balmer and or Paschen/Brackett Stark profiles
! from Lemke will be used until nu=10, otherwise (nu>10) Griem
! BALMER_LEMKE only for specific tests (default false),
! PB=Paschen/Brackett default true, otherwise Griem-broadening
! for all upper levels

  Character :: dc*6, key*6, spec*11, ret*4, spec1*18
  Character :: broad_balmer*6, broad_pb*6, vturb_str*60
! ..
! .. local arrays ..
  Real (dp) :: p(np), profabs(nfobs), profem(nfobs), profile(nfobs), &
    profrot(nfobs), r(ndm), r1(nd1), rho(ndm), temp(ndm), v(ndm), x0(nfobs), &
    xne(nd1), xnefi(ndm), clf(nd1), z(ndm, np), zg(np, nfobs), zray(ltom), &
    opacon(ndm, 2), scont(ndm, 2), opacray(ltom, 2), sconray(ltom, 2), &
    aic(nfobs, 2), clf_test(nd1), fvel(nd1), hpor(nd1), fvol(nd1), v1(nd1), &
    dvdr1(nd1), temp1(nd1)

  Integer (i4b) :: index1(nd1), lmax(np), ltot(np), ncomp(maxno)

  Integer (i4b) :: lineno(maxncomp, maxno), nstarkl(maxncomp, maxno)

  Real (dp) :: xlam_blu(maxno), xlam_red(maxno), xresol(maxno)

  Character :: level(maxncomp, maxno)*6, leveu(maxncomp, maxno)*6, &
    lines(maxno)*20

  Integer (i4b), Allocatable, Dimension (:) :: nl, nu

  Real (dp), Allocatable, Dimension (:) :: gflu, gl, gu, vdop, xmax, xmaxdop, &
    xmaxdop_min, deltaarr
  Real (dp), Allocatable, Dimension (:, :) :: opal, sline
  Real (dp), Allocatable, Dimension (:, :) :: opalray, sray, xcmf
  Real (dp), Allocatable, Dimension (:) :: vturbray


! ..
! .. external subroutines ..
  External :: contin, elscat, ffracc, ffrcdc, ffrkey, ffrncd, ffrnum, formal, &
    model, preformal, staread
! ..
! .. intrinsic functions ..
! ..
  first = .True.

  Print *, ' INPUT CATALOGUE NAME FOR FORMAL CALC.'
  Read (*, Fmt='(A)') file

  Print *, ' Either INPUT of constant VTURB (KM/S)'
  Print *, ' or INPUT VTURBMIN (KM/S), VTURBMAX (KM/S or IN UNITS OF VINF)'
  Read (*, Fmt='(A)') vturb_str
! trick to read unknown number of variables (including description)
  Open (2, Status='SCRATCH', Form='FORMATTED')
  Write (2, *) vturb_str
  Rewind (2)
  Read (2, Fmt=*, End=100, Err=100) vtmi, vtma
  Close (2)
  Go To 110
! in case that only one input value for vturb (old approach)
100 vtma = vtmi

110 Print *, ' INPUT IESCAT (=0, NO; = 1 WITH ESCAT)'

  If (vtmi==0.D0 .And. vtma/=0.D0) Stop &
    ' VTURBMIN = 0 AND VTURBMAX NE 0! NOT ALLOWED'

  Read (*, Fmt=*) iescat
  If (iescat==1 .And. vtmi/=vtma) Stop &
    ' DEPTH DEPENDENT TURBULENCE ONLY FOR NO ESCAT'

  vturb = vtmi*1.D5 !               for photosphere
! VTURBMIN and VTURBMAX will obtain correct units in sub. MODEL

  broad_balmer = 'BUTLER'
! BROAD_BALMER='LEMKE'
! uncomment for tests with Lemke broadening functions
! PRINT *,' IN CASE, SPECIFY STARK BROADENING FOR BALMER LINES:'
! PRINT *,' BUTLER (default) or LEMKE'
! READ (*,FMT='(A)',END=5) BROAD_BALMER

! 5 IF(BROAD_BALMER.EQ.'') BROAD_BALMER ='BUTLER'

! BROAD_PB='GRIEM'
  broad_pb = 'LEMKE'
! uncomment for tests with Lemke broadening functions
! PRINT *,' IN CASE, SPECIFY STARK BROADENING FOR PASCHEN/BRACKETT LINES:'
! PRINT *,' GRIEM(DEFAULT FROM NU=11 ALWAYS) or LEMKE (DEFAULT UNTIL NU=10)'
! READ (*,FMT='(A)',END=6) BROAD_PB

! 6 IF(BROAD_PB.EQ.'') BROAD_PB ='LEMKE'

  If (broad_balmer/='BUTLER' .And. broad_balmer/='LEMKE ') &
    Stop ' WRONG BALMER BROADENING'
  If (broad_pb/='GRIEM' .And. broad_pb/='LEMKE ') Stop &
    ' WRONG PASCHEN/BRACKETT BROADENING'

  balmer_lemke = .False.
  If (broad_balmer=='LEMKE') balmer_lemke = .True.

  Print *, ' BALMER LINE BROADENING FOLLOWING ', broad_balmer

  pb_lemke = .False.
  If (broad_pb=='LEMKE') pb_lemke = .True.

  Print *, ' PASCHEN/BRACKETT LINE BROADENING (UNTIL NU=10) FOLLOWING ', &
    broad_pb
  Print *, ' PASCHEN/BRACKETT LINE BROADENING (FROM  NU=11) FOLLOWING GRIEM'

  dc = ':T'

  Open (1, File='FORMAL_INPUT', Status='OLD', Form='FORMATTED')
  iu = 1
  Call ffracc(iu, ret)
  Select Case (ret)
  Case ('RET1')
    Go To 150
  Case ('RET0')
    Continue
  Case Default
    Stop ' WRONG RETURN IN FORMALSOL'
  End Select

  Call ffrcdc(dc)

! ---  allowing for all characters in string 'file'

  Call ffrncd(ret)
  If (on_error(ret)) Go To 160

  Call ffrnum(realo, integ, ret)
  If (on_error(ret)) Go To 160
  vsini = realo

  If (iescat==0) Then
    escat = .False.
  Else If (iescat==1) Then
    escat = .True.
  Else
    Stop ' ERROR IN ESCAT OPTION'
  End If

  nlin = 0
  Print *, ' LINES TO BE TREATED'

comploop: Do

    nlin = nlin + 1
    If (nlin>maxno) Stop 'CHANGE MAXNO'

    Call ffrkey(lines(nlin), nch, ret)
    Select Case (ret)
    Case ('RET1')
      Go To 120
    Case ('RET2')
      Go To 160
    Case ('RET0')
      Continue
    Case Default
      Stop ' WRONG RETURN IN FORMALSOL'
    End Select

    Print *, lines(nlin)

    If (trim(lines(nlin))/='UV') Then
!     standard treatment
      Call ffrnum(realo, nco, ret)
      If (on_error(ret)) Go To 160

      ncomp(nlin) = nco

      Do i = 1, nco
        Call ffrkey(key, nch, ret)
        If (on_error(ret)) Go To 160
        level(i, nlin) = key

        Call ffrkey(key, nch, ret)
        If (on_error(ret)) Go To 160
        leveu(i, nlin) = key

        Call ffrnum(realo, integ, ret)
        If (on_error(ret)) Go To 160
        lineno(i, nlin) = integ

        Call ffrnum(realo, integ, ret)
        If (on_error(ret)) Go To 160
        nstarkl(i, nlin) = integ

      End Do

    Else
!     spectral range
      Call ffrnum(realo, integ, ret)
      If (on_error(ret)) Go To 160
      xlam_blu(nlin) = realo

      Call ffrnum(realo, integ, ret)
      If (on_error(ret)) Go To 160
      xlam_red(nlin) = realo

      Call ffrnum(realo, integ, ret)
      If (on_error(ret)) Go To 160
      xresol(nlin) = realo

      delta = xlam_red(nlin) - xlam_blu(nlin)
      If (delta<=0.D0) Stop ' LAMBDA_BLUE > LAMBDA_RED, MODIFY INPUT'
!     HERE WAS A BUG (RESOL INSTEAD OF XRESOL(NLIN))
      If (delta<=xresol(nlin)) Stop ' RESOL TOO LARGE, MODIFY INPUT'

    End If

  End Do comploop

! end of file formal-input

120 Continue
  Write (*, Fmt='(A)') ' END OF FORMAL-INPUT, NO MORE LINES '
  Close (1)

  nlin = nlin - 1

  Open (1, File=trim(file)//'/MODEL', Status='OLD', Form='UNFORMATTED')
  Rewind 1
  Read (1) dum, dum, dum, yhe
  Close (1)

all_lines: Do nli = 1, nlin

    If (trim(lines(nli))/='UV') Then
!     standard treatment
      nco = ncomp(nli)
      If (nco>maxncomp) Stop ' TOO MANY COMPONENTS (INCREASE MAXNCOMP)'

      Allocate (inverted(id_ndept,nco))
      Allocate (gammal(nco), gammac(nco), weig(nco), xnue(nco))
      Allocate (nstark(nco), nlevl(nco), nlevu(nco), nws(nco), nts(nco), &
        nes(nco))
      Allocate (nl(nco), nu(nco))
      Allocate (gflu(nco), gl(nco), gu(nco), vdop(nco))
      Allocate (qhalf(nco))
      Allocate (dws(id_maxww,nco))
      Allocate (ts(id_maxtt,nco))
      Allocate (es(id_maxne,nco))
      Allocate (ps(id_maxps,nco))

      Allocate (opal(ndm,nco), sline(ndm,nco))
      Allocate (xmax(nco), xmaxdop(nco), xmaxdop_min(nco), deltaarr(nco))

      Allocate (opalray(ltom,nco), sray(ltom,nco), xcmf(ltom,nco), &
        vturbray(ltom))

      line = lines(nli)
      Do i = 1, nco
        nstark(i) = nstarkl(i, nli)
      End Do

!     ---- begin of calculation for different lines
!     ---- for PREFORMAL, WE USE ONLY VTURB = VTURBMIN

!     ARGUMENT .FALSE. CORRESPONDS TO SRANGE IN PREFORMAL
      Call preformal(file, nco, level(1,nli), leveu(1,nli), lineno(1,nli), &
        nstark, vturb, yhe, balmer_lemke, pb_lemke, .False.)

      Open (1, File=trim(file)//'/IX.DAT', Status='UNKNOWN')
      Rewind 1
      Read (1, Fmt=*) nb
      If (nb/=nco) Stop ' SOMETHING WRONG WITH NO OF COMPONENTS'
      Do i = 1, nb
        Read (1, Fmt=*) xnue(i)
        Read (1, Fmt=*) gl(i), gu(i)
        Read (1, Fmt=*) gflu(i)
        Read (1, Fmt=*) nl(i), nu(i)
        Read (1, Fmt=*) weig(i)
        Read (1, Fmt=*) nstark(i)
        If (nstark(i)==2 .Or. nstark(i)==3) Read (1, Fmt=*) gammal(i), &
          gammac(i)
      End Do
      Close (1)

!     approximate treatment so far only for ESCAT = .FALSE.
      If (at .And. escat) Stop &
        ' Approximate Treatment (AT) AND ESCAT NOT POSSIBLE YET'

!     Read in clumping params from OUTPUT_CLUMPING
      Open (1, File=trim(file)//'/CLUMPING_OUTPUT', Status='OLD')
      Read (1, *) !                 header
      Do i = 1, nd1
        Read (1, Fmt=*) j, dum1, dum1, dum1, clf_test(i), fic(i), fvel(i), &
          hpor(i), fvol(i), tcl_fac_line(i), tcl_fac_cont(i)
      End Do
      Close (1)

      optthick = .True.
      If (maxval(tcl_fac_line)==0. .And. maxval(tcl_fac_cont)==0.) &
        optthick = .False.
      If (optthick) Then
        Print *
        Print *, 'MODEL WITH OPTICALLY THICK CLUMPING'
      End If

!     ------------------------------------
!     output-filename and open output; still named with vturb=vturbmin
!     but profiles with variable vturb denoted by 'VTV' instead of 'VT'

      is = index(line, ' ')
      vturb1 = vturb*1.D-5
      If (vturb1==0.) Then
        If (escat) Then
          spec = 'ESC'
          Open (2, File=trim(file)//'/OUT.'//line(1:is-1)//'_'//spec, &
            Status='UNKNOWN')
        Else
          spec = ' '
          Open (2, File=trim(file)//'/OUT.'//line(1:is-1), Status='UNKNOWN')
        End If
      Else
        iturb = int(vturb1)
        If (iturb>=1000) Stop ' VTURB > 1000 KM/S'
        i1 = iturb/100
        i2 = (iturb-i1*100)/10
        i3 = iturb - i1*100 - i2*10

        If (escat) Then
          spec = 'ESC_VT' // char(48+i1) // char(48+i2) // char(48+i3)
        Else
          If (vtmi==vtma) Then
            spec = 'VT' // char(48+i1) // char(48+i2) // char(48+i3)
          Else
            spec = 'VTV' // char(48+i1) // char(48+i2) // char(48+i3)
          End If
        End If

        Open (2, File=trim(file)//'/OUT.'//line(1:is-1)//'_'//spec, &
          Status='UNKNOWN')
      End If

!     ------------------------------------

!     initialize xnemin (to be checked in perfil)

      xnemin = 0.

      nsum = 0
      Do i = 1, nb
        If (nstark(i)==1) nsum = nsum + 1

!       ------- xnue in cm-1

        xnue(i) = 1.D8/xnue(i)

!       ------- calculation of mass dependent part of VDOP (rest is calculated
!       in subroutine model)

        vdop(i) = sqrt(1./weig(i))
      End Do

!     ---- stark profiles are read

      If (nsum>0) Then
        Call staread(file, nb, nsum, nl, nu)
      End If

      Do i = 1, nb
        If (nstark(i)==2 .Or. nstark(i)==3) nsum = nsum + 10
      End Do

!     ---- hence: first digit in nsum: number of voigt components
!     ----          2nd digit in nsum: number of stark components

      If (escat) Then
        If (1.D8/xnue(1)<uvlimit) Stop &
          ' ELECTRON SCATTERING ONLY POSSIBLE FOR LAM > UVLIMIT'
        Call elscat(xnue, nl, nu, gflu, gl, gu, nb)
      End If

      Write (*, Fmt=170) line
      Print *
      Print *, 'VSINI: ', vsini, ' KM/S, VTURBMIN: ', vturb*1.D-5, ' KM/S'

      Call model(r1, r, v1, v, dvdr1, rho, rmax, ndm, nd, vmax, vdop, index1, &
        srvmax, sr, nb, xne, clf, clf_test, vtmi, vtma)

!     from here on, VDOP is correct (at TEFF, corrected for VTURB), in VMAX
!     VDOP calculated with VTURB=VTURBMIN;

      delta = 0.D0
!     changed from v 6.0, since x-grid w.r.t. xnue(1)
      If (nb>=2) Then
        delta = (xnue(1)-xnue(nb))/xnue(1)*clight/vmax
        If (delta<=0.D0) Stop ' ERROR IN DELTA'
        Do i = 2, nb - 1
          deltaarr(i) = (xnue(1)-xnue(i))/xnue(1)*clight/vmax
          If (deltaarr(i)<=0.D0) Stop ' ERROR IN DELTA'
        End Do
        deltaarr(nb) = delta
      End If

      Print *
      Print *, ' DELTA = ', delta

!     ---- determination of continuum quantities

      Call contin(xnue, nl, nu, gflu, srvmax, index1, r, opal, opacon, sline, &
        scont, gl, gu, nb, sr, xcred, xcblue, tem, tgrad, xmax, xmaxdop, &
        xmaxdop_min, temp, xne, xnefi, clf, escat)

!     no profile will be calculated:
      If (klow==-1 .And. kup==-1) Go To 140


!     ---- determination of the implicit overlapping option


130   deltax = delta - xmax(1) - xmax(nb)
      If (deltax>0.D0) Then
        Write (*, Fmt=180) deltax, delta
        optiov = .False.

      Else
        If (delta==0.D0) Then
          Write (*, Fmt=200) xmax(1)
          optiov = .True.
        Else
          Write (*, Fmt=190) deltax, delta
          optiov = .True.
        End If
      End If

!     -----more than 2 pure Doppler profiles calculated with general approach

      If (nsum==0 .And. nb>2) optiov = .True.

!     -----in present philosophy, stark/voigt broadening requires
!     -----optiov = .true.

      If (nsum>0 .And. .Not. optiov) Then

!       CHEAT TO OBTAIN ALWAYS OVERLAP (=> DELTAX=-0.01): RENORMALIZE XMAX

        relxmax = xmax(nb)/xmax(1)
        xmax(1) = (delta+0.01)/(1.+relxmax)
        xmax(nb) = xmax(1)*relxmax
        Go To 130
      End If

!     -----formal integral

140   Call formal(nd, np, nc, nb, nfobs, rmax, delta, escat, r1, r, v, opal, &
        sline, p, z, lmax, x0, vmax, opalray, sray, zray, vturbray, ltot, &
        xcmf, profile, profabs, profem, zg, aic, xmax, xmaxdop, xmaxdop_min, &
        optiov, nsum, opacon, opacray, scont, sconray, tem, xcred, xcblue, &
        tgrad, xnue(1), vdop, temp, xnefi, profrot, vsini, deltaarr)

      Close (2)

      Deallocate (inverted)
      Deallocate (gammal, gammac, weig, xnue)
      Deallocate (nstark, nlevl, nlevu, nws, nts, nes, nl, nu)
      Deallocate (gflu, gl, gu, vdop)
      Deallocate (qhalf)
      Deallocate (dws, ts, es, ps)

      Deallocate (opal, sline)
      Deallocate (vturbv)
      Deallocate (xmax, xmaxdop, xmaxdop_min, deltaarr)

      Deallocate (opalray, sray, xcmf, vturbray)

    Else
!     spectral range
!     so far, no incoherent e-scat
      lam_blu = xlam_blu(nli)
      lam_red = xlam_red(nli)
      resol = xresol(nli)

      escat_uv = .False.
!     IF(VTURB.LT.5.D0) STOP ' UV AND VTURB < 5 KM/S'
      If (vturb<1.D0) Stop ' UV AND VTURB < 1 KM/S'

      Inquire (File=trim(file)//'/CONT_FORMAL_CMF', Exist=exi)
      If (.Not. exi) Stop &
        ' UV range not possible, since CONT_FORMAL_CMF not existent!'

      nb = 1
      Allocate (vdop(nb))
      Allocate (opalray(ltom,1), sray(ltom,1))

!     Read in clumping params from OUTPUT_CLUMPING
      Open (1, File=trim(file)//'/CLUMPING_OUTPUT', Status='OLD')
      Read (1, *) !                 header
      Do i = 1, nd1
        Read (1, Fmt=*) j, dum1, dum1, dum1, clf_test(i), fic(i), fvel(i), &
          hpor(i), fvol(i), tcl_fac_line(i), tcl_fac_cont(i)
      End Do
      Close (1)

      optthick = .True.
      If (maxval(tcl_fac_line)==0. .And. maxval(tcl_fac_cont)==0.) &
        optthick = .False.
      If (optthick) Then
        Print *
        Print *, 'MODEL WITH OPTICALLY THICK CLUMPING'
      End If


!     ------------------------------------
!     output-filename and open output

      line = lines(nli)
      is = index(line, ' ')
      vturb1 = vturb*1.D-5
      If (resol<10.) Then
        Write (spec1, Fmt='(I5.5,A1,I5.5,A2,F5.3)') int(lam_blu), '_', &
          int(lam_red), '__', resol
      Else If (resol<100.) Then
        Write (spec1, Fmt='(I5.5,A1,I5.5,A1,F6.3)') int(lam_blu), '_', &
          int(lam_red), '_', resol
      Else
        Stop ' RESOL TOO LARGE, CHANGE FORMAT'
      End If

      If (vturb1==0.) Then
        If (escat_uv) Then
          spec = 'ESC'
          Open (2, File=trim(file)//'/OUT.'//line(1:is-1)//'_'//spec1//'_'// &
            spec, Status='UNKNOWN')
        Else
          spec = ' '
          Open (2, File=trim(file)//'/OUT.'//line(1:is-1)//'_'//spec1, &
            Status='UNKNOWN')
        End If
      Else
        iturb = int(vturb1)
        If (iturb>=1000) Stop ' VTURB > 1000 KM/S'
        i1 = iturb/100
        i2 = (iturb-i1*100)/10
        i3 = iturb - i1*100 - i2*10

        If (escat_uv) Then
          spec = 'ESC_VT' // char(48+i1) // char(48+i2) // char(48+i3)
        Else
          If (vtmi==vtma) Then
            spec = 'VT' // char(48+i1) // char(48+i2) // char(48+i3)
          Else
            spec = 'VTV' // char(48+i1) // char(48+i2) // char(48+i3)
          End If
        End If

        Open (2, File=trim(file)//'/OUT.'//line(1:is-1)//'_'//spec1//'_'//spec &
          , Status='UNKNOWN')
      End If

      Write (*, Fmt=170) line
      Print *, ' SPECTRAL RANGE:', lam_blu, ' ', lam_red
      Print *, ' RESOL = ', resol, ' ANGSTROM'
      Print *
      Print *, 'VSINI: ', vsini, ' KM/S, VTURBMIN: ', vturb*1.D-5, ' KM/S'

!     ----mass dependent part of VDOP (rest is calculated in subroutine model)
      vdop(nb) = 0.
      Call model(r1, r, v1, v, dvdr1, rho, rmax, ndm, nd, vmax, vdop, index1, &
        srvmax, sr, nb, xne, clf, clf_test, vtmi, vtma)
!     on output, vdop(nb)=vdop(1)=vturbmin/vinf

      If (first) Then
!       call prince only once! (2nd call will not work)
        Call prince
        first = .False.
      End If

      Call calc_lines(nd1, xne, temp1, clf, r1, v1, dvdr1, sr, srvmax, &
        lam_blu, lam_red, iblu, ired, vturb, yhe, balmer_lemke, pb_lemke, &
        maxncomp)


      Call formal_range(nd, np, nc, vdop(nb), rmax, r, p, z, v, ltot, lmax, &
        r1, xne, temp1, index1, iblu, ired, lam_blu, lam_red, resol, sr, &
        vmax, opacon, scont, zray, opacray(1,1), sconray(1,1), opalray(1,1), &
        sray(1,1))

      Close (2)

      Deallocate (vdop, opalray, sray, vturbv)
      Call dealloc
!     PRINT*,' ALL ARRAYS DEALLOCATED'
      Print *

    End If

  End Do all_lines

  Call cpu_time(time_cpu)
  Print *
  Print *, ' CPU time: ', time_cpu

  Stop

! error in input unit exit

150 Continue
  Write (*, Fmt='(A)') ' ERROR IN INPUT UNIT NUMBER '
  Stop

! error conditions
160 Select Case (ret)
  Case ('RET1') !                   end of file "no esperado" exit
    Write (*, Fmt='(A)') ' EOF NO ESPERADO. ALGO FALLA! '
    Stop
  Case ('RET2') !                   error in ffr subroutines exit
    Write (*, Fmt='(A)') ' ERROR IN FFR SUBROUTINES '
    Stop
  Case Default
    Stop ' WRONG ERROR CONDITION IN FORMALSOL'
  End Select

170 Format ('LINE: ', A)
180 Format (' DELTAX =', F7.2, '  DELTA =', F7.2, ' WELL SEPARATED')
190 Format (' DELTAX =', F7.2, '  DELTA =', F7.2, ' OVERLAP !!')
200 Format (' SINGLE WITH XMAX =', F7.4)
End Program

!***********************************************************************

!subroutines: simple ones

!***********************************************************************

Subroutine contin(xnue, nl, nu, gflu, srvmax, index1, r, opal, opacon, sline, &
  scont, gl, gu, nb, sr, xcred, xcblue, tem, tgrad, xmax, xmaxdop, &
  xmaxdop_min, temp, xne, xnefi, clf, escat)


! includes clumping
! NOTE: continuum opacities already corrected
! line quantities need to be corrected (assumption: resonance zone
! large compared to clumps)

! sline and opal are now calculated from interpolated nl, nu values

! approximate treatment if AT = .TRUE.
! used departures for AT-case are controlled by OPTDEP
! OPTDEP=.FALSE -> DEP = 1, LTE
! OPTDEP=.TRUE  -> DEP_LOW = DEP(19), DEP_UP = DEP(20)
! to obtain similar ratios.

  Use :: nlte_type
  Use :: nlte_dim
  Use :: fund_const, Only: clight
  Use :: formalsol_var, Only: gammal, gammac, xcmfe, sconte, nfcmf, file, &
    inverted, fic, tcl_fac_line, optthick, klow, kup, fcontlow, fcontup, &
    lamtrans, uvlimit

  Use :: preformal_var, Only: at, nlevh21
  Implicit None

! ---- reads and calculates continuum quantities

! .. parameters ..

  Logical, Parameter :: optdep = .False.

  Integer (i4b), Parameter :: nd = id_ndept
  Integer (i4b), Parameter :: ifretot = id_frec1
  Integer (i4b), Parameter :: nd1 = id_depfi, nfmax = id_nfesc

  Real (dp), Parameter :: c1 = 0.02654D0, c2 = 3.97285D-16

! ..
! .. scalar arguments ..
  Real (dp) :: sr, srvmax, tem, tgrad, xcblue, xcred
  Integer (i4b) :: nb
  Logical :: escat
! ..
! .. array arguments ..
  Real (dp) :: gflu(nb), gl(nb), gu(nb), opacon(nd1, 2), opal(nd1, nb), &
    r(nd1), scont(nd1, 2), sline(nd1, nb), temp(nd1), xmax(nb), xmaxdop(nb), &
    xmaxdop_min(nb), xne(nd), xnefi(nd1), clf(nd), xnue(nb), opacnew(nd1, 2)

  Integer (i4b) :: index1(nd), nl(nb), nu(nb)
! ..
! .. local scalars ..
  Real (dp) :: aux, dummy, factor, oc, oc1, oc2, ol, ol1, ol2, oopa, rlj, &
    rlj1, rlj2, sc, sc1, sc2, tl, tl1, tl2, xkl, xl, xl1, xl2, xnl, xnu, xxx, &
    oopa1, taucon, depl, depu, dumfloat, xmed, limlow, limup, limlowt, limupt, &
    dx1, dx2, fmed

  Real (dp) :: tcl, fac2, eff_rat_int

  Integer (i4b) :: i, ik, j, j1, j2, k, k0, ki, l, li, nfre, nrec, nrecet, nx, &
    itauc1, ifre, dumint, imin, imax
! ..
! .. local arrays ..
  Real (dp) :: blevel(id_llevs), alevel(id_llevs), freqc(ifretot), &
    opac(nd, ifretot), scon(nd, ifretot), t(nd), tauross(nd), xnelte(nd), &
    lam(ifretot), logfnu(ifretot), rtau1(ifretot)

  Integer (i4b) :: no(ifretot)

  Real (dp), Dimension (nd, nb) :: eff_rat
  Real (dp), Dimension (nd, ifretot) :: opa_eff_rat

  Integer (i4b), Allocatable, Dimension (:) :: depthl
  Real (dp), Allocatable, Dimension (:) :: xmaxl, xlamb
  Real (dp), Allocatable, Dimension (:, :) :: xnll, xnul

  Logical :: exi

! ..
! .. external subroutines ..
  External :: calcxmax
! ..
! .. intrinsic functions ..
! ..
! cw----xmaxl: maximal 'linestrength'
! cw----depthl: depth of maximal 'linestrength'

! consistency check for approximate treatment (AT)
  If (.Not. at .And. nlevh21/=0) Stop ' NOT AT AND NLEVH21 NE 0'
  If (at) Then
    If (nlevh21==0) Stop ' AT AND NLEVH21 EQ 0'
    If (nl(1)/=id_llevs+1) Stop ' AT: ERROR IN NL'
    If (nu(1)/=id_llevs+2) Stop ' AT: ERROR IN NU'
    If (nb/=1) Stop ' AT: ERROR IN NB'
  End If

  Allocate (depthl(nb), xmaxl(nb), xlamb(nb))
  Allocate (xnll(nd,nb), xnul(nd,nb))

  nrec = id_llevs
  nrecet = 8*nrec

! now including OPA_EFF_RAT
  Open (1, File=trim(file)//'/CONT_FORMAL', Status='OLD', Form='UNFORMATTED')
  Rewind 1
  Read (1) nfre, (freqc(i), i=1, ifretot), ((scon(i,j),i=1,nd), j=1, ifretot), &
    ((opac(i,j),i=1,nd), j=1, ifretot), (t(i), i=1, nd), &
    ((opa_eff_rat(i,j),i=1,nd), j=1, ifretot)
  Close (1)
  Print *, nfre

! when individual lines are compared with those from the
! complete approach (when existent) they should be identical.
! Thus, SCON and OPAC (somewhat different from values above) need to be read
! in
! The file CONT_FORMAL_CMF only exists in the model if the last nlte-run had
! been performed with OPTCMF_FULL = T.
! If the last run was with OPTCMF_FULL = F, CONT_FORMAL_CMF (when previously
! existing) had been deleted.

  Inquire (File=trim(file)//'/CONT_FORMAL_CMF', Exist=exi)
  If (exi) Then
    Print *, ' SCON AND OPAC READ FROM CONT_FORMAL_CMF!'
    Open (1, File=trim(file)//'/CONT_FORMAL_CMF', Status='OLD', &
      Form='UNFORMATTED')
    Rewind 1
    Read (1)((scon(i,j),i=1,nd), j=1, ifretot), ((opac(i,j),i=1,nd), j=1, &
      ifretot)
    Close (1)
  End If

! CHECK CONSISTENCY
  Do i = 1, ifretot
    If (freqc(i)==0.) Exit
  End Do
  ifre = i
! here might be a bug. Sometimes, the frequency files are longer than
! NFRE (because of previous iterations with more freq. points. Valid
! data are only until nfre. Nevertheless, ifre only used here

! do i=1,nd
! do j=1,ifretot
! print*,i,j,opa_eff_rat(i,j)
! enddo
! enddo

  If (.Not. optthick .And. maxval(abs(opa_eff_rat(1:nd,1:ifre)- &
    1.D0))>1.D-6) Then
    Print *, maxval(abs(opa_eff_rat(1:nd,1:ifre)-1.D0))
    Stop ' OPTICALLY THIN CLUMPING AND OPA_EFF_RAT NE 1, subr. contin!!!'
  End If

  tem = t(nd)

  If (at) Then
    Open (1, File=trim(file)//'/TAU_ROS', Status='OLD', Form='FORMATTED')
    Rewind 1
    Do l = 1, nd
      Read (1, Fmt=*) tauross(l), xnelte(l)
    End Do
    Close (1)
  End If



! ---- xnue in cm-1
! ---- it is assumed without any check that xnue(1) > xnue(2) (i.e.,
! ---- xnue(1) corresponds to the blue component)

  If (xnue(1)>freqc(nfre) .Or. xnue(1)<freqc(1)) Stop &
    'FREQ. NOT CALCULATED IN CONTINUUM TRANSFER'

  k0 = nfre
  xxx = freqc(nfre)

  Do While (xnue(1)<xxx)
    k0 = k0 - 1
    xxx = freqc(k0)
  End Do

! check whether pseudo continuum is smooth
! around line, and define suitable cont. frequencies.
  lamtrans = 1.D8/xnue(1)
  Print *
  Print *, ' CHECK CONTINUUM BEHAVIOUR (VIA FLUXCONT)'

  Open (1, File=trim(file)//'/FLUXCONT', Status='OLD')
  Read (1, *) !                     FIRST LINE WITH COMMENTS
  Do i = 1, nfre
    Read (1, *, End=100, Err=100) no(i), lam(i), logfnu(i), dumfloat, &
      dumfloat, dumfloat, dumint, rtau1(i)
  End Do
  Close (1)
  Go To 110

100 Stop ' PROBLEMS IN READING FLUXCONT '

110 Continue
! PRINT*,LAM(1),LOGFNU(1),RTAU1(1)
! PRINT*,LAM(NFRE),LOGFNU(NFRE),RTAU1(NFRE)
! PRINT*,LAM(K0)
  If (k0<11 .Or. nfre-k0<11) Stop &
    ' LINE LOCATED TOO CLOSE TO WAVELENGTH BORDERS'

  Do i = k0, 1, -1
    If (lam(i)-lamtrans>20.) Exit
  End Do
  imin = i
! not for optical lines; otherwise, problem with Balmerjump
  If (lamtrans<uvlimit) imin = min(k0-10, imin)

  Do i = k0 + 1, nfre
    If (lamtrans-lam(i)>20.) Exit
  End Do
  imax = i
! not for optical lines; otherwise, problem with Balmerjump
  If (lamtrans<uvlimit) imax = max(k0+10, imax)

  Print *, ' TRANSITION AT ', lamtrans
  Print *, ' CONSIDERED WAVEL. RANGE FOR MEDIAN = ', lam(imax), lam(imin)

  Call median(logfnu(imin:imax), imax-imin+1, fmed)
  Print *, ' RANGE OF LOG FNUE AROUND LINE:'
  Print *, minval(logfnu(imin:imax)), maxval(logfnu(imin:imax))
  Print *, ' MEDIAN: ', fmed

! ALLOW FOR 0.05 DEX VARIATION IN LOGFNU
  limlow = fmed - 0.05
  limup = fmed + 0.05

  Call median(rtau1(k0-10:k0+10), 21, xmed)
  Print *, ' RANGE OF R(TAU=1) AROUND LINE:'
  Print *, minval(rtau1(k0-10:k0+10)), maxval(rtau1(k0-10:k0+10))
  Print *, ' MEDIAN: ', xmed

! ALLOW FOR 20% VARIATION IN R(TAU=1)
  limlowt = xmed - 0.2D0*(xmed-1.D0)
  limupt = xmed + 0.2D0*(xmed-1.D0)

  If (imax-imin==1) Then
    limlowt = min(limlowt, rtau1(imin), rtau1(imax))
    limupt = max(limupt, rtau1(imin), rtau1(imax))
    Print *, ' ALLOWED RANGE OF RTAU1-VALUES RESET TO ', limlowt, limupt
  End If

! PRINT*,LIMLOW,LIMUP
! FIND UPPER WAVELENGTH POINT KLOW LE K0
  Do i = k0, imin, -1
    If (logfnu(i)>=limlow .And. logfnu(i)<=limup .And. rtau1(i)>=limlowt .And. &
      rtau1(i)<=limupt) Go To 120
!   print*,k0,imin,imax,i
!   print*,logfnu(i),limlow,limup
!   print*,rtau1(i),limlowt,limupt
!   print*
  End Do
  Print *, ' NO SUITABLE CONTINUUM POINT FOUND IN THE UPPER WAVEL. RANGE'
  klow = -1
  kup = -1
  Deallocate (depthl, xmaxl, xlamb)
  Deallocate (xnll, xnul)
  Return

120 klow = i
! PRINT*,K0,KLOW,LAM(KLOW)

! FIND LOWER WAVELENGTH POINT KUP GT K0
  Do i = k0 + 1, imax
    If (logfnu(i)>=limlow .And. logfnu(i)<=limup .And. rtau1(i)>=limlowt .And. &
      rtau1(i)<=limupt) Go To 130
  End Do
  Print *, ' NO SUITABLE CONTINUUM POINT FOUND IN THE LOWER WAVEL. RANGE'
  klow = -1
  kup = -1
  Deallocate (depthl, xmaxl, xlamb)
  Deallocate (xnll, xnul)
  Return

130 kup = i
! PRINT*,K0,KUP,LAM(KUP)
! IF SEPARATION TOO LARGE IN THE UV, PERFORM NO INTERPOLATION,
! BUT USE CLOSER POINT'

  If (lamtrans<uvlimit) Then

    If (lam(klow)-lam(kup)>20.) Then
      dx1 = lam(klow) - lam(k0)
      dx2 = lam(k0) - lam(kup)
      If (dx1<dx2) Then
        kup = klow
      Else
        klow = kup
      End If

      If (min(dx1,dx2)>50.) Then
        Print *, ' NO SUITABLE FREQUENCY POINT WITHIN 50 A FOUND, &
          &LINE CANNOT BE CALCULATED'
        klow = -1
        kup = -1
        Deallocate (depthl, xmaxl, xlamb)
        Deallocate (xnll, xnul)
        Return
      Else If (min(dx1,dx2)>10.) Then
        Print *, ' NO SUITABLE FREQUENCY POINT WITHIN 10 A FOUND, &
          &PROCEED AT OWN RISK'
      End If
    End If

  Else !                            CONTINUUM AND IR SHOULD BEHAVE AS IN
!   PREVIOUS VERSIONS, OTHERWISE RESET
    If (klow/=k0 .Or. kup/=k0+1) Then
      Print *, ' WARNING! WARNING! WARNING!'
      Print *, ' OPTICAL/IR LINE, AND KLOW/KUP DIFFERENT FROM K0/K0+1'
      Print *, ' KLOW, KUP RESET TO K0, K0+1, PROCEED AT OWN RISK'
      Print *
      klow = k0
      kup = k0 + 1
!     KLOW=-1
!     KUP=-1
!     RETURN
    End If
  End If

  Print *, ' TRANSITION AT ', lamtrans
  Print *, ' CONTINUUM TAKEN (INTERPOLATED) IN BETWEEN'
  Print *, lam(kup), lam(klow), kup, klow

  If (klow/=k0 .Or. kup/=k0+1) Then
    Print *, ' WARNING! CONTINUUM POINTS MODIFIED,'
    Print *, ' DUE TO VARIATIONS IN FLUX AND R(TAU=1)'
  End If

  fcontup = logfnu(kup)
  fcontlow = logfnu(klow)

  Print *, ' CORRESPONDING LOG FNU VALUES FROM FLUXCONT', fcontup, fcontlow
  Print *


  Open (1, File=trim(file)//'/NLTE_POP', Status='OLD', Access='DIRECT', &
    Recl=nrecet)
  Open (11, File=trim(file)//'/LTE_POP', Status='OLD', Access='DIRECT', &
    Recl=nrecet)

! ---- calculation of cmf-frequencies for continuum interpolation
! ---- (they corresponds to the limits of the subinterval defined
! ---- by klow and kup)

  xcred = (freqc(klow)-xnue(1))*clight*srvmax/sr/xnue(1)
  xcblue = (freqc(kup)-xnue(1))*clight*srvmax/sr/xnue(1)

! ---- lambda in cm

  Do i = 1, nb
    xlamb(i) = 1.D0/xnue(i)
  End Do

  Do i = 1, nb
    xmaxl(i) = 0.D0
  End Do

  inverted = 0

  taucon = 0.
  itauc1 = 0

  eff_rat = 0.

  Do i = 1, nd

    Read (1, Rec=i)(blevel(j), j=1, nrec)
    Read (11, Rec=i)(alevel(j), j=1, nrec)
!   recheck
    scont(index1(i), 1) = scon(i, klow)
    scont(index1(i), 2) = scon(i, kup)
    opac(i, klow) = opac(i, klow)/opa_eff_rat(i, klow)
    opac(i, kup) = opac(i, kup)/opa_eff_rat(i, kup)
!   JO Jan 2023 not required here, calculated below
!   OPACON(INDEX1(I),1) = OPAC(I,KLOW)*SR
!   OPACON(INDEX1(I),2) = OPAC(I,KUP)*SR
!   OPAC IS ALREADY EFFECTIVE QUANTITY, so needs to be back-corrected.
!   Then follow philosophy  of correcting with largest TCL
!   -- either line or from blended line+cont BG
!   (i.e. opac is corrected again below)


!   ASSUMING RED OPACITY IS HIGHER
!   IF(I.LT.ND.AND.ITAUC1.EQ.0) THEN
!   TAUCON=TAUCON+ &
!   &      .5*(OPACON(INDEX1(I),1)+OPAC(I+1,KLOW)*SR)* &
!   &      (R(INDEX1(I))-R(INDEX1(I+1)))
!   IF(TAUCON.GE.0.66) ITAUC1=I
!   ENDIF
!   Do this calculation AFTER it's been decided which TCL reduction to use

    Do ik = 1, nb

      If (.Not. at) Then
!       standard treatment
        xnl = blevel(nl(ik))
        If (xnl==0.) xnl = 1.D-15*alevel(nl(ik))

        xnu = blevel(nu(ik))
        If (xnu==0.) xnu = 1.D-15*alevel(nu(ik))
      Else
!       approximate treatment
        If (optdep) Then
          aux = alevel(nlevh21)/blevel(nlevh21)*xnelte(i)/xne(i)
!         LAST TWO DEPARTURES, CORRESPONDING TO H119 and H120 in A10HHe
          depl = blevel(nlevh21-2)/alevel(nlevh21-2)*aux
          depu = blevel(nlevh21-1)/alevel(nlevh21-1)*aux
        End If
        Call occup_at(nl(ik), nu(ik), xne(i), t(i), blevel(nlevh21), xnl, xnu)
        If (optdep) Then
          xnl = xnl*depl
          xnu = xnu*depu
        End If
      End If

      xnll(i, ik) = xnl
      xnul(i, ik) = xnu

      xkl = c1*gflu(ik)*(xnl/gl(ik)-xnu/gu(ik))
      sline(index1(i), ik) = c2*(xnue(ik)**3)/(xnl/xnu*gu(ik)/gl(ik)-1.D0)

!     CAN BE MODIFIED FOR TESTS
!     Note: in most cases, this should give a very similar result to the
!     "exact"
!     approach (as long as tau_inv small, i.e., a dominating source term)
!     IN CASE, DON'T FORGET TO INSERT THE ABS BELOW (FOR INTERPOLATED VALUES!
!     XKL = ABS(C1*GFLU(IK)* (XNL/GL(IK)-XNU/GU(IK)) )
!     SLINE(INDEX1(I),IK) = ABS(C2*
!     (XNUE(IK)**3)/(XNL/XNU*GU(IK)/GL(IK)-1.D0))

      If (xkl<0.) Then
!       that was the old formulation; new one allows for a better (though not
!       perfect) treatment of inverted levels
!       XKL = C1*GFLU(IK)* XNL/GL(IK)
!       SLINE(INDEX1(I),IK)=0.
        inverted(i, ik) = 1
      End If
!     here comes the correction (only for lines!)
      xkl = xkl/clf(i)
      opal(index1(i), ik) = xkl*sr
!     correction for thick clumping + consistent inversion
      tcl = abs(xkl)*xlamb(ik)*tcl_fac_line(i)
!     Note: SRVMAX term included in tcl_fac_line multiplier
      fac2 = (1.+fic(i)*tcl)/(1.+tcl)
      eff_rat(i, ik) = min(fac2, opa_eff_rat(i,klow), opa_eff_rat(i,kup))
      If (eff_rat(i,ik)<=0.0D0 .Or. eff_rat(i,ik)>1.0D0) Then
        Print *, eff_rat(i, ik)
        Print *, tcl, fac2, i, ik
        Stop ' OPA_EFF > <OPA>, CONTIN'
      End If
      opal(index1(i), ik) = opal(index1(i), ik)*eff_rat(i, ik)
!     Now effective line opacity, reduce continuum with *same* factor
!     JO JAN 2023: HERE WAS A BUG, CHANGED
      opacnew(i, 1) = opac(i, klow)*eff_rat(i, ik)
!     next statement actually not required, since only used for oopa at KLOW
      opacnew(i, 2) = opac(i, kup)*eff_rat(i, ik)
!     JO JAN 2023: not required here since calculated below
!     OPACON(INDEX1(I),1) = OPAC(I,KLOW)*SR
!     OPACON(INDEX1(I),2) = OPAC(I,KUP)*SR

      oopa = abs(xkl*eff_rat(i,ik))/opacnew(i, 1)

!     cw----search for maximum linestrength

      If (i>=1) Then
        If (oopa>xmaxl(ik)) Then
          xmaxl(ik) = oopa
          depthl(ik) = i
        End If
      End If
    End Do

  End Do

! JO JAN 2023: FOLLOWING DO-LOOP IS NEW, TO CURE PREVIOUS BUG
! also: now opacon is original opac*SR, to be consistent with
! radiative transfer in nlte.
! in previous (buggy) version, opacon was somewhat erroneous

  Do i = 1, nd
!   THE NEXT STATEMENT -- IN LINE WITH ORIGINAL APPROACH -- WOULD BE AN
!   APPROX.,
!   AND MIGHT OVERESTIMATE THE EFFECT (SINCE ONLY VALID CLOSE TO THE STRONGEST
!   LINE).
!   EFF=MINVAL(EFF_RAT(I,1:NB))
!   OPAC(I,KLOW)=OPAC(I,KLOW)*EFF
!   OPAC(I,KUP)=OPAC(I,KUP)*EFF
!   TO BE ON THE SAFE SIDE, WE USE THE SAME CORRECTED CONTINUUM AS IN NLTE
!   NEEDS TO BE CAREFULLY TESTED (CONSISTENCY TEST IS ALREADY DONE).
    opac(i, klow) = opac(i, klow)*opa_eff_rat(i, klow)
    opac(i, kup) = opac(i, kup)*opa_eff_rat(i, kup)
    opacon(index1(i), 1) = opac(i, klow)*sr
    opacon(index1(i), 2) = opac(i, kup)*sr
  End Do

! Need to calculate TAUCON scale here, since might have been modified,
! and OPAC(I+1) only known after end of previous I loop
! ASSUMING RED OPACITY IS HIGHER
  Do i = 1, nd
    If (i<nd .And. itauc1==0) Then
      If (opacon(index1(i+1),1)/=opac(i+1,klow)*sr) Stop &
        ' ERROR IN OPACON/OPAC -- subr. CONTIN'
      taucon = taucon + .5*(opacon(index1(i),1)+opac(i+1,klow)*sr)*(r(index1(i &
        ))-r(index1(i+1)))
      If (taucon>=0.66) itauc1 = i
    End If
  End Do

  Close (1)
  Close (11)

  If (itauc1==0) Stop ' TAUC = 0.66 NOT FOUND'
  Print *
  Print *, ' TAUCON = ', taucon, ' AT LOG NE = ', log10(xne(itauc1))

! cw----calculation of xmax at depth of maximal linestrength

  Do ik = 1, nb

    oopa1 = abs(opal(index1(itauc1),ik))/opacon(index1(itauc1), 1)
!   both opacities (line and cont) corrected

    Print *
    Print *, 'DEPTHL =', depthl(ik)
    Print *, 'LOG(XMAXL)  =', log10(xmaxl(ik))
!   calcxmax remains unaffected, since we assume that profile is
!   formed inside clumps, i.e., at increased densities
    Call calcxmax(xmax(ik), xmaxdop(ik), xmaxdop_min(ik), ik, t(depthl(ik)), &
      xmaxl(ik), xlamb(ik), xne(depthl(ik)), gammal(ik), gammac(ik), &
      t(itauc1), xne(itauc1))
    Print *, 'XMAX  =', xmax(ik)
  End Do

  Do ik = 1, nb
    If (xmax(ik)==0.D0) Stop 'ERROR IN XMAX - CONTIN'
  End Do

! ---- interpolation (linear in log-log) in the micro-grid
! -----interpolate *linearly* in EFF_RAT, to avoid certain problems

  temp(1) = t(1)
  xnefi(1) = xne(1)
iloop: Do i = 2, nd
    j1 = index1(i-1)
    j2 = index1(i)
    temp(j2) = t(i)
    xnefi(j2) = xne(i)

jloop: Do j = j1 + 1, j2 - 1
      rlj = log(r(j))
      rlj1 = log(r(j1))
      rlj2 = log(r(j2))
      factor = (rlj-rlj1)/(rlj2-rlj1)
      tl1 = log(temp(j1))
      tl2 = log(temp(j2))
      tl = tl1 + factor*(tl2-tl1)
      temp(j) = exp(tl)
      xl1 = log(xnefi(j1))
      xl2 = log(xnefi(j2))
      xl = xl1 + factor*(xl2-xl1)
      xnefi(j) = exp(xl)
      Do k = 1, 2
        sc1 = log(scont(j1,k))
        sc2 = log(scont(j2,k))
        sc = sc1 + factor*(sc2-sc1)
        scont(j, k) = exp(sc)
        oc1 = log(opacon(j1,k))
        oc2 = log(opacon(j2,k))
        oc = oc1 + factor*(oc2-oc1)
        opacon(j, k) = exp(oc)
      End Do

      Do k = 1, nb
!       linear in EFF_RAT_INT
        rlj = r(j)
        rlj1 = r(j1)
        rlj2 = r(j2)
        fac2 = (rlj-rlj1)/(rlj2-rlj1)
        ol1 = eff_rat(i-1, k)
        ol2 = eff_rat(i, k)
        eff_rat_int = ol1 + fac2*(ol2-ol1)

        ol1 = log(xnll(i-1,k)/clf(i-1))
        ol2 = log(xnll(i,k)/clf(i))
        ol = ol1 + factor*(ol2-ol1)
        xnl = exp(ol)
        ol1 = log(xnul(i-1,k)/clf(i-1))
        ol2 = log(xnul(i,k)/clf(i))
        ol = ol1 + factor*(ol2-ol1)
        xnu = exp(ol)
        opal(j, k) = c1*gflu(k)*(xnl/gl(k)-xnu/gu(k))*sr*eff_rat_int
!       Now effective
        If (opal(j,k)==0.) Stop ' OPAL = 0 IN INTERPOLATION; CHANGE SCHEME'
        sline(j, k) = c2*(xnue(k)**3)/(xnl/xnu*gu(k)/gl(k)-1.D0)
!       IN CASE (IF ABS HAS BEEN USED ABOVE
!       OPAL(J,K)  = ABS(C1*GFLU(K)* (XNL/GL(K)-XNU/GU(K))*SR)
!       SLINE(J,K) = ABS(C2* (XNUE(K)**3)/(XNL/XNU*GU(K)/GL(K)-1.D0))
        If (sline(j,k)*opal(j,k)<0.) Then
          Print *, r(j), ' ', j, ' ', k, ' OPAL*SL < 0!'
          Stop ' OPAL*SL < 0!'
        End If
      End Do
    End Do jloop

  End Do iloop

  tgrad = abs(t(nd)-t(nd-1))/(r(index1(nd-1))-r(index1(nd)))

  Deallocate (depthl, xmaxl, xlamb)
  Deallocate (xnll, xnul)

  If (.Not. escat) Return

! interpolation of scont (incl. non-coherent el.scat term)
! on micro grid

! !
! open(1,file=TRIM(FILE)//'/exi_'//line,status='unknown',form='formatted')

  Rewind 17
  Read (17, Fmt=*) nx, nfcmf

  If (nx/=nd) Stop 'NX .NE. ND'
  If (nfcmf>nfmax) Stop 'NFCMF > NFMAX'

  Do k = 1, nfcmf
    Read (17, Fmt=*) ki, xcmfe(k), dummy
  End Do

  Do k = 1, nfcmf
    Do l = 1, nd
      Read (17, Fmt=*) ki, li, dummy, sconte(index1(l), k), dummy
    End Do
  End Do

  Close (17)

  Do i = 2, nd
    j1 = index1(i-1)
    j2 = index1(i)
    Do j = j1 + 1, j2 - 1
      rlj = log(r(j))
      rlj1 = log(r(j1))
      rlj2 = log(r(j2))
      factor = (rlj-rlj1)/(rlj2-rlj1)
      Do k = 1, nfcmf
        sc1 = log(sconte(j1,k))
        sc2 = log(sconte(j2,k))
        sc = sc1 + factor*(sc2-sc1)
        sconte(j, k) = exp(sc)
      End Do
    End Do
  End Do

  Return
! WRITE(*,1001) J1,J2,J,OPAL(J1,1),OPAL(J,1),OPAL(J2,1)
! WRITE(*,1001)J1,J2,J,SLINE(J1,1),SLINE(J,1),SLINE(J2,1)
140 Format (3(I4,2X), 3(E9.3,2X))

End Subroutine

!-----------------------------------------------------------------------

Subroutine calcxmax(xmax, xmaxdop, xmaxdop_min, ik, t, opa, xlambda, xne, &
  gammal, gammac, tetauc1, xnetauc1)

  Use :: nlte_type
  Use :: nlte_dim
  Use :: fund_const, Only: akb, xmh => amh, pi, clight
  Use :: formalsol_var, Only: nstark, weig, vmax, vturb, vturbmin, vturbmax
  Implicit None

! ---- calculates xmax for component 'ik'
! ---- now with vturbmax for Doppler;
! ---- for STARK and Voigt, calculate with vturb=vturbmin (photosphere),
! -----and then take maximum (xmax, xmaxdop) except for const. vturb

! .. parameters ..
  Real (dp), Parameter :: wpi1 = 0.56418958D0
! ..
! .. scalar arguments ..
  Real (dp) :: gammac, gammal, opa, t, xlambda, xmax, xmaxdop, xmaxdop_min, &
    xne, tetauc1, xnetauc1
  Integer (i4b) :: ik
! ..
! .. local scalars ..
  Real (dp) :: aux, avoigt, eps, vther, xkl, phii, xi, ximin
  Logical :: flag
! ..
! ..
! .. external functions ..
  Real (dp) :: perfil1
! ..

  If (vturb/=vturbmin) Stop ' VTURB NE VTURBMIN in CALCXMAX'

! minimum for xmax according to dlam = 10. A or 0.3 vinf

  ximin = 10.D-8/xlambda*clight/vmax
  ximin = max(ximin, 0.3D0)

! new: calculated XMAXDOP_MIN as before, with VTURB=VTURBMIN,
! to be used in core region of frequency grid
  vther = sqrt(2.D0*akb*t/xmh/weig(ik)+vturb**2)

! ------- doppler profile(overall or in wind)

  xkl = opa*xlambda/vther
  If (xkl>1.D0) Then

!   ##----  modifiied by jp, assumimg line profile 1/100 at xmax

    aux = xkl*wpi1*100.D0
    aux = log(aux)
    xmaxdop_min = vther/vmax*sqrt(aux)
  Else
    xmaxdop_min = max(1.D-3, 3.D0*vther/vmax)
  End If

! new
! NOW, for XMAX, use VTURBMAX instead of VTURB;
! remains identical to old version for constant vturb,
! since then VTURBMAX=VTURBMIN=VTURB
  vther = sqrt(2.D0*akb*t/xmh/weig(ik)+vturbmax**2)

! ------- doppler profile(overall or in wind)

  xkl = opa*xlambda/vther
  If (xkl>1.D0) Then

!   ##----  modifiied by jp, assumimg line profile 1/100 at xmax

    aux = xkl*wpi1*100.D0
    aux = log(aux)
    xmaxdop = vther/vmax*sqrt(aux)
  Else
    xmaxdop = max(1.D-3, 3.D0*vther/vmax)
  End If

  If (nstark(ik)==0) Then
    xmax = xmaxdop
  Else If (nstark(ik)==1) Then

!   stark profile: find x where stark profile is below threshold
!   note that doppler width is included into perfil

!   THUS: opal*perfil(xmax) le eps*opac gives desired xmax

!   first for tauc = 0.66

    eps = 5.D-3 !                   previous version 1.d-2

!   AUX=EPS/OPA1 !old version
!   in the new version, we use the maximum of the line-strength, to be on the
!   safe side (together with the conditions at tauc=0.66)
    aux = eps/opa

    xi = 0.05D0
    phii = perfil1(ik, xi, tetauc1, xnetauc1, vturb, flag)

    Do While (phii>aux)
      xi = xi + 0.05D0
      phii = perfil1(ik, xi, tetauc1, xnetauc1, vturb, flag)
    End Do

    xmax = max(xi, ximin)
!   JO: might need to be changed
    If (vturbmax>vturbmin) xmax = max(xmax, xmaxdop)
  Else

!   cw------ voigt profil

    If (nstark(ik)/=2 .And. nstark(ik)/=3) Stop ' ERROR IN NSTARK'
    avoigt = (gammal+xne*gammac)*xlambda/vther
    aux = xkl*100.*avoigt/pi
    xmax = vther/vmax*sqrt(aux)

    xmax = max(xmax, ximin)
!   JO: might need to be changed
    If (vturbmax>vturbmin) xmax = max(xmax, xmaxdop)
  End If

  Return
End Subroutine

!-----------------------------------------------------------------------

Subroutine model(r, rprim, v, vprim, dvdr, rhopri, rmax, nd1, nd2, vmax, vdop, &
  index1, srvmax, sr, nb, xne, clf, clf_test, vtmi, vtma)
! !
! enhancement factor clf now read from 'model'
! NOTE: ne, nh and occupation numbers include clf, rho NOT
  Use :: nlte_type
  Use :: nlte_dim
  Use :: fund_const, Only: akb, xmh => amh
  Use :: formalsol_var, Only: vturb, file, vturbmin, vturbmax, vturbv

  Implicit None

! ---- enrique santolaya, 23 de junio de 1993 (noche de san juan)

! ---- modified (5-12-94): no hamman velocity law is used. instead of
! ---- that, the velocity is linearly interpolated in the radial grid to
! ---- obtain the micro-grid.
! ---- any velocity law can then be used (e. santolaya).

! ---- file 'model' is read. it contain all information about stellar
! ---- parameters
! ##   xne only approximative in 'model', actual value from 'enion'
! ##   has to be used

! .. parameters ..
  Integer (i4b), Parameter :: nd = id_ndept, ndm = id_depfi
  Integer (i4b), Parameter :: kel = id_atoms, kis = id_kisat
  Real (dp), Parameter :: rsun = 6.96D10
! ..
! .. scalar arguments ..
  Real (dp) :: rmax, sr, srvmax, vmax, vtmi, vtma
  Integer (i4b) :: nb, nd1, nd2
! ..
! .. array arguments ..
  Real (dp) :: r(nd), rhopri(ndm), rprim(ndm), vdop(nb), vprim(ndm), xne(nd), &
    clf(nd), v(nd), dvdr(nd)
  Real (dp), Dimension (nd) :: clf_test

  Integer (i4b) :: index1(nd)
! ..
! .. local scalars ..
  Real (dp) :: a, aux, beta, cteceq, del, dvmin, dvrho, ggrav, vmin, vvdop, &
    vd, x, xmloss, xmu, yhe, teff
  Integer (i4b) :: i, j, k, nn, nrho, nv
! ..
! .. local arrays ..
  Real (dp) :: enionnd(kel, kis+1, nd), rho(nd), xnh(nd)

  Logical :: range
! ..
! .. intrinsic functions ..
! ..

  Open (1, File=trim(file)//'/MODEL', Status='OLD', Form='UNFORMATTED')
  Rewind 1
  Read (1) teff, ggrav, sr, yhe, xmu, vmax, xmloss, beta, (r(i), i=1, nd), &
    (v(i), i=1, nd), (dvdr(i), i=1, nd), (rho(i), i=1, nd), (xne(i), i=1, nd), &
    (xnh(i), i=1, nd), (clf(i), i=1, nd), (index1(i), i=1, nd)
  Close (1)

  Open (1, File=trim(file)//'/ENION', Status='OLD', Form='UNFORMATTED')
  Read (1) enionnd, xne
  Close (1)

! test consistency
  Do i = 1, nd
!   lower precision in CLUMPING_OUTPUT
    If (abs(clf_test(i)-clf(i))>1.D-3) Stop &
      ' Error in in clumping params!, model_1'
  End Do

  rmax = r(1)

! if vdop(nb) = 0., we are in the "range"-branch (here, nb=1)
  range = vdop(nb) == 0.


! ---- calculation of doppler width(s) using VTURBMIN

  vturbmin = vtmi*1.D5
  If (vturbmin/=vturb) Stop ' SOMETHING WRONG WITH VTURB -- SUBR. MODEL'
  If (vtma<=1.D0) Then
    vturbmax = vtma*vmax
  Else
    vturbmax = vtma*1.D5
  End If

  If (vturbmax<vturbmin) Stop ' VTURBMAX < VTURBMIN, NOT ALLOWED!!'

  Print *, ' MIN(VTURB) = ', vturbmin/1.D5, ' MAX(VTURB) = ', vturbmax/1.D5
  Print *

  vvdop = 1.D50
  a = sqrt(2.D0*akb*teff/xmh)
  Do i = 1, nb
    vd = a*vdop(i)
    vdop(i) = sqrt(vd**2+vturb**2)/vmax ! using vturb=vturbmin
    vvdop = min(vvdop, vdop(i))
  End Do

  srvmax = sr/vmax

  vmin = v(nd)*vmax

  dvmin = min(vvdop/5.D0, 1.D-2)

! -----interpolation on velocity values

! ---- calculation of micro-grid using VTURBMIN (since min(dvmin 0.01 anyway)

  dvrho = 5.D-2
  rprim(1) = r(1)
  vprim(1) = v(1)
  rhopri(1) = rho(1)
  cteceq = rho(1)*v(1)*r(1)*r(1)

  j = 1

  Do i = 2, nd
    nv = int((v(i-1)-v(i))/dvmin) + 1
    nrho = int(log10(rho(i)/rho(i-1))/dvrho) + 1
    nn = max(nv, nrho)
    x = log10(r(i-1)/r(i))/dble(nn+1)
    a = (v(i-1)-v(i))/(r(i-1)-r(i))

    Do k = 1, nn
      j = j + 1
      aux = log10(r(i-1)) - dble(k)*x
      rprim(j) = 10.D0**aux
      del = rprim(j) - r(i)
      vprim(j) = v(i) + a*del
      rhopri(j) = cteceq/vprim(j)/rprim(j)/rprim(j)
    End Do

    j = j + 1
    index1(i) = j
    rprim(j) = r(i)
    vprim(j) = v(i)
    rhopri(j) = rho(i)
  End Do

  Print *
  Print *, 'RSTAR= ', sr/rsun, '    R-STAR/VMAX=', srvmax
  Print *
  Write (*, Fmt=100) xmloss, v(1)*vmax/1.D5, v(nd)*vmax/1.D5, rho(nd)
  Print *
  nd2 = j
  Print *, 'TOTAL NUMBER OF DEPTH POINTS=', nd2

  If (nd2>nd1) Stop 'TOO MANY DEPTH POINTS'

! CALCULATION OF DEPTH DEPENDENT VTURB = VTURBV propto V(I) and VDOPV

  If (.Not. range) Then
!   vturbv defined on micro-grid
    Allocate (vturbv(nd2))
    If (vturbmax/=vturbmin) Then
      Do i = 1, nd2
        vturbv(i) = max(vturbmin, vprim(i)*vturbmax)
!       PRINT*,I,VTURBV(I)/1.D5
      End Do
    Else
      vturbv = vturbmin
    End If

  Else
!   vturbv defined on macro-grid
    Allocate (vturbv(nd))
    If (vturbmax/=vturbmin) Then
      Do i = 1, nd
        vturbv(i) = max(vturbmin, v(i)*vturbmax)
      End Do
    Else
      vturbv = vturbmin
    End If

  End If

  Return

100 Format ('MLOSS/YEAR=', E12.6, '   VMAX(RMAX)=', F7.1, /, '    VMIN=', &
    F11.7, '    RHOPHO=', E13.6)

End Subroutine

!-----------------------------------------------------------------------

Subroutine staread(file, nb, nsum, nl, nu)

  Use :: nlte_type
  Use :: nlte_dim
  Use :: formalsol_var, Only: nstark, nlevl, nlevu, nws, nts, nes, qhalf, dws, &
    ts, es, ps

  Implicit None

! ---- it reads and stores  stark profiles from binary file 'stark.bin'

! .. parameters ..
  Integer (i4b), Parameter :: maxw = id_maxww, maxt = id_maxtt, &
    maxne = id_maxne, maxps = id_maxps
! ..
! .. scalar arguments ..
  Integer (i4b) :: nb, nsum
! ..
! .. array arguments ..
  Integer (i4b) :: nl(nb), nu(nb)
! ..
! .. local scalars ..
  Integer (i4b) :: i, j, nn
! ..
  Character (60) :: file

  Open (1, File=trim(file)//'/STARK.BIN', Status='OLD', Form='UNFORMATTED')
  Rewind 1

  If (nsum==0. .Or. nsum>nb) Stop ' STAREAD CALLED ERRONEOUSLY'

  Do i = 1, nb

    If (nstark(i)/=1) Cycle

    Read (1, End=100, Err=110) nlevl(i), nlevu(i), qhalf(i), nws(i), &
      (dws(j,i), j=1, nws(i)), nts(i), (ts(j,i), j=1, nts(i)), nes(i), &
      (es(j,i), j=1, nes(i)), nn, (ps(j,i), j=1, nn)

!   ---- checking data

    If (nws(i)>maxw) Stop 'TOO MANY FREQUENCIES'
    If (nts(i)>maxt) Stop 'TOO MANY TEMPERATURES'
    If (nes(i)>maxne) Stop 'TOO MANY ELECTRON DENSITIES'
    If (nn>maxps) Stop 'TOO MANY POINTS IN PROFILE'
    If (nlevl(i)/=nl(i)) Stop 'LOWER LEVEL NOT CORRECT'
    If (nlevu(i)/=nu(i)) Stop 'UPPER LEVEL NOT CORRECT'
  End Do

  Close (1)

  Return

! ---- end of file not expected

100 Continue
  Print *, 'END OF FILE STARK.BIN NOT EXPECTED'
  Stop

! ---- error in reading

110 Continue
  Print *, 'ERROR IN READING STARK.BIN'
  Stop

End Subroutine

!***********************************************************************

!subroutines: complex ones
!elscat and related

!***********************************************************************

Subroutine elscat(xnue, nl, nu, gflu, gl, gu, nb)

  Use :: nlte_type
  Use :: nlte_dim
  Use :: fund_const, Only: xmh => amh, clight, sigmae
  Use :: formalsol_var, Only: vmax, xnemin, vturb, file, optthick, fic, &
    tcl_fac_line
  Implicit None

! calculates consistent mean intensities for non-coherent electron
! scattering on the macro grid. transfer in the cmf.
! el. scattering treated as in rybicki and hummer, aa 290, 553.

! CLUMPING INCLUDED (THIN AND THICK)

! .. parameters ..
  Integer (i4b), Parameter :: nd = id_ndept, np = id_npoin
  Integer (i4b), Parameter :: kel = id_atoms, kis = id_kisat
  Integer (i4b), Parameter :: ifretot = id_frec1, nfmax = id_nfesc

  Real (dp), Parameter :: c1 = 0.02654D0, c2 = 3.97285D-16
! ..
! .. scalar arguments ..
  Integer (i4b) :: nb
! ..
! .. array arguments ..
  Real (dp) :: gflu(nb), gl(nb), gu(nb), xnue(nb)
  Integer (i4b) :: nl(nb), nu(nb)
! ..
! .. local scalars ..
  Real (dp) :: a1, a2, aic, aux, b1, b2, beta, betat, const, corrfc, dbdr, dx, &
    dxcoarse, dxfine, err, ggrav, qblue, qred, rel, scnew, scold, sr, srvmax, &
    teff, vthel, xcblue, xcred, xkl, xlambda, xmloss, xmu, xnl, xnu, xxx, yhe, &
    z2, thomson_1, thomson_2, thomson_11, thomson_21, struek, struek1, lammin, &
    lammax, vref, lam_cmf, xmax, xmin, fac1, fac2, xxdiff, phikl, dum1, tcl

  Integer (i4b) :: i, iit, ik, j, jp, k, k0, l, lmax, lz, nc, nfcmf, nfre, &
    nrec, nrecet, k01, idv, inext, ilast, idx, ifre
  Logical :: flag
! ..
! .. local arrays ..
  Real (dp) :: blevel(id_llevs), alevel(id_llevs), dvdr(nd), &
    enionnd(kel, kis+1, nd), freqc(ifretot), p(np), r(nd), z(nd, np), v(nd), &
    t(nd), xne(nd), clf(nd), dum(nd), opac(nd), opacon(nd, ifretot), &
    st(nd, ifretot), strue(nd), thoms(nd, ifretot), thomson(nd), &
    ubluwi(nd, np-1), ucmf(nd), vcmf(nd-1), w0(nd, np-1), w2(nd, np-1), &
    xjall(nd, ifretot), xjold(nd), xx(nfmax), opa_eff_rat(nd, ifretot)

  Real (dp), Allocatable :: aick(:), dbdrk(:), eold(:), x(:), xi(:), fi(:), &
    etac(:, :), aj(:, :), e(:, :), xjnew(:, :), opal_tot(:, :), etal_tot(:, :)

  Real (dp), Allocatable :: opal(:, :), sline(:, :), xlamb(:), xc(:)

! ..
! .. external functions ..
  Real (dp) :: perfil1
  External :: perfil1
! ..
! .. external subroutines ..
  External :: cont1, diffus, ray, scatsol
! ..
! .. intrinsic functions ..
! ..
! .. data statements ..

  Data a1, a2/1.69070371729D0, -.69070371729D0/
  Data b1, b2/1.614249968779D0, 2.154326524957D0/
! ..

  Allocate (opal(nd,nb), sline(nd,nb), xlamb(nb), xc(nb))

  xnemin = 10.5D0
  xnemin = 10.D0**xnemin

  nrec = id_llevs
  nrecet = 8*nrec

! this subroutine assumes in any case and without checking
! that index1(i)=: i

  Open (1, File=trim(file)//'/MODEL', Status='OLD', Form='UNFORMATTED')
  Rewind 1
  Read (1) teff, ggrav, sr, yhe, xmu, vmax, xmloss, beta, (r(i), i=1, nd), &
    (v(i), i=1, nd), (dvdr(i), i=1, nd), (dum(i), i=1, nd), (dum(i), i=1, nd), &
    (dum(i), i=1, nd), (clf(i), i=1, nd)
  Close (1)

  Open (1, File=trim(file)//'/ENION', Status='OLD', Form='UNFORMATTED')
  Read (1) enionnd, xne
  Close (1)

  srvmax = sr/vmax

! ALL QUANTITIES ALREADY CORRECTED FOR CLUMPING
  Open (1, File=trim(file)//'/CONT_FORMAL_ALL', Status='OLD', &
    Form='UNFORMATTED')
  Rewind 1
  Read (1) nfre, (freqc(i), i=1, ifretot), ((opacon(i,j),i=1,nd), j=1, ifretot &
    ), ((thoms(i,j),i=1,nd), j=1, ifretot), ((st(i,j),i=1,nd), j=1, ifretot), &
    ((xjall(i,j),i=1,nd), j=1, ifretot), (t(i), i=1, nd), &
    ((opa_eff_rat(i,j),i=1,nd), j=1, ifretot)
  Close (1)

! CHECK CONSISTENCY
  Do i = 1, ifretot
    If (freqc(i)==0.) Exit
  End Do
  ifre = i
! here might be a bug. Sometimes, the frequency files are longer than
! NFRE (because of previous iterations with more freq. points. Valid
! data are only until nfre. Nevertheless, ifre only used here

  If (.Not. optthick .And. maxval(abs(opa_eff_rat(1:nd,1:ifre)- &
    1.D0))>1.D-6) Stop &
    ' OPTICALLY THIN CLUMPING AND OPA_EFF_RAT NE 1, subr. elscat!!!'


  Open (1, File=trim(file)//'/NLTE_POP', Status='OLD', Access='DIRECT', &
    Recl=nrecet)
  Open (11, File=trim(file)//'/LTE_POP', Status='OLD', Access='DIRECT', &
    Recl=nrecet)

! ---- xnue in cm-1
! ---- it is assumed without any check that xnue(1) > xnue(2) (i.e.,
! ---- xnue(1) corresponds to the blue component)


  If (xnue(1)>freqc(nfre) .Or. xnue(1)<freqc(1)) Stop &
    'FREQ. NOT CALCULATED IN CONTINUUM TRANSFER'

  k0 = nfre
  xxx = freqc(nfre)

  Do While (xnue(1)<xxx)
    k0 = k0 - 1
    xxx = freqc(k0)
  End Do

  k01 = k0 + 1
! commented out by JP (June 2013), since only bluest component decisive
! XXX = FREQC(K01)
! DO WHILE (XNUE(NB).LT.XXX)
! K0=K0-1
! XXX=FREQC(K0)
! END DO

! PRINT*,1.d8/FREQC(K01),1.d8/XNUE(1),1.d8/XNUE(NB),1.d8/XXX

! ---- calculation of cmf-frequencies for continuum interpolation
! ---- (they correspond to the limits of the subinterval defined
! ---- by k0 and k0+1)

! NOTE: from here on, all frequencies displayed as X (i.e., in units of vinf)
! are calculated with respect to the bluemost line frequency (XNUE(1))

  xcred = (freqc(k0)-xnue(1))*clight*srvmax/sr/xnue(1)
  xcblue = (freqc(k01)-xnue(1))*clight*srvmax/sr/xnue(1)

! ---- lambda in cm

  Do i = 1, nb
    xlamb(i) = 1.D0/xnue(i)
  End Do

! up to now linear interpolation to line freq of bluemost componnent (x=0),
! to be consistent with treatment in formacon

  dx = xcblue - xcred
  qblue = 1.D0 - xcblue/dx
  qred = 1.D0 - qblue

  Do i = 1, nd

!   recalculate thomson factor, separating electron-scat. and bg-line
!   contribution
!   re-correct for OPA_EFF_RAT
    thomson_1 = (xne(i)*sigmae*xmh/clf(i))/opacon(i, k0)*opa_eff_rat(i, k0)
    thomson_2 = thoms(i, k0) - thomson_1
    If (thomson_2<-1.D-14) Stop ' ERROR IN THOMSON_1(K0)'
    thomson_2 = max(0.D0, thomson_2)

    thomson_11 = (xne(i)*sigmae*xmh/clf(i))/opacon(i, k01)*opa_eff_rat(i, k01)
    thomson_21 = thoms(i, k01) - thomson_11
    If (thomson_21<-1.D-14) Stop ' ERROR IN THOMSON_11(K1)'
    thomson_21 = max(0.D0, thomson_21)

    opac(i) = (qred*opacon(i,k0)+qblue*opacon(i,k01)) ! already corrected
    dum1 = (qred*opacon(i,k0)/opa_eff_rat(i,k0)+qblue*opacon(i,k01)/ &
      opa_eff_rat(i,k01))
    thomson(i) = (xne(i)*sigmae*xmh/clf(i))/dum1
    opac(i) = opac(i)*sr

!   --------possibility that thoms > 1 due to updated electr. densities
!   in nlte code.

    If (thomson(i)>1.D0) Then
      Print *, 'WARNING THOMSON > 1 ', i, ' ', thomson(i)
      thomson(i) = 1.D0
    End If
!   correct strue for bg-line contribution
!   since THOMSON_2 has been consistently calculated, no additional correction
!   required
    struek = st(i, k0) + thomson_2*xjall(i, k0)
    struek1 = st(i, k01) + thomson_21*xjall(i, k01)

    strue(i) = qred*struek + qblue*struek1
    xjold(i) = qred*xjall(i, k0) + qblue*xjall(i, k01)

    Read (1, Rec=i)(blevel(j), j=1, nrec)
    Read (11, Rec=i)(alevel(j), j=1, nrec)
    Do ik = 1, nb
      xnl = blevel(nl(ik))
      If (xnl==0.) xnl = 1.D-15*alevel(nl(ik))
      xnu = blevel(nu(ik))
      If (xnu==0.) xnu = 1.D-15*alevel(nu(ik))
!     corrected
      xkl = c1*gflu(ik)*(xnl/gl(ik)-xnu/gu(ik))/clf(i)
      tcl = abs(xkl)*xlamb(ik)*tcl_fac_line(i) ! consistent inversion
!     treatment
!     Note: SRVMAX term included in tcl_fac_line multiplier
      fac2 = (1.+fic(i)*tcl)/(1.+tcl)
!     correction for optically thick clumping as in nlte.f90,
!     no back-correction of background
      fac2 = min(fac2, opa_eff_rat(i,k0), opa_eff_rat(i,k01))
      If (fac2<=0.0D0 .Or. fac2>1.0D0) Then
        Print *, fac2
        Print *, tcl, i, ik
        Stop ' OPA_EFF > <OPA>, ELSCAT'
      End If
      opal(i, ik) = xkl*srvmax*xlamb(ik)*fac2
      sline(i, ik) = c2*(xnue(ik)**3)/(xnl/xnu*gu(ik)/gl(ik)-1.D0)
      If (opal(i,ik)<0.) Then
        Print *, 'INVERSION OF COMPONENT', ik, ' AT L = ', i
      End If
    End Do
  End Do

  Close (1)
  Close (11)

! calculation of the p-grid

  nc = np - nd
  Do l = 1, nc
    p(l) = .99D0*dble(l-1)/dble(nc-1)
  End Do

  Do l = nc + 1, np
    p(l) = r(np+1-l)
  End Do

! calculation of the z-grid

  Do l = 1, np
    lmax = min0(np+1-l, nd)
    Do j = 1, lmax
      z2 = r(j)*r(j) - p(l)*p(l)
      If (l>nc .And. j==lmax) Then
        If (z2>1.D-5) Stop 'ERROR IN Z=0'
        z2 = 0.D0
      End If
      z(j, l) = sqrt(z2)
    End Do
  End Do

! set up of cmf freq. grid ----------------------------------

! ---  assumes that turbulence velocity is negligible for electrons

  vthel = 5.5056D5*sqrt(teff)
  Print *, 'vmax = ', vmax/1.D5, ' vth_el = ', vthel/1.D5

! from x = 4*vthel,-2*vinf-3*vthel
! min(lam) from xlamb(1), max(lam) from xlamb(nb)

  lammin = xlamb(1)/(4.*vthel/clight+1.)
  lammax = xlamb(nb)/((-3.-2.*vmax/vthel)*vthel/clight+1.)
! PRINT*,LAMMIN*1.D8,LAMMAX*1.D8,XLAMB*1.D8

  vref = max(vturb, 5.D5)
  vref = vref/3.D0 !                3 points per vth is sufficient
  Print *
  Print *, 'freq. separation inside doppler core = ', vref/1.D5, ' km/s'

! x for lammin/lammax (again: w.r.t. blue component)
  xmax = (xlamb(1)/lammin-1.D0)*clight/vmax
  xmin = (xlamb(1)/lammax-1.D0)*clight/vmax

! transition frequencies in X, w.r.t. blue component
  Do ik = 1, nb
    xc(ik) = (xlamb(1)/xlamb(ik)-1.D0)*clight/vmax
  End Do

  dx = 0.05 !                       5% is typical maximum separation in V (do
! not use larger DX values
! in non core region; otherwise => oscillations in the red
  dx = max(dx, 50.D5/vmax) !        minimum dx corresponding to  50 km/s
  dx = min(dx, 100.D5/vmax) !       maximum dx corresponding to 100 km/s
  dxcoarse = dx
  dxfine = vref/vmax

! account for fine-resolution around line-cores
  idv = 6*3 !                       6 Doppler width * 3*vref; might need to be
! changed to 10

all_comp: Do ik = 1, nb
    If (ik==1) Then
      xx(1) = xc(ik) !              line-center(IK=1)
      inext = 2
      Do i = inext, inext + idv - 1
        xx(i) = xc(ik) + (i+1-inext)*dxfine ! 6 vth towards blue
      End Do
      inext = inext + idv
    End If

    If (ik==nb) Then
      Do i = inext, inext + idv - 1
        xx(i) = xc(ik) - (i+1-inext)*dxfine ! 6 vth towards red
      End Do
      ilast = inext + idv - 1
      Exit all_comp
    End If

    If (xc(ik)-xc(ik+1)>(2*idv*dxfine)) Then
      Do i = inext, inext + idv - 1
        xx(i) = xc(ik) - (i+1-inext)*dxfine ! 6 vth towards red of IK
      End Do
      inext = inext + idv
      xx(inext) = xc(ik+1) !        line-center(IK+1)
      inext = inext + 1
      Do i = inext, inext + idv - 1
        xx(i) = xc(ik+1) + (i+1-inext)*dxfine ! 6 vth towards blue of IK+1
      End Do
      inext = inext + idv
    Else
      If (xc(ik)==xc(ik+1)) Cycle all_comp ! identical frequencies
      idx = (xc(ik)-xc(ik+1))/dxfine + 1
      dx = (xc(ik)-xc(ik+1))/idx
      Do i = inext, inext + idx - 1
        xx(i) = xc(ik) - (i+1-inext)*dx ! towards red, including XC(IK+1)
      End Do
      inext = inext + idx
    End If
  End Do all_comp

  ilast = ilast + 1
  xx(ilast) = xmax
  ilast = ilast + 1
  xx(ilast) = xmin

  xx = -xx
  Call sort(ilast, xx)
  xx = -xx


! check that transition frequencies are included
  Do ik = 1, nb
    Do i = 1, ilast
      If (abs(xx(i)-xc(ik))<1.D-15) Go To 100
    End Do
    Stop ' XC(IK) NOT FOUND IN XX!'
100 Continue
  End Do

! now fill the gaps with DXCOARSE

  inext = ilast + 1
  Do i = 1, ilast - 1
    xxdiff = xx(i) - xx(i+1)
    If (xxdiff>1.2*dxcoarse) Then
      idx = xxdiff/dxcoarse + 1
      dx = xxdiff/idx
      Do k = inext, inext + idx - 2 ! only inside interval
        xx(k) = xx(i) - (k+1-inext)*dx
      End Do
      inext = inext + idx - 1
    End If
  End Do
  ilast = inext - 1

  If (ilast>nfmax) Stop 'NFMAX=ID_NFESC TOO SMALL'

  nfcmf = ilast

  xx = -xx
  Call sort(nfcmf, xx)
  xx = -xx

! for tests
! DO I=1,NFCMF
! PRINT*,I,XX(I),XX(I)-XX(I+1)
! ENDDO

  Print *, 'number of cmf-frequency points = ', nfcmf

  Allocate (aick(nfcmf), dbdrk(nfcmf), eold(nfcmf), x(nfcmf), xi(nfcmf), &
    fi(nfcmf))

  Allocate (etac(nfcmf,nd), aj(nfcmf,nd), e(nfcmf,nd), xjnew(nfcmf,nd), &
    opal_tot(nfcmf,nd), etal_tot(nfcmf,nd))

! now xi = ln (nue)

! nue_0 (corresponding to x=0) taken from first component

  Do i = 1, nfcmf
    x(i) = xx(i)
    lam_cmf = xlamb(1)/(x(i)*vmax/clight+1.D0)
    xnu = clight/lam_cmf
    xi(i) = log(xnu)
  End Do

! end set up of cmf freq. grid ------------------------------


! blue wing boundary (remains fixed during subsequent iteration)

  xlambda = xlamb(1)*1.D8/(vmax*x(1)/clight+1.D0)
  Call diffus(xlambda, t, r, nd, aic, dbdr)
  corrfc = 1.D0

  Call cont1(nd, np, r, p, z, strue, opac, thomson, aic, dbdr, corrfc, &
    xjall(1,1), ubluwi, w0, w2)
! now xjall(i,1) is continuum J at x=0

  err = 0.D0

  Do l = 1, nd - 1
!   last point differs, due to irradiation from x(1) instead of x=0
!   print*,l,' ',xjall(l,1)/xjold(l)
!   write(*,fmt='(i2,2x,f10.5,3(2x,e11.5),2x,f10.5)') &
!   &       l,xjall(l,1)/xjold(l),xjold(l),strue(l),opac(l),thomson(l)
    err = max(err, abs(1.D0-xjall(l,1)/xjold(l)))
  End Do

  Print *
  Print *, ' MAX. INCONSISTENCY IN JNUE(CONT): ', err
  Print *

  If (err>0.01) Stop 'INCONSISTENCY IN JNUE(CONT) TOO LARGE!'

! for tests of stark-broadening
! PRINT*,'lambda, 20000/1D13 40000/1D14 20000/1D14 40000/1D14'
! DO K = 1,NFCMF
! XLAMBDA = XLAMB(2)*1.D8/ (VMAX*X(K)/CLIGHT+1.D0)
! PHIKL1 = PERFIL1(2,X(K),20000.d0,1.D13,VTURB,FLAG)
! PHIKL2 = PERFIL1(2,X(K),40000.d0,1.D13,VTURB,FLAG)
! PHIKL3 = PERFIL1(2,X(K),20000.d0,1.D14,VTURB,FLAG)
! PHIKL4 = PERFIL1(2,X(K),40000.d0,1.D14,VTURB,FLAG)
! WRITE(*,FMT='(5(E12.6,2X))'),XLAMBDA-XLAMB(2)*1.D8,LOG10(PHIKL1), &
! &       LOG10(PHIKL2),LOG10(PHIKL3),LOG10(PHIKL4)
! ENDDO
! STOP


  Do l = 1, nd
    xjnew(:, l) = xjall(l, 1)

    Do k = 1, nfcmf

!     -----perfil with respect to opal*sr instead of opal*srvmax*lambda

      phikl = perfil1(1, x(k), t(l), xne(l), vturb, flag)*vmax/xlamb(1)
      aux = opal(l, 1)*phikl
      opal_tot(k, l) = aux
      etal_tot(k, l) = sline(l, 1)*aux
      Do ik = 2, nb
        phikl = perfil1(ik, x(k)-xc(ik), t(l), xne(l), vturb, flag)*vmax/ &
          xlamb(ik)
        aux = opal(l, ik)*phikl
        opal_tot(k, l) = opal_tot(k, l) + aux
        etal_tot(k, l) = etal_tot(k, l) + sline(l, ik)*aux
      End Do
!     remember the old bug for the bluest frequencies: k1 instead of k
      xlambda = xlamb(1)*1.D8/(vmax*x(k)/clight+1.D0)
      Call diffus(xlambda, t, r, nd, aick(k), dbdrk(k))
    End Do

  End Do


! DELTAX = VREF/VMAX  !old version; now set in subr. RAY

! total opacity (does not change during iteration)
  Do k = 1, nfcmf
    opal_tot(k, :) = opal_tot(k, :) + opac
  End Do
! now, opal_tot is total opacity (lines + cont)

  Do k = 1, nfcmf
    dum = opal_tot(k, :)
    Where (dum<=0.)
!     check for inversion in total opacity;
!     intermed. hack: use only cont. quantities
      opal_tot(k, :) = opac
      etal_tot(k, :) = 0.
    End Where
  End Do

! begin of iteration

itloop: Do iit = 1, 20

!   =======================================
!   cmf solution
!   =======================================


    rel = 0.D0
    aj = 0.D0

    If (iit==1) Then
!     for inversion, etal_tot has been set to zero
      Do k = 1, nfcmf
        etac(k, :) = (strue+thomson*xjnew(k,:))*opac + etal_tot(k, :)
      End Do

    Else
      Do k = 1, nfcmf
        etac(k, :) = (strue+thomson*e(k,:))*opac + etal_tot(k, :)
      End Do
    End If
!   now, etac is total emissivity (lines + cont); changes during iteration

!   for tests
!   do k=1,nfcmf
!   do l=1,nd
!   WRITE (*,FMT='(2(i3,2x),3(E12.6,2x))') k,l,x(k),opal_tot(k,l),etac(k,l)
!   enddo
!   enddo
!   stop

jploop1: Do jp = 1, np - 1
      lmax = min0(np+1-jp, nd)
      lz = lmax - 1

!     bluewing boundary condition

      ucmf(1) = ubluwi(1, jp)

      Do l = 1, lz
        ucmf(l+1) = ubluwi(l+1, jp)
        aux = .5D0*(opac(l)+opac(l+1))
        vcmf(l) = (ucmf(l+1)-ucmf(l))/aux/(z(l,jp)-z(l+1,jp))
      End Do

      Call ray(jp, z, r, v, dvdr, nd, np, ucmf, vcmf, opal_tot, etac, nfcmf, &
        aick, dbdrk, aj, w0(1,jp), x)

    End Do jploop1

!   DO K = 1,NFCMF
!   ERR = 0.D0
!   DO L = 1,ND
!   ERR = MAX(ERR,ABS(1.D0-AJ(K,L)/XJOLD(L)))
!   END DO
!   print*
!   print*,' max. inconsistency in jnue: ',err,' at x(',k,') = ',x(k)
!   END DO

    xjnew = aj
!   hack to account for edge effects (inconsistency
!   cont-transfer/cmf-transfer)
    xjnew(1, :) = xjnew(2, :)

!   normalize J(cmf) to (more accurate) J(cont)
!   in the ideal case, xjnew (cmf) should be equal to j(cont) at blue and red
!   edges
!   (assuming that freq. dependent lower boundary makes no difference expect
!   for lowermost point)
    Do l = 1, nd - 1
      fac1 = xjnew(1, l)/xjold(l)
      fac2 = xjnew(nfcmf, l)/xjold(l)
      If (l==1 .And. (abs(1.-fac1)>0.05 .Or. abs(1.-fac2)>0.05)) Then
        Print *, fac1, fac2
        Print *, &
          'WARNING !!! LARGE INCONSISTENCY BETWEEN J(CMF) AND J(CONT)!!!'
      End If
      If (abs(1.-fac1)>0.1) Stop &
        'LARGE INCONSISTENCY BETWEEN J(CMF) AND J(CONT) AT BLUE EDGE'
      If (abs(1.-fac2)>0.1) Stop &
        'LARGE INCONSISTENCY BETWEEN J(CMF) AND J(CONT) AT RED EDGE'
      If (abs(1.-fac1/fac2)>0.02) Then
        Print *, fac1, fac2
        Stop 'INCONSISTENCY AT BLUE AND RED EDGE FREQUENCY'
      End If
      fac1 = 0.5*(fac1+fac2)
      fac1 = 1./fac1
      xjnew(:, l) = xjnew(:, l)*fac1
    End Do

!   now rybicki and hummer algorithm at cmf frequencies
!   (dipole approx)

!   goto 510
lloop: Do l = 1, nd

      If (iit>=2) eold = e(:, l)

      betat = 1.84D-5*sqrt(t(l))
      const = (betat/b1)**2
!     AJ(:,I),I=1,3 as dummy for TA,TB,TC to save storage
      Call scatsol(xi, xjnew(1,l), fi, const, aj(:,1), aj(:,2), aj(:,3), &
        nfcmf)

      Do k = 1, nfcmf
        e(k, l) = a1*fi(k)
      End Do

      const = (betat/b2)**2
      Call scatsol(xi, xjnew(1,l), fi, const, aj(:,1), aj(:,2), aj(:,3), &
        nfcmf)

      Do k = 1, nfcmf
        e(k, l) = e(k, l) + a2*fi(k)
        If (iit>=2) rel = max(rel, abs(1.D0-e(k,l)/eold(k)))
      End Do

    End Do lloop

    If (iit>=2) Print *, ' ITERATION NO.', iit, ' MAX DEV. IN E(K,L) = ', rel
    If (iit>=2 .And. rel<=1.D-3) Then
      Print *
      Print *, ' CONVERGENCY ACHIEVED'
      Print *
      Go To 110
    End If

  End Do itloop

  Print *
  Print *, ' ITERATION NOT CONVERGED IN 20 ITERATIONS'
  Print *

110 Continue


! change of first and last freq. to be consistent with xgrid
! assuming that pure continuum is reached

  If (x(1)-x(2)<1.) x(1) = x(2) + 1.01D0

  If (x(nfcmf-1)-x(nfcmf)<1.) x(nfcmf) = x(nfcmf-1) - 1.01D0

! !      open(1,file=TRIM(FILE)//'/exi_'//line,
! !     *status='unknown',form='formatted')

  Open (17, Status='SCRATCH', Form='FORMATTED')
  Write (17, Fmt=*) nd, nfcmf
  Do k = 1, nfcmf
    Write (17, Fmt=*) k, ' ', x(k), ' ', xi(k)
  End Do

  Do k = 1, nfcmf
    Do l = 1, nd
      scold = strue(l) + thomson(l)*xjold(l)
      scnew = strue(l) + thomson(l)*e(k, l)
      Write (17, Fmt=*) k, ' ', l, ' ', e(k, l), ' ', scnew, ' ', scnew/scold
!     WRITE (17,FMT=*) K,' ',L,' ',E(K,L),' ',SCOLD,' ', &
!     SCNEW/SCOLD
!     if we replace here scnew by scold, we can test whether the result is
!     consistent with a calculation without elscat

!     WRITE (*,FMT='(E12.6,1x,I3,2(1x,E12.6))') X(K),L,E(K,L),SCNEW/SCOLD
!     WRITE (*,FMT='(E12.6,1x,I3,2(1x,E12.6))')
!     X(K),L,E(K,L),XJNEW(K,L)/XJOLD(L)
    End Do
  End Do

! !      close(1)
! stop

  Deallocate (aick, dbdrk, eold, x, xi, fi)
  Deallocate (etac, aj, e, xjnew, opal_tot, etal_tot)
  Deallocate (opal, sline, xlamb, xc)

  Return
End Subroutine

!-----------------------------------------------------------------------

Subroutine cont1(nd, np, r, p, z, st, opa, thomson, xic1, xic2, corr, xj, u, &
  vp, vp2)

  Use :: nlte_type
  Use :: nlte_dim
  Implicit None

! .......................................................................

! adaptation from nonplusultra. only j is calculated

! attention!!! subroutine cont1 changed for fast and consistent
! continuum transfer in program nonplusultra. so far, only case
! i) can be treated.

! calculates the solution of the following equation applying the
! rybicki algorithm (including the solution of moments' equation,
! if necessary , for case i) )

! xj = lambda (st+thomson*xj) ,

! or, explicitly written

! 1)     j = lamdda (s)
! 2) dj/dt = lambda (ds/dt)

! where s is the usual nlte-source function including electron
! scattering

! upper boundary condition : valid to 3rd order including self-
! consistent calculated i-minus  (dj/dtau neglected)

! lower boundary condition: valid to 2nd order
! diffusion approximation possible : xic1 = bnue
! xic2 = dbnue/dr
! rescaling of lower input flux possible : corr = h-true/h-actual


! case 1) input:

! st = (eta - opa-thomson*j)/opa
! opa = total opacity * stellar radius
! thomson = opa-thomson/opa
! xic1 = i-plus at lower boundary (diff.approx. = bnue(t))
! xic2 = 0. (diff. approx. = dbnue/dr)

! output : xj...mean intensity
! xh...eddington flux (optflux=.true.)


! case 2) input:

! st = d/dt((eta - opa-thomson*j)/opa) + d/dt(thomson) *j
! opa , thomson as  above
! xic1 = d/dt(i-plus)  (diff. approx. = d/dt(bnue(t))
! xic2 = 0. (diff. approx. = d2/drdt (bnue(t)) )

! output : xj = dj/dt


! .......................................................................


! .. parameters ..
  Integer (i4b), Parameter :: nd1 = id_ndept, np1 = id_npoin
  Integer (i4b), Parameter :: imax = 10
! ..
! .. scalar arguments ..
  Real (dp) :: corr, xic1, xic2
  Integer (i4b) :: nd, np
! ..
! .. array arguments ..
  Real (dp) :: opa(nd1), p(np1), r(nd1), st(nd1), thomson(nd1), u(nd1, np1-1), &
    vp(nd1, np1-1), vp2(nd1, np1-1), xj(nd1), z(nd1, np1)
! ..
! .. local scalars ..
  Real (dp) :: e0, e1, edmue, epsmax, heddi, heddo, hin, hout, pi2, ux
  Integer (i4b) :: iii, jp, l, lmax
  Logical :: conv, core
! ..
! .. local arrays ..
  Real (dp) :: akp(nd1), aksum(nd1), f(nd1+3), qq(nd1), ta(nd1), taum(np1), &
    tb(nd1), tb1(nd1), tc(nd1), tc1(nd1), tp(nd1, nd1), tsum(nd1, nd1), &
    up(nd1), wp(nd1+1, np1-1), xjold(nd1), xk(nd1)
! ..
! .. external subroutines ..
  External :: dev, inv, invtri, invtri3, madd, mdmv, mdv, momcont, mvmd, mvv, &
    setup1, vadd, vsub, weigh11
! ..
! .. intrinsic functions ..
! ..

  pi2 = acos(0.D0)
  conv = .False.

  Call weigh11(nd, np, r, z, p, vp, vp2, wp)

  aksum = .0D0
  tsum = .0D0

  Do jp = 1, np - 1

    If (jp==1) Then
      taum(1) = opa(1)*r(1)*(1.D0/3.D0+2.D0*thomson(1)/3.D0)
    Else
      taum(jp) = opa(1)*r(1)*r(1)*(thomson(1)/p(jp)*(pi2-atan(z(1,jp)/p(jp)))+ &
        (1.D0-thomson(1))*r(1)*r(1)/2.D0/p(jp)**3.D0*(pi2-atan(z(1,jp)/ &
        p(jp))-z(1,jp)*p(jp)/r(1)/r(1)))
    End If

    If (taum(jp)<0.D0) Stop 'TAUM NEGATIVE!'

    lmax = min0(np+1-jp, nd)
    core = (lmax==nd)

!   calculation of tp,akp

    Call setup1(lmax, core, z(1,jp), st, opa, thomson, xic1, xic2, corr, up, &
      akp, ta, tb, tc, taum(jp))

!   rybicki algorithm to obtain ajc

    Call invtri3(ta, tb, tc, tp, lmax, nd)
    Call mdmv(tp, vp(1,jp), lmax, nd)
    Call mvv(qq, tp, akp, lmax, lmax, nd)
    Call vadd(aksum, qq, lmax)

!   tp=:(vp*tp-1)*u

    Call mvmd(tp, up, lmax, nd)
    Call madd(tsum, tp, lmax, nd)

  End Do

! addition of unity-matrix to tsum

  Do l = 1, nd
    tsum(l, l) = tsum(l, l) + 1.D0
  End Do

  Call inv(nd, nd, tsum)
  Call mvv(xj, tsum, aksum, nd, nd, nd)

  qq = xj
  xjold = xj

! backsubstitution to obtain u

  iii = 0

100 Continue

  iii = iii + 1

  Do jp = 1, np - 1

    lmax = min0(np+1-jp, nd)
    core = (lmax==nd)

    Call setup1(lmax, core, z(1,jp), st, opa, thomson, xic1, xic2, corr, up, &
      akp, ta, tb, tc, taum(jp))

    Do l = 1, lmax
      ta(l) = -ta(l)
      tc(l) = -tc(l)
    End Do

!   recalculation of total source-function

    Call mdv(qq, up, lmax)
    Call vsub(akp, up, lmax)

    Call invtri(ta, tb, tc, akp, lmax)

    Do l = 1, lmax
      u(l, jp) = akp(l)
    End Do

  End Do

! recalculation of j and calculation of k

  xj = 0.D0
  xk = 0.D0

  Do jp = 1, np - 1
    lmax = min0(np+1-jp, nd)
    Do l = 1, lmax
      ux = u(l, jp)
      xj(l) = xj(l) + vp(l, jp)*ux
      xk(l) = xk(l) + vp2(l, jp)*ux
    End Do
  End Do

  Do l = 1, nd
    f(l) = xk(l)/xj(l)
  End Do

! calculation of inner and outer eddington factors respective to h

! ---- this is the old form for boundary conditions (i.e., achim's form
! ---- instead of gudrun's). i prefer it because the first one takes into
! ---- account the new radiation field in the calculation of h-nue as the
! ---- outer boundary condition.

  hin = 0.D0
  hout = 0.D0
  edmue = 0.D0

  Do jp = 1, np - 1
    lmax = min0(np+1-jp, nd)

!   outer boundary

    If (taum(jp)<100.D0) Then
      e0 = exp(-taum(jp))
      e1 = 1.D0 - e0
    Else
      e1 = 1.D0
    End If

!   dxi=u(1,jp)-e1*(st(1)+thomson(1)*xj(1))
!   ---- this line has been masked (21-10-92). it corresponds to case
!   ---- i_plus<i_minus (at r_max), and now there is no difference
!   between c---- both cases.(also, now is equivalent to gudrun's
!   formulation but c---- the comment stated above).
!   if(dxi.lt.0.d0) e1=u(1,jp)/(st(1)+thomson(1)*xj(1))

    hout = hout + u(1, jp)*wp(1, jp)
    edmue = edmue + e1*wp(1, jp)

!   ---- with the last modification, edmue is always i_minus(1)/s_total(1)
!   ---- integrated over mue

!   inner boundary

    If (lmax==nd) Then
      hin = hin + u(nd, jp)*wp(nd+1, jp)
    End If

  End Do

  heddo = hout/xj(1)
  heddi = hin/xj(nd)
  f(nd+1) = heddo
  f(nd+2) = heddi
  f(nd+3) = edmue

  Call momcont(nd, r, opa, thomson, st, xic1, xic2, corr, qq, ta, tb, tc, tb1, &
    tc1, akp, xj, f)

  If (conv .Or. iii==imax) Then
!   PRINT *,' ACHIEVED ACCURACY IN CONT. TRANSPORT  = ',EPSMAX
!   PRINT *,' IN ',III,' ITERATIONS'
    Return
  End If

! -----------------------------------------------------------------------

  Call dev(r, xj, xjold, nd, epsmax, .False.)

  If (epsmax<1.D-3 .Or. iii==imax) conv = .True.

  qq = xj
  xjold = xj

  Go To 100

End Subroutine

!-----------------------------------------------------------------------

Subroutine dev(r, a1, a, nd, epsmax, opt)

  Use :: nlte_type
  Use :: nlte_dim
  Implicit None


! .. scalar arguments ..
  Real (dp) :: epsmax
  Integer (i4b) :: nd
  Logical :: opt
! ..
! .. array arguments ..
  Real (dp) :: a(nd), a1(nd), r(nd)
! ..
! .. local scalars ..
  Real (dp) :: e1, e2, e3, e4, e5, eps1, eps2, eps3, eps4, eps5
  Integer (i4b) :: l
! ..
! .. intrinsic functions ..
! ..

  eps1 = 0.D0
  eps2 = 0.D0
  eps3 = 0.D0
  eps4 = 0.D0
  eps5 = 0.D0

  Do l = 1, nd
    If (r(l)>=50.D0) Then
      e1 = abs(1.D0-a1(l)/a(l))
      If (e1>eps1) eps1 = e1

    Else If (r(l)<50.D0 .And. r(l)>=10.D0) Then
      e2 = abs(1.D0-a1(l)/a(l))
      If (e2>eps2) eps2 = e2

    Else If (r(l)<10.D0 .And. r(l)>=5.D0) Then
      e3 = abs(1.D0-a1(l)/a(l))
      If (e3>eps3) eps3 = e3

    Else If (r(l)<5.D0 .And. r(l)>=2.D0) Then
      e4 = abs(1.D0-a1(l)/a(l))
      If (e4>eps4) eps4 = e4

    Else If (r(l)<2.D0) Then
      e5 = abs(1.D0-a1(l)/a(l))
      If (e5>eps5) eps5 = e5

    End If
  End Do

  epsmax = max(eps1, eps2, eps3, eps4, eps5)

  If (opt) Then
    Print *
    Print *, ' MAXIMUM DEVIATION ', epsmax
    Print *
    Print *, '  1. < R <   2.  : ', eps5
    Print *, '  2. < R <   5.  : ', eps4
    Print *, '  5. < R <  10.  : ', eps3
    Print *, ' 10. < R <  50.  : ', eps2
    Print *, ' 50. < R < 100.  : ', eps1
  End If

  Return
End Subroutine

!-----------------------------------------------------------------------

Subroutine momcont(nd, r, opa, thomson, st, xic1, xic2, corr, q, ta, tb, tc, &
  tb1, tc1, akp, xj, f)

  Use :: nlte_type
  Use :: nlte_dim
  Implicit None


! -----solves moments equation for continuum

! .. parameters ..
  Integer (i4b), Parameter :: nd1 = id_ndept
! ..
! .. scalar arguments ..
  Real (dp) :: corr, xic1, xic2
  Integer (i4b) :: nd
! ..
! .. array arguments ..
  Real (dp) :: akp(nd1), f(nd1+3), opa(nd1), q(nd1), r(nd1), st(nd1), ta(nd1), &
    tb(nd1), tb1(nd1), tc(nd1), tc1(nd1), thomson(nd1), xj(nd1)
! ..
! .. local scalars ..
  Real (dp) :: dt0, dtm, dtp, edmue, fl, flp, hi, ho, rl, rl2, rlp, rrq
  Integer (i4b) :: l
! ..
! .. external subroutines ..
  External :: invtri
! ..
! .. intrinsic functions ..
! ..

  ho = f(nd+1)
  hi = f(nd+2)
  edmue = f(nd+3)

  q(nd) = 1.D0
  rrq = 1.D0
  fl = 3.D0 - 1.D0/f(nd)

  Do l = nd - 1, 1, -1
    rl = r(l)
    rlp = r(l+1)
    flp = fl
    fl = 3.D0 - 1.D0/f(l)
    rrq = rrq*exp(fl-flp)*(rl/rlp)**((flp*rl-fl*rlp)/(rl-rlp))
    q(l) = rrq/rl/rl
  End Do

! feautrier scheme to solve moments equation :

! (-ta,tb,-tc) * rr * xj = akp

! outer boundary condition

  dtp = 2.D0/((q(1)*opa(1)+q(2)*opa(2))*(r(1)-r(2)))
  tb(1) = (f(1)*q(1)*dtp+ho-edmue*thomson(1))
  tb1(1) = tb(1) + edmue*thomson(1)
  tc(1) = f(2)*q(2)*dtp
  tc1(1) = tc(1)
  akp(1) = r(1)*r(1)*st(1)*edmue

! non boundary points

  Do l = 2, nd - 1
    dtm = dtp
    dtp = 2.D0/((q(l)*opa(l)+q(l+1)*opa(l+1))*(r(l)-r(l+1)))
    dt0 = 2.D0/(1.D0/dtp+1.D0/dtm)
    ta(l) = f(l-1)*q(l-1)*dt0*dtm
    tc(l) = f(l+1)*q(l+1)*dt0*dtp
    tc1(l) = tc(l)
    tb(l) = f(l)*q(l)*dt0*(dtm+dtp) + (1.D0-thomson(l))/q(l)
    tb1(l) = tb(l) + thomson(l)/q(l)
    akp(l) = r(l)*r(l)*st(l)/q(l)
  End Do

! inner boundary condition

  l = nd
  ta(l) = f(l-1)*q(l-1)*dtp
  tb(l) = f(l)*q(l)*dtp + hi
  tb1(l) = tb(l)
  akp(l) = .5D0*xic1 + xic2*corr/(3.D0*opa(l))
  tc1(nd) = akp(l)

  Call invtri(ta, tb, tc, akp, nd)

  Do l = 1, nd
    rl2 = r(l)*r(l)
    xj(l) = akp(l)/rl2
  End Do

  Return
End Subroutine

!-----------------------------------------------------------------------

Subroutine setup1(lmax, core, z, st, opa, thomson, xic1, xic2, corr, up, akp, &
  ta, tb, tc, taum)

  Use :: nlte_type
  Use :: nlte_dim
  Implicit None

! sets up matrix elements for subroutine cont1

! .. parameters ..
  Integer (i4b), Parameter :: nd1 = id_ndept
! ..
! .. scalar arguments ..
  Real (dp) :: corr, taum, xic1, xic2
  Integer (i4b) :: lmax
  Logical :: core
! ..
! .. array arguments ..
  Real (dp) :: akp(nd1), opa(nd1), st(nd1), ta(nd1), tb(nd1), tc(nd1), &
    thomson(nd1), up(nd1), z(nd1)
! ..
! .. local scalars ..
  Real (dp) :: ak, akdz, bb, cc, dt0, dtm, dtp, dz, e0, e1
  Integer (i4b) :: l, lz
! ..
! .. intrinsic functions ..
! ..

  lz = lmax - 1

! outer boundary, 3rd order, corrected for i-minus

  ak = .5D0*(opa(1)+opa(2))
  dz = z(1) - z(2)
  akdz = .5D0*ak*dz
  dtp = 1.D0/ak/dz
  bb = 1.D0/dtp/3.D0
  cc = .5D0*bb

  If (taum<100.D0) Then
    e0 = exp(-taum)
    e1 = 1.D0 - e0
  Else
    e1 = 1.D0
  End If

  If (akdz<.5D0) Then
    up(1) = thomson(1)*(bb+e1) + thomson(2)*cc
    akp(1) = -st(1)*(bb+e1) - st(2)*cc
    tb(1) = -1.D0 - dtp - bb
    tc(1) = dtp - cc
  Else
    up(1) = thomson(1)*e1
    akp(1) = -st(1)*e1
    tb(1) = -1.D0 - dtp
    tc(1) = dtp
  End If

! non boundary points

  Do l = 2, lz
    dtm = dtp
    ak = .5D0*(opa(l)+opa(l+1))
    dz = z(l) - z(l+1)
    dtp = 1.D0/ak/dz
    dt0 = 2.D0/(1.D0/dtm+1.D0/dtp)
    up(l) = thomson(l)
    akp(l) = -st(l)
    ta(l) = dt0*dtm
    tc(l) = dt0*dtp
    tb(l) = -dt0*(dtm+dtp) - 1.D0
  End Do

  l = lmax

! inner boundary,2nd order

  If (core) Then
    up(l) = 0.D0
    ta(l) = dtp
    akp(l) = -xic1 - z(l)*xic2/opa(l)*corr
    tb(l) = -dtp - 1.D0
  Else
    akp(l) = -st(l)
    up(l) = thomson(l)
    ta(l) = 2.D0*dtp*dtp
    tb(l) = -2.D0*dtp*dtp - 1.D0
  End If

  Return
End Subroutine

!-----------------------------------------------------------------------

Subroutine weigh11(nd, np, r, z, p, w0, w2, w1)

  Use :: nlte_type
  Use :: nlte_dim
  Implicit None

! ***  calculation of angular integration weights
! ***  the weights w1  are calculated for intermesh points
! ***  (in this context only for boundary conditions)


! .. parameters ..
  Integer (i4b), Parameter :: nd1 = id_ndept, np1 = id_npoin
! ..
! .. scalar arguments ..
  Integer (i4b) :: nd, np
! ..
! .. array arguments ..
  Real (dp) :: p(np1), r(nd1), w0(nd1, np1-1), w1(nd1+1, np1-1), &
    w2(nd1, np1-1), z(nd1, np1)
! ..
! .. local scalars ..
  Real (dp) :: a, aa, b, bb, c, cc, rl, rl2, rrr12, rz, rzq6, w1lz, ww1
  Integer (i4b) :: jp, l, lmax, lz
! ..
! .. intrinsic functions ..
! ..

jploop: Do jp = 1, np - 1

    lmax = min0(np+1-jp, nd)
    lz = lmax - 1
!   ***
!   ***  0. and 2. moment. the integration is performed in the z variable
    Do l = 1, lmax
      rl = r(l)
      rl2 = rl + rl
      rrr12 = rl*rl*rl*12.D0
!     ***  first step if jp=1
      If (jp==1) Then
        b = z(l, 1)
        a = z(l, 2)
        w0(l, jp) = (b-a)/rl2
        aa = a*a
        bb = b*b
        w2(l, jp) = (b*(3.D0*bb-aa)-a*(bb+aa))/rrr12
      Else
        If (l/=lmax .Or. jp<=(np-nd)) Then
!         ***  intermediate step
          a = z(l, jp+1)
          b = z(l, jp)
          c = z(l, jp-1)
          w0(l, jp) = (c-a)/rl2
          aa = a*a
          bb = b*b
          cc = c*c
          w2(l, jp) = (b*(cc-aa)+c*(cc+bb)-a*(bb+aa))/rrr12
        Else
!         ***  last step, implying z(l,jmax)=0
          b = z(l, jp-1)
          w0(l, jp) = b/rl2
          w2(l, jp) = b*b*b/rrr12
        End If
      End If
    End Do
!   ****
!   **** 1.moment.
!   ****  first step
!   ****
    If (jp==1) Then
      ww1 = p(2)*p(2)
    Else
!     ***  intermediate steps
      a = p(jp-1)
      b = p(jp)
      c = p(jp+1)
      ww1 = (a+b+c)*(c-a)
!     ***  for the last interval (l=lz to the p axis), the next point is an
!     ***  intermesh-point in p
      c = .5D0*(b+c)
      w1lz = (a+b+c)*(c-a)
    End If
!   ***  no weight for the z=0 point is calculated, as v=0 there for symmetry.
!   ***  loop over depth index l

    Do l = 1, lz
      rz = .5D0*(r(l)+r(l+1))
      rzq6 = rz*rz*6.D0
      If (l/=lz .Or. jp<=(np-nd)) Then
        w1(l+1, jp) = ww1/rzq6
      Else
        w1(l+1, jp) = w1lz/rzq6
      End If
    End Do
!   ***  special weights at the outer boundary for h
    rl = r(1)
    w1(1, jp) = ww1/rl/rl/6.D0
!   ***  special weights at the inner boundary
    If (lmax<nd) Cycle jploop

    If (jp>(np-nd)) Then
!     ******  core tangent ray, last step of integration
      ww1 = (b-a)*(2.D0*b+a)
    End If

    w1(nd+1, jp) = ww1/6.D0

  End Do jploop

  Return
End Subroutine

!-----------------------------------------------------------------------

Subroutine diffus(xlambda, t, r, nd, aic, dbdr)

  Use :: nlte_type
  Use :: nlte_dim
  Implicit None

! gives the planck function aic and its radius derivative dbdr
! at the inner boundary from the given temperature-stratification.

! .. parameters ..
  Integer (i4b), Parameter :: nd1 = id_ndept
! ..
! .. scalar arguments ..
  Real (dp) :: aic, dbdr, xlambda
  Integer (i4b) :: nd
! ..
! .. array arguments ..
  Real (dp) :: r(nd1), t(nd1)
! ..
! .. local scalars ..
  Real (dp) :: dtdr
! ..
! .. external functions ..
  Real (dp) :: bnue, dbdt
  External :: bnue, dbdt
! ..

  aic = bnue(xlambda, t(nd))
  dtdr = (t(nd)-t(nd-1))/(r(nd-1)-r(nd))
  dbdr = dbdt(xlambda, t(nd))*dtdr

  Return
End Subroutine

!-----------------------------------------------------------------------

Function bnue(xlambda, t)

  Use :: nlte_type
  Use :: nlte_dim
  Implicit None


! planck function,lambda in angstroem,t in kelvin
! bnue in erg per (cm**2 * sec * hertz )

! constanten: c1=h*c/k,c2=2*h*c


! .. parameters ..
  Real (dp), Parameter :: c1 = 1.4388354967334D8, c2 = 3.972970127D8
! ..
! .. scalar arguments ..
  Real (dp) :: t, xlambda, bnue
! ..
! .. intrinsic functons ..
! ..

  bnue = c2/(exp(c1/xlambda/t)-1.D0)/xlambda/xlambda/xlambda

  Return
End Function

!-----------------------------------------------------------------------

Function dbdt(x, t)

  Use :: nlte_type
  Use :: nlte_dim
  Implicit None


! .. parameters ..
  Real (dp), Parameter :: c1 = 1.4388354967334D8
! ..
! .. scalar arguments ..
  Real (dp) :: t, x, dbdt
! ..
! .. external functions ..
  Real (dp) :: bnue
  External :: bnue
! ..
! .. intrinsic functions ..
! ..

  dbdt = bnue(x, t)*c1/x/t/t/(1.D0-exp(-c1/x/t))

  Return
End Function

!-----------------------------------------------------------------------

Subroutine ray(jp, z, r, velo, gradv, nd, np, u, v, opal, etal, nf, aic, dbdr, &
  aj, w0, x)

  Use :: nlte_type
  Use :: nlte_dim
  Implicit None


! line radiation transfer in the comoving frame from a
! given source function for a given impact parameter jp.
! the integration is carried out in space to
! yield the mean intensity

! j-nue

! note, that here the cont. emissivity is varying


! .. parameters ..
  Integer (i4b), Parameter :: nd1 = id_ndept
! ..
! .. scalar arguments ..
  Integer (i4b) :: jp, nd, nf, np
! ..
! .. array arguments ..
  Real (dp) :: aic(nf), aj(nf, nd), dbdr(nf), opal(nf, nd), etal(nf, nd), &
    u(nd), v(nd-1), z(nd, np), r(nd), velo(nd), gradv(nd), w0(nd), x(nf)
! ..
! .. local scalars ..
  Real (dp) :: rl, xmue, xmue2, deltax
  Integer (i4b) :: k, l, lmax, lz
! ..
! .. local arrays ..
  Real (dp) :: ga(nd1), h(nd1), pp(nd1), pp1(nd1), qq(nd1), s(nd1), ta(nd1), &
    tb(nd1), tc(nd1), ub(nd1), va(nd1), vb(nd1)
! ..
! .. external subroutines ..
  External :: cmfset, gmalu, invtri, mdv, vadd, vmalv
! ..
! .. intrinsic functions ..
! ..

  lmax = min0(np+1-jp, nd)
  lz = lmax - 1

! pp(l)=velocity gradient,projected on the present ray,/deltax

  Do l = 1, lmax
    rl = r(l)
    xmue = z(l, jp)/rl
    xmue2 = xmue*xmue
    pp1(l) = (xmue2*gradv(l)+(1.D0-xmue2)*velo(l)/rl) ! here, without DELTAX
  End Do

! loop for all original frequency points

kloop: Do k = 1, nf
    If (k==1) Go To 100

    deltax = x(k-1) - x(k)
    pp = pp1/deltax
!   THIS IS THE MOST IMPORTANT NEW FEATURE.
!   For K=1, BW values are taken (as usual)
!   But for K=2, we set PP articifically to zero, and decouple K=2 from K=1.
!   In this way, a consistent BW boundary condition is calculated, which
!   does not spoil the overall solution. Otherwise, we always encounter
!   small (0.005) inconsistencies in the pseudo-continuum, because the
!   slight differences between BW (continuum, 2nd order, Rybicki)
!   and CMF solution (1st order, formal) are not damped out if no strong line
!   or optically thick continuum is present.
!   These differences are set to zero by the new method.
!   Note that this trick requires K=2 to be pure continuum!
    If (k==2) pp = 0.

    Call cmfset(z(1,jp), nd, lmax, ta, tb, tc, ub, va, vb, ga, h, s, &
      opal(k,:), etal(k,:), pp, dbdr(k), aic(k))

    Call vmalv(va, vb, v, qq, lmax)
    Call vadd(qq, s, lmax)

    u(1:lmax) = ub(1:lmax)*u(1:lmax)

    Call vadd(u, qq, lmax)
    Call invtri(ta, tb, tc, u, lmax)

!   tested; improves results
    Do l = 1, lmax
      If (u(l)<0.) u(l) = 0.D0
    End Do

!   now u is the new field at index k

    Call mdv(h, v, lz)
    Call gmalu(ga, u, lmax)
    Call vadd(v, ga, lz)

!   now v is the new field at index k

!   adding the new u to 0th moments aj

100 Continue


    Do l = 1, lmax
      aj(k, l) = aj(k, l) + w0(l)*u(l)
!     AJ(K,L) = AJ(K,L) + W0(L)*AMAX1(U(L),0.)
!     print*,k,l,u(l),v(l),aj(k,l)
    End Do

  End Do kloop

  Return
End Subroutine

!-----------------------------------------------------------------------

Subroutine cmfset(z, nd, lmax, ta, tb, tc, ub, va, vb, ga, h, s, opal, etal, &
  pp, dbdr, aic)

  Use :: nlte_type
  Use :: nlte_dim
  Implicit None



! sets up the array elements for the cmf formalism


! ..
! .. scalar arguments ..
  Integer (i4b) :: lmax, nd
  Real (dp) :: aic, dbdr

! ..
! .. array arguments ..
  Real (dp) :: ga(nd), h(nd), s(nd), opal(nd), etal(nd), pp(nd), ta(nd), &
    tb(nd), tc(nd), ub(nd), va(nd), vb(nd), z(nd)
! ..
! .. local scalars ..
  Real (dp) :: ak, az, daz, dazm, dbz, dbzm, dt, dtz, dtzm, dx, dxz, dxz1, &
    dxzm, dxzm1, ek, tau, tauz
  Integer (i4b) :: l, lz
! ..
! .. intrinsic functions ..
! ..

  lz = lmax - 1

! ***  outer boundary condition  -  2nd order

  ak = opal(1)
  az = 0.5D0*(ak+opal(2))
  tauz = z(1) - z(2)
  dx = pp(1)/ak
  s(1) = 0.D0
  tc(1) = 1.D0/ak/tauz
  tb(1) = tc(1) + dx + 1.D0
  ub(1) = dx
  vb(1) = 0.0D0

! ***  for g and h, the matrix element s are not different from inner
! points
  dtzm = 1.D0/az/tauz
  dxzm = (pp(1)+pp(2))/az/2.D0
  dxzm1 = 1.D0/(1.D0+dxzm)
  dazm = dtzm*dxzm1
  dbzm = dxzm*dxzm1
  ga(1) = -dazm
  h(1) = dbzm

! ***  non-boundary points

  Do l = 2, lz
    ak = opal(l)
    ek = etal(l)
    s(l) = ek/ak
    az = 0.5D0*(ak+opal(l+1))
    tau = 0.5D0*(z(l-1)-z(l+1))
    tauz = z(l) - z(l+1)
    dt = 1.D0/ak/tau
    dtz = 1.D0/az/tauz
    dx = pp(l)/ak
    dxz = (pp(l)+pp(l+1))/az/2.D0
    dxz1 = 1.D0/(1.D0+dxz)
    daz = dtz*dxz1
    dbz = dxz*dxz1
    ta(l) = dt*dazm
    tc(l) = dt*daz
    tb(l) = ta(l) + tc(l) + dx + 1.D0
    ub(l) = dx
    va(l) = -dt*dbzm
    vb(l) = dt*dbz
    ga(l) = -daz
    h(l) = dbz
    dazm = daz
    dbzm = dbz

  End Do

  l = lmax

  If (lmax<nd) Go To 100

! ***  inner boundary condition (core rays)  -  only to first order
! ***   diffusion-approximation

  ak = opal(l)
  s(l) = aic + z(l)/ak*dbdr
  tauz = z(l-1) - z(l)
  dx = pp(l)/ak
  dt = 1.D0/tauz/ak
  ta(l) = dt
  tb(l) = dt + dx + 1.D0
  ub(l) = dx
  va(l) = 0.0D0
  vb(l) = 0.0D0

  Return

! ***  inner boundary condition (non-core rays)  -  second order

100 Continue

  ak = opal(l)
  ek = etal(l)
  s(l) = ek/ak
  tauz = z(lz)
  dt = 1.D0/ak/tauz
  dx = pp(l)/ak

  ta(l) = 2.D0*dt*dazm !            INSTEAD OF DAZ
  tb(l) = ta(l) + dx + 1.D0
  ub(l) = dx
  va(l) = -2.D0*dt*dbzm !           INSTEAD OF DBZ

  Return
End Subroutine

!-----------------------------------------------------------------------

Subroutine scatsol(xi, xj, fi, const, ta, tb, tc, nf)

  Use :: nlte_type
  Use :: nlte_dim
  Implicit None

! -----set up of matrix elements (note, that here in reverse order, i.e.
! from high to low frequencies)

! -----note: used scheme here is -ta tb -tc in order to use invtri

! ..
! .. scalar arguments ..
  Real (dp) :: const
  Integer (i4b) :: nf
! ..
! .. array arguments ..
  Real (dp) :: fi(nf), xi(nf), xj(nf), ta(nf), tb(nf), tc(nf)
! ..
! .. local scalars ..
  Real (dp) :: dm, ddp, dq, ta1, tb1
  Integer (i4b) :: l
! ..
! .. external subroutines ..
  External :: invtri
! ..

! first frequency

  ddp = xi(2) - xi(1)
  tb1 = 2.D0*const/ddp**2
  tb(1) = tb1 + 1.D0
  tc(1) = tb1
  fi(1) = xj(1)

  Do l = 2, nf - 1
    dm = ddp
    ddp = xi(l+1) - xi(l)
    dq = .5D0*(ddp+dm)
    ta(l) = const/(dm*dq)
    tb(l) = 2.D0*const/(ddp*dm) + 1.D0
    tc(l) = const/(ddp*dq)
    fi(l) = xj(l)
  End Do

! last freq.

  ta1 = 2.D0*const/ddp**2
  ta(nf) = ta1
  tb(nf) = ta1 + 1.D0
  fi(nf) = xj(nf)

  Call invtri(ta, tb, tc, fi, nf)

  Return
End Subroutine

!***********************************************************************

!subroutines: complex ones
!formal and related

!***********************************************************************

Subroutine formal(nd, np, nc, nb, nfobs, rmax, delta, escat, r1, r, v, opal, &
  sline, p, z, lmax, x0, vmax, opalray, sray, zray, vturbray, ltot, xcmf, &
  profile, profabs, profem, zg, aic, xmax, xmaxdop, xmaxdop_min, optiov, nsum, &
  opacon, opacray, scont, sconray, tem, xcred, xcblue, tgrad, xnue0, vdop, &
  temp, xne, profrot, vsini, deltaarr)

  Use :: nlte_type
  Use :: nlte_dim
  Use :: fund_const, Only: clight, pi
  Use :: formalsol_var, Only: xnemin, inverted, klow, kup, fcontlow, fcontup
  Implicit None


! .. parameters ..
  Real (dp), Parameter :: taucmin = 0.03 ! XNEMIN CHOSEN FROM THIS VALUE

  Integer (i4b), Parameter :: nd1 = id_ndept, ndm = id_depfi, nf = id_nfobs, &
    lto1 = 2*ndm
  Integer (i4b), Parameter :: nc1 = id_cores
  Integer (i4b), Parameter :: np1 = nc1 + 2 + 4*id_nocor
! ..
! .. scalar arguments ..
  Real (dp) :: delta, rmax, tem, tgrad, vmax, vsini, xcblue, xcred, xnue0
  Integer (i4b) :: nb, nc, nd, nfobs, np, nsum
  Logical :: escat, optiov
! ..
! .. array arguments ..
  Real (dp) :: aic(nf, 2), opacon(ndm, 2), opacray(lto1, 2), opal(ndm, nb), &
    opalray(lto1, nb), p(np1), profabs(nf), profem(nf), profile(nf), &
    profrot(nf), r(ndm), r1(nd1), sconray(lto1, 2), scont(ndm, 2), &
    sline(ndm, nb), sray(lto1, nb), vturbray(lto1), temp(ndm), v(ndm), &
    vdop(nb), x0(nf), xcmf(lto1, nb), xmax(nb), xmaxdop(nb), xmaxdop_min(nb), &
    xne(ndm), z(ndm, np1), zg(np1, nf), zray(lto1), deltaarr(nb)
  Integer (i4b) :: lmax(np1), ltot(np1), inver(nd1)
! ..
! .. local scalars ..
  Real (dp) :: aequit, dd, delp, emint, emint1, errmax, relem, tauc, tauk, &
    vplus, vvdop, w, w1, ww, ww1, xcmfmin, xicor, xxlam, xxlam1, xxx0, z2, &
    diff, diffcont
  Integer (i4b) :: i, ijp, ik, irest, istart, izb, izb1, izneb, izner, izr, &
    izr2, j, jj, jp, k, l, lm, lm1, ltaumax, ndelt
  Logical :: core, first, puredop
! ..
! .. local arrays ..
  Real (dp) :: conabs(nf), conem1(nf), conem2(nf), fscon(lto1), profcon(nf), &
    profred(nf), s1pabs(nf), s1pem(nf), s2pem(nf), tau(nf), taucon(lto1), &
    tempray(lto1), xneray(lto1)
! ..
! .. external subroutines ..
  External :: formacon, obsfram, obsfram1, obsfram2, prepray, rconv, xgrid
! ..
! .. intrinsic functions ..
! ..

  If (nd>ndm) Stop ' ND > NDM IN FORMAL'
  If (nf/=nfobs) Stop ' NF .NE. NFOBS IN FORMAL'

! no continuum definition possible, no profile will be calculated
  If (klow==-1 .And. kup==-1) Go To 130

  vvdop = 1.D50

  Do i = 1, nb
    vvdop = min(vvdop, vdop(i)) !   calculated with VTURB=VTURBMIN
  End Do

! ---- open output file for resonant zones with tau=1

! open(20,file=TRIM(FILE)//'/tauuno_'//line,status='unknown',form='formatted')

! -----sets up grid variables

! calculation of the p-grid

  Do jp = 1, nc
    p(jp) = .99D0*dble(jp-1)/dble(nc-1)
  End Do

  p(nc+1) = 1.D0
  p(nc+2) = 1.D0 + 1.D-4

  ndelt = (np-(nc+2))/4

  irest = (4*ndelt+1-nd1)
  If (irest==0) Then
    istart = 1

  Else If (irest==2) Then
    istart = -1

  Else
    Stop ' SOMETHING WRONG WITH NUMBER OF P-RAYS'
  End If

  Do jp = 1, ndelt
    ijp = nc + 2 + 4*jp
    p(ijp) = r1(nd1+1-(istart+4*jp))
    delp = p(ijp) - p(ijp-4)
    dd = delp/4.D0
    Do j = 1, 3
      p(ijp-j) = p(ijp) - j*dd
    End Do
  End Do

  p(np) = rmax

! calculation of the z-grid

  Do jp = 1, np - 1
    Do l = 1, nd
      z2 = r(l)**2 - p(jp)**2

      If (z2>=-1.D-10) Then
        If (z2<0.D0) z2 = 0.D0
        z(l, jp) = sqrt(z2)
        lmax(jp) = l
      Else
        Exit
      End If
    End Do
  End Do

  z(1, np) = 0.D0
  lmax(np) = 1

  vplus = 0.D0
  Do i = 1, nb
    vplus = max(vplus, xmax(i))
  End Do

  Print *
  Print *, ' STARK/DOPPLER BROADENING EFFECTIVE FOR XMAX = ', vplus
  Print *

  Call xgrid(delta, nfobs, xmaxdop, xmaxdop_min, x0, aic, opacon, nd, tem, &
    tgrad, xcred, xcblue, xnue0, vmax, nb, vplus, escat)

  aequit = 0.D0
  errmax = 0.D0

  Do k = 1, nfobs
    tau(k) = 0.D0
    s1pabs(k) = 0.D0
    s1pem(k) = 0.D0
    s2pem(k) = 0.D0
    conabs(k) = 0.D0
    conem1(k) = 0.D0
    conem2(k) = 0.D0
  End Do

! -----find XNEMIN from continuum depth at TAUCMIN (minimum value 10.5)
! -----tauc calculated assuming red opacity is larger

  tauc = 0.
  Do i = 1, nd - 1
    tauc = tauc + 0.5*(opacon(i,1)+opacon(i+1,1))*(r(i)-r(i+1))
    If (tauc>taucmin) Go To 100
  End Do
  Stop ' TAUC > TAUCMIN NOT FOUND!'

100 xnemin = max(3.2D10, xne(i-1))
  xnemin = 3.2D10

110 Continue

jploop: Do jp = 1, np - 1
!   JPLOOP: DO JP = 1,1
    lm = lmax(jp)
    core = lm == nd
    first = .True.

kloop: Do k = 1, nfobs
!     KLOOP: DO K = 81,81
!     KLOOP: DO K = 34,34

      xicor = aic(k, 1) + z(lm, jp)*aic(k, 2)
      zg(jp, k) = -10.D0*rmax

      Call prepray(z(1,jp), p, np, nb, jp, x0(k), ltot, lm, core, w, w1, v, &
        vvdop, r, opal, sline, opalray, sray, zray, vturbray, xcmf, delta, &
        xmaxdop, izr, izr2, izb1, izb, izner, izneb, first, optiov, nsum, &
        scont, sconray, opacon, opacray, temp, tempray, xne, xneray, escat, &
        deltaarr)

!     correction to achim's integration weights is necessary to
!     obtain absolute fluxes

!     ww=w*2.d0*pi*0.5d0/(r(1)*r(1))
!     ww1=w1*2.d0*pi*0.5d0/(r(1)*r(1))

      ww = w*pi/(r(1)*r(1))
      ww1 = w1*pi/(r(1)*r(1))

!     ---- continuum can be included in three differents ways:
!     1) keeping opa and source function constant over the whole
!     frequency range (twice v_inf);
!     2) allowing it to vary over the ray for every observer frequency;
!     3) allowing it to vary over the ray but calculating it only with
!     a certain observer frequency (the central one, for example).
!     to get a faster code, option 3) has been taken. if desired,
!     option 2) can be considered if the 'if(first)' condition is
!     avoided, and the line masked with ccc (in earlier versions)
!     is uncommented.

      xxx0 = 0.D0

      If (escat) xxx0 = x0(k)

      If (first .Or. escat) Call formacon(taucon, ltaumax, fscon, opacray, &
        sconray, zray, ltot(jp), xxx0, xcred, xcblue, escat)
!     check for xnemin

      If (first .And. jp==1) Then
        Do jj = 1, ltot(jp)
          If (taucon(jj)>taucmin) Go To 120
        End Do
        Stop ' TAUC > TAUCMIN  NOT FOUND'

120     Continue

        Print *, ' LOG NE = ', log10(xne(jj)), ' AT TAUC = ', taucon(jj)
        Print *, ' CHOSEN LOG NEMIN = ', log10(xnemin)
        If (xnemin>xne(jj) .And. xnemin/=3.2D10) Then
          Print *, ' SOMETHING WRONG WITH XNEMIN !!!'
!         reset xnemin
          xnemin = 3.2D10
          Go To 110
        Else
          Print *, ' XNEMIN OK!'
        End If
      End If
      first = .False.
      lm1 = ltot(jp)

!     ---- check whether stark zones must me considered

      If (izner==1 .And. izneb==1) Then
        puredop = .True.
      Else
        puredop = .True.
        Do ik = 1, nb
          xcmfmin = min(abs(xcmf(izner,ik)), abs(xcmf(izneb,ik)))
          If (xcmfmin<xmax(ik)) puredop = .False.
        End Do
      End If

      If (((izr==1 .And. izb==1) .Or. (izr==lm1 .And. &
        izb==lm1)) .And. puredop) Then
        If (core) Then
          emint = fscon(1)
          emint1 = xicor*exp(-taucon(lm1))
        Else
          emint = fscon(1)
          emint1 = 0.D0
        End If

      Else If (optiov) Then

        If (puredop) Then

!         -----only doppler-zones in the wind

          Call obsfram1(lm1, nb, core, emint, opalray, sray, zray, vturbray, &
            xcmf, izr, izb, emint1, zg(jp,k), tauk, xicor, taucon, fscon, &
            tempray, xneray, ltaumax)

!         -----doppler zones in the wind + stark broadening zone, or NB > 2

        Else
          Call obsfram(lm1, nb, core, emint, opalray, sray, zray, vturbray, &
            xcmf, izr, izb, izner, izneb, emint1, zg(jp,k), tauk, xicor, &
            taucon, fscon, tempray, xneray, ltaumax)

        End If

        If (core .And. jp==1) tau(k) = tauk
        If (.Not. core .And. jp==nc+2) tau(k) = tauk

      Else

!       -----separated doublets: only doppler broadening!

        Do ik = 1, nb
          If (xmax(ik)/=xmaxdop(ik)) Stop &
            ' SOMETHING WRONG WITH XMAX FOR DOPPLER-DOUBLETS!'
        End Do

        Call obsfram2(lm1, nb, core, emint, opalray, sray, zray, vturbray, &
          xcmf, izr, izr2, izb1, izb, emint1, zg(jp,k), tauk, xicor, taucon, &
          fscon, tempray, xneray, ltaumax)

        If (core .And. jp==1) tau(k) = tauk
        If (.Not. core .And. jp==nc+2) tau(k) = tauk
      End If

!     write(20,*) zg(jp,k)

      s1pabs(k) = s1pabs(k) + ww*emint1
      s1pem(k) = s1pem(k) + ww*emint
      If (core) conabs(k) = conabs(k) + ww*xicor*exp(-taucon(lm1))
      conem1(k) = conem1(k) + ww*fscon(1)

      If (jp>=nc+2) Then
        s2pem(k) = s2pem(k) + ww1*emint
        conem2(k) = conem2(k) + ww1*fscon(1)
      End If

    End Do kloop

  End Do jploop

  diffcont = 0.D0
  Do k = 1, nfobs

    profabs(k) = s1pabs(k)
    profem(k) = s1pem(k) + s2pem(k)

    If (s2pem(k)==0. .And. profem(k)==0.) Then
      relem = 0.D0
    Else
      relem = abs(s2pem(k)/profem(k))
    End If

    errmax = max(errmax, relem)
    profile(k) = profabs(k) + profem(k)
    profcon(k) = conem1(k) + conem2(k) + conabs(k)
!   COMPARE CONTINUUM FLUXES FROM HERE AND FROM FLUXCONT
    diff = min(abs(log10(profcon(k))-fcontlow), abs(log10( &
      profcon(k))-fcontup))
    If (diff>0.2) Then
      Print *, log10(profcon(k)), fcontlow, fcontup
      Print *, ' CONTINUUM FLUXES (OUTPUT VS. FLUXCONT) STRONGLY DIFFERENT'
      klow = -1
      kup = -1
      Go To 130
    End If

    diff = min(abs(1.-profcon(k)/10.**fcontlow), abs(1.-profcon( &
      k)/10.**fcontup))
    diffcont = max(diff, diffcont)
    profred(k) = profile(k)/profcon(1)

!   NOW: PROFILE CONTAINS WAVELENGTHS

    profile(k) = 1.D8/xnue0/(x0(k)*vmax/clight+1.D0)
  End Do

  Print *
  Print *, ' CALCULATED CONTINUUM FLUX (LOG) = ', log10(profcon(1))
  If (log10(profcon(1))>=min(fcontlow,fcontup) .And. log10(profcon( &
    1))<=max(fcontlow,fcontup)) Then
    Print *, ' CONT. FLUXES CONSISTENT WITH FLUXCONT (SEE ABOVE)'
  Else
    Print *, &
      ' MAXIMUM DEVIATION OF CONT. FLUXES (W.R.T. FLUXCONT BOUNDARIES) = ', &
      diffcont
  End If
  Print *


! rot. convolution

  Call rconv(profile, profred, profrot, nfobs, vsini)

! conversion to wavelength in air, xkw is wavenumber in microns,
! formula a la lang, is done in preformal/hallcl

130 Continue

  If (klow==-1 .And. kup==-1) Then
!   error condition

!   pseudo frequency grid
    diff = 2.D0/float(nfobs-1)

    Do k = 1, nfobs
      x0(k) = 1. - diff*(k-1)
      profile(k) = 1.D8/xnue0/(x0(k)*vmax/clight+1.D0)
      xxlam = profile(k)
      profcon(k) = 0.D0
      profred(k) = 0.D0
      profrot(k) = 0.D0

      Write (*, Fmt=150) k, x0(k), xxlam, profcon(k), profred(k), profrot(k)
      Write (2, Fmt=150) k, x0(k), xxlam, profcon(k), profred(k), profrot(k)
    End Do
    Print *
    Print *, ' CONTINUUM AROUND LINE(S) PROBLEMATIC (see output above)'
    Print *, ' NO PROFILE CALCULATED (set to zero)'
    Print *
    aequit = 0.D0
    Write (2, Fmt=*) aequit
  Else


    Do k = 1, nfobs
      xxlam = profile(k)

      If (k/=1) Then
        xxlam1 = profile(k-1)
        aequit = aequit + (.5D0*(profrot(k-1)+profrot(k))-1.D0)*(xxlam-xxlam1)
      End If

      Write (*, Fmt=150) k, x0(k), xxlam, profcon(k), profred(k), profrot(k)
      Write (2, Fmt=150) k, x0(k), xxlam, profcon(k), profred(k), profrot(k)
    End Do

!   #    pessimistic error (from delta i/i with step halfing)

    errmax = 15.D0*errmax
    If (errmax>.03) Print *, ' WARNING!!! WARNING!!! WARNING!!! WARNING!!!'

    Print *
    Print *, ' MAXIMUM ERROR IN ANGULAR INTEGRATION ', errmax
    Print *

    If (errmax>.03) Then
      Print *, ' WARNING!!! WARNING!!! WARNING!!! WARNING!!!'
      Print *
    End If

    Print *, ' OBSERVED AEQUIVALENT WIDTH = ', aequit
    Print *
    Write (2, Fmt=*) aequit


    Do ik = 1, nb
      Do l = 1, nd1
        If (inverted(l,ik)/=0) Go To 140
      End Do
    End Do

    Return

140 Print *, ' WARNING: INVERTED LEVELS FOUND!!!'
    Do ik = 1, nb
      inver = 0
      k = 0
      Do l = 1, nd1
        If (inverted(l,ik)/=0) Then
          k = k + 1
          inver(k) = l
        End If
      End Do
      If (k/=0) Then
        Write (*, 160) ik, (inver(l), l=1, k)
        Write (2, 170) ik, (inver(l), l=1, k)
      End If
    End Do
    Print *

  End If

  Return

150 Format (1X, I3, 5(2X,G14.6))
160 Format ('  FOR COMPONENT ', I2, ' AT L = ', 101(I2))
170 Format ('  INVERSION FOR COMPONENT ', I2, ' AT L = ', 101(I2))

End Subroutine

!-----------------------------------------------------------------------


Subroutine formacon(taucon, ltaumax, fscon, opacray, sconray, zray, ltotal, &
  x0, xcred, xcblue, escat)

  Use :: nlte_type
  Use :: nlte_dim
  Use :: formalsol_var, Only: xcmfp, lamtrans, uvlimit
  Implicit None

! ---- calculates comoving continuum opacity and source function over
! ---- the ray. it calculates also the pure continuum formal solution
! ---- (that is, the so called 'f_s(tau-cont)' for every ray.
! ---- a linear expansion of source function is taken over every
! ---- subinterval
! -----regard, that in the elscat case the source function was already
! -----interpolated to cmf frequencies in prepray
! -----

! .. parameters ..
  Integer (i4b), Parameter :: ndm = id_depfi, lto1 = 2*ndm
  Real (dp), Parameter :: taumax = 10.D0
! ..
! .. scalar arguments ..
  Real (dp) :: x0, xcblue, xcred
  Integer (i4b) :: ltaumax, ltotal
  Logical :: escat
! ..
! .. array arguments ..
  Real (dp) :: fscon(lto1), opacray(lto1, 2), sconray(lto1, 2), taucon(lto1), &
    zray(lto1)
! ..
! .. local scalars ..
  Real (dp) :: dt, dx0, dx00, dxc, e0, e1, w1, w2
  Integer (i4b) :: i
! ..
! .. local arrays ..
  Real (dp) :: opac(lto1), sco(lto1)
! ..
! .. intrinsic functions ..
! ..

  dxc = xcblue - xcred

  If (escat) Then

    Do i = 1, ltotal

!     changed as elscat assumes constant opacities

      If (dxc==0.D0) Stop ' DXC = 0 FOR ESCAT CONDITION'
      dx00 = (0.D0+xcmfp(i)-xcred)/dxc
      opac(i) = opacray(i, 1) + (opacray(i,2)-opacray(i,1))*dx00
      sco(i) = sconray(i, 1)
    End Do

  Else
    If (x0/=0.D0) Then
      Print *, ' X0 NE 0 IN FORMACON'
      Print *, &
        ' NOT POSSIBLE WITH CURRENT SETUP (ONLY ONE FREQ. POINT FOR CONT.)'
      Stop ' X0 NE 0 IN FORMACON'
    End If

    Do i = 1, ltotal

!     changed as elscat assumes constant opacities

      If (dxc==0.D0) Then
!       no interpolation
        If (opacray(i,1)/=opacray(i,2)) Stop &
          ' INCONSISTENCY DXC=0 AND OPACRAY'
        If (sconray(i,1)/=sconray(i,2)) Stop &
          ' INCONSISTENCY DXC=0 AND SCONRAY'
        opac(i) = opacray(i, 1)
        sco(i) = sconray(i, 1)
      Else
!       to avoid cmf-interpolation between (sometimes) strongly varying
!       opacities/
!       source-functions in the UV, we interpolate here only in the observer's
!       frame;
!       note that the freq. shift for the bg-elements is performed only in an
!       approximate way, so that the major part of the pseudo-cont. is in the
!       observer's frame anyway. For consistency with older versions, we keep
!       the old cmf-interpolation for optical/IR lines
        If (lamtrans>=uvlimit) Then
          dx00 = (0.D0+xcmfp(i)-xcred)/dxc
          dx0 = (x0+xcmfp(i)-xcred)/dxc
        Else
          dx00 = (0.D0-xcred)/dxc
          dx0 = (x0-xcred)/dxc
        End If
!       OLDER COMMENT: if strongly varying cont. opacities, the next statement
!       needs to
!       be changed (OPAC with DX0 as well)
!       NOW (Nov. 2021): cured via new approach to find suitable cont. points,
!       and to interpolate, in the UV, only w.r.t. the observer's frame

!       all this only for models without cmf-treatment

        opac(i) = opacray(i, 1) + (opacray(i,2)-opacray(i,1))*dx00
        sco(i) = sconray(i, 1) + (sconray(i,2)-sconray(i,1))*dx0
!       for models with cmf treatment (to be included)
!       OPAC(I) = OPACRAY(I,1)
!       SCO(I) = SCONRAY(I,1)
      End If
    End Do

  End If

  taucon(1) = 0.D0

  Do i = 2, ltotal
    taucon(i) = taucon(i-1) + 0.5D0*(opac(i-1)+opac(i))*(zray(i-1)-zray(i))
!   print*,i,taucon(i),opac(i-1),opac(i),zray(i-1),zray(i)
    If (taucon(i)<taumax) ltaumax = i
  End Do

  ltaumax = i + 1
  fscon(ltotal) = 0.D0

  Do i = ltotal - 1, 1, -1
    dt = taucon(i+1) - taucon(i)

    If (dt<1.D-8) Then
      w1 = .5D0*dt
      w2 = .5D0*dt
    Else
      e0 = exp(-dt)
      e1 = (1.D0-e0)/dt
      w1 = 1.D0 - e1
      w2 = e1 - e0
    End If

    fscon(i) = fscon(i+1) + exp(-taucon(i))*(sco(i)*w1+sco(i+1)*w2)
  End Do
  Return
End Subroutine

!-----------------------------------------------------------------------

Subroutine obsfram(lm, nb, core, emint, opalray, sray, zray, vturbray, xcmf, &
  izr, izb, izner, izneb, emint1, z1, tau, xicor, taucon, fscon, tempray, &
  xneray, ltaumax)

  Use :: nlte_type
  Use :: nlte_dim
  Implicit None


! integration of the emergent intensities in the observers frame
! for singlets and overlapping multiplets:
! pure doppler-case for ne < nemin, stark-broadening inside
! "stark-broadening zone from izner...izneb

! seven different pathes possible:two discontinuous, five conitnuous

! path 1: izr-izb, izner-izneb
! path 2: izner-izneb, izr-izb
! path 3: izr-izner-izb-izneb
! path 4: izr-izner-izneb-izb
! path 5: izner-izr-izb-izneb
! path 6: izner-izr-izneb-izb
! path 7  izner-izneb for (izr and izb) = (1 or lm)

! tauout: resonance zones of tau=tauout stored in z1


! ---- initialization

! .. parameters ..
  Integer (i4b), Parameter :: ndm = id_depfi, ltom = 2*ndm
  Real (dp), Parameter :: taumax = 10.D0, tauout = 1.D0
! ..
! .. scalar arguments ..
  Real (dp) :: emint, emint1, tau, xicor, z1
  Integer (i4b) :: izb, izneb, izner, izr, lm, ltaumax, nb
  Logical :: core
! ..
! .. array arguments ..
  Real (dp) :: fscon(ltom), opalray(ltom, nb), sray(ltom, nb), taucon(ltom), &
    tempray(ltom), xcmf(ltom, nb), xneray(ltom), zray(ltom), vturbray(ltom)
! ..
! .. local scalars ..
  Real (dp) :: dt, fblue, fq, fred, opablue, opai, opais, opared, phii, sblue, &
    si, sis, sred, taublu, taured, xi, zi, zinew
  Integer (i4b) :: iend, iend1, iend2, ista, ista1, ista2, k, l
  Logical :: flag, optout
! ..
! .. external functions ..
  Real (dp) :: expuno, perfil1
  External :: expuno, perfil1
! ..
! .. intrinsic functions ..
! ..

  optout = .True.
  emint = .0D0
  tau = .0D0
  sis = 0.D0

  opais = 0.D0

! ---- decision, whether continuous or discontinuous pathes


  If (izr==1 .And. izb==1 .Or. izr==lm .And. izb==lm) Then

!   --path7
    ista = izner
    iend = izneb

!   --path1
  Else If (izb<izner) Then
    ista1 = izr
    iend1 = izb
    ista2 = izner
    iend2 = izneb
    Go To 120

!   --path2
  Else If (izneb<izr) Then
    ista1 = izner
    iend1 = izneb
    ista2 = izr
    iend2 = izb
    Go To 120

!   -- remaining pathes
  Else
    ista = min(izr, izner)
    iend = max(izb, izneb)
  End If

! -- continous path

  l = ista
  zi = zray(l)

  Do k = 1, nb
    xi = xcmf(l, k)
    phii = perfil1(k, xi, tempray(l), xneray(l), vturbray(l), flag)
    If (flag) Cycle
    opai = phii*opalray(l, k)
    si = opai*sray(l, k)
    opais = opais + opai
    sis = sis + si
  End Do

  opared = opais
  taured = tau

  If (opais/=0.D0) Then
    sred = sis/opais
    fred = sred*exp(-taucon(l)) - fscon(l)
  Else
    fred = -fscon(l)
  End If

  l = l + 1

100 Continue

  sis = 0.D0
  opais = 0.D0
  zinew = zray(l)

  Do k = 1, nb
    xi = xcmf(l, k)
    phii = perfil1(k, xi, tempray(l), xneray(l), vturbray(l), flag)
    If (flag) Cycle
    opai = phii*opalray(l, k)
    si = opai*sray(l, k)
    opais = opais + opai
    sis = sis + si
  End Do

  opablue = opais

  If (opais/=0.D0) Then
    sblue = sis/opais
    fblue = sblue*exp(-taucon(l)) - fscon(l)
  Else
    fblue = -fscon(l)
  End If

  dt = .5D0*(opared+opablue)*(zi-zinew)
  taublu = taured + dt

  If (dt==0.D0) Then
!   if(nb.eq.1) print*,' warning!! error in resonance zone - obsfram'
    Go To 110
  End If

! print*,zi,'  ',zinew,'  ',tau,'  ',dt,'  macrogrid'

  If (dt>1.D0) Then
    fq = fred !                     exact result would be fq=f(dt=1)
  Else
    fq = .5D0*(fblue+fred)
!   almost exact for small dt; when expanded with f=a+b*t,
!   exact result (for small dt) would be
!   fq=a+b*dt/2/(1+dt/2) approx a+b*dt/2 ->
!   fq approx fred+(fblu-fred)/dt * dt/2 = 0.5(fred+fblu)
  End If
! tested by comparison with linear expansion of f and integration
! almost perfect. Note that cancellation effects need to be accounted
! for when f linearly expanded

! This is the test
! IF (DT.GT.0.01D0) THEN ! for smaller dt, cancellation effects
! FQ = (FBLUE-FRED)/DT
! AUX=FRED*(1.D0-EXP(-DT))+FQ*(1.D0-EXP(-DT)*(1.D0+DT))
! EMINT = EMINT + EXP(-TAURED)*AUX
! ELSE
! FQ = .5D0* (FBLUE+FRED)
! EMINT = EMINT + EXP(-TAURED)*FQ*EXPUNO(DT)
! END IF

  emint = emint + exp(-taured)*fq*expuno(dt)

  If (emint>10.D0) Stop ' EMINT > 10 IN OBSFRAM'

  tau = tau + dt

  If (optout .And. tau>=tauout) Then
    z1 = zinew
    optout = .False.
  End If

110 Continue

  If (tau>taumax .Or. l==iend .Or. l==ltaumax) Then
    emint1 = 0.D0
    If (core) emint1 = xicor*exp(-tau-taucon(lm))
    emint = emint + fscon(1)
    Return
  End If

  opared = opablue
  fred = fblue
  taured = taublu
  zi = zinew
  l = l + 1
  Go To 100

! -- discontinous paths

! -- 1st inegration range

120 Continue
  l = ista1
  zi = zray(l)

  Do k = 1, nb
    xi = xcmf(l, k)
    phii = perfil1(k, xi, tempray(l), xneray(l), vturbray(l), flag)
    If (flag) Cycle
    opai = phii*opalray(l, k)
    si = opai*sray(l, k)
    opais = opais + opai
    sis = sis + si
  End Do

  opared = opais
  taured = tau

  If (opais/=0.D0) Then
    sred = sis/opais
    fred = sred*exp(-taucon(l)) - fscon(l)
  Else
    fred = -fscon(l)
  End If

  l = l + 1

130 Continue
  sis = 0.D0
  opais = 0.D0
  zinew = zray(l)

  Do k = 1, nb
    xi = xcmf(l, k)
    phii = perfil1(k, xi, tempray(l), xneray(l), vturbray(l), flag)
    If (flag) Cycle
    opai = phii*opalray(l, k)
    si = opai*sray(l, k)
    opais = opais + opai
    sis = sis + si
  End Do

  opablue = opais

  If (opais/=0.D0) Then
    sblue = sis/opais
    fblue = sblue*exp(-taucon(l)) - fscon(l)
  Else
    fblue = -fscon(l)
  End If

  dt = .5D0*(opared+opablue)*(zi-zinew)
  taublu = taured + dt

  If (dt==0.D0) Then
!   if(nb.eq.1) print*,' warning!! error in resonance zone - obsfram'
    Go To 140
  End If

! print*,zi,'  ',zinew,'  ',tau,'  ',dt,'  macrogrid'

! for tests, see routine obsfram (first incidence of integration)
  If (dt>1.D0) Then
    fq = fred
  Else
    fq = .5D0*(fblue+fred)
  End If


  emint = emint + exp(-taured)*fq*expuno(dt)

  If (emint>10.D0) Stop ' EMINT > 10 IN OBSFRAM'

  tau = tau + dt
  If (optout .And. tau>=tauout) Then
    z1 = zinew
    optout = .False.
  End If

  If (tau>taumax .Or. l==ltaumax) Then
    emint1 = 0.D0
    emint = emint + fscon(1)
    Return
  End If

140 Continue
  If (l==iend1) Go To 150

  opared = opablue
  fred = fblue
  taured = taublu
  zi = zinew
  l = l + 1
  Go To 130

! ------------------------------

! -- 2nd integration range!

150 Continue

  If (ista2==lm) Then
    If (iend2/=lm) Stop ' ERROR IN OBSFRAM'
    emint1 = 0.D0
    If (core) emint1 = xicor*exp(-tau-taucon(lm))
    emint = emint + fscon(1)
    Return
  End If

  taured = tau
  l = ista2
  zi = zray(l)

  Do k = 1, nb
    xi = xcmf(l, k)
    phii = perfil1(k, xi, tempray(l), xneray(l), vturbray(l), flag)
    If (flag) Cycle
    opai = phii*opalray(l, k)
    si = opai*sray(l, k)
    opais = opais + opai
    sis = sis + si
  End Do

  opared = opais
  taured = tau

  If (opais/=0.D0) Then
    sred = sis/opais
    fred = sred*exp(-taucon(l)) - fscon(l)
  Else
    fred = -fscon(l)
  End If

  l = l + 1

160 Continue

  sis = 0.D0
  opais = 0.D0
  zinew = zray(l)

  Do k = 1, nb
    xi = xcmf(l, k)
    phii = perfil1(k, xi, tempray(l), xneray(l), vturbray(l), flag)
    If (flag) Cycle
    opai = phii*opalray(l, k)
    si = opai*sray(l, k)
    opais = opais + opai
    sis = sis + si
  End Do

  opablue = opais

  If (opais/=0.D0) Then
    sblue = sis/opais
    fblue = sblue*exp(-taucon(l)) - fscon(l)
  Else
    fblue = -fscon(l)
  End If

  dt = .5D0*(opared+opablue)*(zi-zinew)
  taublu = taured + dt

  If (dt==0.D0) Then
!   if(nb.eq.1) print*,' warning!! error in resonance zone - obsfram'
    Go To 170
  End If

! print*,zi,'  ',zinew,'  ',tau,'  ',dt,'  macrogrid'

! for tests, see routine obsfram (first incidence of integration)
  If (dt>1.D0) Then
    fq = fred
  Else
    fq = .5D0*(fblue+fred)
  End If


  emint = emint + exp(-taured)*fq*expuno(dt)

  If (emint>10.D0) Stop ' EMINT > 10 IN OBSFRAM'

  tau = tau + dt
  If (optout .And. tau>=tauout) Then
    z1 = zinew
    optout = .False.
  End If

170 Continue

  If (tau>taumax .Or. l==iend2 .Or. l==ltaumax) Then
    emint1 = 0.D0
    If (core) emint1 = xicor*exp(-tau-taucon(lm))
    emint = emint + fscon(1)
    Return
  End If

  opared = opablue
  fred = fblue
  taured = taublu
  zi = zinew
  l = l + 1
  Go To 160

End Subroutine

!-----------------------------------------------------------------------

Subroutine obsfram1(lm, nb, core, emint, opalray, sray, zray, vturbray, xcmf, &
  izr, izb, emint1, z1, tau, xicor, taucon, fscon, tempray, xneray, ltaumax)

  Use :: nlte_type
  Use :: nlte_dim
  Implicit None


! integration of the emergent intensities in the observers frame
! for singlets and overlapping multiplets, pure doppler-case

! tauout: resonance zones of tau=tauout stored in z1


! .. parameters ..
  Integer (i4b), Parameter :: ndm = id_depfi, ltom = 2*ndm
  Real (dp), Parameter :: taumax = 10.D0, tauout = 1.D0
! ..
! .. scalar arguments ..
  Real (dp) :: emint, emint1, tau, xicor, z1
  Integer (i4b) :: izb, izr, lm, ltaumax, nb
  Logical :: core
! ..
! .. array arguments ..
  Real (dp) :: fscon(ltom), opalray(ltom, nb), sray(ltom, nb), taucon(ltom), &
    tempray(ltom), xcmf(ltom, nb), xneray(ltom), zray(ltom), vturbray(ltom)
! ..
! .. local scalars ..
  Real (dp) :: dt, fblue, fq, fred, opablue, opai, opais, opared, phii, sblue, &
    si, sis, sred, taublu, taured, xi, zi, zinew
  Integer (i4b) :: k, l
  Logical :: flag, optout
! ..
! .. external functions ..
  Real (dp) :: expuno, perfil1
  External :: expuno, perfil1
! ..
! .. intrinsic functions ..
! ..

  optout = .True.

  emint = .0D0
  tau = .0D0

  sis = 0.D0
  opais = 0.D0
  l = izr
  zi = zray(l)

  Do k = 1, nb
    xi = xcmf(l, k)
    phii = perfil1(k, xi, tempray(l), xneray(l), vturbray(l), flag)
    If (flag) Cycle
    opai = phii*opalray(l, k)
    si = opai*sray(l, k)
    opais = opais + opai
    sis = sis + si
  End Do

  opared = opais
  taured = tau
  If (opais/=0.D0) Then
    sred = sis/opais
    fred = sred*exp(-taucon(l)) - fscon(l)
  Else
    fred = -fscon(l)
  End If

  l = l + 1

100 Continue

  sis = 0.D0
  opais = 0.D0
  zinew = zray(l)

  Do k = 1, nb
    xi = xcmf(l, k)
    phii = perfil1(k, xi, tempray(l), xneray(l), vturbray(l), flag)
    If (flag) Cycle
    opai = phii*opalray(l, k)
    si = opai*sray(l, k)
    opais = opais + opai
    sis = sis + si
  End Do

  opablue = opais
  If (opais/=0.D0) Then
    sblue = sis/opais
    fblue = sblue*exp(-taucon(l)) - fscon(l)
  Else
    fblue = -fscon(l)
  End If

  dt = .5D0*(opared+opablue)*(zi-zinew)
  taublu = taured + dt

  If (dt==0.D0) Then
!   if(nb.eq.1) print*,' warning!! error in resonance zone - obsfram1'
    Go To 110
  End If

! print*,zi,' ',zinew,' ',tau,' ',dt,'  macrogrid'

! aux1=exp(-taured)-exp(-taublu)
! aux2=taured*exp(-taured)-taublu*exp(-taublu)
! auxa=(fblue-fred)/dt
! auxb=(fred*taublu-fblue*taured)/dt

! emint=emint+aux1*(auxa+auxb)+auxa*aux2
! for tests, see routine obsfram (first incidence of integration)
  If (dt>1.D0) Then
    fq = fred
  Else
    fq = .5D0*(fblue+fred)
  End If


  emint = emint + exp(-taured)*fq*expuno(dt)

  If (emint>10.D0) Stop ' EMINT > 10 IN OBSFRAM1'

! write(*,1001) l,zi,zinew,tau,dt,sblue,fscon(l),taucon(l),fq,emint
! 1001 format(i4,2x,2(f8.3,2x),7(e9.3,2x))

  tau = tau + dt

  If (optout .And. tau>=tauout) Then
    z1 = zinew
    optout = .False.
  End If

110 Continue

  If (tau>taumax .Or. l==izb .Or. l==ltaumax) Then
    emint1 = 0.D0
    If (core) emint1 = xicor*exp(-tau-taucon(lm))
    emint = emint + fscon(1)
!   print*,'emint ',emint
    Return
  End If

  opared = opablue
  fred = fblue
  taured = taublu
  zi = zinew
  l = l + 1
  Go To 100

End Subroutine

!-----------------------------------------------------------------------

Subroutine obsfram2(lm, nb, core, emint, opalray, sray, zray, vturbray, xcmf, &
  izr, izr2, izb1, izb, emint1, z1, tau, xicor, taucon, fscon, tempray, &
  xneray, ltaumax)

  Use :: nlte_type
  Use :: nlte_dim
  Implicit None

! integration of the emergent intensities in the observers frame
! for non-overlapping doublets, pure doppler-case

! tauout: resonance zones of tau=tauout stored in z1


! .. parameters ..
  Integer (i4b), Parameter :: ndm = id_depfi, ltom = 2*ndm
  Real (dp), Parameter :: taumax = 10.D0, tauout = 1.D0
! ..
! .. scalar arguments ..
  Real (dp) :: emint, emint1, tau, xicor, z1
  Integer (i4b) :: izb, izb1, izr, izr2, lm, ltaumax, nb
  Logical :: core
! ..
! .. array arguments ..
  Real (dp) :: fscon(ltom), opalray(ltom, nb), sray(ltom, nb), taucon(ltom), &
    tempray(ltom), xcmf(ltom, nb), xneray(ltom), zray(ltom), vturbray(ltom)
! ..
! .. local scalars ..
  Real (dp) :: dt, fblue, fq, fred, opablue, opared, phii, sblue, sred, &
    taublu, taured, xi, zi, zinew
  Integer (i4b) :: k, l
  Logical :: flag, optout
! ..
! .. external functions ..
  Real (dp) :: expuno, perfil1
  External :: expuno, perfil1
! ..
! .. intrinsic functions ..
! ..

  If (nb>2) Stop ' NB > 2 AND OBSFRAM2!'

  optout = .True.

  emint = .0D0
  tau = .0D0

  If (izr==1 .And. izr2==1) Go To 120

! red component!!!

  taured = tau
  l = izr
  zi = zray(l)

  k = 2
  xi = xcmf(l, 2)
  phii = perfil1(k, xi, tempray(l), xneray(l), vturbray(l), flag)

  If (flag) Then
    fred = -fscon(l)
    opared = 0.D0
  Else
    opared = phii*opalray(l, 2)
    sred = sray(l, 2)
    fred = sred*exp(-taucon(l)) - fscon(l)
  End If


  l = l + 1

100 Continue

  zinew = zray(l)

  xi = xcmf(l, 2)
  phii = perfil1(k, xi, tempray(l), xneray(l), vturbray(l), flag)

  If (flag) Then
    fblue = -fscon(l)
    opablue = 0.D0
  Else
    opablue = phii*opalray(l, 2)
    sblue = sray(l, 2)
    fblue = sblue*exp(-taucon(l)) - fscon(l)
  End If


  dt = .5D0*(opared+opablue)*(zi-zinew)
  taublu = taured + dt

  If (dt==0.D0) Then
!   print*,' warning!! error in resonance zone - obsfram2'
    Go To 110
  End If

! print*,zi,'  ',zinew,'  ',tau,'  ',dt,'  macrogrid'

! aux1=exp(-taured)-exp(-taublu)
! aux2=taured*exp(-taured)-taublu*exp(-taublu)
! auxa=(fblue-fred)/dt
! auxb=(fred*taublu-fblue*taured)/dt

! emint=emint+aux1*(auxa+auxb)+auxa*aux2
! for tests, see routine obsfram (first incidence of integration)
  If (dt>1.D0) Then
    fq = fred
  Else
    fq = 0.5D0*(fred+fblue)
  End If


  emint = emint + exp(-taured)*fq*expuno(dt)

  If (emint>10.D0) Stop ' EMINT > 10 IN OBSFRAM2'

  tau = tau + dt

  If (optout .And. tau>=tauout) Then
    z1 = zinew
    optout = .False.
  End If

  If (tau>taumax .Or. l==ltaumax) Then
    emint1 = 0.D0
    emint = emint + fscon(1)
    Return
  End If

110 Continue

  If (l==izr2) Go To 120

  fred = fblue
  taured = taublu
  opared = opablue
  zi = zinew
  l = l + 1
  Go To 100

! ------------------------------

! blue component!!!

120 Continue

  If (izb1==lm) Then
    If (izb/=lm) Stop ' ERROR IN OBSFRAM2'
    emint1 = 0.D0
    If (core) emint1 = xicor*exp(-tau-taucon(lm))
    emint = emint + fscon(1)
    Return
  End If

  taured = tau
  l = izb1
  zi = zray(l)

  k = 1
  xi = xcmf(l, 1)
  phii = perfil1(k, xi, tempray(l), xneray(l), vturbray(l), flag)

  If (flag) Then
    fred = -fscon(l)
    opared = 0.D0
  Else
    opared = phii*opalray(l, 1)
    sred = sray(l, 1)
    fred = sred*exp(-taucon(l)) - fscon(l)
  End If


  l = l + 1

130 Continue

  zinew = zray(l)

  xi = xcmf(l, 1)
  phii = perfil1(k, xi, tempray(l), xneray(l), vturbray(l), flag)

  If (flag) Then
    fblue = -fscon(l)
    opablue = 0.D0
  Else
    opablue = phii*opalray(l, 1)
    sblue = sray(l, 1)
    fblue = sblue*exp(-taucon(l)) - fscon(l)
  End If


  dt = .5D0*(opared+opablue)*(zi-zinew)
  taublu = taured + dt

  If (dt==0.D0) Then
!   print*,' warning!! error in resonance zone - obsfram2'
    Go To 140
  End If

! print*,zi,'  ',zinew,'  ',tau,'  ',dt,'  macrogrid'

! aux1=exp(-taured)-exp(-taublu)
! aux2=taured*exp(-taured)-taublu*exp(-taublu)
! auxa=(fblue-fred)/dt
! auxb=(fred*taublu-fblue*taured)/dt

! emint=emint+aux1*(auxa+auxb)+auxa*aux2
! for tests, see routine obsfram (first incidence of integration)
  If (dt>1.D0) Then
    fq = fred
  Else
    fq = 0.5D0*(fred+fblue)
  End If


  emint = emint + exp(-taured)*fq*expuno(dt)

  If (emint>10.D0) Stop ' EMINT > 10 IN OBSFRAM2'

  tau = tau + dt
  If (optout .And. tau>=tauout) Then
    z1 = zinew
    optout = .False.
  End If

140 Continue

  If (tau>taumax .Or. l==izb .Or. l==ltaumax) Then
    emint1 = 0.D0
    If (core) emint1 = xicor*exp(-tau-taucon(lm))
    emint = emint + fscon(1)
    Return
  End If

  fred = fblue
  taured = taublu
  opared = opablue
  zi = zinew
  l = l + 1
  Go To 130

End Subroutine

!-----------------------------------------------------------------------

Subroutine prepray(z, p, np, nb, jp, x0, ltot, lmax, core, w, w1, v, vvdop, r, &
  opal, source, opalray, sray, zray, vturbray, xcmf, delta, xmax, izr, izr2, &
  izb1, izb, izner, izneb, first, optiov, nsum, scont, sconray, opacon, &
  opacray, temp, tempray, xne, xneray, escat, deltaarr)

  Use :: nlte_type
  Use :: nlte_dim
  Use :: formalsol_var, Only: xcmfp, vzrp, wp, wp1, nfcmf, xcmfe, sconte, &
    sceray, sce, desce, xnemin, xnue, vturbv

  Implicit None

! defining whole rays including backward half-sphere
! calculation of doppler-(freq. dependent) and stark zones

! note here: xmax corresponds to doppler-broadening

! .. parameters ..
  Integer (i4b), Parameter :: ndm = id_depfi, lto1 = 2*ndm, nc = id_cores
  Integer (i4b), Parameter :: np1 = nc + 2 + 4*id_nocor
! ..
! .. scalar arguments ..
  Real (dp) :: delta, vvdop, w, w1, x0
  Integer (i4b) :: izb, izb1, izneb, izner, izr, izr2, jp, lmax, nb, np, nsum
  Logical :: core, escat, first, optiov
! ..
! .. array arguments ..
  Real (dp) :: opacon(ndm, 2), opacray(lto1, 2), opal(ndm, nb), &
    opalray(lto1, nb), p(np1), r(ndm), sconray(lto1, 2), scont(ndm, 2), &
    source(ndm, nb), sray(lto1, nb), temp(ndm), tempray(lto1), v(ndm), &
    xcmf(lto1, nb), xmax(nb), xne(ndm), xneray(lto1), z(ndm), zray(lto1), &
    vturbray(lto1), deltaarr(nb)

  Real (dp), Dimension (10) :: x01, opai2, delopa2, si2, dels2

  Integer (i4b) :: ltot(np1)
! ..
! .. local scalars ..
  Real (dp) :: dellin, delopa1, dels1, delt1, delvt1, deltaz, delv, delxn1, &
    delz, delzi, deopc1, deopc2, desc1, desc2, ddp, dpold, dpp3, dpp3old, &
    opai1, opc1, opc2, q1, sc1, sc2, si1, sum1, sum2, sum3, sumex, t1, vt1, &
    vdopmax, vdopmin, vzr, vzr1, xcmfl, xmaxi, xn1, z1, zi
  Integer (i4b) :: i, ip, irest, k, kk, l, l0, l1, ll, lm, lm1, ndel
  Logical :: optz0
! ..
! .. intrinsic functions ..
! ..

  If (first) Then

    vdopmin = vvdop/5.D0
    vdopmax = vdopmin*2.D0
!   modified from V6.0 on
    If (nb>=2) Then
      If (xnue(1)/xnue(nb)<=1.D0) Stop ' PREFORMAL: ERROR IN XNUE'
      If (nb>10) Stop ' TOO MANY COMPONENTS (> 10) IN PREFORMAL'
      Do i = 2, nb
        x01(i) = (x0+deltaarr(i))*xnue(1)/xnue(i)
      End Do
    End If

    optz0 = .False.

!   -----core rays and first point of non-core rays without interpolation

    If (core) Then
      ltot(jp) = lmax
      lm1 = lmax
    Else
      lm1 = 1
    End If

    Do l = 1, lm1
      vzr = v(l)*z(l)/r(l)
      xcmfp(l) = -vzr
      xcmf(l, 1) = x0 - vzr
      sray(l, 1) = source(l, 1)
      opalray(l, 1) = opal(l, 1)
      opacray(l, 1) = opacon(l, 1)
      sconray(l, 1) = scont(l, 1)
      opacray(l, 2) = opacon(l, 2)
      sconray(l, 2) = scont(l, 2)
      tempray(l) = temp(l)
      xneray(l) = xne(l)
      zray(l) = z(l)
      vturbray(l) = vturbv(l)
    End Do

    If (escat) Then
      Do k = 1, nfcmf
        Do l = 1, lm1
          sceray(l, k) = sconte(l, k)
        End Do
      End Do
    End If

    If (nb>=2) Then
      Do l = 1, lm1
        Do k = 2, nb
          xcmf(l, k) = x01(k) + xcmfp(l)
          sray(l, k) = source(l, k)
          opalray(l, k) = opal(l, k)
        End Do
      End Do
    End If

!   -----non core rays with interpolation

    If (.Not. core) Then

      If (z(lmax)==0.D0) Then
        optz0 = .True.
        lm1 = lmax
      Else
        lm1 = lmax + 1
      End If

      Do l = 1, lmax
        vzrp(l) = v(l)*z(l)/r(l)
      End Do

      If (.Not. optz0) vzrp(lm1) = 0.D0

      i = 1

lloop1: Do l = 1, lm1 - 1
        l1 = l + 1
        vzr = vzrp(l)
        vzr1 = vzrp(l1)
        delv = vzr - vzr1

        If (delv<vdopmax) Then

!         no interpolation necessary

          i = i + 1
          If (i>lto1) Stop ' MORE THAN LTO1 POINTS USED IN INTERPOLATION!'
          xcmfp(i) = -vzr1

          If (optz0 .Or. l1/=lm1) Then
            xcmf(i, 1) = x0 - vzr1
            sray(i, 1) = source(l1, 1)
            opalray(i, 1) = opal(l1, 1)
            opacray(i, 1) = opacon(l1, 1)
            opacray(i, 2) = opacon(l1, 2)
            sconray(i, 1) = scont(l1, 1)
            sconray(i, 2) = scont(l1, 2)

            If (escat) Then
              Do k = 1, nfcmf
                sceray(i, k) = sconte(l1, k)
              End Do
            End If

            tempray(i) = temp(l1)
            xneray(i) = xne(l1)
            zray(i) = z(l1)
            vturbray(i) = vturbv(l1)

            If (nb>=2) Then
              Do k = 2, nb
                xcmf(i, k) = x01(k) - vzr1
                sray(i, k) = source(l1, k)
                opalray(i, k) = opal(l1, k)
              End Do
            End If

!           different treatment for last point with optz0 = false, assuming
!           constant source function and opacity from z(lmax) to z=0.

          Else
            xcmf(i, 1) = x0 - vzr1
            sray(i, 1) = source(l, 1)
            opalray(i, 1) = opal(l, 1)
            opacray(i, 1) = opacon(l, 1)
            opacray(i, 2) = opacon(l, 2)
            sconray(i, 1) = scont(l, 1)
            sconray(i, 2) = scont(l, 2)

            If (escat) Then
              Do k = 1, nfcmf
                sceray(i, k) = sconte(l, k)
              End Do
            End If

            tempray(i) = temp(l)
            xneray(i) = xne(l)
            vturbray(i) = vturbv(l)
            zray(i) = 0.D0
            If (nb>=2) Then
              Do k = 2, nb
                xcmf(i, k) = x01(k) - vzr1
                sray(i, k) = source(l, k)
                opalray(i, k) = opal(l, k)
              End Do
            End If

          End If

!         interpolation linearly with z

        Else

          ndel = int(delv/vdopmin) + 1
          If (ndel<2) Stop ' NDEL < 2'
          zi = z(l)

          If (optz0 .Or. l1/=lm1) Then
            z1 = z(l1)
            deltaz = zi - z1
            delz = deltaz/dble(ndel)
            deltaz = 1.D0/deltaz
            opai1 = opal(l, 1)
            delopa1 = opal(l1, 1) - opal(l, 1)
            si1 = source(l, 1)
            dels1 = source(l1, 1) - source(l, 1)
            opc1 = opacon(l, 1)
            deopc1 = opacon(l1, 1) - opacon(l, 1)
            sc1 = scont(l, 1)
            desc1 = scont(l1, 1) - scont(l, 1)
            opc2 = opacon(l, 2)
            deopc2 = opacon(l1, 2) - opacon(l, 2)
            sc2 = scont(l, 2)
            desc2 = scont(l1, 2) - scont(l, 2)

            If (escat) Then
              Do k = 1, nfcmf
                sce(k) = sconte(l, k)
                desce(k) = sconte(l1, k) - sconte(l, k)
              End Do
            End If

            t1 = temp(l)
            delt1 = temp(l1) - temp(l)
            xn1 = xne(l)
            delxn1 = xne(l1) - xne(l)
            vt1 = vturbv(l)
            delvt1 = vturbv(l1) - vturbv(l)
            If (nb>=2) Then
              Do k = 2, nb
                opai2(k) = opal(l, k)
                delopa2(k) = opal(l1, k) - opal(l, k)
                si2(k) = source(l, k)
                dels2(k) = source(l1, k) - source(l, k)
              End Do
            End If

          Else

!           different treatment for last point with optz0 = false, assuming
!           constant source function and opacity from z(lmax) to z=0.

            z1 = 0.D0
            deltaz = zi - z1
            delz = deltaz/dble(ndel)
            deltaz = 1.D0/deltaz
            opai1 = opal(l, 1)
            delopa1 = 0.D0
            si1 = source(l, 1)
            dels1 = 0.D0
            opc1 = opacon(l, 1)
            deopc1 = 0.D0
            sc1 = scont(l, 1)
            desc1 = 0.D0
            opc2 = opacon(l, 2)
            deopc2 = 0.D0
            sc2 = scont(l, 2)
            desc2 = 0.D0

            If (escat) Then
              Do k = 1, nfcmf
                sce(k) = sconte(l, k)
                desce(k) = 0.D0
              End Do
            End If

            t1 = temp(l)
            delt1 = 0.D0
            xn1 = xne(l)
            delxn1 = 0.D0
            vt1 = vturbv(l)
            delvt1 = 0.D0
            If (nb>=2) Then
              Do k = 2, nb
                opai2(k) = opal(l, k)
                delopa2(k) = 0.D0
                si2(k) = source(l, k)
                dels2(k) = 0.D0
              End Do
            End If

          End If

kloop:    Do k = 1, ndel
            i = i + 1
            If (i>lto1) Stop ' MORE THAN LTO1 POINTS USED IN INTERPOLATION!'
            delzi = k*delz
            dellin = delzi*deltaz
            zray(i) = zi - delzi

!           attention!! explicitly accounted for sign of delv

            xcmfp(i) = -vzr + delv*dellin
            xcmf(i, 1) = x0 + xcmfp(i)
            sray(i, 1) = si1 + dels1*dellin
            opalray(i, 1) = opai1 + delopa1*dellin
            opacray(i, 1) = opc1 + deopc1*dellin
            opacray(i, 2) = opc2 + deopc2*dellin
            sconray(i, 1) = sc1 + desc1*dellin
            sconray(i, 2) = sc2 + desc2*dellin

            If (escat) Then
              Do kk = 1, nfcmf
                sceray(i, kk) = sce(kk) + desce(kk)*dellin
              End Do
            End If

            tempray(i) = t1 + delt1*dellin
            xneray(i) = xn1 + delxn1*dellin
            vturbray(i) = vt1 + delvt1*dellin

            If (nb>=2) Then
              Do kk = 2, nb
                xcmf(i, kk) = x01(kk) + xcmfp(i)
                sray(i, kk) = si2(kk) + dels2(kk)*dellin
                opalray(i, kk) = opai2(kk) + delopa2(kk)*dellin
              End Do
            End If

          End Do kloop

        End If

      End Do lloop1

!     defining quantities on the back hemisphere

      lm1 = 2*i - 1
      If (lm1>lto1) Stop ' MORE THAN LTO1 POINTS USED IN INTERPOLATION!'
      ltot(jp) = lm1

      Do l = 1, i - 1
        ll = 2*i - l
        vzr = -xcmfp(l)
        xcmfp(ll) = vzr
        xcmf(ll, 1) = x0 + vzr
        sray(ll, 1) = sray(l, 1)
        opalray(ll, 1) = opalray(l, 1)
        opacray(ll, 1) = opacray(l, 1)
        opacray(ll, 2) = opacray(l, 2)
        sconray(ll, 1) = sconray(l, 1)
        sconray(ll, 2) = sconray(l, 2)
        tempray(ll) = tempray(l)
        xneray(ll) = xneray(l)
        zray(ll) = -zray(l)
        vturbray(ll) = vturbray(l)
      End Do

      If (escat) Then
        Do k = 1, nfcmf
          Do l = 1, i - 1
            ll = 2*i - l
            sceray(ll, k) = sceray(l, k)
          End Do
        End Do
      End If

      If (nb>=2) Then
        Do l = 1, i - 1
          ll = 2*i - l
          Do k = 2, nb
            xcmf(ll, k) = x01(k) + xcmfp(ll)
            sray(ll, k) = sray(l, k)
            opalray(ll, k) = opalray(l, k)
          End Do
        End Do
      End If

!     end of non-core ray treatment

    End If

!   w = weights for flux integral

    If (jp==1) Then
      w = p(2)*p(2)/3.D0
      w1 = 0.D0

    Else If (jp<nc+2) Then
      w = (p(jp-1)+p(jp)+p(jp+1))*(p(jp+1)-p(jp-1))/3.D0
      w1 = 0.D0

    Else If (jp==nc+2) Then
      w = (p(jp)-p(jp-1))*(2*p(jp)+p(jp-1))/3.D0
      ddp = p(jp+1) - p(jp)
      dpp3 = p(jp)*ddp*2.D0/3.D0
      w = w + dpp3
      w1 = -dpp3/15.D0

    Else
      ip = jp - (nc+2)
      irest = mod(ip, 4)
      ddp = p(jp+1) - p(jp)
      dpp3 = p(jp)*ddp*2.D0/3.D0

      If (irest==0) Then
        dpold = p(jp) - p(jp-1)
        dpp3old = p(jp)*dpold*2.D0/3.D0
        w = dpp3 + dpp3old
        w1 = -w/15.D0

      Else If (irest==1 .Or. irest==3) Then
        w = 4.D0*dpp3
        w1 = w/15.D0
      Else
        w = 2.D0*dpp3
        w1 = -3.D0*w/15.D0
      End If

    End If

    wp(jp) = w
    wp1(jp) = w1

    If (jp==np-1) Then
      ddp = p(jp+1) - p(jp)
      dpp3 = p(jp+1)*ddp*2.D0/3.D0
      wp(np) = dpp3
      wp1(np) = -dpp3/15.D0

!     -----test of integration weights

      sum1 = 0.D0
      sum2 = 0.D0
      sum3 = 0.D0
      Do ip = 1, np
        sum1 = sum1 + wp(ip)
        sum3 = sum3 + wp1(ip)
        If (ip>=nc+2) sum2 = sum2 + wp1(ip)/p(ip)
      End Do

      If (abs(sum2)>1.D-5) Stop 'ERROR IN SUM2'
      sumex = r(1)**2

      If (abs(sum1+sum3-sumex)>1.D-2) Then
        Print *, sum1 + sum3 - sumex
        Stop 'ERROR IN SUM1'
      End If
    End If

    izner = 1
    izneb = 1

!   definition of "stark-broadening zone"

    If (nsum>0) Then

      If (core) Then
        Do l = 1, ltot(jp)
          If (xneray(l)>=xnemin) Go To 100
        End Do
        Stop 'XNEMIN NOT FOUND IN CORE-RAYS'

100     Continue
        izner = l - 1
        izneb = ltot(jp)

      Else
        l0 = (ltot(jp)-1)/2
        Do l = 1, l0
          If (xneray(l)>=xnemin) Go To 110
        End Do
!       print*,'xnemin not found in non core-ray with p =',p(jp)
        Go To 130

110     Continue

        izner = l - 1
        Do l = ltot(jp), l0, -1
          If (xneray(l)>=xnemin) Go To 120
        End Do
        Stop 'XNEMIN NOT FOUND IN NON CORE-RAY, SOMETHING WRONG!!!'

120     Continue

        izneb = l + 1

130     Continue

      End If
    End If

!   if not first

  Else

    lm = ltot(jp)

    Do k = 1, nb
      If (k==1) x01(1) = x0
!     changed
      If (k>=2) x01(k) = (x0+deltaarr(k))*xnue(1)/xnue(k)
!     #------  do 290 l=izr,lm changed!!!!
      Do l = 1, lm
        xcmf(l, k) = xcmfp(l) + x01(k)
      End Do
    End Do

!   w = weights for flux integral

    w = wp(jp)
    w1 = wp1(jp)

  End If

  lm = ltot(jp)

! defining continuum source function at appropriate cmf frequency
! corresponding to x0. the result is written to sconray(l,1). the
! algorithm exploits the monotonicity of xcmf (decreasing with z)
! note, that sce and xcmfe are defined from blue to red

  If (escat) Then

    xcmfl = xcmfp(lm) + x0
    Do k = 1, nfcmf - 1
      If (xcmfl<=xcmfe(k) .And. xcmfl>=xcmfe(k+1)) Go To 140
    End Do
    Print *, xcmfl
    Print *
    Print *, xcmfe
    Stop ' XCMF(LM) NOT FOUND IN XCMFE'

140 Continue

    q1 = (xcmfl-xcmfe(k+1))/(xcmfe(k)-xcmfe(k+1))
    sconray(lm, 1) = q1*sceray(lm, k) + (1.D0-q1)*sceray(lm, k+1)
    sconray(lm, 2) = 0.D0
    kk = k

    Do l = lm, 1, -1
      xcmfl = xcmfp(l) + x0
      Do k = kk, nfcmf - 1
        If (xcmfl<=xcmfe(k) .And. xcmfl>=xcmfe(k+1)) Go To 150
      End Do

      Do k = 1, kk - 1
        If (xcmfl<=xcmfe(k) .And. xcmfl>=xcmfe(k+1)) Go To 150
      End Do
      Print *, l, ' ', xcmfl
      Print *, xcmfe
      Stop ' XCMF(L) NOT FOUND IN XCMFE'


150   Continue

      q1 = (xcmfl-xcmfe(k+1))/(xcmfe(k)-xcmfe(k+1))
      sconray(l, 1) = q1*sceray(l, k) + (1.D0-q1)*sceray(l, k+1)
      sconray(l, 2) = 0.D0
      kk = k
    End Do

  End If

! location of resonance zone(s) for singlets and overlapping
! multiplets (between comp. 1 and NB)

  If (optiov) Then

    If (first) Then
      izr = 1
      izb = 1
    End If

    k = nb

    xmaxi = max(xmax(1), xmax(nb))

    If (xcmf(1,k)>=-xmaxi) Then
      izr = 1
    Else
      Do l = izr, lm - 1
        If (xcmf(l,k)<=-xmaxi .And. xcmf(l+1,k)>-xmaxi) Exit
      End Do
      izr = l
    End If

    If (xcmf(1,1)>=xmaxi) Then
      izb = 1
    Else If (xcmf(lm,1)<xmaxi) Then
      izb = lm
    Else
      Do l = max(1, izb-1), lm - 1
        If (xcmf(l,1)<xmaxi .And. xcmf(l+1,1)>=xmaxi) Exit
      End Do
      izb = l + 1

      If (izb>lm) Stop ' IZB > LM'
    End If

!   print*,jp,' ',x0,' ',izr,' ',izb,' ',izner,' ',izneb

    Return

  Else

!   location of resonance zone(s) for non-overlapping doublets


    If (first) Then
      izr = 1
      izr2 = 1
      izb1 = 1
      izb = 1
    End If

!   red component(outer one)

    If (xcmf(1,2)>=-xmax(2)) Then
      izr = 1
    Else
      Do l = izr, lm - 1
        If (xcmf(l,2)<=-xmax(2) .And. xcmf(l+1,2)>-xmax(2)) Exit
      End Do
      izr = l
    End If

    If (xcmf(1,2)>=xmax(2)) Then
      izr2 = 1
    Else If (xcmf(lm,2)<xmax(2)) Then
      izr2 = lm
    Else
      Do l = max(1, izr2-1), lm - 1
        If (xcmf(l,2)<xmax(2) .And. xcmf(l+1,2)>=xmax(2)) Exit
      End Do
      izr2 = l + 1

      If (izr2>lm) Stop ' IZR2 > LM'
    End If

!   blue component(inner one)

    If (xcmf(1,1)>=-xmax(1)) Then
      izb1 = 1
    Else
      Do l = izb1, lm - 1
        If (xcmf(l,1)<=-xmax(1) .And. xcmf(l+1,1)>-xmax(1)) Exit
      End Do
      izb1 = l
    End If

    If (xcmf(1,1)>=xmax(1)) Then
      izb = 1
    Else If (xcmf(lm,1)<xmax(1)) Then
      izb = lm
    Else
      Do l = max(1, izb-1), lm - 1
        If (xcmf(l,1)<xmax(1) .And. xcmf(l+1,1)>=xmax(1)) Exit
      End Do
      izb = l + 1

      If (izb>lm) Stop ' IZB > LM'
    End If

!   print*,
!   *jp,' ',x0,' ',izr,' ',izr2,' ',izb1,' ',izb,' ',izner,' ',izneb

    Return

  End If

End Subroutine

!-----------------------------------------------------------------------

Subroutine xgrid(delta, nfobs, xmaxdop, xmaxdop_min, x0, aic, opacon, nd, tem, &
  tgrad, xcred, xcblue, xnue0, vmax, nb, vplus, escat)

! new version, from V7.0 on.

  Use :: nlte_type
  Use :: nlte_dim
  Use :: fund_const, Only: clight
  Use :: formalsol_var, Only: nfcmf, xcmfe, xrest => xnue
  Implicit None


! .. parameters ..
  Integer (i4b), Parameter :: nf = id_nfobs
! ..
! .. scalar arguments ..
  Real (dp) :: delta, tem, tgrad, vmax, vplus, xcblue, xcred, xnue0
  Integer (i4b) :: nb, nd, nfobs
  Logical :: escat
! ..
! .. array arguments ..
  Real (dp) :: aic(nfobs, 2), opacon(nd, 2), x0(nfobs), xmaxdop(nb), &
    xmaxdop_min(nb)
! ..
  Real (dp), Dimension (nb) :: xlam1


! .. local scalars ..
  Real (dp) :: c2, dxs1, dxs2, fmax, fmin, opa, sat, vm, vmm, vmp, xlam, xsat, &
    xx, xx0, dx, fmax1, fmin1
  Integer (i4b) :: i, k, m, m1, nperline
! ..
! .. external functions ..
  Real (dp) :: bplanc
  External :: bplanc
! ..
! .. external subroutines ..
  External :: sort
! ..
! .. intrinsic functions ..
! ..

  If (nf/=nfobs) Stop ' NF NE NFOBS IN XGRID'

  If (nb==1 .And. delta/=0.D0) Stop ' NB=1 AND DELTA NE 0'

! C2 = 1.439426D0/TEM
! changed by JO July 2015
  c2 = 1.4388354967334D0/tem

! definition of satellite components for hei-lines

  xlam = 1.D8/xnue0

  If (abs(xlam-4026.D0)<3.) Then
    Print *
    Print *, 'SATELLITE-COMPONENT AT -0.6 A!!'
    Print *
    sat = -.6D0

  Else If (abs(xlam-4387.D0)<3.) Then
    Print *
    Print *, 'SATELLITE-COMPONENT AT -0.525 A!!'
    Print *
    sat = -0.525D0

  Else If (abs(xlam-4471.D0)<3.) Then
    Print *
    Print *, 'SATELLITE-COMPONENT AT -1.5 A!!'
    Print *
    sat = -1.5D0

  Else If (abs(xlam-4922.D0)<3.) Then
    Print *
    Print *, 'SATELLITE-COMPONENT AT -1.4 A!!'
    Print *
    sat = -1.4D0

  Else
    sat = 0.D0
  End If

  xsat = (xlam/(xlam+sat)-1.D0)*clight/vmax
  If (sat==0.) xsat = 0.D0

! defining maximum observer's frame bandwidth

  vm = 1.D0
  vmp = vm + xmaxdop(1) !           calculated with VDOPMAX
  vmm = vm + xmaxdop(nb) !          calculated with VDOPMAX

  Print *
  Print *, ' VMAX + VD_MAX = ', vmp
  Print *

  fmax = max(vmp, vplus) !          vplus calculated with XMAX
  fmin = max(vmm, vplus) + delta

  fmax1 = fmax
  fmin1 = fmin

  If (escat) Then
!   IN CASE OF VERY BROAD LINES, CHANGE FIRST AND LAST XCMFE-FREQ. POINTS
!   (faked anyway)
    If (fmax1>xcmfe(1)-1.D0) xcmfe(1) = fmax1 + 1.01
    If (fmin1>-xcmfe(nfcmf)-1.D0) xcmfe(nfcmf) = -fmin1 - 1.01

    fmax1 = max(fmax1, xcmfe(1)-1.D0)
    fmin1 = max(fmin1, -xcmfe(nfcmf)-1.D0)
  End If

! so far, equidistant grid, resolved at line center
  nperline = 32 !                   sufficient for 3 lines
  m1 = max(60, nfobs-nb*nperline)

  If (m1==60) Then
    nperline = (nfobs-60)/nb
    If (mod(nperline,2)/=0) nperline = nperline - 1
    m1 = nfobs - nb*nperline
    Print *
    Print *, ' WARNING!!!! WARNING!!!! WARNING!!!!'
    Print *, ' ONLY ', nperline, ' HIGHLY RESOLVED FREQUENCY POINTS PER LINE'
    Print *, ' WARNING!!!! WARNING!!!! WARNING!!!!'
    Print *
  End If
  m1 = m1 - nb

  Do i = 1, nb
    xlam1(i) = 1.D8/xrest(i)
!   XLAM1 in x-units, w.r.t. bluemost component
    xlam1(i) = (xlam/xlam1(i)-1.D0)*clight/vmax ! nu/nu0 corresponds to
!   lam0/lam
  End Do

  If (fmax==fmax1) Then !           either no elscat, or elscat wings narrow
!   equidistant grid with M1 points
    xx = (fmax+fmin)/dble(m1-1)
    Do i = 1, m1
      x0(i) = fmax - (i-1)*xx
    End Do
    m = m1

  Else
    If (.Not. escat) Stop ' ERROR IN ESCAT PHILOSOPHY - SUBR. XGRID'
!   equidistant grid with M1-20 points
    m1 = m1 - 20
    xx = (fmax+fmin)/dble(m1-1)
    Do i = 1, m1
      x0(i) = fmax - (i-1)*xx
    End Do
!   equidistant grid between FMAX1 and FMAX with 10 points
    m = m1
    xx = (fmax1-fmax)/dble(10)
    Do i = 1, 10
      m = m + 1
      x0(m) = fmax1 - (i-1)*xx !    (excl. FMAX)
    End Do
!   equidistant grid between FMIN and FMIN1 with 10 points
    xx = (fmin1-fmin)/dble(10)
    Do i = 1, 10
      m = m + 1
      x0(m) = -fmin - i*xx !        (excl. FMIN)
    End Do
  End If

! higher resolution around line cores, logarithmically distributed
  m = m + 1
  Do k = 1, nb
    dxs1 = alog10(4./0.05)/dble(nperline/2-1) ! from 0.05 to 4 (previously 2)
!   max-doppler units (usually
!   xmaxdop=3 vdop)
    x0(m) = xlam1(k)
    m = m + 1
    Do i = 1, nperline/2
      dx = 0.05*xmaxdop_min(k)*10.**((i-1)*dxs1) ! calculated with
!     XMAXDOP_MIN,
!     to assure resolution of core
      x0(m) = xlam1(k) + dx
!     print*,i,nperline,x0(m),dx,xmaxdop(k),xmax(k)
      m = m + 1
      x0(m) = xlam1(k) - dx
      m = m + 1
    End Do
  End Do

  x0 = -x0
  Call sort(nfobs, x0)
  x0 = -x0

  If (xsat==0.D0) Go To 110

! set in satelite component

  Do k = 1, nfobs - 1
    If (xsat<x0(k) .And. xsat>=x0(k+1)) Go To 100
  End Do

  Stop ' SATELLITE COMPONENT NOT FOUND'

100 Continue

  dxs1 = x0(k) - xsat
  dxs2 = xsat - x0(k+1)

  If (dxs2<dxs1) Then
    x0(k+1) = xsat
  Else
    x0(k) = xsat
  End If

! loop for every observers-frame-frequency

110 Continue

  Do i = 1, nfobs

!   ---- i-core calculated under difussion appr. with limb darkening

!   opa=opacon(nd,1)+(opacon(nd,2)-opacon(nd,1))/(xcblue-
!   *         xcred)*(x0(i)-xcred)

!   changed to be consistent with constant opacity approx. used
!   in elscat

    opa = opacon(nd, 1) + (opacon(nd,2)-opacon(nd,1))/(xcblue-xcred)*(0.D0- &
      xcred)
    xx0 = x0(i)*vmax*xnue0/clight + xnue0

    aic(i, 1) = bplanc(xx0, tem)
    aic(i, 2) = aic(i, 1)*c2*xx0*tgrad/tem/opa/(1.D0-exp(-c2*xx0))

  End Do

  Return

End Subroutine

!-----------------------------------------------------------------------

Subroutine xgrid_old(delta, nfobs, xmaxdop, x0, aic, opacon, nd, tem, tgrad, &
  xcred, xcblue, xnue0, vmax, nb, vplus, escat)

! old version, until V6.1
  Use :: nlte_type
  Use :: nlte_dim
  Use :: fund_const, Only: clight
  Use :: formalsol_var, Only: nfcmf, xcmfe, xrest => xnue
  Implicit None


! .. parameters ..
  Integer (i4b), Parameter :: nf = id_nfobs
! ..
! .. scalar arguments ..
  Real (dp) :: delta, tem, tgrad, vmax, vplus, xcblue, xcred, xnue0
  Integer (i4b) :: nb, nd, nfobs
  Logical :: escat
! ..
! .. array arguments ..
  Real (dp) :: aic(nfobs, 2), opacon(nd, 2), x0(nfobs), xmaxdop(nb)
! ..
! REAL(DP), DIMENSION(NF) :: AUX
  Real (dp), Dimension (nb) :: xlam1


! .. local scalars ..
  Real (dp) :: c2, dxs1, dxs2, fmax, fmin, fxlog, opa, sat, vm, vmm, vmp, &
    xlam, xsat, xx, xx0, yy, dx
  Integer (i4b) :: i, ie, ieb, ier, k, m, m1, mblue, mred, nfobs1, mfirst, &
    nperline
! ..
! .. external functions ..
  Real (dp) :: bplanc
  External :: bplanc
! ..
! .. external subroutines ..
  External :: sort
! ..
! .. intrinsic functions ..
! ..

  If (nf/=nfobs) Stop ' NF NE NFOBS IN XGRID'

  If (nb==1 .And. delta/=0.D0) Stop ' NB=1 AND DELTA NE 0'

! changed by JO July 2015
  c2 = 1.4388354967334D0/tem

! definition of satellite components for hei-lines

  xlam = 1.D8/xnue0

  If (abs(xlam-4026.D0)<3.) Then
    Print *
    Print *, 'SATELLITE-COMPONENT AT -0.6 A!!'
    Print *
    sat = -.6D0

  Else If (abs(xlam-4387.D0)<3.) Then
    Print *
    Print *, 'SATELLITE-COMPONENT AT -0.525 A!!'
    Print *
    sat = -0.525D0

  Else If (abs(xlam-4471.D0)<3.) Then
    Print *
    Print *, 'SATELLITE-COMPONENT AT -1.5 A!!'
    Print *
    sat = -1.5D0

  Else If (abs(xlam-4922.D0)<3.) Then
    Print *
    Print *, 'SATELLITE-COMPONENT AT -1.4 A!!'
    Print *
    sat = -1.4D0

  Else
    sat = 0.D0
  End If

  xsat = (xlam/(xlam+sat)-1.D0)*clight/vmax
  If (sat==0.) xsat = 0.D0

! defining maximum observer's frame bandwidth

  vm = 1.D0
  vmp = vm + xmaxdop(1)
  vmm = vm + xmaxdop(nb)

  Print *
  Print *, ' VMAX + VD_MAX = ', vmp
  Print *

  fmax = max(vmp, vplus)
  fmin = max(vmm, vplus) + delta
  Print *, 'vmp', vmp, vmm
  Print *, 'vplus', vplus


  If (.Not. escat) xcmfe(1) = 0.D0

  If (nb>2) Go To 100

! path for singlets and doublets

  If (.Not. escat .Or. fmax>=xcmfe(1)) Then
    m = nfobs/3

!   ---- blue wing

    m1 = 2*m + 1
    mblue = nfobs - m1
    mblue = mblue/2
    mred = nfobs - (m1+mblue)

    If (vplus<vmp) Then

!     1/3 of points from vmp to .5 and from -.5 to -vmm linear
!     2/3 of points from .5 to -.5 logarithmic for NB=1
!     1/3 of points from .5 to -.5 logarithmic for each component if NB=2


      xx = (fmax-.5D0)/dble(mblue)
      Do i = 1, mblue
        x0(i) = fmax - (i-1)*xx
      End Do
      k = mblue

      fxlog = log10(0.5D0)

      If (nb==1) Then
!       ONE COMPONENT
        xx = (fxlog+4.D0)/dble(m-1)
        Do i = 1, m
          k = k + 1
          x0(k) = 1.D1**(fxlog-dble(i-1)*xx)
        End Do

!       -----one point at x=0.

        k = k + 1
        x0(k) = 0.D0

!       ---- red wing

        fxlog = log10(0.5D0)
        yy = (fxlog+4.D0)/dble(m-1)
        Do i = 1, m
          k = k + 1
          x0(k) = -1.D1**(-4.D0+dble(i-1)*yy)
        End Do

        xx = (fmin-.5D0)/dble(mred)
        Do i = 1, mred
          k = k + 1
          x0(k) = -0.5D0 - i*xx
        End Do

!       -----
      Else
!       TWO COMPONENTS
        If (nb/=2) Stop ' SOMETHING WRONG WITH NB(1)'
        mfirst = m/2
        xx = (fxlog+4.D0)/dble(mfirst-1)
        Do i = 1, mfirst
          k = k + 1
          x0(k) = 1.D1**(fxlog-dble(i-1)*xx)
          k = k + 1
          x0(k) = 1.D1**(fxlog-dble(i-1)*xx) - delta
        End Do

!       -----one point at x=0.

        k = k + 1
        x0(k) = 0.D0
        k = k + 1
        x0(k) = -delta

!       ---- red wing

        fxlog = log10(0.5D0)
        yy = (fxlog+4.D0)/dble(mfirst-1)
        Do i = 1, mfirst
          k = k + 1
          x0(k) = -1.D1**(-4.D0+dble(i-1)*yy)
          k = k + 1
          x0(k) = -1.D1**(-4.D0+dble(i-1)*yy) - delta
        End Do

        If (mred+k/=nfobs-1) Stop ' SOMETHING WRONG IN MRED PHILOSOPHY'
        mred = mred + 1

        xx = (fmin-.5D0)/dble(mred)
        Do i = 1, mred
          k = k + 1
          x0(k) = -0.5D0 - i*xx
        End Do

        Do i = 1, nfobs
          x0(i) = -x0(i)
        End Do

        Call sort(nfobs, x0)

        Do i = 1, nfobs
          x0(i) = -x0(i)
        End Do

        If (k/=nfobs) Stop 'ERROR IN OBSERV. FREQUENCY GRID, K.NE.NF(1)'

      End If !                      ONE OR TWO COMPONENT TREATMENT
    Else

!     2/3 of points from vplus to -vplus logarithmic for NB=1
!     1/3 of points from vplus to -vplus logarithmic for each component if
!     NB=2
!     1/3 of points from vmp to .5 and from -.5 to -vmm
!     linear additionally

      fxlog = log10(fmax)
      If (nb==1) Then
!       ONE COMPONENT
        xx = (fxlog+4.D0)/dble(m-1)
        k = 0
        Do i = 1, m
          k = k + 1
          x0(k) = 1.D1**(fxlog-dble(i-1)*xx)
        End Do

!       -----one point at x=0.

        k = k + 1
        x0(k) = 0.D0

!       ---- red wing

        fxlog = log10(fmin)
        yy = (fxlog+4.D0)/dble(m-1)
        Do i = 1, m
          k = k + 1
          x0(k) = -1.D1**(-4.D0+dble(i-1)*yy)
        End Do

      Else
!       TWO COMPONENTS
        If (nb/=2) Stop ' SOMETHING WRONG WITH NB(2)'
        mfirst = m/2
        xx = (fxlog+4.D0)/dble(mfirst-1)
        k = 0
        Do i = 1, mfirst
          k = k + 1
          x0(k) = 1.D1**(fxlog-dble(i-1)*xx)
          k = k + 1
          x0(k) = 1.D1**(fxlog-dble(i-1)*xx) - delta
        End Do

!       -----one point at x=0.

        k = k + 1
        x0(k) = 0.D0
        k = k + 1
        x0(k) = -delta

!       ---- red wing

        fxlog = log10(fmin)
        yy = (fxlog+4.D0)/dble(mfirst-1)
        Do i = 1, mfirst
          k = k + 1
          x0(k) = -1.D1**(-4.D0+dble(i-1)*yy)
          k = k + 1
          x0(k) = -1.D1**(-4.D0+dble(i-1)*yy) - delta
        End Do

        If (k+mblue+mred/=nfobs-1) Stop &
          ' SOMETHING WRONG IN MBLUE/MRED PHILOSOPHY'
        mred = mred + 1

      End If !                      ONE OR TWO COMPONENT TREATMENT

!     additional points

      xx = (vmp-.5D0)/dble(mblue-1)
      Do i = 1, mblue
        k = k + 1
        x0(k) = vmp - (i-1)*xx
      End Do

      xx = (vmm-.5D0)/dble(mred-1)
      Do i = 1, mred
        k = k + 1
        x0(k) = -0.5D0 - (i-1)*xx
      End Do

      If (k/=nfobs) Stop 'ERROR IN OBSERV. FREQUENCY GRID, K.NE.NF(2)'

      Do i = 1, nfobs
        x0(i) = -x0(i)
      End Do

      Call sort(nfobs, x0)

      Do i = 1, nfobs
        x0(i) = -x0(i)
      End Do
    End If

  Else !                            ELSCAT TREATMENT

    If (xcmfe(1)-xcmfe(2)<1.) Stop ' DELTA X(1) < 1'

    ieb = 0

    Do k = 1, nfcmf
      If (xcmfe(k)>fmax) ieb = ieb + 1
    End Do

    ier = 0

    Do k = 1, nfcmf
      If (xcmfe(k)<-fmin) ier = ier + 1
    End Do

    ie = ieb + ier

!   -----now, ie is the number of freq. points for el. scattering wings

    nfobs1 = nfobs - ie

    m = nfobs1/3

!   ---- blue wing

    m1 = 2*m + 1
    mblue = nfobs1 - m1
    mblue = mblue/2
    mred = nfobs1 - (m1+mblue)

    If (vplus<vmp) Then

!     1/3 of points from vmp to .5 and from -.5 to -vmm linear
!     2/3 of points from .5 to -.5 logarithmic

      xx = (fmax-.5D0)/dble(mblue)
      Do i = 1, mblue
        x0(i) = fmax - (i-1)*xx
      End Do

      k = mblue
      fxlog = log10(0.5D0)
      xx = (fxlog+4.D0)/dble(m-1)

      Do i = 1, m
        k = k + 1
        x0(k) = 1.D1**(fxlog-dble(i-1)*xx)
      End Do

!     -----one point at x=0.

      k = k + 1
      x0(k) = 0.D0

!     ---- red wing

      fxlog = log10(0.5D0)
      yy = (fxlog+4.D0)/dble(m-1)
      Do i = 1, m
        k = k + 1
        x0(k) = -1.D1**(-4.D0+dble(i-1)*yy)
      End Do

      xx = (fmin-.5D0)/dble(mred)
      Do i = 1, mred
        k = k + 1
        x0(k) = -0.5D0 - i*xx
      End Do

    Else

!     2/3 of points from vplus to -vplus logarithmic
!     1/3 of points from vmp to .5 and from -.5 to -vmm
!     linear additionally

      fxlog = log10(fmax)
      xx = (fxlog+4.D0)/dble(m-1)
      k = 0
      Do i = 1, m
        k = k + 1
        x0(k) = 1.D1**(fxlog-dble(i-1)*xx)
      End Do

!     -----one point at x=0.

      k = k + 1
      x0(k) = 0.D0

!     ---- red wing

      fxlog = log10(fmin)
      yy = (fxlog+4.D0)/dble(m-1)
      Do i = 1, m
        k = k + 1
        x0(k) = -1.D1**(-4.D0+dble(i-1)*yy)
      End Do

!     additional points

      xx = (vmp-.5D0)/dble(mblue-1)
      Do i = 1, mblue
        k = k + 1
        x0(k) = vmp - (i-1)*xx
      End Do

      xx = (vmm-.5D0)/dble(mred-1)
      Do i = 1, mred
        k = k + 1
        x0(k) = -0.5D0 - (i-1)*xx
      End Do

    End If

!   additional cmf-points

    k = k + 1
    x0(k) = xcmfe(1) - 1.D0
    Do i = 2, ieb
      k = k + 1
      x0(k) = xcmfe(i)
    End Do

    Do i = 1, ier - 1
      k = k + 1
      x0(k) = xcmfe(nfcmf-ier+i)
    End Do

    k = k + 1
    x0(k) = xcmfe(nfcmf) + 1.D0

    Do i = 1, nfobs
      x0(i) = -x0(i)
    End Do

    Call sort(nfobs, x0)

    Do i = 1, nfobs
      x0(i) = -x0(i)
    End Do

    If (k/=nfobs) Stop 'ERROR IN OBSERV. FREQUENCY GRID, K.NE.NF(3)'

  End If


  Go To 110

! path for three and more components
! so far, equidistant grid, resolved at line center
100 nperline = 32
  m1 = max(60, nfobs-nb*nperline)

  If (m1==60) Then
    nperline = (nfobs-60)/nb
    If (mod(nperline,2)/=0) nperline = nperline - 1
    m1 = nfobs - nb*nperline
    Print *
    Print *, ' WARNING!!!! WARNING!!!! WARNING!!!!'
    Print *, ' ONLY ', nperline, ' HIGHLY RESOLVED FREQUENCY POINTS PER LINE'
    Print *, ' WARNING!!!! WARNING!!!! WARNING!!!!'
    Print *
  End If
  m1 = m1 - nb

  Do i = 1, nb
    xlam1(i) = 1.D8/xrest(i)
    xlam1(i) = (xlam/xlam1(i)-1.D0)*clight/vmax ! nu/nu0 corresponds to
!   lam0/lam
  End Do

! equidistant grid with M1 points
  xx = (fmax+fmin)/dble(m1-1)
  Do i = 1, m1
    x0(i) = fmax - (i-1)*xx
  End Do

! higher resolution around line cores, logarithmically distributed
  m = m1 + 1
  Do k = 1, nb
    dxs1 = alog10(3./0.05)/(nperline/2) ! FROM 0.05 TO 3 DOPPLERUNITS
    x0(m) = xlam1(k)
    m = m + 1
    Do i = 1, nperline/2
      dx = 0.05*xmaxdop(k)*10.**(i*dxs1)
      x0(m) = xlam1(k) + dx
      m = m + 1
      x0(m) = xlam1(k) - dx
      m = m + 1
    End Do
  End Do

  x0 = -x0
  Call sort(nfobs, x0)
  x0 = -x0

110 If (xsat==0.D0) Go To 130

! set in satelite component

  Do k = 1, nfobs - 1
    If (xsat<x0(k) .And. xsat>=x0(k+1)) Go To 120
  End Do

  Stop ' SATELLITE COMPONENT NOT FOUND'

120 Continue

  dxs1 = x0(k) - xsat
  dxs2 = xsat - x0(k+1)

  If (dxs2<dxs1) Then
    x0(k+1) = xsat
  Else
    x0(k) = xsat
  End If

! loop for every observers-frame-frequency

130 Continue

  Do i = 1, nfobs

!   ---- i-core calculated under difussion appr. with limb darkening

!   opa=opacon(nd,1)+(opacon(nd,2)-opacon(nd,1))/(xcblue-
!   *         xcred)*(x0(i)-xcred)

!   changed to be consistent with constant opacity approx. used
!   in elscat

    opa = opacon(nd, 1) + (opacon(nd,2)-opacon(nd,1))/(xcblue-xcred)*(0.D0- &
      xcred)
    xx0 = x0(i)*vmax*xnue0/clight + xnue0

    aic(i, 1) = bplanc(xx0, tem)
    aic(i, 2) = aic(i, 1)*c2*xx0*tgrad/tem/opa/(1.D0-exp(-c2*xx0))

  End Do

  Return

End Subroutine

!***********************************************************************

!subroutines: broadening functions

!***********************************************************************

Function perfil1(i, xi, tem, xne, vturbl, flag)

! no action concerning clumping required, since we solve inside clumps,
! and n_e is the increased value

! new version, allowing for depth dependent vturb (from header)
  Use :: nlte_type
  Use :: nlte_dim
  Use :: fund_const, Only: akb, xmh => amh, clight, rapi => sqpi
  Use :: formalsol_var, Only: gammal, gammac, vmax, xnemin, weig, xnue, &
    nstark, nws, nts, nes, qhalf, dws, ts, es, ps

  Implicit None

! ----it calculates doppler and stark profiles (if neccesary) for the 'i'
! ----component, at frequency 'x', with temperature and electron density
! ----'tem' and 'xne' resp.
! In this version, vturb = vturbmin = const has been folded (PREFORMAL) into
! all
! tabulated profiles (Stark, Iso and Griem)

! ..
! .. scalar arguments ..
  Real (dp) :: tem, xi, xne, perfil1, vturbl
  Integer (i4b) :: i
  Logical :: flag
! ..
! .. local scalars ..
  Real (dp) :: a, aktm, avoigt, deldop, dw, elog, pf, pff, pfl, pl, plf, pll, &
    tlog, vth, wef, wtf, wwf, xd, xlam, xlamb0
  Integer (i4b) :: ifi, j, jeb, jtb, jwb, l1, l2, l3, l4
! ..
! .. external functions ..
  Real (dp) :: voigt
  External :: voigt
! ..
! .. intrinsic functions ..
! ..

  flag = .False.
  xlamb0 = 1.D0/xnue(i)

  If (nstark(i)==0 .Or. xne<xnemin) Then

!   --- pure doppler profile (normalized in frequencies) with depth-dependent
!   vturb

    aktm = 2.D0*akb*tem/xmh/weig(i) + vturbl*vturbl
    a = xi*xi*vmax*vmax/aktm

    If (a>25) Then
      perfil1 = 0.D0
      flag = .True.
    Else
      perfil1 = exp(-a)/sqrt(aktm)*xlamb0/rapi
    End If

  Else If (nstark(i)==2 .Or. nstark(i)==3) Then
!   ---- with VTURBL

!   ------- voigt profile

    vth = sqrt(2.D0*akb*tem/xmh/weig(i)+vturbl*vturbl)
    deldop = vth*xnue(i)

!   --- note: gammal, gammac include 1/(4pi)

    avoigt = (gammal(i)+xne*gammac(i))/deldop
    xd = xi*vmax/vth
    perfil1 = voigt(avoigt, xd)*xlamb0/(vth*rapi)
!   p1=exp(-xd*xd)/vth*xlamb0/rapi
!   print*,log10(xne),' ',xd,' ',avoigt,' ',perfil,' ',p1
    Return

  Else If (nstark(i)==1) Then
!   ---- with VTURB=VTURBMIN (implicit included in tables)

!   --- stark profile (doppler broadening -- temp+Stark -- is included)

    xlam = xlamb0/(1.D0+xi*vmax/clight)

!   --- weights for wavelength

    dw = (xlam-xlamb0)*1.D8
    If (qhalf(i)) dw = abs(dw)

!   --- out of freq. range in table (assuming then phi = 0.)

    If (dw<dws(1,i) .Or. dw>dws(nws(i),i)) Then
      flag = .True.
      perfil1 = 0.D0
      Return
    End If

    Do j = 2, nws(i)
      ifi = j
      If (dws(j,i)>=dw) Go To 100
    End Do

    Stop ' ERROR IN TABLE FREQ. GRID'

100 Continue

    jwb = ifi - 1
    wwf = (dw-dws(ifi-1,i))/(dws(ifi,i)-dws(ifi-1,i))

!   --- weights for t

    tlog = log10(tem)
    tlog = max(ts(1,i), min(ts(nts(i),i),tlog))
    Do j = 2, nts(i)
      ifi = j
      If (ts(j,i)>=tlog) Exit
    End Do

    jtb = ifi - 1
    wtf = (tlog-ts(ifi-1,i))/(ts(ifi,i)-ts(ifi-1,i))

    elog = log10(xne)
    If (elog<es(1,i)) Then
      If (xnemin/=0.) Then
        Print *, elog, es(1, i)
        Stop ' ERROR IN XNEMIN OR STARK PROFILE TABLES'
      End If

!     perfil called by CALCXMAX (xnemin not calculated yet) and
!     very thin atmosphere (SNe)

      elog = es(1, i) !             approximate treatment (just to calculate
!     xmax)
    End If

    Do j = 2, nes(i)
      ifi = j
      If (es(j,i)>=elog) Exit
    End Do

    jeb = ifi - 1
    wef = (elog-es(ifi-1,i))/(es(ifi,i)-es(ifi-1,i))

!   --- interpolation

    l1 = jwb + nws(i)*(jtb+nts(i)*jeb)
    pff = wwf*(ps(l1+1,i)-ps(l1,i)) + ps(l1, i)
    l2 = l1 - nws(i)
    plf = wwf*(ps(l2+1,i)-ps(l2,i)) + ps(l2, i)
    l3 = l1 - nws(i)*nts(i)
    pfl = wwf*(ps(l3+1,i)-ps(l3,i)) + ps(l3, i)
    l4 = l3 - nws(i)
    pll = wwf*(ps(l4+1,i)-ps(l4,i)) + ps(l4, i)
    pf = wtf*(pff-plf) + plf
    pl = wtf*(pfl-pll) + pll

    perfil1 = 10.D0**(wef*(pf-pl)+pl)

  Else

    Stop 'ERROR IN NSTARK  -  PERFIL1'

  End If

  Return
End Function

!-----------------------------------------------------------------------

Subroutine rconv(x, y0, y1, n, vrot)

  Use :: nlte_type
  Use :: nlte_dim
  Implicit None


! -----unrotated profile y0
! -----  rotated profile y1

! -----rotation convolution result returned in y0 y1 is scratch
! -----limb-darkening coefficient 1.5 used
! .. scalar arguments ..
  Real (dp) :: vrot
  Integer (i4b) :: n
! ..
! .. array arguments ..
  Real (dp) :: x(n), y0(n), y1(n)
! ..
! .. local scalars ..
  Real (dp) :: con, div, pb, pbx, pf, pfx, pi1, qad, summ, x0, xb, xf
  Integer (i4b) :: i, j, m0, n0
! ..
! .. intrinsic functions ..
! ..

  If (vrot<0.1) Then
    Do i = 1, n
      y1(i) = y0(i)
    End Do
    Return
  End If

  con = 2.997925D5/vrot
  pi1 = 0.31830D0
  n0 = 2

  Do i = 1, n
    summ = 0.D0
    qad = 0.D0
    m0 = n0
    x0 = x(i)
    div = con/x0
    Do j = m0, n
      xb = (x(j-1)-x0)*div
      xf = (x(j)-x0)*div
      If (xb>1.0) Go To 110
      If (xf<-1.0) Go To 100
      xb = max(-1.0D0, xb)
      xf = min(1.0D0, xf)
      pbx = 1.0D0 - xb*xb
      pfx = 1.0D0 - xf*xf
      pb = pi1*sqrt(pbx) + .375D0*pbx
      pf = pi1*sqrt(pfx) + .375D0*pfx
      summ = summ + (pb*y0(j-1)+pf*y0(j))*(xf-xb)
      qad = qad + (pb+pf)*(xf-xb)
      Exit
100   Continue
      n0 = j
    End Do

110 Continue

    y1(i) = summ/qad
  End Do

  Return
End Subroutine

!-----------------------------------------------------------------------

Function voigt(aa, vv)

  Use :: nlte_type
  Use :: nlte_dim
  Implicit None

! -----computes unnormalized voigt function--i.e. normalized to root pi
! .. scalar arguments ..
  Real (dp) :: aa, vv, voigt
! ..
! .. local scalars ..
  Real (dp) :: a, a2, a4, a6, c1, c2, v, v2, v4, v6, w, x, z, z2
  Integer (i4b) :: i, j
! ..
! .. local arrays ..
  Real (dp) :: h(25)
! ..
! .. external functions ..
  Real (dp) :: dawson
  External :: dawson
! ..
! .. intrinsic functions ..
! ..
! .. save statement ..
  Save
! ..
! .. data statements ..
  Data c1, c2/1.128379167095512D0, 5.64189583547756D-1/
! ..
  v = abs(vv)
  a = aa
  v2 = v*v
  a2 = a*a
  z = a2 + v2
  If (a<=0.5) Go To 110
  If (z<10.) Go To 120
! -----asymptotic expansion for large modulus
100 Continue
  z2 = z*z
  v4 = v2*v2
  v6 = v4*v2
  a4 = a2*a2
  a6 = a4*a2
  voigt = c2*a*(1.D0+((1.875D0*(7.D0*v6-35.D0*a2*v4+21.D0*a4*v2-a6)/z2+ &
    0.75D0*(5.D0*v4-10.D0*a2*v2+a4))/z2+1.5D0*v2-0.5D0*a2)/z2)/z
  Return
! -----harris expansion
110 Continue
  If (v>5.) Go To 100
  w = dawson(v)
  h(1) = exp(-v2)
  h(2) = -c1*(1.D0-2.D0*v*w)
  h(3) = (1.D0-2.D0*v2)*h(1)
  h(4) = -c1*(2.D0*(1.D0-v2)/3.D0-2.D0*v*w*(1.D0-2.D0*v2/3.D0))
! -----higher terms by recursion
  Do i = 5, 11
    x = i - 1
    h(i) = (2.D0*(2.D0*x-3.D0-2.D0*v2)*h(i-2)-4.D0*h(i-4))/(x*(x-1.D0))
  End Do
  voigt = h(11)
  Do i = 1, 10
    j = 11 - i
    voigt = h(j) + a*voigt
  End Do
  Return
! -----gronwall expansion
120 Continue
  x = 1.D0/(1.D0+3.275911D-1*a)
  h(1) = ((((1.061405429D0*x-1.453152027D0)*x+1.421413741D0)*x-2.84496736D-1)* &
    x+2.54829592D-1)*x
  Do i = 2, 25
    x = i - 1
    h(i) = 2.D0*a*(c2-a*h(i-1))/(2.D0*x-1.D0)
  End Do
  voigt = 0.D0
  Do i = 1, 24
    j = 26 - i
    x = j - 1
    voigt = (voigt+h(j))*v2/x
  End Do
  voigt = exp(-v2)*(voigt+h(1))

  Return
End Function

!-----------------------------------------------------------------------

Function dawson(xx)

  Use :: nlte_type
  Use :: nlte_dim
  Implicit None

! -----dawson*s integral using anl algorithm, math comp,1970,171
! .. scalar arguments ..
  Real (dp) :: xx, dawson
! ..
! .. local scalars ..
  Real (dp) :: down, u, up, x
! ..
! .. save statement ..
  Save
! ..
  x = xx
  u = x*x
  If (x<5.0) Go To 100
! -----x greater than 5
  dawson = ((5.0000000167450D-1+7.4999919056701D-1/(-2.5001711668562D0+ &
    u-2.4878765880441D0/(-4.6731202214124D0+u-4.1254406560831D0/( &
    -1.1195216423662D1+u))))/u+1.0D0)/(2.0D0*x)
  Return
! -----x on (3.5,5.0)
100 Continue
  If (x<3.5) Go To 110
  dawson = (5.00001538408193D-1+2.49811162845499D-1/(-1.53672069271915D0+u- &
    6.53419359860764D-1/(-1.77068693717670D1+u+2.04866410976332D2/( &
    7.49584016278357D0+u-2.298758419286D0/(4.02187490205698D1+u+ &
    2.53388006963558D3/(-5.9391591850032D1+u))))))/x
  Return
! -----x on (2.5,3.5)
110 Continue
  If (x<2.5) Go To 120
  dawson = (5.0140106611704D-1+1.8897553014354D-1/(-7.4499050579364D0+u+ &
    7.0204980729194D1/(7.5077816490106D0+u+4.1821806337830D1/( &
    -2.6629001073842D1+u+3.7343084728334D1/(3.0984087863402D1+u+ &
    1.2599323546764D3/(-4.0847391212716D1+u))))))/x
  Return
! -----x less than 2.5
120 Continue
  up = (((((u*2.0846835103886D-2-8.5410681195954D-1)*u+5.4616122556699D1)*u- &
    4.3501160207595D2)*u+9.6696398191665D3)*u-2.9179464300780D4)*u + &
    2.3156975201341D5
  down = (((((u+2.9391995612556D1)*u+4.668490654511D2)*u+ &
    4.7447098440662D3)*u+3.1384620138163D4)*u+1.2520037031851D5)*u + &
    2.3156975201425D5
  dawson = x*(up/down)

  Return
End Function

!***********************************************************************

!subroutines: miscellaneous

!***********************************************************************

Function bplanc(x, t)

  Use :: nlte_type
  Use :: nlte_dim
  Implicit None

! .. scalar arguments ..
  Real (dp) :: t, x, bplanc
! ..
  Real (dp), Parameter :: c1 = 3.972970127D-16, c21 = 1.4388354967334D0
! .. local scalars ..
  Real (dp) :: c2
! ..
! .. intrinsic functions ..
! ..
! C1 = 3.9728D-16
! C2 = 1.439426D0/T
  c2 = c21/t
! changed by JO July 2015
  bplanc = c1*(x**3)/(exp(c2*x)-1.D0)

  Return
End Function

!-----------------------------------------------------------------------

Function expuno(dt)

  Use :: nlte_type
  Use :: nlte_dim
  Implicit None

! ---- calculation of (1-exp(-dt)), expanding the exp for small dt

! .. scalar arguments ..
  Real (dp) :: dt, expuno
! ..
! .. intrinsic functions ..
! ..

  If (abs(dt)<1.D-8) Then
    expuno = dt - 5.D-1*dt*dt
  Else
    expuno = 1.D0 - exp(-dt)
  End If

  Return
End Function

!-----------------------------------------------------------------------

Subroutine sort(n, ra)

  Use :: nlte_type
  Use :: nlte_dim
  Implicit None

! .. scalar arguments ..
  Integer (i4b) :: n
! ..
! .. array arguments ..
  Real (dp) :: ra(n)
! ..
! .. local scalars ..
  Real (dp) :: rra
  Integer (i4b) :: i, ir, j, l
! ..

  l = n/2 + 1
  ir = n
100 Continue

  If (l>1) Then
    l = l - 1
    rra = ra(l)
  Else
    rra = ra(ir)
    ra(ir) = ra(1)
    ir = ir - 1
    If (ir==1) Then
      ra(1) = rra
      Return
    End If
  End If

  i = l
  j = l + l

110 Continue

  If (j<=ir) Then
    If (j<ir) Then
      If (ra(j)<ra(j+1)) j = j + 1
    End If

    If (rra<ra(j)) Then
      ra(i) = ra(j)
      i = j
      j = j + j
    Else
      j = ir + 1
    End If

    Go To 110

  End If

  ra(i) = rra

  Go To 100

End Subroutine

!***********************************************************************

!subroutines: algebraic ones

!***********************************************************************

Subroutine gmalu(ga, u, lmax)

  Use :: nlte_type
  Use :: nlte_dim
  Implicit None

! ***  algebraic routine called from cmfray
! result overwrites ga


! .. scalar arguments ..
  Integer (i4b) :: lmax
! ..
! .. array arguments ..
  Real (dp) :: ga(lmax), u(lmax)
! ..
! .. local scalars ..
  Real (dp) :: a, b
  Integer (i4b) :: l, lz
! ..

  lz = lmax - 1
  b = u(1)

  Do l = 1, lz
    a = b
    b = u(l+1)
    ga(l) = ga(l)*(a-b)
  End Do

  Return
End Subroutine

!-----------------------------------------------------------------------

Subroutine inv(n, ndim, a)

  Use :: nlte_type
  Use :: nlte_dim
  Implicit None

! scaled matrix inversion

! .. scalar arguments ..
  Integer (i4b) :: n, ndim
! ..
! .. array arguments ..
  Real (dp) :: a(ndim, ndim)
! ..
! .. local scalars ..
  Real (dp) :: qi
  Integer (i4b) :: i, ii
! ..
! .. local arrays ..
  Real (dp) :: q(200)
! ..
! .. external subroutines ..
  External :: matinv
! ..
! .. intrinsic functions ..
! ..

  Do i = 1, n
    qi = 0.D0

    Do ii = 1, n
      qi = max(qi, abs(a(ii,i)))
    End Do

    q(i) = qi

    Do ii = 1, n
      a(ii, i) = a(ii, i)/qi
    End Do

  End Do

  Call matinv(a, n, ndim)

  Do ii = 1, n
    qi = q(ii)
    Do i = 1, n
      a(ii, i) = a(ii, i)/qi
    End Do
  End Do

  Return
End Subroutine

!-----------------------------------------------------------------------

Subroutine matinv(a, n, no)

  Use :: nlte_type
  Use :: nlte_dim
  Implicit None


! -----matinv executes matrix inversion by lu decomposition
! -----inversion is accomplished in place,
! -----and the original matrix is replaced by its inverse
! -----note that n must be smaller than no

! .. scalar arguments ..
  Integer (i4b) :: n, no
! ..
! .. array arguments ..
  Real (dp) :: a(no, no)
! ..
! .. local scalars ..
  Real (dp) :: div, summ
  Integer (i4b) :: i, ii, im1, ip1, j, jj, jm1, jp1, k, ko, l
! ..
! .. intrinsic functions ..
! ..

iloop1: Do i = 2, n

    im1 = i - 1
    Do j = 1, im1
      jm1 = j - 1
      If (a(j,j)==0.0D0) a(j, j) = 1.0D-20
      div = a(j, j)
      summ = 0.0D+00
      If (jm1<1) Go To 100

      Do l = 1, jm1
        summ = summ + a(i, l)*a(l, j)
      End Do

100   Continue
      a(i, j) = (a(i,j)-summ)/div
    End Do

    Do j = i, n
      summ = 0.0D+00
      Do l = 1, im1
        summ = summ + a(i, l)*a(l, j)
      End Do
      a(i, j) = a(i, j) - summ
    End Do
  End Do iloop1

iiloop1: Do ii = 2, n
    i = n + 2 - ii
    im1 = i - 1
    If (im1<1) Cycle
    Do jj = 1, im1
      j = i - jj
      jp1 = j + 1
      summ = 0.0D+00
      If (jp1>im1) Go To 110
      Do k = jp1, im1
        summ = summ + a(i, k)*a(k, j)
      End Do
110   Continue
      a(i, j) = -a(i, j) - summ
    End Do
  End Do iiloop1

iiloop2: Do ii = 1, n
    i = n + 1 - ii
    If (a(i,i)==0.0D0) a(i, i) = 1.0D-20
    div = a(i, i)
    ip1 = i + 1
    If (ip1>n) Go To 120
    Do jj = ip1, n
      j = n + ip1 - jj
      summ = 0.0D+00
      Do k = ip1, j
        summ = summ + a(i, k)*a(k, j)
      End Do
      a(i, j) = -summ/div
    End Do
120 Continue
    a(i, i) = 1.D0/a(i, i)
  End Do iiloop2

iloop2: Do i = 1, n
    Do j = 1, n
      ko = max0(i, j)
      If (ko==j) Go To 140
      summ = 0.0D+00

130   Continue

      Do k = ko, n
        summ = summ + a(i, k)*a(k, j)
      End Do
      Go To 150

140   Continue

      summ = a(i, ko)
      If (ko==n) Go To 150
      ko = ko + 1
      Go To 130

150   Continue
      a(i, j) = summ

    End Do
  End Do iloop2

  Return
End Subroutine

!-----------------------------------------------------------------------

Subroutine invtri(a, b, c, q, n)

  Use :: nlte_type
  Use :: nlte_dim
  Implicit None

! solution of tridiagonal equation system with
! tridiag. matrix  -a(l),b(l),-c(l)  ,
! right side: vector q(l),solution overwrites q


! .. scalar arguments ..
  Integer (i4b) :: n
! ..
! .. array arguments ..
  Real (dp) :: a(n), b(n), c(n), q(n)
! ..
! .. local scalars ..
  Real (dp) :: h
  Integer (i4b) :: i
! ..

  c(1) = c(1)/b(1)
  q(1) = q(1)/b(1)

  If (n==1) Return
  If (n==2) Go To 100

  Do i = 2, n - 1
    h = b(i) - a(i)*c(i-1)
    c(i) = c(i)/h
    q(i) = (q(i)+q(i-1)*a(i))/h
  End Do

100 Continue

  q(n) = (q(n)+q(n-1)*a(n))/(b(n)-c(n-1)*a(n))

  Do i = 1, n - 1
    q(n-i) = q(n-i) + c(n-i)*q(n-i+1)
  End Do

  Return
End Subroutine

!-----------------------------------------------------------------------

Subroutine invtri3(ta, tb, tc, a, k, n)

  Use :: nlte_type
  Use :: nlte_dim
  Implicit None


! .. scalar arguments ..
  Integer (i4b) :: k, n
! ..
! .. array arguments ..
  Real (dp) :: a(id_ndept, id_ndept), ta(id_ndept), tb(id_ndept), tc(id_ndept)
! ..
! .. local scalars ..
  Integer (i4b) :: i, j, l
! ..
! .. local arrays ..
  Real (dp) :: d(0:100), dia(100), e(101)
! ..

  If (n>100) Stop 'TOO MANY GRID POINTS IN DIAG'

  d(0) = 0.D0

  Do l = 1, k
    d(l) = -tc(l)/(tb(l)+ta(l)*d(l-1))
  End Do
  e(k+1) = 0.D0

  Do l = k, 1, -1
    e(l) = -ta(l)/(tb(l)+tc(l)*e(l+1))
  End Do

  Do l = 1, k
    dia(l) = 1.D0/((1.D0-d(l)*e(l+1))*(tb(l)+ta(l)*d(l-1)))
    a(l, l) = dia(l)
  End Do

  Do j = 2, k
    Do i = j - 1, 1, -1
      a(i, j) = a(i+1, j)*d(i)
    End Do
  End Do

  Do j = 1, k - 1
    Do i = j + 1, k
      a(i, j) = a(i-1, j)*e(i)
    End Do
  End Do

  Return
End Subroutine

!-----------------------------------------------------------------------

Subroutine madd(a, b, n, ndim)

  Use :: nlte_type
  Use :: nlte_dim
  Implicit None

! a(i,j)=a(i,j)+b(i,j)


! .. scalar arguments ..
  Integer (i4b) :: n, ndim
! ..
! .. array arguments ..
  Real (dp) :: a(ndim, ndim), b(ndim, ndim)
! ..

  a(1:n, 1:n) = a(1:n, 1:n) + b(1:n, 1:n)

  Return
End Subroutine

!-----------------------------------------------------------------------

Subroutine mdmv(a, b, n, ndim)

  Use :: nlte_type
  Use :: nlte_dim
  Implicit None

! matrix a(voll)=: matrix b(diag) * matrix a(voll)


! .. scalar arguments ..
  Integer (i4b) :: n, ndim
! ..
! .. array arguments ..
  Real (dp) :: a(ndim, ndim), b(ndim)
! ..
! .. local scalars ..
  Integer (i4b) :: i, j
! ..

  Do i = 1, n
    Do j = 1, n
      a(i, j) = b(i)*a(i, j)
    End Do
  End Do

  Return
End Subroutine

!-----------------------------------------------------------------------

Subroutine mdv(a, v, n)

  Use :: nlte_type
  Use :: nlte_dim
  Implicit None

! matrix(diagonal) * vector v
! result overwrites v


! .. scalar arguments ..
  Integer (i4b) :: n
! ..
! .. array arguments ..
  Real (dp) :: a(n), v(n)
! ..

  v(1:n) = a(1:n)*v(1:n)

  Return
End Subroutine

!-----------------------------------------------------------------------

Subroutine mvmd(a, b, n, ndim)

  Use :: nlte_type
  Use :: nlte_dim
  Implicit None

! matrix a(voll)=: matrix a(voll) * matrix b(diag)


! .. scalar arguments ..
  Integer (i4b) :: n, ndim
! ..
! .. array arguments ..
  Real (dp) :: a(ndim, ndim), b(ndim)
! ..
! .. local scalars ..
  Integer (i4b) :: i, j
! ..

  Do j = 1, n
    Do i = 1, n
      a(i, j) = a(i, j)*b(j)
    End Do
  End Do

  Return
End Subroutine

!-----------------------------------------------------------------------

Subroutine mvv(wx, b, w, jmax, jmm, jp)

  Use :: nlte_type
  Use :: nlte_dim
  Implicit None

! matrix(voll) b *vektor w= vektor wx
! format: wx(jmax)=b(jmax,jmm)*w(jmm)

! .. scalar arguments ..
  Integer (i4b) :: jmax, jmm, jp
! ..
! .. array arguments ..
  Real (dp) :: b(jp, jp), w(jp), wx(jp)
! ..
! .. local scalars ..
  Real (dp) :: wxi
  Integer (i4b) :: i, k
! ..

  Do i = 1, jmax
    wxi = 0.D0
    Do k = 1, jmm
      wxi = wxi + b(i, k)*w(k)
    End Do
    wx(i) = wxi
  End Do

  Return
End Subroutine

!-----------------------------------------------------------------------

Subroutine vadd(a, b, n)

  Use :: nlte_type
  Use :: nlte_dim
  Implicit None

! vector addition a=a+b


! .. scalar arguments ..
  Integer (i4b) :: n
! ..
! .. array arguments ..
  Real (dp) :: a(n), b(n)
! ..

  a(1:n) = a(1:n) + b(1:n)

  Return
End Subroutine

!-----------------------------------------------------------------------

Subroutine vmalv(va, vb, v, q, lmax)

  Use :: nlte_type
  Use :: nlte_dim
  Implicit None

! ***  algebraic routine called from cmfray


! .. scalar arguments ..
  Integer (i4b) :: lmax
! ..
! .. array arguments ..
  Real (dp) :: q(lmax), v(lmax-1), va(lmax), vb(lmax)
! ..
! .. local scalars ..
  Real (dp) :: a, b
  Integer (i4b) :: l, lz
! ..

  lz = lmax - 1
  b = v(1)
  q(1) = vb(1)*b

  Do l = 2, lz
    a = b
    b = v(l)
    q(l) = va(l)*a + vb(l)*b
  End Do

  q(lmax) = va(lmax)*b

  Return
End Subroutine

!-----------------------------------------------------------------------

Subroutine vsub(a, b, n)

  Use :: nlte_type
  Use :: nlte_dim
  Implicit None

! a(i)=a(i)-b(i)


! .. scalar arguments ..
  Integer (i4b) :: n
! ..
! .. array arguments ..
  Real (dp) :: a(n), b(n)
! ..

  a(1:n) = a(1:n) - b(1:n)

  Return
End Subroutine

!-----------------------------------------------------------------------

Subroutine occup_at(nl, nu, xne, temp, xnhii, xnl, xnu)
  Use :: nlte_type
  Use :: nlte_dim
  Use :: fund_const, Only: saha, hk
  Use :: preformal_var, Only: gl, fl
  Implicit None

! ---- calculates approximate occup. for NL and NU, assuming dep=1.

! ..
! .. scalar arguments ..
  Real (dp) :: xne, temp, xnhii, xnl, xnu, const, ghii = 1.D0
  Integer (i4b) :: nl, nu

  const = 0.5D0*saha/ghii*xne*xnhii/temp**1.5
  xnl = const*gl(nl)*exp(hk*fl(nl)/temp)
  xnu = const*gl(nu)*exp(hk*fl(nu)/temp)

  Return
End Subroutine

!-----------------------------------------------------------------------

Subroutine median(x2, n, xmed)

! Find the median of X2(1), ... , X2(N), using as much of the quicksort
! algorithm as is needed to isolate it.

! Latest revision - 26 November 1996
! From A. Miller repository
! Modified to take real type from module share, C. Allende Prieto 2011
! Changed input array from x to x2, so that it is not changed, CAP 2017

  Use :: nlte_type

  Implicit None

  Integer (i4b), Intent (In) :: n
  Real (dp), Intent (In), Dimension (n) :: x2
  Real (dp), Intent (Out) :: xmed

! Local variables
  Real (dp), Dimension (n) :: x
  Real (dp) :: temp, xhi, xlo, xmax, xmin
  Logical :: odd
  Integer (i4b) :: hi, lo, nby2, nby2p1, mid, i, j, k

! make a copy to avoid altering it
  x = x2

  nby2 = n/2
  nby2p1 = nby2 + 1
  odd = .True.


! HI & LO are position limits encompassing the median.

  If (n==2*nby2) odd = .False.
  lo = 1
  hi = n
  If (n<3) Then
    If (n<1) Then
      xmed = 0.0
      Return
    End If
    xmed = x(1)
    If (n==1) Return
    xmed = 0.5*(xmed+x(2))
    Return
  End If

! Find median of 1st, middle & last values.

100 mid = (lo+hi)/2
  xmed = x(mid)
  xlo = x(lo)
  xhi = x(hi)
  If (xhi<xlo) Then !               Swap xhi & xlo
    temp = xhi
    xhi = xlo
    xlo = temp
  End If
  If (xmed>xhi) Then
    xmed = xhi
  Else If (xmed<xlo) Then
    xmed = xlo
  End If

! The basic quicksort algorithm to move all values <= the sort key (XMED)
! to the left-hand end, and all higher values to the other end.

  i = lo
  j = hi
110 Do
    If (x(i)>=xmed) Exit
    i = i + 1
  End Do
  Do
    If (x(j)<=xmed) Exit
    j = j - 1
  End Do
  If (i<j) Then
    temp = x(i)
    x(i) = x(j)
    x(j) = temp
    i = i + 1
    j = j - 1

!   Decide which half the median is in.

    If (i<=j) Go To 110
  End If

  If (.Not. odd) Then
    If (j==nby2 .And. i==nby2p1) Go To 130
    If (j<nby2) lo = i
    If (i>nby2p1) hi = j
    If (i/=j) Go To 120
    If (i==nby2) lo = nby2
    If (j==nby2p1) hi = nby2p1
  Else
    If (j<nby2p1) lo = i
    If (i>nby2p1) hi = j
    If (i/=j) Go To 120

!   Test whether median has been isolated.

    If (i==nby2p1) Return
  End If
120 If (lo<hi-1) Go To 100

  If (.Not. odd) Then
    xmed = 0.5*(x(nby2)+x(nby2p1))
    Return
  End If
  temp = x(lo)
  If (temp>x(hi)) Then
    x(lo) = x(hi)
    x(hi) = temp
  End If
  xmed = x(nby2p1)
  Return

! Special case, N even, J = N/2 & I = J + 1, so the median is
! between the two halves of the series.   Find max. of the first
! half & min. of the second half, then average.

130 xmax = x(1)
  Do k = lo, j
    xmax = max(xmax, x(k))
  End Do
  xmin = x(n)
  Do k = i, hi
    xmin = min(xmin, x(k))
  End Do
  xmed = 0.5*(xmin+xmax)

  Return
End Subroutine
