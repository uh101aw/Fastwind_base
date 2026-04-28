Module preformal_var

  Use :: nlte_type
  Use :: nlte_dim
  Implicit None

! for STAFIL and LINETAB

  Character (*), Parameter :: fpath = '../inicalc/DATA/'
! ----------------------------------------------------------------------
! comcha

! ID_LLEVS+2, to allow for one additional lower and upper level used
! within approximate treatment (NCOMP = 1!)
  Character (6), Dimension (id_atoms) :: labat
  Character (6), Dimension (id_llevs+2) :: labl
! ----------------------------------------------------------------------
! comnui

  Integer (i4b) :: nalin, nat, nl, nlevh21 ! level no of H21 for AT
  Integer (i4b), Dimension (id_llevs+2) :: le

  Logical :: at, srange1
! ----------------------------------------------------------------------
! comnur

  Real (dp), Dimension (id_atoms) :: weight, zeff
  Real (dp), Dimension (id_llevs+2) :: fl, gl, zl
! ----------------------------------------------------------------------
! charge,griem

  Integer (i4b) :: ns
  Real (dp) :: as, odop, ps, zz
  Real (dp), Dimension (69) :: ss, sx
! ----------------------------------------------------------------------
! inout

  Integer (i4b) :: iu, iou

End Module

!ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
!ccc                                                                cccc
!ccc      preparation of formal solucion                            cccc
!ccc                                                                cccc
!ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc

Subroutine preformal(file, ncomp, level, leveu, lineno, nstark, vturb, yhe, &
  balmer_lemke, pb_lemke, srange)

  Use :: nlte_type
  Use :: nlte_dim
  Use :: fund_const, Only: pi
  Use :: preformal_var, Only: fpath, nat, nl, fl, gl, zl, weight, zeff, iu, &
    iou, at, nlevh21, srange1
  Use :: ffr_error
  Implicit None

! ------ program should be run before the formal solution is performed
! ------ it uses one file as input:

! ATOM_FILE: contains the files of atomic data and stark profiles

! ------ and two files as output:

! i) IX.DAT: will be used by the formal solution as input
! ii) STARK.BIN: binary file containing the stark profile
! of the transition(s) considered

! stark broadening options:
! 0: pure doppler

! 1: stark (from tables or griem)
! for NSTARK = 0 or 1
! if lineno =0, standard approach
! if lineno =1, approximate treatment of transitions between
! high lying levels, and data from LINES.dat

! 2: voigt, with gammal and gammac from LINES.dat

! 3: voigt, with gammal and gammac approximated
! *works only for ncomp=1*:
! lineno   =0: gammal and gammac calculated
! lineno ne 0: gammal calculated, gammac=0.

! NOTE: if GammaL and GammaC not available, can be approximated quite well
! by using one line of the multiplet with lineno=0 and stark=3
! (checked, e.g., by comparing Si line values with given data)
! Works particularly good for lines with low lying lower level.


! ------ enrique santolaya, november 1993, version 1.0

! version 1.1 april 30th 1997 (j.puls)
!   to be used as subroutine in 'formalsol'


! version 1.2 july 4th 1997
!   ix.dat contains now gf-value,
!   broadening option "2" for voigt profiles incorporated


! version 3.0.0 jan 21st 1999
!   changed to f90, inclusion of modification by
!   christian wiethaus: wavelenghts now for air!

! version 3.1 march 11th 1999
!   more f90 changes, obsolete features removed

! version 4.0 april 16th 1999
!   ISO-lines can now be treated,
!   vturb included by additional folding of tabulated data

! version 4.1 may 20th 1999
!   GAMMAL and GAMMAC approximated for lines were no other
!   information given (corresponds to NSTARK = 3)

! version 4.1.1 june 16th 1999
!   error corrected (treatment of lines with two components
!   with the same parent levels). This error was corrected
!   by J. Puls, but new versions (4.2 and 5.0) have this
!   error. I have corrected it again for Carrie's version
!   (Miguel A. Urbaneja, Bremen, january 3rd 2002)

! version 4.2 june 28th 2000 inconsistency in griem (range = 6, &
!   xmax in err = 5) removed

! version 5.0 feb 7th 2001 file paths to STAFIL and LINETAB
!   changed via parameter fpath

! ----lineno  = 0 : line from detail input, assuming that levels are not packed
! ----lineno ne 0 : line from lines.dat, levels are packed

! version 5.1 april 2004 new broadening routine for "iso1" levels
!   according to the data by Dimitrijevic & Sahal-Brechot.
!   So far, proton/HeII contribution only approximate
!   from YHe and HeI input values (assumed to be constant
!   throughout the atmosphere)

! version 5.2 may 2006: stop statements in transic1 removed; allows
!   to have lines in LINES.DAT from levels which are NOT
!   present in the model

! version 6.0 july 2009: generalization to arbitrary NCOMP
!   old bug removed: calculation with pure Doppler-broadening
!   and more than one component possible

! version 6.1 april 2012: approximate treatment (AT) of transitions
!   between high lying levels (mm and submm range), for ALMA
!   coding by NSTARK = 0 or 1 and LINENO=1 (NCOMP=1)

! version 6.2 may 2013: Stark-broadening from Lemke (1997) will be used
!   for H_Balmer lines if BALMER_LEMKE = .true.

! version 7.0 july 2015: routines RDATOM, HALLCL, and IDENT renamed,
!   to RDATOM1, HALLCL1, and IDENT1,
!   because of inclusion and overlap with PRINCESA
!   subr. TRANSIC1: QCL = .FALSE. inserted
!   (old bug, found in March 2015)

! version 7.0.1 july 2016:
!   consistency with version 7.0 and 7.0.1 from standard_v10.4.1
!   inclusion of function IGENIO here (and no longer in formalsol_add.f90)

!   in function HALLCL, correction for
!   atmospheric refraction only down to 2000 AA:
!   Additional IF statement to avoid that UV lines are
!   corrected for air.
!   This had to be done because lines as HEII1640 and
!   Lyman_beta (not present in LINES.DAT) were corrected
!   while specific non-H/He lines such as NV use wavelengths
!   from LINES.DAT as given (here: in vaccum).

! version 7.0.2 march 2017: call ffrver(0) verification level 0,
!   to prevent unwanted output on fort.38

! version 7.0.3 oct 2020: compatible with gfortran

! version 7.0.3.1 nov 2022: potential Lemke broadening for
!   Paschen/Brackett lines (hydrogen) included

! version 7.0.3.2 feb 2023: few changes to allow for gfortran

! version 7.1 dec. 2025: polished

!--------------------------------------------------------------------------
  
!   NOTE on atmospheric refraction

!   1. standard path: SRANGE ('spectral range') = .FALSE.

!   lambda > 2000 A: refraction included
!   a) lines from LINES.DAT: lambda in air
!   b) ISO/ISO1 lines: lambda already in air
!   c) lines with wavelength from level energies:
!   corrected in HALL_CL

!   lambda < 2000 A: vacuum wavelength
!   a) lines from LINES.DAT: lambda in vacuum
!   b) ISO/ISO1 lines: thus far, only optical frequencies,
!   no re-correction required.

!   ATTENTION for future calculations!!!
!   c) lines with wavelength from level energies:
!   NOT corrected in HALL_CL

!   2. spectral ranges: SRANGE = .TRUE.
!   ALL wavelengths in vacuum
!   a) lines from LINES.DAT: not used so far,
!   ATTENTION for future calculations!!!
!   b) ISO/ISO1 lines: lambda in air, re-corrected in PRESTARK
!   c) lines with wavelength from level energies:
!   NOT corrected in HALL_CL

!   initializing data
! .. scalar arguments ..
  Integer (i4b) :: ncomp
  Real (dp) :: vturb, yhe
  Logical :: balmer_lemke, pb_lemke, srange
! ..
! .. array arguments ..
  Integer (i4b), Dimension (ncomp) :: lineno, nstark
  Character, Dimension (ncomp) :: level*6, leveu*6
! ..
! .. local scalars ..
  Integer (i4b) :: i, ns, nsum
  Logical :: rtable
  Character :: dc*6, fichero*32, linetab*32, stafil*32, ret*4, file*60
! ..
! .. local arrays ..
  Real (dp), Allocatable, Dimension (:) :: flu, gammac, gammal, xlamb, sumlo, &
    sumup
  Integer (i4b), Allocatable, Dimension (:) :: natom, nlevl, nlevu, index
  Logical, Allocatable, Dimension (:) :: low
! ..
! .. external subroutines ..
  External :: ffracc, prestark, rdatom1, transic1, transic2, transic_h_at
! ..
! .. intrinsic functions ..
! ..

  Allocate (flu(ncomp), gammac(ncomp), gammal(ncomp), xlamb(ncomp), &
    sumlo(ncomp), sumup(ncomp))
  Allocate (natom(ncomp), nlevl(ncomp), nlevu(ncomp), index(ncomp))
  Allocate (low(ncomp))

  srange1 = srange
  rtable = .False.
  at = .False.
  low = .False.

  nlevh21 = 0

  nl = 0
  nat = 0
  weight = 0.
  zeff = 0.
  fl = 0.
  gl = 0.
  zl = 0.

  natom = 0
  nlevl = 0
  nlevu = 0
  xlamb = 0.
  flu = 0.
  gammal = 0.
  gammac = 0.
  sumlo = 0.
  sumup = 0.

  Do i = 1, ncomp
    If (nstark(i)==2 .Or. lineno(i)/=0) rtable = .True. ! works also for
!   nstark=4
    If (nstark(i)==3 .And. lineno(i)/=0) Stop ' NSTARK = 3 AND LINENO NE 0'
    If (nstark(i)<=1 .And. lineno(i)==1) at = .True.
  End Do

  iu = 37
  iou = 38
  dc = ':T'

  Open (1, File='ATOM_FILE', Status='OLD')
  Rewind 1
  Read (1, Fmt='(A)') fichero
  Read (1, Fmt='(A)') stafil
  If (rtable) Read (1, Fmt='(A)') linetab
  Close (1)

  Open (Unit=iu, File=fichero, Form='FORMATTED', Status='OLD')
  Open (Unit=iou, File='control.dat', Status='UNKNOWN')

! rewind!!!

  iu = -37

  Call ffracc(iu, ret)
  If (on_error(ret)) Go To 110

rloop: Do

    Call rdatom1(ret)
    If (ret/='RET0' .And. ret/='RET3') Go To 120
    If (ret=='RET3') Go To 100

    If (.Not. at) Then
!     standard path
      Call transic1(level, leveu, ncomp, lineno, nlevl, nlevu, natom, xlamb, &
        flu, low, sumlo, sumup, nstark, gammac, ret)
      If (on_error(ret)) Go To 120

    End If
  End Do rloop

! normal exit: end of file

100 Continue

  Write (*, Fmt='(A)') ' DETAIL FILE IS OVER. SUCCESSFUL!! '
  Close (iu)

! JO March 2017 (to prevent unwanted output to file fort.38)'
  Call ffrver(0) !                  set verification level to zero

! now, all contributions to SUMLO, SUMUP are calculated, &
! hence we can calculate GAMMAL.
! the constant below corresponds to 8 pi^2 e^2/(m_e c), if lambda in A

  If (.Not. at) Then
    Do i = 1, ncomp
      If (nstark(i)==3) Then
        gammal(i) = 6.6702082D15*(sumlo(i)+sumup(i))/(4.*pi)
        Print *, gammal(i), gammac(i)
      End If
    End Do
  Else
!   approximate treatment
    Call transic_h_at(level, leveu, ncomp, lineno, nlevl, nlevu, natom, xlamb, &
      flu, nstark)
  End If

  If (rtable) Then
    Open (Unit=iu, File=fpath//linetab, Form='FORMATTED', Status='OLD')

!   rewind!!!

    iu = -37
    Call ffracc(iu, ret)
    If (on_error(ret)) Go To 110

!   note here: flu corresponds to gf/sum(g), see above

    Call transic2(level, leveu, ncomp, nstark, lineno, nlevl, nlevu, natom, &
      xlamb, flu, gammal, gammac, ret)
    If (on_error(ret)) Go To 120

    Write (*, Fmt='(A)') ' LINE FILE IS OVER. SUCCESSFUL!! '
    Close (iu)
  End If

  Go To 130

! error in input unit exit

110 Continue
  Write (*, Fmt='(A)') ' ERROR IN INPUT UNIT NUMBER '
  Stop

! error conditions
120 Select Case (ret)
  Case ('RET1') !                   end of file "no esperado" exit
    Write (*, Fmt='(A)') ' EOF NOT EXPECTED. THERE IS SOMETHING WRONG! '
    Stop
  Case ('RET2') !                   error in ffr subroutines exit
    Write (*, Fmt='(A)') ' ERROR IN FFR SUBROUTINES '
    Stop
  Case Default
    Stop ' WRONG ERROR CONDITION IN PREFORMAL'
  End Select

130 Continue
  Close (iou)

  Do i = 1, ncomp
    If (nlevl(i)==0 .Or. nlevu(i)==0) Stop 'LEVELS NOT FOUND'
    If (xlamb(i)<=0.D0) Stop 'ERROR IN LAMBDA'
    If (flu(i)<=0.D0) Stop 'ERROR IN FLU'
    If (natom(i)==0) Stop 'ERROR IN ATOM'
  End Do

! ---    sorting for wavelengths (blue component at first)

  Call indexx(ncomp, xlamb, index)

! ---    interchange of components

  xlamb = xlamb(index)
  flu = flu(index)
  gammal = gammal(index)
  gammac = gammac(index)
  nlevl = nlevl(index)
  nlevu = nlevu(index)
  natom = natom(index)
  nstark = nstark(index)

! ---    from here on, comp 1 = blue component


! ------ nstark=0 : doppler
! ------ nstark=1 : stark
! ------ nstark=2 : voigt with data from LINETAB (LINES.DAT)
! ------ nstark=3 : voigt with approx. for GAMMAL/GAMMAC (see above comments)


  nsum = 0
  Do i = 1, ncomp
    If (nstark(i)==1) nsum = nsum + 1
  End Do

! ------ nsum=0 : no stark;
! ------ nsum=1 : only one comp. stark (nsum, precs.)
! ------ nsum=n : n comps. stark
  Print *, 'NUMBER OF STARK-COMPONENTS', nsum

  If (nsum/=0) Call prestark(file, nsum, stafil, ncomp, nlevl, nlevu, xlamb, &
    flu, nstark, vturb, yhe, balmer_lemke, pb_lemke)

! ------ output file ix.dat

  Open (1, File=trim(file)//'/IX.DAT', Status='UNKNOWN', Form='FORMATTED')
  Rewind 1

  Write (1, Fmt=*) ncomp
  Do i = 1, ncomp
    Write (1, Fmt=*) xlamb(i)
    Write (1, Fmt=*) gl(nlevl(i)), gl(nlevu(i))
    Write (1, Fmt=*) gl(nlevl(i))*flu(i)
    Write (1, Fmt=*) nlevl(i), nlevu(i)
    Write (1, Fmt=*) weight(natom(i))
    ns = nstark(i)
    Write (1, Fmt=*) ns
    If (ns==2 .Or. ns==3) Write (1, Fmt=*) gammal(i), gammac(i)
  End Do
  Close (1)

  Deallocate (flu, gammac, gammal, xlamb, sumlo, sumup)
  Deallocate (natom, nlevl, nlevu, index, low)

  Return

End Subroutine

!ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc

Subroutine rdatom1(retchar)

  Use :: nlte_type
  Use :: nlte_dim
  Use :: preformal_var, Only: labat, labl, nat, nl, le, fl, gl, zl, weight, &
    zeff
  Use :: ffr_error
  Implicit None

! .. local scalars ..
  Real (dp) :: realo, zef
  Integer (i4b) :: integ, nc
  Character :: label*4, key*6, retchar*4, ret*4
! ..
! .. external subroutines ..
  External :: ffrkey, ffrloc, ffrnum
! ..

  label = 'ATOM'
  Call ffrloc(label, ret)
  Select Case (ret)
  Case ('RET1')
    Go To 140
  Case ('RET2')
    Go To 150
  Case ('RET0')
    Continue
  Case Default
    Stop ' WRONG RETURN IN RDATOM1'
  End Select

  nat = nat + 1

  Call ffrkey(key, nc, ret)
  If (on_error(ret)) Go To 150
  labat(nat) = key

  Call ffrnum(realo, integ, ret)
  If (on_error(ret)) Go To 150
  zeff(nat) = realo

  Call ffrnum(realo, integ, ret)
  If (on_error(ret)) Go To 150
  weight(nat) = realo

  Call ffrnum(realo, integ, ret)
  If (on_error(ret)) Go To 150

  zef = zeff(nat) - 1.D0

100 Continue

  Call ffrkey(key, nc, ret)
  If (on_error(ret)) Go To 150

! ------------- l  levels
  If (key=='L') Then
    zef = zef + 1.D0
110 Continue
    Call ffrkey(key, nc, ret)
    If (on_error(ret)) Go To 150
    If (key=='0') Go To 100
    nl = nl + 1
    labl(nl) = key
    zl(nl) = zef

    Call ffrnum(realo, integ, ret)
    If (on_error(ret)) Go To 150
    gl(nl) = realo

    Call ffrnum(realo, integ, ret)
    If (on_error(ret)) Go To 150
    fl(nl) = realo
    le(nl) = nat
    Call ffrkey(key, nc, ret)
    If (on_error(ret)) Go To 150

    Go To 110

!   ----------------- x  levels
  Else If (key=='X') Then
120 Continue
    Call ffrkey(key, nc, ret)
    If (on_error(ret)) Go To 150

    If (key=='0') Go To 100

    Call ffrnum(realo, integ, ret)
    If (on_error(ret)) Go To 150

    Call ffrnum(realo, integ, ret)
    If (on_error(ret)) Go To 150

    Call ffrkey(key, nc, ret)
    If (on_error(ret)) Go To 150

    Go To 120

!   -------------- s  levels
  Else If (key=='S') Then
130 Continue
    Call ffrkey(key, nc, ret)
    If (on_error(ret)) Go To 150

    If (key=='0') Go To 100

    Call ffrnum(realo, integ, ret)
    If (on_error(ret)) Go To 150

    Call ffrnum(realo, integ, ret)
    If (on_error(ret)) Go To 150

    Call ffrnum(realo, integ, ret)
    If (on_error(ret)) Go To 150

    Call ffrnum(realo, integ, ret)
    If (on_error(ret)) Go To 150

    Call ffrnum(realo, integ, ret)
    If (on_error(ret)) Go To 150

    Call ffrnum(realo, integ, ret)
    If (on_error(ret)) Go To 150

    Call ffrkey(key, nc, ret)
    If (on_error(ret)) Go To 150

    Go To 130

!   --------------------  k  levels
  Else If (key=='K') Then
    nl = nl + 1
    zef = zef + 1
    Call ffrkey(key, nc, ret)
    If (on_error(ret)) Go To 150
    labl(nl) = key

    Call ffrnum(realo, integ, ret)
    If (on_error(ret)) Go To 150
    gl(nl) = realo

    Call ffrnum(realo, integ, ret)
    fl(nl) = realo
    zl(nl) = zef

    Call ffrkey(key, nc, ret)
    If (on_error(ret)) Go To 150

  End If

! the atom is over.
  retchar = 'RET0'
  Return
! no hay mas atomos que leer. ya se ha leido todo el fichero. es
! cuando acaba el programa principal (prince)
140 Continue

  retchar = 'RET3'
  Return

! error handling
150 Select Case (ret)
  Case ('RET1') !                   end of file exit
    Write (*, Fmt='(A)') ' END OF FILE '
    retchar = 'RET1'
    Return
  Case ('RET2') !                   error exit
    Write (*, Fmt='(A)') 'ERROR IN FFR SUBROUTINES '
    retchar = 'RET2'
    Return
  Case Default
    Stop ' WRONG ERROR CONDITION IN RDATOM1'
  End Select

End Subroutine

!ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc

Subroutine transic1(level, leveu, ncomp, lineno, nlevl, nlevu, natom, xlamb, &
  flu, low, sumlo, sumup, nstark, gammac, retchar)

  Use :: nlte_type
  Use :: nlte_dim
  Use :: preformal_var, Only: gl, labl, nl, le, iu
  Use :: ffr_error
  Implicit None

! subrutina que lee los datos referentes a las transiciones, los
! almacena y escribe (para control) en un fichero de salida.
! espero poner suficientes explicaciones a lo largo de la
! subrutina, de modo que este todo lo mas claro posible.

! .. scalar arguments ..
  Integer (i4b) :: ncomp
! ..
! .. array arguments ..
  Real (dp) :: flu(ncomp), xlamb(ncomp), sumlo(ncomp), sumup(ncomp), &
    gammac(ncomp)
  Integer (i4b) :: lineno(ncomp), natom(ncomp), nlevl(ncomp), nlevu(ncomp), &
    nstark(ncomp)
  Logical :: low(ncomp)
  Character :: level(ncomp)*6, leveu(ncomp)*6
! ..
! .. local scalars ..
  Real (dp) :: realo, tl, wlc
  Integer (i4b) :: i, iform, integ, iqi, labl1, labu1, m, nc, ndatm, ndatos, &
    nlong, nn, nnu
  Logical :: qcl, qrt
  Character :: key*6, retchar*4, ret*4
! ..
! .. local arrays ..
  Character :: tabla(8)*6
! ..
! .. external functions ..
  Real (dp) :: hallcl1, calcgam
  External :: hallcl1, calcgam
! ..
! .. external subroutines ..
  External :: ffrkey, ffrnum, ident1
! ..
! .. data statements ..
  Data tabla/'RBB', 'CBB', 'RBX', 'CBX', 'CBS', 'CBF', 'RBF', 'RFF'/
! ..

100 Continue
! JO March 2015: this is a very important statement, and had been forgotten
! in previous versions; if left out, certain reading errors can occur.
  qcl = .False.

  Call ffrkey(key, nc, ret)

  If (on_error(ret)) Go To 230

  If (key=='TY') Then
    Call ffrkey(key, nc, ret)
    If (on_error(ret)) Go To 230

    Call ident1(key, tabla, 8, m)
    Call ffrnum(realo, iform, ret)
    If (on_error(ret)) Go To 230

    Call ffrnum(realo, ndatos, ret)
    If (on_error(ret)) Go To 230
!   son los numeros de formula y de datos, respectivamente
    ndatm = ndatos

    Go To 100

!   muestreo en anchuras doppler. el muestreo queda en dl.
  Else If (key=='DB') Then
    Call ffrnum(realo, nlong, ret)
    If (on_error(ret)) Go To 230
    Do i = 1, nlong
      Call ffrnum(realo, integ, ret)
      If (on_error(ret)) Go To 230
    End Do

    Go To 100

!   muestreo en angstroms
  Else If (key=='DL') Then
    Call ffrnum(realo, nlong, ret)
    If (on_error(ret)) Go To 230
    Do i = 1, nlong
      Call ffrnum(realo, integ, ret)
      If (on_error(ret)) Go To 230
    End Do

    Go To 100

!   temperatura de linea
  Else If (key=='TL') Then
    Call ffrnum(tl, integ, ret)
    If (on_error(ret)) Go To 230

    Go To 100

!   longitud central de la linea
  Else If (key=='CL') Then
    qcl = .True.
    Go To 100

!   tipo de transicion
  Else
    Call ident1(key, tabla, 8, m)
    If (m==0) Then
      Go To 240
    Else
      qrt = .False.
      If ((m==1) .Or. (m==3)) qrt = .True.
      Go To (110, 110, 150, 150, 170, 170, 190, 210) m
!     goto calculado, para que haga una cosa u otra segun el tipo de trans.


    End If

  End If

! -------------transiciones rbb o cbb ------------
110 Continue

  Call ffrkey(key, nc, ret)
  If (on_error(ret)) Go To 230

120 Continue
  Do i = 1, ncomp
    low(i) = .False.
    If (key==level(i) .And. lineno(i)==0) low(i) = .True.
  End Do

  Do iqi = 1, nl
    If (labl(iqi)==key) Then
      nn = iqi
      Go To 130
    End If
  End Do

  Stop 'ERROR IN TRANSIC1 - LABEL L NOT FOUND, BB'
130 Continue

  labl1 = nn

  Call ffrkey(key, nc, ret)
! print*,labl(labl1),key

  If (on_error(ret)) Go To 230

  Do iqi = 1, nl
    If (labl(iqi)==key) Then
      nnu = iqi
      Go To 140
    End If
  End Do
  Stop 'ERROR IN TRANSIC1 - LABEL U NOT FOUND, BB'
140 Continue

  labu1 = nnu
! print*,labl(labl1),' ',labl(labu1),qrt,qcl

! transicion radiativa
  If (qrt) Then
    If (qcl) Then
      Call ffrnum(wlc, integ, ret)
      If (on_error(ret)) Go To 230
    Else
      wlc = hallcl1(labl1, labu1)
    End If

    Do i = 1, ncomp

!     sumlo and sumup contain sum of (gi fij/gj/lambda^2) from all transitions
!     downwards from lower and upper level. Needed to calculate Sum(Aji) for
!     GAMMAL!
!     JO changed Sept. 2015; new version might be still erroneous, partic. for
!     NSTARK=2
!     -according to header, this works only for NSTARK=3 and NCOMP=1
!     -as far as I can see, this works only if in the component list
!      any upper and lower level does not appear more than once
!      (which is usually the case for individual lines complexes).
!     -does no longer work for 'range' calculations
!     -in any case, GAMMAL only calculated for NSTARK=3
!     -thus, seperate reading of oscillator strength
!      (question: what about oscillator strength if NSTARK = 2 in old version???)

!     IF(NSTARK(I).NE.2 .AND. LINENO(I).EQ.0) THEN !old version
      If (nstark(i)==3 .And. lineno(i)==0) Then ! new version (Sept. 2015)
        If (labl(nnu)==level(i)) Then
          Call ffrnum(realo, integ, ret)
          If (on_error(ret)) Go To 230
          ndatm = ndatos - 1
          sumlo(i) = sumlo(i) + gl(nn)*realo/(gl(nnu)*wlc*wlc)
!         PRINT*,'LO ',I,' ',LABL(NN),LABL(NNU),GL(NN),GL(NNU),WLC,REALO
        Else If (labl(nnu)==leveu(i)) Then
          Call ffrnum(realo, integ, ret)
          If (on_error(ret)) Go To 230
          ndatm = ndatos - 1
          sumup(i) = sumup(i) + gl(nn)*realo/(gl(nnu)*wlc*wlc)
!         PRINT*,'UP ',I,' ',LABL(NN),LABL(NNU),GL(NN),GL(NNU),WLC,REALO
        End If
!       new else block, allowing to read osc. strength
      Else
        If (key==leveu(i) .And. low(i)) Then
          Call ffrnum(realo, integ, ret)
          If (on_error(ret)) Go To 230
          ndatm = ndatos - 1
        End If
      End If


      If (key==leveu(i) .And. low(i)) Then
        xlamb(i) = wlc
!       osc. strength already read
        flu(i) = realo
        nlevl(i) = labl1
        nlevu(i) = labu1
        natom(i) = le(labl1)
        If (nstark(i)==3) gammac(i) = calcgam(labl1, labu1, natom(i))
      End If
    End Do

  End If

! lectura de los datos
  Do i = 1, ndatm
    Call ffrnum(realo, integ, ret)
    If (on_error(ret)) Go To 230
  End Do

  ndatm = ndatos

! leer nueva keyword
  Call ffrkey(key, nc, ret)
  If (on_error(ret)) Go To 230

  If (key=='0') Then
    qcl = .False.
    Go To 100
  End If

  Go To 120

! --------------- transiciones rbx o cbx ------------
150 Continue
  Call ffrkey(key, nc, ret)
  If (on_error(ret)) Go To 230

  Call ffrkey(key, nc, ret)
  If (on_error(ret)) Go To 230
160 Continue

  Call ffrkey(key, nc, ret)
  If (on_error(ret)) Go To 230

! transicion radiativa
  If (qrt) Then
    If (qcl) Then
      Call ffrnum(wlc, integ, ret)
      If (on_error(ret)) Go To 230
    End If
  End If

! lectura de datos
  Do i = 1, ndatos
    Call ffrnum(realo, integ, ret)
    If (on_error(ret)) Go To 230
  End Do

! leer nueva keyword
  Call ffrkey(key, nc, ret)
  If (on_error(ret)) Go To 230

  If (key=='0') Then
    qcl = .False.
    Go To 100
  End If

  Go To 160

! ------------- transiciones  cbs o cbf ------------------
170 Continue
  Call ffrkey(key, nc, ret)
  If (on_error(ret)) Go To 230

  Call ffrkey(key, nc, ret)
  If (on_error(ret)) Go To 230

180 Continue

! lectura de datos
  Do i = 1, ndatos
    Call ffrnum(realo, integ, ret)
    If (on_error(ret)) Go To 230
  End Do

! leer nueva keyword
  Call ffrkey(key, nc, ret)
  If (on_error(ret)) Go To 230

  If (key=='0') Go To 100

  Go To 180

! ------------------ transicion  rbf -------------------
190 Continue
  Call ffrkey(key, nc, ret)
  If (on_error(ret)) Go To 230

  Call ffrkey(key, nc, ret)
  If (on_error(ret)) Go To 230

200 Continue
  Call ffrkey(key, nc, ret)
  If (on_error(ret)) Go To 230

  Call ffrkey(key, nc, ret)
  If (on_error(ret)) Go To 230

! leer datos
  Do i = 1, ndatos
    Call ffrnum(realo, integ, ret)
    If (on_error(ret)) Go To 230
  End Do

! leer nueva keyword
  Call ffrkey(key, nc, ret)
  If (on_error(ret)) Go To 230

  If (key=='0') Go To 100

  Go To 200

! ----------------- transiciones rff ------------------------
210 Continue

  Call ffrkey(key, nc, ret)
  If (on_error(ret)) Go To 230

220 Continue

! lectura de datos
  Do i = 1, ndatos
    Call ffrnum(realo, integ, ret)
    If (on_error(ret)) Go To 230
  End Do

! leer nueva keyword
  Call ffrkey(key, nc, ret)
  If (on_error(ret)) Go To 230

  If (key=='0') Go To 100

  Go To 220

! error handling
230 Select Case (ret)
  Case ('RET1') !                   end of file exit
    Write (*, Fmt='(A)') ' END OF FILE '
    retchar = 'RET1'
    Return
  Case ('RET2') !                   error exit
    Write (*, Fmt='(A)') 'ERROR IN FFR SUBROUTINES '
    retchar = 'RET2'
    Return
  Case Default
    Stop ' WRONG ERROR CONDITION IN TRANSIC1'
  End Select

240 Continue
  Backspace iu

  Return

End Subroutine

!ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc

Subroutine transic2(level, leveu, ncomp, nstark, lineno, nlevl, nlevu, natom, &
  xlamb, flu, gammal, gammac, retchar)

  Use :: nlte_type
  Use :: nlte_dim
  Use :: preformal_var, Only: labl, nl, le, at
  Use :: ffr_error
  Implicit None

! finds line-transitions from file lines.dat

! .. scalar arguments ..
  Integer (i4b) :: ncomp
! ..
! .. array arguments ..
  Real (dp) :: flu(ncomp), gammac(ncomp), gammal(ncomp), xlamb(ncomp)
  Integer (i4b) :: lineno(ncomp), natom(ncomp), nlevl(ncomp), nlevu(ncomp), &
    nstark(ncomp)
  Character :: level(ncomp)*6, leveu(ncomp)*6, retchar*4, ret*4
! ..
! .. local scalars ..
  Real (dp) :: realo, sumglo
  Integer (i4b) :: i, integ, iqi, labl1, labu1, nc, ndatos, nline, nn, nnu, &
    iqistart
  Character :: key*6
! ..
! .. local arrays ..
  Logical :: low(ncomp)
! ..
! .. external subroutines ..
  External :: ffrkey, ffrnum
! ..

100 Continue
  Call ffrkey(key, nc, ret)
  If (on_error(ret)) Go To 150
  If (key/='CL') Go To 160

  Call ffrkey(key, nc, ret)
  If (on_error(ret)) Go To 150
  If (key/='TY') Go To 160

  Call ffrkey(key, nc, ret)
  If (on_error(ret)) Go To 150
  If (key/='RBB') Go To 160

  Call ffrnum(realo, integ, ret)
  If (on_error(ret)) Go To 150
  If (integ/=2) Go To 160

  Call ffrnum(realo, integ, ret)
  If (on_error(ret)) Go To 150
  ndatos = integ
  Call ffrkey(key, nc, ret)
  If (on_error(ret)) Go To 150
  If (key/='RBB') Go To 160

  Call ffrkey(key, nc, ret)
  If (on_error(ret)) Go To 150

110 Continue
  Do i = 1, ncomp
    low(i) = .False.
    If (key==level(i) .And. lineno(i)/=0) low(i) = .True.
  End Do

  iqistart = 1
  If (at) Then
    If (nl/=id_llevs+2) Stop ' AT AND NL NE ID_LLEVS+2 IN TRANSIC2'
    iqistart = id_llevs + 1
  End If

  Do iqi = iqistart, nl
    If (labl(iqi)==key) Then
      nn = iqi
      Go To 120
    End If
  End Do
! stop statement removed: not ALL lines in LINES.dat must be present in
! level-list
! STOP 'ERROR IN TRANSIC2 - LABEL L NOT FOUND, BB'
120 Continue


  labl1 = nn

  Call ffrkey(key, nc, ret)
  If (on_error(ret)) Go To 150

  Do iqi = iqistart, nl
    If (labl(iqi)==key) Then
      nnu = iqi
      Go To 130
    End If
  End Do
! stop statement removed: not ALL lines in LINES.dat must be present in
! level-list
! STOP 'ERROR IN TRANSIC2 - LABEL U NOT FOUND, BB'

130 Continue
  labu1 = nnu

  Call ffrnum(realo, integ, ret)
  If (on_error(ret)) Go To 150

  nline = integ

  Do i = 1, ncomp
    If (key==leveu(i) .And. low(i) .And. nline==lineno(i)) Then
      Call ffrnum(realo, integ, ret)
      If (on_error(ret)) Go To 150
      xlamb(i) = realo

      Call ffrnum(realo, integ, ret)
      If (on_error(ret)) Go To 150
      flu(i) = realo

      Call ffrnum(realo, integ, ret)
      If (on_error(ret)) Go To 150
      gammal(i) = realo

      Call ffrnum(realo, integ, ret)
      If (on_error(ret)) Go To 150

      Call ffrnum(realo, integ, ret)
      If (on_error(ret)) Go To 150

      Call ffrnum(realo, integ, ret)
      If (on_error(ret)) Go To 150
      gammac(i) = realo

      Call ffrnum(realo, integ, ret)
      If (on_error(ret)) Go To 150
      sumglo = realo
      nlevl(i) = labl1
      nlevu(i) = labu1
      natom(i) = le(labl1)

      Go To 140

    End If

  End Do

! lectura de los datos
  Do i = 1, ndatos
    Call ffrnum(realo, integ, ret)
    If (on_error(ret)) Go To 150
  End Do

! leer nueva keyword
140 Continue
  Call ffrkey(key, nc, ret)
  If (on_error(ret)) Go To 150

  If (key=='0') Go To 100

  Go To 110

! error handling

150 Select Case (ret)
  Case ('RET1') !                   end of file exit
    Write (*, Fmt='(A)') ' END OF FILE '
    retchar = 'RET1'
    Return
  Case ('RET2') !                   error exit
    Write (*, Fmt='(A)') 'ERROR IN FFR SUBROUTINES '
    retchar = 'RET2'
    Return
  Case Default
    Stop ' WRONG ERROR CONDITION IN TRANSIC2'
  End Select

160 If (key=='THEEND') Then
    retchar = 'RET0'
    Return
  Else
    Stop ' WRONG END CONDITION IN TRANSIC2'
  End If

End Subroutine
!-----------------------------------------------------------------------

Function calcgam(ii1, ii2, natom)

  Use :: nlte_type
  Use :: nlte_dim
  Use :: preformal_var, Only: fl, zl, weight
  Use :: fund_const, Only: clight, pi
  Implicit None

! calculated GAMMAC in approximate weigh, using the relation between
! neff and energy

! ..
  Real (dp) :: calcgam

! .. scalar arguments ..
  Integer (i4b) :: ii1, ii2, natom
! ..
  Real (dp) :: ryd, z

  z = zl(ii1)
  If (z/=zl(ii2)) Stop ' SOMETHING WRONG WITH NET CHARGE'
  z = z + 1.

  ryd = 109737.312*(1.+1./(1836.*weight(natom)))
  calcgam = 4.335D-7/(4.*pi)*z*z*(ryd*clight)**2*(1./(fl(ii2)*fl(ii2))+1./(fl( &
    ii1)*fl(ii1)))

  Return

End Function


!-----------------------------------------------------------------------

Function hallcl1(ii1, ii2)

  Use :: nlte_type
  Use :: nlte_dim
  Use :: preformal_var, Only: fl, srange1
  Implicit None

! calcula la longitud de onda central de una transicion rbb,
! a partir de los datos atomicos de los niveles involucrados

! ..
  Real (dp) :: hallcl1

! .. scalar arguments ..
  Integer (i4b) :: ii1, ii2
! ..
! .. local scalars ..
  Real (dp) :: clight, fc, ff1, ff2, xkw, xn
! ..
! .. intrinsic functions ..

! ..
! .. data statements ..
  Data clight/2.997925E18/
! ..
  ff1 = fl(ii1)
  ff2 = fl(ii2)
  fc = abs(ff1-ff2)
  If (fc==0.D0) fc = 1.D0
  hallcl1 = clight/fc
! ---    primitive way to change hei 4026 line (to avoid
! ---    false wavelength and order of components)
  If (abs(hallcl1-4026.72472212120D0)<1.D-10) hallcl1 = 4027.3D0
! --- same for hei 4388
  If (abs(hallcl1-4388.82518552703D0)<1.D-10) hallcl1 = 4389.16D0

! -----AVOID THAT UV LINES HAVE THEIR WAVELENGTHS CONVERGED TO AIR
! -----CHANGE MADE BY Luiz Carneiro - 25/02/2015
  If (srange1 .Or. hallcl1<=2000.D0) Return

! ------CONVERSION TO WAVELENGTH IN AIR

  xkw = 1.D4/hallcl1
  xn = 1.D0 + 1.D-7*(643.28D0+294981.D0/(146.D0-xkw**2)+2554.D0/(41.D0-xkw**2) &
    )
  hallcl1 = hallcl1/xn
  Return

End Function

!-----------------------------------------------------------------------

Subroutine ident1(key, table, nmax, j)

  Use :: nlte_type
  Use :: nlte_dim
  Implicit None

! compara el caracter key con los de la tabla 'table', y devuelve
! en j el indice que corresponda a key en table, si esta, o 0 en
! otro caso.
! .. scalar arguments ..
  Integer (i4b) :: j, nmax
  Character :: key*6
! ..
! .. array arguments ..
  Character :: table(*)*6
! ..
! .. local scalars ..
  Integer (i4b) :: i
! ..
  i = 1
100 Continue

  If (key==table(i)) Then
    j = i
    Return
  End If

  i = i + 1
  If (i<=nmax) Then
    Go To 100
  End If

  j = 0

  Return

End Subroutine

!ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc

Subroutine prestark(file, nsum, stafil, nb, nlevl, nlevu, xlamb, flu, nstark, &
  vturb, yhe, balmer_lemke, pb_lemke)

  Use :: nlte_type
  Use :: nlte_dim
  Use :: preformal_var, Only: fpath, labl, le, weight, iu, srange1
  Use :: ffr_error
  Implicit None

! ---- preparation of stark solution. stark profiles are read and
! ---- the desired transition(s) is(are) stored in file 'stark.bin'

! .. parameters ..
  Integer (i4b), Parameter :: maxw = id_maxww, maxt = id_maxtt, &
    maxne = id_maxne, maxps = id_maxps, ndatamax = 50
! ..
! .. scalar arguments ..
  Integer (i4b) :: nb, nsum
  Real (dp) :: vturb, yhe
  Logical :: balmer_lemke, pb_lemke
  Character :: stafil*32, file*60
! ..
! .. array arguments ..
  Real (dp) :: xlamb(nb), flu(nb)
  Integer (i4b) :: nlevl(nb), nlevu(nb), nstark(nb), ll, lu
! ..
! .. local scalars ..
  Real (dp) :: aa, bb, ee, fac, osc, peso, realo, tt, xl0, waveair, xkw, xn
  Integer (i4b) :: i, ijk, integ, j, je, jt, mfull, nc, ndata, nin, nn, nn1, &
    nn2, nnat, nnn, nps, ntsiso, nfsiso
  Logical :: qlog
  Character :: key*6, key2*6, ret*4, llchar*6, luchar*6
! ..
! .. local arrays ..
  Real (dp) :: data(ndatamax), ddws(35), sq1(35), eees(11), sq(maxw), ttts(6), &
    prof(maxw)
  Integer (i4b) :: nwsn
  Logical :: q

  Real (dp), Allocatable, Dimension (:, :) :: dws, es, p, ts
  Integer (i4b), Allocatable, Dimension (:) :: nes, nts, nws
  Logical, Allocatable, Dimension (:) :: done, qhalf

! ..
! .. external functions ..
  Integer (i4b) :: nppal
  External :: nppal
! ..
! .. external subroutines ..
  External :: ffracc, ffrkey, ffrloc1, ffrnum, griem
! ..
! .. intrinsic functions ..
! ..
! .. data statements ..

  Data ttts, eees/4., 4.301, 4.602, 4.903, 5.204, 5.505, 10.5, 11.5, 12., &
    12.5, 13., 13.5, 14., 14.5, 15., 16., 17./
! ..

! IF(NSUM.NE.1.AND.NSUM.NE.2) STOP ' MORE THAN TWO STARK-COMPONENTS, CHECK!'

  Allocate (dws(maxw,nb), es(maxne,nb), p(maxps,nb), ts(maxt,nb))
  Allocate (nes(nb), nts(nb), nws(nb))
  Allocate (done(nb), qhalf(nb))

  done = .True.

  Open (1, File=fpath//stafil, Status='OLD', Form='FORMATTED')

! NEEDS TO BE CHANGED FOR MORE THAN TWO STARK LINES
  Do i = 1, nb
    If (nstark(i)==1) done(i) = .False.
  End Do

! BIG_LOOP

big_loop: Do nnn = 1, nb

    If (nstark(nnn)/=1) Cycle


!   rewind!!!

    iu = -1
    Call ffracc(iu, ret)
    Select Case (ret)
    Case ('RET1')
      Go To 120
    Case ('RET0')
      Continue
    Case Default
      Stop ' WRONG RETURN IN PRESTARK'
    End Select

!   ---- read stafil

100 Continue

!   note! 'line' has to be situated within first five characters of
!   input-line

    Call ffrloc1('LINE', ret)
    Select Case (ret)
    Case ('RET1')
      Go To 110
    Case ('RET2')
      Go To 130
    Case ('RET0')
      Continue
    Case Default
      Stop ' WRONG RETURN IN PRESTARK'
    End Select

    Call ffrkey(key, nc, ret)
    Select Case (ret)
    Case ('RET1')
      Go To 130
    Case ('RET2')
      Go To 130
    Case ('RET0')
      Continue
    Case Default
      Stop ' WRONG RETURN IN PRESTARK'
    End Select

!   ----------------------------------------------------------------------

!   frequency grid for ISO lines

    If (key=='DL') Then

      Call ffrnum(realo, nc, ret)
      If (on_error(ret)) Go To 130
      If (nc>maxw) Stop ' NWS(ISO) > MAXW'
      nws(nnn) = nc
      Do i = 1, nc
        Call ffrnum(realo, integ, ret)
        If (on_error(ret)) Go To 130
        dws(i, nnn) = realo
      End Do
      Go To 100

    End If

!   ----------------------------------------------------------------------


!   treat ISO lines

    If (key=='ISO') Then

      Call ffrkey(key, nc, ret)
      If (on_error(ret)) Go To 130

      Call ffrkey(key2, nc, ret)
      If (on_error(ret)) Go To 130
      q = (labl(nlevl(nnn))==key .And. labl(nlevu(nnn))==key2)

      If (.Not. q) Go To 100
      Call ffrnum(realo, integ, ret)
      waveair = realo
      If (srange1) Then
!       convert air to vacuum: ISO wavelength in air (CHECK)
        xkw = 1.D4/xlamb(nnn)
        xn = 1.D0 + 1.D-7*(643.28D0+294981.D0/(146.D0-xkw**2)+2554.D0/(41.D0- &
          xkw**2))
        waveair = xn*waveair !      now, waveair is in vacuum
      End If

      If (abs(waveair-xlamb(nnn))>0.1) Print *, ' CHECK WAVELENGTHS (ISO)! ', &
        waveair, ' ', xlamb(nnn)

!     take value provided (either air or vacuum)

      xlamb(nnn) = waveair
      Call ffrnum(realo, integ, ret)
      osc = realo
      If (abs(log10(osc/flu(nnn)))>0.1) Print *, &
        ' CHECK OSCILLATOR STRENGTH (ISO)! ', osc, ' ', flu(nnn)

!     take value provided

      flu(nnn) = osc

!     read and write data in following format
!     ----OSC,DLP,DENS,NTS,NFS,
!     ----TS(NTS),WS(NTS),DS(NTS),AS(NTS),
!     ----(DLF,FORB(NTS)) FOR EACH NFS

      j = 1
      data(j) = osc
      j = j + 1
      Call ffrnum(realo, integ, ret)
      data(j) = realo !             DLP
      j = j + 1
      Call ffrnum(realo, integ, ret)
      data(j) = realo !             DENS
      j = j + 1
      Call ffrnum(realo, integ, ret)
      ntsiso = integ
      data(j) = integ !             NTS
      j = j + 1
      Call ffrnum(realo, integ, ret)
      nfsiso = integ
      data(j) = integ !             NFS

!     ----NDATA=5+4*NTS+NFS*(1+NTS)

      ndata = 5 + 4*ntsiso + nfsiso*(1+ntsiso)
      If (ndata>ndatamax) Stop ' NDATA > NDATAMAX'

      Do i = 1, 4*ntsiso
        j = j + 1
        Call ffrnum(realo, integ, ret)
        data(j) = realo !           TS,WS,DS,AS
      End Do

      Do i = 1, nfsiso*(1+ntsiso)
        j = j + 1
        Call ffrnum(realo, integ, ret)
        data(j) = realo !           (DLF,FORB) for each forbidden component
      End Do

      nwsn = nws(nnn)
      If (dws(1,nnn)*dws(nwsn,nnn)>=0.D0) Stop ' ERROR IN ISO-FREQ (QHALF)'
      qhalf(nnn) = .False.
      nnat = le(nlevl(nnn))
      peso = weight(nnat)
      xl0 = xlamb(nnn)
      nts(nnn) = 6
      nes(nnn) = 11
      nps = nws(nnn)*nts(nnn)*nes(nnn)
      If (nps>maxps) Stop 'TOO MANY POINTS IN PROFILE'
      Do ijk = 1, 6
        ts(ijk, nnn) = ttts(ijk)
      End Do

      Do ijk = 1, 11
        es(ijk, nnn) = eees(ijk)
      End Do

      nin = 0
jelop: Do je = 1, 11
        Do jt = 1, 6
          tt = 10.D0**ts(jt, nnn)
          ee = 10.D0**es(je, nnn)
          Call iso(data, ndata, tt, ee, peso, vturb, xl0, dws(1:nwsn,nnn), &
            nwsn, sq)
          Do ijk = 1, nwsn
            nin = nin + 1
            p(nin, nnn) = log10(sq(ijk))
          End Do
        End Do
      End Do jelop
      If (nin/=nps) Stop ' ERROR IN NPS(ISO)'

      done(nnn) = .True.
      If (nnn/=nb) Then
        Cycle
      Else
        Exit
      End If

    End If

!   ----------------------------------------------------------------------


!   treat ISO1 lines

    If (key=='ISO1') Then

      Call ffrkey(key, nc, ret)
      If (on_error(ret)) Go To 130

      Call ffrkey(key2, nc, ret)
      If (on_error(ret)) Go To 130
      q = (labl(nlevl(nnn))==key .And. labl(nlevu(nnn))==key2)

      If (.Not. q) Go To 100
      Call ffrnum(realo, integ, ret)
      waveair = realo
      If (srange1) Then
!       convert air to vacuum: ISO1 wavelength in air (CHECK)
        xkw = 1.D4/xlamb(nnn)
        xn = 1.D0 + 1.D-7*(643.28D0+294981.D0/(146.D0-xkw**2)+2554.D0/(41.D0- &
          xkw**2))
        waveair = xn*waveair !      now, waveair is in vacuum
      End If

      If (abs(waveair-xlamb(nnn))>0.1) Print *, ' CHECK WAVELENGTHS (ISO1)! ', &
        waveair, ' ', xlamb(nnn)

!     take value provided (either air or vacuum)

      xlamb(nnn) = waveair
      Call ffrnum(realo, integ, ret)
      osc = realo
      If (abs(log10(osc/flu(nnn)))>0.1) Print *, &
        ' CHECK OSCILLATOR STRENGTH (ISO1)! ', osc, ' ', flu(nnn)

!     take value provided

      flu(nnn) = osc

!     read and write data in following format
!     ----OSC,DLP,DENS,NTS, &
!     ----TS(NTS),WS(NTS),DS(NTS),   electrons
!     ----WS1(NTS),DS1(NTS),WS2(NTS),DS2(NTS) protons, HeII

!     NOTE: tables contain FULL half-widths, in contrast to ISO tables
!     provided by Griem

      j = 1
      data(j) = osc
      j = j + 1
      Call ffrnum(realo, integ, ret)
      data(j) = realo !             DLP
      j = j + 1
      Call ffrnum(realo, integ, ret)
      data(j) = realo !             DENS
      j = j + 1
      Call ffrnum(realo, integ, ret)
      ntsiso = integ
      data(j) = integ !             NTS
      j = j + 1
      Call ffrnum(realo, integ, ret)
      data(j) = realo !             CRIT


!     ----NDATA=5+NTS+6*NTS

      ndata = 5 + 7*ntsiso
      If (ndata>ndatamax) Stop ' NDATA > NDATAMAX'

      Do i = 1, 7*ntsiso
        j = j + 1
        Call ffrnum(realo, integ, ret)
        data(j) = realo !           TS,WS,DS,WS1,DS1,WS2,DS2
      End Do

      nwsn = nws(nnn)
      If (dws(1,nnn)*dws(nwsn,nnn)>=0.D0) Stop ' ERROR IN ISO-FREQ (QHALF)'
      qhalf(nnn) = .False.
      nnat = le(nlevl(nnn))
      peso = weight(nnat)
      xl0 = xlamb(nnn)
      nts(nnn) = 6
      nes(nnn) = 11
      nps = nws(nnn)*nts(nnn)*nes(nnn)
      If (nps>maxps) Stop 'TOO MANY POINTS IN PROFILE'
      Do ijk = 1, 6
        ts(ijk, nnn) = ttts(ijk)
      End Do

      Do ijk = 1, 11
        es(ijk, nnn) = eees(ijk)
      End Do

      nin = 0
jelop1: Do je = 1, 11
        Do jt = 1, 6
          tt = 10.D0**ts(jt, nnn)
          ee = 10.D0**es(je, nnn)
          Call iso1(data, ndata, tt, ee, peso, vturb, xl0, dws(1:nwsn,nnn), &
            nwsn, sq, yhe)
          Do ijk = 1, nwsn
            nin = nin + 1
            p(nin, nnn) = log10(sq(ijk))
          End Do
        End Do
      End Do jelop1
      If (nin/=nps) Stop ' ERROR IN NPS(ISO)'

      done(nnn) = .True.
      If (nnn/=nb) Then
        Cycle
      Else
        Exit
      End If

    End If

!   ----------------------------------------------------------------------

    llchar = labl(nlevl(nnn))
    luchar = labl(nlevu(nnn))
    If (llchar(1:2)=='H1') Then
      Read (llchar(3:4), '(I2)') ll
      Read (luchar(3:4), '(I2)') lu
    End If

!   "normal" lines

!   in case, use Balmer-lines from Lemke
    If (llchar=='H12   ' .And. balmer_lemke) Then
      qhalf(nnn) = .True.
      qlog = .False.
      Call vcs_balmer(xlamb(nnn), ll, lu, nws(nnn), dws(:,nnn), nts(nnn), &
        ts(:,nnn), nes(nnn), es(:,nnn), nps, p(:,nnn))

      If (vturb/=0.) Then
        nwsn = nws(nnn)
        nin = 0
        Do i = 1, nts(nnn)*nes(nnn)
          Do j = 1, nwsn
            nin = nin + 1
            prof(j) = p(nin, nnn)
          End Do

          Call convvturb(prof(1:nwsn), dws(1:nwsn,nnn), nwsn, qlog, &
            qhalf(nnn), i, vturb, xlamb(nnn))
          nin = nin - nwsn
          Do j = 1, nwsn
            nin = nin + 1
            p(nin, nnn) = prof(j)
          End Do
        End Do
        If (nin/=nps) Stop ' ERROR IN NPS(CONVVTURB)'
      End If

!     in case, use Paschen/Brackett-lines from Lemke

    Else If ((llchar=='H13   ' .Or. llchar=='H14   ') .And. pb_lemke .And. &
        lu<=10) Then
      qhalf(nnn) = .True.
      qlog = .False.
      Call vcs_pb(xlamb(nnn), ll, lu, nws(nnn), dws(:,nnn), nts(nnn), &
        ts(:,nnn), nes(nnn), es(:,nnn), nps, p(:,nnn))

      If (vturb/=0.) Then
        nwsn = nws(nnn)
        nin = 0
        Do i = 1, nts(nnn)*nes(nnn)
          Do j = 1, nwsn
            nin = nin + 1
            prof(j) = p(nin, nnn)
          End Do

          Call convvturb(prof(1:nwsn), dws(1:nwsn,nnn), nwsn, qlog, &
            qhalf(nnn), i, vturb, xlamb(nnn))
          nin = nin - nwsn
          Do j = 1, nwsn
            nin = nin + 1
            p(nin, nnn) = prof(j)
          End Do
        End Do
        If (nin/=nps) Stop ' ERROR IN NPS(CONVVTURB)'
      End If

    Else

      Call ffrkey(key, nc, ret)
      If (on_error(ret)) Go To 130

      Call ffrkey(key2, nc, ret)
      If (on_error(ret)) Go To 130

      q = (labl(nlevl(nnn))==key .And. labl(nlevu(nnn))==key2)

      If (.Not. q) Go To 100

      Call ffrnum(realo, integ, ret)
      If (on_error(ret)) Go To 130
      mfull = integ
      qhalf(nnn) = mfull == 0
      Call ffrnum(realo, integ, ret)
      If (on_error(ret)) Go To 130
!     ---------- number of wavelengths
      If (integ>maxw) Stop 'TOO MANY LAMBDAS'
      nws(nnn) = integ
      Do i = 1, integ
        Call ffrnum(realo, integ, ret)
        If (on_error(ret)) Go To 130
        dws(i, nnn) = realo
      End Do

      Call ffrnum(realo, integ, ret)
      If (on_error(ret)) Go To 130
!     ---------- number of temperatures
      If (integ>maxt) Stop 'TOO MANY TEMPS'
      nts(nnn) = integ
      Do i = 1, integ
        Call ffrnum(realo, integ, ret)
        If (on_error(ret)) Go To 130
        ts(i, nnn) = realo
      End Do
      Call ffrnum(realo, integ, ret)
      If (on_error(ret)) Go To 130
!     ---------- number of electron densities
      If (integ>maxne) Stop 'TOO MANY ELECTRON DENS.'
      nes(nnn) = integ
      Do i = 1, integ
        Call ffrnum(realo, integ, ret)
        If (on_error(ret)) Go To 130
        es(i, nnn) = realo
      End Do
      Call ffrnum(realo, integ, ret)
      If (on_error(ret)) Go To 130
      fac = realo
      Call ffrnum(realo, integ, ret)
      If (on_error(ret)) Go To 130
      qlog = integ == 1
!     ---------- number of profile points
      nps = nws(nnn)*nts(nnn)*nes(nnn)
      If (nps>maxps) Stop 'TOO MANY POINTS IN PROFILE'

      If (vturb==0.) Then
        Do i = 1, nps
          Call ffrnum(realo, integ, ret)
          If (on_error(ret)) Go To 130
          p(i, nnn) = realo*fac
          If (qlog) p(i, nnn) = log10(p(i,nnn))
        End Do
      Else
        nwsn = nws(nnn)
        nin = 0
        Do i = 1, nts(nnn)*nes(nnn)
          Do j = 1, nwsn
            Call ffrnum(realo, integ, ret)
            If (on_error(ret)) Go To 130
            prof(j) = realo*fac
          End Do
          Call convvturb(prof(1:nwsn), dws(1:nwsn,nnn), nwsn, qlog, &
            qhalf(nnn), i, vturb, xlamb(nnn))
          If (qlog) Then
            Do j = 1, nwsn
              prof(j) = log10(prof(j))
            End Do
          End If
          Do j = 1, nwsn
            nin = nin + 1
            p(nin, nnn) = prof(j)
          End Do
        End Do
        If (nin/=nps) Stop ' ERROR IN NPS(CONVVTURB)'
      End If

    End If

    done(nnn) = .True.
    If (nnn/=nb) Then
      Cycle
    Else
      Exit
    End If

    Stop ' ERROR IN PATH'


!   ----------------------------------------------------------------------

!   ---- stark profile(s) not in tables. we calculate them in Griem-approx,
!   if hydrogen or HeII

110 Continue

    If (done(nnn)) Stop ' WRONG PATH'

    Print *, 'STARK PROFILES NOT PRESENT; COMPONENT:', nnn, ' , GRIEM USED'
    key = labl(nlevl(nnn))
    If (key(:3)=='HE2' .Or. key(:2)=='H1') Then
      qhalf(nnn) = .True.
      nnat = le(nlevl(nnn))
      peso = weight(nnat)
      xl0 = xlamb(nnn)
      nn1 = nppal(key)
      Print *, 'H or HeII LINE: GRIEM USED  ', xl0, ' ', key, ' ', &
        labl(nlevu(nnn))
      key = labl(nlevu(nnn))
      nn2 = nppal(key)
      aa = dble(nn1)
      bb = dble(nn2)
      nts(nnn) = 6
      nes(nnn) = 11
      nws(nnn) = 35
      nin = 0
      Do ijk = 1, 6
        ts(ijk, nnn) = ttts(ijk)
      End Do

      Do ijk = 1, 11
        es(ijk, nnn) = eees(ijk)
      End Do

jeloop: Do je = 1, 11
        Do jt = 1, 6
          tt = 10.D0**ts(jt, nnn)
          ee = 10.D0**es(je, nnn)
          Call griem(tt, ee, sq1, peso, xl0, aa, bb, ddws, vturb)
          Do ijk = 1, 35
            sq(ijk) = sq1(ijk)
            If (sq(ijk)<0.) Stop ' NEG. PROFILE IN GRIEM'
            nin = nin + 1
            dws(ijk, nnn) = ddws(ijk)
            p(nin, nnn) = log10(sq(ijk))
          End Do
        End Do
      End Do jeloop
    Else
      Print *, 'STARK CALCULATION NOT IMPLEMENTED'
      Print *, 'LINE:', labl(nlevl(nnn)), '-', labl(nlevu(nnn)), ' ', &
        xlamb(nnn)
      If (srange1) Then
        Print *, 'NSTARK RESET TO 0'
        nstark(nnn) = 0
      Else
        Stop ' CHANGE NSTARK IN FORMAL_INPUT'
      End If
    End If

    done(nnn) = .True.
    If (nnn/=nb) Then
      Cycle
    Else
      Exit
    End If

  End Do big_loop

  Close (1)

! ---- output file is written

  Open (1, File=trim(file)//'/STARK.BIN', Status='UNKNOWN', &
    Form='UNFORMATTED')
  Rewind 1

  Do i = 1, nb
    If (nstark(i)/=1) Cycle
    If (.Not. done(i)) Stop ' ERROR IN PRESTARK, BEFORE WRITING'
    nn = nws(i)*nts(i)*nes(i)
    Write (1) nlevl(i), nlevu(i), qhalf(i), nws(i), (dws(j,i), j=1, nws(i)), &
      nts(i), (ts(j,i), j=1, nts(i)), nes(i), (es(j,i), j=1, nes(i)), nn, &
      (p(j,i), j=1, nn)
  End Do
  Close (1)

  Deallocate (dws, es, p, ts)
  Deallocate (nes, nts, nws)
  Deallocate (done, qhalf)

  Return

! ---- error in input unit
120 Continue

  Stop 'ERROR IN INPUT UNIT'

! error handling

130 Select Case (ret)
  Case ('RET1') !                   end of file exit
    Write (*, Fmt='(A)') ' END OF FILE '
    Stop
  Case ('RET2') !                   error exit
    Write (*, Fmt='(A)') 'ERROR IN FFR SUBROUTINES '
    Stop
  Case Default
    Stop ' WRONG ERROR CONDITION IN PRESTARK'
  End Select

End Subroutine

!ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc

Subroutine griem(tt, ee, sq, weight, wavc, aa, bb, dws, vturb)

  Use :: nlte_type
  Use :: nlte_dim
  Use :: fund_const, Only: xmh => amh, akb, pi
  Use :: preformal_var, Only: ns, as, odop, ps, zz, ss, sx

  Implicit None

! -----create a profile table using the griem theory as used by a and m
! -----in their ap j s24 paper
! ..
! .. scalar arguments ..
  Real (dp) :: aa, bb, ee, tt, wavc, weight, vturb
! ..
! .. array arguments ..
  Real (dp) :: dws(35), sq(35)
! ..
! .. local scalars ..
  Real (dp) :: asq, beta, betap, betaw, bsq, clight, conkb, cue, divkab, &
    dopkab, dw, f0, fac, fpw, gam0, gama, gaml, gamw, oba, rft, sre, srt, sss, &
    summ, vmot, wsq, xkabf0
  Integer (i4b) :: i, jm, jp, nbet
! ..
! .. local arrays ..
  Real (dp) :: bet(35), ddws(35)
! ..
! .. external functions ..
  Real (dp) :: dconv, tbg, th
  External :: dconv, tbg, th
! ..
! .. intrinsic functions ..
! ..
! .. save statement ..
  Save
! ..
! .. data statements ..
  Data ddws/0.D0, 1.5849D-3, 2.5119D-3, 3.9811D-3, 6.3096D-3, 1.D-2, &
    1.5849D-2, 2.5119D-2, 3.9811D-2, 6.3096D-2, 1.D-1, 1.585D-1, 1.995D-1, &
    2.512D-1, 3.162D-1, 3.981D-1, 5.012D-1, 6.310D-1, 7.943D-1, 1.00D0, &
    1.259D0, 1.585D0, 1.995D0, 2.512D0, 3.981D0, 6.310D0, 1.0D1, 1.585D1, &
    2.512D1, 3.891D1, 6.310D1, 1.0D2, 1.585D2, 2.512D2, 3.981D2/


  Data clight/2.997925D18/
  Data nbet/35/
! ..

  srt = sqrt(tt)

! ---- calculation of clight/v_thermal
  vmot = clight*1.D-8/sqrt(2.D0*akb*tt/xmh/weight+vturb**2)

! -----doppler width (a)
  rft = wavc/clight
  odop = vmot/wavc

! -----upper and lower level parameters
  asq = aa*aa
  bsq = bb*bb
  wsq = wavc*wavc
  oba = 1.0D0/(bsq-asq)

! -----normal field strength
  sre = sqrt(ee)
  cue = ee**0.3333333333D0
  f0 = 1.25D-9*cue*cue

! -----ratio stark/doppler width
  dopkab = 5.5D-5*(asq*bsq)**2*oba*f0*odop/zz**5

! -----1/stark width     (a)**-1
  divkab = odop/dopkab

! ---- calculation  of beta grid
  xkabf0 = dopkab/odop

  Do i = 1, nbet
    bet(i) = ddws(i)/xkabf0
  End Do

  betaw = 6.951D-9*wsq*zz*tt*divkab/bsq
  betap = 2.995D-15*wsq*sre*divkab
  betap = min(betap, .99D0*betaw)
  gam0 = 3.78D-5*(bsq-3.D0*aa)*bsq*oba*cue*log10(4.D6*zz*tt/bsq/sre)/srt/zz
  gaml = log(betaw/betap)
  gamw = 4.712389D0/sqrt(betaw)
  fpw = gam0/gaml

! -----initialize gama for zero,betap
  gama = gamw + gam0

! -----generate s(beta) on a fixed basis
  jm = nbet
  jp = jm

  Do i = 1, nbet
    beta = bet(i)
    If (beta<betap) Go To 100
    If (beta>betaw) Go To 120

!   -----   assign gama for betap,betaw
    gama = gamw + fpw*log(betaw/beta)

!   -----   look for regions of the gama,beta plane where fast methods are
!   -----   this is not just to save time as the t(beta,gama) analysis
!   -----   is numerically unstable in these same regions
100 Continue
    If (gama>25.) Go To 130
    If (gama<.01D0*beta) Go To 110
    If (gama<.01) Go To 110

!   -----   full case
    sss = tbg(beta, gama)
    Go To 140

!   -----   beta gt 100*gama for beta lt 1
110 Continue
    sss = th(beta)
    If (beta<1.0D0) Go To 140

!   -----   gama lt  .01  or beta gt 100*gama for beta gt 1
    sss = sss + gama/(pi*beta*beta)
    Go To 140

!   -----   beta gt betaw
120 Continue
    sss = 2.0D0*th(beta)
    Go To 140

!   -----   gama gt 25.
130 Continue
    sss = gama/(pi*(gama*gama+beta*beta))

!   -----   fill the symmetric ss,sx vectors
140 Continue
    ss(jm) = sss
    sx(jm) = -beta*dopkab
    ss(jp) = sss
    sx(jp) = beta*dopkab
    jm = jm - 1
    jp = jp + 1
  End Do


! -----make the asymtotic power law constants
  ns = 2*nbet - 1
  ps = log(ss(ns)/ss(ns-1))/log(sx(ns)/sx(ns-1))
  conkb = -2.0D0
  ps = min(ps, conkb)
  as = ss(ns)/sx(ns)**ps

! -----normalize to unit frequency integral
  summ = 0.0D0
  Do i = 2, ns
    summ = summ + (ss(i)+ss(i-1))*(sx(i)-sx(i-1))
  End Do

  fac = 2.0D0*vmot*rft/summ

! -----convolve with doppler function
  Do i = 1, 35
    dw = ddws(i)
    sq(i) = fac*dconv(dw)
    dws(i) = ddws(i)
  End Do

  Return


End Subroutine

!-----------------------------------------------------------------------

Function tbg(bet, gam)

  Use :: nlte_type
  Use :: nlte_dim
  Use :: fund_const, Only: pi
  Implicit None

! .. scalar arguments ..
  Real (dp) :: bet, gam, tbg
! ..
! .. local scalars ..
  Real (dp) :: t
! ..
! .. external functions ..
  Real (dp) :: tf, tg, th
  External :: tf, tg, th
! ..
! .. save statement ..
  Save
! ..

  If (gam>1.D-2) Go To 100

! small gamma
  tbg = th(bet)
  If (bet>1.) tbg = tbg + gam/pi/bet**2
  Return

! normal case
100 Continue
  If (bet/gam>1.D2) Go To 110
  t = gam*(tf(bet,gam)+tf(-bet,gam)+tg(bet,gam)+tg(-bet,gam))/pi
  tbg = t
  Return

110 Continue
  tbg = th(bet) + gam/pi/bet**2

  Return

End Function

!-----------------------------------------------------------------------

Function tg(bet, gam)

  Use :: nlte_type
  Use :: nlte_dim
  Implicit None

! .. scalar arguments ..
  Real (dp) :: bet, gam, tg
! ..
! .. local scalars ..
  Real (dp) :: b, bet2, c, c4, c5, c6, cc, cosa2, d1, d2, d3, d4, d5, d6, d7, &
    d8, d9, g2, p, p2, pi2, q, q2, r, sina, sina2, sum1, sum2, sum3, sum4, x1, &
    x2, x3, x4, x5, y1, y2
! ..
! .. intrinsic functions ..
! ..
! .. save statement ..
  Save
! ..
! .. data statements ..
  Data c4, c5, c6/1.5, 3.4636008D1, -1.3253986D2/
  Data pi2/1.57079632679489D0/
! ..

  g2 = gam*gam
  b = 2.0D0*bet
  c = bet*bet + g2
  bet2 = bet*bet
  cc = c*c
  r = sqrt(c)
  q = sqrt(r)
  q2 = q*q
  p = 1.D0/q
  p2 = p*p
  sina = gam/r
  cosa2 = sqrt(0.5D0*(1.D0-bet/r))
  sina2 = sqrt(0.5D0*(1.D0+bet/r))
  x5 = bet + 4.D0
  y2 = log((x5*x5+g2)/16.D0)
  If (x5/gam>1.0) Go To 100
  y1 = pi2 - atan(x5/gam)

  Go To 110

100 Continue
  y1 = atan(gam/x5)

110 Continue
  d2 = c/192.D0
  d3 = -b/32.D0
  d4 = 0.25D0*(3.D0*bet2-g2)/c
  d5 = b*(g2-bet2)/cc
  d6 = (bet2*(bet2-6.0D0*g2)+g2*g2)/cc
  sum1 = ((d6*y1)/gam+(d4+d3)) + d2
  sum1 = (sum1+d5*y2)*c5/cc
  d1 = c/1024.D0
  d2 = -b/192.D0
  d3 = (3.0D0*bet2-g2)/(32.D0*c)
  d4 = bet*(g2-bet2)/cc
  d5 = 0.5D0*(g2*g2+5.D0*bet2*(bet2-2.0D0*g2))/(c*cc)
  d6 = -bet*(bet2*bet2+5.0D0*g2*(g2-2.0D0*bet2))/(cc*c)
  sum2 = (d5*y2+d3) + d1
  sum2 = (sum2+(d2+d4+(d6*y1)/gam))*c6/cc
  d7 = c4/c
  d8 = d7*(b*b/c-1.D0)/(2.D0*q*q2*sina)
  d9 = d7*(b/(c*c))/(2.D0*p*p2*sina)
  x1 = (4.0D0-q2)/(4.0D0*q*sina2)
  x2 = (q*(q+4.D0*cosa2)+4.D0)/(q*(q-4.D0*cosa2)+4.D0)
  x3 = (0.25D0-p2)/(p*sina2)
  x4 = (p*(p+cosa2)+0.25D0)/(p*(p-cosa2)+0.25D0)
  If (x1>1.) Go To 120
  y1 = pi2 - atan(x1)

  Go To 130

120 Continue
  y1 = atan(1.D0/x1)

130 Continue
  If (x3>-1.) Go To 140
  y2 = -atan(1.D0/x3)

  Go To 150

140 Continue
  y2 = pi2 + atan(x3)

150 Continue
  sum3 = d8*(2.D0*cosa2*y1-sina2*log(x2))
  sum4 = d9*(2.D0*cosa2*y2+sina2*log(x4))
  tg = (sum4+d7*(1.D0/12.D0-b/c)) + sum3
  tg = tg + sum1 + sum2

  Return

End Function

!-----------------------------------------------------------------------

Function tf(b, g)

  Use :: nlte_type
  Use :: nlte_dim
  Implicit None

! .. scalar arguments ..
  Real (dp) :: b, g, tf
! ..
! .. local scalars ..
  Real (dp) :: c0, c1, c2, c3, d, d1, d2, d3, g2, x1
! ..
! .. intrinsic functions ..
! ..
! .. save statement ..
  Save
! ..
! .. data statements ..
  Data c0, c1, c2, c3/1.0007744D-1, 4.93208719D-3, -7.09873526D-3, &
    7.11559325D-4/
! ..
  d3 = c3
  d2 = c2 - 3.0D0*b*c3
  d1 = (3.0D0*c3*b-2.0D0*c2)*b + c1
  d = ((-c3*b+c2)*b-c1)*b + c0
  g2 = g*g
  x1 = b + 4.D0
  tf = 4.0D0*d2 + 4.0D0*d3*(b+2.0D0) + 0.5D0*(d1-g2*d3)*log((x1*x1+g2)/(b*b+g2 &
    )) + (d-g2*d2)*(atan(x1/g)-atan(b/g))/g

  Return

End Function

!-----------------------------------------------------------------------

Function dconv(dlam)

  Use :: nlte_type
  Use :: nlte_dim
  Use :: preformal_var, Only: ns, odop, ss, sx
  Implicit None

! -----
! -----convolution of gaussian profile with s function  6 aug 77
! -----
! .. scalar arguments ..
  Real (dp) :: dlam, dconv
! ..
! .. local scalars ..
  Real (dp) :: con, cr, db, dx, gr, half, range, srtpi, terx, tex, x1, xkb, &
    xn, zero
  Integer (i4b) :: i, j, jm, n0
! ..
! .. local arrays ..
  Real (dp) :: erx(69), ex(69), x(69)
! ..
! .. external functions ..
  Real (dp) :: asint, err
  External :: asint, err
! ..
! .. intrinsic functions ..
! ..
! .. save statement ..
  Save
! ..
! .. data statements ..
  Data half/0.5D0/, zero/0.0D0/, srtpi/5.6418958D-1/
  Data range/6.0D0/
! ..

  db = abs(dlam)*odop
  x1 = sx(1) + db
  If (x1<range) Go To 100

! -----asint -range,range
  dconv = asint(-range, range, db)
  Return

100 Continue
  dconv = zero

! -----asint -range,x1
  If (x1>-range) dconv = dconv + asint(-range, x1, db)

! -----asint xn,range
  xn = sx(ns) + db
  If (xn<range) dconv = dconv + asint(xn, range, db)

! -----set up x in the gaussian frame
! -----store all exponential and error functions once
  n0 = 2
  con = half*srtpi

  Do i = 1, ns
    xkb = sx(i) + db
    x(i) = xkb
    If (xkb<-range) n0 = i + 1
    ex(i) = exp(-xkb*xkb)
    erx(i) = err(xkb, ex(i))
  End Do

! -----by s segment
  terx = zero
  tex = zero

  Do j = n0, ns
    jm = j - 1
    dx = x(j) - x(jm)
    If (dx<0.1) Go To 110
    gr = (ss(j)-ss(jm))/dx
    cr = ss(j) - gr*x(j)
    terx = terx + cr*(erx(j)-erx(jm))
    tex = tex + gr*(ex(j)-ex(jm))

    Go To 120
110 Continue
    dconv = dconv + con*(ss(j)*ex(j)+ss(jm)*ex(jm))*dx
120 Continue
    If (x(j)>range) Exit
  End Do

  dconv = dconv + half*(terx-srtpi*tex)

  Return

End Function

!-----------------------------------------------------------------------

Function asint(x0, x1, db)

  Use :: nlte_type
  Use :: nlte_dim
  Use :: preformal_var, Only: as, ps
  Implicit None

! -----
! -----integral over asymtotic region by simpsons rule
! -----
! .. scalar arguments ..
  Real (dp) :: db, x0, x1, asint
! ..
! .. local scalars ..
  Real (dp) :: c23, dx, h, half, one, srtpi, step, two, x
  Integer (i4b) :: i, n
! ..
! .. intrinsic functions ..
! ..
! .. statement functions ..
  Real (dp) :: f
! ..
! .. save statement ..
  Save
! ..
! .. data statements ..
  Data srtpi/5.6418958D-1/, c23/6.666666667D-1/
  Data half, one, two/0.5D0, 1.0D0, 2.0D0/
  Data step/0.2D0/
! ..
! .. statement function definitions ..
  f(x) = exp(-x*x)*abs(x-db)**ps
! ..

! -----chose h to be sept
  n = int((x1-x0)/step+one)
  dx = (x1-x0)/dble(n)
  h = half*dx
  x = x0 - h
  asint = f(x0)

  Do i = 1, n
    x = x + dx
    asint = asint + two*f(x) + f(x+h)
  End Do

  asint = asint*as*srtpi*h*c23

  Return

End Function

!-----------------------------------------------------------------------

Function err(x, ex)

  Use :: nlte_type
  Use :: nlte_dim
  Implicit None

! -----error function
! .. scalar arguments ..
  Real (dp) :: ex, x, err
! ..
! .. local scalars ..
  Real (dp) :: one, t
! ..
! .. intrinsic functions ..
! ..
! .. save statement ..
  Save
! ..

  one = 1.0D0
  t = abs(x)
  If (t>6.) Go To 100
  t = 1.D0/(1.D0+0.3275911D0*t)
  err = ((((1.061405429D0*t-1.453152027D0)*t+1.421413741D0)*t-0.284496736D0)*t &
    +0.254829592D0)*t
  err = sign(one-ex*err, x)

  Return

100 Continue
  err = sign(one, x)

  Return

End Function

!-----------------------------------------------------------------------

Function th(x)

  Use :: nlte_type
  Use :: nlte_dim
  Implicit None

! -----griem microfield fudge - normalized to unit b integral
! .. scalar arguments ..
  Real (dp) :: x, th
! ..
! .. local scalars ..
  Real (dp) :: b, b2, c0, c1, c2, c3, c4, c5, c6, srt
! ..
! .. intrinsic functions ..
! ..
! .. save statement ..
  Save
! ..
! .. data statements ..
  Data c0, c1, c2, c3/1.0007744D-1, 4.93208719D-3, -7.09873526D-3, &
    7.11559325D-4/
  Data c4, c5, c6/1.5, 3.4636008D1, -1.3253986D2/
! ..

  b = abs(x)
  If (b>4.) Go To 100
  th = ((c3*b+c2)*b+c1)*b + c0

  Return

100 Continue
  srt = sqrt(b)
  b2 = b*b
  th = ((c6/b+c5)/b2+c4/srt)/b2

  Return

End Function

!-----------------------------------------------------------------------

Function nppal(lab)

  Use :: nlte_type
  Use :: nlte_dim
  Use :: preformal_var, Only: zz
  Implicit None

! ---- this function determines the principal quantum number, n,
! ---- of a certain level of an hydrogenic ion, if the level's label
! ---- is written in the "normal detail way"

! ..
  Integer (i4b) :: nppal

! .. scalar arguments ..
  Character :: lab*6
! ..
! .. local scalars ..
  Integer (i4b) :: ij, j
  Character :: car*1, at*2, cara*2
! ..
! .. save statement ..
  Save
! ..

  j = 1

100 Continue
  car = lab(j:j)
  If (car=='0' .Or. car=='1' .Or. car=='2' .Or. car=='3' .Or. car=='4' .Or. &
    car=='5' .Or. car=='6' .Or. car=='7' .Or. car=='8' .Or. car=='9') Then
!   ----------- atom determination for atomic charge
    at = lab(1:j-1)

    If (at=='H') Then
      zz = 1.D0
    Else If (at=='HE') Then
      zz = 2.D0
    Else
      Print *, 'ATOM ', at, ' NOT IMPLEMENTED - NPPAL'
      Stop
    End If

    j = j + 1
    ij = j

110 Continue
    car = lab(j:j)
    If (car=='0' .Or. car=='1' .Or. car=='2' .Or. car=='3' .Or. car=='4' .Or. &
      car=='5' .Or. car=='6' .Or. car=='7' .Or. car=='8' .Or. car=='9') Then
      j = j + 1
      Go To 110
    Else
      j = j - 1
      cara = lab(ij:j)
      Read (cara, Fmt='(I2)') nppal
      If (nppal==0) nppal = 10
    End If

  Else
    j = j + 1
    Go To 100
  End If

  Return
End Function

!-----------------------------------------------------------------------

Subroutine convvturb(prof, ddws, ns, qlog, qhalf, istart, vturb, lambda0)

! convolution with vturb, assuming profile is linear between grid points.
! for istart=1, convolution matrix will be calculated, for istart > 1
! "only" applied


  Use :: nlte_type
  Use :: nlte_dim
  Use :: fund_const, Only: clight
  Implicit None

! Input data
  Integer (i4b) :: ns, istart
  Real (dp) :: prof(ns), ddws(ns), vturb, lambda0
  Logical :: qlog, qhalf

  Integer (i4b), Parameter :: nsmax = id_maxww, nsmax2 = 2*nsmax - 1
  Real (dp), Parameter :: srtpi1 = 0.5641895835D0

  Integer (i4b) :: i, j, jmin, jmax, j1, j2, ns2
  Real (dp) :: const, err, x0, w
  Real (dp) :: q, q1, dx, dx1, dx2, deltaexp, deltaerf, x1, x2, exp1, exp2, &
    erfc1, erfc2
  Real (dp) :: dw(nsmax2), x(0:nsmax2+1), p(0:nsmax2+1), a(nsmax, 0:nsmax2+1), &
    pnew(nsmax)

  Save :: a

! check whether input from blue to red as it is assumed
  If (ddws(1)>ddws(ns)) Stop ' INPUT FROM RED TO BLUE IN CONVVTURB'

! change order from red to blue and complete variables if QHALF

  ns2 = ns

  If (qhalf) Then
    ns2 = 2*ns - 1
    Do i = 1, ns
      j = i + ns - 1
      dw(i) = ddws(ns+1-i)
      p(i) = prof(ns+1-i)
      dw(j) = -ddws(i)
      p(j) = prof(i)
    End Do
  Else
    Do i = 1, ns2
      dw(i) = ddws(ns2+1-i)
      p(i) = prof(ns2+1-i)
    End Do
  End If

  If (.Not. qlog) Then
    Do i = 1, ns2
      p(i) = 10.**p(i)
    End Do
  End If

! calculate matrix once

  If (istart==1) Then

    a = 0.
!   conversion to Doppler units (with respect to vturb)

    const = clight/(vturb*lambda0)
    Do i = 1, ns2
      x(i) = -dw(i)*const
    End Do
    x(0) = x(1) - 5.01D0
    x(ns2+1) = x(ns2) + 5.1D0

!   append first and last element assuming constant profiles
    p(0) = p(1)
    p(ns2+1) = p(ns2)

    jmin = 0
    jmax = 0

!   set up matrix, calculate full (NS2=NS) or only half (NS2 NE NS) range
big: Do i = 1, ns
      x0 = x(i)

!     find range
      Do j = jmin, ns2
        If (x(j)-x0>-5.D0) Exit
      End Do
      jmin = j - 1

      Do j = jmax, ns2 + 1
        If (x(j)-x0>=5.D0) Exit
      End Do
      jmax = j

!     calculate weights
!     first element
      j1 = jmin
      j2 = j1 + 1
      x1 = x(j1)
      x2 = x(j2)
      dx = x2 - x1
      dx1 = x1 - x0
      dx2 = x2 - x0
      q = -dx1/dx
      q1 = 1.D0 - q
      exp1 = exp(-dx1*dx1)
      exp2 = exp(-dx2*dx2)
      deltaexp = exp1 - exp2
      erfc1 = 1.D0 - err(dx1, exp1)
      erfc2 = 1.D0 - err(dx2, exp2)
      deltaerf = erfc1 - erfc2
      w = q1*deltaerf - srtpi1*deltaexp/dx
      a(i, j1) = .5D0*w

!     middle elements
      Do j = jmin + 1, jmax - 1

!       old values
        w = q*deltaerf + srtpi1*deltaexp/dx

!       new values
        j1 = j
        j2 = j1 + 1
        x1 = x(j1)
        x2 = x(j2)
        dx = x2 - x1
        dx1 = x1 - x0
        dx2 = x2 - x0
        q = -dx1/dx
        q1 = 1.D0 - q
        exp1 = exp(-dx1*dx1)
        exp2 = exp(-dx2*dx2)
        deltaexp = exp1 - exp2
        erfc1 = 1.D0 - err(dx1, exp1)
        erfc2 = 1.D0 - err(dx2, exp2)
        deltaerf = erfc1 - erfc2
        w = w + q1*deltaerf - srtpi1*deltaexp/dx
        a(i, j1) = .5D0*w
      End Do

!     last element
      w = q*deltaerf + srtpi1*deltaexp/dx
      a(i, jmax) = .5D0*w

    End Do big

  End If

! convolution
  pnew(1:ns) = matmul(a(1:ns,0:ns2+1), p(0:ns2+1))

! test
  w = 0.D0
  q = 0.D0
  Do i = 1, ns - 1
    dx = .5*(dw(i+1)-dw(i))
    w = w + (p(i)+p(i+1))*dx
    q = q + (pnew(i)+pnew(i+1))*dx
  End Do

! error in equivalent width
  dx = 1. - q/w
  If (abs(dx)>0.1D0) Print *, ' WARNING!!! ERROR IN CONVOLUTION LARGE: ', dx

! renormalization
  Do i = 1, ns
    pnew(i) = pnew(i)*w/q
  End Do

! reversal of order: either pnew is ordered from -xmax to 0 (symmetric case)
! or from -xmax to +xmax (asymmetric case)

  Do i = 1, ns
    prof(i) = pnew(ns+1-i)
  End Do

  If (.Not. qlog) Then
    Do i = 1, ns
      prof(i) = log10(prof(i))
    End Do
  End If

  Return
End Subroutine

!-----------------------------------------------------------------------

Subroutine iso(data, ndata, tt, ee, weight, vturb, lambda0, dws, nws, sq)

  Use :: nlte_type
  Use :: nlte_dim
  Implicit None


! -----ISOLATED AND MULTI COMPONENT HE I LINE PROFILES
! -----INPUT DATA ARE - OSC,DLP,DENS,NTS,NFS,
! -----                 TS(NTS),WS(NTS),DS(NTS),AS(NTS),
! -----                 (DLF,FORB(NTS)) FOR EACH NFS

! -----Updated/corrected by j. puls 13/04/99
! depends majorly on treatment by Barnard, Cooper, Smith, 1974,
! J.Q.S.R.T. 14, 1025


! -----NDATA=5+4*NTS+NFS*(1+NTS)

  Integer (i4b), Parameter :: maxww = id_maxww
  Real (dp), Parameter :: clight = 2.997925D18

! input values

  Integer (i4b) :: ndata, nws
  Real (dp) :: data(ndata), tt, ee, weight, vturb, lambda0, dws(nws), sq(nws)

! local values
  Integer (i4b) :: i, ib, if, j, k, l, nfs, nts
  Real (dp) :: a, alf, con, d, db, dens, dlp, ft, fq(maxww), p, rba, rbhz, &
    rft, rwt, rhom, sigma, srt, va, vb, vth, vmot, vp, w, wf, wq, wt, x, y

  Real (dp) :: ts(4), ws(4), ds(4), as(4), dls(3), fb(6, 3), fos(3)

! .. external functions ..
  Real (dp) :: voigt
  External :: voigt

! -----frequential quantities

  Do i = 1, nws
    wq = lambda0 + dws(i)
    fq(i) = clight/wq
  End Do

  ft = clight/lambda0
  rft = 1.D0/ft

  srt = sqrt(tt)

  vth = sqrt(1.6504D8*tt/weight+vturb**2)
  vmot = clight*1.D-8/vth

  wt = clight*rft
  rwt = ft/clight

! move DATA to local

  i = 1
  dlp = data(i+1)
  dens = data(i+2)
  nts = data(i+3) + 0.5D0
  nfs = data(i+4) + 1.5D0
  k = i + 4

  Do j = 1, nts
    k = k + 1
    ts(j) = data(k)
  End Do

  Do j = 1, nts
    k = k + 1
    ws(j) = data(k)
  End Do

  Do j = 1, nts
    k = k + 1
    ds(j) = data(k)
  End Do

  Do j = 1, nts
    k = k + 1
    as(j) = data(k)
  End Do

  If (nfs>1) Then
    Do l = 2, nfs
      k = k + 1
      dls(l) = data(k)
      Do j = 1, nts
        k = k + 1
        fb(j, l) = data(k)
      End Do
    End Do
  End If


! RBHZ = 1/deltanu-dop
! RBA  = nu0/vth  (for transformation dlambda -> dnue

  rbhz = vmot*rft
  rba = vmot*rwt

  dlp = dlp*rba
  If (nfs>1) Then
    Do l = 2, nfs
      dls(l) = dls(l)*rba
    End Do
  End If

! -----SET UP INTERPOLATION IN T
  Do j = 2, nts
    if = j
    If (ts(j)>=tt) Exit
  End Do

  ib = if - 1
  wf = (tt-ts(ib))/(ts(if)-ts(ib))

  If (tt<ts(ib)) wf = 0.0D0
  If (tt>ts(if)) wf = 1.0D0

! -----PERTURBER QUANTITIES
  y = ee/dens

! VP=8.78D0/SRT !factor cannot be understood, should be

! factor = [(1./A_perturber)+(1./A_radiator)]^(-1/2), A atomic weight
! A_perturber = hydrogen = 1, A_radiatior = helium = 4


  vp = 1.D0/(sqrt(1.+1./weight)*srt) ! inverse of mean vel., without factors

  rhom = 0.62035D0/ee**0.333333D0 ! mean distance, from 4pi/3 Ne rhom^3=1

! -----IMPACT WIDTH WIDTH (A)
  w = (wf*(ws(if)-ws(ib))+ws(ib))*y

! -----RATIO IMPACT SHIFT/WIDTH
  d = wf*(ds(if)-ds(ib)) + ds(ib)

! -----ION BROADENING PARAMETERS
  alf = (wf*(as(if)-as(ib))+as(ib))*y**0.25D0

! factor accounts for transformation dlambda -> dnu and appropriate units
! factor = sqrt(pi*m_hyd/(8 k_botz)) * clight * 1.e8

  sigma = 2.06D14*w*rhom*vp*rwt*rwt ! sig=w*rhom/v (*c/lam0^2)
  x = alf**0.888889D0/sigma**0.333333D0 ! alf^(8/9) sig^(-1/3) propto w=wi/we

! -----TOTAL WIDTH IN DOPPLER UNITS
  a = w*(1.0D0+1.36D0*x)*rba !      we*(1+wi/we) at line-center
  dls(1) = w*d*(1.0D0+2.36D0*x/abs(d))*rba ! de*(1+di/de) at line-center, no
! Debye shielding
! and di/de=(di/we)/(de/we), BCS 6.6
  fos(1) = 1.0D0

! -----SATELLITE COMPONENTS
  x = fos(1)
  If (nfs>1) Then
    Do l = 2, nfs
      fos(l) = wf*(fb(if,l)-fb(ib,l)) + fb(ib, l)
      x = x + fos(l)
    End Do
  End If

! con includes transformation to Doppler units via RBHZ,
! factor 1/sqrt(pi) from normalization of Voigt profile and
! factor 1/9 from using (8:1) weighted Voigt profiles to account
! for fine-structure of lower component, level 2P3.
! Works also for other cases if separation is set to zero

  con = 6.268773D-2*rbhz/x

! -----COMPUTE PROFILE
  Do j = 1, nws
    db = (ft-fq(j))*rbhz
    p = 0.D0
    Do l = 1, nfs
      va = db - dls(l)
      vb = va - dlp
      p = p + fos(l)*(8.0D0*voigt(a,va)+voigt(a,vb))
    End Do
    sq(j) = con*p
  End Do

  Return
End Subroutine

!-----------------------------------------------------------------------

Subroutine iso1(data, ndata, tt, ee, weight, vturb, lambda0, dws, nws, sq, &
  yhe)

  Use :: nlte_type
  Use :: nlte_dim
  Implicit None


! -----ISOLATE HE I LINE PROFILES according to data from
! Dimitrijevic & Sahal-Brechot.
! So far, proton/HeII contribution only approximate
! from YHe input values and HeI (assumed to be constant
! throughout the atmosphere)
! If HeIII major component (corresp. to HeI = 2), same
! values as for HeII used.


! -----INPUT DATA ARE - OSC,DLP,DENS,NTS,
! -----                 TS(NTS),WS(NTS),DS(NTS),  for electrons
! -----                 WS1(NTS),DS1(NTS),WS2(NTS),DS2(NTS) for protons and
! HeII

! -----programmed by j.p. 04/04 in analogy to subroutine iso from above


! -----NDATA=4+7*NTS

  Integer (i4b), Parameter :: maxww = id_maxww
  Real (dp), Parameter :: clight = 2.997925D18

! input values

  Integer (i4b) :: ndata, nws
  Real (dp) :: data(ndata), tt, ee, weight, vturb, lambda0, dws(nws), sq(nws), &
    yhe

! local values
  Integer (i4b) :: i, ib, if, j, k, nts
  Real (dp) :: a, con, crit, d, d1, d2, db, dens, dlp, dls, ft, fq(maxww), p, &
    rba, rbhz, rft, rwt, srt, va, vb, vth, vmot, w, w1, w2, wf, wq, wt, y, np, &
    nhe, hei

  Real (dp) :: ts(6), ws(6), ds(6), ws1(6), ds1(6), ws2(6), ds2(6)

! .. external functions ..
  Real (dp) :: voigt
  External :: voigt

! -----frequential quantities

  Do i = 1, nws
    wq = lambda0 + dws(i)
    fq(i) = clight/wq
  End Do

  ft = clight/lambda0
  rft = 1.D0/ft

  srt = sqrt(tt)

  vth = sqrt(1.6504D8*tt/weight+vturb**2)
  vmot = clight*1.D-8/vth

  wt = clight*rft
  rwt = ft/clight

! move DATA to local

  i = 1
  dlp = data(i+1)
  dens = data(i+2)
  nts = data(i+3) + 0.5D0
  crit = data(i+4)
  k = i + 4

  Do j = 1, nts
    k = k + 1
    ts(j) = data(k)
  End Do

  Do j = 1, nts
    k = k + 1
    ws(j) = data(k)
  End Do

  Do j = 1, nts
    k = k + 1
    ds(j) = data(k)
  End Do

  Do j = 1, nts
    k = k + 1
    ws1(j) = data(k)
  End Do

  Do j = 1, nts
    k = k + 1
    ds1(j) = data(k)
  End Do
  Do j = 1, nts
    k = k + 1
    ws2(j) = data(k)
  End Do

  Do j = 1, nts
    k = k + 1
    ds2(j) = data(k)
  End Do


! RBHZ = 1/deltanu-dop
! RBA  = nu0/vth  (for transformation dlambda -> dnue

  rbhz = vmot*rft
  rba = vmot*rwt

  dlp = dlp*rba

! -----SET UP INTERPOLATION IN T
  Do j = 2, nts
    if = j
    If (ts(j)>=tt) Exit
  End Do

  ib = if - 1
  wf = (tt-ts(ib))/(ts(if)-ts(ib))

  If (tt<ts(ib)) wf = 0.0D0
  If (tt>ts(if)) wf = 1.0D0

! -----PERTURBER QUANTITIES
  y = ee/dens

! -----IMPACT WIDTH WIDTH (A)
  w = (wf*(ws(if)-ws(ib))+ws(ib))*y
  w = .5D0*w !                      tables give full halfwidths

! -----SHIFT (A)
  d = (wf*(ds(if)-ds(ib))+ds(ib))*y

! the same for protons and helium ions (note, that heiii approximated by
! heii!)
! test for crit

! rough estimate
  If (tt>=25000.) Then
    hei = 2.
  Else If (tt>=10000.) Then
    hei = 1.
  Else
    hei = 0.
  End If

  np = ee/(1.+hei*yhe) !            valid also for hei=0
  y = np/dens

  w1 = (wf*(ws1(if)-ws1(ib))+ws1(ib))
  If (w1==0.) Then
    d1 = 0.
  Else
    If (crit/w1<=np) Then
      w1 = .5D0*w1*y !              tables give full halfwidths
!     -----SHIFT (A)
      d1 = (wf*(ds1(if)-ds1(ib))+ds1(ib))*y
    Else !                          NOT RELIABLE
      w1 = 0.
      d1 = 0.
    End If
  End If


  If (hei==0.) Then
!   no contribution from heii/heiii
    w2 = 0.
  Else
    nhe = yhe*np !                  (heii/heiii)
    y = nhe/dens
    w2 = (wf*(ws2(if)-ws2(ib))+ws2(ib))
  End If

  If (w2==0.) Then
    d2 = 0.
  Else
    If (w2/=0 .And. crit/w2<=nhe) Then
      w2 = .5D0*w2*y !              tables give full halfwidths
!     -----SHIFT (A)
      d2 = (wf*(ds2(if)-ds2(ib))+ds2(ib))*y
    Else !                          NOT RELIABLE
      w2 = 0.
      d2 = 0.
    End If
  End If

! -----TOTAL WIDTH IN DOPPLER UNITS
  a = (w+w1+w2)*rba
  dls = (d+d1+d2)*rba

! con includes transformation to Doppler units via RBHZ,
! factor 1/sqrt(pi) from normalization of Voigt profile and
! factor 1/9 from using (8:1) weighted Voigt profiles to account
! for fine-structure of lower component, level 2P3.
! Works also for other cases if separation is set to zero

  con = 6.268773D-2*rbhz

! -----COMPUTE PROFILE

  Do j = 1, nws
    db = (ft-fq(j))*rbhz
    va = db - dls
    vb = va - dlp
    p = 8.0D0*voigt(a, va) + voigt(a, vb)
    sq(j) = con*p
  End Do

  Return
End Subroutine

!-----------------------------------------------------------------------

Subroutine indexx(n, arr, indx)

  Use :: nlte_type
  Use :: nlte_dim
  Implicit None

! .. parameters ..
  Integer (i4b), Parameter :: m = 7, nstack = 50
! ..
! .. scalar arguments ..
  Integer (i4b) :: n
! ..
! .. array arguments ..
  Real (dp) :: arr(n)
  Integer (i4b) :: indx(n)
! ..
! .. local scalars ..
  Real (dp) :: a
  Integer (i4b) :: i, indxt, ir, itemp, j, jstack, k, l
! ..
! .. local arrays ..
  Integer (i4b) :: istack(nstack)
! ..

  Do j = 1, n
    indx(j) = j
  End Do

  jstack = 0
  l = 1
  ir = n

100 Continue

  If (ir-l<m) Then

jloop: Do j = l + 1, ir
      indxt = indx(j)
      a = arr(indxt)
      Do i = j - 1, 1, -1
        If (arr(indx(i))<=a) Go To 110
        indx(i+1) = indx(i)
      End Do
      i = 0

110   Continue

      indx(i+1) = indxt
    End Do jloop

    If (jstack==0) Return
    ir = istack(jstack)
    l = istack(jstack-1)
    jstack = jstack - 2

  Else

    k = (l+ir)/2
    itemp = indx(k)
    indx(k) = indx(l+1)
    indx(l+1) = itemp

    If (arr(indx(l+1))>arr(indx(ir))) Then
      itemp = indx(l+1)
      indx(l+1) = indx(ir)
      indx(ir) = itemp
    End If

    If (arr(indx(l))>arr(indx(ir))) Then
      itemp = indx(l)
      indx(l) = indx(ir)
      indx(ir) = itemp
    End If

    If (arr(indx(l+1))>arr(indx(l))) Then
      itemp = indx(l+1)
      indx(l+1) = indx(l)
      indx(l) = itemp
    End If

    i = l + 1
    j = ir
    indxt = indx(l)
    a = arr(indxt)

120 Continue

    i = i + 1
    If (arr(indx(i))<a) Go To 120

130 Continue

    j = j - 1
    If (arr(indx(j))>a) Go To 130
    If (j<i) Go To 140
    itemp = indx(i)
    indx(i) = indx(j)
    indx(j) = itemp
    Go To 120

140 Continue

    indx(l) = indx(j)
    indx(j) = indxt
    jstack = jstack + 2
    If (jstack>nstack) Stop 'NSTACK TOO SMALL IN INDEXX'
    If (ir-i+1>=j-l) Then
      istack(jstack) = ir
      istack(jstack-1) = i
      ir = j - 1
    Else
      istack(jstack) = j - 1
      istack(jstack-1) = l
      l = i
    End If

  End If

  Go To 100

End Subroutine

!-----------------------------------------------------------------------

Subroutine transic_h_at(level, leveu, ncomp, lineno, nlevl, nlevu, natom, &
  xlamb, flu, nstark)
  Use :: nlte_type
  Use :: nlte_dim
  Use :: fund_const, Only: clight
  Use :: preformal_var, Only: fl, gl, zl, labl, nl, le, nlevh21
  Implicit None

! provides info on high lying H-levels

! .. scalar arguments ..
  Integer (i4b) :: ncomp
! ..
! .. array arguments ..
  Real (dp) :: flu(ncomp), xlamb(ncomp)
  Integer (i4b) :: lineno(ncomp), natom(ncomp), nlevl(ncomp), nlevu(ncomp), &
    nstark(ncomp)
  Character :: level(ncomp)*6, leveu(ncomp)*6

! local variables
  Real (dp), Parameter :: ryd = 109677.58D0

  Integer (i4b) :: nppal, nn1, nn2, i, leh20

! final consistency check
  If (ncomp/=1) Stop ' NCOMP NE 1 IN TRANSIC_H_app'
  If (lineno(ncomp)/=1) Stop ' LINENO NE 1 IN TRANSIC_H_app'
  If (nstark(ncomp)>1) Stop ' NSTARK > 1 IN TRANSIC_H_app'

  If (nl/=id_llevs) Stop ' NL NE ID_LLEVS IN TRANSIC_H_AT'

  nl = nl + 2 !                     two more levels

  labl(id_llevs+1) = level(1)
  labl(id_llevs+2) = leveu(1)

! levels of neutral hydrogen
  zl(id_llevs+1) = 0.
  zl(id_llevs+2) = 0.

! two extra levels
  nlevl(1) = id_llevs + 1
  nlevu(1) = id_llevs + 2

  nn1 = nppal(level(1))
  nn2 = nppal(leveu(1))

  Print *, 'approximate treatment for level ', level(1), '(', nn1, ') and ', &
    leveu(1), '(', nn2, ')'

  gl(id_llevs+1) = 2.*nn1**2
  gl(id_llevs+2) = 2.*nn2**2

  fl(id_llevs+1) = ryd/float(nn1**2)*clight
  fl(id_llevs+2) = ryd/float(nn2**2)*clight

  Do i = 1, id_llevs
    If (labl(i)=='H21') Go To 100
  End Do
  Stop ' H21 NOT FOUND IN LABL'

100 nlevh21 = i
  leh20 = le(nlevh21-1) !           since le(k-level) =0; thus with H120

  le(id_llevs+1) = leh20
  le(id_llevs+2) = leh20

  natom(1) = le(nn1)

! XLAMB and FLU from TRANSIC2

  Return
End Subroutine

!-----------------------------------------------------------------------

Subroutine vcs_balmer(xlamb, ll, lu, nws, dws, nts, ts, nes, es, nps, p)
! ***
! **********************************************************************
! ***
! ***	Read VCS tables, interpolate onto common wavelength grid, &
! ***   and transform to p(nu)d(nu)
! ***
! ***	From A&AS paper `Extended VCS Stark broadening tables for hydrogen
! ***	by Michael Lemke (michael@io.as.utexas.edu,
! ***                     ai26@a400.sternwarte.uni-erlangen.de)
! ***
! ***
! ***    adapted May 2013 by JP
! ***    read according to original data format Nov 2022 by JP
! ***
! ***	 2-AUG-1996 18:00:10.47		ML / Bamberg.
! ***
! **********************************************************************
! ***

  Use :: nlte_type
  Use :: nlte_dim
  Use :: fund_const, Only: clight
  Use :: preformal_var, Only: fpath
  
  Implicit None

  Logical, Parameter :: optout = .False. ! if true, original values will be
! tabulated

  Integer (i4b), Parameter :: maxw = id_maxww, maxt = id_maxtt, &
    maxne = id_maxne, maxps = id_maxps, nnws = 35

  Real (dp) :: xlamb

  Real (dp) :: dws(maxw), ts(maxt), es(maxne), p(maxps), ddws(0:nnws-1)

  Integer (i4b) :: nws1, nts, nes, nps, nws, ll, lu


  Integer (i4b), Parameter :: pmne = 17, pmt = 7, pmp = 65, mline = 20
! max # of n_e
! max # of T
! max # of profile points
! max # of H lines
! VCS profiles

  Real (dp), Parameter :: f0const = 2.5985D0*4.8032D-10

  Real (dp) :: svcs(pmt, pmne, 0:pmp), alpha(pmp), alpha_loc(pmp), &
    svcs_loc(pmt, pmne, 0:nnws-1)
! log Delta_alpha_min

  Real (dp) :: log_alpha0(mline), log_ne0(mline), log_t0(mline)
! log n_e_min
! log T_min
! Delta log Delta_alpha

  Real (dp) :: log_alpha_inc(mline), log_t_inc(mline), log_ne_inc(mline)
! Delta log n_e
! Delta log T
  Integer (i4b) :: mp(mline), mne(mline), mt(mline)
  Integer (i4b), Dimension (mline) :: nlow, nup

  Integer (i4b) :: i, j, k, line, nline, iline, it, kk, nnl, nnu, iu

  Real (dp) :: ne, f0, x, q, q1, renorm

  Character (1) :: null(pmt)

  Character (11) :: charline

!  Character (*), Parameter :: fpath = '../inicalc/DATA/' !in module preformal_var

  Data ddws/0.D0, 1.5849D-3, 2.5119D-3, 3.9811D-3, 6.3096D-3, 1.D-2, &
    1.5849D-2, 2.5119D-2, 3.9811D-2, 6.3096D-2, 1.D-1, 1.585D-1, 1.995D-1, &
    2.512D-1, 3.162D-1, 3.981D-1, 5.012D-1, 6.310D-1, 7.943D-1, 1.00D0, &
    1.259D0, 1.585D0, 1.995D0, 2.512D0, 3.981D0, 6.310D0, 1.0D1, 1.585D1, &
    2.512D1, 3.891D1, 6.310D1, 1.0D2, 1.585D2, 2.512D2, 3.981D2/



! READ IN VCS ARRAYS

  iu = 5 !                          since lu=1,2 already in use
  Open (iu, File=fpath//'balmer_lemke.dat', Status='OLD')

  Read (iu, *) nline
  If (nline>mline) Then
    Write (*, *) 'Table too big.  Not more than ', mline, ' lines.'
    Stop ' Table too big!'
  End If

  Do i = 1, nline
    Read (iu, *) nlow(i), nup(i), log_alpha0(i), log_ne0(i), log_t0(i), &
      log_alpha_inc(i), log_ne_inc(i), log_t_inc(i), mp(i), mne(i), mt(i)
    If (mp(i)>pmp .Or. mne(i)>pmne .Or. mt(i)>pmt) Then
      Write (*, *) 'Table too big in one of these:'
      Write (*, *) 'mp:', mp(i), ' >', pmp
      Write (*, *) 'mne:', mne(i), ' >', pmne
      Write (*, *) 'mt:', mt(i), ' >', pmt
      Stop ' Table too big in one of these!'
    End If
  End Do

  Do line = 1, nline
    Read (iu, Fmt='(A11)') charline
    Read (charline(4:5), '(i2)') nnl
    Read (charline(10:11), '(i2)') nnu

    If (nnl/=nlow(line) .Or. nnu/=nup(line)) Then
      Write (*, *) 'Inconsistency in table for', nlow(line), ' ->', nup(line)
      Stop ' Inconsistency in stark_balmer_lemke!'
    End If
    Read (iu, *)(((svcs(i,j,k),k=0,mp(line)),i=1,mt(line)), j=1, mne(line))

    If (nnl==ll .And. nnu==lu) Go To 100

  End Do
  Stop ' line not found stark_balmer_lemke'

100 Close (iu)
  iline = line

  nws1 = mp(iline)
  nts = mt(iline)
  nes = mne(iline)
  nws = nnws

  nps = nts*nes*nws !               and not nws1

  If (nts>maxt) Stop ' nts > maxt'
  If (nes>maxne) Stop ' nes > maxne'
  If (nps>maxps) Stop ' nps > maxps'

  Do k = 1, nws1
    alpha(k) = log_alpha0(iline) + log_alpha_inc(iline)*(k-1)
  End Do

  Do i = 1, nts
    it = nint(10.D0**(log_t0(iline)+log_t_inc(iline)*(i-1)))
    ts(i) = float(it)
  End Do

  Do j = 1, nes
    es(j) = 10.D0**(log_ne0(iline)+log_ne_inc(iline)*(j-1))
  End Do

  line = iline
  Write (*, *) 'nl = ', nlow(line), ';  nu = ', nup(line), &
    ' Stark broadening according to Lemke (1997)'

  If (optout) Then
    Do j = 1, nes
      ne = 10.D0**(log_ne0(line)+log_ne_inc(line)*(j-1))
      Write (*, *) ne, ' cm^-3'
!     svcs(i,j,0): quality flag (0 = OK)
      Write (*, '(a5,1x,a8,1x,7(:'' ('',i1,'')'',i6))') 'alpha', 'lambda', &
        (abs(int(svcs(i+1,j,0))), nint(10.D0**(log_t0(line)+ &
        log_t_inc(line)*i)), i=0, nts-1)

      f0 = f0const*ne**(2.D0/3.D0)

      Do k = 1, nws1
        x = log_alpha0(line) + log_alpha_inc(line)*(k-1)
        Do i = 1, nts
          null(i) = ' '
        End Do
        Write (*, '(f5.1,1x,1p,e8.3e1,1x,0p,7(f9.5,a))') x, 10.D0**x*f0, &
          (real(svcs(i,j,k)), null(i), i=1, nts)
      End Do
      Write (*, *)
    End Do

  End If

! now we interpolate the profile functions onto a common wavelength grid.
! interpolation is linear in log profile vs. log dalpha

! as a reference wavelength grid, we use the one from subr. GRIEM
! (consistent with the one from thom_new.dat)

! interpolation for all models

! rember: nws1 is number of freq. points in original grid
! : nws is number of freq. points for output grid, with indices(0:nws-1)


  Do j = 1, nes
!   recalculate alpha for given lambda
    f0 = f0const*es(j)**(2.D0/3.D0)

    Do k = 1, nws - 1
      alpha_loc(k) = log10(ddws(k)/f0)
    End Do

    Do k = 1, nws - 1
      If (alpha_loc(k)<alpha(1)) Then
        Do i = 1, nts
          svcs_loc(i, j, k) = svcs(i, j, 1) ! constant for (very) low dw
        End Do
      Else
        Do kk = 1, nws1 - 1
          If (alpha_loc(k)>=alpha(kk) .And. alpha_loc(k)<=alpha(kk+1)) Exit
        End Do
!       end condition
        If (kk==nws1) kk = nws1 - 1
        q = (alpha_loc(k)-alpha(kk))/(alpha(kk+1)-alpha(kk))
        q1 = 1.D0 - q
        Do i = 1, nts
          svcs_loc(i, j, k) = q1*svcs(i, j, kk) + q*svcs(i, j, kk+1)
        End Do
      End If
    End Do

!   value at dlam = 0, taken from first entry at orignal table
    Do i = 1, nts
      svcs_loc(i, j, 0) = svcs(i, j, 1)
    End Do

!   for tests
!   print*,j,log10(es(j))
!   k=0
!   print*,alpha_loc(1)-10.,svcs_loc(1,j,k),svcs_loc(7,j,k)
!   do k=1,nws-1
!   print*,alpha_loc(k),svcs_loc(1,j,k),svcs_loc(7,j,k)
!   enddo

  End Do

! finally, write dws (starting with dlam=0)
  Do k = 0, nws - 1
    dws(k+1) = ddws(k)
  End Do

! ... and profile p, including transformation
! remember that all wavelengths (and alpha) are in Angstrom, &
! and that the original normalization is 0.5 (and not unity) for 0<alpha<inf

  Do j = 1, nes
    f0 = f0const*es(j)**(2.D0/3.D0)
    Do k = 0, nws - 1
      renorm = (xlamb+ddws(k))**2/(f0*clight*1.D8)
      renorm = log10(renorm)
      Do i = 1, nts
        svcs_loc(i, j, k) = svcs_loc(i, j, k) + renorm
      End Do
    End Do
  End Do

! write transformed profile to output p (sequential order)
  kk = 0
  Do j = 1, nes
    Do i = 1, nts
      Do k = 0, nws - 1
        kk = kk + 1
        p(kk) = svcs_loc(i, j, k)
      End Do
    End Do
  End Do


  If (kk/=nps) Stop ' kk ne nps in vcs_balmer!'

  ts(1:nts) = log10(ts(1:nts))
  es(1:nes) = log10(es(1:nes))

  Return

End Subroutine

!-----------------------------------------------------------------------

Subroutine vcs_pb(xlamb, ll, lu, nws, dws, nts, ts, nes, es, nps, p)
! ***
! **********************************************************************
! ***
! ***	Read VCS tables, interpolate onto common wavelength grid, &
! ***   and transform to p(nu)d(nu)
! ***
! ***	From A&AS paper `Extended VCS Stark broadening tables for hydrogen
! ***	by Michael Lemke (michael@io.as.utexas.edu,
! ***                     ai26@a400.sternwarte.uni-erlangen.de)
! ***
! ***
! ***    adapted for Paschen & Brackett series, Nov 2022 by JP
! ***
! ***	 2-AUG-1996 18:00:10.47		ML / Bamberg.
! ***
! **********************************************************************
! ***

  Use :: nlte_type
  Use :: nlte_dim
  Use :: fund_const, Only: clight
  Use :: preformal_var, Only: fpath

  Implicit None

  Logical, Parameter :: optout = .False. ! if true, original values will be
! tabulated

  Integer (i4b), Parameter :: maxw = id_maxww, maxt = id_maxtt, &
    maxne = id_maxne, maxps = id_maxps, nnws = 35

  Real (dp) :: xlamb

! character*6 :: lablow, labup

  Real (dp) :: dws(maxw), ts(maxt), es(maxne), p(maxps), ddws(0:nnws-1)

  Integer (i4b) :: nws1, nts, nes, nps, nws, ll, lu, iu


  Integer (i4b), Parameter :: pmne = 17, pmt = 7, pmp = 65, mline = 20
! max # of n_e
! max # of T
! max # of profile points
! max # of H lines
! VCS profiles

  Real (dp), Parameter :: f0const = 2.5985D0*4.8032D-10

  Real (dp) :: svcs(pmt, pmne, 0:pmp), alpha(pmp), alpha_loc(pmp), &
    svcs_loc(pmt, pmne, 0:nnws-1)
! log Delta_alpha_min

  Real (dp) :: log_alpha0(mline), log_ne0(mline), log_t0(mline)
! log n_e_min
! log T_min
! Delta log Delta_alpha

  Real (dp) :: log_alpha_inc(mline), log_t_inc(mline), log_ne_inc(mline)
! Delta log n_e
! Delta log T
  Integer (i4b) :: mp(mline), mne(mline), mt(mline)
  Integer (i4b), Dimension (mline) :: nlow, nup

  Integer (i4b) :: i, j, k, line, nline, iline, it, kk, nnl, nnu

  Real (dp) :: ne, f0, x, q, q1, renorm

  Character (1) :: null(pmt)

  Character (11) :: charline

!  Character (*), Parameter :: fpath = '../inicalc/DATA/'

  Data ddws/0.D0, 1.5849D-3, 2.5119D-3, 3.9811D-3, 6.3096D-3, 1.D-2, &
    1.5849D-2, 2.5119D-2, 3.9811D-2, 6.3096D-2, 1.D-1, 1.585D-1, 1.995D-1, &
    2.512D-1, 3.162D-1, 3.981D-1, 5.012D-1, 6.310D-1, 7.943D-1, 1.00D0, &
    1.259D0, 1.585D0, 1.995D0, 2.512D0, 3.981D0, 6.310D0, 1.0D1, 1.585D1, &
    2.512D1, 3.891D1, 6.310D1, 1.0D2, 1.585D2, 2.512D2, 3.981D2/



! READ IN VCS ARRAYS

  iu = 5 !                          since lu=1,2 already in use
  Select Case (ll)

  Case (3)
    Open (iu, File=fpath//'paschen_lemke.dat', Status='OLD')
  Case (4)
    Open (iu, File=fpath//'brackett_lemke.dat', Status='OLD')
  Case Default
    Print *, ll
    Stop 'wrong lower level in vcs_pb'

  End Select

  If (lu>10) Stop 'upper level > 10 in vcs_pb'

  Read (iu, *) nline
  If (nline>mline) Then
    Write (*, *) 'Table too big.  Not more than ', mline, ' lines.'
    Stop ' Table too big!'
  End If

  Do i = 1, nline
    Read (iu, *) nlow(i), nup(i), log_alpha0(i), log_ne0(i), log_t0(i), &
      log_alpha_inc(i), log_ne_inc(i), log_t_inc(i), mp(i), mne(i), mt(i)
    If (mp(i)>pmp .Or. mne(i)>pmne .Or. mt(i)>pmt) Then
      Write (*, *) 'Table too big in one of these:'
      Write (*, *) 'mp:', mp(i), ' >', pmp
      Write (*, *) 'mne:', mne(i), ' >', pmne
      Write (*, *) 'mt:', mt(i), ' >', pmt
      Stop ' Table too big in one of these!'
    End If
  End Do

  Do line = 1, nline
    Read (iu, Fmt='(A11)') charline
    Read (charline(4:5), '(i2)') nnl
    Read (charline(10:11), '(i2)') nnu

    If (nnl/=nlow(line) .Or. nnu/=nup(line)) Then
      Write (*, *) 'Inconsistency in table for', nlow(line), ' ->', nup(line)
      Stop ' Inconsistency in table (vcs_pb)!'
    End If
    Read (iu, *)(((svcs(i,j,k),k=0,mp(line)),i=1,mt(line)), j=1, mne(line))

    If (nnl==ll .And. nnu==lu) Go To 100

  End Do
  Stop ' line not found lemke tables'

100 Close (iu)
  iline = line

  nws1 = mp(iline)
  nts = mt(iline)
  nes = mne(iline)
  nws = nnws

  nps = nts*nes*nws !               and not nws1

  If (nts>maxt) Stop ' nts > maxt'
  If (nes>maxne) Stop ' nes > maxne'
  If (nps>maxps) Stop ' nps > maxps'

  Do k = 1, nws1
    alpha(k) = log_alpha0(iline) + log_alpha_inc(iline)*(k-1)
  End Do

  Do i = 1, nts
    it = nint(10.D0**(log_t0(iline)+log_t_inc(iline)*(i-1)))
    ts(i) = float(it)
  End Do

  Do j = 1, nes
    es(j) = 10.D0**(log_ne0(iline)+log_ne_inc(iline)*(j-1))
  End Do

  line = iline
  Write (*, *) 'nl = ', nlow(line), ';  nu = ', nup(line), &
    ' Stark broadening according to Lemke (1997)'

  If (optout) Then
    Do j = 1, nes
      ne = 10.D0**(log_ne0(line)+log_ne_inc(line)*(j-1))
      Write (*, *) ne, ' cm^-3'
!     svcs(i,j,0): quality flag (0 = OK)
      Write (*, '(a5,1x,a8,1x,7(:'' ('',i1,'')'',i6))') 'alpha', 'lambda', &
        (abs(int(svcs(i+1,j,0))), nint(10.D0**(log_t0(line)+ &
        log_t_inc(line)*i)), i=0, nts-1)

      f0 = f0const*ne**(2.D0/3.D0)

      Do k = 1, nws1
        x = log_alpha0(line) + log_alpha_inc(line)*(k-1)
        Do i = 1, nts
          null(i) = ' '
        End Do
        Write (*, '(f5.1,1x,1p,e8.3e1,1x,0p,7(f9.5,a))') x, 10.D0**x*f0, &
          (real(svcs(i,j,k)), null(i), i=1, nts)
      End Do
      Write (*, *)
    End Do

  End If

! now we interpolate the profile functions onto a common wavelength grid.
! interpolation is linear in log profile vs. log dalpha

! as a reference wavelength grid, we use the one from subr. GRIEM
! (consistent with the one from thom_new.dat)

! interpolation for all models

! rember: nws1 is number of freq. points in original grid
! : nws is number of freq. points for output grid, with indices(0:nws-1)


  Do j = 1, nes
!   recalculate alpha for given lambda
    f0 = f0const*es(j)**(2.D0/3.D0)

    Do k = 1, nws - 1
      alpha_loc(k) = log10(ddws(k)/f0)
    End Do

    Do k = 1, nws - 1
      If (alpha_loc(k)<alpha(1)) Then
        Do i = 1, nts
          svcs_loc(i, j, k) = svcs(i, j, 1) ! constant for (very) low dw
        End Do
      Else
        Do kk = 1, nws1 - 1
          If (alpha_loc(k)>=alpha(kk) .And. alpha_loc(k)<=alpha(kk+1)) Exit
        End Do
!       end condition
        If (kk==nws1) kk = nws1 - 1
        q = (alpha_loc(k)-alpha(kk))/(alpha(kk+1)-alpha(kk))
        q1 = 1.D0 - q
        Do i = 1, nts
          svcs_loc(i, j, k) = q1*svcs(i, j, kk) + q*svcs(i, j, kk+1)
        End Do
      End If
    End Do

!   value at dlam = 0, taken from first entry at orignal table
    Do i = 1, nts
      svcs_loc(i, j, 0) = svcs(i, j, 1)
    End Do

!   for tests
!   print*,j,log10(es(j))
!   k=0
!   print*,alpha_loc(1)-10.,svcs_loc(1,j,k),svcs_loc(7,j,k)
!   do k=1,nws-1
!   print*,alpha_loc(k),svcs_loc(1,j,k),svcs_loc(7,j,k)
!   enddo

  End Do

! finally, write dws (starting with dlam=0)
  Do k = 0, nws - 1
    dws(k+1) = ddws(k)
  End Do

! ... and profile p, including transformation
! remember that all wavelengths (and alpha) are in Angstrom, &
! and that the original normalization is 0.5 (and not unity) for 0<alpha<inf

  Do j = 1, nes
    f0 = f0const*es(j)**(2.D0/3.D0)
    Do k = 0, nws - 1
      renorm = (xlamb+ddws(k))**2/(f0*clight*1.D8)
      renorm = log10(renorm)
      Do i = 1, nts
        svcs_loc(i, j, k) = svcs_loc(i, j, k) + renorm
      End Do
    End Do
  End Do

! write transformed profile to output p (sequential order)
  kk = 0
  Do j = 1, nes
    Do i = 1, nts
      Do k = 0, nws - 1
        kk = kk + 1
        p(kk) = svcs_loc(i, j, k)
      End Do
    End Do
  End Do


  If (kk/=nps) Stop ' kk ne nps in vcs_pb!'

  ts(1:nts) = log10(ts(1:nts))
  es(1:nes) = log10(es(1:nes))

  Return

End Subroutine

!-----------------------------------------------------------------------

Function igenio(kk, ii)

  Use :: nlte_type
  Use :: nlte_dim
  Use :: princesa_var, Only: iong, nions

  Implicit None

! ------ this function returns the value of 'general ion' for atom kk and
! ion ii

! .. scalar arguments ..
  Integer (i4b) :: ii, kk, igenio
! ..
! .. local scalars ..
  Integer (i4b) :: i
! ..

  igenio = 0

  Do i = 1, kk - 1
    igenio = igenio + nions(i)
  End Do

  igenio = igenio + ii
  If (igenio>iong) Stop ' error in general ion - igenio '

  Return

End Function
