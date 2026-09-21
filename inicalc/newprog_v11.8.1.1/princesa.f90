!version 1.2 (see below)
Module princesa_var

  Use :: nlte_type
  Use :: nlte_dim
  Implicit None
! ----------------------------------------------------------------------
! atdnin

  Integer (i4b) :: nl = 0, nlx = 0, ns = 0, nat = 0
  Integer (i4b) :: ion, iong = 0
  Integer (i4b), Dimension (id_atoms) :: la0 = 0, la1 = 0
  Integer (i4b), Dimension (id_atoms) :: ixa0 = 0, ixa1 = 0
  Integer (i4b), Dimension (id_atoms) :: isa0 = 0, isa1 = 0
  Integer (i4b), Dimension (id_atoms) :: nions = 0
  Integer (i4b), Dimension (id_iones) :: ifirsl = 1, ifirx = 1, ifirs = 1
  Integer (i4b), Dimension (id_llevs) :: kl = 0, le = 0, li = 0
  Integer (i4b), Dimension (id_xlevs) :: klx = 0, lix = 0
  Integer (i4b), Dimension (id_slevs) :: lis = 0, ns0 = 0, ns1 = 0
! ----------------------------------------------------------------------
! atdnre

  Real (dp), Dimension (id_atoms) :: zeff = 0., weight = 0., abund = 0.
  Real (dp), Dimension (id_llevs) :: gl = 0., fl = 0., zl = 0.
  Real (dp), Dimension (id_xlevs) :: glx = 0., flx = 0.
  Real (dp), Dimension (id_slevs) :: gs0 = 0., gs1 = 0., ryd = 0., qd = 0.
! ----------------------------------------------------------------------
! atomdatc

  Character (6), Dimension (id_atoms) :: labat
  Character (6), Dimension (id_llevs) :: labl, lkl
  Character (6), Dimension (id_xlevs) :: labx, lkx
  Character (6), Dimension (id_slevs) :: labls, kls
! ----------------------------------------------------------------------
! dat

  Integer (i4b) :: ndat = 0
  Real (dp), Dimension (id_ndata) :: data = 0.
! ----------------------------------------------------------------------
! frec

  Integer (i4b) :: nf = 0 !         only needed for perfil
! ----------------------------------------------------------------------
! inout

  Integer (i4b) :: iu, iou
! ----------------------------------------------------------------------
! tranin

  Integer (i4b) :: nlong, indtt = 0, indtr = 0
  Integer (i4b) :: index1 = 0, index2 = 0, index3 = 0, index4 = 0, index5 = 0
  Integer (i4b), Dimension (id_nttrd) :: iaux2
  Integer (i4b), Dimension (id_rcbbt) :: iform1, numd1, indat1, labl1, labu1
  Integer (i4b), Dimension (id_rcbxt) :: iform2, numd2, indat2
  Integer (i4b), Dimension (id_cbsft) :: iform3, numd3, indat3, labl3
  Integer (i4b), Dimension (id_rbftr) :: iform4, numd4, indat4, labl4, &
    drexplicit
  Integer (i4b), Dimension (id_rfftr) :: iform5, numd5, indat5, ipare5
  Integer (i4b), Dimension (id_ntran, 2) :: iaux1
! ----------------------------------------------------------------------
! tranlo

  Logical :: qcl, qdb, qrt
! ----------------------------------------------------------------------
! tranre

  Real (dp) :: ancdo, tl, wlc
  Real (dp), Dimension (id_nlong) :: dl
  Real (dp), Dimension (id_rbftr) :: frecin, frecfn
  Real (dp), Dimension (id_cbsft) :: flcbf
! ----------------------------------------------------------------------
! transc

  Character (6), Dimension (id_rcbxt) :: labl2, labu2, paren2
  Character (6), Dimension (id_cbsft) :: paren3
  Character (6), Dimension (id_rbftr) :: paren4

End Module

!ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
!ccc                                                                cccc
!ccc         read atomic data  subroutines (e. santolaya, 1991)     cccc
!ccc                                                                cccc
!ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc

!15/01/99 change to f90 by j. puls
!03/03/99 more f90 changes
!11/03/99 obsolete fortran features removed (j.puls)


!version 1.1 (march 2006) : initial freq. grid from DETAIL input
!no longer read!!! (can be present or absent in file)
!version 1.2 (july 2009) : new array DREXPLICIT(INDEX4): set to 1 if
!explicit dr data are present (otherwise set to 0)
!version 1.3 (dec. 2016) : new array FLCBF(INDEX3): frequencies
!for CBF transitions

Subroutine prince

  Use :: nlte_type
  Use :: nlte_dim
  Use :: princesa_var, Only: iou, iu, labl
  Use :: ffr_error
  Implicit None

! main program for reading and storaging atomic data
! variables
! ..
! .. local scalars ..

  Integer (i4b) :: i, j

  Character :: dc*6, fichero*32, ret*4
! ..
! .. external subroutines ..
  External :: abu, ffracc, ffrcdc, ffrout, rdatom, transic
! ..

  iu = 37
  iou = 38
  dc = ':T'
  Open (1, File='ATOM_FILE', Status='OLD')
  Rewind 1
  Read (1, Fmt='(A)') fichero

  Close (1)
  Open (Unit=iu, File=fichero, Form='FORMATTED', Status='OLD')
  Open (Unit=iou, File='control.dat', Status='UNKNOWN')

  Call ffrout(iou, ret)
  If (on_error(ret)) Go To 100

  Call ffracc(iu, ret)
  If (on_error(ret)) Go To 110

  Call ffrcdc(dc)

rloop: Do

    Call rdatom(ret)
    If (ret/='RET0' .And. ret/='RET3') Go To 120
    If (ret=='RET3') Go To 130

    Call transic(ret)
    If (on_error(ret)) Go To 120

  End Do rloop

! error in output unit exit
100 Continue
  Write (*, Fmt='(A)') ' ERROR IN OUTPUT UNIT NUMBER '
  Stop

! error in input unit exit
110 Continue
  Write (*, Fmt='(A)') ' ERROR IN INPUT UNIT NUMBER '
  Stop

! error conditions
120 Select Case (ret)
  Case ('RET1') !                   end of file "no esperado" exit
    Write (*, Fmt='(A)') ' EOF NO ESPERADO. ALGO FALLA! '
    Stop
  Case ('RET2') !                   error in ffr subroutines exit
    Write (*, Fmt='(A)') ' ERROR IN FFR SUBROUTINES '
    Stop
  Case Default
    Stop ' WRONG ERROR CONDITION IN PRINCESA'
  End Select

! normal exit: end of file
130 Continue
  Call abu
  Write (*, Fmt='(A)') ' FILE IS OVER. SUCCESSFUL!! '
  Close (iu)

  Close (iou)

! JO FEB 2023: CHECK THAT ALL LEVELS ARE DIFFERENT
  Do i = 1, id_llevs - 1
    Do j = i + 1, id_llevs
      If (labl(i)==labl(j)) Then
        Print *, i, j, labl(i), labl(j)
        Stop ' DIFFERENT LEVELS HAVE IDENTICAL LABELS'
      End If
    End Do
  End Do

  Return


End Subroutine


!ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc

Subroutine rdatom(retchar)

  Use :: nlte_type
  Use :: nlte_dim
  Use :: princesa_var, Only: nl, nlx, ns, nat, ion, iong, la0, la1, ixa0, &
    ixa1, isa0, isa1, nions, ifirsl, ifirx, ifirs, kl, le, li, klx, lix, ns0, &
    ns1, lis, zeff, weight, abund, gl, fl, zl, glx, flx, gs0, gs1, ryd, qd, &
    labat, labl, lkl, labx, lkx, labls, kls
  Use :: ffr_error

  Implicit None

! it reads and storages atomic data
! variables
! nat: atom number (for internal counting)
! ion: number of ion inside each atomo
! iong: number of ion in the general counting
! labat( ): atom label (name)
! zeff( ): charge of the nucleus (=atomic number)
! weight( ): atomic weight
! abund( ): abundance
! la0( ): index of first l level of the atom
! ixa0( ): same with levels x
! isa0( ): same with levels s
! zef: effective charge of the atom
! ifirsl( ): first l level of a given ion
! nl: index for levels l
! kl( ): index of the first level of the next ion
! labl( ): level l label
! gl( ): stadistic weight of level l
! fl( ): ionization frequency of level l
! lkl( ): label of parental level of a given l level
! le( ): number of atom of a level l
! li( ): number of ion of a level l
! zl( ): effective charge of level l
! ifirx( ): first x level of a given ion
! nlx: index for levels x
! labx( ): labels of x levels
! glx( ): stadistic weight of levels x
! flx( ): ionization frequency for levels x
! lkx( ): label of parent level to a x level
! klx( ): index of the first l level of the next ion
! (for x levels)
! lix( ): number of ion of a level x
! ifirs( ): first s level of a given ion
! ns: index for levels s
! labls( ): label of levels s
! gs0( ): stadistic weight of first level of the serial
! gs1( ): same with the last
! ryd( ): rydberg constant of level s
! qd( ): quantum deffect of level s
! ns0( ): index of the first level of the serial
! ns1( ): same with the last
! kls( ): label of the parent level for a s level
! lis( ): ion of level s
! la1( ): number of l levels in the atom (can be avoid)
! ixa1( ): same with levels x (can be avoid)
! isa1( ): same with levels s (can be avoid)
! nions( ): total number of ions in a given atom
! ..
! .. local scalars ..
  Real (dp) :: realvar, zef
  Integer (i4b) :: i, ij, integ, nc
  Character :: label*4, key*6, retchar*4, ret*4
! ..
! .. external subroutines ..
  External :: ffrkey, ffrloc, ffrnum
! ..
! .. common blocks ..
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
    Stop ' WRONG RETURN IN RDATOM'
  End Select

  nat = nat + 1

  ion = 0
  Call ffrkey(key, nc, ret)
  If (on_error(ret)) Go To 150
  labat(nat) = key

  Call ffrnum(realvar, integ, ret)
  If (on_error(ret)) Go To 150
  zeff(nat) = realvar

  Call ffrnum(realvar, integ, ret)
  If (on_error(ret)) Go To 150
  weight(nat) = realvar

  Call ffrnum(realvar, integ, ret)
  If (on_error(ret)) Go To 150
  If (realvar<0.) realvar = 10**realvar

  abund(nat) = realvar
  la0(nat) = nl + 1
  ixa0(nat) = nlx + 1
  isa0(nat) = ns + 1

  zef = zeff(nat) - 1.D0
100 Continue

  Call ffrkey(key, nc, ret)
  If (on_error(ret)) Go To 150

! ------------- l  levels
  If (key=='L') Then
    zef = zef + 1.D0
    ion = ion + 1
    iong = iong + 1

    ifirsl(iong) = nl + 1
110 Continue
    Call ffrkey(key, nc, ret)
    If (on_error(ret)) Go To 150

    If (key=='0') Then
      Do i = ifirsl(iong), nl
        kl(i) = nl + 1
      End Do
      Do ij = iong + 1, id_iones
        ifirsl(ij) = nl + 1
      End Do

      Go To 100
    Else
      nl = nl + 1
      labl(nl) = key
      Call ffrnum(realvar, integ, ret)
      If (on_error(ret)) Go To 150
      gl(nl) = realvar
      Call ffrnum(realvar, integ, ret)
      If (on_error(ret)) Go To 150
      fl(nl) = realvar
      Call ffrkey(key, nc, ret)
      If (on_error(ret)) Go To 150
      lkl(nl) = key
      le(nl) = nat
      li(nl) = ion
      zl(nl) = zef

      Go To 110


    End If

!   ----------------- x  levels
  Else If (key=='X') Then

    ifirx(iong) = nlx + 1
120 Continue
    Call ffrkey(key, nc, ret)
    If (on_error(ret)) Go To 150

    If (key=='0') Then
      Do ij = iong + 1, id_iones
        ifirx(ij) = nlx + 1
      End Do
      Go To 100

    End If
    nlx = nlx + 1
    labx(nlx) = key
    Call ffrnum(realvar, integ, ret)
    If (on_error(ret)) Go To 150
    glx(nlx) = realvar

    Call ffrnum(realvar, integ, ret)
    If (on_error(ret)) Go To 150
    flx(nlx) = realvar

    Call ffrkey(key, nc, ret)
    If (on_error(ret)) Go To 150
    lkx(nlx) = key
    klx(nlx) = nl + 1
    lix(nlx) = ion


    Go To 120

!   -------------- s  levels
  Else If (key=='S') Then
    ifirs(iong) = ns + 1
130 Continue
    Call ffrkey(key, nc, ret)
    If (on_error(ret)) Go To 150

    If (key/='0') Then
      ns = ns + 1
      labls(ns) = key
      Call ffrnum(realvar, integ, ret)
      If (on_error(ret)) Go To 150
      gs0(ns) = realvar

      Call ffrnum(realvar, integ, ret)
      If (on_error(ret)) Go To 150
      gs1(ns) = realvar

      Call ffrnum(realvar, integ, ret)
      If (on_error(ret)) Go To 150
      ryd(ns) = realvar

      Call ffrnum(realvar, integ, ret)
      If (on_error(ret)) Go To 150
      qd(ns) = realvar

      Call ffrnum(realvar, integ, ret)
      If (on_error(ret)) Go To 150
      ns0(ns) = integ

      Call ffrnum(realvar, integ, ret)
      If (on_error(ret)) Go To 150
      ns1(ns) = integ

      Call ffrkey(key, nc, ret)
      If (on_error(ret)) Go To 150
      kls(ns) = key
      lis(ns) = ion

      Go To 130

    End If
    Do ij = iong + 1, id_iones
      ifirs(ij) = ns + 1
    End Do

    Go To 100

!   --------------------  k  levels
  Else If (key=='K') Then
    nl = nl + 1
    zef = zef + 1
    ion = ion + 1
    iong = iong + 1
    ifirsl(iong) = nl
    Call ffrkey(key, nc, ret)
    If (on_error(ret)) Go To 150
    labl(nl) = key

    Call ffrnum(realvar, integ, ret)
    If (on_error(ret)) Go To 150
    gl(nl) = realvar

    Call ffrnum(realvar, integ, ret)
    If (on_error(ret)) Go To 150
    fl(nl) = realvar

    Call ffrkey(key, nc, ret)
    If (on_error(ret)) Go To 150
    lkl(nl) = key
    zl(nl) = zef
    le(nl) = nat
    li(nl) = ion
    kl(nl) = nl

  End If

! the atom is over.
  la1(nat) = nl - la0(nat) + 1
  ixa1(nat) = nlx - ixa0(nat) + 1
  isa1(nat) = ns - isa0(nat) + 1


  nions(nat) = ion

  retchar = 'RET0'
  Return
! no more atoms to read. the main subroutine (prince) ends
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
    Stop ' WRONG ERROR CONDITION IN RDATOM'
  End Select

End Subroutine

!ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc

Subroutine transic(retchar)

  Use :: nlte_type
  Use :: nlte_dim
  Use :: fund_const, Only: clight
  Use :: princesa_var, Only: nl, ns, nat, ifirsl, le, li, weight, fl, labl, &
    labls, data, ndat, iu, nlong, indtt, index1, index2, index3, index4, &
    index5, iform1, iform2, iform3, iform4, iform5, numd1, numd2, numd3, &
    numd4, numd5, indat1, indat2, indat3, indat4, indat5, iaux1, labl4, labl1, &
    labu1, labl3, ipare5, qdb, qcl, qrt, tl, wlc, dl, frecin, frecfn, labl2, &
    paren2, paren3, paren4, labu2, drexplicit, flcbf, ancdo
  Use :: ffr_error

  Implicit None

! it reads and storages all data concerning transitions
! variables:
! tabla: posible transition types
! iform: number of formula
! ndatos: number of data required for the transition
! qdb: flag for the db case (doppler widths sampling)
! nlong: number of points in the profile sampling
! dl( ): points for the sampling in lambda
! tl: line temperatura
! ancdo: doppler width
! qcl: flag for the cl case (central wavelenght)
! qrt: flag for rbb or rbx transition
! indtt: index for the total number of transitions
! iaux1( , ): matrix for storaging the transition type
! and the pointer for other transitions
! counts
! index1( ): index for rbb or cbb transistions
! iform1( ): formula number for rbb or cbb
! numd1( ): same with the number of data
! labl1( ): lower level
! labu1( ): upper level
! indat1( ): index of the first data of this transition,
! storaged in data
! ndat1( ): number of data in data for this transition
! ndat: pointer for data
! data( ): data
! wlc: central wavelength
! for the other transition type, joint in five groups,
! (rbb-cbb; rbx-cbx; cbs-cbf; rbf; rff)
! similar vector and matrices are defined, but using
! the numbers 2, 3, 4 and 5 instead of 1. there are
! only one different:
! paren_( ): parent level label
! frecin( ): actual edge energy  (rbf) [in Kayser]
! frecfn( ): edge energy with respect to next ground
! state [in Kayser]
! flcbf( ): edge frequency with respect to next ground
! state for cbf ionization

! ..
! .. local scalars ..
  Real (dp) :: deltae, frecin2, realvar
  Integer (i4b) :: i, iform, ij, integ, iqi, isi, l0, m, nc, ndatos, nn, nnu
  Character :: key*6, pare*6, retchar*4, ret*4
! ..
! .. local arrays ..
  Character :: tabla(8)*6
! ..
! .. external functions ..
  Real (dp) :: doppler
  Integer (i4b) :: igenio
  External :: doppler, igenio
! ..
! .. external subroutines ..
  External :: ffrkey, ffrnum, hallcl, hallcx, ident, perfil
! ..
! .. intrinsic functions ..
! ..
! .. data statements ..
  Data tabla/'RBB', 'CBB', 'RBX', 'CBX', 'CBS', 'CBF', 'RBF', 'RFF'/
! ..

! JP changed May 3rd 2001
  qcl = .False.

100 Continue

  Call ffrkey(key, nc, ret)
  If (on_error(ret)) Go To 280

  If (key=='TY') Then
    Call ffrkey(key, nc, ret)
    If (on_error(ret)) Go To 280
    Call ident(key, tabla, 8, m)

    Call ffrnum(realvar, iform, ret)
    If (on_error(ret)) Go To 280

    Call ffrnum(realvar, ndatos, ret)
    If (on_error(ret)) Go To 280

!   son los numeros de formula y de datos, respectivamente

    Go To 100

!   sampling in doppler widths. the sampling is storage in dl.
  Else If (key=='DB') Then
    qdb = .True.
    Call ffrnum(realvar, nlong, ret)
    If (on_error(ret)) Go To 280

    Do i = 1, nlong
      Call ffrnum(realvar, integ, ret)
      If (on_error(ret)) Go To 280
      dl(i) = realvar
    End Do


    Go To 100

!   sampling in angstroms
  Else If (key=='DL') Then
    qdb = .False.
    Call ffrnum(realvar, nlong, ret)
    If (on_error(ret)) Go To 280

    Do i = 1, nlong
      Call ffrnum(realvar, integ, ret)
      If (on_error(ret)) Go To 280
      dl(i) = realvar
    End Do


    Go To 100

!   line temperature
  Else If (key=='TL') Then
    Call ffrnum(tl, integ, ret)
    If (on_error(ret)) Go To 280
    ancdo = doppler(tl, weight(nat))

    Go To 100

!   central wavelength
  Else If (key=='CL') Then
    qcl = .True.

    Go To 100

!   transition type
  Else
    Call ident(key, tabla, 8, m)
    If (m==0) Then
      Go To 290

    Else
      qrt = .False.
      If ((m==1) .Or. (m==3)) qrt = .True.
      Go To (110, 110, 150, 150, 170, 170, 210, 250) m
!     goto calculado, para que haga una cosa u otra segun el tipo de trans.


    End If


  End If

! -------------rbb or cbb transitions ------------
110 Continue

  Call ffrkey(key, nc, ret)
  If (on_error(ret)) Go To 280

120 Continue
  Do iqi = 1, nl
    If (labl(iqi)==key) Then
      nn = iqi
      Go To 130
    End If
  End Do

  Stop 'ERROR IN TRANSIC - LABEL L NOT FOUND, BB'

130 Continue
  indtt = indtt + 1
  iaux1(indtt, 1) = m
  index1 = index1 + 1
  iform1(index1) = iform
  numd1(index1) = ndatos + 2
  labl1(index1) = nn

  Call ffrkey(key, nc, ret)
  If (on_error(ret)) Go To 280

  Do iqi = 1, nl
    If (labl(iqi)==key) Then
      nnu = iqi
      Go To 140
    End If
  End Do

  Stop 'ERROR IN TRANSIC - LABEL U NOT FOUND, BB'

140 Continue
  labu1(index1) = nnu
  indat1(index1) = ndat + 1

  iaux1(indtt, 2) = index1

! radiative transition
  If (qrt) Then
    If (qcl) Then
      Call ffrnum(wlc, integ, ret)
      If (on_error(ret)) Go To 280
    Else
      Call hallcl(index1)
    End If

    Call perfil
  Else
    Call hallcl(index1)

  End If

! read data
  ndat = ndat + 1
  data(ndat) = dble(m)
  ndat = ndat + 1
  data(ndat) = wlc

  Do i = 1, ndatos
    ndat = ndat + 1
    Call ffrnum(realvar, integ, ret)
    If (on_error(ret)) Go To 280
    data(ndat) = realvar
!   print*,nn,nnu,realvar
  End Do

! read a new keyword
  Call ffrkey(key, nc, ret)
  If (on_error(ret)) Go To 280

  If (key=='0') Then
    qcl = .False.
    Go To 100
  End If

  Go To 120

! --------------- rbx or cbx transitions -----------
150 Continue
  Call ffrkey(pare, nc, ret)
  If (on_error(ret)) Go To 280

  Call ffrkey(key, nc, ret)
  If (on_error(ret)) Go To 280

160 Continue
  indtt = indtt + 1
  iaux1(indtt, 1) = m
  index2 = index2 + 1
  iform2(index2) = iform
  numd2(index2) = ndatos
  indat2(index2) = ndat + 1
  paren2(index2) = pare
  iaux1(indtt, 2) = index2
  labl2(index2) = key

  Call ffrkey(key, nc, ret)
  If (on_error(ret)) Go To 280

  labu2(index2) = key

! radiative transition
  If (qrt) Then
    If (qcl) Then
      Call ffrnum(wlc, integ, ret)
      If (on_error(ret)) Go To 280
    Else
      Call hallcx(index2)
    End If

    Call perfil

  End If

! read data
  Do i = 1, ndatos
    ndat = ndat + 1
    Call ffrnum(realvar, integ, ret)
    If (on_error(ret)) Go To 280
    data(ndat) = realvar
  End Do

! read new keyword
  Call ffrkey(key, nc, ret)
  If (on_error(ret)) Go To 280

  If (key=='0') Then
    qcl = .False.
    Go To 100
  End If

  Go To 160

! ------------- cbs or cbf transitions  ------------------
170 Continue
  Call ffrkey(pare, nc, ret)
  If (on_error(ret)) Go To 280

  Call ffrkey(key, nc, ret)
  If (on_error(ret)) Go To 280

180 Continue
  indtt = indtt + 1
  iaux1(indtt, 1) = m
  index3 = index3 + 1
  iform3(index3) = iform
  numd3(index3) = ndatos
  indat3(index3) = ndat + 1
  iaux1(indtt, 2) = index3
  paren3(index3) = pare

  Do isi = 1, nl
    If (labl(isi)==key) Then
      nn = isi
      Go To 190
    End If
  End Do

  Do isi = 1, ns
    If (labls(isi)==key) Then
      nn = -isi
      Go To 190
    End If
  End Do

  Stop 'ERROR IN TRANSIC - LABEL L OR S NOT FOUND -CBF-CBS'

190 Continue

  labl3(index3) = nn
  flcbf(index3) = fl(nn)

! calculation of actual edge frequency (if not ground state
! ionization)

  Do isi = 1, nl
    If (labl(isi)==pare) Then
      nn = isi
      Go To 200
    End If
  End Do

  Stop 'ERROR IN TRANSIC - LABEL PARE (CBF) NOT FOUND'

200 Continue
  ij = igenio(le(nn), li(nn))
  l0 = ifirsl(ij)
  If (l0==nn) Then
!   groundstate ionization
    deltae = 0.D0
  Else
!   excited state ionization
    deltae = fl(l0) - fl(nn)
  End If

  flcbf(index3) = flcbf(index3) + deltae

! read data
  Do i = 1, ndatos
    ndat = ndat + 1
    Call ffrnum(realvar, integ, ret)
    If (on_error(ret)) Go To 280
    data(ndat) = realvar
  End Do

! read new keyword
  Call ffrkey(key, nc, ret)
  If (on_error(ret)) Go To 280

  If (key=='0') Go To 100

  Go To 180

! ------------------ rbf transition  -------------------
210 Continue
  Call ffrkey(pare, nc, ret)
  If (on_error(ret)) Go To 280

  Call ffrkey(key, nc, ret)
  If (on_error(ret)) Go To 280

220 Continue
  indtt = indtt + 1
  iaux1(indtt, 1) = m
  index4 = index4 + 1
  iform4(index4) = iform
  numd4(index4) = ndatos
  indat4(index4) = ndat + 1
  iaux1(indtt, 2) = index4
  paren4(index4) = pare
! NEW; CHECK WHETHER DR WITH EXPLICIT STABILIZING RATES
  If (iform==21) Then
    drexplicit(index4) = 1
  Else
    drexplicit(index4) = 0
  End If

  Do isi = 1, nl
    If (labl(isi)==key) Then
      nn = isi
      Go To 230
    End If
  End Do

  Stop 'ERROR IN TRANSIC - LABEL L NOT FOUND'

230 Continue
  labl4(index4) = nn
  frecfn(index4) = fl(nn)

! calculation of actual edge frequency (if not ground state
! ionization)

  Do isi = 1, nl
    If (labl(isi)==pare) Then
      nn = isi
      Go To 240
    End If
  End Do

  Stop 'ERROR IN TRANSIC - LABEL PARE NOT FOUND'

240 Continue
  ij = igenio(le(nn), li(nn))
  l0 = ifirsl(ij)
  If (l0==nn) Then
!   groundstate ionization
    deltae = 0.D0
  Else
!   excited state ionization
    deltae = fl(l0) - fl(nn)
  End If

  Call ffrkey(key, nc, ret)
  If (on_error(ret)) Go To 280

  Call ffrkey(key, nc, ret)
  If (on_error(ret)) Go To 280

  frecin2 = frecfn(index4) + deltae

  frecin(index4) = frecin2/clight
  frecfn(index4) = frecfn(index4)/clight

! a la salida de aqui frecin y frecfn estan en cm-1
! read data
  Do i = 1, ndatos
    ndat = ndat + 1
    Call ffrnum(realvar, integ, ret)
    If (on_error(ret)) Go To 280
    data(ndat) = realvar
  End Do

! read new keyword
  Call ffrkey(key, nc, ret)
  If (on_error(ret)) Go To 280

  If (key=='0') Go To 100

  Go To 220

! ----------------- rff transitions ------------------------
250 Continue

  Call ffrkey(key, nc, ret)
  If (on_error(ret)) Go To 280

260 Continue
  indtt = indtt + 1
  iaux1(indtt, 1) = m
  index5 = index5 + 1
  iform5(index5) = iform
  numd5(index5) = ndatos
  indat5(index5) = ndat + 1
  iaux1(indtt, 2) = index5

  Do isi = 1, nl
    If (labl(isi)==key) Then
      nn = isi
      Go To 270
    End If
  End Do

  Stop 'ERROR IN TRANSIC - LABEL L NOT FOUND'

270 Continue
  ipare5(index5) = nn
! read data
  Do i = 1, ndatos
    ndat = ndat + 1
    Call ffrnum(realvar, integ, ret)
    If (on_error(ret)) Go To 280
    data(ndat) = realvar
  End Do

! read new keyword
  Call ffrkey(key, nc, ret)
  If (on_error(ret)) Go To 280

  If (key=='0') Go To 100

  Go To 260

! --------------------------------------------------------
! error handling
280 Select Case (ret)
  Case ('RET1') !                   end of file exit
    Write (*, Fmt='(A)') ' END OF FILE '
    retchar = 'RET1'
    Return
  Case ('RET2') !                   error exit
    Write (*, Fmt='(A)') 'ERROR IN FFR SUBROUTINES '
    retchar = 'RET2'
    Return
  Case Default
    Stop ' WRONG ERROR CONDITION IN TRANSIC'
  End Select

290 Continue
  Backspace iu

  retchar = 'RET0'
  Return

End Subroutine

!-----------------------------------------------------------------------

Subroutine perfil

  Use :: nlte_type
  Use :: nlte_dim
  Use :: princesa_var, Only: nf, nlong, indtt, indtr, iaux2, qdb, dl

  Implicit None

! calculates the sampling of the atomic profile for radiative
! bound-bound transitions
! ..
! .. local scalars ..
  Integer (i4b) :: i
! ..
! .. local arrays ..
  Integer (i4b) :: ifptr(id_nttrd), nfptr(id_nttrd)
! ..
! ..
! this routine should never be called, since profiles separately treated, THUS

  indtr = indtr + 1
  iaux2(indtr) = indtt
  nfptr(indtr) = nlong

  ifptr(indtr) = nf + 1
! sampling in doppler width units
  If (qdb) Then
    Do i = 1, nlong
!     DL(I) = DL(I)* (ANCDO*WLC)
      dl(i) = 0.
    End Do
  End If
! now dl is in angstroms
! construction of the sampling
! do 20 i=1,nlong
! nf=nf+1
! fq(nf)=clight/(wlc+dl(i))
! 20      continue

  Return

End Subroutine

!-----------------------------------------------------------------------

Subroutine hallcl(index)

  Use :: nlte_type
  Use :: nlte_dim
  Use :: princesa_var, Only: labl1, labu1, fl, wlc
  Implicit None

! it calculates the central wavelength for a rbb transition
! from the atomic data
! character*6 key1,key2
! .. scalar arguments ..
  Integer (i4b) :: index
! ..
! .. local scalars ..
  Real (dp) :: clight, fc, ff1, ff2
  Integer (i4b) :: ii1, ii2
! ..
! .. intrinsic functions ..
! ..
! .. data statements ..
  Data clight/2.997925D18/

! ..
! key1=labl1(index)
! key2=labu1(index)
! nmax=nl
! call ident(key1,labl,nmax,ii1)
! call ident(key2,labl,nmax,ii2)
  ii1 = labl1(index)
  ii2 = labu1(index)
  ff1 = fl(ii1)
  ff2 = fl(ii2)
  fc = abs(ff1-ff2)
  If (fc==0.) Then
    wlc = 1.D20 !                   if coincidentally levels are equal
  Else
    wlc = clight/fc
  End If
  Return

End Subroutine

!-----------------------------------------------------------------------

Subroutine hallcx(index)

  Use :: nlte_type
  Use :: nlte_dim
  Use :: princesa_var, Only: labl2, labu2, labl, nl, labx, nlx, fl, flx, wlc
  Implicit None

! it calculates central wavelength for a rbx transition from the
! atomic data of the levels inolved
! .. scalar arguments ..
  Integer (i4b) :: index
! ..
! ..
! .. local scalars ..
  Real (dp) :: clight, fc, ff1, ff2
  Integer (i4b) :: ii1, ii2
  Character :: key1*6, key2*6
! ..
! .. external subroutines ..
  External :: ident
! ..
! .. intrinsic functions ..
! ..
! .. data statements ..
  Data clight/2.997925D18/
! ..

  key1 = labl2(index)
  key2 = labu2(index)
  Call ident(key1, labl, nl, ii1)

  Call ident(key2, labx, nlx, ii2)
  ff1 = fl(ii1)
  ff2 = flx(ii2)
  fc = abs(ff1-ff2)
  wlc = clight/fc

  Return

End Subroutine

!-----------------------------------------------------------------------

Subroutine ident(key, table, nmax, j)

  Use :: nlte_type
  Use :: nlte_dim
  Implicit None

! it compares the character 'key' with those of table 'table',
! returning in j the value of the corresponding index; if 'key'
! is not in 'table', j=0
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

!-----------------------------------------------------------------------

Function doppler(t, w)

  Use :: nlte_type
  Use :: nlte_dim
  Implicit None

! it calculates the doppler width for a given temperature t
! and atomic weight w.
! the output is an adimensional number that multiplied by the
! central wavelength gives the width in angstroms, or by the
! central frequency, the width in herzts
! .. scalar arguments ..
  Real (dp) :: t, w, doppler
! ..
! .. local scalars ..
  Real (dp) :: amu, cboltz, rccms
! ..
! .. intrinsic functions ..
! ..

  amu = 1.67333D-24

  cboltz = 1.3806D-16

  rccms = 3.3356405D-11

  doppler = rccms*sqrt(2.D0*cboltz*t/(amu*w))

End Function

!ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc

Subroutine abu

  Use :: nlte_type
  Use :: nlte_dim
  Use :: princesa_var, Only: nat, abund, labat
  Implicit None

! this subroutine allows the abundances in the atomic data to be
! written either relative to h or to 1 (always in number of
! particles, not in mass). after this subroutine, abundances are
! always relative to hydrogen, so a certain abundance of h is
! always required (it can be very small, but not 0)
! ..
! .. local scalars ..
  Real (dp) :: abh, summ
  Integer (i4b) :: i, nene
! ..

  summ = 0.D0
  Do i = 1, nat
    summ = summ + abund(i)
  End Do

  If (summ<=1.) Then
    nene = 0
    Do i = 1, nat
      If (labat(i)=='H') nene = i
    End Do

    If (nene==0) Stop 'ERROR. H NOT FOUND - ABU'
    abh = abund(nene)
    Do i = 1, nat
      abund(i) = abund(i)/abh
    End Do
  End If

  Return
End Subroutine
