Program totout

! 03/03/99 changed to f90
! july 2003 number of minor modifications
! march 2006 dimension of enionzc
! sept 2006 reading/output of clf
! oct 2014 length of cat (*50) and related
! sept 2015 recalc. of pressure (to account for update of T and Ne)
! nov 20255 polished
  
  Use :: nlte_type
  Use :: nlte_dim
  Implicit None

! generates complete output for considered model
! in the following form

! changed to treat both old model files (without pressure)
! and new model files (with pressure)

! nd,natom,nkis,nlev,nfreq
! r (dim=nd)
! v (dim=nd)
! dvdr (dim=nd)
! rho (dim=nd)
! p (dim=nd)
! t (dim=nd)
! ne (dim=nd)
! tauross (dim=nd)
! coldens (dim=nd)
! clf(dim=nd)
! level names (dim=nlev)
! departure coeff (dim=nlev,nd)
! ionization fractions (dim=nkis+1,natom,nd)
! lambda(dim=nfreq)
! h-nue (dim=nfreq)
! trad_eff(dim=nfreq)


! .. parameters ..
  Integer (i4b), Parameter :: nd = id_ndept
  Integer (i4b), Parameter :: kel = id_atoms, kis = id_kisat
  Integer (i4b), Parameter :: ifretot = id_frec1
  Integer (i4b), Parameter :: nlevs = id_llevs
  Integer (i4b), Parameter :: natom = 30
! ..
! .. local scalars ..
  Real (dp) :: beta, csound, dum, ggrav, rel, rvmin1, sr, summ, teff, vmax, &
    xmloss, xmu, yhe, summass, sumabu, err
  Integer (i4b) :: i, icorr, ifre, j, k, l, nl, nlev, nu, zeffk, ndef, kismax, &
    iddep
  Character :: cat*50, ltepop*60, nltepo*60, out_tot*60
! ..
! .. local arrays ..
  Real (dp) :: a(nlevs), b(nlevs), dep(nlevs), enionnd(kel, kis+1, nd), &
    fnue(ifretot), fre(ifretot), gv(nd), opac(nd, ifretot), p(nd), r(nd), &
    clf(nd), rho(nd), scont(nd, ifretot), taur(nd), temp(nd), trad(ifretot), &
    v(nd), xlam(ifretot), xm(nd), xne(nd), xnelte(nd), xnh(nd), zeff(kel), &
    abund(natom), aweight(natom), pnew(nd)

  Real (dp), Dimension (:, :, :), Allocatable :: enionzc

  Integer (i4b) :: index(nd), nupper(nlevs)
  Integer (i4b), Dimension (:, :), Allocatable :: defined
  Character :: levnam(nlevs)*6

  Logical :: exist

  Data aweight/1.008, 4.003, 6.941, 9.012, 10.811, 12.011, 14.007, 16.000, &
    18.998, 20.180, 22.990, 24.305, 26.982, 28.085, 30.974, 32.066, 35.453, &
    39.948, 39.098, 40.078, 44.956, 47.88, 50.941, 51.996, 54.938, 55.847, &
    58.933, 58.69, 63.546, 65.39/

! ..
! .. external subroutines ..
  External :: atfile
! ..
! .. intrinsic functions ..
! ..

  Print *, 'GIVE IN MODEL NAME (CATALOG)'
  Read (*, Fmt='(A)') cat

  Print *, 'WITH DEPARTURES = 0, WITH OCCNUM = 1'
! READ*, IDDEP
  iddep = 0

  If (iddep/=0 .And. iddep/=1) Stop ' INPUT MUST BE 0 OR 1!'

  nltepo = trim(cat) // '/NLTE_POP'
  ltepop = trim(cat) // '/LTE_POP'

  If (iddep==0) out_tot = trim(cat) // '/OUT_TOT'
  If (iddep==1) out_tot = trim(cat) // '/OUT_TOT_with_occnum'

  Print *, ' LEVEL STRUCTURE'
  Print *

  Open (10, File=out_tot, Status='UNKNOWN', Form='FORMATTED')

  nlev = nlevs
  nl = 0
  Call atfile(zeff, nupper, nlev, nl, levnam)

  zeffk = maxval(zeff)
  kismax = kis + zeffk
  Allocate (enionzc(kel,kismax+1,nd))

  If (nl>nlevs) Stop ' NL > NLEVS'

  Open (2, File=trim(cat)//'/CONT_FORMAL', Status='OLD', Form='UNFORMATTED')
  Read (2) ifre, (fre(i), i=1, ifretot), ((scont(i,j),i=1,nd), j=1, ifretot), &
    ((opac(i,j),i=1,nd), j=1, ifretot), (temp(i), i=1, nd)
  Close (2)

  Write (10, Fmt=*) nd, kel, kismax, nl, ifre

  Open (2, File=trim(cat)//'/MODEL', Status='OLD', Form='UNFORMATTED')
  Read (2, Err=100) teff, ggrav, sr, yhe, xmu, vmax, xmloss, beta, &
    (r(i), i=1, nd), (v(i), i=1, nd), (gv(i), i=1, nd), (rho(i), i=1, nd), &
    (xne(i), i=1, nd), (xnh(i), i=1, nd), (clf(i), i=1, nd), &
    (index(i), i=1, nd), (p(i), i=1, nd)
! ne from hydro model (no T-corrections)
  Go To 110


100 Continue

! to allow reading of very old version
  Print *
  Print *, ' PRESSURE ONLY FOR MU = CONST'
  Print *
  csound = 1.3806D-16/(xmu*1.673D-24)

  Do l = 1, nd
    p(l) = rho(l)*temp(l)*csound
  End Do

110 Continue

  Close (2)

  Print *
  Print *, ' MODEL-PARAMETERS'
  Print *

  Print *, 'TEFF = ', teff, ' LOG G =', log10(ggrav), &
    ' INNERMOST RADIUS POINT = ', sr/6.96D10
  Print *, ' YHE = ', yhe
  Print *, ' MDOT = ', xmloss*3.1557D7/1.989D33, ' VINF = ', vmax*1.D-5, &
    ' BETA = ', beta
  Print *

  Open (2, File=trim(cat)//'/TAU_ROS', Status='OLD', Form='FORMATTED')
  Do l = 1, nd
    Read (2, Fmt=*) taur(l), xnelte(l)
  End Do
  Close (2)

  xm(1) = 0.D0
  Do l = 2, nd
    xm(l) = xm(l-1) + .5D0*(rho(l-1)+rho(l))*(r(l-1)-r(l))
  End Do

  Do l = 1, nd
    xm(l) = xm(l)*sr
  End Do

  Open (1, File='OUT_LASTMODEL', Status='UNKNOWN', Form='FORMATTED')
  Write (1, Fmt='(A)') ' N      R            V           DV/DR &
    &         RHO            P ''          XNH         TEMP    &
    &    TAUROSS      COL.DENS.     CLF'
! 23456789 123456789 123456789 123456789 123456789 123456789 123456789 12
  Do i = 1, nd
    Write (1, Fmt=130) i, r(i), v(i), gv(i), rho(i), p(i), xnh(i), temp(i), &
      taur(i), xm(i), clf(i)
  End Do

  Open (2, File=trim(cat)//'/ENION', Status='OLD', Form='UNFORMATTED')
  Read (2) enionnd, xne
  Close (2)

! JO new Sept 2015
  Inquire (File=trim(cat)//'/ABUNDAN', Exist=exist)
  If (exist) Then

    Open (2, File=trim(cat)//'/ABUNDAN', Status='OLD')
    Read (2, *)
    Read (2, *)
    Read (2, *)

    Do i = 0, natom - 6, 6
!     NORMALIZED TO NH
      Read (2, 120)(k, abund(j+i), j=1, 6)
    End Do
    Close (2)

    sumabu = 0.
    summass = 0.
!   NOTE THAT SUMMASS HERE IS DIFFERENT FROM nlte_approx,
!   DUE TO DIFFERENT NORMALIZATION
    Do k = 1, natom
      sumabu = sumabu + abund(k)
      summass = summass + abund(k)*aweight(k)
    End Do

  Else
!   APPROXIMATE
    sumabu = 1. + yhe
    summass = 1. + 4.*yhe
  End If

! recalculate pressure
  err = 1.D100
  Do l = 1, nd
    xmu = summass/(sumabu+xne(l)/xnh(l)) ! NOTE THAT XNH ONLY APPROX.
    csound = 1.3806D-16/(xmu*1.673D-24)
    pnew(l) = rho(l)*temp(l)*csound
    err = min(err, abs(1.-pnew(l)/p(l)))
  End Do


  Write (10, Fmt='(3(F16.10,2x))')(r(i), i=1, nd)
  Write (10, Fmt='(3(E17.10,2x))')(v(i), i=1, nd)
  Write (10, Fmt='(3(E17.10,2x))')(gv(i), i=1, nd)
  Write (10, Fmt='(3(E17.10,2x))')(rho(i), i=1, nd)
! JO new Sept 2015
  Write (10, Fmt='(3(E17.10,2x))')(pnew(i), i=1, nd)
  Write (10, Fmt='(3(F16.5,2x))')(temp(i), i=1, nd)
  Write (10, Fmt='(3(E17.10,2x))')(xne(i), i=1, nd)
  Write (10, Fmt='(3(E17.10,2x))')(taur(i), i=1, nd)
  Write (10, Fmt='(3(E17.10,2x))')(xm(i), i=1, nd)
  Write (10, Fmt='(A)')(levnam(l), l=1, nl)

  Open (11, File=ltepop, Status='OLD', Access='DIRECT', Recl=nl*8)
  Open (12, File=nltepo, Status='OLD', Access='DIRECT', Recl=nl*8)

  rel = 0.D0
  Do l = 1, nd
    xnelte(l) = xnelte(l)/xne(l)
    rel = max(rel, abs(1.D0-xnelte(l)))
  End Do
  Print *
  Print *, ' MAX DEV. BETWEEN NE(LTE) AND NE(NLTE) = ', rel
  Print *

  Allocate (defined(nl,nd))

  Do i = 1, nd
    Read (11, Rec=i)(a(j), j=1, nl)
    Read (12, Rec=i)(b(j), j=1, nl)
    Where (b==0)
      defined(:, i) = 0
    Elsewhere
      defined(:, i) = 1
    End Where
    Do j = 1, nl
!     print*,i,j,a(j),b(j)
      nu = nupper(j)
      If (b(nu)/=0.) Then
        If (iddep==0) dep(j) = b(j)/a(j)*a(nu)/b(nu)*xnelte(i)
!       for comparison with Adi
!       DEP(J)=B(J)/A(J)*XNELTE(I)
!       for absolute occupation numbers
        If (iddep==1) dep(j) = b(j)
      Else
        dep(j) = 0.D0
      End If
    End Do

    Write (10, Fmt='(3(E17.10,2x))')(dep(j), j=1, nl)
  End Do

  Do j = 1, nl
    ndef = 0
    Do i = 1, nd
      If (defined(j,i)==0) ndef = 1
    End Do
    If (ndef==1) Then
      Print *, 'LEVEL ', levnam(j), ' ', j, ' NOT DEFINED AT ALL DEPTH-POINTS'
!     ELSE
!     PRINT*,'LEVEL ',LEVNAM(J),' ',J,' DEFINED AT ALL DEPTH-POINTS'
    End If
  End Do

! note: enionzc corrected for effective charge

lloop: Do l = 1, nd
kloop: Do k = 1, kel
      summ = 0.D0
      zeffk = int(zeff(k))
      Do i = 1, kis + 1
        summ = summ + enionnd(k, i, l)
      End Do
      Do i = 1, zeffk
        enionzc(k, i, l) = 0.D0
      End Do
      Do i = zeffk + 1, kismax + 1
        icorr = i - zeffk
        enionzc(k, i, l) = enionnd(k, icorr, l)/summ
      End Do
    End Do kloop
    Write (10, Fmt='(3(E17.10,2x))')((enionzc(k,i,l),i=1,kismax+1), k=1, kel)
  End Do lloop

  Open (2, File=trim(cat)//'/FLUXCONT', Status='OLD', Form='FORMATTED')
  Read (2, Fmt=*)
  Do k = 1, ifre
    Read (2, Fmt=*) i, xlam(i), fnue(i), dum, trad(i)
  End Do

  Print *
  Print *, ' MAX DEV. BETWEEN P AND PNEW = ', err

  Read (2, Fmt=*) rvmin1
  Print *
  Print *, ' EFFECTIVE RADIUS = ', sr*rvmin1/6.96D10
  Close (2)

  Write (10, Fmt='(3(E17.10,2x))')(xlam(i), i=1, ifre)
  Write (10, Fmt='(3(F16.10,2x))')(fnue(i), i=1, ifre)
  Write (10, Fmt='(3(F16.5,2x))')(trad(i), i=1, ifre)
! this is new (appended in order to allow reading of older files)
  Write (10, Fmt='(3(E17.10,2x))')(clf(i), i=1, nd)


  Close (10)
  Close (11)
  Close (12)

  If (iddep==0) Print *, 'FILE OUT_TOT IN CATALOG ', trim(cat), ' WRITTEN'
  If (iddep==1) Print *, 'FILE OUT_TOT_with_occnum IN CATALOG ', trim(cat), &
    ' WRITTEN'
120 Format (1X, 10(I4,1X,E8.2))

130 Format (1X, I2, 10(1X,G12.4))

End Program

!ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc

Subroutine atfile(zeff, nupper, nlev, nl, levnam)

  Use :: nlte_type
  Use :: nlte_dim
  Use :: ffr_error

  Implicit None

! .. scalar arguments ..
  Integer (i4b) :: nl, nlev
! ..
! .. array arguments ..
  Real (dp) :: zeff(id_atoms)
  Integer (i4b) :: nupper(nlev)
  Character :: levnam(nlev)*6
! ..
! .. scalars in common ..
  Integer (i4b) :: iou, iu, nat
! ..
! .. local scalars ..
  Integer (i4b) :: nprev
  Character :: dc*6, fichero*20, ret*4
! ..
! .. external subroutines ..
  External :: ffracc, ffrcdc, ffrout, rdatom
! ..
! .. common blocks ..
  Common /inout/iu, iou, nat
! ..
  iu = 37
  iou = 38
  dc = ':T'

  nat = 0
  Open (1, File='ATOM_FILE', Status='OLD', Form='FORMATTED')
  Read (1, Fmt='(A)') fichero
  Close (1)

  Open (Unit=iu, File=fichero, Form='FORMATTED', Status='OLD')
  Open (Unit=iou, File='control.dat', Status='UNKNOWN')

  Call ffrout(iou, ret)
  If (on_error(ret)) Go To 100

  Call ffracc(iu, ret)
  If (on_error(ret)) Go To 110

  Call ffrcdc(dc)

  nprev = 1

rloop: Do

    Call rdatom(zeff, nprev, nupper, nlev, nl, levnam, ret)
    If (ret/='RET0' .And. ret/='RET3') Go To 120
    If (ret=='RET3') Go To 130

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
    Stop ' WRONG ERROR CONDITION IN DIMENATOM'
  End Select

! normal exit: end of file
130 Continue
  Write (*, Fmt='(A)') ' FILE IS OVER. SUCCESSFUL!! '

  Close (iu)
  Close (iou)

  Return

End Subroutine


!ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc


Subroutine rdatom(zeff, nprev, nupper, nlev, nl, levnam, retchar)

  Use :: nlte_type
  Use :: nlte_dim
  Use :: ffr_error

  Implicit None

! lee y almacena los datos atomicos
! lista de variables (por orden de appearance):
! nat: numero de atomo
! ion: numero de ion dentro del atomo
! iong: numero de ion en el computo general
! .. scalar arguments ..
  Integer (i4b) :: nl, nlev, nprev
! ..
! .. array arguments ..
  Real (dp) :: zeff(id_atoms)
  Integer (i4b) :: nupper(nlev)
  Character :: levnam(nlev)*6
! ..
! .. scalars in common ..
  Integer (i4b) :: iou, iu, nat
! ..
! .. local scalars ..
  Real (dp) :: realvar
  Integer (i4b) :: ii, integ, nc, ns
  Character :: label*4, key*6, retchar*4, ret*4
! ..
! .. external subroutines ..
  External :: ffrkey, ffrloc, ffrnum
! ..
! .. common blocks ..
  Common /inout/iu, iou, nat
! ..

  label = 'ATOM'
  Call ffrloc(label, ret)
  Select Case (ret)
  Case ('RET1')
    Go To 130
  Case ('RET2')
    Go To 140
  Case ('RET0')
    Continue
  Case Default
    Stop ' WRONG RETURN IN RDATOM'
  End Select

  nat = nat + 1
  Call ffrkey(key, nc, ret)
  If (on_error(ret)) Go To 140

  Call ffrnum(realvar, integ, ret)
  If (on_error(ret)) Go To 140
  zeff(nat) = realvar

  Print *, ' ATOM = ', key, ' ZEFF = ', zeff(nat)
  Call ffrnum(realvar, integ, ret)
  If (on_error(ret)) Go To 140

  Call ffrnum(realvar, integ, ret)
  If (on_error(ret)) Go To 140

100 Continue
  Call ffrkey(key, nc, ret)
  If (on_error(ret)) Go To 140

! ------------- l  levels

  If (key=='L') Then
    Print *, 'L LEVELS'
    Do ii = nprev, nl
      nupper(ii) = nl + 1
    End Do
    nprev = nl + 1

110 Continue
    Call ffrkey(key, nc, ret)
    If (on_error(ret)) Go To 140

    If (key=='0') Then
      Go To 100
    Else
      nl = nl + 1
      Print *, ' LEVEL: ', key, ' CORRESPONDING NO.', nl
      levnam(nl) = key
      Call ffrnum(realvar, integ, ret)
      If (on_error(ret)) Go To 140

      Call ffrnum(realvar, integ, ret)
      If (on_error(ret)) Go To 140

      Call ffrkey(key, nc, ret)
      If (on_error(ret)) Go To 140

      Go To 110

    End If

!   ----------------- x  levels

  Else If (key=='X') Then

    Stop 'X-LEVEL FOUND'

!   2000            call ffrkey(key,nc,*9990,*9999)
!   if (key.eq.'0') goto 100
!   nlx=nlx+1
!   call ffrnum(REALVAR, integ,*9990,*9999)
!   call ffrnum(REALVAR,integ,*9990,*9999)
!   call ffrkey(key,nc,*9990,*9999)
!   goto 2000


!   -------------- s  levels

  Else If (key=='S') Then

120 Continue
    Call ffrkey(key, nc, ret)
    If (on_error(ret)) Go To 140

    If (key/='0') Then
      ns = ns + 1
      Call ffrnum(realvar, integ, ret)
      If (on_error(ret)) Go To 140

      Call ffrnum(realvar, integ, ret)
      If (on_error(ret)) Go To 140

      Call ffrnum(realvar, integ, ret)
      If (on_error(ret)) Go To 140

      Call ffrnum(realvar, integ, ret)
      If (on_error(ret)) Go To 140

      Call ffrnum(realvar, integ, ret)
      If (on_error(ret)) Go To 140

      Call ffrnum(realvar, integ, ret)
      If (on_error(ret)) Go To 140

      Call ffrkey(key, nc, ret)
      If (on_error(ret)) Go To 140

      Go To 120

    End If

    Go To 100

!   --------------------  k  levels

  Else If (key=='K') Then
    Print *, 'K LEVEL'
    Do ii = nprev, nl + 1
      nupper(ii) = nl + 1
    End Do

    nprev = nl + 2
    nl = nl + 1
    Call ffrkey(key, nc, ret)
    If (on_error(ret)) Go To 140
    Print *, ' LEVEL: ', key, ' CORRESPONDING NO.', nl
    levnam(nl) = key

    Call ffrnum(realvar, integ, ret)
    If (on_error(ret)) Go To 140

    Call ffrnum(realvar, integ, ret)
    If (on_error(ret)) Go To 140

    Call ffrkey(key, nc, ret)
    If (on_error(ret)) Go To 140

  End If
! the atom is over.


  If (nlev<nl) Stop 'MORE THAN NLEVS LEVELS IN FILE!!'

  retchar = 'RET0'
  Return

! no hay mas atomos que leer. ya se ha leido todo el fichero. es
! cuando acaba el programa principal (prince)
130 Continue

  retchar = 'RET3'
  Return

! error handling
140 Select Case (ret)
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
