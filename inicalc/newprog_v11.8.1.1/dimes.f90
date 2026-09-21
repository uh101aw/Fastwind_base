!version 1.7.2 (Nov 2025):
!  polished

!version 1.7.1 (Sept 2021):
!  number of freq. points changed (rounded + 500)
!  [including add. freq. for Stark-broadened transitions]

!version 1.7 (Oct 2020):
!  compatible with gfortran

!version 1.6 (March 2017):
!  one more freq. point at 303.797 (HeII Ly-alpha)

!version 1.5 (Jan 2015):
!  frequency grid changed between 1600 and 911, consistent with nlte.f90
!  checked July 2016 that v. 1.3.4  from Feb 2015 was already correctly
!  included:
!  FRESCAL_0 CHANGED FOR BETTER RESOLUTION

!version 1.4:
!  allows for X-ray treatment (FRESCAL0 adapted for higher frequencies and
!  K-shells

!version 1.3.3 (July 2013):
!  NFESC set to 500 (for new formalsol)

!version 1.3.2 (April 2013):
!  Improved treatment of freq. spacing between 911 and 925 A, &
!  to ensure a maximum separation of 3 A (minimum box size for nsubinter = 120
!  is 3.6 A, can lead to problems with nlambox_ly when many ions included)

!version 1.3.1 (March 2013):
!  number of freq. points changed (rounded + 200)
!  higher resolution in UV range until Lyman edge
!  new default values for nd1 and np1

!version 1.3 until March 2013

Module nlte_param

  Use :: nlte_type
  Implicit None

! fix parameters for dimes etc.

! nlte
! INTEGER(I4B), PARAMETER :: NF2 = 550             !freq. fine grid (HHE)
  Integer (i4b), Parameter :: nf2 = 4000 ! freq. fine grid (large atoms)
  Integer (i4b), Parameter :: levmax = 4
  Integer (i4b), Parameter :: levmin1 = 3
  Integer (i4b), Parameter :: levmin2 = 2
  Integer (i4b), Parameter :: nmaxion = 5 ! maximum number of ions per
! atom (without highest one)
! INTEGER(I4B), PARAMETER :: ND1 = 47              !depth points (for purely
! optical analysis)
! INTEGER(I4B), PARAMETER :: NP1 = 52              !p-rays
! INTEGER(I4B), PARAMETER :: ND1 = 51              !suggested if IR should be
! calculated
! INTEGER(I4B), PARAMETER :: NP1 = 56              !p-rays
! INTEGER(I4B), PARAMETER :: NP1 = 61              ! high resol. p-rays
! INTEGER(I4B), PARAMETER :: ND1 = 91              !suggested if IR should be
! calculated
! INTEGER(I4B), PARAMETER :: NP1 = 96              !p-rays
! INTEGER(I4B), PARAMETER :: ND1 = 61              !suggested if IR should be
! calculated
! INTEGER(I4B), PARAMETER :: NP1 = 66              !p-rays
! INTEGER(I4B), PARAMETER :: ND1 = 53              !suggested if IR should be
! calculated
! INTEGER(I4B), PARAMETER :: NP1 = 59              !p-rays
  Integer (i4b), Parameter :: nd1 = 67 ! depth points for high resol. of
! photosphere
! INTEGER(I4B), PARAMETER :: NP1 = 72 !p-rays for high resol. of photosphere
  Integer (i4b), Parameter :: np1 = 77 ! high resol. p-rays for high resol. of
! photosphere
  Integer (i4b), Parameter :: nfcmf = 21 ! cmf-freq.
! formal
  Integer (i4b), Parameter :: nfformal = 3000 ! depth points (fine grid)
  Integer (i4b), Parameter :: ncformal = 10 ! core rays
  Integer (i4b), Parameter :: nfobs = 161 ! obs. frame freq.
  Integer (i4b), Parameter :: nfesc = 500 ! cmf frame freq.
! stark
  Integer (i4b), Parameter :: nfstark = 80
  Integer (i4b), Parameter :: ntstark = 50
  Integer (i4b), Parameter :: nestark = 50
  Integer (i4b), Parameter :: ntotstark = 10000
End Module

!----------------------------------------------------------------------

Module dimes_var

  Use :: nlte_type
  Implicit None

! dimensiones

  Integer (i4b) :: nat = 0
  Integer (i4b) :: iong = 0
  Integer (i4b) :: nl = 0
  Integer (i4b) :: nlx = 0
  Integer (i4b) :: ns = 0
  Integer (i4b) :: indtt = 0
  Integer (i4b) :: index1 = 0, index2 = 0, index3 = 0, index4 = 0, index5 = 0
  Integer (i4b) :: indtr = 0
  Integer (i4b) :: ndat = 0
  Integer (i4b) :: nprof = 0

! inout

  Integer (i4b) :: iu
  Integer (i4b) :: iou

! xray-treatment
  Integer (i4b), Parameter :: n_kedges = 59 ! in ATOMDAT_NEW/k_shell_data
  Integer (i4b), Parameter :: k_nmin = 3, k_nmax = 8 ! min/max ion. stage to
! be included in freq.
! grid

  Real (dp), Dimension (n_kedges) :: eth, sigma, s, zeff
  Integer (i4b), Dimension (n_kedges) :: z, n

End Module

!----------------------------------------------------------------------

Module prince_var

  Use :: nlte_type
  Implicit None
! ----------------------------------------------------------------------
! atdnin

  Integer (i4b) :: nl = 0, nlx = 0, ns = 0, nat = 0
  Integer (i4b) :: ion, iong = 0
  Integer (i4b), Dimension (:), Allocatable :: la0, la1
  Integer (i4b), Dimension (:), Allocatable :: ixa0, ixa1
  Integer (i4b), Dimension (:), Allocatable :: isa0, isa1
  Integer (i4b), Dimension (:), Allocatable :: nions
  Integer (i4b), Dimension (:), Allocatable :: ifirsl, ifirx, ifirs
  Integer (i4b), Dimension (:), Allocatable :: kl, le, li
  Integer (i4b), Dimension (:), Allocatable :: klx, lix
  Integer (i4b), Dimension (:), Allocatable :: lis, ns0, ns1
! ----------------------------------------------------------------------
! atdnre

  Real (dp), Dimension (:), Allocatable :: zeff, weight, abund
  Real (dp), Dimension (:), Allocatable :: gl, fl, zl
  Real (dp), Dimension (:), Allocatable :: glx, flx
  Real (dp), Dimension (:), Allocatable :: gs0, gs1, ryd, qd
! ----------------------------------------------------------------------
! atomdatc

  Character (6), Dimension (:), Allocatable :: labat
  Character (6), Dimension (:), Allocatable :: labl, lkl
  Character (6), Dimension (:), Allocatable :: labx, lkx
  Character (6), Dimension (:), Allocatable :: labls, kls
! ----------------------------------------------------------------------
! dat

  Integer (i4b) :: ndat = 0
  Real (dp), Dimension (:), Allocatable :: data
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
  Integer (i4b), Dimension (:), Allocatable :: iaux2
  Integer (i4b), Dimension (:), Allocatable :: iform1, numd1, indat1, labl1, &
    labu1
  Integer (i4b), Dimension (:), Allocatable :: iform2, numd2, indat2
  Integer (i4b), Dimension (:), Allocatable :: iform3, numd3, indat3, labl3
  Integer (i4b), Dimension (:), Allocatable :: iform4, numd4, indat4, labl4
  Integer (i4b), Dimension (:), Allocatable :: iform5, numd5, indat5, ipare5
  Integer (i4b), Dimension (:, :), Allocatable :: iaux1
! ----------------------------------------------------------------------
! tranlo

  Logical :: qcl, qdb, qrt
! ----------------------------------------------------------------------
! tranre

  Real (dp) :: ancdo, tl, wlc
  Real (dp), Dimension (:), Allocatable :: dl
  Real (dp), Dimension (:), Allocatable :: frecin, frecfn
! ----------------------------------------------------------------------
! transc

  Character (6), Dimension (:), Allocatable :: labl2, labu2, paren2
  Character (6), Dimension (:), Allocatable :: paren3
  Character (6), Dimension (:), Allocatable :: paren4

End Module

Program dimenatom

! programa para calcular las dimensiones necesarias para los datos
! atomicos

! 14/01/99 changed to f90 by (j.puls)
! 02/03/99 more f90 changes
! 11/03/99 obsolete fortran features removed (j.puls)

! version 1.1 (march 2006) : initial freq. grid from DETAIL input
! no longer read!!! (can be present or absent in file)
! stop enforced if NMAXION is too low

! calls modified version of princesa to read in all data
! these data are used to calculate a pessimistic guess
! of ID_FREC1 (in a modified version of FRESCAL)
! (in order to save disk space when saving the model-files)
! NOTE: ID_FREC2 is still input data, since only used inside
! FASTWIND
  Use :: nlte_type
  Use :: dimes_var, Only: iu, iou
  Use :: ffr_error

  Implicit None

! ..
! .. local scalars ..
  Integer (i4b) :: nf1

  Character :: dc*6, fichero*32, ret*4
! ..
! .. external subroutines ..
  External :: ffracc, ffrcdc, ffrout, freqcon, outsub, rdatom_0, transic_0
! ..
! ..

  iu = 37
  iou = 38
  dc = ':T'

  Open (1, File='ATOM_FILE', Status='OLD')
  Rewind 1
  Read (1, Fmt='(A)') fichero
  Close (1)

  Open (Unit=iu, File=fichero, Form='FORMATTED', Status='OLD')
  Open (Unit=iou, File='control_'//fichero, Status='UNKNOWN')

  Call ffrout(iou, ret)
  If (on_error(ret)) Go To 100

  Call ffracc(iu, ret)
  If (on_error(ret)) Go To 110

  Call ffrcdc(dc)

rloop: Do

    Call rdatom_0(ret)
    If (ret/='RET0' .And. ret/='RET3') Go To 120
    If (ret=='RET3') Go To 130

    Call transic_0(ret)
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
    Stop ' WRONG ERROR CONDITION IN DIMENATOM'
  End Select

! normal exit: end of file
130 Continue
  Write (*, Fmt='(A)') ' FILE IS OVER. SUCCESSFUL!! '

  Close (iu)
  Close (iou)

! read in all data
  Call prince_run_first

! k-shell data read in nlte.f90
! updated routine (read_kshell_data) in lambda_xray.f90

! calculate pessimistic guess for NF1
  Call frescal_0(nf1)

! write nlte_dim.f90
  Call outsub(fichero, nf1)

End Program


!ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc


Subroutine rdatom_0(retchar)
! lee y almacena los datos atomicos
! lista de variables (por orden de appearance):
! nat: numero de atomo
! ion: numero de ion dentro del atomo
! iong: numero de ion en el computo general
  Use :: nlte_type
  Use :: nlte_param, Only: nmaxion
  Use :: dimes_var, Only: nat, iong, nl, nlx, ns
  Use :: ffr_error
  Implicit None

! .. local scalars ..
  Real (dp) :: realvar
  Integer (i4b) :: integ, ion, nc
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
    Stop ' WRONG RETURN IN RDATOM'
  End Select

  nat = nat + 1

  ion = 0
  Call ffrkey(key, nc, ret)
  If (on_error(ret)) Go To 150

  Call ffrnum(realvar, integ, ret)
  If (on_error(ret)) Go To 150

  Call ffrnum(realvar, integ, ret)
  If (on_error(ret)) Go To 150

  Call ffrnum(realvar, integ, ret)
  If (on_error(ret)) Go To 150

100 Continue

  Call ffrkey(key, nc, ret)
  If (on_error(ret)) Go To 150

! ------------- l  levels
  If (key=='L') Then
    ion = ion + 1
    iong = iong + 1

110 Continue
    Call ffrkey(key, nc, ret)
    If (on_error(ret)) Go To 150

    If (key=='0') Then
      Go To 100
    Else
      nl = nl + 1
      Call ffrnum(realvar, integ, ret)
      If (on_error(ret)) Go To 150

      Call ffrnum(realvar, integ, ret)
      If (on_error(ret)) Go To 150

      Call ffrkey(key, nc, ret)
      If (on_error(ret)) Go To 150

      Go To 110


    End If

!   ----------------- x  levels
  Else If (key=='X') Then
120 Continue
    Call ffrkey(key, nc, ret)
    If (on_error(ret)) Go To 150

    If (key=='0') Go To 100
    nlx = nlx + 1
    Call ffrnum(realvar, integ, ret)
    If (on_error(ret)) Go To 150

    Call ffrnum(realvar, integ, ret)
    If (on_error(ret)) Go To 150

    Call ffrkey(key, nc, ret)
    If (on_error(ret)) Go To 150

    Go To 120

!   -------------- s  levels
  Else If (key=='S') Then
130 Continue
    Call ffrkey(key, nc, ret)
    If (on_error(ret)) Go To 150

    If (key/='0') Then
      ns = ns + 1
      Call ffrnum(realvar, integ, ret)
      If (on_error(ret)) Go To 150

      Call ffrnum(realvar, integ, ret)
      If (on_error(ret)) Go To 150

      Call ffrnum(realvar, integ, ret)
      If (on_error(ret)) Go To 150

      Call ffrnum(realvar, integ, ret)
      If (on_error(ret)) Go To 150

      Call ffrnum(realvar, integ, ret)
      If (on_error(ret)) Go To 150

      Call ffrnum(realvar, integ, ret)
      If (on_error(ret)) Go To 150

      Call ffrkey(key, nc, ret)
      If (on_error(ret)) Go To 150

      Go To 130

    End If

    Go To 100

!   --------------------  k  levels
  Else If (key=='K') Then
    nl = nl + 1
    ion = ion + 1
    iong = iong + 1
    Call ffrkey(key, nc, ret)
    If (on_error(ret)) Go To 150

    Call ffrnum(realvar, integ, ret)
    If (on_error(ret)) Go To 150

    Call ffrnum(realvar, integ, ret)
    If (on_error(ret)) Go To 150

    Call ffrkey(key, nc, ret)
    If (on_error(ret)) Go To 150

  End If

! the atom is over.
  retchar = 'RET0'

  If (ion-1>nmaxion) Then
    Print *, ' INCREASE NMAXION TO ', ion - 1
    Stop ' NMAXION TOO LOW!'
  End If

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
    Stop ' WRONG ERROR CONDITION IN RDATOM'
  End Select

End Subroutine


!ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc


Subroutine transic_0(retchar)
! subrutina que lee los datos referentes a las transiciones, los
! almacena y escribe (para control) en un fichero de salida.
! espero poner suficientes explicaciones a lo largo de la subr.
! de modo que este todo lo mas claro posible.
! variables (por orden de aparicion):
! tabla: posibles tipos de transiciones
! qdb: bandera para el caso db (muestreo en anchuras doppler)
! nlong: numero de puntos en el muestreo del perfil
! qcl: bandera para el caso cl (longitud central)
! qrt: bandera para el caso de transicion rbb o rbx
! indtt: indice para llevar la cuenta del total de transiciones
  Use :: nlte_type
  Use :: dimes_var, Only: nprof, indtt, indtr, ndat, index1, index2, index3, &
    index4, index5, iu
  Use :: ffr_error
  Implicit None

! .. local scalars ..
  Real (dp) :: realvar, tl, wlc
  Integer (i4b) :: i, iform, integ, m, nc, ndatos, nlong
  Logical :: qcl, qdb, qrt
  Character :: key*6, pare*6, retchar*4, ret*4
! ..
! .. local arrays ..
  Character :: tabla(8)*6
! ..
! .. external subroutines ..
  External :: ffrkey, ffrnum, ident_0
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
  If (on_error(ret)) Go To 210

  If (key=='TY') Then
    Call ffrkey(key, nc, ret)
    If (on_error(ret)) Go To 210

    Call ident_0(key, tabla, 8, m)
    Call ffrnum(realvar, iform, ret)
    If (on_error(ret)) Go To 210

    Call ffrnum(realvar, ndatos, ret)
    If (on_error(ret)) Go To 210

!   son los numeros de formula y de datos, respectivamente

    Go To 100

!   muestreo en anchuras doppler. el muestreo queda en dl.
  Else If (key=='DB') Then
    qdb = .True.
    Call ffrnum(realvar, nlong, ret)
    If (on_error(ret)) Go To 210

    nprof = max(nprof, nlong)
    Do i = 1, nlong
      Call ffrnum(realvar, integ, ret)
      If (on_error(ret)) Go To 210
    End Do

    Go To 100

!   muestreo en angstroms
  Else If (key=='DL') Then
    qdb = .False.
    Call ffrnum(realvar, nlong, ret)
    If (on_error(ret)) Go To 210

    Do i = 1, nlong
      Call ffrnum(realvar, integ, ret)
      If (on_error(ret)) Go To 210
    End Do

    Go To 100

!   temperatura de linea
  Else If (key=='TL') Then
    Call ffrnum(tl, integ, ret)
    If (on_error(ret)) Go To 210

    Go To 100

!   longitud central de la linea
  Else If (key=='CL') Then
    qcl = .True.

    Go To 100
!   tipo de transicion
  Else
    Call ident_0(key, tabla, 8, m)
    If (m==0) Then
      Go To 220

    Else
      qrt = .False.
      If ((m==1) .Or. (m==3)) qrt = .True.
      Go To (110, 110, 130, 130, 150, 150, 170, 190) m
!     goto calculado, para que haga una cosa u otra segun el tipo de trans.


    End If


  End If

! -------------transiciones rbb o cbb ------------
110 Continue

  Call ffrkey(key, nc, ret)
  If (on_error(ret)) Go To 210

120 Continue
  indtt = indtt + 1
  index1 = index1 + 1

  Call ffrkey(key, nc, ret)
  If (on_error(ret)) Go To 210

! transicion radiativa
  If (qrt) Then
    If (qcl) Then
      Call ffrnum(wlc, integ, ret)
      If (on_error(ret)) Go To 210
    End If

    indtr = indtr + 1
!   nf=nf+nlong
  End If

! lectura de los datos
  ndat = ndat + 1
  ndat = ndat + 1
  Do i = 1, ndatos
    ndat = ndat + 1
    Call ffrnum(realvar, integ, ret)
    If (on_error(ret)) Go To 210
  End Do

! leer nueva keyword
  Call ffrkey(key, nc, ret)
  If (on_error(ret)) Go To 210

  If (key=='0') Then
    qcl = .False.
    Go To 100
  End If


  Go To 120

! --------------- transiciones rbx o cbx ------------
130 Continue
  Call ffrkey(pare, nc, ret)
  If (on_error(ret)) Go To 210

  Call ffrkey(key, nc, ret)
  If (on_error(ret)) Go To 210

140 Continue
  indtt = indtt + 1
  index2 = index2 + 1


  Call ffrkey(key, nc, ret)
  If (on_error(ret)) Go To 210

! transicion radiativa
  If (qrt) Then
    If (qcl) Then
      Call ffrnum(wlc, integ, ret)
      If (on_error(ret)) Go To 210
    End If

  End If

! lectura de datos
  Do i = 1, ndatos
    ndat = ndat + 1
    Call ffrnum(realvar, integ, ret)
    If (on_error(ret)) Go To 210
  End Do

! leer nueva keyword
  Call ffrkey(key, nc, ret)
  If (on_error(ret)) Go To 210

  If (key=='0') Then
    qcl = .False.
    Go To 100
  End If

  Go To 140

! ------------- transiciones  cbs o cbf ------------------
150 Continue
  Call ffrkey(pare, nc, ret)
  If (on_error(ret)) Go To 210

  Call ffrkey(key, nc, ret)
  If (on_error(ret)) Go To 210

160 Continue
  indtt = indtt + 1
  index3 = index3 + 1

! lectura de datos
  Do i = 1, ndatos
    ndat = ndat + 1
    Call ffrnum(realvar, integ, ret)
    If (on_error(ret)) Go To 210
  End Do

! leer nueva keyword
  Call ffrkey(key, nc, ret)
  If (on_error(ret)) Go To 210

  If (key=='0') Go To 100

  Go To 160

! ------------------ transicion  rbf -------------------
170 Continue
  Call ffrkey(pare, nc, ret)
  If (on_error(ret)) Go To 210

  Call ffrkey(key, nc, ret)
  If (on_error(ret)) Go To 210

180 Continue

  indtt = indtt + 1
  index4 = index4 + 1
  Call ffrkey(key, nc, ret)
  If (on_error(ret)) Go To 210

  Call ffrkey(key, nc, ret)
  If (on_error(ret)) Go To 210

! leer datos
  Do i = 1, ndatos
    ndat = ndat + 1
    Call ffrnum(realvar, integ, ret)
    If (on_error(ret)) Go To 210
  End Do
! leer nueva keyword
  Call ffrkey(key, nc, ret)
  If (on_error(ret)) Go To 210

  If (key=='0') Go To 100

  Go To 180

! ----------------- transiciones rff ------------------------
190 Continue

  Call ffrkey(key, nc, ret)
  If (on_error(ret)) Go To 210

200 Continue

  indtt = indtt + 1

  index5 = index5 + 1
! lectura de datos
  Do i = 1, ndatos
    ndat = ndat + 1
    Call ffrnum(realvar, integ, ret)
    If (on_error(ret)) Go To 210
  End Do

! leer nueva keyword
  Call ffrkey(key, nc, ret)
  If (on_error(ret)) Go To 210

  If (key=='0') Go To 100

  Go To 200

! error handling
210 Select Case (ret)
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

220 Continue

  Backspace iu
  retchar = 'RET0'

  Return

End Subroutine


!------------------------------------------------------------------


Subroutine ident_0(key, table, nmax, j)
! compara el caracter key con los de la tabla 'table', y devuelve
! en j el indice que corresponda a key en table, si esta, o 0 en
! otro caso.
! .. scalar arguments ..
  Use :: nlte_type
  Implicit None

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


Subroutine outsub(fichero, nf1)
  Use :: nlte_type
  Use :: nlte_param
  Use :: dimes_var, Only: nat, iong, nl, nlx, ns, index1, index2, index3, &
    index4, index5, indtt, indtr, ndat, nprof

  Implicit None

! For fix parameters (freq., depth points etc., see module nlte_param)

! .. scalar arguments ..
  Integer (i4b) :: nf1
  Character :: fichero*32
! ..
! .. local scalars ..
  Integer (i4b) :: irest, n, ndep
  Character :: alpha*1
! ..
! .. intrinsic functions ..
! ..

  alpha = ' '

  Open (Unit=3, File='nlte_dim.f90', Status='UNKNOWN')

  Write (3, Fmt='(A)') 'Module nlte_dim'
  Write (3, Fmt='(A)') alpha

  Write (3, Fmt='(A)') 'Use nlte_type'
  Write (3, Fmt='(A)') 'IMPLICIT NONE'
  Write (3, Fmt='(A)') alpha

  Write (3, Fmt='(A,2X,A)') '!  DIMENSIONS FROM FILE', fichero
  Write (3, Fmt='(A)') alpha
  Write (3, Fmt='(A)') alpha

  If (nat==0) Stop ' NAT = 0!'
  Write (3, Fmt='(A)') '!    NUMBER OF ATOMS'
  Write (3, Fmt='(A,I2)') 'INTEGER(I4B), PARAMETER :: ID_ATOMS = ', nat
  Write (3, Fmt='(A)') alpha

  Write (3, Fmt='(A)') &
    '!    NUMBER OF FREQUENCIES (coarse mesh, modify in case)'
  Write (3, Fmt='(A,I5)') 'INTEGER(I4B), PARAMETER :: ID_FREC1 = ', nf1
  Write (3, Fmt='(A)') alpha

  Write (3, Fmt='(A)') '!    NUMBER OF FREQUENCIES (fine mesh, input)'
  Write (3, Fmt='(A,I5)') 'INTEGER(I4B), PARAMETER :: ID_FREC2 = ', nf2
  Write (3, Fmt='(A)') alpha

  Write (3, Fmt='(A)') '!    LEVMAX (see frescal)'
  Write (3, Fmt='(A,I5)') 'INTEGER(I4B), PARAMETER :: ID_LEVMA = ', levmax
  Write (3, Fmt='(A)') alpha

  Write (3, Fmt='(A)') '!    LEVMIN1 (see frescal and fresfin)'
  Write (3, Fmt='(A,I5)') 'INTEGER(I4B), PARAMETER :: ID_LEVM1 = ', levmin1
  Write (3, Fmt='(A)') alpha

  Write (3, Fmt='(A)') '!    LEVMIN2(see frescal)'
  Write (3, Fmt='(A,I5)') 'INTEGER(I4B), PARAMETER :: ID_LEVM2 = ', levmin2
  Write (3, Fmt='(A)') alpha

  If (iong==0) Stop ' IONG = 0!'
  Write (3, Fmt='(A)') '!    NUMBER OF IONS'
  Write (3, Fmt='(A,I3)') 'INTEGER(I4B), PARAMETER :: ID_IONES = ', iong
  Write (3, Fmt='(A)') alpha

  If (nl==0) Stop ' NL = 0!'
  Write (3, Fmt='(A)') '!    NUMBER OF L + K LEVELS'
  Write (3, Fmt='(A,I5)') 'INTEGER(I4B), PARAMETER :: ID_LLEVS = ', nl
  Write (3, Fmt='(A)') alpha

  If (nlx/=0) Then
    Print *, ' WARNING!!! X LEVELS FOUND!'
  Else
    nlx = 1
  End If
  Write (3, Fmt='(A)') '!    NUMBER OF X LEVELS '
  Write (3, Fmt='(A,I5)') 'INTEGER(I4B), PARAMETER :: ID_XLEVS = ', nlx
  Write (3, Fmt='(A)') alpha

  If (ns/=0) Then
    Print *, ' WARNING!!! S LEVELS FOUND!'
  Else
    ns = 1
  End If
  Write (3, Fmt='(A)') '!    NUMBER OF S LEVELS'
  Write (3, Fmt='(A,I5)') 'INTEGER(I4B), PARAMETER :: ID_SLEVS = ', ns
  Write (3, Fmt='(A)') alpha

  If (index1==0) Stop ' INDEX1 = 0!'
  Write (3, Fmt='(A)') '!    NUMBER OF RBB + CBB TRANSITIONS'
  Write (3, Fmt='(A,I5)') 'INTEGER(I4B), PARAMETER :: ID_RCBBT = ', index1
  Write (3, Fmt='(A)') alpha

  If (index2/=0) Stop ' RBX AND/OR CBX TRANSITIONS FOUND!'
  index2 = 1
  Write (3, Fmt='(A)') '!    NUMBER OF RBX + CBX TRANSITIONS'
  Write (3, Fmt='(A,I5)') 'INTEGER(I4B), PARAMETER :: ID_RCBXT = ', index2
  Write (3, Fmt='(A)') alpha

  If (index3==0) Stop ' INDEX3 = 0!'
  Write (3, Fmt='(A)') '!    NUMBER OF CBS + CBF TRANSITIONS'
  Write (3, Fmt='(A,I5)') 'INTEGER(I4B), PARAMETER :: ID_CBSFT = ', index3
  Write (3, Fmt='(A)') alpha

  If (index4==0) Stop ' INDEX4 = 0!'
  Write (3, Fmt='(A)') '!    NUMBER OF RBF TRANSITIONS'
  Write (3, Fmt='(A,I5)') 'INTEGER(I4B), PARAMETER :: ID_RBFTR = ', index4
  Write (3, Fmt='(A)') alpha

  If (index5==0) Stop ' INDEX5 = 0!'
  Write (3, Fmt='(A)') '!    NUMBER OF RFF TRANSITIONS'
  Write (3, Fmt='(A,I5)') 'INTEGER(I4B), PARAMETER :: ID_RFFTR = ', index5
  Write (3, Fmt='(A)') alpha

  If (indtt==0) Stop ' INDTT = 0!'
  Write (3, Fmt='(A)') '!    NUMBER OF TOTAL TRANSITIONS'
  Write (3, Fmt='(A,I5)') 'INTEGER(I4B), PARAMETER :: ID_NTRAN = ', indtt
  Write (3, Fmt='(A)') alpha

  If (indtr==0) Stop ' INDTR = 0!'
  Write (3, Fmt='(A)') '!   NUMBER OF RBB + RBX TRANSITIONS'
  Write (3, Fmt='(A,I5)') 'INTEGER(I4B), PARAMETER :: ID_NTTRD = ', indtr
  Write (3, Fmt='(A)') alpha

  If (ndat==0) Stop ' NDAT = 0!'
  Write (3, Fmt='(A)') '!    NUMBER OF DATA'
  Write (3, Fmt='(A,I5)') 'INTEGER(I4B), PARAMETER :: ID_NDATA = ', ndat
  Write (3, Fmt='(A)') alpha

  If (nprof==0) nprof = 1
  Write (3, Fmt='(A)') '!    NUMBER OF FREQ. PER LINE (from input)'
  Write (3, Fmt='(A,I5)') 'INTEGER(I4B), PARAMETER :: ID_NLONG = ', nprof
  Write (3, Fmt='(A)') alpha

! maximum number of ions/atom excluding the highest one
  Write (3, Fmt='(A)') &
    '!    MAX. NUMBER OF IONS PER ATOM EXCL. LAST ONE (set)'
  Write (3, Fmt='(A,I2)') 'INTEGER(I4B), PARAMETER :: ID_KISAT = ', nmaxion
  Write (3, Fmt='(A)') alpha

  Write (3, Fmt='(A)') '!    NUMBER OF DEPTH POINTS (ND, ND1)'
  ndep = nd1
  If (mod(ndep,2)==0) Stop ' NUMBER OF DEPTH POINTS MUST BE ODD!!!'
  Write (3, Fmt='(A,I2)') 'INTEGER(I4B), PARAMETER :: ID_NDEPT =  ', nd1
  Write (3, Fmt='(A)') alpha

  Write (3, Fmt='(A)') '!    NUMBER OF P-RAYS'
  Write (3, Fmt='(A,I2)') 'INTEGER(I4B), PARAMETER :: ID_NPOIN = ', np1
  Write (3, Fmt='(A)') alpha

  Write (3, Fmt='(A)') '!    NUMBER OF CMF-FREQUENCIES'
  Write (3, Fmt='(A,I2)') 'INTEGER(I4B), PARAMETER :: ID_NFCMF = ', nfcmf
  Write (3, Fmt='(A)') alpha
  Write (3, Fmt='(A)') alpha

  Write (3, Fmt='(A)') '!   DIMENSIONS FOR FORMAL SOLUTION'
  Write (3, Fmt='(A)') alpha

  Write (3, Fmt='(A)') '!    NUMBER OF DEPTH POINTS (FINE MESH)'
  Write (3, Fmt='(A,I4)') 'INTEGER(I4B), PARAMETER :: ID_DEPFI = ', nfformal
  Write (3, Fmt='(A)') alpha

  Write (3, Fmt='(A)') '!    NUMBER OF CORE RAYS'
  Write (3, Fmt='(A,I2)') 'INTEGER(I4B), PARAMETER :: ID_CORES = ', ncformal
  Write (3, Fmt='(A)') alpha

  Write (3, Fmt='(A)') '!    NUMBER OF NON-CORE RAYS MODULO 4 (+1)'

  irest = mod(ndep-1, 4)
  If (irest==0) Then

    n = ndep/4
  Else If (irest==2) Then

    n = ndep/4 + 1
  Else

    Stop ' SOMETHING WAS WRONG IN MY PHILOSOPHY'

  End If
  Write (3, Fmt='(A,I2)') 'INTEGER(I4B), PARAMETER :: ID_NOCOR = ', n
  Write (3, Fmt='(A)') alpha

  Write (3, Fmt='(A)') '!    NUMBER OF OBS. FRAME FREQUENCIES'
  Write (3, Fmt='(A,I3)') ' INTEGER(I4B), PARAMETER :: ID_NFOBS = ', nfobs
  Write (3, Fmt='(A)') alpha

  Write (3, Fmt='(A)') '!    NUMBER OF CMF FREQUENCIES (EL. SCAT.)'
  Write (3, Fmt='(A,I3)') ' INTEGER(I4B), PARAMETER :: ID_NFESC = ', nfesc
  Write (3, Fmt='(A)') alpha
  Write (3, Fmt='(A)') alpha

  Write (3, Fmt='(A)') '!    DIMENSIONS FOR STARK EFFECT TREATMENT'
  Write (3, Fmt='(A)') alpha

  Write (3, Fmt='(A)') '!    MAX. NUMBER OF WAVELENGTHS'
  Write (3, Fmt='(A,I2)') 'INTEGER(I4B), PARAMETER :: ID_MAXWW = ', nfstark
  Write (3, Fmt='(A)') alpha

  Write (3, Fmt='(A)') '!    MAX. NUMBER OF TEMPERATURES'
  Write (3, Fmt='(A,I2)') 'INTEGER(I4B), PARAMETER :: ID_MAXTT = ', ntstark
  Write (3, Fmt='(A)') alpha

  Write (3, Fmt='(A)') '!    MAX. NUMBER OF ELECTRON DENSITIES'
  Write (3, Fmt='(A,I2)') 'INTEGER(I4B), PARAMETER :: ID_MAXNE = ', nestark
  Write (3, Fmt='(A)') alpha

  Write (3, Fmt='(A)') '!    MAX. NUMBER IN PERFIL'
  Write (3, Fmt='(A,I5)') 'INTEGER(I4B), PARAMETER :: ID_MAXPS = ', ntotstark
  Write (3, Fmt='(A)') alpha

  Write (3, Fmt='(A)') 'END MODULE nlte_dim'

  Close (3)

  Return
End Subroutine

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

Subroutine prince_run_first

  Use :: nlte_type
! USE nlte_dim
  Use :: prince_var
  Use :: ffr_error
  Implicit None

! main program for reading and storaging atomic data
! variables
! ..
! .. local scalars ..
  Character :: dc*6, fichero*32, ret*4
! ..
! .. external subroutines ..
  External :: abu, ffracc, ffrcdc, ffrout, rdatom, transic
! ..

  Call alloc

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

  Return


End Subroutine


!ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc

Subroutine rdatom(retchar)

  Use :: nlte_type
  Use :: dimes_var, Only: id_iones => iong
  Use :: prince_var, Only: nl, nlx, ns, nat, ion, iong, la0, la1, ixa0, ixa1, &
    isa0, isa1, nions, ifirsl, ifirx, ifirs, kl, le, li, klx, lix, ns0, ns1, &
    lis, zeff, weight, abund, gl, fl, zl, glx, flx, gs0, gs1, ryd, qd, labat, &
    labl, lkl, labx, lkx, labls, kls
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
  Use :: fund_const, Only: clight
  Use :: prince_var, Only: nl, ns, nat, ifirsl, le, li, weight, fl, labl, &
    labls, data, ndat, iu, nlong, indtt, index1, index2, index3, index4, &
    index5, iform1, iform2, iform3, iform4, iform5, numd1, numd2, numd3, &
    numd4, numd5, indat1, indat2, indat3, indat4, indat5, iaux1, labl4, labl1, &
    labu1, labl3, ipare5, qdb, qcl, qrt, tl, ancdo, wlc, dl, frecin, frecfn, &
    labl2, paren2, paren3, paren4, labu2
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
! frecin( ): actual edge frequency  (rbf)
! frecfn( ): edge frequency with respect to next ground
! state

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
  If (on_error(ret)) Go To 270

  If (key=='TY') Then
    Call ffrkey(key, nc, ret)
    If (on_error(ret)) Go To 270
    Call ident(key, tabla, 8, m)

    Call ffrnum(realvar, iform, ret)
    If (on_error(ret)) Go To 270

    Call ffrnum(realvar, ndatos, ret)
    If (on_error(ret)) Go To 270

!   son los numeros de formula y de datos, respectivamente

    Go To 100

!   sampling in doppler widths. the sampling is storage in dl.
  Else If (key=='DB') Then
    qdb = .True.
    Call ffrnum(realvar, nlong, ret)
    If (on_error(ret)) Go To 270

    Do i = 1, nlong
      Call ffrnum(realvar, integ, ret)
      If (on_error(ret)) Go To 270
      dl(i) = realvar
    End Do


    Go To 100

!   sampling in angstroms
  Else If (key=='DL') Then
    qdb = .False.
    Call ffrnum(realvar, nlong, ret)
    If (on_error(ret)) Go To 270

    Do i = 1, nlong
      Call ffrnum(realvar, integ, ret)
      If (on_error(ret)) Go To 270
      dl(i) = realvar
    End Do


    Go To 100

!   line temperature
  Else If (key=='TL') Then
    Call ffrnum(tl, integ, ret)
    If (on_error(ret)) Go To 270
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
      Go To 280

    Else
      qrt = .False.
      If ((m==1) .Or. (m==3)) qrt = .True.
      Go To (110, 110, 150, 150, 170, 170, 200, 240) m
!     goto calculado, para que haga una cosa u otra segun el tipo de trans.


    End If


  End If

! -------------rbb or cbb transitions ------------
110 Continue

  Call ffrkey(key, nc, ret)
  If (on_error(ret)) Go To 270

120 Continue
  Do iqi = 1, nl
    If (labl(iqi)==key) Then
      nn = iqi
      Go To 130
    End If
  End Do

  Stop ' ERROR IN TRANSIC - LABEL L NOT FOUND, BB'

130 Continue
  indtt = indtt + 1
  iaux1(indtt, 1) = m
  index1 = index1 + 1
  iform1(index1) = iform
  numd1(index1) = ndatos + 2
  labl1(index1) = nn

  Call ffrkey(key, nc, ret)
  If (on_error(ret)) Go To 270

  Do iqi = 1, nl
    If (labl(iqi)==key) Then
      nnu = iqi
      Go To 140
    End If
  End Do

  Stop ' ERROR IN TRANSIC - LABEL U NOT FOUND, BB'

140 Continue
  labu1(index1) = nnu
  indat1(index1) = ndat + 1

  iaux1(indtt, 2) = index1

! radiative transition
  If (qrt) Then
    If (qcl) Then
      Call ffrnum(wlc, integ, ret)
      If (on_error(ret)) Go To 270
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
    If (on_error(ret)) Go To 270
    data(ndat) = realvar
  End Do

! read a new keyword
  Call ffrkey(key, nc, ret)
  If (on_error(ret)) Go To 270

  If (key=='0') Then
    qcl = .False.
    Go To 100
  End If

  Go To 120

! --------------- rbx or cbx transitions -----------
150 Continue
  Call ffrkey(pare, nc, ret)
  If (on_error(ret)) Go To 270

  Call ffrkey(key, nc, ret)
  If (on_error(ret)) Go To 270

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
  If (on_error(ret)) Go To 270

  labu2(index2) = key

! radiative transition
  If (qrt) Then
    If (qcl) Then
      Call ffrnum(wlc, integ, ret)
      If (on_error(ret)) Go To 270
    Else
      Call hallcx(index2)
    End If

    Call perfil

  End If

! read data
  Do i = 1, ndatos
    ndat = ndat + 1
    Call ffrnum(realvar, integ, ret)
    If (on_error(ret)) Go To 270
    data(ndat) = realvar
  End Do

! read new keyword
  Call ffrkey(key, nc, ret)
  If (on_error(ret)) Go To 270

  If (key=='0') Then
    qcl = .False.
    Go To 100
  End If

  Go To 160

! ------------- cbs or cbf transitions  ------------------
170 Continue
  Call ffrkey(pare, nc, ret)
  If (on_error(ret)) Go To 270

  Call ffrkey(key, nc, ret)
  If (on_error(ret)) Go To 270

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

  Stop ' ERROR IN TRANSIC - LABEL L OR S NOT FOUND -CBF-CBS'

190 Continue

  labl3(index3) = nn

! read data
  Do i = 1, ndatos
    ndat = ndat + 1
    Call ffrnum(realvar, integ, ret)
    If (on_error(ret)) Go To 270
    data(ndat) = realvar
  End Do

! read new keyword
  Call ffrkey(key, nc, ret)
  If (on_error(ret)) Go To 270

  If (key=='0') Go To 100

  Go To 180

! ------------------ rbf transition  -------------------
200 Continue
  Call ffrkey(pare, nc, ret)
  If (on_error(ret)) Go To 270

  Call ffrkey(key, nc, ret)
  If (on_error(ret)) Go To 270

210 Continue
  indtt = indtt + 1
  iaux1(indtt, 1) = m
  index4 = index4 + 1
  iform4(index4) = iform
  numd4(index4) = ndatos
  indat4(index4) = ndat + 1
  iaux1(indtt, 2) = index4
  paren4(index4) = pare

  Do isi = 1, nl
    If (labl(isi)==key) Then
      nn = isi
      Go To 220
    End If
  End Do

  Stop ' ERROR IN TRANSIC - LABEL L NOT FOUND'

220 Continue
  labl4(index4) = nn
  frecfn(index4) = fl(nn)

! calculation of actual edge frequency (if not ground state
! ionization)

  Do isi = 1, nl
    If (labl(isi)==pare) Then
      nn = isi
      Go To 230
    End If
  End Do

  Stop ' ERROR IN TRANSIC - LABEL PARE NOT FOUND'

230 Continue
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
  If (on_error(ret)) Go To 270

  Call ffrkey(key, nc, ret)
  If (on_error(ret)) Go To 270

  frecin2 = frecfn(index4) + deltae

  frecin(index4) = frecin2/clight
  frecfn(index4) = frecfn(index4)/clight

! a la salida de aqui frecin y frecfn estan en cm-1
! read data
  Do i = 1, ndatos
    ndat = ndat + 1
    Call ffrnum(realvar, integ, ret)
    If (on_error(ret)) Go To 270
    data(ndat) = realvar
  End Do

! read new keyword
  Call ffrkey(key, nc, ret)
  If (on_error(ret)) Go To 270

  If (key=='0') Go To 100

  Go To 210

! ----------------- rff transitions ------------------------
240 Continue

  Call ffrkey(key, nc, ret)
  If (on_error(ret)) Go To 270

250 Continue
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
      Go To 260
    End If
  End Do

  Stop ' ERROR IN TRANSIC - LABEL L NOT FOUND'

260 Continue
  ipare5(index5) = nn
! read data
  Do i = 1, ndatos
    ndat = ndat + 1
    Call ffrnum(realvar, integ, ret)
    If (on_error(ret)) Go To 270
    data(ndat) = realvar
  End Do

! read new keyword
  Call ffrkey(key, nc, ret)
  If (on_error(ret)) Go To 270

  If (key=='0') Go To 100

  Go To 250

! --------------------------------------------------------
! error handling
270 Select Case (ret)
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

280 Continue
  Backspace iu

  retchar = 'RET0'
  Return

End Subroutine

!-----------------------------------------------------------------------

Subroutine perfil

  Use :: nlte_type

  Use :: prince_var, Only: nlong, indtt, indtr, iaux2, qdb, dl

  Implicit None

! calculates the sampling of the atomic profile for radiative
! bound-bound transitions
! ..
! .. local scalars ..

  Integer (i4b) :: i
! ..
! .. local arrays ..
! unused:  Integer (i4b) :: ifptr(id_nttrd), nfptr(id_nttrd)
! ..
! this routine should never be called, since profiles separately treated, THUS

  indtr = indtr + 1
  iaux2(indtr) = indtt
! unused:  nfptr(indtr) = nlong

! unused:  ifptr(indtr) = nf + 1
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
! Use nlte_dim
  Use :: prince_var, Only: labl1, labu1, fl, wlc
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
! Use nlte_dim
  Use :: prince_var, Only: labl2, labu2, labl, nl, labx, nlx, fl, flx, wlc
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
! Use nlte_dim
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
! Use nlte_dim
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
! Use nlte_dim
  Use :: prince_var, Only: nat, abund, labat
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

    If (nene==0) Stop ' ERROR. H NOT FOUND - ABU'
    abh = abund(nene)
    Do i = 1, nat
      abund(i) = abund(i)/abh
    End Do
  End If

  Return
End Subroutine

!ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc

Subroutine alloc

  Use :: nlte_type
  Use :: prince_var
  Use :: dimes_var, Only: id_atoms => nat, id_iones => iong, id_llevs => nl, &
    id_xlevs => nlx, id_slevs => ns, id_data => ndat, id_nttrd => indtr, &
    id_rcbbt => index1, id_rcbxt => index2, id_cbsft => index3, &
    id_rbftr => index4, id_rfftr => index5, id_ntran => indtt, &
    id_nlong => nprof

  Implicit None
! ----------------------------------------------------------------------

  Allocate (la0(id_atoms), la1(id_atoms))
  la0 = 0
  la1 = 0

  Allocate (ixa0(id_atoms), ixa1(id_atoms))
  ixa0 = 0
  ixa1 = 0

  Allocate (isa0(id_atoms), isa1(id_atoms))
  isa0 = 0
  isa1 = 0

  Allocate (nions(id_atoms))
  nions = 0

  Allocate (ifirsl(id_iones), ifirx(id_iones), ifirs(id_iones))
  ifirsl = 1
  ifirx = 1
  ifirs = 1

  Allocate (kl(id_llevs), le(id_llevs), li(id_llevs))
  kl = 0
  le = 0
  li = 0

  Allocate (klx(id_xlevs), lix(id_xlevs))
  klx = 0
  lix = 0

  Allocate (lis(id_slevs), ns0(id_slevs), ns1(id_slevs))
  lis = 0
  ns0 = 0
  ns1 = 0

  Allocate (zeff(id_atoms), weight(id_atoms), abund(id_atoms))
  zeff = 0.
  weight = 0.
  abund = 0.

  Allocate (gl(id_llevs), fl(id_llevs), zl(id_llevs))
  gl = 0.
  fl = 0.
  zl = 0.

  Allocate (glx(id_xlevs), flx(id_xlevs))
  glx = 0.
  flx = 0.

  Allocate (gs0(id_slevs), gs1(id_slevs), ryd(id_slevs), qd(id_slevs))
  gs0 = 0.
  gs1 = 0.
  ryd = 0.
  qd = 0.

  Allocate (labat(id_atoms))
  Allocate (labl(id_llevs), lkl(id_llevs))
  Allocate (labx(id_xlevs), lkx(id_xlevs))
  Allocate (labls(id_slevs), kls(id_slevs))

  Allocate (data(id_data))
  data = 0.

  Allocate (iaux2(id_nttrd))
  Allocate (iform1(id_rcbbt), numd1(id_rcbbt), indat1(id_rcbbt), &
    labl1(id_rcbbt), labu1(id_rcbbt))
  Allocate (iform2(id_rcbxt), numd2(id_rcbxt), indat2(id_rcbxt))
  Allocate (iform3(id_cbsft), numd3(id_cbsft), indat3(id_cbsft), &
    labl3(id_cbsft))
  Allocate (iform4(id_rbftr), numd4(id_rbftr), indat4(id_rbftr), &
    labl4(id_rbftr))
  Allocate (iform5(id_rfftr), numd5(id_rfftr), indat5(id_rfftr), &
    ipare5(id_rfftr))
  Allocate (iaux1(id_ntran,2))
  Allocate (dl(id_nlong))
  Allocate (frecin(id_rbftr), frecfn(id_rbftr))

  Allocate (labl2(id_rcbxt), labu2(id_rcbxt), paren2(id_rcbxt))
  Allocate (paren3(id_cbsft))
  Allocate (paren4(id_rbftr))

  Return
End Subroutine

!ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc

Function igenio(kk, ii)

  Use :: nlte_type

  Use :: prince_var, Only: iong, nions

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
  If (igenio>iong) Stop ' ERROR IN GENERAL ION - IGENIO '

  Return

End Function

!ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc


Subroutine frescal_0(ifre)

! -----  this subroutine provides a pessimistic guess for ID_FREC1
! in a similar spirit as the original version.
! if major changes are done there, they should be applied also here!!!!

! present philosophy requires that EACH level of a DETAILed ion
! is represented with two corresponding frequencies, namely at
! the corresponding edge and at edge-epslon(edge)

! in this version, we use only LEVMAX and ALL ions

  Use :: nlte_type

  Use :: dimes_var, Only: kel => nat, id_rbftr => index4

  Use :: dimes_var, Only: n_kedges, k_nmin, k_nmax, nionk => n, eth

  Use :: nlte_param, Only: levmax, id_frec2 => nf2, kis => nmaxion

  Use :: fund_const

  Use :: prince_var, Only: nat, nions, le, li, zeff, labat, labl, labl4, &
    frecin, paren4

  Implicit None

! .. parameters ..
  Integer (i4b), Parameter :: nfmin = 300
  Logical, Parameter :: optout = .False.

  Integer (i4b) :: ifre

! ..
! .. local scalars ..
  Real (dp) :: edge, eps, fremin, frenew, aaa, emax, dmin, dmin1, dmin2, dm, &
    dnew, dnew1, xmin, xmax, dfl

  Integer (i4b) :: nopa, i, j, k, l, m, n, ihi, iz, isi, levused, levtot, &
    levporg, ii, ilab, leve, levp, lli, if1, nper, imp, ji, inew, il, io, ll, &
    nk


  Character :: lab*6, pare*6
! ..
! .. local arrays ..
  Real (dp) :: flaux(id_rbftr), freedge(id_rbftr), fre(id_frec2)
  Integer (i4b) :: indedge(id_rbftr), iopa(kis*kel), index(id_frec2)
! ..
! .. statement functions ..
  Real (dp) :: epslon
! ..
! .. statement function definitions ..

  epslon(edge) = 5.D0*10.D0**(dble(int(log10(edge)))-6.D0)
! ..
  If (kel/=nat) Stop ' SOMETHING WRONG WITH KEL,NAT'

  If (id_rbftr>999) Stop ' TOO MANY RBF-TRANSITIONS, CODING OF INDEX AFFECTED'

! most abundant element k1


  nopa = 0

  Do k = 1, nat
    Do m = 1, nions(k) - 1
      nopa = nopa + 1
      iopa(nopa) = k*100 + m
    End Do
  End Do

  If (nopa>kis*kel) Stop ' ERROR IN NOPA'

  Print *
  Print *, ' USED IONS FOR FREQUENCY GRID'
  Print *

  Do n = 1, nopa
    ihi = iopa(n)/100
    iz = int(zeff(ihi))
    isi = iopa(n) - 100*ihi
    If (nions(ihi)==isi) Stop ' ERROR IN NIONS OR ISI'
    levused = levmax
    Write (*, Fmt=130) labat(ihi), isi, iz + isi - 1, levused
  End Do

  Print *

  ifre = 0
  levtot = 0
  Do n = 1, nopa

    k = iopa(n)/100
    i = iopa(n) - 100*k

    If (i==nions(k)) Stop ' ERROR IN NIONS OR ISI'
    If (i<=0) Stop ' ERROR IN ION - FRESCAL_0 '

    levporg = 0
    Do ii = 1, id_rbftr
      ilab = labl4(ii)
      If (le(ilab)==k .And. li(ilab)==i) Then
        levporg = levporg + 1
        flaux(levporg) = frecin(ii)

      End If
    End Do

!   reordering of flaux since fl maybe not ordered by energy
!   example: hei

    levtot = levtot + levporg
    Call sort(levporg, flaux)

!   find which edges shall be resolved

!   always levmax levels (assuming levmax < levporg)

    leve = levmax

    levp = levporg

    Do l = 1, levp
      ifre = ifre + 1
      lli = levporg + 1 - l
      Do ii = 1, id_rbftr
        If (flaux(lli)==frecin(ii)) Go To 100
      End Do
      Stop ' FLAUX NOT FOUND IN FRECIN'

100   Continue
      index(ifre) = 100000*k + i*1000 + ii
      If (l<=leve) Then
        indedge(ifre) = 1
      Else
        indedge(ifre) = 0
      End If

!     note that flaux here is in the "wrong' order

      fre(ifre) = flaux(lli)

    End Do

  End Do

  If (levtot>id_rbftr) Stop ' ERROR IN LEVTOT'

  If (ifre>id_rbftr) Stop ' ERROR IN IFRE' ! THIS IS THE NEW ONE

  Print *, ' NUMBER OF CONSIDERED EDGES = ', ifre
  Print *

  Do i = 1, ifre
    eps = epslon(fre(i))
    fre(i+ifre) = fre(i) - eps
    freedge(i) = fre(i)
  End Do

  if1 = ifre
  ifre = 2*ifre

! K-shell edges
  nk = 0

  Do k = 1, n_kedges
    If (nionk(k)<k_nmin) Cycle
    If (nionk(k)>k_nmax) Cycle
    If (eth(k)>2.D7) Cycle
    nk = nk + 1
    ifre = ifre + 1
    fre(ifre) = eth(k)
    ifre = ifre + 1
    eps = epslon(eth(k))
    fre(ifre) = fre(ifre-1) - eps
  End Do

  Print *, ' NUMBER OF CONSIDERED K-SHELL EDGES = ', nk
  Print *

  Call sort(ifre, fre)

  fremin = fre(1)
  i = log10(fremin) - 1
  i = max(i, 0)
  i = 10**i
  fremin = float(int(fremin/i)*i)

  If (fre(1)<=fremin) Stop ' ERROR IN FREMIN'

! OLD VERSION, WITH FRE(1)=FREMIN
! DO I = IFRE + 1,2,-1
! FRE(I) = FRE(I-1)
! END DO

! FRE(1) = FREMIN
! IFRE = IFRE + 1

! NEW VERSION FOR mm-FLUXES
! 9 points before fremin, starting at 1mm = 10 Kayser
  dfl = log10(fremin/10.)/9. !      10 Kayser, 10 points = 9 intervals

  Do i = ifre + 10, 11, -1
    fre(i) = fre(i-10)
  End Do

  fre(1) = 10.
  Do i = 1, 9
    fre(i+1) = fre(1)*10.**(i*dfl)
  End Do

  ifre = ifre + 10

! AAA = 5.D6  !  ( 20 A) sort of minimum
! new version including x-ray treatment
  aaa = 2.D7 !                      (  5 A)

! HYDROGEN ONLY
! changed March 2012, Valparaiso
  If (nat==1 .And. iopa(1)/100==1) aaa = 4.D5 ! 250 A
! might still lead to problems with ifremax for cool stars

  Print *
  Print *, ' FRE(IFRE): ', 1.0D8/fre(ifre)
  Print *, ' MIN      : ', 1.0D8/aaa
  Print *

  emax = aaa
  emax = max(emax, fre(ifre)*1.2D0)

  If (emax>7.5D5) Then !            in case, resolve HeII edge
    ifre = ifre + 1
    fre(ifre) = 7.5D5 !             (133 A)
  End If

  If (emax>1.2D6) Then
    ifre = ifre + 1
    fre(ifre) = 1.2D6 !             (86 A, half way between 133 A and first
!   K-shell at 39A)
  End If

  If (emax>5.D6) Then !             xray treatment
    ifre = ifre + 1
    fre(ifre) = 5.D6 !              (20 A)
  End If

  ifre = ifre + 1
  fre(ifre) = emax

  Call sort(ifre, fre)

! changed to obtain similar freq. grids in all cases
! dmin=log(fre(ifre)/fre(1))/nfmin

  dmin = log(1.D6/1.D3)/nfmin
  nper = ifre - 1

  Do i = 1, nper

    If (fre(i+1)<1.D4) Then !       >10000 A
      dmin1 = dmin*3
    Else If (fre(i+1)<5.D4) Then !  >2000 A
      dmin1 = dmin*2
    Else If (fre(i+1)<1.D8/1600.) Then ! >1600 A
!     changed Jan 27 2015, for better resol between 2000 and 1600 of non-HHe
!     models
!     DMIN1 = DMIN
      dmin1 = dmin/4.
    Else If (fre(i+1)<1.D8/910.) Then ! > 910 A
      dmin1 = dmin/4.

!     minimum separation for lambda < 227 a  is approx. 2  A

    Else If (fre(i+1)>440528.D0) Then ! < 227 A
      dmin2 = 2.D-8*fre(i+1) !      DELTA LAM / LAM
      If (fre(i+1)>5.D6) dmin2 = 0.3
      dmin1 = dmin2
    Else !                          227 (or min) ... 910 A
      dmin1 = 2000./300000. !       (MIN. RESOL = 2000 KM/S)
    End If

    dm = fre(i+1)/fre(i) - 1.D0

!   ----      find out whether edge should be resolved

    imp = 0
    Do ji = 1, if1
      If (fre(i)==freedge(ji) .And. indedge(ji)==1) imp = 1
    End Do

!   ----      nice trick to obtain resolved edges in case

    If (dm>dmin1/4.D0) Then
      inew = int(dm/dmin1) + 1
      dnew = (fre(i+1)-fre(i))/inew

      Do j = 1, inew
        dnew1 = dnew

!       concentration towards edge for important transitioms

        If (imp==1 .And. j==1) Then
          dnew1 = dnew/4.D0
          Do k = 1, 3
            frenew = fre(i) + k*dnew1
            ifre = ifre + 1
            fre(ifre) = frenew
            If (optout) Write (*, Fmt=160) 1.D8/fre(ifre), dnew1, &
              dnew1/fre(ifre)
          End Do

        Else If (imp==1 .And. j==2) Then
          dnew1 = dnew/2.D0
          frenew = fre(ifre) + dnew1
          ifre = ifre + 1
          fre(ifre) = frenew
          If (optout) Write (*, Fmt=160) 1.D8/fre(ifre), dnew1, &
            dnew1/fre(ifre)
        End If

        If (j/=inew) Then
          frenew = fre(i) + j*dnew
          ifre = ifre + 1
          fre(ifre) = frenew
          If (optout) Write (*, Fmt=160) 1.D8/fre(ifre), dnew1, &
            dnew1/fre(ifre)

        End If

      End Do
    Else If (imp==1) Then

!     important edge, but due to near next edge not resolved
!     assume then that next edge is important and so on, until
!     final resolution

      Do ji = 1, if1
        If (fre(i+2)==freedge(ji)) indedge(ji) = 1
      End Do
    End If

  End Do

  Call sort(ifre, fre)

! check resolution for
! hydrogen resonance lines and hei singlet resonance line


! XMIN=1.D0/1030.D-8
  xmin = 1.D0/1025.D-8
  xmax = 1.D0/912.D-8

! following block changed Feb 2015, to be on the safe side,
! and to include Lyman beta (important for Halpha)
  Do i = 1, ifre - 1

    If (fre(i)<=1.D0/1215.D-8 .And. fre(i+1)>1.D0/1215.D-8) Then ! LYMAN ALPHA
      dmin1 = 0.005D0 !             VERY HIGH RESOLUTION
      dm = fre(i+1)/fre(i) - 1.D0
      If (dm<dmin1) Cycle

!     new treatment
    Else If (fre(i)<=xmin .And. fre(i+1)>xmin) Then ! LYMAN BETA
      dmin1 = 0.001D0 !             EVEN HIGHER RESOLUTION
      dm = fre(i+1)/fre(i) - 1.D0
      If (dm<dmin1) Cycle

    Else If (fre(i)>xmin .And. fre(i+1)<=1.D0/925.D-8) Then ! OTHER HYDROGEN
!     RES. LINES
!     DMIN1=0.01
      dmin1 = 0.001
      dm = fre(i+1)/fre(i) - 1.D0
      If (dm<dmin1) Cycle

    Else If (fre(i)>xmin .And. fre(i)<=1.D0/925.D-8 .And. fre(i+1)<=xmax) Then
      dmin1 = 1.5/925. !            max resol = 1.5 A
      dm = fre(i+1)/fre(i) - 1.D0
      If (dm<dmin1) Cycle

    Else If (fre(i)>1.D0/925.D-8 .And. fre(i+1)<=xmax) Then
!     DMIN1=3./925. ! max resol = 3 A
      dmin1 = 1.5/925. !            max resol = 1.5 A
      dm = fre(i+1)/fre(i) - 1.D0
      If (dm<dmin1) Cycle

    Else If (fre(i)<=xmax .And. fre(i+1)>xmax) Then
!     DMIN1=3./912. ! max resol = 3 A
      dmin1 = 1.5/912. !            max resol = 1.5 A
      dm = fre(i+1)/fre(i) - 1.D0
      If (dm<dmin1) Cycle
!     end new treatment

    Else If (fre(i)<=1.D0/584.D-8 .And. fre(i+1)>1.D0/584.D-8) Then ! HEI RES.
!     LINE
!     DMIN1=0.01
      dmin1 = 0.001
      dm = fre(i+1)/fre(i) - 1.D0
      If (dm<dmin1) Cycle

    Else !                          NO PROBLEMS SO FAR, HEII LYMAN ALPHA
!     RESOLVED ANYWAY
      Cycle
    End If

!   INCLUDE ADDITONAL FREQ. POINTS
    inew = int(dm/dmin1) + 1
    dnew = (fre(i+1)-fre(i))/inew

    Do j = 1, inew - 1
      frenew = fre(i) + j*dnew
      ifre = ifre + 1
      fre(ifre) = frenew
      If (optout) Write (*, Fmt=170) 1.D8/fre(ifre), dnew, dnew/fre(ifre)
    End Do

  End Do

! JO CHANGED March 2017
! ONE ADDITIONAL POINT EXACTLY AT HEII LY-ALPHA
  ifre = ifre + 1
  fre(ifre) = 1.D8/303.797
  If (optout) Write (*, Fmt=180) 1.D8/fre(ifre)


  Call sort(ifre, fre)

  If (optout) Then
    Print *, ' WAVELENGTH(A)     FREQ(HZ)'
    Do i = 1, ifre

      Do j = 1, if1
        If (fre(i)==freedge(j)) Go To 110
      End Do

      Write (*, Fmt=140) 1.D8/fre(i), clight*fre(i)
      Go To 120

110   Continue
      il = index(j)
      k = il/100000
      iz = int(zeff(k))
      il = il - 100000*k
      io = il/1000
      ll = il - 1000*io
      lab = labl(labl4(ll))
      pare = paren4(ll)
      Write (*, Fmt=150) 1.D0/fre(i)*1.D8, clight*fre(i), labat(k), io, &
        iz + io - 1, lab, pare
120 End Do
    Print *
  End If

  Print *, ' TOTAL NUMBER OF FREQUENCIES = ', ifre
  Print *
  i = log10(float(ifre)) !          ROUNDING TO NEXT FULL 100
  i = 10**i
  k = ifre/i
  iz = ifre - k*i
  iz = iz/100

! IFRE=K*I+(IZ+1)*100
! IFRE=K*I+(IZ+2)*100
! changed at March 14th 2013, after request (problem) by roberto
! JO Sept 2021
  ifre = k*i + (iz+2)*100 + 300
  Print *, ' USED MAX. NUMBER OF FREQUENCIES = ', ifre
  Print *


  Return

130 Format (1X, A2, I2, ' Z = ', I2, ' MAX.NO.OF RESOLVED EDGES = ', I2)
140 Format (1X, F12.3, 6X, E12.6)
150 Format (1X, F12.3, 6X, E12.6, 3X, A2, I2, ' Z=', I2, '  FROM ', A6, &
    ' TO ', A6)
160 Format (' ADDIT. POINT AT ', F12.3, 6X, F12.3, ' DELTA E / E = ', F5.3)
170 Format (' ADDIT. POINT AT ', F12.3, 6X, F12.3, ' DELTA E / E = ', F5.3, &
    ' (RESONANCE LINES!)')
180 Format (' ADDIT. POINT AT ', F12.3, ' HEII LY_ALPHA')

End Subroutine

!-----------------------------------------------------------------------

Subroutine sort(n, ra)

  Use :: nlte_type
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
