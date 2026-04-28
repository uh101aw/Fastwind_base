Module ffr_error
Contains
  Logical Function on_error(ret)
    Implicit None
    Character (4) :: ret

    on_error = ret /= 'RET0'
    Return
  End Function

End Module
!ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
!cccc                                                              ccccc
!cccc   free-field read routines : source library                  ccccc
!cccc                                                              ccccc
!ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
!cccc
!cccc   status : 03/03/87
!cccc
!cccc   revison record:
!cccc      date: 02/03/87   programmer: d.husfeld
!cccc            function 'ffrlin' added
!cccc      date: 03/03/87   programmer: d.husfeld
!cccc            subroutines 'ffrkey', 'ffrnch', 'ffrnum' modified
!cccc            to accelerate execution
!cccc      date: 31/03/92   programmer: e. santolaya rey
!cccc            all integer variables changed into implicit type (i-n).
!cccc            implicit double precision (a-h,o-z)
!cccc      date: 14/01/99 polished and changed to f90 (J.Puls)
!cccc      date: 20/01/99 lower case character set included (J.Puls)
!cccc      date: 10/03/99 obsolete fortran features removed (J.Puls)
!cccc      date: march 2006 underscore allowed in key character set (J.Puls)

Subroutine ffracc(iun, retchar)

! docs ffracc
! free-field read : access new input unit
! docc    call ffracc ( iun,retchar )
! iun    : an integer variable or integer array element which
! absolute value identifies the unit number
! to be used for further input if it is in the range
! from 1 to 999. the referenced unit must already
! be opened for sequential and formatted access.
! if the value of 'iun' is negative the referenced
! unit will be rewound; if it is positive the unit
! will be not repositioned.
! if the value of 'iun' is zero when ffracc is
! called, 'iun' receives the current input unit
! number and no further action is taken.
! ret1   : a statement label in the calling program to which
! control is returned from ffracc in case of an error
! docr    latest revision : 05/11/86
! docp    programmer      : d.husfeld
! doce

! free-field read : control variables                        ffrctrl
! latest revision : 10/10/86                                 ffrctrl
! iout = output unit number                                  ffrctrl
! iver = verification level                                  ffrctrl
! nor = no run indicator                                     ffrctrl
! cdi = comment delimiting indicator                         ffrctrl
! sdi = string delimiting indicator                          ffrctrl
! cdc = comment delimiting character                         ffrctrl
! sdc = string delimiting character                          ffrctrl
! int = interrupt character                                  ffrctrl
! kcs = keyword character set                                ffrctrl
! kcl = keyword character set length                         ffrctrl
! free-field read : error status table                       ffrerrs
! latest revision : 26/09/86                                 ffrerrs
! errprg = name of routine that detected the error           ffrerrs
! errmsg = error message                                     ffrerrs
! free-field read : file status variables                    ffrstat
! latest revision : 05/11/86                                 ffrstat
! inp = input unit number                                    ffrstat
! icol = column number in line                               ffrstat
! iseq = line sequence number                                ffrstat
! eof = end-of-file indicator                                ffrstat
! lin = line image                                           ffrstat
! free-field read : unit table                               ffrutab
! latest revision : 05/11/86                                 ffrutab
! nun = number of units in table                             ffrutab
! inpt= input unit number table                              ffrutab
! icolt= column number table                                 ffrutab
! iseqt= line sequence number table                          ffrutab
! lint= line image table                                     ffrutab

! free-field read : error status clearance                   ffreclr
! latest revision : 26/09/86                                 ffreclr
! .. parameters ..
  Use :: nlte_type
  Implicit None

  Character (6) :: prgnam
  Parameter (prgnam='FFRACC')
  Integer (i4b) :: kclm
  Parameter (kclm=64)
  Integer (i4b) :: mun
  Parameter (mun=10)
! ..
! .. scalar arguments ..
  Integer (i4b) :: iun
! ..
! .. scalars in common ..
  Integer (i4b) :: icol, inp, iout, iseq, iver, kcl, nun
  Logical :: cdi, eof, nor, sdi
  Character :: cdc*1, int*1, sdc*1, errprg*6, errmsg*40, lin*80, kcs*(kclm)
! ..
! .. arrays in common ..
  Integer (i4b) :: icolt(mun), inpt(mun), iseqt(mun)
  Character :: lint(mun)*80
! ..
! .. local scalars ..
  Integer (i4b) :: i, ierr
  Logical :: rew, uop
  Character :: acmode*10, fmmode*10, retchar*4
! ..
! .. intrinsic functions ..
! ..
! .. common blocks ..
  Common /ffrct1/iout, iver, nor, cdi, sdi, kcl
  Common /ffrct2/cdc, sdc, int, kcs
  Common /ffrerr/errprg, errmsg
  Common /ffrst1/inp, icol, eof, iseq
  Common /ffrst2/lin
  Common /ffrut1/nun, inpt, icolt, iseqt
  Common /ffrut2/lint
! ..
  errprg = ' '
  errmsg = ' '

! .....return current input unit if un is zero
  If (iun==0) Then
    iun = inp
    retchar = 'RET0'
    Return
  End If

! .....save current file status if unit still exists in table
  Do i = 1, nun
    If (inpt(i)==inp) Then
      icolt(i) = icol
      lint(i) = lin
      iseqt(i) = iseq

    End If
  End Do

! .....check iun
  rew = iun < 0
  iun = iabs(iun)
  If (iun>999) Then
    errprg = prgnam
    Write (errmsg, Fmt='(A,I10)') 'NON-ACCEPTABLE UNIT NUMBER :', iun

    Go To 110

  End If
  Inquire (Unit=iun, Iostat=ierr, Opened=uop, Access=acmode, Form=fmmode)
  If ((.Not. uop) .Or. (ierr/=0) .Or. (acmode/='SEQUENTIAL') .Or. &
    (fmmode/='FORMATTED')) Then
    errprg = prgnam
    Write (errmsg, Fmt='(A,I3,A)') 'UNIT #', iun, ' NOT OPENED PROPERLY'

    Go To 110
  End If

! .....access new unit
  Do i = 1, nun
    If (inpt(i)==iun) Then
!     ...unit already exist in table
      inp = iun
      icol = icolt(i)
      iseq = iseqt(i)
      lin = lint(i)
      eof = .False.

      Go To 100

    End If
  End Do
! ...unit is not contained in table
  If (nun>=mun) Then
    errprg = prgnam
    errmsg = 'UNIT BUFFER IS FULL'

    Go To 110

  End If
  nun = nun + 1
  inpt(nun) = iun
  If (inp/=iun) Then
    inp = iun
    lin = ' '
    icol = 81
    iseq = 0
    eof = .False.

  End If
100 Continue

! .....rewind unit if requested
  If (rew) Then
    Rewind inp
    lin = ' '
    icol = 81
    iseq = 0
    eof = .False.
  End If

! .....verify
  If (iver>=1) Write (iout, Fmt='(A,I3,A)', Iostat=ierr) ' ** UNIT #', iun, &
    ' ACCESSED'
  retchar = 'RET0'
  Return

! .....error exit
110 Continue
  retchar = 'RET1'
  Return

End Subroutine
Block Data ffrbdt

! docs ffrbdt
! free-field read : define initial values
! docr    latest revision : 20/01/99
! docp    programmer      : d.husfeld/chr.wiethaus
! doce

! free-field read : control variables                        ffrctrl
! latest revision : 10/10/86                                 ffrctrl
! iout = output unit number                                  ffrctrl
! iver = verification level                                  ffrctrl
! nor = no run indicator                                     ffrctrl
! cdi = comment delimiting indicator                         ffrctrl
! sdi = string delimiting indicator                          ffrctrl
! cdc = comment delimiting character                         ffrctrl
! sdc = string delimiting character                          ffrctrl
! int = interrupt character                                  ffrctrl
! kcs = keyword character set                                ffrctrl
! kcl = keyword character set length                         ffrctrl
! free-field read : error status table                       ffrerrs
! latest revision : 26/09/86                                 ffrerrs
! errprg = name of routine that detected the error           ffrerrs
! errmsg = error message                                     ffrerrs
! free-field read : file status variables                    ffrstat
! latest revision : 05/11/86                                 ffrstat
! inp = input unit number                                    ffrstat
! icol = column number in line                               ffrstat
! iseq = line sequence number                                ffrstat
! eof = end-of-file indicator                                ffrstat
! lin = line image                                           ffrstat
! free-field read : unit table                               ffrutab
! latest revision : 05/11/86                                 ffrutab
! nun = number of units in table                             ffrutab
! inpt= input unit number table                              ffrutab
! icolt= column number table                                 ffrutab
! iseqt= line sequence number table                          ffrutab
! lint= line image table                                     ffrutab
! .. parameters ..
  Use :: nlte_type
  Implicit None

  Integer (i4b) :: kclm
  Parameter (kclm=64)
  Integer (i4b) :: mun
  Parameter (mun=10)
! ..
! .. scalars in common ..
  Integer (i4b) :: icol, inp, iout, iseq, iver, kcl, nun
  Logical :: cdi, eof, nor, sdi
  Character :: cdc*1, int*1, sdc*1, errprg*6, errmsg*40, lin*80, kcs*(kclm)
! ..
! .. arrays in common ..
  Integer (i4b) :: icolt(mun), inpt(mun), iseqt(mun)
  Character :: lint(mun)*80
! ..
! .. common blocks ..
  Common /ffrct1/iout, iver, nor, cdi, sdi, kcl
  Common /ffrct2/cdc, sdc, int, kcs
  Common /ffrerr/errprg, errmsg
  Common /ffrst1/inp, icol, eof, iseq
  Common /ffrst2/lin
  Common /ffrut1/nun, inpt, icolt, iseqt
  Common /ffrut2/lint
! ..
! .. data statements ..

! ffrctrl

! ffrerrs

! ffrstat

! ffrutab
  Data iout/0/, iver/0/, nor/.True./, cdi/.False./, sdi/.False./, cdc/' '/, &
    sdc/' '/, int/':'/

! fixing the problem when using small letters in formal_input.
! because of limitation of columns in fortran and due to portability
! i had to devide 'kcs' in two pieces.
! 'kcl' is set to '62'.

! 13/3/1998 chr. wiethaus

  Data kcs(1:26)/'ABCDEFGHIJKLMNOPQRSTUVWXYZ'/
  Data kcs(27:63)/'abcdefghijklmnopqrstuvwxyz0123456789_'/
  Data kcl/63/

  Data errprg/' '/, errmsg/' '/
  Data inp/0/, icol/81/, iseq/0/, eof/.False./, lin/' '/
  Data nun/0/, inpt/mun*0/, icolt/mun*81/, iseqt/mun*0/, lint/mun*' '/
! ..

End Block Data
Subroutine ffrcdc(dc)

! docs ffrcdc
! free-field read : define comment delimiting character
! docc    call ffrcdc ( dc )
! dc     : a character variable that defines in its first
! position the new string delimiting character.
! however, if dc contains in its first five positions
! the string 'clear', the delimiting character
! becomes undefined.
! docr    latest revision : 10/10/86
! docp    programmer      : d.husfeld
! doce

! free-field read : control variables                        ffrctrl
! latest revision : 10/10/86                                 ffrctrl
! iout = output unit number                                  ffrctrl
! iver = verification level                                  ffrctrl
! nor = no run indicator                                     ffrctrl
! cdi = comment delimiting indicator                         ffrctrl
! sdi = string delimiting indicator                          ffrctrl
! cdc = comment delimiting character                         ffrctrl
! sdc = string delimiting character                          ffrctrl
! int = interrupt character                                  ffrctrl
! kcs = keyword character set                                ffrctrl
! kcl = keyword character set length                         ffrctrl
! free-field read : error status table                       ffrerrs
! latest revision : 26/09/86                                 ffrerrs
! errprg = name of routine that detected the error           ffrerrs
! errmsg = error message                                     ffrerrs

! free-field read : error status clearance                   ffreclr
! latest revision : 26/09/86                                 ffreclr
! .. parameters ..
  Use :: nlte_type
  Implicit None

  Integer (i4b) :: kclm
  Parameter (kclm=64)
! ..
! .. scalar arguments ..
  Character :: dc*(*)
! ..
! .. scalars in common ..
  Integer (i4b) :: iout, iver, kcl
  Logical :: cdi, nor, sdi
  Character :: cdc*1, int*1, sdc*1, errprg*6, errmsg*40, kcs*(kclm)
! ..
! .. local scalars ..
  Integer (i4b) :: ierr
! ..
! .. intrinsic functions ..
! ..
! .. common blocks ..
  Common /ffrct1/iout, iver, nor, cdi, sdi, kcl
  Common /ffrct2/cdc, sdc, int, kcs
  Common /ffrerr/errprg, errmsg
! ..
  errprg = ' '
  errmsg = ' '

  If ((len(dc)>=5) .And. (dc(1:5)=='CLEAR')) Then
    cdi = .False.
    cdc = ' '

    If (iver>=1) Write (iout, Fmt='(A)', Iostat=ierr) &
      ' ** COMMENT DELIMITER CLEARED'
  Else
    cdi = .True.
    cdc = dc(1:1)
    If (iver>=1) Write (iout, Fmt='(A)', Iostat=ierr) &
      ' ** COMMENT DELIMITER IS SET TO ''' // cdc // ''''

  End If
  Return

End Subroutine
Subroutine ffrdmp

! docs ffrdmp
! free-field read : dump current status
! docc    call ffrdmp
! docr    latest revision : 05/11/86
! docp    programmer      : d.husfeld
! doce

! free-field read : control variables                        ffrctrl
! latest revision : 10/10/86                                 ffrctrl
! iout = output unit number                                  ffrctrl
! iver = verification level                                  ffrctrl
! nor = no run indicator                                     ffrctrl
! cdi = comment delimiting indicator                         ffrctrl
! sdi = string delimiting indicator                          ffrctrl
! cdc = comment delimiting character                         ffrctrl
! sdc = string delimiting character                          ffrctrl
! int = interrupt character                                  ffrctrl
! kcs = keyword character set                                ffrctrl
! kcl = keyword character set length                         ffrctrl
! free-field read : error status table                       ffrerrs
! latest revision : 26/09/86                                 ffrerrs
! errprg = name of routine that detected the error           ffrerrs
! errmsg = error message                                     ffrerrs
! free-field read : file status variables                    ffrstat
! latest revision : 05/11/86                                 ffrstat
! inp = input unit number                                    ffrstat
! icol = column number in line                               ffrstat
! iseq = line sequence number                                ffrstat
! eof = end-of-file indicator                                ffrstat
! lin = line image                                           ffrstat
! .. parameters ..
  Use :: nlte_type
  Implicit None

  Integer (i4b) :: kclm
  Parameter (kclm=64)
! ..
! .. scalars in common ..
  Integer (i4b) :: icol, inp, iout, iseq, iver, kcl
  Logical :: cdi, eof, nor, sdi
  Character :: cdc*1, int*1, sdc*1, errprg*6, errmsg*40, lin*80, kcs*(kclm)
! ..
! .. local scalars ..
  Integer (i4b) :: i, io, nun
  Logical :: lcl, lop
  Character :: word1*7, word2*7
! ..
! .. common blocks ..
  Common /ffrct1/iout, iver, nor, cdi, sdi, kcl
  Common /ffrct2/cdc, sdc, int, kcs
  Common /ffrerr/errprg, errmsg
  Common /ffrst1/inp, icol, eof, iseq
  Common /ffrst2/lin
! ..
! .. data statements ..
  Data lop/.False./, lcl/.False./
! ..

! determine output unit for dump
  If (iout/=0) Then
!   ...if output unit is defined for ffr-routines then use that

    io = iout
  Else
    Inquire (File='OUTPUT', Opened=lop, Number=nun)
    If (lop) Then
!     ...if file 'output' is opened use its unit number

      io = nun
    Else
!     search free unit number and open file 'output' via this unit
      Do i = 1, 999
        Inquire (Unit=i, Opened=lop)
        If (lop) Go To 100
        io = i
        Open (Unit=io, File='OUTPUT', Status='UNKNOWN')
        lcl = .True.

        Go To 110
100   End Do
!     ...output file can not be accessed, return
      Return

110   Continue

    End If
  End If

! .....print error status
  If (errmsg==' ') Then

    Write (io, Fmt='(A)') ' NO ERROR DETECTED'
  Else
    Write (io, Fmt='(A)') ' *** ERROR DETECTED BY ' // errprg // ' : ' // &
      errmsg
  End If

! .....print input file status
  If (eof) Then

    Write (io, Fmt='(A,I3,A,I5,A)') ' CURRENT INPUT UNIT #', inp, &
      ' AT EOF (FOUND AFTER LINE #', iseq, ')'
  Else
    Write (io, Fmt='(A,I3,A,I5)') ' CURRENT INPUT UNIT #', inp, ' AT LINE #', &
      iseq
    If (icol<=40) Then

      Write (io, Fmt='(41A1,A,I2)')(' ', i=0, icol-1), '!', &
        (' ', i=icol+1, 40), ' LINE AT COLUMN #', icol
    Else
      Write (io, Fmt='(A17,I2,21X,50A1)') ' LINE AT COLUMN #', icol, &
        (' ', i=41, icol-1), '!'

    End If
    Write (io, Fmt='(1X,A)') lin

  End If

! .....print control variables
  Write (io, Fmt='(A,I1,A,I3)') ' VERIFICATION LEVEL ', iver, &
    ' , PROTOICOL DIRECTED TO UNIT #', iout
  word1 = 'NOT SET'
  word2 = 'NOT SET'
  If (cdi) word1 = '''' // cdc // ''''
  If (sdi) word2 = '''' // sdc // ''''
  Write (io, Fmt='(A)') ' INTERRUPT CHARACTER ''' // int // &
    ''', STRING DELIMITER ' // word2 // ', COMMENT DELIMITER ' // word1
  Write (io, Fmt='(A)') ' KEYWORD CHARACTER SET ''' // kcs(1:kcl) // ''''

! .....close output file if permitted
  If (lcl) Close (Unit=io)


  Return
End Subroutine
Subroutine ffrint(ic)

! docs ffrint
! free-field read : define interrupt character
! docc    call ffrint ( ic )
! ic     : a character variable that defines in its first
! position the new interrupt character
! docr    latest revision : 26/09/86
! docp    programmer      : d.husfeld
! doce

! free-field read : control variables                        ffrctrl
! latest revision : 10/10/86                                 ffrctrl
! iout = output unit number                                  ffrctrl
! iver = verification level                                  ffrctrl
! nor = no run indicator                                     ffrctrl
! cdi = comment delimiting indicator                         ffrctrl
! sdi = string delimiting indicator                          ffrctrl
! cdc = comment delimiting character                         ffrctrl
! sdc = string delimiting character                          ffrctrl
! int = interrupt character                                  ffrctrl
! kcs = keyword character set                                ffrctrl
! kcl = keyword character set length                         ffrctrl
! free-field read : error status table                       ffrerrs
! latest revision : 26/09/86                                 ffrerrs
! errprg = name of routine that detected the error           ffrerrs
! errmsg = error message                                     ffrerrs

! free-field read : error status clearance                   ffreclr
! latest revision : 26/09/86                                 ffreclr
! .. parameters ..
  Use :: nlte_type
  Implicit None

  Integer (i4b) :: kclm
  Parameter (kclm=64)
! ..
! .. scalar arguments ..
  Character :: ic*(*)
! ..
! .. scalars in common ..
  Integer (i4b) :: iout, iver, kcl
  Logical :: cdi, nor, sdi
  Character :: cdc*1, int*1, sdc*1, errprg*6, errmsg*40, kcs*(kclm)
! ..
! .. local scalars ..
  Integer (i4b) :: ierr
! ..
! .. common blocks ..
  Common /ffrct1/iout, iver, nor, cdi, sdi, kcl
  Common /ffrct2/cdc, sdc, int, kcs
  Common /ffrerr/errprg, errmsg
! ..
  errprg = ' '
  errmsg = ' '

  int = ic(1:1)
  If (iver>=1) Write (iout, Fmt='(A)', Iostat=ierr) &
    ' ** INTERRUPT CHARACTER IS SET TO ''' // int // ''''
  Return

End Subroutine
Subroutine ffripr

! docs ffripr
! free-field read : interrupt processing
! docc    *** for internal use only
! call ffripr
! docr    latest revision : 20/01/99
! docp    programmer      : d.husfeld/ chr. wiethaus
! doce

! free-field read : control variables                        ffrctrl
! latest revision : 10/10/86                                 ffrctrl
! iout = output unit number                                  ffrctrl
! iver = verification level                                  ffrctrl
! nor = no run indicator                                     ffrctrl
! cdi = comment delimiting indicator                         ffrctrl
! sdi = string delimiting indicator                          ffrctrl
! cdc = comment delimiting character                         ffrctrl
! sdc = string delimiting character                          ffrctrl
! int = interrupt character                                  ffrctrl
! kcs = keyword character set                                ffrctrl
! kcl = keyword character set length                         ffrctrl
! free-field read : error status table                       ffrerrs
! latest revision : 26/09/86                                 ffrerrs
! errprg = name of routine that detected the error           ffrerrs
! errmsg = error message                                     ffrerrs
! free-field read : file status variables                    ffrstat
! latest revision : 05/11/86                                 ffrstat
! inp = input unit number                                    ffrstat
! icol = column number in line                               ffrstat
! iseq = line sequence number                                ffrstat
! eof = end-of-file indicator                                ffrstat
! lin = line image                                           ffrstat

! free-field read : error status clearance                   ffreclr
! latest revision : 26/09/86                                 ffreclr
! .. parameters ..
  Use :: nlte_type
  Implicit None

  Integer (i4b) :: kclm
  Parameter (kclm=64)
! ..
! .. scalars in common ..
  Integer (i4b) :: icol, inp, iout, iseq, iver, kcl
  Logical :: cdi, eof, nor, sdi
  Character :: cdc*1, int*1, sdc*1, errprg*6, errmsg*40, lin*80, kcs*(kclm)
! ..
! .. local scalars ..
  Integer (i4b) :: i, ierr, j
  Character :: chr*1, del*1, equ*1
! ..
! .. external subroutines ..
  External :: ffrifc
! ..
! .. common blocks ..
  Common /ffrct1/iout, iver, nor, cdi, sdi, kcl
  Common /ffrct2/cdc, sdc, int, kcs
  Common /ffrerr/errprg, errmsg
  Common /ffrst1/inp, icol, eof, iseq
  Common /ffrst2/lin
! ..
! .. data statements ..
  Data equ/'='/
! ..
  errprg = ' '
  errmsg = ' '


! .....confirm call
  If (lin(1:1)/=int) Return

  chr = lin(2:2)
  del = lin(3:3)

! .....c  comment delimiter
  If (chr=='C') Then
    If (del==equ) Then
      cdi = .True.
      cdc = lin(4:4)

      If (iver>=1) Write (iout, Fmt='(A)', Iostat=ierr) &
        ' ** COMMENT DELIMITER SET TO ''' // cdc // ''''
    Else
      cdi = .False.
      cdc = ' '
      If (iver>=1) Write (iout, Fmt='(A)', Iostat=ierr) &
        ' ** COMMENT DELIMITER CLEARED'
    End If

!   .....d  string delimiter
  Else If (chr=='D') Then
    If (del==equ) Then
      sdi = .True.
      sdc = lin(4:4)

      If (iver>=1) Write (iout, Fmt='(A)', Iostat=ierr) &
        ' ** STRING DELIMITER IS SET TO ''' // sdc // ''''
    Else
      sdi = .False.
      sdc = ' '
      If (iver>=1) Write (iout, Fmt='(A)', Iostat=ierr) &
        ' ** STRING DELIMITER CLEARED'
    End If

!   .....i  interrupt character
  Else If (chr=='I') Then
    int = lin(3:3)
    If (iver>=1) Write (iout, Fmt='(A)', Iostat=ierr) &
      ' ** INTERRUPT CHARACTER IS SET TO ''' // int // ''''

!   .....k  keyword character set
  Else If (chr=='K') Then
    If (del==equ) Then
      kcs(1:1) = lin(4:4)
      kcl = 1
      Do i = 2, kclm
        j = i + 3
        If (lin(j:j)==kcs(kcl:kcl)) Go To 100
        kcs(i:i) = lin(j:j)
        kcl = i
      End Do

100   Continue
    Else
      kcs(1:26) = 'ABCDEFGHIJKLMNOPQRSTUVWXYZ'
      kcs(27:62) = 'abcdefghijklmnopqrstuvwxyz0123456789'
      kcl = 62
    End If
    If (iver>=1) Write (iout, Fmt='(A)', Iostat=ierr) &
      ' ** NEW KEYWORD CHARACTER SET IS :', '    ''' // kcs(1:kcl) // ''''

!   .....s  send line to interface routine 'ffrifc'
  Else If (chr=='S') Then
    Call ffrifc(lin)
    If (iver>=1) Write (iout, Fmt='(A)', Iostat=ierr) &
      ' ** FOLLOWING LINE SENT TO INTERFACE :', '''' // lin // ''''

!   .....v  verification level
  Else If (chr=='V') Then
    Read (lin(3:3), Fmt='(I1)', Iostat=ierr) i
    If (ierr==0) iver = i
    If (iout==0) iver = 0

!   .....any other character : comment line
  Else
    If (iver>=1) Write (iout, Fmt='(A)', Iostat=ierr) ' *' // lin(3:80)
  End If

  Return

End Subroutine
Subroutine ffrkcs(cs)

! docs ffrkcs
! free-field read : defined keyword character set
! docc    call ffrkcs ( cs )
! cs     : a character variable that defines the character set
! acceptable in any keyword following.
! if the length of cs exceeds the maximum length
! allowed, cs will be truncated.
! sl     : a statement label in the calling program to which
! control is returned from ffrkcs in case of an error
! docr    latest revision : 26/09/86
! docp    programmer      : d.husfeld
! doce

! free-field read : control variables                        ffrctrl
! latest revision : 10/10/86                                 ffrctrl
! iout = output unit number                                  ffrctrl
! iver = verification level                                  ffrctrl
! nor = no run indicator                                     ffrctrl
! cdi = comment delimiting indicator                         ffrctrl
! sdi = string delimiting indicator                          ffrctrl
! cdc = comment delimiting character                         ffrctrl
! sdc = string delimiting character                          ffrctrl
! int = interrupt character                                  ffrctrl
! kcs = keyword character set                                ffrctrl
! kcl = keyword character set length                         ffrctrl
! free-field read : error status table                       ffrerrs
! latest revision : 26/09/86                                 ffrerrs
! errprg = name of routine that detected the error           ffrerrs
! errmsg = error message                                     ffrerrs

! free-field read : error status clearance                   ffreclr
! latest revision : 26/09/86                                 ffreclr
! .. parameters ..
  Use :: nlte_type
  Implicit None

  Integer (i4b) :: kclm
  Parameter (kclm=64)
! ..
! .. scalar arguments ..
  Character :: cs*(*)
! ..
! .. scalars in common ..
  Integer (i4b) :: iout, iver, kcl
  Logical :: cdi, nor, sdi
  Character :: cdc*1, int*1, sdc*1, errprg*6, errmsg*40, kcs*(kclm)
! ..
! .. local scalars ..
  Integer (i4b) :: ierr, l
! ..
! .. intrinsic functions ..
! ..
! .. common blocks ..
  Common /ffrct1/iout, iver, nor, cdi, sdi, kcl
  Common /ffrct2/cdc, sdc, int, kcs
  Common /ffrerr/errprg, errmsg
! ..
  errprg = ' '
  errmsg = ' '

  l = len(cs)
  If (l>kclm) Then
    If (iver>=1) Write (iout, Fmt='(A,I3,A)', Iostat=ierr) &
      ' ** NEW KEYWORD CHARACTER SET TRUNCATED TO ', kclm, ' CHARACTERS'
    kcl = kclm

    kcs = cs(1:kclm)
  Else
    kcl = l
    kcs(1:kcl) = cs
  End If

  If (iver>=1) Write (iout, Fmt='(A)', Iostat=ierr) &
    ' ** NEW KEYWORD CHARACTER SET IS :', '    ''' // kcs(1:kcl) // ''''
  Return

End Subroutine
Subroutine ffrkey(key, nc, retchar)

! docs ffrkey
! free-field read : get next keyword
! docc    call ffrkey ( key,nc, retchar )
! key    : a character variable that receives the next keyword
! found in the input stream
! nc     : an integer variable that returns the number
! of characters stored into key
! ret1    : a statement label in the calling program to which
! control is returned from ffrkey when an eof is
! found during search for a keyword.
! note that the keyword will be successfully read
! even when it is terminated by an eof.
! ret2    : a statement label in the calling program to which
! control is returned from ffrkey in case of an error
! docr    latest revision : 03/03/87
! docp    programmer      : d.husfeld
! doce

! free-field read : control variables                        ffrctrl
! latest revision : 10/10/86                                 ffrctrl
! iout = output unit number                                  ffrctrl
! iver = verification level                                  ffrctrl
! nor = no run indicator                                     ffrctrl
! cdi = comment delimiting indicator                         ffrctrl
! sdi = string delimiting indicator                          ffrctrl
! cdc = comment delimiting character                         ffrctrl
! sdc = string delimiting character                          ffrctrl
! int = interrupt character                                  ffrctrl
! kcs = keyword character set                                ffrctrl
! kcl = keyword character set length                         ffrctrl
! free-field read : error status table                       ffrerrs
! latest revision : 26/09/86                                 ffrerrs
! errprg = name of routine that detected the error           ffrerrs
! errmsg = error message                                     ffrerrs

! free-field read : error status clearance                   ffreclr
! latest revision : 26/09/86                                 ffreclr
! .. parameters ..
  Use :: nlte_type
  Implicit None

  Integer (i4b) :: kclm
  Parameter (kclm=64)
! ..
! .. scalar arguments ..
  Integer (i4b) :: nc
  Character :: key*(*)
! ..
! .. scalars in common ..
  Integer (i4b) :: iout, iver, kcl
  Logical :: cdi, nor, sdi
  Character :: cdc*1, int*1, sdc*1, errprg*6, errmsg*40, kcs*(kclm)
! ..
! .. local scalars ..
  Integer (i4b) :: i, ierr, l
  Character :: ch*1, retchar*4, ret*4
! ..
! .. external subroutines ..
  External :: ffrnch
! ..
! .. intrinsic functions ..
! ..
! .. common blocks ..
  Common /ffrct1/iout, iver, nor, cdi, sdi, kcl
  Common /ffrct2/cdc, sdc, int, kcs
  Common /ffrerr/errprg, errmsg
! ..
  errprg = ' '
  errmsg = ' '

  l = len(key)
  key = ' '
  If (iver>=3) Write (iout, Fmt='(A)', Iostat=ierr) ' **** READ NEXT KEYWORD'

! .....scan input for first acceptable character
100 Continue
  Call ffrnch(ch, ret)

  Select Case (ret)
  Case ('RET1')
    Go To 120
  Case ('RET2')
    Go To 130
  Case ('RET0')
    Continue
  Case Default
    Stop ' WRONG RETURN CONDITION'
  End Select

  If (ch==' ') Go To 100
  If (index(kcs(1:kcl),ch)==0) Go To 100

! .....first acceptable character found.
! now icollect characters in key until either key is full
! or next non-acceptable character occurs.
  key(1:1) = ch
  nc = 1
  Do i = 2, l
    Call ffrnch(ch, ret)

    Select Case (ret)
    Case ('RET1')
      Go To 110
    Case ('RET2')
      Go To 130
    Case ('RET0')
      Continue
    Case Default
      Stop ' WRONG RETURN CONDITION'
    End Select

    If (index(kcs(1:kcl),ch)==0) Go To 110
    key(i:i) = ch
    nc = i
  End Do
  Call ffrnch(ch, ret)

  Select Case (ret)
  Case ('RET1')
    Go To 110
  Case ('RET2')
    Go To 130
  Case ('RET0')
    Continue
  Case Default
    Stop ' WRONG RETURN CONDITION'
  End Select

110 Continue

  If (iver>=2) Write (iout, Fmt='(A,A)', Iostat=ierr) ' *** KEYWORD FOUND : ', &
    key(1:nc)

  retchar = 'RET0'
  Return

! .....eof condition exit
120 Continue
  retchar = 'RET1'
  Return

! .....error exit
130 Continue
  retchar = 'RET2'
  Return

End Subroutine
Character (*) Function ffrlin()

! docs ffrlin
! free-field read : supply current input line
! docc    line = ffrlin()
! line   : a character variable, character variable substring,
! character array element or character array element
! substring that receives a copy of the current input
! line.
! docr    latest revision : 02/03/87
! docp     programmmer    : d. husfeld
! doce

! free-field read : file status variables                    ffrstat
! latest revision : 05/11/86                                 ffrstat
! inp = input unit number                                    ffrstat
! icol = column number in line                               ffrstat
! iseq = line sequence number                                ffrstat
! eof = end-of-file indicator                                ffrstat
! lin = line image                                           ffrstat
! .. scalars in common ..
  Use :: nlte_type
  Implicit None

  Integer (i4b) :: icol, inp, iseq
  Logical :: eof
  Character :: lin*80
! ..
! .. common blocks ..
  Common /ffrst1/inp, icol, eof, iseq
  Common /ffrst2/lin
! ..
  ffrlin = lin

  Return
End Function
Subroutine ffrloc(str, retchar)

! docs ffrloc
! free-field read : locate specified string and position read
! pointer behind it.
! docc    call ffrloc ( str,retchar )
! str    : a character expression that specifies the string
! to be located.
! ret1   : a statement label in the calling program to which
! control is returned from ffrloc when an eof is found
! ret2   : a statement label in the calling program to which
! control is returned from ffrloc in case of an error
! docr    latest revision : 05/11/86
! docp    programmer      : d.husfeld
! doce

! free-field read : control variables                        ffrctrl
! latest revision : 10/10/86                                 ffrctrl
! iout = output unit number                                  ffrctrl
! iver = verification level                                  ffrctrl
! nor = no run indicator                                     ffrctrl
! cdi = comment delimiting indicator                         ffrctrl
! sdi = string delimiting indicator                          ffrctrl
! cdc = comment delimiting character                         ffrctrl
! sdc = string delimiting character                          ffrctrl
! int = interrupt character                                  ffrctrl
! kcs = keyword character set                                ffrctrl
! kcl = keyword character set length                         ffrctrl
! free-field read : error status table                       ffrerrs
! latest revision : 26/09/86                                 ffrerrs
! errprg = name of routine that detected the error           ffrerrs
! errmsg = error message                                     ffrerrs
! free-field read : file status variables                    ffrstat
! latest revision : 05/11/86                                 ffrstat
! inp = input unit number                                    ffrstat
! icol = column number in line                               ffrstat
! iseq = line sequence number                                ffrstat
! eof = end-of-file indicator                                ffrstat
! lin = line image                                           ffrstat

! free-field read : error status clearance                   ffreclr
! latest revision : 26/09/86                                 ffreclr
! .. parameters ..
  Use :: nlte_type
  Implicit None

  Integer (i4b) :: kclm
  Parameter (kclm=64)
! ..
! .. scalar arguments ..
  Character :: str*(*)
! ..
! .. scalars in common ..
  Integer (i4b) :: icol, inp, iout, iseq, iver, kcl
  Logical :: cdi, eof, nor, sdi
  Character :: cdc*1, int*1, sdc*1, errprg*6, errmsg*40, lin*80, kcs*(kclm)
! ..
! .. local scalars ..
  Integer (i4b) :: i, ierr, l
  Character :: ch1*1, ch2*1, retchar*4, ret*4
! ..
! .. external subroutines ..
  External :: ffrnch
! ..
! .. intrinsic functions ..
! ..
! .. common blocks ..
  Common /ffrct1/iout, iver, nor, cdi, sdi, kcl
  Common /ffrct2/cdc, sdc, int, kcs
  Common /ffrerr/errprg, errmsg
  Common /ffrst1/inp, icol, eof, iseq
  Common /ffrst2/lin
! ..
  errprg = ' '
  errmsg = ' '

  If (iver>=3) Write (iout, Fmt='(A,A)', Iostat=ierr) &
    ' **** LOCATING STRING : ', str

  l = len(str)
  ch2 = str(1:1)
100 Continue
  Call ffrnch(ch1, ret)

  Select Case (ret)
  Case ('RET1')
    Go To 120
  Case ('RET2')
    Go To 130
  Case ('RET0')
    Continue
  Case Default
    Stop ' WRONG RETURN CONDITION'
  End Select

110 Continue
  If (ch1/=ch2) Go To 100
  Do i = 2, l
    Call ffrnch(ch1, ret)

    Select Case (ret)
    Case ('RET1')
      Go To 120
    Case ('RET2')
      Go To 130
    Case ('RET0')
      Continue
    Case Default
      Stop ' WRONG RETURN CONDITION'
    End Select

    If (ch1/=str(i:i)) Go To 110
  End Do

! .....verify
  If (iver>=2) Write (iout, Fmt='(A,I5)', Iostat=ierr) &
    ' *** SPECIFIED STRING LOCATED IN LINE #', iseq

  retchar = 'RET0'
  Return

! .....eof condition exit
120 Continue
  errprg = 'FFRLOC'
  errmsg = 'EOF WHILE LOCATING ''' // str // ''''
  retchar = 'RET1'
  Return

! .....error exit
130 Continue
  retchar = 'RET2'
  Return

End Subroutine
Subroutine ffrncd(retchar)

! docs ffrncd
! free-field read : load next card into memory
! docc    call ffrncd ( retchar )
! ret1   : a statement label in the calling program to which
! control is returned from ffrncd when an eof is found
! ret2   : a statement label in the calling program to which
! control is returned from ffrncd in case of an error
! docr    latest revision : 05/11/86
! docp    programmer      : d.husfeld
! doce

! free-field read : control variables                        ffrctrl
! latest revision : 10/10/86                                 ffrctrl
! iout = output unit number                                  ffrctrl
! iver = verification level                                  ffrctrl
! nor = no run indicator                                     ffrctrl
! cdi = comment delimiting indicator                         ffrctrl
! sdi = string delimiting indicator                          ffrctrl
! cdc = comment delimiting character                         ffrctrl
! sdc = string delimiting character                          ffrctrl
! int = interrupt character                                  ffrctrl
! kcs = keyword character set                                ffrctrl
! kcl = keyword character set length                         ffrctrl
! free-field read : error status table                       ffrerrs
! latest revision : 26/09/86                                 ffrerrs
! errprg = name of routine that detected the error           ffrerrs
! errmsg = error message                                     ffrerrs
! free-field read : file status variables                    ffrstat
! latest revision : 05/11/86                                 ffrstat
! inp = input unit number                                    ffrstat
! icol = column number in line                               ffrstat
! iseq = line sequence number                                ffrstat
! eof = end-of-file indicator                                ffrstat
! lin = line image                                           ffrstat

! free-field read : error status clearance                   ffreclr
! latest revision : 26/09/86                                 ffreclr
! .. parameters ..
  Use :: nlte_type
  Implicit None

  Character (6) :: prgnam
  Parameter (prgnam='FFRNCD')
  Integer (i4b) :: kclm
  Parameter (kclm=64)
! ..
! .. scalars in common ..
  Integer (i4b) :: icol, inp, iout, iseq, iver, kcl
  Logical :: cdi, eof, nor, sdi
  Character :: cdc*1, int*1, sdc*1, errprg*6, errmsg*40, lin*80, kcs*(kclm)

  Character (4) :: retchar
! ..
! .. external subroutines ..
  External :: ffripr
! ..
! .. common blocks ..
  Common /ffrct1/iout, iver, nor, cdi, sdi, kcl
  Common /ffrct2/cdc, sdc, int, kcs
  Common /ffrerr/errprg, errmsg
  Common /ffrst1/inp, icol, eof, iseq
  Common /ffrst2/lin
! ..
  errprg = ' '
  errmsg = ' '

! .....return immediately if eof flag is already set
  If (eof) Then
    retchar = 'RET1'
    Return
  End If

100 Continue
  Read (inp, Fmt='(A)', End=110, Err=120) lin
  icol = 1
  iseq = iseq + 1
  eof = .False.

! .....check for interrupt
  If (lin(1:1)==int) Then
    Call ffripr
    Go To 100
  End If

! .....normal exit
  retchar = 'RET0'
  Return

! .....eof condition exit
110 Continue
  errprg = prgnam
  Write (errmsg, Fmt='(A,I3)') 'EOF ENCOUNTERED ON UNIT #', inp
  icol = 1
  eof = .True.
  lin = ' '
  retchar = 'RET1'
  Return

! .....read error exit
120 Continue
  errprg = prgnam
  Write (errmsg, Fmt='(A,I3)') 'READ ERROR ON UNIT #', inp
  Go To 130

! .....error exit
130 Continue
  retchar = 'RET2'
  Return

End Subroutine
Subroutine ffrnch(ch, retchar)

! docs ffrnch
! free-field read : get next character
! docc    *** for internal use only
! call ffrnch ( ch,retchar )
! ch     : a character variable that receives in its first po-
! sition the next character found in the input stream
! ret1   : a statement label in the calling program to which
! control is returned from ffrnch when an eof is found
! ret2   : a statement label in the calling program to which
! control is returned from ffrnch in case of an error
! docr    latest revision : 03/03/87
! docp    programmer      : d.husfeld
! doce

! free-field read : control variables                        ffrctrl
! latest revision : 10/10/86                                 ffrctrl
! iout = output unit number                                  ffrctrl
! iver = verification level                                  ffrctrl
! nor = no run indicator                                     ffrctrl
! cdi = comment delimiting indicator                         ffrctrl
! sdi = string delimiting indicator                          ffrctrl
! cdc = comment delimiting character                         ffrctrl
! sdc = string delimiting character                          ffrctrl
! int = interrupt character                                  ffrctrl
! kcs = keyword character set                                ffrctrl
! kcl = keyword character set length                         ffrctrl
! free-field read : error status table                       ffrerrs
! latest revision : 26/09/86                                 ffrerrs
! errprg = name of routine that detected the error           ffrerrs
! errmsg = error message                                     ffrerrs
! free-field read : file status variables                    ffrstat
! latest revision : 05/11/86                                 ffrstat
! inp = input unit number                                    ffrstat
! icol = column number in line                               ffrstat
! iseq = line sequence number                                ffrstat
! eof = end-of-file indicator                                ffrstat
! lin = line image                                           ffrstat

! free-field read : error status clearance                   ffreclr
! latest revision : 26/09/86                                 ffreclr
! .. parameters ..
  Use :: nlte_type
  Implicit None

  Integer (i4b) :: kclm
  Parameter (kclm=64)
! ..
! .. scalar arguments ..
  Character :: ch*(*)
! ..
! .. scalars in common ..
  Integer (i4b) :: icol, inp, iout, iseq, iver, kcl
  Logical :: cdi, eof, nor, sdi
  Character :: cdc*1, int*1, sdc*1, errprg*6, errmsg*40, lin*80, kcs*(kclm)
! ..
! .. local scalars ..
  Integer (i4b) :: i, iicol
  Character (4) :: retchar, ret
! ..
! .. external subroutines ..
  External :: ffrncd
! ..
! .. common blocks ..
  Common /ffrct1/iout, iver, nor, cdi, sdi, kcl
  Common /ffrct2/cdc, sdc, int, kcs
  Common /ffrerr/errprg, errmsg
  Common /ffrst1/inp, icol, eof, iseq
  Common /ffrst2/lin
! ..
  errprg = ' '
  errmsg = ' '

  If (icol>80) Then
    Call ffrncd(ret)
    Select Case (ret)
    Case ('RET1')
      Go To 110
    Case ('RET2')
      Go To 120
    Case ('RET0')
      Continue
    Case Default
      Stop ' WRONG RETURN CONDITION'
    End Select
  End If

  ch(1:1) = lin(icol:icol)
  icol = icol + 1
  retchar = 'RET0'
  If (.Not. cdi) Return
  If (ch(1:1)/=cdc) Return
! ...comment delimiter found. skip to next comment delimiter or end
! of line and return a blank in ch.
  ch(1:1) = ' '
  iicol = icol
  Do i = iicol, 80
    icol = i
    If (lin(i:i)==cdc) Go To 100
  End Do
100 Continue
  icol = icol + 1
  Return

! .....eof condition exit
110 Continue
  retchar = 'RET1'
  Return

! .....error exit
120 Continue
  retchar = 'RET2'
  Return

End Subroutine
Subroutine ffrnum(realvar, integ, retchar)

! docs ffrnum
! free-field read : get next number
! docc    call ffrnum ( realvar,integ,retchar )
! realvar: is a real variable that receives the number read
! integ  : an integer variable that receives the number read
! ret1   : a statement label in the calling program to which
! control is returned from ffrnum when an eof is
! found during search for a valid number.
! note that the number will be correctly read even
! when it is terminated by an eof.
! ret2   : a statement label in the calling program to which
! control is returned from ffrnum in case of an error.
! docr    latest revision : 03/03/87
! docp    programmer      : d.husfeld
! doce

! free-field read : control variables                        ffrctrl
! latest revision : 10/10/86                                 ffrctrl
! iout = output unit number                                  ffrctrl
! iver = verification level                                  ffrctrl
! nor = no run indicator                                     ffrctrl
! cdi = comment delimiting indicator                         ffrctrl
! sdi = string delimiting indicator                          ffrctrl
! cdc = comment delimiting character                         ffrctrl
! sdc = string delimiting character                          ffrctrl
! int = interrupt character                                  ffrctrl
! kcs = keyword character set                                ffrctrl
! kcl = keyword character set length                         ffrctrl
! free-field read : error status table                       ffrerrs
! latest revision : 26/09/86                                 ffrerrs
! errprg = name of routine that detected the error           ffrerrs
! errmsg = error message                                     ffrerrs
! .. parameters ..
  Use :: nlte_type
  Implicit None

  Character (6) :: prgnam
  Parameter (prgnam='FFRNUM')
  Integer (i4b) :: ml
  Parameter (ml=50)
  Integer (i4b) :: kclm
  Parameter (kclm=64)
! ..
! .. scalar arguments ..
  Real (dp) :: realvar
  Integer (i4b) :: integ
! ..
! .. scalars in common ..
  Integer (i4b) :: iout, iver, kcl
  Logical :: cdi, nor, sdi
  Character :: cdc*1, int*1, sdc*1, errprg*6, errmsg*40, kcs*(kclm)
! ..
! .. local scalars ..
  Real (dp) :: intmax
  Integer (i4b) :: i, ierr, jump, l
  Logical :: digit
  Character :: ch*1, form1*10, set*15, num*(ml), retchar*4, ret*4
! ..
! .. external subroutines ..
  External :: ffrnch
! ..
! .. intrinsic functions ..
! ..
! .. common blocks ..
  Common /ffrct1/iout, iver, nor, cdi, sdi, kcl
  Common /ffrct2/cdc, sdc, int, kcs
  Common /ffrerr/errprg, errmsg
! ..
! .. data statements ..
  Data set/'.0123456789-+ED'/, intmax/32767/
! ..
! +                  ,intmax/281474976710655/

! free-field read : error status clearance                   ffreclr
! latest revision : 26/09/86                                 ffreclr
  errprg = ' '
  errmsg = ' '

  If (iver>=3) Write (iout, Fmt='(A)', Iostat=ierr) ' **** READ NEXT NUMBER'

! .....clear
100 Continue
  l = 0
  num = ' '
  digit = .False.

! .....skip over non-acceptable characters
110 Continue
  Call ffrnch(ch, ret)

  Select Case (ret)
  Case ('RET1')
    Go To 220
  Case ('RET2')
    Go To 240
  Case ('RET0')
    Continue
  Case Default
    Stop ' WRONG RETURN CONDITION'
  End Select

  If (ch==' ') Go To 110
  i = index(set, ch)
  If (i==0) Go To 110

! .....sign field
120 Continue
  If (i==12) Then
    jump = 2100
    Go To 200
  End If

! .....integer field
130 Continue
  jump = 2100
  If ((i>=2) .And. (i<=11)) Then
    digit = .True.

    Go To 200
  End If

! .....decimal point
140 Continue
  If (i==1) Then
    jump = 2300
    Go To 200
  End If

! .....fractional field
150 Continue
  jump = 2300
  If ((i>=2) .And. (i<=11)) Then
    digit = .True.
    Go To 200
  End If

! .....restart read if no digit was found yet
  If (.Not. digit) Go To 100

! .....exponent indicator field
160 Continue
  If ((i==14) .Or. (i==15)) Then
    digit = .False.
    jump = 2500
    Go To 200
  End If

! .....exponent sign field
170 Continue
  If ((i==12) .Or. (i==13)) Then
    digit = .False.
    jump = 2600
    Go To 200
  End If

! .....exponent integer field
180 Continue
  jump = 2600
  If ((i>=2) .And. (i<=11)) Then
    digit = .True.
    Go To 200
  End If

! .....store '0' if no exponent integer was read yet
  If (.Not. digit) Then
    If (l>=ml) Go To 210
    l = l + 1
    num(l:l) = '0'
  End If

! .....read number into real variable

190 Continue
  Write (form1, Fmt='(A,I2,A)') '(E', l, '.0)'

! .....high art of fortran: if input has decimal point, location of decimal
! point overrides format specification: thus, only length of field has
! to be correct

  Read (num, Fmt=form1, Err=230) realvar

! .....convert into integer variable
  integ = nint(sign(min(abs(realvar),intmax),realvar))

! .....verify
  If (iver>=2) Write (iout, Fmt='(A,1PG21.14,A,I6)') ' *** NUMBER FOUND :', &
    realvar, ' <-> ', integ

! .....normal exit
  retchar = 'RET0'
  Return

! .....procedure to store current character and to fetch new one
200 Continue
  If (l>=ml) Go To 210
  l = l + 1
  num(l:l) = ch
  Call ffrnch(ch, ret)

  Select Case (ret)
  Case ('RET1')
    Go To 190
  Case ('RET2')
    Go To 240
  Case ('RET0')
    Continue
  Case Default
    Stop ' WRONG RETURN CONDITION'
  End Select

  i = index(set, ch)
  Select Case (jump)
  Case (2000)
    Go To 120
  Case (2100)
    Go To 130
  Case (2200)
    Go To 140
  Case (2300)
    Go To 150
  Case (2400)
    Go To 160
  Case (2500)
    Go To 170
  Case (2600)
    Go To 180
  Case Default
    Stop 'ERROR IN JUMP'
  End Select

! .....error exit when read buffer 'num' is full
210 Continue
  errprg = prgnam
  errmsg = 'BUFFER FULL DURING READ OF NUMBER'
  Go To 240

! .....eof condition exit
220 Continue
  retchar = 'RET1'
  Return

! .....read error exit
230 Continue
  errprg = prgnam
  errmsg = 'ERROR DURING READ OF NUMBER'
  Go To 240

! .....error exit
240 Continue
  retchar = 'RET2'
  Return

End Subroutine
Subroutine ffrout(iun, retchar)

! docs ffrout
! free-field read : define output unit number
! docc    call ffrout ( iun,retchar )
! iun    : an integer variable that identifies the unit number
! to be used for any further output of the free-field
! read routines if its value is in the range from 1
! to 999. the referenced output unit must already be
! be opened for sequential and formatted access.
! if the value of iun is zero when ffrout is called
! un receives the current output unit number and no
! further action is taken.
! ret1   : a statement label in the calling program to which
! control is returned from ffrout in case of an error
! docr    latest revision : 10/10/86
! docp    programmer      : d.husfeld
! doce

! free-field read : control variables                        ffrctrl
! latest revision : 10/10/86                                 ffrctrl
! iout = output unit number                                  ffrctrl
! iver = verification level                                  ffrctrl
! nor = no run indicator                                     ffrctrl
! cdi = comment delimiting indicator                         ffrctrl
! sdi = string delimiting indicator                          ffrctrl
! cdc = comment delimiting character                         ffrctrl
! sdc = string delimiting character                          ffrctrl
! int = interrupt character                                  ffrctrl
! kcs = keyword character set                                ffrctrl
! kcl = keyword character set length                         ffrctrl
! free-field read : error status table
! latest revision : 26/09/86                                 ffrerrs
! errprg = name of routine that detected the error           ffrerrs
! errmsg = error message                                     ffrerrs

! free-field read : error status clearance                   ffreclr
! latest revision : 26/09/86                                 ffreclr
! .. parameters ..
  Use :: nlte_type
  Implicit None

  Integer (i4b) :: kclm
  Parameter (kclm=64)
! ..
! .. scalar arguments ..
  Integer (i4b) :: iun
! ..
! .. scalars in common ..
  Integer (i4b) :: iout, iver, kcl
  Logical :: cdi, nor, sdi
  Character :: cdc*1, int*1, sdc*1, errprg*6, errmsg*40, kcs*(kclm)
! ..
! .. local scalars ..
  Integer (i4b) :: ierr
  Logical :: uop
  Character :: acmode*10, fmmode*10, retchar*4
! ..
! .. common blocks ..
  Common /ffrct1/iout, iver, nor, cdi, sdi, kcl
  Common /ffrct2/cdc, sdc, int, kcs
  Common /ffrerr/errprg, errmsg
! ..
  errprg = ' '
  errmsg = ' '

! .....return current output unit number if iun is zero
  retchar = 'RET0'
  If (iun==0) Then
    iun = iout
    Return
  End If

! .....check iun
  If ((iun<1) .Or. (iun>999)) Then
    If (iver>=1) Write (iout, Fmt='(A,I10)', Iostat=ierr) &
      ' *** NON-ACCEPTABLE UNIT NUMBER ', iun

    Go To 100

  End If
  Inquire (Unit=iun, Iostat=ierr, Opened=uop, Access=acmode, Form=fmmode)
  If ((.Not. uop) .Or. (ierr/=0) .Or. (acmode/='SEQUENTIAL') .Or. &
    (fmmode/='FORMATTED')) Then
    If (iver>=1) Write (iout, Fmt='(A,I3,A)', Iostat=ierr) ' *** UNIT #', iun, &
      ' NOT OPENED PROPERLY'

    Go To 100
  End If

! .....accept unit number
  If (iver>=1) Write (iout, Fmt='(A,I3)', Iostat=ierr) &
    ' ** SWITCHING TO OUTPUT UNIT #', iun
  iout = iun
  Return

! .....error exit
100 Continue
  retchar = 'RET1'
  Return

End Subroutine
Subroutine ffrrtn(iun)

! docs ffrrtn
! free-field read : return specified input unit
! docc    call ffrrtn ( iun )
! iun     : an integer variable which value defines the unit
! to be closed for free-field read
! docr    latest revision : 26/09/86
! docp    programmer      : d.husfeld
! doce

! free-field read : control variables                        ffrctrl
! latest revision : 10/10/86                                 ffrctrl
! iout = output unit number                                  ffrctrl
! iver = verification level                                  ffrctrl
! nor = no run indicator                                     ffrctrl
! cdi = comment delimiting indicator                         ffrctrl
! sdi = string delimiting indicator                          ffrctrl
! cdc = comment delimiting character                         ffrctrl
! sdc = string delimiting character                          ffrctrl
! int = interrupt character                                  ffrctrl
! kcs = keyword character set                                ffrctrl
! kcl = keyword character set length                         ffrctrl
! free-field read : error status table                       ffrerrs
! latest revision : 26/09/86                                 ffrerrs
! errprg = name of routine that detected the error           ffrerrs
! errmsg = error message                                     ffrerrs
! free-field read : file status variables                    ffrstat
! latest revision : 05/11/86                                 ffrstat
! inp = input unit number                                    ffrstat
! icol = column number in line                               ffrstat
! iseq = line sequence number                                ffrstat
! eof = end-of-file indicator                                ffrstat
! lin = line image                                           ffrstat
! free-field read : unit table                               ffrutab
! latest revision : 05/11/86                                 ffrutab
! nun = number of units in table                             ffrutab
! inpt= input unit number table                              ffrutab
! icolt= column number table                                 ffrutab
! iseqt= line sequence number table                          ffrutab
! lint= line image table                                     ffrutab

! free-field read : error status clearance                   ffreclr
! latest revision : 26/09/86                                 ffreclr
! .. parameters ..
  Use :: nlte_type
  Implicit None

  Integer (i4b) :: kclm
  Parameter (kclm=64)
  Integer (i4b) :: mun
  Parameter (mun=10)
! ..
! .. scalar arguments ..
  Integer (i4b) :: iun
! ..
! .. scalars in common ..
  Integer (i4b) :: icol, inp, iout, iseq, iver, kcl, nun
  Logical :: cdi, eof, nor, sdi
  Character :: cdc*1, int*1, sdc*1, errprg*6, errmsg*40, lin*80, kcs*(kclm)
! ..
! .. arrays in common ..
  Integer (i4b) :: icolt(mun), inpt(mun), iseqt(mun)
  Character :: lint(mun)*80
! ..
! .. local scalars ..
  Integer (i4b) :: i, ierr, j
! ..
! .. common blocks ..
  Common /ffrct1/iout, iver, nor, cdi, sdi, kcl
  Common /ffrct2/cdc, sdc, int, kcs
  Common /ffrerr/errprg, errmsg
  Common /ffrst1/inp, icol, eof, iseq
  Common /ffrst2/lin
  Common /ffrut1/nun, inpt, icolt, iseqt
  Common /ffrut2/lint
! ..
  errprg = ' '
  errmsg = ' '

! .....search unit number in table
  Do i = 1, nun
    If (inpt(i)==iun) Then
!     ...shift all following entries one step down
      Do j = i + 1, nun
        inpt(j-1) = inpt(j)
        icolt(j-1) = icolt(j)
        lint(j-1) = lint(j)
      End Do
      nun = nun - 1

      Go To 100

    End If
  End Do
100 Continue
  If (iver>=1) Write (iout, Fmt='(A,I3,A)', Iostat=ierr) ' ** UNIT #', iun, &
    ' DISCARDED'
  Return

End Subroutine
Subroutine ffrsdc(dc)

! docs ffrsdc
! free-field read : define string delimiting character
! docc    call ffrsdc ( dc )
! dc     : a character variable that defines in its first
! position the new string delimiting character.
! however, if dc contains in its first five positions
! the string 'clear', the delimiting character
! becomes undefined.
! docr    latest revision : 26/09/86
! docp    programmer      : d.husfeld
! doce

! free-field read : control variables                        ffrctrl
! latest revision : 10/10/86                                 ffrctrl
! iout = output unit number                                  ffrctrl
! iver = verification level                                  ffrctrl
! nor = no run indicator                                     ffrctrl
! cdi = comment delimiting indicator                         ffrctrl
! sdi = string delimiting indicator                          ffrctrl
! cdc = comment delimiting character                         ffrctrl
! sdc = string delimiting character                          ffrctrl
! int = interrupt character                                  ffrctrl
! kcs = keyword character set                                ffrctrl
! kcl = keyword character set length                         ffrctrl
! free-field read : error status table                       ffrerrs
! latest revision : 26/09/86                                 ffrerrs
! errprg = name of routine that detected the error           ffrerrs
! errmsg = error message                                     ffrerrs

! free-field read : error status clearance                   ffreclr
! latest revision : 26/09/86                                 ffreclr
! .. parameters ..
  Use :: nlte_type
  Implicit None

  Integer (i4b) :: kclm
  Parameter (kclm=64)
! ..
! .. scalar arguments ..
  Character :: dc*(*)
! ..
! .. scalars in common ..
  Integer (i4b) :: iout, iver, kcl
  Logical :: cdi, nor, sdi
  Character :: cdc*1, int*1, sdc*1, errprg*6, errmsg*40, kcs*(kclm)
! ..
! .. local scalars ..
  Integer (i4b) :: ierr
! ..
! .. intrinsic functions ..
! ..
! .. common blocks ..
  Common /ffrct1/iout, iver, nor, cdi, sdi, kcl
  Common /ffrct2/cdc, sdc, int, kcs
  Common /ffrerr/errprg, errmsg
! ..
  errprg = ' '
  errmsg = ' '

  If ((len(dc)>=5) .And. (dc(1:5)=='CLEAR')) Then
    sdi = .False.
    sdc = ' '

    If (iver>=1) Write (iout, Fmt='(A)', Iostat=ierr) &
      ' ** STRING DELIMITER CLEARED'
  Else
    sdi = .True.
    sdc = dc(1:1)
    If (iver>=1) Write (iout, Fmt='(A)', Iostat=ierr) &
      ' ** STRING DELIMITER IS SET TO ''' // sdc // ''''

  End If
  Return

End Subroutine
Subroutine ffrstr(str, nc, retchar)

! docs ffrstr
! free-field read : get string
! docc    call ffrstr ( str,nc,retchar )
! str    : a character variable that receives the string
! nc     : an integer variable that returns the number of
! characters stored into str
! ret1   : a statement label in the calling program to which
! control is returned from ffrstr when an eof is found
! ret2   : a statement label in the calling program to which
! control is returned from ffrstr in case of an error
! docr    latest revision : 05/11/86
! docp    programmer      : d.husfeld
! doce

! free-field read : control variables                        ffrctrl
! latest revision : 10/10/86                                 ffrctrl
! iout = output unit number                                  ffrctrl
! iver = verification level                                  ffrctrl
! nor = no run indicator                                     ffrctrl
! cdi = comment delimiting indicator                         ffrctrl
! sdi = string delimiting indicator                          ffrctrl
! cdc = comment delimiting character                         ffrctrl
! sdc = string delimiting character                          ffrctrl
! int = interrupt character                                  ffrctrl
! kcs = keyword character set                                ffrctrl
! kcl = keyword character set length                         ffrctrl
! free-field read : error status table                       ffrerrs
! latest revision : 26/09/86                                 ffrerrs
! errprg = name of routine that detected the error           ffrerrs
! errmsg = error message                                     ffrerrs

! free-field read : error status clearance                   ffreclr
! latest revision : 26/09/86                                 ffreclr
! .. parameters ..
  Use :: nlte_type
  Implicit None

  Integer (i4b) :: kclm
  Parameter (kclm=64)
! ..
! .. scalar arguments ..
  Integer (i4b) :: nc
  Character :: str*(*)
! ..
! .. scalars in common ..
  Integer (i4b) :: iout, iver, kcl
  Logical :: cdi, nor, sdi
  Character :: cdc*1, int*1, sdc*1, errprg*6, errmsg*40, kcs*(kclm)
! ..
! .. local scalars ..
  Integer (i4b) :: i, ierr, l
  Logical :: cdisav
  Character :: ch*1, retchar*4, ret*4
! ..
! .. external subroutines ..
  External :: ffrnch
! ..
! .. intrinsic functions ..
! ..
! .. common blocks ..
  Common /ffrct1/iout, iver, nor, cdi, sdi, kcl
  Common /ffrct2/cdc, sdc, int, kcs
  Common /ffrerr/errprg, errmsg
! ..
  errprg = ' '
  errmsg = ' '

  If (iver>=3) Write (iout, Fmt='(A)', Iostat=ierr) ' **** READ STRING'

  l = len(str)
  str = ' '
  nc = 0

! .....inhibit temporarily the comment delimiter detection
  cdisav = cdi
  cdi = .False.

! .....icollect in str all characters until either str is full
! or the string delimiter (if any is defined) is encountered
100 Continue
  Do i = 1, l
    Call ffrnch(ch, ret)

    Select Case (ret)
    Case ('RET1')
      Go To 120
    Case ('RET2')
      Go To 130
    Case ('RET0')
      Continue
    Case Default
      Stop ' WRONG RETURN CONDITION'
    End Select

    If ((sdi) .And. (ch==sdc)) Go To 110
    str(i:i) = ch
    nc = i
  End Do
110 Continue
  If (nc==0) Go To 100

! .....restore comment delimiter
  cdi = cdisav

! .....verify
  If (iver>=2) Write (iout, Fmt='(A,A)', Iostat=ierr) ' *** STRING FOUND : ', &
    str(1:nc)

  retchar = 'RET0'
  Return

! .....eof condition exit
120 Continue
  retchar = 'RET1'
  Return

! .....error exit
130 Continue
  retchar = 'RET2'
  Return

End Subroutine
Subroutine ffrver(iv)

! docs ffrver
! free-field read : define verification level
! docc    call ffrver ( iv )
! iv     : an integer variable that defines (after being taken
! to the modulus of 10) the new verification level
! docr    latest revision : 10/10/86
! docp    programmer      : d.husfeld
! doce

! free-field read : control variables                        ffrctrl
! latest revision : 10/10/86                                 ffrctrl
! iout = output unit number                                  ffrctrl
! iver = verification level                                  ffrctrl
! nor = no run indicator                                     ffrctrl
! cdi = comment delimiting indicator                         ffrctrl
! sdi = string delimiting indicator                          ffrctrl
! cdc = comment delimiting character                         ffrctrl
! sdc = string delimiting character                          ffrctrl
! int = interrupt character                                  ffrctrl
! kcs = keyword character set                                ffrctrl
! kcl = keyword character set length                         ffrctrl
! free-field read : error status table                       ffrerrs
! latest revision : 26/09/86                                 ffrerrs
! errprg = name of routine that detected the error           ffrerrs
! errmsg = error message                                     ffrerrs

! free-field read : error status clearance                   ffreclr
! latest revision : 26/09/86                                 ffreclr
! .. parameters ..
  Use :: nlte_type
  Implicit None

  Integer (i4b) :: kclm
  Parameter (kclm=64)
! ..
! .. scalar arguments ..
  Integer (i4b) :: iv
! ..
! .. scalars in common ..
  Integer (i4b) :: iout, iver, kcl
  Logical :: cdi, nor, sdi
  Character :: cdc*1, int*1, sdc*1, errprg*6, errmsg*40, kcs*(kclm)
! ..
! .. intrinsic functions ..
! ..
! .. common blocks ..
  Common /ffrct1/iout, iver, nor, cdi, sdi, kcl
  Common /ffrct2/cdc, sdc, int, kcs
  Common /ffrerr/errprg, errmsg
! ..
  errprg = ' '
  errmsg = ' '

  iver = mod(iv, 10)
  If (iout==0) iver = 0
  Return

End Subroutine
Subroutine ffrifc(line)
! dummy subroutine to satisfy external reference in the free-field
! program package. it may be used to handle interface interrupts
! found in the input field which are presently ignored.
! .. scalar arguments ..
  Use :: nlte_type
  Implicit None

  Character :: line*(*)
! ..

  Return
End Subroutine

!----------------------------------------------------------------------
Subroutine ffrloc1(str, retchar)

! docs ffrloc1
! free-field read : locate specified string and position read
! pointer behind it-only first five characters are checked per
! line for first character of str!!!
! docc    call ffrloc ( str,retchar )
! str    : a character expression that specifies the string
! to be located.
! ret    : a statement label in the calling program to which
! control is returned from ffrloc when an eof is found
! ret2   : a statement label in the calling program to which
! control is returned from ffrloc in case of an error
! docr    latest revision : 30/04/97
! docp    programmer      : j.puls
! doce

! free-field read : control variables                        ffrctrl
! latest revision : 10/10/86                                 ffrctrl
! iout = output unit number                                  ffrctrl
! iver = verification level                                  ffrctrl
! nor = no run indicator                                     ffrctrl
! cdi = comment delimiting indicator                         ffrctrl
! sdi = string delimiting indicator                          ffrctrl
! cdc = comment delimiting character                         ffrctrl
! sdc = string delimiting character                          ffrctrl
! int = interrupt character                                  ffrctrl
! kcs = keyword character set                                ffrctrl
! kcl = keyword character set length                         ffrctrl
! free-field read : error status table                       ffrerrs
! latest revision : 26/09/86                                 ffrerrs
! errprg = name of routine that detected the error           ffrerrs
! errmsg = error message                                     ffrerrs
! free-field read : file status variables                    ffrstat
! latest revision : 05/11/86                                 ffrstat
! inp = input unit number                                    ffrstat
! icol = column number in line                               ffrstat
! iseq = line sequence number                                ffrstat
! eof = end-of-file indicator                                ffrstat
! lin = line image                                           ffrstat

! free-field read : error status clearance                   ffreclr
! latest revision : 26/09/86                                 ffreclr
! .. parameters ..
  Use :: nlte_type
  Implicit None

  Integer (i4b) :: kclm
  Parameter (kclm=64)
! ..
! .. scalar arguments ..
  Character :: str*(*)
! ..
! .. scalars in common ..
  Integer (i4b) :: icol, inp, iout, iseq, iver, kcl
  Logical :: cdi, eof, nor, sdi
  Character :: cdc*1, int*1, sdc*1, errprg*6, errmsg*40, lin*80, kcs*(kclm)
! ..
! .. local scalars ..
  Integer (i4b) :: i, ierr, l
  Character :: ch1*1, ch2*1, retchar*4, ret*4
! ..
! .. external subroutines ..
  External :: ffrnch, ffrnch1
! ..
! .. intrinsic functions ..
! ..
! .. common blocks ..
  Common /ffrct1/iout, iver, nor, cdi, sdi, kcl
  Common /ffrct2/cdc, sdc, int, kcs
  Common /ffrerr/errprg, errmsg
  Common /ffrst1/inp, icol, eof, iseq
  Common /ffrst2/lin
! ..
  errprg = ' '
  errmsg = ' '

  If (iver>=3) Write (iout, Fmt='(A,A)', Iostat=ierr) &
    ' **** LOCATING STRING : ', str

  l = len(str)
  ch2 = str(1:1)
100 Continue
  Call ffrnch1(ch1, ret)

  Select Case (ret)
  Case ('RET1')
    Go To 120
  Case ('RET2')
    Go To 130
  Case ('RET0')
    Continue
  Case Default
    Stop ' WRONG RETURN CONDITION'
  End Select

110 Continue
  If (ch1/=ch2) Go To 100
  Do i = 2, l
    Call ffrnch(ch1, ret)

    Select Case (ret)
    Case ('RET1')
      Go To 120
    Case ('RET2')
      Go To 130
    Case ('RET0')
      Continue
    Case Default
      Stop ' WRONG RETURN CONDITION'
    End Select

    If (ch1/=str(i:i)) Go To 110
  End Do

! .....verify
  If (iver>=2) Write (iout, Fmt='(A,I5)', Iostat=ierr) &
    ' *** SPECIFIED STRING LOCATED IN LINE #', iseq

  retchar = 'RET0'
  Return

! .....eof condition exit
120 Continue
  errprg = 'FFRLOC1'
  errmsg = 'EOF WHILE LOCATING ''' // str // ''''
  retchar = 'RET1'
  Return

! .....error exit
130 Continue
  retchar = 'RET2'
  Return

End Subroutine
!----------------------------------------------------------------------
Subroutine ffrnch1(ch, retchar)

! docs ffrnch
! free-field read : get next character in the first five columns
! docc    *** for internal use only
! call ffrnch ( ch,retchar )
! ch     : a character variable that receives in its first po-
! sition the next character found in the input stream
! ret1   : a statement label in the calling program to which
! control is returned from ffrnch when an eof is found
! ret2   : a statement label in the calling program to which
! control is returned from ffrnch in case of an error
! docr    latest revision : 30/04/93
! docp    programmer      : j.puls
! doce

! free-field read : control variables                        ffrctrl
! latest revision : 10/10/86                                 ffrctrl
! iout = output unit number                                  ffrctrl
! iver = verification level                                  ffrctrl
! nor = no run indicator                                     ffrctrl
! cdi = comment delimiting indicator                         ffrctrl
! sdi = string delimiting indicator                          ffrctrl
! cdc = comment delimiting character                         ffrctrl
! sdc = string delimiting character                          ffrctrl
! int = interrupt character                                  ffrctrl
! kcs = keyword character set                                ffrctrl
! kcl = keyword character set length                         ffrctrl
! free-field read : error status table                       ffrerrs
! latest revision : 26/09/86                                 ffrerrs
! errprg = name of routine that detected the error           ffrerrs
! errmsg = error message                                     ffrerrs
! free-field read : file status variables                    ffrstat
! latest revision : 05/11/86                                 ffrstat
! inp = input unit number                                    ffrstat
! icol = column number in line                               ffrstat
! iseq = line sequence number                                ffrstat
! eof = end-of-file indicator                                ffrstat
! lin = line image                                           ffrstat

! free-field read : error status clearance                   ffreclr
! latest revision : 26/09/86                                 ffreclr
! .. parameters ..
  Use :: nlte_type
  Implicit None

  Integer (i4b) :: kclm
  Parameter (kclm=64)
! ..
! .. scalar arguments ..
  Character :: ch*(*)
! ..
! .. scalars in common ..
  Integer (i4b) :: icol, inp, iout, iseq, iver, kcl
  Logical :: cdi, eof, nor, sdi
  Character :: cdc*1, int*1, sdc*1, errprg*6, errmsg*40, lin*80, kcs*(kclm)
! ..
! .. local scalars ..
  Integer (i4b) :: i, iicol
  Character (4) :: retchar, ret
! ..
! .. external subroutines ..
  External :: ffrncd
! ..
! .. common blocks ..
  Common /ffrct1/iout, iver, nor, cdi, sdi, kcl
  Common /ffrct2/cdc, sdc, int, kcs
  Common /ffrerr/errprg, errmsg
  Common /ffrst1/inp, icol, eof, iseq
  Common /ffrst2/lin
! ..
  errprg = ' '
  errmsg = ' '

  If (icol>5) Then
    Call ffrncd(ret)
    Select Case (ret)
    Case ('RET1')
      Go To 110
    Case ('RET2')
      Go To 120
    Case ('RET0')
      Continue
    Case Default
      Stop ' WRONG RETURN CONDITION'
    End Select
  End If

  ch(1:1) = lin(icol:icol)
  icol = icol + 1
  retchar = 'RET0'
  If (.Not. cdi) Return
  If (ch(1:1)/=cdc) Return
! ...comment delimiter found. skip to next comment delimiter or end
! of line and return a blank in ch.
  ch(1:1) = ' '
  iicol = icol
  Do i = iicol, 5
    icol = i
    If (lin(i:i)==cdc) Go To 100
  End Do
100 Continue
  icol = icol + 1
  Return

! .....eof condition exit
110 Continue
  retchar = 'RET1'
  Return

! .....error exit
120 Continue
  retchar = 'RET2'
  Return

End Subroutine
