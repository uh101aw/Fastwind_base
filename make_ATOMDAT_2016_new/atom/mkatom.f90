    PROGRAM MKATOM

!     Initial Version:  2005, Philip Hultzsch:
!                       generates atomic data files for use by WM-basic,
!                       html catalog of atomic data.
!
!     based on PROGRAM UPDATE by Margie Lennon, September 1995 with
!     modifications by Silvia Becker, 1999
!------------------------------------------------------------------------------
!------------------------------------------------------------------------------
!
!     The following files are created in this program
!
!     "unformatted"
!     nl3i......... indices to the long lines list                    UNIT  81
!     nl3a......... sorted wavelengths of above lines                 UNIT  82
!     nl3g......... corresponding gf values of above lines            UNIT  83
!
!     "ascii"
!     nl3info...... bookkeeping information for the nl3* files        UNIT  10
!     generalinfo.. bookkeeping information for all other files       UNIT  20
!
!     atomcol...... collision strengths                               UNIT  11
!     atomlev...... atomic term designations                          UNIT  12
!     atompho...... photoioniation data (Seaton fits)                 UNIT  14
!     atomnl2...... energy level information                          UNIT  15
!     atomgf3...... transition data for lines used in NLTE calc.      UNIT  16
!     atomnldr..... dielectronic recombination data                   UNIT  17
!     atdiestwnel.. ground state ionisation energies
!                   statistical weight of groundstate of highest ion
!                   completeness of term designations, Seaton fits    UNIT  18
!
!     NOTE: When updating atompho or atomcol, take care that the level
!           numbering corresponds to that already given in atomnl2 and
!           atomlev files
!
!------------------------------------------------------------------------------
!------------------------------------------------------------------------------

      USE M_type_definitions,   ONLY: sp,dp,i4b
      USE M_global_settings,    ONLY: KIS,LDRL,KEL,L3,NUMEL,NCALC
      USE M_physical_constants, ONLY: SUNR,AWEIGHT,SOLAR_ABUND,         &
                                      STRING1,STRING2,HTCSSBF

      IMPLICIT NONE

      INTEGER(i4b), PARAMETER :: NII=200,NKIS=10,NKEL=50,NIOCCUP=30,         &
                            NLDRL=100,NLDR=5000,NL3ARR=9000000,         &
                            NEXC=1000

! Max. energy difference allowed for line packing. This is just a guess value.
! If you lower this value, all lines with energies more different that ERGDIFF
! compared to the difference in level energies are discarded. Currently, only
! a few Cl_I lines are removed. If you lower the value to 25000, many Cu_V and
! Ar_V lines will be removed as well. Details are printed out and can be found
! in the file fort.333.

      REAL(dp),     PARAMETER :: ERGDIFF=50000._dp !

      LOGICAL      :: COLEXT,LEVEXT,NL2EXT,NL3EXT,PHOEXT,NLDEXT,        &
                      INFEXT,ALPHAL

      INTEGER(i4b) :: EXC,GETCWD,I,IE,IELZAHL,II,IL,ILEVCON,            &
                      IMACH,IMAX,IO,IQUELDR,IREC,IU,IZ,J,J1,J2,J3,      &
                      JCOUNT,JJ,K,KK,KK1,KLAUF,L,L1,LDADR,LENGTH,LEN_C, &
                      LEN_F,LEVAT,LEVDR,LEVPHI,LEVX,LINZAHL,LLDADR,     &
                      LLEVDR,LLEVPHI,LLG,LLINZAHL,M1,M2,MCOL,MENE,MGEN, &
                      MGF3,MHTM,MION,MLEV,MLIS,MM,MMA,MMB,MMC,MNI,MNL2, &
                      MNL3,MNL3A,MNL3G,MNL3I,MNL3INFO,MNLD,MPHO,MTRA,   &
                      MYSTAT,N,N1,N2,NA,NCOI,NDR,NEL,NELE,NGF3,NIZ,     &
                      NLEVAT,NLL,NNL2,NNL3,NPHO,NQ,NREC,NZQ,NZQUA,Z,ZN, &
                      ZZ,NDATE,NTEMPDATE
      REAL(dp)     :: ANG,ANGOLD,ALPHA,BETA,SA,FREQ,SLOPG,OMEGG,LEVDIF, &
                      HTMLANG(5)

      LOGICAL,      DIMENSION(NII)    :: LOGMERK

      INTEGER(i4b), DIMENSION(NII)    :: NOEL,NOELD,INOEL,INOELD,       &
                                         IELNAM,IONST,                  &
                                         NN,IATN,IS,NCOL,NLEV,NLEVPHI,  &
                                         NLINZAHL,NLEVDR,NLDADR,NELEV,  &
                                         NSRC
      INTEGER(i4b), DIMENSION(NKIS)   :: NOELT,INOELT,JSTATWKT
      INTEGER(i4b), DIMENSION(NKEL)   :: IELNO,IACTIVE,JSTATWK

      INTEGER(i4b), DIMENSION(NLDR)   :: LLOW,LUP,IQUELLP,LEVUNT,LEVOB, &
                                         IVON,INACH,IQUELLC
      REAL(dp),     DIMENSION(NLDR)   :: ALPH,BET,S,DRWELL,DRGF,        &
                                         SLOP,OMEG

      INTEGER(i4b), DIMENSION(NL3ARR) :: IND,IUOLEV
      REAL(dp)    , DIMENSION(NL3ARR) :: GF,ANGST

      INTEGER(i4b), DIMENSION(NLDRL)  :: NINLEV,NLDADL,ILEV,JLEV,LHUP,  &
                                         MAINQ,MAINQ1,JSTATW,JSTATW1,   &
                                         INLEV,IBGMARK
      REAL(dp),     DIMENSION(NLDRL)  :: XLEVEL,XLEVEL1,ZQUANT

      INTEGER(i4b), DIMENSION(NEXC)   :: EXCZ,EXCI,EXCIL,EXCIU
      REAL(dp),     DIMENSION(NEXC)   :: EXCA

      INTEGER(i4b), DIMENSION(NLDRL,NKIS)      :: ISUM

      REAL(dp), DIMENSION(NLDRL*(NLDRL-1)/2)   :: ANGGF,GFGF,SGFGF
      REAL(dp), DIMENSION(NIOCCUP)             :: IOCC
      REAL(dp), DIMENSION(NKEL,NKIS)           :: ENERGY

      INTEGER(i4b), DIMENSION(NLDRL,NLDRL)     :: GF3INDEX
      REAL(sp),     DIMENSION(NL3ARR), TARGET  :: ANG1,ANG2,GF1,GF2
      INTEGER(i4b), DIMENSION(NL3ARR), TARGET  :: IND1,IND2
      REAL(sp),     DIMENSION(:),      POINTER :: ANGIN,ANGOUT,GFIN,GFOUT
      INTEGER(i4b), DIMENSION(:),      POINTER :: INDIN,INDOUT
      LOGICAL                                  :: INIS1=.TRUE.
      LOGICAL                                  :: OPTHTMLCAT=.TRUE.

      CHARACTER          :: SLASH,CION*7,CIOT*7,ATOM(NLDR)*12,TEST*4,   &
                            DATE*8,COLDEF(4)*7,COLDES(4)*10,COMM*2,     &
                            CHARCON(NLDR)*10,DATEION(NII)*8
      CHARACTER(LEN=260) :: SRCDIR,GENNAM,IONNAM,COLNAM,LEVNAM,PHONAM,  &
                            NL2NAM,NL3NAM,GF3NAM,NLDNAM,FILENAM,NL3INAM,&
                            NL3ANAM,NL3GNAM,NL3INFO,HTMNAM,HTMTIT,      &
                            LISNAM,ENENAM,TRANAM,ATDNAM,ABDIRNAM


      DATA MCOL,MLEV,MPHO,MNL2,MGF3,MNLD,MION,MNL3I,MNL3A,MNL3G,MGEN    &
           / 11,  12,  14,  15,  16,  17,  18,   81,   82,   83,  20 /  &
           MNL3INFO,MNL3,MHTM,MENE,MTRA,MLIS                            &
           /     10, 100, 110, 121, 122, 126 /

      DATA HTMLANG / 0,227,501,911,1800 /

!----------------------------------------------------------------------
!----------------------------------------------------------------------

!     CALL GETENV ('ABDIR',ABDIRNAM)
!     WRITE (*,'(/"ABDIRNAM IS ",A,/)') ABDIRNAM
!     IF (LEN_TRIM(ABDIRNAM).EQ.0.OR. &
!         SCAN(ABDIRNAM,'Ascii-Bin-DataInput',BACK=.TRUE.).GT.0) THEN
!       OPTHTMLCAT=.TRUE.
!     ENDIF
!     IF (OPTHTMLCAT) WRITE (*,'(/"CREATING HTML-CATALOG."/)')

      MYSTAT  = GETCWD(SRCDIR)
      IF(MYSTAT.NE.0) THEN
        WRITE(*,*)' -> Error: cannot get information on directory!'
        STOP
      ENDIF

      LENGTH  = LEN_TRIM(SRCDIR)
      SLASH   = '/' ! SRCDIR(LENGTH-19:LENGTH-19)
!     IF (LEN_TRIM(ABDIRNAM).EQ.0) THEN
!       SRCDIR  = SRCDIR(1:LENGTH-37)//'Ascii-Bin-DataInput'//SLASH// &
!                                      'Atomic-Models-ASCII'//SLASH
!     ELSE
!       SRCDIR  = ABDIRNAM(1:LEN_TRIM(ABDIRNAM))//SLASH// &
!                                      'Atomic-Models-ASCII'//SLASH
!     ENDIF

      CALL GETENV ('ATOM_DIR',SRCDIR)
      LENGTH  = LEN_TRIM(SRCDIR)
      SRCDIR  = SRCDIR(1:LENGTH)//SLASH
      LENGTH  = LEN_TRIM(SRCDIR)
      WRITE (*,'(/"NEW ATOM DIR IS ",/,A,/)') SRCDIR(1:LENGTH)

      GENNAM  = SRCDIR(1:LENGTH)//'generalinfo'
      IONNAM  = SRCDIR(1:LENGTH)//'atdiestwnel'
      COLNAM  = SRCDIR(1:LENGTH)//'atomcol'
      LEVNAM  = SRCDIR(1:LENGTH)//'atomlev'
      PHONAM  = SRCDIR(1:LENGTH)//'atompho'
      NL2NAM  = SRCDIR(1:LENGTH)//'atomnl2'
      GF3NAM  = SRCDIR(1:LENGTH)//'atomgf3'
      NLDNAM  = SRCDIR(1:LENGTH)//'atomnldr'

      SRCDIR  = SRCDIR(1:LENGTH-20)//'Line-List-Bin-Source'//SLASH
      LENGTH  = LEN_TRIM(SRCDIR)
      NL3INAM = SRCDIR(1:LENGTH)//'nl3i'
      NL3ANAM = SRCDIR(1:LENGTH)//'nl3a'
      NL3GNAM = SRCDIR(1:LENGTH)//'nl3g'
      NL3INFO = SRCDIR(1:LENGTH)//'nl3info'

!     MYSTAT  = GETCWD(SRCDIR)
!     LENGTH  = LEN_TRIM(SRCDIR)
!     SRCDIR  = SRCDIR(1:LENGTH-19)//'Atomic-Models-Source'//SLASH
      CALL GETENV ('ATOM_SRC',SRCDIR)
      LENGTH  = LEN_TRIM(SRCDIR)
      SRCDIR  = SRCDIR(1:LENGTH)//SLASH
!     PRINT*,SRCDIR
      LENGTH  = LEN_TRIM(SRCDIR)
      FILENAM = SRCDIR(1:LENGTH)
      ATDNAM  = SRCDIR(1:LENGTH)//'_atdiestwnel'
      IF (OPTHTMLCAT) THEN
        HTMNAM  = SRCDIR(1:LENGTH)//'html'//SLASH//'index.html'
        LISNAM  = SRCDIR(1:LENGTH)//'html'//SLASH//'data'//SLASH//'_list.txt'
!       PRINT*,LISNAM
      ENDIF

!     Read exceptions for gf3-data
      EXC = 0
      OPEN (200,FILE=SRCDIR(1:LENGTH)//'_exceptions.gf3')
      DO WHILE (.NOT. EOF(200))
        EXC = EXC + 1
        READ (200,'(I2,I1,X,I2,I2,X,F10.3)') EXCZ(EXC),EXCI(EXC),       &
                                 EXCIL(EXC),EXCIU(EXC),EXCA(EXC)
      ENDDO
      CLOSE (200)


      N = 0
!     Read file atdiestwnel
      OPEN (200,FILE=ATDNAM)
      DO K=1,3; READ(200,*); ENDDO
!     groundstate ionisation energy
      DO K=1,KEL
        READ (200,'(4X,9(F9.1,1X))') (ENERGY(K,I),I=1,KIS)
      ENDDO

      DO K=1,3; READ (200,*); ENDDO
!     statistical weight of the max+1 ionisation stage, not required
      DO K=1,KEL
        READ (200,'(5X,9(I4,6X))') (JSTATWKT(I),I=1,KIS)
        DO I=1,KIS
          IF (JSTATWKT(I).NE.0) JSTATWK(K) = JSTATWKT(I)
        ENDDO
      ENDDO

      DO K=1,3; READ (200,*); ENDDO
!     number of equivalent electrons per level of the corresp. ionis. st.
      DO K=1,KEL
        READ (200,'(5X,9(I3,I3,4X))') (NOELT(I),INOELT(I),I=1,KIS)
        DO I=1, KIS
          IF (NOELT(I).NE.0 .AND. INOELT(I).NE.0) THEN
            N = N + 1
            NOEL (N) = NOELT (I)
            INOEL(N) = INOELT(I)
          ENDIF
        ENDDO
      ENDDO
      CLOSE (200)


      OPEN (MGEN,FILE=GENNAM)
      OPEN (MCOL,FILE=COLNAM)
      OPEN (MLEV,FILE=LEVNAM)
      OPEN (MPHO,FILE=PHONAM)
      OPEN (MNL2,FILE=NL2NAM)
      OPEN (MGF3,FILE=GF3NAM)
      OPEN (MNLD,FILE=NLDNAM)
      IF (OPTHTMLCAT) THEN
        OPEN (MHTM,FILE=HTMNAM)
        OPEN (MLIS,FILE=LISNAM)
      ENDIF


      WRITE (*,'("* ion     atomnl3   atomgf3   atomnl2   atomlev   ",  &
                "atompho   atomcol  atomnldr *")')

      CALL DATE_AND_TIME(DATE)
      IF (OPTHTMLCAT) THEN
        HTMTIT = 'WMbasic Atomic Database'
        CALL HTMLHEAD(MHTM,HTMTIT)
        CALL COLORR
      ENDIF


      N = 0
      DO IL = 1, NLDRL - 1
        DO IU = IL + 1, NLDRL
          N = N + 1
          GF3INDEX(IL,IU) = N
        ENDDO
      ENDDO


      N   = 0
      ZN  = 0
      LLG = 0

      DO Z = 1, KEL

        IF (OPTHTMLCAT) THEN
          WRITE (MHTM,'("<TR><TD ALIGN=""CENTER"" BGCOLOR=""SILVER"">",I2,        &
                       "</TD><TD></TD>")') Z
        ENDIF

        IELNO(Z)   = Z
        IACTIVE(Z) = 1

        IMAX = MIN(Z,NKIS)
        NELE = 0
        DO I = 1, IMAX

          NPHO = 0
          NGF3 = 0
          NNL2 = 0

          CIOT = STRING1(Z)//' '//STRING2(I)(2:5)
          IF (CIOT(1:1).EQ.' ') CIOT = CIOT(2:7)//' '
          CION = STRING1(Z)//STRING2(I)
          IF (CION(1:1).EQ.' ') CION = CION(2:7)//' '
          LEN_C = LEN_TRIM(CION)
          FILENAM(LENGTH+1:LENGTH+LEN_C+1) = CION(1:LEN_C)
          LEN_F = LENGTH+LEN_C

          INQUIRE (FILE=FILENAM(1:LEN_F)//'.col' ,EXIST=COLEXT)
          INQUIRE (FILE=FILENAM(1:LEN_F)//'.lev' ,EXIST=LEVEXT)
          INQUIRE (FILE=FILENAM(1:LEN_F)//'.nl3' ,EXIST=NL3EXT)
          INQUIRE (FILE=FILENAM(1:LEN_F)//'.pho' ,EXIST=PHOEXT)
          INQUIRE (FILE=FILENAM(1:LEN_F)//'.nldr',EXIST=NLDEXT)
          INQUIRE (FILE=FILENAM(1:LEN_F)//'.info',EXIST=INFEXT)


          IF (NL3EXT) THEN
            N = N + 1
            IF (ZZ .NE. Z) THEN
              ZZ = Z
              ZN = ZN + 1
              IOCC(ZN) = ZZ
            ENDIF
            IELNAM(N) = ZZ
            IONST(N)  = I
          ENDIF

!     COLLISION DATA
          IF (COLEXT) THEN
            NCOL(N) = 0
            OPEN (200,FILE=FILENAM(1:LEN_F)//'.col')
            DO WHILE (.NOT. EOF(200))
              NCOL(N) = NCOL(N) + 1
              READ (200,*)
            ENDDO
            REWIND (200)
            WRITE (MCOL,'(2(2X,I4))') N,NCOL(N)
            DO J = 1, NCOL(N)
              READ (200,*) J1,J2,SLOPG,OMEGG,JJ
              WRITE (MCOL,'(2(2x,I4),2X,E8.2,2X,E8.2,2X,I4)')           &
                           J1,J2,SLOPG,OMEGG,JJ
            ENDDO
            CLOSE (200)
          ENDIF

!     ENERGIES, QUANTUM NOS AND STAT. WT
          IF (LEVEXT) THEN
            NNL2   = 0
            NLEV   = 0
            LHUP   = 0
            ZQUANT = 0
            JLEV   = 0
            ILEV   = 0
            OPEN (200,FILE=FILENAM(1:LEN_F)//'.lev')
            DO WHILE (.NOT. EOF(200))
              NNL2 = NNL2 + 1
              READ (200,'(I4,2x,F12.1,2X,I4,2X,A12)')                   &
                   JLEV(NNL2),XLEVEL(NNL2),JSTATW(NNL2),ATOM(NNL2)
              READ (ATOM(NNL2),'(I2,A)') MAINQ(NNL2)
              IF (LEN_TRIM(ATOM(NNL2)).GT.2) THEN
                NLEV(N) = NLEV(N) + 1
                ILEV(NLEV(N)) = NNL2
                READ (ATOM(NNL2)(4:5),'(I2)') NZQ
                ZQUANT(NLEV(N)) = REAL(NZQ,dp)
                READ (ATOM(NNL2)(10:11),'(I2)') LHUP(NLEV(N))
              ENDIF
            ENDDO
            CLOSE (200)
            IF (NNL2.GT.0) WRITE (MNL2,'(2(2X,I4))') N,NNL2
            IF (NLEV(N).GT.0) WRITE (MLEV,'(2(2X,I4))') N,NLEV(N)
            DO J = 1, NNL2
              WRITE (MNL2,'(F12.1,2X,I2,2X,I4)') XLEVEL(J),MAINQ(J),JSTATW(J)
              IF (LEN_TRIM(ATOM(J)).GT.2)                               &
                WRITE (MLEV,'(I10,2x,A12)') JLEV(J),ATOM(J)
            ENDDO
          ENDIF

!     TRANSITION DATA
          IF (NL3EXT) THEN
            IF(INIS1) THEN
              ANGIN  => ANG1
              GFIN   => GF1
              INDIN  => IND1
              ANGOUT => ANG2
              GFOUT  => GF2
              INDOUT => IND2
            ELSE
              ANGIN  => ANG2
              GFIN   => GF2
              INDIN  => IND2
              ANGOUT => ANG1
              GFOUT  => GF1
              INDOUT => IND1
            ENDIF
            INIS1=.NOT.INIS1

            ANGGF = 0._dp
            OPEN (200,FILE=FILENAM(1:LEN_F)//'.nl3')
            DO J = 1, 12
              READ (200,*)
            ENDDO
            TEST = '    '
            NNL3 = 0
            DO WHILE (TEST.NE.'++++')
              READ (200,'(A4)') TEST
              NNL3 = NNL3 + 1
            ENDDO
            NNL3 = NNL3 - 1
            LLG = LLG + NNL3
            REWIND (200)
            DO J = 1, 12
              IF (J.EQ.3) THEN                   ! Change date to nl3 file
                READ (200,'(21X,A8)') DATEION(N)
              ELSE
                READ (200,*)
              ENDIF
            ENDDO
            K = 1
            ANGOLD = 0.
            DO J = 1, NNL3
              READ (200,*) ANG,GF(J),IND(J)
              IF (ANGOLD.GT.ANG) THEN
                PRINT*,ANGOLD, ANG
                STOP '  ANGOLD .GT. ANG !!!'
              ENDIF  
              ANGOLD = ANG
              IU = IND(J)-(IND(J)/100)*100
              IL = IND(J)/100-(IND(J)/10000)*100
              IF (IL.EQ.0) THEN
                STOP "STOP IL=0"
 !            ELSE IF (IU.EQ.0.AND.IL.LT.NLEV(N)) THEN
 !              LEVDIF = XLEVEL(NLEV(N))-XLEVEL(IL)
 !              IF (LEVDIF.GT.1.E8_dp/ANG+2._dp*ERGDIFF) THEN
 !                WRITE (333,'("NOTICE:  BAD WAVELENGTH FOUND FOR ",I8, &
 !                  " - CHECK  ",F11.1," WITH ",F11.1)') IND(J),ANG,1.E8_dp/LEVDIF
 !                WRITE (*  ,'("NOTICE:  BAD WAVELENGTH FOUND FOR ",I8, &
 !                  " - CHECK  ",F11.1," WITH ",F11.1)') IND(J),ANG,1.E8_dp/LEVDIF
!!                CYCLE
 !              ENDIF
              ELSE IF (IU.GT.0) THEN
                LEVDIF = XLEVEL(IU)-XLEVEL(IL)
                IF (ABS(LEVDIF-1.E8_dp/ANG).GT.ERGDIFF) THEN
                  WRITE (333,'("WARNING: BAD WAVELENGTH FOUND FOR ",I8, &
                    " - IS ",F14.3," SHOULD BE ",F14.3)') IND(J),ANG,1.E8_dp/LEVDIF
                  WRITE (*  ,'("WARNING: BAD WAVELENGTH FOUND FOR ",I8, &
                    " - IS ",F11.1," SHOULD BE ",F11.1)') IND(J),ANG,1.E8_dp/LEVDIF
!                  CYCLE
!                  STOP
                ENDIF
                IF (ANGGF(GF3INDEX(IL,IU)).EQ.0.) THEN
                  ANGGF(GF3INDEX(IL,IU)) = ANG
                  GFGF (GF3INDEX(IL,IU)) = GF(J)
                  SGFGF(GF3INDEX(IL,IU)) = GF(J)
                ELSE
                  SGFGF(GF3INDEX(IL,IU)) = SGFGF(GF3INDEX(IL,IU)) + GF(J)
                  IF (GF(J).GE.GFGF(GF3INDEX(IL,IU))) THEN
                    ANGGF(GF3INDEX(IL,IU)) = ANG
                    GFGF (GF3INDEX(IL,IU)) = GF(J)
                  ENDIF
                ENDIF
              ENDIF
              DO WHILE (ANGIN(K).LT.ANG .AND. ANGIN(K).GT.0.)
                ANGOUT(K+J-1) = ANGIN(K)
                GFOUT (K+J-1) = GFIN (K)
                INDOUT(K+J-1) = INDIN(K)
                K = K + 1
              ENDDO
              ANGOUT(K+J-1) = ANG
              GFOUT (K+J-1) = GF (J)
              INDOUT(K+J-1) = IND(J)
            ENDDO
            CLOSE (200)
            DO WHILE (K.LE.LLG-NNL3)
              ANGOUT(K+NNL3) = ANGIN(K)
              GFOUT (K+NNL3) = GFIN (K)
              INDOUT(K+NNL3) = INDIN(K)
              K = K + 1
            ENDDO
!            print*,'test1',llg,nnl3

            
            NGF3 = 0
            DO IL = 1, NLDRL - 1
              DO IU = IL + 1, NLDRL
                IF (ANGGF(GF3INDEX(IL,IU)).GT.0.) THEN
                  NGF3 = NGF3 + 1
                ENDIF
              ENDDO
            ENDDO
            DO JJ = 1, EXC
              IF (Z.EQ.EXCZ(JJ).AND.I.EQ.EXCI(JJ))                      &
                      ANGGF(GF3INDEX(EXCIL(JJ),EXCIU(JJ))) = EXCA(JJ)
            ENDDO
            IF (NGF3.GT.0) WRITE (MGF3,'(2(2X,I4))') N, NGF3
            DO IL = 1, NLDRL - 1
              DO IU = IL + 1, NLDRL
                IF (ANGGF(GF3INDEX(IL,IU)).GT.0.) THEN
                  WRITE (MGF3,'(2X,I2,I1,"0",2(I2.2),2X,F12.3,2X,E8.2)')&
                        Z,I,IL,IU,ANGGF(GF3INDEX(IL,IU)),               &
                        SGFGF(GF3INDEX(IL,IU))
                ENDIF
              ENDDO
            ENDDO

          ENDIF

!     PHOTOIONISATION DATA
          IF (PHOEXT) THEN
            NPHO = 0
            OPEN (200,FILE=FILENAM(1:LEN_F)//'.pho')
            DO WHILE (.NOT. EOF(200))
              NPHO = NPHO + 1
              READ (200,*) LLOW(NPHO),LUP(NPHO),ALPH(NPHO),BET(NPHO), &
                           S(NPHO),IQUELLP(NPHO)
            ENDDO
            CLOSE (200)
            WRITE (MPHO,'(2(2X,I4))') N,NPHO
            DO J = 1, NPHO
              WRITE (MPHO,'(2(2X,I4),3X,F8.4,1X,F7.2,2X,F7.2,2X,I2)')   &
                LLOW(J),LUP(J),ALPH(J),BET(J),S(J),IQUELLP(J)
            ENDDO
          ENDIF

!     DIELECTRONIC RECOMBINATION
          IF (NLDEXT) THEN
            NLEVDR = 0
            NLDADR = 0
            OPEN (200,FILE=FILENAM(1:LEN_F)//'.nldr')
            DO WHILE (.NOT. EOF(200))
              NLEVDR(N) = NLEVDR(N) + 1
              READ (200,*) NINLEV(NLEVDR(N)),NLDADL(NLEVDR(N))
              NLDADR(N) = NLDADR(N) + NLDADL(NLEVDR(N))
              DO J = 1, NLDADL(NLEVDR(N))
                READ (200,*)
              ENDDO
            ENDDO
            REWIND (200)
            WRITE (MNLD,'(3(2X,I4))') N,NLEVDR(N),NLDADR(N)
            DO J = 1, NLEVDR(N)
              READ (200,*)
              WRITE (MNLD,'(2(2X,I4))') NINLEV(J),NLDADL(J)
              DO JJ = 1, NLDADL(J)
                READ (200,*) DRWELL(JJ),DRGF(JJ)
                WRITE (MNLD,'(2X,F10.1,2X,E10.3)')                      &
                  DRWELL(JJ),DRGF(JJ)
              ENDDO
            ENDDO
            CLOSE (200)
          ENDIF

!     ADDITIONAL INFO: SOURCE, EXCEPTIONS
          IF (INFEXT) THEN
            OPEN (200,FILE=FILENAM(1:LEN_F)//'.info')
            READ (200,*) NSRC(N)
            READ (200,*) NCOI
            CLOSE (200)
          ENDIF


          IF (LEVEXT.OR.NL2EXT.OR.NL3EXT.OR.PHOEXT.OR.NLDEXT) THEN
!     GERNERALINFO
            WRITE (*,'("*",1X,A2,1X,A4,6(2X,I6,2X),2X,I6,1X,"*")')      &
                 STRING1(Z),STRING2(I)(2:5),NNL3,NGF3,NNL2,NLEV(N),     &
                 NPHO,NCOL(N),NLEVDR(N)
            WRITE (MGEN,101) N,Z,I,NCOL(N),NLEV(N),NPHO,NGF3,NLEVDR(N), &
                 NLDADR(N),NNL2,NSRC(N)
101         FORMAT (11I4)

!     HTML-CATALOG
            IF (OPTHTMLCAT) THEN
              CALL HTMLLINK(MHTM,CION,CIOT,COLDEF(NCOI))
            ENDIF

            NELE = NELE + 1

          ENDIF

        ENDDO

        IF (OPTHTMLCAT) THEN
          DO I = NELE, 7
            WRITE(MHTM,'("  <TD ALIGN=""CENTER"" BGCOLOR=""#EDEED9""></TD>")')
          ENDDO
          WRITE(MHTM,'("</TR>")')
        ENDIF

      ENDDO

      NDATE = 0
      DO I=1,N
        READ (DATEION(I),'(I8)') NTEMPDATE
        IF (NTEMPDATE.GT.NDATE) NDATE = NTEMPDATE
      ENDDO
      WRITE (DATE,'(I8)') NDATE
      IF (OPTHTMLCAT) THEN
        CALL HTMLFOOT(MHTM,COLDEF,COLDES,DATE,HTMLANG)
      ENDIF

      WRITE(*,'()')
      WRITE(*,'("WRITING BINARY FILES ...")')
      WRITE(*,'()')

      OPEN (MNL3I,FORM='UNFORMATTED',FILE=NL3INAM,RECL=L3)
      OPEN (MNL3A,FORM='UNFORMATTED',FILE=NL3ANAM,RECL=L3)
      OPEN (MNL3G,FORM='UNFORMATTED',FILE=NL3GNAM,RECL=L3)

      DO I = 1, LLG/L3
        J  = (I - 1) * L3
        WRITE (MNL3I) (INDOUT(J+JJ),JJ=1,L3)
        WRITE (MNL3A) (ANGOUT(J+JJ),JJ=1,L3)
        WRITE (MNL3G) (GFOUT (J+JJ),JJ=1,L3)
        do jj=1,l3-1
        if(angout(j+jj).gt.angout(j+jj+1)) then
          print*,j+jj,angout(j+jj),angout(j+jj+1)
          print*,'record ',i
          stop' error in angout (intermediate record)'
        endif
        enddo   
!        print*,'complete',i,j,l3
      ENDDO

      J = (LLG/L3)*L3+1
      WRITE (MNL3I) (INDOUT(JJ),JJ=J,LLG)
      WRITE (MNL3A) (ANGOUT(JJ),JJ=J,LLG)
      WRITE (MNL3G) (GFOUT (JJ),JJ=J,LLG)
!      print*,'last',j,llg
      do jj=j,llg-1
        if(angout(jj).gt.angout(jj+1)) then
          print*,jj,angout(jj),angout(jj+1)
          stop' error in angout (last record)'
        endif
      enddo   
      
      CLOSE (MNL3I)
      CLOSE (MNL3A)
      CLOSE (MNL3G)

      OPEN  (MNL3INFO,FILE=NL3INFO)
      WRITE (MNL3INFO,'(I8,I4,I8)') LLG, LLG/L3, LLG-(LLG/L3)*L3
      CLOSE (MNL3INFO)

      IF (OPTHTMLCAT) THEN
        WRITE(*,'()')
        WRITE(*,'("WRITING ASCII NL3 FILES ...")')
        WRITE(*,'()')

        WRITE (NL3NAM,'(A,"html",A1,"data",A1,"lines",I1,".nl3")')         &
               SRCDIR(1:LENGTH),SLASH,SLASH,1
        OPEN  (200,FILE=NL3NAM)
        J = 2
        DO I = 1, LLG
          IF (J.LE.5) THEN
            IF (ANGOUT(I).GT.HTMLANG(J)) THEN
              CLOSE (200)
              WRITE (NL3NAM,'(A,"html",A1,"data",A1,"lines",I1,".nl3")')   &
                   SRCDIR(1:LENGTH),SLASH,SLASH,J
              OPEN  (200,FILE=NL3NAM)
              J = J + 1
            ENDIF
          ENDIF
          WRITE (200,'(F12.3,2X,I8.8,2X,E8.2)') ANGOUT(I),INDOUT(I),GFOUT(I)
        ENDDO
        CLOSE (200)

        CLOSE (MHTM)
!----------------------------------------------------------------------
!----------------------------------------------------------------------
!       WRITE HTML-INFO FOR IONS
        WRITE(*,'()')
        WRITE(*,'("WRITING HTML FILES FOR IONS ...")')
        WRITE(*,'()')

        II = N

        REWIND (MGEN)

!       READ THE BOOKKEEPING INFORMATION FOR THE ASCII FILES
        READ (MGEN,101) NN(1),IATN(1),IS(1),NCOL(1),NLEV(1),NLEVPHI(1),   &
                        NLINZAHL(1),NLEVDR(1),NLDADR(1),NELEV(1),NSRC(1)

        DO N = 1, II

          IF (N.LT.II) THEN
            READ (MGEN,101) NN(N+1),IATN(N+1),IS(N+1),NCOL(N+1),NLEV(N+1),&
                            NLEVPHI(N+1),NLINZAHL(N+1),NLEVDR(N+1),       &
                            NLDADR(N+1),NELEV(N+1),NSRC(N+1)
            IF (N+1 .NE. NN(N+1)) STOP 'N+1 .NE. NN(N+1)!!'
          ENDIF

          Z    = IATN(N)
          I    = IS(N)

          CIOT = STRING1(Z)//' '//STRING2(I)(2:5)
          IF (CIOT(1:1).EQ.' ') CIOT = CIOT(2:7)//' '
          CION = STRING1(Z)//STRING2(I)
          IF (CION(1:1).EQ.' ') CION = CION(2:7)//' '
          LEN_C = LEN_TRIM(CION)
          FILENAM(LENGTH+1:LENGTH+LEN_C+1) = CION(1:LEN_C)
          LEN_F = LENGTH+LEN_C

          WRITE (MLIS,'(A4,"-",A2,"-",A2,1X,A7,1X,A7)') DATEION(N)(1:4),  &
            DATEION(N)(5:6),DATEION(N)(7:8),CIOT(1:LEN_C),CION(1:LEN_C)

          ENENAM  = SRCDIR(1:LENGTH)//'html'//SLASH//'data'//SLASH//      &
                    'ener.'//CION(1:LEN_C)
          TRANAM  = SRCDIR(1:LENGTH)//'html'//SLASH//'data'//SLASH//      &
                    'tran.'//CION(1:LEN_C)
          OPEN (MENE,FILE=ENENAM)
          OPEN (MTRA,FILE=TRANAM)

          DO KLAUF = 1, NII
            NOELD(KLAUF)  = NOEL(KLAUF)
            INOELD(KLAUF) = INOEL(KLAUF)
          ENDDO

          REWIND (MCOL)
          REWIND (MLEV)
          REWIND (MPHO)
          REWIND (MNL2)
          REWIND (MGF3)
          REWIND (MNLD)

          LEN_F = LENGTH + LEN_C + 19  !14 if .new is removed
          FILENAM(LENGTH+1:LEN_F) = 'html/cat/'//CION(1:LEN_C)//'.html'
          OPEN (25,FILE=FILENAM(1:LEN_F))
          HTMTIT = CIOT(1:LEN_C)//'           '
          CALL HTMLHEAD(25,HTMTIT)

          IZ      = NCOL(N)
          LEVAT   = NLEV(N)
          LEVPHI  = NLEVPHI(N)
          LINZAHL = NLINZAHL(N)
          LEVDR   = NLEVDR(N)
          LDADR   = NLDADR(N)

          IND     = 0
          LHUP    = 0
          ALPHA   = 0._dp
          BETA    = 0._dp
          SA      = 0._dp

!       COLLISIONAL DATA                           (atomcol)
          IF (IZ.NE.0) CALL COLREAD (IVON,INACH,SLOP,OMEG,IQUELLC)

!       PHOTOIONISATION DATA                       (atompho)
          IF (LEVPHI.NE.0) CALL PHOREAD (LLOW,LUP,ALPH,BET,S,IQUELLP)

!       LEVEL DESIGNATIONS                         (atomlev)
          IF (LEVAT.NE.0) CALL ATREAD (IND,ATOM)

!       ENERGIES, QUANTUM NO, STAT.WT.             (atomnl2)
          IF (NELEV(N).NE.0) THEN
            NREC = N
            CALL READXL (XLEVEL,MAINQ,JSTATW)
          ENDIF

!       TRANSITION DATA                            (atomgf3)
          IF (LINZAHL.NE.0) CALL RLINE (ANGST,GF,IUOLEV)

!       DIELECTRONIC RECOMBINATION                 (atomnldr)
          IF ((LEVDR.NE.0).AND.(LDADR.NE.0))                              &
            CALL RLINEDR (INLEV,IBGMARK,DRWELL,DRGF)

!       ENERGIES, QUANTUM NO, STAT.WT. NEXT LEVEL  (atomnl2)
          IF (NREC.NE.II .AND. NELEV(N+1).NE.0) THEN
            REWIND (MNL2)
            NREC  = N + 1
            CALL READXL (XLEVEL1,MAINQ1,JSTATW1)
          ENDIF


          WRITE (25,900) CION(1:LEN_C),CION(1:LEN_C),CIOT(1:LEN_C),       &
                         STRING1(Z),STRING2(I)(2:5),CION(1:LEN_C),Z,      &
                         NELEV(N),N
!       FORMAT FOR LOCAL CATALOG
 900      FORMAT('<A HREF="../plot/',A,'.pdf"><IMG SRC="',A,           &
                 '.png" ALT="GROTRIAN DIAGRAM OF ',A,'"></A><BR>'/,       &
                 '<PRE>'/,1X,108('*')/,1X,'  ELEMENT = ',A2,1X,A4,1X,     &
                 '(<A HREF="../../',A,'.nl3">NL3-FILE</A>) : '            &
                 'ELEMENT NUMBER = ',I2,'   : NUMBER OF LEVELS = ',I2,    &
                 ' : ION NUMBER = ',I3/,1X,108('*'))
!       FORMAT FOR INTERNET CATALOG
 901      FORMAT(1X,'<A HREF=../plot/',A,'.pdf><IMG SRC=',A,              &
                 '.png ALT="GROTRIAN DIAGRAM OF ',A,'"></A><BR>'/,        &
                 1X,108('*')/,1X,'  ELEMENT = ',A2,1X,A4,1X,              &
                 '(<A HREF=javascript:void(0); onClick=''var a = '        &
                 'window.open("../data/',A,'.nl3","smallnl3file2",'       &
                 '"toolbar=0,location=0,directories=0,status=0,'          &
                 'menubar=0,scrollbars=yes,copyhistory=0,resizable=yes,'  &
                 'width=400"); a.focus();''>NL3-FILE</A>) : '             &
                 'ELEMENT NUMBER = ',I2,'   : NUMBER OF LEVELS = ',I2,    &
                 ' : ION NUMBER = ',I3/,1X,108('*'))
          WRITE (25,105) LEVPHI
 105      FORMAT(20X,'ENERGY-LEVELS',25X,'PHOTOIONISATION: NUMBER = ',I4/)
          WRITE (25,110)
 110      FORMAT(1X,'LEVEL       TERMS       STAT-WEIGHT',3X,             &
              'ENERGY ',10X,'TO',7X,'ENERGY',4X,'ALPHA   BETA  S',4X,     &
              'COMMENT',2X,'SOURCE')
          WRITE (25,120)
 120      FORMAT(9X,'n',1X,'l',1X,'Z',1X,'M',1X,'L',1X,'P',1X,            &
              'I',7X,'g',9X,'(1/CM)')


          IF (LEVAT.NE.0) THEN
            CALL PHINIT (ISUM,LHUP,ZQUANT)
          ENDIF

          DO L=1,NELEV(N)
            J1=0
            J2=0
            IF (LEVAT.NE.0) THEN
              DO KK=1,LEVAT
                IF (L.EQ.IND(KK)) THEN
                  J1=KK
                ENDIF
              ENDDO
            ENDIF

            MM=0
            IF (LEVPHI.NE.0) THEN
              DO KK=1,LEVPHI
                IF (L.EQ.LLOW(KK)) THEN
                  J2=KK
                  MM=MM+1
                ENDIF
              ENDDO
              MM=MM-1
            ENDIF

!       BEGINNING OF OUTPUT

            IF ((J1.NE.0).AND.(J2.NE.0)) THEN
              COMM='SA'
              READ (ATOM(J1)(4:5),'(I2)') IELZAHL
              READ (ATOM(J1)(10:11),'(I2)') ILEVCON
              CALL PHOTOCS (ALPHA,BETA,SA,LEVX,FREQ)
              IF (ALPHA.LT.1.E-10_dp) ALPHA = ALPHA * 1.E+18_dp
              WRITE (MENE,150) L,ATOM(J1)(2:2),ATOM(J1)(3:3),IELZAHL,     &
                   ATOM(J1)(7:7),ATOM(J1)(8:8),ATOM(J1)(9:9),ILEVCON,     &
                   ATOM(J1)(12:12),JSTATW(L),XLEVEL(L),STRING1(Z),        &
                   STRING2(I+1),LEVX,FREQ,ALPHA,BETA,SA,COMM,             &
                   IQUELLP(J2-MM)
              WRITE (25,151) L,ATOM(J1)(2:2),ATOM(J1)(3:3),IELZAHL,       &
                   ATOM(J1)(7:7),ATOM(J1)(8:8),ATOM(J1)(9:9),ILEVCON,     &
                   ATOM(J1)(12:12),JSTATW(L),XLEVEL(L),STRING1(Z),        &
                   STRING2(I+1),LEVX,FREQ,ALPHA,BETA,SA,COMM,             &
                   IQUELLP(J2-MM)
150           FORMAT(3X,I2,4X,A1,1X,A1,I2,1X,A1,1X,A1,1X,A1,I2,A1,4X,     &
                   I3,6X,F9.1,6X,A2,A4,'/',I2,2X,F9.1,2X,F7.3,1X,F6.2,    &
                   1X,F4.2,4X,A2,6X,I2)
151           FORMAT(3X,I2,4X,A1,1X,A1,I2,1X,A1,1X,A1,1X,A1,I2,A1,4X,     &
                   I3,6X,F9.1,6X,A2,A4,'/',I2,2X,F9.1,2X,F7.3,1X,F6.2,    &
                   1X,F4.2,4X,A2,6X,I2)
              IF (MM.NE.0) THEN
                DO KK=1,MM
                  CALL PHOTOCS (ALPHA,BETA,SA,LEVX,FREQ)
                  IF (ALPHA.LT.1.E-10_dp) ALPHA = ALPHA * 1.E+18_dp
                  WRITE (25,200) STRING1(Z),STRING2(I+1),LEVX,FREQ,       &
                       ALPHA,BETA,SA,COMM,IQUELLP(J2-MM+KK)
200               FORMAT(51X,A2,A4,'/',I2,2X,F9.1,2X,F7.3,1X,F6.2,        &
                         1X,F4.2,4X,A2,6X,I2)
                ENDDO
              ENDIF
            ENDIF

            IF ((J1.EQ.0).AND.(J2.NE.0)) THEN
              COMM='SA'
              CALL PHOTOCS (ALPHA,BETA,SA,LEVX,FREQ)
              IF (ALPHA.LT.1.E-10_dp) ALPHA = ALPHA * 1.E+18_dp
              WRITE (MENE,160) L,MAINQ(L),JSTATW(L),XLEVEL(L),            &
                   STRING1(Z),STRING2(I+1),LEVX,FREQ,ALPHA,BETA,SA,       &
                   COMM,IQUELLP(J2-MM)
              WRITE (25,160) L,MAINQ(L),JSTATW(L),XLEVEL(L),STRING1(Z),   &
                   STRING2(I+1),LEVX,FREQ,ALPHA,BETA,SA,COMM,             &
                   IQUELLP(J2-MM)
160           FORMAT(3X,I2,3X,I2,17X,I3,6X,F9.1,6X,A2,A4,'/',             &
                   I2,2X,F9.1,2X,F7.3,1X,F6.2,1X,F4.2,4X,A2,6X,I2)
              IF (MM.NE.0) THEN
                DO KK=1,MM
                  CALL PHOTOCS (ALPHA,BETA,SA,LEVX,FREQ)
                  IF (ALPHA.LT.1.E-10_dp) ALPHA = ALPHA * 1.E+18_dp
                  WRITE (25,200) STRING1(Z),STRING2(I+1),LEVX,FREQ,ALPHA, &
                         BETA,SA,COMM,IQUELLP(J2-MM+KK)
                ENDDO
              ENDIF
            ENDIF

            IF ((J1.NE.0).AND.(J2.EQ.0)) THEN
              COMM='GM'
              READ (ATOM(J1)(4:5),'(I2)') IELZAHL
              READ (ATOM(J1)(10:11),'(I2)') ILEVCON
              CALL PHOTOCS (ALPHA,BETA,SA,LEVX,FREQ)
              ALPHA = ALPHA * 1.E+18_dp
              WRITE (MENE,150) L,ATOM(J1)(2:2),ATOM(J1)(3:3),IELZAHL,     &
                   ATOM(J1)(7:7),ATOM(J1)(8:8),ATOM(J1)(9:9),ILEVCON,     &
                   ATOM(J1)(12:12),JSTATW(L),XLEVEL(L),STRING1(Z),        &
                   STRING2(I+1),LEVX,FREQ,ALPHA,BETA,SA,COMM
              WRITE (25,151) L,ATOM(J1)(2:2),ATOM(J1)(3:3),IELZAHL,       &
                   ATOM(J1)(7:7),ATOM(J1)(8:8),ATOM(J1)(9:9),ILEVCON,     &
                   ATOM(J1)(12:12),JSTATW(L),XLEVEL(L),STRING1(Z),        &
                   STRING2(I+1),LEVX,FREQ,ALPHA,BETA,SA,COMM
            ENDIF

            IF ((J1.EQ.0).AND.(J2.EQ.0)) THEN
              COMM='GA'
              LEVX=1
              CALL PHOTOCS (ALPHA,BETA,SA,LEVX,FREQ)
              ALPHA = ALPHA * 1.E+18_dp
              WRITE (MENE,160) L,MAINQ(L),JSTATW(L),XLEVEL(L),            &
                   STRING1(Z),STRING2(I+1),LEVX,FREQ,ALPHA,BETA,SA,       &
                   COMM
              WRITE (25,160) L,MAINQ(L),JSTATW(L),XLEVEL(L),STRING1(Z),   &
                   STRING2(I+1),LEVX,FREQ,ALPHA,BETA,SA,COMM
            ENDIF


          ENDDO

          CLOSE (MENE)

          WRITE (25,'()')
          WRITE (25,510) STRING1(Z),STRING2(I)(2:5),LINZAHL,IZ,LDADR
510       FORMAT(/1X,108('-'),/1X,'  ELEMENT = ',A2,1X,A4,' : ',          &
               'NUMBER OF LINEDATA-SETS:',I3/,1X,'    LEVEL',             &
               4X,'gf-VALUE',4X,'WAVELENGTH',3X,'COLLISION-STRENGTHS: ',  &
               'NUMBER=',I3,10X,'DR-DATASETS: NUMBER =',I3)
          WRITE (25,520)
520       FORMAT(3X,'LOW   UP',                                           &
               14X,'[ANGSTROEM]',4X,'OMEGA0',4X,'SLOPE',5X,'COMM. SRC.',  &
               5X,'f-VALUE',2X,'WAVELENGTH',2X,'COMM.  SRC.'/1X,108('-'))

          DO L=1,NELEV(N)
            J1=0
            J2=0
            J3=0
            MMA=0
            MMB=0
            MMC=0

            IF(LINZAHL.NE.0) THEN
              DO KK=1,LINZAHL
                WRITE (CHARCON(KK),'(I8)') IUOLEV(KK)
                READ  (CHARCON(KK)(5:6),'(I3)') LEVUNT(KK)
                READ  (CHARCON(KK)(7:8),'(I4)') LEVOB(KK)
                IF (L.EQ.LEVUNT(KK)) THEN
                  J1=KK
                  MMA=MMA+1
                ENDIF
              ENDDO
            ENDIF

            IF(IZ.NE.0) THEN
              DO KK=1,IZ
                IF(L.EQ.IVON(KK)) THEN
                  J2=KK
                  MMB=MMB+1
                ENDIF
              ENDDO
            ENDIF

            IF(LEVDR.NE.0) THEN
              DO KK=1,LEVDR
                IF(L.EQ.INLEV(KK)) THEN
                  IF(KK.NE.LEVDR) THEN
                    J3=KK
                    MMC=(IBGMARK(J3+1)-IBGMARK(J3))
                  ELSE
                    J3=KK
                    MMC=LDADR-IBGMARK(J3)+1
                  ENDIF
                ENDIF
              ENDDO
            ENDIF

!       -----------------------------------------------------------
!       DIE FOLGENDEN AUSGABESTATEMENTS UMFASSEN ALLE MOEGLICHEN
!       FAELLE VON DATENMAECHTIGKEIT DER AUSGABEWERTE
!       -----------------------------------------------------------

            IF ((Z.EQ.7).AND.(I.EQ.3)) THEN
              IQUELDR=2
            ELSE
              IQUELDR=1
            ENDIF

            IF (((J1.NE.0).AND.(J2.NE.0)).AND.(J3.NE.0)) THEN
              IF ((LEVOB(J1-MMA+1)).LT.(INACH(J2-MMB+1))) THEN
                WRITE (MTRA,539) L,LEVOB(J1-MMA+1),GF(J1-MMA+1),          &
                     ANGST(J1-MMA+1),DRGF(IBGMARK(J3)),                   &
                     DRWELL(IBGMARK(J3)),IQUELDR
                WRITE (25,539) L,LEVOB(J1-MMA+1),GF(J1-MMA+1),            &
                     ANGST(J1-MMA+1),DRGF(IBGMARK(J3)),                   &
                     DRWELL(IBGMARK(J3)),IQUELDR
                MMA=MMA-1
                MMC=MMC-1
              ELSE IF (LEVOB(J1-MMA+1).EQ.INACH(J2-MMB+1)) THEN
                WRITE (MTRA,536) L,LEVOB(J1-MMA+1),GF(J1-MMA+1),          &
                     ANGST(J1-MMA+1),OMEG(J2-MMB+1),SLOP(J2-MMB+1),       &
                     IQUELLC(J2-MMB+1),DRGF(IBGMARK(J3)),                 &
                     DRWELL(IBGMARK(J3)),IQUELDR
                WRITE (25,536) L,LEVOB(J1-MMA+1),GF(J1-MMA+1),            &
                     ANGST(J1-MMA+1),OMEG(J2-MMB+1),SLOP(J2-MMB+1),       &
                     IQUELLC(J2-MMB+1),DRGF(IBGMARK(J3)),                 &
                     DRWELL(IBGMARK(J3)),IQUELDR
                MMA=MMA-1
                MMB=MMB-1
                MMC=MMC-1
              ELSE IF (LEVOB(J1-MMA+1).GT.INACH(J2-MMB+1)) THEN
                WRITE (25,540) L,INACH(J2-MMB+1),OMEG(J2-MMB+1),          &
                     SLOP(J2-MMB+1),IQUELLC(J2-MMB+1),DRGF(IBGMARK(J3)),  &
                     DRWELL(IBGMARK(J3)),IQUELDR
                MMB=MMB-1
                MMC=MMC-1
              ENDIF

            ELSE IF (((J1.NE.0).AND.(J2.NE.0)).AND.(J3.EQ.0)) THEN
              IF ((LEVOB(J1-MMA+1)).LT.(INACH(J2-MMB+1))) THEN
                WRITE (MTRA,539) L,LEVOB(J1-MMA+1),GF(J1-MMA+1),          &
                     ANGST(J1-MMA+1)
                WRITE (25,539) L,LEVOB(J1-MMA+1),GF(J1-MMA+1),            &
                     ANGST(J1-MMA+1)
                MMA=MMA-1
              ELSE IF (LEVOB(J1-MMA+1).EQ.INACH(J2-MMB+1)) THEN
                WRITE (MTRA,536) L,LEVOB(J1-MMA+1),GF(J1-MMA+1),          &
                     ANGST(J1-MMA+1),OMEG(J2-MMB+1),SLOP(J2-MMB+1),       &
                     IQUELLC(J2-MMB+1)
                WRITE (25,536) L,LEVOB(J1-MMA+1),GF(J1-MMA+1),            &
                     ANGST(J1-MMA+1),OMEG(J2-MMB+1),SLOP(J2-MMB+1),       &
                     IQUELLC(J2-MMB+1)
                MMA=MMA-1
                MMB=MMB-1
              ELSE IF (LEVOB(J1-MMA+1).GT.INACH(J2-MMB+1)) THEN
                WRITE (25,540) L,INACH(J2-MMB+1),OMEG(J2-MMB+1),          &
                     SLOP(J2-MMB+1),IQUELLC(J2-MMB+1)
                MMB=MMB-1
              ENDIF

            ELSE IF (((J1.NE.0).AND.(J2.EQ.0)).AND.(J3.EQ.0)) THEN
              WRITE (MTRA,536) L,LEVOB(J1-MMA+1),GF(J1-MMA+1),            &
                   ANGST(J1-MMA+1)
              WRITE (25,536) L,LEVOB(J1-MMA+1),GF(J1-MMA+1),              &
                   ANGST(J1-MMA+1)
              MMA=MMA-1

            ELSE IF (((J1.NE.0).AND.(J2.EQ.0)).AND.(J3.NE.0)) THEN
              WRITE (MTRA,539) L,LEVOB(J1-MMA+1),GF(J1-MMA+1),            &
                   ANGST(J1-MMA+1),DRGF(IBGMARK(J3)),DRWELL(IBGMARK(J3)), &
                   IQUELDR
              WRITE (25,539) L,LEVOB(J1-MMA+1),GF(J1-MMA+1),              &
                   ANGST(J1-MMA+1),DRGF(IBGMARK(J3)),DRWELL(IBGMARK(J3)), &
                   IQUELDR
              MMA=MMA-1
              MMC=MMC-1

            ELSE IF(((J1.EQ.0).AND.(J2.NE.0)).AND.(J3.NE.0)) THEN
              WRITE (25,540) L,INACH(J2-MMB+1),OMEG(J2-MMB+1),            &
                   SLOP(J2-MMB+1),IQUELLC(J2-MMB+1),DRGF(IBGMARK(J3)),    &
                   DRWELL(IBGMARK(J3)),IQUELDR
              MMB=MMB-1
              MMC=MMC-1

            ELSE IF(((J1.EQ.0).AND.(J2.NE.0)).AND.(J3.EQ.0)) THEN
              WRITE (25,540) L,INACH(J2-MMB+1),OMEG(J2-MMB+1),            &
                   SLOP(J2-MMB+1),IQUELLC(J2-MMB+1)
              MMB=MMB-1

            ELSE IF(((J1.EQ.0).AND.(J2.EQ.0)).AND.(J3.NE.0)) THEN
              WRITE (25,542) L,DRGF(IBGMARK(J3)),DRWELL(IBGMARK(J3)),     &
                   IQUELDR
              MMC=MMC-1

            ENDIF

536         FORMAT(4X,I2,3X,I2,3X,E8.2,2X,F12.3,4X,G7.2,2X,E8.2,9X,I2,5X, &
                 E9.3,2X,F9.1,10X,I2)
539         FORMAT(4X,I2,3X,I2,3X,E8.2,2X,F12.3,37X,E9.3,2X,F9.1,10X,I2)
540         FORMAT(4X,I2,3X,I2,29X,G7.2,2X,E8.2,9X,I2,5X,E9.3,2X,F9.1,    &
                 10X,I2)
542         FORMAT(4X,I2,67X,E9.3,2X,F9.1,10X,I2)

            KK1=0

570         IF(((MMA.NE.0).AND.(MMB.NE.0)).AND.(MMC.NE.0)) THEN
              KK1=KK1+1
              IF ((LEVOB(J1-MMA+1)).LT.(INACH(J2-MMB+1))) THEN
                WRITE (MTRA,645) L,LEVOB(J1-MMA+1),GF(J1-MMA+1),          &
                     ANGST(J1-MMA+1),DRGF(IBGMARK(J3)+KK1),               &
                     DRWELL(IBGMARK(J3)+KK1),IQUELDR
                WRITE (25,545) LEVOB(J1-MMA+1),GF(J1-MMA+1),              &
                     ANGST(J1-MMA+1),DRGF(IBGMARK(J3)+KK1),               &
                     DRWELL(IBGMARK(J3)+KK1),IQUELDR
                MMA=MMA-1
                MMC=MMC-1
              ELSE IF (LEVOB(J1-MMA+1).EQ.INACH(J2-MMB+1)) THEN
                WRITE (MTRA,644) L,LEVOB(J1-MMA+1),GF(J1-MMA+1),          &
                     ANGST(J1-MMA+1),OMEG(J2-MMB+1),SLOP(J2-MMB+1),       &
                     IQUELLC(J2-MMB+1),DRGF(IBGMARK(J3)+KK1),             &
                     DRWELL(IBGMARK(J3)+KK1),IQUELDR
                WRITE (25,544) LEVOB(J1-MMA+1),GF(J1-MMA+1),              &
                     ANGST(J1-MMA+1),OMEG(J2-MMB+1),SLOP(J2-MMB+1),       &
                     IQUELLC(J2-MMB+1),DRGF(IBGMARK(J3)+KK1),             &
                     DRWELL(IBGMARK(J3)+KK1),IQUELDR
                MMA=MMA-1
                MMB=MMB-1
                MMC=MMC-1
              ELSE IF (LEVOB(J1-MMA+1).GT.INACH(J2-MMB+1)) THEN
                WRITE (25,566) INACH(J2-MMB+1),OMEG(J2-MMB+1),            &
                     SLOP(J2-MMB+1),IQUELLC(J2-MMB+1),                    &
                     DRGF(IBGMARK(J3)+KK1),DRWELL(IBGMARK(J3)+KK1),       &
                     IQUELDR
                MMB=MMB-1
                MMC=MMC-1
              ENDIF
            ENDIF

            IF(((MMA.NE.0).AND.(MMB.NE.0)).AND.(MMC.NE.0)) THEN
              GOTO 570
            ENDIF

575         IF (((MMA.NE.0).AND.(MMB.NE.0)).AND.(MMC.EQ.0)) THEN
              IF ((LEVOB(J1-MMA+1)).LT.(INACH(J2-MMB+1))) THEN
                WRITE (MTRA,645) L,LEVOB(J1-MMA+1),GF(J1-MMA+1),          &
                     ANGST(J1-MMA+1)
                WRITE (25,545) LEVOB(J1-MMA+1),GF(J1-MMA+1),              &
                     ANGST(J1-MMA+1)
                MMA=MMA-1
              ELSE IF (LEVOB(J1-MMA+1).EQ.INACH(J2-MMB+1)) THEN
                WRITE (MTRA,644) L,LEVOB(J1-MMA+1),GF(J1-MMA+1),          &
                     ANGST(J1-MMA+1),OMEG(J2-MMB+1),SLOP(J2-MMB+1),       &
                     IQUELLC(J2-MMB+1)
                WRITE (25,544) LEVOB(J1-MMA+1),GF(J1-MMA+1),              &
                     ANGST(J1-MMA+1),OMEG(J2-MMB+1),SLOP(J2-MMB+1),       &
                     IQUELLC(J2-MMB+1)
                MMA=MMA-1
                MMB=MMB-1
              ELSE IF (LEVOB(J1-MMA+1).GT.INACH(J2-MMB+1)) THEN
                WRITE (25,566) INACH(J2-MMB+1),OMEG(J2-MMB+1),            &
                     SLOP(J2-MMB+1),IQUELLC(J2-MMB+1)
                MMB=MMB-1
              ENDIF
            ENDIF

            IF (((MMA.NE.0).AND.(MMB.NE.0)).AND.(MMC.EQ.0)) THEN
              GOTO 575
            ENDIF

580         IF (((MMA.NE.0).AND.(MMB.EQ.0)).AND.(MMC.NE.0)) THEN
              KK1=KK1+1
              WRITE (MTRA,645) L,LEVOB(J1-MMA+1),GF(J1-MMA+1),            &
                   ANGST(J1-MMA+1),DRGF(IBGMARK(J3)+KK1),                 &
                   DRWELL(IBGMARK(J3)+KK1),IQUELDR
              WRITE (25,545) LEVOB(J1-MMA+1),GF(J1-MMA+1),                &
                   ANGST(J1-MMA+1),DRGF(IBGMARK(J3)+KK1),                 &
                   DRWELL(IBGMARK(J3)+KK1),IQUELDR
              MMA=MMA-1
              MMC=MMC-1
            ENDIF

            IF (((MMA.NE.0).AND.(MMB.EQ.0)).AND.(MMC.NE.0)) THEN
              GOTO 580
            ENDIF

585         IF (((MMA.EQ.0).AND.(MMB.NE.0)).AND.(MMC.NE.0)) THEN
              KK1=KK1+1
              WRITE (25,566) INACH(J2-MMB+1),OMEG(J2-MMB+1),              &
                   SLOP(J2-MMB+1),IQUELLC(J2-MMB+1),DRGF(IBGMARK(J3)+KK1),&
                   DRWELL(IBGMARK(J3)+KK1),IQUELDR
              MMB=MMB-1
              MMC=MMC-1
            ENDIF

            IF (((MMA.EQ.0).AND.(MMB.NE.0)).AND.(MMC.NE.0)) THEN
              GOTO 585
            ENDIF

590         IF (((MMA.NE.0).AND.(MMB.EQ.0)).AND.(MMC.EQ.0)) THEN
              WRITE (MTRA,644) L,LEVOB(J1-MMA+1),GF(J1-MMA+1),            &
                   ANGST(J1-MMA+1)
              WRITE (25,544) LEVOB(J1-MMA+1),GF(J1-MMA+1),                &
                   ANGST(J1-MMA+1)
              MMA=MMA-1
            ENDIF

            IF (((MMA.NE.0).AND.(MMB.EQ.0)).AND.(MMC.EQ.0)) THEN
              GOTO 590
            ENDIF

595         IF (((MMA.EQ.0).AND.(MMB.NE.0)).AND.(MMC.EQ.0)) THEN
              WRITE (25,566) INACH(J2-MMB+1),OMEG(J2-MMB+1),              &
                   SLOP(J2-MMB+1),IQUELLC(J2-MMB+1)
              MMB=MMB-1
            ENDIF

            IF (((MMA.EQ.0).AND.(MMB.NE.0)).AND.(MMC.EQ.0)) THEN
              GOTO 595
            ENDIF

598         IF (((MMA.EQ.0).AND.(MMB.EQ.0)).AND.(MMC.NE.0)) THEN
              KK1=KK1+1
              WRITE(25,561) DRGF(IBGMARK(J3)+KK1),                        &
                   DRWELL(IBGMARK(J3)+KK1),IQUELDR
              MMC=MMC-1
            ENDIF

            IF (((MMA.EQ.0).AND.(MMB.EQ.0)).AND.(MMC.NE.0)) THEN
              GOTO 598
            ENDIF

! ------DIES IST DAS ENDE DER DATENMAECHTIGKEITSABFRAGE

!       FORMAT FOR HTML CATALOG
544         FORMAT(9X,I2,3X,E8.2,2X,F12.3,4X,G7.2,2X,E8.2,9X,I2,5X,        &
                 E9.3,2X,F9.1,10X,I2)
545         FORMAT(9X,I2,3X,E8.2,2X,F12.3,37X,E9.3,2X,F9.1,10X,I2)
561         FORMAT(73X,E9.3,2X,F9.1,10X,I2)
566         FORMAT(9X,I2,29X,G7.2,2X,E8.2,9X,I2,5X,E9.3,2X,F9.1,10X,I2)

!       FORMAT FOR TRANFILE
644         FORMAT(4X,I2,3X,I2,3X,E8.2,2X,F12.3,4X,G7.2,2X,E8.2,9X,I2,5X,  &
                 E9.3,2X,F9.1,10X,I2)
645         FORMAT(4X,I2,3X,I2,3X,E8.2,2X,F12.3,37X,E9.3,2X,F9.1,10X,I2)

!------------------------------------------------------------------------

            WRITE (25,'()')
          ENDDO

          WRITE (25,'()')
          CALL HTMLFOOT(25,COLDEF,COLDES,DATEION(N))
          CLOSE (25)
          CLOSE (MTRA)

        ENDDO
      ENDIF

      CLOSE (MGEN)
      CLOSE (MCOL)
      CLOSE (MLEV)
      CLOSE (MPHO)
      CLOSE (MNL2)
      CLOSE (MGF3)
      CLOSE (MNLD)
      CLOSE (MLIS)

      WRITE(*,'("DONE.")')

    CONTAINS

!----------------------------------------------------------------------
!----------------------------------------------------------------------
    SUBROUTINE COLREAD (IVON,INACH,SLOP,OMEG,IQUELLC)
!     READ COLLISIONAL DATA (atomcol)
      INTEGER(i4b), DIMENSION(NLDR)  :: IVON,INACH,IQUELLC
      REAL(dp),     DIMENSION(NLDR)  :: SLOP,OMEG

221   READ   (MCOL,218) IREC,NIZ
218   FORMAT (2(2X,I4))
      IF ((IREC.EQ.N).AND.(NIZ.NE.IZ)) THEN
        PRINT*,'INCONSISTENCY IN GENERALINFO'
        PRINT*,'FOR RECORD ',N
        PRINT*,'CHECK atomcol DATA ENTRY'
        STOP
      ENDIF
      READ   (MCOL,219) (IVON(J),INACH(J),SLOP(J),OMEG(J),              &
                         IQUELLC(J), J=1,NIZ)
219   FORMAT (2(2X,I4),2X,E8.2,2X,E8.2,2X,I4)
      IF (IREC.EQ.N) RETURN
      GOTO 221
      RETURN
    END SUBROUTINE COLREAD
!----------------------------------------------------------------------
!----------------------------------------------------------------------
    SUBROUTINE ATREAD (IND,ATOM)
!     READ LEVEL INFORMATION (atomlev)
      INTEGER(i4b)                    :: IREC,NLEVAT,J
      INTEGER(i4b), DIMENSION(NL3ARR) :: IND
      CHARACTER,    DIMENSION(NLDR)   :: ATOM*12

221   READ   (MLEV,222) IREC,NLEVAT
222   FORMAT (2(2X,I4))
      IF ((IREC.EQ.N).AND.(NLEVAT.NE.LEVAT)) THEN
        PRINT*,'INCONSISTENCY IN GENERALINFO'
        PRINT*,'FOR RECORD ',N
        PRINT*,'CHECK atomlev DATA ENTRY'
        STOP
      ENDIF
      READ   (MLEV,220) (IND(J),ATOM(J), J=1,NLEVAT)
220   FORMAT (I10,2X,A12)
      IF (IREC.EQ.N) RETURN
      GOTO 221
223   RETURN
    END SUBROUTINE ATREAD
!----------------------------------------------------------------------
!----------------------------------------------------------------------
    SUBROUTINE PHOREAD (LLOW,LUP,ALPH,BET,S,IQUELLP)
!     PHOTOIONISATION DATA (atompho)
      INTEGER(i4b), DIMENSION(NLDRL)  :: LLOW,LUP,IQUELLP
      REAL(dp),     DIMENSION(NLDRL)  :: ALPH,BET,S

221   READ   (MPHO,214) IREC,LLEVPHI
214   FORMAT (2(2X,I4))
      IF ((IREC.EQ.N).AND.(LLEVPHI.NE.LEVPHI)) THEN
        PRINT*,'INCONSISTENCY IN GENERALINFO'
        PRINT*,'FOR RECORD ',N
        PRINT*,'CHECK atompho DATA ENTRY'
        STOP
      ENDIF
      READ   (MPHO,*) (LLOW(J),LUP(J),ALPH(J),BET(J),S(J),IQUELLP(J),   &
                       J=1,LLEVPHI)
      IF (IREC.EQ.N) RETURN
      GOTO 221
      RETURN
    END SUBROUTINE PHOREAD
!----------------------------------------------------------------------
!----------------------------------------------------------------------
    SUBROUTINE READXL (XLEVEL,MAINQ,JSTATW)
!     ENERGIES AND STATISTICAL WEIGHTS (atomnl2)
      INTEGER(i4b), DIMENSION(NLDRL)  :: MAINQ,JSTATW
      REAL(dp),     DIMENSION(NLDRL)  :: XLEVEL

221   READ   (MNL2,242) IREC,NLL
242   FORMAT (2(2x,I4))
      IF ((IREC.EQ.NREC).AND.(NLL.NE.NELEV(NREC))) THEN
        PRINT*,'INCONSISTENCY IN GENERALINFO'
        PRINT*,'FOR RECORD ',NREC
        PRINT*,'CHECK atomnl2 DATA ENTRY'
        STOP
      ENDIF
      READ   (MNL2,243,END=244) (XLEVEL(J),MAINQ(J),JSTATW(J),J=1,NLL)
243   FORMAT (F12.1,2X,I2,2X,I4)
      IF (IREC.EQ.NREC) RETURN
      GOTO 221
244   RETURN
    END SUBROUTINE READXL
!----------------------------------------------------------------------
!----------------------------------------------------------------------
    SUBROUTINE RLINE (ANGST,GF,IUOLEV)
!     DATA FROM SMALL LINE LIST (atomgf3)
      INTEGER(i4b), DIMENSION(NL3ARR)  :: IUOLEV
      REAL(dp),     DIMENSION(NL3ARR)  :: ANGST,GF

221   READ   (MGF3,216) IREC,LLINZAHL
216   FORMAT (2(2X,I4))
      IF ((IREC.EQ.N).AND.(LLINZAHL.NE.LINZAHL)) THEN
        PRINT*,'INCONSISTENCY IN GENERALINFO'
        PRINT*,'FOR RECORD ',N
        PRINT*,'CHECK atomgf3 DATA ENTRY'
        STOP
      ENDIF
      READ   (MGF3,217) (IUOLEV(J),ANGST(J),GF(J), J=1,LLINZAHL)
217   FORMAT (2X,I8,2X,F12.3,2X,E8.2)
      IF (IREC.EQ.N) RETURN
      GOTO 221
      RETURN
    END SUBROUTINE RLINE
!----------------------------------------------------------------------
!----------------------------------------------------------------------
    SUBROUTINE RLINEDR (INLEV,IBGMARK,DRWELL,DRGF)
!     DIELECTRONIC RECOMBINATION DATA (atomnldr)
      INTEGER(i4b), DIMENSION(NLDRL)  :: INLEV,IBGMARK
      REAL(dp),     DIMENSION(NLDR)   :: DRWELL,DRGF

221   READ   (MNLD,227) IREC,LLEVDR,LLDADR
227   FORMAT (3(2x,I4))
      IF ((IREC.EQ.N).AND.((LLEVDR.NE.LEVDR).OR.(LLDADR.NE.LDADR))) THEN
        PRINT*,'INCONSISTENCY IN GENERALINFO'
        PRINT*,'FOR RECORD ',N
        PRINT*,'CHECK atomnldr DATA ENTRY'
        STOP
      ENDIF
      JCOUNT = 0
      DO J = 1,LLEVDR
        READ   (MNLD,271) INLEV(J),NDR
271     FORMAT (2(2X,I4))
        READ   (MNLD,272) (DRWELL(JJ),DRGF(JJ), JJ=JCOUNT+1,JCOUNT+NDR)
272     FORMAT (2X,F10.1,2X,E10.3)
        IBGMARK(J) = JCOUNT + 1
        JCOUNT     = JCOUNT + NDR
      ENDDO
      IF (IREC.EQ.N) RETURN
      GOTO 221
      RETURN
    END SUBROUTINE RLINEDR
!----------------------------------------------------------------------
!----------------------------------------------------------------------
    SUBROUTINE PHINIT (ISUM,LHUP,ZQUANT)
!     PHOTOIONISATION DATA
      INTEGER, PARAMETER                       :: NQ=10
      INTEGER(i4b)                             :: MM,NA,M1,M2,N1,N2
      LOGICAL,      DIMENSION(NII)             :: LOGMERK
      INTEGER(i4b), DIMENSION(NLDRL)           :: LHUP
      REAL(dp),     DIMENSION(NLDRL)           :: ZQUANT
      INTEGER(i4b), DIMENSION(NLDRL,NKIS)      :: ISUM

      DO MM=1,LEVAT
        READ (ATOM(MM)(4:5),'(I2)') NZQ
        ZQUANT(MM) = REAL(NZQ,dp)
        READ (ATOM(MM)(10:11),'(I2)') LHUP(MM)
      ENDDO

      DO MM=1,LEVAT
        LOGMERK(MM) = .FALSE.
      ENDDO
      DO NA=1,NQ
        DO MM=1,LEVAT
          ISUM(MM,NA) = 0
        ENDDO
      ENDDO

      DO MM=1,LEVAT
        IF (.NOT. LOGMERK(MM)) THEN
          M1 = MAINQ(IND(MM))
          M2 = LHUP(MM)
          DO NA=MM,LEVAT
            IF (M1.GT.NQ) STOP 'MAIN QUANTUM NO. GT.10!'
            N1 = MAINQ(IND(NA))
            N2 = LHUP(NA)
            IF (.NOT. LOGMERK(NA) .AND. N1.EQ.M1 .AND.N2.EQ.M2) THEN
              ISUM(N2,N1) = JSTATW(IND(NA)) + ISUM(N2,N1)
              LOGMERK(NA) = .TRUE.
            ENDIF
          ENDDO
        ENDIF
      ENDDO
      RETURN
    END SUBROUTINE PHINIT
!----------------------------------------------------------------------
!----------------------------------------------------------------------
    SUBROUTINE PHOTOCS (ALPHA,BETA,SA,LEVX,FREQ)
!     PARAMETER ALPHA (THRESHOLD) OF THE SEATON FORMULA
      LOGICAL                    :: ALPHAL
      INTEGER(i4b)               :: NZP,NZPL,I1,I2,I3,MNI,LEVX
      REAL(dp)                   :: ALPHATH,F,ALPHA,BETA,SA,FREQ
      DATA NZP,NZPL /1,1/

!-----STATEMENT FUNCTION FOR APPROXIMATION OF ALPHATH
      ALPHATH (I1,I2,I3,F) = HTCSSBF*REAL(I1*I2,dp)/(F*REAL(2*I3**3,dp))

      CALL BETASFR (ALPHAL,BETA,SA,LEVX,FREQ,NZP,NZPL)

      MNI = MAINQ(L)

      IF (ALPHAL) THEN
        ALPHA = ALPH(NZP)
        IF (NZP.EQ.LEVPHI) THEN
          NZP = 1
        ELSE
          NZP = NZP + 1
        ENDIF
      ENDIF


      IF (LEVAT.NE.0 .AND. L.EQ.IND(NZPL)) THEN
!****** THE STATISTICAL FACTOR IN GOULD'S FORMULA IS USED
!****** CORRECTLY.
        IF (.NOT.ALPHAL) THEN
          ALPHA = HTCSSBF / FREQ * ZQUANT(NZPL) / REAL(MNI,dp)
!*Note*
! the statistical factor below underestimates the
! cross section for excited levels
!         ALPHA = ALPHA * REAL(JSTATW(IND(NZPL)),dp) /                  &
!                         REAL(ISUM(LHUP(NZPL),MNI),dp)
          NZQUA = INT(ZQUANT(NZPL))
        ENDIF
        IF (NZPL.EQ.LEVAT) THEN
          NZPL = 1
        ELSE
          NZPL = NZPL + 1
        ENDIF
      ELSE
!****** THE STATISTICAL FACTOR IN GOULD'S FORMULA IS USED
!****** IN A HYDROGENIC APPROXIMATION.
        IF (.NOT.ALPHAL) THEN
          NEL = NELF (L,N,MNI)
          IF (NEL <= 0) THEN
            PRINT*,' ELEMENT ',K,' IN IONIS.STAGE ',IS,' NOEL=',NEL
            STOP ' PHOTOCS: NOEL <= 0'
          END IF
          ALPHA = ALPHATH (JSTATW(L),NEL,MNI,FREQ)
        ENDIF
      ENDIF
      RETURN
    END SUBROUTINE PHOTOCS
!----------------------------------------------------------------------
!----------------------------------------------------------------------
    SUBROUTINE BETASFR (ALPHAL,BETA,SA,LEVX,FREQ,NZP,NZPL)
!     SUPPLY-SUBROUTINE FOR CALCULATION OF THE PARAMETERS BETA AND S
!     IN SEATON'S FORMULA. LEVX IS THE LEVEL OF THE HIGHER ION WHERE
!     IONISATION TAKES PLACE, FREQ IS THE ACCORDING FREQUENCY
!     CALLED BY 'PHOTOCS' - 3.91 ADI
      LOGICAL      :: ALPHAL
      INTEGER(i4b) :: NZP,NZPL,LEVX
      REAL(dp)     :: BETA,SA,FREQ

!     ALPHAL BECOMES TRUE, IF THERE ARE EXTERNAL PHOTODATA FOR LEVEL L
      ALPHAL=.FALSE.


      IF (LEVPHI.NE.0) THEN
        IF (L.EQ.LLOW(NZP)) THEN
          ALPHAL = .TRUE.
          BETA   = BET(NZP)
          SA     = S(NZP)
          IF (I.NE.0) THEN
            LEVX = LUP(NZP)
          ELSE
            LEVX=1
          ENDIF
        ENDIF
      ENDIF

      IF (.NOT.ALPHAL) THEN
        BETA=2._dp
        SA=2._dp
        IF (LEVAT.NE.0 .AND. L.EQ.IND(NZPL) .AND. LHUP(NZPL).NE.0) THEN
          LEVX=LHUP(NZPL)
        ELSE
          LEVX=1
        ENDIF
      ENDIF
      IF (SA.LT.1._dp) THEN
        SA=1._dp
      ELSEIF (SA.GT.5._dp) THEN
        SA=5._dp
      ENDIF
!       SA=REAL(NINT(2.*SA),dp)/2.

!     ZUORDNUNG DES UEBERGANGSLEVELS---------------
      IF (ENERGY(Z,I+1).EQ.0._dp) XLEVEL1(LEVX) = 0._dp
      FREQ = ENERGY(Z,I) - XLEVEL(L) + XLEVEL1(LEVX)
!-------------------------------------------------
      IF (FREQ.LE.0._dp) THEN
        PRINT*,'FREQ IN BFRCF LESS OR EQUAL 0!'
        PRINT*,'Z=',Z,' I=',I,' L=',L
        PRINT*,'ENERGY(Z,I)=',ENERGY(Z,I),' XLEVEL(L)=',XLEVEL(L)
        STOP 'FREQ<=0'
      ENDIF
      RETURN
    END SUBROUTINE BETASFR
!----------------------------------------------------------------------
!----------------------------------------------------------------------
    REAL(dp) FUNCTION NELF (L,N,MNI)
!     CALCULATES THE NUMBER OF ELECTRONS OF ION N
      INTEGER(i4b)  :: L,N,MNI

      IF (INOEL(N).GE.100) THEN
        IF (MNI.EQ.3) THEN
          NELF=NOEL(N)
        ELSE
          NELF=INOEL(N)/100
        ENDIF
      ELSE
        IF (L.GT.INOEL(N)) THEN
          NELF=1
        ELSE
          NELF=NOEL(N)
        ENDIF
      ENDIF
      RETURN
    END FUNCTION NELF
!----------------------------------------------------------------------
!----------------------------------------------------------------------
    SUBROUTINE COLORR

      COLDEF(1)='#88FF88'
      COLDES(1)='excellent'
      COLDEF(2)='#DDFF88'
      COLDES(2)='good'
      COLDEF(3)='#FFCC88'
      COLDES(3)='poor'
      COLDEF(4)='#FF8888'
      COLDES(4)='bad'
      RETURN
    END SUBROUTINE COLORR
!----------------------------------------------------------------------
!----------------------------------------------------------------------

    END PROGRAM MKATOM

!----------------------------------------------------------------------
!----------------------------------------------------------------------
    SUBROUTINE HTMLHEAD(NOUT,TEXT)

      USE M_type_definitions, ONLY: i4b
      INTEGER(i4b) :: NOUT,LEN_T
      CHARACTER    :: TEXT*260

      LEN_T = LEN_TRIM(TEXT)
      WRITE (NOUT,'("<!DOCTYPE HTML PUBLIC ""-//W3C//DTD HTML 4.01 Transitional//EN"">")')
      WRITE (NOUT,'("<HTML><HEAD>")')
      WRITE (NOUT,'("<META http-equiv=""Content-Type"" content=""text/html; charset=US-ASCII"">")')
      WRITE (NOUT,'("<TITLE>",A,"</TITLE>")') TEXT(1:LEN_T)
      WRITE (NOUT,'("</HEAD><BODY BGCOLOR=""#FFF8DC"" LINK=""BLUE"" ALINK=""MAGENTA"" VLINK=""RED"">")')
      IF (NOUT.EQ.110) THEN
        WRITE(NOUT,'("<TABLE BORDER=0>")')
        WRITE(NOUT,'("<TR><TD WIDTH=25></TD><TD></TD>")')
        WRITE(NOUT,'("  <TD WIDTH=55 ALIGN=""CENTER"" BGCOLOR=""SILVER"">I</TD>")')
        WRITE(NOUT,'("  <TD WIDTH=55 ALIGN=""CENTER"" BGCOLOR=""SILVER"">II</TD>")')
        WRITE(NOUT,'("  <TD WIDTH=55 ALIGN=""CENTER"" BGCOLOR=""SILVER"">III</TD>")')
        WRITE(NOUT,'("  <TD WIDTH=55 ALIGN=""CENTER"" BGCOLOR=""SILVER"">IV</TD>")')
        WRITE(NOUT,'("  <TD WIDTH=55 ALIGN=""CENTER"" BGCOLOR=""SILVER"">V</TD>")')
        WRITE(NOUT,'("  <TD WIDTH=55 ALIGN=""CENTER"" BGCOLOR=""SILVER"">VI</TD>")')
        WRITE(NOUT,'("  <TD WIDTH=55 ALIGN=""CENTER"" BGCOLOR=""SILVER"">VII</TD>")')
        WRITE(NOUT,'("  <TD WIDTH=55 ALIGN=""CENTER"" BGCOLOR=""SILVER"">VIII</TD>")')
        WRITE(NOUT,'("</TR><TR><TD COLSPAN=10></TD></TR>")')
      ENDIF
      RETURN

    END SUBROUTINE HTMLHEAD
!----------------------------------------------------------------------
!----------------------------------------------------------------------
    SUBROUTINE HTMLLINK(NOUT,CION,CIOT,COLDEF)

      USE M_type_definitions, ONLY: i4b
      INTEGER(i4b) :: NOUT,LEN1
      CHARACTER    :: CION*7,CIOT*7,COLDEF*7

      LEN1   = LEN_TRIM(CION)
      WRITE  (NOUT,900) COLDEF,CION(1:LEN1),CIOT(1:LEN1)
900   FORMAT ('  <TD ALIGN="CENTER" BGCOLOR="',A,'"><A HREF="cat/',A,   &
              '.html">',A,'</A></TD>')
      RETURN
    END SUBROUTINE HTMLLINK
!----------------------------------------------------------------------
!----------------------------------------------------------------------
    SUBROUTINE HTMLFOOT(NOUT,COLDEF,COLDES,DATE,HTMLANG)

      USE M_type_definitions, ONLY: i4b,dp
      IMPLICIT NONE
      INTEGER(i4b)           :: NOUT,LEN1,JI
      REAL(dp), DIMENSION(5) :: HTMLANG
      CHARACTER              :: COLDEF(4)*7,COLDES(4)*10,DATE*8,        &
                                LAMBDA*(*),LOEQ*(*),LOTH*(*),LINK1*(*), &
                                LINK2*(*)

      PARAMETER(LAMBDA='&lambda;</A>')
      PARAMETER(LOEQ='  &nbsp;&le;&nbsp;')
      PARAMETER(LOTH='&nbsp;&Aring;&nbsp;&lt;&nbsp;')
      PARAMETER(LINK1='  <A HREF="data/lines')
      PARAMETER(LINK2='.nl3">')

      IF (NOUT.EQ.110) THEN
        WRITE (NOUT,'("<TR><TD COLSPAN=10></TD></TR>")')
        WRITE (NOUT,'("<TR><TD COLSPAN=10></TD></TR>")')
        WRITE (NOUT,'("<TR><TD></TD><TD></TD>")')
        WRITE (NOUT,'("<TD COLSPAN=8>&nbsp;atomic data status:</TD></TR>")')
        WRITE (NOUT,'("<TR><TD></TD><TD></TD>")')
        DO JI = 1, 4
          LEN1 = LEN_TRIM(COLDES(JI))
          WRITE (NOUT,'("  <TD COLSPAN=2 BGCOLOR=""",A,""" ALIGN=""CENTER"">", &
                        A,"</TD>")') COLDEF(JI),COLDES(JI)(1:LEN1)
        ENDDO
        WRITE (NOUT,'("</TR>")')
        WRITE (NOUT,'("<TR><TD COLSPAN=10></TD></TR>")')
        WRITE (NOUT,'("<TR><TD COLSPAN=10></TD></TR>")')
        WRITE (NOUT,'("<TR><TD></TD><TD></TD>")')
        WRITE (NOUT,'("<TD COLSPAN=8>&nbsp;line list for all atoms, ",  &
                      "sorted by wavelength:</TD></TR>")')
        WRITE (NOUT,'("<TR><TD></TD><TD></TD>")')
        WRITE (NOUT,'("<TD NOWRAP COLSPAN=8 BGCOLOR=""#EDEED9"" ALIGN=""CENTER"">")')
        WRITE (NOUT,'(A,I1,A,A)') LINK1,1,LINK2,LAMBDA
        DO JI=2,5
          WRITE (NOUT,'(A,I4,A)') LOEQ,INT(HTMLANG(JI)),LOTH
          WRITE (NOUT,'(A,I1,A,A)') LINK1,JI,LINK2,LAMBDA
        ENDDO
        WRITE (NOUT,'("</TD></TR>")')
        WRITE (NOUT,'("</TABLE>")')
      ELSE
        WRITE (NOUT,'("</PRE>")')
      ENDIF
      WRITE (NOUT,'("<BR><FONT COLOR=""SILVER""><SMALL>",               &
                    "&copy; USM Hot Star Group, ",A4,"-",A2,"-",A2)')   &
                    DATE(1:4),DATE(5:6),DATE(7:8)
      WRITE (NOUT,'("</SMALL></FONT></BODY></HTML>")')
      RETURN

    END SUBROUTINE HTMLFOOT
