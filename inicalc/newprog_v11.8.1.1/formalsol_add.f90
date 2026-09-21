Module version_formalsol_add

! additional routines for FORMALSOL, when calculating spectrum in predefined
! range (only if OPTCMF_ALL = .true.), with keyword 'UV'

! NOTES for further modifications

! thus far, Voigt profiles NOT included (only Dopplerbroadening and
! Stark-broadening for HHe considered)

! version 0.1 Aug. 2016: calculation of opacities completely changed.
!   opal_range and sline_range now contain actual
!   (and not log) values  -> beta_p in calc_xmax
!   frequential and spatial interpolations linear

! version 0.2 Sept. 2016: further improvements, particularly lamgrid_cmf
!   still linear interpolation onto radial micro-grid

! version 0.3 Sept. 2016: further improvements, particularly log-log
!   interpolation onto radial micro-grid

! version 0.4 Sept. 2016: now: frequential lin-log interpolation for opacities
!   tested: increasing precision from v0.1 to v0.4
!   tested: v0.4 vs. v0.0 gives almost identical results,
!   but much faster.

! version 0.4.1 end Oct 2016: bug in u0 fixed

! version 0.4.2 May 2017: improved treatment of inversion: opal_plus_range
!   etc.
!   now contains absorption part of opacity
!   version 0.5.0 Oct.2020: compatible with gfortran

! version 0.5.1 Nov.2022: small changes in calc_lines:
!   since LWION1 has changed to taur ge 10,
!   there are now a couple of lines
!   (between LWION and LWION1) with OPAL and SLINE = 0,
!   which lie outside imin, imax. This has been changed:
!   now, also lines outside imin, imax (approximate
!   treatment) with defined occupation numbers
!   (gs,'m','s' above LWION1) are accounted for.
!   NOTE:
!   (i) for selected elements, LINE_LIST_MERGED
!   also includes lines between levels that
!   are NOT gs,'m',or 's'. i.e. with occup=0,
!   outside imin,imax when ll lt lwion
!   (ii)for not-selected = approximate elements,
!   all lines between approximate levels
!   that are NOT gs,'m',or 's' levels have been
!   excluded from LINE_LIST_MERGED (in cmf_all/line_list).
!   Thus, for not-selected elements there should be
!   no lines with occup=0 (except, of course, those
!   with NUP=0)

! version 0.5.2 April 2025: small changes regarding gfortran

! version 1.0 Dec 2025: polished
  
! WRITTEN:
! 07/15 by JP

! HISTORY:
! completely new

Character (5) :: ver_formalsol_add = '1.0.0'
End Module

!--------------------------------------------------------------------------- &

Module formalsol_add_var

  Use :: nlte_type
  Implicit None

  Integer (i4b) :: ntot, icount_hhe

  Real (dp) :: lamb1, lamr1 !       extended edge wavelengths (vmax + vbroad)

  Integer (i4b), Dimension (:), Allocatable :: id1
  Real (sp), Dimension (:), Allocatable :: gf1, xlam1

  Integer (i4b), Dimension (:), Allocatable :: k_range, inverted_range, &
    index_hhe, l23, lstore
  Real (dp), Dimension (:, :), Allocatable :: opal_range, sline_range, &
    opal_plus_range
! note: at first, opal_part_lgrid contains absorption part of opacity.
! at the end of subr. opal_cmfgrid, opal_part_lgrid contains
! stimulated emission part of opacity, and opal_lgrid absorption part,
! whenever opalmin_lgrid(freq) < 0.
  Real (dp), Dimension (:, :), Allocatable :: opal_lgrid, opal_part_lgrid, &
    sl_lgrid

  Real (dp), Dimension (:), Allocatable :: lamgrid, lamgrid_cmf, &
    dlam_max_stark, beta_p, qstore, opalmin_lgrid

  Integer (i4b), Parameter :: natom = 30

! JO: if elements between He and C are implemented, change first 'heavy'
! element (first_hel = 6 so far)
  Integer (i4b), Parameter :: first_hel = 6

! maximum width of Doppler-broadening (in vth, used for non-Stark lines in
! prepray_range)
  Real (dp), Parameter :: xmax = 4.D0

! control parameter for Stark-broadening (see subr. calcxmax_range and
! calcxmax)
  Real (dp), Parameter :: eps_broad = 5.D-3

! maximum integration depth
  Real (dp), Parameter :: taumax = 15.D0
! neglect lines with tau_sob < tauweak in range 1,lwion-1
! to included all lines, set tauweak to very low value (1.d-10 or lower)
  Real (dp), Parameter :: tauweak = 1.D-3 ! compromise; for tests, set to 0.01
  Real (dp), Parameter :: tauout = 1.D0

! maximum width of Stark-broadening (max(vbroad) calculated in
! calcxmax_range),
! to be used in prepray_range
  Real (dp), Parameter :: vmax_broad_safety = 3000.D5 ! (for boundaries in
! calc_lines)
  Real (dp) :: vmax_broad = 0. !    in case that there are no Stark-broadened
! lines

  Character (2), Dimension (natom) :: names1
  Real (dp), Dimension (natom) :: aweight


! for info and checks
  Data names1/'H ', 'HE', 'LI', 'BE', 'B ', 'C ', 'N ', 'O ', 'F ', 'NE', &
    'NA', 'MG', 'AL', 'SI', 'P ', 'S ', 'CL', 'AR', 'K ', 'CA', 'SC', 'TI', &
    'V ', 'CR', 'MN', 'FE', 'CO', 'NI', 'CU', 'ZN'/

  Data aweight/1.008, 4.003, 6.941, 9.012, 10.811, 12.011, 14.007, 16.000, &
    18.998, 20.180, 22.990, 24.305, 26.982, 28.085, 30.974, 32.066, 35.453, &
    39.948, 39.098, 40.078, 44.956, 47.88, 50.941, 51.996, 54.938, 55.847, &
    58.933, 58.69, 63.546, 65.39/


End Module

!--------------------------------------------------------------------------- &

Subroutine calc_lines(nd, xne, te, clfac, r, velo, dvdr, sr, srvmax, lam_blu, &
  lam_red, iblu, ired, vturb, yhe, balmer_lemke, pb_lemke, maxncomp)

! calculates line opacities/source-functions for line-list provided from
! package cmf_all and corresponding occupation numbers

! OPAL corresponds to pi e2/me c * 1.d-8 *SR/vmax * lam * gf * nl/gl / clfac
! = true (mean) opacity * lambda * SR / vmax
! calculated opacity for prepray_range =  ln(true opac
! *SR*lambda/sqrt(pi)/clight)
! =  ln (OPAL* vmax/clight /sqpi)

  Use :: nlte_type
  Use :: nlte_dim, Only: id_ndept, id_llevs, id_frec1
  Use :: fund_const, Only: clight, hc2, sqpi
  Use :: princesa_var, Only: gl, le, labat, labl
  Use :: formalsol_var, Only: modnam => file, nstark, nlevl, nlevu
  Use :: formalsol_add_var, Only: ntot, id1, lamb1, lamr1, xlam1, gf1, names1, &
    opal_range, opal_plus_range, sline_range, k_range, inverted_range, &
    index_hhe, icount_hhe, vmax_broad_safety, tauweak, l23, beta_p

  Implicit None

  Integer (i4b), Parameter :: nrecex = id_llevs
  Integer (i4b), Parameter :: nd1 = id_ndept
  Integer (i4b), Parameter :: ifretot = id_frec1

  Integer (i4b), Parameter :: natom = 30, nion = 9, nlev = 50, nrec = 149

  Real (dp), Parameter :: c_vanreg = 0.5*2.055D-23
  Real (dp), Parameter :: c_forbid = 8.63D-6/(0.5*1.3707D-7)
  Real (dp), Parameter :: gfcut = 1.D-5 ! see sumopal and opacity_cmf

  Real (dp), Parameter :: c1 = 1.4388354967334D8, c2 = 3.972970127D8
! for Planck function, identical with function bnue

  Real (dp), Parameter :: safety = vmax_broad_safety/clight ! extend interval
! to allow for
! red/blue-shifted
! lines

  Integer (i4b), Intent (In) :: nd, maxncomp

  Real (dp), Dimension (nd1), Intent (In) :: xne, clfac, r, velo, dvdr
  Real (dp), Dimension (nd1), Intent (Out) :: te

  Real (dp), Intent (In) :: sr, srvmax, lam_blu, lam_red, vturb, yhe
  Integer (i4b), Intent (Out) :: iblu, ired
  Logical, Intent (In) :: balmer_lemke, pb_lemke

  Integer (i4b), Dimension (nd1) :: lll

  Real (dp), Dimension (nd1) :: vr, tau0, collfac, opal, sline, occlow, occup, &
    xj1, xk1, mustar, eddf, u0, vanreg, vanreg1, eps, opal_plus ! , eps1

  Integer (i4b), Dimension (ifretot) :: lthin

  Real (dp), Dimension (ifretot) :: fre
  Real (dp), Dimension (nd1, ifretot) :: xj, xxk
  Real (dp), Dimension (nd1, ifretot) :: opac, scont, thomson

  Real (dp), Dimension (nrecex, nd1) :: occexpl

  Integer (i4b), Dimension (natom) :: jatom_full
  Integer (i4b), Dimension (natom, nion) :: indrec
  Integer (i4b), Dimension (natom, nd1) :: met_imin, met_imax
  Real (dp), Dimension (nrec, nlev, nd1) :: occng

  Integer (i4b), Dimension (:), Allocatable :: lineno

  Character, Dimension (:), Allocatable :: level*6, leveu*6

  Real (dp), Parameter :: c_tau = 1.D-8*0.02654, hc2_8 = hc2*1.D24 ! conversion
! to
! Angstrom
  Real (dp), Parameter :: const_opal_stark = sqpi*clight*1.D8

  Real (dp) :: c_tau1, const, lam, lam3, occlowi, occupi, ener, q, q1, &
    eddtest, taubi, taui, betacic, beta, bni, epsi, eps1i, vmax, const_opal, &
    maxtau, xlamb, gll, glu, gf, weight, tauc

  Integer (i4b) :: lwion, lwion1, ifre, i, j, ll, ii, index, k, kk, irest, ml, &
    mu, low, up, icount_expl, icount_bg, icount_bg_noup, irec, ifstart, lth, &
    n_inverted, n_weak, ncomp, nst

  Logical :: almost_converged

  If (nd/=nd1) Stop ' nd ne nd1 in calc_lines'

! read merged line list (as calculated and saved in CMF_ALL
  Open (1, File=trim(modnam)//'/LINE_LIST_MERGED', Form='UNFORMATTED')
  Read (1) ntot, lwion, lwion1, indrec
  Print *, 'lwion:', lwion, ' lwion1', lwion1

  Allocate (id1(ntot), xlam1(ntot), gf1(ntot))

  Read (1) id1, xlam1, gf1
  Close (1)


  Open (1, File=trim(modnam)//'/CONT_FORMAL_ALL', Status='unknown', &
    Form='unformatted')

! actually, strue, but not needed
  Read (1) ifre, (fre(i), i=1, ifretot), ((opac(i,j),i=1,nd), j=1, ifretot), &
    ((thomson(i,j),i=1,nd), j=1, ifretot), ((scont(i,j),i=1,nd), j=1, ifretot) &
    , ((xj(i,j),i=1,nd), j=1, ifretot), (te(i), i=1, nd)
  Close (1)

  Open (1, File=trim(modnam)//'/CONT_FORMAL_CMF', Status='unknown', &
    Form='unformatted')

! total source function (CMF)
! actually, opac_nolines (will be used to calculate tauc)
  Read (1)((scont(i,j),i=1,nd), j=1, ifretot), ((opac(i,j),i=1,nd), j=1, &
    ifretot), ((xxk(i,j),i=1,nd), j=1, ifretot), (lthin(j), j=1, ifretot)
  Close (1)

  Do kk = 1, ifre
    If (fre(kk)/=0.) Then
      fre(kk) = log10(1.D8/fre(kk))
      xj(:, kk) = log10(xj(:,kk))
      xxk(:, kk) = log10(xxk(:,kk))
    End If
  End Do

  c_tau1 = c_tau*srvmax

  vr = velo/r
  tau0 = 1.D0/(dvdr/3.+2./3.*vr) !  at mue^2 = 1/3
  collfac = c_vanreg*xne/sqrt(te)

! extend interval
  vmax = sr/srvmax
  lamb1 = lam_blu*(1.-(vmax/clight+safety))
  lamr1 = lam_red*(1.+(vmax/clight+safety))

  If (xlam1(1)>lamr1) Stop ' first line redwards from interval (calc_lines)'
  If (xlam1(ntot)<lamb1) Stop &
    ' last line bluewards from interval (calc_lines)'

  const_opal = vmax/(sqpi*clight) ! for opal_range; division by c, since
! vthray=vth/c

! define range for lines in considered range
  Do i = 1, ntot
    If (xlam1(i)>=lamb1) Exit
  End Do

  iblu = max(1, i-1)

  Do i = iblu + 1, ntot
    If (xlam1(i)>lamr1) Exit
  End Do

  ired = i
  Print *, ' total number of lines in considered range:', ired - iblu + 1

  Allocate (opal_range(nd,iblu:ired), opal_plus_range(nd,iblu:ired), &
    sline_range(nd,iblu:ired))
  Allocate (k_range(iblu:ired), inverted_range(iblu:ired))

  inverted_range = 0

! read NLTE occupation numbers (expl. elements) at all depth points
  Open (Unit=1, File=trim(modnam)//'/NLTE_POP', Status='UNKNOWN', &
    Recl=nrecex*8, Access='DIRECT')

  Do ll = 1, nd
    Read (1, Rec=ll)(occexpl(i,ll), i=1, nrecex)
    lll(ll) = ll !                  auxiliary vector
  End Do
  Close (1)

! read NLTE occupation numbers (bg elements) and considered ions
! at all depth points
  Open (1, File=trim(modnam)//'/OCCNG', Status='unknown', Form='unformatted')
  Read (1) met_imin, met_imax, occng
  Read (1) almost_converged, jatom_full
  Close (1)

  icount_expl = 0 !                 counter for lines from explicit elements
  icount_bg = 0 !                   counter for lines from background elements
  icount_bg_noup = 0 !              counter for lines from explicit elements
! without upper level
  n_inverted = 0 !                  counter for inverted lines
! (inverted_range=1)
  n_weak = 0 !                      counter for weak lines (inverted_range=3)

  icount_hhe = 0 !                  counter for H/He lines

! loop over all considered lines, to calculate opal and etal
all_lines: Do ii = iblu, ired
    index = id1(ii)
    If (index>=1D8) Then
!     explicit element
      kk = index/1D8
      irest = index - kk*1D8
      ml = irest/10000
      mu = irest - ml*10000
      If (kk/=le(ml) .Or. kk/=le(mu)) Stop ' problems with kk (calc_lines)'

      Do i = 1, natom
        If (labat(kk)==names1(i)) Then
          k_range(ii) = i
          Exit
        End If
      End Do

      lam = xlam1(ii)
!     print*,k_range(ii),ml,mu,lam,gf1(ii)
      lam3 = lam**3
!     print*,lam,ml,mu
      const = c_tau1*lam*gf1(ii)
      occlow = occexpl(ml, :)/gl(ml)
      occup = occexpl(mu, :)/gl(mu)
      If (occlow(1)*occup(1)==0.) Then
!       this should never happen (at least for GLOBAL = .true.)
!       JO: check for GLOBAL = .false.
        Print *, index, ' ', lam
        Stop ' occlow or occup eq 0 (expl. elements) in calc_lines'
      End If
!     pi e2/me c * 1.d-8 *SR/vmax * lam * gf * (nl/gl - nu/gu) /clfac
      opal = const*(occlow-occup)/clfac
      opal_plus = const*occlow/clfac
      sline = hc2_8/lam3/(occlow/occup-1.D0)

!     print*,lam,'explicit'
!     do ll=1,nd
!     print*,opal(ll),sline(ll)
!     enddo

      icount_expl = icount_expl + 1

    Else
!     background element
      k = index/1000000
      k_range(ii) = k
      irest = index - k*1000000
      j = irest/100000
      irest = irest - j*100000
      low = irest/100
      up = irest - low*100
      lam = xlam1(ii)
      lam3 = lam**3
      irec = indrec(k, j)
      const = c_tau1*lam*gf1(ii)
!     print*,lam,k,j,low,up,ii
      If (up/=0) Then

!       lower and upper level  present

        If (jatom_full(k)==1) Then
!         selected elements, division lwion1
          Do ll = 1, lwion1 - 1
!           outside
!           JP changed Nov. 2022, to include lines (approximate treatment)
!           outside
!           imin, imax
            occlowi = occng(irec, low, ll) ! occng = n/g
            occupi = occng(irec, up, ll)
!           check "exact" NLTE
            If (j>=met_imin(k,ll) .And. j<=met_imax(k,ll)) Then
              If (occlowi*occupi<=0.) Then
                Print *, index, ' ', lam
                Stop ' selected, outside: occnum = 0 for imin le j le imax'
              End If
            Else
!             approximate solution (outside imin,imax) for selected elements:
!             should happen only for ll < lwion-1
!             NOTE: for selected elements, LINE_LIST_MERGED also includes
!             lines between
!             levels that are NOT gs,'m',or 's'. i.e. with occup=0, outside
!             imin,imax
!             when ll lt lwion
              If (occlowi*occupi==0.) Then
                If (ll>=lwion) Then
                  Print *, index, ' ', lam
                  Stop ' approximate elements, ll ge lwion: occnum=0'
                End If
                opal(ll) = 0.
                opal_plus(ll) = 0.
                sline(ll) = 0.
                Cycle
              End If
            End If
!           either inside imin, imax, or beyond LWION, or lines between
!           gs,'m','s' levels
            opal(ll) = const*(occlowi-occupi)/clfac(ll)
            opal_plus(ll) = const*occlowi/clfac(ll)
            sline(ll) = hc2_8/lam3/(occlowi/occupi-1.D0)
          End Do
!         end change

          Do ll = lwion1, nd
!           inside (all levels should be populated)
            occlowi = occng(irec, low, ll) ! occng = n/g
            occupi = occng(irec, up, ll)
            If (occlowi*occupi==0) Then
              Print *, index, ' ', lam
              Stop ' occlow or occup eq 0 (inner selected elements) &
                &in calc_lines'
            End If
            opal(ll) = const*(occlowi-occupi)/clfac(ll)
            opal_plus(ll) = const*occlowi/clfac(ll)
            sline(ll) = c2/(exp(c1/lam/te(ll))-1.D0)/lam3
          End Do

!         print*,'selected with upper level'
!         print*,lam,k,j,low,up
!         do ll=1,nd
!         print*,opal(ll),sline(ll),met_imin(k,ll),met_imax(k,ll)
!         enddo

        Else
!         other elements (approximate), division lwion
!         note: in cmf_all/line_list, all lines between approximate levels
!         that are NOT gs,'m',or 's' levels have been excluded from
!         LINE_LIST_MERGED.
          Do ll = 1, lwion - 1
            occlowi = occng(irec, low, ll) ! occng = n/g
            occupi = occng(irec, up, ll)
            If (occlowi*occupi==0) Then ! change of ionization
              Print *, index, ' ', lam
              Stop &
                ' occlow or occup eq 0 (outer approx. elements) in calc_lines'
            End If
            opal(ll) = const*(occlowi-occupi)/clfac(ll)
            opal_plus(ll) = const*occlowi/clfac(ll)
            sline(ll) = hc2_8/lam3/(occlowi/occupi-1.D0)
          End Do
          Do ll = lwion, nd
!           inside
            occlowi = occng(irec, low, ll) ! occng = n/g
            occupi = occng(irec, up, ll)
            If (occlowi*occupi==0) Then
              Print *, index, ' ', lam
              Stop ' occlow or occup eq 0 (inner approx.  elements) &
                &in calc_lines'
            End If
            opal(ll) = const*(occlowi-occupi)/clfac(ll)
            opal_plus(ll) = const*occlowi/clfac(ll)
            sline(ll) = c2/(exp(c1/lam/te(ll))-1.D0)/lam3
          End Do

!         print*,'bg approx. with upper level'
!         print*,lam,k,j,low,up
!         do ll=1,nd
!         print*,opal(ll),sline(ll),met_imin(k,ll),met_imax(k,ll)
!         enddo

        End If !                    jatom_full
        icount_bg = icount_bg + 1


!       only lower level  present

      Else

!       prepare two-level approach
!       sline approximated by two level atom; see subr. netmat_met

!       calculate frequency dependent quantities (interpolate)
        ener = log10(lam)
        ifstart = ifre
        Do kk = ifstart, 2, -1
          If (fre(kk)<=ener .And. fre(kk-1)>ener) Go To 100
        End Do
        Stop ' ener not found in fre (calc_lines)'
100     ifstart = kk
        q = (ener-fre(kk))/(fre(kk-1)-fre(kk))
        q1 = 1.D0 - q
        xj1 = q1*xj(:, kk) + q*xj(:, kk-1)
        xk1 = q1*xxk(:, kk) + q*xxk(:, kk-1)
        xj1 = 10.D0**xj1
        xk1 = 10.D0**xk1
        eddtest = xk1(nd)/xj1(nd)
        If (eddtest<0.31 .Or. eddtest>0.35) Stop &
          ' inconsistent Edd.factor in calc_lines'

        lth = lthin(kk)
        Where (lll<lth)
          mustar = 1.D0 - (r(lth)/r)**2 ! correction for lthin (checked)
          mustar = sqrt(mustar)
          eddf = (1.D0-mustar**3)/3.D0/(1.D0-mustar) ! (checked)
        Elsewhere
          eddf = xk1/xj1
        End Where
!       comment: eddf checked in opacity_cmf,
!       smooth transition between outer and inner region,
!       goes to unity for large radii

!       prepare calculation of eps' = Cul/Aul (1-exp(-u0))
!       (see also subroutine opacitl and sumopal)
!       JO (31 Oct. 2016): here was the bug
!       u0=1.4388d8/lam*te
        u0 = 1.4388D8/(lam*te)
        vanreg = collfac*lam**3*(1.D0-exp(-u0)) ! should be OK, but not expl.
!       checked
!       vanreg1 propto lam^2 needs to be finally divided by gf
        vanreg1 = vanreg*c_forbid/lam ! should be OK
        vanreg = vanreg/(1.+vanreg)
!       print*
!       print*,lam,vanreg(1),vanreg1(1)
!       print*,lam,vanreg(43),vanreg1(43)
!       print*,lam,vanreg(nd),vanreg1(nd)

!       line specific quantities, at all depth-points

        If (gf1(ii)<=gfcut) Then
          eps = vanreg1/gf1(ii)
          eps = eps/(1.+eps)
        Else
          eps = vanreg
        End If

        If (jatom_full(k)==1) Then
!         selected elements, division lwion1
          Do ll = 1, lwion1 - 1
!           outside
!           JP changed Nov. 2022, to include lines (approximate treatment)
!           outside
!           imin, imax
            occlowi = occng(irec, low, ll) ! occng = n/g
!           check "exact" NLTE
            If (j>=met_imin(k,ll) .And. j<=met_imax(k,ll)) Then
              If (occlowi<=0.) Then
                Print *, index, ' ', lam
                Stop &
                  ' selected, up=0, outside: occnum = 0 for imin le j le imax'
              End If
            Else
!             approximate solution (outside imin,imax) for selected elements:
!             should happen only for ll < lwion-1
              If (occlowi==0.) Then
                If (ll>=lwion) Then
                  Print *, index, ' ', lam
                  Stop ' approximate elements, up=0, ll ge lwion: occnum=0'
                End If
                opal(ll) = 0.
                opal_plus(ll) = 0.
                sline(ll) = 0.
                Cycle
              End If
            End If
!           either inside imin, imax, or beyond LWION1, or gs,'m' or 's'
!           levels
            opal(ll) = const*occlowi/clfac(ll) ! includes gstat, ind. emission
!           negl.
            opal_plus(ll) = opal(ll)
            taubi = opal(ll)/(eddf(ll)*dvdr(ll)+(1.-eddf(ll))*vr(ll)) ! at
!           mue^2
!           = eddf
!           end change
!           pi e2/me c * 1.d-8 *SR/vmax * lam * gf * nl/gl / clfac /
!           (eddf*dvdr + (1-eddf*v/r))

            If (taubi<1.D-5) Then
              betacic = xj1(ll)
            Else
              betacic = (1.-exp(-taubi))/taubi*xj1(ll)
            End If

            taui = opal(ll)*tau0(ll) ! at mue^2 = 1/3
!           pi e2/me c * 1.d-8 *SR/vmax * lam * gf * nl/gl / clfac / (dvdr/3.
!           + 2./3.*v/r)
            If (taui<1.D-5) Then
              beta = 1.D0
            Else
              beta = (1.-exp(-taui))/taui
            End If

!           now, we can calculate the two-level source-function. We inline
!           Bnue, &
!           and evaluate at actual frequency, to ensure correct thermalization

            epsi = eps(ll)
            eps1i = 1.D0 - epsi
            bni = c2/(exp(c1/lam/te(ll))-1.D0)/lam3
            sline(ll) = (eps1i*betacic+epsi*bni)/(epsi+beta*eps1i)

          End Do
          Do ll = lwion1, nd
!           inside
            occlowi = occng(irec, low, ll) ! occng = n/g
            If (occlowi==0) Then
              Print *, index, ' ', lam
              Stop ' occlow eq 0 (inner selected elements) in opacity_cmf'
            End If
            opal(ll) = const*occlowi/clfac(ll)
            opal_plus(ll) = opal(ll)
            sline(ll) = c2/(exp(c1/lam/te(ll))-1.D0)/lam3
          End Do

!         print*,'selected with lower level only'
!         print*,lam,k,j,low,up
!         do ll=1,nd
!         print*,opal(ll),sline(ll),met_imin(k,ll),met_imax(k,ll)
!         enddo

        Else
!         other elements (approximate), division lwion
          Do ll = 1, lwion - 1
            occlowi = occng(irec, low, ll) ! occng = n/g
            If (occlowi==0) Then
              Print *, index, ' ', lam
              Stop ' occlow eq 0 (outer approx. elements) in opacity_cmf'
            End If
            opal(ll) = const*occlowi/clfac(ll) ! includes gstat, ind. emission
!           negl.
            opal_plus(ll) = opal(ll)
            taubi = opal(ll)/(eddf(ll)*dvdr(ll)+(1.-eddf(ll))*vr(ll)) ! at
!           mue^2
!           = eddf
!           pi e2/me c * 1.d-8 *SR/vmax * lam * gf * nl/gl / clfac /
!           (eddf*dvdr + (1-eddf*v/r))

            If (taubi<1.D-5) Then
              betacic = xj1(ll)
            Else
              betacic = (1.-exp(-taubi))/taubi*xj1(ll)
            End If

            taui = opal(ll)*tau0(ll) ! at mue^2 = 1/3
!           pi e2/me c * 1.d-8 *SR/vmax * lam * gf * nl/gl / clfac / (dvdr/3.
!           + 2./3.*v/r)
            If (taui<1.D-5) Then
              beta = 1.D0
            Else
              beta = (1.-exp(-taui))/taui
            End If

!           now, we can calculate the two-level source-function. We inline
!           Bnue, &
!           and evaluate at actual frequency, to ensure correct thermalization

            epsi = eps(ll)
            eps1i = 1.D0 - epsi
            bni = c2/(exp(c1/lam/te(ll))-1.D0)/lam3
            sline(ll) = (eps1i*betacic+epsi*bni)/(epsi+beta*eps1i)

          End Do
          Do ll = lwion, nd
!           inside
            occlowi = occng(irec, low, ll) ! occng = n/g
            If (occlowi==0) Then
              Print *, index, ' ', lam
              Stop ' occlow eq 0 (inner approx. elements) in opacity_cmf'
            End If
            opal(ll) = const*occlowi/clfac(ll)
            opal_plus(ll) = opal(ll)
            sline(ll) = c2/(exp(c1/lam/te(ll))-1.D0)/lam3
          End Do

!         print*,'bg approx. with lower level only'
!         print*,lam,k,j,low,up
!         do ll=1,nd
!         print*,opal(ll),sline(ll),met_imin(k,ll),met_imax(k,ll)
!         enddo

        End If

        icount_bg_noup = icount_bg_noup + 1

      End If !                      up ne 0
    End If !                        foreground or background

    Do ll = 1, nd
      If (opal(ll)<0.) Then
        If (sline(ll)>0.) Stop ' negative OPAL and positive SLINE'
        inverted_range(ii) = 1
        n_inverted = n_inverted + 1
        Exit
      End If
    End Do

    Do ll = 1, nd
      If (opal(ll)==0.) Then
        If (sline(ll)/=0.) Stop ' OPAL = 0 and SLINE ne 0'
        inverted_range(ii) = 2
        Exit
      End If
    End Do

    maxtau = 0.
    Do ll = 1, lwion - 1
      maxtau = max(maxtau, abs(opal(ll))*tau0(ll))
    End Do
!   JO this is the hack to allow for certain elements only
!   if(maxtau.lt.tauweak  .or. &
!   (k_range(ii).ne.1.and.k_range(ii).ne.2.and.k_range(ii).ne.7)) then
    If (maxtau<tauweak) Then
      inverted_range(ii) = 3
      n_weak = n_weak + 1
    End If

!   if(lam.gt.4471. .and. lam.lt.4473. .and. inverted_range(ii).ne.3) &
!   & print*,index,lam,gf1(ii)

!   WARNING WARNING WARNING
!   lines with inverted_range=2 (opal/sline=0 at some points)
!   can be inverted at different points!!!! Take care!

    If (inverted_range(ii)==0) Then
!     prepare for interpolation
      opal_range(:, ii) = opal(:)*const_opal
      opal_plus_range(:, ii) = opal_plus(:)*const_opal
      sline_range(:, ii) = sline(:)
    Else If (inverted_range(ii)==1) Then
!     special treatment
!     JO previously, we neglected a correct treatment, and used
!     abs(opal)  and abs(sline),
!     since (i) eta = opal*sline is positive and correct for inversion, and
!     (ii) the difference in tau should be negligible in the affected regions.
!     nevertheless, this has to be tested carefully later on
!     opal_range(:,ii)=abs(opal(:))*const_opal
!     sline_range(:,ii)=abs(sline(:))
!     JO May 2017: new treatment, more exact: keep sign
      opal_range(:, ii) = opal(:)*const_opal
      opal_plus_range(:, ii) = opal_plus(:)*const_opal
      sline_range(:, ii) = sline(:)
    Else If (inverted_range(ii)==2) Then
!     keep opal/sline 0 when 0, will be caught in subr. prepray_range
      Do ll = 1, nd
        If (opal(ll)==0) Then
          opal_range(ll, ii) = 0.
          opal_plus_range(ll, ii) = 0.
          sline_range(ll, ii) = 0.
        Else
!         opal_range(ll,ii)=abs(opal(ll))*const_opal
!         sline_range(ll,ii)=abs(sline(ll))
!         JO May 2017: new treatment, more exact: keep sign
          opal_range(ll, ii) = opal(ll)*const_opal
          opal_plus_range(ll, ii) = opal_plus(ll)*const_opal
          sline_range(ll, ii) = sline(ll)
        End If
      End Do
    End If

    If (k_range(ii)<3) icount_hhe = icount_hhe + 1

!   for tests
!   if(inverted_range(ii).ne.3) print*,index,lam,gf1(ii),inverted_range(ii)


  End Do all_lines

! create specific index_file for H/He lines
  Allocate (index_hhe(icount_hhe), l23(icount_hhe), beta_p(icount_hhe))

  icount_hhe = 0
  Do ii = iblu, ired
    If (k_range(ii)<3) Then
      icount_hhe = icount_hhe + 1
      index_hhe(icount_hhe) = ii
    End If
  End Do

  Print *
  Print *, icount_expl + icount_bg + icount_bg_noup, &
    ' line opacities/emissivites calculated'
  Print *, icount_expl, ' lines from explicit elements'
  Print *, icount_hhe, ' lines from H/He'
  Print *, icount_bg, ' lines from background elements with upper level'
  Print *, icount_bg_noup, &
    ' lines from background elements with no upper level'
  Print *, n_weak, ' weak lines with tau_sob < ', tauweak
  Print *, n_inverted, ' inverted lines'

! calculation of Stark broadening profiles etc.
  If (icount_hhe>maxncomp) Then
    Print *, 'NO OF H/HE LINES:', icount_hhe
    Stop ' TOO MANY H/HE LINES FOR STARK BROADENING (INCREASE MAXNCOMP)'
  End If

  Allocate (level(icount_hhe), leveu(icount_hhe), lineno(icount_hhe), &
    nlevl(icount_hhe), nlevu(icount_hhe), nstark(icount_hhe))

  nstark = 1
  lineno = 0

  Do i = 1, icount_hhe
    ii = index_hhe(i)
    index = id1(ii)
    If (index<1D8) Stop ' problems with index (HHe)'
    kk = index/1D8
    irest = index - kk*1D8
    ml = irest/10000
    mu = irest - ml*10000
    level(i) = labl(ml)
    leveu(i) = labl(mu)
  End Do

! calculated with vturbmin=vturb
  Call preformal(modnam, icount_hhe, level, leveu, lineno, nstark, vturb, yhe, &
    balmer_lemke, pb_lemke, .True.)

! test IX.DAT file (should be sorted for frequency) and consistency
  Open (1, File=trim(modnam)//'/IX.DAT', Status='unknown', Form='formatted')
  Rewind 1

  Read (1, *) ncomp
  If (ncomp/=icount_hhe) Stop ' problem with no of components in IX.DAT'

  Do i = 1, ncomp
    ii = index_hhe(i)
    Read (1, *) xlamb
    If (abs(xlamb-xlam1(ii))>0.1) Print *, ' check wavelengths (IX.DAT)! ', &
      xlamb, ' ', xlam1(ii)
    Read (1, *) gll, glu
    Read (1, *) gf
    If (abs(log10(gf/gf1(ii)))>0.1) Print *, &
      ' check oscillator strength (IX.DAT! ', gf, ' ', gf1(ii)
    Read (1, *) nlevl(i), nlevu(i)
    Read (1, *) weight
    Read (1, *) nst
  End Do
  Close (1)

  Print *
  Print *, ' Stark profiles prepared, all data consistent'

  Deallocate (level, leveu, lineno)

! for all Stark-broadened lines, calculate OPAL/OPAC  at tauc=0.66
! further used in calcxmax_range

  ifstart = ifre

  Do i = 1, icount_hhe
    ii = index_hhe(i)
    lam = xlam1(ii)

!   determine (approx) where tauc=2/3
    ener = log10(lam)
    Do kk = ifstart, 2, -1
      If (fre(kk)<=ener .And. fre(kk-1)>ener) Go To 110
    End Do
    Stop ' ener not found (2) in fre (calc_lines)'

110 ifstart = kk

!   use kk to allow for larger opacity in case of edges, to be on the safe
!   side
    tauc = 0.
    Do ll = 2, nd
      tauc = tauc + 0.5*(opac(ll-1,kk)+opac(ll,kk))*(r(ll-1)-r(ll))*sr
      If (tauc>0.66) Go To 120
    End Do
    Stop ' tauc = 2/3 not found (calc_lines)'

120 l23(i) = ll
    beta_p(i) = opal_range(ll, ii)*const_opal_stark/lam/(opac(ll,kk)*sr)
!   print*,i,lam,l23(i),beta_p(i)
  End Do

  Return
End Subroutine

!-----------------------------------------------------------------------

Subroutine formal_range(nd, np, nc, vdop, rmax, r, p, z, v, ltot, lmax, r1, &
  xne1, temp1, index1, iblu, ired, lam_blu, lam_red, resol, sr, vmax, &
  opacon, scont, zray, opacray, sconray, opalray, sray)

! here: nd actual length of microgrid
! ndm = id_depfi maximum length
! nd1 = id_ndept length of standard grid
! variables with ending'1' refer to standard grid

  Use :: nlte_type
  Use :: nlte_dim, Only: id_ndept, id_depfi, id_cores, id_nocor, id_maxww, &
    id_maxtt, id_maxne, id_maxps
  Use :: fund_const, Only: pi
  Use :: formalsol_var, Only: nws, nts, nes, qhalf, dws, ts, es, ps, &
    modnam => file, nstark, nlevl, nlevu
  Use :: formalsol_add_var, Only: lamgrid, icount_hhe, dlam_max_stark, qstore, &
    lstore
  Implicit None


! .. parameters ..
  Integer (i4b), Parameter :: nd1 = id_ndept, ndm = id_depfi, lto1 = 2*ndm
  Integer (i4b), Parameter :: nc1 = id_cores
  Integer (i4b), Parameter :: np1 = nc1 + 2 + 4*id_nocor
! ..
! .. scalar arguments ..
  Real (dp), Intent (In) :: vdop, rmax, lam_blu, lam_red, resol, sr, vmax
  Integer (i4b), Intent (In) :: iblu, ired

  Integer (i4b) :: nd, np, nc

! ..
! .. array arguments ..

  Real (dp) :: r(ndm), v(ndm), p(np1), z(ndm, np1), r1(nd1), &
    xne1(nd1), temp1(nd1)

  Integer (i4b) :: index1(nd1)


  Real (dp), Dimension (lto1) :: opacon, scont ! note, dimension changed from
! (ndm,2) to 2*ndm
  Real (dp), Dimension (lto1) :: zray, opacray, sconray, opalray, sray

  Integer (i4b) :: lmax(np1), ltot(np1)
! ..
! .. local scalars ..

  Real (dp) :: delp, dd, z2, dtdr, xlambda, aic, dbdr, w, w1, ww, ww1, xicor, &
    emint, emint1, errmax, relem, t_start, t_end, time

  Integer (i4b) :: jp, ndelt, irest, istart, ijp, l, j, nf, nfcmf, lm, k, lm1, &
    ltaumax, nsum

  Logical :: core, first
! ..
! .. local arrays ..
  Real (dp) :: temp(ndm), xnefi(ndm), weight_range(ndm), tempray(lto1), &
    xneray(lto1), taucon(lto1), fscon(lto1), tauray(lto1)
! ..
  Real (dp), Dimension (:), Allocatable :: tau, s1pabs, s1pem, s2pem
  Real (dp), Dimension (:), Allocatable :: conabs, conem1, conem2, profcon, &
    profabs, profem, profile, profred
  Integer (i4b), Dimension (:), Allocatable :: nl, nu
  Real (dp), Dimension (:, :), Allocatable :: zg

! .. external functions ..
  Real (dp) :: bnue, dbdt
  External :: bnue, dbdt

  If (nd>ndm) Stop ' nd > ndm in formal_range'

! note: until here, vdop=vturb/vinf

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
    Stop ' something wrong with number of p-rays'
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

! calculation of the z-grid (fine)

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

! vplus = 0.d0
! do i = 1,nb
! vplus = max(vplus,xmax(i))
! end do

! print *
! print *,' stark/doppler broadening effective for xmax = ',vplus
! print *

  Call xgrid_range(lam_blu, lam_red, resol, nf)

  Call interpol_range(r, temp, temp1, xnefi, xne1, index1, weight_range)

  Allocate (tau(nf), s1pabs(nf), s1pem(nf), s2pem(nf), profabs(nf), &
    profem(nf), profile(nf), profred(nf))
  Allocate (conabs(nf), conem1(nf), conem2(nf), profcon(nf))
  Allocate (zg(np1,nf))

  errmax = 0.D0

  tau = 0.D0
! so far, tau not used; if needed, check formalsol how to implement
  s1pabs = 0.D0
  s1pem = 0.D0
  s2pem = 0.D0
  conabs = 0.D0
  conem1 = 0.D0
  conem2 = 0.D0

  dtdr = (temp1(nd1)-temp1(nd1-1))/(r1(nd1-1)-r1(nd1))

! read Stark-profiles
  nsum = sum(nstark)

  If (nsum>0) Then
!   nstar, nlevl, nlevu already allocated in calc_lines
    Allocate (nws(icount_hhe), nts(icount_hhe), nes(icount_hhe), &
      qhalf(icount_hhe))
    Allocate (nl(icount_hhe), nu(icount_hhe), dlam_max_stark(icount_hhe))
    Allocate (dws(id_maxww,icount_hhe))
    Allocate (ts(id_maxtt,icount_hhe))
    Allocate (es(id_maxne,icount_hhe))
    Allocate (ps(id_maxps,icount_hhe))
    nl = nlevl
    nu = nlevu
    Call staread(modnam, icount_hhe, nsum, nl, nu)
    Call calcxmax_range(temp1, xne1)
  End If

  time = 0.

! calculate opacities/emissivities on coarse grid for cmf-frequencies
! will be interpolated onto obs. frame freq. in routine prepray_range
  Call opal_cmfgrid(nfcmf, resol, vdop, vmax, temp1, xne1, iblu, ired)


jploop: Do jp = 1, np - 1
!   jploop: do jp = 1,1
    lm = lmax(jp)
    core = lm == nd
    first = .True.

kloop: Do k = 1, nf
!     kloop: do k = 4,4
      xlambda = lamgrid(k)
      aic = bnue(xlambda, temp1(nd1))
      dbdr = dbdt(xlambda, temp1(nd1))*dtdr

      zg(jp, k) = -10.*rmax

!     print*,k,jp
      Call cpu_time(t_start)
      Call prepray_range(first, nfcmf, xlambda, jp, z(1,jp), p, ltot, np, &
        lm, core, w, w1, vdop, r, v, temp, xnefi, r1, index1, sr, &
        vmax, opacon, scont, weight_range, zray, tempray, xneray, opacray, &
        sconray, opalray, sray, tauray)
      Call cpu_time(t_end)
      time = time + t_end - t_start
      first = .False.

      ww = w*pi/(r(1)*r(1))
      ww1 = w1*pi/(r(1)*r(1))

      lm1 = ltot(jp)

      Call formacon_range(taucon, ltaumax, fscon, opacray, sconray, zray, lm1)
!     ltaumax index where tau_cont > taumax

      If (core) xicor = aic + z(lm1, jp)*dbdr/opalray(lm1) ! for lines + cont

      Call obsfram_range(lm1, core, zray, opalray, sray, tauray, xicor, &
        zg(jp,k), emint, emint1)

!     zg not tested so far

      s1pabs(k) = s1pabs(k) + ww*emint1
      s1pem(k) = s1pem(k) + ww*emint
      If (core) Then
        xicor = aic + z(lm1, jp)*dbdr/opacray(lm1) ! for cont
        conabs(k) = conabs(k) + ww*xicor*exp(-taucon(lm1))
      End If
      conem1(k) = conem1(k) + ww*fscon(1)

      If (jp>=nc+2) Then
        s2pem(k) = s2pem(k) + ww1*emint
        conem2(k) = conem2(k) + ww1*fscon(1)
      End If

    End Do kloop
    Print *, jp, ' from a total of ', np - 1, ' rays done'

  End Do jploop

  Do k = 1, nf
    profcon(k) = conem1(k) + conem2(k) + conabs(k)
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
    profred(k) = profile(k)/profcon(k)
  End Do

  Do k = 1, nf
    xlambda = lamgrid(k)
    Write (*, Fmt=100) k, xlambda, profcon(k), profred(k), profile(k)
    Write (2, Fmt=100) k, xlambda, profcon(k), profred(k), profile(k)
  End Do

  errmax = 15.D0*errmax
  If (errmax>.03) Print *, ' WARNING!!! WARNING!!! WARNING!!! WARNING!!!'

  Print *
  Print *, ' MAXIMUM ERROR IN ANGULAR INTEGRATION ', errmax
  Print *

  If (errmax>.03) Then
    Print *, ' WARNING!!! WARNING!!! WARNING!!! WARNING!!!'
    Print *
  End If

  Deallocate (tau, s1pabs, s1pem, s2pem, profabs, profem, profile, profred)
  Deallocate (conabs, conem1, conem2, profcon)
  Deallocate (zg)
  Deallocate (nlevl, nlevu, nstark)
  Deallocate (qstore, lstore)

  If (nsum>0) Then
    Deallocate (dlam_max_stark, nl, nu)
    Deallocate (nws, nts, nes, qhalf)
    Deallocate (dws, ts, es, ps)
  End If

  Print *, 'time=', time

  Return

100 Format (1X, I5, 4(2X,G14.6))

End Subroutine

!-----------------------------------------------------------------------


Subroutine xgrid_range(lam_blu, lam_red, resol, nf)

  Use :: nlte_type
  Use :: formalsol_add_var, Only: lamgrid

  Implicit None

! .. parameters ..
  Real (dp), Intent (In) :: lam_blu, lam_red, resol
  Integer (i4b), Intent (Out) :: nf

  Integer (i4b) :: i

  nf = int((lam_red-lam_blu)/resol) + 1

  Allocate (lamgrid(nf))

  Do i = 1, nf
    lamgrid(i) = lam_blu + (i-1)*resol
  End Do
  If (abs(lamgrid(nf)-lam_red)>resol) Stop ' error in xgrid_range'
  Print *, ' ACTUAL lambda_red = ', lamgrid(nf)

  Return
End Subroutine

!-----------------------------------------------------------------------

Subroutine interpol_range(r, temp, t, xnefi, xne, index1, weight_range)

  Use :: nlte_type
  Use :: nlte_dim, Only: id_ndept, id_depfi

  Implicit None

  Integer (i4b), Parameter :: nd1 = id_ndept, ndm = id_depfi

  Real (dp), Dimension (ndm), Intent (In) :: r
  Real (dp), Dimension (nd1), Intent (In) :: t, xne
  Real (dp), Dimension (ndm), Intent (Out) :: temp, xnefi, weight_range

  Integer (i4b), Dimension (nd1), Intent (In) :: index1

  Integer (i4b) :: i, j, j1, j2

  Real (dp) :: rlj, rlj1, rlj2, tl, tl1, tl2, factor

! calculating interpolation weights, and
! interpolation of temperature and n_e (and other quantities when needed)
  weight_range = 0.

  temp(1) = t(1)
  xnefi(1) = xne(1)
iloop: Do i = 2, nd1
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
      tl1 = log(xnefi(j1))
      tl2 = log(xnefi(j2))
      tl = tl1 + factor*(tl2-tl1)
      xnefi(j) = exp(tl)
      weight_range(j) = factor
    End Do jloop
  End Do iloop

  Return
End Subroutine

!-----------------------------------------------------------------------

Subroutine opal_cmfgrid(nfcmf, resol, vvdop, vmax, temp, xne, iblu, ired)

! calculates total line opacities and source functions on an appropriate
! cmf frequency grid, for the coarse spatial grid.
! Later on, these will be interpolated onto corresponding obs. frame
! frequencies
! and the fine grid

  Use :: nlte_type
  Use :: fund_const, Only: clight, amh, akb, sqpi
  Use :: nlte_dim, Only: nd1 => id_ndept
  Use :: formalsol_var, Only: nstark, vturbv
  Use :: formalsol_add_var, Only: natom, aweight, lamb1, lamr1, xmax, &
    first_hel, xlam1, k_range, inverted_range, vmax_broad, icount_hhe, &
    index_hhe, dlam_max_stark, opal_range, opal_plus_range, sline_range, &
    opal_lgrid, opal_part_lgrid, sl_lgrid, lamgrid_cmf, opalmin_lgrid

  Implicit None
  Real (dp), Parameter :: const_vth = 2.D0*akb/amh
  Real (dp), Parameter :: xnemin = 10.D0**10.5
  Real (dp), Parameter :: const_opal_stark = sqpi*clight*1.D8


  Real (dp), Intent (In) :: resol, vvdop, vmax
  Real (dp), Dimension (nd1), Intent (In) :: temp, xne

  Integer (i4b), Intent (In) :: iblu, ired
  Integer (i4b), Intent (Out) :: nfcmf

  Real (dp), Allocatable, Dimension (:, :) :: vth

  Real (dp) :: vturb, vturb2, vth2, xlamcmf, lammax, lammax_old, lammin, &
    lammax_hhe, lammax_hhe_old, xl, xk, opal, opal_plus, sl, perf, perfil2, &
    resol_cmf, etal, opalmin, prof
  Integer (i4b) :: i, j, k, l, kel, indstart, indend, indstart_hhe, &
    indend_hhe, loop_count, nlines, nlines_hhe

  Logical :: flag

! thermal speeds for all elements on coarse grid (in units of clight)
! vvdop is vturb/vinf
  Allocate (vth(nd1,natom))
  vturb = vvdop*vmax

  If (abs(1.-vturb/vturbv(nd1))>1.D-12) Stop &
    ' something rotten with vturb/vturbv in opal_cmfgrid'

  Do l = 1, nd1
    vth2 = const_vth*temp(l)
    vturb2 = vturbv(l)**2
    Do kel = 1, natom
      vth(l, kel) = sqrt(vturb2+vth2/aweight(kel))/clight
    End Do
  End Do

! create cmf-frequency grid with minimum vturb
  resol_cmf = vturb/3.D0/clight*lamb1
! JO: test whether resolution with vturb (instead of vturb/3) is sufficient
  resol_cmf = min(resol_cmf, resol)
  nfcmf = (lamr1-lamb1)/resol_cmf + 2

  Allocate (lamgrid_cmf(nfcmf))
  Do k = 1, nfcmf
    lamgrid_cmf(k) = lamb1 + (k-1)*resol_cmf
  End Do
  If (lamgrid_cmf(nfcmf)<=lamr1) Stop ' error in lamgrid_cmf'
  Print *
  Print *, ' cmf frequency grid with resolution ', resol_cmf, ' (A) '
  Print *, ' and ', nfcmf, ' frequency points allocated'

  Allocate (opal_lgrid(nd1,nfcmf), opal_part_lgrid(nd1,nfcmf), &
    sl_lgrid(nd1,nfcmf), opalmin_lgrid(nfcmf))
  opal_lgrid = 0.
  opal_part_lgrid = 0.
  sl_lgrid = 0.
  opalmin_lgrid = 0.

  Do k = 1, nfcmf

    xlamcmf = lamgrid_cmf(k)

    indstart = ired
    indstart_hhe = icount_hhe
    lammax = 1.D20
    lammax_hhe = 1.D20

    Do l = nd1, 1, -1
!     calculate range of cmf-frequencies which affect the given (cmf-) freq.
!     at first for non-HHe lines, because of the lower vtherm
!     order from inside to ouside, since (mostly) vth and thus lammax
!     decreases
      lammax_old = lammax
      lammax = xlamcmf*(1.D0+xmax*vth(l,first_hel))
      lammin = xlamcmf*(1.D0-xmax*vth(l,first_hel))

      loop_count = 0
      If (lammax>lammax_old) indstart = ired ! if vth increasing
!     for tests
!     if(indstart.le.ired-1 .and. xlam1(indstart+1).le.lammax) then
!     print*,l,indstart,lammax,xlam1(indstart+1)
!     stop ' error in indstart (prepray_range)'
!     endif
      Do i = indstart, iblu, -1
        xl = xlam1(i)

        If (xl<=lammax .And. xl>=lammin) Then
          If (loop_count==0) indstart = i
          indend = i
          loop_count = loop_count + 1
        End If
        If (xl<lammin) Exit
      End Do
      If (loop_count/=0) Then
        nlines = indstart - indend + 1
      Else
        indend = indstart + 1 !     hack to ensure that loop below is not
!       executed
        nlines = 0
      End If

!     remember once more: lines with inverted_range=2 (opal/sline=0 at some
!     points)
!     can be inverted at different points!!!! Take care!

      Do i = indend, indstart !     works also for nlines = 0, since then
!       indend > indstart
        kel = k_range(i)
        If (inverted_range(i)==3) Cycle ! weak lines with tau_sob < tauweak
        If (kel<3) Cycle !          special treatment for H/He
        xk = (xlam1(i)/xlamcmf-1.D0)/vth(l, kel) ! v/c / vth/c = v/vth

        If (abs(xk)>4.) Cycle !     only for Doppler
        If (inverted_range(i)==2) Then
!         no interpolation, opal set to zero
          If (opal_range(l,i)==0. .Or. opal_range(l+1,i)==0.) Cycle
        End If
!       opal_range already multiplied by sr*lam*1.d-8/sqpi/c (subr.
!       calc_lines)
!       JO so far, simple Doppler profile. Might be improved by Voigt profile
        prof = exp(-xk*xk)/vth(l, kel)
        opal = opal_range(l, i)*prof
        opal_plus = opal_plus_range(l, i)*prof
        sl = sline_range(l, i)
        opal_lgrid(l, k) = opal_lgrid(l, k) + opal
        opal_part_lgrid(l, k) = opal_part_lgrid(l, k) + opal_plus
!       JO May 2017 for tests, can be skipped later
        etal = opal*sl
        If (etal<=0.) Then
          Print *, kel, xlam1(i), l, opal_range(l, i), sl, inverted_range(i), &
            i
          Stop ' error in etal (metals) - subr. opal_cmfgrid'
        End If
        sl_lgrid(l, k) = sl_lgrid(l, k) + etal
      End Do

      If (icount_hhe/=0) Then
        lammax_hhe_old = lammax_hhe
!       special treatment of H/He lines (Stark-broadened)
!       calculate range of cmf-frequencies which affect the given obs. frame
!       freq.
        If (xne(l)<xnemin) Then
!         note also that xnemin can be met twice (when clumping included);
!         use Doppler-with for hydrogen
          lammax_hhe = xlamcmf*(1.D0+xmax*vth(l,1))
          lammin = xlamcmf*(1.D0-xmax*vth(l,1))
          If (lammax_hhe>lammax_hhe_old) indstart_hhe = icount_hhe ! if vth
!         increasing
        Else
!         reset of indstart (because of broader profile)
          indstart_hhe = icount_hhe
!         use maximum Stark-width
          lammax_hhe = xlamcmf*(1.D0+vmax_broad/clight)
          lammin = xlamcmf*(1.D0-vmax_broad/clight)
        End If

        loop_count = 0
!       for tests
!       if(indstart_hhe.le.icount_hhe-1 .and. &
!       &        xlam1(index_hhe(indstart_hhe+1)).le.lammax_hhe) then
!       print*,l,indstart_hhe,lammax_hhe,xlam1(index_hhe(indstart_hhe+1))
!       stop ' error in indstart_hhe (prepray_range)'
!       endif
        Do i = indstart_hhe, 1, -1
          j = index_hhe(i)
          xl = xlam1(j)
          If (xl<=lammax_hhe .And. xl>=lammin) Then
            If (loop_count==0) indstart_hhe = i
            indend_hhe = i
            loop_count = loop_count + 1
          End If
          If (xl<lammin) Exit
        End Do
        If (loop_count/=0) Then
          nlines_hhe = indstart_hhe - indend_hhe + 1
        Else
          nlines_hhe = 0
        End If

        If (nlines_hhe/=0) Then
          Do i = indend_hhe, indstart_hhe
            j = index_hhe(i)
            kel = k_range(j)
            If (inverted_range(j)==3) Cycle ! weak lines with tau_sob <
!           tauweak
            If (inverted_range(j)==2) Stop ' OPAL/SLINE (H/He) = 0!'
            If (kel>=3) Stop ' problem with index_hhe'

            If (xne(l)<xnemin .Or. nstark(i)==0) Then
!             consistency check
!             JO May 2017: check not possible when all nstark = 0, since then
!             dlam_max_stark not allocated
!             if(nstark(i).eq.0 .and. dlam_max_stark(i).ne.0.) &
!             &             stop ' HHe line with nstark = 0 and dlam_max_stark
!             ne 0'
              xk = (xlam1(j)/xlamcmf-1.D0)/vth(l, kel) ! v/c / vth/c = v/vth

!             Doppler profile (when no other info available)
              If (abs(xk)>4.) Cycle
!             opal_range already multiplied by sr*lam*1.d-8/sqpi/c (subr.
!             calc_lines)
              prof = exp(-xk*xk)/vth(l, kel)
              opal = opal_range(l, j)*prof
              opal_plus = opal_plus_range(l, j)*prof
              sl = sline_range(l, j)
            Else If (nstark(i)==1) Then
              xl = xlam1(j)
              If (abs(xl-xlamcmf)>dlam_max_stark(i)) Cycle
!             Stark broadening, normalized corresponding to
!             exp(-xk*xk)*lambda/sqpi/vdop
!             (or min(vdop) if depth dependent vturb)
              perf = perfil2(i, xl, xlamcmf, temp(l), xne(l), flag)
              If (flag) Cycle
!             opal_range already multiplied by sr*lam*1.d-8/sqpi/c (subr.
!             calc_lines)
!             to account for normalization of perfil2 (see above)
!             we have to multiply OPAL by clight*sqpi*1.d8 and to divide by
!             xlam (in A)
              prof = perf*const_opal_stark/xl
              opal = opal_range(l, j)*prof
              opal_plus = opal_plus_range(l, j)*prof
              sl = sline_range(l, j)
            Else
              Stop ' NSTARK > 1 in prepray_range!'
            End If
            opal_lgrid(l, k) = opal_lgrid(l, k) + opal
            opal_part_lgrid(l, k) = opal_part_lgrid(l, k) + opal_plus
!           JO May 2017 for tests, can be skipped later
            etal = opal*sl
            If (etal<=0.) Then
              Stop ' error in etal (HHe) - subr. opal_cmfgrid'
            End If
            sl_lgrid(l, k) = sl_lgrid(l, k) + etal
          End Do
        End If
      End If !                      HHe treatment

!     for tests
!     print*,k,l,sl_lgrid(l,k)/opal_lgrid(l,k)/bnue(xlamcmf,temp(l))
!     if(opal_lgrid(l,k).eq.0.) print*,l,k,' opal_lgrid=0'
    End Do !                        l-loop
  End Do !                          k-loop

! JO May 2017: new treatment of inversion:
! interpolate absorption and ind. emission part seperately
  Do k = 1, nfcmf
    opalmin = minval(opal_lgrid(:,k))
    If (opalmin<0.D0) Then
      opalmin_lgrid(k) = opalmin
!     print*,k,lamgrid_cmf(k),opalmin
!     print*,opal_lgrid(:,k)
!     print*,sl_lgrid(:,k)
!     print*
!     from here on, opal_lgrid contains absorption, and
!     opal_part_lgrid stim. emssion, whenever opalmin_lgrid <0
      Do l = 1, nd1
        etal = opal_lgrid(l, k) !   abs - em
        opal = opal_part_lgrid(l, k) ! abs
        opal_part_lgrid(l, k) = opal - etal ! now emission part, should be > 0
!       now
        opal_lgrid(l, k) = opal !   now absorption part
        If (opal_lgrid(l,k)<0. .Or. opal_part_lgrid(l,k)<0.) &
          Stop &
          ' absorption or emission part of opacity <= 0 (subr. opal_cmfgrid)!'
      End Do
    End If
  End Do

  Where (opal_lgrid==0.D0)
    opal_lgrid = 1.D-40 !           for log-log interpolation
  End Where

  Where (opal_part_lgrid==0.D0)
    opal_part_lgrid = 1.D-40 !      for log-log interpolation
  End Where

  opal_lgrid = log(opal_lgrid)
  opal_part_lgrid = log(opal_part_lgrid)

  Where (sl_lgrid==0.)
    sl_lgrid = 1.D-40 !             for log-log interpolation
  End Where
  sl_lgrid = log(sl_lgrid)

  Deallocate (vth)

  Print *
  Print *, ' total opacities and source-functions calculated on CMF grid'
  Print *, ' (low spatial resolution)'
  Print *

  Return

End Subroutine

!-----------------------------------------------------------------------

Subroutine prepray_range(first, nfcmf, xlamobs, jp, z, p, ltot, np, lmax, &
  core, w, w1, vvdop, r, v, temp, xne, r1, index1, sr, vmax, &
  opacon, scont, weight_range, zray, tempray, xneray, opacray, sconray, &
  opalray, sray, tauray)

! vvdop is vturb/vinf; thus far, it's no longer used

! remember: lines with inverted_range=2 (opal/sline=0 at some points)
! can be inverted at different points!!!!
! might be few weak lines (inverted_range=3) which have an influence on the
! spectrum, when strong in intermediate/lower photosphere

  Use :: nlte_type
  Use :: fund_const, Only: clight
  Use :: nlte_dim
  Use :: formalsol_var, Only: vzrp, xcmfp, wp, wp1, modnam => file
  Use :: formalsol_add_var, Only: taumax, qstore, lstore, opal_lgrid, &
    opal_part_lgrid, sl_lgrid, lamgrid_cmf, opalmin_lgrid

  Implicit None

! defining whole rays including backward half-sphere

! .. parameters ..
  Integer (i4b), Parameter :: ndm = id_depfi, lto1 = 2*ndm, nc = id_cores
  Integer (i4b), Parameter :: nd1 = id_ndept, ifretot = id_frec1
  Integer (i4b), Parameter :: np1 = nc + 2 + 4*id_nocor

! ..
! .. scalar arguments ..
  Real (dp), Intent (In) :: vvdop, xlamobs, sr, vmax
  Real (dp), Intent (Out) :: w, w1

  Integer (i4b), Intent (In) :: nfcmf, jp, np, lmax

  Logical, Intent (In) :: first, core
! ..
! .. array arguments ..
  Real (dp), Dimension (ndm), Intent (In) :: z, r, v, temp, xne, weight_range
  Real (dp), Dimension (np1), Intent (In) :: p

  Real (dp), Dimension (lto1) :: opacon, scont, opacray, sconray, opalray, &
    sray
  Real (dp), Dimension (lto1) :: zray, tempray, xneray, tauray

  Real (dp), Intent (In) :: r1(nd1)
  Integer (i4b), Intent (In) :: index1(nd1)
  Integer (i4b), Intent (Out) :: ltot(np1)
! ..
! .. local scalars ..

  Real (dp) :: vzr, vzr1, delv, zi, z1, deltaz, delz, t1, delt1, xn1, delxn1, &
    delzi, dellin, ddp, dpp3, dpold, dpp3old, sum1, sum2, sum3, sumex, &
    enerobs, enercmf, opacon_last, scont_last, q, tl1, tl2, factor, rmin, p2, &
    logr, q1, opal, sl, xlamcmf, opali, opali1, sli, sli1, ww, ww1, &
    opali_part, opali1_part, opal_em

  Real (dp), Save :: vc, vdopmin, vdopmax


  Integer (i4b) :: i, j, l, ll, lm1, l1, ndel, ip, irest, kkstart, kk, j1, j2, &
    lmax1, ltot1, jj, jj1, jj2, indstart, lmz0, irstart


  Integer (i4b), Save :: ifre

! .. local arrays ..
  Real (dp), Dimension (nd1), Save :: logr1

  Real (dp), Dimension (lto1), Save :: zcon, qarr
  Integer (i4b), Dimension (lto1), Save :: indexz

  Real (dp), Dimension (ifretot), Save :: fre
  Real (dp), Dimension (nd1, ifretot), Save :: opac, scon

  Logical :: start, optz0
! ..

  Data start/.True./

  If (start) Then

    logr1 = log(r1)
    vc = vmax/clight
    vdopmin = vvdop/5.D0
    vdopmax = vdopmin*2.D0

    Open (1, File=trim(modnam)//'/CONT_FORMAL_ALL', Status='unknown', &
      Form='unformatted')

    Read (1) ifre, (fre(i), i=1, ifretot)
    Close (1)

    Open (1, File=trim(modnam)//'/CONT_FORMAL_CMF', Status='unknown', &
      Form='unformatted')

!   cont source function (CMF)
    Read (1)((scon(i,j),i=1,nd1), j=1, ifretot), ((opac(i,j),i=1,nd1), j=1, &
      ifretot) !                    opac_nolines
    Close (1)

!   for interpolation below
    Do i = 1, nd1
      Do j = 1, ifre
        scon(i, j) = log(scon(i,j))
        opac(i, j) = log(opac(i,j))
      End Do
    End Do

    start = .False.
  End If

  If (first) Then
!   calculate frequency independent quantities along ray
!   (so far, z, temp, xne)

    optz0 = .False.

!   -----core rays and first point of non-core rays without interpolation

!   NOTE: in contrast to standard approach, xcmpf now is +mu(r)*v(r)
!   don't fiddle around with vzrp, since this is  an auxillary array needed
!   (as it is) to do a correct interpolation!!!

    If (core) Then
      ltot(jp) = lmax
      lm1 = lmax
      ltot1 = lmax
    Else
      lm1 = 1
    End If

    Do l = 1, lm1
      vzr = v(l)*z(l)/r(l)
      xcmfp(l) = vzr !              different
      tempray(l) = temp(l)
      xneray(l) = xne(l)
      zray(l) = z(l)
    End Do

!   -----non core rays with interpolation, and definition of zcon (zgrid for
!   opacon and scont)
!   difference between zcon and zray:
!   zcon corresponds to (z,0,-z): thus dim(zcon) = ltot1=2*lm1-1
!   with lm1=lmax (optz0) or lm1=lmax+1 (.not.optz0) before interpol.
!   note that after interpolation lm1 changes value,
!   becomes total number of points including interpolated ones
!   zray includes interpolated values if delv too large; dim(zray)=ltot(jp)

    If (.Not. core) Then
!     consistent with calculation of opacon and scont below
      zcon(1:lmax) = z(1:lmax)

      If (z(lmax)==0.D0) Then
        optz0 = .True.
        lm1 = lmax
      Else
        lm1 = lmax + 1
        zcon(lm1) = 0.D0
      End If

      ltot1 = 2*lm1 - 1

      Do l = lm1 + 1, ltot1
        jj = 2*lm1 - l
        zcon(l) = -z(jj)
      End Do

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
          If (i>lto1) Stop ' more than lto1 points used in interpolation!'
          xcmfp(i) = vzr1 !         different

          If (optz0 .Or. l1/=lm1) Then

            tempray(i) = temp(l1)
            xneray(i) = xne(l1)
            zray(i) = z(l1)

!           different treatment for last point with optz0 = false, assuming
!           constant source function and opacity from z(lmax) to z=0.

          Else
            tempray(i) = temp(l)
            xneray(i) = xne(l)
            zray(i) = 0.D0
          End If

!         interpolation linearly with z

        Else

          ndel = int(delv/vdopmin) + 1
          If (ndel<2) Stop ' ndel < 2'
          zi = z(l)

          If (optz0 .Or. l1/=lm1) Then
            z1 = z(l1)
            deltaz = zi - z1
            delz = deltaz/dble(ndel)
            deltaz = 1.D0/deltaz
            t1 = temp(l)
            delt1 = temp(l1) - temp(l)
            xn1 = xne(l)
            delxn1 = xne(l1) - xne(l)

          Else

!           different treatment for last point with optz0 = false, assuming
!           constant source function and opacity from z(lmax) to z=0.

            z1 = 0.D0
            deltaz = zi - z1
            delz = deltaz/dble(ndel)
            deltaz = 1.D0/deltaz
            t1 = temp(l)
            delt1 = 0.D0
            xn1 = xne(l)
            delxn1 = 0.D0

          End If

kloop:    Do kk = 1, ndel
            i = i + 1
            If (i>lto1) Stop ' more than lto1 points used in interpolation!'
            delzi = kk*delz
            dellin = delzi*deltaz
            zray(i) = zi - delzi

!           attention!! explicitly accounted for sign of delv

            xcmfp(i) = vzr - delv*dellin ! different
            If (kk==ndel) Then
!             can happen due to calculation, and leads to inconsistent grids
              If (abs(zray(i)-z1)>1.D-10) Stop &
                ' error in zray(ndel) -- prepray-range'
              zray(i) = z1
              xcmfp(i) = vzr1
            End If
            tempray(i) = t1 + delt1*dellin
            xneray(i) = xn1 + delxn1*dellin
          End Do kloop

        End If

      End Do lloop1

!     defining quantities on the back hemisphere

      lm1 = 2*i - 1
      If (lm1>lto1) Stop ' more than lto1 points used in interpolation!'
      ltot(jp) = lm1

      Do l = 1, i - 1
        ll = 2*i - l
        vzr = -xcmfp(l) !           change of sign remains
        xcmfp(ll) = vzr
        tempray(ll) = tempray(l)
        xneray(ll) = xneray(l)
        zray(ll) = -zray(l)
      End Do

!     end of non-core ray treatment

    End If

!   for tests
!   print*,'jp=',jp
!   do ll=1,lm1
!   write(*,fmt='(i4,2(2x,f10.4),2x,e10.4)')
!   ll,zray(ll),log10(xneray(ll)),xcmfp(ll)
!   enddo

!   calculation of vtherm (in units of c)
    If (lm1/=ltot(jp)) Stop ' problem with lm1'

!   now, xcmfp is mu*v(r)/c with v(r) in actual units
    xcmfp(1:lm1) = xcmfp(1:lm1)*vc

!   set up index array for interpolation of opacon to opacray etc.
    If (lm1/=ltot1) Then
      kkstart = 1
      indexz = 0
      qarr = 1.D0
      Do i = 1, ltot1
        Do l = kkstart, lm1
          If (zray(l)==zcon(i)) Then
!           print*,jp,l,i,zray(l)
            indexz(l) = i
            kkstart = l + 1
            Go To 100
          End If
          indexz(l) = -i
          qarr(l) = (zray(l)-zcon(i))/(zcon(i-1)-zcon(i))
        End Do
        Print *, 'not found', jp, l, i, zcon(i)
        Stop 'zray not found in zcon (prepray_range)'
100     Continue
      End Do
!     for tests
!     do l=1,lm1
!     i=indexz(l)
!     if(i.ne.0) print*,jp,l,zray(l),i,zcon(abs(i))
!     enddo
    End If

!   calculate w = weights for flux integral

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

      If (abs(sum2)>1.D-5) Stop 'error in sum2'
      sumex = r(1)**2

      If (abs(sum1+sum3-sumex)>1.D-2) Then
        Print *, sum1 + sum3 - sumex
        Stop 'error in sum1'
      End If
    End If

!   if not first

  Else

!   w = weights for flux integral

    w = wp(jp)
    w1 = wp1(jp)

  End If

! -------------------------------------------------------------------------

! calculation of continuum quantities (w.r.t. CMF-freq) on micro-grid
! note: so far, continua 'only' corrected for micro-clumping

! interpolation in frequency and space;
! frequential interpolation only at macro-grid points with respect to
! CMF-freq.
! spatial interpolation to micro-grid only from freq.-interpolated macro-grid
! variables
! more exact method (commented out) also tested, almost identical results

! note: thus far, inversion treated in approx. way (by using abs)
! thus: inverted_range=1 needs no special attention and can be calculated
! in parallel with inverted_range=0

! for both core and non-core rays
! front hemisphere
  enerobs = 1.D8/xlamobs

  kkstart = 1

  j2 = 1
  enercmf = enerobs*(1.D0-z(j2)*v(j2)*vc/r(j2))

  Do kk = kkstart, ifre - 1
    If (fre(kk)<=enercmf .And. fre(kk+1)>enercmf) Go To 110
  End Do
  Stop ' enercmf not found (1) in fre (prepray_range)'

110 kkstart = kk
  q = (enercmf-fre(kk))/(fre(kk+1)-fre(kk))
  q1 = 1.D0 - q
  opacon(j2) = q1*opac(1, kk) + q*opac(1, kk+1)
  scont(j2) = q1*scon(1, kk) + q*scon(1, kk+1)

iloop: Do i = 2, nd1
    j1 = index1(i-1)
    j2 = index1(i)

!   for interpolation below, opacon(j1) and opacon(j2) need to be known
    If (j2<=lmax) Then
      enercmf = enerobs*(1.D0-z(j2)*v(j2)*vc/r(j2))
    Else
      enercmf = enerobs !           approximate assuming z=0 (for .not. optz0)
    End If

    Do kk = kkstart, ifre - 1
      If (fre(kk)<=enercmf .And. fre(kk+1)>enercmf) Go To 120
    End Do
    Stop ' enercmf not found (2) in fre (prepray_range)'

!   20   kkstart=1! since next loop starts !(more exact method, commented out)
120 kkstart = kk
    q = (enercmf-fre(kk))/(fre(kk+1)-fre(kk))
    q1 = 1.D0 - q
!   interpolation w.r.t. index i, even if j2 out of range
!   (accounting for radial interpolation in next loop)
    opacon_last = q1*opac(i, kk) + q*opac(i, kk+1)
    scont_last = q1*scon(i, kk) + q*scon(i, kk+1)
    If (j2<=lmax) Then
      opacon(j2) = opacon_last
      scont(j2) = scont_last
!     otherwise, index j2 out of range
    End If

jloop: Do j = j1 + 1, j2 - 1

      If (j>lmax) Exit iloop !      for non-core rays

!     more exact method commented out
!     enercmf=enerobs*(1.d0-z(j)*v(j)*vc/r(j))

!     do kk=kkstart,ifre-1
!     print*,'two',kk,enercmf,fre(kk),fre(kk+1)
!     if(fre(kk).le.enercmf .and. fre(kk+1).gt.enercmf) goto 30
!     enddo
!     stop ' enercmf not found (2) in fre (prepray)'

!     30       kkstart=kk
!     q=(enercmf-fre(kk))/(fre(kk+1)-fre(kk))
!     tl1=(1.d0-q)*opac(i-1,kk)+q*opac(i-1,kk+1)
!     tl2=(1.d0-q)*opac(i,kk)+q*opac(i,kk+1)
      factor = weight_range(j)
      If (factor==0.D0) Stop ' problems with factor (prepray_range)'
      tl1 = opacon(j1)
      tl2 = opacon_last
      opacon(j) = tl1 + factor*(tl2-tl1)
      tl1 = scont(j1)
      tl2 = scont_last
      scont(j) = tl1 + factor*(tl2-tl1)
    End Do jloop
  End Do iloop

  lmax1 = lmax
  ltot1 = lmax

! back hemisphere for non-core rays
  If (.Not. core) Then

    optz0 = .False.
    If (z(lmax)==0.D0) optz0 = .True.

    If (.Not. optz0) Then
!     assuming constant opacities/source-functions between z(lmax) and z=0.
!     (consistent with standard approach)
      lmax1 = lmax + 1
      opacon(lmax1) = opacon(lmax)
      scont(lmax1) = scont(lmax)
      If (zcon(lmax1)/=0.) Stop ' problem with zcon = 0!'
    End If

    ltot1 = 2*lmax1 - 1

    kkstart = ifre

    j2 = 1
    jj2 = 2*lmax1 - j2
    enercmf = enerobs*(1.D0+z(j2)*v(j2)*vc/r(j2)) ! negative mu gives plus
!   sign

    Do kk = kkstart, 2, -1
      If (fre(kk)>=enercmf .And. fre(kk-1)<enercmf) Go To 130
    End Do
    Stop ' enercmf not found (3) in fre (prepray_range)'

130 kkstart = kk
    q = (enercmf-fre(kk))/(fre(kk-1)-fre(kk))
    q1 = 1.D0 - q
    opacon(jj2) = q1*opac(1, kk) + q*opac(1, kk-1)
    scont(jj2) = q1*scon(1, kk) + q*scon(1, kk-1)

iloop1: Do i = 2, nd1
      j1 = index1(i-1)
      j2 = index1(i)
      jj1 = 2*lmax1 - j1
      jj2 = 2*lmax1 - j2

!     for interpolation below, opacon(j1) and opacon(j2) need to be known
      If (j2<=lmax) Then
        enercmf = enerobs*(1.D0+z(j2)*v(j2)*vc/r(j2)) ! negative mu gives plus
!       sign
      Else
        enercmf = enerobs !         approximate assuming z=0 (for .not. optz0)
      End If

      Do kk = kkstart, 2, -1
        If (fre(kk)>=enercmf .And. fre(kk-1)<enercmf) Go To 140
      End Do
      Stop ' enercmf not found (4) in fre (prepray)_range'

140   kkstart = kk
      q = (enercmf-fre(kk))/(fre(kk-1)-fre(kk))
      q1 = 1.D0 - q
!     interpolation w.r.t. index i, even if out of range
!     (accounting for radial interpolation in next loop)
      opacon_last = q1*opac(i, kk) + q*opac(i, kk-1)
      scont_last = q1*scon(i, kk) + q*scon(i, kk-1)
      If (j2<=lmax) Then
!       consistency test
        If (j2==lmax .And. optz0) Then
          If (jj2/=j2) Stop ' index problems at z=0 (prepray_range)'
          If (abs(opacon_last-opacon(j2))>1.D-10) Stop &
            ' problems at z=0 (prepray_range)'
        End If
        opacon(jj2) = opacon_last
        scont(jj2) = scont_last
!       otherwise, index j2 out of range
      End If

jloop1: Do j = j1 + 1, j2 - 1

        jj = 2*lmax1 - j
        If (j>lmax) Exit iloop1

        factor = weight_range(j)
        tl1 = opacon(jj1)
        tl2 = opacon_last
        opacon(jj) = tl1 + factor*(tl2-tl1)
        tl1 = scont(jj1)
        tl2 = scont_last
        scont(jj) = tl1 + factor*(tl2-tl1)
      End Do jloop1
    End Do iloop1

  End If !                          non core

  opacon(1:ltot1) = exp(opacon(1:ltot1))*sr
  scont(1:ltot1) = exp(scont(1:ltot1))

! for tests
! do l=1,ltot1
! if(k.eq.1) print*,jp,lmax,lmax1,l,scont(l)
! enddo

! -------------------------------------------------------------------------

! now, for all frequencies and angles, we set up opacray and sconray,
! in analogy to procedure from above [if(first) ...]

! -----no interpolation, when number of grid points (zray vs. zcon) equal

  lm1 = ltot(jp)
  If (lm1==ltot1) Then
!   typically, until jp = 30 ... 40

    opacray(1:lm1) = opacon(1:ltot1)
    sconray(1:lm1) = scont(1:ltot1)

  Else
    Do l = 1, lm1
      i = indexz(l)
      q = qarr(l)
      If (i>0.) Then
!       consistency tests, can be skipped later
        If (zray(l)/=zcon(i)) Stop ' error(1) in indexz (prepray_range)'
        If (q/=1.D0) Stop ' error(1) in qarr (prepray_range)'
        opacray(l) = opacon(i)
        sconray(l) = scont(i)
      Else
        i = abs(i)
!       consistency test, can be skipped later
        If (zray(l)<zcon(i-1) .And. zray(l)>zcon(i)) Then
          Continue
        Else
          Stop ' error(2) in indexz (prepray_range)'
        End If
!       finally, the linear interpolation
        q1 = 1.D0 - q
        opacray(l) = q*opacon(i-1) + q1*opacon(i)
        sconray(l) = q*scont(i-1) + q1*scont(i)
      End If
    End Do

  End If

! for tests
! if(k.eq.10) then
! do l=1,lm1
! write(*,fmt='(3(i4,2x),4(F11.3,2x),2(E10.4,2x))')  &
! &
! jp,lm1,l,zray(l),xcmfp(l),tempray(l),log10(xneray(l)),opacray(l),sconray(l)
! enddo
! print*
! endif

! now the most time-consuming loop to calculate the total opacity and
! source-function for a given obs. frame freq. and ray


  p2 = p(jp)**2
! some final consistency tests

  If (core) Then
    rmin = sqrt(zray(lm1)**2+p2)
    If (abs(rmin-1.D0)>1.D-15) Then
      Print *, jp, rmin
      Stop ' rmin(core) ne 1'
    End If
  Else
    lmz0 = (lm1+1)/2.
    If (zray(lmz0)/=0.D0) Then
      Print *, jp, zray(lmz0)
      Stop ' z(lmz0) ne 0'
    End If
  End If


! continuum values (not symmetric)
  opalray(1:lm1) = opacray(1:lm1)
  sray(1:lm1) = opacray(1:lm1)*sconray(1:lm1) ! cont. emissivity
  tauray(1:lm1) = 0.D0

  If (first) Then
!   calculate interpolation weights only once per jp
    If (allocated(qstore)) Deallocate (qstore, lstore)
    Allocate (qstore(lm1), lstore(lm1))

    irstart = 1
    Do l = 1, lm1

      logr = 0.5*log(zray(l)**2+p2)
!     to avoid numerical problems (e.g., model d8v, where
!     logr-logr1(1)=8.d-16)
      If (abs(logr-logr1(1))<1.D-15) logr = logr1(1)
      If (core .And. abs(logr)<1.D-15) logr = 0.D0
      If (core) Then
        Do ll = irstart, nd1 - 1
          If (logr<=logr1(ll) .And. logr>=logr1(ll+1)) Go To 150
        End Do
        Print *, l, irstart, logr
        Print *, logr1
        Stop ' logr not found in logr1 (1)'
150     irstart = ll

      Else If (l<=lmz0) Then
        Do ll = irstart, nd1 - 1
          If (logr<=logr1(ll) .And. logr>=logr1(ll+1)) Go To 160
        End Do
        Print *, l, irstart, logr
        Print *, logr1
        Stop ' logr not found in logr1 (2)'
160     irstart = ll

      Else
        If (zray(l)>=0.) Stop ' error in zray'
        Do ll = irstart, 1, -1
!         do ll=ND1-1,1,-1
          If (logr<=logr1(ll) .And. logr>=logr1(ll+1)) Go To 170
        End Do
!       for tests
!       print*,l,irstart,logr,zray(l)
        Print *, logr1
        Stop ' logr not found in logr1 (3)'
170     irstart = ll

      End If

!     save interpolation weights
      qstore(l) = (logr-logr1(ll+1))/(logr1(ll)-logr1(ll+1))
      lstore(l) = ll
!     print*,jp,l,ll,qstore(l)
    End Do
  End If

! completely new verions, interpolation of total CMF opacities and
! emissivities [sl_lgrid = sum(sl*opal)]


  indstart = nfcmf - 1
lines_l1: Do l = 1, lm1

!   restore interpolation weights
    q = qstore(l)
    q1 = 1.D0 - q
    ll = lstore(l)

    xlamcmf = xlamobs*(1.+xcmfp(l))

!   for tests, can be left out later
    If (xlamcmf<lamgrid_cmf(1) .Or. xlamcmf>lamgrid_cmf(nfcmf)) &
      Stop ' eror in xlamcmf (subr. prepray_range)'


    Do i = indstart, 1, -1
      If (xlamcmf<=lamgrid_cmf(i+1) .And. xlamcmf>=lamgrid_cmf(i)) Then
        indstart = i
        Go To 180
      End If
    End Do
    Stop ' xlamcmf not found in lamgrid_cmf'

180 Continue

!   interpolate (log-log) coarse-grid CMF line opacities/emissivities at i and
!   i+1
!   for core rays and l=lm, ll = nd-1 and  q=0 (=> opali = opal_lgrid(nd)
!   for non-core rays and l=lm, ll = 1 and  q=1 (=> opali = opal_lgrid(1)

    opali = q*opal_lgrid(ll, i) + q1*opal_lgrid(ll+1, i)
    sli = q*sl_lgrid(ll, i) + q1*sl_lgrid(ll+1, i)

    opali1 = q*opal_lgrid(ll, i+1) + q1*opal_lgrid(ll+1, i+1)
    sli1 = q*sl_lgrid(ll, i+1) + q1*sl_lgrid(ll+1, i+1)

!   interpolate micro-grid CMF line opacities/emissivites to find final
!   quantities at 'correct' CMF frequency corresponding to considered
!   observer's frame frequency
!   finally, we use a lin (freq.) vs. log (opacity) interpolation, which
!   is more accurate than lin-lin and saves time (two times exp vs. four times
!   exp)
    ww = (xlamcmf-lamgrid_cmf(i+1))/(lamgrid_cmf(i)-lamgrid_cmf(i+1))
    ww1 = 1.D0 - ww

    sl = exp(ww*sli+ww1*sli1) !     since emissivities (>0), this works always

    If (opalmin_lgrid(i)==0. .And. opalmin_lgrid(i+1)==0.) Then
      opal = exp(ww*opali+ww1*opali1) ! standard
    Else If (opalmin_lgrid(i)/=0. .And. opalmin_lgrid(i+1)/=0.) Then
!     JO May 2017: new treatment of inversion; interpolation for abs. and
!     emission
!     part separately: here, both components inverted
      opali_part = q*opal_part_lgrid(ll, i) + q1*opal_part_lgrid(ll+1, i)
      opali1_part = q*opal_part_lgrid(ll, i+1) + q1*opal_part_lgrid(ll+1, i+1)

!     opali=exp(opali)-exp(opali_part)
!     opali1=exp(opali1)-exp(opali1_part)
!     opal=ww*opali+ww1*opali1
      opal = exp(ww*opali+ww1*opali1) ! abs. part
      opal_em = exp(ww*opali_part+ww1*opali1_part) ! em. part
      opal = (opal-opal_em) !       total opacity
!     if(opal.lt.0.) print*,l,opal,opalray(l)
    Else If (opalmin_lgrid(i)==0.) Then
!     inversion only for i+1
      opali1_part = q*opal_part_lgrid(ll, i+1) + q1*opal_part_lgrid(ll+1, i+1)
      opali = exp(opali) !          already total
      opali1 = exp(opali1) - exp(opali1_part) ! total opacity
!     linear frequency interpolation
      opal = ww*opali + ww1*opali1
    Else If (opalmin_lgrid(i+1)==0.) Then
!     inversion only for i+1
      opali_part = q*opal_part_lgrid(ll, i) + q1*opal_part_lgrid(ll+1, i)
      opali = exp(opali) - exp(opali_part) ! total opacity
      opali1 = exp(opali1) !        already total
!     linear frequency interpolation
      opal = ww*opali + ww1*opali1
    Else
      Stop ' wrong condition in opacity interpolation (subr. opal_cmfgrid)'
    End If

    opalray(l) = opalray(l) + opal ! add to continuum
    sray(l) = sray(l) + sl !        add to continuum
!   print*,l,sray(l)/bnue(xlamobs,temp(l))

!   calculate final source-function
    sray(l) = sray(l)/opalray(l)
!   for tests; can happen
!   if(abs(sray(l)).gt.10.) then
!   print*,'strong emission'
!   print*,k,jp,l,opalray(l),sray(l)
!   stop
!   sray(l)=0.
!   endif

    If (l>1) Then
      tauray(l) = tauray(l-1) + 0.5D0*(opalray(l-1)+opalray(l))*(zray(l-1)- &
        zray(l))
      If (tauray(l)>taumax) Then
!       print*,k,jp,l,lm1,' tau > taumax'
        Exit lines_l1
      End If
    End If

!   if(opalmin_lgrid(i).eq.0. .and. opalmin_lgrid(i+1).eq.0.) then
!   write(*,fmt='(f10.4,2x,i4,3(2x,e10.4))')
!   xlamobs,l,sray(l),opalray(l),tauray(l)
!   else
!   write(*,fmt='(f10.4,2x,i4,3(2x,e10.4),2x,"T")')
!   xlamobs,l,sray(l),opalray(l),tauray(l)
!   endif

  End Do lines_l1

! print*,(xlam1(index_hhe(i)),i=1,icount_hhe)

  Return
End Subroutine

!-----------------------------------------------------------------------

Subroutine formacon_range(taucon, ltaumax, fscon, opacray, sconray, zray, &
  ltotal)

  Use :: nlte_type
  Use :: nlte_dim
  Use :: formalsol_add_var, Only: taumax

  Implicit None
! ---- calculates the pure continuum formal solution for each ray,
! ---- using the opacities and source-functions evaluated at the
! ---- corresponding CMF frequencies in subroutine preformal_range
! ---- in case, a linear expansion of source function is performed over each
! ---- subinterval
! -----

! .. parameters ..
  Integer (i4b), Parameter :: ndm = id_depfi, lto1 = 2*ndm
! ..
! .. scalar arguments ..
  Integer (i4b) :: ltaumax, ltotal
! ..
! .. array arguments ..
  Real (dp), Dimension (lto1) :: fscon, taucon, opacray, sconray, zray
! ..
! .. local scalars ..
  Real (dp) :: dt, e0, e1, w1, w2
  Integer (i4b) :: i
! ..

! in contrast to standard approach, opacray and sconray are already
! evalutated (for each ray and frequency) at the corresponding cmf-frequency

  taucon(1) = 0.D0

  Do i = 2, ltotal
    taucon(i) = taucon(i-1) + 0.5D0*(opacray(i-1)+opacray(i))*(zray(i-1)-zray( &
      i))
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

    fscon(i) = fscon(i+1) + exp(-taucon(i))*(sconray(i)*w1+sconray(i+1)*w2)
  End Do

! for tests
! w1=0.
! do i=2,ltotal
! dt=fscon(i-1)-fscon(i)
! w1=w1+dt
! write(*,fmt='(i4,4(2x,e10.4))') i,w1,dt,sconray(i),opacray(i)
! enddo

  Return

End Subroutine

!-----------------------------------------------------------------------

Subroutine obsfram_range(lm1, core, zray, opalray, sray, tauray, xicor, z1, &
  emint, emint1)

! resonance zones of tau=tauout stored in z1


  Use :: nlte_type
  Use :: nlte_dim
  Use :: formalsol_add_var, Only: taumax, tauout

  Implicit None

! .. parameters ..
  Integer (i4b), Parameter :: ndm = id_depfi, lto1 = 2*ndm

! .. array arguments ..
  Real (dp), Dimension (lto1) :: zray, opalray, sray, tauray

! ..
! .. scalar arguments ..
  Real (dp) :: xicor, z1, emint, emint1
  Integer (i4b) :: lm1

  Logical :: core

! .. local scalars ..
  Real (dp) :: taured, taublu, sred, sblu, dt, aux, opared, opablu, zred, &
    zblu, bb
  Integer (i4b) :: l
  Logical :: optout

  optout = .True.
  emint = .0D0

  taured = tauray(1)
  sred = sray(1)
  opared = opalray(1)
  zred = zray(1)

  Do l = 2, lm1
    taublu = tauray(l)
    If (taublu==0.) Stop ' error in tauray'
    sblu = sray(l)
    opablu = opalray(l)
    dt = taublu - taured
!   JO May 2017: tested, works
    zblu = zray(l)

!   JO May 2017: old integration, problematic when close to (S>0)
!   or in inversion (S<0, OPA < 0).

!   if(dt.gt.1.d0) then
!   sq=sred
!   else if(dt.gt.0.) then
!   else
!   sq=0.5d0*(sred+sblu)
!   endif
!   aux=exp(-taured)*sq*expuno(dt)

!   THUS: integration not over SL, but over eta if dt < 0.01 or sred,sblu<0
!   S exp(-tau) dt => eta exp(-tau) dz corresponding to
!   eta exp(-tau) dt/(0.5(opared+opablu)) [from calculation of tau, dt in
!   subr. prepray_range]
!   Integral then solved via trapez, with ranges [0,dt] and final
!   multiplication
!   with exp(-taured). Factor 0.5 cancels with factor 0.5 from dt/(0.5*opa...)

    If (dt>0.01D0 .And. sblu>0.D0 .And. sred>0.D0) Then
!     linear exansion of S, then integration; only done for larger dt;
!     otherwise
!     cancellation effects
!     see also formalsol.f90, routine obsfram
      bb = (sblu-sred)/dt
      aux = sred*(1.D0-exp(-dt)) + bb*(1.D0-exp(-dt)*(dt+1.D0))
      aux = aux*exp(-taured)
    Else !                          small or negative dt, or inversion, using
!     trapezoidal integration
      aux = (sred*opared+sblu*opablu*exp(-dt))/(opared+opablu)*dt*exp(-taured)
!     aux=0.5d0*(sred*opared+sblu*opablu*exp(-dt))*dz*exp(-taured)
!     gives identical results
    End If
    emint = emint + aux

!   for tests
!   write(*,fmt='(i4,4(2x,e10.4))') l,emint,aux,sray(l),opalray(l)

!   if (emint.gt.10.d0) stop ' emint > 10 in obsfram_range'

    If (optout .And. taublu>=tauout) Then
      z1 = zray(l)
      optout = .False.
    End If

    If (taublu>taumax) Then
      emint1 = 0.D0
!     here we don't reach the core, thus iplus not relevant
!     if we would include the following statement, this could lead to problems
!     when the source function is still low when taumax is reached.
!     in this case, we could have emint1 >> emint
!     THUS, don't active the following statement (unless for tests)
!     if (core) emint1 = xicor*exp(-taublu)
      Return
    End If

    taured = taublu
    sred = sblu
    opared = opablu
    zred = zblu

  End Do

  emint1 = 0.D0
! here we reached the core, thus adding iplus
  If (core) emint1 = xicor*exp(-taublu)


! if finished, test with pure cont, and compare to formacon
  Return
End Subroutine

!-----------------------------------------------------------------------

Subroutine calcxmax_range(temp, xne)

! calculated maximum width/velocity until which Stark etc. broadening is
! required. We use here the same condition as in calcxmax. i.e.
! requiring that opal*perf(dlam)/opac < eps_broad (e.g., 5.d-3)

  Use :: nlte_type
  Use :: nlte_dim
  Use :: fund_const, Only: clight
  Use :: formalsol_var, Only: nstark
  Use :: formalsol_add_var, Only: xlam1, icount_hhe, index_hhe, eps_broad, &
    dlam_max_stark, vmax_broad, vmax_broad_safety, l23, beta_p

  Implicit None

  Integer (i4b), Parameter :: nd = id_ndept

  Real (dp), Dimension (nd), Intent (In) :: temp, xne

  Integer (i4b) :: i, j, l
  Real (dp) :: dlam, perf, vbroad, aux
  Real (dp) :: perfil2
  Logical :: flag

  vbroad = 0.

  Do i = 1, icount_hhe
    If (nstark(i)==0) Then
      dlam_max_stark(i) = 0.
    Else
      l = l23(i)
      dlam = 2.
      aux = eps_broad/beta_p(i)
      perf = perfil2(i, 0.D0, dlam, temp(l), xne(l), flag)
      Do While (perf>aux)
        dlam = dlam + 2.
        perf = perfil2(i, 0.D0, dlam, temp(l), xne(l), flag)
      End Do
      dlam = dlam + 2. !            safety
      dlam_max_stark(i) = dlam
      j = index_hhe(i)
      vbroad = max(vbroad, dlam/xlam1(j)*clight)
!     print*,id1(j),dlam,dlam/xlam1(j)*clight/1.d5
    End If
  End Do

  vmax_broad = vbroad !             for prepray_range
  Print *
  Print *, ' Maximum extent of Stark-broadened lines = ', vmax_broad/1.D5, &
    ' km/s'
  If (vmax_broad>vmax_broad_safety) Then
    Print *, ' WARNING!!!'
    Print *, ' vmax_broad > vmax_broad_safety = ', vmax_broad_safety/1.D5
    Print *, ' might lead to problems at freq. boundaries'
  End If

  Return
End Subroutine

!-----------------------------------------------------------------------

Function perfil2(i, xlamb0, xlam, tem, xne, flag)

! no action concerning clumping required, since we solve inside clumps,
! and n_e is the increased value

  Use :: nlte_type
  Use :: nlte_dim
  Use :: formalsol_var, Only: nstark, nws, nts, nes, qhalf, dws, ts, es, ps

  Implicit None

! ----calculates Stark profiles for the 'i'
! ----component, at wavelength xlam, with temperature and electron density
! ----'tem' and 'xne' resp.
! In this version, vturb = const has been folded (PREFORMAL) into all
! tabulated profiles (Stark, Iso and Griem)

! ..
! .. scalar arguments ..
  Real (dp) :: xlamb0, xlam, tem, xne, perfil2
  Integer (i4b) :: i
  Logical :: flag
! ..
! .. local scalars ..
  Real (dp) :: dw, elog, pf, pff, pfl, pl, plf, pll, tlog, wef, wtf, wwf
  Integer (i4b) :: ifi, j, jeb, jtb, jwb, l1, l2, l3, l4
! ..

  flag = .False.
  If (nstark(i)/=1) Stop ' perfil2 erroneously called'

! --- stark profile (doppler broadening -- temp+stark -- is included)

! --- weights for wavelength

  dw = xlam - xlamb0 !              here in angstrom
  If (qhalf(i)) dw = abs(dw)

! --- out of freq. range in table (assuming then phi = 0.)

  If (dw<dws(1,i) .Or. dw>dws(nws(i),i)) Then
    flag = .True.
    perfil2 = 0.D0
    Return
  End If

  Do j = 2, nws(i)
    ifi = j
    If (dws(j,i)>=dw) Go To 100
  End Do

  Stop ' error in table freq. grid'

100 Continue

  jwb = ifi - 1
  wwf = (dw-dws(ifi-1,i))/(dws(ifi,i)-dws(ifi-1,i))

! --- weights for t

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
!   if(xnemin.ne.0.) then
    Print *, elog, es(1, i)
    Stop ' error in xnemin or stark profile tables'
!   endif
!   elog=es(1,i)  ! for approximate treatment (just to calculate xmax)
  End If

  Do j = 2, nes(i)
    ifi = j
    If (es(j,i)>=elog) Exit
  End Do

  jeb = ifi - 1
  wef = (elog-es(ifi-1,i))/(es(ifi,i)-es(ifi-1,i))

! --- interpolation

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

  perfil2 = 10.D0**(wef*(pf-pl)+pl)

  Return
End Function

!-----------------------------------------------------------------------

Subroutine dealloc

  Use :: formalsol_add_var, Only: id1, xlam1, gf1, k_range, inverted_range, &
    opal_range, opal_plus_range, sline_range, opal_lgrid, opal_part_lgrid, &
    sl_lgrid, lamgrid, lamgrid_cmf, index_hhe, l23, beta_p, opalmin_lgrid

  Implicit None

  Deallocate (id1, xlam1, gf1)
  Deallocate (k_range, inverted_range)
  Deallocate (index_hhe, l23, beta_p)
  Deallocate (opal_range, opal_plus_range, sline_range)
  Deallocate (opal_lgrid, opal_part_lgrid, sl_lgrid, opalmin_lgrid)
  Deallocate (lamgrid, lamgrid_cmf)

  Return
End Subroutine
