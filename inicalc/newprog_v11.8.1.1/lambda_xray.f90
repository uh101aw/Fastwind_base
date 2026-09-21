!routines for x-ray treatment

!Basically, for estimated shock temperatures,
!the Raymond & Smith 1989 code is used to calculate the cooling functions
!(shock emission). Moreover, data for inner-shell absorption and Auger
!ionisation are imported and modified for further use.

!version 1.0   October 2014, programmed by Luiz Carneiro and Jo
!              following Pauldrach+ 1994 and Feldmeier+ 1997

!version 1.1   July 2016, programmed by Jo and Koushik Sen, to
!allow for a better description of radiative shocks
!(see version control in nlte.f90), following Owocki+ 2013
!(see also Puls+ 2020). To switch to that formulation,
!the input variable px needs to be < -999.99

!version 1.2   few changes (mostly, regarding comments), and
!Sept 2020: one bug found by Sarah Brands removed

!version 1.3   compatible with gfortran

!version 1.4   June 2025 obsolet fortran features removed
!(mostly: arithmetic if and do-loops without end do)

!version 1.5 Nov. 2025: polished

!-----------------------------------------------------------------------

Subroutine tshock(ts, v, r, rtau23, xmu, nd)
! calculates stratification of X-ray shock-temperatures, basically following
! Pauldrach et al. 1994.
! Change this routine if you want to test other stratifications

! NOTE: The minimum emission radius is set from the minimum of
! the velocity (V_em > MX*VSOUND) or radius constraint (R_em > RMINX)

  Use :: nlte_type
  Use :: nlte_dim, Only: id_ndept
  Use :: fund_const, Only: amh, akb
  Use :: nlte_var, Only: modnam, vmax, vsound
  Use :: nlte_xrays, Only: fx, px, lxmin, rminx_eff

  Implicit None

  Integer (i4b), Intent (In) :: nd
  Real (dp), Dimension (nd), Intent (In) :: v, r
  Real (dp), Dimension (nd), Intent (Out) :: ts

  Real (dp), Dimension (id_ndept) :: ujump

  Real (dp) :: rtau23, xmu
  Real (dp) :: gamx, mx, uinfx, ts_const, vmin, rminx

  Integer (i4b) :: l, l1, l2

  Logical :: optxray

  Open (1, File=trim(modnam)//'/INXRAY.DAT', Status='OLD')
  Rewind 1
  Read (1, Fmt=*) optxray, fx
  If (.Not. optxray) Then
    Write (999, *) ' STOP:  OPTXRAY = FALSE IN CATALOGUE/INXRAY.DAT'
    Stop ' OPTXRAY = FALSE IN CATALOGUE/INXRAY.DAT'
  End If
  Read (1, Fmt=*) gamx, mx, rminx, uinfx, px
  Close (1)

! for simplicity, we calculate the shock temperature from the approximate
! expression, neglecting terms with O(vsound), i.e., assuming strong shocks
! (also w.r.t. consistency with Feldmeier et al., Eq. 19), and approximate
! v_0=v_pre-v_shock by v_0 approx v_jump=u (assuming an isothermal shock
! w.r.t. v_1=v_post_v_shock << v_0, see Pauldrach et al. 1994)
! Note: the exact expression would be
! Tshock = 3/16 mu mh/k * (u^2 + (14/5 a^2 * (1 - 3/14 a^2/u^2)))
! again using u instead of v_0

  ts_const = 3.D0/16.D0*xmu*amh/akb

  ujump = uinfx*1.D5*v**gamx

! now assuming delmu = 1 (to ensure consistent restart)
  ts = ts_const*ujump**2
! TS=TS_CONST/DELMU*UJUMP**2
! corrected for variable xmu; xmu_actual=xmu/delmu

! print*, 'XMU=',XMU
! print*, 'AMH=',AMH
! print*, 'AKB=',AKB
! print*, 'V= ',V
! print*, 'GAMX =',GAMX
! print*, 'UJUMP=',UJUMP
! print*, 'TS=',TS


  vmin = mx*vsound

  Do l = 1, nd
    If (v(l)*vmax<vmin) Exit
  End Do
  l1 = l

! JO changed July 2016
! RMINX in units of RSTAR (not SR)
  Do l = 1, nd
    If (r(l)/rtau23<rminx) Exit
  End Do
  l2 = l

  lxmin = max(l1, l2) - 1
  rminx_eff = r(lxmin)

  Do l = lxmin + 1, nd
    ts(l) = 0.D0
  End Do

  Write (999, *)
  Write (999, *) ' SHOCK TEMPERATURES CALCULATED, FROM L=1 TO L=', lxmin
  Write (999, *) ' MINIMUM EMISSION RADIUS (in units of RSTAR) = ', &
    rminx_eff/rtau23
  Write (999, *) ' MINIMUM EMISSION RADIUS (in units of SR) = ', rminx_eff
  Write (999, *)
  Write (999, *) ' SHOCK TEMPERATURE AT RMAX = ', ts(1)/1.D6, ' MK'
  Write (999, *) ' MINIMUM SHOCK TEMPERATURE = ', ts(lxmin)/1.D6, ' MK'

  Print *
  Print *, 'SHOCK TEMPERATURES CALCULATED, FROM L=1 TO L=', lxmin
  Print *, 'MINIMUM EMISSION RADIUS (in units of RSTAR) = ', rminx_eff/rtau23
  Print *, 'MINIMUM EMISSION RADIUS (in units of SR) = ', rminx_eff
  Print *
  Print *, 'SHOCK TEMPERATURE AT RMAX = ', ts(1)/1.D6, ' MK'
  Print *, 'MINIMUM SHOCK TEMPERATURE = ', ts(lxmin)/1.D6, ' MK'

  Return
End Subroutine

!---------------------------------------------------------------------------

Subroutine lambda_xray(xmu, ts, r, v, rho, xne, xnh, fcl, nd)

! calculate X-ray cooling functions accounting for cooling zones, both for
! radiative and adiabatic shocks, following Feldmeier et al. 1997
! The basic cooling functions are calculated via Raymond & Smith 1989

! Note: you might play with theta and xne;
! so far, xne = 1.d10 fixed in cooling function
! JO July 2016
! for the electron density used to calculate the cooling time/length, we
! now use the hot-plasma conditions (IHe=2), but check for exceptions

! version 0.1, programmed by Luiz Carneiro (2014)
! version 0.2, modified by Jo (using input from Koushik Sen) (July 2016)

  Use :: nlte_type
  Use :: fund_const, Only: hh, clight, ev
  Use :: nlte_var, Only: yhe => yhein, modnam
  Use :: nlte_xrays, Only: enerx, spec => lambdax, lxmin, fxl, fx, px, &
    rminx_eff

  Implicit None

  Integer (i4b) :: nd
  Real (dp), Dimension (nd), Intent (In) :: ts, r, v, rho, xne, xnh, fcl
  Real (dp) :: xmu

! input for rs-code: variables/data
  Integer (i4b), Parameter :: nf = 1000 ! <= 1000, because required by
! Raymond-Smith

! parameters which may be changed
  Integer (i4b) :: ntdec = 25 !     number of t-points per decade
! if you change the binsize, change also nf
  Real (dp) :: x_binsyz = 2.5D0 !   (eV)


  Integer (i4b) :: i_nbin, i_nbln = 0, i_num = 12, i_iprint = 0, i_jprint = 0, &
    i_iphot = 0, i_idens = 0, i_icx = 0

  Real (dp) :: x_binmin = 1., x_blnmin = 0., x_blnsyz = 0.

  Real (dp) :: logt, tinput

! ***************************************************************************
! NOIE: so far and for computational efficiency,
! we calculate Lambda with constant electron density ne=10^10.
! Tests have shown that this is basically OK, but certain frequencies
! around 10 to 100 ev react weakly on variations
  Real (dp) :: logne = 10.D0

! NOIE: so far we assume theta = T_rev/T_for (post) as unity.
! This is just a guess and needs to be tested
  Real (dp) :: theta = 1.D0
  Real (dp) :: wr = 0.5D0, wf = 0.5D0

  Logical :: isothermal = .False.
! if true, isothermal shocks will be calculated (check that wr+wf=1.)

! regarding the adiabatic double-shock structure:
! if theta = 1, the structure of the cooling zones behind the forward and
! the reverse shocks are identical (but somewhat different from the (single)
! radiative one). Since we now have two shocks per volume, this would require
! a different filling factor (factor two lower), and thus we divide the
! total emission by a factor of two (using weighting factors of 0.5).
! If we would use wr=1 and wf=0, we would only count the contribution from
! the reverse shock (for theta=1 this is identical to 0.5/0.5), and vice versa
! for the forward shock with wr=0 and wf=1.
! Note that reverse shocks are destroyed around 10 Rstar.
! Realistic weighting factors need to be derived from hydro-sims.
! We suggest to keep the above parameterization until further evidence.

! NOTE: As in Feldmeier, we assume that the parameterized
! jump velocity always refers to the reverse-shock temperature
! in the radiative case, but to the forward-shock temperature
! in the adiabatic case. The reason for doing so is that the forward-shock
! structure in the latter is independent of theta, whilst the density of
! the reverse shock changes significantly. Thus, normalizing to the
! forward shock retains a somewhat stable solution if theta is changed

! Note also the inconsistency that we have calculated Tshock from u_jump
! (assuming v_0 approx u_jump), but in the following will use
! v_1 = v_0/4 (adiabatic shock). Nevertheless, this inconsistency is weak,
! since globally most shocks are (quasi-) isothermal with rho_1/rho_0 >> 1
! which controls the shock temperature, whereas locally (close to the jump
! itself), the jump is only a factor of four. Anyway, since the difference is
! only a factor (4/3)^2, there should be no problem here.

! ***************************************************************************


  Integer (i4b), Parameter :: nfeld = 1000 ! number of integration points over
! cooling zone

  Integer (i4b) :: i, k, nt, l

  Real (sp), Dimension (nf) :: energy, lambdanu

  Real (dp) :: tmin, tmax, dt, w1, w2, t1, t2, tt
  Real (dp) :: fract_tc_tf, fract_lc_rs, xne_av, xnh_av, kappa_f, dxi, &
    xne_app_last

  Real (dp), Allocatable, Dimension (:) :: tgrid, lambdainter
  Real (dp), Allocatable, Dimension (:, :) :: lambda

! Declaration of variables in the subroutines

  Real (dp), Dimension (nfeld) :: f_ksi, g_ksi

  Real (dp), Dimension (nfeld) :: f_eta_rev, g_eta_rev, f_eta_for, g_eta_for

  Real (dp), Dimension (nfeld) :: f_rev, g_rev


! JO July 2016
! check consistency of provided ('cool') electron density at lowermost point
! with approximate hot plasma xne-value. If differences, then most likely
! other elements than H/He are dominating, and another approach must be
! provided.
  xne_app_last = xnh(nd)*(1.+2.*yhe)
  If (abs(1.-xne_app_last/xne(nd))>0.05) Then
    Write (999, *) &
      ' approximate electron-density used for x-ray emission erroneous!'
    Write (999, *) ' most likely, h/he not dominating in this model!'
    Write (999, *) ' change approach!'
    Write (999, *) xne(nd), xne_app_last
    Write (999, *) &
      ' STOP: approximate electron-density used for x-ray emission erroneous!'
    Print *, &
      ' approximate electron-density used for X-ray emission erroneous!'
    Print *, ' most likely, H/He not dominating in this model!'
    Print *, ' change approach!'
    Print *, xne(nd), xne_app_last
    Stop ' approximate electron-density used for X-ray emission erroneous!'
  End If


! create grid of models

  tmin = 4.D0
  If (tmin>ts(lxmin)) Then
    Write (999, *) ' STOP: tmin > ts(lxmin). Modify minimum emission radius!'
    Stop ' tmin > ts(lxmin). Modify minimum emission radius!'
  End If

  tmax = max(ts(1), ts(1)*theta) !  to account for the double shock structure
  tmax = int(log10(tmax)) + 1.
  Write (999, *)
  Write (999, *) ' minimum shock temp. (log10) for grid = ', tmin
  Write (999, *) ' maximum shock temp. (log10) for grid = ', tmax
  Write (999, *) ' number of T-points per decade = ', ntdec
  Print *
  Print *, ' minimum shock temp. (log10) for grid = ', tmin
  Print *, ' maximum shock temp. (log10) for grid = ', tmax
  Print *, ' number of T-points per decade = ', ntdec

  dt = tmax - tmin
  nt = dt*ntdec + 2
  dt = 1.D0/float(ntdec)

  Allocate (tgrid(nt))

  Do i = 1, nt
    tgrid(i) = tmin + (i-1)*dt
!   print*,i,tgrid(i)
  End Do

  If (tgrid(nt)<=tmax) Then
    Write (999, *) ' STOP: tgrid(nt) < tmax'
    Stop ' tgrid(nt) < tmax'
  End If

! actually i_nbin can be also smaller than nf, if other freq. grid required,
! e.g.
! i_nbin = int(1000/x_binsyz),
! which would give 500 points for x_binsyz = 2.)

  i_nbin = nf

! print*,nt
  Allocate (lambda(i_nbin,nt))
  Allocate (lambdainter(i_nbin))

  Do i = 1, nt

!   note: RS-code in single precision
    Call raymond_smith(sngl(tgrid(i)), sngl(logne), energy, lambdanu, i_nbin, &
      sngl(x_binmin), sngl(x_binsyz), i_nbln, sngl(x_blnmin), sngl(x_blnsyz), &
      i_num, i_iprint, i_jprint, i_iphot, i_idens, i_icx)
!   ********
!   ********Added this block to make lambdainter = 0
!   ********in case of lamda = 0    05/08/14
!   ********
    Where (lambdanu/=0.D0)
      lambda(:, i) = log10(lambdanu(:))
    Elsewhere
      lambda(:, i) = 0.D0
    End Where
  End Do

! independent of local conditions
  Call radiative_part(g_ksi, f_ksi, nfeld)

  If (isothermal) Then
    g_ksi = 1.D0
    f_ksi = 1.D0
  End If

  Allocate (spec(nf,lxmin))
  spec = 0.D0

  Print *

  If (.Not. allocated(fxl)) Allocate (fxl(nd))
  fxl = 0.D0

  Open (1, File=trim(modnam)//'/OUT_XRAY', Status='unknown')
  Rewind 1

  Print *
  If (px<-999.99) Then
    Write (999, *) ' Old approach with constant filling factor,'
    Write (999, *) ' giving rise to rho^2-dependent emission'
    Write (999, *) ' for radiative and adiabatic shocks!'
    Write (999, *) ' For details, check file OUT_XRAY'
    Write (999, *)
    Print *, ' Old approach with constant filling factor,'
    Print *, ' giving rise to rho^2-dependent emission'
    Print *, ' for radiative and adiabatic shocks!'
    Write (1, *) ' Feldmeier approach with constant filling factor'
    Write (1, *) ' fx = ', fx
    Write (1, *) ' Rmin (in SR) = ', rminx_eff
  Else
    Write (999, *) ' New approach with radius-dependent filling factor,'
    Write (999, *) &
      ' giving rise to rho-dependent emission for radiative shocks,'
    Write (999, *) &
      '          and rho^2-dependent emission for adiabatic ones!'
    Print *, ' New approach with radius-dependent filling factor,'
    Print *, ' giving rise to rho-dependent emission for radiative shocks,'
    Print *, '          and rho^2-dependent emission for adiabatic ones!'
    Write (999, *) ' For details, check file OUT_XRAY'
    Write (999, *)
    Write (1, *) ' New approach with radius dependent filling factor'
    Write (1, *) ' no = dN/dln r = ', fx, ' p = ', px
    Write (1, *) ' Rmin (in SR) = ', rminx_eff
  End If
  Write (999, *)
  Print *

  Do l = 1, lxmin

    logt = log10(ts(l))

    tinput = logt

!   xne_av=xne(l)/fcl(l)
    xnh_av = xnh(l)/fcl(l)
!   modified
    xne_av = xnh_av*(1.+2.*yhe)

!   &
!   (Fract_tc_tf,Fract_Lc_rs,ts(l),r(l),v(l),rho(l),xne_av,xnh_av,xmu/delmu(l),
!   &
!   now assuming delmu=1 to ensure consistent restart
    Call cooling_flowing_time(fract_tc_tf, fract_lc_rs, ts(l), r(l), v(l), &
      rho(l), xne_av, xnh_av, xmu, theta, kappa_f, l)

    If (fract_tc_tf<1.D0) Then
      If (kappa_f/=0.D0) Then
        Write (999, *) ' STOP: radiative and kappa_f ne 0!'
        Stop ' radiative and kappa_f ne 0!'
      End If
      g_rev = g_ksi
      f_rev = f_ksi

    Else
      Call adiabatic_part(theta, kappa_f, f_eta_rev, g_eta_rev, f_eta_for, &
        g_eta_for, nfeld)

      If (isothermal) Then
        f_eta_rev = 1.D0
        g_eta_rev = 1.D0
        f_eta_for = 1.D0
        g_eta_for = 1.D0
      End If

!     g normalized to forward shock, g(eta_f=1)=1
!     rev. shock temperature = : T_for * g_rev =: T_jump * g_rev with T_jump =
!     ts
      g_rev = g_eta_rev
      f_rev = f_eta_rev

!     do i = 1,NFELD
!     print*, i,f_eta_rev(i),g_eta_rev(i),f_eta_for(i),g_eta_for(i)
!     end do
    End If

    Do k = 1, nfeld !               radiative cooling
      logt = tinput

      logt = logt + log10(g_rev(k))

!     The emission of temperature log(ts) < 4.0 is too low and has no big
!     contribution.
      If (logt<4.0D0) Then
        lambdainter = 0.D0
        Go To 110
      End If

!     interpolate the cooling function at logt from the grid
      Do i = nt, 1, -1
        If (logt<=tgrid(i) .And. logt>tgrid(i-1)) Go To 100
      End Do

      Write (999, *) logt
      Write (999, *) tgrid
      Write (999, *) ' logt not found in tgrid - 1'
      Print *, logt
      Print *, tgrid
      Print *, ' logt not found in tgrid - 1'
      Stop


100   t1 = 1.D0/10.D0**tgrid(i-1)
      t2 = 1.D0/10.D0**tgrid(i)
      tt = 1.D0/10.D0**logt

      w1 = (t2-tt)/(t2-t1)
      w2 = 1.D0 - w1

!     print*,w1,w2

!     ********
!     ********Added this block to make lambdainter = 0
!     ********in case of lamda = 0    05/08/14
!     ********
      Where (lambda(:,i-1)/=0.D0 .And. lambda(:,i)/=0.D0)
        lambdainter(:) = 10.D0**(w1*lambda(:,i-1)+w2*lambda(:,i))
      Elsewhere
        lambdainter(:) = 0.D0
      End Where

110   Continue
!     integration weights assumed to be equal (neglecting boudary effects)
      If (fract_tc_tf<1.D0) Then
!       no weight, one shock
        spec(:, l) = spec(:, l) + f_rev(k)**2*lambdainter(:)
      Else

!       contrib. from reverse shock
!       two shocks, weighted
        spec(:, l) = spec(:, l) + wr*f_rev(k)**2*lambdainter(:)


!       ***************************************************************************
!       contrib. from forward shock
        logt = tinput

!       g normalized to forward shock, g(eta_f=1)=1
!       forw. shock temperature = : T_for * g_eta_for =: T_jump * g_eta_for
!       with T_jump = ts
        logt = logt + log10(g_eta_for(k))

        If (logt<4.0D0) Then
          lambdainter = 0.D0
          Go To 130
        End If

        Do i = nt, 1, -1
          If (logt<=tgrid(i) .And. logt>tgrid(i-1)) Go To 120
        End Do

        Write (999, *) ' logt not found in tgrid - 2'
        Print *, ' logt not found in tgrid - 2'
        Stop


120     t1 = 1.D0/10.D0**tgrid(i-1)
        t2 = 1.D0/10.D0**tgrid(i)
        tt = 1.D0/10.D0**logt

        w1 = (t2-tt)/(t2-t1)
        w2 = 1.D0 - w1

!       print*,w1,w2

!       ********
!       ********Added this block to make lambdainter = 0
!       ********in case of lamda = 0    05/08/14
!       ********

        Where (lambda(:,i-1)/=0.D0 .And. lambda(:,i)/=0.D0)
          lambdainter(:) = 10.D0**(w1*lambda(:,i-1)+w2*lambda(:,i))
        Elsewhere
          lambdainter(:) = 0.D0
        End Where

130     Continue

!       two shocks, weighted
        spec(:, l) = spec(:, l) + wf*f_eta_for(k)**2*lambdainter(:)

      End If
!     ***************************************************************************

    End Do !                        cooling zones

  End Do !                          depth loop

  Close (1)

  dxi = 1./float(nfeld)

  spec = spec*dxi*1.D-23 !          since RS codes in these units

  Allocate (enerx(0:nf))

  enerx(0:nf-1) = energy(1:nf)
  enerx(nf) = energy(nf) + x_binsyz

! for tests
! do l=1,lxmin
! do l=9,9
! write(*,'(i4,3es12.4)') (i,enerx(i),spec(i,l),10.d0**lambda(i,66)*1.d-23,
! i=1,i_nbin)
! enddo
  Write (999, *)
  Write (999, *) ' X-ray cooling function calculated for all energies between'
  Write (999, *) energy(1), ' to ', energy(nf) + x_binsyz, ' eV'
  Write (999, *) ' with bin-width = ', x_binsyz, ' eV'
  Print *
  Print *, ' X-ray cooling function calculated for all energies between'
  Print *, energy(1), ' to ', energy(nf) + x_binsyz, ' eV'
  Print *, ' with bin-width = ', x_binsyz, ' eV'

! transform eV grid to frequency grid (in Kayser)
  dxi = hh*clight/ev
  dxi = 1.D0/dxi
  enerx = enerx*dxi

  Return

End Subroutine

!---------------------------------------------------------------------------

Block Data block_data_raymond

  Common /params/nj(12), abunj(12), abund, binmin, binsyz, nbin
  Common /pt/rf(500), tau(150), tauhe(150), tplus(150)

  Data abunj/10.93, 8.52, 7.96, 8.82, 7.92, 7.42, 7.52, 7.2, 6.9, 6.3, 7.6, &
    6.3/

  Data rf/500*1.0/

End Block Data


Subroutine raymond_smith(xlogt, xlogne, hnuvec, binout, i_nbin, x_binmin, &
  x_binsyz, i_nbln, x_blnmin, x_blnsyz, i_num, i_iprint, i_jprint, i_iphot, &
  i_idens, i_icx)

! CALLS XSPCT FOR VARIOUS PARAMETERS:
! Dec 17, 1991 : modify GAUNT5 to improve transition 1 of Li-like ions
! Sep 21, 1993 : add option to choose the solar photospheric and corona
! abundances from Anders and Grevesse (1989, Geochimica et Cosmochimica
! 53, 197) as well as the cosmic abundances from Allen 1973 (SAD)
! Note that the default abundances are Allen 1973 cosmic values.

  Use :: nlte_app, Only: abund_h12
  Use :: fastwind_params, Only: fpath_xrays

  Real, Dimension (1000) :: hnuvec, binout

  Character (1) :: irep, iab, iabsel
  Common /hln/ext(11), half, halfc, alfly, alflyc, hbet, hbetc, twop
  Common /bln/bln(1000), blnmin, blnsyz, nbln
  Common /params/nj(12), abunj(12), abund, binmin, binsyz, nbin
  Common /result/conce(30), gndrec(30), power(220), rhy, heneut, heplus, dne, &
    pcool, pou, pot, re, tu, pm(4)
  Common /pt/rf(500), tau(150), tauhe(150), tplus(150)
  Common /contin/brmev(1000), recev(1000), tufev(1000)
  Common /com/cnc(12, 30), ptot(12, 220), abin(1000), bin(1000)
  Dimension :: abcor(12), abphot(12)
  Data abcor/10.14, 7.90, 7.40, 8.30, 7.46, 7.59, 7.55, 6.93, 5.89, 6.46, &
    7.65, 6.22/, abphot/10.99, 8.56, 8.05, 8.93, 8.09, 7.58, 7.55, 7.21, 6.56, &
    6.36, 7.67, 6.25/

  Dimension :: iabindex(12)
  Data iabindex/2, 6, 7, 8, 10, 12, 14, 16, 18, 20, 26, 28/

  Logical :: start
  Data start/.True./

  External :: block_data_raymond

  Open (Unit=1, File=fpath_xrays//'file1.dat')
  Open (Unit=2, File=fpath_xrays//'atomic.dat')
  Open (Unit=3, File=fpath_xrays//'file2.dat')

! open (unit = 7, file = 'answers.dat')

  nbin = i_nbin
  binmin = x_binmin
  binsyz = x_binsyz

  nbln = i_nbln
  blnmin = x_blnmin
  blnsyz = x_blnsyz


  num = i_num
  iprint = i_iprint
  jprint = i_jprint
  iphot = i_iphot
  idens = i_idens
  icx = i_icx

  Call atread(num)
  If (iphot/=0) Then
    Print *, 'RF1?'
    Read *, rf1
    Call radrd(1)
    Do n = 1, nbin
      abin(n) = rf1*abin(n)
    End Do

  End If

100 Continue

  iab = 'N'
! this needs to be changed

  If (iab/='Y' .And. iab/='y') Go To 140
110 Continue
  Print *, 'WANT CORONAL(C), PHOTOSPHERIC(P), OR 1 by 1(1)?'
  Read (5, 160) iabsel
  If (iabsel=='C' .Or. iabsel=='c') Go To 120
  If (iabsel=='P' .Or. iabsel=='p') Go To 130
  Print *, 'NO, ABUND?'
  Read (5, *) no, abund

  abunj(no) = abund

  Print *, ' ANOTHER ABUNDANCE MODIFICATION?'
  Read (5, 160) iab
  If (iab=='Y' .Or. iab=='y') Go To 110
  Go To 140
120 Continue
  Do i = 1, 12
    abunj(i) = abcor(i)
  End Do

  Go To 140
130 Continue
  Do i = 1, 12
    abunj(i) = abphot(i)
  End Do

140 Continue

! here we overwrite our own abundances
! we do this only once, since ABUNJ in DATA statement
  If (start) Then
    Write (999, *)
    Write (999, *) ' used abundances (consistent)'
    Print *
    Print *, ' used abundances (consistent)'
    Do i = 1, 12
      ki = iabindex(i)
      xnewab = log10(abund_h12(ki)) + 12.D0
      Write (999, *) i, ki, xnewab
      Print *, i, ki, xnewab
      abunj(i) = xnewab
    End Do
    start = .False.
  End If

  t = 10.**xlogt
  dene = 10.**xlogne

! COMPUTE RATIO OF IONIZED TO NEUTRAL HYDROGEN
  cionh = cion(1, 1, 13.6, t)
  If (iphot/=0) cionh = cionh + phot(1, 1, 13.6, t, dene)/dene
  be = 157890./t
  crec = 2.06E-11*(.4288+.5*alog(be)+.469*be**(-.33333))/sqrt(t)

  rhy = cionh/crec

! write (7, * ) ' rhy = ', rhy


! for tests
! print*,nbin,binmin,binsyz
! print*,nbln,blnmin,blnsyz
! print*,num,iprint,jprint,iphot,idens,icx
! print*,t,dene



  Call feline(t, dene, num, iprint, jprint, 0, iphot, idens, icx)
! put hydrogen emission into bins
! st = sqrt(t)
! call hline(rhy,0.,t,st,be,0)
! critn = 100.*sqrt(t)
! p2ph = ext(11)*critn/(critn+dene)
! ext(9) = ext(9) + ext(11)*(1.-critn/(critn+dene))
! imax = (12399./1215.6 - binmin)/binsyz
! imax= min0(imax,nbin)
! if (imax .le. 1) go to 76
! DO 75 ik= 1,imax
! r = (binmin+binsyz*(ik-.5)) * 1215.6/12399.
! tufev(ik) = 12.*p2ph*r*r*(1.-r)*1215.6*binsyz/12399.
! print *,ik,tufev(ik)
! 75 bin(ik) = bin(ik) + tufev(ik)
! 76 continue
! ilmin = (1215.6 - blnmin) / blnsyz + 1
! DO 85 ibln = ilmin,nbln
! wbln = (ibln - 0.5)*blnsyz + blnmin
! ephot = 12399. / wbln
! ibn = (ephot-binmin) / binsyz + 1
! de = ephot - 12399. / (wbln + blnsyz)
! if (ibn .ge. 1 .and. ibn .le. nbin) bln(ibln) = bln(ibln) +
! 1  tufev(ibn) * de / binsyz
! 85   continue
  ibn = (10.2-binmin)/binsyz + 1
  ibl = 0
  If (ibn>=1 .And. ibn<=nbin) bin(ibn) = bin(ibn) + ext(6) + ext(9)
  If (blnsyz>0) ibl = (1215.6-blnmin)/blnsyz + 1

  If (ibl>=1 .And. ibl<=nbln) bln(ibl) = bln(ibl) + ext(6) + ext(9)
  If (nbin==0) Go To 150
  Write (1, 170) t, binmin, binsyz, abunj
  Write (1, 180)(bin(k), k=1, nbin)
  If (nbln==0) Go To 150
  Write (3, 170) t, blnmin, blnsyz, abunj
  Write (3, 180)(bln(k), k=1, nbln)
150 Continue

  irep = 'n'
  If (irep=='Y' .Or. irep=='y') Go To 100

  Do i = 1, nbin
    hnuvec(i) = binmin + i*binsyz - binsyz
    binout(i) = bin(i)
  End Do

  Close (1)
  Close (2)
  Close (3)



  Return


160 Format (A1)
170 Format (' TEMPERATURE=', E10.2, ' BINMIN,   BINSYZ=', 2F7.1, ' ABUNJ', &
    12F6.2)

180 Format (10E12.3)
End Subroutine


Subroutine atread(num)
! CALLS ATRD NUM TIMES TO READ IN ATOMIC DATA
  Common /fel/wj(12, 220), fj(12, 220), e3j(12, 220), pwr(12, 220), &
    c1j(12, 30), c2j(12, 30), ej(12, 30), eaj(12, 30), s2j(12, 30), &
    llj(12, 30), s3j(12, 30), s4j(12, 30), s5j(12, 30)
  Common /params/nj(12), abunj(12), abund, binmin, binsyz, nbin
  Common /dat/e(30), ea(30), s2(30), wave(220), e3(220), f(220), ll(30), &
    s3(30), s4(30), s5(30)

  Do no = 1, num
    Call atrd(no)
    n = nj(no)
    Do j = 1, n
      ej(no, j) = e(j)
      eaj(no, j) = ea(j)
      s2j(no, j) = s2(j)
      llj(no, j) = ll(j)
      s3j(no, j) = s3(j)
      s4j(no, j) = s4(j)
      s5j(no, j) = s5(j)
    End Do
    Do l = 1, 220
      e3j(no, l) = e3(l)
      fj(no, l) = f(l)
      wj(no, l) = wave(l)
    End Do
  End Do
  Return

End Subroutine
Subroutine feline(t, dene, num, ipr, jpr, jcont, iphot, idens, icx)
  Common /bln/bln(1000), blnmin, blnsyz, nbln
  Common /fel/wj(12, 220), fj(12, 220), e3j(12, 220), pwr(12, 220), &
    c1j(12, 30), c2j(12, 30), ej(12, 30), eaj(12, 30), s2j(12, 30), &
    llj(12, 30), s3j(12, 30), s4j(12, 30), s5j(12, 30)
  Common /dat/e(30), ea(30), s2(30), wave(220), e3(220), f(220), ll(30), &
    s3(30), s4(30), s5(30)
  Common /params/nj(12), abunj(12), abund, binmin, binsyz, nbin
  Common /result/conce(30), ca(30), power(220), rhy, heneut, heplus, dne, &
    pcool, pou, pot, re, tu, pm(4)
  Common /aug/caj(12, 30)
  Common /rates/c1(30), c2(30), cth
  Common /com/cnc(12, 30), ptot(12, 220), abin(1000), bin(1000)
  Common /contin/brmev(1000), recev(1000), tufev(1000)
  Dimension :: pc(12)

  Do i = 1, nbin
    recev(i) = 0.
    tufev(i) = 0.
    brmev(i) = 0.
    bin(i) = 0.
  End Do

  nlines = 220
  st = sqrt(t)
  pcool = 0.
  abund = 12.
  conce(2) = rhy/(1.+rhy)
  If (nbin/=0) Call brems(t, 1)
  If (nbin/=0) Call recems(t, 1, 0, 0)

! ELEMENT LOOP
  Do no = 1, num
    re = 0.
    tu = 0.
    abund = abunj(no)
    Do l = 1, nlines
      power(l) = 0.
      wave(l) = wj(no, l)
      f(l) = fj(no, l)
      e3(l) = e3j(no, l)
    End Do

    n = nj(no)
    nn = n + 1
    Do j = 1, n
      e(j) = ej(no, j)
      ea(j) = eaj(no, j)
      s2(j) = s2j(no, j)
      conce(j) = cnc(no, j)
      ll(j) = llj(no, j)
      s3(j) = s3j(no, j)
      s4(j) = s4j(no, j)
      s5(j) = s5j(no, j)
    End Do

    e(n+1) = 0.
    ll(n+1) = 0.
    conce(n+1) = cnc(no, n+1)

    Call single(t, dene, no, io, ipr, jpr, jcont, iphot, idens, icx)

    If (no==1) heneut = conce(1)
    If (no==1) heplus = conce(2)
    Do j = 1, nn
      If (jcont==0) cnc(no, j) = conce(j)
      caj(no, j) = ca(j)
      c1j(no, j) = c1(j)
      c2j(no, j) = c2(j)
    End Do

    If (n==18 .Or. n==20 .Or. n==28) Go To 100
    If (nbin/=0) Call brems(t, n)
    If (nbin/=0) Call recems(t, n, 0, 0)
100 Continue

    Do l = 1, 220
      pwr(no, l) = power(l)
    End Do

    If (nbin/=0) Call bund(no)
    pc(no) = pcool + tu
  End Do
  pcool = 0.
  Do no = 1, num
    pcool = pcool + pc(no)
  End Do

  If (nbin==0) Go To 110
  Do i = 1, nbin
    bin(i) = bin(i) + recev(i) + brmev(i) + tufev(i)
  End Do

  If (jpr==0) Go To 110
  Write (7, 120) t, binmin, binsyz
  Do i = 1, nbin
    hnu = binmin + i*binsyz - binsyz
    Write (7, 130) i, hnu, recev(i), brmev(i), tufev(i), bin(i)
  End Do
110 Continue
  If (nbln/=0) Call blnd(num)
  Return
120 Format ('1', ' TEMP', F10.0, '    BINMIN', F7.2, '    BINSYZ', F7.2, /, &
    ('  BIN    HNU',9X,'RECEMS     BREMS    2 PHOT     TOTAL'), /, /)
130 Format ((1X,I5,F10.1,3X,4E10.3,1X))

End Subroutine
Subroutine atrd(no)
  Common /dat/e(30), ea(30), s2(30), wave(220), e3(220), f(220), ll(30), &
    s3(30), s4(30), s5(30)
  Common /params/nj(12), abunj(12), abund, binmin, binsyz, nbin
  Integer :: aa
! READ ATOMIC DATA
  ix = 0
100 Read (2, 120) n, j, e(j), ea(j), s2(j), s3(j), s4(j), s5(j), ll(j)
  iy = ll(j)
  If (iy<=0) Go To 110
  Do l = 1, iy
    aa = ix + 3*l - 3
    Read (2, 130)(wave(aa+k), e3(aa+k), f(aa+k), k=1, 3)
  End Do
  ix = ix + 3*iy
110 If (n>j) Go To 100
  nj(no) = n
  Return
120 Format (2I5, F5.0, 5F7.0, I5)
130 Format (9F6.0)

End Subroutine
Subroutine single(t, dene, no, io, iprint, jprint, jcont, iphot, idens, icx)
  Integer :: bb
  Common /dat/e(30), s(30), c(30), wave(220), e3(220), f(220), ll(30), &
    sig(30), alf(30), es(30)
  Common /params/nj(12), abunj(12), abund, binmin, binsyz, nbin
  Common /rates/c1(30), c2(30), cth
  Common /result/conce(30), gndrec(30), power(220), rhy, heneut, heplus, dne, &
    pcool, pou, pot, re, tu, pm(4)
  Common /strmet/popm(4, 12)

  Dimension :: jdum(5), ldum(5), wdum(5), pdum(5)

  pcool = 0.
  abund = 12.
  st = sqrt(t)
! IONIZATION STATE OF H MUST BE PASSED FROM OUTSIDE
  re = 0.
  tu = 0.
  abund = abunj(no)
  n = nj(no)
  ix = 0
  nn = n + 1
  e(n+1) = 0.
  ll(n+1) = 0
  Call nquil(t, n, dene, jcont, iphot, idens, icx)
  popm(1, no) = pm(1)
  popm(2, no) = pm(2)
  popm(3, no) = pm(3)
  popm(4, no) = pm(4)
  loc = 0
  If (iprint==0) Go To 100
  Write (7, 170) n, abund, t, dene
  If (n>=6) Write (7, 180) n, pm
  Write (7, 190)(k, conce(k), k=1, nn)
  If (iprint>=3) Write (7, 190)(k, c1(k), k=1, nn)
  If (iprint>=3) Write (7, 190)(k, c2(k), k=1, nn)
100 iik = ix
  ix = 0
  Do il = 1, n
    iy = 3*ll(il)
    If (iy<=0.) Go To 130
    cn = conce(il)
    pw = 0.
    Do l = 1, iy
      bb = ix + l
      If (cn<.00001) Go To 110
      If (wave(bb)<=0.) Go To 120
      g = gaunt_rs(t, e3(bb), n, il, l, dene)
      pw = 10.**(abund+11.)*cn*(alphadi(n,il,l,bb,t)*e3(bb)*1.602E-12+2.72E-15 &
        *f(bb)*exp(-11590*e3(bb)/t)*g/st)*12398./(wave(bb)*e3(bb))
110   power(bb) = pw
      pcool = pcool + pw
120 End Do
    If (il==n-1) Call heseq(n, il, t, dene, ix)
    If (il==n) Call hyseq(n, il, t, dene, ix)
130 ix = ix + iy
  End Do

  If (iprint<=1) Go To 160
  ix = 0
  Do il = 1, n
    iy = 3*ll(il)
    If (conce(il)<=.001) Go To 150
    Do l = 1, iy
      bb = ix + l
      If (wave(bb)<=0.) Go To 140
      loc = loc + 1
      jdum(loc) = il
      ldum(loc) = l
      wdum(loc) = wave(bb)
      pdum(loc) = power(bb)
      If (loc==5) Write (7, 200) n, (jdum(k), ldum(k), wdum(k), pdum(k), k=1, &
        5)
      loc = mod(loc, 5)
140 End Do
150 ix = ix + iy
  End Do

  If (loc/=0) Write (7, 200) n, (jdum(k), ldum(k), wdum(k), pdum(k), k=1, loc)
160 Continue
  If (iprint/=0) Write (7, 210) pcool, re, tu
  Return
170 Format ('0LINES FOR N=', I2, '     ABUND=', F5.2, '     T=', E9.3, &
    '     DENE=', E9.3)
180 Format (' CONCE FOR N=', I3, ' METASTABLES ', 4E10.3)
190 Format (10(I3,E10.2), /, 10(I3,E10.2), /, 10(I3,E10.2), /)
200 Format (I4, 5(I4,I3,F8.2,E10.3))
210 Format (' PCOOL=', E10.3, ' RE=', E10.3, ' TU=', E10.3)
End Subroutine


Function cion(n, j, e, t)
! SM YOUNGER JQSRT 26, 329; 27, 541; 29, 61   WITH MOORES FOR UNDONE
! A0 FOR B-LIKE ION HAS TWICE 2S PLUS ONE 2P  AS IN SUMMERS ET AL
! CHI = kT / I

  Dimension :: a0(30), a1(30), a2(30), a3(30), b0(30), b1(30), b2(30), b3(30), &
    c0(30), c1(30), c2(30), c3(30), d0(30), d1(30), d2(30), d3(30)
  Data a0/13.5, 27.0, 9.07, 11.80, 20.2, 28.6, 37.0, 45.4, 53.8, 62.2, 11.7, &
    38.8, 37.27, 46.7, 57.4, 67.0, 77.8, 90.1, 106., 120.8, 135.6, 150.4, &
    165.2, 180.0, 194.8, 209.6, 224.4, 239.2, 154.0, 268.8/
  Data a1/ -14.2, -60.1, 4.30, 27*0./
  Data a2/40.6, 140., 7.69, 27*0./

  Data a3/ -17.1, -89.8, -7.53, 27*0./
  Data b0/ -4.81, -9.62, -2.47, -3.28, -5.96, -8.64, -11.32, -14.00, -16.68, &
    -19.36, -4.29, -16.7, -14.58, -16.95, -19.93, -23.05, -26.00, -29.45, &
    -34.25, -38.92, -43.59, -48.26, -52.93, -57.60, -62.27, -66.94, -71.62, &
    -76.29, -80.96, -85.63/
  Data b1/9.77, 33.1, -3.78, 27*0./
  Data b2/ -28.3, -82.5, -3.59, 27*0./

  Data b3/11.4, 54.6, 3.34, 27*0./
  Data c0/1.85, 3.69, 1.34, 1.64, 2.31, 2.984, 3.656, 4.328, 5.00, 5.672, &
    1.061, 1.87, 3.26, 5.07, 6.67, 8.10, 9.92, 11.79, 7.953, 8.408, 8.863, &
    9.318, 9.773, 10.228, 10.683, 11.138, 11.593, 12.048, 12.505, 12.96/
  Data c1/0., 4.32, .343, 27*0./
  Data c2/0., -2.527, -2.46, 27*0./

  Data c3/0., .262, 1.38, 27*0./
  Data d0/ -10.9, -21.7, -5.37, -7.58, -12.66, -17.74, -22.82, -27.9, -32.98, &
    -38.06, -7.34, -28.8, -24.87, -30.5, -37.9, -45.3, -53.8, -64.6, -54.54, &
    -61.70, -68.86, -76.02, -83.18, -90.34, -97.50, -104.66, -111.82, -118.98, &
    -126.14, -133.32/
  Data d1/8.90, 42.5, -12.4, 27*0./
  Data d2/ -35.7, -131., -8.09, 27*0./

  Data d3/16.5, 87.4, 1.23, 27*0./

  cion = 0.
  chir = t/(11590.*e)
  If (chir<=.0115) Return
  chi = amax1(chir, 0.1)
  ch2 = chi*chi
  ch3 = ch2*chi
  alpha = (.001193+.9764*chi+.6604*ch2+.02590*ch3)/(1.0+1.488*chi+.2972*ch2+ &
    .004925*ch3)
  beta = (-.0005725+.01345*chi+.8691*ch2+.03404*ch3)/ &
    (1.0+2.197*chi+.2457*ch2+.002503*ch3)
  j2 = j*j
  j3 = j2*j

  iso = n - j + 1
  a = a0(iso) + a1(iso)/j + a2(iso)/j2 + a3(iso)/j3
  b = b0(iso) + b1(iso)/j + b2(iso)/j2 + b3(iso)/j3
  c = c0(iso) + c1(iso)/j + c2(iso)/j2 + c3(iso)/j3

  d = d0(iso) + d1(iso)/j + d2(iso)/j2 + d3(iso)/j3
! FE II EXPERIMENTAL IONIZATION MONTAGUE ET AL: D. NEUFELD FIT
  If (n/=26 .Or. j/=2) Go To 100
  a = -13.825
  b = -11.0395
  c = 21.07262
  d = 0.

100 Continue
  ch = 1./chi
  fchi = 0.3*ch*(a+b*(1.+ch)+(c-(a+b*(2.+ch))*ch)*alpha+d*beta*ch)
  If (iso>=4 .And. iso<=10) fchi = fchi*1.59
! correct Younger JQSRT 27, 541 from table to graphs
  cion = 2.2E-6*sqrt(chir)*fchi*exp(-1./chir)/(e*sqrt(e))
  Return
End Function
Function alphadi(n, j, l, ln, t)
! DIELECTRONIC RECOMBINATION : BURGESS AND TWORKOWSKI FOR H-LIKE

  Common /dat/e(30), s(30), c(30), wave(220), e3(220), f(220), ll(30), &
    sig(30), alf(30), es(30)
! BURGESS AND TWORKOWSKI
  z12 = j - 12.
  twork = .84 + .5/j**2 + .03*z12/(1.+4.5E-5*z12**3)
  If (n/=j) twork = 1.
  alphadi = 0.
  dl = delt(n, j, l)
  If (dl<=0.) Return
  zf = j - 1.
  b = sqrt(zf)*(zf+1.)**2.5/sqrt(zf*zf+13.4)
  x = e3(ln)/((zf+1.)*13.6)
  a = sqrt(x)/(1.+.105*x+.015*x*x)
  ebar = e3(ln)/(1.+.015*zf**3/(zf+1.)**2)
  If (n-j/=1) Go To 100
  b = (0.5*zf/sqrt(zf))*b
  x = .75*j
  a = sqrt(x)*0.5/(1.+0.21*x+0.03*x*x)
! younger jqsrt 29, 67 for he-like ions

100 Continue
  alphadi = .0030*t**(-1.5)*a*b*f(ln)*twork*dl*exp(-ebar*11590./t)
  Return

End Function
Function autoin(n, j, t)
! INNERSHELL EXCITATION FOLLOWED BY AUTOIONIZATION FOR NA,MG,AL,SI
! SEQUENCES.  CRANDALL ET AL PRA 25,143 FOR MG,SI.  MANN FOR FE* .75.
! COMPROMISE BETWEEN COWAN AND MANN AND DOUBLE THAT.
! USE 2S AT 2P ENERGY, ETC FOR B,C,... SEQUENCES
  Dimension :: et(20), om(20)
  Data et/0., 55., 0., 115., 0., 190., 0., 280., 0., 375., 5*0., 802., 0., &
    1002., 0., 0./
  Data om/0., .34, 0., .16, 0., .18, 0., .20, 0., .22, 5*0., .29, 0., .3, 0., &
    0./

  autoin = 0.
  If (n<=10) Return
  i = n - j + 1
  If (i<=10 .Or. i>=15) Return
  autoin = 8.63E-6*om(n-10)*exp(-11590.*et(n-10)/t)/sqrt(t)
! ASSUMES THAT NUMBER OF 32,3P ELECTRONS DOESN'T MATTER
  Return

End Function

Function alflo(n, j, t)
! Low Temperature Dielectronic Recombination from Nussbaumer and
! Storey: C - Si for 1000 - 60,000 K

  Dimension :: iion(8, 16), a(25), b(25), c(25), d(25), f(25)

  Data iion/40*0, 0, 1, 2, 3, 0, 0, 0, 0, 0, 4, 5, 6, 7, 0, 0, 0, 0, 8, 9, 10, &
    11, 12, 0, 0, 8*0, 0, 13, 14, 15, 16, 17, 18, 19, 8*0, 0, 20, 6*0, 0, 21, &
    22, 5*0, 0, 23, 24, 25, 4*0, 8*0, 8*0/

  Data a/.0108, 1.8267, 2.3196, 0.0, 0.0320, -.8806, .4134, 0.0, -.0036, 0.0, &
    .0061, -2.8425, 0.0, .0129, 3.6781, -.0254, -.0141, 19.928, 5.4751, &
    1.2044, .0219, .7086, -.0219, 3.2163, .1203/

  Data b/ -.1075, 4.1012, 10.7328, .6310, -.6624, 11.2406, -4.6319, .0238, &
    .7519, 21.879, .2269, .2283, 0.0, -.1779, 14.1481, 5.5365, 33.8479, &
    235.0536, 203.9751, -4.6836, -.4528, -3.1083, .4364, -12.0571, -2.690/

  Data c/.2810, 4.8443, 6.8830, .1990, 4.3191, 30.7066, 25.9172, .0659, &
    1.5252, 16.2730, 32.1419, 40.4072, 0., .9353, 17.1175, 17.0727, 43.1608, &
    152.5096, 86.9016, 7.6620, 2.5427, 7.0422, .0684, 16.2118, 19.1943/

  Data d/ -.0193, .2261, -.1824, -.0197, .0003, -1.1721, -2.2290, .0349, &
    -.0838, -.7020, 1.9939, -3.4956, 0., -.0682, -.5017, -.7225, -1.6072, &
    9.1413, -7.4568, -.5930, -.1678, .5998, -.0032, -.5886, -.1479/

  Data f/ -.1127, .5960, .4101, .4398, .5946, .6127, .2360, .5334, .2769, &
    1.1899, -.0646, 1.7558, 0., .4516, .2313, .1702, .1942, .1282, 2.5145, &
    1.6260, .2276, .4194, .1342, .5613, .1118/

  alflo = 0.
  If (j==1 .Or. j>8) Return

  If (n>14) Return

  t4 = t*.0001

  If (t4<0.1 .Or. t4>6.) Return
  If (n==6 .And. j==2 .And. t4<0.2) Return
  If (n==8 .And. j==4 .And. t4<0.2) Return
  If (n==8 .And. j==6 .And. t4<0.4) Return
  If (n==10 .And. j==8 .And. t4<0.2) Return
  If (n==12 .And. j==2 .And. t4<0.15) Return

  If (n==14 .And. j==3 .And. t4<0.15) Return
  ij = iion(j, n)

  If (ij==0) Return

  alflo = 1.0E-12*(a(ij)/t4+b(ij)+c(ij)*t4+d(ij)*t4*t4)*exp(-f(ij)/t4)/t4**1.5
  Return



End Function
!FUNCTION ALFLO(N,J,T)
!LOW TEMPERATURE DIELECTRONIC RECOMBINATION FROM NUSSBAUMER AND STOREY
!C THROUGH O AND T < 60000 K
!DIMENSION IION(8,3),A(12),B(12),C(12),D(12),F(12)
!dimension asi(4),bsi(4),csi(4),dsi(4),fsi(4)
!DATA IION/0,1,2,3,0,0,0,0, 0,4,5,6,7,0,0,0, 0,8,9,10,11,12,0,0/
!DATA A/.0108,1.8267,2.3196,0.0,.032,-.8806,.4134,0.,-.0036,0.,
!1  .0061,-2.8425/
!DATA B/-.1075,4.1012,10.7328,.6310,-.6624,11.2406,-4.6319,.0238,
!1  .7519,21.8790,.2269,.2283/
!DATA C/.2810,4.8443,6.883,.199,4.3191,30.7066,25.9172,.0659,
!1  1.5252,16.273,32.1419,40.407/
!DATA D/-.0193,.2261,-.1824,-.0197,.0003,-1.1721,-2.229,.0349,
!1  -.0838,-.7020,1.9939,-3.4956/
!DATA F/-.1127,.5960,.4101,.4398,.5946,.6127,.2360,.5334,.2769,
!1  1.1899,-.0646,1.7558/
!DATA ASI/0.,-.0219,3.2163,0.1203/
!DATA BSI/0.,.4364,-12.0571,-2.690/
!DATA CSI/0.,.0684,16.2118,19.1943/
!DATA DSI/0.,-.0032,-.5886,-.0099/
!DATA FSI/0.,.1342,.5613,.1118/
!ALFLO = 0.
!IF (J .EQ. 1 .OR. N-J .LE. 1) RETURN
!T4 = T *.0001
!IF (T4 .GT. 6.0) RETURN
!IF(T4 .LT. 0.1) RETURN
!IF (N .EQ. 6 .AND. J .EQ. 2 .AND. T4 .LT. 0.2) RETURN
!IF (N .EQ. 8 .AND. J .EQ. 4 .AND. T4 .LT. 0.2) RETURN
!IF (N .EQ. 8 .AND. J .EQ. 6 .AND. T4 .LT. 0.4) RETURN
!IF (N .EQ. 12 .OR. N .EQ. 14) GO TO 900
!IF (N .LT. 6 .OR. N .GT. 8) RETURN
!IJ = IION(J,N-5)
!ALFLO = 1.0E-12 * (A(IJ)/T4 + B(IJ) + C(IJ)*T4 + D(IJ)*T4**2) *
!1  EXP(-F(IJ) / T4) / T4 ** 1.5
!900   IF (N .EQ. 14) GO TO 950
!IF (J .NE. 2) RETURN
!ALFLO = 1.0E-12 * (1.2044/T4 -4.6836 +7.6620*T4 -.5930*T4**2) *
!1  EXP(-1.6260 / T4) / T4 ** 1.5
!RETURN
!950   IF (J .GT. 4) RETURN
!ALFLO = 1.0E-12*(ASI(J)/T4 +BSI(J) +CSI(J)*T4 +DSI(J)*T4**2) *
!1  EXP(-FSI(J) / T4) / T4 ** 1.5
!RETURN
!END
Function dimet(n, j, t, b, dd)
  Dimension :: noj(30), e(12, 3), f(12, 3)
  Data noj/0, 1, 3*0, 2, 3, 4, 0, 5, 0, 6, 0, 7, 0, 8, 0, 9, 0, 10, 5*0, 11, &
    0, 12, 2*0/
  Data e/0, 10.55, 13.4, 16.3, 22.1, 27.9, 40.0, 40.5, 46.7, 52.8, 76.3, 86.2, &
    0., 12.28, 16.06, 19.8, 27.4, 35.0, 42.7, 51.0, 60.9, 73.5, 122., 145., &
    5*0., 4.46, 9.55, 14.5, 19.5, 24.5, 40.7, 51./
  Data f/0., .26, .21, .18, .16, .12, .10, .091, .081, .075, .003, .002, 0., &
    .16, .17, .14, .084, .076, .071, .067, .063, .059, .049, .046, 5*0., .61, &
    .564, 0.458, .408, .366, .30, .25/

  ij = n - j - 2
  If (ij==9) ij = 3
  z = j - 1
  no = noj(n)
  x = e(no, ij)/(j*13.6)
  a = sqrt(x)/(1.+.105*x+.015*x*x)
  ebar = e(no, ij)/(1.+.015*z*z*z/(j*j))
  dimet = .003*t**(-1.5)*a*b*f(no, ij)*exp(-11590.*ebar/t)*dd
  Return

End Function
Subroutine bund(no)
! PUTS EMISSION LINES INTO ARRAY BIN.
  Common /params/nj(12), abunj(12), abund, binmin, binsyz, nbin
  Common /com/cnc(12, 30), ptot(12, 220), abin(1000), bin(1000)
  Common /dat/v(30), s(30), c(30), wave(220), e3(220), f(220), ll(30), &
    sig(30), alf(30), es(30)
  Common /result/conce(30), gndrec(30), power(220), rhy, heneut, heplus, dne, &
    pcool, ppu, pot, re, tu, pm(4)

  ibn(e) = (e-binmin)/binsyz + 1
  Do l = 1, 220
    If (wave(l)<=0) Go To 100
    e = 12399./wave(l)
    ik = ibn(e)
    If (ik<1) Go To 100
    If (ik>nbin) Go To 100

!   MOVE PHOTONS IN SAME BIN AS POPULAR K EDGE BUT HIGHER ENERGY UP
!   TO NEXT HIGHER BIN:  B,C,N,O,F,NE,AL EDGES

    kedge = ibn(187.03)
    If (ik==kedge .And. e>187.03) ik = ik + 1
    kedge = ibn(284.05)
    If (ik==kedge .And. e>284.05) ik = ik + 1
    kedge = ibn(400.07)
    If (ik==kedge .And. e>400.06) ik = ik + 1
    kedge = ibn(532.09)
    If (ik==kedge .And. e>532.09) ik = ik + 1
    kedge = ibn(692.14)
    If (ik==kedge .And. e>692.14) ik = ik + 1
    kedge = ibn(874.16)
    If (ik==kedge .And. e>874.16) ik = ik + 1
    kedge = ibn(1559.3)
    If (ik==kedge .And. e>1559.3) ik = ik + 1
    bin(ik) = bin(ik) + power(l)
100 End Do
  Return

End Subroutine


Function fecool(dene, t, ipr)
! COOLING DUE TO FE II AND FE III ASSUMING THAT CONCE=CONCE(MG)
! IF IRON IS NOT EXPLICITLY CALCULATED.
! NUSSBAUMER AND STOREY FE II OMEGA'S AND ROBB FE III
  Common /params/nj(12), abunj(12), abund, binmin, binsyz, nbin
  Common /com/cnc(12, 30), ptot(12, 220), abin(1000), bin(1000)
  Common /fe/feinf, feperm, feforb, feiii, feiiio

  Common /feir/fe16, fe86, fe25

  st = sqrt(t)
  c2 = cnc(6, 2)
  c3 = cnc(6, 3)
  If (nj(11)==26) c2 = cnc(11, 2)
  If (nj(11)==26) c3 = cnc(11, 3)
  fct = 1.6E11*10.**(abunj(11)-12.)*8.63E-6*c2/st
  feinf = fct*(.31*.048*exp(-553./t)+.42*.99*exp(-11400./t))
  fe16 = fct*.39*.99*exp(-11400./t)*1130.*st/(dene+1130.*st)
  fe86 = fct*.075*1.67*exp(-19400./t)*7050.*st/(dene+7050.*st)
  fe25 = fct*.3*.048*exp(-553./t)*580.*st/(dene+580.*st)
  feinf = feinf*580.*st/(dene+580.*st)
  feperm = fct*3.4*4.47*exp(-55000./t)
  feforb = fct*0.1*2.9*exp(-33000./t)*1.2E6*st/(dene+1.2E6*st)
  fct = 1.6E11*10.**(abunj(11)-12.)*8.63E-6*c3/st
  feiii = fct*.30*.054*exp(-625./t)*1590.*st/(dene+1590.*st)

  feiiio = fct*.43*2.66*exp(-30840./t)*521000.*st/(dene+521000.*st)
  fc = feinf + feperm + feforb + feiii + feiiio
  fecool = fc
  If (ipr/=0) Print 100, fc, feinf, feperm, feforb, feiii, feiiio
  Return
100 Format (' IRON COOLING; FC,II IR, II UV, II OPT, III IR, III   OPT ', &
    6F9.6)

End Function
Function cthr(n, j, t)
! SCOTT'S RATES FOR H I + J TO H II + (J-1)
! THROUGH ARGON : NEUFELT'S FE II - FE III

  Dimension :: a(30), b(30)
  Common /result/conce(30), gndrec(30), power(220), rhy, heneut, heplus, dne, &
    pcool, pou, pot, re, tu, pm(4)
  Data a/0., 0., .00001, 0., 5.9, 6.2, 7.8, 9.4, 11., 13., 14., 16.1, 17.6, &
    19., 20.4, 21.8, 23.3, 13*0./
  Data b/4*0., .25, .18, .13, .1, .08, .05, .04, 19*0./

  cthr = 0.
  If (j==1) Return
  t4 = amax1(t*.0001, 0.1)
  t4 = amin1(t4, 10.)
! NEUFELD'S FE II - FE III, FE III-FE IV, AND NI II - NI III RATES
  If (n==26 .And. j==3) cthr = 1.0E-9*(1.25+.25*alog10(t4))
  If (n==26 .And. j==4) cthr = 3.4E-9*sqrt(t4)
  If (n==28 .And. j==3) cthr = 1.0E-9*(.34+1.86*t4)
  If (n>18) Return
  cthr = a(j)*(1.+b(j)*t4)*1.0E-9
  If (j>=6) Return
  If (n/=2) Go To 100
  cthr = 1.9E-15
  If (j==3) cthr = 1.7E-13
  Return
100 If (n/=6) Go To 110
! ALBERT'S C III TRIPLET P
  If (j==3) cthr = 2.5E-9*pm(1)*exp(-15000./t)
  If (j==4) cthr = 2.9E-9
  If (j==5) cthr = 7.6E-10*t4**1.48
  Return
110 If (n/=7) Go To 120
  cthr = 1.1E-12/(1.+.1*t4)
  If (j==3) cthr = 5.2E-10
  If (j==4) cthr = 2.7E-9*t4**.926
  If (j==5) cthr = 1.7E-10*t4**1.40
! ALBERT'S N VI
  If (j==6) cthr = 6.6E-10
  Return
120 If (n/=8) Go To 130
  cthr = 4.0E-10
  If (j==3) cthr = 7.7E-10*sqrt(t4)
  If (j==4) cthr = 2.1E-9
  If (j==5) cthr = 1.4E-10*(1.+t4)
! ALBERT'S O VII : NEED O VI
  If (j==7) cthr = 5.4E-8
  Return
130 If (n/=10) Go To 140
  If (j==4) cthr = 3.8E-9*sqrt(t4)
  Return
140 If (n/=12) Go To 150
  If (j==4) cthr = 4.4E-9*(1.+.37*t4)
  Return
150 If (n/=14) Go To 160
  If (j==3) cthr = 4.0E-9*t4**.23
  If (j==4) cthr = 4.1E-10
  If (j==5) cthr = 2.2E-9*(1.+.1*t4)
  Return
160 If (n/=16) Go To 170
  If (j==4) cthr = 2.5E-9
  If (j==5) cthr = 7.0E-9
  Return
170 If (n/=18) Go To 180
  If (j==4) cthr = 4.4E-8*t4**.27
  Return
180 Return

End Function
Function cthi(n, j, t)
! SCOTT'S RATES FOR   H II  +  J    TO    HI + J+1
! THROUGH ARGON
  t4 = amax1(.0001*t, 0.1)
  t4 = amin1(t4, 10.)
  cthi = 0.
  If (j>2) Return
  If (j==2) Go To 100
  If (n==6) cthi = 2.5E-15
  If (n==7) cthi = 3.3E-13*t4**.12
  If (n==8) cthi = 3.3E-10
  If (n==16) cthi = 6.7E-10
  Return
100 If (n==12) cthi = 2.6E-11*exp(-6.0/t4)
  If (n==14) cthi = 9.6E-10*exp(-2.74/t4)
! NEUFELD'S FE II - FE III RATES
  If (n==26) cthi = 1.0E-9*1.666*(1.25+.25*alog10(t4))*exp(-3.005/t4)
  Return

End Function
Function cther(n, j, t)
! FITS TO SCOTT'S RATES FOR  HE I + J   TO   HE II + J-1
! 1000  TO 100,000 K THROUGH ARGON
  Dimension :: a(2, 9)
  Data a/4*0., .059, .00001, 1.1, 0.7, .00001, 1.7, .075, 2.2, .096, 1.2, 1.1, &
    .0008, .00001, 1.0/

  cther = 0.
  If (j==1) Return
  t4 = amax1(.0001*t, 0.1)
  t4 = amin1(t4, 10.)
  cr = 5.4E-10*t4
  If (j==2 .Or. n>18) Return
  cther = amax1(cr, (-.1+.4*j)*1.0E-9)
  If (j>=6) Return
  If (j/=3) Go To 100
  cther = 0.
! INCLUDE ALBERT'S RATE TO 1D OF N+
  If (n==7) cther = 3.0E-10*(1.+.26*t4) + 6.2E-11
  If (n==8) cther = 3.3E-10*t4**.7
  If (n==18) cther = 1.3E-10
  Return
100 If (n/=7) Go To 110
  cther = 1.1E-10
  If (j==5) cther = 2.0E-9
  Return
110 Continue
  in = n/2
  cther = a(j-3, in)*1.0E-9
  If (j==5) Go To 120
  If (n==6) cther = 5.1E-11*t4**1.46
  If (n==14) cther = 9.6E-10*t4**.55
  If (n==16) cther = 1.1E-9*t4**.63
  Return
120 If (n==18) cther = 1.0E-9*t4**.91
  Return

End Function
Function cthei(n, j, t, dene)
! SCOTT'S RATES FOR HE+  + J  TO HE0  + J+1
! THROUGH ARGON  0.1 < T4 < 10.
  cthei = 0.
  t4 = 0.0001*t
  If (j==1 .Or. j>=4) Return
  If (j==3) Go To 110
  If (n==6) cthei = 5.3E-10*exp(-13.2/t4)
! ALBERT'S RATES FOR N+ FROM 1D LEVEL
  If (n/=7) Go To 100
  qdw = 5.14E-06*dene/sqrt(t)
  qup = 2.86E-6*exp(-21826./t)
  r21 = qup/(qdw+.0041)
  d1 = r21/(1.+r21)
  p3 = 1. - d1
  cthei = p3*4.1E-10*exp(-10.4/t4) + d1*1.1E-10*exp(-3.59/t4)
100 Continue
  If (n==14) cthei = 2.9E-10*exp(-8.7/t4)
  If (n==18) cthei = 1.1E-10*exp(-3.57/t4)
  Return
110 Continue
  If (n==14) cthei = 3.4E-9*exp(-10.5/t4)
  If (n==16) cthei = 4.8E-10*exp(-15.7/t4)
  If (n==26) cthei = 1.2E-9*sqrt(t4)*exp(-12./t4)
  Return

End Function

Function popmet(n, ij, emet, dene, t)
! RETURNS POPULATION OF METASTABLE LEVEL OF BE,B,C, OR MG-LIKE ION;  CA
! RESONANCE CONTRIBUTIONS CALC FOR HIGH Z WITH SMITH ET AL METHOD
! PASSES REDUC BACK TO GAUNT_RS FOR INTERCOMBINATION LINES

  Real :: mgmm, mgm1, mgm5, mg21, mg5, mg1
  Dimension :: rm5(12), ebe(12), pops(4)
  Dimension :: mgmm(3, 7), mgm1(7), mg5(7), mg21(2, 7), emg(7)
  Dimension :: nj(30), a21(2, 12), b21(3, 12), c21(3, 12), rm1(12)
  Dimension :: rmm(3, 12), bm1(12), bmm(3, 12), cm1(3, 12), cmm(2, 12)
  Dimension :: ecm(2, 12), r(3, 3), x(3), c(3), cpr(3)
  Dimension :: rb5(10), rc5(10), rres(4, 10), eres(4, 10)

  Dimension :: rbe(10), bee(10)
  Common /result/conce(30), gndrec(30), power(220), rhy, heneut, heplus, dne, &
    pcool, pou, pot, re, tu, pm(4)
  Common /interc/reduc(4)

  Common /metln/pln(4, 3, 8), tln(4, 3, 8)
  Data rm5/0., 18.6, 11.7, 6.09, 2.79, 1.54, .93, .61, .41, .31, .13, .10/
  Data ebe/0., 6.20, 7.87, 9.51, 12.81, 16.1, 19.5, 23.11, 26.9, 28.9, 55.5, &
    63.1/
  Data mg21/430., 5.8E-4, 22000., .008, 1.63E+5, .000, 9.E+5, .12, 3.3E+6, &
    .26, 4.5E+7, 1.00, 8.E+7, 2./


  Data mg5/0.2, 0.5, .26, .035, .026, .0058, .0039/
! GILES QUOTED BY MENDOZA FOR S V, BALUJA FOR SI III
  Data mgm1/.045, .35, .072, .026, .016, .0069, .0047/
  Data mgmm/.5, .38, 1.5, 1.4, 1.06, 4.1, .25, .39, 1.32, .29, .21, .82, .17, &
    .13, .51, .033, .25, .43, .02, .13, .30/
  Data emg/1.63, 3.72, 5.50, 7.13, 8.72, 13.9, 17./
  Data nj/0, 1, 0, 0, 0, 2, 3, 4, 0, 5, 0, 6, 0, 7, 0, 8, 0, 9, 0, 10, 5*0, &
    11, 0, 12, 0, 0/
  Data a21/2*0., 77.2, .00282, 471., .00629, 1880., .0118, 16200., .0306, &
    84500., .0629, 320000., .112, 972600., .1823, 2.546E+6, .2781, 5.914E+6, &
    .4029, 7.3E+7, 1.21, 1.4E+8, 1.8/
  Data b21/3*0., 332., 130., 82., 1183., 130., 347., 4210., 499., 1470., &
    28400., 3170., 8910., 152000., 19800., 56200., 705000., 102000., 288000., &
    43790., 11613., 27300., 4.28E+06, 7.08E+05, 2.42E+06, 1.02E+07, 1.55E+06, &
    5.62E+06, 8.E+7, 1.E+7, 5.E+7, 2.E+8, 2.E+7, 1.E+8/


  Data c21/3*0., .00031, .0026, 32., .0041, .0342, 190., .022, .16, 385., &
    .599, 4.64, 6910., 4.67, 37.1, 45200., 26.1, 200., 197000., 115., 856., &
    680000., 424., 3000., 2.0E+06, 1370., 9090., 5.32E+06, 29400., 128000., &
    6.7E+7, 95000., 380000., 1.7E+8/
! REDUCED COLLISION RATES * 10. ** 7


  Data rm1/0., 9.384, 4.944, 2.814, 1.340, .7265, .4718, .3164, .2243, .1732, &
    .08, .062/
! INCLUDES PROTON RATES FROM MUNRO THROUGH SI, EXT F OR S, NONE ABOVE A
  Data rmm/3*0., 29.9, 13.4, 50.8, 17.7, 8.41, 32.2, 11.2, 6.00, 21.8, 6.48, &
    4.19, 13.9, 3.96, 3.04, 9.68, 2.91, 2.67, 7.96, 2.19, 2.42, 6.83, 1.25, &
    0.40, 1.79, .533, .317, 1.40, .046, .16, .67, .020, .12, .52/
  Data bm1/0., 14., 8.1, 5.3, 4.1, 3.01, 2.0, 1.5, 1.1, .95, .67, .58/
  Data bmm/3*0., 23., 7.9, 31., 17.1, 5.73, 20.8, 12.3, 4.17, 14.1, 6.52, &
    2.39, 7.91, 3.84, 1.47, 3.01, 2.65, 1.01, 3.16, 18.6, 13.9, 38., 1.68, &
    .57, 1.9, 1.34, .43, 1.45, .69, .18, .69, .55, .14, .54/
  Data cm1/3*0., 12.1, 12.9, 10., 51., 29.5, 23., 12.0, 8.6, 4.18, 9.1, 6.4, &
    2.92, 6.75, 4.72, 2.1, 5.39, 3.86, 1.55, 4.17, 3.04, 1.19, 2.84, 2.02, &
    1.00, 2.43, 1.57, 0.90, 1.1, .75, .60, 1.0, .56, .54/
  Data cmm/2*0., 12.9, .0005, 32.5, .001, 29., .0015, 17.5, .00205, 11.2, &
    .00173, 7.42, .00173, 5.44, .00518, 4.1, .0069, 3.52, .0138, 1.7, .03, &
    1.5, .05/
  Data ecm/2*0., 1.26, 2.68, 1.90, 4.05, 2.42, 5.34, 3.76, 7.93, 5.09, 10.6, &
    6.58, 13.4, 8.34, 16.5, 10.6, 20.1, 13.5, 24.5, 30.1, 45., 51., 62./
  Data rb5/0., 25., 21.6, 16.9, 8.80, 5.23, 3.64, 15.5, 2*0./
  Data rc5/10*0./
  Data bee/0., 10.5, 13.4, 16.3, 22.1, 28., 34., 40.5, 2*0./
  Data rbe/0., 115., 95.4, 78.4, 47.7, 32.1, 24.6, 18.2, 2*0./
  Data eres/4*0., 12.7, 9.25, 7.94, 0., 16.3, 12.5, 11.5, 0., 19.7, 15.8, &
    14.9, 0., 27.0, 22.3, 21.9, 0., 33.8, 28.8, 28.9, 4.35, 40.9, 35.8, 36.3, &
    10.3, 48.4, 11.6, 44.1, 15.8, 8*0./

  Data rres/4*0., 327., 95., 5.9, 0., 306., 111., 74.9, 0., 229., 97.8, 139., &
    0., 141., 60.7, 111., 0., 96., 41.3, 70.3, 68.0, 71.7, 30.3, 37.1, 956., &
    52.5, 22.0, 26.5, 64.3, 8*0./

  Do i = 1, 3
    Do j = 1, 3
      r(i, j) = 0.
    End Do

  End Do
! avoid overflow for low t
  reduc(ij) = 0.
  popmet = 0.
  If (emet*11590./t>=30.) Return
  st = sqrt(t)
  no = nj(n)
  v = dene*1.E-7/st
  Go To (100, 110, 120, 130) ij
100 Continue
  cn = conce(n-3)
  r(1, 1) = -(rm1(no)+rmm(1,no)*3.+rmm(2,no)*5.)*v
  r(2, 2) = -a21(1, no) - v*(rm1(no)+rmm(1,no)+rmm(3,no)*1.6666667)
  r(3, 3) = -a21(2, no) - v*(rm1(no)+rmm(2,no)+rmm(3,no))
  r5 = rm5(no)*v*exp(-11590.*ebe(no)/t)
  r(1, 1) = r(1, 1) - r5
  r(2, 2) = r(2, 2) - r5
  r(3, 3) = r(3, 3) - r5
  r(2, 1) = v*rmm(1, no)*3.
  r(3, 1) = v*rmm(2, no)*5.
  r(1, 2) = v*rmm(1, no)
  r(1, 3) = v*rmm(2, no)
  r(3, 2) = v*1.666667*rmm(3, no)
  r(2, 3) = v*rmm(3, no)
  c(1) = -rm1(no)*exp(-emet*11590./t)*v
  c(2) = c(1)*3.
  c(3) = c(1)*5.
  Go To 140
110 Continue
! B-LIKE IONS
  cn = conce(n-4)
  r5 = rb5(no)*v*exp(-11590.*(eres(2,no)-emet)/t)
  r(1, 1) = -(bm1(no)+bmm(1,no)*2.+bmm(2,no)*3.)*v - b21(1, no)
  r(2, 2) = -(bm1(no)+bmm(1,no)+bmm(3,no)*1.5)*v - b21(2, no)
  r(3, 3) = -(bm1(no)+bmm(2,no)+bmm(3,no))*v - b21(3, no)
  r(2, 1) = v*bmm(1, no)*2.
  r(3, 1) = v*bmm(2, no)*3.
  r(1, 2) = v*bmm(1, no)
  r(1, 3) = v*bmm(2, no)
  r(3, 2) = v*1.5*bmm(3, no)
  r(2, 3) = v*bmm(3, no)
  c(1) = -bm1(no)*exp(-emet*11590./t)*v/3.
  c(2) = c(1)*2.
  c(3) = c(1)*3.
  Go To 140
120 Continue
! C-LIKE IONS
  cn = conce(n-5)
  r5 = rc5(no)*v*exp(-11590.*(eres(3,no)-emet)/t)
  e1 = exp(-ecm(1,no)*11590./t)
  e2 = exp(-ecm(2,no)*11590./t)
  e3 = exp(-emet*11590./t)
  r(1, 1) = -(cm1(1,no)+cmm(1,no)*.2*e2/e1+cmm(2,no)*e3/e1)*v - c21(1, no)
  r(2, 2) = -(cm1(2,no)+cmm(1,no))*v - c21(2, no)
  r(3, 3) = -(cm1(3,no)+cmm(2,no))*v - c21(3, no)
  r(2, 1) = v*cmm(1, no)*.2*e2/e1
  r(1, 2) = v*cmm(1, no)
  r(3, 1) = v*cmm(2, no)*e3/e1
  r(1, 3) = v*cmm(2, no)
  r(3, 2) = 0.
  r(2, 3) = 0.
  c(1) = -cm1(1, no)*v*5.*e1/9.
  c(2) = -cm1(2, no)*v*e2/9.
  c(3) = -cm1(3, no)*v*e3*5./9.
  Go To 140
130 Continue
  cn = conce(n-11)
  v = 8.63E-6*dene/st
  np = no - 5
  mg1 = v*mgm1(np)
  mgm5 = v*mg5(np)*exp(-11590.*emg(np)/t)
  If (n==14) mgm5 = mgm5/(1.+8.6E-6*t)
  If (n==14) mg1 = mg1/(1.+t*3.5E-6)
  r(1, 1) = -mg1 - mgm5 - mgmm(1, np)*v - mgmm(2, np)*v
  r(2, 2) = -mg1 - mgm5 - mgmm(1, np)*v/3 - mgmm(3, np)*v/3 - mg21(1, np)
  r(3, 3) = -mg1 - mgm5 - mg21(2, np) - mgmm(2, np)*v/5. - mgmm(3, np)*v/5.
  r(2, 1) = v*mgmm(1, np)
  r(3, 1) = v*mgmm(2, np)
  r(1, 2) = v*mgmm(1, np)/3.
  r(1, 3) = v*mgmm(2, np)/5.
  r(3, 2) = v*mgmm(3, np)/3.
  r(2, 3) = v*mgmm(3, np)/5.
  c(1) = -mg1*exp(-11590.*emet/t)
  c(2) = c(1)*3.
  c(3) = c(1)*5.
140 Continue
  nn = noson(r, x, 3, 3)
  If (nn==0) Print 190
  Do i = 1, 3
    cpr(i) = 0.
    Do j = 1, 3
      cpr(i) = cpr(i) + r(i, j)*c(j)
    End Do
  End Do
  sum = 1. + cpr(1) + cpr(2) + cpr(3)
  pops(1) = 1./sum
  pops(2) = cpr(1)/sum
  pops(3) = cpr(2)/sum
  pops(4) = cpr(3)/sum
  popmet = (cpr(1)+cpr(2)+cpr(3))/sum
  If (ij==3) popmet = cpr(3)/sum
  fc = 1.E23*1.6E-12*cn/dene
  Go To (150, 160, 170, 180) ij
! ERGS * 10-23 / S / ATOM
150 pln(1, 1, no) = (pops(1)*rres(1,no)*v*exp(-11590.*eres(1, &
    no)/t)+popmet*r5)*fc*eres(1, no)
  pln(1, 2, no) = popmet*v*rbe(no)*exp(-11590.*bee(no)/t)*fc*bee(no)
  pln(1, 3, no) = pops(3)*a21(1, no)*emet*fc
  sumc = -c(1) - c(2) - c(3)
  reduc(1) = (pops(3)*a21(1,no)+pops(4)*a21(2,no))/sumc
  Return
160 pln(2, 1, no) = (pops(1)*rres(2,no)*v*exp(-11590.*eres(2, &
    no)/t)+popmet*r5)*fc*eres(2, no)
  pln(2, 3, no) = (pops(2)*b21(1,no)+pops(3)*b21(2,no)+pops(4)*b21(3,no))*fc* &
    emet
  sumc = -c(1) - c(2) - c(3)
  reduc(2) = (pops(2)*b21(1,no)+pops(3)*b21(2,no)+pops(4)*b21(3,no))/sumc
  Return
170 pln(3, 1, no) = (pops(1)*rres(3,no)*v*exp(-11590.*eres(3, &
    no)/t)+popmet*r5)*fc*eres(3, no)
  pln(3, 3, no) = pops(4)*c21(3, no)*fc*emet
  reduc(3) = -pops(4)*c21(3, no)/c(3)
  Return
180 pln(4, 1, no) = (pops(1)*rres(4,no)*v*exp(-11590.*eres(4, &
    no)/t)+popmet*mgm5)*fc*eres(4, no)
  pln(4, 3, no) = pops(3)*mg21(1, np)*fc*emet
  sumc = -c(1) - c(2) - c(3)
  reduc(4) = (pops(3)*mg21(1,np)+pops(4)*mg21(2,np))/sumc
  Return
190 Format (' BAD INVERSION')

End Function
Function noson(a, x, l, lmax)
! RUDOLF LOESER, 28 JULY 1965
! NOSON IS A MATRIX INVERSION ROUTINE, USING THE METHOD OUTLINED ON
! P 434 OF HILDEBRAND, INTRODUCTION TO NUMERICAL ANALYSIS,
! (NEW YORK,1956).
! A IS A MATRIX OF ORDER L, WITH COLUMN DIMENSION LMAX, ITS ELEMENTS
! ARE ASSUMED TO BE STORED COLUMNWISE, IN THE USUAL FORTRAN MANNER,
! X IS WORKING STORAGE OF LENGTH L.
! THE INVERSE OF A WILL REPLACE A.
! UPON RETURN, NOSON=1 IF INVERSION WENT OK, =0 IF A DIVISOR IS
! ZERO (IN WHICH CASE A MAY CONTAIN GARBAGE UPON RETURN).
  Dimension :: a(1), x(1)

  n = l - 1
  max = n*lmax + l
  Do i = 1, l
    x(i) = 1.
  End Do
  k1 = -lmax
  Do k = 1, l
    k1 = k1 + lmax
    k2 = k1 + k
    If (a(k2)==0.) Then
      Go To 170
    Else
      Go To 100
    End If
100 Do i = 1, l
      j1 = k1 + i
      If (a(j1)==0.) Then
        Go To 120
      Else
        Go To 110
      End If
110   f = 1./a(j1)
      x(i) = x(i)*f
      Do j1 = i, max, lmax
        a(j1) = a(j1)*f
      End Do
120 End Do
    a(k2) = x(k)
    x(k) = 1.
    Do i = 1, l
      ki = k - i
      If (ki==0) Then
        Go To 150
      Else
        Go To 130
      End If
130   j1 = k1 + i
      If (a(j1)==0.) Then
        Go To 150
      Else
        Go To 140
      End If
140   a(j1) = 0.
      Do j2 = i, max, lmax
        j1 = j2 + ki
        a(j2) = a(j2) - a(j1)
      End Do
150 End Do
  End Do
  Do i = 1, n
    If (x(i)==0.) Then
      Go To 170
    Else
      Go To 160
    End If
160 f = 1./x(i)
    Do j1 = i, max, lmax
      a(j1) = a(j1)*f
    End Do
  End Do
  noson = 1
  Return
170 noson = 0
  Return
End Function

Subroutine nquil(t, n, dene, jcont, iphot, idens, icx)
! ************************  METASTABLE POPULATION DATA
! July 93: CD1 and CD2 for Li-like ions work well for high Z, but
! underestimates density effects for low Z
  Dimension :: emetj(30, 4), cd1(28), cd2(28)
  Dimension :: rat(30), prod(30), c5(30), c6(30), over(30), imeta(30)

  Dimension :: sec(30)
  Common /dat/e(30), ea(30), s2(30), wave(220), e3(220), f(220), ll(30), &
    s3(30), s4(30), s5(30)
  Common /params/nj(12), abunj(12), abund, binmin, binsyz, nbin
  Common /rates/c1(30), c2(30), cth

  Common /result/conce(30), ca(30), power(220), rhy, heneut, heplus, dne, &
    pcool, pou, pot, re, tu, pm(4)
  Data emetj/5*0., 6.5, 8.35, 10.2, 0., 14., 0., 17.6, 0., 21.6, 0., 25.3, 0., &
    29.0, 0., 32.7, 5*0., 47.0, 0., 51.8, 0., 0., 5*0., 5.34, 7.10, 8.87, 0., &
    12.5, 0., 16.0, 0., 19.8, 0., 24.0, 0., 28.6, 0., 34.2, 5*0., 57.1, 0., &
    64.8, 0., 0., 5*0., 4.18, 5.78, 7.48, 0., 11.0, 0., 14.7, 0., 18.6, 0., &
    23.0, 0., 28.0, 0., 33.9, 5*0., 60.4, 0., 69.2, 0., 0., 11*0., 2.72, 0., &
    6.55, 0., 10.4, 0., 14.2, 0., 18.1, 5*0., 29.7, 0., 38.9, 2*0./

  Data imeta/0, 0, 0, 1, 2, 3, 5*0, 4, 18*0/
  Data cd1/.0024, -.0005, .00, .061, .027, .011, .005, .005, 0., .0107, .09, &
    .13, .11, .081, .075, .066, .051, 11*0./
  Data cd2/0., .01485, .30, .108, .107, .024, .075, .051, .054, .0167, .36, &
    17*0./
  Data sec/4*1., .9, .9, 6*1., .6, .5, 1., .7, .9, 13*1./

  ix = 0
  abhe = 10.**(abunj(1)-12.)
  Do i = 1, 4
    pm(i) = 0.
  End Do

  orange = 1.
  c5(1) = 0.
  c5(n+1) = 0.
  jlow = 1
  jmax = n + 1
  jtop = n + 1
  c6(1) = 0
  t6 = t*1.0E-6
  st = sqrt(t)
  Do j = 1, n
!   IONIZATION RATE FROM YOUNGER
    eion = e(j)
    c1(j) = cion(n, j, eion, t)
    If (n-j<=1) Go To 110
!   IONIZATION FROM METASTABLE LEVEL : BE,B,C,MG SEQUENCES
    If (idens==0) Go To 100
    imet = imeta(n-j+1)
    If (imet==0) Go To 100
    emet = emetj(n, imet)
    fmet = popmet(n, imet, emet, dene, t)
    pm(imet) = fmet
    em = e(j) - emet
    c1(j) = c1(j)*(1.-fmet) + fmet*cion(n, j, em, t)
100 Continue
!   INNERSHELL IONIZATION
    If (n-j<=1) Go To 110
    eion = ea(j)
    jp = n - 1
    If (n-j>=10) jp = n - 9
    c1(j) = c1(j) + cion(n, jp, eion, t)
110 Continue
!   IONIZATION FOLLOWING INNERSHELL EXCITATION  :  COWAN AND MANN
!   INCREASED FOR FE
    c1(j) = c1(j) + autoin(n, j, t)
!   ****************      PHOTOIONIZATION
    If (iphot/=0) c1(j) = c1(j) + phot(n, j, e(j), t, dene)/dene
!   ****************      RADIATIVE RECOMBINATION
    zeff = j
    beta = zeff*zeff/(6.34*t6)
    c2(j+1) = 2.07E-11*zeff*zeff*(0.4288+0.5*alog(beta)+0.469*(beta**(-1./3.)) &
      )/sqrt(t)
    apple = c2(j+1)
    plum = 0
    recn = 0.
    rgnd = 0.
    If (n==j) Go To 130
    nval = gnd(n-j) + 0.1
    starn = sqrt(13.6/e(j))*zeff
    Do nqm = 1, nval
      qm = nqm
      be = beta/(qm*qm)
      c2(j+1) = c2(j+1) - 5.197E-14*j*sqrt(be)*seaton(be, qm)
    End Do

    If (c2(j+1)<0) c2(j+1) = 0
    i = n - j
!   RECOMBINATION TO GROUND STATE
    If (i<=17) rgnd = grec(n, j, e(j), t)
    If (i<=1) Go To 130
    If (i>=4 .And. i<=9) Go To 130
    ei = e(j)
    If (i>=17) Go To 120
!   RECOMBINATION TO OTHER STATES IN SAME SHELL
    lion = 1
    If (i==12) lion = 3
    ei = ei - e3(ix+lion)
120 starn = sqrt(13.595/ei)*zeff
    be = beta/(starn*starn)
    recn = effnew(n-j, t, zeff, n)*seaton(be, starn)*2.599E-14*j* &
      sqrt(ei*11590./t)/(starn*starn)
130 c2(j+1) = c2(j+1) + recn + rgnd
    If (beta>=0.5) Go To 140
    chi = 1.2*beta + (0.55+1.04*alog(beta))*beta*beta + &
      (-0.43+1.01*alog(beta))*beta**3
    Go To 150
140 chi = 0.5*(0.735+alog(beta)+1./(3.*beta))
150 Continue
    Do nqm = 1, nval
      qm = nqm
      chi = chi - p(beta, qm)
    End Do
    chi = chi + effn(n-j, zeff, t)*p(beta, starn)/(2.*starn*starn)
    c6(j+1) = e(j)*1.6027E-12*(c2(j+1)+4.22E-22*sqrt(t/1.037E-11)*chi)
!   DIELECTRONIC RECOMBINATION ***** DENSITY DEPENDENCE AND SECONDARY
!   ***** AUTOIONIZATION
    iy = ll(j)*3
    If (j==1) Go To 170
    zeff = j - 1
    rho = (dene/zeff**7.)**.2
    dd = 1.
    If (rho>=3) dd = 1./(1.+(cd1(n-j+1)+cd2(n-j+1)/zeff)*rho)
    Do l = 1, iy
      ln = ix + l
      plum = plum + alphadi(n, j, l, ln, t)*amin1(dd, sec(n-j+1))
    End Do
!   DIELECTRONIC RECOMBINATION FOR METASTABLE LEVELS OF BE,B,C,MG ISO
    imet = imeta(n-j+1)
    If (imet==0) Go To 160
    zf = j - 1
    b = sqrt(zf)*(zf+1.)**2.5/sqrt(zf*zf+13.4)
    If (idens>0) plum = plum*(1.-fmet)
    If (idens>0) plum = plum + fmet*dimet(n, j, t, b, dd)
160 Continue
!   *************** ADD DIELECTRONIC RECOMB AND STOREY'S LOW T DIELEC RE
    c2(j) = c2(j) + plum + alflo(n, j, t)
170 over(j) = plum/orange
    If (j/=1) c5(j) = plum*e(j-1)*10.**(abund-1.)*1.6027
    ix = ix + iy
    orange = apple
  End Do

  If (icx==0) Go To 180
! CHARGE TRANSFER CONTRIBUTION
  fplus = rhy/(1.+rhy)
  denh = 1./(.003+fplus+abhe*heplus+abhe*2.*(1.-heneut-heplus))
  c1(n+1) = 0.
  c2(1) = 0.
  nn = n + 1
! **********************************  S. Butler Charge Transfer Rates
  Do j = 1, n
    c1(j) = c1(j) + cthi(n, j, t)*denh*fplus + cthei(n, j, t, dene)*heplus* &
      abhe
    c2(j+1) = c2(j+1) + cthr(n, j+1, t)*denh*(1.-fplus) + &
      cther(n, j+1, t)*heneut*abhe
  End Do
180 Continue
! AUGER IONIZATION
  If (iphot==2) Call aphot(n, dene, jcont)
  If (jcont==1) Go To 260
  Do j = 1, n
    c2(j+1) = abs(c2(j+1))
    If (c1(j)/c2(j+1)<=0.) Go To 190
    rat(j) = c2(j+1)/c1(j)
    If (j==1 .And. t>=1.0E8) rat(j) = amin1(rat(j), 1.0E4)
    Go To 200
190 rat(j) = 1.0E+6
200 Continue
    If (rat(j)>=1.0E-5) Go To 210
    jlow = j + 1
210 If (rat(j)<1.0E+5) Go To 220
    jmax = j
    Go To 230
220 End Do
  If (jlow==jmax) Go To 240
230 jump = max(jmax-1, 1)
  prod(jump) = rat(jump)
  ktop = jump - jlow
  sum = 1.00000 + prod(jump)
  If (ktop==0) Go To 250
  Do k = 1, ktop
    prod(jump-k) = rat(jump-k)*prod(jmax-k)
    prod(jump-k) = amin1(prod(jump-k), 1.0E30)
    sum = sum + prod(jump-k)
  End Do

  Go To 250
240 sum = 1.0000
250 Continue
  Do j = 1, jtop
    conce(j) = 0
    c5(j) = 0.
  End Do

  conce(jmax) = 1.000/sum
  kmax = jmax - jlow
  If (kmax==0) Go To 260
  Do k = 1, kmax
    conce(jmax-k) = prod(jmax-k)*conce(jmax)
  End Do

260 Continue
  If (n==2) dne = 0.
  worl = 0.
  dei = 0.
  Do j = 1, jtop
    dne = (j-1)*conce(j)*10.**(abund-12.) + dne
    c6(j) = c6(j)*conce(j)*10.0**(11.+abund)
    re = re + c6(j)
    c5(j) = c5(j)*conce(j)
    pcool = pcool + c5(j)
    pcool = pcool + c6(j)
  End Do
  Return
End Subroutine

Subroutine heseq(n, j, t, dene, ix)
! HELIUM-LIKE IONS  :  RECOMBINATION TO EXCITED LEVELS, INNER-SHELL
! IONIZATION OF LI-LIKE IONS  AND DENSITY DEPENDENCE USING PRADHAN,
! MEWE & SCHRIJVER,  BERRINGTON, FON & KINGSTON, DRAKE

  Common /dat/e(30), s(30), c(30), wave(220), e3(220), f(220), ll(30), &
    sig(30), alf(30), es(30)
  Common /result/conce(30), gndrec(30), power(220), rhy, heneut, heplus, dne, &
    pcool, pou, pot, re, tu, pm(4)
  Common /params/nj(12), abunj(12), abund, binmin, binsyz, nbin
  Common /twosv/hetwop

  Dimension :: br(30), cr(3, 30)
  Data br/0., 1.0, 3*0., .89, .78, .71, 0., .67, 0., .63, 0., .58, 0., .50, &
    0., .40, .0, .33, 5*0., .20, 0., .18, 0., 0./
  Data cr/3*0., 7.6E6, 4540., 2210., 9*0., 1.9E10, 1.3E11, 9.6E11, 4.9E10, &
    6.5E11, 2.3E12, 1.25E11, 3.3E12, 5.8E12, 3*0., 2.4E12, 8.4E13, 1.4E14, &
    3*0., 2.3E13, 1.1E15, 1.9E15, 3*0., 2.3E14, 1.5E15, 2.7E16, 3*0., 9.5E14, &
    7.6E16, 1.3E17, 3*0., 3.9E15, 3.9E17, 6.2E17, 3*0., 1.6E16, 2.0E18, &
    3.0E18, 15*0., 4.5E17, 5.3E19, 8.2E19, 3*0., 1.0E18, 1.0E20, 2.0E20, 6*0./

  If (conce(j)<.00001) Go To 110
  abu = 10.**(abund-12.)
  tu = 0.
  en = n
! EXCITATION OF 1S2S 1S LEVEL PRADHAN ; BALUJA,CALLAWAY,HENRY: 1.2 FOR
  y = 11590.*e3(ix+1)/t
  cc = alog((y+1.)/y) - 0.4/(y+1)**2
  gbar = (.06168+.0002522*n) + (-.02404-.001948*n)*y*cc + &
    (-.007604+.002416*n)*(y-y*y*cc)
! GBAR = (.03014+.0009855*N)+(.2392-.007019*N)*Y*CC
! 1  +(-.1292+.003543*N)*(Y-Y*Y*CC)
  If (n==2) gbar = .048*(t/10000.)**(.201)
  om = 1.2*14.5*f(ix+1)*gbar*13.6/e3(ix+1)
  p2ph = (8.63E-6*om/sqrt(t))*exp(-y)*abu*conce(n-1)*e3(ix+1)*1.6E11

  tu = p2ph
! INNERSHELL IONIZATION OF LI-LIKE ION FROM MEWE AND SCHRIJVER
  If (n==2) Go To 100
  eis = 13.6*en*en*3.**(.24-.43/alog10(en))
  cionize = 2.5E-8*sqrt(t)*exp(-eis*11590./t)/eis**2
  addit = conce(j-1)*cionize*e3(ix+2)*1.6E11*abu
  power(ix+2) = power(ix+2) + addit*.75

  p2ph = p2ph + addit*.25
100 Continue
! RADIATIVE RECOMBINATION
  cn = conce(j+1)*abu
  z = n - 1.
  power(ix+2) = power(ix+2) + (5.9E-13*z**1.8*t**(-0.4))*(1.+17*z**0.4*t**( &
    -0.2))*e3(ix+2)*1.6E11*cn
  power(ix+3) = power(ix+3) + (3.6E-11*z**2.4*t**(-0.7)+3.6E-10*z**2.8*t**(-.9 &
    ))*1.6E11*e3(ix+3)*cn
  power(ix+1) = power(ix+1) + 1.6E11*cn*e3(ix+1)*1.4E-11*z**2.52*t**(-.76)*(1. &
    +10.*z**0.4*t**(-0.2))

  p2ph = p2ph + cn*1.6E11*e3(ix+1)*(5.4E-13*z*z*t**(-0.5))*(1.+17.*z**(0.4)*t &
    **(-0.2))

  hetwop = p2ph
  zf = n - 1
  effecn = sqrt(zf*zf*13.6/(e(j)-e3(ix+7)))
  xn = 157890.*zf*zf/(t*effecn*effecn)
  recrat = 5.2E-14*zf*sqrt(xn)*seaton(xn, effecn)
  addit = (5./12.)*recrat*cn*(12399./wave(ix+7))*1.6E11
! HYDROGENIC RECOMBINATION FOR 3D LINES EXCEPT HELIUM FROM ROBBINS
  If (n==2) addit = cn*5.00E-14*(.0001*t)**(-1.333)*1.6E11*12399./wave(ix+7)
  power(ix+7) = power(ix+7) + addit
  power(ix+8) = power(ix+8) + .3333*addit*wave(ix+7)/wave(ix+8)
! DIELECTRONIC RECOMB
  power(ix+1) = power(ix+1) + 1.6E11*e3(ix+1)*cn*adh(1, n, t)
  power(ix+2) = power(ix+2) + 1.6E11*e3(ix+2)*cn*adh(2, n, t)
  power(ix+3) = power(ix+3) + 1.6E11*e3(ix+3)*cn*adh(3, n, t)

  p2ph = p2ph + 1.6E11*e3(ix+1)*cn*adh(4, n, t)
! DENSITY DEPENDENCE :  PRADHAN; BERRINGTON, FON & KINGSTON FOR HE I
  crit1 = cr(1, n)
  crit2 = cr(2, n)
  crit3 = cr(3, n)
  critm = 1./(1./crit1+1./crit2+1./crit3)
  crit = 5000.*z**8*sqrt(t)
  If (n==2) crit = 3.1E5*sqrt(t)
  p2dum = p2ph

  p2ph = (p2ph+power(ix+2)*(critm/crit3)*dene/(dene+critm))*crit/(dene+crit)
  power(ix+1) = power(ix+1) + p2dum*dene/(dene+crit) + &
    power(ix+2)*(critm/crit2)*dene/(dene+critm)
  branch = br(n)
  pdum = power(ix+3)
  power(ix+3) = power(ix+3)*(1.-branch) + power(ix+2)*(critm/crit1)*(dene/( &
    dene+critm))
  fc = branch*(1.-branch)
  If (n==2) fc = 2.03E-5
  power(ix+6) = pdum*branch*wave(ix+3)/wave(ix+6) + power(ix+2)*(wave(ix+2)/ &
    wave(ix+6))*fc*(dene/crit1)/(1.+dene/critm)
  power(ix+2) = (power(ix+2)+branch*pdum)*critm/(dene+critm)

  wav2 = wave(ix+1)
  If (nbin/=0) Call twoph(wav2, p2ph, 0)
110 Continue
  Return

End Subroutine
Function adh(l, n, t)
! DIELECTRONIC RECOMBINATION TO HE-LIKE EXCITED STATES
! USING MEWE&SCHRIJVER
  p = 33*(n-1.)**0.6*t**(-0.3)
  z = n
  z2 = z*z
  z3 = z*z2
  z4 = z*z3
  zpt = (z+.5)*(z+.5)
  a = 6.46E-8*z4*t**(-1.5)
  ex1 = exp(-78900.*zpt/t)
  ex2 = exp(-101800.*z2/t)
  ex3 = exp(-118400.*z2/t)
  Go To (100, 110, 120, 130) l
100 adh = a*(12.*ex1/(1.+6.E-6*z4)+18.*ex2/(1.+3.0E-5*z4)+69.*ex3/(1.+5.0E-3* &
    z3))
  Return
110 adh = a*(9.*ex1/(1.+7.E-5*z4)+27.*ex2/(1.+8.0E-5*z4)+380.*ex3/((1.+p)*(1.+ &
    5.0E-3*z3)))
  Return
120 adh = a*(18.*ex1/9.5+54*ex2/(1.+1.9E-4*z4)+380.*ex3*p/((1.+p)*(1.+ &
    5.0E-3*z3)))
  Return
130 adh = a*(3.*ex1/(1.+3.0E-6*z4)+0.5*ex2/(1.+2.2E-5*z4)+6.3*ex3/(1.+5.0E-3* &
    z3))
  Return

End Function

Subroutine hyseq(n, j, t, dene, ix)
! HYDROGENIC IONS : ADDS RECOMBINATION, DENSITY DEPENDENCE
! HAYES AND SEATON 2S EXCITATION,  HYDROGENIC RECOMBINATION
  Common /dat/e(30), s(30), c(30), wave(220), e3(220), f(220), ll(30), &
    sig(30), alf(30), es(30)
  Common /result/conce(30), gndrec(30), power(220), rhy, heneut, heplus, dne, &
    pcool, pou, pot, re, tu, pm(4)
  Common /params/nj(12), abunj(12), abund, binmin, binsyz, nbin

  If (conce(j)<.00001) Go To 100
  abu = 10.**(abund-12.)
  en = n
  cn = conce(j)*abu
  y = e3(ix+1)*11590./t
  cc = alog((y+1)/y) - 0.4/(y+1)**2
  gbar = 0.047
! GBAR = 0.20 * Y * CC + 0.276 * CC
  If (n==1) gbar = .017 + .1*cc
! 1.2 FOR CASCADES : TU IS DIRECT EXCITATION TWO PHOTON
  om = 1.2*14.5*.416*gbar*13.6/e3(ix+1)
  p2ph = (8.63E-6*om/sqrt(t))*exp(-y)*e3(ix+1)*1.6E11*cn

  tu = tu + p2ph
! RADIATIVE RECOMBINATION : CASE A : exponent revised Aug. '86
  tp = .0001*t/n**2
  cn = conce(j+1)*abu
  rc = n*16.7E-14*tp**(-0.91)
  power(ix+1) = power(ix+1) + rc*cn*e3(ix+1)*1.6E11
  rc = n*3.61E-14*tp**(-0.72)
  power(ix+2) = power(ix+2) + rc*cn*e3(ix+2)*1.6E11
  rc = n*1.40E-14*tp**(-0.757)
  power(ix+3) = power(ix+3) + rc*cn*e3(ix+3)*1.6E11
  rc = n*.693E-14*tp**(-.750)
  power(ix+4) = power(ix+4) + rc*cn*e3(ix+4)*1.6E11
  rc = n*7.84E-14*tp**(-1.04)
  If (n==2) rc = 2*11.8E-14*tp**(-1.02)
! CASE B FOR HELIUM PLUS
  power(ix+5) = power(ix+5) + rc*cn*1.6E11*12399./wave(ix+5)
  rc = n*3.18E-14*tp**(-.610)

  p2ph = p2ph + rc*cn*1.6E11*e3(ix+1)
! DENSITY DEPENDENCE OF METASTABLE LEVEL MEWE AND GRONENSCHILD
  crit = 100.*en**8*sqrt(t)
  power(ix+1) = power(ix+1) + p2ph*dene/(dene+crit)
  p2ph = p2ph*crit/(dene+crit)
  wav2 = wave(ix+1)
  If (nbin/=0) Call twoph(wav2, p2ph, 0)
  power(ix+5) = power(ix+5) + power(ix+2)*.11*wave(ix+2)/wave(ix+5)
  power(ix+2) = .89*power(ix+2)
  power(ix+3) = .86*power(ix+3)
  power(ix+4) = .78*power(ix+4)
100 Continue
  Return
End Subroutine
Subroutine recems(t, n, iprint, jprint)
  Real :: kt

  Logical :: notgrc(25), done(25), switch
  Dimension :: drat(20), lion(2, 20), zeta(2, 20), iedge(26), edge(25), &
    chi(25), esmin3(25), alfx(25), jedge(25), imap(25)

  Dimension :: s2x(25), s3x(25), s4x(25), s5x(25)
  Common /dat/e(30), ea(30), s2(30), wave(220), e3(220), f(220), ll(30), &
    s3(30), s4(30), s5(30)
  Common /params/nj(12), abunj(12), abund, binmin, binsyz, nbin
  Common /rates/c1(30), c2(30), cth
  Common /result/conce(30), gndrec(30), power(220), rhy, heneut, heplus, dne, &
    pcool, pou, pot, re, tu, pm(4)

  Common /contin/brmev(1000), recev(1000), tuf(1000)

  ibn(z) = int((z-binmin)/binsyz+0.5) + 1
  ecen(i) = binmin + (i-0.5)*binsyz

  gfunc(x, y) = 1. + 0.1728*(y*x)**(-2./3.)*(x-2.) - &
    0.0496*(y*x)**(-4./3.)*(x*x-x*2./3.+2./3.)
  Data drat/2., .5, 2., .5, 6., 2.5, 2.22, 2.25, .667, .167, 2., .5, 6., 2.5, &
    2.22, 2.25, .667, .167, 2., .5/
  Data lion/0, 0, 3, 7, 1, 5, 2, 4, 4, 0, 4, 0, 5, 0, 2, 0, 2, 0, 3, 0, 1, 5, &
    1, 4, 3, 5, 1, 3, 1, 4, 1, 3, 1, 3, 1, 2, 4, 0, 6, 0/

  Data zeta/0., 0., 2., 3., 2., 3., 1.75, 3., 3., 0., 3., 0., 3., 0., 3., 0., &
    3., 0., 3., 0., 3., 4., 2.83, 4., 2.67, 4., 2.5, 4., 2.33, 4., 2.17, 4., &
    2., 4., 1.83, 4., 4., 0., 4., 0./
  t6 = t/1.E6
  kt = 86.17*t6
! UPPER LIMIT TO BINNING IS CHOSEN WHERE EXP = 1.E-20
  imax = nbin
  If (imax<1) Return
  ab = 10.**(abund-12.)
  ebk = exp(binsyz/kt)
  q = ab*kt*(ebk-1.)*exp(-binmin/kt)/t6**1.5
  qprime = ab*kt/t6**1.5
  qgrec = q*1.31E-8
  qtk = q*6.52/12399.
  z = n
  qhyd = q*5.28E-4*(z**4)
  qgrecpr = qprime*1.31E-8
  qtkpr = qprime*6.52/12399.
  qhydpr = qprime*5.28E-4*(z**4)
  ebk = 1./ebk
! LIMIT TO SIGNIFICANT CONTRIB FROM AN EDGE; ASSUMES #IONS PRESENT
! GOES LIKE SQRT(Z); COMPARE THRESH RECEMS TO EDGLIM*(CURRENT CONTIN)
  edglim = 0.01/sqrt(z)
! FIRST CONSTRUCT TABLE OF SIGNIFICANT EDGES
  ix = 0
  j = 1
  iso = n
100 If (iso<=20) Go To 110
  ix = ix + 3*ll(j)
  iso = iso - 1
  j = j + 1
  Go To 100
110 iemax = 0
120 If (iso==1) Go To 160
  If (conce(j+1)<1.E-4) Go To 150
  chitry = e(j)
  ignd = ibn(chitry)
  If (ignd>=nbin) Go To 130
  If (ignd<1) ignd = 1
  thresh = s2(j) + s3(j) + s4(j) + s5(j)

  If (chitry/kt>46.) recev(ignd) = recev(ignd) + qgrecpr*conce(j+1)*drat(iso)* &
    (chitry**3)*thresh
  If (chitry/kt>46.) Go To 130

  edgtry = qgrec*conce(j+1)*drat(iso)*(chitry**3)*thresh*exp(chitry/kt)
  en = ecen(ignd)
  x = chitry/en
  x2 = x*x
  edg = edgtry*(ebk**ignd)*(s2(j)/x+s3(j)+s4(j)*x+s5(j)*x2)/thresh
  binnow = brmev(ignd) + recev(ignd)
  If (edg<edglim*binnow) Go To 150
  iemax = iemax + 1
  iedge(iemax) = ignd
  edge(iemax) = edgtry
  chi(iemax) = chitry
  jedge(iemax) = j
  notgrc(iemax) = .False.
  s2x(iemax) = s2(j)/thresh
  s3x(iemax) = s3(j)/thresh
  s4x(iemax) = s4(j)/thresh
  s5x(iemax) = s5(j)/thresh
  If (iemax==25) Go To 180
130 ln = ll(j)*3
  If (ln==0) Go To 150
  Do k = 1, 2
    l = lion(k, iso)
    If (l==0 .Or. l>ln) Go To 140
    e3x = e3(ix+l)
    If (e3x<1.) Go To 140
    ixl = ix + l
    chitry = e(j) - e3(ixl)
    ignd = ibn(chitry)
    If (ignd>=nbin) Go To 140
    If (ignd<1) ignd = 1

    If (chitry/kt>46.) recev(ignd) = recev(ignd) + qtkpr*conce(j+1)*zeta(k, &
      iso)*((chitry/13.6)**2)
    If (chitry/kt>46.) Go To 140

    edgtry = qtk*conce(j+1)*zeta(k, iso)*((chitry/13.6)**2)*exp(chitry/kt)
    en = ecen(ignd)
    edg = edgtry*ebk**ignd
    binnow = brmev(ignd) + recev(ignd)
    If (edg<edglim*binnow) Go To 140
    iemax = iemax + 1
    iedge(iemax) = ignd
    edge(iemax) = edgtry
    chi(iemax) = chitry
    jedge(iemax) = j
    notgrc(iemax) = .True.
    done(iemax) = .True.
    If (iemax==25) Go To 180
140 End Do
150 ix = ix + 3*ll(j)
  j = j + 1
  iso = iso - 1
  Go To 120
160 If (conce(n+1)<1.E-6) Go To 190
  If (conce(n+1)<1.E-4 .And. n/=2) Go To 190
  Do k = 1, 2
    y = k
    chitry = 13.6*(z/y)**2
    ignd = ibn(chitry)
    If (ignd>=nbin) Go To 170
    If (ignd<1) ignd = 1
    en = ecen(ignd)
    x = en/chitry
    gftr = gfunc(x, y)
    If (x>100.) gftr = 7./sqrt(x)
    gmxtry = 1.1 - 1./(10.*y**1.5)
    If (gftr>gmxtry) gftr = gmxtry

    If (chitry/kt>46.) recev(ignd) = recev(ignd) + qhydpr*conce(n+1)*gftr/y**3
    If (chitry/kt>46.) Go To 170

    edgtry = qhyd*conce(n+1)*exp(chitry/kt)/y**3
    edg = edgtry*gftr*ebk**ignd
    binnow = brmev(ignd) + recev(ignd)
    If (edg<edglim*binnow) Go To 170
    iemax = iemax + 1
    iedge(iemax) = ignd
    edge(iemax) = edgtry
    chi(iemax) = chitry
    jedge(iemax) = j
    esmin3(iemax) = gmxtry
    alfx(iemax) = y
    notgrc(iemax) = .True.
    done(iemax) = .False.
    If (iemax==25) Go To 180
170 End Do
  Go To 190
180 Print 310, n
190 If (iprint>0) Print 320, n, t, iemax
  If (iemax==0) Return
! BUBBLE SORT TO PUT IEDGE IN INCREASING ORDER; PARALLEL SORT OF IMAP
  Do ie = 1, 25
    imap(ie) = ie
  End Do

200 switch = .False.
  iedge(iemax+1) = 0
  If (iemax<25) iedge(iemax+2) = 0
  If (iemax==1) Go To 230
  ie = 1
210 iesave = iedge(ie)
  If (iedge(ie+1)>=iesave) Go To 220
  imsave = imap(ie)
  iedge(ie) = iedge(ie+1)
  imap(ie) = imap(ie+1)
  iedge(ie+1) = iesave
  imap(ie+1) = imsave
  If (ie>1) switch = .True.
220 ie = ie + 1
  If (ie<iemax) Go To 210
  If (switch) Go To 200
! END SORT, PRINT TABLE OF EDGES SIGNIFICANT ENOUGH TO BE INCLUDED
230 If (jprint==0) Go To 260
  Print 330
  iep = iemax/2 + 1
  Do ie = 1, iep
    jb = ie*2
    ib = jb - 1
    jp = imap(jb)
    ip = imap(ib)
    If (iedge(ib)==0) Go To 250
    If (iedge(jb)==0) Go To 240
    Print 340, jedge(ip), chi(ip), edge(ip), notgrc(ip), done(ip), jedge(jp), &
      chi(jp), edge(jp), notgrc(jp), done(jp)
    Go To 250
240 Print 340, jedge(ip), chi(ip), edge(ip), notgrc(ip), done(ip)
250 End Do
! BEGIN BINNING, USING LIST IEDGE TO GOVERN NUMBER OF EDGES INCLUDED
260 imin = iedge(1)
  iemax = 1
  fimin = imin - 1
  exprod = ebk**fimin
  Do i = imin, nbin
270 If (i/=iedge(iemax+1)) Go To 280
    iemax = iemax + 1
    Go To 270
280 en = ecen(i)
    binsum = 0.
    Do ie = 1, iemax
      ipnt = imap(ie)
      binn = edge(ipnt)
      If (notgrc(ipnt)) Go To 290
      x = chi(ipnt)/en
      binn = binn*(s2x(ipnt)/x+s3x(ipnt)+s4x(ipnt)*x+s5x(ipnt)*x*x)
      Go To 300
290   If (done(ipnt)) Go To 300
      y = alfx(ipnt)
      x = en/chi(ipnt)
      gftr = 7./sqrt(x)
      If (x<100.) gftr = gfunc(x, y)
      gmax = esmin3(ipnt)
      If (gftr>gmax) gftr = gmax
      binn = binn*gftr
300   binsum = binsum + binn
    End Do

    exprod = exprod*ebk
    recev(i) = recev(i) + binsum*exprod
  End Do
  Return
310 Format (1X, /, /, ' IEMAX OVERRUN ON ELEMENT N =', I3, /, /)
320 Format (1X, /, /, ' RECEMS FOR ELEMENT N =', I3, '  AT T =', 1P, E10.3, &
    ',   WITH', I3, '  SIGNIF. EDGES', /)
330 Format (2(11X,'J     CHI',8X,'EDGE',7X,'N GRC   SHRT'), /)
340 Format (2(10X,I2,0P,F10.1,1P,E12.2,2L8))

End Subroutine
Function fbg(u, gam)
! 39th values of A1, A2, and A3 revised to get rid of 5% bump
! in Kellog et al fit to Karzas and Latter at gamma2 = 0.1,
! u = .3-1. : uses numbers from Larry Molnar, Jan 1988
! and corrections from Jack Hughes

  Double Precision :: t, ai, ak, u4
  Dimension :: a(6, 7, 3), a1(126), gam2(6), gam3(6)
  Equivalence (a, a1)
  Data gam2/.7783, 1.2217, 2.6234, 4.3766, 20., 70./
  Data gam3/1., 1.7783, 3., 5.6234, 10., 30./
  Data (a1(i), i=1, 42)/1.001, 1.004, 1.017, 1.036, 1.056, 1.121, 1.001, &
    1.005, 1.017, 1.046, 1.073, 1.115, .9991, 1.005, 1.030, 1.055, 1.102, &
    1.176, .9970, 1.005, 1.035, 1.069, 1.134, 1.186, .9962, 1.004, 1.042, &
    1.100, 1.193, 1.306, .9874, .9962, 1.047, 1.156, 1.327, 1.485, .9681, &
    .9755, 1.02009, 1.208, 1.525, 1.965/
  Data (a1(i), i=43, 84)/.30290, .16160, .04757, .01300, .00490, -.00320, &
    .49050, .21550, .08357, .02041, .00739, .00029, .65400, .28330, .08057, &
    .03257, .00759, -.00151, 1.0290, .39100, .12660, .05149, .01274, .00324, &
    .95690, .48910, .17640, .05914, .01407, -.00024, 1.2360, .75790, .32600, &
    .10770, .02800, .00548, 1.3270, 1.0170, 0.60166, .20500, .06050, .00187/

  Data (a1(i), i=85, 126)/ -1.3230, -.25400, -.01571, -.001000, -.000184, &
    .00008, -4.7620, -.33860, -.03571, -.001786, -.000300, .00001, -6.3490, &
    -.42060, -.02571, -.003429, -.000234, .00005, -13.231, -.59000, -.04571, &
    -.005714, -.000445, -.00004, -7.6720, -.68520, -.06430, -.005857, &
    -.000420, .00004, -7.1430, -.99470, -.12000, -.010070, -.000851, -.00004, &
    -3.1750, -1.1160, -.22695, -.018210, -.001729, .00023/

  gam1 = gam*1000.
  If (gam1>100.) Go To 120
  u2 = u**2

! *****COMPUTE BORN APPROXIMATION GAUNT_RS FACTOR

  u1 = u/2.
  t = u1/3.75
  u4 = u1/2.
  If (u1>2.) Go To 100
  ai = 1.0 + 3.5156229*t**2 + 3.0899424*t**4 + 1.2067492*t**6 + &
    0.2659732*t**8 + 0.0360768*t**10 + 0.0045813*t**12
  ak = -1.*log(u4)*ai - .57721566 + .42278420*u4**2 + .23069758*u4**4 + &
    .0348859*u4**6 + .00262698*u4**8 + .00010750*u4**10 + .0000074*u4**12
  Go To 110

100 ak = 1.2533141 - .07832358/u4 + .02189568/u4**2 - .01062446/u4**3 + &
    .00587872/u4**4 - .00251540/u4**5 + .00053208/u4**6
  ak = ak/(exp(u1)*sqrt(u1))
110 born = .5513*exp(u1)*ak

! *****COMPUTE POLYMONIAL FACTOR TO MULTIPLY BORN APPROXIMATION

  If (gam1<1.) Go To 130
  If (u<.003) Go To 130
  If (u<=.03) n = 1
  If ((u<=.3) .And. (u>.03)) n = 2
  If ((u<=1.) .And. (u>.3)) n = 3
  If ((u<=5.) .And. (u>1.)) n = 4
  If ((u<=15.) .And. (u>5.)) n = 5
  If (u>15.) n = 6
  If (gam1<=1.7783) m = 1
  If ((gam1<=3.) .And. (gam1>1.7783)) m = 2
  If ((gam1<=5.6234) .And. (gam1>3.)) m = 3
  If ((gam1<=10.) .And. (gam1>5.6234)) m = 4
  If ((gam1<=30.) .And. (gam1>10.)) m = 5
  If ((gam1<=100.) .And. (gam1>30.)) m = 6
  m1 = m + 1
  g1 = (a(n,m,1)+a(n,m,2)*u+a(n,m,3)*u2)*born
  g2 = (a(n,m1,1)+a(n,m1,2)*u+a(n,m1,3)*u2)*born
  p = (gam1-gam3(m))/gam2(m)
  fbg = (1.0-p)*g1 + p*g2
  Return
120 power = -.134/(gam**.2097)
  fbg = 1.5*(3.*u)**power
  Return
130 fbg = born
  Return

End Function
Subroutine brems(t, n)
  Real :: kt
  Common /result/conce(30), gndrec(30), power(220), rhy, heneut, heplus, dne, &
    pcool, pou, pot, re, tu, pm(4)
  Common /params/nj(12), abunj(12), abund, binmin, binsyz, nbin
  Common /contin/brmev(1000), recev(1000), tuf(1000)
  Dimension :: x(30)

  hc = 12399.
  j2 = n + 1
  j1 = 2
  t6 = t/1.E6
  kt = 86.17*t6
  ebk = exp(binsyz/kt)
  q = 10.**(abund-12.)*86.17*sqrt(t6)*20.4*(ebk-1.)*exp(-binmin/kt)/hc
  ebk = 1./ebk
  Do j = j1, j2
    x(j) = 0.
    If (conce(j)>1.E-5) x(j) = conce(j)*q*(j-1)**2
  End Do
  i2 = min0(nbin, int(46.*kt/binsyz)+1)
  exprod = 1.
  Do i = 1, i2
    en2 = binsyz*i + binmin
    en1 = en2 - binsyz
    y = (en2+en1)*.5/kt
    chrgsm = 0.
    Do j = j1, j2
      If (x(j)==0.) Go To 100
      z = (j-1)**2*.158/t6
      chrgsm = chrgsm + x(j)*fbg(y, z)
100 End Do
    exprod = exprod*ebk
    brmev(i) = brmev(i) + chrgsm*exprod
  End Do

  Return

End Subroutine
Subroutine twoph(wav2, p2ph, iii)
  Common /params/nj(12), abunj(12), abund, binmin, binsyz, nbin
  Common /contin/brmev(1000), recev(1000), tufev(1000)

  If (p2ph<=1.E-10) Return
  itmax = (12399./wav2-binmin)/binsyz
  imax = min0(itmax, nbin)
  If (imax<=1) Return
  Do ik = 1, imax
    r = (binmin+binsyz*(ik-.5))*wav2/12399.
    tufev(ik) = tufev(ik) + 12.*p2ph*r*r*(1.-r)*wav2*binsyz/12399.
  End Do
  Return
End Subroutine


Function gaunt_rs(t, e, n, j, l, dene)
! MAY '82 VERSION RELYING HEAVILY ON MANN AND ROBB CALCULATIONS,
! BHATIA FOR DELTA N = 0 AND BHATIA AND MASON DELTA N = 1
! DEc 17, 1991 : Modify to improve trans 1 of Li-like ions
  Dimension :: ah(5), bh(5), ghe(24), gbe(8, 5)
  Dimension :: gbe1(4), gbe2(4), gbe3(11), gbe4(11), gbe5(11), bespec(4, 3)
  Dimension :: gb(8, 6), gc(8, 9), gn(8, 4), gox(8, 8), gf(8, 9), gne(5, 12)
  Dimension :: gna(8, 2), gna4(3, 3), gmgii(12)
  Dimension :: gmg(10), gmgp(8, 2)
  Dimension :: gal(10), gal2(10), a(10), b(10), c(10)
  Dimension :: gmshell(6, 4), gk(10), gca(10)

  Dimension :: a1(12), b1(12), c1(12)
  Common /interc/reduc(4)

  Common /result/conce(30), gndrec(30), power(220), rhy, heneut, heplus, dne, &
    pcool, pou, pot, re, tu, pm(4)
  Data ah/0., .08, .08, .08, 0./
  Data bh/.22, .16, .19, .20, 0./
  Data ghe/.03785, .0002038, .04214, -.001901, -.03263, .004462, .28, 0., &
    .02327, -.0007424, -.06369, .001865, .07711, .000058, 0., 0., -.03168, &
    .0008439, .182, -.007289, -.05873, .009605, 0., 0./
  Data gbe/.7114, .00558, -.9277, -.00153, .2711, .00523, 0.0, 0.0, 0.0, 0.0, &
    -.044, .00735, .482, -.01855, 0.0, 0.0, .44, 0.0, -.010, 5*0.0, .077, &
    .00666, 6*0., -.2943, -.00084, .2574, .01167, .04852, -.00777, .2671, &
    .00555/
  Data gbe1/ -.0264, 0.03, 0.0, 0.20/
  Data gbe2/.0881, 0.10, .117, 0.10/
  Data gbe3/0.0, .064, .069, .048, .035, .0303, .025, .02, .02, .032, .032/
  Data gbe4/0.0, .342, .39, .192, .103, .0836, .063, .043, .043, 0.0, 0.0/
  Data gbe5/0.0, 3*0.0, 55000., 64000., 86000., 97000., 108000., 156000., &
    168000./
  Data bespec/.00069, .122, .0059, .0063, .0016, .13, .0073, .0058, .0017, &
    .127, .0087, .0040/
  Data gb/.2392, .00543, .3723, .0335, .0321, -.0253, .4591, -.00232, .2317, &
    .00654, -.1479, .03186, .3038, -.02130, .3904, -.00302, .0811, .01318, &
    -.2523, .04591, .2122, -.02192, .2304, .00639, -.112, .0063, .050, .011, &
    .0808, -.005, .239, .0023, .126, 0.0, .287, 0.0, 0.0, 0.0, .28, 0.0, .128, &
    -.025, -.46, .133, .138, -.023, -.70, .128/
  Data gc/.3539, .00039, .2314, .02314, .0358, -.01534, .5449, -.00858, .2005, &
    .00395, -.3012, .04332, -.09599, .01606, .3727, -.001517, .3229, -.002883, &
    .0212, .02481, .0783, -.01326, .4671, -.008106, -.112, .0063, .050, .011, &
    .0808, -.005, .239, .0023, .128, -.025, -.46, .133, .138, -.023, -.70, &
    .128, -.269, -.00117, .0318, .0233, -.0102, -.0075, .235, .0172, .22, &
    .0056, 6*0., .198, .0069, 6*0., 8*0./
  Data gn/ -.01876, .0011, .0083, .00182, .0135, -.00083, .040, .00040, .51, &
    -.1, -1.84, .53, .55, -.091, -2.80, .51, .128, -.025, -.46, .133, .138, &
    -.023, -.70, .128, -.07504, .0042, .033, .00726, .054, -.0033, .159, &
    .0015/

  Data gox/.69, 0.0, 4*0.0, .46, 0.0, -.112, .0063, .050, .011, .0808, -.005, &
    .2391, .0023, -.112, .0063, .050, .011, .0808, -.005, .2391, .0023, .128, &
    -.025, -.46, .133, .138, -.023, -.70, .128, .51, -.1, -1.84, .53, .55, &
    -.091, -2.80, .51, .22, .0056, 6*0.0, 8*0., .1261, 0.0, .2869, 0.0, 0.0, &
    0.0, .28, 0.0/
  Data gf/.2246, .0075, .3055, .0062, -.0575, -.0017, .4248, -.0049, -.063, &
    0., .407, 0., -.0238, 0., .478, 0., -.0157, 0., .188, 0., .0195, 0., .283, &
    0., -.001, 0., .134, 0., .123, 0., .158, 0., .51, -.1, -1.84, .53, .55, &
    -.091, -2.8, .51, .039, .0154, 6*0., 8*0., .22, .056, 6*0., 8*0./
  Data gne/ -.72117, 1.2406, 11.746, 8.2169, -7.7772, 1.0227, .70828, 4.5400, &
    4.1450, -3.8657, -1.04290, .84411, 6.7031, 3.1927, -3.1693, -.039582, &
    .26837, .25803, .086929, -.086641, -.017323, .29514, .29301, .13223, &
    -.13071, -.071593, .15504, .12598, -.000449, .000523, .040202, .25113, &
    .14858, .030780, -.030745, .45761, .38477, .52142, .92153, -.91649, &
    -.17862, .32249, .28172, .040677, -.040639, -.062670, .14921, 1.5354, &
    1.0586, -1.0219, -.057871, .030701, .36471, .14784, -.14293, .093106, &
    -.001108, -.067652, -.021663, .021657/

  Data gna/.1586, .007375, .1866, .04156, .02403, -.02416, .3100, -.00098, &
    -.3245, 0.0, .5548, 0.0, -.1562, 0.0, .266, 0./
  Data gna4/1.153, -.0333, 0.0, .001, .0053, .058, .67, -.0133, 0.0/
  Data gmgii/.24, 37.3, 68.3, 96., 96., 2.39, 5.99, 5.99, 0., 0., 0., .142/
  Data gmg/.9, .2, .25, .23, .23, .2, .2, .2, .14, .0094/
  Data gmgp/ -.1333, .0155, -.6758, .0577, .5057, -.0323, .314, 0.0, -.294, &
    0.0, .5043, 0.0, -.1619, 0.0, .2606, 0.0/
  Data gal/0.0139, .01, .06, 0., 0.2, 3*0., 1.64, .11/
  Data gal2/0.0, .28, .28, 0., 0., 3*0., .00, .28/
  Data a/.022, .625, .660, .240, .061, .256, .0368, .271, .833, .696/
  Data b/.0035, .360, .432, .019, .364, .354, .343, .794, .404, .387/
  Data c/.1261, .1261, .1261, .477, .162, .0108, .138, .0719, .1261, .1261/
  Data gmshell/.453, .91, .13, .93, .2, .2, .7, .98, .13, .93, 0., 0., .59, &
    .95, .13, .93, 0., 0., .51, .95, .20, .29, 0., 0./
  Data gk/.35, .35, 1.1, .092, .91, .97, 1.1, 3*0./

  Data gca/.35, .21, 43., .14, .12, .43, 37., .20, .11, 4./
  Data a1/0., .4508, .3957, .5584, .5055, .5055, .5537, .5537, .5108, .5108, &
    1.235, 1.235/
  Data b1/0., -.02787, .4532, .1363, .3461, .3461, .3468, .3468, .3304, .3304, &
    -.6625, -.6625/


  Data c1/0., .1570, -.1219, .1158, -.01872, -.01872, .00516, .00516, .1007, &
    .1007, .8310, .8310/
! SHOULD PUT IN SENSIBLE F'S FOR POTASSIUM ISOSEQUENCE
  t4 = t/10000.
  st = sqrt(t)
  gaunt_rs = 0.2
  y = e*11590./t
  cc = exint1(y, 2)
  If (j==1 .And. n>2) Go To 360

  ii = n - j + 1
  If (ii>=18) Go To 330
  Go To (100, 110, 140, 150, 180, 190, 200, 210, 220, 230, 240, 280, 290, 320, &
    320, 320, 320) ii

100 Continue
! HYDROGENIC IONS: MEWE; HAYES AND SEATON: 3S,3D FROM THOMAS; OTHERS ME
  gaunt_rs = ah(l) + bh(l)*y*cc + .276*cc
  If (l==5) gaunt_rs = .212 - .203*(y-y*y*cc) + .0256*cc
  If (n>2) Return
  If (l==5) gaunt_rs = .1198*y*cc - .1194*(y-y*y*cc) + .0328*cc
  Go To 370
110 Continue
! HELIUM-LIKE IONS:  PRADHAN+CASCADES FOR N=2, MEWE N=3,4
! NEUTRAL HE FROM Berrington et al J Phys B 18,4135; Aggarwal Ap J 278,
! L = 9 HAS GAUNT_RS = 0.0,  SATELLITE LINES : L = 6 DONE BY BRANCH FROM L
  If (n==2) Go To 130
  If (l>3) Go To 120
  gaunt_rs = gncrch(ghe, l, y, cc, n)
  Return
120 If (l==4) gaunt_rs = (1.+7./n)*(.053+.022*(y*y*y*cc-y*y+y)+.276*cc)
  If (l==5) gaunt_rs = (1+1.5/n)*(.053+.022*(y*y*y*cc-y*y+y)+.276*cc)
  If (l==6) gaunt_rs = 0.
  If (l==7) gaunt_rs = .04*(y-y*y*cc)
  If (l==8) gaunt_rs = .02
  If (l==9) gaunt_rs = 0.
  Return
130 If (l==1) gaunt_rs = amin1(.00656*t4**(.89), .30*cc)
  If (l==4 .Or. l==5) gaunt_rs = .30*cc
  If (l==2) gaunt_rs = amin1(.0283*t4**(-.10), .359*t4**(-.569))
  If (l==3) gaunt_rs = amin1(.0105*t4**(.30), .095*t4**(-.479))
  If (l==6) gaunt_rs = 0.0
  If (l==7) gaunt_rs = 4.08E-7*t**.941
  If (l==8) gaunt_rs = 2.46E-8*t**1.22
  If (l==9) gaunt_rs = 0.0
  Go To 370
140 Continue
! LITHIUM-LIKE IONS: MEWE MODIFIED TO FIT CLOSE-COUPLING CALCULATIONS O
! L=1 modified Dec 17, 1991
  i = n - 4
  If (n>=10) i = n/2
  If (n>=26) i = n/2 - 2
  If (l==1) gaunt = a1(i) + (b1(i)*y-c1(i)*y*y+0.279)*cc + c1(i)*y
  x = 1./(n-3.)
! IF (L .EQ. 1) GAUNT_RS = (0.68 + .02*N) * ((.7+.35*X) +
! 1  ((1.-.8*X)*Y+(.5-.5*X)*Y*Y +.28)*CC - (.5-.5*X)*Y)
  If (l==2) gaunt_rs = .053 + .16*cc
  If (l==3 .Or. l==4) gaunt_rs = -(.16+.32*x) + ((.8-.56*x)*y+.2*y*y+.28)*cc - &
    .2*y
  If (l==5) gaunt_rs = (.19+.25*x) + .079*cc
  If (l==6) gaunt_rs = .31 - .1*y*cc
  If (l==7) gaunt_rs = .096 + .32*x
  If (l==8) gaunt_rs = .13
  Go To 370
150 Continue
! BERYLLIUM-LIKE IONS:  QUB FOR N=2 UP THROUGH NE; MANN, ROBB & SAMPSON
! MANN AND MALINOVSKY FOR N = 3 AND QUB O V N=3 FROM WIDING
! N = 4 FROM JOHNSTON & KUNZE, MASON & STOREY (FE XXII) AND LI-LIKE
  i = n - 5
  If (n>=10) i = n/2 - 1
  If (n>=26) i = n/2 - 3
  If (l/=1) Go To 160
  tl = alog10(t)
  gaunt_rs = .54 + .0125*n + .135*cc
  If (n<=10) gaunt_rs = gbe1(i) + gbe2(i)*tl
  Go To 370
160 If (l/=2) Go To 170
  gnt = gbe3(i)/(1.+gbe4(i)/y) + gbe5(i)/t
  If (n==6) gnt = .046/(1.+t*t/6.25E10)
  gaunt_rs = gnt*reduc(1)
  Go To 370
170 Continue
  If (l<=7) gaunt_rs = gncrch(gbe, l-2, y, cc, n)
  If (l==8) gaunt_rs = .1261 + .2869*y*cc + .276*cc
  If (n>8) Go To 370
  em = -4.67 + 1.86*n
  ex = exp(-em*11590./t)
  m = n - 5
  If (l==10) gaunt_rs = bespec(1, m)*(1.-pm(1))*ex + bespec(2, m)*pm(1)
  If (l==11) gaunt_rs = bespec(3, m)*(1.-pm(1))*ex + bespec(4, m)*pm(1)
  Go To 370
180 Continue
! BORON ISOSEQUENCE: N=2 INTERPOLATED FROM O IV AND FE XXII OF ROBB; GO
! NA-S, DERE ET AL AR, MANN C II;  OVEREST BY 36% NEAR THRESHOLD WRT RO
! N = 3 INTERPOLATED C II AND FE XXII FROM MANN; 2S-3L SCALED FROM BE-,
! TO 8% FOR 2S-3P.  N = 4 FROM MASON & STOREY
! L = 9 INTERCOMBINATION LINE
  If (l<=6) gaunt_rs = gncrch(gb, l, y, cc, n)
  If (l==9) gaunt_rs = reduc(2)
  If (l==12) gaunt_rs = 0.
  Go To 370
190 Continue
! CARBON ISOSEQUENCE: N=2 INTERP FROM MANN AND ROBB;  AGREES TO 5-10% W
! N=3,4 GENERAL B-F; AGREES 3-10% WITH MANN 2P-3D; L=6 INTERCOMB
  If (l<=5) gaunt_rs = gncrch(gc, l, y, cc, n)
  If (l==7) gaunt_rs = .198 + .0069*n
  If (l==8) gaunt_rs = .22 + .0056*n
  If (l==9) gaunt_rs = .1261 + (.2869*y+.28)*cc
  If (l==6) gaunt_rs = reduc(3)
  If (l==6 .And. n==6) gaunt_rs = sqrt(.0001*t)*reduc(3)
  If (l==12) gaunt_rs = 0.
  Go To 370
200 Continue
! NITROGEN SEQUENCE:   N=2 FROM MANN FE XX AND MASON & BHATIA MG,SI,S,A
! N=3 AND 4 GENERAL B-F
  If (l==1) gaunt_rs = (1.086-1.71/j)*(.3642+.9358*y*cc-.3758*(y-y*y*cc)+.3586 &
    *cc)
  If (l>=2 .And. l<=5) gaunt_rs = gncrch(gn, l-1, y, cc, n)
  If (l==12) gaunt_rs = 0.
  If (l==8) gaunt_rs = .22 + .0056*n
  If (l==9) gaunt_rs = .1261 + .2869*y*cc + .276*cc
  Go To 370
210 Continue
! OXYGEN ISOSEQUENCE:  N=2 FROM MANN, ROBB FE XIX AND FROM BHATIA,F&D S
! OTHERS GENERIC B-F
  If (l<=8) gaunt_rs = gncrch(gox, l, y, cc, n)
  If (l==10) gaunt_rs = .16 + .0015*n
  If (l==15) gaunt_rs = 0.
  Go To 370
220 Continue
! FLUORINE SEQUENCE:   N=2 FROM MANN AL V AND ROBB FE XVIII;   OTHERS G
  If (l<=8) gaunt_rs = gncrch(gf, l, y, cc, n)
  If (l==9) gaunt_rs = .1261 + (.2869*y+.28)*cc
  Go To 370
230 Continue
! NEON SEQUENCE: SMITH ET AL INCLUDING CASCADES AND RESONANCES
! ASSUME THEY ARE ALL LIKE FE XVII
  If (l<=12) gaunt_rs = gne(1, l) + (gne(2,l)+gne(3,l)*y+gne(4,l)*y*y)*exint1( &
    y, 2) + gne(5, l)*y
  Go To 370
240 Continue
! SODIUM SEQUENCE:  MANN, FLOWER & NUSSBAUMER, AND BLAHA
  If (n/=12) Go To 250
  gaunt_rs = gmgii(l)
  If (l==1) gaunt_rs = .112 + (.0269*y-.0998*y*y+.318)*cc + .0998*y
  If (l==2) gaunt_rs = 141 + (59.3*y-671*y*y+.858)*cc + 671*y
  Return
250 Continue
  If (l>2) Go To 260
  gaunt_rs = gncrch(gna, l, y, cc, n-10)
  If (n==14 .And. l==2) gaunt_rs = -.0172 + (.832*y+.029*y*y+.3513)*cc - &
    .029*y
  Return
260 Continue
  If (l>5) Go To 270
  ll = l - 2
  gaunt_rs = gna4(1, ll) + gna4(2, ll)*n + gna4(3, ll)*cc
  Return
270 If (l==6) gaunt_rs = -.16 + .8*y*cc - .2*(y-y*y*cc) + .276*cc
  If (l==7) gaunt_rs = .44 - 0.1*y*cc
  If (l==8) gaunt_rs = .15 - .05*y*cc
  If (l==9 .Or. l==10) gaunt_rs = 0.03
  If (l==11) gaunt_rs = 0.07
  If (l==12) gaunt_rs = 0.15
  Return
280 Continue
! MAGNESIUM SEQUENCE:   INTERPOLATE SI TO FE FOR 3P, 4P EXCITATIONS;
! USE MANN FE GAUNT_RS FACTORS FOR 3D, 4S, 4F, AND INTERCOMB.  ASSUME
! G FOR 4D = G FOR 4F;  L=10 IS INTERCOMBINATION LINE
  gaunt_rs = gmg(l)
  If (l<=2) gaunt_rs = gncrch(gmgp, l, y, cc, n)
  If (l==10) gaunt_rs = reduc(4)
  If (l==10 .And. n==14) gaunt_rs = reduc(4)*(20000./t)**(.406)
  Return
290 Continue
! ALUMINUM ISO-SEQUENCE:   SI II FROM ROBERTS, S IV 3P FROM MANN, 3D
! FE XIV FROM BLAHA USED FOR NI AND FOR ALL N=4
! RESONANCES INCLUDED FOR METASTABLE USING BHATIA COLLISION STRENGTH
! SI II 1814 LINES FROM BROWN, FERRAZ AND JORDAN
! FE XIV RESONANCES COMPUTED AS IN SMITH ET AL FROM BLAHA OMEGAS
  If (l>=4 .And. l<=8) Go To 310
  If (n>14) Go To 300
  gaunt_rs = gal(l) + gal2(l)*cc
  If (l==1) gaunt_rs = gal(1)*1.9E7*st/(dene+st*1.9E7)
  Return
300 If (n>20) Go To 310
  cn = st*2.7*10.**(7.+j/2.)
  If (l==1) gaunt_rs = ((.128+(.3449*y+.3544*y*y)*cc-.3544*y)*1.8/(1.+.25/y))* &
    .069*cn/(dene+cn)
  If (l==2) gaunt_rs = .0405 + (.2052*y+.0328*y*y+.2311)*cc - .0328*y
  If (l==3) gaunt_rs = .142 + .147*cc
  If (l==9) gaunt_rs = .1330 + (.3833*y+.0934*y*y+.1611)*cc - .0934*y
  If (l==10) gaunt_rs = .1684 + (.2432*y+.0484*y*y+.2638)*cc - .0484*y
  Return
310 Continue
  gaunt_rs = a(l) + b(l)*(1.-1./(1.+c(l)/y))
  If (l==1) gaunt_rs = (.0121+.0036/(1.+.1261/y))*6.0E17/(dene+6.0E17)
  Return

! SILICON THROUGH CHLORINE SEQUENCES
320 gaunt_rs = gmshell(l, ii-13)
  If (l==3 .And. ii/=17) gaunt_rs = gmshell(l, ii-13) - .114/(1.+.138/y)
  Return
330 Continue
  If (ii>18) Go To 340
  Return
340 If (ii>19) Go To 350
  gaunt_rs = gk(l)
  If (n==20) gaunt_rs = gca(l)
  Return
350 If (l<=3) gaunt_rs = amax1(.2, -.606*alog10(y)-.052)
  Return
360 Continue
  gaunt_rs = .01
  If (y<=10.) gaunt_rs = amin1(exp(-.7*alog(y)-2.3), .28*cc)
370 Continue
  Return

End Function
Function gncrch(r, l, y, cc, n)
  Dimension :: r(1)

  il = 8*(l-1)
  a = r(il+1) + r(il+2)*n
  b = r(il+3) + r(il+4)*n
  c = r(il+5) + r(il+6)*n
  d = r(il+7) + r(il+8)*n
  gncrch = a + b*y*cc + c*(y-y*y*cc) + d*cc
  Return

End Function
Function delt(n, j, l)
! OUR VERY OWN MODIFICATIONS TO JACOBS
! DENSITY DEPENDENCE EXTERNAL
! TAKE 0.5 FOR MERTS ET AL   DELTA N .NE. 0 EXCEPT S
! ASSUME 3D OF N,O,F INTERP BETWEEN C AND NE ISOSEQUENCES
! THIS USES YOUNGER'S CLAIM THAT DELT=1 FOR HE-LIKE RESONANCE LINE;
! TRY BELY-DUBAU??????
  delt = 0.
  If (j==1) Return
  i = n - j + 1
  If (l==12 .And. i/=10) Go To 100
  If (l/=1 .Or. i<=2) Go To 110
  If (i==13) Go To 110
100 delt = 1.
  Go To 280
110 Continue
  If (i>=19) Go To 270
  If (i==18) Go To 260
  If (i>=15) Go To 250
  Go To (120, 130, 140, 150, 160, 170, 180, 190, 200, 210, 220, 230, 240, 250) &
    i
120 Continue
! BURGESS & TWORKOWSKI ARE IN ALPHADI
  If (l==1) delt = 1.
  If (l>1) delt = 0.1
  If (l==5) delt = 0.
  Go To 280
130 If (l==4 .Or. l==5) delt = 0.1
  If (l==9) delt = 1.0
  Go To 280
140 If (l>=2 .And. l<=4) delt = .5*(-.35+1.26*n/(n+11.))
  Go To 280
150 If (l==7) delt = amax1(0.5*(-.55+1.50*n/(n+11.)), 0.05)
  Go To 280
160 If (l<=6) delt = 1.
  If (l==4 .Or. l==5) delt = 0.5*(-.54+2.2*n/(n+11.))
  Go To 280
170 If (l<=5) delt = 1.
  If (l==4) delt = amax1(0.5*(-1.28+3.2*n/(n+11)), 0.05)
  Go To 280
180 If (l==4) delt = 1.
  If (l==5) delt = amax1(0.5*(-1.44+3.4*n/(n+11.)), 0.05)
  If (l==7) delt = amax1(0.25*(-1.44+3.4*n/(n+11.)), 0.025)
  Go To 280
190 If (l<=3) delt = 0.5*(-1.60+3.6*n/(n+11.))
  If (l==4 .Or. l==5) delt = 1.
  If (l==7) delt = .25*(-1.6+3.6*n/(n+11.))
  If (l==15) delt = 1.
  Go To 280
200 If (l<=3) delt = 0.5*(-1.75+3.8*n/(n+11.))
  If (l==5) delt = 1.
  If (l==6) delt = 0.25*(-1.75+3.8*n/(n+11.))
  Go To 280
210 If (l<=2) delt = 1.
  If (l>=3 .And. l<=5) delt = 0.5*(-1.9+4.0*n/(n+11.))
  If (l>=6 .And. l<8) delt = .25*(-1.9+4.*n/(n+11.))
  If (l==14) delt = 1.
  Go To 280
220 Continue
  If (l==2 .Or. l==3 .Or. l==6) delt = .25
  Go To 280
230 If (l==2) delt = .25
  Go To 280
240 If (l==2 .Or. l>=9) delt = 1.
  If (l==3 .Or. l==5) delt = .5
  If (l==7) delt = 0.25
  Go To 280
250 If (l<=2) delt = 1.
  If (l==3 .Or. l==4) delt = .25
  Go To 280
260 delt = 0.25
  If (l<=3) delt = 1.
  If (n/=20) Return
  delt = 0.
  If (l==6) delt = 0.05
  Go To 280
270 delt = 0.25
280 Continue
  Return
End Function


Function phot(n, j, e, t, dene)
! USES REILMAN & MANSON FOR INNER SUBSHELL PHOTOIONIZATION:  2S FOR B,C
! SEQUENCES,; 3S FOR AL,SI,P SEQUENCES AND 3P FOR K,CA,SC SEQUENCES
! ASSUMES THAT OUTER SHELL IONIZATION IS IN FOUR TERM POLYNOMIAL,
! S2*X**2 + S3*X**3 + S4*X**4 + S5*X**5
  Common /com/cnc(12, 30), ptot(12, 220), abin(1000), bin(1000)
  Common /hot/he, hea(2)
  Common /hetot/heating
  Common /dat/ee(30), ea(30), s2(30), wave(220), e3(220), f(220), ll(30), &
    s3(30), s4(30), s5(30)
  Common /params/nj(12), abunj(12), abund, binmin, binsyz, nbin
  Common /pt/rf(500), tau(150), tauhe(150), tplus(150)
  Common /result/conce(30), gndrec(30), power(220), rhy, heneut, heplus, dne, &
    pcool, pou, pot, re, tu, pm(4)
  Dimension :: inner(30), ainner(12, 9), binner(12, 9), einner(12, 9)
  Data inner/4*0, 1, 2, 3, 5*0, 4, 5, 6, 3*0, 7, 8, 9, 9*0/
  Data einner/0., 30.9, 55.8, 87.6, 172., 283., 422., 589., 784., 1006., &
    1842., 2178., 0., 15.6, 36.7, 63.8, 139., 241., 371., 528., 713., 925., &
    1731., 2056., 2*0., 20.3, 42.6, 108., 201., 321., 469., 644., 847., 1622., &
    1938., 6*0., 22.9, 57.6, 105.2, 165., 421., 531., 6*0., 13.5, 43.8, 87.6, &
    144., 388., 493., 7*0., 30.7, 70.4, 123., 356., 458., 9*0., 28., 213., &
    296., 9*0., 38., 190., 271., 10*0., 169., 246./
  Data ainner/0., 2.357, 2.745, 1.443, .6024, .3391, .2058, .1385, .09804, &
    .0764, .04443, .0376, 0., 18.26, 6.950, 2.817, .9583, .5001, .2672, .1685, &
    .1184, .0913, .04837, .0407, 2*0., 27.80, 8.139, 1.630, .7415, .3621, &
    .2394, .1478, .112, .05442, .0455, 6*0., 4.961, 2.695, 1.401, .8416, &
    .2929, .2254, 6*0., 10.489, 5.507, 1.575, 1.409, .4378, .3263, 7*0., &
    12.06, 4.726, 2.351, .6321, .4799, 9*0., 29.6, 2.452, 1.788, 9*0., 53.8, &
    3.265, 2.070, 10*0., 4.078, 2.351/
  Data binner/0., -.7352, .0201, .4321, .4094, .2605, .1935, .1423, .111, &
    .0865, .04050, .0343, 0., -6.256, -.2887, 1.463, 1.029, .5726, .4080, &
    .2885, .2132, .164, .07817, .0654, 2*0., -14.93, 1.194, 2.298, 1.106, &
    .7338, .4742, .3510, .267, .1244, .104, 6*0., -4.607, -1.863, -.6435, &
    -.2594, -.0225, -.0111, 6*0., -10.412, -4.613, -1.118, -.4648, -.0221, &
    .0105, 7*0., -11.09, -3.357, -1.014, -.0255, -.0078, 9*0., -28.3, -.386, &
    -.360, 9*0., -45.8, -.610, -.178, 10*0., -.834, .0037/

  phot = 0.
  he = 0.
  Do noj = 1, 12
    no = noj
    If (n<=nj(noj)) Go To 100
  End Do
100 Continue
  If (j/=1) Go To 110
  If (n==12) phot = 1.0E-10
  If (n==16) phot = 2.0E-10
  If (n==26) phot = 3.0E-10
  If (n/=6) Go To 110
! PHOTOIONIZATION OF C I FROM 1D LEVEL
  t4 = amin1(.0001*t, 3.)
  om = 1.16*t4*(1.085-.07507*t4-.02150*t4*t4)
  cn = sqrt(t)*5.*3.26E-4/(8.63E-6*om)
  fd1 = 0.4*exp(-1.45/t4)*dene/(dene+cn)
  ily = (10.2-binmin)/binsyz + 1
  phot = 2.5E-10
  If (ily>=1) phot = phot + fd1*1.03E-17*6.12E10*abin(ily)
110 Continue
  imin = (e-binmin)/binsyz + 1
  If (e<binmin+(imin+0.5)*binsyz) imin = imin + 1
  If (imin>=nbin) Return
  imin = max0(imin, 1)
  If (n/=1) Go To 120
  ss2 = 0.00
  ss3 = 8.44
  ss4 = -2.14
  ss5 = 0.00
  fact = 1.6E-12*1.E23/(1.+rhy)
  Go To 130
120 Continue
  ss2 = s2(j)
  ss3 = s3(j)
  ss4 = s4(j)
  ss5 = s5(j)
  fact = 1.6E-12*1.E23*10.**(abund-12.)*conce(j)
130 Continue
  iso = n - j + 1
  ia = (ea(j)-binmin)/binsyz
  ia = max0(ia, 1)
  imax = nbin
  If (iso>=3) imax = min0(nbin, ia)
  iiso = inner(iso)
  If (iiso==0) Go To 140
  ein = einner(no, iiso)
  iin = (ein-binmin)/binsyz + 1
  If (ein<binmin+(iin+0.5)*binsyz) iin = iin + 1
  iin = max0(iin, 1)
  imax = min0(iin-1, nbin)
140 Continue
  hnu = binmin + (imin-0.5)*binsyz
  Do i = imin, imax
    hnu = hnu + binsyz
    x = e/hnu
    x2 = x*x
    sg = ss2*x2 + ss3*x2*x + ss4*x2*x2 + ss5*x2*x2*x
    phot1 = sg*(1.E-18/1.6E-12)*rf(i)*abin(i)/hnu
    phot = phot + phot1
    he = he + phot1*hnu*fact
  End Do
  If (n==2) hea(j) = he
! INNER SUBSHELL
  If (iiso==0) Go To 150
  If (iin>=nbin) Go To 150
  imax = min0(ia, nbin)
  ss2 = ainner(no, iiso)
  ss3 = binner(no, iiso)
  Do i = iin, imax
    hnu = binmin + binsyz*(i-0.5)
    x = ein/hnu
    sg = ss2*x*x + ss3*x*x*x
    phot1 = sg*(1.E-18/1.6E-12)*rf(i)*abin(i)/hnu
    phot = phot + phot1
    he = he + phot1*hnu*fact
  End Do
150 If (n>=6) heating = heating + he
  Return

End Function

Subroutine aphot(n, dene, icont)
! AUGER PHOTOIONIZATION : N=1 AND N=2 SHELLS
! USES COX & DALTABUIT AND REILMAN & MANSON CROSS SECTIONS
  Common /com/cnc(12, 30), ptot(12, 220), abin(1000), bin(1000)
  Common /hot/he, hea(2)
  Common /hetot/heating
  Common /dat/e(30), ea(30), s2(30), wave(220), e3(220), f(220), ll(30), &
    s3(30), s4(30), s5(30)
  Common /params/nj(12), abunj(12), abund, binmin, binsyz, nbin
  Common /pt/rf(500), tau(150), tauhe(150), tplus(150)
  Common /result/conce(30), ca(30), power(220), rhy, heneut, heplus, dne, &
    pcool, pou, pot, re, tu, pm(4)
  Common /rates/c1(30), c2(30), cth
  Dimension :: a2p(20, 20), b2p(20, 20)
  Data a2p/20*0, 12.60, 20.58, 38*0., 4.520, 3.945, 5.546, 7.213, 36*0., &
    1.286, 1.316, 1.843, 2.556, 3.363, 3.654, 34*0., .6109, .6256, .8755, &
    1.111, 1.254, 1.578, 1.919, 2.307, 32*0., .4075, .4438, .4801, .5164, &
    .6798, .8431, .9660, 1.089, 1.288, 1.488, 110*0., .1339, .1531, .1722, &
    .1909, .2095, .2320, .2545, .3022, .3499, .3628, .3756, .4630, .5504, &
    .5820, .6136, .6679, 24*0., .1145, .1294, .1409, .1525, .1644, .1763, &
    .1797, .1831, .2005, .2179, .2425, .2670, .3393, .4116, .4266, .4415, &
    .5339, .6263, 42*0./
  Data b2p/20*0., .1887, -7.098, 38*0., 1.349, 5.109, 4.291, 6.750, 36*0., &
    2.408, 2.546, 3.116, 3.060, 2.981, 5.273, 34*0., 1.377, 1.456, 1.782, &
    1.930, 2.384, 2.572, 3.023, 3.451, 32*0., .7788, .9714, 1.164, 1.356, &
    1.488, 1.620, 1.869, 2.118, 2.261, 2.405, 110*0., .3454, .3781, .4107, &
    .4565, .5023, .5590, .6157, .6283, .6409, .7392, .8375, .8258, .8142, &
    .9381, 1.062, 1.207, 24*0., .2541, .2787, .3087, .3387, .3770, .4152, &
    .4907, .5661, .6055, .6499, .6729, .7010, .6410, .5809, .6970, .8130, &
    .7768, .7406, 42*0./

  he = 0.
  Do j = 1, 30
    ca(j) = 0.
  End Do

  rn = n
  nn = n - 2
  If (n<=2) Return
  Do j = 1, nn
    fact = 1.6E-12*1.E23*10.**(abund-12.)*conce(j)
    iso = n - j + 1
    riso = iso
    pa1 = 0.
    pa2 = 0.
    eks = 13.6*rn*rn*riso**(.24-.43/alog10(rn))
    ik = (eks-binmin)/binsyz
    imax = nbin
    If (iso>=11) imax = min0(ik, nbin)
    imin = (ea(j)-binmin)/binsyz + 1
    If (ea(j)<binmin+(imin+0.5)*binsyz) imin = imin + 1
    imin = max0(1, imin)
    If (imin>=nbin) Go To 140

    If (iso<=10) Go To 100
    a2 = a2p(iso-10, n-10)
    b2 = b2p(iso-10, n-10)
    Do i = imin, imax
      hnu = binmin + binsyz*(i-0.5)
      x = ea(j)/hnu
      sg = a2*x*x + b2*x*x*x
      phot1 = sg*(1.0E-18/1.6E-12)*rf(i)*abin(i)/hnu
      pa1 = pa1 + phot1
      he = he + phot1*hnu*fact
    End Do
100 Continue
    imin = max0(ik+1, 1)
    If (imin>=nbin) Go To 110
    sss = 1.66*eks**(0.071)
    Do i = imin, nbin
      hnu = binmin + binsyz*(i-0.5)
      x = eks/hnu
      sg = (295./eks)*x**sss
      phot1 = sg*(1.0E-18/1.6E-12)*rf(i)*abin(i)/hnu
      pa2 = pa2 + phot1
      he = he + phot1*hnu*fact
    End Do
110 Continue
    fluor = 0.5*(rn/26.)**4
    ca(j) = (pa1+pa2)*(1.-fluor)/dene
    c1(j) = c1(j) + (pa1+pa2)*fluor/dene
    If (iso/=3) Go To 120
    c1(j) = c1(j) + pa2/dene
    ca(j) = 0.
120 If (iso/=11) Go To 130
    ca(j) = (1.-fluor)*pa2/dene
    c1(j) = c1(j) + (pa2*fluor+pa1)/dene
130 Continue
140 End Do

  heating = heating + he

  If (icont/=0) Return

! AVRETT'S METHOD FOR INCLUDING AUTOIONIZATION IN EQUILIBRIUM

  q = 0.
  x = ca(1)
  nn = n - 1
  Do j = 1, nn
    c1(j) = c1(j) + x
    If (c1(j)<=0.) Go To 150
    q = c2(j+1)*ca(j)/c1(j)
    x = ca(j+1) + q
  End Do
150 Continue
  Return

End Subroutine
Function secont(b, q)
  External :: gamma

  a = b/(q*q)
  r = q**.6666667
  If ((a-50.0)<=0.) Then
    Go To 110
  Else
    Go To 100
  End If

100 secont = (1.0-.1728/r) - 0.0496/(r*r)
  Go To 120
110 Continue
  secont = 1.0 + (.1728/r)*((.33333-2.*a)*gamma(.33333,a)+1.0) - &
    (.0496/(r*r))*(.66667*(1.-a-3.*a*a)*gamma(.66667,a)+(1.+2.*a))
120 Continue
  Return

End Function
Subroutine lnprt(num)
  Common /hln/ext(11), half, halfc, alfly, alflyc, hbet, hbetc, twop
  Common /params/nj(12), abunj(12), abund, binmin, binsyz, nbin
  Common /fln/fwv(320), ftot(320), fout(320), fcool
  Common /fel/wj(12, 220), fj(12, 220), e3j(12, 220), pwr(12, 220), &
    c1j(12, 30), c2j(12, 30), ej(12, 30), sj(12, 30), cj(12, 30), llj(12, 30), &
    alfj(12, 30), esj(12, 30), sigj(12, 30)
  Common /com/cnc(12, 30), ptot(12, 220), abin(1000), bin(1000)
  Integer :: bb
  Common /lnpr/jmx(11), nf(11), ldum(5), jdum(5), wdum(5), pdum(5)

  Write (7, 160) hbet, hbetc, twop
  c = 100./(hbet+hbetc)
  Do no = 1, num
    n = nj(no)
    Write (7, 170) n
    loc = 1
    ix = 0
    jt = amin0(jmx(no), n)
    Do j = 1, jt
      iy = llj(no, j)*3
      If (iy<=0) Then
        Go To 130
      Else
        Go To 100
      End If

100   Continue
      Do l = 1, iy
        bb = ix + l
        If (wj(no,bb)<=0.) Then
          Go To 120
        Else
          Go To 110
        End If

110     Continue
        ldum(loc) = l
        jdum(loc) = j
        wdum(loc) = wj(no, bb)
        pdum(loc) = ptot(no, bb)*c
        If (loc==5) Write (7, 180) n, (jdum(k), ldum(k), wdum(k), pdum(k), k=1 &
          , 5)
        loc = mod(loc, 5) + 1
120   End Do
130   ix = ix + iy
    End Do
    il = loc - 1
    If (loc/=1) Write (7, 180) n, (jdum(k), ldum(k), wdum(k), pdum(k), k=1, il &
      )
    If (n==2) Go To 150
    Write (7, 190) n
    nlow = nf(no-1) + 1
    nhigh = nf(no)
    loc = 1
    Do i = nlow, nhigh
      If (fwv(i)*ftot(i)<=0.) Go To 140
      ldum(loc) = i
      wdum(loc) = fwv(i)
      pdum(loc) = ftot(i)*c
      If (loc==5) Write (7, 200) n, (ldum(k), wdum(k), pdum(k), k=1, 5)
      loc = mod(loc, 5) + 1
140 End Do
    il = loc - 1
    If (loc/=1) Write (7, 200) n, (ldum(k), wdum(k), pdum(k), k=1, il)
150 End Do
  Return
160 Format ('    HBET, HBETC, TWOP    ', 3E10.3)
170 Format (' LINES FOR ELEMENT', I5)
180 Format (I4, 5(I4,I3,F8.2,E10.3))
190 Format ('  ***** FORBIDDEN LINES ***** N=', I3)
200 Format (1X, I2, 5(I6,F9.1,E10.3))

End Subroutine

Subroutine radrd(icont)
! Read in radiation profile.

  Character (100) :: label
  Common /com/cnc(12, 30), power(12, 220), abin(1000), bin(1000)
  Common /params/nj(12), abunj(12), abund, binmin, binsyz, nbin

  If (icont==0) Return
  ird = nbin/10
  Read (3, 100) label
  Do i = 1, ird
    il = 10*i - 10
    Read (3, 110)(abin(il+k), k=1, 10)
  End Do
  Return
100 Format (A100)

110 Format (13X, 10E10.2)

End Subroutine
Function solvx(f, e, h)
  Common /slv/denz, press, pb, dynam
! CALCULATE COMPRESSION
  ax2 = 2.*denz*(h-e)/pb
  ax1 = -5.*press/pb
  ax0 = 4.*dynam/pb
  qx = ax2*ax2/9. - ax1/3.
  rx = (3.*ax0-ax1*ax2)/6. + (ax2**3)/27.
  plink = rx/qx**1.5
  If (plink>0.997) Go To 100
  thetax = acos(plink)/3.
  solvx = sqrt(qx)*(cos(thetax)+1.732*sin(thetax)) - ax2/3.
  Return
100 x = (-ax1+sqrt(ax1*ax1-4.*ax2*ax0))/(2.*ax2)
  Do k = 1, 4
    x = (-ax1+sqrt(ax1*ax1-4.*ax2*(ax0+x**3)))/(2.*ax2)
  End Do

  solvx = x
  Return

End Function
Subroutine hline(rhy, a, t, st, beta, alter)
  Common /hln/ext(11), half, halfc, alfly, alflyc, hbet, hbetc, twop

  rh = rhy/(1.+rhy)
  y = .75*157890./t
  cc = alog((y+1.)/y) - 0.4/(y+1.)**2
! CONTINUUM TO N .GE. 3
  ext(1) = 3.46*rh/st
! CONTINUUM TO N .EQ. 2
  ext(2) = (45.0/st)*secont(beta, 2.)*rh/8.
! CONTINUUM TO GROUND
  ext(3) = 45.0*secont(beta, 1.)*rh/st
  If (alter>=1.) ext(3) = 0.
! TOTAL RECOMBINATION LINES
  chi = (.735+alog(beta)+.3333/beta)/2.
  If (beta<=.5) chi = 1.2*beta + (.55+1.04*alog(beta))*beta*beta + &
    (1.01*alog(beta)-0.43)*beta**3
  chi = chi - p(beta, 1.)
  ex = .4288 + alog(beta)/2. + 0.469*beta**(-.333) - seaton(beta, 1.)
  ext(4) = (2.06*1.607*13.6/st)*(ex+chi/beta-secont(beta,2.)/8.-.077)*rh
! HBETA FROM RECOMBINATION
  ext(5) = .01275*exp(-410./t)*rh/(.0001*t)**.937
! LYMAN ALPHA FROM RECOMBINATION
  ext(6) = 2.06*1.6027*13.6*.75*ex*rh/st
! TOTAL FROM COLLISIONAL EXCITATION
! STILL USING OLD LYMAN ALPHA EXCITATION
  sum = .416*exp(-.75*beta)*1.195 + .07912*exp(-.889*beta)*1.076 + &
    .02899*exp(-.9375*beta)*1.041 + .01394*exp(-.96*beta)*1.026 + &
    .007799*exp(-.972*beta)*1.018
  ext(7) = 4.5E-6*1.6027E-12*(1.E23/sqrt(13.6))*sum/((1.+rhy)*beta**.12)
! HBETA FROM COLLISIONAL EXCITATION : OLD H BETA
  sum = (.02899*.739*3/16.)*exp(-.9375*beta)*1.110 + &
    .0745*.01394*exp(-.96*beta)*1.068 + .007799*.079*exp(-.972*beta)*1.047
  ext(8) = (4.5E-6*1.6027E11/sqrt(13.6))*sum/((1+rhy)*beta**.12)
! LYMAN ALPHA FROM COLLISIONAL EXCITATION
  sum = .0791*exp(-.889*beta)*1.210 + .02899*exp(-.9375*beta)*1.110 + &
    .01394*exp(-.96*beta)*1.068 + .007799*exp(-.972*beta)*1.047
  ext(9) = (4.5E-6*1.6027E-12*1.E23/sqrt(13.6))*sum/(beta**.12*(1.+rhy))
! AGGARWAL FITS TO MCDOWELL AND CALLAWAY WITH OLD CASCADES
  tp = amin1(t, 5.0E5)
  om = .5*(.335+1.45E-5*tp+1.39E-10*tp*tp-5.66E-15*tp*tp*tp)
  If (tp>=25000.) om = .5*(.328+1.43E-5*tp-6.55E-12*tp*tp-2.69E-18*tp*tp*tp)
  If (t>=5.0E5) om = 8.05*.276*cc
  ext(9) = ext(9) + 8.63E-6*om*10.2*1.6E11*exp(-118000./t)/(sqrt(t)*(1.+rhy))
! H ALPHA EXCITATION; AGGARWAL ; CASE B, NO CASCADES
  om = (.175-1.31E-7*tp-4.08E-11*tp*tp+3.10E-15*tp*tp*tp)
  If (tp>=25000.) om = (.138+2.50E-6*tp-4.43E-12*tp*tp+3.59E-18*tp*tp*tp)
  ext(10) = 8.63E-6*om*1.89*1.6E11*exp(-140300./t)/(sqrt(t)*(1.+rhy))
! TWO PHOTON FROM AGGARWAL;  NO CASCADES INCLUDED
  om = .5*(.195+1.28E-5*tp-4.69E-10*tp*tp+6.19E-15*tp*tp*tp)
  If (tp>=25000.) om = .5*(.308+3.58E-7*tp+6.15E-13*tp*tp-1.08E-18*tp*tp*tp)
  ext(11) = 8.63E-6*om*10.2*1.6E11*exp(-118000./t)/(sqrt(t)*(1.+rhy))
  half = half + ext(5)*a*(2.654+1730./t)
  halfc = halfc + ext(10)*a
  hbet = hbet + ext(5)*a
  hbetc = hbetc + ext(8)*a
  alfly = alfly + ext(6)*a
  alflyc = alflyc + ext(9)*a
  twop = twop + ext(11)*a
  Return


End Subroutine

Function hcool(t, rhy, alter)
  Common /hln/ext(11), half, halfc, alfly, alflyc, hbet, hbetc, twop

  beta = 157890./t
  st = sqrt(t)
  Call hline(rhy, 0., t, st, beta, alter)
  hcool = ext(1) + ext(2) + ext(3) + ext(4) + ext(8) + ext(9) + ext(10) + &
    ext(11)
  Return

End Function
Function taup(w, it)
  Common /pt/rf(500), tau(150), tauhe(150), tplus(150)

  taup = 1.E-8
  If (w>912.) Return
  taup = taup + tau(it)*(w/912.)**3
  If (w<504.6) taup = taup + tauhe(it)*(.763*(w/504.6)**1.99+.237*(w/504.6)** &
    2.99)
  If (w<228.) taup = taup + tplus(it)*(w/228.)**3
  Return

End Function
Function expf(x)

  expf = 1.0E-37
  If (x>=-75. .And. x<=75.) expf = exp(x)
  If (x>75.) expf = 1.0E37
  Return

End Function
Function gamma(a, x)
  Double Precision :: can1, can2, cbn1, cbn2, an, bn, can, cbn
! EXP(-X) * X ** A   TAKEN OUTSIDE
  can2 = 1.
  can1 = 0.
  cbn2 = 0.
  cbn1 = 1.
  an = 1.
  bn = x
  n = 1
  can = bn*can1 + an*can2
  cbn = bn*cbn1 + an*cbn2
  fn1 = can/cbn
  Do n = 2, 1000
    If (can<=1.E30) Go To 100
    can1 = can1*1.E-30
    can = can*1.E-30
    cbn1 = cbn1*1.E-30
    cbn = cbn*1.E-30
100 Continue
    in = n
    can2 = can1
    cbn2 = cbn1
    can1 = can
    cbn1 = cbn
    mn = mod(n, 2)
    bn = mn*x + (1.-mn)
    nn = n/2
    an = nn - (1.-mn)*a
    can = bn*can1 + an*can2
    cbn = bn*cbn1 + an*cbn2
    fn = can/cbn
    If (abs((fn-fn1)/fn1)<=.0001 .And. n>=20) Go To 110
    fn2 = fn1
    fn1 = fn
  End Do

  gamma = (fn+fn2)/2.
  Return
110 gamma = fn
  Return

End Function
Function fable(x, y, z)

  If ((x-y)==0.) Then
    Go To 110
  Else
    Go To 100
  End If

100 duck = (expf(-x*z)-expf(-y*z))/(y-x)
  Go To 120
110 duck = z*exp(-x*z)
120 fable = duck
  Continue
  Return

End Function



Subroutine blnd(num)
! COMPUTES EMISSION SPECTRUM IN WAVELENGTH BINS : FIRST BLN EXTENDS FRO
! BLNMIN TO BLNMIN+BLNSYZ
! RESOLUTION AND RANGE OF CONTINUUM ARE ONLY AS GOOD AS THE
! RESOLUTION AND RANGE OF THE ENERGY BIN COMPUTATION.  BE SURE THAT
! NBIN, BINMIN AND BINSYZ PROVIDE THE NECESSARY INTERVAL.  BLN'S
! OUTSIDE THE RANGE OF THE ENERGY BINS WILL HAVE NO CONTINUUM
! CONTRIBUTION, BUT PLACEMENT OF LINES WILL NOT BE AFFECTED.
! SEPARATES RESONANCE DOUBLETS OF LI-LIKE AND NA-LIKE LINES
! FEB. 1983
  Common /bln/bln(1000), blnmin, blnsyz, nbln
  Common /params/nj(12), abunj(12), abund, binmin, binsyz, nbin
  Common /contin/brmev(1000), recev(1000), tufev(1000)
  Common /fel/wj(12, 220), fj(12, 220), e3j(12, 220), pwr(12, 220), &
    c1j(12, 30), c2j(12, 30), ej(12, 30), sj(12, 30), cj(12, 30), llj(12, 30), &
    sigj(12, 30), alfj(12, 30), esj(12, 30)
  Dimension :: l1(12), l2(12), wli(2, 12), wna(2, 12)
  Data l1/0, 28, 37, 43, 52, 85, 103, 127, 88, 148, 184, 166/
  Data l2/5*0., 13, 28, 43, 34, 64, 85, 88/
  Data wli/0., 0., 1548.2, 1550.8, 1238.8, 1242.8, 1031.9, 1037.6, 770.4, &
    780.3, 609.8, 624.9, 499.4, 520.7, 417.6, 445.8, 353.9, 389.1, 302.2, &
    344.8, 192.03, 255.11, 165.42, 234.20/
  Data wna/10*0., 2796.4, 2803.5, 1393.8, 1402.8, 933.4, 944.5, 700.4, 714.0, &
    557.7, 574.0, 335.4, 360.8, 292.0, 320.6/

! CONTINUUM

  Do ibln = 1, nbln
    wbln = (ibln-0.5)*blnsyz + blnmin
    ephot = 12399./wbln
    ibn = (ephot-binmin)/binsyz + 1
    de = ephot - 12399./(wbln+blnsyz)
    If (ibn>=1 .And. ibn<=nbin) bln(ibln) = (brmev(ibn)+recev(ibn)+tufev(ibn)) &
      *de/binsyz
  End Do

! LINES

  Do no = 1, num
    Do l = 1, 220
      If (l==l1(no)) Go To 100
      If (l==l2(no)) Go To 110
      Go To 120
100   ibln = (wli(1,no)-blnmin)/blnsyz + 1
      If (ibln>=1 .And. ibln<=nbln) bln(ibln) = bln(ibln) + &
        0.6666666*pwr(no, l)
      ibln = (wli(2,no)-blnmin)/blnsyz + 1
      If (ibln>=1 .And. ibln<=nbln) bln(ibln) = bln(ibln) + &
        0.333333*pwr(no, l)
      Go To 130
110   ibln = (wna(1,no)-blnmin)/blnsyz + 1
      If (ibln>=1 .And. ibln<=nbln) bln(ibln) = bln(ibln) + &
        0.666666*pwr(no, l)
      ibln = (wna(2,no)-blnmin)/blnsyz + 1
      If (ibln>=1 .And. ibln<=nbln) bln(ibln) = bln(ibln) + 0.33333*pwr(no, l)
      Go To 130
120   Continue
      ibln = (wj(no,l)-blnmin)/blnsyz + 1
      If (ibln>=1 .And. ibln<=nbln) bln(ibln) = bln(ibln) + pwr(no, l)
130   Continue
    End Do
  End Do
  Return
End Subroutine

Function effnew(i, t, z, n)

  t3 = .001*t/z**2
  xx = .4342*log(t3)
  If ((t3-1.)<=0.) Then
    Go To 100
  Else
    Go To 110
  End If

100 f1 = .266
  f2 = .13
  f3 = .13
  Go To 130
110 If ((t3-1.E5)>=0.) Go To 120
  f1 = .266 + .1068*xx - .074*sin(1.2566*xx)
  f2 = .130 + .1160*xx - .074*sin(1.2566*xx)
  f3 = .130 - .012*xx + .05*exp(-(xx-2.)*(xx-2.))
  Go To 130
120 f1 = .80
  f2 = .71
  f3 = .07
130 Continue
  If (i==2 .Or. i==3) effnew = 8.*(1.-f1)
  If (i==10 .Or. i==11) effnew = 18.*(1.-f2)
  If (i>=12 .And. i<=17) effnew = 18.*(1.-3.*f3-f2)
  If (i>=18) effnew = (28.-n+z)*1.8*(1.-3.*f3-f2)
  Return

End Function
Function effn(i, zeff, t)

  t3 = t/(zeff*zeff*1000.0)
  xx = 0.4342*log(t3)
  If ((t3-1.0)<=0.) Then
    Go To 100
  Else
    Go To 110
  End If

100 f1 = 0.266
  f2 = 0.13
  f3 = 0.13
  Go To 140
110 If ((t3-10.**5)<0.) Then
    Go To 120
  Else
    Go To 130
  End If

120 f1 = 0.266 + 0.1068*xx - 0.074*sin(1.2566*xx)
  f2 = 0.130 + 0.1160*xx - 0.074*sin(1.2566*xx)
  f3 = 0.130 - 0.0120*xx + 0.050*exp(-(xx-2.)*(xx-2.))
  Go To 140
130 f1 = 0.80
  f2 = 0.71
  f3 = 0.07
140 Continue
  If (i==0) Go To 160
  If ((i-18)<=0) Then
    Go To 150
  Else
    Go To 330
  End If

150 Go To (160, 170, 180, 190, 200, 210, 220, 230, 240, 250, 260, 270, 280, &
    290, 300, 310, 320, 330) i
160 eye = i
  effn = 2. - eye
  Go To 350
170 effn = 8.0
  Go To 350
180 effn = 8. - (4.*f1)
  Go To 350
190 effn = 8.*(1.-f1)
  Go To 350
200 effn = 6.6667*(1.-f1)
  Go To 350
210 effn = 5.33333*(1.-f1)
  Go To 350
220 effn = 4.*(1.-f1)
  Go To 350
230 effn = 2.6667*(1.-f1)
  Go To 350
240 effn = 1.33333*(1.-f1)
  Go To 350
250 effn = 18.
  Go To 350
260 effn = 18. - (9.*f2)
  Go To 350
270 effn = 18.0*(1.-f2)
  Go To 350
280 effn = 18.*(1.-f2) - 1.*(9.*f3)
  Go To 350
290 effn = 18.*(1.-f2) - 2.*(9.*f3)
  Go To 350
300 effn = 18.*(1.-f2) - 3.*(9.*f3)
  Go To 350
310 effn = 18.*(1.-f2) - 4.*(9.*f3)
  Go To 350
320 effn = 18.*(1.-f2) - 45.0*f3
  Go To 350
330 no = 28 - i
  effn = no*1.8*(1.-3.*f3-f2)
  If (effn<0.) Then
    Go To 340
  Else
    Go To 350
  End If

340 eye = i
  effn = 60. - eye
! GUARD PACKAGE FOR I = 17
! PROBABLY UNNECESSARY
  If (effn<=0) effn = 1.0
350 Continue
  Return

End Function
Function grec(n, j, e, t)
  Common /dat/v(30), ea(30), s2(30), wave(220), e3(220), f(220), ll(30), &
    s3(30), s4(30), s5(30)
  Dimension :: drat(30)
  Data drat/2., .5, 2., .5, 6., 2.5, 2.22, 2.25, .67, .167, 2., .5, 6., 2.5, &
    2.22, 2.25, .67, .167, 12*0./

  deg = drat(n-j+1)
  x1 = e*11590./t
  x12 = x1*x1
  x14 = x12*x12
  g = 5.238E-23*t**1.5*deg*(s2(j)*x12+s4(j)*x12*x1+x12*x1*s5(j)*(0.5-x1/2.)+ &
    exint1(x1,2)*(s3(j)*x12*x1-x14*s4(j)+x1*x14*s5(j)/2.))
  grec = g
  If (x1>=300.) grec = 5.238E-23*t**1.5*deg*(s2(j)+s3(j)+s4(j)+s5(j))*x12
  Return

End Function
Function seaton(x, q)
! R. MOORE OCTOBER 1976
  Implicit Double Precision (A, B)
  Double Precision :: xs1, xs2, x3
  Data a10/8.7469604697013D-4/, b10/ -1.2007397521051D-4/
  Data a11/.20406771231267D0/, b11/ -.055223640825293D0/
  Data a12/ -2.1524482972354D0/, b12/.029171138841798D0/
  Data a13/12.663578339302D0/, b13/.25091093604147D0/
  Data a14/ -43.153566859883D0/, b14/ -.94344532109356D0/
  Data a15/61.113098339262D0/
  Data a20/.011834608639468D0/, b20/ -1.2467753947278D-3/
  Data a21/ -2.9889195903436D-3/, b21/ -.043871158058636D0/
  Data a22/ -.10946027271945D0/, b22/.013617064094285D0/
  Data a23/.097410292577482D0/, b23/ -5.2824884665512D-3/
  Data a24/ -.039676727608179D0/, b24/8.9490487211065D-4/
  Data a25/6.2318768197420D-3/
  Data a30/.018985613015081D0/, b30/ -.010963734653233D0/
  Data a31/ -.064707516794785D0/, b31/ -.028119928428050D0/
  Data a32/6.2172804659938D-3/, b32/1.1971209805431D-3/
  Data a33/ -4.1777961107942D-4/, b33/ -4.8739343472085D-5/
  Data a34/1.5264422582645D-5/, b34/8.4274360230135D-7/
  Data a35/ -2.2708774951499D-7/

  If (x>=20.) Go To 130
  If (x>2.) Go To 120
  If (x>0.2) Go To 110
  If (x>0.02) Go To 100
  xs1 = .4629*x*(1.+4.*x) - 1.0368*x**1.3333333333333*(1.+1.875*x)
  xs2 = -.0672*x*(1.+3.*x) + .1488*x**1.6666666666667*(1.+1.8*x)
  Go To 140
100 xs1 = a10 + (a11+(a12+(a13+(a14+a15*x)*x)*x)*x)*x
  xs2 = b10 + (b11+(b12+(b13+b14*x)*x)*x)*x
  Go To 140
110 xs1 = a20 + (a21+(a22+(a23+(a24+a25*x)*x)*x)*x)*x
  xs2 = b20 + (b21+(b22+(b23+b24*x)*x)*x)*x
  Go To 140
120 xs1 = a30 + (a31+(a32+(a33+(a34+a35*x)*x)*x)*x)*x
  xs2 = b30 + (b31+(b32+(b33+b34*x)*x)*x)*x
  Go To 140
130 x3 = 3.*x
  xs1 = -.1728*x**.33333333333333*(1.+(-8.+(70.+(-800.+11440./x3)/x3)/x3)/x3)
  xs2 = -.0496*x**.66666666666667*(1.+(-3.+(32.-448./x3)/x3)/x3)
140 x3 = x**.33333333333333*q**.66666666666667
  seaton = exint1(x, 3) + (xs1+xs2/x3)/x3
  Return

End Function
Function p(y, q)

  a = y/(q*q)
  If (a<=75) Go To 100
  p = 1./q
  Return
100 p = a*(1.-exint1(a,3))/q
  Return

End Function
Function gnd(num)

  gnd = 4.
  If (num<28) gnd = 3.
  If (num<10) gnd = 2.
  If (num<2) gnd = 1.
  Return

End Function
Function exint1(x, jump)
! R. MOORE OCTOBER 1976
! JUMP=1    EXINT1=E1(X)
! JUMP=2    EXINT1=EXP(X)*E1(X)
! JUMP=3    EXINT1=X*EXP(X)*E1(X)
  If (x>=1.) Go To 120
  exint1 = ((((((((7.122452D-7*x-1.766345D-6)*x+2.928433D-5)* &
    x-.0002335379D0)*x+.001664156D0)*x-.01041576D0)*x+ &
    .05555682D0)*x-.2500001D0)*x+.9999999D0)*x - log(x) - .57721566490153D0
  Go To (150, 100, 110) jump
100 exint1 = exp(x)*exint1
  Return
110 exint1 = x*exp(x)*exint1
  Return
120 x2 = x*x
  x3 = x2*x
  x4 = x3*x
  exint1 = (x4+8.5733287401D0*x3+18.059016973D0*x2+8.6347608925D0*x+ &
    .2677737343D0)/(x4+9.5733223454D0*x3+25.6329561486D0*x2+21.0996530827D0*x+ &
    3.9584969228D0)
  Go To (130, 140, 150) jump
130 exint1 = exint1*exp(-x)/x
  Return
140 exint1 = exint1/x
150 Return

End Function
Function etwo(x)

  etwo = exp(-x) - x*exint1(x, 1)
  Return
End Function

!---------------------------------------------------------------------------

Subroutine radiative_part(g_ksi, f_ksi, nfeld)
! **************************************************************************
! calculates the functions f(ksi) g(ksi), h(ksi)
! being ksi divided in NFELD parts, along Lc (cooling length),
! between 0 and 1. ksi is a parametrization of the radiative cooling
! length as introduced in Feldmeier et al. 1997
! **************************************************************************

  Use :: nlte_type
  Use :: nlte_xrays, Only: a
  Implicit None

  Integer (i4b), Intent (In) :: nfeld
  Real (dp), Dimension (nfeld), Intent (Out) :: g_ksi, f_ksi

  Real (dp), Parameter :: twoseventh = 2.D0/7.D0, third = 1.D0/3.D0

  Real (dp) :: xi, h_ksi

  Integer (i4b) :: i

! **************************************************************************

! integration from eps to 1


  Do i = 1, nfeld
    xi = (float(i)/float(nfeld))**twoseventh ! See that Xi never reachs the
!   value zero.
    h_ksi = a*xi*(1.0D0+(1.0D0/a-1.0D0)*xi)
    f_ksi(i) = 1.0D0/h_ksi
    g_ksi(i) = third*h_ksi*(4.0D0-h_ksi)
!   print*,xi,f_ksi(i),g_ksi(i)
  End Do

End Subroutine

!---------------------------------------------------------------------------

Subroutine adiabatic_part(theta, kappa_f, f_eta_rev, g_eta_rev, f_eta_for, &
  g_eta_for, nfeld)
! ***************************************************************************
! calculates the functions eta_f (for the forward shock)
! and eta_r (reverse shock).
! For details, see Feldmeier et al. 1997
! ***************************************************************************

  Use :: nlte_type
  Use :: nlte_xrays, Only: mfn
  Implicit None

  Integer (i4b), Intent (In) :: nfeld

! theta = T_r/T_f (just after jump), defined in the calling routine
  Real (dp), Intent (In) :: theta, kappa_f ! as defined by Feldmeier
  Real (dp), Dimension (nfeld), Intent (Out) :: f_eta_rev, g_eta_rev, &
    f_eta_for, g_eta_for

! local
  Integer (i4b) :: i

! Contact discontinuity inside the shell
! Location of the foward shock front (= 1.)
! Position of the reverse shock front
  Real (dp) :: eta_c, eta_f, eta_r, h_eta_rev, h_eta_for, xi

  Real (dp) :: eta_rev, eta_for


! ***************************************************************************


  eta_c = 1.D0 - 0.1314D0*kappa_f - 0.02857D0*kappa_f**2

  eta_f = 1.D0 !                    position of forward shock

! JO Sept. 2018: this is (according to Feldmeier) only valued for kappa_f <0.3
  If (kappa_f>=0.4) Then
    Write (999, *) ' kappa_f = ', kappa_f
    Write (999, *) &
      ' STOP: kappa_f > 0.4; either T_shock too large or v(l) to low'
    Print *, ' kappa_f = ', kappa_f
    Stop ' kappa_f > 0.4; either T_shock too large or v(l) to low'
  End If
  eta_r = eta_c - (1.D0-eta_c)*(1.D0-kappa_f+kappa_f**2)*theta**(0.5D0-0.31D0* &
    kappa_f+0.55D0*kappa_f**2)
! print*, eta_c, eta_f, eta_r

! JO Sept. 2018
  If (eta_r<0. .Or. eta_r>eta_c) Then
    Write (999, *) eta_c, eta_r
    Write (999, *) ' STOP: problems in double-shock structure: eta_r'
    Print *, eta_c, eta_r
    Stop ' problems in double-shock structure: eta_r'
  End If

! **************
! REVERSE SHOCKS
! **************
  Do i = 1, nfeld
    xi = float(i)/float(nfeld)
    eta_rev = eta_c - xi*(eta_c-eta_r)
    h_eta_rev = ((eta_c/eta_rev)**3-1.D0)/((eta_c/eta_r)**3-1.D0)
    f_eta_rev(i) = 1.D0/theta*(eta_r/eta_rev)**2*h_eta_rev**mfn
    g_eta_rev(i) = 1.D0/f_eta_rev(i)
!   print*, eta_rev, f_eta_rev(i), g_eta_rev(i)
!   print*, xi, f_eta_rev(i), g_eta_rev(i)
  End Do

! ************
! FOWARD SHOCKS
! ************
  Do i = 1, nfeld
    xi = float(i)/float(nfeld)
    eta_for = eta_c + xi*(1.D0-eta_c)
    h_eta_for = ((eta_c/eta_for)**3-1.D0)/((eta_c/eta_f)**3-1.D0)
    f_eta_for(i) = 1.D0*(eta_f/eta_for)**2*h_eta_for**mfn
    g_eta_for(i) = 1./f_eta_for(i)
!   print*, eta_for, f_eta_for(i), g_eta_for(i)
!   print*, xi, f_eta_for(i), g_eta_for(i)
  End Do

End Subroutine

!---------------------------------------------------------------------------

Subroutine cooling_flowing_time(fract_tc_tf, fract_lc_rs, tsl, rl, vl, rhol, &
  xnel, xnhl, xmu, theta, kappa_f, l)

! **************************************************************************
! calculates t_cool/t_flow and L_cool/r_sh,
! as well as kappa_f (Eq. 19, with vinf replaced by v(r)=vl
! as introduced in Feldmeier et al. 1997

! Note that as in Feldmeier, the jump-temperature as calculated in
! subroutine Tshock refers to the reverse shock (radiative case)
! and to the forward shock (adiabtic case, kappa_f needed)
! Thus, Tsh_rev = Tjump*theta , and Tsh_for = Tjump in the adiabtic case,
! whilst Tsh_rev = Tjump in the radiative case.
! **************************************************************************

  Use :: nlte_type
  Use :: nlte_var, Only: vmax, sr
  Use :: fund_const, Only: akb, amh
  Use :: nlte_xrays, Only: fx, px, fxl, rminx_eff, ar, fac

  Implicit None

  Real (dp), Intent (Out) :: fract_tc_tf, fract_lc_rs, kappa_f
  Real (dp), Intent (In) :: tsl, rl, vl, rhol, xnel, xnhl, xmu, theta
  Integer (i4b), Intent (In) :: l
! xnel and xnhl average values, re-correct for clumping


  Real (dp) :: tfpo

  fract_tc_tf = 40./7./ar*(16./3.)**1.5/4.**4*akb/(xmu*amh)
  fract_tc_tf = fract_tc_tf*(vl*vmax)/(rl*sr)*rhol/(xnel*xnhl)*tsl**1.5

! hack to ensure adiabatic treatment
! Fract_tc_tf=1.1

  fract_lc_rs = fac/ar*(16./3.)**2/4.**5*(akb/(xmu*amh))**1.5
  fract_lc_rs = fract_lc_rs/(rl*sr)*rhol/(xnel*xnhl)*tsl**2



! either use old approach with constant filling factor
! (giving rise to rho^2 dependent emissivities in all cases)
  If (px<-999.99) Then
    fxl(l) = fx
! or calculate new X-ray filling-factor according to work by Koushik Sen + JP,
! unifying the Feldmeier+ and Owocki+ (2013) work (see also Puls+ 2020)
  Else
    fxl(l) = fx*(rminx_eff/rl)**px*min(fract_lc_rs, 1.D0)

! JO: if thin-shell mixing shall be considered, replace amin1() by
!   1.d0/(1.d0+1./Fract_Lc_rs)^m  with mixing exponent m
  End If
! NOTE: this routine is only called for l=1,lxmin,
! and all fxl(l) for l>lxmin are set to zero

  If (fract_tc_tf<1.D0) Then
    kappa_f = 0.D0
    If (l==1) Write (999, 100) l, rl, tsl, fract_tc_tf, fract_lc_rs, fxl(l)
    Write (1, 100) l, rl, tsl, fract_tc_tf, fract_lc_rs, fxl(l)
    Write (*, 100) l, rl, tsl, fract_tc_tf, fract_lc_rs, fxl(l)
  Else
    tfpo = tsl !                    see comment above
    kappa_f = 1.D0 + sqrt(3.D0*xmu*amh/(16.*akb*tfpo))*vl*vmax
    kappa_f = 1.D0/kappa_f
    If (l==1) Write (1, 110) l, rl, tsl, fract_tc_tf, fract_lc_rs, fxl(l), &
      kappa_f
    Write (1, 110) l, rl, tsl, fract_tc_tf, fract_lc_rs, fxl(l), kappa_f
    Write (*, 110) l, rl, tsl, fract_tc_tf, fract_lc_rs, fxl(l), kappa_f
  End If

100 Format ('L =', I3, ' R = ', F9.3, ' Tshock = ', E9.3, ' tcool/tflow = ', &
    E9.3, ' Lcool/rsh = ', E9.3, ' fvol(l) = ', E9.3, ' RADIATIVE!')
110 Format ('L =', I3, ' R = ', F9.3, ' Tshock = ', E9.3, ' tcool/tflow = ', &
    E9.3, ' Lcool/rsh = ', E9.3, ' fvol(l) = ', E9.3, ' ADIABATIC! kappa = ', &
    F9.3)

End Subroutine

!---------------------------------------------------------------------------

Subroutine interpol_xray

! **************************************************************************
! interpolates lambda_xray (per bin) from equidistant energy grid
! onto nlte-freq. grid in such a way that integrals lambda_nu dnu
! remain preserved on output,
! integrals(k,l) contain emitted frequential X-ray energy
! **************************************************************************

  Use :: nlte_type
  Use :: nlte_dim, Only: nd => id_ndept
  Use :: fund_const, Only: clight
  Use :: nlte_xrays, Only: enerx, lambdax, integrals => lambdanu, lxmin, &
    lamred
  Use :: nlte_var, Only: fre, wfre, ifre

  Implicit None


  Real (dp), Allocatable, Dimension (:, :) :: lamint
  Real (dp), Allocatable, Dimension (:, :) :: dum

  Real (dp), Dimension (nd) :: eps0, intlow, intup

  Real (dp) :: dx, ered, x1, q, q1, h, eb, er, error

  Integer (i4b) :: nf, k, i, ired, kp, inext, l, kred, kblue


  If (lbound(enerx,1)/=0) Then
    Write (999, *) ' STOP: lower bound of enerx ne 0'
    Stop ' lower bound of enerx ne 0'
  End If

  nf = ubound(enerx, 1)

  If (1.D8/lamred<enerx(1)) Then
    Write (999, *) ' stop: interpol_xray: lamred outside xray grid'
    Stop ' interpol_xray: lamred outside xray grid'
  End If

  If (fre(ifre)>enerx(nf)) Then
    Write (999, *) ' stop: interpol_xray: fre(ifre) outside xray grid'
    Stop ' interpol_xray: fre(ifre) outside xray grid'
  End If

! test whether equidistant grid
  dx = enerx(1) - enerx(0)

  If (abs(1.-dx/(enerx(nf)-enerx(nf-1)))>1.D-12) Then
    Write (999, *) ' stop: enerx not equidistant'
    Stop ' enerx not equidistant'
  End If

  Allocate (lamint(0:nf,lxmin))
  Allocate (dum(ifre,lxmin))

  If (allocated(integrals)) Then
    Deallocate (integrals)
    Allocate (integrals(ifre,lxmin))
  Else
    Allocate (integrals(ifre,lxmin))
  End If

! indices i correspond to enerx (xray bins), and indices k to fre-grid

! integrate lambdax (sum up, since already energy per bin)

! summation might provide problems for very low numbers.
! tests showed that also quadruple precision and reverse summation does not
! help
! Anyway, problems only when emissivity very low.
  lamint(nf, :) = 0.D0
  integrals = 0.D0
  dum = 0.D0

! integrating from high to low frequencies, to avoid numerical problems
! note indices
  Do i = nf - 1, 0, -1
    lamint(i, :) = lamint(i+1, :) + lambdax(i+1, :)
  End Do

! check for precision
  Do l = 1, lxmin
    Do i = 0, nf - 1
      If (lamint(i,l)<=lamint(i+1,l)) Then
        If (lambdax(i+1,l)/=0.D0) Then
          Write (999, *) i, l, lamint(i, l), lamint(i+1, l), lambdax(i+1, l)
          Write (999, *) ' stop: interpol_xray: problems with precision'
          Print *, i, l, lamint(i, l), lamint(i+1, l), lambdax(i+1, l)
          Stop ' interpol_xray: problems with precision'
        End If
      End If
    End Do
  End Do

  ered = 1.D8/lamred

  Do i = 1, nf
    If (ered>=enerx(i-1) .And. ered<enerx(i)) Go To 100
  End Do
  Write (999, *) ' stop: 1 in interpol_xray'
  Stop ' 1 in interpol_xray'

100 Continue

! ered inside bin i-1,i, corresponding to energy i (see RS manual)
  ired = i

! calculate integrals between fre(k-1) and fre(k)

! lower freq. for first integral = ered
  q = (ered-enerx(ired))/(enerx(ired-1)-enerx(ired))
  q1 = 1.D0 - q

  intlow(1:lxmin) = q*lamint(ired-1, :) + q1*lamint(ired, :)

! define start index (upper freq. for first integral)
  Do k = 1, nf
    If (fre(k)>ered) Go To 110
  End Do
  Write (999, *) ' stop: 2 in interpol_xray'
  Stop ' 2 in interpol_xray'

110 Continue

  kp = k
  inext = ired

  Do k = kp, ifre !                 over all frequencies

    x1 = fre(k)
    Do i = inext, nf
      If (x1>=enerx(i-1) .And. x1<enerx(i)) Go To 120
    End Do
    Write (999, *) ' stop: 3 in interpol_xray'
    Stop ' 3 in interpol_xray'

120 Continue
!   interpolate energy-integrals at fre(k)

    inext = i !                     for next index

    q = (x1-enerx(inext))/(enerx(inext-1)-enerx(inext))
    q1 = 1.D0 - q

    intup(1:lxmin) = q*lamint(inext-1, :) + q1*lamint(inext, :)

!   do l=1,lxmin
!   print*,l,intlow(l),intup(l)
!   enddo

!   reverse order, since integration from high to low frequencies
!   JO Sept 2020: bug found by Sarah Brands removed
!   integrals(k,:)=intlow-intup
    integrals(k, :) = intlow(1:lxmin) - intup(1:lxmin)
    If (minval(integrals(k,:))<0.D0) Then
      Write (999, *) ' stop: interpol_xrays: integrals negative'
      Stop ' interpol_xrays: integrals negative'
    End If

    intlow = intup

  End Do !                          all frequencies

! total energy
  Do l = 1, lxmin
    h = sum(integrals(kp:ifre,l))
    If (h>lamint(1,l)) Then
      Write (999, *) ' stop: interpol_xrays: summ(energies) > lamint!'
      Stop ' interpol_xrays: summ(energies) > lamint!'
    End If
    eps0(l) = h !                   slightly different from lamint, due to
!   edge effects
  End Do


! this is now the frequential energy on the coarse grid w.r.t. staggered
! frequency grid
  k = kp
  integrals(k, :) = integrals(k, :)/(clight*(fre(k)-ered))
  Do k = kp + 1, ifre
    integrals(k, :) = integrals(k, :)/(clight*(fre(k)-fre(k-1)))
  End Do

! now, interpolate frequential energy onto original grid, conserving the
! energy
! (see notes)
! first frequency
  k = kp
  q = (fre(k)-ered)/(fre(k+1)-ered)
  q1 = (fre(k+1)-fre(k))/(fre(k+1)-ered)
  dum(k, :) = q*integrals(k, :) + q1*integrals(k+1, :)

  Do k = kp + 1, ifre - 1
    q = (fre(k)-fre(k-1))/(fre(k+1)-fre(k-1))
    q1 = (fre(k+1)-fre(k))/(fre(k+1)-fre(k-1))
    dum(k, :) = q*integrals(k, :) + q1*integrals(k+1, :)
  End Do

! last frequency (no action required)
  dum(ifre, :) = integrals(ifre, :)

! integrals now contain the frequential energies on the original grid
  integrals = dum


! --------------------------------------------------------------------------
! Now lets test a specific, line contaminated interval, say between 10 to 20 A
! corresponding to 0.6 to 1.2 keV
! can be skipped if everything is running

  Do k = kp, ifre - 1
    If (1.D8/fre(k)>20.D0 .And. 1.D8/fre(k+1)<=20.D0) Exit
  End Do
  kred = k
  Do k = kred, ifre - 1
    If (1.D8/fre(k)>10.D0 .And. 1.D8/fre(k+1)<=10.D0) Exit
  End Do
  kblue = k + 1
  er = fre(kred)
  eb = fre(kblue)

! at first, exact integral
  Do i = 1, nf
    If (er>=enerx(i-1) .And. er<enerx(i)) Exit
  End Do
! er inside bin i-1,i, corresponding to energy i (see RS manual)
  ired = i
! lower freq. for first integral = ered
  q = (er-enerx(ired))/(enerx(ired-1)-enerx(ired))
  q1 = 1.D0 - q
  intlow(1:lxmin) = q*lamint(ired-1, :) + q1*lamint(ired, :)


  Do i = ired, nf
    If (eb>=enerx(i-1) .And. eb<enerx(i)) Exit
  End Do
! eb inside bin i-1,i, corresponding to energy i (see RS manual)
  ired = i
! lower freq. for first integral = ered
  q = (eb-enerx(ired))/(enerx(ired-1)-enerx(ired))
  q1 = 1.D0 - q
  intup(1:lxmin) = q*lamint(ired-1, :) + q1*lamint(ired, :)

  intup = -(intup-intlow)

! now, integral over interpolated values
  intlow(1:lxmin) = 0.D0
  Do k = kred, kblue - 1
    intlow(1:lxmin) = intlow(1:lxmin) + 0.5D0*(integrals(k,:)+integrals(k+1,:) &
      )*(fre(k+1)-fre(k))*clight
  End Do

  error = 0.D0
  Do l = 1, lxmin
    If (intup(l)/=0.D0) error = max(error, abs(1.D0-intlow(l)/intup(l)))
  End Do

  Print *
  Print *, 'test for interval = ', 1.D8/er, 'to', 1.D8/eb, ' A'
  Print *, 'energy not conserved by ', error
  If (error>0.05) Then
    Write (999, *) &
      ' stop: interpol_xrays: energy of interpolated values not conserved(1)!'
    Stop ' interpol_xrays: energy of interpolated values not conserved(1)!'
  End If

! --------------------------------------------------------------------------
! energy conservation in total freq. range, and
! renormalization w.r.t. standard integration weights

  error = 0.D0
  Do l = 1, lxmin
    h = 0.D0
    Do k = kp, ifre
      h = h + wfre(k)*integrals(k, l)
    End Do
    If (eps0(l)/=0.D0) Then
      error = max(error, abs(1.D0-h/eps0(l)))
      integrals(kp:ifre, l) = integrals(kp:ifre, l)*eps0(l)/h
    Else
      If (h/=0.D0) Then
        Write (999, *) ' stop: interpol_xrays: eps0 = 0 and h ne 0!'
        Stop ' interpol_xrays: eps0 = 0 and h ne 0!'
      End If
    End If
  End Do

  Print *
  Print *, 'test for total range = ', 1.D8/ered, 'to', 1.D8/fre(ifre), ' A'
  Print *, 'energy not conserved by ', error
  If (error>0.05) Then
    Write (999, *) &
      ' stop: interpol_xrays: energy of interpolated values not conserved(2)!'
    Stop ' interpol_xrays: energy of interpolated values not conserved(2)!'
  End If

  Write (999, *)
  Write (999, *) ' interpolation of x-ray emission to coarse grid done'
  Write (999, *)
  Print *
  Print *, ' interpolation of x-ray emission to coarse grid done'
  Print *

! can be deleted when everything is running
  Do k = 1, kp - 1
    Do l = 1, lxmin
      If (integrals(k,l)/=0D0) Then
        Write (999, *) ' stop: problem in integrals'
        Stop ' problem in integrals'
      End If
    End Do
  End Do

  Deallocate (lamint, dum)

  Return

! for tests at l
  l = 10
  Do i = 1, nf
    Write (999, *) i, 1.D8/(enerx(i)-0.5*dx), lambdax(i, l)/(clight*dx)
    Print *, i, 1.D8/(enerx(i)-0.5*dx), lambdax(i, l)/(clight*dx)
  End Do

  Do k = kp, ifre
    Write (999, *) k, 1.D8/fre(k), integrals(k, l)
    Print *, k, 1.D8/fre(k), integrals(k, l)
  End Do
  Stop

End Subroutine

!---------------------------------------------------------------------------

Subroutine read_kshell_data

  Use :: nlte_type
  Use :: fund_const, Only: hh, clight, ev
  Use :: fastwind_params, Only: fpath

  Use :: nlte_xrays, Only: n_kedges, z, n, eth, sigma, s, zeff, aug_1_6

  Implicit None

  Integer :: i

  Real (dp) :: zl10, e_thre, cross_s, s_calc, const, eth_kaastra, err, summ

  Open (Unit=2, File=fpath//'k_shell_Auger_data', Status='old')
  Read (2, *)
  Read (2, *)
  Read (2, *)

! aug_1_6 contains fraction of number of electrons ejected
! during complete process (1 K-shell electron + n-1 Auger electrons)

  err = 0.D0
  Do i = 1, n_kedges
!   N is number of electrons
    Read (2, *) z(i), n(i), eth(i), eth_kaastra, sigma(i), s(i), zeff(i), &
      aug_1_6(i, 1), aug_1_6(i, 2), aug_1_6(i, 3), aug_1_6(i, 4), &
      aug_1_6(i, 5), aug_1_6(i, 6)
!   sigma(i)=0.0
    err = max(abs(1.-eth(i)/eth_kaastra), err)
    summ = sum(aug_1_6(i,:))
    If (summ/=10000.) Then
      Write (999, *) ' stop: error in Auger probabilities!'
      Stop ' error in Auger probabilities!'
    End If

!   **************************************************
!   Remember to multiply sigma(i) by 1.E-18
!   **************************************************
!   print*, z(i),n(i),eth(i),sigma(i),s(i),zeff(i)
  End Do
  Close (2)

  Go To 100
! for tests
  i = 12
  zl10 = log10(float(z(i)))
  e_thre = 10.D0**((0.24D0-0.43D0/zl10)*log10(float(n(i)))+2.D0*zl10+ &
    1.133539D0)
  cross_s = 10.D0**(-log10(eth(i))+2.47D0)*1.D-18
  s_calc = 10.D0**(0.07D0*log10(eth(i))+0.221D0)

  Print *, e_thre, cross_s, s_calc
  Print *
  Print *, eth(i), sigma(i)*1.D-18, s(i)


100 Continue
! now calculate ionization stages (astronomical units)
  n = z - n + 1
  sigma = sigma*1.D-18

! transform eV to Kayser
  const = hh*clight/ev
  const = 1.D0/const
  eth = eth*const

  Write (999, *)
  Write (999, *) ' K-shell data read'
  Write (999, *) &
    ' max deviation Daltabuit & Cox vs. Kaastra & Mewe (thresholds) = ', err
  Print *
  Print *, 'K-shell data read'
  Print *, 'max deviation Daltabuit & Cox vs. Kaastra & Mewe (thresholds) = ', &
    err
  If (err>5.D-2) Then
    Write (999, *) ' stop: deviation too large!'
    Stop ' deviation too large!'
  End If
  Write (999, *)
  Print *

  aug_1_6 = aug_1_6/10000.

  Return

End Subroutine

!---------------------------------------------------------------------------

Subroutine opacit_kshell(nd, xne, ilow, imax)
! calculates k-shell opacities (clumping included in occupation numbers),
! but not corrected

  Use :: nlte_type

  Use :: nlte_dim, Only: id_atoms, id_ndept

  Use :: princesa_var, Only: nat, labat, zeff

  Use :: nlte_var, Only: concon, imia, imaa, imianl, imaanl, enionnd, ifre, &
    fre

  Use :: nlte_app, Only: natom, indexel, jatom_full, met_imin, met_imax, &
    abund, fjk, xmuee, summas

  Use :: nlte_xrays, Only: n_kedges, k_nmin, k_nmax, z, nionk => n, eth, &
    sigma, s, name, opa_kshell, alpha_fre

  Implicit None

  Integer (i4b), Parameter :: kel = id_atoms, nd1 = id_ndept

  Integer (i4b), Intent (In) :: nd
  Real (dp), Intent (In) :: xne(nd)

  Integer (i4b), Intent (In) :: ilow(nd1, kel), imax(nd1, kel)

! real(dp), dimension(:), allocatable :: alpha_fre

  Real (dp), Dimension (nd1) :: const, occ

  Integer (i4b) :: nk, zold, k, kex, iex, iz, j, ll, kk, kmin

  Integer (i4b), Dimension (natom) :: iminbg, imaxbg
  Integer (i4b), Dimension (kel) :: iminex, imaxex

  Real (dp) :: eth_min

! minimum and maximum ions for bg-elements
  Do k = 1, natom
    If (jatom_full(k)/=1) Cycle
    iminbg(k) = minval(met_imin(k,:))
    imaxbg(k) = maxval(met_imax(k,:)) + 1 ! accounting for highest ion with
!   ground-state
  End Do

! minimum and maximum ions for explicit elements
  Do k = 1, nat
    iminex(k) = minval(ilow(:,k))
    imaxex(k) = maxval(imax(:,k)) + 1 ! accounting for highest ion with
!   ground-state
    If (.Not. concon) Then
      If (iminex(k)/=imia(k)) Then
        Write (999, *) ' stop: opacit_kshell: imia not consistent'
        Stop ' opacit_kshell: imia not consistent'
      End If
      If (imaxex(k)-1/=imaa(k)) Then
        Write (999, *) ' stop: opacit_kshell: imaa not consistent'
        Stop ' opacit_kshell: imaa not consistent'
      End If
    Else
      If (iminex(k)/=imianl(k)) Then
        Write (999, *) ' stop: opacit_kshell: imianl not consistent'
        Stop ' opacit_kshell: imianl not consistent'
      End If
      If (imaxex(k)-1/=imaanl(k)) Then
        Write (999, *) ' stop: opacit_kshell: imaanl not consistent'
        Stop ' opacit_kshell: imaanl not consistent'
      End If
    End If
  End Do

  eth_min = 1.D99
  Do k = 1, n_kedges
    If (nionk(k)<k_nmin) Cycle !    below k_nmin (typically 3, set in
!   nlte_xrays)
    If (nionk(k)>k_nmax) Cycle !    above k_nmax (typically 8, set in
!   nlte_xrays)
    If (eth(k)>2.D7) Cycle !        below 5 a (maximum energy in fre, set in
!   frescal)
    eth_min = min(eth_min, eth(k))
  End Do

  Do kk = 1, ifre
    If (eth_min==fre(kk)) Go To 100
  End Do
  Write (999, *) ' stop: min edge not found in opacit_kshell'
  Stop ' min edge not found in opacit_kshell'


100 kmin = kk

  If (allocated(opa_kshell)) Then
    Deallocate (opa_kshell, alpha_fre)
    Allocate (opa_kshell(nd1,kmin:ifre), alpha_fre(kmin:ifre))
  Else
    Allocate (opa_kshell(nd1,kmin:ifre), alpha_fre(kmin:ifre))
  End If
  opa_kshell = 0.D0

  nk = 0
  zold = 0
! only edges from ionization stage iii on
  Do k = 1, n_kedges
    If (nionk(k)<k_nmin) Cycle !    below k_nmin (typically 3, set in
!   nlte_xrays)
    If (nionk(k)>k_nmax) Cycle !    above k_nmax (typically 8, set in
!   nlte_xrays)
    If (eth(k)>2.D7) Cycle !        below 5 a (maximum energy in fre, set in
!   frescal)
    zold = z(k)

!   different treatment of explicit and background element

!   bg element
    If (indexel(zold)==-1) Then

      If (jatom_full(zold)/=1) Then
        Write (999, *) &
          ' stop: K-shell absorption from not selected bg-element!'
        Stop ' K-shell absorption from not selected bg-element!'
      End If
      If (name(zold)/=name(zold)) Then
        Write (999, *) ' stop: opacit_kshell: names not consistent (1)!'
        Stop ' opacit_kshell: names not consistent (1)!'
      End If
!     print*,k,z(k),nionk(k),iminbg(zold),imaxbg(zold)
      If (nionk(k)<iminbg(zold) .Or. nionk(k)>imaxbg(zold)) Cycle
      Write (999, *) ' k-shell opacity from bg element', zold, name(zold), &
        nionk(k)
      Print *, ' k-shell opacity from bg element', zold, name(zold), nionk(k)

      const = xmuee*xne/summas*abund(zold) ! depth dependent
      occ = 0.D0
      Do ll = 1, nd
        j = nionk(k)
        If (j<met_imin(zold,ll) .Or. j>met_imax(zold,ll)+1) Cycle
        occ(ll) = const(ll)*fjk(zold, j, ll)
      End Do

!     explicit element
    Else
      kex = indexel(zold)
      If (name(zold)/=labat(kex)) Then
        Write (999, *) ' stop: opacit_kshell: names not consistent (2)!'
        Stop ' opacit_kshell: names not consistent (2)!'
      End If
!     account for zeff
      iz = int(zeff(kex))
      If (float(iz)-zeff(kex)>1.D-12) Then
        Write (999, *) ' stop: opacit_kshell: something wrong with int(zeff)!'
        Stop ' opacit_kshell: something wrong with int(zeff)!'
      End If
      iex = nionk(k) - iz !         this is the index w.r.t. ionization state
!     of expl. elements
      If (iex<iminex(kex) .Or. iex>imaxex(kex)) Cycle
      Write (999, *) ' k-shell opacity from exp. elem.', kex, labat(kex), &
        nionk(k), '(', iex, ')'
      Print *, ' k-shell opacity from exp. elem.', kex, labat(kex), nionk(k), &
        '(', iex, ')'
      occ = 0.D0
      Do ll = 1, nd
        If (iex<ilow(ll,kex) .Or. iex>imax(ll,kex)+1) Cycle
        occ(ll) = enionnd(kex, iex, ll)
      End Do
    End If

    Do kk = 1, ifre
      If (eth(k)==fre(kk)) Go To 110
    End Do
    Write (999, *) ' stop: edge not found in opacit_kshell'
    Stop ' edge not found in opacit_kshell'

110 Continue

    alpha_fre(kk:ifre) = sigma(k)*(eth(k)/fre(kk:ifre))**s(k)
    Do ll = 1, nd
      If (occ(ll)/=0.D0) opa_kshell(ll, kk:ifre) = opa_kshell(ll, kk:ifre) + &
        occ(ll)*alpha_fre(kk:ifre)
    End Do

    nk = nk + 1
  End Do

  Write (999, *)
  Write (999, *) nk, ' k-shell opacities calculated'
  Write (999, *)
  Print *
  Print *, nk, ' k-shell opacities calculated'

  Return
End Subroutine

!-----------------------------------------------------------------------

Subroutine intermepho_auger(x1, itrans, ll)

! calculates photo-integrals for k-shell absorption required for Auger-ionization

  Use :: nlte_type
  Use :: nlte_dim
  Use :: fund_const
  Use :: nlte_var, Only: fre, wfre, ifre, xj
  Use :: nlte_xrays, Only: n_kedges, eth, sigma, s

  Implicit None
  Integer (i4b), Intent (In) :: itrans, ll
  Real (dp), Intent (Out) :: x1

  Integer (i4b) :: k, kk

  Real (dp) :: wf, quant, xint, const

  If (itrans<=0 .Or. itrans>n_kedges) Then
    Write (999, *) ' stop: intermepho_auger: error in itrans'
    Stop ' intermepho_auger: error in itrans'
  End If

  Do kk = 1, ifre
    If (eth(itrans)==fre(kk)) Go To 100
  End Do
  Write (999, *) ' stop: intermepho_auger: edge not found'
  Stop ' intermepho_auger: edge not found'

100 Continue

  xint = 0.

  Do k = kk, ifre

!   ---    note: first integration weight has to be modified

    If (k==kk) Then
      wf = wfre(kk) - .5D0*clight*(fre(kk)-fre(kk-1))
    Else
      wf = wfre(k)
    End If
    quant = sigma(itrans)*(eth(itrans)/fre(k))**s(itrans)/fre(k)*xj(ll, k)
    xint = xint + wf*quant

  End Do

  const = 4.D0*pi/(clight*hh)
  x1 = const*xint

  If (x1<=0.) Then
    Write (999, *) ' stop: INTERMEPHO_AUGER!'
    Stop ' INTERMEPHO_AUGER!'
  End If

  Return
End Subroutine
