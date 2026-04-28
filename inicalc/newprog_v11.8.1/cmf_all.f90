Module version_cmf_all

! NOTES for further modifications:

! NOTE: so far, everything with constant vturb!
! minimum vturb = vturbmin (5 km/s); this value is used even if input lower!
! modify vturbmin if you really need lower values (e.g., for grad)

! IMPORTANT NOTE (referring to subr. line_list)
! lines from background elements de-packed (all multiplet components!)
! example: CIV 1548, 1550 with same index (id,id1)
! vs.
! lines from explicit elements (mostly) packed (from Detail input)
! specific lines are depacked, according to 'linetab' (LINES_xxx.DAT) from
! ATOM_FILE

! IMPORTANT NOTE (referring to subr. PROFILE_CMF)
! so far, only Doppler profiles. At least for H/He, Stark broadening should
! be included to obtain more realistic fluxes at Lyman series limit

! NOTE regarding opacities: for selected (bg-)elements, opacities and
! source functions (OPACITY_CMF) only calculated in between
! met_imin(l), met_imax(l) [controlled by occngold ne 0]
! other opacities set to zero.
! The same also applies for line-rate coefficients incl. ALO (CALC_TLU_MET)
!
! Only for Ng-extrapolation, all defined (resonance) transitions
! are accounted for, by using occng and occngold2
! [however, since rateeq and opacities/source-functions use only
! transitions defined from occngold ne 0, additional info might be not used]

! For non-selected (bg-elements), opacities/source-functions calculated
! from occng, in the spirit of nlte_approx.

! THUS, CAREFUL INSPECTION OF MET_IMIN, MET_IMAX is required in certain cases.
! E.g., the overlap between the OIII/NIII resonance lines might be hampered
! in hotter models, since: outside OIII is present, in the middle wind not
! (OIV-OVI), in the sonic region OIII is present again, and then
! vanishes in the deeper photosphere.

! WRITTEN:
! 04/12 by JP

! HISTORY:
! comletely new

! version 0.0 April 2012
! version 0.1 Oct   2013
! version 1.0 June  2015
! version 1.1 April 2016
! version 1.2 April 2016
! version 1.3 Oct   2016

! version 1.3.1 end Oct 2016: Bug in u0 found and fixed (should play almost no
!   role)

! version 1.3.2 Dec 2016: Sobo-transfer for lines inside range (outside ilow,
!   imax)

! version 1.3.3 Dec 2016: inclusion of depacked components as defined in
!   LINES_xxx.dat

! version 1.3.4 Feb 2017: inclusion of Ng-extrapolation for selected
!   transitions

! version 1.3.5 March/April 2017: lots of convergence control statements

! version 1.3.6 Dec 2017: subr. optcmf_restart changed: interpolation of XJ
!   from previous results, if ifre1 ne ifre (change of ionization,
!   due to different ilow, imax (subr. ilowimax_lte)

! version 1.3.7 March 2018: after including X-rays, two stop statements
!   related to erroneous Edd. factors
!   (subr. interpol_cmf_to_coarse, subr. xj_smooth) changed:
!   check only for lambda > 20 A, otherwise numerical
!   uncertainties too large.
!   Two similar stop statements not manipulated, since
!   in effect only for CMF-range.
!   betacic (no upper level) in OPACITY_CMF now always
!   calculated from consistent Edd. factor, and not
!   approximative for L < LTHIN
!   (might be reconsidered if problems with certain lines transitions)

! version 1.3.8 July 2018 final ALO control:
!   expl. elements (calc_tlu)
!   if ALO < 0.1 then set to zero (calc_tlu)
!   background elements (calc_tlu_met)
!   if ALO < 0.9 then set to zero,
!   except for OIII 374 (1-10)
!   if ALO > 0.99, then set to 0.99

!   subr. line list modified, to allow for changes
!   in ionization structure during convergence

! version 1.3.9 Feb 2019: new treatment of specific transitions with UP=0
!   (see nlte.f90), if OPT_NEW_UPPER = .TRUE.
!   This version has already been implemented in Nov. 2018!

! version 1.4.0 Oct 2020: compatible with gfortran

! version 1.4.1 Dec 2020: first IR version

! version 1.4.2 April 2021: several updates, in particular regarding ALO of
!   HeII 303

! version 1.4.3 March 2022: two stop statements after checking the
!   eddington factor removed and replaced by warnings
!   changes in cmf_complete: interpolate xj_save etc.
!   when frequency grid has changed.
!   no stop for erroneous f-value (Adi) for O3_28 to O3_39

! version 1.4.4 April 2023: few changes in ray_complete when encountering u<0
!   update of chih_params1 to allow for multi-linear
!   regression when "method 2" is of lesser precision

! version 1.4.5 Sept 2023: previous erroneous f-value (Adi) for O3_28 to O3_39
!   corrected, no longer checked.

! version 1.4.6 Nov 2023: smoothing of xj_save now in separate subroutine
!   prep_smooth (before opacity_cmf), instead inside
!   cmf_complete. Irradiation of lines with up=0
!   with smoothed XJ (in opacity_cmf)

! version 1.4.6.1 Feb 2025: Changed back to irradiation with pseudo-continuum,
!   since otherwise too strong peaks around 230 A for hot models

! version 1.4.7 March 2025: completely new philosophy to account for an
!   (almost) consistent ALO after updating J via moments.
!   From now on, no longer any warnings from calc_tlu_met
!   because of inconsistent ALOs ("problems with snew" etc.)
!   Moreover, after zillions of tests, original formulation
!   for the treatment of negative I in ray_complete have
!   shown that the original formulation was OK.
!   Subroutine ray_complete cleaned for unnecessary statements.

! version 1.4.8 July 2025: couple of changes at blue-most frequ. in CMF
!   treatment (k=1), to obtain a consistent description
!   for J and H at the outer boundary. Most changes
!   refer to an approximate Iminus which was set to zero
!   in previous versions. Moreover, to avoid other
!   inconsistencies, now kmom_start = 2, and J and H
!   = J_ray and H_ray for k=1,2

! version 1.4.9 Oct 2025: new approach to handle negative intensities in
!   subr. ray_complete (cmf_all.f90), together with
!   consistent changes in moments_cmf. Now, the number
!   of inconsistent solutions for Jbar from moments_cmf
!   should considerably decrease, as well as the
!   "overall" changes in J_nu from iteration to
!   iteration. Requires new variable add_corr_arr.
!   Note: only tested with moments_cmf1. If moments_cmf
!   itself should be used, check and test.

! version 1.5 Nov 2025: polished
  
  Character (10) :: ver_cmf_all = '1.5'

End Module

!--------------------------------------------------------------------------- &

Module cmf_all_var

  Use :: nlte_type
  Use :: nlte_dim, Only: id_nttrd, id_atoms, id_ndept, id_llevs
  Use :: fastwind_params, Only: natom

  Implicit None

  Integer (i4b) :: ntotnl_all, ntot, ntot_expl_inside, icfirst, iclast
! ntot: all lines in new line list
! ntot_expl_inside:
! lines from expl. elem. inside wavblue, wavred in new line list

  Integer (i4b) :: nftot = 0
! total number of frequency points within wavblue-dlam1, wavred+dlam2

  Integer (i4b) :: no2
! number of overlapping lines with up=0 as calculated in subr. overlap_noup

  Real (dp) :: vref, vturb_cmf

  Real (dp), Dimension (id_nttrd) :: xlamcmf_depacked, sumgf_depacked
! analogous to xlamcmf, but with wavelength of strongest depacked component
  Integer (i4b), Dimension (id_nttrd) :: indxlamc_depacked
! analogous to indxlamc (index for wavelength sorting), refering to
! xlamcmf_depacked

  Integer (i4b), Dimension (id_nttrd) :: indexrbb_inv, indexcmf1, &
    indexcmf1_lam
! indexrbb_inv gives index regarding index1 (rbb+cbb) for j=1,id_nttrd
! indexcmf1 refers to cmf-treatment (1,2,3,4,...)
! indexcmf1_lam points to index of lam_cmf
  Integer (i4b), Dimension (id_atoms) :: indexel_inv
! gives index of explict atom regarding background elements
  Integer (i4b), Dimension (natom) :: elements_cmf, nf_cmf
! provides info on which bg-elements are treated in cmf, and gives maximum
! number
! of freq. points (one sided, without 0) for Doppler profiles w.r.t. vref

  Real (dp), Dimension (natom) :: vth_cmf
! using vturb_cmf=max(vturb,vturbmin)
  Real (dp), Dimension (id_ndept, natom) :: vdoptot

  Real (dp), Dimension (id_llevs, id_ndept) :: occexpl

  Integer (i4b), Dimension (:), Allocatable :: id1, id_rbb, index_id1
  Integer (i4b), Dimension (:), Allocatable :: id1_overlap
  Real (sp), Dimension (:), Allocatable :: gf1, xlam1, opal_ntest

  Integer (i4b), Dimension (:), Allocatable :: index_lam_cmf, indfre_cmf, &
    indfre_staggered_cmf
! indfre_cmf: index regarding cont. grid (for interpolation)
! indfre_staggered_cmf: index regarding cont. grid (staggered)
  Real (dp), Dimension (:), Allocatable :: lam_cmf, h_obs, xhinner_cmf
  Real (dp), Dimension (:, :), Allocatable :: sumopal_cmf, sumetal_cmf, &
    xj_cmf, alo_cmf
  Real (dp), Dimension (:, :), Allocatable :: xh_cmf, fedd, gedd, geddp, &
    add_corr_arr
  Real (dp), Dimension (:), Allocatable :: hbound1, nbound1, hboundnd, &
    nboundnd
  Real (dp), Dimension (:, :), Allocatable :: u_outer
  Real (dp), Dimension (:), Allocatable :: sumopal_outer, sline_outer
  Real (dp), Dimension (:), Allocatable :: sumtaus_outer

  Real (dp), Dimension (:, :), Allocatable :: profdop_h, profdop_he
  Real (dp), Dimension (:, :), Allocatable :: pweightdop_h, pweightdop_he
  Real (dp), Dimension (:, :, :), Allocatable :: profdop
  Real (dp), Dimension (:, :, :), Allocatable :: pweightdop

  Real (dp), Dimension (:), Allocatable :: sum_pweightdop_h, sum_pweightdop_he
  Real (dp), Dimension (:, :), Allocatable :: sum_pweightdop

  Real (dp), Dimension (:), Allocatable :: lamtest, weighttest

  Integer (i4b), Dimension (:, :), Allocatable :: ice
  Real (dp), Dimension (:, :), Allocatable :: icearr

  Integer (i4b) :: lineno

  Type (depacked_data), Dimension (:), Allocatable :: depacked

End Module

!--------------------------------------------------------------------------- &

Module ng_var

  Use :: nlte_type
  Implicit None

! initialization: after restart etc., allow for one iteration without
! saving source functions, to obtain clear old and new occupation numbers
! THUS: itng=-1, itng_met=-1
! to avoid wrong n_ng (number of lines to be extrapolated) = 0 condition,
! initialize as negative
  Integer (i4b) :: itng = -1, n_ng = -1, itng_met = -1, n_ng_met = -1

  Logical :: ng = .False., ng_met = .False.

  Integer (i4b), Dimension (:), Allocatable :: index1_ng, index_ng
  Integer (i4b), Dimension (:, :), Allocatable :: trans_ng
  Real (dp), Dimension (:, :, :), Allocatable :: sl_ng

  Integer (i4b), Dimension (:), Allocatable :: index_ng_met
  Integer (i4b), Dimension (1000) :: trans_ng_met
  Real (dp), Dimension (:, :, :), Allocatable :: sl_ng_met

End Module

!-----------------------------------------------------------------------

Subroutine cmf_all(nd, xne, temp, clfac, r, velo, dvdr, rho, taur, xnh, &
  pressure, ferr, dtfcorr, ilow, imax, rmax, unasol, opt_ray)

! calculates formal integral in the cmf for *all* important lines, i.e.,
! for lines with CALCION(k,j) = 1, including also lines with only the lower
! level present

! all tests so far with GLOBAL = .TRUE.
  Use :: nlte_type
  Use :: nlte_dim, Only: id_ndept, id_atoms
  Use :: nlte_opt, Only: opt_modify_sl_noup
  Use :: fastwind_params, Only: natom, taur_test, test_overlap
  Use :: nlte_var, Only: metals_converged, vsound, vmax, lwion1
  Use :: nlte_app, Only: indexel, lniii_in, lniii_out

  Implicit None

  Integer (i4b), Parameter :: nd1 = id_ndept
  Integer (i4b), Parameter :: kel = id_atoms

! ..
  Integer (i4b), Intent (In) :: nd
  Logical, Intent (In) :: unasol
  Integer (i4b), Dimension (nd1, kel), Intent (In) :: ilow, imax
  Real (dp), Intent (In) :: rmax
! ..
! .. array arguments ..
  Real (dp), Dimension (nd1), Intent (In) :: xne, temp, clfac, r, velo, dvdr, &
    rho, xnh, pressure
! JO, changed April 2016
  Real (dp), Dimension (nd1), Intent (Inout) :: taur
  Real (dp), Dimension (nd1), Intent (Out) :: ferr, dtfcorr

  Integer (i4b) :: l, ntest, lsound

  Logical :: start, opt_ray

  Data start/.True./

  Write (999, *)
  Write (999, *) ' CMF_ALL STARTED'
  Write (999, *)
  Print *
  Print *, ' CMF_ALL STARTED'
  Print *

  If (start) Then

!   JO Sept 2018
!   define region where NIII resonance lines are formed (when important)
    Do l = 1, nd - 1
      If (velo(l)*vmax<vsound) Exit
    End Do
    lsound = l
    lniii_in = lwion1 - 1
    Do l = 1, nd - 1
      If (velo(l)<0.1) Exit
    End Do
    l = l - 1
!   use max(0.1 vinf,v(lsound-10)) as outer boundary
    lniii_out = min(lsound-10, l)
    If (lniii_out<1) Then
      Write (999, *) ' stop: lniii_out < 1'
      Stop ' lniii_out < 1'
    End If
    Write (999, *) 'approx. range of NIII emission line-formation (if any):'
    Write (999, *) 'lniii_out = ', lniii_out, ' lniii_in = ', lniii_in
    Write (999, *)
    Print *, 'approx. range of NIII emission line-formation (if any):'
    Print *, 'lniii_out = ', lniii_out, ' lniii_in = ', lniii_in
    Print *

    Call line_list(test_overlap)
!   merge (explicit and background elements) and modify the original
!   line-lists
!   NOTE: ONLY REQUIRED AS LONG AS IMIN, IMAX ARE CHANGING (start=.true.)
    Call fgrid_cmf
    Call profile_cmf(nd, temp)
    If (metals_converged) start = .False.
  End If

  ntest = 1
  If (test_overlap) Then
    Do l = 1, nd - 1
      If (taur(l)<=taur_test .And. taur(l+1)>taur_test) Exit
    End Do
    ntest = l
    If ((taur(l+1)-taur_test)<abs(taur_test-taur(l))) ntest = l + 1
    Print *
    Print *, ' ntest = ', ntest, ' taur(ntest) = ', taur(ntest)
    Print *
  End If

! JO: Nov. 2023: new subroutine for smoothing, extracted from cmf_complete
  Call prep_smooth

  Call opacity_cmf(nd, xne, temp, clfac, r, velo, dvdr, test_overlap, ntest, 1)

  If (test_overlap) Then
    Call overlap_detector(dvdr(ntest))
    Write (999, *) ' stop: stop after overlap_detector'
    Stop ' stop after overlap_detector'
  End If

! calculate and apply modified source-functions for bg-lines with no upper
! level, but only if non-HHe explicit elements are present
  If (opt_modify_sl_noup .And. maxval(indexel(3:natom))>-1) Then
!   detects lines with up=0 which need to be modified
    Call overlap_noup(clfac)

!   2nd run, source-functions (and sumetal_cmf!) will be modified
!   such 2nd run necessary to avoid a depth dependent array for id1_overlap
    Call opacity_cmf(nd, xne, temp, clfac, r, velo, dvdr, test_overlap, ntest, 2)
  End If

  Call cmf_complete(nd, xne, temp, clfac, r, velo, dvdr, rho, taur, xnh, &
    pressure, ferr, dtfcorr, opt_ray)

! remap radiation field onto coarse grid (conserve freq. integrals),
! and overwrite radiation field (from CONT, coarse grid) by remapped
! cmf-quantities
  Call interpol_cmf_to_coarse(nd)

! calculate obs. frame flux (approximate), write to file 'out_xh_obs',
! and resample onto coarse grid (flux-conservative)
  If (unasol) Then
    Call formal_obs(r, velo(1), rmax)
    Call interpol_hcmf_to_hcoarse
  End If

  Return
End Subroutine

!-----------------------------------------------------------------------

Subroutine version_cmf_all_sub(ver_cmf_all1)

  Use :: version_cmf_all, Only: ver_cmf_all

  Implicit None

  Character (10) :: ver_cmf_all1

  ver_cmf_all1 = ver_cmf_all

End Subroutine

!-----------------------------------------------------------------------

Subroutine optcmf_restart

! provides info for cmf_all

  Use :: nlte_type
  Use :: nlte_dim, Only: id_frec1, id_ndept

  Use :: nlte_opt, Only: optcmf_full
  Use :: nlte_var, Only: modnam, almost_converged, abund_have_changed, ifre, &
    fre, xj_save, kcmf_start, kcmf_end, restart_cmfall

  Use :: nlte_app, Only: met_imin, met_imax, occng
  Use :: fastwind_params, Only: wavblue, wavred
  Implicit None

! .. parameters ..
  Integer (i4b), Parameter :: nd = id_ndept
  Integer (i4b), Parameter :: ifretot = id_frec1

  Integer (i4b) :: ifre1, i, j, k
  Real (dp), Dimension (ifretot) :: fre1
  Real (dp), Dimension (nd, ifretot) :: dummy

  Real (dp) :: eblue, ered, q, q1

  Logical :: newgrid

  If (.Not. optcmf_full) Then
    Write (999, *) ' stop: optcmf_full = F and optcmf_restart called'
    Stop ' optcmf_full = F and optcmf_restart called'
  End If

  Open (1, File=trim(modnam)//'/OCCNG', Status='old', Form='unformatted')
  Rewind (1)
  Read (1) met_imin, met_imax, occng
  Read (1, End=100) almost_converged
  Go To 110

100 Print *, ' Restart: ALMOST_CONVERGED not defined (left at default)'

110 Close (1)


! JO Oct. 25: abund_have_changed policy changed.
! At this point, abund_have_changed should be always true.
! ALMOST_CONVERGED controlled elsewhere
! if(abund_have_changed) almost_converged=.false.
! old version:
! if(abund_have_changed) almost_converged=.false.
! new version:
  If (abund_have_changed) Then
    Write (999, *) ' stop: At this point (optcmf_restart), ABUND_HAVE_CHANGED &
      &should never be .TRUE.!'
    Stop ' At this point (optcmf_restart), ABUND_HAVE_CHANGED should &
      &never be .TRUE.!'
  End If

  Write (999, *) ' Restart: file OCCNG successfully read!'
  Write (999, *)
  Print *, ' Restart: file OCCNG successfully read!'
  Print *

! JO: added April 2016
  Open (1, File=trim(modnam)//'/CONT_FORMAL_ALL', Status='unknown', &
    Form='unformatted')

  Rewind 1

! we read until XJ (first entries after FRE1 are OPAC, THOMSON, STRUE, XJ)
  Read (1) ifre1, (fre1(i), i=1, ifretot), ((dummy(i,j),i=1,nd), j=1, ifretot) &
    , ((dummy(i,j),i=1,nd), j=1, ifretot), ((dummy(i,j),i=1,nd), j=1, ifretot) &
    , ((dummy(i,j),i=1,nd), j=1, ifretot)
  Close (1)

  If (ifre1/=ifre) Then
    Write (999, *) ' old:', ifre1, ' new:', ifre
    Write (999, *) &
      ' WARNING !!!! FREQUENCY GRID HAS CHANGED in optcmf_restart'
    Print *, ' old:', ifre1, ' new:', ifre
    Print *, ' WARNING !!!! FREQUENCY GRID HAS CHANGED in optcmf_restart'
    newgrid = .True.
  Else
    newgrid = .False.
  End If

  If (.Not. newgrid) Then
    Do i = 1, ifre
      If (abs(1.D0-fre1(i)/fre(i))>1.D-12) Then
        Write (999, *) ' stop: fre1 ne fre in optcmf_restart'
        Stop ' fre1 ne fre in optcmf_restart'
      End If
    End Do
!   if (dummy(1,ifre+1).ne.2) stop ' XJ qualifier ne 2 in optcmf_restart'
  Else
    If (abs(1.D0-fre1(1)/fre(1))>1.D-12) Then
      Write (999, *) ' stop: fre1(1) ne fre(1) in optcmf_restart'
      Stop ' fre1(1) ne fre(1) in optcmf_restart'
    End If

    If (abs(1.D0-fre1(ifre1)/fre(ifre))>1.D-12) Then
      Write (999, *) ' stop: fre1(ifre1) ne fre(ifre) in optcmf_restart'
      Stop ' fre1(ifre1) ne fre(ifre) in optcmf_restart'
    End If
    If (dummy(1,ifre1+1)/=2) Then
      Write (999, *) ' stop: newgrid and XJ qualifier ne 2 in optcmf_restart'
      Stop ' newgrid and XJ qualifier ne 2 in optcmf_restart'
    End If
  End If

  Allocate (xj_save(nd,ifre))
  If (.Not. newgrid) Then
    Do i = 1, ifre
      xj_save(:, i) = dummy(:, i)
!     for tests
!     print*,1.d8/fre(i),xj_save(1,i),xj(1,i)
    End Do
    Print *, &
      ' Restart: XJ_SAVE (smoothed) successfully read from CONT_FORMAL_ALL!'
    Print *
  Else !                            interpolation of old results onto new grid
    xj_save(:, 1) = dummy(:, 1)
    xj_save(:, ifre) = dummy(:, ifre1)
outer: Do i = 2, ifre - 1
      Do k = 1, ifre1 - 1
        If (fre(i)>fre1(k) .And. fre(i)<=fre1(k+1)) Then
          q = log10(fre(i)/fre1(k))/log10(fre1(k+1)/fre1(k))
          q1 = 1.D0 - q
          xj_save(:, i) = q1*log10(dummy(:,k)) + q*log10(dummy(:,k+1))
          xj_save(:, i) = 10.D0**xj_save(:, i)
          Cycle outer
        End If
      End Do
      Write (999, *) ' stop: fre not found in fre1 (subr. optcmf_restart)'
      Stop ' fre not found in fre1 (subr. optcmf_restart)'
    End Do outer
    Write (999, *) &
      ' Restart: XJ_SAVE (smoothed) interpolated from CONT_FORMAL_ALL!'
    Write (999, *)
    Print *, ' Restart: XJ_SAVE (smoothed) interpolated from CONT_FORMAL_ALL!'
    Print *
  End If

  If (newgrid) Go To 120 !          kcmf_start and kcmf_end need to be defined

  Open (1, Err=120, File=trim(modnam)//'/KCMF', Status='old', &
    Form='formatted')
  Rewind 1
  Read (1, *) kcmf_start, kcmf_end
  Close (1)
  Go To 130

! file does not exist ('old' models); use wavblue and wavred to calculate
! indices
120 eblue = 1.D8/wavblue
  ered = 1.D8/wavred

  Do i = 1, ifre
    If (fre(i)>ered) Exit
  End Do
  kcmf_start = i

  Do i = ifre, 1, -1
    If (fre(i)<eblue) Exit
  End Do
  kcmf_end = i

  If (abs(1.-fre(kcmf_start)/ered)>0.1 .Or. abs(1.-fre(kcmf_end)/eblue)>0.1) &
    Then
    Write (999, *) &
      ' stop: Problems with kcmf_start or kcmf_end in optcmf_restart'
    Stop ' Problems with kcmf_start or kcmf_end in optcmf_restart'
  End If

130 Continue

  Write (999, *) ' Restart: KCMF_START/END = ', kcmf_start, kcmf_end
  Write (999, *) ' corresponding wavelengths:', 1.D8/fre(kcmf_start), &
    1.D8/fre(kcmf_end)
  Write (999, *)
  Print *, ' Restart: KCMF_START/END = ', kcmf_start, kcmf_end
  Print *, ' corresponding wavelengths:', 1.D8/fre(kcmf_start), &
    1.D8/fre(kcmf_end)
  Print *

  restart_cmfall = .True.
  Print *, ' RESTART_CMFALL set to .true.'
  Print *

  Return
End Subroutine

!-----------------------------------------------------------------------

Subroutine line_list(test_overlap)

! merge (explicit and background elements) and modify the original line-lists

! full coupling between wavblue and wavred
! beyond, cmf transfer with approx. background radiation field
! Thus line list:
! lam < wavblue:   only cmf lines from explicit elements
! wavblue < lam < wavred :   important lines from background,
! plus cmf lines from explicit elements
! wavred < lam < xlam(ired): only cmf lines from explicit elements

! philosophy regarding tests for indexcmf might need to be changed for
! optmixed-approach

! JO -- IMPORTANT NOTE: in current version, also depacked data from
! LINES_xxx.dat included

! NOTE regarding bg-elements:
! please remember that any manipulation here (id1, xlam1, gf1) will affect
! 'only' the lines (opacities and emissivities) as treated in the
! cmf-transport, particularly the resulting J, Jbar and ALO.
! With respect to rates (as calculated in calc_tlu_met),
! the contributing components and data (particularly gf) are taken from
! the original files (id,xlam,gf) and stored in met_rbb. This is done
! in subr. fgrid_cmf. To manipulate any data consistently, you have to
! hack id_rbb created below (and maybe some data in subr. fgrid_cmf)

! 'CaVI-problem' (March 2016)
! During our tests, it turned out that for hot models the convergence of
! the T-struct (and rest) is hampered by the interaction of the OV 1-3 line
! at 629.73 with the Ca VI line 1-4 at 629.49 (this line has three strong
! components, at 629.49, 634.058 and 641.66), Since the cooling/heating
! balance is quite strongly influenced by the CBB transitions from OV 1-3
! (30% of the total heating rate for model d2v), this transition is very
! important. What happens is the following: when Ca VI overlaps with OV,
! the cooling becomes strong because the upper level 3 is no longer as
! populated. Below a certain temperature, Ca VI looses its strength,
! OV can pump, and the third level can heat. Then the cycle starts again
! (more CaVI, less pumping, more cooling, etc.)
! Thus, for specific conditions we exclude the 3rd component of CA VI
! at 629.49 from the line list, and the temp. convergence is no longer
! hampered.

  Use :: nlte_type
  Use :: nlte_dim
  Use :: fund_const, Only: clight
  Use :: fastwind_params, Only: names1, natom, wavblue, wavred, gf_min, &
    deltavmax, fpath_ufunc

  Use :: princesa_var, Only: gl, zl, indtr, data, index1, indat1, le, li, &
    labl, ifirsl, nat, labat

  Use :: ffr_error

  Use :: nlte_var, Only: indxlamc, xlamcmf, indexrbb, mmlow, mmup, indexcmf, &
    restart_cmfall

  Use :: tcorr_var, Only: enatcor, temp_converged

  Use :: nlte_app, Only: jatom, calcion, indexel, metall, jatom_full, indrec, &
    irbb, ilin, met_rbb, indexel, teff

  Use :: nlte_lines, Only: id, gf, xlam, iblue, ired, ntotnl3

  Use :: cmf_all_var, Only: indexrbb_inv, ntotnl_all, id1, gf1, xlam1, ntot, &
    ntot_expl_inside, elements_cmf, indexcmf1, id_rbb, indexel_inv, &
    opal_ntest, index_id1, id1_overlap, lineno, depacked, xlamcmf_depacked, &
    indxlamc_depacked, sumgf_depacked

  Implicit None


  Real (dp) :: xxlam, lam, flu, xlamold, xl, deltav, flu1, sumf, maxf
  Integer (i4b) :: ii, indi, j, k, kk, l, ml, mu, irest, low, up, ic, jc, jce, &
    index, jce_icmf1, jce_red, jce_outside, n, jj, ji, igenio, kold, i, jfull, &
    jfullold, iu, nocomp, ic1, ic2, icomp

  Logical :: first, newj, no_calcium_six, test_overlap

  Character (1) :: met
  Character (4) :: ret
  Character (6) :: lab, dc, level, leveu
  Character (32) :: fichero, stafil, linetab

  Data first/.True./

  Print *, 'number of lines from original nl3-files:', ired - iblue + 1

  If (id_nttrd/=indtr) Then
    Write (999, *) ' stop: something wrong with id_nttdr (subr line_list)'
    Stop ' something wrong with id_nttdr (subr line_list)'
  End If

  If (first) Then

    lineno = 0

!   -----------------------------------------------------------
!   prepare depacking, but only if non-HHe elements are present
    If (maxval(indexel(3:natom))>-1) Then

!     use LINES_xxx.dat (as encoded in linetab) to depack specific transitions
!     from explicit elements
      Open (1, File='ATOM_FILE', Status='old')
      Rewind 1
      Read (1, Fmt='(a)') fichero
      Read (1, Fmt='(a)') stafil
      Read (1, Fmt='(a)') linetab
      Close (1)

      iu = 337
!     iou = 338
      dc = ':T'

      Print *
      Print *, 'file ', trim(linetab), &
        ' used for depacking transitions from explicit elements'

      Open (Unit=iu, File=fpath_ufunc//linetab, Form='formatted', &
        Status='old')
!     open (unit=iou,file='control_lines.dat',status='unknown')

!     rewind!!!

      iu = -337
      Call ffracc(iu, ret)
      If (on_error(ret)) Go To 130

      Call transic_count(lineno, ret)
      Write (999, *) lineno, ' (de-packed) lines present in file'
      Print *, lineno, ' (de-packed) lines present in file'

      If (on_error(ret)) Go To 140

      Allocate (depacked(lineno))

      Rewind (abs(iu))
      Call transic_read(lineno, ret) ! wavelengths (lambda > 2000 A) converted
      If (on_error(ret)) Go To 140

!     for tests
!     do i=1, lineno
!     print*,i,depacked(i)
!     enddo

      Write (999, Fmt='(a)') ' data successfully read! '
      Write (999, Fmt='(a)') ' consistency of statistical weights checked '
      Write (999, Fmt='(a)') &
        ' wavelengths (lambda > 2000 A) converted from air to vacuum'
      Write (999, *)
      Write (*, Fmt='(a)') ' data successfully read! '
      Write (*, Fmt='(a)') ' consistency of statistical weights checked '
      Write (*, Fmt='(a)') &
        ' wavelengths (lambda > 2000 A) converted from air to vacuum'
      Print *
      Close (iu)
!     close (iou)
    Else
!     JP contemporary hack, to ensure that "depacked" is allocated even for
!     HHe as explicit elements
      If (lineno/=0) Then
        Write (999, *) ' stop: HHe and lineno(depacked) ne 0'
        Stop ' HHe and lineno(depacked) ne 0'
      End If
      Allocate (depacked(lineno))
    End If !                        prepare depacking

!   -----------------------------------------------------------

    If (wavblue<xlam(iblue)) Then
!     wavblue=xlam(iblue)+20.
!     changed by JO Aug. 2013; note cooler stars with lines
!     only until 228 A (Teff < 25kK) or even 250 A (Teff < 10kK)
      wavblue = xlam(iblue) + 2.
      Print *, ' wavblue (cmf_all) reset to ', wavblue
    End If

    If (wavred>xlam(ired)) Then
      Write (999, *) ' stop: wavred (cmf_all) > wavred(nlte_approx)'
      Stop ' wavred (cmf_all) > wavred(nlte_approx)'
    End If

    Do j = 1, id_nttrd
      ii = indxlamc(j)
      Do k = 1, index1
        If (indexrbb(k)==ii) Then
          indexrbb_inv(j) = k
          Go To 100
        End If
      End Do
      Write (999, *) ' stop: ii not found in indexrbb (subr. line_list)'
      Stop 'ii not found in indexrbb (subr. line_list)'

100   Continue
    End Do

    ntotnl_all = ntotnl3 + id_nttrd

!   this has to be done only once, since dimension fix
    Allocate (id1(ntotnl_all), gf1(ntotnl_all), xlam1(ntotnl_all))
    Allocate (id1_overlap(ntotnl_all))
    If (test_overlap) Allocate (opal_ntest(ntotnl_all), index_id1(ntotnl_all))

!   index_array providing index of explict atom w.r.t. background elements
!   shifted from subr. profile_cmf
    If (id_atoms/=nat) Then
      Write (999, *) ' stop: id_atoms ne nat (subr. line_list)'
      Stop ' id_atoms ne nat (subr. line_list)'
    End If

    Do k = 1, nat
      Do i = 1, natom
        If (labat(k)==names1(i)) Then
          indexel_inv(k) = i
        End If
      End Do
!     test
      If (indexel(indexel_inv(k))/=k) Then
        Write (999, *) ' stop: something wrong with index_el(inv), &
          &H or He not explicit? (subr. line_list)'
        Stop ' something wrong with index_el(inv), H or He not &
          &explicit? (subr. line_list)'
      End If
    End Do

!   determine the index of met_rbb related to bg line-list within wavblue,
!   wavred.
!   Only for selected NLTE bg-elements (jatom_full=1) and lines with defined
!   upper level,
!   but for ALL ionization stages. Thus, this task needs to be performed only
!   once.
    Allocate (id_rbb(ntotnl3))
    id_rbb = 0

    Write (999, *)
    Write (999, *) ' significant difference (>20%) in wavelength &
      &(bg-elements, line list vs. packed transition)!'
    Print *
    Print *, ' significant difference (>20%) in wavelength (bg-elements, &
      &line list vs. packed transition)!'

lines_rbb: Do l = iblue, ired
      If (xlam(l)<wavblue) Cycle lines_rbb
      If (xlam(l)>wavred) Exit lines_rbb

      k = id(l)/1000000
      If (jatom_full(k)==0) Cycle lines_rbb

      irest = id(l) - k*1000000
      j = irest/100000
      irest = irest - j*100000
      low = irest/100
      up = irest - low*100
      If (up==0) Cycle lines_rbb !  .and. gf(l).lt.gf_min) ! if we want to cut
!     the list

      n = indrec(k, j)
      jj = irbb(n) - 1

      Do ii = 1, ilin(n)
        ji = jj + ii
        If (low==met_rbb(ji)%low .And. up==met_rbb(ji)%lup) Exit
      End Do
      If (abs(1.-met_rbb(ji)%wave/xlam(l))>0.20) Then
        Write (999, Fmt='(4(2x,i2),2(2x,f10.3))') k, j, low, up, xlam(l), &
          met_rbb(ji)%wave
        Write (*, Fmt='(4(2x,i2),2(2x,f10.3))') k, j, low, up, xlam(l), &
          met_rbb(ji)%wave
      End If

      id_rbb(l) = ji
!     hack to exclude certain components
!     if(id(l).eq.8300110 .and. xlam(l).ne.374.073) id_rbb(l)=0
!     if(id(l).eq.8300106 .and. xlam(l).ne.702.337 .and. xlam(l).ne.703.855)
!     id_rbb(l)=0
!     if(id(l).eq.8300105 .and. xlam(l).ne.835.289 .and. xlam(l).ne.832.929)
!     id_rbb(l)=0

    End Do lines_rbb
    Write (999, *)
    Print *

    first = .False.
  End If

! JO July 2018 included to update indexrbb_inv in case
  If (restart_cmfall) Then
    Do j = 1, id_nttrd
      ii = indxlamc(j)
      Do k = 1, index1
        If (indexrbb(k)==ii) Then
          indexrbb_inv(j) = k
          Go To 110
        End If
      End Do
      Write (999, *) ' stop: ii not found in indexrbb (subr. line_list)'
      Stop 'ii not found in indexrbb (subr. line_list)'

110   Continue
    End Do
    Print *
    Print *, 'indexrbb_inv updated'
    Print *
  End If

! local copy of indexcmf:
! lines treated in cmf_complete will obtain an indexcmf1=3 or higher (for
! depacked levels)
! lines outside range will keep their previous status

  indexcmf1 = indexcmf

! initialize
  depacked%index = 0

  ic = 0 !                          counter for new list
  jc = 0 !                          counter refering to explicit elements
  jce = 0 !                         counter for cmf-lines from explicit
! elements in new list
  jce_outside = 0 !                 counter for cmf-lines from explicit
! elements outside wavblue,wavred
  jce_icmf1 = 0 !                   counter for lines from expl. elements with
! index_cmf = 1 (outside imia,imaa)
  jce_red = 0 !                     counter for cmf-lines from explicit
! elements with lam > lam (ired)
  newj = .True.
  elements_cmf = 0

! do k=1,30
! write(*,5) k,(calcion(k,j),j=1,9)
! enddo
! 5 format(10(i2,2x))

! check whether problematic CA VI line should be removed
! only for hot models and when temperature is not converged
! JO: in case, also other lines might need to be removed;
! check with overlap detector at end of this routine

  no_calcium_six = .False.
  If (teff>40000. .And. enatcor .And. .Not. temp_converged) &
    no_calcium_six = .True.

! at first, remove all lines with calcion = 0, wave < wavblue, and explicit
! elements
! then include all lines from explicit elements
lines: Do l = iblue, ired
    If (xlam(l)<wavblue) Cycle lines
!   this is the hack to avoid lines overlapping with HeI sing resonance line
!   if (xlam(l).gt.580.and.xlam(l).lt.590.) cycle lines
!   this is the hack to exclude most lines...
!   if (xlam(l).lt.6000.)  cycle lines
!   print*,id(l),xlam(l)
    k = id(l)/1000000
    If (jatom(k)==0) Cycle lines

!   for tests without non-selected elements
!   if(jatom_full(k).eq.0) cycle lines

    If (indexel(k)>0) Cycle lines ! included in nlte

    irest = id(l) - k*1000000
    j = irest/100000

!   This is the FeV line (374.245) which is also responsible for the NIII
!   emission
!   note: delta v = 37 km/s, i.e., for small turbulence, there should be no
!   effect.
!   if (xlam(l).gt.372.and.xlam(l).lt.376. .and. id(l).eq.26500200) cycle
!   lines


    If (calcion(k,j)==0) Cycle lines ! too low abundance

    irest = irest - j*100000
    low = irest/100
    up = irest - low*100
!   for tests
!   if(k.eq.8 .and. up.eq.0) cycle lines

!   for tests without OIII line overlap
!   if(k.eq.8 .and. j.eq.3 .and.low.eq.1 .and. up.eq.10) then
!   if(xlam(l).ne.374.073) cycle lines
!   print*,'O3 found',xlam(l)
!   endif


!   if(k.eq.8 .and. j.eq.3 .and.low.eq.1 .and. up.eq.6) then
!   if(xlam(l).ne.702.337 .and. xlam(l).ne.703.855) cycle lines
!   print*,'O3 found',xlam(l)
!   endif
!   if(k.eq.8 .and. j.eq.3 .and.low.eq.1 .and. up.eq.5) then
!   if(xlam(l).ne.835.289 .and. xlam(l).ne.832.929) cycle lines
!   print*,'O3 found',xlam(l)
!   endif

    If (up==0 .And. gf(l)<gf_min) Cycle lines ! only stronger resonance lines

    n = indrec(k, j)

!   test overlap with N48
!   if(up.eq.0 .and. gf(l).lt.10. .and. xlam(l).eq.387.371)) cycle lines !
!   overlap with FeV  26500400

!   no check for met_imin, met_imax required, since no problem with
!   opacities missing at certain depth-points

!   if only selected elements, uncomment following line
!   if(jatom_full(k).eq.0) cycle lines

!   overlap between CaVI and OV (629.496 and 629.732) prevents
!   temperature convergence for hot models (e.g., d2v).
!   remove responsible CA VI line in CMF transport
!   JO: when data-base updated, check frequencies
!   JO: other lines might need to be removed as well
    If (no_calcium_six .And. k==20) Then
      If (j==6 .And. low==1 .And. xlam(l)>629. .And. xlam(l)<630.) Then
        Print *, 'CA VI line at', xlam(l), ' removed from cmf transport'
        Print *
        Cycle lines
      End If
    End If

    met = metall(k, j, low)

    If (low/=1 .And. met==' ' .And. jatom_full(k)==0) Cycle lines ! consistent
!   with
!   sumopal

!   exclude lines overlapping with NV 2->3
!   if (xlam(l).gt.250.and.xlam(l).lt.280.) cycle lines


    If (up/=0) Then
      met = metall(k, j, up)
      If (met==' ' .And. jatom_full(k)==0) Cycle lines ! no approx. upper
!     level
    End If

!   for tests
!   if(k.eq.15.and.j.eq.3.) print*,k,j,low,up,xlam(l),gf(l)

    elements_cmf(k) = 1 !           element k treated in cmf

!   JO changed Oct. 2016
!   remove Fe-lines (with uncertain gf values) overlapping with HeI resonance
!   line
!   JO: other lines might need to be removed as well if there are problems
!   in the outer region, in particular the resonance line of MN IV (584.296
!   A).

!   JO: when data-base updated, check frequencies
    If (xlam(l)>=584.23 .And. xlam(l)<=584.43) Then
!     FeIV 584.368/584.397 with up=0
      If (k==26 .And. j==4 .And. (low==6 .Or. low==8)) gf(l) = 1.D-5
    End If

!   for tests
!   if(up.ne.0) then
!   print*,xlam(l),jatom_full(k),id(l),metall(k,j,low),' ',metall(k,j,up)
!   else
!   print*,xlam(l),jatom_full(k),id(l),metall(k,j,low),' 0'
!   endif

!   info on line-transitions from explict elements
120 If (newj) Then
      jc = jc + 1
!     assume that there are lines from explicit elements with lam > xlam(ired)
      If (jc>id_nttrd) Then
        Write (999, *) ' stop: jc > id_nttrd in line_list'
        Stop ' jc > id_nttrd in line_list'
      End If
      ii = indxlamc(jc)
      xxlam = xlamcmf(ii)
      If (xxlam==0.) Then
        If (indexcmf(ii)/=1) Then
          Write (999, *) ' stop: xxlam = 0 and indexcmf(ii) ne 1 in line_list'
          Stop ' xxlam = 0 and indexcmf(ii) ne 1 in line_list'
        End If
        jce_icmf1 = jce_icmf1 + 1
        Go To 120 !                 from ions outside imia,imaa
      End If
      If (indexcmf(ii)/=2 .And. xxlam<xlam(ired)) Then
!       JO: changed Dec. 2016: sobo lines inside range can be present
!       if occup. numbers not defined at ALL depth points
        Print *, ' line at ', xxlam, ' treated in Sobo (outside ilow, imax)'
        jce_icmf1 = jce_icmf1 + 1
        Go To 120 !                 from ions outside imia,imaa
!       stop 'indexcmf(ii) ne 2 in line_list'
      End If
!     this is the end condition
      If (xxlam>=xlam(ired)) Exit lines
      kk = indexrbb_inv(jc)
!     kk is the index w.r.t. index1
      indi = indat1(kk)
      lam = data(indi+1)
      If (lam/=xxlam) Then
        Write (999, *) ' stop: lam ne xxlam (subr. line_list)'
        Stop ' lam ne xxlam (subr. line_list)'
      End If
      flu = data(indi+2)
      ml = mmlow(kk)
      mu = mmup(kk)
!     for tests, works
!     if(lablinelo(ii).ne.labl(ml)) stop ' error in label(low) in line_list'
!     if(lablineup(ii).ne.labl(mu)) stop ' error in label( up) in line_list'
    End If

    If (xlam(l)<xxlam .And. xlam(l)<=wavred) Then
      ic = ic + 1
      id1(ic) = id(l)
      gf1(ic) = gf(l)
      xlam1(ic) = xlam(l)
      newj = .False.
    Else If (xxlam<xlam(ired)) Then
!     lines from explicit elements with lam < wavblue included as first lines
!     lines from explicit elements with lam > wavred included as last lines

!     leading '1' to tell that this is a line from explicit elements
!     (different coding)
      index = le(ml)*1D8 + ml*1D4 + mu
      ic = ic + 1
      jce = jce + 1
      id1(ic) = index
      gf1(ic) = gl(ml)*flu
      xlam1(ic) = lam
      newj = .True.
      If (lam<wavblue .Or. lam>wavred) Then
        jce_outside = jce_outside + 1
      Else
!       considered line will be treated in cmf_complete
        Where (depacked%levlo==labl(ml) .And. depacked%levup==labl(mu))
!         note: only transitions with xlam ne 0 (inside imia,imaa)
          depacked%index = ic
        End Where
        indexcmf1(ii) = 3
      End If

!     print*,ic,l,jc,xlam1(ic-1),xlam1(ic),id1(ic)
      Go To 120

    End If

  End Do lines

! last lines from explicit elements (not used) should be sobo-lines
  Do j = jc, id_nttrd
    ii = indxlamc(j)
    If (indexcmf(ii)/=1) Then
      Write (999, *) ' stop: something wrong in philosophy -- line_list'
      Stop ' something wrong in philosophy -- line_list'
    End If
  End Do
  jce_red = id_nttrd - jc + 1

  ntot = ic
  ntot_expl_inside = jce - jce_outside
  Write (999, *) ' number of lines in merged files: ', ntot
  Write (999, *) ' including ', jce, ' cmf lines from explicit elements'
  Write (999, *) ' number of cmf lines (expl. elem.) inside wavblue, wavred:', &
    ntot_expl_inside
  Write (999, *) ' number of lines (expl. elem.) outside imia,imaa:', &
    jce_icmf1
  Write (999, *) &
    ' number of Sobolev lines (expl. elem.) with lam > lam(ired):', jce_red
  Write (999, *)
  Print *, ' number of lines in merged files: ', ntot
  Print *, ' including ', jce, ' cmf lines from explicit elements'
  Print *, ' number of cmf lines (expl. elem.) inside wavblue, wavred:', &
    ntot_expl_inside
  Print *, ' number of lines (expl. elem.) outside imia,imaa:', jce_icmf1
  Print *, ' number of Sobolev lines (expl. elem.) with lam > lam(ired):', &
    jce_red
  Print *
  If (jce+jce_icmf1+jce_red/=id_nttrd) Then
    Write (999, *) &
      ' stop: not all lines from expl. elements found in subr. line_list'
    Stop ' not all lines from expl. elements found in subr. line_list'
  End If

! test that all explicit lines have been found
  ic = 0
  Do j = 1, id_nttrd
    ii = indxlamc(j)
!   for specific tests
!   print*,j,xlamcmf(ii),indexcmf1(ii)
    If (indexcmf1(ii)==3) ic = ic + 1
  End Do

  If (ic/=ntot_expl_inside) Then
    Write (999, *) ' stop: not all explicit lines inside wavblue, &
      &wavred identified in indexcmf1'
    Stop &
      ' not all explicit lines inside wavblue, wavred identified in indexcmf1'
  End If

! for tests
! do j=1,ntot
! if(xlam1(j).gt.600. .and. xlam1(j).lt.650.) &
! &   print*,j,xlam1(j),id1(j),gf1(j)
! enddo
! stop

! exchange/add depacked components in original data
! works also for lineno = 0
  ic = ntot
  ic1 = 0
  ic2 = 0

! initialize
  xlamcmf_depacked = xlamcmf
  sumgf_depacked = 0.
  depacked%nocomp = 0

  If (lineno>0) Then
    Do j = 1, id_nttrd
      ii = indxlamc(j)
      xxlam = xlamcmf(ii)
      If (xxlam==0.) Cycle !        outside imia,imaa
      kk = indexrbb_inv(j)
!     kk is the index w.r.t. index1
      ml = mmlow(kk)
      lab = labat(le(ml))
      If (lab=='H' .Or. lab=='HE') Cycle ! only non-HHe elements
      If (indexcmf1(ii)/=3) Cycle ! only lines with cmf_all treatment
      mu = mmup(kk)
      level = labl(ml)
      leveu = labl(mu)
      indi = indat1(kk)
      lam = data(indi+1)
      If (lam/=xxlam) Then
        Write (999, *) lam, xxlam
        Write (999, *) ' stop: lam ne xxlam (subr. line_list)'
        Print *, lam, xxlam
        Stop ' lam ne xxlam (subr. line_list)'
      End If
      flu = data(indi+2)
!     print*,lam,level,' ',leveu,flu
!     rewind(iou)

!     check consistency of f-values and provide number of de-packed components
!     thus far, nocomp not required
      sumf = 0.D0
      nocomp = 0
!     note: where statement requires additional array(s), e.g.,
!     idum=0
!     where (depacked/levlo .eq. ...
!     idum=depacked%lineno
!     depacked%nocomp=maxval(idum)
!     endwhere
      maxf = 0.D0
      icomp = 0
      Do i = 1, lineno
        If (depacked(i)%levlo==level .And. depacked(i)%levup==leveu) Then
          nocomp = nocomp + 1
          flu1 = depacked(i)%flu
          sumf = sumf + flu1
!         determine largest f-value of depacked transition
!         NOTE: using amax1 destroys precision
          maxf = max(maxf, flu1)
          If (maxf==flu1) icomp = depacked(i)%lineno
        End If
      End Do

      Where (depacked%levlo==level .And. depacked%levup==leveu)
        depacked%nocomp = nocomp
      End Where

      Do i = 1, lineno
        If (depacked(i)%levlo==level .And. depacked(i)%levup==leveu) Then
          If (icomp==0) Then
            Write (999, *) ' stop:  error in icomp'
            Stop ' error in icomp'
          End If
!         print*,'found
!         ',level,leveu,depacked(i)%lineno,depacked(i)%wave,depacked(i)%nocomp,depacked(i)%index
!         output only once
          If (abs(1.-sumf/flu)>0.2 .And. icomp==depacked(i)%lineno) Then
            Write (999, *) ' inconsistent f-values for ', level, ' ', leveu
            Write (999, *) ' sum(comp):', sumf, ' packed:', flu
            Print *, ' inconsistent f-values for ', level, ' ', leveu
            Print *, ' sum(comp):', sumf, ' packed:', flu
            If (level=='O3_28' .And. leveu=='O3_39') Then
!             JO Sept. 2023: f-value noew corrected in atomic data file(s),
!             should no longer happen
              Write (999, *) 'erroneous f-value (''packed'') in &
                &WM-basic database, needs to be updated'
              Print *, 'erroneous f-value (''packed'') in WM-basic &
                &database, needs to be updated'
            Else
              Write (999, *) ' check components and flu values in ', &
                trim(linetab)
              Print *, 'check components and flu values in ', trim(linetab)
              Write (999, *) ' stop: inconsistent f-values for &
                &transition from explicit element'
              Stop &
                ' inconsistent f-values for transition from explicit element'
            End If
          End If
          If (depacked(i)%lineno==1) Then
!           exchange wavelength and gf value
            k = depacked(i)%index
            xlam1(k) = depacked(i)%wave
            gf1(k) = depacked(i)%flu*gl(ml) ! gf value ; remember:
!           gf(individual) = flu(input) *
!           gl(tot)
            ic1 = ic1 + 1
          Else
!           add info for additional components
            ic = ic + 1
            ic2 = ic2 + 1
            If (ic>ntotnl_all) Then
              Write (999, *) ' stop: ic too large'
              Stop ' ic too large'
            End If
            k = depacked(i)%index
            id1(ic) = id1(k)
            xlam1(ic) = depacked(i)%wave
            gf1(ic) = depacked(i)%flu*gl(ml) ! gf value
          End If

          If (icomp==depacked(i)%lineno) Then
!           change wavelength in xlamcmf_depacked corresponding to wavelength
!           of
!           transition with largest f-value
            xlamcmf_depacked(ii) = depacked(i)%wave
!           if different from zero, indicates that transition depacked/changed
            sumgf_depacked(ii) = sumf*gl(ml)
!           overwrite indexcmf1 with number of components (when different from
!           1):
!           indexcmf1 = 4 means two components, = 5 means three components,
!           and so on
            If (indexcmf1(ii)/=3) Then
              Write (999, *) ' stop: something rotten with indexcmf1'
              Stop ' something rotten with indexcmf1'
            End If
            indexcmf1(ii) = 2 + depacked(i)%nocomp ! works also if only one
!           component
          End If

        End If
      End Do
    End Do

!   for tests
!   do i=1, lineno
!   print*,'testdp',i,depacked(i)
!   enddo



!   create new indexarray for xlamcmf_depacked
    Call indexx(id_nttrd, xlamcmf_depacked, indxlamc_depacked)

    If (ic2/=ic-ntot) Then
      Write (999, *) ' stop: error in exchange/add of components'
      Stop ' error in exchange/add of components'
    End If

    Write (999, *)
    Write (999, *) ic1, ' lines/components exchanged (from ', trim(linetab), &
      ')'
    Write (999, *) ' components added (from ', trim(linetab), ')'
    Write (999, *)
    Print *
    Print *, ic1, ' lines/components exchanged (from ', trim(linetab), ')'
    Print *, ic2, ' components added (from ', trim(linetab), ')'
    Print *

    ntot = ic
    jce = jce + ic2
    ntot_expl_inside = ntot_expl_inside + ic2
    Write (999, *) ' FINAL number of lines in merged files: ', ntot
    Write (999, *) ' including ', jce, ' cmf lines from explicit elements'
    Write (999, *) &
      ' number of cmf lines (expl. elem.) inside wavblue, wavred:', &
      ntot_expl_inside
    Print *, ' FINAL number of lines in merged files: ', ntot
    Print *, ' including ', jce, ' cmf lines from explicit elements'
    Print *, ' number of cmf lines (expl. elem.) inside wavblue, wavred:', &
      ntot_expl_inside

!   final test that %lineno correct (for all used lines, i.e., with %nocomp
!   >0)
    Do i = 1, lineno
      If (depacked(i)%nocomp>0 .And. depacked(i)%lineno>depacked(i)%nocomp) &
        Then
        Write (999, *) depacked(i)%levlo, depacked(i)%levup, ' : ', &
          depacked(i)%lineno
        Write (999, *) ' number of components: ', depacked(i)%nocomp
        Write (999, *) ' stop: something wrong with %lineno in &
          &LINES_xxx.dat file, please correct'
        Print *, depacked(i)%levlo, depacked(i)%levup, ' : ', &
          depacked(i)%lineno
        Print *, ' number of components: ', depacked(i)%nocomp
        Stop ' something wrong with %lineno in LINES_xxx.dat file, &
          &please correct'
      End If
    End Do

!   re-sort line list
    If (size(xlam1)/=ntotnl_all) Then
      Write (999, *) ' stop: error in nsize (line_list)'
      Stop ' error in nsize (line_list)'
    End If

!   sort according to wavelength
    Call sort_line_list(ntot, xlam1(1:ntot), gf1(1:ntot), id1(1:ntot))

  Else If (lineno==0.) Then

    indxlamc_depacked = indxlamc

  Else
    Write (999, *) ' stop: lineno <0!!!'
    Stop ' lineno <0!!!'

  End If

! for tests
! do i=1,ic2
! j=ntot-ic2+i
! print*,id1(j),xlam1(j),gf1(j)
! enddo

! for tests
! do i=1,ntot
! if(xlam1(i).gt.372. .and. xlam1(i).lt.376.)
! print*,'dpack',xlam1(i),gf1(i),id1(i)
! enddo

! do j=1,id_nttrd
! ii=indxlamc_depacked(j)
! i=indxlamc(j)
! print*,j,xlamcmf_depacked(ii),indexcmf1(ii),xlamcmf(i)
! enddo

  Return
! in case, return here

! overlap detector
! detects overlap of important (resonance lines), within prescribed deltav
! E.g., it detects overlap between CaVI and OV (629.496 and 692.732) that
! prevents
! T-conv for model d2v (with 'standard' Delta n = 5 and gfmin = 0.1).
! If Delta n = 10 and gfmin = 0.05, it also detects OIII/NIII overlap
! responsible for OIII emission at late O-types.
! JO: needs to be tested when N/O treated as explicit elements
  kold = 0
  Do ii = 1, ntot
    index = id1(ii)
    If (index>=1D8) Then
!     explicit element
      kk = index/1D8
      irest = index - kk*1D8
      ml = irest/10000
      mu = irest - ml*10000
      ji = igenio(le(ml), li(ml)) ! general counting
      k = indexel_inv(kk)
      j = int(zl(ml)+1.D-6) + 1 !   in astro-definition (1 = neutral)
!     standard DELTA n and gfmin
      If (ml==ifirsl(ji) .And. mu-ml<=5 .And. gf1(ii)>0.1) Then ! important
!       resonance
!       line
        xl = xlam1(ii)
        If (kold==0) Then
          Write (*, Fmt='("explicit   ",4i3,2x,f10.4)') k, j, ml, mu, xl
          xlamold = xl
          kold = k
        Else
          deltav = (xl-xlamold)/xl*clight
          If (k/=kold .And. deltav<deltavmax) Then ! always
            Write (*, Fmt='("explicit overlap  ",4i3,2x,f10.4)') k, j, ml, mu, &
              xl
          Else
            Write (*, Fmt='("explicit   ",4i3,2x,f10.4)') k, j, ml, mu, xl
          End If
          xlamold = xl
          kold = k
        End If
        jfullold = 1
      End If
    Else
      k = index/1000000
      irest = index - k*1000000
      j = irest/100000
      irest = irest - j*100000
      low = irest/100
      up = irest - low*100
!     standard DELTA n and gfmin
      If (low==1 .And. up/=0 .And. up-low<=5 .And. gf1(ii)>0.1) Then ! important
!       resonance
!       line
        xl = xlam1(ii)
        If (kold==0) Then
          Write (*, Fmt='("background ",4i3,2x,f10.4)') k, j, low, up, xl
          xlamold = xl
          kold = k
          jfullold = jatom_full(k)
        Else
          deltav = (xl-xlamold)/xl*clight
          jfull = jatom_full(k)
          If (k/=kold .And. (jfullold==1 .Or. jfull==1) .And. &
            deltav<deltavmax) Then ! only for explict and selected elements
            Write (*, Fmt='("background overlap",4i3,2x,f10.4)') k, j, low, &
              up, xl
          Else
            Write (*, Fmt='("background ",4i3,2x,f10.4)') k, j, low, up, xl
          End If
          xlamold = xl
          kold = k
          jfullold = jfull
        End If
      End If
    End If
  End Do


! error in input unit exit

130 Continue
  Write (*, Fmt='(a)') ' error in input unit number '
  Stop

! error conditions
140 Select Case (ret)
  Case ('RET1') !                   end of file "no esperado" exit
    Write (*, Fmt='(a)') ' eof not expected. there is something wrong! '
    Stop
  Case ('RET2') !                   error in ffr subroutines exit
    Write (*, Fmt='(a)') ' error in ffr subroutines '
    Stop
  Case Default
    Stop ' wrong error condition in transic'
  End Select

End Subroutine

!-----------------------------------------------------------------------

Subroutine fgrid_cmf

! grid for coupled cmf-transfer; equidistant in dlam/lam = vref/clight
! vref=vturb/3, with vturb = max(vturb,vturbmin)
! i.e., if vturb = 0, artificial turb. of vturbmin (5 km/s) is used,
! to avoid too large freq. grid

! after some tests, we decided to use vturb to obtain a reasonable
! resolution for all elements at all depth points.

! To safe time, one might also use vth(O), but then the resolution of
! lines from heavy elements is below optimum

! wavelength shifted to closest grid points.

! Note: central line frequencies can be shifted by +/- 0.5 dlambda,
! corresponding to vref/2 = vturb/6
! For max(vturb) approx. 15 km/s, this corresponds to +/- 2.5 km/s.
! This uncertainty is reflected in any formal integral solution based
! on the freq. grid provided here.

! the scaling of the number of freq. points, nftot, is via
! 1/ln(1+vref/c) approx 1/vref

  Use :: nlte_type
  Use :: nlte_dim
  Use :: fund_const, Only: clight
  Use :: fastwind_params, Only: natom, wavblue, wavred, gf_min, vturbmin

  Use :: nlte_app, Only: met_rbb, jatom_full, jmax, indrec, irbb, ilin, &
    tlu_met, tul_met
  Use :: nlte_var, Only: vturb, fre, ifre, lwion1
  Use :: nlte_lines, Only: gf, xlam, ntotnl3, unpacked_index_cmf, unpacked_gf, &
    unpacked_xlam
  Use :: cmf_all_var, Only: xlam1, ntot, vturb_cmf, vref, lam_cmf, &
    index_lam_cmf, indfre_cmf, indfre_staggered_cmf, nftot, icfirst, iclast, &
    id1, id_rbb, lineno, depacked

  Implicit None

  Integer (i2b) :: max_no_comp
  
  Integer (i4b) :: i, ic, in, ic0, din, din1, lastindex, dinmax, j, istart, &
    ii, ind, inds, nftotnew, k, n, jj, tot_comp, ji, iqual0, &
    iqual1, iqual2, irun, no_comp, index_comp, nsize

  Real (dp) :: vend, dlp, dlm, lmin, lmax, lam, energy, frec, diffe, diffe1, &
    sum_gf, cln, wave, lamax

  Real (dp), Save :: lammin, lammax

  Logical :: first

  Data first/.True./

! --------------------------------------------------------------------

! set up freq. grid and various index files only once
  If (first) Then

    If (nftot/=0) Then
      Write (999, *) ' stop: first and nftot ne 0 in fgrid_cmf'
      Stop ' first and nftot ne 0 in fgrid_cmf'
    End If

!   vturbmin in fastwind_params, typically 5 km/s
    If (vturb<vturbmin) Then
      Write (999, *) ' WARNING!!!!WARNING!!!!WARNING!!!! artificial vturb = ', &
        vturbmin/1.D5, 'km/s used in cmf_all'
      Print *, ' WARNING!!!!WARNING!!!!WARNING!!!! artificial vturb = ', &
        vturbmin/1.D5, 'km/s used in cmf_all'
      vturb_cmf = vturbmin
      vref = vturb_cmf
    Else
      vturb_cmf = vturb
      vref = vturb_cmf
    End If

    vref = vref/3. !                3 points per vth is sufficient
!   note: CMFSING works with NF = 21 points for 6-8 vth
    Write (999, *)
    Write (999, *) ' reference velocity used for freq. grid = ', vref/1.D5, &
      ' km/s'
    Write (999, *) ' max. error due to shift of rest-wavel. is +/- ', &
      vref/2.D5, ' km/s,'
    Write (999, *) ' corresponding to dlambda = ', vref/2./clight*4500, &
      ' A at 4500 A'
    Print *
    Print *, ' reference velocity used for freq. grid = ', vref/1.D5, ' km/s'
    Print *, ' max. error due to shift of rest-wavel. is +/- ', vref/2.D5, &
      ' km/s,'
    Print *, ' corresponding to dlambda = ', vref/2./clight*4500, &
      ' A at 4500 A'

    vend = sqrt(20.D5**2+vturb_cmf**2)*5
!   typical value for vth(H) -- on blue side, even factor 2 smaller - no H,
!   but He --, &
!   assuming strong line 5 vth wide
    lammin = wavblue*(1.-vend/clight)
    lammax = wavred*(1.+vend/clight)

!   verify that nftot does not change
    nftotnew = log(lammax/lammin)/log(1.D0+vref/clight) + 1

    nftot = nftotnew

!   allocate only once, since nftot should not change during iteration
    Allocate (lam_cmf(nftot), index_lam_cmf(nftot), indfre_cmf(nftot), &
      indfre_staggered_cmf(nftot))

    cln = (1.D0+vref/clight)
    Do i = 1, nftot
      lam_cmf(i) = lammin*cln**(i-1)
    End Do

!   define index-file w.r.t. continuum frequency grid

    energy = 1.D8/lam_cmf(1) !      lam_cmf from blue to red, fre vice vera
    If (energy>fre(ifre)) Then
      Write (999, *) ' stop: lam_cmf(blue) outside fre! (subr. fgrid_cmf)'
      Stop ' lam_cmf(blue) outside fre! (subr. fgrid_cmf)'
    End If

    energy = 1.D8/lam_cmf(nftot) !  lam_cmf from blue to red, fre vice vera
    If (energy<fre(1)) Then
      Write (999, *) ' stop: lam_cmf(red) outside fre! (subr. fgrid_cmf)'
      Stop ' lam_cmf(red) outside fre! (subr. fgrid_cmf)'
    End If

!   indfre_staggered_cmf: such that 0.5(fre(ind+1)+fre(ind)) > lam >
!   0.5(fre(ind)+fre(ind-1))
!   indfre_cmf: such that fre(ind) > energy(lam) > fre(ind-1), for
!   interpolation
    istart = 2
    Do ii = nftot, 1, -1
      lam = lam_cmf(ii)
      energy = 1.D8/lam
freloop: Do i = istart, ifre
        frec = fre(i)
        If (energy<=frec) Then
          ind = i
          diffe = abs(energy-frec)
          diffe1 = abs(energy-fre(i-1))
          If (diffe<=diffe1) Then
            inds = i
          Else
            inds = i - 1
          End If
          istart = i
          Exit freloop
        End If
      End Do freloop
!     print*,ii,i,energy,frec,fre(i-1),ind,istart
      indfre_staggered_cmf(ii) = inds
      indfre_cmf(ii) = ind
    End Do


!   provide info on individual line-components for selected bg-elements
    met_rbb%no_comp = 0
    met_rbb%sum_gf = 0.
    max_no_comp = 0

!   NOTE: all data created here refer to the original files (id,xlam,gf)
!   1st run do define no. of unpacked components and to calculate sum(gf)
loop1: Do i = 1, nftot

      If (i==1) Then
        dlp = 0.5D0*(lam_cmf(i+1)-lam_cmf(i))
        lmin = lammin
        lmax = lam_cmf(i) + dlp
!       look for first line in interval
        Do ic = 1, ntotnl3
          If (xlam(ic)>=lmin) Exit
        End Do

      Else If (i==nftot) Then
        dlm = dlp
        lmin = lmax
        lmax = lam_cmf(i)

      Else
        dlm = dlp
        lmin = lmax
        dlp = 0.5D0*(lam_cmf(i+1)-lam_cmf(i))
        lmax = lam_cmf(i) + dlp
      End If

100   Continue

!     might be programmed faster (checking only the upper limit), &
!     but note that certain intervals can be empty
      If (xlam(ic)>=lmin .And. xlam(ic)<lmax) Then
        in = id_rbb(ic)
        If (in/=0) Then
          met_rbb(in)%no_comp = met_rbb(in)%no_comp + 1
          max_no_comp = max(max_no_comp, met_rbb(in)%no_comp)
          met_rbb(in)%sum_gf = met_rbb(in)%sum_gf + gf(ic)
        End If
        ic = ic + 1
!       out of range
        If (xlam(ic)>lammax) Exit loop1
        Go To 100
      End If

    End Do loop1

    Write (999, *)
    Write (999, *) &
      ' max. number of unpacked components of NLTE bg-elements = ', &
      max_no_comp
    Print *
    Print *, ' max. number of unpacked components of NLTE bg-elements = ', &
      max_no_comp

!   test whether sum gf_i (unpacked) = gf (packed) ...
!   note1: this test can be wrong at the borders, since unpacked levels are
!   only considered between WAVBLUE and WAVRED. If a packed level
!   is close to the border, it might contain contributions from inside
!   and outside this interval, where the latter are not considered.
!   note2: With the present database, there are some cases where the
!   gf-values are considerably different, or, even worse, there is
!   no transition in the long list where there should be one according
!   to the short list (packed). This has been improved with the new
!   database (checked together with Tadziu), so that at some point
!   the new database needs to be implemented!

!   ... and define indices where info on unpacked transitions will be located

    met_rbb%quality = -1
    met_rbb%index_comp = -1

    iqual0 = 0
    iqual1 = 0
    iqual2 = 0
    tot_comp = 0

    irun = 1

    Do k = 1, natom
      If (jatom_full(k)==0) Cycle
      Do j = 1, jmax(k)
        n = indrec(k, j)
        jj = irbb(n) - 1

        Do ii = 1, ilin(n)
          ji = jj + ii
          If (met_rbb(ji)%wave>wavblue .And. met_rbb(ji)%wave<wavred) Then
            tot_comp = tot_comp + met_rbb(ji)%no_comp
!           this is the index where info on the unpacked transitions will be
!           located
            met_rbb(ji)%index_comp = irun
!           includes case with no_comp=0
            irun = irun + met_rbb(ji)%no_comp
            If (abs(1.-met_rbb(ji)%sum_gf/met_rbb(ji)%gf)>0.2) Then
!             print*,k,j,ji,met_rbb(ji)%wave,met_rbb(ji)%no_comp,met_rbb(ji)%sum_gf,met_rbb(ji)%gf
              If (met_rbb(ji)%sum_gf==0.) Then
                met_rbb(ji)%quality = 2
                iqual2 = iqual2 + 1
              Else
                met_rbb(ji)%quality = 1
                iqual1 = iqual1 + 1
              End If
            Else
              met_rbb(ji)%quality = 0
              iqual0 = iqual0 + 1
            End If
          End If
        End Do

      End Do
    End Do
    Write (999, *) &
      ' total number of depacked components of NLTE bg-elements = ', tot_comp
    Print *, ' total number of depacked components of NLTE bg-elements = ', &
      tot_comp

    Write (999, *)
    Write (999, *) ' no. of packed transitions of selected elements &
      &(within limits, but all ions)'
    Write (999, *) ' with quality = 0 (OK):', iqual0
    Write (999, *) ' with quality = 1 (large differences in Sum gf_i and gf):' &
      , iqual1
    Write (999, *) ' with quality = 2 (no transition in unpacked data):', &
      iqual2
    Write (999, *)
    Print *
    Print *, ' no. of packed transitions of selected elements (within &
      &limits, but all ions)'
    Print *, ' with quality = 0 (OK):', iqual0
    Print *, ' with quality = 1 (large differences in Sum gf_i and gf):', &
      iqual1
    Print *, ' with quality = 2 (no transition in unpacked data):', iqual2
    Print *

!   iqual=2 implies no_comp=0, sumgf=0.
!   iqual=1,0 implies no_comp ge 1, sumgf gt 0.
!   not considered: iqual=-1, no_comp=0, sumgf=0.

!   these are the files where info on unpacked transitions will be stored
!   unpacked_index_cmf: index of corresponding cmf_freq. (roughly xlam)
!   unpacked_gf: gf-value of individual component
    Allocate (unpacked_index_cmf(tot_comp), unpacked_gf(tot_comp), &
      unpacked_xlam(tot_comp))
    unpacked_index_cmf = -1

!   2nd run do define cmf-index and gf values for individual components
loop2: Do i = 1, nftot

      If (i==1) Then
        dlp = 0.5D0*(lam_cmf(i+1)-lam_cmf(i))
        lmin = lammin
        lmax = lam_cmf(i) + dlp
!       look for first line in interval
        Do ic = 1, ntotnl3
          If (xlam(ic)>=lmin) Exit
        End Do

      Else If (i==nftot) Then
        dlm = dlp
        lmin = lmax
        lmax = lam_cmf(i)

      Else
        dlm = dlp
        lmin = lmax
        dlp = 0.5D0*(lam_cmf(i+1)-lam_cmf(i))
        lmax = lam_cmf(i) + dlp
      End If

110   Continue

!     might be programmed faster (checking only the upper limit), &
!     but note that certain intervals can be empty
      If (xlam(ic)>=lmin .And. xlam(ic)<lmax) Then
        in = id_rbb(ic)

!       bug found by Jon July 2015
        If (in/=0) Then
          If (met_rbb(in)%wave>wavblue .And. met_rbb(in)%wave<wavred) Then
            no_comp = met_rbb(in)%no_comp
            If (no_comp==0 .And. met_rbb(in)%quality/=2) Then
              Write (999, *) ' stop: no_comp = 0 and quality ne 2'
              Stop ' no_comp = 0 and quality ne 2'
            End If
            index_comp = met_rbb(in)%index_comp
!           includes case that no_comp=0
            Do ii = index_comp, index_comp + no_comp - 1
              If (unpacked_index_cmf(ii)==-1) Go To 120
            End Do
            Print *, index_comp, no_comp, ii, met_rbb(in)%wave
            Write (999, *) &
              ' stop: empty index for unpacked component not found'
            Stop ' empty index for unpacked component not found'

120         unpacked_index_cmf(ii) = i ! index of cmf_frequency
            unpacked_gf(ii) = gf(ic) ! gf value of individual component
            unpacked_xlam(ii) = xlam(ic) ! wavelength of individual component
          End If
        End If
        ic = ic + 1
!       out of range
        If (xlam(ic)>lammax) Exit loop2
        Go To 110
      End If

    End Do loop2

!   consistency tests, can be commented out later
    Do k = 1, natom
      If (jatom_full(k)==0) Cycle
      Do j = 1, jmax(k)
        n = indrec(k, j)
        jj = irbb(n) - 1

        Do ii = 1, ilin(n)
          ji = jj + ii
          If (met_rbb(ji)%quality>=0 .And. (met_rbb( &
            ji)%wave<=wavblue .Or. met_rbb(ji)%wave>=wavred)) Then
            Write (999, *) ' stop: quality ge 0 outside interval'
            Stop ' quality ge 0 outside interval'
          End If

          If (met_rbb(ji)%wave>wavblue .And. met_rbb(ji)%wave<wavred) Then
            If (met_rbb(ji)%quality==-1) Then
              Write (999, *) ' stop: quality = -1 inside interval'
              Stop ' quality = -1 inside interval'
            End If

            no_comp = met_rbb(ji)%no_comp
            index_comp = met_rbb(ji)%index_comp
!           includes case that no_comp=0
            sum_gf = 0.
            Do in = index_comp, index_comp + no_comp - 1
              sum_gf = sum_gf + unpacked_gf(in)
            End Do
            If (sum_gf/=met_rbb(ji)%sum_gf) Then
              Print *, k, j, ji, met_rbb(ji)%wave, met_rbb(ji)%no_comp, &
                met_rbb(ji)%sum_gf, sum_gf
            End If
          End If
        End Do

      End Do
    End Do

!   JO shifted out of 'first' block, since lwion1 can change
!   nsize=size(met_rbb%quality)
!   allocate(tlu_met(nsize,lwion1-1),tul_met(nsize,lwion1-1))

!   now we have acquired all necessary info. The rest needs to be calculated
!   as long as calcion varies (since xlam1 varies)

    first = .False.
  End If
! --------------------------------------------------------------------

! calculate indices for depacked levels; this has do be done out of 'first'
! block, since depacked%index changes as long as line_list will be called

  depacked%index = 0 !              reset (data from subr. line_list no longer
! needed)
  cln = log(1.D0+vref/clight)
  lamax = lam_cmf(nftot)
  Do i = 1, lineno
    wave = depacked(i)%wave
    If (wave>=lammin .And. wave<=lamax) Then
      ii = log(depacked(i)%wave/lammin)/cln + 1
      dlp = wave - lam_cmf(ii)
      If (dlp<0.) Then
        Write (999, *) ' stop: dlp < 0'
        Stop 'dlp < 0'
      End If
      dlm = lam_cmf(ii+1) - wave
      If (dlm<0.) Then
        Write (999, *) ' stop: dlm < 0'
        Stop 'dlm < 0'
      End If
      If (dlp<dlm) Then
        depacked(i)%index = ii
      Else
        depacked(i)%index = ii + 1
      End If
!     print*,wave,ii,depacked(i)%index,lam_cmf(ii),lam_cmf(ii+1)
    End If
  End Do

! JO changed Nov. 2015
  nsize = size(met_rbb%quality)
! lwion1 can change
  If (.Not. allocated(tlu_met)) Then
    Allocate (tlu_met(nsize,lwion1-1), tul_met(nsize,lwion1-1))
  Else
    Deallocate (tlu_met, tul_met)
    Allocate (tlu_met(nsize,lwion1-1), tul_met(nsize,lwion1-1))
  End If

  If (nftot==0) Then
    Write (999, *) ' stop: nftot = 0 after first cmf_all iteration'
    Stop ' nftot = 0 after first cmf_all iteration'
  End If

  index_lam_cmf = 0

loop: Do i = 1, nftot

    If (i==1) Then
      dlp = 0.5D0*(lam_cmf(i+1)-lam_cmf(i))
      lmin = lammin
      lmax = lam_cmf(i) + dlp
!     look for first line in interval
      Do ic = 1, ntot
!       if(xlam1(ic).ge.lmin) exit
!       changed by JO March 2015: wavblue is the actual begin
        If (xlam1(ic)>=wavblue) Exit
      End Do
      icfirst = ic
      lastindex = icfirst - 1
!     icfirst is first valid index

    Else If (i==nftot) Then
      dlm = dlp
      lmin = lmax
      lmax = lam_cmf(i)

    Else
      dlm = dlp
      lmin = lmax
      dlp = 0.5D0*(lam_cmf(i+1)-lam_cmf(i))
      lmax = lam_cmf(i) + dlp
    End If

130 Continue

!   might be programmed faster (checking only the upper limit), &
!   but note that certain intervals can be empty
    If (xlam1(ic)>=lmin .And. xlam1(ic)<lmax) Then
      index_lam_cmf(i) = ic
      ic = ic + 1
!     out of range
!     if (xlam1(ic).gt.lammax) exit loop
!     changed by JO March 2015: wavred is the actual end
      If (xlam1(ic)>wavred) Exit loop
      Go To 130
!     no else required; in case there is no line. index_lam_cmf = 0
    End If

!   consistency check after sub-interval finished;
!   might be skipped when everything is working
    in = index_lam_cmf(i)
    If (in/=0) Then
      If (xlam1(in)>lmax .Or. xlam1(in)<lmin) Then
        Write (999, *) ' stop: error1 in fgrid_cmf'
        Stop ' error1 in fgrid_cmf'
      End If
      If (xlam1(lastindex+1)>lmax .Or. xlam1(lastindex+1)<lmin) Then
        Write (999, *) ' stop: error2 in fgrid_cmf'
        Stop ' error2 in fgrid_cmf'
      End If
      lastindex = in
    End If

  End Do loop

  iclast = ic - 1 !                 is the last valid index

! final consistency check
  If (iclast>ntot) Then
    Write (999, *) ' stop: ic-1 > ntot in fgrid_cmf'
    Stop ' ic-1 > ntot in fgrid_cmf'
  End If

! some statistics, might be skipped when everything is working
  ic0 = 0
  din = 0
  dinmax = 0
  lastindex = icfirst - 1

  Do i = 1, nftot

    in = index_lam_cmf(i)
    If (in==0) Then
      ic0 = ic0 + 1
    Else
      din1 = in - lastindex !       actually, in - (lastindex+1) + 1
      din = din + din1
      dinmax = max(din1, dinmax)
      lastindex = in
    End If
  End Do

  If (iclast-icfirst+1/=din) Then
    Write (999, *) ' stop: something wrong with index finder'
    Stop ' something wrong with index finder'
  End If

  Write (999, *) ' considered range for complete coupling: ', lammin, lammax
  Write (999, *) ' number of subintervals                = ', nftot
  Write (999, *) ' corresponding to no. of single lines  = ', nftot/21 + 1
  Write (999, *) ' total number of lines in subintervals = ', din, &
    ' (gf_min = ', gf_min, ')'
  Write (999, *) ' number of subintervals with no lines  = ', ic0
  Write (999, *) ' max. number of lines in subintervals  = ', dinmax
  Write (999, *)
  Print *, ' considered range for complete coupling: ', lammin, lammax
  Print *, ' number of subintervals                = ', nftot
  Print *, ' corresponding to no. of single lines  = ', nftot/21 + 1
  Print *, ' total number of lines in subintervals = ', din, ' (gf_min = ', &
    gf_min, ')'
  Print *, ' number of subintervals with no lines  = ', ic0
  Print *, ' max. number of lines in subintervals  = ', dinmax
  Print *
  Return

! this is an example to find individual lines
  lastindex = icfirst - 1
  Do i = 1, nftot
    in = index_lam_cmf(i)
    If (lam_cmf(i)>4633 .And. lam_cmf(i)<4645) Then
      Print *, lam_cmf(i), in
      If (in/=0) Then
        Do j = lastindex + 1, in
          Print *, j, xlam1(j), id1(j)
        End Do
      End If
    End If
    If (in/=0) lastindex = in
  End Do
  Stop
End Subroutine

!-----------------------------------------------------------------------

Subroutine profile_cmf(nd, temp)

! calculates Doppler profiles (used within subr. opacity_cmf)
! for all required elements as a a function of depth

! IMPORTANT NOTE
! so far, only Doppler profiles. At least for H/He, Stark broadening needs
! to be included to obtain realistic fluxes at Lyman series limit

! JO: presumably, pweightdop is not needed in final version
! test can be skipped after experience with different models has accumulated

  Use :: nlte_type
  Use :: nlte_dim, Only: id_ndept
  Use :: fund_const, Only: clight, wpi => sqpi
  Use :: fastwind_params, Only: names1, natom, xmaxdop
  Use :: princesa_var, Only: nat
  Use :: nlte_var, Only: vmax, vturb, vther
  Use :: nlte_app, Only: teff, indexel, vth
  Use :: cmf_all_var, Only: indexel_inv, nftot, lam_cmf, vturb_cmf, vref, &
    vth_cmf, elements_cmf, nf_cmf, profdop, profdop_h, profdop_he, pweightdop, &
    pweightdop_h, pweightdop_he, sum_pweightdop, sum_pweightdop_h, &
    sum_pweightdop_he, lamtest, weighttest, vdoptot

  Implicit None

  Integer (i4b), Parameter :: nd1 = id_ndept

! ..
! .. scalar arguments ..
  Integer (i4b), Intent (In) :: nd
! ..
! .. array arguments ..
  Real (dp), Dimension (nd1), Intent (In) :: temp

  Integer (i4b) :: k, kk, i, j, l, ixmax, ixmax_h, ixmax_he, maxixmax, kfirst, &
    klast

  Integer (i4b) :: nsafety = 2 !    add. freq. points for wider profile
! (interior)

  Real (dp) :: x, eps, err, vdop, xmax, vdopl, del, del1, xk, ws, dx, wnue1, &
    wnue2, sump, maxerr

  Logical :: first
  Logical, Parameter :: test = .True.

  Data first/.True./

  If (first) Then
!   index_array providing index of explict atom w.r.t. background elements
!   now already performed in subr. line_list
!   if (id_atoms .ne. nat) stop ' id_atoms ne nat (subr. profile_cmf)'
!   do k = 1,nat
!   do i = 1,natom
!   if(labat(k).eq.names1(i)) then
!   indexel_inv(k)=i
!   endif
!   enddo
!   test
!   if (indexel(indexel_inv(k)).ne.k) stop ' something wrong with
!   index_el(inv), H or He not explicit? (subr. profile_cmf)'
!   enddo

!   define vth_cmf
    vth_cmf = vth
    Do k = 1, nat
      i = indexel_inv(k)
      vth_cmf(i) = vther(k) !       overwrite explicit elements
    End Do

!   recalculate in case vturb ne vturb_cmf (usually, vturb is lower then)
    If (vturb/=vturb_cmf) Then
      Do k = 1, natom
        x = sqrt(vth_cmf(k)**2-vturb**2+vturb_cmf**2)
        vth_cmf(k) = x
      End Do
    End If

!   test to ensure that grid is equidistant (almost, to order vshift/c)
!   with respect to !x=(nu-nu0)/delta_nu_dop = c/vth*(1-lambda0/lambda), &
!   where vshift = n*vref
!   'exact' error (2nd order): -n*(n+1)/2*vref/c

!   test at red end, to also account for precision
    j = nftot - 11
    i = j + 10
    x = clight/vref*(1.D0-lam_cmf(j)/lam_cmf(i))
    eps = (x-dble(i-j))
    err = -(i-j)*(i-j+1)/2.*vref/clight
    If (abs(1.-eps/err)>0.01) Then
      Write (999, *) &
        ' stop: problems with equidistance of x (subr. profile_cmf)'
      Stop ' problems with equidistance of x (subr. profile_cmf)'
    End If

!   define maximum number of frequencies, and allocate profile arrays
!   (H and He separately, due to wider grid)


!   NOTE: we need to calculate only half of profile, due to symmetry!
    nf_cmf = elements_cmf

    ixmax_h = 0
    ixmax_he = 0
!   at first, profiles for bg elements
    Do k = 1, natom
!     only for bg-elements which are treated in the cmf
      If (elements_cmf(k)==0) Cycle
      vdop = vth_cmf(k)
      xmax = xmaxdop*vdop/vref
      ixmax = int(xmax) + 1 + nsafety
      If (indexel(k)>0) Then
        Write (999, *) &
          ' stop: something wrong with indexel-1 (subr. profile_cmf)'
        Stop ' something wrong with indexel-1 (subr. profile_cmf)'
      End If
      If (k==1) Then
        ixmax_h = ixmax
      Else If (k==2) Then
        ixmax_he = ixmax
      End If
!     for ALL bg elements (incl. H/He, if present)
      nf_cmf(k) = ixmax
    End Do

!   now, also for explicit elements
    Do k = 1, nat
      i = indexel_inv(k) !          i is the index w.r.t. bg elements
!     (ordered)
!     so far, elements_cmf contain only bg elements
      If (elements_cmf(i)/=0) Then
        Write (999, *) &
          ' stop: something wrong with elements_cmf (subr. profile_cmf)'
        Stop ' something wrong with elements_cmf (subr. profile_cmf)'
      End If
      vdop = vth_cmf(i)
      xmax = xmaxdop*vdop/vref
      ixmax = int(xmax) + 1 + nsafety
      If (indexel(i)<=0) Then
        Write (999, *) &
          ' stop: something wrong with indexel-2 (subr. profile_cmf)'
        Stop ' something wrong with indexel-2 (subr. profile_cmf)'
      End If
      If (i==1) Then
        ixmax_h = ixmax
      Else If (i==2) Then
        ixmax_he = ixmax
      End If
!     for ALL explicit elements (incl. H/He, if present)
      nf_cmf(i) = ixmax
    End Do

    If (ixmax_h/=0) Then
      Allocate (profdop_h(0:ixmax_h,nd))
      Allocate (pweightdop_h(0:ixmax_h,nd))
      Allocate (sum_pweightdop_h(-ixmax_h:ixmax_h))
!     for test
      Allocate (lamtest(0:ixmax_h), weighttest(0:ixmax_h))
      Print *, '(half) Doppler profile for H with ', ixmax_h + 1, &
        ' freq. points'
    End If
    If (ixmax_he/=0) Then
      Allocate (profdop_he(0:ixmax_he,nd))
      Allocate (pweightdop_he(0:ixmax_he,nd))
      Allocate (sum_pweightdop_he(-ixmax_he:ixmax_he))
      Print *, '(half) Doppler profile for He with ', ixmax_he + 1, &
        ' freq. points'
    End If

!   until here no change during iteration
    first = .False.

  Else !                            if not first
!   check elements 3 ... natom, since they may have changed

    nf_cmf(3:natom) = elements_cmf(3:natom)

    Do k = 3, natom
!     for all bg-elements (without H/He) which are treated in the cmf
      If (elements_cmf(k)==0) Cycle
      vdop = vth_cmf(k)
      xmax = xmaxdop*vdop/vref
      ixmax = int(xmax) + 1 + nsafety
      nf_cmf(k) = ixmax
    End Do
    Do k = 1, nat
      i = indexel_inv(k) !          i is the index w.r.t. bg elements
!     (ordered)
!     for all explicit elements (without H/He)
      If (i<=2) Cycle
      vdop = vth_cmf(i)
      xmax = xmaxdop*vdop/vref
      ixmax = int(xmax) + 1 + nsafety
      nf_cmf(i) = ixmax
    End Do

  End If

! for tests
! print*,'elements_cmf'
! print*,elements_cmf
! print*,'nf_cmf'
! print*,nf_cmf

! from here on, changes are possible (elements/temperature)
  maxixmax = 0
  kfirst = 0
  klast = 0

  Do k = 3, natom !                 now, check for ALL elements except for
!   H,He)
!   JO: here was a bug (cured Dec. 2016)
!   if(elements_cmf(k).eq.0) cycle
    If (nf_cmf(k)==0) Cycle
    If (kfirst==0) kfirst = k
    ixmax = nf_cmf(k)
    maxixmax = max0(maxixmax, ixmax)
    klast = k
  End Do
  If (.Not. allocated(profdop)) Then
    Allocate (profdop(0:maxixmax,nd,kfirst:klast))
    Allocate (pweightdop(0:maxixmax,nd,kfirst:klast))
    Allocate (sum_pweightdop(-maxixmax:maxixmax,kfirst:klast))
  Else
    Deallocate (profdop, pweightdop, sum_pweightdop)
    Allocate (profdop(0:maxixmax,nd,kfirst:klast))
    Allocate (pweightdop(0:maxixmax,nd,kfirst:klast))
    Allocate (sum_pweightdop(-maxixmax:maxixmax,kfirst:klast))
  End If

  Print *, '(half) Doppler profiles with max. ', maxixmax + 1, ' freq. points'
  Print *, '... for elements from ', names1(kfirst), ' to ', names1(klast)

! now everything is prepared, and the profiles and weights can be calculated
! as in subr. FGRID
elemloop: Do kk = 1, natom

    ixmax = nf_cmf(kk)
    If (ixmax==0) Cycle elemloop
    vdop = vth_cmf(kk)

    If (kk==1) Then

      Do l = 1, nd
        vdopl = (vdop**2-vturb_cmf**2)*temp(l)/teff
        vdopl = sqrt(vdopl+vturb_cmf**2)
        vdoptot(l, kk) = vdopl
        del = vmax/vdopl
        del1 = vref/vdopl !         x defined on vref-grid
        ws = .0D0
        Do k = 0, ixmax
!         x w.r.t. vref = k, tested above
          xk = k*del1
!         if(l.eq.43) print*,kk,k,xk  ! at l=43, T approx Teff
          pweightdop_h(k, l) = exp(-(xk**2))
          If (k/=0) ws = ws + pweightdop_h(k, l)
          profdop_h(k, l) = pweightdop_h(k, l)*del/wpi ! OPAL includes factor
!         lam*SR/vmax
        End Do

!       renormalization

        ws = 2.D0*ws + 1.D0 !       symmetry + zero point, exp(-0)
        Do k = 0, ixmax
          pweightdop_h(k, l) = pweightdop_h(k, l)/ws
        End Do
      End Do

    Else If (kk==2) Then

      Do l = 1, nd
        vdopl = (vdop**2-vturb_cmf**2)*temp(l)/teff
        vdopl = sqrt(vdopl+vturb_cmf**2)
        vdoptot(l, kk) = vdopl
        del = vmax/vdopl
        del1 = vref/vdopl !         x defined on vref-grid
        ws = .0D0
        Do k = 0, ixmax
          xk = k*del1
!         if(l.eq.43) print*,kk,k,xk
          pweightdop_he(k, l) = exp(-(xk**2))
          If (k/=0) ws = ws + pweightdop_he(k, l)
          profdop_he(k, l) = pweightdop_he(k, l)*del/wpi ! OPAL includes
!         factor lam*SR/vmax
        End Do

!       renormalization

        ws = 2.D0*ws + 1.D0 !       symmetry + zero point, exp(-0)
        Do k = 0, ixmax
          pweightdop_he(k, l) = pweightdop_he(k, l)/ws
        End Do
      End Do

    Else

      Do l = 1, nd
        vdopl = (vdop**2-vturb_cmf**2)*temp(l)/teff
        vdopl = sqrt(vdopl+vturb_cmf**2)
        vdoptot(l, kk) = vdopl
        del = vmax/vdopl
        del1 = vref/vdopl !         x defined on vref-grid
        ws = .0D0
        Do k = 0, ixmax
          xk = k*del1
!         if(l.eq.43) print*,kk,k,xk
          pweightdop(k, l, kk) = exp(-(xk**2))
          If (k/=0) ws = ws + pweightdop(k, l, kk)
          profdop(k, l, kk) = pweightdop(k, l, kk)*del/wpi ! OPAL includes
!         factor
!         lam*SR/vmax
        End Do

!       renormalization

        ws = 2.D0*ws + 1.D0 !       symmetry + zero point, exp(-0)
        Do k = 0, ixmax
          pweightdop(k, l, kk) = pweightdop(k, l, kk)/ws
        End Do
      End Do

    End If

  End Do elemloop

! calculation of sum(pweightdop) for H, He and metals at l=1
! (for optically thick outer boundary), allowing also for a test
  ws = 0.
  If (nf_cmf(1)/=0) Then
    Do k = -nf_cmf(1), nf_cmf(1)
      ws = ws + pweightdop_h(abs(k), 1)
      sum_pweightdop_h(k) = ws
    End Do
    If (abs(ws-1.D0)>1.D-10) Then
      Write (999, *) ' stop: error in pweightdop_H (subr. profile_cmf)'
      Stop ' error in pweightdop_H (subr. profile_cmf)'
    End If
  Else
    Write (999, *) ' WARNING: H absent, not test of pweightdop possible'
    Print *, ' WARNING: H absent, not test of pweightdop possible'
  End If

  ws = 0.
  If (nf_cmf(2)/=0) Then
    Do k = -nf_cmf(2), nf_cmf(2)
      ws = ws + pweightdop_he(abs(k), 1)
      sum_pweightdop_he(k) = ws
    End Do
    If (abs(ws-1.D0)>1.D-10) Then
      Write (999, *) ' stop: error in pweightdop_He (subr. profile_cmf)'
      Stop ' error in pweightdop_He (subr. profile_cmf)'
    End If
  Else
    Write (999, *) ' WARNING: He absent, not test of pweightdop possible'
    Print *, ' WARNING: He absent, not test of pweightdop possible'
  End If


  Do kk = kfirst, klast
    ws = 0.
    If (nf_cmf(kk)/=0) Then
      Do k = -nf_cmf(kk), nf_cmf(kk)
        ws = ws + pweightdop(abs(k), 1, kk)
        sum_pweightdop(k, kk) = ws
      End Do
      If (abs(ws-1.D0)>1.D-10) Then
        Write (999, *) ' stop: error in pweightdop (subr. profile_cmf)'
        Stop ' error in pweightdop (subr. profile_cmf)'
      End If
    End If

  End Do

! test whether profiles correctly normalized and sufficiently broad
  If (test) Then
!   set up test frequency grid, from proto-typical lam0 until
!   max. wavelength for H-profile (this is the broadest one)

    ixmax_h = nf_cmf(1)
    ixmax_he = nf_cmf(2)

    If (ixmax_h==0) Then
      Write (999, *) ' stop: ixmax_H = 0 (no hydrogen) in profile_cmf'
      Stop ' ixmax_H = 0 (no hydrogen) in profile_cmf'
    End If

    lamtest(0) = 1000.D8
    dx = 1.D0 + vref/clight
!   for weights, see subr. cmf_complete
    wnue1 = (1.D0-1.D0/dx)*clight*0.5
    wnue2 = (1.D0-1.D0/(dx**2))*clight*0.5

    Do k = 1, ixmax_h
      lamtest(k) = lamtest(0)*dx**k
    End Do

    Do k = 0, ixmax_h
      If (k==0) Then
        weighttest(k) = wnue1/lamtest(0)
      Else If (k==ixmax_h) Then
        weighttest(k) = wnue1/lamtest(ixmax_h-1)
      Else
        weighttest(k) = wnue2/lamtest(k-1)
      End If
    End Do

    dx = 2.*lamtest(0)/vmax

    Print *
testloop: Do kk = 1, natom
      ixmax = nf_cmf(kk)
      If (ixmax==0) Cycle testloop

      If (kk==1) Then
        If (ixmax/=ixmax_h) Then
          Write (999, *) ' stop: ixmax ne ixmax_H in profile_cmf'
          Stop ' ixmax ne ixmax_H in profile_cmf'
        End If
        maxerr = 0.
        Do l = 1, nd
          sump = 0.
          Do k = 0, ixmax_h
            sump = sump + profdop_h(k, l)*weighttest(k)
          End Do
          sump = sump*dx
          err = abs(1.-sump)
          maxerr = max(err, maxerr)
          If (err>0.01) Then
            Write (999, *) ' H', l, err
            Print *, ' H', l, err
            Write (999, *) ' stop: error in profile function (normalization)'
            Stop ' error in profile function (normalization)'
          End If
        End Do
        Print *, ' max. err in normalization of profile funct. (H)   = ', &
          maxerr
        maxerr = 0.

      Else If (kk==2) Then
        If (ixmax/=ixmax_he) Then
          Write (999, *) ' stop: ixmax ne ixmax_He in profile_cmf'
          Stop ' ixmax ne ixmax_He in profile_cmf'
        End If
        Do l = 1, nd
          sump = 0.
          Do k = 0, ixmax_he
            sump = sump + profdop_he(k, l)*weighttest(k)
          End Do
          sump = sump*dx
          err = abs(1.-sump)
          maxerr = max(err, maxerr)
          If (err>0.01) Then
            Write (999, *) ' He', l, err
            Write (999, *) ' stop: error in profile function (normalization)'
            Print *, ' He', l, err
            Stop ' error in profile function (normalization)'
          End If
        End Do
        Print *, ' max. err in normalization of profile funct. (He)  = ', &
          maxerr
        maxerr = 0.

      Else
        Do l = 1, nd
          sump = 0.
          Do k = 0, ixmax
            sump = sump + profdop(k, l, kk)*weighttest(k)
          End Do
          sump = sump*dx
          err = abs(1.-sump)
          maxerr = max(err, maxerr)
!         if(err.gt.0.001) then
!         Jo changed Nov. 2016
          If (err>0.01) Then
            Write (999, *) kk, l, err
            Write (999, *) ' stop: error in profile function (normalization)'
            Print *, kk, l, err
            Stop ' error in profile function (normalization)'
          End If
        End Do

      End If
    End Do testloop
    Print *, ' max. err in normalization of profile funct. (met) = ', maxerr

  End If !                          test

  Return
End Subroutine

!-----------------------------------------------------------------------

Subroutine opacity_cmf(nd, xne, te, clfac, r, velo, dvdr, test_overlap, ntest, &
  ipath)

! calculates summed line opacities and source-functions for all
! frequency grid points.
! Note that the transition frequencies are slightly shifted (see fgrid_cmf),
! to allow for a correct transport at peak opacity

! only line quantities, since cmf-transport requires separate continuum

! NOTE: I_core and f_edd determined from approx. treatment (CONT1), since
! average quantities required. Since f_edd for L < Lthin is approximated,
! this can lead to small disturbances at the end of convergence if
! LTHIN changes. Since lines affected, this can give rise to moderate
! changes in XJ_CMF (for a typical case, we found changes by 6%).
! To avoid this problem, we now calculate fedd always from XXK/XJ


! included:
! (i) correction of source-function for overlapping lines if up=0
! (ii) calculation of all depacked components according to LINES_xxx.dat
! (iiI) use Ng-extrapolated source-function if ng=.true.

  Use :: nlte_type
  Use :: nlte_dim, Only: id_ndept, id_llevs, id_nttrd
  Use :: nlte_opt, Only: opt_new_upper
  Use :: fund_const, Only: clight, hkl, hc2, c_tau, c_vanreg, c_forbid
  Use :: fastwind_params, Only: xm, minrho, gfcut ! see sumopal

  Use :: princesa_var, Only: gl, indat1, data
  Use :: nlte_var, Only: sr, vmax, fre, xj, xxk, lwion1, no_check, xne_ratio, &
    metals_converged, xj_smoothed => xj_cmf_coarse
  Use :: nlte_app, Only: indrec, occngold, occng, occnglte, jatom_full, &
    met_imin, met_imax, lwion, met_rbb, highest_level, ergion, metall
  Use :: cmf_all_var, Only: vref, lam_cmf, index_lam_cmf, &
    indfre_staggered_cmf, nftot, icfirst, id1, xlam1, gf1, indexel_inv, &
    nf_cmf, occexpl, ntot_expl_inside, profdop, profdop_h, profdop_he, &
    sum_pweightdop, sum_pweightdop_h, sum_pweightdop_he, sumopal_cmf, &
    sumetal_cmf, sumopal_outer, sline_outer, sumtaus_outer, indexcmf1, &
    indexcmf1_lam, opal_ntest, id1_overlap, no2, vdoptot, indxlamc_depacked, &
    xlamcmf_depacked
  Use :: ng_var, Only: ng, n_ng, trans_ng, sl_ng, index_ng, ng_met, n_ng_met, &
    trans_ng_met, sl_ng_met, index_ng_met

  Implicit None

  Integer (i4b), Parameter :: nrec = id_llevs
  Integer (i4b), Parameter :: nd1 = id_ndept

! real(dp), parameter :: c1=1.4388354967334d8,c2=3.972970127d8
! for Planck function, identical with function bnue
  Real (dp), Parameter :: c1 = hkl*1.D8, c2 = hc2*1.D24

  Integer (i4b), Intent (In) :: nd, ipath

  Real (dp), Dimension (nd1), Intent (In) :: xne, te, clfac, r, velo, dvdr

  Logical :: test_overlap

  Integer (i4b), Dimension (nd1) :: lll

  Logical, Dimension (nd1) :: check

  Real (dp), Dimension (nd1) :: opal, sline, occlow, occup, vr, tau0, xxj, &
    eddf, collfac, u0, vanreg, vanreg1, eps, eps1, opalexp, slineexp, &
    xxj_smoothed

  Real (dp), Parameter :: hc2_8 = hc2*1.D24 ! conversion to Angstrom

  Real (dp) :: srvmax, c_tau1, const, lam, lamgrid, test_dl, xlamgrid, opal1, &
    etal1, occlowi, occupi, taubi, taui, betacic, beta, bni, epsi, eps1i, &
    lam3, xxlam, precis, rel, lamexp, constexp, mirho, rho1, deltav, vdopmax, &
    sl, xlamc, fac, erg_fake, occupi_fake

  Integer (i4b) :: i, j, lastindex, in, ii, index, k, kk, ll, irest, ml, mu, &
    low, up, ixmax, irec, icount_expl, icount_bg, icount_bg_noup, inds_old, &
    inds, ik, ic, ntest, iiexp, indexexp, kkexp, kexp, iiexp_old, icount_ov, &
    ix, ifexp, i_ng, indi, up_fake

! print*,'ipath=',ipath

  If (.Not. allocated(sumopal_cmf)) Then
!   only once
    Allocate (sumopal_cmf(nd1,nftot), sumetal_cmf(nd1,nftot))
    Allocate (sumopal_outer(nftot), sline_outer(nftot))
    Allocate (sumtaus_outer(nftot))
  End If
! initialize
  If (ipath==1) sumopal_cmf = 0.
  sumetal_cmf = 0.
  sumopal_outer = 0.
  sline_outer = 0.
  sumtaus_outer = 0.

  If (ipath==1) id1_overlap = 0

  srvmax = sr/vmax
  c_tau1 = c_tau*srvmax
  test_dl = 0.5*vref/clight

  vr = velo/r
  tau0 = 1.D0/(dvdr/3.+2./3.*vr) !  at mue^2 = 1/3
  collfac = c_vanreg*xne/sqrt(te)

  check = .False.
! calculate check-variable (only for lwion (<lwion1)
  If (lwion1<lwion) Then
    Write (999, *) ' stop: something wrong with lwion1'
    Stop ' something wrong with lwion1'
  End If

  Do ll = lwion, nd
    rel = abs(1.D0-xne_ratio(ll))
    check(ll) = .Not. no_check .And. rel < 0.01
!   print*,'check',ll,check(ll),rel
  End Do

! read occupation numbers at all depth points
  Do ll = 1, nd
    Read (17, Rec=ll)(occexpl(i,ll), i=1, nrec)
    lll(ll) = ll !                  auxiliary vector
  End Do

  icount_expl = 0 !                 counter for lines from explicit elements
  icount_bg = 0 !                   counter for lines from background elements
  icount_bg_noup = 0 !              counter for lines from explicit elements
! without upper level
  icount_ov = 0 !                   counter for lines with up=0 and modified
! source-function

! loop over all cmf-frequencies, to calculate opal and etal
  lastindex = icfirst - 1
  inds_old = 0 !                    staggered

  ic = 0
  iiexp_old = 0

  Do i = 1, nftot
    in = index_lam_cmf(i)
    lamgrid = lam_cmf(i)
    inds = indfre_staggered_cmf(i) ! for calculation on staggered grid

!   lines shifted to lam to simplify profile calculations (uniform x)
!   and to obtain 'correct' line centers
    If (in/=0) Then
overlap: Do ii = lastindex + 1, in
        index = id1(ii)
        If (index>=1D8) Then
          If (ipath==2 .And. id1_overlap(ii)/=0) Then
            Write (999, *) ' stop: error(1) id1_overlap'
            Stop ' error(1) id1_overlap'
          End If
!         explicit element
          kk = index/1D8
          irest = index - kk*1D8
          ml = irest/10000
          mu = irest - ml*10000
          k = indexel_inv(kk)
!         in units of vref = vturb/3 (with vturb=max(vturb,vturbmin)
          ixmax = nf_cmf(k)
          lam = xlam1(ii)
          lam3 = lam**3
!         nevertheless, for consistency actual wavelengths used for
!         front-factors
!         if(xlam1(ii).gt.6000 .and. xlam1(ii).lt.7000)
!         print*,index,lam,k,ml,mu,ixmax
!         for tests
          If (abs(1.-lamgrid/lam)>test_dl) Then
            Write (999, *) ' stop: error1 in freq. shift (subr. opacity_cmf)'
            Stop ' error1 in freq. shift (subr. opacity_cmf)'
          End If
          const = c_tau1*lam*gf1(ii)
          occlow = occexpl(ml, :)/gl(ml)
          occup = occexpl(mu, :)/gl(mu)
          If (occlow(1)*occup(1)==0.) Then
!           this should never happen (at least for GLOBAL = .true.)
!           JO: check for GLOBAL = .false.
            Write (999, *) index, ' ', lam
            Print *, index, ' ', lam
            Write (999, *) &
              ' stop: occlow or occup eq 0 (expl. elements) in opacity_cmf'
            Stop ' occlow or occup eq 0 (expl. elements) in opacity_cmf'
          End If
!         pi e2/me c * 1.d-8 *SR/vmax * lam * gf * (nl/gl - nu/gu) /clfac
          opal = const*(occlow-occup)/clfac
          sline = hc2_8/lam3/(occlow/occup-1.D0)
          If (ng) Then
!           no check for almost_converged, since controlled in prep_ng
            Do i_ng = 1, n_ng
              If (ml==trans_ng(i_ng,1) .And. mu==trans_ng(i_ng,2)) Then
                sl = sl_ng(35, i_ng, 4)
                If (sl>0.D0) Then
!                 rescale to actual frequency (extrapolated source function
!                 calc. with xlamc)
                  indi = indat1(index_ng(i_ng))
                  xlamc = data(indi+1)
                  fac = xlamc**3/lam3
!                 print*,'sl_ng:',lam,xlamc,i_ng,sline(35),sl*fac
                  sline = sl_ng(:, i_ng, 4)*fac
                End If
                Exit
              End If
            End Do
          End If

          icount_expl = icount_expl + 1

!         save the freq. index for rateeq
          If (ic==0) Then
            Do ic = 1, id_nttrd
!             first line to be considered
              If (indexcmf1(indxlamc_depacked(ic))>=3) Exit
            End Do
          End If
          ik = indxlamc_depacked(ic)

!         JO Dec. 2016: account for Sobo-lines within range
          Do While (indexcmf1(ik)==1)
            ic = ic + 1
            ik = indxlamc_depacked(ic)
          End Do
          xxlam = xlamcmf_depacked(ik)

          If (abs(1.-lam/xxlam)<1.D-6) Then ! lam in single prec.
            If (indexcmf1(ik)<3) Then
              Write (999, *) ' stop: indexcmf1(ik) lt 3'
              Stop ' indexcmf1(ik) lt 3'
            End If

            indexcmf1_lam(ik) = i
!           print*,i,ic,ik,lam,xxlam,indexcmf1(ik)
!           print*
            ic = ic + 1
!           in case, use extrapolated source-function
          End If

        Else
!         background element
          k = index/1000000
          irest = index - k*1000000
          j = irest/100000
          irest = irest - j*100000
          low = irest/100
          up = irest - low*100
          ixmax = nf_cmf(k)
!         nevertheless, for consistency actual wavelengths used for
!         front-factors
          lam = xlam1(ii)
          lam3 = lam**3
!         if(xlam1(ii).gt.6000 .and. xlam1(ii).lt.7000)
!         print*,index,lam,k,j,low,up,ixmax
          If (abs(1.-lamgrid/lam)>test_dl) Then
            Write (999, *) ' stop: error2 in freq. shift (subr. opacity_cmf)'
            Stop ' error2 in freq. shift (subr. opacity_cmf)'
          End If
          irec = indrec(k, j)
          const = c_tau1*lam*gf1(ii)
          If (up/=0) Then
!           for tests: neglect opacity for corresponding OIII lines
!           if(k.eq.8.and.j.eq.3.and.low.eq.1) then
!           if(up.eq.5 .or. up.eq.6 .or. up.eq.8) const=c_tau1*lam*1.d-10
!           endif
            If (ipath==2 .And. id1_overlap(ii)/=0) Then
              Write (999, *) ' stop: error(2) id1_overlap'
              Stop ' error(2) id1_overlap'
            End If

!           lower and upper level  present

            If (jatom_full(k)==1) Then
!             selected elements, division lwion1
              Do ll = 1, lwion1 - 1
!               outside
                occlowi = occngold(irec, low, ll) ! occng = n/g
                If (j<met_imin(k,ll) .Or. j>met_imax(k,ll) .Or. occlowi==0.) &
                  Then
                  opal(ll) = 0.
                  sline(ll) = 0.
                Else
                  occupi = occngold(irec, up, ll)
                  opal(ll) = const*(occlowi-occupi)/clfac(ll)
                  sline(ll) = hc2_8/lam3/(occlowi/occupi-1.D0)
                End If
              End Do

              If (ng_met) Then
!               should never happen
                If (metals_converged) Then
                  Write (999, *) &
                    ' stop: metals_converged and ng_met in subr. opacity_cmf'
                  Stop ' metals_converged and ng_met in subr. opacity_cmf'
                End If
                Do i_ng = 1, n_ng_met
                  If (index==trans_ng_met(i_ng)) Then
                    sl = sl_ng_met(35, i_ng, 4)
                    If (sl>0.D0) Then
!                     rescale to actual frequency (extrapolated source
!                     function calc. with xlamc)
                      indi = index_ng_met(i_ng)
                      xlamc = met_rbb(indi)%wave
                      fac = xlamc**3/lam3
!                     print*,'sl_ng_met:',lam,xlamc,i_ng,sline(35),sl*fac
                      sline(1:lwion1-1) = sl_ng_met(:, i_ng, 4)*fac
                    End If
                    Exit
                  End If
                End Do
              End If

              Do ll = lwion1, nd
!               inside
                occlowi = occng(irec, low, ll) ! occng = n/g
                If (check(ll) .And. abs(1.-occlowi/occnglte(irec,low, &
                  ll))>0.05) Then
                  Write (999, *) ' lte-error1 ', k, j, irec, low, ll, occlowi, &
                    occnglte(irec, low, ll)
                  Write (999, *) &
                    ' stop: occng(selected, lower atmosphere) ne occnglte'
                  Print *, ' lte-error1 ', k, j, irec, low, ll, occlowi, &
                    occnglte(irec, low, ll)
                  Stop ' occng(selected, lower atmosphere) ne occnglte'
                End If
                occupi = occng(irec, up, ll)
                If (occlowi*occupi==0) Then
                  Write (999, *) index, ' ', lam
                  Write (999, *) ' stop: occlow or occup eq 0 (inner &
                    &selected elements) in opacity_cmf'
                  Print *, index, ' ', lam
                  Stop ' occlow or occup eq 0 (inner selected elements) &
                    &in opacity_cmf'
                End If
                opal(ll) = const*(occlowi-occupi)/clfac(ll)
                sline(ll) = c2/(exp(c1/lam/te(ll))-1.D0)/lam3
              End Do

            Else
!             other elements (approximate), division lwion
              Do ll = 1, lwion - 1
                occlowi = occng(irec, low, ll) ! occng = n/g
                occupi = occng(irec, up, ll)
                If (occlowi*occupi==0) Then ! change of ionization
                  Write (999, *) index, ' ', lam
                  Write (999, *) ' stop: occlow or occup eq 0 (outer &
                    &approx. elements) in opacity_cmf'
                  Print *, index, ' ', lam
                  Stop ' occlow or occup eq 0 (outer approx. elements) &
                    &in opacity_cmf'
                End If
                opal(ll) = const*(occlowi-occupi)/clfac(ll)
                sline(ll) = hc2_8/lam3/(occlowi/occupi-1.D0)
              End Do
              Do ll = lwion, nd
!               inside
                occlowi = occng(irec, low, ll) ! occng = n/g
                If (check(ll) .And. abs(1.-occlowi/occnglte(irec,low, &
                  ll))>0.05) Then
                  Write (999, *) ' lte-error2 ', k, j, irec, low, ll, occlowi, &
                    occnglte(irec, low, ll)
                  Write (999, *) &
                    ' stop: occng(approx., lower atmosphere) ne occnglte'
                  Print *, ' lte-error2 ', k, j, irec, low, ll, occlowi, &
                    occnglte(irec, low, ll)
                  Stop ' occng(approx., lower atmosphere) ne occnglte'
                End If
                occupi = occng(irec, up, ll)
                If (occlowi*occupi==0) Then
                  Write (999, *) index, ' ', lam
                  Write (999, *) ' stop: occlow or occup eq 0 (inner &
                    &approx.  elements) in opacity_cmf'
                  Print *, index, ' ', lam
                  Stop ' occlow or occup eq 0 (inner approx.  elements) &
                    &in opacity_cmf'
                End If
                opal(ll) = const*(occlowi-occupi)/clfac(ll)
                sline(ll) = c2/(exp(c1/lam/te(ll))-1.D0)/lam3
              End Do

            End If
            icount_bg = icount_bg + 1

!           only lower level  present

          Else

!           prepare two-level approach
!           sline approximated by two level atom; see subr. netmat_met
!           note: illuminating radiation field (for betacIc -> beta J)
!           from smooth (pseudo-)continuum
!           Nov. 2023: no longer from pseudo-continuum, but from smoothed
!           J_CMF
!           Feb. 2025: changed back to illumination from smooth pseudo cont.

            If (inds/=inds_old) Then
!             otherwise, already calculated; change approach when
!             interpolating in freq.,
!             and use indfre_cmf instead of indfre_staggered_cmf
              xxj = xj(:, inds)
!             JO Nov. 2023
              xxj_smoothed = xj_smoothed(:, inds)
!             JO: April 2018; now, the calculated Eddington factor is used
!             everywhere
!             (also outside lthin), to avoid oscillations)
!             lth=lthin(inds)
!             where(lll.lt.lth)
!             mustar=1.d0-(r(lth)/r)**2 !correction for lthin (checked)
!             mustar=sqrt(mustar)
!             eddf=(1.d0-mustar**3)/3.d0/(1.d0-mustar)! (checked)
!             elsewhere
              eddf = xxk(:, inds)/xxj
!             endwhere
              If (.Not. no_check .And. (eddf(nd)<0.31 .Or. eddf(nd)>0.35)) &
                Then
                Write (999, *) &
                  ' stop: inconsistent Edd. factor in opacity_cmf'
                Stop ' inconsistent Edd. factor in opacity_cmf'
              End If

!             OLD comment: eddf (approx.) checked, smooth transition between
!             outer and inner region,
!             goes to unity for large radii

!             prepare calculation of eps' = Cul/Aul (1-exp(-u0))
!             (see also subroutine opacitl and sumopal)
!             to save time, frequency at continuum grid
              xlamgrid = 1.D8/fre(inds)
!             JO (31. Oct. 2016): here was the bug
!             u0=1.4388d8/xlamgrid*te
              u0 = 1.4388D8/(xlamgrid*te)
              vanreg = collfac*xlamgrid**3*(1.D0-exp(-u0)) ! checked, OK
!             vanreg1 propto lam^2 needs to be finally divided by gf
              vanreg1 = vanreg*c_forbid/xlamgrid ! checked, OK,
              vanreg = vanreg/(1.+vanreg)
!             print*
!             print*,xlamgrid,vanreg(1),vanreg1(1)
!             print*,xlamgrid,vanreg(43),vanreg1(43)
!             print*,xlamgrid,vanreg(nd),vanreg1(nd)
              inds_old = inds
            End If

!           line specific quantities, at all depth-points

            If (gf1(ii)<=gfcut) Then
              eps = vanreg1/gf1(ii)
              eps = eps/(1.+eps)
            Else
              eps = vanreg
            End If

!           check for overlap with line from explict element,
!           in case modify source-function later on
            If (ipath==2) Then
              iiexp = id1_overlap(ii)
              If (iiexp/=0 .And. iiexp/=iiexp_old) Then
                indexexp = id1(iiexp)
                If (indexexp<1D8) Then
                  Write (999, *) &
                    ' stop: indexexp does not refer to explicit element'
                  Stop ' indexexp does not refer to explicit element'
                End If
                kkexp = indexexp/1D8
                kexp = indexel_inv(kkexp)
                If (kexp==1 .Or. kexp==2) Then
                  Write (999, *) ' stop: indexexp refers to H/He'
                  Stop ' indexexp refers to H/He'
                End If
                irest = indexexp - kkexp*1D8
                ml = irest/10000
                mu = irest - ml*10000
!               in units of vref = vturb/3 (with vturb=max(vturb,vturbmin)
                lamexp = xlam1(iiexp)
!               frequency index
                ifexp = nint(log(lamexp/lam_cmf(1))/log(1.D0+vref/clight)) + 1
                If (abs(1.-lam_cmf(ifexp)/lamexp)>test_dl) Then
                  Write (999, *) ' stop: error in ifexp (subr. opacity_cmf)'
                  Stop ' error in ifexp (subr. opacity_cmf)'
                End If
!               print*
!               print*,ml,mu,lam,lamexp
                constexp = c_tau1*lamexp*gf1(iiexp)
                occlow = occexpl(ml, :)/gl(ml)
                occup = occexpl(mu, :)/gl(mu)
!               pi e2/me c * 1.d-8 *SR/vmax * lam * gf * (nl/gl - nu/gu)
!               /clfac
                opalexp = constexp*(occlow-occup)/clfac
!               JO NOTE: thus far, we don't use extrapolated source-function
!               here,
!               for reasons of comp. time. Should not play any role,
!               since only no_up lines affected. Expected to become
!               consistent in course of iteration
                slineexp = hc2_8/lamexp**3/(occlow/occup-1.D0)
                deltav = abs(lam/lamexp-1.D0)*clight
!               print*,'new',index,indexexp,slineexp(1)
                icount_ov = icount_ov + 1
                iiexp_old = iiexp
              Else If (iiexp/=0 .And. iiexp==iiexp_old) Then
!               print*,'old',index,indexexp,slineexp(1)
                deltav = abs(lam/lamexp-1.D0)*clight
                icount_ov = icount_ov + 1
              End If
            End If !                ipath=2

            If (jatom_full(k)==1) Then
!             selected elements, division lwion1
ll_loop1:     Do ll = 1, lwion1 - 1
!               outside
                occlowi = occngold(irec, low, ll) ! occng = n/g
                If (j<met_imin(k,ll) .Or. j>met_imax(k,ll) .Or. occlowi==0.) &
                  Then
                  opal(ll) = 0.
                  sline(ll) = 0.
                Else
                  opal(ll) = const*occlowi/clfac(ll) ! includes gstat, ind.
!                 emission negl.
                  taubi = opal(ll)/(eddf(ll)*dvdr(ll)+(1.-eddf(ll))*vr(ll))
!                 at mue^2 = eddf
!                 pi e2/me c * 1.d-8 *SR/vmax * lam * gf * nl/gl / clfac /
!                 (eddf*dvdr + (1-eddf*v/r))

                  If (taubi<1.D-5) Then
!                   here, xxj is pseudo-continuum
                    betacic = xxj(ll)
!                   JO: Nov. 2023
!                   here (and in the following), xxj_smoothed is smoothed
!                   J_CMF
!                   JO: Feb. 2025 changed back, since otherwise too strong
!                   emission at 230 A for hot models
!                   betacic=xxj_smoothed(ll)
                  Else
                    betacic = (1.-exp(-taubi))/taubi*xxj(ll)
!                   JO: Nov. 2023/changed back
!                   betacic=(1.-exp(-taubi))/taubi*xxj_smoothed(ll)
                  End If

                  taui = opal(ll)*tau0(ll) ! at mue^2 = 1/3
!                 pi e2/me c * 1.d-8 *SR/vmax * lam * gf * nl/gl / clfac /
!                 (dvdr/3. + 2./3.*v/r)
                  If (taui<1.D-5) Then
                    beta = 1.D0
                  Else
                    beta = (1.-exp(-taui))/taui
                  End If

!                 now, we can calculate the two-level source-function. We
!                 inline Bnue, &
!                 and evaluate at actual frequency, to ensure correct
!                 thermalization

                  epsi = eps(ll)
                  eps1i = 1.D0 - epsi
                  bni = c2/(exp(c1/lam/te(ll))-1.D0)/lam3

                  sline(ll) = (eps1i*betacic+epsi*bni)/(epsi+beta*eps1i)

                  If (opt_new_upper) Then
!                   JO Nov. 2018
                    up_fake = highest_level(irec, low)
!                   if up_fake = 0, no info available, and we keep TLA
!                   approach.
!                   These lines are important, e.g., for hot stars around 1000
!                   A
!                   (e.g., from high lying Fe-levels)
                    If (up_fake/=0 .And. .Not. (low==1 .Or. metall(k,j, &
                      low)=='m')) Then
!                     if(up_fake.ne.0) then
!                     Otherwise, we assume that the upper level is in LTE
!                     w.r.t. up_fake
!                     In other words: the departures of up_fake and up are
!                     assumed to be identical
                      erg_fake = 1.D8/lam - (ergion(k,j,up_fake)-ergion(k,j, &
                        low))

                      If (erg_fake<0.) Then
!                       the upper level lies below up_fake
                        erg_fake = -erg_fake
!                       JO in case, allow for a somewhat larger range
                        If (erg_fake/ergion(k,j,up_fake)>0.1) Then
                          Write (999, *) k, j, low, lam
                          Write (999, *) erg_fake, ergion(k, j, up_fake), &
                            ergion(k, j, low)
                          Write (999, *) &
                            ' stop: upper level too much below up_fake'
                          Print *, k, j, low, lam
                          Print *, erg_fake, ergion(k, j, up_fake), &
                            ergion(k, j, low)
                          Stop ' upper level too much below up_fake'
                        End If
!                       fake upper level now more populated than up_fake
!                       (since below)
                        occupi_fake = occngold(irec, up_fake, ll)* &
                          exp(c1*erg_fake*1.D-8/te(ll))
                      Else
!                       typical situation: upper level less populated than
!                       up_fake
                        occupi_fake = occngold(irec, up_fake, ll)* &
                          exp(-c1*erg_fake*1.D-8/te(ll))
                      End If
                      opal(ll) = const*(occlowi-occupi_fake)/clfac(ll)
                      sline(ll) = hc2_8/lam3/(occlowi/occupi_fake-1.D0)
                    End If
                  End If
!                 end of opt_new_upper treatment


!                 check (similar to subr. overlap_noup) whether strong overlap
!                 with line from explicit element
                  If (ipath==2 .And. iiexp/=0) Then
                    vdopmax = xm*vdoptot(ll, kexp)
                    If (deltav>vdopmax) Cycle ll_loop1 ! outside
!                   depth-dependent range
!                   nothing to be done
                    ix = vdopmax/vref + 1 ! in units of vref/3
                    mirho = 1000.
                    Do ik = 0, ix
                      rho1 = opalexp(ll)*profdop(ik, ll, kexp)/ &
                        sumopal_cmf(ll, ifexp+ik)
!                     if(index.eq.26500400) then
!                     print*,ll,ifexp+ik,rho1,indexexp
!                     print*,opalexp(ll),profdop(ik,ll,kexp),sumopal_cmf(ll,ifexp+ik)
!                     endif
!                     significant or dominating inversion of bg opacities
                      If (rho1>1.1D0 .Or. rho1<0.D0) rho1 = 0.D0
                      mirho = min(mirho, rho1)
                    End Do
                    Do ik = 1, ix
                      rho1 = opalexp(ll)*profdop(ik, ll, kexp)/ &
                        sumopal_cmf(ll, ifexp-ik)
!                     if(index.eq.26500400) then
!                     print*,ll,ifexp-ik,rho1,indexexp
!                     print*,opalexp(ll),profdop(ik,ll,kexp),sumopal_cmf(ll,ifexp-ik)
!                     endif
!                     significant or dominating inversion of bg opacities
                      If (rho1>1.1D0 .Or. rho1<0.D0) rho1 = 0.D0
                      mirho = min(mirho, rho1)
                    End Do
!                   in case, modify source-function (approximating beta_tot by
!                   beta etc)
                    If (mirho>minrho) Then
!                     print*,ml,mu,ll,lamexp,mirho
                      sline(ll) = eps1i*((1.-beta)*slineexp(ll)+betacic) + &
                        epsi*bni
                    End If
                  End If !          ipath=2
                End If !            imin,imax
              End Do ll_loop1
              Do ll = lwion1, nd
!               inside
                occlowi = occng(irec, low, ll) ! occng = n/g
                If (check(ll) .And. abs(1.-occlowi/occnglte(irec,low, &
                  ll))>0.05) Then
                  Write (999, *) ' lte-error3 ', k, j, irec, low, ll, occlowi, &
                    occnglte(irec, low, ll)
                  Write (999, *) &
                    ' stop: occng(selected, lower atmosphere) ne occnglte'
                  Print *, ' lte-error3 ', k, j, irec, low, ll, occlowi, &
                    occnglte(irec, low, ll)
                  Stop ' occng(selected, lower atmosphere) ne occnglte'
                End If
                If (occlowi==0) Then
                  Write (999, *) index, ' ', lam
                  Write (999, *) ' stop: occlow eq 0 (inner selected &
                    &elements) in opacity_cmf'
                  Print *, index, ' ', lam
                  Stop ' occlow eq 0 (inner selected elements) in opacity_cmf'
                End If
                opal(ll) = const*occlowi/clfac(ll)
                sline(ll) = c2/(exp(c1/lam/te(ll))-1.D0)/lam3
              End Do

            Else
!             other elements (approximate), division lwion
ll_loop2:     Do ll = 1, lwion - 1
                occlowi = occng(irec, low, ll) ! occng = n/g
                If (occlowi==0) Then
                  Write (999, *) index, ' ', lam
                  Write (999, *) ' stop: occlow eq 0 (outer approx. &
                    &elements) in opacity_cmf'
                  Print *, index, ' ', lam
                  Stop ' occlow eq 0 (outer approx. elements) in opacity_cmf'
                End If
                opal(ll) = const*occlowi/clfac(ll) ! includes gstat, ind.
!               emission negl.
                taubi = opal(ll)/(eddf(ll)*dvdr(ll)+(1.-eddf(ll))*vr(ll)) ! at
!               mue^2
!               =
!               eddf
!               pi e2/me c * 1.d-8 *SR/vmax * lam * gf * nl/gl / clfac /
!               (eddf*dvdr + (1-eddf*v/r))

                If (taubi<1.D-5) Then
                  betacic = xxj(ll)
!                 JO: Nov. 2023/changed back
!                 betacic=xxj_smoothed(ll)
                Else
                  betacic = (1.-exp(-taubi))/taubi*xxj(ll)
!                 JO: Nov. 2023/changed back
!                 betacic=(1.-exp(-taubi))/taubi*xxj_smoothed(ll)
                End If

                taui = opal(ll)*tau0(ll) ! at mue^2 = 1/3
!               pi e2/me c * 1.d-8 *SR/vmax * lam * gf * nl/gl / clfac /
!               (dvdr/3. + 2./3.*v/r)
                If (taui<1.D-5) Then
                  beta = 1.D0
                Else
                  beta = (1.-exp(-taui))/taui
                End If

!               now, we can calculate the two-level source-function. We inline
!               Bnue, &
!               and evaluate at actual frequency, to ensure correct
!               thermalization

                epsi = eps(ll)
                eps1i = 1.D0 - epsi
                bni = c2/(exp(c1/lam/te(ll))-1.D0)/lam3
                sline(ll) = (eps1i*betacic+epsi*bni)/(epsi+beta*eps1i)

!               check (similar to subr. overlap_noup) whether strong overlap
!               with line from explicit element
                If (ipath==2 .And. iiexp/=0) Then
                  vdopmax = xm*vdoptot(ll, kexp)
                  If (deltav>vdopmax) Cycle ll_loop2 ! outside depth-dependent
!                 range
!                 nothing to be done
                  ix = vdopmax/vref + 1 ! in units of vref/3
                  mirho = 1000.
                  Do ik = 0, ix
                    rho1 = opalexp(ll)*profdop(ik, ll, kexp)/ &
                      sumopal_cmf(ll, ifexp+ik)
!                   significant or dominating inversion of bg opacities
                    If (rho1>1.1D0 .Or. rho1<0.D0) rho1 = 0.D0
                    mirho = min(mirho, rho1)
                  End Do
                  Do ik = 1, ix
                    rho1 = opalexp(ll)*profdop(ik, ll, kexp)/ &
                      sumopal_cmf(ll, ifexp-ik)
!                   significant or dominating inversion of bg opacities
                    If (rho1>1.1D0 .Or. rho1<0.D0) rho1 = 0.D0
                    mirho = min(mirho, rho1)
                  End Do
!                 in case, modify source-function (approximating beta_tot by
!                 beta etc)
                  If (mirho>minrho) Then
!                   print*,ml,mu,ll,lamexp,mirho
                    sline(ll) = eps1i*((1.-beta)*slineexp(ll)+betacic) + &
                      epsi*bni
                  End If
                End If !            ipath=2

              End Do ll_loop2
              Do ll = lwion, nd
!               inside
                occlowi = occng(irec, low, ll) ! occng = n/g
                If (check(ll) .And. abs(1.-occlowi/occnglte(irec,low, &
                  ll))>0.05) Then
                  Write (999, *) ' lte-error4 ', k, j, irec, low, ll, occlowi, &
                    occnglte(irec, low, ll)
                  Write (999, *) &
                    ' stop: occng(approx., lower atmosphere) ne occnglte'
                  Print *, ' lte-error4 ', k, j, irec, low, ll, occlowi, &
                    occnglte(irec, low, ll)
                  Stop ' occng(approx., lower atmosphere) ne occnglte'
                End If
                If (occlowi==0) Then
                  Write (999, *) index, ' ', lam
                  Write (999, *) ' stop: occlow eq 0 (inner approx. &
                    &elements) in opacity_cmf'
                  Print *, index, ' ', lam
                  Stop ' occlow eq 0 (inner approx. elements) in opacity_cmf'
                End If
                opal(ll) = const*occlowi/clfac(ll)
                sline(ll) = c2/(exp(c1/lam/te(ll))-1.D0)/lam3
              End Do
            End If

            icount_bg_noup = icount_bg_noup + 1
          End If !                  up eq or ne 0

        End If !                    explicit or bg element

        If (test_overlap) opal_ntest(ii) = opal(ntest)

!       if(ipath.eq.2 .and. xlam1(ii).ge.387.31 .and. xlam1(ii).le.387.40)
!       then! .and. gf1(ii).ge.0.05) then
!       print*,id1(ii),xlam1(ii)
!       do ll=1,nd
!       print*,ll,opal(ll),sline(ll)
!       enddo
!       endif


!       having calculated opacities and source functions, we can add them
!       (profile-weighted) to all affected frequencies +/- ixmax
!       in parallel, we calculate auxiliary variables for outer boundary,
!       namely
!       sum(opal), sum(sline*opal) and sum(opal*sum(pweight))
!       note that sum(pweight) := sum_pweight is defined from -ixmax...ixmax

!       for outer boundary
        opal1 = opal(1)
        etal1 = opal1*sline(1)

!       if(lam.gt.1545. .and. lam.lt.1555) then
!       print*,lam,k,j,opal1,sline(1)
!       endif
!       for tests
!       if (k.eq.12.and.j.eq.2.and.low.eq.3.and.up.eq.5) then
!       print*,lam
!       do ll=1,nd
!       print*,ll,opal(ll),sline(ll),met_imin(k,ll),met_imax(k,ll)
!       enddo
!       endif


!       sum_tausouter   is sum(opal * int(x...xmax) phi(x)dx) for iminus
!       sum_sline_outer is sum (opal*sline)/sum(opal) (without profile
!       function)
!       for iminus

!       i-ik<0 should not happen for all elements,
!       since 'safety region below wavblue and beyond wavred

!       calculate sumopal_cmf only in first run
        If (k==1) Then
          Do ik = 0, ixmax
            eps = opal*profdop_h(ik, :)
            eps1 = eps*sline
            If (ipath==1) sumopal_cmf(:, i+ik) = sumopal_cmf(:, i+ik) + eps
            sumetal_cmf(:, i+ik) = sumetal_cmf(:, i+ik) + eps1
!           for outer boundary
            sumopal_outer(i+ik) = sumopal_outer(i+ik) + opal1
            sline_outer(i+ik) = sline_outer(i+ik) + etal1
            sumtaus_outer(i+ik) = sumtaus_outer(i+ik) + &
              opal1*sum_pweightdop_h(ik)
            If (ik/=0) Then
              If (i-ik<=0) Then
                Write (999, *) ' stop: i-ik le 0 for H'
                Stop ' i-ik le 0 for H'
              End If
              If (ipath==1) sumopal_cmf(:, i-ik) = sumopal_cmf(:, i-ik) + eps
              sumetal_cmf(:, i-ik) = sumetal_cmf(:, i-ik) + eps1
!             for outer boundary
              sumopal_outer(i-ik) = sumopal_outer(i-ik) + opal1
              sline_outer(i-ik) = sline_outer(i-ik) + etal1
              sumtaus_outer(i-ik) = sumtaus_outer(i-ik) + &
                opal1*sum_pweightdop_h(-ik)
            End If
          End Do
!         print*,k,lam,ixmax

        Else If (k==2) Then
          Do ik = 0, ixmax
            eps = opal*profdop_he(ik, :)
            eps1 = eps*sline
            If (ipath==1) sumopal_cmf(:, i+ik) = sumopal_cmf(:, i+ik) + eps
            sumetal_cmf(:, i+ik) = sumetal_cmf(:, i+ik) + eps1
!           for outer boundary
            sumopal_outer(i+ik) = sumopal_outer(i+ik) + opal1
            sline_outer(i+ik) = sline_outer(i+ik) + etal1
            sumtaus_outer(i+ik) = sumtaus_outer(i+ik) + &
              opal1*sum_pweightdop_he(ik)
            If (ik/=0) Then
              If (i-ik<=0) Then
                Write (999, *) ' stop: i-ik le 0 for He'
                Stop ' i-ik le 0 for He'
              End If
              If (ipath==1) sumopal_cmf(:, i-ik) = sumopal_cmf(:, i-ik) + eps
              sumetal_cmf(:, i-ik) = sumetal_cmf(:, i-ik) + eps1
!             for outer boundary
              sumopal_outer(i-ik) = sumopal_outer(i-ik) + opal1
              sline_outer(i-ik) = sline_outer(i-ik) + etal1
              sumtaus_outer(i-ik) = sumtaus_outer(i-ik) + &
                opal1*sum_pweightdop_he(-ik)
            End If
          End Do
!         print*,k,lam,ixmax

        Else
          Do ik = 0, ixmax
            eps = opal*profdop(ik, :, k)
            eps1 = eps*sline
            If (ipath==1) sumopal_cmf(:, i+ik) = sumopal_cmf(:, i+ik) + eps
            sumetal_cmf(:, i+ik) = sumetal_cmf(:, i+ik) + eps1
!           for outer boundary
            sumopal_outer(i+ik) = sumopal_outer(i+ik) + opal1
            sline_outer(i+ik) = sline_outer(i+ik) + etal1
            sumtaus_outer(i+ik) = sumtaus_outer(i+ik) + &
              opal1*sum_pweightdop(ik, k)
            If (ik/=0) Then
              If (i-ik<=0) Then
                Write (999, *) ' stop: i-ik le 0 for metals'
                Stop ' i-ik le 0 for metals'
              End If
              If (ipath==1) sumopal_cmf(:, i-ik) = sumopal_cmf(:, i-ik) + eps
              sumetal_cmf(:, i-ik) = sumetal_cmf(:, i-ik) + eps1
!             for outer boundary
              sumopal_outer(i-ik) = sumopal_outer(i-ik) + opal1
              sline_outer(i-ik) = sline_outer(i-ik) + etal1
              sumtaus_outer(i-ik) = sumtaus_outer(i-ik) + &
                opal1*sum_pweightdop(-ik, k)
            End If
          End Do

        End If !                    k=1,2,else

      End Do overlap !              all overlapping lines at cmf freq. lamgrid
      lastindex = in
    End If !                        if lines present at lamgrid
  End Do !                          all cmf-frequencies


  If (icount_expl/=ntot_expl_inside) Then
    Write (999, *) icount_expl, ' ', ntot_expl_inside
    Write (999, *) &
      ' stop: not all lines from expl. elements found (subr. opacity_cmf)'
    Print *, icount_expl, ' ', ntot_expl_inside
    Stop ' not all lines from expl. elements found (subr. opacity_cmf)'
  End If

! finally, calculate line source function for outer boundary
  Where (sumopal_outer/=0.)
    sline_outer = sline_outer/sumopal_outer
  Elsewhere
    sline_outer = 0.
  End Where

! test whether no line(s) at first frequency (for blue wing boundary
! condition)
  Do ll = 1, nd
    If (sumopal_cmf(ll,1)/=0.) Then
      Write (999, *) ' stop: line at blue wing (subr. opacity_cmf)'
      Stop ' line at blue wing (subr. opacity_cmf)'
    End If
  End Do

  Write (999, *)
  Write (999, *) icount_expl + icount_bg + icount_bg_noup, &
    ' cmf line opacities/emissivites calculated'
  Write (999, *) icount_expl, ' lines from explicit elements'
  Write (999, *) icount_bg, ' lines from background elements with upper level'
  Write (999, *) icount_bg_noup, &
    ' lines from background elements with no upper level'
  Print *
  Print *, icount_expl + icount_bg + icount_bg_noup, &
    ' cmf line opacities/emissivites calculated'
  Print *, icount_expl, ' lines from explicit elements'
  Print *, icount_bg, ' lines from background elements with upper level'
  Print *, icount_bg_noup, &
    ' lines from background elements with no upper level'
  If (ipath==2) Then
    Write (999, *)
    Write (999, *) icount_ov, &
      ' lines with no upper level and modified source-function'
    Print *
    Print *, icount_ov, &
      ' lines with no upper level and modified source-function'
    If (icount_ov/=no2) Then
      Write (999, *) ' stop: icountov ne no2'
      Stop ' icountov ne no2'
    End If
  End If

! do i=1,nftot
! write(*,100) lam_cmf(i),sumetal_cmf(1,i)/sumopal_cmf(1,i), &
! &   sumetal_cmf(43,i)/sumopal_cmf(43,i),sumetal_cmf(nd,i)/sumopal_cmf(nd,i)
! enddo

! do i=1,nftot
! write(*,100) lam_cmf(i),sumopal_cmf(1,i),sumopal_cmf(43,i),sumopal_cmf(nd,i)
! enddo

! test whether all explicit lines have correct wavelengths
  precis = vref/clight
  Do j = 1, id_nttrd
    ii = indxlamc_depacked(j)
    If (indexcmf1(ii)>=3) Then
!     if(abs(1.-lam_cmf(indexcmf1_lam(ii))/xlamcmf_depacked(ii)).gt.precis) &
!     &     stop ' error in indexcmf1_lam'
      If (abs(1.-lam_cmf(indexcmf1_lam(ii))/xlamcmf_depacked(ii))>precis) &
        Print *, lam_cmf(indexcmf1_lam(ii)), xlamcmf_depacked(ii)
    End If
  End Do

  Return

100 Format (F14.4, 2X, 6(E12.6,2X))

End Subroutine

!----------------------------------------------------------------------------

Subroutine cmf_complete(nd, xne, temp, clfac, r, velo, dvdr, rho, taur, xnh, &
  pressure, ferr, dtfcorr, opt_ray)

! complete cmf transfer for all elements

! here, r, v, and dvdr are dimensionless (normalized),
! and the opacities/emissivities include line and continuum contributions.
! line contributions already profile-weighted (depth dependend) and summed

! both line and continuum quantities multiplied by Rstar
! line quantities weighted with Rstar/ delta_nue_dop_max = Rstar*lambda/vmax
! corrected in profile propto vmax/vdop(r) * exp(-xth(r)^2)/sqrt(pi)

! in these units, deltax = vref/vmax

! flux-like quantities (ah, ah_int, ...) from 1...nd+1 as usual
! will be interpolated to regular grid,
! with quantities xh, xh_int etc. (incl. factor r^2)

! radiative accelerations only on staggered grid, from 1...nd-1

! from Oct. 2025 on, new arrays add_corr (per frequency) and add_corr_ray
! that calculate and store corrections required for ray_complete and
! moments_cmf to avoid negative I_nu and J_nu

  Use :: nlte_type
  Use :: nlte_dim, Only: id_ndept, id_npoin
  Use :: fund_const, Only: amh, sigmae, clight, pi, sigsb, akb
  Use :: fastwind_params, Only: natom, optout_xj_xh, kmom_start

  Use :: run_once, Only: start_cmf_simple
  Use :: nlte_var, Only: sr, vmax, rtau23, corrfc, ifre, fre, wfre, &
    opac_nolines, etat_nolines, xj, xxh, modnam, opac_coarse => opac, &
    q1 => qq1, q2 => qq2, ndiv_calc, xj_smoothed => xj_cmf_coarse, xj_save, &
    kcmf_start, kcmf_end, restart_cmfall
  Use :: nlte_var, Only: z, p, tauz1, tau1, pp1, w0 => vp, w0last, w1 => vp1
  Use :: photstruc, Only: chibar_h
  Use :: nlte_app, Only: teff, abund, summas
  Use :: cmf_all_var, Only: vref, nftot, lam_cmf, indfre_cmf, sumopal_cmf, &
    sumetal_cmf, sumtaus_outer, xj_cmf, alo_cmf, u_outer, sline_outer, &
    xhinner_cmf, fedd, gedd, geddp, add_corr_arr, xh_cmf, hbound1, nbound1, &
    hboundnd, nboundnd

! JO June 2025
  Use :: cmf_multi, Only: iminus_cont, cmf_all_bluwi

  Implicit None
  Integer (i4b), Parameter :: nd1 = id_ndept, np = id_npoin

  Integer (i4b), Intent (In) :: nd

  Real (dp), Dimension (nd1), Intent (In) :: xne, temp, clfac, r, velo, dvdr, &
    rho, xnh, pressure

  Real (dp), Dimension (nd1), Intent (Inout) :: taur

  Real (dp), Dimension (nd1), Intent (Out) :: ferr, dtfcorr

  Real (dp), Allocatable, Dimension (:, :), Save :: xj_old

  Real (dp), Dimension (nd1) :: sigth, opac, etac, scont, x1o, x2o, expo, x1e, &
    x2e, expe, opalk, opakk, slinek, ak1, ucmf, pp, aj, ak, alo, efac, ydum, &
    uoldk, slinetotk, tauross, opaross, oparblue, oparred, pnew, gr, gt, &
    chibar_h_cmf

  Real (dp), Dimension (nd1-1) :: az1

  Real (dp), Dimension (nd1-1) :: vcmf, voldk, v1oldk

  Real (dp), Dimension (nd1-1) :: rad_force, cont_force, el_force, ydum1

  Real (dp), Dimension (nd1+1) :: ah, an, hold, ahweight, ah_int, r2, xh, &
    ah_ray, errh

  Real (dp), Dimension (nd1) :: xh1, xhlast, xh_int, maxchange_xj, aj_ray, &
    errj, jold, add_corr

  Real (dp), Dimension (nd1) :: xh_int_blue, xh_int_red

  Real (dp), Dimension (nd1, np-1) :: ubluwi, uold

  Real (dp), Dimension (nd1-1, np-1) :: vbluwi, vold, v1old

  Real (dp), Dimension (np-1) :: iminusold

  Integer (i4b), Dimension (np-1) :: ncount, nstep

  Logical, Dimension (np-1) :: in_line_old

  Real (dp) :: rl, xmue, xmue2, deltax, xlambda, aic, dbdr, e1, e2, de, dei, &
    aux, opa1, sl1, taus_outer, mubar, vc, e1obs, e2obs, xj1, xj2, summ, &
    summ1, xjapp, mean, dev, mu, sumabu, csound, devp, hbnd, nbnd

  Real (dp) :: y1, rdum, vdum, xfac, errmax, dtdr, nom_flux, fluxdif, &
    fldiff_blue, fldiff_red, planck, planck_blue, planck_red, tw, vc1, vcend, &
    wnue1, wnue2, weight, fre1, fre2, wlast_blue, wlast_red, const, &
    corrfc_cmf, xx1o, xx2o, xexpo, opacnd, opakknd, desh, xhmin

  Real (dp) :: bnue, dbdt, hplus, alo_old

  Integer (i4b) :: jp, lmax, l, lz, ind_old, ind, k, kk, k1, kblue, kred, &
    nsum, counter_u, counter_uneg, counter_v, counter_jneg, counter_alo_zero

  Logical :: inversion, cont_blue_thick, start, opt_ray, jneg

  Integer (i4b) :: icounter

  Character (3) :: counter

  Data icounter/0/

  Data start/.True./

! check that xj processed by obsfram treatment
  If (xj(1,ifre+1)/=1) Then
    Write (999, *) ' stop: XJ = XJ (CMF) IN CMF_COMPLETE'
    Stop ' XJ = XJ (CMF) IN CMF_COMPLETE'
  End If

! recalculate pressure

  sumabu = 0.
! note that summas/sumabu are here different from nlte/totout,
! due to different normalization
  Do k = 1, natom
    sumabu = sumabu + abund(k)/abund(1)
!   change of normalization;
!   abund' = abund/sum(abund) with abund N_k/N_H
!   abund  = abund' * sum(abund)
!   sum_abund = abund/abund' = (e.g.) abund(1)/abund'(1) = 1./abund'(1)
  End Do
  Print *
  Print *, ' sum(abundances, w.r.t. H) = ', sumabu
! recalculate pressure
  devp = 1.D100
  Do l = 1, nd
    mu = summas/(1.D0+xne(l)/xnh(l)/sumabu) ! note that xnh only approx.
    csound = akb/(mu*amh)
    pnew(l) = rho(l)*temp(l)*csound
    devp = min(devp, abs(1.-pnew(l)/pressure(l)))
  End Do

  Print *
  Print *, ' max dev. between p and pnew = ', devp
! now, pnew is the actual pressure
  If (devp>0.05) Then
    Write (999, *) ' stop: deviation in pressure too large'
    Stop ' deviation in pressure too large'
  End If

! restart_cmfall=.false.  !for tests with restart-input XJ(obs. frame)
! JO Nov. 2023: previous preparation/calculation of xj_smoothed shifted to
! subr. prep_smooth (called  before subr. opacity_cmf)
  If (start) Then
    If (restart_cmfall) Then
      start = .False.
    Else
      Continue
!     start set to .false. after end of freloop
    End If
  End If

! line-independent quantities
  If (start_cmf_simple) Then
!   might have already been calculated in cmf_simple (nlte_approx)

!   cmf_simple should have been called before
!   (if restart and metals_converged = T, this has not happened)
!   but: for very low z, no cmf-line might have been calculated.
!   THUS: we comment this stop statement.
!   if (optphotlines .and. .not. metals_converged) &
!   & stop ' optphotlines = T, metals_converged = F and start_cmf_simple = T
!   in cmf_complete'

jploop1: Do jp = 1, np - 1
      lmax = min0(np+1-jp, nd)
      lz = lmax - 1

      Do l = 1, lz
        tauz1(l, jp) = 1.D0/(z(l,jp)-z(l+1,jp))
      End Do

      Do l = 2, lz
        tau1(l, jp) = 2.D0/(z(l-1,jp)-z(l+1,jp))
      End Do

!     pp(l)=velocity gradient,projected on the present ray

      Do l = 1, lmax
        rl = r(l)
        xmue = z(l, jp)/rl
        xmue2 = xmue*xmue
        pp1(l, jp) = xmue2*dvdr(l) + (1.D0-xmue2)*velo(l)/rl
      End Do

    End Do jploop1

    start_cmf_simple = .False.
  End If

  sigth = xne*amh*sigmae/clfac !    corrected

  Do l = 1, nd
    rl = r(l)/rtau23
!   average mu to convert cmf to obs frequencies
!   assuming linear limb-darkening:
!   Int(1 ... mustar) (a+b mu) dmu =: a + b mubar => mubar = 0.5(1+mustar)
    If (rl>=1D0) Then
      mubar = 0.5D0*(1.D0+sqrt(1.D0-1.D0/rl**2))
!     in the outer wind, I_minus mostly zero
    Else
!     mubar=0.5d0
      mubar = 0.
!     in the innermost wind, I_plus and I_minus of same order
!     I_plus comes from bluer, and I_minus from redder frequencies
    End If
    vc = velo(l)*vmax/clight
    efac(l) = 1.D0/(1.D0-mubar*vc)
!   print*,l,mubar,efac(l)
  End Do

! r2 is r^2 on staggered grid
  r2(1) = r(1)*r(1)
  Do l = 1, nd - 1
    rdum = 0.5*(r(l)+r(l+1))
    r2(l+1) = rdum*rdum
  End Do
  r2(nd+1) = r(nd)*r(nd)

  deltax = vref/vmax

! for moments equation
  gr = dvdr/deltax
  gt = velo/r/deltax

  If (.Not. allocated(xj_cmf)) Then
!   only once
    Allocate (xj_cmf(nd1,nftot), alo_cmf(nd1,nftot), xhinner_cmf(nftot))
    Allocate (xj_old(nd1,nftot))
    Allocate (xh_cmf(nd1+1,nftot), fedd(nd1,nftot), gedd(nd1-1,nftot), &
      geddp(nd1-1,nftot))
    Allocate (add_corr_arr(nd1,nftot))
    Allocate (hbound1(nftot), nbound1(nftot), hboundnd(nftot), &
      nboundnd(nftot))
    Allocate (u_outer(np,nftot))
    xj_old = 0.
  End If

! ---------------------------------------------
! calculate CORRFC_CMF, by combining the ranges from the approx. treatment
! (lambda > wavblue and lambda < wavred), and the CMF range

! note: integration weights checked; perfect for blue and red range;
! accuracy for intermediate range O(10^-12)

! preparation of freq. integration weights
  vc1 = lam_cmf(2)/lam_cmf(1) - 1.D0
  vcend = lam_cmf(nftot)/lam_cmf(nftot-1) - 1.D0
  If (abs(vc1-vcend)>1.D-14) Then
    Write (999, *) ' stop: error in precision: lam_cmf'
    Stop ' error in precision: lam_cmf'
  End If

! compact notation
  vc1 = lam_cmf(1)/lam_cmf(2)
  wnue1 = (1.D0-vc1)*clight*0.5D8
  wnue2 = (1.D0-vc1*vc1)*clight*0.5D8

! blue range
  fre1 = 1.D8/lam_cmf(1)
  Do k = ifre, 1, -1
    If (fre(k)<fre1) Then
      kblue = k
      Exit
    End If
  End Do
  If (kblue+1>ifre) Then
    Write (999, *) ' stop: kblue+1 > ifre'
    Stop 'kblue+1 > ifre'
  End If

! renormalize weight for kblue+1
  tw = 0.D0
  Do k = ifre, kblue + 1, -1
    tw = tw + wfre(k)
  End Do
  weight = clight*(fre(ifre)-fre(kblue+1))
  wlast_blue = wfre(kblue+1) + weight - tw
! print*,wlast/wfre(kblue+1)
  If (wlast_blue<0.) Then
    Write (999, *) ' stop: wlast < 0 in blue range'
    Stop ' wlast < 0 in blue range'
  End If

! now integrate over blue range
  planck_blue = 0.D0
  fldiff_blue = 0.D0

  Do k = kblue + 1, ifre
    weight = wfre(k)
    If (k==kblue+1) weight = wlast_blue
    planck_blue = planck_blue + bnue(1.D8/fre(k), teff)*weight
    fldiff_blue = fldiff_blue + dbdt(1.D8/fre(k), temp(nd))/opac_coarse(nd, k) &
      *weight
  End Do

  weight = clight*(fre(kblue+1)-1.D8/lam_cmf(1))
  If (weight<=0.) Then
    Write (999, *) ' stop: weight < 0 (1)'
    Stop ' weight < 0 (1)'
  End If

  rdum = 0.5D0*(bnue(1.D8/fre(kblue+1),teff)+bnue(lam_cmf(1),teff))
  planck_blue = planck_blue + rdum*weight
  rdum = 0.5D0*(dbdt(1.D8/fre(kblue+1),temp(nd))+dbdt(lam_cmf(1),temp(nd)))/ &
    opac_coarse(nd, kblue+1) !      approximation (opakk approx opac_coarse)
  fldiff_blue = fldiff_blue + rdum*weight

! red range
  fre2 = 1.D8/lam_cmf(nftot)
  Do k = 1, ifre
    If (fre(k)>fre2) Then
      kred = k
      Exit
    End If
  End Do
  If (kred-1<1) Then
    Write (999, *) ' stop: kred-1 < 1'
    Stop 'kred-1 < 1'
  End If

! renormalize weight for kred-1
  tw = 0.D0
  Do k = 1, kred - 1
    tw = tw + wfre(k)
  End Do
  weight = clight*(fre(kred-1)-fre(1))
  wlast_red = wfre(kred-1) + weight - tw
! print*,wlast/wfre(kred-1)
  If (wlast_red<0.) Then
    Write (999, *) ' stop: wlast < 0 in blue range'
    Stop ' wlast < 0 in blue range'
  End If

! now integrate over red range
  planck_red = 0.D0
  fldiff_red = 0.D0

  Do k = 1, kred - 1
    weight = wfre(k)
    If (k==kred-1) weight = wlast_red
    planck_red = planck_red + bnue(1.D8/fre(k), teff)*weight
    fldiff_red = fldiff_red + dbdt(1.D8/fre(k), temp(nd))/opac_coarse(nd, k)* &
      weight
  End Do

  weight = clight*(1.D8/lam_cmf(nftot)-fre(kred-1))
  If (weight<=0.) Then
    Write (999, *) ' stop: weight < 0 (2)'
    Stop ' weight < 0 (2)'
  End If

  rdum = 0.5D0*(bnue(1.D8/fre(kred-1),teff)+bnue(lam_cmf(nftot),teff))
  planck_red = planck_red + rdum*weight
  rdum = 0.5D0*(dbdt(1.D8/fre(kred-1),temp(nd))+dbdt(lam_cmf(nftot),temp(nd))) &
    /opac_coarse(nd, kred-1) !      approximation (opakk approx opac_coarse)
  fldiff_red = fldiff_red + rdum*weight

  fluxdif = 0.D0
  planck = 0.D0
  errj = 0.
  errh = 0.

  ind_old = 0 !                     for interpolation

  If (opt_ray) Then
    Write (999, *)
    Write (999, *) ' CMF-transport: RAY-BY-RAY SOLUTION + MOMENTS EQS.'
    Write (999, *)
    Print *
    Print *, ' CMF-transport: RAY-BY-RAY SOLUTION + MOMENTS EQS.'
    Print *
  Else
    Write (999, *)
    Write (999, *) ' CMF-transport: MOMENTS EQS. ONLY'
    Write (999, *)
    Print *
    Print *, ' CMF-transport: MOMENTS EQS. ONLY'
    Print *
  End If

! cmf_range
  Do k = 1, nftot
!   calculate opacity at last depth point, in the same way as for
!   *all* depth points below

    xlambda = lam_cmf(k)
!   calculate continuum opacities in cmf
!   interpolate, to avoid steps in min(opa)
    ind = indfre_cmf(k)
    If (ind/=ind_old) Then
!     prepare interpolation weights (lamgrid always in between ind and ind-1)
      e1 = fre(ind)
      e2 = fre(ind-1)
      de = 1.D0/log10(e1/e2)
      xx1o = opac_nolines(nd, ind)
      xx2o = opac_nolines(nd, ind-1)
      xexpo = log10(xx1o/xx2o)*de
      xx1o = xx1o*sr !              before, sr cancels out
      ind_old = ind
    End If
!   actual interpolation (log-log):
!   x = x1*(nu/nu1)^[(log10(x1/x2)/log10(nu1/nu2)]
    dei = 1.D8/(xlambda*e1) !       here we use the shifted frequencies
    opacnd = xx1o*dei**xexpo !      clumping corrected, includes opath
    opakknd = sumopal_cmf(nd, k) + opacnd

!   no check for inversion

    If (k==1) Then
      weight = wnue1/lam_cmf(1)
    Else If (k==nftot) Then
      weight = wnue1/lam_cmf(nftot-1)
    Else
      weight = wnue2/lam_cmf(k-1)
    End If

!   lower flux in diffusion approx.
    y1 = dbdt(lam_cmf(k), temp(nd))/opakknd
    fluxdif = fluxdif + y1*weight

!   Planck function at Teff
!   y1  = bnue(lam_cmf(k),temp(nd))
    y1 = bnue(lam_cmf(k), teff)
    planck = planck + y1*weight

  End Do

! Perform some checks and calculate corrfc_cmf
  y1 = sigsb/pi*teff**4
  nom_flux = sigsb/(4.*pi)*teff**4*rtau23**2
! corresponds to BNUECON/4 in nlte.f90; correction, because r^2*flux = const
! and we compare to the nominal flux (also w.r.t. errors) at the lower
! boundary
! which is larger by a factor (rstar/sr)^2=rtau23^2
  dtdr = (temp(nd)-temp(nd-1))/(r(nd-1)-r(nd))

! check integration of Planck-function (i.e., freq. weights)
  planck = planck + planck_blue + planck_red
  errmax = 1.D0 - planck/y1
  Print *
  Print *, ' contribution to Planck-function from blue range = ', &
    planck_blue/planck
  Print *, ' contribution to Planck-function from  red range = ', &
    planck_red/planck
  Print *
  Print *, ' Error in integrated Planck-function (total range) = ', errmax
  Print *

  If (abs(errmax)>1.D-4) Then
    Write (999, *) &
      ' stop: something wrong with integration weights (cmf_complete)'
    Stop ' something wrong with integration weights (cmf_complete)'
  End If

! calculate corrfc_cmf

! Note: opac_coarse without SR
  const = corrfc*dtdr/(3.D0*sr)
  fldiff_blue = fldiff_blue*const
  fldiff_red = fldiff_red*const

  const = dtdr/(3.D0)
  fluxdif = fluxdif*const

  corrfc_cmf = (nom_flux-fldiff_blue-fldiff_red)/fluxdif

  Print *, ' CORRFC_CMF = ', corrfc_cmf
  Print *
  fluxdif = fluxdif*corrfc_cmf

  fluxdif = fluxdif + fldiff_blue + fldiff_red
  Print *, ' contribution to fluxdiff from blue range = ', fldiff_blue/fluxdif
  Print *, ' contribution to fluxdiff from  red range = ', fldiff_red/fluxdif
  Print *

! for tests (reformulate, since xxh now on grid)
! do k=1,ifre
! xh_int=xh_int+xxh(:,k)*wfre(k)
! const=(r(nd-1)-r(nd))*0.5d0*(opac_coarse(nd,k)+opac_coarse(nd-1,k))*sr
! const=(r(nd-1)-r(nd))*opac_coarse(nd,k)*sr
! rdum=(xj(nd,k)-xj(nd-1,k))/const
! vdum=(bnue(1.d8/fre(k),temp(nd))-bnue(1.d8/fre(k),temp(nd-1)))/const
! errmax=dbdt(1.d8/fre(k),temp(nd))*dtdr/(opac_coarse(nd,k)*sr)
! fluxdif=errmax/3.d0
! write(*,fmt='(i4,2x,f12.4,2x,3(e10.4,2x),f6.3,2x,e10.4,2x,f6.3)') &
! &  k,1.d8/fre(k),const,rdum,vdum,rdum/vdum,errmax,errmax/vdum
! write(*,fmt='(i4,2x,f12.4,2x,3(e10.4,2x),f6.3,2x,f6.3,2x,f6.3)') &
! &
! k,1.d8/fre(k),const,vdum,errmax,errmax/rdum,errmax/vdum,rdum/(3.*xxh(nd,k))
! enddo

! ---------------------------------------------

! initial values for approximation of Iminus
  iminusold = 0. !                  will be overwritten from values calculated
! in setup2/cont2
  in_line_old = .False. !           correct for pure cont
  ncount = 0 !                      correct for pure cont
  cont_blue_thick = .False.

  Do jp = 1, np - 1
    nstep(jp) = (1.-z(1,jp)/r(1))/deltax
  End Do

  ind_old = 0 !                     for interpolation

  rad_force = 0.D0
  cont_force = 0.D0
  el_force = 0.D0
  ah_int = 0.D0
  xh_int = 0.D0
  chibar_h = 0.D0
  fluxdif = 0.D0
  opaross = 0.D0

  counter_u = 0
  counter_v = 0
  counter_uneg = 0
  counter_jneg = 0
  counter_alo_zero = 0

freloop: Do k = 1, nftot

    xlambda = lam_cmf(k)

    If (mod(k,5000)==0) Print *, k, ' cmf_transport at ', xlambda

!   calculate continuum opacities/emissivities in cmf
!   interpolate, to avoid steps in min(opa,eta)
    ind = indfre_cmf(k)
    If (ind/=ind_old) Then
!     prepare interpolation weights (lamgrid always in between ind and ind-1)
      e1 = fre(ind)
      e2 = fre(ind-1)
      de = 1.D0/log10(e1/e2)
      x1o = opac_nolines(:, ind)
      x2o = opac_nolines(:, ind-1)
      expo = log10(x1o/x2o)*de
      x1o = x1o*sr !                before, sr cancels out

!     if START = .TRUE., i.e., 1st run THEN
!     interpolate xj onto corresponding obs. frame frequency
!     lam_cmf 'sees' radiation from bluer frequencies, lam_obs =
!     lam_cmf(1-mu*v/c)
!     e.g., at the CMF Lyman edge in the outer wind, the radiation field
!     orginates
!     from obs. frame frequencies around 905 A, and has a lower intensity than
!     using xj(lamobs=lamcmf)

!     instead of interpolating directly w.r.t xlambda, we
!     interpolate only the xj grid values w.r.t. the corresponding obs.
!     frequencies
!     (saves a lot of time)
!     could be certainly programmed somewhat more efficient, but programming
!     effort
!     not worth doing so

!     algorithm checked, does not make a lot of difference compared to the
!     approximate approach of using xj(:,ind) and xj(:,ind-1).
!     but one never knows ...

!     instead of xj we work here with xj_smoothed, to approx. account
!     for broadening by Thomson scattering

!     if START = .FALSE., no interpolation of XJ_SAVE required,
!     when line-frequency between KCMF_START and KCMF_END,
!     since the corresponding XJ values are already CMF values.
!     Outside this range, we interpolate w.r.t to the same freq. shift as
!     above,
!     since here the XJ values are obs. frame values.
      If (.Not. start) Then
        If (e2>=fre(kcmf_start) .And. e1<=fre(kcmf_end)) Then
!         .not. start and inside range: take direct values
          Do l = 1, nd
            x1e(l) = etat_nolines(l, ind) + xj_smoothed(l, ind)*sigth(l)
            x2e(l) = etat_nolines(l, ind-1) + xj_smoothed(l, ind-1)*sigth(l)
          End Do
          Go To 120
        End If
      End If
!     start or outside range: interpolate xj onto cmf-frequencies

      Do l = 1, nd
        e2obs = e2*efac(l)
        Do kk = ind - 1, ifre
          If (e2obs>=fre(kk) .And. e2obs<fre(kk+1)) Then
!           here we interpolate only linearly, to save time
            xj2 = xj_smoothed(l, kk) + (xj_smoothed(l,kk+1)-xj_smoothed(l,kk)) &
              /(fre(kk+1)-fre(kk))*(e2obs-fre(kk))
            Go To 100
          End If
        End Do
        Write (999, *) ' stop: not found (interpol. of xj2)'
        Stop ' not found (interpol. of xj2)'

100     e1obs = e1*efac(l)
        Do kk = ind, ifre
          If (e1obs>=fre(kk) .And. e1obs<fre(kk+1)) Then
!           here we interpolate only linearly, to save time
            xj1 = xj_smoothed(l, kk) + (xj_smoothed(l,kk+1)-xj_smoothed(l,kk)) &
              /(fre(kk+1)-fre(kk))*(e1obs-fre(kk))
            Go To 110
          End If
        End Do
        Write (999, *) ' stop: not found (interpol. of xj1)'
        Stop ' not found (interpol. of xj1)'

110     x1e(l) = etat_nolines(l, ind) + xj1*sigth(l)
        x2e(l) = etat_nolines(l, ind-1) + xj2*sigth(l)
      End Do

120   Continue
      expe = log10(x1e/x2e)*de
      x1e = x1e*sr !                before, sr cancels out

      ind_old = ind
    End If
!   actual interpolation (log-log):
!   x = x1*(nu/nu1)^[(log10(x1/x2)/log10(nu1/nu2)]
    dei = 1.D8/(xlambda*e1) !       here we use the shifted frequencies
    opac = x1o*dei**expo !          clumping corrected, includes opath
    etac = x1e*dei**expe !          clumping corrected, includes approx etath

!   for lower boundary
    Call diffus(xlambda, temp, r, nd, aic, dbdr)
!   at very first frequency, there is no line (checked in opacity_cmf)
!   then, continuum transfer can be performed in obs-frame
    scont = etac/opac

!   calculation of continuum radiation field at the bluest frequency

    If (k==1) Then
      Print *
      Print *, ' CORRECTION FACTOR FOR CONTINUUM AT WAVBLUE = ', corrfc
      Print *
      cmf_all_bluwi = .True.
      Call cont2(nd, np, r, p, z, scont, opac, aic, dbdr, corrfc, ubluwi)
      cmf_all_bluwi = .False.
    End If

    opalk = sumopal_cmf(:, k)
    opakk = opalk + opac
    slinek = (sumetal_cmf(:,k)+etac)/opakk

!   if(xlambda.ge.370. .and. xlambda.le.380) then
!   write(*,987) k,xlambda,opalk(43),opakk(43),slinek(43),etac(43)
!   endif
!   987 format('output',i5,f10.5,4(2x,e10.5))

!   total line source function, needed for alo in ray_complete
    Do l = 1, nd
      If (opalk(l)==0.) Then
        If (sumetal_cmf(l,k)/=0) Then
          Write (999, *) &
            ' stop: sumopal_cmf = 0 and sumetal_cmf ne 0 in ray_complete'
          Stop ' sumopal_cmf = 0 and sumetal_cmf ne 0 in ray_complete'
        End If
        slinetotk(l) = 0.
      Else
        slinetotk(l) = sumetal_cmf(l, k)/opalk(l)
      End If
    End Do

!   check for inversion
    inversion = .False.
    Do l = 1, nd
      If (opakk(l)<=0.) Then
        If (slinek(l)>0.) Then
          Write (999, *) &
            ' stop: opakk < 0 and slinek > 0 (subr. cmf_complete)'
          Stop ' opakk < 0 and slinek > 0 (subr. cmf_complete)'
        End If
        inversion = .True.
      End If
      If (slinek(l)<=0.) Then
        If (opakk(l)>0.) Then
          Write (999, *) &
            ' stop: slinek < 0 and opakk > 0 (subr. cmf_complete)'
          Stop ' slinek < 0 and opakk > 0 (subr. cmf_complete)'
        End If
        inversion = .True.
      End If
    End Do

    If (inversion) Then
      Print *, ' Inversion in cmf_complete at ', xlambda
!     stop ' Inversion in cmf_complete'
!     JO  this needs to be treated when occuring
!     hack, needs to be improved
      Where (opakk<0.)
        opalk = 0.
        opakk = opac
        slinek = etac/opac
        slinetotk = 0.
      End Where
    End If

    ak1 = 1./opakk

    Do l = 2, nd
      az1(l-1) = 2.D0/(opakk(l)+opakk(l-1))
    End Do


optray: If (opt_ray .Or. k==1) Then

      aj = 0.D0
      add_corr = 0.D0
      ah = 0.D0
      ak = 0.D0
      an = 0.D0
      alo = 0.D0
      hbnd = 0.D0
      nbnd = 0.D0

      opa1 = opac(1)
      sl1 = sline_outer(k)
      If (sl1<0.) sl1 = 0.
      taus_outer = sumtaus_outer(k)

      If (k==2 .And. opa1*r(1)/3.D0>1.) cont_blue_thick = .True.
!     optically thick continuum only for rho^2 processes

!     if(k.ge.98509 .and. k.le.98511) then
!     do l=1,nd
!     write(*,fmt='(i3,3(2x,e10.4))') l,opalk(l),slinetotk(l),slinek(l)
!     enddo
!     endif


jploop2: Do jp = 1, np - 1
        lmax = min0(np+1-jp, nd)
        lz = lmax - 1

        If (k==1) Then
!         only for bluewing boundary condition
          ucmf(1) = ubluwi(1, jp)
          Do l = 1, lz
            ucmf(l+1) = ubluwi(l+1, jp)
            aux = .5D0*(opac(l)+opac(l+1))
            vcmf(l) = (ucmf(l+1)-ucmf(l))/aux/(z(l,jp)-z(l+1,jp))
          End Do
!         for ALO
          uold(:, jp) = 0.D0
          vold(:, jp) = 0.D0
          v1old(:, jp) = 0.D0
          iminusold(jp) = iminus_cont(jp)
        Else
!         restore u and v
          ucmf(1:lmax) = ubluwi(1:lmax, jp)
          vcmf(1:lz) = vbluwi(1:lz, jp)
!         for ALO
          uoldk(1:lmax) = uold(1:lmax, jp)
          voldk(1:lz) = vold(1:lz, jp)
          v1oldk(1:lz) = v1old(1:lz, jp)
        End If

!       since deltax does not change, this can be calculated in advance
!       (later)
        pp(1:lmax) = pp1(1:lmax, jp)/deltax


!       for tests, to recover static solution
!       pp=1.d-50
!       sl1=0.
!       slinek=0.
!       slinetotk=0.

        Call ray_complete(k, jp, z(:,jp), r, nd, np, ucmf, vcmf, opa1, sl1, &
          taus_outer, iminusold(jp), in_line_old(jp), nstep(jp), ncount(jp), &
          deltax, aic, dbdr, corrfc_cmf, w0(:,jp), w1(:,jp), w0last, aj, ah, &
          ak, an, alo, opalk, slinetotk, slinek, add_corr, ak1, az1, &
          tauz1(:,jp), tau1(:,jp), pp, uoldk, voldk, v1oldk, hbnd, nbnd, &
          xlambda, counter_u, counter_uneg, counter_v)

!       save current u and v
        ubluwi(1:lmax, jp) = ucmf(1:lmax)
        vbluwi(1:lz, jp) = vcmf(1:lz)
        u_outer(jp, k) = ucmf(1)
!       JO     u_outer(jp,k)=ucmf(1)-0.5*iminusold(jp) !This is 0.5*Iplus
!       save current uold, vold, v1old for ALO
        uold(1:lmax, jp) = uoldk(1:lmax)
        vold(1:lz, jp) = voldk(1:lz)
        v1old(1:lz, jp) = v1oldk(1:lz)

      End Do jploop2


!     safety statment. Under very peculiar circumstances, j can be 0
      Do l = 1, nd
        If (aj(l)==0) Then
          aj(l) = xj_cmf(l, k-1)
          Print *, ' reset of aj=0 at', k, xlambda, l, ' to ', aj(l)
        End If
      End Do

!     Jo March 2025: should never happen
      Do l = 2, nd - 1
        If (aj(l)-alo(l)*slinetotk(l)<0.) Then
          Print *, k, xlambda, l, aj(l), alo(l), slinetotk(l)
          Write (999, *) ' stop: j-alo*sl < 0 after ray_complete'
          Stop ' j-alo*sl < 0 after ray_complete'
        End If
      End Do

!     save ALO
      If (maxval(alo)>1.D0+1.D-14 .Or. minval(alo)<0.) Then
        Print *, xlambda, maxval(alo), minval(alo)
        Do l = 1, nd
          Print *, aj(l), ' ', alo(l)
        End Do
        Write (999, *) ' stop: ALO_ray > 1 or ALO_ray= < 0'
        Stop ' ALO_ray > 1 or ALO_ray= < 0'
      End If
      alo_cmf(:, k) = alo
!     for tests
!     if(maxval(alo).le.0) print*,k,'max(alo)=0.'

!     JO check whether clumping correctly accounted for
!     save Eddington factors for moments equations
!     in the diffuse regime (large tau), fedd should reach 1/3 and gedd 3/5.
!     for geddp, this factor becomes depth and frequency dependent,
!     since it is given by 1/5 (dB_nu/dtau_nu)/B_nu  as a function of nu, T
!     and tau_nu
      fedd(:, k) = ak/aj !          grid
      Where (fedd(:,k)<0.05D0) fedd(:, k) = 0.05 ! for saftey: sphericality
!     factor
      gedd(:, k) = an(2:nd)/ah(2:nd) ! staggered grid
      hbound1(k) = ah(1)/aj(1)
      If (hbound1(k)<0.) hbound1(k) = 0.01D0 ! to avoid zero flux
      nbound1(k) = an(1)/aj(1)
      If (nbound1(k)<0.) nbound1(k) = 0.
      hboundnd(k) = hbnd/aj(nd)
      nboundnd(k) = nbnd/aj(nd)
!     geddp defined at staggered grid; thus r^2 J needs to be interpolated
      Do l = 1, nd - 1
        hplus = 0.5*(aj(l)*r(l)*r(l)+aj(l+1)*r(l+1)*r(l+1))/(0.5*(r(l)+r(l+ &
          1)))**2
        geddp(l, k) = an(l+1)/hplus
      End Do
!     this is the new variable required for moments_cmf
      add_corr_arr(:, k) = add_corr

      aj_ray = aj
      ah_ray = ah

    Else !                          not opt_ray

      If (xj_cmf(1,k)==0. .Or. xh_cmf(1,k)==0.) Then
        Write (999, *) ' stop: problems in xj_cmf/xh_cmf in moments iteration'
        Stop ' problems in xj_cmf/xh_cmf in moments iteration'
      End If
!     use previous solution
      If (xj_cmf(1,k)/=xj_old(1,k)) Then
        Write (999, *) ' stop: problems in xj_cmf vs. xj_old'
        Stop ' problems in xj_cmf vs. xj_old'
      End If
      aj_ray = xj_cmf(:, k)
      ah_ray = xh_cmf(:, k)
    End If optray !                 opt_ray


    hplus = 0.5*aic + corrfc_cmf*dbdr/(3.D0*opakk(nd))
!   remember: if inversion=.true., opacity and sourcefunction only continuum
!   values

!   In case (kmom_start > 1), we start later to avoid problems with
!   transition from approx. to exact method.
!   Later on; kmom_start might need to be replaced by k1, or something similar

!   if(k.gt.kmom_start) call
!   moments_cmf(k,nd,aj,aj_ray,ah,ah_ray,jold,hold,r,gr,gt,opakk,slinek,hplus,xlambda,counter_jneg,jneg)


!   working with g' (N/J) instead of g = N/H (see notes)
    If (k>kmom_start) Call moments_cmf1(k, nd, aj, aj_ray, ah, ah_ray, jold, &
      hold, r, gr, gt, opakk, slinek, hplus, xlambda, counter_jneg, jneg)
!   uncomment if moments not called
!   aj=aj_ray
!   ah=ah_ray

!   JO March 2025: adapt/renormalize (approx) alo for J(moments). Assume that
!   J_ray-alo(ray)*slinetotk = J_mom-alo(mom)*slinetotk
!   (if ray performed just before mom)
!   Would be (almost) exact if replacing "alo" (=Lambda*) by complete Lambda

!   in case of a previous mom-solution, alo is very approximate
!   (since slinek for previous alo is different from current slinek)

    If (opt_ray) Then
!     previous ray solution
      If (.Not. jneg) Then
!       no reset of J/H in mom
        Do l = 2, nd - 1
          alo_old = alo_cmf(l, k) ! as calculated in ray
          If (slinetotk(l)==0. .Or. alo_old==0.) Then
!           in certain cases, slinetotk=0 and alo_old ne 0
!           do nothing
            Continue
          Else If (aj(l)-alo_old*slinetotk(l)>0.) Then
!           do nothing, since in this case  alo = alo_old underestimates
!           the approximate value (which is allowed),
!           alo_mom(max)=J_mom/slinetotk
            Continue
          Else
!           in this case, J_mom-alo_old*slinetot < 0 (with alo_old = alo(ray))
!           and J_mom < J_ray, since J_ray-alo_old*slinetotk > 0
!           here, we have to renormalize ALO, resulting in alo_mom < alo_ray
            alo_cmf(l, k) = aj(l)/slinetotk(l)*0.99 ! safety
            If (alo_cmf(l,k)>=alo_old) Then
              Write (999, *) k, xlambda, l, aj(l), aj_ray(l), slinetotk(l), &
                alo_cmf(l, k), alo_old
              Write (999, *) ' stop: alo_mom > alo_ray after renormalization'
              Print *, k, xlambda, l, aj(l), aj_ray(l), slinetotk(l), &
                alo_cmf(l, k), alo_old
              Stop ' alo_mom > alo_ray after renormalization'
            End If
          End If
        End Do
      Else !                        jneg
!       reset of J/H to ray-values in mom, thus alo consistent (do nothing)
        Continue
      End If

    Else !                          .not.opt_ray
!     no current ray solution, previous one (with corresponding slinek
!     different
!     from current one) might be ray+moments or moments alone.
      If (.Not. jneg) Then
!       no reset of J/H in mom
        Do l = 2, nd - 1
          alo_old = alo_cmf(l, k) ! as calculated in previous iteration
          If (slinetotk(l)==0. .Or. alo_old==0.) Then
!           in certain cases, slinetotk=0 and alo_old ne 0
!           do nothing
            Continue
          Else If (aj(l)-alo_old*slinetotk(l)>0.) Then
!           do nothing (similar argument as above)
            Continue
          Else
!           renormalize ALO as above
            alo_cmf(l, k) = aj(l)/slinetotk(l)*0.99 ! safety
            If (alo_cmf(l,k)>=alo_old) Then
              Write (999, *) k, xlambda, l, aj(l), aj_ray(l), slinetotk(l), &
                alo_cmf(l, k), alo_old
              Write (999, *) &
                ' stop: alo_mom > alo_previous after renormalization'
              Print *, k, xlambda, l, aj(l), aj_ray(l), slinetotk(l), &
                alo_cmf(l, k), alo_old
              Stop ' alo_mom > alo_previous after renormalization'
            End If
          End If
        End Do
      Else !                        jneg
!       JO Oct. 2019: this statement is extremely important, since it ensures
!       that no inconsistent alo is used for problematic frequencies where J
!       and H
!       have been reset to the "old" ray-by-ray solution inside moments_cmf
        alo_cmf(:, k) = 0.
      End If
    End If

!   test current ALO
    Do l = 1, nd
      If (alo_cmf(l,k)>1.D0+1.D-14 .Or. alo_cmf(l,k)<0.) Then
        Write (999, *) k, xlambda, l, aj(l), alo_cmf(l, k), slinetotk(l), &
          aj_ray(l), alo_old
        Write (999, *) ' stop: ALO > 1 or ALO < 0 after renormalization'
        Print *, k, xlambda, l, aj(l), alo_cmf(l, k), slinetotk(l), aj_ray(l), &
          alo_old
        Stop ' ALO > 1 or ALO < 0 after renormalization'
      End If

      If (aj(l)-alo_cmf(l,k)*slinetotk(l)<0.) Then
        Write (999, *) k, xlambda, l, aj(l), alo_cmf(l, k), slinetotk(l), &
          aj_ray(l), alo_old
        Write (999, *) aj(l) - alo_cmf(l, k)*slinetotk(l), &
          aj_ray(l) - alo_old*slinetotk(l)
        Write (999, *) &
          ' stop: j-alo*sl < 0 after moments_cmf1 and renormalization'
        Print *, k, xlambda, l, aj(l), alo_cmf(l, k), slinetotk(l), aj_ray(l), &
          alo_old
        Print *, aj(l) - alo_cmf(l, k)*slinetotk(l), &
          aj_ray(l) - alo_old*slinetotk(l)
        Stop ' j-alo*sl < 0 after moments_cmf1 and renormalization'
      End If
    End Do

    errj = errj + aj/aj_ray
    errh = errh + ah/ah_ray

    jold = aj !                     for moments_cmf/next frequency point
    hold = ah !                     for moments_cmf/next frequency point


!   save J, H, H(ND+1)
!   after temp_converged, we use aj_ray to be consistent with ALO and
!   to allow for convergence
!   if(temp_converged) then
!   if(.not.opt_ray) stop ' temp_converged and .not. opt_ray in cmf_all'
!   xj_cmf(:,k)=aj_ray
!   else
!   otherwise
    xj_cmf(:, k) = aj !             note: always used in moments_cmf
!   at next frequency point via jold

!   endif
    xh_cmf(:, k) = ah !             used in moments_cmf at next frequency
!   point via hold
    xhinner_cmf(k) = ah(nd+1)


!   interpolate ah (staggered grid) to radial grid (including r^2)
    xh = ah*r2
    xhmin = minval(xh(2:nd))
    desh = 2.D0*abs(xhmin)

    Do l = 2, nd
      xh(l) = log10(xh(l)+desh)
    End Do

    Do l = 1, nd - 2
      xh(l+1) = q1(l)*xh(l+1) + q2(l)*xh(l+2)
      xh(l+1) = 10.D0**xh(l+1) - desh
    End Do
    xh(nd) = xh(nd+1)

!   now, xh = r^2 * H on mesh points (regular grid)

!   integrate various quantities over frequency within cmf-range

!   integration weights from using trapezoidal rule over delta nu, and
!   delta nu = c/lam(1-1/(1+v/c)) for first and last interval, and
!   delta_nu = (nu_(k-1)-nu_(k+1)) = c/lam(k-1)*(1-1/(1+v/c)^2)
!   for intermediate intervals
!   factor 0.5 (from trapezoidal rule) *clight * 1.d8 included in weights
!   note: 1/(1+v/c)=lam(1)/lam(2)
    If (k==1) Then
      weight = wnue1/lam_cmf(1)
      xh1 = xh(1:nd)
    Else If (k==nftot) Then
      weight = wnue1/lam_cmf(nftot-1)
      xhlast = xh(1:nd)
    Else
      weight = wnue2/lam_cmf(k-1)
    End If

    ahweight = ah*weight

!   rad. acc calculated on staggered grid
!   remember:use ah, ahweight without boundaries, i.e., from 2:ND
!   total rad_acc
!   JO check for clumping
    ydum = opakk/rho
    ydum1 = 0.5D0*(ydum(1:nd-1)+ydum(2:nd))
    rad_force = rad_force + ydum1*ahweight(2:nd)

!   continuum rad_acc
    ydum = opac/rho
    ydum1 = 0.5D0*(ydum(1:nd-1)+ydum(2:nd))
    cont_force = cont_force + ydum1*ahweight(2:nd)

!   thomson_acc
    ydum = sigth/rho
    ydum1 = 0.5D0*(ydum(1:nd-1)+ydum(2:nd))
    el_force = el_force + ydum1*ahweight(2:nd)

!   integrated flux on staggered grid
    ah_int = ah_int + ahweight

!   integrated flux (times r^2) on regular grid
    ydum = xh(1:nd)*weight
    xh_int = xh_int + ydum
    chibar_h = chibar_h + ydum*opakk
!   print*,k,lam_cmf(k),chibar_h(30),opakk(30),xh(30)

!   prepare calculation of opaross (with dummy tauross)
    Do l = 1, nd
      tauross(l) = dbdt(lam_cmf(k), temp(l))/opakk(l)
    End Do
    opaross = opaross + tauross*weight

!   lower flux in diffusion approx.
!   y1  = dbdt(lam_cmf(k),temp(nd))/opakk(nd)
    y1 = tauross(nd)
    fluxdif = fluxdif + y1*weight

!   check for frequencies where alo is zero at ALL depth points
    If (k>1. .And. sum(alo_cmf(:,k))==0.) Then
      If (opt_ray) Then
!       should not happen (alo = 0 only at certain, but not all depth points),
!       except for first cmf iteration (for a couple of first frequencies,
!       alo=0 because uold < 0 and thus set to zero!
        Print *, k, xlambda, ' opt_ray = T and alo = 0 at all depth points'
      End If
      counter_alo_zero = counter_alo_zero + 1
    End If

  End Do freloop
! end of frequency looop

  Write (999, *)
  Write (999, *) ' total number of cmf frequency points = ', nftot
  Write (999, *) ' number of resets J->J_ray = ', counter_jneg
  Write (999, *) ' number of calculated u= ', counter_u
  Print *
  Print *, ' total number of cmf frequency points = ', nftot
  Print *, ' number of resets J->J_ray = ', counter_jneg
  Print *, ' number of calculated u= ', counter_u
! print*,' number of resets u->u_new = ',counter_uneg ! no longer required
! currently, no reset of v
! print*,' number of resets abs(v)->u = ',counter_v
  Write (999, *)
  Write (999, *) ' opt_ray = ', opt_ray
  Print *
  Print *, ' opt_ray = ', opt_ray
! counter_alo_zero should be ge counter_jneg, but not signficantly
! (since for consecutive opt_ray=F the previous zeros are NOT reset)
  Print *, ' number of freq. points with ALO set to zero = ', counter_alo_zero
  Print *
  Print *
  Print *, &
    ' average ratio between moments (moments eq. vs. ray-by-ray solution)'
  Print *
  Print *, '  L   <Jnu(mom)/Jnu(ray)>  <Hnu(mom)/Hnu(ray)>'
  Do l = 1, nd
    Write (*, Fmt='(i3,2(7x,f10.5))') l, errj(l)/nftot, errh(l)/nftot
  End Do


  If (xj_old(1,1)/=0.) Then
    Print *
    Print *, ' max. change in xj_cmf as a function of depth'
    Do l = 1, nd
      maxchange_xj(l) = maxval(abs(1.-xj_cmf(l,:)/xj_old(l,:)))
      Print *, l, ' ', maxchange_xj(l)
    End Do
    Print *, ' overall max. change: ', maxval(maxchange_xj)
    Print *
    Write (999, *)
    Write (999, *) ' overall max. change in xj_xcmf (excluding ND=1) ', &
      maxval(maxchange_xj(2:nd))
    Write (999, *)
  End If

! for tests
! icounter=icounter+1
! write(counter,fmt='(i3)') icounter

! open(1,file=trim(modnam)//'/xj_comp_'//adjustl(counter))
! do k=1,nftot
! write(1,fmt='(i8,4(2x,G12.6))')
! k,lam_cmf(k),xj_cmf(31,k),sumopal_cmf(31,k),sumetal_cmf(31,k)
! enddo
! close(1)
! print*,' file xj_comp written to model'


  xj_old = xj_cmf !                 for next iteration

  start = .False.
! first-run condition ends here

! Test that approximate treatment and cmf transport are consistent
! ----------------------------------------------------------------------------
! (particularly regarding Iminus). Note that there cannot be a perfect
! consistency, since the approx. treatment neglects any freq. shifts.
! Thus, freq. ranges which are dominated by cont. processes should agree
! at *identical* frequencies, xobs = xcmf (i.e., optically thick continua), &
! whilst ranges dominated by photospheric processes should agree when
! xobs=xcmf/(1-v/c) (Lyman or Balmer edges).
! Additionally, there are numerical problems (formal solution for CMF,
! compared
! to moments solution for obsframe.)

! Thus, overall we require a 10% accuracy
! check first and last interval over 3 vinf
  k1 = 3.*vmax/vref

! close to blue edge
  summ = 0.
  summ1 = 0.
  nsum = 0

  Write (999, *)
  Print *
  If (cont_blue_thick) Then
    Write (999, *) ' blue cont. optically thick, no freq. shift in test'
    Print *, ' blue cont. optically thick, no freq. shift in test'
  End If

  Do k = 1, k1
    xlambda = 1.D8/lam_cmf(k) !     Kayser
    If (.Not. cont_blue_thick) xlambda = xlambda*(1.+vmax/clight) ! corrected
!   for shift
    Do l = ifre, 2, -1
      If (fre(l)>=xlambda .And. fre(l-1)<xlambda) Then
        xjapp = xj(1, l)*(xj(1,l)/xj(1,l-1))**((xlambda-fre(l))/(fre(l)-fre(l- &
          1)))
!       log interpol, since app. radiation field on coarse grid
!       print*,fre(l),xlambda,fre(l-1),xjapp,xj_cmf(1,k)
        nsum = nsum + 1
        xjapp = xj_cmf(1, k)/xjapp
        summ = summ + xjapp
        summ1 = summ1 + xjapp**2
        Exit
      End If
    End Do
  End Do

  mean = summ/nsum
  dev = sqrt(summ1/nsum-mean**2)

  Write (999, Fmt= &
    '(''  close to blue edge: <j_cmf(1)/j_app(1)> = '',f7.2,'' +/- '',f8.3)') &
    mean, dev
  Write (*, Fmt= &
    '(''  close to blue edge: <j_cmf(1)/j_app(1)> = '',f7.2,'' +/- '',f8.3)') &
    mean, dev

  If ((mean-3.*dev)>1.1 .Or. (mean+3.*dev)<0.9) Then
    Write (999, *) ' CMF AND APPROX. TRANSFER NOT CONSISTENT AT BLUE EDGE!!!!'
    Write (999, *) ' CMF AND APPROX. TRANSFER NOT CONSISTENT AT BLUE EDGE!!!!'
    Write (999, *) ' CMF AND APPROX. TRANSFER NOT CONSISTENT AT BLUE EDGE!!!!'
    Write (999, *) &
      ' should happen only when continuum becomes optically thick around L=5!'
    Write (999, *) ' check FLUXCONT!!!'
    Print *, ' CMF AND APPROX. TRANSFER NOT CONSISTENT AT BLUE EDGE!!!!'
    Print *, ' CMF AND APPROX. TRANSFER NOT CONSISTENT AT BLUE EDGE!!!!'
    Print *, ' CMF AND APPROX. TRANSFER NOT CONSISTENT AT BLUE EDGE!!!!'
    Print *, &
      ' should happen only when continuum becomes optically thick around L=5!'
    Print *, ' check FLUXCONT!!!'
!   stop ' CMF and approx. transfer not consistent at blue edge'
  End If

! close to red edge
  summ = 0.
  summ1 = 0.
  nsum = 0
  Do k = nftot, nftot - k1, -1
    xlambda = 1.D8/lam_cmf(k)*(1.+vmax/clight) ! Kayser, always corrected for
!   shift
    Do l = 2, ifre
      If (fre(l)>=xlambda .And. fre(l-1)<xlambda) Then
!       here we can afford a linear interpolation (Rayleigh Jeans)
        xjapp = xj(1, l) + (xj(1,l)-xj(1,l-1))/(fre(l)-fre(l-1))*(xlambda-fre( &
          l))
!       print*,fre(l),xlambda,fre(l-1),xjapp,xj_cmf(1,k)
        nsum = nsum + 1
        xjapp = xj_cmf(1, k)/xjapp
        summ = summ + xjapp
        summ1 = summ1 + xjapp**2
        Exit
      End If
    End Do
  End Do
  mean = summ/nsum
  dev = sqrt(summ1/nsum-mean**2)
  Write (999, Fmt= &
    '(''  close to  red edge: <j_cmf(1)/j_app(1)> = '',f7.2,'' +/- '',f8.3)') &
    mean, dev
  Write (999, *)
  Write (*, Fmt= &
    '(''  close to  red edge: <j_cmf(1)/j_app(1)> = '',f7.2,'' +/- '',f8.3)') &
    mean, dev
  Print *


! JO comment out for tests
  If ((mean-3.*dev)>1.1 .Or. (mean+3.*dev)<0.9) Print *, &
    ' CMF and approx. transfer not consistent at red edge!!!!'
! & stop ' CMF and approx. transfer not consistent at red edge'

! ----------------------------------------------------------------------------

! remember: in the following, xh is xh*r^2 on regular grid, contrasted to ah

! save xj(cmf), xj(approx), xh(nd+1) at specific radii
  xhinner_cmf = xhinner_cmf*r2(nd)

  If (optout_xj_xh) Then

    icounter = icounter + 1
    Write (counter, Fmt='(i3)') icounter

!   open(1,file=trim(modnam)//'/out_xj_cmf_'//adjustl(counter))
    Open (1, File=trim(modnam)//'/out_xj_cmf')
    Do k = 1, nftot
      Write (1, 160) lam_cmf(k), xj_cmf(1, k), xj_cmf(43, k), xj_cmf(nd, k), &
        xhinner_cmf(k)
!     write(1,110)
!     lam_cmf(k),xj_cmf(9,k),xj_cmf(35,k),xj_cmf(45,k),xhinner_cmf(k)
    End Do
    Close (1)

!   open(1,file=trim(modnam)//'/out_xj_app_'//adjustl(counter))
    Open (1, File=trim(modnam)//'/out_xj_app')
    Do k = 1, ifre
      Write (1, 160) 1.D8/fre(k), xj(1, k), xj(43, k), xj(nd, k), xxh(nd, k)
!     write(1,110) 1.d8/fre(k),xj(9,k),xj(43,k),xj(nd,k),xxh(nd,k)
    End Do
    Close (1)

!   after tests, comment out next two blocks
    Open (1, File=trim(modnam)//'/out_xj_saved')
    Do k = 1, ifre
      Write (1, 160) 1.D8/fre(k), xj_save(1, k), xj_save(43, k), &
        xj_save(nd, k), xxh(nd, k)
    End Do
    Close (1)

    Open (1, File=trim(modnam)//'/out_xj_smoothed')
    Do k = 1, ifre
      Write (1, 160) 1.D8/fre(k), xj_smoothed(1, k), xj_smoothed(43, k), &
        xj_smoothed(nd, k), xxh(nd, k)
    End Do
    Close (1)

!   if xh_cmf has been calculated
!   open(1,file=trim(modnam)//'/out_xh_cmf')
!   do k=1,nftot
!   write(1,110)
!   lam_cmf(k),xh_cmf(1,k),xh_cmf(43,k),xh_cmf(nd,k),xhinner_cmf(k)
!   enddo
!   close(1)

!   open(1,file=trim(modnam)//'/out_xh_app')
!   do k=1,ifre
!   call diffus (1.d8/fre(k), temp, r, nd, aic, dbdr)
!   write(1,110) 1.d8/fre(k),xxh(1,k),xxh(43,k),xxh(nd,k), &
!   dbdr/(3.d0*opac_coarse(nd,k)*sr)*corrfc
!   enddo
!   close(1)
  End If

! calculate and save rad_forces (within wavblue ... wavred)
  xfac = (4.*pi/clight)/sr

  rad_force = xfac*rad_force
  cont_force = xfac*cont_force
  el_force = xfac*el_force*sr

  chibar_h = chibar_h/sr !          since opakk includes SR

  Print *
  Print *, ' L  grad(stag) grad(grid) grad(nom_flux) chibar_h(cmf-range)'
  Open (1, File=trim(modnam)//'/out_rad_force')
  Do l = 1, nd - 1
!   now on approximate staggered grid:
!   accelerations from 1:ND-1, fluxes (ah*r^2) from 2:ND
    rdum = (r(l)+r(l+1))*0.5
    vdum = velo(l) + (velo(l+1)-velo(l))/(r(l+1)-r(l))*(rdum-r(l))
    ah_int(l+1) = ah_int(l+1)*r2(l+1)
    Write (1, 140) rdum, vdum, rad_force(l), cont_force(l), el_force(l), &
      ah_int(l+1)
!   check consistency:
!   usually, grad(grid)(l) should lie in between grad(stag)(l-1) and
!   grad(stag)(l)
!   differences between grad(grid) and grad(nom_flux) when flux different from
!   nom. flux
!   NOTE: int Hnu dnu here only for CMF range. In lower atmosphere of hot
!   stars,
!   this can be significantly smaller that nominal flux, since blue range
!   (approx. method)
!   might contribute a lot. Thus: H_nom > H_range, and grad(nom_flux) >
!   grad(grid)
    rdum = sigsb*teff**4/clight*(rtau23/r(l))**2
    vdum = 4.*pi/(clight*r(l)*r(l))
    chibar_h_cmf(l) = chibar_h(l)/(xh_int(l)*rho(l))
!   new column chibar_h_cmf included
    Write (*, 150) l, rad_force(l), chibar_h(l)/rho(l)*vdum, &
      chibar_h_cmf(l)*rdum, chibar_h_cmf(l)
  End Do
  chibar_h_cmf(nd) = chibar_h(nd)/(xh_int(nd)*rho(nd))
  Close (1)
  Print *

! outer boundary (so far, ah_int no longer used from here on)
  ah_int(1) = ah_int(1)*r2(1)
  ah_int(nd+1) = ah_int(nd+1)*r2(nd+1)

! now combine freq. ranges: outside wavblue ... wavred, use approx. values
! most quantities have been already calculated above
! xxh is flux*r^2 from approx. solution on regular grid

! integrate over blue range
  xh_int_blue = 0.D0
  oparblue = 0.D0

  Do k = kblue + 1, ifre
    weight = wfre(k)
    xlambda = 1.D8/fre(k)
    If (k==kblue+1) weight = wlast_blue
    xh_int_blue = xh_int_blue + xxh(:, k)*weight
    Do l = 1, nd
      tauross(l) = dbdt(xlambda, temp(l))/opac_coarse(l, k)
    End Do
    oparblue = oparblue + tauross*weight
    chibar_h = chibar_h + xxh(:, k)*opac_coarse(:, k)*weight
!   print*,'blue',k,chibar_h(30),opac_coarse(30,k),xxh(30,k)
  End Do

  weight = clight*(fre(kblue+1)-1.D8/lam_cmf(1))
  If (weight<=0.) Then
    Write (999, *) ' stop: weight < 0 (3)'
    Stop ' weight < 0 (3)'
  End If
  xh_int_blue = xh_int_blue + 0.5D0*(xxh(:,kblue+1)+xh1)*weight
  xlambda = 1.D8/fre(kblue+1)
  Do l = 1, nd
    tauross(l) = dbdt(xlambda, temp(l))/opac_coarse(l, kblue+1)
  End Do
! assuming that opac_coarse is not too different from opakk at lam_cmf(1)
  oparblue = oparblue + tauross*weight
  chibar_h = chibar_h + 0.5D0*(xxh(:,kblue+1)+xh1)*opac_coarse(:, kblue+1)* &
    weight


! integrate over red range
  xh_int_red = 0.D0
  oparred = 0.D0

  Do k = 1, kred - 1
    weight = wfre(k)
    xlambda = 1.D8/fre(k)
    If (k==kred-1) weight = wlast_red
    xh_int_red = xh_int_red + xxh(:, k)*weight
    Do l = 1, nd
      tauross(l) = dbdt(xlambda, temp(l))/opac_coarse(l, k)
    End Do
    oparred = oparred + tauross*weight
    chibar_h = chibar_h + xxh(:, k)*opac_coarse(:, k)*weight
!   print*,'red',k,chibar_h(30),opac_coarse(30,k),xxh(30,k)
  End Do

  weight = clight*(1.D8/lam_cmf(nftot)-fre(kred-1))
  If (weight<=0.) Then
    Write (999, *) ' stop: weight < 0 (3)'
    Stop ' weight < 0 (3)'
  End If
  xh_int_red = xh_int_red + 0.5D0*(xxh(:,kred-1)+xhlast)*weight
  xlambda = 1.D8/fre(kred-1)
  Do l = 1, nd
    tauross(l) = dbdt(xlambda, temp(l))/opac_coarse(l, kred-1)
  End Do
! assuming that opac_coarse is not too different from opakk at lam_cmf(nftot)
  oparred = oparred + tauross*weight
  chibar_h = chibar_h + 0.5D0*(xxh(:,kred-1)+xhlast)*opac_coarse(:, kred-1)* &
    weight

! calculating total opaross and corresponding tauross
! SR included in opakk, but not in opac_coars
  opaross = opaross + (oparblue+oparred)/sr
  tauross = 7.21863D-5*temp**3
  opaross = tauross/opaross

  summ = 0.D0
  tauross(1) = 0.D0
  Do l = 1, nd - 1
    rdum = r(l) - r(l+1)
    summ = summ + rdum*0.5D0*(opaross(l)+opaross(l+1))
    tauross(l+1) = summ
!   print*,l+1,opaross(l+1)
  End Do

  summ = 0.D0
  ferr(1) = 0.D0

! for tests: assuming opaross a power-law in r, we obtain the following result
! note : similar result when assuming opaross exp(-(r-r0)/H)
! in both cases, the derived tauross is a bit *lower* than derived from the
! trapezoidal rule. Since we overestimate the temperature, this does
! not help. We would need a larger tauross (exact) vs. tauross(trapez)
! do l = 1,nd - 1
! rdum = log(opaross(l)/opaross(l+1))/log(r(l)/r(l+1))
! summ = summ + (opaross(l)*r(l)-opaross(l+1)*r(l+1))/(rdum+1.d0)
! ferr(l+1) = summ
! print*,l+1,tauross(l+1),ferr(l+1)
! end do


! recheck and compare with 'approximate' tauross (ferr used as dummy)
  oparred = 0.D0
  Do k = 1, ifre
    weight = wfre(k)
    xlambda = 1.D8/fre(k)
    Do l = 1, nd
      ferr(l) = dbdt(xlambda, temp(l))/opac_coarse(l, k)
    End Do
    oparred = oparred + ferr*weight
  End Do
  ferr = 7.21863D-5*temp**3
  oparred = ferr/oparred*sr

  summ = 0.D0
  ferr(1) = 0.D0
  Do l = 1, nd - 1
    rdum = r(l) - r(l+1)
    summ = summ + rdum*0.5D0*(oparred(l)+oparred(l+1))
!   print*,l+1,oparred(l+1)
    ferr(l+1) = summ
  End Do

  Print *, 'L     taur(CMF)   taur(OPAC) taur(used)'
  Do l = 2, nd
    Write (*, Fmt='(i3,3(2x,f10.4))') l, log10(tauross(l)), log10(ferr(l)), &
      log10(taur(l))
  End Do
  Print *
! JO: so far, it seems that the flux-corrected CMF temperature is too large in
! photospheric regions.
! In any case we need to use the tauross as calculated here (helps a bit),
! but this does not completely solve the problem.
! improving the radial resolution helps, but the grid at higher tau is still
! too coarse (equidistant in LOG(M) and thus LOG(TAUR)
! setting XM(ND) = 1.1 XM(ND-1) also helps, but still:
! -> somewhat erroneous flux -> wrong T(r), though it is conserved
! REAL solution only if moments' equation were solved.

! in new version, this is done (successfully)

  taur = tauross
  Open (1, File=trim(modnam)//'/TAU_ROS', Status='unknown', Form='formatted')
  Rewind 1
  Do l = 1, nd
    Read (1, *) rdum, ferr(l) !     ferr is xne_lte
  End Do
  Rewind 1
  Do l = 1, nd
    Write (1, *) taur(l), ferr(l)
  End Do
  Close (1)

  Write (999, *) ' TAUROSS UPDATED!!!'
  Write (999, *)
  Print *, ' TAUROSS UPDATED!!!'
  Print *

! Perform various checks
! re-check diffusion approx.
  const = corrfc_cmf*dtdr/3.D0 !    SR included in opacities
  fluxdif = fluxdif*const
  errmax = 1. - xh_int(nd)/fluxdif
  Print *, ' Deviation actual flux/diff approx. in CMF-region = ', errmax
  Print *

  fluxdif = fluxdif + fldiff_blue + fldiff_red

  errmax = 1.D0 - fluxdif/nom_flux
  Print *, ' Error in diffusion approx. = ', errmax
  Print *

! if(abs(errmax).gt.1.d-15) stop ' error in diff. approx. > accuracy
! (cmf_complete)'
! allow for a bit larger differences, due to different compilers
  If (abs(errmax)>1.D-14) Then
    Write (999, *) ' stop: error in diff. approx. > accuracy (cmf_complete)'
    Stop ' error in diff. approx. > accuracy (cmf_complete)'
  End If

  Print *, ' L  H_int(blue)   H_int(red)    H_int(CMF)    H_tot/H_nom'
  Do l = 1, nd
    Write (*, 180) l, xh_int_blue(l)/nom_flux, xh_int_red(l)/nom_flux, &
      xh_int(l)/nom_flux, (xh_int(l)+xh_int_blue(l)+xh_int_red(l))/nom_flux
  End Do
  Print *

  xh_int = xh_int + xh_int_blue + xh_int_red

! now calculate flux-errors: again, w.r.t. nominal flux at lower boundary,
! which is larger than the flux at Rstar (cf. ROUTINE CONT)

  errmax = 1.D0 - xh_int(nd)/nom_flux
  Write (999, *) ' Deviation actual flux/nominal flux (lower boundary) = ', &
    errmax
  Write (999, *)
  Print *, ' Deviation actual flux/nominal flux (lower boundary) = ', errmax
  Print *

! check flux-conservation/overall flux error
  ferr = xh_int/nom_flux - 1.D0
  errmax = maxval(abs(ferr))
  Write (999, *) ' Max. absolute flux error = ', errmax
  Write (999, *)
  Print *, ' Max. absolute flux error = ', errmax
  Print *

! calculate flux-weighted opacity
  chibar_h = chibar_h/(xh_int*rho)
  Open (61, File=trim(modnam)//'/CHIBAR_H_CMF', Status='unknown')
  Rewind 61
  Do l = 1, nd
!   here was a bug, cured Sept 12 2019
!   write(61,*) l,chibar_h(l),chibar_h(l)*sigsb*teff**4/clight
    If (abs(1.-chibar_h_cmf(l)/chibar_h(l))>0.1) Then
      Write (999, *) ' WARNING WARNING: chibar_h in CMF and complete &
        &range significantly different:', l, ' ', chibar_h(l), ' ', &
        chibar_h_cmf(l)
      Print *, ' WARNING WARNING: chibar_h in CMF and complete &
        &range significantly different:', l, ' ', chibar_h(l), ' ', &
        chibar_h_cmf(l)
    End If
    Write (61, *) l, chibar_h(l), chibar_h(l)*sigsb*teff**4/clight*(rtau23/r(l &
      ))**2
  End Do
  Close (61)
  Write (999, *)
  Write (999, *) ' file CHIBAR_H_CMF written'
  Write (999, *)

  Print *
  Print *, ' file CHIBAR_H_CMF written'
  Print *

! recheck flux from approx. method (complete freq. range)
  xh_int = 0.D0
  fluxdif = 0.D0

  Do k = 1, ifre
    xh_int = xh_int + xxh(:, k)*wfre(k)
    Call diffus(1.D8/fre(k), temp, r, nd, aic, dbdr)
    fluxdif = fluxdif + dbdr/(3.D0*opac_coarse(nd,k)*sr)*corrfc*wfre(k)
  End Do

  errmax = maxval(abs(1.D0-xh_int/nom_flux))
  Write (999, *) ' Max. absolute flux error (approx. method) = ', errmax
  Write (999, *)
  Print *, ' Max. absolute flux error (approx. method) = ', errmax
  Print *
  errmax = abs(1.D0-fluxdif/nom_flux)
  Print *, ' Error in diff. approx. (approx. method) = ', errmax
  Print *

! temperature correction with respect to flux conservation
  Call tcorr(taur, ferr, teff, temp, r, nd, dtfcorr, ndiv_calc+1)

! for tests
! do l=1,nd
! print*,l,log10(taur(l)),temp(l),xh_int(l),dtfcorr(l)
! enddo
  Write (999, *) ' CMF_COMPLETE DONE'
  Write (999, *)
  Print *, ' CMF_COMPLETE DONE'
  Print *

! call test !calculating specific jbars for some tests
  Return

130 Format (I2, 2X, 6(E9.3,2X))
140 Format (6E15.5)
150 Format (I3, 4(2X,G10.4))
160 Format (F14.4, 4(2X,E10.4))
170 Format (F14.4, 6(2X,E10.4))
180 Format (I3, 4(2X,E12.6))
190 Format (I3, 6(2X,E12.6))

End Subroutine

!-----------------------------------------------------------------------

Subroutine ray_complete(k, jp, z, r, nd, np, u, v, opa1, sl1, taus_outer, &
  iminusold, in_line_old, nstep, ncount, deltax, aic, dbdr, corrfc, w0, w1, &
  w0last, aj, ah, ak, an, alo, opalk, slinetotk, slinek, add_corr, aakblue, &
  aaz1, tauz1, tau1, pp, uold, vold, v1old, hbnd, nbnd, xlambda, counter_u, &
  counter_uneg, counter_v)
! subroutine ray from nlte, adapted for cmf_complete, for freq. k and
! impact parameter k

! Do no longer use I=I+=I- at last ray; leads to problems (spikes at red
! profile edge) due to transition to continuum, and gives almost no
! improvement

! JO: in this new version, u is no longer set to zero for neg. values.
! uold still needs to be set to zero

! after hundreds of tests, we came to the enclosed solution to avoid negative
! u's. For consistency, a new variable add_corr is calculated to be used within
! the moments equations.  
  
  Use :: nlte_type
  Use :: nlte_dim, Only: id_ndept
  Use :: nlte_var, Only: vp2, vp3

  Implicit None


! line radiation transfer in the comoving frame from a
! given source function for a given impact parameter jp.
! the integration is carried out in space and freq. to
! yield the mean line intensity

! new version: inlinen aller algebraischen routinen

! calculated are Jnu and consistent ALO

! opalk, slinetotk: total line opacities and source functions
! slinek: total source function (incl. cont)
! total opacity not needed (since provided via aakblue etc.)

! .. parameters ..
  Integer (i4b), Parameter :: nd1 = id_ndept
! ..
! .. scalar arguments ..
  Real (dp), Intent (In) :: aic, corrfc, dbdr, deltax, w0last, xlambda, opa1, &
    taus_outer, sl1
  Integer (i4b), Intent (In) :: jp, nd, np
! ..
! .. array arguments ..
  Real (dp), Intent (Inout) :: aj(nd), ak(nd), ah(nd+1), an(nd+1), alo(nd), &
    u(nd), v(nd-1), uold(nd), vold(nd-1), v1old(nd-1), add_corr(nd)

  Real (dp), Intent (In) :: pp(nd), r(nd), opalk(nd), slinetotk(nd), &
    slinek(nd), aakblue(nd), aaz1(nd-1), tau1(nd), tauz1(nd), w0(nd), z(nd), &
    w1(nd+1)

  Real (dp), Intent (Inout) :: iminusold, hbnd, nbnd
  Integer (i4b), Intent (Inout) :: ncount, counter_u, counter_uneg, counter_v
  Integer (i4b), Intent (In) :: nstep
  Logical, Intent (Inout) :: in_line_old


! ..
! .. local scalars ..
  Logical :: test

  Real (dp) :: a, akblue, az1, b, daz, dazm, dbz, dbzm, dt, dtz, dtzm, dx, &
    dxz, dxz1, dxzm, dxzm1, h1, help, pw, rmax, tauc, taus, s1, tauz, tb1, &
    ff1, iplus, uo, aux

  Integer (i4b) :: k, l, l1, lmax, lz, ni
! ..
! .. local arrays ..
  Real (dp) :: ga(nd1), h(nd1), qq(nd1), s(nd1), ta(nd1), tb(nd1), tc(nd1), &
    ub(nd1), va(nd1), vb(nd1), uorig(nd1), rhs(nd1), tc1(nd1), delta_u(nd1), &
    delta_b(nd1)

! variables for alo
  Real (dp) :: ff(nd1), ddia(0:nd1), edia(nd1+1)

! ..
! .. intrinsic functions ..
! ..

  lmax = min0(np+1-jp, nd)

  If (k==1) Then
    Do l = 1, lmax
      pw = w0(l)
      aj(l) = aj(l) + u(l)*pw
      ak(l) = ak(l) + u(l)*vp2(l, jp)
!     no action for alo required, since uold=0, and alo=0 set for k=1
    End Do

    Do l = 1, lmax - 1
      pw = w1(l+1)
      ah(l+1) = ah(l+1) + v(l)*pw
      an(l+1) = an(l+1) + v(l)*vp3(l+1, jp)
    End Do
!   outer and inner boundary, approximating Iminus(out)=0.
!   ah(1)=ah(1)+u(1)*w1(1)
!   an(1)=an(1)+u(1)*vp3(1,jp)
!   JO June 2025: now with iminus_cont
    ah(1) = ah(1) + (u(1)-iminusold)*w1(1)
    an(1) = an(1) + (u(1)-iminusold)*vp3(1, jp)
    If (lmax==nd) Then
      iplus = aic + z(nd)*aakblue(nd)*dbdr*corrfc
      ah(nd+1) = ah(nd+1) + (iplus-u(nd))*w1(nd+1)
      hbnd = hbnd + u(nd)*w1(nd+1)
      nbnd = nbnd + u(nd)*vp3(nd+1, jp)
    End If

    add_corr(1:lmax) = 0. !         no correction

    Return
  End If

  lz = lmax - 1
  rmax = r(1)

! uold, vold, v1old initalized outside

! approximate optical depths for outer boundary condition
  tauc = opa1*rmax/3.D0 !           optically thick continuum only for rho^2
! processes
  taus = taus_outer/(pp(1)*deltax) ! since pp divided by deltax


! ***  inlined subroutine cmfset
! ----------------------------------------------------------------------

! ***  outer boundary condition  -  1st order with Iminus

  akblue = aakblue(1)
  az1 = aaz1(1)
  dx = pp(1)*akblue

! JO modified boundary condition; might need to be changed in subr. ray
! as well, but for line-cores each reasonable condition is sufficient

! s1 is source term for outer boundary, corresponds to Imin - P/kappa dImin/dx
! the inherent philosophie and physics for calculating I_minus can be found
! in Puls+ 2020, Appendix  
  Call calc_iminus(s1, xlambda, dx, taus, tauc, sl1, slinek(1), iminusold, &
    in_line_old, nstep, ncount)

! set up of coefficients, starting with outer boundary condition

  s(1) = s1
! s(1)=0. !if not boundary condition for Iminus
! JO changed July 2025: delta u/delta z/kappa_0
! old version
  tc(1) = tauz1(1)*akblue

! might be replaced by delta u/delta z/kappa_1/2, to allow for a more
! precise calculation (e.g., of Eddington factors at outermost
! region). The new formulation results from an expansion of du/dz.

! new version
! tc(1) = tauz1(1)*az1

! after a multitude of tests. it turned out that the old formulation works better
  
! the terms propto dv/dx = d(u-Iminus)/dx ("dx" and dIminus/dx
! included in s1) are still calculated with forefactor P/kappa_0,
! since they refer to the conditions at the outermost point.
! Unfortunatel, the new version does not work (at least) in thick B-star
! winds, thus old formulation still used.
  tb(1) = tc(1) + dx + 1.D0
  ub(1) = dx
  vb(1) = 0.D0
  dtzm = tauz1(1)*az1

! ***  alo for outer boundary

  ff1 = 0.D0
! for first tests, set ff1 to zero
! JO in case, define ALO for optically thick boundary condition
  ff(1) = ff1
! !ALO        if (s(1).eq.sc) ff(1) = 0.d0

! ***  for g and h, the matrix element s are not different from
! ***  inner points

  dxzm = .5D0*(pp(1)+pp(2))*az1
  dxzm1 = 1.D0/(1.D0+dxzm)
  dazm = dtzm*dxzm1
  dbzm = dxzm*dxzm1
  ga(1) = -dazm
  h(1) = dbzm

! ***  non-boundary points

  Do l = 2, lz
    akblue = aakblue(l)
    az1 = aaz1(l)
    s(l) = slinek(l)
    dt = tau1(l)*akblue
    dtz = tauz1(l)*az1
    dx = pp(l)*akblue
    dxz = .5D0*(pp(l)+pp(l+1))*az1
    dxz1 = 1.D0/(1.D0+dxz)
    daz = dtz*dxz1
    dbz = dxz*dxz1
    ta(l) = dt*dazm
    tc(l) = dt*daz
    tb(l) = ta(l) + tc(l) + dx + 1.D0
    ff(l) = opalk(l)*akblue
    ub(l) = dx
    va(l) = -dt*dbzm
    vb(l) = dt*dbz
    ga(l) = -daz
    h(l) = dbz
    dazm = daz
    dbzm = dbz
  End Do

  l = lmax

  If (lmax==nd) Then

!   ***  inner boundary condition (core rays)  -  only to first order
!   ***  diffusion-approximation

    akblue = aakblue(l)
    s(l) = aic + z(l)*akblue*dbdr*corrfc
    iplus = s(l)
    dx = pp(l)*akblue
!   dt = tauz1(l-1)*akblue
!   CHANGED BY JP Sept. 2013, to be consistent with CONT1
!   (if (u(l) linear, then du/dtau = delta_u/delta_tau,
!   and delta_tau = chi_half*dz)
!   this is inconsistent with the dx-term, but the latter does not matter
!   anyway
    dt = tauz1(l-1)*aaz1(l-1)
    ta(l) = dt
    tb(l) = dt + dx + 1.D0
    ff(l) = 0.D0
    ub(l) = dx
    va(l) = 0.0D0
    vb(l) = 0.0D0
  Else

!   ***  inner boundary condition (non-core rays)  -  second order

    akblue = aakblue(l)
    s(l) = slinek(l)
    tauz = z(lz)
    dt = akblue/tauz
    dx = pp(l)*akblue
    ta(l) = 2.D0*dt*dazm !          instead of daz
    tb(l) = ta(l) + dx + 1.D0
    ub(l) = dx
    va(l) = -2.D0*dt*dbzm !         instead of dbz
    ff(l) = opalk(l)*akblue
  End If

! -----------------------------------------------------------------------
! calculation of alo (see Puls 1993)
! ***  continuation of subroutine ray

  qq(1) = 0.D0

  Do l = 2, lz
    qq(l) = va(l)*v1old(l-1) + vb(l)*vold(l)
  End Do

  qq(lmax) = va(lmax)*v1old(lz)

! uold=ub * uold + qq + ff

  Do l = 1, lmax
    uold(l) = ub(l)*uold(l) + qq(l) + ff(l)
  End Do

! ***     now, uold is the solution vector for the inhomgeneous
! ***     equation corresponding to u(sl=1)-u(sl=0)

! inlinen von diag mit option 1


  ddia(0) = 0.D0
  ta(1) = 0.D0
  tc(lmax) = 0.D0

  Do l = 1, lmax
    ddia(l) = tc(l)/(tb(l)-ta(l)*ddia(l-1))
  End Do

  edia(lmax+1) = 0.D0

  Do l = lmax, 1, -1
    edia(l) = ta(l)/(tb(l)-tc(l)*edia(l+1))
  End Do

  Do l = 1, lmax
    uold(l) = uold(l)/((1.D0-ddia(l)*edia(l+1))*(tb(l)-ta(l)*ddia(l-1)))
    uold(l) = max(uold(l), 0.D0)
  End Do

! ***  now,uold is the solution corresponding to u(sl=1)-u(0)

  Do l = 2, lz
    l1 = l - 1
    If (uold(l)==0.) Then
      v1old(l1) = 0.
      vold(l) = 0.
    Else
      v1old(l1) = -ga(l1)*uold(l) + h(l1)*v1old(l1)
      vold(l) = ga(l)*uold(l) + h(l)*vold(l)
    End If
  End Do

  l = lz
  v1old(l) = -ga(l)*uold(lmax) + h(l)*v1old(l)

! ***  now, v1old, vold are the v's corresponding to delta u

! inlinen von vmalv

! now, calculation of J and H
  b = v(1)
  qq(1) = vb(1)*b !                 = 0, since vb(1)=0.
  If (qq(1)/=0.) Then
    Write (999, *) ' stop: qq(1) ne 0 in ray_complete!'
    Stop ' qq(1) ne 0 in ray_complete!'
  End If

  Do l = 2, lz
    a = b
    b = v(l)
    qq(l) = va(l)*a + vb(l)*b
  End Do

  qq(lmax) = va(lmax)*b

! check for diagonal dominance (can be commented out after confirmation)
! JO March 2025
! commented out, since confirmed by zillions of runs
! if(ta(1).ne.0. .or. tc(lmax).ne.0.) stop ' ta or tc ne 0 in ray'
! do l=1,lmax
! if(ta(l).lt.0. .or. tb(l).lt.0. .or. tc(l).lt.0.) then
! print*,'not positive in ray'
! print*,xlambda,jp,l,ta(l),tb(l),tc(l)
! stop ' not positive in ray_complete'
! endif
! if(abs(ta(l))+abs(tc(l)).gt.abs(tb(l))) then
! print*,'not diagonally dominated in ray'
! print*,xlambda,jp,l,ta(l),tb(l),tc(l)
! endif
! enddo

! u = ub * u + qq + s

! JP: new hack March 2019; avoid problems with Iminus, set to zero


  test = .False.

  uorig(1:lmax) = u(1:lmax)

  Go To 100
! -------------------------------------------------------------------
! in case, original version (requires that add_corr remains zero)

  Do l = 1, lmax
    u(l) = ub(l)*u(l) + qq(l) + s(l)
  End Do
  If (u(1)<0) Then
    If (s(1)>=0) Then
      Write (999, *) ' stop: rhs(1) < 0 and s(1) ge 0!'
      Stop ' rhs(1) < 0 and s(1) ge 0!'
    End If
    u(1) = ub(1)*uorig(1)
    test = .True.
  End If

  rhs(1:lmax) = u(1:lmax)
  tc1(1:lmax) = tc(1:lmax)

! ***     inlining of invtri

  tb1 = 1.D0/tb(1)
  tc(1) = tc(1)*tb1
  u(1) = u(1)*tb1

  Do l = 2, lz
    help = tb(l) - ta(l)*tc(l-1)
    h1 = 1.D0/help
    tc(l) = tc(l)*h1
    u(l) = (u(l)+u(l-1)*ta(l))*h1
  End Do

  u(lmax) = (u(lmax)+u(lz)*ta(lmax))/(tb(lmax)-tc(lz)*ta(lmax))

  If (u(lmax)<0.) u(lmax) = 0.

  Do l = 1, lz
    ni = lmax - l
!   u(ni) = max(u(ni),-tc(ni)*u(ni+1)) ! solution vector adapted
    aux = tc(ni)*u(ni+1)
    u(ni) = u(ni) + aux
    If (u(ni)<0.) u(ni) = 0.
  End Do

  Go To 110
! -------------------------------------------------------------------
! new version

100 Do l = 1, lmax
    u(l) = ub(l)*u(l) + qq(l) + s(l)
  End Do
  If (u(1)<0) Then
    If (s(1)>=0) Then
      Write (999, *) ' stop: rhs(1) < 0 and s(1) ge 0!'
      Stop ' rhs(1) < 0 and s(1) ge 0!'
    End If
    u(1) = ub(1)*uorig(1)
    test = .True.
  End If

  rhs(1:lmax) = u(1:lmax)
  tc1(1:lmax) = tc(1:lmax)


! ***     inlining of invtri

  tb1 = 1.D0/tb(1)
  tc(1) = tc(1)*tb1
  u(1) = u(1)*tb1

  Do l = 2, lz
    help = tb(l) - ta(l)*tc(l-1)
    h1 = 1.D0/help
    tc(l) = tc(l)*h1
    u(l) = (u(l)+u(l-1)*ta(l))*h1
  End Do

  u(lmax) = (u(lmax)+u(lz)*ta(lmax))/(tb(lmax)-tc(lz)*ta(lmax))

  Do l = 1, lz
    ni = lmax - l
!   u(ni) = max(u(ni),-tc(ni)*u(ni+1)) ! solution vector adapted
    aux = tc(ni)*u(ni+1)
    u(ni) = u(ni) + aux
  End Do

! test whether T*u=b if u is solution that can become negative
! remember T= -TA TB -TC1; checked, OK

! l=1
! delta_b(l)=tb(l)*u(l)-tc1(l)*u(l+1)
! l=lmax
! delta_b(l)=-ta(l)*u(l-1)+tb(l)*u(l)

! do l=2,lz
! delta_b(l)=-ta(l)*u(l-1)+tb(l)*u(l)-tc1(l)*u(l+1)
! enddo

! rhs is orginal b with first value maybe modified (-> test=true)
! errm=maxval(abs(1.-delta_b(1:lmax)/rhs(1:lmax)))
! precision of 1.d-4 should be always reached (1.d-5 sometimes happens)
! if(errm.gt.1.d-4) then
! print*,'incons',k,jp,errm
! do l=1,lmax
! print*,l,delta_b(l),rhs(l)
! enddo
! endif

  delta_b(1:lmax) = 0.

! this is the new condition: u needs to be improved, such that it
! becomes positive
  If (minval(u(1:lmax))<0.) Then
    delta_u(1:lmax) = 0.
!   delta u should become slightly positive
    Where (u(1:lmax)<0)
      delta_u = -u*(1.D0+1.D-6) !   this should be > 0
    End Where

!   for tests, thus far always OK
!   do l=1,lmax
!   if(u(l).ge.0. .and. delta_u(l).ne.0.) stop ' error in delta_u'
!   enddo

!   calculate delta_b = T * delta_u, needed for moments equations
    l = 1
    If (maxval(delta_u(l:l+1))>0.) delta_b(l) = tb(l)*delta_u(l) - &
      tc1(l)*delta_u(l+1)
    l = lmax
    If (maxval(delta_u(l-1:l))>0.) delta_b(l) = -ta(l)*delta_u(l-1) + &
      tb(l)*delta_u(l)
    Do l = 2, lz
      If (maxval(delta_u(l-1:l+1))>0.) delta_b(l) = -ta(l)*delta_u(l-1) + &
        tb(l)*delta_u(l) - tc1(l)*delta_u(l+1)
    End Do

!   set u to improved valued: u(l)=u_previous(l) for u_previous(l)>0,
!   u(l)=1.d-6*abs(u_previous(l)) for u_previous(l) <0
!   use corresponding delta_b (integrated over dmue) in moments equations
    u = u + delta_u

  End If

110 Continue

! commented out, since confirmed
! do l=1,lmax
! if(tc(l).lt.0.) then
! print*,'tc_prime not positive'
! print*,xlambda,jp,l,tc(l)
! endif
! enddo


! uneghere=.false.
! final test
  Do l = 1, lmax
    If (u(l)<0.) Then
      Write (999, *) ' stop: u < 0'
      Stop ' u < 0'
!     JO March 2025
!     in the following, some ideas to cure u<0 if not explicitely set to zero
!     above. After many tests, currently there is no evidence that such new
!     ideas improve the situation. Thus, we stay with setting u to zero,
!     and comment out everything else

!     uneghere=.true.
!     if(l.eq.nd) stop 'u=0 at lower core boundary'
!     if(ub(l).gt.1.) then
!     ta(l)=0.
!     tc1(l)=0.
!     tb(l)=1.
!     rhs(l)=uorig(l)
!     else
!     if(l.ne.lmax) then
!     dxzm=0.5*(pp(l-1)+pp(l))*aaz1(l-1)
!     dxz=0.5*(pp(l)+pp(l+1))*aaz1(l)
!     ta(l)=ta(l)*(1.+dxzm)
!     tc1(l)=tc1(l)*(1.+dxz)
!     tb(l)=ta(l)+tc1(l)+1.
!     rhs(l)=s(l)
!     else
!     dxzm=0.5*(pp(l-1)+pp(l))*aaz1(l-1)
!     ta(l)=ta(l)*(1.+dxzm)
!     tb(l)=ta(l)+1.
!     rhs(l)=s(l)
!     endif
!     endif
    End If
  End Do

! JO March 2025 commented out
! if(uneghere) then
! tc=tc1
! u=rhs

! tb1 = 1.d0/tb(1)
! tc(1) = tc(1)*tb1
! u(1) = u(1)*tb1

! do l = 2,lz
! help = tb(l) - ta(l)*tc(l-1)
! h1 = 1.d0/help
! tc(l) = tc(l)*h1
! u(l) = (u(l)+u(l-1)*ta(l))*h1
! end do

! u(lmax) = (u(lmax)+u(lz)*ta(lmax))/ (tb(lmax)-tc(lz)*ta(lmax))

! do l = 1,lz
! ni = lmax - l
! u(ni) = u(ni) + tc(ni)*u(ni+1)

! if(u(ni).lt.0.) then
! print*,xlambda,jp,ni,u(ni)
! stop ' u still <0 after correction'
! endif
! end do
! endif

  counter_u = counter_u + 1
! JO March 2025 commented out, since u > 0 everywhere (otherwise, stop before)
! if(uneghere) counter_uneg=counter_uneg+1

! now u is the new field at index k

! inlinen von gmalu und
! v = h * v + gmalu

  b = u(1)
  Do l = 1, lz
    a = b
    b = u(l+1)
    v(l) = h(l)*v(l) + ga(l)*(a-b) ! Note ga(l) = -daz!
  End Do

! JO March 2025 commented out
! goto 1113

! until further evidence, this loop is not executed
! if(uneghere) then
! do l=1,lz
! help=0.5*(u(l)+u(l+1))
! if(abs(v(l)).gt.1.2*help) then
! counter_v=counter_v+1
! if(v(l).gt.0.) then
! v(l)=help
! else
! v(l)=-help
! endif
! endif
! enddo
! endif

! now v is the new field at index k

! adding the new u to aj

! JO March 2025 commented out
! 1113 continue

  Do l = 1, lmax
!   final test once more (just to be sure)
    If (u(l)<0.) Then
      Write (999, *) ' stop: u(l) < 0 in final loop (subr. ray_complete)'
      Stop ' u(l) < 0 in final loop (subr. ray_complete)'
    End If
!   no update required, since u needs to be nelected for Jnu, as well as uold
!   for alo
    pw = w0(l)
    aj(l) = aj(l) + u(l)*pw
    add_corr(l) = add_corr(l) + delta_b(l)*pw ! this is the term used in
!   moments equations which adds
!   to S
    ak(l) = ak(l) + u(l)*vp2(l, jp)
!   JO April 2018: previously, there was a cycle statement; however, such
!   contribution can be essential. Safetyfactor 0.9 anyway!
!   Since uold always < 1., this can happen only if SL > u;
!   in this case, the effective ALO also becomes < 0.9
    If (uold(l)*slinetotk(l)>u(l)) Then

!     indeed,  u = uold *slinetotk + psi (J = alo * slinetotk + psi (since
!     ff is included in alo (uold)), where psi accounts for background and
!     bluewing boundary, but many tests ensured that this test and reset
!     is  sufficient, without accounting for psi. To calculate psi, the above
!     procedure needs to be performed with slinetotk = 0, but including psi
!     in the condition above does not improve the reset.

      uo = 0.9*u(l)/slinetotk(l) !  do not modify uold itself
      If (uo>0.9) Then
        Write (999, *) ' stop: something really wrotten with uold'
        Stop ' something really wrotten with uold'
      End If
      alo(l) = alo(l) + uo*pw
    Else
      alo(l) = alo(l) + uold(l)*pw
    End If
  End Do

! Calculate H, needed for radiative acc, JS, and N, needed for moments eq.
! (JP)
! print*,xlambda
  Do l = 1, lmax - 1
    pw = w1(l+1)
    ah(l+1) = ah(l+1) + v(l)*pw
    an(l+1) = an(l+1) + v(l)*vp3(l+1, jp)
  End Do

! outer and inner boundary
  If (test) Then
!   outer boundary reset to Iminus = 0, thus v=u
    ah(1) = ah(1) + u(1)*w1(1)
    an(1) = an(1) + u(1)*vp3(1, jp)
  Else
    ah(1) = ah(1) + (u(1)-iminusold)*w1(1) ! iminusold corresponds here to
!   actual Iminus
    an(1) = an(1) + (u(1)-iminusold)*vp3(1, jp) ! iminusold corresponds here
!   to actual Iminus
  End If

  If (lmax==nd) Then
    ah(nd+1) = ah(nd+1) + (iplus-u(nd))*w1(nd+1)
!   integrals of u, required for moments equations
    hbnd = hbnd + u(nd)*w1(nd+1)
    nbnd = nbnd + u(nd)*vp3(nd+1, jp)
  End If

  Return

! obsolte
! return if not boundary condition for Iminus
! if (jp.ne.np-1) return

! ***  final value for p = rmax and z=0 in optically thick case!!

! print *,' optically thick outer boundary at',xlambda

! pw = w0last
! aj(1) = aj(1) + iminus*pw
! !ALO             if (s1.eq.sl) alo(1) = alo(1) + alo1*pw
! return

End Subroutine

!----------------------------------------------------------------------------

Subroutine moments_cmf(k, nd, aj, aj_ray, ah, ah_ray, jold, hold, r, gr, gt, &
  opakk, slinek, hplus, xlambda, counter_jneg, jneg)

! solves moments equations for given Eddington factors f and g.
! needs to be updated for g' in case of very small H in g=N/H.
! on output/input, J, H, J_old, and H_old have their usual meaning.
! during calculation, these values and eta are assumed to be multiplied by r^2

! note: whereas H has dimension ND+1 (only elements 2:ND on staggered grid
! will
! be used), gedd has dimension ND-1 (only staggered grid)

! as usual, we solve for (-ta, tb, -tc)J = x

! programmed by J.P. during a trip to Tenerife Feb. 2019
! (based on his diploma work, partly incorporting notes from W.-R. Hamann)

! includes new term add_corr_arr to correct source-function
! thus avoiding (in almost all cases) negative J_nu. Not tested here, only in
! moments_cmf1

  Use :: nlte_type
  Use :: nlte_dim, Only: id_ndept
  Use :: cmf_all_var, Only: fedd, gedd, hbound1, nbound1, hboundnd, nboundnd, &
    add_corr_arr

  Implicit None
  Integer (i4b), Parameter :: nd1 = id_ndept

  Integer (i4b), Intent (In) :: k, nd
  Integer (i4b), Intent (Inout) :: counter_jneg
  Logical, Intent (Out) :: jneg

  Real (dp), Intent (Out) :: aj(nd), ah(nd+1)
  Real (dp), Intent (In) :: jold(nd), hold(nd+1), r(nd), gr(nd), gt(nd), &
    opakk(nd), slinek(nd), aj_ray(nd), ah_ray(nd+1)
  Real (dp), Intent (In) :: hplus, xlambda

  Integer (i4b) :: l, lz, ni

  Real (dp), Dimension (nd1) :: ajold, rr, etak, q
  Real (dp), Dimension (nd1) :: ta, tb, tc, qfp1, qfm1, qfp0, qfm0, blp, blm
  Real (dp), Dimension (nd1-1) :: ahold, r2

  Real (dp) :: rl, rlm, rlp, opa0, opam, opap, gr0, grm, grp, gt0, gtm, gtp, &
    rrq, fl, flp, dr1, denp, denm, dtp, dtm, tb1, help, h1

! converting to quantities as used in the moments equations. Here,
! ahold = r^2* h_nue(k-1) has dimension ND-1, i.e., only staggered grid
! values.
  rr = r*r
  ajold = rr*jold
  etak = rr*opakk*(slinek+add_corr_arr(:,k))

  ahold = hold(2:nd)
  Do l = 1, nd - 1
    rl = 0.5D0*(r(l)+r(l+1))
    r2(l) = rl*rl
  End Do
  ahold = ahold*r2

! sphericality factor
  q(nd) = 1.D0
  rrq = 1.D0
  fl = 3.D0 - 1.D0/fedd(nd, k)
  Do l = nd - 1, 1, -1
    rl = r(l)
    rlp = r(l+1)
    flp = fl
    fl = 3.D0 - 1.D0/fedd(l, k)
    rrq = rrq*exp(fl-flp)*(rl/rlp)**((flp*rl-fl*rlp)/(rl-rlp))
    q(l) = rrq/rl/rl
  End Do


  Do l = 2, nd - 1
    rl = r(l)
    rlm = r(l-1)
    rlp = r(l+1)
    dr1 = 2.D0/(rlm-rlp)
    opa0 = opakk(l)
    opam = .5*(opakk(l-1)+opa0)
    opap = .5*(opakk(l+1)+opa0)
    gr0 = gr(l)
    grm = 0.5*(gr(l-1)+gr0)
    grp = 0.5*(gr(l+1)+gr0)
    gt0 = gt(l)
    gtm = 0.5*(gt(l-1)+gt0)
    gtp = 0.5*(gt(l+1)+gt0)

    denp = opap + (grp-gtp)*gedd(l, k) + gtp
    denm = opam + (grm-gtm)*gedd(l-1, k) + gtm

    dtp = 2.D0/(q(l)+q(l+1))/(rl-rlp)
    dtm = 2.D0/(q(l)+q(l-1))/(rlm-rl)
!   following quantities vectors, to enable re-calculation of H after J
!   in future versions, many of them can be made scalar
    qfp1(l) = dtp*q(l+1)*fedd(l+1, k)/denp
    qfm1(l) = dtm*q(l-1)*fedd(l-1, k)/denm
    qfp0(l) = dtp*q(l)*fedd(l, k)/denp
    qfm0(l) = dtm*q(l)*fedd(l, k)/denm

    blp(l) = ((grp-gtp)*gedd(l,k-1)+gtp)/denp
    blm(l) = ((grm-gtm)*gedd(l-1,k-1)+gtm)/denm

    ta(l) = qfm1(l)*dr1
    tc(l) = qfp1(l)*dr1
    tb(l) = (qfp0(l)+qfm0(l))*dr1 + (gr0-gt0)*fedd(l, k) + gt0 + opa0

    aj(l) = etak(l) + ((gr0-gt0)*fedd(l,k-1)+gt0)*ajold(l) - &
      dr1*(blm(l)*ahold(l-1)-blp(l)*ahold(l))
  End Do

! outer boundary
  l = 1
  dtp = 2.D0/(q(l)+q(l+1))/(r(l)-r(l+1))

  ta(l) = 0.D0
  tc(l) = dtp*q(l+1)*fedd(l+1, k)
  tb(l) = dtp*q(l)*fedd(l, k) + (gr(l)-gt(l))*nbound1(k) + &
    (gt(l)+opakk(l))*hbound1(k)

  aj(l) = ((gr(l)-gt(l))*nbound1(k-1)+gt(l)*hbound1(k-1))*ajold(l)

! inner boundary
  l = nd
  dtm = 2.D0/(q(l)+q(l-1))/(r(l-1)-r(l))
  ta(l) = dtm*q(l-1)*fedd(l-1, k)
  tc(l) = 0.
  tb(l) = dtm*q(l)*fedd(l, k) + (gr(l)-gt(l))*nboundnd(k) + &
    (gt(l)+opakk(l))*hboundnd(k)

  aj(l) = opakk(l)*hplus + ((gr(l)-gt(l))*nboundnd(k-1)+gt(l)*hboundnd(k-1))* &
    ajold(l)


! ***     inlining of invtri; remember: solution vector (aj) overwritten by
! solution


  lz = nd - 1

! same result as when calculated in qp (i.e., precision OK).

  tb1 = 1.D0/tb(1)
  tc(1) = tc(1)*tb1
  aj(1) = aj(1)*tb1

  Do l = 2, lz
    help = tb(l) - ta(l)*tc(l-1)
    h1 = 1.D0/help
    tc(l) = tc(l)*h1
    aj(l) = (aj(l)+aj(l-1)*ta(l))*h1
  End Do

  aj(nd) = (aj(nd)+aj(lz)*ta(nd))/(tb(nd)-tc(lz)*ta(nd))

  Do l = 1, lz
    ni = nd - l
    aj(ni) = aj(ni) + tc(ni)*aj(ni+1)
  End Do

  jneg = .False.
  Do ni = 1, nd
    If (aj(ni)<0.) Then
      If (.Not. jneg) counter_jneg = counter_jneg + 1
      jneg = .True.
    End If
  End Do

  If (jneg) Then
!   reset to ray-by-ray solution or to solution of previous iteration
!   (which might be ray-by-ray solution)
    aj = aj_ray
    ah = ah_ray
    Return
  End If

! now, aj is r^2*J_nu at index k
! let's calculate now r^2 H_nu (staggered grid)
! first defined matrix elements at l=2, corresponding to ah(3), but ahold(2)
! range: 2+1/2 to ND-1+1/2
  Do l = 2, lz
    ah(l+1) = qfp1(l)*aj(l+1) - qfp0(l)*aj(l) + blp(l)*ahold(l)
  End Do

! to check consistency: here, we calculate h from the backward quantities,
! should give identical results. Again, first defined matrix elements at l=2,
! but now corresponding to ah(2) and ahold(1)
! range: 1+1/2 to ND-2+1/2

! test can be commented out for production runs
  Do l = 3, lz !                    (last defined matrix elements at LZ
    help = qfm0(l)*aj(l) - qfm1(l)*aj(l-1) + blm(l)*ahold(l-1)
    If (abs(1.-help/ah(l))>1.D-16) Then
      Write (999, *) ' stop: subr. moments_cmf: problems in forward &
        &and backward matrix elements'
      Stop &
        ' subr. moments_cmf: problems in forward and backward matrix elements'
    End If
  End Do

! this is the point we need
  l = 2
  ah(2) = qfm0(l)*aj(l) - qfm1(l)*aj(l-1) + blm(l)*ahold(l-1)

! now we define H_nu at the outer and innermost points, using the same
! philosophy as above.
  ah(1) = hbound1(k)*aj(1)
  ah(nd+1) = hplus - hboundnd(k)*aj(nd)

! finally, correct for r^2 term
  aj = aj/rr
  ah(2:nd) = ah(2:nd)/r2
  ah(1) = ah(1)/r(1)**2
! nothing to be done at ND, since r(nd)=1.


  Return
End Subroutine

!----------------------------------------------------------------------------

Subroutine moments_cmf1(k, nd, aj, aj_ray, ah, ah_ray, jold, hold, r, gr, gt, &
  opakk, slinek, hplus, xlambda, counter_jneg, jneg)

! solves moments equations for given Eddington factors f and g' = N/J
! (staggered grid!)
! new method, following Hillier and Miller
! on output/input, J, H, J_old, and H_old have their usual meaning.
! during calculation, these values and eta are assumed to be multiplied by r^2

! note: whereas H has dimension ND+1 (only elements 2:ND on staggered grid
! will
! be used), geddp has dimension ND-1 (only staggered grid)

! as usual, we solve for (-ta, tb, -tc)J = x

! programmed by J.P. April 2019
! (based on his diploma work, partly incorporting notes from W.-R. Hamann)

! includes new term add_corr_arr to correct source-function
! thus avoiding (in almost all cases) negative J_nu.
  Use :: nlte_type
  Use :: nlte_dim, Only: id_ndept
  Use :: cmf_all_var, Only: fedd, geddp, hbound1, nbound1, hboundnd, nboundnd, &
    add_corr_arr

  Implicit None
  Integer (i4b), Parameter :: nd1 = id_ndept

  Integer (i4b), Intent (In) :: k, nd
  Integer (i4b), Intent (Inout) :: counter_jneg
  Logical, Intent (Out) :: jneg

  Real (dp), Intent (Out) :: aj(nd), ah(nd+1)
  Real (dp), Intent (In) :: jold(nd), hold(nd+1), r(nd), gr(nd), gt(nd), &
    opakk(nd), slinek(nd), aj_ray(nd), ah_ray(nd+1)
  Real (dp), Intent (In) :: hplus, xlambda

  Integer (i4b) :: l, lz, ni

  Real (dp), Dimension (nd1) :: ajold, rr, etak, q
  Real (dp), Dimension (nd1) :: ta, tb, tc, qfp1p, qfm1p, qfp0p, qfm0p, blpr, &
    blmr, blprr, blmrr
  Real (dp), Dimension (nd1-1) :: ahold, r2, ajold_stag

  Real (dp) :: rl, rlm, rlp, opa0, opam, opap, gr0, grm, grp, gt0, gtm, gtp, &
    rrq, fl, flp, dr1, denpr, denmr, dtp, dtm, tb1, help, h1

! converting to quantities as used in the moments equations. Here,
! ahold = r^2* h_nue(k-1) has dimension ND-1, i.e., only staggered grid
! values.
  rr = r*r
  ajold = rr*jold
  etak = rr*opakk*(slinek+add_corr_arr(:,k))

  ahold = hold(2:nd)
  Do l = 1, nd - 1
    rl = 0.5D0*(r(l)+r(l+1))
    r2(l) = rl*rl
    ajold_stag(l) = 0.5D0*(ajold(l)+ajold(l+1))
  End Do
  ahold = ahold*r2

! sphericality factor
  q(nd) = 1.D0
  rrq = 1.D0
  fl = 3.D0 - 1.D0/fedd(nd, k)
  Do l = nd - 1, 1, -1
    rl = r(l)
    rlp = r(l+1)
    flp = fl
    fl = 3.D0 - 1.D0/fedd(l, k)
    rrq = rrq*exp(fl-flp)*(rl/rlp)**((flp*rl-fl*rlp)/(rl-rlp))
    q(l) = rrq/rl/rl
  End Do


  Do l = 2, nd - 1
    rl = r(l)
    rlm = r(l-1)
    rlp = r(l+1)
    dr1 = 2.D0/(rlm-rlp)
    opa0 = opakk(l)
    opam = .5*(opakk(l-1)+opa0)
    opap = .5*(opakk(l+1)+opa0)
    gr0 = gr(l)
    grm = 0.5*(gr(l-1)+gr0)
    grp = 0.5*(gr(l+1)+gr0)
    gt0 = gt(l)
    gtm = 0.5*(gt(l-1)+gt0)
    gtp = 0.5*(gt(l+1)+gt0)

    denpr = opap + gtp
    denmr = opam + gtm

    dtp = 2.D0/(q(l)+q(l+1))/(rl-rlp)
    dtm = 2.D0/(q(l)+q(l-1))/(rlm-rl)
!   following quantities vectors, to enable re-calculation of H after J
!   in future versions, many of them can be made scalar
    qfp1p(l) = (dtp*q(l+1)*fedd(l+1,k)-(grp-gtp)*geddp(l,k)*0.5D0)/denpr
    qfm1p(l) = (dtm*q(l-1)*fedd(l-1,k)+(grm-gtm)*geddp(l-1,k)*0.5D0)/denmr
    qfp0p(l) = (dtp*q(l)*fedd(l,k)+(grp-gtp)*geddp(l,k)*0.5D0)/denpr
    qfm0p(l) = (dtm*q(l)*fedd(l,k)-(grm-gtm)*geddp(l-1,k)*0.5D0)/denmr

    blpr(l) = ((grp-gtp)*geddp(l,k-1))/denpr
    blprr(l) = gtp/denpr
    blmr(l) = ((grm-gtm)*geddp(l-1,k-1))/denmr
    blmrr(l) = gtm/denmr

    ta(l) = qfm1p(l)*dr1
    tc(l) = qfp1p(l)*dr1
    tb(l) = (qfp0p(l)+qfm0p(l))*dr1 + (gr0-gt0)*fedd(l, k) + gt0 + opa0

    aj(l) = etak(l) + ((gr0-gt0)*fedd(l,k-1)+gt0)*ajold(l) - &
      dr1*(blmrr(l)*ahold(l-1)-blprr(l)*ahold(l)) - dr1*(blmr(l)*ajold_stag(l- &
      1)-blpr(l)*ajold_stag(l))
  End Do

! outer boundary
  l = 1
  dtp = 2.D0/(q(l)+q(l+1))/(r(l)-r(l+1))

  ta(l) = 0.D0
  tc(l) = dtp*q(l+1)*fedd(l+1, k)
  tb(l) = dtp*q(l)*fedd(l, k) + (gr(l)-gt(l))*nbound1(k) + &
    (gt(l)+opakk(l))*hbound1(k)

  aj(l) = ((gr(l)-gt(l))*nbound1(k-1)+gt(l)*hbound1(k-1))*ajold(l)

! inner boundary
  l = nd
  dtm = 2.D0/(q(l)+q(l-1))/(r(l-1)-r(l))
  ta(l) = dtm*q(l-1)*fedd(l-1, k)
  tc(l) = 0.
  tb(l) = dtm*q(l)*fedd(l, k) + (gr(l)-gt(l))*nboundnd(k) + &
    (gt(l)+opakk(l))*hboundnd(k)

  aj(l) = opakk(l)*hplus + ((gr(l)-gt(l))*nboundnd(k-1)+gt(l)*hboundnd(k-1))* &
    ajold(l)


! ***     inlining of invtri; remember: solution vector (aj) overwritten by
! solution


  lz = nd - 1

  Go To 100
! check for diagonal dominance (can be commented out after confirmation)
  If (ta(1)/=0. .Or. tc(nd)/=0.) Stop ' ta or tc ne 0 in moments_cmf1'
  Do l = 1, nd
    If (ta(l)<0. .Or. tb(l)<0. .Or. tc(l)<0.) Then
      Print *, 'not positive in moments_cmf1'
      Print *, xlambda, l, ta(l), tb(l), tc(l)
    End If
    If (abs(ta(l))+abs(tc(l))>abs(tb(l))) Then
      Print *, 'not diagonally dominated in moments_cmf1'
      Print *, xlambda, l, ta(l), tb(l), tc(l)
    End If
  End Do

100 Continue

! same result as when calculated in qp (i.e., precision OK).

  tb1 = 1.D0/tb(1)
  tc(1) = tc(1)*tb1
  aj(1) = aj(1)*tb1

  Do l = 2, lz
    help = tb(l) - ta(l)*tc(l-1)
    h1 = 1.D0/help
    tc(l) = tc(l)*h1
    aj(l) = (aj(l)+aj(l-1)*ta(l))*h1
  End Do

  aj(nd) = (aj(nd)+aj(lz)*ta(nd))/(tb(nd)-tc(lz)*ta(nd))

  Do l = 1, lz
    ni = nd - l
    aj(ni) = aj(ni) + tc(ni)*aj(ni+1)
  End Do

  jneg = .False.
  Do ni = 1, nd
    If (aj(ni)<0.) Then
      If (.Not. jneg) Then
        counter_jneg = counter_jneg + 1
        Print *, 'new jneg at ', k, xlambda
      End If
      jneg = .True.
    End If
  End Do

  If (jneg) Then
!   reset to ray-by-ray solution or to solution of previous iteration
!   (which might be ray-by-ray solution)
    aj = aj_ray
    ah = ah_ray
!   do ni=1,nd
!   write(*,1) xlambda,k,ni,fedd(ni,k),geddp(ni,k)
!   enddo
    Return
  End If

! now, aj is r^2*J_nu at index k
! let's calculate now r^2 H_nu (staggered grid)
! first defined matrix elements at l=2, corresponding to ah(3), but ahold(2)
! range: 2+1/2 to ND-1+1/2
  Do l = 2, lz
    ah(l+1) = qfp1p(l)*aj(l+1) - qfp0p(l)*aj(l) + blpr(l)*ajold_stag(l) + &
      blprr(l)*ahold(l)
  End Do

! to check consistency: here, we calculate h from the backward quantities,
! should give identical results. Again, first defined matrix elements at l=2,
! but now corresponding to ah(2) and ahold(1)
! range: 1+1/2 to ND-2+1/2

! test can be commented out for production runs
  Do l = 3, lz !                    (last defined matrix elements at LZ
    help = qfm0p(l)*aj(l) - qfm1p(l)*aj(l-1) + blmr(l)*ajold_stag(l-1) + &
      blmrr(l)*ahold(l-1)
!   if(abs(1.-help/ah(l)).gt.1.d-16) then
!   allow for a bit larger differences, due to different compilers
    If (abs(1.-help/ah(l))>1.D-14) Then
      Write (999, *) help, ah(l)
      Write (999, *) ' stop: subr. moments_cmf1: problems in forward &
        &and backward matrix elements'
      Print *, help, ah(l)
      Stop ' subr. moments_cmf1: problems in forward and backward &
        &matrix elements'
    End If
  End Do

! this is the point we need
  l = 2
  ah(2) = qfm0p(l)*aj(l) - qfm1p(l)*aj(l-1) + blmr(l)*ajold_stag(l-1) + &
    blmrr(l)*ahold(l-1)

! now we define H_nu at the outer and innermost points, using the same
! philosophy as above.
  ah(1) = hbound1(k)*aj(1)
  ah(nd+1) = hplus - hboundnd(k)*aj(nd)

! finally, correct for r^2 term
  aj = aj/rr
  ah(2:nd) = ah(2:nd)/r2
  ah(1) = ah(1)/r(1)**2
! nothing to be done at ND, since r(nd)=1.

! for tests
! print*,xlambda,ah(1),ah_ray(1),ah(2),ah_ray(2)

  Return
110 Format (F12.5, 2X, I10, 2X, I4, 2(2X,E10.3))
End Subroutine

!----------------------------------------------------------------------------

Subroutine formal_obs(r, velo1, rmax)

! calculates (approx). emergent Eddington flux, normalized to Rstar,
! using u_outer

! JO: in case, needs to be updated for Iminus
! do we need u or 0.5*Iplus? (as calculated in subr. CMF_COMPLTE?
! In any case, u including Iminus needs to be corrected
! also for negative angles!

! JO first tests have shown that using u or u-0.5*Iminus gives similar results
! beyond the HeII edge, but different fluxes below (as to be expected).
! Interestingly, using u is more similar to cmfgen solutions.
! final decision will depend on comparison with solution from FORMALSOL.
! Thus far, we use u and 'live' with the inconsistency with fluxes from CONT
! below 228 A.

  Use :: nlte_type
  Use :: nlte_dim, Only: id_ndept, id_npoin
  Use :: fund_const, Only: clight
  Use :: fastwind_params, Only: optout_xj_xh

  Use :: nlte_var, Only: modnam, vmax, rtau23
  Use :: nlte_var, Only: p, vp1
  Use :: cmf_all_var, Only: nftot, vref, lam_cmf, u_outer, h_obs

  Implicit None
  Integer (i4b), Parameter :: nd = id_ndept, np = id_npoin

  Real (dp), Intent (In) :: velo1, rmax
  Real (dp), Dimension (nd), Intent (In) :: r

  Real (dp), Dimension (np) :: lamfac

  Integer (i4b) :: k, jp, ik, dk, dkmax, kint
  Real (dp) :: h, rmax2, vc, r1, lam1, u_inter

! usally allocated only one, since called only one (if UNASOL)
! but for tests, called more often. Thus
  If (.Not. allocated(h_obs)) Allocate (h_obs(nftot))
! JO Dec. 2018: to account for end-effects
  h_obs = 0.

  r1 = r(1)
  rmax2 = r1/rtau23

! check consistency
  If (abs(1.-rmax2/rmax)>1.D-12) Then
    Write (999, *) ' stop: problem with rtau23 in subr. formal_obs'
    Stop ' problem with rtau23 in subr. formal_obs'
  End If
! print*
! print*,' outermost radius point in units of nominal radius = ',rmax2
  rmax2 = rmax2*rmax2

  vc = velo1*vmax/clight

  dk = vmax/vref

! assuming now that lamcmf is a grid refering to lambda_obs, we calculate
! the corresponding cmf-wavelength where u needs to be evaluated
! according to lam_cmf = lam_obs/(1-mu*v/c) = lam_obs*lamfac(jp)
! Note that lam_obs (for p=Rmax) le lam_cmf le lam_obs/(1-v/c) (for p=0)
  Do jp = 1, np - 1
    lamfac(jp) = 1.D0/(1.D0-sqrt(1.D0-(p(jp)/r1)**2)*vc)
  End Do

  Do k = 1, nftot - (dk+1) !        to account for end effects
    h = 0.
    dkmax = dk
    Do jp = 1, np - 1
!     assuming v(pmax) = 0. because of symmetry
!     JO: modify when outer boundary included
      lam1 = lam_cmf(k)*lamfac(jp)
!     interpolate u(lam1=lam_cmf(lam_obs)) from u(lam)
!     jp=0:    lam1=lam_obs/(1-v/c)
!     jp=rmax: lam1=lam_obs
      Do ik = dkmax, 0, -1
        If (lam_cmf(k+ik+1)>lam1 .And. lam_cmf(k+ik)<lam1) Then
          Go To 100
        End If
      End Do
      Print *, k, jp, lam1, lam_cmf(k+ik+1), lam_cmf(k)
      Write (999, *) ' stop: index not found (subr. formal_obs)'
      Stop 'index not found (subr. formal_obs)'

100   kint = k + ik
!     linear interpolation
      u_inter = u_outer(jp, kint) + (u_outer(jp,kint+1)-u_outer(jp,kint))/( &
        lam_cmf(kint+1)-lam_cmf(kint))*(lam1-lam_cmf(kint))

!     integration
      h = h + u_inter*vp1(1, jp)
      dkmax = ik
    End Do
    h_obs(k) = h*rmax2
  End Do

  Print *
  Print *, ' H-obs calculated (and in case saved to file ''out_xh_obs'')'
  Print *, ' (H-obs corresponds to H_nu * (r(1)/Rstar)^2, contrasted to XXH)'

  If (optout_xj_xh) Then
    Open (1, File=trim(modnam)//'/out_xh_obs')
!   here lam_cmf is actually lam_obs (vacuum)
    Do k = 1, nftot
      Write (1, 110) lam_cmf(k), h_obs(k)
    End Do
    Close (1)
  End If

  Return

110 Format (F14.4, 2X, E10.4)
End Subroutine

!----------------------------------------------------------------------------

Subroutine test

! calculates Jbar for specified transition, and compares with SA
! (e.g., to check outer boundary condition)

  Use :: nlte_type
  Use :: nlte_dim, Only: id_ndept
  Use :: cmf_all_var, Only: icfirst, nftot, index_lam_cmf, lam_cmf, xlam1, &
    id1, indexel_inv, nf_cmf, xj_cmf, pweightdop

  Implicit None
  Integer (i4b), Parameter :: nd = id_ndept

  Integer (i4b) :: i, j, lastindex, in, index, kk, irest, ml, mu, k, ixmax, &
    jj, low, up, ik, l

  Real (dp), Dimension (nd) :: jbar

  lastindex = icfirst - 1
  Do i = 1, nftot
    in = index_lam_cmf(i)
    If (lam_cmf(i)>1540 .And. lam_cmf(i)<1560) Then
!     print*,i,lam_cmf(i),in
      If (in/=0) Then
        Do j = lastindex + 1, in
          index = id1(j)
!         print*,j,index
          If (index>=1D8) Then
!           explicit element
            kk = index/1D8
            irest = index - kk*1D8
            ml = irest/10000
            mu = irest - ml*10000
            k = indexel_inv(kk)
!           in units of vref = vturb/3 (with vturb=max(vturb,vturbmin)
            ixmax = nf_cmf(k)
!           print*,'exp',xlam1(j),kk,k,ml,mu,ixmax
          Else
            k = index/1000000
            irest = index - k*1000000
            jj = irest/100000
            irest = irest - jj*100000
            low = irest/100
            up = irest - low*100
            ixmax = nf_cmf(k)
            If (k==6 .And. low==1 .And. up==2) Then
              Print *, i, lam_cmf(i)
              Print *, ' bg', xlam1(j), k, low, up, ixmax
              jbar(:) = xj_cmf(:, i)*pweightdop(0, :, k) ! ik=0
              Do ik = 1, ixmax
                jbar(:) = jbar(:) + (xj_cmf(:,i+ik)+xj_cmf(:,i-ik))*pweightdop &
                  (ik, :, k)
              End Do
              Do l = 1, nd
                Print *, l, jbar(l)
              End Do
              Print *
            End If
          End If
        End Do
      End If
    End If
    If (in/=0) lastindex = in
  End Do

  Return
End Subroutine

!----------------------------------------------------------------------------

Subroutine calc_iminus(s1, xlambda, dx, taus, tauc, sl1, slinekblue, &
  iminusold, in_line_old, nstep, ncount)

! completely new outer boundary conditions, which are consistent with
! approximate treatment

! s1 is the output quantity, source term for outer boundary condition
! all quantities calculated for ray jp

! 1. Iminus - P/opa dIminus/dx (standard) or approx Iminus - (Iminus - S)
! (alternative)
! old outer boundary condition with Iminus alone resulted in increasing
! flux at long wavelengths. Origin due to neglect of dIminus/dx and
! setting u=0 for negative u
! 2. include opac(1) in total tau for Iminus; otherwise, optically thick
! continuum shortward of HeII edge incorrect

! to test outer boundary condition, compare with Iminus=0. (set s1=0.),
! and compare j(1) with j(2) and xj(1) from approx. treatment.

  Use :: nlte_type

  Implicit None

  Real (dp), Intent (Out) :: s1

  Real (dp), Intent (In) :: xlambda, dx, taus, tauc, sl1, slinekblue
  Real (dp), Intent (Inout) :: iminusold

  Integer (i4b), Intent (In) :: nstep
  Integer (i4b), Intent (Inout) :: ncount

  Logical, Intent (Inout) :: in_line_old

! local variables
  Real (dp) :: iminus, iminus1, sl, diminus

  If (taus>0.05 .And. taus>tauc) Then
!   dominating line processes, and line begins to become optically thick
    sl = sl1 !                      sl1 is more or less constant over line; in
!   this way, we avoid spikes
!   at the edges of the core (setting here slinek(1) is disastrous
    iminus = sl*(1.D0-exp(-taus-tauc)) ! including continuum tau

!   for closely separated lines, iminus from previous line (outer resonance
!   zones) might be larger than iminus from considered line core.
!   In this case, we use the former value for iminus.
!   Exact treatmeat would yield iminus_exact = sl(1-exp(-tau)) +
!   iminus1*exp(-tau)
    If (ncount/=0) Then
!     decaying iminus from previous line
      If (.Not. in_line_old) Then
        Write (999, *) &
          ' stop: calc_iminus: something rotten with in_line_old!'
        Stop ' calc_iminus: something rotten with in_line_old!'
      End If
      ncount = ncount + 1
      iminus1 = iminusold*(1.D0-1.D0/float(nstep-ncount+1))
      If (iminus1>iminus) Then
!       approximate iminus_exact=sl(1-exp(-tau)) + iminus1*exp(-tau)
!       by iminus1 (rather small tau)
        iminus = iminus1
        sl = 1.D99 !                hack to ensure that correct path chosen
!       later on
      Else
!       approximate iminus_exact=sl(1-exp(-tau)) + iminus1*exp(-tau)
!       by iminus=sl(1-exp(-tau)) (large tau)
        ncount = 0
!       end illumination from above resonance zones, since local term
!       dominates
      End If
    Else
      in_line_old = .True.
      ncount = 0 !                  counting begins when line core ends
    End If

  Else If (tauc>1.) Then
!   dominating continuum process, but including line
    sl = slinekblue !               total source function
    iminus = sl*(1.D0-exp(-taus-tauc))
    in_line_old = .True. !          also here, I-minus does not go down
!   instantaneously,
!   when tauc goes to zero, but account for the fact that
!   (optically thin) frequencies are illuminated from
!   bluewards intensities, which originate from higher levels
!   (Iminus(|mu|V,x) = Iminus((|mu|V+deltax,x+deltax)
    ncount = 0 !                    counting begins when continuum becomes
!   optically thin

!   here, we don't apply the upper correction for larger iminus1
!   in case, include it also here

  Else
!   optically thin continuum and line
    If (in_line_old) Then
!     end of optically thick line core (or optically thick continuum)
!     encountered; account for outer resonance zones (continuum):
!     approximated in such a way that Iminus decreases linearly to
!     zero, where the spatial extent is calculated from the condition
!     that irradiation from above must be possible
!     a) in obsframe: outwards resonance zones possible from
!     xobs = -(|mu|V)_last ... -1
!     b) in cmf: from -xmax ... -1 + (|mu|V)_last
!     approximated by 0 ... -1 + |mu|_last = -1 +(z/r)_last
!     THUS: additional NSTEPs = (1-z(1)/r(1))/DELTAX required,
!     where Iminus decreases monotonically to zero (corresponding
!     to the fact that the resonance zones move outwards and SLINE
!     and OPAL (or SCONT and OPAC) decrease
      If (nstep<2) Then
!       for core rays, the irradiation ends abruptly
        sl = 0.
        iminus = 0.
        in_line_old = .False.
        ncount = 0.
      Else

!       first step after end of line core (-xmax)
        If (ncount==0) Then
          ncount = 1
          iminus = iminusold*(1.-1./float(nstep))
        Else
!         next steps
          ncount = ncount + 1
          iminus = iminusold*(1.D0-1.D0/float(nstep-ncount+1))
        End If

        If (iminus==0.D0) Then
!         end of freq. range with outer res. zones encountered
          in_line_old = .False.
          ncount = 0
!         otherwise, leave in_line_old at true
        End If
        sl = 1.D99 !                hack to ensure that correct path chosen
!       later on
      End If

    Else
!     no irradiation from above
      sl = 0.
      iminus = 0.
!     consistency check
      If (ncount/=0) Then
        Write (999, *) &
          ' stop: subroutine calc_iminus: no irradiance, but ncount ne 0'
        Stop 'subroutine calc_iminus: no irradiance, but ncount ne 0'
      End If
!     this path implies that in_line_old is already .false.
    End If
  End If

! now we have calculated iminus, can proceed to calculate dIminus/dx and s1

  If (sl==slinekblue) Then
    diminus = iminus - sl
!   Alternative formulation for a dominating opt. thick continuum
!   diminus calculated according to
!   P/opa dIminus/dx approx Iminus - S -> effective source term is S.
!   Here, the standard formulation gives spurious results in lines (below 228
!   A)
!   or too much absorption at blue edge
  Else
    diminus = dx*(iminusold-iminus)
!   Standard formulation. Inside line, diminus <= 0. Becomes positive
!   when the line core 'ends'. Then, and most important, the effective source
!   term S1 becomes negative, and particularly the intensities for the
!   non-core rays are damped. Without this term (e.g., alternative
!   formulation from above) the intensities remain at the high values
!   from previous S1 terms, and the redward fluxes become to high!
!   This condition also applies when continuum becomes opt. thin
  End If

! This is the final value, Iminus - P/kappa dIminus/dx
  s1 = iminus - diminus

! for tests (include jp in header)
! if(jp.eq.1.and.xlambda.ge.1548.and.xlambda.lt.1555.) then
! write(*,fmt='(f10.3,9(2x,e8.2))') xlambda,s1,iminus,diminus, &
! &   sl,taus,tauc
! endif

! for next frequency
  iminusold = iminus

  Return
End Subroutine

!----------------------------------------------------------------------------

Subroutine calc_tlu(nd, ilow, imax, xne, temp, clf, r, v, dvdr)

! ilow,imax nlte values

! calculates reduced RBB rate_coefficients (tlumat, tulmat) for explicit
! elements
! indexcmf1 = 3,4,5,... : radiation field from cmf_complete
! (number of depacked components is indexcmf1-2
! indexcmf1 = 2 : radiation field from cmf_sing
! indexcmf1 = 1 : Sobolev approach

! TODO:
! attention for xlamc etc. (original from DATA or modified from
! xlamcmf_depacked)
! xlamc etc: attention if SA transport inside range

! IMPORTANT NOTE: all data regarding components as stored in met_rbb refer to
! original data (id,xlam,gf) and have been calculated in subr. fgrid_cmf

! if treat_single=T, corresponding line (defined in upper part of routine)
! will be treated in single line approximation with average background
! in particular, this refers to the HeI resonance line in case opt_hei_sing =
! T

  Use :: nlte_type
  Use :: fund_const, Only: c_tau => cross_class ! lambda in cgs
  Use :: fund_const, Only: e2mc2, hc2, hh, pi

  Use :: nlte_dim, Only: id_ndept, id_atoms, id_nttrd

  Use :: nlte_opt, Only: opt_hei_sing

  Use :: princesa_var, Only: index1, indat1, data, le, ifirsl, iong, gl
  Use :: nlte_var, Only: mmlow, mmup, indexrbb, indexcmf, vdopcmf, vther, &
    blucmf, bulcmf, aulcmf, sr, vmax, corrfc, lablinelo, lablineup, tulmat, &
    tlumat, indexsa, betold, no_alo
  Use :: cmf_multi, Only: indov
  Use :: nlte_app, Only: teff
  Use :: cmf_all_var, Only: indexcmf1, occexpl, indexel_inv, nf_cmf, &
    indexcmf1_lam, xj_cmf, alo_cmf, pweightdop_h, profdop_h, pweightdop_he, &
    profdop_he, pweightdop, profdop, sumopal_cmf, lineno, depacked, &
    indxlamc_depacked, xlamcmf_depacked, sumgf_depacked
  Use :: ng_var, Only: itng, index1_ng, sl_ng, ng
  Use :: nlte_porvor, Only: fic, tcl_fac_line, tcl_fac_cont
  Use :: tcorr_var, Only: temp_converged

  Implicit None
  Integer (i4b), Parameter :: nd1 = id_ndept, kel = id_atoms

! .. arguments ..
  Integer (i4b), Intent (In) :: nd
  Integer (i4b), Dimension (nd1, kel), Intent (In) :: ilow, imax

  Real (dp), Dimension (nd1), Intent (In) :: xne, temp, clf, r, v, dvdr

! .. local variables ..
  Integer (i4b), Dimension (:), Allocatable :: index_noc, iarr, indexx

  Real (dp), Dimension (nd1) :: ajlmean, alo, vr, xmust

  Real (dp) :: er, xlamc, flu, blu, bul, aul, xnl, xnu, xlam, vdop, snew, tlu, &
    tul, c_tau1, lam3, const, occlow, occup, opal, sline, jbar, alobar, rho1, &
    alo1bar, sumgf, lami, fi, consti, opali, jbari, alobari, gfi, lam3i, sl, &
    xlamcmf_dep, deltai, sumrho, const1, fip, lamp, ratio

  Integer (i4b) :: ii, indi, m, ml, mu, nato, indmat, indcmf, ll, ilo, ima, &
    numion, igenio, j, k, ke, il, ixmax, ik, i, noc, inoc, instart, inend, in, &
    iov, ij, il1, ikp

  Real (dp), Dimension (0:60) :: pw1, pw2

  Integer (i4b), Dimension (nd1, kel) :: nfir, nlas

  Logical :: inversion, treat_single, once

  Character (6) :: levlo, levup

! similar to routine cmfprep, but with fixed indices
! assuming that ionization equilibrium does not change any longer
! otherwise, more complex approach required

! calculate nfir and nlas once, to allow for subsequent checks

  If (maxval(nf_cmf)>60) Then
    Write (999, *) ' stop: dimension (pw1, pw2) > 60 in calc_tlu, adjust!'
    Stop ' dimension (pw1, pw2) > 60 in calc_tlu, adjust!'
  End If

  If (nd/=nd1) Then
    Write (999, *) ' stop: nd ne nd1 in calc_tlu'
    Stop ' nd ne nd1 in calc_tlu'
  End If

  Do k = 1, kel
    Do ll = 1, nd
      ilo = ilow(ll, k)
      ima = imax(ll, k)
      If (ima<ilo) Then
        Write (999, *) ' stop: ima < ilo in calc_tlu'
        Stop ' ima < ilo in calc_tlu'
      End If
      numion = igenio(k, ilo)
      nfir(ll, k) = ifirsl(numion)
      numion = igenio(k, ima) + 1
      If (numion>iong) Then
        Write (999, *) ' stop: error in calc_tlu - numion'
        Stop ' error in calc_tlu - numion'
      End If
      nlas(ll, k) = ifirsl(numion)
    End Do
  End Do

! prepara some aux variables
  Do ll = 1, nd
!   generalized blocking factors need to be updated
    Call factinc(r(ll), ll)
  End Do
  vr = v/r
  xmust = sqrt(1.D0-1.D0/r/r)

  c_tau1 = c_tau*sr/vmax


! preparation of cmf- and SA-transport
  Print *
iiloop: Do ii = 1, index1
    indi = indat1(ii)
    er = data(indi)
    m = int(er)
    If (m/=1) Cycle
    ml = mmlow(ii)
    mu = mmup(ii)

    nato = le(ml)

    xlamc = data(indi+1)*1.D-8
    flu = data(indi+2)

    indmat = indexrbb(ii)
    indcmf = indexcmf1(indmat) !    use new index

    levlo = lablinelo(indmat)
    levup = lablineup(indmat)
    xlamcmf_dep = xlamcmf_depacked(indmat)

    treat_single = .False.
!   if (levlo .eq. 'HE21  '.and. levup .eq. 'HE22  ') treat_single=.true.
    If (opt_hei_sing) Then
      If (indcmf==3 .And. (levlo=='HE11S1  ' .And. xlamcmf_dep>584. .And. &
        xlamcmf_dep<585.)) treat_single = .True.
    End If

!   define index-array for depacked lines
    sumgf = sumgf_depacked(indmat)
    If (sumgf/=0) Then
      Do i = 1, lineno
        If (depacked(i)%levlo==levlo .And. depacked(i)%levup==levup) Then
          noc = depacked(i)%nocomp
          Allocate (index_noc(noc), iarr(noc), indexx(noc))
          inoc = 1
          index_noc(inoc) = i
          instart = i
          Exit
        End If
      End Do
      If (noc/=1) Then
        Do i = instart + 1, lineno
          If (depacked(i)%levlo==levlo .And. depacked(i)%levup==levup) Then
            inoc = inoc + 1
            index_noc(inoc) = i
          End If
        End Do
      End If
!     test
      If (inoc/=noc) Then
        Write (999, *) ' stop: something wrong with inoc (subr. calc_tlu)'
        Stop ' something wrong with inoc (subr. calc_tlu)'
      End If
!     sort index_noc according to index (wavelength)
!     most components in LINES_xxx.dat already sorted, but not necessarily
      Do i = 1, noc
        iarr(i) = depacked(index_noc(i))%index
      End Do
      Call indexx_int(noc, iarr, indexx)
!     aux. variable
      iarr = index_noc
      Do i = 1, noc
        index_noc(i) = iarr(indexx(i))
      End Do

    End If !                        sumgf ne 0

!   ---     final check if everything was ok

    once = .True.
    Do ll = 1, nd

      If (ml<nfir(ll,nato) .Or. ml>=nlas(ll,nato)) Then
        If (tlumat(indmat,ll)/=0.) Then
          Write (999, *) ' stop: tlumat ne 0(1)'
          Stop ' tlumat ne 0(1)'
        End If
        Cycle
      End If
      If (mu>nlas(ll,nato)) Then
        If (tlumat(indmat,ll)/=0.) Then
          Write (999, *) ' stop: tlumat ne 0(2)'
          Stop ' tlumat ne 0(2)'
        End If
        Cycle
      End If

      xnl = occexpl(ml, ll)
      xnu = occexpl(mu, ll)
!     all occupation numbers considered here should be defined
      If (xnl*xnu==0.) Then
        Write (999, *) ' stop: occ. num not defined in calc_tlu'
        Stop ' occ. num not defined in calc_tlu'
      End If

      If (xlamcmf_dep==0.D0) Then
        Write (999, *) ' stop: error in index gymnastics (calc_tlu)'
        Stop ' error in index gymnastics (calc_tlu)'
      End If

      If (indcmf==1) Then
!       consistency check
        If (indexcmf(indmat)/=1) Then
          Write (999, *) ' stop: indexcmf1=1 and indexcmf ne 1 in calc_tlu'
          Stop ' indexcmf1=1 and indexcmf ne 1 in calc_tlu'
        End If
!       if (indexsa(indmat).eq.0) stop ' indexsa=0 for SA-lines in calc_tlu'
!       in present philosophy (cmf_all implies GLOBAL=.TRUE.),
!       indexsa should have a value of ND
!       JO Nov. 2019 this might happen for cooler Temp. if He different
!       ilow/imax as a function of depth
!       Meanwhile, hack in routine ILOWIMAX_LTE in nlte.f90
        If (indexsa(indmat)/=nd) Then
          Write (999, *) ' stop: indexsa ne nd for SA-lines in calc_tlu'
          Stop ' indexsa ne nd for SA-lines in calc_tlu'
        End If

!       ---  Sobolev treatment, update of tulmat and tlumat
!       ---  NOTE: Sobolev lines (even inside range) should be never depacked.
!       Thus working with packed data
!       one more consistency check
        If (abs(1.D0-xlamcmf_dep*1.D-8/xlamc)>1.D-10) Then
          Write (999, *) xlamcmf_dep, xlamc*1.D8
          Write (999, *) &
            ' stop: problem with xlamc in SA-path (subr. calc_tlu)'
          Print *, xlamcmf_dep, xlamc*1.D8
          Stop ' problem with xlamc in SA-path (subr. calc_tlu)'
        End If
        betold(indmat, ll) = 0.D0
        If (ll==1) vdopcmf(indmat) = 0.D0
        Call tratrad(ll, ml, mu, flu, xlamc, xnl, xnu, xne(ll), gl, vr(ll), &
          dvdr(ll), tlu, tul, betold(indmat,ll), sr, vmax, xmust(ll), ii, &
          temp(ll), clf(ll), vther(nato), fic(ll), tcl_fac_line(ll), &
          tcl_fac_cont(ll))
!       JO SA path needs to be checked when SA actually calculated
!       (inside range where tlu and tul not set to zero)
!       Roughly checked, but note that I_inc etc. as used within range are now
!       cmf-values with a lot of structure in frequency.
!       Might be not what you want (which are smoothed values)

!       if(ll.eq.1 .or. ll.eq.10 .or. ll.eq.30 .or. ll.eq. 40 .or. ll.eq.50) &
!       &
!       print*,indmat,ll,tlu,tlumat(indmat,ll),tul,tulmat(indmat,ll)

        tlumat(indmat, ll) = tlu
        tulmat(indmat, ll) = tul

        If (once) Write (*, Fmt=130) xlamcmf_dep, levlo, levup, &
          indexsa(indmat)
        once = .False.
!       else if (indcmf.eq.2) then
!       if HeI resonance line should be treated in single line approx.
      Else If (indcmf==2 .Or. treat_single) Then

!       prepare single-line cmf-transport for lines outside range
!       remember: tlumat and tulmat will be overwritten by opal and sline
!       NOTE: outside range, all lines should be packed
!       consistency check
        If (indexcmf(indmat)/=2) Then
          Write (999, *) ' stop: indexcmf1=2 and indexcmf ne 2 in calc_tlu'
          Stop ' indexcmf1=2 and indexcmf ne 2 in calc_tlu'
        End If
        If (abs(1.D0-xlamcmf_dep*1.D-8/xlamc)>1.D-10) Then
          Write (999, *) xlamcmf_dep, xlamc*1.D8
          Write (999, *) &
            ' stop: problem with xlamc in CMFSING-path (subr. calc_tlu)'
          Print *, xlamcmf_dep, xlamc*1.D8
          Stop ' problem with xlamc in CMFSING-path (subr. calc_tlu)'
        End If
!       if(ll.eq.1 .or. ll.eq.10 .or. ll.eq.30 .or. ll.eq. 40 .or. ll.eq.50) &
!       &          print*,indmat,ll,tlumat(indmat,ll),tulmat(indmat,ll)
!       print*,xlamc,indcmf,treat_single
        Call tratcmf(ll, ml, mu, flu, xlamc, xnl, xnu, xne(ll), clf(ll), gl, &
          sr, vmax, ii, fic(ll), tcl_fac_line(ll), tcl_fac_cont(ll))

!       consistency check
!       JO can be skipped when everything runs
        If (ll==1) Then
          If (vdopcmf(indmat)/=vther(nato)) Then
            Write (999, *) ' stop: vdopcmf inconsistent in calc_tlu'
            Stop ' vdopcmf inconsistent in calc_tlu'
          End If
          blu = e2mc2/hh*xlamc*4.D0*pi**2*flu
          bul = gl(ml)*blu/gl(mu)
          aul = hc2/xlamc**3*bul
          If (blucmf(indmat)/=blu) Then
            Write (999, *) ' stop: blu inconsistent in calc_tlu'
            Stop ' blu inconsistent in calc_tlu'
          End If
          If (bulcmf(indmat)/=bul) Then
            Write (999, *) ' stop: bul inconsistent in calc_tlu'
            Stop ' bul inconsistent in calc_tlu'
          End If
          If (aulcmf(indmat)/=aul) Then
            Write (999, *) ' stop: aul inconsistent in calc_tlu'
            Stop ' aul inconsistent in calc_tlu'
          End If
        End If

!       else if (indcmf.eq.3) then
      Else If (indcmf>=3 .And. .Not. treat_single) Then

!       (re-)calculate line opacity and source-function
!       This needs to be done as in subr. opacity_cmf, to ensure consistency
!       of the ALO

!       consistency check
        If (indexcmf(indmat)/=2) Then
          Write (999, *) ' stop: indexcmf1 ge 3 and indexcmf ne 2 in calc_tlu'
          Stop ' indexcmf1 ge 3 and indexcmf ne 2 in calc_tlu'
        End If
        lam3 = xlamc**3 !           in cgs
        const = c_tau1*xlamc*flu*gl(ml)
        occlow = xnl/gl(ml)
        occup = xnu/gl(mu)
!       pi e2/me c * SR/vmax * lam * gf * (nl/gl - nu/gu) /clfac
        opal = const*(occlow-occup)/clf(ll)
        sline = hc2/lam3/(occlow/occup-1.D0)

!       consistency check
!       JO can be skipped when everything runs
        If (ll==1) Then
          If (vdopcmf(indmat)/=vther(nato)) Then
            Write (999, *) ' stop: vdopcmf inconsistent in calc_tlu'
            Stop ' vdopcmf inconsistent in calc_tlu'
          End If
          blu = e2mc2/hh*xlamc*4.D0*pi**2*flu
          bul = gl(ml)*blu/gl(mu)
          aul = hc2/xlamc**3*bul
          If (blucmf(indmat)/=blu) Then
            Write (999, *) ' stop: blu inconsistent in calc_tlu'
            Stop ' blu inconsistent in calc_tlu'
          End If
          If (bulcmf(indmat)/=bul) Then
            Write (999, *) ' stop: bul inconsistent in calc_tlu'
            Stop ' bul inconsistent in calc_tlu'
          End If
          If (aulcmf(indmat)/=aul) Then
            Write (999, *) ' stop: aul inconsistent in calc_tlu'
            Stop ' aul inconsistent in calc_tlu'
          End If
        End If
!       calculate reduced rates
!       ALO for single line is int alo_nu rho1_nu phi dnu, with rho1_nu = opal
!       * phi/opal_tot_nu
!       note that alo_nu already includes rho_nu = chi_nu(line_tot)/chi_nu_tot

        ke = indexel_inv(nato)
        ixmax = nf_cmf(ke)
        il = indexcmf1_lam(indmat) ! central frequency of strongest component
!       test, can be skipped later on
        If (ixmax==0) Then
          Write (999, *) ll, nato, ke, ixmax, il
          Write (999, *) ' stop: ixmax = 0 in calc_tlu'
          Print *, ll, nato, ke, ixmax, il
          Stop 'ixmax = 0 in calc_tlu'
        End If

!       integration over phi
        jbar = 0.
        alobar = 0.

        If (sumgf==0.) Then
!         path for packed components with original frequency or
!         for lineno=0 (when only HHe or no depacked lines present)

!         JO Dec. 2016  this should also apply for lines with no components in
!         long line list when treated as explicit element
!         (corresponding to quality=2 when treated as bg element, see
!         calc_tlu_met)

!         here we have to sum up jbar and alobar (one comp), and
!         to set alo1bar = alobar at the end
!         if(ll.eq.1) write (*,fmt=9005) xlamc*1.d8,levlo,levup
          If (ll==1) Write (*, Fmt=110) xlamcmf_dep, levlo, levup

          If (opal<0. .Or. no_alo) Then ! don't use ALO
!           alobar remains zero
            If (ke==1) Then
              pw1(0:ixmax) = pweightdop_h(0:ixmax, ll)
            Else If (ke==2) Then
              pw1(0:ixmax) = pweightdop_he(0:ixmax, ll)
            Else
              pw1(0:ixmax) = pweightdop(0:ixmax, ll, ke)
            End If

!           red side including center
            Do ik = 0, ixmax
              jbar = jbar + xj_cmf(ll, il+ik)*pw1(ik)
            End Do
!           blue side
            Do ik = 1, ixmax
              jbar = jbar + xj_cmf(ll, il-ik)*pw1(ik)
            End Do
!           standard path, no inversion, use ALO

          Else

            If (ke==1) Then
              pw1(0:ixmax) = pweightdop_h(0:ixmax, ll)
              pw2(0:ixmax) = profdop_h(0:ixmax, ll)
            Else If (ke==2) Then
              pw1(0:ixmax) = pweightdop_he(0:ixmax, ll)
              pw2(0:ixmax) = profdop_he(0:ixmax, ll)
            Else
              pw1(0:ixmax) = pweightdop(0:ixmax, ll, ke)
              pw2(0:ixmax) = profdop(0:ixmax, ll, ke)
            End If

!           red side including center
            Do ik = 0, ixmax
              jbar = jbar + xj_cmf(ll, il+ik)*pw1(ik)
              rho1 = opal*pw2(ik)/sumopal_cmf(ll, il+ik)
!             significant or dominating inversion of bg opacities
              If (rho1>1.1D0 .Or. rho1<0.D0) rho1 = 0.D0
!             very slight difference in lambda (SP vs. DP) or mild inversion
!             of bg
              If (rho1>1.D0) rho1 = 1.D0
              alobar = alobar + alo_cmf(ll, il+ik)*rho1*pw1(ik)
            End Do

!           blue side
            Do ik = 1, ixmax
              jbar = jbar + xj_cmf(ll, il-ik)*pw1(ik)
              rho1 = opal*pw2(ik)/sumopal_cmf(ll, il-ik)
!             significant or dominating inversion of bg opacities
              If (rho1>1.1D0 .Or. rho1<0.D0) rho1 = 0.D0
!             very slight difference in lambda (SP vs. DP) or mild inversion
!             of bg
              If (rho1>1.D0) rho1 = 1.D0
              alobar = alobar + alo_cmf(ll, il-ik)*rho1*pw1(ik)
            End Do
!           if(levlo.eq.'N31' .and. levup.eq.'N311' .and. ll.eq.44) then
!           do ik=-100,100
!           print*,lam_cmf(il+ik),sumopal_cmf(ll,il+ik)
!           enddo
!           endif


          End If !                  inversion or not
          alo1bar = alobar !        in this case (packed transitions),
!         no difference between alobar and alo1bar
!         if no_alo, alobar=0, and thus alo1bar=0

        Else
!         path for depacked components or single component with modified
!         frequency
          If (ke==1 .Or. ke==2) Then
            Write (999, *) ' stop: H/He line in calc_tlu (depacked)'
            Stop ' H/He line in calc_tlu (depacked)'
          End If

!         here we have to sum up jbari and alobari for individual components,
!         and then to calculate jbar, alobar and alo1bar for the packed
!         transition
!         within the rate equations by weighting with gfi
          alo1bar = 0.

          Do in = 1, noc

            i = index_noc(in)
!           print*,i,xlamcmf_dep,depacked(i)%wave
            If (depacked(i)%levlo/=levlo .Or. depacked(i)%levup/=levup) Then
              Write (999, *) ' stop: something wrong with index_noc'
              Stop ' something wrong with index_noc'
            End If
            lami = depacked(i)%wave
            fi = depacked(i)%flu

            If (ll==1 .And. lami==xlamcmf_dep) Write (*, Fmt=120) xlamcmf_dep, &
              levlo, levup

            lami = lami*1.D-8
            consti = lami*fi/(xlamc*flu)
            opali = opal*consti
            il = depacked(i)%index
!           only one component
            If (noc==1) Then
              If (abs(il-indexcmf1_lam(indmat))>1) Then
!               under specific circumstances, there might be an
!               index-difference of 1
!               must not be more!
                Write (999, *) ll, il, indexcmf1_lam(indmat)
                Write (999, *) ' stop: error in index for lamcmf (ncomp=1)'
                Print *, ll, il, indexcmf1_lam(indmat)
                Stop ' error in index for lamcmf (ncomp=1)'
              End If
              instart = in
              inend = in
            Else
!             overlapping components within profile; remember: index_noc
!             wavelength-ordered
              Do iov = in, 1, -1
                ij = index_noc(iov)
                deltai = depacked(ij)%index - il
!               for tests
                If (deltai>0) Then
                  Stop 'deltai > 0 (calc_tlu)'
                  Stop 'deltai > 0 (calc_tlu)'
                End If
                If (abs(deltai)<=ixmax) instart = iov
              End Do
              Do iov = in, noc
                ij = index_noc(iov)
                deltai = depacked(ij)%index - il
!               for tests:
                If (deltai<0) Then
                  Write (999, *) ' stop: deltai < 0 (calc_tlu)'
                  Stop ' deltai < 0 (calc_tlu)'
                End If
                If (deltai<=ixmax) inend = iov
              End Do
!             if(ll.eq.1) print*,levlo,levup,depacked(i)%lineno,instart,inend
            End If
            jbari = 0.
            alobari = 0.

            If (opali<0. .Or. no_alo) Then ! don't use ALO
!             alobari remains zero
              pw1(0:ixmax) = pweightdop(0:ixmax, ll, ke)

!             red side including center
              Do ik = 0, ixmax
                jbari = jbari + xj_cmf(ll, il+ik)*pw1(ik)
              End Do
!             blue side
              Do ik = 1, ixmax
                jbari = jbari + xj_cmf(ll, il-ik)*pw1(ik)
              End Do
!             standard path, no inversion, use ALO

            Else

              pw1(0:ixmax) = pweightdop(0:ixmax, ll, ke)
              pw2(0:ixmax) = profdop(0:ixmax, ll, ke)

!             red side including center
              Do ik = 0, ixmax
                jbari = jbari + xj_cmf(ll, il+ik)*pw1(ik)
                If (instart==inend) Then
                  rho1 = opali*pw2(ik)/sumopal_cmf(ll, il+ik)
                Else
                  sumrho = opali*pw2(ik)
                  const1 = 0.
                  Do iov = instart, inend
                    If (iov==in) Cycle
                    ij = index_noc(iov)
                    il1 = depacked(ij)%index
                    ikp = abs(ik+il-il1) ! other component's profile index
                    If (ikp>ixmax) Cycle ! outside other's component range
                    fip = depacked(ij)%flu
                    lamp = depacked(ij)%wave*1.D-8 ! in cgs
!                   correction for (slightly) different frequencies in SL:
!                   as calculated ALO refers to SL(iov);
!                   however, we want to apply it to SL(in). Thus:
!                   const1=const1+(fip/fi)*(lamp/lami)*pw2(ikp)*(lami/lamp)**3
!                   shorter
                    const1 = const1 + (fip/fi)*pw2(ikp)*(lami/lamp)**2
                  End Do
                  sumrho = sumrho + const1*opali
                  rho1 = sumrho/sumopal_cmf(ll, il+ik)
                End If

!               significant or dominating inversion of bg opacities
                If (rho1>1.1D0 .Or. rho1<0.D0) rho1 = 0.D0
!               very slight difference in lambda (SP vs. DP) or mild inversion
!               of bg
                If (rho1>1.D0) rho1 = 1.D0
                alobari = alobari + alo_cmf(ll, il+ik)*rho1*pw1(ik)
              End Do

!             blue side
              Do ik = 1, ixmax
                jbari = jbari + xj_cmf(ll, il-ik)*pw1(ik)
                If (instart==inend) Then
                  rho1 = opali*pw2(ik)/sumopal_cmf(ll, il-ik)
                Else
                  sumrho = opali*pw2(ik)
                  const1 = 0.
                  Do iov = instart, inend
                    If (iov==in) Cycle
                    ij = index_noc(iov)
                    il1 = depacked(ij)%index
                    ikp = abs(-ik+il-il1) ! other component's profile index
                    If (ikp>ixmax) Cycle ! outside other's component range
                    fip = depacked(ij)%flu
                    lamp = depacked(ij)%wave*1.D-8 ! in cgs
!                   correction for (slightly) different frequencies in SL:
!                   as calculated ALO refers to SL(iov);
!                   however, we want to apply it to SL(in). Thus:
!                   const1=const1+(fip/fi)*(lamp/lami)*pw2(ikp)*(lami/lamp)**3
!                   shorter
                    const1 = const1 + (fip/fi)*pw2(ikp)*(lami/lamp)**2
                  End Do
                  sumrho = sumrho + const1*opali
                  rho1 = sumrho/sumopal_cmf(ll, il-ik)
                End If

!               significant or dominating inversion of bg opacities
                If (rho1>1.1D0 .Or. rho1<0.D0) rho1 = 0.D0
!               very slight difference in lambda (SP vs. DP) or mild inversion
!               of bg
                If (rho1>1.D0) rho1 = 1.D0
                alobari = alobari + alo_cmf(ll, il-ik)*rho1*pw1(ik)
              End Do

            End If !                inversion or not

!           sum up the contributions from individual components
            gfi = fi*gl(ml)
            lam3i = lam3/lami**3 !  both quantities in cgs
            jbar = jbar + jbari*gfi
            alobar = alobar + alobari*gfi
            alo1bar = alo1bar + alobari*gfi*lam3i ! includes correction for
!           sline_i
!           print*,ll,depacked(i)%wave,depacked(i)%levlo,depacked(i)%levup
!           print*,jbari,alobari,gfi

          End Do !                  loop over in (number of components)

!         if(levlo.eq.'N31' .and. levup.eq.'N311' .and. ll.eq.44) then
!         do ik=-100,100
!         print*,lam_cmf(il+ik),sumopal_cmf(ll,il+ik)
!         enddo
!         endif

!         now, calculate actual jbar and alo by normalizing with sumgf and
!         sumgflam3 = sum (gfi * lam^3/lam_i^3), respecively
!         JO Dec 2016: sumgflam3 does not exist, and was not used in previous
!         version
!         of calc_tlu_met (approach changed from early on, and forgot to
!         delete comment?)
          jbar = jbar/sumgf
          alobar = alobar/sumgf
          alo1bar = alo1bar/sumgf
!         print*,jbar,alobar,sumgf
!         print*

        End If !                    packed or not
!       here, all required info should be present;
!       for only one transition (packed), alo1bar = alobar
        If (itng==2) Then
!         actually, jbar corresponds to J(n-1) corresponding itng=3,
!         but itng is only updated (in subr. prep_ng) after present routine
!         has been called. Thus, itng=2 still.
          ik = index1_ng(indmat)
          If (ik/=0) sl_ng(ll, ik, 5) = jbar
        End If

!       in case, use Ng-extrapolated source-function.
!       no rescaling of lambda required here (calculated with same xlamc)
        If (ng) Then
          ik = index1_ng(indmat)
          If (ik/=0) Then
            sl = sl_ng(ll, ik, 4)
            If (sl>0.D0) sline = sl
          End If
        End If

        If (abs(alobar)>1.D0) Then
!         print*,' warning!!!! something wrong with line-alo in calc_tlu'
          Write (999, *) &
            ' stop: warning!!!! something wrong with line-alo in calc_tlu'
          Stop ' warning!!!! something wrong with line-alo in calc_tlu'
          alobar = 0.D0
        End If

!       here we have to consider the differences between packed und individual
!       source function, thus working with alo1_tot
!       the alo itself, however, is independent of source function and thus
!       does
!       not need to be corrected
        snew = jbar - alo1bar*sline

!       for tests (no ALO)
!       snew = jbar
!       alobar = 0.

        If (levlo=='N31' .And. levup=='N311') Then
          Write (18, *) 'N', ll, occlow, occup
!         print*,xlamc,ll,indmat
!         write(18,*) jbar,alobar,alo1bar,sline,snew
        End If

!       if(levlo .eq. 'HE11S1  '.and. levup .eq. 'HE12P1  ') &
!       print*,ll,jbar,sline,snew
!       if(levlo .eq. 'HE21    '.and. levup .eq. 'HE22    ') &
!       write(*,222) ll,jbar,alobar,sline,snew
!       222 format('calcHE2',i3,4(2x,e12.6))


        If (snew<0.D0) Then
!         JO March 2025
          Write (999, *) ' stop: problems (1) with snew in calc_tlu. &
            &Should never happen!!!'
          Stop ' problems (1) with snew in calc_tlu. Should never happen!!!'

!         JO Dec. 2016: in case, modify according to calc_tlu_met
          Print *, ' warning!!!! something wrong with snew in calc_tlu'
          Print *, xlamc, ll, indmat
          Print *, jbar, alobar, sline, snew
          If (abs(snew)<sline/1000.) Then
!           snew=0.
!           print*,'numerical problems, snew set to zero!'
          Else
!           stop ' warning!!!! something wrong with snew in calc_tlu'
            Write (999, *) &
              ' warning!!!! something wrong with snew in calc_tlu'
            Print *, ' warning!!!! something wrong with snew in calc_tlu'
            snew = jbar
            alobar = 0.
          End If
        End If

!       hack to avoid oscillations in strongly coupled lines
        If (alobar<0.1) Then
          snew = jbar
          alobar = 0.D0
        End If


        If (temp_converged) Then
!         this is important for OB-stars: HeII 303 needs to converge as
!         reliable as possible,
!         thus do not prohibit ALOs close to unity
          If (alobar>0.99 .And. .Not. (levlo=='HE21    ' .And. levup== &
            'HE22    ')) Then
            ratio = alo1bar/alobar
            alobar = 0.99D0 !       underestimate allowed
            alo1bar = ratio*0.99D0
            snew = jbar - alo1bar*sline
          End If
        End If

!       if(levup.eq.'N54   ' .and.ll.le.15) &
!       & write(*,fmt='(A6,2x,i2,2x,e10.4,2x,f10.4,2x,2(e10.4,2x))') &
!       & levlo,ll,jbar,alobar,sline,snew

!       JO check whether this is also OK for inversion and alo=0
        tlu = blucmf(indmat)*snew
        tul = bulcmf(indmat)*snew + aulcmf(indmat)*(1.D0-alobar)
!       if(xlamc.gt.955.33d-8 .and. xlamc.lt.955.4d-8) then
!       print*,xlamc,ll,jbar,alobar,sline,snew,tlu,tul
!       endif


!       for tests
!       if((ll .le. 40) .and. levlo .eq. 'N51   ' .and. levup .eq. 'N52   ')
!       THEN
!       !           if(ll.eq.1 .or. ll.eq.10 .or. ll.eq.30 .or. ll.eq. 40 .or.
!       ll.eq.50) &
!       write(*,fmt='(f10.3,3(2x,i4),2(2x,e8.2))') xlamc*1.d8,ii,indmat,ll, &
!       &            1.-tlu/tlumat(indmat,ll),1.-tul/tulmat(indmat,ll)
!       print*,tlu,tul,jbar,alobar,sline,snew
!       endif
!       if((ll .le. 40) .and. levlo .eq. 'N44   ' .and. levup .eq. 'N429  ')
!       THEN
!       if(ll.eq.1 .or. ll.eq.10 .or. ll.eq.30 .or. ll.eq. 40 .or. ll.eq.50) &
!       write(*,fmt='(f10.3,3(2x,i4),2(2x,e8.2))') xlamc*1.d8,ii,indmat,ll, &
!       &            1.-tlu/tlumat(indmat,ll),1.-tul/tulmat(indmat,ll)
!       print*,tlu,tul,jbar,alobar,sline,snew
!       endif

        tlumat(indmat, ll) = tlu
        tulmat(indmat, ll) = tul

      Else
        Write (999, *) ' stop: error in indexcmf1(calc_tlu)'
        Stop ' error in indexcmf1(calc_tlu)'

      End If

    End Do !                        depthloop

    If (sumgf/=0.) Deallocate (index_noc, iarr, indexx)

  End Do iiloop
  Print *

! perform cmf-transport for lines outside range, update tulmat and tlumat
  Do j = 1, id_nttrd

    ii = indxlamc_depacked(j)
    xlam = xlamcmf_depacked(ii)

    treat_single = .False.
!   if (lablinelo(ii) .eq. 'HE21  '.and. lablineup(ii) .eq. 'HE22  ')
!   treat_single=.true.
    If (opt_hei_sing) Then
      If (indexcmf1(ii)==3 .And. (lablinelo(ii)=='HE11S1 ' .And. xlam>584. &
        .And. xlam<585.)) treat_single = .True.
    End If

!   if (indexcmf1(ii).eq.2) then
!   if HeI should be treated in single line approx.
    If (indexcmf1(ii)==2 .Or. (treat_single .And. indexcmf1(ii)/=1)) Then

!     similar as in routine cmf (nlte.f90)
!     consistency check
      If (indexcmf(ii)/=2) Then
        Write (999, *) ' stop: indexcmf1=2 and indexcmf ne 2 in calc_tlu'
        Stop ' indexcmf1=2 and indexcmf ne 2 in calc_tlu'
      End If
      vdop = vdopcmf(ii)
      If (vdop==0.D0) Then
        Write (999, *) ' stop: error in vthcmf'
        Stop ' error in vthcmf'
      End If
!     no overlap for lines outside overlap range
      If (indov(j)/=0) Then
        Write (999, *) ' stop: indov ne 0 in calc_tlu. Set optsing = .true.'
        Stop ' indov ne 0 in calc_tlu. Set optsing = .true.'
      End If

      Call cmfsing(ii, r, v, dvdr, temp, xlam, vdop, vmax, teff, corrfc, &
        ajlmean, alo, inversion)

      Write (*, Fmt=100) xlam, lablinelo(ii), lablineup(ii)

      Do ll = 1, nd
        snew = ajlmean(ll) - alo(ll)*tulmat(ii, ll)
        If (snew<0.D0 .Or. abs(alo(ll))>1.D0) Then
          If (.Not. inversion) Then
            Write (999, *) ' warning!!!! something wrong with line-alo'
            Write (999, *) ' stop: problems (2) with snew in calc_tlu. &
              &Should never happen!!!'
            Print *, ' warning!!!! something wrong with line-alo'
            Stop ' problems (2) with snew in calc_tlu. Should never happen!!!'
          End If
          If (inversion) Then
            Write (999, *) &
              ' warning!!!! jbar corrected (calc_tlu, inversion) at', ll
            Print *, ' warning!!!! jbar corrected (calc_tlu, inversion) at', &
              ll
!           JO at some point, this needs to be improved (also in nlte.f90)
            ajlmean(ll) = 0.
          End If
          alo(ll) = 0.D0
          snew = ajlmean(ll)
        End If

!       if (lablinelo(ii) .eq. 'HE21  '.and. lablineup(ii) .eq. 'HE22  ') &
!       write(*,222) ll,ajlmean(ll),alo(ll),tulmat(ii,ll),snew

!       if(lablinelo(ii) .eq. 'HE11S1  '.and.lablineup(ii) .eq. 'HE12P1  ') &
!       print*,ll,ajlmean(ll),tulmat(ii,ll),snew


        tlumat(ii, ll) = blucmf(ii)*snew
        tulmat(ii, ll) = bulcmf(ii)*snew + aulcmf(ii)*(1.D0-alo(ll))
!       if(ll.eq.1 .or. ll.eq.10 .or. ll.eq.30 .or. ll.eq. 40 .or. ll.eq.50) &
!       &              print*,ii,ll,tlumat(ii,ll),tulmat(ii,ll)
      End Do

    End If
  End Do

  Print *
  If (.Not. no_alo) Then
    Print *, ' RBB-RATES FOR EXPLICIT ELEMENTS PREPARED, ALO NE 0 (CALC_TLU)'
  Else
    Print *, ' RBB-RATES FOR EXPLICIT ELEMENTS PREPARED, ALO = 0 (CALC_TLU)'
  End If
  Print *

  Return

100 Format ('CMFT AT', F12.3, ' A ', A6, ' TO ', A6, ' SINGLE LINE')
110 Format ('CMFT AT', F12.3, ' A ', A6, ' TO ', A6, ' MULTI LINE (PACKED)')
120 Format ('CMFT AT', F12.3, ' A ', A6, ' TO ', A6, &
    ' MULTI LINE (DEPACKED, modified wavelengths)')
130 Format ('SA-T AT', F12.3, ' A ', A6, ' TO ', A6, ' FOR ', I3, ' DEPTH', &
    ' POINTS')
End Subroutine

!----------------------------------------------------------------------------

Subroutine calc_tlu_met(nd, clf)

! calculates reduced RBB rate_coefficients for selected bg elements
! packing of individual components as follows (with slinei the
! source-function of the individual components and sline the source function
! of the packed transition: slinei=lam^3/lami^3*sline

! jbartot = sum gfi*jbari / sum gfi (see thesis Jo)

! =>

! tlu=snew = jbartot - alotot*sline = sum gfi*(jbari-aloi*slinei) / sum gfi =
! jbartot - sum (aloi*lam^3/lam_i^3 sline) / sum gfi

! =>

! snue=jbartot - alo1tot*sline    and  tul=alotot
! with jbartot as above, alo1tot = sum (gfi*aloi*lam^3/lami^3)/sum gfi,
! alotot= sum (gfi*alo)/sum gfi


  Use :: nlte_type
  Use :: nlte_dim, Only: id_ndept
  Use :: fund_const, Only: hc2, c_tau
  Use :: fastwind_params, Only: natom
  Use :: nlte_var, Only: lwion1, sr, vmax, no_alo_metals
  Use :: nlte_lines, Only: unpacked_index_cmf, unpacked_gf, unpacked_xlam
  Use :: nlte_app, Only: met_imin, met_imax, met_rbb, jatom_full, indexel, &
    indrec, irbb, ilin, occngold, tlu_met, tul_met
  Use :: cmf_all_var, Only: nf_cmf, elements_cmf, pweightdop, profdop, xj_cmf, &
    alo_cmf, sumopal_cmf
  Use :: ng_var, Only: ng_met, itng_met, sl_ng_met
  Use :: tcorr_var, Only: temp_converged

  Implicit None

  Integer (i4b), Intent (In) :: nd
  Real (dp), Dimension (nd), Intent (In) :: clf

  Real (dp), Parameter :: hc2_8 = hc2*1.D24 ! conversion to Angstrom

  Real (dp), Dimension (id_ndept) :: opalfac, sline, jbar_tot, alo_tot, &
    alo1_tot
! required field-length only LWION1-1

  Real (dp), Dimension (0:60) :: pw1, pw2

  Real (dp) :: srvmax, c_tau1, xlamc, lam3, occlowi, occupi, gfi, lam, const, &
    opal, jbar, alobar, sumgf, snew, rho1, lam3i, ratio, dum, sl

  Real (dp) :: sumrho, const1, gfip, lamp

  Integer (i4b) :: k, j, ilow, imax, n, jj, ii, ji, iqual, no_comp, &
    index_comp, in, lw1, ixmax, m1, m2, ll, il, ik

  Integer (i4b) :: iov, deltai, instart, inend, il1, ikp

  Logical, Dimension (id_ndept) :: calc, inv

! required field-length only LWION1-1

  If (maxval(nf_cmf)>60) Then
    Write (999, *) ' stop: dimension (pw1, pw2) > 60 in calc_tlu_met, adjust!'
    Stop ' dimension (pw1, pw2) > 60 in calc_tlu_met, adjust!'
  End If

  srvmax = sr/vmax
  c_tau1 = c_tau*srvmax

  lw1 = lwion1 - 1

  tlu_met = -9999.D0
  tul_met = -9999.D0

  Do k = 1, natom
    If (jatom_full(k)==0) Cycle

!   explicit elements will be treated by conventional (simplified) approach
!   in netmat_met: since atomic data are different, particularly line-
!   frequencies, and the cmf-transport is performed with the 'explicit' data,
!   Jbars calculated using bg-data would be erroneous anyway.
    If (indexel(k)>0) Cycle

    If (k<3) Then
      Write (999, *) ' stop: k < 3 in calc_tlu_met'
      Stop ' k < 3 in calc_tlu_met'
    End If

    ixmax = nf_cmf(k)

!   if(ixmax.eq.0) stop 'ixmax = 0 in calc_tlu_met'
!   change July2014. Consider only elements which are included in elements_cmf
    If (ixmax==0) Then
      If (elements_cmf(k)/=0) Then
        Write (999, *) &
          ' stop: ixmax = 0 and elements_cmf(k) ne 0 in calc_tlu_met'
        Stop ' ixmax = 0 and elements_cmf(k) ne 0 in calc_tlu_met'
      End If
      Cycle
    End If

    ilow = minval(met_imin(k,1:lw1))
    imax = maxval(met_imax(k,1:lw1))

!   we consider all ionization stages present between 1 to LWION1-1
    Do j = ilow, imax

      n = indrec(k, j)
      jj = irbb(n) - 1

      Do ii = 1, ilin(n)
        ji = jj + ii
        iqual = met_rbb(ji)%quality
        no_comp = met_rbb(ji)%no_comp
        If (iqual==2 .And. no_comp/=0) Then
          Write (999, *) ' stop: iqual = 2 and no_comp ne 0'
          Stop ' iqual = 2 and no_comp ne 0'
        End If
!       outside range or no individual component
!       will be treated by conventional (simplified) approach in netmat_met
        If (iqual==-1 .Or. iqual==2) Cycle

        m1 = met_rbb(ji)%low
        m2 = met_rbb(ji)%lup
!       test_oxy=.false.
!       if(k.eq.8.and.j.eq.3.and.m1.eq.1.and.(m2.eq.10.or.m2.eq.5.or.m2.eq.6))
!       test_oxy=.true.
        index_comp = met_rbb(ji)%index_comp
        xlamc = met_rbb(ji)%wave !  packed transition, in A
        lam3 = xlamc**3

!       depth dependent quantities which are identical for all components
        Do ll = 1, lw1

          occlowi = occngold(n, m1, ll) ! occng = n/g
          If (j<met_imin(k,ll) .Or. j>met_imax(k,ll) .Or. occlowi==0.) Then
!           needn't to be considered or new ion
            calc(ll) = .False.
          Else
            occupi = occngold(n, m2, ll)
!           following two quantities are identical for all components
            opalfac(ll) = (occlowi-occupi)/clf(ll)
            sline(ll) = hc2_8/lam3/(occlowi/occupi-1.D0)
!           note: calculated with packed wavelength, corrected below
            If (opalfac(ll)<=0.) Then
              inv(ll) = .True.
            Else
              inv(ll) = .False.
            End If
            calc(ll) = .True.
          End If
        End Do !                    first depth loop

        jbar_tot = 0.
        alo_tot = 0.
        alo1_tot = 0.

!       loop over all components
        Do in = index_comp, index_comp + no_comp - 1
!         this is the depth-independent part

!         print*,k,j,ii,no_comp,'comp=',in
          il = unpacked_index_cmf(in)
          gfi = unpacked_gf(in)
!         lam=lam_cmf(il) ! in A (rounded)
          lam = unpacked_xlam(in) ! in A
!         correction for difference in sline_i and sline
!         sline_i = lam^3/lami^3*sline
          lam3i = lam3/lam**3

!         overlapping components within profile
          Do iov = in, index_comp, -1
            deltai = unpacked_index_cmf(iov) - il
!           for tests:
            If (deltai>0) Then
              Write (999, *) ' stop: deltai > 0 (calc_tlu_met)'
              Stop ' deltai > 0 (calc_tlu_met)'
            End If
            If (abs(deltai)<=ixmax) instart = iov
          End Do
          Do iov = in, index_comp + no_comp - 1
            deltai = unpacked_index_cmf(iov) - il
!           for tests:
            If (deltai<0) Then
              Write (999, *) ' stop: deltai < 0 (calc_tlu_met)'
              Stop ' deltai < 0 (calc_tlu_met)'
            End If
            If (deltai<=ixmax) inend = iov
          End Do
!         hack to avoid coupled ALOs
!         instart=in
!         inend=in


!         for tests; note that there are certain ions that do not have
!         overlapping lines according to line-list, e.g., PVI and SVII
!         if(instart.ne.inend) then
!         print*,k,j,ii,lam
!         print*,no_comp,'comp=',in,instart,inend
!         endif

!         if(lam.gt.623.and.lam.lt.634.and.gfi.gt.0.05)
!         print*,lam,k,j,m1,m2,gfi

          const = c_tau1*lam*gfi
!         this is the depth-dependent part

          Do ll = 1, lw1
            If (.Not. calc(ll)) Cycle

!           pi e2/me c * SR/vmax * lam * gf * (nl/gl - nu/gu) /clf
            opal = const*opalfac(ll)

!           calculate reduced rates
!           ALO for single line is int alo_nu rho1_nu phi dnu, with rho1_nu =
!           opal * phi/opal_tot_nu
!           note that alo_nu already includes rho_nu =
!           chi_nu(line_tot)/chi_nu_tot

!           integration over phi
            jbar = 0.
            alobar = 0.

            If (inv(ll)) Then !     inversion, don't use alo
!             alobar remains zero
              pw1(0:ixmax) = pweightdop(0:ixmax, ll, k)

!             red side including center
              Do ik = 0, ixmax
                jbar = jbar + xj_cmf(ll, il+ik)*pw1(ik)
              End Do
!             blue side
              Do ik = 1, ixmax
                jbar = jbar + xj_cmf(ll, il-ik)*pw1(ik)
              End Do
!             if(test_oxy) print*,'alobari_O',lam,jbar,alobar

!             standard path, no inversion
            Else
              pw1(0:ixmax) = pweightdop(0:ixmax, ll, k)
              pw2(0:ixmax) = profdop(0:ixmax, ll, k)

!             red side including center
              Do ik = 0, ixmax
                jbar = jbar + xj_cmf(ll, il+ik)*pw1(ik)
                dum = sumopal_cmf(ll, il+ik)
                If (dum==0.D0) Then
!                 new ion, no sumopal has been calculated before
                  rho1 = 0.D0
                Else

                  If (instart==inend) Then ! one component only
                    rho1 = opal*pw2(ik)/dum
                  Else
                    sumrho = opal*pw2(ik)
                    const1 = 0.D0
                    Do iov = instart, inend
                      If (iov==in) Cycle
                      il1 = unpacked_index_cmf(iov)
                      ikp = abs(ik+il-il1) ! other component's profile index
                      If (ikp>ixmax) Cycle ! outside other's component range
                      gfip = unpacked_gf(iov)
                      lamp = unpacked_xlam(iov) ! in A
!                     correction for (slightly) different frequencies in SL:
!                     as calculated ALO refers to SL(iov);
!                     however, we want to apply it to SL(in). Thus:
                      const1 = const1 + gfip*lamp*pw2(ikp)*(lam/lamp)**3
                    End Do
                    sumrho = sumrho + const1*c_tau1*opalfac(ll)
                    rho1 = sumrho/dum
                  End If

                End If
!               if(test_oxy) print*,lam,gfi,ik,rho1
!               if(k.eq.18.and.j.eq.5.and.m1.eq.2.and.m2.eq.10 .and. ll.ge.9
!               .and. ll.le.13) then
!               print*,ll,ik,il+ik,lam,xlamc,rho1,opal*pw2(ik),sumopal_cmf(ll,il+ik),
!               &
!               &              xj_cmf(ll,il+ik),alo_cmf(ll,il+ik),pw1(ik)
!               endif

!               significant or dominating inversion of bg opacities
                If (rho1>1.1D0 .Or. rho1<0.D0) rho1 = 0.D0
!               slight difference in lambda (grid vs. exact) or mild inversion
!               of bg
                If (rho1>1.D0) rho1 = 1.D0
                alobar = alobar + alo_cmf(ll, il+ik)*rho1*pw1(ik)
              End Do

!             blue side
              Do ik = 1, ixmax
                jbar = jbar + xj_cmf(ll, il-ik)*pw1(ik)
                dum = sumopal_cmf(ll, il-ik)
                If (dum==0.D0) Then
!                 new ion, no sumopal has been calculated before
                  rho1 = 0.D0
                Else

                  If (instart==inend) Then ! one component only
                    rho1 = opal*pw2(ik)/dum
                  Else
                    sumrho = opal*pw2(ik)
                    const1 = 0.D0
                    Do iov = instart, inend
                      If (iov==in) Cycle
                      il1 = unpacked_index_cmf(iov)
                      ikp = abs(-ik+il-il1) ! other component's profile index
                      If (ikp>ixmax) Cycle ! outside other's component range
                      gfip = unpacked_gf(iov)
                      lamp = unpacked_xlam(iov) ! in A
!                     correction for (slightly) different frequencies in SL:
!                     as calculated ALO refers to SL(iov);
!                     however, we want to apply it to SL(in). Thus:
                      const1 = const1 + gfip*lamp*pw2(ikp)*(lam/lamp)**3
                    End Do
                    sumrho = sumrho + const1*c_tau1*opalfac(ll)
                    rho1 = sumrho/dum
                  End If
                End If
!               if(test_oxy) print*,lam,gfi,ik,rho1
!               if(k.eq.18.and.j.eq.5.and.m1.eq.2.and.m2.eq.10 .and. ll.ge.9
!               .and. ll.le.13) then
!               print*,ll,ik,il-ik,lam,xlamc,rho1,opal*pw2(ik),sumopal_cmf(ll,il-ik),
!               &
!               &              xj_cmf(ll,il-ik),alo_cmf(ll,il-ik),pw1(ik)
!               endif


!               significant or dominating inversion of bg opacities
                If (rho1>1.1D0 .Or. rho1<0.D0) rho1 = 0.D0
!               slight difference in lambda (grid vs. exact) or mild inversion
!               of bg
                If (rho1>1.D0) rho1 = 1.D0
                alobar = alobar + alo_cmf(ll, il-ik)*rho1*pw1(ik)
              End Do
!             if(test_oxy) print*,'alobari_O',lam,jbar,alobar

            End If !                inversion or not

!           sum up jbar and alo for all components
            jbar_tot(ll) = jbar_tot(ll) + jbar*gfi
            alo_tot(ll) = alo_tot(ll) + alobar*gfi
            alo1_tot(ll) = alo1_tot(ll) + alobar*gfi*lam3i ! includes
!           correction for
!           sline_i

          End Do !                  depth-loop for component in

        End Do !                    component loop for component in with index
!       il

!       now, calculate actual jbar and alo by normalizing with sumgf and
!       sumgflam3 = sum (gfi * lam^3/lam_i^3), respecively
!       JO Dec 2016: sumgflam3 does not exist, and was not used in previous
!       version
!       of calc_tlu_met (approach changed from early on, and forgot to delete
!       comment?)

!       sum(gfi) = sumgf already checked in subr. fgrid_cmf
        sumgf = met_rbb(ji)%sum_gf
!       for all depth points
        jbar_tot = jbar_tot/sumgf
        alo_tot = alo_tot/sumgf
        alo1_tot = alo1_tot/sumgf

        If (itng_met==2) Then
!         actually, jbar corresponds to J(n-1) corresponding itng_met=3,
!         but itng_met is only updated (in subr. prep_ng_met) after present
!         routine
!         has been called. Thus, itng_met=2 still.
          ik = met_rbb(ji)%index1_ng
          If (ik/=0) sl_ng_met(:, ik, 5) = jbar_tot(1:lw1)
        End If

!       in case, use Ng-extrapolated source-function.
!       no rescaling of lambda required here (calculated with same xlamc)
        If (ng_met) Then
!         here it might happen that ng_met and metals_converged, if latter
!         was set to true in subr. Trad just before call to nlte_approx
!         Thus, NO similar check as in subr. opacity_cmf
          ik = met_rbb(ji)%index1_ng
          If (ik/=0) Then
            sl = sl_ng_met(35, ik, 4)
            If (sl>0.D0) sline(1:lw1) = sl_ng_met(:, ik, 4)
          End If
        End If

!       now, perform some test and, in case, redefine alo!
        Do ll = 1, lw1
          If (.Not. calc(ll)) Then
            If (jbar_tot(ll)/=0. .Or. alo_tot(ll)/=0.) Then
              Write (999, *) &
                ' stop: calc = .false. but jbar or alo ne 0 in calc_tlu_met'
              Stop ' calc = .false. but jbar or alo ne 0 in calc_tlu_met'
            End If
            tlu_met(ji, ll) = 0.
            tul_met(ji, ll) = 0.

          Else

            If (abs(alo_tot(ll))>1.D0) Then
!             print*,' warning!!!! something wrong with line-alo inq calc_tlu'
              Write (999, *) ' stop: warning!!!! something wrong &
                &with line-alo in calc_tlu_met'
              Stop &
                ' warning!!!! something wrong with line-alo in calc_tlu_met'
              alo_tot(ll) = 0.D0
            End If

!           here we have to consider the differences between packed und
!           individual
!           source function, thus working with alo1_tot
!           the alo itself, however, is independent of source function and
!           thus does
!           not need to be corrected
            snew = jbar_tot(ll) - alo1_tot(ll)*sline(ll)

!           for tests (no ALO)
!           snew = jbar_tot(ll)
!           alo_tot(ll) = 0.


!           if(test_oxy) then
!           print*,xlamc,ll
!           print*,jbar_tot(ll),alo_tot(ll),alo1_tot(ll),sline(ll),snew
!           endif


!           final test on snew (with potential correction of alo)
            If (snew<0.D0 .And. .Not. no_alo_metals) Then
!             JO March 2025 in the current version, this should never happen
!             JO May 2025: indeed, it *very* seldom happens, needs to be
!             investigated when time
!             stop ' problems with snew in calc_tlu_met. Should never happen!'
              If (no_comp==1 .Or. alo_tot(ll)<=0.99D0) Then
                Write (999, *)
                Write (999, *) &
                  ' warning!!!! something wrong with snew in calc_tlu_met'
                Write (999, Fmt='(f10.3,6(2x,i3))') lam, k, j, m1, m2, ll, &
                  no_comp
                Write (999, *) snew, jbar_tot(ll), alo1_tot(ll), sline(ll), &
                  alo_tot(ll)
                Print *
                Print *, &
                  ' warning!!!! something wrong with snew in calc_tlu_met'
                Write (*, Fmt='(f10.3,6(2x,i3))') lam, k, j, m1, m2, ll, &
                  no_comp
                Print *, snew, jbar_tot(ll), alo1_tot(ll), sline(ll), &
                  alo_tot(ll)
!               JO April 2017: stop statement removed: problems can happen
!               right after
!               hydro update (d10v with ND=51)
                snew = jbar_tot(ll)
                alo_tot(ll) = 0.
              Else
                Write (999, *)
                Write (999, *) &
                  ' warning!!!! problems with snew in calc_tlu_met'
                Write (999, Fmt='(f10.3,6(2x,i3))') lam, k, j, m1, m2, ll, &
                  no_comp
                Write (999, *) snew, jbar_tot(ll), alo1_tot(ll), sline(ll), &
                  alo_tot(ll)
                Print *
                Print *, ' warning!!!! problems with snew in calc_tlu_met'
                Write (*, Fmt='(f10.3,6(2x,i3))') lam, k, j, m1, m2, ll, &
                  no_comp
                Print *, snew, jbar_tot(ll), alo1_tot(ll), sline(ll), &
                  alo_tot(ll)
!               try to correct (might happen if alo very close to one)
                ratio = alo1_tot(ll)/alo_tot(ll)
                alo_tot(ll) = 0.99D0 ! underestimate allowed
                alo1_tot(ll) = ratio*0.99D0
                snew = jbar_tot(ll) - alo1_tot(ll)*sline(ll)
                If (snew<0.D0) Then
!                 print*,xlamc,k,j,ii,m1,m2,no_comp,ll,inv(ll)
!                 print*,jbar_tot(ll),alo_tot(ll),alo1_tot(ll),sline(ll),snew
!                 JO Sept. 2019: stop statement removed; problems can occur in
!                 first
!                 iterations of moments.
!                 stop ' problems cannot be corrected'
                  snew = jbar_tot(ll)
                  alo_tot(ll) = 0.
                End If
              End If
            End If

!           for tests
!           if (k.eq.8 .and. j.eq.3 .and. m1.eq.1 .and. m2.eq.10) then
!           n=indrec(8,3)
!           write(18,*) 'O',ll,occngold(n,1,ll),occngold(n,10,ll)
!           write(18,*) jbar,alobar,alo1bar,sline,snew
!           write(18,*) jbar_tot(ll),alo_tot(ll),alo1_tot(ll),sline(ll),snew
!           print*,'O3/1-10',ll,sline(ll),snew,alo_tot(ll)
!           endif


!           hack to avoid oscillations in strongly coupled lines
!           JO April 2018: changed from 0.8 to 0.9, since also for 0.9, there
!           is a reasonable convergence without ALO (n propto 1/beta =
!           1/(1-ALO)),
!           but we avoid certain convergence issues due to impact of
!           non-diagonal
!           elements which can be still substantial for ALO = 0.8
!           (remember Ar V in the outer wind of zeta Pup)

            If (no_alo_metals) Then
              snew = jbar_tot(ll)
              alo_tot(ll) = 0.D0
            Else If (alo_tot(ll)<0.9) Then
!             JO July 2018: Nevertheless, for O3 1-10 and N3 1-11 we keep the
!             ALO
              If (k==8 .And. j==3 .And. m1==1 .And. m2==10) Then
                Continue
              Else If (k==7 .And. j==3 .And. m1==1 .And. m2==11) Then
                Continue
              Else
                snew = jbar_tot(ll)
                alo_tot(ll) = 0.D0
              End If
            End If

!           JO June 2018; do not allow for too large ALOs. Can lead to
!           oscillations,
!           if slightly different from iteration to iteration, and ALO at
!           0.999...
!           but use only if temperature already converged
            If (temp_converged) Then
              If (alo_tot(ll)>0.99D0) Then
                ratio = alo1_tot(ll)/alo_tot(ll)
                alo_tot(ll) = 0.99D0 ! underestimate allowed
                alo1_tot(ll) = ratio*0.99D0
                snew = jbar_tot(ll) - alo1_tot(ll)*sline(ll)
              End If
            End If

!           for tests
            If (k==8 .And. j==3 .And. m1==1 .And. m2==10) Then
              n = indrec(8, 3)
              Write (18, *) 'O', ll, occngold(n, 1, ll), occngold(n, 10, ll)
!             write(18,*) jbar_tot(ll),alo_tot(ll),alo1_tot(ll),sline(ll),snew
!             print*,'O3/1-10',ll,sline(ll),snew,alo_tot(ll)
            End If

!           if(alo_tot(ll).gt.0.95.and.instart.ne.inend)
!           print*,k,j,ll,ii,alo_tot(ll),alo1_tot(ll)

!           tlu = blucmf(indmat)*snew
!           tul = bulcmf(indmat)*snew + aulcmf(indmat)* (1.d0-alobar)
            tlu_met(ji, ll) = snew
            tul_met(ji, ll) = alo_tot(ll)

!           for tests
!           if(k.eq.8.and.j.eq.5.and.m1.eq.1.and.m2.eq.3 .and. ll.ge.17 .and.
!           ll.le.25) then
!           print*,k,j,ii,m1,m2,no_comp,ll,inv(ll),jbar_tot(ll),alo_tot(ll),sline(ll),snew
!           endif

          End If !                  calc or .not. calc

        End Do !                    depth loop for packed line

      End Do !                      line loop

    End Do !                        ion loop

  End Do !                          atom loop

  Print *
  Print *, ' RBB-RATES FOR BACKGROUND ELEMENTS PREPARED (CALC_TLU_MET)'
  Print *

  If (ng_met) Then
!   in case that ALMOST_CONVERGED was set to T just before
    ng_met = .False.
    Print *, ' ng_met set to false'
    Print *
  End If

  Return

End Subroutine

!----------------------------------------------------------------------------

Subroutine interpol_cmf_to_coarse(nd)

! interpolates cmf-quantities (J,K) onto coarse grid, under the condition
! that integral values remain conserved.
! Note that integrations are performed over energy (Kayser) = 1/lambda
! IMPORTANT: inter- and extrapolation works only well when freq. grid
! is equidistant. To this end, subr. FRESCAL and FRESFIN have been
! modified at important edges (IMP=1).

  Use :: nlte_type
  Use :: nlte_dim, Only: id_ndept
  Use :: fastwind_params, Only: wavblue, wavred, optout_xj_xh
  Use :: nlte_var, Only: ifre, fre, xj, xj_cmf_coarse, xxk, alo, modnam, &
    kcmf_start, kcmf_end, no_check, xj_save
  Use :: cmf_all_var, Only: lam_cmf, nftot, xj_cmf

  Implicit None

  Integer (i4b), Intent (In) :: nd

  Integer (i4b) :: i, nmax, is, nmax_ind, il, ig, k, n1, n2, ks
  Real (dp) :: eblue, ered, l1, l2, wnue1, weight, xn0, xn1, xn2, xn3, deh, &
    de, errmax, eddftest

  Integer (i4b), Allocatable, Dimension (:) :: index
  Real (dp), Allocatable, Dimension (:) :: ener_cmf
  Real (dp), Allocatable, Dimension (:, :) :: integ_cmf, xint
! real (dp), allocatable, dimension(:,:), save :: xj_old

  Real (dp), Dimension (id_ndept) :: di, di1, i0, i1, i2, i3, xiquart, xihalf, &
    diff1, diff2

! interpol_lin included in contain

  eblue = 1.D8/wavblue
  ered = 1.D8/wavred

! for tests
! eblue=1.d8/240.
! ered=1.d8/8264.

! note that fre is in rydberg, and ordered from red to blue
! the index-file and lam_cmf is ordered from blue to red

  Do i = 1, ifre - 1
    If (fre(i)>=ered) Exit
  End Do

  nmax = i
  Do i = ifre - 1, nmax, -1
    If (fre(i)<=eblue) Exit
  End Do
  is = i

  If (abs(1.-fre(is)/eblue)>0.1 .Or. abs(1.-fre(nmax)/ered)>0.1) Then
    Write (999, *) wavblue, 1.D8/fre(is), wavred, 1.D8/fre(nmax)
    Write (999, *) ' stop: something wrong with freq. boundaries &
      &(subr. interpol_cmf_to_coarse)'
    Print *, wavblue, 1.D8/fre(is), wavred, 1.D8/fre(nmax)
    Stop &
      ' something wrong with freq. boundaries (subr. interpol_cmf_to_coarse)'
  End If

  nmax_ind = is - nmax + 1

  Allocate (index(0:nmax_ind-1))
  Allocate (integ_cmf(nd,nftot), ener_cmf(nftot))
  Allocate (xint(nd,ifre))

! xj_cmf_coarse should have been allocated in CMF_COMPLETE, under XJ_SMOOTHED
! in parallel, xj_save has been allocated.
  If (.Not. allocated(xj_cmf_coarse)) Then
    Write (999, *) &
      ' stop: xj_cmf_coarse not allocated in INTERPOL_CMF_TO_COARSE'
    Stop ' xj_cmf_coarse not allocated in INTERPOL_CMF_TO_COARSE'
  End If

  xj_cmf_coarse = 0.D0 !            overwrite previous values from XJ_SMOOTHED
  ener_cmf = 1.D8/lam_cmf

  index(0) = is
  index(1) = is
  ig = 2

  Do il = is - 1, nmax, -1
    l1 = fre(il+1)
    l2 = fre(il)

    If (abs(1.-l2/l1)<1.D-5) Then
      If (il==is-1) Then
        index(0) = is
        index(1) = is - 1
        ig = 2
      Else
!       if 3 frequency points very narrow, make new range
!       e.g., instead of (479, 480),(480,481) -> (479,480),(481,482)
        If (index(ig-1)/=il+1) Then
          index(ig) = il + 1
          index(ig+1) = il
!         write(*,fmt='(i4,2(2x,f10.4),3(2x,i4))') &
!         &             il,1.d8/l1,1.d8/l2,ig,index(ig),index(ig+1)
          ig = ig + 2
        End If
      End If
    End If
  End Do

  ig = ig - 1
  If (index(ig)/=nmax) Then
    index(ig+1) = nmax
    index(ig+2) = nmax
    ig = ig + 2
  End If

  If (modulo(ig,2)/=1) Then
    Write (999, *) ' ig even in subroutine interpol_cmf_to_coarse'
    Print *, ' ig even in subroutine interpol_cmf_to_coarse'
    Stop
  End If

! for tests
! do il=0,ig-1,2
! write(*,fmt='(3(i4,2x),2(f10.4),2x)') &
! il,index(il),index(il+1),1.d8/fre(index(il)),1.d8/fre(index(il+1))
! enddo

! preparation of energy integration weights (already checked in cmf_complete)
! integration weights from using trapezoidal rule over delta E, and
! delta E = 1/lam*(1-1/(1.+v/c))
! factor 0.5 (from trapezoidal rule) * 1.d8 included in weights
! Basically
! vc1=lam_cmf(2)/lam_cmf(1)-1.d0
! wnue1=(1.d0-1.d0/(1.d0+vc1))*0.5d8
! leads to
  wnue1 = (1.D0-lam_cmf(1)/lam_cmf(2))*0.5D8

  integ_cmf(:, 1) = 0.

! calculate integral over cmf-quantities (as a function of freq.)
  Do k = 2, nftot

    weight = wnue1/lam_cmf(k-1) !   corresponds to 0.5*(nu(k-1)-nu(k))
    di = (xj_cmf(:,k-1)+xj_cmf(:,k))*weight
    integ_cmf(:, k) = integ_cmf(:, k-1) + di
!   test for accuracy (truncation errors if J varies considerably with nu)
    di1 = integ_cmf(:, k) - integ_cmf(:, k-1)
!   JO July 2018 accuracy changed (1.d-5 to 1.d-3)
    If (maxval(abs(1.-di1/di))>1.D-3) Then
      Write (999, *) ' problem with accuracy of J-integration'
      Write (999, *) k, lam_cmf(k)
      Write (999, *) di
      Write (999, *)
      Write (999, *) di1
      Write (999, *)
      Write (999, *) 1. - di1/di
      Write (999, *) ' stop: problem with accuracy of J-integration'
      Print *, ' problem with accuracy of J-integration'
      Print *, k, lam_cmf(k)
      Print *, di
      Print *
      Print *, di1
      Print *
      Print *, 1. - di1/di
      Stop ' problem with accuracy of J-integration'
    End If

  End Do

! now, approximate quantities for coarse grid, so that integrals
! remain conserved
  ks = 1

  Do i = 1, ig - 2, 2
    n1 = index(i)
    n2 = index(i+1)

!   piecewise, from edge to edge
    Do k = n1, n2, -1

      If (k==n1) Then
!       standard approach and extrapolation (assuming linear function) at
!       first point
        xn1 = fre(k)
        xn2 = 0.5D0*(fre(k)+fre(k-1))
        xn3 = fre(k-1)
        i1 = interpol_lin(integ_cmf, ener_cmf, xn1, ks, nd, nftot)
        i2 = interpol_lin(integ_cmf, ener_cmf, xn2, ks, nd, nftot)
        i3 = interpol_lin(integ_cmf, ener_cmf, xn3, ks, nd, nftot)
        deh = xn1 - xn2
        xiquart = (i2-i1)/deh
        xj_cmf_coarse(:, k) = xiquart ! standard
        de = xn1 - xn3
        xihalf = (i3-i1)/de
        xint(:, k) = 2.D0*xiquart - xihalf ! linear extrapolation
!       print*,'k=n1',k,1.d8/xn1,1.d8/xn2,1.d8/xn3
!       print*,k,1.d8/fre(k),xj_cmf_coarse(1,k),xint(1,k)

      Else If (k==n2) Then
!       standard approach and extrapolation (assuming linear function) at last
!       point
        xn0 = fre(k+1)
        xn1 = xn2
        xn2 = fre(k)
!       restart of ks
        ks = 1
        i0 = interpol_lin(integ_cmf, ener_cmf, xn0, ks, nd, nftot)
        i1 = i2 !                   previously calculated
!       i1=interpol_lin(integ_cmf,ener_cmf,xn1,ks,nd,nftot)
        i2 = interpol_lin(integ_cmf, ener_cmf, xn2, ks, nd, nftot)
        deh = xn1 - xn2
        xiquart = (i2-i1)/deh
        xj_cmf_coarse(:, k) = xiquart ! standard
        de = xn0 - xn2
        xihalf = (i2-i0)/de
        xint(:, k) = 2.*xiquart - xihalf ! linear extrapolation
!       print*,'k=n2',k,1.d8/xn0,1.d8/xn1,1.d8/xn2
!       print*,k,1.d8/fre(k),xj_cmf_coarse(1,k),xint(1,k)

      Else
!       almost exact at intermediate points
        xn1 = xn2
        xn2 = 0.5*(fre(k)+fre(k-1))
        de = xn1 - xn2
!       i1=interpol_lin(integ_cmf,ener_cmf,xn1,ks,nd,nftot)
        i1 = i2 !                   previously calculated
!       ks continued
        i2 = interpol_lin(integ_cmf, ener_cmf, xn2, ks, nd, nftot)
        xj_cmf_coarse(:, k) = (i2-i1)/de
!       print*,k,fre(k+1),fre(k),fre(k-1),i1(1),i2(1),xj_cmf_coarse(1,k)
!       xint(:,k)=xj_cmf_coarse(k) ; not needed
!       print*,'else',k,1.d8/xn1,1.d8/xn2
!       print*,k,1.d8/fre(k),xj_cmf_coarse(1,k),' 0.'
      End If
    End Do
  End Do

! --------
! check of accuracy: until here, total integrals should be identical
  ered = fre(index(ig-1))
  Do i = nftot, 1, -1
    If (ener_cmf(i)>ered) Exit
  End Do
  nmax = i + 1
! print*,1.d8/ered,lam_cmf(nmax)

  eblue = fre(index(1))
  Do i = 1, nftot
    If (ener_cmf(i)<eblue) Exit
  End Do
  is = i - 1
! print*,1.d8/eblue,lam_cmf(is)

! this is the energy integral over the cmf-value (slightly larger)
  xiquart = (integ_cmf(:,nmax)-integ_cmf(:,is))

! now we calculate the energy integral over the exactly interpolated values,
! using the trapezoidal rule (wfre contains different weights, particularly at
! end)
  xihalf = 0.
  Do i = index(ig) + 1, index(0)
    xihalf = xihalf + 0.5D0*(xj_cmf_coarse(:,i)+xj_cmf_coarse(:,i-1))*(fre(i)- &
      fre(i-1))
  End Do

  diff1 = abs(1.-xihalf/xiquart)
  errmax = maxval(diff1)

! do i=1,nd
! print*,i,xiquart(i),xihalf(i)
! enddo

  Print *
  Print *, ' interpolation of J-CMF values onto coarse grid'
  Print *, ' maximum deviation between total integrals = ', errmax
  Print *

! if(errmax.gt.1.d-4) stop ' total integrals (cmf vs coarse) before
! extrapolation not conserved'
  If (errmax>1.D-2) Then
    Write (999, *) ' stop: total integrals (cmf vs coarse) before &
      &extrapolation not conserved'
    Stop ' total integrals (cmf vs coarse) before extrapolation not conserved'
  End If
! --------


! now we check which condition (standard or extrapol.) is better.
! standardwise, we use extrapol at the boundaries and standard else,
! except for the case that diff(extrapol.) < diff(standard) over the edges
! this allows for a smooth transition for smooth intensities
  Where (xint(:,index(1))>0.)
    xj_cmf_coarse(:, index(1)) = xint(:, index(1))
  End Where
  Where (xint(:,index(ig-1))>0.)
    xj_cmf_coarse(:, index(ig-1)) = xint(:, index(ig-1))
  End Where

  Do i = 3, ig - 2, 2
    n1 = index(i)
    diff1 = abs(xint(:,n1)-xint(:,n1+1))
    diff2 = abs(xj_cmf_coarse(:,n1)-xj_cmf_coarse(:,n1+1))
    Where (diff1<diff2 .And. xint(:,n1)>0. .And. xint(:,n1+1)>0.)
      xj_cmf_coarse(:, n1+1) = xint(:, n1+1)
      xj_cmf_coarse(:, n1) = xint(:, n1)
    End Where
  End Do

! specific treatment in case that first/last points are edges
! Check whether there are edges in XJ (Jnu(approx)).
! If so, modify xj_cmf_coarse(index(0)) and/or xj_cmf_coarse(index(ig)),
! by using the corresponding values from xj.
! Otherwise, use interpolated values at neighbouring frequencies
  If (index(0)/=index(1)) Then
    diff1 = abs(1.-xj(:,index(0))/xj(:,index(1)))
    Where (diff1>0.05)
      xj_cmf_coarse(:, index(0)) = xj(:, index(0))
    Elsewhere
      xj_cmf_coarse(:, index(0)) = xj_cmf_coarse(:, index(1))
    End Where
  End If

  If (index(ig)/=index(ig-1)) Then
    diff1 = abs(1.-xj(:,index(ig))/xj(:,index(ig-1)))
    Where (diff1>0.05)
      xj_cmf_coarse(:, index(ig)) = xj(:, index(ig))
    Elsewhere
      xj_cmf_coarse(:, index(ig)) = xj_cmf_coarse(:, index(ig-1))
    End Where
  End If

! 2nd check of accuracy, now for final values (not as exact, since linear
! extrapolation)
  ered = fre(index(ig))
  Do i = nftot, 1, -1
    If (ener_cmf(i)>ered) Exit
  End Do
  nmax = i + 1
! print*,1.d8/ered,lam_cmf(nmax)

  eblue = fre(index(0))
  Do i = 1, nftot
    If (ener_cmf(i)<eblue) Exit
  End Do
  is = i - 1
! print*,1.d8/eblue,lam_cmf(is)

! this is the energy integral over the cmf-value from eblue to ered (slightly
! larger)
  xiquart = (integ_cmf(:,nmax)-integ_cmf(:,is))

! now we calculate the energy integral over the interpolated value,
! using the trapezoidal rule (wfre contains different weights, particularly at
! end)
  xihalf = 0.
  Do i = index(ig) + 1, index(0)
    xihalf = xihalf + 0.5D0*(xj_cmf_coarse(:,i)+xj_cmf_coarse(:,i-1))*(fre(i)- &
      fre(i-1))
  End Do

  diff1 = abs(1.-xihalf/xiquart)
  errmax = maxval(diff1)

! do i=1,nd
! print*,i,xiquart(i),xihalf(i)
! enddo

  Print *, &
    ' maximum deviation between total integrals (incl. extrapolation) = ', &
    errmax
  Print *

! if(errmax.gt.0.01) stop ' too large deviation between total integrals
! (final)'
! JO changed Nov. 2019, to exclude problems in few specific cases
  If (errmax>0.05) Then
    Write (999, *) &
      ' stop: too large deviation between total integrals (final)'
    Stop ' too large deviation between total integrals (final)'
  End If

! finally, &
! we approximate XK in the cmf-region (from the coarse Eddington-factors),
! use approx. values for XJ in those regions where lam < wavblue or lam >
! wavred, ...

! JO: at some point, we might use XK values from the cmf-solution
! ... and overwrite present values (from CONT) with new ones (from CMF)

! for tests (if photo-rates shall be constant)
! if(.not.allocated(xj_old)) then
! allocate(xj_old(nd,ifre))
! xj_old=xj_cmf_coarse
! print*,'xj_old allocated'
! else
! xj_cmf_coarse=xj_old
! print*,'xj_old used'
! endif

! JO changed
! no where statement, since xj and xj_cmf_coarse of different dimensions:
! xj(nd,id_frec1),xj_cmf_coarse(nd,ifre)
  kcmf_start = 0
  Do k = 1, ifre
    If (xj_cmf_coarse(1,k)/=0.D0) Then
      If (kcmf_start==0) kcmf_start = k
      kcmf_end = k
!     print*,k,1.d8/fre(k),xxk(nd,k),xj(nd,k),xxk(nd,k)/xj(nd,k),xj_cmf_coarse(nd,k)
      xxk(:, k) = xxk(:, k)/xj(:, k)*xj_cmf_coarse(:, k) ! approximate and
!     overwrite
      alo(:, k) = 0.D0 !            reset alo in line region
      xj(:, k) = xj_cmf_coarse(:, k) ! overwrite
    End If
!   JO March 2018
!   test only for lambda > 20 (for higher freq. -- X-rays --), inconsistencies
!   possible
    If (1.D8/fre(k)>20.) Then
      eddftest = xxk(nd, k)/xj(nd, k)
!     print*,k,1.d8/fre(k),xxk(nd,k),xj(nd,k),eddftest
!     if(.not. no_check .and. (eddftest.lt.0.31 .or. eddftest.gt.0.35)) &
!     &     stop ' inconsistent Edd. factor in interpol_cmf_to_coarse'
      If (.Not. no_check .And. (eddftest<0.31 .Or. eddftest>0.35)) Then
        Write (999, *) ' Warning!!! Warning!!! Warning!!!'
        Write (999, *) 1.D8/fre(k), ' ', eddftest
        Write (999, *) ' inconsistent Edd. factor in interpol_cmf_to_coarse'
        Write (999, *)
        Print *, ' Warning!!! Warning!!! Warning!!!'
        Print *, 1.D8/fre(k), ' ', eddftest
        Print *, ' inconsistent Edd. factor in interpol_cmf_to_coarse'
        Print *
      End If
    End If
  End Do


! JO changed; xj_cmf_coarse only filled in CMF_COMPLETE range

  xj(1, ifre+1) = 2 !               indicates cmf treatment

! needed for cmf_complete
  Do k = 1, ifre
    xj_save(:, k) = xj(:, k)
  End Do

  Print *, ' XJ and XK (coarse grid, from CONT) overwritten by CMF-values!,'
  Print *, ' in the range ', 1.D8/fre(kcmf_end), ' to ', 1.D8/fre(kcmf_start)
  Print *

  If (1.D8/fre(kcmf_end+1)>wavblue) Then
    Write (999, *) ' stop: problems in cmf-range (blue)'
    Stop ' problems in cmf-range (blue)'
  End If

  If (1.D8/fre(kcmf_start-1)<wavred) Then
    Write (999, *) ' stop: problems in cmf-range (red)'
    Stop ' problems in cmf-range (red)'
  End If

  If (optout_xj_xh) Then
    Open (1, File=trim(modnam)//'/out_xj_cmf_coarse')
    Do k = 1, ifre
      Write (1, 100) 1.D8/fre(k), xj(1, k), xj(43, k), xj(nd, k)
    End Do
    Close (1)
  End If

! do k=1,ifre
! errmax=maxval(alo(:,k))
! write(*,110) 1.d8/fre(k),xj(43,k),xj_cmf_coarse(43,k),alo(43,k),errmax
! enddo

! for tests
! do i=ifre,1,-1
! write(*,fmt='(i3,2x,f20.4,2x,4(e12.6,2x))')
! i,1.d8/fre(i),xj(1,i),xj_cmf_coarse(1,i),xj(51,i),xj_cmf_coarse(51,i)
! enddo

  Deallocate (index, integ_cmf, ener_cmf, xint)

  Print *, ' J/K-CMF interpolated onto coarse grid (integrals conserved),'
  Print *, ' and previous values (from CONT) replaced'
  Print *

  Return

100 Format (F14.4, 3(2X,E10.4))
110 Format (F14.4, 2(2X,E10.4), 2(2X,F12.4))


Contains

  Function interpol_lin(y, x, x1, ks, nd, nf)
!   linear interpolation at x1 in y(:,i),y(:,i+1)
!   abscissa values of y are x, located in descending order

    Use :: nlte_type
    Use :: nlte_dim, Only: id_ndept

    Implicit None

    Real (dp), Dimension (id_ndept) :: interpol_lin

    Integer (i4b), Intent (In) :: nd, nf
    Integer (i4b), Intent (Inout) :: ks

    Real (dp), Dimension (id_ndept, nf), Intent (In) :: y
    Real (dp), Dimension (nf), Intent (In) :: x
    Real (dp), Intent (In) :: x1

    Integer (i4b) :: i
    Real (dp) :: q, q1

    Do i = ks, nf - 1
      If (x1<=x(i) .And. x1>x(i+1)) Go To 100
!     print*,i,x(i),x1,x(i+1)
    End Do

    Write (999, *) ' stop: x1 not found in x (interpol_lin)'
    Stop ' x1 not found in x (interpol_lin)'

100 ks = i !                        for next index
    q = (x1-x(i+1))/(x(i)-x(i+1))
    q1 = 1.D0 - q

    interpol_lin = q*y(:, i) + q1*y(:, i+1)

    Return

  End Function

End Subroutine

!----------------------------------------------------------------------------

Subroutine interpol_hcmf_to_hcoarse

! interpolates h_obs (high resol) onto coarse grid, under the condition
! that integral values remain conserved.
! Note that integrations are performed over energy (Kayser) = 1/lambda
! (programmed in analogy to INTERPOL_CMF_TO_COARSE

! Note that lam_cmf and integ_cmf are actually obs. frame quantities
  Use :: nlte_type
  Use :: fund_const, Only: clight, sigsb, pi
  Use :: fastwind_params, Only: wavblue, wavred, optout_xj_xh

  Use :: nlte_app, Only: teff
  Use :: nlte_var, Only: modnam, ifre, fre, xxh, xh_obs_coarse, kcmf_start, &
    kcmf_end, rtau23
  Use :: cmf_all_var, Only: lam_cmf, nftot, h_obs

  Implicit None

  Integer (i4b) :: i, nmax, is, nmax_ind, il, ig, k, n1, n2, ks, kcmf_start1, &
    kcmf_end1
  Real (dp) :: eblue, ered, l1, l2, wnue1, weight, xn0, xn1, xn2, xn3, deh, &
    de, fluxtot

  Integer (i4b), Allocatable, Dimension (:) :: index
  Real (dp), Allocatable, Dimension (:) :: ener_cmf
  Real (dp), Allocatable, Dimension (:) :: integ_cmf, xint, xh

  Real (dp) :: di, di1, i0, i1, i2, i3, xiquart, xihalf, diff1, diff2

! interpol_lin included in contain

  eblue = 1.D8/wavblue
  ered = 1.D8/wavred

! for tests
! eblue=1.d8/240.
! ered=1.d8/8264.

! note that fre is in rydberg, and ordered from red to blue
! the index-file and lam_cmf is ordered from blue to red

  Do i = 1, ifre - 1
    If (fre(i)>=ered) Exit
  End Do

  nmax = i
  Do i = ifre - 1, nmax, -1
    If (fre(i)<=eblue) Exit
  End Do
  is = i

  If (abs(1.-fre(is)/eblue)>0.1 .Or. abs(1.-fre(nmax)/ered)>0.1) Then
    Write (999, *) wavblue, 1.D8/fre(is), wavred, 1.D8/fre(nmax)
    Write (999, *) ' stop: something wrong with freq. boundaries &
      &(subr. interpol_hcmf_to_hcoarse)'
    Print *, wavblue, 1.D8/fre(is), wavred, 1.D8/fre(nmax)
    Stop ' something wrong with freq. boundaries (subr. interpol_hcmf_to_hcoars&
      &e)'
  End If

  nmax_ind = is - nmax + 1

  Allocate (index(0:nmax_ind-1))
  Allocate (integ_cmf(nftot), ener_cmf(nftot))
  Allocate (xint(ifre), xh(ifre))

! allocated only once, since routine only called once (UNASOL)
! usally allocated only one, since called only one (if UNASOL)
! but for tests, called more often. Thus
  If (.Not. allocated(xh_obs_coarse)) Allocate (xh_obs_coarse(nftot))

  xh = xxh(1, 1:ifre)
! rescale xh to be consistent with h_obs
! xxh = hnu*(r(1)/sr)^2
! hobs = hnu*(r(1)/srnom)^2 = hnu*rmax^2
! xh -> xxh * (sr/srnom)^2 = xxh/rtau23^2
  xh = xh/rtau23**2

  xh_obs_coarse = 0.D0
  ener_cmf = 1.D8/lam_cmf

  index(0) = is
  index(1) = is
  ig = 2

  Do il = is - 1, nmax, -1
    l1 = fre(il+1)
    l2 = fre(il)

    If (abs(1.-l2/l1)<1.D-5) Then
      If (il==is-1) Then
        index(0) = is
        index(1) = is - 1
        ig = 2
      Else
!       if 3 frequency points very narrow, make new range
!       e.g., instead of (479, 480),(480,481) -> (479,480),(481,482)
        If (index(ig-1)/=il+1) Then
          index(ig) = il + 1
          index(ig+1) = il
!         write(*,fmt='(i4,2(2x,f10.4),3(2x,i4))') &
!         &             il,1.d8/l1,1.d8/l2,ig,index(ig),index(ig+1)
          ig = ig + 2
        End If
      End If
    End If
  End Do

  ig = ig - 1
  If (index(ig)/=nmax) Then
    index(ig+1) = nmax
    index(ig+2) = nmax
    ig = ig + 2
  End If

  If (modulo(ig,2)/=1) Then
    Write (999, *) ' ig even in subroutine interpol_hcmf_to_hcoarse'
    Print *, ' ig even in subroutine interpol_hcmf_to_hcoarse'
    Stop
  End If

! for tests
! do il=0,ig-1,2
! write(*,fmt='(3(i4,2x),2(f10.4),2x)') &
! il,index(il),index(il+1),1.d8/fre(index(il)),1.d8/fre(index(il+1))
! enddo

! preparation of energy integration weights (already checked in cmf_complete)
! integration weights from using trapezoidal rule over delta E, and
! delta E = 1/lam*(1-1/(1.+v/c))
! factor 0.5 (from trapezoidal rule) * 1.d8 included in weights
! Basically
! vc1=lam_cmf(2)/lam_cmf(1)-1.d0
! wnue1=(1.d0-1.d0/(1.d0+vc1))*0.5d8
! leads to
  wnue1 = (1.D0-lam_cmf(1)/lam_cmf(2))*0.5D8

  integ_cmf(1) = 0.

! calculate integral over cmf-quantities (as a function of freq.)
  Do k = 2, nftot

    weight = wnue1/lam_cmf(k-1) !   corresponds to 0.5*(nu(k-1)-nu(k))
    di = (h_obs(k-1)+h_obs(k))*weight
    integ_cmf(k) = integ_cmf(k-1) + di
!   test for accuracy (truncation errors if H varies considerably with nu)
    di1 = integ_cmf(k) - integ_cmf(k-1)
    If (di==0.D0) Then
      If (di1/=0.D0) Then
        Print *, 'di = 0 and di1 ne 0 in subr. interpol_hcmf_to_hcoarse'
      End If
    Else
!     JO July 2018 accuracy changed (1.d-5 to 1.d-3)
      If (abs(1.-di1/di)>1.D-3) Then
        Write (999, *) ' problem with accuracy of H-integration'
        Write (999, *) i, lam_cmf(k)
        Write (999, *) di
        Write (999, *) di1
        Write (999, *) ' stop: problem with accuracy of H-integration'
        Print *, 'problem with accuracy of H-integration'
        Print *, i, lam_cmf(k)
        Print *, di
        Print *, di1
        Stop ' problem with accuracy of H-integration'
      End If
    End If
  End Do

! now, approximate quantities for coarse grid, so that integrals
! remain conserved
  ks = 1

  Do i = 1, ig - 2, 2
    n1 = index(i)
    n2 = index(i+1)

!   piecewise, from edge to edge
    Do k = n1, n2, -1

      If (k==n1) Then
!       standard approach and extrapolation (assuming linear function) at
!       first point
        xn1 = fre(k)
        xn2 = 0.5D0*(fre(k)+fre(k-1))
        xn3 = fre(k-1)
        i1 = interpol_lin(integ_cmf, ener_cmf, xn1, ks, nftot)
        i2 = interpol_lin(integ_cmf, ener_cmf, xn2, ks, nftot)
        i3 = interpol_lin(integ_cmf, ener_cmf, xn3, ks, nftot)
        deh = xn1 - xn2
        xiquart = (i2-i1)/deh
        xh_obs_coarse(k) = xiquart ! standard
        de = xn1 - xn3
        xihalf = (i3-i1)/de
        xint(k) = 2.D0*xiquart - xihalf ! linear extrapolation

      Else If (k==n2) Then
!       standard approach and extrapolation (assuming linear function) at last
!       point
        xn0 = fre(k+1)
        xn1 = xn2
        xn2 = fre(k)
!       restart of ks
        ks = 1
        i0 = interpol_lin(integ_cmf, ener_cmf, xn0, ks, nftot)
        i1 = i2 !                   previously calculated
!       i1=interpol_lin(integ_cmf,ener_cmf,xn1,ks,nftot)
        i2 = interpol_lin(integ_cmf, ener_cmf, xn2, ks, nftot)
        deh = xn1 - xn2
        xiquart = (i2-i1)/deh
        xh_obs_coarse(k) = xiquart ! standard
        de = xn0 - xn2
        xihalf = (i2-i0)/de
        xint(k) = 2.*xiquart - xihalf ! linear extrapolation

      Else
!       almost exact at intermediate points
        xn1 = xn2
        xn2 = 0.5*(fre(k)+fre(k-1))
        de = xn1 - xn2
!       i1=interpol_lin(integ_cmf,ener_cmf,xn1,ks,nftot)
        i1 = i2 !                   previously calculated
!       ks continued
        i2 = interpol_lin(integ_cmf, ener_cmf, xn2, ks, nftot)
        xh_obs_coarse(k) = (i2-i1)/de
      End If
    End Do
  End Do

! --------
! check of accuracy: until here, total integrals should be identical
  ered = fre(index(ig-1))
  Do i = nftot, 1, -1
    If (ener_cmf(i)>ered) Exit
  End Do
  nmax = i + 1
! print*,1.d8/ered,lam_cmf(nmax)

  eblue = fre(index(1))
  Do i = 1, nftot
    If (ener_cmf(i)<eblue) Exit
  End Do
  is = i - 1
! print*,1.d8/eblue,lam_cmf(is)

! this is the energy integral over the obs-value (slightly larger)
  xiquart = (integ_cmf(nmax)-integ_cmf(is))

! now we calculate the energy integral over the exactly interpolated values,
! using the trapezoidal rule (wfre contains different weights, particularly at
! end)
  xihalf = 0.
  Do i = index(ig) + 1, index(0)
    xihalf = xihalf + 0.5D0*(xh_obs_coarse(i)+xh_obs_coarse(i-1))*(fre(i)-fre( &
      i-1))
  End Do

  diff1 = abs(1.-xihalf/xiquart)

! print*,xiquart,xihalf

  Print *
  Print *, ' interpolation of H-obs values onto coarse grid'
  Print *, ' maximum deviation between total integrals = ', diff1
  Print *

  If (diff1>1.D-4) Then
    Write (999, *) ' stop: total H-integrals (fine vs coarse) before &
      &extrapolation not conserved'
    Stop &
      ' total H-integrals (fine vs coarse) before extrapolation not conserved'
  End If
! --------


! now we check which condition (standard or extrapol.) is better.
! standardwise, we use extrapol at the boundaries and standard else,
! except for the case that diff(extrapol.) < diff(standard) over the edges
! this allows for a smooth transition for smooth intensities
  If (xint(index(1))>0.) xh_obs_coarse(index(1)) = xint(index(1))
  If (xint(index(ig-1))>0.) xh_obs_coarse(index(ig-1)) = xint(index(ig-1))

  Do i = 3, ig - 2, 2
    n1 = index(i)
    diff1 = abs(xint(n1)-xint(n1+1))
    diff2 = abs(xh_obs_coarse(n1)-xh_obs_coarse(n1+1))
    If (diff1<diff2 .And. xint(n1)>0. .And. xint(n1+1)>0.) Then
      xh_obs_coarse(n1+1) = xint(n1+1)
      xh_obs_coarse(n1) = xint(n1)
    End If
  End Do

! specific treatment in case that first/last points are edges
! Check whether there are edges in XH (Hnu(approx)).
! If so, modify xh_obs_coarse(index(0)) and/or xh_cmf_coarse(index(ig)),
! by using the corresponding values from XH.
! Otherwise, use interpolated values at neighbouring frequencies
  If (index(0)/=index(1)) Then
    diff1 = abs(1.-xh(index(0))/xh(index(1)))
    If (diff1>0.05) Then
      xh_obs_coarse(index(0)) = xh(index(0))
    Else
      xh_obs_coarse(index(0)) = xh_obs_coarse(index(1))
    End If
  End If

  If (index(ig)/=index(ig-1)) Then
    diff1 = abs(1.-xh(index(ig))/xh(index(ig-1)))
    If (diff1>0.05) Then
      xh_obs_coarse(index(ig)) = xh(index(ig))
    Else
      xh_obs_coarse(index(ig)) = xh_obs_coarse(index(ig-1))
    End If
  End If

! 2nd check of accuracy, now for final values (not as exact, since linear
! extrapolation)
  ered = fre(index(ig))
  Do i = nftot, 1, -1
    If (ener_cmf(i)>ered) Exit
  End Do
  nmax = i + 1
! print*,1.d8/ered,lam_cmf(nmax)

  eblue = fre(index(0))
  Do i = 1, nftot
    If (ener_cmf(i)<eblue) Exit
  End Do
  is = i - 1
! print*,1.d8/eblue,lam_cmf(is)

! this is the energy integral over the cmf-value from eblue to ered (slightly
! larger)
  xiquart = (integ_cmf(nmax)-integ_cmf(is))

! now we calculate the energy integral over the interpolated value,
! using the trapezoidal rule (wfre contains different weights, particularly at
! end)
  xihalf = 0.
  Do i = index(ig) + 1, index(0)
    xihalf = xihalf + 0.5D0*(xh_obs_coarse(i)+xh_obs_coarse(i-1))*(fre(i)-fre( &
      i-1))
  End Do

  diff1 = abs(1.-xihalf/xiquart)

! print*,xiquart,xihalf

  Print *, &
    ' maximum deviation between total H-integrals (incl. extrapolation) = ', &
    diff1
  Print *

  If (diff1>0.01) Then
    Write (999, *) &
      ' stop: too large deviation between total H-integrals (final)'
    Stop ' too large deviation between total H-integrals (final)'
  End If

! now, add approximate fluxes in region outside CMF transfer
  kcmf_start1 = 0
  Do k = 1, ifre
    If (xh_obs_coarse(k)/=0.D0) Then
      If (kcmf_start1==0) kcmf_start1 = k
      kcmf_end1 = k
    Else
      xh_obs_coarse(k) = xh(k) !    remember, xh rescaled to (r(1)/srnom)^2
    End If
  End Do

! JO Nov 2025: under specific conditions, kcmf_start1 (from H) and kcmf_start
! (from J) can be (slightly) different, because kcmf_start1 is calculated
! after a transformation from cmf to obs! Thus, H_obs might become zero since
! the obs. frame freq. might not be covered by the corresponding cmf freq.
! Similar, for kcmf_end.
  If (kcmf_start1/=kcmf_start) Then
    Write (999, *) kcmf_start, kcmf_start1
    Write (999, *) ' WARNING!!! problems in kcmf_start1'
    Print *, kcmf_start, kcmf_start1
    Print *, ' WARNING!!! problems in kcmf_start1'
    If (abs(kcmf_start-kcmf_start1)>=5) Then
      Write (999, *) ' stop: problems in kcmf_start1'
      Stop ' problems in kcmf_start1'
    End If
  End If

  If (kcmf_end1/=kcmf_end) Then
    Write (999, *) kcmf_end, kcmf_end1
    Write (999, *) ' WARNING!!! problems in kcmf_end1'
    Print *, kcmf_end, kcmf_end1
    Print *, ' WARNING!!! problems in kcmf_end1'
    If (abs(kcmf_end-kcmf_end1)>=5) Then
      Write (999, *) ' stop: problems in kcmf_end1'
      Stop ' problems in kcmf_end1'
    End If
  End If

  fluxtot = sigsb*teff**4

  xihalf = 0.
! do i=ifre-1,1,-1
! JO changed May 2025
  Do i = ifre, 2, -1
    xihalf = xihalf + 0.5D0*(xh_obs_coarse(i)+xh_obs_coarse(i-1))*(fre(i)-fre( &
      i-1))
  End Do

  xihalf = xihalf*4.D0*pi*clight

  diff1 = (1.-xihalf/fluxtot)
  Print *, ' maximum deviation flux-integral at Rmax (coarse grid)'
  Print *, ' and nominal flux = ', diff1, ' (positive: lower than nominal)'
  Print *

  If (optout_xj_xh) Then
    Open (1, File=trim(modnam)//'/out_xh_obs_coarse')
    Do k = 1, ifre
      Write (1, 100) 1.D8/fre(k), xh_obs_coarse(k)
    End Do
    Close (1)
  End If

  Deallocate (index, integ_cmf, ener_cmf, xint, xh)

  Print *, ' H-obs remapped onto coarse grid,'
  Print *, ' in the range ', 1.D8/fre(kcmf_end), ' to ', 1.D8/fre(kcmf_start)
  Print *

  Return

100 Format (F14.4, 3(2X,E10.4))


Contains

  Function interpol_lin(y, x, x1, ks, nf)
!   linear interpolation at x1 in y(i),y(i+1)
!   abscissa values of y are x, located in descending order

    Use :: nlte_type

    Implicit None

    Real (dp) :: interpol_lin

    Integer (i4b), Intent (In) :: nf
    Integer (i4b), Intent (Inout) :: ks

    Real (dp), Dimension (nf), Intent (In) :: y
    Real (dp), Dimension (nf), Intent (In) :: x
    Real (dp), Intent (In) :: x1

    Integer (i4b) :: i
    Real (dp) :: q, q1

    Do i = ks, nf - 1
      If (x1<=x(i) .And. x1>x(i+1)) Go To 100
!     print*,i,x(i),x1,x(i+1)
    End Do

    Write (999, *) ' stop: x1 not found in x (interpol_lin)'
    Stop ' x1 not found in x (interpol_lin)'

100 ks = i !                        for next index
    q = (x1-x(i+1))/(x(i)-x(i+1))
    q1 = 1.D0 - q

    interpol_lin = q*y(i) + q1*y(i+1)

    Return

  End Function

End Subroutine

!----------------------------------------------------------------------------

Subroutine prep_smooth

  Use :: nlte_type

  Use :: nlte_dim, Only: id_ndept

  Use :: fastwind_params, Only: wavblue, wavred

  Use :: nlte_var, Only: ifre, ifreold1, fre, freold, xj_save, &
    xj_smoothed => xj_cmf_coarse, restart_cmfall, kcmf_start, kcmf_end

  Use :: cmf_all_var, Only: ice, icearr

  Implicit None
  Integer (i4b), Parameter :: nd1 = id_ndept

  Real (dp), Allocatable, Dimension (:, :) :: xj_save_old

  Real (dp) :: eblue, ered, qx, qx1

  Integer (i4b) :: ifretest, i, k

  Logical :: newfreq, start

  Data start/.True./

  If (start) Then
    Call prep_convol
    Write (999, *)
    Write (999, *) ' restart_cmfall', restart_cmfall
    Print *
    Print *, ' restart_cmfall', restart_cmfall
    If (restart_cmfall) Then
!     restart
      If (.Not. allocated(xj_save)) Then
        Write (999, *) &
          ' stop: restart and xj_save not allocated in prep_smooth!'
        Stop ' restart and xj_save not allocated in prep_smooth!'
      End If
      Allocate (xj_smoothed(nd1,ifre))
!     make sure that frequency grid has not changed! this should never happen
      ifretest = ubound(xj_save, 2)
      If (ifretest/=ifre) Then
        Write (999, *) ' old:', ifretest, ' new:', ifre
        Write (999, *) ' stop: IFRE HAS CHANGED IN PREP_SMOOTH (1)!!!!'
        Print *, 'old:', ifretest, ' new:', ifre
        Stop ' IFRE HAS CHANGED IN PREP_SMOOTH (1)!!!!'
      End If
!     smoothing not required, since input value has been smoothed before
!     saving
      xj_smoothed = xj_save
      Print *
      Print *, ' restart of CMF with KCMF_START/END = ', kcmf_start, kcmf_end
    Else
      Allocate (xj_smoothed(nd1,ifre), xj_save(nd1,ifre))
!     allocate(xj_smoothed(nd1,ifre)) !for tests with restart-input XJ(obs.
!     frame)
!     smooth XJ for further use in Thomson emissivity (cmf_complete)
!     here: XJ contains only obs. frame quantities (after CONT)
!     smoothed XJ overwrites XJ_SMOOTH -> XJ_CMF_COARSE
      Call xj_smooth(2)
!     start set to .false. after end of freloop
      Write (999, *)
      Write (999, *) ' after xj_smooth'
      Print *
      Print *, ' after xj_smooth'
    End If
  Else
!   make sure that frequency grid has not changed
!   within iteration, frequency grid can change, due to temp-update
    newfreq = .False.
    ifretest = ubound(xj_smoothed, 2)
    If (ifretest/=ifreold1) Then
      Write (999, *) ifretest, ifreold1, ifre
      Write (999, *) ' stop: ifretest ne ifreold1'
      Print *, ifretest, ifreold1, ifre
      Stop ' ifretest ne ifreold1'
    End If
!   JO March 2022
!   when, e.g., Helium has a very low abundance, the main ions of other
!   elements might change during the T-iteration, and thus the freq. grid
!   might change. In this case, xj_save and xj_smooth need to be adapted
    If (ifretest/=ifre) Then
      Print *, ' old:', ifretest, ' new:', ifre
      Print *, ' IFRE HAS CHANGED IN PREP_SMOOTH (2)!!!!'
      Print *
      newfreq = .True.
    Else
      Do k = 1, ifre
        If (abs(1.D0-fre(k)/freold(k))>1.D-15) Then
          Print *, k, ' ', fre(k), ' ', freold(k)
          Print *, ' FREQUENCY GRID HAS CHANGED IN PREP_SMOOTH (2)!!!!'
          Print *
          newfreq = .True.
          Exit
        End If
      End Do
    End If
    If (newfreq) Then
      Allocate (xj_save_old(nd1,ifreold1))
      xj_save_old = xj_save
      Deallocate (xj_smoothed, xj_save)
      Allocate (xj_smoothed(nd1,ifre), xj_save(nd1,ifre))
!     as in optcmf_restart
      xj_save(:, 1) = xj_save_old(:, 1)
      xj_save(:, ifre) = xj_save_old(:, ifretest)
outer: Do i = 2, ifre - 1
        Do k = 1, ifreold1 - 1
          If (fre(i)>freold(k) .And. fre(i)<=freold(k+1)) Then
            qx = log10(fre(i)/freold(k))/log10(freold(k+1)/freold(k))
            qx1 = 1.D0 - qx
            xj_save(:, i) = qx1*log10(xj_save_old(:,k)) + &
              qx*log10(xj_save_old(:,k+1))
            xj_save(:, i) = 10.D0**xj_save(:, i)
            Cycle outer
          End If
        End Do
        Write (999, *) ' stop: fre not found in freold (subr. prep_smooth)'
        Stop ' fre not found in freold (subr. prep_smooth)'
      End Do outer

      ifreold1 = ifre
      freold = fre
      Print *, &
        ' prep_smooth: XJ_SAVE (smoothed) interpolated from previous values!'
      Print *

      Deallocate (xj_save_old)

      eblue = 1.D8/wavblue
      ered = 1.D8/wavred

      Do i = 1, ifre
        If (fre(i)>ered) Exit
      End Do
      kcmf_start = i

      Do i = ifre, 1, -1
        If (fre(i)<eblue) Exit
      End Do
      kcmf_end = i

      If (abs(1.-fre(kcmf_start)/ered)>0.1 .Or. abs(1.-fre( &
        kcmf_end)/eblue)>0.1) Then
        Write (999, *) &
          ' stop: Problems with kcmf_start or kcmf_end in cmf_complete'
        Stop ' Problems with kcmf_start or kcmf_end in cmf_complete'
      End If
      Print *, ' xj_save updated (re-dimensioned) in prep_smooth'
      Print *

      Deallocate (ice, icearr)
      Call prep_convol

    End If

!   smooth XJ_SAVE for further use in Thomson emissivity (cmf_complete)
!   here: XJ_SAVE contains CMF + obs. frame quantities
!   (from previous call of INTERPOL_CMF_TO_COARSE)
!   smoothed XJ_SAVE overwrites XJ_SMOOTH -> XJ_CMF_COARSE
    Call xj_smooth(3)
  End If

  start = .False.

  Return
End Subroutine

!----------------------------------------------------------------------------

Subroutine prep_convol

! prepares weights and corresponding array for approx. convolution with
! electron-scattering profile

! slightly improved by JP Nov 2023: for large broadening, boundary
! conditions insufficient in previous verions

  Use :: nlte_type
  Use :: fund_const, Only: clight, akb, emass
  Use :: nlte_var, Only: ifre, fre
  Use :: nlte_app, Only: teff
  Use :: cmf_all_var, Only: ice, icearr
  Implicit None

  Real (dp) :: vthec, vthecm, e, de, ex, dx, summ

  Integer (i4b) :: k, kk, il, im, nmax

! assuming thermal broadening at Teff (to save space and time);
! treatment approximate anyway)
! broadening considered until 3*delta_nu
  vthec = sqrt(2.D0*akb*teff/emass)/clight
  vthecm = 3.*vthec

  Allocate (ice(2,ifre))

  ice(1, 1) = 1
  ice(2, 1) = 1
  ice(1, ifre) = ifre
  ice(2, ifre) = ifre

  Do k = 2, ifre - 1
    e = fre(k)
    de = e*vthecm
!   JP Nov 2023: next line inserted
    ice(1, k) = 1
    Do kk = k - 1, 1, -1
      If (e-fre(kk)>de) Then
        ice(1, k) = kk
        Exit
      End If
    End Do
!   JP Nov 2023: next line inserted
    ice(2, k) = ifre
    Do kk = k + 1, ifre
      If (fre(kk)-e>de) Then
        ice(2, k) = kk
        Exit
      End If
    End Do
!   print*,k,e,ice(1,k),ice(2,k)
  End Do

  nmax = maxval(ice(2,:)-ice(1,:)) + 1

  Allocate (icearr(nmax,ifre))

! first and last weight unity
  icearr(1, 1) = 1.
  icearr(1, ifre) = 1.

! factor 0.5 and 1/(sqrt(pi*delta nu) left out, since renormalized anyway
! relative integration weights from trapez
  Do k = 2, ifre - 1
!   do k=2,3
    e = fre(k)
    de = e*vthec
    il = ice(1, k)
    im = ice(2, k)

    If (im-il<2) Then
      Write (999, *) ' stop: IM-IL < 2 in PREP_CONVOL'
      Stop ' IM-IL < 2 in PREP_CONVOL'
    End If

    kk = il
    ex = (e-fre(kk))/de
    dx = fre(kk+1) - fre(kk)
    icearr(kk+1-il, k) = exp(-ex*ex)*dx
!   print*,il,ex,dx,icearr(kk+1-il,k)

    kk = im
    ex = (e-fre(kk))/de
    dx = fre(kk) - fre(kk-1)
    icearr(kk+1-il, k) = exp(-ex*ex)*dx
!   print*,im,ex,dx,icearr(kk+1-il,k)

    Do kk = il + 1, im - 1
      ex = (e-fre(kk))/de
      dx = fre(kk+1) - fre(kk-1)
      icearr(kk+1-il, k) = exp(-ex*ex)*dx
!     print*,kk,ex,dx,icearr(kk+1-il,k)
    End Do
!   renormalization
    summ = 0.
!   note change in indices
    im = im - il + 1
    il = 1
    Do kk = il, im
      summ = summ + icearr(kk, k)
    End Do
    summ = 1./summ
    icearr(il:im, k) = icearr(il:im, k)*summ

!   for typical tests
!   if(k.eq.990) then
!   print*,1.d8/fre(k)
!   print*,il,im
!   print*,1.d8/fre(ice(1,k):ice(2,k))
!   print*,icearr(il:im,k)
!   endif
  End Do

End Subroutine

!----------------------------------------------------------------------------

Subroutine xj_smooth(iopt)

! smooth xj for further use in Thomson emissivity
! if iopt=1, then xj (and xj_smoothed) will be overwritten;
! additionally, xxk will be smoothed.
! so far, only to be called in formal (nlte.f90)
! if iopt=2, then xj will be smoothed and xj_smoothed overwritten
! so far, only to be called in cmf_complete (cmf_all.f90)
! if iopt=3, then xj_save will smoothed and xj_smoothed overwritten
! so far, only to be called in cmf_complete (cmf_all.f90)


  Use :: nlte_type
  Use :: nlte_dim, Only: id_ndept
  Use :: nlte_var, Only: ifre, xj, xj_smoothed => xj_cmf_coarse, xxk, xj_save, &
    fre
  Use :: cmf_all_var, Only: ice, icearr
  Implicit None

  Integer (i4b), Intent (In) :: iopt

  Integer (i4b) :: k, kk, il, im
  Real (dp), Dimension (id_ndept) :: xj_aux

  Real (dp) :: eddftest

  If (iopt==1) Then

!   smooth xxk
    Do k = 1, ifre
      il = ice(1, k)
      im = ice(2, k)
      xj_aux = 0.D0
      Do kk = il, im
        xj_aux = xj_aux + xxk(:, kk)*icearr(kk+1-il, k)
      End Do
      xj_smoothed(:, k) = xj_aux !  xj_smoothed used as auxiliary variable
    End Do
!   overwrite xxk
    Do k = 1, ifre
      xxk(:, k) = xj_smoothed(:, k)
    End Do

!   smooth xj
    Do k = 1, ifre
      il = ice(1, k)
      im = ice(2, k)
      xj_aux = 0.D0
      Do kk = il, im
        xj_aux = xj_aux + xj(:, kk)*icearr(kk+1-il, k)
      End Do
      xj_smoothed(:, k) = xj_aux !  xj_smoothed used as auxiliary variable
    End Do

!   overwrite xj
    Do k = 1, ifre
      xj(:, k) = xj_smoothed(:, k)
      eddftest = xxk(id_ndept, k)/xj(id_ndept, k)
!     JO March 2018
      If (1.D8/fre(k)>20. .And. (eddftest<0.31 .Or. eddftest>0.35)) Then
        Write (999, *) ' Warning!!! Warning!!! Warning!!!'
        Write (999, *) 1.D8/fre(k), ' ', eddftest
        Write (999, *) ' inconsistent Edd. factor in xj_smooth'
        Write (999, *)
        Print *, ' Warning!!! Warning!!! Warning!!!'
        Print *, 1.D8/fre(k), ' ', eddftest
        Print *, ' inconsistent Edd. factor in xj_smooth'
        Print *
      End If
    End Do

  Else If (iopt==2) Then

!   smooth only xj, and do not overwrite
    Do k = 1, ifre
      il = ice(1, k)
      im = ice(2, k)
      xj_aux = 0.D0
      Do kk = il, im
        xj_aux = xj_aux + xj(:, kk)*icearr(kk+1-il, k)
      End Do
      xj_smoothed(:, k) = xj_aux !  overwrite xj_smoothed
    End Do

  Else If (iopt==3) Then

!   smooth only xj_save, and do not overwrite
    Do k = 1, ifre
      il = ice(1, k)
      im = ice(2, k)
      xj_aux = 0.D0
      Do kk = il, im
        xj_aux = xj_aux + xj_save(:, kk)*icearr(kk+1-il, k)
      End Do
      xj_smoothed(:, k) = xj_aux !  overwrite xj_smoothed
    End Do

  Else
    Print *, iopt
    Write (999, *) ' stop: wrong iopt in xj_smooth'
    Stop ' wrong iopt in xj_smooth'
  End If

  Write (999, *)
  Write (999, *) ' subr. xj_smooth(', iopt, ') called'
  Print *
  Print *, ' subr. xj_smooth(', iopt, ') called'

  Return
End Subroutine

!----------------------------------------------------------------------------

Subroutine save_line_list
  Use :: nlte_type
  Use :: nlte_var, Only: modnam, lwion1
  Use :: nlte_app, Only: lwion, indrec
  Use :: cmf_all_var, Only: ntot, id1, xlam1, gf1

  Open (1, File=trim(modnam)//'/LINE_LIST_MERGED', Form='UNFORMATTED')
  Write (1) ntot, lwion, lwion1, indrec
  Write (1) id1(1:ntot), xlam1(1:ntot), gf1(1:ntot)
  Close (1)

! print*,ntot
! print*,xlam1(1),xlam1(ntot)

  Return
End Subroutine

!-----------------------------------------------------------------------

Subroutine chih_params1(chi, p, t, delchi, sigemin, const, expo, expo_p, nd, &
  ns, nsig)

! parameterizes chibar_h, in case by optimizing sigemin
! changed April 2023 to include multi-linear fit as in chih_params

  Use :: nlte_type
  Use :: nlte_dim, Only: id_ndept
  Implicit None

  Integer (i4b), Parameter :: nd1 = id_ndept
  Integer (i4b), Intent (In) :: nd, ns, nsig
  Real (dp), Dimension (nd), Intent (In) :: chi, p, t
  Real (dp), Dimension (nd), Intent (Out) :: delchi

  Real (dp), Intent (Inout) :: sigemin ! different from original version
  Real (dp), Intent (Out) :: const, expo, expo_p

  Real (dp), Dimension (nd1) :: fit, sig, chi1, t1, p1, sig1, fit_lin, &
    delchi_lin

  Integer (i4b) :: imin, i, l, minimin, minimin_lin, ntot
  Real (dp) :: chi2min, expoi, consti, corr, chi2, maxchi, sigemini, chisq, &
    sigemin_orig
  Real (dp) :: expo_lin, expo_p_lin, const_lin, sigemin_lin

  Integer (i4b), Dimension (3) :: ia
  Real (dp), Dimension (3, 3) :: covar
  Real (dp), Dimension (3) :: a
  External :: lfit_funcs

! at first, old method (with updated calculation of CHI2), i.e, with given
! SIGEMIN
  chi2min = 1.D5
  expo = 0.
  minimin = 0
! test for best interval to perform the fit

  Do imin = nsig, nd - 4
!   regression from imin on
    Call linreg(log10(t(imin:nd)), log10((chi(imin:nd)- &
      sigemin)/p(imin:nd)), nd+1-imin, expoi, consti, corr)

!   fit quality for "all" points (from nsig on)
    fit(nsig:nd) = (10.**consti*t(nsig:nd)**expoi*p(nsig:nd)) + sigemin
    chi2 = sum((fit(nsig:nd)/chi(nsig:nd)-1.D0)**2)
    If (chi2<chi2min) Then
      minimin = imin
      chi2min = chi2
      expo = expoi
      expo_p = 1.D0 !               JO: new parameter
      const = consti
    End If
  End Do

  If (expo==0.) Then
    Write (999, *) ' stop: regression not successful (chih_params1#1)'
    Stop ' regression not successful (chih_params1#1)'
  End If

  const = 10.**const
! for all photospheric points (from NS on)
  fit(ns:nd) = (const*t(ns:nd)**expo*p(ns:nd)) + sigemin
  delchi(ns:nd) = chi(ns:nd)/fit(ns:nd)
  Write (999, *) ' linear method 1'
  Write (999, *) ' deviations in range ', minval(delchi(ns:nd)), ' ', &
    maxval(delchi(ns:nd))
  Write (999, *)
  Print *, ' linear method 1'
  Print *, ' deviations in range ', minval(delchi(ns:nd)), ' ', &
    maxval(delchi(ns:nd))
  Print *
! do i=ns,nd
! write(*,10) i,chi(i),fit(i),delchi(i)
! enddo

! new condition
! if(maxval(delchi(ns:nd)).le.1.2 .and. minval(delchi(ns:nd)).ge.0.8) then
  If (maxval(delchi(ns:nd))<=1.05 .And. minval(delchi(ns:nd))>=0.95) Then
    Write (999, *)
    Write (999, *) ' routine CHIH_PARAMS1: linear method 1 used'
    Write (999, *) ' ns = ', ns, ' nsig = ', nsig, ' imin = ', minimin
    Write (999, *) ' sigemin_orig: ', sigemin
    Write (999, *) ' parameters: ', sigemin, ' ', const, ' ', expo, ' ', &
      expo_p
    Write (999, *) ' deviations in range ', minval(delchi(ns:nd)), &
      maxval(delchi(ns:nd))
    Write (999, *)
    Print *
    Print *, ' routine CHIH_PARAMS1: linear method 1 used'
    Print *, ' ns = ', ns, ' nsig = ', nsig, ' imin = ', minimin
    Print *, ' sigemin_orig: ', sigemin
    Print *, ' parameters: ', sigemin, ' ', const, ' ', expo, ' ', expo_p
    Print *, ' deviations in range ', minval(delchi(ns:nd)), &
      maxval(delchi(ns:nd))
    Print *
    Do i = ns, nd
      Write (999, 110) i, chi(i), fit(i), t(i), p(i)
      Write (*, 110) i, chi(i), fit(i), t(i), p(i)
    End Do
    Print *
    Return
  End If

! if DELCHI problematic (e.g., for first models of d10v), optimize sigemin

  chi2min = 1.D5
  expo = 0.
  minimin = 0
  sigemin_orig = sigemin
  sigemini = sigemin
  sigemini = int(sigemini*10.)/10.D0 ! round SIGEMIN
  maxchi = maxval(chi(nsig:nd))

  Do While (sigemini<maxchi)

loop: Do imin = nsig, nd - 4
!     regression from imin on
      i = 0
      Do l = imin, nd
        If ((chi(l)-sigemini)>0.) Then
          i = i + 1
          chi1(i) = chi(l)
          t1(i) = t(l)
          p1(i) = p(l)
        End If
      End Do
      If (i<2) Cycle loop

      Call linreg(log10(t1(1:i)), log10((chi1(1:i)-sigemini)/p1(1:i)), i, &
        expoi, consti, corr)

!     fit quality for "all" points (from nsig on)
      fit(nsig:nd) = (10.**consti*t(nsig:nd)**expoi*p(nsig:nd)) + sigemini
!     though we finally calculate chi/fit, we optimize w.r.t. fit/chi, to
!     allow
!     for larger errors if chi is small.
      chi2 = sum((fit(nsig:nd)/chi(nsig:nd)-1.D0)**2)
!     print*,sigemini,imin,chi2
      If (chi2<chi2min) Then
        minimin = imin
        chi2min = chi2
        expo = expoi
        expo_p = 1.D0 !             J0: new parameter
        const = consti
        sigemin = sigemini
      End If
    End Do loop

    sigemini = sigemini + 0.1D0

  End Do

  If (expo==0.) Then
    Write (999, *) ' stop: regression not successful (chih_params1#2)'
    Stop ' regression not successful (chih_params1#2)'
  End If

! print*,chi2min,expo,const,sigemin
  const = 10.**const
! for all photospheric points (from NS on)
  fit(ns:nd) = (const*t(ns:nd)**expo*p(ns:nd)) + sigemin
  delchi(ns:nd) = chi(ns:nd)/fit(ns:nd)
  Write (999, *) ' linear method2'
  Write (999, *) ' deviations in range ', minval(delchi(ns:nd)), ' ', &
    maxval(delchi(ns:nd))
  Write (999, *)
  Print *, ' linear method2'
  Print *, ' deviations in range ', minval(delchi(ns:nd)), ' ', &
    maxval(delchi(ns:nd))
  Print *
! do i=ns,nd
! write(*,10) i,chi(i),fit(i),delchi(i)
! enddo

! new condition
! if(maxval(delchi(ns:nd)).le.1.2 .and. minval(delchi(ns:nd)).ge.0.8) then
  If (maxval(delchi(ns:nd))<=1.05 .And. minval(delchi(ns:nd))>=0.95) Then
    Write (999, *)
    Write (999, *) ' routine CHIH_PARAMS1: linear method 2 used'
    Write (999, *) ' sigemin_orig: ', sigemin_orig
    Write (999, *) ' ns = ', ns, ' nsig = ', nsig, ' imin = ', minimin
    Write (999, *) ' parameters: ', sigemin, ' ', const, ' ', expo, ' ', &
      expo_p
    Write (999, *) ' deviations in range ', minval(delchi(ns:nd)), ' ', &
      maxval(delchi(ns:nd))
    Write (999, *)
    Print *
    Print *, ' routine CHIH_PARAMS1: linear method 2 used'
    Print *, ' sigemin_orig: ', sigemin_orig
    Print *, ' ns = ', ns, ' nsig = ', nsig, ' imin = ', minimin
    Print *, ' parameters: ', sigemin, ' ', const, ' ', expo, ' ', expo_p
    Print *, ' deviations in range ', minval(delchi(ns:nd)), ' ', &
      maxval(delchi(ns:nd))
    Print *
    Do i = ns, nd
      Write (999, 110) i, chi(i), fit(i), t(i), p(i)
      Write (*, 110) i, chi(i), fit(i), t(i), p(i)
    End Do
    Write (999, *)
    Print *
    Return
  Else
    Write (999, *) ' linear fit in chih_params1 not good enough!'
    Write (999, *) ' sigemin_orig: ', sigemin_orig
    Write (999, *) ' parameters: ', sigemin, ' ', const, ' ', expo, ' ', &
      expo_p
    Write (999, *) ' trying multi-linear regression now'
    Write (999, *)
    Print *, ' linear fit in chih_params1 not good enough!'
    Print *, ' sigemin_orig: ', sigemin_orig
    Print *, ' parameters: ', sigemin, ' ', const, ' ', expo, ' ', expo_p
    Print *, ' trying multi-linear regression now'
    Print *
!   save values from linear regression, might be used below in case
    minimin_lin = minimin
    const_lin = const
    expo_lin = expo
    expo_p_lin = expo_p
    fit_lin = fit
    delchi_lin = delchi
    sigemin_lin = sigemin
  End If

! multi-linear regression including P (see subr. LFIT_FUNCS)
! test for best interval to perform the fit, and proceed as for linear fit
  sig = 1. !                        no error
  ia = 1 !                          all parameters

! NOTE: sigemin adapted above
  chi2min = 1.D5
  expo = 0.
  minimin = 0
  sigemini = sigemin_orig
  sigemini = int(sigemini*10.)/10.D0 ! round SIGEMIN
! maxchi already defined

  Do While (sigemini<maxchi)

loop_multi: Do imin = nsig, nd - 4
!     regression from imin on
      i = 0
      Do l = imin, nd
        If ((chi(l)-sigemini)>0.) Then
          i = i + 1
          chi1(i) = chi(l)
          t1(i) = t(l)
          p1(i) = p(l)
          sig1(i) = sig(l)
        End If
      End Do
      If (i<2) Cycle loop_multi
      ntot = i

      Call lfit(log10(t1(1:i)), log10(chi1(1:i)-sigemini), sig1(1:i), p1(1:i), &
        ntot, a, ia, 3, covar, 3, chisq, lfit_funcs)
!     fit parameters contained in a

!     fit quality for "all" points (from nsig on)
      fit(nsig:nd) = (10.**a(1)*t(nsig:nd)**a(2)*p(nsig:nd)**a(3)) + sigemin
!     in case, regarding log-values
!     chi2=sum((log10(chi(nsig:nd))-log10(fit(nsig:nd)))**2/log10(chi(nsig:nd))**2)
!     chi2=sum((chi(nsig:nd)-fit(nsig:nd))**2/chi(nsig:nd)**2)=
      chi2 = sum((fit(nsig:nd)/chi(nsig:nd)-1.D0)**2)
!     print*,' parameters: ',imin,' ',sigemin, chi2
!     print*,' parameters: ',a(1),' ',a(2),' ',a(3)
!     print*

      If (chi2<chi2min) Then
        minimin = imin
        chi2min = chi2
        expo = a(2)
        expo_p = a(3)
        const = a(1)
        sigemin = sigemini
      End If
    End Do loop_multi

    sigemini = sigemini + 0.1D0

  End Do

  If (expo==0.) Then
    Write (999, *) &
      ' stop: MULTI-LINEAR REGRESSION NOT SUCCESSFUL (CHIH_PARAMS)'
    Stop ' MULTI-LINEAR REGRESSION NOT SUCCESSFUL (CHIH_PARAMS)'
  End If

  const = 10.**const
! for all photospheric points (from ns on)
  fit(ns:nd) = (const*t(ns:nd)**expo*p(ns:nd)**expo_p) + sigemin
  delchi(ns:nd) = chi(ns:nd)/fit(ns:nd)
  Write (999, *) ' multi-linear'
  Write (999, *) ' deviations in range ', minval(delchi(ns:nd)), ' ', &
    maxval(delchi(ns:nd))
  Write (999, *)
  Print *, ' multi-linear'
  Print *, ' deviations in range ', minval(delchi(ns:nd)), ' ', &
    maxval(delchi(ns:nd))
  Print *
! do i=ns,nd
! write(*,10) i,chi(i),fit(i),delchi(i)
! enddo

  If (maxval(delchi(ns:nd))>2. .Or. minval(delchi(ns:nd))<0.2) Then
!   still old condition
    Do i = ns, nd
      Write (999, 100) i, chi(i), fit(i), delchi(i)
      Write (*, 100) i, chi(i), fit(i), delchi(i)
    End Do
    Write (999, *) ' no multi-linear regression possible (chih_params1)!'
    Print *, ' no multi-linear regression possible (chih_params1)!'
  Else If (expo>0. .Or. expo_p<0.) Then
!   unphysical parameters, may appear if locally close to eddington limit'
    Write (999, *) ' multi-linear regression parameters not appropriate!'
    Write (999, *) ' expo for t (should be negative) = ', expo
    Write (999, *) ' expo for p (should be in between 0...1) = ', expo_p
    Print *, ' multi-linear regression parameters not appropriate!'
    Print *, ' expo for t (should be negative) = ', expo
    Print *, ' expo for p (should be in between 0...1) = ', expo_p
  Else
    Write (999, *) ' multi-linear fit in chih_params1 successful'
    Print *, ' multi-linear fit in chih_params1 successful'
    If (expo_p>1.) Then
      Write (999, *) ' expo for p = ', expo_p, ' > 1, proceed at own risk!!!'
      Print *, ' expo for p = ', expo_p, ' > 1, proceed at own risk!!!'
    End If
    Write (999, *) ' ns = ', ns, ' nsig = ', nsig, ' imin = ', minimin
    Write (999, *) ' sigemin_orig: ', sigemin_orig
    Write (999, *) ' parameters: ', sigemin, ' ', const, ' ', expo, ' ', &
      expo_p
    Write (999, *) ' deviations in range ', minval(delchi(ns:nd)), ' ', &
      maxval(delchi(ns:nd))
    Write (999, *)
    Print *, ' ns = ', ns, ' nsig = ', nsig, ' imin = ', minimin
    Print *, ' sigemin_orig: ', sigemin_orig
    Print *, ' parameters: ', sigemin, ' ', const, ' ', expo, ' ', expo_p
    Print *, ' deviations in range ', minval(delchi(ns:nd)), ' ', &
      maxval(delchi(ns:nd))
    Print *

!   compare with results from linear fit (method2)
    If (maxval(delchi_lin(ns:nd))<=maxval(delchi(ns:nd)) .And. minval( &
      delchi_lin(ns:nd))>=minval(delchi(ns:nd))) Then
      Write (999, *) ' linear method 2 better suited than multi-linear'
      Print *, ' linear method 2 better suited than multi-linear'
      minimin = minimin_lin
      const = const_lin
      expo = expo_lin
      expo_p = expo_p_lin
      fit = fit_lin
      delchi = delchi_lin
      sigemin = sigemin_lin
      Write (999, *) ' ns = ', ns, ' nsig = ', nsig, ' imin = ', minimin
      Write (999, *) ' sigemin_orig: ', sigemin_orig
      Write (999, *) ' parameters: ', sigemin, ' ', const, ' ', expo, ' ', &
        expo_p
      Write (999, *) ' deviations in range ', minval(delchi(ns:nd)), ' ', &
        maxval(delchi(ns:nd))
      Print *, ' ns = ', ns, ' nsig = ', nsig, ' imin = ', minimin
      Print *, ' sigemin_orig: ', sigemin_orig
      Print *, ' parameters: ', sigemin, ' ', const, ' ', expo, ' ', expo_p
      Print *, ' deviations in range ', minval(delchi(ns:nd)), ' ', &
        maxval(delchi(ns:nd))
    End If
    Write (999, *) ' finally adopted values (chi,fit,delchi)'
    Print *, ' finally adopted values (chi,fit,delchi)'
    Do i = ns, nd
      Write (999, 100) i, chi(i), fit(i), delchi(i)
      Write (*, 100) i, chi(i), fit(i), delchi(i)
    End Do
    Write (999, *)
    Write (999, *) ' finally adopted values (chi,fit,t,p)'
    Print *
    Print *, ' finally adopted values (chi,fit,t,p)'
    Do i = ns, nd
      Write (999, 110) i, chi(i), fit(i), t(i), p(i)
      Write (*, 110) i, chi(i), fit(i), t(i), p(i)
    End Do
    Write (999, *)
    Print *
    Return
  End If

! final possibility: if linear fit not too bad, use these values
  If (maxval(delchi_lin(ns:nd))<=2 .And. minval(delchi_lin(ns:nd))>=0.2) Then
    Write (999, *) ' linear regression used instead'
    Print *, ' linear regression used instead'
    minimin = minimin_lin
    const = const_lin
    expo = expo_lin
    expo_p = expo_p_lin
    fit = fit_lin
    delchi = delchi_lin
    sigemin = sigemin_lin
    Write (999, *) ' ns = ', ns, ' nsig = ', nsig, ' imin = ', minimin
    Write (999, *) ' sigemin_orig: ', sigemin_orig
    Write (999, *) ' parameters: ', sigemin, ' ', const, ' ', expo, ' ', &
      expo_p
    Write (999, *) ' deviations in range ', minval(delchi(ns:nd)), ' ', &
      maxval(delchi(ns:nd))
    Write (999, *) ' finally adopted values (chi,fit,delchi)'
    Print *, ' ns = ', ns, ' nsig = ', nsig, ' imin = ', minimin
    Print *, ' sigemin_orig: ', sigemin_orig
    Print *, ' parameters: ', sigemin, ' ', const, ' ', expo, ' ', expo_p
    Print *, ' deviations in range ', minval(delchi(ns:nd)), ' ', &
      maxval(delchi(ns:nd))
    Print *, ' finally adopted values (chi,fit,delchi)'
    Do i = ns, nd
      Write (999, 100) i, chi(i), fit(i), delchi(i)
      Write (*, 100) i, chi(i), fit(i), delchi(i)
    End Do
    Write (999, *)
    Write (999, *) ' finally adopted values (chi,fit,t,p)'
    Print *
    Print *, ' finally adopted values (chi,fit,t,p)'
    Do i = ns, nd
      Write (999, 110) i, chi(i), fit(i), t(i), p(i)
      Write (*, 110) i, chi(i), fit(i), t(i), p(i)
    End Do
    Print *
    Write (999, *)
    Return
  Else
!   problems cannot be cured
    Write (999, *) ' stop: NO REGRESSION POSSIBLE (CHIH_PARAMS1)! &
      &RE-TRY WITH UPDATE_STRUCT = .FALSE.'
    Stop ' NO REGRESSION POSSIBLE (CHIH_PARAMS1)! RE-TRY WITH UPDATE_STRUCT &
      &= .FALSE.'
  End If

100 Format (I3, 2X, 3(E10.3,2X))
110 Format (I3, 2X, 4(E10.3,2X))

End Subroutine

!-----------------------------------------------------------------------

Subroutine overlap_detector(dvdr)

! detects strong overlapping lines for lines from explicit and selected
! elements
! note: uses reduced line-list

  Use :: nlte_type
  Use :: fund_const, Only: clight
  Use :: fastwind_params, Only: gf_cut, opal_cut, tau_cut, fac_vturb
  Use :: princesa_var, Only: zl, labl
  Use :: nlte_app, Only: jatom_full
  Use :: cmf_all_var, Only: id1, xlam1, opal_ntest, vturb_cmf, icfirst, &
    iclast, indexel_inv, index_id1, gf1

  Implicit None

  Real (dp), Intent (In) :: dvdr

  Integer (i4b) :: ii, iexp, isel, isel1, index, kk, irest, k, j, low, up, i, &
    indexi, ki, ji, lowi, upi, nsize, nsize1, in, kold, jold, lowold

  Real (sp) :: lam, opalii, lami

  Real (dp) :: fac, red, blue, dxlam, lamm, lamp, tauii

  Logical :: expl, new, expli

  nsize = size(id1)
  nsize1 = size(index_id1)
  If (nsize/=nsize1) Then
    Write (999, *) ' stop: error in nsize (overlap_detector)'
    Stop ' error in nsize (overlap_detector)'
  End If
! sort according to index;
! after sorting, selected elements first, explicit elements at the end
  Call indexx_int(nsize, id1, index_id1)

  fac = fac_vturb*vturb_cmf/clight
  blue = xlam1(icfirst)
  red = xlam1(iclast)
  Print *
  Print *, ' overlap-detector between', blue, red

  iexp = 0
  isel = 0
  isel1 = 0

  kold = 0
  jold = 0
  lowold = 0

iiloop: Do in = 1, nsize
!   sorted according to element, ionization, lower and upper level
    ii = index_id1(in)
    If (gf1(ii)<gf_cut) Cycle iiloop
    index = id1(ii)
    If (index==0) Cycle iiloop
    lam = xlam1(ii)
    If (lam<blue .Or. lam>red) Cycle iiloop
    dxlam = fac*lam
    lamm = lam - dxlam
    lamp = lam + dxlam
    If (index>=1D8) Then
!     explicit element
      kk = index/1D8
      irest = index - kk*1D8
      low = irest/10000
      up = irest - low*10000
      k = indexel_inv(kk)
      j = int(zl(low)+1.D-6) + 1
      opalii = opal_ntest(ii)
      tauii = opalii/dvdr
      If (tauii<tau_cut) Cycle iiloop
      If (opalii==0.) Then
        Write (999, *) ' stop: opalii = 0 for explicit element'
        Stop ' opalii = 0 for explicit element'
      End If
      expl = .True.
      iexp = iexp + 1
    Else
      k = index/1000000
      If (jatom_full(k)/=1) Cycle iiloop
      irest = index - k*1000000
      j = irest/100000
      irest = irest - j*100000
      low = irest/100
      up = irest - low*100
      If (up==0) Cycle iiloop
      opalii = opal_ntest(ii)
      tauii = opalii/dvdr
      If (tauii<tau_cut) Cycle iiloop
      isel1 = isel1 + 1
      If (opalii==0.) Cycle
      isel = isel + 1
      expl = .False.
    End If

    new = .True.

iloop1: Do i = ii - 1, icfirst, -1
      lami = xlam1(i)
      If (lami<lamm) Exit iloop1
      If (opal_ntest(i)<opal_cut*opalii) Cycle iloop1
      indexi = id1(i)
      If (indexi>=1D8) Then
!       explicit element
        kk = indexi/1D8
        irest = indexi - kk*1D8
        lowi = irest/10000
        upi = irest - lowi*10000
        ki = indexel_inv(kk)
        ji = int(zl(lowi)+1.D-6) + 1
        expli = .True.
      Else
        ki = indexi/1000000
        irest = indexi - ki*1000000
        ji = irest/100000
        irest = irest - ji*100000
        lowi = irest/100
        upi = irest - lowi*100
        expli = .False.
      End If

      If (k/=kold .Or. j/=jold .Or. low/=lowold) Then
        Print *
        kold = k
        jold = j
        lowold = low
      End If

      If (new) Then
        If (expl) Then
          Write (*, 100) k, j, low, up, lam, adjustr(labl(low)), &
            adjustr(labl(up))
        Else
          Write (*, 110) k, j, low, up, lam
        End If
        new = .False.
      End If
      If (expli) Then
        Write (*, 120) ki, ji, lowi, upi, lami, opal_ntest(i)/opalii, &
          adjustr(labl(lowi)), adjustr(labl(upi))
      Else
        Write (*, 130) ki, ji, lowi, upi, lami, opal_ntest(i)/opalii
      End If

    End Do iloop1

iloop2: Do i = ii + 1, iclast
      lami = xlam1(i)
      If (lami>lamp) Exit iloop2
      If (opal_ntest(i)<opal_cut*opalii) Cycle iloop2
      indexi = id1(i)
      If (indexi>=1D8) Then
!       explicit element
        kk = indexi/1D8
        irest = indexi - kk*1D8
        lowi = irest/10000
        upi = irest - lowi*10000
        ki = indexel_inv(kk)
        ji = int(zl(lowi)+1.D-6) + 1
        expli = .True.
      Else
        ki = indexi/1000000
        irest = indexi - ki*1000000
        ji = irest/100000
        irest = irest - ji*100000
        lowi = irest/100
        upi = irest - lowi*100
        expli = .False.
      End If

      If (k/=kold .Or. j/=jold .Or. low/=lowold) Then
        Print *
        kold = k
        jold = j
        lowold = low
      End If

      If (new) Then
        If (expl) Then
          Write (*, 100) k, j, low, up, lam, adjustr(labl(low)), &
            adjustr(labl(up))
        Else
          Write (*, 110) k, j, low, up, lam
        End If
        new = .False.
      End If
      If (expli) Then
        Write (*, 120) ki, ji, lowi, upi, lami, opal_ntest(i)/opalii, &
          adjustr(labl(lowi)), adjustr(labl(upi))
      Else
        Write (*, 130) ki, ji, lowi, upi, lami, opal_ntest(i)/opalii
      End If

    End Do iloop2

  End Do iiloop

  Print *
  Print *, ' overlap interval +/- ', fac_vturb*vturb_cmf/1.D5, ' km/s'
  Print *, ' all lines until gf > ', gf_cut
  Print *, ' all lines until tau_Sob > ', tau_cut
  Print *, ' all overlaps until a ratio > ', opal_cut, &
    ' (overlapping to considered line)'
  Print *
  Print *, ' number of lines from explicit elements ', iexp
  Print *, ' number of lines from selected elements (up ne 0) ', isel
  Print *, ' number of lines from selected elements (up ne 0, opal ne 0) ', &
    isel1
  Print *
  Return

100 Format (I2, 2X, I2, 2X, I5, 2X, I5, 2X, F12.4, 2X, A6, 2X, A6)
110 Format (I2, 2X, I2, 2X, I5, 2X, I5, 2X, F12.4)

120 Format ('overlap with: ', I2, 2X, I2, 2X, I5, 2X, I5, 2X, F12.4, 2X, &
    E10.4, 2X, A6, 2X, A6)
130 Format ('overlap with: ', I2, 2X, I2, 2X, I5, 2X, I5, 2X, F12.4, 2X, &
    E10.4)


End Subroutine

!-----------------------------------------------------------------------

Subroutine indexx_int(n, arr, indx)

! as subr. indexx (nlte.f90), but for integer array

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
  Integer (i4b) :: arr(n), indx(n)
! ..
! .. local scalars ..
! JO changed March 2017
! real(dp) ::  a
  Integer (i4b) :: a
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
    If (jstack>nstack) Then
      Write (999, *) ' stop: nstack too small in indexx_int'
      Stop ' nstack too small in indexx_int'
    End If
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

Subroutine overlap_noup(clfac)

! detects lines with noup overlapping with lines from explicit
! (non-HHe) elements, and prepares special treatment for corresponding
! source-functions in subr. opacity_cmf
! Note(1): only if the opacity of the explicit-element line is dominating,
! special treatment will be applied. In this way, also the impact of
! two or more narrowly spaced lines will be  treated correctly:
! a) in case the 2nd line is from background-elements, this line will either
! have a marginal influence on the (strong) explicit line, and the special
! treatment takes place. Or 2nd line is dominating, and nothing needs to be
! done, since the first line will be not dominating.
! b) in case the 2nd line is from an explicit element, and this line is
! strong, its impact is correctly accounted for, since any previous effect
! from the first line will be overwritten for those lines where the 2nd one
! affects them stronger than the first line.

! Note(2): uses reduced line-list

  Use :: nlte_type
  Use :: nlte_dim, Only: id_ndept
  Use :: fund_const, Only: clight, c_tau

  Use :: fastwind_params, Only: xm, minrho

  Use :: princesa_var, Only: gl
  Use :: nlte_var, Only: sr, vmax, lwion1
  Use :: cmf_all_var, Only: id1, xlam1, gf1, icfirst, iclast, nftot, &
    index_lam_cmf, lam_cmf, indexel_inv, nf_cmf, occexpl, vdoptot, vref, &
    profdop, sumopal_cmf, id1_overlap, no2

  Implicit None

  Integer (i4b), Parameter :: nd1 = id_ndept

  Real (dp), Dimension (nd1), Intent (In) :: clfac

  Real (dp), Dimension (nd1) :: opal, occlow, occup

  Integer (i4b) :: lastindex, i, in, ii, index, kk, k, irest, ml, mu, ixmax, &
    ix, l, ik, index1, up, no1

  Real (dp) :: lamgrid, lam, srvmax, c_tau1, const, rho1, mirho, lam1, deltav, &
    dw1, dw2

! id_overlap allocated in subr. line-list
! id_overlap initialized (=0) in subr. opacity_cmf for ipath=1

  srvmax = sr/vmax
  c_tau1 = c_tau*srvmax

  lastindex = icfirst - 1

  Do i = 1, nftot
    in = index_lam_cmf(i)
    lamgrid = lam_cmf(i)

    If (in/=0) Then
overlap: Do ii = lastindex + 1, in
        index = id1(ii)
        If (index<1D8) Cycle overlap ! not explicit element
        kk = index/1D8
        k = indexel_inv(kk)
        If (k==1 .Or. k==2) Cycle overlap ! H or He
        irest = index - kk*1D8
        ml = irest/10000
        mu = irest - ml*10000
!       in units of vref = vturb/3 (with vturb=max(vturb,vturbmin)
        ixmax = nf_cmf(k)
        lam = xlam1(ii)
!       print*
!       print*,ml,mu,lam
        const = c_tau1*lam*gf1(ii)
        occlow = occexpl(ml, :)/gl(ml)
        occup = occexpl(mu, :)/gl(mu)
!       pi e2/me c * 1.d-8 *SR/vmax * lam * gf * (nl/gl - nu/gu) /clfac
        opal = const*(occlow-occup)/clfac
        If (opal(1)==0.) Then
          Write (999, *) &
            ' stop: opal = 0 for expl. element, subr. overlap_noup)'
          Stop ' opal = 0 for expl. element, subr. overlap_noup)'
        End If
!       below LTE for bg-elements; note: lwion1 (selected) > lwion
!       (non-selected)
        Do l = 1, lwion1 - 1
          ix = xm*vdoptot(l, k)/vref + 1 ! in units of vref/3
          If (ix>ixmax) Then
            Write (999, *) ' stop: ix > ixmax in subr. overlap_nopup'
            Stop ' ix > ixmax in subr. overlap_nopup'
          End If
          mirho = 1000.
          Do ik = 0, ix
            rho1 = opal(l)*profdop(ik, l, k)/sumopal_cmf(l, i+ik)
!           significant or dominating inversion of bg opacities
            If (rho1>1.1D0 .Or. rho1<0.D0) rho1 = 0.D0
            mirho = min(mirho, rho1)
          End Do
          Do ik = 1, ix
            rho1 = opal(l)*profdop(ik, l, k)/sumopal_cmf(l, i-ik)
!           significant or dominating inversion of bg opacities
            If (rho1>1.1D0 .Or. rho1<0.D0) rho1 = 0.D0
            mirho = min(mirho, rho1)
          End Do
          If (mirho>minrho) Then
!           print*,ml,mu,l,lam,mirho
            Do ik = ii + 1, iclast
              index1 = id1(ik)
              If (index1>=1D8) Cycle
              irest = index1/100
              up = index1 - irest*100
              If (up/=0) Cycle
              lam1 = xlam1(ik)
              deltav = (lam1/lam-1.D0)*clight
              If (deltav>xm*vdoptot(l,k)) Exit
!             print*,'overlap with'
!             print*,lam1,deltav/1.d5,index1
              If (id1_overlap(ik)==0) Then
                id1_overlap(ik) = ii
              Else
                If (id1_overlap(ik)/=ii) Then
!                 in this case, line ik is closer to new line ii (since line
!                 ik is longwards
!                 of new line ii  (i.e., even more longwards from old line ii)
                  id1_overlap(ik) = ii
!                 print*,'2nd line(+)',l,id1(id1_overlap(ik)),id1(ii),index1
!                 print*,'2nd line(+)',xlam1(id1_overlap(ik)),xlam1(ii),lam1
                End If
              End If

            End Do

            Do ik = ii - 1, icfirst, -1
              index1 = id1(ik)
              If (index1>=1D8) Cycle
              irest = index1/100
              up = index1 - irest*100
              If (up/=0) Cycle
              lam1 = xlam1(ik)
              deltav = (1.D0-lam1/lam)*clight
              If (deltav>xm*vdoptot(l,k)) Exit
!             print*,'overlap with'
!             print*,lam1,deltav/1.d5,index1
              If (id1_overlap(ik)==0) Then
                id1_overlap(ik) = ii
              Else
                If (id1_overlap(ik)/=ii) Then
!                 here, line ik is in between old and new line ii. For reasons
!                 of simplicity
!                 we choose the closer one as being of major impact. This is
!                 justified, since
!                 both lines have a signficant impact (min(rho1) > minrho

!                 ik can be even shortwards of old ii
                  dw1 = abs(lam1-xlam1(id1_overlap(ik)))
!                 ik must be shortwards of new ii
                  dw2 = lam - lam1
                  If (dw2<0.) Then
                    Write (999, *) ' stop: dw2 < 0 in subr. overlap_noup'
                    Stop ' dw2 < 0 in subr. overlap_noup'
                  End If
!                 line ik closer to new line ii
                  If (dw2<dw1) id1_overlap(ik) = ii
!                 print*,'2nd line(-)',l,id1(id1_overlap(ik)),id1(ii),index1
!                 print*,'2nd line(-)',xlam1(id1_overlap(ik)),xlam1(ii),lam1
!                 JO keep this stop statement; test for more explicit elements
!                 if stop occurs, design a 2nd array id2_overlap and include
!                 the 2nd line
!                 in opacity_cmf, then test which line affects what
!                 stop
                End If
              End If
            End Do
          End If !                  mirho
        End Do !                    depth-loop
      End Do overlap !              all overlapping lines at cmf freq. lamgrid
      lastindex = in
    End If !                        if lines present at lamgrid
  End Do !                          all cmf-frequencies

  no1 = 0
  no2 = 0
  Do ii = icfirst, iclast
    index = id1(ii)
    If (index>=1D8) Cycle
    irest = index/100
    up = index - irest*100
    If (up/=0) Cycle
    no1 = no1 + 1
    If (id1_overlap(ii)/=0) no2 = no2 + 1
  End Do

  Print *
  Print *, no1, ' lines from background elements with no upper level'
  Print *, no2, ' from those with modified source-function'

  Return
End Subroutine

!-----------------------------------------------------------------------

Subroutine transic_count(lineno, retchar)

  Use :: nlte_type
  Use :: nlte_dim
  Use :: ffr_error
  Implicit None

! provides number of lines in file LINES_xxx.dat

! ..
  Integer (i4b) :: lineno
  Character :: retchar*4, ret*4
! ..
! .. local scalars ..
  Real (dp) :: realo
  Integer (i4b) :: i, integ, nc, ndatos
  Character :: key*6

! ..
! .. external subroutines ..
  External :: ffrkey, ffrnum
! ..

  lineno = 0

100 Continue
  Call ffrkey(key, nc, ret)
  If (on_error(ret)) Go To 120
  If (key/='CL') Go To 130

  Call ffrkey(key, nc, ret)
  If (on_error(ret)) Go To 120
  If (key/='TY') Go To 130

  Call ffrkey(key, nc, ret)
  If (on_error(ret)) Go To 120
  If (key/='RBB') Go To 130

  Call ffrnum(realo, integ, ret)
  If (on_error(ret)) Go To 120
  If (integ/=2) Go To 130

  Call ffrnum(realo, integ, ret)
  If (on_error(ret)) Go To 120
  ndatos = integ
  Call ffrkey(key, nc, ret)
  If (on_error(ret)) Go To 120
  If (key/='RBB') Go To 130

  Call ffrkey(key, nc, ret)
  If (on_error(ret)) Go To 120

110 Continue

  lineno = lineno + 1

  Call ffrkey(key, nc, ret)
  If (on_error(ret)) Go To 120

  Call ffrnum(realo, integ, ret) !  nline
  If (on_error(ret)) Go To 120


! lectura de los datos
  Do i = 1, ndatos
    Call ffrnum(realo, integ, ret)
    If (on_error(ret)) Go To 120
  End Do

! leer nueva keyword
  Call ffrkey(key, nc, ret)
  If (on_error(ret)) Go To 120

  If (key=='0') Go To 100

  Go To 110

! error handling

120 Select Case (ret)
  Case ('RET1') !                   end of file exit
    Write (*, Fmt='(a)') ' end of file '
    retchar = 'RET1'
    Return
  Case ('RET2') !                   error exit
    Write (*, Fmt='(a)') 'error in ffr subroutines '
    retchar = 'RET2'
    Return
  Case Default
    Stop ' wrong error condition in transic_count'
  End Select

130 If (key=='THEEND') Then
    retchar = 'RET0'
    Return
  Else
    Stop ' wrong end condition in transic_count'
  End If

End Subroutine

!-----------------------------------------------------------------------

Subroutine transic_read(lineno, retchar)

  Use :: nlte_type
  Use :: nlte_dim
  Use :: ffr_error
  Use :: princesa_var, Only: labl, gl
  Use :: cmf_all_var, Only: depacked
  Implicit None

! provides number of lines in file LINES_xxx.dat

! ..
  Integer (i4b), Intent (In) :: lineno
  Character :: retchar*4, ret*4
! ..
! .. local scalars ..
  Real (dp) :: realo, sumglo, xkw, xn, vac
  Integer (i4b) :: integ, nc, ndatos, nline, l
  Character :: key*6

! ..
! .. external subroutines ..
  External :: ffrkey, ffrnum
! ..

  nline = 0

100 Continue
  Call ffrkey(key, nc, ret)
  If (on_error(ret)) Go To 130
  If (key/='CL') Go To 140

  Call ffrkey(key, nc, ret)
  If (on_error(ret)) Go To 130
  If (key/='TY') Go To 140

  Call ffrkey(key, nc, ret)
  If (on_error(ret)) Go To 130
  If (key/='RBB') Go To 140

  Call ffrnum(realo, integ, ret)
  If (on_error(ret)) Go To 130
  If (integ/=2) Go To 140

  Call ffrnum(realo, integ, ret)
  If (on_error(ret)) Go To 130
  ndatos = integ
  Call ffrkey(key, nc, ret)
  If (on_error(ret)) Go To 130
  If (key/='RBB') Go To 140

  Call ffrkey(key, nc, ret)
  If (on_error(ret)) Go To 130

110 Continue

  nline = nline + 1
  depacked(nline)%levlo = key

  Call ffrkey(key, nc, ret)
  If (on_error(ret)) Go To 130
  depacked(nline)%levup = key

  Call ffrnum(realo, integ, ret) !  nline
  If (on_error(ret)) Go To 130
  depacked(nline)%lineno = integ

! read remaining data
  Call ffrnum(realo, integ, ret)
  If (on_error(ret)) Go To 130

! convert air to vacuum (for lambda > 2000)
  If (realo>2000.) Then
!   two iterations are sufficient (until 3rd digit in Angstrom)
    xkw = 1.D4/realo
    xn = 1.D0 + 1.D-7*(643.28D0+294981.D0/(146.D0-xkw**2)+2554.D0/(41.D0-xkw** &
      2))
    vac = realo*xn !                first guess for vacuum wavelength;
!   actually, xn needs to be calculated from vacuum
    xkw = 1.D4/vac
    xn = 1.D0 + 1.D-7*(643.28D0+294981.D0/(146.D0-xkw**2)+2554.D0/(41.D0-xkw** &
      2))
    realo = realo*xn !              2nd iteration
  End If

  depacked(nline)%wave = realo

  Call ffrnum(realo, integ, ret)
  If (on_error(ret)) Go To 130
  depacked(nline)%flu = realo

  Call ffrnum(realo, integ, ret)
  If (on_error(ret)) Go To 130
! gammal = realo

  Call ffrnum(realo, integ, ret)
  If (on_error(ret)) Go To 130

  Call ffrnum(realo, integ, ret)
  If (on_error(ret)) Go To 130

  Call ffrnum(realo, integ, ret)
  If (on_error(ret)) Go To 130
! gammac = realo

  Call ffrnum(realo, integ, ret)
  If (on_error(ret)) Go To 130
  sumglo = realo
! consistency test
  Do l = 1, id_llevs
    If (depacked(nline)%levlo==labl(l)) Then
      If (sumglo/=gl(l)) Then
        Write (999, *) labl(l), sumglo, gl(l)
        Write (999, *) &
          ' stop: inconsistent statistical weights in LINES_xxx.dat'
        Print *, labl(l), sumglo, gl(l)
        Stop ' inconsistent statistical weights in LINES_xxx.dat'
      End If
      Go To 120
!     NOTE: there might be more levels in LINES_xxx.dat than in model atom
    End If
  End Do
! for tests
! print*,' level ',depacked(nline)%levlo,' not found in model atom'

! leer nueva keyword
120 Call ffrkey(key, nc, ret)
  If (on_error(ret)) Go To 130

  If (key=='0') Go To 100

  Go To 110

! error handling

130 Select Case (ret)
  Case ('RET1') !                   end of file exit
    Write (*, Fmt='(a)') ' end of file '
    retchar = 'RET1'
    Return
  Case ('RET2') !                   error exit
    Write (*, Fmt='(a)') 'error in ffr subroutines '
    retchar = 'RET2'
    Return
  Case Default
    Stop ' wrong error condition in transic_read'
  End Select

140 If (key=='THEEND') Then
    retchar = 'RET0'
    If (nline/=lineno) Stop ' nline ne lineno in subr. transic_read'
    Return
  Else
    Stop ' wrong end condition in transic_read'
  End If

End Subroutine

!-----------------------------------------------------------------------

Subroutine sort_line_list(n, xl, gf, id)

  Use :: nlte_type
  Implicit None


! sorts for absolute values


! .. scalar arguments ..
  Integer (i4b) :: n
! ..
! .. array arguments ..
  Real (sp) :: xl(n), gf(n)
  Integer (i4b) :: id(n)
! ..
! .. local scalars ..
  Real (sp) :: xxl, xxg
  Integer (i4b) :: i, ir, j, l, xxi
! ..

  l = n/2 + 1
  ir = n

100 Continue

  If (l>1) Then
    l = l - 1
    xxl = xl(l)
    xxg = gf(l)
    xxi = id(l)
  Else
    xxl = xl(ir)
    xxg = gf(ir)
    xxi = id(ir)
    xl(ir) = xl(1)
    gf(ir) = gf(1)
    id(ir) = id(1)
    ir = ir - 1

    If (ir==1) Then
      xl(1) = xxl
      gf(1) = xxg
      id(1) = xxi
      Return
    End If
  End If

  i = l
  j = l + l

110 Continue

  If (j<=ir) Then

    If (j<ir) Then
      If (xl(j)<xl(j+1)) j = j + 1
    End If

    If (xxl<xl(j)) Then
      xl(i) = xl(j)
      gf(i) = gf(j)
      id(i) = id(j)
      i = j
      j = j + j
    Else
      j = ir + 1
    End If

    Go To 110

  End If

  xl(i) = xxl
  gf(i) = xxg
  id(i) = xxi

  Go To 100

End Subroutine

!-----------------------------------------------------------------------

Subroutine prep_ng(ilow, imax, r, flag_ng)

! prepare and perform Ng-extrapolation for explicit elements

! note: occ_expl refers to 'old' occupation numbers, calculated
! one iteration before last call of rateeq.
! thus, at itng=1 we read S(n-3) from occexpl
! at itng=2 we read S(n-2) from occexpl
! at itng=3 we read S(n-1) from occexpl
! AND                     S(n)   from unit 17


  Use :: nlte_type
  Use :: nlte_dim, Only: id_nttrd, id_ndept, id_atoms, id_llevs
  Use :: fund_const, Only: hc2
  Use :: fastwind_params, Only: gfmin_ng, no_ng_el, ng_el

  Use :: princesa_var, Only: ifirsl, iong, indat1, data, labat, le, li, gl, &
    labl
  Use :: nlte_var, Only: mmlow, mmup, indexrbb, indxlamc, almost_converged

  Use :: cmf_all_var, Only: indexcmf1, occexpl, indexrbb_inv
  Use :: ng_var, Only: itng, ng, n_ng, index1_ng, index_ng, sl_ng, trans_ng
  Implicit None

  Integer (i4b), Parameter :: nd1 = id_ndept, kel = id_atoms, nrec = id_llevs

  Integer (i4b), Dimension (nd1, kel), Intent (In) :: ilow, imax
  Real (dp), Dimension (nd1), Intent (In) :: r

  Integer (i4b), Dimension (kel) :: nfir, nlas

  Integer (i4b) :: k, ilo, ima, numion, igenio, j, ii, indi, ml, mu, nato, i, &
    indmat, indcmf, ion, ll

  Real (dp), Dimension (nd1) :: jbar, sl1n, sl1n1, sl1n2, sl1n3, w, vec

  Real (dp) :: gf, xlamc, lam3, glow, gup, const, xnl, xnu, errmax

  Real (dp) :: a1, a2, b1, b2, c1, c2, aa, bb, denom, nom1, nom2

  Logical :: first, flag_ng

  Data first/.True./

  ng = .False.

  If (first) Then
    Allocate (index1_ng(id_nttrd))
    first = .False.
  Else
    If (n_ng==0) Then
      If (itng/=1) Then
        Write (999, *) &
          ' stop: error in philosophy: n_ng = 0 and itng ne 1 (prep_ng)'
        Stop ' error in philosophy: n_ng = 0 and itng ne 1 (prep_ng)'
      End If
      Return
    End If
  End If


  itng = itng + 1

  If (itng<1) Return

  If (itng==1) Then
!   set up index-array for resonance lines treated with Ng

!   only for ions which are present everywhere
    Do k = 1, kel
      ilo = maxval(ilow(:,k))
      ima = minval(imax(:,k))
      If (ima<ilo) Then
        Write (999, *) ' stop: ima < ilo in prep_ng, itng=1'
        Stop ' ima < ilo in prep_ng, itng=1'
      End If
      numion = igenio(k, ilo)
      nfir(k) = ifirsl(numion)
      numion = igenio(k, ima) + 1
      If (numion>iong) Then
        Write (999, *) ' stop: error in prep_ng - numion'
        Stop 'error in prep_ng - numion'
      End If
      nlas(k) = ifirsl(numion)
!     print*,k,nfir(k),nlas(k)
    End Do

!   find resonance transitions in ng-selected elements (ng_el)
    Print *
    Print *, &
      ' Resonance lines (explicit elements) potentially extrapolated by Ng'
    n_ng = 0
iiloop: Do j = 1, id_nttrd
      ii = indexrbb_inv(j)
!     ii is the index w.r.t. index1
      indi = indat1(ii)
      ml = mmlow(ii)
      mu = mmup(ii)
      nato = le(ml)
      Do i = 1, no_ng_el
        If (labat(nato)==ng_el(i)) Go To 100
      End Do
      Cycle

100   indmat = indxlamc(j)
      indcmf = indexcmf1(indmat) !  use new index
      If (indcmf<3) Cycle !         NOTE INCLUDE AFTER TESTS

      If (ml<nfir(nato) .Or. ml>=nlas(nato)) Cycle
      If (mu>nlas(nato)) Cycle !    should be read 'ge' instead of 'gt'
      ion = li(ml)
      numion = igenio(nato, ion)
!     print*,ml,mu,numion,ifirsl(ion)
      If (ml/=ifirsl(numion)) Cycle ! only resonance lines
      gf = data(indi+2)*gl(ml)
      If (gf<gfmin_ng) Cycle !      no forbidden lines
!     if HeII 303 should be treated, uncomment following line
!     if(labat(nato).eq.'HE' .and. lablineup(indmat).ne.'HE22  ') cycle
      Print *, labat(nato), ion, ml, mu, data(indi+1)
      n_ng = n_ng + 1
      index1_ng(n_ng) = ii
    End Do iiloop

    If (n_ng==0) Then
      Print *, ' No resonance lines (explicit) to be extrapolated by Ng found'
      Print *
      Return
    End If

    Print *
!   for each extrapolation cycle, new allocation, since n_ng might have
!   changed (new ions)
    If (allocated(index_ng)) Then
      Deallocate (index_ng)
      Deallocate (sl_ng)
      Deallocate (trans_ng)
    End If
    Allocate (index_ng(n_ng))
!   4 entries for SL, the last one for Jbar from calc_tlu
    Allocate (sl_ng(nd1,n_ng,5))
    Allocate (trans_ng(n_ng,2))
    index_ng = index1_ng(1:n_ng)
    sl_ng = 0.D0

    index1_ng = 0
    Do i = 1, n_ng
      ii = index_ng(i)
      indmat = indexrbb(ii)
      index1_ng(indmat) = i
    End Do
  End If

  If (itng>=1 .And. itng<=3) Then
!   calculate and save source-functions from occexpl (see above)
ngloop: Do i = 1, n_ng
      ii = index_ng(i)
      indi = indat1(ii)
      xlamc = data(indi+1)*1.D-8
      lam3 = xlamc**3 !             in cgs

      ml = mmlow(ii)
      mu = mmup(ii)
      glow = gl(ml)
      gup = gl(mu)

      const = hc2/lam3
      Do ll = 1, nd1
        xnl = occexpl(ml, ll)/glow
        xnu = occexpl(mu, ll)/gup
        sl_ng(ll, i, itng) = const/(xnl/xnu-1.D0)
      End Do
      If (itng==1) Then
        trans_ng(i, 1) = ml
        trans_ng(i, 2) = mu
      End If
    End Do ngloop
  End If

  If (itng==3) Then
    If (almost_converged) Then
      Do i = 1, n_ng
!       do not use extrapolated value
!       in subr. opacity_cmf and calc_tlu, update of
!       sline -> sl_ng only performed when sl_ng(35,i,4) > 0.
        sl_ng(:, i, 4) = 0.
      End Do
      Print *
      Print *, ' Source-functions for selected resonance lines &
        &(expl.) NOT extrapolated,'
      Print *, ' since ALMOST_CONVERGED previously set to T!'
      Print *
      itng = -2 !                   next two iterations without Ng
      ng = .True.
      Return

    Else
!     at first, we read the actual SL from last call of rateeq
!     here, we have to change the order (outer loop = depth)
      Do ll = 1, nd1
!       note that occexpl (cmf_all_var) is overwritten; this should not
!       matter though, since at this state it's no longer used!
!       for next iteration (opacities etc., it's re-read again anyhow)
        Read (17, Rec=ll)(occexpl(ii,ll), ii=1, nrec)

ngloop1: Do i = 1, n_ng
          ii = index_ng(i)
          indi = indat1(ii)
          xlamc = data(indi+1)*1.D-8
          lam3 = xlamc**3 !         in cgs

          ml = mmlow(ii)
          mu = mmup(ii)
          glow = gl(ml)
          gup = gl(mu)

          const = hc2/lam3

          xnl = occexpl(ml, ll)/glow
          xnu = occexpl(mu, ll)/gup
          sl_ng(ll, i, 4) = const/(xnl/xnu-1.D0) ! this is the actual SL
!         (index=4)
        End Do ngloop1
      End Do

!     extrapolate latest 4 source-functions to new one using Ng-extrapolation
!     inlined to save time
      Do i = 1, n_ng
        ii = index_ng(i)
        ml = mmlow(ii)
        mu = mmup(ii)
        jbar = sl_ng(:, i, 5)
        sl1n = sl_ng(:, i, 4)
        sl1n1 = sl_ng(:, i, 3)
        sl1n2 = sl_ng(:, i, 2)
        sl1n3 = sl_ng(:, i, 1)
!       round to 4 digits, to avoid numerical problems
        Call round(sl1n, 4, nd1)
        Call round(sl1n1, 4, nd1)
        Call round(sl1n2, 4, nd1)
        Call round(sl1n3, 4, nd1)
!       not calculated or inversion, set extrapolated value to zero
!       in subr. opacity_cmf and calc_tlu, update of
!       sline -> sl_ng only performed when sl_ng(35,i,4) > 0.

        If (minval(sl1n)<=0.D0) Then
          Write (*, Fmt='("sl1n le 0, no ng-extrapolation:",i3,2(2x,a6))') i, &
            labl(ml), labl(mu)
          sl_ng(:, i, 4) = 0.
          Cycle
        End If
        If (minval(sl1n1)<=0.D0) Then
          Write (*, Fmt='("sl1n1 le 0, no ng-extrapolation:",i3,2(2x,a6))') i, &
            labl(ml), labl(mu)
          sl_ng(:, i, 4) = 0.
          Cycle
        End If
        If (minval(sl1n2)<=0.D0) Then
          Write (*, Fmt='("sl1n2 le 0, no ng-extrapolation:",i3,2(2x,a6))') i, &
            labl(ml), labl(mu)
          sl_ng(:, i, 4) = 0.
          Cycle
        End If
        If (minval(sl1n3)<=0.D0) Then
          Write (*, Fmt='("sl1n3 le 0, no ng-extrapolation:",i3,2(2x,a6))') i, &
            labl(ml), labl(mu)
          sl_ng(:, i, 4) = 0.
          Cycle
        End If

        Do ll = 1, nd1
          If (jbar(ll)<=0.D0) Then
            Write (999, *) i, ll
            Write (999, *) ' stop: jbar le 0 in prep_ng (itng=3)'
            Print *, i, ll
            Stop ' jbar le 0 in prep_ng (itng=3)'
          End If
        End Do
!       weights according to tests (read_it_ng.pro)
        w = 1./(jbar*r*r)**2
!       w=1./(jbar*r*r)
!       for tests
!       w(:)=1.d0

        vec = (sl1n-2.D0*sl1n1+sl1n2)**2*w
        a1 = sum(vec)
        vec = (sl1n-sl1n1-sl1n2+sl1n3)*(sl1n-2.D0*sl1n1+sl1n2)*w
        b1 = sum(vec)
        a2 = b1
        vec = (sl1n-sl1n1-sl1n2+sl1n3)**2*w
        b2 = sum(vec)
        vec = (sl1n-2.D0*sl1n1+sl1n2)*(sl1n-sl1n1)*w
        c1 = sum(vec)
        vec = (sl1n-sl1n1-sl1n2+sl1n3)*(sl1n-sl1n1)*w
        c2 = sum(vec)

!       modified for denom=0 by Jo Oct. 2025
        denom = (a1*b2-a2*b1)
        nom1 = (c1*b2-c2*b1)
        nom2 = (c2*a1-c1*a2)
        If (denom==0.) Then
          If (abs(nom1)>1.D-5 .Or. abs(nom2)>1.D-5) Then
            Print *, a1, a2, b1, b2, c1, c2
            Print *, denom, nom1, nom2
            Stop ' denom = 0 and nom1 or nom2 ne 0 in prep_ng'
          End If
          sl_ng(:, i, 4) = sl1n
        Else
          aa = nom1/denom
          bb = nom2/denom
          sl_ng(:, i, 4) = (1.D0-aa-bb)*sl1n + aa*sl1n1 + bb*sl1n2
        End If

        If (minval(sl_ng(:,i,4))<=0.D0) Then
          Print *, i, ' extrapolated source function le 0 (prep_ng)'
          Print *, 'NO EXTRAPOLATION!!!'
!         stop ' extrapolated source function le 0 (prep_ng)'
          sl_ng(:, i, 4) = sl1n
        End If
        errmax = maxval(abs(1.-sl_ng(:,i,4)/sl1n(:)))
        Write (*, Fmt='("Ng-extrapolation (expl.): ",i3,2(2x,a6),2x,e10.4)') &
          i, labl(ml), labl(mu), errmax
!       do ll=1,nd1
!       write(*,fmt='(i5,5(2x,e10.4))')
!       ll,sl1n3(ll),sl1n2(ll),sl1n1(ll),sl1n(ll),sl_ng(ll,i,4)
!       enddo
      End Do
      Print *
      Print *, &
        ' Source-functions for selected resonance lines (expl.) extrapolated!'
      Print *
      itng = -2 !                   next two iterations without Ng
      ng = .True.
      flag_ng = .True.

    End If !                        .NOT. ALMOST_CONVERGED
  End If !                          it_ng = 3

  Return
End Subroutine

!-----------------------------------------------------------------------

Subroutine prep_ng_met(r, flag_ng_met)

! prepare and perform Ng-extrapolation for selected background elements

! note: occngold2 refers to 'old' occupation numbers, calculated
! one iteration before last call of rateeq.
! thus, at itng_met=1 we read S(n-3) from occngold2
! at itng_met=2 we read S(n-2) from occngold2
! at itng_met=3 we read S(n-1) from occngold2
! AND                         S(n)   from occng
  Use :: nlte_type
  Use :: fund_const, Only: hc2
  Use :: nlte_dim, Only: id_ndept
  Use :: fastwind_params, Only: names1, natom, gfmin_ng, no_ng_el, ng_el
  Use :: nlte_var, Only: lwion1, metconv2
  Use :: nlte_app, Only: met_imin, met_imax, met_rbb, jatom_full, indexel, &
    indrec, irbb, ilin, occng, occngold2
  Use :: ng_var, Only: itng_met, ng_met, n_ng_met, index_ng_met, sl_ng_met, &
    trans_ng_met
  Use :: nlte_lines, Only: metconv1

  Implicit None

  Integer (i4b), Parameter :: nd1 = id_ndept

  Real (dp), Dimension (nd1), Intent (In) :: r

  Integer (i4b) :: lw1, ilow, imax, k, j, i, n, jj, ii, ji, iqual, ml, mu, &
    ind, irest, ll

  Real (dp), Dimension (:), Allocatable, Save :: jbar, sl1n, sl1n1, sl1n2, &
    sl1n3, w, vec

  Real (dp) :: gf, xlamc, lam3, const, xnl, xnu, errmax

  Real (dp) :: a1, a2, b1, b2, c1, c2, aa, bb

  Logical :: first, flag_ng_met

  Data first/.True./

  ng_met = .False.
  lw1 = lwion1 - 1

  If (first) Then
    Allocate (jbar(lw1), sl1n(lw1), sl1n1(lw1), sl1n2(lw1), sl1n3(lw1), &
      w(lw1), vec(lw1))
    first = .False.
  Else
    If (n_ng_met==0) Then
      If (itng_met/=1) Then
        Write (999, *) ' stop: error in philosophy: n_ng_met = &
          &0 and itng_met ne 1 (prep_ng_met)'
        Stop &
          ' error in philosophy: n_ng_met = 0 and itng_met ne 1 (prep_ng_met)'
      End If
      Return
    End If
  End If

  itng_met = itng_met + 1

  If (itng_met<1) Return

  If (itng_met==1) Then
!   set up index-array for resonance lines treated with Ng

    Print *
    Print *, &
      ' Resonance lines (background elements) potentially extrapolated by Ng'
    n_ng_met = 0
    Do k = 1, natom
      If (indexel(k)>0) Cycle !     explicit element
      If (jatom_full(k)==0) Cycle ! only selected elements

!     find resonance transitions in ng-selected elements (ng_el)
      Do i = 1, no_ng_el
        If (names1(k)==ng_el(i)) Go To 100
      End Do
      Cycle

!     only for ions which are present everywhere
100   ilow = maxval(met_imin(k,1:lw1))
      imax = minval(met_imax(k,1:lw1))

      Do j = ilow, imax
!       hack
!       if(j.ne.3) cycle
        n = indrec(k, j)
        jj = irbb(n) - 1

        Do ii = 1, ilin(n)
          ji = jj + ii
          iqual = met_rbb(ji)%quality
!         outside range or no transition in line-list,
!         will be treated by conventional (simplified) approach in netmat_met,
!         see calc_tlu_met
          If (iqual==-1 .Or. iqual==2) Cycle

          ml = met_rbb(ji)%low
          If (ml/=1) Cycle !        only resonance lines
          gf = met_rbb(ji)%gf
          If (gf<gfmin_ng) Cycle !  no forbidden lines
          mu = met_rbb(ji)%lup
!         hack
!         if(mu.ne.10) cycle
          xlamc = met_rbb(ji)%wave ! packed transition, in A
          ind = k*1000000 + j*100000 + ml*100 + mu
          Print *, names1(k), ind, xlamc
          n_ng_met = n_ng_met + 1
          If (n_ng_met>1000) Then
            Write (999, *) &
              ' stop: n_ng_met > 1000! Enlarge array trans_ng_met!'
            Stop ' n_ng_met > 1000! Enlarge array trans_ng_met!'
          End If
          met_rbb(n_ng_met)%index1_ng = ji ! used as dummy-variable
          trans_ng_met(n_ng_met) = ind
        End Do
      End Do
    End Do

    If (n_ng_met==0) Then
      Print *, &
        ' No resonance lines (background) to be extrapolated by Ng found'
      Print *
      Return
    End If

    Print *
!   for each extrapolation cycle, new allocation, since n_ng_met might have
!   changed (new ions)
    If (allocated(index_ng_met)) Then
      Deallocate (index_ng_met)
      Deallocate (sl_ng_met)
    End If
    Allocate (index_ng_met(n_ng_met))
!   4 entries for SL, the last one for Jbar from calc_tlu
    Allocate (sl_ng_met(lw1,n_ng_met,5))
    index_ng_met = met_rbb(1:n_ng_met)%index1_ng
    sl_ng_met = 0.D0

    met_rbb%index1_ng = 0
    Do i = 1, n_ng_met
      ji = index_ng_met(i)
      met_rbb(ji)%index1_ng = i
    End Do
  End If

  If (itng_met>=1 .And. itng_met<=3) Then
!   calculate and save source-functions from occexpl (see above)
ngloop: Do i = 1, n_ng_met
      ji = index_ng_met(i)
      xlamc = met_rbb(ji)%wave*1.D-8
      lam3 = xlamc**3 !             in cgs

      ind = trans_ng_met(i)

      k = ind/1000000
      irest = ind - k*1000000
      j = irest/100000
      irest = irest - j*100000
      ml = irest/100
      mu = irest - ml*100

      n = indrec(k, j)
      const = hc2/lam3

      Do ll = 1, lw1
!       if(ll.eq.1.or.ll.eq.lw1) &
!       &    write(*,fmt='(i3,3(2x,e10.4))')
!       i,1.-occngold2(n,mu,ll)/occng(n,mu,ll), &
!       &    occngold2(n,mu,ll),occng(n,mu,ll)
        xnl = occngold2(n, ml, ll) ! occng = n/g
        xnu = occngold2(n, mu, ll)
        sl_ng_met(ll, i, itng_met) = const/(xnl/xnu-1.D0)
      End Do

      If (itng_met==3) Then
!       here, we additionally use the most recent occupation numbers to obtain
!       sl1n
        Do ll = 1, lw1
          xnl = occng(n, ml, ll) !  occng = n/g
          xnu = occng(n, mu, ll)
          sl_ng_met(ll, i, 4) = const/(xnl/xnu-1.D0)
        End Do
      End If

    End Do ngloop
  End If

  If (itng_met==3) Then
!   JO Jan. 2021; changed from metconv1 to (metconv1 .or metconv2);
!   tests have shown that for metconv2 = T the convergence is better without
!   extrapolation;
!   in those cases where metconv2 does not appear (opt_oiii_iter=F), metconv1
!   is still used as the decisive critirion.
!   if(metconv1) then
    If (metconv1 .Or. metconv2) Then
      Do i = 1, n_ng_met
!       do not use extrapolated value
!       in subr. opacity_cmf and calc_tlu, update of
!       sline -> sl_ng only performed when sl_ng(35,i,4) > 0.
        sl_ng_met(:, i, 4) = 0.
      End Do
      Print *
      Print *, ' Source-functions for selected resonance lines &
        &(backgr.) NOT extrapolated,'
      Print *, ' since METCONV1 or METCONV2 previously set to T!'
      Print *
      itng_met = -2 !               next two iterations without Ng
      ng_met = .True.
      Return

    Else
!     extrapolate latest 4 source-functions to new one using Ng-extrapolation
!     inlined to save time
      Do i = 1, n_ng_met
        ji = index_ng_met(i)
        ind = trans_ng_met(i)
        k = ind/1000000
        irest = ind - k*1000000
        j = irest/100000
        irest = irest - j*100000
        ml = irest/100
        mu = irest - ml*100

        jbar = sl_ng_met(:, i, 5)
        sl1n = sl_ng_met(:, i, 4)
        sl1n1 = sl_ng_met(:, i, 3)
        sl1n2 = sl_ng_met(:, i, 2)
        sl1n3 = sl_ng_met(:, i, 1)
!       round to 4 digits, to avoid numerical problems
        Call round(sl1n, 4, lw1)
        Call round(sl1n1, 4, lw1)
        Call round(sl1n2, 4, lw1)
        Call round(sl1n3, 4, lw1)
!       not calculated or inversion, set extrapolated value to zero
!       in subr. opacity_cmf and calc_tlu_met, update of
!       sline -> sl_ng_met only performed when sl_ng_met(35,i,4) > 0.

        If (minval(sl1n)<=0.D0) Then
          Write (*, Fmt= &
            '("sl1n le 0, no ng-extrapolation:",i3,2x,a2,3(2x,i3))') i, &
            names1(k), j, ml, mu
          sl_ng_met(:, i, 4) = 0.
          Cycle
        End If
        If (minval(sl1n1)<=0.D0) Then
          Write (*, Fmt= &
            '("sl1n1 le 0, no ng-extrapolation:",i3,2x,a2,3(2x,i3))') i, &
            names1(k), j, ml, mu
          sl_ng_met(:, i, 4) = 0.
          Cycle
        End If
        If (minval(sl1n2)<=0.D0) Then
          Write (*, Fmt= &
            '("sl1n2 le 0, no ng-extrapolation:",i3,2x,a2,3(2x,i3))') i, &
            names1(k), j, ml, mu
          sl_ng_met(:, i, 4) = 0.
          Cycle
        End If
        If (minval(sl1n3)<=0.D0) Then
          Write (*, Fmt= &
            '("sl1n3 le 0, no ng-extrapolation:",i3,2x,a2,3(2x,i3))') i, &
            names1(k), j, ml, mu
          sl_ng_met(:, i, 4) = 0.
          Cycle
        End If

!       JO changed Oct. 2025
!       old version commented out
!       do ll=1,lw1
!       if(jbar(ll).le.0.d0) then
!       print*,i,ll
!       stop ' jbar le 0 in prep_ng (itng_met=3)'
!       endif
!       enddo
!       can happen, if select provides ions different ('adapted') from those
!       which were present at set-up
        If (minval(jbar(1:lw1))<=0.D0) Then
          Write (*, Fmt= &
            '("jbar le 0 (check previous run of select):",i3,2x,a2,3(2x,i3))') &
            i, names1(k), j, ml, mu
          Cycle
        End If
!       weights according to tests (read_it_ng.pro)
        w = 1./(jbar*r(1:lw1)*r(1:lw1))**2
!       w=1./(jbar*r(1:lw1)*r(1:lw1))
!       for tests
!       w(:)=1.d0

        vec = (sl1n-2.D0*sl1n1+sl1n2)**2*w
        a1 = sum(vec)
        vec = (sl1n-sl1n1-sl1n2+sl1n3)*(sl1n-2.D0*sl1n1+sl1n2)*w
        b1 = sum(vec)
        a2 = b1
        vec = (sl1n-sl1n1-sl1n2+sl1n3)**2*w
        b2 = sum(vec)
        vec = (sl1n-2.D0*sl1n1+sl1n2)*(sl1n-sl1n1)*w
        c1 = sum(vec)
        vec = (sl1n-sl1n1-sl1n2+sl1n3)*(sl1n-sl1n1)*w
        c2 = sum(vec)
!       Jo Sept 18
!       if identical source-functions (e.g., when bg elements kept constant)
!       somewhat different compared to prep_ng (for explicit elements), but
!       similar result (all source functions are equal)
        If (a1*b2-a2*b1==0.) Then
          If (c1*b2-c2*b1/=0) Stop ' ng-extrapolation of bg: problems with aa'
          If (c2*a1-c1*a2/=0) Stop ' ng-extrapolation of bg: problems with bb'
          aa = 0.5D0
          bb = 0.5D0
        Else
          aa = (c1*b2-c2*b1)/(a1*b2-a2*b1)
          bb = (c2*a1-c1*a2)/(a1*b2-a2*b1)
        End If
        sl_ng_met(:, i, 4) = (1.D0-aa-bb)*sl1n + aa*sl1n1 + bb*sl1n2
        errmax = maxval(abs(1.-sl_ng_met(:,i,4)/sl1n(:)))
        If (minval(sl_ng_met(:,i,4))<=0.D0) Then
          Write (999, *) i, ' extrapolated source function le 0 (prep_ng_met)'
          Write (999, *) ' NO EXTRAPOLATION!!!'
          Print *, i, ' extrapolated source function le 0 (prep_ng_met)'
          Print *, ' NO EXTRAPOLATION!!!'
!         stop ' extrapolated source function le 0 (prep_ng)'
          sl_ng_met(:, i, 4) = sl1n ! rounded
        End If
        errmax = maxval(abs(1.-sl_ng_met(:,i,4)/sl1n(:)))
        Write (*, Fmt= &
          '("Ng-extrapolation (bg.): ",i3,2x,a2,3(2x,i3),2x,e10.4)') i, &
          names1(k), j, ml, mu, errmax
!       do ll=1,lw1
!       write(*,fmt='(i5,5(2x,e10.4))'),ll,sl1n3(ll),sl1n2(ll),sl1n1(ll),sl1n(ll),sl_ng_met(ll,i,4)
!       enddo
      End Do
      Print *
      Print *, ' Source-functions for selected resonance lines &
        &(backgr.) extrapolated!'
      Print *
      itng_met = -2 !               next two iterations without Ng
      ng_met = .True.
      flag_ng_met = .True.

    End If !                        .NOT. METCONV1

  End If !                          it_ng_met=3

  Return
End Subroutine

!-----------------------------------------------------------------------

Subroutine round(x, prec, nd)

! rounds vector x of dimension nd to prec
  Use :: nlte_type

  Integer (i4b), Intent (In) :: nd, prec
  Real (dp), Dimension (nd) :: x

  Integer (i4b), Dimension (nd) :: expo ! automatic array
  Real (dp), Dimension (nd) :: rexp ! automatic array


! take care: difference, since INT works different for positive and negative
! arguments
  Where (x>=1.D0)
    expo = -int(log10(x)) - 1 + prec
  Elsewhere
    expo = -int(log10(x)) + prec
  End Where
  rexp = 10.D0**expo

  x = nint(x*rexp)/rexp
  Return
End Subroutine
