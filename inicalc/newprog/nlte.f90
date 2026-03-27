Module version_nlte
! -----------------------------------------------------------------------

! program history

! version  0.0  june 1995
!   original version by quique

! version  1.0  may 30th 1996
!   new hydro, cmf, number of bugs corrected

! version: 1.1  july 10th 1996
!   modified for mixed cmf/sa ,
!   outer boundary (sa-lines),
!   density inversion, maximum corrections

! version: 1.2  oct  11th 1996
!   modified for dq'/dtau'-term in hydro
!   and hopf-parameter input

! version: 1.3  nov  13th 1996
!  hydrogen collision rates new

! version: 1.4  dec  4th 1996
!   consistent treatment of p vs. rho including mu(r)
!   storage of pressure in model-file
!   check of consistencies of tnew and tact
!   ground1 (iteration scheme)

! version: 1.5  march 11th 1997
!   changed location of 2nd radius point (wrs!)
!   now allowed to start from converged hydro-models
!   if itlast.ne.0

!   !!! note that due to a slightly different updating
!   of xlevel = (ni/nk)*, minor differences in
!   the occupation numbers will arise if one
!   compares the results of a calculation with
!   and without an interupt!!!

! version: 1.6  march 17th 1997
!   z-integration for optically thick lines in cmf
!   changed (z=0 point at p = rmax included).
!   most important!!! consistent alo at r = rmax
!   for optically thick lines developed and included.

! version: 1.7  march 21st 1997
!   first version of mixed sa/cmf approach according
!   to the requirements specified in subroutine
!   cmfprep. this possibility is used if the parameter
!   optmixed in the input (the former parameter df,
!   which has no meaning any longer) is set to unity.
!   if all lines are to be treated in the cmf, optmixed
!   must be zero.

! version: 2.0  may 23rd 1997
!   completely revised, can calculate now all atoms
!   provided by detail.
!   most important changes refer to photo-integrals,
!   which are calculated by a completely new algorithm,
!   and frequency meshes, which are now two!

! version: 2.0.1 june 4th 1997

!   qualifiers for problematic levels introduced
!   occupation numbers and ionization set to zero
!   where outside ilow, imax

! version: 2.0.2 june 5th 1997

!   rybicki scheme in continuum transport only after
!   major changes of occ.no., else pure iteration of
!   eddington factors


! version: 2.0.3 june 17th 1997

!   checked for conistency with old version.
!   found two additional errors in old one:
!   update of electron density for ni* (in netmat,
!   coll. terms) forgotten. inconsistent treatment of
!   gaunt-factor in opacity (all freq.) and photo-
!   integrals (only at edge) leads to inappropriate alo.

!   correction of bug in intermepho (wrong end condition
!   in some cases)

! version: 2.0.4 dec 4th 1997

!   density inversions now allowed!
!   only warning on output!


! version: 2.0.5 jan 9th 1998
!   copy of indat.dat into directory

! version: 2.0.6 feb 10th 1998
!   check for photospheric mass
!   problems for too fast convergence removed

! version: 2.1.0 feb 13th 1998
!   polished, i4 compilation now possible

! version: 2.2.0 march 9th 1998
!   accuracy for matrix solution improved
!   left out line is now that with largest occ.num.
!   in sobo: inversion treated differently

! version: 2.2.1 nov 18th 1998
!   yhe has no longer to be updated in the atomic
!   data input file. abund(he), if existent, is
!   overwritten by yhein(from indat.dat) in "start"

! version: 3.0.0 jan 15th 1999  conversion to f90

! version: 3.0.1 jan 20th 1999  inclusion of work by christian wiethaus
!   micro-turbulence included
!   mg i,ii (simple version) can be treated now

! version: 3.0.2 feb 2nd 1999
!   correction of  (very old) error in ray(=> cmfset)

! version: 3.1 march 4th 1999  almost real f90 now

! version: 3.2 may 6th 1999  line quantities for inverted levels
!   resulting in total negative opacity (line + cont)
!   approximated in a VERY poor way to avoid difficulties
!   in convergence

! version: 4.0 feb 2000 including nlte_approx
!   (Note: xmet added to input)

!   bug in control of gff(-> i) fixed.
!   (note: no effect for all models calculated so far,
!   due to same charge of HII and HeI)

! version: 4.1 may 2000 final version including approx. metal
!   background. Valid for O/B stars.
!   to rund old version, simply set xmet = 0.

! version: 4.2 june 7th 2000 catalog name (modnam) extended up
!   to 25 characters

! version: 4.2.1 oct 18th 2000 improvement of modnam extension (30 char)

! version: 4.3 jan 24th 2001 minimum freq. in frescal/fresfin
!   automatically calculated.
!   Sobo_tranport modified with respect to calc. of bcIc
!   in cases of taus ne O(1)
!   Test version of SA transport for inverted levels

! version: 5.0 feb 7th 2001 internal beta-Iteration in SA-Transport
!   abolished
!   All tables for U and U2 functions moved to DATA

! version: 5.1 feb 20th 2001 Sobo_transport for inverted levels finalized
!   dilfac1 modified

! version: 5.1.1 march 1st 2001 : Limits for FASTBCIC in TRATRAD changed!
!   If convergence fails for some models with pure SOBO- &
!   transport, check limits again.

! version: 6.0 minor modifications to allow for inclusion of averaged
!   metal line background

! version: 6.1 april 2nd 2001: finalized version for metallic line
!   backgroud: input in INDAT.DAT modified, requires now
!   also input of LINES and LINES_IN_MODEL (cf. docu.txt)

! version: 6.2 april 10th 2001: cmf-treatment of inverted levels
!   now "exact", smaller corrections to allow for -nan comp.

! version: 6.2.1 april 23rd 2001: taumax for inversion changed,
!   CMF-lines until 150000 A to include H5->6 transition

! version: 6.3 april 24th 2001: inclusion of H- opacity, higher resolution
!   in photosphere possible if dimes.f90 changed (nd->67),
!   meaning of first logical variable in INDAT.DAT changed:
!   OPTNEUPDATE = .TRUE., as usual
!   OPTNEUPDATE = .FALSE., use always LTE electron density (e.g.,
!   for high gravity A-stars)
!   Rayleigh scattering can be included if lines including
!   "rayleigh" are uncommented

! version: 6.4 may 2nd 2001: can be used now to use detailled
!   H/He atomic model A10HHe: all required formulas implemented

! version: 6.5 may 14th 2001: error in line transfer (update of
!   cont. source function) found and corrected

! version: 6.6 may 18th 2001: freq. grid changed to resolve line-blocked
!   regions

! version: 6.7 june 21st 2001: Two Hopfparameter files, in dependence of xmet

! version: 7.0 july 25th 2001: line overlap included/checked

! version: 7.1  nov  7th 2001: log q in Hopfpara-files assumed

! version: 7.1.1 april 18th 2002: subroutine hopfpara changed
!   change of subroutine prepov, in order to allow
!   for Sobolev lines inside overlap complex

! version: 7.1.2 april 26th 2002: INDAT.DAT used from cat

! version: 7.2  july 3rd 2002: vdiv may adapt now in MODVL to
!   prevent inversion
!   additional file FASTWIND.LOG

! version: 7.2.1 march 2003, Tcorr almost finalized
!   more metals (see below)

! version: 7.2.2 april 2003, Tcorr finalized
!   including line-blocking (approx.treatment)tcorrvar

! version: 7.3   alternate path: modification for new nlte_approx
!   october 10th 2002: xxk (2nd moment) calculated
!   always (req. for eddf)
!   different treatment for LINES_IN_MODEL = .FALSE.
!   model calculated without background-lines, however
!   final temperature with lines (allows for a comparsion
!   with WM-basic); Restart capabilities not tested so far
!   used only by Jo.

! version: 7.4   june 2003: Vers. 7.2.2 and 7.3 coupled

! version: 8.0   july 1st 2003: most important metals now "exactly treated, &
!   T-correction from corresponding transitions
!   needs nlte_approx Vers. 5.0,
!   (-> new atomic data file-path: ATOMDAT_NEW!)

! version: 8.1   july 30th: small changes after testing for Bstars;
!   most important changes: check for TAURLIM and increase
!   if strong winds; imin, imax for Si.
!   check also for TMIN1 in module fastwind_params!

! version: 8.2   inclusion of approximate treatment of clumping
!   different layout of INDAT.DAT required
!   use nlte_approx version 5.1 and higher

!   NOTE1: ne, nh and occup. numbers include enhancement-factor,
!   opacities are corrected
!   (except FFOPA, which is needed only for T-structure)
!   to account for spatial averages in RT calculations.
!   PAY ATTENTION if line processes are considered, 
!   particularly in SA, and
!   whenever local opacities are needed
!   WARNING: since the energy balance is calculated inside
!   the clumps (alternative argument: enhancement-factor
!   cancels out), pay attention if adiabatic cooling is
!   included, since then the rates MUST be corrected!
!   NOTE2: not rigorously tested so far. Use at own risk!

! version: 8.3   changes in INTERMEPHO (calculate heating/cooling rates
!   without ALO)
!   XNE -> XNELTE whenever LTE values are used
!   (to clarify programming)
!   heating/cooling rates of explicit elements calculated
!   with NEW occupation numbers (which is more consistent)
!   START of Sobolev Treatment from it = 4 on!

! version: 8.3.1 improvement of T-convergence, hundreds of tests

! version: 8.3.2 May 2004
!   inclusion of DELPARA in calculation of hydro-structure
!   to account for deviations from power-law in case of
!   cool stars (Teff < 12000). Now the actual and used
!   opacities should agree perfectly.
!   Note that in this case vdiv is NOT updated, &
!   thus density inversions are possible.
!   Bug in logic concerning OPTNEUPDATE=.FALSE. and T-corr
!   removed. Ne is now updated to actual LTE density.
!   New hack in the calculation of T(hydro). If Tact
!   is still Tmin, dT is set to zero. This improves the
!   consistency between Tact and Thydro significantly.
!   The deviation between both (-> FASTWIND.LOG) should be
!   now of order of 1%.

! version: 8.4   August 2004
!   inclusion of formula 46 plus function prhcol
!   (n, np, z, e), as referred to by
!   Przybilla and Butler, A&A 2004.
!   New formula 113 (old H-coll, cf. Burke et al. 1967),
!   corresponds to the data given in A10HHe if the
!   modification as outlined by formula 13 (see below)
!   should NOT be applied.
!   This version is applicable to both the optical and the
!   infrared provided that the correct DETAIL DATASET
!   is taken; i.e., either
!   A10HHe (old H-coll. data from Giovanardi et al., &
!   subr. hcol, formula 13)
!   or
!   A11HHe (new H-coll. data from Butler, see above)

! version: 8.5   Nov. 2004
!   change in calculation of dqdt partials for explicit
!   elements to be consistent with nlte_approx.
!   smaller changes in treatment of inverted levels

!   changes in frescal and fresfin to resolve H and HeI
!   resonance lines

! version: 8.5.1 April 2005: Bug in Tratcol (formula 113) identified

! version: 8.5.2 Sept. 2005: Rel diff in cooling changed: at first
!   differences > 1e-3 are checked, then differences > 1e-4
!   abundance changes for explicit and background elements
!   now from INDAT.DAT (see comments around line 870 and docu.txt
!   NEW solar abundances from file abundances_solar

! version: 8.6   Oct. 2005: Careful tests concerning clumping performed;
!   no programming error in original approach from
!   version 8.2 detected.
!   Input of clumping-parameters generalized:
!   Variable number of parameters (up to ten) possible
!   (e.g., if user changes stratification in subroutine clumping)
!   parameters now contained in allocated field clf_field)

!   IMPORTANT: first parameter has to be REPRESENTATIVE
!   clumping factor which will be ONLY used to find
!   appropriate HOPF-parameters for calculating the
!   T-structure in the start-model.

!   Correspondance with old approach:

!   clf_field(1)=clf, clf_field(2)=vclstart, clf_field(3)=vclmax

! version: 8.6.1 Oct. 2005: distribution of radial points changed,
!   to better sample the transonic region (for IR!)
!   xnh precision in GROUND1 changed, to allow for restart
!   with different (explicit) metal abundances

! version: 8.6.2 Feb. 2006: formula 17 - CBF (needed, e.g., for ground
!   state ionization of NIII) changed. Note that
!   there was a bug in the old version!

!   precision of inversion of pure continuum rate eq.
!   increased, in analogy to line case

!   pure continuum iteration only for H/He


! version: 8.7 March 2006: inclusion of Miguel's version 8.5.2_alpha
!   treatment of the OP rbf data improved: a new
!   subroutine (op_setup) takes care of the initial
!   set up of the OP data. Minor modifications to account for
!   changes in the frequency grid (see LTE subroutine).
!   It also checks for differences in the ionization energies
!   defined in the input DETAIL file and the OP energies, for
!   consistency.
!   Bug detected and corrected in TRATCOL formula (26), so far
!   used only by CIII model atom.
!   Bug in INTERMECOL formula (17) corrected (Achim),
!   only affects NIII model atom, so far
!   Integration limit DEX (subroutine INTERMEPHO) changed to
!   allow for very detailed model atoms (mainly Norbert's neutral
!   ions) that have levels with fairly similar energies.

! version: 8.7.1 March 2006: smoothing of OP-DATA RBF cross-sections

! version: 9.0 March 2006: REQUIRES NEW nlte_type.f90
!   new formualae for Adi's data set
!   completely new approach to calculate ILOW, IMAX
!   new routines ilowimax_lte and ilowimax_nlte
!   in case, update ONLY ILOWIMAX_NLTE (contained in nlte_approx.f90)
!   The old treatment (which is used also for z=0) is
!   kept in routine ILOWIMAX_OLD (requires "old" detail files)
!   Now, we make a division between
!   ILOW, IMAX, IMIA, IMAA
!   which refer to LTE values and define, e.g., the
!   frequency grid, and
!   ILOWNL, IMAXNL, IMIANL, IMAANL
!   which refer to NLTE values and are used to solve the
!   rate equations. For H and He, both quantities are
!   identical and calculated in OCCLTE.
!   Dimensions of xl, partit and occnum corrected.

!   Note: z=0 treatment not checked so far

!   default value for nd (id_ndepth) is now 51. One point
!   close to photosphere (just below) included, to allow
!   for appropriate thermalization of background elements

!   iterative improvement of solution of rate-equations
!   calculated with residuum in qp precision
!   (possible with recent versions of intel compiler)

! version: 9.1 August 2006: subr. ILOWIMAX_NLTE moved to nlte_approx in
!   order to allow compilation from scratch (uses nlte_app)
!   VDIV allowed to INCREASE only (old "bug")

! version: 9.2 August 2006: subr. ILOWIMAX_NLTE moved to nlte_approx

! version: 9.2.1 Nov. 2006 ensure convergence of oscillating cool/high gravity
!   models.
!   Call to TLUCY in MODVL (with e-scat only) only for first
!   Hydro-Interation

! version: 9.2.2 Dec. 2006: old bug in optsing = .false.
!   (multi-line-approach) fixed

! version: 9.2.3 Jan. 2007: improved calc. of dvdr

! version: 9.2.4 March 2007: crossbf formula 101: read tabulated
!   (OP-) photocross-sections provided by Norbert P.
!   and convolve with Gaussian
!   small changes in 'smooth' to account for this

! version: 9.2.5 July 2007: photospheric grad written to GRAD.OUT

! version: 9.2.6 October 2007: width of zones with resonances increased
!   from 10 to 20 times emin (in subroutine smooth)
!   to allow for correct treatment of SiII bf cross-sections

! version: 9.3 Sept 2008: detail/kurucz models can be read in if
!   OPTMOD = TRUE in the first iteration.
!   Filename has to be Modelname+'.struct'
!   New module photstruc

! version: 9.4 Feb 2009: tlusty models can be read in, in the same
!   spirit as above. Here, we calculate P from rho instead
!   rho from P.
!   Filename has to be Modelname+'.tlusty'
!   module photstruc updated

! !!! NOTE: SO FAR, EXTERNAL PHOTOSPHERIC STRUCTURES
!   ONLY FOR UNCLUMPED MODELS !!!!
!   Warning: if you get the message
!   ERROR IN R(TAU-LUCY=2/3),
!   then most probably Rstar is too small!

! version: 9.5 June 2009: lots of smaller changes to allow a consistent
!   treatment for explicit elements with MORE than three ions
!   in particular, all occupation numbers outside the local
!   values of ILOWNL, IMAXNL are set to zero. Only at restart
!   or Temp-update, these numbers are set to LTE*1.D-10, to
!   allow for calculation of meaningful bf/ff opacities.

!   NOTE: NEW parameter in subroutine ILOWIMAX_NLTE (package
!   nlte_approx), to control globally constant or locally
!   consistent values of ILOWNL, IMAXNL. Default is global.

! version: 9.6 July 2009: inclusion of dielectronic recombination for
!   explicit ions following the method by Mihalas, Hummer etc.
!   RBF formula 21: PhotoData described by smooth cross-sections
!   in Seaton parameterization, resonances described as
!   stabilizing lines from the double excited configuration
!   (fromulated with respect to the ground-state of the next higher ion)

! version: 9.7 July 2009: new formulas 60...65 in tratcol (for NIII)
!   interpolation of OP-data where zero cross-section
!   Hack in subroutine THERMAL BALANCE (to avoid stop because
!   of inaccurate flux-integration - done in Sept. 2009)

! version: 9.8 deprecated

! version: 9.9 Sept. 2009: all global options put into module nlte_opt
!   needs nlte_type v1.2

!   inclusion from work done in March 2009:
!   small changes required to use
!   nlte_approx v7.0 and higher (photospheric line transfer)
!   weighting scheme routine tempcorr slightly changed,
!   towards larger influence of dT(flux-correction)

! version: 10.0 Oct. 2009: update of photo-structure
!   new module RUN_ONCE with new variables
!   (replacing START and FIRST in the routines
!   WEIGH11, CMFSING, CMFMULTI (nlte.f90) and
!   jbar_metals, cmf_simple (nlte_approx,f90) introduced
!   calling sequence of WEIGH11 changed

!   number of smaller topics changed
!   GLOBAL (controlling ILOWNL, IMAXNL) put to nlte_opt
!   EXPIN: argument zero explictly treated
!   check in ILOWIMAX_LTE: if blevel (HE) is zero, &
!   update to a small number

! version : 10.1 Dec. 2009: formal solution (jbar and alo) of RT
!   for photospheric lines changed (nlte_approx, V7.1)
!   new method uses differential formulation
!   switch between new (differential) and old (integral)
!   formulation with OPTPHOTLINES_DIFF = TRUE/FALSE
!   in module nlte_opt.

! subroutine TCORR: calculation of dFerr/dtau improved
!   for closely separated grid-points

! subroutine MODVL and MODVL_UPDATE: check that point
!   just below SRNOM lies not too close to previous point.
!   in case, cure this problem by skipping the previous point.

!   new file CONVERG_METALS: controls convergence of metals

!   if UPDATE_STRUCT = TRUE,
!   convergence criterium for first hydro-iterations changed:
!   if ITHYD > 20 and abs(corr) > 0.97, model considered as
!   converged (usually, < 10 iterations required).
!   Should not matter since model updated anyway.

!   test of photospheric clumping with rho_phot = rho_phot/fcl
!   (constant clumping factor) enabled if OPT_PHOTCLUMP:
!   two models (unclumped model 0 and clumped model 1 with
!   FCL > 1 but constant) should be identical with respect
!   to model structure, ionization fraction etc
!   (plotted as a function of M or TAUR) and profiles
!   if Mdot1 = fcl Mdot0 and Rstar1 = fcl Rstar0

!   in MODVL, accuracy w.r.t. location of v < vs2 changed

!   new formula 66 in tratcol (for NIV)
!   old bug in FRESFIN corrected
!   for restarted models with different metallicity,
!   metallic background now recalculated (see nlte_approx)

!   no update of phot. structure performed for dense winds
!   where radacc(start) > radacc(actual)
!   at the outer photosphere (wind acceleration becomes of importance)

!   subroutine ILOWIMAX_OLD updated w.r.t. new N atom

! version : 10.1.1 April 2010: Also cmfgen input can be used (with ending
!   .tlusty), see subroutine MODVL (extrapolation of
!   input structure now allowed.)

! version : 10.1.2 Jan 2011: three changes of precision (for intel 11.0)

! version : 10.1.3 March 2012: can use also Michel Cure's input.
!   NO update of phot. structure
!   OPTMOD = TRUE in the first iteration.
!   Filename has to be Modelname+'.michel'
!   NOTE: The input value of beta has to be representative
!   for Michel's v-field
!   Bug in MODVL_UPDATE fixed: SIGMA_THOMSON needs still
!   to be calculated from XNE(LTE) to be consistent with TAUR
!   (problems occured for cool stars).
!   New philosophy in MODVL and MODVL_UPDATE to cure problems
!   with point close to SRNOM (2nd try with safety-factor 1.2)
!   version : 10.1.4 April 2012: inclusion of Paco's clumping law with
!   4 parameters (CL1 ... CL4)
!   small bug in FRESCAL removed, regarding LAM_LY

! version : 10.1.5 April 2013: Improved treatment of freq. spacing between
!   911 and 925 A, to ensure a maximum separation of 3 A
!   (minimum box size for nsubinter = 120  is 3.6 A, can lead
!   to problems with nlambox_ly when many ions included)
!   -> changes in FRESCAL (also in dimes.f90)

! version : 10.1.6 July 2013: Update of routines that convolve
!   OP-cross-sections with a Gaussian (mainly SMOOTH and
!   CONVOLVE, but also OP_RBFSETUP, used by formula 20), &
!   to allow for a correct treatment of Miguel's data
!   (e.g., MgI and MgII)

! version : 10.2.0 October 2014: inclusion of wind-embedded X-rays
!   new module nlte_xrays, new package lambda_xray.f90,
!   new data file k_shell_Auger_data
!   X-ray emission following Raymond & Krolik,
!   shock cooling zones following Feldmeier et al. (1998)
!   requires additional input with keyword XRAYs (see docu.txt)
!   K-shell Auger ionization if OPTAUGER = .TRUE. (default)
!   L/M-shell Auger ionization for higher elements missing ..,
!   L_x/L_bol calculated and saved (file XLUM),
!   in range 100 eV to 2.5 keV
!   Tested: explicit and bg elements (if consistent) give very
!   similar results for ionization fractions
!   Tested: HHe and other models give very similar results
!   still to do: change ne and nh to account for hot plasma

! version : 10.3.0 October 2014: inclusion of optically thick clumping
!   (porosity and vorosity), programmed by Jon Sundqvist
!   included into version 10.2.0
!   still to do: improve function CLUMPING
!   (input, and factor (1-(1-fvol)fic for cont. porosity)

! inclusion of v10.1.7 and v10.1.7.1
!   Update of FRESFIN and FRESCAL to improve
!   resolution between Lyman beta and Lyman-edge
!   some bug-fixes due to compiler issues


! FROM HERE ON, DIFFERENT PATHES FOR "old = standard" and "new" version.

! version : 10.4.0 Jan 2015: New approach with complete CMF solution
!   included from previous versions

!   NOTE: porosity and vorosity still to be included in cmf_all

! version : 10.4.1 Dec 2015: inclusion of v10.1.7.2 (cure bug in v10.1.7.1)

! version : 10.4.2 Jan 2016: couple of smaller changes to improve diffusion
!   approx. at lower boundary;
!   most important: XM(ND) = XM(ND-1)*1.1 in MODVL and MODVL_UPDATE

! version : 10.4.3 April 2016: more smaller changes, in particular
!   TAUR (CMF) now used for everything when OPTCMF_ALL=.TRUE.
!   new lower boundary only if ND ge 61 -- for lower number
!   of grid-points, deep photosphere badly resolved, and
!   close distance of last two points lead to problematic
!   temperature corrections from TCORR;
!   smoothing of velocity etc. at NS

! version : 10.4.4 April 2016: use XJ(CMF) for Thomson emissivity
!   instead of XJ(cont). Change restart accordingly.
!   all numbers refering to itlast=30 changed
!   to iitcmf_full.
!   variable RESTART_CMFALL introduced

! NOTE: version 10.4.0 to 10.4.4 only in this path, no to be confused
!   with standard_v10.4 (the latter should be consistent with v11.0)

! version: 11.0 July 2016: inclusion of all changes from standard_v10.4
!   from here on, standard_v10.4 and v11.0 should be
!   consistent if OPTCMF_FULL=.FALSE.

!   COPYOCC did not need to be rechanged, already OK in v10.4.4

!   inclusion of new parameterization of X-ray filling
!   factor, as resulting from the work by Koushik Sen/Jo
!   when combining the results from Feldmeier et al. 1997
!   and Owocki et al. 2013.
!   The present input allows to use the old parameterization
!   (rho^2 dependence of emissivity in radiative and
!   adiabatic shocks), and the new one, where the
!   emissivity in radiative shocks is rho-dependent.
!   When the parameter PX is set to -1000 or lower, the old
!   parameterization is used. If p > -1000, the new one is
!   used, with PX the exponent describing the distribution
!   of shocks (see Owocki+). For the old description,
!   FX corresponds to the volume-filling factor, for the
!   new one to the normalization of the shock-number, n_0.
!   The combined, radius-dependent filling factor is now
!   called FXL
!   Also updated: electron density for X-ray calculations;
!   since this density refers to the hot plasma, the
!   electron density has to refer to completely ionized
!   Helium (error in Feldmeier et al., but also used in
!   previous x-ray versions).
!   Also updated: the input value Rmin refers to units
!   of Rstar and not to the lowermost radius.

!   new formula 101 for CBB-transitions implemented,
!   to allow for an effective reading for tabulated collision-strengths

! version: 11.1 Aug 2016: updated routine CLUMPING from Jon
!   optically thick clumping generalized,
!   fvel normalized,
!   input generalized.
!   still missing: stop statement to prevent clumping below *sonic* point
!   consistent with standard v10.4.2

! version: 11.2 Sept 2016: changes regarding OP-data (ionization to
!   excited states) as in standard v10.4.3:
!   routine OP_RBFSETUP/CONVOL modified;
!   improved OP-cross sections when ionization to excited levels
!   STOP for RBF-formula 101 enforced
!   related changes in FRESFIN (INDFRQ1), NETMAT and OPACITC.
!   new routine INTERMEPHO1.
!   in brief: split of photo-integrals. For FRE ge FRECIN,
!   rates towards excited (and, if FRECIN = FRECFN, towards
!   ground state). Calculated in INTERMEPHO.
!   For FRECFN1 le FREC le FRECIN, corresponding photo-rates
!   towards ground-state. Calculated in INTERMEPHO1
!   Corresponds to explicit DR-approch,
!   rates are called RBF1 in megasalida.
!   Approach should be a good approximation. It assumes
!   that largest part of cross-sections (ge FRECIN) relates to
!   excited or ground-state, whilst lower frequency part
!   (if present)  corresponds to resonances and relates
!   to ground-state.
!  Somewhat wrong if ionization to higher than 2nd stage.

! version: 11.2.1 Oct. 2016: HeI hack no longer necessary, since FeIV lines
!   overlapping with HeI 584.33 manipulated in a similar way
!   as done by Paco (in LINE_LIST/cmf_all.f90).
!   Might need to be updated when additional problems with
!   HeI singlet lines. Thus far: OPT_HEI_SING -> .false.
!   New option OPT_DAMP_FG: damping of level corrections for
!   foreground elements only when this option is .true.
!   Default now (contrasted to previous versions) is .false.
!   (no damping)
!   Routine MODVL_UPDATE slightly modified
!   (check for min(chibar > min(sig_th) improved)

! version: 11.2.2 Oct. 2016: if OPTCMF_ALL, corrections of fore-ground
!   elements different from H/He always damped if MEANERR < -2.

! version: 11.2.3 Dec. 2016: new option OPT_MODIFY_SL_NOUP, to allow for
!   modified source-function of bg-lines with no upper level
!   that are overlapping with non-HHe elements.
!   Sobo transport (indexcmf/indexcmf1 = 1) for lines
!   from transitions outside ilow,imax

! version: 11.2.4 Dec. 2016: CBF rates for ionization towards excited
!   states improved. required updated princesa.f90,
!   consistent with standard 10.4.4

! version: 11.3.0 Feb. 2017: inclusion of Ng-extrapolation,
!   new option OPT_NG

! version: 11.3.1 March 2017: Inclusion at one additional freq. point
!   at HeII 303.797 (corresponding to A10HHe-energies),
!   to enforce identical solutions for different elements
!   (different freq. grids). Driven by a first discrepancy
!   when calculating HHe, HHeC and HHeN and HHeNC models
!   (in the HHeC models, HeII did not recombine for d40,
!   contrasted to the other models).
!   corresponds to v10.4.5 of standard-version
!   NOTE: change wavelength when HeII energies are changed
!   in data-file. Consistent with standard 10.4.5

!   Improved ALO for those cases were fine-structure
!   components from same transition overlap. Very
!   important for convergence!

! version: 11.3.2 July 2017: subr. OPACITC modified; calculate
!   only those bf-opacities below ILOW where ENIONND NE 0
!   (otherwise problems in RESTART for OPTXRAY (oxygen!)
!   Consistent with standard 10.4.6

! version: 11.3.3 Dec 2017: subr. ILOWIMAX_LTE only called when OPTCMF_ALL =
!  .FALSE.

! version: 11.3.4 March 2018: subr. CLUMPING -- previous bug identified
!   by Jon (1-fic^2) -> (1-fic)^2

! version: 11.3.5 April 2018: new variable for CONVERG_METALS:
!   error_met_mean (average error for all metals --occ1--
!   and all considered depth points).
!   MET_IMIN and MET_IMAX for restart (e.g., if OPTXRAY = T)
!   now transfered to subr. RESTART

! version: 11.3.6 July 2018: two additions (after call to ILOWIMAX_NLTE)
!   to allow for changes in ionization structure
!   (index-files for explicit elements need to be updated)


! version: 11.3.7 Feb. 2019: smaller changes, and introduction of
!   new approach for UP=0 transitions if the lower level
!   is not the ground or a metastable state (new option OPT_NEW_UPPER)

! version: 11.4.0 Oct. 2020: compatible with gfortran

! version: 11.4.1 Dec. 2020: in previous versions, approx. bg elements were
!   treated in LTE after T-update, disturbing the T-convergence.
!   (zig-zag in metal-convergence).
!   Now, we use the NLTE values from the previous NLTE iteration
!   for the outer part (the inner one is in LTE w.r.t. the new
!   T-structure), which gives a much smoother convergence.
!   Moreover, in this way we avoid problems in the CMF-transfer,
!   which might become (stronly) perturbed by the LTE source-functions,
!   initiating problems in the red part of the spectrum,
!   due to erroneous I-minus after the last strong line.
!   This, e.g., happened for model s10a because of an erroneous
!   AlIII resonance line, which was in extreme emission when using
!   the current LTE source function instead of (now) using the
!   NLTE occup. numbers from the previous iteration,
!   for the outer part.
!   Note: opacities for the CMF-transfer are calculated from
!   (i) occngold (NLTE, and .ne. zero only in between ilow, imax)
!       in the outer part for selected elements, and from
!  (ii) occng (LTE) in the inner part,
! (iii) whilst the approximate bg elements use the current
!   occng-values  (in LTE after T-update) for all depth points.
!   By overwriting with the NLTE occng-values from the outer part,
!   everything can converge much better.

! version: 11.4.2 Jan. 2021: additional changes (mostly in nlte_approx)
!   to improve the convergence.
!   If approximate approach (OPTCMF_FULL = .FALSE.),
!   delete specific files, which might have survived from
!   previous models, in particular CONT_FORMAL_CMF,
!   whose existence is checked in FORMALSOL

! version: 11.5.0 Sept. 2021: inclusion of updates made within 10.4.7, 10.5.1
!   and 10.6.0
!   for path with OPTCMF_FULL=.FALSE.

!   subr. MODVL_UPDATE modified; allow for
!   somewhat larger inconsistencies, but particularly
!   re-define xm0 exactly. Consistent with 10.4.7

!   small inconsistencies in WM-Basic Data Base
!   regarding DR data for OII (re-) detected. For those
!   transitions, stop statement in subr. intermedr de-activated
!   Consistent with 10.5.1

!   Change in FRESCAL AND FRESFIN, to allow for a consistent
!   frequency grid around Ly_gam and Ly_delta, in the same
!   spirit as previously for Ly_alpha and Ly_beta (important
!   for B-stars).
!   Specific problem in detail input detected:
!   for high frequencies (ga 100 A), at least RBF formula 12
!   might give erroneous values, when out of region where
!   fit is valid (detected in OIII data from Norbert/Miguel).
!   In future, range of validity needs to be provided.
!   Until then, RBF formula-12-cross-sections for lambda < 100
!   are approximated via
!   sigma(nu) = sigma_old(100A)*(nu(100A)/nu)^2
!   whenever sigma_old(nu) > sigma_old(100A)

!   Symmetric extension of Stark-broadening
!   to neigbouring "boxes" in sumopal (lte and nlte).
!   Addition (in FRESCAL) of freq. points for most important
!   Stark-broadened transitions (defined in sumopal_nlte)
!   saved in subr. SAVE_ADD_FREQ_POINTS (new file
!   ADD_FREQUENCIES), read in (or deleted, if starting
!   with iteration "0") after INDAT.DAT has been read.
!   FRESCAL updated (new variable NEWFREQ, which indicates
!   that freq. grid has changed even if IFRE did not change.
!   (can happen). Important to reset OP_FLAG, OP_FLAG2.
!   MEANERR <= -4.5 as additional convergence criterium

!   From here on, the OPTCMF_FULL = .FALSE. path
!   should be fully consistent with standard v10.6.0

!   This version required nlte_approx.f90 v8.4.0 and
!   dimes.f90 v1.3.2

! version: 11.5.1 March 2022: few changes in FRESCAL (definition of
!   freold, ifreold1), and increase of LWION1 to
!   if OPTDR=F, also fg-elements *without* DR for explicit
!   DR-transitions
!   tau_Ross approx. 10 (sometimes, 2 is not enough)

! version: 11.5.2 Nov 2022: few changes, to account for updates
!   from 10.6 to 10.6.1:
!   current hack, needs to be improved:
!   no longer stop, but only print if
!   ' RESTART, BUT FILE ADD_FREQUENCIES MISSING'

!   no longer stop, but only warning, when external
!   structure and
!   'xm not found at outer boundary of external struct. '

!   for further changes, see nlte_approx (sumopal)
!   and formalsol (Lemke broadening for Paschen + Brackett
!   lines until n=10)

! versions 11.6.0 to 11.6.2 do not exist. Next version is 11.6.3,
!   to make numbering consistent with v10.

! version: 11.6.3 Feb 2023: few changes in cmf_all/ray_complete.
!   (for u<0), to improve convergence.
!   ADDITIONAL CHANGES TO BE CONSISTENT WITH V10.6.3
!   (frequency grid: edges close to Lyman lines now allowed).
!   Multi-linear regression (including P^EXPO_P) for
!   calculating chibar_h in case linear regression is not
!   sufficient. The latter fit is only accepted if the
!   resulting exponents are physically meaningful
!   (EXPO < 0 AND EXPO_P > 0). If this condition is not met,
!   either linear fit results are used, or program stops.
!   Also: Additional checks in MODVL_UPDATE, whether
!   current radiative acceleration (for L > NDIV) is
!   locally beyond (then STOP) or close to Eddington limit.

!   Small changes in princesa, totout, and preformal
!   changes in formalsol/contin for optical thick
!   clumping -- remember: thus far allowed only in standard
!   approximation)

!   Still: no line blocking in sumopal_nlte/sumopal_lte
!   (nlte_approx.f90) for k=8 AND up=0
!   in "standard" path (needs to be checked)

! version: 11.6.4 May 2023: update of op_rbfsetup: couple of changes
!   to allow for additional rbf-files (Nahar etc.)
!   subr. smooth: since in new rbf files (e.g., CII)
!   the ratio of the maximum/minimum can be quite large,
!   the double precision fft used for the convolution
!   with a Gaussian suffers strongly from numerical
!   inaccuracies. Therefore now a test on this ratio:
!   if lower than 1.d8, then old approach, if larger,
!   then fft and convolution in quadruple precision.
!   in case, lower this limiting ratio

! version: 11.6.5 Sept 2023: Very old bug in subr. TRATRAD
!   cured (if TAUSF very low, expand exponential,
!   to prevent BCIC becoming zero!

!   Major change: Additional frequencies from fgrid (central
!   values) included, to allow for similar line-irradiation
!   when different foreground elements are used. To this
!   end, FGRID_LINES already called in subr. FRESCAL (once).
!   preparation for this routine (xion, vth(8)) in subr.
!   START.
!   Two interpolation schemes for continuum quantities
!   (J_nu etc.) now possible in subr. TRATCMF:
!   path=1 (old scheme, log log, suggested), and
!   path=1,2,3 (new scheme, constant values between
!   bg-frequency grid, not suggested, since large changes
!   compared to conventional approach)

!   almost perfect consistency with 10.6.5 if same NP

! version: 11.6.5.1 Nov 2023: now with modified cmf_all v1.4.6
!   (levels with up=0 illuminated by smoothed XJ)

! version: 11.6.5.2 Feb 2025: in subr. opacity_cmf, illumination
!   of levels with up=0 by pseudo-continuum (old approach).
!   Subr. ray_complete changed with respect to final check
!   of ALO (inclusion of ff).

! version: 11.6.5.3 March 2025: small changes in ILOWIMAX_LTE
!   for teff ge 15 kK, always HeI/II
!   for teff  < 15 kK, only HeI
!   MIGHT NEED TO BE ADAPTED

! version: 11.6.5.4 March 2025: small changes in calculation of T-struct
!   when dt(flux) < 0 and dt(th-balance) > 0, to avoid
!   artificial temperature minima
!   update of alo after being calculated in cmf_all
!   (after ray-complete and moments), to ensure a
!   (more or less) correct snew in calc_tlu_metals.
!   from now on, warning messages such as "problems with snew"
!   shall no longer apppear.

! version: 11.6.5.5 March 2025: small changes regarding gfortran,
!   and in nlte_approx.f90 to allow for higher vmic until
!   vmic=40 km/s. However, changes for the latter value
!   are severe, and effects need to be tested.
!   Thus far, we suggest a maximum value of 25!

! version: 11.6.5.6 May 2025: few bug fixes

! version: 11.7.0 May 2025: obsolete fortran features removed

! version: 11.7.1 June 2025: FRESCAL and FRESFIN slightly adapted
!   to allow for equidistant freq. steps behind
!   important edges for lambda < 227 A; only in
!   case of OPTCMF_FULL and OPTXRAY (see corresponding
!   routines. Equidist. steps are important for
!   a reasonable inter- and extrapolation of xj_cmf_coarse
!   in subr. interpol_cmf_to_coarse (package cmf_all.f90)
!   Currently, this works reasonably well for
!   wavblue = 200. (A) in module cmf_all_var (package cmf_all.f90)
!   but might need to be improved for lower values (but see
!   notes in that module),

!   couple of changes at blue-most frequ. in CMF treatment
!   (k=1), to obtain a consistent description
!   for J and H at the outer boundary. See comments
!   in cmf_all.f90, v 1.4.8

!   [slight modification of outer boundary condition (cmf_complete)
!   to allow for more accurate Eddington factors of the outermost
!   grid points. Tested with R(2) = 0.9 R(1).
!   In the default case of R(2) =0.98 R(1), only small differences.]
!   Though this modification works well in the O-star grid,
!   it doesn't work for strong B-star winds. Thus, re-modified

!   modification of met_imin and met_imax in nlte_approx.f90
!   (now version 8.5.7), subr. select and jbar_metals.
!   Now, if required, 4 ions in parallel

! version: 11.8.0 Sept...Dec 2025:
!   changes of subr. fgrid_lines in nlte_approx.f90
!   (details described there). Basically, the wavelength boxes
!   for calculating the averages opacities required for
!   line-blocking/blanketing are now the same for all models
!   (independent of Teff), avoiding potential discontinuities
!   when Teff is changed by a bit (consistent with FW 10.7.0

!   changes in subr. tcorr, to stabilize the lower boundary

!   new approach to handle negative intensities in
!   subr. ray_complete (cmf_all.f90), together with
!   consistent changes in moments. Now, the number
!   of inconsistent solutions for Jbar from moments
!   should considerably decrease, as well as the
!   "overall" changes in J_nu from iteration to iteration

!   ABUND_HAVE_CHANGED variable only checked after model-restart,
!   after first check reset to .false.

!   Still, certain transitions oscillate (e.g., NV1240, NIV4-29).
!   Carefully checked, no numerical, but physical problem.
!   Damping of foreground elements changed to geometrical
!   mean, helps a lot. Interestingly, shutting off the ALO
!   results in similar solutions

!   To ensure meaningful damping (see above), we now
!   artificially enforce ALMOST_CONVERGED = .true.
!   after iteration IIT=IIT_TEMP_CONV+20, where
!   IIT_TEMP_CONV is the iteration number of temperature
!   convergence. Moreover, currently we check (if
!   OPT_OIII_ITER = .true.) for a stabilization below 1.d-3
!   (before 5.d-4)

! version: 11.8.1 Dec 2025:
!   polished with nag-tools

! !! WARNING !!!!
! for future calculations, remember that opac_nolines is mean, not
! effective opacity
! !! WARNING !!!!

  Character (10) :: ver_nlte = '11.8.1'

End Module

!----------------------------------------------------------------------

Module fastwind_params

! important control parameters; change only if you know what you are doing
! (other control parameters included in specific models).

  Use :: nlte_type
  Use :: fund_const
  Implicit None

! I. Parameters in nlte.f90

! U and U2 functions, and much more
  Character (*), Parameter :: fpath_ufunc = '../inicalc/DATA/'

! specific data independent of used line list
Character*(*), Parameter :: fpath = '../inicalc/ATOMDAT_NEW/' 

!NOTE: atomic data for bg elements moved to specific sub-directories 
! standard data (old) with      4157124          83        7124 lines
!Character*(*), Parameter :: fpath1 = '../inicalc/ATOMDAT_NEW/ATOMDAT_NEW_old/' 

! or

! standard data (new) with      4157110          83        7110 lines
! repaired for N_VI, F_I, Cl_I, Al_III
Character*(*), Parameter :: fpath1 = '../inicalc/ATOMDAT_NEW/ATOMDAT_NEW_new/' 

!or

! standard data 2016 (new) with  6642198         132       42198 lines
! repaired for Cl_I, Al_III
! NOTE: UV-resonance line SIV at 1067 A erroneous!!!
!Character*(*), Parameter :: fpath1='../inicalc/ATOMDAT_NEW/ATOMDAT_2016_new/' 

! explicit RBF cross sections (Opacity Project Data etc.)
  Character (*), Parameter :: oppath = '../inicalc/OP_DATA_NEW/'

! data for xray treatment
  Character (*), Parameter :: fpath_xrays = '../inicalc/RaymondSmith/'

! controls velocity range where overlap is accounted for
  Real (dp), Parameter :: tauscrit = 0.01D0

! Limit for the Flux-Tcorr
  Real (dp), Parameter :: taurlim_default = .5D0 ! standard value, is adapted
! in case of thick winds

  Real (dp) :: taurlim !            actual value

  Real (dp) :: tmin1 = .4D0 !       minimum value for T(r) in units of Teff
! from comparison with Adi's models, a value of .4 to .5
! seems to be reasonable.

  Real (dp) :: tminabs = 6000. !    absolute minimum temp. which can be
! handeled

! range for X-ray luminosities (0.1, 0.15. 0.35 to 2.5 keV); Rosat and XMM
  Real (dp), Parameter :: const_ev = ev/(hh*clight)
  Real (dp), Parameter :: xred1 = const_ev*100.D0, xred2 = const_ev*150.D0, &
    xred3 = const_ev*350.D0, xblu = const_ev*2450.D0

  Real (dp), Parameter :: resol = 100. ! AVERAGE resolution of Gaussian, &
! corresponding to 3000 km/s (FWHM)
! for OP_RBFSETUP

  Real (dp), Parameter :: rsmpl_fac = 0.7D0, max_resol = 30000. ! for CONVOL

! --------------------
! params for nlte_approx and related

! ndebug= 1 minimum output (major ionization stages)
! ndebug= 2 including ionization fractions in subr. 'ionout', file 'METAL_IDL'
! generated
! ndebug= 3 including partition functions
  Integer (i4b), Parameter :: ndebug = 2

! allow or forbid inclusion of H for extended Starkwings in sumopal_lte and
! sumopal_nlte: default = F
  Logical, Parameter :: opt_extend_h_wings = .False.

! cutoff : cutoff energy (in units of ionization energy) for
! excited levels ('s' and 'm') to be included for opacities
  Real (sp), Parameter :: cutoff = 0.75 ! (1.0 = all)

! cutoff value to consider a transition as allowed (used in atomdat)
  Real (sp), Parameter :: gfcut_allowed = 1.D-4 ! to be on the save side

! minimum number fraction njk/ntot to be present if to be included in
! opacities
  Real (sp), Parameter :: fracmin = 1.D-12

! maximum level numbers for Adi's database (background)
  Integer (i4b), Parameter :: natom = 30, nion = 9, nlev = 50, nrec = 149, &
    levmax = 5000

  Character (2), Dimension (natom), Parameter :: names1 = (/ 'H ', 'HE', 'LI', &
    'BE', 'B ', 'C ', 'N ', 'O ', 'F ', 'NE', 'NA', 'MG', 'AL', 'SI', 'P ', &
    'S ', 'CL', 'AR', 'K ', 'CA', 'SC', 'TI', 'V ', 'CR', 'MN', 'FE', 'CO', &
    'NI', 'CU', 'ZN' /)

  Character (2), Dimension (natom), Parameter :: name = (/ ' H', 'He', 'Li', &
    'Be', ' B', ' C', ' N', ' O', ' F', 'Ne', 'Na', 'Mg', 'Al', 'Si', ' P', &
    ' S', 'Cl', 'Ar', ' K', 'Ca', 'Sc', 'Ti', ' V', 'Cr', 'Mn', 'Fe', 'Co', &
    'Ni', 'Cu', 'Zn' /)

  Integer (i4b), Parameter :: nbunch = 50000

  Integer (i4b), Parameter :: nion1 = nion - 1

! number of channels for nominal vturb = 10km/s and Teff = 30 kK
! Don't fiddle around with this number, if Stark-broadening should be used:
! Hydrogen Voigt parameters fitted for this value;
  Integer (i4b), Parameter :: sample_width = 120

! switch for radiative allowed/forbidden transitions, used to calculate
! collisional strengths (in sumopal_nlte, netmat_met)
  Real (dp), Parameter :: gfcut = 1.D-5

! minimum ionization fraction for met_imax (might need to be adapted)
! present value taken from zeta Pup model with X-rays (otherwise problems
! with Si at l=29 (used in subr. select)
  Real (dp), Parameter :: fracmin_imax = 1.D-6

! --------------------
! params for cmf_all and related

! adapt in case
! JO tested: for wavblue = 130 (NV->VI edge), no essential difference

! if wavblue set to lower values, then subr. frescal and fresfin should
! be modified to allow for equidistant steps. Otherwise, there is a certain
! zigzag around the true cmf-solution when J(CMF) is interpolated onto the
! coarse grid. This might be irrelevant though, since all J_nu integrals use
! trapezoidal rule, and the integrals should be OK.

  Real (dp) :: wavblue = 200.D0  ! this is not a parameter, since might be modified
! JO: at high temperatures, lower wavblue?
! (NIV reson. lines around 196 A, NV reson. lines around 140 A)
! seems to be OK; tested for D45 model, no changes

  Real (dp), Parameter :: wavred = 10000.D0
! for tests:  
! Real (dp), Parameter :: wavred = 300.D0
! "allowed" values 10000,15000,20000,25000,50000
! certain values (e.g., 45000) can lead to problems in xh_obs_coarse
! or even to stops (kcmf_start1 ne kcmf_start), because of
! unsuited frequency grid


! impact of gf_min tested; if set to 1.d-5, no change
! might be different if improved Fe-model will be used
  Real (dp), Parameter :: gf_min = 1.D-3 ! all lines with up=0 and gf < gf_min
! neglected
  Real (dp), Parameter :: vturbmin = 5.D5
  Real (dp), Parameter :: xmaxdop = 3.

! parameters for detecting significant line-overlap with explicit lines
  Real (dp), Parameter :: xm = 1.D0 ! overlap interval in units of vdop
  Real (dp), Parameter :: minrho = 0.3D0 ! minimum opacity contribution of
! considered explicit line

! if optout_xj_xh set to true, following files will be written:
! out_xj_cmf
! out_xj_cmf_coarse
! out_xj_app
! out_xh_obs
! out_xh_obs_coarse
! NOTE: xh_obs (fine and coarse) multiplied by (r(1)/Rstar)^2 = RMAX^2
! in contrast, XXH (from approx. solution) multiplied by (r(1)/SR)^2
  Logical, Parameter :: optout_xj_xh = .True.

! cut-off
  Real (dp), Parameter :: gfmin_ng = 1.D-3
! elements (explicit and background) where Ng-acceleration should be performed
! for tests with Ng also for HeII resonance line(s)
! integer(i4b), parameter :: no_ng_el = 4
! character*2, dimension(no_ng_el), parameter :: ng_el=(/'HE','C','N','O'/)
! integer(i4b), parameter :: no_ng_el = 1
! character*2, dimension(no_ng_el), parameter :: ng_el=(/'O'/)
  Integer (i4b), Parameter :: no_ng_el = 3
  Character (2), Dimension (no_ng_el), Parameter :: ng_el = (/ 'C', 'N', 'O' /)

! for testing overlaps (subr. overlap_detector)
  Real (dp), Parameter :: taur_test = 0.1
  Logical, Parameter :: test_overlap = .False.
! maximum separation of lines (in cm/s) for overlap detector
  Real (dp), Parameter :: deltavmax = 150.D5
  Real (dp), Parameter :: gf_cut = 1.D-3 ! not to be confused with gfcut!
  Real (dp), Parameter :: opal_cut = 0.1
  Real (dp), Parameter :: tau_cut = 0.1
  Real (dp), Parameter :: fac_vturb = 1.


  Integer (i4b), Parameter :: kmom_start = 2 ! first kmom_start frequ.
! point(s) only ray-by-ray

End Module

!----------------------------------------------------------------------

Module nlte_opt

! control options

  Implicit None

! to be consistent with older versions (before V9.6), set the next 5 options
! to .false.

! Enable full cmf-treatment, use cmf_all
  Logical, Parameter :: optcmf_full = .true.

! Enables update of photo-struct. by using flux-mean around iteration 10
! or higher (for OPTCMF_FULL)
! Only one update is performed

! default .true.
  Logical :: update_struct = .True.

! Exact treatment of photospheric lines from selected background elements
! default .true.
  Logical :: optphotlines = .True.

! Stark-broadening for (quasi-) resonance lines of metals and H at Lyman-jump
! default .true.  if approx. treatment
! default .false. if full cmf treatment
! LOGICAL :: OPTSTARK=.false.
  Logical :: optstark = .Not. optcmf_full

! controls whether DR (dielectronic recomb) for background elements enabled
! default .true.
! NOW (from 11.5.1 on): if OPTDR=F, not only bg-, but also fg-elements
! with explicit DR-data are treated without DR

  Logical, Parameter :: optdr = .True.

! ----------------------------------------------------------------------------
! the following options should be only changed if you know what you are doing

! new option
! if set to .true., the uppermost levels of bg-elements connected by
! specific transitions (see below) will be approximated in a different way
! than previously. In previous versions (and for OPT_NEW_UPPER = .FALSE.),
! ALL source functions of transitions with UP=0 were calculated in a TLA
! approach. If OPT_NEW_UPPER is set to .TRUE. the upper levels of transitions
! with UP=0 and lower levels NOT being the ground or a metastable state
! are approximated differently (UP = 0 and lower level = ground or
! metastable level remain in TLA): Here, we assume that the upper level is
! in LTE to the uppermost known level of the corresponding series, in other
! words, we assume that the departure coefficients remain constant.
! This affects particularly the Fe lines below 400 A, and the radiative
! acceleration, becoming lower in the transonic region.
! For a final solution, a complete Fe model atom (using Superlevels) is
! required
  Logical, Parameter :: opt_new_upper = .False.


! new option
! fast convergence if .false., but NIII emission lines will be wrong.
! If set to .true., more iterations until OIII level 10 (the pumped one)
! is (more or less) converged
  Logical, Parameter :: opt_oiii_iter = .True.

! new option
! if set to true, Ng extrapolation of source-functions for selected elements
! will be performed. For data and setup, see module ng_var in cmf_all.f90
  Logical, Parameter :: opt_ng = .True.

! new option, applies to OPTCMF_FULL
! if set to true, HeI singlet resonance line (584 A) will
! be treated as single line via CMF_SING (to avoid contamination be Fe/Ni
! lines)
! from v 11.2.1 default = .false., since corresponding changes in
! line_list/cmf_all.f90
  Logical, Parameter :: opt_hei_sing = .False.

! new option, applies to OPTCMF_FULL, default = .true.
! allows for modified source-function of bg-lines with no upper level
! that are overlapping with non-HHe elements
! change to .false. if original approach (two-level atom for all lines
! with no upper level) aimed at.
  Logical, Parameter :: opt_modify_sl_noup = .True.

! new parameter, from V6.1 on
! if set to .true., ILOWNL, IMAXNL will be constant throughout complete
! atmosphere, &
! controlled by the average in between NS-5 to NS+5
! if set to .false., ILOWNL, IMAXNL will be locally consistent (but with
! disadvantages regarding the CMF lines, see below)
! Thus far, GLOBAL must be set to .true. if OPTCMF_FULL
  Logical, Parameter :: global = .True.

! new parameter, from V11.2.1 on
! damping of level corrections for ALL foreground elements (default = F)
! if ALMOST_CONVERGED = .TRUE.
! Note: starting at version 11.2.2.,
! if OPTCMF_ALL, corrections of fore-ground elements different from H/He
! always damped if MEANERR < -2.
  Logical, Parameter :: opt_damp_fg = .False.

! ----------------------------------------------------------------------------
! treatment of photospheric lines via differential or integral method
! (integral method needs additional tables (k2atable, l2atable, l2table,
! l3atable) in DATA)
! default .true.
  Logical :: optphotlines_diff = .True.

! enhanced output for control of temperature correction (additional files!)
! default .false.
  Logical, Parameter :: outthb = .False.

! enhanced output for optically thick clumping (additional files!)
! default .false.
  Logical, Parameter :: outthick = .False.

! lot of enhanced output for opacity data cross-sections and convolution
! (additional files!)
! default .false.
  Logical, Parameter :: opdebug = .False.

! use Kurucz abundances for background elements
! default .false.
  Logical :: optkur = .False.

! simple temperature correction (out-dated). If enabled, check
! nlte_approx.f90,
! such that rates are not calculated (and added) twice -> subr. NETMET_MET
! default .false.
  Logical :: opttcor_simple = .False.

! test of photospheric clumping with rho_phot=rho_phot/fcl
! default .false.
  Logical :: opt_photclump = .False.

! to force single line treatment for all lines in (old) standard approach
! (v10.x), set OPTSING = .TRUE.
! default .true.
  Logical, Parameter :: optsing = .True.

! to forbid treatment of wind induced overlap (=>tests), set IOV_ONLY = .TRUE.
! default .false.
  Logical, Parameter :: iov_only = .False.

End Module

!----------------------------------------------------------------------

Module nlte_xrays

  Use :: nlte_type
  Use :: nlte_dim
  Use :: fund_const, Only: pi
  Use :: fastwind_params, Only: name => names1 ! NOTE: in fastwind_params,
! NAMES1 (here used as NAME) in upper case, NAME in mixed case
  Implicit None

! enable xray-treatment; controlled by content of INDAT.DAT
  Logical :: optxray
! enable Auger (K-shell) ionization;
!!!!! Auger L and M shell ionization for heavy elements still missing !!!!
  Logical :: optauger = .True.

  Real (dp), Allocatable, Dimension (:) :: enerx, fxl, alpha_fre
  Real (dp), Allocatable, Dimension (:, :) :: lambdax, lambdanu, opa_kshell

  Integer (i4b), Parameter :: n_kedges = 59 ! in ATOMDAT_NEW/k_shell_Auger_data
  Integer (i4b), Parameter :: k_nmin = 3, k_nmax = 8 ! min/max ion. stage to
! be included in freq. grid

  Real (dp) :: fx !                 VOLUME FILLING FACTOR or SHOCK NUMBER
! (corresponding to n_0)            FOR SHOCK EMISSION
  Real (dp) :: px !                 PARAMETER FOR dN/dlnr = n_0*(RMIN/r)^p
  Real (dp) :: rminx_eff !          MINIMUM EMISSION RADIUS (in units of SR)
! (slightly different from input RMINX (in units of Rstar))

  Integer (i4b) :: lxmin = 1 !      INDEX FOR MINIMUM EMISSION RADIUS

  Real (dp), Dimension (n_kedges) :: eth, sigma, s, zeff

  Real (dp), Dimension (n_kedges, 6) :: aug_1_6

  Integer (i4b), Dimension (n_kedges) :: z, n

! ---------------------------
! SPECIFIC PARAMETERS, MOSTLY FOLLOWING FELDMEIER+ 1997

! BEGIN OF X-RAY TREATMENT (RED SIDE, IN A)
  Real (dp), Parameter :: lamred = 3500.

! CONSTANTS FOR FREQ. INTEGRATED COOLING FUNCTION, FELDMEIER EQ. 7
  Real (dp), Parameter :: ar = 1.64D-19
  Real (dp), Parameter :: fac = (93.D0*sqrt(3.D0)-40.D0*pi)/10.D0

! CONSTANTS FOR ADIABATIC PART
  Real (dp), Parameter :: mfn = -4.D0/9.D0

! CONSTANTS FOR RADIATIVE PART
  Real (dp), Parameter :: a = 0.8722544D0 ! SEE FELDMEIER EQ. 11

! for info and checks
! JO Oct. 2025: now USEd from fastwind_params
! CHARACTER*2, DIMENSION(30) :: NAME
! DATA NAME /'H ','HE','LI','BE','B ','C ', &
! &          'N ','O ','F ','NE','NA','MG', &
! &          'AL','SI','P ','S ','CL','AR',&
! &          'K ','CA','SC','TI','V ','CR',&
! &          'MN','FE','CO','NI','CU','ZN'/

End Module

!----------------------------------------------------------------------

Module run_once

! global options

  Implicit None
  Logical :: start_weigh11 = .True.
  Logical :: start_cmfsing = .True.
  Logical :: start_cmfmulti = .True.
  Logical :: start_jbar_metals = .True.
  Logical :: start_cmf_simple = .True.

End Module

!----------------------------------------------------------------------

Module nlte_var

  Use :: nlte_type
  Use :: nlte_dim
  Implicit None
! ----------------------------------------------------------------------
! wavelength (A) for Hydrogen Lyman-jump (red side)

  Real (dp) :: lam_ly
! ----------------------------------------------------------------------
! transition point for selected elements

  Integer (i4b) :: lwion1

! ----------------------------------------------------------------------
! correction for lowermost input flux dB/dtau

  Real (dp) :: corrfc

! ----------------------------------------------------------------------
! various velocities

  Real (dp) :: vmax, vsound, vdiv
  Integer (i4b) :: ndiv = 0, ndiv_calc ! (OUTERMOST POINT OF QUASI-HYDROSTATIC
! REGIME

! ----------------------------------------------------------------------
! variables for clumping parameters

  Integer (i4b) :: npara_clf = 11 ! maximum number of parameters,
! will be overwritten
  Real (dp), Dimension (11) :: para_clf_input = -1.D99
  Real (dp), Allocatable, Dimension (:) :: clf_field

! ----------------------------------------------------------------------
! variables for selected Stark-broadened transitions in subr. sumopal

  Integer (i4b) :: nstark
  Integer (i4b), Parameter :: nstarkmax = 500
  Integer (i4b), Dimension (nstarkmax) :: idstark
  Real (dp), Dimension (nstarkmax) :: wavestark

! ----------------------------------------------------------------------
! variables for ff gaunt-factors

! local parameter IONMAX

  Integer (i4b), Parameter :: ionmax = 9
  Integer (i4b) :: izmin, izmax

  Real (dp), Dimension (ionmax, id_ndept, id_frec1) :: gffnlte

! ----------------------------------------------------------------------
! alphafs


  Real (dp), Dimension (id_rbftr, id_frec1) :: alphafs

! ----------------------------------------------------------------------
! blocf

  Real (dp), Dimension (id_frec1) :: as, a2s, a3s, bs, b2s, b3s
! ----------------------------------------------------------------------

! ccom1 (without alevel)

  Real (dp), Dimension (id_ndept) :: xnelte, xne_ratio
  Integer (i4b), Dimension (id_atoms, id_ndept) :: ionis
  Integer (i4b), Dimension (id_atoms) :: ionisst
! ----------------------------------------------------------------------
! ccom4 dimension corrected (from kis to kis+1)

  Real (dp), Dimension (id_atoms, id_kisat+1) :: xl
  Real (dp), Dimension (id_atoms, id_kisat+1, id_ndept) :: partit
! ----------------------------------------------------------------------
! ccomoc

  Real (dp), Dimension (id_atoms, id_kisat+1, id_ndept) :: occnum
  Integer (i4b), Dimension (id_ndept, id_atoms) :: mainion
! ----------------------------------------------------------------------
! cffine

  Real (dp), Dimension (id_frec2) :: fref, wfref
  Integer (i4b), Dimension (id_frec2) :: indfin
  Integer (i4b), Dimension (id_rbftr) :: indfrq, indfrq1
  Integer (i4b) :: ifref
! ----------------------------------------------------------------------
! cfgrid

  Real (dp), Dimension (id_frec1) :: fre, wfre, hc, r1k
  Real (dp), Dimension (id_frec1) :: freold
  Integer (i4b), Dimension (id_frec1) :: lthin
  Integer (i4b) :: ifre = 0 !       initialize
  Integer (i4b) :: ifreold1
! ----------------------------------------------------------------------
! new

  Real (dp), Dimension (:), Allocatable :: fre_lines
  Integer (i4b) :: ifre_lines
! ----------------------------------------------------------------------
! clevel

  Real (dp), Dimension (id_llevs) :: alevel, blevel, nist, xnkk
  Real (dp), Dimension (id_atoms, id_kisat+1, id_ndept) :: enionnd, &
    enionnd_lte
! ----------------------------------------------------------------------
! cmms

  Integer (i4b), Dimension (id_rcbbt) :: mmlow, mmup
  Integer (i4b), Dimension (id_rbftr) :: mlow, mup
  Integer (i4b), Dimension (id_cbsft) :: mclow, mcup
! ----------------------------------------------------------------------
! comcmf1,comcmf2,comcfm3,comcmfs

  Integer (i4b), Dimension (id_nttrd) :: indxlamc
  Character (6), Dimension (id_nttrd) :: lablinelo, lablineup

  Real (dp), Dimension (id_nttrd) :: xlamcmf, vdopcmf, blucmf, bulcmf, aulcmf, &
    xmaxdop
  Real (dp), Dimension (id_ndept, id_npoin) :: tauz1, tau1, pp1 = 0.D0
  Real (dp), Dimension (id_nfcmf, id_ndept) :: ak1, az1, slinek
! ----------------------------------------------------------------------
! comcmfa,comcmfb,comtemp

  Logical :: optneupdate, optcmf, optmodel, optlucy, optmet, &
    optlines = .False., lines, lines_in_model, metals_converged = .False., &
    metconv2 = .False., lte_update = .False., concon = .False., &
    optcmf_all = .False., almost_converged = .False., no_alo = .False., &
    no_alo_metals = .False., no_check = .False., after_update = .False., &
    restart_cmfall = .False.

  Real (dp) :: optmixed

! ----------------------------------------------------------------------
! comcont

  Real (dp), Dimension (id_ndept, id_frec1) :: opac, strue, xj, xxk, alo, &
    strue_m_old, strue_m_new, opat_m_old, opat_m_new, opac_nolines, &
    opat_m_nolines, etat_nolines, etat_m_nolines, thomson_lines
! , RAYLEIGH = 0.D0
  Real (dp), Dimension (:, :), Allocatable :: xj_save, xj_cmf_coarse ! also
! used as
! XJ_SMOOTHED
! in
! CMF_COMPLETE
  Real (dp), Dimension (:), Allocatable :: xh_obs_coarse
  Real (dp), Dimension (id_ndept, id_frec1) :: tradj, tradj1
  Real (dp), Dimension (id_ndept, id_frec1) :: xxh
  Real (dp), Dimension (id_ndept-2) :: qq1, qq2

  Integer (i4b) :: kcmf_start, kcmf_end

! ----------------------------------------------------------------------

! comiact
  Integer (i4b), Dimension (id_atoms) :: imia, imaa, imianl, imaanl
! ----------------------------------------------------------------------

! comopac
  Real (dp), Dimension (id_nttrd, id_ndept) :: opaclin
! ----------------------------------------------------------------------
! compz

  Real (dp), Dimension (id_npoin) :: p
  Real (dp), Dimension (id_ndept, id_npoin) :: z
! ----------------------------------------------------------------------
! comtlucy

  Integer (i4b) :: nsdiv
  Real (dp) :: sr, srvmin, srnom, constt, rtau23
! ----------------------------------------------------------------------
! comtrad, comsa

  Integer (i4b), Dimension (id_rcbbt) :: indexrbb
  Integer (i4b), Dimension (id_nttrd) :: indexcmf, indexsa
  Real (dp), Dimension (id_nttrd, id_ndept) :: tlumat, tulmat, opacm, scontm, &
    betold
! ----------------------------------------------------------------------
! comvelo

  Real (dp) :: h1
! ----------------------------------------------------------------------
! comweigh

  Real (dp) :: w0last
  Real (dp), Dimension (id_ndept, id_npoin-1) :: vp, vp2
! Added for line-force, JS and moments equations (JP)
  Real (dp), Dimension (id_ndept+1, id_npoin-1) :: vp1, vp3
! ----------------------------------------------------------------------
! der,der1

  Integer (i4b) :: idum, itmin = 1
  Real (dp) :: ggrav, teff, sigem, ckappa, xt, a2te
  Real (dp), Dimension (id_ndept) :: delsig = 1.D0, dqdtau = 0.D0, &
    delmu = 1.D0, delpara = 1.D0

! ----------------------------------------------------------------------
! cqual

  Integer (i4b), Dimension (id_llevs, id_ndept) :: iqual
! ----------------------------------------------------------------------
! hopf,mic_turb,precision

  Real (dp) :: qinf, q0, gamhopf, vturb, precis
! ----------------------------------------------------------------------
! ithist

  Real (dp), Dimension (30) :: exith, xkapith, corrith
! ----------------------------------------------------------------------
! modnam

  Character (50) :: modnam
! ----------------------------------------------------------------------
! otra

  Integer (i4b) :: indr, indc
! ----------------------------------------------------------------------
! output

  Integer (i4b) :: iconver
! ----------------------------------------------------------------------
! path (ode)

  Integer (i4b) :: kmax, kount
  Real (dp) :: dxsav, xp(200), yp(10, 200)
! ----------------------------------------------------------------------
! soluna
  Logical :: unasol = .False., megas
! ----------------------------------------------------------------------
! store

  Real (dp), Dimension (id_frec2) :: bn
  Integer (i4b), Dimension (id_frec2) :: index
! ----------------------------------------------------------------------
! ufuncg,ufunck

  Real (dp) :: ufung(9, 17), ufunk(21, 76)
! ----------------------------------------------------------------------

! .. parameters ..
  Integer (i4b), Parameter :: n1_b = 21, n2_b = 6
  Integer (i4b), Parameter :: n1_l1 = 61, n2_l1 = 71
  Integer (i4b), Parameter :: n1_l2 = 12, n2_l2 = 18
  Integer (i4b), Parameter :: n1_s = 192, n2_s = 200
! ..
! U2-values from tables
  Real (dp), Dimension (n1_b, n2_b) :: u2b
  Real (dp), Dimension (n1_l1, n2_l1) :: u2l1
  Real (dp), Dimension (n1_l2, n2_l2) :: u2l2
  Real (dp), Dimension (n1_s, n2_s) :: u2s
! ----------------------------------------------------------------------
! vthermal, abundnd, comyhe

  Real (dp) :: yhein
  Real (dp), Dimension (id_atoms) :: vther
  Real (dp), Dimension (id_atoms, id_ndept) :: abundnd

! prespecified elements
  Integer (i4b) :: nat_tot, nat_bg
  Real (dp), Dimension (30) :: abundbg
  Real (dp), Dimension (31) :: abchanged ! might include "FX"
  Character (6), Dimension (30) :: names_bg
  Character (6), Dimension (31) :: names_changed  ! might include "XRAYS"

  Logical :: abund_have_changed = .False.
! ----------------------------------------------------------------------
! weights

  Real (dp) :: w5(5), w9(9)
! ----------------------------------------------------------------------
  Real (dp), Dimension (id_rbftr, id_frec1) :: op_data
! FRECFN1 analogous to FRECFN (ground-state ionization), but for OP-data
! (rbf-formula 20)
! only until minimum energy where sigma ne 0 (in between FRECIN AND FRECFN)
! if sigma ne 0 for energies lower than FRECFN, then FRECFN1=FRECFN
! calculated in routine RBF_SETUP. This is a certain approximation
! (sometimes, sig ne 0 even below ground-state edge), but used
! for consistency with the remaining approach.
  Real (dp), Dimension (id_rbftr) :: eopdiff, frecfn1
  Logical, Dimension (id_rbftr) :: op_flag, op_flag2

! ----------------------------------------------------------------------

End Module

!***********************************************************************

Module cmf_multi

  Use :: nlte_type
  Use :: nlte_dim
  Implicit None

  Type :: infotable
    Integer (i4b) :: index
    Integer (i4b) :: indbw
    Real (dp) :: xov
    Real (dp) :: if
  End Type

  Integer (i4b), Parameter :: lto = 2*id_ndept - 1

  Integer (i4b), Dimension (id_nttrd) :: indov

  Integer (i4b), Dimension (:), Allocatable :: indoverlap_order

  Type (infotable), Dimension (:), Allocatable :: infotab

  Real (dp), Dimension (:, :, :), Allocatable :: usave, vsave
  Real (dp), Dimension (:, :), Allocatable :: iplus, iminus

  Integer (i4b), Dimension (id_npoin-1) :: lmaxjp
  Real (dp), Dimension (id_npoin-1) :: fmin, fmax
  Real (dp), Dimension (id_npoin-1) :: iminus_cont
  Real (dp), Dimension (id_npoin-id_ndept+1) :: iplus_cont

  Real (dp), Dimension (lto, id_npoin-1) :: uv, zray

  Logical :: cmf_all_bluwi = .False.

End Module

!***********************************************************************

Module tcorr_var

  Use :: nlte_type
  Use :: nlte_dim
  Implicit None

! Limit for the Flux-Tcorr
  Real (dp) :: taurlim ! actual value (overrules TAURLIM_DEFAULT)

  Real (dp) :: emaxtc

! new variable: convergence control of metals
  Real (dp) :: error_met, error_met_mean
! .
! Free-free quantities
  Real (dp), Dimension (id_ndept, id_frec1) :: ffopa, dtffopa, ffopa_m = 0., &
    dtffopa_m = 0.

  Real (dp), Dimension (id_ndept) :: qffh, qffc, dtqffh, dtqffc
! .
! Continuum quantities
  Real (dp), Dimension (id_ndept) :: qcbfr, qcbfi, qrbfr, qrbfi, dtqrbfr, &
    dtqrbfi, dtqcbfh, dtqcbfc

  Real (dp), Dimension (id_ndept) :: qcbbu, qcbbd, dtqcbbu, dtqcbbd

! Quantities for metal-background (no coll. ionization, since not present
! in our approx.)
  Real (dp), Dimension (id_ndept) :: qcbfr_m = 0., qcbfi_m = 0., &
    dtqcbfr_m = 0., dtqcbfi_m = 0., qrbfr_m = 0., qrbfi_m = 0., &
    dtqrbfr_m = 0., dtqrbfi_m = 0., qcbbu_m = 0., qcbbd_m = 0., &
    dtqcbbu_m = 0., dtqcbbd_m = 0.


! Logical controls
  Logical :: enatcor, opttcor, temp_converged = .False.
! .

End Module

!***********************************************************************

Module photstruc

  Use :: nlte_type
  Use :: nlte_dim
  Implicit None


  Integer (i4b) :: ndmod
  Real (dp) :: xt, const, expo, expo_p, sigemin, clf_const = 1.
  Real (dp), Dimension (:), Allocatable :: xmext, text
  Real (dp), Dimension (:), Allocatable :: pkur, xnekur, gradkur
  Real (dp), Dimension (:), Allocatable :: rhotlu, xnetlu
  Real (dp), Dimension (:), Allocatable :: rmich, rmich1, vmich, rhomich
  Real (dp), Dimension (id_ndept) :: xm, chibar_h, delchi
  Character (6) :: modtype
End Module

!***********************************************************************

Module nlte_porvor

  Use :: nlte_type
  Use :: nlte_dim, Only: id_ndept, id_frec1
  Implicit None

  Real (dp), Parameter :: epsi = 1.D-12 ! FOR PRECISION OF OPA_EFF_RAT IN
! OPTICALLY THIN APPROACH
  Real (dp), Parameter :: epsi1 = 1.D0 + 1.D-12 ! FOR PRECISON OF OPA_EFF_RAT

  Logical :: optthick

  Real (dp) :: w_bg_red, w_bg_blue ! will be set to WAVBLUE and WAVRED in
! nlte_approx
  Real (dp) :: conv_eff_opa_rat = 1.0D0

  Real (dp), Allocatable, Dimension (:) :: fic_field, fvel_field, hpor_field

  Real (dp), Dimension (id_ndept) :: fic, tcl_fac_line, tcl_fac_cont, fvel, &
    fvol, hpor

  Real (dp), Dimension (id_ndept, id_frec1) :: opa_eff_rat = 1.D0, &
    opa_eff_rat_old = 1.D0
! these are the crucial quantities which correct the mean opacities
! <opa> = opa/clf to effecitve opacities opa_eff = <opa> * opa_eff_rat with
! opa_eff_rat = (1 + fic * tau_cl)/(1 + tau_cl)
! in the optically thin case, opa_eff_rat has to be unity at all frequency
! points, by forcing tau_cl to 0.

End Module

!***********************************************************************

Module nlte_wavcon

! new (Sept 2023) module for transporting WAVBLUE, WAVRED, WAVCON
! inside nlte; transfer from nlte_approx via subr. FRESCAL_EXPL_LINES

  Use :: nlte_type
  Implicit None

  Real (dp) :: wwavblue, wwavred, wwavcon

End Module

!***********************************************************************

!main program

!***********************************************************************

Program fastwind

! -----this program will - hopefully - solve the nlte line-formation
! -----problem for unified atmospheres in an optimum  way

! ----- V11 : complete cmf-solution, for explicit and background (bg) elements,
!             use OPTCMF_FULL = T
! until V10 : cmf-solution for explicit elements, approximate method
!             for bg ones: use OPTCMF_FULL = F
  
! program can currently be used for calculation of hydrogen and helium,
! silicon, magnesium, C, N, O (DETAIL versions)
! still, certain improvements might be required

! and for all elements with Adi's data (if you have a corresponding
! input file). Note, however, that using elements with a large
! number of levels (and line number) is inconsistent with present assumptions.
! At some point, super-level approach required!
! Currently: simplified packing, assuming that for sublevels n_i/g_i = const.
! in super-level approach, the latter estimate needs to be replace by
! LTE conditions. Implies that current approach is correct as long as sublevels
! are separated by h delta nu/kT << 1  
  
! developing/programming: j.puls, e. santolaya-rey (version 0.0), 
! miguel a. urbaneja (tcorr and more ions, version 7.2.1) and
! few others, incl. Jorge Rivero Gonzalez.

! The user can include her/his own clumping laws by manipulating subroutine
! clumping (max. number of parameters: 10), Thus far, optically thick clumping
! in the approx. of Sundqvist+ 2018 can only be treated in the v10
! (optcmf_full = F) approach. In the v11 approach, only optically thin clumping
! can be considered  

! note: if new atoms are implemented, provide new formulae in
! subroutines crossbf (photo cross sections),
! intermecol (coll. cross sections (ionisation))
! and tratcol (coll. excitation)
! if necessary!!!
! also required is update of ilow, imax in routine ilowimax
! (_old, _lte, _nlte (latter in nlte_approx.f90))

! missing points (ONLY when OPTCMF_FULL = .FALSE.)

! i)calculating taur(nlte)
! ii)updating of hydro-structure with respect to tau'(nlte)
! iii)(overlapping) lines included in continuum and ionization integrals
! maybe important at edges
! iv) external photospheric structures only for unclumped models (generally)
! v) RBF-formula 101 (currently not allowed) needs to be updated and checked.

! -----------------------------------------------------------------------

! IF OPTCMF_ALL = .TRUE., (which happens for OPTCMF_FULL = T with begin
! of iteration 30), then ...

! 1.
! ... XJ corresponds to XJ_CMF_COARSE
! (thus far, just from end of CONT until begin of
! CMF_COMPLETE, XJ is 'standard' approximate obs. frame quantity)

! In the former case (complete cmf-transfer), the qualifyer XJ(1,IFRE+1) #
! is set to 2 (otherwise -- approximate treatment, it's set to 1).
! But note that only those frequencies that are treated by CMF_COMPLETE
! provide actual (coarse grid) CMF intensities.
! Outside this range (< wavblue or > wavred), the old (approx. values)
! are present, and used for those quantities (ionization integrals,
! Sobo-beta's, dielectronic rates, line rates ...) which are needed in the
! rate equations (both for explicit and background elements).
! In this case, the frequencies and mean intensities are observer's frame
! quantities (whilst the rest are cmf-quantities)!!!

! 2.
! ... for XXK there is the same situation, though XXK is approximated
! by assuming the same eddington-factor as in the approx. treatment

! 3.
! ... the corresponding "continuum" ALO (in the CMF_COMPLETE range,
! when calculating the photo-integrals) has been set to zero.
! The line-ALOs are (almost) always considered
  

! 4.
! ... all important rates are calculated in the CMF, since the
! frequencies and mean intensities are CMF quantities. Only outside the
! CMF_COMPLETE range, the rates are (somewhat inconsistently) calculated
! in the obs. frame, and the frequencies are also interpreted as being
! obs. frame frequencies. In particular, this refers to Sobo-lines in
! subr. TRATRAD that uses generalized limb-darkening coefficients
! from FACTINC (obs. frame radiation field stored in CONT1, and converted
! in CONVERSION). Note that for OPTCMF_ALL = .TRUE. Sobolev line transfer
! inside the CMF_COMPLETE range is forbidden (checked in TRATRAD).

! 5.
! TAUR is read from TAU_ROS after each restart
! -----------------------------------------------------------------------

! XNH is hydrogen density (total one, HI+HII), independent of LTE/NLTE

! -----this program uses three values for the stellar radius:
! --   1) the input value rstar, which is the nominal radius where
! --      tau(-lucy) = 2/3; note that log g  and teff are
! --      also defined with respect to this radius. srnom = rstar * rsun.
! --   2) the lowest value of r(nd) = sr < srnom which is the scaling
! --      radius of the poblem.
! --   3) the radius where the beta-field and the photosphere merge
! --      ("transition point"): srvmin
! --
! --   rtau23 is the ratio of srnom/sr
! --

! -- Basic scheme of calculation:
! In tlucy, calculate r23= r(taup=2/3). From this value, and using
! r(ns) (from modvl) as well as srvmin/srnom = r(ns)/r23,
! srvmin = = srnom*r(ns)/r23can be calculated.

! This latter value is iterated in parallel with modvl:
! input srvmin -> r(ns) etc.
! This defines finally srvmin and rtau23= srnom/sr
! Note: all this relies on r23, which is calculated from taup = 2/3
! in the hydro iteration scheme. Since taup bases on LTE opacities,
! the position of the actual value of taup=2/3 (defining srnom) is
! somewhat approximative (but see next paragraph).
! For a check, subsequent values T(taur=2/3) are printed out.

! IF the photospheric structure is NOT updated, srvmin and rtau23 will remain
! constant throughout all NLTE iterations.
! If the structure IS updated, then the division point srvmin
! (in absolute values) will remain unchanged (assuming that taup has not
! changed much in the wind), and a new SR can be calculated using the
! updated grad-values. From this quantity, rtau23=srnom/sr(new)
! is updated as well

! WARNING!!! program nlte uses

! abund = n(k)/n(H), i.e. normalized to hyd.

! program nlte_approx uses
! abund = n(k)/Sum(n(k)), i.e. normalized to ntot

! Thus, also the value of summass = Sum(abund(i)) in both cases differ!!!

! END WARNING!!!

! !! WARNING !!!!
! for future calculations, remember that opac_nolines is mean, not
! effective opacity
! !! END WARNING !!!!

!IMPORTANT THINGS TO DO
!data: check completeness of Fe (background), and include missing data
! (e.g., FeIV 4f series and most FeII lines in the optical
!code: implement superlevel approach, level dissolution, and optically thick
! clumping for v11

! LU (FILE) numbers (except no 1 to 7, which are used only locally

! *** in nlte.f90
! 8 CONVERG
! 9 CONVERG_METALS
! 10 MAXTCORR.DAT
! 11 MAX_CHANGE
! 13 CLUMPING_OUTPUT (used only locally)
! 14 XLUM_ITERATION

! 15 LTE_POP
! 16 NLTE_POP
! 17 BLEVEL (scratch)
! 18 OCCNO

! 20 MODEL
! 21 TEMP
! 23 FLUXCONT (used only locally)

! UNIT 30-33 files for ulecture (used only locally)

! UNIT 40-43 (DIRECT ACCESS, scratch):
! J,H,K,N (cont), function of L, index=energy/lambda
! UNIT 50-53 (DIRECT ACCESS,scratch):
! J,H,K,N (cont), function of energy/lambda, index=L

! 60 CHIBAR_H
! 61 GRAD.OUT, GRAD_START.OUT, and CHIBAR_H_CMF (both used only locally)

! 70, 71: TABLES FOR ULECTURE (used only locally)

! 77 cross-sections for bf and related quanties 'opori'
! 771 as 77, used only locally

! 100 megasalida

! 999 FASTWIND.LOG

! *** in nlte_approx.f90

! 111 nl3i_all
! 112 nl3a_all
! 113 nl3g_all
! 114 nl3info_all

! 200 OUT.IDL
! 201 BG_LINES.CMF

! 211 generalinfo, atomgf3
! 212 partit, atomcol
! 213 atomnl2_meta_v2.0
! 214 atomnldr

! 287 METAL_IDL_LTE/METAL_IDL
! 287 METAL_IONS_LTE/METAL_IONS

! 291-294: kl_tables

! *** in cmf_all.f90

! 337 linetab = LINES_xxx.dat in inicalc/DATA

! *** no important files used in lambda_xray.f90

  Use :: version_nlte
  Use :: nlte_type
  Use :: nlte_dim
  Use :: fund_const, Only: amh, akb
  Use :: fastwind_params, Only: taurlim_default

  Use :: run_once !                 all start variables

  Use :: nlte_opt, Only: optcmf_full, update_struct, opt_oiii_iter, outthb, &
    opt_photclump, outthick, global, opt_ng, optstark

  Use :: nlte_xrays, Only: optxray, fx

  Use :: princesa_var, Only: labl, nat

  Use :: nlte_var, Only: enionnd, optneupdate, optcmf, optmet, optmixed, &
    optmodel, optlucy, lines, lines_in_model, metals_converged, metconv2, &
    optcmf_all, almost_converged, no_alo, no_alo_metals, no_check, &
    after_update, iqual, qinf, q0, gamhopf, vturb, precis, xkapith, exith, &
    corrith, modnam, iconver, unasol, megas, lte_update, nat_tot, &
    names_changed, abchanged, npara_clf, para_clf_input, clf_field, imia, &
    imaa, imianl, imaanl, ionis, concon, lwion1, ndiv, ndiv_calc, fre, ifre, &
    xne_ratio, restart_cmfall, nstarkmax, nstark, wavestark

  Use :: nlte_var, Only: ifre, xj

  Use :: nlte_var, Only: corrfc, sr, vmax, vsound, vdiv, rtau23

  Use :: nlte_var, Only: xnelte_common => xnelte

  Use :: tcorr_var, Only: enatcor, opttcor, taurlim, temp_converged, emaxtc, &
    error_met, error_met_mean

  Use :: photstruc, Only: xt

  Use :: nlte_porvor, Only: optthick, fic, tcl_fac_line, tcl_fac_cont, fvel, &
    fvol, hpor, fic_field, fvel_field, hpor_field, opa_eff_rat, &
    opa_eff_rat_old, conv_eff_opa_rat

  Implicit None

! .. parameters ..
  Integer (i4b), Parameter :: nd1 = id_ndept, np1 = id_npoin
  Integer (i4b), Parameter :: kel = id_atoms ! , kis = id_kisat
! ..
! .. local scalars ..

  Character (256) :: clf_para_str
  Character (10) :: ver_nlte_app, ver_cmf_all
  Character (3) :: ncname
  Real (dp) :: beta, corvm, emax, exkap, ggrav, hei, relex, relhyd, relkap, &
    rmax, srnom, teff, tmin, vmin, vmin1, xmet, xkap, xmloss, xmmax, xmu, yhe, &
    yhein, meanerr

  Real (dp) :: gamx, mx, rminx, uinfx, px

  Real (dp) :: dtmean

! not used:  Integer (i4b) :: ngstart

  Integer (i4b) :: i, iend, iitcmf, iitsobo, iitupdate, in, iq, istart, ithyd, &
    ithyds, itlast, iit, itmore, k, j, l, nrec, ntemp, ntemp1, nd, np, ittcor, &
    ncor, set_step, set_first, lwion, iit_temp, iitcmf_full, opt_ray_counter, &
    delta_iit, iit_start_cmf, iit_updated, iit_temp_conv

  Logical :: accel, fastsobo, frenew, grey, helone, hopfself, optcmf1, optcv, &
    optlte, optmod, rybicki, expansion, restart_met_lte, lastlte, &
    first_metals, updated, update_done, exi, start_xray, first_rateli, &
    update_taur, flag_ng, flag_ng_met, opt_ray, first_cmf
! ..
! .. local arrays ..
  Real (dp) :: dilfac(nd1), dilfac1(nd1), dvdr(nd1), err(nd1), errold(nd1), &
    r(nd1), rho(nd1), pressure(nd1), rtemp(nd1), specmat(5, nd1), taue(nd1), &
    taur(nd1), temp(nd1), velo(nd1), xne(nd1), xnh(nd1), dtfcorr(nd1), &
    fluxerr(nd1), clfac(nd1), xnelte(nd1), xnelte_taur(nd1), taur_taur(nd1), &
    tsho(nd1), ferr_cmf(nd1), dtfcorr_cmf(nd1)

  Real (dp) :: co_dum, time_cpu

  Integer (i4b) :: ilow(nd1, kel), imax(nd1, kel), ilownl(nd1, kel), &
    imaxnl(nd1, kel), index(nd1), opt_ray_arr1(4), opt_ray_arr2(4), &
    optrayarr1, optrayarr2

! ..
! .. external subroutines ..
  External :: checkmass, cont, copyocc, formal, fresfin, hopfpara, lte, model, &
    opacitc, pzgrid, rateeq, rateli, start, ulecture
! ..
! .. intrinsic functions ..
! INTRINSIC ABS,MAX,MIN
! ..
! .. data statements ..
  Data grey/.True./
  Data optlte/.True./
  Data first_metals/.True./
  Data start_xray/.True./
  Data first_rateli/.True./
! currently not used:  Data before_hydro_update/.True./
  Data iit_updated/0/
  Data update_done/.False./

! opt_ray before hyd. update
  opt_ray_arr1 = (/ 0, 1, 3, 7 /)
! opt_ray after hyd. update
  opt_ray_arr2 = (/ 0, 1, 3, 5 /)

! to set up model
  nrec = id_llevs
  nd = nd1
  np = np1

  meanerr = 99.

! the following quantities might be changed during iteration
  iitsobo = 4 !                     start of line-treatment with SOBO
  iitcmf = 21 !                     typical start of cmf
  iitupdate = 11 !                  hydro-update if .not. OPTCMF_FULL

! remains fixed
  If (optcmf_full) Then
    iitcmf_full = 30 !              start of cmf_full
    iitupdate = 41 !                hydro-update if OPTCMF_FULL
    first_cmf = .True.
  End If

  fastsobo = .False.

  precis = 1.D-6

! INITIALIZE ENIONND
  enionnd = 0.D0

! now it begins...

! part1------------------------------------------------------------------------

! read input data from indat.dat, copy to catalogue,
! and prepare specific quantities when reading from existing model

  Open (1, File='INDAT.DAT', Status='OLD')
  Rewind 1
  Read (1, Fmt=*) modnam


! JO August 2025 till Dec. 2025 changed.
! FASTWIND.LOG now contains important warnings and messages
! OPEN (999,FILE=TRIM(MODNAM)//'/FASTWIND.LOG',STATUS='UNKNOWN', &
! &                                            POSITION='APPEND')
  Open (999, File=trim(modnam)//'/FASTWIND.LOG', Status='UNKNOWN')
  Rewind 999

  Write (999, *) ' MODEL = ', trim(modnam)
  Write (999, *)
  Print *, ' MODEL = ', trim(modnam)
  Print *

  Call version_nlte_app(ver_nlte_app)
  Call version_cmf_all_sub(ver_cmf_all)

  Write (999, *) ' FASTWIND VERSION ', ver_nlte, ' (nlte)'
  Write (999, *) '          VERSION ', ver_nlte_app, ' (nlte_approx)'
  Write (999, *) '          VERSION ', ver_cmf_all, ' (cmf_all)'
  Write (999, *)
  Write (*, *) ' FASTWIND VERSION ', ver_nlte, ' (nlte)'
  Write (*, *) '          VERSION ', ver_nlte_app, ' (nlte_approx)'
  Write (*, *) '          VERSION ', ver_cmf_all, ' (cmf_all)'
  Write (*, *)

! CONTINUE READING INDAT.DAT
  Read (1, Fmt=*) optneupdate, helone, itlast, itmore
  Read (1, Fmt=*) optmixed
  Read (1, Fmt=*) teff, ggrav, srnom
  Read (1, Fmt=*) rmax, tmin
  Read (1, Fmt=*) xmloss, vmin, vmax, beta, vdiv
  Read (1, Fmt=*) yhein, hei
  Read (1, Fmt=*) optmod, optlucy, megas, accel, optcmf1
  Read (1, Fmt=*) vturb, xmet, lines, lines_in_model
  Read (1, Fmt=*) enatcor, expansion, set_first, set_step

  If (optcmf_full .And. itlast==iitcmf_full-1) Then
    Print *, ' Restart and ITLAST = ', itlast, &
      ' CORRESPONDING TO IITCMF_FULL-1'
    Print *, ' To ensure a reliable restart, use either ITLAST = ', &
      iitcmf_full - 2
    Print *, ' or ITLAST = ', iitcmf_full, &
      ' (in case at least one CMF_FULL iteration has already been performed)'
    Write (999, *) ' STOP: OPTCMF_FULL AND ITLAST = IITCMF_FULL-1'
    Stop ' OPTCMF_FULL AND ITLAST = IITCMF_FULL-1'
  End If

! no update of phot-struct if external model
  If (optmod) update_struct = .False.

!  If (optcmf_full) Then
!   start of Ng-extrapolation -- currently not used
!   If (update_struct) Then
!   ngstart = iitupdate + 6
!   Else
!   ngstart = iitcmf_full + 4
!   End If
!  End If

  Write (999, *) ' IMPORTANT CONTROL OPTIONS'
  Write (999, *) ' OPTCMF_FULL = ', optcmf_full, 'UPDATE_STRUCT = ', &
    update_struct, ' OPT_OIII_ITER = ', opt_oiii_iter
  Write (999, *)
  Write (*, *) ' IMPORTANT CONTROL OPTIONS'
  Write (*, *) ' OPTCMF_FULL = ', optcmf_full, 'UPDATE_STRUCT = ', &
    update_struct, ' OPT_OIII_ITER = ', opt_oiii_iter
  Write (*, *)


  If (optmod .And. opt_photclump) Then
    Write (999, *) ' STOP: external photospheric model .AND. photospheric &
      &clumping not implemented'
    Stop ' external photospheric model .AND. photospheric clumping &
      &not implemented'
  End If

! JO changed Dec. 2015
! no update if no temp. correction
  If (.Not. enatcor) update_struct = .False.

  If (optmod .And. itlast==0) Then
!   IF (ENATCOR) STOP ' Kurucz/Detail model to be used, but ENATCOR=.TRUE.'
    Inquire (File=trim(modnam)//'.michel', Exist=exi)
    If (.Not. optlucy .And. .Not. exi) Then
      Write (999, *) &
        ' STOP: external photospheric model to be used, but OPTLUCY=.FALSE.'
      Stop ' external photospheric model to be used, but OPTLUCY=.FALSE.'
    End If
  End If

! ------------------------
! input of clumping prescription (thin or thick)

  Read (1, Fmt='(A)') clf_para_str

  If (trim(adjustl(clf_para_str))=='THICK') Then
    optthick = .True.
    Read (1, Fmt='(A)') clf_para_str
!   read clumping factor
  Else
    optthick = .False.
!   clumping factor already read, no further action required
  End If

  If (optthick .And. optcmf_full) Then
    Write (999, *) ' STOP: OPTICALLY THICK CLUMPING NOT YET IMPLEMENTED &
      &IN FULL CMF TREATMENT'
    Stop ' OPTICALLY THICK CLUMPING NOT YET IMPLEMENTED IN FULL CMF TREATMENT'
  End If
! trick to read unknown number of variables (including description)
  Open (2, Status='SCRATCH', Form='FORMATTED')
  Write (2, *) clf_para_str
  Rewind (2)
  Read (2, Fmt=*, End=100, Err=100)(para_clf_input(i), i=1, npara_clf)
  Close (2)

100 Do i = 1, npara_clf
    If (para_clf_input(i)==-1.D99) Exit
  End Do

  npara_clf = i - 1 !               (one for the exit)

  Allocate (clf_field(npara_clf))

  Do i = 1, npara_clf
    clf_field(i) = para_clf_input(i)
  End Do

! check for TYPICAL clumping factor
  If (clf_field(1)<1.D0) Then
    Write (999, *) ' STOP: FIRST CLUMPING PARAMETER < 1!'
    Stop ' FIRST CLUMPING PARAMETER < 1!'
  End If
! Below parameters controlling if optically thick clumping should
! be considered, or if old assumption of only optically thin clumps
! should be used
  Allocate (fic_field(npara_clf))
  Allocate (fvel_field(npara_clf))
  Allocate (hpor_field(npara_clf))
  If (optthick) Then
    Read (1, Fmt=*) fic_field
    Read (1, Fmt=*) fvel_field
    Read (1, Fmt=*) hpor_field
  Else
    fic_field = 0.0
!   JO changed Aug 2016, since now normalized
    fvel_field = 1.0
    hpor_field = 0.0
  End If

! for tests
! DO I=1,NPARA_CLF
! PRINT*,FIC_FIELD(I),FVEL_FIELD(I),HPOR_FIELD(I)
! ENDDO
! STOP ' JS-TESTING!'

! end of clumping input
! ------------------------

  optmet = .False.
  If (xmet/=0.) optmet = .True.

! check for consistency
  If (.Not. optmet .And. lines) Then
    Write (999, *) ' STOP: XMET = 0., BUT LINE-BLOCKING REQUIRED!'
    Stop ' XMET = 0., BUT LINE-BLOCKING REQUIRED!'
  End If

  If (.Not. lines .And. lines_in_model) Then
    Write (999, *) ' STOP: NO LINE-BLOCKING, BUT REQUIRED IN MODEL!'
    Stop ' NO LINE-BLOCKING, BUT REQUIRED IN MODEL!'
  End If

  If (enatcor) Then
    If (lines) Then
      If (set_step<2) Then
        Write (999, *) ' STOP: METAL LINE BACKGROUND AND TCORR: &
          &MINIMUM VALUE FOR SET_STEP = 2!'
        Stop &
          ' METAL LINE BACKGROUND AND TCORR: MINIMUM VALUE FOR SET_STEP = 2!'
      End If
    End If
    If (lines .And. .Not. lines_in_model) Then
      Write (999, *) ' STOP: TCORR AND LINES: SET LINES_IN_MODEL TO .TRUE.!'
      Stop ' TCORR AND LINES: SET LINES_IN_MODEL TO .TRUE.!'
    End If
  End If

! ------------------------
! output of input data, both to model-dir (INDAT.DAT), to output-file, and to
! FASTWIND.LOG

  Write (999, *)
  Write (999, *) ' NUMBER OF PARAMETERS FOR CLUMPING FACTOR:', npara_clf
  Write (999, *) ' CORRESPONDING PARAMETERS:'
  Write (999, *) clf_field
  Write (999, *)
  Print *
  Print *, ' NUMBER OF PARAMETERS FOR CLUMPING FACTOR:', npara_clf
  Print *, ' CORRESPONDING PARAMETERS:'
  Print *, clf_field
  Print *
  If (optthick) Then
    Print *, ' OPTICALLY THICK CLUMPING'
    Print *
    Write (999, *) ' OPTICALLY THICK CLUMPING'
    Write (999, *)
  Else
    Print *, ' OPTICALLY THIN CLUMPING'
    Print *
    Write (999, *) ' OPTICALLY THIN CLUMPING'
    Write (999, *)
  End If

  Open (2, File=trim(modnam)//'/INDAT.DAT', Status='UNKNOWN')
  Rewind 2
  Write (2, Fmt=*) '''', modnam, ''''
  Write (2, Fmt=*) optneupdate, helone, itlast, itmore
  Write (2, Fmt=*) optmixed
  Write (2, Fmt=*) teff, ggrav, srnom
  Write (2, Fmt=*) rmax, tmin
  Write (2, Fmt='(E12.6,4(F10.4))') xmloss, vmin, vmax, beta, vdiv
  Write (2, Fmt=*) yhein, hei
  Write (2, Fmt=*) optmod, optlucy, megas, accel, optcmf1
  Write (2, Fmt=*) vturb, xmet, lines, lines_in_model
  Write (2, Fmt=*) enatcor, expansion, set_first, set_step
  If (optthick) Write (2, Fmt=*) 'THICK'
  Write (2, Fmt='(10(G10.4,2x))')(clf_field(i), i=1, npara_clf)

  If (optthick) Then
    Write (2, Fmt='(10(G10.4,2x))')(fic_field(i), i=1, npara_clf)
    Write (2, Fmt='(10(G10.4,2x))')(fvel_field(i), i=1, npara_clf)
    Write (2, Fmt='(10(G10.4,2x))')(hpor_field(i), i=1, npara_clf)
  End If

! end output of input data
! ------------------------

  If (optcmf_full) Open (18, File=trim(modnam)//'/OCC_NO', Status='UNKNOWN')
! occupation numbers of specified N and O levels


! ------------------------
! data related to temp. and temp. correction

! setup of TC
  ittcor = itlast + set_first
  opttcor = .False.

  If (enatcor) Then
!   TOTAL NUMBER OF T-CORRECTIONS
    ncor = -1
    emaxtc = 1.D6
    If (itlast==0) Then
      Open (10, File=trim(modnam)//'/MAXTCORR.dat', Status='REPLACE')
    Else
      Open (10, File=trim(modnam)//'/MAXTCORR.dat', Status='UNKNOWN')
      Rewind 10
      Do
!       IF RESTART FROM A MODEL WITH T-CORR, READ T-CONVERGENCE FROM PREVIOUS
!       RUN
        Read (10, *, End=110) ncor, emaxtc
      End Do
    End If

!   position file just before EOF
110 Backspace 10
    ncor = ncor + 1

!   CHECK WHETHER PREVIOUS RUN HAS BEEN CONVERGED
    If (.Not. optcmf_full) Then
      If (emaxtc<3.D-3) Then
        enatcor = .False.
        temp_converged = .True.
      End If
    Else
      If (emaxtc<3.D-3 .And. itlast>=iitcmf_full) Then
        enatcor = .False.
        temp_converged = .True.
      End If
    End If
  End If

  Write (999, *)
  Write (999, *) ' >>> VTURB =', vturb, ' km/s'
  Write (999, *)
  Print *
  Print *, ' >>> VTURB =', vturb, ' km/s'
  Print *

! ----vturb in cm/s

  vturb = vturb*1.0D5

! definition of hopfparameters

  If (.Not. optlucy) Then
    Read (1, Fmt=*) hopfself
    Write (2, Fmt=*) hopfself
    If (hopfself) Then
      Call hopfpara(teff, ggrav, yhein, xmet, xmloss, vmax, srnom, &
        clf_field(1))
    Else
      Read (1, Fmt=*) qinf, q0, gamhopf
      Write (2, Fmt=*) qinf, q0, gamhopf
    End If
  End If
! end temperature data
! ------------------------


! ------------------------
! data related to abundances

! prespecified abundances (explicit and background elements), EXCEPT FOR He
! if explicit element, abundance from DETAIL input file will be overwritten
! if background element, rescaled solar abundance will be overwritten

! input RELATIVE TO HYDROGEN and ALWAYS LOGARITHMICALLY
! either with respect to H=12 (positive number)
! or directly, log10(abund/H), i.e., negative

! if log10(abund/H) gt 0 (e.g., for Y > 1), use only first possibility
! (with respect to H=12), since positive numbers always interpreted in this
! way

! NOTE THAT THESE NUMBERS WILL BE USED DIRECTLY, i.e.,
! NOT RENORMALIZED IN ROUTINE ABUNDAN

  optxray = .False.
  nat_tot = 0
  Do i = 1, 31
    Read (1, Fmt=*, End=120) names_changed(i), abchanged(i)
!   TAKE CARE: if all 30 elements would be changed, the line containing
!   the x-ray parameters would not be read, thus I_max = 31
    If (names_changed(i)=='XRAYS') Then
      optxray = .True.
      fx = abchanged(i)
      Go To 120
    End If
    nat_tot = nat_tot + 1
  End Do

120 Continue

  If (nat_tot>0) Then
    Print *
    Print *, 'NUMBER OF ELEMENTS WITH SPECIFIED ABUNDANCES = ', nat_tot
    Print *
    Do i = 1, nat_tot
      Write (2, Fmt=*) names_changed(i), abchanged(i)
      If (abchanged(i)>0.) abchanged(i) = abchanged(i) - 12.D0
      abchanged(i) = 10.D0**abchanged(i)
      Print *, 'ELEMENT = ', names_changed(i), &
        'prespec. abundance (relative to H) = ', abchanged(i)
    End Do
    Print *
  End If

! end abundance data
! ------------------------


! ------------------------
! data related to X-rays

  Open (3, File=trim(modnam)//'/INXRAY.DAT', Status='UNKNOWN')
  If (optxray) Then
    If (.Not. optmet) Then
      Write (999, *) ' STOP: XMET = 0., BUT XRAY-TREATMENT'
      Stop ' XMET = 0., BUT XRAY-TREATMENT'
    End If

    If (teff<25000.) Then
      Write (999, *) ' STOP: TOO LOW TEFF FOR XRAY-TREATMENT'
      Stop ' TOO LOW TEFF FOR XRAY-TREATMENT'
    End If
!   UINFX IN KM/S
    Read (1, Fmt=*) gamx, mx, rminx, uinfx, px
    Rewind 3
    Write (3, Fmt=*) optxray, fx
    Write (3, Fmt='(5(G10.4,2x))') gamx, mx, rminx, uinfx, px
    Write (2, Fmt=*) 'XRAYS', fx
    Write (2, Fmt='(5(G10.4,2x))') gamx, mx, rminx, uinfx, px
  Else
    Write (3, Fmt=*) 'NO XRAY TREATMENT'
  End If

  Close (1)
  Close (2)
  Close (3)

! end X-ray data
! ------------------------

! end of input

! end
! part1--------------------------------------------------------------------

! prepare and check additional files (particularly if restart)


! if approximate approach (OPTCMF_FULL = .FALSE.), delete specific files
! which might have survived from previous models
  If (.Not. optcmf_full) Call delete_old_model_files

  If (optxray) Call read_kshell_data

  If (optcmf1 .And. optmixed/=1. .And. optmixed/=0.) Then
    Write (999, *) ' STOP: ERROR IN OPTMIXED'
    Stop ' ERROR IN OPTMIXED'
  End If


! ***  check for consistency: file indexcmf has to present under certain
! ***  conditions

  If (optmixed==1 .And. itlast>=20) Then
    Open (1, Err=270, File=trim(modnam)//'/INDEXCMF', Status='OLD')
    Close (1)
  End If

  Open (8, File=trim(modnam)//'/CONVERG', Status='UNKNOWN', Form='FORMATTED')

  If (optxray) Open (14, File=trim(modnam)//'/XLUM_ITERATION', &
    Status='UNKNOWN', Form='FORMATTED')
  Open (9, File=trim(modnam)//'/CONVERG_METALS', Status='UNKNOWN', &
    Form='FORMATTED')

! **** initializing of iqual (for non-treated levels)

  iqual = 0

! JO Sept 2023 subr. START modified
  Call start(nrec, teff, yhein, optmet)

! JO Sept. 2021
! check whether file ADD_FREQUENCIES exist, and delete if starting with
! iteration 0
  nstark = 0 !                      initialize
  If (optstark) Then

    If (optcmf_full) Then
      Write (999, *) ' STOP: OPTCMF_FULL .AND. OPTSTARK'
      Stop ' OPTCMF_FULL .AND. OPTSTARK'
    End If

    Inquire (File=trim(modnam)//'/ADD_FREQUENCIES', Exist=exi)
    If (exi) Then
      If (itlast==0) Then
        Open (1, File=trim(modnam)//'/ADD_FREQUENCIES', Status='OLD')
        Close (1, Status='DELETE')
        Print *
        Print *, 'old file ADD_FREQUENCIES deleted'
      Else
        Open (1, File=trim(modnam)//'/ADD_FREQUENCIES', Status='OLD')
        Read (1, *) nstark
        If (nstark>nstarkmax) Then
          Write (999, *) ' STOP: NSTARK > NSTARKMAX (after input)'
          Stop ' NSTARK > NSTARKMAX (after input)'
        End If
        Do i = 1, nstark
          Read (1, *) wavestark(i)
        End Do
        Close (1)
        Print *
        Print *, nstark, ' additional freq. points (Stark) read'
      End If
    Else
!     current hack, needs to be improved
      If (itlast/=0) Print *, ' RESTART, BUT FILE ADD_FREQUENCIES MISSING'
!     IF(ITLAST.NE.0) STOP ' RESTART, BUT FILE ADD_FREQUENCIES MISSING'
    End If

  End If

! part2------------------------------------------------------------------------

! .... establish photospheric structure with approximate grad
! (see Santolaya-Rey+ 1993)
! OR

! ... enable restart
! OR

! .... update photospheric structure (with actual grad)
! OR

! ...  no change in hydro structure (ICONVER=1) when T has been corrected

! PLUS -- in all above cases -- update of LTE occupation numbers (fg and bg)
! to account for changes in T and/or rho

  vmin1 = 0.D0
  xmmax = 200.D0

  iit = itlast + 1

! ***** new LTE/NLTE-cycle (first iteration or after T-update) *****
130 Continue


! establish different paths, in dependence of ITLAST and OPTTCOR
! basic control via ICONVER (structure converged) and
! ITHYDS (number of iterations to be performed)

! OPTMOD: HERE: model is present -- or (see modvl) -- external structure shall
! be used
! OPTMODEL: external structure shall be used

! if ENATCOR and NLTE iteration no. (IIT-ITLAST) > 1, this point should be
! only reached
! if OPTTCOR=T
  If (.Not. opttcor .And. enatcor .And. iit-itlast>1) Then
    Write (999, *) ' STOP: SOMETHING WRONG IN OPTTCOR/HYDRO-STRUCT PHILOSOPHY'
    Stop ' SOMETHING WRONG IN OPTTCOR/HYDRO-STRUCT PHILOSOPHY'
  End If

  If (itlast==0 .And. .Not. opttcor) Then
!   start from scratch, and temperature not corrected in iteration before
!   (IIT-1)
    iconver = 0
    ithyds = 1
    optmodel = .False.
    If (optmod) Then
      Print *, ' External photospheric structure will be used'
      optmodel = .True.
      optmod = .False.
    End If
  Else
!   either restart, or temperature corrected

!   .....start from previous model!  ... or T-corr!
!   .....check whether file 'temp' was previously written
!   .....to allow restart in old philosophy

    Rewind 21
    Read (21, *, Err=140) rtemp(nd)
    If (rtemp(nd)/=1.D0) Then
      Write (999, *) ' STOP: SOMETHING ROTTEN WITH RTEMP(ND) IN FILE TEMP'
      Stop ' SOMETHING ROTTEN WITH RTEMP(ND) IN FILE TEMP'
    End If
!   only one "iteration" to initialize LTE values
!   using old (restart) or new (corrected) temperature
    iconver = 1
    ithyds = 29
    optmod = .True.
    Go To 150

140 Continue
!   file temp does not exist, thus start from scratch
    Print *, ' HYDRO MODEL STARTED FROM SCRATCH, SINCE NO FILE "TEMP"'
    iconver = 0
    ithyds = 1
    optmodel = .False.
    If (optmod) Then
      Print *, ' External photospheric structure will be used'
      optmodel = .True.
      optmod = .False.
    End If
150 Continue
  End If

! photospheric
! struct.---------------------------------------------------------

ithydro: Do ithyd = ithyds, 30

!   for tests
!   print*,'iit =',iit,ithyds,iconver,ithyd,itlast,opttcor,optmod

!   either converge photosph. structure (ITHYDS=1), or
!   allow for restart,  update structure,
!   or account for an update of T-struct (all latter three ITHYDS=29)
!   NOTE: update of phot. struct. implies OPTTCOR=T (= ICONVER=1 and
!   ITHYDS=29),
!   independent of ITLAST

    If (iconver==1) Then
      If (itlast==0 .And. .Not. opttcor) Then
        Write (999, *) ' STOP: ERROR IN START MODEL'
        Stop ' ERROR IN START MODEL'
      End If
      Go To 170
    End If

    relkap = 1.D0
    relex = 1.D0
    If (ithyd>3) Then
      relkap = abs(1.D0-xkapith(ithyd-1)/xkapith(ithyd-2))
      relex = abs(1.D0-exith(ithyd-1)/exith(ithyd-2))
      relhyd = max(relkap, relex)
      If (relhyd<0.01D0) iconver = 1
!     this hack pretends convergence if H-recombination oscillates
      If (teff<12000. .And. ithyd>=20) iconver = 1

      If (update_struct .And. ithyd>=20) Then
!       this hack pretends convergence if oscillations of different origin
!       present
!       (but only if correlation high enough)
!       limit of 20 and corr > 0.97 from inspection of ob_grid
!       should not matter since structure updated later on anyway
        If (abs(corrith(ithyd-1))>=0.97) iconver = 1
      End If
    End If

!   ---- corvm added to correct (if neccessary) vmin to get proper
!   ---- tau-ross maximum (e. santolaya, 11/2/94)
    optcv = .False.
160 Continue

    Print *
    Print *, ' -------------- START VALUES FOR NE,TEMP,TAU-E-------------'
    Print *, ' ---------------ITHYD = ', ithyd
    Print *

!   input/output: ACTUAL XNE (lte), temp values (ne,nh corrected for clumping)
!   this routine is called only at iteration 0 AND the first time when
!   ICONVER=1

    Call model(taue, vmin1, .True., ithyd, xkap, exkap, temp, xmmax, rtau23, &
      xnh, xnelte, vdiv, ndiv)
!   ---- here's where parameters get defined (incl. clfac)
    Print *
    Print *, ' -------------------------------------------------------------'
    Print *

    If (iconver==1) Then
      Print *
      Print *, &
        ' +++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++'
      Print *, '            FINAL ITERATION NO . ', ithyd, ' FOR HYDRO MODEL'
      Print *, &
        ' +++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++'
      Print *
    End If

!   ---- file 'model' is read, particularly when ICONVER=1 (see above)

170 Continue
    Rewind 20

!   XNE(lte),nh corrected for clumping
!   last entry is now NDIV

    Print *
    Print *, ' MODEL FILE READ'
!   xnelte dummy except for first iteration(s)
    Read (20, Err=180) teff, ggrav, sr, yhe, xmu, vmax, xmloss, beta, &
      (r(i), i=1, nd), (velo(i), i=1, nd), (dvdr(i), i=1, nd), &
      (rho(i), i=1, nd), (xnelte(i), i=1, nd), (xnh(i), i=1, nd), &
      (clfac(i), i=1, nd), (index(i), i=1, nd), (pressure(i), i=1, nd), &
      ndiv_calc, updated, xt

    Go To 190
!   OLD VERSION OF MODEL FILE, SET UPDATE_STRUCTURE TO FALSE
180 update_struct = .False.
    Print *, ' WARNING WARNING WARNING!'
    Print *, ' OLD VERSION OF MODEL'
    Print *, ' UPDATE-STRUCT SET TO FALSE!'
    Print *

    Write (999, *)
    Write (999, *) ' ********************************************'
    Write (999, *) ' WARNING WARNING WARNING!'
    Write (999, *) ' OLD VERSION OF MODEL'
    Write (999, *) ' UPDATE-STRUCT SET TO FALSE!'
    Write (999, *) ' ********************************************'
    Write (999, *)

!   ---- comment upper lines and uncomment following line if UPDATE_STRUCT for
!   models
!   ---- with old version of MODEL-file should be performed
!   80   UPDATED=.FALSE.

190 Continue

!   ---- in case, read in clumping params
!   ---- standardwise, not required, since quantities stored in  module
!   nlte_porvor
!   Note that when skipping this read, the solution will be slightly
!   different,
!   since CLUMPING_OUTPUT is written with a lower precision than present in
!   the module
!   reading in the parameters is required after restart!

    If (optthick .And. itlast/=0 .And. iit==itlast+1) Then

      Open (13, File=trim(modnam)//'/CLUMPING_OUTPUT', Status='UNKNOWN')
      Read (13, *) !                header
      Do l = 1, nd1
        Read (13, Fmt=*) i, co_dum, co_dum, co_dum, co_dum, fic(l), fvel(l), &
          hpor(l), fvol(l), tcl_fac_line(l), tcl_fac_cont(l)
      End Do
      Close (13)
      Print *
      Print *, ' OPTTHICK AND RESTART: CLUMPING PARAMETERS READ'
    End If

    Print *
    vsound = sqrt(teff*akb/(xmu*amh))
    Print *, ' VSOUND = ', vsound/1.D5, ' KM/S'


!   ---- set up of p-z-grid

    Call pzgrid(nd, np, r)

!   ---- set up of electron temperature

!   ---- file 'temp' is read
!   ---- Temp updated in LTE (hydro-cycle) or in TEMPCORR_MAIN
!   ---- xnelte updated in LTE

    Rewind 21
    Read (21, *)(rtemp(i), i=nd, 1, -1), (temp(i), i=nd, 1, -1)
    Do l = 1, nd
      If (abs(rtemp(l)-r(index(l)))>1.D-10) Then
        Write (999, *) ' STOP: RTEMP NOT COMPATIBLE WITH R'
        Stop ' RTEMP NOT COMPATIBLE WITH R'
      End If
    End Do

!   .....input of XNE(lte) from converged hydro-model (in case of No
!   T-correction)
!   (actually, value from next to last iteration, to ensure restart with
!   identical numbers)
!   and rtau23.

!   in the case of T-correction, xnelte has been updated in TEMPCORR_MAIN
!   and refers to the last T-iteration
    If (iconver==1 .And. optmod) Then
      Read (21, *)(xnelte(i), i=nd, 1, -1) ! includes clumping
      Read (21, *) rtau23
    End If


!   -----in case, update of photo-stratification

    If (lte_update .And. update_struct) Then
      If (.Not. updated .And. iit>=iitupdate .And. emaxtc<0.03) Then
!       -----here, SRNOM is Rstar (in Rsun), and GGRAV is 10^logg
!       -----needs to be transfered into module nlte_var inside routine
        Call modvl_update(nd, r, velo, dvdr, rho, xne, xnelte, xnh, temp, &
          taue, taur, clfac, teff, ggrav, srnom, xmloss, beta, ndiv_calc, &
          yhein, yhe, hei, rtau23, update_done)

!       -----update p-z-grid
        If (update_done) Then
!         -----after actual update (for OPTCMF_FULL, since larger
!         corrections),
!         the next NLTE iteration should be done with LTE electron densities,
!         since the previous NLTE values are out of order.
!         This will be controlled by the following variable
          If (optcmf_full) after_update = .True.
          Call pzgrid(nd, np, r)
!         -----relax for 4 additional iterations
!         immediate update destroys lower boundary
          ittcor = ittcor + 4
!         -----reset all start variables
          start_weigh11 = .True.
          start_cmfsing = .True.
          start_cmfmulti = .True.
          start_jbar_metals = .True.
          start_cmf_simple = .True.
!         comment: recalculate lambda
        End If
      End If
    End If


!   -----set up of lte numbers and spherical grey temperature (if grey =
!   -----.true.)

    grey = .Not. optmod

!   to obtain consistent temperature stratification
!   this is the new statement from ver 7.3 (restart capab. not tested)
    If (iconver==1 .And. lines) lines_in_model = .True.

    restart_met_lte = .False.
!   if METAL-RESTART, to obtain consistent metal-opacities
    If (itlast/=0 .And. iit==itlast+1 .And. optmet) restart_met_lte = .True.
!   xnelte is updated and transfered to "common block" in module nlte_var and
!   written to TAU_ROS
!   JO April 2016: if OPTCMF_ALL or RESTART, LTE value of taur is replaced by
!   CMF/NLTE value (from previous content of TAU_ROS)
    If (optcmf_all .Or. (optcmf_full .And. itlast>=iitcmf_full)) Then
      update_taur = .True.
    Else
      update_taur = .False.
    End If
    Call lte(teff, sr, taur, taue, r, velo, dvdr, rho, xnelte, xnh, temp, &
      dilfac, clfac, ilow, imax, nd, grey, corrfc, ntemp, corvm, optcv, ithyd, &
      xkap, exkap, xmmax, rtau23, xmet, restart_met_lte, lte_update, &
      update_taur)

!   ---- from here on, temp = actual non-lte temperature

    If (optcv .And. ithyd==1) Then
      vmin1 = velo(nd)*corvm*vmax*1.D-5
      Print *, ' '
      Print *, 'CORRECTION TO VMIN; NEW VMIN (KM/S)=', vmin1
      Print *, ' '
      Go To 160
    End If

!   check whether all operations concerning XNELTE where successful

    Do i = 1, nd
      If (xnelte(i)/=xnelte_common(i)) Then
        Write (999, *) ' STOP: XNELTE .NE. XNELTE (common)'
        Stop ' XNELTE .NE. XNELTE (common)'
      End If
    End Do
    Open (1, File=trim(modnam)//'/TAU_ROS', Status='UNKNOWN', &
      Form='FORMATTED')
    Rewind 1
    Do i = 1, nd
      Read (1, *) taur_taur(i), xnelte_taur(i)
      If (taur(i)/=0.D0 .And. abs(1.-taur(i)/taur_taur(i))>1.D-13) Then
        Write (999, *) ' STOP: TAUR .NE. TAUR IN TAU_ROS'
        Stop ' TAUR .NE. TAUR IN TAU_ROS'
      End If
      If (abs(1.-xnelte(i)/xnelte_taur(i))>1.D-13) Then
        Write (999, *) ' STOP: XNELTE .NE. XNELTE IN TAU_ROS' ! not exact,
!       since
!       formattedly
!       written/read
        Stop ' XNELTE .NE. XNELTE IN TAU_ROS' ! not exact, since formattedly
!       written/read
      End If
    End Do
    Close (1)

    If (iconver==1) Go To 200
  End Do ithydro

! **** from here on, hydro is converged or pretended to be converged

  Print *, ' ITERATION HISTORY OF KRAMER OPACITY REGR. COEFF.'
  Print *
  Do l = 1, ithyd - 1
    Write (*, 280) l, xkapith(l), exith(l), corrith(l)
    Print *
  End Do

  Write (999, *) ' ITERATION HISTORY OF KRAMER OPACITY REGR. COEFF.'
  Write (999, *)
  Do l = 1, ithyd - 1
    Write (999, *)
    Write (999, 280) l, xkapith(l), exith(l), corrith(l)
  End Do

  Write (999, *) ' STOP: WARNING!!!! HYDRO NOT CONVERGED IN 30 ITERATIONS'
  Stop ' WARNING!!!! HYDRO NOT CONVERGED IN 30 ITERATIONS'

! -----radiative transfer to obtain start values for j

200 Continue


! end
! part2---------------------------------------------------------------------

! all hydro-related stuff now done:
! approx. strat. when calculated from scratch,
! stratification read in when restart,
! or updated when required
! Moreover, LTE-values with current T- and rho-stratification calculated:

! final check for photospheric mass < fraction of total mass

  Call checkmass(sr, ggrav, rtau23, r, rho, nd, ndiv, ndiv_calc)

  frenew = .True.

! ---- u-function tabular values are read and stored (only once)

  Call ulecture

! ---- u2-function tabular values are read and stored (only once)

  Call u2lecture

  errold = precis

  specmat = 1.D0

! calculation of fine frequency grid

  Call fresfin(nd, xnh, temp(nd), .True., teff)

! redefinition of dilfac to ensure thermalization in nlte_approx

  If (optmet) Then
    Where (dilfac<0.5)
      dilfac1 = dilfac
    Elsewhere
      dilfac1 = 1.
    End Where
    Do l = 1, nd
      If (dilfac1(l)==1.) Then
        dilfac1(l) = 0.75 !         smooth transition
        Exit
      End If
    End Do
  End If

  lwion = l + 1
  taurlim = taurlim_default
! check whether TAURLIM is outside metal-lte range or too small
  If (enatcor .And. optmet) Then
    If (taurlim<taur(ndiv_calc+1)) taurlim = taur(ndiv_calc+1) ! STRONG WIND
    If (taurlim>2.) taurlim = 2. !  MAXIMUM VALUE EVEN FOR STRONG WINDS
    If (taurlim>taur(lwion)) taurlim = taur(lwion)
    Print *, ' TAURLIM        = ', taurlim, '   (DIVISION OF T-CORR SCHEMES)'
    Write (999, *)
    Write (999, *) ' TAURLIM        = ', taurlim, &
      '   (DIVISION OF T-CORR SCHEMES)'
  End If
  Print *, ' TAUROSS(LWION) = ', taur(lwion), &
    ' (DIVISION OF LTE/NLTE for non-seclected METALS)'
  Write (999, *) ' TAUROSS(LWION) = ', taur(lwion), &
    ' (DIVISION OF LTE/NLTE for non-seclected METALS)'

! specify transition point for selected ions
  Do l = lwion, nd
!   IF(TAUR(L).GT.2.D0) EXIT
!   JO March 2022 (11.5.1): increased in new version (taur=2 sometimes too
!   low)
    If (taur(l)>10.D0) Exit
  End Do

  lwion1 = l
  Print *, ' TAUROSS(LWION1) = ', taur(lwion1), &
    ' (DIVISION OF LTE/NLTE for seclected METALS)'
  Print *, ' FROM NLTE: LWION = ', lwion
  Print *, ' FROM NLTE: LWION1 = ', lwion1
  Write (999, *) ' TAUROSS(LWION1) = ', taur(lwion1), &
    ' (DIVISION OF LTE/NLTE for seclected METALS)'
  Write (999, *) ' FROM NLTE: LWION = ', lwion
  Write (999, *) ' FROM NLTE: LWION1 = ', lwion1
  Write (999, *)

  If (itlast/=0 .Or. opttcor) Go To 210

  Print *, '               ++++++++++++++++++++++++++'
  Print *, '               +     ITERATION NO 0     +'
  Print *, '               ++++++++++++++++++++++++++'
  Print *
  Write (999, *)
  Write (999, *) '               ++++++++++++++++++++++++++'
  Write (999, *) '               +     ITERATION NO 0     +'
  Write (999, *) '               ++++++++++++++++++++++++++'
  Write (999, *)

! XNE is input as LTE
  If (opttcor) Then
    Write (999, *) ' STOP: SOMETHING WRONG WITH OPTTCOR'
    Stop ' SOMETHING WRONG WITH OPTTCOR'
  End If
  lastlte = .True.
  rybicki = .True.

! JO Oct. 2021: here was a bug (IIT missing in call)
  Call cont(nd, np, r, rho, xnelte, temp, clfac, teff, sr, corrfc, optlte, &
    concon, taur, rtau23, rybicki, dtfcorr, fluxerr, ndiv_calc+1, iit)

  If (optmet) Then
    Print *, ' BEGIN OF APPROXIMATE NLTE FOR METALS (It 0, xne(LTE))'
    Write (999, *) ' BEGIN OF APPROXIMATE NLTE FOR METALS (It 0, xne(LTE))'
    Call trad(nd, temp, dilfac1, xnelte, clfac, iit, optcmf1, lastlte)

!   CALL NLTE_APPROX(TEFF,YHE,XMET,RHO,XNELTE,DILFAC1,TEMP, &
!   &       DVDR*VMAX/SR,VELO/R*VMAX/SR,VELO*VMAX,CLFAC,ND, &
!   &       FRENEW,.FALSE.,1,.TRUE.)
!   dummy arguments changed (Jan 2015)
    Call nlte_approx(teff, yhe, xmet, rho, xnelte, dilfac1, temp, r, velo, &
      dvdr, clfac, ilownl, imaxnl, nd, frenew, .False., 1, .True.)

    Call opacitm(nd, xnelte, temp, clfac, .False., 1, .True., frenew, .True.)

    Print *, ' END OF APPROXIMATE NLTE FOR METALS'
    Print *
    Write (999, *) ' END OF APPROXIMATE NLTE FOR METALS'
    Write (999, *)
!   JO Feb. 2017
    Call update_occng(.False.) !    since LTEOPT=.FALSE.
  End If

! INPUT: LTE, OUTPUT NLTE (updated in UPDATE!)
  xne = xnelte
! print*,'xnecheck: iteration zero'
! print*,xnelte
! print*,xne

! here, we do not need ilow and imax for the metals, since they
! are calculated only from the line-iteration on (CONCON = T).
! ilow and imax (H/He) from ILOWIMAX_LTE

  Call rateeq(xnh, xne, temp, clfac, ilow, imax, nd, err, concon, r, velo, &
    errold, specmat, accel, optmet, meanerr)

  frenew = .False.

  errold = err

210 Continue

  If (itlast/=0 .Or. opttcor) Then

!   if start from existing model, then update of electron density
!   and (ni/nk)*, which are normally updated in UPDATE, but at this
!   stage defined by the lte-electron density.
!   input XNE(lte), output XNE(nlte)

    Print *, ' RESTART FROM PREVIOUS ITERATION'
    Write (999, *) ' RESTART FROM PREVIOUS ITERATION'
    xne = xnelte
!   ilow and imax only dummy arguments here, to preserve RESTART
!   locate new levels for .not.(concon.and.optmet)
    Call restart(nd, xne, dilfac1, optmet, ilow, imax)
!   here was the bug; update nlte electron density by actual lte e-density
!   could be done also in subroute copyocc, but here it's more clear
    If (.Not. optneupdate) xne = xnelte
!   use lte electron densities just after update (for OPTCMF_FULL)
!   variable will be reset to .false. in OPACITM
    If (after_update) xne = xnelte
!   print*,xnelte
!   print*,xne

!   call to select now in subr. RESTART
  End If

  optlte = .False.
  ntemp1 = 1
  emax = 1000.D0

! part
! 3-----------------------------------------------------------------------

! Begin of the NLTE-cycle iteration

220 If (iit==itlast+itmore+1) Go To 230
! from here on, only NLTE values for XNE as input

  Print *, '               ++++++++++++++++++++++++++'
  Print *, '               +  ITERATION NO ', iit, '  +'
  Print *, '               ++++++++++++++++++++++++++'
  Print *
  Write (999, *)
  Write (999, *) '               ++++++++++++++++++++++++++'
  Write (999, Fmt='(A,I6,A)') '                +  ITERATION NO ', iit, '   +'
  Write (999, *) '               ++++++++++++++++++++++++++'
  Write (999, *)

  optcmf = optcmf1
  If (iit<=20 .And. .Not. fastsobo) Then
    optcmf = .False.

!   if fast convergence of sobo

    If (concon .And. optcmf1 .And. emax<=.1 .And. (.Not. update_struct .Or. &
      updated)) Then
      fastsobo = .True.
      optcmf = .True.
      iitcmf = iit
    End If
  End If

  If (iit>3) concon = .True.

! if fast convergence of continuum

  If (iit<=3 .And. emax<=.1) Then
    concon = .True.
    iitsobo = iit
  End If

! ---- complete output for the last iteration

  If (iit==itlast+itmore) unasol = .True.

! ---- file for complete output

  If (unasol .And. megas) Open (100, File=trim(modnam)//'/megasalida', &
    Status='UNKNOWN', Form='FORMATTED')

! print*,'xnecheck: normal iterations'
! print*,xne

! part3a-----------------------------------------------------------------------

! specification of ILOWNL and IMAXNL

  If (.Not. concon) Then
!   only H/He values for ilownl/imaxnl are needed,
!   as well as LTE values for ilow/imax.
!   These have been calculated already in ILOWIMAX_LTE
    ilownl = ilow
    imaxnl = imax
  End If
! if CONCON, metals are switched on; thus we need ILOWNL and IMAXNL
! for all routines, including OPACITC.
! Usually, we obtain them from our approximate NLTE solution
! (met_imin, met_imax, previous iteration).
! Only in case z=0, we have already calculated them in ILOWIMAX_OLD.
! ILOWNL, IMAXNL must be updated only when the temperature has changed,
! and we are not close to convergence. THUS

! IF (CONCON .AND. (FIRST_METALS.OR. (LTE_UPDATE .AND. EMAXTC.GE.0.05))) THEN
! JO changed Sept. 2014, to allow for better guess
! IF (CONCON .AND. (FIRST_METALS.OR. (LTE_UPDATE .AND. EMAXTC.GE.0.03))) THEN
! JO CHANGED (v10.4): one final update of ILOWNL etc. after temperature has
! converged, &
! to allow a consistent restart
  If (concon .And. (first_metals .Or. (lte_update .And. &
    emaxtc>=0.03) .Or. (lte_update .And. temp_converged))) Then
    If (first_metals) first_metals = .False.
    If (optmet) Then
!     --- this is the subroutine to fiddle around if you are not satisfied
!     with
!     --- ILOW and/or IMAX !!! (in nlte_approx.f90)
      Call ilowimax_nlte(ilownl, imaxnl, ilow, imax, nd, ndiv_calc)
!     JO July 2018: after potential change in ionization, indices for
!     cmf transfer (expl. elements) have to be prepared in RATELI,
!     using restart-option
      If (optcmf_all) Then
        first_rateli = .True.
        restart_cmfall = .True.
      End If
    Else
      ilownl = ilow
      imaxnl = imax
      imianl = imia
      imaanl = imaa
    End If
  Else
!   Usually, ILOWIMAX_NLTE ensures that IMAX(NLTE) always LE IONIS(LTE).
!   Under very conspicous circumstances, however, IONIS might decrease
!   even if the the temperature correction becomes low
!   (remember the first CNOP test).
!   In this case then, IMAXNL would be too large since no LTE info would
!   be present any longer (leads to NISTAR = infinity).
!   Thus have to perform the following hack
    Do k = 1, nat
      Do l = 1, nd
        If (imaxnl(l,k)>ionis(k,l)) Then
          imaxnl(l, k) = ionis(k, l)
          If (imaxnl(l,k)<ilownl(l,k)) Then
            Write (999, *) ' STOP: IMAXNL < ILOWNL'
            Stop ' IMAXNL < ILOWNL'
          End If
          Print *, ' CORRECTION OF IMAX(NLTE) REQUIRED!!!'
          Print *, k, l, ' IMAX(NLTE) = ', imaxnl(l, k), ' (ZEFF CORRECTED)'
          Write (999, *) ' CORRECTION OF IMAX(NLTE) REQUIRED!!!'
          Write (999, *) k, l, ' IMAX(NLTE) = ', imaxnl(l, k), &
            ' (ZEFF CORRECTED)'
        End If
      End Do
    End Do
  End If

! part3b-----------------------------------------------------------------------

! calculate run of shock-temperatures
  If (optxray) Then
    If (start_xray) Then
!     needs only to be iterated when the change in ne should be accounted for
!     in this case, LXMIN might change as well
!     JO: comment changed July 2016
!     BUT: ne is not changing, since it refers to the hot plasma;
!     and the total wind-structure is NOT changing even after hydro-update,
!     since only photosphere affected by this update.
!     THUS: Indeed, only one call of TSHOCK and LAMBDA_XRAY required
!     inconsistently, we use XMU from the cool wind. Does not matter!
      Call tshock(tsho, velo, r, rtau23, xmu, nd)
!     calculate X-ray cooling functions accounting for cooling zones, both for
!     radiative and adiabatic shocks, following Feldmeier et al. 1997,
!     OR, in case, Owocki+ 2013
!     same situation as for TSHOCK
      Call lambda_xray(xmu, tsho, r, velo, rho, xnelte, xnh, clfac, nd)
      Call interpol_xray !          usually, in the first call frenew = false
      start_xray = .False.
    End If

!   this needs to be reevaluated when frequency grid changes
    If (frenew) Call interpol_xray
  End If

! part3c-----------------------------------------------------------------------

! continuum opacities incl. K-shell ones, and continuum transport
  Call opacitc(nd, xne, xnh, temp, clfac, ilownl, imaxnl, optlte, ntemp1, &
    .True., frenew, optmet)
! STOP 'TESTING-JS!'
! THICK: HERE OPAC GETS UPDATED; it is the critical part of the effective
! opacity inclusion!

! rybicki solution only if larger changes in continuum to be ex-
! pected, i.e. after first nlte-solution of sobolev and cmf case

  rybicki = .False.
  If (iit==itlast+1) rybicki = .True.
  If (iit==iitsobo+1) rybicki = .True.
  If (iit==iitcmf+1) rybicki = .True.
  If (opttcor) rybicki = .True.

  Print *
  If (rybicki) Then
    Print *, ' CONTINUUM WITH RYBICKI !!!! '
    Write (999, *) ' CONTINUUM WITH RYBICKI !!!! '
  Else
    Print *, ' CONTINUUM WITHOUT RYBICKI '
    Write (999, *) ' CONTINUUM WITHOUT RYBICKI '
  End If
  Print *

  Call cont(nd, np, r, rho, xne, temp, clfac, teff, sr, corrfc, optlte, &
    concon, taur, rtau23, rybicki, dtfcorr, fluxerr, ndiv_calc+1, iit)

  lastlte = .False.
  If (opttcor) lastlte = .True. !   previous LTE-update

  opttcor = .False.
! no T-correction in final iteration (to obtain consistent temp. and
! LTE-occ.num.)
  If (enatcor .And. iit==ittcor .And. .Not. unasol) Then
    opttcor = .True. !              enables the update of the LTE numbers and
!   calculation of partial derivatives
    Call tempcorr_start !           starting up all the ThB quantities
  End If

! part3d-----------------------------------------------------------------------

! ****  begin CMF_ALL treatment

! start/restart

  If (optcmf_full) Then
    If (optcmf .And. optmet .And. global) Then
!     restart: read OCCNG, MET_IMIN, MET_IMAX from file OCCNG (only once)
!     JO: at some point, the input value of Jnu (for escat) has to be updated
!     in restart (but: thus far, no problems)
      If (.Not. optcmf_all .And. itlast>=iitcmf_full) Call optcmf_restart
      If (iit>=iitcmf_full) optcmf_all = .True.
      If (itlast<iitcmf_full) Then
        If (.Not. update_struct .And. iit<=iitcmf_full+5) &
          metals_converged = .False.
        If (update_struct .And. iit<=iitupdate+1) metals_converged = .False.
      End If
    End If
  End If

! ensure that rateli is called before cmf_all after restart (to prepare index
! files etc.)
  If (optcmf_all .And. first_rateli) Then
    If (xj(1,ifre+1)/=1) Then
      Write (999, *) ' STOP: Wrong XJ-flag before first call of RATELI'
      Stop ' Wrong XJ-flag before first call of RATELI'
    End If
    Call rateli(xne, temp, clfac, ilownl, imaxnl, nd, concon, r, velo, dvdr, &
      sr, vmax, teff, corrfc)
    Print *
    Print *, &
      ' RBB-RATES FOR EXPLICIT ELEMENTS PREPARED AFTER RESTART (RATELI)'
    Print *
    Write (999, *)
    Write (999, *) &
      ' RBB-RATES FOR EXPLICIT ELEMENTS PREPARED AFTER RESTART (RATELI)'
    Write (999, *)
  End If

  Print *, ' OPTIONS - CURRENT STATE: IIT,OPTCMF_ALL,METALS_CONVERGED,METCONV2,&
    &ALMOST_CONVERGED,NO_ALO_METALS,NO_ALO'
  Print *, iit, optcmf_all, metals_converged, metconv2, almost_converged, &
    no_alo_metals, no_alo
  Write (999, *) ' OPTIONS - CURRENT STATE: IIT,OPTCMF_ALL,METALS_CONVERGED,MET&
    &CONV2,ALMOST_CONVERGED,NO_ALO_METALS,NO_ALO'
  Write (999, *) iit, optcmf_all, metals_converged, metconv2, &
    almost_converged, no_alo_metals, no_alo

  If (optcmf_all) Then
!   perform complete CMF transport in a predefined range (wavblue < lambda
!   < wavred), for all significant elements.
!   line opacities from occngold(1...LWION1-1)/occnglte(LWION1...ND) for
!   selected elements and occng for approx. elements,

!   remap radiation field onto coarse grid (conserve freq. integrals),
!   and overwrite radiation field (from CONT, coarse grid) by remapped
!   cmf-quantities

!   at first: control of cmf-solution:
!   OPT_RAY = .TRUE. -> ray by ray + moments eqs.
!   OPT_RAY = .FALSE. -> only moments eqs., with
!   Eddington factors from previous ray-by-ray
    If (temp_converged) Then
      opt_ray = .True. !            to allow for smooth convergence
    Else
!     some preparations
      If (set_step/=2 .Or. set_first/=1) Then
        Write (999, *) &
          ' STOP: CHANGE SET_STEP/SET_FIRST, OR OPT_RAY-PHILOSOPHY!'
        Stop ' CHANGE SET_STEP/SET_FIRST, OR OPT_RAY-PHILOSOPHY!'
      End If
      If (first_cmf) Then
        iit_start_cmf = iit
        opt_ray_counter = 1
      End If
      If (update_done .And. iit_updated==0) Then
        iit_updated = iit
        opt_ray_counter = 1
      End If
!     now, set OPT_RAY
!     be aware that there might be some interference between alo and ng
!     (if ng not avoided before temp_converged).
      opt_ray = .False. !           except for following cases
      If (.Not. update_done) Then ! before hydro-update
        delta_iit = iit - iit_start_cmf
        If (delta_iit<0) Then
          Write (999, *) ' STOP: ERROR IN OPT_RAY-PHILOSOPHY(1)'
          Stop ' ERROR IN OPT_RAY-PHILOSOPHY(1)'
        End If
        If (opt_ray_counter<=4) Then
          optrayarr1 = opt_ray_arr1(opt_ray_counter)
        Else
          optrayarr1 = opt_ray_arr1(4) + (opt_ray_counter-4)*4
        End If
        If (delta_iit==optrayarr1) Then
          opt_ray = .True.
          opt_ray_counter = opt_ray_counter + 1
!         JO Aug 2025: obsolete
!         IF(OPT_RAY_COUNTER.GT.40) STOP ' OPT_RAY_COUNTER TOO LARGE (1)'
        End If
!       just before hydro_update
!       IF (UPDATE_STRUCT.AND..NOT.UPDATED.AND.IIT.GE.IITUPDATE &
!       &                .AND.EMAXTC.LT.0.03.AND.BEFORE_HYDRO_UPDATE) THEN
!       PRINT*,'BEFORE HYDRO-UPDATE'
!       OPT_RAY=.TRUE.
!       BEFORE_HYDRO_UPDATE=.FALSE.
!       ENDIF
      Else !                        after actual hydro-update
        delta_iit = iit - iit_updated
        If (delta_iit<0) Then
          Write (999, *) ' STOP: ERROR IN OPT_RAY-PHILOSOPHY(2)'
          Stop ' ERROR IN OPT_RAY-PHILOSOPHY(2)'
        End If
        If (opt_ray_counter<=4) Then
          optrayarr2 = opt_ray_arr2(opt_ray_counter)
        Else
          optrayarr2 = opt_ray_arr2(4) + (opt_ray_counter-4)*4
        End If
        If (delta_iit==optrayarr2) Then
          opt_ray = .True.
          opt_ray_counter = opt_ray_counter + 1
!         JO Aug 2025: obsolete
!         IF(OPT_RAY_COUNTER.GT.40) STOP ' OPT_RAY_COUNTER TOO LARGE (2)'
        End If
      End If
    End If

!   As it turned out, during temp-iteration it is of utmost importance to keep
!   ALO; with and without ALO the heating/cooling rates may change sign, and
!   a different (and wrong) T-structure (from electron thermal balance) may
!   arise in intermediate regions.

    Print *
    Print *, ' OPT_RAY CONTROL: ', iit, ' ', opt_ray, 'UPDATE_DONE = ', &
      update_done
    Print *
    Write (999, *)
    Write (999, *) ' OPT_RAY CONTROL: ', iit, ' ', opt_ray, 'UPDATE_DONE = ', &
      update_done
    Write (999, *)

!   Jo Dec 2020 additional check
    If (optlte) Then
      Write (999, *) ' LTEOPT = T in NLTE-path'
      Stop ' LTEOPT = T in NLTE-path'
    End If

    Call cmf_all(nd, xne, temp, clfac, r, velo, dvdr, rho, taur, xnh, &
      pressure, ferr_cmf, dtfcorr_cmf, ilownl, imaxnl, rmax, unasol, opt_ray)
    restart_cmfall = .False.
    first_cmf = .False.


!   for tests
!   DO L=1,ND
!   WRITE(*,FMT='(I3,3(2X,F10.4))'),L,LOG10(TAUR(L)),FERR_CMF(L),DTFCORR_CMF(L)
!   ENDDO
  End If

! **** end CMF_ALL treatment

! end
! part3d--------------------------------------------------------------------

! print fluxes and flux errors from subr. CONT and/or CMF_ALL to file FLUXCONT
  If (unasol) Then
    Call print_fluxes(r)
    Call print_fluxerror(r, taur, fluxerr, ferr_cmf)
  End If

! part3e----------------------------------------------------------------------

  If (optmet .And. .Not. metals_converged) Then
!   IF (OPTMET.AND.(.NOT.METALS_CONVERGED.OR.OPTCMF_FULL)) THEN
    Print *, ' BEGIN OF APPROXIMATE NLTE FOR METALS'
    Write (999, *) ' BEGIN OF APPROXIMATE NLTE FOR METALS'
    Call trad(nd, temp, dilfac1, xne, clfac, iit, optcmf1, lastlte)

!   xne_ratio required for opacity_cmf:
!   background occup. numbers are calculated with xne as input
!   but, before next call to cmf_all, xne has been updated,
!   so this needs to be saved here
    Do l = 1, nd
!     one final check
      If (xnelte(l)/=xnelte_common(l)) Then
        Write (999, *) &
          ' STOP: XNELTE NE XNE_COMMON (before call to nlte_approx)'
        Stop ' XNELTE NE XNE_COMMON (before call to nlte_approx)'
      End If
      xne_ratio(l) = xne(l)/xnelte(l)
!     PRINT*,XNE(L),XNELTE(L),XNE_RATIO(L)
    End Do

!   Jo Dec 2020 additional check
    If (optlte) Then
      Write (999, *) ' STOP: LTEOPT = T in NLTE-path'
      Stop ' LTEOPT = T in NLTE-path'
    End If

!   dummy arguments changed (Jan. 2015)
    Call nlte_approx(teff, yhe, xmet, rho, xne, dilfac1, temp, r, velo, dvdr, &
      clfac, ilownl, imaxnl, nd, frenew, optlte, ntemp1, .True.)

    Call opacitm(nd, xne, temp, clfac, optlte, ntemp1, .True., frenew, &
      .False.)
    Write (999, *) ' END OF APPROXIMATE NLTE FOR METALS'
    Write (999, *)
    Print *, ' END OF APPROXIMATE NLTE FOR METALS'
    Print *
  End If

! end
! part3e----------------------------------------------------------------------

! ---  calculation of rbb-rates (explicit elements) from old occup.numbers
! ---  for all depth points (cmf or sobolev)

  If (concon) Then
    first_rateli = .False. !        to control restart
    If (.Not. optcmf_all) Then
!     prepare/perform SA/CMF transport for individual lines,
!     and calculate TLU/TUL for explicit elements
      Call rateli(xne, temp, clfac, ilownl, imaxnl, nd, concon, r, velo, dvdr, &
        sr, vmax, teff, corrfc)
      Print *
      Print *, &
        ' RBB-RATES (SOBO/CMF_SING) FOR EXPLICIT ELEMENTS PREPARED (RATELI)'
      Print *
      Write (999, *)
      Write (999, *) &
        ' RBB-RATES (SOBO/CMF_SING) FOR EXPLICIT ELEMENTS PREPARED (RATELI)'
      Write (999, *)
    Else
!     calculate TLU/TUL for explicit elements, from complete cmf-solution
!     (and for cmfsing/SA transitions)
      Call calc_tlu(nd, ilownl, imaxnl, xne, temp, clfac, r, velo, dvdr)
    End If

  End If

! part3f----------------------------------------------------------------------

! --- now, full rateeq including actual values ilownl, imaxnl

  Call rateeq(xnh, xne, temp, clfac, ilownl, imaxnl, nd, err, concon, r, velo, &
    errold, specmat, accel, optmet, meanerr)

! --- prepare and perform Ng-extrapolation

! IF(OPT_NG .AND. OPTCMF_ALL .AND. IIT.GE.NGSTART) THEN
  flag_ng = .False.
  flag_ng_met = .False.
  If (opt_ng .And. optcmf_all .And. temp_converged) Then
!   --- extrapolation (affecting OPACITY_CMF, CALC_TLU and CALC_TLU_MET)
!   only performed when in corresponding iteration ALMOST_CONVERGED = F
!   No general statement for explicit elements, since ALMOST_CONVERGED
!   can oscillate
    Call prep_ng(ilownl, imaxnl, r, flag_ng)
    If (.Not. metals_converged) Call prep_ng_met(r, flag_ng_met)
  End If

  If (optmet .And. .Not. metals_converged .And. .Not. optlte) Then

!   --- reset occngold2 = occng (bg-elements)

    Call update_occng(optlte)
  End If

! end
! part3f-------------------------------------------------------------------

! part3g-----------------------------------------------------------------------

  If (enatcor .And. iit==ittcor .And. .Not. unasol) Then
    If (optcmf_all) Then
      dtfcorr = dtfcorr_cmf
      fluxerr = ferr_cmf
    End If

    Call tempcorr_main(temp, xne, taur, dtfcorr, teff, ncor, emaxtc, &
      dtmean, expansion, outthb)
    Print *
    Print *, ' MEAN RELATIVE CHANGE DT/T (LOG) = ', dtmean
    Print *
    Write (999, *)
    Write (999, *) ' MEAN RELATIVE CHANGE DT/T (LOG) = ', dtmean
    Write (999, *)
!   JO CHANGED Feb 2017: 21,31,41 changed to IITCMF, IITCMF_FULL+1, IITUPDATE
    iit_temp = iitcmf
    If (optcmf_full .And. .Not. update_struct) iit_temp = iitcmf_full + 1
    If (optcmf_full .And. update_struct) iit_temp = iitupdate
!   JO Feb 2016: convergence only after update
    If (.Not. update_struct .Or. update_struct .And. updated) Then
      If (emaxtc<3.D-3 .And. iit>iit_temp) Then ! prevent convergence in Sobo
!       cycle
        Print *, ' !!!!!!!!!!!!!!!!!!!!!'
        Print *, ' TEMPERATURE CONVERGED'
        Print *, ' !!!!!!!!!!!!!!!!!!!!!'
        Write (999, *) ' !!!!!!!!!!!!!!!!!!!!!'
        Write (999, *) ' TEMPERATURE CONVERGED'
        Write (999, *) ' !!!!!!!!!!!!!!!!!!!!!'
        enatcor = .False. !         no more temperature-corrections required
        temp_converged = .True.
!       JO June 2018: in case of bad t-convergence, allow for weaker condition
      Else If (dtmean<-3.3 .And. iit>max(iit_temp,80)) Then ! only if many
!       iterations have
!       been done
        Print *, ' !!!!!!!!!!!!!!!!!!!!!'
        Print *, ' TEMPERATURE adopted as CONVERGED, DTMEAN = ', dtmean
        Print *, ' !!!!!!!!!!!!!!!!!!!!!'
        Write (999, *) ' !!!!!!!!!!!!!!!!!!!!!'
        Write (999, *) ' TEMPERATURE adopted as CONVERGED, DTMEAN = ', dtmean
        Write (999, *) ' !!!!!!!!!!!!!!!!!!!!!'
        enatcor = .False. !         no more temperature-corrections required
        temp_converged = .True.
        emaxtc = 2.9D-3 !           reset, to be consistent with present
!       philosophy
      End If
    End If

    If (temp_converged) iit_temp_conv = iit

    Write (10, Fmt='(1X,I3,1X,G12.5,1X,I3)') ncor, emaxtc, iit

    ittcor = ittcor + set_step !    iteration for the next T-Corr
    ncor = ncor + 1 !               number of performed T-Corr
!   At this point, T has change and the TEMP file has been also written
    Print *, ' '
!   JO Jan 2016: reset after first temp. correction following photospheric
!   update
    If (no_check) no_check = .False.
  End If

! end
! part3g-------------------------------------------------------------------

! part3h-----------------------------------------------------------------------

! perform some final calculations (for extrapolations), check status of
! photospheric update, check status of convergence, and set specific
! control statements
  If (concon) Then
    Do j = 1, 4
      specmat(j, :) = specmat(j+1, :)
    End Do

    Do l = 1, nd
      If (err(l)>3.D-3 .And. errold(l)/=0.) Then
        specmat(5, l) = err(l)/errold(l)
      Else
        specmat(5, l) = 0.D0
      End If
    End Do
  End If

  errold = err

  frenew = .False.

  emax = 0.D0

  Do in = 1, nd
    emax = max(emax, err(in))
  End Do

  Write (999, *) ' CORR. MAX: ', emax
  Print *, ' CORR. MAX: ', emax
  precis = min(precis, emax/10.D0)

  meanerr = sum(log10(err(1:lwion)))/lwion

! JO Oct. 2025: this is the new hack: 20 iterations after temperature
! convergence, MEANERR is set at least to -2.6, to enforce ALMOST_CONVERGED =
! T
! and the damping control of the occupation numbers, to prevent oscillations

  If (optcmf_full .And. temp_converged) Then
    If (iit>=iit_temp_conv+20) meanerr = min(meanerr, -2.6D0)
  End If


! JO changed Dec. 2015/Jan. 2016
! JO further changes April 2023
  If (update_struct .And. .Not. updated) Then
    almost_converged = .False.
    Print *, ' MEAN LOGERR(1:LWION) = ', meanerr
    Write (999, *) ' MEAN LOGERR(1:LWION) = ', meanerr
    If (no_check) Then
!     UPDATED (SINCE NO_CHECK), BUT CORRESPONDING VARIABLE NOT READ (FROM
!     MODEL FILE)
      If (.Not. update_done) Then
        Write (999, *) &
          ' STOP: INCONSISTENCY BETWEEN NO_CHECK AND UPDATE_DONE(1)'
        Stop ' INCONSISTENCY BETWEEN NO_CHECK AND UPDATE_DONE(1)'
      End If
      Print *, ' PHOTOSPHERIC STRUCTURE (HYDRO) UPDATED, NO T-CORR'
      Print *
      Write (999, *) ' PHOTOSPHERIC STRUCTURE (HYDRO) UPDATED, NO T-CORR'
      Write (999, *)
    Else !                          w.r.t. no_check
      If (update_done) Then
!       UPDATED AND NO-CHECK = F (ONE ITERATION BEFORE MODEL FILE WITH
!       UPDATED=T WILL BE READ)
        Print *, ' PHOTOSPHERIC STRUCTURE (HYDRO) UPDATED, NO T-CORR'
        Write (999, *) ' PHOTOSPHERIC STRUCTURE (HYDRO) UPDATED, NO T-CORR'
      Else If (itlast==0) Then
!       EITHER BEFORE UPDATE, OR NON-SUCCESSFUL UPDATE
        Print *, ' PHOTOSPHERIC STRUCTURE (HYDRO) NOT UPDATED SO FAR'
        Write (999, *) ' PHOTOSPHERIC STRUCTURE (HYDRO) NOT UPDATED SO FAR'
      Else
!       NO INFO PRESENT, SINCE RESTART AFTER UPDATE (UPDATE_DONE REMAINS F)
        Print *, ' ITLAST NE 0, NO INFO ON PHOTOSPERIC UPDATE AVAILABLE'
        Write (999, *) ' ITLAST NE 0, NO INFO ON PHOTOSPERIC UPDATE AVAILABLE'
      End If
      Print *
      Write (999, *)
    End If !                        w.r.t. no_check
  Else !                            w.r.t.update_struct.and..not.updated
!   UPDATED (MODEL FILE READ IN)
    If (.Not. optcmf_full) Then
!     UPDATED AND NOT OPTCMF_FULL (FAST SOLUTION)
!     METCONV2 = .FALSE. REMAINS UNCHANGED
!     JO changed Jan 2016 from -2 to -2.5
      If (meanerr<=-2.5D0) Then
        almost_converged = .True.
        Print *, ' MEAN LOGERR(1:LWION) = ', meanerr, ': ALMOST_CONVERGED'
        Write (999, *) ' MEAN LOGERR(1:LWION) = ', meanerr, &
          ': ALMOST_CONVERGED'
      Else
        almost_converged = .False.
        Print *, ' MEAN LOGERR(1:LWION) = ', meanerr
        Write (999, *) ' MEAN LOG ERR(1:LWION) = ', meanerr
      End If

    Else !                          w.r.t. .not. optcmf_full
      If (.Not. optcmf_all) Then
!       UPDATED AND OPTCMF_FULL AND NOT OPTCMF_ALL (FIRST IITCMF_FULL, usually
!       0...29, ITERATIONS)
        almost_converged = .False.
        Print *, ' MEAN LOGERR(1:LWION) = ', meanerr, ' OPTCMF_ALL = .FALSE.'
        Write (999, *) ' MEAN LOGERR(1:LWION) = ', meanerr, &
          ' OPTCMF_ALL = .FALSE.'

      Else !                        w.r.t. .not. optcmf_all
!       UPDATED AND OPTCMF_FULL AND OPTCMF_ALL (AFTER BEGIN OF CMFALL))
!       JO Sept 2018
!       until further evidence, keep ALO_METALS
!       JO MARCH 2017: keep ALO always, otherwise oscillations and/or
!       lambda-iteration
!       JO Sept. 2018 completely (re-) changed, due to METCONV2
        If (meanerr<=-2.5) Then
          almost_converged = .True.
          no_alo_metals = .False.
          no_alo = .False.
          Print *, ' MEAN LOGERR(1:LWION) = ', meanerr, ': ALMOST_CONVERGED'
          Print *, ' ALOs USED IN NEXT ITERATION'
          Write (999, *) ' MEAN LOGERR(1:LWION) = ', meanerr, &
            ': ALMOST_CONVERGED'
          Write (999, *) ' ALOs USED IN NEXT ITERATION'
        Else
          almost_converged = .False.
          no_alo_metals = .False.
          no_alo = .False.
          Print *, ' MEAN LOGERR(1:LWION) = ', meanerr
          Write (999, *) ' MEAN LOGERR(1:LWION) = ', meanerr
        End If

      End If !                      w.r.t. .not.optcmf_all
    End If !                        .w.r.t. .not. optcmf_full

    If (update_struct) Then
      If (update_done) Then
        Print *, ' PHOTOSPHERIC STRUCTURE (HYDRO) UPDATED'
        Write (999, *) ' PHOTOSPHERIC STRUCTURE (HYDRO) UPDATED'
      Else If (itlast==0) Then
!       only then, info on UPDATE_DONE is present
        Print *, ' PHOTOSPHERIC STRUCTURE (HYDRO) COULD NOT BE UPDATED'
        Write (999, *) ' PHOTOSPHERIC STRUCTURE (HYDRO) COULD NOT BE UPDATED'
      Else
!       no info on UPDATE_DONE present
        Print *, ' ITLAST NE 0, NO INFO ON PHOTOSPERIC UPDATE AVAILABLE'
        Write (999, *) ' ITLAST NE 0, NO INFO ON PHOTOSPERIC UPDATE AVAILABLE'
      End If
    Else !                          .not. update_structure
      Print *, &
        ' UPDATE_STRUCTURE = F, NO UPDATE OF PHOTOSPH. STRUCTURE (HYDRO)'
      Write (999, *) &
        ' UPDATE_STRUCTURE = F, NO UPDATE OF PHOTOSPH. STRUCTURE (HYDRO)'
    End If !                        w.r.t. update_struct
    Print *
    Write (999, *)
  End If !                          outermost if clause, w.r.t.
! update_struct.and..not.updated

  If (flag_ng) Then
    Write (8, Fmt=*) iit, emax, ' ', meanerr, ' NG'
  Else
    Write (8, Fmt=*) iit, emax, ' ', meanerr
  End If

! ---- convergence output also for CONV_EFF_OPA_RAT
! NOTE: This tests the convergence *ONLY* within the range set by
! wavblue wavred. *If* we want to test complete convergence, we might
! include this in OPACITC, just before setting opa_eff_rat_old=opa_eff_rat
! ------------------------------------------------------------------
  conv_eff_opa_rat = maxval(abs((opa_eff_rat-opa_eff_rat_old)/opa_eff_rat))
  Print *, ' CORR. MAX FOR OPA_EFF_RAT: ', conv_eff_opa_rat
  Print *

  If (flag_ng_met) Then
    Write (9, Fmt='(I3,3(2X,G12.6),2X,A2)') iit, error_met, error_met_mean, &
      conv_eff_opa_rat, 'NG'
  Else
    Write (9, Fmt='(I3,3(2X,G12.6))') iit, error_met, error_met_mean, &
      conv_eff_opa_rat
  End If

  If (unasol) Go To 230 !           finalize calculation below

  If (.Not. enatcor .And. .Not. opttcor) Then

!   to prevent too early convergence
    If (optcmf_full .And. iit<=iitcmf_full) Then
      If (emax<3.D-3) emax = 3.1D-3
    End If

!   JO MARCH 2017: use MEANERR AS ADDITIONAL CONVERGENCE CRITERIUM
    If (optcmf_full) Then
      If (opt_ray .And. (emax<3.D-3 .Or. meanerr<-3.3)) unasol = .True.
!     IF (OPT_RAY .AND. (EMAX.LT.3.D-3 .OR. MEANERR.LT.-3.3)) UNASOL = .false.
    Else
      If (emax<3.D-3) unasol = .True.
    End If
  End If

! new NLTE iteration cycle
  iit = iit + 1

  lte_update = .False.

! ---- output OPA_EFF_RAT (for the total pseudo-continuum)
! ---- output created also for optically thin clumping, to check consistency
! ---- analogous to OUTTHB, for additional printouts when testing
  If (outthick) Then
    Call ncor_name(iit, ncname)
  Else
    ncname = 'LA'
  End If

  Open (1, File=trim(modnam)//'/OPA_EFF_RAT_'//trim(ncname), Status='UNKNOWN')
  Write (1, Fmt=*) nd, ifre
  Do i = 1, nd
    Do j = 1, ifre
      Write (1, 290) opa_eff_rat(i, j), 1.D8/fre(j), r(i), taur(i), i
    End Do
  End Do
  Close (1)


! UPDATE LTE-NUMBERS at begin of next iteration,
! save current NLTE numbers (old temp. struct.) for next NLTE iteration
  If (opttcor) Then
    If (optmet) Call save_metal_info
    lte_update = .True.
    Go To 130 !                     goto begin of LTE-update cycle
  End If

! NEW CYCLE
  Go To 220 !                       new NLTE cycle, without LTE-update


230 Continue

! part4--------------------------------------------------------------------
! ***finalization

! ---- output for formal solution

! if OPTCMF_ALL, then XJ = XJ_CMF_COARSE in CMF-COMPLETE range, and
! XJ = XJ(approx.) outside this range
  Call formal(xne, temp, clfac, r, velo, vmax, rtau23)
! after FORMAL, XJ has been smoothed if OPTCMF_ALL

  If (optmet) Then
    Call save_metal_info
!   JO Sept. 2021
    If (nstark>0) Then
      If (optcmf_full) Then
        Write (999, *) ' STOP: OPTCMF_FULL .AND. NSTARK > O'
        Stop ' OPTCMF_FULL .AND. NSTARK > O'
      End If
      Call save_add_freq_points
    End If
  End If

! output of problematic levels

  iq = 0
  Do l = 1, nd
    Do i = 1, id_llevs
      iq = max(iq, iqual(i,l))
    End Do
  End Do

  If (conv_eff_opa_rat>1.0D-3) Then
    Print *, ' WARNING!!! -- Maximum last relative correction in &
      &OPA_EFF_RAT > 1.d-3!'
    Print *, ' View by 2nd column in CONVERG_METALS'
    Print *
    Write (999, *) ' WARNING!!! -- Maximum last relative correction &
      &in OPA_EFF_RAT > 1.d-3!'
    Write (999, *) ' View by 2nd column in CONVERG_METALS'
    Write (999, *)
  End If

  If (enatcor .Or. temp_converged) Then
    Print *, ' TOTAL NUMBER OF APPLIED T CORRECTION(S) ', ncor
    Print *, ' MAXIMUM RELATIVE CORRECTION IN THE LAST TC ', emaxtc
    Print *
    Write (999, *) ' TOTAL NUMBER OF APPLIED T CORRECTION(S) ', ncor
    Write (999, *) ' MAXIMUM RELATIVE CORRECTION IN THE LAST TC ', emaxtc
    Write (999, *)
    Close (10)
  End If

! ---- printing also cpu time, for testings
  Call cpu_time(time_cpu)
  Print *
  Print *, ' CPU TIME: ', time_cpu
  Write (999, *)
  Write (999, *) ' CPU TIME: ', time_cpu

  If (iq==0) Then
    Print *
    Print *, ' ALL LEVELS OK!!!'
    Write (999, *)
    Write (999, *) ' ALL LEVELS OK!!!'
  Else
    Write (999, *)

    Do i = 1, id_llevs
      Do l = 1, nd
        If (iqual(i,l)/=0) Then
          istart = l
          Go To 240
        End If
      End Do
      Go To 260

240   Continue
      Do l = nd, 1, -1
        If (iqual(i,l)/=0) Then
          iend = l
          Go To 250
        End If
      End Do

      Write (999, *) ' STOP: WRONG END CONDITION'
      Stop ' WRONG END CONDITION'
250   Continue
      Write (*, Fmt=300) labl(i), i, istart, iend, r(iend)
      Write (999, Fmt=300) labl(i), i, istart, iend, r(iend)

260 End Do
  End If

! finally...

  Close (8) !                       CONVERG
  Close (9) !                       CONVERG_METALS
! 10 MAXTCORR already closed
  Close (11) !                      MAX_CHANGE
  If (optxray) Close (14) !         XLUM_ITERATION

  Close (15) !                      LTE_POT
  Close (16) !                      NLTE_POP
  Close (17) !                      SCRATCH (actual pop)
  If (optcmf_full) Close (18) !     OCCUP_NO

  Close (20) !                      MODEL
  Close (21) !                      TEMP
  Close (40) !                      SCRATCH (jnue,nf)
  Close (41) !                      SCRATCH (hnue,nf)
  Close (42) !                      SCRATCH (knue,nf)
  Close (43) !                      SCRATCH (nnue,nf)
  Close (50) !                      SCRATCH (jnue,nd)
  Close (51) !                      SCRATCH (hnue,nd)
  Close (52) !                      SCRATCH (knue,nd)
  Close (53) !                      SCRATCH (nnue,nd)
  Close (60) !                      CHIBAR_H
  Close (100) !                     megasalida (if not open, close has no
! effect)

  Print *
  Print *, '! ESTO ES EL ACABOSE !'
  Write (999, *)
  Write (999, *) '! STOP: ESTO ES EL ACABOSE !'
  Stop '! ESTO ES EL ACABOSE !'

  Close (999) !                     FASTWIND.LOG

270 Continue

  Write (999, *) ' STOP: OPTMIXED, ITLAST > 19 AND NO FILE INDEXCMF'
  Stop ' OPTMIXED, ITLAST > 19 AND NO FILE INDEXCMF'
280 Format (' IT.NO. ', I2, ' XKAP = ', E12.6, ' EX = ', F10.5, ' CORR = ', &
    F10.5)
290 Format (F15.5, E20.10, F15.5, E20.10, I10)
300 Format (' PROBLEMS WITH LEVEL ', A6, ' (NO.=', I3, ') FROM L = ', I2, &
    ' TO ', I2, '(R = ', F10.5, ')')

End Program

!***********************************************************************

!subroutines: simple ones

!***********************************************************************

Subroutine checkmass(sr, ggrav, rtau23, r, rho, nd, ndiv, ns)
! calculates photospheric mass and provides ns
! new approach, NS saved in MODEL_FILE

  Use :: nlte_type
  Use :: nlte_dim
  Use :: fund_const
  Implicit None

! note: sr is lowest radius in cgs, stellar radius is sr*rtau23

! ..
! .. scalar arguments ..
  Real (dp) :: ggrav, rtau23, sr
  Integer (i4b) :: nd
  Integer (i4b), Intent (In) :: ndiv ! TRANSITION POINT (NE 0 IF START FROM
! SCRATCH)
  Integer (i4b), Intent (In) :: ns ! TRANSITION POINT (READ IN FROM
! MODEL_FILE)


! ..
! .. array arguments ..
  Real (dp) :: r(nd), rho(nd), xm(id_ndept)
! ..
! .. local scalars ..
  Real (dp) :: dr, frac, stm1, stmass, summ, xmphot ! not used: dm, dmold
  Integer (i4b) :: i

  stmass = ggrav*(sr*rtau23)**2/gconst
  stm1 = stmass/xmsun


  summ = 0.D0
  Do i = 1, nd - 1
    dr = r(i) - r(i+1)
    summ = summ + .5D0*(rho(i)*r(i)*r(i)+rho(i+1)*r(i+1)*r(i+1))*dr
  End Do

  xm(1) = 0.D0
  Do i = 2, nd
    xm(i) = xm(i-1) + .5D0*(rho(i-1)+rho(i))*(r(i-1)-r(i))
!   for tests
!   print*,i,rho(i),r(i),log10(xm(i))
  End Do

! FROM VERSION 9.0 ON, THE OLD METHOD IS NO LONGER WORKING
! (since point close to photosphere might be close to sonic point)
! Thus: NS read from file

! RECALCULATE TRANSITION POINT
! first region with constant dlog_m
! DMOLD=LOG10(XM(ND-1)/XM(ND))
! DO I=ND-1,2,-1
! DM=LOG10(XM(I-1)/XM(I))
! print*,i,dmold,dm
! IF(ABS(1.-DMOLD/DM).GT.0.15) EXIT
! DMOLD=DM
! ENDDO
! 2nd region with constant dlog_m/3 (2 times) and dlog_m/6 (2 times)
! NS=I-4 ! cf. SUBROUTINE MODVL

  If (ndiv/=0) Then
!   CHECK CONSISTENCY IN CASE OF START FROM SCRATCH
    If (ns/=ndiv) Then
      Print *, ndiv, ' ', ns
      Write (999, *) &
        ' STOP: SOMETHING WRONG WITH TRANSITION POINT (CHECKMASS)!'
      Stop ' SOMETHING WRONG WITH TRANSITION POINT (CHECKMASS)!'
    End If
  End If

  Print *
  Print *, ' TRANSITION POINT (WIND/PHOTOSPHERE) AT L = ', ns

! note: r in units sr

  xmphot = 4.D0*pi*summ*sr**3
  frac = xmphot/stmass

  Print *
  Print *, ' STELLAR MASS = ', stm1, ' SOLAR MASSES'
  Print *, ' PHOTOSPH. MASS = ', frac, ' STELLAR MASSES'
  Print *
  Print *
  Print *, '-------------------------------------------------------------'
  Print *

  If (frac>=.01) Then
    Write (999, *) ' STOP: TOO MUCH MASS IN PHOTOSPHERE'
    Stop ' TOO MUCH MASS IN PHOTOSPHERE'
  End If

  Return

End Subroutine

!-----------------------------------------------------------------------

Subroutine copyocc(xne, nd)

  Use :: nlte_type
  Use :: nlte_dim
  Use :: nlte_var, Only: xnelte, blevel, enionnd, modnam, alevel, optmet, &
    concon

  Implicit None

! update of ne-lte to ne-nlte, copying old nlte-numbers to new ones
! changed from version 9.0 on

  Integer (i4b), Parameter :: nd1 = id_ndept, kel = id_atoms
! ..
! .. scalar arguments ..
  Integer (i4b) :: nd
! ..
! .. array arguments ..
  Real (dp) :: xne(nd1)

! .. local scalars ..
  Integer (i4b) :: i, ll, nrec
! ..
! .. local arrays ..
  Real (dp) :: xnenlte(nd1)
! ..

  nrec = id_llevs

  Open (1, File=trim(modnam)//'/ENION', Form='UNFORMATTED', Status='OLD')

  Rewind 1
  Read (1) enionnd, xnenlte
  Close (1)

  Do i = 1, nrec
    If (blevel(i)/=alevel(i)) Then
      Write (999, *) ' STOP: SOMETHING WRONG IN COPYOCC-START'
      Stop ' SOMETHING WRONG IN COPYOCC-START'
    End If

  End Do

  Print *
  Do ll = 1, nd
    Read (16, Rec=ll)(blevel(i), i=1, nrec)

    If (.Not. optmet .Or. .Not. concon) Then ! check for "new levels"
      Read (15, Rec=ll)(alevel(i), i=1, nrec) ! forgotten in versions earlier
!     than 9.0
!     NOT REQUIRED FROM VER 9.5 ON
!     above statement wrong. still required as long as ILOWIMAX_NLTE not
!     called, i.e., in very first iterations (example: CV/CVI)

      Do i = 1, nrec
        If (blevel(i)==0.) Then
!         VERY SMALL OCCUP. NUMBER
          blevel(i) = alevel(i)*1.D-10
          If (alevel(i)/=0.) Print *, 'CONCON = .FALSE., new level', ll, i
        End If
      End Do

    End If !                        always: copy nlte numbers from last nlte
!   iteration to local ones
    Write (17, Rec=ll)(blevel(i), i=1, nrec)
    xnelte(ll) = xne(ll)
    xne(ll) = xnenlte(ll)
  End Do

! print*,'used in copyocc,'
! do ll=1,nd
! print*,xnelte(ll), xne(ll)
! enddo
! print*,'xnechek: ne read in copyocc'
! print*,xnelte
! print*,xne

  Return

End Subroutine

!-----------------------------------------------------------------------

Subroutine formal(xne, temp, clf, r, velo, vmax, rtau23)

! UPDATED FOR CLUMPING, THOMSON_LINES ALREADY CORRECTED
! OPA_EFF_RAT added to output, since needed in formalsol

! if optcmf_all, then scont/xj/xk values in CONT_FORMAL, CONT_FORMAL_ALL,
! and CONT_FORMAL_CMF are smoothed!!!

  Use :: nlte_type
  Use :: nlte_dim
  Use :: fund_const, Only: sigmae, xmh => amh, clight
  Use :: nlte_var, Only: fre, ifre, opac, strue, xj, modnam, optlines, &
    thomson_lines, optcmf_all, kcmf_start, kcmf_end, etat_nolines, &
    opac_nolines, xxk, lthin
! ,RAYLEIGH
  Use :: nlte_xrays, Only: optxray, opa_kshell

  Use :: nlte_porvor, Only: opa_eff_rat_old

  Implicit None


! ---- it writes the file used in formal solution


! .. parameters ..
  Integer (i4b), Parameter :: nd = id_ndept
  Integer (i4b), Parameter :: ifretot = id_frec1
! ..
! .. array arguments ..
  Real (dp) :: temp(nd), xne(nd), clf(nd), r(nd), velo(nd)
! .. scalar arguments ..
  Real (dp) :: vmax, rtau23

! ..
  Integer (i4b) :: i, j, l, k, kk
! ..
! .. local arrays ..
  Real (dp) :: scont(nd, ifretot), thomson(nd, ifretot), efac(nd), &
    eobs(ifretot), xjlog(ifretot), sigth(nd)

  Real (dp) :: rl, mubar, vc, e1, xj1
! ..
! Just the Format of OPACITY_LP
  Logical :: out_lp = .True.
  Character (30) :: forma

! ..
  If (optcmf_all) Then

    If (xj(1,ifre+1)/=2) Then
      Write (999, *) &
        ' STOP: OPTCMF_ALL = TRUE, BUT APPROX. XJ (OBSFRAM) USED IN FORMAL'
      Stop ' OPTCMF_ALL = TRUE, BUT APPROX. XJ (OBSFRAM) USED IN FORMAL'
    End If

!   smooth xj for further use in emissivity
    Call xj_smooth(1)

!   changed July 2015 by JO

!   the following transformation is not necessary,
!   since Jnue is evaluated at CMF-frequency in FORMALSOL
!   by interpolation
!   (e.g., OPAC(nu_obs) = OPAC(nu_cmf = nu_obs*(1-mu v(r)/c))
!   thus,

    Go To 110

!   we leave this here because it might be required later on
!   transform XJ to obs. frame in CMF_COMPLETE RANGE, in analogy to
!   CMF_COMPLETE
    Do l = 1, nd
      rl = r(l)/rtau23
!     average mu to convert obs to cmf frequencies
!     assuming linear limb-darkening:
!     int(1 ... mustar) (a+b mu) dmu =: a + b mubar => mubar = 0.5(1+mustar)
      If (rl>=1.D0) Then
        mubar = 0.5D0*(1.D0+sqrt(1.D0-1.D0/rl**2))
!       in the outer wind, i_minus mostly zero
      Else
        mubar = 0.
!       in the innermost wind, i_plus and i_minus of same order
!       i_plus comes from bluer, and i_minus from redder frequencies
      End If
      vc = velo(l)*vmax/clight
      efac(l) = 1.D0 - mubar*vc
    End Do

    Do l = 1, nd
      If (efac(l)==1.D0) Exit
      eobs = fre/efac(l)
!     BUG FOUND BY JON July 2015
!     XJLOG=LOG10(XJ(L,:))
      Where (xj(l,:)>0.D0)
        xjlog = log10(xj(l,:))
      Elsewhere
        xjlog = -50.D0
      End Where

      Do k = kcmf_start, kcmf_end
        e1 = fre(k)
        Do kk = k, 2, -1
!         PRINT*,E1,EOBS(KK),EOBS(KK-1)
          If (e1<=eobs(kk) .And. e1>eobs(kk-1)) Then
            xj1 = xjlog(kk) + (xjlog(kk-1)-xjlog(kk))/(eobs(kk-1)-eobs(kk))*( &
              e1-eobs(kk))
            xj1 = 10.D0**xj1
            Go To 100
          End If
        End Do
        Write (999, *) ' STOP: NOT FOUND (INTERPOL IN FORMAL)'
        Stop ' NOT FOUND (INTERPOL IN FORMAL)'
100     Continue
        xj(l, k) = xj1
      End Do
    End Do

110 Continue

    sigth = xne*sigmae*xmh/clf !    CORRECTED
  End If


  Open (1, File=trim(modnam)//'/CONT_FORMAL', Status='UNKNOWN', &
    Form='UNFORMATTED')

  Rewind 1


  If (optlines) Then

!   if not OPTCMF_ALL, this is correct
!   else: certain inconsistency, since J_nu should be
!   from approx. treatment and observer's frame quantity,
!   but it is the CMF quantity at CMF frequency
!   => profiles using CONT_FORMAL slightly wrong (to be checked)

    Do i = 1, nd
      Do j = 1, ifre
!       ---- OPAC is effective opacity, thus ratio needs to be corrected
!       ---- Since OPAC is one iteration behind, we use OPA_EFF_RAT_OLD here
!       ---- (As long as we're converged, it of course makes no difference)
!       (XNE(I)*SIGMAE*XMH/CLF(I)+THOMSON_LINES(I,J)+RAYLEIGH(I,J)/CLF(I))/OPAC(I,J)
        thomson(i, j) = (xne(i)*sigmae*xmh/clf(i)+thomson_lines(i,j))/ &
          opac(i, j)*opa_eff_rat_old(i, j)
      End Do
    End Do
  Else

!   correct in any case, since no line treatment
    Do i = 1, nd
      Do j = 1, ifre
!       ----      As above
!       THOMSON(I,J) = (XNE(I)*SIGMAE*XMH+RAYLEIGH(I,J))/CLF(I)/OPAC(I,J)
        thomson(i, j) = (xne(i)*sigmae*xmh)/clf(i)/opac(i, j)* &
          opa_eff_rat_old(i, j)
      End Do
    End Do
  End If

  Do i = 1, nd
    Do j = 1, ifre
      scont(i, j) = strue(i, j) + thomson(i, j)*xj(i, j)
    End Do
  End Do

  Write (1) ifre, (fre(i), i=1, ifretot), ((scont(i,j),i=1,nd), j=1, ifretot), &
    ((opac(i,j),i=1,nd), j=1, ifretot), (temp(i), i=1, nd), &
    ((opa_eff_rat_old(i,j),i=1,nd), j=1, ifretot)

  Close (1)


  Open (1, File=trim(modnam)//'/CONT_FORMAL_ALL', Status='UNKNOWN', &
    Form='UNFORMATTED')

  Rewind 1

! now, all continuum quantities are stored in such a way that
! scont = strue + thomson * xj


! SAME CAVEATS AS ABOVE

  Write (1) ifre, (fre(i), i=1, ifretot), ((opac(i,j),i=1,nd), j=1, ifretot), &
    ((thomson(i,j),i=1,nd), j=1, ifretot), ((strue(i,j),i=1,nd), j=1, ifretot) &
    , ((xj(i,j),i=1,nd), j=1, ifretot), (temp(i), i=1, nd), &
    ((opa_eff_rat_old(i,j),i=1,nd), j=1, ifretot)

  Close (1)

! from Luiz
  If (optxray .And. out_lp) Then
    Open (1, File=trim(modnam)//'/OPACITY_LP', Status='UNKNOWN')
    forma = '(E13.6,F11.6,I7.3,E14.6,E14.6)'
    Write (1, *) &
      '  OPACITY      LAMBDA      dp  OPACITY-OPA_KSHELL   OPA_KSHELL'
!   WRITE (1,FMT=*) ((OPAC(I,J),I=1,ND),J=1,IFRETOT)
    Do i = 1, nd
      Do j = 1, ifre
        If (1.D8/fre(j)>0. .And. 1.D8/fre(j)<50.) Then
          Write (1, forma) opac(i, j), 1.D8/fre(j), i, &
            (opac(i,j)-opa_kshell(i,j)), opa_kshell(i, j)
        End If
      End Do
    End Do
!   WRITE (1,FMT=*) ((OPAC(I,J),I=1,ND),J=1,IFRETOT)
    Close (1)
  End If

  If (optcmf_all) Then

!   indices dividing CMF/obs. frame treatment of XJ
    Open (1, File=trim(modnam)//'/KCMF', Status='UNKNOWN', Form='FORMATTED')
    Rewind 1
    Write (1, *) kcmf_start, kcmf_end
    Close (1)

!   calculate 'new' scont from opac_nolines and etat_nolines and
!   write additional file saving theses quantities (for FORMALSOL)

!   profiles calculated with these quantities should be correct
    Open (1, File=trim(modnam)//'/CONT_FORMAL_CMF', Status='UNKNOWN', &
      Form='UNFORMATTED')

    Rewind 1

!   JO: etat_nolines, opac_nolines and sigth (from above) only corrected for
!   CLF
!   might need to be corrected for optically thick clumping
    Do i = 1, nd
      Do j = 1, ifre
        scont(i, j) = (etat_nolines(i,j)+sigth(i)*xj(i,j))/opac_nolines(i, j)
      End Do
    End Do

!   NOTE: XJ and XK (as well as OPAC and SCONT) are CMF quantities at
!   CMF-frequencies
    Write (1)((scont(i,j),i=1,nd), j=1, ifretot), ((opac_nolines(i,j),i=1,nd), &
      j=1, ifretot), ((xxk(i,j),i=1,nd), j=1, ifretot), &
      (lthin(j), j=1, ifretot)

!   JO: maybe OPA_EFF_RAT needs to be adapted and saved

    Close (1)

    Call save_line_list

  End If


  Return

End Subroutine

!-----------------------------------------------------------------------

Subroutine fresfin(nd, xnh, tmax, optout, teff)

  Use :: nlte_type
  Use :: nlte_dim
  Use :: fund_const
  Use :: princesa_var, Only: nl, nat, nions, le, li, zeff, abund, labat, labl, &
    labl4, frecin, frecfn, paren4

  Use :: nlte_var, Only: mainion, fref, wfref, indfin, indfrq, indfrq1, ifref, &
    fre, ifre, enionnd, imia, imaa, iconver, bn, index, op_flag, frecfn1

  Use :: nlte_opt, Only: optcmf_full

  Use :: nlte_xrays, Only: optxray, n_kedges, k_nmin, k_nmax, z, nionk => n, &
    eth, name

  Implicit None

! -------this subroutine calculates the fine frequency grid
! accounting for all edges. otherwise, this routine
! has to be consistent with frescal

! parameter las frecuencias del fichero de datos atomicos

! hay que definir tambien las matrices de /atomdat/

! improved for edges close to Lyman lines in Feb 2023 (JO)

! .. parameters ..
  Integer (i4b), Parameter :: kel = id_atoms, kis = id_kisat
  Integer (i4b), Parameter :: nfmin = 300
! INTEGER(I4B), PARAMETER :: LEVMIN1=ID_LEVM1
! Integer (i4b), Parameter :: ifretot = id_frec1
  Integer (i4b), Parameter :: ifrefin = id_frec2
  Integer (i4b), Parameter :: nd1 = id_ndept

! JO Dec. 2015
! FREQUENCIES FOR Lya and Lyb which should be used for consistency
! precision 1.d-9
  Real (dp), Parameter :: lya_lo = 1.D8/1217.00043262610D0, &
    lya_hi = 1.D8/1212.85710645286D0, lyb_lo = 1.D8/1027.46136261866D0, &
    lyb_hi = 1.D8/1024.50655646359D0, lyg_lo = 1.D8/973.196165978500D0, &
    lyg_hi = 1.D8/972.310779289066D0, lyd_lo = 1.D8/949.843120035809D0, &
    lyd_hi = 1.D8/948.999697008538D0
! ..
! .. scalar arguments ..
  Real (dp) :: tmax, teff
  Integer (i4b) :: nd
  Logical :: optout
! ..
! .. array arguments ..
  Real (dp) :: xnh(nd1)
! ..
! .. local scalars ..
  Real (dp) :: aaa, dm, dmin, dmin1, dmin2, dnew, dnew1, edge, emax, eps, &
    errp, frenew, pf, pfe, sum1, sum2, x, x1, xm1, xm2, xnh1, fremin, dfl
  Integer (i4b) :: i, if1, ii, ij, il, ilab, inew, inewd, io, iold, iz, j, k, &
    k1, k2, l, levp, levporg, levtot, ll, lli, m, m1, m2, mm1, mm2, n, nopa, &
    nper, nk, istart, iend, nskip

  Integer (i4b) :: nedge_ly

  Character :: lab*6, pare*6
! ..
! .. local arrays ..
  Real (dp) :: flaux(id_rbftr), freedge(id_rbftr)
  Integer (i4b) :: iopa(kis*kel)

  Integer (i4b), Dimension (50) :: arr_edge_ly

! ..
! .. external functions ..
  Real (dp) :: bnue, xintfre
  External :: bnue, xintfre
! ..
! .. external subroutines ..
  External :: sort
! ..
! .. intrinsic functions ..
! INTRINSIC ABS,DBLE,INT,LOG,LOG10,MAX,MIN
! ..
! .. statement functions ..

  Real (dp) :: epslon
! ..
! .. statement function definitions ..

  epslon(edge) = 5.D0*10.D0**(dble(int(log10(edge)))-6.D0)
! ..

  nedge_ly = 0
  arr_edge_ly = 0

  If (iconver/=1) Then
    Write (999, *) ' STOP: INCORRECT CALL OF FRESFIN'
    Stop ' INCORRECT CALL OF FRESFIN'
  End If

  If (nl>999) Then
    Write (999, *) ' STOP: TOO MANY LEVELS, CODING OF INDEX AFFECTED'
    Stop ' TOO MANY LEVELS, CODING OF INDEX AFFECTED'
  End If

! most abundant element k1

  k1 = 1
  x = abund(1)
  Do k = 2, nat
    x1 = max(x, abund(k))
    If (x1==x) Cycle
    x = x1
    k1 = k
  End Do

! 2nd abundant element if present k2

  If (nat>1) Then
    k2 = 1
    If (k1==1) k2 = 2
    x = abund(k2)
    Do k = 1, nat
      If (k==k1) Cycle
      x1 = max(x, abund(k))
      If (x1==x) Cycle
      x = x1
      k2 = k
    End Do

  End If

! dominant ionisation stages of k1, k2 (if present)

  sum1 = 0.D0
  sum2 = 0.D0

  Do l = 1, nd
    sum1 = sum1 + mainion(l, k1)
    If (nat>1) sum2 = sum2 + mainion(l, k2)
  End Do

  xm1 = sum1/nd
  xm2 = sum2/nd
  m1 = int(xm1)
  m2 = int(xm2)
  If (xm1-m1>.5D0) m1 = m1 + 1
  If (xm2-m2>.5D0) m2 = m2 + 1

  iopa(1) = k1*100 + m1

  If (nat>1) Then
    iopa(2) = k2*100 + m2
  Else
    k2 = 0
    If (m2/=0) Then
      Write (999, *) ' STOP: ERROR IN M2'
      Stop ' ERROR IN M2'
    End If
    iopa(2) = 0
  End If

  If (imaa(k1)>1) Then

!   2nd important ionization stage of k1 (if present)

    If (imaa(k1)==2) Then

!     e.g. helium

      mm1 = 3 - m1
      iopa(3) = k1*100 + mm1
    Else

      If (m1==imaa(k1)) Then
        mm1 = m1 - 1
        Go To 100
      End If

!     above or below m1?

      sum1 = 0.D0
      sum2 = 0.D0
      Do l = 1, nd
        xnh1 = 1.D0/xnh(l)
        sum1 = sum1 + enionnd(k1, m1+1, l)*xnh1
        sum2 = sum2 + enionnd(k1, m1-1, l)*xnh1
      End Do
      sum1 = sum1/dble(nd)
      sum2 = sum2/dble(nd)

      If (sum1>sum2) Then
        mm1 = m1 + 1
      Else
        mm1 = m1 - 1
      End If
100   Continue
      iopa(3) = k1*100 + mm1

    End If

    nopa = 3
  Else

!   most abundant elememt has only one ionization stage

    mm1 = m1
    nopa = 2
  End If

  If (nat>1) Then

!   2nd important ionization stage of k2, as for k1

    If (imaa(k2)>1) Then

      If (imaa(k2)==2) Then
        mm2 = 3 - m2
        iopa(nopa+1) = k2*100 + mm2
      Else
        If (m2==imaa(k2)) Then
          mm2 = m2 - 1
          Go To 110
        End If
        sum1 = 0.D0
        sum2 = 0.D0

        Do l = 1, nd
          xnh1 = 1.D0/xnh(l)
          sum1 = sum1 + enionnd(k2, m2+1, l)*xnh1
          sum2 = sum2 + enionnd(k2, m2-1, l)*xnh1
        End Do

        sum1 = sum1/dble(nd)
        sum2 = sum2/dble(nd)
        If (sum1>sum2) Then
          mm2 = m2 + 1
        Else
          mm2 = m2 - 1
        End If
110     Continue
        iopa(nopa+1) = k2*100 + mm2

      End If

      nopa = nopa + 1
    Else
      mm2 = m2

!     may happen if only one ion per element present
!     if(nopa.eq.2) stop ' error in nopa'

    End If

  End If

! total number of important ions (with levmax,levmin1)

  Do k = 1, nat
    If (abund(k)/abund(k1)<=1.D-7) Go To 120
    Do m = imia(k), imaa(k)
      If (k==k1 .And. (m==m1 .Or. m==mm1)) Cycle
      If (k==k2 .And. (m==m2 .Or. m==mm2)) Cycle
      nopa = nopa + 1
      iopa(nopa) = k*100 + m
    End Do
120 End Do

  If (nopa>kis*kel) Then
    Write (999, *) ' STOP: ERROR IN NOPA'
    Stop ' ERROR IN NOPA'
  End If

  ifref = 0
  levtot = 0

  Do n = 1, nopa

    If (iopa(n)==0) Go To 140
    k = iopa(n)/100
    i = iopa(n) - 100*k

    If (i==nions(k)) Then
      Write (999, *) ' STOP: ERROR IN NIONS OR ISI'
      Stop ' ERROR IN NIONS OR ISI'
    End If

    If (i<=0) Then
      Write (999, *) ' STOP: ERROR IN ION - FRESFIN '
      Stop ' ERROR IN ION - FRESFIN '
    End If

    levporg = 0
    Do ii = 1, id_rbftr
      ilab = labl4(ii)
      If (le(ilab)==k .And. li(ilab)==i) Then

!       JO Dec. 2015: avoid edges around Lya and Lyb
!       JO July 2021: and around Lyg and Lyd
!       JO Feb 2023: new approach: this is now allowed, will be accounted for
!       below
        If (frecin(ii)>lya_lo .And. frecin(ii)<lya_hi) Then
          Print *, ' WARNING!!!! WARNING!!!! WARNING!!!!'
          Print *, labl(ilab), &
            ' EDGE CLOSE TO LY-ALPHA NEGLECTED IN FREQ. GRID', k, ' ', i
          Write (999, *) ' WARNING!!!! WARNING!!!! WARNING!!!!'
          Write (999, *) labl(ilab), &
            ' EDGE CLOSE TO LY-ALPHA NEGLECTED IN FREQ. GRID', k, ' ', i
!         STOP ' WHEN YOU SEE THIS STATEMENT, PLEASE CONTACT J.P.'
!         after tests, the stop statement can be exchanged with the cycle
!         statement
          nedge_ly = nedge_ly + 1
          If (nedge_ly>50) Then
            Write (999, *) ' STOP: CHECK, AND INCREASE NEDGE_LY!!!'
            Stop ' CHECK, AND INCREASE NEDGE_LY!!!'
          End If
          arr_edge_ly(nedge_ly) = ii
          Cycle
        End If

        If (frecin(ii)>lyb_lo .And. frecin(ii)<lyb_hi) Then
          Print *, ' WARNING!!!! WARNING!!!! WARNING!!!!'
          Print *, labl(ilab), &
            ' EDGE CLOSE TO LY-BETA NEGLECTED IN FREQ. GRID', k, ' ', i
          Write (999, *) ' WARNING!!!! WARNING!!!! WARNING!!!!'
          Write (999, *) labl(ilab), &
            ' EDGE CLOSE TO LY-BETA NEGLECTED IN FREQ. GRID', k, ' ', i
!         STOP ' WHEN YOU SEE THIS STATEMENT, PLEASE CONTACT J.P.'
!         after tests, the stop statement can be exchanged with the cycle
!         statement
          nedge_ly = nedge_ly + 1
          If (nedge_ly>50) Then
            Write (999, *) ' STOP: CHECK, AND INCREASE NEDGE_LY!!!'
            Stop ' CHECK, AND INCREASE NEDGE_LY!!!'
          End If
          arr_edge_ly(nedge_ly) = ii
          Cycle
        End If

        If (frecin(ii)>lyg_lo .And. frecin(ii)<lyg_hi) Then
          Print *, ' WARNING!!!! WARNING!!!! WARNING!!!!'
          Print *, labl(ilab), &
            ' EDGE CLOSE TO LY-GAMMA NEGLECTED IN FREQ. GRID', k, ' ', i
          Write (999, *) ' WARNING!!!! WARNING!!!! WARNING!!!!'
          Write (999, *) labl(ilab), &
            ' EDGE CLOSE TO LY-GAMMA NEGLECTED IN FREQ. GRID', k, ' ', i
!         STOP ' WHEN YOU SEE THIS STATEMENT, PLEASE CONTACT J.P.'
!         after tests, the stop statement can be exchanged with the cycle
!         statement
          nedge_ly = nedge_ly + 1
          If (nedge_ly>50) Then
            Write (999, *) ' STOP: CHECK, AND INCREASE NEDGE_LY!!!'
            Stop ' CHECK, AND INCREASE NEDGE_LY!!!'
          End If
          arr_edge_ly(nedge_ly) = ii
          Cycle
        End If

        If (frecin(ii)>lyd_lo .And. frecin(ii)<lyd_hi) Then
          Print *, ' WARNING!!!! WARNING!!!! WARNING!!!!'
          Print *, labl(ilab), &
            ' EDGE CLOSE TO LY-DELTA NEGLECTED IN FREQ. GRID', k, ' ', i
          Write (999, *) ' WARNING!!!! WARNING!!!! WARNING!!!!'
          Write (999, *) labl(ilab), &
            ' EDGE CLOSE TO LY-DELTA NEGLECTED IN FREQ. GRID', k, ' ', i
!         STOP ' WHEN YOU SEE THIS STATEMENT, PLEASE CONTACT J.P.'
!         after tests, the stop statement can be exchanged with the cycle
!         statement
          nedge_ly = nedge_ly + 1
          If (nedge_ly>50) Then
            Write (999, *) ' STOP: CHECK, AND INCREASE NEDGE_LY!!!'
            Stop ' CHECK, AND INCREASE NEDGE_LY!!!'
          End If
          arr_edge_ly(nedge_ly) = ii
          Cycle
        End If

        levporg = levporg + 1
        flaux(levporg) = frecin(ii)

      End If
    End Do

!   reordering of flaux since fl maybe not ordered by energy
!   example: hei

    levtot = levtot + levporg
    Call sort(levporg, flaux)

!   ---  all levels, up to the highest ones consistent with frescal

    levp = levporg
    If (levp==0) Then
      Write (999, *) ' STOP: ERROR IN LEVP - FRESMIN '
      Stop ' ERROR IN LEVP - FRESMIN '
    End If

!   levp1 = min(levp,levmin1)

    Do l = 1, levp
      ifref = ifref + 1
      lli = levporg + 1 - l
      Do ii = 1, id_rbftr
        If (flaux(lli)==frecin(ii)) Go To 130
      End Do
      Write (999, *) ' STOP: FLAUX NOT FOUND IN FRECIN'
      Stop ' FLAUX NOT FOUND IN FRECIN'

130   Continue

      index(ifref) = 100000*k + i*1000 + ii

!     note that flaux here is in the "wrong' order

      fref(ifref) = flaux(lli)

    End Do
140 End Do

  If (levtot>id_rbftr) Then
    Write (999, *) ' STOP ERROR IN LEVTOT'
    Stop ' ERROR IN LEVTOT'
  End If

! IF (IFREF.GT.ID_LLEVS) STOP  ! this was the old statement with one edge per
! level
  If (ifref>id_rbftr) Then
    Write (999, *) ' STOP: ERROR IN IFREF' ! THIS IS THE NEW ONE
    Stop ' ERROR IN IFREF' !        THIS IS THE NEW ONE
  End If

  Print *, ' NUMBER OF CONSIDERED EDGES = ', ifref
  Print *

  Do i = 1, ifref
    eps = epslon(fref(i))
    fref(i+ifref) = fref(i) - eps
    freedge(i) = fref(i)
  End Do

  if1 = ifref
  ifref = 2*ifref

! K-shell edges
  nk = 0

  If (optxray) Then
!   only edges from ionization stage III on
    Do k = 1, n_kedges
      If (nionk(k)<k_nmin) Cycle
      If (nionk(k)>k_nmax) Cycle
      If (eth(k)>2.D7) Cycle
      nk = nk + 1
      ifref = ifref + 1
      fref(ifref) = eth(k)
      ifref = ifref + 1
      eps = epslon(eth(k))
      fref(ifref) = fref(ifref-1) - eps
    End Do
    Print *, ' NUMBER OF CONSIDERED K-SHELL EDGES = ', nk
    Print *
  End If

  Call sort(ifref, fref)

  fremin = fref(1)
  i = log10(fremin) - 1
  i = max(i, 0)
  i = 10**i
  fremin = float(int(fremin/i)*i)

  If (fref(1)<=fremin) Then
    Write (999, *) ' STOP: ERROR IN FREMIN'
    Stop ' ERROR IN FREMIN'
  End If

! OLD VERSION, WITH FRE(1)=FREMIN
! DO I = IFREF + 1,2,-1
! FREF(I) = FREF(I-1)
! END DO

! FREF(1) = FREMIN
! IFREF = IFREF + 1

! NEW VERSION FOR mm-FLUXES
! 9 points before fremin, starting at 1mm = 10 Kayser

  dfl = log10(fremin/10.)/9. !      10 Kayser, 10 points = 9 intervals

  Do i = ifref + 10, 11, -1
    fref(i) = fref(i-10)
  End Do

  fref(1) = 10.
  Do i = 1, 9
    fref(i+1) = fref(1)*10.**(i*dfl)
  End Do

  ifref = ifref + 10

  If (optxray) Then
    aaa = 2.D7 !                    (5 A)

  Else
    emax = 10.D0*tmax/hkl

    If (nat==1 .And. k1==1 .Or. teff<10000.) Then
      aaa = 4.D5 !                  (250 A)
    Else If (teff<20000.) Then
      aaa = 1.D8/180.D0 !           (180 A)
    Else If (teff<35000.) Then
      aaa = 1.D6 !                  (100 A)
    Else
      aaa = 5.D6 !                  ( 20 A)
    End If

  End If

  Print *
  Print *, ' INTERMEDIATE FREF(IFREF): ', 1.0D8/fref(ifref)
  Print *, ' MIN      : ', 1.0D8/aaa
  Print *


  emax = max(emax, aaa)
  emax = max(emax, fref(ifref)*1.2D0)

  If (emax>7.5D5) Then !            in case, resolve HeII edge
    ifref = ifref + 1
    fref(ifref) = 7.5D5 !           (133 A)
  End If


  If (optxray) Then
    ifref = ifref + 1
    fref(ifref) = 1.2D6 !           (86 A, half way between 133 A and first
!   K-shell at 39A)
    ifref = ifref + 1
    fref(ifref) = 5.D6 !            ( 20 A)
  End If

  ifref = ifref + 1
  fref(ifref) = emax

! here was a bug, missing sort statement until Oct. 2014
  Call sort(ifref, fref)

! changed to obtain similar freq. grids in all cases
! dmin=log(fre(ifre)/fre(1))/nfmin

  dmin = log(1.D6/1.D3)/nfmin
  nper = ifref - 1
  Do i = 1, nper

    If (fref(i+1)<1.D4) Then !      >10000 A
      dmin1 = dmin*3
    Else If (fref(i+1)<5.D4) Then ! >2000 A
      dmin1 = dmin*2
    Else If (fref(i+1)<1.D8/1600.) Then ! >1600 A
!     changed Feb 2015, for better resol between 2000 and 1600 of non-HHe
!     models
!     DMIN1 = DMIN
      dmin1 = dmin/4.
!     here was the  bug (corrected from V10.1 on): fre instead of fref)
    Else If (fref(i+1)<1.D8/910.) Then ! > 910 A
      dmin1 = dmin/4.

!     minimum separation behind HeII edge should be approx. 10  A or 2 A (if
!     X-rays)

    Else If (fref(i+1)>440528.D0) Then ! < 227 A
      dmin2 = 10.D-8*fref(i+1) !    DELTA LAM / LAM
      If (optxray) dmin2 = 2.D-8*fref(i+1) ! DELTA LAM / LAM
      If (optxray .And. fref(i+1)>5.D6) dmin2 = 0.3
      dmin1 = dmin2
    Else !                          227 (or min) ... 910 A
      dmin1 = 2000./300000. !       (MIN. RESOL = 2000 KM/S)
    End If

    dm = fref(i+1)/fref(i) - 1.D0

!   nice trick to obtain resolved edges in any case

    If (dm>dmin1/4.D0) Then
      inew = int(dm/dmin1) + 1
      dnew = (fref(i+1)-fref(i))/inew

      Do j = 1, inew
        dnew1 = dnew

!       concentration towards edges for all transitions

        If (j==1) Then
          dnew1 = dnew/4.D0
          Do k = 1, 3
            frenew = fref(i) + k*dnew1
            ifref = ifref + 1
            fref(ifref) = frenew
          End Do

        Else If (j==2) Then
!         Jo June 2025: this block has changed. Also here, we use the
!         same stepsize (DNEW/4.) as for J=1, in case of OPTCMF_FULL AND .NOT.
!         OPTXRAY
!         for edges below the HeII one, to obtain equidistant grid-points
!         which allows for a better interpolation of xj_cmf_coarse.
          If (optcmf_full .And. .Not. optxray .And. fref(i+1)>440528.D0) Then
!           < 227 A
            dnew1 = dnew/4.D0
            Do k = 1, 3
              frenew = fref(ifref) + dnew1
              ifref = ifref + 1
              fref(ifref) = frenew
            End Do
!           old version
          Else
            dnew1 = dnew/2.D0
            frenew = fref(ifref) + dnew1
            ifref = ifref + 1
            fref(ifref) = frenew
          End If
        End If

        If (j/=inew) Then
          frenew = fref(i) + j*dnew
          ifref = ifref + 1
          fref(ifref) = frenew
        End If
      End Do

    End If
  End Do

! JO Dec. 2015
! modify frequency points around Lya and Lyb so that consistency
! with pure HHe freq. grid (if HHe only, no action performed)
! points chosen in such a way that Lya and Lyb not completely centered,
! to allow for a compromise (otherwise self-shadowing too strong)

! in case, add points
  Do i = 1, ifref
    If (abs(fref(i)-lya_lo)<1.D-9) Go To 150
  End Do
  ifref = ifref + 1
  fref(ifref) = lya_lo
  If (optout .And. iconver==1) Write (*, Fmt=340) 1.D8/fref(ifref)

150 Do i = 1, ifref
    If (abs(fref(i)-lya_hi)<1.D-9) Go To 160
  End Do
  ifref = ifref + 1
  fref(ifref) = lya_hi
  If (optout .And. iconver==1) Write (*, Fmt=340) 1.D8/fref(ifref)

160 Do i = 1, ifref
    If (abs(fref(i)-lyb_lo)<1.D-9) Go To 170
  End Do
  ifref = ifref + 1
  fref(ifref) = lyb_lo
  If (optout .And. iconver==1) Write (*, Fmt=340) 1.D8/fref(ifref)

170 Do i = 1, ifref
    If (abs(fref(i)-lyb_hi)<1.D-9) Go To 180
  End Do
  ifref = ifref + 1
  fref(ifref) = lyb_hi
  If (optout .And. iconver==1) Write (*, Fmt=340) 1.D8/fref(ifref)

180 Continue

  Call sort(ifref, fref)

! -----------------------------------------------
! in case, remove points in between lya_lo,lya_hi

  Do i = 1, ifref
    If (abs(fref(i)-lya_lo)<1.D-9) Go To 190
  End Do
  Write (999, *) ' STOP: LYA_LO NOT FOUND IN FRESFIN'
  Stop ' LYA_LO NOT FOUND IN FRESFIN'


190 istart = i
  Do i = istart, ifref
    If (abs(fref(i)-lya_hi)<1.D-9) Go To 200
  End Do
  Write (999, *) ' STOP: LYA_HI NOT FOUND IN FRESFIN'
  Stop ' LYA_HI NOT FOUND IN FRESFIN'

200 iend = i

  If (optout .And. iconver==1) Then
    Do i = istart + 1, iend - 1
      Write (*, Fmt=350) 1.D8/fref(i)
    End Do
  End If

  nskip = iend - istart - 1
  If (nskip<0) Then
    Write (999, *) ' STOP: ERROR IN NSKIP (FRESFIN)'
    Stop ' ERROR IN NSKIP (FRESFIN)'
  End If

  If (nskip/=0) Then
    ifref = ifref - nskip
    Do i = istart + 1, ifref
      fref(i) = fref(i+nskip)
    End Do
  End If

! -----------------------------------------------
! in case, remove points in between lyb_lo,lyb_hi

  Do i = 1, ifref
    If (abs(fref(i)-lyb_lo)<1.D-9) Go To 210
  End Do
  Write (999, *) ' STOP: LYB_LO NOT FOUND IN FRESFIN'
  Stop ' LYB_LO NOT FOUND IN FRESFIN'

210 istart = i
  Do i = istart, ifref
    If (abs(fref(i)-lyb_hi)<1.D-9) Go To 220
  End Do
  Write (999, *) ' STOP: LYB_HI NOT FOUND IN FRESFIN'
  Stop ' LYB_HI NOT FOUND IN FRESFIN'

220 iend = i

  If (optout .And. iconver==1) Then
    Do i = istart + 1, iend - 1
      Write (*, Fmt=350) 1.D8/fref(i)
    End Do
  End If

  nskip = iend - istart - 1
  If (nskip<0) Then
    Write (999, *) ' STOP: ERROR IN NSKIP (FRESFIN)'
    Stop ' ERROR IN NSKIP (FRESFIN)'
  End If

  If (nskip/=0) Then
    ifref = ifref - nskip
    Do i = istart + 1, ifref
      fref(i) = fref(i+nskip)
    End Do
  End If

! include high resolution points from frescal
! (hydrogen resonance lines (in particular, Lyg and Lyd), hei singlet
! resonance line and heii 303)

outer: Do i = 1, ifre
    Do ij = 1, ifref
      If (abs(1.-fref(ij)/fre(i))<1.D-14) Cycle outer ! CAN BE EQUAL BY CHANCE
    End Do
    ifref = ifref + 1
    fref(ifref) = fre(i)
    Write (*, Fmt=330) 1.D8/fref(ifref)
  End Do outer


  Call sort(ifref, fref)

  If (ifref>ifrefin) Then
    Write (999, *) ' STOP: TOO MANY FREQUENCIES IN FRESFIN!'
    Stop ' TOO MANY FREQUENCIES IN FRESFIN!'
  End If

  Write (999, *)
  Write (999, *) ' TOTAL NUMBER OF FINE GRID FREQUENCIES = ', ifref
  Write (999, *)
  Write (999, *)
  Write (999, *) ' FINAL FREF(1): ', 1.0D8/fref(1)
  Write (999, *) ' FINAL FREF(IFREF): ', 1.0D8/fref(ifref)

  Print *
  Print *, ' TOTAL NUMBER OF FINE GRID FREQUENCIES = ', ifref
  Print *
  Print *
  Print *, ' FINAL FREF(IFREF): ', 1.0D8/fref(ifref)


  If (optout .And. iconver==1) Then

!   WRITE(999,*) ' NO    WAVELENGTH(A)     FREQ(HZ)'
    Print *, ' NO    WAVELENGTH(A)     FREQ(HZ)'
    Do i = 1, ifref
      Do j = 1, if1
        If (fref(i)==freedge(j)) Go To 230
      End Do

      Do ij = 1, ifre
        If (abs(1.-fref(i)/fre(ij))<1.D-14) Then
          indfin(i) = ij

          If (optxray) Then
            Do j = 1, n_kedges
              If (fref(i)==eth(j)) Then
!               WRITE (999,FMT=9025) INDFIN(I),1.D8/FREF(I), &
!               CLIGHT*FREF(I),NAME(Z(J)),NIONK(J)-1
                Write (*, Fmt=310) indfin(i), 1.D8/fref(i), clight*fref(i), &
                  name(z(j)), nionk(j) - 1
                Go To 240
              End If
            End Do
          End If

!         WRITE (999,FMT=9000) INDFIN(I),1.D8/FREF(I), &
!         CLIGHT*FREF(I)
          Write (*, Fmt=280) indfin(i), 1.D8/fref(i), clight*fref(i)
          Go To 240
        End If
      End Do

      If (i==1) Then
        Write (999, *) ' STOP: ERROR IN PHILOSOPHY - FRESFIN '
        Stop ' ERROR IN PHILOSOPHY - FRESFIN '
      End If
      indfin(i) = -abs(indfin(i-1))
!     WRITE (999,FMT=9010) INDFIN(I),1.D8/FREF(I),CLIGHT*FREF(I)
      Write (*, Fmt=290) indfin(i), 1.D8/fref(i), clight*fref(i)
      Go To 240

230   Continue
      il = index(j)
      k = il/100000
      iz = int(zeff(k))
      il = il - 100000*k
      io = il/1000
      ll = il - 1000*io
      lab = labl(labl4(ll))
      pare = paren4(ll)
      Do ij = 1, ifre
        If (fref(i)==fre(ij)) Then
          indfin(i) = ij
!         WRITE (999,FMT=9020) INDFIN(I),1.D0/FREF(I)*1.D8, &
!         CLIGHT*FREF(I),LABAT(K),IO,IZ + IO - 1,LAB,PARE
          Write (*, Fmt=300) indfin(i), 1.D0/fref(i)*1.D8, clight*fref(i), &
            labat(k), io, iz + io - 1, lab, pare
          Go To 240
        End If
      End Do

      If (i==1) Then
        Write (999, *) ' STOP: ERROR IN PHILOSOPHY - FRESFIN '
        Stop ' ERROR IN PHILOSOPHY - FRESFIN '
      End If
      indfin(i) = -abs(indfin(i-1))
!     WRITE (999,FMT=9030) 1.D0/FREF(I)*1.D8,CLIGHT*FREF(I), &
!     LABAT(K),IO,IZ + IO - 1,LAB,PARE
      Write (*, Fmt=320) 1.D0/fref(i)*1.D8, clight*fref(i), labat(k), io, &
        iz + io - 1, lab, pare

240 End Do

    Print *

  End If

! ---    check of indfin

  iold = indfin(1)

  Do i = 2, ifref
    inew = abs(indfin(i))
    inewd = inew - iold
    If (inewd/=0 .And. inewd/=1) Then
      Print *, ifref, i, inew, iold, inewd
      Write (999, *) ' STOP: ERROR IN INDFIN'
      Stop ' ERROR IN INDFIN'
    End If
    iold = inew
  End Do

! ---   final check

  If (abs(indfin(ifref))/=ifre) Then
    Write (999, *) ' STOP: ERROR IN INDFIN: NUMBER OF FREQ. POINTS'
    Stop ' ERROR IN INDFIN: NUMBER OF FREQ. POINTS'
  End If

! index array with indices of edge frequencies
! according to input (i.e., ground or excited states)

  indfrq = 0
  indfrq1 = 0

  Do ii = 1, id_rbftr

    Do ij = 1, ifref
      If (fref(ij)==frecin(ii)) Then
        indfrq(ii) = ij
        Go To 250
      End If
    End Do

!   may happen if not all ions present
!   stop ' edge not found in fref'

250 End Do

! JO Feb 2023: couple of changes, to allow edges close to Lyman lines

! count edged that are not included in FRECIN (close to Lyman, or if not
! all ions present)
  ij = 0
  Do ii = 1, id_rbftr
    If (indfrq(ii)==0) ij = ij + 1
  End Do

! lower than, since additional edges might be absent if not all ions present
! minimum IJ=NEDGE_LY
  If (ij<nedge_ly) Then
    Print *, ij, nedge_ly
    Write (999, *) ' STOP: SOMETHING ROTTEN WITH INDFRQ, 1'
    Stop ' SOMETHING ROTTEN WITH INDFRQ, 1'
  End If

  If (nedge_ly/=0) Then
    Write (999, *) ' EDGES CLOSE TO LYMAN LINES FOUND. NUMBER OF EDGES = ', &
      nedge_ly
    Write (999, *) ' LOWER LEVELS AND WAVELENGTHS ARE:'
    Print *, ' EDGES CLOSE TO LYMAN LINES FOUND. NUMBER OF EDGES = ', nedge_ly
    Print *, ' LOWER LEVELS AND WAVELENGTHS ARE:'
    Do ij = 1, nedge_ly
      ii = arr_edge_ly(ij)
      Write (999, *) ' ', labl(labl4(ii)), ' ', 1.D8/frecin(ii)
      Print *, ' ', labl(labl4(ii)), ' ', 1.D8/frecin(ii)
      If (indfrq(ii)/=0) Then
        Write (999, *) ' STOP: SOMETHING ROTTEN WITH INDFRQ, 2'
        Stop ' SOMETHING ROTTEN WITH INDFRQ, 2'
      End If
    End Do
    Write (999, *)
    Print *
  End If

  Do i = 1, nedge_ly
    ii = arr_edge_ly(i)
    Do ij = 1, ifref - 1
      If (fref(ij)<frecin(ii) .And. fref(ij+1)>=frecin(ii)) Then
        indfrq(ii) = ij + 1
!       INDFIN SHOULD BE POSITIVE, SINCE PROBLEMATIC EDGES BETWEEN COARSE GRID
!       FREQUENCY POINTS
        Print *, ' EDGE FOR ', labl(labl4(ii)), ' SHIFTED TO ', &
          1.D8/fref(indfrq(ii))
        If (indfin(ij+1)<=0) Then
          Write (999, *) ' STOP: INDFIN <= 0 FOR EDGE CLOSE TO LYMAN LINES'
          Stop ' INDFIN <= 0 FOR EDGE CLOSE TO LYMAN LINES'
        End If
        Go To 260
      End If
    End Do
    Write (999, *) ' STOP: EDGE CLOSE TO LYMAN LINES NOT FOUND IN FREF'
    Stop ' EDGE CLOSE TO LYMAN LINES NOT FOUND IN FREF'
260 End Do

! index array with indices of edge frequencies for OP-data
! (range between ground and excited states)


! JO Feb. 2023, checked: the additional path with INDFRQ1 is NOT affected
! by any edge in the Lyman range, since the corresponding FRECFN1 is adapted
! from the coarse freq. grid, which has been already modified.
  Do ii = 1, id_rbftr

    If (.Not. op_flag(ii)) Then
      indfrq1(ii) = 0
      Cycle
    Else
      Do ij = 1, ifref
!       JO Sept. 2016: changed to account for FRECFN1 (OP-DATA)
        If (frecfn1(ii)/=frecfn(ii)) Then
          If (fref(ij)==frecfn1(ii)) Then
            indfrq1(ii) = ij
            If (indfin(ij)<=0) Then
              Write (999, *) ' STOP: INDFIN(IJ) <= 0 FOR INDFRQ1 (PATH1)'
              Stop ' INDFIN(IJ) <= 0 FOR INDFRQ1 (PATH1)'
            End If
!           PRINT*,'1',II,LABL(LABL4(II)),FRECFN1(II),INDFIN(IJ),FREF(IJ)
            Go To 270
          End If
        Else
          If (fref(ij)>=frecfn1(ii) .And. indfin(ij)>0) Then
!           HERE: GE, since in this case FRECFN1 (=FRECFN) not necessarily a
!           frequency grid point;
!           (it is definitely a grid point if either FRECFN1 > FRECFN --
!           previous
!           condition or FRECFN1=FRECFN=FRECIN -- this condition)
!           INDFIN has to be positive, since edge needs to be on the coarse
!           grid
            indfrq1(ii) = ij
!           PRINT*,'2',II,LABL(LABL4(II)),FRECFN1(II),INDFIN(IJ),FREF(IJ)
            Go To 270
          End If
        End If
      End Do
!     PRINT*,II,LABL(LABL4(II)),OP_FLAG(II)
!     STOP 'NOT FOUND'

!     may happen if not all ions present
!     stop ' edge not found in fre'

    End If

270 End Do

  Do k = 1, ifref
    bn(k) = bnue(1.D8/fref(k), tmax)
  End Do

  pf = xintfre(bn, wfref, 1, ifref, ifref, fref, 4, 'NEW')
  pfe = sigsb/pi*tmax**4
  errp = abs(1.D0-pf/pfe)
  Print *, ' ERROR IN PLANCK-FUNCTION WITH FINE FREQ. GRID = ', errp
  Print *

  Return

280 Format ('*', I4, 2X, F12.3, 6X, E12.6)
290 Format (I5, 2X, F12.3, 6X, E12.6)
300 Format ('*', I4, 2X, F12.3, 6X, E12.6, 3X, A2, I2, ' Z=', I2, '  FROM ', &
    A6, ' TO ', A6)
310 Format ('*', I4, 2X, F12.3, 6X, E12.6, 3X, A2, ' Z=', I2, &
    ' K-SHELL ABSORPTION')
320 Format (I5, 2X, F12.3, 6X, E12.6, 3X, A2, I2, ' Z=', I2, '  FROM ', A6, &
    ' TO ', A6)
330 Format (' ADDIT. POINT FROM FRESCAL AT ', F12.3)
340 Format (' ADDIT. POINT AT ', F12.3, ' LY_ALPHA/LY_BETA')
350 Format ('        POINT AT ', F12.3, ' DROPPED (LY_ALPHA/LY_BETA)')
End Subroutine

!-----------------------------------------------------------------------

Subroutine hopfpara(teff, ggrav, yhe, xmet, xmloss, vmax, srnom, clf)

  Use :: nlte_type
  Use :: nlte_dim
  Use :: nlte_var, Only: qinf, q0, gam => gamhopf

  Implicit None

! assumes files 'hopfpara_all_HHe/met' in catalog above local directory

! .. scalar arguments ..
  Real (dp) :: ggrav, teff, yhe, xmet, xmloss, vmax, srnom, clf
! ..
! .. local scalars ..
  Real (dp) :: gin, tin, yin, qqin, qqq
! ..
! ..

! INCLUDING MAX. CLF
  qqq = log10(sqrt(clf)*xmloss/(vmax*srnom)**1.5)
  Print *, 'log Q = (for Hopfpara-file)', qqq

  If (abs(xmet)<0.2D0) Then
    Open (3, File='../HOPFPARA_ALL_HHe', Status='OLD')
  Else
    Open (3, File='../HOPFPARA_ALL_met', Status='OLD')
  End If

! first line is description

  Read (3, Fmt=*)

readfile: Do

!   READ (1,FMT=*,END=20) TIN,GIN,QQIN,YIN,QINF,Q0,GAM
!   IF (TIN.EQ.TEFF .AND. GIN.EQ.GGRAV .AND. YIN.EQ.YHE .AND. &
!   &    abs(QQIN-QQQ).lt.1.d-3) THEN

    Read (3, Fmt=*, End=100) tin, gin, yin, qinf, q0, gam
    If (tin==teff .And. gin==ggrav .And. yin==yhe) Then
      Print *
      Print *, ' QINF = ', qinf, ' Q0 = ', q0, ' GAMMA = ', gam
      Print *
      Close (3)
      Return
    End If

  End Do readfile

100 Close (3)


  If (abs(xmet)<0.2D0) Then
    Open (3, File='../HOPFPARA_ALL_HHe', Status='OLD')
  Else
    Open (3, File='../HOPFPARA_ALL_met', Status='OLD')
  End If

  Read (3, Fmt=*)

  qqq = -14.

readfile_1: Do

    Read (3, Fmt=*, End=110) tin, gin, qqin, yin, qinf, q0, gam

    If (tin==teff .And. gin==ggrav .And. yin==yhe .And. qqin==qqq) Then
      Print *
      Print *, ' QINF = ', qinf, ' Q0 = ', q0, ' GAMMA = ', gam
      Print *
      Close (3)
      Return
    End If

  End Do readfile_1

110 Continue
  Write (999, *) ' STOP: HOPF-PARAMETERS NOT FOUND'
  Stop ' HOPF-PARAMETERS NOT FOUND'

End Subroutine

!-----------------------------------------------------------------------

Subroutine pzgrid(nd, np, r)

  Use :: nlte_type
  Use :: nlte_dim
  Use :: nlte_var, Only: p, z
  Implicit None

! .. parameters ..
  Integer (i4b), Parameter :: nd1 = id_ndept, np1 = id_npoin
! ..
! .. scalar arguments ..
  Integer (i4b) :: nd, np
! ..
! .. array arguments ..
  Real (dp) :: r(nd1)
! ..
! .. local scalars ..
  Real (dp) :: z2
  Integer (i4b) :: j, l, lmax, nc
! ..
! .. intrinsic functions ..
! INTRINSIC DBLE,MIN0,SQRT
! ..

  If (nd/=nd1 .Or. np/=np1) Then
    Write (999, *) ' STOP: WRONG DIMENSIONS IN PZGRID'
    Stop ' WRONG DIMENSIONS IN PZGRID'
  End If

  nc = np - nd

! calculation of the p-grid

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
        If (z2>1.D-5) Then
          Write (999, *) ' STOP: ERROR IN Z=0'
          Stop ' ERROR IN Z=0'
        End If
        z2 = 0.
      End If
      z(j, l) = sqrt(z2)
    End Do
  End Do

  Return

End Subroutine

!-----------------------------------------------------------------------

Subroutine start(nrec, teff, yhein, optmet)

! JO Sept 2023: additional call of start_expl_lines

  Use :: nlte_type
  Use :: nlte_dim
  Use :: fund_const, Only: akb, amh
  Use :: princesa_var, Only: nat, labat, abund, weight, frecfn
  Use :: nlte_var, Only: vturb, vther, nat_tot, names_changed, abchanged, &
    nat_bg, names_bg, abundbg, frecfn1

  Implicit None

! ----- read atomic data

! .. scalar arguments ..
  Real (dp) :: teff, yhein
  Integer (i4b) :: nrec
  Logical :: optmet
! ..
! .. local scalars ..
  Real (dp) :: yhe
  Integer (i4b) :: i, isi, nene
! ..
! .. external subroutines ..
  External :: openfs, prince
! ..
! .. intrinsic functions ..
! INTRINSIC SQRT
! ..

  Call prince
! check for identical rbf frequencies and modify them
  Call rbf_check
! JO Sept 2016 initialize FRECFN1
  frecfn1 = frecfn


  nene = 0
  Do isi = 1, nat
    If (labat(isi)=='HE') nene = isi
  End Do

  If (nene/=0) Then
    yhe = abund(nene)
    If (abs(1.D0-yhe/yhein)>1.D-5) Then
      Print *
      Print *, ' ********************************************'
      Print *, ' ** WARNING: YHE AND YHEIN ARE DIFFERENT  ***'
      Print *, ' ** MODEL ATOM HAS YHE= ', yhe, ' INPUT IS ', yhein
      Print *, ' ** ADOPTED : ', yhein
      Print *, ' ********************************************'
      Print *
      Write (999, *)
      Write (999, *) ' ********************************************'
      Write (999, *) ' ** WARNING: YHE AND YHEIN ARE DIFFERENT  ***'
      Write (999, *) ' ** MODEL ATOM HAS YHE= ', yhe, ' INPUT IS ', yhein
      Write (999, *) ' ** ADOPTED : ', yhein
      Write (999, *) ' ********************************************'
      Write (999, *)
      abund(nene) = yhein
    End If
  End If
! prespecified abundances of explicit and background elements
  nat_bg = 0

  Do i = 1, nat_tot
    If (names_changed(i)=='HE') Then
      Write (999, *) ' STOP: HE FOUND IN CHANGED ABUNDANCES'
      Stop ' HE FOUND IN CHANGED ABUNDANCES'
    End If

    nene = 0
    Do isi = 1, nat
      If (labat(isi)==names_changed(i)) nene = isi
    End Do

    If (nene/=0) Then
      Print *
      Print *, ' ********************************************'
      Print *, ' ** EXPLICIT ABUNDANCE CHANGED ***'
      Print *, ' ** MODEL ATOM ', labat(nene), ' CHANGED TO', abchanged(i)
      Print *, ' ********************************************'
      Print *
      Write (999, *)
      Write (999, *) ' ********************************************'
      Write (999, *) ' ** EXPLICIT ABUNDANCE CHANGED ***'
      Write (999, *) ' ** MODEL ATOM ', labat(nene), ' CHANGED TO', &
        abchanged(i)
      Write (999, *) ' ********************************************'
      Write (999, *)
      abund(nene) = abchanged(i)
    Else
!     must be background element
      nat_bg = nat_bg + 1
      abundbg(nat_bg) = abchanged(i)
      names_bg(nat_bg) = names_changed(i)
      Print *
      Print *, ' ********************************************'
      Print *, ' ** BACKGROUND ABUNDANCE CHANGED ***'
      Print *, ' ** MODEL ATOM ', names_bg(nat_bg), ' CHANGED TO', &
        abundbg(nat_bg)
      Print *, ' ********************************************'
      Print *
      Write (999, *)
      Write (999, *) ' ********************************************'
      Write (999, *) ' ** BACKGROUND ABUNDANCE CHANGED ***'
      Write (999, *) ' ** MODEL ATOM ', names_bg(nat_bg), ' CHANGED TO', &
        abundbg(nat_bg)
      Write (999, *) ' ********************************************'
      Write (999, *)
    End If
  End Do

  Do i = 1, nat

!   finally changed

    vther(i) = sqrt(2.D0*akb*teff/weight(i)/amh+vturb**2)
  End Do

! -----open direct access-files

  Call openfs(nrec)

! new call from Sept 2023 on

  If (optmet) Call start_expl_lines(teff)

  Return

End Subroutine

!-----------------------------------------------------------------------

Subroutine openfs(nrecet)

  Use :: nlte_type
  Use :: nlte_dim
  Use :: nlte_var, Only: modnam
  Implicit None

! ------ this subroutine opens all direct access (and some other files.)
! the files have a status of 'scratch', so they will be dismissed
! when the program ends

! .. parameters ..
  Integer (i4b), Parameter :: nd = id_ndept
  Integer (i4b), Parameter :: ifretot = id_frec1
! ..
! .. scalar arguments ..
  Integer (i4b) :: nrecet
! ..
! .. local scalars ..
  Integer (i4b) :: infre, nde, nrec
! ..
  nrec = nrecet
  nrec = nrec*8
  nde = 8*nd

  infre = ifretot*8

! -----------------------------------------------------------------------
! info file for changes in occupation numbers and extrapolations
! (acceleration)

  Open (Unit=11, File=trim(modnam)//'/MAX_CHANGE', Status='UNKNOWN')

! file for lte population (alevel(i))

  Open (Unit=15, File=trim(modnam)//'/LTE_POP', Status='UNKNOWN', Recl=nrec, &
    Access='DIRECT')

! -----------------------------------------------------------------------
! file for occnlte_old (blevel(i))

  Open (Unit=16, File=trim(modnam)//'/NLTE_POP', Status='UNKNOWN', Recl=nrec, &
    Access='DIRECT')

! -----------------------------------------------------------------------
! file for actual population (blevel(i))

  Open (Unit=17, Status='SCRATCH', Recl=nrec, Access='DIRECT')

! -----------------------------------------------------------------------
! file 'model'

  Open (Unit=20, File=trim(modnam)//'/MODEL', Status='UNKNOWN', &
    Form='UNFORMATTED')

! -----------------------------------------------------------------------
! file 'temp'

  Open (Unit=21, File=trim(modnam)//'/TEMP', Status='UNKNOWN', &
    Form='FORMATTED')

! -----------------------------------------------------------------------
! from V10 on: file 'CHIBAR_H': flux-weighted opacities

  Open (Unit=60, File=trim(modnam)//'/CHIBAR_H', Status='UNKNOWN', &
    Form='FORMATTED')

! -----------------------------------------------------------------------
! files used in sobolev+continuum approach (lines)
! file for j_nue values at different depth points

  Open (Unit=50, Status='SCRATCH', Recl=nde, Access='DIRECT')

! -----------------------------------------------------------------------
! file for h_nue values at different depth points, no longer req.

  Open (Unit=51, Status='SCRATCH', Recl=nde, Access='DIRECT')

! -----------------------------------------------------------------------
! file for k_nue values at different depth points

  Open (Unit=52, Status='SCRATCH', Recl=nde, Access='DIRECT')

! -----------------------------------------------------------------------
! file for n_nue values at different depth points, no longer req.

  Open (Unit=53, Status='SCRATCH', Recl=nde, Access='DIRECT')

! same as before but at a given depth points, for different
! nue values. 40: j; 41: h; 42: k; 43: n.

  Open (Unit=40, Status='SCRATCH', Recl=infre, Access='DIRECT')
  Open (Unit=42, Status='SCRATCH', Recl=infre, Access='DIRECT')

  Open (Unit=41, Status='SCRATCH', Recl=infre, Access='DIRECT')
  Open (Unit=43, Status='SCRATCH', Recl=infre, Access='DIRECT')

  Return

End Subroutine

!-----------------------------------------------------------------------

Subroutine ulecture

  Use :: nlte_type
  Use :: nlte_dim
  Use :: fastwind_params, Only: fpath_ufunc

  Use :: nlte_var, Only: ufung, ufunk
  Implicit None

! ---- it reads and stores u-function values from files utablak.dat and
! ---- utablag.dat. values are stores into module

! .. local scalars ..
  Real (dp) :: beta
  Integer (i4b) :: ib, it
! ..
! .. local arrays ..
  Real (dp) :: betag(17), betak(76)
! ..
! .. intrinsic functions ..
! INTRINSIC LOG10
! ..
! .. data statements ..
  Data betak/1.0D-10, 2.8D-10, 4.6D-10, 6.4D-10, 8.2D-10, 1.0D-9, 2.8E-9, &
    4.6D-9, 6.4D-9, 8.2D-9, 1.0D-8, 2.8D-8, 4.6D-8, 6.4D-8, 8.2D-8, 1.0D-7, &
    2.8D-7, 4.6D-7, 6.4D-7, 8.2D-7, 1.0D-6, 2.8D-6, 4.6D-6, 6.4D-6, 8.2D-6, &
    1.0D-5, 2.8D-5, 4.6D-5, 6.4D-5, 8.2D-5, 1.0D-4, 2.8D-4, 4.6D-4, 6.4D-4, &
    8.2D-4, 1.0D-3, 2.8D-3, 4.6D-3, 6.4D-3, 8.2D-3, 1.0D-2, 2.8D-2, 4.6D-2, &
    6.4D-2, 8.2D-2, 1.0D-1, 2.8D-1, 4.6D-1, 6.4D-1, 8.2D-1, 1.0D+0, 2.8D+0, &
    4.6D+0, 6.4D+0, 8.2D+0, 1.0D+1, 2.8D+1, 4.6D+1, 6.4D+1, 8.2D+1, 1.0D+2, &
    2.8D+2, 4.6D+2, 6.4D+2, 8.2D+2, 1.0D+3, 2.8D+3, 4.6D+3, 6.4D+3, 8.2D+3, &
    1.0D+4, 2.8D+4, 4.6D+4, 6.4D+4, 8.2D+4, 1.0D+5/

  Data betag/1.D-10, 1.D-9, 1.D-8, 1.D-7, 1.D-6, 1.D-5, 1.D-4, 1.D-3, 1.D-2, &
    1.D-1, 5.D-1, 1.D+0, 1.D+1, 1.D+2, 1.D+3, 1.D+4, 1.D+5/
! ..

  Open (Unit=70, File=fpath_ufunc//'utablak.dat', Status='OLD', &
    Form='FORMATTED')

  Open (Unit=71, File=fpath_ufunc//'utablag.dat', Status='OLD', &
    Form='FORMATTED')

  Do it = 1, 9
    Do ib = 1, 17
      Read (71, Fmt=*) ufung(it, ib)
    End Do
  End Do

  Do it = 1, 21
    Do ib = 1, 76
      Read (70, Fmt=*) ufunk(it, ib)
    End Do
  End Do

! ---- u-g function is changed into ug/beta if beta<1,

  Do it = 1, 9
    Do ib = 1, 11
      beta = betag(ib)
      ufung(it, ib) = ufung(it, ib)/beta
    End Do
  End Do

! ---- u-k function is changed into log10(uk) (log10(uk/beta) if beta<1)

  Do it = 1, 21
    Do ib = 1, 50
      beta = betak(ib)
      ufunk(it, ib) = log10(ufunk(it,ib)/beta)
    End Do

    Do ib = 51, 76
      ufunk(it, ib) = log10(ufunk(it,ib))
    End Do

  End Do

  Close (70)
  Close (71)
  Return

End Subroutine

!-----------------------------------------------------------------------

Subroutine u2lecture
  Use :: nlte_type
  Use :: fastwind_params, Only: fpath_ufunc

  Use :: nlte_var, Only: u2b, u2s, u2l1, u2l2
  Implicit None


! ..
! .. local scalars ..
  Integer (i4b) :: i, j, k
! ..
  Open (Unit=30, File=fpath_ufunc//'U2B.dat', Status='old')
  Open (Unit=31, File=fpath_ufunc//'U2S.dat', Status='old')
  Open (Unit=32, File=fpath_ufunc//'U2L1.dat', Status='old')
  Open (Unit=33, File=fpath_ufunc//'U2L2.dat', Status='old')

  Do j = 1, 21
    Read (30, Fmt=*)(u2b(j,k), k=1, 3)
    Read (30, Fmt=*)(u2b(j,k), k=4, 6)
  End Do

  Do i = 1, 192
    k = 0
    Do j = 1, 66
      Read (31, Fmt=*) u2s(i, j+k), u2s(i, j+k+1), u2s(i, j+k+2)
      k = k + 2
    End Do
    Read (31, Fmt=*) u2s(i, 199), u2s(i, 200)
  End Do

  Do i = 1, 61
    k = 0
    Do j = 1, 23
      Read (32, Fmt=*) u2l1(i, j+k), u2l1(i, j+k+1), u2l1(i, j+k+2)
      k = k + 2
    End Do
    Read (32, Fmt=*) u2l1(i, 70), u2l1(i, 71)
  End Do

  Do i = 1, 12
    k = 0
    Do j = 1, 6
      Read (33, Fmt=*) u2l2(i, j+k), u2l2(i, j+k+1), u2l2(i, j+k+2)
      k = k + 2
    End Do
  End Do

  Close (30)
  Close (31)
  Close (32)
  Close (33)

  Return

End Subroutine

!***********************************************************************

!subroutines: complex ones
!cont and related

!***********************************************************************

Subroutine cont(nd, np, r, rho, xne, temp, clf, teff, sr, corrfc, optlte, &
  concon, taur, rtau23, rybicki, dtfcorr, fltot, ndiv, iit)

! updated for clumping; all opacities correspond to average values!
! only change: redefinition of THOMSON

  Use :: nlte_type
  Use :: nlte_dim
  Use :: fund_const
  Use :: fastwind_params, Only: const_ev, xred1, xred2, xred3, xblu

  Use :: nlte_var, Only: fre, wfre, hc, ifre, lthin, r1k, opac, st => strue, &
    xj, alo, z, p, modnam, precis, unasol, optlines, thomson_lines, xxh, qq1, &
    qq2, optcmf_all
! ,RAYLEIGH

  Use :: photstruc, Only: chibar_h

  Use :: nlte_xrays, Only: optxray

  Use :: nlte_porvor, Only: opa_eff_rat ! , tcl_fac_line, tcl_fac_cont

  Implicit None

! .. parameters ..

! Integer (i4b), Parameter :: ifretot = id_frec1
  Integer (i4b), Parameter :: nd1 = id_ndept, np1 = id_npoin
! ..
! .. scalar arguments ..
  Real (dp) :: corrfc, rtau23, sr, teff
  Integer (i4b) :: nd, np, ndiv, iit
  Logical :: concon, optlte, rybicki
! ..
! .. array arguments ..
  Real (dp), Intent (In) :: r(nd1), rho(nd1), taur(nd1), temp(nd1), xne(nd1), &
    clf(nd1)
  Real (dp), Intent (Out) :: dtfcorr(nd1), fltot(nd1)

! ..
! .. local scalars ..
  Real (dp) :: aic, bnuecon, dbdr, dbdtf, dtdr, ferri, ferrmax, ferro, &
    fluxcin, fluxcon, fluxdif, op, opand, rthin, sig, tradic, tradic1, wfr, &
    xic, xicin, xicr, xicv1, xlambda, xhmin, desh, convlx, tau, dt, r1, hck

! not used  Real (dp) :: tcl

  Integer (i4b) :: imod, k, l, ixr1, ixr2, ixr3, ixb

! ..
! .. local arrays ..
  Real (dp) :: opa(nd1), q1(nd1-2), q2(nd1-2), thomson(nd1), wp(nd1+1, np1-1), &
    wp2(nd1+1, np1-1), xh(nd1+1), htot(nd1), lumx1(nd1+1), lumx2(nd1+1), &
    lumx3(nd1+1), aux(nd1+1)
! ..
! .. external functions ..
  Real (dp) :: dbdt !               , bnue
  External :: dbdt !                , bnue
! ..
! .. external subroutines ..
  External :: cont1, conversion, diffus
! ..
! .. intrinsic functions ..
! INTRINSIC ABS,LOG,LOG10,MAX,MOD
! ..

  sig = amh*sigmae
  fluxcin = 0.D0
  fluxcon = 0.D0
  fluxdif = 0.D0

  fltot = 0.D0
  htot = 0.D0
  chibar_h = 0.D0


! -----same argument as in subr. lte
! -----flux has to be corrected with respect to the lower radius according to
! rstar^2*flux(teff) = rmin^2*flux(rmin)
! -> flux(rmin)=(rstar/rmin)^2 flux(teff)=RTAU23^2 flux(teff)

  bnuecon = sigsb/pi*teff**4*rtau23**2 ! corresponds to 4 H (astrophys. flux)
! * (rstar/sr)^2 =
! total flux AT LOWER BOUNDARY
  If (.Not. optlte) Then
    dtdr = (temp(nd)-temp(nd-1))/(r(nd-1)-r(nd))
    Do k = 1, ifre
      wfr = wfre(k)
      dbdtf = dbdt(1.D8/fre(k), temp(nd))
      opand = opac(nd, k)
      fluxdif = fluxdif + wfr*dbdtf/opand
    End Do

    fluxdif = fluxdif*4.D0*dtdr/(3.D0*sr)
    corrfc = bnuecon/fluxdif
    Print *
    Print *, ' NEW CORRECTION FACTOR = ', corrfc
    Print *
    fluxdif = fluxdif*corrfc
  End If

! boundaries for X-ray treatment
  If (.Not. optlte .And. optxray) Then
    lumx1 = 0.D0
    lumx2 = 0.D0
    lumx3 = 0.D0

    Do k = 1, ifre
      If (fre(k)>xred1) Exit
    End Do
    ixr1 = k - 1

    Do k = ixr1, ifre
      If (fre(k)>xred2) Exit
    End Do
    ixr2 = k - 1

    Do k = ixr2, ifre
      If (fre(k)>xred3) Exit
    End Do
    ixr3 = k - 1

    Do k = ixr3, ifre
      If (fre(k)>=xblu) Exit
    End Do
    ixb = k

    If (fre(ixb)<=xblu) Then
      Write (999, *) ' STOP: CONT: IXB NOT FOUND IN FRE!'
      Stop ' CONT: IXB NOT FOUND IN FRE!'
    End If
!   missing RTAU23^2 corrected (Oct. 16, 2014)
    convlx = 4.D0*pi/(sigsb*teff**4*rtau23**2) ! for conversion, see below
    If (abs(4.D0-convlx*bnuecon)>1.D-6) Then
      Write (999, *) ' STOP: ERROR IN CONVLX'
      Stop ' ERROR IN CONVLX'
    End If
  End If

  If (ifre<50) Then
    imod = 1
  Else
    imod = ifre/50
  End If

  Print *
  Print *, ' NOMINAL RADIUS = ', rtau23
  Print *


  Write (*, Fmt='(2A)') '   NO        LAMBDA    LOG F-NUE(app) TRAD(a', &
    'pp)  TRAD1(app)     RTHIN   LTHIN  R(TAU=1)'
  Print *

! output of fluxes etc. to file FLUXCONT moved to routine PRINT_FLUXES

  xj(1, ifre+1) = 1 !               INDICATES OBSFRAM-TREATMENT

freloop: Do k = 1, ifre

    xlambda = 1.D8/fre(k)
    Call diffus(xlambda, temp, r, nd, aic, dbdr)

    Do l = 1, nd
      op = opac(l, k)
!     ----  OPAC already effective, from OPACITC
!     thus just rescale, and we're fine with ratio
!     (since thomson_lines etc. are mean quantities)
!     ----  NOTE: This routine called directly after OPACITC,
!     so OK to scale with opa_eff_rat (rather than old values opa_eff_rat_old)
      If (optlines) Then
!       THOMSON(L) =
!       (SIG*XNE(L)/CLF(L)+THOMSON_LINES(L,K)+RAYLEIGH(L,K)/CLF(L))/OP
        thomson(l) = (sig*xne(l)/clf(l)+thomson_lines(l,k))/op* &
          opa_eff_rat(l, k)
      Else
!       THOMSON(L) = (SIG*XNE(L)+RAYLEIGH(L,K))/CFL(L)/OP
        thomson(l) = (sig*xne(l))/clf(l)/op*opa_eff_rat(l, k)
      End If
      If (thomson(l)>=1.D0) Then
        Print *, ' WARNING: ERROR IN THOMSON', l, thomson(l)
        Print *, k, xlambda, op, sig, xne(l)/clf(l), thomson_lines(l, k)
        Write (999, *) ' WARNING: ERROR IN THOMSON', l, thomson(l)
        Write (999, *) k, xlambda, op, sig, xne(l)/clf(l), thomson_lines(l, k)
        thomson(l) = 1. !           should happen only between interpolation
!       of freq. grids
      End If
      opa(l) = op*sr
    End Do

    If (mod(k,imod)==0 .Or. unasol) Then
!     calculate r(tau_k=1) for output frequencies
      tau = taur(2)/10.
      Do l = 2, nd
        dt = 0.5*(opa(l-1)+opa(l))*(r(l-1)-r(l))
        If (tau+dt>1.) Go To 100
        tau = tau + dt
      End Do
      Print *, ' TAU = ', tau, ' (< 1!) AT ', xlambda
      Write (999, *) ' STOP: TAU(K) < 1 AT PHOTOSPHERE'
      Stop ' TAU(K) < 1 AT PHOTOSPHERE'

!     assuming log tau linear in r
100   r1 = r(l-1) - log10(tau)*(r(l)-r(l-1))/log10(1.D0+dt/tau)
      If (r1>r(l-1) .Or. r1<r(l)) Then
        Write (999, *) ' STOP: SOMETHING WRONG WITH R(TAU=1)!'
        Stop ' SOMETHING WRONG WITH R(TAU=1)!'
      End If
      r1k(k) = r1
    End If

!   -----------------------------------------------------------------------

!   calculation of continuum radiation field

!   if(xlambda.lt.8.79) then
!   print*,xlambda,aic,dbdr
!   do l=1,nd
!   print*,l,opa(l),thomson(l),st(l,k)
!   enddo
!   print*
!   endif

!   opa=1.d4
!   thomson=0.1
!   where(r .gt. 1.2) thomson=0.99
!   DO L=1,ND
!   BN(L)=BNUE(XLAMBDA,TEMP(L))
!   ENDDO
!   DTDR = (TEMP(ND)-TEMP(ND-1))/ (R(ND-1)-R(ND))
!   DT=BN(ND)/TEMP(ND)*DTDR
!   print*,'call',aic,dbdr
!   CALL CONT1(ND,NP,R,P,Z,BN,OPA,THOMSON,AIC,DBDR,CORRFC, &
!   WP,WP2,ALO(:,K),XJ(:,K),XH,.TRUE.,XLAMBDA,Q1,Q2,K,CONCON,RYBICKI)

!   ---- note: at the end of CONT1, XJ to XN (1,ND) are saved at each
!   frequency,
!   -----and further processed in subr. CONVERSION and FACTINC, to be used
!   -----in TRATRAD for Sobo-lines
    Call cont1(nd, np, r, p, z, st(:,k), opa, thomson, aic, dbdr, corrfc, wp, &
      wp2, alo(:,k), xj(:,k), xh, .True., xlambda, q1, q2, k, concon, rybicki)

    qq1 = q1
    qq2 = q2

!   ---- note that xh is not the eddington flux but (h*r**2). to obtain
!   ---- emergent absolute fluxes it has to be multiplied by 4*pi/r**2

!   for tests
!   IF(XLAMBDA.GT.912..AND.XLAMBDA.LT.1230.) XH=XH*0.7

    xic = 4.D0*xh(1)
    xicin = 4.D0*xh(nd+1)
    hc(k) = xic*pi/r(1)/r(1)

!   DBDTF = DBDT(1.D8/FRE(K),TEMP(ND))
!   OPAND = DBDTF/OPAC(ND,K)*DTDR/(3.D0*SR)*CORRFC
!   print*,k,xlambda,opand/xh(nd)


!   for test purposes

!   alo(:,k)=0.
!   do 150 l=1,nd
!   if(k.ge.1118.and.k.le.1150) then
!   write(*,151) xlambda,l,opa(l),st(l,k),thomson(l),xj(l,k)
!   print*,xlambda,l,xj(l,k),st(l,k),alo(l,k),thomson(l),opa(l)
!   endif
!   aux=1./(1.-alo(l,k)*thomson(l))
!   alolk=alo(l,k)*aux
!   if(xj(l,k)/st(l,k).lt.alolk) stop ' alo > j/s'
!   150  continue
!   print*
!   151  format(f8.3,1x,i2,1x,4(e8.3,1x))
!   -----------------------------------------------------------------------

!   for testing LTE
!   do l=1,nd
!   xj(l,k)=bnue(xlambda,temp(l))
!   alo(l,k) =0.
!   enddo

    Do l = 2, nd
      xicr = 4.D0*xh(l)
      If (abs(1.D0-xicr/xic)>.1D0) Go To 110
    End Do

110 Continue

    rthin = r(l-1)
    lthin(k) = l - 1

    hck = hc(k)
    If (xic>0.D0) Then
      tradic = 1.4388D8/xlambda/(log(3.97D8/xlambda**3/xic+1.D0))
      xicv1 = xic/(rtau23**2)
      tradic1 = 1.4388D8/xlambda/(log(3.97D8/xlambda**3/xicv1+1.D0))
    Else
      tradic = -1.4388D8/xlambda/(log(3.97D8/xlambda**3/abs(xic)+1.D0))
      xicv1 = xic/(rtau23**2)
      tradic1 = -1.4388D8/xlambda/(log(3.97D8/xlambda**3/abs(xicv1)+1.D0))
    End If

!   ---- normal running output (reduced)

!   print*,xlambda,opa(nd),xh(nd)
    If (mod(k,imod)==0) Write (*, Fmt=120) k, xlambda, log10(abs(hck)), &
      tradic, tradic1, rthin, lthin(k), r1k(k)

!   ---- absolute fluxes at r=rmax are displayed

!   ---- converged output (full) -> subr. PRINT_FLUXES

    wfr = wfre(k)
    fluxcon = fluxcon + wfr*xic
    fluxcin = fluxcin + wfr*xicin

!   ---- calculate X-ray luminosity in three bands (staggered grid)
    If (.Not. optlte .And. optxray .And. k>=ixr1 .And. k<=ixb) Then
      aux = wfr*xh
      lumx1 = lumx1 + aux
      If (k>=ixr2) lumx2 = lumx2 + aux
      If (k>=ixr3) lumx3 = lumx3 + aux
    End If

!   here was a bug (see below):
!   fltot defined with respect to 'regular' grid and not to staggered grid
!   DO L = 2, ND
!   FLTOT(L) = FLTOT(L) + 4.D0*XH(L)*WFR
!   END DO
    If (optlte) fluxdif = fluxdif + wfr*4.D0*dbdr*corrfc/(3.D0*opa(nd))

!   ---- calculation of flux-weighted opacity

!   ---- at first, interpolate H onto radial grid
!   ---- note that xh here is h*r**2
    xhmin = 0.D0
    Do l = 2, nd
      xhmin = min(xhmin, xh(l))
    End Do

    desh = 2.D0*abs(xhmin)

    Do l = 2, nd
      xh(l) = log10(xh(l)+desh)
    End Do

    Do l = 1, nd - 2
      xh(l+1) = q1(l)*xh(l+1) + q2(l)*xh(l+2)
      xh(l+1) = 10.D0**xh(l+1) - desh
    End Do
    xh(nd) = xh(nd+1)
!   new: XXH now on mesh points
    xxh(:, k) = xh(1:nd)

!   ---- now, calculate chibar_h, corrected for clumping
!   NOTE: dp/dm=4pi/(c rho) H chibar; dp=4pi/(c rho) H chibar rho dr
!   = 4pi/c H chibar dr;
!   to obtain correct result, volume filling factor is needed
!   ALREADY INCLUDED IN OPAC
!   JO check whether optically thick case consistent as well
!   Thus, no additional correction required if average density rho is used
!   ---- new version: standard formulation with chi/rho


    Do l = 1, nd
      op = 4.D0*xh(l)*wfr
      htot(l) = htot(l) + op
      chibar_h(l) = chibar_h(l) + op*opac(l, k)/rho(l)
    End Do
!   corrected in v10.4 (Jan 2015): fltot has to be calculated HERE
!   (after XH defined on regular mesh)
    Do l = 2, nd
      fltot(l) = fltot(l) + 4.D0*xh(l)*wfr
    End Do

  End Do freloop

! ---- transformation between frequency stored binary files and depth
! ---- ones

  If (concon) Call conversion

  If (abs(1.D0-fluxdif/bnuecon)>precis) Then
    Write (999, *) ' STOP: ERROR IN DBDR'
    Stop ' ERROR IN DBDR'
  End If

  Print *
  Print *, ' USED CORRECTION FACTOR = ', corrfc
  Print *
  Print *, '  SIGSB/PI * TEFF**4 * RTAU23^2 = ', bnuecon
  Print *, ' INTEGRAL(4*HNUE*DNUE)           = ', fluxcin, ' (R_IN)'
  Print *, ' INTEGRAL(4*HNUE*DNUE)* RMAX^2   = ', fluxcon, ' (R_OUT)'

  chibar_h = chibar_h/htot
! ---- write actual CHIBAR_H (reg. grid) and grad to file CHIBAR_H
  Rewind 60
  Do l = 1, nd
    Write (60, *) l, chibar_h(l), chibar_h(l)*sigsb/clight*teff**4
  End Do

! note: following quantities flux-errors;
! since H values actually (r/sr)^2*H and H_nom = sigsb/(4pi)*Teff^4 at r=rstar
! fltot/bnuecon=4*(r/sr)^2 * H_int(r)/(sigsb/pi*teff^4*(rstar/sr)^2 =
! (r/rstar)^2 * H_int(r) /(sigsb/(4pi)*teff^4) = 1 if no error!
  ferri = fluxcin/bnuecon - 1.D0
  ferro = fluxcon/bnuecon - 1.D0

  fltot(1) = ferro
  Do l = 2, nd
    fltot(l) = fltot(l)/bnuecon - 1.D0
  End Do
! new check
  If (fltot(nd)/=ferri) Then
    Write (999, *) ' STOP: ERROR IN FLTOT'
    Stop ' ERROR IN FLTOT'
  End If

  ferrmax = maxval(abs(fltot(2:nd)))

  Print *
  Print *, ' DIFFUSION APPROX. NOT VALID BY ', ferri
  Print *, ' FLUX (RMAX)  NOT  CONSERVED BY ', ferro
  Print *, ' MAXIMUM FLUX ERROR = ', ferrmax
  Print *

  If (.Not. optlte .And. optxray) Then
    lumx1 = lumx1*convlx
    lumx2 = lumx2*convlx
    lumx3 = lumx3*convlx
    Write (14, 130) iit, r(1), lumx1(1), lumx2(1), lumx3(1)

    Write (999, *)
    Write (999, *) '**************************************************'
    Write (999, *) 'radius (staggered)  L_x/L_bol'
    Write (999, *) '**************************************************'
    Write (999, *) '(range =', 1.D8/fre(ixr1), ' to', 1.D8/fre(ixb), ' A'
    Write (999, *) '       =', fre(ixr1)/const_ev, ' to', fre(ixb)/const_ev, &
      ' eV'

    Print *
    Write (*, *) '**************************************************'
    Write (*, *) 'radius (staggered)  L_x/L_bol'
    Write (*, *) '**************************************************'
    Write (*, *) '(range =', 1.D8/fre(ixr1), ' to', 1.D8/fre(ixb), ' A'
    Write (*, *) '       =', fre(ixr1)/const_ev, ' to', fre(ixb)/const_ev, &
      ' eV'
    Write (*, *) ' 2nd and 3rd range until', fre(ixr2)/const_ev, ' and', &
      fre(ixr3)/const_ev, ' eV)'
    l = 1
    Write (999, 130) l, r(l), lumx1(l), lumx2(l), lumx3(l)
    Write (999, *)
    Write (*, 130) l, r(l), lumx1(l), lumx2(l), lumx3(l)
    Do l = 2, nd
      desh = 0.5D0*(r(l-1)+r(l))
      Write (*, 130) l, desh, lumx1(l), lumx2(l), lumx3(l)
    End Do
    l = nd + 1
    Write (*, 130) l, r(l-1), lumx1(l), lumx2(l), lumx3(l)
    Print *
  End If

  If (unasol) Then
!   output of flux error now in routine PRINT_FLUXERROR

    If (.Not. optlte .And. optxray) Then
!     CONVERT TO Lx/Lbol
!     here was a (moderate) bug (corrected by JO, Oct. 16, 2014)
!     LX = 4 pi Rstar^2 * 4 pi (r/Rstar)^2 H = 4 pi Rstar^2 * 4 Pi (r/SR)^2 H
!     * (SR/Rstar)^2
!     => Lx/Lbol = 4 pi Int(Hc)/(sigsb * Teff^4 * RTAU23^2)
!     WITH RTAU23=RSTAR/SR
      Open (Unit=1, File=trim(modnam)//'/XLUM', Status='UNKNOWN')
      Write (1, *) '**************************************************'
      Write (1, *) 'radius (staggered)  L_x/L_bol'
      Write (1, *) '**************************************************'
      l = 1
      Write (1, 130) l, r(l), lumx1(l), lumx2(l), lumx3(l)
      Do l = 2, nd
        desh = 0.5D0*(r(l-1)+r(l))
        Write (1, 130) l, desh, lumx1(l), lumx2(l), lumx3(l)
      End Do
      l = nd + 1
      Write (1, 130) l, r(l-1), lumx1(l), lumx2(l), lumx3(l)
      Close (1)
    End If
  End If

! TEMPERATURE CORRECTION WITH RESPECT TO FLUX CONSERVATION
! JO SEPT 2025: ONLY FOR .NOT. OPTCMF_ALL
  If (.Not. optcmf_all) Call tcorr(taur, fltot, teff, temp, r, nd, dtfcorr, &
    ndiv)

! for tests
! xh=0.
! do k=1,ifre
! xh(1:nd)=xh(1:nd)+4.d0*xxh(1:nd,k)*wfre(k)
! enddo
! ferrmax=maxval(abs(1.d0-xh(1:nd)/bnuecon))
! print*,' max error (xxh) = ',ferrmax

! do l=1,nd
! print*,l,xh(l),(fltot(l)+1.d0)*bnuecon
! enddo
! stop

  Return

120 Format (I5, 2X, F15.5, 2X, F10.5, 2(2X,F10.0), 3X, F8.3, 3X, I3, 3X, F8.3)
130 Format (I4, '     ', F10.4, ' ', 3(E10.4,2X))

End Subroutine

!-----------------------------------------------------------------------

Subroutine print_fluxes(r)

! ---- fluxes etc. written to FLUXCONT (only in last iteration)
! different treatment in dependence of OPTCMF_ALL


  Use :: nlte_type
  Use :: nlte_dim
  Use :: fund_const
  Use :: nlte_var, Only: modnam, optcmf_all, fre, ifre, rtau23, xxh, &
    xh_obs_coarse, lthin, r1k

  Implicit None

  Integer (i4b), Parameter :: nd1 = id_ndept

  Real (dp), Dimension (nd1), Intent (In) :: r

  Integer (i4b) :: imod, k

  Real (dp) :: r1, xlambda, xic, xicv1, xicv2, hck, tradic, tradic1, rthin


  If (ifre<50) Then
    imod = 1
  Else
    imod = ifre/50
  End If

  r1 = r(1)

  If (optcmf_all) Then
    Write (*, Fmt='(2A)') '   NO        LAMBDA    LOG F-NUE(CMF) TRAD1(', &
      'app) TRAD1(CMF)     RTHIN   LTHIN  R(TAU=1)'
  Else
    Write (*, Fmt='(2A)') '   NO        LAMBDA    LOG F-NUE(app) TRAD(a', &
      'pp)  TRAD1(app)     RTHIN   LTHIN  R(TAU=1)'
  End If
  Print *

! ---- output file for continuum fluxes

  Open (23, File=trim(modnam)//'/FLUXCONT', Status='UNKNOWN', &
    Form='FORMATTED')

  If (optcmf_all) Then
    Write (23, Fmt='(2A)') '   NO        LAMBDA    LOG F-NUE(CMF) TRAD1(', &
      'app) TRAD1(CMF)     RTHIN   LTHIN  R(TAU=1)'
  Else
    Write (23, Fmt='(2A)') '   NO        LAMBDA    LOG F-NUE(app) TRAD(a', &
      'pp)  TRAD1(app)     RTHIN   LTHIN  R(TAU=1)'
  End If


freloop: Do k = 1, ifre

    xlambda = 1.D8/fre(k)
    xic = 4.D0*xxh(1, k)

!   different output, in dependence of OPTCMF_ALL

    If (optcmf_all) Then
!     note the different scaling of XH_OBS_COARSE and XH/XXH
      xicv1 = xic/(rtau23**2) !     for approx. solution
      xicv2 = 4.*xh_obs_coarse(k) ! corresponds to XICV1 for approx. method,
!     since different scaling
      hck = pi*xicv2/(r1/rtau23)**2 ! actual FLUX
      If (xicv1>0.D0) Then
        tradic = 1.4388D8/xlambda/(log(3.97D8/xlambda**3/xicv1+1.D0))
      Else
        tradic = -1.4388D8/xlambda/(log(3.97D8/xlambda**3/abs(xicv1)+1.D0))
      End If

      If (xicv2>0.D0) Then
        tradic1 = 1.4388D8/xlambda/(log(3.97D8/xlambda**3/xicv2+1.D0))
      Else
        tradic1 = -1.4388D8/xlambda/(log(3.97D8/xlambda**3/abs(xicv2)+1.D0))
      End If

    Else
      hck = xic*pi/r1/r1
      If (xic>0.D0) Then
        tradic = 1.4388D8/xlambda/(log(3.97D8/xlambda**3/xic+1.D0))
        xicv1 = xic/(rtau23**2)
        tradic1 = 1.4388D8/xlambda/(log(3.97D8/xlambda**3/xicv1+1.D0))
      Else
        tradic = -1.4388D8/xlambda/(log(3.97D8/xlambda**3/abs(xic)+1.D0))
        xicv1 = xic/(rtau23**2)
        tradic1 = -1.4388D8/xlambda/(log(3.97D8/xlambda**3/abs(xicv1)+1.D0))
      End If

    End If

    rthin = r(lthin(k))


!   ---- normal running output (reduced)

!   print*,xlambda,opa(nd),xh(nd)
    If (mod(k,imod)==0) Write (*, Fmt=100) k, xlambda, log10(abs(hck)), &
      tradic, tradic1, rthin, lthin(k), r1k(k)

!   ---- absolute fluxes at r=rmax are displayed


!   ---- converged output (full)

    Write (23, Fmt=100) k, xlambda, log10(abs(hck)), tradic, tradic1, rthin, &
      lthin(k), r1k(k)

  End Do freloop

  Write (23, Fmt=*) rtau23
  Close (23)

  Return

100 Format (I5, 2X, F15.5, 2X, F10.5, 2(2X,F10.0), 3X, F8.3, 3X, I3, 3X, F8.3)

End Subroutine

!-----------------------------------------------------------------------

Subroutine print_fluxerror(r, taur, fluxerr, fluxerr_cmf)

  Use :: nlte_type
  Use :: nlte_dim
  Use :: nlte_var, Only: optcmf_all, modnam

  Implicit None

  Integer (i4b), Parameter :: nd1 = id_ndept

! ..
! .. array arguments ..
  Real (dp), Dimension (nd1), Intent (In) :: r, taur, fluxerr, fluxerr_cmf

  Integer (i4b) :: l

  Open (23, File=trim(modnam)//'/FLUXCONT', Position='APPEND', &
    Form='FORMATTED')

  Write (23, Fmt=*)
  Write (23, Fmt=*) '  L        R        TAUR       FLUX-ERROR  '

  If (.Not. optcmf_all) Then

    l = 1
    Write (23, Fmt=100) l, r(l), taur(l), fluxerr(l)

    Do l = 2, nd1
      Write (23, Fmt=100) l, r(l), log10(taur(l)), fluxerr(l)
    End Do

  Else

    l = 1
    Write (23, Fmt=110) l, r(l), taur(l), fluxerr_cmf(l), fluxerr(l)

    Do l = 2, nd1
      Write (23, Fmt=110) l, r(l), log10(taur(l)), fluxerr_cmf(l), fluxerr(l)
    End Do

  End If

  Close (23)

  Return

100 Format (I3, 3(2X,F10.5))
110 Format (I3, 4(2X,F10.5))

End Subroutine

!-----------------------------------------------------------------------

Subroutine tcorr(taur, ferr, teff, temp, r, nd, dt, ndiv)
! temperature correction with respect to flux conservation

! to be used in lower part of atmosphere

! modified Lucy Unsoeld procedure (Mihalas page 174), &
! assuming all averaged opacities are equal.

! modified in V10.1 (check for closely separated grid-points, &
! and improve calculation of df/dt)

! Sept 2025: lower boundary changed, to stabilize iteration

  Use :: nlte_type
  Use :: nlte_var, Only: rtau23, optcmf_all

  Implicit None


  Integer (i4b), Intent (In) :: nd, ndiv

  Real (dp), Dimension (nd), Intent (In) :: taur, ferr, temp, r
  Real (dp), Intent (In) :: teff
  Real (dp), Dimension (nd), Intent (Out) :: dt

  Integer (i4b) :: i, is, ii
  Real (dp), Dimension (nd) :: df, dfdt, taup, dt1
  Real (dp) :: taulim, aux, dlogm, dlogp, taurm, taurp, ferrm, ferrp

  Logical :: test = .False.


  taulim = taur(ndiv)
  taulim = taulim*2. !              safety factor, to be well away from the
! transition zone
  If (taulim>0.1) taulim = 0.1 !    the old standard value;
! otherwise (larger), outer boundary for integrals insufficient
  If (taulim<0.001) taulim = 0.001

  Do i = 1, nd
    If (taur(i)>=taulim) Exit
  End Do

  is = i

  Print *, 'TCORR'
  df(is-1) = 0.
  Do i = is, nd

!   integral f(tau) dtau, expressed with respect to dlntau
    df(i) = df(i-1) + .5*(ferr(i-1)*taur(i-1)+ferr(i)*taur(i))*(log(taur(i))- &
      log(taur(i-1)))
    If (i==nd) Then
      dfdt(i) = (ferr(i-1)-ferr(i))/(taur(i-1)-taur(i))
    Else
!     derivative df(tau)/dtau, expressed with respect to dlntau
!     check for closely separated grid-points
      dlogm = log(taur(i-1)) - log(taur(i))
      dlogp = log(taur(i+1)) - log(taur(i))
!     print*,i,dlogm,dlogp

      If (abs(dlogm)<0.1 .Or. abs(dlogp)<0.1) Then ! factor 1.1
        taurm = taur(i)/1.1
        taurp = taur(i)*1.1
        Do ii = i - 1, 1, -1
          If (taur(ii)<taurm .And. taur(ii+1)>taurm) Exit
        End Do
        ferrm = ferr(ii) + (ferr(ii+1)-ferr(ii))/(taur(ii+1)-taur(ii))*(taurm- &
          taur(ii))
!       JO changed Jan 2016 (old version: loop till ND, no reset of II)
        Do ii = i, nd - 1
          If (taur(ii)<taurp .And. taur(ii+1)>taurp) Go To 100
        End Do
        ii = nd - 1 !               in case no point is found
100     ferrp = ferr(ii) + (ferr(ii+1)-ferr(ii))/(taur(ii+1)-taur(ii))*(taurp- &
          taur(ii))

        dfdt(i) = .5/taur(i)*((ferrm-ferr(i))/(log(taurm)-log(taur(i)))+(ferrp &
          -ferr(i))/(log(taurp)-log(taur(i))))
      Else
        dfdt(i) = .5/taur(i)*((ferr(i-1)-ferr(i))/dlogm+(ferr(i+ &
          1)-ferr(i))/dlogp)
      End If

      aux = abs(dfdt(i))
!     to avoid too large corrections in case of very closely separated
!     grid points; idea max df/dlntau from 0.03/(.33*ln10) (3 points per
!     decade) plus safety factor 1.5
      If (aux*taur(i)>0.06) Then
!       PRINT*,'orig',I,DFDT(I)
        dfdt(i) = aux/dfdt(i)*0.06/taur(i)
!       PRINT*,'corr',I,DFDT(I)
      End If
    End If
  End Do

  dt(1:is-1) = 0.D0
  Do i = is, nd
    dt(i) = 1.D0/16.D0*(teff/temp(i))**4*temp(i)*(-3.*df(i)+dfdt(i))
!   The following formulation is more exact than the one from above,
!   which results from a plane-parallel approach (see notes in Diplom-folder),
!   but gives identical results (tested!), since in both cases DT -> 0 for DF
!   -> 0
!   The lower formulation acounts for two issues
!   a) actually, H0 has been calculated at ND, i.e., H0 propto Teff^4*RTAU23
!   b) in spherical symmetry there is a factor 1/r^2 in front of DFDT,
!   and a factor 1/r'^2 inside the integral, which has been taken out as 1/r^2
!   DT(I)=1.D0/16.D0*(TEFF/TEMP(I))**4*TEMP(I)*(-3.*DF(I)+DFDT(I))*(RTAU23/R(I))**2
    aux = abs(dt(i))
!   restrict correction to 7.5%
    If (aux>0.075*temp(i)) dt(i) = 0.075*aux/dt(i)*temp(i) ! correct sign
  End Do

! JO March 2025
! following insert (commented out) testing a temperature correction
! calculated from
! T(taup) = Teff * (3/4(taup+2/3))^0.25
! works, but not faster than original code, and problems at transition
! around taup=3

! JO Sept 2025

  If (optcmf_all) Then
!   IF(test) THEN
    dt1 = dt

!   accurate (tested) approximation for taup (spherical generalisation of
!   taur)
    taup(1) = taur(1)
    Do i = 2, nd
      taup(i) = taup(i-1) + (taur(i)-taur(i-1))*rtau23**2/(r(i-1)*r(i))
    End Do

!   approximate value for DT at large optical depths
!   WHERE(TAUP.GE.3.)
!   DT=TEFF*(0.75*(TAUP+2./3.))**(0.25)-TEMP
!   ENDWHERE
    If (maxval(abs(dt(nd-1:nd)))>500.) Then
      Do i = nd - 1, nd
        dt(i) = teff*(0.75*(taup(i)+2./3.))**(0.25) - temp(i)
      End Do
    End If

  End If

  Do i = is, nd
!   for tests
!   PRINT*,DF(I),DFDT(I),-3.*DF(I)+DFDT(I)
    If (optcmf_all) Then
!     IF(test) THEN
      Write (*, Fmt='(1X,I2,5(1X,G14.7))') i, taur(i), ferr(i), temp(i), &
        dt(i), dt1(i)
    Else
      Write (*, Fmt='(1X,I2,4(1X,G14.7))') i, taur(i), ferr(i), temp(i), dt(i)
    End If
  End Do

  Return
End Subroutine

!-----------------------------------------------------------------------

Subroutine cont1(nd, np, r, p, z, st, opa, thomson, xic1, xic2, corr, wp, wp2, &
  alo, xj, xh, optflux, xlam, q1, q2, k, concon, rybicki)

! .......................................................................

! attention!!! subroutine cont1 changed for fast and consistent
! continuum transfer in program FASTWIND. so far, only case
! i) can be treated. in this new version, both j and h are taken
! from the solution of the moments' equation( here, h corresponds
! to r**2 * h-nue).
! additionally, the alo is calculated consistent with j.



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

! optflux = .true. : flux = xh is calculated, where
! xh(1) , xh(nd+1) boundary values consistent
! with actual i-plus,i-minus
! xh(i),i=2...nd   intermesh values at .5*(r(i-1)+r(i))


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

  Use :: nlte_type
  Use :: nlte_dim
  Use :: run_once, Only: start_weigh11
  Use :: nlte_var, Only: precis, vp, vp2, xxk

  Implicit None


  Integer (i4b), Parameter :: nd1 = id_ndept, np1 = id_npoin
  Integer (i4b), Parameter :: imax = 10
! ..
! .. scalar arguments ..
  Real (dp) :: corr, xic1, xic2, xlam
  Integer (i4b) :: k, nd, np
  Logical :: concon, optflux, rybicki
! ..
! .. array arguments ..
  Real (dp) :: alo(nd1), opa(nd1), p(np1), q1(nd1-2), q2(nd1-2), r(nd1), &
    st(nd1), thomson(nd1), wp(nd1+1, np1-1), wp2(nd1+1, np1-1), xh(nd1+1), &
    xj(nd1), z(nd1, np1)
! ..
! .. local scalars ..
  Real (dp) :: dtaug, dxi, e0, e1, edmue, epsmax, heddi, heddo, hin, hout, &
    pi2, ux, v, vi, xiplus
  Integer (i4b) :: iii, ij, jp, l, ll, lmax, lstart, lz
  Logical :: conv, core
! ..
! .. local arrays ..
  Real (dp) :: akp(nd1), akpjp(nd1*np1-1), aksum(nd1), f(nd1+3), qq(nd1), &
    ta(nd1), tajp(nd1*np1-1), taum(np1), tb(nd1), tb1(nd1), tbjp(nd1*np1-1), &
    tc(nd1), tc1(nd1), tcjp(nd1*np1-1), tp(nd1, nd1), tsum(nd1, nd1), &
    u(nd1, np1-1), up(nd1), upjp(nd1*np1-1), xhs(nd1+1), xjold(nd1), xk(nd1), &
    xns(nd1+1)
! ..
! .. external subroutines ..
  External :: dev, interpol, inv, invtri, invtri3, madd, mdmv, mdv, momalo, &
    momcont, mvmd, mvv, setup1, vadd, vsub, weigh11, weigh5, weigh9
! ..
! .. intrinsic functions ..
! INTRINSIC ACOS,ATAN,EXP,MIN,MIN0
! ..

  pi2 = acos(0.D0)
  conv = .False.

! calculation of integration weights for xj, xh only once
! ---- also calculation of interpolation weights for xh and xn

  If (start_weigh11) Then
    Call weigh11(nd, np, r, z, p, wp, wp2, q1, q2)
    Call weigh9
    Call weigh5
    start_weigh11 = .False.
  End If


  If (rybicki) Then
    aksum = .0D0
    tsum = .0D0
  End If

  Do jp = 1, np - 1

    If (jp==1) Then
      taum(1) = opa(1)*r(1)*(1.D0/3.D0+2.D0*thomson(1)/3.D0)
    Else
      taum(jp) = opa(1)*r(1)*r(1)*(thomson(1)/p(jp)*(pi2-atan(z(1,jp)/p(jp)))+ &
        (1.D0-thomson(1))*r(1)*r(1)/2.D0/p(jp)**3*(pi2-atan(z(1,jp)/ &
        p(jp))-z(1,jp)*p(jp)/r(1)/r(1)))
    End If

    If (taum(jp)<0.D0) Then
      Write (999, *) ' STOP: TAUM NEGATIVE!'
      Stop ' TAUM NEGATIVE!'
    End If

    lmax = min0(np+1-jp, nd)
    core = (lmax==nd)

!   calculation of tp,akp

    Call setup1(lmax, core, z(:,jp), st, opa, thomson, xic1, xic2, corr, up, &
      akp, ta, tb, tc, taum(jp))
    ta(1) = 0.D0
    tc(lmax) = 0.D0

!   ---  storage for further use

!   ---  start index

    lstart = nd1*(jp-1)
    Do l = 1, lmax
      ll = lstart + l
      upjp(ll) = up(l)
      akpjp(ll) = akp(l)
      tajp(ll) = ta(l)
      tbjp(ll) = tb(l)
      tcjp(ll) = tc(l)
    End Do

    If (.Not. rybicki) Go To 100

!   rybicki algorithm to obtain ajc

    Call invtri3(ta, tb, tc, tp, lmax, nd)
    Call mdmv(tp, vp(:,jp), lmax, nd)
    Call mvv(qq, tp, akp, lmax, lmax, nd)
    Call vadd(aksum, qq, lmax)

!   tp=:(vp*tp-1)*u

    Call mvmd(tp, up, lmax, nd)
    Call madd(tsum, tp, lmax, nd)
100 End Do

! addition of unity-matrix to tsum


  If (rybicki) Then
    Do l = 1, nd
      tsum(l, l) = tsum(l, l) + 1.D0
    End Do
    Call inv(nd, nd, tsum)
    Call mvv(xj, tsum, aksum, nd, nd, nd)
  End If

  qq = xj
  xjold = xj

  If (.Not. optflux .And. rybicki) Return

! backsubstitution to obtain u

  iii = 0

110 Continue

  iii = iii + 1

jploop: Do jp = 1, np - 1

    lmax = min0(np+1-jp, nd)
    core = (lmax==nd)

!   ---  start index

    lstart = nd1*(jp-1)
    Do l = 1, lmax
      ll = lstart + l
      up(l) = upjp(ll)
      akp(l) = akpjp(ll)
      ta(l) = -tajp(ll)
      tb(l) = tbjp(ll)
      tc(l) = -tcjp(ll)
    End Do

!   recalculation of total source-function

    Call mdv(qq, up, lmax)
    Call vsub(akp, up, lmax)

    Call invtri(ta, tb, tc, akp, lmax)

    Do l = 1, lmax
      u(l, jp) = akp(l)
    End Do

  End Do jploop

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

!   DXI = U(1,JP) - E1* (ST(1)+THOMSON(1)*XJ(1))

!   ---- this line has been masked (21-10-92). it corresponds to case
!   ---- i_plus<i_minus (at r_max), and now there is no difference
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

  Call momcont(nd, r, opa, thomson, st, xic1, xic2, corr, qq, xh, ta, tb, tc, &
    tb1, tc1, akp, alo, xj, f, .False.)

! -----------------------------------------------------------------------


  If (conv .Or. iii==imax) Then

    Do l = 1, nd
      xk(l) = f(l)*xj(l)
      xxk(l, k) = xk(l)
    End Do

    If (.Not. concon) Go To 120


!   --------------------------------------------------------------------------
!   ---- here we start the specific treatment for lines. h_nue and n_nue
!   are calculated so that ic can be expressed as a combination of
!   blocking factors and mue's. from gudrun taresch's program.
!   a new iteration step has been made, so that u's values be more
!   accurate

    xns = 0.D0
    xhs = 0.D0

    Do jp = 1, np - 1

      lmax = min(np+1-jp, nd)
      lz = lmax - 1

!     ---------- outer boundary for h,n

      e1 = 1.D0 - exp(-taum(jp))
      dxi = u(1, jp) - e1*(st(1)+thomson(1)*xj(1))
      xns(1) = xns(1) + dxi*wp2(1, jp)
      xhs(1) = xhs(1) + dxi*wp(1, jp)

!     ---------- intermesh points of h,n

      Do ll = 1, lz
        dtaug = 2.D0/((opa(ll)+opa(ll+1))*(z(ll+1,jp)-z(ll,jp)))
        v = (u(ll,jp)-u(ll+1,jp))*dtaug
        xhs(ll+1) = xhs(ll+1) + v*wp(ll+1, jp)
        xns(ll+1) = xns(ll+1) + v*wp2(ll+1, jp)
      End Do

!     ---------- inner boundary for h, n

      If (lmax==nd) Then
        xiplus = xic1 + z(nd, jp)*xic2/opa(nd)*corr
        vi = xiplus - u(lmax, jp)
        xhs(nd+1) = xhs(nd+1) + vi*wp(nd+1, jp)
        xns(nd+1) = xns(nd+1) + vi*wp2(nd+1, jp)
      End If

    End Do

!   ------- interpolation of h and n to get points on the grid

    Call interpol(xhs, xns, nd, q1, q2)

!   ------- copy of momenta to a binary file


    Write (50, Rec=k) xj
    Write (52, Rec=k) xk
    Write (51, Rec=k)(xhs(ij), ij=1, nd)
    Write (53, Rec=k)(xns(ij), ij=1, nd)
!   ----------------------------------------------------------------------------

120 Continue

    If (iii==imax .And. epsmax>1.D-5) Print *, &
      ' ACHIEVED ACCURACY IN CONT. TRANSPORT AT ', xlam, ' = ', epsmax

    Call momalo(nd, r, opa, thomson, st, ta, tb1, tc1, qq, f, akp, alo, xh, &
      xic2, corr)

    Return

  End If

! -----------------------------------------------------------------------

  Call dev(r, xj, xjold, nd, epsmax, .False.)

  If (epsmax<precis .Or. iii==imax) conv = .True.

  qq = xj
  xjold = xj

  Go To 110

End Subroutine

!----------------------------------------------------------------

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
! INTRINSIC ABS,MAX
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

Subroutine interpol(xhs, xns, nd, q1, q2)

  Use :: nlte_type
  Use :: nlte_dim
  Implicit None

! ---- this subroutine interpolates h and n to get points on the r-grid
! ---- changed (26/3/93) so that negative fluxes does not affect the
! ---- logarithmic interpolation
! ----
! ---- comment from jo: this is an elegant, but not very exact way,
! ---- since it is more appropriate for linear interpolation. however,
! ---- since log (c*r^n + desh) = log(c'*r^n) with c' = c + desh*r^(-n),
! ---- the real assumption underlying this approach is that c' instead
! ---- of c is constant over the considered interval. this is not too
! ---- bad and was checked by myself that it works fairly well.

! ---- if somebody has doubt about the above statement or there are
! ---- conspicous cases, uncomment the commented lines below and try!

! .. parameters ..
  Integer (i4b), Parameter :: nd1 = id_ndept
! ..
! .. scalar arguments ..
  Integer (i4b) :: nd
! ..
! .. array arguments ..
  Real (dp) :: q1(nd1-2), q2(nd1-2), xhs(nd1+1), xns(nd1+1)
! ..
! .. local scalars ..
  Real (dp) :: desh, desn, xhmin, xnmin
  Integer (i4b) :: l
! ..
! .. intrinsic functions ..
! INTRINSIC ABS,LOG10,MIN
! ..

  xhmin = 0.D0
  xnmin = 0.D0
  Do l = 2, nd
    xhmin = min(xhmin, xhs(l))
    xnmin = min(xnmin, xns(l))
  End Do

  desh = 2.D0*abs(xhmin)
  desn = 2.D0*abs(xnmin)

  Do l = 2, nd
    xhs(l) = log10(xhs(l)+desh)
    xns(l) = log10(xns(l)+desn)
  End Do

  Do l = 1, nd - 2
    xhs(l+1) = q1(l)*xhs(l+1) + q2(l)*xhs(l+2)
    xns(l+1) = q1(l)*xns(l+1) + q2(l)*xns(l+2)
    xhs(l+1) = 10.D0**xhs(l+1) - desh
    xns(l+1) = 10.D0**xns(l+1) - desn
!   xhmin=min(xhold(l+1),xhs(l+1),xhold(l+2))
!   xhmax=max(xhold(l+1),xhs(l+1),xhold(l+2))
!   xnmin=min(xnold(l+1),xns(l+1),xnold(l+2))
!   xnmax=max(xnold(l+1),xns(l+1),xnold(l+2))
!   if(xhmin.eq.xhs(l+1)) stop ' error in interpolation'
!   if(xhmax.eq.xhs(l+1)) stop ' error in interpolation'
!   if(xnmin.eq.xns(l+1)) stop ' error in interpolation'
!   if(xnmax.eq.xns(l+1)) stop ' error in interpolation'
  End Do

  xhs(nd) = xhs(nd+1)
  xns(nd) = xns(nd+1)

  Return

End Subroutine

!-----------------------------------------------------------------------

Subroutine momalo(nd, r, opa, thomson, st, ta, tb1, tc1, q, f, akp, alo, xh, &
  xic2, corr)

! AKP on input is here r^2 *J

  Use :: nlte_type
  Use :: nlte_dim
  Implicit None


! .. parameters ..
  Integer (i4b), Parameter :: nd1 = id_ndept
! ..
! .. scalar arguments ..
  Integer (i4b) :: nd
  Real (dp) :: xic2, corr
! ..
! .. array arguments ..
  Real (dp) :: akp(nd1), alo(nd1), f(nd1+3), opa(nd1), q(nd1), r(nd1), &
    st(nd1), ta(nd1), tb1(nd1), tc1(nd1), thomson(nd1), xh(nd1+1)
! ..
! .. local scalars ..
  Real (dp) :: dtau1, r2
  Integer (i4b) :: l
! ..
! .. external subroutines ..
  External :: diag
! ..

  Call diag(ta, tb1, tc1, alo, nd, 1)

  r2 = r(1)*r(1)
  xh(1) = akp(1)*(f(nd+1)-f(nd+3)*thomson(1)) - f(nd+3)*st(1)*r2

  Do l = 2, nd
    dtau1 = 2.D0/((r(l-1)-r(l))*(opa(l-1)+opa(l)))
    xh(l) = (f(l)*q(l)*akp(l)-f(l-1)*q(l-1)*akp(l-1))*dtau1/(.5D0*(q(l- &
      1)+q(l)))
!   LEAVE THIS STATEMENT HERE; OTHERWISE AN OPTIMIZATION ERROR WITH THE
!   RESULT XH = 0 MIGHT OCCUR!
!   ALSO, THIS CAN HAPPEN INDEED (FOR X-RAYS, IN FIRST ITERATION)
!   and also in iteration zero
    If (xh(l)==0.D0) Then
!     PRINT*,L,DTAU1,XH(L)
!     PRINT*,F(L),Q(L),AKP(L)
!     PRINT*,F(L-1),Q(L-1),AKP(L-1)
!     PRINT*,'WARNING!!!! XH = 0 IN MOMALO!!!!!'
      xh(l) = 1.D-40
!     STOP ' XH = 0 IN MOMALO!'
    End If
  End Do

  xh(nd+1) = tc1(nd) - akp(nd)*f(nd+2)
! PRINT*,'test',XH(ND+1),TC1(ND),AKP(ND),F(ND+2)
! JO Jan 2016: might happen after update of hyd-structure
  If (xh(nd+1)<=0.D0) Then
    xh(nd+1) = xic2*corr/(3.D0*opa(nd))
!   PRINT*,'CORR',XH(ND+1)
  End If

  Return

End Subroutine

!-----------------------------------------------------------------------

Subroutine diag(a, b, c, dia, n, iopt)

! diagonal of inverse of tridiag matrix
! (acc. to Rybicky & Hummer, 1991, A&A 245, 171)

  Use :: nlte_type
  Use :: nlte_dim
  Implicit None

! .. scalar arguments ..
  Integer (i4b) :: iopt, n
! ..
! .. array arguments ..
  Real (dp) :: a(n), b(n), c(n), dia(n)
! ..
! .. local scalars ..
  Integer (i4b) :: l
! ..
! .. local arrays ..
  Real (dp) :: d(0:100), e(101)
! ..

  If (n>100) Then
    Write (999, *) ' STOP: TOO MANY GRID POINTS IN DIAG'
    Stop ' TOO MANY GRID POINTS IN DIAG'
  End If

  d(0) = 0.D0

  Do l = 1, n
    d(l) = c(l)/(b(l)-a(l)*d(l-1)) ! trick: although a(1) undefined, d(0)=0.
  End Do
  e(n+1) = 0.D0

  Do l = n, 1, -1
    e(l) = a(l)/(b(l)-c(l)*e(l+1)) ! trick: although c(n) undefined, e(n+1)=0.
  End Do

  If (iopt==0) Then
    Do l = 1, n
      dia(l) = 1.D0/((1.D0-d(l)*e(l+1))*(b(l)-a(l)*d(l-1)))
    End Do

  Else If (iopt==1) Then
    Do l = 1, n
      dia(l) = dia(l)/((1.D0-d(l)*e(l+1))*(b(l)-a(l)*d(l-1)))
    End Do

  Else
    Write (999, *) ' STOP: WRONG OPTION IN DIAG'
    Stop ' WRONG OPTION IN DIAG'

  End If

  Return

End Subroutine

!-----------------------------------------------------------------------

Subroutine momcont(nd, r, opa, thomson, st, xic1, xic2, corr, q, xh, ta, tb, &
  tc, tb1, tc1, akp, alo, xj, f, optflux)

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
  Logical :: optflux
! ..
! .. array arguments ..
  Real (dp) :: akp(nd1), alo(nd1), f(nd1+3), opa(nd1), q(nd1), r(nd1), &
    st(nd1), ta(nd1), tb(nd1), tb1(nd1), tc(nd1), tc1(nd1), thomson(nd1), &
    xh(nd1+1), xj(nd1)
! ..
! .. local scalars ..
  Real (dp) :: dt0, dtm, dtp, edmue, fl, flp, hi, ho, rl, rl2, rlp, rrq
  Integer (i4b) :: l
! ..
! .. external subroutines ..
  External :: invtri, momalo
! ..
! .. intrinsic functions ..
! INTRINSIC EXP
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
  alo(1) = edmue

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
    alo(l) = 1.D0/q(l)
    akp(l) = r(l)*r(l)*st(l)/q(l)
  End Do

! inner boundary condition

  l = nd
  ta(l) = f(l-1)*q(l-1)*dtp
  tb(l) = f(l)*q(l)*dtp + hi
  tb1(l) = tb(l)
  alo(l) = 0.D0
  akp(l) = .5D0*xic1 + xic2*corr/(3.D0*opa(l))
! print*,'hplus',.5*xic1,xic2*corr/(3.d0*opa(l)),akp(l)
  tc1(nd) = akp(l)

  Call invtri(ta, tb, tc, akp, nd)
! now, AKP is r^2 J

  Do l = 1, nd
    rl2 = r(l)*r(l)
    xj(l) = akp(l)/rl2
  End Do

  If (.Not. optflux) Return

! calculation of alo and h

  Call momalo(nd, r, opa, thomson, st, ta, tb1, tc1, q, f, akp, alo, xh, xic2, &
    corr)

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
! INTRINSIC EXP
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

  If (lz<2) Go To 100

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

100 Continue
  l = lmax

! inner boundary, 2nd order

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

Subroutine weigh11(nd, np, r, z, p, w1, w3, q1, q2)

  Use :: nlte_type
  Use :: nlte_dim
  Use :: nlte_var, Only: w0 => vp, w2 => vp2, w0last, vp1, vp3


  Implicit None

! ***  calculation of angular integration weights
! ***  the weights w1  are calculated for intermesh points
! ***  and interpolation weights q1,q2 for h,n

! ***  changed from V10 on: everything calculated, &
! ***  called after update of structure


! .. parameters ..
  Integer (i4b), Parameter :: nd1 = id_ndept, np1 = id_npoin
! ..
! .. scalar arguments ..
  Integer (i4b) :: nd, np
! ..
! .. array arguments ..
  Real (dp) :: p(np1), q1(nd1-2), q2(nd1-2), r(nd1), w1(nd1+1, np1-1), &
    w3(nd1+1, np1-1), z(nd1, np1)
! ..
! .. local scalars ..
  Real (dp) :: a, aa, b, bb, c, cc, rd1, rd2, rl, rl2, rrr12, rrrr20, rz, &
    rzq6, w1lz, w3lz, ww1, ww3, xlog2

! for tests
  Real (dp) :: sum1, sum2, sum3, sum4, asum1, asum2, asum3, asum4
  Integer (i4b) :: jp, l, lmax, lz
! ..
! .. intrinsic functions ..
! INTRINSIC LOG10,MIN0
! ..

jploop1: Do jp = 1, np - 1

    lmax = min0(np+1-jp, nd)
    lz = lmax - 1
!   ***
!   ***  0. and 2. moment. the integration is performed in the z variable
!   ***
lloop1: Do l = 1, lmax
      rl = r(l)
      rl2 = rl + rl
      rrr12 = rl*rl*rl*12.D0
!     ***
!     ***  first step if jp=1
!     ***
      If (jp==1) Then
        b = z(l, 1)
        a = z(l, 2)
        w0(l, jp) = (b-a)/rl2
        aa = a*a
        bb = b*b
        w2(l, jp) = (b*(3.D0*bb-aa)-a*(bb+aa))/rrr12
      Else
        If (l/=lmax .Or. jp<=(np-nd)) Then
!         ***
!         ***  intermediate step
!         ***
          a = z(l, jp+1)
          b = z(l, jp)
          c = z(l, jp-1)
          w0(l, jp) = (c-a)/rl2
          aa = a*a
          bb = b*b
          cc = c*c
          w2(l, jp) = (b*(cc-aa)+c*(cc+bb)-a*(bb+aa))/rrr12
        Else
!         ***
!         ***  last step, implying z(l,jmax)=0
!         ***
          b = z(l, jp-1)
          w0(l, jp) = b/rl2
          w2(l, jp) = b*b*b/rrr12
        End If
      End If

    End Do lloop1

  End Do jploop1
! ***
! *** integration weight for p = rmax and z=0
! ***
  w0last = z(1, np-1)/(r(1)+r(1))

jploop2: Do jp = 1, np - 1

    lmax = min0(np+1-jp, nd)
    lz = lmax - 1
!   ***
!   *** 1.moment.
!   *** first step
!   ***
    If (jp==1) Then
      ww1 = p(2)*p(2)
      ww3 = ww1*ww1
    Else
!     ***
!     ***  intermediate steps
!     ***
      a = p(jp-1)
      b = p(jp)
      c = p(jp+1)
      ww1 = (a+b+c)*(c-a)
      aa = a*a
      bb = b*b
      cc = c*c
      ww3 = (b-a)*(aa*(a+2.D0*b)+bb*(3.D0*a+4.D0*b)) + &
        (c-b)*(cc*(c+2.D0*b)+bb*(3.D0*c+4.D0*b))
!     ***
!     ***  for the last interval (l=lz to the p axis), the next point is an
!     ***  intermesh-point in p
!     ***
      c = .5D0*(b+c)
      w1lz = (a+b+c)*(c-a)
      cc = c*c
      w3lz = (b-a)*(aa*(a+2.D0*b)+bb*(3.D0*a+4.D0*b)) + &
        (c-b)*(cc*(c+2.D0*b)+bb*(3.D0*c+4.D0*b))
    End If
!   ***
!   ***  no weight for the z=0 point is calculated, as v=0 there for
!   ***  symmetry.
!   ***
!   ***  loop over depth index l
!   ***

lloop2: Do l = 1, lz
      rz = .5D0*(r(l)+r(l+1))
      rzq6 = rz*rz*6.D0
      rrrr20 = rz*rz*rz*rz*20.D0
      If (l/=lz .Or. jp<=(np-nd)) Then
        w1(l+1, jp) = ww1/rzq6
        w3(l+1, jp) = w1(l+1, jp) - ww3/rrrr20

      Else
        w1(l+1, jp) = w1lz/rzq6
        w3(l+1, jp) = w1(l+1, jp) - w3lz/rrrr20
      End If
    End Do lloop2

!   ***
!   ***  special weights at the outer boundary for h
!   ***
    rl = r(1)
    w1(1, jp) = ww1/rl/rl/6.D0
    w3(1, jp) = w1(1, jp) - ww3/rl/rl/rl/rl/20.D0
!   ***
!   ***  special weights at the inner boundary
!   ***
    If (lmax<nd) Cycle
    If (jp>(np-nd)) Then
!     ***  core tangent ray, last step of integration
      ww1 = (b-a)*(2.D0*b+a)
      ww3 = (b-a)*(aa*(a+2.D0*b)+bb*(3.D0*a+4.D0*b))
    End If
    w1(nd+1, jp) = ww1/6.D0
    w3(nd+1, jp) = w1(nd+1, jp) - ww3/20.D0

  End Do jploop2

! ***
! ***  interpolation weights
! ***

  xlog2 = log10(2.D0)
! weights rechecked Feb. 2015, Ok.
! Q2=(log r1'-log r)/(log r1' - log r2') = log(r1'/r)/log(r1'/r2') with
! r  = r(l+1)
! r1'=1/2(r(l)+r(l+1))
! r2'=1/2(r(l+1)+r(l+2)
  Do l = 1, nd - 2
    rd1 = (r(l)+r(l+1))/r(l+1)
    rd2 = (r(l)+r(l+1))/(r(l+1)+r(l+2))
    q2(l) = (log10(rd1)-xlog2)/log10(rd2)
    q1(l) = 1.D0 - q2(l)
  End Do

  vp1 = w1
  vp3 = w3
! Weights for rad_force saved also in common block, JS

  Return

! for tests of H and N at R(L=1,1.5,.5*(ND-1+ND),and ND)
  sum1 = 0.
  sum2 = 0.
  sum3 = 0.
  sum4 = 0.
  asum1 = 0.
  asum2 = 0.
  asum3 = 0.
  asum4 = 0.
  Do jp = 1, np - 1

    lmax = min0(np+1-jp, nd)
    lz = lmax - 1
    sum1 = sum1 + w1(67, jp)
    sum2 = sum2 + w1(60, jp)
    sum3 = sum3 + w1(2, jp)
    sum4 = sum4 + w1(1, jp)
    asum1 = asum1 + w3(67, jp)
    asum2 = asum2 + w3(60, jp)
    asum3 = asum3 + w3(2, jp)
    asum4 = asum4 + w3(1, jp)
    Print *, jp, lmax, lz
    Print *, 0.5*(r(66)+r(67))
    Print *, '66.5', w1(67, jp), sum1, w3(67, jp), asum1
    Print *, 0.5*(r(59)+r(60))
    Print *, '59.5', w1(60, jp), sum2, w3(60, jp), asum2
    Print *, 0.5*(r(1)+r(2))
    Print *, '1.5', w1(2, jp), sum3, w1(3, jp), asum3
    Print *, r(1)
    Print *, '1.0', w1(1, jp), sum4, w3(1, jp), asum4


  End Do

End Subroutine

!-----------------------------------------------------------------------

Subroutine weigh5

  Use :: nlte_type
  Use :: nlte_dim
  Use :: nlte_var, Only: w5

  Implicit None

! ---- calculates integration weights for five points integration
! ---- (simpson modified)

  w5(1) = 14.D0/45.D0
  w5(2) = 64.D0/45.D0
  w5(3) = 24.D0/45.D0
  w5(4) = w5(2)
  w5(5) = w5(1)

  Return

End Subroutine

!-----------------------------------------------------------------------

Subroutine weigh9

  Use :: nlte_type
  Use :: nlte_dim
  Use :: nlte_var, Only: w9

  Implicit None

! ---- calculates integration weights for nine points integrations

  w9(1) = 14.D0/45.D0
  w9(2) = 64.D0/45.D0
  w9(3) = 24.D0/45.D0
  w9(4) = w9(2)
  w9(5) = 28.D0/45.D0
  w9(6) = w9(2)
  w9(7) = w9(3)
  w9(8) = w9(2)
  w9(9) = w9(1)

  Return

End Subroutine

!-----------------------------------------------------------------------

Subroutine conversion

  Use :: nlte_type
  Use :: nlte_dim
  Use :: nlte_var, Only: ifre

  Implicit None

! ---- it converts binary files stored by frecuency into binary files
! ---- stores by depth

! .. parameters ..
  Integer (i4b), Parameter :: ifretot = id_frec1
  Integer (i4b), Parameter :: nd = id_ndept
! ..
! .. local scalars ..
  Integer (i4b) :: ifr, ind, j
! ..
! .. local arrays ..
  Real (dp) :: a(ifretot, nd)
! ..
! ..
  Do ifr = 1, ifre
    Read (50, Rec=ifr)(a(ifr,j), j=1, nd)
  End Do

  Do ind = 1, nd
    Write (40, Rec=ind)(a(j,ind), j=1, ifre)
  End Do


  Do ifr = 1, ifre
    Read (52, Rec=ifr)(a(ifr,j), j=1, nd)
  End Do

  Do ind = 1, nd
    Write (42, Rec=ind)(a(j,ind), j=1, ifre)
  End Do

  Do ifr = 1, ifre
    Read (51, Rec=ifr)(a(ifr,j), j=1, nd)
  End Do

  Do ind = 1, nd
    Write (41, Rec=ind)(a(j,ind), j=1, ifre)
  End Do


  Do ifr = 1, ifre
    Read (53, Rec=ifr)(a(ifr,j), j=1, nd)
  End Do

  Do ind = 1, nd
    Write (43, Rec=ind)(a(j,ind), j=1, ifre)
  End Do

  Return
End Subroutine

!-----------------------------------------------------------------------

Function dbdt(x, t)

  Use :: nlte_type
  Use :: nlte_dim
  Use :: fund_const, Only: hkl
  Implicit None


! .. parameters ..
! REAL(DP), PARAMETER :: C1=1.4388354967334D8
  Real (dp), Parameter :: c1 = hkl*1.D8
! ..
! .. scalar arguments ..
  Real (dp) :: t, x, dbdt
! ..
! .. external functions ..
  Real (dp) :: bnue
  External :: bnue
! ..
! .. intrinsic functions ..
! INTRINSIC EXP
! ..

  dbdt = bnue(x, t)*c1/x/t/t/(1.D0-exp(-c1/x/t))

  Return

End Function

!-----------------------------------------------------------------------

Function bnue(xlambda, t)

  Use :: nlte_type
  Use :: nlte_dim
  Use :: fund_const, Only: hc2, hkl
  Implicit None


! planck function,lambda in angstroem,t in kelvin
! bnue in erg per (cm**2 * sec * hertz )

! constanten: c1=h*c/k,c2=2*h*c


! .. parameters ..
! REAL(DP), PARAMETER :: C1=1.4388354967334D8,C2=3.972970127D8
  Real (dp), Parameter :: c1 = hkl*1.D8, c2 = hc2*1.D24
! ..
! .. scalar arguments ..
  Real (dp) :: t, xlambda, bnue
! ..

  Real (dp) :: x, expo
! .. intrinsic functions ..
! INTRINSIC EXP
! ..
! MODIFIED TO ACCOUNT FOR LARGE EXPONENTS (X-RAY REGIME)

  x = c1/xlambda/t

  If (x<200.D0) Then
    bnue = c2/(exp(x)-1.D0)/xlambda/xlambda/xlambda
  Else
    expo = log(c2/xlambda**3) - x
    bnue = exp(expo)
  End If

  Return

End Function

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

!***********************************************************************

!subroutines: complex ones
!lte and related

!***********************************************************************

Subroutine lte(teff, sr, taur, taue, r, velo, dvdr, rho, xne, xnh, temp, &
  dilfac, clf, ilow, imax, nd, grey, corrfc, ntemp, corvm, optcv, ithyd, xkap, &
  ex, xmmax, rtau23, xmet, restart_met_lte, lte_update, update_taur)
! clf is clumping factor (overdensity in clumps)

  Use :: nlte_type
  Use :: nlte_dim
  Use :: fund_const
  Use :: nlte_var, Only: fre, wfre, ifre, opac, sigem, delsig, xkapith, exith, &
    corrith, modnam, iconver, optmet, yhein, lines_in_model, delpara, op_flag, &
    op_flag2, frecfn1
  Use :: princesa_var, Only: frecfn

  Use :: nlte_porvor, Only: fic, tcl_fac_cont


  Implicit None

! .. parameters ..
  Integer (i4b), Parameter :: kel = id_atoms
! Integer (i4b), Parameter :: ifretot = id_frec1
  Integer (i4b), Parameter :: nd1 = id_ndept
! ..
! .. scalar arguments ..
  Real (dp) :: corrfc, corvm, ex, rtau23, sr, teff, xkap, xmmax
  Integer (i4b) :: ithyd, nd, ntemp
  Logical :: grey, optcv, restart_met_lte, lte_update, update_taur
! ..
! .. array arguments ..
  Real (dp) :: dilfac(nd1), dvdr(nd1), r(nd1), rho(nd1), taue(nd1), taur(nd1), &
    temp(nd1), velo(nd1), xne(nd1), xnh(nd1), clf(nd1)

  Integer (i4b) :: ilow(nd1, kel), imax(nd1, kel)
! ..
! ..
! .. local scalars ..
  Real (dp) :: a, bt, corr, del, deltmax, epsmax, flux, fluxross, op, opth, &
    rs, sig, srconst, summ, tmin, xlam, xmet, x, ttaur, tcl
  Integer (i4b) :: i, iii, iit, k, l, ll, lross, n, nt1, ifreold
  Logical :: conv, opteta, optout1, opttemp, newfreq
! ..
! .. local arrays ..
! TAUPTEST(ND1),
  Real (dp) :: arho(nd1), brho(nd1), crho(nd1), oparold(nd1), opaross(nd1), &
    sumopar(nd1), temp1(nd1), xnesav(nd1), dummy(nd1), taur_nlte(nd1)
  Data ifreold/0/

! ..
! .. external functions ..
  Real (dp) :: dbdt
  External :: dbdt
! ..
! .. external subroutines ..
  External :: dev, frescal, integra, linreg, nistar, occlte, opacitc, prmodel, &
    spline, tlucy
! ..
! .. intrinsic functions ..
! INTRINSIC ABS,LOG10,MAX,SQRT
! ..

  opteta = .False.
! not used:  writeoc = .False.
  optout1 = .False.

  Do l = 1, nd
    rs = rtau23/r(l)
    If (rs>1.D0) Then
      dilfac(l) = .5D0
    Else
      dilfac(l) = .5D0*(1.D0-sqrt(1.D0-rs**2))
    End If
  End Do

  If (grey) Then
    conv = .False.
    iit = 20
  Else
    conv = .True.
    iit = 1
  End If

  Do l = 2, nd
!   IF (TEMP(L).GT.TEMP(L-1)) EXIT
!   CHANGED
    If (temp(l)/=temp(l-1)) Exit
  End Do

  ntemp = l - 1

  If (ntemp==1) Then
    opttemp = .True.
    tmin = 0.D0
  Else
    opttemp = .False.
    tmin = temp(ntemp)
  End If

  Print *, ' TMIN = ', tmin
  Print *

! calculation of grey temperature

  If (grey) Then
    Print *
    Print *, ' TEMPERATURE CORRECTION '
    Print *
  End If

  ntemp = 1

! ===================================================================

ittemp: Do iii = 1, iit

    If (conv) Then
      opteta = .True.
!     not used: writeoc = .False.
      optout1 = .True.
    End If

    If (.Not. opttemp .And. iii>=2) Then
      Do l = 1, nd
        If (temp(l)>tmin) Then
          ntemp = l - 1
          Go To 100
        End If
      End Do
    End If

100 Continue

    If (conv) ntemp = 1

    Call occlte(rho, xne, xnh, clf, temp, ilow, imax, nd, ntemp, lte_update)
!   all occup. numbers include clf
!   NOT changed in case of OPTTHICK, since occupation numbers remain

    Call nistar

    If (ifreold/=0) ifreold = ifre
    Call frescal(nd, xnh, temp(nd), optout1, teff, newfreq)

!   IF(IFRE.NE.IFREOLD) THEN
!   JO Sept. 2021: New condition
    If (newfreq) Then
!     IF FREQ. GRID HAS CHANGED, NTEMP = 1 REQUIRED TO UPDATE ALL DEPTH POINTS
      ntemp = 1
      tmin = temp(1)
      ifreold = ifre
!     A change in the frequency grid requires to read again the OP rbf
!     data file, to map the cross sections in the new frequency grid; miguel
!     Sept. 2016: in this case, reset FRECFN1
      op_flag(:) = .False.
      op_flag2(:) = .False.
      frecfn1 = frecfn
      Print *, ' OP_FLAG AND FRECFN1 RESET'
    End If

    If (optmet) Then
      Print *, ' LTE FOR METALS'

!     note: dvdr only dummy argument here

!     CALL NLTE_APPROX(TEFF,YHEIN,XMET,RHO,XNE,DILFAC,TEMP,DUMMY,DUMMY,DUMMY,
!     &
!     &        CLF,ND,.TRUE.,.TRUE.,NTEMP,OPTOUT1)
!     dummy arguments changed (Jan. 2015)
      Call nlte_approx(teff, yhein, xmet, rho, xne, dilfac, temp, dummy, &
        dummy, dummy, clf, ilow, imax, nd, .True., .True., ntemp, optout1)
!     NOTE: no call of update_occng required, since lteopt=.true.

      Call opacitm(nd, xne, temp, clf, .True., ntemp, opteta, .True., .False.)

    End If

    Call opacitc(nd, xne, xnh, temp, clf, ilow, imax, .True., ntemp, opteta, &
      .True., optmet)
!   remember, that all opacities have been corrected for clumping
!   to allow for correctly taken spatial averages

!   if restart and very first LTE run, 2nd call in order
!   to provide correct OPAC_NOLINES

    If (restart_met_lte .And. lines_in_model) Then
      Print *, ' RESTART AND LTE: 2ND CALL OF OPACITL TO OBTAIN &
        &CONSISTENT LTE-OPACITIES'
      Call opacitl(xne, .True., ntemp, dummy, dummy, dummy, clf, .True.)

      Call opacitm(nd, xne, temp, clf, .True., ntemp, opteta, .True., .False.)

      Call opacitc(nd, xne, xnh, temp, clf, ilow, imax, .True., ntemp, opteta, &
        .True., optmet)
    End If


    Do k = 1, ifre
      xlam = 1.D8/fre(k)
      Do l = ntemp, nd
        opaross(l) = dbdt(xlam, temp(l))/opac(l, k)
!       ---- OPAC is effective, thus effective rosseland mean as well!
      End Do
      Call integra(opaross, sumopar, nd, ifre, wfre(k), k)
    End Do

!   print*,'renormalization for opaross '

    fluxross = opaross(nd)*(temp(nd)-temp(nd-1))/(r(nd-1)-1.D0)/sr/3.D0

!   only used if ntemp>1, i.e. if iii.gt.1

    If (iii==1 .And. ntemp/=1) Then
      Write (999, *) ' STOP: ERROR IN CHOSEN PHILOSOPHY -- OPAROLD'
      Stop ' ERROR IN CHOSEN PHILOSOPHY -- OPAROLD'
    End If

    Do l = 1, ntemp - 1
      opaross(l) = oparold(l)
    End Do

    Do l = ntemp, nd
      bt = 7.21863D-5*temp(l)**3
      opaross(l) = bt/opaross(l)/r(l)/r(l)
    End Do

    oparold = opaross

!   note that opaross is corrected for clumping

!   ---- calculation of tau-lucy and tau-ross
!   ---- arho: tau-lucy
!   ---- brho: tau-ross
!   ---- crho:(effective or mean) opa-ross (optthick = .true./false.)

!   -----x with respect to srnom, but r with respect to sr
!   hence: constant in front of integral srnom*rtau23=sr*rtau23^2

    srconst = sr*rtau23**2

    Do l = 1, nd
      crho(l) = opaross(l)*srconst
    End Do

    summ = 0.D0
    arho(1) = 0.D0
    Do l = 1, nd - 1
      del = r(l+1) - r(l)
      summ = summ - del*5.D-1*(crho(l)+crho(l+1))
      arho(l+1) = summ
    End Do

    If (conv) Go To 110

    Do l = 1, nd
      crho(l) = opaross(l)*r(l)*r(l)*sr
    End Do

    summ = 0.D0
    brho(1) = 0.D0

    Do l = 1, nd - 1
      del = r(l+1) - r(l)
      summ = summ - del*5.D-1*(crho(l)+crho(l+1))
      brho(l+1) = summ
    End Do

    Call tlucy(r, arho, brho, temp1, nd, tmin, teff, rtau23, ithyd)
!   print*,'from hydro,tlucy',temp1
    Call dev(r, temp1, temp, nd, epsmax, .False.)
    Print *

    temp = temp1

    If (epsmax<0.01D0) Then
      conv = .True.

!     .....storage of electron density
!     (from next to last hydro-iteration, to ensure precise restart)

      xnesav = xne
    End If

  End Do ittemp

! ===================================================================

! -----calculation of rosseland-optical depth

110 Continue

  sig = sigmae*amh

  Do l = 1, nd
    sumopar(l) = opaross(l)*r(l)*r(l)*sr
  End Do

  summ = 0.D0
  taur(1) = 0.D0
  Do l = 1, nd - 1
    del = r(l+1) - r(l)
    summ = summ - del*5.D-1*(sumopar(l)+sumopar(l+1))
    taur(l+1) = summ
  End Do

! ---- if needed, vmin is corrected to obtain a large enough value of
! ---- tau-ross at the inner point (e. santolaya, 11/2/94)

  If (taur(nd)<50.D0) Then
    optcv = .True.
    corvm = taur(nd)/60.D0
  Else
    optcv = .False.
  End If

  optcv = optcv .And. grey

  If (optcv .And. ithyd==1) Return

  If (taur(nd)>80.D0) xmmax = xmmax*70.D0/taur(nd)

  If (taur(nd)<60.D0) xmmax = xmmax*70.D0/taur(nd)

  Do l = 1, nd - 1
    If (taur(l)<.66D0 .And. taur(l+1)>=.66D0) Go To 120
  End Do

  Print *, ' WARNING!!!! TAUROSS = 2/3 NOT FOUND'
  Print *
  Write (999, *) ' WARNING!!!! TAUROSS = 2/3 NOT FOUND'
  Write (999, *)

  lross = 0
  Go To 130

120 Continue

  lross = l + 1

130 Continue

  Print *
  Print *, ' FOLLOWING TAU-ROSS VALUES ARE LTE-VALUES'
  Print *, ' TAU-ROSS AT RSTAR = ', taur(nd)
  Print *

  If (lross/=0) Then
    Print *, ' TAU-ROSS = ', taur(lross), ' AT R = ', r(lross)/rtau23, &
      ' NOM. RADII'
    ttaur = temp(lross)**4 + (temp(lross)**4-temp(lross-1)**4)/(taur(lross)- &
      taur(lross-1))*(0.66-taur(lross))
    ttaur = ttaur**0.25
    Print *, ' T(TAUROSS=2/3) = ', ttaur
    Print *
  End If

  Print *
  Print *, ' FOLLOWING TABLE WITH EFFECTIVE OPACITIES (BUT STILL LTE)'
  Print *
  Print *, ' NO    OPAROSS     OPAR/OPATH    OPTH/RHO      TAUROSS', &
    '       TAUP         TEMP'
  Print *

! JO March 2025: testing an approximate calculations for taup. WORKS!!!
! (with high accuracy for values la 0.2)
! TAUPTEST(1)=ARHO(1)
! PRINT*,RTAU23,R(1),R(ND)
! DO L=2,ND
! TAUPTEST(L)=TAUPTEST(L-1)+(TAUR(L)-TAUR(L-1))*RTAU23**2/(R(L-1)*R(L))
! ENDDO

  Do l = 1, nd
    op = opaross(l)*r(l)*r(l)
    opth = sig*xne(l)/clf(l) !      since op has been corrected
    tcl = opth*tcl_fac_cont(l)
    opth = opth*(1.+fic(l)*tcl)/(1.+tcl)
!   ---- printing effective thomson opacity, as if the atmosphere
!   would have only this as opacity source.
!   ARHO(L),TAUPTEST(L),TEMP(L)/TEFF
    Write (*, Fmt=190) l, log10(op), op/opth, opth/rho(l), taur(l), arho(l), &
      temp(l)/teff
  End Do

! -----this is the Eddington flux with respect to the nom. radius rstar(teff)

  flux = 4.511693D-6*teff**4 !      sigma_B/(4pi)

! -----it has to be corrected with respect to the lower radius according to
! rstar^2*flux(teff) = rmin^2*flux(rmin)
! -> flux(rmin)=(rstar/rmin)^2 flux(teff)=RTAU23^2 flux(teff)

  flux = flux*rtau23**2

  corrfc = flux/fluxross

  Print *
  Print *, ' CORRFC FOR LOWER BOUNDARY  = ', corrfc
  Print *
  If (corrfc<.8D0 .Or. corrfc>1.2D0) Then
    Print *, ' WARNING!!! CORRFC EXTREMELY LARGE'
    Print *
    Write (999, *) ' WARNING!!! CORRFC EXTREMELY LARGE'
    Write (999, *)
  End If

! JO, April 2016: following block included; taur will be overwritten by NLTE
! value
  If (update_taur) Then
    Open (1, File=trim(modnam)//'/TAU_ROS', Status='UNKNOWN', &
      Form='FORMATTED')
    Rewind 1
    Do l = 1, nd
      Read (1, Fmt=*) taur_nlte(l), x ! xne includes clf
    End Do
    Close (1)

    taur = taur_nlte

    Print *
    Print *, ' TAU-ROSS UPDATED BY CMF/NLTE VALUE'

    Do l = 1, nd - 1
      If (taur(l)<.66D0 .And. taur(l+1)>=.66D0) Exit
    End Do

    lross = l + 1

    Print *
    Print *, ' TAU-ROSS(NLTE) AT RSTAR = ', taur(nd)
    Print *

    Print *, ' TAU-ROSS = ', taur(lross), ' AT R = ', r(lross)/rtau23, &
      ' NOM. RADII'
    ttaur = temp(lross)**4 + (temp(lross)**4-temp(lross-1)**4)/(taur(lross)- &
      taur(lross-1))*(0.66-taur(lross))
    ttaur = ttaur**0.25
    Print *, ' T(TAUROSS_NLTE=2/3) = ', ttaur
    Print *

  End If

! ne(lte) written to tau_ros for calculation of dep. coeffiecients in
! totout.f90

  Open (1, File=trim(modnam)//'/TAU_ROS', Status='UNKNOWN', Form='FORMATTED')
  Rewind 1
  Do l = 1, nd
    Write (1, Fmt=*) taur(l), xne(l) ! xne includes clf
  End Do
  Close (1)
! print*,'xnecheck: tau_ros written'
! print*,taur
! print*,xne

! .....only for models which start from given hydro-structure
! (.not.grey corresponds to optmod)

  If (.Not. grey) Go To 150

  Do l = 1, nd
    sumopar(l) = sig*xne(l)*sr/clf(l) ! corrected
!   ---- corrected here as above, assuming the atmosphere has only e-
!   scattering
    tcl = sumopar(l)*tcl_fac_cont(l)/sr ! SR already included
    sumopar(l) = sumopar(l)*(1.+fic(l)*tcl)/(1.+tcl)
!   ---- effective opacity
    delsig(l) = sumopar(l)/(sigem*rho(l)*sr)
  End Do

! regression coefficients for kramer opacities used for photospheric
! correction in wind model

  Do l = 1, nd
    If (taur(l)>.1D0) Go To 140
  End Do
  Write (999, *) ' STOP: NT1 NOT FOUND'
  Stop ' NT1 NOT FOUND'

140 Continue

! FROM HERE ON, WE LEAVE EVERYTHING AS IN THE UNCLUMPED CASE, &
! SINCE IT APPLIES ONLY TO THE PHOTOSPHERE

  nt1 = l
  If (ithyd>=15) nt1 = nt1 + 1

! all quantities now refer to average density

! correction for delpara (if .ne. 1)
  Do l = nt1, nd
    ll = l - nt1 + 1
    arho(ll) = log10((opaross(l)*r(l)**2*sr/sumopar(l)-1.D0)/(rho(l)*delpara( &
      l)))
    brho(ll) = log10(temp(l))
  End Do

  n = nd - nt1 + 1

  Call linreg(brho, arho, n, ex, a, corr)
  xkap = 10.D0**a
  ex = -ex
  xkapith(ithyd) = xkap
  exith(ithyd) = ex
  corrith(ithyd) = corr
  Print *
  Print *, '   LINEAR REGRESSION OF FLUX-MEAN OF TYPE :'
  Print *
  Print *, &
    '   -- CHI(R) = OPA-TH(R)*(1.+ XKAP*DELPARA(L)*RHO(L)/T(L)**EX) -- '
  Print *
  Print *, '   FROM R = 1. UP TO  R = ', r(nt1)
  Print *
  Print *, '   YIELDS XKAP = ', xkap, '   EX = ', ex
  Print *, '   CORR. COEFF = ', corr
  Print *


  If (abs(corr)<.8D0) Then
    Write (999, *) ' STOP: CORR IN KRAMER OPACITIES TOO LARGE'
    Stop ' CORR IN KRAMER OPACITIES TOO LARGE'
  End If

! ACCOUNTING FOR DIFFERENCES BETWEEN PARAMETERIZED AND ACTUAL OPACITY
  If (teff<12000. .And. ithyd>=5 .And. ithyd<15) Then
!   if H-recom is oscillating, freeze DELPARA
    Do l = 1, nd
      arho(l) = opaross(l)*r(l)**2*sr/sumopar(l) - 1.D0
      brho(l) = xkap*rho(l)*temp(l)**(-ex)
      delpara(l) = arho(l)/brho(l)
      If (delpara(l)<0.D0) delpara(l) = 0.D0
    End Do
  End If

  If (iconver/=1) Go To 160

  Print *, ' ITERATION HISTORY OF KRAMER OPACITY REGR. COEFF.'
  Print *

  Do l = 1, ithyd
    Write (*, 180) l, xkapith(l), exith(l), corrith(l)
    Print *
  End Do
  Print *

  Write (999, *) ' ITERATION HISTORY OF KRAMER OPACITY REGR. COEFF.'
  Write (999, *)

  Do l = 1, ithyd
    Write (999, *)
    Write (999, 180) l, xkapith(l), exith(l), corrith(l)
  End Do

  Write (999, *)

  Rewind 21
  Read (21, *)(brho(i), i=nd, 1, -1), (temp1(i), i=nd, 1, -1)
  deltmax = 0.D0
  Do l = 1, nd
    deltmax = max(deltmax, abs(1.D0-temp1(l)/temp(l)))
  End Do
! print*,'from hydro,readfrom temp',temp1

  Rewind 21
  Write (21, *)(brho(i), i=nd, 1, -1), (temp(i), i=nd, 1, -1)
  Write (21, *)(xnesav(i), i=nd, 1, -1)
  Write (21, *) rtau23
! print*,'from hydro,writeto temp',temp

  Do l = 1, nd
    arho(l) = sumopar(l)/sr*(1.D0+xkap*delpara(l)*rho(l)*temp(l)**(-ex))
  End Do

  Print *
  Print *, '   CHECK VALUES , DELTA-MAX(T_HYDRO,T_ACT) = ', deltmax
  Print *
  Print *, '   ALL OPACITIES CORRESPOND TO AVERAGE DENSITY'
  Print *
  Print *, &
    ' NO   OPA(CALC)    OPA(APPROX)   T(HYDRO)       T(ACT)      DELPARA'

  Do l = 1, nd
    Write (*, Fmt=200) l, log10(opaross(l)*r(l)**2), log10(arho(l)), temp1(l), &
      temp(l), delpara(l)
  End Do

  Write (999, *)
  Write (999, *) '   CHECK VALUES , DELTA-MAX(T_HYDRO,T_ACT) = ', deltmax
  Write (999, *)
  Write (999, *) '   ALL OPACITIES CORRESPOND TO AVERAGE DENSITY'
  Write (999, *)
  Write (999, *) &
    ' NO   OPA(CALC)    OPA(APPROX)   T(HYDRO)   T(ACT)   DELPARA'

  Do l = 1, nd
    Write (999, Fmt=200) l, log10(opaross(l)*r(l)**2), log10(arho(l)), &
      temp1(l), temp(l), delpara(l)
  End Do

150 Continue

! JO changed Jan. 2016
! recalc taue in micro_clumping approx.
! (above, calculated including optically thick clumping)
  sumopar = sig*xne*sr/clf !        corrected
  Call spline(r, sumopar, arho, brho, crho, nd)
  summ = 0.D0
  taue(1) = 0.D0

  Do l = 1, nd - 1
    del = r(l+1) - r(l)
    summ = summ - del*(sumopar(l)+del*(arho(l)/2.D0+del*(brho(l)/3.D0+ &
      del*crho(l)/4.D0)))
    taue(l+1) = summ
  End Do

  Print *
  Print *, ' -------------- ACTUAL WIND MODEL--------------------------'
  Print *

! REMEMBER: TAUE=TAU-TH ONLY IN MICRO-CLUMPING APPROX.
  Call prmodel(r, velo, dvdr, rho, xne, xnh, temp, taue, clf, teff, nd)

  Print *
  Print *, ' -------------------------------------------------------------'
  Print *

160 Continue

  Do l = 2, nd
    If (temp(l)>temp(l-1)) Go To 170
  End Do

170 Continue

  ntemp = l - 1

  Return
180 Format (' IT.NO. ', I2, ' XKAP = ', E12.6, ' EX = ', F10.5, ' CORR = ', &
    F10.5)

! 9000 FORMAT (I3,8 (3X,F10.5))  !incl. tauptest
190 Format (I3, 7(3X,F10.5))
200 Format (I3, 3X, 2(F10.5,3X), 2(F10.0,3X), G12.6)

End Subroutine

!-----------------------------------------------------------------------

Subroutine frescal(nd, xnh, tmax, optout, teff, newfreq)

! -----  this subroutine calculates the frequency grid, nh corrected for
! clumping

! present philosophy requires that EACH level of a DETAILed ion
! is represented with two corresponding frequencies, namely at
! the corresponding edge and at edge-epslon(edge)

! levmax : maximum number of resolved edges for dominant ions of
! the two most abundant atoms
! levmin1: maximum number of resolved edges for 2nd dominant ions
! of the two most abundant atoms
! levmin2: maximum number of resolved edges for all other ions

! NOTE: we always count from the highest transition frequency.
! Two resolved edges means using the 2 highest frequencies.
! This does not necessarily correspond to the two lowest levels,
! since transitions to excited levels might lie at higher freqs.
! Example: SiIII (usually 2nd ion). Highest transitions are
! Si301->401 and Si304->402

! resolved means here, that we use a sufficient number of freq.
! points above the corresponding edges.

  Use :: nlte_type
  Use :: nlte_dim
  Use :: fund_const

  Use :: princesa_var, Only: nat, nions, le, li, zeff, abund, labat, labl, &
    labl4, frecin, paren4

  Use :: nlte_var, Only: mainion, fre, wfre, ifre, enionnd, imia, imaa, &
    iconver, bn, index, lam_ly, nstark, wavestark, freold, ifreold1, optmet, &
    fre_lines, ifre_lines

  Use :: nlte_opt, Only: optcmf_full

  Use :: nlte_xrays, Only: optxray, n_kedges, k_nmin, k_nmax, z, nionk => n, &
    eth, name

  Use :: nlte_wavcon, Only: wwavblue, wwavred, wwavcon

  Implicit None

! .. parameters ..
  Integer (i4b), Parameter :: kel = id_atoms, kis = id_kisat
  Integer (i4b), Parameter :: nfmin = 300
  Integer (i4b), Parameter :: levmax = id_levma, levmin1 = id_levm1, &
    levmin2 = id_levm2
  Integer (i4b), Parameter :: ifretot = id_frec1
  Integer (i4b), Parameter :: nd1 = id_ndept

! JO Dec. 2015
! FREQUENCIES FOR Lya and Lyb which should be used for consistency
! precision 1.d-9
! JO July 2021
! extended for Lyg and Lyd to avoid problems for B-stars
  Real (dp), Parameter :: lya_lo = 1.D8/1217.00043262610D0, &
    lya_hi = 1.D8/1212.85710645286D0, lyb_lo = 1.D8/1027.46136261866D0, &
    lyb_hi = 1.D8/1024.50655646359D0, lyg_lo = 1.D8/973.196165978500D0, &
    lyg_hi = 1.D8/972.310779289066D0, lyd_lo = 1.D8/949.843120035809D0, &
    lyd_hi = 1.D8/948.999697008538D0, lam_ciii = 977.0200D0 ! CIII RESONANCE
! LINE (1-3)
! ..
! .. scalar arguments ..
  Real (dp) :: tmax, teff
  Integer (i4b) :: nd
  Logical :: optout, newfreq, start
! ..
! .. array arguments ..
  Real (dp) :: xnh(nd1)
! ..
! .. local scalars ..
  Real (dp) :: aaa, dm, dmin, dmin1, dmin2, dnew, dnew1, edge, emax, eps, &
    errp, frenew, pf, pfe, sum1, sum2, x, x1, xm1, xm2, xnh1, fremin, xmin, &
    xmax, dfl, kayser
  Integer (i4b) :: i, if1, ihi, ii, il, ilab, imp, inew, io, isi, iz, j, ji, &
    k, k1, k2, l, leve, leve1, leve2, levp, levporg, levtot, levused, ll, lli, &
    m, m1, m2, mm1, mm2, n, nimp, nopa, nper, nk, istart, iend, nskip, ifreold
  Character :: lab*6, pare*6
! ..
! .. local arrays ..
  Real (dp) :: flaux(id_rbftr), freedge(id_rbftr)
  Integer (i4b) :: indedge(id_rbftr), iopa(kis*kel)
! ..
! .. external functions ..
  Real (dp) :: bnue, xintfre
  External :: bnue, xintfre
! ..
! .. external subroutines ..
  External :: sort
! ..
! .. intrinsic functions ..
! INTRINSIC ABS,DBLE,INT,LOG,LOG10,MAX,MIN
! ..
! .. statement functions ..
  Real (dp) :: epslon
! ..
! .. statement function definitions ..

  epslon(edge) = 5.D0*10.D0**(dble(int(log10(edge)))-6.D0)
! ..
! optout=.true.
! iconver=1
  Data start/.True./

  If (id_rbftr>999) Then
    Write (999, *) ' STOP: TOO MANY RBF-TRANSITIONS, CODING OF INDEX AFFECTED'
    Stop ' TOO MANY RBF-TRANSITIONS, CODING OF INDEX AFFECTED'
  End If

! JO Sept. 2021
  ifreold = ifre
  If (ifreold/=0) Then
    ifreold1 = ifre
    Print *, 'IFREOLD1 1:', ifreold1
    Do i = 1, ifre
      freold(i) = fre(i)
    End Do
  End If

! most abundant element k1

  k1 = 1
  x = abund(1)
  Do k = 2, nat
    x1 = max(x, abund(k))
    If (x1==x) Cycle
    x = x1
    k1 = k
  End Do

! 2nd abundant element if present k2

  If (nat>1) Then
    k2 = 1
    If (k1==1) k2 = 2
    x = abund(k2)
    Do k = 1, nat
      If (k==k1) Cycle
      x1 = max(x, abund(k))
      If (x1==x) Cycle
      x = x1
      k2 = k
    End Do
  End If

! dominant ionisation stages of k1, k2 (if present)

  sum1 = 0.D0
  sum2 = 0.D0

  Do l = 1, nd
    sum1 = sum1 + mainion(l, k1)
    If (nat>1) sum2 = sum2 + mainion(l, k2)
  End Do

  xm1 = sum1/nd
  xm2 = sum2/nd
  m1 = int(xm1)
  m2 = int(xm2)
  If (xm1-m1>.5D0) m1 = m1 + 1
  If (xm2-m2>.5D0) m2 = m2 + 1

  iopa(1) = k1*100 + m1

  If (nat>1) Then
    iopa(2) = k2*100 + m2
  Else
    k2 = 0
    If (m2/=0) Then
      Write (999, *) ' STOP: ERROR IN M2'
      Stop ' ERROR IN M2'
    End If
    iopa(2) = 0
  End If

  If (imaa(k1)>1) Then

!   2nd important ionization stage of k1 (if present)

    If (imaa(k1)==2) Then

!     e.g. helium

      mm1 = 3 - m1
      iopa(3) = k1*100 + mm1
    Else

      If (m1==imaa(k1)) Then
        mm1 = m1 - 1
        Go To 100
      End If

!     above or below m1?

      sum1 = 0.D0
      sum2 = 0.D0
      Do l = 1, nd
        xnh1 = 1.D0/xnh(l)
        sum1 = sum1 + enionnd(k1, m1+1, l)*xnh1
        sum2 = sum2 + enionnd(k1, m1-1, l)*xnh1
      End Do
      sum1 = sum1/dble(nd)
      sum2 = sum2/dble(nd)
      If (sum1>sum2) Then
        mm1 = m1 + 1
      Else
        mm1 = m1 - 1
      End If

100   Continue
      iopa(3) = k1*100 + mm1
    End If
    nopa = 3
  Else

!   most abundant elememt has only one ionization stage

    mm1 = m1
    nopa = 2
  End If

  If (nat>1) Then

!   2nd important ionization stage of k2, as for k1

    If (imaa(k2)>1) Then
      If (imaa(k2)==2) Then
        mm2 = 3 - m2
        iopa(nopa+1) = k2*100 + mm2
      Else
        If (m2==imaa(k2)) Then
          mm2 = m2 - 1
          Go To 110
        End If
        sum1 = 0.D0
        sum2 = 0.D0
        Do l = 1, nd
          xnh1 = 1.D0/xnh(l)
          sum1 = sum1 + enionnd(k2, m2+1, l)*xnh1
          sum2 = sum2 + enionnd(k2, m2-1, l)*xnh1
        End Do
        sum1 = sum1/dble(nd)
        sum2 = sum2/dble(nd)
        If (sum1>sum2) Then
          mm2 = m2 + 1
        Else
          mm2 = m2 - 1
        End If

110     Continue
        iopa(nopa+1) = k2*100 + mm2
      End If
      nopa = nopa + 1
    Else
      mm2 = m2

!     may happen if only one ion per element present
!     if(nopa.eq.2) stop ' error in nopa'

    End If
  End If

! total number of important ions (with levmax,levmin1)

  nimp = nopa

  Do k = 1, nat
!   commented out, since we might consider also depleted elements
!   IF (ABUND(K)/ABUND(K1).LE.1.D-7) GO TO 90
    Do m = imia(k), imaa(k)
      If (k==k1 .And. (m==m1 .Or. m==mm1)) Cycle
      If (k==k2 .And. (m==m2 .Or. m==mm2)) Cycle
      nopa = nopa + 1
      iopa(nopa) = k*100 + m
    End Do
  End Do

  If (nopa>kis*kel) Then
    Write (999, *) ' STOP: ERROR IN NOPA'
    Stop ' ERROR IN NOPA'
  End If

  Print *
  Print *, ' USED IONS FOR FREQUENCY GRID'
  Print *

  Do n = 1, nopa
    If (iopa(n)==0) Cycle
    ihi = iopa(n)/100
    iz = int(zeff(ihi))
    isi = iopa(n) - 100*ihi
    If (nions(ihi)==isi) Then
      Write (999, *) ' STOP: ERROR IN NIONS OR ISI'
      Stop ' ERROR IN NIONS OR ISI'
    End If
    levused = levmax
    If (n>2) levused = levmin1
    If (n>nimp) Then !              NEW
      If (imaa(ihi)-imia(ihi)>=2) Then ! at least 3 ions present
        levused = levmin2
      Else
        levused = levmin1
      End If
    End If
    Write (*, Fmt=320) labat(ihi), isi, iz + isi - 1, levused
  End Do

  Print *

  ifre = 0
  levtot = 0
  Do n = 1, nopa

    If (iopa(n)==0) Go To 130
    k = iopa(n)/100
    i = iopa(n) - 100*k

    If (i==nions(k)) Then
      Write (999, *) ' STOP: ERROR IN NIONS OR ISI'
      Stop ' ERROR IN NIONS OR ISI'
    End If
    If (i<=0) Then
      Write (999, *) ' STOP: ERROR IN ION - FRESCAL '
      Stop ' ERROR IN ION - FRESCAL '
    End If

    levporg = 0
    Do ii = 1, id_rbftr
      ilab = labl4(ii)
      If (le(ilab)==k .And. li(ilab)==i) Then

!       JO Dec. 2015: avoid edges around Lya and Lyb
!       JO July 2021: and around Lyg and Lyd
!       JO Feb 2023: new approach: this is now allowed, will be accounted for
!       in FRESFIN
        If (frecin(ii)>lya_lo .And. frecin(ii)<lya_hi) Then
          Print *, ' WARNING!!!! WARNING!!!! WARNING!!!!'
          Print *, labl(ilab), &
            ' EDGE CLOSE TO LY-ALPHA NEGLECTED IN FREQ. GRID', k, ' ', i
          Write (999, *) ' WARNING!!!! WARNING!!!! WARNING!!!!'
          Write (999, *) labl(ilab), &
            ' EDGE CLOSE TO LY-ALPHA NEGLECTED IN FREQ. GRID', k, ' ', i
!         STOP ' WHEN YOU SEE THIS STATEMENT, PLEASE CONTACT J.P.'
!         after tests, the stop statement can be exchanged with the cycle
!         statement
          Cycle
        End If

        If (frecin(ii)>lyb_lo .And. frecin(ii)<lyb_hi) Then
          Print *, ' WARNING!!!! WARNING!!!! WARNING!!!!'
          Print *, labl(ilab), &
            ' EDGE CLOSE TO LY-BETA NEGLECTED IN FREQ. GRID', k, ' ', i
          Write (999, *) ' WARNING!!!! WARNING!!!! WARNING!!!!'
          Write (999, *) labl(ilab), &
            ' EDGE CLOSE TO LY-BETA NEGLECTED IN FREQ. GRID', k, ' ', i
!         STOP ' WHEN YOU SEE THIS STATEMENT, PLEASE CONTACT J.P.'
!         after tests, the stop statement can be exchanged with the cycle
!         statement
          Cycle
        End If

        If (frecin(ii)>lyg_lo .And. frecin(ii)<lyg_hi) Then
          Print *, ' WARNING!!!! WARNING!!!! WARNING!!!!'
          Print *, labl(ilab), &
            ' EDGE CLOSE TO LY-GAMMA NEGLECTED IN FREQ. GRID', k, ' ', i
          Write (999, *) ' WARNING!!!! WARNING!!!! WARNING!!!!'
          Write (999, *) labl(ilab), &
            ' EDGE CLOSE TO LY-GAMMA NEGLECTED IN FREQ. GRID', k, ' ', i
!         STOP ' WHEN YOU SEE THIS STATEMENT, PLEASE CONTACT J.P.'
!         after tests, the stop statement can be exchanged with the cycle
!         statement
          Cycle
        End If

        If (frecin(ii)>lyd_lo .And. frecin(ii)<lyd_hi) Then
          Print *, ' WARNING!!!! WARNING!!!! WARNING!!!!'
          Print *, labl(ilab), &
            ' EDGE CLOSE TO LY-DELTA NEGLECTED IN FREQ. GRID', k, ' ', i
          Write (999, *) ' WARNING!!!! WARNING!!!! WARNING!!!!'
          Write (999, *) labl(ilab), &
            ' EDGE CLOSE TO LY-DELTA NEGLECTED IN FREQ. GRID', k, ' ', i
!         STOP ' WHEN YOU SEE THIS STATEMENT, PLEASE CONTACT J.P.'
!         after tests, the stop statement can be exchanged with the cycle
!         statement
          Cycle
        End If

        levporg = levporg + 1
        flaux(levporg) = frecin(ii)

      End If
    End Do

!   reordering of flaux since fl maybe not ordered by energy
!   example: hei

    levtot = levtot + levporg
    Call sort(levporg, flaux)

!   find which edges shall be resolved

!   levmax levels for main ions. stage of elem. k1,k2

    leve = min(levporg, levmax)
    leve1 = min(leve, levmin1)

    If (imaa(k)-imia(k)>=2) Then !  NEW, at least 3 ions present
      leve2 = min(leve, levmin2)
    Else
      leve2 = leve1
    End If


!   levmin1 levels for 2nd important ions. stage of elem. k1,k2

    If (n>2) leve = leve1

!   levmin2 levels for rest

    If (n>nimp) leve = leve2

    levp = levporg
    If (levp==0) Then
      Write (999, *) ' STOP: ERROR IN LEVP - FRESMIN '
      Stop ' ERROR IN LEVP - FRESMIN '
    End If
!   levp1 = min(levmin1, levp) !    not used in the following

    Do l = 1, levp
      ifre = ifre + 1
      lli = levporg + 1 - l
      Do ii = 1, id_rbftr
        If (flaux(lli)==frecin(ii)) Go To 120
      End Do
      Write (999, *) ' STOP: FLAUX NOT FOUND IN FRECIN'
      Stop ' FLAUX NOT FOUND IN FRECIN'

120   Continue
!     K number of element
!     I number of ion
!     II number of transition
!     all w.r.t. DATA file (MAINION and ENIONND w.r.t enumeration as in data
!     file)
      index(ifre) = 100000*k + i*1000 + ii

!     here L counts from the highest transition freq. on.
!     If there is only ground state ionization, L corresponds to level
!     (starting at 1).
!     For ionization to excited states, L does *not* correspond to level
!     NOTE: INDEDGE modified below, if required
      If (l<=leve) Then
        indedge(ifre) = 1
      Else
        indedge(ifre) = 0
      End If

!     note that flaux here is in the "wrong' order

      fre(ifre) = flaux(lli)

    End Do

130 End Do

  If (levtot>id_rbftr) Then
    Write (999, *) ' STOP: ERROR IN LEVTOT'
    Stop ' ERROR IN LEVTOT'
  End If

! IF (IFRE.GT.ID_LLEVS) STOP  ! this was the old statement with one edge per
! level
  If (ifre>id_rbftr) Then
    Write (999, *) ' STOP: ERROR IN IFRE' ! THIS IS THE NEW ONE
    Stop ' ERROR IN IFRE' !         THIS IS THE NEW ONE
  End If

  Print *, ' NUMBER OF CONSIDERED EDGES = ', ifre
  Print *

  lam_ly = 0.
  Do i = 1, ifre
    eps = epslon(fre(i))
    fre(i+ifre) = fre(i) - eps
!   RED SIDE OF LYMAN-JUMP
!   IF(INDEX(I).EQ.101001) LAM_LY=1.D8/FRE(I+IFRE)
!   above statement not general, new version
    il = index(i)
    k = il/100000
    il = il - 100000*k
    io = il/1000
    ll = il - 1000*io
    lab = labl(labl4(ll))
    pare = paren4(ll)
    If (lab=='H11' .And. pare=='H21') lam_ly = 1.D8/fre(i+ifre)
    freedge(i) = fre(i)
  End Do

  if1 = ifre
  ifre = 2*ifre

! K-shell edges
  nk = 0

  If (optxray) Then
!   only edges from ionization stage III on
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
  End If

  Call sort(ifre, fre)

  fremin = fre(1)
  i = log10(fremin) - 1
  i = max(i, 0)
  i = 10**i
  fremin = float(int(fremin/i)*i)

  If (fre(1)<=fremin) Then
    Write (999, *) ' STOP: ERROR IN FREMIN'
    Stop ' ERROR IN FREMIN'
  End If

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

  If (optxray) Then
    aaa = 2.D7 !                    (5 A)

  Else
    emax = 10.D0*tmax/hkl

    If (nat==1 .And. k1==1 .Or. teff<10000.) Then
      aaa = 4.D5 !                  (250 A)
    Else If (teff<20000.) Then
      aaa = 1.D8/180.D0 !           (180 A)
    Else If (teff<35000.) Then
      aaa = 1.D6 !                  (100 A)
    Else
      aaa = 5.D6 !                  ( 20 A)
    End If

  End If

  Print *
  Print *, ' INTERMEDIATE FRE(IFRE): ', 1.0D8/fre(ifre)
  Print *, ' MIN      : ', 1.0D8/aaa
  Print *

  emax = max(emax, aaa)
  emax = max(emax, fre(ifre)*1.2D0)

  If (emax>7.5D5) Then !            in case, resolve HeII edge
    ifre = ifre + 1
    fre(ifre) = 7.5D5 !             (133 A)
  End If

  If (optxray) Then
    ifre = ifre + 1
    fre(ifre) = 1.2D6 !             (86 A, half way between 133 A and first
!   K-shell at 39A)
    ifre = ifre + 1
    fre(ifre) = 5.D6 !              ( 20 A)
  End If

  ifre = ifre + 1
  fre(ifre) = emax

! here was a bug, missing sort statement until Oct. 2014
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
!     changed Feb 2015, for better resol between 2000 and 1600 of non-HHe
!     models
!     DMIN1 = DMIN*2
      dmin1 = dmin/4.
    Else If (fre(i+1)<1.D8/910.) Then ! > 910 A
      dmin1 = dmin/4.

!     minimum separation behind HeII edge should be approx. 10  A or 2 A (if
!     X-rays)

    Else If (fre(i+1)>440528.D0) Then ! < 227 A
      dmin2 = 10.D-8*fre(i+1) !     DELTA LAM / LAM
      If (optxray) dmin2 = 2.D-8*fre(i+1) ! DELTA LAM / LAM
      If (optxray .And. fre(i+1)>5.D6) dmin2 = 0.3
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
            If (optout .And. iconver==1) Write (*, Fmt=370) 1.D8/fre(ifre), &
              dnew1, dnew1/fre(ifre)
          End Do

        Else If (imp==1 .And. j==2) Then
!         Jo June 2025: this block has changed. Also here, we use the
!         same stepsize (DNEW/4.) as for J=1, in case of OPTCMF_FULL AND .NOT.
!         OPTXRAY
!         for edges below the HeII one, to obtain equidistant grid-points
!         which allows for a better interpolation of xj_cmf_coarse.
          If (optcmf_full .And. .Not. optxray .And. fre(i+1)>440528.D0) Then
!           < 227 A
            dnew1 = dnew/4.D0
            Do k = 1, 3
              frenew = fre(ifre) + dnew1
              ifre = ifre + 1
              fre(ifre) = frenew
              If (optout .And. iconver==1) Write (*, Fmt=370) 1.D8/fre(ifre), &
                dnew1, dnew1/fre(ifre)
            End Do
!           old version
          Else
            dnew1 = dnew/2.D0
            frenew = fre(ifre) + dnew1
            ifre = ifre + 1
            fre(ifre) = frenew
            If (optout .And. iconver==1) Write (*, Fmt=370) 1.D8/fre(ifre), &
              dnew1, dnew1/fre(ifre)
          End If
        End If

        If (j/=inew) Then
          frenew = fre(i) + j*dnew
          ifre = ifre + 1
          fre(ifre) = frenew
          If (optout .And. iconver==1) Write (*, Fmt=370) 1.D8/fre(ifre), &
            dnew1, dnew1/fre(ifre)
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


! JO Dec. 2015
! modify frequency points around Lya and Lyb so that consistency
! with pure HHe freq. grid (if HHe only, no action performed)
! points chosen in such a way that Lya and Lyb not completely centered,
! to allow for a compromise (otherwise self-shadowing too strong)

! in case, add points
  Do i = 1, ifre
    If (abs(fre(i)-lya_lo)<1.D-9) Go To 140
  End Do
  ifre = ifre + 1
  fre(ifre) = lya_lo
  If (optout .And. iconver==1) Write (*, Fmt=390) 1.D8/fre(ifre)

140 Do i = 1, ifre
    If (abs(fre(i)-lya_hi)<1.D-9) Go To 150
  End Do
  ifre = ifre + 1
  fre(ifre) = lya_hi
  If (optout .And. iconver==1) Write (*, Fmt=390) 1.D8/fre(ifre)

150 Do i = 1, ifre
    If (abs(fre(i)-lyb_lo)<1.D-9) Go To 160
  End Do
  ifre = ifre + 1
  fre(ifre) = lyb_lo
  If (optout .And. iconver==1) Write (*, Fmt=390) 1.D8/fre(ifre)

160 Do i = 1, ifre
    If (abs(fre(i)-lyb_hi)<1.D-9) Go To 170
  End Do
  ifre = ifre + 1
  fre(ifre) = lyb_hi
  If (optout .And. iconver==1) Write (*, Fmt=390) 1.D8/fre(ifre)

170 Continue

  Call sort(ifre, fre)

! -----------------------------------------------
! in case, remove points in between lya_lo,lya_hi
  Do i = 1, ifre
    If (abs(fre(i)-lya_lo)<1.D-9) Go To 180
  End Do
  Write (999, *) ' STOP: LYA_LO NOT FOUND IN FRESCAL'
  Stop ' LYA_LO NOT FOUND IN FRESCAL'

180 istart = i
  Do i = istart, ifre
    If (abs(fre(i)-lya_hi)<1.D-9) Go To 190
  End Do
  Write (999, *) ' STOP: LYA_HI NOT FOUND IN FRESCAL'
  Stop ' LYA_HI NOT FOUND IN FRESCAL'

190 iend = i

  If (optout .And. iconver==1) Then
    Do i = istart + 1, iend - 1
      Write (*, Fmt=400) 1.D8/fre(i)
    End Do
  End If

  nskip = iend - istart - 1
  If (nskip<0) Then
    Write (999, *) ' STOP: ERROR IN NSKIP (FRESCAL)'
    Stop ' ERROR IN NSKIP (FRESCAL)'
  End If

  If (nskip/=0) Then
    ifre = ifre - nskip
    Do i = istart + 1, ifre
      fre(i) = fre(i+nskip)
    End Do
  End If

! -----------------------------------------------
! in case, remove points in between lyb_lo,lyb_hi

  Do i = 1, ifre
    If (abs(fre(i)-lyb_lo)<1.D-9) Go To 200
  End Do
  Write (999, *) ' STOP: LYB_LO NOT FOUND IN FRESCAL'
  Stop ' LYB_LO NOT FOUND IN FRESCAL'

200 istart = i
  Do i = istart, ifre
    If (abs(fre(i)-lyb_hi)<1.D-9) Go To 210
  End Do
  Write (999, *) ' STOP: LYB_HI NOT FOUND IN FRESCAL'
  Stop ' LYB_HI NOT FOUND IN FRESCAL'

210 iend = i

  If (optout .And. iconver==1) Then
    Do i = istart + 1, iend - 1
      Write (*, Fmt=400) 1.D8/fre(i)
    End Do
  End If

  nskip = iend - istart - 1
  If (nskip<0) Then
    Write (999, *) ' STOP: ERROR IN NSKIP (FRESCAL)'
    Stop ' ERROR IN NSKIP (FRESCAL)'
  End If

  If (nskip/=0) Then
    ifre = ifre - nskip
    Do i = istart + 1, ifre
      fre(i) = fre(i+nskip)
    End Do
  End If

! check resolution for
! hydrogen resonance lines and hei singlet resonance line


! XMIN=1.D0/1030.D-8
  xmin = 1.D0/1025.D-8 !            Lyman-beta
  xmax = 1.D0/912.D-8

! following block changed Feb 2015 to improve resolution
! and to include Lyman beta (important for Halpha)
! JO Dec. 2015: around Lya and Lyb, nothing should happen now
! (predefined points used for all atomic models (hopefully)
  Do i = 1, ifre - 1

    If (fre(i)<=1.D0/1215.D-8 .And. fre(i+1)>1.D0/1215.D-8) Then ! LYMAN ALPHA
      dmin1 = 0.005D0 !             VERY HIGH RESOLUTION
!     leave this value here;
!     usually, not completely centered, to allow for a compromise (otherwise
!     self-shadowing too strong)
!     if freq. points at Lya as given used, no additional point should be
!     included
      dm = fre(i+1)/fre(i) - 1.D0
      If (dm<dmin1) Cycle

!     new treatment
    Else If (fre(i)<=xmin .And. fre(i+1)>xmin) Then ! LYMAN BETA
!     that was the condition in the 'erroneous' version 10.7-1
!     DMIN1=0.001D0  !EVEN HIGHER RESOLUTION
      dmin1 = 0.005D0 !             VERY HIGH RESOLUTION
!     leave this value here;
!     usually, not completely centered, to allow for a compromise (otherwise
!     self-shadowing too strong)
!     if freq. points at Lyb as given used, no additional point should be
!     included
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
      If (optout .And. iconver==1) Write (*, Fmt=380) 1.D8/fre(ifre), dnew, &
        dnew/fre(ifre)
    End Do

  End Do

! JO CHANGED March 2017
! ONE ADDITIONAL POINT EXACTLY AT HEII LY-ALPHA
! comment out next 3 lines when working with older models!!!!!
  ifre = ifre + 1
  fre(ifre) = 1.D8/303.797
  If (optout .And. iconver==1) Write (*, Fmt=410) 1.D8/fre(ifre)

  Call sort(ifre, fre)

! JO July 2021
! modify frequency points around Lyg and Lyd so that consistency
! with pure HHe freq. grid (if HHe only, no action performed)

  Do i = 1, ifre
    If (abs(fre(i)-lyg_lo)<1.D-9) Go To 220
  End Do
  ifre = ifre + 1
  fre(ifre) = lyg_lo
  If (optout .And. iconver==1) Write (*, Fmt=390) 1.D8/fre(ifre)

220 Do i = 1, ifre
    If (abs(fre(i)-lyg_hi)<1.D-9) Go To 230
  End Do
  ifre = ifre + 1
  fre(ifre) = lyg_hi
  If (optout .And. iconver==1) Write (*, Fmt=390) 1.D8/fre(ifre)

230 Do i = 1, ifre
    If (abs(fre(i)-lyd_lo)<1.D-9) Go To 240
  End Do
  ifre = ifre + 1
  fre(ifre) = lyd_lo
  If (optout .And. iconver==1) Write (*, Fmt=390) 1.D8/fre(ifre)

240 Do i = 1, ifre
    If (abs(fre(i)-lyd_hi)<1.D-9) Go To 250
  End Do
  ifre = ifre + 1
  fre(ifre) = lyd_hi
  If (optout .And. iconver==1) Write (*, Fmt=390) 1.D8/fre(ifre)

250 Continue

  Call sort(ifre, fre)

! -----------------------------------------------
! in case, remove points in between lyg_lo,lyg_hi
  Do i = 1, ifre
    If (abs(fre(i)-lyg_lo)<1.D-9) Go To 260
  End Do
  Write (999, *) ' STOP: LYG_LO NOT FOUND IN FRESCAL'
  Stop ' LYG_LO NOT FOUND IN FRESCAL'

260 istart = i
  Do i = istart, ifre
    If (abs(fre(i)-lyg_hi)<1.D-9) Go To 270
  End Do
  Write (999, *) ' STOP: LYG_HI NOT FOUND IN FRESCAL'
  Stop ' LYG_HI NOT FOUND IN FRESCAL'

270 iend = i

  If (optout .And. iconver==1) Then
    Do i = istart + 1, iend - 1
      Write (*, Fmt=400) 1.D8/fre(i)
    End Do
  End If

  nskip = iend - istart - 1
  If (nskip<0) Then
    Write (999, *) ' STOP: ERROR IN NSKIP (FRESCAL)'
    Stop ' ERROR IN NSKIP (FRESCAL)'
  End If

  If (nskip/=0) Then
    ifre = ifre - nskip
    Do i = istart + 1, ifre
      fre(i) = fre(i+nskip)
    End Do
  End If

! -----------------------------------------------
! in case, remove points in between lyd_lo,lyd_hi
  Do i = 1, ifre
    If (abs(fre(i)-lyd_lo)<1.D-9) Go To 280
  End Do
  Write (999, *) ' STOP: LYD_LO NOT FOUND IN FRESCAL'
  Stop ' LYD_LO NOT FOUND IN FRESCAL'

280 istart = i
  Do i = istart, ifre
    If (abs(fre(i)-lyd_hi)<1.D-9) Go To 290
  End Do
  Write (999, *) ' STOP: LYD_HI NOT FOUND IN FRESCAL'
  Stop ' LYG_HI NOT FOUND IN FRESCAL'

290 iend = i

  If (optout .And. iconver==1) Then
    Do i = istart + 1, iend - 1
      Write (*, Fmt=400) 1.D8/fre(i)
    End Do
  End If

  nskip = iend - istart - 1
  If (nskip<0) Then
    Write (999, *) ' STOP: ERROR IN NSKIP (FRESCAL)'
    Stop ' ERROR IN NSKIP (FRESCAL)'
  End If

  If (nskip/=0) Then
    ifre = ifre - nskip
    Do i = istart + 1, ifre
      fre(i) = fre(i+nskip)
    End Do
  End If

  Call sort(ifre, fre)

! -----------------------------------------------
! in case, add points at center of important Stark-broadened lines (k ge 6)
  If (nstark/=0) Then

outer: Do ii = 1, nstark !          not executed for nstark=0

!     DO NOT TOUCH INTERVALS FOR HYDROGEN LY-LINES
      kayser = 1.D8/wavestark(ii)
      If ((kayser<=lya_hi .And. kayser>=lya_lo) .Or. &
        (kayser<=lyb_hi .And. kayser>=lyb_lo) .Or. (kayser<=lyg_hi .And. &
        kayser>=lyg_lo) .Or. (kayser<=lyd_hi .And. kayser>=lyd_lo)) &
        Cycle outer

      Do i = 1, ifre
        If (abs(fre(i)-kayser)<fre(i)*0.5D-3 .And. abs(wavestark( &
          ii)-lam_ciii)>1.D-4) Then ! corresponding to 0.5 A at 1000 A
!         PRINT*,'NEGLECTED',1.D8/FRE(I),WAVESTARK(II)
          Cycle outer
        End If
      End Do
      ifre = ifre + 1
      If (ifre>ifretot) Then
        Write (999, *) ' STOP: TOO MANY FREQUENCIES!'
        Stop ' TOO MANY FREQUENCIES!'
      End If
      fre(ifre) = kayser
      If (optout .And. iconver==1) Write (*, Fmt=420) 1.D8/fre(ifre)
    End Do outer

    Call sort(ifre, fre)
  End If

! JO Sept 2023 couple of additional frequencies from fgrid (central values)
! included, to allow for similar line-irradiation when different foreground
! elements are used. To this end, FGRID_LINES already called here (once),
! preparation for this routine (xion, vth(8)) in subr. START.
  If (optmet) Then
    If (start) Then
      Call fgrid_lines(teff)
      start = .False.
    End If

!   Add additional frequencies corresponding to rest wavelengths of explicit
!   lines,
!   and transfer wavblue,wavred,wavcon from nlte_approx to nlte
    Call frescal_expl_lines(wwavblue, wwavred, wwavcon)
    If (ifre+ifre_lines>ifretot) Then
      Write (999, *) ' STOP: TOO MANY FREQUENCIES (IFRE_LINES)'
      Stop ' TOO MANY FREQUENCIES (IFRE_LINES)'
    End If
    Do i = 1, ifre_lines
      ifre = ifre + 1
      fre(ifre) = fre_lines(i)
      If (optout .And. iconver==1) Write (*, Fmt=430) 1.D8/fre(ifre)
    End Do
    Deallocate (fre_lines)
    Call sort(ifre, fre)
  End If

  Print *, ' TOTAL NUMBER OF FREQUENCIES = ', ifre
  Print *
  Print *, ' FINAL FRE(IFRE): ', 1.0D8/fre(ifre)
  Print *


  If (optout .And. iconver==1) Then
    Print *, ' WAVELENGTH(A)     FREQ(HZ)'
    Do i = 1, ifre

      Do j = 1, if1
        If (fre(i)==freedge(j)) Go To 300
      End Do

      If (optxray) Then
        Do j = 1, n_kedges
          If (fre(i)==eth(j)) Then
            Write (*, Fmt=360) 1.D8/fre(i), clight*fre(i), name(z(j)), &
              nionk(j) - 1
            Go To 310
          End If
        End Do
      End If

      Write (*, Fmt=330) 1.D8/fre(i), clight*fre(i)
      Go To 310

300   Continue
      il = index(j)
      k = il/100000
      iz = int(zeff(k))
      il = il - 100000*k
      io = il/1000
      ll = il - 1000*io
      lab = labl(labl4(ll))
      pare = paren4(ll)
      If (indedge(j)/=1) Then
        Write (*, Fmt=340) 1.D0/fre(i)*1.D8, clight*fre(i), labat(k), io, &
          iz + io - 1, lab, pare
      Else
        Write (*, Fmt=350) 1.D0/fre(i)*1.D8, clight*fre(i), labat(k), io, &
          iz + io - 1, lab, pare
      End If

310 End Do
    Print *
  End If

  Do k = 1, ifre
    bn(k) = bnue(1.D8/fre(k), tmax)
  End Do

  pf = xintfre(bn, wfre, 1, ifre, ifre, fre, 4, 'NEW')
  pfe = sigsb/pi*tmax**4
  errp = abs(1.D0-pf/pfe)
  Print *, ' ERROR IN PLANCK-FUNCTION WITH ACTUAL FREQ. GRID = ', errp
  Print *

! JO Sept. 2021
! check whether frequency grid has changed (points can be the same,
! but frequencies different). Important to reset OP_FLAG
  newfreq = .False.

  If (ifreold==0) Then
    newfreq = .True.
    ifreold1 = ifre
    Print *, 'IFREOLD1 2:', ifreold1
    Do i = 1, ifre
      freold(i) = fre(i)
    End Do

  Else If (ifre/=ifreold) Then
    newfreq = .True.
  Else
    Do k = 1, ifre
      If (abs(1.D0-fre(k)/freold(k))>1.D-15) Then
        newfreq = .True.
        Exit
      End If
    End Do
  End If

  Print *
  If (newfreq) Then
    Print *, ' MESSAGE FROM FRESCAL: FREQ. GRID HAS CHANGED!'
  Else
    Print *, ' MESSAGE FROM FRESCAL: FREQ. GRID HAS NOT CHANGED!'
  End If
  Print *

! stop
  Return

320 Format (1X, A2, I2, ' Z = ', I2, ' MAX.NO.OF RESOLVED EDGES = ', I2)
330 Format (1X, F12.3, 6X, E12.6)
340 Format (1X, F12.3, 6X, E12.6, 3X, A2, I2, ' Z=', I2, '  FROM ', A6, &
    ' TO ', A6)
350 Format (1X, F12.3, 6X, E12.6, 3X, A2, I2, ' Z=', I2, '  FROM ', A6, &
    ' TO ', A6, ' IMPORTANT')
360 Format (1X, F12.3, 6X, E12.6, 3X, A2, ' Z=', I2, ' K-SHELL ABSORPTION')
370 Format (' ADDIT. POINT AT ', F12.3, 6X, F12.3, ' DELTA E / E = ', F5.3)
380 Format (' ADDIT. POINT AT ', F12.3, 6X, F12.3, ' DELTA E / E = ', F5.3, &
    ' (RESONANCE LINES!)')
390 Format (' ADDIT. POINT AT ', F12.3, ' LY_ALPHA/LY_BETA/LY_GAM/LY_DELTA')
400 Format ('        POINT AT ', F12.3, &
    ' DROPPED (LY_ALPHA/LY_BETA/LY_GAM/LY_DELTA)')
410 Format (' ADDIT. POINT AT ', F12.3, ' HEII LY_ALPHA')
420 Format (' ADDIT. POINT AT ', F12.3, &
    ' TRANSITION WITH STARK TREATMENT IN SUMOPAL')
430 Format (' ADDIT. POINT AT ', F12.3, ' TRANSITION FROM EXPLICIT ELEMENT')
End Subroutine

!-----------------------------------------------------------------------

Function xintfre(quan, w, ilow, imax, n, x, iopt, char)

  Use :: nlte_type
  Use :: nlte_dim
  Implicit None


! integrates quan over frequency from ilow to imax

! char = 'new'..... new integration weights calculated in frewe
! char = 'old'..... old integration weights used

! x.... input for frewe as described there
! opt.. input for frewe describing x


! .. scalar arguments ..
  Integer (i4b) :: ilow, imax, iopt, n
  Character :: char*3
! ..
! .. array arguments ..
  Real (dp) :: quan(n), w(n), x(n), xintfre
! ..
! .. local scalars ..
  Real (dp) :: summ
  Integer (i4b) :: i
! ..
! .. external subroutines ..
  External :: frewe
! ..

  If (char=='NEW') Then
    Call frewe(x, w, n, iopt)
  Else If (char/='OLD') Then
    Write (999, *) ' STOP: WRONG OPTION FOR CHAR IN INTFRE'
    Stop ' WRONG OPTION FOR CHAR IN INTFRE'
  End If

  summ = 0.D0
  Do i = ilow, imax
    summ = summ + quan(i)*w(i)
  End Do

  xintfre = summ

  Return

End Function

!-----------------------------------------------------------------------

Subroutine frewe(x, w, n, iopt)

  Use :: nlte_type
  Use :: nlte_dim
  Use :: fund_const, Only: clight
  Use :: nlte_var, Only: iconver, precis

  Implicit None

! calculates integration weights: input x, dim(x) = n  < 2500

! iopt=1 frequencies
! iopt=2 wavelengths (a)
! iopt=3 wavelengths (cm)
! iopt=4 kayser (1/cm)

! x can be ordered from high to low or vice versa

! output w, frequency integration weights such that

! integral(x-nue dnue) = sum over i(x-i*w-i)

! .. scalar arguments ..
  Integer (i4b) :: iopt, n
! ..
! .. array arguments ..
  Real (dp) :: w(n), x(n)
! ..
! .. local scalars ..
  Real (dp) :: delnue, delnue1, summ
  Integer (i4b) :: i, ii, k, kint, kk, nint, nint2, nint3, nint4, nint5, &
    nint6, nint7, nint8
  Logical :: optinv, optout
! ..
! .. local arrays ..
  Real (dp) :: del(4000), xnue(4000)
  Integer (i4b) :: index(4000)
! ..
! .. external subroutines ..
  External :: frewe2, frewe3, frewe4, frewe5, frewe6, frewe7, frewe8
! ..
! .. intrinsic functions ..
! INTRINSIC ABS,MOD
! ..

  If (n>4000) Then
    Write (999, *) ' STOP: N > 4000 IN FREWE'
    Stop ' N > 4000 IN FREWE'
  End If

  nint2 = 0
  nint3 = 0
  nint4 = 0
  nint5 = 0
  nint6 = 0
  nint7 = 0
  nint8 = 0

  optout = .True.

  optinv = .False.

  If (iopt==1 .Or. iopt==4) Then
    If (x(1)<x(n)) Then
      optinv = .True.
      Do i = 1, n
        xnue(i) = x(n+1-i)
      End Do
    Else
      Do i = 1, n
        xnue(i) = x(i)
      End Do
    End If
    If (iopt==4) Then
      Do i = 1, n
        xnue(i) = clight*xnue(i)
      End Do
    End If

  Else If (iopt==2 .Or. iopt==3) Then
    If (x(1)>x(n)) Then
      optinv = .True.
      Do i = 1, n
        xnue(i) = x(n+1-i)
      End Do
    Else
      Do i = 1, n
        xnue(i) = x(i)
      End Do
    End If
    If (iopt==2) Then
      Do i = 1, n
        xnue(i) = clight*1.D8/xnue(i)
      End Do
    Else
      Do i = 1, n
        xnue(i) = clight/xnue(i)
      End Do
    End If

  Else
    Write (999, *) ' STOP: WRONG INPUT OPTION'
    Stop ' WRONG INPUT OPTION'

  End If


  k = 1
  i = 1
  index(1) = 1

100 Continue
  delnue = xnue(i) - xnue(i+1)
  If (i+1==n) Then
    del(k) = delnue
    k = k + 1
    index(k) = n
    Go To 120
  End If

  Do ii = i + 1, n - 1
    delnue1 = xnue(ii) - xnue(ii+1)
    If (abs(1.D0-delnue/delnue1)>1.D-5) Go To 110
  End Do

  del(k) = delnue
  k = k + 1
  index(k) = n
  Go To 120

110 Continue
  del(k) = delnue
  k = k + 1
  index(k) = ii
  i = ii
  Go To 100

120 Continue
  kint = k
  If (optout .And. iconver==1) Then
    Print *
    Print *, ' NUMBER OF EQUIDISTANT SUBINTERVALS = ', kint - 1
    Print *
  End If

  w(1) = 0.D0
  Do kk = 2, kint
    nint = index(kk) - index(kk-1) + 1

    If (nint==1) Then
      Write (999, *) ' STOP: WRONG NUMBER OF SUBINTERVALS'
      Stop ' WRONG NUMBER OF SUBINTERVALS'

    Else If (nint==2) Then
      nint2 = nint2 + 1
      Call frewe2(index(kk-1), index(kk), del(kk-1), w, n)

    Else If (nint==3) Then
      nint3 = nint3 + 1
      Call frewe3(index(kk-1), index(kk), del(kk-1), w, n)

    Else If (nint==4) Then
      nint4 = nint4 + 1
      Call frewe4(index(kk-1), index(kk), del(kk-1), w, n)

    Else If (mod(nint,4)==1) Then
      nint5 = nint5 + 1
      Call frewe5(index(kk-1), index(kk), del(kk-1), w, n)

    Else If (mod(nint,4)==2) Then
      nint6 = nint6 + 1
      Call frewe6(index(kk-1), index(kk), del(kk-1), w, n)

    Else If (mod(nint,4)==3) Then
      nint7 = nint7 + 1
      Call frewe7(index(kk-1), index(kk), del(kk-1), w, n)

    Else If (mod(nint,4)==0) Then
      nint8 = nint8 + 1
      Call frewe8(index(kk-1), index(kk), del(kk-1), w, n)

    End If
  End Do

  If (optout .And. iconver==1) Then
    Print *
    Print *, ' NUMBER OF              TRAPEZ INTEGRATIONS = ', nint2
    Print *, ' NUMBER OF             SIMPSON INTEGRATIONS = ', nint3
    Print *, ' NUMBER OF                 3/8 INTEGRATIONS = ', nint4
    Print *, ' NUMBER OF        5-PT-SIMPSON INTEGRATIONS = ', nint5
    Print *, ' NUMBER OF 5-PT-SIMPSON+TRAPEZ INTEGRATIONS = ', nint6
    Print *, ' NUMBER OF        7-PT-SIMPSON INTEGRATIONS = ', nint7
    Print *, ' NUMBER OF 7-PT-SIMPSON+TRAPEZ INTEGRATIONS = ', nint8
    Print *
  End If

  summ = 0.D0
  Do k = 1, n
    summ = summ + w(k)
  End Do

  delnue = xnue(1) - xnue(n)
  If (abs(1.D0-delnue/summ)>precis) Then
    Write (999, *) ' STOP: ERROR IN FREWE!'
    Stop ' ERROR IN FREWE!'
  End If

  If (.Not. optinv) Return

  Do k = 1, n
    xnue(k) = w(n+1-k)
  End Do

  Do k = 1, n
    w(k) = xnue(k)
  End Do

  Return

End Subroutine

!-----------------------------------------------------------------------

Subroutine frewe2(i1, i2, del, w, n)

  Use :: nlte_type
  Use :: nlte_dim
  Implicit None

! trapezoidal rule


! .. scalar arguments ..
  Real (dp) :: del
  Integer (i4b) :: i1, i2, n
! ..
! .. array arguments ..
  Real (dp) :: w(n)
! ..

  w(i1) = w(i1) + .5D0*del
  w(i2) = .5D0*del

  Return

End Subroutine

!-----------------------------------------------------------------------

Subroutine frewe3(i1, i2, del, w, n)

  Use :: nlte_type
  Use :: nlte_dim
  Implicit None

! simpson's rule


! .. scalar arguments ..
  Real (dp) :: del
  Integer (i4b) :: i1, i2, n
! ..
! .. array arguments ..
  Real (dp) :: w(n)
! ..

  w(i1) = w(i1) + del/3.D0
  w(i1+1) = 4.D0*del/3.D0
  w(i2) = del/3.D0

  Return

End Subroutine

!-----------------------------------------------------------------------

Subroutine frewe4(i1, i2, del, w, n)

  Use :: nlte_type
  Use :: nlte_dim
  Implicit None

! 3/8 rule


! .. scalar arguments ..
  Real (dp) :: del
  Integer (i4b) :: i1, i2, n
! ..
! .. array arguments ..
  Real (dp) :: w(n)
! ..

  w(i1) = w(i1) + 3.D0*del/8.D0
  w(i1+1) = 9.D0*del/8.D0
  w(i1+2) = 9.D0*del/8.D0
  w(i2) = 3.D0*del/8.D0

  Return

End Subroutine

!-----------------------------------------------------------------------

Subroutine frewe5(i1, i2, del, w, n)

  Use :: nlte_type
  Use :: nlte_dim
  Implicit None

! modified simpson's rule


! .. scalar arguments ..
  Real (dp) :: del
  Integer (i4b) :: i1, i2, n
! ..
! .. array arguments ..
  Real (dp) :: w(n)
! ..
! .. local scalars ..
  Integer (i4b) :: k, nint, nk
! ..
! .. intrinsic functions ..
! INTRINSIC MOD
! ..

  nint = i2 - i1 + 1
  Do k = i1, i2
    nk = k - i1 + 1

    If (nk==1) Then
      w(k) = w(k) + (1.D0-1.D0/15.D0)*del/3.D0

    Else If (nk==nint) Then
      w(k) = (1.D0-1.D0/15.D0)*del/3.D0

    Else If (mod(nk,2)==0) Then
      w(k) = 4.D0*(1.D0+1.D0/15.D0)*del/3.D0

    Else If (mod(nk,4)==1) Then
      w(k) = 2.D0*(1.D0-1.D0/15.D0)*del/3.D0

    Else
      w(k) = 1.6D0*del/3.D0

    End If
  End Do

  Return

End Subroutine

!-----------------------------------------------------------------------

Subroutine frewe6(i1, i2, del, w, n)

  Use :: nlte_type
  Use :: nlte_dim
  Implicit None

! modified simpson's rule + trapezoidal rule for last interval


! .. scalar arguments ..
  Real (dp) :: del
  Integer (i4b) :: i1, i2, n
! ..
! .. array arguments ..
  Real (dp) :: w(n)
! ..
! .. local scalars ..
  Integer (i4b) :: k, nint, nk
! ..
! .. intrinsic functions ..
! INTRINSIC MOD
! ..

  nint = i2 - i1 + 1

  Do k = i1, i2
    nk = k - i1 + 1

    If (nk==1) Then
      w(k) = w(k) + (1.D0-1.D0/15.D0)*del/3.D0

    Else If (nk==nint-1) Then
      w(k) = (1.D0-1.D0/15.D0)*del/3.D0 + del/2.D0

    Else If (nk==nint) Then
      w(k) = del/2.D0

    Else If (mod(nk,2)==0) Then
      w(k) = 4.D0*(1.D0+1.D0/15.D0)*del/3.D0

    Else If (mod(nk,4)==1) Then
      w(k) = 2.D0*(1.D0-1.D0/15.D0)*del/3.D0

    Else
      w(k) = 1.6D0*del/3.D0

    End If
  End Do

  Return

End Subroutine

!-----------------------------------------------------------------------

Subroutine frewe7(i1, i2, del, w, n)

  Use :: nlte_type
  Use :: nlte_dim
  Implicit None

! simpson's rule


! .. scalar arguments ..
  Real (dp) :: del
  Integer (i4b) :: i1, i2, n
! ..
! .. array arguments ..
  Real (dp) :: w(n)
! ..
! .. local scalars ..
  Integer (i4b) :: k, nint, nk
! ..
! .. intrinsic functions ..
! INTRINSIC MOD
! ..

  nint = i2 - i1 + 1

  Do k = i1, i2
    nk = k - i1 + 1
    If (nk==1) Then
      w(k) = w(k) + del/3.D0

    Else If (nk==nint) Then
      w(k) = del/3.D0

    Else If (mod(nk,2)==0) Then
      w(k) = 4.D0*del/3.D0

    Else
      w(k) = 2.D0*del/3.D0

    End If
  End Do

  Return

End Subroutine

!-----------------------------------------------------------------------

Subroutine frewe8(i1, i2, del, w, n)

  Use :: nlte_type
  Use :: nlte_dim
  Implicit None

! simpson's rule + trapezoidal rule for last interval


! .. scalar arguments ..
  Real (dp) :: del
  Integer (i4b) :: i1, i2, n
! ..
! .. array arguments ..
  Real (dp) :: w(n)
! ..
! .. local scalars ..
  Integer (i4b) :: k, nint, nk
! ..
! .. intrinsic functions ..
! INTRINSIC MOD
! ..

  nint = i2 - i1 + 1

  Do k = i1, i2
    nk = k - i1 + 1

    If (nk==1) Then
      w(k) = w(k) + del/3.D0

    Else If (nk==nint-1) Then
      w(k) = del/3.D0 + del/2.D0

    Else If (nk==nint) Then
      w(k) = del/2.D0

    Else If (mod(nk,2)==0) Then
      w(k) = 4.D0*del/3.D0

    Else
      w(k) = 2.D0*del/3.D0

    End If
  End Do

  Return

End Subroutine

!-----------------------------------------------------------------------

Subroutine nistar

  Use :: nlte_type
  Use :: nlte_dim


  Use :: princesa_var, Only: nl, nat, index1, index3, index4, ixa1, kl, labl, &
    labl1, labl3, labl4, labu1, le, paren3, paren4

  Use :: nlte_var, Only: mlow, mup, mmlow, mmup, mclow, mcup, indr, indc

  Implicit None

! ------this subroutine calculates (nlow/nup)*
! hay que definir tambien las matrices de /atomdat/

! .. parameters ..
! Integer (i4b), Parameter :: kel = id_atoms, kis = id_kisat
! ..
! .. local scalars ..
  Integer (i4b) :: i, ihi, ilo, iup, jj, k, mclowa, mcupa, mlowa, mmlowa, mupa
  Character :: lupp*6
! ..

  Do i = 1, index4

    lupp = paren4(i)
    ilo = labl4(i)
    mlowa = 0
    k = le(ilo)
    Do jj = 1, k - 1
      mlowa = mlowa + ixa1(jj)
    End Do

    mlow(i) = mlowa + ilo
    ihi = kl(ilo)

100 Continue

    If (labl(ihi)==lupp) Then
      iup = ihi
      Go To 110
    End If

    If (ihi==nl) Then
      Write (999, *) ' STOP: UPPER LEVEL NOT FOUND IN RBF '
      Stop ' UPPER LEVEL NOT FOUND IN RBF '
    End If
    ihi = ihi + 1

    Go To 100

110 Continue

    k = le(iup)
    mupa = 0

    Do jj = 1, k - 1
      mupa = mupa + ixa1(jj)
    End Do

    mup(i) = mupa + iup

  End Do

! ------ introduced (4-oct-92) for bound-bound transitions treatment

  Do i = 1, index1

    ilo = labl1(i)
    mmlowa = 0
    k = le(ilo)

    Do jj = 1, k - 1
      mmlowa = mmlowa + ixa1(jj)
    End Do

    mmlow(i) = mmlowa + ilo
    iup = labu1(i)
    mmup(i) = mmlowa + iup
  End Do

! ------ for cbf transitions

  Do i = 1, index3

    ilo = labl3(i)
    If (ilo<0) Then
      Write (999, *) ' STOP: ILO < 0 IN LABL3'
      Stop ' ILO < 0 IN LABL3'
!     mclow(i)=0
!     mcup(i)=0
    Else
      mclowa = 0
      k = le(ilo)

      Do jj = 1, k - 1
        mclowa = mclowa + ixa1(jj)
      End Do

      mclow(i) = mclowa + ilo
      lupp = paren3(i)
      ihi = kl(ilo)

120   Continue

      If (labl(ihi)==lupp) Then
        iup = ihi
        Go To 130
      End If

      If (ihi==nl) Then
        Write (999, *) ' STOP: UPPER LEVEL NOT FOUND - CBF'
        Stop ' UPPER LEVEL NOT FOUND - CBF'
      End If
      ihi = ihi + 1
      Go To 120

130   Continue

      k = le(iup)
      mcupa = 0
      Do jj = 1, k - 1
        mcupa = mcupa + ixa1(jj)
      End Do

      mcup(i) = mcupa + iup

    End If

  End Do

  indr = index4
  indc = index3
  If (indr<(nl-nat)) Then
    Write (999, *) ' STOP: LEVELS NOT CONSIDERED IN PHOTO. '
    Stop ' LEVELS NOT CONSIDERED IN PHOTO. '
  End If

  Return

End Subroutine

!----------------------------------------------------------------------

Subroutine occlte(rho, xne, xnh, clf, te, ilow, imax, nd, ntemp, lte_update)
! ---- unaffected by OPTTHICK

  Use :: nlte_type
  Use :: nlte_dim
  Use :: fund_const
  Use :: princesa_var, Only: nat, ixa1, nions, ifirsl, ifirx, gl, fl, glx, &
    flx, labat

  Use :: nlte_var, Only: ionis, xnelte, occnum, alevel, blevel, enionnd, &
    abundnd, yhein, iconver, modnam, partit, optmet, optcmf_all

  Use :: tcorr_var, Only: enatcor, emaxtc, temp_converged

  Implicit None

! -------this subroutine calculates lte occupation numbers, &
! accounting for clumping with enhancement-factor clf and
! corrected values for xne and xnh

! .. parameters ..
  Integer (i4b), Parameter :: kel = id_atoms, kis = id_kisat
  Integer (i4b), Parameter :: nd1 = id_ndept
! ..
! .. scalar arguments ..
  Integer (i4b) :: nd, ntemp
! ..
! .. array arguments ..
  Real (dp) :: rho(nd1), te(nd1), xne(nd1), xnh(nd1), clf(nd1)
  Integer (i4b) :: ilow(nd1, kel), imax(nd1, kel)
! ..
! .. local scalars ..

  Real (dp) :: deltae, dumreal, e0, occ, optmixed, rstnom, teff, xlogg
  Integer (i4b) :: i, ihi, isi, j, k, k1, l, leo, levl, levx, ll, m, nene, &
    nnn, nrec
  Logical :: helone, optout, optneupdate, lte_update, dumlogic, optcmf1
  Character :: dumchar*30
! ..
! .. external subroutines ..
  External :: ground1
! ..
! .. intrinsic functions ..
! INTRINSIC ABS,EXP,MAX,MAX0,MIN,MIN0
! ..
  optout = .True.

  nrec = id_llevs

! check consistency
  If (lte_update .And. .Not. enatcor .And. emaxtc>0.003D0) Then
    Write (999, *) ' STOP: INCONSISTENCY IN OCCLTE'
    Stop ' INCONSISTENCY IN OCCLTE'
  End If

  Open (1, File=trim(modnam)//'/INDAT.DAT', Status='OLD')
  Rewind 1
  Read (1, Fmt=*) dumchar
  Read (1, Fmt=*) optneupdate, helone, nnn, nnn
  Read (1, Fmt=*) optmixed
  Read (1, Fmt=*) teff, xlogg, rstnom
  Read (1, Fmt=*) dumreal, dumreal
  Read (1, Fmt=*) dumreal, dumreal, dumreal
  Read (1, Fmt=*) yhein, dumreal
  Read (1, Fmt=*) dumlogic, dumlogic, dumlogic, dumlogic, optcmf1
  Close (1)

  If (optout .And. iconver==1) Then
    Write (*, Fmt=110) teff, xlogg, yhein
    Print *
  End If

! -----------------------------------------------------------------------

  nene = 0
  Do isi = 1, nat
    If (labat(isi)=='HE') nene = isi
    If (ixa1(isi)/=0) Then
      Write (999, *) ' STOP: X-LEVELS NOT YET IMPLEMENTED'
      Stop ' X-LEVELS NOT YET IMPLEMENTED'
    End If
  End Do

  Call ground1(te, rho, xne, xnh, clf, nd, optout, ntemp, yhein, nene)

! occmain maximum fraction of all ions EXCLUDING the highest one without
! transitions

! IF TEMP. ALMOST CONVERGED, DON'T CHANGE IMIN, IMAX ETC.
! IF(LTE_UPDATE.AND.EMAXTC.LT.0.05) THEN
! JO CHANGED: one final update after temperature has converged, to allow
! successful restart
  If (lte_update .And. emaxtc<0.03 .And. .Not. temp_converged) Then
!   at least, check for changes in IONIS
    Do k = 1, nat
      Do l = 1, nd
        If (imax(l,k)>ionis(k,l)) Then
          imax(l, k) = ionis(k, l)
          If (imax(l,k)<ilow(l,k)) Then
            Write (999, *) ' STOP: IMAX < ILOW'
            Stop ' IMAX < ILOW'
          End If
          Print *, 'CORRECTION OF IMAX(LTE) REQUIRED!!!'
          Print *, k, l, ' IMAX(NLTE) = ', imax(l, k), ' (ZEFF CORRECTED)'
        End If
      End Do
    End Do
    Go To 100
  End If

! JO Dec. 2017: when complete cmf-treatment has been switched on,
! ilow/imax(LTE) must not be changed to keep freq. grid (IFRE).
! Otherwise: stop 'IFRE HAS CHANGED IN CMF_COMPLETE!!!!' might occur.
! needs to be tested.

  If (optcmf_all) Go To 100

! ----------- determination of ilow and imax etc, for LTE calculations
! and H/He (always) in the old spirit
! Note the major change from version 9.0 on

  If (optmet) Then
    Call ilowimax_lte(ilow, imax, te, teff, nd, ntemp, helone, optcmf1, &
      optout, lte_update)
  Else
!   this is the very old version kept for consistency
!   in case optmet=.false., i.e., no background metals
!   (which are used to predict the correct ionization
!   equilibrium in the standard approach)
    Call ilowimax_old(ilow, imax, te, teff, nd, ntemp, helone, optcmf1, &
      optout)
  End If

! -------------calculation of n(i,j,k) = g_i *(n1/g1) * exp(-dE/kT)
! n1/g1=occnum at this stage (see ground1)

100 Continue
depthloop: Do ll = ntemp, nd

    j = 0
    m = 0

kloop: Do k = 1, nat

      k1 = nions(k) !               excitation for ALL levels
!     thus: excited levels of IONIS(K,L) populated
!     rest = 0.

ionloop: Do i = 1, k1 - 1
        j = j + 1
        levl = ifirsl(j+1) - ifirsl(j)
        levx = ifirx(j+1) - ifirx(j)
!       levs=ifirs(j+1)-ifirs(j)
        e0 = fl(ifirsl(j))

        Do leo = 0, levl - 1
          m = m + 1
          ihi = ifirsl(j) + leo
          deltae = abs(e0-fl(ihi))
          occ = occnum(k, i, ll)*exp(-1.D0*hk*deltae/te(ll))
          alevel(m) = occ*gl(ihi)
          blevel(m) = alevel(m)
        End Do

        Do leo = 0, levx - 1
          m = m + 1
          ihi = ifirx(j) + leo
          deltae = (e0-flx(ihi))
          occ = occnum(k, i, ll)*exp(-1.D0*hk*deltae/te(ll))
          alevel(m) = occ*glx(ihi)
          blevel(m) = alevel(m)
        End Do

!       do leo=o,levs
!       ihi=ifirs(j)+leo
!       do isi=ns0(ihi),ns1(ihi)
!       m=m+1
!       deltae=e0-ryd(ihi)/(isi-qd(ihi))**2
!       occ=occnum(k,i,ll)*exp(-1.d0*deltae/akb/te(ll))
!       alevel(m)=occ
!       blevel(m)=alevel(m)
!       end do
!       end do

      End Do ionloop

!     -----------k level

      j = j + 1
      m = m + 1


!     here was the bug (forgotten to multiply with partition function of
!     uppermost ion)!

      If (partit(k,k1,ll)/=gl(m)) Then
        Write (999, *) ' STOP: SOMETHING WRONG WITH LAST G-VALUE!'
        Stop ' SOMETHING WRONG WITH LAST G-VALUE!'
      End If
      occ = occnum(k, k1, ll)*gl(m)
      alevel(m) = occ
      blevel(m) = alevel(m)

    End Do kloop

    Write (15, Rec=ll)(alevel(i), i=1, nrec)
    Write (17, Rec=ll)(blevel(i), i=1, nrec)

    xnelte(ll) = xne(ll)

    If (iconver==1) Then
      Print *, ' L= ', ll, '   TOTAL NUMBER OF LEVELS= ', m

!     ---    calculation of depth dependent abundance (including clumping)

      Do k = 1, kel
        abundnd(k, ll) = 0.D0
        Do i = 1, kis + 1
          abundnd(k, ll) = abundnd(k, ll) + enionnd(k, i, ll)
        End Do
      End Do
    End If

  End Do depthloop

  Print *
  Print *, '  FILE OCCLTE WRITTEN '
  Print *

  Return

110 Format (10X, 'TEFF= ', F7.0, 'LOGG= ', F5.2, 'YHE= ', F5.2)

End Subroutine

!-----------------------------------------------------------------------

Subroutine ground1(te, rho, xne, xnhhyd, clf, nd, optout, ntemp, yhein, nene)
! ---- unaffected by OPTTHICK

  Use :: nlte_type
  Use :: nlte_dim
  Use :: fund_const
  Use :: princesa_var, Only: nat, nions, ifirsl, ifirx, ifirs, ns0, ns1, zeff, &
    weight, abund, gl, fl, zl, glx, flx, gs0, gs1, ryd, qd

  Use :: nlte_var, Only: ionis, partit, xl, occnum, mainion, enionnd, &
    enionnd_lte, iconver, precis

  Implicit None

! this subroutine calculates the lte population of ground level
! (divided by g1), accounting for clumping with clumping factor clf
! xne and xnh are corrected, rho is average (uncorrected) density


! ------- calculation of partition function

! .. parameters ..
  Integer (i4b), Parameter :: kel = id_atoms, kis = id_kisat
  Integer (i4b), Parameter :: nd1 = id_ndept
! ..
! .. scalar arguments ..
  Real (dp) :: yhein
  Integer (i4b) :: nd, nene, ntemp
  Logical :: optout
! ..
! .. array arguments ..
  Real (dp) :: rho(nd1), te(nd1), xne(nd1), xnhhyd(nd1), clf(nd1)
! ..
! .. local scalars ..
  Real (dp) :: ag, al, dev, e1, edenh, el, g, occmain, summass, sumxmue, xneh, &
    xnh, zi, zi1
  Integer (i4b) :: i, ihi, ii, iit, ij, ik, is, iz, j, jj, k, k1, l, l1, l2, &
    ll, lx, n, nionmax
  Logical :: conv, non
! ..
! .. local arrays ..
  Real (dp) :: enion(kel), u(nd1)
! ..
! .. external functions ..
  Integer (i4b) :: igenio
  External :: igenio
! ..
! .. intrinsic functions ..
! INTRINSIC ABS,DBLE,EXP,INT,MAX
! ..

  summass = 0.D0
  nionmax = 0

  Do i = 1, nat
    nionmax = max(nionmax, nions(i))
    summass = summass + weight(i)*abund(i)
  End Do

  If (nene==0 .And. yhein/=0.) Then
    summass = summass + 4.D0*yhein
  End If

  If (nionmax-1>kis) Then
    Write (999, *) ' STOP: KIS TOO SMALL, RE-DIMENSION IN DIMES!'
    Stop ' KIS TOO SMALL, RE-DIMENSION IN DIMES!'
  End If

  j = 0

kloop1: Do k = 1, nat
    ik = nions(k)
ionloop: Do jj = 1, ik - 1
      j = j + 1
      Do lx = ntemp, nd
        u(lx) = 0.D0
      End Do

      n = ifirsl(j)
      e1 = fl(n)
      ihi = ifirsl(j+1) - 1

      Do ll = ifirsl(j), ihi
        el = abs(fl(ll)-e1)
        Do lx = ntemp, nd
          u(lx) = u(lx) + gl(ll)*exp(-1.D0*hk*el/te(lx))
        End Do
      End Do

      ihi = ifirx(j+1) - 1
      Do ll = ifirx(j), ihi
        el = abs(flx(ll)-e1)
        Do lx = ntemp, nd
          u(lx) = u(lx) + glx(ll)*exp(-1.D0*hk*el/te(lx))
        End Do
      End Do

      ihi = ifirs(j+1) - 1

      Do ll = ifirs(j), ihi
        Do is = ns0(ll), ns1(ll)
          el = abs(abs(ryd(ll))/(dble(is)-qd(ll))**2) - e1*hh
          g = gs0(ll) + gs1(ll)*dble(is)**2
          Do lx = ntemp, nd
            u(lx) = u(lx) + g*exp(-1.D0*el/akb/te(lx))
          End Do
        End Do
      End Do

      partit(k, jj, ntemp:nd) = u(ntemp:nd)

    End Do ionloop

    j = j + 1
    ihi = ifirsl(j)

    Do lx = ntemp, nd
      partit(k, ik, lx) = gl(ihi)
    End Do

  End Do kloop1

! -------calculation of n(j,k)/n(h) !

  sumxmue = 0.D0

  iit = 0

depthloop: Do ll = ntemp, nd

    conv = .False.
    dev = 1.D0

100 Continue

    iit = iit + 1

!   -------at first, we calculate of l(j,k)

    i = 0
kloop2: Do k = 1, nat
      k1 = nions(k)
      non = .False.

      Do ii = 1, k1 - 1
        i = i + 1
        n = ifirsl(i)
        al = fl(n)
        ag = al*hk/te(ll)

        If ((ag>150D0) .And. .Not. non) Then
          ionis(k, ll) = ii - 1
          non = .True.
        End If
!       xne includes clumping
        xl(k, ii) = .5D0*saha*partit(k, ii, ll)*xne(ll)*exp(ag)/ &
          partit(k, ii+1, ll)/te(ll)**1.5D0
      End Do

      i = i + 1
      xl(k, k1) = 1.D0
      If (.Not. non) ionis(k, ll) = k1 - 1
    End Do kloop2

!   ------calculation of n(j,k)/n(h) = occnum(k,i,ll)

    xnh = rho(ll)*clf(ll)/amh/summass ! corrected for clumping
!   IF (ABS(1.D0-XNHHYD(LL)/XNH).GT.PRECIS) STOP 'ERROR IN XNH'
    If (abs(1.D0-xnhhyd(ll)/xnh)>0.01) Then
      Write (999, *) ' STOP: ERROR IN XNH' ! changed
      Stop ' ERROR IN XNH' !        changed
    End If

kloop3: Do k = 1, nat
      k1 = ionis(k, ll)
      al = 1.D0

      Do l = 1, k1
        ag = 1.D0

        Do l1 = l, k1
          ag = ag*xl(k, l1)
        End Do

        al = al + ag
      End Do

      enion(k) = abund(k)/al
      occnum(k, k1, ll) = xl(k, k1)*enion(k)

      Do l2 = k1 - 1, 1, -1
        occnum(k, l2, ll) = xl(k, l2)*occnum(k, l2+1, ll)
      End Do

    End Do kloop3

!   ------------- new electronic density, accounting for zeff

    edenh = 0.D0

kloop4: Do k = 1, nat
      k1 = ionis(k, ll)
      Do i = 1, k1
        zi = zeff(k) + i - 1

!       check of charge

        ij = igenio(k, i)
        zi1 = zl(ifirsl(ij))
        If (zi/=zi1) Then
          Write (999, *) ' STOP: ERROR IN CHARGE'
          Stop ' ERROR IN CHARGE'
        End If
        edenh = edenh + zi*occnum(k, i, ll)
      End Do

      zi = zeff(k) + k1
      ij = igenio(k, k1+1)
      zi1 = zl(ifirsl(ij))
      If (zi/=zi1) Then
        Write (999, *) ' STOP: ERROR IN CHARGE-HIGHEST ION'
        Stop ' ERROR IN CHARGE-HIGHEST ION'
      End If
      edenh = edenh + zi*enion(k)

    End Do kloop4

    If (conv) sumxmue = sumxmue + rho(ll)*clf(ll)/(edenh*xnh)
    xneh = xne(ll)/xnh
    dev = abs(xneh-edenh)/xneh
    xne(ll) = edenh*xnh

!   --------------calculation of mainion and output

    If (conv) Then

      Do k = 1, nat
        occmain = 0.D0
        k1 = ionis(k, ll)

        Do l = 1, k1
          enionnd(k, l, ll) = occnum(k, l, ll)*xnh
          occmain = max(occmain, occnum(k,l,ll))
          If (occmain==occnum(k,l,ll)) mainion(ll, k) = l
        End Do

        enionnd(k, k1+1, ll) = enion(k)*xnh
        enionnd(k, k1+2:kis, ll) = 0.D0 ! in case ionization changes
      End Do
      enionnd_lte = enionnd !       for metals in iterations before CONCON


      If (optout .And. iconver==1) Then
        Print *
        Print *, '-----------------------------------------'
        Print *
        Print *, '  L=', ll, '    TE=', te(ll), '    NE/NH=', xneh
        Print *
        Print *, ' NO(ELEM) NO(MAIN ION) Z(MAIN ION) // NO  Z  N(J,K)/N(H)  '
        Do k = 1, nat
          k1 = ionis(k, ll)
          iz = int(zeff(k))
          Write (*, Fmt=110) k, mainion(ll, k), (iz+mainion(ll,k)-1), &
            (l, (iz+l-1), occnum(k,l,ll), l=1, k1), k1 + 1, iz + k1, enion(k)
        End Do
      End If

    End If

!   ---------------calculation of n(j,k)/u(j,k)=occnum(k,i,ll) = n1/g1

kloop5: Do k = 1, nat
      k1 = ionis(k, ll)

      Do i = 1, k1
        occnum(k, i, ll) = occnum(k, i, ll)/partit(k, i, ll)*xnh
      End Do

      l1 = nions(k)
      occnum(k, k1+1, ll) = enion(k)/partit(k, k1+1, ll)*xnh

      Do i = k1 + 2, l1
        occnum(k, i, ll) = 0.D0
      End Do

    End Do kloop5

    If (conv) Cycle
    If (dev<precis) conv = .True.
    Go To 100

  End Do depthloop

! -----------------------------------------------------------------------

  sumxmue = sumxmue/amh/(nd+1-ntemp)
  iit = iit/(nd+1-ntemp)
  Print *
  Print *, ' MUE-ELECTRON= ', sumxmue
  Print *, ' AVERAGE NUMBER OF NE-ITERATIONS WAS ', iit
  Print *

  Return
110 Format (/, I3, 2(3X,'/',I3), /, 9(I3,3X,I3,3X,E12.6,/))

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
  If (igenio>iong) Then
    Write (999, *) ' STOP: ERROR IN GENERAL ION - IGENIO '
    Stop ' ERROR IN GENERAL ION - IGENIO '
  End If

  Return

End Function

!-----------------------------------------------------------------------

Subroutine opacitc(nd, xne, xnh, temp, clf, ilow, imax, lteopt, ntemp, opteta, &
  frenew, optmet)


! -------this subroutine calculates the free-free and bound-free
! opacity for every frequency at every depth point

! ---    note: standard detail input provides only ff opacities
! for h and he. e.g., if calculating wr's, take care for cno!!!

! CHANGED FROM VER 9.5
! ---    warning: opacities calculated for all transitions with
! n_l ne 0 (bf) or ENIONND ne 0 (ff).
! Occupation/ionic numbers with zero value should
! occur only  in the range outside ILOW / IMAX

! ilow and imax can be lte or nlte values, dependent of call

! THE FINAL OPACITIES, EMISSIVITES AND SOURCE-FUNCTIONS
! HAVE BEEN CORRECTED FOR CLUMPING, EXCEPT FFOPA WHICH IS NEEDED
! ONLY FOR THE CALCULATION OF THE T-STRUCTURE IN ORDER TO ENABLE
! CORRECTLY EVALUATED SPATIAL AVERAGES

! OPAC CHANGED FOR EFFECTS FROM OPTICALLY THICK CLUMPING, I.E.,
! OPAC IS NOW EFFECTIVE OPACITY
! ALL OTHER OPACITIES ETC. SHOULD REMAIN AS BEFORE (OPTICALLY THIN
! CLUMPING), I.E., SHOULD BE MEAN OPACITIES

! includes emissivities from hot plasma (wind-embedded shocks)

  Use :: nlte_type
  Use :: nlte_dim
  Use :: fund_const

  Use :: princesa_var, Only: kl, labl, labl4, le, li, zl, ifirsl, index4, &
    index5, iform4, iform5, frecin, frecfn, ipare5

  Use :: nlte_var, Only: alphafs, gffnlte, ionmax, izmin, izmax, xnelte, fre, &
    ifre, alevel, blevel, enionnd, mlow, mup, opac, strue, opat_m_old, &
    opat_m_new, strue_m_old, strue_m_new, opac_nolines, opat_m_nolines, &
    etat_nolines, etat_m_nolines, precis, yhe => yhein, op_flag, frecfn1
! ,RAYLEIGH

  Use :: nlte_xrays, Only: optxray, lambdanu, fxl, lxmin, opa_kshell

  Use :: tcorr_var, Only: ffopa, dtffopa

  Use :: nlte_porvor, Only: epsi, epsi1, optthick, fic, tcl_fac_cont, &
    opa_eff_rat, opa_eff_rat_old, w_bg_red, w_bg_blue

  Implicit None

! .. parameters ..
  Integer (i4b), Parameter :: kel = id_atoms ! , kis = id_kisat
! Integer (i4b), Parameter :: ifretot = id_frec1
  Integer (i4b), Parameter :: nd1 = id_ndept

! ..
! .. scalar arguments ..
  Integer (i4b) :: nd, ntemp
  Logical :: frenew, lteopt, opteta, optmet
! ..
! .. array arguments ..
  Real (dp) :: temp(nd1), xne(nd1), xnh(nd1), clf(nd1)

  Integer (i4b) :: ilow(nd1, kel), imax(nd1, kel)
! ..
! .. local scalars ..
  Real (dp) :: a, alpha, bn, bn1, chiff, chiffc, depart, ehnuekt, etab, etabf, &
    etaffc, freq, gf, hcfre3, opab, opabf, opaff, opaff1, opath, sig, xnerat, &
    xni, xnistar, xnk, z, opahminus, en, constx_ll, x, expo, tcl, lam, xnel

  Integer (i4b) :: i, ii, ij, iz, inato, kk, ll, nato, nlas, nm, nn, nrec, &
    kmin
! ..
! .. external functions ..
  Real (dp) :: gff
  Integer (i4b) :: igenio
  External :: gff, igenio
! ..
! .. external subroutines ..
  External :: crossbf
! ..
! .. intrinsic functions ..
! INTRINSIC ABS,EXP,SQRT
! ..

  nrec = id_llevs
  sig = amh*sigmae

  If (.Not. lteopt .And. optxray) Then
    Call opacit_kshell(nd, xne, ilow, imax)
    kmin = lbound(opa_kshell, dim=2)
    If (ubound(opa_kshell,dim=2)/=ifre) Then
      Write (999, *) ' STOP: ERROR IN UBOUND OF OPA_KSHELL'
      Stop ' ERROR IN UBOUND OF OPA_KSHELL'
    End If
  End If

  ffopa = 0.D0
  dtffopa = 0.D0

  If (frenew) Then !                IN LTE-ITERATION, FRENEW IS ALWAYS TRUE
    izmin = 100
    izmax = 0

!   Find out all possible charges to be treated

    Do i = 1, index5
      If (iform5(i)/=1) Then
!       ------        has to be implemented if x-levels are present
        Write (999, *) ' STOP: RFF NOT IMPLEMENTED; ', i, '!'
        Print *, 'RFF NOT IMPLEMENTED; ', i, '!'
        Stop
      End If

      ii = ipare5(i)
      z = zl(ii)
      izmin = min(izmin, int(z+1.D-6))
      izmax = max(izmax, int(z+1.D-6))
    End Do

    If (izmax>ionmax) Then
      Write (999, *) ' STOP: IONMAX > 9, REDIMENSION GFFNLTE!'
      Stop ' IONMAX > 9, REDIMENSION GFFNLTE!'
    End If
!   corresponds to NION in nlte_approx

!   prepare ff gaunt-factors, for all possible charges, depth points and
!   freq., &
!   if not already done in opacitm

    If (.Not. optmet) Then
      Do ll = ntemp, nd
        Do kk = 1, ifre
          freq = fre(kk)
          Do iz = izmin, izmax
            z = float(iz)
            gffnlte(iz, ll, kk) = gff(freq*clight, temp(ll), z, ll, kk, iz)*z* &
              z
          End Do
        End Do
      End Do
    End If

  End If

depthloop: Do ll = ntemp, nd

    If (.Not. lteopt .And. optxray) Then
!     JO modified July 2016
      xnel = xnh(ll)*(1.+2.*yhe)
      constx_ll = xnel/clf(ll)*xnh(ll)/clf(ll)*fxl(ll)/(4.*pi)
!     ---- No changes for optically thick clumping, since emission from hot
!     plasma
!     ---- with independent filling factor FXL (either constant or
!     radius-dependent)
    End If

    If (lteopt) Then
      If (xne(ll)/=xnelte(ll)) Then
        Write (999, *) ' STOP: ERROR IN NE -- LTEOPT'
        Stop ' ERROR IN NE -- LTEOPT'
      End If
      xnerat = 1.D0
    Else
!     IF (OPTNEUPDATE.AND.XNE(LL).EQ.XNELTE(LL)) STOP ' ERROR IN NE --
!     NLTEOPT'
      xnerat = xne(ll)/xnelte(ll)
    End If

    Read (17, Rec=ll)(blevel(i), i=1, nrec)
    Read (15, Rec=ll)(alevel(i), i=1, nrec)

!   check for clumping if included (so far, assumed that RAYLEIGH itself not
!   corrected)
!   IF(TEMP(LL).LT.12000.) THEN
!   CALL RAYLEIGH_SCAT(LL)
!   ELSE IF(RAYLEIGH(LL,1).NE.0.D0) THEN
!   RAYLEIGH(LL,:) = 0.D0
!   ENDIF

kloop: Do kk = 1, ifre
      freq = fre(kk)
      If (opteta) hcfre3 = hc2*freq**3
      ehnuekt = exp(-hkl*freq/temp(ll))
      opabf = 0.D0
      etabf = 0.D0

iloop: Do i = 1, index4
!       JO Sept. 2016: for N=20, lower frequencies allowed
        If (iform4(i)/=20) Then
          a = frecin(i)
        Else
          a = frecfn1(i)
!         JO can be skipped later on
          If (a<frecfn(i)) Then
            Write (999, *) ' STOP: FRECFN1 LT FRECFN'
            Stop ' FRECFN1 LT FRECFN'
          End If
        End If
        If (a>freq) Cycle iloop

        If (frenew .And. ll==ntemp) Then ! ONLY ONCE
          Call crossbf(i, alpha, freq)
          If (alpha==0.) Then
            If (op_flag(i) .And. freq<frecin(i)) Then
!             sigma = 0 in OP-data shortwards (Jo: should read longwards???)
!             from FRECFN, no further problems
              Cycle iloop
            Else
              Write (999, *) ' STOP: ERROR IN ALPHA(1)'
              Stop ' ERROR IN ALPHA(1)'
            End If
!           should not happen
          End If
          alphafs(i, kk) = alpha

        Else
          alpha = alphafs(i, kk)
          If (alpha==0.) Then
            Print *, trim(labl(labl4(i))), ll
            Print *, freq, i, op_flag(i), frecin(i), frecfn1(i)
            Write (999, *) ' STOP: ERROR IN ALPHA(2)'
            Stop ' ERROR IN ALPHA(2)'
!           something rotten, since FRECFN1 should be OK now (defined in
!           previous call
!           of OP_RBFSETUP
          End If
        End If

        nn = mlow(i)
        inato = li(nn)
        nato = le(nn)
        xni = blevel(nn)

        If (xni==0.) Then
!         can happen when restart
!         IF (INATO.GE.ILOW(LL,NATO).AND.INATO.LE.IMAX(LL,NATO)) THEN
!         PRINT*,LL,NATO,INATO,ILOW(LL,NATO),IMAX(LL,NATO)
!         STOP ' INCONSISTENCY IN XNI=0 (OPACITC)'
!         ENDIF
          Cycle iloop
        End If

        If (inato>imax(ll,nato)) Cycle iloop
!       JO July 2017: hopefully cures problems in restart with X-rays,
!       when new levels with low occ (but large source-function) are included
        If (inato<ilow(ll,nato)) Then
          If (enionnd(nato,inato,ll)==0.) Cycle iloop
        End If

!       either to ground state (for OP_FLAG and low frequencies)
        If (op_flag(i) .And. freq<frecin(i)) Then
          nm = kl(nn)
!         IF(LL.EQ.1) PRINT*,LABL(LABL4(I)),FREQ,FRECIN(I),NM
!         or to ground/excited state as defined in DETAIL input
        Else
          nm = mup(i)
        End If

!       --- outside ilow, imax!


        ij = igenio(nato, imax(ll,nato)) + 1
        nlas = ifirsl(ij)

!       --- ionization to excited level of imax+1, consistent with netmat
!       completely changed, see notes.
!       idea: assume excited levels to be in LTE to ground-state of nk = nk1
!       then: nistar = nki * (ni/nki)* = nk1 * (ni/nk1)*

!       Note also that there was a bug here in former versions:
!       We used an alternative formulation, but only if xnk_i = 0,
!       although xnk_i was at its (absolute) LTE value. Since the ground-state
!       was in NLTE, the ratio was wrong, and nistar was wrong as well.


        If (nm>nlas) Then !         this is the condition
          If (op_flag(i) .And. freq<frecin(i)) Stop &
            ' NM > NLAS for low freq. OP-data!'
          nm = nlas
        End If
!       THAT'S ALL!!!
        xnk = blevel(nm)
        xnistar = xnk*alevel(nn)/alevel(nm)*xnerat
        If (xnistar==0.) Then
          Print *, ll, nn, nm, xnk, alevel(nn), alevel(nm), xnerat
          Print *, nato, ilow(ll, nato), imax(ll, nato)
          Print *, blevel
          Write (999, *) ' STOP: XNISTAR = 0!'
          Stop ' XNISTAR = 0!'
        End If
        depart = xni/xnistar

        If (lteopt .And. abs(depart-1.D0)>precis) Then
          Write (999, *) ' STOP: ERROR IN DEPART - LTE'
          Stop ' ERROR IN DEPART - LTE'
        End If

        opabf = opabf + xnistar*alpha*(depart-ehnuekt)
!       if(kk.ge.739.and.kk.le.743.and.mod(ll,6).eq.0 .and. &
!       &                nato.eq.3.and.inato.ge.2) then
!       print*,ll,kk,XNISTAR*ALPHA*
!       (DEPART-EHNUEKT),i,nato,inato,nm,nn,alpha,xni,xnk,xnistar,depart
!       endif

        If (opteta) etabf = etabf + xnistar*alpha

      End Do iloop

!     -------- now we have calculated bf opacity for all transitions
!     ------             we start ff opacity

      opaff = 0.D0
      chiff = 1.3695D-23*xne(ll)/(freq**3*sqrt(temp(ll)))
      chiffc = chiff*(1.D0-ehnuekt)

      If (opteta) etaffc = chiff

      Do i = 1, index5
        ii = ipare5(i)
        nato = le(ii)
        inato = li(ii)
        en = enionnd(nato, inato, ll)
        If (en==0.) Then
!         can happen when restart
!         IF (INATO.GE.ILOW(LL,NATO).AND.INATO.LE.IMAX(LL,NATO)+1) THEN
!         PRINT*,LL,NATO,INATO,ILOW(LL,NATO),IMAX(LL,NATO)
!         STOP ' INCONSISTENCY IN ENIONND = 0 (OPACITC)'
!         ENDIF
          Cycle
        End If
        z = zl(ii)
        iz = int(z+1.D-6)
        gf = gffnlte(iz, ll, kk)
!       FOR TESTS
!       IF(GF.LE.0.) STOP 'ERROR IN FF GAUNT-FACTORS'
        opaff = opaff + gf*en
      End Do

      opaff1 = opaff*chiffc
!     the following two quanties remain uncorrected for clumping, since
!     they are needed only for the T-structure and all other quantities
!     entering the energy balance remain uncorrected as well

      ffopa(ll, kk) = opaff*chiff ! maup
      dtffopa(ll, kk) = (-1.D0)*ffopa(ll, kk)/(2.D0*temp(ll))
!     considering that the gaunt fact. does not depend on temperature

      opath = xne(ll)*sig/clf(ll) ! xne includes cl =>  mean opacity

      If (temp(ll)<=10000. .And. 1.D8/freq>200.D0) Then
!       opahminus yields uncorrected values
        Call hminus(xne(ll), temp(ll), 1.D8/freq, ehnuekt, ll, opahminus)
      Else
        opahminus = 0.D0
      End If
      opab = opabf + opaff1 + opahminus
!     still uncorrected

!     this is the only place where the k-shell opacities should appear
      If (.Not. lteopt .And. optxray .And. kk>=kmin) opab = opab + &
        opa_kshell(ll, kk)

      opab = opab/clf(ll) !         now corrected => mean opacity

!     OPAC(LL,KK) = OPAB + OPATH + RAYLEIGH(LL,KK)/CLF(LL)
!     print*,ll,1.d8/freq,opab,rayleigh(ll,kk)

      opac(ll, kk) = opab + opath ! corrected => mean opacity

!     ---- Now, OPAC is changed below to its EFFECTIVE value
      If (optmet) Then
!       opat_m_new etc corrected (mean opacities)

        opat_m_old(ll, kk) = opat_m_new(ll, kk) ! already corrected in OPACITM
        strue_m_old(ll, kk) = strue_m_new(ll, kk)
        opac(ll, kk) = opac(ll, kk) + opat_m_new(ll, kk) ! so far, mean

!       if(mod(ll,6).eq.0.and..not.lteopt) &
!       print*,'old',kk,ll,1.d8/fre(kk), &
!       &
!       (opat_m_new(ll,kk)-opat_m_nolines(ll,kk))/opac_nolines(ll,kk)
        opac_nolines(ll, kk) = opab + opath + opat_m_nolines(ll, kk)
!       ---- Remains as mean (not effective) for next iteration
!       NOTE: K-Shell absorption included in both (OPAC and OPAC_NO_LINES)

!       if(mod(ll,6).eq.0.and..not.lteopt) &
!       print*,'new',kk,ll,1.d8/fre(kk), &
!       &              (opac(ll,kk)/opac_nolines(ll,kk)-1.)

!       ---- Inside or outside line-blocking range?
        lam = 1.D8/fre(kk)
        If (lam>w_bg_red .Or. lam<w_bg_blue) Then
!         JO test for transition effects (are all frequencies treated?)
!         ---- Just continuum bg opacity
          tcl = opac_nolines(ll, kk)*tcl_fac_cont(ll)
          opa_eff_rat(ll, kk) = (1.+fic(ll)*tcl)/(1.+tcl)
        End If

      Else
!       ---- for models with only explicit elements (e.g., pure HHe),
!       ---- OPAC consists of continuum only. THUS
        tcl = opac(ll, kk)*tcl_fac_cont(ll)
        opa_eff_rat(ll, kk) = (1.+fic(ll)*tcl)/(1.+tcl)
      End If

      If (.Not. optthick) Then
        If (abs(opa_eff_rat(ll,kk)-1.D0)>epsi) Then
!         ---- this is a good test, since also in thin clumping opa_eff_grid
!         is
!         ---- calculated in sumopal. Test proves appropriate averaging
          Print *, 'opa', ll, kk, 1.D8/freq, opa_eff_rat(ll, kk) - 1.D0
          Write (999, *) &
            ' STOP: OPTICALLY THIN CLUMPING AND OPA_EFF_RAT NE 1 IN OPACITC'
          Stop ' OPTICALLY THIN CLUMPING AND OPA_EFF_RAT NE 1 IN OPACITC'
        End If
!       FINALLY, SET TO UNITY
        opa_eff_rat(ll, kk) = 1.D0
      End If

      opa_eff_rat_old(ll, kk) = opa_eff_rat(ll, kk)

      If (opa_eff_rat(ll,kk)>epsi1 .Or. opa_eff_rat(ll,kk)<=0.D0) Then
        Print *, opa_eff_rat(ll, kk), ll, kk
        Write (999, *) ' STOP: OPA_EFF_RAT OUT OF BOUNDS IN OPACITC_1'
        Stop ' OPA_EFF_RAT OUT OF BOUNDS IN OPACITC_1'
      End If
!     ---- This is the crucial correction!
      opac(ll, kk) = opac(ll, kk)*opa_eff_rat(ll, kk)


      If (opteta) Then

        If (ehnuekt/=0.D0) Then
          etab = ehnuekt*hcfre3*(etabf+opaff*etaffc) ! still uncorrected
        Else
          x = hkl*freq/temp(ll)
          If (x<200.D0) Then
            Write (999, *) ' STOP: OPACITC: PROBLEM AT HIGH ENERGIES'
            Stop ' OPACITC: PROBLEM AT HIGH ENERGIES'
          End If
          expo = log(hcfre3*(etabf+opaff*etaffc)) - x
          etab = exp(expo) !        still uncorrected
        End If

        If (temp(ll)<=10000. .And. 1.D8/freq>200.D0) etab = etab + &
          opahminus*hcfre3/(1.D0/ehnuekt-1.D0) ! still uncorrected

        etab = etab/clf(ll) !       now corrected => mean emissivity
!       JO: here we have to think about correction (mean or effective
!       emissivity???)
        etat_nolines(ll, kk) = etab + etat_m_nolines(ll, kk)

        etab = etab*opa_eff_rat(ll, kk) ! now effective emissivity

!       add x-ray emissivity
        If (.Not. lteopt .And. optxray .And. ll<=lxmin) Then
          x = constx_ll*lambdanu(kk, ll)
          etab = etab + x !         2nd term corrected, since constx_ll
!         includes fx
!         JO May 2025: following statement forgotten in previous versions
          etat_nolines(ll, kk) = etat_nolines(ll, kk) + x
        End If
!       for tests
!       IF(.NOT.LTEOPT.AND.OPTXRAY .AND. LL.LE.LXMIN) THEN
!       IF(KK.GE.KMIN) THEN
!       WRITE(*,5) LL,KK,OPA_KSHELL(LL,KK),CONSTX_LL,LAMBDANU(KK,LL)
!       ELSE
!       WRITE(*,5) LL,KK,ZERO,CONSTX_LL,LAMBDANU(KK,LL)
!       ENDIF
!       ENDIF
!       5 FORMAT('XRCF',I3,2X,I5,2X,3(E10.4,2X))

!       ---- we do NOT correct X-ray emission for optically thick effects
!       (assuming X-ray emission/absorption dispatched, see notes etc.)

        If (.Not. optmet) Then
          strue(ll, kk) = etab/opac(ll, kk) ! correction cancels out anyway
        Else
!         STRUE(LL,KK) = &
!         &                (ETAB + OPAT_M_NEW(LL,KK)*STRUE_M_NEW(LL,KK)) /
!         OPAC(LL,KK)
!         ---- needed to make opat_m_new effective quantity here!
          strue(ll, kk) = (etab+opa_eff_rat(ll,kk)*opat_m_new(ll,kk)* &
            strue_m_new(ll,kk))/opac(ll, kk)
        End If

!       IF(LL.EQ.1.AND.OPTXRAY.AND..NOT.LTEOPT) THEN
!       PRINT*
!       IF(KK.GE.KMIN) THEN
!       PRINT*,1.D8/FRE(KK),OPAC(LL,KK),OPA_KSHELL(LL,KK)
!       ELSE
!       PRINT*,1.D8/FRE(KK),OPAC(LL,KK),'0.'
!       ENDIF
!       PRINT*,ETAB,CONSTX_LL*LAMBDANU(KK,LL),STRUE(LL,KK)
!       PRINT*,OPAT_M_NEW(KK,LL),STRUE_M_NEW(LL,KK)
!       ENDIF


!       print*,ll,kk,1.d8/fre(kk),strue(ll,kk)
        If (strue(ll,kk)==0.) Then
          If (1.D8/freq>20.D0) Then
            Write (999, *) ' STOP: LAMBDA > 20 A and STRUE = 0!'
            Stop ' LAMBDA > 20 A and STRUE = 0!'
          End If
        End If

        If (lteopt) Then
!         JO MODIFIED NOV. 2025
          x = hkl*freq/temp(ll)
!         IF(EHNUEKT.NE.0.D0) THEN
          If (x<200.D0) Then !      see also subr. opacitm
            bn = hcfre3/(1.D0/ehnuekt-1.D0)
          Else
            expo = log(hcfre3) - x
            bn = exp(expo)
          End If
!         ---- since ETAB has been corrected, but not OPAB we need to adjust
!         for that here
!         BN1 = ETAB/OPAB
          bn1 = etab/(opab*opa_eff_rat(ll,kk))
          If (1.D8/freq>20.D0) Then
            If (abs(1.D0-bn1/bn)>precis) Then
              Print *, 'PROBLEMS IN LTE:', bn, bn1 ! !!
              If (bn>1.D-100) Then
                Write (999, *) ' STOP: ERROR IN LTE!!!'
                Stop ' ERROR IN LTE!!!'
              End If
            End If
          End If
        End If
      End If

    End Do kloop

  End Do depthloop

! for tests
! DO KK = 1,IFRE
! FREQ = 1.d8/FRE(KK)
! print*,freq,opac(54,kk),opac_nolines(54,kk)
! enddo


! finally, all opacities (and source-functions) have been corrected w.r.t. CLF
! i.e., are mean quantities, in particular OPAC_NOLINES,
! except for FFOPA (which remains uncorrected at all),
! and ETAB and OPAC(L,K), which have been corrected for macro-clumping.

! SUMMARY: on output, OPAC is effective, whilst all other quantities remain
! as in the optically thin case.
  Return

End Subroutine

!-----------------------------------------------------------------------

Subroutine crossbf(i, alpha, freq)
! JO Sept 2021 new provisional treatment for formula 12
  Use :: nlte_type
  Use :: nlte_dim
  Use :: fund_const

  Use :: princesa_var, Only: le, li, zl, weight, data, iform4, indat4, labat, &
    labl, labl4, frecin, frecfn

  Use :: nlte_var, Only: op_flag, op_data, ifre, fre, modnam, frecfn1

  Implicit None

! ------ this subroutine calculates photoionization cross-sections in the
! ------ way described by 'detail'. the cross-section is return in
! ------ 'alpha'

! ------ used (and checked so far)

! ------ formula 4: mg i, ii
! ------ formula 5: h i, he i, he ii, si ii,iii,iv
! ------ formula 6: he i
! ------ formula 7: he i
! ------ formula 2: he i
! ------ formula 3: he i
! ------ formala 12: si ii,iii,iv. attention: assuming here that only
! ------             first three coefficients are different from zero!!!
! ------             o ii, o iii: all coefficients might be different from
! zero
! ------             provisional treatment for lambda < 100 A
! ------ formula 13: o ii (rbf)
! ------ formula 16: he i (rbf)
! ------ formula 17: he i (rbf)
! ------ formula 20: OPACITY PROJECT data, convolved with Gaussian of average
! width E/Delta E = 100
! ------ formula 101: rbf data (new detail format), convolved with Gaussian
! of average width E/Delta E = 100

! ------ formula 51: Seaton approximation
! ------ formula 21: Seaton approximation, followed by DR-data for the
! stabilizing transitions

! hay que definir tambien las matrices de /atomdat/


! .. scalar arguments ..

  Integer (i4b), Parameter :: mmm = 1
  Real (dp), Dimension (mmm) :: alpham, xnum

  Real (dp) :: alpha, freq
  Integer (i4b) :: i
! ..
! .. local scalars ..
  Real (dp) :: alp0, alp1, alp2, g, summ, x, xip, xiq, xlog, xnu, z, z2, aw, &
    alp3, alp4, alp5, x100, alpha100
  Integer (i4b) :: ii1, in, n, natomo, ndatos, nivll, numi, ll, nn, il, is, iz
  Character :: key*6

! for the use of the OP data
  Integer (i4b) :: ndatosmax, nrecord, k
  Real (dp) :: peq = 0.D0, lop1 = 0.D0, lop2 = 0.D0 ! hawaii

  Character (30) :: namefileop
  Character (80) :: dummychar
! ..
! ..
! .. external functions ..
  Real (dp) :: agau, gaunt
  External :: agau, gaunt
! ..
! .. intrinsic functions ..
! INTRINSIC EXP,INT,LOG
! ..
! ..

  n = iform4(i)

! JO Sept. 2016: for N=20, frequencies down to ground-state ionization allowed
! NOTE: at first call, FRECFN1 is equal to FRECFN. After later calls
! (op_flag=.true.) FRECFN1 (with FRECFN LE FRECFN1 LE FRECIN) present

  If (n/=20) Then
    x = frecin(i)/freq
  Else
    x = frecfn1(i)/freq
  End If

  If (x>1.D0) Then
    Write (999, *) ' STOP: ERROR IN THRESHOLD '
    Stop ' ERROR IN THRESHOLD '
  End If

  xnu = freq*clight
  ii1 = indat4(i)

  If (n==1) Then
    alp0 = data(ii1)
    alpha = alp0*x**3
    Return

  Else If (n==2) Then
    alp0 = data(ii1)
    xip = data(ii1+1)
    alpha = alp0*x**(-xip)
    Return

  Else If (n==3) Then
    alp0 = data(ii1)
    xip = data(ii1+1)
    alp1 = data(ii1+2)
    xiq = data(ii1+3)
    alpha = alp0*x**(-xip) + alp1*x**(-xiq)
    Return

  Else If (n==4) Then
    alp0 = data(ii1)
    alp1 = data(ii1+1)
    alp2 = data(ii1+2)
    alpha = alp0*(alp1+(1.D0-alp1)*x)*x**alp2
    Return

  Else If (n==5) Then
    in = int(data(ii1)+0.5D0)
    nivll = labl4(i)
    z = zl(nivll) + 1.D0
    z2 = z*z
    natomo = le(nivll)
    key = labat(natomo)
    numi = li(nivll)
    If (key=='H' .Or. key=='HE') Then

!     artemios gaunt faktoren

      g = agau(in, xnu, key, numi)
    Else

!     mihalas gaunt faktoren

      g = gaunt(in, xnu/z2)
    End If
    alpha = 2.815D29*z2*z2/((in**5)*(xnu**3))*g
    Return

  Else If (n==6) Then
    alp0 = data(ii1)
    alp1 = data(ii1+1)
    alp2 = data(ii1+2)
    alpha = alp0*exp(alp1+alp2*xnu)
    Return

  Else If (n==7) Then
    alp0 = data(ii1)
    alp1 = data(ii1+1)
    alp2 = data(ii1+2)
    alpha = exp(alp0+log(xnu)*(alp1+alp2*log(xnu)))
    Return

  Else If (n==12) Then
!   JO new provisional version from Sept. 2021
    alp0 = data(ii1)
    alp1 = data(ii1+1)
    alp2 = data(ii1+2)
    alp3 = data(ii1+3)
    alp4 = data(ii1+4)
    alp5 = data(ii1+5)

    xlog = log(x)
    If (alp3==0. .And. alp4==0. .And. alp5==0.) Then
      summ = alp0 + xlog*(alp1+xlog*alp2)
    Else
      summ = alp0 + xlog*(alp1+xlog*(alp2+xlog*(alp3+xlog*(alp4+xlog*alp5))))
    End If
    alpha = exp(summ)

    If (freq<=1.D6) Return !        100 A

    x100 = frecin(i)/1.D6
    xlog = log(x100)
    If (alp3==0. .And. alp4==0. .And. alp5==0.) Then
      summ = alp0 + xlog*(alp1+xlog*alp2)
    Else
      summ = alp0 + xlog*(alp1+xlog*(alp2+xlog*(alp3+xlog*(alp4+xlog*alp5))))
    End If
    alpha100 = exp(summ)
    If (alpha>alpha100) Then !      PRESUMABLY OUT OF VALIDITY RANGE
!     PRINT*,' INCREASING ALPHA: ',FREQ,I,ALPHA,ALPHA100
      alpha = alpha100*(1.D6/freq)**2
!     PRINT*,' new approximate value: ',ALPHA
    End If
    Return

  Else If (n==13) Then !            not checked by J.Puls
    alp0 = data(ii1)
    alp1 = data(ii1+1)
    alp2 = data(ii1+2)
    alp3 = data(ii1+3)

    summ = alp0*((alp1+(1.D0-alp1)*x)*x**alp2*(1.D0-exp(alp3*(1.D0/x-1.D0)))+ &
      exp(1.D0-1.D0/x))
    alpha = summ
    Return

  Else If (n==16) Then
    nivll = labl4(i)
    natomo = le(nivll)
    aw = weight(natomo)
    nn = int(data(ii1)+0.5)
    ll = int(data(ii1+1)+0.5)
    z = zl(nivll) + 1.D0
    xnum(1) = xnu
    Call pix11(z, aw, nn, ll, mmm, xnum, alpham)
    alpha = alpham(1)
    Return

  Else If (n==17) Then
    nn = int(data(ii1)+0.5)
    is = int(data(ii1+1)+0.5)
    il = int(data(ii1+2)+0.5)
    Call pix21(is, il, nn, xnu, alpha)
    Return

  Else If (n==20) Then !            OP data

!   JO changed Sept. 2016

    If (.Not. op_flag(i)) Then !    this is only required once
      alp0 = int(data(ii1)+0.5D0)
      ndatosmax = int(data(ii1+1)+0.5D0) ! maximum number of data points
      nrecord = int(data(ii1+2)+0.5D0) ! record to be read
      ndatos = int(data(ii1+3)+0.5D0) ! number of OP data points in this
!     record
      Call op_rbfsetup(i, alp0, ndatosmax, nrecord, ndatos, namefileop)
      op_flag(i) = .True.
      If (freq<frecfn1(i)) Then
!       in first call(s) of CROSSBF (OP_FLAG=.FALSE.), FREQ can be located
!       between
!       FRECFN and FRECFN1, since FRECFN1 updated only in OP_RBFSETUP
        If (freq<frecfn(i)) Then
          Write (999, *) ' STOP: ERROR IN FREQ (OP_RBFSETUP)'
          Stop ' ERROR IN FREQ (OP_RBFSETUP)'
        End If
        alpha = 0.
        Return
      End If
    End If

!   After the first time that is called, the smoothed cross sections
!   are calculated by interpolation of the OP data stored in OP_DATA matrix, &
!   one point per each of the frequencies in the (coarse) grid.

    Do k = 2, ifre
      If (freq<=fre(k)) Go To 100
    End Do
    Write (999, *) ' STOP: freq not found in fre(k) -- crossbf'
    Stop ' freq not found in fre(k) -- crossbf'

!   JO CHANGED SEPT 2016
!   10   PEQ = FREQ - FRECIN(I)
100 peq = freq - frecfn1(i)

!   fre(k) > freq > fre(k-1)

    If (abs(peq)<=1.D-5 .Or. (fre(k)-freq)<=1.D-5) Then
      alpha = op_data(i, k) !       exactly at the edge OR freq. point
    Else
      lop1 = log10(freq/fre(k))/log10(fre(k-1)/fre(k))
      lop2 = 1.D0 - lop1
      alpha = lop1*log10(op_data(i,k-1)) + lop2*log10(op_data(i,k))
      alpha = 10.D0**alpha
    End If

    If (alpha<=0.D0) Then
      nrecord = int(data(ii1+2)+0.5D0) ! record to be read
      Print *
      Print *, ' ERROR in OP data file ', trim(namefileop), &
        ' cross-section ZERO OR NEGATIVE!!!!!', alpha
      Print *, ' Transition numer: ', i, ' at frequency: ', freq
      Print *, ' OP data record: ', nrecord, ' corresponds to level: ', &
        labl(labl4(i))
      Print *, ' level''s ionization edge : ', frecfn(i)
      Print *, ' frequencies : ', fre(k), ' ', fre(k+1), ' ', op_data(i, k), &
        ' ', op_data(i, k+1)
      Print *
      dummychar = trim(modnam) // trim('/crossbf-error-'//trim(labl(labl4(i))) &
        //'.dat')
      Open (77, File=dummychar, Status='UNKNOWN')
      Do k = 1, ifre
        Write (77, *) fre(k), op_data(i, k)
      End Do
      Close (77)
      Write (999, *) ' STOP!'
      Stop
    End If
    Return

  Else If (n==21) Then
!   first three data points for Seaton approx
    alp0 = data(ii1)
    alp1 = data(ii1+1)
    alp2 = data(ii1+2)
    alpha = alp0*(alp1*x**alp2+(1.D0-alp1)*x**(alp2+1.D0))
    If (alpha<0.) Then
      If (alp1>=0.) Then
        Write (999, *) ' STOP: ALPHA NEGATIVE AND BETA POSITIV'
        Stop ' ALPHA NEGATIVE AND BETA POSITIV'
      End If
      alpha = 1.D-40
    End If
    Return

  Else If (n==33) Then
    Print *, ' FORMULA NO. ', n, ' NOT CHECKED YET'
    Write (999, *) ' STOP!'
    Stop
!   alpha0=data(ii1)*1.d-18
!   beta=data(ii1+1)
!   ese=data(ii1+2)
!   nivel=data(ii1+3)
!   nivll=labl4(i)
!   natomo=le(nivll)
!   key=labat(natomo)
!   numi=li(nivll)
!   c          facgau=agau(nivel,xnu,key,numi)
!   c          alpha0=alpha0*facgau
!   alpha=alpha0*(x**ese)*(beta+(1.d0-beta)*x)
!   return

  Else If (n==51) Then
    alp0 = data(ii1)
    alp1 = data(ii1+1)
    alp2 = data(ii1+2)
    alpha = alp0*(alp1*x**alp2+(1.D0-alp1)*x**(alp2+1.D0))
    If (alpha<0.) Then
      If (alp1>=0.) Then
        Write (999, *) ' STOP: ALPHA NEGATIVE AND BETA POSITIV'
        Stop ' ALPHA NEGATIVE AND BETA POSITIV'
      End If
      alpha = 1.D-40
    End If
    Return

!   comment JO Sept 2021: modify according to formula 20
  Else If (n==101) Then !           rbf cross-sections in new DETAIL format
    Print *, &
      ' Cross sections according to FORMULA 101 need to be adapted/tested'
    Print *, ' whenever ionization to excited level occurs'
    Print *, ' (effective edge is ground-state edge due to dielectr. recomb.)'
    Print *, ' Proceed in analogy to subr. OP_RBFSETUP'
    Write (999, *) ' STOP: Cross sections according to FORMULA &
      &101 need to be adapted/tested'
    Stop ' Cross sections according to FORMULA 101 need to be adapted/tested'
    nivll = labl4(i)
    natomo = le(nivll)
    key = labat(natomo)
    z = zl(nivll) + 1.D0
    iz = int(z+0.5D0)
    If (z/=int(data(ii1))) Then
      Write (999, *) ' STOP: ERROR IN RBF CODING (FORMULA101)'
      Stop ' ERROR IN RBF CODING (FORMULA101)'
    End If
    namefileop = trim(key) // achar(48+iz) // 'PIC' ! assumes specific
!   labelling of data-files
    nrecord = int(data(ii1+1)) !    corresponds to label
    If (.Not. op_flag(i)) Then !    this is only required once
      Call op_rbfsetup1(i, key, iz, nrecord, namefileop)
      op_flag(i) = .True.
    End If

!   After the first time that is called, the smoothed cross sections
!   are calculated by interpolation of the OP data stored in OP_DATA matrix, &
!   one point per each of the frequencies in the (coarse) grid.

    Do k = 2, ifre
      If (freq<=fre(k)) Go To 110
    End Do
    Write (999, *) ' STOP: freq not found in fre(k) -- crossbf'
    Stop ' freq not found in fre(k) -- crossbf'

110 peq = freq - frecin(i)

!   fre(k) > freq > fre(k-1)

    If (abs(peq)<=1.D-5 .Or. (fre(k)-freq)<=1.D-5) Then
      alpha = op_data(i, k) !       exactly at the edge OR freq. point
    Else
      lop1 = log10(freq/fre(k))/log10(fre(k-1)/fre(k))
      lop2 = 1.D0 - lop1
      alpha = lop1*log10(op_data(i,k-1)) + lop2*log10(op_data(i,k))
      alpha = 10.D0**alpha
    End If

    If (alpha<=0.D0) Then
      Print *
      Print *, ' ERROR in OP data file ', trim(namefileop), &
        ' cross-section ZERO OR NEGATIVE!!!!!', alpha
      Print *, ' Transition numer: ', i, ' at frequency: ', frecin(i)
      Print *, ' OP data record: ', nrecord, ' corresponds to level: ', &
        labl(labl4(i))
      Print *, ' level''s ionization edge : ', frecin(i)
      Print *, ' frequencies : ', fre(k), ' ', fre(k+1), ' ', op_data(i, k), &
        ' ', op_data(i, k+1)
      Print *
      dummychar = trim(modnam) // trim('/crossbf-error-'//trim(labl(labl4(i))) &
        //'.dat')
      Open (77, File=dummychar, Status='UNKNOWN')
      Do k = 1, ifre
        Write (77, *) fre(k), op_data(i, k)
      End Do
      Close (77)
      Write (999, *) ' STOP!'
      Stop
    End If
    Return


  Else
!   ihi=labl4(i)
!   call alphas(ihi,alpha,freq)
!   return
    Print *, ' FORMULA NO. ', n, ' NOT IMPLEMENTED YET'
    Write (999, *) ' STOP!'
    Stop
  End If

End Subroutine

!-----------------------------------------------------------------------

Subroutine op_codex(num, name, checkname, rydberg)

! It translates the codification of the file name for the OP photionization
! cross section data

  Use :: nlte_type
  Use :: nlte_dim
  Implicit None

! .. scalar arguments ..
  Integer (i4b) :: num
  Character (30) :: name, checkname
  Real (dp) :: rydberg
! ..
! .. local scalars ..
  Integer (i4b) :: i, j, div
  Integer (i4b) :: limit, indice
! integer(I4B) :: num0
! ..
  Real (dp), Dimension (26) :: auxrydberg
  Character (6), Dimension (26) :: atomo
  Character (5), Dimension (10) :: estado
  Character (8), Dimension (3) :: xsec_data_name, xsec_control_name
  Character (20) :: auxname

  Data atomo/'H', 'He', 'X', 'X', 'X', 'C', 'N', 'O', 'F', 'Ne', 'Na', 'Mg', &
    'Al', 'Si', 'P', 'S', 'Cl', 'Ar', 'K', 'Ca', 'Sc', 'Ti', 'V', 'Cr', 'Mn', &
    'Fe'/
  Data estado/'I', 'II', 'III', 'IV', 'V', 'VI', 'VII', 'VIII', 'IX', 'X'/
  Data auxrydberg/109677.62, 109722.28, 109728.64, 109730.64, 109731.75, &
    109732.30, 109733.02, 109733.55, 109734.15, 109734.33, 109734.70, &
    109734.84, 109735.08, 109735.17, 109735.37, 109735.44, 109735.62, &
    109735.81, 109735.78, 109735.81, 109735.98, 109736.06, 109736.13, &
    109736.16, 109736.22, 109736.24/
  Data xsec_data_name/'_rbf_dat', '_xop_dat', '_nah_dat'/
  Data xsec_control_name/'_rbf_lis', '_xop_lis', '_nah_lis'/


! DIV=1000
! IF (NUM .LT. 1000) DIV=10

! num0 = num

  limit = 1000
  If (num<limit) Then
    div = 10
    indice = 1
  Else
    indice = int(num/limit+0.5)
    num = num - limit*indice
    div = 10
    indice = indice + 1
  End If

  i = int(num/div+0.5)
  j = int(num-i*div)

  If (j>i) Then
    Print *, 'I: ', i, '  J: ', j
    Write (999, *) ' STOP: ERROR IN OP DATA CODEX'
    Stop ' ERROR IN OP DATA CODEX'
  End If
! AUXNAME=TRIM(ATOMO(I))//TRIM(ESTADO(J))//'_rbf_dat'
  auxname = trim(atomo(i)) // trim(estado(j)) // trim(xsec_data_name(indice))
  name = trim(auxname)
! CHECKNAME=TRIM(ATOMO(I))//TRIM(ESTADO(J))//'_rbf_lis'
  checkname = trim(atomo(i)) // trim(estado(j)) // trim(xsec_control_name( &
    indice))

! for debugging purposes only
! print*,NUM,LIMIT,INDICE
! print*,auxname
! print*,checkname

  rydberg = auxrydberg(i)
  Return
End Subroutine

!-----------------------------------------------------------------------

Subroutine prmodel(r, v, gradv, rho, xne, xnh, temp, tau, clfac, teff, n)

  Use :: nlte_type
  Use :: nlte_dim
  Implicit None

! druckt r in sternradien,v in einheiten von v-unendlich,
! log rho, log ne, log nh, t-lucy/teff  und tau-thomson aus


! .. parameters ..
  Integer (i4b), Parameter :: nd1 = id_ndept
! ..
! .. scalar arguments ..
  Real (dp) :: teff
  Integer (i4b) :: n
! ..
! .. array arguments ..
  Real (dp) :: gradv(nd1), r(nd1), rho(nd1), tau(nd1), temp(nd1), v(nd1), &
    xne(nd1), xnh(nd1), clfac(nd1)
! ..
! .. local scalars ..
  Integer (i4b) :: i
! ..
! .. intrinsic functions ..
! INTRINSIC LOG10
! ..

  Write (*, Fmt=100)

  Do i = 1, n
    Write (*, Fmt=110) i, r(i), v(i), gradv(i), log10(rho(i)), log10(xnh(i)), &
      log10(xne(i)), temp(i)/teff, tau(i), clfac(i)
  End Do

! do i=1,n
! if(gradv(i).le.0.d0) then
! print*,' negative vel. gradient -- density inversion! at n =',i
! print*,' increase log g!!!!'
! stop ' density inversion'
! endif
! END DO


  Return

100 Format (' NO   RADIUS     VELOCITY   GRADIENT  DENSITY(AV) &
    &N-HYD(ENH) N-EL(ENH)   T/TEFF TAU-TH(MIC) ENH-FAC', /, /)
110 Format (I3, 2X, 8(F9.5,2X), F6.2)

End Subroutine

!-----------------------------------------------------------------------

Subroutine tlucy(r, taul, tau, t, n, tmin, teff, rtau23, ithyd)

  Use :: nlte_type
  Use :: nlte_dim
! & IDUM,GGRAV,CKAPPA,XT,A2TE,DELSIG,DELMU,DQDTAU, &
! JO: fixed July 2015 (XT not needed here, overlap with XT from photstruc)
  Use :: nlte_var, Only: opt => optlucy, optmodel, srvmin, srnom, constt, &
    nsdiv, dqdtau, qinf, q0, gam => gamhopf
  Use :: photstruc

  Implicit None

! berechnet die temperatur nach dem verfahren von
! lucy bei vorgegebener dichte. (opt=.true.)

! constt already with respect to rtau23 = srnom/sr

! uses nlte hopf function (opt=.false.)

! for optmodel=.true., input T-struct. is used

! .. parameters ..
  Integer (i4b), Parameter :: nd1 = id_ndept
! ..
! .. scalar arguments ..
  Real (dp) :: rtau23, teff, tmin
  Integer (i4b) :: ithyd, n
! ..
! .. array arguments ..
  Real (dp) :: r(nd1), t(nd1), tau(nd1), taul(nd1)
! ..
! .. local scalars ..
  Real (dp) :: dtaudr, q, qh, qold, qprime, r23, ri, rold, rp, tauold, taup, &
    taur, w, w1, x, dtdm, tminkur
  Integer (i4b) :: k, i
  Logical :: optfound
! ..
! .. intrinsic functions ..
! INTRINSIC EXP,SQRT
! ..

  optfound = .False.

  If (opt) Then

    tauold = 0.D0
    If (.Not. optmodel .Or. ithyd==1) Then

      w = .5D0*(1.D0-sqrt(1.D0-(rtau23/r(1))**2))
      t(1) = teff*(w+constt)**.25D0

      Do k = 1, n - 1
        If (t(k)<tmin) t(k) = tmin
        w1 = 1.D0 - (rtau23/r(k+1))**2
        If (w1<0.D0) w1 = 0.D0
        w = .5D0*(1.D0-sqrt(w1))
        taup = taul(k+1)

        If (taup>=.66D0 .And. tauold<.66D0) Then
          optfound = .True.
          rold = r(k)
          rp = r(k+1)
          dtaudr = (taup-tauold)/(rp-rold)
          r23 = rold + (.66D0-tauold)/dtaudr
          If (r23>rold .Or. r23<rp) Then
            Write (999, *) ' STOP: ERROR IN R(TAU-LUCY=2/3)'
            Stop ' ERROR IN R(TAU-LUCY=2/3)'
          End If
          srvmin = srnom*r(nsdiv)/r23
          Print *
          Print *, ' MESSAGE FROM TLUCY: SRVMIN/SRNOM = ', srvmin/srnom
          Print *
        End If
        tauold = taup
        t(k+1) = teff*(w+.75D0*taup+constt)**.25D0
      End Do

    Else
!     interpolation of input structure with respect to m
      t(1) = tmin
      tminkur = 0.

      Do k = 1, n - 1

        taup = taul(k+1)

        If (taup>=.66D0 .And. tauold<.66D0) Then
          optfound = .True.
          rold = r(k)
          rp = r(k+1)
          dtaudr = (taup-tauold)/(rp-rold)
          r23 = rold + (.66D0-tauold)/dtaudr
          If (r23>rold .Or. r23<rp) Then
            Write (999, *) ' STOP: ERROR IN R(TAU-LUCY)=2/3)'
            Stop ' ERROR IN R(TAU-LUCY)=2/3)'
          End If
          srvmin = srnom*r(nsdiv)/r23
          Print *
          Print *, ' MESSAGE FROM TLUCY: SRVMIN/SRNOM = ', srvmin/srnom
          Print *
        End If
        tauold = taup

        x = xm(k+1)

        If (x==0.) Then
          t(k+1) = tmin
        Else

          Do i = 1, ndmod - 1
            If (x>xmext(i) .And. x<xmext(i+1)) Exit
          End Do
          If (i==ndmod) Then
            If (x<xmext(1)) Then
              i = 1
            Else
              i = ndmod - 1
!             STOP ' xm not found in external phot. struct. (tlucy)'
            End If
          End If

          dtdm = log10(text(i+1)/text(i))/log10(xmext(i+1)/xmext(i))
          t(k+1) = log10(text(i)) + dtdm*log10(x/xmext(i))
          t(k+1) = 10.**t(k+1)
!         IF (TMINKUR.EQ.0.) TMINKUR=T(K+1)
          If (tminkur==0.) tminkur = max(tmin, t(k+1))
        End If
      End Do

      Where (t==tmin) t = tminkur
      tmin = tminkur
!     until here, we have only the othermost points set to
!     max(tmin_input,tmin(kur)
      Where (t<tmin) t = tmin
!     the latter statement ensures to have all temperatures above tmin7tminkur
!     in case the t-struct. is not montonic, this feature will be removed
!     comment out the line if you want to preserve this feature
    End If

    If (ithyd>=2 .And. .Not. optfound) Then
      Write (999, *) ' STOP: TAU-LUCY = 2/3 NOT FOUND'
      Stop ' TAU-LUCY = 2/3 NOT FOUND'
    End If
    Return

  End If

! ---  preliminary version at tau = 0.
! ---  should not matter since t(l) anyway tmin

  ri = r(1)/rtau23
  qh = q0/(ri**2)
  t(1) = teff*(.75D0*qh)**.25D0
  tauold = 0.D0
  qold = 0.D0
  dqdtau(1) = 0.D0

  Do k = 1, n - 1
    If (t(k)<tmin) t(k) = tmin
    taup = taul(k+1)
    taur = tau(k+1)

    If (taup>=.66D0 .And. tauold<.66D0) Then
      optfound = .True.
      rold = r(k)
      rp = r(k+1)
      dtaudr = (taup-tauold)/(rp-rold)
      r23 = rold + (.66D0-tauold)/dtaudr
      If (r23>rold .Or. r23<rp) Then
        Write (999, *) ' STOP: ERROR IN R(TAU-LUCY=2/3)'
        Stop ' ERROR IN R(TAU-LUCY=2/3)'
      End If
      srvmin = srnom*r(nsdiv)/r23
      Print *
      Print *, ' MESSAGE FROM TLUCY: SRVMIN/SRNOM = ', srvmin/srnom
      Print *
    End If

    q = qinf - (qinf-q0)*exp(-taur*gam)
    qprime = taup/taur*q
    qh = taup/taur*(q+taur)
    dqdtau(k+1) = (qold-qprime)/(tauold-taup)
    t(k+1) = teff*(.75D0*qh)**.25D0
    tauold = taup
    qold = qprime
  End Do

  If (ithyd>=2 .And. .Not. optfound) Then
    Write (999, *) ' STOP: TAU-LUCY = 2/3 NOT FOUND'
    Stop ' TAU-LUCY = 2/3 NOT FOUND'
  End If

  Return

End Subroutine

!***********************************************************************

!subroutines: complex ones
!model and related

!***********************************************************************

Subroutine model(tau, vmin1, optout, ithyd, xkap, ex, tact, xmmax, rtau23, &
  xnh, xne, vdiv, ndiv)

  Use :: nlte_type
  Use :: nlte_dim
  Use :: nlte_var, Only: iconver, modnam, optmodel

  Implicit None


! calculates atmospheric structure for analytical velocity laws
! and approx. photospheric structure


! analytical velocity law
! calculated values of ne, nh, tau, t  assume pure h/he atmosphere

! v(r) = vmax * (1.-b/r)**beta   ( r > srnom, i.e,
! v > vdiv * vsound)
! v(r) = vmin * exp((r-sr/hl)    ( r < srnom )

! clumping is now accounted for according to parameterization
! provided by input parameters (clf_field) and parameterization
! given in subroutine clumping
! -----------------------------------------------------------------------
! -----------------------------------------------------------------------

! input
! local input files: indat

! indat:
! teff, log g, rstar/rsun(nominal value)
! rmax,tmin
! mloss,v0,vmax,beta,vdiv
! hei
! optmod
! **
! rmax.....in stellar radii
! tmin.....minimum temperature for lucy's temp.stratification
! mloss....mass loss in m-sun/year
! v0,vmax..in km/s
! optmod...if optmod=.true., Kurucz/detail structure will be used
! to create structure (implies optmodel=.true.)

! ********************************************************************


! .. parameters ..
  Integer (i4b), Parameter :: nd = id_ndept
! ..
! .. scalar arguments ..
  Real (dp) :: ex, rtau23, vmin1, xkap, xmmax, vdiv
  Integer (i4b) :: ithyd, ndiv
  Logical :: optout
! ..
! .. array arguments ..
  Real (dp) :: tact(nd), tau(nd), xne(nd), xnh(nd)
! ..
! .. local scalars ..
  Real (dp) :: beta, hei, optmixed, rmax, rstarnom, teff, tmin, vdivold, vmax, &
    vmin, xlogg, xmloss, yhein
  Integer (i4b) :: nd1, nd2, nd3, nd4, nnn
  Logical :: optmod
! ..
! .. local arrays ..
  Real (dp) :: gradv(nd), r(nd), r1(nd), rho(nd), rho1(nd), td(nd), temp(nd), &
    temp1(nd), v(nd), v1(nd), wne(nd), clfac(nd), xne1(nd), xnh1(nd)

! ..
! .. external subroutines ..

  Character :: dumchar*30
  External :: modvl, prmodel
! ..

  If (nd<47) Then
    Write (999, *) ' STOP: ND TOO LOW'
    Stop ' ND TOO LOW'
  End If

  Open (1, File=trim(modnam)//'/INDAT.DAT', Status='OLD')
  Rewind 1
  Read (1, Fmt=*) dumchar
  Read (1, Fmt=*) optmod, optmod, nnn, nnn
  Read (1, Fmt=*) optmixed
  Read (1, Fmt=*) teff, xlogg, rstarnom
  Read (1, Fmt=*) rmax, tmin
  Read (1, Fmt=*) xmloss, vmin, vmax, beta, vdivold
  Read (1, Fmt=*) yhein, hei
  Read (1, Fmt=*) optmod
  Close (1)

  If (optmod .And. .Not. optmodel) Then
    Write (999, *) ' STOP: ERROR: OPTMOD.AND..NOT.OPTMODEL'
    Stop ' ERROR: OPTMOD.AND..NOT.OPTMODEL'
  End If

  If (tmin<2.D0) tmin = tmin*teff
  If (vmin1/=0.D0) vmin = vmin1

  nd1 = 10
  nd2 = 13
  nd3 = 16 !                        inclusive one extra point at .98 Rmax
! ND2 = 20
! ND3 = 9  !inclusive one extra point at .98 Rmax
! (one less than in versions before 8.6.1)
! ND1 = 10
! ND2 = 10
! ND3 = 10

! save 7 points for N5 + one extra point close to R=1)

  If (nd<51) Then
    nd4 = (nd-(nd1+nd2+nd3+7)) !    old version, no particular resolution of
!   transonic region
  Else
!   save additonal 4 points for transonic region (if vdiv not too large)
    nd4 = (nd-(nd1+nd2+nd3+11))
  End If

  If (nd4<0) Then
    Write (999, *) ' STOP: ERROR IN ND4'
    Stop ' ERROR IN ND4'
  End If

  Call modvl(teff, xlogg, rstarnom, rmax, tmin, xmloss, vmin, vmax, beta, &
    vdiv, yhein, hei, optmod, nd, nd1, nd2, nd3, nd4, r, v, gradv, rho, xne, &
    xnh, temp, tau, td, wne, r1, v1, rho1, xne1, xnh1, temp1, ithyd, xkap, ex, &
    tact, xmmax, rtau23, ndiv, clfac)

! REMEMBER: TAU=TAU-TH ONLY IN MICRO-CLUMPING APPROX.
  If (optout .And. iconver==1) Call prmodel(r, v, gradv, rho, xne, xnh, temp, &
    tau, clfac, teff, nd)

  Return

End Subroutine

!-----------------------------------------------------------------------

Subroutine modvl(teff, xlogg, rstarnom, rmax, tmin, xmloss, vmin, vmax, beta, &
  vdiv, yhein, hei, optmod, nd, nd1, nd2, nd3, nd4, r, v, gradv, rho, xne, &
  xnh, temp, tau, td, wne, r1, v1, rho1, xne1, xnh1, temp1, ithyd, xkap, ex, &
  tact, xmmax, rtau23, ns, clfac)

  Use :: nlte_var, Only: iconver
  Use :: nlte_type
  Use :: nlte_dim
  Use :: fund_const, Only: xmh => amh, akb, sigmae => sigmae_cross, &
    xmsu => xmsun, rsu => rsun
  Use :: nlte_opt, Only: opt_photclump

  Use :: princesa_var, Only: nat, abund, labat, weight

  Use :: nlte_var, Only: modnam, sr, srvmin, srnom, constt, nsdiv, h1, idum, &
    itmin, ggrav, te => teff, sigem, ckappa, xt1 => xt, a2te, delsig, delmu, &
    delpara

  Use :: photstruc !                includes XM

  Implicit None

! calculation of atmospheric structures for velocity law
! as described in Santolaya-Rey et al.,
! accounting now for clumping (see below)

! possibility to use Kurucz/Detail or Tlusty structure for photosphere

! output: first iteration: ne, temp1=temp as start values
! other iterations:    temp1=(tact,tnew),
! connected at sonic point
! input : in hydro iteration cycle: ne, tact

! in hydro: tnew, p (from diff eq., with approximated opacity)
! on output p is recalculated from tact to be
! intrinsically consistent

! after convergence, tnew has to be a good approximation of tact
! this is checked in routine lte by "check"


! !!!! recheck carefully, if more than h/he present!!!

! .. parameters ..
  Integer (i4b), Parameter :: nnde = id_ndept

! ..
! .. scalar arguments ..
  Real (dp) :: beta, ex, hei, rmax, rstarnom, rtau23, teff, tmin, vdiv, vmax, &
    vmin, xkap, xlogg, xmloss, xmmax, yhein

  Integer (i4b) :: ithyd, nd, nd1, nd2, nd3, nd4, ns
  Logical :: optmod, found, first
  Logical :: updated = .False.
! ..
! .. array arguments ..
  Real (dp) :: gradv(nnde), r(nnde), r1(nnde), rho(nnde), rho1(nnde), &
    tact(nd), tau(nnde), td(nnde), temp(nnde), temp1(nnde), v(nnde), v1(nnde), &
    wne(nnde), xne(nnde), xne1(nnde), xnh(nnde), xnh1(nnde), clfac(nnde)
! ..
! .. local scalars ..
  Real (dp) :: a2, auxi, b, c1, c2, del, del1, deltam, deltam1, deltar, diffv, &
    divrho, dm, dpx, eps1, gamma, h2, h3, hl, hs, p0, pi, r0, rc, rdiv2, rho0, &
    rhopho, rnew1, rpar, rper, rrr, rvmax, sige, srconst, summ, sumatom, &
    summass, t0, taumax, taumin, vdiv2, vdst, vs2, vsonic, x1, x2, xi, xm0, &
    xmu, xmu1, yhe, term, grad, pl, tl, rl, dmdr, dummy, rhol, vns, rns, &
    safety, integ, tau_e, rr1, rr2, det, v1test, v2test, rr0, vsmooth, &
    rhosmooth, vnsold

  Integer (i4b) :: i, ii, isi, ivs2, j, k, l, ll, m, n5, n6, nbad, nene, nok, &
    vinv, istart, iend, nr
! ..
! .. local arrays ..
  Real (dp) :: arho(id_ndept), brho(id_ndept), crho(id_ndept), p(id_ndept), &
    tnew(id_ndept), y(3), rhocl(id_ndept), xnecl(id_ndept), radacc(id_ndept), &
    zz1(4), zz2(4), vmat(4, 4), yy(4), gradvnew(id_ndept)

  Integer (i4b) :: index(id_ndept), indlud(4)
! ..
! .. external functions ..
  Real (dp) :: rhofunc, vfvl
  External :: rhofunc, vfvl
! ..
! .. external subroutines ..
  External :: derivs, derivs_kur, derivs_tlu, odeint, rkqc, spline, tlucy, &
    vindex
! ..
! .. intrinsic functions ..
! INTRINSIC ABS,ACOS,DBLE,EXP,INT,LOG,LOG10,SQRT
! ..
! .. data statements ..


  Data eps1/.5D0/
  Data first/.True./

! in order to account for clumping, only the last part of this
! routine is affected,
! a) since tau-scale of first part refers to
! tau_e alone (which remains unaffected -- clf cancels out)

! NOTE: THIS IS NO LONGER THE CASE, SINCE POSSIBILITY FOR POROSITY
! in e-scattering. Thus this scale might be modified accordingly.
! Not considered yet, since it is only a start-model -- and SO FAR
! it has not been any problems at all.
! Note by JO: tau_e only used as reference scale for setting up the grid
! MIGHT NEED TO BE UPDATED IN FUTURE.

! b) since photosphere is not affected by clumping

! ---- variable optstw incorporated (e. santolaya, 6-11-94) to allow the
! ---- input of a given atmosphere model


  nene = 0

  Do isi = 1, nat
    If (labat(isi)=='HE') nene = isi
  End Do

  If (nene==0) Then
    Print *, 'HE NOT FOUND - MODVL13'
    If (hei/=0.D0) Then
      Print *, ' NO HELIUM BUT ELSCAT-CONTRIBUTION!!!'
      Print *, ' EITHER SET HEI TO ZERO OR CALCULATE WITH HELIUM!!!'
!     comment the following line out for specific tests
      Write (999, *) ' STOP: INCONSISTENT TREATMENT OF HELIUM'
      Stop ' INCONSISTENT TREATMENT OF HELIUM'

    End If
    yhe = yhein
  Else
    yhe = abund(nene)
    If (abs(1.D0-yhe/yhein)>1.D-5) Then
      Write (999, *) ' STOP: ERROR IN YHE - MODVL13'
      Stop ' ERROR IN YHE - MODVL13'
    End If
  End If


  Write (999, Fmt=170) teff, xlogg, rstarnom, yhe
  Write (*, Fmt=170) teff, xlogg, rstarnom, yhe

  taumax = 20.D0
  te = teff
  ggrav = 10.D0**xlogg
  vdst = sqrt(2.D0*akb*teff/xmh)

  srnom = rstarnom*rsu
  vmax = vmax*1.D5
  vmin = vmin*1.D5
  h1 = xmloss*xmsu/3.9657D8

  Write (999, Fmt=180) xmloss, vmax/1.D5, vmin/1.D5
  Write (*, Fmt=180) xmloss, vmax/1.D5, vmin/1.D5

  sumatom = 0.D0
  summass = 0.D0

  Do k = 1, nat
    sumatom = sumatom + abund(k)
    summass = summass + weight(k)*abund(k)
  End Do

  If (nene==0 .And. yhe/=0.) Then
    sumatom = sumatom + yhe
    summass = summass + 4.D0*yhe
  End If

  summass = summass*xmh
  c1 = summass

! ---- for sigem, only contribution of helium and hydrogen

  c2 = (1.D0+hei*yhe)/c1
  sigem = sigmae*c2

  taumin = h1*sigem/(vmax*rmax*srnom)
  constt = taumin/(4.D0*rmax**2)


  If (.Not. optmod) Then

    modtype = 'FASTWD'

  Else !                            (i.e., if optmod)
!   IF THIS FEATURE IS ENABLED, TAKE CARE CONCERNING CORRECT TREATMENT OF
!   CLUMPING!!!
!   to be considered later; so far assumed that no photospheric clumping

!   ------- the model is read from structure file

    If (first) Then
      Open (4, File=trim(modnam)//'.struct', Status='OLD', Form='FORMATTED', &
        Err=100)
      modtype = 'KURUCZ'
      Print *
      Print *, ' ** STRUCTURE = KURUCZ **'
      Go To 130

100   Open (4, File=trim(modnam)//'.tlusty', Status='OLD', Form='FORMATTED', &
        Err=110)
      modtype = 'TLUSTY'
      Print *
      Print *, ' ** STRUCTURE = TLUSTY **'
      Go To 130


110   Open (4, File=trim(modnam)//'.michel', Status='OLD', Form='FORMATTED', &
        Err=120)
      modtype = 'MICHEL'
      Print *
      Print *, ' ** STRUCTURE = MICHEL **'
      Go To 130
120   Continue
      Write (999, *) ' STOP: Model structure (*.struct/*.tlusty/*.michel) &
        &not present - please create!'
      Stop ' Model structure (*.struct/*.tlusty/*.michel) not present &
        &- please create!'

130   Continue

      If (modtype=='KURUCZ') Then
        Read (4, *)
        Read (4, *) ndmod
        Allocate (xmext(ndmod), text(ndmod), pkur(ndmod), xnekur(ndmod), &
          gradkur(ndmod))
        Do i = 1, ndmod
          Read (4, *) xmext(i), text(i), pkur(i), xnekur(i), gradkur(i)
        End Do

      Else If (modtype=='TLUSTY') Then ! also for cmfgen input
        Read (4, *) ndmod
        Allocate (xmext(ndmod), text(ndmod), rhotlu(ndmod), xnetlu(ndmod))
        Do i = 1, ndmod
          Read (4, *) ii, xmext(i), dummy, dummy, text(i), xnetlu(i), &
            rhotlu(i)
        End Do
        If (ii/=ndmod) Then
          Write (999, *) ' STOP: WRONG INPUT FILE - CHECK!'
          Stop ' WRONG INPUT FILE - CHECK!'
        End If

      Else If (modtype=='MICHEL') Then
        Read (4, *) ndmod
        Allocate (rmich(ndmod), rmich1(ndmod), vmich(ndmod), rhomich(ndmod))
        Do i = 1, ndmod
          Read (4, *) rmich(i), dummy, vmich(i), dummy, dummy, rhomich(i)
          vmich(i) = vmich(i)*1.D5
!         PRINT*,I,RMICH(I),VMICH(I),RHOMICH(I)
!         consistency check for Mdot/(4.*Pi)
          dummy = (srnom*rmich(i))**2*vmich(i)*rhomich(i)
          Print *, dummy, h1
!         if this line is commented out, one can use mdots inconsistent
!         with the density stratification
          If (abs(1.-dummy/h1)>0.01) Then
            Write (999, *) &
              ' STOP: DENSITY(MICHEL) AND MDOT(INPUT) INCONSISTENT'
            Stop ' DENSITY(MICHEL) AND MDOT(INPUT) INCONSISTENT'
          End If
        End Do
!       consistency check for VINF
        If (abs(1.-vmich(ndmod)/vmax)>0.01) Then
          Write (999, *) ' STOP: VMAX(MICHEL) AND VINF(INPUT) INCONSISTENT'
          Stop ' VMAX(MICHEL) AND VINF(INPUT) INCONSISTENT'
        End If
        vmax = vmich(ndmod)


      Else
        Write (999, *) ' STOP: WRONG MODEL-TYPE'
        Stop ' WRONG MODEL-TYPE'
      End If

      Close (4)

      first = .False.
    End If

  End If

! ------- the model is calculated

  xmu = (1.D0+4.D0*yhe)/(2.D0+yhe*(1.D0+hei))

  If (ithyd/=1) Then
    Do l = 1, nd
!     clumping cancels out
      xmu1 = summass/(sumatom+xne(l)/xnh(l))/xmh

!     ---  delmu defined in order to correct a2te (vs(iso)**2/t)
!     ---  for ionisation effects

      delmu(l) = xmu/xmu1
    End Do
  End If

  gamma = sigem*1.8913D-15*te**4/ggrav

  If (gamma>=1.D0) Then
    Write (999, *) ' STOP: GAMMA .GE. 1'
    Stop ' GAMMA .GE. 1'
  End If
! ---- DO NOT ALLOW POROSITY-MEDIATED MASS LOSS FROM PHOTOSPHERE YET
! thus keep this condition

  hl = akb*te/xmu/xmh/(ggrav*(1.D0-gamma))

  Write (999, *)
  Write (999, *) ' GAMMA = ', gamma
  Write (999, *) ' PRESS. SCALE HEIGHT(TEFF) = ', hl/srnom, ' (IN NOM. RADII)'
  Write (999, *)
  Print *
  Print *, ' GAMMA = ', gamma
  Print *, ' PRESS. SCALE HEIGHT(TEFF) = ', hl/srnom, ' (IN NOM. RADII)'
  Print *

! a=(log(vmax/vmin)+sr/hl)/beta
! -----approximate formula
! yy=a-log(a)
! rc=beta*hl*yy
! r0=rc/(1.d0/yy+1.d0)

! -----division point at vdiv*vsound corresponding to srnom

  vsonic = sqrt(akb*te/xmh/xmu)
  vs2 = vdiv*vsonic
! check for too large Mdot which cannot handled by FASTWIND
  b = 1.D0 - (vsonic/vmax)**(1.D0/beta)

! integral is int_1^inf dr/[r^2*(1-b/r)^beta]

  If (beta/=1.D0) Then
    integ = (1.-b)**(-beta)*(b-1.+(1.-b)**beta)
    integ = integ/(b*(1.-beta))
  Else
    integ = -log(1.-b)/b
  End If
  If (integ<0.) Then
    Write (999, *) ' STOP: ERROR IN INTEG (MODVL)'
    Stop ' ERROR IN INTEG (MODVL)'
  End If

! TAU_E is optical depth in e-scat (micro) at sonic point
  tau_e = sigem*h1/(vmax*srnom)*integ
  Print *, ' TAU_E (MICRO) AT VSOUND = ', tau_e
  Print *
  If (tau_e>2.) Then
    Print *, &
      ' MASS LOSS TOO LARGE OR VINF TOO LOW TO BE HANDLED BY FASTWIND!!!'
    Print *, ' REDUCE MDOT AT LEAST TO ', xmloss*2./tau_e
    Print *, ' TO BE ON THE SAFE SIDE, REDUCE AT LEAST TO ', xmloss/tau_e
    Print *
    Print *, ' --------------------'
    Print *, ' NO MODEL CALCULATED'
    Print *, ' --------------------'
    Write (999, *) &
      ' MASS LOSS TOO LARGE OR VINF TOO LOW TO BE HANDLED BY FASTWIND!!!'
    Write (999, *) ' REDUCE MDOT AT LEAST TO ', xmloss*2./tau_e
    Write (999, *) ' TO BE ON THE SAFE SIDE, REDUCE AT LEAST TO ', &
      xmloss/tau_e
    Write (999, *)
    Write (999, *) ' --------------------'
    Write (999, *) ' NO MODEL CALCULATED'
    Write (999, *) ' --------------------'
    Write (999, *) ' STOP!'
    Stop
  End If

  If (vs2<vmin) Then
    Write (999, *) ' STOP: VDIV < VMIN'
    Stop ' VDIV < VMIN'
  End If

  If (ithyd==1) Then
    rc = srnom
  Else
    rc = srvmin
  End If

  sr = rc - hl*log(vs2/vmin)
  r0 = rc*(1.D0-(vs2/vmax)**(1.D0/beta))
! here and in the following R0 corrsponds to b, and RC -> R(NS) at VDIV
! in the 2nd part of this routine (photosphere), R0 corresponds to R(NS) at
! VDIV

  Print *, ' LOWER (SCALING RADIUS) AT ', sr/srnom, ' NOM. RADII'
  Print *, ' SRVMIN = ', rc/srnom, ' NOM. RADII'


  rhopho = h1/sr**2/vmin

! aufstellen des r-rasters ueber nd1 stuetzstellen,die
! bezueglich log tau aequidistant verteilt sind

  r(1) = rmax*srnom
  v(1) = vfvl(r(1), rc, r0, sr, vmax, vmin, beta, hl)
  rho(1) = rhofunc(r(1), v(1))

  r(nd1) = sr
  v(nd1) = vmin
  rho(nd1) = rhopho

140 Continue

  del1 = log10(taumax/taumin)/dble(nd1-1)
  h2 = 1.D0 - 10.D0**(-del1)

  Do k = 2, nd1 - 1
    sige = sigem*rho(k-1)
    r(k) = r(k-1) - taumin*10.D0**(del1*(k-1))*h2/sige

    If (r(k)<sr) Then
      taumax = 4.D0*taumax/5.D0
      Go To 140
    End If

    v(k) = vfvl(r(k), rc, r0, sr, vmax, vmin, beta, hl)
    rho(k) = rhofunc(r(k), v(k))
  End Do

! einschieben von nd2 punkten,so dass v "gleichmaessig
! ueber das r-raster verteilt wird.

  Do j = 0, nd2 - 1
    diffv = (v(1)-v(2))
    i = 2

    Do k = 3, nd1 + j
      h3 = (v(k-1)-v(k))
      If (h3>diffv) Then
        diffv = h3
        i = k
      End If
    End Do

    m = nd1 + j - i
    Call vindex(r, v, rho, m, i)
    r(i) = .5D0*(r(i)+r(i-1))
    v(i) = vfvl(r(i), rc, r0, sr, vmax, vmin, beta, hl)
    rho(i) = rhofunc(r(i), v(i))
  End Do

! einschieben von nd3-1 punkten,sodass rho "gleichmaessig"
! ueber das r-gitter verteilt wird

  Do j = 0, nd3 - 2
    divrho = rho(2)/rho(1)
    i = 2

    Do k = 3, nd1 + nd2 + j
      h3 = rho(k)/rho(k-1)
      If (h3>divrho) Then
        divrho = h3
        i = k
      End If
    End Do

    m = nd1 + nd2 + j - i
    Call vindex(r, v, rho, m, i)
    r(i) = .5D0*(r(i)+r(i-1))
    v(i) = vfvl(r(i), rc, r0, sr, vmax, vmin, beta, hl)
    rho(i) = rhofunc(r(i), v(i))
  End Do

! ---- einschieben von einem punkt bei .98*rmax

  i = 2
  m = nd1 + nd2 + (nd3-1) - i
  Call vindex(r, v, rho, m, i)
  r(i) = .98D0*r(1)
! for tests
! R(I) = .9D0*R(1)

  If (r(3)>r(i)) Then
    Write (999, *) ' STOP: R(3) > R(2)'
    Stop ' R(3) > R(2)'
  End If

  v(i) = vfvl(r(i), rc, r0, sr, vmax, vmin, beta, hl)
  rho(i) = rhofunc(r(i), v(i))

! ND4 points exclusively for the photosphere

  If (nd4>0) Then
    rpar = r(nd1+nd2+nd3-1)
    rper = r(nd1+nd2+nd3)
    deltar = (rpar-rper)/dble(nd4+1)
    Do ii = 1, nd4
      rnew1 = rpar - dble(ii)*deltar

      Do j = 1, nd1 + nd2 + nd3 + ii
        If (r(j)>rnew1) k = j + 1
      End Do

      m = nd1 + nd2 + nd3 + ii - k
      Call vindex(r, v, rho, m, k)
      r(k) = rnew1
      v(k) = vfvl(r(k), rc, r0, sr, vmax, vmin, beta, hl)
      rho(k) = rhofunc(r(k), v(k))
    End Do
  End If

! ---- addition of points at photospheric division
! -old:(between 2.*vs2 and .2*vs2)
! changed March 2016; form
! -----(between 2.*vs2 and .3*vs2)

  If (3.*vsonic>vmax) Then
    Write (999, *) ' STOP: 3 VSOUND > VMAX! CHANGE MODVL!'
    Stop ' 3 VSOUND > VMAX! CHANGE MODVL!'
  End If

  If (nd<51 .Or. (3.*vsonic-2.*vs2)<vsonic) Then
    n6 = 0
    n5 = nd - (nd1+nd2+nd3+nd4+1) ! old version or large vdiv
    If (nd<51) Then
      If (n5/=6) Then
        Write (999, *) ' STOP: ERROR IN N5'
        Stop ' ERROR IN N5'
      End If
    Else
      If (n5/=10) Then
        Write (999, *) ' STOP: ERROR IN N5'
        Stop ' ERROR IN N5'
      End If
    End If
  Else
    n6 = 4
    n5 = nd - (nd1+nd2+nd3+nd4+1) - n6 ! save 4 points for transonic region
    If (n5/=6) Then
      Write (999, *) ' STOP: ERROR IN N5'
      Stop ' ERROR IN N5'
    End If
  End If

  auxi = (2.D0*vs2/vmax)**(1.D0/beta)
  rdiv2 = r0/(1.D0-auxi)

  If (rdiv2<=rc) Then
    Write (999, *) ' STOP: ERROR IN RDIV2'
    Stop ' ERROR IN RDIV2'
  End If

  vdiv2 = 0.3D0*vs2

  If (vdiv2<=vmin) Then
    Write (999, *) &
      ' STOP: VMIN TOO LARGE - MODIFY INPUT (EITHER VMIN OR VDIV)'
    Stop ' VMIN TOO LARGE - MODIFY INPUT (EITHER VMIN OR VDIV)'
  End If
! modified March 2016 (nearest index not w.r.t. to v, but to r)
! otherwise, points become to close to each other

! nearest index to rc (corresponding to vs2)

! XI = (2.D0-1.D0)* (N5-1)/ (2.D0-.2D0)
! IVS2 = INT(XI)
! IF (XI-IVS2.GT..5D0) IVS2 = IVS2 + 1


  rvmax = sr + hl*log(vdiv2/vmin)

  If (rvmax>=rc) Then
    Write (999, *) ' STOP: ERROR IN 0.3 * VS2'
    Stop ' ERROR IN 0.3 * VS2'
  End If

  rpar = rdiv2
  rper = rvmax

! modified
  xi = (rpar-rc)*(n5-1)/(rpar-rper)
  ivs2 = int(xi)
  If (xi-ivs2>.5D0) ivs2 = ivs2 + 1

  deltar = (rpar-rper)/dble(n5-1)

! ----------------------------------------------------------------------------
! MARCH 2016
! remove points between rpar+deltar and rc (will be filled up with N5/2 points
! below),
! and redistribute them w.r.t. velocity
  istart = 0
  iend = 0

  nbad = nd1 + nd2 + nd3 + nd4
  Do i = 1, nbad
    If (r(i)<=rpar+deltar .And. r(i)>=rc) Then
      istart = i
      Exit
    End If
  End Do

  Do i = nbad, 1, -1
    If (r(i)<=rpar+deltar .And. r(i)>=rc) Then
      iend = i
      Exit
    End If
  End Do

  If (istart==0 .And. iend/=0) Then
    Write (999, *) ' STOP: ERROR IN ISTART, IEND (REMOVE IN MODVL)'
    Stop ' ERROR IN ISTART, IEND (REMOVE IN MODVL)'
  End If

  If (istart/=0) Then
    nr = iend - istart + 1
    If (nr<1) Then
      Write (999, *) ' STOP: ERROR IN NR (REMOVE IN MODVL)'
      Stop ' ERROR IN NR (REMOVE IN MODVL)'
    End If
!   remove points
    Do j = iend + 1, nbad
      r(j-nr) = r(j)
      v(j-nr) = v(j)
      rho(j-nr) = rho(j)
    End Do
!   re-distribute them w.r.t. velocity (in this case, in the outer regions)
    Do j = 0, nr - 1
      diffv = (v(1)-v(2))
      i = 2

      Do k = 3, nbad - nr + j
        h3 = (v(k-1)-v(k))
        If (h3>diffv) Then
          diffv = h3
          i = k
        End If
      End Do

      m = nbad - nr + j - i
      Call vindex(r, v, rho, m, i)
      r(i) = .5D0*(r(i)+r(i-1))
      v(i) = vfvl(r(i), rc, r0, sr, vmax, vmin, beta, hl)
      rho(i) = rhofunc(r(i), v(i))
    End Do
  End If
! ----------------------------------------------------------------------------


  Do ii = 0, n5 - 1
    rnew1 = rpar - dble(ii)*deltar

!   one point exactly at rc

    If (ii==ivs2) rnew1 = rc !      usually, at II=2 (for N5=6)

    Do j = 1, nd1 + nd2 + nd3 + nd4 + ii
      If (r(j)>rnew1) k = j + 1
    End Do

    m = nd1 + nd2 + nd3 + nd4 + ii - k
    Call vindex(r, v, rho, m, k)
    r(k) = rnew1
    v(k) = vfvl(r(k), rc, r0, sr, vmax, vmin, beta, hl)
    rho(k) = rhofunc(r(k), v(k))
  End Do


  If (n6>0) Then
!   N6 points to resolve transonic region from 2*VDIV2 to 3*VSOUND

    Do i = 1, nd - 1
      If (v(i)<3.*vsonic) Exit
    End Do

!   DIFFV=ABS(V(I)-3.*VSONIC)

!   IF(DIFFV.LT.ABS(V(I-1)-3.*VSONIC)) THEN
!   ISTART=I
!   ELSE
    istart = i - 1
!   ENDIF

    rpar = r(istart)

    Do i = istart, nd - 1
!     CHANGED IN V10.1 BECAUSE OF ACCURACY
      If (v(i)<(2.D0-1.D-12)*vs2) Exit
    End Do

    iend = i - 1
    rper = r(iend)

    If (rpar-rper<=0.) Then
      Write (999, *) ' STOP: RPAR < RPER IN TRANSONIC POINTS'
      Stop ' RPAR < RPER IN TRANSONIC POINTS'
    End If

    Do j = 0, n6 - 1

      diffv = log10(v(istart)/v(istart+1))
      i = istart + 1
      Do k = istart + 2, iend + j
        h3 = log10(v(k-1)/v(k))
        If (h3>diffv) Then
          diffv = h3
          i = k
        End If
      End Do

      m = nd1 + nd2 + nd3 + nd4 + n5 + j - i
      Call vindex(r, v, rho, m, i)
      r(i) = .5D0*(r(i)+r(i-1))
      v(i) = vfvl(r(i), rc, r0, sr, vmax, vmin, beta, hl)
!     PRINT*,RPAR,RPER,V(istart)/1.e5,v(iend)/1.e5,R(I),V(I)/1.E5
      rho(i) = rhofunc(r(i), v(i))
    End Do

  End If

! einschieben eines punktes nahe der photosphaere

  Call vindex(r, v, rho, 0, nd-1)
  r(nd-1) = r(nd) + eps1*(r(nd-2)-r(nd))
  v(nd-1) = vfvl(r(nd-1), rc, r0, sr, vmax, vmin, beta, hl)
  rho(nd-1) = rhofunc(r(nd-1), v(nd-1))

! DO K=1,ND
! PRINT*,K,R(K)/SR,V(K)/1.E5
! ENDDO

  Do k = 2, nd
    If (v(k)>v(k-1)) Then
      v(k) = v(k-1)
      rho(k) = rhofunc(r(k), v(k))
    End If
  End Do

  Do i = 1, nd
    rrr = r(i)
    If (rrr>=rc) Then
      gradv(i) = vmax*beta*r0/rrr/rrr*(1.D0-r0/rrr)**(beta-1.D0)
    Else
      gradv(i) = vmin/hl*exp((rrr-sr)/hl)
    End If
  End Do

! ---- end of calculation of the analytical model

! correction for photopheric density stratification
! warning!!! warning!!! warning!!!
! up to now neglected: advection term (minor)
! approximation of sigma-thomson in vsound (minor)

! So far, no optically thick clumping in photosphere allowed!

  Do l = 1, nd
    If (r(l)==rc) Go To 150
  End Do
  Write (999, *) ' STOP: ERROR: RC NOT FOUND'
  Stop ' ERROR: RC NOT FOUND'

150 Continue

  ns = l
  nsdiv = ns

  If (ithyd==1) Then
    rtau23 = srnom/sr
    rmax = r(1)/sr
  Else
    Print *
    Print *, ' PHOTOSPHERIC CORRECTION, M(RSTAR) = ', xmmax
    Print *
    xt1 = ex !                      for module nlte_var
    xt = ex !                       for module photstruc
    xmloss = xmloss*1.989D33/3.1557D7
    pi = acos(-1.D0)
    a2 = vsonic**2
    a2te = a2/te
    ckappa = sigem*xkap/a2te
    vmin = v(ns)
    r0 = r(ns)
    t0 = tact(ns)
    b = 1.D0 - (vmin/vmax)**(1.D0/beta)

    If (beta==1.D0) Then
      h2 = -log(1.D0-b)
    Else
      h2 = (1.D0-(vmin/vmax)**((1.D0-beta)/beta))/(1.D0-beta)
    End If

    rho0 = rho(ns)
    p0 = a2te*delmu(ns)*rho0*t0
    xm0 = xmloss/(4.D0*pi*r0*vmax)/b*h2

    Print *
    Print *, ' DIVISION AT ', vmin*1.D-5, ' KM/S (L = ', ns, ')'
    Print *, ' CORRESPONDING TO ', r0/srnom, ' NOM. RADII'
    Print *, ' AND LG(M0) = ', log10(xm0)
    Print *

    deltam = log10(xmmax/xm0)/(nd-ns-5)
    xm(ns) = log10(xm0)

    radacc = 0.

    Do l = ns + 1, nd
      deltam1 = deltam
      If (l<=ns+8) deltam1 = deltam/2.D0
      If (l<=ns+4) deltam1 = deltam/3.D0
      If (l<=ns+2) deltam1 = deltam/6.D0
      xm(l) = xm(l-1) + deltam1
    End Do

    Do l = ns, nd
      xm(l) = 10.D0**xm(l)
    End Do
!   JO Jan/April 2016: last point close to previous one, to improve diffusion
!   approx.
    If (nd>=61) xm(nd) = xm(nd-1)*1.1

    x1 = xm0
    y(1) = p0
    If (modtype=='TLUSTY') y(1) = rho0
    y(2) = t0
    y(3) = r0
    p(ns) = p0
    tnew(ns) = t0
!   RHO(NS) ALREADY DEFINED

!   THIS IS NEW, ACCOUNTING FOR TMIN
    Do l = ns + 1, nd
      If (tact(l)==t0) Cycle
      Exit
    End Do
    itmin = l - 1

!   PRINT*,'ITMIN=',ITMIN


!   LOOP OVER ALL XM
!   FROM V10.1 ON: WHILE LOOP, TO ALLOW FOR RESET OF L
!   -------------------------------------------------------------------------
    l = ns + 1

    Do While (l<nd+1)

      pl = y(1)
      If (modtype=='TLUSTY') rhol = y(1)
      tl = y(2)
      rl = y(3)
      x2 = xm(l)
      hs = (x2-x1)/20.D0
      idum = l

      If (modtype=='FASTWD' .Or. modtype=='MICHEL') Then
        Call odeint(y, 3, x1, x2, 1.D-5, hs, 0.D0, nok, nbad, derivs, rkqc)
      Else If (modtype=='KURUCZ') Then
        Call odeint(y, 3, x1, x2, 1.D-5, hs, 0.D0, nok, nbad, derivs_kur, &
          rkqc)
      Else If (modtype=='TLUSTY') Then
        Call odeint(y, 3, x1, x2, 1.D-5, hs, 0.D0, nok, nbad, derivs_tlu, &
          rkqc)
      Else
        Write (999, *) ' STOP: WRONG MODEL TYPE'
        Stop ' WRONG MODEL TYPE'
      End If

      p(l) = y(1)
      tnew(l) = y(2)
      r(l) = y(3)
      If (modtype=='TLUSTY') Then
        rho(l) = y(1)
        p(l) = rho(l)*(a2te*delmu(l)*tnew(l))
      End If
!     PRINT*,'TLUSTY',L,log10(XM(L)),R(L)/SR,log10(RHO(L)),TNEW(L),P(L)

!     one point just below SRNOM
!     always for standard calculation (OPT_PHOTCLUMP=.FALSE.) or from ITHYD >
!     3
!     if photospheric clumping (otherwise, photosphere too close to division
!     point)
      If (r(l)<srnom .And. r(l-1)>srnom .And. (.Not. opt_photclump .Or. ithyd> &
        3)) Then
        If (l<=ns+4) Then
          Print *
          Print *, ' WARNING! PHOTOSPHERE AT ', l, &
            ' VERY CLOSE TO DIVISION POINT AT', ns
          Print *
          Write (999, *)
          Write (999, *) ' WARNING! PHOTOSPHERE AT ', l, &
            ' VERY CLOSE TO DIVISION POINT AT', ns
          Write (999, *)
        End If
!       new try
        safety = 1.1D0
160     dmdr = (xm(l)-xm(l-1))/(r(l-1)-r(l))
        deltam = safety*dmdr*(r(l-1)-srnom) ! including safety factor
        If (deltam<0.) Then
          Write (999, *) ' STOP: ERROR IN DELTAM < 0'
          Stop ' ERROR IN DELTAM < 0'
        End If
!       FROM V10.1 ON: AVOID TOO SMALL SEPARATIONS (except L - 1 = NS, i.e.,
!       last point = starting point)
        If (deltam/xm(l-1)<0.1 .And. l-1/=ns) Then ! FACTOR 1.1
!         SKIP LAST POINT, OVERWRITE IT
          xm(l-1) = xm(l-1) + deltam
!         RESET INDEX, INITIALIZE ODE
          l = l - 1
          x1 = xm(l-1)
          y(1) = p(l-1)
          If (modtype=='TLUSTY') y(1) = rho(l-1)
          y(2) = tnew(l-1)
          y(3) = r(l-1)
          idum = l
        Else
!         ACCEPT LAST POINT, RECALCULATE NEW POINT
          xm(l) = xm(l-1) + deltam
          y(1) = pl
          If (modtype=='TLUSTY') y(1) = rhol
          y(2) = tl
          y(3) = rl
        End If
        x2 = xm(l)
        hs = (x2-x1)/20.D0

        If (modtype=='FASTWD' .Or. modtype=='MICHEL') Then
          Call odeint(y, 3, x1, x2, 1.D-5, hs, 0.D0, nok, nbad, derivs, rkqc)
        Else If (modtype=='KURUCZ') Then
          Call odeint(y, 3, x1, x2, 1.D-5, hs, 0.D0, nok, nbad, derivs_kur, &
            rkqc)
        Else If (modtype=='TLUSTY') Then
          Call odeint(y, 3, x1, x2, 1.D-5, hs, 0.D0, nok, nbad, derivs_tlu, &
            rkqc)
        Else
          Write (999, *) ' STOP: WRONG MODEL TYPE'
          Stop ' WRONG MODEL TYPE'
        End If

        p(l) = y(1)
        tnew(l) = y(2)
        r(l) = y(3)
        If (modtype=='TLUSTY') Then
          rho(l) = y(1)
          p(l) = rho(l)*(a2te*delmu(l)*tnew(l))
        End If
        If (r(l)/srnom>=1.0D0) Then
!         NEW PHILOSOPHY FROM V10.1.3 ON
          If (safety==1.1D0) Then
            Print *, ' WARNING!!!  PROBLEMS WITH POINT CLOSE TO PHOTOSPHERE!'
            Write (999, *) &
              ' WARNING!!!  PROBLEMS WITH POINT CLOSE TO PHOTOSPHERE!'
            safety = 1.2
            Go To 160
          End If

          Print *, xm(l-2), xm(l-1), xm(l)
          Print *, r(l-2)/srnom, r(l-1)/srnom, r(l)/srnom
          Print *, r(l)/srnom
          Write (999, *) &
            ' STOP: PROBLEMS WITH POINT CLOSE TO PHOTOSPHERE CANNOT BE CURED!'
          Stop ' PROBLEMS WITH POINT CLOSE TO PHOTOSPHERE CANNOT BE CURED!'
        End If
!       define new xm grid (equidistant)
        deltam = log10(xmmax/x2)/(nd-l)
        deltam = 10.D0**deltam
        Do ll = l + 1, nd
          xm(ll) = xm(ll-1)*deltam
        End Do
!       JO Jan/April 2016: last point close to previous one, to improve
!       diffusion approx.
        If (nd>=61) xm(nd) = xm(nd-1)*1.1

      End If
      x1 = x2
!     re-calculate R^2*GRAD (OK, checked!!!)
      term = sigem + ckappa*delpara(l)/delmu(l)*p(l)*y(2)**(-xt-1.D0)
!     note that ckappa includes sigmae
      term = term*delsig(l)
      grad = 5.67D-5/2.9979D10*teff**4*term
      radacc(l) = grad

      l = l + 1
    End Do

!   ---- while loop finished

    Do l = ns + 1, nd
      If (p(l)<0.D0) Then
        Write (999, *) ' STOP: PRESSURE NEGATIVE!'
        Stop ' PRESSURE NEGATIVE!'
      End If
      If (modtype/='TLUSTY') rho(l) = p(l)/(a2te*delmu(l)*tnew(l))
!     for tests; note that differences in grad should scale with
!     differences in chiross/rho, where chiross corresponds to TERM below
!     TERM = SIGEM*(1.D0+XKAP*DELPARA(L)*RHO(L)*TNEW(L)**(-XT))
!     TERM = TERM*DELSIG(L)*RHO(L)
!     GRAD=5.67D-5*TEFF**4/2.9979D10*TERM/RHO(l)*(SRNOM/R(L))**2
!     print*,l,' ',p(l),' ',tnew(l),' ',tact(l),' ',delmu(l),' ', &
!     &           rho(l),' ',LOG10(TERM)
!     print*,l,' ',grad,' ',radacc(l)
    End Do

    sr = r(nd)
    rmax = r(1)/sr
    rtau23 = srnom/sr
    srvmin = r0
    If (1.D0-1.D0/rtau23>.1D0) Print *, ' WARNING!!! EXTENDED PHOTOSPHERE'
    If (1.D0-1.D0/rtau23>.1D0) Write (999, *) &
      ' WARNING!!! EXTENDED PHOTOSPHERE'

    Print *, ' LOWER (PHOT.) RADIUS AT ', 1.D0/rtau23, ' NOM. RADII'
    Print *, ' MAX. RADIUS (RMAX) = ', rmax, ' (IN SR)'

    Do l = ns + 1, nd
      v(l) = xmloss/(4.D0*pi*r(l)**2*rho(l))
    End Do

!   modified to ensure smooth transition

!   new version
    Do l = ns, nd - 1
      dm = log(v(l)/v(l-1))/(r(l)-r(l-1))
      dpx = log(v(l)/v(l+1))/(r(l)-r(l+1))
      gradv(l) = .5D0*(dm+dpx)*v(l)
    End Do

    l = nd
    gradv(l) = log(v(l)/v(l-1))/(r(l)-r(l-1))*v(l)

  End If

! -------------
  If (modtype=='MICHEL') Then

!   OVERWRITE THE WIND-STRUCTURE WITH MICHEL'S VALUES, AFTER RENORMALIZING
    vns = v(ns)
    Do l = 1, ndmod - 1
      If (vmich(l)<vns .And. vmich(l+1)>=vns) Exit
    End Do
    rns = rmich(l) + (rmich(l+1)-rmich(l))/(vmich(l+1)-vmich(l))*(vns-vmich(l) &
      )
    Print *, vns, rns, r(ns), srnom
!   THIS POINT NOW SHOULD CORRESPOND TO R(NS) = SRNOM in the first iteration,
!   and to SRVMIN later on (until convergence)

!   print*,'renormalization!'
!   print*,rns,rns/srnom,l

!   RENORMALIZATION
    Do l = 1, ndmod
      rmich1(l) = rmich(l)/rns*r(ns)
    End Do

    If (rmich1(ndmod)<=r(1)) Then
      Write (999, *) &
        ' STOP: OUTER RADIUS (MICHEL) SMALLER THAN RMAX, DECREASE RMAX'
      Stop ' OUTER RADIUS (MICHEL) SMALLER THAN RMAX, DECREASE RMAX'
    End If

    Do k = 1, ns - 1
      Do l = 1, ndmod - 1
        If (rmich1(l)<r(k) .And. rmich1(l+1)>=r(k)) Exit
      End Do
      vns = vmich(l) + (vmich(l+1)-vmich(l))/(rmich1(l+1)-rmich1(l))*(r(k)- &
        rmich1(l))
      v(k) = vns
      rho(k) = h1/(r(k)**2*v(k))
      Print *, 'michel', k, r(k)/srnom, v(k), rho(k)
    End Do
!   DVDR
    Call deriv_3p(gradv, r, v, nd)

  End If
! -------------

  Do k = 1, nd
    r(k) = r(k)/sr
    v(k) = v(k)/vmax
    gradv(k) = gradv(k)*sr/vmax
    If (ithyd>1) Then
      p(k) = rho(k)*(a2te*delmu(k)*tact(k))
    End If
  End Do

  r(nd) = 1.D0

! -------------
! JO Apri16: now we 'smooth' the velocity at the transition point,
! requiring a 3rd degree polynom and prescribed v, gradv above and below
! transition point (4 unknowns for the parabola, four conditions to fix them)

  If (ithyd>1) Then
    rr1 = r(ns-1)
    rr2 = r(ns+1)

    zz1 = (/ rr1**3, rr1**2, rr1, 1.D0 /)
    zz2 = (/ rr2**3, rr2**2, rr2, 1.D0 /)
    vmat(1, :) = zz1
    vmat(2, :) = zz2
    vmat(3, :) = (/ 3.D0*rr1**2, 2.D0*rr1, 1.D0, 0.D0 /)
    vmat(4, :) = (/ 3.D0*rr2**2, 2.D0*rr2, 1.D0, 0.D0 /)

    yy = (/ v(ns-1), v(ns+1), gradv(ns-1), gradv(ns+1) /)

    Call ludcmp(vmat, 4, 4, indlud, det)
    Call lubksb(vmat, 4, 4, indlud, yy)

!   test precision
    v1test = dot_product(yy, zz1)
    If (abs(1.D0-v1test/v(ns-1))>1.D-4) Then
      Write (999, *) ' STOP: PRECISION OF V1 NOT SUFFICIENT (ROUTINTE MODVL)'
      Stop ' PRECISION OF V1 NOT SUFFICIENT (ROUTINTE MODVL)'
    End If
    v2test = dot_product(yy, zz2)
    If (abs(1.D0-v2test/v(ns+1))>1.D-4) Then
      Write (999, *) ' STOP: PRECISION OF V2 NOT SUFFICIENT (ROUTINTE MODVL)'
      Stop ' PRECISION OF V2 NOT SUFFICIENT (ROUTINTE MODVL)'
    End If

!   calculate smoothed velocity at NS
    rr0 = r(ns)
    zz1 = (/ rr0**3, rr0**2, rr0, 1.D0 /)
    vsmooth = dot_product(yy, zz1)
    Print *
    Print *, 'OLD/NEW V-VALUES AROUND NS:'
    Write (*, Fmt='(3(F12.6,2X))') v(ns-1), v(ns), v(ns+1)
    Write (*, Fmt='(3(F12.6,2X))') v(ns-1), vsmooth, v(ns+1)

!   recalculate rho and p, assuming r^2*rho*v = const
    rhosmooth = rho(ns)*v(ns)/vsmooth
    vnsold = v(ns)
    v(ns) = vsmooth
    p(ns) = p(ns)/rho(ns)*rhosmooth
    rho(ns) = rhosmooth

    If (modtype=='MICHEL') Then
      Call deriv_3p(gradvnew, r, v, nd)
    Else
      gradvnew = gradv
      Do l = ns - 1, ns + 1
        dm = log(v(l)/v(l-1))/(r(l)-r(l-1))
        dpx = log(v(l)/v(l+1))/(r(l)-r(l+1))
        gradvnew(l) = .5D0*(dm+dpx)*v(l)
      End Do
    End If
!   for tests
    Print *, 'OLD/NEW DVDR-VALUES AROUND NS:'
    Write (*, Fmt='(5F12.6,2X)') gradv(ns-2:ns+2)
    Write (*, Fmt='(5F12.6,2X)') gradvnew(ns-2:ns+2)

    gradv = gradvnew
  End If

! ----------------------------------------------------------------------------

! ----from here on, we have to consider clumping
! RHOCL is the overdense rho (in the average clump), whilst
! RHO is the mean rho, consistent with Mdot
  Call clumping(r, v, rho, clfac, rhocl, rtau23, nd, ns, vmax, gradv)

  Do l = 1, nd
    If (ithyd==1) xne(l) = c2*rhocl(l) ! corrected
    xnh(l) = rhocl(l)/c1 !          corrected
  End Do

! assumed that xne for ithyd ne 1 has been corrected

  xnecl = xne/clfac !               corrected for integration of optical depth
! in the following and until further evidence, we use the micro-clumped
! tau_e since this is presumably used for start-values only
! JO check

  Call spline(r, xnecl, arho, brho, crho, nd)
  tau(1) = taumin
  summ = taumin

  Do l = 1, nd - 1
    del = r(l+1) - r(l)
    summ = summ - sigmae*sr*del*(xnecl(l)+del*(arho(l)/2.D0+del*(brho(l)/3.D0+ &
      del*crho(l)/4.D0)))
    tau(l+1) = summ
  End Do

  Print *
  Print *, ' TAU-THOMSON (MICRO-CLUMPED) AT RSTAR = ', tau(nd)
  Print *


! ---- arho: tau-lucy(electron)
! ---- brho: tau-electron

  srconst = sr*rtau23**2

  Do l = 1, nd
    crho(l) = srconst*sigmae*xnecl(l)/r(l)/r(l)
  End Do

  summ = 0.D0
  arho(1) = 0.D0
  Do l = 1, nd - 1
    del = r(l+1) - r(l)
    summ = summ - del*5.D-1*(crho(l)+crho(l+1))
    arho(l+1) = summ
  End Do

  Do l = 1, nd
    crho(l) = sr*sigmae*xnecl(l)
  End Do

  summ = 0.D0
  brho(1) = 0.D0

  Do l = 1, nd - 1
    del = r(l+1) - r(l)
    summ = summ - del*5.D-1*(crho(l)+crho(l+1))
    brho(l+1) = summ
  End Do

! changed: here, tlucy uses tau_e only. Thus, rtau23 is no longer found
! for recombined models. Necessary only in first iteration anyway
  If (ithyd<2) Call tlucy(r, arho, brho, temp, nd, tmin, teff, rtau23, ithyd)


  Do l = 1, nd
    index(l) = l
  End Do

! in case, adapt vdiv to prevent density inversion
  If (teff>=12000. .And. ithyd>=3 .And. vdiv<0.9) Then

    vinv = 0
!   JO July 2016: this test has to be performed with unsmoothed v(ns)=VNSOLD
    v(ns) = vnsold
    Do l = nd - 1, 1, -1
      If (v(l)<v(l+1)) Then
        vinv = l + 1
        Exit
      End If
    End Do

    If (vinv/=0) Then
      Do l = vinv, 1, -1
        If (v(l)>v(vinv)) Then
          auxi = vdiv !             old value
          vdiv = 1.2*v(l)*vmax/vsonic ! including safety factor 1.2
          vdiv = max(vdiv, 1.2*auxi) ! this is the new version; ensure that
!         VDIV
!         becomes always larger (not smaller)
          vdiv = min(vdiv, 0.9D0) ! not beyond vsound
          iconver = 0
          Exit
        End If
      End Do
      Print *, ' WARNING! WARNING! WARNING!'
      Print *, ' VDIV ADAPTED, NEW VALUE = ', vdiv
      Print *
      Write (999, *) ' WARNING! WARNING! WARNING!'
      Write (999, *) ' VDIV ADAPTED, NEW VALUE = ', vdiv
      Write (999, *)
    End If
!   V(NS)=VNSOLD
!   JO August 2016: I think this was a typo
    v(ns) = vsmooth
  End If

  found = .False.
  Do l = 1, nd
    If (gradv(l)<=0.D0) Then
      found = .True.
      gradv(l) = -gradv(l)
      Print *, ' NEGATIVE VEL. GRADIENT -- DENSITY INVERSION! AT N =', l
      Print *, ' ACTUAL VDIV = ', vdiv
      If (iconver==1) Then
        Print *, ' NEGATIVE VEL. GRADIENT -- DENSITY INVERSION! AT N =', l
        Write (999, *) ' NEGATIVE VEL. GRADIENT -- DENSITY INVERSION! AT N =', &
          l
      End If
    End If
  End Do
  If (found .And. iconver==1) Then
    Print *, ' ACTUAL VDIV = ', vdiv
    Print *, ' INCREASE LOG G???'
    Write (999, *) ' ACTUAL VDIV = ', vdiv
    Write (999, *) ' INCREASE LOG G???'
    Write (999, *) ' PROCEED AT OWN RISK!!!!'
  End If


! ---- file 'model' is written: XNE, XNH include CLFAC
! (RHO average density without CLFAC)
! from version 9.0 on: NS as last entry
! from version 10.0 on: UPDATED and XT as last entries

  Rewind 20
! UPDATED always .false. here

  Write (20) teff, ggrav, sr, yhe, xmu, vmax, xmloss, beta, (r(i), i=1, nd), &
    (v(i), i=1, nd), (gradv(i), i=1, nd), (rho(i), i=1, nd), &
    (xne(i), i=1, nd), (xnh(i), i=1, nd), (clfac(i), i=1, nd), &
    (index(i), i=1, nd), (p(i), i=1, nd), ns, updated, xt

! store photospheric rad. acceleration
  Open (61, File=trim(modnam)//'/GRAD.OUT', Status='UNKNOWN')
  Do l = 1, nd
    Write (61, *) 'RADACC ', l, ' ', radacc(l)
  End Do
  Close (61)


  If (ithyd>1) Then

    Do l = 1, ns
      temp(l) = tact(l)
    End Do

    Do l = ns + 1, nd
      temp(l) = tnew(l)
    End Do

  End If

  Do ll = 1, nd
    l = index(nd+1-ll)
    temp1(ll) = temp(l)
    td(ll) = 1.2D0*sigmae*xnh(l)*vdst/(gradv(l)*vmax/sr)
    wne(ll) = xne(l)/(.5D0*(1.D0-sqrt(1.D0-1.D0/r(l)/r(l))))
    r1(ll) = r(l)
    v1(ll) = v(l)*vmax*1.D-5
    rho1(ll) = rho(l)
    xnh1(ll) = xnh(l)
    xne1(ll) = xne(l)/xnh(l)
!   write(*,190) l,td(ll),wne(ll),r1(ll),v1(ll),rho1(ll),xne1(ll)
  End Do

! open(1,file='stwithd',form='unformatted')
! rewind 1
! write(1) (td(i),i=1,nx),(wne(i),i=1,nx),(r1(i),i=1,nx),
! *         (v1(i),i=1,nx),(rho1(i),i=1,nx)
! close(1)

! ---- file 'temp' is written

  Rewind 21
  Write (21, *) r1, temp1

  Return

170 Format ('  TEFF = ', F10.1, '  LOG G = ', F5.2, '  RSTAR(NOMINAL)', ' = ', &
    F5.1, ' YHE = ', F6.2, /)
180 Format ('  MLOSS/YEAR = ', E12.6, '  VMAX = ', F7.1, '  VMIN = ', F11.7)
190 Format (I3, 6(2X,E10.5))

End Subroutine

!-----------------------------------------------------------------------

Subroutine clumping(r, v, rho, clfac, rhocl, rtau23, nd, ns, vmax, dvdr)


! calculates clumping factor clfac, CAN BE MODIFIED BY USER
! parameters (up to ten) provided in array clf_field
! number of parameters given by NPARA_CLF

! first parameter MUST be representative clumping factor, i.e. in any case

! MIN(CLFAC) <= CLF_FIELD(1) <= MAX(CLFAC)

! can be used, of course, do indicate that no clumping is present
! (this is NOT a must, however)

! User can change parameterization
! (as a function of V = v(r)/vinf or R = r/SR, with SR INNER radius)
! RTAU23 is ratio of SRNOM/SR, i.e., if you want to refer to
! the nominal radius, all values of r have to be divided by RTAU23
! (in this case then, r(nd) < 1 finally) => available in array RSCAL

! NS is the index of the transition point between wind and photosphere, can
! be used as well

! output values are clfac, rhocl=rho*clfac (overdensity in clumps)

! adapted for optically thick clumping by JS
! NOTE: variables HPOR, FVEL, FVOL *ONLY* needed for printing oputputs,
! ALL calculations work on TCL_FAC_LINE, TCL_FAC_CONT, and FIC -- see notes
! NOTE2: In Sundqvist+ 2016, TCL_FAC_LIN = G_l , TCL_FAC_CONT = G_c

! bug in conversion from fcl/fic to fvol identified March 2018

  Use :: nlte_type
  Use :: nlte_var, Only: clf_field, npara_clf, sr, modnam
  Use :: nlte_opt, Only: opt_photclump
  Use :: photstruc, Only: clf_const

  Use :: nlte_porvor, Only: optthick, fic_field, fvel_field, hpor_field, fic, &
    tcl_fac_line, tcl_fac_cont, fvel, fvol, hpor

  Implicit None

  Integer (i4b), Intent (In) :: nd, ns
  Real (dp), Intent (In) :: rtau23, vmax
  Real (dp), Dimension (nd), Intent (In) :: r, v, rho, dvdr
  Real (dp), Dimension (nd), Intent (Out) :: clfac, rhocl

  Real (dp) :: clf_representative ! always
  Real (dp) :: clf, vclstart, vclmax, clfmin, clfmax ! for NPARA_CLF=3
  Real (dp) :: dummy !              for NPARA_CLF=4
  Real (dp) :: clfmid, rmid, clfout, rout ! for NPARA_CLF=5
  Real (dp) :: cl1, cl2, cl3, cl4 ! for NPARA_CLF=6
  Real (dp) :: rin, rfar, clfin, clffar ! additionally for NPARA_CLF=9, 11
  Real (dp) :: rphot, clfphot !     additionally for NPARA_CL=11
  Real (dp) :: fvor

  Real (dp), Dimension (nd) :: rscal ! automatic array, in case we need it

  Integer (i4b) :: i, imax, istart, imid
  Integer (i4b) :: iphot, iin, iout ! additionally for NPARA_CLF=11

  clf_representative = clf_field(1)
  rscal = r/rtau23 !                scaled to nominal radius Rstar

  If (clf_representative<1.D0) Then
    Write (999, *) ' STOP: REPRESENT. VALUE OF CLF < 1!'
    Stop ' REPRESENT. VALUE OF CLF < 1!'
  End If

! ------------------------------------------------------------
! CHANGE IN CASE
  If (opt_photclump .And. npara_clf/=3) Then
    Write (999, *) ' STOP: PHOTOSPHERIC CLUMPING .AND. NPARA NE 3'
    Stop ' PHOTOSPHERIC CLUMPING .AND. NPARA NE 3'
  End If

  If (npara_clf==3) Then

    clf = clf_representative
    vclstart = clf_field(2)
    vclmax = clf_field(3)
!   "old" approach as used from version 8.2 on (linear increase from 1 to CLF
!   in between vclstart and vclmax
!   vor v < vclstart, clfac = 1, for v > vclmax, clfac = clf

!   clf = 1 gives unclumped model
!   clf = x and vclstart = vclmax = 0. gives constant clumping factor = x
!   EVERYWHERE

    If (vclstart>vclmax) Then
      Write (999, *) ' STOP: VCLSTART > VCLMAX!'
      Stop ' VCLSTART > VCLMAX!'
    End If

    If (clf==1.D0) Then

      clfac = 1.D0

    Else

      Do i = 1, nd
        If (v(i)<vclmax) Exit
      End Do

      imax = max0(1, i-1)
      Do i = imax, nd
        If (v(i)<vclstart) Exit
      End Do
      istart = i - 1 !              works also for vclstart = 0. (completely
!     clumped atmosphere)

      If (istart<imax) Then
        Print *, istart, imax
        Write (999, *) ' STOP: ISTART < IMAX IN CLUMPING'
        Stop ' ISTART < IMAX IN CLUMPING'
      End If

!     photospheric clumping now allowed except for tests (OPT_PHOTCLUMP)
      If (.Not. opt_photclump) Then
        If (istart>=ns) Then
          Print *, istart, ns, v(ns)
          Write (999, *) ' STOP: ISTART < NS IN CLUMPING'
          Stop ' ISTART < NS IN CLUMPING'
        End If
      End If

      clfac(1:imax) = clf

      Do i = imax + 1, istart
!       linear increase
        clfac(i) = 1. + (clf-1.)/(vclmax-vclstart)*(v(i)-vclstart)
      End Do

      clfac(istart+1:nd) = 1.D0
    End If

!   except for tests, no photospheric clumping
    If (.Not. opt_photclump) Then
      clfac(ns:nd) = 1.D0
    Else
      clfmin = minval(clfac)
      clfmax = maxval(clfac)
      If (clfmin/=clfmax) Then
        Write (999, *) ' STOP: PHOTOSPHERIC CLUMPING, BUT CLF NOT CONSTANT'
        Stop ' PHOTOSPHERIC CLUMPING, BUT CLF NOT CONSTANT'
      End If
      clf_const = clfmin
    End If

!   for tests
!   PRINT*
!   PRINT*,' CLUMPING FACTOR'
!   DO I=1,ND
!   PRINT*,I,V(I),CLFAC(I)
!   ENDDO

!   ------------------------------------------------------------
  Else If (npara_clf==4) Then
    dummy = clf_field(4)

    clfout = clf_field(2)
    vclstart = clf_field(3) !       in km/s

    If (dummy==0.) Then
!     either Exponential law a la Hillier; input filling factor, vcl, dummy=0.
      clfac = clfout + (1.-clfout)*exp(-(v*vmax*1.D-5/vclstart))

    Else
!     Paco's law returning to an unclumped outer wind
      clfac = clfout + (1.-clfout)*exp(-(v*vmax*1.D-5/vclstart)) + &
        (1.-clfout)*exp((v-v(1))*vmax*1.D-5/dummy)
    End If

    clfac = 1./clfac
!   no photospheric clumping
    clfac(ns:nd) = 1.D0

!   for tests
    Print *
    Print *, ' CLUMPING FACTOR'
    Do i = 1, nd
      Print *, i, rscal(i), v(i)*vmax*1.D-5, clfac(i)
    End Do

!   ------------------------------------------------------------
  Else If (npara_clf==5) Then
!   preliminary parameterization of clfac starting from unity (at NS)
!   to clmid (at rmid), reaching clout (at rout).
!   For r> rout, clfac is set to clout, in between linearly interpolated
!   idea: maximum value (around 5 to 9 at clmid (roughly at Rstar = 2),
!   then declining, following the results by Jo, Nevy and Salvo from the
!   combined IR/radio analyis

!   rmid and rout in units of nominal radius

    clfmid = clf_field(2)
    rmid = clf_field(3)
    clfout = clf_field(4)
    rout = clf_field(5)
    If (clfmid<1.D0 .Or. clfout<1.D0) Then
      Write (999, *) ' STOP: CLFMID OR CLFOUT < 1'
      Stop ' CLFMID OR CLFOUT < 1'
    End If

    If (rmid>rout) Then
      Write (999, *) ' STOP: RMID > ROUT!'
      Stop ' RMID > ROUT!'
    End If

    If (clfmid==1.D0 .And. clfout==1.D0) Then

      clfac = 1.D0

    Else

      Do i = 1, nd
        If (rscal(i)<rout) Exit
      End Do
      imax = max0(1, i-1)

      Do i = imax, nd
        If (r(i)<rmid) Exit
      End Do
      imid = i - 1

      istart = ns
!     PHOTOSPHERIC CLUMPING NOT ALLOWED
      If (imid>=istart) imid = istart - 1

      If (imid<=imax) Then
        Print *, imid, imax
        Write (999, *) ' STOP: IMID LE IMAX IN CLUMPING'
        Stop ' IMID LE IMAX IN CLUMPING'
      End If

      clfac(1:imax) = clfout

      Do i = imax + 1, imid
!       linear interpol.
        clfac(i) = clfout + (clfmid-clfout)/(rmid-rout)*(rscal(i)-rout)
      End Do

      Do i = imid + 1, istart - 1
!       linear increase
        clfac(i) = 1. + (clfmid-1.)/(rmid-rscal(ns))*(rscal(i)-rscal(ns))
      End Do

      clfac(istart:nd) = 1.D0
    End If
!   for tests
!   PRINT*
!   PRINT*,' CLUMPING FACTOR'
!   DO I=1,ND
!   PRINT*,I,RSCAL(I),CLFAC(I)
!   ENDDO

!   ------------------------------------------------------------
  Else If (npara_clf==6) Then

!   Paco's law with CL1 to CL4; one more (dummy) parameter than
!   actually needed, to allow for 6 parameters


    cl1 = clf_field(2)
    cl2 = clf_field(3)
    cl3 = clf_field(4)
    cl4 = clf_field(5)
    dummy = clf_field(6)


    If (dummy/=0.) Then
      Write (999, *) ' STOP: CLF_FIELD(6) NE 0, MODIFY!'
      Stop ' CLF_FIELD(6) NE 0, MODIFY!'
    End If

    clfac = cl1 + (1.-cl1)*exp(-(v*vmax*1.D-5/cl2)) + &
      (cl4-cl1)*exp((v-v(1))*vmax*1.D-5/cl3)

    clfac = 1./clfac
!   no photospheric clumping
    clfac(ns:nd) = 1.D0

!   for tests
    Print *
    Print *, ' CLUMPING FACTOR'
    Do i = 1, nd
      Print *, i, rscal(i), v(i)*vmax*1.D-5, clfac(i), 1./clfac(i)
    End Do

!   ------------------------------------------------------------
  Else If (npara_clf==11) Then
!   preliminary parameterization of clfac starting from unity (at NS)
!   over clin (at rin), clmid(rmid), clout(rout) reaching clfar (at rfar).

!   this is the linear formulation with respect to the histogram solution
!   provided by Puls et al. (2006) (see npara_clf = 9 version below)

!   all radii in units of nominal radius

    rphot = clf_field(2)

!   PHOTOSPHERIC CLUMPING NOT ALLOWED
    If (rphot<rscal(ns)) Then
      Write (999, *) ' STOP: PHOTOSPHERIC CLUMPING NOT ALLOWED(1)!!!'
      Stop ' PHOTOSPHERIC CLUMPING NOT ALLOWED(1)!!!'
    End If

    rin = clf_field(3)
    rmid = clf_field(4)
    rout = clf_field(5)
    rfar = clf_field(6)


    clfphot = clf_field(7)
    If (clfphot/=1.) Then
      Write (999, *) ' STOP: PHOTOSPHERIC CLUMPING NOT ALLOWED(2)!!!'
      Stop ' PHOTOSPHERIC CLUMPING NOT ALLOWED(2)!!!'
    End If

    clfin = clf_field(8)
    clfmid = clf_field(9)
    clfout = clf_field(10)
    clffar = clf_field(11)

    Do i = 1, nd
      If (rscal(i)<rfar) Exit
    End Do
    imax = max0(1, i-1)

    Do i = imax, nd
      If (rscal(i)<rout) Exit
    End Do
    iout = i - 1

    Do i = iout, nd
      If (rscal(i)<rmid) Exit
    End Do
    imid = i - 1

    Do i = imid, nd
      If (rscal(i)<rin) Exit
    End Do
    iin = i - 1

    Do i = iin, nd
      If (rscal(i)<rphot) Exit
    End Do
    iphot = i - 1

    istart = ns
!   PHOTOSPHERIC CLUMPING NOT ALLOWED
    If (iphot>=istart) iin = istart - 1

    clfac(1:imax) = clffar

    Do i = imax + 1, iout
!     linear interpol.
      clfac(i) = clffar + (clfout-clffar)/(rout-rfar)*(rscal(i)-rfar)
    End Do

    Do i = iout + 1, imid
!     linear interpol.
      clfac(i) = clfout + (clfmid-clfout)/(rmid-rout)*(rscal(i)-rout)
    End Do

    Do i = imid + 1, iin
!     linear interpol.
      clfac(i) = clfmid + (clfin-clfmid)/(rin-rmid)*(rscal(i)-rmid)
    End Do

    Do i = iin + 1, iphot
!     linear interpol.
      clfac(i) = clfin + (clfphot-clfin)/(rphot-rin)*(rscal(i)-rin)
    End Do

    Do i = iphot + 1, istart - 1
!     linear increase in case clfphot ne 1
      clfac(i) = 1. + (clfphot-1.)/(rphot-rscal(ns))*(rscal(i)-rscal(ns))
    End Do

    clfac(istart:nd) = 1.D0

!   for tests
    Print *
    Print *, ' CLUMPING FACTOR'
    Do i = 1, nd
      Print *, i, rscal(i), clfac(i)
    End Do

!   ------------------------------------------------------------
  Else If (npara_clf==9) Then
!   parameterization according to Puls et al. (2006)
!   rin, rmid, rout and rfar in units of nominal radius

    rin = clf_field(2)

!   PHOTOSPHERIC CLUMPING NOT ALLOWED
    If (rin<rscal(ns)) Then
      Write (999, *) ' STOP: PHOTOSPHERIC CLUMPING NOT ALLOWED!!!'
      Stop ' PHOTOSPHERIC CLUMPING NOT ALLOWED!!!'
    End If

    rmid = clf_field(3)
    rout = clf_field(4)
    rfar = clf_field(5)

    clfin = clf_field(6)
    clfmid = clf_field(7)
    clfout = clf_field(8)
    clffar = clf_field(9)

    clfac = 1.D0

    Do i = 1, nd

      If (rscal(i)>=rfar) Then
        clfac(i) = clffar

      Else If (rscal(i)<rfar .And. rscal(i)>=rout) Then
        clfac(i) = clfout

      Else If (rscal(i)<rout .And. rscal(i)>=rmid) Then
        clfac(i) = clfmid

      Else If (rscal(i)<rmid .And. rscal(i)>=rin) Then
        clfac(i) = clfin
      End If

    End Do

!   for tests
!   PRINT*
!   PRINT*,' CLUMPING FACTOR'
!   DO I=1,ND
!   PRINT*,I,RSCAL(I),CLFAC(I)
!   ENDDO

!   ------------------------------------------------------------
  Else
    Write (999, *) ' NUMBER OF CLUMPING PARAMETERS = ', npara_clf
    Write (999, *) &
      ' STOP: CORRESPONDING IMPLEMENTATION OF CLFAC NOT AVAILABLE'
    Print *, ' NUMBER OF CLUMPING PARAMETERS = ', npara_clf
    Stop ' CORRESPONDING IMPLEMENTATION OF CLFAC NOT AVAILABLE'
  End If
! ------------------------------------------------------------

! at this stage, clf-field (both for optically thin and optically
! thin clumping) has been processed.
! now we turn to optically thick clumping alone, and process
! fic_field, fvel_field, and hpor_field

  If (optthick) Then

!   POSSIBLE TO INCLUDE SEVERAL OPTION FOR THICK CLUMPING AS WELL
!   NOTE: MUST BE CONSISTENT WITH CORRESPNDING OPTION FOR THIN CLUMPING!
!   (I.E. CLFAC IS HERE *NOT* RECALCULATED)

    If (npara_clf==9) Then
!     Puls+ 2006 paramterization

!     FIC first
      rin = fic_field(2)
      rmid = fic_field(3)
      rout = fic_field(4)
      rfar = fic_field(5)

      clfin = fic_field(6)
      clfmid = fic_field(7)
      clfout = fic_field(8)
      clffar = fic_field(9)
      fic = 1.0D0
      fvol = 1.0D0
      Do i = 1, nd


!       bug created by JON, cured March 2018
!       (1-fic^2) -> (1-fic)^2
        If (rscal(i)>=rfar) Then
          fic(i) = clffar
          fvol(i) = (1.-fic(i))**2/(clfac(i)-2.*fic(i)+fic(i)**2)
        Else If (rscal(i)<rfar .And. rscal(i)>=rout) Then
          fic(i) = clfout
          fvol(i) = (1.-fic(i))**2/(clfac(i)-2.*fic(i)+fic(i)**2)
        Else If (rscal(i)<rout .And. rscal(i)>=rmid) Then
          fic(i) = clfmid
          fvol(i) = (1.-fic(i))**2/(clfac(i)-2.*fic(i)+fic(i)**2)
        Else If (rscal(i)<rmid .And. rscal(i)>=rin) Then
          fic(i) = clfin
          fvol(i) = (1.-fic(i))**2/(clfac(i)-2.*fic(i)+fic(i)**2)
        End If
!       volume filling factor, just as output for viewing
        If (fvol(i)>1.0) Then
          Write (999, *) ' STOP: FVOL > 1 -- Check FCL and FIC input'
          Stop ' FVOL > 1 -- Check FCL and FIC input'
        End If
      End Do

!     TCL FOR lines (multiplication factor) next
      rin = fvel_field(2)
      rmid = fvel_field(3)
      rout = fvel_field(4)
      rfar = fvel_field(5)
!     FVEL
      clfin = fvel_field(6)
      clfmid = fvel_field(7)
      clfout = fvel_field(8)
      clffar = fvel_field(9)
!     NOTE: JS -- CHANGED INPUT TO NORMALIZED VELOCITY FILLING FRACTION, cf
!     Sundqvist+ 2014
!     FVEL = 1.0D30
      fvel = 1.0
      tcl_fac_line = 0.0
      Do i = 1, nd

        If (rscal(i)>=rfar) Then
          fvel(i) = clffar
          If (fvel(i)>=1.0D0) Then
            fvel(i) = 1.0D0
            tcl_fac_line(i) = 0.0D0
          Else
            fvor = fvel(i)/(1.-fvel(i))
            tcl_fac_line(i) = (1.-(1.-fvol(i))*fic(i))/(abs(dvdr(i))*fvor)*sr/ &
              vmax
          End If
        Else If (rscal(i)<rfar .And. rscal(i)>=rout) Then
          fvel(i) = clfout
          If (fvel(i)>=1.0D0) Then
            fvel(i) = 1.0D0
            tcl_fac_line(i) = 0.0D0
          Else
            fvor = fvel(i)/(1.-fvel(i))
            tcl_fac_line(i) = (1.-(1.-fvol(i))*fic(i))/(abs(dvdr(i))*fvor)*sr/ &
              vmax
          End If
        Else If (rscal(i)<rout .And. rscal(i)>=rmid) Then
          fvel(i) = clfmid
          If (fvel(i)>=1.0D0) Then
            fvel(i) = 1.0D0
            tcl_fac_line(i) = 0.0D0
          Else
            fvor = fvel(i)/(1.-fvel(i))
            tcl_fac_line(i) = (1.-(1.-fvol(i))*fic(i))/(abs(dvdr(i))*fvor)*sr/ &
              vmax
          End If
        Else If (rscal(i)<rmid .And. rscal(i)>=rin) Then
          fvel(i) = clfin
          If (fvel(i)>=1.0D0) Then
            fvel(i) = 1.0D0
            tcl_fac_line(i) = 0.0D0
          Else
            fvor = fvel(i)/(1.-fvel(i))
            tcl_fac_line(i) = (1.-(1.-fvol(i))*fic(i))/(abs(dvdr(i))*fvor)*sr/ &
              vmax
          End If
        End If


!       Note: this includes a factor SR/VMAX!
!       Also we don't do it in one wash-up, since default does not become
!       tcl=0. then
      End Do

!     HPOR LAST
      rin = hpor_field(2)
      rmid = hpor_field(3)
      rout = hpor_field(4)
      rfar = hpor_field(5)
!     HPOR
      clfin = hpor_field(6)
      clfmid = hpor_field(7)
      clfout = hpor_field(8)
      clffar = hpor_field(9)

      hpor = 0.0
      tcl_fac_cont = 0.0

      Do i = 1, nd

        If (rscal(i)>=rfar) Then
          hpor(i) = clffar

        Else If (rscal(i)<rfar .And. rscal(i)>=rout) Then
          hpor(i) = clfout

        Else If (rscal(i)<rout .And. rscal(i)>=rmid) Then
          hpor(i) = clfmid

        Else If (rscal(i)<rmid .And. rscal(i)>=rin) Then
          hpor(i) = clfin
        End If
        tcl_fac_cont(i) = hpor(i)*v(i)*sr
!       Jo times SR or times SRNOM?
        tcl_fac_cont(i) = tcl_fac_cont(i)*(1.-(1.-fvol(i))*fic(i))

!       JS-CHANGE: Jun -16, now included non-void ic correction properly
!       Assume velocity stretch porosity law as default
!       Include scale-factor SR (gives porosity-length h 'real' units)
!       and also results in consistent tau_cl
      End Do

    Else If (npara_clf==3) Then

!     'standard' simple option, linearly increase to max/min value from onset
!     in v/vinf
!     NOTE: Is good e.g. in B/A-supergiant cases where velocity very low
!     beyond tau=2/3,
!     and will also be useful for consideration of optically thick winds

!     FIC FIRST-----------------------------------------------
      clf = fic_field(1)
      vclstart = fic_field(2)
      vclmax = fic_field(3)
      If (vclstart>vclmax) Then
        Write (999, *) ' STOP: VCLSTART > VCLMAX!'
        Stop ' VCLSTART > VCLMAX!'
      End If
      Do i = 1, nd
        If (v(i)<vclmax) Exit
      End Do
      imax = max0(1, i-1)
      Do i = imax, nd
        If (v(i)<vclstart) Exit
      End Do
      istart = i - 1 !              works also for vclstart = 0. (completely
!     clumped atmosphere)
      If (istart<imax) Then
        Print *, istart, imax
        Write (999, *) ' STOP: ISTART < IMAX IN CLUMPING'
        Stop ' ISTART < IMAX IN CLUMPING'
      End If
!     photospheric clumping not allowed except for tests (OPT_PHOTCLUMP)
      If (.Not. opt_photclump) Then
        If (istart>=ns) Then
          Print *, istart, ns, v(ns)
          Write (999, *) ' STOP: ISTART < NS IN CLUMPING'
          Stop ' ISTART < NS IN CLUMPING'
        End If
      End If

!     value reached
      fic(1:imax) = clf
!     linear increase
      Do i = imax + 1, istart
        fic(i) = 1. + (clf-1.)/(vclmax-vclstart)*(v(i)-vclstart)
      End Do
!     Unclumped inner atmosphere (FIC=1.0)
      fic(istart+1:nd) = 1.D0
!     -----------------------------------------------------
!     Volume filling factor
      Do i = 1, nd
        If (clfac(i)<=1.0 .Or. fic(i)>=1.0) Then
          fvol(i) = 1.0
        Else
!         bug created by JON, cured Sept 2018
          fvol(i) = (1.-fic(i))**2/(clfac(i)-2.*fic(i)+fic(i)**2)
!         print*,fic(i),clfac(i),fvol(i)
          If (fvol(i)>1.0) Then
            Write (999, *) ' STOP: FVOL >1, subroutine CLUMPING'
            Stop ' FVOL >1, subroutine CLUMPING'
          End If
        End If
      End Do

!     TCL_FAC_LINE AND FVEL
!     NEXT-----------------------------------------------
      clf = fvel_field(1)
      vclstart = fvel_field(2)
      vclmax = fvel_field(3)
      If (vclstart>vclmax) Then
        Write (999, *) ' STOP: VCLSTART > VCLMAX!'
        Stop ' VCLSTART > VCLMAX!'
      End If
      Do i = 1, nd
        If (v(i)<vclmax) Exit
      End Do
      imax = max0(1, i-1)
      Do i = imax, nd
        If (v(i)<vclstart) Exit
      End Do
      istart = i - 1 !              works also for vclstart = 0. (completely
!     clumped atmosphere)
      If (istart<imax) Then
        Print *, istart, imax
        Write (999, *) ' STOP: ISTART < IMAX IN CLUMPING'
        Stop ' ISTART < IMAX IN CLUMPING'
      End If
!     photospheric clumping not allowed except for tests (OPT_PHOTCLUMP)
      If (.Not. opt_photclump) Then
        If (istart>=ns) Then
          Print *, istart, ns, v(ns)
          Write (999, *) ' STOP: ISTART < NS IN CLUMPING'
          Stop ' ISTART < NS IN CLUMPING'
        End If
      End If
!     value reached
      fvel(1:imax) = clf
!     linear increase
      Do i = imax + 1, istart
        fvel(i) = 1. + (clf-1.)/(vclmax-vclstart)*(v(i)-vclstart)
      End Do
!     Unclumped inner atmosphere (FIC=1.0)
      fvel(istart+1:nd) = 1.D0

      Do i = 1, nd
        If (clfac(i)<=1.0 .Or. fic(i)>=1.0) Then
          tcl_fac_line(i) = 0.0D0
        Else
          fvor = fvel(i)/(1.-fvel(i))
          tcl_fac_line(i) = (1.-(1.-fvol(i))*fic(i))/(abs(dvdr(i))*fvor)*sr/ &
            vmax
        End If
      End Do
!     -----------------------------------------------------

!     HPOR and TCL_FAC_CONT last---------------------------
      clf = hpor_field(1)
      vclstart = hpor_field(2)
      vclmax = hpor_field(3)
      If (vclstart>vclmax) Then
        Write (999, *) ' STOP: VCLSTART > VCLMAX!'
        Stop ' VCLSTART > VCLMAX!'
      End If
      Do i = 1, nd
        If (v(i)<vclmax) Exit
      End Do
      imax = max0(1, i-1)
      Do i = imax, nd
        If (v(i)<vclstart) Exit
      End Do
      istart = i - 1 !              works also for vclstart = 0. (completely
!     clumped atmosphere)
      If (istart<imax) Then
        Print *, istart, imax
        Write (999, *) ' STOP: ISTART < IMAX IN CLUMPING'
        Stop ' ISTART < IMAX IN CLUMPING'
      End If
!     photospheric clumping not allowed except for tests (OPT_PHOTCLUMP)
      If (.Not. opt_photclump) Then
        If (istart>=ns) Then
          Print *, istart, ns, v(ns)
          Write (999, *) ' STOP: ISTART < NS IN CLUMPING'
          Stop ' ISTART < NS IN CLUMPING'
        End If
      End If
!     value reached
      hpor(1:imax) = clf
!     Unclumped inner atmosphere (FIC=1.0)
      hpor(istart+1:nd) = 0.D0
!     linear increase between
      Do i = imax + 1, istart
        hpor(i) = hpor(istart+1) + (hpor(imax)-hpor(istart+1))/(vclmax- &
          vclstart)*(v(i)-vclstart)
      End Do

      Do i = 1, nd
        If (clfac(i)<=1.0 .Or. fic(i)>=1.0) Then
          tcl_fac_cont(i) = 0.0D0
        Else
!         Jo times SR or times SRNOM?
          tcl_fac_cont(i) = hpor(i)*v(i)*sr
          tcl_fac_cont(i) = tcl_fac_cont(i)*(1.-(1.-fvol(i))*fic(i))
        End If
      End Do
!     -----------------------------------------------------

    Else
      Write (999, *) ' STOP: THIS OPTHICK CLUMPING-OPTION NOT ALLOWED!, &
        &subroutine clumping'
      Stop ' THIS OPTHICK CLUMPING-OPTION NOT ALLOWED!, subroutine clumping'
    End If

  Else

    Do i = 1, nd
      fic(i) = 1.0D0
      fvel(i) = 1.0D0
      hpor(i) = 0.0D0
      fvol(i) = 1.0D0
      tcl_fac_line(i) = 0.0D0
      tcl_fac_cont(i) = 0.0D0
!     Switch off by setting all clump optical depths to zero
    End Do

  End If

  Open (13, File=trim(modnam)//'/CLUMPING_OUTPUT', Status='UNKNOWN')
! WRITE(13,fmt='(A)') '            r/Rstar          velo      abs(dvdr)
! f_cl       f_ic           f_vel              hpor           f_vol
! tcl_fac_line             tcl_fac_cont'
  Write (13, Fmt='(A)') '            r/Rstar          velo     &
    &    abs(dvdr)         f_cl       f_ic           f_vel     &
    &         hpor           f_vol           tcl_fac_line      &
    &       tcl_fac_cont'
  Do i = 1, nd
    Write (13, 100) i, rscal(i), v(i), abs(dvdr(i)), clfac(i), fic(i), &
      fvel(i), hpor(i), fvol(i), tcl_fac_line(i), tcl_fac_cont(i)
  End Do
  Close (13)

! stop ' testing'
! ALWAYS!!!
! check for representative value of clf
  If (clf_representative<minval(clfac) .Or. clf_representative>maxval(clfac)) &
    Then
    Print *, ' REPRESENTATIVE VALUE OF CLUMPING FACTOR OUTSIDE ACTUAL RUN'
    Print *, clf_representative
    Do i = 1, nd
      Print *, i, r(i), v(i), clfac(i)
    End Do
    Write (999, *) &
      ' STOP: REPRESENTATIVE VALUE OF CLUMPING FACTOR OUTSIDE ACTUAL RUN'
    Stop ' REPRESENTATIVE VALUE OF CLUMPING FACTOR OUTSIDE ACTUAL RUN'
  End If

! finally, calculate overdensity
  rhocl = rho*clfac

  Return
100 Format (I3, 2X, 4F15.5, 2E15.5, 2F15.5, 2E25.10)

End Subroutine

!-----------------------------------------------------------------------

Subroutine derivs(x, y, dy)

  Use :: nlte_type
  Use :: nlte_dim
  Use :: nlte_opt, Only: opt_photclump
  Use :: nlte_var, Only: srnom, idum, ggrav, teff, sigmae => sigem, ckappa, &
    xt, a2te, delsig, delmu, dqdtau, delpara, itmin

  Use :: photstruc, Only: clf_const

  Implicit None

! delsig corrects for constant sigmae
! delmu  corrects for constant mue
! deli   corrects for recombination

! .. scalar arguments ..
  Real (dp) :: x
! ..
! .. array arguments ..
  Real (dp) :: dy(3), y(3)
! ..
! .. local scalars ..
  Real (dp) :: delm, dels, dy1, dy2, term0, term1, term2, deli
! ..

  dels = .5D0*(delsig(idum-1)+delsig(idum))
  delm = .5D0*(delmu(idum-1)+delmu(idum))
  deli = .5D0*(delpara(idum-1)+delpara(idum))

  term0 = sigmae + ckappa*deli/delm*y(1)*y(2)**(-xt-1.D0)

! note that ckappa includes sigmae
  term1 = term0*dels

  dy1 = ggrav - 5.67D-5/2.9979D10*teff**4*term1
! this is the ONLY hack to ensure that rho_phot => rho_phot/fcl
  If (opt_photclump) dy1 = dy1/clf_const

! with this statement we can check for grad
! if(idum.ge.38) print*,idum,5.67D-5/2.9979D10*TEFF**4*TERM1
  dy(1) = dy1*(srnom/y(3))**2

  term2 = dels/y(2)**3*term0

  dy2 = 3.D0/16.D0*teff**4*term2*(1.D0+dqdtau(idum))
  dy(2) = dy2*(srnom/y(3))**2

! THIS IS NEW
  If (idum<=itmin) dy(2) = 0.D0

  dy(3) = -y(2)/y(1)*a2te*delm

  Return

End Subroutine

!-----------------------------------------------------------------------

Subroutine derivs_kur(x, y, dy)

  Use :: nlte_type
  Use :: nlte_dim
  Use :: nlte_var, Only: idum, a2te, delmu, itmin
  Use :: photstruc

  Implicit None

! delmu  corrects for constant mue

! .. scalar arguments ..
  Real (dp) :: x
! ..
! .. array arguments ..
  Real (dp) :: dy(3), y(3)
! ..
! .. local scalars ..
  Real (dp) :: delm, dpdm, dtdm
! ..
  Integer (i4b) :: i

  Do i = 1, ndmod - 1
    If (x>xmext(i) .And. x<xmext(i+1)) Exit
  End Do
! re-changed Nov. 2022
  If (i==ndmod) Then
    If (x<xmext(1)) Then
      Print *, 'xm not found at outer boundary of external struct. ', x, &
        xmext(1)
      i = 1
    Else
      Print *, x, xmext(1), xmext(ndmod)
      If (x<xmext(1)) Print *, ' INCREASE MDOT!!!'
      Write (999, *) ' STOP: XM NOT FOUND IN EXTERNAL PHOT. STRUCT.'
      Stop ' XM NOT FOUND IN EXTERNAL PHOT. STRUCT.'
    End If
  End If

  dpdm = (pkur(i+1)-pkur(i))/(xmext(i+1)-xmext(i))
  dtdm = (text(i+1)-text(i))/(xmext(i+1)-xmext(i))

  delm = .5D0*(delmu(idum-1)+delmu(idum))


  dy(1) = dpdm


  dy(2) = dtdm
! THIS IS NEW
  If (idum<=itmin) dy(2) = 0.D0

  dy(3) = -y(2)/y(1)*a2te*delm

  Return

End Subroutine

!-----------------------------------------------------------------------

Subroutine derivs_tlu(x, y, dy)

  Use :: nlte_type
  Use :: nlte_dim
  Use :: nlte_var, Only: idum, itmin
  Use :: photstruc

  Implicit None

! delmu  corrects for constant mue

! .. scalar arguments ..
  Real (dp) :: x
! ..
! .. array arguments ..
  Real (dp) :: dy(3), y(3)
! ..
! .. local scalars ..
  Real (dp) :: drhodm, dtdm
! ..
  Integer (i4b) :: i

  Do i = 1, ndmod - 1
    If (x>xmext(i) .And. x<xmext(i+1)) Exit
  End Do
  If (i==ndmod) Then
!   STOP ' xm not found in external phot. struct.'
!   extraplolation at own risk
    If (x>xmext(ndmod)) Then
      i = ndmod - 1
    Else If (x<xmext(1)) Then
      i = 1
    Else
      Write (999, *) &
        ' STOP: SOMETHING WRONG WITH XM IN EXTERNAL PHOT. STRUCT.'
      Stop ' SOMETHING WRONG WITH XM IN EXTERNAL PHOT. STRUCT.'
    End If
  End If

  drhodm = (rhotlu(i+1)-rhotlu(i))/(xmext(i+1)-xmext(i))
  dtdm = (text(i+1)-text(i))/(xmext(i+1)-xmext(i))

! no longer used  delm = .5D0*(delmu(idum-1)+delmu(idum))

  dy(1) = drhodm

  dy(2) = dtdm
! THIS IS NEW
  If (idum<=itmin) dy(2) = 0.D0

  dy(3) = -1.D0/y(1)

  Return

End Subroutine

!-----------------------------------------------------------------------

Function rhofunc(r, v)

  Use :: nlte_type
  Use :: nlte_dim
  Use :: nlte_var, Only: h1
  Implicit None

! berechnet die dichte aus mloss und v

! .. scalar arguments ..
  Real (dp) :: r, v, rhofunc
! ..

  rhofunc = h1/(r**2*v)

  Return

End Function

!-----------------------------------------------------------------------

Function vfvl(r, rc, r0, sr, vinf, vmin, beta, hl)

  Use :: nlte_type
  Use :: nlte_dim
  Implicit None

! calculates velocities for analytical vel.law (hamann, 1985)

! .. scalar arguments ..
  Real (dp) :: beta, hl, r, r0, rc, sr, vinf, vmin, vfvl
! ..
! .. intrinsic functions ..
! INTRINSIC EXP
! ..

  If (r<rc) Then
    vfvl = vmin*exp((r-sr)/hl)
  Else
    vfvl = vinf*(1.D0-r0/r)**beta
  End If

  Return

End Function

!-----------------------------------------------------------------------

Subroutine vindex(r, v, rho, m, i)

  Use :: nlte_type
  Use :: nlte_dim
  Implicit None

! verschiebt m indices ab index i um 1 nach oben

! .. parameters ..
  Integer (i4b), Parameter :: nd1 = id_ndept
! ..
! .. scalar arguments ..
  Integer (i4b) :: i, m
! ..
! .. array arguments ..
  Real (dp) :: r(nd1), rho(nd1), v(nd1)
! ..
! .. local scalars ..
  Integer (i4b) :: j, k
! ..

  Do j = 0, m
    k = m + i + 1 - j
    r(k) = r(k-1)
    v(k) = v(k-1)
    rho(k) = rho(k-1)
  End Do

  Return

End Subroutine

!***********************************************************************

!subroutines: complex ones
!rateeq and related

!***********************************************************************

Subroutine rateeq(xnh, xne, temp, clf, ilow, imax, nd, err, concon, r, v,  &
  errold, specmat, accel, optmet, meanerr)

! UPATED FOR CLUMPING
! ILOW AND IMAX ARE EITHER LTE OR NLTE VALUES, DEPENDENT ON CALL

  Use :: nlte_type
  Use :: nlte_dim

  Use :: princesa_var, Only: nat, ifirsl, labat, labl, le

  Use :: nlte_var, Only: alevel, blevel, enionnd, enionnd_lte, iqual, modnam, &
    almost_converged, optcmf_all

  Use :: nlte_opt, Only: opt_damp_fg

  Use :: tcorr_var, Only: qcbfr, qcbfi, qrbfr, qrbfi, dtqrbfr, dtqrbfi, &
    dtqcbfh, dtqcbfc, qcbbu, qcbbd, dtqcbbu, dtqcbbd, opttcor

  Implicit None

! .. parameters ..
  Integer (i4b), Parameter :: kel = id_atoms, kis = id_kisat
  Integer (i4b), Parameter :: nrec = id_llevs
  Integer (i4b), Parameter :: nd1 = id_ndept

! ..
! .. scalar arguments ..
  Real (dp) :: meanerr
  Integer (i4b) :: nd
  Logical :: accel, concon, optmet
! ..
! .. array arguments ..
  Real (dp) :: err(nd1), errold(nd1), r(nd1), specmat(5, nd1), &
    temp(nd1), v(nd1), xne(nd1), xnh(nd1), clf(nd1)

  Integer (i4b) :: ilow(nd1, kel), imax(nd1, kel)
! ..
! .. local scalars ..
  Real (dp) :: ampl, binf, en, errn, occ, occold, occoldmax, specmax, specmin, &
    specnew, velr, xmust
  Integer (i4b) :: i, i1, i2, iilo, iima, j, k, ll, numion, nato, jerr1
! ..
! .. local arrays ..
  Real (dp), Dimension (2) :: err1
  Integer (i4b), Dimension (2) :: jerr

  Real (dp) :: blev0(id_llevs)

  Real (dp), Dimension (nrec, 2) :: qrion, qrrec, qcion, qcrec, qlu, qld

  Real (dp) :: summ, summ1, bold, error, error1
  Logical :: damped

! ..
! .. external functions ..
  Integer (i4b) :: igenio
  External :: igenio
! ..
! .. external subroutines ..
  External :: factinc, netmat, update
! ..
! .. intrinsic functions ..
! INTRINSIC ABS,MAX,MIN,SQRT
! ..

depthloop: Do ll = 1, nd

    err1 = 0.D0
    Read (17, Rec=ll)(blevel(i), i=1, nrec)
    Read (15, Rec=ll)(alevel(i), i=1, nrec)

    Do i = 1, nrec

!     !!! value from previous iteration
      nato = le(i)

      blev0(i) = blevel(i)

!     !!! set to zero to obtain only actually calculated values in
!     dependence of ilow,imax
!     ONLY FOR H,He if pure continuum

      If (.Not. concon) Then
        If (labat(nato)=='H' .Or. labat(nato)=='HE') Then
          blevel(i) = 0.D0
        Else
!         ---- use UPDATED LTE value for metals!!!
!         otherwise, blevel and enionnd reamain at their LTE value from it. 0,
!         and the boundary condition is deteriorated (=> opacont must be
!         consistent)
          blevel(i) = alevel(i)
        End If
      End If
    End Do

!   ---- calculation of limb-darkening coefficients

    If (concon) Then
      Call factinc(r(ll), ll)
!     THIS IS NEW FROM V9.5 ON
      blevel = 0.D0
    End If

    velr = v(ll)/r(ll)
    xmust = sqrt(1.D0-1.D0/r(ll)/r(ll))

    If (opttcor) Then
      qrion = 0.
      qrrec = 0.
      qcion = 0.
      qcrec = 0.
      If (concon) Then
        qlu = 0.
        qld = 0.
      End If
    End If

kloop: Do k = 1, nat
      iilo = ilow(ll, k)
      iima = imax(ll, k)
      If (iima<iilo) Then
        Write (999, *) ' STOP: IIMA < IILO'
        Stop ' IIMA < IILO'
      End If
      If (.Not. concon .And. labat(k)/='H' .And. labat(k)/='HE') Then
!       -----use UPDATED enionnd for metals (see above)
        enionnd(k, :, :) = enionnd_lte(k, :, :)
        Go To 110
      End If

!     -----calculation of net bound-free rates and derivatives
!     -----setup of nlte equations and newton-raphson iteration

      Call netmat(k, iilo, iima, ll, xne(ll), temp(ll), clf(ll), blev0, &
        concon, errold(ll), optmet, qrion, &
        qrrec, qcion, qcrec, qlu, qld)

!     for tests
!     IF(LL.EQ.9) THEN
!     DO I = IILO, IIMA + 1
!     NUMION = IGENIO(K,I)
!     I1 = IFIRSL(NUMION)
!     I2 = IFIRSL(NUMION+1) - 1
!     IF(I.EQ.IIMA+1) I2=I1

!     DO J=I1,I2
!     BOLD=BLEV0(J)
!     PRINT*,K,I,J,BOLD,BLEVEL(J)
!     ENDDO
!     ENDDO
!     ENDIF

      summ = sum(blevel)
!     JO changed Oct. 2016
!     tested again April 2017: this damping at -2 better than no damping
!     JO comment from June 2018: this needs to be kept, otherwise quite
!     strong oscillations are possible
      If ((almost_converged .And. opt_damp_fg) .Or. (optcmf_all .And. meanerr< &
        -2. .And. labat(k)/='H' .And. labat(k)/='HE')) Then
!       LABAT(K).eq.'XX') ) THEN
!       damp changes. Typically, these are levels which oscillate due
!       to numerical problems of bb-rates, or due to overlapping lines.
!       Could be cured by more depth points.
        damped = .False.
        Do i = iilo, iima + 1
          numion = igenio(k, i)
          i1 = ifirsl(numion)
          i2 = ifirsl(numion+1) - 1
          If (i==iima+1) i2 = i1

          Do j = i1, i2
            bold = blev0(j)
            If (bold==0.D0) Then
              Print *, ll, k, i, j
              Write (999, *) ' STOP:  BOLD = 0 IN RATEEQ'
              Stop ' BOLD = 0 IN RATEEQ'
            End If
            error1 = 1. - bold/blevel(j)
            error = abs(error1)
            If (error>0.05) Then !  log = -1.3
!             large changes damped to 2.5%
!             IF(ERROR1.GT.0) THEN
!             BLEVEL(J)=BOLD/0.975
!             ELSE
!             BLEVEL(J)=BOLD/1.025
!             ENDIF
!             JO Oct. 2025: take geometric mean
              blevel(j) = sqrt(blevel(j)*bold)
              Write (*, 120) ll, k, i, j, error1
              damped = .True.

            Else If (error>0.025) Then
!             smaller changes damped to 1.5% (note: convergence at 1%)
!             IF(ERROR1.GT.0) THEN
!             BLEVEL(J)=BOLD/0.985
!             ELSE
!             BLEVEL(J)=BOLD/1.015
!             ENDIF
              blevel(j) = sqrt(blevel(j)*bold)
!             WRITE(*,15) LL,K,J,I,ERROR1
!             15              FORMAT(' L = ',I3,': UPDATE OF LEVEL
!             ',3(2X,I3),2X,F10.5,' DAMPED')
              damped = .True.
            End If

          End Do
        End Do
!       renormalize level occupation
        If (damped) Then
          summ1 = sum(blevel)
          error = summ/summ1
          If (error>1.05 .Or. error<0.95) Then
            Print *, error
            Write (999, *) &
              ' STOP: TOO LARGE RENORMALIZATION OF DAMPED CORRECTIONS'
            Stop ' TOO LARGE RENORMALIZATION OF DAMPED CORRECTIONS'
          End If
          blevel = blevel*error
!         PRINT*,' LEVEL OCCUP. RENORMALIZED BY FACTOR',ERROR
        End If

      End If !                      damping

!     for tests, to keep occup const.
!     blevel=blev0


!     -----calculation of enionnd(k,i,ll)

iloop: Do i = iilo, iima
        numion = igenio(k, i)
        i1 = ifirsl(numion)
        i2 = ifirsl(numion+1) - 1
        en = 0.D0

!       ---  maxmimum occ.numbers (old)

        occoldmax = 0.D0

        Do j = i1, i2 + 1
          occoldmax = max(occoldmax, blev0(j))
        End Do

        occoldmax = 1.D0/occoldmax

jloop:  Do j = i1, i2

!         check for negative occup. no.

          If (blevel(j)<=0.D0) Then
            Print *, ll, ' WARNING: NEGATIVE OCCUP. NO AT LEVEL ', j
            Write (999, *) ll, ' WARNING: NEGATIVE OCCUP. NO AT LEVEL ', j
            blevel(j) = blev0(j)
          End If

          occ = blevel(j)
          occold = blev0(j)
!         ---
!         ---     only those errors are considered, which are not contaminated
!         ---     by insufficient numerics
!         ---
          If (occold*occoldmax<1.D-11) Go To 100
          If (iqual(j,ll)/=0) Go To 100
          errn = abs(1.D0-occ/occold)
          If (errn>err1(1)) Then
            err1(1) = errn
            jerr(1) = j
          Else If (errn>err1(2)) Then
            err1(2) = errn
            jerr(2) = j
          End If

100       Continue

          en = en + blevel(j)

        End Do jloop
        enionnd(k, i, ll) = en
      End Do iloop

      i1 = ifirsl(numion+1)

      occ = blevel(i1)

      If (blevel(i1)<=0.D0) Then
        Print *, ll, ' WARNING: NEGATIVE OCCUP. NO AT IMAX ', i1
        Write (999, *) ll, ' WARNING: NEGATIVE OCCUP. NO AT IMAX ', i1
        blevel(i1) = blev0(i1)
      End If

      occold = blev0(i1)

      If (occold*occoldmax>1.D-11 .And. iqual(i1,ll)==0) Then
        errn = abs(1.D0-occ/occold)
        If (errn>err1(1)) Then
          err1(1) = errn
          jerr(1) = i1
        Else If (errn>err1(2)) Then
          err1(2) = errn
          jerr(2) = i1
        End If
      End If

      jerr1 = jerr(1)
      If (jerr1/=i1 .And. accel .And. .Not. almost_converged) Then

!       extrapolation of occupation number

        If (errold(ll)>6.D-3) Then
          specnew = err1(1)/errold(ll)
        Else
          specnew = 0.D0
        End If

        specmax = max(specmat(1,ll), specmat(2,ll), specmat(3,ll), &
          specmat(4,ll), specmat(5,ll), specnew)
        specmin = min(specmat(1,ll), specmat(2,ll), specmat(3,ll), &
          specmat(4,ll), specmat(5,ll), specnew)
        If (specmax<1.D0 .And. specmin>.95D0) Then
          ampl = 1.D0/(1.D0-specmin)
          Print *
          Print *, 'EXTRAPOLATION OF LEVEL ', jerr1
          Print *, 'EXTRAPOLATION OF LEVEL ', jerr1
          Print *, 'EXTRAPOLATION OF LEVEL ', jerr1
          Print *
          binf = blev0(jerr1) + (blevel(jerr1)-blev0(jerr1))*ampl

          If (binf<.1D0*blevel(jerr1)) binf = .1D0*blevel(jerr1)

          If (binf>10.D0*blevel(jerr1)) binf = 10.D0*blevel(jerr1)

          blevel(jerr1) = binf
          Write (11, Fmt=*) 'EXTRAPOLATED VALUE', binf
          err1(1) = abs(1.D0-binf/blev0(jerr1))
        End If
      End If

!     ---  last ionization stage assumed always as single level

      enionnd(k, iima+1, ll) = blevel(i1)

!     not calculated ions explicitly set to zero

110   Do i = 1, iilo - 1
        enionnd(k, i, ll) = 0.D0
      End Do

      Do i = iima + 2, kis + 1
        enionnd(k, i, ll) = 0.D0
      End Do

    End Do kloop

!   update heating/cooling rates: at this place, since consistent
!   nlte numbers have to be used
!   (rate-coefficients from actual lte numbers!)

!   NOTE: since only those coefficients are populated which have
!   transitions within (imin,imax), we can multiply with TOTAL vector blevel
    If (opttcor) Then
      qrbfr(ll) = dot_product(qrrec(:,1), blevel)
      dtqrbfr(ll) = dot_product(qrrec(:,2), blevel)
      qrbfi(ll) = dot_product(qrion(:,1), blevel)
      dtqrbfi(ll) = dot_product(qrion(:,2), blevel)
      qcbfr(ll) = dot_product(qcrec(:,1), blevel)
      dtqcbfh(ll) = dot_product(qcrec(:,2), blevel)
      qcbfi(ll) = dot_product(qcion(:,1), blevel)
      dtqcbfc(ll) = dot_product(qcion(:,2), blevel)
      If (concon) Then
        qcbbu(ll) = dot_product(qlu(:,1), blevel)
        dtqcbbu(ll) = dot_product(qlu(:,2), blevel)
        qcbbd(ll) = dot_product(qld(:,1), blevel)
        dtqcbbd(ll) = dot_product(qld(:,2), blevel)
      End If
    End If

    Call update(xnh, xne, temp, ilow, imax, ll, nat)

    err(ll) = err1(1)
!   nnn=ll/6
!   nnuu=nnn*6
!   if(ll.eq.1.or.ll.eq.47) nnuu=ll
!   if(nnuu.eq.ll) then

    Print *
    Print *, ' MAXIMUM CHANGE IN OCCUPATION NUMBERS = ', err(ll), ' ', &
      labl(jerr(1))
    Write (11, Fmt=*) ll, ' MAX. CHANGE IN OCC. NUMBERS = ', err(ll), ' ', &
      labl(jerr(1))
    Write (11, Fmt=*) ll, ' 2ND  CHANGE IN OCC. NUMBERS = ', err1(2), ' ', &
      labl(jerr(2))
    Print *
    Print *

!   endif

    Write (16, Rec=ll)(blevel(i), i=1, nrec)
    Write (17, Rec=ll)(blevel(i), i=1, nrec)

  End Do depthloop

  Open (1, File=trim(modnam)//'/ENION', Form='UNFORMATTED', Status='UNKNOWN')
  Rewind 1
  Write (1) enionnd, xne
  Close (1)
! print*,'xnecheck: enion written'
! print*,xne
  Return
120 Format (' L = ', I3, ': UPDATE OF LEVEL ', 3(2X,I3), 2X, F10.5, ' DAMPED')

End Subroutine

!-----------------------------------------------------------------------

Subroutine factinc(rr, l)

  Use :: nlte_type
  Use :: nlte_dim
  Use :: nlte_var, Only: as, a2s, a3s, bs, b2s, b3s, ifre

  Implicit None

! below, find the description of the old approach no longer supported

! ---- calculation of limb darkening factors ('blocking factors') for
! ---- i_inc using momenta j,k,h,n
! i_inc = as +bs * mue    for mu_*<mu<1
! i_inc = a2s+b2s* mue    for  0 <mu<mu_*
! i_inc = a3s+b3s* mue    for  -1<mue<0
! ---- it has been changed (e. santolaya, october 92) so that a new
! ---- phylosophy is applied: j, h, k and n are now vectors wich elements
! ---- contain the radiation field at the different frecuency points
! ---- at a certain depth point.

! .. parameters ..
  Integer (i4b), Parameter :: ifretot = id_frec1
! ..
! .. scalar arguments ..
  Real (dp) :: rr
  Integer (i4b) :: l
! ..
! .. local scalars ..
  Real (dp) :: aa, bb, cc, dd, smue, smued, smueq, smuez, xxhn
  Integer (i4b) :: i, lll
! ..
! .. local arrays ..
  Real (dp) :: xxh(ifretot), xxj(ifretot), xxk(ifretot), xxn(ifretot)
! ..
! .. intrinsic functions ..
! INTRINSIC SQRT
! ..

  Read (40, Rec=l)(xxj(i), i=1, ifre)
  Read (42, Rec=l)(xxk(i), i=1, ifre)
  Read (41, Rec=l)(xxh(i), i=1, ifre)
  Read (43, Rec=l)(xxn(i), i=1, ifre)

! This was the old approach

  smue = sqrt(1.D0-1.D0/rr/rr)

! ---- smueq=(1-(mu_*)**2)**2

  smueq = 1.D0/rr**4
  smued = 3.D0 - smue**2
  smuez = smue**2 - 2.D0

  Do lll = 1, ifre

    xxhn = 24.D0*xxh(lll) - 40.D0*xxn(lll)

    If (smueq/=0.D0) Then
      bb = xxhn/smueq
    Else
      bb = xxhn*rr**4
    End If

    cc = 6.D0*xxh(lll) - bb*.5D0*smued
    dd = bb*smue*smuez + 24.D0*xxk(lll) - 8.D0*xxj(lll)
    aa = 0.5D0*smue*bb*smued + 6.D0*xxj(lll) - 12.D0*xxk(lll)

    as(lll) = (aa+bb)/2.D0
    a2s(lll) = (aa-bb)/2.D0
    a3s(lll) = a2s(lll)
    bs(lll) = (cc+dd)/2.D0

    If (smue/=0.D0) Then
      b2s(lll) = (as(lll)-a2s(lll))/smue + bs(lll)
    Else
      b2s(lll) = 0.D0
    End If

    b3s(lll) = (cc-dd)/2.D0

  End Do

  Return

End Subroutine

!-----------------------------------------------------------------------

Subroutine netmat(nato, ilo, ima, ll, xne, temp, clf, blev0, concon, &
     errold, optmet, qrion, qrrec, qcion, qcrec, qlu, qld)

! UPDATED FOR CLUMPING
! ADDITIONAL CYCLE FOR EXPLICIT DR TRANSITIONS INCLUDED:
! ASSUMES THAT PARENTAL LEVEL IS GROUNDSTATE OF NEXT ION (KL!!!)

  Use :: nlte_type
  Use :: nlte_dim
  Use :: fund_const
! Use :: princesa_var, Only: nl, nlx, ns, nat, ion, iong, la0, la1, ixa0, &
! ixa1, isa0, isa1, nions, ifirsl, ifirx, ifirs, kl, le, li, klx, lix, ns0, &
! ns1, lis, zeff, weight, abund, gl, fl, zl, glx, flx, gs0, gs1, ryd, qd, &
! drexplicit, lkl, labl, labat, labl3, labl4, frecin, frecfn, flcbf

  Use :: princesa_var, Only: le, li, iong, kl, fl, ifirsl, labat, zeff, zl, &
    labl, labl3, labl4, drexplicit, frecin, frecfn, flcbf


  Use :: nlte_var, Only: abundnd, xnelte, alevel, blevel, mlow, mup, mclow, &
    mcup, index4 => indr, index3 => indc, precis, unasol, megas, nist, xnkk, &
    fre, ifre, op_flag, frecfn1, indfrq1

  Use :: nlte_xrays, Only: optxray, optauger, name, n_kedges, k_nmin, k_nmax, &
    eth, z_k => z, n_k => n, aug_1_6

  Use :: tcorr_var, Only: opttcor

  Use :: nlte_opt, Only: optdr

  Implicit None

! -------this subroutine calculates the net rates, setups the rate
! -------equation and solves them iteratively in a newton-raphson scheme
! -------starting with the populations obtained by solving 'a pelo' the
! -------system.

! .. parameters ..
! Integer (i4b), Parameter :: kel = id_atoms, kis = id_kisat
  Integer (i4b), Parameter :: ndim = id_llevs
  Integer (i4b), Parameter :: nnrec = id_llevs
! ..
! .. scalar arguments ..
  Real (dp) :: errold, temp, xne, clf
  Integer (i4b) :: ilo, ima, ll, nato
  Logical :: concon, optmet
! ..
! .. array arguments ..
  Real (dp), Dimension (nnrec) :: blev0
  Real (dp), Dimension (nnrec, 2) :: qrion, qrrec, qcion, qcrec, qlu, qld
! ..
! ..
! .. local scalars ..
  Real (dp) :: bmax, casta1, casta2, cik, depart, det, errmax, x1, x2, xnerat, &
    xni, xnistar, xnk, xxlev, x3, x4, dtcik, x5, x6, hhnu, cik1, prob
  Integer (i4b) :: i, iii, itm, j, jmax, ml, mn, mu, mu1, n1, n2, nfir, nlas, &
    numion, k, l, itrans, jeff, nfm, nlm, jj, ij, lo
! ..
! .. local arrays ..
  Real (dp) :: bold(ndim), deltan(ndim), f(ndim), flud(ndim), &
    rat2(ndim, ndim), ratmat(ndim, ndim)

  Real (qp) :: aux(ndim)

  Integer (i4b) :: indlud(ndim)
! ..
! .. external functions ..
  Integer (i4b) :: igenio
  External :: igenio
! ..
! .. external subroutines ..
  External :: intermecol, intermepho, linesob, lubksb, ludcmp
! ..
! .. intrinsic functions ..
! INTRINSIC ABS,MAX
! ..

  itm = 10

  nist = 0
  xnerat = xne/xnelte(ll)
  numion = igenio(nato, ilo)
  nfir = ifirsl(numion)

  numion = igenio(nato, ima) + 1
  If (numion>iong) Then
    Write (999, *) ' STOP: ERROR IN NETMAT - NUMION'
    Stop ' ERROR IN NETMAT - NUMION'
  End If
  nlas = ifirsl(numion)

  mn = nlas - nfir + 1

! ------ initialitation of rate matrix and solution vector

  f(1:mn) = 0.D0
  ratmat(1:mn, 1:mn) = 0.D0

! ------ radiative b-f transitions

rbfloop: Do i = 1, index4

!   standard treatment, all rbf-data, integration in INTERMEPHO
!   starting at FRECIN (gs or excited state)

    ml = mlow(i)
    mu = mup(i)

    If (ml<nfir .Or. ml>=nlas) Cycle

    xni = blev0(ml)

    If (mu>nlas) Then
!     this path refers to ionization to excited levels of imax+1
!     upper level assumed to be in lte with respect to nlas (groundstate of
!     imax+1)
!     new formulation, see notes
      mu = nlas
    End If

    xnk = blev0(mu)
    xxlev = alevel(ml)/alevel(mu)*xnerat

    xnistar = xnk*xxlev !           n_u *(n_i/n_u)*
    nist(ml) = xxlev !              witout nk
    xnkk(ml) = xnk
    depart = xni/xnistar

!   -------   photo-rates (integration starting at INDFRQ(I),
!   corresponding to FRECIN)

    Call intermepho(i, depart, xnistar, ll, x1, x2, xne, temp, clf, optmet, &
      x3, x4, x5, x6)

!   all rates entering the energy equation are with respect to overdense
!   clumps
!   and according occupation numbers (NOT AVERAGED!)

    If (opttcor) Then
      hhnu = hh*fl(ml) !            to be consistent with LTE
      If (hhnu<=0.D0) Then
        Write (999, *) ' STOP: ERROR IN HNU - RBF'
        Stop ' ERROR IN HNU - RBF'
      End If

!     dtrbfc(n1) = XNK*XXLEV*X6 !does presently not work at low temp. Why?
      qrrec(mu, 1) = qrrec(mu, 1) + x4*xxlev ! (*n_k)
      qrion(ml, 1) = qrion(ml, 1) + x3 ! (*n_i)
      qrrec(mu, 2) = qrrec(mu, 2) + xxlev*(x6-x4*(1.5D0+hhnu/akb/temp)/temp)
!     (*n_k)
      qrion(ml, 2) = qrion(ml, 2) - x3*(1.5D0+hhnu/akb/temp)/temp ! (*n_i),
!     minus!
    End If

    n1 = ml - nfir + 1
    n2 = mu - nfir + 1

    If (n2>mn .Or. n1<0) Then
      Write (999, *) ' STOP: ERROR IN NETMAT, N'
      Stop ' ERROR IN NETMAT, N'
    End If

!   ----------output (last iteration)

    If (unasol .And. megas) Then
      casta1 = xni*x1
      casta2 = xnk*x2*xxlev
!     WRITE (100,FMT=9000) LL,'RBF',ML,MU,x1,x2*xxlev, &
      Write (100, Fmt=140) ll, 'RBF', ml, mu, casta1, casta2, casta2 - casta1
!     NIST*RIK, NIST*RKI, NIST*RKI-DEP*NIST*RIK
!     WRITE (100,FMT=9000) LL,'RBF',ML,MU,CASTA1*NIST(ML)/XNI,CASTA2/XNK, &
!     (CASTA2 - CASTA1)/XNK
    End If

!   ----------set up of matrix

    If (x1==0.D0 .Or. x2*xxlev==0.D0) Then
      Write (999, *) ' STOP: INTERMEPHO'
      Stop ' INTERMEPHO'
    End If

!   if(nato.eq.3)print*,ll,n1,n2,x1,x2,depart,xxlev
    ratmat(n1, n1) = ratmat(n1, n1) - x1
    ratmat(n1, n2) = ratmat(n1, n2) + x2*xxlev ! ok, since with respect to n_k
!   (and not n_u)

!   ---makes only sense for intermediate ions, couples lower to upper ones

    ratmat(n2, n1) = ratmat(n2, n1) + x1
    ratmat(n2, n2) = ratmat(n2, n2) - x2*xxlev

    mu1 = kl(ml)
    If (.Not. op_flag(i) .Or. mu==mu1) Cycle rbfloop

!   additional treatment for OP-data with ionization to excited states
!   (resonances between FRECFN and FRECIN)

!   consistency tests
    If (frecfn(i)==frecin(i)) Then
      Print *, labl(labl4(i)), frecin(i), frecfn(i), frecfn1(i)
      Write (999, *) &
        ' STOP: PROBLEM IN ADDITIONAL PATH OF RBF-TREATMENT (OP-DATA)'
      Stop ' PROBLEM IN ADDITIONAL PATH OF RBF-TREATMENT (OP-DATA)'
    End If

    If (indfrq1(i)==0) Then
      Write (999, *) &
        ' STOP: INDFRQ1 = 0 IN ADDITIONAL PATH OF RBF-TREATMENT (OP-DATA)'
      Stop ' INDFRQ1 = 0 IN ADDITIONAL PATH OF RBF-TREATMENT (OP-DATA)'
    End If

!   upper state now ground-state of next ion
!   (this is an approximation, particularly  when excited state = 3 or higher
!   -- in the latter case, upper level could be the 2nd state etc.)
    mu = mu1
    If (mu>nlas) Then
      Write (999, *) &
        ' STOP: MU > NLAS IN ADDITIONAL PATH OF RBF-TREATMENT (OP-DATA)'
      Stop ' MU > NLAS IN ADDITIONAL PATH OF RBF-TREATMENT (OP-DATA)'
    End If

    xnk = blev0(mu)
    xxlev = alevel(ml)/alevel(mu)*xnerat

    xnistar = xnk*xxlev !           n_u *(n_i/n_u)*
    nist(ml) = xxlev !              witout nk
    xnkk(ml) = xnk
    depart = xni/xnistar

!   -------   additional photo-rates (integration starting at INDFRQ1(I),
!   corresponding to FRECFN1, until FRECIN; analogous to INTERMEDR)

!   no rates for heating and cooling (already treated above)
    Call intermepho1(i, depart, xnistar, ll, x1, x2, xne, temp, clf, optmet)

    n2 = mu - nfir + 1

    If (n2>mn) Then
      Write (999, *) ' STOP: ERROR IN NETMAT, N2 (additional rbf-path)'
      Stop ' ERROR IN NETMAT, N2 (additional rbf-path)'
    End If


!   ----------output (last iteration)

    If (unasol .And. megas) Then
      casta1 = xni*x1
      casta2 = xnk*x2*xxlev
!     WRITE (100,FMT=9001) LL,'RBF1',ML,MU,x1,x2*xxlev, &
      Write (100, Fmt=150) ll, 'RBF1', ml, mu, casta1, casta2, casta2 - casta1
!     NIST*RIK, NIST*RKI, NIST*RKI-DEP*NIST*RIK
!     WRITE (100,FMT=9000) LL,'RBF',ML,MU,CASTA1*NIST(ML)/XNI,CASTA2/XNK, &
!     (CASTA2 - CASTA1)/XNK
    End If

!   ----------set up of matrix

    If (x1==0.D0 .Or. x2*xxlev==0.D0) Then
      Write (999, *) ' STOP: INTERMEPHO (additional path)'
      Stop ' INTERMEPHO (additional path)'
    End If

!   now additional rates (de-) populating gs of next ion (analogous to
!   INTERMEDR)
!   if(nato.eq.3)print*,ll,n1,n2,x1,x2,depart,xxlev
    ratmat(n1, n1) = ratmat(n1, n1) - x1
    ratmat(n1, n2) = ratmat(n1, n2) + x2*xxlev ! ok, since with respect to n_k
!   (and not n_u)

!   ---makes only sense for intermediate ions, couples lower to upper ones

    ratmat(n2, n1) = ratmat(n2, n1) + x1
    ratmat(n2, n2) = ratmat(n2, n2) - x2*xxlev

  End Do rbfloop


! ------ explicit DR-transitions

  If (.Not. optdr) Go To 100
drloop: Do i = 1, index4
    If (drexplicit(i)==0) Cycle
    ml = mlow(i)
    mu = kl(ml) !                   w.r.t. IONIC GROUND-STATE!

    If (ml<nfir .Or. ml>=nlas) Cycle

    If (mu>nlas) Then
!     SHOULD NEVER HAPPEN. THUS
      Write (999, *) ' STOP: UPPER STATE OF DR TRANSITION NOT GROUND-STATE'
      Stop ' UPPER STATE OF DR TRANSITION NOT GROUND-STATE'
    End If

    xxlev = alevel(ml)/alevel(mu)*xnerat

!   -------   DR-rates

    Call intermedr(i, ll, x1, x2, temp)

!   HEATING COOLING RATES NOT INCLUDED SO FAR

    n1 = ml - nfir + 1
    n2 = mu - nfir + 1

    If (n2>mn .Or. n1<0) Then
      Write (999, *) ' STOP: ERROR IN NETMAT, N'
      Stop ' ERROR IN NETMAT, N'
    End If

!   ----------output (last iteration)

    If (unasol .And. megas) Then
      xni = blev0(ml)
      xnk = blev0(mu)
      casta1 = xni*x1
      casta2 = xnk*x2*xxlev
!     WRITE (100,FMT=9000) LL,'RDR',ML,MU,x1,x2*xxlev, &
      Write (100, Fmt=140) ll, 'RDR', ml, mu, casta1, casta2, casta2 - casta1
!     NIST*RIK, NIST*RKI, NIST*RKI-DEP*NIST*RIK
!     WRITE (100,FMT=9000) LL,'RDR',ML,MU,CASTA1*XXLEV/XNI,CASTA2/XNK, &
!     (CASTA2 - CASTA1)/XNK
    End If

!   ----------set up of matrix

    If (x1==0.D0 .Or. x2*xxlev==0.D0) Then
      Write (999, *) ' STOP: INTERMEDR'
      Stop ' INTERMEDR'
    End If

    ratmat(n1, n1) = ratmat(n1, n1) - x1
    ratmat(n1, n2) = ratmat(n1, n2) + x2*xxlev ! ok, since with respect to n_k
!   (and not n_u)

!   ---makes only sense for intermediate ions, couples lower to upper ones

    ratmat(n2, n1) = ratmat(n2, n1) + x1
    ratmat(n2, n2) = ratmat(n2, n2) - x2*xxlev
  End Do drloop

100 Continue
! ---
! --- Auger-ionization
! --- in its present form (see table), only for K-shell ionization
! --- NOT only for two ejected electrons, but generalized

  If (optxray .And. optauger) Then
    Do k = 1, 30
      If (name(k)==labat(nato)) Go To 110
    End Do
    Write (999, *) ' STOP: LABAT(NATO) NOT FOUND IN NAME'
    Stop ' LABAT(NATO) NOT FOUND IN NAME'

110 Continue

augerloop: Do j = ilo, ima
      jeff = j + zeff(nato)

      If (jeff<k_nmin .Or. jeff>k_nmax) Cycle augerloop

      Do l = 1, n_kedges
        If (z_k(l)==k .And. n_k(l)==jeff) Go To 120
      End Do
!     ELEMENT (OR ION) NOT PRESENT IN K_SHELL_AUGER_DATA
      Cycle augerloop

120   itrans = l
      If (eth(itrans)>=fre(ifre)) Cycle augerloop

      Call intermepho_auger(x1, itrans, ll)

      numion = igenio(nato, j)
      nfm = ifirsl(numion)
      nlm = ifirsl(numion+1) - 1


      Do ml = nfm, nlm

!       LOOP OVER ALL PROBABILITIES (1,2,3 ETC. EJECTED ELECTRONS)
proba:  Do jj = j + 1, ima + 1
          mu = ifirsl(igenio(nato,jj))
          n1 = ml - nfir + 1
          n2 = mu - nfir + 1
!         E.G., IF JJ=J+1, WE HAVE JJ-J = 1 EJECTED ELECTRONS, FOR
!         JJ=J+2  WE HAVE JJ-J = 2 EJECTED ELECTRONS (TYPICAL CASE), AND SO ON
          prob = aug_1_6(itrans, jj-j)
          If (prob==0.D0) Cycle proba

!         FOR TESTS
!         IF(LL.EQ.1) PRINT*,K,J,ITRANS,JJ,ML,MU,N1,N2,PROB,X1

!         ----------output (last iteration)

          If (unasol .And. megas) Then
            xni = blev0(ml)
            casta1 = xni*x1*prob
            casta2 = 0.
            Write (100, Fmt=140) ll, 'AUG', ml, mu, casta1, casta2, &
              casta2 - casta1
          End If
!         only upwards rates
          ratmat(n1, n1) = ratmat(n1, n1) - x1*prob
          ratmat(n2, n1) = ratmat(n2, n1) + x1*prob
        End Do proba
      End Do
    End Do augerloop
  End If !                          AUGER


! ------- collisional b-f transitions

cbfloop: Do i = 1, index3
    ml = mclow(i)
    mu = mcup(i)

    If (ml<nfir .Or. ml>=nlas) Cycle

    xni = blev0(ml)

    If (mu>nlas) Then
!     this path refers to ionization to excited levels of imax+1
!     upper level assumed to be in lte with respect to nlas (groundstate of
!     imax+1)
!     see above
      mu = nlas
    End If
    xnk = blev0(mu)
    xxlev = alevel(ml)/alevel(mu)*xnerat

!   for tests
    ij = igenio(le(mu), li(mu))
    lo = ifirsl(ij)
    If (lo==mu) Then
!     frequencies can be different if actual upper level > NLAS
!     in this case, use actual frequency (FLCBF(I)),
!     but add rate ground-state of imax+1
      If (flcbf(i)/=fl(ml) .And. mu==mcup(i)) Then
        Print *, flcbf(i), fl(ml), labl3(i), ml, mu, ilo, ima, ll
        Write (999, *) &
          ' STOP: PROBLEMS WITH GROUND-STATE IONIZATION FREQUENCY (CBF)'
        Stop ' PROBLEMS WITH GROUND-STATE IONIZATION FREQUENCY (CBF)'
      End If
    End If

    Call intermecol(i, temp, cik, fl, zl)
    cik = cik*xne !                 includes clumping

    If (opttcor) Then
      Call intermecol(i, temp*1.01, cik1, fl, zl)
!     partial derivative
      dtcik = (cik1*xne-cik)/(0.01*temp)
!     JO CHANGED DEC. 2016
!     HHNU = HH*FL(ML) ! to be consistent with LTE
      hhnu = hh*flcbf(i)
      If (hhnu<=0.D0) Then
        Write (999, *) ' STOP: ERROR IN HNU - CBF'
        Stop ' ERROR IN HNU - CBF'
      End If

      qcrec(mu, 1) = qcrec(mu, 1) + hhnu*xxlev*cik ! (*n_k)
      qcion(ml, 1) = qcion(ml, 1) + hhnu*cik ! (*n_i)
!     DTQCBFH(LL,NATO) = DTQCBFH(LL,NATO) + HHNU*XNK*XLEVEL*DTCIK !does not
!     work
      qcrec(mu, 2) = qcrec(mu, 2) + hhnu*xxlev*cik*(dtcik/cik-(3.D0/2.D0+hhnu/ &
        akb/temp)/temp) !           (*n_k)
!     Here comes the question! Why does the first formulation give a better
!     flux, however leads to a destabilization, whereas the 2nd one give a
!     stable
!     convergence?
!     DTQCBFC(LL,NATO) = DTQCBFC(LL,NATO) + HHNU*XNI*DTCIK
      qcion(ml, 2) = qcion(ml, 2) + hhnu*cik*(dtcik/cik-(3.D0/2.D0+hhnu/akb/ &
        temp)/temp) !               (*n_i)!

    End If

    n1 = ml - nfir + 1
    n2 = mu - nfir + 1
    If (n2>mn .Or. n1<0) Then
      Write (999, *) ' STOP: ERROR IN NETMAT, N'
      Stop ' ERROR IN NETMAT, N'
    End If

!   ----------output (last iteration)

    If (unasol .And. megas) Then
      casta1 = xni*cik
      casta2 = xnk*cik*xxlev
!     WRITE (100,FMT=9000) LL,'CBF',ML,MU,Cik,Cik*xxlev, &
      Write (100, Fmt=140) ll, 'CBF', ml, mu, casta1, casta2, casta2 - casta1
!     NIST*CIK, NIST*CIK, NIST*CIK-DEP*NIST*CIK, with CKI=(NI/NK)_star*CIK
!     WRITE (100,FMT=9000) LL,'CBF',ML,MU,CASTA1*NIST(ML)/XNI,CASTA2/XNK, &
!     (CASTA2 - CASTA1)/XNK
    End If

!   ----------set up of matrix

    If (cik==0.D0 .Or. xxlev==0.D0) Stop ' intermecoll'
!   if(nato.eq.3)print*,ll,n1,n2,cik,depart,xxlev
    ratmat(n1, n1) = ratmat(n1, n1) - cik
    ratmat(n1, n2) = ratmat(n1, n2) + cik*xxlev

!   ---makes only sense for intermediate ions, couples lower to upper ones

    ratmat(n2, n1) = ratmat(n2, n1) + cik
    ratmat(n2, n2) = ratmat(n2, n2) - cik*xxlev

  End Do cbfloop

! -----  line treatment

  If (concon) Go To 130

! changed, improved approach (see linesob)
! ------ last row of ratmat is explicitly overwritten (pure cont.)

! DO I = 1,MN
! RATMAT(MN,I) = 1.D0
! END DO

! F(MN) = ABUNDND(NATO,LL)

  jmax = mn
  bmax = 0.

  Do i = 1, mn
    j = i + nfir - 1
    bold(i) = blev0(j)
    If (bold(i)>bmax) Then
      bmax = bold(i)
      jmax = i
    End If
  End Do

  Do i = 1, mn
    ratmat(jmax, i) = 1.D0
  End Do

  f(jmax) = abundnd(nato, ll)


  flud(1:mn) = f(1:mn)
  rat2(1:mn, 1:mn) = ratmat(1:mn, 1:mn)


! -----solution of rate equations

  Call ludcmp(ratmat, mn, ndim, indlud, det)

  Call lubksb(ratmat, mn, ndim, indlud, f)

! -----------------------------------------------------------------
! --------newton raphson improvement-cycle-------------------------
! -----------------------------------------------------------------

errit: Do iii = 1, itm

!   old version (if no qp precision available)
!   DO I = 1,MN
!   SDP = -FLUD(I)

!   high precision summation

!   DO J=1,MN
!   AUX(J)=RAT2(I,J)*F(J)
!   END DO

!   CALL SORT1(MN,AUX)

!   DO J=1,MN
!   SDP=SDP+AUX(J)
!   END DO

!   DELTAN(I) = SDP
!   END DO

!   new version (residuum calculated in qp precision
    aux(1:mn) = matmul(real(rat2(1:mn,1:mn),qp), real(f(1:mn),qp)) - real(flud &
      (1:mn), qp)
    deltan(1:mn) = aux(1:mn)
    Call lubksb(ratmat, mn, ndim, indlud, deltan)

    errmax = 0.D0

    Do i = 1, mn
      f(i) = f(i) - deltan(i)
      errmax = max(errmax, abs(deltan(i)/f(i)))
    End Do

!   print*,'errmax: ',ll,errmax

    If (errmax<precis) Then
      Do i = 1, mn
        j = i + nfir - 1
        blevel(j) = f(i)
      End Do
      Go To 130
    End If

  End Do errit

  Do i = 1, mn
    j = i + nfir - 1
    blevel(j) = f(i)
  End Do

  If (errmax>1.D-3) Then
    Print *, ' WARNING!! N-R CONVERGENCE NOT ACHIEVED IN ', iii, &
      ' ITERATIONS BY', errmax, ' (NETMAT)'
    Print *, 'ATOM: ', nato, '   DEPTH: ', ll
    Write (999, *) ' WARNING!! N-R CONVERGENCE NOT ACHIEVED IN ', iii, &
      ' ITERATIONS BY', errmax, ' (NETMAT)'
    Write (999, *) 'ATOM: ', nato, '   DEPTH: ', ll
!   EMERGENCY
    If (errmax>0.5) Then
      Do i = 1, mn
        j = i + nfir - 1
        blevel(j) = blev0(i)
      End Do
    End If
  End If

130 Continue

  If (concon) Call linesob(ll, nfir, nlas, blev0, ratmat, ndim, &
    xne, nato, temp, errold, qlu, qld)

  Return

140 Format (1X, I2, 3X, A3, 2(3X,I3), 3(2X,G16.6))
150 Format (1X, I2, 3X, A4, 2X, I3, 3X, I3, 3(2X,G16.6))

End Subroutine

!-----------------------------------------------------------------------

Subroutine intermecol(ntr, temp, cik, fl, zl)

  Use :: nlte_type
  Use :: nlte_dim
  Use :: fund_const
  Use :: princesa_var, Only: data, iform3, indat3, labl3, flcbf

  Implicit None

! ---- it calculates collisional bound-free coeffs.
! ---- used (and checked so far)

! ------ formula 1 : simple cbf coefficient, a a la Mihalas (for Adi's input)
! as formula 15, but with GQ as input
! ------ formula 10: h i (all levels without level 2)
! ------ formula 11: h i (level 2), he ii
! ------ formula 12: he ii
! ------ formula 14: he i
! ------ formula 15: he ii, si ii,iii,iv
! ------ formula 17: n iii  (Feb. 2006: bug found by JO)
! ------ formula 20: mg ii
! ------ formula 21: he i

! JO: changed Dec. 2016, to account for ionization to excited levels

! .. scalar arguments ..
  Real (dp) :: cik, temp
  Integer (i4b) :: ntr
! ..
! .. array arguments ..
  Real (dp) :: fl(id_llevs), zl(id_llevs)

! ..
! .. local scalars ..
  Real (dp) :: aa1, anon, apar, gamma, gq, tlog, usu1, usu2, usub0, xnu0, z, &
    xxx, e(10)
  Integer (i4b) :: ida, iz, nnf, nnl, nn, n, i
! ..
! .. external functions ..
  Real (dp) :: expin1
  External :: expin1
! ..
! .. intrinsic functions ..
! INTRINSIC EXP,INT,LOG10,SQRT
! ..

  aa1 = 5.465727564D-11
  nnl = labl3(ntr)
! XNU0 = FL(NNL)
  xnu0 = flcbf(ntr)
  nnf = iform3(ntr)


! thus far, we only allow for ionization to excited states for formula 1 and
! 15
  If (nnf/=1 .And. nnf/=15) Then
    If (xnu0/=fl(nnl)) Then
      Write (999, *) ' STOP: CBF FORMULA 17 AND IONIZATION TO EXCITED LEVEL'
      Stop ' CBF FORMULA 17 AND IONIZATION TO EXCITED LEVEL'
    End If
  End If

  usub0 = hh*xnu0/akb/temp
  ida = indat3(ntr)


  z = zl(nnl) + 1.D0

  iz = int(z+0.5D0)

  If (nnf==1) Then
    cik = 1.55D13*data(ida)*data(ida+1)*exp(-usub0)/(sqrt(temp)*usub0)

  Else If (nnf==11) Then
    gamma = data(ida)/(temp**2) + data(ida+1)/temp + data(ida+2) + &
      data(ida+3)*temp + data(ida+4)*temp**2

    cik = aa1*gamma*exp(-usub0)*sqrt(temp)

  Else If (nnf==12) Then
    tlog = log10(temp)
    gamma = data(ida)/(tlog**2) + data(ida+1)/tlog + data(ida+2) + &
      data(ida+3)*tlog + data(ida+4)*tlog**2

    cik = aa1*gamma*exp(-usub0)*sqrt(temp)

  Else If (nnf==13) Then
    gamma = (data(ida)/z)**4

    cik = aa1*gamma*exp(-usub0)*sqrt(temp)

  Else If (nnf==14) Then
    usu1 = usub0 + 0.27D0
    usu2 = usub0 + 1.43D0
    apar = usub0*expin1(usub0)*exp(-usub0) - 0.728D0*usub0**2/usu1*expin1(usu1 &
      )*exp(-usu1)
    anon = 0.189D0*usub0**2*exp(-usub0)*(2.D0+usu2)/usu2**3

    cik = (apar-anon)*data(ida)*aa1/sqrt(temp)

  Else If (nnf==15) Then

    If (iz==1) Then
      gq = .1D0
    Else If (iz==2) Then
      gq = .2D0
    Else
      gq = .3D0
    End If

    cik = 1.55D13*exp(-usub0)*data(ida)*gq/usub0/sqrt(temp)

  Else If (nnf==17) Then !          checked by JO, bug identified, changed
!   there is still an inconsistency between
!   the Detail program and the manual.
!   original source by Moore 72 need to be checked
    n = int(data(ida)+0.5D0)
    gq = 3.039652336D-16*xnu0 !     nu/nu_hyd, d.h. Ei in Ryd
    usu1 = hh*clight/akb/temp*1.097373177D5 ! = h/kt * nu_hyd = u_ion for hyd.
    e(1) = expin1(usub0)*exp(-usub0)
    Do i = 2, n
      e(i) = (exp(-usub0)-usub0*e(i-1))/(float(i-1))
    End Do
    gamma = 0.D0
    Do i = 1, n
      gamma = gamma + e(i)*data(ida+2+i)
    End Do
    gamma = gamma*usub0
    gamma = gamma + data(ida+1)*e(1) + data(ida+2)*exp(-usub0)
    gamma = gamma*usu1/gq !         usu1/gq = u_hyd/(nu/nu_hyd) = uh/ei

    cik = aa1*sqrt(temp)*gamma

  Else If (nnf==20) Then !          taken from DETAIL, not checked by JO

    gamma = hh*clight/akb/temp*1.097373177D5/usub0
    usu1 = usub0 + data(ida+2)
    gamma = gamma**2*(data(ida)*usub0*expin1(usub0)*exp(-1.D0*usub0)+data(ida+ &
      1)*(usub0/usu1)**2*(expin1(usu1)*exp(-1.D0*usu1)+exp(-1.D0*usu1)))
    cik = aa1*sqrt(temp)*gamma

  Else If (nnf==21) Then
    nn = int(data(ida)+5.D-1)
!   ------- gamma is calculated with electron density equals to unity
!   ------- because is scaled by it in netmat
    xxx = 1.D0
    Call hecion(temp, xxx, nn, gamma)
    cik = gamma

  Else
    Print *, 'TRANSITION CBF NOT IMPLEMENTED, NFOR=', nnf
    Write (999, *) 'STOP!'
    Stop

  End If

  Return

End Subroutine

!-----------------------------------------------------------------------

Subroutine intermepho(ntr, depart, xnistar, ll, x1, x2, xne, temp, clf, &
  optmet, x3, x4, x5, x6)

! ----  several quantities 'back-corrected' here, since
! ----  OPAC variable EFFECTIVE opacity. See comments.

! ----  comment by JO (Jan. 2016):
! the scattering component of ALO (pure Thomson so far)
! might need to be updated regarding the line pseudo continuum
! (note that this refers to the 'old' approach, i.e.,
! it is an old problem)

  Use :: nlte_type
  Use :: nlte_dim
  Use :: fund_const

  Use :: nlte_var, Only: alphafs, fre1 => fref, wfre1 => wfref, indfin, &
    indfrq, ifre1 => ifref, fre, ifre, opac, xj, alo, opat_m_old, strue_m_old, &
    strue_m_new, metals_converged, optcmf_all

  Use :: nlte_porvor, Only: opa_eff_rat, opa_eff_rat_old

  Implicit None

! ---  new algorithm
! ---  subroutine intemepho (see netmat). it calculates photo-integrals
! ---  for some the formulae implemented in detail.

! .. parameters ..
! Integer (i4b), Parameter :: nd1 = id_ndept
! Integer (i4b), Parameter :: ifretot = id_frec1
! Integer (i4b), Parameter :: ifrefin = id_frec2
! ..
! .. scalar arguments ..
  Real (dp) :: depart, temp, x1, x2, xne, x3, x4, x5, x6, edgenu, xnistar, clf
  Integer (i4b) :: ll, ntr
  Logical :: optmet, optmet1
! ..
! .. local scalars ..
  Real (dp) :: alovth, alpha, aux1, aux2, dex, dfreq, dfx, dr, exk, f, f1, f2, &
    facf, facg, fgconst, freq, freqii, freqii1, g, g1, g2, hcfre3, hckt, &
    hnue1, opa, pi4, sig, sum1, sum2, thom, thomalo, wf, xalo, xf, xf1, xf2, &
    xfreq, xfreq1, xg, xg1, xg2, xxj, betar, sum3, sum4, sum5, sum6, xxjj, &
    aux11, ff, gg, ff1, ff2, gg1, gg2, xff1, xgg1, xff2, xgg2, facff, facgg, &
    xff, xgg
  Integer (i4b) :: ii, ii1, iiold, istart, kk
! ..
! .. external subroutines ..
  External :: crossbf
! ..
! .. intrinsic functions ..
! INTRINSIC ABS,EXP,LOG10,MAX
! ..

  If (optcmf_all) Then
    If (xj(1,ifre+1)/=2) Then
      Write (999, *) ' STOP: OPTCMF_ALL = TRUE, BUT APPROX. XJ &
        &(OBSFRAM) USED IN INTERMEPHO'
      Stop ' OPTCMF_ALL = TRUE, BUT APPROX. XJ (OBSFRAM) USED IN INTERMEPHO'
    End If
!   remember: ALO=0 for cmf_all frequencies
  End If


  optmet1 = optmet
  If (optmet1 .And. metals_converged) optmet1 = .False.

  hckt = hkl/temp

! begin of integration of the four photo-integrals
! (for definition, see the santolaya et al. paper)
! note: integration over frequency, fre in kayser

  iiold = 0
! not used:  rerr = 1.D-5
  sig = amh*sigmae
  istart = indfrq(ntr)
  If (istart<=0) Then
    Print *, 'TRANSITION', ntr, ' AT', ll
    Print *, 'INDFRQ NOT DEFINED: CHECK ILOW OR IMAX (LTE AND NLTE!!!)'
    Write (999, *) ' STOP: INDFRQ NOT DEFINED: CHECK ILOW OR IMAX &
      &(LTE AND NLTE!!!) - FOR TESTS'
    Stop &
      ' INDFRQ NOT DEFINED: CHECK ILOW OR IMAX (LTE AND NLTE!!!) - FOR TESTS'
  End If

  sum1 = 0.D0
  sum2 = 0.D0
  ff = 0.D0 !                       maup

! ---    note: first integration weight has to be modified

  kk = istart
  wf = wfre1(kk) - .5D0*clight*(fre1(kk)-fre1(kk-1))
  freq = fre1(kk)

! ---    ii index corresponding to coarse freq. grid

  edgenu = freq !                   maup
! IF(LL.EQ.1. .AND. OP_FLAG(NTR)) PRINT*,LABL(LABL4(NTR)), EDGENU

  ii = indfin(kk)

  If (ii>0.) Then
    alpha = alphafs(ntr, ii)
  Else
    Write (999, *) &
      ' STOP: II MUST NOT BE LT ZERO IN INDFIN AT FIRST FREQ. POINT'
    Stop ' II MUST NOT BE LT ZERO IN INDFIN AT FIRST FREQ. POINT'
  End If


! OPA = OPAC(LL,II)*CLF !corrected
  opa = opac(ll, ii)*clf/opa_eff_rat_old(ll, ii) ! back-corrected
! ---- since OPAC is already an effective quantity, back-corrected here
! (and thus also below).
! NOTE: Use OLD values, since opat_m_old is used below
  thom = xne*sig/opa
  thomalo = 1.D0 - alo(ll, ii)*thom
  alovth = alo(ll, ii)/thomalo
  xxj = xj(ll, ii)
  If (optmet1) Then
    betar = opat_m_old(ll, ii)*clf/opa ! corrected (OPAT_M_OLD is mean
!   opacity)
    xxj = xxj + alovth*betar*(strue_m_new(ll,ii)-strue_m_old(ll,ii))
  End If

  exk = exp(-hckt*freq)
  hcfre3 = hc2*freq**3
  hnue1 = 2.D0/(hc2*freq)

! ---    f corresponds to i3-i4, g to i1-i2

  f = alpha*hnue1*(xxj-hcfre3*exk*alovth*alpha*xnistar/opa)
  g = alpha*hnue1*exk*(xxj+hcfre3*(1.D0-alovth*alpha*depart*xnistar/opa))

  sum1 = sum1 + f*wf
  sum2 = sum2 + g*wf
! SUM3,SUM4,SUM5 and SUM6 = 0. due to (1-EDGENU/nu) factor
  sum3 = 0.D0
  sum4 = 0.D0
  sum5 = 0.D0
  sum6 = 0.D0

! ----   rest of integration, with all necessary interpolations

kkloop: Do kk = istart + 1, ifre1

    wf = wfre1(kk)
    freq = fre1(kk)
    ii = indfin(kk)

    If (ii>0.) Then
      alpha = alphafs(ntr, ii)
!     OPA = OPAC(LL,II)*CLF !corrected
      opa = opac(ll, ii)*clf/opa_eff_rat_old(ll, ii) ! back-corrected
      thom = xne*sig/opa
      xalo = alo(ll, ii)
      thomalo = 1.D0 - xalo*thom
      alovth = xalo/thomalo
      xxj = xj(ll, ii)
      xxjj = xxj
      If (optmet1) Then
        betar = opat_m_old(ll, ii)*clf/opa ! corrected
        xxj = xxj + alovth*betar*(strue_m_new(ll,ii)-strue_m_old(ll,ii))
      End If
      exk = exp(-hckt*freq)
      hcfre3 = hc2*freq**3
      hnue1 = 2.D0/(hc2*freq)
      f = alpha*hnue1*(xxj-hcfre3*exk*alovth*alpha*xnistar/opa)
      g = alpha*hnue1*exk*(xxj+hcfre3*(1.D0-alovth*alpha*depart*xnistar/opa))
      ff = alpha*hnue1*xxjj
      gg = alpha*hnue1*exk*(xxjj+hcfre3)
!     FF=ALPHA*HNUE1*HCFRE3*EXK*ALOVTH*ALPHA*XNISTAR/OPA


      If (f<0. .Or. g<0.) Then
        Print *, ntr, fre1(istart), kk, f, g, xxj, alovth, depart
        Print *, ii, freq, opa, xnistar
        Print *, alphafs(ntr, ii)
        Print *, 'EFFECTIVE OPACITY ISSUE?:', opa_eff_rat(ll, ii), &
          opa_eff_rat_old(ll, ii)
        Print *, opac(ll, ii), clf
        Print *, 'NEGATIVE CONTRIBUTION TO PHOTO-RATES1'
        Write (999, *) ' STOP!'
        Stop
      End If

    Else If (ii<0) Then

!     logarithmic interpolation for desired quantities

      xf = log10(freq)
      ii = -ii
      If (ii==iiold) Go To 100
      iiold = ii
      ii1 = ii + 1
!     OPA = OPAC(LL,II)*CLF !corrected
      opa = opac(ll, ii)*clf/opa_eff_rat_old(ll, ii) ! back-corrected
      thom = xne*sig/opa
      xalo = alo(ll, ii)
      thomalo = 1.D0 - xalo*thom
      alovth = xalo/thomalo
      xxj = xj(ll, ii)
      xxjj = xxj
      If (optmet1) Then
        betar = opat_m_old(ll, ii)*clf/opa ! corrected
        xxj = xxj + alovth*betar*(strue_m_new(ll,ii)-strue_m_old(ll,ii))
      End If
      freqii = fre(ii)
      exk = exp(-hckt*freqii)
      hcfre3 = hc2*freqii**3

      alpha = alphafs(ntr, ii)
      aux1 = xxj/(hcfre3*alpha)
      aux11 = xxjj/(hcfre3*alpha)
      aux2 = alovth*xnistar/opa
!     new evalution scheme, to allow for very low EXK (x-ray regime)
!     F1 = AUX1/EXK - AUX2
      f1 = aux1 - aux2*exk
      g1 = 1.D0/alpha + aux1 - depart*aux2
!     FF1 = AUX11/EXK
      ff1 = aux11
      gg1 = 1.D0/alpha + aux11


!     OPA = OPAC(LL,II1)*CLF  !corrected
      opa = opac(ll, ii1)*clf/opa_eff_rat_old(ll, ii1) ! back-corrected
      thom = xne*sig/opa
      xalo = alo(ll, ii1)
      thomalo = 1.D0 - xalo*thom
      alovth = xalo/thomalo
      xxj = xj(ll, ii1)
      xxjj = xxj
      If (optmet1) Then
        betar = opat_m_old(ll, ii1)*clf/opa ! corrected
        xxj = xxj + alovth*betar*(strue_m_new(ll,ii1)-strue_m_old(ll,ii1))
      End If
      freqii1 = fre(ii1)
      exk = exp(-hckt*freqii1)
      hcfre3 = hc2*freqii1**3

      alpha = alphafs(ntr, ii1)
      aux1 = xxj/(hcfre3*alphafs(ntr,ii1))
      aux11 = xxjj/(hcfre3*alpha)
      aux2 = alovth*xnistar/opa
!     F2 = AUX1/EXK - AUX2
      f2 = aux1 - aux2*exk
      g2 = 1.D0/alpha + aux1 - depart*aux2
!     FF2 = AUX11/EXK
      ff2 = aux11
      gg2 = 1.D0/alpha + aux11

      If (f1<0. .Or. g1<0. .Or. f2<0. .Or. g2<0.) Then
        Print *, ntr, fre1(istart), kk, f1, g1, f2, g2, xxj, alovth, depart
        Print *, ii, ii1, freqii, freqii1
        Print *, alphafs(ntr, ii), alphafs(ntr, ii1)
        Print *, 'effective issue?:', opa_eff_rat(ll, ii), &
          opa_eff_rat_old(ll, ii)
        Print *, opac(ll, ii), clf
        Print *, 'NEGATIVE CONTRIBUTION TO PHOTO-RATES1'
        Print *, 'NEGATIVE CONTRIBUTION TO PHOTO-RATES2'
        Write (999, *) ' STOP!'
        Stop
      End If

      xfreq = log10(freqii)
      xfreq1 = log10(freqii1)
!     XF1 = LOG10(F1)
!     XF2 = LOG10(F2)
!     ADD log10(exp(hnu/kt))
      xf1 = log10(f1) + hckt*freqii*log10e
      xf2 = log10(f2) + hckt*freqii1*log10e
      xg1 = log10(g1)
      xg2 = log10(g2)

!     XFF1 = LOG10(FF1)
!     XFF2 = LOG10(FF2)
!     ADD log10(exp(hnu/kt))
      xff1 = log10(ff1) + hckt*freqii*log10e
      xff2 = log10(ff2) + hckt*freqii1*log10e
      xgg1 = log10(gg1)
      xgg2 = log10(gg2)


      dfreq = 1.D0/(xfreq1-xfreq)
      facf = (xf2-xf1)*dfreq
      facg = (xg2-xg1)*dfreq
      facff = (xff2-xff1)*dfreq
      facgg = (xgg2-xgg1)*dfreq

!     xff1 = LOG10(FF1)
!     xff2 = LOG10(FF2)
!     facff = (xff2-xff1)*DFREQ

100   Continue

      dfx = xf - xfreq
      xf = xf1 + facf*dfx
      xg = xg1 + facg*dfx
!     XF = 10.D0**XF
!     XG = 10.D0**XG

      xff = xff1 + facff*dfx
      xgg = xgg1 + facgg*dfx
!     XFF = 10.D0**XFF
!     XGG = 10.D0**XGG
      Call crossbf(ntr, alpha, freq)
      hnue1 = 2.D0/(hc2*freq)
      exk = exp(-hckt*freq)
      hcfre3 = hc2*freq**3
      fgconst = log10(alpha*hnue1*hcfre3*alpha) - hckt*freq*log10e
      f = 10.D0**(fgconst+xf)
      g = 10.D0**(fgconst+xg)
      ff = 10.D0**(fgconst+xff)
      gg = 10.D0**(fgconst+xgg)

!     XFF = XFF1 + FACF*DFX
!     XFF = 10.D0**XFF
!     FF = XFF*FGCONST

    Else
      Write (999, *) ' STOP: II MUST NOT BE ZERO IN INDFIN'
      Stop ' II MUST NOT BE ZERO IN INDFIN'
    End If

    x1 = f*wf
    x2 = g*wf

    x3 = ff*wf*(1.D0-edgenu/freq)*hh*clight*freq
    x4 = gg*wf*(1.D0-edgenu/freq)*hh*clight*freq
    x5 = x3*(hckt*freq)/temp
    x6 = x4*(hckt*freq)/temp

    sum1 = sum1 + x1
    sum2 = sum2 + x2
    sum3 = sum3 + x3
    sum4 = sum4 + x4
    sum5 = sum5 + x5
    sum6 = sum6 + x6

    dr = abs(x1/sum1)
    dr = max(dr, abs(x2/sum2))
    dr = max(dr, abs(x3/sum3))
    dr = max(dr, abs(x4/sum4))
    dr = max(dr, abs(x6/sum6))

    dex = 1.D0 - fre1(kk-1)/freq

!   IF (DEX.LT.1.D-7) STOP ' ERROR IN DEX'

!   The limit has been changed to allow for close energy terms
!   in detailed models for neutral ions (energy terms with almost
!   the same ionization energy). Note that this will increase the
!   computation time, but I would not expect a huge effect
    If (dex<1.D-9) Then
      Write (999, *) ' STOP: ERROR IN DEX'
      Stop ' ERROR IN DEX'
    End If
!   IF (DR.LT.RERR .AND. DEX.GT..01D0) EXIT
  End Do kkloop

  pi4 = 4.D0*pi
  x1 = pi4*sum1
  x2 = pi4*sum2

  x3 = pi4*sum3
  x4 = pi4*sum4
  x5 = pi4*sum5
  x6 = pi4*sum6

! print*,ntr,x1,' ',x2,' ',dr,' ',istart,' ',kk

  Return

End Subroutine

!-----------------------------------------------------------------------

Subroutine intermepho1(ntr, depart, xnistar, ll, x1, x2, xne, temp, clf, &
  optmet)

! ----  as subroutine INTERMEPHO, but only calculating X1 and X2, in
! the range INDFRQ1(NTR) ... FRECIN(NTR), assuming ionization to ground-state


! X1 and X2 are the additional rbf rates for OP-data, when ionizing
! to excited levels, but resonances below corresponding edge is present
! (analogous to INTERMEDR)

! ----  OPAC variable EFFECTIVE opacity. See comments.

! ----  comment by JO (Jan. 2016):
! the scattering component of ALO (pure Thomson so far)
! might need to be updated regarding the line pseudo continuum
! (note that this refers to the 'old' approach, i.e.,
! it is an old problem)

! ----  remember: ALO=0 for cmf_all frequencies

  Use :: nlte_type
  Use :: nlte_dim
  Use :: fund_const
  Use :: princesa_var, Only: frecin
  Use :: nlte_var, Only: alphafs, fre1 => fref, wfre1 => wfref, indfin, &
    indfrq1, ifre1 => ifref, fre, opac, xj, alo, opat_m_old, strue_m_old, &
    strue_m_new, metals_converged, op_flag, frecfn1

  Use :: nlte_porvor, Only: opa_eff_rat, opa_eff_rat_old

  Implicit None

! ---  new algorithm
! ---  subroutine intemepho1 (see netmat). it calculates photo-integrals
! ---  for some the formulae implemented in detail.

! .. parameters ..
! Integer (i4b), Parameter :: nd1 = id_ndept
! Integer (i4b), Parameter :: ifretot = id_frec1
! Integer (i4b), Parameter :: ifrefin = id_frec2
! ..
! .. scalar arguments ..
  Real (dp) :: depart, temp, x1, x2, xne, edgenu, xnistar, clf
  Integer (i4b) :: ll, ntr
  Logical :: optmet, optmet1
! ..
! .. local scalars ..
  Real (dp) :: alovth, alpha, aux1, aux2, dex, dfreq, dfx, dr, exk, f, f1, f2, &
    facf, facg, fgconst, freq, freqii, freqii1, g, g1, g2, hcfre3, hckt, &
    hnue1, opa, pi4, sig, sum1, sum2, thom, thomalo, wf, xalo, xf, xf1, xf2, &
    xfreq, xfreq1, xg, xg1, xg2, xxj, betar, lastfreq
  Integer (i4b) :: ii, ii1, iiold, istart, kk
! ..
! .. external subroutines ..
  External :: crossbf
! ..
! .. intrinsic functions ..
! INTRINSIC ABS,EXP,LOG10,MAX
! ..

  optmet1 = optmet
  If (optmet1 .And. metals_converged) optmet1 = .False.

  hckt = hkl/temp

! begin of integration of the two photo-integrals
! (for definition, see the santolaya et al. paper)
! note: integration over frequency, fre in kayser

  iiold = 0
! not used: rerr = 1.D-5
  sig = amh*sigmae
  istart = indfrq1(ntr)
  If (istart<=0) Then
    Print *, 'TRANSITION', ntr, ' AT', ll
    Print *, 'INDFRQ1 NOT DEFINED: CHECK ILOW OR IMAX (LTE AND NLTE!!!)'
  End If

  sum1 = 0.D0
  sum2 = 0.D0

! ---    note: first integration weight has to be modified

  kk = istart
  wf = wfre1(kk) - .5D0*clight*(fre1(kk)-fre1(kk-1))
  freq = fre1(kk)

! ---    ii index corresponding to coarse freq. grid

  edgenu = freq !                   maup
  lastfreq = frecin(ntr)

! some tests
  If (ll==1) Then
    If (.Not. op_flag(ntr)) Then
      Write (999, *) ' STOP: INTERMEPHO1 CALLED AND OP_FLAG = F'
      Stop ' INTERMEPHO1 CALLED AND OP_FLAG = F'
    End If
    If (edgenu>=lastfreq .Or. edgenu<frecfn1(ntr)) Then
      Write (999, *) ' STOP: SOMETHING ROTTEN WITH EDGENU -- INTERMEPHO1!'
      Stop ' SOMETHING ROTTEN WITH EDGENU -- INTERMEPHO1!'
    End If
  End If

  ii = indfin(kk)

  If (ii>0.) Then
    alpha = alphafs(ntr, ii)
  Else
    Write (999, *) ' STOP: II MUST NOT BE LT ZERO IN INDFIN AT &
      &FIRST FREQ. POINT (INTERMEPHO1)'
    Stop &
      ' II MUST NOT BE LT ZERO IN INDFIN AT FIRST FREQ. POINT (INTERMEPHO1)'
  End If


! OPA = OPAC(LL,II)*CLF !corrected
  opa = opac(ll, ii)*clf/opa_eff_rat_old(ll, ii) ! back-corrected
! ---- since OPAC is already an effective quantity, back-corrected here
! (and thus also below).
! NOTE: Use OLD values, since opat_m_old is used below
  thom = xne*sig/opa
  thomalo = 1.D0 - alo(ll, ii)*thom
  alovth = alo(ll, ii)/thomalo
  xxj = xj(ll, ii)
  If (optmet1) Then
    betar = opat_m_old(ll, ii)*clf/opa ! corrected (OPAT_M_OLD is mean
!   opacity)
    xxj = xxj + alovth*betar*(strue_m_new(ll,ii)-strue_m_old(ll,ii))
  End If

  exk = exp(-hckt*freq)
  hcfre3 = hc2*freq**3
  hnue1 = 2.D0/(hc2*freq)

! ---    f corresponds to i3-i4, g to i1-i2

  f = alpha*hnue1*(xxj-hcfre3*exk*alovth*alpha*xnistar/opa)
  g = alpha*hnue1*exk*(xxj+hcfre3*(1.D0-alovth*alpha*depart*xnistar/opa))

  sum1 = sum1 + f*wf
  sum2 = sum2 + g*wf

! ----   rest of integration, with all necessary interpolations

kkloop: Do kk = istart + 1, ifre1

    wf = wfre1(kk)
    freq = fre1(kk)
    If (freq>=lastfreq) Exit kkloop
    ii = indfin(kk)

    If (ii>0.) Then
      alpha = alphafs(ntr, ii)
!     OPA = OPAC(LL,II)*CLF !corrected
      opa = opac(ll, ii)*clf/opa_eff_rat_old(ll, ii) ! back-corrected
      thom = xne*sig/opa
      xalo = alo(ll, ii)
      thomalo = 1.D0 - xalo*thom
      alovth = xalo/thomalo
      xxj = xj(ll, ii)
!     not used: xxjj = xxj
      If (optmet1) Then
        betar = opat_m_old(ll, ii)*clf/opa ! corrected
        xxj = xxj + alovth*betar*(strue_m_new(ll,ii)-strue_m_old(ll,ii))
      End If
      exk = exp(-hckt*freq)
      hcfre3 = hc2*freq**3
      hnue1 = 2.D0/(hc2*freq)
      f = alpha*hnue1*(xxj-hcfre3*exk*alovth*alpha*xnistar/opa)
      g = alpha*hnue1*exk*(xxj+hcfre3*(1.D0-alovth*alpha*depart*xnistar/opa))

      If (f<0. .Or. g<0.) Then
        Print *, ntr, fre1(istart), kk, f, g, xxj, alovth, depart
        Print *, ii, freq, opa, xnistar
        Print *, alphafs(ntr, ii)
        Print *, 'EFFECTIVE OPACITY ISSUE?:', opa_eff_rat(ll, ii), &
          opa_eff_rat_old(ll, ii)
        Print *, opac(ll, ii), clf
        Print *, 'NEGATIVE CONTRIBUTION TO PHOTO-RATES1'
        Write (999, *) ' STOP!'
        Stop
      End If

    Else If (ii<0) Then

!     logarithmic interpolation for desired quantities

      xf = log10(freq)
      ii = -ii
      If (ii==iiold) Go To 100
      iiold = ii
      ii1 = ii + 1
!     OPA = OPAC(LL,II)*CLF !corrected
      opa = opac(ll, ii)*clf/opa_eff_rat_old(ll, ii) ! back-corrected
      thom = xne*sig/opa
      xalo = alo(ll, ii)
      thomalo = 1.D0 - xalo*thom
      alovth = xalo/thomalo
      xxj = xj(ll, ii)
!     not used: xxjj = xxj
      If (optmet1) Then
        betar = opat_m_old(ll, ii)*clf/opa ! corrected
        xxj = xxj + alovth*betar*(strue_m_new(ll,ii)-strue_m_old(ll,ii))
      End If
      freqii = fre(ii)
      exk = exp(-hckt*freqii)
      hcfre3 = hc2*freqii**3

      alpha = alphafs(ntr, ii)
      aux1 = xxj/(hcfre3*alpha)
!     not used: aux11 = xxjj/(hcfre3*alpha)
      aux2 = alovth*xnistar/opa
!     new evalution scheme, to allow for very low EXK (x-ray regime)
!     F1 = AUX1/EXK - AUX2
      f1 = aux1 - aux2*exk
      g1 = 1.D0/alpha + aux1 - depart*aux2

!     OPA = OPAC(LL,II1)*CLF  !corrected
      opa = opac(ll, ii1)*clf/opa_eff_rat_old(ll, ii1) ! back-corrected
      thom = xne*sig/opa
      xalo = alo(ll, ii1)
      thomalo = 1.D0 - xalo*thom
      alovth = xalo/thomalo
      xxj = xj(ll, ii1)
!     not used:  xxjj = xxj
      If (optmet1) Then
        betar = opat_m_old(ll, ii1)*clf/opa ! corrected
        xxj = xxj + alovth*betar*(strue_m_new(ll,ii1)-strue_m_old(ll,ii1))
      End If
      freqii1 = fre(ii1)
      exk = exp(-hckt*freqii1)
      hcfre3 = hc2*freqii1**3

      alpha = alphafs(ntr, ii1)
      aux1 = xxj/(hcfre3*alphafs(ntr,ii1))
      aux2 = alovth*xnistar/opa
!     F2 = AUX1/EXK - AUX2
      f2 = aux1 - aux2*exk
      g2 = 1.D0/alpha + aux1 - depart*aux2

      If (f1<0. .Or. g1<0. .Or. f2<0. .Or. g2<0.) Then
        Print *, ntr, fre1(istart), kk, f1, g1, f2, g2, xxj, alovth, depart
        Print *, ii, ii1, freqii, freqii1
        Print *, alphafs(ntr, ii), alphafs(ntr, ii1)
        Print *, 'effective issue?:', opa_eff_rat(ll, ii), &
          opa_eff_rat_old(ll, ii)
        Print *, opac(ll, ii), clf
        Print *, 'NEGATIVE CONTRIBUTION TO PHOTO-RATES1'
        Print *, 'NEGATIVE CONTRIBUTION TO PHOTO-RATES2'
        Write (999, *) ' STOP!'
        Stop
      End If

      xfreq = log10(freqii)
      xfreq1 = log10(freqii1)
!     XF1 = LOG10(F1)
!     XF2 = LOG10(F2)
!     ADD log10(exp(hnu/kt))
      xf1 = log10(f1) + hckt*freqii*log10e
      xf2 = log10(f2) + hckt*freqii1*log10e
      xg1 = log10(g1)
      xg2 = log10(g2)

      dfreq = 1.D0/(xfreq1-xfreq)
      facf = (xf2-xf1)*dfreq
      facg = (xg2-xg1)*dfreq

100   Continue

      dfx = xf - xfreq
      xf = xf1 + facf*dfx
      xg = xg1 + facg*dfx
!     XF = 10.D0**XF
!     XG = 10.D0**XG

      Call crossbf(ntr, alpha, freq)
      hnue1 = 2.D0/(hc2*freq)
      exk = exp(-hckt*freq)
      hcfre3 = hc2*freq**3
      fgconst = log10(alpha*hnue1*hcfre3*alpha) - hckt*freq*log10e
      f = 10.D0**(fgconst+xf)
      g = 10.D0**(fgconst+xg)

    Else
      Write (999, *) ' STOP: II MUST NOT BE ZERO IN INDFIN'
      Stop ' II MUST NOT BE ZERO IN INDFIN'
    End If

    x1 = f*wf
    x2 = g*wf

    sum1 = sum1 + x1
    sum2 = sum2 + x2

    dr = abs(x1/sum1)
    dr = max(dr, abs(x2/sum2))

    dex = 1.D0 - fre1(kk-1)/freq

!   IF (DEX.LT.1.D-7) STOP ' ERROR IN DEX'

!   The limit has been changed to allow for close energy terms
!   in detailed models for neutral ions (energy terms with almost
!   the same ionization energy). Note that this will increase the
!   computation time, but I would not expect a huge effect
    If (dex<1.D-9) Then
      Write (999, *) ' STOP: ERROR IN DEX'
      Stop ' ERROR IN DEX'
    End If
!   IF (DR.LT.RERR .AND. DEX.GT..01D0) EXIT
  End Do kkloop

  pi4 = 4.D0*pi
  x1 = pi4*sum1
  x2 = pi4*sum2

! print*,ntr,x1,' ',x2,' ',dr,' ',istart,' ',kk

  Return

End Subroutine


!-----------------------------------------------------------------------

Subroutine intermedr(i, ll, x1, x2, temp)

! CALCULATES DR RATES, XNUE=CM^-1

  Use :: nlte_type
  Use :: nlte_dim
  Use :: fund_const
  Use :: princesa_var, Only: data, iform4, indat4, numd4, kl, frecin, frecfn, &
    labl, labl4
  Use :: nlte_var, Only: ifre, fre, xj, optcmf_all, mlow, mup

  Implicit None

  Real (dp), Parameter :: const1 = e2mc2/hh*4.D0*pi**2

! .. scalar arguments ..
  Real (dp) :: temp, x1, x2
  Integer (i4b) :: i, ll

! .. local scalars ..
  Integer (i4b) :: n, ii1, ndata, ni, ntran, index, k, ml, mu, mu1
  Real (dp) :: xnue, flu, xlamc, blu, afac, ehnuekt, xxx, wilog, j1, j2, jnue
  Character (5) :: lab

  If (optcmf_all) Then
    If (xj(1,ifre+1)/=2) Then
      Write (999, *) &
        ' STOP: OPTCMF_ALL = TRUE, BUT APPROX. XJ (OBSFRAM) USED IN INTERMEDR'
      Stop ' OPTCMF_ALL = TRUE, BUT APPROX. XJ (OBSFRAM) USED IN INTERMEDR'
    End If
  End If

  n = iform4(i)
  If (n/=21) Then
    Write (999, *) ' STOP: WRONG FORMULA IN INTERMEDR'
    Stop ' WRONG FORMULA IN INTERMEDR'
  End If

  ii1 = indat4(i) + 3 !             OFFSET
  ndata = numd4(i)

  ntran = (ndata-3)/2

  x1 = 0.
  x2 = 0.

! LOOP OVER ALL RESONANCES
  k = 2

  Do ni = 1, ntran

    index = ii1 + 2*(ni-1)
    xnue = data(index)
!   --- for tests of Adi's DR data:
!   --- if ionization to excited states, first resonance should lie in between
!   --- ground-state (FRECFN) and excited state (FRECIN) edge.
    ml = mlow(i)
    mu = mup(i)
    mu1 = kl(ml)
    If (ni==1 .And. mu/=mu1) Then
      If (xnue<frecfn(i) .Or. xnue>frecin(i)) Then
        lab = labl(labl4(i))
!       JP May 2021: neglect small inconsistencies for OII in WM-Basic data
!       base
        If ((lab=='O2_5 ' .Or. lab=='O2_17' .Or. lab=='O2_15') .And. &
          xnue/frecfn(i)>0.95 .And. xnue<=frecin(i)) Then
          Continue
        Else
          Print *, labl(labl4(i)), xnue, frecfn(i), frecin(i)
!         after tests, this stop statement might be removed
          Write (999, *) ' STOP: PROBLEMS IN RESONANCES!!!'
          Stop ' PROBLEMS IN RESONANCES!!!'
        End If
      End If
    End If

    flu = data(index+1)

    If (xnue>fre(ifre) .Or. xnue<fre(1)) Cycle

    xlamc = 1.D0/xnue

    blu = const1*xlamc*flu
    afac = hc2/xlamc**3
    ehnuekt = exp(-hkl*xnue/temp)


!   ---------    log-log interpolation for JNUE

    k = k - 1
!   note: not all resonances in monotonic order
    If (k/=1.) Then
      If (xnue<fre(k-1)) k = 1
    End If

    xxx = fre(k)
    Do While (xnue>xxx)
      k = k + 1
      xxx = fre(k)
    End Do

    If (k>ifre) Then
      Write (999, *) ' STOP: INDEX NOT FOUND IN INTERMEDR'
      Stop ' INDEX NOT FOUND IN INTERMEDR'
    End If
!   can be commented after tested
    If (xnue>fre(k) .Or. xnue<fre(k-1)) Then
      Write (999, *) ' STOP: ERROR IN DR FRE-INDEX'
      Stop ' ERROR IN DR FRE-INDEX'
    End If

    wilog = log10(fre(k)/xnue)/log10(fre(k-1)/fre(k))

    j1 = xj(ll, k)
    j2 = xj(ll, k-1)
    jnue = j1*10.D0**(log10(j1/j2)*wilog)
!   -------------------------------------------

!   SUM UP RATES, ASSUMING BROAD RESONANCES (-> JNUE CAN BE USED)
    x1 = x1 + blu*jnue
    x2 = x2 + ehnuekt*blu*(afac+jnue)

  End Do

  Return
End Subroutine


!-----------------------------------------------------------------------

Subroutine linesob(ll, nfir, nlas, blev0, rat3, ndim, xnel, &
  nato, tempe, errold, qlu, qld)

! so far, this routine should be indendent on clumping, as long as
! all rates for energy equation remain unaffected


! new qualifiers for levels:
! iqual = 0 ok
! iqual = 1 problems, no high precision
! iqual = 2 negative occup. numbers achieved, reset to old value!


  Use :: nlte_type
  Use :: nlte_dim
  Use :: fund_const !               maup

  Use :: princesa_var, Only: le, data, index1, iform1, indat1, labat, labl1, &
    labu1

  Use :: nlte_var, Only: abundnd, alevel, blevel, mmlow, mmup, tlumat, tulmat, &
    indexrbb, iqual, unasol, megas, nist

  Use :: tcorr_var, Only: opttcor

  Implicit None

! .. parameters ..
  Integer (i4b), Parameter :: ndim1 = id_llevs
  Integer (i4b), Parameter :: nnrec = id_llevs
! Integer (i4b), Parameter :: kel = id_atoms, kis = id_kisat
! ..
! .. scalar arguments ..
  Real (dp) :: errold, tempe, xnel
  Integer (i4b) :: ll, nato, ndim, nfir, nlas
! ..
! .. array arguments ..
  Real (dp) :: blev0(nnrec), rat3(ndim1, ndim1)
  Real (dp), Dimension (nnrec, 2) :: qlu, qld

! ..
! .. local scalars ..
  Real (dp) :: bmax, casta1, casta2, clu, det, er, errmax, fiold, precis1, &
    tlu, tul, xlamc, xnl, xnlstr, xnu, xnustr, dtclu, hhnu, clu1
  Integer (i4b) :: i, ii, iii, indi, indmat, itm, j, jmax, m, ml, mn, mu, n1, &
    n2, nfor, nlo, nu
  Logical :: unados
! ..
! .. local arrays ..
  Real (dp) :: bold(ndim1), deltan(ndim1), erri(id_llevs), f(ndim1), &
    flud(ndim1), rat2(ndim1, ndim1), ratmat(ndim1, ndim1)

  Real (qp) :: aux(ndim1)

  Integer (i4b) :: indlud(ndim1), iq(id_llevs)
! ..
! .. external subroutines ..
  External :: lubksb, ludcmp, tratcol, tratrad
! ..
! .. intrinsic functions ..
! INTRINSIC ABS,INT,MAX,MIN
! ..

  If (ndim/=ndim1) Then
    Write (999, *) ' STOP: ERROR IN NDIM NE NDIM1'
    Stop ' ERROR IN NDIM NE NDIM1'
  End If

  precis1 = min(errold/10.D0, 1.D-3)
  precis1 = max(1.D-5, precis1)
  unados = unasol

  itm = 100

  mn = nlas - nfir + 1


  f(1:mn) = 0.D0
  iq(1:mn) = 0
  ratmat(1:mn, 1:mn) = rat3(1:mn, 1:mn)

iiloop: Do ii = 1, index1
    ml = mmlow(ii)
    mu = mmup(ii)
    If (ml<nfir .Or. ml>=nlas) Cycle
    If (mu>nlas) Cycle

    xnl = blev0(ml)
    xnu = blev0(mu)

    xnlstr = alevel(ml)
    xnustr = alevel(mu)
    If (nist(ml)==0.) Then
      Write (999, *) ' STOP: NIST(ML=0)'
      Stop ' NIST(ML=0)'
    End If
    If (nist(mu)==0.) Then
      Write (999, *) ' STOP: NIST(MU=0)'
      Stop ' NIST(MU=0)'
    End If
    indi = indat1(ii)
    er = data(indi)
    m = int(er)

    If (m==1) Then

!     ---  rbb transitions with occ.no. from previous iteration

      indmat = indexrbb(ii)
      tlu = tlumat(indmat, ll)
      tul = tulmat(indmat, ll)

      If (unados .And. megas) Then
        casta1 = xnl*tlu
        casta2 = xnu*tul
!       WRITE (100,FMT=9000) LL,'RBB',ML,MU,tlu,tul, &
        Write (100, Fmt=110) ll, 'RBB', ml, mu, casta1, casta2, &
          casta2 - casta1
!       WRITE (100,FMT=9000) LL,'RBB',ML,MU,CASTA1*NIST(ML)/XNL, &
!       &                CASTA2*NIST(MU)/XNU, (CASTA2 - CASTA1)/XNKK(ML)
      End If
    Else If (m==2) Then

      nlo = labl1(ii)
      nu = labu1(ii)
      xlamc = data(indi+1)*1.D-8
      nfor = iform1(ii)

      Call tratcol(ii, nlo, nu, xlamc, tempe, labl1, clu, nfor, indat1, &
        data, ll, 'A')
      tlu = clu*xnel

!     ---note: no correction for ne necessary, since transition in same ion

      tul = tlu*xnlstr/xnustr

      If (opttcor) Then
         Call tratcol(ii, nlo, nu, xlamc, tempe*1.01, labl1, clu1, nfor, indat1, &
           data, ll, 'B')
        dtclu = (clu1-clu)/(0.01*tempe)*xnel
!       There is no difference in the energies by using XLAMC or the
!       FL values
        hhnu = hh*clight/xlamc
!       DTCBBD(II) = XNU*(XNLSTR/XNUSTR)*DTCLU*HHNU !does not work
        qlu(ml, 1) = qlu(ml, 1) + tlu*hhnu ! (*x_l)
        qld(mu, 1) = qld(mu, 1) + tul*hhnu ! (*x_u)
        qlu(ml, 2) = qlu(ml, 2) + dtclu*hhnu ! (*x_l)

        If (labat(nato)=='H' .Or. labat(nato)=='HE') Then
!         standard treatment for H/He
          qld(mu, 2) = qld(mu, 2) + (xnlstr/xnustr)*dtclu*hhnu - &
            (tul*hhnu)*hhnu/akb/tempe/tempe ! (*x_u)
!         other elements
        Else !                      to be consistent with the background
!         treatment, see nlte_approx package
          If (abs(1.D0-tlu*blev0(ml)/(tul*blev0(mu)))<0.001) Then
            qld(mu, 2) = qld(mu, 2) + (xnlstr/xnustr)*dtclu*hhnu ! to allow
!           for
!           cancellation
          Else
            qld(mu, 2) = qld(mu, 2) + (xnlstr/xnustr)*dtclu*hhnu - &
              (tul*hhnu)*hhnu/akb/tempe/tempe ! (*x_u)
          End If
        End If
      End If

      If (unados .And. megas) Then
        casta1 = xnl*tlu
        casta2 = xnu*tul
!       WRITE (100,FMT=9000) LL,'CBB',ML,MU,tlu,tul, &
        Write (100, Fmt=110) ll, 'CBB', ml, mu, casta1, casta2, &
          casta2 - casta1
!       WRITE (100,FMT=9000) LL,'CBB',ML,MU,CASTA1*NIST(ML)/XNL, &
!       &                CASTA2*NIST(MU)/XNU, (CASTA2 - CASTA1)/XNKK(ML)
      End If

    Else

      Write (999, *) ' STOP: ERROR IN TRANSITION TYPE RBB/CBB'
      Stop ' ERROR IN TRANSITION TYPE RBB/CBB'
    End If

    n1 = ml - nfir + 1
    n2 = mu - nfir + 1

    If (n2>mn .Or. n1<=0) Then
      Write (999, *) ' STOP: ERROR IN LINESOB'
      Stop ' ERROR IN LINESOB'
    End If

!   ------- coeficients matrix

    ratmat(n1, n1) = ratmat(n1, n1) - tlu
    ratmat(n1, n2) = ratmat(n1, n2) + tul
    ratmat(n2, n1) = ratmat(n2, n1) + tlu
    ratmat(n2, n2) = ratmat(n2, n2) - tul


  End Do iiloop

  unados = .False.

! ---- newton-raphson for rate equations cont+lines

! ---- changed: leave out line with largest occ.num

  jmax = mn
  bmax = 0.

  Do i = 1, mn
    j = i + nfir - 1
    bold(i) = blev0(j)
    If (bold(i)>bmax) Then
      bmax = bold(i)
      jmax = i
    End If
  End Do

  Do i = 1, mn
    ratmat(jmax, i) = 1.D0
  End Do

  f(jmax) = abundnd(nato, ll)

  flud(1:mn) = f(1:mn)
  rat2(1:mn, 1:mn) = ratmat(1:mn, 1:mn)


! -----solution of rate equations

  Call ludcmp(ratmat, mn, ndim, indlud, det)

  Call lubksb(ratmat, mn, ndim, indlud, f)

! -----------------------------------------------------------------
! --------newton raphson improvement-cycle-------------------------
! -----------------------------------------------------------------

itloop: Do iii = 1, itm

!   old version (if no qp precision available)
!   DO I = 1,MN
!   SDP = -FLUD(I)

!   high precision summation

!   DO J=1,MN
!   AUX(J)=RAT2(I,J)*F(J)
!   END DO

!   CALL SORT1(MN,AUX)

!   DO J=1,MN
!   SDP=SDP+AUX(J)
!   END DO

!   DELTAN(I) = SDP
!   END DO

!   new version (residuum calculated in qp precision
    aux(1:mn) = matmul(real(rat2(1:mn,1:mn),qp), real(f(1:mn),qp)) - real(flud &
      (1:mn), qp)
    deltan(1:mn) = aux(1:mn)

    Call lubksb(ratmat, mn, ndim, indlud, deltan)

    errmax = 0.D0

    Do i = 1, mn

!     changed order (at first update, then error) to allow for f=0
!     after first iteration

      fiold = f(i)
      f(i) = f(i) - deltan(i)
      erri(i) = abs(deltan(i)/f(i))

      If (erri(i)<precis1) Then
        iq(i) = 0
      Else
        iq(i) = 1
      End If

      If (f(i)<0.D0) Then
        If (fiold<0.D0) Then
          f(i) = bold(i)
        Else
          f(i) = fiold
        End If
        iq(i) = 2
      End If

      errmax = max(errmax, erri(i))
    End Do

!   print*,'errmax: ',errmax

    If (errmax<precis1) Then

      Do i = 1, mn
        j = i + nfir - 1
        blevel(j) = f(i)

!       for this end condition, all levels have to be 0

        iqual(j, ll) = iq(i)
      End Do

      Go To 100
    End If

  End Do itloop

! at this point, we have problems with our occup. numbers
! only those errors are counted, which have "perfect" occup. numbers

  Print *, 'ATOM: ', nato, '   DEPTH: ', ll
  errmax = 0.D0

  Do i = 1, mn
    j = i + nfir - 1
    blevel(j) = f(i)
    iqual(j, ll) = iq(i)
    If (iq(i)==0) Then
      errmax = max(errmax, erri(i))
    Else
      Write (*, Fmt=120) j, iq(i)
    End If

  End Do

100 Continue

  Return

110 Format (1X, I2, 3X, A3, 2(3X,I3), 3(2X,G16.6))
120 Format (' LEVEL(ABS. NUMBER)', I3, ' WITH IQUAL = ', I1)

End Subroutine

!-----------------------------------------------------------------------

Subroutine tratcol(ii, nlo, nu, xlamc, tempe, labl1, clu, nfor, indat1, &
  data, ll, char)

! ---- subroutine for treating cbb transitions.
! ---- used (and partly checked so far)

! ------ formula 1:  si
! ------ formula 11: mg
! ------ formula 13/113: h i (either Giovanardi et al. or Burke et al.)
! ------ formula 14: he ii
! ------ formula 15: he i (singlett)
! ------ formula 16: he i (singlett, triplett)
! ------ formula 16: he i (singlett)
! ------ formula 17: he i (singlett, forbidden, intercomb.)
! ------ formula 19: he i (triplett, forbidden)
! ------ formula 20: he i (intercomb.)
! ------ formula 21: he i (intercomb.)
! ------ formula 22: mg (see formula 1) !------ formula 23: n iii
! ------ formula 24: si
! ------ formula 25: si, general
! ------ formula 26: c (bug identified and corrected, 05 Feb 2005)
! ------ formula 32: mg
! ------ formula 35: he i (singlett, forbidden)
! formula 37, 40: he i
! formula 45: h i
! formula 46: h i
! formula 51: for Adi's input, omega with slope
! formula 52: for Adi's input, omega = 1
! formula 60-65 : n iii (fine structure fits for gamma) from Stafford (1993)
! formula 66: niv (fine structure fits for gamma) from Ramsbottom (1994)

! formula 100: h i/ he ii (fits to exactly calculated cross-sections, &
! either from Butler (2004) or from Aggarwal (1991)

! JO July 2016
! formula 101: tabulated collision-strengths (e.g. CII/CIV from Aggarwal, CIII
! from Mitnik (2003)),
! in combination with formula 100. The first transition is read by formula 100
! (temperatures and coll. strengths), for the others formula 101 is used (only
! collision-strengths) when the T-grid is the same.
! NOTE!!!!! DIMENSION OF TGRID SET TO 50 (i.e., maximum number of grid-points
! = 50)
! NOTE!!!!! for reasons of comp. time, it's NOT checked whether NPT>50

! note: A10HHe(hydrogen) uses formula 13 for all transitions
! (from Giovanardi et al. for nl le 14, nu le 15,
! otherwise according to Burke
! A11HHe(hydrogen)
! uses formula 100 (Butler) for all trans. with nl, nu  le 7
! uses formula 113 (Burke)  for all trans. with nl le 4, nl ge 8
! uses formula  46 (P & R)  for all trans. with nl ge 5, nl ge 8

  Use :: nlte_type
  Use :: nlte_dim
  Use :: fund_const

  Use :: princesa_var, Only: fl, gl, labl, zl

  Implicit None
! .. parameters ..
  Integer (i4b), Parameter :: ngau = 6, nsp = 10
  Integer (i4b), Parameter :: nd1 = id_ndept
! ..
! .. scalar arguments ..
  Real (dp) :: clu, tempe, xlamc
  Integer (i4b) :: ii, nfor, nlo, nu, ll
  Character (1) :: char
! ..
! .. array arguments ..
  Real (dp) :: data(id_ndata)
  Integer (i4b) :: indat1(id_rcbbt), labl1(id_rcbbt)

  Real (dp) :: wlag(ngau), xlag(ngau), xsp(nsp)


  Real (dp), Dimension (nd1), Save :: tempa, tempb ! TEMPB CORRESPONDS TO
! 1.01*TEMPA
  Real (dp), Dimension (0:8, nd1), Save :: tsa, tsb
  Real (dp), Dimension (0:8) :: ts
  Real (dp), Dimension (0:23) :: ai
  Real (dp), Dimension (0:9) :: asp, asp2
  Real (dp), Dimension (50), Save :: tgrid ! for formula 100/101
! ..
! .. local scalars ..
  Real (dp) :: acon, acon1, adn, aifit, an, an2, au, au1, au2, au3, c0, c1, &
    c_1, c_2, e02, er, flu, gamma, gammae, gbar, omega, scale, t1, u0, usu1, &
    x, xlamio, eps1, eps2, fac, f, e1, e2, p1, p2, z, zzr, rydo, xnulo, xnuup, &
    ulo, uup, rn, n2, np2, rx, rx2, twon2, twon2rx, logtwnx, ann, bn, rfac, &
    bnn, con2, yex, yex2, ry, ry2, z0, zex, zex2, rz, rz2, upsfac, delta, &
    upsln, con1, te10d4, ee, tlog, slope, om, xgamma, lgamma1, lgamma2, &
    lgamma3, lgamma4, lgamma5, lgamma6, gamma1, gamma2, gamma3, gamma4, &
    gamma5, gamma6, ltempe, yp1, ypn

  Integer (i4b) :: ij, indi, j, lev, ndelta, nfit, nivl, nivu, nnlo, nnup, &
    npt, indx, ipt, k, i

  Integer (i4b), Save :: nptstore

  Character :: car*6, lab*6
  Logical :: found
! ..
! .. external functions ..
  Real (dp) :: erfcm, expin, expin1, hcol, pvr, pftn, prhcol
  External :: erfcm, expin, expin1, hcol, pvr, pftn, prhcol
! ..
! .. intrinsic functions ..
! INTRINSIC DBLE,EXP,INT,LOG,MAX,MIN,SQRT
! ..
  Data xlag, wlag/0.222846604D0, 1.1889321D0, 2.99273633D0, 5.77514357D0, &
    9.83746742D0, 15.982874D0, 0.458964674D0, 0.417000831D0, 0.113373382D0, &
    0.0103991975D0, 0.000261017203D0, 8.98547906D-07/

  Data xsp/3.50515, 3.90309, 4.20412, 4.50515, 4.90309, 5.20412, 5.50515, &
    5.90309, 6.20412, 6.50515/

  Data tempa/nd1*0./
  Data tempb/nd1*0./

  u0 = hh*clight/xlamc/akb/tempe
  acon1 = 5.465727564D-11
  xlamio = 9.117532285D-6

! ---- atom determination

  lev = labl1(ii)
! not used:  na = le(lev)
  indi = indat1(ii)

! --- prepare chebychev polynomials for formula 60 to 65


  If (nfor>=60 .And. nfor<=65) Then

    If (char=='A') Then
      If (tempe==tempa(ll)) Then
        ts = tsa(:, ll)
      Else
        tempa(ll) = tempe
        xgamma = (0.542519D0*log(tempe)) - 5.620738D0
        tsa(0, ll) = 1.D0
        tsa(1, ll) = xgamma
        Do i = 2, 8
          tsa(i, ll) = 2.D0*xgamma*tsa(i-1, ll) - tsa(i-2, ll)
        End Do
        ts = tsa(:, ll)
      End If

    Else If (char=='B') Then
      If (tempe==tempb(ll)) Then
        ts = tsb(:, ll)
      Else
        tempb(ll) = tempe
        xgamma = (0.542519D0*log(tempe)) - 5.620738D0
        tsb(0, ll) = 1.D0
        tsb(1, ll) = xgamma
        Do i = 2, 8
          tsb(i, ll) = 2.D0*xgamma*tsb(i-1, ll) - tsb(i-2, ll)
        End Do
        ts = tsb(:, ll)
      End If

    Else
      Write (999, *) ' STOP: WRONG CHAR IN TRATCOL'
      Stop ' WRONG CHAR IN TRATCOL'
    End If

  End If


! ---- principal quantum number, only for hi and heii

  If (nfor==13 .Or. nfor==113 .Or. nfor==14 .Or. nfor==45 .Or. nfor==46) Then
    lab = labl(nlo)
    j = 1

100 Continue

    car = lab(j:j)

    If (car=='0' .Or. car=='1' .Or. car=='2' .Or. car=='3' .Or. car=='4' .Or. &
      car=='5' .Or. car=='6' .Or. car=='7' .Or. car=='8' .Or. car=='9') Then
      j = j + 1
      ij = j

110   Continue

      car = lab(j:j)

      If (car=='0' .Or. car=='1' .Or. car=='2' .Or. car=='3' .Or. &
        car=='4' .Or. car=='5' .Or. car=='6' .Or. car=='7' .Or. car=='8' .Or. &
        car=='9') Then
        j = j + 1
        Go To 110
      Else
        j = j - 1
        car = lab(ij:j)
        Read (car, Fmt='(I2)') nivl
!       #             if(nivl.eq.0) nivl=10
        If (nivl==0) Then
          Write (999, *) ' STOP: NIVL = 0'
          Stop ' NIVL = 0'
        End If
      End If
    Else
      j = j + 1
      Go To 100
    End If

    lab = labl(nu)
    j = 1

120 Continue

    car = lab(j:j)

    If (car=='0' .Or. car=='1' .Or. car=='2' .Or. car=='3' .Or. car=='4' .Or. &
      car=='5' .Or. car=='6' .Or. car=='7' .Or. car=='8' .Or. car=='9') Then
      j = j + 1
      ij = j
130   Continue

      car = lab(j:j)

      If (car=='0' .Or. car=='1' .Or. car=='2' .Or. car=='3' .Or. &
        car=='4' .Or. car=='5' .Or. car=='6' .Or. car=='7' .Or. car=='8' .Or. &
        car=='9') Then
        j = j + 1
        Go To 130
      Else
        j = j - 1
        car = lab(ij:j)
        Read (car, Fmt='(I2)') nivu
!       #             if(nivu.eq.0) nivu=10
        If (nivu==0) Then
          Write (999, *) ' STOP: NIVU = 0'
          Stop ' NIVU = 0'
        End If
      End If
    Else
      j = j + 1
      Go To 120
    End If
  End If

! ---

  If (nfor==1 .Or. nfor==22) Then

!   ------- van regemorter

    flu = data(indi+2)
    gbar = data(indi+3)
    If (gbar/=.2D0 .And. gbar/=.7D0) Then
      Write (999, *) ' STOP: ERROR IN GBAR'
      Stop ' ERROR IN GBAR'
    End If

!   ------- note that expin1(u0) is not e_1(u0) but e_1(u0)*exp(u0)

    gammae = max(gbar, .276D0*expin1(u0))
    acon = acon1*(xlamc/xlamio)**2*flu*u0
    clu = 14.5D0*acon*sqrt(tempe)*exp(-u0)*gammae

    Return

  Else If (nfor==11) Then

!   warning: not checked by jp so far

    c_2 = data(indi+2)
    c_1 = data(indi+3)
    c0 = data(indi+4)
    c1 = data(indi+5)
    gamma = c_2/tempe**2 + c_1/tempe + c0 + c1*tempe
    clu = acon1*sqrt(tempe)*gamma*exp(-1.D0*u0)

!   print *
!   print *, ' >>> c_2  : ', c_2
!   print *, ' >>> c_1  : ', c_1
!   print *, ' >>> c0   : ', c0
!   print *, ' >>> c1   : ', c1
!   print *, ' >>> tempe: ', tempe
!   print *, ' >>> u0   : ', u0
!   print *, ' >>> xlamc: ', xlamc
!   print *
!   print *, ' >>> gamma: ', gamma
!   print *, ' >>> clu  : ', clu
!   print *

    Return

  Else If (nfor==13) Then

!   ------- hydrogen atom

    If (nivl>14 .Or. nivu>15) Then

      ndelta = nivu - nivl

!     ------- note that expin1(u0) is not e_1(u0) but e_1(u0)*exp(u0)
!     ------- idem with expin(5,u0)

      flu = data(indi+2)
      acon = acon1*(xlamc/xlamio)**2*flu*u0
      clu = acon*4.D0*(expin1(u0)+0.148D0*u0*expin(5,u0))*exp(-u0)
      If (ndelta/=1) Then
        au1 = 3.D0 - 1.2D0/dble(nivl)
        au2 = 1.8D0 - 0.4D0/dble(nivl)**2
        au3 = au1 + 2.D0*(au2-au1)/dble(ndelta)
        clu = clu*au3
      End If
      clu = clu*sqrt(tempe)

      Return
    Else
      clu = hcol(nivl, nivu, tempe)*exp(-u0)

!     here, we have included the correct temperature dependence,
!     hence:

      Return

    End If


  Else If (nfor==113) Then

!   ------- hydrogen atom, very old treatment (Burke et al. 1967)

    ndelta = nivu - nivl

!   ------- note that expin1(u0) is not e_1(u0) but e_1(u0)*exp(u0)
!   ------- idem with expin(5,u0)

    flu = data(indi+2)
    acon = acon1*(xlamc/xlamio)**2*flu*u0
    clu = acon*4.D0*(expin1(u0)+0.148D0*u0*expin(5,u0))*exp(-u0)
    If (ndelta/=1) Then
      au1 = 3.D0 - 1.2D0/dble(nivl)
      au2 = 1.8D0 - 0.4D0/dble(nivl)**2
      au3 = au1 + 2.D0*(au2-au1)/dble(ndelta)
      clu = clu*au3
    End If
    clu = clu*sqrt(tempe)

    Return

  Else If (nfor==14) Then

    flu = data(indi+2)
    acon = acon1*(xlamc/xlamio)**2*flu*u0
    clu = acon*exp(1.D0)*(log(2.D0)+expin1(u0))*exp(-u0)
    ndelta = nivu - nivl
    an = dble(nivl)
    adn = dble(ndelta)
    an2 = an - (an-1.D0)/adn
    au = min(adn, an2)
    If (nivl>1) au = au*1.1D0
    clu = clu*au*sqrt(tempe)

    Return

  Else If (nfor==15) Then

    flu = data(indi+2)
    acon = acon1*(xlamc/xlamio)**2*flu*u0
    clu = acon*4.D0*expin1(u0)*exp(-u0)*sqrt(tempe)

    Return

  Else If (nfor==16) Then

    e02 = .818730753D0
    flu = data(indi+2)
    acon = acon1*(xlamc/xlamio)**2*flu*u0
    usu1 = u0 + 0.2D0
    clu = acon*4.D0*(expin1(u0)*exp(-u0)-e02*u0/usu1*expin1(usu1)*exp(-usu1))
    clu = clu*sqrt(tempe)

    Return

  Else If (nfor==17) Then

    x = u0 + 1.D0/data(indi+5)
    clu = data(indi+3)*exp(-data(indi+4)*u0)*(x*data(indi+4)+2.D0)/x**3
    x = u0 + 1.D0/data(indi+8)
    clu = clu + data(indi+6)*exp(-data(indi+7)*u0)*(x*data(indi+7)+2.D0)/x**3
    clu = data(indi+2)*u0*(data(indi+9)*expin1(u0)*exp(-u0)+u0*clu)
    clu = clu*acon1*sqrt(tempe)

    Return

  Else If (nfor==19) Then

    clu = acon1*data(indi+2)*u0*(1.D0-u0*expin1(u0))*exp(-u0)*sqrt(tempe)

    Return

  Else If (nfor==20) Then

    er = erfcm(u0)*exp(-u0**2)
    clu = 2.D0/3.D0*(1.D0+u0)*exp(-u0) - sqpi*sqrt(u0)*er*(1.D0+2.D0/3.D0*u0)
    clu = clu*acon1*2.D0*data(indi+2)*u0*u0*sqrt(tempe)

    Return

  Else If (nfor==21) Then

    clu = acon1*data(indi+2)*u0*u0*(expin(2,u0)-expin(3,u0))*exp(-u0)* &
      sqrt(tempe)

    Return

  Else If (nfor==23) Then !         taken from DETAIL, not checked by  JO
!   Hummer distorted wave

    flu = data(indi+12)/data(indi+13)
    If (data(indi+3)>0.D0) flu = flu*data(indi+2)/data(indi+3)
    tlog = log10(tempe) - 4.D0
    gamma = data(indi+11)
    Do k = 1, 7
      j = indi + 11 - k
      gamma = gamma*tlog + data(j)
    End Do

    clu = gamma*exp(-1.D0*u0)*flu/sqrt(tempe)
    Return


  Else If (nfor==24) Then

!   ------- allen's cross section

    omega = data(indi+2)
    gamma = omega/gl(lev)*u0*(xlamc/xlamio)
    clu = acon1*sqrt(tempe)*exp(-u0)*gamma

    Return

  Else If (nfor==25) Then

!   ------- fit of effective collision strength

    t1 = data(indi+2)
    scale = data(indi+3)
    x = (tempe-t1)/scale
    nfit = int(data(indi+4))
    gamma = 0.D0

    Do j = 1, nfit
      aifit = data(indi+4+j)
      gamma = gamma + aifit*x**(j-1)
    End Do

    clu = 8.631D-6/(sqrt(tempe)*gl(lev))*exp(-u0)*gamma

    Return

  Else If (nfor==26) Then !         taken from DETAIL, not checked by  JO
!   powers of Log(T)

    t1 = data(indi+2)
    nfit = int(data(indi+4)+0.5D0)
    tlog = log10(tempe)
    gamma = 0.D0

    Do j = 1, nfit !                bug detected/corrected Feb 2005, miguel
      aifit = data(nfit-j+3+indi)
      gamma = aifit + gamma*tlog
    End Do

    clu = 8.631D-6/(sqrt(tempe)*gl(lev))*exp(-u0)*gamma

    Return

  Else If (nfor==32) Then !         taken from DETAIL, not checked by  JO

    acon = 4.D0*acon1*(109734.84D0*xlamc/clight)**2*data(indi+2)*u0
    clu = acon*sqrt(tempe)*(data(indi+3)*expin1(u0)+data(indi+4))
    clu = clu*exp(-1.D0*u0)
    Return


  Else If (nfor==35) Then

    nnlo = int(data(indi+2)+0.5D0)
    nnup = int(data(indi+3)+0.5D0)
!   CLU = HECOL(NNLO,NNUP,TEMPE)*EXP(-U0)  !old version for A7HHe
!   is replaced now by
    Call coll21(tempe, nnlo, nnup, gamma)
    clu = gamma*exp(-u0)

!   here, we have included the correct temperature dependence,
!   hence:

    Return

  Else If (nfor==37) Then


!   -------- helium i

    gamma = pvr(u0)
    acon = acon1*u0*(xlamc/xlamio)**2
    clu = acon*gamma*1.45D1*data(indi+2)*exp(-u0)*sqrt(tempe)

    Return

  Else If (nfor==40) Then

!   --------- helium i
    rydo = 3.29D15
    z = zl(nlo) + 1.D0
    zzr = z*z*rydo
    ulo = fl(nlo)/zzr
    uup = fl(nu)/zzr
    xnulo = 1.0D0/sqrt(ulo)
    xnuup = 1.0D0/sqrt(uup)

!   classical energy band

    eps1 = ulo - 1.0D0/((xnuup-0.5D0)**2)
    eps2 = ulo - 1.0D0/((xnuup+0.5D0)**2)
    p1 = pftn(tempe, ulo, eps1)
    p2 = pftn(tempe, ulo, eps2)
    e1 = exp(-1.579D5*eps1/tempe)
    e2 = exp(-1.579D5*eps2/tempe)
    If (eps1>0.0D0 .And. eps2>0.0D0) Then
      f = e2*p2 - e1*p1
    Else If (eps1<0.0D0 .And. eps2>0.0D0) Then
      f = e2*p2 - p1
    Else If (eps1<0.0D0 .And. eps2<0.0D0) Then
      f = p2 - p1
    Else
      Write (6, 140) tempe, eps1, eps2
      Write (999, *) ' STOP!'
      Stop
    End If
    clu = 8.69D-8*f*(1.579D5)**1.5D0
    fac = data(indi+2)*data(indi+2)*4.D0
    clu = clu/fac*gl(nu)/tempe/tempe*sqrt(tempe)

    Return

  Else If (nfor==45) Then
    Write (999, *) ' STOP: formula 45 not in use'
    Stop ' formula 45 not in use'

!   -------- hydrogen atom
!   -------- Johnson 1972

    rn = 1.D0/nivl
    n2 = nivl*nivl
    np2 = nivu*nivu
    x = 1.D0 - n2/np2
    rx = 1.D0/x
    rx2 = rx*rx
    twon2 = 2.D0*n2
    twon2rx = twon2*rx
    logtwnx = log(twon2rx)
    ann = data(indi+2)*twon2rx
    If (nivl==1) Then
      bn = -6.03D-1
      rfac = 4.5D-1
    Else
      bn = (((-2.809D1*rn+3.624D1)*rn-1.863D1)*rn+4.D0)*rn
      rfac = 1.94D0*rn**1.57D0
    End If
    rfac = rfac*x
    bnn = twon2*twon2*rx2/(np2*nivu)*(1.D0+4.D0/3.D0*rx+bn*rx2)
    con1 = 5.436D-11*twon2rx
    con2 = bnn - ann*logtwnx
    yex = expin1(u0)*exp(-u0)
    yex2 = expin(2, u0)*exp(-u0)
    ry = 1.D0/u0
    ry2 = ry + 0.5D0
    z0 = u0 + rfac
    zex = expin1(z0)*exp(-z0)
    zex2 = expin(2, z0)*exp(-z0)
    rz = 1.D0/z0
    rz2 = rz + 0.5D0
    clu = con1*u0*u0*sqrt(tempe)*(ann*(ry2*yex-rz2*zex)+con2*(ry*yex2-rz*zex2) &
      )


    Return

  Else If (nfor==46) Then

!   -------- hydrogen atom
!   -------- Percival & Richards 1978


    z = zl(nlo) + 1.D0
    acon = nivl*nivl*nivl*nivl*gl(nu)/(gl(nlo)*z*z)
    gamma = 0.D0
    upsfac = 8.631D-6/sqrt(tempe)
    Do j = 1, ngau
      ee = 6.33D-6*tempe*xlag(j)
      gamma = gamma + wlag(j)*prhcol(nivl, nivu, z, ee)
    End Do
    clu = gamma*acon*upsfac*exp(-u0)

    Return

  Else If (nfor==51) Then
!   omega and slope as input
    omega = data(indi+2)
    slope = data(indi+3)
    om = omega + (tempe-2.D4)/1.D4*slope ! FROM ADI
    om = max(omega/10., om)
    clu = 8.631D-6/(sqrt(tempe)*gl(lev))*exp(-u0)*om

    Return

  Else If (nfor==52) Then
!   omega = 1
    clu = 8.631D-6/(sqrt(tempe)*gl(lev))*exp(-u0)

    Return

!   New formulas included by Jorge
  Else If (nfor==60) Then

    Do i = 0, 3
      ai(i) = data(indi+2+i)
    End Do

    lgamma1 = 0.5D0*ai(0) + dot_product(ai(1:3), ts(1:3))

    gamma = exp(lgamma1)

    clu = 8.631D-6/(sqrt(tempe)*gl(lev))*exp(-u0)*gamma

    Return

  Else If (nfor==61) Then

    Do i = 0, 7
      ai(i) = data(indi+2+i)
    End Do

    lgamma1 = 0.5D0*ai(0) + dot_product(ai(1:3), ts(1:3))
    lgamma2 = 0.5D0*ai(4) + dot_product(ai(5:7), ts(1:3))

    gamma1 = exp(lgamma1)
    gamma2 = exp(lgamma2)

    gamma = gamma1 + gamma2

    clu = 8.631D-6/(sqrt(tempe)*gl(lev))*exp(-u0)*gamma

    Return

  Else If (nfor==62) Then

    Do i = 0, 11
      ai(i) = data(indi+2+i)
    End Do

    lgamma1 = 0.5D0*ai(0) + dot_product(ai(1:3), ts(1:3))
    lgamma2 = 0.5D0*ai(4) + dot_product(ai(5:7), ts(1:3))
    lgamma3 = 0.5D0*ai(8) + dot_product(ai(9:11), ts(1:3))

    gamma1 = exp(lgamma1)
    gamma2 = exp(lgamma2)
    gamma3 = exp(lgamma3)

    gamma = gamma1 + gamma2 + gamma3

    clu = 8.631D-6/(sqrt(tempe)*gl(lev))*exp(-u0)*gamma

    Return

  Else If (nfor==63) Then

    Do i = 0, 15
      ai(i) = data(indi+2+i)
    End Do

    lgamma1 = 0.5D0*ai(0) + dot_product(ai(1:3), ts(1:3))
    lgamma2 = 0.5D0*ai(4) + dot_product(ai(5:7), ts(1:3))
    lgamma3 = 0.5D0*ai(8) + dot_product(ai(9:11), ts(1:3))
    lgamma4 = 0.5D0*ai(12) + dot_product(ai(13:15), ts(1:3))

    gamma1 = exp(lgamma1)
    gamma2 = exp(lgamma2)
    gamma3 = exp(lgamma3)
    gamma4 = exp(lgamma4)

    gamma = gamma1 + gamma2 + gamma3 + gamma4

    clu = 8.631D-6/(sqrt(tempe)*gl(lev))*exp(-u0)*gamma

    Return

  Else If (nfor==64) Then

    Do i = 0, 23
      ai(i) = data(indi+2+i)
    End Do

    lgamma1 = 0.5D0*ai(0) + dot_product(ai(1:3), ts(1:3))
    lgamma2 = 0.5D0*ai(4) + dot_product(ai(5:7), ts(1:3))
    lgamma3 = 0.5D0*ai(8) + dot_product(ai(9:11), ts(1:3))
    lgamma4 = 0.5D0*ai(12) + dot_product(ai(13:15), ts(1:3))
    lgamma5 = 0.5D0*ai(16) + dot_product(ai(17:19), ts(1:3))
    lgamma6 = 0.5D0*ai(20) + dot_product(ai(21:23), ts(1:3))

    gamma1 = exp(lgamma1)
    gamma2 = exp(lgamma2)
    gamma3 = exp(lgamma3)
    gamma4 = exp(lgamma4)
    gamma5 = exp(lgamma5)
    gamma6 = exp(lgamma6)

    gamma = gamma1 + gamma2 + gamma3 + gamma4 + gamma5 + gamma6

    clu = 8.631D-6/(sqrt(tempe)*gl(lev))*exp(-u0)*gamma

    Return

  Else If (nfor==65) Then

    Do i = 0, 17
      ai(i) = data(indi+2+i)
    End Do


    lgamma1 = 0.5D0*ai(0) + dot_product(ai(1:8), ts(1:8))
    lgamma2 = 0.5D0*ai(9) + dot_product(ai(10:17), ts(1:8))

    gamma1 = exp(lgamma1)
    gamma2 = exp(lgamma2)

    gamma = gamma1 + gamma2

    clu = 8.631D-6/(sqrt(tempe)*gl(lev))*exp(-u0)*gamma

    Return

  Else If (nfor==66) Then

    Do i = 0, 9
      asp(i) = data(indi+2+i)
    End Do

    ltempe = log10(tempe)

    Call devp(xsp, asp, nsp, yp1, ypn)

    Call spline_col(xsp, asp, nsp, yp1, ypn, asp2)

    Call splinter(xsp, asp, asp2, nsp, ltempe, gamma)

    clu = 8.631D-6/(sqrt(tempe)*gl(lev))*exp(-u0)*gamma

    Return
!   ***********

  Else If (nfor==100) Then
!   fits to "exact" data, used for hydrogen (Butler, 2004) or helium II

    te10d4 = log10(tempe)

    upsfac = 8.631D-6/sqrt(tempe)

    npt = int(data(indi+2)+.5)

    found = .False.

    tgrid = 0.
    nptstore = npt
    tgrid(1:npt) = data(indi+3:indi+2+npt)

    If (tgrid(1)>te10d4) Then
      indx = 0
      delta = 0.D0
      found = .True.
    Else If (tgrid(npt)<te10d4) Then
      indx = npt
      delta = 0.D0
      found = .True.
    Else
      Do ipt = 1, npt - 1
        If ((.Not. found) .And. (tgrid(ipt)<te10d4) .And. (tgrid(ipt+ &
          1)>=te10d4)) Then
          found = .True.
          indx = ipt
          delta = te10d4 - tgrid(ipt)
          delta = delta/(tgrid(ipt+1)-tgrid(ipt))
        End If
      End Do
    End If

    If (indx==0) Then
      upsln = data(indi+2+npt+1)
    Else If (indx==npt) Then
      upsln = data(indi+2+npt+npt)
    Else
      upsln = data(indi+2+npt+indx+1)*delta
      upsln = upsln + (1.-delta)*data(indi+2+npt+indx)
    End If

    con1 = upsln/gl(nlo)

    clu = con1*upsfac*exp(-u0)

    Return
!   ***********

  Else If (nfor==101) Then
!   tabulated collision-strengths, in combination with formula 100
!   the temperature grid is stored from the last call to formula 100

    te10d4 = log10(tempe)

    upsfac = 8.631D-6/sqrt(tempe)

    npt = int(data(indi+2)+.5)
    If (npt/=nptstore) Then
      Print *, npt, nptstore
      Write (999, *) ' STOP: SOMETHING WRONG WITH FORMULA 101'
      Stop ' SOMETHING WRONG WITH FORMULA 101'
    End If

    found = .False.

    If (tgrid(1)>te10d4) Then
      indx = 0
      delta = 0.D0
      found = .True.
    Else If (tgrid(npt)<te10d4) Then
      indx = npt
      delta = 0.D0
      found = .True.
    Else
      Do ipt = 1, npt - 1
        If ((.Not. found) .And. (tgrid(ipt)<te10d4) .And. (tgrid(ipt+ &
          1)>=te10d4)) Then
          found = .True.
          indx = ipt
          delta = te10d4 - tgrid(ipt)
          delta = delta/(tgrid(ipt+1)-tgrid(ipt))
        End If
      End Do
    End If

    If (indx==0) Then
      upsln = data(indi+3)
    Else If (indx==npt) Then
      upsln = data(indi+2+npt)
    Else
      upsln = data(indi+2+indx+1)*delta
      upsln = upsln + (1.-delta)*data(indi+2+indx)
    End If

    con1 = upsln/gl(nlo)

    clu = con1*upsfac*exp(-u0)

    Return
!   ***********

  Else

    Print *, 'FORMULA NOT IMPLEMENTED YET'
    Print *, nfor
    Write (999, *) ' STOP!'
    Stop

  End If
140 Format (' TRATCOL : EPS CONDITION FAILS FOR TEM=', E10.3, '  EPS1=', &
    E10.3, '  EPS2=', E10.3, /)

End Subroutine

!----------------------------------------------------------------------- &

Function prhcol(n, np, z, e)
  Use :: nlte_type
  Use :: nlte_dim
  Implicit None
  Real (dp) :: prhcol
  Real (dp) :: l, z, e, e2, z2, twoznpe, rnp2, rs, fnnp, rfnnp, facn, a, d, f, &
    g, fac, xp, xm, y, h, eghtthd, twothd
  Integer (i4b) :: s, twsp1, n, np, n2
  Parameter (eghtthd=8.D0/3.D0, twothd=2.D0/3.D0)

  e2 = e*e
  z2 = z*z
  n2 = n*n
  twoznpe = 2.D0*z/(n2*e)
  rnp2 = 1.D0/(np*np)
  s = np - n
  twsp1 = 2*s + 1
  rs = 1.D0/s
  fnnp = np*n
  rfnnp = 1.D0/fnnp
  facn = np*rs/n
  a = eghtthd*rs*facn**3*(0.184D0-0.04D0*rs**twothd)* &
    (1.D0-0.20D0*s*rfnnp)**twsp1
  d = exp(-z2*rfnnp/e2)
  l = log((1.D0+0.53D0*e2*fnnp/z2)/(1.D0+0.4D0*e))
  f = (1.D0-0.3D0*s*d*rfnnp)**twsp1
  g = 0.5D0*(e*n2/(z*np))**3
  fac = sqrt(2.D0-n2*rnp2)
  xp = twoznpe/(fac+1.D0)
  xm = twoznpe/(fac-1.D0)
  y = 1.D0 - 0.25D0*d*log(18.D0*s)*rs
  y = 1.D0/y
  h = c2(xm, y) - c2(xp, y)
  prhcol = a*d*l + f*g*h
  Return
Contains
  Real (dp) Function c2(x, y)
    Real (dp) :: x, y

    c2 = x*x*log(1.D0+twothd*x)/(2.D0*y+1.5D0*x)
  End Function
End Function

!-----------------------------------------------------------------------

Subroutine tratrad(l, nlo, nu, flu, xlamc, xnl, xnu, xne, gl, velr, gvel, tlu, &
  tul, beto, sr, vmax, xmust, ii, clf, vtherm, fic, tcl_fac_line)

! UPDATED FOR CLUMPING

  Use :: nlte_type
  Use :: nlte_dim
  Use :: fund_const
  Use :: nlte_var, Only: as, a2s, a3s, bs, b2s, b3s, fre, ifre, opac, strue, &
    xj, xxk, opaclin, indexrbb, weig9 => w9, optlines, thomson_lines, &
    optcmf_all, restart_cmfall
! ,RAYLEIGH

  Use :: nlte_porvor, Only: opa_eff_rat_old

  Implicit None

! ---- radiative treatment. all elements needed for rate equation
! ---- coefficientes are calculated (in the framework of sobolev plus!
! ---- continuum absortion approach)

! ---- calculation of bcIc acclerated

! UPDATED FOR CLUMPING: present assumption: clumps small compared
! to Sobolev length; thus, "usual" averaging

! .. parameters ..
! Integer (i4b), Parameter :: nd1 = id_ndept
! Integer (i4b), Parameter :: ifretot = id_frec1

! maximum allowed tau_sob for inversion, to avoid too large corrections
  Real (dp), Parameter :: taumax = 2.3D0, sqthird = 0.5773502692D0, &
    prec = 1.D-8
  Real (dp), Parameter :: const1 = e2mc2/hh*4.D0*pi**2
! ..
! .. scalar arguments ..
  Real (dp) :: beto, flu, gvel, sr, tlu, tul, velr, vmax, vtherm, &
    xlamc, xmust, xne, xnl, xnu, clf
  Real (dp) :: fic, tcl_fac_line
  Integer (i4b) :: ii, l, nlo, nu
! ..
! .. array arguments ..
  Real (dp) :: gl(id_llevs)
! ..
! .. local scalars ..
  Real (dp) :: aul, bb1, bb2, bcic, betap, betas, betbig, betbim, betlem, &
    betles, blu, bul, delm1, delm2, dn, opacon, s1, s2, scont, taus, taus1, &
    taus2, ubar, wi1, wi2, wilog, xkli, xkline, xmue1, xmue2, xnue, xxx, &
    tausf, f1, f2, eddf, sl, jbar

  Real (dp) :: aas, aa2s, aa3s, bbs, bb2s, bb3s

  Real (dp) :: tcl, fac1

  Integer (i4b) :: i, indl, k
  Logical :: fastbeta, fastbcic, inversion, warning
! ..
! .. external functions ..
  Real (dp) :: betmue, betinv, qfunc, ups, u2_inv
  External :: betmue, betinv, qfunc, ups, u2_inv
! ..
! .. intrinsic functions ..
! INTRINSIC ABS,DBLE,LOG10
! ..

  xnue = 1.D0/xlamc

  If (optcmf_all) Then
    If (restart_cmfall) Then
      If (xj(1,ifre+1)/=1) Then
        Write (999, *) &
          ' STOP: RESTART_CMFALL, BUT XJ (CMF) INDICATED IN TRATRAD'
        Stop ' RESTART_CMFALL, BUT XJ (CMF) INDICATED IN TRATRAD'
      End If
!     at restart, all frequencies possible
    Else
      If (xj(1,ifre+1)/=2) Then
        Write (999, *) &
          ' STOP: OPTCMF_ALL = TRUE, BUT APPROX. XJ (OBSFRAM) USED IN TRATRAD'
        Stop ' OPTCMF_ALL = TRUE, BUT APPROX. XJ (OBSFRAM) USED IN TRATRAD'
      End If
!     JO Dec. 2016. In previous versions,
!     only evaluated at non-cmf_all frequencies
!     now: sobo transport also for lines inside range (if outside ilow,imax)
!     thus commented out
!     IF(XNUE.LT.FRE(KCMF_END) .AND. XNUE.GT.FRE(KCMF_START) )THEN
!     PRINT*,XLAMC
!     STOP ' OPTCMF_ALL = TRUE, BUT SOBO-LINE INSIDE CMF-RANGE (TRATRAD)'
!     ENDIF
    End If
  End If

  indl = indexrbb(ii)

! ---- einstein coefficients and line opacity
! ---- (warning! xlamc= central lambda in cm)


  If (xnue>fre(ifre) .Or. xnue<fre(1)) Then
    tlu = 0.D0
    tul = 0.D0
    Return
  End If

  blu = const1*xlamc*flu
  bul = gl(nlo)*blu/gl(nu)
  aul = hc2/xlamc**3*bul
  dn = (xnl/gl(nlo)-xnu/gl(nu))/clf ! corrected

  inversion = .False.
  If (dn<0.) inversion = .True.

  xkli = pi*e2mc2*clight*gl(nlo)*flu*dn
  xkline = xkli*xlamc*sr/vmax

  k = ifre
  xxx = fre(ifre)
  Do While (xnue<xxx)
    k = k - 1
    xxx = fre(k)
  End Do

! ---- as discussed by JS & JP, we use MAX(tcl_line,tcl_bg) =
! MIN(fac1,opa_eff_rat_old)
! (depth variable L, freq. variable K)
! ---- inversion consistent with the non-inversion case
  tcl = abs(xkline)*tcl_fac_line*vmax/sr
! ---- back-correct, since tcl_fac_line already includes SR/VMAX
  fac1 = (1.+fic*tcl)/(1.+tcl)
! IF (TCL.LT.0.0) FAC1=1.0D0   !old treatment of inversion (poor man's way..)
  fac1 = min(fac1, opa_eff_rat_old(l,k))
  If (fac1<=0.0D0 .Or. fac1>1.0D0) Then
    Print *, fac1, opa_eff_rat_old(l, k)
    Print *, tcl, l, k
    Write (999, *) ' STOP: OPA_EFF > <OPA>, TRATRAD_1'
    Stop ' OPA_EFF > <OPA>, TRATRAD_1'
  End If
! IF (XLAMC*1.d8.gt.1000.and.XLAMC*1.d8.lt.3000.) THEN
! Print*,'test-JS-TRATRAD:'
! Print*,L,K,XLAMC*1.d8
! Print*,XKLINE,TCL_FAC_LINE
! Print*,FAC1,OPA_EFF_RAT_OLD(L,K),L
! Print*,'-----------'
! ENDIF

  xkline = xkline*fac1
! now effective quantity

! changed: linear interpolation for generalized blocking factors,
! log-log interpolation for scont and opac

  wilog = log10(fre(k)/xnue)/log10(fre(k+1)/fre(k))
  wi1 = (fre(k+1)-xnue)/(fre(k+1)-fre(k))
  wi2 = 1.D0 - wi1

! at first, check whether tausf is in critical range
! to do so, calculate at first eddington factor

  f1 = xxk(l, k)/xj(l, k)
  f2 = xxk(l, k+1)/xj(l, k+1)
  eddf = f1*wi1 + f2*wi2

  tausf = xkline/qfunc(sqrt(eddf), velr, gvel)

  fastbcic = .False.
  fastbeta = .False.
  If (.Not. inversion) Then

!   slowbeta range has to include slowbcic range, otherwise betles to betbim
!   required for slowbeta are undefined

    If (tausf<=0.01D0 .Or. tausf>8.D0) fastbcic = .True.
    If (tausf<=0.0001D0 .Or. tausf>10.D0) fastbeta = .True.
  Else
    If (abs(tausf)>taumax) fastbcic = .True.
  End If

  If (fastbcic .Or. inversion) Then
    aas = xj(l, k)*10.D0**(log10(xj(l,k)/xj(l,k+1))*wilog)
  Else
    aas = as(k)*wi1 + as(k+1)*wi2
    bbs = bs(k)*wi1 + bs(k+1)*wi2
    aa2s = a2s(k)*wi1 + a2s(k+1)*wi2
    aa3s = a3s(k)*wi1 + a3s(k+1)*wi2
    bb2s = b2s(k)*wi1 + b2s(k+1)*wi2
    bb3s = b3s(k)*wi1 + b3s(k+1)*wi2
  End If

! ---- s_cont and beta_puls


  If (optlines) Then
!   S1 = STRUE(L,K) + XJ(L,K)* &
!   &
!   (XNE*SIGMAE*AMH/CLF+THOMSON_LINES(L,K)+RAYLEIGH(L,K)/CLF)/OPAC(L,K)
!   S2 = STRUE(L,K+1) + XJ(L,K+1)* &
!   &
!   (XNE*SIGMAE*AMH/CLF+THOMSON_LINES(L,K+1)+RAYLEIGH(L,K+1)/CLF)/OPAC(L,K+1)

!   ---- here only OPAC is effective, so we have to back-correct again
!   ---- BUT again OPAC HAS *NOT* been updated in this loop yet
!   ---- (i.e. NLTE_APPROX and OPACITM called, and so OPA_EFF_RAT updated,
!   but OPACITC NOT called yet) Thus OPA_EFF_RAT_OLD.
    s1 = strue(l, k) + xj(l, k)*(xne*sigmae*amh/clf+thomson_lines(l,k))/opac(l &
      , k)*opa_eff_rat_old(l, k) !  back-corrected
    s2 = strue(l, k+1) + xj(l, k+1)*(xne*sigmae*amh/clf+thomson_lines(l,k+1))/ &
      opac(l, k+1)*opa_eff_rat_old(l, k+1) ! back-corrected

  Else
    s1 = strue(l, k) + xne*sigmae*amh*xj(l, k)/clf/opac(l, k)*opa_eff_rat_old( &
      l, k) !                       back-corrected
    s2 = strue(l, k+1) + xne*sigmae*amh*xj(l, k+1)/clf/opac(l, k+1)* &
      opa_eff_rat_old(l, k+1) !     back-corrected
  End If

  scont = s1*10.D0**(log10(s1/s2)*wilog)

  opacon = opac(l, k)*10.D0**(log10(opac(l,k)/opac(l,k+1))*wilog) ! already OK
! ---- OPACON already effective
  opaclin(indl, l) = opacon

  If (.Not. inversion) Then

!   --------------------------------------------------------------------

!   this is the path for no inversion

    betap = opacon/xkline*sr*vtherm/vmax ! ratio remains unaffected
!   ---- same goes here -- both OPACON and XKLINE are already effective

    If (fastbeta) Then

!     beta evaluated at MUE^2 = 1/3

      taus = xkline/qfunc(sqthird, velr, gvel)
      betas = betmue(taus)
      ubar = ups(taus, betap)
    Else

!     ---- angular integration

      betles = 0.D0
      betbig = 0.D0
      betlem = 0.D0
      betbim = 0.D0
      ubar = 0.D0


!     ------- nine points for both intervals

      delm1 = xmust/8.D0
      delm2 = (1.D0-xmust)/8.D0
      Do i = 1, 9
        xmue1 = dble(i-1)*delm1
        xmue2 = xmust + dble(i-1)*delm2
        taus1 = xkline/qfunc(xmue1, velr, gvel)
        taus2 = xkline/qfunc(xmue2, velr, gvel)
        bb1 = betmue(taus1)*delm1
        bb2 = betmue(taus2)*delm2
        betles = betles + bb1*weig9(i)
        betbig = betbig + bb2*weig9(i)
        betlem = betlem + bb1*xmue1*weig9(i)
        betbim = betbim + bb2*xmue2*weig9(i)
        ubar = ubar + (ups(taus1,betap)*delm1+ups(taus2,betap)*delm2)*weig9(i)
      End Do

      betas = betles + betbig
    End If

    If (fastbcic) Then

!     betac_ic at MUE^2 = EDDF

!     JO: improved Sept 2023
      If (tausf<prec) Then
        bcic = (1.D0-0.5D0*tausf)*aas
      Else
        bcic = (1.D0-exp(-tausf))/tausf*aas
      End If
    Else

!     ---- more accurate integration: beta*i_inc integrated in mue and
!     beta_sob

      bcic = 0.5D0*(aas+aa3s)*betbig + .5D0*(bbs-bb3s)*betbim + &
        .5D0*(aa2s+aa3s)*betles + .5D0*(bb2s-bb3s)*betlem
    End If

    If (bcic<0.D0) Then
      If (bcic<-1.D-12) Then
        Print *, ' WARNING: BETAC IC NEGATIVE AT TRANSITION:', ii, xlamc*1.D8
        Write (999, *) ' WARNING: BETAC IC NEGATIVE AT TRANSITION:', ii, &
          xlamc*1.D8
      End If
      bcic = 0.
    End If

    beto = betas

!   ---- terms for rate equations

    tlu = blu*(bcic+scont*ubar)
    tul = bul*(bcic+scont*ubar) + aul*(betas+ubar)

!   if (xlamc.lt.305.d-8.and.xlamc.gt.303.d-8) then
!   print*,l,bcic,betas,betap,ubar,scont,tlu,tul
!   print*
!   endif

    Return

!   --------------------------------------------------------------------

  Else !                            path for inversion

!   note: provide always the possibility for a cancellation of Jnue and Scont

    warning = .False.
    xkline = abs(xkline) !          already corrected
    tausf = abs(tausf)

    betap = opacon/xkline*sr*vtherm/vmax

    If (fastbcic) Then

!     strong inversion, approx. treatment ok

      If (tausf<taumax) Then
        Write (999, *) ' STOP: ERROR IN FASTBCIC'
        Stop ' ERROR IN FASTBCIC'
      End If
      tausf = taumax
      warning = .True.
      betas = betinv(tausf)

!     ubar is now u2
      ubar = u2_inv(tausf, betap)

!     bcic now corresponds to Bc(Ic-Sc) = beta(sqrt(f))*(J-Sc)

    Else

!     ---- angular integration, J drawn out of integral

      betles = 0.D0
      betbig = 0.D0
      ubar = 0.D0

!     ------- nine points for both intervals

      delm1 = xmust/8.D0
      delm2 = (1.D0-xmust)/8.D0
      Do i = 1, 9
        xmue1 = dble(i-1)*delm1
        xmue2 = xmust + dble(i-1)*delm2
        taus1 = xkline/qfunc(xmue1, velr, gvel)
        If (taus1>taumax) Then
          taus1 = taumax
          warning = .True.
        End If
        taus2 = xkline/qfunc(xmue2, velr, gvel)
        If (taus2>taumax) Then
          taus2 = taumax
          warning = .True.
        End If
        bb1 = betinv(taus1)*delm1
        bb2 = betinv(taus2)*delm2
        betles = betles + bb1*weig9(i)
        betbig = betbig + bb2*weig9(i)
        ubar = ubar + (u2_inv(taus1,betap)*delm1+u2_inv(taus2,betap)*delm2)* &
          weig9(i)
      End Do

      betas = betles + betbig

    End If

    bcic = betas*(aas-scont)

    sl = -aul*xnu/(xnl*blu-xnu*bul) ! here, CLF cancels out (thin and thick
!   clumping)

    jbar = bcic + scont*(1.+ubar) + sl*ubar

    beto = betas

!   ---- terms for rate equations

    If (jbar>0.D0) Then

      tlu = blu*(bcic+scont*(1.D0+ubar))
      tul = bul*(bcic+scont*(1.D0+ubar)) + aul*(1.D0+ubar)

    Else !                          can only happen if Jnue < Scont
      Print *, 'JBAR NEGATIVE AT TRANSITION ', ii, xlamc*1.D8

!     jbar set to zero, use explicit rates

      tlu = 0.
      tul = aul
    End If

    If (warning) Then
      Print *, ' STRONG INVERSION AT', ii, xlamc*1.D8, l
!     occurs too often, thus only in output
!     WRITE(999,*) ' STRONG INVERSION AT',II,XLAMC*1.D8,L
    End If

    Return

  End If

End Subroutine

!-----------------------------------------------------------------------

Function betmue(taus)

  Use :: nlte_type
  Use :: nlte_dim
  Implicit None

! .. scalar arguments ..
  Real (dp) :: taus, betmue
! ..
! .. local scalars ..
  Real (dp) :: prec
! ..
! .. intrinsic functions ..
! INTRINSIC EXP
! ..

  prec = 1.D-8

  If (taus<prec) Then
    betmue = 1.D0 - 0.5D0*taus
  Else
    betmue = (1.D0-exp(-taus))/taus
  End If

  Return

End Function

!-----------------------------------------------------------------------

Function betinv(taus)

  Use :: nlte_type
  Use :: nlte_dim
  Implicit None

! .. scalar arguments ..
  Real (dp) :: taus, betinv
! ..
! .. local scalars ..
  Real (dp) :: prec
! ..
! .. intrinsic functions ..
! INTRINSIC EXP
! ..

  prec = 1.D-8

  If (taus<prec) Then
    betinv = 1.D0 + 0.5D0*taus
  Else
    betinv = (exp(taus)-1.D0)/taus
  End If

  Return

End Function

!-----------------------------------------------------------------------

Function qfunc(xmue, velr, gvel)

  Use :: nlte_type
  Use :: nlte_dim
  Implicit None

! .. scalar arguments ..
  Real (dp) :: gvel, velr, xmue, qfunc
! ..

  qfunc = velr*(1.D0-xmue**2) + gvel*xmue**2

  Return

End Function

!-----------------------------------------------------------------------

Function ups(tau, beta)

  Use :: nlte_type
  Use :: nlte_dim
  Use :: nlte_var, Only: ufung, ufunk
  Implicit None

! ---- approximation of function u (used in sobolev+continuum theory)
! ---- all beta's and tau's are possible.
! ---- e. santolaya, 1992 (inspired in gudrun taresch's)

! .. scalar arguments ..
  Real (dp) :: beta, tau, ups
! ..
! .. local scalars ..
  Real (dp) :: aux, b1, b2, betai, betal, betas, betau2, t1, t2, taui, taul, &
    taus, u, ui, uii, uis, us, usi, uss
  Integer (i4b) :: icb, ict, nb, ne, nn, nt
! ..
! .. local arrays ..
  Real (dp) :: betag(17), betak(76)
! ..
! .. external functions ..
  Real (dp) :: erfcm
  External :: erfcm
! ..
! .. intrinsic functions ..
! INTRINSIC DBLE,INT,LOG10,SQRT
! ..
! .. data statements ..
  Data betak/1.0D-10, 2.8D-10, 4.6D-10, 6.4D-10, 8.2D-10, 1.0D-9, 2.8D-9, &
    4.6D-9, 6.4D-9, 8.2D-9, 1.0D-8, 2.8D-8, 4.6D-8, 6.4D-8, 8.2D-8, 1.0D-7, &
    2.8D-7, 4.6D-7, 6.4D-7, 8.2D-7, 1.0D-6, 2.8D-6, 4.6D-6, 6.4D-6, 8.2D-6, &
    1.0D-5, 2.8D-5, 4.6D-5, 6.4D-5, 8.2D-5, 1.0D-4, 2.8D-4, 4.6D-4, 6.4D-4, &
    8.2D-4, 1.0D-3, 2.8D-3, 4.6D-3, 6.4D-3, 8.2D-3, 1.0D-2, 2.8D-2, 4.6D-2, &
    6.4D-2, 8.2D-2, 1.0D-1, 2.8D-1, 4.6D-1, 6.4D-1, 8.2D-1, 1.0D+0, 2.8D+0, &
    4.6D+0, 6.4D+0, 8.2D+0, 1.0D+1, 2.8D+1, 4.6D+1, 6.4D+1, 8.2D+1, 1.0D+2, &
    2.8D+2, 4.6D+2, 6.4D+2, 8.2D+2, 1.0D+3, 2.8D+3, 4.6D+3, 6.4D+3, 8.2D+3, &
    1.0D+4, 2.8D+4, 4.6D+4, 6.4D+4, 8.2D+4, 1.0D+5/


  Data betag/1.D-10, 1.D-9, 1.D-8, 1.D-7, 1.D-6, 1.D-5, 1.D-4, 1.D-3, 1.D-2, &
    1.D-1, 5.D-1, 1.D+0, 1.D+1, 1.D+2, 1.D+3, 1.D+4, 1.D+5/
! ..

  taul = log10(tau)

  betal = log10(beta)

! ---- casuistry begins...

  If (taul>2.D0) Then

    If (betal>5.D0) Then
      ups = 1.D0
      Go To 100
    End If

!   ------- interpolation indexes for tau and beta

    nt = int(taul) - 1
    nb = int(betal+10.D0) + 1

    If (beta>=0.5D0) nb = nb + 1

!   ------- we continue

    If (betal<-10.D0) Then

      If (taul<10.D0) Then
        taui = dble(nt) + 1.D0
        taus = taui + 1.D0
        ui = ufung(nt, 1)
        us = ufung(nt+1, 1)
        t1 = taul - taui
        t2 = taus - taul
        u = ui*t2 + us*t1
!       #               ups=u*1.d10*beta**2
        ups = u*beta
      Else
        ups = ufung(9, 1)*beta
      End If

    Else

!     ---------- weights for tau and beta (ufung region)

      taui = dble(nt) + 1.D0
      taus = taui + 1.D0

      If (betal==5.D0) Then
        icb = -1
      Else
        icb = 1
      End If

      betai = log10(betag(nb))
      betas = log10(betag(nb+icb))
      b1 = (betal-betai)/(betas-betai)
      b2 = (betas-betal)/(betas-betai)

!     ---------- cases again...

      If (taul>=10.D0) Then
        ui = ufung(9, nb)
        us = ufung(9, nb+icb)
        ups = ui*b2 + us*b1

        If (betal<0.D0) ups = ups*beta

      Else
        t1 = taul - taui
        t2 = taus - taul
        uii = ufung(nt, nb)
        uis = ufung(nt, nb+icb)
        usi = ufung(nt+1, nb)
        uss = ufung(nt+1, nb+icb)
        ups = uii*t2*b2 + uis*t2*b1 + usi*t1*b2 + uss*t1*b1

        If (betal<0.D0) ups = ups*beta

      End If
    End If

  Else

    If (taul>=-2.D0) Then

!     ---------- weights for tau and beta (ufunk region)

      nn = int(betal+10.D0) - 10
      aux = beta/10.D0**dble(nn)
      ne = int((aux-1.D0)/1.8D0) + 1
      nb = (nn+10)*5 + ne
      nt = int((taul+2.D0)/.2D0) + 1

      If (taul==2.D0) Then
        ict = 0
      Else
        ict = 1
      End If

      taui = (dble(nt)-1.D0)*.2D0 - 2.D0
      taus = taui + .2D0
      t1 = (taul-taui)/.2D0
      t2 = (taus-taul)/.2D0

!     ---------- we continue...

      If (betal<-10.D0) Then
        ui = ufunk(nt, 1)
        us = ufunk(nt+ict, 1)
        u = ui*t2 + us*t1
        ups = 10.D0**u
!       #               ups=ups*1.d10*beta**2
        ups = ups*beta

      Else If (betal<5.D0) Then
        betai = log10(betak(nb))
        betas = log10(betak(nb+1))
        b1 = (betal-betai)/(betas-betai)
        b2 = (betas-betal)/(betas-betai)
        uii = ufunk(nt, nb)
        uis = ufunk(nt, nb+1)
        usi = ufunk(nt+ict, nb)
        uss = ufunk(nt+ict, nb+1)
        u = uii*t2*b2 + uis*t2*b1 + usi*t1*b2 + uss*t1*b1
        ups = 10.D0**u

        If (beta<1.D0) ups = ups*beta

      Else
        ui = ufunk(nt, 76)
        us = ufunk(nt+ict, 76)
        u = ui*t2 + us*t1
        ups = 10.D0**u
      End If

    Else If (betal<=0.D0) Then

      ups = beta*.398942D0*tau**2

    Else If (beta*tau<10.D0) Then

      betau2 = beta*tau/sqrt(2.D0)
      ups = tau/2.D0*(1.D0-erfcm(betau2))

    Else

      ups = tau/2.D0 - .398942D0/beta
    End If

  End If

! ---- we finish!

100 Continue

End Function

!-----------------------------------------------------------------------

Function u2_inv(tau, beta)

! calculates u2 function for inverted levels;
! tau and beta are absolute values

  Use :: nlte_type
  Use :: nlte_var, Only: n1_b, n2_b, n1_s, n2_s, n1_l1, n2_l1, n1_l2, n2_l2, &
    u2b, u2s, u2l1, u2l2
  Implicit None

  Real (dp) :: u2_inv


! note, that the following number is not exactly sqrt(2 pi), &
! however the numerical analogon used here to be consistent with the
! numerical integrations performed to obtain u2.
! don't change it!

! .. parameters ..
  Real (dp), Parameter :: rdpi = 0.39908365D0
! ..
! .. scalar arguments ..
  Real (dp) :: beta, tau
! ..
! ..
! .. local scalars ..
  Real (dp) :: b1, b2, betax, bxl, t1, t2, taux, txl, u2, xb, xx, xxx, y1, y2, &
    y3, y4, yb, yy, yyy
  Integer (i4b) :: ib, ibl, ibl2, ibl2p1, iblp1, ibp1, ibs, ibsp1, it, itl, &
    itl2, itl2p1, itlp1, itp1, its, itsp1
! ..
! .. local arrays ..
  Real (dp), Dimension (n2_b) :: betab
  Real (dp), Dimension (n2_l1) :: betal1
  Real (dp), Dimension (n2_l2) :: betal2
  Real (dp), Dimension (n2_s) :: betas

  Real (dp), Dimension (n1_b) :: taub
  Real (dp), Dimension (n1_l1) :: taul1
  Real (dp), Dimension (n1_l2) :: taul2
  Real (dp), Dimension (n1_s) :: taus

! ..
! .. external functions ..
  Real (dp) :: dinte, zinter
  External :: dinte, zinter
! ..
! .. intrinsic functions ..
! INTRINSIC EXP, INT, LOG10
! ..
! .. data statements ..
  Data taub/.10000000000000D-05, .17782794100389D-05, .31622776601684D-05, &
    .56234132519035D-05, .10000000000000D-04, .17782794100389D-04, &
    .31622776601684D-04, .56234132519035D-04, .10000000000000D-03, &
    .17782794100389D-03, .31622776601684D-03, .56234132519035D-03, &
    .10000000000000D-02, .17782794100389D-02, .31622776601684D-02, &
    .56234132519035D-02, .10000000000000D-01, .17782794100389D-01, &
    .31622776601684D-01, .56234132519035D-01, .10000000000000D+00/
  Data betab/.1000D-05, .1000D-04, .1000D-03, .1000D-02, .1000D-01, .1000D+00/

  Data taus/ -.100D+01, -.980D+00, -.960D+00, -.940D+00, -.920D+00, -.900D+00, &
    -.880D+00, -.860D+00, -.840D+00, -.820D+00, -.800D+00, -.780D+00, &
    -.760D+00, -.740D+00, -.720D+00, -.700D+00, -.680D+00, -.660D+00, &
    -.640D+00, -.620D+00, -.600D+00, -.580D+00, -.560D+00, -.540D+00, &
    -.520D+00, -.500D+00, -.480D+00, -.460D+00, -.440D+00, -.420D+00, &
    -.400D+00, -.380D+00, -.360D+00, -.340D+00, -.320D+00, -.300D+00, &
    -.280D+00, -.260D+00, -.240D+00, -.220D+00, -.200D+00, -.180D+00, &
    -.160D+00, -.140D+00, -.120D+00, -.100D+00, -.800D-01, -.600D-01, &
    -.400D-01, -.200D-01, 0.000D-00, 0.200D-01, 0.400D-01, 0.600D-01, &
    0.800D-01, 0.100D+00, 0.120D+00, 0.140D+00, 0.160D+00, 0.180D+00, &
    0.200D+00, 0.220D+00, 0.240D+00, 0.260D+00, 0.280D+00, 0.300D+00, &
    0.320D+00, 0.340D+00, 0.360D+00, 0.380D+00, 0.400D+00, 0.420D+00, &
    0.440D+00, 0.460D+00, 0.480D+00, 0.500D+00, 0.520D+00, 0.540D+00, &
    0.560D+00, 0.580D+00, 0.600D+00, 0.620D+00, 0.640D+00, 0.660D+00, &
    0.680D+00, 0.700D+00, 0.720D+00, 0.740D+00, 0.760D+00, 0.780D+00, &
    0.800D+00, 0.820D+00, 0.840D+00, 0.860D+00, 0.880D+00, 0.900D+00, &
    0.920D+00, 0.940D+00, 0.960D+00, 0.980D+00, 0.100D+01, 0.102D+01, &
    0.104D+01, 0.106D+01, 0.108D+01, 0.110D+01, 0.112D+01, 0.114D+01, &
    0.116D+01, 0.118D+01, 0.120D+01, 0.122D+01, 0.124D+01, 0.126D+01, &
    0.128D+01, 0.130D+01, 0.132D+01, 0.134D+01, 0.136D+01, 0.138D+01, &
    0.140D+01, 0.142D+01, 0.144D+01, 0.146D+01, 0.148D+01, 0.150D+01, &
    0.152D+01, 0.154D+01, 0.156D+01, 0.158D+01, 0.160D+01, 0.162D+01, &
    0.164D+01, 0.166D+01, 0.168D+01, 0.170D+01, 0.172D+01, 0.174D+01, &
    0.176D+01, 0.178D+01, 0.180D+01, 0.182D+01, 0.184D+01, 0.186D+01, &
    0.188D+01, 0.190D+01, 0.192D+01, 0.194D+01, 0.196D+01, 0.198D+01, &
    0.200D+01, 0.202D+01, 0.204D+01, 0.206D+01, 0.208D+01, 0.210D+01, &
    0.212D+01, 0.214D+01, 0.216D+01, 0.218D+01, 0.220D+01, 0.222D+01, &
    0.224D+01, 0.226D+01, 0.228D+01, 0.230D+01, 0.232D+01, 0.234D+01, &
    0.236D+01, 0.238D+01, 0.240D+01, 0.242D+01, 0.244D+01, 0.246D+01, &
    0.248D+01, 0.250D+01, 0.252D+01, 0.254D+01, 0.256D+01, 0.258D+01, &
    0.260D+01, 0.262D+01, 0.264D+01, 0.266D+01, 0.268D+01, 0.270D+01, &
    0.272D+01, 0.274D+01, 0.276D+01, 0.278D+01, 0.280D+01, 0.282D+01/
  Data betas/ -.300D+01, -.298D+01, -.296D+01, -.294D+01, -.292D+01, &
    -.290D+01, -.288D+01, -.286D+01, -.284D+01, -.282D+01, -.280D+01, &
    -.278D+01, -.276D+01, -.274D+01, -.272D+01, -.270D+01, -.268D+01, &
    -.266D+01, -.264D+01, -.262D+01, -.260D+01, -.258D+01, -.256D+01, &
    -.254D+01, -.252D+01, -.250D+01, -.248D+01, -.246D+01, -.244D+01, &
    -.242D+01, -.240D+01, -.238D+01, -.236D+01, -.234D+01, -.232D+01, &
    -.230D+01, -.228D+01, -.226D+01, -.224D+01, -.222D+01, -.220D+01, &
    -.218D+01, -.216D+01, -.214D+01, -.212D+01, -.210D+01, -.208D+01, &
    -.206D+01, -.204D+01, -.202D+01, -.200D+01, -.198D+01, -.196D+01, &
    -.194D+01, -.192D+01, -.190D+01, -.188D+01, -.186D+01, -.184D+01, &
    -.182D+01, -.180D+01, -.178D+01, -.176D+01, -.174D+01, -.172D+01, &
    -.170D+01, -.168D+01, -.166D+01, -.164D+01, -.162D+01, -.160D+01, &
    -.158D+01, -.156D+01, -.154D+01, -.152D+01, -.150D+01, -.148D+01, &
    -.146D+01, -.144D+01, -.142D+01, -.140D+01, -.138D+01, -.136D+01, &
    -.134D+01, -.132D+01, -.130D+01, -.128D+01, -.126D+01, -.124D+01, &
    -.122D+01, -.120D+01, -.118D+01, -.116D+01, -.114D+01, -.112D+01, &
    -.110D+01, -.108D+01, -.106D+01, -.104D+01, -.102D+01, -.100D+01, &
    -.980D+00, -.960D+00, -.940D+00, -.920D+00, -.900D+00, -.880D+00, &
    -.860D+00, -.840D+00, -.820D+00, -.800D+00, -.780D+00, -.760D+00, &
    -.740D+00, -.720D+00, -.700D+00, -.680D+00, -.660D+00, -.640D+00, &
    -.620D+00, -.600D+00, -.580D+00, -.560D+00, -.540D+00, -.520D+00, &
    -.500D+00, -.480D+00, -.460D+00, -.440D+00, -.420D+00, -.400D+00, &
    -.380D+00, -.360D+00, -.340D+00, -.320D+00, -.300D+00, -.280D+00, &
    -.260D+00, -.240D+00, -.220D+00, -.200D+00, -.180D+00, -.160D+00, &
    -.140D+00, -.120D+00, -.100D+00, -.800D-01, -.600D-01, -.400D-01, &
    -.200D-01, 0.000D-00, 0.200D-01, 0.400D-01, 0.600D-01, 0.800D-01, &
    0.100D+00, 0.120D+00, 0.140D+00, 0.160D+00, 0.180D+00, 0.200D+00, &
    0.220D+00, 0.240D+00, 0.260D+00, 0.280D+00, 0.300D+00, 0.320D+00, &
    0.340D+00, 0.360D+00, 0.380D+00, 0.400D+00, 0.420D+00, 0.440D+00, &
    0.460D+00, 0.480D+00, 0.500D+00, 0.520D+00, 0.540D+00, 0.560D+00, &
    0.580D+00, 0.600D+00, 0.620D+00, 0.640D+00, 0.660D+00, 0.680D+00, &
    0.700D+00, 0.720D+00, 0.740D+00, 0.760D+00, 0.780D+00, 0.800D+00, &
    0.820D+00, 0.840D+00, 0.860D+00, 0.880D+00, 0.900D+00, 0.920D+00, &
    0.940D+00, 0.960D+00, 0.980D+00/

  Data taul1/ -.600D+01, -.590D+01, -.580D+01, -.570D+01, -.560D+01, &
    -.550D+01, -.540D+01, -.530D+01, -.520D+01, -.510D+01, -.500D+01, &
    -.490D+01, -.480D+01, -.470D+01, -.460D+01, -.450D+01, -.440D+01, &
    -.430D+01, -.420D+01, -.410D+01, -.400D+01, -.390D+01, -.380D+01, &
    -.370D+01, -.360D+01, -.350D+01, -.340D+01, -.330D+01, -.320D+01, &
    -.310D+01, -.300D+01, -.290D+01, -.280D+01, -.270D+01, -.260D+01, &
    -.250D+01, -.240D+01, -.230D+01, -.220D+01, -.210D+01, -.200D+01, &
    -.190D+01, -.180D+01, -.170D+01, -.160D+01, -.150D+01, -.140D+01, &
    -.130D+01, -.120D+01, -.110D+01, -.100D+01, -.900D+00, -.800D+00, &
    -.700D+00, -.600D+00, -.500D+00, -.400D+00, -.300D+00, -.200D+00, &
    -.100D+00, 0.000D+00/
  Data betal1/ -.100D+01, -.900D+00, -.800D+00, -.700D+00, -.600D+00, &
    -.500D+00, -.400D+00, -.300D+00, -.200D+00, -.100D+00, 0.000D+00, &
    0.100D+00, 0.200D+00, 0.300D+00, 0.400D+00, 0.500D+00, 0.600D+00, &
    0.700D+00, 0.800D+00, 0.900D+00, 0.100D+01, 0.110D+01, 0.120D+01, &
    0.130D+01, 0.140D+01, 0.150D+01, 0.160D+01, 0.170D+01, 0.180D+01, &
    0.190D+01, 0.200D+01, 0.210D+01, 0.220D+01, 0.230D+01, 0.240D+01, &
    0.250D+01, 0.260D+01, 0.270D+01, 0.280D+01, 0.290D+01, 0.300D+01, &
    0.310D+01, 0.320D+01, 0.330D+01, 0.340D+01, 0.350D+01, 0.360D+01, &
    0.370D+01, 0.380D+01, 0.390D+01, 0.400D+01, 0.410D+01, 0.420D+01, &
    0.430D+01, 0.440D+01, 0.450D+01, 0.460D+01, 0.470D+01, 0.480D+01, &
    0.490D+01, 0.500D+01, 0.510D+01, 0.520D+01, 0.530D+01, 0.540D+01, &
    0.550D+01, 0.560D+01, 0.570D+01, 0.580D+01, 0.590D+01, 0.600D+01/

  Data taul2/0.0000D+00, 0.2500D+00, 0.5000D+00, 0.7500D+00, 1.0000D+00, &
    1.2500D+00, 1.5000D+00, 1.7500D+00, 2.0000D+00, 2.2500D+00, 2.5000D+00, &
    2.7500D+00/
  Data betal2/0.7500D+00, 1.0000D+00, 1.2500D+00, 1.5000D+00, 1.7500D+00, &
    2.0000D+00, 2.2500D+00, 2.5000D+00, 2.7500D+00, 3.0000D+00, 3.2500D+00, &
    3.5000D+00, 3.7500D+00, 4.0000D+00, 4.2500D+00, 4.5000D+00, 4.7500D+00, &
    5.0000D+00/
! ..
! r
! r    reset values of betax/taux if they can stop the program
! r
  taux = tau

  betax = beta
  If (betax<.1D-5) Then
    betax = .1D-05
  End If

! COMMENTED OUT, SINCE CAUGHT IN CALLING ROUTINE
! IF (TAUX.GT.0.65D3) THEN
! TAUX = 0.65D3
! ENDIF
! r

  txl = log10(taux)
  bxl = log10(betax)
! r
! r    set the values of u2 out of the table
! r    -------------------------------------

  If (txl<-.6D1) Then
    u2_inv = 0.5D0*taux

    Return
  End If

  If (bxl>.5D1 .And. txl>0.) Then
    u2_inv = rdpi/betax

    Return
  End If
! r
  If (bxl>.6D1) Then
    u2_inv = rdpi/betax

    Return
  End If

! r    "extrapolation" out of little2 (tau=last value in table)

  If (txl>=.275D1 .And. bxl>=.1D1) Then
    itl2 = n1_l2
    yyy = (bxl-1.D0)*4.D0 + 1.D0
    ibl2 = int(yyy)
    ibl2p1 = ibl2 + 1
    y1 = u2l2(itl2, ibl2)
    y2 = u2l2(itl2, ibl2p1)
    u2 = dinte(y1, y2, bxl, betal2(ibl2), betal2(ibl2p1))
    u2_inv = 10.D0**u2

    Return
  End If
! r
  If (bxl<=-.6D1) Then
    u2_inv = (exp(taux)-1.D0)/taux - 1.D0

    Return
  End If
! r
! r
  If (txl>.278D1 .And. bxl<-.1D1) Then
    u2_inv = (exp(taux)-1.D0)/taux - 1.D0

    Return
  End If
! r
! r    "extrapolation" out of small (same as above)
! r
  If (txl>.278D1 .And. bxl>-.1D1 .And. bxl<.1D1) Then
    yy = (bxl+3.D0)*50.D0 + 1.D0
    ibs = int(yy)
    its = n1_s
    ibsp1 = ibs + 1
    y1 = u2s(its, ibs)
    y2 = u2s(its, ibsp1)
    u2 = dinte(y1, y2, bxl, betas(ibs), betas(ibsp1))
    u2_inv = 10.D0**u2

    Return
  End If
! r
! r    interpolation regions
! r    --------------------------

! small table - interpolation of log(u2) in log(t) & log(b)
! range: -1>=log(tau)>=2.82 & -3>=log(beta)>=0.98
! r
  If (txl>=-.1D1 .And. bxl>=-.3D1 .And. bxl<=.98D0) Then
!   r
    xx = (txl+1.D0)*50.D0 + 1.D0
    its = int(xx)
    itsp1 = its + 1

    yy = (bxl+3.D0)*50.D0 + 1.D0
    ibs = int(yy)
    ibsp1 = ibs + 1
!   r
    y1 = u2s(its, ibsp1)
    y2 = u2s(itsp1, ibsp1)
    y3 = u2s(itsp1, ibs)
    y4 = u2s(its, ibs)
    t1 = taus(its)
    t2 = taus(itsp1)
    b1 = betas(ibs)
    b2 = betas(ibsp1)
    u2 = zinter(y1, y2, y3, y4, t1, t2, b1, b2, txl, bxl)
    u2_inv = 10.D0**u2
    Return

  End If
! r
! little1 table - interpolation of log(u2) in log(t) & log(b)
! range: -6>=log(tau)>=0 & -1>=log(beta)>=6

  If (bxl>=-1.D0 .And. txl<=0.D0) Then
!   r
    xx = (txl+6.D0)*10.D0 + 1.D0
    itl = int(xx)
    itlp1 = itl + 1

    yy = (bxl+1.D0)*10.D0 + 1.D0
    ibl = int(yy)
    iblp1 = ibl + 1
!   r
    y1 = u2l1(itl, iblp1)
    y2 = u2l1(itlp1, iblp1)
    y3 = u2l1(itlp1, ibl)
    y4 = u2l1(itl, ibl)
    t1 = taul1(itl)
    t2 = taul1(itlp1)
    b1 = betal1(ibl)
    b2 = betal1(iblp1)
    u2 = zinter(y1, y2, y3, y4, t1, t2, b1, b2, txl, bxl)
    u2_inv = 10.D0**u2
    Return

  End If
! r
! little2 table - interpolation of log(u2) in log(t) & log(b)
! range: 0>=log(tau)>=2.75 & 0.75>=log(beta)>=5

  If (bxl>0.75D0 .And. txl>0.D0) Then
!   r
    xxx = txl*4.D0 + 1.D0
    itl2 = int(xxx)
    itl2p1 = itl2 + 1

    yyy = (bxl-0.75D0)*4.D0 + 1.D0
    ibl2 = int(yyy)
    ibl2p1 = ibl2 + 1
!   r
    y1 = u2l2(itl2, ibl2p1)
    y2 = u2l2(itl2p1, ibl2p1)
    y3 = u2l2(itl2p1, ibl2)
    y4 = u2l2(itl2, ibl2)
    t1 = taul2(itl2)
    t2 = taul2(itl2p1)
    b1 = betal2(ibl2)
    b2 = betal2(ibl2p1)
    u2 = zinter(y1, y2, y3, y4, t1, t2, b1, b2, txl, bxl)
    u2_inv = 10.D0**u2
    Return

  End If
! r
! big table - interpolation of u2 in tau & beta
! range: 1d-6>=tau>=0.1 & 1d-6>=beta>=0.1

  If (txl<=-1.D0 .And. bxl<-1.D0) Then
    xb = (txl+6.D0)*4.D0 + 1.D0
    it = int(xb)

    itp1 = it + 1
!   r
    yb = bxl + 7.D0
    ib = int(yb)

    ibp1 = ib + 1

    y1 = u2b(it, ibp1)/(.5D0*taub(it))
    y2 = u2b(itp1, ibp1)/(.5D0*taub(itp1))
    y3 = u2b(itp1, ib)/(.5D0*taub(itp1))
    y4 = u2b(it, ib)/(.5D0*taub(it))
    t1 = taub(it)
    t2 = taub(itp1)
    b1 = betab(ib)
    b2 = betab(ibp1)
    u2 = zinter(y1, y2, y3, y4, t1, t2, b1, b2, taux, betax)
    u2_inv = .5D0*taux*u2

    Return
  End If

! region of possible big numbers
! good aprox. u2=(exp(tau)-1)/tau-1

  If (txl>-1.D0 .And. bxl<-3.D0) Then
    Print *, 'Possible big number, no interpolation'
    u2_inv = (exp(taux)-1.D0)/taux - 1.D0

    Return
  End If

End Function
!r
!r    ***************************************************************
!r    function zinter performs the interpolation
!r    ***************************************************************
Function zinter(y1, y2, y3, y4, t1, t2, b2, b1, taux, betax)
  Use :: nlte_type
  Implicit None

  Real (dp) :: zinter
! .. scalar arguments ..
  Real (dp) :: b1, b2, betax, t1, t2, taux, y1, y2, y3, y4
! ..
! .. local scalars ..
  Real (dp) :: t, u
! ..
  t = (taux-t1)/(t2-t1)
  u = (betax-b1)/(b2-b1)
  zinter = (1.D0-t)*(1.D0-u)*y1 + t*(1.D0-u)*y2 + t*u*y3 + (1.D0-t)*u*y4

  Return
End Function
!r
!r    ***************************************************************
!r    function dinte interpolates only beta values
!r    *************************************************************** &
Function dinte(y1, y2, betax, b, bp1)
  Use :: nlte_type
  Implicit None

  Real (dp) :: dinte
! .. scalar arguments ..
  Real (dp) :: b, betax, bp1, y1, y2
! ..
! .. local scalars ..
  Real (dp) :: u
! ..
  u = betax - b
  u = u/(bp1-b)
  dinte = u*(y2-y1) + y1

  Return
End Function

!-----------------------------------------------------------------------

Subroutine update(xnh, xne, te, ilow, imax, ll, nat)

  Use :: nlte_type
  Use :: nlte_dim

  Use :: princesa_var, Only: labat, zeff

  Use :: nlte_var, Only: enionnd, optneupdate

  Implicit None

! .. parameters ..
  Integer (i4b), Parameter :: kel = id_atoms ! , kis = id_kisat
  Integer (i4b), Parameter :: nd1 = id_ndept
! ..
! .. scalar arguments ..
  Integer (i4b) :: ll, nat
! ..
! .. array arguments ..
  Real (dp) :: te(nd1), xne(nd1), xnh(nd1)
  Integer (i4b) :: ilow(nd1, kel), imax(nd1, kel)
! ..
! .. local scalars ..
  Real (dp) :: aux, err, rel, xneh, xnenew, xneold, zi
  Integer (i4b) :: i, iz, k, nnn, nnuu
  Integer (i4b), Dimension (1) :: main
! ..
! .. intrinsic functions ..
! INTRINSIC ABS,INT
! ..

  xneold = xne(ll)
  xnenew = 0.D0

  Do k = 1, nat

    Do i = ilow(ll, k), imax(ll, k) + 1

      zi = zeff(k) + i - 1
      xnenew = xnenew + enionnd(k, i, ll)*zi
    End Do
!   check for consistency
    Do i = 1, ilow(ll, k) - 1
      If (enionnd(k,i,ll)/=0.) Then
        Write (999, *) ' STOP: ENIONND NE 0 FOR I < ILOW '
        Stop ' ENIONND NE 0 FOR I < ILOW '
      End If
    End Do

    Do i = imax(ll, k) + 2, id_kisat + 1
      If (enionnd(k,i,ll)/=0.) Then
        Write (999, *) ' STOP: ENIONND NE 0 FOR I > IMAX +1 '
        Stop ' ENIONND NE 0 FOR I > IMAX +1 '
      End If
    End Do
  End Do

  If (.Not. optneupdate) xnenew = xneold ! LEAVE NE AT ITS LTE VALUE
  xne(ll) = xnenew
! print*,'xnecheck: in update',optneupdate
! if(ll.eq.nd1) print*,xne


  rel = xnenew/xneold
  err = abs(1.D0-rel)

  xneh = xnenew/xnh(ll)
  nnn = ll/6
  nnuu = nnn*6
  If (ll==1 .Or. ll==47) nnuu = ll

  If (nnuu==ll) Then
    Print *
    Print *, &
      '----------------------------------------------------------------'
    Print *, '   L=', ll, '   TE = ', te(ll), ' NE/NH = ', xneh
    Print *
    Print *, ' ATOM // Z(MAIN ION) // N(J,K)/N(H)'

    Do k = 1, nat

      iz = int(zeff(k))
      iz = iz - 1
      aux = 0.D0
      main = maxloc(enionnd(k,1:imax(ll,k),ll))
      Write (*, Fmt=100) labat(k), iz + main(1), (aux, i=1, iz+ilow(ll,k)), &
        (enionnd(k,i,ll)/xnh(ll), i=ilow(ll,k), imax(ll,k)+1)
    End Do

    Print *
    Print *, ' CHANGE IN ELECTRON DENSITY = ', err
    Print *

  End If

  Return

100 Format (/, A4, '/', I3, 2X, '/', 9(E10.4,2X))

End Subroutine

!***********************************************************************

!subroutines: complex ones
!rateli and related (esp. cmf)

!***********************************************************************

Subroutine rateli(xne, temp, clf, ilow, imax, nd, concon, r, v, gradv, sr, &
  vmax, teff, corrfc)

! UPDATED FOR CLUMPING
! ILOW, IMAX, IMIANL AND IMAANL: NLTE VALUES ONLY!!!!

  Use :: nlte_type
  Use :: nlte_dim

  Use :: princesa_var, Only: nat, iong, ifirsl

  Use :: nlte_var, Only: blevel, indxlamc, xlamcmf, optcmf, optmixed, imianl, &
    imaanl, indexcmf, indexsa, modnam, lte_update

  Use :: nlte_porvor, Only: fic, tcl_fac_line, tcl_fac_cont

  Implicit None

! .. parameters ..
  Integer (i4b), Parameter :: kel = id_atoms ! , kis = id_kisat
  Integer (i4b), Parameter :: nd1 = id_ndept

! ..
! .. scalar arguments ..
  Real (dp) :: corrfc, sr, teff, vmax
  Integer (i4b) :: nd
  Logical :: concon
! ..
! .. array arguments ..
  Real (dp) :: gradv(nd1), r(nd1), temp(nd1), v(nd1), xne(nd1), clf(nd1)

  Integer (i4b) :: ilow(nd1, kel), imax(nd1, kel)
! ..
! .. local scalars ..
  Real (dp) :: velr, xlam, xmust

  Integer (i4b) :: i, icmf, icount, ii, ilo, ima, j, ll, lx, nato, nfi, nfir, &
    nfir1, nfirm, nla, nlas, nlas1, nlasm, nrec, numion
  Integer (i4b), Save :: ipres

  Logical :: startcmf
! ..
! .. external functions ..
  Integer (i4b) :: igenio

! .. data statements ..

  Data startcmf/.True./
! ..

! AFTER LTE-UPDATE, CREATE NEW CMF-LIST (CHANGE IN IMIN,IMAX)
  If (.Not. startcmf .And. lte_update) startcmf = .True.

  nrec = id_llevs

  If (.Not. concon) Then
    Write (999, *) ' STOP: ERROR IN RATELI'
    Stop ' ERROR IN RATELI'
  End If

! -----sobolev case

  If (.Not. optcmf) Then

    Do ll = 1, nd

      Read (17, Rec=ll)(blevel(i), i=1, nrec)
      Call factinc(r(ll), ll)
      velr = v(ll)/r(ll)
      xmust = sqrt(1.D0-1.D0/r(ll)/r(ll))

      Do nato = 1, nat
        ilo = ilow(ll, nato)
        ima = imax(ll, nato)
!       print*,nato,'  ',ll,'  ',ilo,'  ',ima
        If (ima<ilo) Cycle
        numion = igenio(nato, ilo)
        nfir = ifirsl(numion)
        numion = igenio(nato, ima) + 1
        If (numion>iong) Then
          Write (999, *) ' STOP: ERROR IN RATELI - NUMION'
          Stop ' ERROR IN RATELI - NUMION'
        End If
        nlas = ifirsl(numion)
        Call sobprep(ll, velr, gradv(ll), nfir, nlas, blevel, xne(ll), sr, &
          vmax, nato, xmust, clf(ll), fic(ll), tcl_fac_line(ll))
      End Do

    End Do

    Return

!   -----cmf-treatment for some lines

  Else

!   reset of indexsa

    indexsa(1:id_nttrd) = 0
!   JO Dec 2016: initialize xlamcmf
    If (startcmf) xlamcmf = 0.

depthloop: Do ll = 1, nd

      Read (17, Rec=ll)(blevel(i), i=1, nrec)

!     ---  if all lines in cmf, the next statement can be skipped

      Call factinc(r(ll), ll)
      velr = v(ll)/r(ll)
      xmust = sqrt(1.D0-1.D0/r(ll)/r(ll))

natoloop: Do nato = 1, nat
        ilo = ilow(ll, nato)
        ima = imax(ll, nato)
!       print*,nato,'  ',ll,'  ',ilo,'  ',ima

        If (ima<ilo) Then
          If (nato==1 .And. ll==1) Then
            Print *, nato, '  ', ll, '  ', ilo, '  ', ima
            Print *, ' NO LINES FOR NATO =1 AND LL=1'
            Print *, ' IN THIS CASE, CMFPREP IS NOT CALLED AND INDEXCMF'
            Print *, ' MIGHT BE NOT CORRECTLY UPDATED'
            Write (999, *) ' STOP: NO LINES FOR NATO =1 AND LL=1'
            Stop ' NO LINES FOR NATO =1 AND LL=1'
          End If
          Cycle natoloop
        End If

        numion = igenio(nato, ilo)
        nfir = ifirsl(numion)
        numion = igenio(nato, ima) + 1
        If (numion>iong) Then
          Write (999, *) ' STOP: ERROR IN RATELI - NUMION'
          Stop ' ERROR IN RATELI - NUMION'
        End If
        nlas = ifirsl(numion)

!       ---     calculation of min/max levels for cmf-transfer: Here we
!       prepare the decision whether SA even in general CMF situation, &
!       in cases where we have not everywhere occupation numbers

        If (startcmf) Then
          nfirm = 0
          nlasm = 10000
          Do lx = 1, nd
            ilo = ilow(lx, nato)
            ima = imax(lx, nato)
            numion = igenio(nato, ilo)
            nfir1 = ifirsl(numion)
            numion = igenio(nato, ima) + 1
            nlas1 = ifirsl(numion)
            nfirm = max(nfirm, nfir1)
            nlasm = min(nlasm, nlas1)
          End Do

          numion = igenio(nato, imianl(nato))
          nfi = ifirsl(numion)
          numion = igenio(nato, imaanl(nato)) + 1

          If (numion>iong) Then
            Write (999, *) ' STOP: ERROR IN RATELI - NUMION'
            Stop ' ERROR IN RATELI - NUMION'
          End If
          nla = ifirsl(numion)

        End If

!       note: nfirm, nlasm only needed for startcmf = .true.,
!       i.e., at first depth point of first cmf iteration
!       the same is true for nfi,nla

!       print*,nato,'  ',ll,'  ',ilo,'  ',ima
!       print*,nato,'  ',ll,'  ',nfi,'  ',nla
!       print*,nato,'  ',ll,'  ',nfir,'  ',nlas
!       print*,nato,'  ',ll,'  ',nfirm,'  ',nlasm
        Call cmfprep(ll, velr, gradv(ll), nfir, nlas, blevel, xne(ll), sr, &
          vmax, nato, xmust, clf(ll), nfirm, nlasm, nfi, nla, &
          startcmf, optmixed, ipres, r, v, fic(ll), tcl_fac_line(ll), &
          tcl_fac_cont(ll))

      End Do natoloop
      If (startcmf) Then

!       ***  index array so that cmf-line wavelengths are in order

!       JO Dec. 2016: calculate index-array only once
        icount = id_nttrd
        Call indexx(icount, xlamcmf, indxlamc)
        Print *
jloop:  Do j = 1, id_nttrd
          ii = indxlamc(j)
          xlam = xlamcmf(ii)
          icmf = indexcmf(ii)
          If (xlam==0.D0 .And. icmf/=1) Then
            Write (999, *) 'STOP: ERROR IN XLAM=0(1)!'
            Stop ' ERROR IN XLAM=0(1)!'
          End If
          If (icmf==1) Then

            If (xlam/=0.D0) Print *, ' LINE AT ', xlam, &
              ' TREATED IN SOBOLEV APPROX.'

          Else If (icmf==2) Then
            Print *, ' LINE AT ', xlam, ' TREATED IN CMF '

          Else
            Write (999, *) ' STOP: ERROR IN INDEXCMF(RATELI)'
            Stop ' ERROR IN INDEXCMF(RATELI)'
          End If
        End Do jloop

        Print *
        If (ipres==1) Then

!         ***     file indexcmf exists, no further action required

          Continue

        Else If (ipres==2) Then

          If (optmixed==1.) Then
            Write (999, *) ' STOP: OPTMIXED NOT IMPLEMENTED YET'
            Stop ' OPTMIXED NOT IMPLEMENTED YET'
            Open (1, Err=100, File=trim(modnam)//'/INDEXCMF', Status='NEW')
            Write (1, Fmt=*)(indexcmf(ii), ii=1, id_nttrd)
            Close (1)
          End If

        Else
          Write (999, *) ' STOP: ERROR IN IPRES (RATELI)'
          Stop ' ERROR IN IPRES (RATELI)'
        End If
      End If

!     ***  startcmf only for first depth point

      startcmf = .False.

    End Do depthloop

!   perform cmf-calc.

!   so far, no necessity to update subroutine CMF, since all quantities
!   already corrected
!   this statement MUST be executed also after restart of CMF_ALL
    Call cmf(r, v, gradv, temp, vmax, teff, corrfc, nd)

    Print *

    Return

  End If

100 Continue

  Write (999, *) ' STOP: ERROR IN OPENING FILE INDEXCMF'
  Stop ' ERROR IN OPENING FILE INDEXCMF'

End Subroutine

!-----------------------------------------------------------------------

Subroutine cmf(r, v, gradv, temp, vmax, teff, corrfc, nd)

! performs cmf calculations, both for the single line as well as the
! overlap case

! note that all overlap calc. are done with respect to the "def"
! delta nu = nu_o * v/c and
! x = (nu/nu_o - 1) c/v

  Use :: nlte_type
  Use :: fund_const, Only: clight
  Use :: nlte_dim

  Use :: nlte_var, Only: indxlamc, indexcmf, indexsa, xlamcmf, vdopcmf, &
    xmaxdop, tlumat, tulmat, blucmf, bulcmf, aulcmf, lablinelo, lablineup

  Use :: nlte_opt, Only: optsing
  Use :: cmf_multi, Only: indov, indoverlap_order, infotab, usave, vsave, &
    iplus, iminus


  Implicit None

  Integer (i4b), Parameter :: nd1 = id_ndept, np = id_npoin
  Integer (i4b), Parameter :: nf = id_nfcmf


  Integer (i4b) :: nd
  Real (dp) :: vmax, teff, corrfc
  Real (dp), Dimension (nd1) :: r, v, gradv, temp

  Real (dp), Dimension (nd1) :: ajlmean, alo


  Integer (i4b), Dimension (:), Allocatable :: index, indexov
  Real (dp), Dimension (:), Allocatable :: lamblue
  Real (dp), Dimension (:, :), Allocatable :: tluaux, tulaux

  Integer (i4b) :: ii, ii1, j, jj, jj1, jk, l, jstart, jend, nov, nov1, iov

  Real (dp) :: xlam, vdop, snew, xov, xnf, xmax, xlam1, dl, dl1

  Logical :: flag, inversion
! prepare overlap calc.

  If (optsing) Then
!   forces single line treatment
    indov = 0
    Call calcxmax(vmax, teff, temp, nd)
  Else
    Call prepov(r, v, gradv, temp, vmax, teff, nd)
  End If

! do the cmf transport for all lines, either in single or overlap approach

  j = 1
! OPEN(333,FILE='OX') !obsolete
cmfloop: Do While (j<=id_nttrd)

    ii = indxlamc(j)
    xlam = xlamcmf(ii)

    If (indexcmf(ii)==2) Then

      vdop = vdopcmf(ii)
      If (vdop==0.D0) Then
        Write (999, *) ' STOP: ERROR IN VTHCMF'
        Stop ' ERROR IN VTHCMF'
      End If

      If (indov(j)==0) Then
!       no overlap

        Call cmfsing(ii, r, v, gradv, temp, xlam, vdop, vmax, teff, corrfc, &
          ajlmean, alo, inversion)
        If (j==1 .Or. inversion) Then
          If (j==1) Write (999, *) ' ONLY FIRST EXPLICIT LINE HERE &
            &(FASTWIND.LOG), FOR REST SEE OUTPUT'
          Write (999, Fmt=100) xlam, lablinelo(ii), lablineup(ii)
        End If
        Write (*, Fmt=100) xlam, lablinelo(ii), lablineup(ii)

        Do l = 1, nd
          snew = ajlmean(l) - alo(l)*tulmat(ii, l)
!         if(lablinelo(ii).eq.'HE11S1  '.and.lablineup(ii).eq.'HE12P1  ') &
!         print*,'test_HE11S1 ',l,ajlmean(l),alo(l),tulmat(ii,l),snew

          If (snew<0.D0 .Or. abs(alo(l))>1.D0) Then
            If (.Not. inversion) Then
              Print *, ' WARNING!!!! SOMETHING WRONG WITH LINE-ALO'
              Write (999, *) ' WARNING!!!! SOMETHING WRONG WITH LINE-ALO'
            End If
            If (inversion) Then
              Print *, ' WARNING!!!! Jbar corrected at', l
              Write (999, *) ' WARNING!!!! Jbar corrected at', l
!             JO at some point, this needs to be improved (also in
!             cmf_all.f90)
              ajlmean(l) = 0.
            End If
            alo(l) = 0.D0
            snew = ajlmean(l)
          End If

          tlumat(ii, l) = blucmf(ii)*snew
          tulmat(ii, l) = bulcmf(ii)*snew + aulcmf(ii)*(1.D0-alo(l))
        End Do


      Else

!       line overlap (intrinsic + wind induced)


        If (indov(j)/=1) Then
          Write (999, *) ' STOP: ERROR IN INDOV'
          Stop ' ERROR IN INDOV'
        End If
        jstart = j
        j = j + 1
        nov = 1
        Do While (indov(j)/=3)
          j = j + 1
          If (indexcmf(indxlamc(j))==2) nov = nov + 1
        End Do
        jend = j

        nov = nov + 1
        nov1 = jend - jstart + 1

        If (nov/=nov1) Then
          Print *, jstart, jend, nov, nov1
          Write (999, *) &
            ' STOP: SA LINE IN OVERLAP COMPLEX, SHOULD NEVER HAPPEN'
          Stop ' SA LINE IN OVERLAP COMPLEX, SHOULD NEVER HAPPEN'
        End If

!       prepare info concering blue-wing boundary condition
!       ------------------------------------------------------------------------

!       locals
        If (allocated(lamblue)) Deallocate (lamblue, index, indexov, tluaux, &
          tulaux)
!       common
        If (allocated(indoverlap_order)) Deallocate (indoverlap_order, &
          infotab, usave, vsave, iplus, iminus)

        Allocate (lamblue(nov), index(nov), indoverlap_order(nov), &
          infotab(nov), usave(nd1,np-1,nov), vsave(nd1-1,np-1,nov), &
          iplus(np-nd1+1,nov), iminus(np-1,nov), indexov(nov), tluaux(nd,nov), &
          tulaux(nd,nov))

        infotab%index = 0
        infotab%indbw = 0
        infotab%xov = 0.
        infotab%if = 0.

        Do jj = jstart, jend
          ii = indxlamc(jj)
          vdop = vdopcmf(ii)
          If (vdop==0.D0) Then
            Write (999, *) ' STOP: ERROR IN VTHCMF'
            Stop ' ERROR IN VTHCMF'
          End If
!         blue-wing freq., from def. of dnu
          lamblue(jj-jstart+1) = xlamcmf(ii)/(1.D0+xmaxdop(ii)*vdop/clight)
        End Do

!       sort lines for decreasing blue-wing frequencies
        Call indexx(nov, lamblue, index)

        Do jj = 1, nov
          jj1 = index(jj) + jstart - 1
!         JJ1 is now the index corresponding to blue-wing order in INDXLAMC
          ii = indxlamc(jj1)
!         array INDOVERLAP_ORDER does not require any rearrangement,
!         but can be used directly to find indices corresponding to
!         blue wing order in overlap complex
          indoverlap_order(jj) = ii
        End Do

!       Find out which line (II1) gets the b.w. condition from what line (II)
!       and calculate frequency (XOV) where intensities from line (II) have
!       to be saved.
!       INFOTAB%IF then gives the corresponding freq. nr.

!       If NO intrinsic line overlap, find out which line (II1) gets the bw.
!       cond. from the red wing of line (II). In this case, XOV is larger
!       than XMAX(II), and
!       INFOTAB%IF is set to NF


        Do jj = 1, nov - 1
          ii = indoverlap_order(jj)
          If (jj==1) infotab(jj)%index = ii
search:   Do jj1 = jj + 1, nov
            ii1 = indoverlap_order(jj1)
!           defined with opposite sign, from def. of x
            xov = (1.D0-xlamcmf(ii)/lamblue(index(jj1)))*clight/vdopcmf(ii)
!           since also xmax is positive
            xmax = xmaxdop(ii)
            If (xov>xmax) Then
!             for narrow line, xov can be larger than xmax, altough iov was
!             possible for prior line
              If (infotab(jj1)%if==0. .Or. infotab(jj1)%if==dble(nf)) Then
!               no intrinsic line-overlap, use closest line (this is a MUST),
!               since otherwise the range in between two lines would not be
!               line-free!

                infotab(jj1)%index = ii1
                infotab(jj1)%indbw = ii
                infotab(jj1)%xov = -xov
                infotab(jj1)%if = dble(nf)
              End If
              Exit search
            End If
            If (infotab(jj1)%if==0. .Or. infotab(jj1)%if==dble(nf)) Then
!             intrinsic line-overlap, use line with most blueward freq.
!             (numerically most exact)
              infotab(jj1)%index = ii1
              infotab(jj1)%indbw = ii
              infotab(jj1)%xov = -xov
              xnf = 1.D0 + (xmax+xov)*(nf-1)*.5D0/xmax ! (-(-XOV)
              If (xnf<=1. .Or. xnf>=dble(nf)) Then
                Write (999, *) ' STOP: ERROR IN XNF'
                Stop ' ERROR IN XNF'
              End If
              infotab(jj1)%if = xnf
            End If
          End Do search
        End Do

!       perform some tests
        If (infotab(1)%indbw/=0) Then
          Write (999, *) ' STOP: ERROR 1 IN INFOTABLE'
          Stop ' ERROR 1 IN INFOTABLE'
        End If

        Do jj = 2, nov
          If (infotab(jj)%if==dble(nf)) Then
            If (infotab(jj)%indbw/=infotab(jj-1)%index) Then
              Write (999, *) ' STOP: ERROR 2 IN INFOTABLE'
              Stop ' ERROR 2 IN INFOTABLE'
            End If
          End If
!         for tests only
!         PRINT*,JJ,INFOTAB(JJ)%INDEX,INFOTAB(JJ)%INDBW, &
!         INFOTAB(JJ)%XOV,INFOTAB(JJ)%IF
        End Do

!       ----------------------------------------------------------------------



!       perform cmf calc. for overlap interval, in order of descending
!       bw. frequencies

        Do jj = 1, nov
          jj1 = indxlamc(jj+jstart-1)
          ii = indoverlap_order(jj)
          xlam = xlamcmf(ii)
          vdop = vdopcmf(ii)
          dl = xlam*xmaxdop(ii)*vdop/clight

!         ----------------------------------------------------------------------

!         Find out intrinsically overlapping lines to allow for calc. of
!         total background opacities/em.
          indexov = 0
          iov = 0
          Do jk = 1, jj - 1
            ii1 = indoverlap_order(jk)
            xlam1 = xlamcmf(ii1)
            dl1 = xlam1*xmaxdop(ii1)*vdopcmf(ii1)/clight
            If (abs(xlam1-xlam)<dl+dl1) Then
              iov = iov + 1
              indexov(iov) = ii1
            End If
          End Do
          Do jk = jj + 1, nov
            ii1 = indoverlap_order(jk)
            xlam1 = xlamcmf(ii1)
            dl1 = xlam1*xmaxdop(ii1)*vdopcmf(ii1)/clight
            If (abs(xlam1-xlam)<dl+dl1) Then
              iov = iov + 1
              indexov(iov) = ii1
            End If
          End Do
!         test
          If (iov==0) Then
            If (jj/=1. .And. infotab(jj)%if/=dble(nf)) Then
              Write (999, *) ' STOP: ERROR IN INDEXOV'
              Stop ' ERROR IN INDEXOV'
            End If
!           PRINT*,'LINE AT ',XLAMCMF(II),':NO !INTRINSIC OVERLAP'
          Else
            If (jj/=1 .And. infotab(jj)%if/=dble(nf)) Then
              flag = .False.
              Do jk = 1, iov
                If (indexov(jk)==infotab(jj)%indbw) flag = .True.
              End Do
              If (.Not. flag) Then
                Print *, infotab(jj)%index, infotab(jj)%indbw, &
                  infotab(jj)%xov, infotab(jj)%if
                Print *, indexov(1:iov)
                Write (999, *) ' STOP: BW INDEX NOT FOUND IN INDEXOV'
                Stop ' BW INDEX NOT FOUND IN INDEXOV'
              End If
            End If
!           PRINT*,'LINE AT ',XLAMCMF(II),':!INTRINSIC OVERLAP FOR ', &
!           INDEXOV(1:IOV)
          End If

!         ----------------------------------------------------------------------

          Call cmfmulti(jj, ii, r, v, gradv, temp, xlam, vdop, vmax, teff, &
            corrfc, ajlmean, alo, nov, indexov, iov)
          If (ii==jj1) Then
            Write (999, Fmt=110) xlam, lablinelo(ii), lablineup(ii), iov
            Write (*, Fmt=110) xlam, lablinelo(ii), lablineup(ii), iov
          Else
            Write (999, Fmt=120) xlam, lablinelo(ii), lablineup(ii), iov
            Write (*, Fmt=120) xlam, lablinelo(ii), lablineup(ii), iov
          End If

          Do l = 1, nd
            snew = ajlmean(l) - alo(l)*tulmat(ii, l)
            If (snew<0.D0 .Or. abs(alo(l))>1.D0) Then
              Print *, ' WARNING!!!! SOMETHING WRONG WITH LINE-ALO'
              Write (999, *) ' WARNING!!!! SOMETHING WRONG WITH LINE-ALO'
              alo(l) = 0.D0
              snew = ajlmean(l)
            End If

!           avoid to overwrite TLUMAT (->OPALM) and TULMAT(->SLINEM)

            tluaux(l, jj) = blucmf(ii)*snew
            tulaux(l, jj) = bulcmf(ii)*snew + aulcmf(ii)*(1.D0-alo(l))
          End Do
        End Do
!       overlap complex finished, can overwrite TLUMAT and TULMAT

        Do jj = 1, nov
          ii = indoverlap_order(jj)
          Do l = 1, nd
            tlumat(ii, l) = tluaux(l, jj)
            tulmat(ii, l) = tulaux(l, jj)
          End Do
        End Do

        j = jend

      End If !                      SINGLE OR OVERLAP

    Else If (indexcmf(ii)==1) Then

      If (xlam==0.D0) Then
        If (indexsa(ii)/=0) Then
          Write (999, *) ' STOP: ERROR IN XLAM=0(2)!'
          Stop ' ERROR IN XLAM=0(2)!'
        End If
      Else
        Write (999, Fmt=130) xlam, lablinelo(ii), lablineup(ii), indexsa(ii)
        Write (*, Fmt=130) xlam, lablinelo(ii), lablineup(ii), indexsa(ii)
      End If
    End If

    j = j + 1

  End Do cmfloop
! CLOSE(333)  !obsolete

  Return

100 Format ('CMFT AT', F12.3, ' A ', A6, ' TO ', A6, ' SINGLE LINE')
110 Format ('CMFT AT', F12.3, ' A ', A6, ' TO ', A6, ' OVERLAP, IOV =', I2)
120 Format ('CMFT AT', F12.3, ' A ', A6, ' TO ', A6, &
    ' OVERLAP, ORDER CHANGED, IOV = ', I2)
130 Format ('SA-T AT', F12.3, ' A ', A6, ' TO ', A6, ' FOR ', I3, ' DEPTH', &
    ' POINTS')

End Subroutine

!-----------------------------------------------------------------------

Subroutine prepov(r, v, gradv, temp, vmax, teff, nd)

! prepares consideration of line overlap (intrinsic and wind induced)

  Use :: nlte_type
  Use :: fund_const, Only: clight
  Use :: nlte_dim
  Use :: fastwind_params, Only: tauscrit

  Use :: nlte_var, Only: indxlamc, xlamcmf, vdopcmf, xmaxdop, indexcmf, &
    tlumat, z

  Use :: nlte_opt, Only: iov_only
  Use :: cmf_multi, Only: indov, lmaxjp, fmax, fmin, uv, zray


  Implicit None

  Integer (i4b), Parameter :: nd1 = id_ndept, np = id_npoin

  Integer (i4b) :: nd
  Real (dp) :: vmax, teff
  Real (dp), Dimension (nd1) :: r, v, gradv, temp

  Integer (i4b) :: ii, j, jj, jlast, jstart, jold, ll, ja, ia
  Integer (i4b) :: jp, lmax, l

  Real (dp) :: delta, ddop, xlstart, xlam, xlamold, taus, nextlam, nextlamold
  Real (dp) :: fm, zl

  Logical :: first, core
  Logical, Save :: do_once = .True.

  nextlamold = 0.
  indov = 0
  first = .True.
  jlast = 0
  Do j = id_nttrd, 1, -1
    ii = indxlamc(j)
    If (indexcmf(ii)==2) Then
      jlast = j
      Exit
    End If
  End Do
  If (jlast==0) Then
    Write (999, *) ' STOP: JLAST NOT FOUND!'
    Stop ' JLAST NOT FOUND!'
  End If

! calc. xmax
  Call calcxmax(vmax, teff, temp, nd)

! NOTE: if SA line inside overlap complex, end complex

! INDOV = 0 SINGLE LINE OR SOBO
! INDOV = 1 BEGIN OF OVERLAP COMPLEX
! INDOV = 2 INSIDE OF OVERLAP COMPLEX
! INDOV = 3 END OF OVERLAP COMPLEX

! note: INDOV ordered in actual sequence of restwavelengths

  Do j = 1, id_nttrd

    ii = indxlamc(j)

    xlam = xlamcmf(ii)
    If (indexcmf(ii)==2) Then
      Do ll = 1, nd
        taus = tlumat(ii, ll)/gradv(ll)
!       use velocities only where significant contribution
        If (taus>tauscrit) Exit
      End Do
      If (ll>nd) ll = nd

!     maximum influence due to expansion and thermal line width
      ddop = xlam*xmaxdop(ii)*vdopcmf(ii)/clight
      delta = xlam*(2.*v(ll)*vmax)/clight + ddop
!     note that previous line can have wider range of influence
      nextlam = max(xlam+delta, nextlamold)
      If (first) Then
        xlstart = xlam
        jstart = j
        nextlamold = nextlam
        xlamold = xlam
        jold = j
        first = .False.
        Cycle
      End If

      If (j==jlast) Then

!       for lines inside overlap complex, check whether SA-lines present;
!       if so, break up overlap complex into subcomplexes in order
!       to be able have SA-lines as single (correction of previous and
!       actual lines required)

        If (xlstart/=xlamold) Then ! OVERLAP
!         first line
          indov(jstart) = 1
          Do jj = jstart + 1, j - 1 ! all lines inside complex
            ia = indxlamc(jj)
            If (indexcmf(ia)==1) Then ! SA-LINE
              indov(jj) = 0
              ja = indov(jj-1)
              If (ja==1) indov(jj-1) = 0 ! correction of previous index
              If (ja==2) indov(jj-1) = 3 ! similar
            Else
              indov(jj) = 2 !       standard, if no SA lines in between
              ja = indov(jj-1)
              If (ja==0) indov(jj) = 1 ! correction of actual index
            End If
          End Do
!         last line J (end of complex)
          ia = indxlamc(j)
          If (indexcmf(ia)==1) Then ! SA-LINE
            indov(j) = 0
            ja = indov(j-1)
            If (ja==1) indov(j-1) = 0 ! correction of previous index
            If (ja==2) indov(j-1) = 3 ! similar
          Else
            indov(j) = 3 !          standard, if not SA line
            ja = indov(j-1)
            If (ja==0) indov(j) = 0 ! correction of actual index
          End If

        End If

!       corrected for thermal width
      Else If (xlam-ddop>nextlamold) Then
!       same philosophy as above
        If (xlstart/=xlamold) Then ! OVERLAP
          indov(jstart) = 1
          Do jj = jstart + 1, jold - 1
            ia = indxlamc(jj)
            If (indexcmf(ia)==1) Then ! SA-LINE
              indov(jj) = 0
              ja = indov(jj-1)
              If (ja==1) indov(jj-1) = 0 ! correction of previous index
              If (ja==2) indov(jj-1) = 3 ! similar
            Else
              indov(jj) = 2 !       standard
              ja = indov(jj-1)
              If (ja==0) indov(jj) = 1 ! correction of actual index
            End If
          End Do
          ia = indxlamc(jold)
          If (indexcmf(ia)==1) Then ! SA-LINE
            indov(jold) = 0
            ja = indov(jold-1)
            If (ja==1) indov(jold-1) = 0 ! correction of previous index
            If (ja==2) indov(jold-1) = 3 ! similar
          Else
            indov(jold) = 3 !       standard
            ja = indov(jold-1)
            If (ja==0) indov(jold) = 0 ! correction of actual index
          End If

        End If
        xlstart = xlam
        jstart = j
      End If

      nextlamold = nextlam
      xlamold = xlam
      jold = j
    End If
  End Do

! for tests
! DO J = 1,ID_NTTRD
! II = INDXLAMC(J)
! XLAM = XLAMCMF(II)
! PRINT*,J,XLAM,INDEXCMF(II),INDOV(J)
! ENDDO

  If (iov_only .Or. .Not. do_once) Return

  do_once = .False.

  Return

! prepare some variables required for wind induced overlap, old version
! NOT REQUIRED for new version
  Write (999, *) ' STOP: SHOULD NEVER STOP HERE (PREPOV)'
  Stop ' SHOULD NEVER STOP HERE (PREPOV)'

  Do jp = 1, np - 1
    lmax = min0(nd, np+1-jp)
    lmaxjp(jp) = lmax
    core = lmax == nd
    fm = z(1, jp)/r(1)*v(1)
    fmax(jp) = fm

    If (.Not. core) Then
      fmin(jp) = -fm
    Else
      fmin(jp) = z(nd, jp)*v(nd)
    End If

    uv(1, jp) = fm
    zray(1, jp) = z(1, jp)
    Do l = 2, lmax
      zl = z(l, jp)
      uv(l, jp) = zl/r(l)*v(l)
      zray(l, jp) = zl
    End Do

    If (.Not. core) Then
      Do l = 1, lmax - 1
        ll = 2*lmax - l
        uv(ll, jp) = -uv(l, jp)
        zray(ll, jp) = -zray(l, jp)
      End Do
    End If

  End Do

  Return

End Subroutine


!-----------------------------------------------------------------------

Subroutine calcxmax(vmax, teff, temp, nd)

! calculates xmax as function of beta

  Use :: nlte_type
  Use :: fund_const, Only: wpi => sqpi
  Use :: nlte_dim

  Use :: nlte_var, Only: indxlamc, indexcmf, vdopcmf, xmaxdop, &
    opalm => tlumat, opacm

  Implicit None
  Integer (i4b), Parameter :: nd1 = id_ndept

  Integer (i4b) :: nd
  Real (dp) :: teff, vmax
  Real (dp), Dimension (nd1) :: temp

  Integer (i4b) :: ii, j, ltemp
  Real (dp) :: betap, beta1, xmax0

! index for calculating beta_puls

  Do j = 1, nd
    If (temp(j)>teff) Go To 100
  End Do
  Write (999, *) ' STOP: ERROR IN TEMP.STRAT.'
  Stop ' ERROR IN TEMP.STRAT.'

100 Continue

  ltemp = j - 1
  Print *

  Do j = 1, id_nttrd

    ii = indxlamc(j)
    If (indexcmf(ii)==2) Then

      betap = opacm(ii, ltemp)*vdopcmf(ii)/(opalm(ii,ltemp)*vmax)

      If (betap<-1.D0) Then
        Print *, ' WARNING: BETA_P NEGATIVE IN CALCXMAX'
!       occurs too often, thus only in output
!       WRITE(999,*) ' WARNING: BETA_P NEGATIVE IN CALCXMAX'
      End If

      xmax0 = 3.D0
      beta1 = abs(betap)

      If (beta1<.01D0) xmax0 = max(3.D0, sqrt(-log(.01D0*beta1*wpi)))
      xmaxdop(ii) = xmax0
    End If

  End Do

  Return

End Subroutine

!-----------------------------------------------------------------------

Subroutine cmfsing(ii, r, v, dvdr, t, xlambda, vdop, vmax, teff, corrfc, &
  ajlmean, alo, inversion)

! prepares and carries out cmf transfer for single lines

  Use :: nlte_type
  Use :: nlte_dim
  Use :: run_once, Only: start_cmfsing
  Use :: nlte_var, Only: tauz1, tau1, pp1, ak1, az1, slinek, z, p, &
    opalm => tlumat, slinem => tulmat, opacm, scontm, w0 => vp, w0last, &
    xmaxdop

  Implicit None

! .. parameters ..
  Integer (i4b), Parameter :: nd = id_ndept, np = id_npoin
  Integer (i4b), Parameter :: nf = id_nfcmf
! ..
! .. scalar arguments ..
  Real (dp) :: corrfc, teff, vdop, vmax, xlambda
  Integer (i4b) :: ii
! ..
! .. array arguments ..
  Real (dp) :: ajlmean(nd), alo(nd), dvdr(nd), r(nd), t(nd), v(nd)
! ..
! .. local scalars ..
  Real (dp) :: aic, aux, dbdr, deltax, etak, opakk, opalk, rl, xmue, xmue2
  Integer (i4b) :: jp, k, l, lmax, lz
  Logical :: inversion
! ..
! .. local arrays ..
  Real (dp) :: etac(nd), opac(nd), opal(nd), phi(nf, nd), pp(nd), &
    pweight(nf, nd), scont(nd), sline(nd), ubluwi(nd, np-1), ucmf(nd), &
    vcmf(nd-1)
! ..
! .. external subroutines ..
  External :: cont2, diffus, fgrid, ray
! ..
! .. intrinsic functions ..
! INTRINSIC MIN0
! ..
! .. data statements ..

! ..

  If (start_cmfsing) Then

!   ***  might have already been calculated in cmf_simple (nlte_approx)
!   ***  anyway, recalulate, since previous update

!   ***  line-independent quantities

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
        pp1(l, jp) = (xmue2*dvdr(l)+(1.D0-xmue2)*v(l)/rl)
      End Do

    End Do jploop1

    start_cmfsing = .False.
  End If

  Do l = 1, nd
    opal(l) = opalm(ii, l)
    sline(l) = slinem(ii, l)
    opac(l) = opacm(ii, l)
    scont(l) = scontm(ii, l)
    etac(l) = scont(l)*opac(l)
  End Do

! blue wing boundary

  Call diffus(xlambda, t, r, nd, aic, dbdr)

  Call cont2(nd, np, r, p, z, scont, opac, aic, dbdr, corrfc, ubluwi)

  Call fgrid(nf, nd, phi, pweight, deltax, vdop, vmax, teff, t, xmaxdop(ii))

  inversion = .False.


! check whether inversion


invcheck: Do l = 1, nd
    opakk = phi((nf+1)/2, l)*opal(l) + opac(l)
!   if(xlambda.gt.4686..and.xlambda.lt.4689.) then
!   print*,l,opakk,PHI((NF+1)/2,L),opal(l),opac(l)
!   endif
    If (opakk<0.D0 .And. .Not. inversion) Then
      Write (999, *) ' INVERTED LEVELS FOR FOLLOWING TRANSITION '
      Print *, ' INVERTED LEVELS FOR FOLLOWING TRANSITION '
      inversion = .True.
      Exit invcheck
    End If
  End Do invcheck


  If (.Not. inversion) Then
!   STANDARD TREATMENT
    Do l = 1, nd
      Do k = 1, nf
        opalk = phi(k, l)*opal(l)
        etak = sline(l)*opalk + etac(l)
        opakk = opalk + opac(l)
        slinek(k, l) = etak/opakk
!       if(xlambda.gt.4686..and.xlambda.lt.4689..and.k.eq.(nf+1)/2) then
!       print*,l,opakk,opalk,opac(l)
!       print*,l,etak,sline(l)*opalk,etac(l)
!       print*,l,slinek(k,l)
!       print*
!       endif
        ak1(k, l) = 1.D0/opakk
        If (l/=1) az1(k, l-1) = 2.D0/(opakk+1.D0/ak1(k,l-1))
      End Do
    End Do

  Else
!   INVERSION
    Do l = 1, nd
      Do k = 1, nf
        opalk = phi(k, l)*opal(l)
        etak = sline(l)*opalk + etac(l)
        opakk = opalk + opac(l)

!       ATTENTION: for inversion, ak1 is now the actual opacity,
!       and slinek is just the emisivity

        ak1(k, l) = opakk
        slinek(k, l) = etak
        If (l/=1) az1(k, l-1) = (opakk+ak1(k,l-1))/2.D0
      End Do
    End Do

  End If


  ajlmean = 0.D0
  alo = 0.D0

jploop2: Do jp = 1, np - 1
    lmax = min0(np+1-jp, nd)
    lz = lmax - 1

!   bluewing boundary condition

    ucmf(1) = ubluwi(1, jp)

    Do l = 1, lz
      ucmf(l+1) = ubluwi(l+1, jp)
      aux = .5D0*(opac(l)+opac(l+1))
      vcmf(l) = (ucmf(l+1)-ucmf(l))/aux/(z(l,jp)-z(l+1,jp))
    End Do

    Do l = 1, lmax
      pp(l) = pp1(l, jp)/deltax
    End Do

    If (.Not. inversion) Then
!     STANDARD SOLUTION OVER TAU
      Call ray(.True., 1, 1, jp, z(:,jp), r, nd, np, ucmf, vcmf, sline, phi, &
        pweight, deltax, opac, etac, opal, aic, dbdr, corrfc, w0(:,jp), &
        w0last, ajlmean, alo, tauz1(:,jp), tau1(:,jp), pp, xlambda)

    Else
!     NEW FORMULATION OVER Z
      Call ray_inv(.True., 1, 1, jp, z(:,jp), r, nd, np, ucmf, vcmf, &
        pweight, deltax, opac, etac, opal, aic, dbdr, corrfc, w0(:,jp), &
        ajlmean, tauz1(:,jp), tau1(:,jp), pp, xlambda)
    End If

  End Do jploop2

! NOTE: SINCE ALO IS NOT MODIFIED IN THE INVERTED CASE, WE SIMPLY DON'T TOUCH
! IT

! DO L=1,ND
! print*,l,' ',sline(l),' ',ajlmean(l),'  ',alo(l),'  ',opal(l)
! print*,l,' ',ajlmean(l),'  ',alo(l)
! END DO

  Return

End Subroutine

!-----------------------------------------------------------------------

Subroutine cmfmulti(jj, ii, r, v, dvdr, t, xlambda, vdop, vmax, teff, corrfc, &
  ajlmean, alo, nov, indexov, iov)

! prepares and carries out cmf transfer for given line in overlap complex

  Use :: nlte_type
  Use :: fund_const, Only: clight, wpi => sqpi
  Use :: nlte_dim
  Use :: run_once, Only: start_cmfmulti
  Use :: nlte_var, Only: tauz1, tau1, pp1, ak1, az1, slinek, z, p, &
    opalm => tlumat, slinem => tulmat, opacm, scontm, w0 => vp, w0last, vturb, &
    vdopcmf, xlamcmf, xmaxdop

  Use :: nlte_opt, Only: iov_only
  Use :: cmf_multi, Only: infotab, usave, vsave

  Implicit None

! .. parameters ..
  Integer (i4b), Parameter :: nd = id_ndept, np = id_npoin
  Integer (i4b), Parameter :: nf = id_nfcmf
! ..
! .. scalar arguments ..
  Real (dp) :: corrfc, teff, vdop, vmax, xlambda
  Integer (i4b) :: jj, ii, nov, iov
! ..
! .. array arguments ..
  Integer (i4b) :: indexov(iov)
  Real (dp) :: ajlmean(nd), alo(nd), dvdr(nd), r(nd), t(nd), v(nd)
! ..
! .. local scalars ..
  Real (dp) :: aic, aux, dbdr, deltax, etak, opakk, opalk, rl, xmue, xmue2, &
    opabl, sbl, etabl, xmax, xmax1, vdopl, del, del1, cv, phikl, xlami, w, w1, &
    xk, xki, lred, lblu, xmaxr, xmaxb, xov, dew, delta, deltax1

  Integer (i4b) :: i, jp, k, l, lmax, lz, ii1, iir, iib
  Logical :: inversion, flag
! ..
! .. local arrays ..
  Real (dp) :: etac(nd), opac(nd), opal(nd), phi(nf, nd), pp(nd), &
    pweight(nf, nd), scont(nd), sline(nd), ubluwi(nd, np-1), &
    vbluwi(nd-1, np-1), ucmf(nd), vcmf(nd-1), opabg(nf, nd), etabg(nf, nd)

! ..
! .. external subroutines ..
  External :: cont2, diffus, fgrid, ray
! ..
! .. intrinsic functions ..
! INTRINSIC MIN0
! ..
! .. data statements ..

  If (infotab(jj)%index/=ii) Then
    Write (999, *) ' STOP: ERROR IN JJ <=> II'
    Stop ' ERROR IN JJ <=> II'
  End If

  If (start_cmfmulti) Then

!   ***  line-independent quantities

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
        pp1(l, jp) = (xmue2*dvdr(l)+(1.D0-xmue2)*v(l)/rl)
      End Do

    End Do jploop1

    start_cmfmulti = .False.
  End If

  Do l = 1, nd
    opal(l) = opalm(ii, l)
    sline(l) = slinem(ii, l)
  End Do

! continuum background
  Do l = 1, nd
    opabl = opacm(ii, l)
    sbl = scontm(ii, l)
    etabl = sbl*opabl
    opac(l) = opabl
    scont(l) = sbl
    etac(l) = etabl
    Do k = 1, nf
      opabg(k, l) = opabl
      etabg(k, l) = etabl
    End Do
  End Do

  If (iov>0) Then

!   include overlapping lines

    del1 = vdop/vmax
    xmax = xmaxdop(ii)*del1
    deltax = 2.D0*xmax/dble(nf-1)
    cv = clight/vmax

    Do i = 1, iov
      ii1 = indexov(i)
      xlami = xlamcmf(ii1)
      xmax1 = xmaxdop(ii1)*vdopcmf(ii1)/vmax
      w = xlami/xlambda
      w1 = (w-1.D0)*cv

      flag = .False.
      Do l = 1, nd
        opabl = opalm(ii1, l)
        etabl = slinem(ii1, l)*opabl

        vdopl = (vdopcmf(ii1)**2-vturb**2)*t(l)/teff
        vdopl = sqrt(vdopl+vturb**2)
        del = vmax/vdopl

        Do k = 1, nf
          xk = xmax - (k-1)*deltax
          xki = w*xk + w1
          If (abs(xki)<xmax1) flag = .True.
          phikl = exp(-((xki*del)**2))*del/wpi
          opabg(k, l) = opabg(k, l) + opabl*phikl
          etabg(k, l) = etabg(k, l) + etabl*phikl
        End Do

      End Do
      If (.Not. flag) Then
        Print *, xlambda, xlami, ' SOMETHING WRONG IN !INTRINSIC OVERLAP'
      End If
    End Do
  End If


! lower boundary
  Call diffus(xlambda, t, r, nd, aic, dbdr)

! freq. grid and integr. weights
  Call fgrid(nf, nd, phi, pweight, deltax, vdop, vmax, teff, t, xmaxdop(ii))

  inversion = .False.

! check whether inversion, assuming worst case scenario:
! background pure cont

invcheck: Do l = 1, nd
    opakk = phi((nf+1)/2, l)*opal(l) + opac(l)
    If (opakk<0.D0 .And. .Not. inversion) Then
      Print *, ' INVERTED LEVELS FOR FOLLOWING TRANSITION '
      inversion = .True.
      Exit invcheck
    End If
  End Do invcheck



! blue wing boundary
! a)  first line, pure cont

  If (jj==1) Then

!   check consistency
    If (infotab(jj)%if/=0.) Then
      Write (999, *) ' STOP: ERROR IN INFOTAB(1)'
      Stop ' ERROR IN INFOTAB(1)'
    End If

    Call cont2(nd, np, r, p, z, scont, opac, aic, dbdr, corrfc, ubluwi)

!   b)  no intrinsic line overlap, only wind induced
  Else If (infotab(jj)%if==dble(nf)) Then
!   check consistency
    If (usave(nd,np-1,jj)/=dble(nf)) Then
      Write (999, *) ' STOP: ERROR IN U, V BW(1)'
      Stop ' ERROR IN U, V BW(1)'
    End If

    If (iov_only) Then
      Call cont2(nd, np, r, p, z, scont, opac, aic, dbdr, corrfc, ubluwi)
    Else

!     old version
!     CALL CONT2(ND,NP,R,P,Z,SCONT,OPAC,AIC,DBDR,CORRFC,UBLUWI)
!     CALL OVERLAP(JJ,VMAX,ND,UBLUWI,VBLUWI,OPAC,SCONT)
!     goto 100
!     ____________________________________________________________________________
!     new version, cmf transport from red wing (blue) to blue wing (red)

      cv = clight/vmax

      iir = infotab(jj)%index
      iib = infotab(jj)%indbw

      lred = xlamcmf(iir)
      lblu = xlamcmf(iib)

      xmaxr = xmaxdop(iir)*vdopcmf(iir)/vmax
      xmaxb = xmaxdop(iib)*vdopcmf(iib)/vmax
!     for tests, remember that XOV with respct to blue comp.
      xov = -infotab(jj)%xov*vdopcmf(iib)/vmax - xmaxb

      dew = lblu/lred
!     DELTA also with respect to blue comp.
      delta = (1.-dew)*cv - xmaxb - xmaxr*dew

      If (abs(1.-delta/xov)>1.D-6) Then
        Write (999, *) ' STOP:  DELTA AND XOV DIFFERENT'
        Stop '  DELTA AND XOV DIFFERENT'
      End If

      deltax1 = delta/dble(nf-1)

      Do l = 1, nd
        Do k = 1, nf
          opakk = opac(l)
          slinek(k, l) = scont(l)
          ak1(k, l) = 1.D0/opakk
          If (l/=1) az1(k, l-1) = 2.D0/(opakk+1.D0/ak1(k,l-1))
        End Do
      End Do

      Do jp = 1, np - 1
        lmax = min0(np+1-jp, nd)
        lz = lmax - 1

!       bluewing boundary condition from previous line

        ucmf(1:lmax) = usave(1:lmax, jp, jj)
        vcmf(1:lz) = vsave(1:lz, jp, jj)

        Do l = 1, lmax
          pp(l) = pp1(l, jp)/deltax1
        End Do

        Call raycont(jp, z(:,jp), r, nd, np, ucmf, vcmf, opac, etac, &
          aic, dbdr, corrfc, tauz1(:,jp), tau1(:,jp), pp)

        ubluwi(1:lmax, jp) = ucmf
        vbluwi(1:lz, jp) = vcmf

      End Do
!     100 continue
!     _____________________________________________________________________________

    End If

!   c)  intrinsic line overlap
  Else
!   check consistency
    w = int(infotab(jj)%if) + 1.D0 - infotab(jj)%if
    If (usave(nd,np-1,jj)/=w) Then
      Write (999, *) ' STOP: ERROR IN U, V BW(2)'
      Stop ' ERROR IN U, V BW(2)'
    End If
!   W=1.D0-W
!   numerically consistent to last digits (calculated as in subroutine RAY)
    w = infotab(jj)%if + 1.D0 - (int(infotab(jj)%if)+1)
    If (vsave(nd-1,np-1,jj)/=w) Then
      Write (999, *) ' STOP: ERROR IN U, V BW(3)'
      Stop ' ERROR IN U, V BW(3)'
    End If
  End If


  If (.Not. inversion) Then
!   STANDARD TREATMENT
    Do l = 1, nd
      Do k = 1, nf
        opalk = phi(k, l)*opal(l)
        etak = sline(l)*opalk + etabg(k, l)
        opakk = opalk + opabg(k, l)
        slinek(k, l) = etak/opakk
        ak1(k, l) = 1.D0/opakk
        If (l/=1) az1(k, l-1) = 2.D0/(opakk+1.D0/ak1(k,l-1))
      End Do
    End Do

  Else
!   INVERSION

    Do l = 1, nd
      Do k = 1, nf
        opalk = phi(k, l)*opal(l)
        etak = sline(l)*opalk + etabg(k, l)
        opakk = opalk + opabg(k, l)

!       ATTENTION: for inversion, ak1 is now the actual opacity,
!       and slinek is just the emisivity

        ak1(k, l) = opakk
        slinek(k, l) = etak
        If (l/=1) az1(k, l-1) = (opakk+ak1(k,l-1))/2.D0
      End Do
    End Do

  End If

  ajlmean = 0.D0
  alo = 0.D0

jploop2: Do jp = 1, np - 1
    lmax = min0(np+1-jp, nd)
    lz = lmax - 1

    If (jj==1 .Or. (iov_only .And. infotab(jj)%if==dble(nf))) Then

!     bluewing boundary condition from cont

      ucmf(1) = ubluwi(1, jp)

      Do l = 1, lz
        ucmf(l+1) = ubluwi(l+1, jp)
        aux = .5D0*(opac(l)+opac(l+1))
        vcmf(l) = (ucmf(l+1)-ucmf(l))/aux/(z(l,jp)-z(l+1,jp))
      End Do

    Else If (infotab(jj)%if==dble(nf)) Then

!     bluewing boundary condition from routine OVERLAP

      ucmf(1) = ubluwi(1, jp)
      Do l = 1, lz
        ucmf(l+1) = ubluwi(l+1, jp)
        vcmf(l) = vbluwi(l, jp)
      End Do
    Else

!     bluewing boundary condition from previous line

      ucmf(1:lmax) = usave(1:lmax, jp, jj)
      vcmf(1:lz) = vsave(1:lz, jp, jj)
    End If

    Do l = 1, lmax
      pp(l) = pp1(l, jp)/deltax
    End Do

    If (.Not. inversion) Then
!     STANDARD SOLUTION OVER TAU
      Call ray(.False., ii, nov, jp, z(:,jp), r, nd, np, ucmf, vcmf, sline, &
        phi, pweight, deltax, opac, etac, opal, aic, dbdr, corrfc, w0(:,jp), &
        w0last, ajlmean, alo, tauz1(:,jp), tau1(:,jp), pp, xlambda)

    Else
!     NEW FORMULATION OVER Z
      Call ray_inv(.False., ii, nov, jp, z(:,jp), r, nd, np, ucmf, vcmf, &
        pweight, deltax, opac, etac, opal, aic, dbdr, corrfc, &
        w0(:,jp), ajlmean, tauz1(:,jp), tau1(:,jp), pp, xlambda)
    End If

  End Do jploop2

! NOTE: SINCE ALO IS NOT MODIFIED IN THE INVERTED CASE, WE SIMPLY DON'T TOUCH
! IT


! DO L=1,ND
! PRINT*,l,' ',sline(l),' ',ajlmean(l),'  ',alo(l),'  ',opal(l)
! END DO

  Return

End Subroutine

!-----------------------------------------------------------------------

Subroutine overlap(jj, vmax, nd, ubluwi, vbluwi, opac, scont)

! transforms red wing radiation from line IIB to blue wing radiation for
! line IIR
! UBLUWI and OPAC are CONTINUUM quantities at line IIR

  Use :: nlte_type
  Use :: fund_const, Only: clight
  Use :: nlte_dim

  Use :: nlte_var, Only: xlamcmf, vdopcmf, xmaxdop

  Use :: cmf_multi, Only: lto, infotab, lmaxjp, fmax, fmin, uv, zray, usave, &
    vsave, iplus, iminus, iplus_cont, iminus_cont

  Implicit None

  Integer (i4b), Parameter :: nd1 = id_ndept, np = id_npoin

  Integer (i4b) :: jj, nd, it
  Real (dp) :: vmax
  Real (dp), Dimension (nd, np-1) :: ubluwi
  Real (dp), Dimension (nd-1, np-1) :: vbluwi
  Real (dp), Dimension (nd) :: opac, scont

  Real (dp) :: u(nd1), v(nd1-1), zp(lto), xip(lto), xim(nd1), uvj(lto), &
    ucont(nd1), vcont(nd1-1), xipc(nd1), ximc(nd1), aux(nd1), vaux(nd1), &
    ximold(nd1), opap(lto), sconp(lto), dtau(lto-1), fsp(lto-1), fsm(2:nd)

  Integer (i4b) :: iir, iib, jp, l, ll, lmax, ltot, lz, lminj, lj
  Real (dp) :: cv, lred, lblu, xmaxr, xmaxb, dew, delta, xov, xiplus, ximinus, &
    vnew, fminj, fmaxj, uvr, f, del, auxc, cplus, cminus, umin, sigm, sig1, &
    dt, e0, e1, w1, w2


  Logical :: core

  cv = clight/vmax

  iir = infotab(jj)%index
  iib = infotab(jj)%indbw

  lred = xlamcmf(iir)
  lblu = xlamcmf(iib)

  xmaxr = xmaxdop(iir)*vdopcmf(iir)/vmax
  xmaxb = xmaxdop(iib)*vdopcmf(iib)/vmax
! for tests, remember that XOV with respct to blue comp.
  xov = -infotab(jj)%xov*vdopcmf(iib)/vmax - xmaxb


  dew = lblu/lred
! DELTA also with respect to blue comp.
  delta = (1.-dew)*cv - xmaxb - xmaxr*dew

  If (abs(1.-delta/xov)>1.D-6) Then
    Write (999, *) ' STOP: DELTA AND XOV DIFFERENT'
    Stop '  DELTA AND XOV DIFFERENT'
  End If

  Do jp = 1, np - 1
    lmax = lmaxjp(jp)
    lz = lmax - 1
    core = lmax == nd
    ltot = 2*lmax - 1
    If (core) ltot = nd
    uvj(1:ltot) = uv(1:ltot, jp)
    zp(1:ltot) = zray(1:ltot, jp)

!   ------------------------------------------
!   prepare cont-transport for given p-ray
    If (core) Then
      opap(1:ltot) = opac(1:ltot)
      sconp(1:ltot) = scont(1:ltot)
    Else
      opap(1:lmax) = opac(1:lmax)
      sconp(1:lmax) = scont(1:lmax)
      Do l = lmax + 1, ltot
        ll = 2*lmax - l
        opap(l) = opac(ll)
        sconp(l) = scont(ll)
      End Do
    End If

    Do l = 1, ltot - 1
      dtau(l) = 0.5D0*(opap(l)+opap(l+1))*(zp(l)-zp(l+1))
      dt = dtau(l)
      If (dt<1.D-8) Then
        w1 = .5D0*dt
        w2 = .5D0*dt
      Else
        e0 = exp(-dt)
        e1 = (1.D0-e0)/dt
        w1 = 1.D0 - e1
        w2 = e1 - e0
      End If
      fsp(l) = sconp(l)*w1 + sconp(l+1)*w2
      If (l<=lmax-1) fsm(l+1) = sconp(l+1)*w1 + sconp(l)*w2
    End Do
!   ------------------------------------------

!   calculate I_plus and I_minus (CONT) at blue wing of red line from
!   UBLUWI, accounting for consistent boundary conditions

    ucont(1) = ubluwi(1, jp)

    Do l = 1, lz
      ucont(l+1) = ubluwi(l+1, jp)
      auxc = .5D0*(opac(l)+opac(l+1))
      vcont(l) = (ucont(l+1)-ucont(l))/auxc/(zp(l)-zp(l+1))
    End Do

    ximinus = iminus_cont(jp)
    xipc(1) = 2.D0*ucont(1) - ximinus
    ximc(1) = ximinus

    If (core) Then
      xiplus = iplus_cont(jp)
      xipc(lmax) = xiplus
      vnew = 2.D0*ucont(lmax) - xiplus
      If (vnew>=0.) Then
        ximc(lmax) = vnew
      Else
        ximc(lmax) = 0.
      End If
    Else
      ximc(lmax) = ucont(lmax)
      xipc(lmax) = ximc(lmax)
    End If

!   iterative solution for I_plus and I minus, &
!   assuming I(l+1/2)=sqrt(I(l)*I(l+1)), see notes.

    If (lz>=2) Then
!     with this choice, we obtain the best results in cases of numerical
!     problems
      umin = 1.D-6*minval(ucont(1:lmax))

!     start values
      Do l = 2, lz
        xipc(l) = max(ucont(l)+vcont(l), umin)
        ximc(l) = max(ucont(l)-vcont(l), umin)
      End Do
      ximold(2:lz) = ximc(2:lz)

!     iteration
      Do it = 1, 20
        Do l = 2, lz
          cplus = sqrt(xipc(l+1)/xipc(l))
          cminus = sqrt(ximc(l+1)/ximc(l))
          auxc = 2.D0/(cplus+cminus)
          xipc(l) = auxc*(cminus*ucont(l)+vcont(l))
          xipc(l) = max(xipc(l), umin)
          ximc(l) = auxc*(cplus*ucont(l)-vcont(l))
          ximc(l) = max(ximc(l), umin)
        End Do
        sigm = sqrt(dot_product((ximc(2:lz)-ximold(2:lz)),(ximc(2:lz)- &
          ximold(2:lz))))
        If (it==2) Then
          sig1 = sigm
        Else If (it>2) Then
          If (sigm<.05*sig1) Exit
        End If
        ximold(2:lz) = ximc(2:lz)
      End Do
    End If

!   test for continuum: solution should be "perfect" in u
    Do l = 1, lmax
      ximold(l) = .5D0*(xipc(l)+ximc(l))
    End Do
    sig1 = maxval(abs(1.-ximold(1:lmax)/ucont(1:lmax)))
    If (sig1>.003D0) Then
      Print *, ' PROBLEM IN CONT-INTERPOLATION AT JP =', jp, sig1
      Do l = 1, lz
        vaux(l) = .5D0*(sqrt(xipc(l)*xipc(l+1))-sqrt(ximc(l)*ximc(l+1)))
        Print *, l, ximold(l), ucont(l), vaux(l), vcont(l)
      End Do
      l = lmax
      Print *, l, ximold(l), ucont(l)
    End If

!   ---
!   calculate I_plus and I_minus at red wing of blue line from u and v, &
!   accounting for consistent boundary conditions

    u(1:lmax) = usave(1:lmax, jp, jj)
    v(1:lz) = vsave(1:lz, jp, jj)
!   with this choice, we obtain the best results in cases of numerical
!   problems.
!   Note, that U may have been set to zero in CMF.
    umin = 1.D-6*minval(u(1:lmax), mask=u(1:lmax)>0.D0)
    If (umin==0.D0) Then
      Write (999, *) ' STOP: UMIN = 0 IN OVERLAP'
      Stop ' UMIN = 0 IN OVERLAP'
    End If
    Where (u(1:lmax)==0D0) u(1:lmax) = umin

!   boundary conditions
    ximinus = iminus(jp, jj)
    xip(1) = 2.D0*u(1) - ximinus
    xim(1) = ximinus

    If (core) Then
      xiplus = iplus(jp, jj)
      xip(lmax) = xiplus
      vnew = 2.D0*u(lmax) - xiplus
      If (vnew>=0.) Then
        xim(lmax) = vnew
      Else
        xim(lmax) = 0.
      End If
    Else
      xim(lmax) = u(lmax)
      xip(lmax) = xim(lmax)
    End If

!   iterative solution for I_plus and I minus for non-boundary points, &
!   assuming I(l+1/2)=sqrt(I(l)*I(l+1)), see notes.
    If (lz>=2) Then

!     start values
      Do l = 2, lz
        xip(l) = max(u(l)+v(l), umin)
        xim(l) = max(u(l)-v(l), umin)
      End Do
      ximold(2:lz) = xim(2:lz)

!     iteration
      Do it = 1, 20
        Do l = 2, lz
          cplus = sqrt(xip(l+1)/xip(l))
          cminus = sqrt(xim(l+1)/xim(l))
!         print*,it,jp,l,xip(l+1),xip(l),xim(l+1),xim(l),umin !for tests
          auxc = 2.D0/(cplus+cminus)
          xip(l) = auxc*(cminus*u(l)+v(l))
          xip(l) = max(xip(l), umin)
          xim(l) = auxc*(cplus*u(l)-v(l))
          xim(l) = max(xim(l), umin)
        End Do
        sigm = sqrt(dot_product((xim(2:lz)-ximold(2:lz)),(xim(2:lz)- &
          ximold(2:lz))))
        If (it==2) Then
          sig1 = sigm
        Else If (it>2) Then
          If (sigm<.05*sig1) Exit
        End If
        ximold(2:lz) = xim(2:lz)
      End Do
    End If

!   test for lines: solution should be "OK" in u
    Do l = 1, lmax
      ximold(l) = .5D0*(xip(l)+xim(l))
    End Do
    sig1 = maxval(abs(1.-ximold(1:lmax)/u(1:lmax)), mask=u(1:lmax)/=umin)
    If (sig1>.1D0) Then
      Print *, ' PROBLEM IN LINE-INTERPOLATION AT JP =', jp, sig1
!     for tests only
!     DO L=1,LZ
!     VAUX(L)=.5D0*(SQRT(XIP(L)*XIP(L+1))-SQRT(XIM(L)*XIM(L+1)))
!     PRINT*,L,XIMOLD(L),U(L),VAUX(L),VSAVE(L,JP,JJ)
!     ENDDO
!     L=LMAX
!     PRINT*,L,XIMOLD(L),USAVE(L)
    End If

    If (.Not. core) Then
      Do l = 1, lmax - 1
        ll = 2*lmax - l
        xip(ll) = xim(l)
      End Do
    End If

!   ---
!   find out which intensities at red wing (blue comp.) correspond to
!   blue wing intensities for red comp.)
!   In case of no interfering continuum, both intensities are equal

    lminj = 1
    fminj = fmin(jp)
    fmaxj = fmax(jp)

plus: Do l = 1, lmax

      uvr = uvj(l)
      f = uvr - delta
      lj = l
      If (f<fminj) Then

        Do ll = lj, lmax
          aux(ll) = xipc(ll)
!         in this case, Iplus is totally determined by the local blue-wing
!         continuum
!         PRINT*,'PLUS',JP,LL,XIPC(LL)
        End Do
        Exit plus

      Else

        Do ll = lminj, ltot - 1

          If (f<uvj(ll) .And. f>=uvj(ll+1)) Then
            del = (f-uvj(ll))/(uvj(ll+1)-uvj(ll))
!           del corresponds to both uv and z, since assumed as linear in both
!           cases
            If (xip(ll+1)==0.D0) Then
              If (ll+1==ltot) Then ! outer boundary OK
                xiplus = xip(ll)*(1.D0-del) + xip(ll+1)*del
              Else
                Write (999, *) ' STOP: XIP EQ 0 AND LL+1 NE LTOT!'
                Stop ' XIP EQ 0 AND LL+1 NE LTOT!'
              End If
            Else
              xiplus = xip(ll)**(1.D0-del)*xip(ll+1)**del
            End If
!           xiplus input and output!
            Call formacon_plus(jp, l, ll, del, xiplus, xipc(l), zp, &
              opap, sconp, dtau, fsp)
            aux(l) = xiplus
!           PRINT*,'PLUS',JP,L,UVR,LL,LL+1,AUX(L)
!           This is the incident intensity(plus) which is additionally
!           modified by continuum procecesses between zmin(f) and z(l)
!           to yield the final blue wing intensity(plus) at z(l)
            lminj = ll
            Cycle plus
          End If

        End Do

        Print *, ' AT XBLUE = ', lblu, ' AND DELTA = ', delta
        Write (999, *) ' STOP: UVBLUE(XIP) NOT FOUND!'
        Stop ' UVBLUE(XIP) NOT FOUND!'

      End If

    End Do plus


    u(1:lmax) = .5D0*aux(1:lmax) !  0.5 Iplus
    vaux(1:lmax) = aux(1:lmax) !    Iplus

minus: Do l = 1, lmax

      uvr = uvj(l)
      f = uvr + delta
      lj = l
      If (f>=fmaxj) Then
        aux(l) = ximc(l)
!       in this case, Iminus is totally determined by the local blue-wing
!       continuum
!       PRINT*,'MINUS',JP,L,XIMC(L)
        Cycle minus
      End If

      Do ll = 1, lj - 1

        If (f<uvj(ll) .And. f>=uvj(ll+1)) Then
          del = (f-uvj(ll))/(uvj(ll+1)-uvj(ll))
          If (xim(ll)==0.D0) Then
            If (ll==1) Then !       outer boundary OK
              ximinus = xim(ll)*(1.D0-del) + xim(ll+1)*del
            Else
              Write (999, *) ' STOP: XIM EQ 0 AND LL NE 1!'
              Stop ' XIM EQ 0 AND LL NE 1!'
            End If
          Else
            ximinus = xim(ll)**(1.D0-del)*xim(ll+1)**del
          End If
!         ximinus input and output!
          Call formacon_minus(jp, l, ll, del, ximinus, ximc(l), zp, &
            opap, sconp, dtau, fsm)
          aux(l) = ximinus
!         PRINT*,'MINUS',JP,L,UVR,LL,LL+1,AUX(L)
!         This is the incident intensity(minus) which is additionally
!         modified by continuum procecesses between zmax(f) and z(l)
!         to yield the final blue wing intensity(minus) at z(l)
          Cycle minus
        End If

      End Do

      Print *, ' AT XBLUE = ', lblu, ' AND DELTA = ', delta
      Write (999, *) ' STOP: UVBLUE(XIM) NOT FOUND!'
      Stop ' UVBLUE(XIM) NOT FOUND!'

    End Do minus

    ubluwi(1:lmax, jp) = u(1:lmax) + .5D0*aux(1:lmax) ! .5 (Iplus + Iminus)

!   from construction (see above); here, we simply interpolate between
!   continuum quantities and line quantities at the border of influence
    Do l = 1, lz
      vbluwi(l, jp) = .5D0*(sqrt(vaux(l)*vaux(l+1))-sqrt(aux(l)*aux(l+1)))
    End Do

  End Do

  Return
End Subroutine

!-----------------------------------------------------------------------

Subroutine formacon_plus(jp, l, ll, del, xiplus, xipcont, zp, opap, &
  sconp, dtau, fsp)

  Use :: nlte_type
  Use :: nlte_dim
  Use :: cmf_multi, Only: lto

  Implicit None

! Integer (i4b), Parameter :: nd = id_ndept

  Integer (i4b) :: jp, l, ll
  Real (dp) :: del, xiplus, xipcont
  Real (dp), Dimension (lto) :: zp, opap, sconp
  Real (dp), Dimension (lto-1) :: dtau, fsp

  Integer (i4b) :: i
  Real (dp), Dimension (l:ll+1) :: taucon
  Real (dp) :: zlast, opalast, slast, ipc, dt, e0, e1, w1, w2

  If (l>=ll+1) Then
    Write (999, *) ' STOP: ERROR IN FORMACON_PLUS'
    Stop ' ERROR IN FORMACON_PLUS'
  End If

! last, shorter step
  opalast = opap(ll)*(1.D0-del) + opap(ll+1)*del
  slast = sconp(ll)*(1.D0-del) + sconp(ll+1)*del
  zlast = zp(ll)*(1.D0-del) + zp(ll+1)*del

  If (zlast>zp(ll) .Or. zlast<zp(ll+1)) Then
    Print *, jp, ll, zlast, zp(ll), zp(ll+1)
    Write (999, *) ' STOP: ERROR IN ZPLUS'
    Stop ' ERROR IN ZPLUS'
  End If

  taucon(l) = 0.D0
  Do i = l + 1, ll
    taucon(i) = taucon(i-1) + dtau(i-1) ! note DTAU defined at "head"
  End Do
  taucon(ll+1) = taucon(ll) + .5D0*(opap(ll)+opalast)*(zp(ll)-zlast)

  If (taucon(ll+1)>=10.D0) Then
!   interfering continuum optically thick, use continuum value for consistency
    xiplus = xipcont
    Return
  End If

  ipc = exp(-taucon(ll+1))*xiplus ! incident intensity


! first, shorter step
  dt = taucon(ll+1) - taucon(ll)

  If (dt<1.D-8) Then
    w1 = .5D0*dt
    w2 = .5D0*dt
  Else
    e0 = exp(-dt)
    e1 = (1.D0-e0)/dt
    w1 = 1.D0 - e1
    w2 = e1 - e0
  End If
  ipc = ipc + exp(-taucon(ll))*(sconp(ll)*w1+slast*w2)


  Do i = ll - 1, l, -1
    ipc = ipc + exp(-taucon(i))*fsp(i)
  End Do

! consistency check
! IF(TAUCON(LL+1).LT.1.D-3) THEN
! IF(ABS(1.D0-IPC/XIPLUS).GT..02) STOP ' PROBLEMS IN FORMACON_PLUS'
! ENDIF

  xiplus = ipc
  Return
End Subroutine


!-----------------------------------------------------------------------

Subroutine formacon_minus(jp, l, ll, del, ximinus, ximcont, zp, opap, &
  sconp, dtau, fsm)

  Use :: nlte_type
  Use :: nlte_dim
  Use :: cmf_multi, Only: lto

  Implicit None

  Integer (i4b), Parameter :: nd = id_ndept

  Integer (i4b) :: jp, l, ll
  Real (dp) :: del, ximinus, ximcont
  Real (dp), Dimension (lto) :: zp, opap, sconp
  Real (dp), Dimension (lto-1) :: dtau
  Real (dp), Dimension (2:nd) :: fsm

  Integer (i4b) :: i
  Real (dp), Dimension (ll:l) :: taucon
  Real (dp) :: zlast, opalast, slast, imc, dt, e0, e1, w1, w2

  If (l<=ll) Then
    Write (999, *) ' STOP: ERROR IN FORMACON_MINUS'
    Stop ' ERROR IN FORMACON_MINUS'
  End If

! last, shorter step
  opalast = opap(ll)*(1.D0-del) + opap(ll+1)*del
  slast = sconp(ll)*(1.D0-del) + sconp(ll+1)*del
  zlast = zp(ll)*(1.D0-del) + zp(ll+1)*del

  If (zlast>zp(ll) .Or. zlast<zp(ll+1)) Then
    Print *, jp, ll, zlast, zp(ll), zp(ll+1)
    Write (999, *) ' STOP: ERROR IN ZMINUS'
    Stop ' ERROR IN ZMINUS'
  End If

  taucon(l) = 0.D0
  Do i = l - 1, ll + 1, -1
    taucon(i) = taucon(i+1) + dtau(i) ! note DTAU defined at "head"
  End Do
  taucon(ll) = taucon(ll+1) + .5D0*(opap(ll+1)+opalast)*(zlast-zp(ll+1))

  If (taucon(ll)>=10.D0) Then
!   interfering continuum optically thick, use continuum value for consistency
    ximinus = ximcont
    Return
  End If

  imc = exp(-taucon(ll))*ximinus !  incident intensity

! first, shorter step
  dt = taucon(ll) - taucon(ll+1)

  If (dt<1.D-8) Then
    w1 = .5D0*dt
    w2 = .5D0*dt
  Else
    e0 = exp(-dt)
    e1 = (1.D0-e0)/dt
    w1 = 1.D0 - e1
    w2 = e1 - e0
  End If
  imc = imc + exp(-taucon(ll+1))*(sconp(ll+1)*w1+slast*w2)

  Do i = ll + 2, l
    imc = imc + exp(-taucon(i))*fsm(i)
  End Do

! consistency check not performed, since I_minus oscillatory (and
! sometimes much smaller than scont, so that IMC different from
! XIMINUS even for low TAUCON(LL)

  ximinus = imc
  Return
End Subroutine

!-----------------------------------------------------------------------


Subroutine cont2(nd, np, r, p, z, scont, opa, xic1, xic2, corr, u)

  Use :: nlte_type
  Use :: nlte_dim

  Use :: cmf_multi, Only: cmf_all_bluwi

  Implicit None

! ......................................................................

! blue-wing intensities from pure formal solution,
! only u and j are calculated

! .. parameters ..
  Integer (i4b), Parameter :: nd1 = id_ndept, np1 = id_npoin
! ..
! .. scalar arguments ..
  Real (dp) :: corr, xic1, xic2
  Integer (i4b) :: nd, np
! ..
! .. array arguments ..
  Real (dp) :: opa(nd1), p(np1), r(nd1), scont(nd1), u(nd1, np1-1), &
    z(nd1, np1)
! ..
! .. local scalars ..
  Real (dp) :: pi2, thom1
  Integer (i4b) :: jp, l, lmax
  Logical :: core
! ..
! .. local arrays ..
  Real (dp) :: akp(nd1), ta(nd1), taum(np1), tb(nd1), tc(nd1)
! ..
! .. external subroutines ..
  External :: invtri, setup2
! ..
! .. intrinsic functions ..
! INTRINSIC ACOS,ATAN,MIN0
! ..

  pi2 = acos(0.D0)
  thom1 = 1.D0

jploop: Do jp = 1, np - 1

!   ---  implicit assumption: thomson(1):=1

    If (cmf_all_bluwi) Then
      taum(jp) = opa(1)*r(1)/3.D0

    Else

      If (jp==1) Then
        taum(1) = opa(1)*r(1)*(1.D0/3.D0+2.D0*thom1/3.D0)
      Else
        taum(jp) = opa(1)*r(1)*r(1)*(thom1/p(jp)*(pi2-atan(z(1,jp)/p(jp)))+( &
          1.D0-thom1)*r(1)*r(1)/2.D0/p(jp)**3*(pi2-atan(z(1,jp)/p(jp))-z(1, &
          jp)*p(jp)/r(1)/r(1)))
      End If

    End If
    If (taum(jp)<0.D0) Then
      Write (999, *) ' STOP: TAUM NEGATIVE!'
      Stop ' TAUM NEGATIVE!'
    End If

    lmax = min0(np+1-jp, nd)
    core = (lmax==nd)

    Call setup2(jp, lmax, core, z(:,jp), scont, opa, xic1, xic2, corr, akp, &
      ta, tb, tc, taum(jp))

    Call invtri(ta, tb, tc, akp, lmax)

    Do l = 1, lmax
      u(l, jp) = akp(l)
    End Do

  End Do jploop

  Return

End Subroutine

!-----------------------------------------------------------------------

Subroutine setup2(jp, lmax, core, z, scont, opa, xic1, xic2, corr, akp, ta, &
  tb, tc, taum)

  Use :: nlte_type
  Use :: nlte_dim

  Use :: nlte_opt, Only: iov_only
  Use :: cmf_multi, Only: iplus_cont, iminus_cont



  Implicit None

! sets up matrix elements for subroutine cont2
! note: ta,tc negative compared to setup1


! .. parameters ..
  Integer (i4b), Parameter :: nd1 = id_ndept
! ..
! .. scalar arguments ..
  Real (dp) :: corr, taum, xic1, xic2
  Integer (i4b) :: jp, lmax
  Logical :: core
! ..
! .. array arguments ..
  Real (dp) :: akp(nd1), opa(nd1), scont(nd1), ta(nd1), tb(nd1), tc(nd1), &
    z(nd1)
! ..
! .. local scalars ..
  Real (dp) :: ak, akdz, bb, cc, dt0, dtm, dtp, dz, e0, e1
  Integer (i4b) :: l, lz
! ..
! .. intrinsic functions ..
! INTRINSIC EXP
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
    akp(1) = -scont(1)*(bb+e1) - scont(2)*cc
    tb(1) = -1.D0 - dtp - bb
    tc(1) = -dtp + cc
  Else
    akp(1) = -scont(1)*e1
    tb(1) = -1.D0 - dtp
    tc(1) = -dtp
  End If

  If (.Not. iov_only) iminus_cont(jp) = scont(1)*e1

! non boundary points

  Do l = 2, lz
    dtm = dtp
    ak = .5D0*(opa(l)+opa(l+1))
    dz = z(l) - z(l+1)
    dtp = 1.D0/ak/dz
    dt0 = 2.D0/(1.D0/dtm+1.D0/dtp)
    akp(l) = -scont(l)
    ta(l) = -dt0*dtm
    tc(l) = -dt0*dtp
    tb(l) = -dt0*(dtm+dtp) - 1.D0
  End Do

  l = lmax

! inner boundary,2nd order

  If (core) Then
    ta(l) = -dtp
    akp(l) = -xic1 - z(l)*xic2/opa(l)*corr
    tb(l) = -dtp - 1.D0
    If (.Not. iov_only) iplus_cont(jp) = -akp(l)
  Else
    akp(l) = -scont(l)
    ta(l) = -2.D0*dtp*dtp
    tb(l) = -2.D0*dtp*dtp - 1.D0
  End If

  Return

End Subroutine

!-----------------------------------------------------------------------

Subroutine fgrid(nf, nd, phi, pweight, deltax, vdop, vmax, teff, t, xmax0)

  Use :: nlte_type
  Use :: nlte_dim
  Use :: fund_const, Only: wpi => sqpi
  Use :: nlte_var, Only: vturb
  Implicit None


! calculates frequency scale,lineprofile phi and
! renormalized integration weights pweight


! ..
! .. scalar arguments ..
  Real (dp) :: deltax, teff, vdop, vmax, xmax0
  Integer (i4b) :: nd, nf
! ..
! .. array arguments ..
  Real (dp) :: phi(nf, nd), pweight(nf, nd), t(nd)

! .. local scalars ..
  Real (dp) :: del, del1, vdopl, ws, xk, xmax
  Integer (i4b) :: k, l
! ..

! print*,' xmax0 = ',xmax0

  del1 = vdop/vmax
  xmax = xmax0*del1
  deltax = 2.D0*xmax/dble(nf-1)

  Do l = 1, nd
    vdopl = (vdop**2-vturb**2)*t(l)/teff
    vdopl = sqrt(vdopl+vturb**2)
    del = vmax/vdopl
    ws = .0D0
    Do k = 1, nf
      xk = xmax - (k-1)*deltax
      pweight(k, l) = exp(-((xk*del)**2))
      ws = ws + pweight(k, l)
      phi(k, l) = pweight(k, l)*del/wpi
    End Do

!   renormalization

    Do k = 1, nf
      pweight(k, l) = pweight(k, l)/ws
    End Do

  End Do

  Return

End Subroutine

!-----------------------------------------------------------------------

Subroutine ray(sing, ii, nov, jp, z, r, nd, np, u, v, sline, phi, pweight, &
  deltax, opa, eta, opal, aic, dbdr, corrfc, w0, w0last, ajlmean, alo, tauz1, &
  tau1, pp, xlambda)

  Use :: nlte_type
  Use :: nlte_dim
  Use :: nlte_var, Only: aak1 => ak1, aaz1 => az1, slinek
  Use :: cmf_multi, Only: infotab, usave, vsave, iplus, iminus

  Implicit None


! line radiation transfer in the comoving frame from a
! given source function for a given impact parameter jp.
! the integration is carried out in space and freq. to
! yield the mean line intensity

! new version: inlinen aller algebraischen routinen

! j-bar

! .. parameters ..
  Integer (i4b), Parameter :: nd1 = id_ndept, nf = id_nfcmf
! ..
! .. scalar arguments ..
  Real (dp) :: aic, corrfc, dbdr, deltax, w0last, xlambda
  Integer (i4b) :: ii, nov, jp, nd, np
! ..
! .. array arguments ..
  Real (dp) :: ajlmean(nd), alo(nd), opa(nd), eta(nd), opal(nd), phi(nf, nd), &
    pp(nd1), pweight(nf, nd), r(nd), sline(nd), tau1(nd), tauz1(nd), u(nd), &
    v(nd-1), w0(nd), z(nd)

  Logical :: sing
! ..
! .. local scalars ..
  Real (dp) :: a, ak1, alo1, az1, b, daz, dazm, dbz, dbzm, dt, dtz, dtzm, dx, &
    dxz, dxz1, dxzm, dxzm1, erfck, erfcko, ff1, h1, help, pw, rmax, s1, sc, &
    sl, tauc, taus, tauz, tb1, w
  Integer (i4b) :: i, k, l, l1, lmax, lz, ni
! ..
! .. local arrays ..
  Real (dp) :: ddia(0:nd1), edia(nd1+1), ff(nd1), ga(nd1), h(nd1), qq(nd1), &
    s(nd1), ta(nd1), tb(nd1), tc(nd1), ub(nd1), uold(nd1), v1old(nd1), &
    va(nd1), vb(nd1), vold(nd1)
! ..
! .. intrinsic functions ..
! INTRINSIC EXP,MAX,MIN0
! ..

  lmax = min0(np+1-jp, nd)
  lz = lmax - 1
  rmax = r(1)

  erfcko = pweight(1, 1)

  uold(1:lmax) = 0.D0
  vold(1:lz) = 0.D0
  v1old(1:lz) = 0.D0

! until further evidence, use optically thick cont. without line background
! in case of multi-line transport

  tauc = opa(1)*rmax/3.D0
  sc = 0.D0
  alo1 = 0.D0

  If (tauc>=1.D0) sc = eta(1)/opa(1)

  taus = opal(1)/(pp(1)*deltax)
  sl = 0.D0


! loop for all original frequency points


kloop: Do k = 1, nf

    If (k==1) Go To 100

    erfck = erfcko + pweight(k, 1)

!   if(k.eq.nf.and.abs(erfck-1.d0).gt.1.d-10)
!   *   stop 'error in cmf: erfck'

!   ***  inlined subroutine cmfset
!   ----------------------------------------------------------------------

!   ***  outer boundary condition  -  2nd order

    ak1 = aak1(k, 1)
    az1 = aaz1(k, 1)
    dx = pp(1)*ak1

!   approximate treatment for optically thick continuum
!   (rho^2 processes)
!   respectively optically thick line at outer boundary


    ff1 = 0.D0
    If (taus>=.5D0) Then
      ff1 = 1.D0 - exp(-taus*erfck)
      sl = sline(1)*ff1

!     slp=sline(1)/(taus*deltax)*(exp(-taus*erfck)-exp(-taus*erfcko))
!     sl=sl-slp

    End If

    s(1) = max(sl, sc)
    tc(1) = tauz1(1)*ak1
    tb(1) = tc(1) + dx + 1.D0
    ub(1) = dx
    vb(1) = 0.D0
    dtzm = tauz1(1)*az1

!   ***  alo for outer boundary

    ff(1) = ff1
    If (s(1)==sc) ff(1) = 0.D0

!   ***  for g and h, the matrix element s are not different from
!   ***  inner points

    dxzm = .5D0*(pp(1)+pp(2))*az1
    dxzm1 = 1.D0/(1.D0+dxzm)
    dazm = dtzm*dxzm1
    dbzm = dxzm*dxzm1
    ga(1) = -dazm
    h(1) = dbzm

!   ***  non-boundary points

    Do l = 2, lz
      ak1 = aak1(k, l)
      az1 = aaz1(k, l)
      s(l) = slinek(k, l)
      dt = tau1(l)*ak1
      dtz = tauz1(l)*az1
      dx = pp(l)*ak1
      dxz = .5D0*(pp(l)+pp(l+1))*az1
      dxz1 = 1.D0/(1.D0+dxz)
      daz = dtz*dxz1
      dbz = dxz*dxz1
      ta(l) = dt*dazm
      tc(l) = dt*daz
      tb(l) = ta(l) + tc(l) + dx + 1.D0
      ff(l) = phi(k, l)*opal(l)*ak1
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

!     ***  inner boundary condition (core rays)  -  only to first order
!     ***  diffusion-approximation

      ak1 = aak1(k, l)
      s(l) = aic + z(l)*ak1*dbdr*corrfc
      dx = pp(l)*ak1
      dt = tauz1(l-1)*ak1
      ta(l) = dt
      tb(l) = dt + dx + 1.D0
      ff(l) = 0.D0
      ub(l) = dx
      va(l) = 0.0D0
      vb(l) = 0.0D0
    Else

!     ***  inner boundary condition (non-core rays)  -  second order

      ak1 = aak1(k, l)
      s(l) = slinek(k, l)
      tauz = z(lz)
      dt = ak1/tauz
      dx = pp(l)*ak1
      ta(l) = 2.D0*dt*dazm !        instead of daz
      tb(l) = ta(l) + dx + 1.D0
      ub(l) = dx
      va(l) = -2.D0*dt*dbzm !       instead of dbz
      ff(l) = phi(k, l)*opal(l)*ak1
    End If

!   -----------------------------------------------------------------------
!   ***  continuation of subroutine ray

    erfcko = erfck

    qq(1) = 0.D0

    Do l = 2, lz
      qq(l) = va(l)*v1old(l-1) + vb(l)*vold(l)
    End Do

    qq(lmax) = va(lmax)*v1old(lz)

!   uold=ub * uold + qq + ff

    Do l = 1, lmax
      uold(l) = ub(l)*uold(l) + qq(l) + ff(l)
    End Do

!   ***     now, uold is the solution vector for the inhomgeneous
!   ***     equation corresponding to u(sl=1)-u(sl=0)

!   inlinen von diag mit option 1


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

!   ***  now,uold is the solution corresponding to u(sl=1)-u(0)

    Do l = 2, lz
      l1 = l - 1
      v1old(l1) = -ga(l1)*uold(l) + h(l1)*v1old(l1)
      vold(l) = ga(l)*uold(l) + h(l)*vold(l)
    End Do

    l = lz
    v1old(l) = -ga(l)*uold(lmax) + h(l)*v1old(l)

!   ***  now, v1old, vold are the v's corresponding to delta u

!   inlinen von vmalv

    b = v(1)
    qq(1) = vb(1)*b


    Do l = 2, lz
      a = b
      b = v(l)
      qq(l) = va(l)*a + vb(l)*b
    End Do

    qq(lmax) = va(lmax)*b

!   u = ub * u + qq + s

    Do l = 1, lmax
      u(l) = ub(l)*u(l) + qq(l) + s(l)
    End Do

!   ***     inlining of invtri

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
      u(ni) = u(ni) + tc(ni)*u(ni+1)
    End Do

    Do l = 1, lmax
      uold(l) = max(uold(l), 0.D0)

!     i don't like the following statement at all, but as one says
!     in germany:  man hat schon pferde kotzen sehen.
!     actually, i found one case where due to a very strange behaviour
!     of occupation numbers v became strongly negative and in the
!     course also qq, such that u became negative. once u is negative,
!     this can no longer be compensated. hence:

      u(l) = max(u(l), 0.D0)
    End Do

!   now u is the new field at index k

!   inlinen von gmalu und
!   v = h * v + gmalu

    b = u(1)
    Do l = 1, lz
      a = b
      b = u(l+1)
      v(l) = h(l)*v(l) + ga(l)*(a-b)
    End Do

!   now v is the new field at index k

!   adding the new u to ajlmean and to the 0.and 2nd moments aj,ak


    u(lmax+1:nd) = 0.D0
    v(lz+1:nd-1) = 0.D0

100 Continue

    Do l = 1, lmax
      pw = pweight(k, l)*w0(l)
      ajlmean(l) = ajlmean(l) + u(l)*pw
      alo(l) = alo(l) + uold(l)*pw
    End Do

    If (sing) Cycle kloop

!   store u,v if necessary

    Do i = 1, nov

      If (infotab(i)%indbw==ii) Then
        If (int(infotab(i)%if)==nf .And. k==nf) Then
          usave(1:lmax, jp, i) = u(1:lmax)
          vsave(1:lz, jp, i) = v(1:lz)
          usave(nd, np-1, i) = nf
          iminus(jp, i) = s(1)
          If (lmax==nd) iplus(jp, i) = s(nd)
!         IF(JP.EQ.1) PRINT*,II,'NF STORED',INFOTAB(I)%IF,' FOR L<INE ',I,'
!         ',NF
        Else If (int(infotab(i)%if)==k) Then
          w = k + 1.D0 - infotab(i)%if
          usave(1:lmax, jp, i) = w*u(1:lmax)
          vsave(1:lz, jp, i) = w*v(1:lz)
          usave(nd, np-1, i) = w
!         IF(JP.EQ.1) PRINT*,II,K,' STORED (FIRST)',INFOTAB(I)%IF,' FOR LINE
!         ',I,' ',W
        Else If (int(infotab(i)%if)==k-1) Then
          w = infotab(i)%if + 1.D0 - k
          usave(1:lmax, jp, i) = usave(1:lmax, jp, i) + w*u(1:lmax)
          vsave(1:lz, jp, i) = vsave(1:lz, jp, i) + w*v(1:lz)
          vsave(nd-1, np-1, i) = w
!         IF(JP.EQ.1) PRINT*,II,K,' STORED (2ND)',INFOTAB(I)%IF,' FOR LINE
!         ',I,' ',W
        End If
      End If

    End Do


  End Do kloop

  If (jp/=np-1) Return

! ***  final value for p = rmax and z=0 in optically thick case!!

  taus = opal(1)/(pp(1)*deltax)

  If (sc/=0.D0 .Or. taus>=.5D0) Then
    Print *, ' OPTICALLY THICK OUTER BOUNDARY AT', xlambda
    erfcko = pweight(1, 1)

    Do k = 2, nf
      erfck = erfcko + pweight(k, 1)

!     if(k.eq.nf.and.abs(erfck-1.d0).gt.1.d-10)
!     *stop 'error in cmf: erfck'

      alo1 = 1.D0 - exp(-taus*erfck)
      sl = sline(1)*alo1

!     slp=sline(1)/(taus*deltax)*(exp(-taus*erfck)-exp(-taus*erfcko))
!     sl=sl-slp

      s1 = max(sl, sc)
      pw = pweight(k, 1)*w0last
      ajlmean(1) = ajlmean(1) + s1*pw
      If (s1==sl) alo(1) = alo(1) + alo1*pw
      erfcko = erfck
    End Do
  End If

! if(xlambda.gt.303..and.xlambda.lt.304.) then
! do l=1,nd
! write(*,100) l,alo(l),ajlmean(l),ajlmean(l)/sline(l),eta(l)/opa(l),opa(l)
! enddo
! endif
! if(xlambda.gt.256..and.xlambda.lt.257.) then
! do l=1,nd
! write(*,100) l,alo(l),ajlmean(l),ajlmean(l)/sline(l),eta(l)/opa(l),opa(l)
! enddo
! endif
! if(xlambda.gt.243..and.xlambda.lt.244.) then
! do l=1,nd
! write(*,100) l,alo(l),ajlmean(l),ajlmean(l)/sline(l),eta(l)/opa(l),opa(l)
! enddo
! endif
! if(xlambda.gt.237..and.xlambda.lt.238.) then
! do l=1,nd
! write(*,100) l,alo(l),ajlmean(l),ajlmean(l)/sline(l),eta(l)/opa(l),opa(l)
! enddo
! endif

! 100 format(i3,2x,f12.6,2x,e12.6,2x,f10.4,2x,e12.6,2x,e10.4)

  Return

End Subroutine

!-----------------------------------------------------------------------

Subroutine raycont(jp, z, r, nd, np, u, v, opa, eta, aic, dbdr, &
  corrfc, tauz1, tau1, pp)

  Use :: nlte_type
  Use :: nlte_dim
  Use :: nlte_var, Only: aak1 => ak1, aaz1 => az1, slinek

  Implicit None


! line radiation transfer in the comoving frame from a
! given source function for a given impact parameter jp.
! the integration is carried out in space and freq. to
! yield the mean line intensity

! new version: inlinen aller algebraischen routinen

! j-bar

! .. parameters ..
  Integer (i4b), Parameter :: nd1 = id_ndept, nf = id_nfcmf
! ..
! .. scalar arguments ..
  Real (dp) :: aic, corrfc, dbdr
  Integer (i4b) :: jp, nd, np
! ..
! .. array arguments ..
  Real (dp) :: opa(nd), eta(nd), pp(nd1), r(nd), tau1(nd), tauz1(nd), u(nd), &
    v(nd-1), z(nd)

! ..
! .. local scalars ..
  Real (dp) :: a, ak1, az1, b, daz, dazm, dbz, dbzm, dt, dtz, dtzm, dx, dxz, &
    dxz1, dxzm, dxzm1, h1, help, rmax, sc, tauc, tauz, tb1

  Integer (i4b) :: k, l, lmax, lz, ni
! ..
! .. local arrays ..
  Real (dp) :: ga(nd1), h(nd1), qq(nd1), s(nd1), ta(nd1), tb(nd1), tc(nd1), &
    ub(nd1), va(nd1), vb(nd1)
! ..
! .. intrinsic functions ..
! INTRINSIC EXP,MAX,MIN0
! ..

  lmax = min0(np+1-jp, nd)
  lz = lmax - 1
  rmax = r(1)

! until further evidence, use optically thick cont. without line background
! in case of multi-line transport

  tauc = opa(1)*rmax/3.D0
  sc = 0.D0

  If (tauc>=1.D0) sc = eta(1)/opa(1)

kloop: Do k = 2, nf

!   ***  outer boundary condition  -  2nd order

    ak1 = aak1(k, 1)
    az1 = aaz1(k, 1)
    dx = pp(1)*ak1

!   approximate treatment for optically thick continuum
!   (rho^2 processes)
!   respectively optically thick line at outer boundary

    s(1) = sc
    tc(1) = tauz1(1)*ak1
    tb(1) = tc(1) + dx + 1.D0
    ub(1) = dx
    vb(1) = 0.D0
    dtzm = tauz1(1)*az1

!   ***  for g and h, the matrix element s are not different from
!   ***  inner points

    dxzm = .5D0*(pp(1)+pp(2))*az1
    dxzm1 = 1.D0/(1.D0+dxzm)
    dazm = dtzm*dxzm1
    dbzm = dxzm*dxzm1
    ga(1) = -dazm
    h(1) = dbzm

!   ***  non-boundary points

    Do l = 2, lz
      ak1 = aak1(k, l)
      az1 = aaz1(k, l)
      s(l) = slinek(k, l)
      dt = tau1(l)*ak1
      dtz = tauz1(l)*az1
      dx = pp(l)*ak1
      dxz = .5D0*(pp(l)+pp(l+1))*az1
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

    If (lmax==nd) Then

!     ***  inner boundary condition (core rays)  -  only to first order
!     ***  diffusion-approximation

      ak1 = aak1(k, l)
      s(l) = aic + z(l)*ak1*dbdr*corrfc
      dx = pp(l)*ak1
      dt = tauz1(l-1)*ak1
      ta(l) = dt
      tb(l) = dt + dx + 1.D0
      ub(l) = dx
      va(l) = 0.0D0
      vb(l) = 0.0D0
    Else

!     ***  inner boundary condition (non-core rays)  -  second order

      ak1 = aak1(k, l)
      s(l) = slinek(k, l)
      tauz = z(lz)
      dt = ak1/tauz
      dx = pp(l)*ak1
      ta(l) = 2.D0*dt*dazm !        instead of daz
      tb(l) = ta(l) + dx + 1.D0
      ub(l) = dx
      va(l) = -2.D0*dt*dbzm !       instead of dbz
    End If

!   -----------------------------------------------------------------------
!   ***  continuation of subroutine ray


!   inlinen von vmalv

    b = v(1)
    qq(1) = vb(1)*b


    Do l = 2, lz
      a = b
      b = v(l)
      qq(l) = va(l)*a + vb(l)*b
    End Do

    qq(lmax) = va(lmax)*b

!   u = ub * u + qq + s

    Do l = 1, lmax
      u(l) = ub(l)*u(l) + qq(l) + s(l)
    End Do

!   ***     inlining of invtri

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
      u(ni) = u(ni) + tc(ni)*u(ni+1)
    End Do

    Do l = 1, lmax

!     i don't like the following statement at all, but as one says
!     in germany:  man hat schon pferde kotzen sehen.
!     actually, i found one case where due to a very strange behaviour
!     of occupation numbers v became strongly negative and in the
!     course also qq, such that u became negative. once u is negative,
!     this can no longer be compensated. hence:

      u(l) = max(u(l), 0.D0)
    End Do

!   now u is the new field at index k

!   inlinen von gmalu und
!   v = h * v + gmalu

    b = u(1)
    Do l = 1, lz
      a = b
      b = u(l+1)
      v(l) = h(l)*v(l) + ga(l)*(a-b)
    End Do

!   now v is the new field at index k

!   adding the new u to ajlmean and to the 0.and 2nd moments aj,ak


    u(lmax+1:nd) = 0.D0
    v(lz+1:nd-1) = 0.D0
!   print*,jp,k,U
!   print*,jp,k,v

  End Do kloop

  Return

End Subroutine

!-----------------------------------------------------------------------

Subroutine ray_inv(sing, ii, nov, jp, z, r, nd, np, u, v, pweight, &
  deltax, opa, eta, opal, aic, dbdr, corrfc, w0, ajlmean, tauz1, tau1, &
  pp, xlambda)

  Use :: nlte_type
  Use :: nlte_dim
  Use :: nlte_var, Only: aak1 => ak1, aaz1 => az1, slinek
  Use :: cmf_multi, Only: infotab, usave, vsave, iplus, iminus

  Implicit None


! line radiation transfer in the comoving frame from a
! given source function for a given impact parameter jp.
! the integration is carried out in space and freq. to
! yield the mean line intensity, all for inverted transitions!

! NOTE: here, we solve over z, which has been shown to
! be stable (if implicit scheme used) by
! Mihalas, Kunasz and Hummer, 1975, and has
! been checked carefully by R. Venero (March 2001)


! NOTE ALSO: no ALO calculated
! j-bar

! .. parameters ..
  Integer (i4b), Parameter :: nd1 = id_ndept, nf = id_nfcmf
! ..
! .. scalar arguments ..
  Real (dp) :: aic, corrfc, dbdr, deltax, xlambda
  Integer (i4b) :: ii, nov, jp, nd, np
! ..
! .. array arguments ..
  Real (dp) :: ajlmean(nd), eta(nd), opa(nd), opal(nd), pp(nd1), &
    pweight(nf, nd), r(nd), tau1(nd), tauz1(nd), u(nd), v(nd-1), &
    w0(nd), z(nd)

  Logical :: sing
! ..
! .. local scalars ..
  Real (dp) :: a, ak1, az1, b, daz, dazm, dbz, dbzm, dt, dtz, dtzm, dx, dxz, &
    dxz1, dxzm, dxzm1, h1, help, pw, rmax, sc, tauc, taus, tauz, tb1, w
  Integer (i4b) :: i, k, l, lmax, lz, ni
! ..
! .. local arrays ..
  Real (dp) :: ga(nd1), h(nd1), qq(nd1), s(nd1), ta(nd1), tb(nd1), tc(nd1), &
    ub(nd1), va(nd1), vb(nd1)
! ..
! .. intrinsic functions ..
! INTRINSIC MIN0
! ..


  lmax = min0(np+1-jp, nd)
  lz = lmax - 1
  rmax = r(1)

  sc = 0.D0
  tauc = opa(1)*rmax/3.D0
  If (tauc>=1.D0) sc = eta(1)

  taus = opal(1)/(pp(1)*deltax)
  If (taus<-0.5) Then
    Print *, ' WARNING!!!!: optically thick INVERTED line!!! LAM = ', xlambda, &
      ' JP = ', jp
    Write (999, *) ' WARNING!!!!: optically thick INVERTED line!!! LAM = ', &
      xlambda, ' JP = ', jp
  End If

! loop for all original frequency points

kloop: Do k = 1, nf

    If (k==1) Go To 100


!   ***  inlined subroutine cmfset for inverted levels


!   ***  outer boundary condition  -  2nd order

    ak1 = aak1(k, 1)
    az1 = aaz1(k, 1)
    dx = pp(1)
    s(1) = sc
    tc(1) = tauz1(1)
    tb(1) = tc(1) + dx + ak1
    ub(1) = dx
    vb(1) = 0.D0
    dtzm = tauz1(1)

!   ***  for g and h, the matrix element s are not different from
!   ***  inner points

    dxzm = .5D0*(pp(1)+pp(2))
    dxzm1 = 1.D0/(az1+dxzm)
    dazm = dtzm*dxzm1
    dbzm = dxzm*dxzm1
    ga(1) = -dazm
    h(1) = dbzm

!   ***  non-boundary points

    Do l = 2, lz
      ak1 = aak1(k, l)
      az1 = aaz1(k, l)
!     r for inverted levels S(L) is just the emisivity
      s(l) = slinek(k, l)
      dt = tau1(l)
      dtz = tauz1(l)
      dx = pp(l)
      dxz = .5D0*(pp(l)+pp(l+1))
      dxz1 = 1.D0/(az1+dxz)
      daz = dtz*dxz1
      dbz = dxz*dxz1
      ta(l) = dt*dazm
      tc(l) = dt*daz
      tb(l) = ta(l) + tc(l) + dx + ak1
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

!     ***  inner boundary condition (core rays)  -  only to first order
!     ***  diffusion-approximation

      ak1 = aak1(k, l)
      s(l) = aic*ak1 + z(l)*dbdr*corrfc
      dx = pp(l)
      dt = tauz1(l-1)
      ta(l) = dt
      tb(l) = dt + dx + ak1
      ub(l) = dx
      va(l) = 0.0D0
      vb(l) = 0.0D0
    Else

!     ***  inner boundary condition (non-core rays)  -  second order

      ak1 = aak1(k, l)
!     r   same as above
      s(l) = slinek(k, l)
      tauz = z(lz)
      dt = 1.D0/tauz
      dx = pp(l)
      ta(l) = 2.D0*dt*dazm
      tb(l) = ta(l) + dx + ak1
      ub(l) = dx
      va(l) = -2.D0*dt*dbzm
    End If

!   -----------------------------------------------------------------------
!   ***  continuation of subroutine ray


!   inlinen von vmalv

    b = v(1)
    qq(1) = vb(1)*b


    Do l = 2, lz
      a = b
      b = v(l)
      qq(l) = va(l)*a + vb(l)*b
    End Do

    qq(lmax) = va(lmax)*b

!   u = ub * u + qq + s

    Do l = 1, lmax
      u(l) = ub(l)*u(l) + qq(l) + s(l)
    End Do

!   ***     inlining of invtri

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
      u(ni) = u(ni) + tc(ni)*u(ni+1)
    End Do

    Do l = 1, lmax
!     see comment in subroutine Ray
      u(l) = max(u(l), 0.D0)
    End Do


    b = u(1)
    Do l = 1, lz
      a = b
      b = u(l+1)
      v(l) = h(l)*v(l) + ga(l)*(a-b)
    End Do

!   now v is the new field at index k

!   adding the new u to ajlmean and to the 0.and 2nd moments aj,ak

    u(lmax+1:nd) = 0.D0
    v(lz+1:nd-1) = 0.D0

100 Continue

    Do l = 1, lmax
      pw = pweight(k, l)*w0(l)
      ajlmean(l) = ajlmean(l) + u(l)*pw
    End Do

    If (sing) Cycle kloop

!   store u,v if necessary

    Do i = 1, nov

      If (infotab(i)%indbw==ii) Then
        If (int(infotab(i)%if)==nf .And. k==nf) Then
          usave(1:lmax, jp, i) = u(1:lmax)
          vsave(1:lz, jp, i) = v(1:lz)
          usave(nd, np-1, i) = nf
          iminus(jp, i) = s(1)
          If (lmax==nd) iplus(jp, i) = s(nd)/aak1(nf, nd)
!         IF(JP.EQ.1) PRINT*,II,'NF STORED',INFOTAB(I)%IF,' FOR LINE ',I,'
!         ',NF
        Else If (int(infotab(i)%if)==k) Then
          w = k + 1.D0 - infotab(i)%if
          usave(1:lmax, jp, i) = w*u(1:lmax)
          vsave(1:lz, jp, i) = w*v(1:lz)
          usave(nd, np-1, i) = w
!         IF(JP.EQ.1) PRINT*,II,K,' STORED (FIRST)',INFOTAB(I)%IF,' FOR LINE
!         ',I,' ',W
        Else If (int(infotab(i)%if)==k-1) Then
          w = infotab(i)%if + 1.D0 - k
          usave(1:lmax, jp, i) = usave(1:lmax, jp, i) + w*u(1:lmax)
          vsave(1:lz, jp, i) = vsave(1:lz, jp, i) + w*v(1:lz)
          vsave(nd-1, np-1, i) = w
!         IF(JP.EQ.1) PRINT*,II,K,' STORED (2ND)',INFOTAB(I)%IF,' FOR LINE
!         ',I,' ',W
        End If
      End If

    End Do

  End Do kloop

  Return

End Subroutine

!-----------------------------------------------------------------------

Subroutine cmfprep(ll, velr, vgrad, nfir, nlas, blevel, xnel, sr, vmax, nato, &
  xmust, clf, nfirm, nlasm, nfi, nla, startcmf, optmixed, ipres, r, &
  velo, fic, tcl_fac_line, tcl_fac_cont)

! UPDATED FOR CLUMPING

  Use :: nlte_type
  Use :: nlte_dim
  Use :: fund_const

  Use :: princesa_var, Only: data, gl, index1, indat1, labl

  Use :: nlte_var, Only: mmlow, mmup, xlamcmf, vdopcmf, blucmf, bulcmf, &
    aulcmf, lablinelo, lablineup, opaclin, tlumat, tulmat, betold, indexrbb, &
    indexcmf, indexsa, modnam, vther, fre, lte_update

  Implicit None

! --- calculates rbb-rates from old occupation numbers for sob. treatment
! --- prepares cmf calculations and decides which lines are treated in
! --- cmf. (optmixed = 0.: all lines; optmixed = 1., lines according to
! --- the criteria specified below.)


! .. parameters ..
  Integer (i4b), Parameter :: nnrec = id_llevs, nd = id_ndept
! ..
! .. scalar arguments ..
  Real (dp) :: optmixed, sr, velr, vgrad, vmax, xmust, xnel, clf
  Real (dp) :: fic, tcl_fac_line, tcl_fac_cont
  Integer (i4b) :: ipres, ll, nato, nfi, nfir, nfirm, nla, nlas, nlasm
  Logical :: startcmf
! ..
! .. array arguments ..
  Real (dp) :: blevel(nnrec), r(id_ndept), velo(id_ndept)
! ..
! .. local scalars ..
  Real (dp) :: blu, bul, er, flu, tauc, tlu, tul, vcont, vsobo, vthmax, xlam, &
    xlamc, xnl, xnu, maxlam
  Integer (i4b) :: i, icount, ii, indcmf, indi, indmat, l, m, ml, mu
  Logical :: start
! ..
! .. local arrays ..
  Integer (i4b) :: indexold(id_nttrd)
! ..
! .. external subroutines ..
  External :: tratcmf, tratrad
! ..
! .. intrinsic functions ..
! INTRINSIC INT,MAX
! ..

! ---  for cmf-lines, we have the following
! ---- correspondance (before having performed the cmf treatment)
! tlumat -- opal
! tulmat -- sline

! .. data statements ..
  Data start/.True./
! ..

! AFTER LTE-UPDATE, CREATE NEW CMF-LIST (CHANGE IN IMIN,IMAX)
  If (.Not. start .And. lte_update .And. nato==1 .And. ll==1) start = .True.

  If (start) Then
!   THIS BLOCK MUST BE CREATED ONLY ONCE AT NATO=1 AND LL=1, &
!   SINCE ALL LINES (FROM ALL ELEMENTS) ARE INITIALIZED IN THE I-LOOP'
    Open (1, Err=100, File=trim(modnam)//'/INDEXCMF', Status='OLD')
    Rewind 1
    Read (1, Fmt=*)(indexold(ii), ii=1, id_nttrd)
    Close (1)
    ipres = 1
    Go To 110

100 Continue

    ipres = 2

110 Continue
    icount = 0
    Print *, ipres
    Do i = 1, index1
      indi = indat1(i)
      er = data(indi)
      m = int(er)

      If (m==1) Then
        icount = icount + 1
        indexrbb(i) = icount

!       ---  at first, all lines are assumed as sobolev ones

        indexcmf(icount) = 1
      Else If (m==2) Then
        indexrbb(i) = 0
      Else
        Write (999, *) ' STOP: ERROR IN RBB OR CBB'
        Stop ' ERROR IN RBB OR CBB'
      End If

    End Do

    If (icount/=id_nttrd) Then
      Write (999, *) ' STOP: NOT ALL RBB TRANSITIONS FOUND'
      Stop ' NOT ALL RBB TRANSITIONS FOUND'
    End If
    start = .False.
  End If

! ---  now, according to the specified criteria, lines are changed from
! ---  sobolev to cmf treatment
! ---  only at first depth point, sucessively for all elements

  If (startcmf) Then

!   ***  this is the border velocity we require for the mixed treatment.

    vthmax = 0.2D0*vther(nato)/vmax
    maxlam = 1.D8/fre(1)

iloop: Do i = 1, index1
      ml = mmlow(i)
      mu = mmup(i)

!     ---  note: since all ions are to be inspected,
!     we have             nfi   and nla   (depth indep., maximum set)
!     different from      nfir  and nlas  (actual, depth dependent)
!     and                 nfirm and nlasm (depth indep., minimum set)

      If (ml<nfi .Or. ml>=nla) Cycle
      If (mu>nla) Cycle
      indmat = indexrbb(i)
      If (indmat==0) Cycle

      If (indexcmf(indmat)/=1) Then
        Print *, startcmf, start, lte_update, ipres
        Print *, i, ml, mu, indmat, indexcmf(indmat)
        Print *, nfi, nla
        Write (999, *) ' STOP: ERROR IN INDCMF-PHILOSOPHY'
        Stop ' ERROR IN INDCMF-PHILOSOPHY'
      End If
      indi = indat1(i)
      xlam = data(indi+1)
      xlamcmf(indmat) = xlam
      lablinelo(indmat) = labl(ml)
      lablineup(indmat) = labl(mu)

!     ***  this is the present criterium to prespecify the maximum number
!     of cmf lines: lambda < 100000 (150000 to included H5->6 transition)
!     and to account only for ions which are everywhere present

      If (xlam<=maxlam) Then
!       in this way, one can simulate the optmixed approach
!       IF (XLAM.LE.MAXLAM .and. nato.le.2) THEN
        If (ml>=nfirm .And. mu<=nlasm) Then
          indexcmf(indmat) = 2
        End If
      End If

!     ***  now, we specify the cmf lines according to the previous sa-results
!     in case of optmixed

      If (optmixed==1. .And. indexcmf(indmat)==2) Then
!       -----
!       -----note: if this path needs to be continued, take care to update
!       corresponding quantities if model (phot.struct.) has been updated
        Write (999, *) ' STOP: OPTMIXED NOT IMPLEMENTED YET'
        Stop ' OPTMIXED NOT IMPLEMENTED YET'
        If (ipres==2) Then
          vsobo = 0.D0
          vcont = vsobo
          tauc = 0.D0

          Do l = 1, nd

!           ***  if the lines become optically thick (corresponding to a mean
!           tau_s
!           ***  of unity) below the border velocity, they can be treated in
!           sa.
!           ***  they can be treated also in sa, if the continuum becomes
!           optically
!           ***  thick prior to the considered line.

            If (betold(indmat,l)<0.632D0) vsobo = max(vsobo, velo(l))
            If (l/=1) tauc = tauc + .5D0*(opaclin(indmat,l-1)+opaclin(indmat,l &
              ))*(r(l-1)-r(l))*sr
            If (tauc>=1.D0) vcont = max(vcont, velo(l))
          End Do

          If (vsobo<vthmax .Or. vsobo<vcont) indexcmf(indmat) = 1

        Else If (ipres==1) Then
          indexcmf(indmat) = indexold(indmat)
        Else
          Write (999, *) ' STOP: ERROR IN IPRES (CMFPREP)'
          Stop ' ERROR IN IPRES (CMFPREP)'
        End If

      End If

    End Do iloop

  End If

! **** general path

! ---  setup of rbb rates with old occ.no.

iiloop: Do ii = 1, index1
    indi = indat1(ii)
    er = data(indi)
    m = int(er)
    If (m/=1) Cycle
    ml = mmlow(ii)
    mu = mmup(ii)

    If (ml<nfir .Or. ml>=nlas) Cycle
    If (mu>nlas) Cycle

    xnl = blevel(ml)
    xnu = blevel(mu)
!   ALL OCCUPATION NUMBERS CONSIDERED HERE SHOULD BE DEFINED
    If (xnl*xnu==0.) Then
      Write (999, *) ' STOP: OCC. NUM NOT DEFINED IN CMFPREP'
      Stop ' OCC. NUM NOT DEFINED IN CMFPREP'
    End If

    xlamc = data(indi+1)*1.D-8
    flu = data(indi+2)

    indmat = indexrbb(ii)
    indcmf = indexcmf(indmat)

!   ---     final check if everything was ok

    If (xlamcmf(indmat)==0.D0) Then
      Write (999, *) ' STOP: ERROR IN INDEX GYMNASTICS'
      Stop ' ERROR IN INDEX GYMNASTICS'
    End If

!   ---  sobolev treatment

    If (indcmf==1) Then
      betold(indmat, ll) = 0.D0
      If (ll==1) vdopcmf(indmat) = 0.D0
      Call tratrad(ll, ml, mu, flu, xlamc, xnl, xnu, xnel, gl, velr, vgrad, &
        tlu, tul, betold(indmat,ll), sr, vmax, xmust, ii, clf, &
        vther(nato), fic, tcl_fac_line)

      tlumat(indmat, ll) = tlu
      tulmat(indmat, ll) = tul

      indexsa(indmat) = indexsa(indmat) + 1

!     ---  cmf preparation

    Else If (indcmf==2) Then
      Call tratcmf(ll, ml, mu, flu, xlamc, xnl, xnu, xnel, clf, gl, sr, vmax, &
        ii, fic, tcl_fac_line)

      If (ll==1) Then
        vdopcmf(indmat) = vther(nato)
        blu = e2mc2/hh*xlamc*4.D0*pi**2*flu
        bul = gl(ml)*blu/gl(mu)
        blucmf(indmat) = blu
        bulcmf(indmat) = bul
        aulcmf(indmat) = hc2/xlamc**3*bul
      End If

    Else
      Write (999, *) ' STOP: ERROR IN INDEXCMF(CMFPREP)'
      Stop ' ERROR IN INDEXCMF(CMFPREP)'
    End If

  End Do iiloop

  Return

End Subroutine

!-----------------------------------------------------------------------

Subroutine tratcmf(l, nlo, nu, flu, xlamc, xnl, xnu, xne, clf, gl, sr, vmax, &
     ii, fic, tcl_fac_line)

! UPDATED FOR CLUMPING: same assumption as in TRATRAD, i.e., thin clumps small
! optically thick clumps as in TRATRAD
! JO Sept 2023: updated w.r.t. position of continuum quantities

! after many tests, old implemenation preferred. To switch to newer one (not
! preferred)
! (see standard_v10.6.5_old), disable "PATH=1"  and "GOTO 100" statement

  Use :: nlte_type
  Use :: nlte_dim
  Use :: fund_const
  Use :: nlte_var, Only: fre, ifre, opac, strue, xj, tlumat, tulmat, opacm, &
    scontm, indexrbb, optlines, thomson_lines, optcmf_all, kcmf_start, &
    kcmf_end, restart_cmfall
! ,RAYLEIGH
  Use :: nlte_opt, Only: opt_hei_sing

  Use :: nlte_porvor, Only: opa_eff_rat_old

  Use :: nlte_wavcon, Only: wwavblue, wwavcon


  Implicit None

! ---  calculation of line/continuum quantities for cmf-treatment

! .. parameters ..
! Integer (i4b), Parameter :: nd1 = id_ndept
! Integer (i4b), Parameter :: ifretot = id_frec1
! ..
! .. scalar arguments ..
  Real (dp) :: flu, sr, vmax, xlamc, xne, xnl, xnu, clf
  Real (dp) :: fic, tcl_fac_line
  Integer (i4b) :: ii, l, nlo, nu
! ..
! .. array arguments ..
  Real (dp) :: gl(id_llevs)
! ..
! .. local scalars ..
  Real (dp) :: deln, opacon, s1, s2, scont, sline, wilog, xkli, xkline, xnue, &
    xnugl, xxx
  Real (dp) :: tcl, fac1
  Integer (i4b) :: indl, k, nlambox, nlambox1, nlambox2, path, kk
! ..
! .. intrinsic functions ..
! INTRINSIC LOG10
! ..
! ..

  indl = indexrbb(ii)

! ---- line opacity
! ---- (warning! xlamc= central lambda in cm)

  xnue = 1.D0/xlamc

  If (xnue>fre(ifre) .Or. xnue<fre(1)) Then
    Write (999, *) ' STOP: LINE OUT OF FREQ. RANGE'
    Stop ' LINE OUT OF FREQ. RANGE'
  End If

  xnugl = xnu/gl(nu)
  deln = (xnl/gl(nlo)-xnugl)

  xkli = pi*e2mc2*clight*gl(nlo)*flu*deln/clf ! corrected
  xkline = xkli*xlamc*sr/vmax

  sline = 3.973D-16*xnue**3*xnugl/deln ! remains
! TLUMAT(INDL,L) = XKLINE
! TULMAT(INDL,L) = SLINE

  k = ifre
  xxx = fre(ifre)
  Do While (xnue<xxx)
    k = k - 1
    xxx = fre(k)
  End Do

  path = 1
  Go To 100
! to switch to alternative description, comment out previous two statements

  path = 0
  If (xlamc*1.D8>=910. .Or. xlamc*1.D8<=wwavblue .Or. .Not. optlines) Then
    path = 1 !                      interpolation as before (for XLAMC > 910 A
!   OR XLAMC < WWAVBLUE
!   OR .NOT.OPTLINES)
  Else
    nlambox = 1 + log(xlamc*1.D8/wwavblue)/wwavcon
    nlambox1 = 1 + log(1.D8/(fre(k)*wwavblue))/wwavcon
    nlambox2 = 1 + log(1.D8/(fre(k+1)*wwavblue))/wwavcon

    If (nlambox==nlambox1 .And. nlambox==nlambox2) Then
!     line in between K,K+1 of same NLAMBOX -> interpolation as before
      path = 1
    Else If (nlambox==nlambox1) Then
!     line in same NLAMBOX interval as FRE(K) -> no interpolation
      path = 2
    Else If (nlambox==nlambox2) Then
!     line in same NLAMBOX interval as FRE(K+1) -> no interpolation
      path = 3
    Else
      Print *, wwavblue, wwavcon
      Print *, xlamc, k, fre(k), k + 1, fre(k+1)
      Print *, nlambox, nlambox1, nlambox2
      Write (999, *) &
        ' STOP: NLAMBOX NE NLAMBOX1/NLAMBOX2 for wavblue < lambda < 910 A'
      Stop ' NLAMBOX NE NLAMBOX1/NLAMBOX2 for wavblue < lambda < 910 A'
    End If
  End If

100 Continue
  If (path==0) Then
    Write (999, *) ' STOP: SOMETHING ROTTEN WITH PATH (subr. TRATCMF)'
    Stop ' SOMETHING ROTTEN WITH PATH (subr. TRATCMF)'
  End If


! ---- effective quantity calculated as in TRATRAD (incl. inversion)
  tcl = abs(xkline)*tcl_fac_line*vmax/sr ! SR/VMAX term included in
! tcl_fac_line
  fac1 = (1.+fic*tcl)/(1.+tcl)
! IF (FAC1.LT.0.0) FAC1=1.0D0   ! old treatment of inversion
  fac1 = min(fac1, opa_eff_rat_old(l,k))
  If (fac1>1.0D0 .Or. fac1<=0.0D0) Then
    Print *, fac1, opa_eff_rat_old(l, k)
    Print *, l, k, tcl
    Write (999, *) ' STOP: OPA_EFF > <OPA>, TRATCMF_1'
    Stop ' OPA_EFF > <OPA>, TRATCMF_1'
  End If
! if(xlamc*1.d8.gt.237. .and. xlamc*1.d8.lt.238.) fac1=1.
! if(xlamc*1.d8.gt.243. .and. xlamc*1.d8.lt.244.) fac1=1.
! if(xlamc*1.d8.gt.256. .and. xlamc*1.d8.lt.257.) fac1=1.

! if(xlamc*1.d8.gt.237. .and. xlamc*1.d8.lt.238.) &
! print*,'237',l,tcl,fac1,opa_eff_rat_old(l,k)
! if(xlamc*1.d8.gt.243. .and. xlamc*1.d8.lt.244.) &
! print*,'243',l,tcl,fac1,opa_eff_rat_old(l,k)
! if(xlamc*1.d8.gt.256. .and. xlamc*1.d8.lt.257.) &
! print*,'256',l,tcl,fac1,opa_eff_rat_old(l,k)
! if(xlamc*1.d8.gt.303. .and. xlamc*1.d8.lt.304.) &
! print*,'303',l,tcl,fac1,opa_eff_rat_old(l,k)

! IF (XLAMC*1.d8.gt.1000.and.XLAMC*1.d8.lt.3000.) THEN
! Print*,'test-JS-CMF:'
! Print*,L,K,XLAMC*1.d8
! Print*,XKLINE,TCL_FAC_LINE
! Print*,FAC1,OPA_EFF_RAT_OLD(L,K),L
! Print*,'-----------'
! ENDIF
  xkline = xkline*fac1
! ---- now effective

  tlumat(indl, l) = xkline
  tulmat(indl, l) = sline

  If (path==1) Then

!   ---- interpolation weights for frecuency magnitudes. interpolation
!   ---- carries out magnitude versus log(nue) (linear interp.)

    wilog = log10(fre(k)/xnue)/log10(fre(k+1)/fre(k))

!   ---- s_cont and opac*sr

    If (optlines) Then

      If (optcmf_all) Then
        If (restart_cmfall) Then
          If (xj(1,ifre+1)/=1) Then
            Write (999, *) &
              ' STOP: RESTART_CMFALL, BUT XJ (CMF) INDICATED IN TRATCMF'
            Stop ' RESTART_CMFALL, BUT XJ (CMF) INDICATED IN TRATCMF'
          End If
!         at restart, all frequencies possible
        Else
          If (xj(1,ifre+1)/=2) Then
            Write (999, *) ' STOP: OPTCMF_ALL = TRUE, BUT APPROX. &
              &XJ (OBSFRAM) USED IN TRATRAD'
            Stop &
              ' OPTCMF_ALL = TRUE, BUT APPROX. XJ (OBSFRAM) USED IN TRATRAD'
          End If
!         only evaluated at non-cmf_all frequencies
          If (xnue<fre(kcmf_end) .And. xnue>fre(kcmf_start)) Then
            If (opt_hei_sing .And. abs(xlamc*1.D8-584.)<1.) Then
              Continue
            Else
              Print *, 'single line treatment at', xlamc*1.D8
!             JO comment out when specific lines are treated in single line
!             approximation
              Write (999, *) ' STOP: OPTCMF_ALL = TRUE, BUT SINGLE &
                &CMF-LINE INSIDE CMF-RANGE (TRATCMF)'
              Stop ' OPTCMF_ALL = TRUE, BUT SINGLE CMF-LINE INSIDE &
                &CMF-RANGE (TRATCMF)'
            End If
          End If
        End If
      End If

!     S1 = STRUE(L,K) + XJ(L,K)* &
!     &
!     (XNE*SIGMAE*AMH/CLF+THOMSON_LINES(L,K)+RAYLEIGH(L,K)/CLF)/OPAC(L,K)
!     S2 = STRUE(L,K+1) + XJ(L,K+1)* &
!     &
!     (XNE*SIGMAE*AMH/CLF+THOMSON_LINES(L,K+1)+RAYLEIGH(L,K+1)/CLF)/OPAC(L,K+1)

!     ---- same as in TRATRAD
      s1 = strue(l, k) + xj(l, k)*(xne*sigmae*amh/clf+thomson_lines(l,k))/opac &
        (l, k)*opa_eff_rat_old(l, k) ! back-corrected
      s2 = strue(l, k+1) + xj(l, k+1)*(xne*sigmae*amh/clf+thomson_lines(l,k+1) &
        )/opac(l, k+1)*opa_eff_rat_old(l, k+1) ! back-corrected
    Else
!     ---- as above
      s1 = strue(l, k) + xne*sigmae*amh*xj(l, k)/clf/opac(l, k)* &
        opa_eff_rat_old(l, k) !     back-corrected
      s2 = strue(l, k+1) + xne*sigmae*amh*xj(l, k+1)/clf/opac(l, k+1)* &
        opa_eff_rat_old(l, k+1) !   back-corrected
    End If

    scont = s1*10.D0**(log10(s1/s2)*wilog)

    opacon = opac(l, k)*10.D0**(log10(opac(l,k)/opac(l,k+1))*wilog)

    opacm(indl, l) = opacon*sr !    already OK
!   ---- OPACM now effective quantity
    scontm(indl, l) = scont

  Else
!   JO Sept 2023: NOTE THAT FOLLOWING PATH=2 AND PATH=3 BRANCHES NOT TESTED
!   (SINCE PATH=1 PER DEFAULT). IF ENABLED, PRESUMABLY NEED TO BE UPDATED
!   SIMILAR TO PATH=1 (IN CASE OF OPTCMF_ALL)

    If (.Not. optlines) Then
      Write (999, *) ' STOP: OPTLINES = .FALSE. FOR PATH > 1 (subr. TRATCMF)'
      Stop ' OPTLINES = .FALSE. FOR PATH > 1 (subr. TRATCMF)'
    End If
!   JO Sept 2023. NO INTERPOLATION, BUT VALUE CORRESPONDING TO K OR K+1 WITH
!   SAME NLAMBOX AS TRANSITION.
!   NOTE: EVEN FOR WEAK LINE-BG., THIS SHOULD WORK, SINCE CONTINUUM VARIATIONS
!   NOT LARGE OVER ONE LINE-INTERVAL (AND, IF EDGE, THE CORRESPONDING
!   INDEX SHOULD BE CORRECT AS WELL).
    If (path==2) Then
      kk = k
    Else If (path==3) Then
      kk = k + 1
    Else
      Write (999, *) ' STOP: WRONG PHILOSOPHY IN subr. CMFPREP'
      Stop ' WRONG PHILOSOPHY IN subr. CMFPREP'
    End If

    s1 = strue(l, kk) + xj(l, kk)*(xne*sigmae*amh/clf+thomson_lines(l,kk))/ &
      opac(l, kk)*opa_eff_rat_old(l, kk) ! back-corrected


    scont = s1

    opacon = opac(l, kk)

    opacm(indl, l) = opacon*sr !    already OK
!   ---- OPACM now effective quantity
    scontm(indl, l) = scont
  End If

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
    If (jstack>nstack) Then
      Write (999, *) ' STOP: NSTACK TOO SMALL IN INDEXX'
      Stop ' NSTACK TOO SMALL IN INDEXX'
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

Subroutine sobprep(ll, velr, vgrad, nfir, nlas, blevel, xnel, sr, vmax, nato, &
  xmust, clf, fic, tcl_fac_line)

! UPDATED FOR CLUMPING

  Use :: nlte_type
  Use :: nlte_dim

  Use :: princesa_var, Only: data, gl, index1, indat1

  Use :: nlte_var, Only: mmlow, mmup, tlumat, tulmat, betold, indexrbb, vther

  Implicit None

! ---  calculates rbb-rates from old occupation numbers


! .. parameters ..
  Integer (i4b), Parameter :: nnrec = id_llevs
! ..
! .. scalar arguments ..
  Real (dp) :: sr, velr, vgrad, vmax, xmust, xnel, clf
! THICK
  Real (dp) :: fic, tcl_fac_line
  Integer (i4b) :: ll, nato, nfir, nlas
! ..
! .. array arguments ..
  Real (dp) :: blevel(nnrec)
! ..
! .. local scalars ..
  Real (dp) :: er, flu, tlu, tul, xlamc, xnl, xnu
  Integer (i4b) :: i, icount, ii, indi, indmat, m, ml, mu

  Logical :: start
! ..
! .. external subroutines ..
  External :: tratrad
! ..
! .. intrinsic functions ..
! INTRINSIC INT
! ..
! .. data statements ..

  Data start/.True./
! ..

  If (start) Then
    icount = 0

    Do i = 1, index1
      indi = indat1(i)
      er = data(indi)
      m = int(er)

      If (m==1) Then
        icount = icount + 1
        indexrbb(i) = icount
      Else If (m==2) Then
        indexrbb(i) = 0
      Else
        Write (999, *) ' STOP: ERROR IN RBB OR CBB'
        Stop ' ERROR IN RBB OR CBB'
      End If
    End Do

    If (icount/=id_nttrd) Then
      Write (999, *) ' STOP: NOT ALL RBB TRANSITIONS FOUND'
      Stop ' NOT ALL RBB TRANSITIONS FOUND'
    End If
    start = .False.
  End If

! ---  setup of rbb rates with old occ.no.

  Do ii = 1, index1
    indi = indat1(ii)
    er = data(indi)
    m = int(er)
    If (m/=1) Cycle
    ml = mmlow(ii)
    mu = mmup(ii)
    If (ml<nfir .Or. ml>=nlas) Cycle
    If (mu>nlas) Cycle

    xnl = blevel(ml)
    xnu = blevel(mu)
    If (xnl*xnu==0.) Then
      Write (999, *) ' STOP: OCC. NUM NOT DEFINED IN SOBPREP'
      Stop ' OCC. NUM NOT DEFINED IN SOBPREP'
    End If

    xlamc = data(indi+1)*1.D-8
    flu = data(indi+2)
    indmat = indexrbb(ii)
    betold(indmat, ll) = 0.D0

    Call tratrad(ll, ml, mu, flu, xlamc, xnl, xnu, xnel, gl, velr, vgrad, tlu, &
      tul, betold(indmat,ll), sr, vmax, xmust, ii, clf, vther(nato), &
      fic, tcl_fac_line)

    tlumat(indmat, ll) = tlu
    tulmat(indmat, ll) = tul

  End Do

  Return

End Subroutine

!***********************************************************************

!subroutines: atomic quantities

!***********************************************************************

Subroutine rayleigh_scat(ll)

! Set up the hydrogen/helium I cross-sections according to fits
! lifted from the ATLAS program (Kurucz, 1970, SAO Special Report No.309).
! NOTE however that here the cutoff for hydrogen is set to
! 2.40e15 Hz = 1250 A, whereas Kurucz goes down to 1026 A (where
! the whole approximation is invalid because that's inmidst the Ly lines)

! This subroutine is presently not used, test calculations (not rigorous)
! for Vega have indicated that the influence is only weak.
! Note that for a final use the high frequency parts have to be modified!

  Use :: nlte_type
  Use :: nlte_dim
  Use :: princesa_var, Only: labl
  Use :: nlte_var, Only: fre, ifre
! ,RAYLEIGH  uncomment if subroutine used

  Implicit None

  Integer (i4b), Intent (In) :: ll

  Integer (i4b), Save :: levh = 0, levhe = 0
  Integer (i4b) :: i

  Real (dp) :: lam, lam2

  Logical :: start

  Data start/.True./

  If (start) Then
    start = .False.
    Do i = 1, id_llevs
      If (labl(i)=='H11') Then
        levh = i
        Exit
      End If
    End Do

    Do i = 1, id_llevs
      If (labl(i)=='HE11S1') Then
        levhe = i
        Exit
      End If
    End Do
  End If

  If (levh==0) Then
    Write (999, *) ' STOP: NO HYDROGEN IN ATMOSPHERE!'
    Stop ' NO HYDROGEN IN ATMOSPHERE!'
  End If

  Do i = 1, ifre

!   LAMBDA IN A
    lam = 1.D8/fre(i)
    lam = max(lam, 1250.D0)
    lam2 = lam*lam
!   RAYLEIGH(LL,I)= &    UNCOMMENT!!!!
!   BLEVEL(LEVH)*(5.799D-13/LAM2**2+1.422D-6/LAM2**3+2.784/LAM2**4)
  End Do

  If (levhe==0) Return

  Do i = 1, ifre

    lam = 1.D8/fre(i)
    lam = max(lam, 582.D0)
    lam2 = lam*lam
!   RAYLEIGH(LL,I)=RAYLEIGH(LL,I)+ BLEVEL(LEVHE)*5.484D-14/LAM2**2 * &
!   &   ((5.94D10/(LAM2-2.90D5)+2.44D5)/LAM2+1.D0)**2   UNCOMMENT!!!
  End Do

  Return
End Subroutine

!***********************************************************************

Subroutine hminus(xne, temp, lam, ehnuekt, ll, opahminus)

! calculates bf and ff opacities for H-minus

! lambda in A

! formulas and fits according to the references in:
! D.G. Gray, "The observation and analysis of stellar photospheres",
! 2nd Edition (1992), Cambridge Astrophysics Series 20, Cambridge
! University Press, p. 135-137

  Use :: nlte_type
  Use :: nlte_dim
  Use :: fund_const, Only: akb
  Use :: princesa_var, Only: nat, labat
  Use :: nlte_var, Only: enionnd

  Implicit None

  Real (dp), Intent (In) :: xne, temp, lam, ehnuekt
  Real (dp), Intent (Out) :: opahminus
  Integer (i4b), Intent (In) :: ll


  Real (dp) :: theta, thetalog, nh, pe, opabf, opaff, stim, xlam, f0, f1, f2

  Integer (i4b) :: k

! fit coefficients for bf opacity
  Real (dp) :: a0 = 1.99654D0, a1 = -1.18267D-5, a2 = 2.64243D-6, &
    a3 = -4.40524D-10, a4 = 3.23992D-14, a5 = -1.39568D-18, a6 = 2.78701D-23

! fit coefficients for ff opacity
  Real (dp) :: b0 = -2.2763D0, b1 = -1.685D0, b2 = .76661D0, b3 = -.0533464D0
  Real (dp) :: c0 = 15.2827D0, c1 = -9.2846D0, c2 = 1.99381D0, c3 = -.142631D0
  Real (dp) :: d0 = -197.789D0, d1 = 190.266D0, d2 = -67.9775D0, &
    d3 = 10.6913D0, d4 = -0.625151D0

  If (lam<2000.) Then
    opahminus = 0.D0
    Return
  End If

  Do k = 1, nat
    If (labat(k)=='H') Go To 100
  End Do
  Write (999, *) ' STOP: HYDROGEN NOT FOUND IN HMINUS'
  Stop ' HYDROGEN NOT FOUND IN HMINUS'

100 nh = enionnd(k, 1, ll) !        neutral hydrogen density
  pe = xne*akb*temp !               electron pressure

  theta = 5040.D0/temp

! bf opacity

  opabf = 0.
  If (lam<=16000.) Then
    opabf = a0 + lam*(a1+lam*(a2+lam*(a3+lam*(a4+lam*(a5+lam*a6)))))
    opabf = opabf*1.D-18
    stim = 1. - ehnuekt
    opabf = 4.158D-10*opabf*pe*theta**2.5D0*10.D0**(0.754D0*theta)*stim*nh
    If (opabf<0.) opabf = 0.
  End If

  thetalog = log10(theta)
  xlam = log10(lam)

  f0 = b0 + xlam*(b1+xlam*(b2+xlam*b3))
  f1 = c0 + xlam*(c1+xlam*(c2+xlam*c3))
  f2 = d0 + xlam*(d1+xlam*(d2+xlam*(d3+xlam*d4)))

  opaff = 1.D-26*pe*10.D0**(f0+f1*thetalog+f2*thetalog**2)*nh

  opahminus = opabf + opaff
  Return

End Subroutine

!-----------------------------------------------------------------------


Function agau(nivel, xnu, key, numi)

  Use :: nlte_type
  Use :: nlte_dim
  Implicit None

! this function returns the bound-free gaunt factor for h and he
! with the approximation used by herrero (tesis doctoral).

! .. scalar arguments ..
  Real (dp) :: xnu, agau
  Integer (i4b) :: nivel, numi
  Character :: key*6
! ..
! .. local arrays ..

  Real (dp) :: a(3), b(3), c(3)
! ..
! .. data statements ..
  Data a/0.9916D0, 1.105D0, 1.101D0/
  Data b/2.71852D13, -2.37496D14, -9.863186D13/
  Data c/ -2.26846D30, 4.07677D28, 1.03537D28/
! ..

  If (key=='H') Then

    If (numi/=1) Then
      Write (999, *) ' STOP: ERROR IN NUMI -H- GBF'
      Stop ' ERROR IN NUMI -H- GBF'
    End If

    If (nivel<=3) Then
      agau = a(nivel) + b(nivel)/xnu + c(nivel)/xnu/xnu
      Return
    Else
      agau = 1.D0
      Return
    End If

  Else If (key=='HE') Then

    If (numi==1) Then
      If (nivel<=3) Then
        agau = a(nivel) + b(nivel)/xnu + c(nivel)/xnu/xnu
        Return
      Else
        agau = 1.D0
        Return
      End If

    Else If (numi==2) Then
      If (nivel<=3) Then
        agau = a(nivel) + 4.D0*b(nivel)/xnu + 16.D0*c(nivel)/xnu/xnu
        Return
      Else
        agau = 1.D0
        Return
      End If
    Else
      Write (999, *) ' STOP: ERROR IN NUMI -HE- GBF'
      Stop ' ERROR IN NUMI -HE- GBF'
    End If

  Else
    Write (999, *) ' STOP: AGAU INCORRECTLY CALLED, ONLY FOR H, HE'
    Stop ' AGAU INCORRECTLY CALLED, ONLY FOR H, HE'
  End If

End Function

!-----------------------------------------------------------------------

Function gaunt(n, f)

  Use :: nlte_type
  Use :: nlte_dim
  Implicit None

! -----bound free gaunt factors for hydrogen like atoms, f = nue/z**2

! .. scalar arguments ..
  Real (dp) :: f, gaunt
  Integer (i4b) :: n
! ..
! .. local scalars ..
  Real (dp) :: of, one
! ..
! .. local arrays ..
  Real (dp) :: a1(6), a2(6), b1(6), b2(6), c1(6), c2(6)
! ..
! .. save statement ..
  Save
! ..
! .. data statements ..
  Data a1/1.078D0, 1.0926D0, 1.0983D0, 1.0954D0, 1.0912D0, 1.0925D0/
  Data b1/ -8.754D14, -2.019D14, -9.45D13, -5.188D13, -3.2D13, -2.331D13/
  Data c1/ -1.791D29, 1.836D28, 9.177D27, 3.552D27, 1.576D27, 9.325D26/
  Data a2/0.798D0, 0.768D0, 0.793D0, 0.831D0, 0.758D0, 0.79D0/
  Data b2/5.358D15, 6.242D15, 5.48D15, 4.094D15, 6.633D15, 5.808D15/
  Data c2/ -3.484D31, -3.208D31, -2.318D31, -1.43D31, -3.32D31, -2.844D31/
  Data one/1.0D0/
! ..

  If (n>6) Go To 110
  of = one/f

  If (f>1.0D16) Go To 100

! -----less than 1.0e16 hz

  gaunt = a1(n) + (b1(n)+c1(n)*of)*of
  Return

! -----greater than 1.0e16 hz

100 Continue

  gaunt = a2(n) + (b2(n)+c2(n)*of)*of
  Return

! -----default

110 Continue

  gaunt = one
  Return

End Function

!-----------------------------------------------------------------------

Function gff(frq, t, z, nt, l, i)

  Use :: nlte_type
  Use :: nlte_dim
  Use :: nlte_var, Only: kis => ionmax
  Implicit None

! -----generates thermally averaged free-free non-relativistic gaunt
! -----factor for a hydrogenic ion of charge z, with a maximum relative
! -----error of .oo7, (rms fitting error = .001) for temperature and
! -----frequency in intervals:
! 10**-4 le u le 10**1.5,
! 10**-3 le gams le 10**3,
! -----where
! gams = z**2*ryd/k*t = z**2*1.579e5/t
! max          t=1.6e8k
! min          t=5000.k

! u = h*nu/k*t = 1.44*fr(cm-1)/t
! max  .17      t=2400.k
! 10.      t=140ooo.k
! min  3.2e7    t=1.4e6k
! 1.e6     t=45000.k
! 1.1e5    t=5000.k

! -----to obtain the
! -----stated accuracy, the full number of significant figures in the
! -----coefficients must be retained.

! -----this subroutine uses a two-dimensional chebyshev expansion
! -----computed from expressions given by karzas and latter (ap.j.
! -----suppl., v.6, p.167, 1961) augmented by various limiting forms
! -----of energy-specific gaunt factors.
! d.g.hummer, jila, may 1987.

! subroutine arguments are:
! z    =  nuclear charge
! t    =  temperature in degrees kelvin
! frq  =  frequency

! called by opacitc and opacitm

! adapted by A. Pauldrach and JP

! .. parameters ..
  Integer (i4b), Parameter :: nd1 = id_ndept, kisnd = kis*nd1, &
    kisnum = 8*kisnd
! ..
! .. scalar arguments ..
  Real (dp) :: frq, t, z, gff
  Integer (i4b) :: i, l, nt
! ..
! .. local scalars ..
  Real (dp) :: txu, xlf
  Integer (i4b) :: ir, j
! ..
! .. local arrays ..
  Real (dp) :: b(8), c(8, kis, nd1), con(kis, nd1)
! ..
! .. external functions ..
  Real (dp) :: gfree
  External :: gfree
! ..
! .. external subroutines ..
  External :: gffinit
! ..
! .. intrinsic functions ..
! INTRINSIC ABS,LOG10
! ..
! .. data statements ..

  Data con, c/kisnd*0.D0, kisnum*0.D0/
! ..

! for first freq., calculate freq. independent part of gaunt-factor

  If (l==1) Call gffinit(z, t, c(:,i,nt), con(i,nt))

  If (con(i,nt)==100.D0) Then
    gff = gfree(frq, t, z)
    Return

  Else

!   -------------------------------------------------------------
!   h [erg/s] / 1 ryd = 3.0397e-16 ; 1 ryd = 2.18e-11 erg

    xlf = log10(frq*3.0397D-16)
!   -------------------------------------------------------------
    txu = 0.72727273D0*xlf + con(i, nt)

    If (abs(txu)>2.0D0) Then
      gff = gfree(frq, t, z)
      Return
    End If

    ir = 6
    b(8) = c(8, i, nt)
    b(7) = txu*b(8) + c(7, i, nt)
    Do j = 1, 6
      b(ir) = txu*b(ir+1) - b(ir+2) + c(ir, i, nt)
      ir = ir - 1
    End Do

    gff = b(1) - b(3)
  End If

  Return

End Function

!-----------------------------------------------------------------------

Subroutine gffinit(z, t, c, con)

  Use :: nlte_type
  Use :: nlte_dim
  Implicit None

! -----initiation of the nue-independent part of the f-f gaunt-factor
! -----called by 'gff' - m 4.91 adi

! .. scalar arguments ..
  Real (dp) :: con, t, z
! ..
! .. array arguments ..
  Real (dp) :: c(8)
! ..
! .. local scalars ..
  Real (dp) :: txg, xlrkt
  Integer (i4b) :: i, ir, j
! ..
! .. local arrays ..
  Real (dp) :: b(11), d(8, 11)
! ..
! .. intrinsic functions ..
! INTRINSIC ABS,LOG10
! ..
! .. data statements ..

  Data ((d(i,j),i=1,8), j=1, 6)/8.986940175D+00, -4.009515855D+00, &
    8.808871266D-01, 2.640245111D-02, -4.580645915D-02, -3.568055702D-03, &
    2.827798067D-03, 3.365860195D-04, -8.006936989D-01, 9.466021705D-01, &
    9.043402532D-02, -9.608451450D-02, -1.885629865D-02, 1.050313890D-02, &
    2.800889961D-03, -1.078209202D-03, -3.781305103D-01, 1.102726332D-01, &
    -1.543619180D-02, 8.310561114D-03, 2.179620525D-02, 4.259726289D-03, &
    -4.181588794D-03, -1.770208330D-03, 1.877213132D-02, -1.004885705D-01, &
    -5.483366378D-02, -4.520154409D-03, 8.366530426D-03, 3.700273930D-03, &
    6.889320423D-04, 9.460313195D-05, 7.300158392D-02, 3.576785497D-03, &
    -4.545307025D-03, -1.017965604D-02, -9.530211924D-03, -3.450186162D-03, &
    1.040482914D-03, 1.407073544D-03, -1.744671550D-03, 2.864013856D-02, &
    1.903394837D-02, 7.091074494D-03, -9.668371391D-04, -2.999107465D-03, &
    -1.820642230D-03, -3.874082085D-04/
  Data ((d(i,j),i=1,8), j=7, 11)/ -1.707268366D-02, -4.694254776D-03, &
    1.311691517D-03, 5.316703136D-03, 5.178193095D-03, 2.451228935D-03, &
    -2.277321615D-05, -8.182359057D-04, 2.567331664D-04, -9.155339970D-03, &
    -6.997479192D-03, -3.571518641D-03, -2.096101038D-04, 1.553822487D-03, &
    1.509584686D-03, 6.212627837D-04, 4.098322531D-03, 1.635218463D-03, &
    -5.918883504D-04, -2.333091048D-03, -2.484138313D-03, -1.359996060D-03, &
    -5.371426147D-05, 5.553549563D-04, 3.837562402D-05, 2.938325230D-03, &
    2.393747064D-03, 1.328839809D-03, 9.135013312D-05, -7.137252303D-04, &
    -7.656848158D-04, -3.504683798D-04, -8.491991820D-04, -3.615327726D-04, &
    3.148015257D-04, 8.909207650D-04, 9.869737522D-04, 6.134671184D-04, &
    1.068883394D-04, -2.046080100D-04/
! ..

! -----xlrxt is log(ryd/kt)

  xlrkt = 5.1983649D0 - log10(t)
  txg = 0.66666667D0*(2.0D0*log10(z)+xlrkt)

  If (abs(txg)>2.0D0) Then
    con = 100.D0
    Return
  End If

  con = 0.72727273D0*xlrkt + 0.90909091D0

  Do j = 1, 8
    ir = 9
    b(11) = d(j, 11)
    b(10) = txg*b(11) + d(j, 10)

    Do i = 1, 9
      b(ir) = txg*b(ir+1) - b(ir+2) + d(j, ir)
      ir = ir - 1
    End Do

    c(j) = 0.25D0*(b(1)-b(3))
  End Do

  Return

End Subroutine

!-----------------------------------------------------------------------

Function gfree(frq, t, z)

  Use :: nlte_type
  Use :: nlte_dim
  Implicit None

! -----free-free-gaunt-factors
! -----called by 'gff'

! -----------------------------------------------------------------------
! .. scalar arguments ..
  Real (dp) :: frq, t, z, gfree
! ..
! .. local scalars ..
  Real (dp) :: c1, c2, c3, c4, thet, x
! ..

  thet = 5.0404D3/t
  thet = thet*z*z
  If (thet<4.0D-2) thet = 4.0D-2
  x = frq/2.997929D14
  x = x/z/z
  If (x>1.0D0) Go To 100
  If (x<0.2D0) x = 0.2D0

  gfree = (1.0823D0+2.98D-2/thet) + (6.7D-3+1.12D-2/thet)/x

  Return

100 Continue

  c1 = (3.9999187D-3-7.8622889D-5/thet)/thet + 1.070192D0
  c2 = (6.4628601D-2-6.1953813D-4/thet)/thet + 2.6061249D-1
  c3 = (1.3983474D-5/thet+3.7542343D-2)/thet + 5.7917786D-1
  c4 = 3.4169006D-1 + 1.1852264D-2/thet

  gfree = ((c4/x-c3)/x+c2)/x + c1

  Return

End Function

!-----------------------------------------------------------------------

Function hcol(nlo, nup, t)

  Use :: nlte_type
  Use :: nlte_dim
  Implicit None

! generates hydrogen collision rate coefficients
! using fits given by giovanardi, natta, and
! palla,  aa. suppl., 70, 269, 1987
! coefficients have been summed over all states

! nlo  =  lower principal quantum number (le.14)
! nup  =  upper principal quantum number (le.15)
! temp =  temperature  (2000 < t < 500000)

! d.g.hummer, jila, nov 1989
! modified k. butler, munich, jan 1990
! extended to higher temperatures and modified for use
! in detail (not checked so far)
! modified by j. puls for different data
! used by david for lowest transitions
! ly-alpha,ly-beta, balmer-alpha according to callaway
! (at.data nucl.data tables 57, 9-20


! .. scalar arguments ..
  Real (dp) :: t, hcol
  Integer (i4b) :: nlo, nup
! ..
! .. local scalars ..
  Real (dp) :: cutoff, gamfac, glo, gup, ryd, tr, tx, xnlos, xnups
  Integer (i4b) :: i, k
! ..
! .. local arrays ..
  Real (dp) :: bba(0:6), bla(0:5, 2), blb(0:7, 3), c(420), d(420), damp(2)
  Integer (i4b) :: n0(14)
! ..
! .. intrinsic functions ..
! INTRINSIC DBLE,EXP,LOG,SQRT
! ..
! .. data statements ..

! n=1 lower state
! n=2 lower state
! n=3 lower state
! n=4 lower state
! n=5 lower state
! n=6 lower state
! n=7 lower state
! n=8 lower state
! n=9 lower state
! n=10 lower state
! n=11 lower state
! n=12 lower state
! n=13 lower state
! n=14 lower state
! for each successive lower state, data start at n0:
  Data ryd/6.336D-6/
  Data damp/28.056D0, 7.2945D0/
  Data bla/0.36105477D0, 1.73103797D0, 0.37297254D0, -0.43240114D0, &
    0.14296095D0, -0.01670572D0, 0.25501041D0, 0.09973588D0, -0.05404749D0, &
    0.02520311D0, -0.00405306D0, 0.03508995D0/
  Data blb/0.09895325D0, 0.23240113D0, 0.20112555D0, 0.87813331D0, &
    -3.34239518D0, 3.70755260D0, -1.37320827D0, 0.D0, 0.06461453D0, &
    -0.11160826D0, 0.76957269D0, -2.18918875D0, 3.14216481D0, -2.23498961D0, &
    0.62651998D0, 0.01125150D0, 0.05742134D0, -0.00609531D0, 0.93848281D0, &
    -3.40393911D0, 5.50742903D0, -4.25196238D0, 1.26871040D0, 0.00724136D0/
  Data bba/21.6781D0, 207.5215D0, -47.581D0, 5.0449D0, -1.0262D0, 0.4085D0, &
    -5.93D-2/
  Data (c(i), i=1, 56)/5.732D-01, 1.829D-05, -1.158D-10, 9.429D-16, 2.122D-01, &
    -9.190D-07, 8.772D-11, -5.679D-16, 6.323D-03, 2.237D-06, -1.620D-11, &
    8.955D-17, 2.035D-02, 6.076D-07, -2.175D-13, -2.495D-18, 1.136D-02, &
    3.428D-07, -1.467D-13, -1.300D-18, 6.999D-03, 2.126D-07, -9.963D-14, &
    -7.672D-19, 4.624D-03, 1.410D-07, -6.969D-14, -4.927D-19, 3.217D-03, &
    9.836D-08, -5.031D-14, -3.361D-19, 2.329D-03, 7.135D-08, -3.737D-14, &
    -2.400D-19, 1.741D-03, 5.342D-08, -2.845D-14, -1.775D-19, 1.336D-03, &
    4.103D-08, -2.213D-14, -1.351D-19, 1.048D-03, 3.220D-08, -1.754D-14, &
    -1.053D-19, 8.369D-04, 2.574D-08, -1.413D-14, -8.368D-20, 6.791D-04, &
    2.090D-08, -1.154D-14, -6.763D-20/
  Data (d(i), i=1, 56)/5.856D-01, 1.551D-05, -9.669D-12, 5.716D-19, 1.537D-01, &
    3.548D-06, -3.224D-12, -7.626D-18, 2.400D-02, 1.419D-06, -2.008D-12, &
    1.356D-18, 2.002D-02, 6.325D-07, -7.070D-13, 4.096D-19, 1.123D-02, &
    3.549D-07, -3.998D-13, 2.331D-19, 6.940D-03, 2.194D-07, -2.483D-13, &
    1.453D-19, 4.593D-03, 1.453D-07, -1.648D-13, 9.667D-20, 3.199D-03, &
    1.012D-07, -1.150D-13, 6.758D-20, 2.318D-03, 7.334D-08, -8.349D-14, &
    4.910D-20, 1.727D-03, 5.493D-08, -6.270D-14, 3.695D-20, 1.326D-03, &
    4.218D-08, -4.821D-14, 2.844D-20, 1.040D-03, 3.310D-08, -3.786D-14, &
    2.236D-20, 8.305D-04, 2.645D-08, -3.028D-14, 1.790D-20, 6.740D-04, &
    2.147D-08, -2.460D-14, 1.455D-20/
! N=2 LOWER STATE
  Data (c(i), i=57, 108)/1.515D+01, 1.268D-03, 5.808D-09, -5.831D-14, &
    7.816D-01, 5.413D-04, -1.827D-09, 5.100D-17, 1.459D+00, 2.858D-04, &
    -2.207D-09, 9.028D-15, 7.172D-01, 1.440D-04, -1.139D-09, 4.755D-15, &
    4.107D-01, 8.360D-05, -6.699D-10, 2.823D-15, 2.591D-01, 5.319D-05, &
    -4.293D-10, 1.819D-15, 1.747D-01, 3.608D-05, -2.925D-10, 1.243D-15, &
    1.237D-01, 2.567D-05, -2.087D-10, 8.891D-16, 9.097D-02, 1.893D-05, &
    -1.543D-10, 6.585D-16, 6.896D-02, 1.438D-05, -1.174D-10, 5.017D-16, &
    5.356D-02, 1.119D-05, -9.150D-11, 3.913D-16, 4.247D-02, 8.887D-06, &
    -7.272D-11, 3.112D-16, 3.425D-02, 7.176D-06, -5.877D-11, 2.516D-16/
  Data (d(i), i=57, 108)/1.710D+01, 1.530D-03, -2.553D-09, 1.924D-15, &
    8.237D+00, 3.554D-04, -7.566D-10, 6.420D-16, 5.932D+00, 1.301D-04, &
    -2.912D-10, 2.535D-16, 2.991D+00, 6.419D-05, -1.444D-10, 1.260D-16, &
    1.733D+00, 3.689D-05, -8.324D-11, 7.267D-17, 1.102D+00, 2.333D-05, &
    -5.273D-11, 4.605D-17, 7.472D-01, 1.576D-05, -3.567D-11, 3.116D-17, &
    5.312D-01, 1.118D-05, -2.532D-11, 2.212D-17, 3.919D-01, 8.232D-06, &
    -1.865D-11, 1.630D-17, 2.977D-01, 6.245D-06, -1.416D-11, 1.237D-17, &
    2.312D-01, 4.855D-06, -1.101D-11, 9.622D-18, 1.838D-01, 3.851D-06, &
    -8.734D-12, 7.635D-18, 1.484D-01, 3.108D-06, -7.050D-12, 6.164D-18/
! N=3 LOWER STATE
  Data (c(i), i=109, 156)/ -1.290D+01, 2.059D-02, 5.469D-08, -9.082D-13, &
    3.562D+02, 7.337D-03, -9.622D-08, 5.596D-13, 5.744D+00, 3.570D-03, &
    -3.259D-08, 1.452D-13, 2.968D+00, 1.813D-03, -1.703D-08, 7.744D-14, &
    1.756D+00, 1.065D-03, -1.016D-08, 4.667D-14, 1.135D+00, 6.865D-04, &
    -6.601D-09, 3.053D-14, 7.802D-01, 4.713D-04, -4.558D-09, 2.116D-14, &
    5.615D-01, 3.390D-04, -3.292D-09, 1.532D-14, 4.189D-01, 2.528D-04, &
    -2.461D-09, 1.148D-14, 3.213D-01, 1.939D-04, -1.891D-09, 8.833D-15, &
    2.523D-01, 1.522D-04, -1.487D-09, 6.953D-15, 2.018D-01, 1.218D-04, &
    -1.192D-09, 5.576D-15/
  Data (d(i), i=109, 156)/1.940D+02, 1.949D-02, -3.832D-08, 3.137D-14, &
    4.729D+02, 1.927D-03, -4.171D-09, 3.628D-15, 6.741D+01, 1.315D-03, &
    -3.145D-09, 2.814D-15, 3.444D+01, 6.477D-04, -1.560D-09, 1.399D-15, &
    2.031D+01, 3.744D-04, -9.054D-10, 8.130D-16, 1.311D+01, 2.388D-04, &
    -5.789D-10, 5.203D-16, 9.007D+00, 1.629D-04, -3.955D-10, 3.556D-16, &
    6.484D+00, 1.166D-04, -2.835D-10, 2.550D-16, 4.837D+00, 8.666D-05, &
    -2.108D-10, 1.896D-16, 3.711D+00, 6.631D-05, -1.614D-10, 1.452D-16, &
    2.914D+00, 5.194D-05, -1.265D-10, 1.138D-16, 2.332D+00, 4.150D-05, &
    -1.011D-10, 9.100D-17/
! N=4 LOWER STATE
  Data (c(i), i=157, 200)/4.139D+03, 4.645D-01, -7.097D-06, 4.388D-11, &
    1.794D+03, 4.443D-02, -6.484D-07, 3.936D-12, 1.536D+01, 2.042D-02, &
    -2.065D-07, 9.734D-13, 8.730D+00, 1.033D-02, -1.074D-07, 5.161D-13, &
    5.434D+00, 6.084D-03, -6.423D-08, 3.116D-13, 3.628D+00, 3.938D-03, &
    -4.196D-08, 2.048D-13, 2.554D+00, 2.718D-03, -2.914D-08, 1.428D-13, &
    1.873D+00, 1.967D-03, -2.119D-08, 1.041D-13, 1.418D+00, 1.476D-03, &
    -1.594D-08, 7.843D-14, 1.102D+00, 1.138D-03, -1.232D-08, 6.075D-14, &
    8.744D-01, 8.987D-04, -9.746D-09, 4.809D-14/
  Data (d(i), i=157, 200)/7.204D+03, 1.627D-01, -5.181D-07, 5.605D-13, &
    2.507D+03, 9.370D-03, -2.091D-08, 1.842D-14, 3.823D+02, 6.480D-03, &
    -1.600D-08, 1.452D-14, 1.950D+02, 3.161D-03, -7.869D-09, 7.157D-15, &
    1.154D+02, 1.823D-03, -4.561D-09, 4.154D-15, 7.486D+01, 1.165D-03, &
    -2.924D-09, 2.665D-15, 5.178D+01, 7.977D-04, -2.006D-09, 1.830D-15, &
    3.752D+01, 5.737D-04, -1.444D-09, 1.318D-15, 2.816D+01, 4.283D-04, &
    -1.080D-09, 9.858D-16, 2.174D+01, 3.293D-04, -8.307D-10, 7.587D-16, &
    1.717D+01, 2.592D-04, -6.544D-10, 5.978D-16/
! N=5 LOWER STATE
  Data (c(i), i=201, 240)/ -9.122D+02, 1.260D+00, -1.070D-05, 4.290D-11, &
    3.959D+01, 2.108D-01, -2.162D-06, 1.020D-11, 3.691D+01, 7.806D-02, &
    -8.485D-07, 4.166D-12, 2.352D+01, 3.911D-02, -4.365D-07, 2.179D-12, &
    1.542D+01, 2.296D-02, -2.601D-07, 1.310D-12, 1.062D+01, 1.487D-02, &
    -1.699D-07, 8.608D-13, 7.642D+00, 1.029D-02, -1.183D-07, 6.014D-13, &
    5.695D+00, 7.464D-03, -8.621D-08, 4.394D-13, 4.368D+00, 5.617D-03, &
    -6.508D-08, 3.323D-13, 3.430D+00, 4.348D-03, -5.051D-08, 2.583D-13/
  Data (d(i), i=201, 240)/2.166D+04, 4.690D-01, -1.122D-06, 1.008D-12, &
    3.874D+03, 6.443D-02, -1.596D-07, 1.452D-13, 1.465D+03, 2.207D-02, &
    -5.556D-08, 5.082D-14, 7.410D+02, 1.062D-02, -2.698D-08, 2.476D-14, &
    4.374D+02, 6.096D-03, -1.556D-08, 1.430D-14, 2.841D+02, 3.889D-03, &
    -9.962D-09, 9.167D-15, 1.969D+02, 2.663D-03, -6.838D-09, 6.297D-15, &
    1.431D+02, 1.918D-03, -4.935D-09, 4.547D-15, 1.078D+02, 1.436D-03, &
    -3.698D-09, 3.409D-15, 8.353D+01, 1.107D-03, -2.854D-09, 2.632D-15/
! N=6 LOWER STATE
  Data (c(i), i=241, 276)/ -3.431D+03, 4.116D+00, -3.853D-05, 1.679D-10, &
    4.397D+01, 6.434D-01, -7.008D-06, 3.431D-11, 8.927D+01, 2.325D-01, &
    -2.667D-06, 1.350D-11, 6.153D+01, 1.152D-01, -1.354D-06, 6.957D-12, &
    4.165D+01, 6.729D-02, -8.024D-07, 4.156D-12, 2.923D+01, 4.349D-02, &
    -5.232D-07, 2.724D-12, 2.130D+01, 3.008D-02, -3.641D-07, 1.902D-12, &
    1.603D+01, 2.185D-02, -2.656D-07, 1.391D-12, 1.239D+01, 1.647D-02, &
    -2.008D-07, 1.054D-12/
  Data (d(i), i=241, 276)/7.146D+04, 1.379D+00, -3.346D-06, 3.023D-12, &
    1.187D+04, 1.794D-01, -4.501D-07, 4.118D-13, 4.380D+03, 5.990D-02, &
    -1.527D-07, 1.405D-13, 2.192D+03, 2.846D-02, -7.324D-08, 6.759D-14, &
    1.288D+03, 1.621D-02, -4.197D-08, 3.881D-14, 8.351D+02, 1.031D-02, &
    -2.678D-08, 2.480D-14, 5.790D+02, 7.050D-03, -1.837D-08, 1.702D-14, &
    4.213D+02, 5.079D-03, -1.326D-08, 1.229D-14, 3.179D+02, 3.804D-03, &
    -9.944D-09, 9.226D-15/
! N=7 LOWER STATE
  Data (c(i), i=277, 308)/ -9.280D+03, 1.116D+01, -1.122D-04, 5.167D-10, &
    6.658D+01, 1.651D+00, -1.884D-05, 9.487D-11, 2.172D+02, 5.833D-01, &
    -6.977D-06, 3.615D-11, 1.535D+02, 2.858D-01, -3.499D-06, 1.838D-11, &
    1.049D+02, 1.660D-01, -2.060D-06, 1.090D-11, 7.412D+01, 1.070D-01, &
    -1.339D-06, 7.118D-12, 5.428D+01, 7.389D-02, -9.304D-07, 4.963D-12, &
    4.103D+01, 5.366D-02, -6.786D-07, 3.629D-12/
  Data (d(i), i=277, 308)/1.954D+05, 3.426D+00, -8.397D-06, 7.624D-12, &
    3.057D+04, 4.266D-01, -1.080D-06, 9.917D-13, 1.103D+04, 1.392D-01, &
    -3.582D-07, 3.309D-13, 5.458D+03, 6.530D-02, -1.696D-07, 1.572D-13, &
    3.189D+03, 3.693D-02, -9.653D-08, 8.966D-14, 2.062D+03, 2.338D-02, &
    -6.136D-08, 5.707D-14, 1.429D+03, 1.595D-02, -4.200D-08, 3.910D-14, &
    1.039D+03, 1.148D-02, -3.029D-08, 2.282D-14/
! N=8 LOWER STATE
  Data (c(i), i=309, 336)/ -2.069D+04, 2.637D+01, -2.802D-04, 1.342D-09, &
    2.055D+02, 3.731D+00, -4.420D-05, 2.276D-10, 5.123D+02, 1.292D+00, &
    -1.599D-05, 8.442D-11, 3.578D+02, 6.265D-01, -7.922D-06, 4.235D-11, &
    2.438D+02, 3.616D-01, -4.633D-06, 2.494D-11, 1.721D+02, 2.322D-01, &
    -3.000D-06, 1.622D-11, 1.260D+02, 1.601D-01, -2.081D-06, 1.129D-11/
  Data (d(i), i=309, 336)/4.651D+05, 7.527D+00, -1.859D-05, 1.694D-11, &
    6.930D+04, 9.038D-01, -2.302D-06, 2.121D-12, 2.450D+04, 2.891D-01, &
    -7.487D-07, 6.939D-13, 1.200D+04, 1.340D-01, -3.505D-07, 3.260D-13, &
    6.970D+03, 7.523D-02, -1.981D-07, 1.846D-13, 4.493D+03, 4.741D-02, &
    -1.254D-07, 1.170D-13, 3.106D+03, 3.226D-02, -8.559D-08, 7.997D-14/
! N=9 LOWER STATE
  Data (c(i), i=337, 360)/ -4.032D+04, 5.614D+01, -6.231D-04, 3.073D-09, &
    6.989D+02, 7.655D+00, -9.352D-05, 4.903D-10, 1.141D+03, 2.605D+00, &
    -3.313D-05, 1.777D-10, 7.755D+02, 1.250D+00, -1.624D-05, 8.808D-11, &
    5.234D+02, 7.175D-01, -9.437D-06, 5.153D-11, 3.677D+02, 4.590D-01, &
    -6.087D-06, 3.338D-11/
  Data (d(i), i=337, 360)/9.956D+05, 1.506D+01, -3.741D-05, 3.418D-11, &
    1.425D+05, 1.754D+00, -4.489D-06, 4.146D-12, 4.949D+04, 5.510D-01, &
    -1.435D-06, 1.333D-12, 2.401D+04, 2.526D-01, -6.645D-07, 6.196D-13, &
    1.386D+04, 1.408D-01, -3.729D-07, 3.485D-13, 8.904D+03, 8.835D-02, &
    -2.350D-07, 2.200D-13/
! N=10 LOWER STATE
  Data (c(i), i=361, 380)/ -7.097D+04, 1.101D+02, -1.266D-03, 6.390D-09, &
    2.018D+03, 1.455D+01, -1.824D-04, 9.708D-10, 2.383D+03, 4.875D+00, &
    -6.348D-05, 3.449D-10, 1.569D+03, 2.319D+00, -3.081D-05, 1.691D-10, &
    1.046D+03, 1.323D+00, -1.779D-05, 9.830D-11/
  Data (d(i), i=361, 380)/1.961D+06, 2.798D+01, -6.982D-05, 6.394D-11, &
    2.715D+05, 3.175D+00, -8.158D-06, 7.551D-12, 9.279D+04, 9.821D-01, &
    -2.567D-06, 2.390D-12, 4.460D+04, 4.458D-01, -1.178D-06, 1.100D-12, &
    2.561D+04, 2.468D-01, -6.566D-07, 6.150D-13/
! N=11 LOWER STATE
  Data (c(i), i=381, 396)/ -1.150D+05, 2.020D+02, -2.392D-03, 1.231D-08, &
    4.988D+03, 2.601D+01, -3.334D-04, 1.797D-09, 4.675D+03, 8.595D+00, &
    -1.142D-04, 6.273D-10, 2.986D+03, 4.054D+00, -5.491D-05, 3.046D-10/
  Data (d(i), i=381, 396)/3.613D+06, 4.898D+01, -1.227D-04, 1.125D-10, &
    4.861D+05, 5.434D+00, -1.401D-05, 1.299D-11, 1.638D+05, 1.658D+00, &
    -4.348D-06, 4.054D-12, 7.810D+04, 7.456D-01, -1.976D-06, 1.850D-12/
! N=12 LOWER STATE
  Data (c(i), i=397, 408)/ -1.737D+05, 3.511D+02, -4.263D-03, 2.227D-08, &
    1.094D+04, 4.419D+01, -5.774D-04, 3.146D-09, 8.673D+03, 1.442D+01, &
    -1.950D-04, 1.082D-09/
  Data (d(i), i=397, 408)/6.300D+06, 8.163D+01, -2.051D-04, 1.884D-10, &
    8.271D+05, 8.881D+00, -2.296D-05, 2.131D-11, 2.753D+05, 2.676D+00, &
    -7.037D-06, 6.571D-12/
! N=13 LOWER STATE
  Data (c(i), i=409, 416)/ -2.459D+05, 5.829D+02, -7.233D-03, 3.830D-08, &
    2.191D+04, 7.194D+01, -9.561D-04, 5.259D-09/
  Data (d(i), i=409, 416)/1.049D+07, 1.305D+02, -3.288D-04, 3.025D-10, &
    1.348D+06, 1.396D+01, -3.617D-05, 3.361D-11/
! N=14 LOWER STATE
  Data (c(i), i=417, 420)/ -3.273D+05, 9.312D+02, -1.178D-02, 6.306D-08/
  Data (d(i), i=417, 420)/1.680D+07, 2.016D+02, -5.089D-04, 4.687D-10/
! FOR EACH SUCCESSIVE LOWER STATE, DATA START AT N0:

  Data (n0(i), i=1, 14)/1, 57, 109, 157, 201, 241, 277, 309, 337, 361, 381, &
    397, 409, 417/
! ..

! check quantum number and return zero if no data available

  If (nlo==nup .Or. nlo>14 .Or. nup>15) Then
    Write (999, *) ' STOP: ERROR IN H-COLL RATES'
    Stop ' ERROR IN H-COLL RATES'
  End If

! set pointer for start of coefficients

  xnlos = dble(nlo**2)
  xnups = dble(nup**2)
  glo = 2.0D0*xnlos
  gup = 2.0D0*xnups

  If (nlo==1 .And. nup==2) Then

!   ly-alpha

    tr = t*ryd
    tx = 1.D0
    gamfac = bla(0, 1) + bla(0, 2)

    Do k = 1, 4
      tx = tx*tr
      gamfac = gamfac + (bla(k,1)+bla(k,2))*tx
    End Do

    cutoff = bla(5, 2)*log(damp(1)*tr)*exp(-damp(2)*tr)
    tx = tx*tr
    gamfac = gamfac + bla(5, 1)*tx + cutoff
    hcol = gamfac*8.631D-6/(sqrt(t)*glo)
    Return

  Else If (nlo==1 .And. nup==3) Then

!   ly-beta

    tr = t*ryd
    tx = 1.D0
    gamfac = blb(0, 1) + blb(0, 2) + blb(0, 3)

    Do k = 1, 6
      tx = tx*tr
      gamfac = gamfac + (blb(k,1)+blb(k,2)+blb(k,3))*tx
    End Do

    cutoff = (blb(7,2)+blb(7,3))*log(damp(1)*tr)*exp(-damp(2)*tr)
    gamfac = gamfac + cutoff
    hcol = gamfac*8.631D-6/(sqrt(t)*glo)
    Return

  Else If (nlo==2 .And. nup==3) Then

!   balmer-alpha (packed in advance)

    tr = t*ryd
    tx = 1.D0
    gamfac = bba(0)

    Do k = 1, 6
      tx = tx*tr
      gamfac = gamfac + bba(k)*tx
    End Do

    hcol = gamfac*8.631D-6/(sqrt(t)*glo)
    Return

  Else

!   other transitions according to italian mafia

!   index of 1st coef.
    k = n0(nlo) + (nup-nlo-1)*4

!   evaluate rates at specified temperatures

    If (t>60000) Then
      gamfac = 8.631D-6*(d(k)+t*(d(k+1)+t*(d(k+2)+t*d(k+3))))/sqrt(t)
    Else
      gamfac = 8.631D-6*(c(k)+t*(c(k+1)+t*(c(k+2)+t*c(k+3))))/sqrt(t)
    End If

    hcol = gamfac/glo
    Return

  End If

End Function

!-----------------------------------------------------------------------

Function hecol(nl, nu, t)

  Use :: nlte_type
  Use :: nlte_dim
  Implicit None

! for specified t, 1000.le.t.le.30000, gives de-excitation
! rate constants in array c(iu,il), where 11.gt.iu.ge.1.
! certain transitions are missing: (8,6), (9,6), (10,7),
! (11,7), and (9,8).  the states are coded as follows:
! 1: 1 sing s
! 2: 2 trip s
! 3: 2 sing s
! 4: 2 trip p
! 5: 2 sing p
! 6: 3 trip s
! 7: 3 sing s
! 8: 3 trip p
! 9: 3 trip d
! 10: 3 sing d
! 11: 3 sing p

! based on calculations by berrington and kingston (1987)
! d.g.hummer, july, 1987
! modified k. butler. for use in detail july 1987.
! range extended beyond 30000k by setting collision strength to
! value at 30000k.


! .. scalar arguments ..
  Real (dp) :: t, hecol
  Integer (i4b) :: nl, nu
! ..
! .. local scalars ..
  Real (dp) :: c0, cm, cp, d1, d2, d3, f1, f2, f3, st, x, x0, xm, xp
  Integer (i4b) :: i, i0, im, ip, k
! ..
! .. local arrays ..
  Real (dp) :: g(8, 51), g1(80), g2(72), g3(64), g4(56), g5(48), g6(40), &
    g7(16), g8(16), g9(16), tn(8), w(11)
  Integer (i4b) :: iumax(9), iumin(9), nstart(9)
! ..
! .. intrinsic functions ..
! INTRINSIC SQRT
! ..
! .. equivalences ..
  Equivalence (g(1,1), g1), (g(1,11), g2), (g(1,20), g3), (g(1,28), g4), &
    (g(1,35), g5), (g(1,41), g6), (g(1,46), g7), (g(1,48), g8), (g(1,50), g9)
! ..
! .. save statement ..

  Save
! ..
! .. data statements ..
  Data iumin/2, 3, 4, 5, 6, 7, 8, 10, 10/
  Data iumax/11, 11, 11, 11, 11, 11, 9, 11, 11/
  Data nstart/0, 10, 19, 27, 34, 40, 45, 47, 49/
  Data tn/1.0D3, 2.0D3, 5.0D3, 1.0D4, 1.5D4, 2.0D4, 2.5D4, 3.0D4/
  Data w/1.0D0, 3.0D0, 1.0D0, 9.0D0, 3.0D0, 3.0D0, 1.0D0, 9.0D0, 1.5D1, 5.0D0, &
    3.0D0/
  Data g1/3.09D-2, 4.95D-2, 6.82D-2, 7.24D-2, 7.22D-2, 7.16D-2, 7.10D-2, &
    7.03D-2, 1.60D-2, 2.36D-2, 3.31D-2, 3.83D-2, 4.08D-2, 4.26D-2, 4.39D-2, &
    4.49D-2, 6.37D-3, 1.01D-2, 1.71D-2, 2.44D-2, 2.92D-2, 3.30D-2, 3.61D-2, &
    3.84D-2, 3.02D-3, 5.29D-3, 1.05D-2, 1.63D-2, 2.09D-2, 2.52D-2, 2.91D-2, &
    3.24D-2, 1.82D-2, 1.87D-2, 1.82D-2, 1.81D-2, 1.88D-2, 1.95D-2, 2.00D-2, &
    2.03D-2, 9.36D-3, 9.42D-3, 9.18D-3, 9.62D-3, 1.06D-2, 1.16D-2, 1.23D-2, &
    1.28D-2, 5.86D-3, 6.19D-3, 7.01D-3, 8.40D-3, 1.01D-2, 1.17D-2, 1.30D-2, &
    1.38D-2, 8.94D-4, 1.67D-3, 2.34D-3, 2.60D-3, 2.80D-3, 2.93D-3, 2.99D-3, &
    3.01D-3, 4.49D-3, 5.05D-3, 5.19D-3, 5.42D-3, 5.93D-3, 6.41D-3, 6.75D-3, &
    6.95D-3, 6.56D-4, 1.10D-3, 2.25D-3, 4.24D-3, 6.74D-3, 9.21D-3, 1.13D-2, &
    1.29D-2/
  Data g2/1.01D+0, 1.60D+0, 2.24D+0, 2.40D+0, 2.33D+0, 2.21D+0, 2.09D+0, &
    1.98D+0, 3.30D+0, 5.59D+0, 1.36D+1, 2.62D+1, 3.69D+1, 4.58D+1, 5.30D+1, &
    5.86D+1, 3.17D-1, 4.79D-1, 7.60D-1, 9.71D-1, 1.06D+0, 1.09D+0, 1.10D+0, &
    1.09D+0, 2.73D+0, 2.89D+0, 2.67D+0, 2.47D+0, 2.47D+0, 2.54D+0, 2.61D+0, &
    2.65D+0, 4.81D-1, 4.99D-1, 4.62D-1, 4.03D-1, 3.67D-1, 3.41D-1, 3.19D-1, &
    3.01D-1, 1.79D+0, 1.84D+0, 1.84D+0, 1.83D+0, 1.90D+0, 1.98D+0, 2.04D+0, &
    2.07D+0, 4.09D-1, 7.45D-1, 1.37D+0, 2.14D+0, 2.87D+0, 3.51D+0, 4.00D+0, &
    4.35D+0, 2.42D-1, 2.73D-1, 2.98D-1, 3.16D-1, 3.34D-1, 3.47D-1, 3.53D-1, &
    3.53D-1, 4.34D-2, 8.14D-2, 1.42D-1, 1.69D-1, 1.80D-1, 1.84D-1, 1.85D-1, &
    1.84D-1/
  Data g3/9.96D-1, 1.22D+0, 1.51D+0, 1.71D+0, 1.76D+0, 1.76D+0, 1.73D+0, &
    1.68D+0, 1.10D+0, 2.97D+0, 9.16D+0, 1.81D+1, 2.53D+1, 3.11D+1, 3.56D+1, &
    3.89D+1, 8.73D-1, 8.94D-1, 7.48D-1, 5.88D-1, 5.01D-1, 4.44D-1, 4.03D-1, &
    3.70D-1, 6.20D-1, 6.17D-1, 6.23D-1, 6.47D-1, 7.11D-1, 7.86D-1, 8.53D-1, &
    9.03D-1, 5.90D-1, 6.03D-1, 6.12D-1, 5.60D-1, 5.14D-1, 4.78D-1, 4.48D-1, &
    4.22D-1, 9.61D-2, 1.67D-1, 2.79D-1, 3.56D-1, 4.08D-1, 4.42D-1, 4.60D-1, &
    4.68D-1, 5.29D-1, 6.44D-1, 8.73D-1, 1.25D+0, 1.67D+0, 2.06D+0, 2.36D+0, &
    2.59D+0, 1.14D-1, 1.88D-1, 2.93D-1, 3.60D-1, 4.11D-1, 4.56D-1, 4.90D-1, &
    5.14D-1/
  Data g4/1.92D+0, 2.60D+0, 3.57D+0, 4.38D+0, 4.64D+0, 4.75D+0, 4.72D+0, &
    4.72D+0, 7.43D+0, 7.11D+0, 6.25D+0, 5.93D+0, 6.13D+0, 6.44D+0, 6.69D+0, &
    6.86D+0, 6.51D-1, 8.04D-1, 1.07D+0, 1.06D+0, 9.68D-1, 8.91D-1, 8.26D-1, &
    7.70D-1, 8.59D+0, 1.05D+1, 1.23D+1, 1.37D+1, 1.53D+1, 1.68D+1, 1.81D+1, &
    1.89D+1, 3.53D+0, 5.71D+0, 9.81D+0, 1.54D+1, 2.13D+1, 2.67D+1, 3.10D+1, &
    3.42D+1, 1.33D+0, 1.50D+0, 1.60D+0, 1.62D+0, 1.67D+0, 1.70D+0, 1.70D+0, &
    1.68D+0, 3.41D-1, 5.22D-1, 7.96D-1, 9.14D-1, 9.45D-1, 9.47D-1, 9.33D-1, &
    9.10D-1/
  Data g5/1.21D+0, 1.46D+0, 1.38D+0, 1.16D+0, 1.03D+0, 9.37D-1, 8.67D-1, &
    8.09D-1, 6.55D-1, 7.39D-1, 8.04D-1, 8.84D-1, 1.01D+0, 1.14D+0, 1.24D+0, &
    1.32D+0, 1.55D+0, 1.77D+0, 1.87D+0, 1.76D+0, 1.66D+0, 1.58D+0, 1.50D+0, &
    1.42D+0, 5.71D-1, 8.47D-1, 1.24D+0, 1.52D+0, 1.70D+0, 1.80D+0, 1.85D+0, &
    1.85D+0, 3.23D+0, 3.96D+0, 5.19D+0, 7.23D+0, 9.68D+0, 1.20D+1, 1.39D+1, &
    1.54D+1, 6.12D-1, 9.70D-1, 1.57D+0, 2.22D+0, 2.90D+0, 3.53D+0, 4.03D+0, &
    4.41D+0/
  Data g6/3.60D+0, 4.03D+0, 3.68D+0, 2.92D+0, 2.44D+0, 2.12D+0, 1.88D+0, &
    1.70D+0, 0.00D+0, 0.00D+0, 0.00D+0, 0.00D+0, 0.00D+0, 0.00D+0, 0.00D+0, &
    0.00D+0, 0.00D+0, 0.00D+0, 0.00D+0, 0.00D+0, 0.00D+0, 0.00D+0, 0.00D+0, &
    0.00D+0, 3.11D+0, 3.26D+0, 3.10D+0, 2.78D+0, 2.55D+0, 2.36D+0, 2.20D+0, &
    2.06D+0, 5.54D-1, 9.44D-1, 1.40D+0, 1.50D+0, 1.45D+0, 1.38D+0, 1.31D+0, &
    1.23D+0/
  Data g7/3.36D+0, 3.80D+0, 3.70D+0, 3.20D+0, 2.81D+0, 2.52D+0, 2.28D+0, &
    2.09D+0, 4.00D-1, 8.69D-1, 1.54D+0, 1.76D+0, 1.77D+0, 1.72D+0, 1.64D+0, &
    1.56D+0/
  Data g8/1.22D+1, 1.27D+1, 1.20D+1, 1.06D+1, 9.58D+0, 8.75D+0, 8.06D+0, &
    7.46D+0, 2.72D+0, 3.71D+0, 4.85D+0, 5.27D+0, 5.26D+0, 5.08D+0, 4.84D+0, &
    4.58D+0/
  Data g9/6.08D+0, 8.78D+0, 1.19D+1, 1.32D+1, 1.37D+1, 1.38D+1, 1.36D+1, &
    1.32D+1, 1.29D+0, 2.64D+0, 4.66D+0, 5.55D+0, 5.75D+0, 5.70D+0, 5.54D+0, &
    5.33D+0/
! ..

! check that temperature is in bounds

  If (nl<1 .Or. nl>9) Then
    Write (*, Fmt=*) 'NL OUT OF RANGE IN HECOL.  NU = ', nu, ' NL = ', nl
    Write (999, *) ' STOP!'
    Stop
  End If

  If (nu<iumin(nl) .Or. nu>iumax(nl)) Then
    Write (*, Fmt=*) 'NU OUT OF RANGE IN HECOL.  NU = ', nu, ' NL = ', nl
    Write (999, *) ' STOP!'
    Stop
  End If

  i = nstart(nl) + nu - nl
  st = sqrt(t)

  If (t<tn(1)) Then
    hecol = 8.63D-6*g(1, i)/(w(nl)*st)
    Return
  End If

  If (t>tn(8)) Then
    hecol = 8.63D-6*g(8, i)/(w(nl)*st)
    Return
  End If

! find interpolation region

  k = 2

100 Continue

  If (tn(k)<t) Then
    k = k + 1
    Go To 100
  End If

  If (k==2) Then
    im = 1
    i0 = 2
    ip = 3
  Else
    im = k - 2
    i0 = k - 1
    ip = k
  End If

! calculate interpolation weights

  x = t
  xm = tn(im)
  x0 = tn(i0)
  xp = tn(ip)
  d1 = x0 - xm
  d2 = xp - x0
  d3 = xp - xm
  f1 = x - xm
  f2 = x - x0
  f3 = x - xp
  cm = f2*f3/(d1*d3)
  c0 = -f1*f3/(d1*d2)
  cp = f1*f2/(d2*d3)

! interpolate averaged collision strengths and
! evaluate downwards collision rate coefficients

  hecol = g(im, i)*cm + g(i0, i)*c0 + g(ip, i)*cp
  hecol = 8.63D-6*hecol/(w(nl)*st)

  Return

End Function

!***********************************************************************

!subroutines: miscellaneous

!***********************************************************************

Function erfcm(x)

  Use :: nlte_type
  Use :: nlte_dim
  Implicit None

! ---- complementary error function of a given value (x)
! ---- modified (e. santolaya, 23/7/93) so that it returns not
! ---- erfc(x) but erfc(x)*exp(x**2)

! .. scalar arguments ..
  Real (dp) :: x, erfcm
! ..
! .. local scalars ..
  Real (dp) :: a1, a2, a3, a4, a5, p, t
! ..
! .. data statements ..

  Data p/.3275911D0/
  Data a1, a2, a3, a4, a5/.254829592D0, -.284496736D0, 1.421413741D0, &
    -1.453152027D0, 1.061405429D0/
! ..

  t = 1.D0/(1.D0+p*x)
  erfcm = t*(a1+t*(a2+t*(a3+t*(a4+t*a5))))

  If (erfcm<0.D0) Then
    Write (999, *) ' STOP: ERROR IN ERFCM'
    Stop ' ERROR IN ERFCM'
  End If

  Return

End Function

!-----------------------------------------------------------------------

Function expin(n, x)

  Use :: nlte_type
  Use :: nlte_dim
  Implicit None

! returns exponential integral(degree n>1) up to factor exp(-x)

! warning!!! this function does not return e_n(x) but
! e_n(x)*exp(x). the result must be multiplied by exp(-x) to
! get e_n(x)

! .. scalar arguments ..
  Real (dp) :: x, expin
  Integer (i4b) :: n
! ..
! .. local scalars ..
  Real (dp) :: summ, xn, xnq
  Integer (i4b) :: k
! ..
! .. external functions ..
  Real (dp) :: expin1, fakul
  External :: expin1, fakul
! ..

  If (x==0.) Then
    expin = 1.D0/(float(n)-1.)
    Return
  End If

  If (n>=20 .Or. (n>=8 .And. x>5.D0)) Then
    xn = x + n
    xnq = xn*xn
    summ = 1.D0 + n/xnq + n*(n-2.D0*x)/xnq/xnq + n*(6.D0*x*x-8.D0*n*x+n*n)/xnq &
      /xnq/xnq
    expin = summ/xn
    Return
  End If

  summ = 0.D0

  Do k = 0, n - 2
    summ = summ + (-1)**k*x**k*fakul(n-k-2)
  End Do

  summ = summ + (-1)**(n-1)*x**(n-1)*expin1(x)
  expin = summ/fakul(n-1)

  Return

End Function

!-----------------------------------------------------------------------

Function expin1(x)

  Use :: nlte_type
  Use :: nlte_dim
  Implicit None

! returns first exponential integral up to factor exp(-x)

! warning!!! this function does not return e_1(x) but
! e_1(x)*exp(x). the result must be multiplied by exp(-x) to
! get e_1(x)


! .. scalar arguments ..
  Real (dp) :: x, expin1
! ..
! .. local scalars ..
  Real (dp) :: a0, a1, a2, a3, a4, a5, b1, b2, b3, b4, c1, c2, c3, c4, sum1, &
    sum2
! ..
! .. intrinsic functions ..

! INTRINSIC EXP,LOG
! ..
! .. data statements ..
  Data a0, a1, a2, a3, a4, a5/ -.57721566D0, .99999193D0, -.24991055D0, &
    .05519968D0, -.00976004D0, .00107857D0/
  Data b1, b2, b3, b4/8.5733287401D0, 18.0590169730D0, 8.6347608925D0, &
    .2677737343D0/
  Data c1, c2, c3, c4/9.5733223454D0, 25.6329561486D0, 21.0996530827D0, &
    3.9584969228D0/
! ..

  If (x<=1.D0) Then
    expin1 = a0 + x*(a1+x*(a2+x*(a3+x*(a4+x*a5)))) - log(x)
    expin1 = expin1*(exp(x))
    Return
  End If

  sum1 = b4 + x*(b3+x*(b2+x*(b1+x)))
  sum2 = c4 + x*(c3+x*(c2+x*(c1+x)))
  expin1 = sum1/(sum2*x)

  Return

End Function

!-----------------------------------------------------------------------

Function fakul(n)

  Use :: nlte_type
  Use :: nlte_dim
  Implicit None

! .. scalar arguments ..
  Real (dp) :: fakul
  Integer (i4b) :: n
! ..
! .. local scalars ..
  Real (dp) :: xj
  Integer (i4b) :: i
! ..
  If (n<0) Then
    Write (999, *) ' STOP: ARGUMENT OF FAKULTAET NEGATIVE!'
    Stop ' ARGUMENT OF FAKULTAET NEGATIVE!'
  End If

  If (n<=1) Then
    fakul = 1.D0
    Return
  End If

  xj = 1.D0

  Do i = 1, n
    xj = xj*i
  End Do

  fakul = xj

  Return
End Function

!-----------------------------------------------------------------------

Subroutine integra(a, summ, nd, nn, del, ii)

  Use :: nlte_type
  Use :: nlte_dim
  Implicit None

! .. parameters ..
  Integer (i4b), Parameter :: nd1 = id_ndept
! ..
! .. scalar arguments ..
  Real (dp) :: del
  Integer (i4b) :: ii, nd, nn
! ..
! .. array arguments ..
  Real (dp) :: a(nd1), summ(nd1)
! ..
! .. local scalars ..
  Integer (i4b) :: l
! ..

  Do l = 1, nd

    If (ii==1) Then
      summ(l) = a(l)*del
    Else
      summ(l) = summ(l) + a(l)*del
    End If

  End Do

  If (ii/=nn) Return

  a(1:nd) = summ(1:nd)

  Return

End Subroutine

!-----------------------------------------------------------------------

Subroutine linreg(x, y, n, xm, b, r)

  Use :: nlte_type
  Use :: nlte_dim
  Implicit None


! -----one-dimensional linear regression


! .. scalar arguments ..
  Real (dp) :: b, r, xm
  Integer (i4b) :: n
! ..
! .. array arguments ..
  Real (dp) :: x(n), y(n)
! ..
! .. local scalars ..
  Real (dp) :: sigx, sigy, sumx, sumxx, sumxy, sumy, sumyy, xl, yl
  Integer (i4b) :: l
! ..
! .. intrinsic functions ..
! INTRINSIC DBLE,SQRT
! ..

  sumx = 0.D0
  sumy = 0.D0
  sumxx = 0.D0
  sumxy = 0.D0
  sumyy = 0.D0

  Do l = 1, n
    xl = x(l)
    yl = y(l)
    sumx = sumx + xl
    sumy = sumy + yl
    sumxx = sumxx + xl*xl
    sumyy = sumyy + yl*yl
    sumxy = sumxy + xl*yl
  End Do

  xm = sumxy - sumx*sumy/dble(n)
  xm = xm/(sumxx-sumx*sumx/dble(n))
  b = (sumy-xm*sumx)/dble(n)
  sigx = sqrt((sumxx-sumx*sumx/dble(n))/dble(n-1))
  sigy = sqrt((sumyy-sumy*sumy/dble(n))/dble(n-1))
  r = xm*sigx/sigy

  Return

End Subroutine

!-----------------------------------------------------------------------

Subroutine odeint(ystart, nvar, x1, x2, eps, h1, hmin, nok, nbad, derivs, &
  rkqc)

! FROM V10.0 ON: HERE AND IN CORRESPONDING ROUTINES, NMAX CHANGED TO 4
  Use :: nlte_type
  Use :: nlte_dim
  Use :: nlte_var, Only: kmax, kount, dxsav, xp, yp

  Implicit None

! .. parameters ..
  Integer (i4b) :: maxstp, nmax
  Real (dp) :: two, zero, tiny
  Parameter (maxstp=10000, nmax=4, two=2.D0, zero=0.D0, tiny=1.D-30)
! ..
! .. scalar arguments ..
  Real (dp) :: eps, h1, hmin, x1, x2
  Integer (i4b) :: nbad, nok, nvar
! ..
! .. array arguments ..
  Real (dp) :: ystart(nvar)
! ..
! .. subroutine arguments ..
  External :: derivs, rkqc
! ..
! .. local scalars ..
  Real (dp) :: h, hdid, hnext, x, xsav
  Integer (i4b) :: i, nstp
! ..
! .. local arrays ..
  Real (dp) :: dydx(nmax), y(nmax), yscal(nmax)
! ..
! .. intrinsic functions ..
! INTRINSIC ABS,SIGN
! ..

  x = x1
  h = sign(h1, x2-x1)
  nok = 0
  nbad = 0
  kount = 0

  Do i = 1, nvar
    y(i) = ystart(i)
  End Do

  xsav = x - dxsav*two

itloop: Do nstp = 1, maxstp

    Call derivs(x, y, dydx)

    Do i = 1, nvar
      yscal(i) = abs(y(i)) + abs(h*dydx(i)) + tiny
    End Do

    If (kmax>0) Then

      If (abs(x-xsav)>abs(dxsav)) Then
        If (kount<kmax-1) Then
          kount = kount + 1
          xp(kount) = x
          Do i = 1, nvar
            yp(i, kount) = y(i)
          End Do

          xsav = x
        End If
      End If
    End If

    If ((x+h-x2)*(x+h-x1)>zero) h = x2 - x

    Call rkqc(y, dydx, nvar, x, h, eps, yscal, hdid, hnext, derivs)

    If (hdid==h) Then
      nok = nok + 1
    Else
      nbad = nbad + 1
    End If

    If ((x-x2)*(x2-x1)>=zero) Then
      Do i = 1, nvar
        ystart(i) = y(i)
      End Do

      If (kmax/=0) Then
        kount = kount + 1
        xp(kount) = x
        Do i = 1, nvar
          yp(i, kount) = y(i)
        End Do
      End If

      Return
    End If

    If (abs(hnext)<hmin) Then
      Write (999, *) ' STOP: STEPSIZE SMALLER THAN MINIMUM.'
      Stop ' STEPSIZE SMALLER THAN MINIMUM.'
    End If

    h = hnext

  End Do itloop

  Write (999, *) ' STOP: TOO MANY STEPS.'
  Stop 'TOO MANY STEPS.'

End Subroutine

!-----------------------------------------------------------------------

Subroutine rkqc(y, dydx, n, x, htry, eps, yscal, hdid, hnext, derivs)

  Use :: nlte_type
  Use :: nlte_dim
  Implicit None

! .. parameters ..
  Integer (i4b) :: nmax
  Real (dp) :: fcor, one, safety, errcon
  Parameter (nmax=4, fcor=.0666666667D0, one=1.D0, safety=0.9D0, errcon=6.D-4)
! ..
! .. scalar arguments ..
  Real (dp) :: eps, hdid, hnext, htry, x
  Integer (i4b) :: n
! ..
! .. array arguments ..
  Real (dp) :: dydx(n), y(n), yscal(n)
! ..
! .. subroutine arguments ..
  External :: derivs
! ..
! .. local scalars ..
  Real (dp) :: errmax, h, hh, pgrow, pshrnk, xsav
  Integer (i4b) :: i
! ..
! .. local arrays ..
  Real (dp) :: dysav(nmax), ysav(nmax), ytemp(nmax)
! ..
! .. external subroutines ..
  External :: rk4
! ..
! .. intrinsic functions ..
! INTRINSIC ABS,MAX
! ..

  pgrow = -0.20D0
  pshrnk = -0.25D0
  xsav = x

  Do i = 1, n
    ysav(i) = y(i)
    dysav(i) = dydx(i)
  End Do

  h = htry

100 Continue

  hh = 0.5D0*h

  Call rk4(ysav, dysav, n, xsav, hh, ytemp, derivs)

  x = xsav + hh

  Call derivs(x, ytemp, dydx)
  Call rk4(ytemp, dydx, n, x, hh, y, derivs)

  x = xsav + h

  If (x==xsav) Print *, 'STEPSIZE NOT SIGNIFICANT IN RKQC.'

  Call rk4(ysav, dysav, n, xsav, h, ytemp, derivs)
  errmax = 0.D0

  Do i = 1, n
    ytemp(i) = y(i) - ytemp(i)
    errmax = max(errmax, abs(ytemp(i)/yscal(i)))
  End Do

  errmax = errmax/eps

  If (errmax>one) Then
    h = safety*h*(errmax**pshrnk)
    Go To 100
  Else
    hdid = h
    If (errmax>errcon) Then
      hnext = safety*h*(errmax**pgrow)
    Else
      hnext = 4.D0*h
    End If
  End If

  Do i = 1, n
    y(i) = y(i) + ytemp(i)*fcor
  End Do

  Return
End Subroutine

!-----------------------------------------------------------------------

Subroutine rk4(y, dydx, n, x, h, yout, derivs)

  Use :: nlte_type
  Use :: nlte_dim
  Implicit None

! .. parameters ..
  Integer (i4b) :: nmax
  Parameter (nmax=4)
! ..
! .. scalar arguments ..
  Real (dp) :: h, x
  Integer (i4b) :: n
! ..
! .. array arguments ..
  Real (dp) :: dydx(n), y(n), yout(n)
! ..
! .. subroutine arguments ..
  External :: derivs
! ..
! .. local scalars ..
  Real (dp) :: h6, hh, xh
  Integer (i4b) :: i
! ..
! .. local arrays ..
  Real (dp) :: dym(nmax), dyt(nmax), yt(nmax)
! ..

  hh = h*0.5D0
  h6 = h/6.D0
  xh = x + hh

  Do i = 1, n
    yt(i) = y(i) + hh*dydx(i)
  End Do

  Call derivs(xh, yt, dyt)

  Do i = 1, n
    yt(i) = y(i) + hh*dyt(i)
  End Do

  Call derivs(xh, yt, dym)


  Do i = 1, n
    yt(i) = y(i) + h*dym(i)
    dym(i) = dyt(i) + dym(i)
  End Do


  Call derivs(x+h, yt, dyt)


  Do i = 1, n
    yout(i) = y(i) + h6*(dydx(i)+dyt(i)+2.D0*dym(i))
  End Do

  Return

End Subroutine

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

!-----------------------------------------------------------------------

Subroutine sort1(n, ra)

  Use :: nlte_type
  Use :: nlte_dim
  Implicit None


! sorts for absolute values


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
      If (abs(ra(j))<abs(ra(j+1))) j = j + 1
    End If

    If (abs(rra)<abs(ra(j))) Then
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

!-----------------------------------------------------------------------

Subroutine spline(x, y, a, b, c, n)

  Use :: nlte_type
  Use :: nlte_dim
  Implicit None

! -----computation of cubic spline coefficients a,b,c
! -----for approximation of a data array (x,y) such that
! -----f(x)=y(i)+del*(a(i)+del*(b(i)+c(i)*del)) in the
! -----interval x(i)-x-x(i+1), where del=x-x(i).
! -----vanishing curvature at the boundaries is assumed:
! -----f"(x(i))=f"(x(n))=0

! .. scalar arguments ..
  Integer (i4b) :: n
! ..
! .. array arguments ..
  Real (dp) :: a(n), b(n), c(n), x(n), y(n)
! ..
! .. local scalars ..
  Real (dp) :: x1, z1
  Integer (i4b) :: i
! ..
! .. local arrays ..
  Real (dp) :: xx(1000)
! ..

  If (n>1000) Then
    Write (999, *) ' STOP: N IN SPLINE TOO LARGE'
    Stop ' N IN SPLINE TOO LARGE'
  End If

  Do i = 1, n
    xx(i) = x(i)
  End Do

  x1 = x(1)

  Do i = 2, n
    z1 = x(i)
    x(i) = z1 - x1
    x1 = z1
  End Do

  Do i = 2, n - 1
    a(i) = x(i+1)/(x(i)+x(i+1))
    b(i) = 1.D0 - a(i)
    c(i) = 6.D0*((y(i+1)-y(i))/x(i+1)-(y(i)-y(i-1))/x(i))/(x(i)+x(i+1))
  End Do

  b(2) = 2.D0

  Do i = 3, n - 1
    x1 = b(i)/b(i-1)
    b(i) = 2.D0 - x1*a(i-1)
    c(i) = c(i) - x1*c(i-1)
  End Do

  c(n-1) = c(n-1)/b(n-1)

  Do i = n - 2, 2, -1
    c(i) = (c(i)-a(i)*c(i+1))/b(i)
  End Do

  c(1) = 0.D0
  c(n) = 0.D0

  Do i = 1, n - 1
    a(i) = (y(i+1)-y(i))/x(i+1) - x(i+1)*(2.D0*c(i)+c(i+1))/6.D0
    b(i) = .5D0*c(i)
    c(i) = (c(i+1)-c(i))/x(i+1)/6.D0
  End Do

  Do i = 1, n
    x(i) = xx(i)
  End Do

  Return

End Subroutine

!***********************************************************************

!subroutines: algebraic ones

!***********************************************************************

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
! INTRINSIC ABS,MAX
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
! INTRINSIC MAX0
! ..
iloop: Do i = 2, n
    im1 = i - 1

jloop1: Do j = 1, im1
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

    End Do jloop1

    Do j = i, n
      summ = 0.0D+00

      Do l = 1, im1
        summ = summ + a(i, l)*a(l, j)
      End Do
      a(i, j) = a(i, j) - summ
    End Do

  End Do iloop

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

back: Do i = 1, n
jloop2: Do j = 1, n
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
    End Do jloop2

  End Do back

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

! for this program package we assume without proving that n > 1


! .. parameters ..
  Integer (i4b), Parameter :: nd = id_ndept
! ..
! .. scalar arguments ..
  Integer (i4b) :: n
! ..
! .. array arguments ..
  Real (dp) :: a(nd), b(nd), c(nd), q(nd)
! ..
! .. local scalars ..
  Real (dp) :: h
  Integer (i4b) :: i, ni
! ..

  c(1) = c(1)/b(1)
  q(1) = q(1)/b(1)

! if (n.eq.1) return

  Do i = 2, n - 1
    h = b(i) - a(i)*c(i-1)
    c(i) = c(i)/h
    q(i) = (q(i)+q(i-1)*a(i))/h
  End Do

  q(n) = (q(n)+q(n-1)*a(n))/(b(n)-c(n-1)*a(n))

  Do i = 1, n - 1
    ni = n - i
    q(ni) = q(ni) + c(ni)*q(ni+1)
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

  If (n>100) Then
    Write (999, *) ' STOP: TOO MANY GRID POINTS IN DIAG'
    Stop ' TOO MANY GRID POINTS IN DIAG'
  End If

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

Subroutine lubksb(a, n, np, indx, b)

  Use :: nlte_type
  Use :: nlte_dim
  Implicit None

! .. scalar arguments ..
  Integer (i4b) :: n, np
! ..
! .. array arguments ..
  Real (dp) :: a(np, np), b(n)
  Integer (i4b) :: indx(n)
! ..
! .. local scalars ..
  Real (dp) :: summ
  Integer (i4b) :: i, ii, j, ll
! ..

  ii = 0

  Do i = 1, n
    ll = indx(i)
    summ = b(ll)
    b(ll) = b(i)

    If (ii/=0) Then
      Do j = ii, i - 1
        summ = summ - a(i, j)*b(j)
      End Do
    Else If (summ/=0.D0) Then
      ii = i
    End If

    b(i) = summ
  End Do

  Do i = n, 1, -1
    summ = b(i)

    Do j = i + 1, n
      summ = summ - a(i, j)*b(j)
    End Do

    b(i) = summ/a(i, i)
  End Do

  Return

End Subroutine

!-----------------------------------------------------------------------

Subroutine ludcmp(a, n, np, indx, d)

  Use :: nlte_type
  Use :: nlte_dim
  Implicit None

! .. parameters ..
  Integer (i4b) :: nmax
  Parameter (nmax=id_llevs+1)
! ..
! .. scalar arguments ..
  Real (dp) :: d
  Integer (i4b) :: n, np
! ..
! .. array arguments ..
  Real (dp) :: a(np, np)
  Integer (i4b) :: indx(n)
! ..
! .. local scalars ..
  Real (dp) :: aamax, dum, summ
  Integer (i4b) :: i, imax, j, k
! ..
! .. local arrays ..
  Real (dp) :: vv(nmax)
! ..
! .. intrinsic functions ..
! INTRINSIC ABS
! ..

  d = 1.D0

  Do i = 1, n
    aamax = 0.D0
    Do j = 1, n
      If (abs(a(i,j))>aamax) aamax = abs(a(i,j))
    End Do
    vv(i) = 1.D0/aamax
  End Do

jloop: Do j = 1, n

    Do i = 1, j - 1
      summ = a(i, j)

      Do k = 1, i - 1
        summ = summ - a(i, k)*a(k, j)
      End Do

      a(i, j) = summ
    End Do

    aamax = 0.D0

    Do i = j, n
      summ = a(i, j)

      Do k = 1, j - 1
        summ = summ - a(i, k)*a(k, j)
      End Do

      a(i, j) = summ
      dum = vv(i)*abs(summ)

      If (dum>=aamax) Then
        imax = i
        aamax = dum
      End If
    End Do

    If (j/=imax) Then

      Do k = 1, n
        dum = a(imax, k)
        a(imax, k) = a(j, k)
        a(j, k) = dum
      End Do

      d = -d
      vv(imax) = vv(j)
    End If

    indx(j) = imax

    If (a(j,j)==0.) Then
      Write (999, *) ' STOP: MATRIX SINGULAR'
      Stop ' MATRIX SINGULAR'
    End If

    If (j/=n) Then
      dum = 1.D0/a(j, j)
      Do i = j + 1, n
        a(i, j) = a(i, j)*dum
      End Do
    End If

  End Do jloop

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
! .. local scalars ..
  Integer (i4b) :: i
! ..

  Do i = 1, n
    v(i) = a(i)*v(i)
  End Do

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

  a = a + b

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

  a = a - b

  Return

End Subroutine
!ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
!cc                                                                  ccc
!cc      subroutines and functions for new atomic data (partly not used)
!cc      (NOT ADAPTED FOR FAST F90)                                  ccc
!ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc

Function pvr(y)

  Use :: nlte_type
  Use :: nlte_dim
  Implicit None

! generates van regemorter function for neutrals

! .. scalar arguments ..
  Real (dp) :: y, pvr
! ..
! .. local scalars ..
  Real (dp) :: ppl, yyl
  Integer (i4b) :: i, il
! ..
! .. local arrays ..
  Real (dp) :: pl(11), yl(11)
! ..
! .. intrinsic functions ..
! INTRINSIC LOG,LOG10,SQRT
! ..
! .. save statement ..
  Save
! ..
! .. data statements ..
  Data yl/1.0D0, 0.60206D0, 0.30103D0, 0.0D0, -0.39794D0, -0.69897D0, -1.0D0, &
    -1.39794D0, -1.69897D0, -2.0D0, -2.30103D0/
  Data pl/ -1.6383D0, -1.398D0, -1.201D0, -1.0D0, -0.68D0, -0.4802D0, &
    -0.3072D0, -0.1203D0, -0.01954D0, 0.06446D0, 0.11442D0/
! ..
  If (y>1.0D1) Then

    pvr = 0.066D0/sqrt(y)
  Else If (y<=5.0D-3) Then

    pvr = 0.2757D0*(-0.57722D0-log(y))
  Else
    yyl = log10(y)
    Do i = 2, 11
      If (yyl>=yl(i)) Then
        il = i
        Exit
      End If
    End Do

    ppl = pl(il) + (yyl-yl(il))*(pl(il-1)-pl(il))/(yl(il-1)-yl(il))
    pvr = 10.0D0**ppl

  End If

  Return
End Function

!-----------------------------------------------------------------------

Function pftn(t, u, eps)

  Use :: nlte_type
  Use :: nlte_dim
  Implicit None

! basic function for classical collisional excitation

! .. scalar arguments ..
  Real (dp) :: eps, t, u, pftn
! ..
! .. local scalars ..
  Real (dp) :: rat1, rat2, res1, res2
! ..
! .. external functions ..
  Real (dp) :: expino
  External :: expino
! ..
! .. intrinsic functions ..
! INTRINSIC ABS,EXP
! ..
! .. save statement ..
  Save
! ..
  If (abs(eps)<=1.0E-35) Then
    rat1 = u*1.579D5/t
    res1 = expino(rat1)

    pftn = (2.0D0+rat1)*exp(rat1)*res1 - (1.0D0+(1.0D0/rat1))
  Else
    rat1 = (u+abs(eps))*1.579D5/t
    rat2 = u*1.579D5/t
    res1 = expino(rat1)
    res2 = expino(rat2)
    pftn = ((u/abs(eps))+1.0D0)*exp(rat1)*res1 - ((u/abs(eps))-1.0D0)*exp(rat2 &
      )*res2 - 1.0D0/rat2

  End If

  Return
End Function

!-----------------------------------------------------------------------

Subroutine coll21(t, kl, ku, c)

  Use :: nlte_type
  Use :: nlte_dim
  Implicit None

! generates collisional rate coefficients among the 19 states of
! helium with n = 1, 2, 3, and 4, using rates evaluated from the
! cross sections calculated by berrington and kingston (j. phys.b.
! 20, 6631(1987)). collisional rate coefficients have been evaluated
! numerically by p.j.storey from the unpublished computer output
! files of berrington and kingston.

! the states included in the calculation are labelled sequentially
! in order of increasing energy:
! 1          1 sing s
! 2          2 trip s
! 3          2 sing s
! 4          2 trip p
! 5          2 sing p
! 6          3 trip s
! 7          3 sing s
! 8          3 trip p
! 9          3 trip d
! 10          3 sing d
! 11          3 sing p
! 12          4 trip s
! 13          4 sing s
! 14          4 trip p
! 15          4 trip d
! 16          4 sing d
! 17          4 trip f
! 18          4 sing f
! 19          4 sing p

! this ordering differs slightly from that of berrington and
! kingston, in which 15 and 16, and 17 and 18, were interchanged.

! .. scalar arguments ..
  Real (dp) :: c, t
  Integer (i4b) :: kl, ku
! ..
! .. local scalars ..
  Real (dp) :: c1, c2, half, one, temp, tfac, two, xxx, zero
  Integer (i4b) :: i, ir, j, jj, n1, nf, nt, ntm2
! ..
! .. local arrays ..
  Real (dp) :: a(929), b(10), ener(19), s(19), stwt(19)
  Integer (i4b) :: key(4, 4, 2), l(19), n(19), nstart(172)
  Character :: iden(19)*14
! ..
! .. intrinsic functions ..
! INTRINSIC LOG10,SQRT
! ..
! .. save statement ..
  Save
! ..
! .. data statements ..
  Data n/1, 2, 2, 2, 2, 3, 3, 3, 3, 3, 3, 4, 4, 4, 4, 4, 4, 4, 4/
  Data l/0, 0, 0, 1, 1, 0, 0, 1, 2, 2, 1, 0, 0, 1, 2, 2, 3, 3, 1/
  Data s/0, 1, 0, 1, 0, 1, 0, 1, 1, 0, 0, 1, 0, 1, 1, 0, 1, 0, 0/
  Data key/1, 3, 7, 13, 0, 5, 11, 19, 0, 0, 10, 16, 0, 0, 0, 18, 0, 2, 6, 12, &
    0, 4, 8, 14, 0, 0, 9, 5, 0, 0, 0, 17/
  Data ener/0.0D0, 19.8198D0, 20.6160D0, 20.96432D0, 21.2182D0, 22.7187D0, &
    22.9206D0, 23.00731D0, 23.0739D0, 23.0743D0, 23.0873D0, 23.5942D0, &
    23.6738D0, 23.7081D0, 23.7363D0, 23.7366D0, 23.7373D0, 23.7373D0, &
    23.7423D0/
  Data stwt/1.0D0, 3.0D0, 1.0D0, 9.0D0, 3.0D0, 3.0D0, 1.0D0, 9.0D0, 1.5D1, &
    5.0D0, 3.0D0, 3.0D0, 1.0D0, 9.0D0, 1.5D1, 5.0D0, 2.1D1, 7.0D0, 3.0D0/
  Data iden/'HE I 1 SING S', 'HE I 2 TRIP S', 'HE I 2 SING S', &
    'HE I 2 TRIP P', 'HE I 2 SING P', 'HE I 3 TRIP S', 'HE I 3 SING S', &
    'HE I 3 TRIP P', 'HE I 3 TRIP D', 'HE I 3 SING D', 'HE I 3 SING P', &
    'HE I 4 TRIP S', 'HE I 4 SING S', 'HE I 4 TRIP P', 'HE I 4 TRIP D', &
    'HE I 4 SING D', 'HE I 4 TRIP F', 'HE I 4 SING F', 'HE I 4 SING P'/
  Data nstart/1, 6, 11, 16, 20, 28, 32, 40, 44, 52, 57, 62, 67, 72, 77, 82, &
    88, 92, 98, 104, 110, 114, 120, 125, 129, 135, 139, 147, 151, 157, 164, &
    170, 177, 183, 190, 195, 202, 208, 213, 220, 225, 232, 236, 243, 247, 251, &
    260, 266, 273, 278, 285, 290, 300, 304, 309, 316, 320, 324, 329, 333, 338, &
    343, 347, 352, 357, 362, 367, 372, 376, 382, 386, 391, 395, 401, 405, 410, &
    414, 421, 425, 431, 435, 440, 445, 449, 454, 459, 465, 470, 475, 480, 487, &
    491, 497, 503, 508, 515, 520, 525, 530, 536, 542, 547, 552, 559, 564, 571, &
    576, 581, 587, 592, 598, 603, 608, 613, 617, 623, 630, 635, 642, 646, 650, &
    655, 660, 666, 671, 677, 683, 689, 695, 702, 707, 713, 718, 723, 728, 732, &
    737, 741, 745, 750, 754, 759, 765, 771, 777, 782, 789, 796, 801, 805, 810, &
    815, 819, 824, 831, 837, 844, 850, 856, 861, 868, 873, 877, 882, 890, 895, &
    905, 909, 913, 920, 925, 930/

  Data (a(i), i=1, 95)/1.7339D-07, 2.7997D-08, -1.3812D-08, 2.6639D-09, &
    1.7776D-09, 2.9820D-07, 7.5210D-08, -3.5975D-09, 3.2270D-09, 1.5245D-09, &
    1.5601D-05, 1.5340D-06, -2.2122D-06, -1.1073D-07, 1.9249D-07, 2.3682D-08, &
    1.0638D-08, 2.0959D-09, 2.8381D-10, 3.0497D-05, 1.9252D-05, 6.3109D-06, &
    6.9098D-07, -2.8039D-07, -2.1128D-07, -1.2192D-07, -4.4417D-08, &
    1.3896D-06, 1.5715D-07, -8.4358D-08, -2.8800D-08, 5.9599D-08, 3.4756D-08, &
    1.2183D-08, 3.7999D-09, 8.5500D-10, -3.9428D-10, -5.3999D-10, -2.2962D-10, &
    2.2510D-06, 6.1436D-07, -1.2437D-07, -7.1718D-08, 5.9026D-05, 3.8150D-05, &
    1.1426D-05, 9.2886D-07, -6.4827D-07, -4.4270D-07, -1.8611D-07, &
    -5.6403D-08, 5.8752D-06, 2.5167D-06, 2.0787D-07, -2.3353D-07, -8.9900D-08, &
    5.6334D-08, 2.9313D-09, 1.7775D-09, 5.5494D-10, -5.8914D-10, 8.1939D-06, &
    2.4014D-07, 3.8681D-07, 2.8446D-07, -6.1936D-08, 1.8173D-06, -4.7530D-07, &
    -6.6432D-08, 5.4898D-08, -1.2377D-08, 2.0732D-05, 6.5991D-07, 1.5840D-06, &
    4.9920D-07, -1.8332D-07, 3.2273D-06, -4.9880D-07, -2.4929D-07, 1.1964D-07, &
    -2.1996D-08, 9.7096D-08, 1.3557D-08, 7.4404D-09, 1.3858D-09, -1.3778D-09, &
    -8.4885D-10, 3.5068D-06, -5.0675D-07, -1.2252D-07, 6.0514D-08, 7.0524D-06, &
    1.4454D-06, 8.3966D-07, 2.7203D-07/
  Data (a(i), i=96, 190)/ -2.3854D-08, -8.6693D-08, 7.1193D-06, 8.3111D-09, &
    -9.3916D-07, -2.7944D-08, 1.7803D-07, -3.9216D-08, 9.2760D-06, 2.4761D-06, &
    1.0095D-06, 3.4039D-07, -8.6900D-08, -1.1156D-07, 2.5036D-05, -5.7791D-06, &
    -1.8197D-06, 1.0630D-06, 9.8746D-09, 3.1048D-09, 1.2172D-09, 1.8411D-10, &
    -1.5835D-10, -1.4443D-10, 1.9240D-06, 1.5012D-07, 7.9741D-08, 2.9323D-08, &
    -1.2796D-08, 5.0457D-07, -5.5146D-08, -2.7506D-08, 1.4341D-09, 1.3321D-05, &
    2.6960D-06, 2.8059D-07, 1.7548D-08, -2.3623D-07, -1.4619D-07, 1.5294D-06, &
    -9.4925D-08, -1.1347D-07, 8.1980D-09, 1.3324D-04, 7.8068D-05, 2.9238D-05, &
    4.2718D-06, -1.6556D-06, -1.4529D-06, -8.6046D-07, -2.1062D-07, &
    2.8510D-06, -4.7936D-07, -2.4042D-07, 6.2333D-08, 1.3493D-09, 3.5113D-10, &
    -4.7269D-11, 2.4872D-11, -9.7136D-12, -1.4217D-11, 1.5680D-06, 8.9272D-07, &
    2.8313D-07, 5.2456D-08, -2.3404D-08, -2.7304D-08, -8.4539D-09, 1.7967D-07, &
    6.2133D-08, -3.2814D-09, -3.7299D-09, -1.9587D-09, -1.6685D-09, &
    1.2707D-05, 7.5029D-06, 2.8330D-06, 6.7478D-07, -1.5012D-07, -2.3916D-07, &
    -8.9594D-08, 7.6160D-07, 2.0422D-07, -2.5881D-08, -1.7624D-08, &
    -8.3518D-09, -4.1743D-09, 4.6044D-05, 2.2425D-05, 3.3079D-06, 3.6752D-07, &
    -8.4476D-08, -4.8455D-07, -2.5391D-07, 7.0824D-07/
  Data (a(i), i=191, 285)/1.4917D-07, -1.2005D-07, -1.6983D-08, 1.1545D-08, &
    7.2866D-04, 3.8907D-04, 9.6365D-05, 2.3153D-05, 1.8462D-07, -9.0627D-06, &
    -5.4332D-06, 1.0458D-08, 1.6430D-09, 5.6453D-10, 2.0276D-10, -1.6501D-10, &
    -1.0656D-10, 5.2198D-07, 4.2898D-08, -1.0332D-08, -5.6754D-09, &
    -7.1960D-09, 3.0390D-06, 1.5385D-06, 6.2110D-07, 1.4111D-07, -3.7415D-08, &
    -4.8395D-08, -1.7219D-08, 2.7227D-06, 1.4381D-07, -5.5030D-08, &
    -2.2192D-10, -2.8583D-08, 1.8417D-05, 9.3860D-06, 3.9999D-06, 1.0424D-06, &
    -2.1337D-07, -3.3271D-07, -1.2684D-07, 4.6063D-06, -6.7185D-07, &
    -2.5830D-07, 2.9976D-08, 7.8012D-05, 3.4921D-05, 8.3012D-06, 1.5643D-06, &
    -4.3892D-07, -9.5318D-07, -4.6534D-07, 1.7584D-05, -2.8817D-06, &
    -1.0081D-06, 1.3259D-07, 1.8561D-05, 2.6602D-06, -1.7243D-06, -3.5533D-07, &
    2.3315D-08, 1.6239D-08, 7.3739D-09, 1.9570D-09, -2.5667D-10, -5.8101D-10, &
    -2.3947D-10, 4.5821D-12, 4.5024D-11, 3.8796D-07, 1.1239D-07, -2.9033D-08, &
    -5.7237D-09, 2.4186D-09, -2.9670D-09, 1.0887D-06, 4.3737D-07, 8.3753D-08, &
    3.6771D-08, 1.6545D-09, -1.3849D-08, -5.8855D-09, 2.1028D-06, 4.3725D-07, &
    -1.7621D-07, -3.0743D-08, 1.7119D-08, 8.4698D-06, 4.5395D-06, 1.5774D-06, &
    4.1647D-07, -7.9865D-08, -1.6003D-07, -6.0171D-08, 3.1692D-06/
  Data (a(i), i=286, 380)/3.6965D-07, -4.9333D-07, -2.8856D-08, 4.2819D-08, &
    2.0492D-04, 1.3705D-04, 4.0972D-05, 2.9565D-06, -1.4061D-06, -1.8775D-06, &
    -1.6306D-06, -3.6380D-07, 3.1213D-07, 2.0912D-07, 1.1965D-05, 9.9088D-07, &
    -1.3565D-06, -2.3128D-07, 1.1379D-05, 2.6632D-06, -1.6877D-06, &
    -4.0541D-07, 9.8921D-08, 5.9262D-03, 3.2332D-03, 9.2954D-04, 1.6597D-04, &
    -3.7972D-05, -7.1646D-05, -4.0073D-05, 2.5789D-08, 5.2124D-09, &
    -2.2517D-10, -7.9400D-10, 3.3181D-06, 1.9117D-07, 1.0299D-07, -7.7443D-08, &
    2.2953D-06, -1.0879D-06, 3.4368D-07, -1.1230D-07, 1.6826D-08, 6.9774D-06, &
    9.2721D-07, -1.8983D-08, -1.6068D-07, 1.5843D-06, -3.5616D-08, &
    -1.1262D-07, -4.3051D-08, 8.8140D-09, 3.4244D-05, 8.0163D-06, 2.9703D-06, &
    -1.6480D-07, -6.3282D-07, 2.3791D-06, -4.9305D-07, -1.7122D-07, &
    1.2147D-07, 7.2131D-05, 1.7699D-05, 8.4281D-06, -9.3966D-07, -6.8048D-07, &
    5.3785D-05, 6.9745D-06, -2.3554D-06, -7.2141D-07, 3.8501D-07, 4.4794D-06, &
    -6.2435D-07, -2.7340D-07, 2.3127D-08, 4.2314D-08, 3.3784D-06, -5.1356D-07, &
    -2.4321D-07, -6.2698D-09, 3.0812D-08, 5.6419D-08, 1.6889D-08, 2.7037D-09, &
    -1.7433D-09, -9.4507D-10, 1.9714D-06, 4.2743D-08, -7.4114D-08, &
    -2.9304D-08, 2.8973D-06, 1.0079D-06, 2.9561D-07, -6.2977D-09, -4.5782D-08/
  Data (a(i), i=381, 475)/ -2.4331D-08, 4.3284D-06, 2.8398D-07, -3.2705D-07, &
    -1.5079D-07, 5.2686D-06, 1.6539D-06, 3.6079D-07, -1.1589D-07, -5.4904D-08, &
    7.0939D-06, -1.9534D-06, 3.2391D-08, 5.5702D-08, 3.9072D-05, 1.7299D-05, &
    5.1530D-06, -6.0911D-07, -1.2652D-06, -4.6507D-07, 1.1506D-05, &
    -2.0422D-06, -6.1195D-07, 5.8641D-08, 1.0150D-05, -4.6977D-07, &
    -6.9446D-07, 2.2516D-08, 1.1109D-07, 3.6715D-05, 5.6186D-06, -3.2309D-06, &
    -1.6403D-06, 4.1051D-05, 2.3335D-05, 7.8106D-06, 2.0279D-07, -8.1139D-07, &
    -3.7619D-07, -1.4982D-07, 2.5621D-05, -1.0798D-05, 1.4607D-06, 3.2421D-07, &
    6.5478D-09, 2.7233D-09, 3.8646D-10, -1.9143D-10, -1.0483D-10, -4.8664D-11, &
    1.0698D-06, 2.7752D-07, -2.6636D-08, -3.7583D-08, 2.6989D-07, 3.0325D-08, &
    -2.4613D-08, -1.0828D-08, 2.2522D-09, 8.0777D-06, 2.2558D-06, 7.4760D-08, &
    -2.4140D-07, -5.4292D-08, 8.8362D-07, 9.5045D-08, -5.0543D-08, &
    -1.9134D-08, 1.1284D-05, 4.0659D-06, 7.6246D-07, -9.0394D-08, -8.7408D-08, &
    1.2192D-06, -2.0362D-07, -1.0700D-07, 2.1990D-08, 1.2519D-08, 5.2758D-05, &
    2.0175D-05, 3.6690D-06, -4.9014D-07, -7.0474D-07, -3.9931D-07, 4.7638D-05, &
    1.5546D-05, 6.2294D-07, -1.3428D-06, -1.8177D-07, 4.1203D-06, -4.9340D-07, &
    -3.0198D-07, -1.5902D-08, 2.8567D-08, 2.2397D-06/
  Data (a(i), i=476, 570)/ -2.1306D-08, -1.7897D-07, 2.8178D-08, 4.1722D-08, &
    2.0594D-04, 1.2838D-04, 5.8219D-05, 1.1735D-05, -6.9778D-06, -6.2486D-06, &
    -1.3800D-06, 2.9170D-06, -7.0833D-07, -1.2673D-07, 4.5069D-08, 1.0974D-09, &
    4.7908D-10, 3.0294D-11, -7.0815D-11, -1.2039D-11, 5.6357D-12, 8.1274D-07, &
    4.5158D-07, 1.1655D-07, -2.1481D-08, -2.2493D-08, -6.5884D-09, 9.8696D-08, &
    3.6499D-08, 1.4768D-09, -8.0141D-09, -2.8147D-09, 5.8287D-06, 3.2469D-06, &
    9.4927D-07, -7.2550D-08, -1.4231D-07, -5.8408D-08, -1.4468D-08, &
    4.6070D-07, 1.5956D-07, -3.2311D-09, -3.3487D-08, -7.9013D-09, 4.0092D-06, &
    1.2195D-06, 1.1192D-07, -1.3870D-07, -6.2194D-08, 5.4918D-07, -5.1133D-08, &
    -5.4844D-08, -3.1054D-09, 6.1049D-09, 2.4833D-05, 1.1280D-05, 2.2680D-06, &
    -4.4637D-07, -2.8942D-07, -1.2382D-07, 6.8223D-05, 3.6505D-05, 7.8792D-06, &
    -1.7946D-06, -1.4925D-06, -4.6018D-07, 2.5677D-06, -7.2440D-08, &
    -2.2441D-07, -4.1769D-08, 2.3833D-08, 1.1885D-06, 4.0457D-09, -1.2273D-07, &
    -2.3271D-08, 1.5149D-08, 9.1912D-05, 5.6621D-05, 1.4338D-05, -2.6072D-06, &
    -2.3688D-06, -6.2490D-07, -1.4386D-07, 1.0821D-06, -1.6304D-07, &
    -7.9756D-08, -4.9881D-09, 9.9527D-09, 1.5288D-03, 9.8517D-04, 3.4966D-04, &
    3.4238D-05, -4.3610D-05, -3.2734D-05, -8.3676D-06/
  Data (a(i), i=571, 665)/9.6630D-09, 3.2976D-09, 2.1545D-10, -3.3938D-10, &
    -1.4397D-10, 3.4140D-07, 7.6536D-08, -2.2485D-08, -1.8244D-08, &
    -2.0611D-09, 1.3128D-06, 6.4165D-07, 1.4728D-07, -2.7600D-08, -3.0125D-08, &
    -1.0320D-08, 1.6295D-06, 3.3697D-07, -5.7547D-08, -7.3134D-08, &
    -1.5301D-08, 8.0375D-06, 4.0695D-06, 1.1441D-06, -5.8423D-08, -1.6488D-07, &
    -7.5920D-08, 2.1173D-06, -1.9184D-07, -2.4986D-07, 7.2656D-09, 2.7560D-08, &
    7.9040D-06, 3.3963D-06, 4.8082D-07, -3.1619D-07, -1.6501D-07, 5.8017D-06, &
    -5.5173D-07, -6.3613D-07, 1.9633D-08, 7.8681D-08, 9.4568D-06, 9.8227D-08, &
    -7.8983D-07, -2.0343D-07, 7.0351D-05, 3.6017D-05, 8.0504D-06, -1.7276D-06, &
    -1.5329D-06, -4.7799D-07, 3.5263D-05, 1.8639D-05, 4.2070D-06, -4.2643D-07, &
    -5.2726D-07, -2.5481D-07, -1.0144D-07, 4.4613D-06, -8.7802D-07, &
    -3.9633D-07, 6.7953D-08, 4.1539D-08, 1.8780D-04, 1.0547D-04, 2.0983D-05, &
    -3.7076D-06, -3.7090D-06, -1.5934D-06, -2.6146D-07, 1.6384D-05, &
    -3.3765D-06, -7.8390D-07, 1.2382D-07, 1.5748D-05, -9.4352D-07, &
    -4.8936D-07, -4.9967D-07, 3.8177D-10, 9.6781D-11, -1.9622D-11, &
    -1.7859D-11, -4.4781D-12, 3.0310D-07, 1.1193D-07, -5.4059D-09, &
    -1.4249D-08, -1.9549D-09, 4.7664D-08, 8.5684D-09, -5.0948D-09, &
    -3.1313D-09, 4.9557D-10, 4.2702D-10/
  Data (a(i), i=666, 760)/2.3433D-06, 8.3824D-07, 2.6585D-08, -7.5944D-08, &
    -1.5093D-08, 1.8775D-07, 2.7614D-08, -2.1193D-08, -1.0264D-08, 1.9660D-09, &
    1.1637D-09, 7.7890D-06, 3.2901D-06, -2.9753D-07, -5.3153D-07, -5.2098D-08, &
    3.3947D-08, 5.5937D-07, 5.1363D-08, -8.5390D-08, -2.2580D-08, 1.0857D-08, &
    3.5157D-09, 4.1303D-05, 2.1947D-05, 3.9480D-06, -1.4234D-06, -9.0220D-07, &
    -2.4485D-07, 1.1031D-04, 6.6726D-05, 1.9581D-05, -4.3637D-07, -2.8911D-06, &
    -1.4332D-06, -3.2332D-07, 3.0922D-06, 2.0624D-07, -4.0771D-07, &
    -1.0402D-07, 4.8807D-08, 1.2491D-06, 2.1187D-07, -1.5951D-07, -7.2537D-08, &
    1.4035D-08, 9.6514D-09, 2.5146D-05, -9.4490D-08, -3.4103D-06, 2.4020D-07, &
    2.9564D-07, 6.6566D-07, -3.0772D-08, -1.0055D-07, -3.1082D-09, 1.4897D-08, &
    7.7189D-04, 3.3710D-04, 3.9258D-05, -1.5271D-05, -6.0652D-06, 2.1986D-04, &
    -3.3135D-05, -8.5682D-06, -9.2988D-07, 3.8085D-06, -1.3259D-07, &
    -3.8517D-07, -5.7276D-08, 3.7226D-08, 2.0825D-09, 1.6707D-10, -1.2443D-10, &
    -3.9117D-11, 1.0068D-07, 6.2278D-09, -5.9169D-09, -3.9535D-09, 5.3334D-07, &
    2.2022D-07, 1.4522D-08, -2.8440D-08, -8.1406D-09, 6.0187D-07, 5.8318D-08, &
    -2.7772D-08, -2.6639D-08, 3.0929D-06, 1.0648D-06, 1.2317D-07, -9.9233D-08, &
    -4.0796D-08, 1.5217D-06, 8.3095D-08/
  Data (a(i), i=761, 855)/ -1.9819D-07, -6.0417D-08, 2.9553D-08, 9.0905D-09, &
    8.6651D-06, 3.9125D-06, -2.2871D-07, -6.4651D-07, -7.4440D-08, 4.3543D-08, &
    4.4505D-06, 2.7350D-07, -4.6428D-07, -2.2293D-07, 5.5275D-08, 3.6505D-08, &
    8.3471D-06, 5.0215D-07, -1.1189D-06, -2.5404D-07, 1.3175D-07, 1.1687D-04, &
    7.0674D-05, 2.0219D-05, -1.2136D-06, -3.2299D-06, -1.3981D-06, &
    -2.8271D-07, 3.7298D-05, 2.2054D-05, 5.5510D-06, -9.1515D-07, -1.0832D-06, &
    -3.6657D-07, -4.8581D-08, 2.2506D-06, -2.9469D-07, -2.1641D-07, &
    -3.9630D-08, 5.2879D-08, 2.7894D-05, -3.9312D-06, -2.6413D-06, 5.5848D-07, &
    8.7829D-06, -1.4800D-06, -7.0982D-07, -2.1645D-08, 1.3258D-07, 1.2248D-05, &
    -1.4244D-06, -7.6696D-07, -1.5029D-07, 6.8467D-08, 1.9392D-04, &
    -2.6676D-05, -8.0536D-06, -6.3609D-07, 2.3074D-05, 7.2856D-07, &
    -2.3751D-06, -5.8464D-07, 1.3685D-07, 1.9646D-08, 1.2918D-08, 4.6924D-09, &
    5.1593D-10, -4.4715D-10, -3.4091D-10, -9.7254D-11, 3.4117D-07, 1.2731D-07, &
    -2.2300D-08, -2.5297D-08, -1.1877D-10, 2.3436D-09, 9.4463D-07, 4.9464D-07, &
    8.9138D-08, -2.2158D-08, -1.2008D-08, -4.3715D-09, -2.8719D-09, &
    1.4091D-06, 4.8184D-07, -9.6135D-08, -1.0945D-07, -9.3174D-09, 9.0166D-09, &
    4.3747D-06, 2.1053D-06, 1.7702D-07, -3.0956D-07, -1.2902D-07, -8.9418D-09/
  Data (a(i), i=856, 929)/1.4297D-06, 7.7612D-08, -1.7188D-07, 6.4348D-10, &
    1.1572D-08, 7.4828D-06, 4.0177D-06, 1.0890D-06, 8.6393D-08, -5.4701D-08, &
    -5.7656D-08, -3.3030D-08, 5.0379D-06, 1.2861D-07, -7.1313D-07, &
    -3.3703D-09, 7.9367D-08, 6.2235D-06, 9.8196D-07, -4.6416D-07, -2.3278D-07, &
    2.9519D-05, 1.5275D-05, 3.1448D-06, -6.4571D-07, -3.4445D-07, 3.6508D-05, &
    2.0524D-05, 3.1203D-06, -2.3606D-06, -1.5012D-06, -2.5807D-07, 1.3807D-07, &
    9.9656D-08, 3.1628D-06, -5.6361D-08, -3.8861D-07, -1.0385D-09, 2.9930D-08, &
    3.1351D-04, 2.3774D-04, 1.0342D-04, 1.2429D-05, -1.4107D-05, -9.2597D-06, &
    -1.7338D-06, 8.4903D-07, 7.3795D-07, 2.2065D-07, 1.2134D-05, 4.6651D-07, &
    -9.3659D-07, -2.7485D-07, 1.2462D-05, 6.2559D-07, -6.1634D-07, &
    -5.1080D-07, 1.2493D-02, 8.2057D-03, 2.9562D-03, 3.3090D-04, -3.7349D-04, &
    -2.8886D-04, -7.4606D-05, 1.2726D-05, -1.6480D-07, -1.5006D-06, &
    -1.1181D-07, 1.4846D-07, 8.7160D-03, 4.2652D-03, 5.2455D-04, -2.6363D-04, &
    -7.9836D-05/

  Data zero, one, two, half/0.D0, 1.D0, 2.D0, 0.5D0/
! ..
  c1 = half*log10(5.0D7)
  c2 = half*log10(5.0D1)
  temp = t
! if(t.gt.30000.) temp = 30000.
! if(t.lt.2000.) temp = 2000.

! evaluate gammas for requested levels

  xxx = two*(log10(temp)-c1)/c2
  tfac = one/sqrt(temp)
! mapping in cheby array
  j = ((ku*ku-3*ku+4)/2) + kl - 1
! location of first cheby coef
  n1 = nstart(j)
! location of last cheby coef
  nf = nstart(j+1) - 1
  nt = nf - n1 + 1
  ntm2 = nt - 2

! clenshaw summation

  b(nt) = a(nf)
  b(nt-1) = xxx*b(nt) + a(nf-1)
  ir = ntm2
  jj = nf - 2
  Do j = 1, ntm2
    b(ir) = xxx*b(ir+1) - b(ir+2) + a(jj)
    ir = ir - 1
    jj = jj - 1
  End Do
  c = (b(1)-b(3))*tfac*stwt(ku)/stwt(kl)

  Return
End Subroutine

!-----------------------------------------------------------------------

Subroutine pix11(zed, aw, n, lp, nfr, fq, anl)

! pix11 (formula RBF 16) checked and inserted on 30.01.01
! ---- from detail_v201
! ---- modified by e. santolaya (18-11-94) so that only one
! ---- frequency is calculated each time
! ---- fq in hz

  Use :: nlte_type
  Implicit None

! .. scalar arguments ..
  Real (dp) :: aw, zed
  Integer (i4b) :: lp, n, nfr
! ..
! .. array arguments ..
  Real (dp), Dimension (nfr) :: anl, fq
! ..
! .. local scalars ..
  Real (dp) :: a, ai, al, alfac, clight, con1, con2, con3, e, ee, fac, fll, &
    flm, fn, four, frfac, fth, fthhz, g11, g12, g21, g22, g31, g32, gn0, gn1e, &
    gne, one, p, p1, p2, ryd, se, sl, sl4, sm, sm4, sn, sn4, sum, t1, t2, ten, &
    term, two, x, zero
  Integer (i4b) :: i, if, il, ilmax, j, jf, k, l, ll, llk, lll, llm, lm, &
    lmax1, m, mulp, mulr, muls
  Logical :: done
! ..
! .. local arrays ..
  Real (dp), Dimension (1) :: ab, freq
  Real (dp), Dimension (1, 2) :: g2, g3
  Real (dp), Dimension (1500) :: alo, fal

  Integer (i4b), Dimension (1) :: mm
! ..
! .. external functions ..
  Real (dp) :: exp1
  External :: exp1
! ..
! .. intrinsic functions ..
! INTRINSIC ABS, ATAN, DBLE, LOG10, MIN, SQRT
! ..
! .. save statement ..
  Save
! ..
! .. data statements ..
  Data done/.False./
! ..
  If (.Not. done) Then

!   double precision constants

    zero = 0.0D0
    one = 1.0D0
    two = 2.0D0
    four = 4.0D0
    ten = 1.0D1
    con1 = 8.5594D-19
    clight = 2.997925D10


!   evaluate log factorials

    fal(1) = zero
    Do i = 2, 1500
      ai = dble(i-1)
      alo(i-1) = log10(ai)
      fal(i) = alo(i-1) + fal(i-1)
    End Do

    done = .True.

  End If
! rydberg constant

  ryd = 109737.312D0/(one+one/(1822.84D0*aw))

! frequency scaling

  frfac = one/(zed**2*ryd*clight)

! calculation of hydrogenic photoionization cross section


! initialization for n

  fn = dble(n)
  sn = fn*fn
  sn4 = four*sn
  con2 = con1*(fn/zed)**2
  fth = one/sn
  fthhz = fth/frfac*10.D0
  gn0 = 2.3052328943D0 - 2.302585093D0*fal(n+n) - fn*0.61370563888D0 + &
    alo(n)*(fn+one)*2.30258093D0
  lmax1 = min(lp, n-1)
  ilmax = n - lp
  If (lp>100) Then
    lmax1 = n - 1
    ilmax = n
  End If

! reverse frequencies and keep the range down to avoid numerical
! problems

  jf = nfr
  Do if = 1, nfr
    anl(if) = 0.D0
    If (fq(if)<fthhz) Then
      freq(jf) = fq(if)
    Else
      freq(jf) = fthhz
    End If
    jf = jf - 1
  End Do

! initialize g's

  Do i = 1, nfr
    mm(i) = 0
    Do j = 1, 2
      g2(i, j) = zero
      g3(i, j) = zero
    End Do
  End Do

! l loop

lloop: Do il = 1, ilmax
    l = n - il
    m = 0
    al = dble(l)
    k = n - l - 2
    con3 = con2/(two*al+one)

!   frequency loop (freq units after multiplication by frfac
!   are ryd*zed**2)

    jf = nfr
ifloop: Do if = 1, nfr
      If ((frfac*freq(if))<fth) Then
        If (l<=lmax1) anl(jf) = anl(jf) + zero
        jf = jf - 1
        Cycle ifloop
      End If
      g11 = zero
      m = m + 1
      se = (frfac*freq(if)) - fth
      e = sqrt(se)
      x = one + sn*se
      If (k<0) Then

        If (e>=0.314E0) Then
          ee = 6.2831853D0/e
          p = 0.5D0*log10(exp1(-ee))
        Else
          p = zero
        End If

        If (e>1.0D-6) Then
          a = two*(fn-atan(fn*e)/e)
        Else
          a = zero
        End If

        ab(m) = (gn0+a)/2.302585D0 - p - (fn+two)*log10(x)
        gne = 0.1D0
        gn1e = x*gne/(fn+fn)
        g3(m, 2) = gne
        g3(m, 1) = gn1e
        g2(m, 2) = gne*fn*x*(fn+fn-one)
        g2(m, 1) = gn1e*(fn+fn-one)*(four+(fn-one)*x)
      End If
      g22 = g2(m, 2)
      g32 = g3(m, 2)
      g21 = g2(m, 1)
      g31 = g3(m, 1)
      muls = mm(m)
      If (k<0) Then
        Go To 140
      Else If (k>0) Then
        Go To 100
      Else
        If (k/=0) Then
          Write (999, *) ' STOP: K NE 0 IN SUBR. PIX11'
          Stop ' K NE 0 IN SUBR. PIX11'
        End If
        Go To 150
      End If

!     l.lt.n-2

100   Continue
      If (k>1) Go To 110
      ll = n - 1
      lm = n - 2
110   Continue
      sl = dble(ll*ll)
      sl4 = four*sl
      fll = dble(ll)
      g12 = (sn4-sl4+(two*sl-fll)*x)*g22 - sn4*(sn-sl)*(one+(fll+one)**2*se)* &
        g32
      If (l==0) Go To 120
      sm = dble(lm*lm)
      sm4 = four*sm
      flm = dble(lm)
      g11 = (sn4-sm4+(two*sm+flm)*x)*g21 - sn4*(sn-(flm+one)**2)*(one+sm*se)* &
        g31
      g31 = g21
      g21 = g11
120   Continue
      g32 = g22
      g22 = g12
      If (if/=nfr) Go To 130
      ll = ll - 1
      If (l==0) Go To 130
      lm = lm - 1
130   Continue
      If (g12<1.D20) Go To 160
      muls = muls + 35
      g22 = g22*1.D-35
      g32 = g32*1.D-35
      g12 = g12*1.D-35
      If (l==0) Go To 160
      g11 = g11*1.D-35
      g21 = g21*1.D-35
      g31 = g31*1.D-35
      Go To 160

!     l.eq.n-1

140   Continue
      g11 = g31
      If (l==0) g11 = zero
      g12 = g32
      Go To 160

!     l.eq.n-2

150   Continue
      g11 = g21
      If (l==0) g11 = zero
      g12 = g22
160   Continue
      mm(m) = muls
      g2(m, 2) = g22
      g3(m, 2) = g32
      g2(m, 1) = g21
      g3(m, 1) = g31
      alfac = fal(n+l+1) - fal(n-l) + two*(al-fn)*alo(2*n)
      p1 = one
      lll = l + 1
      llm = l - 1
      mulr = 0
      If (llm<1) Go To 170
      Do i = 1, llm
        ai = dble(i)
        p1 = p1*(one+ai*ai*se)
        If (p1>=1.D20) Then
          p1 = p1*1.D-10
          mulr = mulr + 10
        End If
      End Do
170   Continue
      p2 = p1
      llk = llm + 1
      If (llk<1) llk = 1
      Do i = llk, lll
        ai = dble(i)
        p2 = p2*(one+ai*ai*se)
      End Do
      mulp = 0
180   Continue
      If (g12<one) Go To 190
      mulp = mulp - 10
      g12 = g12*1.D-10
      If (l==0) Go To 180
      g11 = g11*1.D-10

      Go To 180
190   Continue
      sum = alfac + dble(mulr) + two*(ab(m)+dble(muls-mulp+1))
      fac = zero
      If (abs(sum)<50.D0) fac = ten**sum
      If (l==0) Go To 200
      g11 = g11*p1*g11
      t1 = fac*g11*x*con3
      If (t1==zero) Go To 210
200   Continue
      g12 = g12*p2*g12
      t2 = fac*g12*x*con3
      If (t2==zero) Go To 210

      If (l<=lmax1) Then
        term = t1*al + t2*(al+one)

        If (lp<100) Then
          anl(jf) = term
        Else
          fac = (2.D0*al+one)/sn
          anl(jf) = anl(jf) + fac*term
        End If
      End If
      jf = jf - 1
    End Do ifloop
  End Do lloop

  Do if = 1, nfr

    If (fq(if)>=fthhz) Then
      x = fthhz/fq(if)
      x = x**3
      anl(if) = anl(if)*x
    End If

  End Do

  Return

210 Continue
  Write (6, Fmt=220) n, l, if

  Write (999, *) ' STOP!'
  Stop
220 Format (1X, 'EXPONENT UNDERFLOW FOR N =', I4, '  L =', I3, '  FREQUEN', &
    'CY NUMBER =', I4)
End Subroutine

!-----------------------------------------------------------------------

Subroutine pix21(is, il, n, xnu, alpha)

! subroutine used only for formula 17 (rbf).
! cross section is placed in alpha.
! xnu is the frecuency, in hertzs

  Use :: nlte_type
  Implicit None


! .. scalar arguments ..
  Real (dp) :: alpha, xnu
  Integer (i4b) :: il, is, n
! ..
! .. local scalars ..
  Real (dp) :: fl, p, x
  Integer (i4b) :: i, ill, iss, j, k, nsl0
! ..
! .. local arrays ..
  Real (dp), Dimension (53) :: a, b, fl0, xfitm
  Real (dp), Dimension (4, 53) :: coef

  Integer (i4b), Dimension (3, 2) :: ist, n0
! ..
! .. intrinsic functions ..

! INTRINSIC LOG10
! ..
! .. data statements ..
  Data ist/1, 36, 20, 11, 45, 28/
  Data n0/1, 2, 3, 2, 2, 3/
  Data fl0/15.7742D0, 14.9831D0, 14.6054D0, 14.3443D0, 14.1440D0, 13.9814D0, &
    13.8445D0, 13.7262D0, 13.6223D0, 13.5294D0, 15.0618D0, 14.6550D0, &
    14.3806D0, 14.1726D0, 14.0050D0, 13.8646D0, 13.7438D0, 13.6378D0, &
    13.5433D0, 14.5634D0, 14.3134D0, 14.1195D0, 13.9611D0, 13.8272D0, &
    13.7111D0, 13.6088D0, 13.5173D0, 14.5636D0, 14.3135D0, 14.1196D0, &
    13.9612D0, 13.8273D0, 13.7112D0, 13.6089D0, 13.5174D0, 14.9110D0, &
    14.5597D0, 14.3105D0, 14.1171D0, 13.9591D0, 13.8254D0, 13.7096D0, &
    13.6075D0, 13.5172D0, 14.9426D0, 14.5822D0, 14.3277D0, 14.1310D0, &
    13.9707D0, 13.8354D0, 13.7184D0, 13.6152D0, 13.5231D0/
  Data xfitm/3.262D-01, 6.135D-01, 9.233D-01, 8.438D-01, 1.020D+00, 1.169D+00, &
    1.298D+00, 1.411D+00, 1.512D+00, 1.602D+00, 7.228D-01, 1.076D+00, &
    1.206D+00, 1.404D+00, 1.481D+00, 1.464D+00, 1.581D+00, 1.685D+00, &
    1.777D+00, 9.586D-01, 1.187D+00, 1.371D+00, 1.524D+00, 1.740D+00, &
    1.854D+00, 1.955D+00, 2.046D+00, 9.585D-01, 1.041D+00, 1.371D+00, &
    1.608D+00, 1.739D+00, 1.768D+00, 1.869D+00, 1.803D+00, 7.360D-01, &
    1.041D+00, 1.272D+00, 1.457D+00, 1.611D+00, 1.741D+00, 1.855D+00, &
    1.870D+00, 1.804D+00, 9.302D-01, 1.144D+00, 1.028D+00, 1.210D+00, &
    1.362D+00, 1.646D+00, 1.761D+00, 1.863D+00, 1.954D+00/
  Data a/6.95319D-01, 1.13101D+00, 1.36313D+00, 1.51684D+00, 1.64767D+00, &
    1.75643D+00, 1.84458D+00, 1.87243D+00, 1.85628D+00, 1.90889D+00, &
    9.01802D-01, 1.25389D+00, 1.39033D+00, 1.55226D+00, 1.60658D+00, &
    1.65930D+00, 1.68855D+00, 1.62477D+00, 1.66726D+00, 1.83599D+00, &
    2.50403D+00, 3.08564D+00, 3.56545D+00, 4.25922D+00, 4.61346D+00, &
    4.91417D+00, 5.19211D+00, 1.74181D+00, 2.25756D+00, 2.95625D+00, &
    3.65899D+00, 4.04397D+00, 4.13410D+00, 4.43538D+00, 4.19583D+00, &
    1.79027D+00, 2.23543D+00, 2.63942D+00, 3.02461D+00, 3.35018D+00, &
    3.62067D+00, 3.85218D+00, 3.76689D+00, 3.49318D+00, 1.16294D+00, &
    1.86467D+00, 2.02110D+00, 2.24231D+00, 2.44240D+00, 2.76594D+00, &
    2.93230D+00, 3.08109D+00, 3.21069D+00/
  Data b/ -1.29300D+00, -2.15771D+00, -2.13263D+00, -2.10272D+00, &
    -2.10861D+00, -2.11507D+00, -2.11710D+00, -2.08531D+00, -2.03296D+00, &
    -2.03441D+00, -1.85905D+00, -2.04057D+00, -2.02189D+00, -2.05930D+00, &
    -2.03403D+00, -2.02071D+00, -1.99956D+00, -1.92851D+00, -1.92905D+00, &
    -4.58608D+00, -4.40022D+00, -4.39154D+00, -4.39676D+00, -4.57631D+00, &
    -4.57120D+00, -4.56188D+00, -4.55915D+00, -4.41218D+00, -4.12940D+00, &
    -4.24401D+00, -4.40783D+00, -4.39930D+00, -4.25981D+00, -4.26804D+00, &
    -4.00419D+00, -4.47251D+00, -3.87960D+00, -3.71668D+00, -3.68461D+00, &
    -3.67173D+00, -3.65991D+00, -3.64968D+00, -3.48666D+00, -3.23985D+00, &
    -2.95758D+00, -3.07110D+00, -2.87157D+00, -2.83137D+00, -2.82132D+00, &
    -2.91084D+00, -2.91159D+00, -2.91336D+00, -2.91296D+00/
  Data ((coef(i,j),i=1,4), j=1, 10)/8.734D-01, -1.545D+00, -1.093D+00, &
    5.918D-01, 9.771D-01, -1.567D+00, -4.739D-01, -1.302D-01, 1.174D+00, &
    -1.638D+00, -2.831D-01, -3.281D-02, 1.324D+00, -1.692D+00, -2.916D-01, &
    9.027D-02, 1.445D+00, -1.761D+00, -1.902D-01, 4.401D-02, 1.546D+00, &
    -1.817D+00, -1.278D-01, 2.293D-02, 1.635D+00, -1.864D+00, -8.252D-02, &
    9.854D-03, 1.712D+00, -1.903D+00, -5.206D-02, 2.892D-03, 1.782D+00, &
    -1.936D+00, -2.952D-02, -1.405D-03, 1.845D+00, -1.964D+00, -1.152D-02, &
    -4.487D-03/
  Data ((coef(i,j),i=1,4), j=11, 19)/7.377D-01, -9.327D-01, -1.466D+00, &
    6.891D-01, 9.031D-01, -1.157D+00, -7.151D-01, 1.832D-01, 1.031D+00, &
    -1.313D+00, -4.517D-01, 9.207D-02, 1.135D+00, -1.441D+00, -2.724D-01, &
    3.105D-02, 1.225D+00, -1.536D+00, -1.725D-01, 7.191D-03, 1.302D+00, &
    -1.602D+00, -1.300D-01, 7.345D-03, 1.372D+00, -1.664D+00, -8.204D-02, &
    -1.643D-03, 1.434D+00, -1.715D+00, -4.646D-02, -7.456D-03, 1.491D+00, &
    -1.760D+00, -1.838D-02, -1.152D-02/
  Data ((coef(i,j),i=1,4), j=20, 27)/1.258D+00, -3.442D+00, -4.731D-01, &
    -9.522D-02, 1.553D+00, -2.781D+00, -6.841D-01, -4.083D-03, 1.727D+00, &
    -2.494D+00, -5.785D-01, -6.015D-02, 1.853D+00, -2.347D+00, -4.611D-01, &
    -9.615D-02, 1.955D+00, -2.273D+00, -3.457D-01, -1.245D-01, 2.041D+00, &
    -2.226D+00, -2.669D-01, -1.344D-01, 2.115D+00, -2.200D+00, -1.999D-01, &
    -1.410D-01, 2.182D+00, -2.188D+00, -1.405D-01, -1.460D-01/
  Data ((coef(i,j),i=1,4), j=28, 35)/1.267D+00, -3.417D+00, -5.038D-01, &
    -1.797D-02, 1.565D+00, -2.781D+00, -6.497D-01, -5.979D-03, 1.741D+00, &
    -2.479D+00, -6.099D-01, -2.227D-02, 1.870D+00, -2.336D+00, -4.899D-01, &
    -6.616D-02, 1.973D+00, -2.253D+00, -3.972D-01, -8.729D-02, 2.061D+00, &
    -2.212D+00, -3.072D-01, -1.060D-01, 2.137D+00, -2.189D+00, -2.352D-01, &
    -1.171D-01, 2.205D+00, -2.186D+00, -1.621D-01, -1.296D-01/
  Data ((coef(i,j),i=1,4), j=36, 44)/1.129D+00, -3.149D+00, -1.910D-01, &
    -5.244D-01, 1.431D+00, -2.511D+00, -3.710D-01, -1.933D-01, 1.620D+00, &
    -2.303D+00, -3.045D-01, -1.391D-01, 1.763D+00, -2.235D+00, -1.829D-01, &
    -1.491D-01, 1.879D+00, -2.215D+00, -9.003D-02, -1.537D-01, 1.978D+00, &
    -2.213D+00, -2.066D-02, -1.541D-01, 2.064D+00, -2.220D+00, 3.258D-02, &
    -1.527D-01, 2.140D+00, -2.225D+00, 6.311D-02, -1.455D-01, 2.208D+00, &
    -2.229D+00, 7.977D-02, -1.357D-01/

  Data ((coef(i,j),i=1,4), j=45, 53)/1.204D+00, -2.809D+00, -3.094D-01, &
    1.100D-01, 1.455D+00, -2.254D+00, -4.795D-01, 6.872D-02, 1.619D+00, &
    -2.109D+00, -3.357D-01, -2.532D-02, 1.747D+00, -2.065D+00, -2.317D-01, &
    -5.224D-02, 1.853D+00, -2.058D+00, -1.517D-01, -6.647D-02, 1.943D+00, &
    -2.055D+00, -1.158D-01, -6.081D-02, 2.023D+00, -2.070D+00, -6.470D-02, &
    -6.800D-02, 2.095D+00, -2.088D+00, -2.357D-02, -7.250D-02, 2.160D+00, &
    -2.107D+00, 1.065D-02, -7.542D-02/
! ..

! check input parameters

  If (is/=1 .And. is/=3) Then
    Write (*, Fmt=100) is
    Write (999, *) ' STOP!'
    Stop
  End If

  If (il<0 .Or. il>2) Then
    Write (*, Fmt=110) il
    Write (999, *) ' STOP!'
    Stop
  End If

  If (n>10) Then
    Write (*, Fmt='(A)') 'N TOO LARGE IN PIX21. GREATER THAN 10'
    Write (999, *) ' STOP!'
    Stop
  End If

! selected beginning and end of coefficientes

  iss = (is+1)/2
  ill = il + 1
  nsl0 = n0(ill, iss)

  If (n<nsl0) Then
    Write (*, Fmt='(A)') 'N TOO SMALL IN PIX21 '
    Write (999, *) ' STOP!'
    Stop
  End If

  i = ist(ill, iss) + n - nsl0

! calculation of cross section for frec. xnu

  fl = log10(xnu)

  x = fl - fl0(i)
  If (x>=-0.001D0) Then

!   allows freq slightly below threshold

    If (x<xfitm(i)) Then
      p = coef(4, i)
      Do k = 1, 3
        p = x*p + coef(4-k, i)
      End Do
      alpha = 10.D0**p*1.D-18

    Else
      alpha = 10.D0**(a(i)+b(i)*x)*1.D-18
    End If

  Else
    alpha = 0.D0
  End If

  Return
100 Format ('  MULTIPLICITY S=', I3, '  NOT ALLOWED FOR HE I')
110 Format ('   ANGULAR MOMENTUM L=', I3, '   NOT ALLOWED FOR HE I')
End Subroutine

!-----------------------------------------------------------------------

Subroutine hecion(t, yne, n, c)

  Use :: nlte_type
  Use :: nlte_dim
  Implicit None

! .. scalar arguments ..
  Real (dp) :: c, t, yne
  Integer (i4b) :: n
! ..
! .. local scalars ..
  Real (dp) :: gam, tfac, u, ui, ul, y
  Integer (i4b) :: i
! ..
! .. local arrays ..
  Real (dp) :: a(6), alpha(73), beta(3), e(73)
! ..
! .. intrinsic functions ..
! INTRINSIC EXP,LOG,LOG10,SQRT
! ..
! .. save statement ..
  Save
! ..
! .. data statements ..

! data for the helium i atom

! ionization energies in cm-1


! threshold photoionization cross sections *10**18


! data for ground state ionization

  Data (e(i), i=1, 73)/198310.76D0, 38454.69D0, 32083.21D0, 29223.69D0, &
    27175.76D0, 15073.87D0, 13445.82D0, 12746.07D0, 12209.11D0, 12205.70D0, &
    12101.29D0, 8012.55D0, 7370.43D0, 7093.62D0, 6866.17D0, 6864.20D0, &
    6858.78D0, 6858.77D0, 6817.94D0, 4963.67D0, 4647.13D0, 4509.95D0, &
    4393.51D0, 4392.37D0, 4389.58D0, 4389.57D0, 4389.03D0, 4368.19D0, &
    3374.53D0, 3195.76D0, 3117.85D0, 3050.59D0, 3049.90D0, 3048.27D0, &
    3048.26D0, 3047.89D0, 3035.72D0, 2442.41D0, 2331.72D0, 2283.36D0, &
    2241.03D0, 2240.54D0, 2239.49D0, 2239.48D0, 2239.24D0, 2231.52D0, &
    1849.34D0, 1775.88D0, 1743.94D0, 1715.58D0, 1715.22D0, 1714.59D0, &
    1714.59D0, 1714.41D0, 1709.25D0, 1448.72D0, 1397.78D0, 1375.34D0, &
    1355.48D0, 1355.24D0, 1354.72D0, 1354.72D0, 1354.60D0, 1350.97D0, &
    1165.48D0, 1128.59D0, 1112.42D0, 1097.88D0, 1097.69D0, 1097.33D0, &
    1097.33D0, 1097.22D0, 1094.52D0/
  Data (alpha(i), i=1, 73)/7.5D0, 5.5D0, 9.5D0, 16.0D0, 13.5D0, 8.0D0, 14.9D0, &
    28.5D0, 18.5D0, 18.1D0, 27.0D0, 10.7D0, 21.1D0, 41.6D0, 36.7D0, 35.7D0, &
    2*19.1D0, 41.7D0, 13.6D0, 27.9D0, 55.8D0, 55.1D0, 55.3D0, 2*41.0D0, &
    17.5D0, 53.9D0, 16.8D0, 35.2D0, 71.3D0, 74.1D0, 71.3D0, 2*71.7D0, 26.7D0, &
    75.7D0, 20.0D0, 43.2D0, 87.7D0, 94.0D0, 90.2D0, 2*86.0D0, 35.6D0, 95.1D0, &
    23.6D0, 51.5D0, 105.4D0, 115.1D0, 110.0D0, 2*109.0D0, 44.4D0, 115.9D0, &
    27.2D0, 60.5D0, 124.5D0, 137.1D0, 130.3D0, 2*132.6D0, 52.9D0, 138.0D0, &
    31.0D0, 70.0D0, 144.5D0, 160.3D0, 152.1D0, 2*156.9D0, 59.9D0, 161.4D0/
  Data a/1.4999D-8, 5.6656D-10, -6.0821D-9, -3.589D-9, 1.5529D-9, 1.3207D-9/
  Data gam/3.1373D-8/
  Data beta/4.7893D-8, -7.7359D-7, 3.7366D-6/
! ..

! ground state - fits to experimental data

  If (n==1) Then
!   =kt/eion
    u = .69501D0*t/e(1)
    ui = 1.0D0/u
    If (u<=1.0D0) Then
      ul = log10(u)

      c = (a(1)+ul*(a(2)+ul*(a(3)+ul*(a(4)+ul*(a(5)+ul*a(6))))))*sqrt(u)*exp( &
        -ui)
    Else
      c = sqrt(u)*(gam*log(u)+beta(1)+ui*(beta(2)+ui*beta(3)))

    End If

    c = yne*c
  Else

!   excited states

    tfac = 1.0D0/sqrt(t)
    y = 1.4388D0*e(n)/t
    c = yne*1.1D-6*alpha(n)*tfac*exp(-y)/y
  End If


  Return
End Subroutine

!-----------------------------------------------------------------------

Function exp1(x)

  Use :: nlte_type
  Use :: nlte_dim
  Implicit None

! .. scalar arguments ..
  Real (dp) :: x, exp1
! ..
! .. local scalars ..
  Real (dp) :: dx
! ..
! .. intrinsic functions ..
! INTRINSIC ABS,EXP
! ..
! .. save statement ..
  Save
! ..
  dx = abs(x)
  If (dx<1.0D-9) Go To 100
  If (dx<1.0D-5) Go To 110
  If (dx<1.0D-3) Go To 120
  exp1 = 1.0D0 - exp(x)

  Return
100 Continue
  exp1 = -x

  Return
110 Continue
  exp1 = ((-x*0.D0)-1.0D0)*x

  Return
120 Continue
  exp1 = (((-x*0.1666666667D0)-0.5D0)*x-1.0D0)*x

  Return
End Function


!-----------------------------------------------------------------------

Function expino(x)

  Use :: nlte_type
  Use :: nlte_dim
  Implicit None

! exponential integral for positive arguments after cody and
! thacher, math. of comp.,22,641(1968)
! .. scalar arguments ..
  Real (dp) :: x, expino
! ..
! .. local scalars ..
  Real (dp) :: a0, a1, a2, a3, a4, a5, b0, b1, b2, b3, b4, c0, c1, c2, c3, c4, &
    c5, c6, d1, d2, d3, d4, d5, d6, e0, e1, e2, e3, e4, e5, e6, ex, ex1, f1, &
    f2, f3, f4, f5, f6, x1
! ..
! .. intrinsic functions ..
! INTRINSIC EXP,LOG
! ..
! .. save statement ..
  Save
! ..
! .. data statements ..
  Data x1/ -1.D20/
  Data a0, a1, a2, a3, a4, a5, b0, b1, b2, b3, b4/ -44178.5471728217D0, &
    57721.7247139444D0, 9938.31388962037D0, 1842.11088668000D0, &
    101.093806161906D0, 5.03416184097568D0, 76537.3323337614D0, &
    32597.1881290275D0, 6106.10794245759D0, 635.419418378382D0, &
    37.2298352833327D0/
  Data c0, c1, c2, c3, c4, c5, c6, d1, d2, d3, d4, d5, d6/4.65627107975096D-7, &
    .999979577051595D0, 9.04161556946329D0, 24.3784088791317D0, &
    23.0192559391333D0, 6.90522522784444D0, .430967839469389D0, &
    10.0411643829054D0, 32.4264210695138D0, 41.2807841891424D0, &
    20.4494785013794D0, 3.31909213593302D0, .103400130404874D0/
  Data e0, e1, e2, e3, e4, e5, e6, f1, f2, f3, f4, f5, &
    f6/ -.999999999998447D0, -26.6271060431811D0, -241.055827097015D0, &
    -895.927957772937D0, -1298.85688746484D0, -545.374158883133D0, &
    -5.66575206533869D0, 28.6271060422192D0, 292.310039388533D0, &
    1332.78537748257D0, 2777.61949509163D0, 2404.01713225909D0, &
    631.657483280800D0/
! ..
  If (x==x1) Go To 130
  ex = exp(-x)
  x1 = x
  If (x>4.E0) Go To 100
  If (x>1.E0) Go To 110
  If (x>0.E0) Go To 120
  ex1 = 0.D0

  Go To 130
100 Continue
  ex1 = (ex+ex*(e0+(e1+(e2+(e3+(e4+(e5+e6/x)/x)/x)/x)/x)/x)/(x+f1+(f2+(f3+(f4+ &
    (f5+f6/x)/x)/x)/x)/x))/x

  Go To 130
110 Continue
  ex1 = ex*(c6+(c5+(c4+(c3+(c2+(c1+c0*x)*x)*x)*x)*x)*x)/(d6+(d5+(d4+ &
    (d3+(d2+(d1+x)*x)*x)*x)*x)*x)

  Go To 130
120 Continue
  ex1 = (a0+(a1+(a2+(a3+(a4+a5*x)*x)*x)*x)*x)/(b0+(b1+(b2+ &
    (b3+(b4+x)*x)*x)*x)*x) - log(x)
130 Continue
  expino = ex1

  Return
End Function

!***********************************************************************

!subroutines: temperature correction

!***********************************************************************

Subroutine tempcorr_start

  Use :: tcorr_var, Only: qffh, qffc, qrbfi, qrbfr, qcbbd, qcbbu, dtqffh, &
    dtqffc, dtqcbfh, dtqcbfc, dtqcbbd, dtqcbbu, dtqrbfr, dtqrbfi, qcbfr, qcbfi

! Energies
  qffh = 0.D0
  qffc = 0.D0
  qcbbd = 0.D0
  qcbbu = 0.D0
  qrbfr = 0.D0
  qrbfi = 0.D0
  qcbfr = 0.D0
  qcbfi = 0.D0
! Derivatives
  dtqffh = 0.D0
  dtqffc = 0.D0
  dtqcbfh = 0.D0
  dtqcbfc = 0.D0
  dtqcbbu = 0.D0
  dtqcbbd = 0.D0
  dtqrbfr = 0.D0
  dtqrbfi = 0.D0

  Return
End Subroutine

!-------------------------------------------------------------------------

Subroutine tempcorr_main(t, xne, taur, dtfcorr, teff, ncor, emaxtc, &
  dtmean, expansion, outthb)

! no update for clumping required, since all rates refer to clumps (incl.
! OPAFF)
! warning: update the algorithm if adiabatic EXPANSION is accounted for

  Use :: princesa_var
  Use :: nlte_type
  Use :: nlte_dim
  Use :: fund_const
  Use :: fastwind_params, Only: tmin1, tminabs

  Use :: nlte_var, Only: modnam, xnelte, optcmf_all
  Use :: tcorr_var, Only: taurlim

  Implicit None

  Integer (i4b), Parameter :: nd = id_ndept
! .
! .. scalar arguments
  Logical :: expansion, outthb
  Integer (i4b) :: ncor
  Real (dp) :: emaxtc, teff, dtmean
! .
! .. array arguments
  Real (dp) :: t(nd), dtfcorr(nd), xne(nd), taur(nd)
! .
! .. local scalars

  Character (3) :: ncname
  Real (dp) :: rtau23, emaxne, tlim, taulast, w1, w2, wtot, dtout, dtin, ddt, &
    dtmax, dtmin
  Integer (i4b) :: i, irel, itaurlim, ilast, icount
! .
! .. local arrays
  Real (dp) :: deltat(nd), dummy(nd), r(nd), told(nd), dtold(nd), &
    xnesaved(nd), xnelte_taur(nd), taur_taur(nd), dt(nd), rel_err(nd)

  Logical :: flag

  Data dtold/nd*0.D0/

  deltat = 0.D0
  told = t
  Call thermal_balance(taur, xne, t, deltat, expansion, ncor+1, outthb, dtold, &
    rel_err, irel)

  Do i = nd, 1, -1
    If (taur(i)<taurlim) Exit
  End Do
  itaurlim = i + 1

! check for outermost point where dtfcorr is provided, including one safety
! point
  Do i = nd, 1, -1
    If (dtfcorr(i)==0.D0) Go To 100
  End Do
  Write (999, *) ' STOP: DTFCORR EQ 0 NOT FOUND'
  Stop ' DTFCORR EQ 0 NOT FOUND'

100 ilast = i + 1
  taulast = taur(ilast)
  If (optcmf_all) taulast = max(0.05D0, taulast) ! OTHERWISE, TOO LARGE
! CORRECTIONS

  Do i = nd, 1, -1
    If (taur(i)<3.*taulast) Exit !  ANOTHER SAFETY FACTOR, TO ALLOW
!   AT LEAST A FACTOR OF 3 FOR BORDER EFFECTS
  End Do
  ilast = i + 1
  If (taur(ilast)>taur(itaurlim)) ilast = itaurlim

! check for irel and modify it if necessary

  Print *, ' TAURLIM = ', taurlim
  Print *, ' TRANSITION REGION AT MAXIMUM UNTIL TAUR = ', taur(ilast)
  Print *, ' TAU(REL_ENER <  1.D-3) = ', taur(irel)

! itaurlim: minimum tau until we do flux-conservation in any case (inside)
! irel    : maximum tau from where on we do thermal balance in any case
! (outside)

  If (irel<ilast) irel = ilast !    THIS IS THE MAXIMUM WE ALLOW

  If (itaurlim-irel<2) Then
    irel = itaurlim - 2 !           at least two points for transition region
    tlim = taur(irel)
  Else
    tlim = taur(irel)
  End If

  If (optcmf_all) Then
!   JO March 2025
!   avoid to strong depression of T in case DTFCORR is negative and DELTAT
!   positive.
    If (dtfcorr(itaurlim)<0. .And. deltat(itaurlim)>0) Then
      Do i = itaurlim, nd
        If (dtfcorr(i)>=0.) Exit
        If (dtfcorr(i)<0. .And. deltat(i)<0.) Exit
      End Do
      itaurlim = i
      If (taur(i)>=2.) Then
        Print *, ' WARNING! WARNING! WARNING! TAURLIM MIGHT BE TOO LARGE!'
        Write (999, *) &
          ' WARNING! WARNING! WARNING! TAURLIM MIGHT BE TOO LARGE!'
      End If
      taurlim = taur(i)
      Print *
      Print *, ' TAURLIM CHANGED, SINCE OTHERWISE DT(FLUX) < 0 &
        &AND DT (THERMAL BALANCE) > 0'
      Print *
      Write (999, *)
      Write (999, *) ' TAURLIM CHANGED, SINCE OTHERWISE DT(FLUX) &
        &< 0 AND DT (THERMAL BALANCE) > 0'
      Write (999, *)
    End If
  End If

  Print *, ' FLUX-CORRECTION FROM FLUX-ERROR PERFORMED UNTIL TAUR > ', taurlim
  Print *, ' FLUX-CORRECTION FROM TH.BALANCE PERFORMED FROM  TAUR < ', tlim
  Print *, ' TRANSITION REGION FROM ', taur(itaurlim-1), ' UNTIL ', taur(irel)

  emaxtc = 0.D0
  Print *, 'TEMPERATURE CORRECTION!!!'
  Write (*, *) &
    ' ND   Tau Ross.      Temp.         dT (ThB)    dT (FlCorr.)    dT (used)'

  dtin = dtfcorr(itaurlim)
  dtout = deltat(irel-1)
  dtmin = min(dtin, dtout)
  dtmax = max(dtin, dtout)
  ddt = taur(itaurlim) - taur(irel-1)

  flag = .False.

  Do i = nd, 1, -1
    If (taur(i)<tlim) Then
!     outside
      dt(i) = deltat(i)
      t(i) = t(i) + dt(i)
      w1 = 1.
      w2 = 0.

    Else If (taur(i)>taurlim) Then
!     inside
      dt(i) = dtfcorr(i)
      t(i) = t(i) + dt(i)
      w1 = 0.
      w2 = 1.
    Else
!     intermediate range: interpolation
      wtot = taurlim - tlim
      If (wtot<=0.) Then
        Write (999, *) ' STOP: ERROR IN TAURLIM,TLIM'
        Stop ' ERROR IN TAURLIM,TLIM'
      End If
      w1 = 1. - exp(-(taurlim-taur(i))/wtot) ! this gives larger weights to
!     flux-corr
!     W1=1.-EXP(-(TAURLIM-TAUR(I))/WTOT/10.)
      w2 = 1. - w1
      dt(i) = w1*deltat(i) + w2*dtfcorr(i)
!     JO April 2017 (to avoid too large corrections, in particular if DTFCORR
!     is no
!     longer reliable): in this case, interpolate between inner (flux) and
!     outer (thermal balance) correction
      If (optcmf_all) Then
!       if alternative correction has been done once (FLAG=T), do following
!       interpolations
!       also with same method, to ensure a smooth transition
        If (dt(i)<dtmin .Or. dt(i)>dtmax .Or. flag) Then
          w1 = (taur(itaurlim)-taur(i))/ddt
          If (w1>1. .Or. w1<0.) Then
            Write (999, *) ' STOP: ERROR IN W1 (TEMPCORR_MAIN)'
            Stop ' ERROR IN W1 (TEMPCORR_MAIN)'
          End If
          w2 = 1. - w1
!         print*,'too large correction'
!         print*,i,dt(i),dtin,dtout,w1,w2
          dt(i) = w1*dtout + w2*dtin
          flag = .True.
        End If
      End If
!     for tests
!     IF(TAUR(I).LE.1. .AND. TAUR(I).GT. 0.01) DT(I)=DTFCORR(I)
      t(i) = t(i) + dt(i)
    End If

!   to avoid problems, we fix the minimum temperature that can be reached
!   at a given value
    If (t(i)<tmin1*teff) Then
      t(i) = tmin1*teff
      dt(i) = 0.
    End If
    If (t(i)<tminabs) Then
      t(i) = tminabs
      dt(i) = 0.
    End If
    emaxtc = max(emaxtc, abs(dt(i)/t(i)))
    Write (*, Fmt='(1X,I2,5(1X,G14.7),2(1X,F4.2),1X,L1)') i, taur(i), t(i), &
      deltat(i), dtfcorr(i), dt(i), w1, w2, flag
  End Do
  Print *, 'Maximum relative T correction ', emaxtc
  dtold = t - told

  If (outthb) Then
    Call ncor_name(ncor+1, ncname)
  Else
    ncname = '_LA'
  End If
  Open (7, File=trim(modnam)//'/DIFF_TC'//trim(ncname)//'.dat', &
    Status='UNKNOWN')
  Rewind 7
  Write (7, *) ' ND   Tau Ross.    Temp.       dT (ThB)    dT (FlCorr.) &
    &  dT (used)    XNE '
  Print *, ' '
  Do i = nd, 1, -1
    Write (7, Fmt='(1X,I3,1X,6(G12.5,1X))') i, taur(i), t(i), deltat(i), &
      dtfcorr(i), dt(i), xne(i)
  End Do
  Close (7)

! JO june 2018: calculate mean correction (used in main program)
  dtmean = 0.
  icount = 0
  Do i = 1, nd
    If (taur(i)<tlim .And. dt(i)/=0.) Then
      icount = icount + 1
      dtmean = dtmean + log10(abs(dt(i)/t(i)))
!     PRINT*,I,LOG10(ABS(DT(I)/T(I))),DTMEAN
    End If
  End Do
  dtmean = dtmean/float(icount)
! PRINT*,ICOUNT,DTMEAN


! updating the TEMP file (including XNELTE from last T-iteration)

  Rewind 21
  Read (21, *)(r(i), i=nd, 1, -1), (dummy(i), i=nd, 1, -1)
  Read (21, *)(xnesaved(i), i=nd, 1, -1)
  Read (21, *) rtau23

! just to be sure
  Open (1, File=trim(modnam)//'/TAU_ROS', Status='UNKNOWN', Form='FORMATTED')
  Rewind 1
  Do i = 1, nd
    Read (1, *) taur_taur(i), xnelte_taur(i)
    If (taur(i)/=0.D0 .And. abs(1.-taur(i)/taur_taur(i))>1.D-13) Then
      Write (999, *) ' STOP: TAUR .NE. TAUR IN TAU_ROS'
      Stop ' TAUR .NE. TAUR IN TAU_ROS'
    End If
    If (abs(1.-xnelte(i)/xnelte_taur(i))>1.D-13) Then
      Write (999, *) ' STOP: XNELTE .NE. XNELTE IN TAU_ROS' ! not exact, since
!     formattedly
!     written/read
      Stop ' XNELTE .NE. XNELTE IN TAU_ROS' ! not exact, since formattedly
!     written/read
    End If
  End Do
  Close (1)

  emaxne = 0.
  Do i = 1, nd
    emaxne = max(emaxne, abs(1.-xnelte(i)/xnesaved(i)))
  End Do

  Print *
  Print *, ' TEMPERATURE AND LTE ELECTRON DENSITY (LAST ITERATION) &
    &UPDATED IN FILE TEMP!'
  Print *, ' MAX. CHANGE IN XNELTE = ', emaxne
  Print *

  Rewind 21
  Write (21, *)(r(i), i=nd, 1, -1), (t(i), i=nd, 1, -1)
  Write (21, *)(xnelte(i), i=nd, 1, -1)
  Write (21, *) rtau23

  Return
End Subroutine

!________________________________________________________________________


Subroutine meanmolweight(meanmw, elecmw)

! check if EXPANSION is incorporated

  Use :: nlte_type
  Use :: nlte_dim
  Use :: princesa_var, Only: zeff, weight, abund
  Use :: nlte_var, Only: enionnd
  Use :: fund_const
  Implicit None

  Integer (i4b), Parameter :: nd = id_ndept, kel = id_atoms, kis = id_kisat

  Integer (i4b) :: i, i1, i2, aux, aux2
  Real (dp) :: ion_par(kel, kis+1, nd), meanmw(nd), elecmw(nd), sumaele, &
    sumapar

  ion_par(:, :, :) = 0.D0
  meanmw(:) = 0.D0
  elecmw(:) = 0.D0

  Do i = nd, 1, -1
    Do i1 = 1, kel
      sumapar = 0.D0
      aux = int(zeff(i1))
      Do i2 = 1, kis + 1
        sumapar = sumapar + enionnd(i1, i2, i)
      End Do
      Do i2 = aux + 1, kis + 1
        aux2 = i2 - aux
        ion_par(i1, i2, i) = enionnd(i1, aux2, i)/sumapar
      End Do
    End Do
    Do i1 = 1, kel
      sumapar = 0.D0
      sumaele = 0.D0
      aux = int(zeff(i1))
      Do i2 = 1, kis + 1
        sumapar = sumapar + ion_par(i1, i2, i)*(1.D0+aux)
        sumaele = sumaele + ion_par(i1, i2, i)*aux
        aux = aux + 1
      End Do
      sumapar = abund(i1)*sumapar/weight(i1)
      sumaele = abund(i1)*sumaele/weight(i1)
      meanmw(i) = meanmw(i) + sumapar
      elecmw(i) = elecmw(i) + sumaele
    End Do
    meanmw(i) = 1.D0/meanmw(i)
    elecmw(i) = 1.D0/elecmw(i)
  End Do

  Return

End Subroutine
!_________________________________________________________________


Subroutine thermal_balance(taur, xne, t, deltat, expansion, ncor, outthb, &
  dtold, rel_err, irel)

! It computes all the terms needed to check the thermal balance of the
! electrons, giving the T-Corr.
! If OUTTHB is true, then the numbers are save into a number of ascii
! files on every TCorr cycle. WARNING: if a model is restarted to perform
! some new T corrections, then the output files will be overwritten


! rates always formulated via heating-cooling:
! a final positive term means that the net process heats the medium

  Use :: nlte_type
  Use :: nlte_dim
  Use :: fund_const

  Use :: nlte_var, Only: ifre, fre, wfre, xj, modnam, optcmf_all

  Use :: tcorr_var, Only: ffopa, dtffopa, qffh, qffc, qcbfr, qcbfi, qrbfr, &
    qrbfi, dtqrbfr, dtqrbfi, dtqffc, dtqffh, dtqcbfh, dtqcbfc, qcbbu, qcbbd, &
    dtqcbbu, dtqcbbd, ffopa_m, dtffopa_m, qcbfr_m, qcbfi_m, dtqcbfr_m, &
    dtqcbfi_m, qrbfr_m, qrbfi_m, dtqrbfr_m, dtqrbfi_m, qcbbu_m, qcbbd_m, &
    dtqcbbu_m, dtqcbbd_m


  Implicit None

  Integer (i4b), Parameter :: nd = id_ndept, nf = id_frec1
! .
! .. scalar arguments
  Integer (i4b) :: ncor, irel
  Logical :: expansion, outthb
! .
! .. array arguments
  Real (dp) :: xne(nd), deltat(nd), t(nd), taur(nd), dtold(nd), rel_err(nd)
! .
! .. local scalars
  Character (3) :: ncname
  Integer (i4b) :: i, j, iifreq
  Real (dp) :: aux, x, expo
! .
! .. local arrays
  Real (dp), Dimension (nd, nf) :: x1, x2, x3, x4, x5
  Real (dp), Dimension (nd) :: jint, bint
  Real (dp) :: energies(nd, 9), dtq(nd, 9), taux1(nd, 6), taux2(nd, 6), &
    taux3(nd, 6), meanmw(nd), dummy(nd), rel_ene(nd, 6), enerdiffold(nd), &
    dtqold(nd), auxold(nd)
! ,REL_ENE2(ND,5),RADIATIVEEQ(ND,9),AUX1(ND,NF)

  Data enerdiffold/nd*0.D0/
! ..
! .. external subroutines ..
  External :: meanmolweight, ncor_name
! ..
! .. external functions
  Real (dp) :: xintfre, fftot, dtfftot
  External :: xintfre
! ..
! .. IntrInsIc functIons ..
! INTRINSIC EXP, MAX, MIN, ABS


  If (optcmf_all) Then
    If (xj(1,ifre+1)/=2) Then
      Write (999, *) ' STOP: OPTCMF_ALL = TRUE, BUT APPROX. XJ &
        &(OBSFRAM) USED IN THERMAL_BALANCE'
      Stop &
        ' OPTCMF_ALL = TRUE, BUT APPROX. XJ (OBSFRAM) USED IN THERMAL_BALANCE'
    End If
  End If

  aux = 0.D0
  energies = 0.D0
  dtq = 0.D0
  x1 = 0.D0
  x2 = 0.D0
  x3 = 0.D0
  x4 = 0.D0
  x5 = 0.D0
! RADIATIVEEQ=0.D0

  Do i = 1, nd !                    depth loop
!   FREE-FREE, other terms are calculated In NETMAT, RATEEQ and LINESOB
    Do iifreq = 1, ifre !           frequency loop
      aux = (xj(i,iifreq)+hc2*fre(iifreq)**3)*exp(-hkl*fre(iifreq)/t(i))
      fftot = ffopa(i, iifreq) + ffopa_m(i, iifreq)
      dtfftot = dtffopa(i, iifreq) + dtffopa_m(i, iifreq)

      x1(i, iifreq) = fftot*xj(i, iifreq)
      x3(i, iifreq) = fftot*aux
!     derivatives vs. temperature
      x2(i, iifreq) = dtfftot*xj(i, iifreq)
      x4(i, iifreq) = aux*(dtfftot+fftot*hkl*fre(iifreq)/(t(i)**2))
!     Planck function, to check consistency
      x = hkl*fre(iifreq)/t(i)
      If (x<200.D0) Then
        x5(i, iifreq) = hc2*fre(iifreq)**3/(exp(x)-1.D0)
      Else
        expo = log(hc2*fre(iifreq)**3) - x
        x5(i, iifreq) = exp(expo)
      End If
    End Do

    qffh(i) = 4.D0*pi*xintfre(x1(i,:), wfre, 1, ifre, ifre, fre, 4, 'OLD')
    qffc(i) = 4.D0*pi*xintfre(x3(i,:), wfre, 1, ifre, ifre, fre, 4, 'OLD')
    dtqffh(i) = 4.D0*pi*xintfre(x2(i,:), wfre, 1, ifre, ifre, fre, 4, 'OLD')
    dtqffc(i) = 4.D0*pi*xintfre(x4(i,:), wfre, 1, ifre, ifre, fre, 4, 'OLD')
    jint(i) = xintfre(xj(i,:), wfre, 1, ifre, ifre, fre, 4, 'OLD')
    bint(i) = xintfre(x5(i,:), wfre, 1, ifre, ifre, fre, 4, 'OLD')
  End Do

! TEST OF FREQ. GRID AND INTEGRATION, JUST IN CASE
  aux = sigsb/pi*t(1)**4
  If (abs(1.-bint(1)/aux)>1.D-2) Then
    Print *, aux, bint(1)
    Write (999, *) &
      ' STOP: THERMAL BALANCE: SOMETHING WRONG IN FREQ. INTEGRATION (OUTSIDE)'
    Stop ' THERMAL BALANCE: SOMETHING WRONG IN FREQ. INTEGRATION (OUTSIDE)'
  End If
  aux = sigsb/pi*t(nd)**4
  If (abs(1.-bint(nd)/aux)>1.D-2) Then
    Print *, aux, bint(nd)
!   STOP ' THERMAL BALANCE: SOMETHING WRONG IN FREQ. INTEGRATION (INSIDE)'
    Print *, &
      ' WARNING!! THERMAL BALANCE: FREQ. INTEGRATION (INSIDE) INACCURATE'
    Print *, &
      ' WARNING!! THERMAL BALANCE: FREQ. INTEGRATION (INSIDE) INACCURATE'
    Print *, &
      ' WARNING!! THERMAL BALANCE: FREQ. INTEGRATION (INSIDE) INACCURATE'
    Write (999, *) &
      ' WARNING!! THERMAL BALANCE: FREQ. INTEGRATION (INSIDE) INACCURATE'
    Write (999, *) &
      ' WARNING!! THERMAL BALANCE: FREQ. INTEGRATION (INSIDE) INACCURATE'
    Write (999, *) &
      ' WARNING!! THERMAL BALANCE: FREQ. INTEGRATION (INSIDE) INACCURATE'
  End If

  Do i = 1, nd !                    depth loop
!   G heating
!   L cooling
!   negative rates=cooling
    energies(i, 1) = qffh(i) !      Gamma_FF
    energies(i, 2) = qffc(i) !      Lambda_FF
    dtq(i, 1) = dtqffh(i) !         dGamma_FF/dT
    dtq(i, 2) = dtqffc(i) !         dLambda_FF/dT

!   radiative equilibrium terms
!   RADIATIVEEQ(I,1) = RADIATIVEEQ(I,1) + QFFC(I)	!REQ heatIng
!   RADIATIVEEQ(I,2) = RADIATIVEEQ(I,2) + QFFH(I)	!REQ coolIng


    energies(i, 3) = qcbfr(i) !     G CBF
    energies(i, 4) = qcbfi(i) !     L CBF
    energies(i, 5) = qrbfi(i) !     G RBF
    energies(i, 6) = qrbfr(i) !     L RBF
    energies(i, 7) = qcbbd(i) !     G CBB
    energies(i, 8) = qcbbu(i) !     L CBB

    dtq(i, 3) = dtqcbfh(i) !        dG/dT CBF
    dtq(i, 4) = dtqcbfc(i) !        dL/dT CBF
    dtq(i, 5) = dtqrbfi(i) !        dG/dT RBF
    dtq(i, 6) = dtqrbfr(i) !        dL/dT RBF
    dtq(i, 7) = dtqcbbd(i) !        dG/dT CBB
    dtq(i, 8) = dtqcbbu(i) !        dL/dT CBB
!   radiative equilibrium terms
!   do IINAT=1,NAT ! atom loop
!   RADIATIVEEQ(I,3) = RADIATIVEEQ(I,3) + RERBFH(I,IINAT)
!   RADIATIVEEQ(I,4) = RADIATIVEEQ(I,4) + RERBFC(I,IINAT)
!   enddo
    dtq(i, 3) = dtq(i, 3) + dtqcbfr_m(i) ! dG/dT CBF_METALS
    dtq(i, 4) = dtq(i, 4) + dtqcbfi_m(i) ! dL/dT CBF_METALS
    dtq(i, 5) = dtq(i, 5) + dtqrbfi_m(i) ! dG/dT RBF_METALS
    dtq(i, 6) = dtq(i, 6) + dtqrbfr_m(i) ! dL/dT RBF_METALS
    dtq(i, 7) = dtq(i, 7) + dtqcbbd_m(i) ! dG/dT CBB_METALS
    dtq(i, 8) = dtq(i, 8) + dtqcbbu_m(i) ! dL/dT CBB_METALS
!   in case, print out the different contributions
!   print*,i,energies(i,3)-energies(i,4),qcbfr_m(i)-qcbfi_m(i)
!   print*,i,energies(i,5)-energies(i,6),qrbfi_m(i)-qrbfr_m(i)
!   print*,i,energies(i,7)-energies(i,8),qcbbd_m(i)-qcbbu_m(i)
!   print*
    energies(i, 3) = energies(i, 3) + qcbfr_m(i) ! G CBF_METALS
    energies(i, 4) = energies(i, 4) + qcbfi_m(i) ! L CBF_METALS
    energies(i, 5) = energies(i, 5) + qrbfi_m(i) ! G RBF_METALS
    energies(i, 6) = energies(i, 6) + qrbfr_m(i) ! L RBF_METALS
    energies(i, 7) = energies(i, 7) + qcbbd_m(i) ! G CBB_METALS
    energies(i, 8) = energies(i, 8) + qcbbu_m(i) ! L CBB_METALS
!   for tests
!   PRINT*,'THERMAL_BAL',I,QCBBD(I),QCBBU(I),QCBBD_M(I),QCBBU_M(I)

    energies(i, 9) = (energies(i,1)-energies(i,2)) + &
      (energies(i,3)-energies(i,4)) + (energies(i,5)-energies(i,6)) + &
      (energies(i,7)-energies(i,8)) ! Sum (G - L)
    dtq(i, 9) = (dtq(i,1)-dtq(i,2)) + (dtq(i,3)-dtq(i,4)) + &
      (dtq(i,5)-dtq(i,6)) + (dtq(i,7)-dtq(i,8)) ! Sum (dG/dT - dL/dT)

!   RADIATIVEEQ(I,7) = RADIATIVEEQ(I,1) + RADIATIVEEQ(I,3) + &
!   RADIATIVEEQ(I,5)
!   RADIATIVEEQ(I,8) = RADIATIVEEQ(I,2) + RADIATIVEEQ(I,4) + &
!   RADIATIVEEQ(I,6)
!   RADIATIVEEQ(I,9) = RADIATIVEEQ(I,7) - RADIATIVEEQ(I,8)


!   relative differences
    rel_ene(i, 1) = 1.D0 - energies(i, 1)/energies(i, 2) ! 1-g/l FF
    rel_ene(i, 2) = 1.D0 - energies(i, 3)/energies(i, 4) ! 1-g/l CBF

    If (energies(i,6)/=0.D0) Then
      rel_ene(i, 3) = 1.D0 - energies(i, 5)/energies(i, 6) ! 1-g/l rbf
    Else
      If (energies(i,5)/=0.D0) Then
        Write (999, *) ' STOP: SOMETHING WRONG WITH RBF ENERGIES'
        Stop ' SOMETHING WRONG WITH RBF ENERGIES'
      End If
      rel_ene(i, 3) = 0.D0
    End If

    If (energies(i,8)/=0.D0) Then
      rel_ene(i, 4) = 1.D0 - energies(i, 7)/energies(i, 8) ! 1-G/L CBB
    Else
      If (energies(i,7)/=0.D0) Then
        Write (999, *) ' STOP: SOMETHING WRONG WITH CBB ENERGIES'
        Stop ' SOMETHING WRONG WITH CBB ENERGIES'
      End If
      rel_ene(i, 4) = 0.D0
    End If

    rel_ene(i, 5) = 1.D0 - (energies(i,1)+energies(i,3)+energies(i,5)+energies &
      (i,7))/(energies(i,2)+energies(i,4)+energies(i,6)+energies(i,8)) ! 1-G/L
!   ALL

    rel_ene(i, 6) = 1.D0 - jint(i)/bint(i) ! 1-J/B

!   relative differences for partials (not required presently)
!   rel_ene2(I,1)=1.D0-ABS(DTQ(I,1))/ABS(DTQ(I,2))	! 1-abs(Dg)/abs(Dl) FF
!   rel_ene2(I,2)=1.D0-ABS(DTQ(I,3))/ABS(DTQ(I,4))	! 1-abs(Dg)/abs(Dl) CBF
!   rel_ene2(I,3)=1.D0-ABS(DTQ(I,5))/ABS(DTQ(I,6))	! 1-abs(Dg)/abs(Dl) RBF

!   IF (DTQ(I,8) .NE. 0.D0) THEN
!   rel_ene2(I,4)=1.D0-ABS(DTQ(I,7))/ABS(DTQ(I,8)) ! 1-g/l CBB
!   ELSE
!   IF (DTQ(I,7) .ne. 0.D0) STOP ' something wrong with cbb partials'
!   rel_ene2(I,4)=0.D0
!   ENDIF
!   REL_ENE2(I,5) = 0.D0

!   absolute differences
    taux1(i, 1) = energies(i, 1) - energies(i, 2) ! G - L FF
    taux1(i, 2) = energies(i, 3) - energies(i, 4) ! G - L CBF
    taux1(i, 3) = energies(i, 5) - energies(i, 6) ! G - L RBF
    taux1(i, 4) = energies(i, 7) - energies(i, 8) ! G - L CBB
    taux1(i, 5) = taux1(i, 1) + taux1(i, 2) + taux1(i, 3) + taux1(i, 4)
    taux1(i, 6) = rel_ene(i, 5) !   1-G/L ALL

    taux2(i, 1) = dtq(i, 1) - dtq(i, 2) ! dG/dT - dL/dT FF
    taux2(i, 2) = dtq(i, 3) - dtq(i, 4) ! dG/dT - dL/dT CBF
    taux2(i, 3) = dtq(i, 5) - dtq(i, 6) ! dG/dT - dL/dT RBF
    taux2(i, 4) = dtq(i, 7) - dtq(i, 8) ! dG/dT - dL/dT CBB
    taux2(i, 5) = taux2(i, 1) + taux2(i, 2) + taux2(i, 3) + taux2(i, 4)

  End Do !                          depth loop


  deltat = 0.D0 !                   T-Corr
! call EXPAN_COOLING(R,XNE,V,T,DELTATE,ENERGIES(:,9),DTQ(:,9),VDIV,OPT, &
! &                  IND2,RHO,XMLOSS,ESCALARADIO,ESCALAVELO)
  Call meanmolweight(meanmw, dummy)
  Do i = 1, nd
    taux3(i, 1) = meanmw(i)
    taux3(i, 2) = t(i)
    taux3(i, 3) = taur(i)
    aux = -energies(i, 9)/dtq(i, 9) ! delta T = (G-L)/(dL/dT-dG/dT)
!   in case the temperature iteration (controlled by the partials) does
!   not work, consider the following possibility (set taur(i) to 1.d-2).
    If (enerdiffold(i)/=0. .And. dtold(i)/=0. .And. taur(i)<-100.) Then
      dtqold(i) = (energies(i,9)-enerdiffold(i))/dtold(i)
      If (enerdiffold(i)*energies(i,9)<0.D0) Then
        auxold(i) = -energies(i, 9)/dtqold(i) ! BISECT
      Else If (abs(energies(i,9))<abs(enerdiffold(i))) Then
        auxold(i) = (.5D0*energies(i,9)-energies(i,9))/dtqold(i)
      Else
        auxold(i) = (.5D0*enerdiffold(i)-energies(i,9))/dtqold(i)
      End If
    Else
      auxold(i) = aux
    End If
    aux = auxold(i)
!   restrict correction to 7.5%
    If (abs(aux)>0.075*t(i)) Then
      taux3(i, 4) = aux/abs(aux)*0.075*t(i) ! correct sign
    Else
      taux3(i, 4) = aux
    End If
    taux3(i, 5) = taux3(i, 4)/t(i) ! delta T / T
    taux3(i, 6) = xne(i) !          electron density
!   TAUX3(I,5)=DELTAT(I)		     ! delta T from EXPAN_COOLING
    deltat(i) = taux3(i, 4) !       ThB  deltaT-Corr
  End Do
  enerdiffold = energies(:, 9)

  If (outthb) Then
    Call ncor_name(ncor, ncname)
  Else
    ncname = '_LA'
  End If
  Open (7, File=trim(modnam)//'/ENERGY_TC'//trim(ncname)//'.dat', &
    Status='UNKNOWN')
  Rewind 7
  Do i = 1, nd
    Write (7, Fmt='(1X,I3,1X,9(G12.5,1X))') i, (energies(i,j), j=1, 9)
  End Do
  Close (7)
  Open (7, File=trim(modnam)//'/DERIVATIVES_TC'//trim(ncname)//'.dat', &
    Status='UNKNOWN')
  Rewind 7
  Do i = 1, nd
    Write (7, Fmt='(1X,I3,1X,9(G12.5,1X))') i, (dtq(i,j), j=1, 9)
  End Do
  Close (7)
  Open (7, File=trim(modnam)//'/ENERDIFF_TC'//trim(ncname)//'.dat', &
    Status='UNKNOWN')
  Rewind 7
  Write (7, *) ' ND   FF (G-L)     CBF (G-L)     RBF (G-L)     &
    &CBB (G-L)    TOTAL     1-G/L ALL'
  Do i = 1, nd
    Write (7, Fmt='(1X,I3,1X,6(G12.5,1X))') i, (taux1(i,j), j=1, 6)
  End Do
  Close (7)
  Open (7, File=trim(modnam)//'/DERIVADIFF_TC'//trim(ncname)//'.dat', &
    Status='UNKNOWN')
  Rewind 7
  Write (7, *) 'dG/dT-dL/dT   FF      CBF           RBF        &
    &CBB          TOTAL       DUMMY  '
  Do i = 1, nd
    Write (7, Fmt='(1X,I3,1X,5(G12.5,1X))') i, (taux2(i,j), j=1, 5)
  End Do
  Close (7)
  Open (7, File=trim(modnam)//'/DELTAT_TC'//trim(ncname)//'.dat', &
    Status='UNKNOWN')
  Rewind 7
! write(7,*) ' ND  MeanMolW      T (K)        Tau Ross      dT (K)
! dT/T         XNE '
  Do i = 1, nd
    Write (7, Fmt='(1X,I3,1X,6(G12.5,1X))') i, (taux3(i,j), j=1, 6)
  End Do
  Close (7)
  Open (7, File=trim(modnam)//'/REL_ENERGY_TC'//trim(ncname)//'.dat', &
    Status='UNKNOWN')
  Rewind 7
  Write (7, *) ' ND   1-g/l FF    1-g/l CBF      1-g/l RBF    1-g/l &
    &CBB    1-G/L ALL   1-J/B '
  Do i = 1, nd
    Write (7, Fmt='(1X,I3,1X,6(G12.5,1X))') i, (rel_ene(i,j), j=1, 6)
  End Do
  Close (7)

! open (7,FILE=TRIM(MODNAM)//'/RADEQ_ENERGIES.dat',STATUS='UNKNOWN')
! rewind 7
! do I=1,ND
! write(7,FMT='(1X,I3,1X,9(G12.5,1X))') I,(radIatIveeq(I,tecla),&
! &		tecla=1,9)
! enddo
! close (7)

! find point where relative differences are too low to allow for successful
! correction

  rel_err = taux1(:, 6)

  Do i = nd - 1, 1, -1
    If (abs(taux1(i,6))>1.D-3) Go To 100
  End Do

  Do i = nd - 1, 1, -1
    If (abs(taux1(i,6))>1.D-4) Go To 100
  End Do

  Write (999, *) ' STOP: REL DIFF IN HEATING/COOLING NEVER LARGER THAN 1.D-4!'
  Stop ' REL DIFF IN HEATING/COOLING NEVER LARGER THAN 1.D-4!'

100 irel = i + 1


  Return
End Subroutine

!-----------------------------------------------------------------

Subroutine ncor_name(n, nname)

  Use :: nlte_type
  Use :: nlte_dim
  Use :: fund_const

  Implicit None
  Integer (i4b) :: n
  Character (3) :: nname
  Integer (i4b) :: i, j


  i = n/10
  If (i<1) Then
    nname = trim(char(n+48))
  Else
    j = n - i*10
    nname = trim(char(i+48)) // trim(char(j+48))
  End If
  nname = trim(nname)

  Return
End Subroutine

!***********************************************************************

!subroutines: recent ones (from vers. 9.0 on)

!***********************************************************************

Subroutine op_rbfsetup(ii, nunname, ndatosmax, nrecord, ndatos, fileop)

! for formula 20
! changed Sept. 2016. most important, FRECFN1 is calculated and overwrites
! the start value (=FRECFN)

  Use :: nlte_type
  Use :: nlte_dim
  Use :: fund_const
  Use :: fastwind_params, Only: oppath, resol

  Use :: princesa_var, Only: labl, labl4, frecin, frecfn, fl
  Use :: nlte_var, Only: ifre, fre, eopdiff, op_flag2, op_data, modnam, &
    frecfn1
  Use :: nlte_opt, Only: opdebug


  Implicit None
! ..
! .. parameters ..
! INTEGER(I4B), PARAMETER :: IMAX=4000
  Integer (i4b), Parameter :: imax = 6000 ! MAU (Apr 2023)

! JP (July 2013)
! INTEGER(I4B), PARAMETER :: NMAX_LOW = 2000 ! for structured part
! INTEGER(I4B), PARAMETER :: NMAX_HIGH = 200 ! for smooth high freq. tail
! INTEGER(I4B), PARAMETER :: NMAX_HIGH = 2000 ! for smooth high freq. tail
! MAU (Apr 2023)
  Integer (i4b), Parameter :: nmax_low = 3000
  Integer (i4b), Parameter :: nmax_high = 3000


! ..
! .. scalar arguments ..

  Integer (i4b) :: ii, ndatosmax, nrecord, ndatos
  Real (dp) :: nunname
  Character (30) :: fileop

! ..
! .. local scalars

  Integer (i4b) :: j, estado, k, l1, l2, ivoid, nout, l2start
  Real (dp) :: xnaux, atmryd, dvoid, eryd, reldiff, q, q1, alpha, xnauxmin

  Character (30) :: checkfile
  Character (80) :: dummychar

! ..
! .. local arrays ..

  Real (dp) :: fgs(imax), pics(imax), aux2(imax), fgsold(2)
  Real (dp), Dimension (nmax_low+nmax_high) :: eout, picsout

  If (ndatosmax>imax) Then
    Write (999, *) ' STOP: NDATOSMAX > IMAX'
    Stop ' NDATOSMAX > IMAX'
  End If

  fgs(:) = 0.D0
  pics(:) = 0.D0
  atmryd = 0.D0
  estado = 1
  eryd = 0.D0
  alpha = 0.D0
  reldiff = 0.D0

  op_flag2(ii) = .False.

  Call op_codex(int(nunname), fileop, checkfile, atmryd) ! filename decoder
! (CODEX)

  Open (771, File=oppath//trim(fileop), Status='OLD', Access='DIRECT', &
    Recl=2*ndatosmax*8, Iostat=estado)

  If (estado/=0) Then
    Print *, 'FILE OP_RBF ', oppath // trim(fileop), ' NOT FOUND'
    Write (999, *) ' STOP!'
    Stop
  End If

  Read (771, Rec=nrecord)(fgs(k), k=1, ndatos), (pics(j), j=1, ndatos)
  Close (771)

  If (opdebug) Then
    dummychar = trim(modnam) // '/opori-crossbf-' // trim(labl(labl4(ii))) // &
      '.ascii'
    Open (771, File=trim(dummychar), Status='UNKNOWN')
    Write (771, Fmt='(1X,A10,3(1X,G18.8))') trim(labl(labl4(ii))), eryd, &
      fl(labl4(ii))/clight, frecin(ii)
    Do k = 1, ndatos
      xnaux = fgs(k)/clight
      Write (771, Fmt='(1X,2(G18.8,1X))') xnaux, pics(k)
    End Do
    Close (771)
  End If

! interpolate OPACITY data when zero values are present

  Do k = 1, ndatos
    If (pics(k)==0.) Then
      Go To 100
    End If
  End Do

  Go To 110

100 Continue
  If (opdebug) Then
    Print *
    Print *, ' WARNING!!! ZEROES FOUND IN OPACITY DATA FOR LEVEL ', &
      trim(labl(labl4(ii)))
    Write (999, *)
    Write (999, *) ' WARNING!!! ZEROES FOUND IN OPACITY DATA FOR LEVEL ', &
      trim(labl(labl4(ii)))
  End If

  Call op_fit(fgs, pics, imax, ndatos)

  If (opdebug) Then
    dummychar = trim(modnam) // '/opori-mod-crossbf-' // &
      trim(labl(labl4(ii))) // '.ascii'
    Open (771, File=trim(dummychar), Status='UNKNOWN')
    Write (771, Fmt='(1X,A10,3(1X,G18.8))') trim(labl(labl4(ii))), eryd, &
      fl(labl4(ii))/clight, frecin(ii)
    Do k = 1, ndatos
      xnaux = fgs(k)/clight
      Write (771, Fmt='(1X,2(G18.8,1X))') xnaux, pics(k)
    End Do
    Close (771)
  End If
! *******************


! Is there any offset in the (energy) level definition between the OP data
! and DETAIL?
! This is only checked once.

! We need to compare the FL variable (ionization energy w.r.t. ground-state)
! with the OP energy (w.r.t. ground-state), not FRECIN (actual ionization
! energy),
! because an ionization transition to an excited level is possible.
! Remember that [FL] = s^{-1};
! NOTE: variable FRECFN consistent with FL (though in Kayser)

110 If (.Not. op_flag2(ii)) Then
    estado = 1
    Open (771, File=oppath//trim(checkfile), Status='OLD', Iostat=estado)
    If (estado/=0) Then
      Print *, ' CHECK FILE OP_RBF ', trim(checkfile), ' NOT FOUND'
      Write (999, *) ' STOP!'
      Stop
    End If
    Read (771, *) dummychar
    Read (771, *) dummychar

    Do l1 = 1, nrecord
      Read (771, *) ivoid, ivoid, dvoid, dvoid, dvoid, eryd, dvoid
    End Do
    Close (771)
    eryd = eryd*(-1.D0)*atmryd
!   eopdiff(ii) = eryd - frecin(ii)
!   reldiff=eopdiff(ii)/frecin(ii)

    xnaux = fl(labl4(ii))/clight
    If (abs(1.-xnaux/frecfn(ii))>1.D-15) Then
      Write (999, *) ' STOP: INCONSISTENCY BETWEEN FL AND FRECFN!'
      Stop ' INCONSISTENCY BETWEEN FL AND FRECFN!'
    End If
    eopdiff(ii) = eryd - xnaux
    reldiff = eopdiff(ii)/xnaux

!   WRITE(*,*) ' Relative energy difference at the edge ', &
!   TRIM(LABL(LABL4(II))),' ', RELDIFF

!   Just a warning, proceed at user's risk
    If (abs(reldiff)>0.2) Then
      Print *, ' RBF OP-DATA WARNING -- level ', trim(labl(labl4(ii))), &
        ' shows a large (rel) energy difference ', reldiff
      Print *, eryd, xnaux, eopdiff(ii)
      Write (999, *) &
        ' STOP: TOO LARGE DIFFERENCE IN IONIZATION EDGES OP-DATA VS. DETAIL'
      Stop ' TOO LARGE DIFFERENCE IN IONIZATION EDGES OP-DATA VS. DETAIL'
    End If

    If (abs(reldiff)>1.D-3) Then
      Print *
      Print *, 'RBF OP-DATA WARNING -- shifting cross section for ', &
        trim(labl(labl4(ii))), ' from ', eryd, ' cm^(-1) to ', &
        eryd - eopdiff(ii), ' cm^(-1)'
      fgsold(1) = fgs(ndatos)/clight
      fgsold(2) = fgs(1)/clight
      fgs(:) = fgs(:) - (eopdiff(ii)*clight)
    End If

    op_flag2(ii) = .True.
    If (opdebug) Then
      dummychar = trim(modnam) // '/op-crossbf-' // trim(labl(labl4(ii))) // &
        '.ascii'
      Open (771, File=trim(dummychar), Status='UNKNOWN')
      Write (771, Fmt='(1X,A10,3(1X,G18.8))') trim(labl(labl4(ii))), eryd, &
        fl(labl4(ii))/clight, frecin(ii)
      Do k = 1, ndatos
        xnaux = fgs(k)/clight
        Write (771, Fmt='(1X,2(G18.8,1X))') xnaux, pics(k)
      End Do
      Close (771)

      dummychar = trim(modnam) // '/detsrf-crossbf-' // &
        trim(labl(labl4(ii))) // '.ascii'
      Open (771, File=trim(dummychar), Status='UNKNOWN')
      Write (771, Fmt='(1X,A10,3(1X,G18.8))') trim(labl(labl4(ii))), eryd, &
        fl(labl4(ii)), frecin(ii)*clight
      Write (771, Fmt='(1X,A2,1X,G18.8,1X,A6)') 'F1', fgs(1), &
        trim(labl(labl4(ii)))
      Write (771, Fmt='(1X,A2,1X,G18.8,1X,A6)') 'F2', fgs(ndatos), &
        trim(labl(labl4(ii)))

      Close (771)
    End If

  End If

  aux2(:) = fgs(:)/clight

  If (aux2(ndatos)>frecin(ii)) Then
    Print *, &
      ' OP STOP  --- rbf cross-section does not include the ionization edge'
    Print *, '   Level ', trim(labl(labl4(ii))), '  -----   record ', nrecord
    Print *, '   Level defined at ', fl(labl4(ii))/clight, ' cm^{-1}'
    Print *, '   Ionization edge located at ', frecin(ii), ' cm^{-1}'
    Print *, '   while the available range is ', aux2(ndatos), ' -- ', &
      aux2(1), ' cm^{-1} '
    Print *, '   <--> before correction was: ', fgsold(1), ' -- ', fgsold(2)
    Print *, ' CHECK THE UPPER LEVEL for this transition !!!! '
    Write (999, *) ' STOP!'
    Stop
  End If

! smoothing with Gaussian: energy units used!
! the logical variable controls output (no output if set to .false.)

! CALL
! SMOOTH(AUX2,PICS,NDATOS,RESOL,NMAX_LOW,NMAX_HIGH,FRECIN(II),EOUT,PICSOUT,NOUT,.false.)
! JO: changed Sept. 2016, to allow for frequencies lower than the
! standard ionization energy. Note that for ionizations to excited levels,
! resonances are possible in the range down or even below the ground-state
! edge


  Call smooth(aux2, pics, ndatos, resol, nmax_low, nmax_high, aux2(ndatos), &
    eout, picsout, nout, .False.)


  If (eout(1)/=fgs(1)/clight) Then
    Write (999, *) ' STOP: ERROR IN FGS(1)' ! highest energy, MUST be
!   identical
    Stop ' ERROR IN FGS(1)' !       highest energy, MUST be identical
  End If
  If (abs(1.D0-eout(nout)/(fgs(ndatos)/clight))>1.D-14) Then
    Print *, eout(nout), fgs(ndatos)/clight
    Write (999, *) ' STOP: ERROR IN FGS(NDATOS)' ! lowest energy, identical
!   within precision
    Stop ' ERROR IN FGS(NDATOS)' !  lowest energy, identical within precision
  End If
  eout(nout) = fgs(ndatos)/clight ! to remain consistent

! Once the cross-section has been read and convolved, &
! we map the cross-section on FASTWIND's frequency grid

! units: fre, frecin, eout in cm^-1; xnaux, fgs Hz
! energy grid fre  ordered from red to blue
! energy grid eout ordered from blue to red

  l2start = nout - 1

! JO Sept2016; for some transitions, data not down to groundstate edge
  If (eout(nout)>frecfn(ii)) Then
    If (eout(nout)>frecin(ii)) Then
      Write (999, *) ' STOP: SOMETHING WRONG WITH OP-DATA (SUBR. OP_RBFSETUP)'
      Stop ' SOMETHING WRONG WITH OP-DATA (SUBR. OP_RBFSETUP)'
    End If
  End If

  xnauxmin = 0.
  Do l1 = 1, ifre
    xnaux = fre(l1)
!   IF (FRE(L1) .LE. EOUT(1) .AND. FRE(L1) .GE. FRECIN(II) ) THEN
!   JO: changed Sept. 2016, to allow for frequencies lower than the
!   standard ionization energy. Note that for ionizations to excited levels,
!   resonances are possible in the range down or even below the ground-state
!   edge

    If (xnaux<=eout(1) .And. xnaux>=eout(nout)) Then
      If (xnauxmin==0.) xnauxmin = xnaux
!     emax > fre(l1) > edge, starting close to edge
!     we do not allow lower frequencies than the effective edge

      Do l2 = l2start, 1, -1
        If (xnaux<=eout(l2) .And. xnaux>eout(l2+1)) Then
          j = l2
          k = l2 + 1
          If (picsout(j)<=0.D0) Then
            Write (999, *) ' STOP: SOUT(J) LE 0.'
            Stop ' SOUT(J) LE 0.'
          End If
          If (picsout(k)<=0.D0) Then
            Write (999, *) ' STOP: SOUT(K) LE 0.'
            Stop ' SOUT(K) LE 0.'
          End If
          Go To 120
        End If
      End Do
      Write (999, *) ' STOP: FRE NOT FOUND IN EOUT'
      Stop ' FRE NOT FOUND IN EOUT'
!     logarithmic interpolation
120   l2start = l2
      q = log10(xnaux/eout(k))/log10(eout(j)/eout(k))
      q1 = 1.D0 - q
      alpha = q*log10(picsout(j)) + q1*log10(picsout(k))
      alpha = 10.D0**alpha

      If (alpha<=0.D0) Then
        Print *, 'ERROR IN ', trim(oppath//trim(fileop)), &
          ' alpha zero or negative'
        Print *, '   record ', nrecord, '  level ', trim(labl(labl4(ii)))

        dummychar = trim(modnam) // '/crossbf-error-' // &
          trim(labl(labl4(ii))) // '.dat'
        Open (77, File=trim(dummychar), Status='UNKNOWN')
        Do k = 1, ifre
          Write (77, *) fre(k), op_data(ii, k)
        End Do
        Close (77)
        Write (999, *) ' STOP!'
        Stop
      End If

      op_data(ii, l1) = alpha

    Else If (xnaux>eout(1)) Then !  nu**3 extrapolation for high energies,
!     just in case
      op_data(ii, l1) = picsout(1)*(eout(1)/xnaux)**3

    Else !                          just in case, we will never reach this
!     frequency
      op_data(ii, l1) = 0.D0
    End If

  End Do

! redefine FRECFN1
  frecfn1(ii) = max(frecfn(ii), xnauxmin)

! PRINT*,LABL(LABL4(II)),FRECIN(II),FRECFN(II),FRECFN1(II)

  If (opdebug) Then
    dummychar = trim(modnam) // '/fstwnd-crossbf-' // trim(labl(labl4(ii))) // &
      '.ascii'
    Open (771, File=trim(dummychar), Status='UNKNOWN')
    Write (771, Fmt='(1X,A10,2(1X,G18.8))') trim(labl(labl4(ii))), &
      fl(labl4(ii))/clight, frecin(ii)
    Do k = 1, ifre
      xnaux = fre(k)*clight
      Write (771, Fmt='(1X,2(G18.8,1X))') fre(k), op_data(ii, k)
    End Do
    Close (771)
  End If

  Print *, 'TRANSITION FROM ', labl(labl4(ii)), ' DONE'

  Return
End Subroutine

!-----------------------------------------------------------------

Subroutine op_rbfsetup1(ii, key, iz, nrecord, fileop)

! for formula 101
! similar to OP_RBFSETUP, adapted for new DETAIL format

  Use :: nlte_type
  Use :: nlte_dim
  Use :: fund_const
  Use :: fastwind_params, Only: oppath, resol

  Use :: princesa_var, Only: labl, labl4, frecin, fl

  Use :: nlte_var, Only: modnam, ifre, fre, op_data

  Use :: nlte_opt, Only: opdebug

  Implicit None
! ..
! .. parameters ..
  Integer (i4b), Parameter :: nmax_low = 1000 ! for structured part
! INTEGER(I4B), PARAMETER :: NMAX_HIGH =200 ! for smooth high freq. tail
  Integer (i4b), Parameter :: nmax_high = 1000 ! for smooth high freq. tail


! ..
! .. scalar arguments ..

  Integer (i4b) :: ii, nrecord, iz
  Character (30) :: fileop
  Character (6) :: key

! ..
! .. local scalars

  Integer (i4b) :: estado, icode, nlab, j, k, ndatos, nout, l1, l2, l2start

  Real (dp) :: xnaux, eryd, q, q1, alpha

  Character (30) :: label, label1

  Character (80) :: dummychar

  Logical :: lab0 !                 indicates that alpha = 0 for a certain
! range close
! to the edge

! ..
! .. local arrays ..

  Real (dp), Dimension (:), Allocatable :: fgs, pics, aux2

  Real (dp), Dimension (nmax_low+nmax_high) :: eout, picsout

  Character (6), Dimension (30) :: names

  Data names/'H ', 'HE', 'LI', 'BE', 'B ', 'C ', 'N ', 'O ', 'F ', 'NE', 'NA', &
    'MG', 'AL', 'SI', 'P ', 'S ', 'CL', 'AR', 'K ', 'CA', 'SC', 'TI', 'V ', &
    'CR', 'MN', 'FE', 'CO', 'NI', 'CU', 'ZN'/

  Write (999, *) ' STOP: SUBR. OP_RBFSETUP1 SHOULD NOT BE CALLED'
  Stop ' SUBR. OP_RBFSETUP1 SHOULD NOT BE CALLED'

  Write (label, *) nrecord
  label = 'LABEL' // adjustl(label)

  eryd = 0.D0
  lab0 = .False.

  Open (771, File=oppath//trim(fileop), Status='OLD', Iostat=estado)

  If (estado/=0) Then
    Print *, 'FILE OP_RBF ', oppath // trim(fileop), ' NOT FOUND'
    Write (999, *) ' STOP!'
    Stop
  End If

  Read (771, *) icode, nlab
  j = icode/100
  If (trim(names(j))/=trim(key)) Then
    Write (999, *) ' STOP: OP_RBF FILE INCONSISTENT (ELEMENT)'
    Stop ' OP_RBF FILE INCONSISTENT (ELEMENT)'
  End If
  k = icode - j*100
  If (k+1/=iz) Then
    Write (999, *) ' STOP: OP_RBF FILE INCONSISTENT (ION)'
    Stop ' OP_RBF FILE INCONSISTENT (ION)'
  End If

  Do
    Read (771, *, End=100) label1, ndatos
    If (trim(label1)==trim(label)) Go To 110
    Read (771, *)(xnaux, j=1, 2*ndatos)
  End Do

100 Continue
  Write (999, *) ' STOP: LABEL NOT FOUND IN OP_RBF FILE'
  Stop ' LABEL NOT FOUND IN OP_RBF FILE'

110 Allocate (fgs(ndatos+1), pics(ndatos+1), aux2(ndatos+1)) ! one additional
! point may be
! required

  Read (771, *)(fgs(j), j=1, ndatos), (pics(j), j=1, ndatos)

  Close (771)

  If (opdebug) Then
    dummychar = trim(modnam) // '/opori-crossbf-' // trim(labl(labl4(ii))) // &
      '.ascii'
    Open (771, File=trim(dummychar), Status='UNKNOWN')
    Write (771, Fmt='(1X,A10,3(1X,G18.8))') trim(labl(labl4(ii))), eryd, &
      fl(labl4(ii))/clight, frecin(ii)
    Do k = 1, ndatos
      xnaux = fgs(k)/clight
      Write (771, Fmt='(1X,2(G18.8,1X))') xnaux, pics(k)
    End Do
    Close (771)
  End If

! Is there any offset in the (energy) level definition between the OP data
! and DETAIL?
! In contrast to OP_RBFSETUP, this is NOT checked here, since we
! do not have the OP energy info. Thus, we do not shift the energy,
! and "believe" that everything is OK
! (Norbert P. has promised that he has done the control by hand, &
! and shifted, if necessary).

! If such an info will become available, we suggest to proceed as in the
! alternative routine


  aux2 = fgs/clight

  If (aux2(ndatos)>frecin(ii)) Then
    If (opdebug) Then
      Print *, &
        ' OP --- rbf cross-section does not include the ionization edge'
      Print *, '   Level ', trim(labl(labl4(ii))), '  -----  ', label
      Print *, '   Level defined at ', fl(labl4(ii))/clight, ' cm^{-1}'
      Print *, '   Ionization edge located at ', frecin(ii), ' cm^{-1}'
      Print *, '   while the available range is ', aux2(ndatos), ' -- ', &
        aux2(1), ' cm^{-1} '
      Print *, '   Assumed that cross-section is low in missing area !!!! '
    End If
    lab0 = .True.
  End If

  If (lab0) Then !                  one additional point at low cross-section,
!   to prevent problems
!   in the photo integration routine (integral does not change)
    ndatos = ndatos + 1
    aux2(ndatos) = frecin(ii)
    fgs(ndatos) = frecin(ii)*clight
    pics(ndatos) = pics(ndatos-1)/1000.
  End If

! smoothing with Gaussian: energy units used!
! the logical variable controls output (no output if set to .false.)

  Call smooth(aux2, pics, ndatos, resol, nmax_low, nmax_high, frecin(ii), &
    eout, picsout, nout, .False.)

  If (eout(1)/=fgs(1)/clight) Then
    Write (999, *) ' STOP: error in fgs(1)' ! highest energy, MUST be
!   identical
    Stop ' error in fgs(1)' !       highest energy, MUST be identical
  End If
! print*,abs(1.-eout(nout)/(fgs(ndatos)/clight))
  If (abs(1.D0-eout(nout)/(fgs(ndatos)/clight))>5.D-14) Then
    Write (999, *) ' STOP: error in fgs(ndatos)' ! lowest energy, identical
!   within precision
    Stop ' error in fgs(ndatos)' !  lowest energy, identical within precision
  End If
  eout(nout) = fgs(ndatos)/clight ! to remain consistent

  Deallocate (fgs, pics, aux2)

! Once the cross-section has been read and convolved, &
! we map the cross-section on FASTWIND's frequency grid

! units: fre, frecin, eout in cm^-1; xnaux, fgs Hz
! energy grid fre  ordered from red to blue
! energy grid eout ordered from blue to red

  l2start = nout - 1

  Do l1 = 1, ifre
    xnaux = fre(l1)
    If (fre(l1)<=eout(1) .And. fre(l1)>=frecin(ii)) Then
!     emax > fre(l1) > edge, starting close to edge
!     we do not allow lower frequencies than

      Do l2 = l2start, 1, -1 !      the edge one
        If (xnaux<=eout(l2) .And. xnaux>=eout(l2+1)) Then
          j = l2
          k = l2 + 1
!         can happen under bad conditions (high resol) due to numerical
!         problems of FFT
          If (picsout(j)<=0.D0) Then
            Write (999, *) ' STOP: sout(j) le 0.'
            Stop ' sout(j) le 0.'
          End If
          If (picsout(k)<=0.D0) Then
            Write (999, *) ' STOP: sout(k) le 0.'
            Stop ' sout(k) le 0.'
          End If
          Go To 120
        End If
      End Do
!     IF(LAB0) THEN ! in this case, everything should be OK
!     ALPHA=0.
!     GOTO 60
!     ENDIF
      Write (999, *) ' STOP: fre not found in eout'
      Stop ' fre not found in eout'
!     logarithmic interpolation
120   l2start = l2
      q = log10(xnaux/eout(k))/log10(eout(j)/eout(k))
      q1 = 1.D0 - q
      alpha = q*log10(picsout(j)) + q1*log10(picsout(k))
      alpha = 10.D0**alpha

      If (alpha<=0.D0) Then
        Print *, 'ERROR IN ', trim(oppath//trim(fileop)), &
          ' alpha zero or negative'
        Print *, '   label ', label, '  level ', trim(labl(labl4(ii)))

        dummychar = trim(modnam) // '/crossbf-error-' // &
          trim(labl(labl4(ii))) // '.dat'
        Open (77, File=trim(dummychar), Status='UNKNOWN')
        Do k = 1, ifre
          Write (77, *) fre(k), op_data(ii, k)
        End Do
        Close (77)
        Write (999, *) ' STOP!'
        Stop
      End If

      op_data(ii, l1) = alpha

    Else If (xnaux>eout(1)) Then !  nu**3 extrapolation for high energies,
!     just in case
      op_data(ii, l1) = picsout(1)*(eout(1)/xnaux)**3

    Else !                          just in case, we will never reach this
!     frequency
      op_data(ii, l1) = 0.D0
    End If

  End Do

  If (opdebug) Then
    dummychar = trim(modnam) // '/fstwnd-crossbf-' // trim(labl(labl4(ii))) // &
      '.ascii'
    Open (771, File=trim(dummychar), Status='UNKNOWN')
    Write (771, Fmt='(1X,A10,2(1X,G18.8))') trim(labl(labl4(ii))), &
      fl(labl4(ii))/clight, frecin(ii)
    Do k = 1, ifre
      xnaux = fre(k)*clight
      Write (771, Fmt='(1X,2(G18.8,1X))') fre(k), op_data(ii, k)
    End Do
    Close (771)
  End If

  Print *, ' PHOTO-CROSS SECTION FOR LEVEL', trim(labl(labl4(ii))), &
    ' CONVOLVED WITH GAUSSIAN'
  Return
End Subroutine

!-----------------------------------------------------------------

Subroutine smooth(e, s, n, resol, nmax_low, nmax_high, edge, eout, sout, nout, &
  message)

  Use :: nlte_type
  Use :: fund_const

  Implicit None

  Integer (i4b), Intent (In) :: n, nmax_low, nmax_high
  Real (dp), Dimension (n), Intent (In) :: e, s
  Real (dp), Intent (In) :: resol, edge
  Logical :: message, low

  Real (dp), Dimension (nmax_low+nmax_high), Intent (Out) :: eout, sout

  Integer (i4b), Intent (Out) :: nout

  Integer (i4b) :: nmax, i, istart

  Real (dp) :: dsx, rdst_max, ratio

! integer(i4b) :: index_miguel !miguel2023

  nmax = nmax_low + nmax_high
  eout = 0.
  sout = 0.

  If (4*resol>nmax_low) Then
    Write (999, *) ' STOP: nmax too low (change resol?)'
    Stop ' nmax too low (change resol?)'
  End If
  nout = 4*nint(resol)

! define smooth and non-smooth regime(s)
  rdst_max = 0.

! ordered from high to low energies
  If (e(1)<e(n)) Then
    Write (999, *) ' STOP: wrong order of energy in crossbf file'
    Stop ' wrong order of energy in crossbf file'
  End If

  Do i = 1, n - 1
    dsx = (1.-e(i+1)/e(i))
    rdst_max = max(rdst_max, dsx)
  End Do

! JP (July 2013)
  low = (edge<10000.) !             for low ionisation stages with edge >
! 10000 A (here: 10000 Kayser = 10000 A)if
! (low) then
  If (low) Then
    ratio = 40.
  Else
    ratio = 20.
  End If
! print*,edge,'low=',low


  istart = 1

100 Do i = istart, n - 1
    dsx = (1.-e(i+1)/e(i))
!   print*,istart,i,dsx,rdst_max,e(i),s(i)
    If (dsx<0.05*rdst_max) Exit
  End Do

  i = max(i, 2) !                   in case i=1
! changed from 10 to 20 in vers. 9.2.5, to allow for Si2 cross-sect.
! additional change by JP (July 2013) to allow for low-enery edges (e.g., MgI)
  If (e(i-1)/max(e(n),edge)>ratio) Then ! inclusion of e(i) to obtain better
!   overlap
!   but might be "normal" edge
!   "real" edge (which differs from e(n), particularly if ionisation to
!   excited state) included
    istart = i + 1
    If (istart>=n-1) Then
      Write (999, *) &
        ' STOP: resonance part too broad, further reduction required'
      Stop ' resonance part too broad, further reduction required'
    End If
    Go To 100
!   JO might be commented out for tests. Think about ratio
  End If

! print*,i,istart,n,e(i),dsx,edge,e(n)

! in case that everything smooth
  If (i==n) Then
    If (message) Then
      Print *
      Print *, ' smooth distribution from', e(1), 'to', e(i)
      Print *, ' no convolution required'
      Print *
    End If
    eout(1:i) = e
    sout(1:i) = s
    nout = n
    Return
  End If

  If (message) Then
    Print *
    Print *, ' smooth energy tail from', e(1), 'to', e(i)
    Print *, '     resonance part from', e(i+1), 'to', e(n)
    Print *
  End If

  Call convol(e(i-1:n), s(i-1:n), n+2-i, eout, sout, nout, resol, message, &
    low)

! add high freq. tail
  If (sout(1)/=s(i-1)) Then
    Write (999, *) ' STOP: error in signal at transition freq.'
    Stop ' error in signal at transition freq.'
  End If

  nout = nout + i - 2
  If (nout>nmax) Then
    Print *, nmax_low, nmax_high
    Print *, nout, nmax
    Write (999, *) ' STOP: nout > nmax'
    Stop ' nout > nmax'
  End If

  eout = cshift(eout, -(i-2))
  sout = cshift(sout, -(i-2))
  eout(1:i-2) = e(1:i-2)
  sout(1:i-2) = s(1:i-2)

! for tests
! open(1,file='cross_out')
! do i=1,nout
! write(1,*),eout(i),sout(i)
! enddo
! close(1)

  Return

End Subroutine

!-----------------------------------------------------------------

Subroutine convol(e, s, n, eout, sout, nout, resol, message, low)

! convolution with gaussian of width resol (e/de, FWHM, see below)
! (assumed at the center of the given range)

! always resample input, since padding on both sides necessary

! n is original dimension (e,s)
! n1 refers to orginal plus padded data

! convolved signal s on modified, equidistant grid e is returned as sout(eout)
! (typically, 4 points per average resol. element)

  Use :: nlte_type
  Use :: fund_const
  Use :: fastwind_params, Only: rsmpl_fac, max_resol

  Implicit None

  Integer (i4b), Parameter :: nwmax = 6 ! maximum width of gaussian, in
! doppler units

  Integer (i4b) :: n, nout
  Real (dp), Dimension (n), Intent (In) :: e, s
  Real (dp), Dimension (nout), Intent (Out) :: eout, sout
  Real (dp), Allocatable, Dimension (:) :: x, y, sxx, syy, px, py, logx, logy
  Real (dp) :: resol
  Logical :: message, swap, low

  Integer (i4b) :: n1, nn, pow, k, i, j, kk, min_cnr, max_cnr, cnr
  Real (dp) :: s1, sn, x1, xn, rdst, rdst_min, rdst_max, dsx, meanw, eps, &
    diff, xpow, q, q1, fac, emean, vdop, dedop, percent, pxunit, errmax


  If (max_resol<10.*resol) Then
    Write (999, *) ' STOP: increase max_resol'
    Stop ' increase max_resol'
  End If
  n1 = n + 2*nwmax !                account for additional points for padding

  Allocate (x(n1), y(n1))

  If (e(1)>e(n)) Then !             change order (assumed from low to high)
    swap = .True.
    Do k = 1, n
      i = n + 1 - k
      x(i) = e(k)
      y(i) = s(k)
    End Do
  Else
    swap = .False.
    x(1:n) = e
    y(1:n) = s
  End If

  s1 = y(1)
  sn = y(n)
  x1 = x(1)
  xn = x(n)

  fac = 2.*sqrt(alog(2.))

! JP (July 2013)
! in most cases, resol assumed to be valid at center of range
  If (.Not. low) Then
    emean = 0.5*(x1+xn)
  Else
    emean = 0.75*x1 + 0.25*xn !     x1+0.25(xn-x1)
  End If

  vdop = clight/(resol*fac)
  dedop = vdop*(emean/clight)

! end padding, assuming symmetric response function (refering to no. of
! points)
  Do k = 1, nwmax
    x(n+k) = xn + k*dedop
    y(n+k) = sn !                   constant end
  End Do

  x = cshift(x, -nwmax) !           shift by nwmax to the right
  y = cshift(y, -nwmax)

! padding at begin, again assuming symmteric response
  Do k = 1, nwmax
    x(k) = x1 - (nwmax+1-k)*dedop
    y(k) = s1 !                     constant beginning
  End Do

  If (x(nwmax+1)/=x1) Then
    Write (999, *) ' STOP: error in x1'
    Stop ' error in x1'
  End If
  If (x(n+nwmax)/=xn) Then
    Write (999, *) ' STOP: error in xn'
    Stop ' error in xn'
  End If

  If (message) Then
    Print *, ' min(e) = ', x1, ' incl. padded reg.', x(1)
    Print *, ' max(e) = ', xn, ' incl. padded reg.', x(n1)
    Print *
  End If

  Allocate (logx(n1), logy(n1))
  logx = log10(x)
  logy = log10(y)

  rdst = (x(n1)-x(1))/(n1-1) !      resampling distance

! resample x-vector always

! dx from actual minimum * rsmpl_fac or max. resol

  rdst_min = 1.D10
  rdst_max = 0.
  Do k = 1, n1 - 1
    dsx = x(k+1) - x(k)
    rdst_min = min(rdst_min, dsx)
    rdst_max = max(rdst_max, dsx)
  End Do

! JP (July 2013)
  If (.Not. low) Then
    meanw = 0.5*(x(1)+x(n1))
  Else
    meanw = 0.75*x(1) + 0.25*x(n1) ! x1+0.25(xn1-x1)
  End If

  eps = meanw/max_resol
  diff = rdst_max - rdst_min

  If (message) Then
    Print *, ' Max(delta dE): ', diff
    Print *
    Print *, ' RDMIN = ', rdst_min, ', RDMAX = ', rdst_max, ' EPS = ', eps
    Print *
  End If

  rdst = max(rdst_min*rsmpl_fac, eps)
  nn = (x(n1)-x(1))/rdst + 1

  xpow = log(dble(nn))/log(2.D0)
  pow = int(xpow+1.E-15)

! print*, xpow,pow,nn,low

  If (pow<=2) Then
    Write (999, *) ' STOP: too few points in convol!'
    Stop ' too few points in convol!'
  End If

  If (.Not. low .And. pow>15) Then
    If (message) Then
      Print *, ' pow = ', pow, ': too many points in convol!'
      Print *, ' reset to pow =15'
    End If
    pow = 15
  End If

  If (low .And. pow>16) Then
    If (message) Then
      Print *, ' pow = ', pow, ': too many points in convol!'
      Print *, ' reset to pow =16'
    End If
    pow = 16
  End If

  pow = pow + 1
  nn = 2**pow

  rdst = (x(n1)-x(1))/(nn-1)

  If (rdst>rdst_min .And. message) Then
    Print *, ' WARNING WARNING WARNING WARNING WARNING WARNING WARNING!'
    Print *, ' WARNING WARNING WARNING WARNING WARNING WARNING WARNING!'
    Print *, ' WARNING WARNING WARNING WARNING WARNING WARNING WARNING!'
    Print *
    Print *, &
      ' rdst > rdst_min: extreme narrow spikes may be not well sampled!'
    Print *
    Write (999, *) ' WARNING WARNING WARNING WARNING WARNING WARNING WARNING!'
    Write (999, *) ' WARNING WARNING WARNING WARNING WARNING WARNING WARNING!'
    Write (999, *) ' WARNING WARNING WARNING WARNING WARNING WARNING WARNING!'
    Write (999, *)
    Write (999, *) &
      ' rdst > rdst_min: extreme narrow spikes may be not well sampled!'
    Write (999, *)
  End If

  Allocate (sxx(nn), syy(nn))

  Do k = 0, nn - 1
    sxx(k+1) = x(1) + k*rdst
  End Do

  If (abs(1.D0-sxx(1)/x(1))>1.D-14) Then
    Write (999, *) ' STOP: error in sxx(1)'
    Stop ' error in sxx(1)'
  End If
  If (abs(1.D0-sxx(nn)/x(n1))>1.D-14) Then
    Write (999, *) ' STOP: error in sxx(nn)'
    Stop ' error in sxx(nn)'
  End If

! first part, padded
  Do j = 1, nn
    If (sxx(j)>=x1) Exit
    syy(j) = s1
  End Do

  kk = 1

! interpolation of orginal signal onto fine grid
  Do k = j, nn - 1

    If (sxx(k)>xn) Exit !           padded region

    Do i = kk, n1 - 1
      If (sxx(k)>x(i) .And. sxx(k)<=x(i+1)) Go To 100
    End Do
    Write (999, *) ' STOP: sxx(k) not found in x'
    Stop ' sxx(k) not found in x'
100 kk = i

!   logarithmic interpolation
    q = (log10(sxx(k))-logx(i+1))/(logx(i)-logx(i+1))
    q1 = 1.D0 - q
    syy(k) = q*logy(i) + q1*logy(i+1)
    syy(k) = 10.**syy(k)
  End Do

! padded end region
  Do i = k, nn
    syy(i) = sn
  End Do

  If (message) Then
    Print *, ' resampled grid with resampling distance: ', rdst
    Print *, ' resampled spectrum has ', nn, ' points'
    Print *
  End If

! open(1,file='resampled')
! do i=1,nn
! write(1,*),sxx(i),syy(i)
! enddo
! close(1)


! Gauss profile

  min_cnr = 5
  max_cnr = nn*0.8

  If (message) Then
    Write (*, Fmt='('' Vdop='',F7.0,'' km/s  <=>'',F7.0,'' Kayser'')') vdop/ &
      1.E5, dedop
    Write (*, Fmt='('' FWHM='',F7.0,'' km/s  <=>'',F7.0,'' Kayser'')') fac* &
      vdop/1.E5, fac*dedop
    Print *
  End If

  cnr = (4.*dedop/rdst) + 1

  If (max_cnr<=cnr) Then
    percent = cnr/float(nn)*100.
    Write (*, Fmt='('' Gauss profile is very broad:'', i2,     &
      &     '' % of spectrum points are forming the profile'')') percent
  End If

  If (min_cnr>=cnr) Then
    Write (*, Fmt='('' Gauss profile is very narrow: only '', i2, &
      &         '' data points are forming the profile'')') cnr
  End If

  Allocate (px(nn), py(nn))

  pxunit = rdst/dedop
  Do i = 1, nn
    k = i - 1 - nn/2
    px(i) = k*pxunit
  End Do

  py = 0.

  Where (abs(px)<=float(nwmax))
    py = exp(-(px**2))/(sqrt(pi)*dedop)
  End Where

  py = cshift(py, -nn/2) !          wrap around zero


! check whether zero padding was successful:
! on the left  side, signal has to be e1 where py has finite values
! on the right side, signal has to be en where py has finite values

  Do k = 1, nn
    If (py(k)==0.D0) Exit
    If (syy(k)/=s1) Then
      Write (999, *) ' STOP: error in padding on right side'
      Stop ' error in padding on right side'
    End If
  End Do

  Do k = nn - 1, 1, -1
    If (py(k)==0.D0) Exit
    If (syy(k)/=sn) Then
      Write (999, *) ' STOP: error in padding on right side'
      Stop ' error in padding on right side'
    End If
  End Do

! for tests
! do k=1,nn
! print*,k,sxx(k),syy(k),py(k)
! enddo

  errmax = maxval(syy(1:nn))/minval(syy(1:nn))
  If (message) Then
    Print *, ' ratio largest/lowest cross-section in res. part: ', errmax
    Print *
  End If

! Jo 2023: for large differences in amplitudes, use FFT in quadruple precision
! currently (from experience with CII), we use a ratio of 1.d8.
! If something is rotten, use a lower value

  If (errmax<1.D8) Then
    Call conv_spec(sxx, syy, py, nn)
  Else
    Call conv_spec_quad(sxx, syy, py, nn)
  End If

! do k=1,nn
! print*,k,sxx(k),syy(k)
! enddo

  Deallocate (px, py)

! interpolation onto appropriate freq. grid
! (4 points per average resol. element)
  rdst = (xn-x1)/(nout-1)

  If (swap) Then !                  change order

    Do k = 1, nout
      eout(k) = xn - (k-1)*rdst
    End Do

    kk = nn

    Do k = 1, nout
      Do i = kk, 2, -1
        If (eout(k)<=sxx(i) .And. eout(k)>=sxx(i-1)) Go To 110
      End Do
      Write (999, *) ' STOP: eout(k) not found in sxx'
      Stop ' eout(k) not found in sxx'
110   kk = i

!     linear interpolation, since sxx well resolved
      q = (eout(k)-sxx(i-1))/(sxx(i)-sxx(i-1))
      q1 = 1.D0 - q
      sout(k) = q*syy(i) + q1*syy(i-1)
!     print*,k,eout(k),sout(k)
    End Do

!   check consistency
!   JP (July 2013)
    If (.Not. low) Then
      errmax = 2.
    Else
      errmax = 10.
    End If

    If (abs(1.-sout(1)/sn)>0.1) Then
      If (message) Then
        Print *, ' WARNING!!!! Problem at begin of high energy tail', sout(1), &
          sn
        Write (999, *) ' WARNING!!!! Problem at begin of high energy tail', &
          sout(1), sn
      End If
!     JO Sept. 2016 commented out
!     if (abs(1.-sout(1)/sn).gt.errmax)  stop ' error at begin of high energy
!     tail (swap)'
!     might need to be commented out
    End If

    sout(1) = sn
    sout(nout) = s1 !               to obtain similar alpha0
!   (sout(nout) might not be equal to s1 due to convolution)

  Else

    Do k = 1, nout
      eout(k) = x1 + (k-1)*rdst
    End Do

    kk = 1

    Do k = 1, nout
      Do i = kk, nn - 1
        If (eout(k)>sxx(i) .And. eout(k)<=sxx(i+1)) Go To 120
      End Do
      Write (999, *) ' STOP: eout(k) not found in sxx'
      Stop ' eout(k) not found in sxx'
120   kk = i

!     linear interpolation, since sxx well resolved
      q = (eout(k)-sxx(i+1))/(sxx(i)-sxx(i+1))
      q1 = 1.D0 - q
      sout(k) = q*syy(i) + q1*syy(i+1)
    End Do

!   check consistency
    If (abs(1.-sout(nout)/sn)>0.1) Then
      If (message) Then
        Print *, ' WARNING!!!! Problem at begin of high energy tail', &
          sout(nout), sn
        Write (999, *) ' WARNING!!!! Problem at begin of high energy tail', &
          sout(nout), sn
      End If
!     JO Sept. 2016 commented out
!     if (abs(1.-sout(nout)/sn).gt.2.)  stop ' error at begin of high energy
!     tail (noswap)'
!     might need to be commented out
    End If

    sout(nout) = sn
    sout(1) = s1 !                  to obtain identical alpha0
!   (sout(1) might not be equal to s1 due to convolution)

  End If

! for tests
! open(1,file='cross_out1')
! do i=1,nn
! if(sxx(i).lt.x1) cycle
! if(sxx(i).gt.xn) exit
! write(1,*),sxx(i),syy(i)
! enddo
! close(1)
  Deallocate (x, y, logx, logy, sxx, syy)

  Return
End Subroutine

!------------------------------------------------------------------------

Subroutine conv_spec(sx, sy, py, n)
! double precision

! note1:
! in contrast to the idl routine (convol.pro), we wrap the response function
! around zero, to enable padding in a simple way.
! This is necessary here since we deal with functions which are completely
! non-periodical (s(1) and s(N) have completely different orders),
! in contrast to convolving rectified spectra which begin and end at unity.

! note2:
! this routine uses the "slow" way. By exploiting the fact that
! signal and response are real functions, one could save a factor of
! roughly two. Since the convolution(s) have to be done only once per
! frequ. grid, we use here the "standard" way for clarity.

  Use :: nlte_type
  Use :: nlte_dim
  Use :: fund_const

  Implicit None
  Integer (i4b) :: n
  Real (dp), Dimension (n) :: sx, sy, py
  Complex (dp), Allocatable, Dimension (:) :: sy_t, py_t

  Allocate (sy_t(n), py_t(n))

  sy_t = sy
  py_t = py

  Call fft(sy_t, 1, n) !            forward
  Call fft(py_t, 1, n) !            forward
  sy_t = sy_t*py_t !                convolution
  Call fft(sy_t, -1, n) !           backwards
  sy = dble(sy_t) !                 real part
  sy = ((sx(n)-sx(1))/(n-1))*sy !   Physical domain (hn = Hn*Delta)

  Deallocate (sy_t, py_t)

  Return

End Subroutine

!------------------------------------------------------------------------

Subroutine conv_spec_quad(sx, sy, py, n)
! quadruple precision

! note1:
! in contrast to the idl routine (convol.pro), we wrap the response function
! around zero, to enable padding in a simple way.
! This is necessary here since we deal with functions which are completely
! non-periodical (s(1) and s(N) have completely different orders),
! in contrast to convolving rectified spectra which begin and end at unity.

! note2:
! this routine uses the "slow" way. By exploiting the fact that
! signal and response are real functions, one could save a factor of
! roughly two. Since the convolution(s) have to be done only once per
! frequ. grid, we use here the "standard" way for clarity.

  Use :: nlte_type
  Use :: nlte_dim
  Use :: fund_const

  Implicit None
  Integer (i4b) :: n
  Real (dp), Dimension (n) :: sx, sy, py
  Complex (qp), Allocatable, Dimension (:) :: sy_t, py_t

  Allocate (sy_t(n), py_t(n))

  sy_t = sy
  py_t = py

  Call fft_quad(sy_t, 1, n) !       forward
  Call fft_quad(py_t, 1, n) !       forward
  sy_t = sy_t*py_t !                convolution
  Call fft_quad(sy_t, -1, n) !      backwards
  sy = dble(sy_t) !                 real part
  sy = ((sx(n)-sx(1))/(n-1))*sy !   Physical domain (hn = Hn*Delta)

  Deallocate (sy_t, py_t)

  Return

End Subroutine

!------------------------------------------------------------------------

Subroutine fft(x, dir, length)
! fast-fourier transformation. length has to be a power of 2.
  Use :: nlte_type
  Use :: fund_const, Only: pi
  Implicit None

  Integer (i4b), Intent (In) :: dir, length
  Complex (dp), Dimension (0:length-1) :: x

  Complex (dp) :: wp, w, ws, temp
  Real (dp) :: theta, pi1
  Integer (i4b) :: nob, i, j, mmax, istep, m, bit_reverse

  nob = nint(log(real(length))/log(2.0))
  If (2**nob/=length) Then
    Write (999, *) ' STOP: input data length must be a power of 2'
    Stop ' input data length must be a power of 2'
  End If

  Do i = 1, length - 2
    j = bit_reverse(i, nob)
    If (j>i) Call swap(x(i), x(j))
  End Do

  If (dir==1) Then
    pi1 = -pi
  Else If (dir==-1) Then
    pi1 = pi
  Else
    Write (999, *) ' STOP: wrong dir in fft'
    Stop ' wrong dir in fft'
  End If

  mmax = 1
  Do
    If (length<=mmax) Exit
    istep = 2*mmax

    theta = pi1/mmax

    wp = cmplx(-2.0D0*sin(.5D0*theta)**2, sin(theta), kind=dp)
    w = cmplx(1.D0, 0.D0, kind=dp)

    Do m = 1, mmax
      ws = w
      Do i = m - 1, length - 1, istep
        j = i + mmax
        temp = ws*x(j)
!       if(dir.eq.-1) print*,mmax,i,j,ws,wp,temp,x(i),x(j)
        x(j) = x(i) - temp
        x(i) = x(i) + temp
!       if(dir.eq.-1) print*,mmax,i,j,x(i),x(j)
      End Do
      w = w + w*wp
    End Do
    mmax = istep
  End Do

  If (dir==-1) x = x/length

End Subroutine

!------------------------------------------------------------------------

Subroutine fft_quad(x, dir, length)
! fast-fourier transformation. length has to be a power of 2.
! quadruple precision
  Use :: nlte_type
  Use :: fund_const, Only: pi => piqp
  Implicit None

  Integer (i4b), Intent (In) :: dir, length
  Complex (qp), Dimension (0:length-1) :: x

  Complex (qp) :: wp, w, ws, temp
  Real (qp) :: theta, pi1
  Integer (i4b) :: nob, i, j, mmax, istep, m, bit_reverse

  nob = nint(log(real(length))/log(2.0))
  If (2**nob/=length) Then
    Write (999, *) ' STOP: input data length must be a power of 2'
    Stop ' input data length must be a power of 2'
  End If

  Do i = 1, length - 2
    j = bit_reverse(i, nob)
    If (j>i) Call swap_quad(x(i), x(j))
  End Do

  If (dir==1) Then
    pi1 = -pi
  Else If (dir==-1) Then
    pi1 = pi
  Else
    Write (999, *) ' STOP: wrong dir in fft'
    Stop ' wrong dir in fft'
  End If

  mmax = 1
  Do
    If (length<=mmax) Exit
    istep = 2*mmax

    theta = pi1/mmax

    wp = cmplx(-2.0_qp*sin(.5_qp*theta)**2, sin(theta), kind=qp)
    w = cmplx(1.0_qp, 0._qp, kind=qp)

    Do m = 1, mmax
      ws = w
      Do i = m - 1, length - 1, istep
        j = i + mmax
        temp = ws*x(j)
!       if(dir.eq.-1) print*,mmax,i,j,ws,wp,temp,x(i),x(j)
        x(j) = x(i) - temp
        x(i) = x(i) + temp
!       if(dir.eq.-1) print*,mmax,i,j,x(i),x(j)
      End Do
      w = w + w*wp
    End Do
    mmax = istep
  End Do

  If (dir==-1) x = x/length

End Subroutine

!------------------------------------------------------------------------

Function bit_reverse(i, size)
  Use :: nlte_type
  Implicit None
  Integer (i4b) :: i, size, bit_reverse
  Integer (i4b) :: length, temp

  temp = i
  Do length = size, 2, -1
    temp = ishftc(temp, -1, length)
  End Do

  bit_reverse = temp
End Function

!------------------------------------------------------------------------

Subroutine swap(a, b)
  Use :: nlte_type
  Implicit None
  Complex (dp) :: a, b, temp

  temp = a
  a = b
  b = temp
End Subroutine

!------------------------------------------------------------------------

Subroutine swap_quad(a, b)
! quadruple precision
  Use :: nlte_type
  Implicit None
  Complex (qp) :: a, b, temp

  temp = a
  a = b
  b = temp
End Subroutine

!------------------------------------------------------------------------

Subroutine ilowimax_lte(ilow, imax, te, teff, nd, ntemp, helone, optcmf1, &
  optout, lte_update)

! definition of ilow and imax for LTE calculations (first iteration cycle)
! and H/He always (in the old spirit, which worked)

! most of the input parameters refer to the older version, and might
! come into handy in future again

  Use :: nlte_type
  Use :: nlte_dim
  Use :: princesa_var, Only: nat, labat, nions, le, li, ifirsl

  Use :: nlte_var, Only: mainion, ionis, enionnd, imia, imaa, blevel, alevel

  Implicit None

! .. parameters ..
  Integer (i4b), Parameter :: kel = id_atoms
  Integer (i4b), Parameter :: nrec = id_llevs
  Integer (i4b), Parameter :: nd1 = id_ndept
! ..
! .. scalar arguments ..
  Integer (i4b) :: nd, ntemp
  Real (dp) :: teff
  Logical :: helone, optcmf1, optout, lte_update
! ..
! .. array arguments ..
  Real (dp) :: te(nd1)
  Integer (i4b) :: ilow(nd1, kel), imax(nd1, kel)

  Integer (i4b) :: k, l, i, main, ioni, imi, ima, i1, i2, deltai, nato, inato, &
    ij, nlas, igenio

  Character :: at*6

! all ionisation stages relative (1 is lowest ion in DETAIL)

  Do k = 1, nat
    Do l = ntemp, nd
      main = mainion(l, k)
      ioni = ionis(k, l)
      at = labat(k)

      If (at=='H') Then
        i1 = 1
        i2 = 1

      Else If (at=='HE') Then
!       exception for He, since option HELONE might be set to .false.
!       (for very hot stars)

        If (helone) Then !          this is the standard
          i1 = 1
        Else
          i1 = 2
        End If

        If (main>=2) Then
          If (main==3) Then
            Write (999, *) ' STOP: CHECK IONIZATION INDEXING!'
            Stop ' CHECK IONIZATION INDEXING!'
          End If
          i2 = 2
          If (enionnd(k,3,l)/enionnd(k,2,l)<1.D-12) i2 = 1
        Else If (enionnd(k,3,l)/enionnd(k,1,l)>=1.D-12) Then
          i2 = 2
        Else
          If (.Not. helone) Then
            Write (999, *) ' STOP: NO HE I, BUT DOMINANT ION!!!'
            Stop ' NO HE I, BUT DOMINANT ION!!!'
          End If
          i2 = 1 !                  for very cool stars
        End If
!       JO Nov. 2019: to avoid different ilow, imax for cool stars, following
!       hack
!       March 2025: include 15kK
        If (teff>=15000.) Then
          If (i2/=2) Then
            Print *, l, ' IMAX(HE) SET TO 2 (ILOWIMAX_LTE)'
!           force to include HeII
            i2 = 2
          End If
        Else
          If (.Not. helone) Then
            Write (999, *) ' STOP: HELONE = F FOR Teff < 15 kK (ILOWIMAX_LTE)'
            Stop ' HELONE = F FOR Teff < 15 kK (ILOWIMAX_LTE)'
          End If
          If (i2/=1) Then
!           changed Oct. 2021, to allow consistency with standard version
!           changed back, since otherwise SA treatment problematic
            If (i2/=2) Then
              Write (999, *) ' STOP: IMAX(HE) ne (1 or 2)'
              Stop ' IMAX(HE) ne (1 or 2)'
            End If
!           force to exclude HeII
            i2 = 1
            Print *, l, ' IMAX(HE) SET TO 1 (ILOWIMAX_LTE)'
!           PRINT*,L,' IMAX(HE) NE 1 (ILOWIMAX_LTE) for Teff lt 15 kK'
          End If
        End If
      Else !                        ALL OTHER ELEMENTS
        If (main==ioni) Then
          i1 = main - 2
          i2 = main
        Else
          i1 = main - 1
          i2 = main + 1
        End If
      End If

      ilow(l, k) = max(i1, 1)
      imax(l, k) = min(i2, ioni)
      deltai = imax(l, k) - ilow(l, k)
      If (deltai<0) Then
        Write (999, *) ' STOP: ERROR IN ILOW,IMAX'
        Stop ' ERROR IN ILOW,IMAX'
      End If
      If (labat(k)/='H' .And. labat(k)/='HE' .And. deltai<1) Then
        Print *, ' WARNING!!! K=', k, ' L=', l, ' ILOW = ', ilow(l, k), &
          ' IMAX = ', imax(l, k)
        Write (999, *) ' WARNING!!! K=', k, ' L=', l, ' ILOW = ', ilow(l, k), &
          ' IMAX = ', imax(l, k)
      End If
    End Do
  End Do

! ------------ determination of imia, imaa

  Do k = 1, nat
    imi = nions(k)
    ima = 0

    Do l = 1, nd
      imi = min0(ilow(l,k), imi)
      ima = max0(imax(l,k), ima)
    End Do

    If (ima==0) Then
      Write (999, *) ' STOP: SOMETHING WRONG WITH IMA'
      Stop ' SOMETHING WRONG WITH IMA'
    End If
    imia(k) = imi
    imaa(k) = ima
!   IF (OPTOUT .AND. ICONVER.EQ.1)
    Write (*, Fmt=100) k, imi, ima
    Print *
  End Do

  If (.Not. lte_update) Return

! from v10 on: check for new levels
! ONLY HERE, ILOW and IMAX can change regarding He.
! If they change to higher stages,
! subroutine OPACITC has no info how to deal with this.

! Example: in a certain iteration, ILOW=1, IMAX=1 at a certain L.
! after ILOWIMAX_LTE has been called (either due to restart, or because
! of temperature update, ILOW/IMAX might change to 2/2. Although the
! occupation
! numbers of the ground-level of i=2 is known, all other levels as well
! as the bf-transition remain unkown. THUS:


  Do l = 1, nd

    Read (17, Rec=l)(blevel(i), i=1, nrec)
    Read (15, Rec=l)(alevel(i), i=1, nrec)

    Do i = 1, nrec
      nato = le(i)
      If (labat(nato)/='HE') Cycle
      inato = li(i)
      ij = igenio(nato, imax(l,nato)) + 1
      nlas = ifirsl(ij)
!     slightly modified from V6.1 on

      If (blevel(i)==0. .And. ((inato>=ilow(l,nato) .And. &
        inato<=imax(l,nato)) .Or. i==nlas)) Then
!       VERY SMALL OCCUP. NUMBER
        blevel(i) = alevel(i)*1.D-10
        Print *, 'NEW LEVEL ', i, ' AT L = ', l, &
          ': OCC MODIFIED TO OCC-LTE*1.D-10'
!       print*,i,nato,inato,ilow(l,nato),imax(l,nato),blevel(i)
      End If
    End Do
!   so far, we assume that this happens only after T-update.
!   Thereafter, restart is invoked, which reads from file 27 (instead of 17)
!   If something goes wrong, try to uncomment the following line
!   WRITE (17,REC=L) (BLEVEL(I),I=1,NREC)
    Write (16, Rec=l)(blevel(i), i=1, nrec)
  End Do

  Return

100 Format (' ELEM. NO. ', I2, ' WITH ILOW= ', I1, ' AND IMAX= ', I1, &
    ' (LTE!!!, ZEFF CORRECTED)')
End Subroutine

!------------------------------------------------------------------------

Subroutine ilowimax_old(ilow, imax, te, teff, nd, ntemp, helone, optcmf1, &
  optout)

! definition of ilow and imax as done in versions until 8.7.1,
! kept for consistency to allow for a treatment of z=0 or comparison
! calculations

  Use :: nlte_type
  Use :: nlte_dim
  Use :: princesa_var, Only: nat, labat, nions

  Use :: nlte_var, Only: mainion, ionis, enionnd, imia, imaa, iconver

  Implicit None

! .. parameters ..
  Integer (i4b), Parameter :: kel = id_atoms
  Integer (i4b), Parameter :: nd1 = id_ndept
! ..
! .. scalar arguments ..
  Integer (i4b) :: nd, ntemp
  Real (dp) :: teff
  Logical :: helone, optcmf1, optout
! ..
! .. array arguments ..
  Real (dp) :: te(nd1)
  Integer (i4b) :: ilow(nd1, kel), imax(nd1, kel)

  Integer (i4b) :: k, l, main, ioni, i1, i2, imi, ima, imi1, ima1

  Real (dp) :: en, rel

  Character :: at*6

! all ionisation stages relative (1 is lowest ion in DETAIL)
! occmain maximum fraction of all ions EXCLUDING the highest one without
! transitions

  Do k = 1, nat

    Do l = ntemp, nd
      main = mainion(l, k)
      ioni = ionis(k, l)
      at = labat(k)
      If (at=='H') Then
        i1 = 1
        i2 = 1

      Else If (at=='HE') Then
        If (helone) Then
          i1 = 1
        Else
          i1 = 2
        End If

        If (main>=2) Then
          If (main==3) Then
            Write (999, *) ' STOP: CHECK IONIZATION INDEXING!'
            Stop ' CHECK IONIZATION INDEXING!'
          End If
          i2 = 2
          If (enionnd(k,3,l)/enionnd(k,2,l)<1.D-12) i2 = 1
        Else If (enionnd(k,3,l)/enionnd(k,1,l)>=1.D-12) Then
          i2 = 2
        Else
          If (.Not. helone) Then
            Write (999, *) ' STOP: NO HE I, BUT DOMINANT ION!!!'
            Stop ' NO HE I, BUT DOMINANT ION!!!'
          End If
          i2 = 1
        End If

!       ELSE IF (AT.EQ.'SI' .or. at.eq.'N') THEN
!       if old N-atom (only NII to NIV), use this block
      Else If (at=='SI') Then
!       NOTE THAT PRESENTLY ONLY SiII to SiIV are available from detail
!       thus, ilow,imax = 1 correspond to SiII, ilow,imax=3 to SiIV and so on
        If (at=='SI' .And. ioni/=3) Then
          Print *, ' EITHER: si atomic model has changed, update &
            &the following blocks'
          Print *, '     OR: IONIS (Si) ne 3'
          Write (999, *) ' STOP: PROBLEMS WITH HIGHEST ION FROM SILICON'
          Stop ' PROBLEMS WITH HIGHEST ION FROM SILICON'
        End If
        If (main==1 .Or. te(l)<=15000.) Then
          i1 = 1
          i2 = 2
        Else If (main==ioni) Then
          en = enionnd(k, main, l)
          rel = enionnd(k, main-2, l)/en
          If (rel>1.D-10) Then
            i1 = main - 2
            i2 = main
          Else
            i1 = main - 1
            i2 = main
          End If
        Else
          en = enionnd(k, main, l)
          rel = enionnd(k, main+2, l)/en
          If (rel>1.D-7) Then
            i1 = main - 1
            i2 = main + 1
          Else
            i1 = main - 1
            i2 = main
          End If
        End If

!       --- C, O, Mg for present Detail models (only TWO ions!)

      Else If (at=='C' .Or. at=='O' .Or. at=='MG') Then
        If (main==3) Then
          Write (999, *) ' STOP: CHECK IONIZATION INDEXING!'
          Stop ' CHECK IONIZATION INDEXING!'
        End If
        i1 = 1
        i2 = 2
        If (enionnd(k,1,l)/enionnd(k,2,l)<1.D-12) Then
          i1 = 2
          i2 = 2
        End If
        If (enionnd(k,2,l)/enionnd(k,1,l)<1.D-12) Then
          i1 = 1
          i2 = 1
        End If


      Else If (at=='N') Then
!       for new model atom from NII to NV
        If (main==ioni) Then
          i1 = main - 2
          i2 = main
        Else
          i1 = main - 1
          i2 = main + 1
        End If

      Else
        Write (999, *) ' STOP: ATOM NOT IMPLEMENTED YET'
        Stop ' ATOM NOT IMPLEMENTED YET'
      End If

      ilow(l, k) = max(i1, 1)
      imax(l, k) = min(i2, ioni)
    End Do

  End Do

! ------------ determination of imia, imaa

  Do k = 1, nat

    at = labat(k)

!   ------------ want to have at least SiII/III or SiIII/IV at
!   ALL depth points (minimum set) if CMF treatment aimed at

    If ((at=='SI' .Or. at=='N') .And. optcmf1) Then
      imi1 = maxval(ilow(:,k))
      ima1 = minval(imax(:,k))
      If ((ima1-imi1)<0) Then
        Write (999, *) ' STOP: Si: SOMETHING WRONG WITH IMI1,IMA1 IN OCCLTE'
        Stop ' Si: SOMETHING WRONG WITH IMI1,IMA1 IN OCCLTE'
      End If
      If ((ima1-imi1)==0) Then !    ONLY ONE ION PRESENT AT ALL DEPTH PONTS
        If (ima1==3) Then
          Do l = 1, nd
            If (imax(l,k)/=3) Then
              Write (999, *) ' STOP: Si: INCONSISTENCY IN IMA1'
              Stop ' Si: INCONSISTENCY IN IMA1'
            End If
            ilow(l, k) = min0(ilow(l,k), 2)
          End Do
        Else If (ima1==1) Then
          Do l = 1, nd
            If (ilow(l,k)/=1) Then
              Write (999, *) ' STOP: Si: INCONSISTENCY IN IMI1'
              Stop ' Si: INCONSISTENCY IN IMI1'
            End If
            imax(l, k) = max0(imax(l,k), 2)
          End Do
        Else
          If (ima1/=2 .Or. imi1/=2) Then
            Write (999, *) ' STOP: Si: INCONSISITENCY IN IMA1,IMI1'
            Stop ' Si: INCONSISITENCY IN IMA1,IMI1'
          End If
!         decide wether II/III or III/IV
          If (teff<=20000.) Then !  might be changed
            Do l = 1, nd
              ilow(l, k) = 1
            End Do
          Else
            Do l = 1, nd
              imax(l, k) = 3
            End Do
          End If
        End If

      End If
!     check once more
      imi1 = maxval(ilow(:,k))
      ima1 = minval(imax(:,k))
      If ((ima1-imi1)<=0) Then
        Write (999, *) ' STOP: Si: unsuccessful def. of ilow, imax'
        Stop ' Si: unsuccessful def. of ilow, imax'
      End If
      Print *
      Print *, ' MINIMUM IONIC SET FOR SILICON:', imi1, ' TO ', ima1
      Print *
    End If
!   ------------ Si treatment

    imi = nions(k)
    ima = 0

    Do l = 1, nd
      imi = min0(ilow(l,k), imi)
      ima = max0(imax(l,k), ima)
    End Do

    If (ima==0) Then
      Write (999, *) ' STOP: SOMETHING WRONG WITH IMA'
      Stop ' SOMETHING WRONG WITH IMA'
    End If
    imia(k) = imi
    imaa(k) = ima
    If (optout .And. iconver==1) Write (*, Fmt=100) k, imi, ima
    Print *
  End Do

  Return

100 Format (' ELEM. NO. ', I2, ' WITH ILOW= ', I1, ' AND IMAX= ', I1, &
    ' (LTE!!!, ZEFF CORRECTED)')
End Subroutine

!------------------------------------------------------------------------

Subroutine rbf_check

! checks for unique rbf transition frequencies and modifies them in case
  Use :: nlte_type
  Use :: nlte_dim
  Use :: princesa_var, Only: nf => id_rbftr, frecin, labl, labl4

  Implicit None

  Integer (i4b) :: i, j, ilog
  Real (dp) :: f1, f2

! first run: change
  Do i = 1, nf - 1
    f1 = frecin(i)
    Do j = i + 1, nf
      f2 = frecin(j)
      If (f2==f1) Then
!       change 7th digit (such a frequency is in the MHz domain and not
!       present here)
        ilog = log10(f2) - 7
        Print *, ' RBF transition frequencies changed, level ', labl(labl4(j))
        Print *, ' old frequency = ', f2
        f2 = f2 + 10.**ilog
        frecin(j) = f2
        Print *, ' new frequency = ', f2
        Print *
      End If
    End Do
  End Do

! 2nd run: check! due to bad luck, new conincidence might be created
  Do i = 1, nf - 1
    f1 = frecin(i)
    Do j = i + 1, nf
      f2 = frecin(j)
      If (f2==f1) Then
        Print *, ' still identical frequencies', f1
        Write (999, *) ' STOP: INDENTICAL RBF FREQUENCIES'
        Stop ' INDENTICAL RBF FREQUENCIES'
      End If
    End Do
  End Do

  Return
End Subroutine

!------------------------------------------------------------------------

Subroutine op_fit(fgs, pics, imax, ndatos)

! interpolation for original OPACITY data where sigma = 0.
! additionally, frequency grid is adapted

  Use :: nlte_type

  Implicit None
! .. scalar arguments ..
  Integer (i4b) :: ndatos, imax

! .. array arguments ..
  Real (dp), Dimension (imax) :: fgs, pics
! ..
! .. local scalars

  Integer (i4b) :: ia, ib, ia1, ib1, k
  Real (dp) :: x1, x2, y1, y2, dsde, df

! ....
! .. automatic array ..
  Real (dp), Dimension (ndatos) :: yop

  yop(1:ndatos) = pics(1:ndatos)

  If (yop(1)==0. .Or. yop(ndatos)==0.) Then
    Write (999, *) ' STOP: sig = 0 at boundaries'
    Stop ' sig = 0 at boundaries'
  End If

bigloop: Do

    Do k = 1, ndatos - 1
      If (yop(k)==0.) Then
        ia = k
        Go To 100
      End If
    End Do

    If (k/=ndatos) Then
      Write (999, *) ' STOP: ERROR1 IN OP-FIT'
      Stop ' ERROR1 IN OP-FIT'
    End If
    Exit bigloop

100 Do k = ia, ndatos
      If (yop(k)/=0.) Then
        ib = k - 1
        Go To 110
      End If
    End Do

    Write (999, *) ' STOP: OP_FIT: SHOULD NOT STOP HERE'
    Stop ' OP_FIT: SHOULD NOT STOP HERE'

!   LOG-LOG INTERPOLATION
110 ia1 = ia - 1
    ib1 = ib + 1

    x1 = fgs(ia1)
    x2 = fgs(ib1)

    y1 = yop(ia1)
    y2 = yop(ib1)

    dsde = log10(y2/y1)/log10(x2/x1)
    df = (x1-x2)/dble(ia1-ib1)

!   PRINT*,X1,X2,Y1,Y2
    Do k = ia, ib
      fgs(k) = fgs(ia1) + (k-ia+1)*df ! works for both increasing and
!     decreasing grid
      x2 = log10(y1) + dsde*log10(fgs(k)/x1)
      yop(k) = 10.D0**x2
!     PRINT*,K,FGS(K),YOP(K)
    End Do
!   PRINT*

  End Do bigloop

  pics(1:ndatos) = yop(1:ndatos)

  Return
End Subroutine

!------------------------------------------------------------------------

Subroutine modvl_update(nd, r, velo, dvdr, rho, xne, xnelte, xnh, temp, tau, &
  taur, clf, teff, ggrav, rstarnom, xmloss, beta, ndiv, yhein, yhe, hei, &
  rtau23, update_done)

! BASIC ASSUMPTION: NO CLUMPING IN PHOTOSPHERE (I.E., BELOW TRANSITION POINT)

! UPDATE OF PHOTOSPHERIC STRATIFICATION WITH GIVEN TEMP, CHIBAR_H AND NE
! CLUMPING FACTOR ASSUMED TO BE UNITY IN PHOTOSPHERIC REGION

! NOTE: CHIBAR_H calculated from NLTE (with XNE).
! JO, changed April 2016
! a) for standard calc (OPTCMF_FULL = .FALSE.)
! TAUR calculated from LTE (with XNELTE). Can lead to inconsistencies
! for cooler stars when H recombines
! b) for OPTCMF_FULL = .TRUE.
! TAUR calculated from NLTE (with XNE).

! JO: re-tested (to best of knowledge) April 2016

  Use :: nlte_type
  Use :: nlte_dim
  Use :: fund_const, Only: pi, xmh => amh, akb, rsu => rsun, &
    sigmae => sigmae_cross
  Use :: nlte_opt, Only: opt_photclump, optcmf_full

  Use :: princesa_var, Only: nat, labat, abund, weight

  Use :: nlte_var, Only: modnam, vmax, sr, srvmin, srnom, idum, gg => ggrav, &
    te => teff, a2te, delmu, no_check

  Use :: photstruc, Only: ndmod, xt, const, expo, expo_p, sigemin, xmext, &
    text, xm, chibar_h, delchi, pold => pkur, taurold => xnekur, &
    sig => gradkur

  Use :: nlte_porvor, Only: optthick

  Implicit None

! .. parameters ..

  Integer (i4b), Parameter :: nd1 = id_ndept
! ..
! .. scalar arguments ..
  Real (dp) :: teff, ggrav, xmloss, beta, yhe
  Real (dp), Intent (In) :: rstarnom, yhein, hei
  Real (dp), Intent (Inout) :: rtau23
  Integer (i4b), Intent (In) :: nd, ndiv
! ..
! .. array arguments ..
  Real (dp), Intent (In) :: taur(nd1)

  Real (dp), Intent (Inout) :: r(nd1), velo(nd1), dvdr(nd1), rho(nd1), &
    temp(nd1), clf(nd1), xne(nd1), xnelte(nd1), xnh(nd1), tau(nd1)
! ..
! .. local scalars ..
  Real (dp) :: vmin, b, h2, xm0, alpha, dm, r0, yhe1, sumatom, summass, xmu, &
    xmu1, vsonic, a2, xmmax, deltam, deltam1, x1, x2, p0, t0, rho0, taur0, pl, &
    tl, rl, taurl, hs, dmdr, srnew, rmax, dpx, sigem, taumin, summ, del, &
    dchidm, delchibar, chibar, average, xtin, safety, average_crit, rr1, rr2, &
    det, v1test, v2test, rr0, vsmooth, rhosmooth, dum


  Integer (i4b) :: l, ns, nene, isi, k, nok, nbad, ll, i, nsig
! ..
! .. local arrays ..
  Real (dp) :: radacc(nd1), p(nd1), tnew(nd1), y(4), rnew(nd1), rhonew(nd1), &
    v(nd1), gradv(nd1), rhocl(nd1), rhonewcl(nd1), xnecl(nd1), arho(nd1), &
    brho(nd1), crho(nd1), taurnew(nd1), rin(nd1), veloin(nd1), dvdrin(nd1), &
    rhoin(nd1), xneltein(nd1), xnhin(nd1), clfin(nd1), zz1(4), zz2(4), &
    vmat(4, 4), yy(4), gradvnew(nd1)


  Integer (i4b) :: index(nd1), indlud(4)

  Logical :: updated
  Logical, Intent (Out) :: update_done

  Character (7) :: char

  External :: derivs_upd, rkqc
! ..

! CHECK WHETHER UPDATE NEEDS TO BE PERFORMED
! ARHO: START VALUE OF RADACC
  Open (61, File=trim(modnam)//'/GRAD.OUT', Status='UNKNOWN')
  Do l = 1, nd
    Read (61, *) char, i, arho(l)
  End Do
  Close (61)

! CRHO: ACTUAL VALUE OF RADACC
  If (.Not. optcmf_full) Then
    Rewind (60)
    Do l = 1, nd
      Read (60, *) i, brho(l), crho(l)
    End Do
  Else
    Open (61, File=trim(modnam)//'/CHIBAR_H_CMF', Status='UNKNOWN')
    Rewind 61
    Do l = 1, nd
      Read (61, *) i, brho(l), crho(l)
    End Do
    Close (61)
  End If

! TEST CONSISTENCY
  Do l = 1, nd
    If (abs(1.-brho(l)/chibar_h(l))>1.D-13) Then
      Print *, brho(l), chibar_h(l)
      Write (999, *) ' STOP: INCONSISTENCY IN CHIBAR_H'
      Stop ' INCONSISTENCY IN CHIBAR_H'
    End If
  End Do

! JO Feb 2023: check whether close to Eddington limit
  dum = ggrav/(5.67D-5/2.9979D10*teff**4)
  Do l = ndiv, nd
    brho(l) = chibar_h(l)/dum
  End Do
  dum = maxval(brho(ndiv:nd))
  If (dum>=1.) Then
    Print *
    Print *, ' MAX. RATIO g_rad/g_grav (for L GE NDIV) = ', dum
    Write (999, *)
    Write (999, *) ' MAX. RATIO g_rad/g_grav (for L GE NDIV) = ', dum
    Write (999, *) ' STOP: LOCALLY BEYOND EDDINGTON LIMIT IN PHOTOSPHERE; &
      &INCREASE LOG G!!!!'
    Stop ' LOCALLY BEYOND EDDINGTON LIMIT IN PHOTOSPHERE; INCREASE LOG G!!!!'
  Else If (dum>0.9) Then
    Print *
    Print *, ' WARNING!!! WARNING!!! WARNING!!!'
    Print *, ' MAX. RATIO g_rad/g_grav (for L GE NDIV) = ', dum
    Print *, ' LOCALLY CLOSE (> 0.9) TO EDDINGTOT LIMIT'
    Print *, ' EITHER INCREASE LOG G OR PROCEED AT OWN RISK'
    Print *
    Write (999, *)
    Write (999, *) ' WARNING!!! WARNING!!! WARNING!!!'
    Write (999, *) ' MAX. RATIO g_rad/g_grav (for L GE NDIV) = ', dum
    Write (999, *) ' LOCALLY CLOSE (> 0.9) TO EDDINGTOT LIMIT'
    Write (999, *) ' EITHER INCREASE LOG G OR PROCEED AT OWN RISK'
    Write (999, *)
  End If

  Do l = 1, nd
    If (arho(l)/=0.D0) Exit
  End Do
  i = l

  If (i>=nd-5) Then
    Write (999, *) ' STOP: SOMETHING WRONG WITH RADACC_START NE 0'
    Stop ' SOMETHING WRONG WITH RADACC_START NE 0'
  End If

  average = 0.
  Do l = i, i + 4
    average = average + crho(l)/arho(l)
  End Do
  average = average/5.
  Print *, ' AVERAGE DEVIATION IN RADACC AT OUTER PHOTOSPHERE =', average
  Write (999, *) ' AVERAGE DEVIATION IN RADACC AT OUTER PHOTOSPHERE =', &
    average

! JO Jan 2016
! for standard path, update only done if pressure/density becomes
! lower than by very first solution. This is checked by AVERAGE_CRIT.
! IF AVERAGE_CRIT. > 1, this means that actual flux-mean opacities are
! larger than in first solution (based on rosseland mean); thus, rad.
! acceleration will be larger, and density lower.
! for OPTCMF_FULL, update is done later, and radacc "exactly" calculated
! thus update should be donw even if flux-mean opa has decreased.
  If (.Not. optcmf_full) Then
    average_crit = 1.
  Else
    average_crit = 0.8
  End If

  If (average<average_crit) Then
!   NO UPDATE PERFORMED
    Print *
    Print *, ' NO UPDATE OF PHOTOSPHERIC STRUCTURE PERFORMED: AVERAGE &
      &< AVERAGE_CRIT = ', average_crit
    Print *
    Write (999, *)
    Write (999, *) ' NO UPDATE OF PHOTOSPHERIC STRUCTURE PERFORMED: &
      &AVERAGE < AVERAGE_CRIT = ', average_crit
    Write (999, *)

    Rewind 20
    Read (20) teff, ggrav, sr, yhe, xmu, vmax, xmloss, beta, &
      (rin(i), i=1, nd), (veloin(i), i=1, nd), (dvdrin(i), i=1, nd), &
      (rhoin(i), i=1, nd), (xneltein(i), i=1, nd), (xnhin(i), i=1, nd), &
      (clfin(i), i=1, nd), (index(i), i=1, nd), (p(i), i=1, nd), ns, updated, &
      xtin

!   CHECK CONSISTENCY
    Do i = 1, nd
!     because of lazyness, no stop-output to LU 999
      If (r(i)/=rin(i)) Stop ' MODVL_UPDATE: R NE RIN'
      If (velo(i)/=veloin(i)) Stop ' MODVL_UPDATE: VELO NE VELOIN'
      If (dvdr(i)/=dvdrin(i)) Stop ' MODVL_UPDATE: DVDR NE DVDRIN'
      If (rho(i)/=rhoin(i)) Stop ' MODVL_UPDATE: RHO NE RHOIN'
      If (xnh(i)/=xnhin(i)) Stop ' MODVL_UPDATE: XNH NE XNHIN'
      If (clf(i)/=clfin(i)) Stop ' MODVL_UPDATE: CLF NE CLFIN'
    End Do

    If (xt/=xtin) Then
      Write (999, *) ' STOP: MODVL_UPDATE: XT NE XTIN'
      Stop ' MODVL_UPDATE: XT NE XTIN'
    End If


!   NOTE: XNELTEIN NE XNELTE!!!
!   XNELTEIN -> REWRITTEN TO MODEL
!   XNE, XNELTE (ARGUMENTS) REMAIN UNCHANGED

    updated = .True.
    update_done = .False.

    Rewind 20
    Write (20) teff, ggrav, sr, yhe, xmu, vmax, xmloss, beta, (r(i), i=1, nd), &
      (velo(i), i=1, nd), (dvdr(i), i=1, nd), (rho(i), i=1, nd), &
      (xneltein(i), i=1, nd), (xnh(i), i=1, nd), (clf(i), i=1, nd), &
      (index(i), i=1, nd), (p(i), i=1, nd), ns, updated, xt
    Return
  End If


  Print *
  Print *, ' =========== UPDATE OF PHOTOSPHERIC STRUCTURE ============'
  Print *
  Write (999, *)
  Write (999, *) ' =========== UPDATE OF PHOTOSPHERIC STRUCTURE ============'
  Write (999, *)


! check that no clumping except for tests

  If (.Not. opt_photclump) Then
    Do l = ndiv, nd
      If (clf(l)/=1.D0) Then
        Write (999, *) &
          ' STOP: PHOTOSPHERIC CLUMPING FOUND (ROUTINE MODVL_UPDATE)'
        Stop ' PHOTOSPHERIC CLUMPING FOUND (ROUTINE MODVL_UPDATE)'
      End If
    End Do
  End If

! transfer to module nlte_var (needed in DERIVS)
  srnom = rstarnom*rsu

  rhocl = rho*clf !                 OLD QUANTITY

  ndmod = nd !                      for module photstruc

  rnew = r*sr !                     in cgs, will contain updated radius
  v = velo*vmax !                   in cgs, will contain updated velocity
  rhonew = rho !                    will contain updated density
  tnew = temp !                     will contain updated temperature
  taurnew = taur !                  will contain updated taur (guess-value)
  gradv = dvdr*vmax/sr

  Allocate (xmext(ndmod), text(ndmod), pold(ndmod), taurold(ndmod), &
    sig(ndmod))
  xmext = 0. !                      standard grid (needed for dT/dm)
  text = temp

! integration assuming rho(r)=r^alpha*b
! => int = 1/(alpha+1)*(rho(r2)*r2-rho(r1)*r1)
! with
! alpha=delta log rho/delta log r and
! (previous version)
! b=rho(r)/r^alpha  => rho(x) = rho(r)*(x/r)^alpha with x on staggered grid
  Do l = 2, nd
    alpha = log10(rho(l-1)/rho(l))/log10(r(l-1)/r(l))
    xmext(l) = xmext(l-1) + (r(l-1)*rho(l-1)-r(l)*rho(l))/(alpha+1.) ! on
!   normal
!   grid
  End Do
  xmext = xmext*sr

! recalculate XM0 = M(NDIV) and check consistency
  vmin = v(ndiv)

  b = 1.D0 - (vmin/vmax)**(1.D0/beta)

  If (beta==1.D0) Then
    h2 = -log(1.D0-b)
  Else
    h2 = (1.D0-(vmin/vmax)**((1.D0-beta)/beta))/(1.D0-beta)
  End If

  r0 = rnew(ndiv)

  xm0 = xmloss/(4.D0*pi*r0*vmax)/b*h2

  dm = xm0/xmext(ndiv)

! IF(ABS(LOG10(DM)).GT.0.05) STOP ' PROBLEMS IN COLUMN DENSITY (ROUTINE
! MODVL_UPDATE)'
! JO Oct. 2017
  If (abs(log10(dm))>0.20) Then
    Do l = 1, ndiv
      Print *, rho(l), r(l), xmext(l)
    End Do
    Print *
    Print *, xm0, xmext(ndiv), dm, log10(dm)
    Write (999, *) ' STOP: PROBLEMS IN COLUMN DENSITY (ROUTINE MODVL_UPDATE)'
    Stop ' PROBLEMS IN COLUMN DENSITY (ROUTINE MODVL_UPDATE)'
  End If

  xmext = xmext*dm

! JO Sept 2021: next statements commented out, approach from version 10
! with modified precision used
! IF(ABS(XMEXT(NDIV)-XM0).GT.1.D-12) THEN
! PRINT*,XMEXT(NDIV),XM0
! STOP ' PRECISION PROBLEM: XMEXT(NDIV) NE XM0 (ROUTINE MODVL_UPDATE)'
! ENDIF
! JO Oct. 2017
! IF(ABS(1.-XMEXT(NDIV)/XM0).GT.1.D-15) THEN !version 10
  If (abs(1.-xmext(ndiv)/xm0)>1.D-12) Then ! modified precision
    Print *, xmext(ndiv), xm0
    Write (999, *) ' STOP: XMEXT(NDIV) NE XM0 (ROUTINE MODVL_UPDATE)'
    Stop ' XMEXT(NDIV) NE XM0 (ROUTINE MODVL_UPDATE)'
  End If
! JO Sept 2021: this is the new statement
  xmext(ndiv) = xm0 !               This is important for index-finding (IT)
! in routine DERIVE_UPD



  Print *, 'CHIBAR_H'
  Do l = 2, nd1 !                   XMEXT(1)=0.
    Print *, l, ' ', log10(xmext(l)), ' ', chibar_h(l)
  End Do

  ns = ndiv

  Print *
  Print *, ' DIVISION AT ', vmin*1.D-5, ' KM/S (L = ', ns, ')'
  Print *, ' CORRESPONDING TO ', r0/srnom, ' NOM. RADII'
  Print *, ' AND LG(M0) = ', log10(xm0)
  Print *

! to be consistent to routine MODVL
  nene = 0

  Do isi = 1, nat
    If (labat(isi)=='HE') nene = isi
  End Do

  If (nene==0) Then
    Print *, 'HE NOT FOUND - MODVL_UPDATE'
    If (hei/=0.D0) Then
      Print *, ' NO HELIUM BUT ELSCAT-CONTRIBUTION!!!'
      Print *, ' EITHER SET HEI TO ZERO OR CALCULATE WITH HELIUM!!!'
!     comment the following line out for specific tests
      Write (999, *) ' STOP: INCONSISTENT TREATMENT OF HELIUM'
      Stop ' INCONSISTENT TREATMENT OF HELIUM'
    End If
    yhe1 = yhein
  Else
    yhe1 = abund(nene)
    If (abs(1.D0-yhe1/yhein)>1.D-5) Then
      Write (999, *) ' STOP: ERROR IN YHE - MODVL_UPDATE'
      Stop ' ERROR IN YHE - MODVL_UPDATE'
    End If
  End If

  If (yhe1/=yhe) Then
    Write (999, *) ' STOP: ERROR IN YHE - MODVL_UPDATE'
    Stop ' ERROR IN YHE - MODVL_UPDATE'
  End If

! transfer to module nlte_var (needed in DERIVS)
  te = teff
  gg = ggrav

  sumatom = 0.D0
  summass = 0.D0

  Do k = 1, nat
    sumatom = sumatom + abund(k)
    summass = summass + weight(k)*abund(k)
  End Do

  If (nene==0 .And. yhe/=0.) Then
    sumatom = sumatom + yhe
    summass = summass + 4.D0*yhe
  End If

  summass = summass*xmh

  xmu = (1.D0+4.D0*yhe)/(2.D0+yhe*(1.D0+hei))

  Do l = 1, nd
!   clumping cancels out
    xmu1 = summass/(sumatom+xne(l)/xnh(l))/xmh

!   ---  delmu defined in order to correct a2te (vs(iso)**2/t)
!   ---  for ionisation effects

    delmu(l) = xmu/xmu1
!   PRINT*,L,DELMU(L)
  End Do

  xm = xmext !                      XM new grid, as calculated here

  vsonic = sqrt(akb*te/xmh/xmu)
  a2 = vsonic**2
  a2te = a2/te

  t0 = tnew(ns)
  rho0 = rhonew(ns)
  p0 = a2te*delmu(ns)*rho0*t0
  taur0 = taurnew(ns)

! P0 AND POLD (FROM 'MODEL') DIFFERENT, DUE TO DIFFERENT TEMPERATURE
! THUS, RECALCULATE P0LD


  Do l = ns, nd
    pold(l) = a2te*delmu(l)*rho(l)*temp(l)
!   JO changed April 2016: for OPTCMF_FULL, use NLTE values (both for ne and
!   taur)
    If (optcmf_full) Then
      sig(l) = xne(l)*sigmae/rhocl(l)
!     CHANGED FROM v10.1.3, SINCE USED FOR TAUR(=LTE)
    Else
      sig(l) = xnelte(l)*sigmae/rhocl(l)
    End If
  End Do
  taurold = taur

! JO April 2016
! For tests. Calculate such grad which is consistent with the present
! structure, and then do the update. In this case, the calculated grad
! has to be very similar with the (faked) input value.
! do l=ns+1,nd
! dpdm=(pold(l-1)-pold(l))/(xmext(l-1)-xm(l))
! grad=ggrav-(rnew(l)/srnom)**2*dpdm
! chibar_h(l)=grad/5.67d-5*2.9979d10/teff**4
! enddo
! print*,chibar_h

  If (.Not. optthick .And. .Not. optcmf_full) Then
!   chibar_h calculated with xne(NLTE)
    sigemin = minval(chibar_h)
!   IF (SIGEMIN.LT.MINVAL(SIG)) &
!   CHANGED FROM v10.1.3, SINCE NOW SIG REFERS TO LTE, WHILST SIGEMIN TO NLTE
    If (sigemin<minval(sig*xne/xnelte)) Then
      Write (999, *) &
        ' STOP: MIN(CHIBAR_H) < SIG(THOMSON). MAYBE ERROR IN WIND-OPACITIES'
      Stop ' MIN(CHIBAR_H) < SIG(THOMSON). MAYBE ERROR IN WIND-OPACITIES'
    End If
  Else
!   check only in unclumped/photospheric part (for OPTCMF_FULL, there can be a
!   narrow region outside NS with negative CHIBAR_H due to negative fluxes).
!   But note that later on SIG(L) will be used.
!   let's hope that there are no effects from wind conditions that influence
!   the region below NS
    sigemin = minval(chibar_h(ns:nd))
!   JO, April 2016: now, xne(NLTE) used, thus no correction
    If (optcmf_full) Then !         THICK OR NOT THICK
      If (sigemin<minval(sig(ns:nd))) Then
!       changed Oct. 2016; if transition point close to sonic point (vdiv >
!       0.1),
!       or negative fluxes close to NS (e.g., d10v), negative or low values
!       possible.
!       THUS
        nsig = min(ns+5, nd)
        sigemin = minval(chibar_h(nsig:nd))
        If (sigemin<minval(sig(ns:nd))) Then
          Write (999, *) ' STOP: MIN(CHIBAR_H) < SIG(THOMSON) BELOW NS+5'
          Stop ' MIN(CHIBAR_H) < SIG(THOMSON) BELOW NS+5'
        End If
      End If
    Else !                          NOT FULL
      If (sigemin<minval(sig*xne/xnelte)) Then
        Write (999, *) ' STOP: MIN(CHIBAR_H) < SIG(THOMSON) BELOW NS'
        Stop ' MIN(CHIBAR_H) < SIG(THOMSON) BELOW NS'
      End If
    End If
  End If

! JO Jan 2016: SIGEMIN redefined. Otherwise log10(CHI-SIGEMIN) (in
! CHIH_PARAMS)
! becomes infinite at that point where CHI = SIGEMIN
  If (optthick .Or. optcmf_full) Then
    If (optcmf_full) Then
      sigemin = minval(sig(ns:nd))
    Else
      sigemin = minval(sig(ns:nd)*xne(ns:nd)/xnelte(ns:nd))
    End If
  End If

  Do l = ns, nd
    If (abs(chibar_h(l)-sigemin)>0.05) Exit
  End Do
  nsig = max(ns, l-1)
! IN CASES WHERE ONLY THOMSON SCATTERING OR CONSTANT GRAD
  If (nd-nsig<=4) Then
    Do l = ns, nd
      If (abs(chibar_h(l)-sigemin)>0.01) Exit
    End Do
    nsig = max(ns, l-1)
    If (nd-nsig<=4) Then
      Write (999, *) ' STOP: TOO FEW POINTS FOR REGRESSION OF CHIBAR_H'
      Stop ' TOO FEW POINTS FOR REGRESSION OF CHIBAR_H'
    End If
  End If

! PARAMETERIZE CHIBAR_H IN THE RANGE NSIG...ND, &
! AND EVALUATE CHIBAR_H(FIT)/CHIBAR_H(ACTUAL) IN THE RANGE NS...ND
  If (.Not. optcmf_full) Then
    Call chih_params(chibar_h, pold, temp, delchi, sigemin, const, expo, &
      expo_p, nd, ns, nsig)
  Else
!   BECAUSE OF MORE 'DIFFICULT' RUN OF CHIBAR, MODIFIED CALCULATION OF
!   PARAMETERS, NOW INCLUDING SIGEMIN, REQUIRED
    Call chih_params1(chibar_h, pold, temp, delchi, sigemin, const, expo, &
      expo_p, nd, ns, nsig)
  End If

  xmmax = .995*xm(nd) !             SAFETY FACTOR, TO REMAIN INSIDE GRID

  deltam = log10(xmmax/xm0)/(nd-ns-5)

  If (abs(1.D0-xm(ns)/xm0)>1.D-14) Then
    Write (999, *) ' STOP: ERROR IN XM0 - MODVL_UPDATE'
    Stop ' ERROR IN XM0 - MODVL_UPDATE'
  End If
  xm(ns) = log10(xm0)

  radacc = 0.

  Do l = ns + 1, nd
    deltam1 = deltam
    If (l<=ns+8) deltam1 = deltam/2.D0
    If (l<=ns+4) deltam1 = deltam/3.D0
    If (l<=ns+2) deltam1 = deltam/6.D0
    xm(l) = xm(l-1) + deltam1
  End Do

  Do l = ns, nd
    xm(l) = 10.D0**xm(l)
  End Do
! JO Jan/April 2016: last point close to previous one, to improve diffusion
! approx.
  If (nd>=61) xm(nd) = xm(nd-1)*1.1

  x1 = xm0
  y(1) = p0
  y(2) = t0
  y(3) = r0
  y(4) = taur0
  p(ns) = p0
! TNEW(NS),RHONEW(NS) DO NOT CHANGE

! LOOP OVER ALL XM
! FROM V10.1 ON: WHILE LOOP, TO ALLOW FOR RESET OF L
! -------------------------------------------------------------------------
  l = ns + 1

  Do While (l<nd+1)

    pl = y(1)
    tl = y(2)
    rl = y(3)
    taurl = y(4)
    x2 = xm(l)
    hs = (x2-x1)/20.D0
    idum = l

!   new photosheric structure, &
!   using grad(m) via CHIBAR_H and approximate update of temperature
    Call odeint(y, 4, x1, x2, 1.D-5, hs, 0.D0, nok, nbad, derivs_upd, rkqc)

    p(l) = y(1)
    tnew(l) = y(2)
    rnew(l) = y(3)
    taurnew(l) = y(4)

!   one point just below SRNOM
    If (rnew(l)<srnom .And. rnew(l-1)>srnom) Then
      If (l<=ns+4) Then
        Print *
        Print *, ' WARNING! PHOTOSPHERE AT ', l, &
          ' VERY CLOSE TO DIVISION POINT AT', ns
        Print *
        Write (999, *)
        Write (999, *) ' WARNING! PHOTOSPHERE AT ', l, &
          ' VERY CLOSE TO DIVISION POINT AT', ns
        Write (999, *)
      End If
!     new try
      safety = 1.1D0
100   dmdr = (xm(l)-xm(l-1))/(rnew(l-1)-rnew(l))
      deltam = safety*dmdr*(rnew(l-1)-srnom) ! including safety factor
      If (deltam<0.) Then
        Write (999, *) ' STOP: ERROR IN DELTAM < 0'
        Stop ' ERROR IN DELTAM < 0'
      End If

!     FROM V10.1 ON: AVOID TOO SMALL SEPARATIONS (except L - 1 = NS, i.e.,
!     last point = starting point)
      If (deltam/xm(l-1)<0.1 .And. l-1/=ns) Then ! FACTOR 1.1
!       SKIP LAST POINT
        xm(l-1) = xm(l-1) + deltam
        l = l - 1
        x1 = xm(l-1)
        y(1) = p(l-1)
        y(2) = tnew(l-1)
        y(3) = rnew(l-1)
        y(4) = taurnew(l-1)
        idum = l
      Else
!       ACCEPT LAST POINT, RECALCULATE NEW POINT
        xm(l) = xm(l-1) + deltam
        y(1) = pl
        y(2) = tl
        y(3) = rl
        y(4) = taurl
      End If

      x2 = xm(l)
      hs = (x2-x1)/20.D0
      Call odeint(y, 4, x1, x2, 1.D-5, hs, 0.D0, nok, nbad, derivs_upd, rkqc)

      p(l) = y(1)
      tnew(l) = y(2)
      rnew(l) = y(3)
      taurnew(l) = y(4)

!     stop statement removed. shall occur only very seldom, &
!     and should be of no harm
      If (r(l)/srnom>=1.0D0) Then
!       NEW PHILOSOPHY FROM V10.1.3 ON
        If (safety==1.1D0) Then
          Print *, ' WARNING!!!  PROBLEMS WITH POINT CLOSE TO PHOTOSPHERE &
            &(MODVL_UPDATE)!'
          Write (999, *) ' WARNING!!!  PROBLEMS WITH POINT CLOSE &
            &TO PHOTOSPHERE (MODVL_UPDATE)!'
          safety = 1.2
          Go To 100
        End If

        Print *, xm(l-2), xm(l-1), xm(l)
        Print *, r(l-2)/srnom, r(l-1)/srnom, r(l)/srnom
        Print *, r(l)/srnom
        Write (999, *) ' STOP: PROBLEMS WITH POINT CLOSE TO PHOTOSPHERE &
          &CANNOT BE CURED (MODVL_UPDATE)!'
        Stop ' PROBLEMS WITH POINT CLOSE TO PHOTOSPHERE CANNOT &
          &BE CURED (MODVL_UPDATE)!'
      End If

!     define new xm grid (equidistant)
      deltam = log10(xmmax/x2)/(nd-l)
      deltam = 10.D0**deltam
      Do ll = l + 1, nd
        xm(ll) = xm(ll-1)*deltam
      End Do
!     JO Jan/April 2016: last point close to previous one, to improve
!     diffusion approx.
      If (nd>=61) xm(nd) = xm(nd-1)*1.1

    End If

    x1 = x2
!   re-calculate R^2*GRAD

    Do i = 1, ndmod - 1
      If (x1>=xmext(i) .And. x1<xmext(i+1)) Exit
    End Do
    If (i==ndmod) Then
      Write (999, *) ' STOP: x1 not found in XM - MODVL_UPDATE'
      Stop ' x1 not found in XM - MODVL_UPDATE'
    End If

    dchidm = (delchi(i+1)-delchi(i))/(xmext(i+1)-xmext(i))
    delchibar = delchi(i) + dchidm*(x1-xmext(i))
    chibar = delchibar*(const*p(l)**expo_p*tnew(l)**expo+sigemin)
    radacc(l) = 5.67D-5/2.9979D10*teff**4*chibar
    If (radacc(l)<=0.) Then
      Write (999, *) ' STOP: RADACC = 0 IN MODVL_UPDATE'
      Stop ' RADACC = 0 IN MODVL_UPDATE'
    End If

    l = l + 1

  End Do
! -------------------------------------------------------------------------

  Do l = ns + 1, nd
    If (p(l)<0.D0) Then
      Write (999, *) ' STOP: PRESSURE NEGATIVE!'
      Stop ' PRESSURE NEGATIVE!'
    End If
    rhonew(l) = p(l)/(a2te*delmu(l)*tnew(l))
!   PRINT*,L,XM(L),XMEXT(L)
!   PRINT*,RNEW(L)/SRNOM,R(L)*SR/SRNOM
!   PRINT*,RHONEW(L),RHO(L)
!   PRINT*,TNEW(L), TEMP(L)
!   PRINT*,TAURNEW(L), TAUR(L)
!   PRINT*,RADACC(L),5.67D-5/2.9979D10*TEFF**4*CHIBAR_H(L)
!   WRITE(*,FMT='(I3,4(2X,E10.4))')
!   L,XM(L),RADACC(L),XMEXT(L),5.67D-5/2.9979D10*TEFF**4*CHIBAR_H(L)
!   PRINT*
  End Do

  srnew = rnew(nd)

  If (abs(1.-sr/srnew)>0.3) Then
    Write (999, *) ' STOP: TOO LARGE CHANGE IN SR -- MODVL_UPDATE'
    Stop ' TOO LARGE CHANGE IN SR -- MODVL_UPDATE'
  End If

  sr = srnew

  rmax = rnew(1)/sr
  rtau23 = srnom/sr
  srvmin = r0
  If (1.D0-1.D0/rtau23>.1D0) Then
    Print *, ' WARNING!!! EXTENDED PHOTOSPHERE'
    Write (999, *) ' WARNING!!! EXTENDED PHOTOSPHERE'
  End If
  Print *, ' LOWER (PHOT.) RADIUS AT ', 1.D0/rtau23, ' NOM. RADII'
  Print *, ' MAX. RADIUS (RMAX) = ', rmax, ' (IN SR)'

  Do l = ns + 1, nd
    v(l) = xmloss/(4.D0*pi*rnew(l)**2*rhonew(l))
  End Do

! modified to ensure smooth transition

! new version
  Do l = ns, nd - 1
    dm = log(v(l)/v(l-1))/(rnew(l)-rnew(l-1))
    dpx = log(v(l)/v(l+1))/(rnew(l)-rnew(l+1))
    gradv(l) = .5D0*(dm+dpx)*v(l)
  End Do

  l = nd
  gradv(l) = log(v(l)/v(l-1))/(rnew(l)-rnew(l-1))*v(l)

! NORMALIZATION
  Do k = 1, nd
    rnew(k) = rnew(k)/sr
    v(k) = v(k)/vmax
    gradv(k) = gradv(k)*sr/vmax
    p(k) = rhonew(k)*(a2te*delmu(k)*tnew(k))
  End Do

  rnew(nd) = 1.D0

! -------------
! JO Apri16: now we 'smooth' the velocity at the transition point,
! requiring a 3rd degree polynom and prescribed v, gradv above and below
! transition point (4 unknowns for the parabola, four conditions to fix them)

  rr1 = rnew(ns-1)
  rr2 = rnew(ns+1)

  zz1 = (/ rr1**3, rr1**2, rr1, 1.D0 /)
  zz2 = (/ rr2**3, rr2**2, rr2, 1.D0 /)
  vmat(1, :) = zz1
  vmat(2, :) = zz2
  vmat(3, :) = (/ 3.D0*rr1**2, 2.D0*rr1, 1.D0, 0.D0 /)
  vmat(4, :) = (/ 3.D0*rr2**2, 2.D0*rr2, 1.D0, 0.D0 /)

  yy = (/ v(ns-1), v(ns+1), gradv(ns-1), gradv(ns+1) /)

  Call ludcmp(vmat, 4, 4, indlud, det)
  Call lubksb(vmat, 4, 4, indlud, yy)

! test precision
  v1test = dot_product(yy, zz1)
  If (abs(1.D0-v1test/v(ns-1))>1.D-4) Then
    Write (999, *) &
      ' STOP: PRECISION OF V1 NOT SUFFICIENT (ROUTINTE MODVL_UPDATE)'
    Stop ' PRECISION OF V1 NOT SUFFICIENT (ROUTINTE MODVL_UPDATE)'
  End If
  v2test = dot_product(yy, zz2)
  If (abs(1.D0-v2test/v(ns+1))>1.D-4) Then
    Write (999, *) &
      ' STOP: PRECISION OF V2 NOT SUFFICIENT (ROUTINTE MODVL_UPDATE)'
    Stop ' PRECISION OF V2 NOT SUFFICIENT (ROUTINTE MODVL_UPDATE)'
  End If

! calculate smoothed velocity at NS
  rr0 = rnew(ns)
  zz1 = (/ rr0**3, rr0**2, rr0, 1.D0 /)
  vsmooth = dot_product(yy, zz1)
  Print *
  Print *, 'OLD/NEW V-VALUES AROUND NS:'
  Write (*, Fmt='(3(F12.6,2X))') v(ns-1), v(ns), v(ns+1)
  Write (*, Fmt='(3(F12.6,2X))') v(ns-1), vsmooth, v(ns+1)

! recalculate rho and p, assuming r^2*rho*v = const
  rhosmooth = rhonew(ns)*v(ns)/vsmooth
  v(ns) = vsmooth
  p(ns) = p(ns)/rhonew(ns)*rhosmooth
  rhonew(ns) = rhosmooth

  gradvnew = gradv
  Do l = ns - 1, ns + 1
    dm = log(v(l)/v(l-1))/(r(l)-r(l-1))
    dpx = log(v(l)/v(l+1))/(r(l)-r(l+1))
    gradvnew(l) = .5D0*(dm+dpx)*v(l)
  End Do

! for tests
  Print *, 'OLD/NEW DVDR-VALUES AROUND NS:'
  Write (*, Fmt='(5F12.6,2X)') gradv(ns-2:ns+2)
  Write (*, Fmt='(5F12.6,2X)') gradvnew(ns-2:ns+2)

  gradv = gradvnew

! ----------------------------------------------------------------------------

  Do l = 1, nd
    If (gradv(l)<=0.D0) Then
      gradv(l) = -gradv(l)
      Print *, &
        ' MODVL_UPDATE: NEGATIVE VEL. GRADIENT -- DENSITY INVERSION! AT N =', &
        l
      Write (999, *) &
        ' MODVL_UPDATE: NEGATIVE VEL. GRADIENT -- DENSITY INVERSION! AT N =', &
        l
    End If
  End Do


! ----from here on, we have to consider clumping
  Call clumping(rnew, v, rhonew, clf, rhonewcl, rtau23, nd, ns, vmax, gradv)

! scaling of ne and nh
  Do l = 1, nd
    dm = rhonewcl(l)/rhocl(l)
    xne(l) = xne(l)*dm
    xnelte(l) = xnelte(l)*dm
    xnh(l) = xnh(l)*dm
  End Do

! assumed that xne for ithyd ne 1 has been corrected

  xnecl = xne/clf !                 corrected for integration of tau_e (see
! subr. MODVL)

  sigem = sigmae*(1.D0+hei*yhe)/summass

  taumin = sigem*xmloss/(4.*pi*vmax*rnew(1)*sr)

  Call spline(rnew, xnecl, arho, brho, crho, nd)
  tau(1) = taumin
  summ = taumin

  Do l = 1, nd - 1
    del = rnew(l+1) - rnew(l)
    summ = summ - sigmae*sr*del*(xnecl(l)+del*(arho(l)/2.D0+del*(brho(l)/3.D0+ &
      del*crho(l)/4.D0)))
    tau(l+1) = summ
  End Do

  Print *
  Print *, ' TAU-THOMSON (MICRO-CLUMPED) AT RSTAR = ', tau(nd)
  Print *

  Do l = 1, nd
    index(l) = l
  End Do


! ---- file 'model' is written: XNE, XNH include CLFAC
! (RHO average density without CLFAC)
! from version 9.0 on: NS as last entry
! from version 10.0 on UPDATED as last entry


! UPDATE done, thus
  updated = .True.
  update_done = .True.
! JO JAN 2016
  no_check = .True.

! note: (hopefully) xnelte (and xnh) from model will not be used during
! restart;
! otherwise, they might need to be replaced by better estimates
  Rewind 20
  Write (20) teff, ggrav, sr, yhe, xmu, vmax, xmloss, beta, &
    (rnew(i), i=1, nd), (v(i), i=1, nd), (gradv(i), i=1, nd), &
    (rhonew(i), i=1, nd), (xnelte(i), i=1, nd), (xnh(i), i=1, nd), &
    (clf(i), i=1, nd), (index(i), i=1, nd), (p(i), i=1, nd), ns, updated, xt

! read photospheric rad. acceleration (start values)
  Open (61, File=trim(modnam)//'/GRAD.OUT', Status='UNKNOWN')
  Do l = 1, nd
    Read (61, *) char, i, arho(l)
  End Do

! store current photospheric rad. acceleration (start values)
  Rewind 61
  Do l = 1, nd
    Write (61, *) 'RADACC ', l, ' ', radacc(l)
  End Do
  Close (61)


! store photospheric rad. acceleration (start)
  Open (61, File=trim(modnam)//'/GRAD_START.OUT', Status='UNKNOWN')
  Do l = 1, nd
    Write (61, *) 'RADACC ', l, ' ', arho(l)
  End Do
  Close (61)

  Do ll = 1, nd
    l = index(nd+1-ll)
    arho(ll) = rnew(l)
    brho(ll) = tnew(l)
  End Do

! ---- file 'temp' is written

  Rewind 21
  Write (21, *) arho, brho
  Write (21, *)(xnelte(i), i=nd, 1, -1)
  Write (21, *) rtau23

  Print *, ' MODEL UPDATED'
  Print *
  Write (999, *) ' MODEL UPDATED'
  Write (999, *)

! REMEMBER: TAU=TAU-TH ONLY IN MICRO-CLUMPING APPROX.
  Call prmodel(rnew, v, gradv, rhonew, xne, xnh, tnew, tau, clf, teff, nd)

! update of variables
  r = rnew
  velo = v
  dvdr = gradv
  rho = rhonew
  temp = tnew

  Return

End Subroutine

!-----------------------------------------------------------------------

Subroutine derivs_upd(x, y, dy)
! basic idea: T(taur) remains preserved, we have only to calculate
! the new taur from the old one and using
! dtau/dm approx se + C rho T^(-x) approx se + C' p T^(-x-1)
! thus, dtaunew/dm approx [(dtauold/dm -se)*pnew/pold (tnew/told)^(-x-1)]+se

! with se = ne * sigmae/rho := sig

! Note that this assumes that the parameterization of Chi_Ross did not change,
! i.e., C' and x (from first hydro model) remain unchanged.

! In the standard approx. (OPTCMF_FULL = .FALSE.), that's no problem.
! However, for OPTCMF_FULL = .TRUE., it is, since now taur includes a
! multitude of CMF lines, and T(taur) has already been updated, so that the
! parameterization (w.r.t. the first hydro-solution) is certainly different.

! However, (few) tests have shown that also this is not a real problem, since
! both taur and T(taur) change in a reasonable way. (Basically, the change in
! taur (when actual values are used) is not so large, so that also T(m(taur))
! does not change significantly.)

! What we really want to know is how the density changes due to changes
! in grad alone, and this requires that we know (in addition to grad)
! P and T as a function of m. As long as T(m) is OK, everything is fine,
! and the latter is warrented because we assume that T(taur) is preserved,
! and only taur needs to be updated (with respect to the actual (CMF) taur.)

! Y(1) --- pressure
! Y(2) --- temperature
! Y(3) --- radius
! Y(4) --- taur_new

  Use :: nlte_type
  Use :: nlte_dim
  Use :: nlte_opt, Only: opt_photclump
  Use :: nlte_var, Only: srnom, idum, ggrav, teff, a2te, delmu
  Use :: photstruc, Only: ndmod, xt, const, expo, expo_p, sigemin, xmext, &
    text, delchi, pold => pkur, taurold => xnekur, sig => gradkur, clf_const

  Implicit None

! delmu  corrects for constant mue

! .. scalar arguments ..
  Real (dp) :: x
! ..
! .. array arguments ..
  Real (dp) :: dy(4), y(4)
! ..
! .. local scalars ..
  Real (dp) :: delchibar, dchidm, chibar, delm, dy1, dpdm, po, dtdm, to, &
    dtaudm, dum
! ..
  Integer (i4b) :: it

  If (ndmod/=id_ndept) Then
    Write (999, *) ' STOP: ERROR IN NDMOD - DERIVS_UPD'
    Stop ' ERROR IN NDMOD - DERIVS_UPD'
  End If

! CHIBAR_H, POLD, TOLD  on regular grid
  Do it = 1, ndmod - 1
    If (x>=xmext(it) .And. x<xmext(it+1)) Exit
  End Do
  If (it==ndmod) Then
    Write (999, *) ' STOP: x not found in XM - DERIVS_UPD'
    Stop ' x not found in XM - DERIVS_UPD'
  End If

  dchidm = (delchi(it+1)-delchi(it))/(xmext(it+1)-xmext(it))
  delchibar = delchi(it) + dchidm*(x-xmext(it))

  dpdm = (pold(it+1)-pold(it))/(xmext(it+1)-xmext(it))
  po = pold(it) + dpdm*(x-xmext(it))

  dtdm = (text(it+1)-text(it))/(xmext(it+1)-xmext(it))
  to = text(it) + dtdm*(x-xmext(it))

  dtaudm = (taurold(it+1)-taurold(it))/(xmext(it+1)-xmext(it))

! NEW (former version: < 1.*SIG(IT)
  If (dtaudm<=.95*sig(it)) Then
    Print *, x, it, dtaudm, sig(it)
    Print *
    Do it = 1, ndmod
      Print *, it, xmext(it), delchi(it), pold(it), text(it), taurold(it), &
        sig(it)
    End Do
    Write (999, *) ' STOP: DTAUDM < SIG_TH'
    Stop ' DTAUDM < SIG_TH'
  End If

! FROM REGRESSION
! JO Feb. 2023: new exponent EXPO_P (=1 if linear regression successful)
  chibar = delchibar*(const*y(1)**expo_p*y(2)**expo+sigemin)

  dy1 = ggrav - 5.67D-5/2.9979D10*teff**4*chibar
  dum = ggrav/(5.67D-5/2.9979D10*teff**4)

  If (dy1<0.) Then
    Print *, ' DELTA P < 0 IN DERIVS_UPD (LOCALLY BEYOND EDDINGTON LIMIT)'
    Print *, ' ALLOWED MAXIMUM VALUE = ', dum, ' ACTUAL VALUE = ', dy1
    Write (999, *) &
      ' DELTA P < 0 IN DERIVS_UPD (LOCALLY BEYOND EDDINGTON LIMIT)'
    Write (999, *) ' ALLOWED MAXIMUM VALUE = ', dum, ' ACTUAL VALUE = ', dy1
    Write (999, *) ' STOP: DELTA P < 0 IN DERIVS_UPD'
    Stop ' DELTA P < 0 IN DERIVS_UPD'
  End If

! this is the ONLY hack to ensure that rho_phot => rho_phot/fcl
  If (opt_photclump) dy1 = dy1/clf_const
  dy(1) = dy1*(srnom/y(3))**2

! DTAUDM FROM OLD TAUR, POLD AND TOLD AND NEW P AND T, + EXPONENT
  dy(4) = ((dtaudm-sig(it))*(y(1)/po)*(to/y(2))**(xt+1.)) + sig(it)

! FIND INTERVAL WITH RESPECT TO NEW TAUR
  Do it = 1, ndmod - 1
    If (y(4)>=taurold(it) .And. y(4)<taurold(it+1)) Exit
  End Do
  If (it==ndmod) it = ndmod - 1

! DTDM = DTDTAU|new * DTAUDM
  dy(2) = (text(it+1)-text(it))/(taurold(it+1)-taurold(it))*dy(4)

  delm = .5D0*(delmu(idum-1)+delmu(idum))

  dy(3) = -y(2)/y(1)*a2te*delm

  Return

End Subroutine

!-----------------------------------------------------------------------

Subroutine chih_params(chi, p, t, delchi, sigemin, const, expo, expo_p, nd, &
  ns, nsig)

! parameterizes CHIBAR_H

  Use :: nlte_type
  Use :: nlte_dim, Only: id_ndept
  Implicit None

  Integer (i4b), Parameter :: nd1 = id_ndept

  Integer (i4b), Intent (In) :: nd, ns, nsig
  Real (dp), Dimension (nd), Intent (In) :: chi, p, t
  Real (dp), Dimension (nd), Intent (Out) :: delchi

  Real (dp), Intent (In) :: sigemin
  Real (dp), Intent (Out) :: const, expo, expo_p ! JO Feb. 2023 new
! fit-parameter

  Real (dp), Dimension (nd1) :: fit, sig, fit_lin, delchi_lin

  Integer (i4b) :: imin, i, ntot, minimin, minimin_lin
  Real (dp) :: chi2min, expoi, consti, corr, chi2, chisq
  Real (dp) :: expo_lin, expo_p_lin, const_lin


  Integer (i4b), Dimension (3) :: ia
  Real (dp), Dimension (3, 3) :: covar
  Real (dp), Dimension (3) :: a
  External :: lfit_funcs

  chi2min = 1.D5
  expo = 0.
  minimin = 0


! JO Feb. 2023: perform at first linear test (old approach), and
! accept values if fit OK. If not, perform multi-linear fit
! in the former case, new exponent for P, EXPO_P, is set to unity

! test for best interval to perform the fit
  Do imin = nsig, nd - 4
!   regression from IMIN on
    Call linreg(log10(t(imin:nd)), log10((chi(imin:nd)- &
      sigemin)/p(imin:nd)), nd+1-imin, expoi, consti, corr)

!   fit quality for "all" points (from NSIG on)
    fit(nsig:nd) = (10.**consti*t(nsig:nd)**expoi*p(nsig:nd)) + sigemin

!   old formulation (including error for consistency, see below)
    chi2 = sum((chi(nsig:nd)-fit(nsig:nd))**2/chi(nsig:nd))
!   JO Jan 2016 -- Note:
!   actually, the denominator should be CHI^2 because we want to minimize
!   (fit/chi-1)^2 = (fit-chi)^2/chi^2 = (chi-fit)^2/chi^2,
!   but for consistency with older versions we leave this here.
!   (This routine is only called for OPTCMF_FULL = F, and in this case
!   the parameterization makes no problem anyway, independent of the choice of
!   chi2
    If (chi2<chi2min) Then
      minimin = imin
      chi2min = chi2
      expo = expoi
      expo_p = 1.D0 !               J0: new parameter
      const = consti
    End If
  End Do

  If (expo==0.) Then
    Write (999, *) ' STOP: LINEAR REGRESSION NOT SUCCESSFUL (CHIH_PARAMS)'
    Stop ' LINEAR REGRESSION NOT SUCCESSFUL (CHIH_PARAMS)'
  End If

  const = 10.**const
! for all photospheric points (from NS on)
  fit(ns:nd) = (const*t(ns:nd)**expo*p(ns:nd)) + sigemin
  delchi(ns:nd) = chi(ns:nd)/fit(ns:nd)

! new condition
  If (maxval(delchi(ns:nd))<=1.2 .And. minval(delchi(ns:nd))>=0.8) Then
!   for tests, if multi-linear fit shall be enforced
!   IF(MAXVAL(DELCHI(NS:ND)).LE.1.02 .AND. MINVAL(DELCHI(NS:ND)).GE.0.98) THEN
    Print *, ' LINEAR FIT IN CHIH_PARAMS SUCCESSFUL'
    Print *, ' NS = ', ns, ' NSIG = ', nsig, ' IMIN = ', minimin
    Print *, ' PARAMETERS: ', sigemin, ' ', const, ' ', expo, ' ', expo_p
    Print *, ' DEVIATIONS IN RANGE ', minval(delchi(ns:nd)), ' ', &
      maxval(delchi(ns:nd))
    Print *
    Write (999, *) ' LINEAR FIT IN CHIH_PARAMS SUCCESSFUL'
    Write (999, *) ' NS = ', ns, ' NSIG = ', nsig, ' IMIN = ', minimin
    Write (999, *) ' PARAMETERS: ', sigemin, ' ', const, ' ', expo, ' ', &
      expo_p
    Write (999, *) ' DEVIATIONS IN RANGE ', minval(delchi(ns:nd)), ' ', &
      maxval(delchi(ns:nd))
    Write (999, *)

    Do i = ns, nd
      Write (*, 110) i, chi(i), fit(i), t(i), p(i)
    End Do
    Print *
    Return
  Else
!   DO I=NS,ND
!   WRITE(*,10) I,CHI(I),FIT(I),DELCHI(I)
!   ENDDO
    Print *, ' LINEAR FIT IN CHIH_PARAMS NOT GOOD ENOUGH!'
    Print *, ' PARAMETERS: ', sigemin, ' ', const, ' ', expo, ' ', expo_p
    Print *, ' TRYING MULTI-LINEAR REGRESSION NOW'
    Print *
    Write (999, *) ' LINEAR FIT IN CHIH_PARAMS NOT GOOD ENOUGH!'
    Write (999, *) ' PARAMETERS: ', sigemin, ' ', const, ' ', expo, ' ', &
      expo_p
    Write (999, *) ' TRYING MULTI-LINEAR REGRESSION NOW'
    Write (999, *)
!   save values from linear regression, might be used below in case
    minimin_lin = minimin
    const_lin = const
    expo_lin = expo
    expo_p_lin = expo_p
    fit_lin = fit
    delchi_lin = delchi
  End If

! multi-linear regression including P (see subr. LFIT_FUNCS)

! test for best interval to perform the fit
  sig = 1. !                        no error
  ia = 1 !                          all parameters


  chi2min = 1.D5
  expo = 0.
  minimin = 0

  Do imin = nsig, nd - 4
!   regression from IMIN on

    ntot = nd - imin + 1
    Call lfit(log10(t(imin:nd)), log10(chi(imin:nd)-sigemin), sig(imin:nd), &
      p(imin:nd), ntot, a, ia, 3, covar, 3, chisq, lfit_funcs)
!   fit parameters contained in A

!   fit quality for "all" points (from NSIG on)
    fit(nsig:nd) = (10.**a(1)*t(nsig:nd)**a(2)*p(nsig:nd)**a(3)) + sigemin
!   in case, regarding log-values
!   CHI2=SUM((LOG10(CHI(NSIG:ND))-LOG10(FIT(NSIG:ND)))**2/LOG10(CHI(NSIG:ND))**2)

!   regarding linear deviation, correct formulation (denominator!)
    chi2 = sum((chi(nsig:nd)-fit(nsig:nd))**2/chi(nsig:nd)**2)
!   PRINT*,' PARAMETERS: ',IMIN,' ',SIGEMIN, CHI2
!   PRINT*,' PARAMETERS: ',A(1),' ',A(2),' ',A(3)
!   PRINT*

    If (chi2<chi2min) Then
      minimin = imin
      chi2min = chi2
      expo = a(2)
      expo_p = a(3)
      const = a(1)
    End If
  End Do

  If (expo==0.) Then
    Write (999, *) &
      ' STOP: MULTI-LINEAR REGRESSION NOT SUCCESSFUL (CHIH_PARAMS)'
    Stop ' MULTI-LINEAR REGRESSION NOT SUCCESSFUL (CHIH_PARAMS)'
  End If

  const = 10.**const
! for all photospheric points (from NS on)
  fit(ns:nd) = (const*t(ns:nd)**expo*p(ns:nd)**expo_p) + sigemin
  delchi(ns:nd) = chi(ns:nd)/fit(ns:nd)

  If (maxval(delchi(ns:nd))>2. .Or. minval(delchi(ns:nd))<0.2) Then
!   still old condition
    Do i = ns, nd
      Write (*, 100) i, chi(i), fit(i), delchi(i)
    End Do
    Print *, ' NO MULTI-LINEAR REGRESSION POSSIBLE (CHIH_PARAMS)!'
    Write (999, *) ' NO MULTI-LINEAR REGRESSION POSSIBLE (CHIH_PARAMS)!'
  Else If (expo>0. .Or. expo_p<0.) Then
!   unphysical parameters, may appear if locally close to Eddington limit'
    Print *, ' MULTI-LINEAR REGRESSION PARAMETERS NOT APPROPRIATE!'
    Print *, ' EXPO FOR T (SHOULD BE NEGATIVE) = ', expo
    Print *, ' EXPO FOR P (SHOULD BE IN BETWEEN 0...1) = ', expo_p
    Write (999, *) ' MULTI-LINEAR REGRESSION PARAMETERS NOT APPROPRIATE!'
    Write (999, *) ' EXPO FOR T (SHOULD BE NEGATIVE) = ', expo
    Write (999, *) ' EXPO FOR P (SHOULD BE IN BETWEEN 0...1) = ', expo_p
  Else
    Print *, ' MULTI-LINEAR FIT IN CHIH_PARAMS SUCCESSFUL'
    Write (999, *) ' MULTI-LINEAR FIT IN CHIH_PARAMS SUCCESSFUL'
    If (expo_p>1.) Then
      Print *, ' EXPO FOR P = ', expo_p, ' > 1, PROCEED AT OWN RISK!!!'
      Write (999, *) ' EXPO FOR P = ', expo_p, ' > 1, PROCEED AT OWN RISK!!!'
    End If
    Print *, ' NS = ', ns, ' NSIG = ', nsig, ' IMIN = ', minimin
    Print *, ' PARAMETERS: ', sigemin, ' ', const, ' ', expo, ' ', expo_p
    Print *, ' DEVIATIONS IN RANGE ', minval(delchi(ns:nd)), ' ', &
      maxval(delchi(ns:nd))
    Print *
    Write (999, *) ' NS = ', ns, ' NSIG = ', nsig, ' IMIN = ', minimin
    Write (999, *) ' PARAMETERS: ', sigemin, ' ', const, ' ', expo, ' ', &
      expo_p
    Write (999, *) ' DEVIATIONS IN RANGE ', minval(delchi(ns:nd)), ' ', &
      maxval(delchi(ns:nd))
    Write (999, *)
    Do i = ns, nd
      Write (*, 110) i, chi(i), fit(i), t(i), p(i)
    End Do
    Print *
    Return
  End If

! final possibility: if linear fit not too bad, use these values
  If (maxval(delchi_lin(ns:nd))<=2 .And. minval(delchi_lin(ns:nd))>=0.2) Then
    Print *, ' LINEAR REGRESSION USED INSTEAD'
    Write (999, *) ' LINEAR REGRESSION USED INSTEAD'
    minimin = minimin_lin
    const = const_lin
    expo = expo_lin
    expo_p = expo_p_lin
    fit = fit_lin
    delchi = delchi_lin
    Print *, ' NS = ', ns, ' NSIG = ', nsig, ' IMIN = ', minimin
    Print *, ' PARAMETERS: ', sigemin, ' ', const, ' ', expo, ' ', expo_p
    Print *, ' DEVIATIONS IN RANGE ', minval(delchi(ns:nd)), ' ', &
      maxval(delchi(ns:nd))
    Write (999, *) ' NS = ', ns, ' NSIG = ', nsig, ' IMIN = ', minimin
    Write (999, *) ' PARAMETERS: ', sigemin, ' ', const, ' ', expo, ' ', &
      expo_p
    Write (999, *) ' DEVIATIONS IN RANGE ', minval(delchi(ns:nd)), ' ', &
      maxval(delchi(ns:nd))
    Do i = ns, nd
      Write (*, 110) i, chi(i), fit(i), t(i), p(i)
    End Do
    Print *
    Return
  Else
!   problems cannot be cured
    Write (999, *) ' STOP: NO REGRESSION POSSIBLE (CHIH_PARAMS)! &
      &RE-TRY WITH UPDATE_STRUCT = .FALSE.'
    Stop ' NO REGRESSION POSSIBLE (CHIH_PARAMS)! RE-TRY WITH UPDATE_STRUCT &
      &= .FALSE.'
  End If

  Return

100 Format (I3, 2X, 3(E10.3,2X))
110 Format (I3, 2X, 4(E10.3,2X))

End Subroutine

!-----------------------------------------------------------------------

!next three subroutines for TRATCOL, NFOR=66
Subroutine spline_col(x, y, n, yp1, ypn, y2)

! Cubic spline interpolation (Numerical Recipes p.88)
! Given arrays X and Y of legth N containing a tabulated function
! y(i)=f(x(i))
! and given values YP1 and YPN for the first derivative of the
! interpolating function at points 1 and N, returns an array Y2
! of length N which contains second derivatives of the interpolating
! function at tabulated points x(i)

! uses natural splines only if derivatives at boundary are extremely large

  Use :: nlte_type
  Implicit None

  Integer (i4b), Parameter :: nmax = 100
  Integer (i4b) :: i, n, k

  Real (dp), Dimension (n) :: x, y, y2
  Real (dp), Dimension (nmax) :: u

  Real (dp) :: yp1, ypn, sig, p, qn, un


  If (yp1>.99D30) Then
    y2(1) = 0.
    u(1) = 0.
  Else
    y2(1) = -0.5
    u(1) = (3./(x(2)-x(1)))*((y(2)-y(1))/(x(2)-x(1))-yp1)
  End If

  Do i = 2, n - 1

    sig = (x(i)-x(i-1))/(x(i+1)-x(i-1))
    p = sig*y2(i-1) + 2
    y2(i) = (sig-1)/p
    u(i) = (6.*((y(i+1)-y(i))/(x(i+1)-x(i))-(y(i)-y(i-1))/(x(i)-x(i-1)))/(x(i+ &
      1)-x(i-1))-sig*u(i-1))/p

  End Do

  If (ypn==.99D30) Then
    qn = 0
    un = 0
  Else
    qn = 0.5
    un = (3./(x(n)-x(n-1)))*(ypn-(y(n)-y(n-1))/(x(n)-x(n-1)))
  End If

  y2(n) = (un-qn*u(n-1))/(qn*y2(n-1)+1.)

  Do k = n - 1, 1, -1
    y2(k) = y2(k)*y2(k+1) + u(k)
  End Do

  Return

End Subroutine

!-----------------------------------------------------------------------

Subroutine splinter(xa, ya, y2a, n, x, y)

! Cubic spline interpolation (Numerical Recipes p.89)
! Given arrays xa aand ya of length n, which tabulated function and given
! array y2a (output from spline_col) and given a value of x, returns a
! cubic spline interpolated value y

  Use :: nlte_type

  Implicit None

  Integer (i4b) :: k, klo, khi, n
  Real (dp), Dimension (n) :: xa, ya, y2a

  Real (dp) :: h, a, b, y, x

  klo = 1
  khi = n

100 If (khi-klo>1) Then
    k = (khi+klo)/2

    If (xa(k)>x) Then
      khi = k
    Else
      klo = k
    End If
    Go To 100
  End If

  h = xa(khi) - xa(klo)

  If (h==0.D0) Then
    Write (999, *) ' STOP: Bad xa input in routine SPLINTER'
    Stop ' Bad xa input in routine SPLINTER'
  End If

  a = (xa(khi)-x)/h
  b = (x-xa(klo))/h
  y = a*ya(klo) + b*ya(khi) + ((a**3-a)*y2a(klo)+(b**3-b)*y2a(khi))*(h**2)/6.

  Return
End Subroutine

!-----------------------------------------------------------------------

Subroutine devp(x, y, n, yp1, ypn)

! Given an array of length N, returns first derivative at points 1 and N

  Use :: nlte_type

  Implicit None

  Integer (i4b) :: n
  Real (dp), Dimension (n) :: x, y
  Real (dp) :: h1, h2, yp1, ypn

  h1 = x(2) - x(1)
  h2 = x(n) - x(n-1)

  yp1 = ((y(2))-y(1))/h1
  ypn = ((y(n)-y(n-1)))/h2

  Return
End Subroutine

!-----------------------------------------------------------------------

Subroutine deriv_3p(d, x, y, n)

! Perform numerical differentiation using 3-point, Lagrangian interpolation.
! df/dx =
! y0*(2x-x1-x2)/(x01*x02)+y1*(2x-x0-x2)/(x10*x12)+y2*(2x-x0-x1)/(x20*x21)
! Where: x01 = x0-x1, x02 = x0-x2, x12 = x1-x2, etc.

! taken from IDL; note that cshift (f90) and shift(idl) use different sign
! conventions

  Use :: nlte_type

  Implicit None

  Integer (i4b) :: n, n2
  Real (dp), Dimension (n) :: x, y, d

  Real (dp), Dimension (n) :: x12, x01, x02

  If (n<3) Then
    Write (999, *) ' STOP: N LT 3 IN DERIV_3P'
    Stop ' N LT 3 IN DERIV_3P'
  End If

  x12 = x - cshift(x, 1) !          x1 - x2
  x01 = cshift(x, -1) - x !         x0 - x1
  x02 = cshift(x, -1) - cshift(x, 1) ! x0 - x2

! Middle points
  d = cshift(y, -1)*(x12/(x01*x02)) + y*(1./x12-1./x01) - &
    cshift(y, 1)*(x01/(x02*x12))

! Formulae for the first and last points:
  d(1) = y(1)*(x01(2)+x02(2))/(x01(2)*x02(2)) - y(2)*x02(2)/(x01(2)*x12(2)) + &
    y(3)*x01(2)/(x02(2)*x12(2))

  n2 = n - 1
  d(n) = -y(n-2)*x12(n2)/(x01(n2)*x02(n2)) + y(n-1)*x02(n2)/(x01(n2)*x12(n2)) &
    - y(n)*(x02(n2)+x12(n2))/(x02(n2)*x12(n2))

  Return
End Subroutine

Subroutine delete_old_model_files

! the following files might be left from previous models with OPTCMF_FULL = T.
! To avoid confusion, and since CONT_FORMAL_CMF must not be present for
! formal if OPTCMF_FULL = F, all these files are deleted.
  Use :: nlte_var, Only: modnam

  Implicit None

  Logical :: exi

  Inquire (File=trim(modnam)//'/LINE_LIST_MERGED', Exist=exi)
  If (exi) Then
    Open (1, File=trim(modnam)//'/LINE_LIST_MERGED', Status='OLD')
    Close (1, Status='DELETE')
    Print *, 'old file LINE_LIST_MERGED deleted'
    Write (999, *) 'old file LINE_LIST_MERGED deleted'
  End If

  Inquire (File=trim(modnam)//'/CONT_FORMAL_CMF', Exist=exi)
  If (exi) Then
    Open (1, File=trim(modnam)//'/CONT_FORMAL_CMF', Status='OLD')
    Close (1, Status='DELETE')
    Print *, 'old file CONT_FORMAL_CMF deleted'
    Write (999, *) 'old file CONT_FORMAL_CMF deleted'
  End If

  Inquire (File=trim(modnam)//'/KCMF', Exist=exi)
  If (exi) Then
    Open (1, File=trim(modnam)//'/KCMF', Status='OLD')
    Close (1, Status='DELETE')
    Print *, 'old file KCMF deleted'
    Write (999, *) 'old file KCMF deleted'
  End If

  Inquire (File=trim(modnam)//'/CHIBAR_H_CMF', Exist=exi)
  If (exi) Then
    Open (1, File=trim(modnam)//'/CHIBAR_H_CMF', Status='OLD')
    Close (1, Status='DELETE')
    Print *, 'old file CHIBAR_H_CMF deleted'
    Write (999, *) 'old file CHIBAR_H_CMF deleted'
  End If

  Inquire (File=trim(modnam)//'/OCCNG', Exist=exi)
  If (exi) Then
    Open (1, File=trim(modnam)//'/OCCNG', Status='OLD')
    Close (1, Status='DELETE')
    Print *, 'old file OCCNG deleted'
    Write (999, *) 'old file OCCNG deleted'
  End If

  Inquire (File=trim(modnam)//'/OCC_NO', Exist=exi)
  If (exi) Then
    Open (1, File=trim(modnam)//'/OCC_NO', Status='OLD')
    Close (1, Status='DELETE')
    Print *, 'old file OCC_NO deleted'
    Write (999, *) 'old file OCC_NO deleted'
  End If

  Inquire (File=trim(modnam)//'/out_xj_cmf', Exist=exi)
  If (exi) Then
    Open (1, File=trim(modnam)//'/out_xj_cmf', Status='OLD')
    Close (1, Status='DELETE')
    Print *, 'old file out_xj_cmf deleted'
    Write (999, *) 'old file out_xj_cmf deleted'
  End If

  Inquire (File=trim(modnam)//'/out_xj_app', Exist=exi)
  If (exi) Then
    Open (1, File=trim(modnam)//'/out_xj_app', Status='OLD')
    Close (1, Status='DELETE')
    Print *, 'old file out_xj_app deleted'
    Write (999, *) 'old file out_xj_app deleted'
  End If

  Inquire (File=trim(modnam)//'/out_xj_cmf_coarse', Exist=exi)
  If (exi) Then
    Open (1, File=trim(modnam)//'/out_xj_cmf_coarse', Status='OLD')
    Close (1, Status='DELETE')
    Print *, 'old file out_xj_cmf_coarse deleted'
    Write (999, *) 'old file out_xj_cmf_coarse deleted'
  End If

  Inquire (File=trim(modnam)//'/out_xh_obs', Exist=exi)
  If (exi) Then
    Open (1, File=trim(modnam)//'/out_xh_obs', Status='OLD')
    Close (1, Status='DELETE')
    Print *, 'old file out_xh_obs deleted'
    Write (999, *) 'old file out_xh_obs deleted'
  End If

  Inquire (File=trim(modnam)//'/out_xh_obs_coarse', Exist=exi)
  If (exi) Then
    Open (1, File=trim(modnam)//'/out_xh_obs_coarse', Status='OLD')
    Close (1, Status='DELETE')
    Print *, 'old file out_xh_obs_coarse deleted'
    Write (999, *) 'old file out_xh_obs_coarse deleted'
  End If

  Inquire (File=trim(modnam)//'/out_rad_force', Exist=exi)
  If (exi) Then
    Open (1, File=trim(modnam)//'/out_rad_force', Status='OLD')
    Close (1, Status='DELETE')
    Print *, 'old file out_rad_force deleted'
    Write (999, *) 'old file out_rad_force deleted'
  End If

  Print *
  Write (999, *)

  Return
End Subroutine

!-----------------------------------------------------------------------

!JO Feb. 2023
!program package to fit chih according to
!chih_fit = (const * t**exp1 * p**exp2)+ chimin
!parameters contained in a=[const,exp1,exp2]

Subroutine lfit_funcs(x, p, afunc, ma)
! adapted for fit (see above), including P
  Use :: nlte_type
  Implicit None

  Integer (i4b) :: ma
  Real (dp) :: x, p

  Real (dp), Dimension (ma) :: afunc

  afunc(1) = 1.
  afunc(2) = x
  afunc(3) = log10(p)

  Return
End Subroutine

!-----------------------------------------------------------------------

Subroutine lfit(x, y, sig, p, ndat, a, ia, ma, covar, npc, chisq, funcs)
! multi-linear fit, from NumRec

  Use :: nlte_type
  Implicit None

  Integer (i4b), Parameter :: mmax = 3
  Integer (i4b) :: ma, npc, ndat
  Integer (i4b), Dimension (ma) :: ia

  Real (dp) :: chisq

  Real (dp), Dimension (ma) :: a
  Real (dp), Dimension (npc, npc) :: covar
  Real (dp), Dimension (ndat) :: x, y, sig, p

  Integer (i4b) :: i, j, k, l, m, mfit
  Real (dp) :: sig2i, summ, wt, ym

  Real (dp), Dimension (mmax) :: afunc, beta

  External :: funcs
! U    USES COVSRT,GAUSSJ


  mfit = 0
  Do j = 1, ma
    If (ia(j)/=0) mfit = mfit + 1
  End Do
  If (mfit==0) Then
    Write (999, *) ' STOP: LFIT: NO PARAMETERS TO BE FITTED'
    Stop ' LFIT: NO PARAMETERS TO BE FITTED'
  End If
  Do j = 1, mfit
    Do k = 1, mfit
      covar(j, k) = 0.
    End Do
    beta(j) = 0.
  End Do
  Do i = 1, ndat
    Call funcs(x(i), p(i), afunc, ma)
    ym = y(i)
    If (mfit<ma) Then
      Do j = 1, ma
        If (ia(j)==0) ym = ym - a(j)*afunc(j)
      End Do
    End If
    sig2i = 1./sig(i)**2
    j = 0
    Do l = 1, ma
      If (ia(l)/=0) Then
        j = j + 1
        wt = afunc(l)*sig2i
        k = 0
        Do m = 1, l
          If (ia(m)/=0) Then
            k = k + 1
            covar(j, k) = covar(j, k) + wt*afunc(m)
          End If
        End Do
        beta(j) = beta(j) + ym*wt
      End If
    End Do
  End Do
  Do j = 2, mfit
    Do k = 1, j - 1
      covar(k, j) = covar(j, k)
    End Do
  End Do
  Call gaussj(covar, mfit, npc, beta, 1, 1)
  j = 0
  Do l = 1, ma
    If (ia(l)/=0) Then
      j = j + 1
      a(l) = beta(j)
    End If
  End Do
  chisq = 0.
  Do i = 1, ndat
    Call funcs(x(i), p(i), afunc, ma)
    summ = 0.
    Do j = 1, ma
      summ = summ + a(j)*afunc(j)
    End Do
    chisq = chisq + ((y(i)-summ)/sig(i))**2
  End Do

! uncomment if covariance matrix required
! CALL COVSRT(COVAR,NPC,MA,IA,MFIT)
  Return
End Subroutine

!-----------------------------------------------------------------------

Subroutine gaussj(a, n, np, b, m, mp)
! used for lfit, from NumRec

  Use :: nlte_type
  Implicit None

  Integer (i4b), Parameter :: nmax = 50
  Integer (i4b) :: m, mp, n, np

  Real (dp), Dimension (np, np) :: a
  Real (dp), Dimension (np, mp) :: b

  Integer (i4b) :: i, icol, irow, j, k, l, ll
  Integer (i4b), Dimension (nmax) :: indxc, indxr, ipiv

  Real (dp) :: big, dum, pivinv

  If (np>nmax) Then
    Write (999, *) ' STOP: NP > NMAX in GAUSSJ'
    Stop ' NP > NMAX in GAUSSJ'
  End If

  Do j = 1, n
    ipiv(j) = 0
  End Do
  Do i = 1, n
    big = 0.
    Do j = 1, n
      If (ipiv(j)/=1) Then
        Do k = 1, n
          If (ipiv(k)==0) Then
            If (abs(a(j,k))>=big) Then
              big = abs(a(j,k))
              irow = j
              icol = k
            End If
          Else If (ipiv(k)>1) Then
            Write (999, *) ' STOP: SINGULAR MATRIX IN GAUSSJ 1'
            Stop ' SINGULAR MATRIX IN GAUSSJ 1'
          End If
        End Do
      End If
    End Do
    ipiv(icol) = ipiv(icol) + 1
    If (irow/=icol) Then
      Do l = 1, n
        dum = a(irow, l)
        a(irow, l) = a(icol, l)
        a(icol, l) = dum
      End Do
      Do l = 1, m
        dum = b(irow, l)
        b(irow, l) = b(icol, l)
        b(icol, l) = dum
      End Do
    End If
    indxr(i) = irow
    indxc(i) = icol
    If (a(icol,icol)==0.) Then
      Write (999, *) ' STOP: SINGULAR MATRIX IN GAUSSJ 2'
      Stop ' SINGULAR MATRIX IN GAUSSJ 2'
    End If
    pivinv = 1./a(icol, icol)
    a(icol, icol) = 1.
    Do l = 1, n
      a(icol, l) = a(icol, l)*pivinv
    End Do
    Do l = 1, m
      b(icol, l) = b(icol, l)*pivinv
    End Do
    Do ll = 1, n
      If (ll/=icol) Then
        dum = a(ll, icol)
        a(ll, icol) = 0.
        Do l = 1, n
          a(ll, l) = a(ll, l) - a(icol, l)*dum
        End Do
        Do l = 1, m
          b(ll, l) = b(ll, l) - b(icol, l)*dum
        End Do
      End If
    End Do
  End Do
  Do l = n, 1, -1
    If (indxr(l)/=indxc(l)) Then
      Do k = 1, n
        dum = a(k, indxr(l))
        a(k, indxr(l)) = a(k, indxc(l))
        a(k, indxc(l)) = dum
      End Do
    End If
  End Do

  Return
End Subroutine


!-----------------------------------------------------------------------

Subroutine covsrt(covar, npc, ma, ia, mfit)
! in case, used for lfit, from NumRec

  Use :: nlte_type
  Implicit None

  Integer (i4b) :: ma, mfit, npc
  Integer (i4b), Dimension (ma) :: ia

  Real (dp), Dimension (npc, npc) :: covar

  Integer (i4b) :: i, j, k
  Real (dp) :: swap

  Do i = mfit + 1, ma
    Do j = 1, i
      covar(i, j) = 0.
      covar(j, i) = 0.
    End Do
  End Do
  k = mfit
  Do j = ma, 1, -1
    If (ia(j)/=0) Then
      Do i = 1, ma
        swap = covar(i, k)
        covar(i, k) = covar(i, j)
        covar(i, j) = swap
      End Do
      Do i = 1, ma
        swap = covar(k, i)
        covar(k, i) = covar(j, i)
        covar(j, i) = swap
      End Do
      k = k - 1
    End If
  End Do

  Return
End Subroutine
