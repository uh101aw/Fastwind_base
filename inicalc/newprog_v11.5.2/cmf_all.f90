module version_cmf_all
!
! NOTES for further modifications
! NOTE: so far, everything with constant vturb!
! minimum vturb = vturbmin (5 km/s); this value is used even if input lower!
!                 modify vturbmin if you really need lower values (e.g., for grad)  
!
! IMPORTANT NOTE (referring to sub. line_list)
!  lines from background elements de-packed (all multiplet components!)
!  example: CIV 1548, 1550 with same index (id,id1)
!  vs.
!  lines from explicit elements (mostly) packed (from Detail input)
!  specific lines are depacked, according to 'linetab' (LINES_xxx.DAT) from ATOM_FILE
!  Might need to be de-packed in final version (idea: use info from master line-list)
!  example: NIII triplet, only one packed transition at 4639
!
! IMPORTANT NOTE (referring to sube. PROFILE_CMF)
!  so far, only Doppler profiles. At least for H/He, Stark broadening needs
!  to be included to obtain realistic fluxes at Lyman series limit
!
! NOTE regarding opacities: for selected (bg-)elements, opacities and
!  source functions(OPACITY_CMF) only calculated in between
!  met_imin(l), met_imax(l) [controlled by occngold ne 0]
!  other opacities set to zero.
!  Similar applies for line-rate coefficients incl. ALO (CALC_TLU_MET)
!  Only for Ng-extrapolation, all defined (resonance) transitions  
!  are accounted for, by using occng and occngold2
!  [however, since rateeq and opacities/source-functions use only
!  transitions defined from occngold ne 0, additional info might be not used]
!  For non-selected (bg-elements), opacities/source-functions calculated
!  from occng, in the spirit of nlte_approx.
!
!  THUS, CAREFUL INSPECTION OF MET_IMIN, MET_IMAX in certain cases is required
!  E.g., the overlap between the OIII/NIII resonance lines might be hampered
!  in hotter models, since: outside OIII is present, in the middle wind not
!  (OIV-OVI), in the sonic region OIII is present again, and then
!  vanishes in the deeper photosphere.
!
! WRITTEN: 
!        04/12 by JP
!
! HISTORY:
!         comletely new
!
! version 0.0 April 2012
! version 0.1 Oct   2013
! version 1.0 June  2015
! version 1.1 April 2016
! version 1.2 April 2016
! version 1.3 Oct   2016
! version 1.3.1 end Oct 2016: Bug in u0 found and fixed (should play almost no role)
! version 1.3.2 Dec 2016: Sobo-transfer for lines inside range (outside ilow, imax)
! version 1.3.3 Dec 2016: inclusion of depacked components as defined in LINES_xxx.dat
! version 1.3.4 Feb 2017: inclusion of Ng-extrapolation for selected transitions
! version 1.3.5 March/April 2017: lots of convergence control statements
! version 1.3.6 Dec 2017: subr. optcmf_restart changed: interpolation of XJ
!                         from previous results, if ifre1 ne ifre (change of ionization,
!                         due to different ilow, imax (subr. ilowimax_lte)
! version 1.3.7 March 2018: after including X-rays, two stop statements
!                         related to erroneous Edd. factors
!                         (subr. interpol_cmf_to_coarse, subr. xj_smooth) changed:
!                         check only for lambda > 20 A, otherwise numerical
!                         uncertainties too large. 
!                         Two similar stop statements not manipulated, since
!                         in effect only for CMF-range.  
!                         betacic (no upper level) in OPACITY_CMF now always
!                         calculated from consistent Edd. factor, and not
!                         approximative for L < LTHIN  
!                         (might be reconsidered if problems with certain lines transitions)
! 
! version 1.3.8 July 2018 final ALO control:
!                         expl. elements (calc_tlu)
!                                if ALO < 0.1 then set to zero (calc_tlu)
!                         background elements (calc_tlu_met)
!                                if ALO < 0.9 then set to zero,
!                                          except for OIII 374 (1-10)
!                                if ALO > 0.99, then set to 0.99
!
!                         subr. line list modified, to allow for changes
!                         in ionization structure during convergence  
!
! version 1.3.9 Feb 2019: new treatment of specific transitions with UP=0
!                         (see nlte.f90), if OPT_NEW_UPPER = .TRUE.
!                         This version has already been implemented in Nov. 2018!  
!
! version 1.4.0 Oct 2020: compatible with gfortran
!
! version 1.4.1 Dec 2020: first IR version
!
! version 1.4.2 April 2021: several updates, in particular regarding ALO of HeII 303
!
! version 1.4.3 March 2022: two stop statements after checking the
!                         eddington factor removed and replaced by warnings
!                         changes in cmf_complete: interpolate xj_save etc.
!                         when frequency grid has changed.  
!                         no stop for erroneous f-value (Adi) for O3_28 to O3_39
!  
character*5 :: ver_cmf_all='1.4.3' 
end module version_cmf_all
!
!--------------------------------------------------------------------------- &
!
module cmf_all_var
!
Use nlte_type      
Use nlte_dim, only: id_nttrd, id_atoms, id_ndept, id_llevs
Use nlte_app, only: natom
implicit none

! this needs to be adapted finally
!JO tested: for wavblue = 130 (NV->VI edge), no essential difference
real(dp) :: wavblue = 200.d0
!JO: at high temperatures, lower wavblue?
!    (NIV reson. lines around 196 A, NV reson. lines around 140 A)
!    seems to be OK; tested for D45 model, no changes

real(dp), parameter :: wavred = 10000.d0
!"allowed" values 10000,15000,20000,25000,50000
! certain values (e.g., 45000) can lead to problems in xh_obs_coarse
! or even to stops (kcmf_start1 ne kcmf_start), because of
! unsuited frequency grid


! impact of gf_min tested; if set to 1.d-5, no change
! might be different if improved Fe-model will be used
real(dp), parameter :: gf_min = 1.d-3 ! all lines with up=0 and gf < gf_min neglected
real(dp), parameter :: vturbmin = 5.d5
real(dp), parameter :: xmaxdop = 3.

! parameters for detecting significant line-overlap with explicit lines
real(dp), parameter :: xm = 1.d0      !overlap interval in units of vdop 
real(dp), parameter :: minrho = 0.3d0 !minimum opacity contribution of
                                      !considered explicit line 

! if set to true, following files will be written:
! out_xj_cmf
! out_xj_cmf_coarse
! out_xj_app
! out_xh_obs
! out_xh_obs_coarse
! NOTE: xh_obs (fine and coarse) multiplied by (r(1)/Rstar)^2 = RMAX^2
!       in contrast, XXH (from approx. solution) multiplied by (r(1)/SR)^2
logical, parameter :: optout_xj_xh = .true.

integer(i4b) :: ntotnl_all, ntot, ntot_expl_inside, icfirst, iclast
!ntot: all lines in new line list
!ntot_expl_inside: lines from expl. elem. inside wavblue, wavred in new line list

integer(i4b) :: nftot=0
!total number of frequency points within wavblue-dlam1, wavred+dlam2

integer(i4b) :: no2
!number of overlapping lines with up=0 as calculated in subr. overlap_noup

real(dp) :: vref,vturb_cmf

real(dp), dimension(id_nttrd) :: xlamcmf_depacked, sumgf_depacked
! analogous to xlamcmf, but with wavelength of strongest depacked component
integer(i4b), dimension(id_nttrd) :: indxlamc_depacked
! analogous to indxlamc (index for wavelength sorting), refering to xlamcmf_depacked

integer(i4b), dimension(id_nttrd) :: indexrbb_inv, indexcmf1, indexcmf1_lam
! indexrbb_inv gives index regarding index1 (rbb+cbb) for j=1,id_nttrd
! indexcmf1 refers to cmf-treatment (1,2,3,4,...)
! indexcmf1_lam points to index of lam_cmf
integer(i4b), dimension(id_atoms)  :: indexel_inv
!gives index of explict atom regarding background elements
integer(i4b), dimension(natom)  :: elements_cmf, nf_cmf
!provides info on which bg-elements are treated in cmf, and gives maximum number
!of freq. points (one sided, without 0) for Doppler profiles w.r.t. vref

real(dp), dimension(natom)  :: vth_cmf
!using vturb_cmf=max(vturb,vturbmin)
real(dp), dimension(id_ndept,natom)  :: vdoptot

real(dp), dimension(id_llevs,id_ndept) :: occexpl

integer(i4b), dimension(:), allocatable:: id1, id_rbb, index_id1
integer(i4b), dimension(:), allocatable:: id1_overlap
real(sp), dimension(:), allocatable :: gf1, xlam1, opal_ntest

integer(i4b), dimension(:), allocatable :: index_lam_cmf, indfre_cmf, indfre_staggered_cmf
!indfre_cmf: index regarding cont. grid (for interpolation)
!indfre_staggered_cmf: index regarding cont. grid (staggered)
real(dp), dimension(:), allocatable :: lam_cmf, h_obs, xhinner_cmf
real(dp), dimension(:,:), allocatable :: sumopal_cmf, sumetal_cmf, xj_cmf, alo_cmf
real(dp), dimension(:,:), allocatable :: xh_cmf,fedd, gedd, geddp
real(dp), dimension(:), allocatable :: hbound1,nbound1,hboundnd,nboundnd
real(dp), dimension(:,:), allocatable :: u_outer
real(dp), dimension(:), allocatable ::  sumopal_outer, sline_outer 
real(dp), dimension(:), allocatable ::  sumtaus_outer 

real(dp), dimension(:,:), allocatable :: profdop_H, profdop_He
real(dp), dimension(:,:), allocatable :: pweightdop_H, pweightdop_He
real(dp), dimension(:,:,:), allocatable :: profdop
real(dp), dimension(:,:,:), allocatable :: pweightdop

real(dp), dimension(:), allocatable :: sum_pweightdop_H, sum_pweightdop_He
real(dp), dimension(:,:), allocatable :: sum_pweightdop

real(dp), dimension(:), allocatable :: lamtest, weighttest

integer(i4b), dimension(:,:), allocatable :: ice
real(dp), dimension(:,:), allocatable :: icearr

integer(i4b) :: lineno

type depacked_data
  character*6 :: levlo
  character*6 :: levup
  integer(i4b) :: lineno
  integer(i4b) :: nocomp
  integer(i4b) :: index ! in subr. line_list, points to index in new line-list
                        ! in subr. fgrid and later, points to index of lamcmf
  real(dp) :: wave
  real(dp) :: flu
end type

type(depacked_data), dimension(:), allocatable :: depacked 

end module cmf_all_var
!
!--------------------------------------------------------------------------- &
!
module ng_var
!
Use nlte_type      
implicit none

! cut-off 
real(dp), parameter :: gfmin_ng = 1.d-3 
! elements (explicit and background) where Ng-acceleration should be performed
integer(i4b), parameter :: no_ng_el = 3
character*2, dimension(no_ng_el), parameter :: ng_el=(/'C','N','O'/)
!for tests with Ng also for HeII resonance line(s)
!integer(i4b), parameter :: no_ng_el = 4
!character*2, dimension(no_ng_el), parameter :: ng_el=(/'HE','C','N','O'/)
!integer(i4b), parameter :: no_ng_el = 1
!character*2, dimension(no_ng_el), parameter :: ng_el=(/'O'/)

! initialization: after restart etc., allow for one iteration without
! saving source functions, to obtain clear old and new occupation numbers
! THUS: itng=-1, itng_met=-1
! to avoid wrong n_ng (number of lines to be extrapolated) = 0 condition,
! initialize as negative
integer(i4b) :: itng=-1, n_ng=-1, itng_met=-1,n_ng_met=-1 

logical :: ng = .false., ng_met = .false.

integer(i4b), dimension(:), allocatable :: index1_ng, index_ng
integer(i4b), dimension(:,:), allocatable :: trans_ng
real(dp), dimension(:,:,:), allocatable :: sl_ng

integer(i4b), dimension(:), allocatable :: index_ng_met
integer(i4b), dimension(1000) :: trans_ng_met
real(dp), dimension(:,:,:), allocatable :: sl_ng_met



end module ng_var
!
!-----------------------------------------------------------------------
!
subroutine cmf_all(nd,xne,temp,clfac,r,velo,dvdr,rho,taur,xnh,pressure, &
                   ferr,dtfcorr,ilow,imax,rmax,unasol,opt_ray)
! 
! calculates formal integral in the cmf for *all* important lines, i.e.,
! for lines with CALCION(k,j) = 1, including also lines with only the lower level present
!
! all tests so far with GLOBAL = .TRUE.
USE nlte_type      
USE nlte_dim, only: id_ndept, id_atoms
USE nlte_opt, only: opt_modify_sl_noup
USE nlte_var, only: metals_converged, vsound, vmax, lwion1
USE nlte_app, only: natom, indexel,lniii_in,lniii_out

implicit none

integer(i4b), parameter :: nd1=id_ndept  
integer(i4b), parameter :: kel=id_atoms  

real(dp), parameter :: taur_test=0.1
logical, parameter :: test_overlap=.false.


!     ..
integer(i4b), intent(in) ::  nd
logical, intent(in) :: unasol
integer(i4b), dimension(nd1,kel), intent(in) :: ilow, imax
real(dp), intent(in) :: rmax
!     ..
!     .. array arguments ..
real(dp), dimension(nd1), intent(in) :: xne,temp,clfac,r,velo,dvdr,rho,xnh,pressure
!JO, changed April 2016
real(dp), dimension(nd1), intent(inout) :: taur
real(dp), dimension(nd1), intent(out) :: ferr,dtfcorr

integer(i4b) :: l,ntest,lsound

logical start, opt_ray

data start/.true./

print*
print*,'cmf_all started'

if (start) then

!JO Sept 2018  
! define region where NIII resonance lines are formed (when important)
  do l=1,nd-1
    if(velo(l)*vmax.lt.vsound) exit
  enddo
  lsound=l
  lniii_in=lwion1-1
  do l=1,nd-1
    if(velo(l).lt.0.1) exit
  enddo
  l=l-1
! use max(0.1 vinf,v(lsound-10)) as outer boundary 
  lniii_out=min(lsound-10,l)
  if(lniii_out.lt.1) stop ' lniii_out < 1'
  print*,'approx. range of NIII emission line-formation (if any):'
  print*,'lniii_out = ',lniii_out,' lniii_in = ',lniii_in
  print*
  
  call line_list(test_overlap)
! merge (explicit and background elements) and modify the original line-lists
! NOTE: ONLY REQUIRED AS LONG AS IMIN, IMAX ARE CHANGING (start=.true.)
  call fgrid_cmf
  call profile_cmf(nd,temp)
  if (metals_converged) start=.false.
endif

ntest=1
if(test_overlap) then
  do l=1,nd-1
    if (taur(l).le.taur_test .and. taur(l+1).gt.taur_test) exit
  enddo
  ntest=l
  if((taur(l+1)-taur_test).lt. abs(taur_test-taur(l))) ntest=l+1 
  print*
  print*,' ntest = ',ntest,' taur(ntest) = ',taur(ntest)
  print*
endif  
  
call opacity_cmf(nd,xne,temp,clfac,r,velo,dvdr,test_overlap,ntest,1)

if(test_overlap) then
  call overlap_detector(dvdr(ntest))
  stop ' stop after overlap_detector'
endif

!calculate and apply modified source-functions for bg-lines with no upper
!level, but only if non-HHe explicit elements are present
if (opt_modify_sl_noup .and. maxval(indexel(3:natom)).gt.-1) then
!detects lines with up=0 which need to be modified
  call overlap_noup(nd,clfac)

!2nd run, source-functions (and sumetal_cmf!) will be modified
!such 2nd run necessary to avoid a depth dependent array for id1_overlap
  call  opacity_cmf(nd,xne,temp,clfac,r,velo,dvdr,test_overlap,ntest,2)
endif

call cmf_complete(nd,xne,temp,clfac,r,velo,dvdr,rho,taur,xnh, &
                  pressure,ferr,dtfcorr,opt_ray)

! remap radiation field onto coarse grid (conserve freq. integrals),
! and overwrite radiation field (from CONT, coarse grid) by remapped
! cmf-quantities       
call interpol_cmf_to_coarse(nd)

! calculate obs. frame flux (approximate), write to file 'out_xh_obs',
! and resample onto coarse grid (flux-conservative)
if(unasol) then
  call formal_obs(r,velo(1),rmax)
  call interpol_hcmf_to_hcoarse(nd)
endif
  
return
end
!
!-----------------------------------------------------------------------
!
subroutine optcmf_restart
!
! provides info for cmf_all
!  
USE nlte_type
USE nlte_dim, ONLY: id_frec1, id_ndept

USE nlte_opt, ONLY: optcmf_full
USE nlte_var, ONLY: modnam, almost_converged,abund_have_changed, &
                    ifre,fre,xj_save,kcmf_start, kcmf_end, restart_cmfall, xj
                    
USE nlte_app, ONLY: met_imin, met_imax, occng
USE cmf_all_var, ONLY: wavblue, wavred
implicit none

!     .. parameters ..
integer(i4b), parameter :: nd=id_ndept  
integer(i4b), parameter :: ifretot=id_frec1  

integer(i4b) :: ifre1,i,j,k
real(dp), dimension(ifretot) :: fre1
real(dp), dimension(nd,ifretot) :: dummy

real(dp) :: eblue, ered, q, q1

logical :: newgrid

if (.not.optcmf_full) stop ' optcmf_full = F and optcmf_restart called'

open (1,file=trim(modnam)//'/OCCNG',status='old',form='unformatted')
rewind (1)
read(1) met_imin,met_imax,occng
read(1,end=10) almost_converged
goto 20

10 print*,' Restart: ALMOST_CONVERGED not defined (left at default)'

20 close(1)


if(abund_have_changed) almost_converged=.false.

print*,' Restart: file OCCNG successfully read!'
print*

!JO: added April 2016
open (1,file=trim(modnam)//'/CONT_FORMAL_ALL',status='unknown', &
& form='unformatted')

rewind 1  

! we read until XJ (first entries after FRE1 are OPAC, THOMSON, STRUE, XJ)
read (1) ifre1, (fre1(i),i=1,ifretot),((dummy(i,j),i=1,nd),j=1,ifretot), &
&  ((dummy(i,j),i=1,nd),j=1,ifretot),((dummy(i,j),i=1,nd),j=1,ifretot), &
&  ((dummy(i,j),i=1,nd),j=1,ifretot)
close (1)  

if(ifre1.ne.ifre) then
  print*,'old:',ifre1,' new:',ifre
  print*,'WARNING !!!! FREQUENCY GRID HAS CHANGED in optcmf_restart'
  newgrid=.true.
else
  newgrid=.false.
endif

if(.not.newgrid) then
  do i=1,ifre
    if(abs(1.d0-fre1(i)/fre(i)).gt.1.d-12) stop 'fre1 ne fre in optcmf_restart'
  enddo
!  if (dummy(1,ifre+1).ne.2) stop ' XJ qualifier ne 2 in optcmf_restart'
else
  if(abs(1.d0-fre1(1)/fre(1)).gt.1.d-12) stop 'fre1(1) ne fre(1) in optcmf_restart'
  if(abs(1.d0-fre1(ifre1)/fre(ifre)).gt.1.d-12) stop 'fre1(ifre1) ne fre(ifre) in optcmf_restart'
  if (dummy(1,ifre1+1).ne.2) stop ' newgrid and XJ qualifier ne 2 in optcmf_restart'
endif

allocate(xj_save(nd,ifre))
if(.not.newgrid) then
  do i=1,ifre
    xj_save(:,i)=dummy(:,i)
!for tests
!  print*,1.d8/fre(i),xj_save(1,i),xj(1,i)
  enddo
  print*,' Restart: XJ_SAVE (smoothed) successfully read from CONT_FORMAL_ALL!'
  print*
else ! interpolation of old results onto new grid
  xj_save(:,1)=dummy(:,1)
  xj_save(:,ifre)=dummy(:,ifre1)
outer: do i=2,ifre-1
    do k=1,ifre1-1
    if (fre(i).gt.fre1(k) .and. fre(i).le.fre1(k+1)) then
      q=log10(fre(i)/fre1(k))/log10(fre1(k+1)/fre1(k))
      q1=1.d0-q
      xj_save(:,i)=q1*log10(dummy(:,k))+q*log10(dummy(:,k+1))
      xj_save(:,i)=10.d0**xj_save(:,i)
      cycle outer                      
    endif                        
    enddo                        
    stop ' fre not found in fre1 (subr. optcmf_restart)'
  enddo outer
  print*,' Restart: XJ_SAVE (smoothed) interpolated from CONT_FORMAL_ALL!'
  print*
endif

if(newgrid) goto 40 ! kcmf_start and kcmf_end need to be defined

open(1,err=40,file=trim(modnam)//'/KCMF',status='old',form='formatted')
rewind 1
read(1,*) kcmf_start, kcmf_end
close (1)
goto 50

!file does not exist ('old' models); use wavblue and wavred to calculate indices
40 eblue=1.d8/wavblue
ered=1.d8/wavred

do i=1,ifre
  if(fre(i).gt.ered) exit
enddo
kcmf_start=i

do i=ifre,1,-1
  if(fre(i).lt.eblue) exit
enddo
kcmf_end=i

if(abs(1.-fre(kcmf_start)/ered).gt.0.1 .or. abs(1.-fre(kcmf_end)/eblue).gt.0.1) &
  stop ' Problems with kcmf_start or kcmf_end in optcmf_restart'


50 print*,' Restart: KCMF_START/END = ',kcmf_start,kcmf_end
print*,' corresponding wavelengths:',1.d8/fre(kcmf_start),1.d8/fre(kcmf_end)
print*

restart_cmfall = .true.
print*,' RESTART_CMFALL set to .true.'
print*

return
end
!
!-----------------------------------------------------------------------
!
subroutine line_list(test_overlap)
!
! merge (explicit and background elements) and modify the original line-lists
!
! full coupling between wavblue and wavred
! beyond, cmf transfer with approx. background radiation field
! Thus line list:
!           lam < wavblue:   only cmf lines from explicit elements 
! wavblue < lam < wavred :   important lines from background,
!                            plus cmf lines from explicit elements 
!  wavred < lam < xlam(ired):only cmf lines from explicit elements 
!
!
! philosophy regarding tests for indexcmf might need to be changed for
! optmixed-approach
!
!JO -- IMPORTANT NOTE: in current version, also depacked data from
!  LINES_xxx.dat included 
!
!NOTE regarding bg-elements:
! please remember that any manipulation here (id1, xlam1, gf1) will affect
! 'only' the lines (opacities and emissivities) as treated in the
! cmf-transport, particularly the resulting J, Jbar and ALO.
! With respect to rates (as calculated in calc_tlu_met),
! the contributing components and data (particularly gf) are taken from
! the original files (id,xlam,gf) and stored in met_rbb. This is done
! in subr. fgrid_cmf. To manipulate any data consistently,you have to
! hack id_rbb created below (and maybe some data in subr. fgrid_cmf)
!
!
!  'CaVI-problem' (March 2016)
!  During our tests, it turned out that for hot models the convergence of
!  the T-struct (and rest) is hampered by the interaction of the OV 1-3 line
!  at 629.73 with the Ca VI line 1-4 at 629.49 (this line has three strong
!  components, at 629.49, 634.058 and 641.66), Since the cooling/heating
!  balance is quite strongly influenced by the CBB transitions from OV 1-3
!  (30% of the total heating rate for model d2v), this transition is very
!  important. What happens is the following: when Ca VI overlaps with OV,
!  the cooling becomes strong because the upper level 3 is no longer as
!  populated. Below a certain temperature, Ca VI looses its strength,
!  OV can pump, and the third level can heat. Then the cycle starts again
!  (more CaVI, less pumping, more cooling, etc.)
!  Thus, for specific conditions we exclude the 3rd component of CA VI
!  at 629.49 from the line list, and the temp. convergence is no longer
!  hampered.
!
USE nlte_type      
USE nlte_dim
USE fund_const, only: clight
USE princesa_var, only: zeff,weight,abund,gl,fl,zl,glx,flx,gs0,gs1,ryd,qd, &
&       data,ndat,le, &
&       nlong,indtt,index1,index2,index3,index4,index5, &
&       indtr,iaux2,iform1,iform2,iform3,iform4,iform5,numd1,numd2, &
&       numd3,numd4,numd5,indat1,indat2,indat3,indat4,indat5,iaux1, &
&       labl,labl4,labl1,labu1,labl3,ipare5, &
&       li,ifirsl,nat,labat
USE ffr_error

USE nlte_var, only: indxlamc,xlamcmf,indexrbb,mmlow,mmup,lablinelo,lablineup, &
&                   indexcmf,restart_cmfall

USE tcorr_var, only: enatcor, temp_converged

USE nlte_app, only: jatom,calcion,indexel,metall,jatom_full,indrec,irbb,ilin, &
                    met_rbb,natom,names1,indexel,teff
USE nlte_lines, only: id, gf, xlam, iblue, ired, ntotnl3
USE cmf_all_var, only: wavblue, wavred, gf_min, indexrbb_inv, ntotnl_all, &
&                   id1, gf1, xlam1, ntot, ntot_expl_inside, elements_cmf, &
&                   indexcmf1, id_rbb, indexel_inv, opal_ntest, index_id1, &
&                   id1_overlap, lineno, &
&                   depacked, xlamcmf_depacked, indxlamc_depacked, sumgf_depacked

implicit none

!maximum separation of lines (in cm/s) for overlap detector
real(dp), parameter :: deltavmax=150.d5 

real(dp) :: xxlam,lam,flu,xlamold,xl,deltav,flu1,sumf,maxf
integer(i4b) :: ii,indi,j,k,kk,l,ml,mu,irest,low,up,ic,jc,jce,index, &
&               jce_icmf1,jce_red,jce_outside,n,jj,ji, &
&               igenio, kold, i, jfull, jfullold, &
&               iu, iou, nocomp, ic1, ic2, icomp

logical first, newj, no_calcium_six, test_overlap, found

character*1 :: met
character*4 :: ret
character*6 :: lab, dc, level, leveu
character*32 :: fichero, stafil, linetab

data first/.true./

print*,'number of lines from original nl3-files:',ired-iblue+1

if (id_nttrd .ne. indtr) stop ' something wrong with id_nttdr (subr line_list)'

if (first) then
!
lineno=0
!
! -----------------------------------------------------------
! prepare depacking, but only if non-HHe elements are present
if (maxval(indexel(3:natom)).gt.-1) then

! use LINES_xxx.dat (as encoded in linetab) to depack specific transitions
! from explicit elements
   open (1,file='ATOM_FILE',status='old')  
   rewind 1  
   read (1,fmt='(a)') fichero  
   read (1,fmt='(a)') stafil  
   read (1,fmt='(a)') linetab  
   close (1)  

   iu = 37  
!   iou = 38  
   dc = ':T'  

   print* 
   print*,'file ',trim(linetab),' used for depacking transitions from explicit elements'   
     
   open (unit=iu,file='../inicalc/DATA/'//linetab,form='formatted',status='old')  
!   open (unit=iou,file='control_lines.dat',status='unknown')  
!
!   rewind!!!
!
   iu = -37  
   call ffracc(iu,ret)  
   if (on_error(ret))  goto 60

   call transic_count(lineno,ret)
   print*,lineno,' (de-packed) lines present in file'
   if (on_error(ret))  goto 70

   allocate(depacked(lineno))
   
   rewind(abs(iu))
   call transic_read(lineno,ret) !wavelengths (lambda > 2000 A) converted
   if (on_error(ret))  goto 70
! for tests
!   do i=1, lineno
!     print*,depacked(i)%levlo, depacked(i)%levup, depacked(i)%wave, depacked(i)%lineno, depacked(i)%flu
!   enddo
   write (*,fmt='(a)') ' data successfully read! '
   write (*,fmt='(a)') ' consistency of statistical weights checked '
   write (*,fmt='(a)') ' wavelengths (lambda > 2000 A) converted from air to vacuum'
   print*
   close (iu)  
!   close (iou)  
else
!JP contemporary hack, to ensure that "depacked" is allocated even for
! HHe as explicit elements
  if (lineno.ne.0) stop 'HHe and lineno(depacked) ne 0'
   allocate(depacked(lineno))
endif !prepare depacking
!
! -----------------------------------------------------------
!
if(wavblue.lt.xlam(iblue)) then
! wavblue=xlam(iblue)+20.
! changed by JO Aug. 2013; note cooler stars with lines
! only until 228 A (Teff < 25kK) or even 250 A (Teff < 10kK)
  wavblue=xlam(iblue)+2.
  print*,' wavblue (cmf_all) reset to ',wavblue
endif
  
if( wavred.gt.xlam(ired))  stop ' wavred (cmf_all) > wavred(nlte_approx)'

do j=1,id_nttrd
  ii=indxlamc(j)
  do k = 1,index1
    if(indexrbb(k).eq.ii) then
      indexrbb_inv(j)=k
      goto 10
    endif  
  enddo
  stop 'ii not found in indexrbb (subr. line_list)'
  
10 continue
enddo

ntotnl_all = ntotnl3 + id_nttrd

! this has to be done only once, since dimension fix
allocate(id1(ntotnl_all),gf1(ntotnl_all),xlam1(ntotnl_all))
allocate(id1_overlap(ntotnl_all))
if(test_overlap) allocate(opal_ntest(ntotnl_all),index_id1(ntotnl_all))

! index_array providing index of explict atom w.r.t. background elements
! shifted from subr. profile_cmf
if (id_atoms .ne. nat) stop ' id_atoms ne nat (subr. line_list)'
do k = 1,nat
  do i = 1,natom
    if(labat(k).eq.names1(i)) then
      indexel_inv(k)=i
    endif  
  enddo
!test
  if (indexel(indexel_inv(k)).ne.k) stop ' something wrong with index_el(inv), H or He not explicit? (subr. line_list)'
enddo

! determine the index of met_rbb related to bg line-list within wavblue, wavred.
! Only for selected NLTE bg-elements (jatom_full=1) and lines with defined upper level, 
! but for ALL ionization stages. Thus, this task needs to be performed only once.
allocate(id_rbb(ntotnl3))
id_rbb=0

print*
print*,'significant difference (>20%) in wavelength (bg-elements, line list vs. packed transition)!'

lines_rbb: do l=iblue,ired
   if (xlam(l).lt.wavblue) cycle lines_rbb
   if (xlam(l).gt.wavred) exit lines_rbb

   k = id(l)/1000000
   if(jatom_full(k).eq.0) cycle lines_rbb

   irest = id(l) - k*1000000
   j = irest/100000
   irest = irest - j*100000
   low = irest/100
   up  = irest-low*100
   if(up.eq.0) cycle lines_rbb    !.and. gf(l).lt.gf_min) ! if we want to cut the list

   n=indrec(k,j)
   jj=irbb(n)-1

   do ii=1,ilin(n)
     ji=jj+ii
     if (low .eq. met_rbb(ji)%low .and. up .eq. met_rbb(ji)%lup) exit
   enddo
   if (abs(1.-met_rbb(ji)%wave/xlam(l)) .gt. 0.20) then
     write(*,fmt='(4(2x,i2),2(2x,f10.3))') k,j,low,up,xlam(l),met_rbb(ji)%wave
   endif  

   id_rbb(l)=ji
! hack to exclude certain components
!   if(id(l).eq.8300110 .and. xlam(l).ne.374.073) id_rbb(l)=0
!   if(id(l).eq.8300106 .and. xlam(l).ne.702.337 .and. xlam(l).ne.703.855) id_rbb(l)=0
!   if(id(l).eq.8300105 .and. xlam(l).ne.835.289 .and. xlam(l).ne.832.929) id_rbb(l)=0

enddo lines_rbb
print*

first=.false.
endif

!JO July 2018 included to update indexrbb_inv in case
if(restart_cmfall) then
do j=1,id_nttrd
  ii=indxlamc(j)
  do k = 1,index1
    if(indexrbb(k).eq.ii) then
      indexrbb_inv(j)=k
      goto 15
    endif  
  enddo
  stop 'ii not found in indexrbb (subr. line_list)'
  
15 continue
enddo
print* 
print*,'indexrbb_inv updated'  
print*
endif

! local copy of indexcmf:
! lines treated in cmf_complete will obtain an indexcmf1=3 or higher (for depacked levels)
! lines outside range will keep their previous status
!
indexcmf1=indexcmf

! initialize
depacked%index=0

ic=0 ! counter for new list
jc=0 ! counter refering to explicit elements
jce=0! counter for cmf-lines from explicit elements in new list
jce_outside=0 !counter for cmf-lines from explicit elements outside wavblue,wavred
jce_icmf1=0   !counter for lines from expl. elements with index_cmf = 1 (outside imia,imaa)
jce_red=0     !counter for cmf-lines from explicit elements with lam > lam (ired)
newj=.true.
elements_cmf=0

!do k=1,30
!  write(*,5) k,(calcion(k,j),j=1,9)
!enddo  
!5 format(10(i2,2x))

!check whether problematic CA VI line should be removed
!only for hot models and when temperature is not converged
!JO: in case, also other lines might need to be removed;
!    check with overlap detector at end of this routine

no_calcium_six=.false.
if(teff .gt. 40000.  .and. enatcor .and. .not. temp_converged) no_calcium_six=.true.

!at first, remove all lines with calcion = 0, wave < wavblue, and explicit elements
!then include all lines from explicit elements
lines: do l=iblue,ired
   if (xlam(l).lt.wavblue) cycle lines
! this is the hack to avoid lines overlapping with HeI sing resonance line
!   if (xlam(l).gt.580.and.xlam(l).lt.590.) cycle lines
! this is the hack to exclude most lines...
!   if (xlam(l).lt.6000.)  cycle lines
!   print*,id(l),xlam(l)
   k = id(l)/1000000
   if(jatom(k).eq.0) cycle lines

!for tests without non-selected elements
!   if(jatom_full(k).eq.0) cycle lines

   if (indexel(k).gt.0) cycle lines  ! included in nlte
   
   irest = id(l) - k*1000000
   j = irest/100000

! This is the FeV line (374.245) which is also responsible for the NIII emission
! note: delta v = 37 km/s, i.e., for small turbulence, there should be no effect.   
!   if (xlam(l).gt.372.and.xlam(l).lt.376. .and. id(l).eq.26500200) cycle lines

   
   if(calcion(k,j).eq.0) cycle lines !too low abundance

   irest = irest - j*100000
   low = irest/100
   up  = irest-low*100
! for tests
!   if(k.eq.8 .and. up.eq.0) cycle lines
   
!for tests without OIII line overlap
!   if(k.eq.8 .and. j.eq.3 .and.low.eq.1 .and. up.eq.10) then
!     if(xlam(l).ne.374.073) cycle lines
!     print*,'O3 found',xlam(l)
!   endif  

!   if(k.eq.8 .and. j.eq.3 .and.low.eq.1 .and. up.eq.6) then
!     if(xlam(l).ne.702.337 .and. xlam(l).ne.703.855) cycle lines
!     print*,'O3 found',xlam(l)
!   endif  
!   if(k.eq.8 .and. j.eq.3 .and.low.eq.1 .and. up.eq.5) then
!     if(xlam(l).ne.835.289 .and. xlam(l).ne.832.929) cycle lines
!     print*,'O3 found',xlam(l)
!   endif  

   if(up.eq.0 .and. gf(l).lt.gf_min) cycle lines ! only stronger resonance lines

   n=indrec(k,j)
     
! test overlap with N48
!   if(up.eq.0 .and. gf(l).lt.10. .and. xlam(l).eq.387.371)) cycle lines ! overlap with FeV  26500400
   
! no check for met_imin, met_imax required, since no problem with
! opacities missing at certain depth-points

!  if only selected elements, uncomment following line
!   if(jatom_full(k).eq.0) cycle lines

! overlap between CaVI and OV (629.496 and 629.732) prevents
! temperature convergence for hot models (e.g., d2v).
! remove responsible CA VI line in CMF transport 
! JO: when data-base updated, check frequencies 
! JO: other lines might need to be removed as well
   if (no_calcium_six .and. k.eq. 20) then
       if (j.eq.6 .and. low.eq.1  .and. xlam(l).gt.629. .and. xlam(l).lt.630.) then
       print*,'CA VI line at',xlam(l),' removed from cmf transport'    
       print*
       cycle lines
     endif  
   endif
     
   met=metall(k,j,low)

   if(low.ne.1.and.met.eq.' '.and.jatom_full(k).eq.0) cycle lines ! consistent with sumopal

! exclude lines overlapping with NV 2->3
!   if (xlam(l).gt.250.and.xlam(l).lt.280.) cycle lines

   
   if(up.ne.0) then
     met=metall(k,j,up)
     if(met.eq.' '.and.jatom_full(k).eq.0) cycle lines ! no approx. upper level
   endif

!  for tests
!   if(k.eq.15.and.j.eq.3.) print*,k,j,low,up,xlam(l),gf(l)
!
   elements_cmf(k)=1 ! element k treated in cmf

! JO changed Oct. 2016
! remove Fe-lines (with uncertain gf values) overlapping with HeI resonance line
! JO: other lines might need to be removed as well if there are problems
! in the outer region, in particular the resonance line of MN IV (584.296 A). 
   
! JO: when data-base updated, check frequencies 
   if(xlam(l).ge.584.23 .and. xlam(l).le.584.43) then
! FeIV 584.368/584.397 with up=0
     if(k.eq.26 .and. j.eq.4 .and. (low.eq.6 .or. low.eq.8)) gf(l)=1.d-5
   endif  

! for tests 
!   if(up.ne.0) then
!     print*,xlam(l),jatom_full(k),id(l),metall(k,j,low),' ',metall(k,j,up)
!   else
!      print*,xlam(l),jatom_full(k),id(l),metall(k,j,low),' 0'
!   endif
!
! info on line-transitions from explict elements
20 if (newj) then
     jc=jc+1
! assume that there are lines from explicit elements with lam > xlam(ired)
     if (jc.gt.id_nttrd) stop ' jc > id_nttrd in line_list'
     ii=indxlamc(jc)
     xxlam=xlamcmf(ii)
     if (xxlam.eq.0.) then
       if (indexcmf(ii).ne.1) stop ' xxlam = 0 and indexcmf(ii) ne 1 in line_list'
       jce_icmf1=jce_icmf1+1
       goto 20 !from ions outside imia,imaa
     endif
     if (indexcmf(ii).ne.2 .and. xxlam .lt. xlam(ired)) then
!JO: changed Dec. 2016: sobo lines inside range can be present
!    if occup. numbers not defined at ALL depth points 
       print*,' line at ',xxlam,' treated in Sobo (outside ilow, imax)'
       jce_icmf1=jce_icmf1+1
       goto 20 !from ions outside imia,imaa
!       stop 'indexcmf(ii) ne 2 in line_list'
     endif  
! this is the end condition
     if (xxlam .ge. xlam(ired)) exit lines
     kk=indexrbb_inv(jc)
! kk is the index w.r.t. index1
     indi=indat1(kk)
     lam=data(indi+1)
     if(lam .ne. xxlam) stop ' lam ne xxlam (subr. line_list)'
     flu=data(indi+2)
     ml=mmlow(kk)
     mu=mmup(kk)
!for tests, works
!   if(lablinelo(ii).ne.labl(ml)) stop ' error in label(low) in line_list'
!   if(lablineup(ii).ne.labl(mu)) stop ' error in label( up) in line_list'
   endif

   if (xlam(l).lt.xxlam .and. xlam(l).le. wavred) then
     ic=ic+1
     id1(ic)=id(l)
     gf1(ic)=gf(l)
     xlam1(ic)=xlam(l)
     newj=.false.
   else if (xxlam .lt. xlam(ired)) then
! lines from explicit elements with lam < wavblue included as first lines
! lines from explicit elements with lam > wavred included as last lines
!
! leading '1' to tell that this is a line from explicit elements (different coding) 
     index=le(ml)*1d8+ml*1d4+mu
     ic=ic+1
     jce=jce+1
     id1(ic)=index
     gf1(ic)=gl(ml)*flu
     xlam1(ic)=lam
     newj=.true.
     if(lam.lt.wavblue .or. lam.gt.wavred) then
       jce_outside=jce_outside+1
     else
!    considered line will be treated in cmf_complete
       where(depacked%levlo .eq. labl(ml) .and. depacked%levup .eq. labl(mu))
! note: only transitions with xlam ne 0 (inside imia,imaa)
         depacked%index=ic
       endwhere
       indexcmf1(ii)=3
     endif

!     print*,ic,l,jc,xlam1(ic-1),xlam1(ic),id1(ic)
     goto 20

   endif  

enddo lines

! last lines from explicit elements (not used) should be sobo-lines
do j=jc,id_nttrd
  ii=indxlamc(j)
  if (indexcmf(ii) .ne. 1) stop ' something wrong in philosophy -- line_list'
enddo
jce_red=id_nttrd-jc+1

ntot=ic
ntot_expl_inside=jce-jce_outside
print*,'number of lines in merged files: ',ntot
print*,'including ',jce,' cmf lines from explicit elements'
print*,'number of cmf lines (expl. elem.) inside wavblue, wavred:',ntot_expl_inside
print*,'number of lines (expl. elem.) outside imia,imaa:', jce_icmf1
print*,'number of Sobolev lines (expl. elem.) with lam > lam(ired):',jce_red
print*
if (jce+jce_icmf1+jce_red.ne.id_nttrd) &
& stop ' not all lines from expl. elements found in subr. line_list'

!test that all explicit lines have been found
ic=0
do j=1,id_nttrd
  ii=indxlamc(j) 
! for specific tests
!  print*,j,xlamcmf(ii),indexcmf1(ii)
  if (indexcmf1(ii).eq.3) ic=ic+1
enddo

if(ic.ne.ntot_expl_inside) &
& stop ' not all explicit lines inside wavblue, wavred identified in indexcmf1'

!for tests
!do j=1,ntot
!  if(xlam1(j).gt.600. .and. xlam1(j).lt.650.) &
!&   print*,j,xlam1(j),id1(j),gf1(j)
!enddo
!stop

!exchange/add depacked components in original data
!works also for lineno = 0
ic=ntot
ic1=0
ic2=0

! initialize
xlamcmf_depacked = xlamcmf
sumgf_depacked = 0.
depacked%nocomp = 0

if(lineno .gt. 0) then
do j=1,id_nttrd
     ii=indxlamc(j) 
     xxlam=xlamcmf(ii)
     if(xxlam.eq.0.) cycle !outside imia,imaa
     kk=indexrbb_inv(j)
! kk is the index w.r.t. index1
     ml=mmlow(kk)
     lab=labat(le(ml))
     if(lab.eq.'H' .or. lab.eq.'HE') cycle !only non-HHe elements
     if(indexcmf1(ii).ne.3) cycle !only lines with cmf_all treatment
     mu=mmup(kk)
     level=labl(ml)
     leveu=labl(mu)
     indi=indat1(kk)
     lam=data(indi+1)
     if(lam .ne. xxlam) then
       print*,lam,xxlam
       stop ' lam ne xxlam (subr. line_list)'
     endif
     flu=data(indi+2)
!     print*,lam,level,' ',leveu,flu
!     rewind(iou)

! check consistency of f-values and provide number of de-packed components 
! thus far, nocomp not required
     sumf=0.d0
     nocomp=0
! note: where statement requires additional array(s), e.g.,
! idum=0
! where (depacked/levlo .eq. ...
!  idum=depacked%lineno
!  depacked%nocomp=maxval(idum)
! endwhere
     maxf=0.d0
     icomp=0
     do i=1,lineno
       if(depacked(i)%levlo .eq. level .and. depacked(i)%levup .eq. leveu) then
         nocomp=nocomp+1
         flu1=depacked(i)%flu
         sumf=sumf+flu1
! determine largest f-value of depacked transition
! NOTE: using amax1 destroys precision
         maxf=max(maxf,flu1)
         if(maxf.eq.flu1) icomp=depacked(i)%lineno
       endif
     enddo  
     
     where(depacked%levlo .eq. level .and. depacked%levup .eq. leveu)
       depacked%nocomp=nocomp
     endwhere

     do i=1,lineno
       if(depacked(i)%levlo .eq. level .and. depacked(i)%levup .eq. leveu) then
         if(icomp.eq.0) stop ' error in icomp'
!         print*,'found ',level,leveu,depacked(i)%lineno,depacked(i)%wave,depacked(i)%nocomp,depacked(i)%index
! output only once
         if(abs(1.-sumf/flu).gt.0.2 .and. icomp.eq.depacked(i)%lineno) then
           print*,'inconsistent f-values for ',level,' ',leveu
           print*,'sum(comp):',sumf,' packed:',flu
           if(level.eq.'O3_28' .and. leveu.eq.'O3_39') then
             print*,"erroneous f-value ('packed') in WM-basic database, needs to be updated"
           else
             print*,'check components and flu values in ',trim(linetab)
             stop ' inconsistent f-values for transition from explicit element'
           endif
           endif
         if(depacked(i)%lineno .eq. 1) then
! exchange wavelength and gf value
           k=depacked(i)%index             
           xlam1(k)=depacked(i)%wave
           gf1(k)=depacked(i)%flu*gl(ml) !gf value ; remember: gf(individual) = flu(input) * gl(tot)
           ic1=ic1+1
         else  
! add info for additional components
           ic=ic+1
           ic2=ic2+1
           if(ic .gt. ntotnl_all) stop ' ic too large' 
           k=depacked(i)%index             
           id1(ic)=id1(k)
           xlam1(ic)=depacked(i)%wave
           gf1(ic)=depacked(i)%flu*gl(ml) !gf value
         endif  

         if(icomp.eq.depacked(i)%lineno) then
! change wavelength in xlamcmf_depacked corresponding to wavelength of
! transition with largest f-value
           xlamcmf_depacked(ii)=depacked(i)%wave
! if different from zero, indicates that transition depacked/changed
           sumgf_depacked(ii)=sumf*gl(ml)
! overwrite indexcmf1 with number of components (when different from 1): 
! indexcmf1 = 4 means two components, = 5 means three components, and so on
           if(indexcmf1(ii).ne.3) stop ' something rotten with indexcmf1'
           indexcmf1(ii)=2+depacked(i)%nocomp ! works also if only one component
         endif

       endif
     enddo
enddo

!create new indexarray for xlamcmf_depacked
call indexx(id_nttrd,xlamcmf_depacked,indxlamc_depacked)

if(ic2 .ne. ic-ntot) stop 'error in exchange/add of components'
print*
print*,ic1,' lines/components exchanged (from ',trim(linetab),')'
print*,ic2,' components added (from ',trim(linetab),')'
print*

ntot=ic
jce=jce+ic2
ntot_expl_inside=ntot_expl_inside+ic2
print*,'FINAL number of lines in merged files: ',ntot
print*,'including ',jce,' cmf lines from explicit elements'
print*,'number of cmf lines (expl. elem.) inside wavblue, wavred:',ntot_expl_inside

!final test that %lineno correct (for all used lines, i.e., with %nocomp >0)
do i=1,lineno
  if(depacked(i)%nocomp .gt. 0 .and. depacked(i)%lineno.gt.depacked(i)%nocomp) then
    print*,depacked(i)%levlo,depacked(i)%levup,' : ',depacked(i)%lineno
    print*,'number of components: ',depacked(i)%nocomp
    stop ' something wrong with %lineno in LINES_xxx.dat file, please correct'
  endif
enddo  

!re-sort line list
if(size(xlam1).ne.ntotnl_all) stop ' error in nsize (line_list)'
! sort according to wavelength
call sort_line_list(ntot,xlam1(1:ntot),gf1(1:ntot),id1(1:ntot))

else if(lineno.eq.0.) then

indxlamc_depacked=indxlamc  

else
stop ' lineno <0!!!'
  
endif

!for tests
!do i=1,ic2
!  j=ntot-ic2+i
!  print*,id1(j),xlam1(j),gf1(j)
!enddo  

!for tests
!do i=1,ntot
!  if(xlam1(i).gt.3400.) print*,xlam1(i),gf1(i),id1(i)
!enddo
!do j=1,id_nttrd
!     ii=indxlamc_depacked(j)
!     i=indxlamc(j) 
!     print*,j,xlamcmf_depacked(ii),indexcmf1(ii),xlamcmf(i)
!enddo     

return
!in case, return here

!overlap detector
!detects overlap of important (resonance lines), within prescribed deltav
!E.g., it detects overlap between CaVI and OV (629.496 and 692.732) that prevents
!T-conv for model d2v (with 'standard' Delta n = 5 and gfmin = 0.1).
!If Delta n = 10 and gfmin = 0.05, it also detects OIII/NIII overlap
!responsible for OIII emission at late O-types.
!JO: needs to be tested when N/O treated as explicit elements 
kold=0
do ii=1,ntot
        index=id1(ii)
        if (index.ge.1d8) then
! explicit element
          kk = index/1d8
          irest = index - kk*1d8
          ml = irest/10000
          mu = irest - ml*10000
          ji=igenio(le(ml),li(ml)) !general counting
          k=indexel_inv(kk)
          j=int(zl(ml)+1.d-6)+1 ! in astro-definition (1 = neutral)
! standard DELTA n and gfmin
          if(ml.eq.ifirsl(ji).and.mu-ml.le.5.and.gf1(ii).gt.0.1) then !important resonance line
            xl=xlam1(ii)
            if(kold.eq.0) then
              write(*,fmt='("explicit   ",4i3,2x,f10.4)') k,j,ml,mu,xl
              xlamold=xl
              kold=k
            else 
              deltav=(xl-xlamold)/xl*clight
              if(k.ne.kold .and. deltav.lt.deltavmax) then ! always 
                write(*,fmt='("explicit overlap  ",4i3,2x,f10.4)') k,j,ml,mu,xl
              else
                write(*,fmt='("explicit   ",4i3,2x,f10.4)') k,j,ml,mu,xl
              endif
              xlamold=xl
              kold=k              
            endif
            jfullold=1
          endif
        else        
          k = index/1000000
          irest = index - k*1000000
          j = irest/100000
          irest = irest - j*100000
          low = irest/100
          up  = irest-low*100
! standard DELTA n and gfmin
          if(low.eq.1.and.up.ne.0.and.up-low.le.5.and.gf1(ii).gt.0.1) then !important resonance line
            xl=xlam1(ii)
            if(kold.eq.0) then
              write(*,fmt='("background ",4i3,2x,f10.4)') k,j,low,up,xl
              xlamold=xl
              kold=k
              jfullold=jatom_full(k)
            else 
              deltav=(xl-xlamold)/xl*clight
              jfull=jatom_full(k)
              if(k.ne.kold .and. (jfullold.eq.1 .or. jfull.eq.1) .and. &
                 deltav.lt.deltavmax) then ! only for explict and selected elements 
                write(*,fmt='("background overlap",4i3,2x,f10.4)') k,j,low,up,xl
              else
                write(*,fmt='("background ",4i3,2x,f10.4)') k,j,low,up,xl
              endif
              xlamold=xl
              kold=k              
              jfullold=jfull
            endif
          endif
        endif
enddo

!
!       error in input unit exit
!
60 continue  
write (*,fmt='(a)') ' error in input unit number '  
stop

!       error conditions
70 select case(ret)
case ('RET1') !       end of file "no esperado" exit  
  write (*,fmt='(a)') ' eof not expected. there is something wrong! '  
  stop
case ('RET2') !       error in ffr subroutines exit  
  write (*,fmt='(a)') ' error in ffr subroutines '  
  stop
case default
  stop ' wrong error condition in transic'
end select

end
!
!-----------------------------------------------------------------------
!
subroutine fgrid_cmf
!
! grid for coupled cmf-transfer; equidistant in dlam/lam = vref/clight
! vref=vturb/3, with vturb = max(vturb,vturbmin)
! i.e., if vturb = 0, artificial turb. of vturbmin (5 km/s) is used, 
! to avoid too large freq. grid  
!
! after some tests, we decided to use vturb to obtain a reasonable
! resolution for all elements at all depth points.
!
! To safe time, one might also use vth(O), but then the resolution of
! lines from heavy elements is below optimum  
!
! wavelength shifted to closest grid points.
!
! Note: central line frequencies can be shifted by +/- 0.5 dlambda,
! corresponding to vref/2 = vturb/6
! For max(vturb) approx. 15 km/s, this corresponds to +/- 2.5 km/s.
! This uncertainty is reflected in any formal integral solution based
! on the freq. grid provided here.

! the scaling of the number of freq. points, nftot, is via
! 1/ln(1+vref/c) approx 1/vref  

USE nlte_type      
USE nlte_dim
USE fund_const, only: clight
USE nlte_app, only: vth, met_rbb, natom,jatom_full,jmax,indrec,irbb,ilin, &
&                tlu_met, tul_met
USE nlte_var, only: vturb,fre,ifre,lwion1
USE nlte_lines, only: id, gf, xlam, ntotnl3, unpacked_index_cmf, unpacked_gf, &
&                unpacked_xlam
USE cmf_all_var, only: wavblue, wavred, gf_min, xlam1, ntot, &
&                vturbmin, vturb_cmf, vref, &
&                lam_cmf, index_lam_cmf, indfre_cmf, indfre_staggered_cmf, nftot, &
&                icfirst, iclast, id1, id_rbb, lineno, depacked

implicit none

integer(i4b) :: i,ic,in,ic0,din,din1,lastindex,dinmax,j,istart,ii, &
&               ind,inds,nftotnew,max_no_comp,k,n,jj,tot_comp,ji, &
&               iqual0,iqual1,iqual2,irun,no_comp,index_comp,nsize

real(dp) :: vend,dlp,dlm,lmin,lmax,lam,energy,frec,diffe,diffe1,sum_gf, &
                cln,wave,lamax

real(dp), save :: lammin,lammax

logical first

data first/.true./

!--------------------------------------------------------------------
!
!set up freq. grid and various index files only once
if(first) then

if(nftot.ne.0) stop ' first and nftot ne 0 in fgrid_cmf'

if(vturb.lt.vturbmin) then
!if(vturb.ne.vturbmin) then
print*,'WARNING!!!!WARNING!!!!WARNING!!!! artificial vturb = ',vturbmin/1.d5, &
& 'km/s used in cmf_all'
vturb_cmf=vturbmin
!vref=sqrt(vth(8)**2-vturb**2+vturb_cmf**2) !gives less freq. and profile points
vref=vturb_cmf
else
!vref=vth(8)
vturb_cmf=vturb
vref=vturb_cmf
endif

vref=vref/3. ! 3 points per vth is sufficient
! note: CMFSING works with NF = 21 points for 6-8 vth
print*
print*,'reference velocity used for freq. grid = ',vref/1.d5,' km/s'
print*,'max. error due to shift of rest-wavel. is +/- ',vref/2.d5,' km/s,' 
print*,'corresponding to dlambda = ',vref/2./clight*4500,' A at 4500 A'

vend=sqrt(20.d5**2+vturb_cmf**2)*5
             ! typical value for vth(H) -- on blue side, even factor 2 smaller - no H, but He --, &
             ! assuming strong line 5 vth wide
lammin=wavblue*(1.-vend/clight)
lammax=wavred* (1.+vend/clight)

! verify that nftot does not change
nftotnew=log(lammax/lammin)/log(1.d0+vref/clight)+1

nftot=nftotnew

! allocate only once, since nftot should not change during iteration
allocate(lam_cmf(nftot), index_lam_cmf(nftot), indfre_cmf(nftot),indfre_staggered_cmf(nftot))

cln=(1.d0+vref/clight)
do i=1,nftot
  lam_cmf(i)=lammin*cln**(i-1)
enddo

! define index-file w.r.t. continuum frequency grid

energy=1.d8/lam_cmf(1)  ! lam_cmf from blue to red, fre vice vera
if (energy.gt.fre(ifre)) stop ' lam_cmf(blue) outside fre! (subr. fgrid_cmf)'
energy=1.d8/lam_cmf(nftot)  ! lam_cmf from blue to red, fre vice vera
if (energy.lt.fre(1)) stop ' lam_cmf(red) outside fre! (subr. fgrid_cmf)'


!indfre_staggered_cmf: such that 0.5(fre(ind+1)+fre(ind)) > lam > 0.5(fre(ind)+fre(ind-1)) 
!indfre_cmf: such that fre(ind) > energy(lam) > fre(ind-1), for interpolation 
istart=2
do ii=nftot,1,-1
   lam=lam_cmf(ii) 
   energy=1.d8/lam
   freloop: do i=istart,ifre
      frec=fre(i)
      if(energy.le.frec) then
         ind=i
         diffe =abs(energy-frec)
         diffe1=abs(energy-fre(i-1))
         if (diffe.le.diffe1) then
           inds=i
         else
           inds=i-1
         endif  
         istart=i
         exit freloop
      endif
   enddo freloop
!print*,ii,i,energy,frec,fre(i-1),ind,istart 
indfre_staggered_cmf(ii)=inds
indfre_cmf(ii)=ind
enddo


! provide info on individual line-components for selected bg-elements
met_rbb%no_comp=0
met_rbb%sum_gf=0.
max_no_comp=0

!NOTE: all data created here refer to the original files (id,xlam,gf)
!1st run do define no. of unpacked components and to calculate sum(gf)
loop1: do i=1,nftot

if (i.eq.1) then
dlp=0.5d0*(lam_cmf(i+1)-lam_cmf(i))
lmin=lammin 
lmax=lam_cmf(i)+dlp  
! look for first line in interval
do ic=1,ntotnl3
  if(xlam(ic).ge.lmin) exit
enddo

else if (i.eq.nftot) then
dlm=dlp
lmin=lmax 
lmax=lam_cmf(i)  

else
dlm=dlp
lmin=lmax
dlp=0.5d0*(lam_cmf(i+1)-lam_cmf(i))
lmax=lam_cmf(i)+dlp  
endif

5 continue

! might be programmed faster (checking only the upper limit), &
! but note that certain intervals can be empty
if(xlam(ic).ge.lmin .and. xlam(ic).lt.lmax) then
  in=id_rbb(ic)
  if(in.ne.0) then
    met_rbb(in)%no_comp=met_rbb(in)%no_comp+1  
    max_no_comp=max0(max_no_comp,met_rbb(in)%no_comp)
    met_rbb(in)%sum_gf=met_rbb(in)%sum_gf+gf(ic)  
  endif  
  ic=ic+1
! out of range
  if (xlam(ic).gt.lammax) exit loop1 
  goto 5
endif

enddo loop1
print*,'max. number of unpacked components of NLTE bg-elements = ',max_no_comp

!test whether sum gf_i (unpacked) = gf (packed) ...
!note1: this test can be wrong at the borders, since unpacked levels are
!       only considered between WAVBLUE and WAVRED. If a packed level
!       is close to the border, it might contain contributions from inside
!       and outside this interval, where the latter are not considered.
!note2: With the present database, there are some cases where the
!       gf-values are considerably different, or, even worse, there is
!       no transition in the long list where there should be one according
!       to the short list (packed). This has been improved with the new
!       database (checked together with Tadziu), so that at some point
!       the new database needs to be implemented!
!
! ... and define indices where info on unpacked transitions will be located

met_rbb%quality=-1
met_rbb%index_comp=-1

iqual0=0
iqual1=0
iqual2=0
tot_comp=0

irun=1

do k=1,natom
 if (jatom_full(k).eq.0) cycle
 do j=1,jmax(k)
 n=indrec(k,j)
 jj=irbb(n)-1

   do ii=1,ilin(n)
     ji=jj+ii
     if(met_rbb(ji)%wave.gt.wavblue .and. met_rbb(ji)%wave.lt.wavred) then
     tot_comp=tot_comp+met_rbb(ji)%no_comp
! this is the index where info on the unpacked transitions will be located
     met_rbb(ji)%index_comp=irun
! includes case with no_comp=0
     irun=irun+met_rbb(ji)%no_comp
       if (abs(1.-met_rbb(ji)%sum_gf/met_rbb(ji)%gf).gt.0.2) then
!         print*,k,j,ji,met_rbb(ji)%wave,met_rbb(ji)%no_comp,met_rbb(ji)%sum_gf,met_rbb(ji)%gf
         if(met_rbb(ji)%sum_gf .eq. 0.) then
           met_rbb(ji)%quality=2      
           iqual2=iqual2+1
         else
           met_rbb(ji)%quality=1  
           iqual1=iqual1+1
        endif       
       else
         met_rbb(ji)%quality=0  
         iqual0=iqual0+1
       endif
     endif
   enddo

 enddo
enddo 
print*,'total number of depacked components of NLTE bg-elements = ',tot_comp

print*
print*,' no. of packed transitions of selected elements (within limits, but all ions)'
print*,' with quality = 0 (OK):',iqual0
print*,' with quality = 1 (large differences in Sum gf_i and gf):',iqual1
print*,' with quality = 2 (no transition in unpacked data):',iqual2
print*

!iqual=2 implies no_comp=0, sumgf=0.
!iqual=1,0 implies no_comp ge 1, sumgf gt 0.
!not considered: iqual=-1, no_comp=0, sumgf=0.

!these are the files where info on unpacked transitions will be stored
!unpacked_index_cmf: index of corresponding cmf_freq. (roughly xlam) 
!unpacked_gf: gf-value of individual component 
allocate(unpacked_index_cmf(tot_comp),unpacked_gf(tot_comp),unpacked_xlam(tot_comp))
unpacked_index_cmf=-1

!2nd run do define cmf-index and gf values for individual components
loop2: do i=1,nftot

if (i.eq.1) then
dlp=0.5d0*(lam_cmf(i+1)-lam_cmf(i))
lmin=lammin 
lmax=lam_cmf(i)+dlp  
! look for first line in interval
do ic=1,ntotnl3
  if(xlam(ic).ge.lmin) exit
enddo

else if (i.eq.nftot) then
dlm=dlp
lmin=lmax 
lmax=lam_cmf(i)  

else
dlm=dlp
lmin=lmax
dlp=0.5d0*(lam_cmf(i+1)-lam_cmf(i))
lmax=lam_cmf(i)+dlp  
endif

6 continue

! might be programmed faster (checking only the upper limit), &
! but note that certain intervals can be empty
if(xlam(ic).ge.lmin .and. xlam(ic).lt.lmax) then
  in=id_rbb(ic)

! bug found by Jon July 2015
  if(in.ne.0) then
    if(met_rbb(in)%wave.gt.wavblue .and. met_rbb(in)%wave.lt.wavred) then
      no_comp=met_rbb(in)%no_comp
      if(no_comp.eq.0 .and. met_rbb(in)%quality.ne.2) stop ' no_comp = 0 and quality ne 2'
      index_comp=met_rbb(in)%index_comp 
! includes case that no_comp=0
      do ii=index_comp,index_comp+no_comp-1
        if(unpacked_index_cmf(ii).eq.-1) goto 7
      enddo  
      print*,index_comp,no_comp,ii,met_rbb(in)%wave
      stop ' empty index for unpacked component not found'

7     unpacked_index_cmf(ii)=i  ! index of cmf_frequency
      unpacked_gf(ii)=gf(ic)    ! gf value of individual component
      unpacked_xlam(ii)=xlam(ic)! wavelength of individual component
    endif
  endif  
  ic=ic+1
! out of range
  if (xlam(ic).gt.lammax) exit loop2 
  goto 6
endif

enddo loop2

!consistency tests, can be commented out later
do k=1,natom
 if (jatom_full(k).eq.0) cycle
 do j=1,jmax(k)
 n=indrec(k,j)
 jj=irbb(n)-1

   do ii=1,ilin(n)
     ji=jj+ii
     if(met_rbb(ji)%quality.ge.0 .and. &
&      (met_rbb(ji)%wave.le.wavblue .or. met_rbb(ji)%wave.ge.wavred)) &
&    stop ' quality ge 0 outside interval'        
     if(met_rbb(ji)%wave.gt.wavblue .and. met_rbb(ji)%wave.lt.wavred) then
       if(met_rbb(ji)%quality.eq.-1) stop ' quality = -1 inside interval'        

       no_comp=met_rbb(ji)%no_comp
       index_comp=met_rbb(ji)%index_comp 
! includes case that no_comp=0
       sum_gf=0.
       do in=index_comp,index_comp+no_comp-1
         sum_gf=sum_gf+unpacked_gf(in)
       enddo
       if (sum_gf .ne. met_rbb(ji)%sum_gf) then
         print*,k,j,ji,met_rbb(ji)%wave,met_rbb(ji)%no_comp,met_rbb(ji)%sum_gf,sum_gf
       endif
     endif
   enddo

 enddo
enddo

!JO shifted out of 'first' block, since lwion1 can change
!nsize=size(met_rbb%quality)
!allocate(tlu_met(nsize,lwion1-1),tul_met(nsize,lwion1-1))

!now we have acquired all necessary info. The rest needs to be calculated
!as long as calcion varies (since xlam1 varies) 

first=.false.
endif
!--------------------------------------------------------------------
!
! calculate indices for depacked levels; this has do be done out of 'first'
! block, since depacked%index changes as long as line_list will be called
!
depacked%index=0 !reset (data from subr. line_list no longer needed)
cln=log(1.d0+vref/clight)
lamax=lam_cmf(nftot)
do i=1,lineno
  wave=depacked(i)%wave
  if(wave.ge.lammin .and. wave.le.lamax) then
    ii=log(depacked(i)%wave/lammin)/cln+1
    dlp=wave-lam_cmf(ii)
    if(dlp.lt.0.) stop 'dlp < 0'
    dlm=lam_cmf(ii+1)-wave
    if(dlm.lt.0.) stop 'dlm < 0'
    if(dlp.lt.dlm) then
      depacked(i)%index=ii
    else  
      depacked(i)%index=ii+1
    endif
!    print*,wave,ii,depacked(i)%index,lam_cmf(ii),lam_cmf(ii+1)
  endif 
enddo

!JO changed Nov. 2015
nsize=size(met_rbb%quality)
!lwion1 can change
if (.not.allocated(tlu_met)) then
  allocate(tlu_met(nsize,lwion1-1),tul_met(nsize,lwion1-1))
else
  deallocate(tlu_met,tul_met)
  allocate(tlu_met(nsize,lwion1-1),tul_met(nsize,lwion1-1))
endif
  
if(nftot.eq.0) stop ' nftot = 0 after first cmf_all iteration'

index_lam_cmf=0

loop: do i=1,nftot

if (i.eq.1) then
dlp=0.5d0*(lam_cmf(i+1)-lam_cmf(i))
lmin=lammin 
lmax=lam_cmf(i)+dlp  
! look for first line in interval
do ic=1,ntot
!  if(xlam1(ic).ge.lmin) exit
! changed by JO March 2015: wavblue is the actual begin
  if(xlam1(ic).ge.wavblue) exit
enddo
icfirst=ic
lastindex=icfirst-1
! icfirst is first valid index

else if (i.eq.nftot) then
dlm=dlp
lmin=lmax 
lmax=lam_cmf(i)  

else
dlm=dlp
lmin=lmax
dlp=0.5d0*(lam_cmf(i+1)-lam_cmf(i))
lmax=lam_cmf(i)+dlp  
endif

10 continue

! might be programmed faster (checking only the upper limit), &
! but note that certain intervals can be empty
if(xlam1(ic).ge.lmin .and. xlam1(ic).lt.lmax) then
  index_lam_cmf(i)=ic
  ic=ic+1
! out of range
!  if (xlam1(ic).gt.lammax) exit loop 
! changed by JO March 2015: wavred is the actual end
  if (xlam1(ic).gt.wavred) exit loop 
  goto 10
! no else required; in case there is no line. index_lam_cmf = 0
endif

! consistency check after sub-interval finished;
! might be skipped when everything is working
in=index_lam_cmf(i)
if (in.ne.0) then
  if(xlam1(in).gt.lmax .or. xlam1(in) .lt. lmin) stop ' error1 in fgrid_cmf'
  if(xlam1(lastindex+1).gt.lmax .or. xlam1(lastindex+1) .lt. lmin) stop ' error2 in fgrid_cmf'
  lastindex=in
endif

enddo loop

iclast=ic-1 !is the last valid index

! final consistency check
if (iclast.gt.ntot) stop ' ic-1 > ntot in fgrid_cmf'

! some statistics, might be skipped when everything is working
ic0=0
din=0
dinmax=0
lastindex=icfirst-1

do i=1,nftot

  in=index_lam_cmf(i)
  if (in.eq.0) then
    ic0=ic0+1
  else
    din1=in-lastindex !actually, in - (lastindex+1) + 1
    din=din+din1
    dinmax=max(din1,dinmax)
    lastindex=in
  endif
enddo

if(iclast-icfirst+1 .ne. din) stop ' something wrong with index finder'

print*,'considered range for complete coupling: ',lammin,lammax
print*,'number of subintervals                = ',nftot
print*,'corresponding to no. of single lines  = ',nftot/21 +1
print*,'total number of lines in subintervals = ',din,' (gf_min = ',gf_min,')'
print*,'number of subintervals with no lines  = ',ic0
print*,'max. number of lines in subintervals  = ',dinmax

return

! this is an example to find individual lines
lastindex=icfirst-1
do i=1,nftot
  in=index_lam_cmf(i)
  if (lam_cmf(i).gt.4633 .and. lam_cmf(i).lt.4645) then
    print*,lam_cmf(i),in
    if(in.ne.0) then
      do j=lastindex+1,in
        print*,j,xlam1(j),id1(j)
      enddo  
    endif
  endif
  if(in.ne.0) lastindex=in
enddo
stop
end
!
!-----------------------------------------------------------------------
!
subroutine profile_cmf(nd,temp)
!
! calculates Doppler profiles (used within subr. opacity_cmf)
! for all required elements as a a function of depth
!
! IMPORTANT NOTE
!  so far, only Doppler profiles. At least for H/He, Stark broadening needs
!  to be included to obtain realistic fluxes at Lyman series limit
!
!JO: presumably, pweightdop is not needed in final version
!test can be skipped after experience with different models has accumulated

USE nlte_type      
USE nlte_dim,only: id_ndept, id_atoms
USE fund_const, only: clight,wpi=>sqpi
USE princesa_var, only: nat, labat
USE nlte_var, only: vmax, vturb, vther
USE nlte_app, only: teff, natom, names1, indexel, vth
USE cmf_all_var, only: indexel_inv, nftot, lam_cmf, vturb_cmf, vref, &
&                vth_cmf, elements_cmf, nf_cmf, xmaxdop, &
&                profdop, profdop_H, profdop_He, &
&                pweightdop, pweightdop_H, pweightdop_He, &
&                sum_pweightdop, sum_pweightdop_H, sum_pweightdop_He, &
&                lamtest, weighttest, vdoptot

implicit none

integer(i4b), parameter :: nd1=id_ndept  

!     ..
!     .. scalar arguments ..
integer(i4b), intent(in) ::  nd  
!     ..
!     .. array arguments ..
real(dp), dimension(nd1), intent(in) :: temp

integer(i4b) k,kk,i,j,l,ixmax,ixmax_H,ixmax_He,maxixmax,kfirst,klast

integer(i4b) :: nsafety=2 !add. freq. points for wider profile (interior)

real(dp) :: x,eps,err,vdop,xmax,vdopl,del,del1,xk,ws,dx,wnue1,wnue2,sump,maxerr

logical first
logical, parameter :: test=.true.

data first/.true./

if (first) then
! index_array providing index of explict atom w.r.t. background elements
! now already performed in subr. line_list
!if (id_atoms .ne. nat) stop ' id_atoms ne nat (subr. profile_cmf)'
!do k = 1,nat
!  do i = 1,natom
!    if(labat(k).eq.names1(i)) then
!      indexel_inv(k)=i
!    endif  
!  enddo
!test
!  if (indexel(indexel_inv(k)).ne.k) stop ' something wrong with index_el(inv), H or He not explicit? (subr. profile_cmf)'
!enddo

!define vth_cmf  
vth_cmf=vth
do  k = 1,nat
  i=indexel_inv(k)
  vth_cmf(i)=vther(k) !overwrite explicit elements
enddo  

!recalculate in case vturb ne vturb_cmf (usually, vturb is lower then)
if (vturb.ne.vturb_cmf) then
 do k=1,natom
    x=sqrt(vth_cmf(k)**2-vturb**2+vturb_cmf**2)
    vth_cmf(k)=x  
 enddo
endif

!test to ensure that grid is equidistant (almost, to order vshift/c)
!with respect to !x=(nu-nu0)/delta_nu_dop = c/vth*(1-lambda0/lambda), &
!where vshift = n*vref
!'exact' error (2nd order): -n*(n+1)/2*vref/c

! test at red end, to also account for precision
j=nftot-11
i=j+10
x=clight/vref*(1.d0-lam_cmf(j)/lam_cmf(i))
eps=(x-dble(i-j))
err=-(i-j)*(i-j+1)/2.*vref/clight
if (abs(1.-eps/err).gt.0.01) stop ' problems with equidistance of x (subr. profile_cmf)'
!
! define maximum number of frequencies, and allocate profile arrays
! (H and He separately, due to wider grid)
! 

!NOTE: we need to calculate only half of profile, due to symmetry!
nf_cmf=elements_cmf

ixmax_H=0
ixmax_He=0
! at first, profiles for bg elements
do k=1,natom
! only for bg-elements which are treated in the cmf
  if(elements_cmf(k).eq.0) cycle
  vdop = vth_cmf(k)
  xmax=xmaxdop*vdop/vref
  ixmax=int(xmax)+1+nsafety
  if(indexel(k).gt.0) stop ' something wrong with indexel-1 (subr. profile_cmf)'
  if (k.eq.1) then
    ixmax_H=ixmax
  else if (k.eq.2) then
    ixmax_He=ixmax
  endif
! for ALL bg elements (incl. H/He, if present)
  nf_cmf(k)=ixmax
enddo

! now, also for explicit elements
do k=1,nat
  i=indexel_inv(k)  ! i is the index w.r.t. bg elements (ordered)
! so far, elements_cmf contain only bg elements
  if(elements_cmf(i) .ne. 0) stop ' something wrong with elements_cmf (subr. profile_cmf)'
  vdop = vth_cmf(i)
  xmax=xmaxdop*vdop/vref
  ixmax=int(xmax)+1+nsafety
  if(indexel(i).le.0) stop ' something wrong with indexel-2 (subr. profile_cmf)'
  if (i.eq.1) then
    ixmax_H=ixmax
  else if (i.eq.2) then
    ixmax_He=ixmax
  endif
! for ALL explicit elements (incl. H/He, if present)
  nf_cmf(i)=ixmax
enddo

if (ixmax_H.ne.0) then
  allocate(profdop_H(0:ixmax_H,nd))
  allocate(pweightdop_H(0:ixmax_H,nd))
  allocate(sum_pweightdop_H(-ixmax_H:ixmax_H))
! for test
  allocate(lamtest(0:ixmax_H),weighttest(0:ixmax_H))
  print*,'(half) Doppler profile for H with ',ixmax_H+1,' freq. points'
endif
if (ixmax_He.ne.0) then
  allocate(profdop_He(0:ixmax_He,nd))
  allocate(pweightdop_He(0:ixmax_He,nd))
  allocate(sum_pweightdop_He(-ixmax_He:ixmax_He))
  print*,'(half) Doppler profile for He with ',ixmax_He+1,' freq. points'
endif

! until here no change during iteration
first=.false.

else ! if not first
! check elements 3 ... natom, since they may have changed

nf_cmf(3:natom)=elements_cmf(3:natom)

do k=3,natom
! for all bg-elements (without H/He) which are treated in the cmf
  if(elements_cmf(k).eq.0) cycle
  vdop = vth_cmf(k)
  xmax=xmaxdop*vdop/vref
  ixmax=int(xmax)+1+nsafety
  nf_cmf(k)=ixmax
enddo
do k=1,nat
  i=indexel_inv(k)  ! i is the index w.r.t. bg elements (ordered)
! for all explicit elements (without H/He)
  if(i.le.2) cycle
  vdop = vth_cmf(i)
  xmax=xmaxdop*vdop/vref
  ixmax=int(xmax)+1+nsafety
  nf_cmf(i)=ixmax
enddo

endif

!for tests
!print*,'elements_cmf'
!print*,elements_cmf
!print*,'nf_cmf'
!print*,nf_cmf

! from here on, changes are possible (elements/temperature)
maxixmax=0
kfirst=0
klast=0

do k=3,natom ! now, check for ALL elements except for H,He)
!JO: here was a bug (cured Dec. 2016)
!  if(elements_cmf(k).eq.0) cycle
  if(nf_cmf(k).eq.0) cycle
  if (kfirst .eq. 0) kfirst=k
  ixmax=nf_cmf(k)
  maxixmax=max0(maxixmax,ixmax)
  klast=k
enddo
if (.not.allocated (profdop)) then
allocate(profdop(0:maxixmax,nd,kfirst:klast))
allocate(pweightdop(0:maxixmax,nd,kfirst:klast))
allocate(sum_pweightdop(-maxixmax:maxixmax,kfirst:klast))
else
deallocate(profdop, pweightdop, sum_pweightdop)
allocate(profdop(0:maxixmax,nd,kfirst:klast))
allocate(pweightdop(0:maxixmax,nd,kfirst:klast))
allocate(sum_pweightdop(-maxixmax:maxixmax,kfirst:klast))
endif

print*,'(half) Doppler profiles with max. ',maxixmax+1,' freq. points'
print*,'... for elements from ',names1(kfirst),' to ',names1(klast)

! now everything is prepared, and the profiles and weights can be calculated as in subr. FGRID
elemloop: do kk=1,natom
    
ixmax=nf_cmf(kk)
if (ixmax.eq.0) cycle elemloop
vdop=vth_cmf(kk)

if (kk.eq.1) then

do l = 1,nd  
     vdopl = (vdop**2-vturb_cmf**2)*temp(l)/teff  
     vdopl=sqrt(vdopl+vturb_cmf**2)
     vdoptot(l,kk)=vdopl
     del =  vmax/vdopl
     del1 = vref/vdopl ! x defined on vref-grid
     ws = .0d0  
     do k = 0,ixmax  
! x w.r.t. vref = k, tested above
          xk = k*del1  
!          if(l.eq.43) print*,kk,k,xk  ! at l=43, T approx Teff
          pweightdop_H(k,l) = exp(- (xk**2))  
          if (k.ne.0) ws = ws + pweightdop_H(k,l)  
          profdop_H(k,l) = pweightdop_H(k,l)*del/wpi ! OPAL includes factor lam*SR/vmax
     end do
!
!    renormalization
!
     ws=2.d0*ws+1.d0 ! symmetry + zero point, exp(-0) 
     do k = 0,ixmax  
          pweightdop_H(k,l) = pweightdop_H(k,l)/ws  
     end do
end do  

else if (kk.eq.2) then

do l = 1,nd  
     vdopl = (vdop**2-vturb_cmf**2)*temp(l)/teff  
     vdopl=sqrt(vdopl+vturb_cmf**2)
     vdoptot(l,kk)=vdopl
     del =  vmax/vdopl
     del1 = vref/vdopl ! x defined on vref-grid
     ws = .0d0  
     do k = 0,ixmax  
          xk = k*del1  
!          if(l.eq.43) print*,kk,k,xk
          pweightdop_He(k,l) = exp(- (xk**2))  
          if (k.ne.0) ws = ws + pweightdop_He(k,l)  
          profdop_He(k,l) = pweightdop_He(k,l)*del/wpi ! OPAL includes factor lam*SR/vmax
     end do
!
!    renormalization
!
     ws=2.d0*ws+1.d0 ! symmetry + zero point, exp(-0) 
     do k = 0,ixmax  
          pweightdop_He(k,l) = pweightdop_He(k,l)/ws  
     end do
end do  

else

do l = 1,nd  
     vdopl = (vdop**2-vturb_cmf**2)*temp(l)/teff  
     vdopl=sqrt(vdopl+vturb_cmf**2)
     vdoptot(l,kk)=vdopl
     del =  vmax/vdopl
     del1 = vref/vdopl ! x defined on vref-grid
     ws = .0d0  
     do k = 0,ixmax  
          xk = k*del1  
!          if(l.eq.43) print*,kk,k,xk
          pweightdop(k,l,kk) = exp(- (xk**2))  
          if(k.ne.0) ws = ws + pweightdop(k,l,kk)  
          profdop(k,l,kk) = pweightdop(k,l,kk)*del/wpi ! OPAL includes factor lam*SR/vmax
     end do
!
!     renormalization
!
     ws=2.d0*ws+1.d0 ! symmetry + zero point, exp(-0) 
     do k = 0,ixmax  
          pweightdop(k,l,kk) = pweightdop(k,l,kk)/ws  
     end do
end do  

endif

enddo elemloop

! calculation of sum(pweightdop) for H, He and metals at l=1
! (for optically thick outer boundary), allowing also for a test
ws=0.
if (nf_cmf(1).ne.0) then 
  do k=-nf_cmf(1),nf_cmf(1)
   ws = ws + pweightdop_H(abs(k),1)
   sum_pweightdop_H(k)=ws
  enddo
  if (abs(ws-1.d0).gt.1.d-10) stop ' error in pweightdop_H (subr. profile_cmf)'
else
  print*,'WARNING: H absent, not test of pweightdop possible' 
endif

ws=0.
if (nf_cmf(2).ne.0) then 
  do k=-nf_cmf(2),nf_cmf(2)
   ws = ws + pweightdop_He(abs(k),1)
   sum_pweightdop_He(k)=ws
  enddo
  if (abs(ws-1.d0).gt.1.d-10) stop ' error in pweightdop_He (subr. profile_cmf)'
else
  print*,'WARNING: He absent, not test of pweightdop possible' 
endif
  

do kk=kfirst,klast
ws=0.
if (nf_cmf(kk).ne.0) then 
  do k=-nf_cmf(kk),nf_cmf(kk)
   ws = ws + pweightdop(abs(k),1,kk)
   sum_pweightdop(k,kk)=ws
  enddo
  if (abs(ws-1.d0).gt.1.d-10) stop ' error in pweightdop (subr. profile_cmf)'
endif

enddo

! test whether profiles correctly normalized and sufficiently broad
if(test) then
! set up test frequency grid, from proto-typical lam0 until
!max. wavelength for H-profile (this is the broadest one)
!
ixmax_H=nf_cmf(1)
ixmax_He=nf_cmf(2)

if (ixmax_H.eq.0) stop ' ixmax_H = 0 (no hydrogen) in profile_cmf' 

lamtest(0)=1000.d8
dx=1.d0+vref/clight
! for weights, see subr. cmf_complete
wnue1=(1.d0-1.d0/dx)*clight*0.5
wnue2=(1.d0-1.d0/(dx**2))*clight*0.5

do k=1,ixmax_H
  lamtest(k)=lamtest(0)*dx**k  
enddo

do k=0,ixmax_H
  if(k.eq.0) then
     weighttest(k)=wnue1/lamtest(0)
  else if (k.eq.ixmax_H) then
     weighttest(k)=wnue1/lamtest(ixmax_H-1)
  else
     weighttest(k)=wnue2/lamtest(k-1)
  endif
enddo

dx=2.*lamtest(0)/vmax

print*
testloop: do kk=1,natom
  ixmax=nf_cmf(kk)
  if (ixmax.eq.0) cycle testloop

  if(kk.eq.1) then
    if (ixmax.ne.ixmax_H) stop ' ixmax ne ixmax_H in profile_cmf'
    maxerr=0.
    do l=1,nd
    sump=0.
    do k=0,ixmax_H
      sump=sump+profdop_H(k,l)*weighttest(k)
    enddo
    sump=sump*dx
    err=abs(1.-sump)
    maxerr=amax1(err,maxerr)
    if(err.gt.0.01) then
      print*,'H',l,err
      stop ' error in profile function (normalization)'
      endif
    enddo 
    print*,' max. err in normalization of profile funct. (H)   = ',maxerr   
    maxerr=0.
    
  else if(kk.eq.2) then
    if (ixmax.ne.ixmax_He) stop ' ixmax ne ixmax_He in profile_cmf'
    do l=1,nd
    sump=0.
    do k=0,ixmax_He
      sump=sump+profdop_He(k,l)*weighttest(k)
    enddo
    sump=sump*dx
    err=abs(1.-sump)
    maxerr=amax1(err,maxerr)
    if(err.gt.0.01) then
      print*,'He',l,err
      stop ' error in profile function (normalization)'
      endif
    enddo
    print*,' max. err in normalization of profile funct. (He)  = ',maxerr   
    maxerr=0.

  else
    do l=1,nd
    sump=0.
    do k=0,ixmax
      sump=sump+profdop(k,l,kk)*weighttest(k)
    enddo
    sump=sump*dx
    err=abs(1.-sump)
    maxerr=amax1(err,maxerr)
!    if(err.gt.0.001) then
! Jo changed Nov. 2016
    if(err.gt.0.01) then
      print*,kk,l,err
      stop ' error in profile function (normalization)'
      endif
    enddo

  endif
enddo testloop
print*,' max. err in normalization of profile funct. (met) = ',maxerr   

endif !test

return
end
!
!-----------------------------------------------------------------------
!
subroutine opacity_cmf(nd,xne,te,clfac,r,velo,dvdr,test_overlap,ntest,ipath)
!
! calculates summed line opacities and source-functions for all
! frequency grid points.
! Note that the transition frequencies are slightly shifted (see fgrid_cmf),
! to allow for a correct transport at peak opacity
!
! only line quantities, since cmf-transport requires separate continuum
!
! NOTE: I_core and f_edd determined from approx. treatment (CONT1), since
!       average quantities required. Since f_edd for L < Lthin is approximated,
!       this can lead to small disturbances at the end of convergence if
!       LTHIN changes. Since lines affected, this can give rise to moderate
!       changes in XJ_CMF (for a typical case, we found changes by 6%).
!       To avoid this problem, we now calculate fedd always from XXK/XJ

!
! included:
!   (i) correction of source-function for overlapping lines if up=0
!  (ii) calculation of all depacked components according to LINES_xxx.dat
! (iiI) use Ng-extrapolated source-function if ng=.true. 
  
USE nlte_type      
USE nlte_dim, only: id_ndept, id_llevs, id_nttrd
USE nlte_opt, only: opt_new_upper
USE fund_const, only: clight,hc2
USE princesa_var, only: gl, indat1, data
USE nlte_var, only: sr,vmax,blevel,fre,xj,xxk,lthin,lwion1, &
                 no_check,xne_ratio,metals_converged   
USE nlte_app, only: indrec, occngold, occng, occnglte, &
&                jatom_full, met_imin, met_imax, lwion, met_rbb, highest_level, &
                 ergion, metall
USE cmf_all_var, only: vref, lam_cmf, index_lam_cmf, indfre_staggered_cmf, &
&                nftot, icfirst, iclast, &
&                id1, xlam1, gf1, indexel_inv, nf_cmf, &
&                occexpl, ntot_expl_inside, &
&                profdop, profdop_H, profdop_He, & 
&                sum_pweightdop, sum_pweightdop_H, sum_pweightdop_He, & 
&                sumopal_cmf, sumetal_cmf, &
&                sumopal_outer, sline_outer, &
&                sumtaus_outer, &
&                indexcmf1,indexcmf1_lam,opal_ntest,id1_overlap, &
&                xm,minrho,no2,vdoptot, &
&                indxlamc_depacked,xlamcmf_depacked
USE ng_var, only: ng, n_ng, trans_ng, sl_ng, index_ng, &
                  ng_met, n_ng_met, trans_ng_met, sl_ng_met, index_ng_met

implicit none

integer(i4b), parameter :: nrec=id_llevs
integer(i4b), parameter :: nd1=id_ndept

real(dp), parameter:: c_vanreg= 0.5 * 2.055d-23
real(dp), parameter:: c_forbid= 8.63d-6/(0.5*1.3707d-7)
real(dp), parameter :: gfcut=1.d-5 ! see sumopal

real(dp), parameter :: c1=1.4388354967334d8,c2=3.972970127d8
! for Planck function, identical with function bnue


integer(i4b), intent(in) :: nd, ipath

real(dp), dimension(nd1), intent(in) :: xne,te,clfac,r,velo,dvdr

logical :: test_overlap

integer(i4b), dimension(nd1) :: lll

logical, dimension(nd1) :: check

real(dp), dimension(nd1) :: opal, sline, occlow, occup, &
&                           vr, tau0, xxj, mustar, eddf, &
&                           collfac, u0, vanreg, vanreg1, eps, eps1, &
&                           opalexp, slineexp

real(dp), parameter :: c_tau=1.d-8*0.02654, hc2_8=hc2*1.d24 !conversion to Angstrom

real(dp) :: srvmax, c_tau1, const, lam, lamgrid, test_dl, xlamgrid, opal1, etal1, &
&           occlowi, occupi, taubi, taui, betacic, beta, bni, epsi, eps1i, lam3, &
&           xxlam, precis, rel, lamexp, constexp, mirho, rho1, deltav,vdopmax, sl, &
&           xlamc, fac, erg_fake, occupi_fake

integer(i4b) :: i,j,lastindex,in,ii,index,k,kk,ll,irest,ml,mu,low,up,ixmax, &
&               irec,icount_expl,icount_bg, icount_bg_noup, &
&               inds_old,inds,lth,ik,ic,ntest, &
&               iiexp,indexexp,kkexp,kexp,iiexp_old,icount_ov,ix,ifexp,i_ng, &
&               indi,up_fake

!print*,'ipath=',ipath

if (.not.allocated(sumopal_cmf)) then
! only once
 allocate(sumopal_cmf(nd1,nftot),sumetal_cmf(nd1,nftot))
 allocate(sumopal_outer(nftot), sline_outer(nftot)) 
 allocate(sumtaus_outer(nftot)) 
endif
! initialize
if(ipath.eq.1) sumopal_cmf=0.  
sumetal_cmf=0.  
sumopal_outer=0.
sline_outer=0.
sumtaus_outer=0. 

if(ipath.eq.1) id1_overlap=0

srvmax=sr/vmax
c_tau1=c_tau*srvmax
test_dl=0.5*vref/clight

vr=velo/r
tau0=1.d0/(dvdr/3.+2./3.*vr)  ! at mue^2 = 1/3
collfac=c_vanreg*xne/sqrt(te)

check=.false.
!calculate check-variable (only for lwion (<lwion1) 
if(lwion1.lt.lwion) stop ' something wrong with lwion1'
do ll=lwion,nd
  rel=abs(1.d0-xne_ratio(ll))
  check(ll)=.not.no_check .and. rel.lt.0.01
!  print*,'check',ll,check(ll),rel
enddo

! read occupation numbers at all depth points
do ll=1,nd
  read (17,rec=ll) (occexpl(i,ll),i=1,nrec)  
  lll(ll)=ll !auxiliary vector
enddo

icount_expl=0 !counter for lines from explicit elements
icount_bg =0 !counter for lines from background elements
icount_bg_noup =0 !counter for lines from explicit elements without upper level
icount_ov =0!counter for lines with up=0 and modified source-function

! loop over all cmf-frequencies, to calculate opal and etal
lastindex=icfirst-1
inds_old=0 ! staggered

ic=0
iiexp_old=0

do i=1,nftot
  in=index_lam_cmf(i)
  lamgrid=lam_cmf(i)
  inds=indfre_staggered_cmf(i) ! for calculation on staggered grid
  
! lines shifted to lam to simplify profile calculations (uniform x)
!  and to obtain 'correct' line centers
  if(in.ne.0) then
      overlap: do ii=lastindex+1,in
        index=id1(ii)
        if (index.ge.1d8) then
        if(ipath.eq.2.and.id1_overlap(ii).ne.0) stop ' error(1) id1_overlap'
! explicit element
        kk = index/1d8
        irest = index - kk*1d8
        ml = irest/10000
        mu = irest - ml*10000
        k=indexel_inv(kk)
! in units of vref = vturb/3 (with vturb=max(vturb,vturbmin)
        ixmax=nf_cmf(k)
        lam=xlam1(ii)
        lam3=lam**3
! nevertheless, for consistency actual wavelengths used for front-factors
!        if(xlam1(ii).gt.6000 .and. xlam1(ii).lt.7000) print*,index,lam,k,ml,mu,ixmax
! for tests
        if (abs(1.-lamgrid/lam).gt. test_dl) stop ' error1 in freq. shift (subr. opacity_cmf)'
        const=c_tau1*lam*gf1(ii)
        occlow=occexpl(ml,:)/gl(ml)
         occup=occexpl(mu,:)/gl(mu)
        if(occlow(1)*occup(1).eq.0.) then
! this should never happen (at least for GLOBAL = .true.)
!JO: check for GLOBAL = .false.
          print*,index,' ',lam
          stop ' occlow or occup eq 0 (expl. elements) in opacity_cmf'
        endif
! pi e2/me c * 1.d-8 *SR/vmax * lam * gf * (nl/gl - nu/gu) /clfac
        opal=const*(occlow-occup)/clfac
        sline=hc2_8/lam3/(occlow/occup-1.d0)
        if(ng) then
! no check for almost_converged, since controlled in prep_ng 
            do i_ng=1,n_ng
              if(ml.eq.trans_ng(i_ng,1).and.mu.eq.trans_ng(i_ng,2)) then
                sl=sl_ng(35,i_ng,4)
                if(sl.gt.0.d0) then
! rescale to actual frequency (extrapolated source function calc. with xlamc)
                  indi = indat1(index_ng(i_ng))  
                  xlamc = data(indi+1)  
                  fac=xlamc**3/lam3
!                  print*,'sl_ng:',lam,xlamc,i_ng,sline(35),sl*fac
                  sline=sl_ng(:,i_ng,4)*fac
                endif
                exit
              endif  
            enddo
        endif  
        
        icount_expl=icount_expl+1
!
!       save the freq. index for rateeq
        if (ic.eq.0) then
          do ic=1,id_nttrd
! first line to be considered
            if (indexcmf1(indxlamc_depacked(ic)).ge.3) exit
          enddo
        endif  
        ik=indxlamc_depacked(ic)

!JO Dec. 2016: account for Sobo-lines within range
        do while (indexcmf1(ik).eq.1)
          ic=ic+1
          ik=indxlamc_depacked(ic)
        enddo  
        xxlam=xlamcmf_depacked(ik)
        
        if (abs(1.-lam/xxlam).lt.1.d-6) then !lam in single prec.
          if(indexcmf1(ik).lt.3) stop ' indexcmf1(ik) lt 3'
          indexcmf1_lam(ik)=i 
!          print*,i,ic,ik,lam,xxlam,indexcmf1(ik)
!          print*
          ic=ic+1
! in case, use extrapolated source-function
        endif
        
        else
! background element
        k = index/1000000
        irest = index - k*1000000
        j = irest/100000
        irest = irest - j*100000
        low = irest/100
        up  = irest-low*100
        ixmax=nf_cmf(k)
! nevertheless, for consistency actual wavelengths used for front-factors
        lam=xlam1(ii)
        lam3=lam**3
!        if(xlam1(ii).gt.6000 .and. xlam1(ii).lt.7000) print*,index,lam,k,j,low,up,ixmax
        if (abs(1.-lamgrid/lam).gt. test_dl) stop ' error2 in freq. shift (subr. opacity_cmf)'
        irec=indrec(k,j)
        const=c_tau1*lam*gf1(ii)
        if(up.ne.0) then
! for tests: neglect opacity for corresponding OIII lines
!        if(k.eq.8.and.j.eq.3.and.low.eq.1) then
!          if(up.eq.5 .or. up.eq.6 .or. up.eq.8) const=c_tau1*lam*1.d-10
!        endif  
        if(ipath.eq.2.and.id1_overlap(ii).ne.0) stop ' error(2) id1_overlap'
!
! lower and upper level  present
!
        if(jatom_full(k).eq.1) then
! selected elements, division lwion1
        do  ll=1,lwion1-1
!outside
          occlowi=occngold(irec,low,ll) ! occng = n/g
          if(j.lt.met_imin(k,ll).or.j.gt.met_imax(k,ll).or.occlowi.eq.0.) then
            opal(ll)=0.
            sline(ll)=0.
          else
            occupi=occngold(irec,up,ll)
            opal(ll)=const*(occlowi-occupi)/clfac(ll)
            sline(ll)=hc2_8/lam3/(occlowi/occupi-1.d0)
          endif
        enddo 

        if(ng_met) then
! should never happen
          if(metals_converged) stop ' metals_converged and ng_met in subr. opacity_cmf'
            do i_ng=1,n_ng_met
              if(index.eq.trans_ng_met(i_ng)) then
                sl=sl_ng_met(35,i_ng,4)
                if(sl.gt.0.d0) then
! rescale to actual frequency (extrapolated source function calc. with xlamc)
                  indi = index_ng_met(i_ng)  
                  xlamc = met_rbb(indi)%wave  
                  fac=xlamc**3/lam3
!                  print*,'sl_ng_met:',lam,xlamc,i_ng,sline(35),sl*fac
                  sline(1:lwion1-1)=sl_ng_met(:,i_ng,4)*fac
                endif
                exit
              endif  
            enddo
        endif  

        do ll=lwion1,nd 
!inside
          occlowi=occng(irec,low,ll) ! occng = n/g
          if (check(ll) .and. abs(1.-occlowi/occnglte(irec,low,ll)).gt.0.05) then
            print*,'lte-error1 ',k,j,irec,low,ll,occlowi,occnglte(irec,low,ll)
            stop ' occng(selected, lower atmosphere) ne occnglte'
          endif
          occupi=occng(irec,up,ll)
          if(occlowi*occupi.eq.0) then
           print*,index,' ',lam
           stop ' occlow or occup eq 0 (inner selected elements) in opacity_cmf'
          endif
          opal(ll)=const*(occlowi-occupi)/clfac(ll)
          sline(ll)=c2/(exp(c1/lam/te(ll))-1.d0)/lam3  
        enddo

        else
! other elements (approximate), division lwion
        do  ll=1,lwion-1
           occlowi=occng(irec,low,ll) ! occng = n/g
           occupi=occng(irec,up,ll)
          if(occlowi*occupi.eq.0) then !change of ionization
           print*,index,' ',lam
           stop ' occlow or occup eq 0 (outer approx. elements) in opacity_cmf'
          endif
          opal(ll)=const*(occlowi-occupi)/clfac(ll)
          sline(ll)=hc2_8/lam3/(occlowi/occupi-1.d0)
        enddo 
        do ll=lwion,nd 
!inside
          occlowi=occng(irec,low,ll) ! occng = n/g
          if (check(ll) .and. abs(1.-occlowi/occnglte(irec,low,ll)).gt.0.05) then
            print*,'lte-error2 ',k,j,irec,low,ll,occlowi,occnglte(irec,low,ll)
            stop ' occng(approx., lower atmosphere) ne occnglte'
          endif
          occupi=occng(irec,up,ll)
          if(occlowi*occupi.eq.0) then
           print*,index,' ',lam
           stop ' occlow or occup eq 0 (inner approx.  elements) in opacity_cmf'
          endif
          opal(ll)=const*(occlowi-occupi)/clfac(ll)
          sline(ll)=c2/(exp(c1/lam/te(ll))-1.d0)/lam3  
        enddo

        endif
        icount_bg=icount_bg+1
!
! only lower level  present
!
        else
!
! prepare two-level approach
! sline approximated by two level atom; see subr. netmat_met
! note: illuminating radiation field (for betacIc -> beta J)
!       from smooth (pseudo-)continuum

        if (inds.ne.inds_old) then
! otherwise, already calculated; change approach when interpolating in freq.,
! and use indfre_cmf instead of indfre_staggered_cmf
          xxj=xj(:,inds)
!JO: April 2018; now, the calculated Eddington factor is used everywhere
!    (also outside lthin), to avoid oscillations) 
!          lth=lthin(inds)
!          where(lll.lt.lth)
!           mustar=1.d0-(r(lth)/r)**2 !correction for lthin (checked)
!           mustar=sqrt(mustar)
!           eddf=(1.d0-mustar**3)/3.d0/(1.d0-mustar)! (checked)
!          elsewhere          
           eddf=xxk(:,inds)/xxj
!          endwhere
          if(.not. no_check .and. (eddf(nd).lt.0.31 .or. eddf(nd).gt.0.35)) &
&           stop ' inconsistent Edd. factor in opacity_cmf'
          
! OLD comment: eddf (approx.) checked, smooth transition between outer and inner region,
!          goes to unity for large radii 

! prepare calculation of eps' = Cul/Aul (1-exp(-u0)) 
! (see also subroutine opacitl and sumopal)
! to save time, frequency at continuum grid
          xlamgrid=1.d8/fre(inds)
!JO (31. Oct. 2016): here was the bug
!          u0=1.4388d8/xlamgrid*te
          u0=1.4388d8/(xlamgrid*te)
          vanreg=collfac*xlamgrid**3*(1.d0-exp(-u0)) !checked, OK
! vanreg1 propto lam^2 needs to be finally divided by gf 
          vanreg1=vanreg*c_forbid/xlamgrid !checked, OK, 
          vanreg=vanreg/(1.+vanreg)
!          print*
!          print*,xlamgrid,vanreg(1),vanreg1(1)
!          print*,xlamgrid,vanreg(43),vanreg1(43)
!          print*,xlamgrid,vanreg(nd),vanreg1(nd)
          inds_old=inds
        endif
!
! line specific quantities, at all depth-points
!
        if(gf1(ii).le.gfcut) then
          eps=vanreg1/gf1(ii)
          eps=eps/(1.+eps)           
        else
          eps=vanreg
        endif
!
!check for overlap with line from explict element,
!in case modify source-function later on
        if(ipath.eq.2) then
          iiexp=id1_overlap(ii)          
          if(iiexp.ne.0 .and.iiexp.ne.iiexp_old) then
            indexexp=id1(iiexp)
            if (indexexp.lt.1d8) stop ' indexexp does not refer to explicit element'
            kkexp = indexexp/1d8
            kexp=indexel_inv(kkexp)
            if(kexp.eq.1 .or. kexp.eq.2) stop ' indexexp refers to H/He'
            irest = indexexp - kkexp*1d8
            ml = irest/10000
            mu = irest - ml*10000
! in units of vref = vturb/3 (with vturb=max(vturb,vturbmin)
            lamexp=xlam1(iiexp)
!frequency index
            ifexp=nint(log(lamexp/lam_cmf(1))/log(1.d0+vref/clight))+1
            if (abs(1.-lam_cmf(ifexp)/lamexp).gt. test_dl) &
&              stop ' error in ifexp (subr. opacity_cmf)'
!        print*
!        print*,ml,mu,lam,lamexp
            constexp=c_tau1*lamexp*gf1(iiexp)
            occlow=occexpl(ml,:)/gl(ml)
            occup=occexpl(mu,:)/gl(mu)
! pi e2/me c * 1.d-8 *SR/vmax * lam * gf * (nl/gl - nu/gu) /clfac
            opalexp=constexp*(occlow-occup)/clfac
!JO NOTE: thus far, we don't use extrapolated source-function here,
!         for reasons of comp. time. Should not play any role,
!         since only no_up lines affected. Expected to become
!         consistent in course of iteration 
            slineexp=hc2_8/lamexp**3/(occlow/occup-1.d0)
            deltav=abs(lam/lamexp-1.d0)*clight
!            print*,'new',index,indexexp,slineexp(1)
            icount_ov=icount_ov+1
            iiexp_old=iiexp
          else if(iiexp.ne.0 .and.iiexp.eq.iiexp_old) then
!            print*,'old',index,indexexp,slineexp(1)
            deltav=abs(lam/lamexp-1.d0)*clight
            icount_ov=icount_ov+1
          endif 
        endif !ipath=2
!        
        if(jatom_full(k).eq.1) then
! selected elements, division lwion1
ll_loop1: do  ll=1,lwion1-1
!outside
          occlowi=occngold(irec,low,ll) ! occng = n/g
          if(j.lt.met_imin(k,ll).or.j.gt.met_imax(k,ll).or.occlowi.eq.0.) then
            opal(ll)=0.
            sline(ll)=0.
          else
            opal(ll)=const*occlowi/clfac(ll) ! includes gstat, ind. emission negl.
            taubi=opal(ll)/(eddf(ll)*dvdr(ll)+(1.-eddf(ll))*vr(ll)) ! at mue^2 = eddf
! pi e2/me c * 1.d-8 *SR/vmax * lam * gf * nl/gl / clfac / (eddf*dvdr + (1-eddf*v/r))
        
            if(taubi.lt.1.d-5) then
            betacic=xxj(ll)
            else
            betacic=(1.-exp(-taubi))/taubi*xxj(ll)
            endif
            
            taui=opal(ll)*tau0(ll)  ! at mue^2 = 1/3
! pi e2/me c * 1.d-8 *SR/vmax * lam * gf * nl/gl / clfac / (dvdr/3. + 2./3.*v/r)
            if(taui.lt.1.d-5) then
            beta=1.d0
            else
            beta=(1.-exp(-taui))/taui
            endif
!
! now, we can calculate the two-level source-function. We inline Bnue, &
! and evaluate at actual frequency, to ensure correct thermalization
!
            epsi=eps(ll)
            eps1i=1.d0-epsi
            bni = c2/ (exp(c1/lam/te(ll))-1.d0)/lam3  

            sline(ll)=(eps1i*betacic+epsi*bni)/(epsi+beta*eps1i)

            if(opt_new_upper) then
!JO Nov. 2018 
            up_fake=highest_level(irec,low)
! if up_fake = 0, no info available, and we keep TLA approach.
! These lines are important, e.g., for hot stars around 1000 A
! (e.g., from high lying Fe-levels) 
            if(up_fake.ne.0 .and. .not.(low.eq.1 .or. metall(k,j,low).eq.'m')) then
!            if(up_fake.ne.0) then
! Otherwise, we assume that the upper level is in LTE w.r.t. up_fake
! In other words: the departures of up_fake and up are assumed to be identical 
              erg_fake=1.d8/lam-(ergion(k,j,up_fake)-ergion(k,j,low))

              if(erg_fake.lt.0.) then
! the upper level lies below up_fake
                erg_fake=-erg_fake
!JO in case, allow for a somewhat larger range 
                if(erg_fake/ergion(k,j,up_fake).gt.0.1) then
                   print*,k,j,low,lam
                   print*,erg_fake,ergion(k,j,up_fake),ergion(k,j,low)
                   stop ' upper level too much below up_fake'
                endif
! fake upper level now more populated than up_fake (since below)
                occupi_fake=occngold(irec,up_fake,ll)*exp(c1*erg_fake*1.d-8/te(ll))
              else
! typical situation: upper level less populated than up_fake
                occupi_fake=occngold(irec,up_fake,ll)*exp(-c1*erg_fake*1.d-8/te(ll))
              endif
              opal(ll)=const*(occlowi-occupi_fake)/clfac(ll)
              sline(ll)=hc2_8/lam3/(occlowi/occupi_fake-1.d0)
            endif
            endif
! end of opt_new_upper treatment
            
!              
! check (similar to subr. overlap_noup) whether strong overlap
! with line from explicit element
            if(ipath.eq.2.and.iiexp.ne.0) then
              vdopmax=xm*vdoptot(ll,kexp)
              if(deltav.gt.vdopmax) cycle ll_loop1 !outside depth-dependent range
                                                      !nothing to be done
              ix=vdopmax/vref+1 !in units of vref/3
              mirho=1000.
              do ik=0,ix
                rho1=opalexp(ll)*profdop(ik,ll,kexp)/sumopal_cmf(ll,ifexp+ik)
!                if(index.eq.26500400) then
!                  print*,ll,ifexp+ik,rho1,indexexp
!                  print*,opalexp(ll),profdop(ik,ll,kexp),sumopal_cmf(ll,ifexp+ik)
!               endif
!  significant or dominating inversion of bg opacities
                if(rho1.gt.1.1d0 .or.  rho1.lt.0.d0) rho1=0.d0
                mirho=amin1(mirho,rho1)
              enddo   
              do ik=1,ix
                rho1=opalexp(ll)*profdop(ik,ll,kexp)/sumopal_cmf(ll,ifexp-ik)
!                if(index.eq.26500400) then
!                  print*,ll,ifexp-ik,rho1,indexexp
!                  print*,opalexp(ll),profdop(ik,ll,kexp),sumopal_cmf(ll,ifexp-ik)
!                endif
!  significant or dominating inversion of bg opacities
                if(rho1.gt.1.1d0 .or.  rho1.lt.0.d0) rho1=0.d0
                mirho=amin1(mirho,rho1)
              enddo   
! in case, modify source-function (approximating beta_tot by beta etc)
              if(mirho.gt.minrho) then
!                 print*,ml,mu,ll,lamexp,mirho
                  sline(ll)=eps1i*((1.-beta)*slineexp(ll)+betacic)+epsi*bni
              endif
            endif !ipath=2 
          endif !imin,imax
        enddo ll_loop1
        do ll=lwion1,nd 
!inside
          occlowi=occng(irec,low,ll) ! occng = n/g
          if (check(ll) .and. abs(1.-occlowi/occnglte(irec,low,ll)).gt.0.05) then
            print*,'lte-error3 ',k,j,irec,low,ll,occlowi,occnglte(irec,low,ll)
            stop ' occng(selected, lower atmosphere) ne occnglte'
          endif
          if(occlowi.eq.0) then
           print*,index,' ',lam
           stop ' occlow eq 0 (inner selected elements) in opacity_cmf'
          endif
          opal(ll)=const*occlowi/clfac(ll)
          sline(ll)=c2/(exp(c1/lam/te(ll))-1.d0)/lam3
        enddo

        else
! other elements (approximate), division lwion
ll_loop2:   do  ll=1,lwion-1
            occlowi=occng(irec,low,ll) ! occng = n/g
            if(occlowi.eq.0) then
             print*,index,' ',lam
             stop ' occlow eq 0 (outer approx. elements) in opacity_cmf'
            endif
            opal(ll)=const*occlowi/clfac(ll) ! includes gstat, ind. emission negl.
            taubi=opal(ll)/(eddf(ll)*dvdr(ll)+(1.-eddf(ll))*vr(ll)) ! at mue^2 = eddf
! pi e2/me c * 1.d-8 *SR/vmax * lam * gf * nl/gl / clfac / (eddf*dvdr + (1-eddf*v/r))
        
            if(taubi.lt.1.d-5) then
            betacic=xxj(ll)
            else
            betacic=(1.-exp(-taubi))/taubi*xxj(ll)
            endif
 
            taui=opal(ll)*tau0(ll)  ! at mue^2 = 1/3
! pi e2/me c * 1.d-8 *SR/vmax * lam * gf * nl/gl / clfac / (dvdr/3. + 2./3.*v/r)
            if(taui.lt.1.d-5) then
            beta=1.d0
            else
            beta=(1.-exp(-taui))/taui
            endif
!
! now, we can calculate the two-level source-function. We inline Bnue, &
! and evaluate at actual frequency, to ensure correct thermalization
!
            epsi=eps(ll)
            eps1i=1.d0-epsi
            bni = c2/ (exp(c1/lam/te(ll))-1.d0)/lam3  
            sline(ll)=(eps1i*betacic+epsi*bni)/(epsi+beta*eps1i)
!
! check (similar to subr. overlap_noup) whether strong overlap
! with line from explicit element
            if(ipath.eq.2.and.iiexp.ne.0) then
              vdopmax=xm*vdoptot(ll,kexp)
              if(deltav.gt.vdopmax) cycle ll_loop2 !outside depth-dependent range
                                                      !nothing to be done
              ix=vdopmax/vref+1 !in units of vref/3
              mirho=1000.
              do ik=0,ix
                rho1=opalexp(ll)*profdop(ik,ll,kexp)/sumopal_cmf(ll,ifexp+ik)
!  significant or dominating inversion of bg opacities
                if(rho1.gt.1.1d0 .or.  rho1.lt.0.d0) rho1=0.d0
                mirho=amin1(mirho,rho1)
              enddo   
              do ik=1,ix
                rho1=opalexp(ll)*profdop(ik,ll,kexp)/sumopal_cmf(ll,ifexp-ik)
!  significant or dominating inversion of bg opacities
                if(rho1.gt.1.1d0 .or.  rho1.lt.0.d0) rho1=0.d0
                mirho=amin1(mirho,rho1)
              enddo   
! in case, modify source-function (approximating beta_tot by beta etc)
              if(mirho.gt.minrho) then
!                print*,ml,mu,ll,lamexp,mirho
                sline(ll)=eps1i*((1.-beta)*slineexp(ll)+betacic)+epsi*bni
              endif
            endif !ipath=2 
!        
        enddo ll_loop2
        do ll=lwion,nd 
!inside
          occlowi=occng(irec,low,ll) ! occng = n/g
          if (check(ll) .and. abs(1.-occlowi/occnglte(irec,low,ll)).gt.0.05) then
            print*,'lte-error4 ',k,j,irec,low,ll,occlowi,occnglte(irec,low,ll)
            stop ' occng(approx., lower atmosphere) ne occnglte'
          endif
          if(occlowi.eq.0) then
           print*,index,' ',lam
           stop ' occlow eq 0 (inner approx. elements) in opacity_cmf'
          endif
          opal(ll)=const*occlowi/clfac(ll)
          sline(ll)=c2/(exp(c1/lam/te(ll))-1.d0)/lam3  
        enddo
        endif
!        
        icount_bg_noup=icount_bg_noup+1
      endif ! up eq or ne 0

      endif ! explicit or bg element

      if(test_overlap) opal_ntest(ii)=opal(ntest)

!      if(ipath.eq.2 .and. xlam1(ii).ge.387.31 .and. xlam1(ii).le.387.40) then! .and. gf1(ii).ge.0.05) then     
!        print*,id1(ii),xlam1(ii)
!        do ll=1,nd
!          print*,ll,opal(ll),sline(ll)
!        enddo  
!      endif  

!
! having calculated opacities and source functions, we can add them
! (profile-weighted) to all affected frequencies +/- ixmax 
! in parallel, we calculate auxiliary variables for outer boundary, namely
! sum(opal), sum(sline*opal) and sum(opal*sum(pweight))
! note that sum(pweight) := sum_pweight is defined from -ixmax...ixmax 
!
! for outer boundary
      opal1=opal(1)
      etal1=opal1*sline(1)

!      if(lam.gt.1545. .and. lam.lt.1555) then
!        print*,lam,k,j,opal1,sline(1)
!      endif
!for tests
!      if (k.eq.12.and.j.eq.2.and.low.eq.3.and.up.eq.5) then        
!        print*,lam
!        do ll=1,nd
!          print*,ll,opal(ll),sline(ll),met_imin(k,ll),met_imax(k,ll)
!        enddo
!      endif  

!
!sum_tausouter   is sum(opal * int(x...xmax) phi(x)dx) for iminus 
!sum_sline_outer is sum (opal*sline)/sum(opal) (without profile function)
!                for iminus 

!i-ik<0 should not happen for all elements,
!since 'safety region below wavblue and beyond wavred 
!
!calculate sumopal_cmf only in first run 
      if (k.eq.1) then
      do ik=0,ixmax
        eps=opal*profdop_H(ik,:) 
        eps1=eps*sline
        if(ipath.eq.1) sumopal_cmf(:,i+ik)=sumopal_cmf(:,i+ik)+eps
        sumetal_cmf(:,i+ik)=sumetal_cmf(:,i+ik)+eps1
! for outer boundary
        sumopal_outer(i+ik)=sumopal_outer(i+ik)+opal1
        sline_outer(i+ik)=sline_outer(i+ik)+etal1
        sumtaus_outer(i+ik)=sumtaus_outer(i+ik)+opal1*sum_pweightdop_H(ik) 
        if(ik.ne.0) then
        if(i-ik.le.0) stop ' i-ik le 0 for H'
        if(ipath.eq.1) sumopal_cmf(:,i-ik)=sumopal_cmf(:,i-ik)+eps
        sumetal_cmf(:,i-ik)=sumetal_cmf(:,i-ik)+eps1
! for outer boundary
        sumopal_outer(i-ik)=sumopal_outer(i-ik)+opal1
        sline_outer(i-ik)=sline_outer(i-ik)+etal1
        sumtaus_outer(i-ik)=sumtaus_outer(i-ik)+opal1*sum_pweightdop_H(-ik) 
        endif
      enddo
!      print*,k,lam,ixmax
        
      else if (k.eq.2) then
      do ik=0,ixmax
        eps=opal*profdop_He(ik,:) 
        eps1=eps*sline
        if(ipath.eq.1) sumopal_cmf(:,i+ik)=sumopal_cmf(:,i+ik)+eps
        sumetal_cmf(:,i+ik)=sumetal_cmf(:,i+ik)+eps1
! for outer boundary
        sumopal_outer(i+ik)=sumopal_outer(i+ik)+opal1
        sline_outer(i+ik)=sline_outer(i+ik)+etal1
        sumtaus_outer(i+ik)=sumtaus_outer(i+ik)+opal1*sum_pweightdop_He(ik) 
        if(ik.ne.0) then
        if(i-ik.le.0) stop ' i-ik le 0 for He'
        if(ipath.eq.1) sumopal_cmf(:,i-ik)=sumopal_cmf(:,i-ik)+eps
        sumetal_cmf(:,i-ik)=sumetal_cmf(:,i-ik)+eps1
! for outer boundary
        sumopal_outer(i-ik)=sumopal_outer(i-ik)+opal1
        sline_outer(i-ik)=sline_outer(i-ik)+etal1
        sumtaus_outer(i-ik)=sumtaus_outer(i-ik)+opal1*sum_pweightdop_He(-ik) 
        endif
      enddo
!      print*,k,lam,ixmax

      else  
      do ik=0,ixmax
        eps=opal*profdop(ik,:,k) 
        eps1=eps*sline
        if(ipath.eq.1) sumopal_cmf(:,i+ik)=sumopal_cmf(:,i+ik)+eps
        sumetal_cmf(:,i+ik)=sumetal_cmf(:,i+ik)+eps1
! for outer boundary
        sumopal_outer(i+ik)=sumopal_outer(i+ik)+opal1
        sline_outer(i+ik)=sline_outer(i+ik)+etal1
        sumtaus_outer(i+ik)=sumtaus_outer(i+ik)+opal1*sum_pweightdop(ik,k) 
        if(ik.ne.0) then
        if(i-ik.le.0) stop ' i-ik le 0 for metals'
        if(ipath.eq.1) sumopal_cmf(:,i-ik)=sumopal_cmf(:,i-ik)+eps
        sumetal_cmf(:,i-ik)=sumetal_cmf(:,i-ik)+eps1
! for outer boundary
        sumopal_outer(i-ik)=sumopal_outer(i-ik)+opal1
        sline_outer(i-ik)=sline_outer(i-ik)+etal1
        sumtaus_outer(i-ik)=sumtaus_outer(i-ik)+opal1*sum_pweightdop(-ik,k) 
        endif
      enddo
!
      endif !k=1,2,else

      enddo overlap! all overlapping lines at cmf freq. lamgrid
      lastindex=in
  endif ! if lines present at lamgrid
enddo ! all cmf-frequencies 


if (icount_expl.ne.ntot_expl_inside) then
  print*,icount_expl,' ',ntot_expl_inside
  stop ' not all lines from expl. elements found (subr. opacity_cmf)'
endif
  
! finally, calculate line source function for outer boundary
where (sumopal_outer.ne.0.)
  sline_outer=sline_outer/sumopal_outer
elsewhere
  sline_outer=0.
endwhere

!test whether no line(s) at first frequency (for blue wing boundary condition)
do ll=1,nd
  if(sumopal_cmf(ll,1).ne.0.) stop ' line at blue wing (subr. opacity_cmf)'
enddo  

print*
print*,icount_expl+icount_bg+icount_bg_noup,' cmf line opacities/emissivites calculated'
print*,icount_expl,' lines from explicit elements' 
print*,icount_bg,' lines from background elements with upper level' 
print*,icount_bg_noup,' lines from background elements with no upper level' 
if(ipath.eq.2) then
  print*,icount_ov,' lines with no upper level and modified source-function' 
  if(icount_ov.ne.no2) stop ' icountov ne no2'
endif  
  
!do i=1,nftot
!  write(*,100) lam_cmf(i),sumetal_cmf(1,i)/sumopal_cmf(1,i), &
!&   sumetal_cmf(43,i)/sumopal_cmf(43,i),sumetal_cmf(nd,i)/sumopal_cmf(nd,i)
!enddo

!do i=1,nftot
!  write(*,100) lam_cmf(i),sumopal_cmf(1,i),sumopal_cmf(43,i),sumopal_cmf(nd,i)
!enddo

! test whether all explicit lines have correct wavelengths 
precis=vref/clight
do j=1,id_nttrd
  ii=indxlamc_depacked(j) 
  if (indexcmf1(ii).ge.3) then
!    if(abs(1.-lam_cmf(indexcmf1_lam(ii))/xlamcmf_depacked(ii)).gt.precis) &
!&     stop ' error in indexcmf1_lam'
    if(abs(1.-lam_cmf(indexcmf1_lam(ii))/xlamcmf_depacked(ii)).gt.precis) &
&      print*,lam_cmf(indexcmf1_lam(ii)),xlamcmf_depacked(ii)
  endif  
enddo

return

100 format(f14.4,2x,6(e12.6,2x))

end
!
!----------------------------------------------------------------------------
!
subroutine cmf_complete(nd,xne,temp,clfac,r,velo,dvdr,rho,taur,xnh, &
                        pressure,ferr,dtfcorr,opt_ray)
!
! complete cmf transfer for all elements 
!
! here, r, v, and dvdr are dimensionless (normalized), 
! and the opacities/emissivities include line and continuum contributions.
! line contributions already profile-weighted (depth dependend) and summed  
!
! both line and continuum quantities multiplied by Rstar
! line quantities weighted with Rstar/ delta_nue_dop_max = Rstar*lambda/vmax
! corrected in profile propto vmax/vdop(r) * exp(-xth(r)^2)/sqrt(pi)
! 
! in these units, deltax = vref/vmax
!
! flux-like quantities (ah, ah_int, ...) from 1...nd+1 as usual
! will be interpolated to regular grid,
! with quantities xh, xh_int etc. (incl. factor r^2)
  
! radiative accelerations only on staggered grid, from 1...nd-1

USE nlte_type
USE nlte_dim, only: id_ndept, id_npoin
USE fund_const, only: amh, sigmae, clight, pi, sigsb, akb
USE nlte_opt, only: optphotlines
USE run_once, only: start_cmf_simple
USE nlte_var, only : sr, vmax, rtau23, corrfc, ifre, fre, wfre, &
&                    opac_nolines, etat_nolines, xj, xxh, modnam, &
&                    opac_coarse=>opac, q1=>qq1, q2=>qq2, ndiv_calc, &
&                    xj_smoothed=>xj_cmf_coarse, xj_save, &
&                    kcmf_start, kcmf_end, restart_cmfall, &
&                    freold, ifreold1
USE nlte_var, only : z,p,tauz1,tau1,pp1,w0=>vp,w0last,w1=>vp1, metals_converged, &
                     kcmf_start, kcmf_end
USE photstruc, only : chibar_h
USE nlte_app, only : teff, natom, abund, summas
USE cmf_all_var, only: vref, nftot, lam_cmf, indfre_cmf, &
&                    sumopal_cmf, sumetal_cmf, sumtaus_outer, &
&                    xj_cmf, alo_cmf, u_outer, sline_outer, xhinner_cmf, &
&                    fedd, gedd, geddp, xh_cmf, hbound1, nbound1, hboundnd, nboundnd, &
&                    optout_xj_xh, wavblue, wavred, ice, icearr
USE tcorr_var, only: temp_converged
 
implicit none
integer(i4b), parameter :: nd1=id_ndept, np=id_npoin

integer(i4b), parameter :: kmom_start=1 !first kmom_start frequ. point(s) only ray-by-ray

integer(i4b), intent(in) :: nd

real(dp), dimension(nd1), intent(in) :: xne,temp,clfac,r,velo,dvdr,rho,xnh,pressure

real(dp), dimension(nd1), intent(inout) :: taur

real(dp), dimension(nd1), intent(out) :: ferr, dtfcorr

real(dp), allocatable, dimension(:,:), save :: xj_old

real(dp), allocatable, dimension(:,:) :: xj_save_old

real(dp), dimension(nd1) :: sigth, opac, etac, scont, &
&                           x1o, x2o, expo, x1e, x2e, expe, &
&                           opalk, opakk, slinek, ak1, ucmf, pp, aj, ak, alo, &
&                           efac, ydum, uoldk, slinetotk, tauross, &
&                           opaross, oparblue, oparred, pnew, gr, gt, chibar_h_cmf

real(dp), dimension(nd1-1) :: az1

real(dp), dimension(nd1-1) :: vcmf, voldk, v1oldk

real(dp), dimension(nd1-1) :: rad_force, cont_force, el_force, ydum1

real(dp), dimension(nd1+1) :: ah, an, hold, ahweight, ah_int, r2, xh, ah_ray, errh

real(dp), dimension(nd1) ::  xh1, xhlast, xh_int, maxchange_xj, aj_ray, errj, jold

real(dp), dimension(nd1) :: xh_int_blue, xh_int_red

real(dp), dimension(nd1,np-1) :: ubluwi, uold

real(dp), dimension(nd1-1,np-1) :: vbluwi, vold, v1old

real(dp), dimension(np-1) :: iminusold

integer(i4b), dimension(np-1) ::  ncount, nstep

logical, dimension(np-1) ::  in_line_old

real(dp) :: rl, xmue, xmue2, deltax, xlambda, aic, dbdr, e1, e2, de, dei, aux, &
&           opa1, sl1, taus_outer, mubar, vc, e1obs, e2obs, xj1, xj2, &
&           summ, summ1, xjapp, mean, dev, mu, sumabu, csound, devp, hbnd, nbnd

real(dp) :: y1, rdum, vdum, xfac, errmax, dtdr, &
&           nom_flux, fluxdif, fldiff_blue, fldiff_red, &
&           planck, planck_blue, planck_red, tw, &
&           vc1, vcend, wnue1, wnue2, weight, &
&           fre1,fre2, wlast_blue, wlast_red, const, corrfc_cmf, &
&           xx1o, xx2o, xexpo, opacnd, opakknd, desh, xhmin

real(dp) :: bnue, dbdt, hplus, eblue, ered, qx, qx1

integer(i4b) :: jp, lmax, l, lz, ind_old, ind, k, kk, k1, kblue, kred, nsum, &
                ifretest, i, &
                counter_u, counter_uneg, counter_v, counter_jneg, &
                counter_alo_zero

logical :: inversion, cont_blue_thick, start, opt_ray, jneg, newfreq

integer(i4b) :: icounter

character (len=3) :: counter

data icounter/0/

data start/.true./

!check that xj processed by obsfram treatment
if(xj(1,ifre+1).ne.1) stop ' XJ = XJ (CMF) IN CMF_COMPLETE'

!recalculate pressure

sumabu=0.
!note that summas/sumabu are here different from nlte/totout,
!due to different normalization
do k=1,natom
    sumabu=sumabu+abund(k)/abund(1) 
! change of normalization;
! abund' = abund/sum(abund) with abund N_k/N_H
! abund  = abund' * sum(abund)
! sum_abund = abund/abund' = (e.g.) abund(1)/abund'(1) = 1./abund'(1)
enddo
print* 
print*,' sum(abundances, w.r.t. H) = ',sumabu  
!recalculate pressure
devp=1.d100
do l=1,nd
  mu=summas/(1.d0+xne(l)/xnh(l)/sumabu) !note that xnh only approx.
  csound = akb / (mu*amh)  
  pnew(l) = rho(l)*temp(l)*csound   
  devp=min(devp,abs(1.-pnew(l)/pressure(l)))
enddo

print*          
print*,' max dev. between p and pnew = ',devp
!now, pnew is the actual pressure
if(devp.gt.0.05) stop ' deviation in pressure too large'

!restart_cmfall=.false.  !for tests with restart-input XJ(obs. frame)
if(start) then
  call prep_convol
  if(restart_cmfall) then
! restart 
    if(.not.allocated(xj_save)) stop ' restart and xj_save not allocated in cmf_complete!'
    allocate(xj_smoothed(nd1,ifre))
! make sure that frequency grid has not changed! this should never happen
    ifretest=ubound(xj_save,2) 
    if(ifretest.ne.ifre) then
      print*,'old:',ifretest,' new:',ifre
      stop ' IFRE HAS CHANGED IN CMF_COMPLETE (1)!!!!'
    endif  
! smoothing not required, since input value has been smoothed before saving
    xj_smoothed=xj_save
    print*
    print*,' restart of CMF_COMPLETE with KCMF_START/END = ',kcmf_start, kcmf_end
    start=.false.
  else
    allocate(xj_smoothed(nd1,ifre), xj_save(nd1,ifre))
!    allocate(xj_smoothed(nd1,ifre)) !for tests with restart-input XJ(obs. frame)
! smooth XJ for further use in Thomson emissivity
! here: XJ contains only obs. frame quantities (after CONT)
!       smoothed XJ overwrites XJ_SMOOTH -> XJ_CMF_COARSE 
    call xj_smooth(2)
! start set to .false. after end of freloop
  endif 
else
! make sure that frequency grid has not changed
! within iteration, frequency grid can change, due to temp-update  
  newfreq=.false.
  ifretest=ubound(xj_smoothed,2)
  if(ifretest.ne.ifreold1) then
    print*,ifretest, ifreold1, ifre
    stop ' ifretest ne ifreold1'
  endif
! JO March 2022
! when, e.g., Helium has a very low abundance, the main ions of other
! elements might change during the T-iteration, and thus the freq. grid
! might change. In this case, xj_save and xj_smooth need to be adapted  
  if(ifretest.ne.ifre) then
    print*,' old:',ifretest,' new:',ifre
    print*,' IFRE HAS CHANGED IN CMF_COMPLETE (2)!!!!'
    print*
    newfreq=.true.
  else 
    do k=1,ifre
      if(abs(1.d0-fre(k)/freold(k)).gt.1.d-15) then
        print*,k,' ',fre(k),' ',freold(k)
        print*,' FREQUENCY GRID HAS CHANGED IN CMF_COMPLETE (2)!!!!'
        print*
        newfreq=.true.
        exit    
      endif
    enddo  
  endif 
  if(newfreq) then
    allocate(xj_save_old(nd1,ifreold1))
    xj_save_old=xj_save         
    deallocate(xj_smoothed, xj_save)
    allocate(xj_smoothed(nd1,ifre), xj_save(nd1,ifre))
! as in optcmf_restart
    xj_save(:,1)=xj_save_old(:,1)
    xj_save(:,ifre)=xj_save_old(:,ifretest)
    outer: do i=2,ifre-1
    do k=1,ifreold1-1
    if (fre(i).gt.freold(k) .and. fre(i).le.freold(k+1)) then
      qx=log10(fre(i)/freold(k))/log10(freold(k+1)/freold(k))
      qx1=1.d0-qx
      xj_save(:,i)=qx1*log10(xj_save_old(:,k))+qx*log10(xj_save_old(:,k+1))
      xj_save(:,i)=10.d0**xj_save(:,i)
      cycle outer                      
    endif                        
    enddo                        
    stop ' fre not found in freold (subr. cmf_complete)'
    enddo outer
    ifreold1=ifre
    freold=fre
    print*,' cmf_complete: XJ_SAVE (smoothed) interpolated from previous values!'
    print*

    eblue=1.d8/wavblue
    ered=1.d8/wavred

    do i=1,ifre
      if(fre(i).gt.ered) exit
    enddo
    kcmf_start=i

    do i=ifre,1,-1
      if(fre(i).lt.eblue) exit
    enddo
    kcmf_end=i

    if(abs(1.-fre(kcmf_start)/ered).gt.0.1 .or. abs(1.-fre(kcmf_end)/eblue).gt.0.1) &
    stop ' Problems with kcmf_start or kcmf_end in cmf_complete'

    print*,' xj_save updated (re-dimensioned) in cmf_complete'
    print*

    deallocate(ice,icearr)
    call prep_convol
    
  endif           
             
! smooth XJ_SAVE for further use in Thomson emissivity
! here: XJ_SAVE contains CMF + obs. frame quantities
!       (from previous call of INTERPOL_CMF_TO_COARSE)
!       smoothed XJ_SAVE overwrites XJ_SMOOTH -> XJ_CMF_COARSE
  call xj_smooth(3)
endif  
                    
! line-independent quantities
if (start_cmf_simple) then
! might have already been calculated in cmf_simple (nlte_approx)
  
! cmf_simple should have been called before
! (if restart and metals_converged = T, this has not happened)
! but: for very low z, no cmf-line might have been calculated.
! THUS: we comment this stop statement.  
!if (optphotlines .and. .not. metals_converged) &
!& stop ' optphotlines = T, metals_converged = F and start_cmf_simple = T in cmf_complete'

jploop1: do jp = 1,np - 1  
          lmax = min0(np+1-jp,nd)  
          lz = lmax - 1  

          do l = 1,lz  
               tauz1(l,jp) = 1.d0/ (z(l,jp)-z(l+1,jp))  
          end do

          do l = 2,lz  
               tau1(l,jp) = 2.d0/ (z(l-1,jp)-z(l+1,jp))  
          end do
!
!     pp(l)=velocity gradient,projected on the present ray
!
          do l = 1,lmax  
               rl = r(l)  
               xmue = z(l,jp)/rl  
               xmue2 = xmue*xmue  
               pp1(l,jp) = xmue2*dvdr(l)+ (1.d0-xmue2)*velo(l)/rl
          end do

end do jploop1

start_cmf_simple=.false.
endif

sigth=xne*amh*sigmae/clfac !corrected

do l=1,nd
rl=r(l)/rtau23
! average mu to convert cmf to obs frequencies
! assuming linear limb-darkening:
! Int(1 ... mustar) (a+b mu) dmu =: a + b mubar => mubar = 0.5(1+mustar)
if (rl.ge.1d0) then
  mubar=0.5d0*(1.d0+sqrt(1.d0-1.d0/rl**2))
! in the outer wind, I_minus mostly zero
else
!  mubar=0.5d0
  mubar=0.
! in the innermost wind, I_plus and I_minus of same order 
! I_plus comes from bluer, and I_minus from redder frequencies
endif
vc=velo(l)*vmax/clight
efac(l)=1.d0/(1.d0-mubar*vc)
!print*,l,mubar,efac(l)
enddo

! r2 is r^2 on staggered grid
r2(1)=r(1)*r(1)
do l=1,nd-1
 rdum=0.5*(r(l)+r(l+1)) 
 r2(l+1)=rdum*rdum
enddo
r2(nd+1)=r(nd)*r(nd)
  
deltax=vref/vmax

! for moments equation
gr=dvdr/deltax
gt=velo/r/deltax

if (.not.allocated(xj_cmf)) then
! only once
  allocate(xj_cmf(nd1,nftot), alo_cmf(nd1,nftot),xhinner_cmf(nftot))
  allocate(xj_old(nd1,nftot))
  allocate(xh_cmf(nd1+1,nftot),fedd(nd1,nftot),gedd(nd1-1,nftot),geddp(nd1-1,nftot))
  allocate(hbound1(nftot),nbound1(nftot),hboundnd(nftot),nboundnd(nftot))
  allocate(u_outer(np,nftot))
  xj_old=0.
endif

! ---------------------------------------------
! calculate CORRFC_CMF, by combining the ranges from the approx. treatment
! (lambda > wavblue and lambda < wavred), and the CMF range
!
! note: integration weights checked; perfect for blue and red range;
!       accuracy for intermediate range O(10^-12)
!
!preparation of freq. integration weights
vc1=lam_cmf(2)/lam_cmf(1)-1.d0
vcend=lam_cmf(nftot)/lam_cmf(nftot-1)-1.d0
if(abs(vc1-vcend).gt.1.d-14) stop ' error in precision: lam_cmf'
! compact notation
vc1=lam_cmf(1)/lam_cmf(2)
wnue1=(1.d0-vc1)*clight*0.5d8
wnue2=(1.d0-vc1*vc1)*clight*0.5d8

! blue range
fre1=1.d8/lam_cmf(1)
do k=ifre,1,-1
  if(fre(k).lt.fre1) then
    kblue=k
    exit
  endif
enddo
if(kblue+1.gt.ifre) stop 'kblue+1 > ifre'

! renormalize weight for kblue+1
tw=0.d0
do k=ifre,kblue+1,-1
  tw=tw+wfre(k)
enddo  
weight=clight*(fre(ifre)-fre(kblue+1))
wlast_blue=wfre(kblue+1)+weight-tw
!  print*,wlast/wfre(kblue+1)
if(wlast_blue.lt.0.) stop ' wlast < 0 in blue range'

! now integrate over blue range
planck_blue=0.d0
fldiff_blue=0.d0

do k=kblue+1,ifre
  weight=wfre(k)
  if(k.eq.kblue+1) weight=wlast_blue
  planck_blue=planck_blue + bnue(1.d8/fre(k),teff)*weight
  fldiff_blue=fldiff_blue + dbdt(1.d8/fre(k),temp(nd))/opac_coarse(nd,k)*weight
enddo

weight= clight*(fre(kblue+1)-1.d8/lam_cmf(1))  
if(weight.le.0.) stop ' weight < 0 (1)'
rdum=0.5d0*(bnue(1.d8/fre(kblue+1),teff)+bnue(lam_cmf(1),teff))
planck_blue=planck_blue + rdum*weight 
rdum=0.5d0*(dbdt(1.d8/fre(kblue+1),temp(nd)) + dbdt(lam_cmf(1),temp(nd)))/ &
&           opac_coarse(nd,kblue+1) !approximation (opakk approx opac_coarse)
fldiff_blue=fldiff_blue + rdum*weight

! red range
fre2=1.d8/lam_cmf(nftot)
do k=1,ifre
  if(fre(k).gt.fre2) then
    kred=k
    exit
  endif
enddo
if(kred-1.lt.1) stop 'kred-1 < 1'

! renormalize weight for kred-1
tw=0.d0
do k=1,kred-1
  tw=tw+wfre(k)
enddo  
weight=clight*(fre(kred-1)-fre(1))
wlast_red=wfre(kred-1)+weight-tw
!  print*,wlast/wfre(kred-1)
if(wlast_red.lt.0.) stop ' wlast < 0 in blue range'

! now integrate over red range
planck_red=0.d0
fldiff_red=0.d0

do k=1,kred-1
  weight=wfre(k)
  if(k.eq.kred-1) weight=wlast_red
  planck_red=planck_red + bnue(1.d8/fre(k),teff)*weight
  fldiff_red=fldiff_red + dbdt(1.d8/fre(k),temp(nd))/opac_coarse(nd,k)*weight
enddo

weight= clight*(1.d8/lam_cmf(nftot)-fre(kred-1))  
if(weight.le.0.) stop ' weight < 0 (2)'
rdum=0.5d0*(bnue(1.d8/fre(kred-1),teff)+bnue(lam_cmf(nftot),teff))
planck_red=planck_red + rdum*weight
rdum=0.5d0*(dbdt(1.d8/fre(kred-1),temp(nd)) + dbdt(lam_cmf(nftot),temp(nd)))/ &
&           opac_coarse(nd,kred-1) !approximation (opakk approx opac_coarse)
fldiff_red=fldiff_red + rdum*weight

fluxdif = 0.d0
planck = 0.d0
errj=0.
errh=0.

ind_old=0 ! for interpolation

if(opt_ray) then
  print*
  print*,' CMF-transport: RAY-BY-RAY SOLUTION + MOMENTS EQS.'
  print*
else
  print*,' CMF-transport: MOMENTS EQS. ONLY'
  print*
endif  

! cmf_range
do k=1,nftot
! calculate opacity at last depth point, in the same way as for
! *all* depth points below

  xlambda=lam_cmf(k)
! calculate continuum opacities in cmf
! interpolate, to avoid steps in min(opa) 
  ind=indfre_cmf(k)
  if(ind.ne.ind_old) then
!prepare interpolation weights (lamgrid always in between ind and ind-1)
    e1=fre(ind)
    e2=fre(ind-1)  
    de=1.d0/log10(e1/e2)
    xx1o=opac_nolines(nd,ind)   
    xx2o=opac_nolines(nd,ind-1)   
    xexpo=log10(xx1o/xx2o)*de
    xx1o=xx1o*sr ! before, sr cancels out
    ind_old=ind
  endif    
! actual interpolation (log-log):
! x = x1*(nu/nu1)^[(log10(x1/x2)/log10(nu1/nu2)]
  dei=1.d8/(xlambda*e1) ! here we use the shifted frequencies
  opacnd=xx1o*dei**xexpo ! clumping corrected, includes opath
  opakknd=sumopal_cmf(nd,k)+opacnd

! no check for inversion

  if(k.eq.1) then
     weight=wnue1/lam_cmf(1)
  else if (k.eq.nftot) then
     weight=wnue1/lam_cmf(nftot-1)
  else
     weight=wnue2/lam_cmf(k-1)
  endif

! lower flux in diffusion approx.
  y1  = dbdt(lam_cmf(k),temp(nd))/opakknd  
  fluxdif = fluxdif + y1 * weight 

! Planck function at Teff
!y1  = bnue(lam_cmf(k),temp(nd))  
  y1  = bnue(lam_cmf(k),teff)  
  planck = planck + y1 * weight 
    
enddo

! Perform some checks and calculate corrfc_cmf
y1=sigsb/pi*teff**4
nom_flux=sigsb/(4.*pi) * teff**4 * rtau23**2
! corresponds to BNUECON/4 in nlte.f90; correction, because r^2*flux = const
! and we compare to the nominal flux (also w.r.t. errors) at the lower boundary
! which is larger by a factor (rstar/sr)^2=rtau23^2
dtdr = (temp(nd)-temp(nd-1))/ (r(nd-1)-r(nd))  

! check integration of Planck-function (i.e., freq. weights)
planck = planck+planck_blue+planck_red
errmax = 1.d0-planck/y1
print *
print *,' contribution to Planck-function from blue range = ',planck_blue/planck
print *,' contribution to Planck-function from  red range = ',planck_red/planck
print *
print *,' Error in integrated Planck-function (total range) = ',errmax
print *

if(abs(errmax).gt.1.d-4) stop ' something wrong with integration weights (cmf_complete)'

!calculate corrfc_cmf

! Note: opac_coarse without SR
const=corrfc*dtdr/(3.d0*sr)
fldiff_blue=fldiff_blue*const
fldiff_red=fldiff_red*const

const=dtdr/(3.d0)
fluxdif=fluxdif*const

corrfc_cmf=(nom_flux-fldiff_blue-fldiff_red)/fluxdif

print*,' CORRFC_CMF = ',corrfc_cmf 
print*
fluxdif=fluxdif*corrfc_cmf

fluxdif=fluxdif+fldiff_blue+fldiff_red
print *,' contribution to fluxdiff from blue range = ',fldiff_blue/fluxdif
print *,' contribution to fluxdiff from  red range = ',fldiff_red/fluxdif
print *

!for tests (reformulate, since xxh now on grid)
!do k=1,ifre
!  xh_int=xh_int+xxh(:,k)*wfre(k)
!  const=(r(nd-1)-r(nd))*0.5d0*(opac_coarse(nd,k)+opac_coarse(nd-1,k))*sr 
!  const=(r(nd-1)-r(nd))*opac_coarse(nd,k)*sr 
!  rdum=(xj(nd,k)-xj(nd-1,k))/const
!  vdum=(bnue(1.d8/fre(k),temp(nd))-bnue(1.d8/fre(k),temp(nd-1)))/const
!  errmax=dbdt(1.d8/fre(k),temp(nd))*dtdr/(opac_coarse(nd,k)*sr)
!  fluxdif=errmax/3.d0
!  write(*,fmt='(i4,2x,f12.4,2x,3(e10.4,2x),f6.3,2x,e10.4,2x,f6.3)') &
!&  k,1.d8/fre(k),const,rdum,vdum,rdum/vdum,errmax,errmax/vdum   
!  write(*,fmt='(i4,2x,f12.4,2x,3(e10.4,2x),f6.3,2x,f6.3,2x,f6.3)') &
!&  k,1.d8/fre(k),const,vdum,errmax,errmax/rdum,errmax/vdum,rdum/(3.*xxh(nd,k))   
!enddo  

! ---------------------------------------------

! initial values for approximation of Iminus 
iminusold=0.
in_line_old=0
ncount=0
cont_blue_thick=.false.

do jp=1,np-1
  nstep(jp)=(1.-z(1,jp)/r(1))/deltax
enddo

ind_old=0 ! for interpolation

rad_force = 0.d0
cont_force = 0.d0
el_force = 0.d0
ah_int = 0.d0
xh_int = 0.d0
chibar_h = 0.d0
fluxdif = 0.d0
opaross = 0.d0

counter_u=0
counter_v=0
counter_uneg=0
counter_jneg=0
counter_alo_zero=0

freloop: do k=1,nftot

  xlambda=lam_cmf(k)

  if(mod(k,5000).eq.0) print*,k,' cmf_transport at ',xlambda

! calculate continuum opacities/emissivities in cmf
! interpolate, to avoid steps in min(opa,eta) 
  ind=indfre_cmf(k)
  if(ind.ne.ind_old) then
!prepare interpolation weights (lamgrid always in between ind and ind-1)
  e1=fre(ind)
  e2=fre(ind-1)  
  de=1.d0/log10(e1/e2)
  x1o=opac_nolines(:,ind)   
  x2o=opac_nolines(:,ind-1)   
  expo=log10(x1o/x2o)*de
  x1o=x1o*sr ! before, sr cancels out

! if START = .TRUE., i.e., 1st run THEN
! interpolate xj onto corresponding obs. frame frequency
! lam_cmf 'sees' radiation from bluer frequencies, lam_obs = lam_cmf(1-mu*v/c) 
! e.g., at the CMF Lyman edge in the outer wind, the radiation field orginates
! from obs. frame frequencies around 905 A, and has a lower intensity than
! using xj(lamobs=lamcmf)  

! instead of interpolating directly w.r.t xlambda, we
! interpolate only the xj grid values w.r.t. the corresponding obs. frequencies
! (saves a lot of time) 
! could be certainly programmed somewhat more efficient, but programming effort
! not worth doing so

! algorithm checked, does not make a lot of difference compared to the
! approximate approach of using xj(:,ind) and xj(:,ind-1).
! but one never knows ...

! instead of xj we work here with xj_smoothed, to approx. account
! for broadening by Thomson scattering

! if START = .FALSE., no interpolation of XJ_SAVE required,
! when line-frequency between KCMF_START and KCMF_END,
! since the corresponding XJ values are already CMF values.
! Outside this range, we interpolate w.r.t to the same freq. shift as above,
! since here the XJ values are obs. frame values.  

  if (.not.start .and. e2.ge.fre(kcmf_start) .and. e1.le.fre(kcmf_end)) then 
! .not. start and inside range: take direct values
    do l=1,nd
     x1e(l)=etat_nolines(l,ind)+xj_smoothed(l,ind)*sigth(l)
     x2e(l)=etat_nolines(l,ind-1)+xj_smoothed(l,ind-1)*sigth(l)
   enddo
  else
! start or (.not. start and outside range)
! interpolate xj onto cmf-frequencies
    do l=1,nd 
      e2obs=e2*efac(l)
      do kk=ind-1,ifre 
        if(e2obs.ge.fre(kk).and.e2obs.lt.fre(kk+1)) then
! here we interpolate only linearly, to save time 
         xj2=xj_smoothed(l,kk)+ &
           (xj_smoothed(l,kk+1)-xj_smoothed(l,kk))/(fre(kk+1)-fre(kk))*(e2obs-fre(kk))
         goto 10
        endif    
      enddo
      stop ' not found (interpol. of xj2)'

10    e1obs=e1*efac(l)
      do kk=ind,ifre 
        if(e1obs.ge.fre(kk).and.e1obs.lt.fre(kk+1)) then
! here we interpolate only linearly, to save time 
         xj1=xj_smoothed(l,kk)+ &
           (xj_smoothed(l,kk+1)-xj_smoothed(l,kk))/(fre(kk+1)-fre(kk))*(e1obs-fre(kk))
         goto 20
        endif    
      enddo
      stop ' not found (interpol. of xj1)'
    
20    x1e(l)=etat_nolines(l,ind)+xj1*sigth(l)
      x2e(l)=etat_nolines(l,ind-1)+xj2*sigth(l)
    enddo   
  endif
    
  expe=log10(x1e/x2e)*de
  x1e=x1e*sr ! before, sr cancels out
  
  ind_old=ind
  endif    
! actual interpolation (log-log):
! x = x1*(nu/nu1)^[(log10(x1/x2)/log10(nu1/nu2)]
  dei=1.d8/(xlambda*e1) ! here we use the shifted frequencies
  opac=x1o*dei**expo ! clumping corrected, includes opath
  etac=x1e*dei**expe ! clumping corrected, includes approx etath

! for lower boundary
  call diffus (xlambda, temp, r, nd, aic, dbdr)  
! at very first frequency, there is no line (checked in opacity_cmf)
! then, continuum transfer can be performed in obs-frame
  scont=etac/opac
!
!     calculation of continuum radiation field at the bluest frequency
!
  if(k.eq.1) then
    print*
    print* ,' CORRECTION FACTOR FOR CONTINUUM AT WAVBLUE = ',corrfc 
    print*
    call cont2 (nd, np, r, p, z, scont, opac, aic, dbdr, corrfc, ubluwi)
  endif

  opalk=sumopal_cmf(:,k)
  opakk=opalk+opac
  slinek=(sumetal_cmf(:,k)+etac)/opakk

! total line source function, needed for alo in ray_complete
  do l=1,nd
    if (opalk(l) .eq. 0.) then
      if (sumetal_cmf(l,k).ne.0) &
&       stop ' sumopal_cmf = 0 and sumetal_cmf ne 0 in ray_complete'
      slinetotk(l)=0.
    else
      slinetotk(l)=sumetal_cmf(l,k)/opalk(l)
    endif
  enddo  

! check for inversion
  inversion=.false.
  do l=1,nd
    if (opakk(l).le.0.) then
      if(slinek(l).gt.0.) stop ' opakk < 0 and slinek > 0 (subr. cmf_complete)' 
      inversion=.true.
    endif
    if (slinek(l).le.0.) then
      if(opakk(l).gt.0.)  stop ' slinek < 0 and opakk > 0 (subr. cmf_complete)' 
      inversion=.true.
    endif
  enddo  

  if(inversion) then
    print*,' Inversion in cmf_complete at ',xlambda
!  stop ' Inversion in cmf_complete'
!JO  this needs to be treated when occuring
!  hack, needs to be improved
    where (opakk .lt. 0.)
      opalk=0.
      opakk=opac
      slinek=etac/opac
      slinetotk=0.
    endwhere
  endif

  ak1=1./opakk

  do l = 2,nd
    az1(l-1) = 2.d0/ (opakk(l)+opakk(l-1))  
  end do


  optray: if(opt_ray .or. k.eq.1) then
  
  aj = 0.d0  
  ah = 0.d0 
  ak = 0.d0  
  an = 0.d0
  aj = 0.d0  
  alo = 0.d0  
  hbnd = 0.d0 
  nbnd = 0.d0

  opa1=opac(1)
  sl1=sline_outer(k)
  if(sl1.lt.0.) sl1=0.
  taus_outer=sumtaus_outer(k)

  if(k.eq.2 .and. opa1*r(1)/3.d0.gt.1.) cont_blue_thick=.true.
! optically thick continuum only for rho^2 processes

  jploop2: do jp = 1,np - 1  
       lmax = min0(np+1-jp,nd)  
       lz = lmax - 1  

       if(k.eq.1) then
! only for bluewing boundary condition
         ucmf(1) = ubluwi(1,jp)  
         do l = 1,lz  
           ucmf(l+1) = ubluwi(l+1,jp)  
           aux = .5d0* (opac(l)+opac(l+1))  
           vcmf(l) = (ucmf(l+1)-ucmf(l))/aux/ (z(l,jp)-z(l+1,jp))  
         end do
! for ALO
         uold(:,jp)=0.d0
         vold(:,jp)=0.d0
         v1old(:,jp)=0.d0
       else
! restore u and v
         ucmf(1:lmax)=ubluwi(1:lmax,jp)
         vcmf(1:lz)=vbluwi(1:lz,jp)
! for ALO
         uoldk(1:lmax)=uold(1:lmax,jp)
         voldk(1:lz)=vold(1:lz,jp)
         v1oldk(1:lz)=v1old(1:lz,jp)
       endif

! since deltax does not change, this can be calculated in advance (later)
       pp(1:lmax) = pp1(1:lmax,jp)/deltax

       
       call ray_complete(k,jp,z(:,jp),r,nd,np,ucmf,vcmf,opa1,sl1,taus_outer, &
&         iminusold(jp),in_line_old(jp),nstep(jp),ncount(jp), &
&         deltax,aic,dbdr,corrfc_cmf,w0(:,jp),w1(:,jp),w0last, &
&         aj,ah,ak,an,alo,opalk,slinetotk,slinek,ak1,az1,tauz1(:,jp),tau1(:,jp),pp, &
&         uoldk,voldk,v1oldk,hbnd,nbnd,xlambda, &
&         counter_u,counter_uneg,counter_v)

! save current u and v
       ubluwi(1:lmax,jp)=ucmf(1:lmax)
       vbluwi(1:lz,jp)=vcmf(1:lz)
       u_outer(jp,k)=ucmf(1)
!JO     u_outer(jp,k)=ucmf(1)-0.5*iminusold(jp) !This is 0.5*Iplus
! save current uold, vold, v1old for ALO
       uold(1:lmax,jp)=uoldk(1:lmax)
       vold(1:lz,jp) =voldk(1:lz)
       v1old(1:lz,jp)=v1oldk(1:lz)

  end do jploop2

! save ALO
  if(maxval(alo).gt.1.d0+1.d-14 .or. minval(alo) .lt.0.) then
    print*,xlambda,maxval(alo),minval(alo)
    do l=1,nd
      print*,aj(l),' ',alo(l)
    enddo  
    stop ' ALO > 1 or ALO < 0'
  endif
  alo_cmf(:,k)=alo
! for tests
!  if(maxval(alo).le.0) print*,k,'max(alo)=0.'
  
!JO check whether clumping correctly accounted for
! save Eddington factors for moments equations
  fedd(:,k)=ak/aj  ! grid
  where(fedd(:,k).lt.0.05d0) fedd(:,k)=0.05 !for saftey: sphericality factor
  gedd(:,k)=an(2:nd)/ah(2:nd) !staggered grid
  hbound1(k)=ah(1)/aj(1)
  if(hbound1(k).lt.0.) hbound1(k)=0.01d0 !to avoid zero flux
  nbound1(k)=an(1)/aj(1)
  if(nbound1(k).lt.0.) nbound1(k)=0.
  hboundnd(k)=hbnd/aj(nd)
  nboundnd(k)=nbnd/aj(nd)
! geddp defined at staggered grid; thus r^2 J needs to be interpolated
  do l=1,nd-1
    hplus=0.5*(aj(l)*r(l)*r(l)+aj(l+1)*r(l+1)*r(l+1))/(0.5*(r(l)+r(l+1)))**2
    geddp(l,k)=an(l+1)/hplus
  enddo  
  
  aj_ray=aj
  ah_ray=ah
 
  else !not opt_ray

  if(xj_cmf(1,k).eq.0. .or. xh_cmf(1,k).eq.0.) stop ' problems in xj_cmf/xh_cmf in moments iteration'
! use previous solution
  if(xj_cmf(1,k).ne.xj_old(1,k)) stop ' problems in xj_cmf vs. xj_old'
  aj_ray=xj_cmf(:,k)
  ah_ray=xh_cmf(:,k)
  endif optray !opt_ray

  
  hplus=0.5*aic+corrfc_cmf*dbdr/(3.d0*opakk(nd))
! remember: if inversion=.true., opacity and sourcefunction only continuum values

! In case (kmom_start > 1), we start later to avoid problems with
! transition from approx. to exact method. 
! Later on; kmom_start might need to be replaced by k1, or something similar

!  if(k.gt.kmom_start) call moments_cmf(k,nd,aj,aj_ray,ah,ah_ray,jold,hold,r,gr,gt,opakk,slinek,hplus,xlambda,counter_jneg,jneg)

! 
! working with g' (N/J) instead of g = N/H (see notes)
  if(k.gt.kmom_start) call moments_cmf1(k,nd,aj,aj_ray,ah,ah_ray,jold,hold,r,gr,gt,opakk,slinek,hplus,xlambda,counter_jneg,jneg)

  if(.not.opt_ray) then
! JO Oct. 2019: this statement is extremely important, since it ensures
! that no inconsistent alo is used for problematic frequencies where J and H
! have been reset to the "old" ray-by-ray solution inside moments_cmf(1)
  if(jneg) alo_cmf(:,k)=0.
!  do l=1,nd
!    if(abs(1.-aj(l)/aj_ray(l)).gt.0.5) then
!    print*,xlambda,l,aj(l),aj_ray(l),ah(l),ah_ray(l)
!    endif                                               
!  enddo  
  endif
  
  errj=errj+aj/aj_ray
  errh=errh+ah/ah_ray

  jold=aj ! for moments_cmf/next frequency point 
  hold=ah ! for moments_cmf/next frequency point 
  
     
! save J, H, H(ND+1)
! after temp_converged, we use aj_ray to be consistent with ALO and
!to allow for convergence  
!  if(temp_converged) then
!    if(.not.opt_ray) stop ' temp_converged and .not. opt_ray in cmf_all'
!    xj_cmf(:,k)=aj_ray
!  else
!otherwise
    xj_cmf(:,k)=aj  ! note: always used in moments_cmf
                    ! at next frequency point via jold

!  endif  
  xh_cmf(:,k)=ah  ! used in moments_cmf at next frequency point via hold
  xhinner_cmf(k)=ah(nd+1)

  
!interpolate ah (staggered grid) to radial grid (including r^2)
  xh=ah*r2
  xhmin = minval(xh(2:nd))
  desh = 2.d0*abs(xhmin)  

  do l = 2,nd  
    xh(l) = log10(xh(l)+desh)  
  end do

  do l = 1,nd - 2  
    xh(l+1) = q1(l)*xh(l+1) + q2(l)*xh(l+2)  
    xh(l+1) = 10.d0**xh(l+1) - desh  
  end do  
  xh(nd) = xh(nd+1)  

! now, xh = r^2 * H on mesh points (regular grid)
  
! integrate various quantities over frequency within cmf-range

! integration weights from using trapezoidal rule over delta nu, and 
! delta nu = c/lam(1-1/(1+v/c)) for first and last interval, and
! delta_nu = (nu_(k-1)-nu_(k+1)) = c/lam(k-1)*(1-1/(1+v/c)^2)
!            for intermediate intervals
! factor 0.5 (from trapezoidal rule) *clight * 1.d8 included in weights 
! note: 1/(1+v/c)=lam(1)/lam(2) 
  if(k.eq.1) then
     weight=wnue1/lam_cmf(1)
     xh1=xh(1:nd)
  else if (k.eq.nftot) then
     weight=wnue1/lam_cmf(nftot-1)
     xhlast=xh(1:nd)
  else
     weight=wnue2/lam_cmf(k-1)
  endif

  ahweight=ah*weight

! rad. acc calculated on staggered grid
! remember:use ah, ahweight without boundaries, i.e., from 2:ND
! total rad_acc
! JO check for clumping
  ydum=opakk/rho
  ydum1=0.5d0*(ydum(1:nd-1)+ydum(2:nd))
  rad_force=rad_force + ydum1*ahweight(2:nd)

! continuum rad_acc
  ydum=opac/rho
  ydum1=0.5d0*(ydum(1:nd-1)+ydum(2:nd))
  cont_force=cont_force + ydum1*ahweight(2:nd)

! thomson_acc
  ydum=sigth/rho
  ydum1=0.5d0*(ydum(1:nd-1)+ydum(2:nd))
  el_force=el_force + ydum1*ahweight(2:nd)

! integrated flux on staggered grid
  ah_int = ah_int+ahweight

! integrated flux (times r^2) on regular grid
  ydum=xh(1:nd) * weight
  xh_int = xh_int + ydum
  chibar_h=chibar_h + ydum * opakk
!  print*,k,lam_cmf(k),chibar_h(30),opakk(30),xh(30)
  
! prepare calculation of opaross (with dummy tauross)
  do l=1,nd
    tauross(l) = dbdt(lam_cmf(k),temp(l))/opakk(l)
  enddo
  opaross = opaross + tauross* weight
    
! lower flux in diffusion approx.
!  y1  = dbdt(lam_cmf(k),temp(nd))/opakk(nd)  
  y1 = tauross(nd)
  fluxdif = fluxdif + y1 * weight 

! check for frequencies where alo is zero at ALL depth points
  if(k.gt.1. .and. sum(alo_cmf(:,k)).eq.0.) then
    if(opt_ray) then
! should not happen (alo = 0 only at certain, but not all depth points),
! except for first cmf iteration (for a couple of first frequencies,
! alo=0 because uold < 0 and thus set to zero! 
      print*,k,xlambda,' opt_ray = T and alo = 0 at all depth points' 
    endif
    counter_alo_zero=counter_alo_zero+1
  endif
enddo freloop
! end of frequency looop

print* 
print*,' total number of cmf frequency points = ',nftot
print*,' number of resets J->J_ray = ',counter_jneg
print*,' number of calculated u= ',counter_u
print*,' number of resets u->u_new = ',counter_uneg
!currently, no reset of v
!print*,' number of resets abs(v)->u = ',counter_v
print* 
print*,' opt_ray = ',opt_ray  
! counter_alo_zero should be ge counter_jneg, but not signficantly
!(since for consecutive opt_ray=F the previous zeros are NOT reset)
print*,' number of freq. points with ALO set to zero = ',counter_alo_zero
print*  
print*
print*,' average ratio between moments (moments eq. vs. ray-by-ray solution)'
print*
print*,'  L   <Jnu(mom)/Jnu(ray)>  <Hnu(mom)/Hnu(ray)>'
do l=1,nd
  write(*,fmt='(i3,2(7x,f10.5))') l,errj(l)/nftot,errh(l)/nftot
enddo  


if(xj_old(1,1).ne.0.) then
  print*
  print*,' max. change in xj_cmf as a function of depth'
  do l=1,nd
    maxchange_xj(l)=maxval(abs(1.-xj_cmf(l,:)/xj_old(l,:)))
    print*,l,' ',maxchange_xj(l)
  enddo
  print*,' overall max. change: ',maxval(maxchange_xj)
  print* 
!for tests  
!  open(1,file=trim(modnam)//'/xj_comp')  
!  do k=1,nftot
!    write(1,fmt='(i8,3(2x,G12.6))') k,lam_cmf(k),xj_old(33,k),xj_cmf(33,k)
!  enddo
!  close(1)
!  print*,' file xj_comp written to model'
endif

xj_old=xj_cmf !for next iteration

start=.false. 
!first-run condition ends here

! Test that approximate treatment and cmf transport are consistent
!----------------------------------------------------------------------------
! (particularly regarding Iminus). Note that there cannot be a perfect
! consistency, since the approx. treatment neglects any freq. shifts.
! Thus, freq. ranges which are dominated by cont. processes should agree
! at *identical* frequencies, xobs = xcmf (i.e., optically thick continua), &
! whilst ranges dominated by photospheric processes should agree when
! xobs=xcmf/(1-v/c) (Lyman or Balmer edges).
! Additionally, there are numerical problems (formal solution for CMF, compared
! to moments solution for obsframe.)

! Thus, overall we require a 10% accuracy
! check first and last interval over 3 vinf
k1=3.*vmax/vref

!close to blue edge
summ=0.
summ1=0.
nsum=0

print*
if(cont_blue_thick) print*,' blue cont. optically thick, no freq. shift in test'

do k=1,k1
  xlambda=1.d8/lam_cmf(k) ! Kayser
  if(.not.cont_blue_thick) xlambda=xlambda*(1.+vmax/clight) ! corrected for shift 
  do l=ifre,2,-1
    if (fre(l).ge.xlambda .and. fre(l-1) .lt. xlambda) then
      xjapp=xj(1,l)*(xj(1,l)/xj(1,l-1))**((xlambda-fre(l))/(fre(l)-fre(l-1)))
! log interpol, since app. radiation field on coarse grid
!      print*,fre(l),xlambda,fre(l-1),xjapp,xj_cmf(1,k)
      nsum=nsum+1
      xjapp=xj_cmf(1,k)/xjapp
      summ=summ+xjapp
      summ1=summ1+xjapp**2
      exit
    endif  
  enddo
enddo

mean=summ/nsum
dev=sqrt(summ1/nsum-mean**2)

write(*,fmt='(''  close to blue edge: <j_cmf(1)/j_app(1)> = '',f7.2,'' +/- '',f8.3)') &
& mean, dev

if((mean-3.*dev).gt.1.1  .or. (mean+3.*dev).lt.0.9) then
print*,' CMF AND APPROX. TRANSFER NOT CONSISTENT AT BLUE EDGE!!!!'
print*,' CMF AND APPROX. TRANSFER NOT CONSISTENT AT BLUE EDGE!!!!'
print*,' CMF AND APPROX. TRANSFER NOT CONSISTENT AT BLUE EDGE!!!!'
print*,' should happen only when continuum becomes optically thick around L=5!'
print*,' check FLUXCONT!!!'
! stop ' CMF and approx. transfer not consistent at blue edge'
endif

!close to red edge
summ=0.
summ1=0.
nsum=0
do k=nftot,nftot-k1,-1
  xlambda=1.d8/lam_cmf(k)*(1.+vmax/clight) ! Kayser, always corrected for shift
  do l=2,ifre
    if (fre(l).ge.xlambda .and. fre(l-1) .lt. xlambda) then
! here we can afford a linear interpolation (Rayleigh Jeans)
      xjapp=xj(1,l)+(xj(1,l)-xj(1,l-1))/(fre(l)-fre(l-1))*(xlambda-fre(l))
!      print*,fre(l),xlambda,fre(l-1),xjapp,xj_cmf(1,k)
      nsum=nsum+1
      xjapp=xj_cmf(1,k)/xjapp
      summ=summ+xjapp
      summ1=summ1+xjapp**2
      exit
    endif  
  enddo
enddo
mean=summ/nsum
dev=sqrt(summ1/nsum-mean**2)
write(*,fmt='(''  close to  red edge: <j_cmf(1)/j_app(1)> = '',f7.2,'' +/- '',f8.3)') &
& mean,dev
print*


!JO comment out for tests
if((mean-3.*dev).gt.1.1  .or. (mean+3.*dev).lt.0.9) &
& stop ' CMF and approx. transfer not consistent at red edge'

!----------------------------------------------------------------------------

! remember: in the following, xh is xh*r^2 on regular grid, contrasted to ah

! save xj(cmf), xj(approx), xh(nd+1) at specific radii
xhinner_cmf=xhinner_cmf*r2(nd)

if(optout_xj_xh) then

icounter=icounter+1
write(counter,fmt='(i3)') icounter

!open(1,file=trim(modnam)//'/out_xj_cmf_'//adjustl(counter))
open(1,file=trim(modnam)//'/out_xj_cmf')
do k=1,nftot
  write(1,110) lam_cmf(k),xj_cmf(1,k),xj_cmf(43,k),xj_cmf(nd,k),xhinner_cmf(k)
!  write(1,110) lam_cmf(k),xj_cmf(20,k),xj_cmf(35,k),xj_cmf(45,k),xhinner_cmf(k)
enddo  
close(1)

open(1,file=trim(modnam)//'/out_xj_app')
do k=1,ifre
  write(1,110) 1.d8/fre(k),xj(1,k),xj(43,k),xj(nd,k),xxh(nd,k)
enddo
close(1)

!open(1,file=trim(modnam)//'/out_xj_smoothed')
!do k=1,ifre
!  write(1,110) 1.d8/fre(k),xj_smoothed(1,k),xj_smoothed(43,k),xj_smoothed(nd,k),xxh(nd,k)
!enddo
!close(1)


!if xh_cmf has been calculated
!open(1,file=trim(modnam)//'/out_xh_cmf')
!do k=1,nftot
!  write(1,110) lam_cmf(k),xh_cmf(1,k),xh_cmf(43,k),xh_cmf(nd,k),xhinner_cmf(k)
!enddo  
!close(1)

!open(1,file=trim(modnam)//'/out_xh_app')
!do k=1,ifre
!  call diffus (1.d8/fre(k), temp, r, nd, aic, dbdr)  
!  write(1,110) 1.d8/fre(k),xxh(1,k),xxh(43,k),xxh(nd,k), &
!    dbdr/(3.d0*opac_coarse(nd,k)*sr)*corrfc
!enddo  
!close(1)
endif

! calculate and save rad_forces (within wavblue ... wavred)
xfac=(4.*pi/clight)/sr

rad_force  = xfac*rad_force
cont_force = xfac*cont_force
el_force =   xfac*el_force*sr

chibar_h=chibar_h/sr ! since opakk includes SR

print*
print*,' L  grad(stag) grad(grid) grad(nom_flux) chibar_h(cmf-range)'  
open(1,file=trim(modnam)//'/out_rad_force')
  do l=1,nd-1
! now on approximate staggered grid:
! accelerations from 1:ND-1, fluxes (ah*r^2) from 2:ND 
    rdum = (r(l)+r(l+1) )*0.5
    vdum = velo(l) + (velo(l+1)-velo(l))/(r(l+1)-r(l)) * (rdum-r(l))
    ah_int(l+1)=ah_int(l+1)*r2(l+1)
    write(1,105) rdum,vdum,rad_force(l),cont_force(l),el_force(l),ah_int(l+1) 
! check consistency:
! usually, grad(grid)(l) should lie in between grad(stag)(l-1) and grad(stag)(l)
! differences between grad(grid) and grad(nom_flux) when flux different from nom. flux
! NOTE: int Hnu dnu here only for CMF range. In lower atmosphere of hot stars,
! this can be significantly smaller that nominal flux, since blue range (approx. method)
! might contribute a lot. Thus: H_nom > H_range, and grad(nom_flux) > grad(grid)
    rdum=sigsb*teff**4/clight*(rtau23/r(l))**2
    vdum=4.*pi/(clight*r(l)*r(l))
    chibar_h_cmf(l)=chibar_h(l)/(xh_int(l)*rho(l))
! new column chibar_h_cmf included
    write(*,106),l,rad_force(l),chibar_h(l)/rho(l)*vdum,chibar_h_cmf(l)*rdum,chibar_h_cmf(l)
  enddo
  chibar_h_cmf(nd)=chibar_h(nd)/(xh_int(nd)*rho(nd))
  close(1)
print*

! outer boundary (so far, ah_int no longer used from here on)
ah_int(1)=ah_int(1)*r2(1)
ah_int(nd+1)=ah_int(nd+1)*r2(nd+1)

! now combine freq. ranges: outside wavblue ... wavred, use approx. values
! most quantities have been already calculated above
! xxh is flux*r^2 from approx. solution on regular grid
!
! integrate over blue range
xh_int_blue=0.d0
oparblue=0.d0

do k=kblue+1,ifre
  weight=wfre(k)
  xlambda=1.d8/fre(k)
  if(k.eq.kblue+1) weight=wlast_blue
  xh_int_blue=xh_int_blue + xxh(:,k)*weight
  do l=1,nd
   tauross(l) = dbdt(xlambda,temp(l))/opac_coarse(l,k)
  enddo 
  oparblue=oparblue+tauross*weight
  chibar_h=chibar_h+xxh(:,k)*opac_coarse(:,k)*weight
!  print*,'blue',k,chibar_h(30),opac_coarse(30,k),xxh(30,k)
enddo

weight= clight*(fre(kblue+1)-1.d8/lam_cmf(1))  
if(weight.le.0.) stop ' weight < 0 (3)'
xh_int_blue=xh_int_blue + 0.5d0*(xxh(:,kblue+1)+xh1)*weight
xlambda=1.d8/fre(kblue+1)
do l=1,nd
 tauross(l) = dbdt(xlambda,temp(l))/opac_coarse(l,kblue+1)
enddo 
! assuming that opac_coarse is not too different from opakk at lam_cmf(1)
oparblue=oparblue+tauross*weight
chibar_h=chibar_h + 0.5d0*(xxh(:,kblue+1)+xh1)*opac_coarse(:,kblue+1)*weight


! integrate over red range
xh_int_red=0.d0
oparred=0.d0

do k=1,kred-1
  weight=wfre(k)
  xlambda=1.d8/fre(k)
  if(k.eq.kred-1) weight=wlast_red
  xh_int_red=xh_int_red + xxh(:,k)*weight
  do l=1,nd
   tauross(l) = dbdt(xlambda,temp(l))/opac_coarse(l,k)
  enddo 
  oparred=oparred+tauross*weight
  chibar_h=chibar_h+xxh(:,k)*opac_coarse(:,k)*weight
!  print*,'red',k,chibar_h(30),opac_coarse(30,k),xxh(30,k)
enddo

weight= clight*(1.d8/lam_cmf(nftot)-fre(kred-1))  
if(weight.le.0.) stop ' weight < 0 (3)'
xh_int_red=xh_int_red + 0.5d0*(xxh(:,kred-1)+xhlast)*weight
xlambda=1.d8/fre(kred-1)
do l=1,nd
 tauross(l) = dbdt(xlambda,temp(l))/opac_coarse(l,kred-1)
enddo 
! assuming that opac_coarse is not too different from opakk at lam_cmf(nftot)
oparred=oparred+tauross*weight
chibar_h=chibar_h + 0.5d0*(xxh(:,kred-1)+xhlast)*opac_coarse(:,kred-1)*weight

! calculating total opaross and corresponding tauross
! SR included in opakk, but not in opac_coars
opaross=opaross+(oparblue+oparred)/sr
tauross = 7.21863d-5*temp**3  
opaross = tauross / opaross

summ = 0.d0  
tauross(1) = 0.d0  
do l = 1,nd - 1  
     rdum =  r(l) - r(l+1)  
     summ = summ + rdum * 0.5d0 * (opaross(l)+opaross(l+1))  
     tauross(l+1) = summ  
!     print*,l+1,opaross(l+1)
end do  

summ = 0.d0  
ferr(1) = 0.d0  

!for tests: assuming opaross a power-law in r, we obtain the following result
!    note : similar result when assuming opaross exp(-(r-r0)/H)
!in both cases, the derived tauross is a bit *lower* than derived from the
!   trapezoidal rule. Since we overestimate the temperature, this does
!   not help. We would need a larger tauross (exact) vs. tauross(trapez)
!do l = 1,nd - 1  
!     rdum = log(opaross(l)/opaross(l+1))/log(r(l)/r(l+1))  
!     summ = summ + (opaross(l)*r(l)-opaross(l+1)*r(l+1))/(rdum+1.d0)  
!     ferr(l+1) = summ 
!     print*,l+1,tauross(l+1),ferr(l+1)                                                     
!end do  


! recheck and compare with 'approximate' tauross (ferr used as dummy)
oparred=0.d0
do k=1,ifre
  weight=wfre(k)
  xlambda=1.d8/fre(k)
  do l=1,nd
   ferr(l) = dbdt(xlambda,temp(l))/opac_coarse(l,k)
  enddo 
  oparred=oparred+ferr*weight
enddo
ferr = 7.21863d-5*temp**3  
oparred = ferr / oparred  *sr

summ = 0.d0  
ferr(1) = 0.d0  
do l = 1,nd - 1  
     rdum =  r(l) - r(l+1)  
     summ = summ + rdum * 0.5d0 * (oparred(l)+oparred(l+1))  
!     print*,l+1,oparred(l+1)
     ferr(l+1) = summ  
end do  

print*,'L     taur(CMF)   taur(OPAC) taur(used)'
do l=2,nd
  write(*,fmt='(i3,3(2x,f10.4))') l,log10(tauross(l)),log10(ferr(l)),log10(taur(l))
enddo
print*
!JO: so far, it seems that the flux-corrected CMF temperature is too large in
!photospheric regions.
!In any case we need to use the tauross as calculated here (helps a bit),
!but this does not completely solve the problem.
!improving the radial resolution helps, but the grid at higher tau is still
!too coarse (equidistant in LOG(M) and thus LOG(TAUR)
!setting XM(ND) = 1.1 XM(ND-1) also helps, but still:
!-> somewhat erroneous flux -> wrong T(r), though it is conserved
!REAL solution only if moments' equation were solved.

!in new version, this is done (successfully)

taur=tauross
open(1,file=trim(modnam)//'/TAU_ROS',status='unknown',form='formatted')  
rewind 1
do l=1,nd
  read(1,*) rdum, ferr(l) ! ferr is xne_lte
enddo
rewind 1
do l=1,nd
  write(1,*) taur(l),ferr(l)
enddo
close(1)

print*,' TAUROSS UPDATED!!!'
print*

! Perform various checks
! re-check diffusion approx.
const=corrfc_cmf*dtdr/3.d0 ! SR included in opacities
fluxdif=fluxdif*const
errmax = 1.-xh_int(nd)/fluxdif 
print *,' Deviation actual flux/diff approx. in CMF-region = ',errmax  
print *

fluxdif=fluxdif+fldiff_blue+fldiff_red

errmax = 1.d0-fluxdif/nom_flux  
print *,' Error in diffusion approx. = ',errmax  
print * 

if(abs(errmax).gt.1.d-15) stop ' error in diff. approx. > accuracy (cmf_complete)'

print*,' L  H_int(blue)   H_int(red)    H_int(CMF)    H_tot/H_nom'
do l=1,nd
    write(*,120) l,xh_int_blue(l)/nom_flux,xh_int_red(l)/nom_flux, &
&                xh_int(l)/nom_flux,(xh_int(l)+xh_int_blue(l)+xh_int_red(l))/nom_flux
enddo
print*

xh_int=xh_int+xh_int_blue+xh_int_red

! now calculate flux-errors: again, w.r.t. nominal flux at lower boundary,
! which is larger than the flux at Rstar (cf. ROUTINE CONT)

errmax =  1.d0-xh_int(nd)/nom_flux  
print *,' Deviation actual flux/nominal flux (lower boundary) = ',errmax  
print * 

! check flux-conservation/overall flux error
ferr=xh_int/nom_flux-1.d0
errmax = maxval(abs(ferr))
print*,' Max. absolute flux error = ',errmax
print*

!calculate flux-weighted opacity
chibar_h = chibar_h/(xh_int*rho)
open(61,file=trim(modnam)//'/CHIBAR_H_CMF',status='unknown')
rewind 61
do l=1,nd
! here was a bug, cured Sept 12 2019
!  write(61,*) l,chibar_h(l),chibar_h(l)*sigsb*teff**4/clight
  if(abs(1.-chibar_h_cmf(l)/chibar_h(l)).gt.0.1) then
     print*,' WARNING WARNING: chibar_h in CMF and complete range significantly different:'
     print*,l,' ',chibar_h(l),' ',chibar_h_cmf(l)
  endif   
  write(61,*) l,chibar_h(l),chibar_h(l)*sigsb*teff**4/clight*(rtau23/r(l))**2
enddo
close(61)
print*
print*,'file CHIBAR_H_CMF written'
print*

! recheck flux from approx. method (complete freq. range)
xh_int=0.d0
fluxdif=0.d0

do k=1,ifre
  xh_int=xh_int+xxh(:,k)*wfre(k)
  call diffus (1.d8/fre(k), temp, r, nd, aic, dbdr)  
  fluxdif=fluxdif+dbdr/(3.d0*opac_coarse(nd,k)*sr)*corrfc*wfre(k)
enddo  

errmax = maxval(abs(1.d0-xh_int/nom_flux))
print*,' Max. absolute flux error (approx. method) = ',errmax
print*
errmax = abs(1.d0-fluxdif/nom_flux)
print*,' Error in diff. approx. (approx. method) = ',errmax
print*

! temperature correction with respect to flux conservation
call tcorr(taur,ferr,teff,temp,nd,dtfcorr,ndiv_calc+1)

!for tests
!do l=1,nd
!  print*,l,log10(taur(l)),temp(l),xh_int(l),dtfcorr(l)
!enddo  
print*,' cmf_complete done'
print* 

!call test !calculating specific jbars for some tests 
return

100 format(i2,2x,6(e9.3,2x))
105 format(6e15.5)
106 format(i3,4(2x,g10.4))
110 format(f14.4,4(2x,e10.4))
115 format(f14.4,6(2x,e10.4))
120 format(i3,4(2x,e12.6))
125 format(i3,6(2x,e12.6))

end
!
!-----------------------------------------------------------------------
!
subroutine ray_complete(k,jp,z,r,nd,np,u,v,opa1,sl1,taus_outer, &
&           iminusold,in_line_old,nstep,ncount, &
&           deltax,aic,dbdr,corrfc,w0,w1,w0last, &
&           aj,ah,ak,an,alo,opalk,slinetotk,slinek,aakblue,aaz1,tauz1,tau1,pp, &
&           uold,vold,v1old,hbnd,nbnd,xlambda, &
&           counter_u,counter_uneg,counter_v)
! subroutine ray from nlte, adapted for cmf_complete, for freq. k and
! impact parameter k 

! Do no longer use I=I+=I- at last ray; leads to problems (spikes at red
! profile edge) due to transition to continuum, and gives almost no improvement

!JO: in this new version, u is no longer set to zero for neg. values.
!    uold still needs to be set to zero
!    This might need to be re-changed again
!    Had been re-changed

use nlte_type
use nlte_dim, only: id_ndept
use nlte_var, only: vp2, vp3

implicit none
!
!
!     line radiation transfer in the comoving frame from a
!     given source function for a given impact parameter jp.
!     the integration is carried out in space and freq. to
!     yield the mean line intensity
!
!     new version: inlinen aller algebraischen routinen
!
!     calculated are Jnu and consistent ALO
!
!     opalk, slinetotk: total line opacities and source functions
!     slinek: total source function (incl. cont)
!     total opacity not needed (since provided via aakblue etc.)
!
!     .. parameters ..
integer(i4b), parameter :: nd1=id_ndept 
!     ..
!     .. scalar arguments ..
real(dp), intent(in) ::  aic,corrfc,dbdr,deltax,w0last,xlambda, &
&                        opa1,taus_outer,sl1
integer(i4b), intent(in) ::  jp,nd,np  
!     ..
!     .. array arguments ..
real(dp), intent(inout) ::  aj(nd), ak(nd), ah(nd+1), an(nd+1), alo(nd), &
                            u(nd), v(nd-1), uold(nd), vold(nd-1), v1old(nd-1)

real(dp), intent(in) ::  pp(nd),r(nd),opalk(nd),slinetotk(nd),slinek(nd),aakblue(nd), &
&                        aaz1(nd-1),tau1(nd),tauz1(nd),w0(nd),z(nd),w1(nd+1)

real(dp), intent(inout) :: iminusold, hbnd, nbnd
integer(i4b), intent(inout) ::  ncount,counter_u,counter_uneg,counter_v
integer(i4b), intent(in) ::  nstep
logical, intent(inout) ::  in_line_old


!     ..
!     .. local scalars ..
logical :: test, uneghere

real(dp) ::  a,akblue,az1,b,daz,dazm,dbz,dbzm,dt,dtz,dtzm, &
&                 dx,dxz,dxz1,dxzm,dxzm1,h1,help, &
&                 pw,rmax,tauc,taus,s1,tauz,tb1,w,ff1,iplus,uo

integer(i4b) ::  i,k,l,l1,lmax,lz,ni,lmax1
!     ..
!     .. local arrays ..
real(dp) ::       ga(nd1),h(nd1), &
&                 qq(nd1),s(nd1),ta(nd1),tb(nd1),tc(nd1),ub(nd1), &
&                 va(nd1),vb(nd1),uorig(nd1),rhs(nd1),tc1(nd1)

!variables for alo
real(dp) :: ff(nd1),ddia(0:nd1),edia(nd1+1)

!     ..
!     .. intrinsic functions ..
intrinsic exp,max,min0  
!     ..

lmax = min0(np+1-jp,nd)  

if (k.eq.1) then
   do l = 1,lmax  
      pw = w0(l)  
      aj(l) = aj(l) + u(l)*pw  
      ak(l) = ak(l) + u(l)*vp2(l,jp)  
! no action for alo required, since uold=0, and alo=0 set for k=1
   end do

   do l = 1,lmax-1  
      pw = w1(l+1)  
      ah(l+1) = ah(l+1) + v(l)*pw  
      an(l+1) = an(l+1) + v(l)*vp3(l+1,jp)  
   end do
! outer and inner boundary, approximating Iminus(out)=0. 
   ah(1)=ah(1)+u(1)*w1(1) 
   an(1)=an(1)+u(1)*vp3(1,jp) 
   if(lmax.eq.nd) then
      iplus = aic+z(nd)*aakblue(nd)*dbdr*corrfc  
      ah(nd+1)=ah(nd+1)+(iplus-u(nd))*w1(nd+1)
      hbnd=hbnd+u(nd)*w1(nd+1)
      nbnd=nbnd+u(nd)*vp3(nd+1,jp)
   endif  
   return
endif   
  
lz = lmax - 1  
rmax = r(1)  

! uold, vold, v1old initalized outside

! approximate optical depths for outer boundary condition
tauc = opa1*rmax/3.d0  ! optically thick continuum only for rho^2 processes
taus = taus_outer/ (pp(1)*deltax) !since pp divided by deltax 

!
!***  inlined subroutine cmfset
!----------------------------------------------------------------------
!
!***  outer boundary condition  -  1st order with Iminus
!
     akblue = aakblue(1)  
     az1 = aaz1(1)  
     dx = pp(1)*akblue  

! JO modified boundary condition; might need to be changed in subr. ray
! as well, but for line-cores each reasonable condition is sufficient
!
! s1 is source term for outer boundary, corresponds to Imin - P/kappa dImin/dx
     call calc_iminus(s1,xlambda,dx,taus,tauc,sl1,slinek(1), &
&          iminusold,in_line_old,nstep,ncount)
!
! set up of coefficients, starting with outer boundary condition
! 
     s(1) = s1
!     s(1)=0. !if not boundary condition for Iminus
     tc(1) = tauz1(1)*akblue  
     tb(1) = tc(1) + dx + 1.d0  
     ub(1) = dx  
     vb(1) = 0.d0  
     dtzm = tauz1(1)*az1  
!
!***  alo for outer boundary
!
     ff1=0.d0
! for first tests, set ff1 to zero
!JO in case, define ALO for optically thick boundary condition
     ff(1) = ff1  
!!ALO        if (s(1).eq.sc) ff(1) = 0.d0  
!
!***  for g and h, the matrix element s are not different from
!***  inner points
!
     dxzm = .5d0* (pp(1)+pp(2))*az1  
     dxzm1 = 1.d0/ (1.d0+dxzm)  
     dazm = dtzm*dxzm1  
     dbzm = dxzm*dxzm1  
     ga(1) = -dazm  
     h(1) = dbzm  
!
!***  non-boundary points
!
     do l = 2,lz  
          akblue = aakblue(l)  
          az1 = aaz1(l)  
          s(l) = slinek(l)  
          dt = tau1(l)*akblue  
          dtz = tauz1(l)*az1  
          dx = pp(l)*akblue  
          dxz = .5d0* (pp(l)+pp(l+1))*az1  
          dxz1 = 1.d0/ (1.d0+dxz)  
          daz = dtz*dxz1  
          dbz = dxz*dxz1  
          ta(l) = dt*dazm  
          tc(l) = dt*daz  
          tb(l) = ta(l) + tc(l) + dx + 1.d0  
          ff(l) = opalk(l)*akblue  
          ub(l) = dx  
          va(l) = -dt*dbzm  
          vb(l) = dt*dbz  
          ga(l) = -daz  
          h(l) = dbz  
          dazm = daz  
          dbzm = dbz  
     end do

     l = lmax  

     if (lmax.eq.nd) then  
!
!***  inner boundary condition (core rays)  -  only to first order
!***  diffusion-approximation
!
        akblue = aakblue(l)  
        s(l) = aic+z(l)*akblue*dbdr*corrfc  
        iplus=s(l)
        dx = pp(l)*akblue  
!        dt = tauz1(l-1)*akblue  
! CHANGED BY JP Sept. 2013, to be consistent with CONT1
! (if (u(l) linear, then du/dtau = delta_u/delta_tau, 
! and delta_tau = chi_half*dz)
! this is inconsistent with the dx-term, but the latter does not matter anyway
        dt = tauz1(l-1)*aaz1(l-1)  
        ta(l) = dt  
        tb(l) = dt + dx + 1.d0  
        ff(l) = 0.d0  
        ub(l) = dx  
        va(l) = 0.0d0  
        vb(l) = 0.0d0  
     else  
!
!***  inner boundary condition (non-core rays)  -  second order
!
        akblue = aakblue(l)  
        s(l) = slinek(l)  
        tauz = z(lz)   
        dt = akblue/tauz  
        dx = pp(l)*akblue  
        ta(l) = 2.d0*dt*dazm !instead of daz  
        tb(l) = ta(l) + dx + 1.d0  
        ub(l) = dx  
        va(l) = -2.d0*dt*dbzm !instead of dbz
        ff(l) = opalk(l)*akblue  
      endif
!      
!-----------------------------------------------------------------------
!***  continuation of subroutine ray
!
        qq(1) = 0.d0  

        do l = 2,lz  
             qq(l) = va(l)*v1old(l-1) + vb(l)*vold(l)  
        end do

        qq(lmax) = va(lmax)*v1old(lz)  
!
!     uold=ub * uold + qq + ff
!
        do l = 1,lmax
          uold(l) = ub(l)*uold(l) + qq(l) + ff(l)  
        end do
!
!***     now, uold is the solution vector for the inhomgeneous
!***     equation corresponding to u(sl=1)-u(sl=0)
!
!     inlinen von diag mit option 1
!
!
        ddia(0) = 0.d0  
        ta(1) = 0.d0
        tc(lmax) = 0.d0

        do l = 1,lmax  
             ddia(l) = tc(l)/ (tb(l)-ta(l)*ddia(l-1))  
        end do

        edia(lmax+1) = 0.d0  

        do l = lmax,1,-1  
             edia(l) = ta(l)/ (tb(l)-tc(l)*edia(l+1))  
        end do

        do l = 1,lmax  
             uold(l) = uold(l)/ ((1.d0-ddia(l)*edia(l+1))* &
&              (tb(l)-ta(l)*ddia(l-1)))
             uold(l) = max(uold(l),0.d0)  
        end do
!
!***  now,uold is the solution corresponding to u(sl=1)-u(0)
!
        do l = 2,lz  
             l1 = l - 1  
             if (uold(l).eq.0.) then
               v1old(l1) = 0.
               vold(l)=0.
             else  
             v1old(l1) = -ga(l1)*uold(l) + h(l1)*v1old(l1)  
             vold(l) = ga(l)*uold(l) + h(l)*vold(l)  
             endif
        end do

        l = lz  
        v1old(l) = -ga(l)*uold(lmax) + h(l)*v1old(l)  
!
!***  now, v1old, vold are the v's corresponding to delta u
!
!     inlinen von vmalv
!
     b = v(1)  
     qq(1) = vb(1)*b  


     do l = 2,lz  
          a = b  
          b = v(l)  
          qq(l) = va(l)*a + vb(l)*b  
     end do

     qq(lmax) = va(lmax)*b  
!
!     u = ub * u + qq + s
!
!JP: new hack March 2019; avoid problems with Iminus, set to zero
!
     test=.false.     

     uorig=u
     
     do l = 1,lmax  
        u(l) = ub(l)*u(l) + qq(l) + s(l)
     end do
     if(u(1).lt.0) then
        u(1)=ub(1)*uorig(1)
        test=.true.
     endif   

     rhs=u
     tc1=tc

     uneghere=.false.
!
!***     inlining of invtri
!
     tb1 = 1.d0/tb(1)  
     tc(1) = tc(1)*tb1  
     u(1) = u(1)*tb1  

      do l = 2,lz  
          help = tb(l) - ta(l)*tc(l-1)  
          h1 = 1.d0/help  
          tc(l) = tc(l)*h1  
          u(l) = (u(l)+u(l-1)*ta(l))*h1  
     end do

     u(lmax) = (u(lmax)+u(lz)*ta(lmax))/ (tb(lmax)-tc(lz)*ta(lmax))
            
     do l = 1,lz  
          ni = lmax - l  
          u(ni) = u(ni) + tc(ni)*u(ni+1)
     end do


     do l=1,lmax
       if(u(l).lt.0.) then
         uneghere=.true.
         if(l.eq.nd) stop 'u=0 at lower core boundary' 
         if(ub(l).gt.1.) then
           ta(l)=0.
           tc1(l)=0.
           tb(l)=1.
           rhs(l)=uorig(l)           
         else
           if(l.ne.lmax) then
             dxzm=0.5*(pp(l-1)+pp(l))*aaz1(l-1)
             dxz=0.5*(pp(l)+pp(l+1))*aaz1(l)
             ta(l)=ta(l)*(1.+dxzm)
             tc1(l)=tc1(l)*(1.+dxz)
             tb(l)=ta(l)+tc1(l)+1.
             rhs(l)=s(l)
           else
             dxzm=0.5*(pp(l-1)+pp(l))*aaz1(l-1)
             ta(l)=ta(l)*(1.+dxzm)
             tb(l)=ta(l)+1.
             rhs(l)=s(l)
           endif         
         endif
       endif  
     enddo  

     if(uneghere) then
       tc=tc1
       u=rhs
       
       tb1 = 1.d0/tb(1)  
       tc(1) = tc(1)*tb1  
       u(1) = u(1)*tb1  

       do l = 2,lz  
          help = tb(l) - ta(l)*tc(l-1)  
          h1 = 1.d0/help  
          tc(l) = tc(l)*h1  
          u(l) = (u(l)+u(l-1)*ta(l))*h1  
       end do

       u(lmax) = (u(lmax)+u(lz)*ta(lmax))/ (tb(lmax)-tc(lz)*ta(lmax))
            
       do l = 1,lz  
          ni = lmax - l  
          u(ni) = u(ni) + tc(ni)*u(ni+1)

          if(u(ni).lt.0.) then
            print*,xlambda,jp,ni,u(ni)
            stop ' u still <0 after correction'
          endif  
       end do
     endif
          
     counter_u=counter_u+1
     if(uneghere) counter_uneg=counter_uneg+1
!
!     now u is the new field at index k
!
!     inlinen von gmalu und
!     v = h * v + gmalu
!
     b = u(1)  
     do l = 1,lz  
          a = b  
          b = u(l+1)  
          v(l) = h(l)*v(l) + ga(l)* (a-b)  !Note ga(l) = -daz!
     end do

     goto 1113

! until further evidence, this loop is not executed 
     if(uneghere) then
     do l=1,lz
       help=0.5*(u(l)+u(l+1))
       if(abs(v(l)).gt.1.2*help) then
         counter_v=counter_v+1
         if(v(l).gt.0.) then
           v(l)=help
         else
           v(l)=-help
         endif  
       endif
     enddo
     endif
!
!     now v is the new field at index k
!
!     adding the new u to aj
!

1113 continue
       
     do l = 1,lmax  
          if(u(l).le.0.) stop ' u(l) < 0 in final loop (subr. ray_complete)'
! no update required, since u needs to be nelected for Jnu, as well as uold for alo
          pw = w0(l)  
          aj(l) = aj(l) +   u(l)*pw  
          ak(l) = ak(l) +   u(l)*vp2(l,jp)  
!JO April 2018: previously, there was a cycle statement; however, such
!         contribution can be essential. Safetyfactor 0.9 anyway! 
!         Since uold always < 1., this can happen only if SL > u;
!         in this case, the effective ALO also becomes < 0.9 
          if(uold(l)*slinetotk(l) .gt. u(l)) then
            uo=0.9*u(l)/slinetotk(l) !do not modify uold itself
            if(uo.gt.0.9) stop ' something really wrotten with uold'
            alo(l) = alo(l) + uo*pw
          else
            alo(l) = alo(l) + uold(l)*pw
          endif  
     end do
     
!Calculate H, needed for radiative acc, JS, and N, needed for moments eq. (JP)
!     print*,xlambda
     do l = 1,lmax-1  
          pw = w1(l+1)
          ah(l+1) = ah(l+1) + v(l)*pw  
          an(l+1) = an(l+1) + v(l)*vp3(l+1,jp)  
     enddo

! outer and inner boundary
     if(test) then
! outer boundary reset to Iminus = 0, thus v=u
      ah(1)=ah(1)+u(1)*w1(1) 
      an(1)=an(1)+u(1)*vp3(1,jp) 
     else     
      ah(1)=ah(1)+(u(1)-iminusold)*w1(1) !iminusold corresponds here to actual Iminus
      an(1)=an(1)+(u(1)-iminusold)*vp3(1,jp) !iminusold corresponds here to actual Iminus
     endif

     if(lmax.eq.nd) then
       ah(nd+1)=ah(nd+1)+(iplus-u(nd))*w1(nd+1)
! integrals of u, required for moments equations
       hbnd=hbnd+u(nd)*w1(nd+1)
       nbnd=nbnd+u(nd)*vp3(nd+1,jp)
     endif  

return

! obsolte
!return if not boundary condition for Iminus
!if (jp.ne.np-1) return  
!
!***  final value for p = rmax and z=0 in optically thick case!!
!
!     print *,' optically thick outer boundary at',xlambda  
!
!     pw = w0last  
!     aj(1) = aj(1) + iminus*pw  
!!ALO             if (s1.eq.sl) alo(1) = alo(1) + alo1*pw    
!return  

end
!
!----------------------------------------------------------------------------
!
subroutine moments_cmf(k,nd,aj,aj_ray,ah,ah_ray,jold,hold,r,gr,gt,opakk,slinek,hplus, &
                       xlambda,counter_jneg,jneg)
!
! solves moments equations for given Eddington factors f and g.
! needs to be updated for g' in case of very small H in g=N/H.
! on output/input, J, H, J_old, and H_old have their usual meaning.
! during calculation, these values and eta are assumed to be multiplied by r^2
!
! note: whereas H has dimension ND+1 (only elements 2:ND on staggered grid will
! be used), gedd has dimension ND-1 (only staggered grid)

! as usual, we solve for (-ta, tb, -tc)J = x

!programmed by J.P. during a trip to Tenerife Feb. 2019
!(based on his diploma work, partly incorporting notes from W.-R. Hamann) 
!
USE nlte_type
USE nlte_dim, only: id_ndept
USE cmf_all_var, only: fedd, gedd, hbound1, nbound1, hboundnd, nboundnd

implicit none
integer(i4b), parameter :: nd1=id_ndept

integer(i4b), intent(in) :: k,nd
integer(i4b), intent(inout) :: counter_jneg
logical, intent(out) :: jneg

real(dp), intent(out) :: aj(nd),ah(nd+1)
real(dp), intent(in) ::  jold(nd),hold(nd+1),r(nd),gr(nd),gt(nd), &
&                        opakk(nd),slinek(nd),aj_ray(nd),ah_ray(nd+1)
real(dp), intent(in) :: hplus, xlambda

integer(i4b) :: l,lz,ni

real(dp), dimension(nd1) :: ajold, rr, etak, q
real(dp), dimension(nd1) :: ta,tb,tc,qfp1,qfm1,qfp0,qfm0,blp,blm
real(dp), dimension(nd1-1) :: ahold, r2

real(dp) :: rl,rlm,rlp,dr,opa0,opam,opap,gr0,grm,grp,gt0,gtm,gtp, &
&           rrq,fl,flp,dr1,denp,denm,dtp,dtm,tb1,help,h1

! converting to quantities as used in the moments equations. Here,
! ahold = r^2* h_nue(k-1) has dimension ND-1, i.e., only staggered grid values.
rr=r*r
ajold=rr*jold
etak=rr*opakk*slinek

ahold=hold(2:nd)
do l=1,nd-1
  rl=0.5d0*(r(l)+r(l+1))
  r2(l)=rl*rl
enddo
ahold=ahold*r2

!sphericality factor
q(nd)=1.d0
rrq=1.d0
fl=3.d0-1.d0/fedd(nd,k)
do l=nd-1,1,-1
  rl=r(l)
  rlp=r(l+1)
  flp=fl
  fl=3.d0-1.d0/fedd(l,k)
  rrq=rrq*exp(fl-flp)*(rl/rlp)**((flp*rl-fl*rlp)/(rl-rlp))
  q(l)=rrq/rl/rl
enddo


do l=2,nd-1
  rl=r(l)
  rlm=r(l-1)
  rlp=r(l+1)
  dr1=2.d0/(rlm-rlp)
  opa0=opakk(l)
  opam=.5*(opakk(l-1)+opa0)
  opap=.5*(opakk(l+1)+opa0)
  gr0=gr(l)
  grm=0.5*(gr(l-1)+gr0)
  grp=0.5*(gr(l+1)+gr0)
  gt0=gt(l)
  gtm=0.5*(gt(l-1)+gt0)
  gtp=0.5*(gt(l+1)+gt0)

  denp=opap+(grp-gtp)*gedd(l,k)+gtp
  denm=opam+(grm-gtm)*gedd(l-1,k)+gtm

  dtp=2.d0/(q(l)+q(l+1))/(rl-rlp)
  dtm=2.d0/(q(l)+q(l-1))/(rlm-rl)
! following quantities vectors, to enable re-calculation of H after J
! in future versions, many of them can be made scalar
  qfp1(l)=dtp*q(l+1)*fedd(l+1,k)/denp
  qfm1(l)=dtm*q(l-1)*fedd(l-1,k)/denm
  qfp0(l)=dtp*q(l)*fedd(l,k)/denp
  qfm0(l)=dtm*q(l)*fedd(l,k)/denm

  blp(l)=((grp-gtp)*gedd(l,k-1)+gtp)/denp
  blm(l)=((grm-gtm)*gedd(l-1,k-1)+gtm)/denm
    
  ta(l)=qfm1(l)*dr1
  tc(l)=qfp1(l)*dr1
  tb(l)=(qfp0(l)+qfm0(l))*dr1+(gr0-gt0)*fedd(l,k)+gt0+opa0

  aj(l)=etak(l)+((gr0-gt0)*fedd(l,k-1)+gt0)*ajold(l)- &
&      dr1*(blm(l)*ahold(l-1)-blp(l)*ahold(l))
enddo  

!outer boundary
l=1
dtp=2.d0/(q(l)+q(l+1))/(r(l)-r(l+1))

ta(l)=0.d0
tc(l)=dtp*q(l+1)*fedd(l+1,k)
tb(l)=dtp*q(l)*fedd(l,k)+(gr(l)-gt(l))*nbound1(k)+(gt(l)+opakk(l))*hbound1(k)

aj(l)=((gr(l)-gt(l))*nbound1(k-1)+gt(l)*hbound1(k-1))*ajold(l)

!inner boundary
l=nd
dtm=2.d0/(q(l)+q(l-1))/(r(l-1)-r(l))
ta(l)=dtm*q(l-1)*fedd(l-1,k)
tc(l)=0.
tb(l)=dtm*q(l)*fedd(l,k)+(gr(l)-gt(l))*nboundnd(k)+(gt(l)+opakk(l))*hboundnd(k)

aj(l)=opakk(l)*hplus+((gr(l)-gt(l))*nboundnd(k-1)+gt(l)*hboundnd(k-1))*ajold(l)

!
!***     inlining of invtri; remember: solution vector (aj) overwritten by solution
!
  
lz=nd-1

! same result as when calculated in qp (i.e., precision OK). 

tb1 = 1.d0/tb(1)  
tc(1) = tc(1)*tb1  
aj(1) = aj(1)*tb1  

do l = 2,lz  
     help = tb(l) - ta(l)*tc(l-1)  
     h1 = 1.d0/help  
     tc(l) = tc(l)*h1  
     aj(l) = (aj(l)+aj(l-1)*ta(l))*h1  
end do

aj(nd) = (aj(nd)+aj(lz)*ta(nd))/ (tb(nd)-tc(lz)*ta(nd))

do l = 1,lz  
     ni = nd - l  
     aj(ni) = aj(ni) + tc(ni)*aj(ni+1)  
end do

jneg=.false.
do ni=1,nd
  if(aj(ni).lt.0.)then
    if (.not.jneg) counter_jneg=counter_jneg+1
    jneg=.true.
  endif  
enddo

if(jneg) then
! reset to ray-by-ray solution or to solution of previous iteration
! (which might be ray-by-ray solution)
  aj=aj_ray
  ah=ah_ray
  return
endif  

!now, aj is r^2*J_nu at index k
!let's calculate now r^2 H_nu (staggered grid)
!first defined matrix elements at l=2, corresponding to ah(3), but ahold(2)
! range: 2+1/2 to ND-1+1/2
do l=2,lz
  ah(l+1)=qfp1(l)*aj(l+1)-qfp0(l)*aj(l)+blp(l)*ahold(l)
enddo

!to check consistency: here, we calculate h from the backward quantities,
!should give identical results. Again, first defined matrix elements at l=2,
!but now corresponding to ah(2) and ahold(1)
! range: 1+1/2 to ND-2+1/2

! test can be commented out for production runs
do l=3,lz !(last defined matrix elements at LZ
  help=qfm0(l)*aj(l)-qfm1(l)*aj(l-1)+blm(l)*ahold(l-1)
  if(abs(1.-help/ah(l)).gt.1.d-16) &
&  stop ' subr. moments_cmf: problems in forward and backward matrix elements'
enddo

!this is the point we need
l=2
ah(2)=qfm0(l)*aj(l)-qfm1(l)*aj(l-1)+blm(l)*ahold(l-1)

!now we define H_nu at the outer and innermost points, using the same
!philosophy as above.
ah(1)=hbound1(k)*aj(1)
ah(nd+1)=hplus-hboundnd(k)*aj(nd)

!finally, correct for r^2 term
aj=aj/rr
ah(2:nd)=ah(2:nd)/r2
ah(1)=ah(1)/r(1)**2
! nothing to be done at ND, since r(nd)=1.


return
end
!
!----------------------------------------------------------------------------
!
subroutine moments_cmf1(k,nd,aj,aj_ray,ah,ah_ray,jold,hold,r,gr,gt,opakk,slinek,hplus, &
                       xlambda,counter_jneg,jneg)
!
! solves moments equations for given Eddington factors f and g' = N/J (staggered grid!)
! new method, following Hillier and Miller  
! on output/input, J, H, J_old, and H_old have their usual meaning.
! during calculation, these values and eta are assumed to be multiplied by r^2
!
! note: whereas H has dimension ND+1 (only elements 2:ND on staggered grid will
! be used), geddp has dimension ND-1 (only staggered grid)

! as usual, we solve for (-ta, tb, -tc)J = x

!programmed by J.P. April 2019
!(based on his diploma work, partly incorporting notes from W.-R. Hamann) 
!
USE nlte_type
USE nlte_dim, only: id_ndept
USE cmf_all_var, only: fedd, geddp, hbound1, nbound1, hboundnd, nboundnd

implicit none
integer(i4b), parameter :: nd1=id_ndept

integer(i4b), intent(in) :: k,nd
integer(i4b), intent(inout) :: counter_jneg
logical, intent(out) :: jneg

real(dp), intent(out) :: aj(nd),ah(nd+1)
real(dp), intent(in) ::  jold(nd),hold(nd+1),r(nd),gr(nd),gt(nd), &
&                        opakk(nd),slinek(nd),aj_ray(nd),ah_ray(nd+1)
real(dp), intent(in) :: hplus, xlambda

integer(i4b) :: l,lz,ni

real(dp), dimension(nd1) :: ajold, rr, etak, q
real(dp), dimension(nd1) :: ta,tb,tc,qfp1p,qfm1p,qfp0p,qfm0p,blpr,blmr,blprr,blmrr
real(dp), dimension(nd1-1) :: ahold, r2, ajold_stag

real(dp) :: rl,rlm,rlp,dr,opa0,opam,opap,gr0,grm,grp,gt0,gtm,gtp, &
&           rrq,fl,flp,dr1,denpr,denmr,dtp,dtm,tb1,help,h1

! converting to quantities as used in the moments equations. Here,
! ahold = r^2* h_nue(k-1) has dimension ND-1, i.e., only staggered grid values.
rr=r*r
ajold=rr*jold
etak=rr*opakk*slinek

ahold=hold(2:nd)
do l=1,nd-1
  rl=0.5d0*(r(l)+r(l+1))
  r2(l)=rl*rl
  ajold_stag(l)=0.5d0*(ajold(l)+ajold(l+1))
enddo
ahold=ahold*r2

!sphericality factor
q(nd)=1.d0
rrq=1.d0
fl=3.d0-1.d0/fedd(nd,k)
do l=nd-1,1,-1
  rl=r(l)
  rlp=r(l+1)
  flp=fl
  fl=3.d0-1.d0/fedd(l,k)
  rrq=rrq*exp(fl-flp)*(rl/rlp)**((flp*rl-fl*rlp)/(rl-rlp))
  q(l)=rrq/rl/rl
enddo


do l=2,nd-1
  rl=r(l)
  rlm=r(l-1)
  rlp=r(l+1)
  dr1=2.d0/(rlm-rlp)
  opa0=opakk(l)
  opam=.5*(opakk(l-1)+opa0)
  opap=.5*(opakk(l+1)+opa0)
  gr0=gr(l)
  grm=0.5*(gr(l-1)+gr0)
  grp=0.5*(gr(l+1)+gr0)
  gt0=gt(l)
  gtm=0.5*(gt(l-1)+gt0)
  gtp=0.5*(gt(l+1)+gt0)

  denpr=opap+gtp
  denmr=opam+gtm

  dtp=2.d0/(q(l)+q(l+1))/(rl-rlp)
  dtm=2.d0/(q(l)+q(l-1))/(rlm-rl)
! following quantities vectors, to enable re-calculation of H after J
! in future versions, many of them can be made scalar
  qfp1p(l)=(dtp*q(l+1)*fedd(l+1,k)-(grp-gtp)*geddp(l,k)*0.5d0)/denpr
  qfm1p(l)=(dtm*q(l-1)*fedd(l-1,k)+(grm-gtm)*geddp(l-1,k)*0.5d0)/denmr
  qfp0p(l)=(dtp*q(l)*fedd(l,k)+(grp-gtp)*geddp(l,k)*0.5d0)/denpr
  qfm0p(l)=(dtm*q(l)*fedd(l,k)-(grm-gtm)*geddp(l-1,k)*0.5d0)/denmr

  blpr(l)=((grp-gtp)*geddp(l,k-1))/denpr
  blprr(l)=gtp/denpr
  blmr(l)=((grm-gtm)*geddp(l-1,k-1))/denmr
  blmrr(l)=gtm/denmr
    
  ta(l)=qfm1p(l)*dr1
  tc(l)=qfp1p(l)*dr1
  tb(l)=(qfp0p(l)+qfm0p(l))*dr1+(gr0-gt0)*fedd(l,k)+gt0+opa0

  aj(l)=etak(l)+((gr0-gt0)*fedd(l,k-1)+gt0)*ajold(l)- &
&      dr1*(blmrr(l)*ahold(l-1)-blprr(l)*ahold(l))- &
&      dr1*(blmr(l)*ajold_stag(l-1)-blpr(l)*ajold_stag(l))
enddo  

!outer boundary
l=1
dtp=2.d0/(q(l)+q(l+1))/(r(l)-r(l+1))

ta(l)=0.d0
tc(l)=dtp*q(l+1)*fedd(l+1,k)
tb(l)=dtp*q(l)*fedd(l,k)+(gr(l)-gt(l))*nbound1(k)+(gt(l)+opakk(l))*hbound1(k)

aj(l)=((gr(l)-gt(l))*nbound1(k-1)+gt(l)*hbound1(k-1))*ajold(l)

!inner boundary
l=nd
dtm=2.d0/(q(l)+q(l-1))/(r(l-1)-r(l))
ta(l)=dtm*q(l-1)*fedd(l-1,k)
tc(l)=0.
tb(l)=dtm*q(l)*fedd(l,k)+(gr(l)-gt(l))*nboundnd(k)+(gt(l)+opakk(l))*hboundnd(k)

aj(l)=opakk(l)*hplus+((gr(l)-gt(l))*nboundnd(k-1)+gt(l)*hboundnd(k-1))*ajold(l)

!
!***     inlining of invtri; remember: solution vector (aj) overwritten by solution
!
  
lz=nd-1

! same result as when calculated in qp (i.e., precision OK). 

tb1 = 1.d0/tb(1)  
tc(1) = tc(1)*tb1  
aj(1) = aj(1)*tb1  

do l = 2,lz  
     help = tb(l) - ta(l)*tc(l-1)  
     h1 = 1.d0/help  
     tc(l) = tc(l)*h1  
     aj(l) = (aj(l)+aj(l-1)*ta(l))*h1  
end do

aj(nd) = (aj(nd)+aj(lz)*ta(nd))/ (tb(nd)-tc(lz)*ta(nd))

do l = 1,lz  
     ni = nd - l  
     aj(ni) = aj(ni) + tc(ni)*aj(ni+1)  
end do

jneg=.false.
do ni=1,nd
  if(aj(ni).lt.0.)then
    if (.not.jneg) counter_jneg=counter_jneg+1
    jneg=.true.
  endif  
enddo

if(jneg) then
! reset to ray-by-ray solution or to solution of previous iteration
! (which might be ray-by-ray solution)
  aj=aj_ray
  ah=ah_ray
  return
endif  

!now, aj is r^2*J_nu at index k
!let's calculate now r^2 H_nu (staggered grid)
!first defined matrix elements at l=2, corresponding to ah(3), but ahold(2)
! range: 2+1/2 to ND-1+1/2
do l=2,lz
  ah(l+1)=qfp1p(l)*aj(l+1)-qfp0p(l)*aj(l)+blpr(l)*ajold_stag(l)+blprr(l)*ahold(l)
enddo

!to check consistency: here, we calculate h from the backward quantities,
!should give identical results. Again, first defined matrix elements at l=2,
!but now corresponding to ah(2) and ahold(1)
! range: 1+1/2 to ND-2+1/2

! test can be commented out for production runs
do l=3,lz !(last defined matrix elements at LZ
  help=qfm0p(l)*aj(l)-qfm1p(l)*aj(l-1)+blmr(l)*ajold_stag(l-1)+blmrr(l)*ahold(l-1)
  if(abs(1.-help/ah(l)).gt.1.d-16) then
    print*,help,ah(l)
    stop ' subr. moments_cmf1: problems in forward and backward matrix elements'
  endif
enddo

!this is the point we need
l=2
ah(2)=qfm0p(l)*aj(l)-qfm1p(l)*aj(l-1)+blmr(l)*ajold_stag(l-1)+blmrr(l)*ahold(l-1)

!now we define H_nu at the outer and innermost points, using the same
!philosophy as above.
ah(1)=hbound1(k)*aj(1)
ah(nd+1)=hplus-hboundnd(k)*aj(nd)

!finally, correct for r^2 term
aj=aj/rr
ah(2:nd)=ah(2:nd)/r2
ah(1)=ah(1)/r(1)**2
! nothing to be done at ND, since r(nd)=1.


return
end
!
!----------------------------------------------------------------------------
!
subroutine formal_obs(r,velo1,rmax)
!
! calculates (approx). emergent Eddington flux, normalized to Rstar, 
! using u_outer
!
!JO: in case, needs to be updated for Iminus  
! do we need u or 0.5*Iplus? (as calculated in subr. CMF_COMPLTE?
! In any case, u including Iminus needs to be corrected
! also for negative angles!

!JO first tests have shown that using u or u-0.5*Iminus gives similar results
! beyond the HeII edge, but different fluxes below (as to be expected).
! Interestingly, using u is more similar to cmfgen solutions.
! final decision will depend on comparison with solution from FORMALSOL.
! Thus far, we use u and 'live' with the inconsistency with fluxes from CONT
! below 228 A.  
  
USE nlte_type
USE nlte_dim, only: id_ndept, id_npoin
USE fund_const, only: clight
USE nlte_var, only : modnam, vmax, rtau23
USE nlte_var, only : p, vp1
USE cmf_all_var, only: nftot, vref, lam_cmf, u_outer, h_obs, optout_xj_xh

implicit none
integer(i4b), parameter :: nd=id_ndept, np=id_npoin

real(dp), intent(in) :: velo1, rmax
real(dp), dimension(nd), intent(in) :: r

real(dp), dimension(np) :: lamfac

integer(i4b) :: k, jp, ik, dk, dkmax
real(dp) :: h, rmax2, vc, r1, lam1, u_inter, kint

! usally allocated only one, since called only one (if UNASOL)
! but for tests, called more often. Thus
if(.not.allocated(h_obs)) allocate(h_obs(nftot))
! JO Dec. 2018: to account for end-effects
h_obs=0.

r1=r(1)
rmax2=r1/rtau23

! check consistency
if(abs(1.-rmax2/rmax).gt.1.d-12) stop ' problem with rtau23 in subr. formal_obs'
!print*
!print*,' outermost radius point in units of nominal radius = ',rmax2
rmax2=rmax2*rmax2

vc=velo1*vmax/clight

dk=vmax/vref 

! assuming now that lamcmf is a grid refering to lambda_obs, we calculate
! the corresponding cmf-wavelength where u needs to be evaluated
! according to lam_cmf = lam_obs/(1-mu*v/c) = lam_obs*lamfac(jp)
! Note that lam_obs (for p=Rmax) le lam_cmf le lam_obs/(1-v/c) (for p=0)
do jp=1,np-1
  lamfac(jp)=1.d0/(1.d0-sqrt(1.d0-(p(jp)/r1)**2)*vc)
enddo

do k=1,nftot-(dk+1) ! to account for end effects
  h=0.
  dkmax=dk
  do jp=1,np-1
! assuming v(pmax) = 0. because of symmetry
!JO: modify when outer boundary included
    lam1=lam_cmf(k)*lamfac(jp)
! interpolate u(lam1=lam_cmf(lam_obs)) from u(lam)
! jp=0:    lam1=lam_obs/(1-v/c)
! jp=rmax: lam1=lam_obs
    do ik=dkmax,0,-1
      if(lam_cmf(k+ik+1).gt.lam1 .and. lam_cmf(k+ik).lt.lam1) then
       goto 10    
      endif
    enddo
    print*,k,jp,lam1,lam_cmf(k+ik+1),lam_cmf(k)
    stop 'index not found (subr. formal_obs)'

10  kint=k+ik
! linear interpolation
    u_inter=u_outer(jp,kint) + (u_outer(jp,kint+1)-u_outer(jp,kint))/ &
&           (lam_cmf(kint+1)-lam_cmf(kint))*(lam1-lam_cmf(kint))

! integration
    h=h+u_inter*vp1(1,jp)    
    dkmax=ik  
  enddo  
  h_obs(k)=h*rmax2 
enddo

print* 
print*," H-obs calculated (and in case saved to file 'out_xh_obs')"
print*,' (H-obs corresponds to H_nu * (r(1)/Rstar)^2, contrasted to XXH)'

if(optout_xj_xh) then
  open(1,file=trim(modnam)//'/out_xh_obs')
! here lam_cmf is actually lam_obs (vacuum)
    do k=1,nftot
      write(1,100) lam_cmf(k),h_obs(k)
    enddo  
  close(1)
endif
  
100 format(f14.4,2x,e10.4)

return
end
!
!----------------------------------------------------------------------------
!
subroutine test
!
! calculates Jbar for specified transition, and compares with SA
! (e.g., to check outer boundary condition)  

USE nlte_type
USE nlte_dim, only: id_ndept
USE cmf_all_var, only: icfirst, nftot, index_lam_cmf, lam_cmf, xlam1, id1, &
&                      indexel_inv, nf_cmf,xj_cmf,pweightdop

implicit none
integer(i4b), parameter :: nd=id_ndept

integer(i4b) :: i, j, lastindex, in, index, kk, irest, ml, mu, k, ixmax, &
&               jj, low, up, ik, l

real(dp), dimension(nd) :: jbar

lastindex=icfirst-1
do i=1,nftot
  in=index_lam_cmf(i)
  if (lam_cmf(i).gt.1540 .and. lam_cmf(i).lt.1560) then
!    print*,i,lam_cmf(i),in
    if(in.ne.0) then
      do j=lastindex+1,in
        index=id1(j)
!        print*,j,index
        if (index.ge.1d8) then
! explicit element
          kk = index/1d8
          irest = index - kk*1d8
          ml = irest/10000
          mu = irest - ml*10000
          k=indexel_inv(kk)
! in units of vref = vturb/3 (with vturb=max(vturb,vturbmin)
          ixmax=nf_cmf(k)
!          print*,'exp',xlam1(j),kk,k,ml,mu,ixmax
        else
          k = index/1000000
          irest = index - k*1000000
          jj = irest/100000
          irest = irest - jj*100000
          low = irest/100
          up  = irest-low*100
          ixmax=nf_cmf(k)
          if(k.eq.6.and.low.eq.1.and.up.eq.2) then
            print*,i,lam_cmf(i)
            print*,' bg',xlam1(j),k,low,up,ixmax
            jbar(:)=xj_cmf(:,i)*pweightdop(0,:,k) ! ik=0
            do ik=1,ixmax
               jbar(:)=jbar(:)+(xj_cmf(:,i+ik)+xj_cmf(:,i-ik))*pweightdop(ik,:,k)
            enddo
            do l=1,nd
               print*,l,jbar(l)
            enddo
            print*
          endif
        endif
      enddo  
    endif
  endif
  if(in.ne.0) lastindex=in
enddo

return
end
!
!----------------------------------------------------------------------------
!
subroutine calc_iminus(s1,xlambda,dx,taus,tauc,sl1,slinekblue, &
&   iminusold,in_line_old,nstep,ncount)
!
! completely new outer boundary conditions, which are consistent with
! approximate treatment
!
! s1 is the output quantity, source term for outer boundary condition
! all quantities calculated for ray jp 

! 1. Iminus - P/opa dIminus/dx (standard) or approx Iminus - (Iminus - S) (alternative)
!    old outer boundary condition with Iminus alone resulted in increasing
!    flux at long wavelengths. Origin due to neglect of dIminus/dx and
!    setting u=0 for negative u 
! 2. include opac(1) in total tau for Iminus; otherwise, optically thick
!    continuum shortward of HeII edge incorrect  

! to test outer boundary condition, compare with Iminus=0. (set s1=0.), 
! and compare j(1) with j(2) and xj(1) from approx. treatment. 

USE nlte_type
USE nlte_dim, only: id_ndept

implicit none

real(dp), intent(out) :: s1

real(dp), intent(in) :: xlambda, dx, taus, tauc, sl1, slinekblue
real(dp), intent(inout) :: iminusold

integer(i4b), intent(in) ::  nstep
integer(i4b), intent(inout) ::  ncount

logical, intent(inout) ::  in_line_old

! local variables
real(dp) :: iminus, iminus1, sl, diminus

if(taus.gt.0.05 .and. taus.gt.tauc) then
! dominating line processes, and line begins to become optically thick
sl=sl1 !sl1 is more or less constant over line; in this way, we avoid spikes
       !at the edges of the core (setting here slinek(1) is disastrous
iminus = sl*(1.d0 - exp(-taus-tauc)) ! including continuum tau

! for closely separated lines, iminus from previous line (outer resonance
! zones) might be larger than iminus from considered line core.
! In this case, we use the former value for iminus. 
! Exact treatmeat would yield iminus_exact = sl(1-exp(-tau)) + iminus1*exp(-tau)
  if(ncount.ne.0) then
! decaying iminus from previous line
    if(.not.in_line_old) stop ' calc_iminus: something rotten with in_line_old!'
    ncount=ncount+1
    iminus1=iminusold*(1.d0-1.d0/float(nstep-ncount+1))
    if (iminus1.gt.iminus) then
! approximate iminus_exact=sl(1-exp(-tau)) + iminus1*exp(-tau)
!             by iminus1 (rather small tau)
      iminus=iminus1
      sl=1.d99  ! hack to ensure that correct path chosen later on
    else
! approximate iminus_exact=sl(1-exp(-tau)) + iminus1*exp(-tau)
!             by iminus=sl(1-exp(-tau)) (large tau)
      ncount=0 
! end illumination from above resonance zones, since local term dominates
    endif
  else 
    in_line_old=.true.
    ncount=0  !counting begins when line core ends
  endif

else if(tauc.gt.1.) then
! dominating continuum process, but including line
sl=slinekblue !total source function
iminus = sl*(1.d0 - exp(-taus-tauc))
in_line_old=.true. ! also here, I-minus does not go down instantaneously, 
                       ! when tauc goes to zero, but account for the fact that
                       ! (optically thin) frequencies are illuminated from
                       ! bluewards intensities, which originate from higher levels
                       ! (Iminus(|mu|V,x) = Iminus((|mu|V+deltax,x+deltax)
ncount=0 !counting begins when continuum becomes optically thin

! here, we don't apply the upper correction for larger iminus1
! in case, include it also here

else
! optically thin continuum and line
    if (in_line_old) then
!   end of optically thick line core (or optically thick continuum)
!   encountered; account for outer resonance zones (continuum):
!   approximated in such a way that Iminus decreases linearly to
!   zero, where the spatial extent is calculated from the condition
!   that irradiation from above must be possible
!   a) in obsframe: outwards resonance zones possible from
!      xobs = -(|mu|V)_last ... -1
!   b) in cmf: from -xmax ... -1 + (|mu|V)_last
!      approximated by 0 ... -1 + |mu|_last = -1 +(z/r)_last 
!      THUS: additional NSTEPs = (1-z(1)/r(1))/DELTAX required, 
!      where Iminus decreases monotonically to zero (corresponding
!      to the fact that the resonance zones move outwards and SLINE
!      and OPAL (or SCONT and OPAC) decrease
      if (nstep.lt.2) then
!     for core rays, the irradiation ends abruptly
        sl=0.
        iminus = 0.
        in_line_old=.false.
        ncount=0.
      else

!     first step after end of line core (-xmax)
        if(ncount.eq.0) then
          ncount=1
          iminus=iminusold*(1.-1./float(nstep))
        else
!     next steps
          ncount=ncount+1
          iminus=iminusold*(1.d0-1.d0/float(nstep-ncount+1))
        endif

        if(iminus.eq.0.d0) then
! end of freq. range with outer res. zones encountered
          in_line_old=.false.
          ncount=0
! otherwise, leave in_line_old at true
        endif
        sl=1.d99  ! hack to ensure that correct path chosen later on
      endif

    else
! no irradiation from above
      sl=0.
      iminus = 0.
! consistency check
      if(ncount.ne.0) &
&       stop 'subroutine calc_iminus: no irradiance, but ncount ne 0'
! this path implies that in_line_old is already .false. 
    endif
endif

! now we have calculated iminus, can proceed to calculate dIminus/dx and s1

if (sl.eq.slinekblue) then
        diminus=iminus-sl
! Alternative formulation for a dominating opt. thick continuum
! diminus calculated according to 
! P/opa dIminus/dx approx Iminus - S -> effective source term is S.
! Here, the standard formulation gives spurious results in lines (below 228 A) 
! or too much absorption at blue edge
     else
        diminus=dx*(iminusold-iminus)
! Standard formulation. Inside line, diminus <= 0. Becomes positive 
! when the line core 'ends'. Then, and most important, the effective source
! term S1 becomes negative, and particularly the intensities for the
! non-core rays are damped. Without this term (e.g., alternative
! formulation from above) the intensities remain at the high values
! from previous S1 terms, and the redward fluxes become to high! 
! This condition also applies when continuum becomes opt. thin
     endif
!
! This is the final value, Iminus - P/kappa dIminus/dx
     s1=iminus-diminus

!for tests (include jp in header)
!if(jp.eq.1.and.xlambda.ge.1548.and.xlambda.lt.1555.) then
!  write(*,fmt='(f10.3,9(2x,e8.2))') xlambda,s1,iminus,diminus, &
!&   sl,taus,tauc
!endif

! for next frequency
     iminusold=iminus

return
end
!
!----------------------------------------------------------------------------
!
subroutine calc_tlu(nd,ilow,imax,xne,temp,clf,r,v,dvdr)

! ilow,imax nlte values

! calculates reduced RBB rate_coefficients (tlumat, tulmat) for explicit elements 
! indexcmf1 = 3,4,5,... : radiation field from cmf_complete
!                 (number of depacked components is indexcmf1-2
! indexcmf1 = 2 : radiation field from cmf_sing
! indexcmf1 = 1 : Sobolev approach

! TODO:
! attention for xlamc etc. (original from DATA or modified from xlamcmf_depacked)
! xlamc etc: attention if SA transport inside range
!
! IMPORTANT NOTE: all data regarding components as stored in met_rbb refer to
! original data (id,xlam,gf) and have been calculated in subr. fgrid_cmf  
!
! if treat_single=T, corresponding line (defined in upper part of routine)
! will be treated in single line approximation with average background
! in particular, this refers to the HeI resonance line in case opt_hei_sing = T  
  
USE nlte_type
USE fund_const

USE nlte_dim, only: id_ndept, id_atoms, id_nttrd

USE nlte_opt, only: opt_hei_sing

USE princesa_var, only: index1,indat1, data, le, ifirsl, iong, gl
USE nlte_var, only: mmlow, mmup, indexrbb, indexcmf, vdopcmf, vther, &
&             blucmf, bulcmf, aulcmf, sr, vmax, corrfc, &
&             lablinelo, lablineup, tulmat, tlumat, indexsa, betold, no_alo
USE cmf_multi, only: indov
USE nlte_app, only: teff
USE cmf_all_var, only: indexcmf1, occexpl, indexel_inv, nf_cmf, &
&                indexcmf1_lam, xj_cmf, alo_cmf, pweightdop_H, profdop_H, &
&                pweightdop_He, profdop_He, pweightdop, profdop, sumopal_cmf, &
&                sumetal_cmf, lam_cmf, lineno, depacked, &
&                indxlamc_depacked, xlamcmf_depacked, sumgf_depacked
USE ng_var, only: itng, index1_ng, sl_ng, ng
USE nlte_porvor, only: fic,tcl_fac_line,tcl_fac_cont
USE tcorr_var, only: temp_converged

implicit none
integer(i4b), parameter :: nd1=id_ndept, kel=id_atoms
real(dp), parameter :: c_tau=0.02654 !lambda in cgs

! .. arguments ..
integer(i4b), intent(in) :: nd
integer(i4b), dimension(nd1,kel), intent(in) :: ilow, imax

real(dp), dimension(nd1), intent(in) :: xne,temp,clf,r,v,dvdr

! .. local variables ..
integer(i4b), dimension(:), allocatable :: index_noc, iarr, indexx

real(dp), dimension(nd1) :: ajlmean, alo, vr, xmust

real(dp) :: er, xlamc, flu, blu, bul, aul, xnl, xnu, xlam, vdop, snew, &
&           tlu, tul, c_tau1, lam3, const, occlow, occup, opal, sline, &
&           jbar, alobar, rho1, alo1bar, sumgf, &
&           lami, fi, consti, opali, jbari, alobari, gfi, lam3i, sl, &
&           xlamcmf_dep, deltai, sumrho, const1, fip, lamp, ratio

integer(i4b) :: ii, indi, m, ml, mu, nato, indmat, indcmf, ll, ilo, ima, &
&               numion, igenio, j, k, ke, il, ixmax, ik, i, &
&               noc, inoc, instart, inend, in, iov, ij, il1, ikp

real(dp), dimension(0:60) :: pw1, pw2

integer(i4b), dimension(nd1,kel) :: nfir, nlas

logical :: inversion, treat_single, once

character*6 :: levlo, levup

!similar to routine cmfprep, but with fixed indices 
!assuming that ionization equilibrium does not change any longer
!otherwise, more complex approach required

! calculate nfir and nlas once, to allow for subsequent checks

if(maxval(nf_cmf).gt.60) stop ' dimension (pw1, pw2) > 60 in calc_tlu, adjust!' 

if(nd .ne. nd1) stop ' nd ne nd1 in calc_tlu'

do k=1,kel
  do ll=1,nd
    ilo = ilow(ll,k)  
    ima = imax(ll,k)  
    if (ima.lt.ilo) stop ' ima < ilo in calc_tlu'
    numion = igenio(k,ilo)  
    nfir(ll,k) = ifirsl(numion)  
    numion = igenio(k,ima) + 1  
    if (numion.gt.iong) stop 'error in calc_tlu - numion'  
    nlas(ll,k) = ifirsl(numion)  
  enddo
enddo

!prepara some aux variables
do ll=1,nd
! generalized blocking factors need to be updated
   call factinc(r(ll),ll)  
enddo
vr = v/r  
xmust = sqrt(1.d0-1.d0/r/r)  

c_tau1=c_tau*sr/vmax


! preparation of cmf- and SA-transport
print*
iiloop: do ii = 1,index1  
     indi = indat1(ii)  
     er = data(indi)  
     m = int(er)  
     if (m.ne.1) cycle  
     ml = mmlow(ii)  
     mu = mmup(ii)  

     nato=le(ml)

     xlamc = data(indi+1)*1.d-8  
     flu = data(indi+2)  

     indmat = indexrbb(ii)  
     indcmf = indexcmf1(indmat) !use new index

     levlo=lablinelo(indmat)  
     levup=lablineup(indmat)  
     xlamcmf_dep=xlamcmf_depacked(indmat)
     
     treat_single=.false.
!     if (levlo .eq. 'HE21  '.and. levup .eq. 'HE22  ') treat_single=.true.  
     if (opt_hei_sing) then
       if(indcmf.eq.3 .and. (levlo .eq. 'HE11S1  ' &
        .and. xlamcmf_dep.gt.584. .and. xlamcmf_dep.lt.585.)) treat_single=.true.
     endif
     
!    define index-array for depacked lines 
     sumgf=sumgf_depacked(indmat)
     if(sumgf.ne.0) then
       do i=1,lineno
          if(depacked(i)%levlo .eq. levlo .and. depacked(i)%levup .eq. levup) then
            noc=depacked(i)%nocomp
            allocate(index_noc(noc),iarr(noc),indexx(noc))
            inoc=1
            index_noc(inoc)=i
            instart=i
            exit
          endif  
       enddo
       if(noc.ne.1) then
         do i=instart+1,lineno
            if(depacked(i)%levlo .eq. levlo .and. depacked(i)%levup .eq. levup) then
              inoc=inoc+1
              index_noc(inoc)=i
            endif
         enddo   
       endif  
! test       
       if(inoc.ne.noc) stop 'something wrong with inoc (subr. calc_tlu)'          
! sort index_noc according to index (wavelength)
! most components in LINES_xxx.dat already sorted, but not necessarily
       do i=1,noc
         iarr(i)=depacked(index_noc(i))%index
       enddo
       call indexx_int(noc,iarr,indexx)
! aux. variable
       iarr=index_noc
       do i=1,noc
         index_noc(i)=iarr(indexx(i))
       enddo       

     endif !sumgf ne 0
!
!---     final check if everything was ok
!
     once=.true.
     do ll=1,nd

       if (ml.lt.nfir(ll,nato) .or. ml.ge.nlas(ll,nato)) then
         if(tlumat(indmat,ll).ne.0.) stop ' tlumat ne 0(1)' 
         cycle
       endif
       if (mu.gt.nlas(ll,nato)) then
         if(tlumat(indmat,ll).ne.0.) stop ' tlumat ne 0(2)' 
         cycle
       endif  
       
       xnl = occexpl(ml,ll)  
       xnu = occexpl(mu,ll)  
! all occupation numbers considered here should be defined
       if (xnl*xnu.eq.0.) stop ' occ. num not defined in calc_tlu'

       if (xlamcmf_dep.eq.0.d0) stop ' error in index gymnastics (calc_tlu)'

       if (indcmf.eq.1) then
! consistency check
          if (indexcmf(indmat).ne.1) stop ' indexcmf1=1 and indexcmf ne 1 in calc_tlu'
!          if (indexsa(indmat).eq.0) stop ' indexsa=0 for SA-lines in calc_tlu'
! in present philosophy (cmf_all implies GLOBAL=.TRUE.),
! indexsa should have a value of ND
! JO Nov. 2019 this might happen for cooler Temp. if He different ilow/imax as a function of depth 
! Meanwhile, hack in routine ILOWIMAX_LTE in nlte.f90
          if (indexsa(indmat).ne.nd) stop ' indexsa ne nd for SA-lines in calc_tlu'          
!
!---  Sobolev treatment, update of tulmat and tlumat
!---  NOTE: Sobolev lines (even inside range) should be never depacked.
!     Thus working with packed data
!     one more consistency check
          if(abs(1.d0-xlamcmf_dep*1.d-8/xlamc).gt.1.d-10) then
            print*,xlamcmf_dep,xlamc*1.d8
            stop ' problem with xlamc in SA-path (subr. calc_tlu)'
          endif
          betold(indmat,ll) = 0.d0  
          if (ll.eq.1) vdopcmf(indmat) = 0.d0  
          call tratrad(ll,ml,mu,flu,xlamc,xnl,xnu,xne(ll),gl,vr(ll), &
           dvdr(ll),tlu,tul,betold(indmat,ll),sr, vmax, xmust(ll),ii, &
           temp(ll),clf(ll),vther(nato), &
           fic(ll),tcl_fac_line(ll),tcl_fac_cont(ll))
! JO SA path needs to be checked when SA actually calculated
!    (inside range where tlu and tul not set to zero)
! Roughly checked, but note that I_inc etc. as used within range are now
! cmf-values with a lot of structure in frequency.
! Might be not what you want (which are smoothed values) 
          
!          if(ll.eq.1 .or. ll.eq.10 .or. ll.eq.30 .or. ll.eq. 40 .or. ll.eq.50) &
!&          print*,indmat,ll,tlu,tlumat(indmat,ll),tul,tulmat(indmat,ll)

          tlumat(indmat,ll) = tlu  
          tulmat(indmat,ll) = tul  
  
          if(once) write (*,fmt=9010) &
&           xlamcmf_dep,levlo,levup,indexsa(indmat)   
          once=.false.  
!       else if (indcmf.eq.2) then  
! if HeI resonance line should be treated in single line approx.
       else if (indcmf.eq.2  .or. treat_single) then 

! prepare single-line cmf-transport for lines outside range
! remember: tlumat and tulmat will be overwritten by opal and sline
! NOTE: outside range, all lines should be packed
! consistency check
          if (indexcmf(indmat).ne.2) stop ' indexcmf1=2 and indexcmf ne 2 in calc_tlu'
          if(abs(1.d0-xlamcmf_dep*1.d-8/xlamc).gt.1.d-10) then
            print*,xlamcmf_dep,xlamc*1.d8
            stop ' problem with xlamc in CMFSING-path (subr. calc_tlu)'
          endif
!           if(ll.eq.1 .or. ll.eq.10 .or. ll.eq.30 .or. ll.eq. 40 .or. ll.eq.50) &
!&          print*,indmat,ll,tlumat(indmat,ll),tulmat(indmat,ll)
!           print*,xlamc,indcmf,treat_single
          call tratcmf(ll,ml,mu,flu,xlamc,xnl,xnu,xne(ll),clf(ll),gl,sr,vmax,ii, &
                       fic(ll),tcl_fac_line(ll),tcl_fac_cont(ll))

! consistency check
! JO can be skipped when everything runs
          if (ll.eq.1) then  
            if (vdopcmf(indmat) .ne. vther(nato)) stop ' vdopcmf inconsistent in calc_tlu'  
            blu = e2mc2/hh*xlamc*4.d0*pi**2*flu
            bul = gl(ml)*blu/gl(mu)  
            aul = hc2/xlamc**3*bul
            if (blucmf(indmat) .ne. blu)  stop ' blu inconsistent in calc_tlu'
            if (bulcmf(indmat) .ne. bul)  stop ' bul inconsistent in calc_tlu'
            if (aulcmf(indmat) .ne. aul ) stop ' aul inconsistent in calc_tlu'
          end if  

!       else if (indcmf.eq.3) then
       else if (indcmf.ge.3  .and. .not.treat_single) then

! (re-)calculate line opacity and source-function
! This needs to be done as in subr. opacity_cmf, to ensure consistency of the ALO
!          
! consistency check
          if (indexcmf(indmat).ne.2) stop ' indexcmf1 ge 3 and indexcmf ne 2 in calc_tlu'
          lam3=xlamc**3 !in cgs
          const=c_tau1*xlamc*flu*gl(ml)
          occlow=xnl/gl(ml)
          occup =xnu/gl(mu)
! pi e2/me c * SR/vmax * lam * gf * (nl/gl - nu/gu) /clfac
          opal=const*(occlow-occup)/clf(ll)
          sline=hc2/lam3/(occlow/occup-1.d0)

! consistency check
! JO can be skipped when everything runs
          if (ll.eq.1) then  
            if (vdopcmf(indmat) .ne. vther(nato)) stop ' vdopcmf inconsistent in calc_tlu'  
            blu = e2mc2/hh*xlamc*4.d0*pi**2*flu
            bul = gl(ml)*blu/gl(mu)  
            aul = hc2/xlamc**3*bul
            if (blucmf(indmat) .ne. blu)  stop ' blu inconsistent in calc_tlu'
            if (bulcmf(indmat) .ne. bul)  stop ' bul inconsistent in calc_tlu'
            if (aulcmf(indmat) .ne. aul ) stop ' aul inconsistent in calc_tlu'
          end if  
! calculate reduced rates
! ALO for single line is int alo_nu rho1_nu phi dnu, with rho1_nu = opal * phi/opal_tot_nu
! note that alo_nu already includes rho_nu = chi_nu(line_tot)/chi_nu_tot

          ke=indexel_inv(nato)
          ixmax=nf_cmf(ke)
          il=indexcmf1_lam(indmat) !central frequency of strongest component
!test, can be skipped later on
          if(ixmax.eq.0) then
            print*,ll,nato,ke,ixmax,il
            stop 'ixmax = 0 in calc_tlu'
          endif  

! integration over phi
          jbar=0.
          alobar=0.
          
          if(sumgf.eq.0.) then
! path for packed components with original frequency or
! for lineno=0 (when only HHe or no depacked lines present) 
!
!JO Dec. 2016  this should also apply for lines with no components in
!         long line list when treated as explicit element
!         (corresponding to quality=2 when treated as bg element, see calc_tlu_met) 
!
! here we have to sum up jbar and alobar (one comp), and
! to set alo1bar = alobar at the end
!           if(ll.eq.1) write (*,fmt=9005) xlamc*1.d8,levlo,levup 
          if(ll.eq.1) write (*,fmt=9005) xlamcmf_dep,levlo,levup  
          
          if(opal.lt.0. .or. no_alo) then !don't use ALO
! alobar remains zero 
             if (ke.eq.1) then
               pw1(0:ixmax)=pweightdop_H(0:ixmax,ll)
             else if(ke.eq.2) then
               pw1(0:ixmax)=pweightdop_He(0:ixmax,ll)
             else
               pw1(0:ixmax)=pweightdop(0:ixmax,ll,ke)
             endif

! red side including center
             do ik=0,ixmax
               jbar=jbar+xj_cmf(ll,il+ik)*pw1(ik)
             enddo
! blue side
             do ik=1,ixmax
               jbar=jbar+xj_cmf(ll,il-ik)*pw1(ik)
             enddo
! standard path, no inversion, use ALO
          
          else  

             if (ke.eq.1) then
               pw1(0:ixmax)=pweightdop_H(0:ixmax,ll)
               pw2(0:ixmax)=profdop_H(0:ixmax,ll)
             else if(ke.eq.2) then
               pw1(0:ixmax)=pweightdop_He(0:ixmax,ll)
               pw2(0:ixmax)=profdop_He(0:ixmax,ll)
             else
               pw1(0:ixmax)=pweightdop(0:ixmax,ll,ke)
               pw2(0:ixmax)=profdop(0:ixmax,ll,ke)
             endif

! red side including center
             do ik=0,ixmax
               jbar=jbar+xj_cmf(ll,il+ik)*pw1(ik)
               rho1=opal*pw2(ik)/sumopal_cmf(ll,il+ik)
!  significant or dominating inversion of bg opacities
               if(rho1.gt.1.1d0 .or.  rho1.lt.0.d0) rho1=0.d0
!  very slight difference in lambda (SP vs. DP) or mild inversion of bg
               if(rho1.gt.1.d0) rho1=1.d0
               alobar=alobar+alo_cmf(ll,il+ik)*rho1*pw1(ik)
             enddo  

! blue side
             do ik=1,ixmax
               jbar=jbar+xj_cmf(ll,il-ik)*pw1(ik)
               rho1=opal*pw2(ik)/sumopal_cmf(ll,il-ik)
!  significant or dominating inversion of bg opacities
               if(rho1.gt.1.1d0 .or.  rho1.lt.0.d0) rho1=0.d0
!  very slight difference in lambda (SP vs. DP) or mild inversion of bg
               if(rho1.gt.1.d0) rho1=1.d0
               alobar=alobar+alo_cmf(ll,il-ik)*rho1*pw1(ik)
             enddo
!             if(levlo.eq.'N31' .and. levup.eq.'N311' .and. ll.eq.44) then
!               do ik=-100,100
!                 print*,lam_cmf(il+ik),sumopal_cmf(ll,il+ik)
!               enddo  
!             endif
  
             
          endif !inversion or not
          alo1bar=alobar ! in this case (packed transitions),
                         ! no difference between alobar and alo1bar
                         ! if no_alo, alobar=0, and thus alo1bar=0
          
          else
! path for depacked components or single component with modified frequency
          if(ke.eq.1 .or. ke.eq.2) stop ' H/He line in calc_tlu (depacked)'

! here we have to sum up jbari and alobari for individual components,
! and then to calculate jbar, alobar and alo1bar for the packed transition
! within the rate equations by weighting with gfi
          alo1bar=0.
          
          do in=1,noc

            i=index_noc(in) 
!            print*,i,xlamcmf_dep,depacked(i)%wave
            if(depacked(i)%levlo .ne. levlo .or. depacked(i)%levup .ne. levup) &
&             stop ' something wrong with index_noc'              
            lami=depacked(i)%wave
            fi=depacked(i)%flu

            if(ll.eq.1 .and. lami.eq. xlamcmf_dep) &
&              write (*,fmt=9006) xlamcmf_dep,levlo,levup  

              lami=lami*1.d-8
              consti=lami*fi/(xlamc*flu)
              opali=opal*consti
              il=depacked(i)%index
! only one component
              if(noc.eq.1) then
                if(abs(il-indexcmf1_lam(indmat)).gt.1) then
!under specific circumstances, there might be an index-difference of 1
!must not be more!
                  print*,ll,il,indexcmf1_lam(indmat)
                  stop ' error in index for lamcmf (ncomp=1)'
                endif  
                instart=in
                inend=in
              else
!overlapping components within profile; remember: index_noc wavelength-ordered
                do iov=in,1,-1
                  ij=index_noc(iov)
                  deltai=depacked(ij)%index-il
! for tests
                  if(deltai.gt.0) stop 'deltai > 0 (calc_tlu)'
                  if (abs(deltai).le.ixmax) instart=iov
                enddo  
                do iov=in,noc
                  ij=index_noc(iov)
                  deltai=depacked(ij)%index-il
! for tests:
                  if(deltai.lt.0) stop 'deltai < 0 (calc_tlu)'
                  if (deltai.le.ixmax) inend=iov
                enddo  
!                if(ll.eq.1) print*,levlo,levup,depacked(i)%lineno,instart,inend
              endif                
              jbari=0.
              alobari=0.
              
              if(opali.lt.0. .or. no_alo) then !don't use ALO
! alobari remains zero 
                pw1(0:ixmax)=pweightdop(0:ixmax,ll,ke)
 
! red side including center
                do ik=0,ixmax
                  jbari=jbari+xj_cmf(ll,il+ik)*pw1(ik)
                enddo
! blue side
                do ik=1,ixmax
                  jbari=jbari+xj_cmf(ll,il-ik)*pw1(ik)
                enddo
! standard path, no inversion, use ALO
          
              else  

                pw1(0:ixmax)=pweightdop(0:ixmax,ll,ke)
                pw2(0:ixmax)=profdop(0:ixmax,ll,ke)

! red side including center
                do ik=0,ixmax
                  jbari=jbari+xj_cmf(ll,il+ik)*pw1(ik)
                  if(instart.eq.inend) then
                    rho1=opali*pw2(ik)/sumopal_cmf(ll,il+ik)
                  else
                    sumrho=opali*pw2(ik)
                    const1=0.
                    do iov=instart,inend
                      if(iov.eq.in) cycle
                      ij=index_noc(iov)
                      il1=depacked(ij)%index
                      ikp=abs(ik+il-il1) ! other component's profile index
                      if(ikp.gt.ixmax) cycle ! outside other's component range
                      fip=depacked(ij)%flu
                      lamp=depacked(ij)%wave*1.d-8 ! in cgs
! correction for (slightly) different frequencies in SL:
! as calculated ALO refers to SL(iov);
! however, we want to apply it to SL(in). Thus: 
!                      const1=const1+(fip/fi)*(lamp/lami)*pw2(ikp)*(lami/lamp)**3
! shorter
                      const1=const1+(fip/fi)*pw2(ikp)*(lami/lamp)**2
                    enddo
                    sumrho=sumrho+const1*opali
                    rho1=sumrho/sumopal_cmf(ll,il+ik)
                  endif
                    
!  significant or dominating inversion of bg opacities
                  if(rho1.gt.1.1d0 .or.  rho1.lt.0.d0) rho1=0.d0
!  very slight difference in lambda (SP vs. DP) or mild inversion of bg
                  if(rho1.gt.1.d0) rho1=1.d0
                  alobari=alobari+alo_cmf(ll,il+ik)*rho1*pw1(ik)
                enddo  

! blue side
                do ik=1,ixmax
                  jbari=jbari+xj_cmf(ll,il-ik)*pw1(ik)
                  if(instart.eq.inend) then
                    rho1=opali*pw2(ik)/sumopal_cmf(ll,il-ik)
                  else
                    sumrho=opali*pw2(ik)
                    const1=0.
                    do iov=instart,inend
                      if(iov.eq.in) cycle
                      ij=index_noc(iov)
                      il1=depacked(ij)%index
                      ikp=abs(-ik+il-il1) ! other component's profile index
                      if(ikp.gt.ixmax) cycle ! outside other's component range
                      fip=depacked(ij)%flu
                      lamp=depacked(ij)%wave*1.d-8 ! in cgs
! correction for (slightly) different frequencies in SL:
! as calculated ALO refers to SL(iov);
! however, we want to apply it to SL(in). Thus: 
!                      const1=const1+(fip/fi)*(lamp/lami)*pw2(ikp)*(lami/lamp)**3
! shorter
                      const1=const1+(fip/fi)*pw2(ikp)*(lami/lamp)**2
                    enddo
                    sumrho=sumrho+const1*opali
                    rho1=sumrho/sumopal_cmf(ll,il-ik)
                  endif

!  significant or dominating inversion of bg opacities
                  if(rho1.gt.1.1d0 .or.  rho1.lt.0.d0) rho1=0.d0
!  very slight difference in lambda (SP vs. DP) or mild inversion of bg
                  if(rho1.gt.1.d0) rho1=1.d0
                  alobari=alobari+alo_cmf(ll,il-ik)*rho1*pw1(ik)
                enddo

              endif ! inversion or not

! sum up the contributions from individual components
              gfi=fi*gl(ml)
              lam3i=lam3/lami**3 !both quantities in cgs
              jbar=jbar+jbari*gfi
              alobar=alobar+alobari*gfi
              alo1bar=alo1bar+alobari*gfi*lam3i !includes correction for sline_i
!              print*,ll,depacked(i)%wave,depacked(i)%levlo,depacked(i)%levup  
!              print*,jbari,alobari,gfi
                  
          enddo ! loop over in (number of components)

!          if(levlo.eq.'N31' .and. levup.eq.'N311' .and. ll.eq.44) then
!               do ik=-100,100
!                 print*,lam_cmf(il+ik),sumopal_cmf(ll,il+ik)
!               enddo  
!          endif

! now, calculate actual jbar and alo by normalizing with sumgf and
! sumgflam3 = sum (gfi * lam^3/lam_i^3), respecively
!JO Dec 2016: sumgflam3 does not exist, and was not used in previous version
!of calc_tlu_met (approach changed from early on, and forgot to delete comment?)
          jbar=jbar/sumgf
          alobar=alobar/sumgf
          alo1bar=alo1bar/sumgf
!          print*,jbar,alobar,sumgf
!          print*
          
          endif ! packed or not 
! here, all required info should be present; 
! for only one transition (packed), alo1bar = alobar 
          if(itng.eq.2) then
! actually, jbar corresponds to J(n-1) corresponding itng=3,
! but itng is only updated (in subr. prep_ng) after present routine
! has been called. Thus, itng=2 still.
            ik=index1_ng(indmat)
            if(ik.ne.0) sl_ng(ll,ik,5)=jbar
          endif  

! in case, use Ng-extrapolated source-function.
! no rescaling of lambda required here (calculated with same xlamc) 
          if(ng) then
            ik=index1_ng(indmat)
            if(ik.ne.0) then
              sl=sl_ng(ll,ik,4)
              if(sl.gt.0.d0) sline = sl
            endif
          endif  
          
          if (abs(alobar).gt.1.d0) then
!              print*,' warning!!!! something wrong with line-alo in calc_tlu'
              stop ' warning!!!! something wrong with line-alo in calc_tlu'
              alobar = 0.d0  
          end if 
               
! here we have to consider the differences between packed und individual
! source function, thus working with alo1_tot
! the alo itself, however, is independent of source function and thus does
! not need to be corrected 
          snew = jbar - alo1bar*sline  

! for tests (no ALO)
!          snew = jbar
!          alobar = 0.
          
          if(levlo.eq.'N31' .and. levup.eq.'N311') then
             write(1111,*) 'N',ll,occlow,occup
!              print*,xlamc,ll,indmat
!             write(1111,*) jbar,alobar,alo1bar,sline,snew
          endif             
          
!           if(levlo .eq. 'HE11S1  '.and. levup .eq. 'HE12P1  ') &
!             print*,ll,jbar,sline,snew
!            if(levlo .eq. 'HE21    '.and. levup .eq. 'HE22    ') &
!             write(*,222) ll,jbar,alobar,sline,snew
!222 format('calcHE2',i3,4(2x,e12.6))           

           
          if (snew.lt.0.d0) then
!JO Dec. 2016: in case, modify according to calc_tlu_met
              print*,' warning!!!! something wrong with snew in calc_tlu'
              print*,xlamc,ll,indmat
              print*,jbar,alobar,sline,snew
              if(abs(snew).lt.sline/1000.) then
!                 snew=0.
!                 print*,'numerical problems, snew set to zero!'
              else
!                  stop ' warning!!!! something wrong with snew in calc_tlu'
                  print*,' warning!!!! something wrong with snew in calc_tlu'
                  snew=jbar
                  alobar=0.
              endif
          end if  

!hack to avoid oscillations in strongly coupled lines
          if (alobar.lt.0.1) then
             snew=jbar
             alobar=0.d0
          endif


          if(temp_converged) then
! this is important for OB-stars: HeII 303 needs to converge as reliable as possible,
! thus do not prohibit ALOs close to unity 
             if(alobar.gt.0.99 .and. &
                .not.(levlo .eq. 'HE21    '.and. levup .eq. 'HE22    ')) then
               ratio=alo1bar/alobar
               alobar=0.99d0 ! underestimate allowed
               alo1bar=ratio*0.99d0
               snew = jbar - alo1bar*sline  
             endif
          endif
          
!if(levup.eq.'N54   ' .and.ll.le.15) &
!& write(*,fmt='(A6,2x,i2,2x,e10.4,2x,f10.4,2x,2(e10.4,2x))') &
!& levlo,ll,jbar,alobar,sline,snew

! JO check whether this is also OK for inversion and alo=0
          tlu = blucmf(indmat)*snew  
          tul = bulcmf(indmat)*snew + aulcmf(indmat)* (1.d0-alobar)  
!if(xlamc.gt.955.33d-8 .and. xlamc.lt.955.4d-8) then
!  print*,xlamc,ll,jbar,alobar,sline,snew,tlu,tul
!endif


!           if(ll.eq.1 .or. ll.eq.10 .or. ll.eq.30 .or. ll.eq. 40 .or. ll.eq.50) &
!           write(*,fmt='(f10.3,3(2x,i4),2(2x,e8.2))') xlamc*1.d8,ii,indmat,ll, &
!&            1.-tlu/tlumat(indmat,ll),1.-tul/tulmat(indmat,ll)
          tlumat(indmat,ll) = tlu
          tulmat(indmat,ll) = tul 

       else
          stop 'error in indexcmf1(calc_tlu)'  

       end if  

     enddo  !depthloop  

     if(sumgf.ne.0.) deallocate(index_noc,iarr,indexx)

end do iiloop 
print*
             
!perform cmf-transport for lines outside range, update tulmat and tlumat
do j=1,id_nttrd  

    ii = indxlamc_depacked(j)  
    xlam = xlamcmf_depacked(ii)

     treat_single=.false.
!     if (lablinelo(ii) .eq. 'HE21  '.and. lablineup(ii) .eq. 'HE22  ') treat_single=.true.  
     if (opt_hei_sing) then
       if(indexcmf1(ii).eq.3 .and. (lablinelo(ii) .eq. 'HE11S1 ' &
          .and. xlam.gt.584. .and. xlam.lt.585.)) treat_single=.true.
     endif
    
!    if (indexcmf1(ii).eq.2) then
! if HeI should be treated in single line approx.
    if (indexcmf1(ii).eq.2 .or. (treat_single .and. indexcmf1(ii).ne.1)) then  

! similar as in routine cmf (nlte.f90)
! consistency check
        if (indexcmf(ii).ne.2) stop ' indexcmf1=2 and indexcmf ne 2 in calc_tlu'
        vdop = vdopcmf(ii)  
        if (vdop.eq.0.d0) stop ' error in vthcmf'  
! no overlap for lines outside overlap range
        if(indov(j).ne.0) stop ' indov ne 0 in calc_tlu. Set optsing = .true.' 

        call cmfsing &
&            (ii,r,v,dvdr,temp,xlam,vdop,vmax,teff,corrfc,ajlmean,alo,inversion)

        write (*,fmt=9000) xlam,lablinelo(ii),lablineup(ii)  

        do ll = 1,nd  
             snew = ajlmean(ll) - alo(ll)*tulmat(ii,ll)  
             if (snew.lt.0.d0 .or. abs(alo(ll)).gt.1.d0) then
               if(.not.inversion) print*,' warning!!!! something wrong with line-alo'
               if(inversion) then
                 print*,' warning!!!! jbar corrected at',ll
!JO at some point, this needs to be improved (also in nlte.f90)
                 ajlmean(ll)=0.                  
               endif  
               alo(ll) = 0.d0  
               snew = ajlmean(ll)  
             end if  

!             if (lablinelo(ii) .eq. 'HE21  '.and. lablineup(ii) .eq. 'HE22  ') &
!               write(*,222) ll,ajlmean(ll),alo(ll),tulmat(ii,ll),snew
             
!           if(lablinelo(ii) .eq. 'HE11S1  '.and.lablineup(ii) .eq. 'HE12P1  ') &
!             print*,ll,ajlmean(ll),tulmat(ii,ll),snew

             
             tlumat(ii,ll) = blucmf(ii)*snew  
             tulmat(ii,ll) = bulcmf(ii)*snew + aulcmf(ii)* (1.d0-alo(ll))  
!             if(ll.eq.1 .or. ll.eq.10 .or. ll.eq.30 .or. ll.eq. 40 .or. ll.eq.50) &
!&              print*,ii,ll,tlumat(ii,ll),tulmat(ii,ll)
         end do 

    endif
enddo

print *  
if(.not.no_alo) then
  print *,' RBB-RATES FOR EXPLICIT ELEMENTS PREPARED, ALO NE 0 (CALC_TLU)'  
else
  print *,' RBB-RATES FOR EXPLICIT ELEMENTS PREPARED, ALO = 0 (CALC_TLU)'  
endif
print *  

return

9000 format ('CMFT AT',F12.3,' A ',A6,' TO ',A6,' SINGLE LINE')  
9005 format ('CMFT AT',F12.3,' A ',A6,' TO ',A6,' MULTI LINE (PACKED)')  
9006 format ('CMFT AT',F12.3,' A ',A6,' TO ',A6,' MULTI LINE (DEPACKED, modified wavelengths)')  
9010 format ('SA-T AT',F12.3,' A ',A6,' TO ',A6,' FOR ',I3,' DEPTH', &
&       ' POINTS')
end
!
!----------------------------------------------------------------------------
!
subroutine calc_tlu_met(nd,clf)

! calculates reduced RBB rate_coefficients for selected bg elements 
! packing of individual components as follows (with slinei the
! source-function of the individual components and sline the source function
! of the packed transition: slinei=lam^3/lami^3*sline  

!jbartot = sum gfi*jbari / sum gfi (see thesis Jo)  

!=>

!tlu=snew = jbartot - alotot*sline = sum gfi*(jbari-aloi*slinei) / sum gfi =  
!jbartot - sum (aloi*lam^3/lam_i^3 sline) / sum gfi 

!=>

!snue=jbartot - alo1tot*sline    and  tul=alotot
!with jbartot as above, alo1tot = sum (gfi*aloi*lam^3/lami^3)/sum gfi, 
!     alotot= sum (gfi*alo)/sum gfi
  

USE nlte_type
USE nlte_dim, only: ID_NDEPT
USE fund_const, only: hc2
USE nlte_var, only: lwion1,sr,vmax,no_alo_metals
USE nlte_lines, only: unpacked_index_cmf, unpacked_gf,unpacked_xlam
USE nlte_app, only: met_imin, met_imax, met_rbb, natom, jatom_full, &
&                   indexel, indrec, irbb, ilin, occngold, tlu_met, tul_met
USE cmf_all_var, only: nf_cmf, elements_cmf, lam_cmf, pweightdop, &
&                   profdop, xj_cmf, alo_cmf, sumopal_cmf
USE ng_var, only: ng_met, itng_met, sl_ng_met
USE tcorr_var, only: temp_converged

implicit none

integer(i4b), intent(in) :: nd
real(dp), dimension(nd), intent(in) :: clf

real(dp), parameter :: c_tau=1.d-8*0.02654, hc2_8=hc2*1.d24 !conversion to Angstrom

real(dp), dimension(ID_NDEPT) :: opalfac, sline, jbar_tot, alo_tot, alo1_tot 
! required field-length only LWION1-1

real(dp), dimension(0:60) :: pw1, pw2

real(dp) :: srvmax, c_tau1, xlamc, lam3, occlowi, occupi, gfi, lam, const, &
&           opal, jbar, alobar, sumgf, snew, rho1, lam3i, ratio, dum, sl

real(dp) :: sumrho, const1, gfip, lamp

integer(i4b) :: k, j, ilow, imax, n, jj, ii, ji, iqual, no_comp, index_comp, &
&               in, lw1, ixmax, m1, m2, ll, il, ik

integer(i4b) :: iov, deltai, instart, inend, il1, ikp

logical, dimension(ID_NDEPT) :: calc, inv
logical :: test_oxy


! required field-length only LWION1-1

if(maxval(nf_cmf).gt.60) stop ' dimension (pw1, pw2) > 60 in calc_tlu_met, adjust!' 

srvmax=sr/vmax
c_tau1=c_tau*srvmax

lw1=lwion1-1

tlu_met=-9999.d0
tul_met=-9999.d0

do k=1,natom
 if (jatom_full(k).eq.0) cycle

! explicit elements will be treated by conventional (simplified) approach
! in netmat_met: since atomic data are different, particularly line- 
! frequencies, and the cmf-transport is performed with the 'explicit' data, 
! Jbars calculated using bg-data would be erroneous anyway. 
 if (indexel(k).gt.0) cycle
 
 if(k.lt.3) stop ' k < 3 in calc_tlu_met'

 ixmax=nf_cmf(k)

! if(ixmax.eq.0) stop 'ixmax = 0 in calc_tlu_met'
! change July2014. Consider only elements which are included in elements_cmf
 if (ixmax.eq.0) then
   if(elements_cmf(k).ne.0) &
     stop ' ixmax = 0 and elements_cmf(k) ne 0 in calc_tlu_met'
   cycle  
 endif
 
 ilow=minval(met_imin(k,1:lw1))
 imax=maxval(met_imax(k,1:lw1))

! we consider all ionization stages present between 1 to LWION1-1
do j=ilow,imax

   n=indrec(k,j)
   jj=irbb(n)-1

   do ii=1,ilin(n)
     ji=jj+ii
     iqual=met_rbb(ji)%quality
     no_comp=met_rbb(ji)%no_comp
     if (iqual.eq.2.and.no_comp.ne.0) stop ' iqual = 2 and no_comp ne 0'
! outside range or no individual component
! will be treated by conventional (simplified) approach in netmat_met
     if(iqual.eq.-1 .or.iqual.eq.2) cycle  
     
     m1=met_rbb(ji)%low  
     m2=met_rbb(ji)%lup
!     test_oxy=.false.
!     if(k.eq.8.and.j.eq.3.and.m1.eq.1.and.(m2.eq.10.or.m2.eq.5.or.m2.eq.6)) test_oxy=.true.
     index_comp=met_rbb(ji)%index_comp
     xlamc=met_rbb(ji)%wave !packed transition, in A
     lam3=xlamc**3
     
! depth dependent quantities which are identical for all components
     do ll=1,lw1

       occlowi=occngold(n,m1,ll) ! occng = n/g
       if(j.lt.met_imin(k,ll).or.j.gt.met_imax(k,ll).or.occlowi.eq.0.) then
! needn't to be considered or new ion
         calc(ll)=.false.
       else
         occupi=occngold(n,m2,ll)
! following two quantities are identical for all components
         opalfac(ll)=(occlowi-occupi)/clf(ll)  
         sline(ll)=hc2_8/lam3/(occlowi/occupi-1.d0)
! note: calculated with packed wavelength, corrected below
         if(opalfac(ll).le.0.) then
           inv(ll)=.true.
         else
           inv(ll)=.false.
         endif  
         calc(ll)=.true.
       endif
     enddo ! first depth loop 

     jbar_tot=0.
     alo_tot=0.
     alo1_tot=0.
     
! loop over all components
     do in=index_comp,index_comp+no_comp-1
! this is the depth-independent part

!        print*,k,j,ii,no_comp,'comp=',in
       il=unpacked_index_cmf(in)
       gfi=unpacked_gf(in)
!       lam=lam_cmf(il) ! in A (rounded)
       lam=unpacked_xlam(in) ! in A
! correction for difference in sline_i and sline
! sline_i = lam^3/lami^3*sline
       lam3i=lam3/lam**3 

!      overlapping components within profile
       do iov=in,index_comp,-1
         deltai=unpacked_index_cmf(iov)-il
! for tests:
         if(deltai.gt.0) stop 'deltai > 0 (calc_tlu_met)'
         if (abs(deltai).le.ixmax) instart=iov
       enddo  
       do iov=in,index_comp+no_comp-1
         deltai=unpacked_index_cmf(iov)-il
! for tests:
         if(deltai.lt.0) stop 'deltai < 0 (calc_tlu_met)'
         if (deltai.le.ixmax) inend=iov
       enddo  
! hack to avoid coupled ALOs
!       instart=in
!       inend=in
       

! for tests; note that there are certain ions that do not have
! overlapping lines according to line-list, e.g., PVI and SVII
!       if(instart.ne.inend) then
!         print*,k,j,ii,lam
!         print*,no_comp,'comp=',in,instart,inend
!       endif       
       
!       if(lam.gt.623.and.lam.lt.634.and.gfi.gt.0.05) print*,lam,k,j,m1,m2,gfi
           
       const=c_tau1*lam*gfi          
! this is the depth-dependent part       

       do ll=1,lw1
          if(.not.calc(ll)) cycle
          
! pi e2/me c * SR/vmax * lam * gf * (nl/gl - nu/gu) /clf
          opal=const*opalfac(ll)

! calculate reduced rates
! ALO for single line is int alo_nu rho1_nu phi dnu, with rho1_nu = opal * phi/opal_tot_nu
! note that alo_nu already includes rho_nu = chi_nu(line_tot)/chi_nu_tot

! integration over phi
          jbar=0.
          alobar=0.

          if(inv(ll)) then !inversion, don't use alo
! alobar remains zero 
            pw1(0:ixmax)=pweightdop(0:ixmax,ll,k)

! red side including center
            do ik=0,ixmax
              jbar=jbar+xj_cmf(ll,il+ik)*pw1(ik)
            enddo
! blue side
            do ik=1,ixmax
              jbar=jbar+xj_cmf(ll,il-ik)*pw1(ik)
            enddo
!            if(test_oxy) print*,'alobari_O',lam,jbar,alobar
             
! standard path, no inversion
          else  
            pw1(0:ixmax)=pweightdop(0:ixmax,ll,k)
            pw2(0:ixmax)=profdop(0:ixmax,ll,k)

! red side including center
            do ik=0,ixmax
              jbar=jbar+xj_cmf(ll,il+ik)*pw1(ik)
              dum=sumopal_cmf(ll,il+ik)              
              if(dum.eq.0.d0) then
! new ion, no sumopal has been calculated before
                rho1=0.d0
              else

                if(instart.eq.inend) then ! one component only
                  rho1=opal*pw2(ik)/dum
                else
                  sumrho=opal*pw2(ik)
                  const1=0.d0
                  do iov=instart,inend
                    if(iov.eq.in) cycle                     
                    il1=unpacked_index_cmf(iov)
                    ikp=abs(ik+il-il1) ! other component's profile index
                    if(ikp.gt.ixmax) cycle ! outside other's component range
                    gfip=unpacked_gf(iov)
                    lamp=unpacked_xlam(iov) ! in A
! correction for (slightly) different frequencies in SL:
! as calculated ALO refers to SL(iov);
! however, we want to apply it to SL(in). Thus: 
                    const1=const1+gfip*lamp*pw2(ikp)*(lam/lamp)**3
                  enddo
                  sumrho=sumrho+const1*c_tau1*opalfac(ll)
                  rho1=sumrho/dum
                endif

              endif
!              if(test_oxy) print*,lam,gfi,ik,rho1
!           if(k.eq.18.and.j.eq.5.and.m1.eq.2.and.m2.eq.10 .and. ll.ge.9 .and. ll.le.13) then
!             print*,ll,ik,il+ik,lam,xlamc,rho1,opal*pw2(ik),sumopal_cmf(ll,il+ik), &
!&              xj_cmf(ll,il+ik),alo_cmf(ll,il+ik),pw1(ik)              
!           endif

!  significant or dominating inversion of bg opacities
              if(rho1.gt.1.1d0 .or.  rho1.lt.0.d0) rho1=0.d0
!  slight difference in lambda (grid vs. exact) or mild inversion of bg
              if(rho1.gt.1.d0) rho1=1.d0
              alobar=alobar+alo_cmf(ll,il+ik)*rho1*pw1(ik)
            enddo  

! blue side
            do ik=1,ixmax
              jbar=jbar+xj_cmf(ll,il-ik)*pw1(ik)
              dum=sumopal_cmf(ll,il-ik)              
              if(dum.eq.0.d0) then
! new ion, no sumopal has been calculated before
                rho1=0.d0
              else

                if(instart.eq.inend) then ! one component only
                  rho1=opal*pw2(ik)/dum
                else
                  sumrho=opal*pw2(ik)
                  const1=0.d0
                  do iov=instart,inend
                    if(iov.eq.in) cycle                     
                    il1=unpacked_index_cmf(iov)
                    ikp=abs(-ik+il-il1) ! other component's profile index
                    if(ikp.gt.ixmax) cycle ! outside other's component range
                    gfip=unpacked_gf(iov)
                    lamp=unpacked_xlam(iov) ! in A
! correction for (slightly) different frequencies in SL:
! as calculated ALO refers to SL(iov);
! however, we want to apply it to SL(in). Thus: 
                    const1=const1+gfip*lamp*pw2(ikp)*(lam/lamp)**3
                  enddo
                  sumrho=sumrho+const1*c_tau1*opalfac(ll)
                  rho1=sumrho/dum
                endif
              endif
!              if(test_oxy) print*,lam,gfi,ik,rho1
!           if(k.eq.18.and.j.eq.5.and.m1.eq.2.and.m2.eq.10 .and. ll.ge.9 .and. ll.le.13) then
!             print*,ll,ik,il-ik,lam,xlamc,rho1,opal*pw2(ik),sumopal_cmf(ll,il-ik), &
!&              xj_cmf(ll,il-ik),alo_cmf(ll,il-ik),pw1(ik)
!           endif


!  significant or dominating inversion of bg opacities
              if(rho1.gt.1.1d0 .or.  rho1.lt.0.d0) rho1=0.d0
!  slight difference in lambda (grid vs. exact) or mild inversion of bg
              if(rho1.gt.1.d0) rho1=1.d0
              alobar=alobar+alo_cmf(ll,il-ik)*rho1*pw1(ik)
            enddo
!              if(test_oxy) print*,'alobari_O',lam,jbar,alobar

          endif !inversion or not

! sum up jbar and alo for all components
          jbar_tot(ll)=jbar_tot(ll)+jbar*gfi
          alo_tot(ll)=alo_tot(ll)+alobar*gfi
          alo1_tot(ll)=alo1_tot(ll)+alobar*gfi*lam3i !includes correction for sline_i
 
        enddo !depth-loop for component in

     enddo ! component loop for component in with index il

! now, calculate actual jbar and alo by normalizing with sumgf and
! sumgflam3 = sum (gfi * lam^3/lam_i^3), respecively
!JO Dec 2016: sumgflam3 does not exist, and was not used in previous version
!of calc_tlu_met (approach changed from early on, and forgot to delete comment?)
!
! sum(gfi) = sumgf already checked in subr. fgrid_cmf
       sumgf=met_rbb(ji)%sum_gf
! for all depth points
       jbar_tot=jbar_tot/sumgf
       alo_tot=alo_tot/sumgf
       alo1_tot=alo1_tot/sumgf

       if(itng_met.eq.2) then
! actually, jbar corresponds to J(n-1) corresponding itng_met=3,
! but itng_met is only updated (in subr. prep_ng_met) after present routine
! has been called. Thus, itng_met=2 still.
          ik=met_rbb(ji)%index1_ng
          if(ik.ne.0) sl_ng_met(:,ik,5)=jbar_tot(1:lw1)
       endif  

! in case, use Ng-extrapolated source-function.
! no rescaling of lambda required here (calculated with same xlamc) 
       if(ng_met) then
! here it might happen that ng_met and metals_converged, if latter
! was set to true in subr. Trad just before call to nlte_approx
! Thus, NO similar check as in subr. opacity_cmf
         ik=met_rbb(ji)%index1_ng
          if(ik.ne.0) then
              sl=sl_ng_met(35,ik,4)
              if(sl.gt.0.d0) sline(1:lw1) = sl_ng_met(:,ik,4)
          endif
       endif  

! now, perform some test and, in case, redefine alo!
       do ll=1,lw1
         if(.not.calc(ll)) then 
           if(jbar_tot(ll).ne.0. .or. alo_tot(ll).ne.0.) &
&            stop ' calc = .false. but jbar or alo ne 0 in calc_tlu_met'
           tlu_met(ji,ll)=0.
           tul_met(ji,ll)=0.

         else

           if (abs(alo_tot(ll)).gt.1.d0) then
!           print*,' warning!!!! something wrong with line-alo inq calc_tlu'
             stop ' warning!!!! something wrong with line-alo in calc_tlu_met'
             alo_tot(ll) = 0.d0  
           end if  

! here we have to consider the differences between packed und individual
! source function, thus working with alo1_tot
! the alo itself, however, is independent of source function and thus does
! not need to be corrected 
           snew = jbar_tot(ll) - alo1_tot(ll)*sline(ll)  

! for tests (no ALO)
!           snew = jbar_tot(ll)
!           alo_tot(ll) = 0.

           
!          if(test_oxy) then 
!              print*,xlamc,ll
!              print*,jbar_tot(ll),alo_tot(ll),alo1_tot(ll),sline(ll),snew
!          endif             
           

! final test on snew (with potential correction of alo)
           if (snew.lt.0.d0 .and. .not.no_alo_metals) then
             if(no_comp.eq.1 .or. alo_tot(ll).le.0.99d0) then 
               print*,lam,xlamc,k,j,m1,m2,no_comp
               print*,snew,jbar_tot(ll),alo1_tot(ll),sline(ll),alo_tot(ll)
!JO April 2017: stop statement removed: problems can happen right after
!                hydro update (d10v with ND=51)
               print*,'warning!!!! something wrong with snew in calc_tlu_met'
                  snew=jbar_tot(ll)
                  alo_tot(ll)=0.             
             else   
               print*,' warning!!!! problems with snew in calc_tlu_met'
               write(*,fmt='(f10.3,6(2x,i3))') lam,k,j,m1,m2,ll,no_comp
!try to correct (might happen if alo very close to one)
               ratio=alo1_tot(ll)/alo_tot(ll)
               alo_tot(ll)=0.99d0 ! underestimate allowed
               alo1_tot(ll)=ratio*0.99d0
               snew = jbar_tot(ll) - alo1_tot(ll)*sline(ll)  
               if (snew.lt.0.d0) then
!                 print*,xlamc,k,j,ii,m1,m2,no_comp,ll,inv(ll)
!                 print*,jbar_tot(ll),alo_tot(ll),alo1_tot(ll),sline(ll),snew
!JO Sept. 2019: stop statement removed; problems can occur in first
!               iterations of moments.
!                 stop ' problems cannot be corrected'
                  snew=jbar_tot(ll)
                  alo_tot(ll)=0.
               endif
             endif
           end if  

!for tests 
!           if (k.eq.8 .and. j.eq.3 .and. m1.eq.1 .and. m2.eq.10) then
!             n=indrec(8,3)
!             write(1111,*) 'O',ll,occngold(n,1,ll),occngold(n,10,ll)              write(1111,*) jbar,alobar,alo1bar,sline,snew
!             write(1111,*) jbar_tot(ll),alo_tot(ll),alo1_tot(ll),sline(ll),snew
!             print*,'O3/1-10',ll,sline(ll),snew,alo_tot(ll)
!           endif  


!hack to avoid oscillations in strongly coupled lines
!JO April 2018: changed from 0.8 to 0.9, since also for 0.9, there
!is a reasonable convergence without ALO (n propto 1/beta = 1/(1-ALO)),
!but we avoid certain convergence issues due to impact of non-diagonal
!elements which can be still substantial for ALO = 0.8
!(remember Ar V in the outer wind of zeta Pup)

           if (no_alo_metals) then
               snew=jbar_tot(ll)
               alo_tot(ll)=0.d0
           else if (alo_tot(ll).lt.0.9) then
!JO July 2018: Nevertheless, for O3 1-10 and N3 1-11 we keep the ALO 
             if (k.eq.8 .and. j.eq.3 .and. m1.eq.1 .and. m2.eq.10) then
               continue
             else if (k.eq.7 .and. j.eq.3 .and. m1.eq.1 .and. m2.eq.11) then
               continue
             else
               snew=jbar_tot(ll)
               alo_tot(ll)=0.d0
             endif
           endif
           
!JO June 2018; do not allow for too large ALOs. Can lead to oscillations,
!if slightly different from iteration to iteration, and ALO at 0.999...
!but use only if temperature already converged
           if(temp_converged) then
             if(alo_tot(ll).gt.0.99d0) then
               ratio=alo1_tot(ll)/alo_tot(ll)
               alo_tot(ll)=0.99d0 ! underestimate allowed
               alo1_tot(ll)=ratio*0.99d0
               snew = jbar_tot(ll) - alo1_tot(ll)*sline(ll)  
             endif
           endif

!for tests 
           if (k.eq.8 .and. j.eq.3 .and. m1.eq.1 .and. m2.eq.10) then
             n=indrec(8,3)
             write(1111,*) 'O',ll,occngold(n,1,ll),occngold(n,10,ll)
!             write(1111,*) jbar_tot(ll),alo_tot(ll),alo1_tot(ll),sline(ll),snew
!             print*,'O3/1-10',ll,sline(ll),snew,alo_tot(ll)
           endif  
           
!           if(alo_tot(ll).gt.0.95.and.instart.ne.inend) print*,k,j,ll,ii,alo_tot(ll),alo1_tot(ll)
           
!           tlu = blucmf(indmat)*snew  
!           tul = bulcmf(indmat)*snew + aulcmf(indmat)* (1.d0-alobar)  
           tlu_met(ji,ll)=snew
           tul_met(ji,ll)=alo_tot(ll)

!for tests
!           if(k.eq.8.and.j.eq.5.and.m1.eq.1.and.m2.eq.3 .and. ll.ge.17 .and. ll.le.25) then
!             print*,k,j,ii,m1,m2,no_comp,ll,inv(ll),jbar_tot(ll),alo_tot(ll),sline(ll),snew
!           endif
           
         endif! calc or .not. calc  

       enddo ! depth loop for packed line 

   enddo !line loop

 enddo! ion loop

enddo! atom loop

print *  
print *,' RBB-RATES FOR BACKGROUND ELEMENTS PREPARED (CALC_TLU_MET)'  
print *  

if(ng_met) then
! in case that ALMOST_CONVERGED was set to T just before
  ng_met=.false. 
  print *,' ng_met set to false'  
  print *  
endif
  
return

end
!
!----------------------------------------------------------------------------
!
subroutine interpol_cmf_to_coarse(nd)

! interpolates cmf-quantities (J,K) onto coarse grid, under the condition
! that integral values remain conserved.
! Note that integrations are performed over energy (Kayser) = 1/lambda 
  

USE nlte_type
USE nlte_dim, only: ID_NDEPT
USE fund_const, only: clight
USE nlte_var, only : ifre, fre, xj, xj_cmf_coarse, xxk, alo, modnam, &
                     kcmf_start, kcmf_end, no_check, xj_save
USE cmf_all_var, only: wavblue, wavred, lam_cmf, nftot, xj_cmf, optout_xj_xh

implicit none

integer(i4b), intent(in) :: nd

integer(i4b) :: i, nmax, is, nmax_ind, il, ig, k, n1, n2, ks
real(dp) :: eblue, ered, l1, l2, wnue1, weight, xn0, xn1, xn2, xn3, &
&           deh, de, errmax, eddftest

integer(i4b), allocatable, dimension(:) :: index
real(dp), allocatable, dimension(:) :: ener_cmf
real (dp), allocatable, dimension(:,:) :: integ_cmf, xint
!real (dp), allocatable, dimension(:,:), save :: xj_old

real(dp), dimension(id_ndept) :: di, di1, i0, i1, i2, i3, &
&         xiquart, xihalf, diff1, diff2

! interpol_lin included in contain

eblue=1.d8/wavblue
 ered=1.d8/wavred

!for tests
!eblue=1.d8/240.
! ered=1.d8/8264.

!note that fre is in rydberg, and ordered from red to blue
!the index-file and lam_cmf is ordered from blue to red

do i=1,ifre-1   
  if(fre(i).ge.ered) exit 
enddo

nmax=i
do i=ifre-1,nmax,-1   
  if(fre(i).le.eblue) exit 
enddo
is=i

if(abs(1.-fre(is)/eblue) .gt. 0.1 .or. abs(1.-fre(nmax)/ered) .gt. 0.1) then
  print*,wavblue,1.d8/fre(is),wavred,1.d8/fre(nmax)
  stop ' something wrong with freq. boundaries (subr. interpol_cmf_to_coarse)'  
endif

nmax_ind=is-nmax+1

allocate(index(0:nmax_ind-1))
allocate(integ_cmf(nd,nftot),ener_cmf(nftot))
allocate(xint(nd,ifre))

! xj_cmf_coarse should have been allocated in CMF_COMPLETE, under XJ_SMOOTHED
! in parallel, xj_save has been allocated.
if (.not. allocated(xj_cmf_coarse)) &
  stop ' xj_cmf_coarse not allocated in INTERPOL_CMF_TO_COARSE'

xj_cmf_coarse=0.d0 !overwrite previous values from XJ_SMOOTHED
ener_cmf=1.d8/lam_cmf

index(0)=is
index(1)=is
ig=2
  
do il=is-1,nmax,-1
      l1=fre(il+1)
      l2=fre(il)

      if(abs(1.-l2/l1) .lt. 1.d-5) then 
         if(il .eq. is-1) then 
           index(0)=is          	 
           index(1)=is-1
	   ig=2
         else
! if 3 frequency points very narrow, make new range
! e.g., instead of (479, 480),(480,481) -> (479,480),(481,482)
           if(index(ig-1) .ne. il+1) then
             index(ig)=il+1
             index(ig+1)=il
!             write(*,fmt='(i4,2(2x,f10.4),3(2x,i4))') &
!&             il,1.d8/l1,1.d8/l2,ig,index(ig),index(ig+1)
             ig=ig+2
	   endif
         endif      
       endif
enddo  

ig=ig-1
if (index(ig) .ne. nmax) then
  index(ig+1)=nmax   
  index(ig+2)=nmax
  ig=ig+2   
endif
   
if(modulo(ig,2) .ne. 1) then
  print*,'ig even in subroutine interpol_cmf_to_coarse'
  stop
endif  

!for tests
!do il=0,ig-1,2
!  write(*,fmt='(3(i4,2x),2(f10.4),2x)') &
!    il,index(il),index(il+1),1.d8/fre(index(il)),1.d8/fre(index(il+1))
!enddo

! preparation of energy integration weights (already checked in cmf_complete)
! integration weights from using trapezoidal rule over delta E, and 
! delta E = 1/lam*(1-1/(1.+v/c))
! factor 0.5 (from trapezoidal rule) * 1.d8 included in weights
! Basically
! vc1=lam_cmf(2)/lam_cmf(1)-1.d0
! wnue1=(1.d0-1.d0/(1.d0+vc1))*0.5d8
! leads to
wnue1=(1.d0-lam_cmf(1)/lam_cmf(2))*0.5d8

integ_cmf(:,1)=0.

!calculate integral over cmf-quantities (as a function of freq.)
do k=2,nftot

    weight=wnue1/lam_cmf(k-1) ! corresponds to 0.5*(nu(k-1)-nu(k))
    di=(xj_cmf(:,k-1)+xj_cmf(:,k))*weight
    integ_cmf(:,k)=integ_cmf(:,k-1)+di
!test for accuracy (truncation errors if J varies considerably with nu)
    di1=integ_cmf(:,k)-integ_cmf(:,k-1)
!JO July 2018 accuracy changed (1.d-5 to 1.d-3)
    if(maxval(abs(1.-di1/di)) .gt. 1.d-3) then
      print*,'problem with accuracy of J-integration'
      print*,i,lam_cmf(k)
      print*,di
      print*
      print*,di1
      print*
      print*,1.-di1/di
      stop ' problem with accuracy of J-integration'
    endif  

enddo

!now, approximate quantities for coarse grid, so that integrals
! remain conserved  
ks=1

do i=1,ig-2,2 
    n1=index(i)
    n2=index(i+1)
    
! piecewise, from edge to edge
    do k=n1,n2,-1

     if(k .eq. n1) then 
! standard approach and extrapolation (assuming linear function) at first point
      xn1=fre(k)
      xn2=0.5d0*(fre(k)+fre(k-1))
      xn3=fre(k-1)
      i1=interpol_lin(integ_cmf,ener_cmf,xn1,ks,nd,nftot)
      i2=interpol_lin(integ_cmf,ener_cmf,xn2,ks,nd,nftot)
      i3=interpol_lin(integ_cmf,ener_cmf,xn3,ks,nd,nftot)
      deh=xn1-xn2
      xiquart=(i2-i1)/deh
      xj_cmf_coarse(:,k)=xiquart ! standard
      de =xn1-xn3
      xihalf=(i3-i1)/de
      xint(:,k)=2.d0*xiquart-xihalf !linear extrapolation
    
     else if(k .eq. n2) then
! standard approach and extrapolation (assuming linear function) at last point
      xn0=fre(k+1) 
      xn1=xn2
      xn2=fre(k)
!restart of ks
      ks=1
      i0=interpol_lin(integ_cmf,ener_cmf,xn0,ks,nd,nftot)
      i1=i2 ! previously calculated
!      i1=interpol_lin(integ_cmf,ener_cmf,xn1,ks,nd,nftot)
      i2=interpol_lin(integ_cmf,ener_cmf,xn2,ks,nd,nftot)
      deh=xn1-xn2
      xiquart=(i2-i1)/deh
      xj_cmf_coarse(:,k)=xiquart ! standard
      de =xn0-xn2
      xihalf=(i2-i0)/de
      xint(:,k)=2.*xiquart-xihalf !linear extrapolation

     else
! almost exact at intermediate points
      xn1=xn2
      xn2=0.5*(fre(k)+fre(k-1))
      de=xn1-xn2
!      i1=interpol_lin(integ_cmf,ener_cmf,xn1,ks,nd,nftot)
      i1=i2 ! previously calculated
! ks continued
      i2=interpol_lin(integ_cmf,ener_cmf,xn2,ks,nd,nftot)
      xj_cmf_coarse(:,k)=(i2-i1)/de
!     xint(:,k)=xj_cmf_coarse(k) ; not needed
     endif
   enddo 
enddo


!--------
! check of accuracy: until here, total integrals should be identical
ered=fre(index(ig-1))
do i=nftot,1,-1   
  if(ener_cmf(i).gt.ered) exit 
enddo
nmax=i+1
!print*,1.d8/ered,lam_cmf(nmax)

eblue=fre(index(1))
do i=1,nftot   
  if(ener_cmf(i).lt.eblue) exit 
enddo
is=i-1
!print*,1.d8/eblue,lam_cmf(is)

! this is the energy integral over the cmf-value (slightly larger)
xiquart=(integ_cmf(:,nmax)-integ_cmf(:,is))

! now we calculate the energy integral over the exactly interpolated values,
! using the trapezoidal rule (wfre contains different weights, particularly at end)
xihalf=0.
do i=index(ig)+1,index(0)
  xihalf=xihalf+0.5d0*(xj_cmf_coarse(:,i)+xj_cmf_coarse(:,i-1))*(fre(i)-fre(i-1))
enddo

diff1=abs(1.-xihalf/xiquart)
errmax=maxval(diff1)

!do i=1,nd
!  print*,i,xiquart(i),xihalf(i)
!enddo  

print*
print*,' interpolation of J-CMF values onto coarse grid'
print*,' maximum deviation between total integrals = ',errmax  
print*

if(errmax.gt.1.d-4) stop ' total integrals (cmf vs coarse) before extrapolation not conserved'
!--------


! now we check which condition (standard or extrapol.) is better.
! standardwise, we use extrapol at the boundaries and standard else, 
! except for the case that diff(extrapol.) < diff(standard) over the edges
! this allows for a smooth transition for smooth intensities  
where(xint(:,index(1)).gt.0.)
  xj_cmf_coarse(:,index(1))=xint(:,index(1))
endwhere
where(xint(:,index(ig-1)).gt.0.)
  xj_cmf_coarse(:,index(ig-1))=xint(:,index(ig-1))
endwhere

do i=3,ig-2,2
    n1=index(i)
    diff1=abs(xint(:,n1)-xint(:,n1+1))
    diff2=abs(xj_cmf_coarse(:,n1)-xj_cmf_coarse(:,n1+1))
    where (diff1.lt.diff2 .and. xint(:,n1).gt.0. .and. xint(:,n1+1).gt.0.)
      xj_cmf_coarse(:,n1+1)=xint(:,n1+1)
      xj_cmf_coarse(:,n1)=xint(:,n1)
    endwhere
enddo
    
! specific treatment in case that first/last points are edges
! Check whether there are edges in XJ (Jnu(approx)).
! If so, modify xj_cmf_coarse(index(0)) and/or xj_cmf_coarse(index(ig)), 
! by using the corresponding values from xj.
! Otherwise, use interpolated values at neighbouring frequencies
if(index(0) .ne. index(1)) then
  diff1=abs(1.-xj(:,index(0))/xj(:,index(1)))
  where(diff1.gt.0.05)  
   xj_cmf_coarse(:,index(0))=xj(:,index(0))  
  elsewhere
   xj_cmf_coarse(:,index(0))=xj_cmf_coarse(:,index(1))  
  endwhere
endif

if(index(ig) .ne. index(ig-1)) then
  diff1=abs(1.-xj(:,index(ig))/xj(:,index(ig-1)))
  where(diff1.gt.0.05)  
   xj_cmf_coarse(:,index(ig))=xj(:,index(ig))  
  elsewhere
   xj_cmf_coarse(:,index(ig))=xj_cmf_coarse(:,index(ig-1))  
  endwhere
endif

!2nd check of accuracy, now for final values (not as exact, since linear extrapolation)
ered=fre(index(ig))
do i=nftot,1,-1   
  if(ener_cmf(i).gt.ered) exit 
enddo
nmax=i+1
!print*,1.d8/ered,lam_cmf(nmax)

eblue=fre(index(0))
do i=1,nftot   
  if(ener_cmf(i).lt.eblue) exit 
enddo
is=i-1
!print*,1.d8/eblue,lam_cmf(is)

! this is the energy integral over the cmf-value from eblue to ered (slightly larger)
xiquart=(integ_cmf(:,nmax)-integ_cmf(:,is))

! now we calculate the energy integral over the interpolated value,
! using the trapezoidal rule (wfre contains different weights, particularly at end)
xihalf=0.
do i=index(ig)+1,index(0)
  xihalf=xihalf+0.5d0*(xj_cmf_coarse(:,i)+xj_cmf_coarse(:,i-1))*(fre(i)-fre(i-1))
enddo

diff1=abs(1.-xihalf/xiquart)
errmax=maxval(diff1)

!do i=1,nd
!  print*,i,xiquart(i),xihalf(i)
!enddo  

print*,' maximum deviation between total integrals (incl. extrapolation) = ',errmax  
print*

!if(errmax.gt.0.01) stop ' too large deviation between total integrals (final)'
!JO changed Nov. 2019, to exclude problems in few specific cases
if(errmax.gt.0.05) stop ' too large deviation between total integrals (final)'

!finally, &
!we approximate XK in the cmf-region (from the coarse Eddington-factors), 
!use approx. values for XJ in those regions where lam < wavblue or lam > wavred, ...
!
!JO: at some point, we might use XK values from the cmf-solution
!... and overwrite present values (from CONT) with new ones (from CMF)

!for tests (if photo-rates shall be constant)
!if(.not.allocated(xj_old)) then
!  allocate(xj_old(nd,ifre))
!  xj_old=xj_cmf_coarse
!  print*,'xj_old allocated'
!else
!  xj_cmf_coarse=xj_old
!  print*,'xj_old used'
!endif  

!JO changed
!no where statement, since xj and xj_cmf_coarse of different dimensions:
!xj(nd,id_frec1),xj_cmf_coarse(nd,ifre)
kcmf_start=0
do k=1,ifre
  if (xj_cmf_coarse(1,k) .ne. 0.d0) then
  if(kcmf_start.eq.0) kcmf_start=k
  kcmf_end=k
!  print*,k,1.d8/fre(k),xxk(nd,k),xj(nd,k),xxk(nd,k)/xj(nd,k),xj_cmf_coarse(nd,k)
  xxk(:,k)=xxk(:,k)/xj(:,k)*xj_cmf_coarse(:,k) ! approximate and overwrite 
  alo(:,k)=0.d0 !reset alo in line region
  xj(:,k)=xj_cmf_coarse(:,k) ! overwrite
  endif
! JO March 2018
! test only for lambda > 20 (for higher freq. -- X-rays --), inconsistencies possible
  if(1.d8/fre(k).gt.20.) then
    eddftest=xxk(nd,k)/xj(nd,k)
!    print*,k,1.d8/fre(k),xxk(nd,k),xj(nd,k),eddftest
!    if(.not. no_check .and. (eddftest.lt.0.31 .or. eddftest.gt.0.35)) &
!&     stop ' inconsistent Edd. factor in interpol_cmf_to_coarse'
    if(.not. no_check .and. (eddftest.lt.0.31 .or. eddftest.gt.0.35)) then
     print*,' Warning!!! Warning!!! Warning!!!'
     print*,1.d8/fre(k),' ',eddftest
     print*,' inconsistent Edd. factor in interpol_cmf_to_coarse'
     print*
    endif   
  endif
enddo


!JO changed; xj_cmf_coarse only filled in CMF_COMPLETE range

xj(1,ifre+1)=2 !indicates cmf treatment

! needed for cmf_complete
do k=1,ifre
  xj_save(:,k)=xj(:,k)
enddo  

print*,' XJ and XK (coarse grid, from CONT) overwritten by CMF-values!,'
print*,' in the range ',1.d8/fre(kcmf_end),' to ',1.d8/fre(kcmf_start)
print*

if(1.d8/fre(kcmf_end+1).gt.wavblue) stop ' problems in cmf-range (blue)'
if(1.d8/fre(kcmf_start-1).lt.wavred) stop ' problems in cmf-range (red)'

if (optout_xj_xh) then
  open(1,file=trim(modnam)//'/out_xj_cmf_coarse')
    do k=1,ifre
      write(1,100) 1.d8/fre(k),xj(1,k),xj(43,k),xj(nd,k)
    enddo  
  close(1)
endif
  
!do k=1,ifre
!errmax=maxval(alo(:,k))
!  write(*,110) 1.d8/fre(k),xj(43,k),xj_cmf_coarse(43,k),alo(43,k),errmax
!enddo

!for tests
!do i=ifre,1,-1
!  write(*,fmt='(i3,2x,f20.4,2x,4(e12.6,2x))') i,1.d8/fre(i),xj(1,i),xj_cmf_coarse(1,i),xj(51,i),xj_cmf_coarse(51,i)
!enddo

deallocate(index,integ_cmf,ener_cmf,xint)

print*,' J/K-CMF interpolated onto coarse grid (integrals conserved),'
print*,' and previous values (from CONT) replaced' 
print*

return

100 format(f14.4,3(2x,e10.4))
110 format(f14.4,2(2x,e10.4),2(2x,f12.4))


contains

function interpol_lin(y,x,x1,ks,nd,nf)
! linear interpolation at x1 in y(:,i),y(:,i+1)
! abscissa values of y are x, located in descending order  

USE nlte_type

implicit none

real(dp), dimension(nd) :: interpol_lin

integer(i4b), intent(in) :: nd, nf
integer(i4b), intent(inout) :: ks

real(dp), dimension(nd,nf), intent(in) :: y
real(dp), dimension(nf), intent(in) :: x
real(dp), intent(in) :: x1

integer(i4b) :: i
real(dp) :: q, q1

do i=ks,nf-1
  if(x1 .le. x(i) .and. x1.gt. x(i+1)) goto 10
!  print*,i,x(i),x1,x(i+1)
enddo

stop ' x1 not found in x (interpol_lin)'

10 ks=i !for next index
q=(x1-x(i+1))/(x(i)-x(i+1))
q1=1.d0-q

interpol_lin=q*y(:,i)+q1*y(:,i+1)

return

end function

end
!
!----------------------------------------------------------------------------
!
subroutine interpol_hcmf_to_hcoarse(nd)

! interpolates h_obs (high resol) onto coarse grid, under the condition
! that integral values remain conserved.
! Note that integrations are performed over energy (Kayser) = 1/lambda 
! (programmed in analogy to INTERPOL_CMF_TO_COARSE  

! Note that lam_cmf and integ_cmf are actually obs. frame quantities 
USE nlte_type
USE nlte_dim, only: ID_NDEPT
USE fund_const, only: clight, sigsb, pi
USE nlte_app, only: teff
USE nlte_var, only : modnam, ifre, fre, xxh, xh_obs_coarse, &
  kcmf_start, kcmf_end, rtau23
USE cmf_all_var, only: wavblue, wavred, lam_cmf, nftot, h_obs, optout_xj_xh

implicit none

integer(i4b), intent(in) :: nd

integer(i4b) :: i, nmax, is, nmax_ind, il, ig, k, n1, n2, ks, &
                kcmf_start1, kcmf_end1
real(dp) :: eblue, ered, l1, l2, wnue1, weight, xn0, xn1, xn2, xn3, &
&           deh, de, fluxtot

integer(i4b), allocatable, dimension(:) :: index
real(dp), allocatable, dimension(:) :: ener_cmf
real (dp), allocatable, dimension(:) :: integ_cmf, xint, xh

real(dp) :: di, di1, i0, i1, i2, i3, &
&         xiquart, xihalf, diff1, diff2

! interpol_lin included in contain

eblue=1.d8/wavblue
 ered=1.d8/wavred

!for tests
!eblue=1.d8/240.
! ered=1.d8/8264.

!note that fre is in rydberg, and ordered from red to blue
!the index-file and lam_cmf is ordered from blue to red

do i=1,ifre-1   
  if(fre(i).ge.ered) exit 
enddo

nmax=i
do i=ifre-1,nmax,-1   
  if(fre(i).le.eblue) exit 
enddo
is=i

if(abs(1.-fre(is)/eblue) .gt. 0.1 .or. abs(1.-fre(nmax)/ered) .gt. 0.1) then
  print*,wavblue,1.d8/fre(is),wavred,1.d8/fre(nmax)
  stop ' something wrong with freq. boundaries (subr. interpol_hcmf_to_hcoarse)'  
endif

nmax_ind=is-nmax+1

allocate(index(0:nmax_ind-1))
allocate(integ_cmf(nftot),ener_cmf(nftot))
allocate(xint(ifre),xh(ifre))

! allocated only once, since routine only called once (UNASOL)
! usally allocated only one, since called only one (if UNASOL)
! but for tests, called more often. Thus
if(.not.allocated(xh_obs_coarse)) allocate(xh_obs_coarse(nftot))

xh=xxh(1,1:ifre)
! rescale xh to be consistent with h_obs
! xxh = hnu*(r(1)/sr)^2 
! hobs = hnu*(r(1)/srnom)^2 = hnu*rmax^2
! xh -> xxh * (sr/srnom)^2 = xxh/rtau23^2 
xh=xh/rtau23**2

xh_obs_coarse=0.d0
ener_cmf=1.d8/lam_cmf

index(0)=is
index(1)=is
ig=2
  
do il=is-1,nmax,-1
      l1=fre(il+1)
      l2=fre(il)

      if(abs(1.-l2/l1) .lt. 1.d-5) then 
         if(il .eq. is-1) then 
           index(0)=is          	 
           index(1)=is-1
	   ig=2
         else
! if 3 frequency points very narrow, make new range
! e.g., instead of (479, 480),(480,481) -> (479,480),(481,482)
           if(index(ig-1) .ne. il+1) then
             index(ig)=il+1
             index(ig+1)=il
!             write(*,fmt='(i4,2(2x,f10.4),3(2x,i4))') &
!&             il,1.d8/l1,1.d8/l2,ig,index(ig),index(ig+1)
             ig=ig+2
	   endif
         endif      
       endif
enddo  

ig=ig-1
if (index(ig) .ne. nmax) then
  index(ig+1)=nmax   
  index(ig+2)=nmax
  ig=ig+2   
endif
   
if(modulo(ig,2) .ne. 1) then
  print*,'ig even in subroutine interpol_hcmf_to_hcoarse'
  stop
endif  

!for tests
!do il=0,ig-1,2
!  write(*,fmt='(3(i4,2x),2(f10.4),2x)') &
!    il,index(il),index(il+1),1.d8/fre(index(il)),1.d8/fre(index(il+1))
!enddo

! preparation of energy integration weights (already checked in cmf_complete)
! integration weights from using trapezoidal rule over delta E, and 
! delta E = 1/lam*(1-1/(1.+v/c))
! factor 0.5 (from trapezoidal rule) * 1.d8 included in weights
! Basically
! vc1=lam_cmf(2)/lam_cmf(1)-1.d0
! wnue1=(1.d0-1.d0/(1.d0+vc1))*0.5d8
! leads to
wnue1=(1.d0-lam_cmf(1)/lam_cmf(2))*0.5d8

integ_cmf(1)=0.

!calculate integral over cmf-quantities (as a function of freq.)
do k=2,nftot

    weight=wnue1/lam_cmf(k-1) ! corresponds to 0.5*(nu(k-1)-nu(k))
    di=(h_obs(k-1)+h_obs(k))*weight
    integ_cmf(k)=integ_cmf(k-1)+di
!test for accuracy (truncation errors if H varies considerably with nu)
    di1=integ_cmf(k)-integ_cmf(k-1)
!JO July 2018 accuracy changed (1.d-5 to 1.d-3)
    if(abs(1.-di1/di) .gt. 1.d-3) then
      print*,'problem with accuracy of H-integration'
      print*,i,lam_cmf(k)
      print*,di
      print*,di1
      stop ' problem with accuracy of H-integration'
    endif  

enddo

!now, approximate quantities for coarse grid, so that integrals
!remain conserved  
ks=1

do i=1,ig-2,2 
    n1=index(i)
    n2=index(i+1)
    
! piecewise, from edge to edge
    do k=n1,n2,-1

     if(k .eq. n1) then 
! standard approach and extrapolation (assuming linear function) at first point
      xn1=fre(k)
      xn2=0.5d0*(fre(k)+fre(k-1))
      xn3=fre(k-1)
      i1=interpol_lin(integ_cmf,ener_cmf,xn1,ks,nftot)
      i2=interpol_lin(integ_cmf,ener_cmf,xn2,ks,nftot)
      i3=interpol_lin(integ_cmf,ener_cmf,xn3,ks,nftot)
      deh=xn1-xn2
      xiquart=(i2-i1)/deh
      xh_obs_coarse(k)=xiquart ! standard
      de =xn1-xn3
      xihalf=(i3-i1)/de
      xint(k)=2.d0*xiquart-xihalf !linear extrapolation
    
     else if(k .eq. n2) then
! standard approach and extrapolation (assuming linear function) at last point
      xn0=fre(k+1) 
      xn1=xn2
      xn2=fre(k)
!restart of ks
      ks=1
      i0=interpol_lin(integ_cmf,ener_cmf,xn0,ks,nftot)
      i1=i2 ! previously calculated
!      i1=interpol_lin(integ_cmf,ener_cmf,xn1,ks,nftot)
      i2=interpol_lin(integ_cmf,ener_cmf,xn2,ks,nftot)
      deh=xn1-xn2
      xiquart=(i2-i1)/deh
      xh_obs_coarse(k)=xiquart ! standard
      de =xn0-xn2
      xihalf=(i2-i0)/de
      xint(k)=2.*xiquart-xihalf !linear extrapolation

     else
! almost exact at intermediate points
      xn1=xn2
      xn2=0.5*(fre(k)+fre(k-1))
      de=xn1-xn2
!      i1=interpol_lin(integ_cmf,ener_cmf,xn1,ks,nftot)
      i1=i2 ! previously calculated
! ks continued
      i2=interpol_lin(integ_cmf,ener_cmf,xn2,ks,nftot)
      xh_obs_coarse(k)=(i2-i1)/de
     endif
   enddo 
enddo


!--------
! check of accuracy: until here, total integrals should be identical
ered=fre(index(ig-1))
do i=nftot,1,-1   
  if(ener_cmf(i).gt.ered) exit 
enddo
nmax=i+1
!print*,1.d8/ered,lam_cmf(nmax)

eblue=fre(index(1))
do i=1,nftot   
  if(ener_cmf(i).lt.eblue) exit 
enddo
is=i-1
!print*,1.d8/eblue,lam_cmf(is)

! this is the energy integral over the obs-value (slightly larger)
xiquart=(integ_cmf(nmax)-integ_cmf(is))

! now we calculate the energy integral over the exactly interpolated values,
! using the trapezoidal rule (wfre contains different weights, particularly at end)
xihalf=0.
do i=index(ig)+1,index(0)
  xihalf=xihalf+0.5d0*(xh_obs_coarse(i)+xh_obs_coarse(i-1))*(fre(i)-fre(i-1))
enddo

diff1=abs(1.-xihalf/xiquart)

!print*,xiquart,xihalf 

print*
print*,' interpolation of H-obs values onto coarse grid'
print*,' maximum deviation between total integrals = ',diff1 
print*

if(diff1.gt.1.d-4) stop ' total H-integrals (fine vs coarse) before extrapolation not conserved'
!--------


! now we check which condition (standard or extrapol.) is better.
! standardwise, we use extrapol at the boundaries and standard else, 
! except for the case that diff(extrapol.) < diff(standard) over the edges
! this allows for a smooth transition for smooth intensities  
if (xint(index(1)).gt.0.) xh_obs_coarse(index(1))=xint(index(1))
if (xint(index(ig-1)).gt.0.) xh_obs_coarse(index(ig-1))=xint(index(ig-1))

do i=3,ig-2,2
    n1=index(i)
    diff1=abs(xint(n1)-xint(n1+1))
    diff2=abs(xh_obs_coarse(n1)-xh_obs_coarse(n1+1))
    if (diff1.lt.diff2 .and. xint(n1).gt.0. .and. xint(n1+1).gt.0.) then
      xh_obs_coarse(n1+1)=xint(n1+1)
      xh_obs_coarse(n1)=xint(n1)
    endif
enddo
    
! specific treatment in case that first/last points are edges
! Check whether there are edges in XH (Hnu(approx)).
! If so, modify xh_obs_coarse(index(0)) and/or xh_cmf_coarse(index(ig)), 
! by using the corresponding values from XH.
! Otherwise, use interpolated values at neighbouring frequencies
if(index(0) .ne. index(1)) then
  diff1=abs(1.-xh(index(0))/xh(index(1)))
  if (diff1.gt.0.05) then  
   xh_obs_coarse(index(0))=xh(index(0))  
  else
   xh_obs_coarse(index(0))=xh_obs_coarse(index(1))  
  endif
endif

if(index(ig) .ne. index(ig-1)) then
  diff1=abs(1.-xh(index(ig))/xh(index(ig-1)))
  if (diff1.gt.0.05) then  
   xh_obs_coarse(index(ig))=xh(index(ig))  
  else
   xh_obs_coarse(index(ig))=xh_obs_coarse(index(ig-1))  
  endif
endif

!2nd check of accuracy, now for final values (not as exact, since linear extrapolation)
ered=fre(index(ig))
do i=nftot,1,-1   
  if(ener_cmf(i).gt.ered) exit 
enddo
nmax=i+1
!print*,1.d8/ered,lam_cmf(nmax)

eblue=fre(index(0))
do i=1,nftot   
  if(ener_cmf(i).lt.eblue) exit 
enddo
is=i-1
!print*,1.d8/eblue,lam_cmf(is)

! this is the energy integral over the cmf-value from eblue to ered (slightly larger)
xiquart=(integ_cmf(nmax)-integ_cmf(is))

! now we calculate the energy integral over the interpolated value,
! using the trapezoidal rule (wfre contains different weights, particularly at end)
xihalf=0.
do i=index(ig)+1,index(0)
  xihalf=xihalf+0.5d0*(xh_obs_coarse(i)+xh_obs_coarse(i-1))*(fre(i)-fre(i-1))
enddo

diff1=abs(1.-xihalf/xiquart)

!  print*,xiquart,xihalf

print*,' maximum deviation between total H-integrals (incl. extrapolation) = ',diff1  
print*

if(diff1.gt.0.01) stop ' too large deviation between total H-integrals (final)'

!now, add approximate fluxes in region outside CMF transfer 
kcmf_start1=0
do k=1,ifre
  if (xh_obs_coarse(k) .ne. 0.d0) then
    if(kcmf_start1.eq.0) kcmf_start1=k
    kcmf_end1=k
  else
    xh_obs_coarse(k)=xh(k) !remember, xh rescaled to (r(1)/srnom)^2
  endif
enddo

if(kcmf_start1.ne.kcmf_start) then
  print*,kcmf_start,kcmf_start1
  print*,kcmf_end,kcmf_end1
  stop ' problems in kcmf_start1'
endif
if(kcmf_end1.ne.kcmf_end) stop ' problems in kcmf_end1'


fluxtot=sigsb * teff**4

xihalf=0.
do i=ifre-1,1,-1
  xihalf=xihalf+0.5d0*(xh_obs_coarse(i)+xh_obs_coarse(i-1))*(fre(i)-fre(i-1))
enddo

xihalf=xihalf*4.d0*pi*clight

diff1=(1.-xihalf/fluxtot)
print*,' maximum deviation flux-integral at Rmax (coarse grid)'
print*,' and nominal flux = ',diff1, ' (positive: lower than nominal)'
print*

if (optout_xj_xh) then
  open(1,file=trim(modnam)//'/out_xh_obs_coarse')
    do k=1,ifre
      write(1,100) 1.d8/fre(k),xh_obs_coarse(k) 
    enddo  
  close(1)
endif

deallocate(index,integ_cmf,ener_cmf,xint,xh)

print*,' H-obs remapped onto coarse grid,'
print*,' in the range ',1.d8/fre(kcmf_end),' to ',1.d8/fre(kcmf_start)
print*

return

100 format(f14.4,3(2x,e10.4))


contains

function interpol_lin(y,x,x1,ks,nf)
! linear interpolation at x1 in y(i),y(i+1)
! abscissa values of y are x, located in descending order  

USE nlte_type

implicit none

real(dp) :: interpol_lin

integer(i4b), intent(in) :: nf
integer(i4b), intent(inout) :: ks

real(dp), dimension(nf), intent(in) :: y
real(dp), dimension(nf), intent(in) :: x
real(dp), intent(in) :: x1

integer(i4b) :: i
real(dp) :: q, q1

do i=ks,nf-1
  if(x1 .le. x(i) .and. x1.gt. x(i+1)) goto 10
!  print*,i,x(i),x1,x(i+1)
enddo

stop ' x1 not found in x (interpol_lin)'

10 ks=i !for next index
q=(x1-x(i+1))/(x(i)-x(i+1))
q1=1.d0-q

interpol_lin=q*y(i)+q1*y(i+1)

return

end function

end
!
!----------------------------------------------------------------------------
!
subroutine prep_convol

! prepares weights and corresponding array for approx. convolution with
! electron-scattering profile
  
USE nlte_type
USE nlte_dim, only: ID_NDEPT
USE fund_const, only: clight, akb
USE nlte_var, only : ifre, fre, xj
USE nlte_app, only: teff
USE cmf_all_var, only: ice, icearr
implicit none

real(dp), parameter :: emass = 9.1095d-28

real(dp) :: vthec, vthecm, e, de, ex, dx, summ

integer(i4b) :: k, kk, il, im, nmax

!assuming thermal broadening at Teff (to save space and time);
!treatment approximate anyway)
!broadening considered until 3*delta_nu
vthec=sqrt(2.d0*akb*teff/emass)/clight
vthecm=3.*vthec

allocate(ice(2,ifre))

ice(1,1)=1
ice(2,1)=1
ice(1,ifre)=ifre
ice(2,ifre)=ifre

do k=2,ifre-1
  e=fre(k)
  de=e*vthecm
  do kk=k-1,1,-1
    if(e-fre(kk).gt.de) then
      ice(1,k)=kk
      exit
    endif
  enddo
  do kk=k+1,ifre
    if(fre(kk)-e.gt.de) then
      ice(2,k)=kk
      exit
    endif
  enddo
!  print*,k,e,ice(1,k),ice(2,k)
enddo  

nmax=maxval(ice(2,:)-ice(1,:))+1

allocate(icearr(nmax,ifre))

!first and last weight unity 
icearr(1,1)=1.
icearr(1,ifre)=1.

! factor 0.5 and 1/(sqrt(pi*delta nu) left out, since renormalized anyway
! relative integration weights from trapez
do k=2,ifre-1
!do k=2,3
  e=fre(k)
  de=e*vthec
  il=ice(1,k)
  im=ice(2,k)

  if(im-il.lt.2) stop ' IM-IL < 2 in PREP_CONVOL'
  
  kk=il
  ex=(e-fre(kk))/de
  dx=fre(kk+1)-fre(kk)
  icearr(kk+1-il,k)=exp(-ex*ex)*dx
!  print*,il,ex,dx,icearr(kk+1-il,k)
  
  kk=im
  ex=(e-fre(kk))/de
  dx=fre(kk)-fre(kk-1)
  icearr(kk+1-il,k)=exp(-ex*ex)*dx
!  print*,im,ex,dx,icearr(kk+1-il,k)
  
  do kk=il+1,im-1
    ex=(e-fre(kk))/de
    dx=fre(kk+1)-fre(kk-1)
    icearr(kk+1-il,k)=exp(-ex*ex)*dx
!    print*,kk,ex,dx,icearr(kk+1-il,k)
  enddo  
!renormalization
  summ=0.
! note change in indices
  im=im-il+1
  il=1
  do kk=il,im
    summ=summ+icearr(kk,k)
  enddo
  summ=1./summ
  icearr(il:im,k)=icearr(il:im,k)*summ
enddo

end  
!
!----------------------------------------------------------------------------
!
subroutine xj_smooth(iopt)

!smooth xj for further use in Thomson emissivity
!if iopt=1, then xj (and xj_smoothed) will be overwritten;
!           additionally, xxk will be smoothed.
!           so far, only to be called in formal (nlte.f90)
!if iopt=2, then xj will be smoothed and xj_smoothed overwritten
!           so far, only to be called in cmf_complete (cmf_all.f90)
!if iopt=3, then xj_save will smoothed and xj_smoothed overwritten
!           so far, only to be called in cmf_complete (cmf_all.f90)

  
USE nlte_type
USE nlte_dim, only: ID_NDEPT
USE nlte_var, only : ifre, xj, xj_smoothed=>xj_cmf_coarse, xxk, xj_save, fre
USE cmf_all_var, only: ice, icearr
implicit none

integer(i4b), intent(in) :: iopt

integer(i4b) :: k,kk,il,im
real(dp), dimension(id_ndept) :: xj_aux

real(dp) :: eddftest

if(iopt.eq.1) then

! smooth xxk
  do k=1,ifre
    il=ice(1,k)
    im=ice(2,k)         
    xj_aux=0.d0
    do kk=il,im
      xj_aux=xj_aux+xxk(:,kk)*icearr(kk+1-il,k)
    enddo
    xj_smoothed(:,k)=xj_aux  !xj_smoothed used as auxiliary variable
  enddo
! overwrite xxk
  do k=1,ifre
    xxk(:,k)=xj_smoothed(:,k)
  enddo
  
! smooth xj  
  do k=1,ifre
    il=ice(1,k)
    im=ice(2,k)         
    xj_aux=0.d0
    do kk=il,im
      xj_aux=xj_aux+xj(:,kk)*icearr(kk+1-il,k)
    enddo
    xj_smoothed(:,k)=xj_aux  !xj_smoothed used as auxiliary variable
  enddo

! overwrite xj
  do k=1,ifre
    xj(:,k)=xj_smoothed(:,k)
    eddftest=xxk(id_ndept,k)/xj(id_ndept,k)
!JO March 2018
    if (1.d8/fre(k).gt. 20. .and. (eddftest.lt.0.31 .or. eddftest .gt. 0.35)) then
     print*,' Warning!!! Warning!!! Warning!!!'
     print*,1.d8/fre(k),' ',eddftest
     print*,' inconsistent Edd. factor in xj_smooth'
     print*
    endif   
  enddo
  
else if(iopt.eq.2) then

! smooth only xj, and do not overwrite  
  do k=1,ifre
    il=ice(1,k)
    im=ice(2,k)         
    xj_aux=0.d0
    do kk=il,im
      xj_aux=xj_aux+xj(:,kk)*icearr(kk+1-il,k)
    enddo
    xj_smoothed(:,k)=xj_aux  ! overwrite xj_smoothed
  enddo

else if(iopt.eq.3) then

! smooth only xj_save, and do not overwrite  
  do k=1,ifre
    il=ice(1,k)
    im=ice(2,k)         
    xj_aux=0.d0
    do kk=il,im
      xj_aux=xj_aux+xj_save(:,kk)*icearr(kk+1-il,k)
    enddo
    xj_smoothed(:,k)=xj_aux  ! overwrite xj_smoothed
  enddo

else 
    print*,iopt      
    stop ' wrong iopt in xj_smooth'      
endif

return
end          
!
!----------------------------------------------------------------------------
!
subroutine save_line_list
USE nlte_type
USE nlte_var, ONLY: modnam,lwion1
USE nlte_app, ONLY: lwion,indrec
USE cmf_all_var, ONLY: ntot,id1,xlam1,gf1

open(1,file=trim(modnam)//'/LINE_LIST_MERGED',form='UNFORMATTED')
write(1) ntot,lwion,lwion1,indrec
write(1) id1(1:ntot),xlam1(1:ntot),gf1(1:ntot)
close(1)

!print*,ntot
!print*,xlam1(1),xlam1(ntot)

return
end
!
!-----------------------------------------------------------------------
!
subroutine chih_params1(chi,p,t,delchi,sigemin,const,expo,nd,ns,nsig)
!
! parameterizes chibar_h, in case by optimizing sigemin
!
use nlte_type
use nlte_dim, ONLY: id_ndept
implicit none

integer(i4b), parameter:: nd1=id_ndept
integer(i4b), intent(in) :: nd, ns, nsig
real(dp), dimension(nd), intent(in) :: chi,p,t
real(dp), dimension(nd), intent(out) :: delchi

real(dp), intent(inout)  :: sigemin !different from original version
real(dp), intent(out) :: const,expo

real(dp), dimension(nd1) :: fit,chi1,t1,p1

integer(i4b) :: k,imin,i,l
real(dp) :: chi2min, expoi, consti, corr, chi2, maxchi, sigemini

!at first, old method (with updated calculation of CHI2), i.e, with given SIGEMIN
chi2min=1.d5
expo=0.

! test for best interval to perform the fit

do imin=nsig,nd-4
! regression from imin on
  call linreg(log10(t(imin:nd)),log10((chi(imin:nd)-sigemin)/p(imin:nd)), &
&   nd+1-imin,expoi,consti,corr)

! fit quality for "all" points (from nsig on)  
  fit(nsig:nd)=(10.**consti *t(nsig:nd)**expoi *p(nsig:nd))+sigemin
  chi2=sum((fit(nsig:nd)/chi(nsig:nd)-1.d0)**2)
  if (chi2 .lt. chi2min) then
    chi2min=chi2
    expo=expoi
    const=consti
  endif  
enddo

!if (expo.eq.0.) stop ' regression not successful (chih_params1#1)'

const=10.**const
! for all photospheric points (from NS on)
fit(ns:nd)=(const *t(ns:nd)**expo *p(ns:nd))+sigemin
delchi(ns:nd)=chi(ns:nd)/fit(ns:nd)
!do i=ns,nd
!  print*,i,' ',chi(i),' ',fit(i),' ',delchi(i) 
!enddo

if(maxval(delchi(ns:nd)).le.2. .and. minval(delchi(ns:nd)).ge.0.2) then
  print*
  print*,' routine CHIH_PARAMS1: method 1 used'
  print*,' Min/Max(DELCHI) = ',minval(delchi(ns:nd)),maxval(delchi(ns:nd))
  print*
  return
endif
  
! if DELCHI problematic (e.g., for first models of d10v), optimize sigemin

chi2min=1.d5
expo=0.
sigemini=sigemin
sigemini=int(sigemini*10.)/10.d0 !round SIGEMIN
maxchi=maxval(chi(nsig:nd))

do while (sigemini.lt.maxchi) 

loop: do imin=nsig,nd-4
! regression from imin on
  i=0
  do l=imin,nd
    if((chi(l)-sigemini).gt.0.) then
      i=i+1
      chi1(i)=chi(l)
      t1(i)=t(l)
      p1(i)=p(l)
    endif  
  enddo
  if(i.lt.2) cycle loop
  
  call linreg(log10(t1(1:i)),log10((chi1(1:i)-sigemini)/p1(1:i)), &
&              i,expoi,consti,corr)
  
! fit quality for "all" points (from nsig on)  
  fit(nsig:nd)=(10.**consti *t(nsig:nd)**expoi *p(nsig:nd))+sigemini
! though we finally calculate chi/fit, we optimize w.r.t. fit/chi, to allow
! for larger errors if chi is small.  
!  chi2=sum((fit(nsig:nd)/chi(nsig:nd)-1.d0)**2)
  chi2=sum((fit(nsig:nd)/chi(nsig:nd)-1.d0)**2)
!  print*,sigemini,imin,chi2
  if (chi2 .lt. chi2min) then
    chi2min=chi2
    expo=expoi
    const=consti
    sigemin=sigemini
  endif  
enddo loop

sigemini=sigemini+0.1d0

enddo

if (expo.eq.0.) stop ' regression not successful (chih_params1#2)'

print*,chi2min,expo,const,sigemin
const=10.**const
! for all photospheric points (from NS on)
fit(ns:nd)=(const *t(ns:nd)**expo *p(ns:nd))+sigemin
delchi(ns:nd)=chi(ns:nd)/fit(ns:nd)
!do i=ns,nd
!  print*,i,' ',chi(i),' ',fit(i),' ',delchi(i) 
!enddo

if(maxval(delchi(ns:nd)).gt.2. .or. minval(delchi(ns:nd)).lt.0.2) then
  do i=ns,nd
    print*,i,' ',chi(i),' ',fit(i),' ',delchi(i) 
  enddo
  stop ' no regression possible (chih_params1)! re-try with update_struct = .false.'
endif

print*
print*,' routine CHIH_PARAMS1: method 2 used'
print*,' Min/Max(DELCHI) = ',minval(delchi(ns:nd)),maxval(delchi(ns:nd))
print*

return
end
!
!-----------------------------------------------------------------------
!
subroutine overlap_detector(dvdr)
!
! detects strong overlapping lines for lines from explicit and selected elements 
! note: uses reduced line-list
!  
use nlte_type
use fund_const, ONLY: clight
use princesa_var, ONLY: zl,labl
use nlte_app, ONLY: jatom_full
use cmf_all_var, ONLY: id1,xlam1,opal_ntest,vturb_cmf,icfirst,iclast,indexel_inv, &
  index_id1,gf1

implicit none

real(dp), parameter :: gf_cut = 1.d-3
real(dp), parameter :: opal_cut = 0.1
real(dp), parameter :: tau_cut = 0.1
real(dp), parameter :: fac_vturb = 1.

real(dp), intent(in) :: dvdr

integer(i4b) :: ii,iexp,isel,isel1,index,kk,irest,k,j,low,up, &
                i,indexi,ki,ji,lowi,upi,nsize,nsize1,in, &
                kold,jold,lowold

real(sp) :: lam, opalii, lami

real(dp) :: fac, red, blue, dxlam, lamm, lamp, tauii

logical :: expl, new, expli

nsize=size(id1)
nsize1=size(index_id1)
if (nsize.ne.nsize1) stop ' error in nsize (overlap_detector)'
! sort according to index;
! after sorting, selected elements first, explicit elements at the end
call indexx_int(nsize,id1,index_id1)

fac=fac_vturb*vturb_cmf/clight
blue=xlam1(icfirst)
red=xlam1(iclast)
print* 
print*,' overlap-detector between',blue,red

iexp=0
isel=0
isel1=0

kold=0
jold=0
lowold=0

iiloop: do in=1,nsize
!sorted according to element, ionization, lower and upper level
   ii=index_id1(in)
   if(gf1(ii).lt.gf_cut) cycle iiloop
   index=id1(ii)
   if(index.eq.0) cycle iiloop
   lam=xlam1(ii)
   if(lam.lt.blue .or. lam.gt.red) cycle iiloop
   dxlam=fac*lam
   lamm=lam-dxlam
   lamp=lam+dxlam
   if(index.ge.1d8) then
! explicit element
      kk = index/1d8
      irest = index - kk*1d8
      low = irest/10000
      up = irest - low*10000
      k=indexel_inv(kk)
      j=int(zl(low)+1.d-6)+1
      opalii=opal_ntest(ii)
      tauii=opalii/dvdr
      if(tauii.lt.tau_cut) cycle iiloop
      if(opalii.eq.0.) stop ' opalii = 0 for explicit element'
      expl=.true.
      iexp=iexp+1
   else
      k = index/1000000
      if(jatom_full(k).ne.1) cycle iiloop
      irest = index - k*1000000
      j = irest/100000
      irest = irest - j*100000
      low = irest/100
      up  = irest-low*100
      if(up.eq.0) cycle iiloop
      opalii=opal_ntest(ii)      
      tauii=opalii/dvdr
      if(tauii.lt.tau_cut) cycle iiloop
      isel1=isel1+1
      if(opalii.eq.0.) cycle
      isel=isel+1
      expl=.false.
   endif    

   new=.true.

iloop1: do i=ii-1,icfirst,-1
      lami=xlam1(i)
      if(lami.lt.lamm) exit iloop1
      if(opal_ntest(i).lt.opal_cut*opalii) cycle iloop1 
      indexi=id1(i)
      if(indexi.ge.1d8) then
! explicit element
         kk = indexi/1d8
         irest = indexi - kk*1d8
         lowi = irest/10000
         upi = irest - lowi*10000
         ki=indexel_inv(kk)
         ji=int(zl(lowi)+1.d-6)+1
         expli=.true.
      else
         ki = indexi/1000000
         irest = indexi - ki*1000000
         ji = irest/100000
         irest = irest - ji*100000
         lowi = irest/100
         upi  = irest-lowi*100
         expli=.false.
      endif

      if (k.ne.kold .or. j.ne.jold .or. low.ne.lowold) then
        print*
        kold=k
        jold=j
        lowold=low
      endif

      if(new) then
        if(expl) then          
          write(*,9000) k,j,low,up,lam,adjustr(labl(low)),adjustr(labl(up))      
        else
          write(*,9001) k,j,low,up,lam     
        endif
        new=.false.
      endif
      if (expli) then
         write(*,9002) ki,ji,lowi,upi,lami,opal_ntest(i)/opalii,adjustr(labl(lowi)),adjustr(labl(upi))      
      else
         write(*,9003) ki,ji,lowi,upi,lami,opal_ntest(i)/opalii
      endif
      
   enddo iloop1 

iloop2: do i=ii+1,iclast
      lami=xlam1(i)
      if(lami.gt.lamp) exit iloop2
      if(opal_ntest(i).lt.opal_cut*opalii) cycle iloop2 
      indexi=id1(i)
      if(indexi.ge.1d8) then
! explicit element
         kk = indexi/1d8
         irest = indexi - kk*1d8
         lowi = irest/10000
         upi = irest - lowi*10000
         ki=indexel_inv(kk)
         ji=int(zl(lowi)+1.d-6)+1
         expli=.true.
      else
         ki = indexi/1000000
         irest = indexi - ki*1000000
         ji = irest/100000
         irest = irest - ji*100000
         lowi = irest/100
         upi  = irest-lowi*100
         expli=.false.
      endif    

      if (k.ne.kold .or. j.ne.jold .or. low.ne.lowold) then
        print*
        kold=k
        jold=j
        lowold=low
      endif

      if(new) then
        if(expl) then          
          write(*,9000) k,j,low,up,lam,adjustr(labl(low)),adjustr(labl(up))      
        else
          write(*,9001) k,j,low,up,lam     
        endif
        new=.false.
      endif
      if (expli) then
         write(*,9002) ki,ji,lowi,upi,lami,opal_ntest(i)/opalii,adjustr(labl(lowi)),adjustr(labl(upi))      
      else
         write(*,9003) ki,ji,lowi,upi,lami,opal_ntest(i)/opalii
      endif
            
   enddo iloop2 

enddo iiloop

print*
print*,' overlap interval +/- ',fac_vturb*vturb_cmf/1.d5,' km/s'
print*,' all lines until gf > ',gf_cut
print*,' all lines until tau_Sob > ',tau_cut
print*,' all overlaps until a ratio > ',opal_cut,' (overlapping to considered line)'
print*
print*,' number of lines from explicit elements ',iexp
print*,' number of lines from selected elements (up ne 0) ',isel
print*,' number of lines from selected elements (up ne 0, opal ne 0) ',isel1
print*
return

9000 format(i2,2x,i2,2x,i5,2x,i5,2x,f12.4,2x,a6,2x,a6)
9001 format(i2,2x,i2,2x,i5,2x,i5,2x,f12.4)

9002 format('overlap with: ',i2,2x,i2,2x,i5,2x,i5,2x,f12.4,2x,e10.4,2x,a6,2x,a6)
9003 format('overlap with: ',i2,2x,i2,2x,i5,2x,i5,2x,f12.4,2x,e10.4)


end
!
!-----------------------------------------------------------------------
!
subroutine indexx_int(n,arr,indx)  
!
! as subr. indexx (nlte.f90), but for integer array
!  
USE nlte_type
USE nlte_dim
implicit none
!
!     .. parameters ..
integer(i4b), parameter :: m=7,nstack=50  
!     ..
!     .. scalar arguments ..
integer(i4b) ::  n  
!     ..
!     .. array arguments ..
integer(i4b) ::  arr(n), indx(n)  
!     ..
!     .. local scalars ..
!JO changed March 2017
!real(dp) ::  a  
integer(i4b) :: a
integer(i4b) ::  i,indxt,ir,itemp,j,jstack,k,l  
!     ..
!     .. local arrays ..
integer(i4b) ::  istack(nstack)  
!     ..

do j = 1,n  
     indx(j) = j  
end do  

jstack = 0  
l = 1  
ir = n  

   20 continue  

if (ir-l.lt.m) then  

jloop: do j = l + 1,ir  
          indxt = indx(j)  
          a = arr(indxt)  
          do i = j - 1,1,-1  
               if (arr(indx(i)).le.a) go to 40  
               indx(i+1) = indx(i)  
          end do
          i = 0  

   40           continue  

          indx(i+1) = indxt  
     end do jloop

     if (jstack.eq.0) return  
     ir = istack(jstack)  
     l = istack(jstack-1)  
     jstack = jstack - 2  

else  

     k = (l+ir)/2  
     itemp = indx(k)  
     indx(k) = indx(l+1)  
     indx(l+1) = itemp  

     if (arr(indx(l+1)).gt.arr(indx(ir))) then  
          itemp = indx(l+1)  
          indx(l+1) = indx(ir)  
          indx(ir) = itemp  
     end if  

     if (arr(indx(l)).gt.arr(indx(ir))) then  
          itemp = indx(l)  
          indx(l) = indx(ir)  
          indx(ir) = itemp  
     end if  

     if (arr(indx(l+1)).gt.arr(indx(l))) then  
          itemp = indx(l+1)  
          indx(l+1) = indx(l)  
          indx(l) = itemp  
     end if  

     i = l + 1  
     j = ir  
     indxt = indx(l)  
     a = arr(indxt)  

   60      continue  

     i = i + 1  
     if (arr(indx(i)).lt.a) go to 60  

   70      continue  

     j = j - 1  
     if (arr(indx(j)).gt.a) go to 70  
     if (j.lt.i) go to 80  
     itemp = indx(i)  
     indx(i) = indx(j)  
     indx(j) = itemp  
     go to 60  

   80      continue  

     indx(l) = indx(j)  
     indx(j) = indxt  
     jstack = jstack + 2  
     if (jstack.gt.nstack) stop 'nstack too small in indexx_int'  
     if (ir-i+1.ge.j-l) then  
          istack(jstack) = ir  
          istack(jstack-1) = i  
          ir = j - 1  
     else  
          istack(jstack) = j - 1  
          istack(jstack-1) = l  
          l = i  
     end if  

end if  

go to 20  

end
!
!-----------------------------------------------------------------------
!
subroutine overlap_noup(nd, clfac)
!
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
!  
! Note(2): uses reduced line-list  
!  
use nlte_type
USE nlte_dim, only: id_ndept
use fund_const, ONLY: clight
USE princesa_var, only: gl
USE nlte_var, only: sr,vmax,lwion1
use cmf_all_var, ONLY: id1,xlam1,gf1,icfirst,iclast,nftot,index_lam_cmf,lam_cmf, &
&           indexel_inv,nf_cmf,occexpl,vdoptot,vref,profdop,sumopal_cmf, &
&           id1_overlap,xm,minrho,no2          

implicit none

integer(i4b),intent(in) :: nd

real(dp), parameter :: c_tau=1.d-8*0.02654! , hc2_8=hc2*1.d24 !conversion to Angstrom

integer(i4b), parameter :: nd1=id_ndept

real(dp), dimension(nd1), intent(in) :: clfac

real(dp), dimension(nd1) :: opal, occlow, occup

integer(i4b) :: lastindex, i, in, ii, index, kk, k, irest, ml, mu, ixmax, &
                ix, l, ik, index1, up, no1      

real(dp) :: lamgrid, lam, srvmax, c_tau1, const, rho1, mirho, lam1, deltav, &
            dw1, dw2

!id_overlap allocated in subr. line-list
!id_overlap initialized (=0) in subr. opacity_cmf for ipath=1

srvmax=sr/vmax
c_tau1=c_tau*srvmax

lastindex=icfirst-1

do i=1,nftot
  in=index_lam_cmf(i)
  lamgrid=lam_cmf(i)  

  if(in.ne.0) then
      overlap: do ii=lastindex+1,in
        index=id1(ii)
        if (index.lt.1d8) cycle overlap ! not explicit element
        kk = index/1d8
        k=indexel_inv(kk)
        if(k.eq.1 .or. k.eq.2) cycle overlap ! H or He
        irest = index - kk*1d8
        ml = irest/10000
        mu = irest - ml*10000
! in units of vref = vturb/3 (with vturb=max(vturb,vturbmin)
        ixmax=nf_cmf(k)
        lam=xlam1(ii)
!        print*
!        print*,ml,mu,lam
        const=c_tau1*lam*gf1(ii)
        occlow=occexpl(ml,:)/gl(ml)
         occup=occexpl(mu,:)/gl(mu)
! pi e2/me c * 1.d-8 *SR/vmax * lam * gf * (nl/gl - nu/gu) /clfac
        opal=const*(occlow-occup)/clfac
        if(opal(1).eq.0.) stop ' opal = 0 for expl. element, subr. overlap_noup)'
! below LTE for bg-elements; note: lwion1 (selected) > lwion (non-selected) 
        do l=1,lwion1-1
          ix=xm*vdoptot(l,k)/vref+1 !in units of vref/3
          if(ix.gt.ixmax) stop ' ix > ixmax in subr. overlap_nopup'
          mirho=1000.
          do ik=0,ix
             rho1=opal(l)*profdop(ik,l,k)/sumopal_cmf(l,i+ik)
!  significant or dominating inversion of bg opacities
             if(rho1.gt.1.1d0 .or.  rho1.lt.0.d0) rho1=0.d0
             mirho=amin1(mirho,rho1)
          enddo   
          do ik=1,ix
             rho1=opal(l)*profdop(ik,l,k)/sumopal_cmf(l,i-ik)
!  significant or dominating inversion of bg opacities
             if(rho1.gt.1.1d0 .or.  rho1.lt.0.d0) rho1=0.d0
             mirho=amin1(mirho,rho1)
          enddo   
          if(mirho.gt.minrho) then
!            print*,ml,mu,l,lam,mirho
            do ik=ii+1,iclast
                index1=id1(ik)
                if(index1.ge.1d8) cycle
                irest=index1/100
                up=index1-irest*100
                if(up.ne.0) cycle
                lam1=xlam1(ik)
                deltav=(lam1/lam-1.d0)*clight
                if(deltav.gt.xm*vdoptot(l,k)) exit
!                print*,'overlap with'
!                print*,lam1,deltav/1.d5,index1
                if(id1_overlap(ik).eq.0) then
                  id1_overlap(ik)=ii
                else
                  if(id1_overlap(ik).ne.ii) then
!in this case, line ik is closer to new line ii (since line ik is longwards
!of new line ii  (i.e., even more longwards from old line ii)
                    id1_overlap(ik)=ii
!                    print*,'2nd line(+)',l,id1(id1_overlap(ik)),id1(ii),index1
!                    print*,'2nd line(+)',xlam1(id1_overlap(ik)),xlam1(ii),lam1
                  endif
                endif  
                
              enddo   

             do ik=ii-1,icfirst,-1
                index1=id1(ik)
                if(index1.ge.1d8) cycle
                irest=index1/100
                up=index1-irest*100
                if(up.ne.0) cycle
                lam1=xlam1(ik)
                deltav=(1.d0-lam1/lam)*clight
                if(deltav.gt.xm*vdoptot(l,k)) exit
!                print*,'overlap with'
!                print*,lam1,deltav/1.d5,index1
                if(id1_overlap(ik).eq.0) then
                  id1_overlap(ik)=ii
                else
                  if(id1_overlap(ik).ne.ii) then
!here, line ik is in between old and new line ii. For reasons of simplicity
!we choose the closer one as being of major impact. This is justified, since
!both lines have a signficant impact (min(rho1) > minrho 

!ik can be even shortwards of old ii
                    dw1=abs(lam1-xlam1(id1_overlap(ik))) 
!ik must be shortwards of new ii
                    dw2=lam-lam1
                    if(dw2.lt.0.) stop ' dw2 < 0 in subr. overlap_noup' 
!line ik closer to new line ii
                    if(dw2.lt.dw1) id1_overlap(ik)=ii
!                    print*,'2nd line(-)',l,id1(id1_overlap(ik)),id1(ii),index1
!                    print*,'2nd line(-)',xlam1(id1_overlap(ik)),xlam1(ii),lam1
! JO keep this stop statement; test for more explicit elements
! if stop occurs, design a 2nd array id2_overlap and include the 2nd line
! in opacity_cmf, then test which line affects what 
!                    stop
                  endif
                endif  
             enddo   
          endif !mirho
        enddo !depth-loop
      enddo overlap! all overlapping lines at cmf freq. lamgrid
      lastindex=in
  endif ! if lines present at lamgrid
enddo ! all cmf-frequencies 

no1=0
no2=0
do ii=icfirst,iclast
  index=id1(ii)
  if(index.ge.1d8) cycle
  irest=index/100
  up=index-irest*100
  if(up.ne.0) cycle
  no1=no1+1
  if(id1_overlap(ii).ne.0) no2=no2+1
enddo

print* 
print*,no1,' lines from background elements with no upper level'
print*,no2,' from those with modified source-function'

return
end
!
!-----------------------------------------------------------------------
!
subroutine transic_count(lineno,retchar)
!
USE nlte_type
USE nlte_dim
USE ffr_error
IMPLICIT NONE
!
!        provides number of lines in file LINES_xxx.dat
!
!     ..
integer(i4b) ::  lineno
character retchar*4,ret*4  
!     ..
!     .. local scalars ..
real(dp) ::  realo
integer(i4b) ::  i,integ,nc,ndatos
character key*6  

!     ..
!     .. external subroutines ..
external ffrkey,ffrnum  
!     ..

lineno=0

   10 continue  
call ffrkey(key,nc,ret)  
if (on_error(ret))  goto 110
if (key.ne.'CL') go to 120  

call ffrkey(key,nc,ret)  
if (on_error(ret))  goto 110
if (key.ne.'TY') go to 120  

call ffrkey(key,nc,ret)  
if (on_error(ret))  goto 110
if (key.ne.'RBB') go to 120  

call ffrnum(realo,integ,ret)  
if (on_error(ret))  goto 110
if (integ.ne.2) go to 120

call ffrnum(realo,integ,ret)  
if (on_error(ret))  goto 110
ndatos = integ  
call ffrkey(key,nc,ret)  
if (on_error(ret))  goto 110
if (key.ne.'RBB') go to 120  

call ffrkey(key,nc,ret)  
if (on_error(ret))  goto 110

   20 CONTINUE  

lineno=lineno+1

call ffrkey(key,nc,ret)  
if (on_error(ret))  goto 110

call ffrnum(realo,integ,ret)  !nline  
if (on_error(ret))  goto 110


!       lectura de los datos
do i = 1,ndatos  
     call ffrnum(realo,integ,ret)  
     if (on_error(ret))  goto 110
end do  

!       leer nueva keyword
call ffrkey(key,nc,ret)  
if (on_error(ret))  goto 110

if (key.eq.'0') go to 10  

go to 20  
!
! error handling
!
110 select case(ret)
case('RET1') !       end of file exit 
  write (*,fmt='(a)') ' end of file '  
  retchar='RET1'
  return
case('RET2') !       error exit 
  write (*,fmt='(a)') 'error in ffr subroutines '  
  retchar='RET2'
  return
case default
  stop ' wrong error condition in transic_count'
end select

120  if (key.eq.'THEEND') then
    retchar='RET0'
    return
else
    stop ' wrong end condition in transic_count'
endif

end
!
!-----------------------------------------------------------------------
!
subroutine transic_read(lineno,retchar)
!
USE nlte_type
USE nlte_dim
USE ffr_error
USE princesa_var, ONLY: labl, gl
USE cmf_all_var, ONLY: depacked
IMPLICIT NONE
!
!        provides number of lines in file LINES_xxx.dat
!
!     ..
integer(i4b), intent(in) ::  lineno
character retchar*4,ret*4  
!     ..
!     .. local scalars ..
real(dp) ::  realo, sumglo, xkw, xn, vac
integer(i4b) ::  i,integ,nc,ndatos,nline,l
character key*6  

!     ..
!     .. external subroutines ..
external ffrkey,ffrnum  
!     ..

nline=0

   10 continue  
call ffrkey(key,nc,ret)  
if (on_error(ret))  goto 110
if (key.ne.'CL') go to 120  

call ffrkey(key,nc,ret)  
if (on_error(ret))  goto 110
if (key.ne.'TY') go to 120  

call ffrkey(key,nc,ret)  
if (on_error(ret))  goto 110
if (key.ne.'RBB') go to 120  

call ffrnum(realo,integ,ret)  
if (on_error(ret))  goto 110
if (integ.ne.2) go to 120

call ffrnum(realo,integ,ret)  
if (on_error(ret))  goto 110
ndatos = integ  
call ffrkey(key,nc,ret)  
if (on_error(ret))  goto 110
if (key.ne.'RBB') go to 120  

call ffrkey(key,nc,ret)  
if (on_error(ret))  goto 110

   20 CONTINUE  

nline=nline+1
depacked(nline)%levlo = key

call ffrkey(key,nc,ret)  
if (on_error(ret))  goto 110
depacked(nline)%levup = key

call ffrnum(realo,integ,ret)  !nline  
if (on_error(ret))  goto 110
depacked(nline)%lineno = integ

!       read remaining data
   call ffrnum(realo,integ,ret)  
   if (on_error(ret))  goto 110

! convert air to vacuum (for lambda > 2000)
   if(realo.gt.2000.) then
! two iterations are sufficient (until 3rd digit in Angstrom)
     xkw=1.d4/realo      
     xn=1.d0+1.d-7*(643.28d0+294981.d0/(146.d0-xkw**2)+2554.d0/(41.d0-xkw**2))
     vac=realo*xn ! first guess for vacuum wavelength;
                  ! actually, xn needs to be calculated from vacuum
     xkw=1.d4/vac
     xn=1.d0+1.d-7*(643.28d0+294981.d0/(146.d0-xkw**2)+2554.d0/(41.d0-xkw**2))
     realo=realo*xn ! 2nd iteration   
   endif
     
   depacked(nline)%wave = realo

   call ffrnum(realo,integ,ret)  
   if (on_error(ret))  goto 110
   depacked(nline)%flu = realo

   call ffrnum(realo,integ,ret)  
   if (on_error(ret))  goto 110
!   gammal = realo  

   call ffrnum(realo,integ,ret)  
   if (on_error(ret))  goto 110

   call ffrnum(realo,integ,ret)  
   if (on_error(ret))  goto 110

   call ffrnum(realo,integ,ret)  
   if (on_error(ret))  goto 110
!   gammac = realo  

   call ffrnum(realo,integ,ret)  
   if (on_error(ret))  goto 110
   sumglo = realo
! consistency test
   do l=1,id_llevs
     if(depacked(nline)%levlo .eq. labl(l)) then
       if(sumglo .ne. gl(l)) then
         print*,labl(l),sumglo,gl(l)
         stop ' inconsistent statistical weights in LINES_xxx.dat'
       endif  
       goto 100
! NOTE: there might be more levels in LINES_xxx.dat than in model atom
     endif 
   enddo  
! for tests
!   print*,' level ',depacked(nline)%levlo,' not found in model atom'

!       leer nueva keyword
 100 call ffrkey(key,nc,ret)  
if (on_error(ret))  goto 110

if (key.eq.'0') go to 10  

go to 20  
!
! error handling
!
110 select case(ret)
case('RET1') !       end of file exit 
  write (*,fmt='(a)') ' end of file '  
  retchar='RET1'
  return
case('RET2') !       error exit 
  write (*,fmt='(a)') 'error in ffr subroutines '  
  retchar='RET2'
  return
case default
  stop ' wrong error condition in transic_read'
end select

120  if (key.eq.'THEEND') then
    retchar='RET0'
    if(nline .ne. lineno) stop ' nline ne lineno in subr. transic_read'    
    return
else
    stop ' wrong end condition in transic_read'
endif

end
!
!-----------------------------------------------------------------------
!
subroutine sort_line_list(n,xl,gf,id)  

use nlte_type
implicit none
!
!
!     sorts for absolute values
!                                                                       
!
!     .. scalar arguments ..
integer(i4b) ::  n  
!     ..
!     .. array arguments ..
real(sp) ::  xl(n), gf(n)  
integer(i4b) :: id(n)
!     ..
!     .. local scalars ..
real(sp) ::  xxl,xxg  
integer(i4b) ::  i,ir,j,l,xxi  
!     ..

l = n/2 + 1  
ir = n  

   10 continue  

if (l.gt.1) then  
     l = l - 1  
     xxl = xl(l)  
     xxg = gf(l)  
     xxi = id(l)  
else  
     xxl = xl(ir)  
     xxg = gf(ir)  
     xxi = id(ir)  
     xl(ir) = xl(1)  
     gf(ir) = gf(1)  
     id(ir) = id(1)  
     ir = ir - 1  

     if (ir.eq.1) then  
          xl(1) = xxl  
          gf(1) = xxg  
          id(1) = xxi  
          return  
     end if  
end if  

i = l  
j = l + l  

   20 continue  

if (j.le.ir) then  

     if (j.lt.ir) then  
          if (xl(j).lt.xl(j+1)) j = j + 1  
     end if  

     if (xxl.lt.xl(j)) then  
          xl(i) = xl(j)  
          gf(i) = gf(j)  
          id(i) = id(j)  
          i = j  
          j = j + j  
     else  
          j = ir + 1  
     end if  

     go to 20  

end if  

xl(i) = xxl  
gf(i) = xxg  
id(i) = xxi  

go to 10  

end
!
!-----------------------------------------------------------------------
!
subroutine prep_ng(ilow,imax,r,flag_ng)  
! 
! prepare and perform Ng-extrapolation for explicit elements
!
! note: occ_expl refers to 'old' occupation numbers, calculated
! one iteration before last call of rateeq.
! thus, at itng=1 we read S(n-3) from occexpl 
!       at itng=2 we read S(n-2) from occexpl
!       at itng=3 we read S(n-1) from occexpl 
! AND                     S(n)   from unit 17 

  
use nlte_type
use nlte_dim, only: id_nttrd, id_ndept, id_atoms, id_llevs
use fund_const, only: hc2
use princesa_var, only: ifirsl, iong, indat1, data, labat, le, li, gl, labl
use nlte_var, only: mmlow, mmup, indexrbb, indxlamc, almost_converged, lablineup
use cmf_all_var, only: indexcmf1, occexpl, indexrbb_inv
use ng_var, only: itng, ng, no_ng_el, ng_el, n_ng, gfmin_ng, &
&   index1_ng, index_ng, sl_ng, trans_ng
implicit none

integer(i4b), parameter :: nd1=id_ndept, kel=id_atoms, nrec=id_llevs

integer(i4b), dimension(nd1,kel), intent(in) :: ilow, imax
real(dp), dimension(nd1), intent(in) :: r

integer(i4b), dimension(kel) :: nfir, nlas

integer(i4b) :: k, ilo, ima, numion, igenio, j, ii, indi, ml, mu, nato, &
&               i, indmat, indcmf, ion, ll

real(dp), dimension(nd1) :: jbar, sl1n, sl1n1, sl1n2, sl1n3, w, vec

real(dp) :: gf, xlamc, lam3, glow, gup, const, xnl, xnu, errmax

real(dp) :: a1, a2, b1, b2, c1, c2, aa, bb

logical first, flag_ng

data first/.true./

ng=.false.

if (first) then
  allocate(index1_ng(id_nttrd))
  first=.false.
else
  if(n_ng.eq.0) then
    if(itng.ne.1) stop ' error in philosophy: n_ng = 0 and itng ne 1 (prep_ng)' 
    return
  endif
endif  


itng = itng + 1

if(itng .lt. 1) return

if(itng.eq.1) then
! set up index-array for resonance lines treated with Ng

! only for ions which are present everywhere  
  do k=1,kel
    ilo = maxval(ilow(:,k))  
    ima = minval(imax(:,k))  
    if (ima.lt.ilo) stop ' ima < ilo in prep_ng, itng=1'
    numion = igenio(k,ilo)  
    nfir(k) = ifirsl(numion)  
    numion = igenio(k,ima) + 1  
    if (numion.gt.iong) stop 'error in prep_ng - numion'  
    nlas(k) = ifirsl(numion)  
!    print*,k,nfir(k),nlas(k)
  enddo

! find resonance transitions in ng-selected elements (ng_el)
  print*
  print*,' Resonance lines (explicit elements) potentially extrapolated by Ng' 
  n_ng=0
  iiloop: do j = 1,id_nttrd  
   ii=indexrbb_inv(j)
! ii is the index w.r.t. index1
   indi = indat1(ii)  
   ml = mmlow(ii)  
   mu = mmup(ii)  
   nato=le(ml)
   do i=1,no_ng_el
     if(labat(nato).eq.ng_el(i)) goto 10
   enddo
   cycle

10 indmat = indxlamc(j)
   indcmf = indexcmf1(indmat) !use new index
   if(indcmf.lt.3) cycle !NOTE INCLUDE AFTER TESTS
   
   if (ml.lt.nfir(nato) .or. ml.ge.nlas(nato)) cycle
   if (mu.gt.nlas(nato)) cycle ! should be read 'ge' instead of 'gt'
   ion=li(ml)
   numion = igenio(nato,ion)  
!   print*,ml,mu,numion,ifirsl(ion)
   if(ml.ne.ifirsl(numion)) cycle !only resonance lines
   gf = data(indi+2) * gl(ml)  
   if (gf.lt.gfmin_ng) cycle ! no forbidden lines
! if HeII 303 should be treated, uncomment following line
!   if(labat(nato).eq.'HE' .and. lablineup(indmat).ne.'HE22  ') cycle
   print*,labat(nato),ion,ml,mu,data(indi+1)   
   n_ng=n_ng+1   
   index1_ng(n_ng)=ii
  enddo iiloop

  if(n_ng.eq.0) then
  print*,' No resonance lines (explicit) to be extrapolated by Ng found' 
  print*
  return
  endif

  print*
! for each extrapolation cycle, new allocation, since n_ng might have
! changed (new ions)
  if (allocated(index_ng)) then
    deallocate(index_ng)
    deallocate(sl_ng)
    deallocate(trans_ng)
  endif  
  allocate(index_ng(n_ng))
! 4 entries for SL, the last one for Jbar from calc_tlu
  allocate(sl_ng(nd1,n_ng,5))
  allocate(trans_ng(n_ng,2))
  index_ng=index1_ng(1:n_ng)
  sl_ng=0.d0  

  index1_ng=0
  do i = 1,n_ng
    ii=index_ng(i)
    indmat=indexrbb(ii)
    index1_ng(indmat)=i
  enddo
endif
  
if(itng.ge.1 .and. itng.le.3) then
! calculate and save source-functions from occexpl (see above)
ngloop: do i = 1,n_ng
  ii=index_ng(i)
  indi = indat1(ii)  
  xlamc = data(indi+1)*1.d-8  
  lam3=xlamc**3 !in cgs
  
  ml = mmlow(ii)  
  mu = mmup(ii)  
  glow=gl(ml)
  gup=gl(mu)

  const=hc2/lam3
  do ll=1,nd1
    xnl = occexpl(ml,ll)/glow  
    xnu = occexpl(mu,ll)/gup  
    sl_ng(ll,i,itng)=const/(xnl/xnu-1.d0)
  enddo
  if(itng.eq.1) then
    trans_ng(i,1)=ml
    trans_ng(i,2)=mu
  endif
  enddo ngloop
endif

if(itng.eq.3) then
   if(almost_converged) then
     do i = 1, n_ng
! do not use extrapolated value 
! in subr. opacity_cmf and calc_tlu, update of
! sline -> sl_ng only performed when sl_ng(35,i,4) > 0.
       sl_ng(:,i,4) = 0.
     enddo
     print*
     print*,' Source-functions for selected resonance lines (expl.) NOT extrapolated,'
     print*,' since ALMOST_CONVERGED previously set to T!'
     print*
     itng=-2 ! next two iterations without Ng
     ng=.true.
     return

   else      
! at first, we read the actual SL from last call of rateeq
! here, we have to change the order (outer loop = depth)
   do ll=1,nd1
! note that occexpl (cmf_all_var) is overwritten; this should not
! matter though, since at this state it's no longer used!
! for next iteration (opacities etc., it's re-read again anyhow)
    read (17,rec=ll) (occexpl(ii,ll),ii=1,nrec)  

    ngloop1: do i = 1,n_ng
      ii=index_ng(i)
      indi = indat1(ii)  
      xlamc = data(indi+1)*1.d-8  
      lam3=xlamc**3 !in cgs
  
      ml = mmlow(ii)  
      mu = mmup(ii)  
      glow=gl(ml)
      gup=gl(mu)

      const=hc2/lam3
    
      xnl = occexpl(ml,ll)/glow  
      xnu = occexpl(mu,ll)/gup  
      sl_ng(ll,i,4)=const/(xnl/xnu-1.d0) !this is the actual SL (index=4)
    enddo ngloop1
  enddo

! extrapolate latest 4 source-functions to new one using Ng-extrapolation
! inlined to save time
  do i=1,n_ng  
    ii=index_ng(i)
    ml = mmlow(ii)  
    mu = mmup(ii)  
    jbar=sl_ng(:,i,5)
    sl1n=sl_ng(:,i,4)
    sl1n1=sl_ng(:,i,3)
    sl1n2=sl_ng(:,i,2)
    sl1n3=sl_ng(:,i,1)
! round to 4 digits, to avoid numerical problems
    call round(sl1n,4,nd1)
    call round(sl1n1,4,nd1)
    call round(sl1n2,4,nd1)
    call round(sl1n3,4,nd1)
! not calculated or inversion, set extrapolated value to zero    
! in subr. opacity_cmf and calc_tlu, update of
! sline -> sl_ng only performed when sl_ng(35,i,4) > 0.
!    
    if(minval(sl1n).le.0.d0) then
      write(*,fmt='("sl1n le 0, no ng-extrapolation:",i3,2(2x,a6))') &
&      i,labl(ml),labl(mu)
      sl_ng(:,i,4) = 0.
      cycle
    endif  
    if(minval(sl1n1).le.0.d0) then
      write(*,fmt='("sl1n1 le 0, no ng-extrapolation:",i3,2(2x,a6))') &
&      i,labl(ml),labl(mu)
      sl_ng(:,i,4) = 0.
      cycle
    endif  
    if(minval(sl1n2).le.0.d0) then
      write(*,fmt='("sl1n2 le 0, no ng-extrapolation:",i3,2(2x,a6))') &
&      i,labl(ml),labl(mu)
      sl_ng(:,i,4) = 0.
      cycle
    endif  
    if(minval(sl1n3).le.0.d0) then
      write(*,fmt='("sl1n3 le 0, no ng-extrapolation:",i3,2(2x,a6))') &
&      i,labl(ml),labl(mu)
      sl_ng(:,i,4) = 0.
      cycle
    endif  

    do ll=1,nd1
      if(jbar(ll).le.0.d0) then
        print*,i,ll
        stop ' jbar le 0 in prep_ng (itng=3)'
      endif
    enddo  
! weights according to tests (read_it_ng.pro) 
    w=1./(jbar*r*r)**2
!    w=1./(jbar*r*r)
!for tests
!    w(:)=1.d0

    vec=(sl1n-2.d0*sl1n1+sl1n2)**2*w
    a1=sum(vec)
    vec=(sl1n-sl1n1-sl1n2+sl1n3)*(sl1n-2.d0*sl1n1+sl1n2)*w
    b1=sum(vec)
    a2=b1
    vec=(sl1n-sl1n1-sl1n2+sl1n3)**2*w
    b2=sum(vec)
    vec=(sl1n-2.d0*sl1n1+sl1n2)*(sl1n-sl1n1)*w
    c1=sum(vec)
    vec=(sl1n-sl1n1-sl1n2+sl1n3)*(sl1n-sl1n1)*w
    c2=sum(vec)
    aa=(c1*b2-c2*b1)/(a1*b2-a2*b1)
    bb=(c2*a1-c1*a2)/(a1*b2-a2*b1)
    sl_ng(:,i,4)=(1.d0-aa-bb)*sl1n+aa*sl1n1+bb*sl1n2
    if(minval(sl_ng(:,i,4)).le.0.d0) then
      print*,i,' extrapolated source function le 0 (prep_ng)'      
      print*,'NO EXTRAPOLATION!!!' 
!      stop ' extrapolated source function le 0 (prep_ng)'
      sl_ng(:,i,4) = sl1n
    endif  
    errmax=maxval(abs(1.-sl_ng(:,i,4)/sl1n(:)))
    write(*,fmt='("Ng-extrapolation (expl.): ",i3,2(2x,a6),2x,e10.4)') &
&      i,labl(ml),labl(mu), errmax
!    do ll=1,nd1
!      write(*,fmt='(i5,5(2x,e10.4))'),ll,sl1n3(ll),sl1n2(ll),sl1n1(ll),sl1n(ll),sl_ng(ll,i,4)
!    enddo     
  enddo
  print*
  print*,' Source-functions for selected resonance lines (expl.) extrapolated!'
  print*
  itng=-2 ! next two iterations without Ng
  ng=.true.
  flag_ng=.true.
  
  endif ! .NOT. ALMOST_CONVERGED
endif !it_ng = 3 

return
end
!
!-----------------------------------------------------------------------
!
subroutine prep_ng_met(r, flag_ng_met)  
! 
! prepare and perform Ng-extrapolation for selected background elements
!
! note: occngold2 refers to 'old' occupation numbers, calculated
! one iteration before last call of rateeq.
! thus, at itng_met=1 we read S(n-3) from occngold2 
!       at itng_met=2 we read S(n-2) from occngold2
!       at itng_met=3 we read S(n-1) from occngold2
! AND                         S(n)   from occng
use nlte_type
use fund_const, only: hc2
use nlte_dim, only: id_ndept
use nlte_var, only: lwion1, metconv2
USE nlte_app, only: met_imin, met_imax, names1, met_rbb, natom, jatom_full, &
&                   indexel, indrec, irbb, ilin, occng, occngold2
use ng_var, only: itng_met, ng_met, n_ng_met, no_ng_el, ng_el, gfmin_ng, &
                  index_ng_met, sl_ng_met, trans_ng_met
use nlte_lines, only: metconv1

implicit none

integer(i4b), parameter :: nd1=id_ndept

real(dp), dimension(nd1), intent(in) :: r

integer(i4b) :: lw1, ilow, imax, k, j, i, n, jj, ii, ji, iqual, ml, mu, &
&               ind, irest, ll

real(dp), dimension(:), allocatable, save :: jbar, sl1n, sl1n1, sl1n2, sl1n3, w, vec

real(dp) :: gf, xlamc, lam3, const, xnl, xnu, errmax

real(dp) :: a1, a2, b1, b2, c1, c2, aa, bb

logical first, flag_ng_met

data first/.true./

ng_met=.false.
lw1=lwion1-1

if (first) then
  allocate(jbar(lw1), sl1n(lw1), sl1n1(lw1), sl1n2(lw1), sl1n3(lw1), w(lw1), vec(lw1))
  first=.false.
else
  if(n_ng_met.eq.0) then
    if(itng_met.ne.1) stop ' error in philosophy: n_ng_met = 0 and itng_met ne 1 (prep_ng_met)' 
    return
  endif
endif  

itng_met = itng_met + 1

if(itng_met .lt. 1) return

if(itng_met.eq.1) then
! set up index-array for resonance lines treated with Ng

  print*
  print*,' Resonance lines (background elements) potentially extrapolated by Ng' 
  n_ng_met=0
  do k=1,natom
    if (indexel(k).gt.0) cycle    !explicit element
    if (jatom_full(k).eq.0) cycle !only selected elements

! find resonance transitions in ng-selected elements (ng_el)
    do i=1,no_ng_el
      if(names1(k).eq.ng_el(i)) goto 10
    enddo
    cycle

! only for ions which are present everywhere  
10  ilow=maxval(met_imin(k,1:lw1))
    imax=minval(met_imax(k,1:lw1))

    do j=ilow,imax
!hack      
!      if(j.ne.3) cycle
      n=indrec(k,j)
      jj=irbb(n)-1

      do ii=1,ilin(n)
        ji=jj+ii
        iqual=met_rbb(ji)%quality
! outside range or no transition in line-list,
! will be treated by conventional (simplified) approach in netmat_met,
! see calc_tlu_met
        if(iqual.eq.-1 .or.iqual.eq.2) cycle  

        ml=met_rbb(ji)%low  
        if(ml.ne.1) cycle !only resonance lines
        gf=met_rbb(ji)%gf  
        if (gf.lt.gfmin_ng) cycle ! no forbidden lines
        mu=met_rbb(ji)%lup
!hack      
!        if(mu.ne.10) cycle
        xlamc=met_rbb(ji)%wave !packed transition, in A
        ind=k*1000000+j*100000+ml*100+mu
        print*,names1(k),ind,xlamc   
        n_ng_met=n_ng_met+1   
        if(n_ng_met.gt.1000) stop ' n_ng_met > 1000! Enlarge array trans_ng_met!'
        met_rbb(n_ng_met)%index1_ng=ji !used as dummy-variable
        trans_ng_met(n_ng_met)=ind
      enddo 
    enddo
  enddo
  
  if(n_ng_met.eq.0) then
    print*,' No resonance lines (background) to be extrapolated by Ng found' 
    print*
    return
  endif

  print*
! for each extrapolation cycle, new allocation, since n_ng_met might have
! changed (new ions)
  if (allocated(index_ng_met)) then
    deallocate(index_ng_met)
    deallocate(sl_ng_met)
  endif  
  allocate(index_ng_met(n_ng_met))
! 4 entries for SL, the last one for Jbar from calc_tlu
  allocate(sl_ng_met(lw1,n_ng_met,5))
  index_ng_met=met_rbb(1:n_ng_met)%index1_ng
  sl_ng_met=0.d0  

  met_rbb%index1_ng=0
  do i = 1,n_ng_met
    ji=index_ng_met(i)
    met_rbb(ji)%index1_ng=i
  enddo
endif

if(itng_met.ge.1 .and. itng_met.le.3) then
! calculate and save source-functions from occexpl (see above)
ngloop: do i = 1,n_ng_met
  ji=index_ng_met(i)
  xlamc = met_rbb(ji)%wave*1.d-8  
  lam3=xlamc**3 !in cgs
  
  ind=trans_ng_met(i)

  k = ind/1000000
  irest = ind - k*1000000
  j = irest/100000
  irest = irest - j*100000
  ml = irest/100
  mu  = irest-ml*100
   
  n=indrec(k,j)
  const=hc2/lam3

  do ll=1,lw1
!    if(ll.eq.1.or.ll.eq.lw1) &
!&    write(*,fmt='(i3,3(2x,e10.4))') i,1.-occngold2(n,mu,ll)/occng(n,mu,ll), &
!&    occngold2(n,mu,ll),occng(n,mu,ll)
    xnl = occngold2(n,ml,ll)  !occng = n/g  
    xnu = occngold2(n,mu,ll)  
    sl_ng_met(ll,i,itng_met)=const/(xnl/xnu-1.d0)
  enddo

  if(itng_met.eq.3) then
! here, we additionally use the most recent occupation numbers to obtain
! sl1n  
  do ll=1,lw1
    xnl = occng(n,ml,ll)  !occng = n/g  
    xnu = occng(n,mu,ll)  
    sl_ng_met(ll,i,4)=const/(xnl/xnu-1.d0)
  enddo
  endif

enddo ngloop
endif

if(itng_met.eq.3) then
! JO Jan. 2021; changed from metconv1 to (metconv1 .or metconv2);
! tests have shown that for metconv2 = T the convergence is better without
! extrapolation;
! in those cases where metconv2 does not appear (opt_oiii_iter=F), metconv1
! is still used as the decisive critirion.
!  if(metconv1) then
  if(metconv1.or.metconv2) then
     do i = 1, n_ng_met
! do not use extrapolated value 
! in subr. opacity_cmf and calc_tlu, update of
! sline -> sl_ng only performed when sl_ng(35,i,4) > 0.
       sl_ng_met(:,i,4) = 0.
     enddo
     print*
     print*,' Source-functions for selected resonance lines (backgr.) NOT extrapolated,'
     print*,' since METCONV1 or METCONV2 previously set to T!'
     print*
     itng_met=-2 ! next two iterations without Ng
     ng_met=.true.
     return

   else      
! extrapolate latest 4 source-functions to new one using Ng-extrapolation
! inlined to save time
   do i=1,n_ng_met  
    ji=index_ng_met(i)
    ind=trans_ng_met(i)
    k = ind/1000000
    irest = ind - k*1000000
    j = irest/100000
    irest = irest - j*100000
    ml = irest/100
    mu  = irest-ml*100

    jbar=sl_ng_met(:,i,5)
    sl1n=sl_ng_met(:,i,4)    
    sl1n1=sl_ng_met(:,i,3)
    sl1n2=sl_ng_met(:,i,2)
    sl1n3=sl_ng_met(:,i,1)
! round to 4 digits, to avoid numerical problems
    call round(sl1n,4,lw1)
    call round(sl1n1,4,lw1)
    call round(sl1n2,4,lw1)
    call round(sl1n3,4,lw1)
! not calculated or inversion, set extrapolated value to zero
! in subr. opacity_cmf and calc_tlu_met, update of
! sline -> sl_ng_met only performed when sl_ng_met(35,i,4) > 0.
!    
    if(minval(sl1n).le.0.d0) then
      write(*,fmt='("sl1n le 0, no ng-extrapolation:",i3,2x,a2,3(2x,i3))') &
&      i,names1(k),j,ml,mu
      sl_ng_met(:,i,4) = 0.
      cycle
    endif  
    if(minval(sl1n1).le.0.d0) then
      write(*,fmt='("sl1n1 le 0, no ng-extrapolation:",i3,2x,a2,3(2x,i3))') &
&      i,names1(k),j,ml,mu
      sl_ng_met(:,i,4) = 0.
      cycle
    endif  
    if(minval(sl1n2).le.0.d0) then
      write(*,fmt='("sl1n2 le 0, no ng-extrapolation:",i3,2x,a2,3(2x,i3))') &
&      i,names1(k),j,ml,mu
      sl_ng_met(:,i,4) = 0.
      cycle
    endif  
    if(minval(sl1n3).le.0.d0) then
      write(*,fmt='("sl1n3 le 0, no ng-extrapolation:",i3,2x,a2,3(2x,i3))') &
&      i,names1(k),j,ml,mu
      sl_ng_met(:,i,4) = 0.
      cycle
    endif  

    do ll=1,lw1
      if(jbar(ll).le.0.d0) then
        print*,i,ll
        stop ' jbar le 0 in prep_ng (itng_met=3)'
      endif
    enddo  
! weights according to tests (read_it_ng.pro) 
    w=1./(jbar*r(1:lw1)*r(1:lw1))**2
!    w=1./(jbar*r(1:lw1)*r(1:lw1))
!for tests
!    w(:)=1.d0

    vec=(sl1n-2.d0*sl1n1+sl1n2)**2*w
    a1=sum(vec)
    vec=(sl1n-sl1n1-sl1n2+sl1n3)*(sl1n-2.d0*sl1n1+sl1n2)*w
    b1=sum(vec)
    a2=b1
    vec=(sl1n-sl1n1-sl1n2+sl1n3)**2*w
    b2=sum(vec)
    vec=(sl1n-2.d0*sl1n1+sl1n2)*(sl1n-sl1n1)*w
    c1=sum(vec)
    vec=(sl1n-sl1n1-sl1n2+sl1n3)*(sl1n-sl1n1)*w
    c2=sum(vec)
! Jo Sept 18
! if identical source-functions (when bg elements kept constant)
    if(a1*b2-a2*b1 .eq. 0.) then
      if(c1*b2-c2*b1 .ne.0) stop ' ng-extrapolation of bg: problems with aa' 
      if(c2*a1-c1*a2 .ne.0) stop ' ng-extrapolation of bg: problems with bb'
      aa=0.5d0
      bb=0.5d0
    else  
      aa=(c1*b2-c2*b1)/(a1*b2-a2*b1)
      bb=(c2*a1-c1*a2)/(a1*b2-a2*b1)
    endif
    sl_ng_met(:,i,4)=(1.d0-aa-bb)*sl1n+aa*sl1n1+bb*sl1n2
    errmax=maxval(abs(1.-sl_ng_met(:,i,4)/sl1n(:)))
    if(minval(sl_ng_met(:,i,4)).le.0.d0) then
      print*,i,' extrapolated source function le 0 (prep_ng_met)'      
      print*,'NO EXTRAPOLATION!!!' 
!      stop ' extrapolated source function le 0 (prep_ng)'
      sl_ng_met(:,i,4) = sl1n !rounded
    endif  
    errmax=maxval(abs(1.-sl_ng_met(:,i,4)/sl1n(:)))
    write(*,fmt='("Ng-extrapolation (bg.): ",i3,2x,a2,3(2x,i3),2x,e10.4)') &
&      i,names1(k),j,ml,mu,errmax
!    do ll=1,lw1
!      write(*,fmt='(i5,5(2x,e10.4))'),ll,sl1n3(ll),sl1n2(ll),sl1n1(ll),sl1n(ll),sl_ng_met(ll,i,4)
!    enddo     
   enddo
   print*
   print*,' Source-functions for selected resonance lines (backgr.) extrapolated!'
   print*
   itng_met=-2 ! next two iterations without Ng
   ng_met=.true.
   flag_ng_met=.true.

   endif ! .NOT. METCONV1

endif !it_ng_met=3      

return
end
!
!-----------------------------------------------------------------------
!
subroutine round(x,prec,nd)  
! 
! rounds vector x of dimension nd to prec
use nlte_type

integer(i4b), intent(in) :: nd, prec
real(dp), dimension(nd) :: x

integer(i4b), dimension(nd) :: expo !automatic array
real(dp), dimension(nd) :: rexp !automatic array


! take care: difference, since INT works different for positive and negative
! arguments
where(x.ge.1.d0)
expo=-int(log10(x))-1+prec
elsewhere
expo=-int(log10(x))+prec
endwhere
rexp=10.d0**expo

x=nint(x*rexp)/rexp
return
end
