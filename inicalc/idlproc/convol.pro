PRO CONVOL,sx,sy,cx,cy,VDOP=vdop,AINST=ain,RESOL=resol,VSINI=vrot,BETA=beta, $
	   vmacro=vmacro, $
           MESSAGE=message, ORIGINAL=original
  
common convol_all, clight  
common share_rot, vrot1, beta1
common share_vmac, vmacro1, ccr1, spi
common share_resol, vdop1, fwhm, dlam_fwhm
  
max_resol=60000
  
; ALL UNITS IN ANGSTROM OR km/s!
  
; original: remap finally on original wavelength grid
; rt-macroturbulence according to Gray, with vmac_t = vmac_r


; keywords vdop, ain, resol, vmacro mutually exclusive.
; rot. folding if present done at first, of course  
 
 
; ain = FWHM
; resolution = lambda/delta lambda
  
;modified version of older routines by Reetz, Fuhrmann and Haser  
;small changes applied and checked for performance. 
 
; version 2.0: allows convolution for large ranges, now accounting
; for lambda-dependent delta lambda (using overlap add method)
; also: smaller changes regarding resampling, and other stuff

IF NOT KEYWORD_SET(beta) THEN beta = 1.5 
IF NOT KEYWORD_SET(message) THEN message = 0
IF NOT KEYWORD_SET(original) THEN original = 0 

clight    = 299792.5
rsmpl_fac = 0.7

sxx       = sx
syy       = sy
n         = N_ELEMENTS(sxx) 
nn        = n
rdst      = (sxx(n-1)-sxx(0))/(n-1)             ; resampling distance

if rdst lt 0. then begin
  print,'decreasing wavelengths, wrong order!'
  return
endif  

intp      = 0

;
; handle not equidistant wavelength vectors
;

dsx        = ABS(sxx(1:*)-sxx(0:N_ELEMENTS(sxx)-2))
rdst_min   = MIN(dsx,MAX=rdst_max)
meanw      = 0.5*(sxx(0)+sxx(n-1))
eps        = meanw/max_resol
diff       = ABS(rdst_max-rdst_min) 
not_eqdist = diff GT eps

IF message THEN PRINT,'Max(dlam): ', diff
IF not_eqdist THEN BEGIN
   IF message THEN $
     PRINT,'Spectrum is not equidistant: RDMIN= ',rdst_min,', RDMAX= ',rdst_max
   rdst = MAX([(rdst_min * rsmpl_fac),eps]) 
   intp = 1
   nn = (sxx(n-1)-sxx(0))/rdst
ENDIF          

xpow = ALOG(nn)/ALOG(2.)
pow = LONG(xpow)

IF ((pow LE 14 AND pow GT 2) OR intp) AND ((xpow-pow) NE 0.) THEN BEGIN
   pow = pow + 1
   nn_x = LONG(2.^pow)
   IF (FLOAT(nn)/FLOAT(nn_x) GT rsmpl_fac) AND (NOT not_eqdist) THEN BEGIN
      nn = nn_x * 2L
   ENDIF ELSE BEGIN
      nn = nn_x
   ENDELSE
   
   rdst = (sxx(n-1)-sxx(0))/(nn-1)
   sxx = LINDGEN(nn)*rdst + sx(0)
   syy = INTERPOL(sy,sx,sxx)
   res=finite(syy)
   bres=where(res eq 0, countb)
   if countb ne 0 then begin

     for i=0,nn-1 do begin
      b=where(sx eq sxx(i),count)
      if count gt 1 then begin
       sym=mean(sy(b))
       syy(i)=sym
      endif
     endfor
   endif

   IF message THEN BEGIN
    PRINT,'Resampled grid with resampling distance: ',rdst
    PRINT, nn, FORMAT= "('Resampled spectrum has ',I6,' points')"
   ENDIF

ENDIF ELSE BEGIN

 IF (DIFF/RDST_MAX GT 0.1) THEN BEGIN

   if rdst ne (sxx(n-1)-sxx(0))/(n-1) or nn ne n then begin
     print,'something rotten in resampling'
     stop
   endif   

   sxx = LINDGEN(nn)*rdst + sx(0)
   syy = INTERPOL(sy,sx,sxx)
   res=finite(syy)
   bres=where(res eq 0, countb)
   if countb ne 0 then begin

     for i=0,nn-1 do begin
      b=where(sx eq sxx(i),count)
      if count gt 1 then begin
       sym=mean(sy(b))
       syy(i)=sym
      endif
     endfor
   endif

   IF message THEN BEGIN
    PRINT,'Resampled grid (with old n) and resampling distance: ',rdst
    PRINT,'previous rdst_min/max = ',rdst_min,rdst_max
    PRINT, nn, FORMAT= "('Resampled spectrum has ',I6,' points')"
   ENDIF

 ENDIF
   
ENDELSE

;if rdst gt eps then begin
;  print,' actual resolution',rdst,' larger than max. resolution', eps
;  stop
;endif  
  
rat=sxx(nn-1)/sxx(0)
if rat gt 1.1 then goto, piecewise

; STANDARD PATH FOR NOT TOO LONG WAVELENGTH RANGES

px = (FINDGEN(nn)-LONG(nn/2))*rdst

; at first, rotational profile, to allow for final Gauss-folding

IF KEYWORD_SET(vrot) THEN BEGIN
 IF vrot ge 5. THEN BEGIN 
   lambda0 = (sxx(nn-1)+sxx(0))/2.
   dlam_max = vrot*(lambda0/clight)
   xx = px / dlam_max
   one = WHERE(ABS(xx) GT 1)
   xx(one) = 1.
   xx_2 = xx*xx
   unsoeld =  (2.*     SQRT(1-xx_2)/!PI + (1-xx_2)*beta/2.) / (1 + 2*beta/3.)
;  gray = (2.*(1-beta)*SQRT(1-xx_2)/!PI + (1-xx_2)*beta/2.) / (1 -   beta/3.) 
   py = unsoeld
   py = py / dlam_max
   CONV_SPEC,sxx,syy,px,py,cx,cy
   cy=float(cy)
   syy = cy
 ENDIF ELSE BEGIN
   print,' vsini too low, no rotational convolution performed'
   cx=sxx
   cy=syy
 ENDELSE  
ENDIF

; *********************************************************
; ------- radial-tangential macroturbulence profile -------
; *********************************************************
; assuming vt=vr
; remember: then only vr necessary, gives same result as with separate calculations
; (otherwise, see routine by Sergio S., don't forget to modify corresponding statements)

IF KEYWORD_SET(vmacro) THEN BEGIN

    AAR=1.

;old    PX2 = PX(NN/2:NN-1)
;new formulation; use all freq. points to avoid asymmetries 

    spi=sqrt(!dpi)
    
    LAMBDA0  = (SXX(NN-1) + SXX(0))/2.
    
    MR=vmacro*lambda0/clight
    CCR=2.*AAR/spi/MR
    
    PXMR=abs(PX)/MR
; from resubstitution y=1/u and partial integration
;    PY2T=CCT*PXMT*[-spi+exp(-PXMT^2)/PXMT+spi*errorf(PXMT)]
;    PY2R=CCR*PXMR*[-spi+exp(-PXMR^2)/PXMR+spi*errorf(PXMR)]
; new
         
;old PY2 = PY2R + PY2T 
;old    PY = [REVERSE(PY2),PY2(0),PY2(0),PY2]
    PY=CCR*[exp(-PXMR^2)+spi*pxmr*(errorf(PXMR)-1.)]
    
    CONV_SPEC,SXX,SYY,PX,PY,CX,CY                  
    CY=FLOAT(CY)
    SYY = CY 
endif   
    
; *********************************************************
; Gauss profile
; *********************************************************


IF KEYWORD_SET(vdop) or KEYWORD_SET(ain) or KEYWORD_SET(resol) THEN BEGIN
; all relations re-checked, OK! (Jo, August 2015)

  fac     = 2.*sqrt(alog(2.))
  lambda0 = (sxx(nn-1)+sxx(0))/2.
  
  if KEYWORD_SET(ain) then begin
    dlam_dop = ain/fac
    vdop     = dlam_dop*clight/lambda0

  endif else if KEYWORD_SET(resol) then begin
    vdop=clight/(resol*fac)
    dlam_dop = vdop*(lambda0/clight)
 
  endif else begin  
    dlam_dop = vdop*(lambda0/clight)
  endelse

  resol=clight/(vdop*fac)
  IF(resol gt max_resol) then begin
;  IF message THEN print,' HiRes anyway, not convolution necessary'
  print,' HiRes anyway, not convolution necessary'
  goto, no_resol
  endif

  
  min_cnr = 5                  
  max_cnr = LONG(nn*0.8)

   IF message THEN BEGIN 
      PRINT,    vdop,    dlam_dop, FORMAT="('Vdop=',F7.3,'km/s  <=>',F7.3,'A')"
      PRINT,fac*vdop,fac*dlam_dop, FORMAT="('FWHM=',F7.3,'km/s  <=>',F7.3,'A')"
   ENDIF
   cnr = LONG(4.*dlam_dop/rdst) + 1       
   IF max_cnr LE cnr THEN BEGIN
      percent = cnr/FLOAT(nn)*100.
      PRINT,'Gauss profile is very large:' + STRING(percent,'(F6.2)') $
                               + '% of spectrum points are forming the profile'
    ENDIF
   IF min_cnr GE cnr THEN BEGIN
   PRINT,'Gauss profile is very small: only ' + STRING(cnr,'(I4)') $  
                                       + ' data points are forming the profile'
   ENDIF
   xx = px / dlam_dop
   six = WHERE(ABS(xx) lt 6.)
   py = fltarr(nn)
   py(six) = EXP(-(xx(six)^2))/(SQRT(!PI)*dlam_dop)
   CONV_SPEC,sxx,syy,px,py,cx,cy       
   cy=float(cy)
   syy = cy
ENDIF

no_resol:

if original then begin
   cy = INTERPOL(syy,sxx,sx)
   cx = sx
endif

return

;--------------------------------------------------------------------------

piecewise:
; PIECEWISE CONVOLUTION 
; such that delta nu (propto 1/lambda) does not vary too much
; (at maximum, 10% change in lambda allowed).
; algorithm according to 'overlap-add-method', 
; and normalized in such a way as to provide identical results
; with the standard method (using CONV_SPEC)

nseg=(rat-1.)/0.1     ; minimum number of segments (float)
lmax=long(nn/nseg)+1  ; maximum length of intervals 

imax=nint(alog(nn)/alog(2.))

; at first, rotational profile, to allow for final Gauss-folding
IF KEYWORD_SET(vrot) THEN BEGIN
 vrot1=vrot
 beta1=beta
 IF vrot ge 5. THEN BEGIN 
   if vrot lt 20. then begin
     print,'vrot < 20 km/s'
     print,'PIECEWISE CONVOLUTION MIGHT PRODUCE ERRONEOUS RESULTS!!!'
     print,'proceed at own risk'
   endif  

; calculate M(max) from highest frequency; will be decreased later on
   lambda0 = sxx(nn-1)
   dlam_max = vrot*(lambda0/clight)
; factor 2, since two times dlam_max
   m=2*long(dlam_max/rdst)+1

; calculate l such that n is a power of two
   for i=1,imax-1 do begin
     l=2l^i-m+1
     if l le lmax then iout=i
   endfor   

   l=2l^iout-m+1
   n=l+m-1
   print,'lambda min/max = ',sxx(0),sxx(nn-1),' nn = ',nn
   print,'piecewise convolution (vsini) with',nn/l+1,' segments'
   print,'L = ',l,' and N = ',n
   if(l gt nn/2) then begin
     print,l,' > ', nn/2
     stop
   endif  
   if(l lt m) then begin
; segment length lower than length of response function 
     print,l,' < ', m,' segment length lower than response'
     print,'reduce vsini'
     stop
   endif  
   if (sxx(l-1)-sxx(0)) gt 0.1*sxx(0) then begin
     print,'interval too large'
     return
   endif  
   
   px = (FINDGEN(n)-LONG(n/2))*rdst

   y=fltarr(nn+m-1)   

   tty = (syy(0)+syy(nn-1))/2.           ; average offset
   delta_n=(sxx(nn-1)-sxx(0))/(nn-1) * n ; renormalization constant 

; calculation of shift   
   xx = px / dlam_max
   core = WHERE(ABS(xx) LE 1.)
;   print,i,m,n_elements(core),dlam_max
   if(n_elements(core) ne m) then begin
     print,'dimension of response function ne M'
     stop
   endif  

   one = WHERE(ABS(xx) GT 1.)
   xx(one) = 1.
   xx_2 = xx*xx
   unsoeld =  (2.*     SQRT(1-xx_2)/!PI + (1-xx_2)*beta/2.) / (1 + 2*beta/3.)
;  gray = (2.*(1-beta)*SQRT(1-xx_2)/!PI + (1-xx_2)*beta/2.) / (1 -   beta/3.) 
   py = unsoeld
   py = py / dlam_max
; calculating shift for broadest response function
; to apply overlap-add method, response function has to be located at the begin
; of the array (otherwise, overlap at both sides of segments)
; (leads to small error at most blueward points)
   offset=where(py ne 0)
   shft=offset(0)

; perform convolution   
   convol_piecewise,y,sxx,syy,px,py,nn,l,m,n,shft,tty,delta_n,'rot'
   
; for next convolution (if there is any)
   syy=y(0:nn-1)

; for output
   cy=syy
   cx=sxx
 endif else begin
   print,' vsini too low, no rotational convolution performed'
   cy=syy
   cx=sxx
 endelse  

endif

; *********************************************************
; ------- radial-tangential macroturbulence profile -------
; *********************************************************
; assuming vt=vr
; remember: then only vr necessary, gives same result as with separate calculations
; (otherwise, see routine by Sergio S., don't forget to modify corresponding statements)

IF KEYWORD_SET(vmacro) THEN BEGIN
   
   vmacro1 = vmacro

   aar=1.
    
   spi=sqrt(!dpi)
    
; calculate M(max) from largest wavelength; will be decreased later on
   lambda0  = sxx(nn-1)
   mr=vmacro*lambda0/clight

; factor 12, since 2*6 times mr (function goes down to \sim 1.d-16)
   
   m=2*long(6.*mr/rdst)+1

; calculate l such that n is a power of two
   for i=1,imax-1 do begin
     l=2l^i-m+1
     if l le lmax then iout=i
   endfor   

   l=2l^iout-m+1
   n=l+m-1
   print,'lambda min/max = ',sxx(0),sxx(nn-1),' nn = ',nn
   print,'piecewise convolution (vmacro) with',nn/l+1,' segments'
   print,'L = ',l,' and N = ',n
   if(l gt nn/2) then begin
     print,l,' > ', nn/2
     stop
   endif  
   if(l lt m) then begin
; segment length lower than length of response function 
     print,l,' < ', m,' segment length lower than response'
     print,'reduce vmacro'
     stop
   endif  
   if (sxx(l-1)-sxx(0)) gt 0.1*sxx(0) then begin
     print,'interval too large'
     stop
   endif  

   px = (FINDGEN(n)-LONG(n/2))*rdst

   y=fltarr(nn+m-1)   

   tty = (syy(0)+syy(nn-1))/2.           ; average offset
   delta_n=(sxx(nn-1)-sxx(0))/(nn-1) * n ; renormalization constant 

; calculation of shift
   pxmr=abs(px)/mr

   core = WHERE(pxmr LE 6.)
;   print,i,m,n_elements(core),mr
   if(n_elements(core) ne m) then begin
     print,'dimension of response function ne M'
     stop
   endif  

; from resubstitution y=1/u and partial integration
;    py2t=cct*pxmt*[-spi+exp(-pxmt^2)/pxmt+spi*errorf(pxmt)]
;    py2r=ccr*pxmr*[-spi+exp(-pxmr^2)/pxmr+spi*errorf(pxmr)]

   ccr1=2.*aar/spi

   py=fltarr(n)
   xx = pxmr(core)
   py(core)=ccr1/mr*[exp(-xx^2)+spi*xx*(errorf(xx)-1.)]

   offset=where(py ne 0)
   shft=offset(0)

; perform convolution   
   convol_piecewise,y,sxx,syy,px,py,nn,l,m,n,shft,tty,delta_n,'vmac'
   
; for next convolution (if there is any)
   syy=y(0:nn-1)

; for output
   cy=syy
   cx=sxx

endif

; *********************************************************
; Gauss profile
; *********************************************************


IF KEYWORD_SET(vdop) or KEYWORD_SET(ain) or KEYWORD_SET(resol) THEN BEGIN
; all relations re-checked, OK! (Jo, August 2015)
  
  fac     = 2.*sqrt(alog(2.))

; for output  
  lambda0 = (sxx(nn-1)+sxx(0))/2.

  if KEYWORD_SET(ain) then begin
    dlam_dop = ain/fac
    vdop     = dlam_dop*clight/lambda0; wavelength dependent, but not required
    fwhm = 1
    dlam_fwhm=dlam_dop; wavelength independent
    
  endif else if KEYWORD_SET(resol) then begin
    vdop=clight/(resol*fac)   ; wavelength independent
    dlam_dop = vdop*(lambda0/clight)
    fwhm = 0
    
  endif else begin  
    dlam_dop = vdop*(lambda0/clight)
    fwhm = 0
  endelse

  vdop1 = vdop; either wavelength independent or not required later on (ain)
  resol=clight/(vdop*fac)
  IF(resol gt max_resol) then begin
;  IF message THEN print,' HiRes anyway, not convolution necessary'
  print,' HiRes anyway, not convolution necessary'
  goto, no_resol_piece
  endif
  
  min_cnr = 5                  
  max_cnr = LONG(nn*0.8)

   IF message THEN BEGIN 
      print,'LONG INTERVAL!'
      print,'following variables evaluated at',lambda0
      print,'during calculation, lambda-dependent dlam_dop accounted for'
      PRINT,    vdop,    dlam_dop, FORMAT="('Vdop=',F7.3,'km/s  <=>',F7.3,'A')"
      PRINT,fac*vdop,fac*dlam_dop, FORMAT="('FWHM=',F7.3,'km/s  <=>',F7.3,'A')"
   ENDIF
   cnr = LONG(4.*dlam_dop/rdst) + 1       
   IF max_cnr LE cnr THEN BEGIN
      percent = cnr/FLOAT(nn)*100.
      PRINT,'Gauss profile is very large:' + STRING(percent,'(F6.2)') $
                               + '% of spectrum points are forming the profile'
    ENDIF
   IF min_cnr GE cnr THEN BEGIN
   PRINT,'Gauss profile is very small: only ' + STRING(cnr,'(I4)') $  
                                       + ' data points are forming the profile'
   ENDIF

; for actual calculation:
; calculate M(max) from largest wavelength; will be decreased later on
  lambda0  = sxx(nn-1)

  if KEYWORD_SET(ain) then begin
; dlam_dop wavelength independent, 
    dlam_dop = ain/fac; dlam_fwhm already calculated above 
    
  endif else if KEYWORD_SET(resol) then begin
    vdop=clight/(resol*fac)
    if abs(1.-vdop1/vdop) gt 1.d-5  then begin
      print,'vdop1 ne vdop'
      stop
    endif  
    dlam_dop = vdop*(lambda0/clight)
 
  endif else begin  
    dlam_dop = vdop*(lambda0/clight)
  endelse

; factor 12, since 2*6 times dlam_dop (function goes down to \sim 1.d-16)
   
   m=2*long(6.*dlam_dop/rdst)+1

; calculate l such that n is a power of two
   for i=1,imax-1 do begin
     l=2l^i-m+1
     if l le lmax then iout=i
   endfor   

   l=2l^iout-m+1
   n=l+m-1
   print,'lambda min/max = ',sxx(0),sxx(nn-1),' nn = ',nn
   print,'piecewise convolution (vmacro) with',nn/l+1,' segments'
   print,'L = ',l,' and N = ',n
   if(l gt nn/2) then begin
     print,l,' > ', nn/2
     stop
   endif  
   if(l lt m) then begin
; segment length lower than length of response function 
     print,l,' < ', m,' segment length lower than response'
     print,'reduce fwhm/vdop or increase resolving power'
     stop
   endif  
   if (sxx(l-1)-sxx(0)) gt 0.1*sxx(0) then begin
     print,'interval too large'
     stop
   endif  

   px = (FINDGEN(n)-LONG(n/2))*rdst

   y=fltarr(nn+m-1)   

   tty = (syy(0)+syy(nn-1))/2.           ; average offset
   delta_n=(sxx(nn-1)-sxx(0))/(nn-1) * n ; renormalization constant 

; calculation of shift
   pxmr = abs(px)/ dlam_dop
; dlam_dop wavelength dependent (vdop, resol) or constant (ain)

   core = WHERE(pxmr LE 6.)
;   print,i,m,n_elements(core),dlam_dop
   if(n_elements(core) ne m) then begin
     print,'dimension of response function ne M'
     stop
   endif  

   py=fltarr(n)
   xx = pxmr(core)
   py(core)=exp(-xx^2)/(SQRT(!PI)*dlam_dop)

   offset=where(py ne 0)
   shft=offset(0)

; perform convolution   
   convol_piecewise,y,sxx,syy,px,py,nn,l,m,n,shft,tty,delta_n,'resol'
   
; for next convolution (if there is any)
   syy=y(0:nn-1)

; for output
   cy=syy
   cx=sxx

endif

no_resol_piece:

if original then begin
   cy = INTERPOL(syy,sxx,sx)
   cx = sx
endif

return
end

;
;------------------------------------------------------------------------
;

PRO CONV_SPEC,sx,sy,px,py,cx,cy

; note:
; no zero padding applied, however sufficient precision also at end points
; note also, that response function centered, not wrapped arround, 
; since referred to central frequency of signal.    
  
n = N_ELEMENTS(sx)
IF n NE N_ELEMENTS(px) THEN PRINT,'ERROR: Different array dimensions'
shft = LONG(n/2)
tty = FLTARR(n) + (sy(0)+sy(n-1))/2. 
sy_t = sy - tty
cx = sx       
cy = SHIFT(FFT(FFT(sy_t,-1)*FFT(py,-1),1),shft)
cy = ((cx(n-1)-cx(0))/(n-1)) * n * cy         ; renorming result of convolution
cy = cy + tty
END

;--------------------------------------------------------------------------

pro convol_piecewise,y,sxx,syy,px,py,nn,l,m,n,shft,tty,delta_n,type
   
; start index (first segment)
   i=0L
; go over all segments
   while i le nn-1 do begin

   il=min([i+l-1,nn-1])
    k=min([i+n-1,m+nn-2])

    
   if type eq 'rot'   then prep_rot,sxx,i,il,m,px,py,shft
   if type eq 'vmac'  then prep_vmac,sxx,i,il,m,n,px,py,shft
   if type eq 'resol' then prep_resol,sxx,i,il,m,n,px,py,shft

   sig=fltarr(n) ; set to zero
; subtract offset (enhances accuracy)
   sig(0:il-i)=syy(i:il)-tty ; zeros from il+1...n-1 padded

   if k eq m+nn-2 then begin
; in last segment, last value of signal padded in between il=nn-1 and il1
; to ensure roughly correct endpoints  
; usually, il1= k = nn + m -1 - 1 (the last -1 for the index).
; under certain circumstances, however, this leads to problems, since
; the last M points in the signal-vector (size N) need to be zero.
; (at least this is my feeling, and tests have shown that this is OK.
; If further problems occur, this needs to be checked more carefully though.)
; Thus, N - (k-i) ge M, and ... 
     il1=min([k,n-m+i])
;     il1=k
     if(il1 lt il) then begin
       print,'il1<il, stop'
       stop
     endif
     sig(il-i+1:il1-i)=syy(nn-1)-tty; padding with last value of signal - offset
   endif else begin
     il1=il
   endelse
; acutal convolution  
   cy = FFT(FFT(sig,-1)*FFT(py,-1),1)
; renormalize result of convolution via delta * n
   cy = delta_n * cy 
; correct for offset
   cy(0:il1-i)=cy(0:il1-i)+tty
;   print,sxx(i),i,il,il1,k
 
; add overlapping convolved signal
   y(i:k)=y(i:k)+float(cy(0:k-i))
;   stop
   i=i+l
   endwhile
   y=shift(y,-m/2) ; shift back final signal (since response not centered)

; in case, interpolate in 'featureless' overlap regions.
; In those regions, small oscillations with an amplitude approx 0.002
; are created, due to the abrupt change of delta_nu.
; when overlap region is featureless, these become visible and can be
; avoided by interpolation   
   i=l
   while i le nn-1 do begin
     iminus=i-m/2
     iplus=min([i+m/2,nn-2])
     ym=y(iminus-1)
     yp=y(iplus+1)
     av=mean(y(iminus-1:iplus+1))
;     print,sxx(i),max(abs(y(iminus-1:iplus+1)-av))
     if max(abs(y(iminus-1:iplus+1)-av)) lt 0.005 then begin
       yv=[ym,yp]
       xv=[sxx(iminus-1),sxx(iplus+1)]
       y(iminus:iplus)=interpol(yv,xv,sxx(iminus:iplus))
;       print,'values around',sxx(i),' interpolated'
     endif
   i=i+l
   endwhile

   return
end
   
;--------------------------------------------------------------------------
   
pro prep_rot,sxx,i,il,m,px,py,shft

common convol_all, clight  
common share_rot, vrot, beta
  
   lambda0 = (sxx(i)+sxx(il))/2.
   dlam_max = vrot*(lambda0/clight)
   xx = px / dlam_max
   core = WHERE(ABS(xx) LE 1.)
;   print,i,m,n_elements(core),dlam_max
   if(n_elements(core) gt m) then begin
     print,'dimension of response function > M'
     stop
   endif  

   one = WHERE(ABS(xx) GT 1.)
   xx(one) = 1.
   xx_2 = xx*xx
   unsoeld =  (2.*     SQRT(1-xx_2)/!PI + (1-xx_2)*beta/2.) / (1 + 2*beta/3.)
;  gray = (2.*(1-beta)*SQRT(1-xx_2)/!PI + (1-xx_2)*beta/2.) / (1 -   beta/3.) 
   py = unsoeld
   py = py / dlam_max
; to apply overlap-add method, response function has to be located at the begin
; of the array (otherwise, overlap at both sides of segments)
; (leads to small error at most blueward points)
; use offset of broadest function
   py=shift(py,-shft)

   return
   
end

;--------------------------------------------------------------------------
   
pro prep_vmac,sxx,i,il,m,n,px,py,shft

common convol_all, clight  
common share_vmac, vmacro, ccr1, spi
  
   lambda0 = (sxx(i)+sxx(il))/2.
   mr = vmacro*(lambda0/clight)

   pxmr=abs(px)/mr

   core = WHERE(pxmr LE 6.)
;   print,i,m,n_elements(core),mr
   if(n_elements(core) gt m) then begin
     print,'dimension of response function gt M'
     stop
   endif  

; from resubstitution y=1/u and partial integration
;    py2t=cct*pxmt*[-spi+exp(-pxmt^2)/pxmt+spi*errorf(pxmt)]
;    py2r=ccr*pxmr*[-spi+exp(-pxmr^2)/pxmr+spi*errorf(pxmr)]

   py=fltarr(n)
   xx = pxmr(core)
   py(core)=ccr1/mr*[exp(-xx^2)+spi*xx*(errorf(xx)-1.)]

; to apply overlap-add method, response function has to be located at the begin
; of the array (otherwise, overlap at both sides of segments)
; (leads to small error at most blueward points)
; use offset of broadest function
   py=shift(py,-shft)
   return
   
end


;--------------------------------------------------------------------------
   
pro prep_resol,sxx,i,il,m,n,px,py,shft

common convol_all, clight  
common share_resol, vdop, fwhm, dlam_fwhm

   if fwhm then begin
     dlam_dop=dlam_fwhm
; wavelength independent
   endif else begin  
     lambda0 = (sxx(i)+sxx(il))/2.
     dlam_dop = vdop*(lambda0/clight)
   endelse
     
   pxmr=abs(px)/dlam_dop

   core = WHERE(pxmr LE 6.)
;   print,i,m,n_elements(core),mr
   if(n_elements(core) gt m) then begin
     print,'dimension of response function gt M'
     stop
   endif  

   py=fltarr(n)
   xx = pxmr(core)
   py(core)=exp(-xx^2)/(SQRT(!PI)*dlam_dop)

; to apply overlap-add method, response function has to be located at the begin
; of the array (otherwise, overlap at both sides of segments)
; (leads to small error at most blueward points)
; use offset of broadest function
   py=shift(py,-shft)
   return
   
end
