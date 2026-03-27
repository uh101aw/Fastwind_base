pro irprof_new,mods,star,ext=extin, $
    vsini=vsini,vmacro=vmacro,resol=resol, $
    obscomp=obscomp, vrad=vrad, vrl=vrl, vacuum=vacuum,$
    lamcomp=lamcomp,profcomp=profcomp, $
    lamc1=lamc1,profc1=profc1, $
    layout=layout,oldnames=oldnames, $
    nocolor=nocolor,ps=ps, $
    comp=comp,e1=e1,e2=e2, $
    lb_resol=lb_resol	      

if not keyword_set(ps) then window,0,xsize=1600,ysize=800
  
if keyword_set(lb_resol) then begin
resol_lband=lb_resol
endif else begin
;resol_lband=6700. ;VLT-ISAAC
resol_lband=2500. ; SPEX
endelse

if keyword_set(nocolor) then begin
  color=0
endif else begin
  color=1
endelse  

models=mods    
; last three keywords for compatibility with older versions
if keyword_set(comp) then begin
  dim1=(size(models,/dim))(0)
  dim2=(size(comp,/dim))(0)
  if dim1 ne 0 or dim2 ne 0 then begin
    print,'Inconsistent input for models and/or comp!!!'
    return
  endif  
  models=[models,comp]    
endif

if keyword_set(e1) then begin
  dim1=(size(e1,/dim))(0)
  if dim1 ne 0 then begin
    print,'Inconsistent input for e1!!!'
    return
  endif  
  extin=e1    
endif

if keyword_set(e2) then begin
  dim2=(size(e2,/dim))(0)
  if dim2 ne 0 or n_elements(extin) eq '' then begin
    print,'Inconsistent input for e2!!!'
    return
  endif  
  extin=[extin,e2]    
endif

;maximum data set (Nov. 2022): 24 lines
lines1=['PGAMMA','PBETA','PALPHA', $
       'BR14','BR13','BR12','BR11','BR10', $
       'BRGAMMA','BRBETA','BRALPHA', $
       'PF10','PF9','PFGAMMA', $
       'HEI170','HEI205','HEI211','HEI370', $
       'HEII57','HEII611','HEII713','HEII712','HEII218','HEII67']

npmax=n_elements(lines1)
if npmax ne 24 then stop,' more than 24 lines'

titles1=['Pa_gamma','Pa_beta','Pa_alpha', $
	'Br14','Br13','Br12','Br11','Br10', $
	'Br_gamma','Br_beta','Br_alpha', $
        'Pf_10','Pf_9','Pf_gamma', $
	'HeI1.70','HeI2.05','HeI2.11','HeI3.70', $
	'HeII5_7','HeII6_11','HeII7_13','HeII1.69','HeII2.18','HEII6_7']

if n_elements(titles1) ne npmax then stop,'dim(lines1) ne dim(titles1)'

xranges1=fltarr(2,npmax); set to zero for all non-specified line-ranges

xranges1[*,3]=[1.578,1.598]    ;Br14
xranges1[*,4]=[1.6009,1.6209]  ;Br13
xranges1[*,5]=[1.63,1.65]      ;Br12
xranges1[*,6]=[1.67,1.69]      ;Br11
xranges1[*,7]=[1.726,1.746]    ;Br10
xranges1[*,8]=[2.15,2.18]      ;Br_gamma
xranges1[*,10]=[4.02,4.08]     ;Br_alpha 
xranges1[*,13]=[3.7,3.78]      ;Pf_gamma
xranges1[*,14]=[1.695,1.706]   ;HeI 1.7
xranges1[*,15]=[2.05,2.065]    ;HeI 2.05
xranges1[*,16]=[2.107,2.118]   ;HeI 2.11
xranges1[*,21]=[1.688,1.697]   ;HeII 1.69
xranges1[*,22]=[2.183,2.195]   ;HeII 2.18

nspec=24
if keyword_set(layout) then nspec=layout
;different layouts
case nspec of
24: begin
    lines=['PGAMMA','PBETA','PALPHA', $
       'BR14','BR13','BR12','BR11','BR10', $
       'BRGAMMA','BRBETA','BRALPHA', $
       'PF10','PF9','PFGAMMA', $
       'HEI170','HEI205','HEI211','HEI370', $
       'HEII57','HEII611','HEII713','HEII712','HEII218','HEII67']
    !p.multi=[0,8,3]
    end

21: begin
lines=['PGAMMA','PBETA','PALPHA', $
       'BR12','BR11','BR10', $
       'BRGAMMA','BRBETA','BRALPHA', $
       'PF10','PF9','PFGAMMA', $
       'HEI170','HEI205','HEI211', $
       'HEII57','HEII611','HEII713','HEII712','HEII218','HEII67']
    !p.multi=[0,7,3]
    end

12: begin
    lines=['BR12','BR11','BR10','BRGAMMA','PFGAMMA','BRALPHA', $
       'HEI170','HEI205','HEI211','HEI370','HEII712','HEII218']
    !p.multi=[0,4,3]
    end

8:  begin
    lines=['BR11','BR10','BRGAMMA', $
        'HEI170','HEI205','HEI211','HEII712','HEII218']
    !p.multi=[0,4,2]
    end

6:  begin
    lines=['BR14','BR13','BR12','BR11','BR10','BRGAMMA']
    !p.multi=[0,3,2]
    end

1:  begin
;    lines=['BRGAMMA']
    lines=['BRALPHA']
    !p.multi=0
    end

else: stop,'wrong layout, allowed are 24,21,12,8,6,1'
endcase  

np=n_elements(lines)
titles=strarr(np)
xranges=fltarr(2,np)

for i=0,n_elements(lines)-1 do begin
  j=where(lines(i) eq lines1,count)
  if count ne 1 then stop,' lines not found in lines1'
  titles(i)=titles1(j)
  xranges(*,i)=xranges1(*,j)
endfor

if keyword_set(oldnames) then begin
  ind=where(lines eq 'BR10')
  lines(ind)='BRZETA'
endif

;if !p.charsize eq 1. then !p.charsize=1.5
!p.charsize=2.

no_lin=(size(lines))(1)
no_tit=(size(lines))(1)
if no_tit ne no_lin then begin
  print,'no_lin ne no_tit!' 
  return
endif  

; make and check arrays for model(s) and extension(s)
cat=models  
no_mod=(size(cat,/dimensions))(0)
if no_mod eq 0 then cat=[cat]

if keyword_set(extin) then begin
  ext=extin
  no_ext=(size(ext,/dimensions))(0)
  if no_ext ne no_mod then begin
    print,'Dimension of models and extensions not consistent!'
    return
  endif
  if no_ext eq 0 then ext=[ext]
endif

if no_mod eq 0 then no_mod=1

; color definition
if keyword_set(color) then begin
  loadct,12
  col=[0,200,100,80,150]
  if no_mod gt 4 then begin
    print,'More colors needed!'
    return
  endif  
endif else begin
  col=[-1] 
endelse


; modify extensions so that appropriate
if not keyword_set(extin) then begin
  ext=strarr(no_mod)
  for i=0,no_mod-1 do begin
    print,'model ',cat(i)
    print,'Which extensions shall be processed?'
    print,'standard [RET]  --  _ESC [e]  --  _VTxxx [xxx]  --  _ESC_VTxxx [exxx]'
    extension=''
    read,extension
    extension=strtrim(extension,2)
    dummy=strlen(extension)
    if (dummy eq 1) then extension='_ESC'
    if (dummy eq 3) then extension='_VT'+extension
    if (dummy eq 4) then extension='_ESC_VT'+strmid(extension,1,3)
    ext(i)=extension
  endfor
  endif else begin
  for i=0,no_mod-1 do begin
  if ext(i) eq '0' then begin
  extension=''
  endif else begin
  extension='_'+strtrim(ext(i),2)  
  endelse
  ext(i)=extension 
  endfor
endelse

print,cat(0)
rindat,cat(0),teff
if teff gt 35000. then begin
  ind=where(lines eq 'BRGAMMA')
  lines(ind)='BRGAMMA1'
  file=cat(0)+'/OUT.BRGAMMA1'+ext(0)
  if not file_test(file) then lines(ind)='BRGAMMA'
  ind=where(lines eq 'BRALPHA')
  lines(ind)='BRALPHA1'
  file=cat(0)+'/OUT.BRALPHA1'+ext(0)
  if not file_test(file) then lines(ind)='BRALPHA'
endif 

; observed/comparison profiles
if keyword_set(obscomp) then begin
  if(obscomp eq 'oplot') then begin
   wobs=lamcomp
   pobs=profcomp
  endif else begin

;  obsfile='$HOME/Observations/IR/HiRes/'+obscomp
;  obsfile='$HOME/Observations/IR/B-SG/'+obscomp
;  obsfile='$HOME/Observations/IR/GalCenOB/'+obscomp+'.dat'
  obsfile='$HOME/Observations/IR/O-SG/combined/'+obscomp
  
  if not file_test(obsfile) then begin
    print,obsfile,' does not exist'
    return
  endif  
  rspec,obsfile,wobs,pobs
  wobs=wobs/10000.
  endelse

  if keyword_set(vacuum) then begin
    if not keyword_set(lamcomp) then begin
      print,' keyword vacuum only in connection with "oplot"!'
      return
    endif  
    refrac,wobs
    wobs=wobs/10000.
  endif  

  if keyword_set(vrad) then wobs=wobs*(1.-vrad/2.99792e5)
  wobs1=wobs
  if keyword_set(vrl)  then wobs1=wobs*(1.-vrl/2.99792e5)

endif else begin
   wobs=[0.]
   pobs=[0.]
   wobs1=[0.]
endelse 

dat={star:'',  $
     cat:strarr(no_mod), ext:strarr(no_mod), $
     vsini:0., vmacro:0., resol:0., $
     lines: strarr(no_lin), titles: strarr(no_lin)}

if nspec eq 12 then begin
  nobs=(size(wobs1))(1)
  lband={resol: 0., wobs1: fltarr(nobs)}
  lband.resol=resol_lband
  lband.wobs1=wobs1
endif else begin
  lband=0
endelse  

if n_params() ne 1 then dat.star=star+': '
dat.cat=cat
dat.ext=ext
dat.lines=lines
dat.titles=titles

if keyword_set(vsini)  then dat.vsini=vsini
if keyword_set(vmacro) then dat.vmacro=vmacro
if keyword_set(resol)  then dat.resol=resol


if keyword_set(ps) then begin
set_plot,'ps'
;print,' Give in Output-filename'
outfile=ps
;read,outfile
if (nspec ne 1) then begin 
device,filename=outfile,ysize=20.,xoffset=1.,yoffset=6.5,/color
endif else begin
device,filename=outfile,/color,/landscape
endelse  
endif

; plot all profiles
;plot_profiles_old, dat, wobs, pobs, col, xranges=xranges,/ir,lband=lband
plot_profiles, dat, wobs, pobs, col, xranges=xranges, lamc1=lamc1, profc1=profc1,/ir;,lband=lband

if keyword_set(fcl) then begin
  readmod,dat.cat(0)+'/OUT_TOT',' ',a
  maxclf=max(a.fcl)
  for i=1,no_mod-1 do begin
    readmod,dat.cat(i)+'/OUT_TOT',' ',b
    maxclf=max([maxclf,b.fcl])
  endfor
  !y.range=[0,maxclf+1]  
  plfcl,a,/plv
  for i=1,no_mod-1 do begin
    readmod,dat.cat(i)+'/OUT_TOT',' ',a
    if keyword_set(color) then begin
      oplot,a.v,a.fcl,color=col(i)
    endif else begin
      oplot,a.v,a.fcl,line=i
    endelse
  endfor
endif  

if keyword_set(color) then loadct,0
!p.title=''
!p.charsize=1.
!p.multi=0
!x.range=0.
!y.range=0.

if keyword_set(ps) then begin
device,/close
set_plot,'x'
endif else begin
char=''
read,char
wdelete,0
endelse

return
end
