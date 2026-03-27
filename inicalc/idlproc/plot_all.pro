pro plot_all,models,star,obs=obs,ext=extin,vrad=vrad,vsini=vsini,vmacro=vmacro,resol=resol,hcopy=hcopy, $
	      lamcomp=lamcomp,profcomp=profcomp,lamc1=lamc1,profc1=profc1


if keyword_set(obs) then begin
  if obs ne 'oplot' then begin
; test for complete file_path
    print,file_test(obs)
    if file_test(obs) then goto, cont
;  obsfile='$HOME/Observations/optical/'+obs
;  obsfile='$HOME/Observations/optical/Gal_Bsg_Nevy/'+obs+'.sp'
  if file_test(obsfile) then begin
    obs=obsfile
    goto, cont
  endif else begin
    print,obsfile,' does not exist'
    return
  endelse
  endif  
endif
  
cont:  
if keyword_set(hcopy) then begin

!x.thick=4.
!y.thick=4.
!p.charthick=4.
!p.thick=4.

set_plot,'ps'
device,file=star+'.ps',/color,xsize=19.7,ysize=25.,xoff=0.,yoff=0.

hprof,models,star,ext=extin,obscomp=obs,vrad=vrad,vsini=vsini,vmacro=vmacro,resol=resol,/color,/fcl, $
       	      lamcomp=lamcomp,profcomp=profcomp,lamc1=lamc1,profc1=profc1

;heiprof,models,star,ext=extin,obs=obs,vrad=vrad,vsini=vsini,vmacro=vmacro,resol=resol,/color, $
;	 	      lamcomp=lamcomp,profcomp=profcomp,lamc1=lamc1,profc1=profc1
heiprof,models,star,ext=extin,obs=obs,vrad=vrad,vsini=vsini,vmacro=vmacro,resol=resol,/color, $
	 	      lamcomp=lamcomp,profcomp=profcomp,lamc1=lamc1,profc1=profc1,layout=8

heiiprof,models,star,ext=extin,obs=obs,vrad=vrad,vsini=vsini,vmacro=vmacro,resol=resol,/color, $
	  	      lamcomp=lamcomp,profcomp=profcomp,lamc1=lamc1,profc1=profc1


device,/close
set_plot,'x'

!x.thick=1.
!y.thick=1.
!p.charthick=1.
!p.thick=1.
endif

;-------------------------------------

window,0
hprof,models,star,ext=extin,obscomp=obs,vrad=vrad,vsini=vsini,vmacro=vmacro,resol=resol,/color,/fcl, $
       	      lamcomp=lamcomp,profcomp=profcomp,lamc1=lamc1,profc1=profc1

window,1
;heiprof,models,star,ext=extin,obs=obs,vrad=vrad,vsini=vsini,vmacro=vmacro,resol=resol,/color, $
;	 	      lamcomp=lamcomp,profcomp=profcomp,lamc1=lamc1,profc1=profc1
heiprof,models,star,ext=extin,obs=obs,vrad=vrad,vsini=vsini,vmacro=vmacro,resol=resol,/color, $
	 	      lamcomp=lamcomp,profcomp=profcomp,lamc1=lamc1,profc1=profc1,layout=8

window,2
heiiprof,models,star,ext=extin,obs=obs,vrad=vrad,vsini=vsini,vmacro=vmacro,resol=resol,/color, $
	  	      lamcomp=lamcomp,profcomp=profcomp,lamc1=lamc1,profc1=profc1

a=' '
read,a
wdelete,0
wdelete,1
wdelete,2
return
end
