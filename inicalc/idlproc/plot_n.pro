pro plot_n,models,star,obs=obs,ext=extin,vrad=vrad,vsini=vsini,vmacro=vmacro,resol=resol,hcopy=hcopy, $
	    lamcomp=lamcomp,profcomp=profcomp,lamc1=lamc1,profc1=profc1


if keyword_set(hcopy) then begin

!x.thick=4.
!y.thick=4.
!p.charthick=4.
!p.thick=4.

set_plot,'ps'
if star eq '' or star eq ' ' then begin
  psfile='idl.ps'
endif else begin
  psfile=star+'.ps'
endelse
  
device,file=psfile,/color,xsize=19.7,ysize=25.,xoff=0.,yoff=0.

hprof,models,star,ext=extin,obscomp=obs,vrad=vrad,vsini=vsini,vmacro=vmacro,resol=resol,/color,/fcl, $
	    lamcomp=lamcomp,profcomp=profcomp,lamc1=lamc1,profc1=profc1
       
heiprof,models,star,ext=extin,obs=obs,vrad=vrad,vsini=vsini,vmacro=vmacro,resol=resol,/color, $
	 	    lamcomp=lamcomp,profcomp=profcomp,lamc1=lamc1,profc1=profc1

heiiprof,models,star,ext=extin,obs=obs,vrad=vrad,vsini=vsini,vmacro=vmacro,resol=resol,/color, $
	  	    lamcomp=lamcomp,profcomp=profcomp,lamc1=lamc1,profc1=profc1


niiprof,models,star,ext=extin,obscomp=obs,vrad=vrad,vsini=vsini,vmacro=vmacro,resol=resol,/color, $
	    lamcomp=lamcomp,profcomp=profcomp,lamc1=lamc1,profc1=profc1

ntriplet,models,star,obs=obs,ext=extin,vrad=vrad,vsini=vsini,vmacro=vmacro,resol=resol,hcopy=hcopy, $
	  	    lamcomp=lamcomp,profcomp=profcomp,lamc1=lamc1,profc1=profc1

nivprof,models,star,ext=extin,obscomp=obs,vrad=vrad,vsini=vsini,vmacro=vmacro,resol=resol,/color, $
	 	    lamcomp=lamcomp,profcomp=profcomp,lamc1=lamc1,profc1=profc1

nvprof,models,star,ext=extin,obscomp=obs,vrad=vrad,vsini=vsini,vmacro=vmacro,resol=resol,/color, $
		    lamcomp=lamcomp,profcomp=profcomp,lamc1=lamc1,profc1=profc1



device,/close
set_plot,'x'

!x.thick=1.
!y.thick=1.
!p.charthick=1.
!p.thick=1.
endif else begin

;-------------------------------------

window,2,xsize=800,ysize=800
niiprof,models,star,ext=extin,obscomp=obs,vrad=vrad,vsini=vsini,vmacro=vmacro,resol=resol,/color, $
  	    lamcomp=lamcomp,profcomp=profcomp,lamc1=lamc1,profc1=profc1


; window 4 and 6
ntriplet,models,star,obs=obs,ext=extin,vrad=vrad,vsini=vsini,vmacro=vmacro,resol=resol,hcopy=hcopy, $
	  	    lamcomp=lamcomp,profcomp=profcomp,lamc1=lamc1,profc1=profc1


window,5,xsize=800,ysize=800
nivprof,models,star,ext=extin,obscomp=obs,vrad=vrad,vsini=vsini,vmacro=vmacro,resol=resol,/color, $
	 	    lamcomp=lamcomp,profcomp=profcomp,lamc1=lamc1,profc1=profc1


window,7,xsize=800,ysize=800
nvprof,models,star,ext=extin,obscomp=obs,vrad=vrad,vsini=vsini,vmacro=vmacro,resol=resol,/color, $
		    lamcomp=lamcomp,profcomp=profcomp,lamc1=lamc1,profc1=profc1


;window,0
;hprof,models,star,ext=extin,obscomp=obs,vrad=vrad,vsini=vsini,vmacro=vmacro,resol=resol,/color,/fcl, $
;	    lamcomp=lamcomp,profcomp=profcomp,lamc1=lamc1,profc1=profc1


;window,1
;heiprof,models,star,ext=extin,obs=obs,vrad=vrad,vsini=vsini,vmacro=vmacro,resol=resol,/color, $
;	    lamcomp=lamcomp,profcomp=profcomp,lamc1=lamc1,profc1=profc1


;window,3
;heiiprof,models,star,ext=extin,obs=obs,vrad=vrad,vsini=vsini,vmacro=vmacro,resol=resol,/color.
;	    lamcomp=lamcomp,profcomp=profcomp,lamc1=lamc1,profc1=profc1



a=' '
read,a

wdelete,2
wdelete,4
wdelete,6
wdelete,5
wdelete,7
;wdelete,0
;wdelete,1
;wdelete,3

endelse

return
end
