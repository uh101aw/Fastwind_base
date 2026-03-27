pro ntriplet,models,star,obs=obs,ext=extin,vrad=vrad,vsini=vsini,vmacro=vmacro,resol=resol,hcopy=hcopy, $
	      lamcomp=lamcomp,profcomp=profcomp,lamc1=lamc1,profc1=profc1

window,4,xsize=800,ysize=800
;niiiprof,models,star,ext=extin,obscomp=obs,vrad=vrad,vsini=vsini,vmacro=vmacro,resol=resol,restn=1,/color, $
niiiprof,models,star,ext=extin,obscomp=obs,vrad=vrad,vsini=vsini,vmacro=vmacro,resol=resol,restk=1,/color, $
	  	      lamcomp=lamcomp,profcomp=profcomp,lamc1=lamc1,profc1=profc1

window,6,xsize=800,ysize=800
niiiprof,models,star,ext=extin,obscomp=obs,vrad=vrad,vsini=vsini,vmacro=vmacro,resol=resol,trip=1,/color, $
	  	      lamcomp=lamcomp,profcomp=profcomp,lamc1=lamc1,profc1=profc1
return
end
