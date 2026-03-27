pro plot_range,model,resoltheo,ext,modelcomp=modelcomp,resol=resol, $
		vsini=vsini,vmacro=vmacro
  readcol,model+'/OUT.UV_00900_02000__'+resoltheo+'_'+ext,idum,lam1,xhdum1,prof1,xhdum1

  if keyword_set(modelcomp) then begin
    readcol,modelcomp+'/OUT.UV_00900_02000__'+resoltheo+'_'+ext,idum,lam1c,xhdum1,prof1c,xhdum1
  endif

  if not keyword_set(resol) then resol=10000.
  if not keyword_set(vsini) then vsini=100.
  if not keyword_set(vmacro) then vmacro=20.

  convol,lam1,prof1,l1,p1,vsini=vsini,vmac=vmacro,resol=resol
  lam1=l1
  prof1=p1
  
  if keyword_set(modelcomp) then begin
    convol,lam1c,prof1c,l1,p1,vsini=vsini,vmac=vmacro,resol=resol
    lam1c=l1
    prof1c=p1
  endif

  loadct,12
  window,1,xsize=1200,ysize=800
  !p.multi=[0,1,2]
  plot,lam1,prof1,xrange=[900,1500],xs=1,yrange=[0,2],ys=1
  if keyword_set(modelcomp) then oplot,lam1c,prof1c,col=200

  plot,lam1,prof1,xrange=[1500,2000],xs=1,yrange=[0,2],ys=1
  if keyword_set(modelcomp) then oplot,lam1c,prof1c,col=200

  readcol,model+'/OUT.UV_03400_07000__'+resoltheo+'_'+ext,idum,lam2,xhdum1,prof2,xhdum1

  if keyword_set(modelcomp) then begin
    readcol,modelcomp+'/OUT.UV_03400_07000__'+resoltheo+'_'+ext,idum,lam2c,xhdum1,prof2c,xhdum1
  endif

  convol,lam2,prof2,l1,p1,vsini=vsini,vmac=vmacro,resol=resol
  lam2=l1
  prof2=p1
  
  if keyword_set(modelcomp) then begin
    convol,lam2c,prof2c,l1,p1,vsini=vsini,vmac=vmacro,resol=resol
    lam2c=l1
    prof2c=p1
  endif

  window,0,xsize=1200,ysize=800

  !p.multi=[0,1,3]

  y1=min(prof2)
  y2=max(prof2)
  plot,lam2,prof2,xrange=[3400,5200],xs=1,yrange=[y1,y2],ys=1
  if keyword_set(modelcomp) then oplot,lam2c,prof2c,col=200

  plot,lam2,prof2,xrange=[5200,7000],xs=1,yrange=[y1,y2],ys=1
  if keyword_set(modelcomp) then oplot,lam2c,prof2c,col=200

  plot,lam2,prof2,xrange=[4600,4750],xs=1,yrange=[y1,y2],ys=1
  if keyword_set(modelcomp) then oplot,lam2c,prof2c,col=200


  
  !p.multi=0
return
end
