pro rindat1,file,model
  openr,1,file
  a=' '
  readf,1,a
  stri=str_sep(a,' ')
  b=where(stri ne '')
  str=stri(b)
  model=str(0)
  
  readf,1,a
  readf,1,a
  readf,1,a
  
  stri=str_sep(a,' ')
  b=where(stri ne '')
  str=stri(b)
  teff=float(str(0))
  logg=float(str(1))
  rstar=float(str(2))
  
  readf,1,a
  stri=str_sep(a,' ')
  b=where(stri ne '')
  str=stri(b)
  tmin=float(str(1))

  readf,1,a
  stri=str_sep(a,' ')
  b=where(stri ne '')
  str=stri(b)
  mdot=float(str(0))
  vinf=float(str(2))
  beta=float(str(3))

  readf,1,a
  stri=str_sep(a,' ')
  b=where(stri ne '')
  str=stri(b)
  yhe=float(str(0))

  readf,1,a
  readf,1,a
  stri=str_sep(a,' ')
  b=where(stri ne '')
  str=stri(b)
  vturb=float(str(0))
  z=float(str(1))

  q=alog10(mdot/(rstar*vinf)^1.5)
  char=['A','B','C','D','E','F','G','H','larger than H']
  qgrid=[-14.,-13.5,-13.15,-12.8,-12.45,-12.1,-11.75,-11.4]
  qgridm=qgrid
  for i=1,7 do begin
    qgridm(i-1)=qgrid(i)+.5*(qgrid(i-1)-qgrid(i))
  endfor
  
  for i=0,7 do begin
    if q le qgridm(i) then goto, lab
  endfor
  i=8

  lab: qchar=char(i)
  print
  print,' model= ',model
  print
  print,' Teff = ',strtrim(string(teff),2),' , logg = ',strtrim(string(logg),2), $
	 ' , YHe =',strtrim(string(yhe),2),' , Rstar = ',strtrim(string(rstar),2)
  print,' ',qchar,': log Q = ',strtrim(string(q),2)
  print,' vturb = ',strtrim(string(vturb),2),' , z = ',strtrim(string(z),2), $
	 ' , Tmin = ',strtrim(string(tmin),2)	 
  print,' Mdot = ',strtrim(string(mdot),2),' , vinf = ',strtrim(string(vinf),2), $
	 ' , beta = ',strtrim(string(beta),2)  
  print
  
  close,1
  return
end  
