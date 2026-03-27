pro plifrac_met,file,dat,metal,comp=dat1,cmod=cmod, $
		 atom=atom,imin=imin,imax=imax,taur=taur

loadct,12

char=' '
lab: 	print,' give in x=axis: m or taur'
if not keyword_set(taur) then begin
  read,char
endif else begin
  char='taur'
endelse

if char eq 'taur' then begin
x=dat.taur
xtit='tau_Ross'
endif else if char eq 'm' then begin
x=dat.m
xtit='m'
endif else begin
print,' Wrong Input'
goto, lab
endelse

ndim=size(x)
nd=ndim(1)

file=file+'/METAL_IDL'
openr,1,file

metal=dblarr(10,30,nd)
readf,1,metal
close,1

lab0: print,' '
print,' give in No. of atom'
if not keyword_set(atom) then begin
  read,na
endif else begin
  na=atom
endelse  
lab1: print,' '
print,' give in min, max  ionization stage (astron. convention)'
if not keyword_set(atom) then begin
  read,imin,imax
endif 
  
if imax gt 9 then begin
print,' Wrong input'
goto, lab1
endif

xmin=x(1)

y=metal(imin,na-1,*)
ymin=min(y)

;for specific tests
;y1=metal(imin:imax,na-1,*)
;ymin=min(y1)

plot_oo,x,y,title=dat.name,xtitle=xtit,ytitle='ioniz. fraction', $
xrange=[x(nd-1),xmin],yrange=[ymin,2]
;plot_oi,x,y,title=dat.name,xtitle=xtit,ytitle='ioniz. fraction', $
;xrange=[x(nd-1),xmin],yrange=[0,1.1]

lines=0
for i=imin+1,imax do begin
lines=lines+1
if lines gt 5 then lines=0
y=metal(i,na-1,*)
oplot,x,y,linestyle=lines
endfor

if keyword_set(dat1) then begin
file1=''
if keyword_set(cmod) then begin
  file1=cmod
endif else begin  
  print,'input comparison model'
  read,file1
endelse
  
file1=file1+'/METAL_IDL'
openr,1,file1

metal1=fltarr(10,30,nd)
readf,1,metal1
close,1

if char eq 'taur' then begin
x1=dat1.taur
endif else if char eq 'm' then begin
x1=dat1.m
endif

y=metal1(imin,na-1,*)
ymin=min(y)
oplot,x1,y,color=50

lines=0
for i=imin+1,imax do begin
lines=lines+1
if lines gt 5 then lines=0
y=metal1(i,na-1,*)
oplot,x1,y,linestyle=lines,col=50
endfor
endif

if keyword_set(taur) then return

lab3: print,' '
print,' new plot = 0 '
print,' end      = 2 '
read,ov
if ov eq 2 then return
goto, lab0

return
end
