pro rd_fitsferos,filename,xx,yy,time
i=long(1.)
yy=readfits(filename,hd)
yy=yy*1.d0
nd=long(sxpar(hd,'naxis1'))
nplo=long(sxpar(hd,'naxis'))
xla1=double(sxpar(hd,'crval1'))
delt1=double(sxpar(hd,'cdelt1'))
;if keyword_set(time) then time=double(sxpar(hd,'exptime'))
time=double(sxpar(hd,'exptime'))
xx=dblarr(nd)
for i=0l,nd-1L do begin 
  xx(i)=xla1 + i * delt1
endfor
return
end


