pro read_iue,name,wave,flux
;
;to find out the content (name) of all entries, have a look into stars.dat
;name has to be a character, usually, but not always, the HD number
;example: read_iue,'66811' or read_iue,'Sk 80' 

openr,1,'~uh101aw/Observations/IUE/Walborn85/spectra.dat'

hd=''
while not eof(1) do begin
readf,1,format='(A10,2x,I4,4(1x,f9.5,2x,f4.2))',hd,wave1,flux1, $
       qual1,flux2,qual2,flux3,qual3,flux4,qual4
       if(strtrim(hd,2) eq name) then goto, start
endwhile
print,'entry not found'
close,1
return

start:
wave=fltarr(4*800)
flux=wave

wave(0)=wave1
wave(1)=wave1+0.25
wave(2)=wave1+0.5
wave(3)=wave1+0.75
flux(0)=flux1
flux(1)=flux2
flux(2)=flux3
flux(3)=flux4

for i=1,799 do begin
readf,1,format='(A10,2x,I4,4(1x,f9.5,2x,f4.2))',hd,wave1,flux1, $
       qual1,flux2,qual2,flux3,qual3,flux4,qual4
       if(strtrim(hd,2) ne name) then begin
        print,'something wrotten with number of entries for object'
	close,1
        stop
	return
       endif
       j=4*i
       wave(j)=wave1
       flux(j)=flux1
       j=j+1
       wave(j)=wave1
       flux(j)=flux1
       j=j+1
       wave(j)=wave1
       flux(j)=flux1
       j=j+1
       wave(j)=wave1
       flux(j)=flux1       
endfor
close,1

!p.multi=[0,1,3]
plot,wave,flux,xrange=[1000,1250],xs=1,yrange=[0,2]
plot,wave,flux,xrange=[1250,1500],xs=1,yrange=[0,2]
plot,wave,flux,xrange=[1500,1750],xs=1,yrange=[0,2]

!p.multi=0

return
end
