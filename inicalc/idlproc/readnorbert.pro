pro readnorbert,file,x,y,plot=plot
   readcol,file,format='F,F,F,F',x1,y1,x2,y2
   i1=n_elements(x1)
   i2=2*i1

   x=fltarr(i2)
   y=x
   for i=0,i1-1 do begin
     j=2*i
     j1=2*i+1
     x(j)=x1(i)
     y(j)=y1(i)
     x(j1)=x2(i)
     y(j1)=y2(i)
   endfor  
   
   if keyword_set(plot) then begin
     plot,x,y,xrange=[4000,5000]
   endif
   
return
end
