pro profcolor,no_mod,col

  loadct,12
; standard
  col=[0,200,80,80,150,300,50,200]
;  col=[0,200,80,80,150,300,50,80]
;  col=[0,200,200,80,150,300,50,80]
  if no_mod gt 5 then begin
    stop,'More colors needed!'
  endif

  return
end  
