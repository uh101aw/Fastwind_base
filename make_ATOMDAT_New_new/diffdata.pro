pro diffdata, type

; input types: '*.lev','*.nldr','*.col','*.pho','*.nl3','*.info'
  
    
  files=file_search('fw_data/'+type)
  
  nf=n_elements(files)
  var=strarr(2*nf)
  
  for i=0,nf-1 do begin
    var(2*i:2*i+1) = strsplit(files(i),'/',/extract)
  endfor

  files=var(1:2*nf-1:2)

  nf=n_elements(files)
;  dir1='/home/helix/hoffmann/work/sci/codes/Data-Source-Utils/Atomic-Models-Source/'
  dir1='tadziu_data/'
  for i = 0,nf-1 do begin
    print,files(i)
; -Z: neglect trailing blancs (in *.lev for labname)    
    spawn,['diff','-Z','-y','--suppress-common-lines','fw_data/'+files(i),dir1+files(i)],/noshell
;    spawn,'diff -s fw_data/'+files(i)+' '+dir1+files(i)   
    
  endfor

  print,'comparison for ',nf,' files completed' 

  
return

end
