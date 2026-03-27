pro refrac, lam

;input  lam(vacuum) [Angst]
;output lam(air)    [Angst]

;only for lambda > 2000
b=where(lam gt 2000, count)
if count eq 0 then return

lam1=lam(b)

sig2=(10000.d0/lam1)^2
n=643.28d0+294981.d0/(146.d0-sig2)+2554.d0/(41.d0-sig2)
n=1.d0+n*1.d-7
sig2=double(lam1)
lam1=sig2/n
lam(b)=lam1

return
end
