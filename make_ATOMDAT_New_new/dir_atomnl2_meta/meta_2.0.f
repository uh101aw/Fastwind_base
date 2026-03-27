      program metast
*----------------------------------------------------------------------
*     identifies metastable levels from Adi Pauldrach's packed
*     level list 'atomlev' and writes new file 'atomnl_meta'
*     which contains all level information needed for program 
*     'dice3.f'
*
*     For statistical purposes, selected lines can be printed
*     in SR lines (file 'testlist'). This file can be further
*     processed by program 'linstat.f'.
*
*     AUTHOR: 
*        U. Springmann, 03/1996
*     
*     INPUT:
*        files in fpath, atomic data files
*        includefile 'metacom', common blocks
*
*     OUTPUT: 
*        file 'atomnl2_meta_v2.0', file 'testlist'
*        recent changes in atomnl2_meta_v2.0
*        for all ' ' levels which are coupled to the ground-state multiplet,
*        we provide also the parent-level,the Aji value and the transition
*        frequency for the strongest transition. 
*
*     COMPILATION:
*        if subroutine 'lines' is used, the program 
*        HAS TO BE COMPILED WITH OPTION -align dcommons
*
*     HISTORY:
*        03/96  correction for level designations
*               modular design for different puposes
*               preprocessing line list
*        05/96  different file paths (in metacom)
*        11/96  extract relevant lines in SR lines -> much smaller
*               data set
*        09/02  changed by Jo: block data and file pathes (-> old)
*               modified bug with respect to forbidden transitions
*               (no separation of different spin systems, ->
*                more meta-stable levels, less transitions to ground state,
*                more transitions to meta-stable states
*                compare, eg., NIII and NIV)
*     version:1.0   12/04  uses Aji instead of gf to define 's' vs. ' ' levels
*     version:2.0   include upper level of bf transitions (new format of atomnl2_meta_v2.0)
*                   (note that term design. of V has changed in atomlev)
*                   for bf transitions with more than one upper level, we use
*                   the one with the highest cross-section as the representative
*                   one, with alpha as the sum of all possible transitions
*----------------------------------------------------------------------
      include 'metacom'

      call atomin               ! read in atomic data
      call sortlev              ! sort levels for upper level
      call ident                ! identifies metastable levels
      call outnl2               ! write atomnl2_meta_v2.0
      call lines                ! write out line list 
!      call forbout              ! output of forbidden transitions

      end

      subroutine atomin
*----------------------------------------------------------------------
*     read atomic data
*----------------------------------------------------------------------
      include 'metacom'
      parameter (nmx=natom*nion*50)

      open( 9,file=fpath//'atomlev',status='old',err=100,iostat=ios)
      open(10,file=fpath//'atomgf3',status='old',err=101,iostat=ios)
      open(11,file=fpath//'generalinfo',status='old',err=102,
     &     iostat=ios)
      open(12,file=fpath//'atomnl2',status='old',err=103,iostat=ios)
      open(15,file=fpath//'ionisen',status='old',err=104,iostat=ios)
      open(16,file=fpath//'atompho',status='old',err=105,iostat=ios)
      open(17,file=fpath//'partit',status='old',err=106,iostat=ios)

* read in general information
* ---------------------------

      do i=1,nrec
         read(11,*)irec,kel(i),iel(i),idum,ilev(i),levphi(i),ilin(i),
     +             idum,idum,ielev(i)
      enddo

* read stat. weight for last ionis. stage and number of stages per atom
* ---------------------------------------------------------------------

      read(17,*)(jmax(i),i=1,natom)
      read(17,*)(upart(i),i=1,natom)


* read ionisation energies
* ------------------------

      do i=1,natom
         read(15,*)(ergion(i,j),j=1,jmax(i))
      enddo

* read level energies
* ------------------------

      do i=1,nrec
        k=kel(i)
        j=iel(i)
        read(12,*) id1,id2
        if(id2.ne.ielev(i)) stop' error in ielev'
        do jj=1,ielev(i)
           read(12,*) ener
           energy(k,j,jj)=ener
        enddo
      enddo
      rewind(12)

* read level term designations from file 'atomlev'
* -----------------------------------------------

      ien=0

      do i=1,nrec
         if(ilev(i) .gt. 0)then
            read(9,*)n,level

            if(level .ne. ilev(i))then
               write(*,*) 'error in atomlev', i,level,ilev(i)
            endif

            do k=1,ilev(i)
               read(9,20)l,atom(ien+l)
            enddo
         endif
         ien=ien+ielev(i)
      enddo
 20   format(i10,2x,a11)

* read level transition information from file 'atomgf3'
* -----------------------------------------------------
c#a
      gfcut = 1.e-4
c      print*,' gf_cut?'
c      read*,gfcut
c#e
      j = 0
      ien=0
      n=1

 50   continue
         read(10,*) nr,lines
 
         if(lines .ne. ilin(nr))then
            write(*,*) 'error in atomgf3:', nr,lines,ilin(nr)
         endif

* for missing transition information, update cum. level index
* (currently element 23I ... 23III)
* -----------------------------------------------------------

         if(nr .gt. n)then
            do l=n,nr-1
               ien=ien+ielev(l)
            enddo
            n=nr
         endif

         do i=1,ilin(nr)
            j = j+1
            read(10,*) ind(j),wave(j),gf(j)
* reverse lower and upper level
* -----------------------------

            npart = ind(j)/10000
            nrest = ind(j) - npart*10000
            low = nrest/100
            lup = nrest - 100*low

* test for forbidden transitions 
* ------------------------------
* either all transitions with gf < gfcut are considered forbidden,
* or test according to LS-coupling rules (also serves to test for
* incorrect atomic level designations)

            if(gf(j) .lt. gfcut)then
            call forbid(atom,nrec,ilev,n,ien,low,lup,ndec)
c#a
c               ndec=1
c#e
               if(ndec.eq.1)then
                  forb(j)='f'
               endif
            endif

            kind(j) = npart*10000 + lup*100 + low
         enddo
         ien=ien+ielev(n)
         n=n+1
      if(n.le.nrec)goto 50
      number = j

      return
 100  write(*,*) ' error in opening atomlev, iostat= ',ios
      stop
 101  write(*,*) ' error in opening atomgf3, iostat= ',ios
      stop
 102  write(*,*) ' error in opening generalinfo, iostat= ',ios
      stop
 103  write(*,*) ' error in opening atomnl2, iostat= ',ios
      stop
 104  write(*,*) ' error in opening ionisen, iostat= ',ios
      stop
 105  write(*,*) ' error in opening atompho, iostat= ',ios
      stop
 106  write(*,*) ' error in opening partit, iostat= ',ios
      stop
      end

      block data init
      include 'metacom'
      parameter (nmx=natom*nion*50)

      data name/' H', 'He', 'Li', 'Be', ' B', ' C', ' N', ' O', ' F', 
     +          'Ne', 'Na', 'Mg', 'Al', 'Si', ' P', ' S', 'Cl', 'Ar',
     +          ' K', 'Ca', 'Sc', 'Ti', ' V', 'Cr', 'Mn', 'Fe', 'Co', 
     +          'Ni', 'Cu', 'Zn'/,
     +    stage/'_I  ', '_II ', '_III', '_IV ', '_V  ', '_VI ', '_VII',
     +          'VIII', '_IX ', '_X  ', '_XI '/

      data meta/5000*'m'/,forb/nlev*' '/,source/5000*' '/,
     +     cross/100*0./,gcross/100*0./,ground/5000*0./,
     +     ocross/5000*0./,konlev/nmx*0/,ilow/5000*0/,
     +     atom/5000*' '/,metall/nmx*'m'/,ipop/5000*0/,
     +     maxgf/5000*0./,gs_multi/5000*0/,ajimax/5000*0./,
     +     maxajisub/5000*0./,ilowsub/5000*0/

      end

      subroutine sortlev
*----------------------------------------------------------------------
*     sort for upper level; make corresponding changes in wave, gf
*----------------------------------------------------------------------
      include 'metacom'
      character rforb(nlev)*1
      integer wksp(nlev),iwksp(nlev)
      real rwksp(nlev)
      
      data wksp/nlev*0/,iwksp/nlev*0/

      call indexx(number,kind,iwksp)
* resorting for upper level: ind -> kind
      do i=1,number
         wksp(i)=kind(i)
      enddo
   
      do i=1,number
         kind(i)=wksp(iwksp(i))
      enddo
* resorting wave
      do i=1,number
         rwksp(i)=wave(i)
      enddo
      do i=1,number
         wave(i)=rwksp(iwksp(i))
      enddo
* resorting gf
      do i=1,number
         rwksp(i)=gf(i)
      enddo
      do i=1,number
         gf(i)=rwksp(iwksp(i))
      enddo
* resorting forb
      do i=1,number
         rforb(i)=forb(i)
      enddo
      do i=1,number
         forb(i)=rforb(iwksp(i))
      enddo

      end

      subroutine ident
*----------------------------------------------------------------------
*     make connection to levels in 'atomnl2' and identify transitions
*----------------------------------------------------------------------

*     CHANGED: use aji- instead of gf-values

*    ajimax:  maximum gf/gi/lami^2 (propto Aji) to ANY lower level 
*     gfmax:  maximum gf to ANY lower level (not needed)

*    maxaji:  maximum Aji to ilow  (1 or 'm')
*     maxgf:  corresponding gf-value
      
*     maxaji2: 2nd strongest Aji-value to low=1 or 'm'
*      maxgf2: corresponding gf-value

*     ilow: lower level for strongest transition
*     ilow2: lower level for 2nd strongest transition 

*     maxajisub : max Aji value for those transitons with up = ' ',
*                 low ne 1 and m and both levels belonging to gs-system
*        ilowsub: lower level of such a transition
      
      include 'metacom'
      integer mgf(-14:3),ngf(-14:3)
      real gfmax(5000),maxgf2(5000),maxaji(5000),maxaji2(5000)
      integer ilow2(5000)
      data mgf/18*0/,ngf/18*0/
      data gfmax/5000*0./,maxgf2/5000*0./,ilow2/5000*0/,
     +     maxaji/5000*0./maxaji2/5000*0./

      k=0
      ien=0

      do i=1,nrec
         meta(ien+1)=' '               ! ground states 
         metall(kel(i),iel(i),1) = ' '
         do j=1,ilin(i)
            k=k+1
            npart = kind(k)/10000
            lup = kind(k) - npart*10000
            lup = lup/100
            if(forb(k).ne.'f') then
               meta(lup+ien) =  ' '    ! all states except metastable ones 
               metall(kel(i),iel(i),lup) = ' '
            endif
         enddo
         ien=ien+ielev(i)
      enddo

* identify unique lower level belonging to subordinate lines
* ----------------------------------------------------------

      k=0
      ien=0
      n=1

 45   continue                  ! record loop

      lold=0

      read(12,*)
      do jj=1,ielev(n)
         read(12,*) rlam(jj),nquant,jstat(jj)
      enddo

      do j=1,ilin(n)
         k=k+1
         npart = kind(k)/10000
         nrest = kind(k) - npart*10000
         lup = nrest/100
         low = nrest - lup*100
* in case of multiple downward transitions, look for strongest
* transition; keep it if lower level is ground or metastable state

* find out also whether metastable is populated via excited level
* required that strong transition also to ground state

* with present philosophy, population of metastable level via
* excited level FROM lower lying metastable level not possible
* example: CI (3rd,'m') level from (2nd,'m') over (6th,'s'). 

C test output for HeI and NIV

c         if(n.eq.3.or. n.eq.12) then
c         if(lup.ne.lold) print*
c         print*,kind(k),lup,low,meta(ien+low),gf(k),jstat(low),forb(k)
c         endif

c in the following, we use Einstein-coefficients (because they are the decisive quantities)

         ajinew=gf(k)/float(jstat(lup))/wave(k)**2
c required for calculating aji_max
         gfnew=gf(k)
         if(lup.eq.lold) then
            if(((low.eq.1.) .or. (meta(ien+low) .eq. 'm')).and.
c this is the new condition: don't consider forbidden lines
     +      forb(k).ne.'f') then
              if(ajinew.gt.maxaji(ien+lup)) then
                maxaji2(ien+lup)=maxaji(ien+lup)
                maxgf2(ien+lup)=maxgf(ien+lup)
                ilow2(ien+lup)=ilow(ien+lup)

                maxaji(ien+lup)=ajinew
                maxgf(ien+lup)=gfnew
                ilow(ien+lup)=low

              else if(ajinew.gt.maxaji2(ien+lup)) then
                  maxaji2(ien+lup)=ajinew
                  maxgf2(ien+lup)=gfnew
                  ilow2(ien+lup)=low
              endif
            endif

         else 
            if(((low.eq.1.) .or. (meta(ien+low) .eq. 'm')).and.
     +      forb(k).ne.'f')then
               ilow(ien+lup)=low
               maxaji(ien+lup)=ajinew
               maxgf(ien+lup)=gfnew
            endif
            lold=lup
         endif
! maximum Aji and gf-value to ANY lower level
         ajimax(ien+lup)=max(ajimax(ien+lup),ajinew)
         gfmax(ien+lup)=max(gfmax(ien+lup),gfnew)

c         if(n.eq.3 .or. n.eq.12)
c     &  print*,ajimax(ien+lup),maxaji(ien+lup),maxaji2(ien+lup),
c     &        ilow(ien+lup),ilow2(ien+lup)
      enddo

* konlev = 1: considered level for line transitions (i.e. ground and
* metastable levels and upper levels connected to these by strong(est)
* (within factor 0.1) downward transition 

      do low=1,ielev(n) 
      if((low.eq.1.) .or. (meta(ien+low). eq. 'm'))then
            konlev(kel(n),iel(n),low)=1
      endif
      
      if(ilow(ien+low).ne.0) then
        if (maxaji(ien+low)/ajimax(ien+low).gt.0.1) then
c     &   .or.maxgf(ien+low).gt.0.5) then
              konlev(kel(n),iel(n),low)=1
              metall(kel(n),iel(n),low)= 's'
              meta(ien+low)='s'
        else
          ilow(ien+low)=0      
        endif
c     for further processing, we convert to gf-values
      endif      
      enddo

c loop to identify transitions belonging to ground-state multiplet
      k=k-ilin(n)
      do j=1,ilin(n)
         k=k+1
         if(forb(k).eq.'f') goto 30
         npart = kind(k)/10000
         nrest = kind(k) - npart*10000
         lup = nrest/100
         low = nrest - lup*100

         if(meta(ien+lup).eq.'s'.and.ilow(ien+lup).eq.1)
     +      gs_multi(ien+lup)=1

         if(meta(ien+lup).eq.' ') then
c         print*,low,lup,meta(ien+lup),gs_multi(ien+low)
           if (low.eq.1.or.gs_multi(ien+low).eq.1) then
             gs_multi(ien+lup)=1
           endif
         endif
30    continue
      enddo

c for tests
c      print*,kel(n),iel(n)
c      do jj=1,ielev(n)
c        print*,jj,meta(ien+jj),gs_multi(ien+jj),gfmax(ien+jj),
c     &    ajimax(ien+jj)
c      enddo      


c loop to identify maxajisub and ilowsub (belonging to gs-system)
      k=k-ilin(n)
      do j=1,ilin(n)
         k=k+1
         npart = kind(k)/10000
         nrest = kind(k) - npart*10000
         lup = nrest/100
         if(gs_multi(ien+lup).ne.1) goto 40
         low = nrest - lup*100

         ajinew=gf(k)/float(jstat(lup))/wave(k)**2
         if((low.ne.1.) .and. (meta(ien+low) .ne. 'm').and.
     +    gs_multi(ien+low).eq.1) then
           maxajisub(ien+lup)=max(maxajisub(ien+lup),ajinew)
           if(maxajisub(ien+lup).eq.ajinew) ilowsub(ien+lup)=low
         endif
40     continue
      enddo

      
      ien=ien+ielev(n)
      n=n+1
      if(n.le.nrec)goto 45

* gf-statistic for transitions

      rewind(12)

      open(30,file='gfdens.dat',status='unknown')
      ien=0
      do i=1,nrec
c         write(30,*)kel(i),iel(i)
         do l=2,ielev(i)
c            write(30,35)l,ilow(ien+l),gfmax(ien+l),meta(ien+l),
c     &           konlev(kel(i),iel(i),l)
c
c  identify populating level for those meta-stable levels fed by higher levels
           if(meta(ien+l).eq.'s' .and. ilow(ien+l) .ne.1 .and.
     &      ilow2(ien+l) .eq. 1 .and .maxgf2(ien+l).ge.0.1) then
             ip=ien+ilow(ien+l)
             if(ipop(ip).eq.0) then  ! take the lowest one
               ipop(ip)=l
c               print*,'case1',i,l,ilow(ien+l),ilow2(ien+l),ip,ipop(ip)
             endif  
           endif

           if(meta(ien+l).eq.'s' .and. ilow(ien+l) .eq.1 .and.
     &      ilow2(ien+l) .ne. 1 .and .maxgf2(ien+l).ge.0.1) then
             ip=ien+ilow2(ien+l)
             if(ipop(ip).eq.0) then  ! take the lowest one
               ipop(ip)=l
c               print*,'case2',i,l,ilow(ien+l),ilow2(ien+l),ip,ipop(ip)
             endif  
           endif
c
c  for further line-statistics, the following has to be modified
c
           if(gfmax(ien+l).ne.0.)then
               indgf = log10(gfmax(ien+l))
               if(gfmax(ien+l).lt.0.)indgf=indgf-1
               if(indgf.ge.-14 .and. indgf.lt.3)then
                  if(meta(ien+l).eq.'m')then
                     mgf(indgf)=mgf(indgf)+1
                  endif
                  ngf(indgf)=ngf(indgf)+1
               endif
           endif
         enddo
         ien=ien+ielev(i)
      enddo
 35   format(2x,2(i2,2x),e10.3,2x,a1,2x,i2)

      write(30,50)'#','dN/dgf','ngf','mgf'
      do i = -14,2
         del = 10.**(i+1) - 10.**i
         y = ngf(i)/del
         if(ngf(i).ne.0)then
            x = log10(float(ngf(i)))
         else
            x = 0.
         endif
         write(30,55)i,y,ngf(i),mgf(i)
      enddo

 50   format(a4,3(1x,a10))
 55   format(i4,1x,e10.4,2(1x,i10))
 19   format(i4,2(2x,i5))
      close(30)
  
      end

      subroutine outnl2
*----------------------------------------------------------------------
*     output of new file atomnl2_meta_v2.0: marks metastable levels and
*     has photoionisation cross sections added
*----------------------------------------------------------------------
      include 'metacom'
      real trans(5000)
      data trans/5000*0./
      
      open(13,file='atomnl2_meta_v2.0',status='unknown')

      k=0
      kold=0

      do i=1,nrec
         n1=kel(i)
         n2=iel(i)
         read(12,*)
         erg=ergion(n1,n2)

         do l=1,ielev(i)
            ibfup(l)=0
            beta(l)=0.
            sval(l)=0.
            cross(l)=0.
            gcross(l)=0.
         enddo

* read exact cross sections where available; add up channels to ex. levels
         if(levphi(i).ne.0)then
            read(16,*)i1,i2
            if(i2.ne.levphi(i))stop 'error in atompho'
         endif
         
         do l=1,levphi(i)
           read(16,*)ile,ib,alfa,b,s
           ener=erg-energy(n1,n2,ile)+energy(n1,n2+1,ib)
c renormalization
           if(b .gt. s+1.) then
              xmax=-b*s/(1.-b)/(1.+s)
              if(xmax.lt.0.or.xmax.gt.1.) then
                print*,b,s,xmax
                stop' something wrong in renormalization'
              endif  
              alfa=alfa*(b*xmax**s+(1.-b)*xmax**(s+1.))
              b=s+1.
           endif

c here comes the new hack: if more upper levels available, use only
c those transitions which have energies within 20% of the lowest one
c for those transitions, sum up the cross-sections

c note that in certain cases the population of the higher ion might
c become erroneous (if closely neightboured spin systems (e.g. OV),
c but at least the depopulation is OK then.
           
           if(ibfup(ile) .eq. 0) then
c first instance of ile
             enerfirst=ener
             alfamax=alfa
             ibfup(ile)=ib
             beta(ile)=b
             sval(ile)=s
             cross(ile)=alfa
           else  
             if (ener.lt.enerfirst) then
               print*,i,n1,n2,ile,ib,enerfirst,ener
               print*,' something wrong with bf to excited states'
               if(n1.eq.9 .and. n2.eq.2) then
                 print*,' (Fluor), needs to be updated, but OK thus far'
               else
                 stop ' should happen only for F_II'
               endif  
c happens for fluor
             else if (ener/enerfirst .lt. 1.2) then
               if (alfa.gt.alfamax) then
c in case, replace transition seaton parameters
                 alfamax=alfa
                 beta(ile)=b
                 sval(ile)=s
               endif
c used alpha value corresponds to the sum of contributing cross-sections
               cross(ile)=cross(ile)+alfa
             else
c otherwise, keep old transition
                continue
             endif         
           endif
         enddo
   
         do j=1,ielev(i)
            k=kold+j
            read(12,*) rlam(j),nquant,jstat(j)
            ground(k) = rlam(j)

* calculate cross section according to Gould (1978, ApJ 219, 250)
            ergi=erg-rlam(j)
            ergh=ergion(1,1)
            
            if(atom(k)(1:3).ne.'   ')then
               read(atom(k)(4:5),'(i2)') nel
            else
               nel=1
            endif

            xne=float(nel)/float(nquant)

* classical value (Kramers 1923, Phil. Mag. 46, 836)     
* this is an upper limit, since quant.-mech. a specific final state
* (ground state term L,S) is selected from a set of possible states 
* L'S' -> LS + (nl), where (L'S') ground state term of ionized atom
            gcross(j) = 7.907*ergh/ergi*xne

* correction is gweight = M(LS)/Sum(M(LS), where M(LS) is 
* (2S+1)(2L+1), sum is over possible final states. But we add over
* all target states:
            gweight = 1.
            gcross(j) = gweight*gcross(j)

            if(cross(j).eq.0.)then
               source(k)='G'
               cross(j)=gcross(j)
               sval(j)=2.
               beta(j)=2.
            endif
         enddo

* calculate zeta and eta values for each species
         summet=0.
         sumres=0.

         do j=1,ielev(i)
            k=kold+j
            part=float(jstat(j))/float(jstat(1))*cross(j)/cross(1)*
     +           ((erg-rlam(j))/erg)**2
            ocross(k) = part
            if(meta(k).eq.'m')then
               summet=summet+part
            else
               sumres=sumres+part
            endif
         enddo
         
         zeta=1./(summet+sumres)
         eta = summet*zeta

         write(13,14)i,kel(i),iel(i),ielev(i),
     +               ergion(kel(i),iel(i)),zeta,eta

         do j=1,ielev(i)
            k=kold+j
            if(meta(k) .ne. metall(kel(i),iel(i),j))then
               write(*,*)kel(i),iel(i),j
               stop
            endif
                         
            il=ilow(k)
            if(meta(k).eq.'s'.and.il.eq.1) then
              trans(k)=rlam(j)-rlam(1)
            endif  
            if(meta(k).eq.'s'.and.il.gt.1) then
              trans(k)=rlam(j)-rlam(il)
            endif  
            if(meta(k).eq.'m'.and.ipop(k).ne.0) then
              trans(k)=rlam(ipop(k))-rlam(1)
            endif  
            
            if(meta(k).eq.'s') then
              if(ipop(il+kold).ne.0) ipop(k)=il
            endif

c   identify ground-state system and provide ajimax
c   aji=6.6702082e15*ajimax

c   for those ' '-levels with known lower levels 's' or 'm'(ipop ne 0)
c   ilow is set to the actual value,
c   whereas for unknown lower levels (only few)
c   we provide only the identification of belonging to the gs-system
            
            if(meta(k).eq.' '.and.gs_multi(k).eq.1) then
              maxgf(k)=maxajisub(k)
              ils=ilowsub(k)
              ilsk=kold+ils
              if(gs_multi(ilsk).eq.1.and.
     &             meta(ilsk).eq.'s'.or.meta(ilsk).eq.'m') then
c check consistency
                if(meta(ilsk).eq.'m') then
                  print*,kel(i),iel(i),j,ils,maxgf(k),meta(ilsk)
                  stop' considered level should be s-level!' 
                endif
c                print*,'known',kel(i),iel(i),j,ils,maxgf(k)
                ilow(k)=ils
                trans(k)=rlam(j)-rlam(ils)
              else  
c                print*,'unkno',kel(i),iel(i),j,ils,maxgf(k)
                ilow(k)=1
              endif
            endif  

            if(ibfup(j) .eq. 0) then
              if (atom(k) .eq. ' ') then
                ibfup(j)=1
              else
                read (atom(k)(10:11),'(i2)') ibfup(j)
              endif
            endif  
      write(13,15)j,ibfup(j),ilow(k),rlam(j),jstat(j),meta(k),ipop(k),
     &        trans(k),maxgf(k),cross(j),ocross(k)/ocross(kold+1),
     &        beta(j),sval(j),source(k),gcross(j),atom(k)
         enddo

         kold=kold+ielev(i)
      enddo                     ! record loop

 14   format(i3,1x,i2,i4,2x,i4,11x,f9.1,2x,f4.3,3x,f4.3)
 15   format(i6,2x,i2,2x,i4,2x,f10.1,2x,i4,2x,'''',a1,'''',1x,i3,
     &  2x,f10.2,1x,
     &  e10.3,1x,f8.3,2x,f10.5,1x,f6.2,1x,f5.2,1x,a1,1x,f8.3,2x,a)

      return
      end

      subroutine forbout
*----------------------------------------------------------------------
*     output of all transitions, forbidden ones are flagged
*----------------------------------------------------------------------
      include 'metacom'

      open(14,file='gftest',status='unknown')

      k=0
      ien=0
      n=1
 55   continue
         lold=0
         gfold=0.

         if(ilin(n) .gt. 0)then
            write(14,'(2(2x,i4))')n,ilin(n)
            do j=1,ilin(n)
               k=k+1
               npart = kind(k)/10000
               nrest = kind(k) - npart*10000
               lup = nrest/100
               low = nrest - lup*100
c#a
               if(forb(k).eq.'f')then
                  write(14,18)kind(k),wave(k),gf(k),forb(k),
     &                 meta(ien+lup),atom(ien+low),atom(ien+lup)
               endif
c#e               
            enddo
         endif
         ien=ien+ielev(n)
         n=n+1
      if(n.le.nrec)goto 55
 18   format(2x,i8,2x,f12.3,2x,e8.3,4(2x,a))
     
      end

      subroutine forbid(atom,nrec,ilev,irec,ien,low,lup,ndec)
*----------------------------------------------------------------------
*     test whether suspicious transition satisfies electric dipole
*     transition selection rule (LS-coupling)      
*
*     INPUT :  
*        atom: level designations
*        irec: record number      
*         ien: cumulative level number up to present record
*         low: lower energy level
*         lup: upper energy level 
* 
*     OUTPUT:
*        ndec: 0 allowed, 1 forbidden      
*----------------------------------------------------------------------
      character atom(5000)*11,atlow*11,atlup*11,atl*1,atu*1
      dimension ilev(nrec)
     
* perform test if data available      
* -------------------------------
      ndec=0
      
      atlow=atom(ien+low)
      atlup=atom(ien+lup)

      if(atlow(1:3).eq.'   ' .or. atlup(1:3).eq.'   ')then
         return
      endif
      
* first test for parity inequality (Laporte's rule)
* -------------------------------------------------
* should not occur; occurrence is spurious due to false level
* designations by Margie Lennon in Kurucz data     

      if(atlow(9:9) .eq. atlup(9:9))then
         ndec = 1
         return
      endif
      
* then test for LS-coupling selection rules
* -----------------------------------------
      if(atlow(7:7) .ne. atlup(7:7))then    ! no intercombination
         ndec = 1
         return
      endif

      atl = atlow(8:8)
      if(atl .eq. 'S')then
         nlow = 0
      else if(atl .eq. 'P')then
         nlow = 1
      else if(atl .eq. 'D')then
         nlow = 2
      else if(atl .eq. 'F')then
         nlow = 3
      else if(atl .eq. 'G')then
         nlow = 4
      else if(atl .eq. 'H')then
         nlow = 5
      else if(atl .eq. 'I')then
         nlow = 6
      endif

      atu = atlup(8:8)
      if(atu .eq. 'S')then
         nlup = 0
      else if(atu .eq. 'P')then
         nlup = 1
      else if(atu .eq. 'D')then
         nlup = 2
      else if(atu .eq. 'F')then
         nlup = 3
      else if(atu .eq. 'G')then
         nlup = 4
      else if(atu .eq. 'H')then
         nlup = 5
      else if(atu .eq. 'I')then
         nlup = 6
      endif  

      if((abs(nlow-nlup) .gt. 1).or.(nlow.eq.0 .and. nlup.eq.0))then
         ndec = 1
         return
      endif
      
      end
      
      subroutine lines
*----------------------------------------------------------------------
*     preprocessing line list from nl3-files nl3*, *=i,a,g. For the
*     master line list, (quasi-)resonance lines and '1st order'
*     subordinate lines are selected (see Abbott & Lucy 1985,
*     ApJ 288, 679)
*
*     OUTPUT: 
*          file 'testlist' (ascii)
*
*     CHANGES:
*        11/96 convert Adi's real*8 format to real*4
*----------------------------------------------------------------------
      include 'metacom'
*     changed, since wmbasic files used
      integer*4 lineid(50000)
      real*4 ang(50000),gfval(50000)
      integer*4 lidneu(50000)
      real*4 angneu(50000),gfneu(50000)

      open(20,file=fpath2//'nl3i',form='unformatted',status='old',
     +     recl=4*50000,err=100,iostat=ios)
      open(21,file=fpath2//'nl3a',form='unformatted',status='old',
     +     recl=4*50000,err=101,iostat=ios)
      open(22,file=fpath2//'nl3g',form='unformatted',status='old',
     +     recl=4*50000,err=102,iostat=ios)
* take care to convert from dos to unix 
      open(23,file=fpath2//'nl3info',status='old',err=103,iostat=ios)

c      open(24,file = 'testlist', status = 'unknown')

* output files
      open(30,file = 'nl3i_all_ascii')
      open(31,file = 'nl3a_all_ascii')
      open(32,file = 'nl3g_all_ascii')

c      open(30,file = 'nl3i_all',form = 'unformatted', 
c     &     status = 'unknown')
c      open(31,file = 'nl3a_all',form = 'unformatted', 
c     &     status = 'unknown')
c      open(32,file = 'nl3g_all',form = 'unformatted', 
c     &     status = 'unknown')

      open(33,file = 'nl3info_all', status = 'unknown')
      read(23,*) ntotnl3,nfullrec,nlast
      print*,ntotnl3,nfullrec,nlast

* read recordwise and test selection criteria
      mlin = 0
      mtot = 0
      nmeta = 0
      nreso = 0
      nsub = 0
      nend = 50000
      neurec = 0
      neulast = 0

      do i=1,nfullrec+1
         if(i .eq. nfullrec+1)nend=nlast
         read(20)(lineid(j),j=1,nend)
         read(21)(ang(j),j=1,nend)
         read(22)(gfval(j),j=1,nend)
         
         do l=1,nend
            k=lineid(l)/1000000
            nrest=lineid(l)-1000000*k
            j=nrest/100000
            nrest=nrest-100000*j
            low=nrest/100
            lup=nrest-100*low

* test input data
            if((k.lt.1 .or. k.gt.natom) .or. (j.lt.1 .or. j.gt.8)
     &           .or. (low.lt.1))then
               print*,l,ang(l),lineid(l),k,j,low
            endif

c            if(gfval(l) .ge. 10.)then
c               print*,ang(l),lineid(l),gfval(l)
c            endif

* output lines according to some criteria
* ---------------------------------------
C-Jo        all wavelengths
C-Jo          if(ang(l).lt.1.e4.and.ang(l).ge.228.) then
c             if(ang(l).gt.1000..and.ang(l).lt.2000.) then
               mtot=mtot+1

C               write(24,10)ang(l),gfval(l),lineid(l)
C 10            format(f12.3,2x,e10.4,2x,i10)

* resonance + metastable as lower level
* use all lines!
*               if(konlev(k,j,low).eq.1)then 

                  if(low.eq.1)nreso = nreso + 1
                  if(metall(k,j,low) .eq. 'm')nmeta = nmeta + 1
                  if(metall(k,j,low) .eq. 's')nsub  = nsub  + 1
                  
                  mlin = mlin + 1
                  neulast = neulast + 1
                  lidneu(neulast) = lineid(l)
                  angneu(neulast) = ang(l)
                  gfneu(neulast)  = gfval(l)

* write reduced line lists
                  if(mod(mlin,50000) .eq. 0)then
                     neurec = neurec + 1
                     neulast = 0
                     write(30,*)(lidneu(n),n=1,50000)
                     write(31,*)(angneu(n),n=1,50000)
                     write(32,*)(gfneu(n),n=1,50000)
c                     write(30)(lidneu(n),n=1,50000)
c                     write(31)(angneu(n),n=1,50000)
c                     write(32)(gfneu(n),n=1,50000)
                     print*
     &,' record ',neurec,' written from ',angneu(1),' to ',
     &angneu(50000)
                  endif
                  
C-Jp           endif            ! selection for transition
C-Jo        endif               ! selection for wavelength and stage
         enddo                  ! block loop
      enddo                     ! record loop

      nsubo=mlin-(nreso+nmeta+nsub)
      write(*,*) 'mtot=',mtot,' mlin=',mlin
      write(*,*) 'nreso=',nreso,' nmeta=',nmeta,' nsub=',nsub
c      if (nsubo.ne. nsub) stop" Something's wrong with 's'-levels!"
      write(*,*) 'rest=',nsubo

      write(30,*)(lidneu(l),l=1,neulast)
      write(31,*)(angneu(l),l=1,neulast)
      write(32,*)(gfneu(l),l=1,neulast)
c      write(30)(lidneu(l),l=1,neulast)
c      write(31)(angneu(l),l=1,neulast)
c      write(32)(gfneu(l),l=1,neulast)
      print*
     &,' last record ',neurec+1,' written from ',angneu(1),' to ',
     &angneu(neulast)
c      
      write (33,*) mlin,neurec,neulast
      
      close(20)
      close(21)
      close(22)
      close(23)
c      close(24)

      close(30)
      close(31)
      close(32)
      close(33)

      return
 100  write(*,*) ' error in reading file nl3i, iostat=',ios
      stop
 101  write(*,*) ' error in reading file nl3a, iostat=',ios
      stop
 102  write(*,*) ' error in reading file nl3g, iostat=',ios
      stop
 103  write(*,*) ' error in reading file nl3info, iostat=',ios
      stop

      end

      SUBROUTINE indexx(n,arr,indx)
      INTEGER n,indx(n),M,NSTACK
      INTEGER arr(n)
      PARAMETER (M=7,NSTACK=50)
      INTEGER i,indxt,ir,itemp,j,jstack,k,l,istack(NSTACK)
      INTEGER a
      do 11 j=1,n
        indx(j)=j
11    continue
      jstack=0
      l=1
      ir=n
1     if(ir-l.lt.M)then
        do 13 j=l+1,ir
          indxt=indx(j)
          a=arr(indxt)
          do 12 i=j-1,1,-1
            if(arr(indx(i)).le.a)goto 2
            indx(i+1)=indx(i)
12        continue
          i=0
2         indx(i+1)=indxt
13      continue
        if(jstack.eq.0)return
        ir=istack(jstack)
        l=istack(jstack-1)
        jstack=jstack-2
      else
        k=(l+ir)/2
        itemp=indx(k)
        indx(k)=indx(l+1)
        indx(l+1)=itemp
        if(arr(indx(l+1)).gt.arr(indx(ir)))then
          itemp=indx(l+1)
          indx(l+1)=indx(ir)
          indx(ir)=itemp
        endif
        if(arr(indx(l)).gt.arr(indx(ir)))then
          itemp=indx(l)
          indx(l)=indx(ir)
          indx(ir)=itemp
        endif
        if(arr(indx(l+1)).gt.arr(indx(l)))then
          itemp=indx(l+1)
          indx(l+1)=indx(l)
          indx(l)=itemp
        endif
        i=l+1
        j=ir
        indxt=indx(l)
        a=arr(indxt)
3       continue
          i=i+1
        if(arr(indx(i)).lt.a)goto 3
4       continue
          j=j-1
        if(arr(indx(j)).gt.a)goto 4
        if(j.lt.i)goto 5
        itemp=indx(i)
        indx(i)=indx(j)
        indx(j)=itemp
        goto 3
5       indx(l)=indx(j)
        indx(j)=indxt
        jstack=jstack+2
        if(jstack.gt.NSTACK)pause 'NSTACK too small in indexx'
        if(ir-i+1.ge.j-l)then
          istack(jstack)=ir
          istack(jstack-1)=i
          ir=j-1
        else
          istack(jstack)=j-1
          istack(jstack-1)=l
          l=i
        endif
      endif
      goto 1
      END
C  (C) Copr. 1986-92 Numerical Recipes Software ,2:-5K#R..
