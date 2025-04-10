C
C  This file is part of EPOS4
C  Copyright (C) 2022 research institutions and authors (See CREDITS file)
C  This file is distributed under the terms of the GNU General Public License version 3 or later
C  (See COPYING file for the text of the licence)
C

c-----------------------------------------------------------------------
      subroutine emsaa(nfr,iret)
c-----------------------------------------------------------------------
c  energy-momentum sharing
c-----------------------------------------------------------------------

#include "aaa.h"
#include "ems.h"
#include "par.h"
#include "sem.h"
      integer ioTestFact
      common /ctestfact/ ioTestFact
      integer noptTestFact
      common /cnopttestfact/ noptTestFact
      common/cwzero/wzero,wzerox
      double precision omega,omlog,oma,omb,wab,wba,wmatrix,wzero,nbar
     *,wzerox,rrr,eps,xprem,xmrem,om1intgck,PhiExpok,drangen,om1intbc
     *,w(-2:9),xhrem,xxrem,om1intbyk
      parameter(eps=1.d-30)
      common/col3/ncol,kolpt,ncoli
      data jkoevtsum/0/ ikoevtsum/0/
      save jkoevtsum,ikoevtsum
c      logical modu
      double precision plc,s
      common/cems5/plc,s
      common/ems6/ivp0,iap0,idp0,isp0,ivt0,iat0,idt0,ist0
      common/emsptl/nppr(npommx,kollmx),npproj(mamx),nptarg(mamx)
      common/nucl3/phi,bimp
      common/ckopjtg/kopj(mamx),kotg(mamx)
      logical laccept
      dimension ishuff(2*mamx,2),icp(2),ict(2),jcp(nflav,2),jct(nflav,2)
     &          ,nishuff(2)
      dimension q2kk(2)
      double precision iutime(5),tiu3,tiu4,tidi
      logical ltime  
      data ncount/0/
      save ncount

      integer getAccumJerr

      call utpri('emsaa ',ish,ishini,4)

      call timer(iutime)
      tiu3=iutime(3)
      tiu4=iutime(4)
      ltime=.false.

      irea=iret

      do j=1,2
        do i=1,nflav
          jcp(i,j)=0
          jct(i,j)=0
        enddo
      enddo

      iret=0
      iret2=0
      iomegasave=iomega
      

c     initialize
c     ----------

      call GfunParK(0)   !nprmx(k) fixed in there so should be defined before emsigr
      call emsipt   !initialize projectile and target
      call emsigr   !initialize grid
      !call xBinary(2.)

      kd=0
      if(iomega.lt.2.and.koll.le.100)then
        k=0
        do while (k.lt.koll.and.kd.ge.0)
          k=k+1
          dprob=om1intbyk(k)
          if(drangen(dble(bk(1))).le.dprob)then
            if(kd.eq.0)then
              kd=k
            else            !only one diffractive pair allowed with Nuclei
              kd=-1
            endif
          endif
        enddo
      endif

      !~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
      if(koll.eq.1.and.kd.eq.1)then         !diffraction (for pp or very peripheral collisions) done first to save time
      !~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

      wzero=0d0                 ! to be sure to always produce a Pomeron here
      k=kd
      n=1
      call ProPo(k,n)
      call ProXY(k,n)
      if(ish.ge.4)write(ifch,*)'Diffractive process'

      !~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
      else                               !non-diffractive event
      !~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
        
      iomega=2
      call DefXminDf(s)
      call GfunParK(-1)
      
      call GfunParP !initialize pair position dependent parameters for omega calculation

      if(koll.eq.1.and.ikolmx.lt.1000000)then
       if(nprmx(1).lt.ikolmx)then
        write(ifmt,*)'lattice too small:',nprmx(1),ikolmx
        iret=1
        goto 1000
       endif
      endif

      !~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
      if(iokoll.eq.0)then !normal case
      !~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

c     initialize 2
c     ------------

      nSprmx=0
      do k=1,koll
        nSprmx=nSprmx+nprmx(k)
      enddo

      omlog=0
      nemsi=nemsi+1
      if(nemsi.le.4.and.iemsi1.eq.1)call xEmsI1(1,0,omlog)
      if(ish.ge.6)write (ifch,*)'after xEmsI1'
      if(nemsi.le.4.and.iemsi2.eq.1)call xEmsI2(1,0)
      if(ish.ge.6)write (ifch,*)'after xEmsI2'
      if(ish.ge.6)call XPrint('Before Markov:&')

      if(ish.ge.4)write(ifch,*)'Markov Process'
      kint=int(max(15.,10.*engy**0.2))
      if(koll.gt.50)kint=3*kint/int(log(float(koll)))
      kmcmx=nSprmx*kint        !50*kint  !100*kint

      laccept=.false.

c     start Metropolis
c     ----------------

      kmc=0
      do while (kmc.lt.kmcmx)
       kmc=kmc+1

       knprmx=0
       rrr=drangen(dble(kmc))
       do ik=1,koll
         knprmx=knprmx+nprmx(ik)
         if(rrr.le.dble(knprmx)/dble(nSprmx))then ! k-th pair
           k=ik
           goto 10
         endif
       enddo
 10    continue

       ip=iproj(k)
       it=itarg(k)
       n=1+int(rangen()*float(nprmx(k)))  ! n-th spot for k-th pair
       nbar=dble(npr(0,k))
       if(idpr(n,k).eq.0)nbar=nbar-1d0

c Alternative: generation follow om1*(1-xp)**alplea*(1-xm)**alplea
       xprem=1d0
       xmrem=1d0
c Alternative: generation follow om1*(xpr-xp)*(xmr-xm)
c       xprem=xpp(ip)
c       xmrem=xmt(it)
c       if(idpr(n,k).gt.0)then
c         xprem=xprem+xppr(n,k)
c         xmrem=xmrem+xmpr(n,k)
c       endif
c Alternative: generation follow om1
c       xprem=xpp(ip)
c       xmrem=xmt(it)
c       if(idpr(n,k).gt.0)then
c         xprem=xprem+xppr(n,k)
c         xmrem=xmrem+xmpr(n,k)
c       endif

       wzerox=(nbar+1d0)
       wzero=wzerox    / ( wzerox
     &                    +om1intc(k)*gammaV(k) )    !if xrem=1, faster to use om1intc
c     &                    +om1intgck(n,k,xprem,xmrem)*gammaV(k) )

       if(ish.ge.8)write(ifch,*)'wzero',k,n,wzero,wzerox,gammaV(k)
     &                      ,om1intc(k),om1intgck(n,k,xprem,xmrem)
       if(ish.ge.1.and.100000*(kmc/100000).eq.kmc)
     & write(ifmt,*)'kmc',kmc,kmcmx

       call StoCon(1,k,n)
       if(laccept.and.idpr(n,k).ne.0.and.nprt(1).le.ikolmn)then
        xpp(ip)=xpp(ip)+xppr(n,k)
        xmt(it)=xmt(it)+xmpr(n,k)
       else
        call RemPom(k,n)
        call ProPo(k,n)
       endif
       call ProXY(k,n)

       call StoCon(2,k,n)

       if(ish.ge.9)write(ifch,*)'propo',xppr(n,k),xmpr(n,k)
     &                                 ,idx0,idpr(n,k)
       if(idpr(n,k).eq.0.and.idx0.eq.0)then
         accept=accept+1.
       else
         omb=omega(n,k)
         if(omb.le.0.d0)then
           reject=reject+1.
           call RemPom(k,n)
           call StoCon(-1,k,n)
c           if(ish.ge.9)write(ifch,*)'omb',omb
         else

           wab=wmatrix(k,n)
           if(ish.ge.8)write(ifch,*)'wab',wab,omb
           if(wab.le.0.d0.or.omb.ne.omb)then
             write (ifch,*)'wab,kmc',wab,omb,kmc,k,n,xpr(n,k),ypr(n,k)
     &  ,xppr(n,k),xmpr(n,k),xpp(ip),xmt(it),ip,it,idpr(n,k)
             write(ifmt,'(a,i12,d25.15)')'ems,seedf',nrevt+1,seedc
             iret=1
             goto 1000
           endif
           call RemPom(k,n)
           call StoCon(-1,k,n)
           oma=omega(n,k)
           wba=wmatrix(k,n)
           if(oma.ge.0.d0.and.oma.le.eps*omb*wba/wab)then
             accept=accept+1.
             call RemPom(k,n)
             call StoCon(-2,k,n)
             omlog=omlog+dlog(omb)
             goto 500
           elseif(oma.le.1.d-300.or.oma.ne.oma.or.omb.ne.omb)then
             write (ifch,*)'oma,kmc',oma,omb,kmc,k,n,xpr(n,k),ypr(n,k)
     &  ,xppr(n,k),xmpr(n,k),idpr(n,k),npr(1,k),xpp(ip),xmt(it),ip,it
             write(ifmt,'(a,i12,d25.15)')'ems,seedf',nrevt+1,seedc
             iret=1
             goto 1000
           endif

           z=sngl(omb/oma*wba/wab)
           if(ish.ge.8)write(ifch,*)'z,oma',z,oma,wba,k,n
           if(rangen().gt.z)then
             reject=reject+1.
           else
             accept=accept+1.
             call RemPom(k,n)
             call StoCon(-2,k,n)
             omlog=omlog-dlog(oma)+dlog(omb)
           endif

 500       continue

         endif

       endif

       if(nemsi.le.4)then
         kplot=int(float(kmc)/float(kmcmx)*100.)
         if(iemsi1.eq.1)call xEmsI1(1,kplot,omlog)
         if(iemsi2.eq.1)call xEmsI2(1,kplot)
       endif

       if(kmc.eq.kmcmx.and.nprt(1).gt.0
     .  .and.(nprt(1).lt.ikolmn.or.nprt(1).gt.ikolmx)
     .  .and.nfr.gt.0)then
         kmc=0
         laccept=.true.
       endif

      enddo                     !-----> end Metropolis

c     second chance for diffraction for nuclei will 1<koll<100
      if(koll.gt.1.and.kd.gt.0)then

c Count interactions to check if collision is not non-diff

        ncol=0
        do k=1,koll
          if(nprt(k).gt.0)ncol=ncol+1
        enddo
        if(ncol.eq.0)then
          iomega=iomegasave
          call DefXminDf(s)
          call GfunParK(-1)
          wzero=0d0             ! to be sure to always produce a Pomeron here
          k=kd
          n=1
          call ProPo(k,n)
          call ProXY(k,n)
          if(ish.ge.4)write(ifch,*)'Diffractive nuclear process',kd
        endif
        
      endif

      
      !~~~~~~~~~~~~~~~~~~~~~~
      else !iokoll.ne.0
      !~~~~~~~~~~~~~~~~~~~~~~

        n=1

        rrr=drangen(dble(kmc))
        if(iokoll.lt.0)then
          wzero=om1intbc(bb)
          wzero=1d0-om1intbc(bk(1))/wzero
        else
          wzero=PhiExpoK(1,1d0,1d0)
        endif

        if(rrr.gt.wzero)then
          wzero=0d0          ! to be sure to always produce a Pomeron here
          do k=1,koll
           call ProPo(k,n)
           call ProXY(k,n)
          enddo
        endif

      !~~~~~~~~~~~~~~~~~~~~~~
      endif !iokoll
      !~~~~~~~~~~~~~~~~~~~~~~

      !~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
      endif  !diff / non-diff
      !~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

c --- check ---

      call timer(iutime)
      tidi=iutime(3)-tiu3 + (iutime(4)-tiu4)*1d-3
      if(ltime)write(ifmt,'(a,f5.2)')'cpu time emsaa 1/5',tidi
      tiu3=iutime(3)
      tiu4=iutime(4)

      if(ish.ge.3)then
        do k=1,koll
          kk=min(80,nprt(k))
          mm=0
          if(bk(k).lt.1.5)then
          mm=mm+1
          write(ifch,'(a,i4,f7.2,2i4,3x,90a1)')'POMS ',mm,bk(k)
     .     ,nprt(k),nprmx(k)
     .      ,('*',j=1,kk)
          endif
        enddo
      endif

c --- Initialize Pomeron x-distributions

      if(iemspx.eq.1)then
       q2kk(1)=0.
       q2kk(2)=0.
       call xxEmsPx(-1,0,0,0.,0.,0.,0.,0,0.,0.,q2kk)
      endif
      if(iemspbx.eq.1)then
       q2kk(1)=0.
       q2kk(2)=0.
       call xxEmsP2(1,-1,0, 0,  0.,0.,0.,0.,0.,0.,0.,0.,0.,0.,q2kk)
      endif
      if(iemspdf.eq.1)then
       q2kk(1)=0.
       q2kk(2)=0.
       call xxEmsP3(1,-1,0,0,0,0.,0.,0.,0.,q2kk)
      endif

c --- Plot Pomeron b-distributions ---

      if(ish.ge.6)call XPrint('After Markov :&')

      if(iemsb.eq.1)then ! plot
       do k=1,koll
        call xEmsB(1,1,k)
       enddo
      endif

      if(iemsbg.eq.1)then ! plot
        call xEmsBg(3,0,0,0)
        do k=1,koll
          call xEmsBg(1,0,k,0)
          if(nprt(k).gt.0)then
            call xEmsBg(1,-1,k,0)
          endif
        enddo
      endif


c --- Count all interactions ---

      ncol=0
      ncoltot=0
      do k=1,koll
        if(nprt(k).gt.0)then
c          if(iscreen.ne.0)call FuseEnhPom(k) !Pomeron fusion to mimic enhanced diagram behavior
          ncol=ncol+1
          ncoltot=ncoltot+nprt(k)
          ip=iproj(k)
          it=itarg(k)
          kolp(ip)=kolp(ip)+nprt(k) !number of cut Pomerons
          kolt(it)=kolt(it)+nprt(k) !on remnants
        endif
      enddo

c --- Plot distr of pomeron number ---


      if(iemspm.eq.1)then
       do k=1,koll
         call xEmsPm(1,k,nprt(k),nprmx(k))
       enddo
      endif

c --- Count

      ncoli=0     
      do k=1,koll
        q2kk(1)=q2nmin !q2kmin(1,0,k)
        q2kk(2)=q2nmin          !q2kmin(2,0,k)
        ip=iproj(k)
        it=itarg(k)
        idh=0
        xprem=xpp(ip)
        xmrem=xmt(it)
        xhrem=0d0
        xxrem=0d0
        do n=1,nprmx(k)
          if(idpr(n,k).ne.0)then
            q2kmin(1,n,k)=q2kk(1) 
            q2kmin(2,n,k)=q2kk(2) 
            idh=idh+1
            xxrem=max(xxrem,xpr(n,k))
            xhrem=xhrem+xpr(n,k)
            xprem=xprem+xppr(n,k)
            xmrem=xmrem+xmpr(n,k)
            if(xpr(n,k).gt.dble(esatur)/sqrt(dble(kolp(ip)*kolt(it)))
     .      )then
              npr(4,k)=npr(4,k)+1 !count only scattering significant for MPI
            endif
          endif
        enddo
        if(npr(4,k).gt.0)ncoli=ncoli+1
        xnpo(k)=1d0
        xnpo(k)=max(xnpo(k),dble(npr(4,k)))  
        xnpp(k)=dble(max(0,kolp(ip)-nprt(k))) !-dsatur*log(xprem)
        xnpt(k)=dble(max(0,kolt(it)-nprt(k))) !-dsatur*log(xmrem)
      enddo

c -- Fix Pomeron type (idpr=1 or idpr=3) ---

      do k=1,koll
      ip=iproj(k)
      it=itarg(k)
      if(nprt(k).gt.0)then
        do n=1,nprmx(k)
          if(idpr(n,k).ne.0.and.ivpr(n,k).eq.1)then
            call ProPoTy(k,n,0)
          endif
        enddo
      endif
      enddo

c --- Count, fix pair type and remnant excitations iep(),iet() ---

      ncol=0
      ncoli=0      !count hard scattering only here (temporary for ProPoTy and Q2s)
      do k=1,koll
        itpr(k)=0
        idh=0
        do n=1,nprmx(k)
          !================!
          call ProPoDif(k,n)
          !================! 
          if(idpr(n,k).eq.3.and.ivpr(n,k).eq.1)idh=idh+1
        enddo
        if(idh.gt.0)ncoli=ncoli+1
        if(itpr(k).lt.0)then       !inelastic
          ncol=ncol+1
c          ncoli=ncoli+1
          itpr(k)=-1
        elseif(itpr(k).gt.0)then    !diffractive
          ncol=ncol+1
          itpr(k)=-2
        endif
      enddo
      if(ish.ge.5)write(ifch,*)'ncol:',ncol
      !--------------------------------------------------------+
      !  itpr=  -2  only mass diagrams or exclusive CD         |
      !  itpr = -1  soft or hard exchange                      |
      !--------------------------------------------------------+

      call timer(iutime)
      tidi=iutime(3)-tiu3 + (iutime(4)-tiu4)*1d-3
      if(ltime)write(ifmt,'(a,f5.2)')'cpu time emsaa 2/5',tidi
      tiu3=iutime(3)
      tiu4=iutime(4)

c -- Fix Final value of Q2 ---

      ncolh=0    
      do k=1,koll
      ip=iproj(k)
      it=itarg(k)
      if(nprt(k).gt.0)then
        idh=0
        do n=1,nprmx(k)
          if(idpr(n,k).ne.0.and.ivpr(n,k).eq.1)then
            call ProPoTy(k,n,1) !calls WomTy
            if(idpr(n,k).eq.3.and.ivpr(n,k).eq.1)idh=idh+1  !count hard
          endif
        enddo
        if(idh.gt.0)ncolh=ncolh+1
      endif
      enddo

      call timer(iutime)
      tidi=iutime(3)-tiu3 + (iutime(4)-tiu4)*1d-3
      if(ltime)write(ifmt,'(a,f5.2)')'cpu time emsaa 3/5',tidi
      tiu3=iutime(3)
      tiu4=iutime(4)

      if(ioTestFact.eq.0)then
      ncoli=ncolh      !# of real hard scattering (temporary for WomTy and Q2s)
      do k=1,koll
      if(nprt(k).gt.0)then
        do n=1,nprmx(k)
          if(idpr(n,k).eq.3.and.ivpr(n,k).eq.1)then
            if(npr(3,k).eq.0)call utstop("Should never happen !!!&")
            call WomTy(w,n,k,2)
          endif
        enddo
      endif
      enddo
      endif

      call timer(iutime)
      tidi=iutime(3)-tiu3 + (iutime(4)-tiu4)*1d-3
      if(ltime)write(ifmt,'(a,f5.2)')'cpu time emsaa 4/5',tidi
      tiu3=iutime(3)
      tiu4=iutime(4)

c -- Write sum of Q2s into zzremn to be used for remnant excitation and mass

      call WriteZZ

c -- set npom dependent parameters  

      call setParams2() 

c -- Update Diffractive Pomeron type ---

      !---------------------------------------------------+
      !  mass diagr have idpr > 3 before call ProItfpr    |
      !  and then are suppressed                          |
      !---------------------------------------------------+
      do k=1,koll
        if(itpr(k).ne.0)call ProItfpr(k)       !define itpr
      enddo

      !----------------------------------------------------------------------+
      !  itpr=  -2  (soft diffraction) if only mass Poms exchanged           |
      !  itpr=  2   (soft central diffraction) if mass Poms are missing      |
      !  itpr=  1   (hard diffraction) if mass Poms are missing              |
      !                                (one or both side)                    |
      !  itpr = -1  (nondiffractive)      else                               |
      !----------------------------------------------------------------------+

      if(iemsbg.eq.1)then ! plot
        do k=1,koll
          if(itpr(k).ne.0)then
            do n=1,nprmx(k)
              if(idpr(n,k).ne.0.and.ivpr(n,k).eq.1.)then
                if(xpr(n,k).lt.1.d-3.and.idpr(n,k).eq.1)then  !low x soft pom
                  call xEmsBg(1,2,k,n)
                else
                  call xEmsBg(1,min(idpr(n,k),6),k,n)
                endif
                if(idpr(n,k).eq.3)then
                call xEmsBg(1,7,k,n) !plot q2kmn(0,k,1) for b of pairs
                call xEmsBg(1,8,k,n) !plot q2kmn(0,k,2) for b of pairs
                call xEmsBg(1,9,k,n) !plot q2kmn(0,k,1) for b of events
                call xEmsBg(1,10,k,n) !plot q2kmn(0,k,2) for b of events
                endif
              endif
            enddo
          endif
        enddo
      endif


c --- fix all variables


      if(ish.ge.4)write(ifch,*)'fix all variables'


c --- Update remnant excitation iep(),iet() and suppress mass diagram accordingly

      do ip=1,maproj
       if(lproj(ip).ne.0)then
         !=================!
         call ProReEx( 1,ip)
         !=================!
         if(iremn.ge.2)call UpdateFlav(ip,jcp,0) !reset jcpref to 0
       endif
      enddo
      do it=1,matarg
       if(ltarg(it).ne.0)then
         !=================!
         call ProReEx(-1,it)
         !=================!
         if(iremn.ge.2)call UpdateFlav(it,jct,0) !reset jctref to 0
       endif
      enddo


c --- Plot MC pomeron number ---

      if(nemsi.le.4.and.irea.ge.0)then
       if(iemsi1.eq.1)call xEmsI1(1,100,omlog)
       if(iemsi2.eq.1)call xEmsI2(1,100)
       if(iemsi1.eq.1.and.ncol.gt.0)call xEmsI1(2,0,omlog)
       if(iemsi2.eq.1.and.ncol.gt.0)call xEmsI2(2,0)
       if((iemsi1.eq.1.or.iemsi2.eq.1).and.ncol.eq.0)nemsi=nemsi-1
      endif

      if(iemsb.eq.1)then        ! plot
        do k=1,koll
          if(itpr(k).ne.0)call xEmsB(1,2,k) !something
          if(itpr(k).eq.0)call xEmsB(1,3,k) !nothing
          if(itpr(k).eq.-1)call xEmsB(1,4,k) !cut
          if(abs(itpr(k)).eq.2.or.itpr(k).eq.1)call xEmsB(1,5,k) !diffr
          if(itpr(k).eq.1)call xEmsB(1,6,k) !diffr cut
        enddo
      endif

c --- Overwrite coordpr 

      !-------------------------------------------------------
      ! Spatial parton distribution, depending on Pomeron configuration
      ! Bigger value of facpos gives smoother en density distribution. 
      !-------------------------------------------------------

      facpevt=0.
      facpcnt=0
      do k=1,koll
        call setfacpos(k)
        if(nprt(k).gt.0)then
          facpevt=facpevt+0.5*(facpos(1)+facpos(2))
          facpcnt=facpcnt+1.
        endif
        do n=1,nprmx(k)
          call conigr(n,k)
        enddo
      enddo
      if(facpcnt.gt.0)facpevt=facpevt/facpcnt

c --- Set string end type and pt, x for soft ones, check Pom mass

      iret0=0
      iret2=0
      do k=1,koll
 88     iret1=0
        do n=1,nprmx(k)
          iret=0
          if(idpr(n,k).gt.0.and.ivpr(n,k).eq.1)then
            call StringEnd(k,n,iret,ierr) !not for backup pomerons which will get all from their hard copy
            if(ierr.ne.0)then
              iret=1
              goto 1000 !redo event
            endif 
            if(iret.eq.0)iret0=iret0+1
            if(iret.eq.1)iret1=iret1+1
            if(iret.eq.2)iret2=iret2+1
          endif
        enddo
        if(iret1.gt.0)goto 88   !if some hard Pom deleted, define SE for backup
      enddo
      if(iret0.eq.0.and.iret2.gt.0)then
        if(ish.ge.2)then
          write(ifch,*)'All CD Pomeron lost, redo event !'
          write(ifmt,*)'All CD Pomeron lost, redo event !'
        endif
        iret=1
        goto 1000 !redo event
      endif

c --- Write /cevt/ ---

      call emszz
      if(ncol.eq.0)then
        iret=0
        goto 1000
      endif


      do k=1,koll
       if(itpr(k).ne.0)call emswrpom(k,iproj(k),maproj+itarg(k))
      enddo

      !call getMcentr !compute M centrality  
      !print*,'EMS M centr : ',mcentrf() 

c --- Determine event typ  typevt

      call EventType
      !--------------------------------------------------------------
      !   itpr   = 0    elastic | - typevt= 0   |                   |
      !   itpr   =-1      inel  | - typevt=-1   | both remn excit   |
      !   itpr   = 1       DPE  | - typevt=-2   | high mass diff    |
      !                   SDpro | - typevt=-3   | high mass diff    |
      !                   SDtar | - typevt=-4   | high mass diff    |
      !                         | - |typevt|>1  | NSD               |
      !   itpr   =-2      DD    | - typevt= 1   |     soft          |
      !                   SDpro | - typevt= 3   |     diff          |
      !                   SDtar | - typevt= 4   |     no pom        |
      !   itpr = 2        CD    | - typevt= 2   | exclusive central |
      !--------------------------------------------------------------

c --- Treat hard Pomerons
      if(ish.ge.8)call blist('list before psahot&',1,nptl)

      jkoevt=0
      do k=1,koll
        !if(itpr(k).ne.0)
        !.print*,'k itpr iep iet M M',k,itpr(k),iep(iproj(k)),iet(itarg(k)),
        !.sqrt(xpp(iproj(k))*xmp(iproj(k))*engy**2),
        !.sqrt(xpt(itarg(k))*xmt(itarg(k))*engy**2)
        estpom=0 !multiplicity estimate all poms per collision
        do n=1,nprmx(k)
         jsplitpom(n,k)=1
         if(idpr(n,k).ne.0.and.ivpr(n,k).eq.1)then
          estpom=estpom+log(xpr(n,k)*plc*plc) !deltaRap
         endif 
        enddo 
        do n=1,nprmx(k)
          if(idpr(n,k).eq.3.and.ivpr(n,k).eq.1)then
            !kkkkkkkkkkk!KW1909
            if(ioTestFact.ne.0)then
              q2kmin(1,n,k)=q2sft !q2nmin
              q2kmin(2,n,k)=q2sft !q2nmin
            endif
            !kkkkkkkkkkk
            q2cmin(1)=q2kmin(1,n,k)
            q2cmin(2)=q2kmin(2,n,k)
            call ipoBasics(q2cmin) !x dependent
            call ipoCSzeroTables(q2cmin)
            if(ishpom.eq.1)then
              if(ipsahot.eq.1)call testpsahot(k,n) !bg
              nptl0=nptl
              !====================
              call psahot(k,n,iret)
              !====================
              !kkkkkkkkkkk!KW1909
              if(iret.eq.55)goto 55
              if(ioTestFact.ne.0)goto 55
              !kkkkkkkkkkk
              if(ish.ge.8)write(ifch,*)'iret from psahot',iret
              if(iret.eq.1)then
                nptl=nptl0   !reset nptl to be sure no new particle is added
                ivi=16
                call VirPom(k,n,ivi)
                if(ivi.lt.0)then
                  call setAccumJerr(7,getAccumJerr(7)+1)
                  iret=1
                  goto 1000
                endif
                istptl(nppr(n,k))=32
              elseif(nbkpr(n,k).ne.0)then
                jkoevt=jkoevt+1 !successful hard Pom
                nn=nbkpr(n,k)
                ivi=17
                call VirPom(k,nn,ivi)
                if(ivi.lt.0)then
                  call setAccumJerr(7,getAccumJerr(7)+1)
                  iret=1
                  goto 1000
                endif
                istptl(nppr(nn,k))=32
                nbkpr(n,k)=0
              endif
              iret=0
            else
              istptl(nppr(n,k))=32
              if(nbkpr(n,k).ne.0)then
                nn=nbkpr(n,k)
                istptl(nppr(nn,k))=32
              endif
            endif
            if(ish.ge.8)call blist('list after psahot&',1,nptl)
          endif
        enddo
      enddo

c --- Set string end flavor for soft ---

      do k=1,koll
        do n=1,nprmx(k)
          if(ivpr(n,k).eq.1)then
            if(isopom.eq.1)then
              call ProSeF(k,n,iret)
              if(iret.eq.1)then
                ivi=18
                call VirPom(k,n,ivi)
                if(ivi.lt.0)then
                  call setAccumJerr(7,getAccumJerr(7)+1)
                  iret=1
                  goto 1000
                endif
                istptl(nppr(n,k))=32
              endif
              iret=0
            elseif(idpr(n,k).eq.1)then
              istptl(nppr(n,k))=32
            endif
          endif
        enddo
      enddo
      if(ish.ge.8)call blist('list after soft&',1,nptl)

c --- Diffractive Pt

      iret=1
      do k=1,koll
        call ProDiPt(k,1,iret)
      enddo
      if(iret.ne.0)then
        call setAccumJerr(8,getAccumJerr(8)+1)
        ivi=99
        if(ish.ge.2)then
          write(ifch,*)'All Pomeron lost, redo event !'
          write(ifmt,*)'All Pomeron lost, redo event !'
        endif
        iret=1
        goto 1000
      endif

      if(iremn.ge.2)then
c --- Add valence quark to jcpref and jctref for remnants ---
        do ip=1,maproj
          if(iep(ip).ne.-1)then
            call UpdateFlav(ip,jcp,10)
            do nnn=1,nrflav
              jcpval(nnn,1,ip)=jcp(nnn,1)
            enddo
            do nnn=1,nrflav
              jcpval(nnn,2,ip)=jcp(nnn,2)
            enddo
          else
            icp(1)=icproj(1,ip)
            icp(2)=icproj(2,ip)
            call iddeco(icp,jcp)
            do nnn=1,nrflav
              jcpval(nnn,1,ip)=jcp(nnn,1)
            enddo
            do nnn=1,nrflav
              jcpval(nnn,2,ip)=jcp(nnn,2)
            enddo
          endif
        enddo
        do it=1,matarg
          if(iet(it).ne.-1)then
            call UpdateFlav(it,jct,20)
            do nnn=1,nrflav
              jctval(nnn,1,it)=jct(nnn,1)
            enddo
            do nnn=1,nrflav
              jctval(nnn,2,it)=jct(nnn,2)
            enddo
          else
            ict(1)=ictarg(1,it)
            ict(2)=ictarg(2,it)
            call iddeco(ict,jct)
            do nnn=1,nrflav
              jctval(nnn,1,it)=jct(nnn,1)
            enddo
            do nnn=1,nrflav
              jctval(nnn,2,it)=jct(nnn,2)
            enddo
          endif
        enddo
      endif

      do ip=1,maproj
        !------------------------------------------------------------
        ! Here and later "kolp(ip).ne.0" replaced by "iep(ip).ne.-1"
        ! to count projectile and target nucleons which are counted in
        ! paires but are not used in collision (no diffractive or
        ! inelastic interaction) as slow particles at the end.
        ! Then we can use them in ProRem to give mass to all other
        ! nucleons and avoid energy conservation violation that utrescl
        ! can not treat (and it gives a reasonnable number of grey
        ! particles even if distributions are not really reproduced).
        !       if(kolp(ip).ne.0)call ProCop(ip,ip)
        !------------------------------------------------------------
        if(iep(ip).ne.-1)call ProCop(ip,ip)
      enddo
      do it=1,matarg
        if(iet(it).ne.-1)call ProCot(it,maproj+it)
        !if(kolt(it).ne.0)call ProCot(it,maproj+it)
      enddo


c --- Fix Pion Exchange in diffractive excited remnants

c      do ip=1,maproj
c       if(iep(ip)/10.eq.2)call ProReEx( 2,ip)
c      enddo
c      do it=1,matarg
c       if(iet(it)/10.eq.2)call ProReEx( -2,it)
c      enddo

c ---- Remnant Masses (ProReM)

      if(ish.ge.6)call XPrint('Before  ProReM:&')
      ntry=0
      iret=0
      iretshu=0
      call StoRe(1)             !Store Remnant configuration
 123  ntry=ntry+1
      nishuff(1)=0
      nishuff(2)=0
      do ip=1,maproj
        if(iep(ip).eq.0)then
          nishuff(1)=nishuff(1)+1
          ishuff(nishuff(1),1)=ip      !positive for non excited projectile
        elseif(iep(ip).gt.0)then
          nishuff(2)=nishuff(2)+1
          ishuff(nishuff(2),2)=ip      !positive for excited projectile
        endif
      enddo
      do it=1,matarg
        if(iet(it).eq.0)then
          nishuff(1)=nishuff(1)+1
          ishuff(nishuff(1),1)=-it !negative for non excited  target
        elseif(iet(it).gt.0)then
          nishuff(2)=nishuff(2)+1
          ishuff(nishuff(2),2)=-it !negative for excited  target
        endif
      enddo

      do while(nishuff(1)+nishuff(2).gt.0)

        if(nishuff(1).gt.0.and.nishuff(2).gt.0)then
          ir=1+int(rangen()+0.5)
        elseif(nishuff(1).gt.0)then
          ir=1
        else
          ir=2
        endif

        indx=1+int(rangen()*float(nishuff(ir)))
        if(ishuff(indx,ir).gt.0)then
          ip=ishuff(indx,ir)
          !=====================!
          call ProReM( 1,ip,iret)
          !=====================!
        else
          it=-ishuff(indx,ir)
          !=====================!
          call ProReM(-1,it,iret)
          !=====================!
        endif
        if(ish.ge.10)call XPrint('In  ProReM:&')

        if(iret.eq.1)then
          !----------------------------------------
          !If there is a problem, try again shuffle (30 times),
          !if it doesn't work, for pp, try 10 times with the same type
          !of event and if doesn't work redo event;
          !for pA redo event ; and for AB (with A or B >10)
          !continue with some ghosts ...
          !----------------------------------------
          if(ntry.lt.50)then
            if(ish.ge.3)write(ifch,*)'shuffle, try again',ntry
            call StoRe(-1)         !Restore Remnant configuration
            iret=0
            goto 123
          elseif(ntry.lt.100)then
            if(ish.ge.3)write(ifch,*)'shuffle, try again',ntry
            call StoRe(-1)         !Restore Remnant configuration
            iret=5
            goto 123
          elseif(maproj.le.20.or.matarg.le.20)then
            if(ish.ge.1)then
              write(ifch,*)'ProRem, redo event ! ntry=',ntry
              write(ifmt,*)'ProRem, redo event ! ntry=',ntry
            endif
            iret=1
            goto 1000
          else
            call StoRe(-1)         !Restore Remnant configuration
            iret=10
            iretshu=10
            goto 123
          endif
        endif

        ishuff(indx,ir)=ishuff(nishuff(ir),ir)
        nishuff(ir)=nishuff(ir)-1

      enddo
      !print*,'M M',
      !.sqrt(xpp(1)*xmp(1)*engy**2),sqrt(xpt(1)*xmt(1)*engy**2),
      !.0.5*log(xpp(1)/xmp(1)),0.5*log(xpt(1)*xmt(1)),vparam()/vparam()

c --- Correction for Diffractive Pt (from Ralph but seems to be less good for NA49)

c      do k=1,koll
c        call ProDiPt(k,2,idum)
c      enddo


      iret=0
      if(ish.ge.6)call XPrint('After ProReM:&')


c --- Write Remnants


      do ip=1,maproj
c       if(kolp(ip).ne.0)call emswrp(ip,ip)
       if(iep(ip).ne.-1)call emswrp(ip,ip)
      enddo

      do it=1,matarg
c       if(kolt(it).ne.0)call emswrt(it,maproj+it)
       if(iet(it).ne.-1)call emswrt(it,maproj+it)
      enddo

c --- Remnant Flavors (ProReF)


      do ip=1,maproj
        call ProReF(1,ip,iret)
        if(iret.ne.0)goto 1000
      enddo
      do it=1,matarg
        call ProReF(-1,it,iret)
        if(iret.ne.0)goto 1000
      enddo

c --- Correct energy

      if(abs(iokoll).le.0)then
      ihit=0
      do loo=1,2
      einit=maproj*engy/2+matarg*engy/2
      esu=0
      do i=1,nptl
        if(mod(istptl(i),10).eq.0)then
          esu=esu+pptl(4,i)
        endif
      enddo
      idf=abs(esu-einit)/einit*1000
      if(idf.gt.10.or.iretshu.eq.10.or.ihit.eq.1)then
        if(nfr.eq.0.and.loo.eq.2)then
          ncount=ncount+1 
          if(ncount.eq.1.or.ncount.eq.10.or.ncount.eq.100
     .      .or.ncount.eq.1000)then
            lcount=log10(1.*ncount)+1
            write(ifmt,'(a,i1,a,3x,3f10.1,2i6,i6)')
     .      'WARNING ',lcount,' ems Energy check'
     .      ,einit,esux,esu,idfx,idf,iretshu
          endif
          if(ish.ge.3)then
            write(ifch,'(a,5x,3f10.1,2i6,2i7)')
     .      'ems - energy check',einit,esux,esu,idfx,idf,iretshu,nrevt
            !call blist('list after ems energy check issue&',1,nptl)
            write(ifch,'(a)')
     .      '"list after ems energy check issue" not activated' 
          endif
        endif
        if(loo.eq.1.and.irescl.gt.0)then
          call utrescxx(iret,0,103) 
          !call utresc(iret,2.0)  !taking the default 1.02 will reject almost all
                                 !events for AuAu at 7 GeV
                                 !With 2.0 doable, but still many rejections
          !*****************************************************
          ! IMPORTANT when updating utresc:
          !*****************************************************
          ! In subroutine utresc: errlim changed to errlim=min(0.005,...)
          !   otherwise too big error accepted at low energies
          !*****************************************************
          ! VERY IMPORTANT: In case of many rejections, a bias is introduced,
          !  pretty much uncontrolled (probably high Npom is more affected etc)
          !    SHOULD BE AVOIDED (to be discussed!!!)
          !*****************************************************
          if(iret.gt.0)then
            iret=1
            goto 1000
          endif
        endif
        ihit=1
        esux=esu
        idfx=idf
      endif
      enddo
      endif

c     plot
c     ----

  55  continue

      call getMcentr   ! -> mcentrf() -> mmxcentrf()

      !call xEmsTest
      !call xEmsTest2
      !call xEmsTest3

      if(iemspx.eq.1)then
       do ko=1,koll
        if(itpr(ko).ne.0)then
         do np=1,nprmx(ko)
          if(idpr(np,ko).ge.1.and.idpr(np,ko).le.3
     .     .and.ivpr(np,ko).ne.2)then
            q2kk(1)=q2kmin(1,np,ko)
            q2kk(2)=q2kmin(2,np,ko)
            do jij=1,jsplitpom(np,ko)
              call xpprbor_get(jij,np,ko, xpprbor_ )
              call xmprbor_get(jij,np,ko, xmprbor_ )
              call ptprboo_get(1,jij,np,ko, ptprboo1_ )
              call ptprboo_get(2,jij,np,ko, ptprboo2_ )
              jot=1+max(-1,idhpr(np,ko)) !0=sat >0=nor
              call xxEmsPx(1,jot,jij,sngl(xpr(np,ko)),sngl(ypr(np,ko))
     *        ,xpprbor_,xmprbor_,idpr(np,ko),ptprboo1_,ptprboo2_,q2kk)
            enddo
          endif
         enddo
        endif
       enddo
      endif

      if(iemspbx.eq.1.or.iemspdf.eq.1)then
       do ko=1,koll
        if(nprt(ko).gt.0)then
         do np=1,nprmx(ko)
          if(idpr(np,ko).gt.0.and.ivpr(np,ko).ne.0)then
            q2kk(1)=q2kmin(1,np,ko)
            q2kk(2)=q2kmin(2,np,ko)
            je1=min(1,nemispr(1,np,ko))
            je2=min(1,nemispr(2,np,ko))
            jex=1+je1+2*je2
            if(idpr(np,ko).eq.3)then
             do jij=1,jsplitpom(np,ko)
              call xpprbor_get(jij,np,ko, xpprbor_ )
              call xmprbor_get(jij,np,ko, xmprbor_ )
              call ptprboo_get(1,jij,np,ko, ptprboo1_ )
              call ptprboo_get(2,jij,np,ko, ptprboo2_ )
              call rapprboo_get(1,jij,np,ko, rapprboo1_ )
              call rapprboo_get(2,jij,np,ko, rapprboo2_ )
              call gbpom_get(1,jij,np,ko, gbpom1_ )
              call gbpom_get(2,jij,np,ko, gbpom2_ )
              call idbor_get(1,jij,np,ko, idbor1__ )
              call idbor_get(2,jij,np,ko, idbor2__ )
              idbor1_=abs(idbor1__)
              idbor2_=abs(idbor2__)
              if(iemspbx.eq.1)call xxEmsP2(1,1+max(-1,idhpr(np,ko)),jex
     *            ,jij,sngl(xppr(np,ko)),sngl(xmpr(np,ko))
     *            ,xpprbor_,xmprbor_,ptprboo1_,ptprboo2_
     *            ,rapprboo1_,rapprboo2_,gbpom1_,gbpom2_,q2kk)
              if(iemspdf.eq.1)call xxEmsP3(1,1+max(-1,idhpr(np,ko)),jex
     *            ,idbor1_,idbor2_,xpprbor_,xmprbor_
     *            ,q2bor(np,ko),q2bor(np,ko),q2kk)
             enddo
            endif
          endif
         enddo
        endif
       enddo
      endif

      if(iret.eq.55)goto 1055

      if(iemsse.eq.1)then
       do ko=1,koll
        if(nprt(ko).gt.0)then
         do np=1,nprmx(ko)
          if(idpr(np,ko).gt.0.and.ivpr(np,ko).ne.0)then
           ptp1=sngl(xxp1pr(np,ko)**2+xyp1pr(np,ko)**2)
           ptp2=sngl(xxp2pr(np,ko)**2+xyp2pr(np,ko)**2)
           ptm1=sngl(xxm1pr(np,ko)**2+xym1pr(np,ko)**2)
           ptm2=sngl(xxm2pr(np,ko)**2+xym2pr(np,ko)**2)
           call xEmsSe(1,sngl(xp1pr(np,ko)),ptp1,1,1)
           call xEmsSe(1,sngl(xp2pr(np,ko)),ptp2,1,1)
           call xEmsSe(1,sngl(xm1pr(np,ko)),ptm1,-1,1)
           call xEmsSe(1,sngl(xm2pr(np,ko)),ptm2,-1,1)
           call xEmsSe(1,sngl(xp1pr(np,ko)),sngl(xm1pr(np,ko)),1,2)
           call xEmsSe(1,sngl(xm2pr(np,ko)),sngl(xp2pr(np,ko)),1,2)
          endif
         enddo
        endif
       enddo
      endif

      if(iemsdr.eq.1)then
       do i=maproj+matarg+1,nptl
        if(istptl(iorptl(i)).eq.41)then
          xpdr=(pptl(4,i)+pptl(3,i))/sngl(plc)
          xmdr=(pptl(4,i)-pptl(3,i))/sngl(plc)
          if(ityptl(i).eq.41)call xEmsDr(1,xpdr,xmdr,1)
          if(ityptl(i).eq.51)call xEmsDr(1,xpdr,xmdr,2)
          if(ityptl(i).eq.42)call xEmsDr(1,xpdr,xmdr,3)
          if(ityptl(i).eq.52)call xEmsDr(1,xpdr,xmdr,4)
        endif
       enddo
      endif

      if(iemsrx.eq.1)then
       if(maproj.eq.1.and.matarg.eq.1.and.idproj.eq.idtarg)then
         i=1
         j=1
         k=1
         if(iep(i)/10.ne.0.and.iet(j)/10.ne.0.and.nprt(k).gt.2)
     .   call xEmsRx(1,1,sngl(xpp(i)),sngl(xmp(i)))
         if(iep(i)/10.eq.0.and.iet(j)/10.ne.0.and.nprt(k).le.2)
     .   call xEmsRx(1,2,sngl(xmt(j)),sngl(xpt(j)))
       else
         do i=1,maproj
           if(kolp(i).gt.0)call xEmsRx(1,1,sngl(xpp(i)),sngl(xmp(i)))
         enddo
         do j=1,matarg
           if(kolt(j).gt.0)call xEmsRx(1,2,sngl(xmt(j)),sngl(xpt(j)))
         enddo
       endif
      endif

      if(ixbDens.eq.1)call xbDens(1)

c     exit
c     ----

 1000 continue
      iomega=iomegasave
      call DefXminDf(s)
c      write(*,*)'emsaa-iret',iret
      if(ish.ge.2.and.iret.ne.0)write(ifch,*)'iret not 0 (ems)=> redo'
     &                                       ,iret,ivi
 1055 continue
      call timer(iutime)
      tidi=iutime(3)-tiu3 + (iutime(4)-tiu4)*1d-3
      if(ltime)write(ifmt,'(a,f5.2)')'cpu time emsaa 5/5',tidi
      call utprix('emsaa ',ish,ishini,4)
      return
      end


c----------------------------------------------------------------------
      subroutine StoCon(mode,k,n)
c----------------------------------------------------------------------
c store or restore configuration
c   mode = 1 (store) or -1 (restore)
c   k = collision index
c   n = pomeron index
c----------------------------------------------------------------------

#include "aaa.h"
#include "ems.h"

      ip=iproj(k)
      it=itarg(k)

      if(mode.eq.1)then

       do i=0,3
        nprx0(i)=npr(i,k)
       enddo
       nprtx0=nprt(k)
       idx0=idpr(n,k)
       xxpr0=xpr(n,k)
       yx0=ypr(n,k)
       xxppr0=xppr(n,k)
       xxmpr0=xmpr(n,k)
c       nppx0=npp(ip)
c       nptx0=npt(it)
       xppx0=xpp(ip)
       xppstx0=xppmx(ip)
       xmpstx0=xppmn(ip)
       xmtx0=xmt(it)
       xptstx0=xmtmx(it)
       xmtstx0=xmtmn(it)

      elseif(mode.eq.2)then

       do i=0,3
        nprx(i)=npr(i,k)
       enddo
       nprtx=nprt(k)
       idx=idpr(n,k)
       xxpr=xpr(n,k)
       yx=ypr(n,k)
       xxppr=xppr(n,k)
       xxmpr=xmpr(n,k)
c       nppx=npp(ip)
c       nptx=npt(it)
       xppx=xpp(ip)
       xppstx=xppmx(ip)
       xmpstx=xppmn(ip)
       xmtx=xmt(it)
       xptstx=xmtmx(it)
       xmtstx=xmtmn(it)

      elseif(mode.eq.-1)then

       do i=0,3
        npr(i,k)=nprx0(i)
       enddo
       nprt(k)=nprtx0
       idpr(n,k)=idx0
       xpr(n,k)=xxpr0
       ypr(n,k)=yx0
       xppr(n,k)=xxppr0
       xmpr(n,k)=xxmpr0
c       npp(ip)=nppx0
c       npt(it)=nptx0
       xpp(ip)=xppx0
       xppmx(ip)=xppstx0
       xppmn(ip)=xmpstx0
       xmt(it)=xmtx0
       xmtmx(it)=xptstx0
       xmtmn(it)=xmtstx0

      elseif(mode.eq.-2)then

       do i=0,3
        npr(i,k)=nprx(i)
       enddo
       nprt(k)=nprtx
       idpr(n,k)=idx
       xpr(n,k)=xxpr
       ypr(n,k)=yx
       xppr(n,k)=xxppr
       xmpr(n,k)=xxmpr
c       npp(ip)=nppx
c       npt(it)=nptx
       xpp(ip)=xppx
       xppmx(ip)=xppstx
       xppmn(ip)=xmpstx
       xmt(it)=xmtx
       xmtmx(it)=xptstx
       xmtmn(it)=xmtstx

      else
      call utstop('mode should integer from -2 to 2 (without 0)&')
      endif
      return
      end

c-------------------------------------------------------------------------
      subroutine RemPom(k,n)
c-------------------------------------------------------------------------
c remove pomeron
c-------------------------------------------------------------------------
#include "aaa.h"
#include "ems.h"

      ip=iproj(k)
      it=itarg(k)
      npr(idpr(n,k),k)=npr(idpr(n,k),k)-1  !nr of pomerons
      nprt(k)=npr(1,k)+npr(2,k)+npr(3,k)
      if(idpr(n,k).gt.0)then
c       npp(ip)=npp(ip)-1                     !nr of pomerons per proj
c       npt(it)=npt(it)-1                     !nr of pomerons per targ
       idpr(n,k)=0
       xpp(ip)=xpp(ip)+xppr(n,k)
       xmt(it)=xmt(it)+xmpr(n,k)
       xpr(n,k)=0.d0
       ypr(n,k)=0.d0
       xppr(n,k)=0.d0
       xmpr(n,k)=0.d0



      endif

      end

c-------------------------------------------------------------------------
      subroutine ProPo(k,n)
c-------------------------------------------------------------------------
c propose pomeron type = idpr(n,k
c-------------------------------------------------------------------------
#include "aaa.h"
#include "ems.h"
      double precision wzero,wzerox,drangen
      common/cwzero/wzero,wzerox

      ip=iproj(k)
      it=itarg(k)

      idpr(n,k)=0

      if(drangen(wzero).gt.wzero)then
        idpr(n,k)=1


c nbr of pomerons per proj
c       npp(ip)=npp(ip)+1
c nbr of pomerons per targ
c       npt(it)=npt(it)+1

      endif

      npr(idpr(n,k),k)=npr(idpr(n,k),k)+1 !nr of pomerons
      nprt(k)=npr(1,k)+npr(2,k)+npr(3,k)


      end


c-------------------------------------------------------------------------
      subroutine ProXY(k,n)
c-------------------------------------------------------------------------
c propose pomeron x,y
c the Pomeron x distribution is dominated by this proposal
c if proposal for pomeron is consistent then x distribution is correct
c and event profile cross-section) but number of pomeron is always low
c and cannot exceed 5 even with different parameters. 
c If we use xrem=1 consistently then not enough pomeron neither. 
c If we use xrem=1 to propose pom only then too many pomerons are accepted
c If we propose pomerons with xrem=1 and propose x with xrem=sqt(xrem) then 
c everything perfect if npom not too large ???
c-------------------------------------------------------------------------

#include "aaa.h"
#include "par.h"
#include "ems.h"
#include "sem.h"
      double precision xp,xm,om1xprk,om1xmrk,axmt,axpp!,eps,drangen
     &,xprem,xmrem,om1x,om1xI,om1xo,om1xoI!,om1xrk,om1yrk
      external om1x,om1xI,om1xo,om1xoI
c      parameter (eps=1.d-30)


      ip=iproj(k)
      it=itarg(k)
      ntry=0
 1    ntry=ntry+1


      xpr(n,k)=0.d0
      ypr(n,k)=0.d0

      if(idpr(n,k).ne.0)then
          xprem=xpp(ip)
          xmrem=xmt(it)
c Generation follow om1*(xmr-xm)**alplea*(xpr-xp)**alplea
c because of fom, it's not symetric any more if we choose always xp first
c and then xm ... so choose it randomly.
c Alternative: generation follows om1*(xpr-xp)**alplea*(xmr-xm)**alplea between xmin and xr
        if(alplea(iclpro).ge.1..or.alplea(icltar).ge.1.)then
          if(rangen().lt.0.5)then
            xp=om1xprk(om1x,om1xI,n,k,xprem,xmrem,1)    !smaller x
            xm=om1xmrk(om1x,om1xI,n,k,xp,xprem,xmrem,1)    !smaller x
          else
            xm=om1xprk(om1x,om1xI,n,k,xmrem,xprem,-1)    !smaller x
            xp=om1xmrk(om1x,om1xI,n,k,xm,xmrem,xprem,-1)    !smaller x
          endif
          xpr(n,k)=xp*xm
          ypr(n,k)=0.d0
        else
c Alternative: generation follows om1 and x between xmin and xr
c and reject large x with probability (xr-x)**alplea so it follows om1*(xpr-xp)**alplea*(xmr-xm)**alplea
          if(rangen().lt.0.5)then
            xp=om1xprk(om1xo,om1xoI,n,k,xprem,xmrem,10)    !smaller x
            xm=om1xmrk(om1xo,om1xoI,n,k,xp,xprem,xmrem,10)    !smaller x
          else
            xm=om1xprk(om1xo,om1xoI,n,k,xmrem,xprem,-10)    !smaller x
            xp=om1xmrk(om1xo,om1xoI,n,k,xm,xmrem,xprem,-10)    !smaller x
          endif
          xpr(n,k)=xp*xm
          ypr(n,k)=0.d0
        endif
c          xpr(n,k)=om1xrk(n,k,xprem,xmrem)
c          ypr(n,k)=om1yrk(n,k,xpr(n,k),xprem,xmrem)
c          r=rangen()
c          if(r.lt.0.33)then
c            ypr(n,k)=-2d0
c          elseif(r.gt.0.66)then
c            ypr(n,k)=2d0
c          else
c            ypr(n,k)=0d0
c          endif
c          xp=sqrt(xpr(n,k))*exp(ypr(n,k))
c          xm=sqrt(xpr(n,k))*exp(-ypr(n,k))
          if(xpr(n,k).gt.XminDf)then
            ypr(n,k)=0.5D0*log(xp/xm)
            xppr(n,k)=xp
            xmpr(n,k)=xm
          elseif(ntry.lt.1000)then
            goto 1
          else
            if(ish.ge.10)write(ifch,*)
c            print *,
     .              'Warning in ProXY ',xp,xm,xpr(n,k),XminDf
            npr(idpr(n,k),k)=npr(idpr(n,k),k)-1
            idpr(n,k)=0
            npr(idpr(n,k),k)=npr(idpr(n,k),k)+1
            xpr(n,k)=0.d0
            ypr(n,k)=0.d0
            xppr(n,k)=0.d0
            xmpr(n,k)=0.d0
            nprt(k)=npr(1,k)+npr(2,k)+npr(3,k)
c            npp(ip)=npp(ip)-1   !nr of pomerons per proj
c            npt(it)=npt(it)-1   !nr of pomerons per targ
            return
          endif

c Update xp and xm of remnant, and change the limit to have big enought mass.

        xpp(ip)=xpp(ip)-xppr(n,k)
        xmt(it)=xmt(it)-xmpr(n,k)
        axmt=0d0
        do li=1,lproj(ip)
          kk=kproj(ip,li)
          itl=itarg(kk)
          if(xmt(itl).gt.0d0)axmt=axmt+xmt(itl)
        enddo
        xppmn(ip)=xpmn(ip)/min(2d0*axmt,xmpmx(ip)) !limit by momentum available on the other side
        axpp=0d0
        do li=1,ltarg(it)
          kk=ktarg(it,li)
          ipl=iproj(kk)
          if(xpp(ipl).gt.0d0)axpp=axpp+xpp(ipl)
        enddo
        xmtmn(it)=xtmn(it)/min(2d0*axpp,xptmx(it)) !limit by momentum available on the other side

      endif

      end

c-------------------------------------------------------------------------
      double precision function wmatrix(k,n)
c-------------------------------------------------------------------------
c proposal matrix w(a->b), considering pomeron type, x, y
c the number of Pomeron is dominated by this proposal
c so alternatice one give the same <npom> than Poisson but then x distribution
c are changed and cross-section a bit higher than from calculation
c other alternatice give correct x but lower npom
c-------------------------------------------------------------------------

#include "ems.h"
      double precision wzero,wzerox,Womegak,xprem,xmrem!,om1intgck
      common/cwzero/wzero,wzerox


      ip=iproj(k)
      it=itarg(k)

      if(idpr(n,k).eq.0)then
        wmatrix=wzero
      else
c Alternative: generation follow om1*(1-xp)**alplea*(1-xm)**alplea but fix xr=1. otherwise not enought MPI
          xprem=1d0
          xmrem=1d0
c Alternative: generation follow om1*(xpr-xp)*(xmr-xm)
c          xprem=xpp(ip)+xppr(n,k)
c          xmrem=xmt(it)+xmpr(n,k)
c Alternative: generation follow om1
c          xprem=xpp(ip)+xppr(n,k)
c          xmrem=xmt(it)+xmpr(n,k)
c IF SOMETHING CHANGED PLEASE CHECK IN OMG.F THAT WOMEGAK IS CHANGED TOO !

          wmatrix=(1d0-wzero)/om1intc(k)!om1intgck(n,k,xprem,xmrem) !if xrem=1., faster to use om1intc
     *           *Womegak(xppr(n,k),xmpr(n,k),xprem,xmrem,n,k)
      endif


      end

c-------------------------------------------------------------------------
      double precision function omega(n,k)
c-------------------------------------------------------------------------
c calculates partial omega for spot (k,n)
c-------------------------------------------------------------------------

#include "aaa.h"
#include "ems.h"
#include "sem.h"
      common/cwzero/wzero,wzerox
      double precision wzero,wzerox,eps
      parameter(eps=1.d-15)
      double precision PhiExpoK,omGamk,xp,xm
      double precision plc,s
      common/cems5/plc,s
      common/nucl3/phi,bimp
      common/crangrid/irangrid(npommx,kollmx)

      omega=0.d0

      ip=iproj(k)
      it=itarg(k)

c      write(ifch,*)'omega',xppmn(ip),xpp(ip),xppmx(ip)
c     &                    ,xmtmn(it),xmt(it),xmtmx(it)

      if(xpp(ip).lt.xppmn(ip)+eps.or.xpp(ip).gt.xppmx(ip)+eps)goto 1001
      if(xmt(it).lt.xmtmn(it)+eps.or.xmt(it).gt.xmtmx(it)+eps)goto 1001

      omega=xpp(ip)**dble(alplea(iclpro))
     &     *xmt(it)**dble(alplea(icltar))

c      zprj=zparpro(1,k)
c      ztgt=zpartar(1,k)
      if(idpr(n,k).eq.0)then
        omega=omega*wzerox
      else
        xp=xppr(n,k)
        xm=xmpr(n,k)
        omega=omega*omGamk(n,k,xp,xm,ntymin,ntymax)
     &             *gammaV(k)
      endif

      omega=omega*PhiExpoK(k,xpp(ip),xmt(it))

      if(omega.le.0.d0)goto 1001

      if(koll.gt.1)then
        do li=1,lproj(ip)
          kk=kproj(ip,li)
          if(itarg(kk).ne.it)then
            ipl=iproj(kk)
            itl=itarg(kk)
            omega=omega*PhiExpoK(kk,xpp(ipl),xmt(itl))
            if(omega.le.0.d0)goto 1001
          endif
        enddo
        do li=1,ltarg(it)
          kk=ktarg(it,li)
          if(iproj(kk).ne.ip)then
            ipl=iproj(kk)
            itl=itarg(kk)
            omega=omega*PhiExpoK(kk,xpp(ipl),xmt(itl))
            if(omega.le.0.d0)goto 1001
          endif
        enddo
      endif

      if(omega.lt.1.d-100)then
        if(ish.ge.6)write(*,*)'omega-exit',omega
        omega=0.d0
      elseif(omega.gt.1.d100)then
        if(ish.ge.6)write(*,*)'omega-exit',omega
        omega=0.d0
      endif

      return

 1001 continue

      omega=0.d0
      return

      end

c-------------------------------------------------------------------------
      subroutine ProPoTy(k,n,iqq)
c-------------------------------------------------------------------------
c propose pomeron type
c To keep Factorization formula, the exact value of om51p should be used
c and not the relative one. A Pomeron has been accepted using om1 (fit)
c defining the maximum amplitude, but to select Pomeron type we use om51p
c directly. The order for selection is important because if sum(om51p)>om1
c some processes may be reduced (and if sum(om51p)<om1 a process can be
c enhanced.
c - Diffractive component fixed here on top of Pomeron type selection
c since a hard Pomeron can be present in diff exchange (diff part incl in pdf)
c - Here we define q-q, q-g, g-q to get valence quarks as much
c right as possible, than we complete the hard component with g-g. The hard
c part is most important to get jets correctly.
c - Finally since the soft contribution contributes only to the multiplicity,
c it is used to complete ww (+ or -) to be consistent with cross-section.
c - idfpr (connection to remnant) is defined. It is used not to
c have valence quark interaction without remnant connection.
c iqq=0 diff/non-diff is defined to know how many Nhard there is
c iqq=1 gg/gq/qg/qq is defined with the proper value of Q2s taking into 
c account Nhard 
c-------------------------------------------------------------------------

#include "aaa.h"
#include "ems.h"
#include "sem.h"
#include "par.h"
      common/cems5/plc,s
      double precision s,plc,xvpr
      common/cems13/xvpr(0:3)
      double precision ww,w0,w1,w2,w3,w4,w5,w(-2:9),aks,drangen,om5Jk
      double precision eps
      parameter(eps=1.d-10)
      integer ioTestFact
      common /ctestfact/ ioTestFact
      integer idhprTestFact
      common /cidhprtestfact/ idhprTestFact

      if(idpr(n,k).ne.1.and.idpr(n,k).ne.3)return

      ip=iproj(k)
      it=itarg(k)

      if(ish.ge.4)write(ifch,*)'ProPoTy:k,n,ip,it,idpr,x,iqq',k,n,ip
     *                               ,it,idpr(n,k),xpr(n,k),iqq
      if(idpr(n,k).ne.1.and.iqq.eq.0)
     .call utstop('ProPoTy: should not happen&')

      if(iqq.eq.0)then
        idfpr(n,k)=1
        idpr(n,k)=1
      endif

      if(ioTestFact.ne.0)then
        npr(idpr(n,k),k)=npr(idpr(n,k),k)-1
        idpr(n,k)=3
        npr(3,k)=npr(3,k)+1
        return
      endif

      ww=0.d0
      w0=0.d0
      w1=0.d0
      w2=0.d0
      w3=0.d0
      w4=0.d0
      w5=0.d0
      do i=-2,9
        w(i)=0.d0
      enddo

      if(iqq.eq.0.and.iomega.lt.2)then  ! ----------------diffractif event-----------------

        do i=6,9
          w(i)=om5Jk(n,k,-10*i)
          ww=ww+w(i)
        enddo
        aks=drangen(ww)*ww
        if(aks.le.w(6))then     !Y+ (target excitation)
          idfpr(n,k)=3          !connexion on target
          wwnk(1,n,k)=w(6)
        elseif(aks.le.w(6)+w(7))then !Y- (projectile excitation)
          idfpr(n,k)=2          !connexion on projectile
          wwnk(1,n,k)=w(7)
        elseif(aks.le.w(6)+w(7)+w(8))then !.or.wdh.le.0d0)then !X (no excitation)
          idfpr(n,k)=4          !no connection to proj or targ
          wwnk(1,n,k)=w(8)
        else                    !DD
          idfpr(n,k)=0          !no Pomeron anymore 
          idpr(n,k)=7
        endif

        q2kmin(1,n,k)=q2nmin    !if diff set Q2s to minimum
        q2kmin(2,n,k)=q2nmin
c Warning : if idpr for diff is changed, please change behavior of Virpom for idpr>3
            
        if(ish.ge.4)write(ifch,*)'ProPoTy:diff ',idpr(n,k),idfpr(n,k)

      elseif(xpr(n,k).gt.xvpr(1))then   ! -------------require minimum mass otherwise to DD------------

        if(iqq.eq.0)then
          ww=om5Jk(n,k,-1)        !G for MC : diff. is included
        elseif(idpr(n,k).eq.3)then
          ww=wwnk(2,n,k)
        elseif(idpr(n,k).eq.1)then
          ww=wwnk(1,n,k)
        else
          call utstop("Should not happen in ProPoTy !&")
        endif
        aks=drangen(ww)*ww

        if(iregge.ne.0.and.iqq.eq.0)w(9)=om5Jk(n,k,5)

        if(aks.lt.w(9))then ! ---Reggeon---

          idpr(n,k)=2
          npr(2,k)=npr(2,k)+1
          npr(1,k)=npr(1,k)-1
          if(ish.ge.4)write(ifch,*)
     .     'ProPoTy:reggeon',idpr(n,k),idfpr(n,k)
       
        else ! ---no Reggeon---

          if(ishpom.eq.1)then
            call WomTy(w,n,k,iqq)
          else
            w(0)=1d0
            ww=1d0
          endif
          if(w(0).gt.0.d0)w0=w(0) !soft if iqq=1, SD+CD if iqq=0
          if(w(1).gt.0.d0)w1=w(1) !hard (incl. diff)
          if(w(2).gt.0.d0)w2=w(2)
          if(w(3).gt.0.d0)w3=w(3)
          if(w(4).gt.0.d0)w4=w(4)
          if(w(5).gt.0.d0)w5=w(5) !hard sea-sea no evolution
        
          wdh=w1+w2+w3+w4+w5
          if(iqq.eq.0)then
            aks=aks-w(9)
            if(iomega.lt.2)then   !fix SD and CD
              do i=6,8
                w(i)=om5Jk(n,k,-10*i) 
              enddo
              if(aks.le.w(6))then !Y+ (target excitation)
                idfpr(n,k)=3      !connexion on target
              elseif(aks.le.w(6)+w(7))then !Y- (projectile excitation)
                idfpr(n,k)=2      !connexion on projectile
              elseif(aks.le.w(6)+w(7)+w(8))then !.or.wdh.le.0d0)then !X (no excitation)
                idfpr(n,k)=4      !no connection to proj or targ
              endif
            endif
            ww=ww-w(9)
            aks=drangen(ww)*ww    !decide if it is soft or hard or DD (by deduction) for iqq=0
          endif
      
          if(ish.ge.6)write(ifch,*)'ProPoTy:ww,aks,ww_i(%),idfpr'
c         print*,'ProPoTy:ww,ww_i'
     *       ,ww,aks/ww*100.d0,w0/ww*100.d0,w1/ww*100.d0,w2/ww*100.d0
     *       ,w3/ww*100.d0,w4/ww*100.d0,w5/ww*100.d0
     *       ,w(6)/ww*100.d0,w(7)/ww*100.d0,w(8)/ww*100.d0
     *       ,wdh/ww,(w0+wdh)/ww,w(9)/ww,idfpr(n,k)

          
          if(ww.gt.eps.and.aks.le.wdh)then !hard pomeron (iqq=1) with soft for iqq=0

            npr(idpr(n,k),k)=npr(idpr(n,k),k)-1
            idpr(n,k)=3
            npr(3,k)=npr(3,k)+1
c            if(iqq.eq.0.and.w1.gt.w0)npr(4,k)=npr(4,k)+1  !count only scattering with high chance to be counted as hard (to fix Q2s)
c            if(iqq.eq.0.and.xpr(n,k).gt.dble(dsatur)/dble(nprt(k)))
c     .      npr(4,k)=npr(4,k)+1   !count only scattering significant for MPI
            if(ish.ge.4)write(ifch,*)'ProPoTy:idpr',idpr(n,k)
c            aks=wdh      !for tests
            if(iqq.gt.0)then
              if(aks.gt.w1+w2+w3+w4)then !soft+gg-pomeron
                wwnk(2,n,k)=w(-1) !w5
                idhpr(n,k)=-1
c               if(aks.gt.w1+w2+w3+w4+w(-1))idhpr(n,k)=-2 !bottom contribution
              elseif(aks.gt.w1+w2+w3)then !qq-pomeron
                wwnk(2,n,k)=w4 !cKW2108
                idhpr(n,k)=3
                ivp(ip)=ivp(ip)-1
                ivt(it)=ivt(it)-1
              elseif(aks.gt.w1)then !qg or gq-pomeron
                if(rangen().gt.0.5)then !randomize not to always pick-up the same if ww is too small
                  if(aks.gt.w1+w3)then !qg-pomeron
                    wwnk(2,n,k)=w2 !cKW2108
                    idhpr(n,k)=1
                    ivp(ip)=ivp(ip)-1
                  else            !gq-pomeron
                    wwnk(2,n,k)=w3 !cKW2108
                    idhpr(n,k)=2
                    ivt(it)=ivt(it)-1
                  endif
                else
                  if(aks.gt.w1+w2)then !gq-pomeron
                    wwnk(2,n,k)=w3 !cKW2108
                    idhpr(n,k)=2
                    ivt(it)=ivt(it)-1
                  else            !qg-pomeron
                    wwnk(2,n,k)=w2 !cKW2108
                    idhpr(n,k)=1
                    ivp(ip)=ivp(ip)-1
                  endif
                endif
              else                !gg-pomeron
                wwnk(2,n,k)=w1 !cKW2108
                idhpr(n,k)=0
              endif
              if(ish.ge.5)write(ifch,*)'ProPoTy:idhpr'
c             print*,'ProPoTy:idhpr'
     &           ,idhpr(n,k),idfpr(n,k),' |',ip,ivp(ip),' |',it,ivt(it)
            endif

          elseif(ww.gt.eps.and.(aks.le.wdh+w0.or.idfpr(n,k).ne.1
     &                            .or.iqq.ge.1.or.iomega.eq.2))then !soft pomeron + SD + CD
            q2kmin(1,n,k)=q2nmin 
            q2kmin(2,n,k)=q2nmin 
            npr(idpr(n,k),k)=npr(idpr(n,k),k)-1
            idpr(n,k)=1
            npr(1,k)=npr(1,k)+1
            if(iqq.eq.1)vparam(80)=vparam(80)+1
c           if(iqq.eq.0.and.xpr(n,k).gt.dble(dsatur)/dble(nprt(k)))
c     .     npr(4,k)=npr(4,k)+1   !count only scattering significant for MPI
            if(ish.ge.4)write(ifch,*)'ProPoTy:soft',idpr(n,k),idfpr(n,k)
            
          else                    !DD

            aks=aks-wdh-w0
            if(iokoll.eq.0.or.aks.gt.0d0)then !remaining contri. is low mass diffraction
              q2kmin(1,n,k)=q2nmin !if diff set Q2s to minimum
              q2kmin(2,n,k)=q2nmin
c Warning : if idpr for diff is changed, please change behavior of Virpom for idpr>3
              idfpr(n,k)=0        !no Pomeron anymore 
              idpr(n,k)=7
              if(ish.ge.4)write(ifch,*)'ProPoTy:DD',idpr(n,k),idfpr(n,k)
            else
              itpr(k)=0           !interaction should be completely deleted when iokoll.ne.0
              call VirPom(k,n,0)
              return 
            endif

          endif

        endif  ! ---no Reggeon---

      else ! ------------------ low mass --------------------

        idfpr(n,k)=0            !no Pomeron anymore 
        idpr(n,k)=7

        q2kmin(1,n,k)=q2nmin    !if diff set Q2s to minimum
        q2kmin(2,n,k)=q2nmin
c Warning : if idpr for diff is changed, please change behavior of Virpom for idpr>3
            
        if(ish.ge.4)write(ifch,*)'ProPoTy:DD low mass '
     .                           ,idpr(n,k),idfpr(n,k)

        
      endif  ! ------------------------------------------------

      if(idpr(n,k).gt.0)then
        if(iqq.eq.0)then
          antot=antot+1
          antotf=antotf+1
        else
          if(abs(idpr(n,k)).eq.1)then
            ansf=ansf+1
            ansff=ansff+1
          endif
          if(abs(idpr(n,k)).eq.3)then
            ansh=ansh+1
            anshf=anshf+1
          endif
        endif
      endif

c      if(abs(idpr(n,k)).eq.3.and.xpr(n,k).lt.xggfit)then
      if(abs(idpr(n,k)).eq.3.and.iqq.gt.0)then

                       !Backup soft Pomeron if sh not possible later

          kb=k
          nb=n
          ip=iproj(kb)
          it=itarg(kb)
          do nn=1,nprmx(kb)
            if(idpr(nn,kb).eq.0)then !empty spot
              nbkpr(nb,kb)=nn
              nvpr(nn,kb)=nb
              idpr(nn,kb)=1
              ivpr(nn,kb)=2
              xpr(nn,kb)=xpr(nb,kb)
              ypr(nn,kb)=ypr(nb,kb)
              xppr(nn,kb)=xppr(nb,kb)
              xmpr(nn,kb)=xmpr(nb,kb)
              idfpr(nn,kb)=idfpr(nb,kb)
c              bhpr(nn,kb)=bhpr(nb,kb)
              idp1pr(nn,kb)=0
              idp2pr(nn,kb)=0
              idm1pr(nn,kb)=0
              idm2pr(nn,kb)=0
              xm1pr(nn,kb)=0.d0
              xp1pr(nn,kb)=0.d0
              xm2pr(nn,kb)=0.d0
              xp2pr(nn,kb)=0.d0
              xxm1pr(nn,kb)=0.d0
              xym1pr(nn,kb)=0.d0
              xxp1pr(nn,kb)=0.d0
              xyp1pr(nn,kb)=0.d0
              xxm2pr(nn,kb)=0.d0
              xym2pr(nn,kb)=0.d0
              xxp2pr(nn,kb)=0.d0
              xyp2pr(nn,kb)=0.d0
              goto 10
            endif
          enddo
          if(ish.ge.2)write(ifmt,*)'no empty lattice site, backup lost'
c          if(ish.ge.2)write(ifch,*)'no empty lattice site, backup lost'

 10       continue
      endif

      return
      end

c-----------------------------------------------------------------------
      subroutine WomTy(w,n,k,iqq)
c-----------------------------------------------------------------------
c k - pair indice; n - position indice;
c Computes w(...)
c w(6) used as input for Y+ diff MC contribution
c w(7) used as input for Y- diff MC contribution
c w(8) used as input for X diff MC contribution
c w(9) used as input for Regge contribution
c    ity = 0   - soft
c    ity = 1   - gg
c    ity = 2   - qg
c    ity = 3   - gq
c    ity = 4   - qq
c    ity = 5   - sat
c iqq = 0,1,2 (iteration, womty is called several times)
c-----------------------------------------------------------------------
#include "aaa.h"
#include "par.h"
#include "ems.h"
#include "sem.h"
#include "tab.h"
      double precision plc,s
      common/nucl3/phi,bimp
      common/cems5/plc,s
      common/cbglaub/bglaub
      common/cnglndi/nglndik
      common/col3/ncol,kolpt,ncoli
      double precision om5Jk,w(-2:9),xp,xm,xh,yp,ww,omt,wsh,w0
     & ,cgcQ2s,eps,xvpr,fxxx,wdd
     & ,w1
     & ,PomIncXExactk,PomIncPExactk,PomIncMExactk,PomIncExactk,PomInc
     & ,PomIncExact
      real binfac,binfun
      common/cems13/xvpr(0:3)
      common/geom/rmproj,rmtarg,bmax,bkmx
      logical idf(-1:4),lndif(2)
      integer noptTestFact
      common /cnopttestfact/ noptTestFact
      logical lpass,lsat
      external PomIncXExactk,PomIncExactk
      dimension xnpoma(0:3),nrevtxx(0:1),nrevtxxx(0:2)
      data xnpoma/0.,0.,0.,0./,nrevtxx/1,0/,nrevtxxx/0,0,0/
      save xnpoma,nrevtxx,nrevtxxx
      data ncount/0/
      save ncount

      dummy=PomIncExact(1d0,1d0,0.)
      call  getSystemType(isys,amassMax,amassAsy)

      call setOptionsWom() 
      call deformoptget(1,iodeform)  !2 = new deformation function method (DFM2)
      if(iodeform.eq.1)stop'ERROR iodeform = 1 not supported any more' !removed in 3444m
      !iodeform was used to compute 'binfac' and 'binfun' and later for plots in xEmsPx, using old/new
      ip=iproj(k)
      it=itarg(k)
      xh=xpr(n,k)
      yp=ypr(n,k)
      spom=xh*plc**2
      nmx=n                    
      xp=xppr(nmx,k)
      xm=xmpr(nmx,k)
      rho=zparpnx(2,k)+zpartnx(2,k)
      cgcq2s=1d0
      binfun=1.0
      con1=xnpp(k)
      con2=xnpt(k)
      xpome=xnpo(k)             !float(max(1,npr(4,k)))
      conn1=con1+xpome
      conn2=con2+xpome
      zpom=0
      do kk=1,koll
        zpom=zpom+nprt(kk)
      enddo
      jval=igetJval(maproj,matarg,engy)
      binfac=1.0   
      
      if(iqq.eq.0)then  ! ================================= iqq=0 ================================

        q2pmin(1)=q2nmin        !start at minimum Q2s
        q2pmin(2)=q2pmin(1)
        call ipoOm5Tab(xp,xm,bk(k))
        wdd=0d0                 !if 0, soft is not rescaled, 1. to scale soft

        if(jval.eq.1)then
          zval=(conn1+conn2)/2
        elseif(jval.eq.2)then
          zval=log(zpom)
        else
          stop'ERROR 09062022a'
        endif
        binfun=1/deform(maproj,matarg,engy,zval,sngl(xp*xm)) 
        if(iomega.ne.2)then
          ww=max(0d0,om5Jk(n,k,-1)-w(9)) !Reggeon already taken out
          wwnk(0,n,k)=max(0d0,min(ww,om5Jk(n,k,-3)))
          wwnk(1,n,k)=max(0d0,om5Jk(n,k,101)+om5Jk(n,k,102)
     .         +om5Jk(n,k,103)+om5Jk(n,k,104))
          w(1)=max(0d0,wwnk(1,n,k)*binfac*binfun)
          w(1)=min(w(1),wwnk(0,n,k))
          w1=max(0.d0,ww-w(1))  !diff diagram with high x should only be soft, not used for scaling
          w(0)=max(0d0,w1-om5Jk(n,k,-90)) !don't count DD in possible soft Pom
        else !normal case iomega=2
          ww=max(0d0,om5Jk(n,k,-1))    !fit all
          wwnk(0,n,k)=max(0d0,min(ww,om5Jk(n,k,-3)))    !fit hard only
          wwnk(1,n,k)=max(0d0,om5Jk(n,k,101)+om5Jk(n,k,102)
     .         +om5Jk(n,k,103)+om5Jk(n,k,104))!true hard with special screening
          w(1)=max(0d0,wwnk(1,n,k)*binfac*binfun)
          w(1)=min(w(1),wwnk(0,n,k))
          w(0)=max(0d0,ww-w(1))
        endif
        if(isopom.eq.0)then     !no soft
          w(1)=w(1)+w(0)
          w(0)=0.d0
        endif
        if(ishpom.eq.0)then     !no hard
          w(0)=w(1)+w(0)
          w(1)=0.d0
        endif
        omt=w(1)
        wwnk(2,n,k)=w(1) !hard
        wwnk(1,n,k)=w(0) !soft
        return

      elseif(iqq.eq.1)then  ! ============================= iqq=1 ============================    !for iqq=1 hard/soft allready fixed

        !- - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - 
        !The total weight ww is determined in the call before (iqq=0) via wwnk, essentially:
        !  ww (hard) =  om5Jk(101-104)   =  sum of tabu-hard with hard screening from Zpair(-1)
        !  ww (soft) =  ALL - om5Jk(101-104)   
        !  ALL = om5Jk(-1) = fit-all = omGamk(0-3) 0=DD,1=ND,2-3=SD (different epsilons)
        !- - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - 
        !In the following ww is the hard or the soft one, depending on idpr
        !One conputes ww = wsh(q2s) = (1-wdd)*om5Jk(0) + binfun*om5Jk(1-4) + om5Jk(11)
        ! taking only hard or soft/sat contributions depending on idpr
        !- - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - 
        if(idpr(n,k).eq.1)then    !Soft or Sat
          q2pmin(1)=q2mnval(2)      !do not start at minimum Q2s but maximum Pomsat
          q2pmin(2)=q2pmin(1)
          ww=wwnk(1,n,k) 
          if(factsat.le.0.)then  !to avoid large Q2s due to small ww
            q2pmin(1)=q2mnval(1)
            q2pmin(2)=q2mnval(1)
            w(0)=1d0
            goto 9999
          endif
          do i=-1,4
            idf(i)=.false.
          enddo
          idf(0)=.true.
          idf(-1)=.true.
          wdd=1d0
        else ! hard (idpr=3)
          q2pmin(1)=q2nmin      !start at minimum Q2s
          q2pmin(2)=q2pmin(1)
          do i=-1,4
            idf(i)=.true.
          enddo
          if(factsat.le.0.)idf(-1)=.false.
          idf(0)=.false.        !soft not included
          wdd=1d0               !if 0, soft is not rescaled, 1. to scale soft
          ww=wwnk(2,n,k)        !probability used to select hard pomeron
        endif
        lndif(1)=idfpr(n,k).eq.1.or.idfpr(n,k).eq.2
        lndif(2)=idfpr(n,k).eq.1.or.idfpr(n,k).eq.3
        if(iremn.ge.2)then
          if(idf(2))idf(2)=(ivp(ip).gt.0).and.lndif(1)
          if(idf(3))idf(3)=(ivt(it).gt.0).and.lndif(2)
          if(idf(4))idf(4)=(ivp(ip)*ivt(it).gt.0)
     .                     .and.lndif(1).and.lndif(2)
        endif
        call ipoli(log(max(1e-6,q2pmin(1))),nq1,frac) 
        nq2=nq1
        inh=4       

      else   ! ============================= iqq=2 ============================= !Pom type already defined, recompute Q2s with Nhard

        inh=3    !hard only
        idf(0)=.false.
        wdd=1d0                 !not to count soft, it should be like that
        idf(-1)=idhpr(n,k).eq.-1
        do i=1,4
          idf(i)=idhpr(n,k).eq.i-1
        enddo
        q2pmin(1)=q2nmin        !start at minimum Q2s
        if(idf(-1))q2pmin(1)=q2mnval(2)!2.*qcmass**2 !q2mnval(2)      !do not start at minimum Q2s but maximum Pomsat
        q2pmin(2)=q2pmin(1)
        call ipoli(log(max(1e-6,max(q2pmin(1),q2pmin(2)))),nq1,frac) 
        nq2=nq1
        ww=wwnk(2,n,k)            !probability used to select hard pomeron

      endif  ! =================================== iqq ====================================

     
      if(iq2sat.ne.0.and.iscreen.ne.0)then 
      !===================================== iq2sat NE 0 ========================  !possiblity to use q2nmin only (as before)

      !###################!
      !     konn conn     !
      !###################!
      konn1=0
      konn2=0
      conn1=0.
      conn2=0.
      xnpom=0.
      do kxy=1,koll           !ckw21
       if(iproj(kxy).eq.ip)conn1=conn1+npr(inh,kxy)
       if(itarg(kxy).eq.it)conn2=conn2+npr(inh,kxy)
       if(iproj(kxy).eq.ip.and.npr(inh,kxy).gt.0)konn1=konn1+1
       if(itarg(kxy).eq.it.and.npr(inh,kxy).gt.0)konn2=konn2+1
      enddo 
      konn1=max(1,konn1)
      konn2=max(1,konn2)      
      conn12=max(1.,0.5*(conn1+conn2))

      xnpomn=float(max(1,npr(inh,k)))
      if(abs(iokoll).gt.0)then
         do kk=1,koll
           if(kk.ne.k)then
             nup=npr(inh,kk)
             xnpomn=xnpomn+float(nup)
           endif
         enddo
      endif
      if(.not.idf(0))then
        if(jval.eq.1)then
          zval=(conn1+conn2)/2
        elseif(jval.eq.2)then
          zval=log(zpom)
        else
          stop'ERROR 09062022b'
        endif
        binfun=1/deform(maproj,matarg,engy,zval,sngl(xh)) 
      else
        binfun=1d0
      endif
      fxxx=1.d0

      if(abs(iq2sat).ge.2)then
        con1=conn1
        con2=conn2
        xpome=npr(inh,k)
      endif  !abs(iq2sat).ge.2

      !+++++++++++++++++++++++++++++
      call defParamsWom(conn12,fuuu) 
      !+++++++++++++++++++++++++++++
      fxxx=fuuu*fxxx

      if(.not.idf(0))then
      fxxx=dble(csatur)*fxxx
      endif

      binfac=binfac*fxxx

      call eventvariset(20,binfac) !z20zevt
      !call getBinfacFit(binfacFit) 
      !call eventvariset(21,binfacFit) !z21zevt
      call eventvariset(21,binfac) !temporary, remove when binfacFit is defined

      fasym=1.
      cgcQ2s=dble(fasym)
     .    *PomIncMExactk(xm,float(k))/PomIncPExactk(xp,float(k))
     .     *(PomInc(PomIncPExactk,xp,float(k),con1,con2,xpome,2)
     .      /PomInc(PomIncMExactk,xm,float(k),con2,con1,xpome,2))

      cgcQ2s=cgcQ2s**asatur
      if(con1.eq.con2.and.rangen().gt.0.5)then
        cgcQ2s=1d0/cgcQ2s
      endif

      if(cgcQ2s.le.1d0)then  !projectile Q2s is larger
        iqp=1
        iqt=2
      else                     !target Q2s is larger
        cgcQ2s=1d0/cgcQ2s
        iqp=2
        iqt=1
      endif

      !##################!
      ! G=Geff procedure !
      !##################! 

c Initial conditions

      call ipoOm5Tab(xp,xm,bk(k))
      w(0)=om5Jk(nmx,k,0)
      w0=w(0)
      omt=0d0
      if(idf(1))omt=omt+om5Jk(nmx,k,1)*binfun
      if(idf(2))omt=omt+om5Jk(nmx,k,2)*binfun
      if(idf(3))omt=omt+om5Jk(nmx,k,3)*binfun
      if(idf(4))omt=omt+om5Jk(nmx,k,4)*binfun
      if(idf(-1))omt=omt+om5Jk(nmx,k,11)
      w1=omt
      eps=0.05
      if(iqq.eq.2.and.idf(-1))eps=0.25     !high precision not needed for Pom Sat and might lead to extreme Q2s values
      if(iqq.le.1)then
        if(w0.le.0d0)then
          w(0)=1d0
          goto 9999        
        endif
        eps=eps*max(1d0,w0/(omt+w0)) !set precision at level of 5% of soft contribution (without screening)
      endif
      if(idf(0))omt=omt+wdd*w0
      wsh=binfac*omt+(1d0-wdd)*w0

      oooooo=om5Jk(n,k,11)

c look for Q2s which fit best the used omega

      if(ish.ge.8)
     .  write(ifch,*)'0',iqq,iqp,iqt,q2pmin(1),q2pmin(2),wsh,ww,cgcQ2s
     .               ,oooooo,nq1,nq2,binfac,xnpom,idf
      lpass=.false.
      lsat=.false.
      if(wsh.lt.ww)then
        if(iqq.eq.1)then
          q2pmin(1)=q2nmin
          q2pmin(2)=q2pmin(1)
        endif
        if(ish.ge.8)
     .  write(ifch,*)'1 skip'
        call incrementNskipNpassWomTy(1)
      else ! wsh.ge.ww
       do while (wsh.gt.ww*(1d0-eps).and.nq2.lt.maxq2mx)
        lpass=.true.
        call incrementNskipNpassWomTy(2)
        nq1=nq2
        nq2=nq2+1
        q2pmin(iqp)=max(q2nmin,q2mnval(nq2))
        q2pmin(iqt)=max(q2nmin,sngl(cgcQ2s)*q2pmin(iqp))
        call ipoOm5Tab(xp,xm,bk(k))
        omt=0d0
        if(idf(1))omt=omt+om5Jk(nmx,k,1)*binfun  !---------------------------------------!
        if(idf(2))omt=omt+om5Jk(nmx,k,2)*binfun  !om5Jk depends on global variable q2pmin!
        if(idf(3))omt=omt+om5Jk(nmx,k,3)*binfun  !---------------------------------------!
        if(idf(4))omt=omt+om5Jk(nmx,k,4)*binfun
        w1234=omt*binfac
        if(idf(-1))omt=omt+om5Jk(nmx,k,11) 
        w1=omt
        w0=om5Jk(nmx,k,0)
        if(idf(0))omt=omt+w0*wdd
        w5a=0
        if(idf(-1))w5a=om5Jk(nmx,k,11)*binfac
        w0a=0
        if(idf(0))w0a=w0*wdd*binfac
        w0a=w0a+(1d0-wdd)*w0
         wsh=binfac*omt+(1d0-wdd)*w0
cKW2108        if(iqq.eq.1)then!+++++
cKW2108          write(*,'(a,2i2,1x,2i4,$)')'G=Geff-iter++',iqq,idpr(n,k),k,n
cKW2108          write(*,'(f8.2,3x,3f8.2,$)')q2pmin(iqp),w0a,w5a,w1234
cKW2108          write(*,'(f16.2,f8.2,3x,f8.2)')wsh,ww,binfac 
cKW2108        endif !+++++++++++++++
         if(ish.ge.8)
     .    write(ifch,*)'1',iqq,q2pmin(iqp),q2pmin(iqt),wsh,ww,nq1,nq2
     .           ,binfac*om5Jk(nmx,k,0),lndif(iqp),lndif(iqt),idf
       enddo
       if(wsh.gt.ww)then
        nq1=maxq2mx
       endif
cKW2108        if(iqq.eq.1)then!+++++
cKW2108          write(*,'(a)')'G=Geff-iter++'
cKW2108        endif !+++++++++++++++
      endif

c fine tune value of Q2s between q2mnval(nq1) and q2mnval(nq2)

      if(lpass.and.nq2.gt.1.and.nq1.lt.maxq2mx)then
       q2m1=max(q2nmin,q2mnval(nq1))
       q2m2=max(q2nmin,q2mnval(nq2))
       do while ((abs(ww-wsh)/ww.gt.eps.and.abs(q2m1-q2m2).gt.0.5)
     .          .or.(wsh.gt.ww*(1d0+eps).and.nq2.le.maxq2mx))
        q2pmin(iqp)=max(q2nmin,0.5*(q2m2+q2m1))
        q2pmin(iqt)=max(q2nmin,sngl(cgcQ2s)*q2pmin(iqp))
        call ipoOm5Tab(xp,xm,bk(k))
        omt=0d0
        if(idf(1))omt=omt+om5Jk(nmx,k,1)*binfun  !---------------------------------------!
        if(idf(2))omt=omt+om5Jk(nmx,k,2)*binfun  !om5Jk depends on global variable q2pmin!
        if(idf(3))omt=omt+om5Jk(nmx,k,3)*binfun  !---------------------------------------!
        if(idf(4))omt=omt+om5Jk(nmx,k,4)*binfun
        if(idf(-1))omt=omt+om5Jk(nmx,k,11) 
        w1=omt
        w0=om5Jk(nmx,k,0)
        if(idf(0))omt=omt+w0*wdd
        wsh=binfac*omt+(1d0-wdd)*w0
        if(ish.ge.8)
     .  write(ifch,*)'2',q2pmin(iqp),q2m1,q2m2,wsh,ww
     .  ,om5Jk(nmx,k,11)
        if(wsh.gt.ww)then
          q2m1=q2pmin(iqp)
          if(abs(q2m2-q2m1).le.1e-4)then   !to avoid infinite loop
            nq2=nq2+1
            if(nq2.le.maxq2mx)q2m2=max(q2nmin,q2mnval(nq2))
          endif
        else
          q2m2=q2pmin(iqp)
        endif
       enddo
      else
        q2m1=0.
        q2m2=0.
      endif
      if(max(nq1,nq2).eq.maxq2mx.and.wsh.gt.ww.and.w0.lt.ww)then
        if(iqq.eq.2.and.idf(-1))then
          return
        else
          ncount=ncount+1 
          if(ncount.eq.1.or.ncount.eq.10.or.ncount.eq.100
     .    .or.ncount.eq.1000)then
            lcount=log10(1.*ncount)+1
            write(ifmt,
     .      '(a,i1,a,i4,3f7.2,2f8.0,2f7.2,2f7.4,2f8.0,f7.4,6l2)')
     .      'WARNING ',lcount,' ems WomTy maxQ2s too low'
     .      ,iqq,xpr(nmx,k),binfac,cgcq2s,q2pmin(1),q2pmin(2),wsh/ww,ww
     .      ,w0,om5Jk(nmx,k,11),q2m1,q2m2,eps,idf
          endif
          if(maxq2mx.lt.maxq2mn)then
            write(ifmt,*)'WARNING : maxq2mx is too low in setBasics ! '
          endif
        endif
      endif
      
      else !===================================== iq2sat EQ 0 =====================================

        q2pmin(1)=csatur*q2nmin
        q2pmin(2)=csatur*q2nmin
        xp=xppr(n,k)
        xm=xmpr(n,k)
        call ipoOm5Tab(xp,xm,bk(k))
        omt=0d0
        if(idf(1))omt=omt+om5Jk(n,k,1)*binfun
        if(idf(2))omt=omt+om5Jk(n,k,2)*binfun
        if(idf(3))omt=omt+om5Jk(n,k,3)*binfun
        if(idf(4))omt=omt+om5Jk(n,k,4)*binfun
        if(idf(-1))omt=omt+om5Jk(n,k,11) 
        w1=omt
        w0=om5Jk(n,k,0)             
        if(idf(0))omt=omt+w0*wdd
          wsh=omt+(1d0-wdd)*w0
          binfac=ww/wsh

      endif    !======================================= iq2sat ====================================               

      if(iqq.eq.1)then !================================= iqq=1 ========================================
        
      omt=0d0
      do i=-1,4           !only hard part
        w(i)=0
        if(idf(i))then
          ii=i
          if(i.eq.-1)ii=11
          if(i.gt.0)then                      !---------------------------------------!
            w(i)=binfac*om5Jk(n,k,ii)*binfun  !om5Jk depends on global variable q2pmin!
            omt=omt+w(i)                      !---------------------------------------!
          elseif(i.eq.0)then
            w(0)=om5Jk(n,k,0)
            if(wdd.gt.0d0)w(0)=binfac*w(0)
            omt=omt+w(i)
          else ! -1
            w(i)=om5Jk(n,k,ii) 
          endif
        endif
      enddo
      w(-1)=binfac*w(-1)
      w(5)=w(-1)
      wsum=w(0)+w(1)+w(2)+w(3)+w(4)+w(5)
      ichoi=0 ! 1,2,3
      if(ichoi.eq.0)then
        iw=0
        if(w(5).gt.1d-10)iw=5 
        wsum=wsum-w(iw)
      elseif(ichoi.eq.1)then
        iw=0
        if(w(5).gt.1d10)iw=5 
        wsum=wsum-w(iw)
      elseif(ichoi.eq.2)then
        iw=0
        if(w(5).gt.1d-10)iw=5 
        wsum=omt
      elseif(ichoi.eq.3)then
        iw=0
        if(w(5).gt.1d+10)iw=5 
        wsum=omt
      else
        stop'ERROR 22072021b'
      endif
      w0old=w(0)
      w5old=w(5)
      if(iw.eq.5)then
        w(5)=max(0d0,ww-wsum)    
      elseif(iw.eq.0)then
        w(0)=max(0d0,ww-wsum)    
      else
        stop'ERROR 22072021'
      endif

      if(ish.ge.7)write(ifch,*)"w(i)",w,om5Jk(n,k,0)
     . ,om5Jk(n,k,11),ww

      endif !================================= iqq=1 end =========================================

      if(ish.ge.6)write(ifch,*)'WomTy : '
     .,iqq,n,k,xp,xm,binfac,cgcq2s,xpome,con1,con2
     .,q2pmin(1),q2pmin(2),wsh/ww,bk(k),idf

 9999 continue
      w(-1)=w(5)
c      if(iqq.eq.1)then 
c       if(idpr(n,k).eq.1)then  !.not.lpass)then
c       write(*,'(a,i8,4f8.2,i8,4x,l7,f9.2)')'TEST-WomTy'
c     . ,idpr(n,k),w(0),w(5),w(1)+w(2)+w(3)+w(4),ww,nint(conn1+conn2)
c     . ,lpass,q2pmin(iqp)
c       endif
c      endif 
      q2kmin(1,n,k)=q2pmin(1)
      q2kmin(2,n,k)=q2pmin(2)

      return
      end


c###############################################################################################
c###############################################################################################
c          New deformation function method (New DF method)
c###############################################################################################
c###############################################################################################


c----------------------------------------------------------------------
      function deform(maproj,matarg,engy,z,x) !fit of deformation function (x,y)
c----------------------------------------------------------------------
      u= -log(x)/20
      w= defoW(
     .piDeform(1,z,maproj,matarg,engy),piDeform(2,z,maproj,matarg,engy),
     .piDeform(3,z,maproj,matarg,engy),piDeform(4,z,maproj,matarg,engy),
     .piDeform(5,z,maproj,matarg,engy),u)
      y= exp(20*w)
      deform=y
      end
      !------------------+------------------+
      !  u = -log(x)/20  ! x = exp(-20*u)   !
      !  w = log(y)/20   ! y = exp(20*w)    !
      !------------------+------------------+ 
      function defoW(a1,a2,a3,a4,a5,u) !fit of deformation function (u,w)
      v=abs(u-a1)/0.5
      w=a2-v**a3*a4
      if(u.gt.a1)w=max(0.,w) 
      w=max(a5,w)
      defoW=w
      end

c----------------------------------------------------------------------
      integer function igetJval(maproj,matarg,engy)
c----------------------------------------------------------------------
      parameter (KDIM=8,NSYSDIM=100,MDIM=5)  
      character*500 path
      character *80 line
      data ncount/0/
      save ncount,jval
      ncount=ncount+1
      path=' ' 
      if(ncount.eq.1)then   
        jval=0
        call getEposPath(path,length)
        open(105,file= path(1:length)//'/src/KWt/deform.dt'
     .  ,status='old')
        rewind(105)
        it=0
        do 
          read(105,*,end=777)maprojR,matargR,engyR,jvalR
          if(((maproj.eq.maprojR.and.matarg.eq.matargR)
     .    .or.(maproj.eq.matargR.and.matarg.eq.maprojR))
     .    .and.nint(engy*1000).eq.nint(engyR*1000))then
            jval=jvalR
          endif
          do K=KDIM,1,-1                
            read(105,*)line
          enddo
        enddo
 777    close(105)
      endif
      if(jval.eq.0)stop'ERROR System not found in deform.dt'
      igetJval=jval
      end

c----------------------------------------------------------------------
      function piDeform(M,z,maproj,matarg,engy) !polynomial interpolation (linear)
c----------------------------------------------------------------------
      parameter (KDIM=8,NSYSDIM=100,MDIM=5)  
      real zar(KDIM),war(KDIM)
      real zx(NSYSDIM,KDIM),wx(NSYSDIM,KDIM,MDIM),engyx(NSYSDIM)
      integer maprojx(NSYSDIM),matargx(NSYSDIM)
      character*500 path
      save zx,wx,maprojx,matargx,engyx
      data ncount/0/
      save ncount
      ncount=ncount+1
      if(ncount.eq.1)then   
        call getEposPath(path,length)
        write(ifmt,'(2a)')'load ',path(1:length)//'/src/KWt/deform.dt'
        open(105,file= path(1:length)//'/src/KWt/deform.dt'
     .  ,status='old')
        !open(105,file='../../util/deform.dt',status='old')
          rewind(105)
          it=0
          do 
            read(105,*,end=777)maprojR,matargR,engyR
            it=it+1
            if(it.gt.NSYSDIM )stop'NSYSDIM too small'
            maprojx(it)=maprojR
            matargx(it)=matargR
            engyx(it)  =engyR
            do K=KDIM,1,-1    
              read(105,*)Kx,zx(it,K),(wx(it,K,Mx),Mx=1,MDIM)
            enddo
          enddo
 777    close(105)
      endif
      nsys=0
      do it=1,NSYSDIM
        if(
     .     (maproj.eq.maprojx(it).and.matarg.eq.matargx(it).or.
     .      matarg.eq.maprojx(it).and.maproj.eq.matargx(it) )
     .     .and.nint(engy*1000).eq.nint(engyx(it)*1000)
     .    )nsys=it
      enddo 
      if(nsys.eq.0)stop'System not in deform table' !eventually we interpolate here
      do K=1,KDIM    
        zar(K)=zx(nsys,K)
        war(K)=wx(nsys,K,M)
      enddo
      if(z.lt.zar(1))then
        w=war(1)
      else
        iz=1
        do while(iz.lt.KDIM-1.and.zar(iz+1).gt.0..and.zar(iz+1).lt.z)
          iz=iz+1
        enddo
        if(iz.lt.KDIM.and.zar(iz+1).eq.0.)iz=iz-1
        z1=zar(iz)
        z2=zar(iz+1)
        w1=war(iz)
        w2=war(iz+1)
        w=w1+(z-z1)/(z2-z1)*(w2-w1)
      endif
      piDeform=w
      end 


c###############################################################################################
c###############################################################################################


c-----------------------------------------------------------------------
      subroutine StoRe(imod)
c-----------------------------------------------------------------------
c Store Remnant configuration (imod=1) before shuffle  to restore the
c initial configuration (imod=-1) in case of problem.
c-----------------------------------------------------------------------

#include "aaa.h"
#include "ems.h"

      if(imod.eq.1)then

c       initialize projectile

        do i=1,maproj
          iepst(i)=iep(i)
          xppst(i)=xpp(i)
          xmpst(i)=xmp(i)
          xposst(i)=xpos(i)
        enddo

c       initialize target

        do j=1,matarg
          ietst(j)=iet(j)
          xmtst(j)=xmt(j)
          xptst(j)=xpt(j)
          xtosst(j)=xtos(j)
        enddo

      elseif(imod.eq.-1)then

c       restore projectile

        do i=1,maproj
          iep(i)=iepst(i)
          xpp(i)=xppst(i)
          xmp(i)=xmpst(i)
          xpos(i)=xposst(i)
        enddo

c       restore target

        do j=1,matarg
          iet(j)=ietst(j)
          xmt(j)=xmtst(j)
          xpt(j)=xptst(j)
          xtos(j)=xtosst(j)
        enddo

      else

        call utstop('Do not know what to do in StoRe.&')

      endif

      return
      end

c-----------------------------------------------------------------------
      subroutine ProSeF(k,n,iret)
c-----------------------------------------------------------------------
c starting from string properties as already determined in EMS,
c one determines string end flavors
c by checking compatibility with remnant masses.
c strings are written to /cems/ and then to /cptl/
c remnant ic is updated (icproj,ictarg)
c------------------------------------------------------------------------

#include "aaa.h"
#include "ems.h"
#include "sem.h"

      double precision plc,s,pstg,pend
      common/cems5/plc,s
      common/cems/pstg(5,2),pend(4,4),idend(4)
      common/emsptl/nppr(npommx,kollmx),npproj(mamx),nptarg(mamx)
      integer icp(2),ict(2),ic(2),icp1(2),icp2(2),icm1(2),icm2(2)
      integer jcp(nflav,2),jct(nflav,2),jcpv(nflav,2),jctv(nflav,2)
      integer jcp1(nflav,2),jcp2(nflav,2),jcm1(nflav,2),jcm2(nflav,2)
      common/col3/ncol,kolpt,ncoli /cfacmss/facmss /cts/its

      integer getAccumJerr
c     entry
c     -----

      iret=0

      if(idpr(n,k).eq.0.or.ivpr(n,k).eq.0)return
      if(idpr(n,k).eq.3)return
      call utpri('ProSeF',ish,ishini,5)
      np=nppr(n,k)
      idptl(np)=abs(idptl(np))
      ip=iproj(k)
      it=itarg(k)
      if(ish.ge.5)then
          write(ifch,*)'soft Pomeron'
          write(ifch,*)'k:',k,'  n:',n,'  ip:',ip,'  it:',it
     .              ,'  id:',idpr(n,k)
      endif


c         string ends

          pend(1,1)=xxp1pr(n,k)
          pend(2,1)=xyp1pr(n,k)
          pend(3,1)=xp1pr(n,k)*plc/2d0
          pend(4,1)=dsqrt(pend(1,1)**2+pend(2,1)**2+pend(3,1)**2)
          pend(1,2)=xxp2pr(n,k)
          pend(2,2)=xyp2pr(n,k)
          pend(3,2)=xp2pr(n,k)*plc/2d0
          pend(4,2)=dsqrt(pend(1,2)**2+pend(2,2)**2+pend(3,2)**2)
          pend(1,4)=xxm1pr(n,k)
          pend(2,4)=xym1pr(n,k)
          pend(3,4)=-xm1pr(n,k)*plc/2d0
          pend(4,4)=dsqrt(pend(1,4)**2+pend(2,4)**2+pend(3,4)**2)
          pend(1,3)=xxm2pr(n,k)
          pend(2,3)=xym2pr(n,k)
          pend(3,3)=-xm2pr(n,k)*plc/2d0
          pend(4,3)=dsqrt(pend(1,3)**2+pend(2,3)**2+pend(3,3)**2)

c         strings

          pstg(1,1)=xxp1pr(n,k)+xxm2pr(n,k)
          pstg(2,1)=xyp1pr(n,k)+xym2pr(n,k)
          pstg(3,1)=(xp1pr(n,k)-xm2pr(n,k))*plc/2d0
          pstg(4,1)=(xp1pr(n,k)+xm2pr(n,k))*plc/2d0
          pstg(5,1)=dsqrt((pstg(4,1)-pstg(3,1))*(pstg(4,1)+pstg(3,1))
     &                   -pstg(1,1)**2-pstg(2,1)**2)
          pstg(1,2)=xxp2pr(n,k)+xxm1pr(n,k)
          pstg(2,2)=xyp2pr(n,k)+xym1pr(n,k)
          pstg(3,2)=(xp2pr(n,k)-xm1pr(n,k))*plc/2d0
          pstg(4,2)=(xp2pr(n,k)+xm1pr(n,k))*plc/2d0
          pstg(5,2)=dsqrt((pstg(4,2)-pstg(3,2))*(pstg(4,2)+pstg(3,2))
     &                   -pstg(2,2)**2-pstg(1,2)**2)

c         initialize

          ntry=0
  777     ntry=ntry+1
          if(ntry.gt.100)goto 1001

          if(iremn.ge.2)then    !uses precalculated flavors
            do i=1,2
              icp(i)=icproj(i,ip)
              ict(i)=ictarg(i,it)
            enddo
            call iddeco(icp,jcpv)
            call iddeco(ict,jctv)
            do j=1,2
              do i=1,nrflav
                jcp(i,j)=jcpref(i,j,ip)
                jct(i,j)=jctref(i,j,it)
              enddo
              do i=nrflav+1,nflav
                jcp(i,j)=0
                jct(i,j)=0
                jcpv(i,j)=0
                jctv(i,j)=0
               enddo
            enddo
            do i=1,2
              icp(i)=icproj(i,ip)
              ict(i)=ictarg(i,it)
            enddo
          else
            do i=1,2
              icp(i)=icproj(i,ip)
              ict(i)=ictarg(i,it)
            enddo
            call iddeco(icp,jcp)
            call iddeco(ict,jct)
            do j=1,2
              do i=1,nflav
                jcpv(i,j)=0
                jctv(i,j)=0
              enddo
            enddo
          endif
          do i=1,2
           icp1(i)=0
           icp2(i)=0
           icm1(i)=0
           icm2(i)=0
           do j=1,nflav
            jcp1(j,i)=0
            jcp2(j,i)=0
            jcm1(j,i)=0
            jcm2(j,i)=0
           enddo
          enddo
          idpj0=idtr2(icp)
          idtg0=idtr2(ict)
          do j=1,4
           idend(j)=0
          enddo

          if(ish.ge.7)then
            write(ifch,'(a,3x,6i3,3x,6i3,i9)')' proj: '
     *     ,jcp,idpj0
            write(ifch,'(a,6i3,3x,6i3)')' proj val:  ',jcpv
          endif
          if(ish.ge.7)then
            write(ifch,'(a,3x,6i3,3x,6i3,i9)')' targ: '
     *    ,jct,idtg0
            write(ifch,'(a,6i3,3x,6i3)')' targ val:  ',jctv
          endif

c         determine string flavors

          call fstrfl(jcp,jct,jcpv,jctv,icp1,icp2,icm1,icm2
     *                ,idp1pr(n,k),idp2pr(n,k),idm1pr(n,k),idm2pr(n,k)
     *                                   ,idsppr(n,k),idstpr(n,k),iret)
          if(iret.ne.0)goto 1002

c         check mass string 1

          ic(1)=icp1(1)+icm2(1)
          ic(2)=icp1(2)+icm2(2)
          if(ic(1).gt.0.or.ic(2).gt.0)then
           am=sngl(pstg(5,1))
           idr=0
c           if(idpr(n,k).eq.2)then    !test resonance case for Reggeon
c             id=idtra(ic,0,0,0)
c             call idres(id,am,idr,iadj,1)
c           endif
           if(idr.eq.0)then     !string
             call iddeco(icp1,jcp1)
             call iddeco(icm2,jcm2)
             ammns=utamnx(jcp1,jcm2)
             if(ish.ge.7)write(ifch,'(a,2i7,2e12.3)')
     *           ' string 1 - ic,mass,min.mass:',ic,am,ammns
             if(am.lt.ammns*facmss)then
               goto 777         !avoid virpom
             endif
             idend(1)=idtra(icp1,0,0,0)
             idend(3)=idtra(icm2,0,0,0)
             if(ish.ge.7)write(ifch,'(a,2i6)') ' string 1 - SE-ids:'
     *            ,idend(1),idend(3)
           else      !resonance
             idend(1)=idr
             idend(3)=0.
           endif
          endif

c         check mass string 2

          if(idpr(n,k).eq.6)then       !low mass CD

          if(ish.ge.5)then
          write(ifch,'(a,i10)')' reg:   '
     *    ,idptl(np)
          write(ifch,'(a,2i5)')' str 1: '
     *    ,idend(1),idend(3)
          endif

c         write strings to /cptl/

          its=idp1pr(n,k)+idm2pr(n,k)
          call fstrwr(1,1,3,k,n)

          else                        !Pomerons

          ic(1)=icp2(1)+icm1(1)
          ic(2)=icp2(2)+icm1(2)
          if(ic(1).gt.0.or.ic(2).gt.0)then
           am=sngl(pstg(5,2))
           call iddeco(icp2,jcp2)
           call iddeco(icm1,jcm1)
           ammns=utamnx(jcp2,jcm1)
           if(ish.ge.7)write(ifch,'(a,2i7,2e12.3)')
     *           ' string 2 - ic,mass,min.mass:',ic,am,ammns
           if(am.lt.ammns*facmss)then
             goto 777  !avoid virpom
           endif
           idend(2)=idtra(icp2,0,0,0)
           idend(4)=idtra(icm1,0,0,0)
           if(ish.ge.7)write(ifch,'(a,2i6)') ' string 2 - SE-ids:'
     *      ,idend(2),idend(4)
          endif

          if(ish.ge.5)then
          write(ifch,'(a,i10)')' pom:   '
     *    ,idptl(np)
          write(ifch,'(a,2i5)')' str 1: '
     *    ,idend(1),idend(3)
          write(ifch,'(a,2i5)')' str 2: '
     *    ,idend(2),idend(4)
          endif

c         write strings to /cptl/

          its=idp1pr(n,k)+idm2pr(n,k)
          call fstrwr(1,1,3,k,n)
          its=idp2pr(n,k)+idm1pr(n,k)
          call fstrwr(2,2,4,k,n)

          endif

c         update remnant ic

c determine icp,ict
c Similar process for hard pomeron in epos-rsh !!!!

          if(iremn.ge.2)then    !uses precalculated flavors

            do j=1,2
              do i=1,nrflav
                jcpref(i,j,ip)=jcp(i,j)
                jctref(i,j,it)=jct(i,j)
              enddo
            enddo
            call idenco(jcpv,icp,iret)
            if(iret.ne.0)goto 1002
            call idenco(jctv,ict,iret)
            if(iret.ne.0)goto 1002
            do i=1,2
              icproj(i,ip)=icp(i)
              ictarg(i,it)=ict(i)
            enddo
            if(ish.ge.5)then
              write(ifch,'(a,6i3,3x,6i3)')' proj:  ',jcp
              write(ifch,'(a,6i3,3x,6i3)')' proj val:  ',jcpv
              write(ifch,'(a,6i3,3x,6i3)')' targ:  ',jct
              write(ifch,'(a,6i3,3x,6i3)')' targ val:  ',jctv
            endif

          else

            call idenco(jcp,icp,iret)
            if(iret.ne.0)goto 1002
            call idenco(jct,ict,iret)
            if(iret.ne.0)goto 1002
            do i=1,2
              icproj(i,ip)=icp(i)
              ictarg(i,it)=ict(i)
            enddo
            if(ish.ge.5)then
              write(ifch,'(a,2i7,1x,a)')' proj:  '
     *             ,(icp(l),l=1,2)
              write(ifch,'(a,2i7,1x,a)')' targ:  '
     *             ,(ict(l),l=1,2)
            endif

          endif

c     exit
c     ----

1000  continue
      call utprix('ProSeF',ish,ishini,5)
      return

 1002 call setAccumJerr(1,getAccumJerr(1)+1)         ! > 9 quarks per flavor attempted.
 1001 iret=1
      if(ish.ge.5)write(ifch,'(a)')'Problem in ProSeF ... '
      goto 1000

      end

c-----------------------------------------------------------------------
      subroutine fstrwr(j,ii,jj,k,n)
c-----------------------------------------------------------------------
c take pstg(5,j),pend(4,ii),idend(ii),pend(4,jj),idend(jj)  (/cems/)
c and write it to /cptl/
c-----------------------------------------------------------------------
c  j:     string 1 or 2
c  ii,jj: string end (1,2: proj; 3,4: targ)
c  k:     current collision
c  n:     current pomeron
c-----------------------------------------------------------------------

#include "aaa.h"
#include "ems.h"

      double precision pstg,pend,ptt3!,utpcmd
      common/cems/pstg(5,2),pend(4,4),idend(4)
      common/emsptl/nppr(npommx,kollmx),npproj(mamx),nptarg(mamx)
      double precision  pp(4)
      common/cts/its

      call utpri('fstrwr',ish,ishini,7)

      if(idend(ii).ne.0.and.idend(jj).ne.0)then

        am1=0.
        am2=0.
        ptt3=0.5d0*pstg(5,j)

       call utlob2(1,pstg(1,j),pstg(2,j),pstg(3,j),pstg(4,j),pstg(5,j)
     * ,pend(1,ii),pend(2,ii),pend(3,ii),pend(4,ii),20)
       pp(1)=0d0
       pp(2)=0d0
       pp(3)=ptt3!.5d0*pstg(5,j)
       pp(4)=sqrt(ptt3*ptt3+dble(am1*am1))!.5d0*pstg(5,j)
       call utrot2
     * (-1,pend(1,ii),pend(2,ii),pend(3,ii),pp(1),pp(2),pp(3))
       call utlob2(-1,pstg(1,j),pstg(2,j),pstg(3,j),pstg(4,j),pstg(5,j)
     * ,pp(1),pp(2),pp(3),pp(4),21)

       npom=nppr(n,k)
       if(ifrptl(1,npom).eq.0)ifrptl(1,npom)=nptl+1
       ifrptl(2,npom)=nptl+2
       istptl(npom)=31

       nptl=nptl+1
       pptl(1,nptl)=sngl(pp(1))
       pptl(2,nptl)=sngl(pp(2))
       pptl(3,nptl)=sngl(pp(3))
       pptl(4,nptl)=sngl(pp(4))
       pptl(5,nptl)=am1 !0.
       istptl(nptl)=20
       iorptl(nptl)=npom
       jorptl(nptl)=0
       ifrptl(1,nptl)=0
       ifrptl(2,nptl)=0
       xorptl(1,nptl)=coordpr(1,n,k)
       xorptl(2,nptl)=coordpr(2,n,k)
       xorptl(3,nptl)=coord(3,k)
       xorptl(4,nptl)=coord(4,k)
       tivptl(1,nptl)=xorptl(4,nptl)
       tivptl(2,nptl)=xorptl(4,nptl)
       idptl(nptl)=idend(ii)
       ityptl(nptl)=ityptl(npom)+j
       itsptl(nptl)=its
       rinptl(nptl)=-9999
       qsqptl(nptl)=0.
       zpaptl(1,nptl)=zpaptl(1,npom)
       zpaptl(2,nptl)=zpaptl(2,npom)

       nptl=nptl+1
       do i=1,4
        pptl(i,nptl)=sngl(pstg(i,j))-pptl(i,nptl-1)
       enddo
       pptl(5,nptl)=am2!0.

       istptl(nptl)=20
       iorptl(nptl)=nppr(n,k)
       jorptl(nptl)=0
       ifrptl(1,nptl)=0
       ifrptl(2,nptl)=0
       xorptl(1,nptl)=coordpr(1,n,k)
       xorptl(2,nptl)=coordpr(2,n,k)
       xorptl(3,nptl)=coord(3,k)
       xorptl(4,nptl)=coord(4,k)
       tivptl(1,nptl)=xorptl(4,nptl)
       tivptl(2,nptl)=xorptl(4,nptl)
       idptl(nptl)=idend(jj)
       ityptl(nptl)=ityptl(npom)+j
       itsptl(nptl)=its
       rinptl(nptl)=-9999
       qsqptl(nptl)=0.
       zpaptl(1,nptl)=zpaptl(1,npom)
       zpaptl(2,nptl)=zpaptl(2,npom)

       if(ish.ge.7)then
        write(ifch,100)' kink:',(pptl(l,nptl-1),l=1,4),idptl(nptl-1)
        write(ifch,100)' kink:',(pptl(l,nptl),l=1,4),idptl(nptl)
       endif

      elseif(idend(ii).ne.0.and.idend(jj).eq.0)then

c resonance

       npom=nppr(n,k)
       if(ifrptl(1,npom).eq.0)ifrptl(1,npom)=nptl+1
       ifrptl(2,npom)=nptl+1
       istptl(npom)=31

       nptl=nptl+1
       idptl(nptl)=idend(ii)
       pptl(1,nptl)=sngl(pstg(1,j))
       pptl(2,nptl)=sngl(pstg(2,j))
       pptl(3,nptl)=sngl(pstg(3,j))
       pptl(4,nptl)=sngl(pstg(4,j))
       pptl(5,nptl)=sngl(pstg(5,j))
       istptl(nptl)=0
       iorptl(nptl)=npom
       jorptl(nptl)=0
       ifrptl(1,nptl)=0
       ifrptl(2,nptl)=0
       xorptl(1,nptl)=coordpr(1,n,k)
       xorptl(2,nptl)=coordpr(2,n,k)
       xorptl(3,nptl)=coord(3,k)
       xorptl(4,nptl)=coord(4,k)
       tivptl(1,nptl)=coord(4,k)
       call idtau(idptl(nptl),pptl(4,nptl),pptl(5,nptl),taugm)
       tivptl(2,nptl)=tivptl(1,nptl)+taugm*(-alog(rangen()))
       ityptl(nptl)=ityptl(npom)+2+j
       itsptl(nptl)=its
       rinptl(nptl)=-9999
       qsqptl(nptl)=0.
       zpaptl(1,nptl)=zpaptl(1,npom)
       zpaptl(2,nptl)=zpaptl(2,npom)
       if(ish.ge.7)then
        write(ifch,100)'  res:',(pptl(l,nptl),l=1,4),idptl(nptl)
       endif
      elseif(idend(ii).eq.0.and.idend(jj).eq.0)then
       goto 1000
      else
       write(ifmt,*)'----->  j,ii,jj,idend(ii),idend(jj):'
       write(ifmt,*)'----->',j,ii,jj,idend(ii),idend(jj)
       stop'####### ERROR in fstrwr #######'
      endif

  100 format(a,4e9.3,i5)

1000  continue
      call utprix('fstrwr',ish,ishini,7)
      return
      end

c-----------------------------------------------------------------------
      function idtr2(ic)
c-----------------------------------------------------------------------
c transforms ic to id such that only hadrons have nonzero id
c-----------------------------------------------------------------------
      parameter (nidt=30)
      integer idt(3,nidt),ic(2)
      data idt/
     * 100000,100000, 110   ,100000,010000, 120   ,010000,010000, 220
     *,100000,001000, 130   ,010000,001000, 230   ,001000,001000, 330
     *,100000,000100, 140   ,010000,000100, 240   ,001000,000100, 340
     *,000100,000100, 440
     *,300000,000000,1111   ,210000,000000,1120   ,120000,000000,1220
     *,030000,000000,2221   ,201000,000000,1130   ,111000,000000,1230
     *,021000,000000,2230   ,102000,000000,1330   ,012000,000000,2330
     *,003000,000000,3331   ,200100,000000,1140   ,110100,000000,1240
     *,020100,000000,2240   ,101100,000000,1340   ,011100,000000,2340
     *,002100,000000,3340   ,100200,000000,1440   ,010200,000000,2440
     *,001200,000000,3440   ,000300,000000,4441/

      idtr2=0
      if(ic(1).eq.0.and.ic(2).eq.0)then
       if(rangen().ge.0.5)then
        idtr2=110
        ic(1)=100000
        ic(2)=100000
       else
        idtr2=220
        ic(1)=10000
        ic(2)=10000
       endif
       return
      endif
      do 1 i=1,nidt
       if(ic(2).eq.idt(1,i).and.ic(1).eq.idt(2,i))idtr2=-idt(3,i)
       if(ic(1).eq.idt(1,i).and.ic(2).eq.idt(2,i))idtr2=idt(3,i)
1     continue
      return
      end

c----------------------------------------------------------------------
      subroutine emsini(e,idpji,idtgi,iset)
c----------------------------------------------------------------------
c  energy-momentum sharing initializations
c----------------------------------------------------------------------
#include "aaa.h"
#include "ems.h"
#include "sem.h"
#include "par.h"
      double precision d,a,b,plc,s,amd,xvpr,xdm,alptr,xdm2
      common/cemsr5/alptr(0:1,0:5)
      common/cems5/plc,s
      common/cems10/a(0:ntypmx),b(0:ntypmx),d(0:ntypmx)
      common/ems6/ivp0,iap0,idp0,isp0,ivt0,iat0,idt0,ist0
      common/cems13/xvpr(0:3)
      common/geom/rmproj,rmtarg,bmax,bkmx


c parameter test

      if(nflavems.lt.nrflav)
     &   call utstop("nflavems<nrflav : change it in epos-ems !&")


c abreviations

      plc=dble(e)
      s=plc**2
      amd=dble(exmass) !large enough in case of strangeness in string end

c alpha (0=0, 1=s, 2=v, 4=d, 8=f)

      a(0)=0d0
      a(1)=dble(alpsea)
      a(2)=dble(alpval)
      a(3)= 0.0d0
      a(4)=dble(alpdiq)
      a(5)=dble(a(4))
      a(6)= 0.0d0
      a(7)= 0.0d0
      a(8)=dble(a(2))
      a(9)= 0.0d0

c beta (0=0, 1=s, 2=v, 4=d, 8=f) !to be updated if needed again

c      b(0)=0.0d0
c      b(1)=dble(-alpqua)
c      b(2)=dble(-alpqua)
c      b(3)=0.0d0
c      b(4)=0.0d0
c      b(5)=0.0d0
c      b(6)=0.0d0
c      b(7)=0.0d0
c      b(8)=dble(-alpqua)
c      b(9)=0.0d0


c alpha_trailing and beta_trailing (0=meson, 1=baryon;
c                                   0=no excit, 1=nondiffr, 2=diffr,
c                                   3=nondiffr split, 5=diffr split)

      alptr(0,0)=0.0d0
      alptr(0,1)=dble(alpdi(1))
      alptr(0,2)=dble(alppom)
      alptr(0,3)=dble(alpndi(1))
      alptr(0,4)=dble(0.5*(sloreg(1)+sloreg(3)))
      alptr(0,5)=dble(alphigh(1))
      alptr(1,0)=0.0d0
      alptr(1,1)=dble(alpdi(2))
      alptr(1,2)=dble(alppom)
      alptr(1,3)=dble(alpndi(2))
      alptr(1,4)=dble(sloreg(2))
      alptr(1,5)=dble(alphigh(2))

c minimal string masses ( i+j, each one: 0=0, 1=s, 2=v, 4=d, 5=d, 8=f)

      ammn(0)=0d0
      ammn(1)=0d0
      ammn(2)=dble(ammsqq)+amd
      ammn(3)=dble(ammsqq)
      ammn(4)=dble(ammsqq)
      ammn(5)=dble(ammsqd)+amd
      ammn(6)=dble(ammsqd)+amd
      ammn(7)=0d0
      ammn(8)=dble(ammsdd)+amd
      ammn(9)=dble(ammsqd)+amd
      ammn(10)=dble(ammsqd)+amd
      ammn(12)=dble(ammsqd)+amd
      ammn(16)=0.14d0

c minimal pomeron masses (0=soft or gg, 1=qg, 2=gq, 3=qq)

      amprmn(0)=ammsqq
      !amprmn(1/2/3) defined later

c cutoff for virtual pomeron (0=0, 1=soft Pom, 2=regge, 3=hard)

      xvpr(0)=0d0
      xvpr(1)=min(1d0,dble(cumpox**2)/s)
      xvpr(2)=min(1d0,dble(ammsqq**2)/s)
      xvpr(3)=min(1d0,(2d0*ammn(2))**2/s)

c minimal remnant masses (0=meson, 1=baryon)

      idpj=idpji
      xdm=0.35d0                  !<pt>
      call idmass(idpj,ampj)
      if(iabs(idpj).gt.1000)then
       ampmn(0)=0.14d0+xdm
       ampmn(1)=dble(ampj)+xdm
      else
       ampmn(0)=dble(ampj)+xdm
       ampmn(1)=0.94d0+xdm
      endif
      idtg=idtgi
      if(idtg.eq.0)idtg=1120
      call idmass(idtg,amtg)
      if(iabs(idtg).gt.1000)then
       amtmn(0)=0.14d0+xdm
       amtmn(1)=dble(amtg)+xdm
      else
       amtmn(0)=dble(amtg)+xdm
       amtmn(1)=0.94d0+xdm
      endif

c minimal excitation masses (0=meson, 1=baryon
c                            0=no excit, 1=nondiffr, 2=diffr,
c                                   6=nondiffr but no pomeron)

      xdm2=0.15d0
      amemn(0,0)=0.d0
      amemn(1,0)=0.d0
      amemn(0,4)=0.d0
      amemn(1,4)=0.d0
      amemn(0,6)=0.d0
      amemn(1,6)=0.d0

      amemn(0,1)=xdm2!+dble(delrex)
      amemn(0,2)=xdm2+dble(delrex)
      amemn(0,3)=xdm2+dble(delrex)
      amemn(0,5)=xdm2+dble(delrex+amdrmin) !remnant excited with MPI

      amemn(1,1)=xdm2!+dble(delrex)
      amemn(1,2)=xdm2+dble(delrex)
      amemn(1,3)=xdm2+dble(delrex)
      amemn(1,5)=xdm2+dble(delrex+amdrmin) !remnant excited with MPI

c maximal excitation masses (0=no excit, 1=nondiffr, 2=diffr)

      amemx(0)=2d0*xdm
      amemx(1)=plc
      amemx(2)=plc

      if(idpj.gt.1000)then     ! baryon

c initial quark configuration
       ivp0=3
       iap0=0
       idp0=1
       isp0=1

c no val quark for exotic projectile
       if(iremn.ge.2.and.(idpj.ne.1120.and.idpj.ne.1220))ivp0=0

      elseif(idpj.lt.-1000)then     ! antibaryon

c initial quark configuration
       ivp0=0
       iap0=3
       idp0=1
       isp0=1

c no val quark for exotic projectile
       if(iremn.ge.2.and.(idpj.ne.-1120.and.idpj.ne.-1220))iap0=0

      else      ! meson

c initial quark configuration
       ivp0=1
       iap0=1
       idp0=0
       isp0=0

c no val quark for exotic projectile
       if(iremn.ge.2.and.(mod(abs(idpj/100),10).gt.4
     &                 .or.mod(abs(idpj/10),10).gt.4
     &    .or.mod(abs(idpj/100),10)/mod(abs(idpj/10),10).eq.1))then
         ivp0=0
         iap0=0
       endif
      endif

      if(idtg.gt.1000)then    ! baryon

c initial quark configuration
       ivt0=3
       iat0=0
       idt0=1
       ist0=1

c no val quark for exotic target
       if(iremn.ge.2.and.(idtg.ne.1120.and.idtg.ne.1220))ivt0=0

      elseif(idtg.lt.-1000)then   ! antibaryon

c initial quark configuration
       ivt0=0
       iat0=3
       idt0=1
       ist0=1

c no val quark for exotic target
       if(iremn.ge.2.and.(idtg.ne.-1120.and.idtg.ne.-1220))iat0=0

      else       ! meson

c initial quark configuration
       ivt0=1
       iat0=1
       idt0=0
       ist0=0

c no val quark for exotic target
       if(iremn.ge.2.and.(mod(abs(idtg/100),10).gt.4
     &                 .or.mod(abs(idtg/10),10).gt.4
     &    .or.mod(abs(idtg/100),10)/mod(abs(idtg/10),10).eq.1))then
         ivt0=0
         iat0=0
       endif

      endif

c for parton position density

c       hsigp=1d0/dble(4.*0.0389*r2har(iclpro)) !0.0389 for conversion from GeV-2 to fm2 : (hbar*c)**2=0.0389 GeV^2.fm^2, and r2had in GeV^-2
c       hnormp=dble(sqrt(0.0389*hsigp/pi))     !(in Omega with use GeV not fm)
c       hsigt=1d0/dble(4.*0.0389*r2har(icltar))
c       hnormt=dble(sqrt(0.0389*hsigt/pi))
c       ssigp=1d0/dble(4.*0.0389*(max(gwidth*rexres(iclpro),1.)
c     &          +0.5*max(slopom,slopoms)*log(s)))
c       snormp=dble(sqrt(0.0389*ssigp/pi))
c       ssigt=1d0/dble(4.*0.0389*(max(gwidth*rexres(icltar),1.)
c     &          +0.5*max(slopom,slopoms)*log(s)))
c       snormt=dble(sqrt(0.0389*ssigt/pi))

       if(iset.eq.0)then
c counters

       antot=0.
       ansh=0.
       ansf=0.
       antotf=0.
       anshf=0.
       ansff=0.
       pp4max=0.
       pp4ini=0.
       andropl=0.
       anstrg0=0.
       anstrg1=0.
       anreso0=0.
       anreso1=0.
       anghadr=0.
       antotre=0.
       anintdiff=0.
       anintsdif=0.
       anintddif=0.
       anintine=0.

c Fit parameters

       if(isetcs.gt.-2)then
         q2pmin(1)=q2nmin
         q2pmin(2)=q2nmin
         call ipoOm5Tables(1)   !used for psvin and thus for om51p
         call paramini(1)
         if(ish.ge.4)then
           write(ifch,'(a,2f6.2,a)')
     *          '****** q2pmin =',q2pmin(1),q2pmin(2),' **********'
           do i=idxD0,idxD1
             write(ifch,'(9(a,f8.4))')
     *   'AlpD:',alpD(i,iclpro,icltar)
     * ,' AlpDp:',alpDp(i,iclpro,icltar)
     * ,' AlpDpp:',alpDpp(i,iclpro,icltar)
     * ,' BetD:',betD(i,iclpro,icltar)
     * ,' BetDp:',betDp(i,iclpro,icltar)
     * ,' BetDpp:',betDpp(i,iclpro,icltar)
     * ,' GamD:',gamD(i,iclpro,icltar)
     * ,' DelD:',delD(i,iclpro,icltar)
           enddo
         endif
c         bkmxndif=conbmxndif()
         bkmx=conbmx()
c         if(ish.ge.3)write(ifch,*)'bkmx,bkmxndif',bkmx,bkmxndif
         if(ish.ge.3)write(ifch,*)'bkmx',bkmx
         if(iokoll.eq.0.and.(maproj.gt.1.or.matarg.gt.1))then
           bmax=rmproj+rmtarg
         else
           bmax=bkmx
         endif
         call xsigma
         !update xmxrem only for nuclei to keep constiency between xs and MC
         if(iokoll.eq.0)then
           if(maproj.gt.1)then
c maximum momentum fraction for mass (m*<<sqrt(s/(2.m.R/hbar/c)) limit due to coherence lenght in diffraction (Feinberg and Pomeranchuk))
             xmxrem(1)=min(1.,  !max(1.01*amhadr(icltar)/e,
     .            xmxmas/(10.135*amproj*rmproj)) !projectile (10.135(GeV-1.fm-1)=2./hbar/c=2/0.197(GeV.fm))
           endif
           if(matarg.gt.1)then
             xmxrem(2)=min(1.,  !max(1.01*amhadr(icltar)/e,
     .            xmxmas/(10.135*amtarg*rmtarg)) !target (10.135(GeV-1.fm-1)=2./hbar/c=2/0.197(GeV.fm))
           endif
c           print *,xmxrem
         elseif(iokoll.lt.0)then !factorized hard cross section
           sigine=SigGFF()/float(abs(iokoll))
           write(ifmt,*)'Cross-section for iokoll=',iokoll,sigine
         endif

       else
         xmxrem(1)=1.
         xmxrem(2)=1.
         xnpomcol=1.
         sigine=xsectionpar()
       endif
c        sigine=xsectionpar()
      endif

      return
      end

c-----------------------------------------------------------------------
      subroutine emsigr
c-----------------------------------------------------------------------
c     initialize grid
c     nprmx, npr and nprt are intialized in GfunparK (omg.f)
c-----------------------------------------------------------------------

#include "aaa.h"
#include "ems.h"

      common/emsptl/nppr(npommx,kollmx),npproj(mamx),nptarg(mamx)
      common/crangrid/irangrid(npommx,kollmx)

      call utpri('emsigr',ish,ishini,5)

      do k=1,koll  !----k-loop---->

c initial value for interaction type

       itpr(k)=0

c initial value for nuclear splitting

c       do ir=1,2
c         knucnt(ir,k)=0
c         do ncon=1,mamx
c           npnuc(ncon,ir,k)=0
c           irnuc(ncon,ir,k)=0
c           xxnuc(ncon,ir,k)=0d0
c         enddo
c       enddo

c initialize grid

       btest=0.
       do n=1,nprmx(k)
        call conigr(n,k)
        btest=btest+bk(k)       !bhpr(n,k)
        idpr(n,k)=0
        idfpr(n,k)=0
        ivpr(n,k)=1
        nppr(n,k)=0
        nbkpr(n,k)=0
        nvpr(n,k)=0
        idsppr(n,k)=0
        idstpr(n,k)=0
        idrpr(n,k)=0
        idhpr(n,k)=0
        xpr(n,k)=0d0
        ypr(n,k)=0d0
        xppr(n,k)=0d0
        xmpr(n,k)=0d0
        xp1pr(n,k)=0d0
        xp2pr(n,k)=0d0
        xm1pr(n,k)=0d0
        xm2pr(n,k)=0d0
        xp1pr(n,k)=0d0
        xp2pr(n,k)=0d0
        xm1pr(n,k)=0d0
        xm2pr(n,k)=0d0
        idp1pr(n,k)=0
        idp2pr(n,k)=0
        idm1pr(n,k)=0
        idm2pr(n,k)=0
        xxp1pr(n,k)=0d0
        xyp1pr(n,k)=0d0
        xxp2pr(n,k)=0d0
        xyp2pr(n,k)=0d0
        xxm1pr(n,k)=0d0
        xym1pr(n,k)=0d0
        xxm2pr(n,k)=0d0
        xym2pr(n,k)=0d0
        do jij=1,mxsplit
        call xpprbor_set(jij,n,k, 0. )
        call xmprbor_set(jij,n,k, 0. )
        call ptprboo_set(1,jij,n,k, 0. )
        call ptprboo_set(2,jij,n,k, 0. )
        call rapprboo_set(1,jij,n,k, 0. )
        call rapprboo_set(2,jij,n,k, 0. )
        call gbpom_set(1,jij,n,k, 0. )
        call gbpom_set(2,jij,n,k, 0. )
        enddo
        q2kmin(1,n,k)=0.
        q2kmin(2,n,k)=0.
        nemispr(1,n,k)=0
        nemispr(2,n,k)=0
c        npp(k)=0
c        npt(k)=0
        irangrid(n,k)=1
        if(rangen().lt.0.5)then
        irangrid(n,k)=-1
        endif
       enddo
       if(ish.ge.5)then
         if(nprmx(k).gt.0)btest=btest/float(nprmx(k))
         write(ifch,'(a,i6,a,i2,a,i2,a,f5.2,a,f5.2)')'Pair ',k,' ('
     &         ,iproj(k),',',itarg(k),') : bk=',bk(k),' <b>=',btest
       endif

      enddo !  <----k-loop-----

      call utprix('emsigr',ish,ishini,5)
      return
      end

c-----------------------------------------------------------------------
      subroutine emsipt
c-----------------------------------------------------------------------
c initialize projectile and target
c-----------------------------------------------------------------------

#include "aaa.h"
#include "ems.h"

      common/emsptl/nppr(npommx,kollmx),npproj(mamx),nptarg(mamx)
      common/cems5/plc,s
      common/ems6/ivp0,iap0,idp0,isp0,ivt0,iat0,idt0,ist0
      common /cncl/xproj(mamx),yproj(mamx),zproj(mamx)
     *            ,xtarg(mamx),ytarg(mamx),ztarg(mamx)

      double precision s,plc

      call utpri('emsipt',ish,ishini,7)
c initialize projectile

      do i=1,maproj
       idp(i)=idp0
       ivp(i)=ivp0+iap0
       iap(i)=iap0
       isp(i)=isp0
       iep(i)=-1
       ifp(i)=0
       kolp(i)=0
c       npp(i)=0
       npproj(i)=0
       xxp(i)=0d0
       xyp(i)=0d0
       xpmn(i)=dble(xminremn)*(amemn(idp(i),0)+ampmn(isp(i)))**2/s  !minimum xp=xpp*xmp given by the rest mass
       xpmx(i)=min(0.25d0,amemx(1)**2/s)  !maximum xp=xpp*xmp given by the maximum mass  (at most particle beam energy sqrt(s)/2)
       xpos(i)=0.9d0*(amemx(0)+ampmn(isp(i)))**2/s  !refined in ProCop befor being used
       xppmx(i)=0.5d0*(1d0+sqrt((1d0+0.5d0*sqrt(xpmn(i)))
     .                         *(1d0-0.5d0*sqrt(xpmn(i))))) !xppmx=0.5*(1+sqrt(1-0.25*M**2/s))
       xmpmn(i)=0.5d0*(1d0-sqrt((1d0+0.5d0*sqrt(xpmn(i)))
     .                         *(1d0-0.5d0*sqrt(xpmn(i))))) !xmpmn=0.5*(1-sqrt(1-0.25*M**2/s))
       xmpmx(i)=xpmx(i)/xppmx(i)      !maximum mass with maximum xp momentum
       xppmn(i)=xpmn(i)/xmpmx(i)      !used in metropolis
       xpp(i)=xppmx(i)                !take into the mass at the beginning
       xmp(i)=xmpmn(i)
       xppst(i)=0.d0
       xmpst(i)=0.d0
       xposst(i)=0.d0
      enddo

c initialize target

      do j=1,matarg
       idt(j)=idt0
       ivt(j)=ivt0+iat0
       iat(j)=iat0
       ist(j)=ist0
       iet(j)=-1
       ift(j)=0
       kolt(j)=0
c       npt(j)=0
       nptarg(j)=0
       xxt(j)=0d0
       xyt(j)=0d0
       xtmn(j)=dble(xminremn)*(amemn(idt(j),0)+amtmn(ist(j)))**2/s  !minimum xt=xpt*xmt given by the rest mass
       xtmx(j)=min(0.25d0,amemx(1)**2/s)  !maximum xt=xpt*xmt given by the maximum mass  (at most particle beam energy sqrt(s)/2)
       xtos(j)=0.9d0*(amemx(0)+amtmn(ist(j)))**2/s  !refined in ProCot befor being used
       xmtmx(j)=0.5d0*(1d0+sqrt((1d0+0.5d0*sqrt(xtmn(j)))
     .                         *(1d0-0.5d0*sqrt(xtmn(j))))) !xppmx=0.5*(1+sqrt(1-0.25*M**2/s))
       xptmn(j)=0.5d0*(1d0-sqrt((1d0+0.5d0*sqrt(xtmn(j)))
     .                         *(1d0-0.5d0*sqrt(xtmn(j))))) !xmpmn=0.5*(1-sqrt(1-0.25*M**2/s))
       xptmx(j)=xtmx(j)/xmtmx(j)      !maximum mass with maximum xt momentum
       xmtmn(j)=xtmn(j)/xptmx(j)      !used in metropolis
       xmt(j)=xmtmx(j)                !take into the mass at the beginning
       xpt(j)=xptmn(j)
       xmtst(j)=0.d0
       xptst(j)=0.d0
       xtosst(j)=0.d0
      enddo

      if(ish.ge.7)then
      do i=1,maproj
        write(ifch,115)i,xpmn(i),xpmx(i),xpos(i)
     .                  ,xppmx(i),xmpmn(i),xmpmx(i),xppmn(i)
      enddo
      do j=1,matarg
        write(ifch,116)j,xtmn(j),xtmx(j),xtos(j)
     .                  ,xmtmx(j),xptmn(j),xptmx(j),xmtmn(j)
      enddo
  115 format(1x,'/proj/',1i6,1p,7(e10.3,1x))
  116 format(1x,'/targ/',1i6,1p,7(e10.3,1x))
      endif

      return
      end


c-----------------------------------------------------------------------
      subroutine emszz
c-----------------------------------------------------------------------
c     completes /cptl/ for nucleons, checks for no interaction
c     writes   /cevt/
c-----------------------------------------------------------------------
#include "aaa.h"
#include "ems.h"
      common/nucl3/phi,bimp
      common/col3/ncol,kolpt,ncoli
      integer kolpz(mamx),koltz(mamx)

      call utpri('emszz ',ish,ishini,6)

c     write /cptl/
c     ------------

c      if(iokoll.eq.1)then   ! precisely matarg collisions

c nothing to do

c      else

c determine ncol

       ncolx=ncol
       ncolix=ncoli
       ncol=0
       ncoli=0
       do 8 k=1,koll
      if(ish.ge.10)write(ifch,*)'k,itpr,ncol,ncolx',k,itpr(k),ncol,ncolx
        if(itpr(k).eq.0)goto 8
        ncol=ncol+1
        i=iproj(k)
        j=itarg(k)
        if(abs(itpr(k)).eq.1)then
          ncoli=ncoli+1
          if(ish.ge.5)
     .    write(ifch,'(a,i7,3x,2i6)')' ++++ems-111++++',ncoli,i,maproj+j
        else
          if(ish.ge.5)
     .    write(ifch,'(a,i7,3x,2i6)')' ++++ems-222++++',ncol,i,maproj+j
        endif
        istptl(i)=1
        iorptl(i)=-1
        tivptl(2,i)=coord(4,k)
        istptl(maproj+j)=1
        iorptl(maproj+j)=-1
        tivptl(2,maproj+j)=coord(4,k)
8      continue
       if(ncolx.ne.ncol)write(6,*)'ncolx,ncol:', ncolx,ncol
       if(ncolx.ne.ncol)call utstop('********ncolx.ne.ncol********&')
c       if(ncolix.ne.ncoli)call utstop('*******ncolix.ne.ncoli*******&')
       if(ncol.eq.0)goto 1001

c determine npj, ntg

       do ip=1,maproj
        kolpz(ip)=0
       enddo
       do it=1,matarg
        koltz(it)=0
       enddo
      do k=1,koll
       if(itpr(k).ne.0.and.abs(itpr(k)).ne.3)then
        ip=iproj(k)
        it=itarg(k)
        kolpz(ip)=kolpz(ip)+1
        koltz(it)=koltz(it)+1
       endif
      enddo
      npj=0
      do ip=1,maproj
       if(kolpz(ip).gt.0.or.iep(ip).ge.30)npj=npj+1
      enddo
      ntg=0
      do it=1,matarg
       if(koltz(it).gt.0.or.iet(it).ge.30)ntg=ntg+1
      enddo
c     write(6,*)'npj,ntg,npj+ntg:',npj,ntg,npj+ntg

c       endif

c     write /cevt/
c     ------------

      nevt=1
      bimevt=bimp
      phievt=phi
      kolevt=ncol
      koievt=ncoli
      npjevt=npj
      ntgevt=ntg
      pmxevt=pnll
      egyevt=engy
c      print*,' ===== ',kolevt,koievt,nglevt,' ====='

c     exit
c     ----

      if(ish.ge.7)then
      do n=1,nptl
      write(ifch,115)iorptl(n),jorptl(n),n,istptl(n)
     *,tivptl(1,n),tivptl(2,n)
      enddo
  115 format(1x,'/cptl/',2i6,2i10,2(e10.3,1x))
      endif

1000  continue
      call utprix('emszz ',ish,ishini,6)
      return

1001  continue
      if(ish.ge.3)then
      write(ifch,*)
      write(ifch,*)'   ***** no interaction!!!'
      write(ifch,*)'   ***** ncol=0 detected in emszz'
      write(ifch,*)
      endif
      goto 1000

      end

c-----------------------------------------------------------------------
      subroutine ProCop(i,ii)
c-----------------------------------------------------------------------
c Propose Coordinates of remnants from active projectile nucleons
c-----------------------------------------------------------------------

#include "aaa.h"
#include "ems.h"
#include "sem.h"

      double precision xmptmp,aproj
      common/cems5/plc,s
      common/emsptl/nppr(npommx,kollmx),npproj(mamx),nptarg(mamx)
      integer icrmn(2),jc(nflav,2),icini(2)
      double precision s,plc

      nptl=nptl+1
      npproj(i)=nptl
      idptl(nptl)=idptl(ii)*100+99  !100*10**idp(i)+iep(i)
      istptl(nptl)=40
      ityptl(nptl)=40
      iorptl(nptl)=ii
      jorptl(nptl)=0
      ifrptl(1,nptl)=0
      ifrptl(2,nptl)=0
      do j=1,2
        do k=1,nflav
          jc(k,j)=0
        enddo
      enddo

      istptl(ii)=1

c     determine kolz

      if(lproj(i).gt.1)then
        zmax=-ainfin
        kolz=0
        do l=1,lproj(i)
          k=kproj(i,l)
          z=coord(3,k)
          if(itpr(k).ne.0.and.z.gt.zmax)then
            zmax=z
            kolz=k
          endif
        enddo
      else
        kolz=1
      endif
c      if(kolz.eq.0)call utstop(' kolz=0 (proj)&')
      if(kolz.eq.0)then
        t=0.
      else
        t=coord(4,kolz)
      endif

      xorptl(1,nptl)=xorptl(1,ii)
      xorptl(2,nptl)=xorptl(2,ii)
      xorptl(3,nptl)=xorptl(3,ii)
      xorptl(4,nptl)=t
      tivptl(1,nptl)=t
      tivptl(2,nptl)=t
c save q2s min and max attached to this remnant
      qsqptl(nptl)=0.
      zpaptl(1,nptl)=0.
      zpaptl(2,nptl)=0

      nqu=0
      naq=0
      if(iremn.ge.2)then   !update icproj
        idp(i)=abs(idp(i))
        k=1
        do n=1,nrflav
          jc(n,k)=jcpref(n,k,i)
          nqu=nqu+jc(n,k)
        enddo
        k=2
        do n=1,nrflav
          jc(n,k)=jcpref(n,k,i)
          naq=naq+jc(n,k)
        enddo
        isum=nqu+naq
        call idenco(jc,icrmn,iret)
        if(iret.eq.0.and.(isum.le.3.or.iremn.ne.3))then
          icproj(1,i)=icrmn(1)
          icproj(2,i)=icrmn(2)
        elseif(iremn.eq.3)then
      write(ifch,*)'Problem in projectile flavor :',i,' ->',jc,' :',isum
          call utstop('Procop: Problem in projectile flavor !&')
        else     !for iremn=2 and large number of quark define icproj=999999
          icproj(1,i)=999999
          icproj(2,i)=999999
        endif
      endif

      icrmn(1)=icproj(1,i)
      icrmn(2)=icproj(2,i)

      if(iremn.ge.1)then      !excited remnant ?
        call idtr4(idptl(ii),icini)
        if(ish.ge.5)write(ifch,*)'Procop icini proj',i,icini,' ->',icrmn
        if((icrmn(1)-icini(1))+(icrmn(2)-icini(2)).ne.0)then
          if(iep(i).eq.60)then
            write(ifch,'(a,d25.15)')
     &'Flavor problem in proj for pseudo-inelastic collision !',seedc
c          elseif(iep(i).eq.40)then !can not be low mass excitation
c            iep(i)=11
          endif
c        elseif(iep(i).eq.10)then !non excited inelastic remnant
c          iep(i)=0
        endif
      endif

c      if(iremn.eq.2)then
c          if(.not.((nqu.eq.3.and.naq.eq.0).or.(nqu.eq.0.and.naq.eq.3)
c     &               .or.(nqu.eq.1.and.naq.eq.1)).and.iet(j)/10.eq.1)
c     &  iep(i)=30+mod(iep(i),10)
c      endif

      if(ish.ge.5)write(ifch,'(a,i3,a,i3,a,i2)')
     &            'Procop part ',ii,', iep(',i,'): ',iep(i)

      if(iremn.le.1)call iddeco(icrmn,jc)
      if(iep(i).ge.1.and.iep(i).ne.60)then
        aproj=dble(max(amproj,fremnux(jc)))
      else
        aproj=dble(max(amproj,fremnux2(jc)))
      endif
c      aprojex=max(ampmn(isp(i))+amemn(idp(i),iep(i)/10)
c     &           ,dble(fremnux(jc)))
      xmptmp=(aproj**2+xxp(i)*xxp(i)+xyp(i)*xyp(i))
     &       /(xpp(i)*s)
      xpos(i)=xpp(i)*xmptmp
      if(ish.ge.5)write(ifch,*)'Procop mass : ',aproj,xpos(i)*s
      if(xmptmp.gt.1.d0)then
        xmptmp=0.d0
      if(ish.ge.1)write(ifmt,*)'Warning in ProCop, Remnant mass too low'
      endif

      pptl(1,nptl)=sngl(xxp(i))
      pptl(2,nptl)=sngl(xyp(i))
      pptl(3,nptl)=sngl((xpp(i)-xmptmp)*plc/2d0)
      pptl(4,nptl)=sngl((xpp(i)+xmptmp)*plc/2d0)
      pptl(5,nptl)=aproj


c      write(ifmt,*)'ProCop',i,nptl

      return

      end

c-----------------------------------------------------------------------
      subroutine ProCot(j,jj)
c-----------------------------------------------------------------------
c Propose Coordinates of remnants from active targets nucleons
c-----------------------------------------------------------------------

#include "aaa.h"
#include "ems.h"
#include "sem.h"

      double precision xpttmp,atarg
      common/cems5/plc,s
      common/emsptl/nppr(npommx,kollmx),npproj(mamx),nptarg(mamx)
      integer icrmn(2),jc(nflav,2),icini(2)
      double precision s,plc

      nptl=nptl+1
      nptarg(j)=nptl

      idptl(nptl)=idptl(jj)*100+99    !100*10**idt(j)+iet(j)
      istptl(nptl)=40
      ityptl(nptl)=50
      iorptl(nptl)=jj
      jorptl(nptl)=0
      ifrptl(1,nptl)=0
      ifrptl(2,nptl)=0
      do k=1,2
        do i=1,nflav
          jc(i,k)=0
        enddo
      enddo

      istptl(jj)=1

c     determine kolz

      if(ltarg(j).gt.1)then
        zmin=ainfin
        kolz=0
        do l=1,ltarg(j)
          k=ktarg(j,l)
          z=coord(3,k)
          if(itpr(k).ne.0.and.z.lt.zmin)then
            zmin=z
            kolz=k
          endif
        enddo
      else
        kolz=1
      endif
c      if(kolz.eq.0)call utstop(' kolz=0 (targ)&')
      if(kolz.eq.0)then
        t=0.
      else
        t=coord(4,kolz)
      endif

      xorptl(1,nptl)=xorptl(1,jj)
      xorptl(2,nptl)=xorptl(2,jj)
      xorptl(3,nptl)=xorptl(3,jj)
      xorptl(4,nptl)=t
      tivptl(1,nptl)=t
      tivptl(2,nptl)=t
c save q2s min and max attached to this remnant
      qsqptl(nptl)=0.
      zpaptl(1,nptl)=0.
      zpaptl(2,nptl)=0.

      nqu=0
      naq=0
      if(iremn.ge.2)then   !update ictarg
        idt(j)=abs(idt(j))
        k=1
        do n=1,nrflav
          jc(n,k)=jctref(n,k,j)
          nqu=nqu+jc(n,k)
        enddo
        k=2
        do n=1,nrflav
          jc(n,k)=jctref(n,k,j)
          naq=naq+jc(n,k)
        enddo
        isum=nqu+naq
        call idenco(jc,icrmn,iret)
        if(iret.eq.0.and.(isum.le.3.or.iremn.ne.3))then
          ictarg(1,j)=icrmn(1)
          ictarg(2,j)=icrmn(2)
        elseif(iremn.eq.3)then
      write(ifch,*)'Problem in projectile flavor :',j,' ->',jc,' :',isum
          call utstop('Procot: Problem in target flavor !&')
        else     !for iremn=2 and large number of quark define ictarg=999999
          ictarg(1,j)=999999
          ictarg(2,j)=999999
        endif
      endif

      icrmn(1)=ictarg(1,j)
      icrmn(2)=ictarg(2,j)

      if(iremn.ge.1)then      !excited remnant ?
        call idtr4(idptl(jj),icini)
        if(ish.ge.5)write(ifch,*)'Procot icini targ',j,icini,' ->',icrmn
        if((icrmn(1)-icini(1))+(icrmn(2)-icini(2)).ne.0)then
          if(iet(j).eq.60)then
            write(ifch,'(a,d25.15)')
     &'Flavor problem in targ for pseudo-inelastic collision !',seedc
c          elseif(iet(j).eq.40)then !can not be low mass excitation
c            iet(j)=11
          endif
c        elseif(iet(j).eq.10)then   !non excited inelastic remnant
c          iet(j)=0
        endif
      endif

c      if(iremn.eq.2)then
c          if(.not.((nqu.eq.3.and.naq.eq.0).or.(nqu.eq.0.and.naq.eq.3)
c     &               .or.(nqu.eq.1.and.naq.eq.1)).and.iet(j)/10.eq.1)
c     &  iet(j)=30+mod(iet(j),10)
c      endif
      if(ish.ge.5)write(ifch,'(a,i3,a,i3,a,i2)')
     &            'Procot part ',jj,', iet(',j,'): ',iet(j)



      if(iremn.le.1)call iddeco(icrmn,jc)
      if(iet(j).ge.1.and.iet(j).ne.60)then
        atarg=dble(max(amtarg,fremnux(jc)))
      else
        atarg=dble(max(amtarg,fremnux2(jc)))
      endif
c      atargex=max(amtmn(ist(j))+amemn(idt(j),iet(j)/10)
c     &           ,dble(fremnux(jc)))
      xpttmp=(atarg**2+xxt(j)*xxt(j)+xyt(j)*xyt(j))
     &       /(xmt(j)*s)
      xtos(j)=xpttmp*xmt(j)
      if(ish.ge.5)write(ifch,*)'Procot mass : ',atarg,xtos(j)*s
      if(xpttmp.gt.1.d0)then
        xpttmp=0.d0
      if(ish.ge.1)write(ifch,*)'Warning in ProCot, Remnant mass too low'
      endif

      pptl(1,nptl)=sngl(xxt(j))
      pptl(2,nptl)=sngl(xyt(j))
      pptl(3,nptl)=sngl((xpttmp-xmt(j))*plc/2d0)
      pptl(4,nptl)=sngl((xpttmp+xmt(j))*plc/2d0)
      pptl(5,nptl)=atarg

c      write(ifmt,*)'ProCot',j,nptl

      return
      end

c-----------------------------------------------------------------------
      subroutine emswrp(i,ii)
c-----------------------------------------------------------------------
c Update values from Procop
c-----------------------------------------------------------------------

#include "aaa.h"
#include "ems.h"

      double precision p5sq
      common/cems5/plc,s
      common/emsptl/nppr(npommx,kollmx),npproj(mamx),nptarg(mamx)
      double precision s,plc
      parameter(eps=1.e-5)

      if(npproj(i).eq.0)then
        write(*,*)'emswrp i ii',i,ii
        call utstop('emswrp with npproj=0 should never happen !&')

c        t=xorptl(4,kolp(i))
c        istptl(ii)=1
c        iorptl(ii)=-1
c        tivptl(2,ii)=t
c        nptl=nptl+1
c        npproj(i)=nptl
c        idptl(nptl)=idptl(ii)*100+99 !100*10**idp(i)+iep(i)
c        istptl(nptl)=40
c        ityptl(nptl)=40
c        iorptl(nptl)=ii
c        jorptl(nptl)=kolp(i)
c        ifrptl(1,nptl)=0
c        ifrptl(2,nptl)=0
c        xorptl(1,nptl)=xorptl(1,ii)
c        xorptl(2,nptl)=xorptl(2,ii)
c        xorptl(3,nptl)=xorptl(3,ii)
c        xorptl(4,nptl)=t
c        tivptl(1,nptl)=t
c        tivptl(2,nptl)=t
c        mm=nptl
c        kolp(i)=1
      else
        mm=npproj(i)
      endif
      pptl(1,mm)=sngl(xxp(i))
      pptl(2,mm)=sngl(xyp(i))
      pptl(3,mm)=sngl((xpp(i)-xmp(i))*plc/2d0)
      pptl(4,mm)=sngl((xpp(i)+xmp(i))*plc/2d0)
      if(pptl(4,mm).lt.-eps)call utstop('E pro<0 !&')
      p5sq=xpp(i)*xmp(i)*s-xxp(i)*xxp(i)-xyp(i)*xyp(i)
      if(p5sq.gt.1.d-10)then
        pptl(5,mm)=sngl(sqrt(p5sq))
      elseif(iep(i).eq.0)then
        pptl(5,mm)=pptl(5,ii)
      else
        if(ish.ge.2)then
          write(ifch,*)'problem with mass for projectile, '
     &         ,'continue with zero mass'
          write(ifch,*)i,mm,xxp(i),xyp(i),xpp(i),xmp(i),p5sq
        endif
        pptl(5,mm)=0.
      endif

      do l=1,4
       ibptl(l,mm)=0
      enddo

      return

      end

c-----------------------------------------------------------------------
      subroutine emswrt(j,jj)
c-----------------------------------------------------------------------
c Update values from Procot
c-----------------------------------------------------------------------

#include "aaa.h"
#include "ems.h"

      double precision p5sq
      common/cems5/plc,s
      common/emsptl/nppr(npommx,kollmx),npproj(mamx),nptarg(mamx)
      double precision s,plc
      parameter(eps=1.e-5)

      if(nptarg(j).eq.0)then

        write(*,*)'emswrt j jj',j,jj
        call utstop('emswrt with nptarg=0 should never happen !&')

      else
        mm=nptarg(j)
      endif
      pptl(1,mm)=sngl(xxt(j))
      pptl(2,mm)=sngl(xyt(j))
      pptl(3,mm)=sngl((xpt(j)-xmt(j))*plc/2d0)
      pptl(4,mm)=sngl((xpt(j)+xmt(j))*plc/2d0)
      if(pptl(4,mm).lt.-eps)call utstop('E targ<0 !&')
      p5sq=xpt(j)*xmt(j)*s-xxt(j)*xxt(j)-xyt(j)*xyt(j)
      if(p5sq.gt.1.d-10)then
        pptl(5,mm)=sngl(sqrt(p5sq))
      elseif(iet(j).eq.0)then
        pptl(5,mm)=pptl(5,jj)
      else
        if(ish.ge.2)then
          write(ifch,*)'problem with mass for target, '
     &            ,'continue with zero mass'
          write(ifch,*)j,mm,xxt(j),xyt(j),xpt(j),xmt(j),p5sq
        endif
        pptl(5,mm)=0.
      endif

      do l=1,4
       ibptl(l,mm)=0
      enddo

      return
      end

c-----------------------------------------------------------------------
      subroutine emswrpom(k,i,j)
c-----------------------------------------------------------------------
c should be call for diffractive pair too because of CD with Reggeon
c-----------------------------------------------------------------------

#include "aaa.h"
#include "ems.h"

      common/cems5/plc,s
      common/emsptl/nppr(npommx,kollmx),npproj(mamx),nptarg(mamx)
      double precision s,px,py,plc


      do 30 n=1,nprmx(k)
       if(idpr(n,k).eq.0.or.ivpr(n,k).eq.0)goto 30
       if(idpr(n,k).eq.2)goto 30 !reggeon with remnant mass only for the moment ???????????????
       nptl=nptl+1
       nppr(n,k)=nptl
       px=xxp1pr(n,k)+xxp2pr(n,k)+xxm1pr(n,k)+xxm2pr(n,k)
       py=xyp1pr(n,k)+xyp2pr(n,k)+xym1pr(n,k)+xym2pr(n,k)
       pptl(1,nptl)=sngl(px)
       pptl(2,nptl)=sngl(py)
       pptl(3,nptl)=sngl(dsqrt(xpr(n,k))*dsinh(ypr(n,k))*plc)
       pptl(4,nptl)=sngl(dsqrt(xpr(n,k))*dcosh(ypr(n,k))*plc)
       pptl(5,nptl)=sngl(dsqrt(xpr(n,k)*s-px*px-py*py))
   !    print*,pptl(5,nptl)/plc
       idptl(nptl)=idpr(n,k)*10000
     &     +idp1pr(n,k)*1000
     &     +idp2pr(n,k)*100
     &     +idm1pr(n,k)*10
     &     +idm2pr(n,k)
       idptl(nptl)=idptl(nptl)*100+99
       if(ivpr(n,k).eq.2)idptl(nptl)=-idptl(nptl)
       istptl(nptl)=30
       iorptl(nptl)=i
       jorptl(nptl)=j
       ifrptl(1,nptl)=0
       ifrptl(2,nptl)=0
       xorptl(1,nptl)=coordpr(1,n,k)
       xorptl(2,nptl)=coordpr(2,n,k)
       xorptl(3,nptl)=coord(3,k)
       xorptl(4,nptl)=coord(4,k)
       tivptl(1,nptl)=coord(4,k)
       tivptl(2,nptl)=coord(4,k)
       qsqptl(nptl)=max(q2kmin(1,n,k),q2kmin(2,n,k))
       zpaptl(1,nptl)=q2kmin(1,n,k)
       zpaptl(2,nptl)=q2kmin(2,n,k)
       if(idpr(n,k).eq.1)then
        ityptl(nptl)=20
        if(itpr(k).gt.0)ityptl(nptl)=25
       elseif(idpr(n,k).eq.6)then
        ityptl(nptl)=21
        if(itpr(k).gt.0)ityptl(nptl)=26
       elseif(idpr(n,k).eq.2)then
        ityptl(nptl)=22
        if(itpr(k).gt.0)ityptl(nptl)=27
       elseif(idpr(n,k).eq.3)then
        ityptl(nptl)=30+idhpr(n,k)+1
        if(itpr(k).gt.0)ityptl(nptl)=35
       else
        call utstop('emswrpom: unknown id&')
       endif
       do l = 1,4
        ibptl(l,nptl)=0
       enddo
30    continue

      return
      end

c-----------------------------------------------------------------------
      subroutine xxEmsPx(iii,jot,jij,xmc,ymc,xpl,xml,npos,pt1,pt2,q2kk) !formerly in xEmsPx(
c-----------------------------------------------------------------------
c plot  x-distribution and y-distribution of Pomerons
c-----------------------------------------------------------------------

#include "aaa.h"
#include "ems.h"
      common/geom/rmproj,rmtarg,bmax,bkmx
      parameter(nbix=100,nbib=51,nbiy=50)
      common/cemspxwxx/
     *   wxmcxx(2,nbix),wxmcQpxx(2,nbix),wxmcQtxx(2,nbix),wwmcxx(2,nbix)
     *  ,wxpxx(2,nbix),wxmxx(2,nbix),wxpQpxx(2,nbix)
     *  ,wxmQpxx(2,nbix),wxpQtxx(2,nbix),wxmQtxx(2,nbix)
     *  ,wymcxx(nbiy),wymcQtxx(2,nbiy),wymcQpxx(2,nbiy),nypxx,nymxx
     *  ,womcxx(2,nbix)
      common/cemspx/xu,xo,yu,yo,dy,xlu,xlo,bb,nn,db,mm,nm,nt
      common/cemspxq2/q2kkk1(2),nq2kkk1
      dimension q2kk(2)

      if(iemspx.eq.0)call utstop('ERROR in XemsPx: iemspx = 0&')

      if(iii.eq.-1)then

        do k=1,2
        do i=1,nbix
          wxmcxx(k,i) =0
          womcxx(k,i) =0
          wxmcQpxx(k,i)=0
          wxmcQtxx(k,i)=0
          wwmcxx(k,i) =0
          wxpxx(k,i)  =0
          wxpQpxx(k,i) =0
          wxpQtxx(k,i) =0
          wxmxx(k,i)  =0
          wxmQpxx(k,i) =0
          wxmQtxx(k,i) =0
        enddo
        enddo
        do i=1,nbiy
          wymcxx(i)    =0
          wymcQpxx(1,i)   =0
          wymcQtxx(1,i)   =0
          wymcQpxx(2,i)   =0
          wymcQtxx(2,i)   =0
        enddo
        nypxx  =0
        nymxx  =0
        q2kkk1(1)=0
        q2kkk1(2)=0
        nq2kkk1=0


      elseif(iii.eq.1)then

       xp=sqrt(xmc)*exp(ymc)
       xm=sqrt(xmc)*exp(-ymc)
       mm=mm+1

c       if(npos.ne.3)goto 11
       if(xmc.lt.xu)goto11
       i=1+int(alog(xmc/xu)/alog(xo/xu)*nbix)
       if(i.gt.nbix)goto1
       if(i.lt.1)goto1
       if(jij.eq.1)             wxmcxx(1,i)=wxmcxx(1,i)+1. !avoid double counting
       if(jij.eq.1.and.jot.ge.1)womcxx(1,i)=womcxx(1,i)+1. !avoid double counting
       if(npos.eq.3)then      !Qs vs x_pom
       wxmcQpxx(1,i)=wxmcQpxx(1,i)+q2kk(1)
       wxmcQtxx(1,i)=wxmcQtxx(1,i)+q2kk(2)
       wxmcQpxx(2,i)=wxmcQpxx(2,i)+1.
       wxmcQtxx(2,i)=wxmcQtxx(2,i)+1.
       endif
       if(pt1.gt.100.or.pt2.gt.100)
     .  wwmcxx(1,i)=wwmcxx(1,i)+1.
1      continue
       i=1+int((xmc-xu)/(xo-xu)*nbix)
       if(i.gt.nbix)goto11
       if(i.lt.1)goto11
       if(jij.eq.1)             wxmcxx(2,i)=wxmcxx(2,i)+1. !avoid double counting
       if(jij.eq.1.and.jot.ge.1)womcxx(2,i)=womcxx(2,i)+1. !avoid double counting
       if(pt1.gt.100.or.pt2.gt.100)
     .  wwmcxx(2,i)=wwmcxx(2,i)+1.
11     continue

       if(jij.gt.1)goto 12         !avoid double counting
       if(xp.lt.xlu)goto12
       i=1+int(alog(xp/xlu)/alog(xlo/xlu)*nbix)
       if(i.gt.nbix)goto2
       if(i.lt.1)goto2
       wxpxx(1,i)=wxpxx(1,i)+1.
2      continue
       i=1+int((xp-xlu)/(xlo-xlu)*nbix)
       if(i.gt.nbix)goto12
       if(i.lt.1)goto12
       wxpxx(2,i)=wxpxx(2,i)+1.
12     continue

       if(jij.gt.1)goto 13         !avoid double counting
       if(xm.lt.xlu)goto13
       i=1+int(alog(xm/xlu)/alog(xlo/xlu)*nbix)
       if(i.gt.nbix)goto3
       if(i.lt.1)goto3
       wxmxx(1,i)=wxmxx(1,i)+1.
3      continue
       i=1+int((xm-xlu)/(xlo-xlu)*nbix)
       if(i.gt.nbix)goto13
       if(i.lt.1)goto13
       wxmxx(2,i)=wxmxx(2,i)+1.
13     continue

       if(npos.eq.3)then      !Qs vs xp and xm
       if(xp.lt.xlu)goto 14
       i=1+int(alog(xp/xlu)/alog(xlo/xlu)*nbix)
       if(i.gt.nbix)goto 14
       if(i.lt.1)goto 14
       wxpQpxx(1,i)=wxpQpxx(1,i)+q2kk(1)
       wxpQpxx(2,i)=wxpQpxx(2,i)+1.
 14    continue

       if(xm.lt.xlu)goto 15
       i=1+int(alog(xm/xlu)/alog(xlo/xlu)*nbix)
       if(i.gt.nbix)goto 15
       if(i.lt.1)goto 15
       wxmQtxx(1,i)=wxmQtxx(1,i)+q2kk(2)
       wxmQtxx(2,i)=wxmQtxx(2,i)+1.
 15    continue

       if(xpl.lt.xlu)goto 16
       i=1+int(alog(xpl/xlu)/alog(xlo/xlu)*nbix)
       if(i.gt.nbix)goto 16
       if(i.lt.1)goto 16
       wxpQtxx(1,i)=wxpQtxx(1,i)+max(q2kk(1),q2kk(2))
       wxpQtxx(2,i)=wxpQtxx(2,i)+1.
 16    continue

       if(xml.lt.xlu)goto 17
       i=1+int(alog(xml/xlu)/alog(xlo/xlu)*nbix)
       if(i.gt.nbix)goto 17
       if(i.lt.1)goto 17
       wxmQpxx(1,i)=wxmQpxx(1,i)+max(q2kk(1),q2kk(2))
       wxmQpxx(2,i)=wxmQpxx(2,i)+1.
 17    continue
       endif

       if(ymc.lt.yu)return
       i=int((ymc-yu)/dy)+1
       if(i.gt.nbiy)i=nbiy !return
       if(i.lt.1)i=1       !return
       wymcxx(i)=wymcxx(i)+1
       if(npos.eq.3)then
       wymcQpxx(1,i)=wymcQpxx(1,i)+q2kk(1)
       wymcQtxx(1,i)=wymcQtxx(1,i)+q2kk(2)
       wymcQpxx(2,i)=wymcQpxx(2,i)+1.
       wymcQtxx(2,i)=wymcQtxx(2,i)+1.
       endif
       if(ymc.gt.0)nypxx=nypxx+1
       if(ymc.lt.0)nymxx=nymxx+1

       if(npos.eq.3)then
         q2kkk1(1)=q2kkk1(1)+q2kk(1)
         q2kkk1(2)=q2kkk1(2)+q2kk(2)
         nq2kkk1=nq2kkk1+1
       endif

       else

          write(ifch,*)'\n\n  ERROR 27092012 \n\n'
                   stop'\n\n  ERROR 27092012 \n\n'
       endif
       end   !xxEmsPx

c-----------------------------------------------------------------------
      subroutine xEmsPx(iii,xmc,ymc,xpl,xml,npos,pt1,pt2)
c-----------------------------------------------------------------------
c plot  x-distribution and y-distribution of Pomerons
c-----------------------------------------------------------------------

#include "aaa.h"
#include "par.h"
#include "ems.h"
      common/geom/rmproj,rmtarg,bmax,bkmx

      parameter(nbix=100,nbib=51,nbiy=50,mmxcentrx=9)
      common/cxw/x(2,nbix),dx(2,nbix),wxmc(0:mmxcentrx+1,2,nbix)
     . ,wxmcQp(0:mmxcentrx+1,2,nbix),xl(2,nbix),dxl(2,nbix)
     . ,wxp(0:mmxcentrx+1,2,nbix),wxm(0:mmxcentrx+1,2,nbix)
     . ,wxpQp(0:mmxcentrx+1,2,nbix),wxmQp(0:mmxcentrx+1,2,nbix)
     . ,wxpQt(0:mmxcentrx+1,2,nbix),wxmQt(0:mmxcentrx+1,2,nbix)
     . ,wxmcQt(0:mmxcentrx+1,2,nbix),wwmc(0:mmxcentrx+1,2,nbix)
     . ,womc(0:mmxcentrx+1,2,nbix)
      common/cyw/y(nbiy),wymc(0:mmxcentrx+1,nbiy)
     .,wymcQt(0:mmxcentrx+1,2,nbiy),wymcQp(0:mmxcentrx+1,2,nbiy),nyp,nym
      double precision PomIncXExact,PomIncPExact,PomIncMExact,Pominc
      double precision PomIncXIExact,PomIncPIExact,PomIncMIExact
      external PomIncXExact,PomIncPExact,PomIncMExact
      external PomIncXIExact,PomIncPIExact,PomIncMIExact
      common/cncoll/wncoll(3,0:mmxcentrx+1),wnevt(3,0:mmxcentrx+1)
     .             ,wbim(0:mmxcentrx+1),wpom(0:7,0:mmxcentrx+1)
      common/cwq2/wq2(3,0:mmxcentrx+1,3)
      common/cemspx/xu,xo,yu,yo,dy,xlu,xlo,bb,nn,db,mm,nm,nt
      character modu*5, imod*5, txtxm*6
      common/cicentrality/icentrality
      real q2amin(2),zu(nbix)
      logical lppb
      real ff(5)

      dummy=xmc
      dummy=ymc
      dummy=xpl
      dummy=xml
      dummy=npos
      dummy=pt1
      dummy=pt2
      iG100=0
      if(iemspx.eq.0)call utstop('ERROR in XemsPx: iemspx = 0&')
      if(mmxcentrf().gt.mmxcentrx)stop'\n\n ERROR 27092012c \n\n '

      if(iii.eq.0)then

       !write(ifch,*)'initialize xEmsPx tables for'
       !. ,mmxcentrf(),' centralities'
       xu=15./engy**2 !15=empirical value, seems to work for all energies
       xo=1.
       xlu=15./engy**2
       xlo=1.
       yu=-alog(engy**2)
       yo=alog(engy**2)
       dy=(yo-yu)/nbiy
       do i=1,nbix
        x(1,i)=xu*(xo/xu)**((i-0.5)/nbix)
        x(2,i)=xu+(xo-xu)*((i-0.5)/nbix)
        dx(1,i)=xu*(xo/xu)**(1.*i/nbix)*(1.-(xo/xu)**(-1./nbix))
        dx(2,i)=(xo-xu)/nbix
        do n=0,mmxcentrf()
        wxmc(n,1,i)=0.
        wxmc(n,2,i)=0.
        womc(n,1,i)=0.
        womc(n,2,i)=0.
        wxmcQp(n,1,i)=0.
        wxmcQp(n,2,i)=0.
        wxmcQt(n,1,i)=0.
        wxmcQt(n,2,i)=0.
        wwmc(n,1,i)=0.
        wwmc(n,2,i)=0.
        enddo
       enddo
       do i=1,nbix
        xl(1,i)=xlu*(xlo/xlu)**((i-0.5)/nbix)
        xl(2,i)=xlu+(xlo-xlu)*((i-0.5)/nbix)
        dxl(1,i)=xlu*(xlo/xlu)**(1.*i/nbix)*(1.-(xlo/xlu)**(-1./nbix))
        dxl(2,i)=(xlo-xlu)/nbix
        do n=0,mmxcentrf()
        wxp(n,1,i)=0.
        wxp(n,2,i)=0.
        wxm(n,1,i)=0.
        wxm(n,2,i)=0.
        wxpQp(n,1,i)=0.
        wxpQp(n,2,i)=0.
        wxmQp(n,1,i)=0.
        wxmQp(n,2,i)=0.
        wxpQt(n,1,i)=0.
        wxpQt(n,2,i)=0.
        wxmQt(n,1,i)=0.
        wxmQt(n,2,i)=0.
        enddo
       enddo
       do i=1,nbiy
        y(i)=yu+dy/2.+float(i-1)*dy
        do n=0,mmxcentrf()
        wymc(n,i)=0.
        wymcQp(n,1,i)=0.
        wymcQp(n,2,i)=0.
        wymcQt(n,1,i)=0.
        wymcQt(n,2,i)=0.
        enddo
       enddo
       mm=0
       nt=0
       nyp=0
       nym=0
       db=bkmx*2./float(nbib-1)
       do n=0,mmxcentrf()
        wncoll(1,n)=0.
        wnevt(1,n)=0.
        wbim(n)=0.
        do j=0,7
          wpom(j,n)=0.
        enddo
        do j=1,3
          wq2(1,n,j)=0
        enddo
       enddo

      elseif(iii.eq.2)then

       lppb=.false.

       if(mmxcentrf().gt.0) then
        ncy1=1
       else
        ncy1=0
       endif
       ncy2=mmxcentrf()

       if(nint(xpar7).eq.1)then
         ncy1=0
         ncy2=0
       endif
       kollini=koll

       do ncy=ncy1,ncy2  !-------------ncy------------->

       nstat=nfull
       if(nstat.lt.0)nstat=max(1,nevent)
       hw=nstat
       zavpoa=1.
       zavpom=1.
       zavco1=0.
       zavco2=0.
       zavco5=0.
       zavco6=0.
       bb=0.
       q2amin(1)=0
       q2amin(2)=0
       if(wnevt(1,ncy).gt.0.)then
         bb=wbim(ncy)/wnevt(1,ncy)
       endif
       if(wpom(0,ncy).gt.0.)then
         zavpoa=wpom(4,ncy)/wpom(0,ncy)
         zavpom=wpom(3,ncy)/wpom(0,ncy)
         zavco1=wpom(1,ncy)/wpom(0,ncy)
         zavco2=wpom(2,ncy)/wpom(0,ncy)
         zavco5=wpom(5,ncy)/wpom(0,ncy)
         zavco6=wpom(6,ncy)/wpom(0,ncy)
         zavco7=wpom(7,ncy)/wpom(0,ncy)
       endif
       if(wq2(1,ncy,3).gt.0.)then
         q2amin(1)=wq2(1,ncy,1)/wq2(1,ncy,3)  !nuclear component of pair averaged Q2s (no x)
         q2amin(2)=wq2(1,ncy,2)/wq2(1,ncy,3)
       endif

c       q2min1=0.     !should be 0 for proper use of PomInc
c       q2min2=0.
       !cKW22 fix "add" problem
       !former ff is now ff(1)*ff(2)/ff(3), ff(i) put into histo
       ff(1)=0. !NOT USED
       ff(2)=0. !NOT USED
       ff(3)=0.
       ff(4)=0. !NOT USED
       ff(5)=0.
       if(maproj.eq.1.and.matarg.eq.1.and.bminim.eq.bmaxim)then
        mmmm=1
        bb=bmaxim
        ff(3)=float(ntevt)
        imod='   dn'
       elseif(maproj.eq.1.and.matarg.eq.1)then
        mmmm=3
        ff(3)=wnevt(1,ncy)
        if(ncy.gt.0)hw=wnevt(1,ncy)
        lppb=jzmode(3).eq.1 !centrality variable = bim
        imod='   dn'
        !q2pmin fixed in function direclty (averaged over b)
       else
        mmmm=2
        if(ncy.gt.0)then
          if(wnevt(1,ncy).gt.0.)then
            ff(3)=wncoll(1,ncy) 
          endif
          hw=wnevt(1,ncy)
        elseif(bminim.lt.0.001.and.bmaxim.gt.20)then     !min bias
          area=pi*(rmproj+rmtarg)**2
          ff(3)=wncoll(1,ncy)
        else
          write(ifmt,*)'ignored'
          return
        endif
        imod='   dn'
       endif
        !only averaged nuclear component has to be added in functions
        q2pmin(1)=max(q2nmin,q2amin(1))
        q2pmin(2)=max(q2nmin,q2amin(2))
        call ipoOm5Tables(1)

       if(ncy1.eq.1.and.wnevt(1,ncy).lt.0.0001)then
         hww=0
         ihww=0
       else
         hww=1
         ihww=1
       endif

       koll=1

       kk=1
       modu=' log '

       write(ifhi,'(a)') "!-----------"  !-------------------------------------------------------------
       write(ifhi,'(a)') "! MC Norm   "  !cKW22 Bugfix: ff must go to histogram (to use "add" properly) 
       write(ifhi,'(a)') "!-----------"  !-------------------------------------------------------------
       write(ifhi,'(a,i1)')'openhisto name PxNorm3L'//modu(3:4),ncy
       write(ifhi,'(a)')   'histoweight -1' !,dble(hw)
       write(ifhi,'(a)')  'array 2'
       write(ifhi,'(4e11.3)')1.,ff(3)
       write(ifhi,'(a)')  'endarray'
       write(ifhi,'(a)')  'closehisto'
       ff(5)=0
       do i=1,nbix
        ff(5)=ff(5)+wxmc(ncy,kk,i)
       enddo
       write(ifhi,'(a,i1)')'openhisto name PxNorm5L'//modu(3:4),ncy
       write(ifhi,'(a)')   'histoweight -1' !,dble(hw)
       write(ifhi,'(a)')  'array 2'
       write(ifhi,'(4e11.3)')1.,ff(5)
       write(ifhi,'(a)')  'endarray'
       write(ifhi,'(a)')  'closehisto'

       write(ifhi,'(a)')       '!----------------------------------'
       write(ifhi,'(a)')       '!   Q2s x distribution    '//modu
       write(ifhi,'(a)')       '!----------------------------------'

       write(ifhi,'(a,i1)')  'openhisto name xQ2pSimuL'//modu(3:4),ncy
       write(ifhi,'(a)')  'htyp lru xmod'//modu//'ymod lin'
       write(ifhi,'(a,2e11.3)')'xrange',xu,xo
       write(ifhi,'(a)')    'text 0 0 "xaxis x?PE!"'
       write(ifhi,'(a)') 'text 0 0 "yaxis Q2s? pro / tar!"'
       write(ifhi,'(a,i2,a)')'text 0.6 0.85  "M',ncy,'"'
       write(ifhi,'(a,e22.14)')'histoweight ',dble(hw)
       write(ifhi,'(a)')       ' array 2'
       do i=1,nbix
       u=x(kk,i)
       z=q2nmin
       if(wxmcQp(ncy,2,i).gt.0.)z=wxmcQp(ncy,kk,i)/wxmcQp(ncy,2,i)
        write(ifhi,'(2e11.3)')u,z
       enddo
       write(ifhi,'(a)')    '  endarray'
       write(ifhi,'(a)')    'closehisto plot 0-'

       write(ifhi,'(a,i1)')  'openhisto name xQ2tSimuL'//modu(3:4),ncy
       write(ifhi,'(a)')  'htyp lba xmod'//modu//'ymod log'
       write(ifhi,'(a,2e11.3)')'xrange',xu,xo
       write(ifhi,'(a,i2,a)')'text 0.6 0.85  "M',ncy,'"'
       write(ifhi,'(a,e22.14)')'histoweight ',dble(hw)
       write(ifhi,'(a)')       ' array 2'
       do i=1,nbix
       u=x(kk,i)
       z=q2nmin
       if(wxmcQt(ncy,2,i).gt.0.)z=wxmcQt(ncy,kk,i)/wxmcQt(ncy,2,i)
        write(ifhi,'(2e11.3)')u,z
       enddo
       write(ifhi,'(a)')    '  endarray'
       write(ifhi,'(a)')    'closehisto plot 0'

       write(ifhi,'(a)')           '!--------------------------------'
       write(ifhi,'(a)')           '!   y Q2s distribution   '//modu
       write(ifhi,'(a)')           '!--------------------------------'

       write(ifhi,'(a,i1)')    'openhisto name yQ2pSimuL'//modu(3:4),ncy
       write(ifhi,'(a)')       'htyp lru xmod lin ymod lin'
       write(ifhi,'(a,2e11.3)')'xrange',yu,yo
       write(ifhi,'(a)')    'text 0 0 "xaxis y?PE!"'
       write(ifhi,'(a)') 'text 0 0 "yaxis  Q2s? pro / tar! "'
       write(ifhi,'(a,e22.14)')'histoweight ',dble(hw)
       write(ifhi,'(a)')       ' array 2'
       s1=0
       do i=1,nbiy
       u=y(i)
       z=q2nmin
       if(wymcQp(ncy,2,i).gt.0.)z=wymcQp(ncy,kk,i)/wymcQp(ncy,2,i)
       s1=s1+z*dy
        write(ifhi,'(2e11.3)')u,z
       enddo
       write(ifhi,'(a)')    '  endarray'
       write(ifhi,'(a,f7.3,a)')
     *                       'text .1 .15 "Ip= ',s1/(y(nbiy)-y(1)),'"'
       write(ifhi,'(a)')    'closehisto plot 0-'

       write(ifhi,'(a,i1)')    'openhisto name yQ2tSimuL'//modu(3:4),ncy
       write(ifhi,'(a)')       'htyp lba xmod lin ymod'//modu
       write(ifhi,'(a,2e11.3)')'xrange',yu,yo
       write(ifhi,'(a,e22.14)')'histoweight ',dble(hw)
       write(ifhi,'(a)')       ' array 2'
       s1=0
       do i=1,nbiy
       u=y(i)
       z=q2nmin
       if(wymcQt(ncy,2,i).gt.0.)z=wymcQt(ncy,kk,i)/wymcQt(ncy,2,i)
       s1=s1+z*dy
        write(ifhi,'(2e11.3)')u,z
       enddo
       write(ifhi,'(a)')    '  endarray'
       write(ifhi,'(a,f7.3,a)')
     *                      'text .65 .15 "It= ',s1/(y(nbiy)-y(1)),'"'
       write(ifhi,'(a)')    'closehisto plot 0'

       write(ifhi,'(a)')       '!----------------------------------'
       write(ifhi,'(a)')       '!   Q2s x+ distribution    '//modu
       write(ifhi,'(a)')       '!----------------------------------'

       write(ifhi,'(a,i1)')'openhisto name xpQ2SimuL'//modu(3:4),ncy
       write(ifhi,'(a)')   'htyp lru xmod'//modu//'ymod log'
       write(ifhi,'(a,2e11.3)')'xrange',xlu,xlo
       write(ifhi,'(a)')  'text 0 0 "xaxis x+?PE! or x-?IB!" '
       write(ifhi,'(2a)') 'text 0 0 "yaxis Q2s?pro! or Q2s?tar!"'
       write(ifhi,'(a,e22.14)')'histoweight ',dble(hw)
       write(ifhi,'(a)')       ' array 2'
       do i=1,nbix
       u=xl(kk,i)
       z=q2nmin
       if(wxpQp(ncy,2,i).gt.0.)z=wxpQp(ncy,kk,i)/wxpQp(ncy,2,i)
       write(ifhi,'(2e11.3)')u,z
       enddo
       write(ifhi,'(a)')    '  endarray'
       write(ifhi,'(a)')    'closehisto plot 0-'


       write(ifhi,'(a)')       '!----------------------------------'
       write(ifhi,'(a)')       '!   Q2s x- distribution    '//modu
       write(ifhi,'(a)')       '!----------------------------------'

       write(ifhi,'(a,i1)')'openhisto name xmQ2SimuL'//modu(3:4),ncy
       write(ifhi,'(a)')   'htyp lba xmod'//modu//'ymod log'
       write(ifhi,'(a,2e11.3)')'xrange',xlu,xlo
       write(ifhi,'(a,e22.14)')'histoweight ',dble(hw)
       write(ifhi,'(a)')       ' array 2'
       do i=1,nbix
       u=xl(kk,i)
       z=q2nmin
       if(wxmQp(ncy,2,i).gt.0.)z=wxmQp(ncy,kk,i)/wxmQp(ncy,2,i)
        write(ifhi,'(2e11.3)')u,z
       enddo
       write(ifhi,'(a)')    '  endarray'
       write(ifhi,'(a)')    'closehisto plot 0'

       kk1=nint(xpar1)
       kk2=nint(xpar2)
       mode=nint(xpar3)

       do kk=kk1,kk2

       if(kk.eq.1)modu=' log '
       if(kk.eq.2)modu=' lin '


       write(ifhi,'(a)')       '!----------------------------------'
       write(ifhi,'(a)')       '!   Pomeron x distribution    '//modu
       write(ifhi,'(a)')       '!----------------------------------'

       write(ifhi,'(a)')  'openhisto'
       write(ifhi,'(a)')  'xmod'//modu//'ymod log'
       write(ifhi,'(a,2e11.3)')'xrange',xu,xo
       write(ifhi,'(a)')  'yrange 1e-4 1e7'
       write(ifhi,'(a)')    'text 0 0 "xaxis x?PE!"'
       write(ifhi,'(a)') 'text 0 0 "yaxis'//imod//'?Pom! / dx?PE!"'
       write(ifhi,'(a,i2,a)')'text 0.03 0.89  "M',ncy,'"'
       write(ifhi,'(a)')'text 0.0 1.03 "red:MC"'  

       write(ifhi,'(a,i1)')    'name xNorSimuL'//modu(3:4),ncy
       write(ifhi,'(a)')       'htyp lyv' 
       write(ifhi,'(a)')       'histoweight -1' !,dble(hw)
       write(ifhi,'(a)')       ' array 2'
       do i=1,nbix
        u=x(kk,i)
        z=womc(ncy,kk,i)/dx(kk,i) !---Nor--- (no Sat)
        write(ifhi,'(2e11.3)')u,z
       enddo
       write(ifhi,'(a)')    '  endarray'
       !-----------------
       !multiply with ff
       !-----------------
       write(ifhi,'(a,i1,$)')
     . 'calc *0 / PxNorm5L'//modu(3:4),ncy
       write(ifhi,'(a)')' ; '
       !---------------
       write(ifhi,'(a)')    'closehisto plot 0-'

       write(ifhi,'(a)')'openhisto' 
       write(ifhi,'(a,i1)')    'name xPomSimuL'//modu(3:4),ncy
       write(ifhi,'(a)')       'htyp lru xmod'//modu//'ymod log'
       write(ifhi,'(a,2e11.3)')'xrange',xu,xo
       write(ifhi,'(a)')       'histoweight -1' !,dble(hw)
       write(ifhi,'(a)')       'array 2'
       do i=1,nbix
       u=x(kk,i)
       z=wxmc(ncy,kk,i)/dx(kk,i) !----All----
        write(ifhi,'(2e11.3)')u,z
       enddo
       write(ifhi,'(a)')    '  endarray'
       !-----------------
       !multiply with ff
       !-----------------
       write(ifhi,'(a,i1,$)')
     . 'calc *0 / PxNorm5L'//modu(3:4),ncy
       write(ifhi,'(a)')' ; '
       !---------------
       write(ifhi,'(a)')    'closehisto plot 0-'


       if(iG100.eq.1)then!~~~~~~~~~~
       write(ifhi,'(a,i1,a)')
     . 'openhisto name xPomSimuL'//modu(3:4),ncy,'G100 set factor 1e5'
       write(ifhi,'(a)') 'htyp lkv xmod'//modu//'ymod log'
       write(ifhi,'(a,2e11.3)')'xrange',xu,xo
       write(ifhi,'(a)')       'histoweight -1' !,dble(hw)
       write(ifhi,'(a)')       'array 2'
       do i=1,nbix
       u=x(kk,i)
       z=wwmc(ncy,kk,i)/dx(kk,i)
        write(ifhi,'(2e11.3)')u,z
       enddo
       write(ifhi,'(a)')    '  endarray'
       !-----------------
       !multiply with ff
       !-----------------
       write(ifhi,'(a,i1,$)')
     . 'calc *0 / PxNorm3L'//modu(3:4),ncy
       write(ifhi,'(a)')' ; '
       !---------------
       write(ifhi,'(a)')    'closehisto plot 0-'
       endif!~~~~~~~~~~

       write(ifhi,'(a,i1)')'openhisto name xPomTheoL'//modu(3:4),ncy
       write(ifhi,'(a)')  'htyp lba xmod'//modu//'ymod log'
       write(ifhi,'(a,2e11.3)')'xrange',xu,xo
       write(ifhi,'(a)')'text 0.35 1.03 "blue:Theo"'  
       write(ifhi,'(a)')    'text 0 0 "xaxis x?PE!"'
       write(ifhi,'(a)') 'text 0 0 "yaxis'//imod//'?Pom! / dx?PE!"'
       write(ifhi,'(a,e22.14)')'histoweight ',dble(hww)
       write(ifhi,'(a)')       'array 2'
       s2=0.
       do i=1,nbix
        u=x(kk,i)
        zu(i)=0
        if(ihww.gt.0)then
          if(ncy.gt.0.and.lppb)then
            zu(i)=sngl(PomIncXExact(dble(u),bb))
          else
            if(mmmm.eq.1)zu(i)=sngl(PomIncXExact(dble(u),bb))
            if(mmmm.eq.2)zu(i)=sngl(PomIncXIExact(dble(u)))/sigine*10
            if(mmmm.eq.3)zu(i)=sngl(PomIncXIExact(dble(u)))/sigine*10
          endif
          s2=s2+dx(kk,i)*zu(i)
        endif
       enddo
       do i=1,nbix
         if(ncy.gt.0.and.s2.gt.0.)then
           zu(i)=zu(i)/s2    !to have the same normalization than in MC (per Pom)
         endif
         write(ifhi,'(2e11.3)')x(kk,i),zu(i)
       enddo
       write(ifhi,'(a)')    '  endarray'
       write(ifhi,'(a,f7.3,a,f7.3,a)')
     *                       '!text .1 .75 "I=',s1,' (',s2,')"'
       write(ifhi,'(a,f7.3,a)')'!text .1 .85 "x+ : av=',avxp,'"'
c       if(ncy.gt.0.or.iokoll.ne.0)then
       write(ifhi,'(a)')    'closehisto plot 0-'

       write(ifhi,'(a,i1)')   'openhisto name xPomFitL'//modu(3:4),ncy
       write(ifhi,'(a)')   'htyp lgi xmod'//modu//'ymod log'
       write(ifhi,'(a,2e11.3)')'xrange',xlu,xlo
       write(ifhi,'(a)')'text 0.8 1.03 "green:Fit"'  
       write(ifhi,'(a)')    'text 0 0 "xaxis x+?PE!"'
       write(ifhi,'(a)') 'text 0 0 "yaxis'//imod//'?Pom! / dx?PE!"'
       write(ifhi,'(a,e22.14)')'histoweight ',dble(hww)
       write(ifhi,'(a)')       'array 2'
       s3=0
       call  getSystemType(isys,amassMax,amassAsy)
       call deformoptget(1,iodeform)  
       if(iodeform.eq.1)stop'ERROR iodeform = 1 not supported any more'
       conn1=zavco5
       conn2=zavco6
       jval=igetJval(maproj,matarg,engy)
       if(jval.eq.1)then
         zval=(conn1+conn2)/2
       elseif(jval.eq.2)then
         zval=zavco7
       else
         stop'ERROR 09062022d'
       endif
       do i=1,nbix
        u=x(kk,i)
        defrm=deform(maproj,matarg,engy,zval,u) 
        if(defrm.gt.100.)then
          print*,'ERROR defrm',defrm,engy,zval,u
          stop'forced STOP'
        endif
        zu(i)=0.
        if(ihww.gt.0)then
          if(ncy.gt.0.and.lppb)then !pp & centrality variable = bim
            zu(i)=sngl(PomIncXExact(dble(u),bb))*defrm 
          else
            if(mmmm.eq.1)then!pp bminim EQ bmaxim
              zu(i)=sngl(PomIncXExact(dble(u),bb))*defrm
            elseif(mmmm.eq.2)then!NOTpp
              zu(i)=sngl(PomIncXIExact(dble(u)))/sigine*10*defrm 
            elseif(mmmm.eq.3)then!pp bminim NE bmaxim
              zu(i)=sngl(PomIncXIExact(dble(u)))/sigine*10*defrm 
            endif
          endif
          s3=s3+dx(kk,i)*zu(i)
         endif
       enddo
       !print*,'xEmsPx PomInc normalization:',s3
       s3=s2 ! Use same normalization as Theo 
       do i=1,nbix
         if(s3.gt.0.)then
           if(abs(iokoll).gt.0)then
             zu(i)=zu(i)/s3*float(abs(iokoll)) !to have the same normalization than in MC
           elseif(ncy.gt.0)then
             zu(i)=zu(i)/s3     !to have the same normalization than in MC (per Pom)
           endif
         endif
         write(ifhi,'(2e11.3)')x(kk,i),zu(i)
       enddo
       if(ncy.gt.0)s3=1.    !by definition
       write(ifhi,'(a)')    '  endarray'
       write(ifhi,'(a,f7.3,a)')'!text .6 .65 "(',s3,')"'
c       endif
       write(ifhi,'(a)')    'closehisto plot 0'

       if(mode.eq.1)goto 777 !xpar3

       if(kk.eq.kk1)then
       write(ifhi,'(a)')           '!--------------------------------'
       write(ifhi,'(a)')           '!   Pomeron y distribution   '//modu
       write(ifhi,'(a)')           '!--------------------------------'

       write(ifhi,'(a,i1)')  'openhisto name yPomSimuL'//modu(3:4),ncy
       write(ifhi,'(a)')       'htyp lru xmod lin ymod'//modu
       write(ifhi,'(a,2e11.3)')'xrange',yu,yo
       write(ifhi,'(a)')    'text 0 0 "xaxis y?PE!"'
       write(ifhi,'(a)') 'text 0 0 "yaxis'//imod//'?Pom!/dy?PE!"'
       write(ifhi,'(a)')       'histoweight -1' !,dble(hw)
       write(ifhi,'(a)')       ' array 2'
       s1=0
       do i=1,nbiy
       u=y(i)
       z=wymc(ncy,i)/dy
       s1=s1+z*dy
        write(ifhi,'(2e11.3)')u,z
       enddo
       write(ifhi,'(a)')    '  endarray'
       !-----------------
       !multiply with ff
       !-----------------
       write(ifhi,'(a,i1,$)')
     . 'calc *0 / PxNorm3L'//modu(3:4),ncy
       write(ifhi,'(a)')' ; '
       !---------------
       write(ifhi,'(a,f7.3,a)')
     *                       'text .1 .15 "I= ',s1,'"'
       write(ifhi,'(a)')    'closehisto plot 0'
       endif

       write(ifhi,'(a)')       '!----------------------------------'
       write(ifhi,'(a)')       '!   Pomeron x+ distribution    '//modu
       write(ifhi,'(a)')       '!----------------------------------'

       write(ifhi,'(a,i1)')'openhisto name xpPomSimuL'//modu(3:4),ncy
       write(ifhi,'(a)')   'htyp lru xmod'//modu//'ymod log'
       write(ifhi,'(a,2e11.3)')'xrange',xlu,xlo
       write(ifhi,'(a)')    'text 0 0 "xaxis x+?PE!"'
       write(ifhi,'(a)') 'text 0 0 "yaxis'//imod//'?Pom! / dx+?PE!"'
       write(ifhi,'(a,f7.3,a)')'text 0.1 0.85 "Con. av',zavco1,'"'
       write(ifhi,'(a)')       'histoweight -1' !,dble(hw)
       write(ifhi,'(a)')       ' array 2'
       s1=0
       do i=1,nbix
       u=xl(kk,i)
       z=wxp(ncy,kk,i)/dxl(kk,i)
       s1=s1+z*dxl(kk,i)
        write(ifhi,'(2e11.3)')u,z
       enddo
       write(ifhi,'(a)')    '  endarray'
       !-----------------
       !multiply with ff
       !-----------------
       write(ifhi,'(a,i1,$)')
     . 'calc *0 / PxNorm3L'//modu(3:4),ncy
       write(ifhi,'(a)')' ; '
       !---------------
       write(ifhi,'(a)')    'closehisto plot 0-'

       write(ifhi,'(a,i1)')   'openhisto name xpPomUnitL'//modu(3:4),ncy
       write(ifhi,'(a)')   'htyp lba xmod'//modu//'ymod log'
       write(ifhi,'(a,2e11.3)')'xrange',xlu,xlo
       write(ifhi,'(a)')    'text 0 0 "xaxis x+?PE!"'
       write(ifhi,'(a)') 'text 0 0 "yaxis'//imod//'?Pom! / dx+?PE!"'
       write(ifhi,'(a,e22.14)')'histoweight ',dble(hww)
       write(ifhi,'(a)')       'array 2'
       s2=0
       do i=1,nbix
        u=xl(kk,i)
        zu(i)=0
        if(ihww.gt.0)then
          if(ncy.gt.0.and.lppb)then
            zu(i)=sngl(PomIncPExact(dble(u),bb))
          else
            if(mmmm.eq.1)zu(i)=sngl(PomIncPExact(dble(u),bb))
            if(mmmm.eq.2)zu(i)=sngl(PomIncPIExact(dble(u)))/sigine*10
            if(mmmm.eq.3)zu(i)=sngl(PomIncPIExact(dble(u)))/sigine*10
          endif
          s2=s2+dxl(kk,i)*zu(i)
        endif
       enddo
       do i=1,nbix
         if(ncy.gt.0.and.s2.gt.0.)then
           zu(i)=zu(i)/s2    !to have the same normalization than in MC (per Pom)
         endif
         write(ifhi,'(2e11.3)')x(kk,i),zu(i)
       enddo
       if(ncy.gt.0)s2=1.    !by definition
       write(ifhi,'(a)')    '  endarray'
       write(ifhi,'(a,f7.3,a,f7.3,a)')
     *                       'text .1 .15 "I= ',s1,' (',s2,')"'
       
       write(ifhi,'(a)')    'closehisto plot 0-'

       write(ifhi,'(a,i1)')   'openhisto name xpPomSatL'//modu(3:4),ncy
       write(ifhi,'(a)')   'htyp lgi xmod'//modu//'ymod log'
       write(ifhi,'(a,2e11.3)')'xrange',xlu,xlo
       write(ifhi,'(a)')    'text 0 0 "xaxis x+?PE!"'
       write(ifhi,'(a)') 'text 0 0 "yaxis'//imod//'?Pom! / dx+?PE!"'
       write(ifhi,'(a,e22.14)')'histoweight ',dble(hww)
       write(ifhi,'(a)')       'array 2'
       s3=0
       do i=1,nbix
        u=xl(kk,i)
        zu(i)=0.
        if(ihww.gt.0)then
          if(ncy.gt.0.and.lppb)then
      zu(i)=sngl(PomInc(PomIncPExact,dble(u),bb,zavco1,zavco2,zavpom,2))
          else
            if(mmmm.eq.1)zu(i)=
     .    sngl(PomInc(PomIncPExact,dble(u),bb,zavco1,zavco2,zavpom,2))
            if(mmmm.eq.2)zu(i)=10./sigine
     .   *sngl(PomInc(PomIncPIExact,dble(u),0.,zavco1,zavco2,zavpom,-2))
            if(mmmm.eq.3)zu(i)=10./sigine
     .   *sngl(PomInc(PomIncPIExact,dble(u),0.,zavco1,zavco2,zavpom,-2))
          endif
          s3=s3+dxl(kk,i)*zu(i)
         endif
       enddo
       do i=1,nbix
         if(s3.gt.0.)then
           if(abs(iokoll).gt.0)then
             zu(i)=zu(i)/s3*float(abs(iokoll)) !to have the same normalization than in MC
           elseif(ncy.gt.0)then
             zu(i)=zu(i)/s3     !to have the same normalization than in MC (per Pom)
           endif
         endif
         write(ifhi,'(2e11.3)')xl(kk,i),zu(i)
       enddo
       if(ncy.gt.0)s3=1.    !by definition
       write(ifhi,'(a)')    '  endarray'
       write(ifhi,'(a,f7.3,a)')'text .1 .05 "(',s3,')"'
       write(ifhi,'(a)')    'closehisto plot 0'

       write(ifhi,'(a)')       '!----------------------------------'
       write(ifhi,'(a)')       '!   x-?PE! distribution    '//modu
       write(ifhi,'(a)')       '!----------------------------------'

       write(ifhi,'(a,i1)')   'openhisto name xmPomSimuL'//modu(3:4),ncy
       write(ifhi,'(a)')   'htyp lru xmod'//modu//'ymod log'
       write(ifhi,'(a,2e11.3)')'xrange',xlu,xlo
       write(ifhi,'(a)')    'text 0 0 "xaxis x-?PE!"'
       write(ifhi,'(a)') 'text 0 0 "yaxis'//imod//'?Pom! / dx-?PE!"'
       write(ifhi,'(a,f7.3,a)')'text 0.1 0.85 "Con. av',zavco2,'"'
       write(ifhi,'(a)')       'histoweight -1' !,dble(hw)
       write(ifhi,'(a)')       ' array 2'
       s1=0
       do i=1,nbix
       u=xl(kk,i)
       z=wxm(ncy,kk,i)/dxl(kk,i)
       s1=s1+z*dxl(kk,i)
        write(ifhi,'(2e11.3)')u,z
       enddo
       write(ifhi,'(a)')    '  endarray'
       !-----------------
       !multiply with ff
       !-----------------
       write(ifhi,'(a,i1,$)')
     . 'calc *0 / PxNorm3L'//modu(3:4),ncy
       write(ifhi,'(a)')' ; '
       !---------------
       write(ifhi,'(a)')    'closehisto plot 0-'

       write(ifhi,'(a,i1)')'openhisto name xmPomUnitL'//modu(3:4),ncy
       write(ifhi,'(a)')   'htyp lba xmod'//modu//'ymod log'
       write(ifhi,'(a,2e11.3)')'xrange',xlu,xlo
       write(ifhi,'(a)')    'text 0 0 "xaxis x-?PE!"'
       write(ifhi,'(a)') 'text 0 0 "yaxis'//imod//'?Pom! / dx-"'
       write(ifhi,'(a,e22.14)')'histoweight ',dble(hww)
       write(ifhi,'(a)')       'array 2'
       s2=0
       do i=1,nbix
        u=xl(kk,i)
        zu(i)=0
        if(ihww.gt.0)then
          if(ncy.gt.0.and.lppb)then
            zu(i)=sngl(PomIncMExact(dble(u),bb))
          else
            if(mmmm.eq.1)zu(i)=sngl(PomIncMExact(dble(u),bb))
            if(mmmm.eq.2)zu(i)=sngl(PomIncMIExact(dble(u))/sigine*10)
            if(mmmm.eq.3)zu(i)=sngl(PomIncMIExact(dble(u))/sigine*10)
          endif
          s2=s2+dxl(kk,i)*zu(i)
        endif
       enddo
       do i=1,nbix
         if(ncy.gt.0.and.s2.gt.0.)then
           zu(i)=zu(i)/s2    !to have the same normalization than in MC (per Pom)
         endif
         write(ifhi,'(2e11.3)')xl(kk,i),zu(i)
       enddo
       if(ncy.gt.0)s2=1.    !by definition
       write(ifhi,'(a)')    '  endarray'
       write(ifhi,'(a,f7.3,a,f7.3,a)')
     *                       'text .1 .15 "I= ',s1,' (',s2,')"'
       
       write(ifhi,'(a)')    'closehisto plot 0-'

       write(ifhi,'(a,i1)')   'openhisto name xmPomSatL'//modu(3:4),ncy
       write(ifhi,'(a)')   'htyp lgi xmod'//modu//'ymod log'
       write(ifhi,'(a,2e11.3)')'xrange',xlu,xlo
       write(ifhi,'(a)')    'text 0 0 "xaxis x+?PE!"'
       write(ifhi,'(a)') 'text 0 0 "yaxis'//imod//'?Pom! / dx+?PE!"'
       write(ifhi,'(a,e22.14)')'histoweight ',dble(hww)
       write(ifhi,'(a)')       'array 2'
       s3=0
       do i=1,nbix
        u=xl(kk,i)
        zu(i)=0
        if(ihww.gt.0)then
          if(ncy.gt.0.and.lppb)then
      zu(i)=sngl(PomInc(PomIncMExact,dble(u),bb,zavco2,zavco1,zavpom,2))
          else
            if(mmmm.eq.1)zu(i)=
     .    sngl(PomInc(PomIncMExact,dble(u),bb,zavco2,zavco1,zavpom,2))
            if(mmmm.eq.2)zu(i)=10./sigine
     .   *sngl(PomInc(PomIncMIExact,dble(u),0.,zavco2,zavco1,zavpom,-2))
            if(mmmm.eq.3)zu(i)=10./sigine
     .   *sngl(PomInc(PomIncMIExact,dble(u),0.,zavco2,zavco1,zavpom,-2))
          endif
          s3=s3+dxl(kk,i)*zu(i)
         endif
       enddo
       do i=1,nbix
         if(s3.gt.0.)then
           if(abs(iokoll).gt.0)then
             zu(i)=zu(i)/s3*float(abs(iokoll)) !to have the same normalization than in MC
           elseif(ncy.gt.0)then
             zu(i)=zu(i)/s3     !to have the same normalization than in MC (per Pom)
           endif
         endif
         write(ifhi,'(2e11.3)')xl(kk,i),zu(i)
       enddo
       if(ncy.gt.0)s3=1.    !by definition
       write(ifhi,'(a)')    '  endarray'
       write(ifhi,'(a,f7.3,a)')'text .1 .05 "(',s3,')"'
       write(ifhi,'(a)')    'closehisto plot 0'

  !................................................................
       xm=-1. !xm integration
       txtxm='xm int'
       if(mode.eq.2) goto 111
       do jjb=0,3
       b=jjb*0.5
       do jj=0,2

       write(ifhi,'(a)')       '!----------------------------------'
       write(ifhi,'(a,3i1)')   '!   ffom11    '//modu,jjb,jj
       write(ifhi,'(a)')       '!----------------------------------'

       write(ifhi,'(a,3i1)')'openhisto name ffom11L'//modu(3:4),jjb,jj+8
     .                                                        ,ncy
       write(ifhi,'(a)')    'htyp lin xmod'//modu//'ymod log'
       write(ifhi,'(a,2e11.3)')'xrange ',xlu,xlo
       write(ifhi,'(a)')'txt "xaxis  x+?PE!"'
       write(ifhi,'(a)')'txt "yaxis dn?Pom! / dx+?PE! "'
       write(ifhi,'(a)')'text 0.05 0.1  "fit and exact, all contrib."'
       if(jjb.lt.3)write(ifhi,'(a,f4.1,3a)')
     *             'txt "title ffom11   b =',b,'   ',txtxm,'"'
       if(jjb.ge.3)write(ifhi,'(3a)')
     *             'txt "title ffom11   b over   ',txtxm,'"'
       write(ifhi,'(a,e22.14)')'histoweight ',dble(hww)
       write(ifhi,'(a)')       'array 2'
       do i=1,nbix
       u=xl(kk,i)
       z=0
       if(ihww.gt.0)then
       if(jjb.lt.3.and.jj.eq.0)z= ffom11(u,xm,b,-1,-1)
       if(jjb.lt.3.and.jj.eq.1)z= ffom11(u,xm,b,0,5)
       if(jjb.lt.3.and.jj.eq.2)z= ffom11(u,xm,b,0,4)
       if(jjb.eq.3.and.jj.eq.0)z=ffom11a(u,xm,-1,-1)
       if(jjb.eq.3.and.jj.eq.1)z=ffom11a(u,xm,0,5)
       if(jjb.eq.3.and.jj.eq.2)z=ffom11a(u,xm,0,4)
       endif
       write(ifhi,'(2e11.3)')u,z
       enddo
       write(ifhi,'(a)')    '  endarray'
       if(jj.le.1)write(ifhi,'(a)')    'closehisto plot 0-'
       if(jj.eq.2)write(ifhi,'(a)')    'closehisto plot 0'

       enddo
       enddo

 111   continue
       do jjb=0,3
       b=jjb*0.5
       do jjj=1,6
       jj=jjj
       if(jjj.eq.6)jj=0

       write(ifhi,'(a)')       '!----------------------------------'
       write(ifhi,'(a,3i1)')   '!   ffom11    '//modu,jjb,jj
       write(ifhi,'(a)')       '!----------------------------------'

       do itt=1,1,-1           !test om51p tabulation in ffom11
         if(itt.eq.1)then
           isgn=1                   !tables
c           print *,'ici',jjb,b
         else
           isgn=-1                  !function (slow)
         endif

       write(ifhi,'(a,4i1)')'openhisto name om1ffL'//modu(3:4)
     .  ,jjb,jj,ncy
       if(itt.eq.1)then
       if(jj.eq.1)write(ifhi,'(a)')    'htyp lbu xmod'//modu//'ymod log'
       if(jj.eq.2)write(ifhi,'(a)')    'htyp lru xmod'//modu//'ymod log'
       if(jj.eq.3)write(ifhi,'(a)')    'htyp lyu xmod'//modu//'ymod log'
       if(jj.eq.4)write(ifhi,'(a)')    'htyp lgu xmod'//modu//'ymod log'
       if(jj.eq.5)write(ifhi,'(a)')    'htyp lfa xmod'//modu//'ymod log'
       if(jj.eq.0)write(ifhi,'(a)')    'htyp lro xmod'//modu//'ymod log'
       else
       if(jj.eq.1)write(ifhi,'(a)')    'htyp pbc xmod'//modu//'ymod log'
       if(jj.eq.2)write(ifhi,'(a)')    'htyp prc xmod'//modu//'ymod log'
       if(jj.eq.3)write(ifhi,'(a)')    'htyp pyc xmod'//modu//'ymod log'
       if(jj.eq.4)write(ifhi,'(a)')    'htyp pgc xmod'//modu//'ymod log'
       if(jj.eq.5)write(ifhi,'(a)')    'htyp pfq xmod'//modu//'ymod log'
       if(jj.eq.0)write(ifhi,'(a)')    'htyp prs xmod'//modu//'ymod log'
       endif
       write(ifhi,'(a,2e11.3)')'xrange ',xlu,xlo
       write(ifhi,'(a)')'yrange 1.e-6 auto'
       if(jj.eq.1)then
       write(ifhi,'(a)') 'txt "xaxis  x+?PE!"'
       write(ifhi,'(a)') 'txt "yaxis  dn?Pom! / dx+?PE!  "'
       if(kk.eq.2)then
        write(ifhi,'(a)') 'text 0.1 0.2  "soft sea-sea"'
        write(ifhi,'(a)') 'text 0.1 0.1  "val-sea sea-val val-val"'
       else
        write(ifhi,'(a)') 'text 0.05 0.8  "soft"'
        write(ifhi,'(a)') 'text 0.05 0.7  "diff"'
        write(ifhi,'(a)') 'text 0.05 0.6  "sea-sea"'
        write(ifhi,'(a)') 'text 0.05 0.5  "val-sea"'
        write(ifhi,'(a)') 'text 0.05 0.4  "sea-val"'
        write(ifhi,'(a)') 'text 0.05 0.3  "val-val"'
       endif
       if(jjb.lt.3)write(ifhi,'(a,f4.1,3a)')
     *             'txt "title ffom11   b =',b,'  ',txtxm,'"'
       if(jjb.ge.3)write(ifhi,'(3a)')
     *             'txt "title ffom11   b over  ',txtxm,'"'
       endif
       write(ifhi,'(a,e22.14)')'histoweight ',dble(hww)
       write(ifhi,'(a)')       'array 2'
c       print *,'ici'
       do i=1,nbix
       u=xl(kk,i)
       z=0
       if(ihww.gt.0)then
       if(jjb.lt.3)z= ffom11(u,xm,b,isgn*jj,isgn*jj)
       if(jjb.eq.3)z=ffom11a(u,xm,isgn*jj,isgn*jj)
       endif
       write(ifhi,'(2e11.3)')u,z
       enddo
       write(ifhi,'(a)')    '  endarray'
       if(jjj.ne.6.or.itt.eq.2)then
         write(ifhi,'(a)')    'closehisto plot 0-'
       else
         write(ifhi,'(a)')    'closehisto plot 0'
       endif

       enddo
       enddo
       enddo

 777   continue

      enddo

      enddo ! <-------------ncy----------------

      koll=kollini

      else

          write(ifch,*)'\n\n  ERROR 27092012b \n\n'
                   stop'\n\n  ERROR 27092012b \n\n'
      endif

      return
      end

c-----------------------------------------------------------------------
      subroutine emsaaFin  !calls xEmsP2fill
c-----------------------------------------------------------------------
      call xEmsPxfill
      call xEmsP2fill
      call xEmsP3fill
      end

c-----------------------------------------------------------------------
      subroutine emsaaDeb  !some counting even if no event
c-----------------------------------------------------------------------
      parameter(mmxcentrx=9)
      common/cncoll/wncoll(3,0:mmxcentrx+1),wnevt(3,0:mmxcentrx+1)
     .             ,wbim(0:mmxcentrx+1),wpom(0:7,0:mmxcentrx+1)
#include "aaa.h"
      common/nucl3/phi,bimp

      if(jzmode(3).ne.1.or.maproj+matarg.gt.2)return !only for bimevt (pp)
      call getMcentr   ! -> mcentrf() -> mmxcentrf()
      ncyi=mcentrf()
      if(ncyi.lt.0.or.ncyi.gt.mmxcentrf())stop'\n\n ERROR 27092012d\n\n'
      if(ncyi.eq.0)ncyi=mmxcentrx+1

c      print *,bimp,ncyi,ncol
      do no=1,1 !2  !not for MB
      if(no.eq.1)then
        ncy=ncyi
      else
        ncy=0
      endif
      wnevt(1,ncy)=wnevt(1,ncy)+1
      wnevt(2,ncy)=wnevt(2,ncy)+1
      wbim(ncy)=wbim(ncy)+bimp

      enddo

      end

c-----------------------------------------------------------------------
      subroutine xEmsPxfill
c-----------------------------------------------------------------------
      parameter(nbix=100,nbiy=50,mmxcentrx=9)
      common/cxw/x(2,nbix),dx(2,nbix),wxmc(0:mmxcentrx+1,2,nbix)
     . ,wxmcQp(0:mmxcentrx+1,2,nbix),xl(2,nbix),dxl(2,nbix)
     . ,wxp(0:mmxcentrx+1,2,nbix),wxm(0:mmxcentrx+1,2,nbix)
     . ,wxpQp(0:mmxcentrx+1,2,nbix),wxmQp(0:mmxcentrx+1,2,nbix)
     . ,wxpQt(0:mmxcentrx+1,2,nbix),wxmQt(0:mmxcentrx+1,2,nbix)
     . ,wxmcQt(0:mmxcentrx+1,2,nbix),wwmc(0:mmxcentrx+1,2,nbix)
     . ,womc(0:mmxcentrx+1,2,nbix)
      common/cyw/y(nbiy),wymc(0:mmxcentrx+1,nbiy)
     .,wymcQt(0:mmxcentrx+1,2,nbiy),wymcQp(0:mmxcentrx+1,2,nbiy),nyp,nym
      common/cemspxwxx/
     *   wxmcxx(2,nbix),wxmcQpxx(2,nbix),wxmcQtxx(2,nbix),wwmcxx(2,nbix)
     *  ,wxpxx(2,nbix),wxmxx(2,nbix),wxpQpxx(2,nbix)
     *  ,wxmQpxx(2,nbix),wxpQtxx(2,nbix),wxmQtxx(2,nbix)
     *  ,wymcxx(nbiy),wymcQtxx(2,nbiy),wymcQpxx(2,nbiy),nypxx,nymxx
     *  ,womcxx(2,nbix)
      common/cncoll/wncoll(3,0:mmxcentrx+1),wnevt(3,0:mmxcentrx+1)
     .             ,wbim(0:mmxcentrx+1),wpom(0:7,0:mmxcentrx+1)
      common/cwq2/wq2(3,0:mmxcentrx+1,3)
      common/cemspxq2/q2kkk1(2),nq2kkk1
      logical lppb
#include "aaa.h"
#include "ems.h"
      ncyi=mcentrf()
      if(ncyi.lt.0.or.ncyi.gt.mmxcentrf())stop'\n\n ERROR 27092012d\n\n'
      if(ncyi.eq.0)ncyi=mmxcentrx+1

      nko=kolevt !nglevt   !to compare MC to theo, kolevt matters
      lppb=.false.
      xpoma=0.
      xpome=0.
      con1=0.
      con2=0.
      nconn=0
      conn1=0
      conn2=0
      if(maproj.eq.1.and.matarg.eq.1)then
        nko=ikoevt              !Npom (all)
        k=1
        xpoma=nprt(k) !all poms
        xpome=npr(3,k) !xnpo(k)!=npr(4,k)          !high x Pom
        con1=xnpp(k)
        con2=xnpt(k)
        if(jzmode(3).eq.1)then  !trigger on b
          lppb=.true.
        endif
      else
        do k=1,koll
          if(itpr(k).ne.0)then
            xpoma=xpoma+nprt(k)
            xpome=xpome+npr(3,k) !xnpo(k)!=npr(4,k)
            con1=con1+xnpp(k)
            con2=con2+xnpt(k)
          endif
        enddo
        if(nko.gt.0)then
          xpoma=xpoma/float(nko)
          xpome=xpome/float(nko)
          con1=con1/float(nko)
          con2=con2/float(nko)
        endif
      endif
      do k=1,koll
        ip=iproj(k)
        it=itarg(k)
        if(npr(3,k).gt.0)then
          nconn=nconn+1
          do kxy=1,koll           !ckw21
            if(iproj(kxy).eq.ip)conn1=conn1+npr(3,kxy)
            if(itarg(kxy).eq.it)conn2=conn2+npr(3,kxy)
          enddo 
        endif
      enddo
      if(nconn.gt.0.)conn1=conn1/nconn
      if(nconn.gt.0.)conn2=conn2/nconn

      if(abs(iokoll).gt.0)nko=abs(iokoll)

      do no=1,2
      if(no.eq.1)then
        ff=1
        if(nko.ne.0.and..not.lppb)then     !not with b trigger to get absolute value (not per Pom)
         ff=1.  ! /nko done later!!!
        endif
        ncy=ncyi
      else
        ff=1
        ncy=0
      endif

      !write(ifmt,*)'nglevt = ',nglevt,'   ff = ',ff,'  (xEmsPx)'
      do k=1,2
      do i=1,nbix
        wxmc(ncy,k,i) = wxmc (ncy,k,i)+wxmcxx(k,i) *ff
        womc(ncy,k,i) = womc (ncy,k,i)+womcxx(k,i) *ff
        wxmcQp(ncy,k,i)= wxmcQp(ncy,k,i)+wxmcQpxx(k,i)
        wxmcQt(ncy,k,i)= wxmcQt(ncy,k,i)+wxmcQtxx(k,i)
        wwmc(ncy,k,i) = wwmc(ncy,k,i) +wwmcxx(k,i)*ff
        wxp(ncy,k,i)  = wxp (ncy,k,i) +wxpxx(k,i) *ff
        wxpQp(ncy,k,i) = wxpQp(ncy,k,i) +wxpQpxx(k,i)
        wxpQt(ncy,k,i) = wxpQt(ncy,k,i) +wxpQtxx(k,i)
        wxm(ncy,k,i)  = wxm(ncy,k,i)  +wxmxx(k,i) *ff
        wxmQp(ncy,k,i) = wxmQp(ncy,k,i) +wxmQpxx(k,i)
        wxmQt(ncy,k,i) = wxmQt(ncy,k,i) +wxmQtxx(k,i)
      enddo
      enddo
      do i=1,nbiy
        wymc(ncy,i) =wymc(ncy,i)+wymcxx(i) *ff
        wymcQp(ncy,1,i)=wymcQp(ncy,1,i)+wymcQpxx(1,i)
        wymcQt(ncy,1,i)=wymcQt(ncy,1,i)+wymcQtxx(1,i)
        wymcQp(ncy,2,i)=wymcQp(ncy,2,i)+wymcQpxx(2,i)
        wymcQt(ncy,2,i)=wymcQt(ncy,2,i)+wymcQtxx(2,i)
      enddo
      nyp=nyp+nypxx
      nym=nym+nymxx
      wncoll(1,ncy)=wncoll(1,ncy)+nko
      nrpom=igetNpom()
      wpom(0,ncy)=wpom(0,ncy)+1.
      wpom(1,ncy)=wpom(1,ncy)+con1
      wpom(2,ncy)=wpom(2,ncy)+con2
      wpom(3,ncy)=wpom(3,ncy)+xpome
      wpom(4,ncy)=wpom(4,ncy)+xpoma
      wpom(5,ncy)=wpom(5,ncy)+conn1
      wpom(6,ncy)=wpom(6,ncy)+conn2
      wpom(7,ncy)=wpom(7,ncy)+nrpom
      if(no.eq.1)then 
        call eventvariset(10,conn1) !z10zevt
        call eventvariset(11,conn2) !z11zevt
      endif
      if(ncy.eq.0.or..not.lppb)then !for bimevt, nglevt, filled in emsaaDeb
      wnevt(1,ncy)=wnevt(1,ncy)+1
      wbim(ncy)=wbim(ncy)+bimevt
      endif
      wq2(1,ncy,1)=wq2(1,ncy,1)+q2kkk1(1)
      wq2(1,ncy,2)=wq2(1,ncy,2)+q2kkk1(2)
      wq2(1,ncy,3)=wq2(1,ncy,3)+nq2kkk1

      enddo

      end

c-----------------------------------------------------------------------
      subroutine xEmsP2fill
c-----------------------------------------------------------------------
c 
c-----------------------------------------------------------------------
      parameter(nbixp=25,nbixm=25,nbipt=77,nbish=200,mmxcentrx=9)
      parameter(nbigb=25)
      common/cxb/xlp(2,nbixp),dxlp(2,nbixp),gbval(nbigb)
     *          ,xlm(2,nbixm),dxlm(2,nbixm)
     *          ,wxb(0:mmxcentrx+1,2,0:4,4,nbixp,nbixm)
     *          ,wxe(0:mmxcentrx+1,2,0:4,4,nbixp,nbixm)
     *          ,wxi(0:mmxcentrx+1,2,0:4,4,nbixp,nbixm)
     *          ,wxs(0:mmxcentrx+1,2,0:4,4,nbixp,nbixm)
      common/cxb3/wgb(0:mmxcentrx+1,2,0:4,4,nbigb)
      common/cxb4/wptob(0:mmxcentrx+1,0:4,4,nbipt)
      common/cptb/ptu,pto,ptob(nbipt),rapu,rapo
      common/cncoll/wncoll(3,0:mmxcentrx+1),wnevt(3,0:mmxcentrx+1)
     .             ,wbim(0:mmxcentrx+1),wpom(0:7,0:mmxcentrx+1)
      common/cemsp2q2/q2kkk2(2),nq2kkk2
      common/cwq2/wq2(3,0:mmxcentrx+1,3)
      common/cemswxx/
     *           wxbxx(2,0:4,4,nbixp,nbixm)
     *          ,wxexx(2,0:4,4,nbixp,nbixm)
     *          ,wxixx(2,0:4,4,nbixp,nbixm)
     *          ,wxsxx(2,0:4,4,nbixp,nbixm)
     *        ,wptobxx(0:4,4,nbipt)
     *        ,wgbxx(2,0:4,4,nbigb)
      logical lppb
#include "aaa.h"
      ncyi=mcentrf()
      if(ncyi.lt.0.or.ncyi.gt.mmxcentrf())stop'\n\n ERROR 27092012g\n\n'
      if(ncyi.eq.0)ncyi=mmxcentrx+1

      nko=nglevt
      lppb=.false.
      if(maproj.eq.1.and.matarg.eq.1)then
        nko=ikhevt
        if(jzmode(3).eq.1)lppb=.true.
      endif

      do no=1,2

      if(no.eq.1)then
        ff=1
        if(nko.ne.0.and..not.lppb)then
         ff=1./nko
        endif
        ncy=ncyi
      else
        ff=1
        ncy=0
      endif

c      write(ifmt,*)ncy,'nglevt = ',nglevt,nko,'   ff = ',ff,'  (xEmsP2)'
      do j=1,nbixm
      do i=1,nbixp
      do jexi=1,4
      do jaai=0,4
      do ij=1,2
        wxb(ncy,ij,jaai,jexi,i,j)=wxb(ncy,ij,jaai,jexi,i,j)
     .                     +wxbxx(ij,jaai,jexi,i,j)*ff
        wxe(ncy,ij,jaai,jexi,i,j)=wxe(ncy,ij,jaai,jexi,i,j)
     .                     +wxexx(ij,jaai,jexi,i,j)*ff
        wxi(ncy,ij,jaai,jexi,i,j)=wxi(ncy,ij,jaai,jexi,i,j)
     .                     +wxixx(ij,jaai,jexi,i,j)*ff
        wxs(ncy,ij,jaai,jexi,i,j)=wxs(ncy,ij,jaai,jexi,i,j)
     .                     +wxsxx(ij,jaai,jexi,i,j)*ff
      enddo
      enddo
      enddo
      enddo
      enddo
      do i=1,nbipt
      do jexi=1,4
      do jaai=0,4
      wptob(ncy,jaai,jexi,i)=wptob(ncy,jaai,jexi,i)
     .                +wptobxx(jaai,jexi,i)*ff
      enddo
      enddo
      enddo
      do i=1,nbigb
      do jexi=1,4
      do jaai=0,4
      do ij=1,2
      wgb(ncy,ij,jaai,jexi,i)=wgb(ncy,ij,jaai,jexi,i)
     .                +wgbxx(ij,jaai,jexi,i)*ff
      enddo
      enddo
      enddo
      enddo

c      nko=nglevt
c      if(maproj.eq.1.and.matarg.eq.1)nko=ikhevt
      wncoll(2,ncy)=wncoll(2,ncy)+nko
      if(ncy.eq.0.or..not.lppb) !for bimevt filled in emsaaDeb
     .wnevt(2,ncy)=wnevt(2,ncy)+1
      wq2(2,ncy,1)=wq2(2,ncy,1)+q2kkk2(1)
      wq2(2,ncy,2)=wq2(2,ncy,2)+q2kkk2(2)
      wq2(2,ncy,3)=wq2(2,ncy,3)+nq2kkk2
      !print*,'xEmsPcheckB',ncy,nq2kkk2,q2kkk2,wq2(2,ncy,2),wq2(2,ncy,1)

      enddo

      end

c-----------------------------------------------------------------------
      subroutine xxEmsP2(iii,jai,jex,jij,xpd,xmd,xpb,xmb,pt1,pt2
     . ,rap1,rap2,gb1,gb2,q2kk)
c-----------------------------------------------------------------------
c       xxEmsP2(): formerly in xEmsP2(iii=1)
c-----------------------------------------------------------------------
c plot  x+ distributions of Pomeron ends (PE) (xpd)
c          and Pomeron's in Born (IB) partons (xpb),
c     and pt dist of Pomeron's out Born (OB) partons
c       integrated over x- bins (xmd,xmb)
c     and F2 taken from pp collision
c  iii=1: fill arrays
c  jai: type of semihard Pomeron
c         0= sea-sea no evolution,
c         1= sea-sea, 2= val=sea, 3= sea-val, 4= val-val
c         5= all hard  for iii=2
c  jex: emission type
c         1= no emission, 2= proj emis, 3= targ emis, 4= both sides
c         5= all for iii=2
c-----------------------------------------------------------------------

#include "aaa.h"
#include "sem.h"
#include "ems.h"
      common/geom/rmproj,rmtarg,bmax,bkmx
      parameter(nbixp=25,nbixm=25,nbipt=77,nbish=200,mmxcentrx=9)
      parameter(nbigb=25)
      common/cptb/ptu,pto,ptob(nbipt),rapu,rapo
      common/cemswxx/
     *           wxbxx(2,0:4,4,nbixp,nbixm)
     *          ,wxexx(2,0:4,4,nbixp,nbixm)
     *          ,wxixx(2,0:4,4,nbixp,nbixm)
     *          ,wxsxx(2,0:4,4,nbixp,nbixm)
     *        ,wptobxx(0:4,4,nbipt)
     *        ,wgbxx(2,0:4,4,nbigb)
      common/cemspbx/xlub1,xlub2,xlob,gbug,gbog
      common /ar3/   x1(7),a1(7)
      common/shat_common/shat
      common/cemsp2q2/q2kkk2(2),nq2kkk2
      dimension q2kk(2)
      real shat
      double precision psjgrv08x,psjcteq6x,psjvrx,psjwox
      external psjgrv08x,psjcteq6x,psjvrx,psjwox
c      sdummy=sh_min
c      sdummy=sh_max

      iG100=0
      
      if(iemspbx.eq.0)call utstop('ERROR in xEmsP2: iemspbx = 0&')
      
      if(iii.eq.1.and.jai.eq.-1)then

       do ij=1,2
       do i=1,nbixp
       do j=1,nbixm
       do jaai=0,4
       do jexi=1,4
        wxbxx(ij,jaai,jexi,i,j)=0.
        wxexx(ij,jaai,jexi,i,j)=0.
        wxixx(ij,jaai,jexi,i,j)=0.
        wxsxx(ij,jaai,jexi,i,j)=0.
       enddo
       enddo
       enddo
       enddo
       enddo
       do i=1,nbipt
       do jaai=0,4
       do jexi=1,4
       wptobxx(jaai,jexi,i)=0
       enddo
       enddo
       enddo
       do ij=1,2
       do i=1,nbigb
       do jaai=0,4
       do jexi=1,4
       wgbxx(ij,jaai,jexi,i)=0
       enddo
       enddo
       enddo
       enddo
       q2kkk2(1)=0
       q2kkk2(2)=0
       nq2kkk2=0


      elseif(iii.eq.1)then

       jaa=jai
       if(jaa.eq.1.or.jaa.gt.10)then !select subprocess in sea-sea contribution
         if(iseahot.ge.0.and.jaa.ne.1+iseahot*10)then
           return
c         else
c           print *,'ici',jaa,1+iseahot*10
         endif
         jaa=1
       endif
        
       xp=xpb
       xm=xmb
       if(xp.lt.xlub1)goto 2
       if(xm.lt.xlub1)goto 2
       i=1+int(alog(xp/xlub1)/alog(xlob/xlub1)*nbixp)
       if(i.gt.nbixp)goto 2
       if(i.lt.1)goto 2
       j=1+int(alog(xm/xlub1)/alog(xlob/xlub1)*nbixm)
       if(j.gt.nbixm)goto 2
       if(j.lt.1)goto 2
       wxbxx(1,jaa,jex,i,j)=wxbxx(1,jaa,jex,i,j)+1.
2      continue

       if(xp.lt.xlub2)goto 12
       if(xm.lt.xlub2)goto 12
       i=1+int((xp-xlub2)/(xlob-xlub2)*nbixp)
       if(i.gt.nbixp)goto 12
       if(i.lt.1)goto 12
       j=1+int((xm-xlub2)/(xlob-xlub2)*nbixm)
       if(j.gt.nbixm)goto 12
       if(j.lt.1)goto 12
       wxbxx(2,jaa,jex,i,j)=wxbxx(2,jaa,jex,i,j)+1.
12     continue

       xp=xpd
       xm=xmd
       if(jij.gt.1)goto 32 !avoid double counting

       if(xp.lt.xlub1)goto 22
       if(xm.lt.xlub1)goto 22
       i=1+int(alog(xp/xlub1)/alog(xlob/xlub1)*nbixp)
       if(i.gt.nbixp)goto 22
       if(i.lt.1)goto 22
       j=1+int(alog(xm/xlub1)/alog(xlob/xlub1)*nbixm)
       if(j.gt.nbixm)goto 22
       if(j.lt.1)goto 22
       wxexx(1,jaa,jex,i,j)=wxexx(1,jaa,jex,i,j)+1.
       if(pt1.gt.100.or.pt2.gt.100)
     . wxixx(1,jaa,jex,i,j)=wxixx(1,jaa,jex,i,j)+1.
  22   continue

       if(xp.lt.xlub2)goto 32
       if(xm.lt.xlub2)goto 32
       i=1+int((xp-xlub2)/(xlob-xlub2)*nbixp)
       if(i.gt.nbixp)goto 32
       if(i.lt.1)goto 32
       j=1+int((xm-xlub2)/(xlob-xlub2)*nbixm)
       if(j.gt.nbixm)goto 32
       if(j.lt.1)goto 32
       wxexx(2,jaa,jex,i,j)=wxexx(2,jaa,jex,i,j)+1.
       if(pt1.gt.100.or.pt2.gt.100)
     . wxixx(2,jaa,jex,i,j)=wxixx(2,jaa,jex,i,j)+1.
  32   continue

       do m=1,2
       if(m.eq.1)then
         pt=pt1
         rap=rap1
       elseif(m.eq.2)then
         pt=pt2
         rap=rap2
       endif
       if(rap.lt.rapu)goto 42
       if(rap.gt.rapo)goto 42
       if(pt.lt.ptu)goto 42 !KW1908: this command is needed !!!
       delptLOG=(log(pto)-log(ptu))/nbipt  !KW1908: log binning needed at low pt !!!
       i=1+int( (log(pt)-log(ptu)) / delptLOG )
       !i=1+int((pt-ptu)/(pto-ptu)*nbipt)
       if(i.lt.1)goto 42
       if(i.gt.nbipt)goto 42
       wptobxx(jaa,jex,i)=wptobxx(jaa,jex,i)+1.
   42  continue
       enddo

       if(gb1.le.0.)goto 3
       gb=log10(gb1)
       if(gb.lt.gbug)goto 3
       i=1+int((gb-gbug)/(gbog-gbug)*nbigb)
       if(i.gt.nbigb)goto 3
       if(i.lt.1)goto 3
       wgbxx(1,jaa,jex,i)=wgbxx(1,jaa,jex,i)+1.
3      continue

       if(gb2.le.0.)goto 4
       gb=log10(gb2)
       if(gb.lt.gbug)goto 4
       i=1+int((gb-gbug)/(gbog-gbug)*nbigb)
       if(i.gt.nbigb)goto 4
       if(i.lt.1)goto 4
       wgbxx(2,jaa,jex,i)=wgbxx(2,jaa,jex,i)+1.
4      continue

       xp=xpd !??? not very  
       xm=xpb !??? useful

       if(xp.lt.xlub1)goto 52
       if(xm.lt.xlub1)goto 52
       i=1+int(alog(xp/xlub1)/alog(xlob/xlub1)*nbixp)
       if(i.gt.nbixp)goto 52
       if(i.lt.1)goto 52
       j=1+int(alog(xm/xlub1)/alog(xlob/xlub1)*nbixm)
       if(j.gt.nbixm)goto 52
       if(j.lt.1)goto 52
       wxsxx(1,jaa,jex,i,j)=wxsxx(1,jaa,jex,i,j)+1.
  52   continue

       if(xp.lt.xlub2)goto 62
       if(xm.lt.xlub2)goto 62
       i=1+int((xp-xlub2)/(xlob-xlub2)*nbixp)
       if(i.gt.nbixp)goto 62
       if(i.lt.1)goto 62
       j=1+int((xm-xlub2)/(xlob-xlub2)*nbixm)
       if(j.gt.nbixm)goto 62
       if(j.lt.1)goto 62
       wxsxx(2,jaa,jex,i,j)=wxsxx(2,jaa,jex,i,j)+1.
  62   continue

c       if(jaa.gt.0)then
         q2kkk2(1)=q2kkk2(1)+q2kk(1)
         q2kkk2(2)=q2kkk2(2)+q2kk(2)
         nq2kkk2=nq2kkk2+1
c       endif
       !print*,'xEmsPcheckA',nq2kkk2,q2kkk2

      else

          write(ifch,*)'\n\n  ERROR 27092012e \n\n'
                   stop'\n\n  ERROR 27092012e \n\n'
      endif

      return
      end

c-----------------------------------------------------------------------
      subroutine xEmsP2(iii,jaa,jex,xpd,xmd,xpb,xmb,pt1,pt2,rap1,rap2)
c-----------------------------------------------------------------------
c plot  x+ distributions of Pomeron ends (PE) (xpd)
c          and Pomeron's in Born (IB) partons (xpb),
c     and pt dist of Pomeron's out Born (OB) partons
c       integrated over x- bins (xmd,xmb)
c  iii=0: initialize
c  iii>=2: make histogram
c           (2 - Pomeron end PE, 3 - in Born IB, 4 - out Born OB - 5 sh
c            6 - GBmin         , 7 - GBmax     , 8 - xp_Born vs xp_PE)
c  jaa: type of semihard Pomeron
c         0= soft+gluon,
c         1= sea-sea, 2= val=sea, 3= sea-val, 4= val-val
c         5= all hard for iii=2
c  jex: emission type
c         1= no emission, 2= proj emis, 3= targ emis, 4= both sides
c         5= all  for iii=2
c-----------------------------------------------------------------------
c    xpar2 .... choices of (up to three) models
c                  1 = ours + MC
c                  2 = ours + CTEQ + MC
c                  3 = ours + CTEQ
c                  4 = ours(charm) + CTEQ(charm)
c                 41 = CTEQ(charm)
c                  5 = ours(charm)
c                  6 = ours(bottom) + CTEQ(bottom)
c    xpar3 .... 1=log 2=lin
c    xpar4 .... xmax                             for iii eq 9
c    xpar5 .... ihq=nint(xpar5)                  for iii ne 4
c    xpar6 .... xmax (not for iii=4)
c    xpar7 .... 1=MB
c    xpar8 .... xmin for theo function (if 1: adds 'plot 0')
c    xpar9 .... number of bins
c    xpar10 ... factor for theo fuction
c    xpar11 ... q2cmin(1 and 2)=xpar11
c    xpar12 ... pdfparamset(1 and 3, xpar12 )
c    xpar13 ... nlow (lowest bin considered)
c-----------------------------------------------------------------------

#include "aaa.h"
#include "sem.h"
#include "ems.h"
      common/geom/rmproj,rmtarg,bmax,bkmx
      parameter(nbixp=25,nbixm=25,nbipt=77,nbish=200,mmxcentrx=9)
      parameter(nbigb=25)
      common/cxb/xlp(2,nbixp),dxlp(2,nbixp),gbval(nbigb)
     *          ,xlm(2,nbixm),dxlm(2,nbixm)
     *          ,wxb(0:mmxcentrx+1,2,0:4,4,nbixp,nbixm)
     *          ,wxe(0:mmxcentrx+1,2,0:4,4,nbixp,nbixm)
     *          ,wxi(0:mmxcentrx+1,2,0:4,4,nbixp,nbixm)
     *          ,wxs(0:mmxcentrx+1,2,0:4,4,nbixp,nbixm)
      common/cxb3/wgb(0:mmxcentrx+1,2,0:4,4,nbigb)
      common/cxb4/wptob(0:mmxcentrx+1,0:4,4,nbipt)
      common/cptb/ptu,pto,ptob(nbipt),rapu,rapo
      common/cxb2/uxlp(2,nbixp),oxlp(2,nbixp)
     .           ,uxlm(2,nbixm),oxlm(2,nbixm)
      common/cncoll/wncoll(3,0:mmxcentrx+1),wnevt(3,0:mmxcentrx+1)
     .             ,wbim(0:mmxcentrx+1),wpom(0:7,0:mmxcentrx+1)
      common/cwq2/wq2(3,0:mmxcentrx+1,3)
      common/cemspbx/xlub1,xlub2,xlob,gbug,gbog
      common /ar3/   x1(7),a1(7)
      common/shat_common/shat
      integer ioTestFact
      common /ctestfact/ ioTestFact
      integer laddTestFact
      common /claddtestfact/ laddTestFact
      integer modeDeltaFzero
      common /cmodeDeltaFzero/ modeDeltaFzero
      integer iffsigiEfficiency
      common /ciffsigiEfficiency/ iffsigiEfficiency
         !double precision xCTEQ2,pifpartone!TEST-CTEQ-EPOS
      real shat
      character cii*1,dii*6,c3*3,modu*5
      integer ii,jj
      real xpmin0
      real ff(3)
      common/cicentrality/icentrality
      integer jjselect(10)
      real q2amin(2)
      double precision psjgrv08x,psjcteq6x,psjvrx,psjwox
      external psjgrv08x,psjcteq6x,psjvrx,psjwox
      data ncntEmsP21/0/ncntEmsP22/0/ncntEmsP23/0/ncntEmsP24/0/
     .ncntEmsP25/0/ncntEmsP26/0/ncntEmsP27/0/ncntEmsP28/0/ncntEmsP29/0/
c      save sh_min,sh_max
      save ncntEmsP21,ncntEmsP22,ncntEmsP23,ncntEmsP24
     .,ncntEmsP25,ncntEmsP26,ncntEmsP27,ncntEmsP28,ncntEmsP29
c      sdummy=sh_min
c      sdummy=sh_max
      dummy=xpd
      dummy=xmd
      dummy=xpb
      dummy=xmb
      dummy=pt1
      dummy=pt2
      dummy=rap1
      dummy=rap2

      iG100=0

      do jj=1,10
        jjselect(jj)=0
      enddo
      jjlow=1 

      ipar2=nint(xpar2)
      if(ipar2.eq.0.or.ipar2.eq.1)then
        jjup=1 
        jjselect(1)=1 !ours
        iSKIP_MC=0    !MC
      elseif(ipar2.eq.2)then
        jjup=2
        jjselect(1)=1 !ours    
        jjselect(2)=6 !CTEQ14  
        iSKIP_MC=0    !MC
      elseif(ipar2.eq.3)then
        jjup=2
        jjselect(1)=1 !ours    
        jjselect(2)=6 !CTEQ14  
        iSKIP_MC=1    !skip MC
      elseif(ipar2.eq.4)then
        jjup=2
        jjselect(1)=7 !ours charm    
        jjselect(2)=8 !CTEQ14 charm 
        iSKIP_MC=1    !skip MC
      elseif(ipar2.eq.41)then
        jjup=1
        jjselect(1)=8 !CTEQ14 charm 
        iSKIP_MC=1    !skip MC
      elseif(ipar2.eq.5)then
        jjup=1
        jjselect(1)=7 !ours charm    
        iSKIP_MC=1    !skip MC
      elseif(ipar2.eq.6)then
        jjup=2
        jjselect(1)=9 !ours bottom   
        jjselect(2)=10 !CTEQ14 bottom 
        iSKIP_MC=1    !skip MC
      else
        stop'ERROR 14042022'
      endif
     
      ifacpdf4=0
      if(iii.eq.4)then
        if(xpar12.gt.0.00001)then
          call pdfparamget(1,facpdf4)
          call pdfparamget(3,facpdf5)
          facpdf4_save=facpdf4
          facpdf5_save=facpdf5
          ifacpdf4=1
          call pdfparamset(1, xpar12 )
          call pdfparamset(3, xpar12 )
        endif 
      endif

      if(iemspbx.eq.0)call utstop('ERROR in xEmsP2: iemspbx = 0&')
      if(mmxcentrf().gt.mmxcentrx)stop'\n\n ERROR 27092012h \n\n '

      !---------------
      if(iii.eq.0)then
      !---------------

       !write(ifch,*)'initialize xEmsP2 tables for'
       !. ,mmxcentrf(),' centralities'
       xlub1=1./engy**2
       xlub2=0.01!/engy
       xlob=1.
       gbug=-6
       gbog=2
       do i=1,nbixp
        xlp(1,i)=xlub1*(xlob/xlub1)**((i-0.5)/nbixp)
        xlp(2,i)=xlub2+(xlob-xlub2)*((i-0.5)/nbixp)
        dxlp(1,i)=xlub1*(xlob/xlub1)**(1.*i/nbixp)
     *             *(1.-(xlob/xlub1)**(-1./nbixp))
        dxlp(2,i)=(xlob-xlub2)/nbixp
        uxlp(1,i)=xlub1*(xlob/xlub1)**((i-1.)/nbixp)
        uxlp(2,i)=xlub2+(xlob-xlub2)*((i-1.)/nbixp)
        oxlp(1,i)=xlub1*(xlob/xlub1)**((i-0.)/nbixp)
        oxlp(2,i)=xlub2+(xlob-xlub2)*((i-0.)/nbixp)
       enddo
       do i=1,nbixm
        xlm(1,i)=xlub1*(xlob/xlub1)**((i-0.5)/nbixm)
        xlm(2,i)=xlub2+(xlob-xlub2)*((i-0.5)/nbixm)
        dxlm(1,i)=xlub1*(xlob/xlub1)**(1.*i/nbixm)
     *             *(1.-(xlob/xlub1)**(-1./nbixm))
        dxlm(2,i)=(xlob-xlub2)/nbixm
        uxlm(1,i)=xlub1*(xlob/xlub1)**((i-1.)/nbixm)
        uxlm(2,i)=xlub2+(xlob-xlub2)*((i-1.)/nbixm)
        oxlm(1,i)=xlub1*(xlob/xlub1)**((i-0.)/nbixm)
        oxlm(2,i)=xlub2+(xlob-xlub2)*((i-0.)/nbixm)
       enddo
       do i=1,nbigb
        gbval(i)=gbug+(gbog-gbug)*((i-0.5)/nbigb)
       enddo
       do ij=1,2
       do i=1,nbixp
       do j=1,nbixm
       do jaai=0,4
       do jexi=1,4
       do n=0,mmxcentrf()
        wxb(n,ij,jaai,jexi,i,j)=0.
        wxe(n,ij,jaai,jexi,i,j)=0.
        wxi(n,ij,jaai,jexi,i,j)=0.
        wxs(n,ij,jaai,jexi,i,j)=0.
       enddo
       enddo
       enddo
       enddo
       enddo
       enddo
       ptu=2
       pto=2500.
       rapo= 0.5
       if(modeDeltaFzero.eq.1)then !Fzero~delta(x-1) 
         if(laddTestFact.eq.1)then !dn/dy in ffsig is delta function
           rapo=+100
         else
           rapo=+0.5
         endif
       endif
       rapu=-rapo 
       do i=1,nbipt
       ptob(i)=ptu*(pto/ptu)**((i-0.5)/nbipt) !log binning needed at low pt !!!
       !ptob(i)=ptu+(pto-ptu)*(i-0.5)/nbipt
       do jaai=0,4
       do jexi=1,4
       do n=0,mmxcentrf()
       wptob(n,jaai,jexi,i)=0
       enddo
       enddo
       enddo
       enddo
       do ij=1,2
       do i=1,nbigb
       do jaai=0,4
       do jexi=1,4
       do n=0,mmxcentrf()
       wgb(n,ij,jaai,jexi,i)=0
       enddo
       enddo
       enddo
       enddo
       enddo
       do n=0,mmxcentrf()
        wncoll(2,n)=0.
        wnevt(2,n)=0.
        do j=1,3
          wq2(2,n,j)=0
        enddo
       enddo

      !-------------------
      elseif(iii.ge.2)then
      !-------------------

       if(mmxcentrf().gt.0)then
         ncy1=1
       else
         ncy1=0
       endif
       ncy2=mmxcentrf()

       if(nint(xpar7).eq.1)then
         ncy1=0
         ncy2=0
       endif
       !zavcollmb--> ff4/ff5
       ff4=wncoll(2,0)
       ff5=wnevt(2,0)
       jff45=1

       do ncy=ncy1,ncy2

       nstat=nfull
       if(nstat.lt.0)nstat=max(1,nevent)
       hw=nstat
       zavcoll=0
       q2amin(1)=q2nmin
       q2amin(2)=q2nmin
       q2cmin(1)=q2nmin
       q2cmin(2)=q2nmin
       if(wnevt(2,ncy).gt.0..and.wq2(2,ncy,3).gt.0.)then
         zavcoll=wncoll(2,ncy)/wnevt(2,ncy)
c used in iii=2 and iii=5 only (q2cmin redefined in ffom12)
         q2amin(1)=wq2(2,ncy,1)/wq2(2,ncy,3) !nuclear component of pair averaged Q2s (no x)
         q2amin(2)=wq2(2,ncy,2)/wq2(2,ncy,3)
c used in iii=3 and iii=4 only
c        q2cmin updated later in fom function to take into account x evolution taken from the data
         if(iii.eq.4)then !for jet xs, better to use the reference line always
           q2cmin(1)=q2nmin
           q2cmin(2)=q2nmin
         else
           q2cmin(1)=max(q2nmin,q2amin(1))
           q2cmin(2)=max(q2nmin,q2amin(2))
         endif
         q2pmin(1)=q2cmin(1)
         q2pmin(2)=q2cmin(2)
       elseif(xpar11.gt.0.)then
         q2cmin(1)=xpar11
         q2cmin(2)=xpar11
       endif

       !cKW21 fix "add" problem
       !former ff is now ff(1)*ff(2)/ff(3), ff(i) put into histo
       ff(1)=1.
       ff(2)=1.
       ff(3)=1.
       if(maproj.eq.1.and.matarg.eq.1.and.bminim.eq.bmaxim)then
        ff(2)=float(nstat)
        ff(3)=float(ntevt)
        bb=bmaxim
       elseif(maproj.eq.1.and.matarg.eq.1)then
        if(ncy.gt.0.and.wnevt(2,ncy).gt.0.)then
          ff(2)=ff(2)*nstat
          ff(3)=ff(3)*wnevt(2,ncy)
        endif 
        if(ncy.gt.0.and.zavcoll.gt.0..and.jzmode(3).eq.1)then
          ff(2)=ff(2)*wnevt(2,ncy)    !in case of b trigger, no normalization by 1/npom in xEmsP2fill
          ff(3)=ff(3)*wncoll(2,ncy)
        endif
        if(ncy.gt.0)hw=wnevt(2,ncy)
        !q2pmin fixed in function direclty (averaged over b)
       else
         ! 1./nglevt is already included (in xEmsP2fill) for ncy.ne.0
         jff45=0 !no zavcollmb factor   !this is to correct by <Npom> in pp only
         if(ncy.gt.0)then
           if(wnevt(2,ncy).gt.0.)then
             ff(2)=ff(2)*nstat
             ff(3)=ff(3)*wnevt(2,ncy)
           endif 
           hw=wnevt(2,ncy)
         elseif(bminim.lt.0.001.and.bmaxim.gt.20)then !min bias
           area=pi*(rmproj+rmtarg)**2
           ff(1)=area/(maproj*matarg)/sigine*10 !Ncol=A*B included
           ff(2)=float(nstat)
           ff(3)=float(ntevt) 
         else !make plot nevertheless to avoid "histo not found"
           area=1
           ff(1)=1
           ff(2)=1
           ff(3)=1
         endif
       endif

c       print*,'xEmsPcheckC',ncy,q2amin,wq2(2,ncy,2),wq2(2,ncy,1),ff

       if(wnevt(2,ncy).lt.0.0001.or.wq2(2,ncy,3).lt.0.0001)then
         hww=0
         ihww=0
       else
         hww=1
         ihww=1
       endif
       
       ihwfac=1 
       if(iffsigiEfficiency.eq.1)then
         if(iopcnt.gt.1)then
           ihwfac=0
           if(iii.eq.4)hww=0       !overwrite weight not to use ffsig
         else        !first file only (batch mode)
           ihww=0
           if(iii.eq.4)then
             hww=1              !overwrite weight always to use ffsig
           else                 !other simu should not be used
             hww=0
             ihww=0
           endif
         endif
       endif

       j1=1        !first xminus bin
       j2=nbixm    !last xminus bin
       if(iii.eq.4)j2=1
       kkk=2
       if(nint(xpar3).gt.0)kkk=nint(xpar3)  !1=log 2=lin
       if(kkk.eq.1)then
         xlub=xlub1
       elseif(kkk.eq.2)then
         xlub=xlub2
       endif
       if(iii.eq.9)then
         if(kkk.eq.1)then
           j=1+int(log(xpar4/xlub1)/alog(xlob/xlub1)*nbixp)
         else
           j=1+int((xpar4-xlub2)/(xlob-xlub2)*nbixp)
         endif
         j1=min(max(1,j),nbixp)
         j2=j1
       endif

       jaa1=jaa
       jaa2=jaa
       jex1=jex
       jex2=jex
       if(jaa.eq.5)then
       jaa1=1
       jaa2=4
       endif
       if(jaa.eq.0)then
       jaa1=11 !5
       jaa2=11 !5
       endif
       if(jex.eq.5)then
       jex1=1
       jex2=4
       endif

       if(jex.eq.1)then
        je1=0
        je2=0
       elseif(jex.eq.2)then
        je1=1
        je2=0
       elseif(jex.eq.3)then
        je1=0
        je2=1
       elseif(jex.eq.4)then
        je1=1
        je2=1
       elseif(jex.eq.5)then
        je1=2
        je2=2
       endif

       if(iii.eq.2)cii='i'
       if(iii.eq.5)cii='j'
       if(iii.eq.2)dii='x+?PE!'
       if(iii.eq.5)dii='x?PE! '

c Initialize some needed variables for function F()
c       if(iscreen.ne.0)then
c         s=engy*engy
c         b=0.
c         do i=idxDmin(iomega),idxDmax(iomega)
c           zp=0.
c           zt=0.
c         call Gfunpar(zp,zt,0.,0.,1,i,b,s,alpx,betx,betpx,epsp,epst,epss,gamv)
c         enddo
c       endif

       if(kkk.eq.1)modu=' log '
       if(kkk.eq.2)modu=' lin '
       
       ihq    = 0               !all processes in born
       if(iii.eq.4.and.abs(xpar5).gt.0.)ihq=nint(xpar5)

       call ipoBasics(q2cmin)
       call ipoCSzeroTables(q2cmin)

       !################## Theoretical functions ##################

       if(iii.eq.2.or.iii.eq.5)then

        write(ifhi,'(a)')       '!----------------------------------'
        write(ifhi,'(a,3i1)')   '!   PE    ',jaa,jex
        write(ifhi,'(a)')       '!----------------------------------'

        sum=0
        if(ihww.gt.0.and.jaa1.ne.5)then
        sum=ffom12aii(max(1,jaa1),jaa2,je1,je2)
        endif
        write(ifhi,'(a,3i1)')'openhisto name ffom12a'//cii,jaa,jex,ncy
        write(ifhi,'(a)')'htyp lin xmod '//modu//' ymod log'
        write(ifhi,'(a,2e11.3)')'xrange ',xlub,xlob
        if(kkk.eq.1)write(ifhi,'(a)')'yrange 1e-5 auto'
        write(ifhi,'(a)')'txt "xaxis  '//dii//'"'
        write(ifhi,'(a)')'txt "yaxis dn?semi! / d'//dii//'    "'
        if(jaa1.eq.5.or.jaa1.eq.11.and.jex.eq.5)then
        write(ifhi,'(a,2i1,a)')
     .  'txt "title ffom11a'//cii//'  MC ('
     .  ,jaa,jex,')"'
        else
        write(ifhi,'(a,2i1,a)')
     .  'txt "title ffom12a'//cii//'  ffom11a'//cii//'  MC ('
     .  ,jaa,jex,')"'
        endif
        write(ifhi,'(a,e22.14)')'histoweight ',dble(hww)
        write(ifhi,'(a)')'array 2'
        do i=1,nbixp
         u=xlp(kkk,i)
         z=0
         if(ihww.gt.0..and.jaa1.ne.5)then
         if(iii.eq.2)z=ffom12ai(u,xlub,max(1,jaa1),jaa2,je1,je2)
         if(iii.eq.5)z=ffom12aj(u,max(1,jaa1),jaa2,je1,je2)
         endif
         write(ifhi,'(2e11.3)')u,z
        enddo
        write(ifhi,'(a)')    '  endarray'
        if(jaa1.eq.5.or.jaa1.eq.11.or.jex.eq.5)then
          call ipoOm5Tables(1)
          write(ifhi,'(a)')    'closehisto plot 0-'
          write(ifhi,'(a,3i1)')'openhisto name ffom11a'//cii,jaa,jex,ncy
          write(ifhi,'(a)')'htyp lba'
          write(ifhi,'(a,2i1,a)')'text 0.83 0.75 "(',jaa,jex,')"'
          write(ifhi,'(a,e22.14)')'histoweight ',dble(hww)
          write(ifhi,'(a)')'array 2'
          do i=1,nbixp
           u=xlp(kkk,i)
           z=0
           if(ihww.gt.0)then
           if(iii.eq.2)z=ffom11ai(u,xlub,jaa1,jaa2)
           if(iii.eq.5)z=ffom11aj(u,jaa1,jaa2)
           endif
           write(ifhi,'(2e11.3)')u,z
          enddo
          write(ifhi,'(a)')    '  endarray'
        endif

       elseif(iii.eq.3)then

        jaaa=jaa
        if(jaa.eq.0)jaaa=11
        write(ifhi,'(a)')       '!----------------------------------'
        write(ifhi,'(a,3i1)')   '!   IB    ',jaa,jex
        write(ifhi,'(a)')       '!----------------------------------'

    !.......total integral
        sum=0
        ei=1./engy
        xpmin0 = ei**2
        if(jaaa.eq.11)then
          s2min=4*q2sft
        else
          s2min=4*max(q2cmin(1),q2cmin(2))
        endif
        zmin=s2min/engy**2
        if(ihww.gt.0)then
        zmax=1
        xpmax=1
        ig1=7
        ig2=7
        r1=0
        do i1=1,ig1
        do m1=1,2
          z=zmin*(zmax/zmin)**(.5+tgss(ig1,i1)*(m1-1.5))
          xpmin=max(z,xpmin0)
          r2=0
          if(xpmin.lt.xpmax)then
          do i2=1,ig2
          do m2=1,2
            xp=xpmin*(xpmax/xpmin)**(.5+tgss(ig2,i2)*(m2-1.5))
            xm=z/xp
            r2=r2+wgss(ig2,i2)*ffsigiut(xp,xm,max(1,jaaa),je1,je2)
          enddo
          enddo
          endif
          r2=r2*0.5*log(xpmax/xpmin)
          r1=r1+wgss(ig1,i1)*r2*z
        enddo
        enddo
        r1=r1*0.5*log(zmax/zmin)
        res=  r1 * factk /sigine * .0389 *10  !((hbar*c)**2=0.389 GeV^2.mbarn, and ffsig in GeV^-2)
        sum=res
        endif
   !.......plot
        xx2max = 1
        xx1min = sqrt(max(zmin,xpmin0))
        xx1max = 1.
        nbins  = nbixp

        write(ifhi,'(a,3i1)') 'openhisto name ffsig',jaa,jex,ncy
        write(ifhi,'(a,2e11.3)')'xrange ',xx1min,xx1max
        write(ifhi,'(a)') 'yrange auto auto htyp lin '
        write(ifhi,'(a)')'htyp lin xmod '//modu//' ymod log'
        write(ifhi,'(a)') 'txt "xaxis x+?IB!         "              '
        write(ifhi,'(a)') 'txt "yaxis dn?semi! / dx+?IB!  "'
        write(ifhi,'(a,2i1,a)')'txt "title ffsig + MC   (',jaa,jex,')"'
        write(ifhi,'(a,e22.14)')'histoweight ',dble(hww)
        write(ifhi,'(a)') 'array 2'
        del=(xx1max-xx1min)/nbins
        do ii=1,nbins
          if(kkk.eq.1)then
            xx1=xx1min*(xx1max/xx1min)**((ii-0.5)/float(nbins))
          else
            xx1=xx1min+(ii-0.5)*del
          endif
          sig=0
          if(ihww.gt.0)then
            if(jaaa.eq.11)then
              xx2min = (4.*q2sft)*xpmin0/xx1
            else
              xx2min = (4.*max(q2cmin(1),q2cmin(2)))*xpmin0/xx1
            endif
            xx2min=max(xx2min,xlub)
            ig2=7
            r2=0
            do i2=1,ig2
            do m2=1,2
            xx2=xx2min*(xx2max/xx2min)**(.5+tgss(ig2,i2)*(m2-1.5))
            r2=r2+wgss(ig2,i2)*ffsigiut(xx1,xx2,max(1,jaaa),je1,je2)*xx2
            enddo
            enddo
           sig=r2*0.5*log(xx2max/xx2min)
           sig   = 2. * sig * factk /sigine *10* .0389 !((hbar*c)**2=0.389 GeV^2.mbarn, and ffsig in GeV^-2) * 2. (for number of jets)
          endif
          write(ifhi,'(2e12.4)')xx1,sig
        enddo
        write(ifhi,'(a)')  '  endarray'


       elseif(iii.eq.4)then

        write(ifhi,'(a)')       '!----------------------------------'
        write(ifhi,'(a,3i1)')   '!   OB    ',jaa,jex
        write(ifhi,'(a)')       '!----------------------------------'

        y2     = rapo
        nbins  = nbipt
        nup=nbins
        if(nint(xpar9).gt.0.)nup=nint(xpar9)
        if(nint(xpar13).gt.0.)nlow=nint(xpar13)
        if(nlow.gt.nup)stop'ERROR xpar13 > xpar9'
        sx=engy**2
        xup=ptob(nup)
        ncntEmsP21=ncntEmsP21+1
        !jjmx=5
        !if(ihq.ne.0)jjmx=1
        !++++++++++++++++++++++++++++++++++++++++++
        !bbb=pifpartone(1,dble(0.1),5.,0,2,0)
        !print*,'TEST-CTEQ-EPOS'
        !do i=4,1,-1
        !  x=1./10**i
        !  aaa=xCTEQ2(1,dble(x),5.,0,2,0)
        !  bbb=pifpartone(1,dble(x),5.,0,2,0)
        !  aaa1=xCTEQ2(1,dble(x),5.,1,2,0)
        !  bbb1=pifpartone(1,dble(x),5.,1,2,0)
        !  print*,'TEST-CTEQ-EPOS',aaa,bbb,aaa1,bbb1
        !enddo
        !++++++++++++++++++++++++++++++++++++++++++
        do jjx=jjlow,jjup
        jj=jjselect(jjx) 
        if(jjx.le.9)then
         write(ifhi,'(a,i1,$)')'openhisto name jet',jjx
        else
         write(ifhi,'(a,i2,$)')'openhisto name jet',jjx
        endif 
        if(ncntEmsP21.le.9)then
          write(ifhi,'(i1)')ncntEmsP21
        else
          write(ifhi,'(i2)')ncntEmsP21
        endif 
        write(ifhi,'(a,e12.6)')'xrange 0 ',xup
        if(xpar8.gt.0.001)write(ifhi,'(a,e12.6)')'xmin ',xpar8
        write(ifhi,'(a)')'yrange auto auto xmod lin ymod log '
        write(ifhi,'(a)') 'txt "xaxis pt?OB!         "           '
        write(ifhi,'(a)') 'txt "yaxis dn?ptn! / dpt?OB!  "'
        if(jj.eq.2)write(ifhi,'(a)')'htyp lkv'
        if(jj.eq.3)write(ifhi,'(a)')'htyp lkb'
        if(jj.eq.4)write(ifhi,'(a)')'htyp lgv'
        if(jj.eq.5)write(ifhi,'(a)')'htyp lgb'
        if(jj.ge.30)write(ifhi,'(a)')'htyp liv'
        if(ipar2.eq.3.or.ipar2.eq.4.or.ipar2.eq.6
     .  .or.ipar2.eq.31.or.ipar2.eq.41)then
          if(jj.eq.1)write(ifhi,'(a)')'htyp lin'
          if(jj.eq.6)write(ifhi,'(a)')'htyp lin'
          if(jj.eq.7)write(ifhi,'(a)')'htyp lin'
          if(jj.eq.8)write(ifhi,'(a)')'htyp lin'
          if(jj.eq.9)write(ifhi,'(a)')'htyp lin'
          if(jj.eq.10)write(ifhi,'(a)')'htyp lin'
        else
          if(jj.eq.1)write(ifhi,'(a)')'htyp lrv'
          if(jj.eq.6)write(ifhi,'(a)')'htyp lwb'
          if(jj.eq.7)write(ifhi,'(a)')'htyp lrv'
          if(jj.eq.8)write(ifhi,'(a)')'htyp lwb'
          if(jj.eq.9)write(ifhi,'(a)')'htyp lrv'
          if(jj.eq.10)write(ifhi,'(a)')'htyp lwb'
        endif

        write(ifhi,'(a,e22.14)')'histoweight ',dble(hww)
        if(xpar10.ne.0.)write(ifhi,'(a,2x,e12.6)')'set factor',xpar10
        write(ifhi,'(a)')'array 2'
        sum=0.
        if(ncy.eq.0.and.ihq.eq.0.and.ioTestFact.eq.0)then
          qmin=max(q2amin(1),q2amin(2))
          dqq=(0.25*sx/qmin)**(0.5/float(nup))
          do i=0,nup
            qq=qmin*dqq**(2.*i+1.)
            sum=sum+ffsigi(qq,50.,0,1)*(dqq**(i+1)-dqq**i)
          enddo
        endif
        do i=nlow,nup
         pt=ptu*(pto/ptu)**((i-0.5)/nbipt) !log binning needed at low pt !!!
         !pt=ptu+(pto-ptu)*(i-0.5)/nbipt
         if(pt.le.pto)then
          sig=0
c          if(ihww.gt.0)then
          sig=1.
          if(jj.eq.1.and.ihwfac.eq.1)then
            sig=ffsigi(pt**2,y2,ihq,1) ! our stuff
          elseif(jj.eq.2)then
            if(engy.ge.10.)sig=psjpdf(psjvrx,pt**2,sx,y2) !grv
          elseif(jj.eq.3)then
            if(engy.ge.10.)sig=psjpdf(psjwox,pt**2,sx,y2) !duke-owens
          elseif(jj.eq.4)then
            if(engy.ge.10.)sig=psjpdf(psjgrv08x,pt**2,sx,y2)  !GRV08
          elseif(jj.eq.5)then
            if(engy.ge.10.)sig=psjpdf(psjcteq6x,pt**2,sx,y2)  !CTEQ6
          elseif(jj.eq.6)then
            if(engy.ge.10.)then
              !if(i.eq.1)write(ifmt,'(a)')'CTEQ start'
              !if(i.eq.1)call clop(3)
              sig=ffsigi(pt**2,y2,ihq,6)   !CTEQ14
              !sigour=ffsigi(pt**2,y2,ihq,1)
              !write(ifmt,'(a,3e11.2)')'CTEQ',pt,sig,sigour
              !call clop(3)
            endif 
          elseif(jj.eq.7)then
            if(engy.ge.10.)sig=ffsigi(pt**2,y2,4,1) !our stuff charm
          elseif(jj.eq.8)then
            if(engy.ge.10.)sig=ffsigi(pt**2,y2,4,6) !CTEQ14 charm
          elseif(jj.eq.9)then
            if(engy.ge.10.)sig=ffsigi(pt**2,y2,5,1) !our stuff bottom
          elseif(jj.eq.10)then
            if(engy.ge.10.)sig=ffsigi(pt**2,y2,5,6) !CTEQ14 bottom
          endif
          !factors: ((hbar*c)**2=0.389 GeV^2.mbarn, and ffsig in GeV^-2)
          !  division by delta_rap: 1/2 already done in ffsigi and psjpdf (integration from 0 to y2)
          !*2. CS for 1 hard process (dijet) but we want CS for a single jet to compare to data (so 2 jets for 1 dijet)
          sig=sig*factk/sigine*10*.0389/y2*2.
c          sig=sig*factk/xsectionpar()*10*.0389/y2*2.  !use parametrized cross section to avoid sigine dependence
c          endif
          write(ifhi,'(2e12.4)')pt,sig
         endif
        enddo
        write(ifhi,'(a)')       '  endarray'
        if(ncy.eq.0.and.ihq.eq.0)then
          sum=sum*factk/sigine*10*.0389*4. !number of hard pom from theo (int over y (*2) )
         write(ifhi,'(a,f6.3,a)')'text 0.15 0.15 "N?hard! theo ',sum,'"'
        endif
        if(jjx.ne.jjup)write(ifhi,'(a)')       'closehisto'
        if(jjx.ne.jjup)write(ifhi,'(a)')  'plot 0-'
        enddo

        if(ifacpdf4.eq.1)then
          call pdfparamset(1, facpdf4_save )
          call pdfparamset(3, facpdf5_save )
        endif

       elseif(iii.eq.6.or.iii.eq.7)then

        write(ifhi,'(a)')       '!----------------------------------'
        write(ifhi,'(a,3i1)')   '!   GB    ',jaa,jex
        write(ifhi,'(a)')       '!----------------------------------'
        write(ifhi,'(a,3i1)') 'openhisto name gb',jaa,jex,ncy
        write(ifhi,'(a,2f12.3)') 'xrange ',gbug,gbog
        write(ifhi,'(a)') 'yrange auto auto htyp lka xmod lin ymod log'
        write(ifhi,'(a)') 'txt "xaxis lg(gb)         "              '
        write(ifhi,'(a)') 'txt "yaxis dn/d_lg(gb) "'
        if(iii.eq.6)
     .  write(ifhi,'(a,2i1,a)')'txt "title  lg gb min  (',jaa,jex,')"'
        if(iii.eq.7)
     .  write(ifhi,'(a,2i1,a)')'txt "title  lg gb max  (',jaa,jex,')"'
        write(ifhi,'(a,e22.14)')'histoweight ',dble(hww)
        write(ifhi,'(a)') 'array 2'
        do ii=1,2
          xx1=(2*ii-3)*0.0001
          sig=0
          if(ihww.gt.0)then
          sig   = 0.001+(ii-1)*10
          endif
          write(ifhi,'(2e12.4)')xx1,sig
        enddo
        write(ifhi,'(a)')  '  endarray'

       elseif(iii.eq.8)then

        write(ifhi,'(a)')       '!----------------------------------'
        write(ifhi,'(a,3i1)')   '!   xb vs xe    ',jaa,jex
        write(ifhi,'(a)')       '!----------------------------------'
        write(ifhi,'(a,3i1)') 'openhisto name xbvsxe',jaa,jex,ncy
        write(ifhi,'(a)')'htyp lba xmod '//modu//' ymod log'
        write(ifhi,'(a,2e11.3)')'xrange ',xlub,xlob
        write(ifhi,'(a)')'txt "xaxis x+?PE!"'
        write(ifhi,'(a)')'txt "yaxis x+?IB! "'
        write(ifhi,'(a,2i1,a)')'txt "title  correl PE-IB (',jaa,jex,')"'
        write(ifhi,'(a,e22.14)')'histoweight ',dble(hww)
        write(ifhi,'(a)') 'array 2'
        do i=j1,j2
          u=xlp(kkk,i)
          z=0.1*u
c          if(u.gt.0.)z=z/u**0.2
          write(ifhi,'(2e11.3)')u,z
        enddo
        write(ifhi,'(a)')  '  endarray'

       elseif(iii.eq.9)then

        write(ifhi,'(a)')       '!----------------------------------'
        write(ifhi,'(a,3i1)')   '!   xb vs xe    ',jaa,jex
        write(ifhi,'(a)')       '!----------------------------------'
        write(ifhi,'(a,3i1)') 'openhisto name xbxe',jaa,jex,ncy
        write(ifhi,'(a)')'htyp lba xmod '//modu//' ymod lin'
        write(ifhi,'(a,2e11.3)')'xrange ',xlub,xlob
        write(ifhi,'(a,2e11.3)')'yrange ',1e-6,0.25
        write(ifhi,'(a)')'txt "xaxis x+?IB!"'
        write(ifhi,'(a)')'txt "yaxis P(x+?IB!) "'
        write(ifhi,'(a,2i1,a,e11.3,a)')'txt "title  x+?PE! ('
     .            ,jaa,jex,')=',xlp(kkk,j1),'"'
        write(ifhi,'(a,e22.14)')'histoweight ',dble(hww)
        write(ifhi,'(a)') 'array 2'
        do i=max(1,j1-1),min(j2+1,nbixp)
          u=xlp(kkk,i)
          z=0.
          if(i.eq.j1)z=0.2
          write(ifhi,'(2e11.3)')u,z
        enddo
        write(ifhi,'(a)')  '  endarray'

       endif

       x=0.1+(min(3,iii)-2)*0.30
       y=0.2+(min(3,iii)-2)*0.50

       if(ihq.eq.0)then

       write(ifhi,'(a)')  'closehisto'

       !########################## MC ######################### 

       if(iSKIP_MC.eq.1)goto 789 !SKIP_MC

       write(ifhi,'(a)')  'plot 0-'
       if(iii.eq.2)then
        ncntEmsP22=ncntEmsP22+1
        ncntxx=ncntEmsP22
       elseif(iii.eq.3)then
        ncntEmsP23=ncntEmsP23+1
        ncntxx=ncntEmsP23
       elseif(iii.eq.4)then
        ncntEmsP24=ncntEmsP24+1
        ncntxx=ncntEmsP24
       elseif(iii.eq.5)then
        ncntEmsP25=ncntEmsP25+1
        ncntxx=ncntEmsP25
       elseif(iii.eq.6)then
        ncntEmsP26=ncntEmsP26+1
        ncntxx=ncntEmsP26
       elseif(iii.eq.7)then
        ncntEmsP27=ncntEmsP27+1
        ncntxx=ncntEmsP27
       elseif(iii.eq.8)then
        ncntEmsP28=ncntEmsP28+1
        ncntxx=ncntEmsP28
       elseif(iii.eq.9)then
        ncntEmsP29=ncntEmsP29+1
        ncntxx=ncntEmsP29
       endif

       if(iii.ne.8.and.iii.ne.9)then !multiply with ff via calc
       write(ifhi,'(a)') "!-----------------------------"
       write(ifhi,'(a)') "! MC Norm   "  !cKW21 Bugfix: ff must go to histogram (to use "add" properly) 
       write(ifhi,'(a)') "!-----------------------------"
       write(ifhi,'(a,4i1)')'openhisto name P2Norm1',iii,jaa,jex,ncntxx
       write(ifhi,'(a,e22.14)')'histoweight ',dble(hw)
       write(ifhi,'(a)')  'array 2'
       write(ifhi,'(4e11.3)')1.,ff(1)
       write(ifhi,'(a)')  'endarray'
       write(ifhi,'(a)')  'closehisto'
       write(ifhi,'(a,4i1)')'openhisto name P2Norm2',iii,jaa,jex,ncntxx
       write(ifhi,'(a,e22.14)')'histoweight ',dble(hw)
       write(ifhi,'(a)')  'array 2'
       write(ifhi,'(4e11.3)')1.,ff(2)
       write(ifhi,'(a)')  'endarray'
       write(ifhi,'(a)')  'closehisto'
       write(ifhi,'(a,4i1)')'openhisto name P2Norm3',iii,jaa,jex,ncntxx
       write(ifhi,'(a,e22.14)')'histoweight ',dble(hw)
       write(ifhi,'(a)')  'array 2'
       write(ifhi,'(4e11.3)')1.,ff(3)
       write(ifhi,'(a)')  'endarray'
       write(ifhi,'(a)')  'closehisto'
       write(ifhi,'(a,4i1)')'openhisto name P2Norm4',iii,jaa,jex,ncntxx
       write(ifhi,'(a,e22.14)')'histoweight ',dble(hw)
       write(ifhi,'(a)')  'array 2'
       write(ifhi,'(4e11.3)')1.,ff4
       write(ifhi,'(a)')  'endarray'
       write(ifhi,'(a)')  'closehisto'
       write(ifhi,'(a,4i1)')'openhisto name P2Norm5',iii,jaa,jex,ncntxx
       write(ifhi,'(a,e22.14)')'histoweight ',dble(hw)
       write(ifhi,'(a)')  'array 2'
       write(ifhi,'(4e11.3)')1.,ff5
       write(ifhi,'(a)')  'endarray'
       write(ifhi,'(a)')  'closehisto'
       endif

       write(ifhi,'(a)') "!-----------------------------"
       write(ifhi,'(a)') "! MC   "
       write(ifhi,'(a)') "!-----------------------------"

       if(iii.eq.2)then
        write(ifhi,'(a,3i1)')'openhisto name dndxPE'//cii,jaa,jex
     .   ,ncntEmsP22
       elseif(iii.eq.3)then
        write(ifhi,'(a,3i1)')'openhisto name dndxIB',jaa,jex
     .    ,ncntEmsP23
       elseif(iii.eq.4)then
        write(ifhi,'(a,3i1)')'openhisto name dndptOB',jaa,jex
     .    ,ncntEmsP24
       elseif(iii.eq.5)then
        write(ifhi,'(a,3i1)')'openhisto name dndxPE'//cii,jaa,jex
     .   ,ncntEmsP25
       elseif(iii.eq.6)then
        write(ifhi,'(a,3i1)')'openhisto name dndlgGBmin',jaa,jex
     .   ,ncntEmsP26
       elseif(iii.eq.7)then
        write(ifhi,'(a,3i1)')'openhisto name dndlgGBmax',jaa,jex
     .   ,ncntEmsP27
       elseif(iii.eq.8)then
        write(ifhi,'(a,3i1)')'openhisto name xbvsxPEmc',jaa,jex
     .   ,ncntEmsP28
       elseif(iii.eq.9)then
        write(ifhi,'(a,3i1)')'openhisto name xbxPEmc',jaa,jex
     .   ,ncntEmsP29
       endif
       write(ifhi,'(a,e22.14)')'histoweight ',dble(hw)
       if(xpar8.gt.0.001)write(ifhi,'(a,e12.6)')'xmin ',xpar8
       if(iii.eq.4)then
         write(ifhi,'(a)')     'htyp prc '
       elseif(iii.eq.6.or.iii.eq.7.or.iii.eq.9)then
         write(ifhi,'(a)')     'htyp lru '
       else
         write(ifhi,'(a)')     'htyp prs '
       endif
       write(ifhi,'(a)')     'array 2'
       sum=0
       imax=nbixp
       if(iii.eq.4)imax=nup
       if(iii.eq.6)imax=nbigb
       !!!!!1if(iii.eq.5)imax=nbish ????? now 5 used for other stuff
       jaa11=jaa1
       jaa12=jaa2
       if(jaa11.eq.11)then
         jaa11=0
         jaa12=0
       endif
       !integral for normalization to get x+IB distribution for a given x+PE (j1)
       z1=0.
       if(iii.eq.9)then
         do i=1,imax
           do jaai=jaa11,jaa12
             do jexi=jex1,jex2
               z1=z1+wxs(ncy,kkk,jaai,jexi,j1,i)
             enddo
           enddo
         enddo
       endif
       iff45=0
       iff6=0
       do i=1,imax
         if(iii.eq.4)then
           u=ptob(i)
         elseif(iii.eq.6.or.iii.eq.7)then
           u=gbval(i)
         else
           u=xlp(kkk,i)
         endif
         z=0.
         if(iii.eq.4)then
           !del=(pto-ptu)!/nbipt
           del=ptu*(pto/ptu)**((i-0.0)/nbipt) !log binning needed at low pt !!!
     .        -ptu*(pto/ptu)**((i-1.0)/nbipt)
           dely=rapo-rapu    !to normalize by rap range
         elseif(iii.eq.6.or.iii.eq.7)then
           del=gbval(2)-gbval(1)
           dely=1.
         else
           del=dxlp(kkk,i)
           dely=1.
         endif
         if(iii.eq.2.or.iii.eq.3.or.iii.eq.8.or.iii.eq.9)then
           if(iii.eq.8)z1=0.
           do j=j1,j2
             do jaai=jaa11,jaa12
               do jexi=jex1,jex2
                 !integral over xm
                 if(iii.eq.2)z=z+wxe(ncy,kkk,jaai,jexi,i,j)
                 if(iii.eq.3)z=z+wxb(ncy,kkk,jaai,jexi,i,j)
                 if(iii.eq.8)then
                 !mean value of x+IB for a given x+PE
                   z1=z1+wxs(ncy,kkk,jaai,jexi,i,j) !normalization
                   z=z+wxs(ncy,kkk,jaai,jexi,i,j)*xlm(kkk,j) !x+IB value
                 endif
                 if(iii.eq.9)z=z+wxs(ncy,kkk,jaai,jexi,j1,i)
               enddo
             enddo
           enddo
           if(iii.le.3.and.ncy.gt.0)iff45=1  !multiply by <Ncol>_minbias to be compared to inclusive cross-section
         elseif(iii.eq.5)then
           do k=1,nbixp
           do j=1,nbixm
             xp=xlp(kkk,k)
             xm=xlm(kkk,j)
             x=xp*xm
             if(x.gt.uxlp(kkk,i).and.x.lt.oxlp(kkk,i))then
              do jaai=jaa11,jaa12
               do jexi=jex1,jex2
                 z=z+wxe(ncy,kkk,jaai,jexi,k,j)
               enddo
              enddo
             endif
           enddo
           enddo
           if(ncy.gt.0)iff45=1 !multiply by <Ncol>_minbias to be compared to inclusive cross-section
         elseif(iii.eq.4)then
           if(jaa.eq.5)then
             jaa11=0
             jaa12=4
           endif
           do jaai=jaa11,jaa12
             do jexi=jex1,jex2
               z=z+wptob(ncy,jaai,jexi,i)
             enddo
           enddo
           if(ncy.gt.0)iff45=1 !multiply by <Ncol>_minbias to be compared to inclusive cross-section
         elseif(iii.eq.6.or.iii.eq.7)then
           do jaai=jaa11,jaa12
             do jexi=jex1,jex2
               if(iii.eq.6)z=z+wgb(ncy,1,jaai,jexi,i)
               if(iii.eq.7)z=z+wgb(ncy,2,jaai,jexi,i)
             enddo
           enddo
         endif
         if(iii.eq.8.or.iii.eq.9)then
           if(z1.gt.0.)then
             z=z/z1
           else
             z=0.
           endif
         elseif(iii.ne.8.and.iii.ne.9)then !use calc to multiply with ff
           z=z/del/max(nstat,1)/dely
         else
           stop'ERROR 09072021'
         endif
         if(iii.ne.4.or.xpar6.le.0.001.or.u.le.xpar6)
     .   write(ifhi,'(2e11.3)')u,z
         sum=sum+z*del
       enddo
       write(ifhi,'(a)')    'endarray'
       !---------------
       if(iii.ne.8.and.iii.ne.9)then !multiply with ff
         write(ifhi,'(a)')    'closehisto'
         write(ifhi,'(a)')    'openhisto'
         if(iii.eq.4)then
           write(ifhi,'(a)')     'htyp prc '
         elseif(iii.eq.6.or.iii.eq.7)then
           write(ifhi,'(a)')     'htyp lru '
         else
           write(ifhi,'(a)')     'htyp prs '
         endif
         write(ifhi,'(3(a,4i1),$)')
     .   'calc *1 * P2Norm1',iii,jaa,jex,ncntxx
     .   ,      ' * P2Norm2',iii,jaa,jex,ncntxx
     .   ,      ' / P2Norm3',iii,jaa,jex,ncntxx
         if(iff45.eq.1.and.jff45.eq.1)then
           write(ifhi,'(2(a,4i1),$)')
     .          ' * P2Norm4',iii,jaa,jex,ncntxx
     .     ,    ' / P2Norm5',iii,jaa,jex,ncntxx
         endif
         write(ifhi,'(a)')' ; '
       endif
       !---------------
       iiix=iii
       if(iii.eq.5)iiix=2
       yrf=(min(3,iiix)-2)*0.50
       x=0.1+(min(3,iiix)-2)*0.25
       y=0.4+yrf
c       if(engy.gt.100)then
c       write(ifhi,'(a,2f5.2,a,f6.3,a)')'text',x,y,' "simu ',sum,'"'
c       else
c       write(ifhi,'(a,2f5.2,a,f6.5,a)')'text',x,y,' "simu ',sum,'"'
c       endif
       y=0.15+yrf
       c3='7.1'
       if(zavcoll.lt.1000)c3='6.1'
       if(zavcoll.lt.100)c3='6.2'
       if(zavcoll.lt.10)c3='6.3'
       if(iii.eq.4)
     . write(ifhi,'(a,2f5.2,a,f'//c3//',a)')'text'
     . ,x,y+0.15,' "Ncoll av',zavcoll,'"'
       y=0.05+yrf
       c3='7.1'
       if(max(q2cmin(1),q2cmin(2)).lt.1000)c3='6.1'
       if(max(q2cmin(1),q2cmin(2)).lt.100)c3='6.2'
       if(max(q2cmin(1),q2cmin(2)).lt.10)c3='6.3'
       yi4=0
       if(iii.eq.4)then
        yi4=0.15
        y=y+yi4
        write(ifhi,'(a,2f5.2,a,f'//c3//',a)')'text'
     .  ,x,y,' "Q2sat av',q2cmin(1),'"'
        write(ifhi,'(a,2f5.2,a,f'//c3//',a)')'text'
     .  ,x+0.33,y-0.10,' "',q2cmin(2),'"'
        write(ifhi,'(a,2f5.2,a,f'//c3//',a)')'text'
     .  ,x+0.33,y-0.20,' "',q2amin(1),'"'
        write(ifhi,'(a,2f5.2,a,f'//c3//',a)')'text'
     .  ,x+0.33,y-0.30,' "',q2amin(2),'"'
       endif
       y=0.25+yrf+yi4
       write(ifhi,'(a,2f5.2,a,i2,a)')'text',x,y,' "M',ncy,'"'
       !if(iii.eq.4)write(ifhi,'(a)')'text 0.4 0.5 "P2 simu"'
       endif                     !ihq
       write(ifhi,'(a)')    'closehisto'

 789   continue !SKIP_MC

       if(iii.eq.2)then
        if(iG100.eq.1)then!~~~~~~~~~~
        write(ifhi,'(a)')      'plot 0-'
        write(ifhi,'(a,3i1,a)')
     .  'openhisto name dndxPE',jaa,jex,ncy,'G100 set factor 1e5'
        write(ifhi,'(a)')     'htyp lkv'
        write(ifhi,'(a,e22.14)')'histoweight ',dble(hw)
        write(ifhi,'(a)')     'array 2'
        imax=nbixp
        do i=1,imax
         u=xlp(kkk,i)
         z=0.
         do j=j1,j2
           do jaai=jaa1,jaa2
             do jexi=jex1,jex2
               z=z+wxi(ncy,kkk,jaai,jexi,i,j)
              enddo
           enddo
         enddo
         del=dxlp(kkk,i)
         dely=1.
         z=z/del*ff(1)*ff(2)/ff(3)/max(nstat,1)/dely   !ff should be put into histo to do "add" properly
         write(ifhi,'(2e11.3)')u,z
        enddo
        write(ifhi,'(a)')    '  endarray'
        write(ifhi,'(a)')    'closehisto'
        endif
       endif

       if(ncy.lt.ncy2.or.nint(xpar8).eq.1)write(ifhi,'(a)')'plot 0'

      enddo                     ! ncy

      else

          write(ifch,*)'\n\n  ERROR 27092012f \n\n'
                   stop'\n\n  ERROR 27092012f \n\n'
      endif

      return
      end

c----------------------------------------------------------------------
      subroutine xEmsTest3
c----------------------------------------------------------------------
#include "aaa.h"
#include "sem.h"
      double precision fzeroGlu
      write(ifmt,'(14a)')('#####',i=1,14)
      write(ifmt,'(a,$)')'xEmsTest3'
      pt=10
      y2=1
      !fv=ffsigi(pt**2,y2,0,1)
      !print*,'ffsigi=',fv

      fz1=fzeroGlu(.1d0,q2nmin, 1)
      fz2=fzeroGlu(.1d-2,q2nmin, 1)
      fz3=fzeroGlu(.1d-4,q2nmin, 1)
      fz4=fzeroGlu(.1d-6,q2nmin, 1)
      print*,'fzeroGlu=',fz1,fz2,fz3,fz4

      write(ifmt,'(14a)')('#####',i=1,14)
      stop'normal STOP in xEmsTest3'
      end

c----------------------------------------------------------------------
      subroutine xEmsTest2
c----------------------------------------------------------------------
#include "aaa.h"
#include "sem.h"
      data ncount/0/
      save ncount
      ncount=ncount+1
      if(ncount.gt.1)return

      write(ifmt,'(14a)')('#####',i=1,14)
      write(ifmt,'(a,$)')'xEmsTest2  '
c         return

      !call ipoAllTables(6.25)
      !q2pmin(1)=q2nmin
      !q2pmin(2)=q2pmin(1)
      !call ipoBasics(q2pmin)
      !call ipoOm5Tables(1)
      sypp=engy*engy
      !--------set q2----------------
      q2pmin(1)=q2nmin
      q2pmin(2)=q2pmin(1)
      q2cmin(1)=q2pmin(1)
      q2cmin(2)=q2pmin(2)
      q1=q2cmin(1)
      q2=q2cmin(2)
      qmx=max(q1,q2)
      print*,'    q2',q2cmin
      call ipoBasics(q2pmin)
      call ipoCSzeroTables(q2cmin)
      call ipoOm5Tables(1)
      dmy=pijet( 2 ,sypp,0,0) !to make csjet table 
      print*,' '
      !------start tests--------------
      do i=1,3
      xm=0.8
      xp=0.8
      xh=xm*xp
      sy=sypp*xh
      b=0.5*(i-1)
      rp=r2had(iclpro)+r2had(icltar)+slopom*log(max(1.,sy))
      zb=exp(-b**2/(4.*.0389*rp)) 
      print*,'b,x = ',b,xh
      !------------comparing om52pp and om11pp--------------- 
      x1=om52pp(sy,1.,1.,b,0,2,2) / ( factk / sigine * 10 )
      x2=om11pp(sy,1d0,zb,0)* sy**delh
      print*,'G(om52pp),G(om11pp) = ', x1,x2
      !------------testing interpolation--------------- 
      m2=0
      l2=0   
      call eolpib( m2 , l2 , ml ) !KW1808b in pijet
      call eolpi( ml , m1 , l1 ) !KW1808c     
      call psjti0(sy,CSipoX,sdy,m2,l2) 
      CSipo = pijet( 2 ,sy,m2,l2) + pijet( 1 ,sy,m2,l2)  !used in 
     .     +pijet( -1,sy,m2,l2) + pijet( 0 ,sy,m2,l2)    !om52pp
      CSexact = psjet(q1,q2,qmx,sy,m1,l1,0)
     .        + psjet1(q1,q2,qmx,sy,m1,l1,0)
     .        + psjet1(q2,q1,qmx,sy,m1,l1,0) 
     .        + psborn(q1,q2,qmx,sy,m1,l1,0,0) 
      print*,'CSipoX,CSipo,CSexact = ',CSipoX,CSipo,CSexact
      enddo 
      write(ifmt,'(14a)')('#####',i=1,14)
      print*,'return from xEmsTest2'
      return

      !***************************************
      !
      ! for these tests:  use "protect inirj"
      !
      !***************************************
      !x1=om52pp(sy,xp, 1., b ,  1 ,  2,2)
      !x2=om11pp(sy,xp,     zb,  1       )
      !x1=om52pp(sy,1., xm, b ,  2 ,  2,2)
      !x2=om11pp(sy,xm,     zb,  2      )
      !x1=om52pp(sy,xp, xm, b ,  3 ,  2,2)
      !x2=om11pp3(sy,xp,xm              )
      !***************************************
      !
      ! for these tests:  do NOT use "protect inirj"
      !
      !     because ffom11 uses the table inirj
      !
      !***************************************
      x1=ffom12(xp,xm,  b, 1,1    ,2,2)
      call ipoOm5Tables(1)
      x2=ffom11(xp,xm,  b, 1,1   )
      x3=ffom11(xp,xm,  b, -1,-1   )
      print*,' '
      print*,x1,x2,x3
      x1=ffom12(xp,xm,  b, 2,2    ,2,2)
      call ipoOm5Tables(1)
      x2=ffom11(xp,xm,  b, 2,2   )
      print*,x1,x2
      x1=ffom12(xp,xm,  b, 3,3    ,2,2)
      call ipoOm5Tables(1)
      x2=ffom11(xp,xm,  b, 3,3   )
      print*,x1,x2
      x1=ffom12(xp,xm,  b, 4,4    ,2,2)
      call ipoOm5Tables(1)
      x2=ffom11(xp,xm,  b, 4,4   )
      print*,x1,x2
      write(ifmt,'(14a)')('#####',i=1,14)
      call ipoCSzeroTables(q2cmin)
      x1=om52pp(sy,1.,1.,b,0,2,2) / ( factk / sigine * 10 )
      x2=om11pp(sy,1d0,zb,0)* sy**delh
      print*,x1,x2
c      x1=om52pp(sy,1.,1.,b,10,2,2) / ( factk / sigine * 10 )
c      x2=om11pp4(sy,zb,10)* sy**delh
c      print*,x1,x2
      write(ifmt,'(14a)')('#####',i=1,14)
      stop'normal STOP in xEmsTest2'
      end

c----------------------------------------------------------------------
      subroutine xEmsTest
c----------------------------------------------------------------------
c  Not updated with Q2s
c----------------------------------------------------------------------
#include "aaa.h"
#include "sem.h"
      double precision Fgfactor,xp,xm
      sypp=engy*engy
      !call ipoAllTables(6.25)
      xm=0.1d0
      xp=0.1d0
      xh=xm*xp
      sy=sypp*xh
      yp=0.5*log(xp/xm)
      b=0.5
      rp=r2had(iclpro)+r2had(icltar)+slopom*log(max(1.,sy))
      rh=r2had(iclpro)+r2had(icltar)-slopom*log(0.3)
      zb=exp(-b**2/(4.*.0389*rp))
      do i=idxDmin(iomega),idxDmax(iomega)
       zzp=xpar98
       zzt=xpar99
       call Gfunpar(zzp,zzt,0.,0.,1,i,beff,sypp,alpx,betx,betpx
     . ,epsp,epst,epss,gamv)
      enddo
      xxx1=Fgfactor(xp,xm,1)
      epsG1=epssUni(1)
      epsGp1=epspUni(1)
      epsGt1=epstUni(1)
      do i=idxDmin(iomega),idxDmax(iomega)
       zzp=xpar98
       zzt=xpar99
       call Gfunpar(zzp,zzt,0.,0.,1,i,0.5,sypp,alpx,betx,betpx
     . ,epsp,epst,epss,gamv)
      enddo
      xxx3=Fgfactor(xp,xm,1)
      epsG3=epssUni(1)
      epsGp3=epspUni(1)
      epsGt3=epstUni(1)
      write(ifmt,'(14a)')('#####',i=1,14)
      write(*,'(a,$)')'xEmsTest 111111111111111111111111111111111111'
      x1=ffom12a( sngl(xp), sngl(xm),         1,1,     2,2)
      print*,' '
      call ipoOm5Tables(1)
      x2=ffom11a( sngl(xp), sngl(xm), 1,1)
      x3=b*ffom11(sngl(xp),sngl(xm),b, 1,1)
      print*,'111 ffom12a =  ',x1
      print*,'222 ffom11a =  ',x2
      print*,'333 b*ffom11 =', x3
      print*,'- - - - - - - - - - - - - - - - - - - - - -'
      x1a=om52pp(sy,1.,1.,0.,0,2,2)
      x3a=psvin(sy,xp,xm,zb,1)
   !G_basic:
      x1a1=om52pp(sy,1.,1.,0.,0,2,2)
     . / ( factk * .0389 / sigine * 10 / 2. )
      x3a1=om11pp(sy,1.d0,zb,0)            !!
     . / ( .0389 / sy**delh )
     . / ( 1./(4.*pi)  *zb**(rp/rh)/rh/.0389 )
   !other factors
      x1a2=factk * .0389
      x1a3=1./ sigine * 10 /2.
      x1b=2*engy**epsG1*engy**epsG1
      x1c=gamhad(iclpro)*gamhad(icltar)*xp**(-alppar)*xm**(-alppar)
      x1d=(1-xp)**alplea(iclpro)*(1-xm)**alplea(icltar)
      x1e=xp**(epsGp1)*xm**(epsGt1)

      x3a2=1./(4.*pi)  *zb**(rp/rh)/rh/.0389  * .0389 /sy**delh
      x3a3=factk*sy**delh                              !G
      x3a4=gamhad(iclpro)*gamhad(icltar)*xp**(-alppar)*xm**(-alppar)
      x3b=b*.5*2.
      x3c=(1-xm)**alplea(icltar)*(1-xp)**alplea(iclpro)
      x3d=xp**epsGp3*xm**epsGt3
      print*,'111 om52pp =',x1a,'  ',x1a*x1b*x1c*x1d*x1e*xxx1
      print*,'333 psvin =',x3a,'  ',x3a*x3b*x3c*x3d*xxx3
      print*,'- - - - - - - - - - - - - - - - - - - - - -'
      print*,'111 G_bas(om52pp)',x1a1,  x1a1*x1a2,   x1a1*x1a2*x1a3
      print*,'333 G_bas(om11pp)',x3a1,x3a1*x3a2*x3a3,x3a1*x3a2*x3a3*x3a4
      print*,'- - - - - - - - - - - - - - - - - - - - - -'
      write(*,'(a,f9.3,$)')' 111',x1a1*x1a2*x1a3*x1b*x1c*x1d*x1e*xxx1
      write(*,'(10f9.3)')        x1a1,x1a2,x1a3,x1b,x1c,x1d,x1e,xxx1
      write(*,'(a,f9.3,$)')' 333',x3a1*x3a2*x3a3*x3a4*x3b*x3c*x3d*xxx3
      write(*,'(10f9.3)')        x3a1,x3a2,x3a3,x3a4,x3b,x3c,x3d,xxx3

      print*,'111',xxx1*x1a1*factk*.0389
     .*gamhad(iclpro)*gamhad(icltar)*xp**(-alppar)*xm**(-alppar)
     .*(1-xp)**alplea(iclpro)*(1-xm)**alplea(icltar)
     .,xp**(epsGp1)*xm**(epsGt1)
     .,1./ sigine * 10
     ., engy**epsG1*engy**epsG1       !missing in 3,  but value = 1

      print*,'333',xxx3*x3a1*.0389*factk
     .*gamhad(iclpro)*gamhad(icltar)*xp**(-alppar)*xm**(-alppar)
     .*(1-xm)**alplea(icltar)*(1-xp)**alplea(iclpro)
     .,xp**epsGp3*xm**epsGt3
     .,  b*.5*2.
     .,  1./(4.*pi)  *zb**(rp/rh)/rh/.0389
       !factor 1./ sigine * 10 included in ffom11a
      write(ifmt,'(14a)')('#####',i=1,14)
      stop'normal STOP in xEmsTest'
      end

