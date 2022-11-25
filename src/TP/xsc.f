C
C  This file is part of EPOS4
C  Copyright (C) 2022 research institutions and authors (See CREDITS file)
C  This file is distributed under the terms of the GNU General Public License version 3 or later
C  (See COPYING file for the text of the licence)
C

!###############################################################################
!###############################################################################
!           screening related stuff
!###############################################################################
!###############################################################################

cc-----------------------------------------------------------------------
c      subroutine CalcScrPair(b)
cckw NOT used any more, no real theoretical justification
cc-----------------------------------------------------------------------
cc  remove screened pairs from lproj3-kproj3-ltarg3-ktarg3-list
cc-----------------------------------------------------------------------
c
c#include "aaa.h"
c#include "ems.h"
c#include "par.h"
c      
c      bdummy=b !unused!!
c      do i=1,maproj
c        if(lproj3(i).gt.1)then
c          do l=1,lproj3(i)
c            nbproj3(i,l)=1
c          enddo
c        endif    
c      enddo          
c      do i=1,matarg
c        if(ltarg3(i).gt.1)then
c          do l=1,ltarg3(i)
c            nbtarg3(i,l)=1
c          enddo
c        endif    
c      enddo          
c
c      do i=1,maproj
c        if(lproj3(i).gt.1)then
c          do l=1,lproj3(i)
c          npair=1
c          do m=1,lproj3(i)
c            if(m.ne.l)then
cc update saturation radius according to number of connected pairs
c              bglaubxx=min(zbrmax,zbrads*sqrt(bglaubx
c     &                                   +b2xscr*log(float(npair)**2)))
c              kl=kproj3(i,l)
c              km=kproj3(i,m)
c              if((bkx(kl)-bkx(km))**2+(bky(kl)-bky(km))**2
c     .        .lt.bglaubxx**2)nbproj3(i,l)=nbproj3(i,l)+1
c              if((bkx(kl)-bkx(km))**2+(bky(kl)-bky(km))**2
c     .        .lt.zbrmax**2)npair=npair+1
c            endif
c          enddo
c          enddo
c        endif
c      enddo
c      do i=1,matarg
c        if(ltarg3(i).gt.1)then
c          do l=1,ltarg3(i)
c          npair=1
c          do m=1,ltarg3(i)
c            if(m.ne.l)then
cc update saturation radius according to number of connected pairs
c              bglaubxx=min(zbrmax,zbrads*sqrt(bglaubx
c     &                                   +b2xscr*log(float(npair)**2)))
c              kl=ktarg3(i,l)
c              km=ktarg3(i,m)
c              if((bkx(kl)-bkx(km))**2+(bky(kl)-bky(km))**2
c     .        .lt.bglaubxx**2)nbtarg3(i,l)=nbtarg3(i,l)+1
c              if((bkx(kl)-bkx(km))**2+(bky(kl)-bky(km))**2
c     .        .lt.zbrmax**2)npair=npair+1
c            endif
c          enddo
c          enddo
c        endif
c      enddo
c      
c      do i=1,maproj
c        if(lproj3(i).gt.1)then
c          do l=1,lproj3(i)
c            if(rangen().gt.1./float(nbproj3(i,l)))kproj3(i,l)=0
c          enddo
c        endif
c      enddo
c      do i=1,matarg
c        if(ltarg3(i).gt.1)then
c          do l=1,ltarg3(i)
c            if(rangen().gt.1./float(nbtarg3(i,l)))ktarg3(i,l)=0
c          enddo
c        endif
c      enddo
c
c      do i=1,maproj
c        if(lproj3(i).gt.1)then
c          ll=0
c          do l=1,lproj3(i)
c            if(kproj3(i,l).ne.0)then
c              ll=ll+1
c              kproj3(i,ll)=kproj3(i,l)
c            endif
c          enddo
c          lproj3(i)=ll
c        endif
c      enddo
c      do i=1,matarg
c        if(ltarg3(i).gt.1)then
c          ll=0
c          do l=1,ltarg3(i)
c            if(ktarg3(i,l).ne.0)then
c              ll=ll+1
c              ktarg3(i,ll)=ktarg3(i,l)
c            endif
c          enddo
c          ltarg3(i)=ll
c        endif
c      enddo
c
c      end
c

cc-----------------------------------------------------------------------
c      subroutine FuseEnhPom(k)
cc-----------------------------------------------------------------------
cc  Fuse some Pomeron together to mimic enhanced diagram behavior and
cc  and increase fluctuations for charm production and multiplicity
cc  distribution
cc-----------------------------------------------------------------------
c
c#include "aaa.h"
c#include "ems.h"
c#include "par.h"
c      double precision om1intgk,xp,xm,xps,xms,Penh,Pfus,GCorr,drangen
c
cc Probability to have enhanced diagram in this pair
c      ip=iproj(k)
c      it=itarg(k)
c      Penh=om1intgk(-k,xpp(ip),xmt(it))/om1intgk(k,xpp(ip),xmt(it))
c      if(ish.ge.4)write(ifch,*)'FuseEnhPom',Penh,nprmx(k)
c                
c      if(drangen(Penh).gt.Penh)then
cc Do some Pomeron fusion
c
cc start from a random position
c        n=int(rangen()*float(nprmx(k))) ! n-th spot for k-th pair
c        nn=1
cc Main loop
c        do while(nn.le.nprmx(k))
c          nn=nn+1
c          n=n+1
c          if(n.gt.nprmx(k))n=1
cc Pomeron to fuse with others
c          if(idpr(n,k).ne.0)then
c            xps=xppr(n,k)
c            xms=xmpr(n,k)
c            xpp(ip)=xpp(ip)+xps
c            xmt(it)=xmt(it)+xms
cc start again from a random position
c            n1=int(rangen()*float(nprmx(k))) ! n-th spot for k-th pair
c            nf=1
cc second loop
c            do while(nf.le.nprmx(k))
c              nf=nf+1
c              n1=n1+1
c              if(n1.gt.nprmx(k))n1=1
c              if(n1.ne.n.and.idpr(n1,k).ne.0.and.idfpr(n1,k).eq.0)then
cc candidate for fusion with selected Pom
c                xp=xps+xppr(n1,k)
c                xm=xms+xmpr(n1,k)
c                Pfus=1d0-Gcorr(xp,xm,bk(k),k,1.)/Penh
c                if(drangen(Pfus).lt.Pfus)then
cc fusion accepted
c                  xps=xp
c                  xms=xm
cc suppress merged Pom
c                  call RemPom(k,n1)
c                  if(ish.ge.4)write(ifch,*)'merged',Pfus,n1,'->',n
c               endif
c              endif
c            enddo
c            xpp(ip)=xpp(ip)-xps
c            xmt(it)=xmt(it)-xms
c            xppr(n,k)=xps
c            xmpr(n,k)=xms
c            idfpr(n,k)=-1      !not to reuse it (idfpr will be reinitialize in ProPoTy)
c          endif
c
c        enddo
c
c      endif
c
c      end
c

c-------------------------------------------------------------
      subroutine getZparEpsiK(zpav,ztav)
c-------------------------------------------------------------      
#include "aaa.h"
#include "sem.h"
#include "ems.h"
#include "par.h"
      common /cncl/xproj(mamx),yproj(mamx),zproj(mamx)
     *            ,xtarg(mamx),ytarg(mamx),ztarg(mamx)
      parameter(betxx=-0.9)


! Z_parton_projectile (zparpro) and Z_parton_target (zpartar)-----------

        call getZpar(zpav,ztav)

cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc


      
      ! calculation of epsilongp epsilongt
      
          !ip=0
       do k=1,koll
          !ipp=ip
        do i=-1,ntymx
          epsilongs(i,k)=0.
          epsilongp(i,k)=0.
          epsilongt(i,k)=0.
        enddo
        gammaV(k)=1.d0

        if(iscreen.eq.1)then

        ip=iproj(k)             !...........projectile
        if(lproj(ip).gt.0)then
          x0=zparpro(0,k)
          x1=zparpro(1,k)
          x2=zparpro(2,k)
          x3=zparpro(-1,k)
c          x0=x1              !to have same screening for diff and ND
          epsilongp(0,k)=abs(epscrd)*x0
          if(epscrd.le.0.)then
            epsilongs(0,k)=epsilongs(0,k)-epsilongp(0,k)
            betx=max(0.,betDp(0,iclpro,icltar)-alppar-betxx
     &          +min(0.,gamD(0,iclpro,icltar))*bk(k)**2)
            epsilongp(0,k)=0. !-min(betx,epsilongp(0,k))
          endif
          epsilongp(1,k)=x1!*abs(betDp(1,iclpro,icltar))
          if(epscrh.le.0.)then
            if(x1.gt.0.)
     &      epsilongs(1,k)=epsilongs(1,k)-epsilongp(1,k)*x2/x1
            betx=max(0.,betDp(1,iclpro,icltar)-alppar-betxx
     &          +min(0.,gamD(1,iclpro,icltar))*bk(k)**2)
            epsilongp(1,k)=-min(betx,epsilongp(1,k))
          endif
          epsilongp(-1,k)=abs(epscrh)*x3  !*abs(betDp(1,iclpro,icltar))
          if(epscrh.le.0.)then
            epsilongs(-1,k)=epsilongs(-1,k)-epsilongp(-1,k)
            betx=max(0.,betDp(1,iclpro,icltar)-alppar-betxx
     &          +min(0.,gamD(1,iclpro,icltar))*bk(k)**2)
            epsilongp(-1,k)=-min(betx,epsilongp(-1,k))
          endif
          if(ntymax.gt.ntynd)then
          epsilongp(2,k)=abs(epscrd)*x0
          epsilongp(3,k)=abs(epscrd)*x0
          if(epscrd.le.0.)then
            epsilongs(2,k)=epsilongs(2,k)-epsilongp(2,k)
            betx=max(0.,betDp(2,iclpro,icltar)-alppar-betxx
     &          +min(0.,gamD(2,iclpro,icltar))*bk(k)**2)
            epsilongp(2,k)=0. !-min(betx,epsilongp(2,k))
            epsilongs(3,k)=epsilongs(3,k)-epsilongp(3,k)
            betx=max(0.,betDp(3,iclpro,icltar)-alppar-betxx
     &          +min(0.,gamD(3,iclpro,icltar))*bk(k)**2)
            epsilongp(3,k)=0. !-min(betx,epsilongp(3,k))
          endif
          endif
ckw          gammaV(k)=1.d0
        endif
        it=itarg(k)             !...........target
        if(ltarg(it).gt.0)then
          x0=zpartar(0,k)
          x1=zpartar(1,k)
          x2=zpartar(2,k)
          x3=zpartar(-1,k)
c          x0=x1              !to have same screening for diff and ND
          epsilongt(0,k)=abs(epscrd)*x0
          if(epscrd.le.0.)then
            epsilongs(0,k)=epsilongs(0,k)-epsilongt(0,k)
            betx=max(0.,betDpp(0,iclpro,icltar)-alppar-betxx
     &          +min(0.,gamD(0,iclpro,icltar))*bk(k)**2)
            epsilongt(0,k)=0. !-min(betx,epsilongt(0,k))
          endif
          epsilongt(1,k)=x1!*abs(betDpp(1,iclpro,icltar))
          if(epscrh.le.0.)then
            if(x1.gt.0.)
     &      epsilongs(1,k)=epsilongs(1,k)-epsilongt(1,k)*x2/x1
            betx=max(0.,betDpp(1,iclpro,icltar)-alppar-betxx
     &          +min(0.,gamD(1,iclpro,icltar))*bk(k)**2)
            epsilongt(1,k)=-min(betx,epsilongt(1,k))
          endif
          epsilongt(-1,k)=abs(epscrh)*x3!*abs(betDp(1,iclpro,icltar))
          if(epscrh.le.0.)then
            epsilongs(-1,k)=epsilongs(-1,k)-epsilongt(-1,k)
            betx=max(0.,betDpp(1,iclpro,icltar)-alppar-betxx
     &          +min(0.,gamD(1,iclpro,icltar))*bk(k)**2)
            epsilongt(-1,k)=-min(betx,epsilongt(-1,k))
          endif
          if(ntymax.gt.ntynd)then
          epsilongt(2,k)=abs(epscrd)*x0
          epsilongt(3,k)=abs(epscrd)*x0
          if(epscrd.le.0.)then
            epsilongs(2,k)=epsilongs(2,k)-epsilongt(2,k)
            betx=max(0.,betDpp(2,iclpro,icltar)-alppar-betxx
     &          +min(0.,gamD(2,iclpro,icltar))*bk(k)**2)
            epsilongt(2,k)=0. !-min(betx,epsilongt(2,k))
            epsilongs(3,k)=epsilongs(3,k)-epsilongt(3,k)
            betx=max(0.,betDpp(3,iclpro,icltar)-alppar-betxx
     &          +min(0.,gamD(3,iclpro,icltar))*bk(k)**2)
            epsilongt(3,k)=0. !-min(betx,epsilongt(3,k))
          endif
          endif
ckw          gammaV(k)=1.d0
        endif
        do i=-1,ntymax
          if(i.eq.0)then
            if(epscrd.le.0.)epsilongs(i,k)=0.5*epsilongs(i,k)
          elseif(abs(i).eq.1)then
            if(epscrh.le.0.)epsilongs(i,k)=0.5*epsilongs(i,k)
          else
            if(epscrd.le.0.)epsilongs(i,k)=0.5*epsilongs(i,k)
          endif
        enddo

        if(iomega.ge.2)then
          epsilongs(0,k)=epsilongs(1,k)
          epsilongp(0,k)=epsilongp(1,k)
          epsilongt(0,k)=epsilongt(1,k)
        endif

        endif                     ! no screening

      enddo


      end        

c-------------------------------------------------------------
      subroutine xEpsilonPrep(irea,epsG,epsGp,epsGt,gamb,zpav,ztav)
c-------------------------------------------------------------      
#include "aaa.h"
#include "sem.h"
#include "ems.h"
#include "par.h"
      parameter(nxeps=22,nyeps=32)
      common/cxeps1/w(0:nxeps,nyeps),y1(nyeps),y2(nyeps)
      common/cxeps2/dbp,b1p,b2p
      common/nucl3/phi,bimp
      
      do k=1,koll
        do i=idxDmin(iomega),idxDmax(iomega)
      
   !...........prepare plot in xepsilon
          if(irea.eq.0)then
           kk=min(nyeps,max(1,int((bk(k)-b1p)/dbp)+1))
           if(i.lt.2)then
             if(i.eq.1)w(0,kk)=w(0,kk)+1.
             if(i.eq.1)w(1,kk)=w(1,kk)+abs(epsGp)
             if(i.eq.1)w(2,kk)=w(2,kk)+abs(epsGt)
c                     w(3+i,kk)=w(3+i,kk)+abs(epsG)
         !...5-8 soft ... 9-12 semi
                     w(5+4*i,kk)=w(5+4*i,kk)
     *              +betDp(i,iclpro,icltar)   !prj eff
     *              +epsGp+gamb
             w(6+4*i,kk)=w(6+4*i,kk)
     *              +betDpp(i,iclpro,icltar)  !tgt eff
     *              +epsGt+gamb
             w(7+4*i,kk)=w(7+4*i,kk)
     *              +betDp(i,iclpro,icltar)  !prj unscr
     *              +gamb
             w(8+4*i,kk)=w(8+4*i,kk)
     *              +betDpp(i,iclpro,icltar) !tgt unscr
     *              +gamb
             if(i.eq.0)w(13,kk)=w(13,kk)+zparpro(1,k)
             if(i.eq.0)w(14,kk)=w(14,kk)+zpartar(1,k)
c             if(i.eq.0)w(15,kk)=w(15,kk)+zparpnx(k)
c             if(i.eq.0)w(16,kk)=w(16,kk)+zpartnx(k)
           else
             if(epscrd.lt.0.)then
               w(3,kk)=w(3,kk)+epsG
               w(4,kk)=0.
             else
               w(3,kk)=w(3,kk)+epsGp
               w(4,kk)=w(4,kk)+epsGt
             endif
           endif
          endif
   !................
        enddo
      enddo
      
           zppevt=zpav/koll
           zptevt=ztav/koll
           if(irea.eq.0)then
             ktot=int(2.*bimp)+1
             if(ktot.le.nyeps)then
               w(17,ktot)=w(17,ktot)+zppevt
               w(18,ktot)=w(18,ktot)+zptevt
               w(19,ktot)=w(19,ktot)+1
             endif
             n=1+int(float(nglevt)/(0.1*maproj*matarg))*(nyeps-1)
             if(nglevt.ge.1.and.n.ge.1.and.n.le.nyeps)then
               w(20,n)=w(20,n)+zppevt
               w(21,n)=w(21,n)+zptevt
               w(22,n)=w(22,n)+1
             endif
           endif
      
      
      end
      
c-------------------------------------------------------------
      subroutine printGfunParK
c-------------------------------------------------------------      
#include "aaa.h"
#include "sem.h"
#include "ems.h"
#include "par.h"
      double precision atildg,btildgp,btildgpp
      common/cgtilde/atildg(idxD0:idxD1,kollmx)
     *,btildgp(idxD0:idxD1,kollmx),btildgpp(idxD0:idxD1,kollmx)

      if(ish.ge.10)then
      do k=1,koll
        ip=iproj(k)
        it=itarg(k)
        write(ifch,*)' k,b,ip,it,',k,bk(k),ip,it
        write(ifch,*)' zparpro,zpartar'
     *      ,zparpro(0,k),zparpro(1,k),zpartar(0,k),zpartar(1,k)
        write(ifch,*)' gammaV,epsilonGP0-2,epsilonGT0-2,epsilonGs0-2'
     *   ,gammaV(k)
     *   ,epsilongp(0,k),epsilongp(1,k),epsilongp(2,k),epsilongp(3,k)
     *   ,epsilongt(0,k),epsilongt(1,k),epsilongt(2,k),epsilongt(3,k)
     *   ,epsilongs(0,k),epsilongs(1,k),epsilongs(2,k),epsilongs(3,k)
        write(ifch,*)'*******************************************'
c        write(ifch,*)" atilde,btildep,btildepp "
c        do i=ntymin,ntymax
c        write(ifch,*)i,atilde(i,k),btildep(i,k),btildepp(i,k)
c        enddo
c        write(ifch,*)" atildg,btildgp,btildgpp "
c        do i=ntymin,ntymax
c        write(ifch,*)i,atildg(i,k),btildgp(i,k),btildgpp(i,k)
c        enddo
        sy=engy*engy
        zp=0.
        zt=0.
        call GfunPar(zp,zt,0.,0.,1,0,bk(k),sy,alp,bet,betp
     &                  ,epsp,epst,epss,gamvv)
        zp=0.
        zt=0.
        call GfunPar(zp,zt,0.,0.,1,1,bk(k),sy,alp,bet,betp
     &                  ,epsp,epst,epss,gamvv)
      enddo
      endif

      end
      
c-------------------------------------------------------------
      subroutine getZparEpsi(m,ii,b,zzip,zzit,eflp,eflt
     .                      ,epsG,epsGp,epsGt)
c-------------------------------------------------------------   
ctp 21.09.2017 Metropolis works very well, so if some difference
ctp appears between MC and theory (for Npom, profile function, 
ctp cross section) this is probably due to a difference in screening
ctp definition between here (theory) and Zpair/getZpar (MC)
c-------------------------------------------------------------   
#include "aaa.h"
#include "sem.h"
#include "par.h"
#include "ems.h"
      parameter(eps=1.d-20,betxx=-0.9)

      i=abs(ii)
      epsG=0.
      epsGp=0.
      epsGt=0.
      zzp=0.
      zzt=0.
      zb=0.
      zb0=0.
      zzz=0.
      bc=0.
      bx=0.
      if(iscreen.ne.0.and.isetcs.ge.-1)then
        if(ish.ge.10)write(ifch,*)'Eps in:',m,i,zzip,zzit
        b2=b*b
        cfalpro=ucfpro
        cfaltar=ucftar
        gamb=gamD(i,iclpro,icltar)*b2
c        b2x=bglaub2
c        b2x=b2x*epscrp
        b2x=b2xscr
c        b2x=epscrp**2           !max(0.5*b2x,b2x+epscrp*fegy)
        if(i.eq.1.or.iomega.ge.2)then
          epsmx=epscrx
          if(ii.gt.0.and.m.gt.0)then
            b2x=2.*b2x
            bc=sign(1.,epsmx)
            bx=sqrt(b2x)
          else !for Q2s
            b2x=2.*b2x  !bsatur
          endif
          b2y=0.
        else
          b2x=epscrp**2         !max(0.5*b2x,b2x+epscrp*fegy)
          epsmx=0. !epscrs
c          b2y=exp(-epscrb/(slopom*log(max(1.,engy*engy))))
        endif
        epsmx=abs(epsmx)
        if(abs(bc).gt.0.)bc=bc*zbcutx !bcutz(fegy,b2x,epsmx)
          zzz=zzip+zzit
          dzzz=0.
          if(zzz.gt.0.)dzzz=abs(zzip-zzit)
        do iii=1,2
        if(iii.eq.1)then
          eflu=eflp
        else
          eflu=eflt
        endif
c        dzzz=0.
c        if(zzz.gt.0.)dzzz=abs(zzip-zzit)
c        zzz=zzz+dzzz
c        if(m.lt.0)dzzz=0.
c        b2x=b2xscr
        fmn=0.
        if(i.eq.1.or.iomega.ge.2)then
c          fmn=1.
          epsmx=abs(epscrx)+znurho*zzz
          fegy=epsmx            !max(epsmx,epscrw) !+zrhoinc*zzz)
c          if(ii.gt.0)then
c            fegy=epsmx
c          else
c            fegy=1.
c          endif
        else
          fegy=fegypp
c          fegy=max(epscrw,fegy)
          b2y=exp(-epscrb/fegy)
        endif
c        b2b=max(min(0.,bc),abs(b)-abs(bc))
c        b2a=max(0.,b2b)**2
        b2a=max(0.,abs(b)-bx)**2/b2x!*4.
        b2c=b**2/b2xscr!/epscrs
        zb01=0.
        zb02=0.
        if(i.eq.1.or.iomega.ge.2)then
c          zb=fegy/(1.+epscrw*abs(b))
c          b2b=max(0.,abs(b)-abs(bc))**2/(bx-abs(bc))**2*5.
          zb=fegy!+epscrw*exp(-b2b)
          if(bc.lt.0)then
            b2b=(abs(b)/abs(bc))**2
            zb00=exp(-b2b/(bc*bc)*4.)
            zb01=epsmx*epscrg+zrhoinc*dzzz-eflu    !change slope
            zb02=epsmx*epscrs+zrhoinc*dzzz-eflu         !change normalization
            zb01=zb01*zb00
            zb02=zb02*zb00
          endif
c        if(abs(bc).gt.0.)then
c          fegy=epsmx
c        else
c          fegy=max(epsmx,fegy)
c        endif
          if(abs(b).ge.bx)zb=zb*exp(-b2a)
        else
          b2b=max(0.,abs(b)+0.5*sqrt(b2xscr))**2/b2x*4. !b2/b2x*5.
c          b2b=b2/b2x*5.
c        if(abs(b).lt.abs(bc))then
c          fmn=1.
c          if(i.eq.1)fmn=epscrs          !*max(0.,1.-zrhoinc*zzz)
c        elseif(abs(bc).gt.0)then
c          b2a=min(b2a,sqrt(b2a))
c        else
c         do not use a gaussian but an exponential to have smooth profile function
c         b2a=b2a**epscrb
          if(b2y.gt.0.)then
c            b2b=b2b/b2y
            b2b=min(b2b,b2b**b2y)
          else
            b2b=0.
          endif
c          b2a=min(b2a/b2x,(b2a/epscrw)**epscrb)
c        endif
c        zb=fegy*max(fmn,exp(-(1.-fmn**0.5)*b2a))
          zb=fegy*exp(-b2b)
c          if(abs(b).ge.sqrt(b2x))zb=0.
          zb=zb+zzz
        endif
c        if(abs(b).ge.sqrt(b2x))zb=zb*exp(-b2a)
c        if(i.eq.1.and.abs(b).ge.abs(bc))then
c          zb=min(zb+zzz,epsmx)
c        endif
        zz=max(0.,zb+zb01)
        zz0=max(0.,zb+zb02)
         
        if(iii.eq.1)then
          zzp=zz
          zzp0=zz0
        else
          zzt=zz
          zzt0=zz0
        endif
        enddo

        x=zzp
        x2=zzp0
        epsG=0.
        betx=max(0.,betDp(i,iclpro,icltar)-alppar-betxx+min(0.,gamb))
        if(ii.eq.-1)then   !screening for hard component (no saturation)
          epsGp=abs(epscrh)*x
          if(epscrh.le.0.)then
            epsG=epsG-epsGp
            epsGp=-min(betx,epsGp)
          endif
        elseif((i.gt.idxD0.and.i.le.idxD).or.iomega.ge.2)then
          epsGp=x!*abs(betDp(i,iclpro,icltar))
          if(epscrh.le.0.)then
            if(x.gt.0.)epsG=epsG-epsGp*x2/x
            epsGp=-min(betx,epsGp)
          endif
        elseif(i.gt.idxD)then
          epsGp=abs(epscrd)*x
          if(epscrd.le.0.)then
            epsG=epsG-epsGp
            epsGp=0.!-min(betx,epsGp)
          endif
        else
          epsGp=abs(epscrd)*x
          if(epscrd.lt.0)then
            epsG=epsG-epsGp
            epsGp=0.!-min(betx,epsGp)
          endif
        endif
c        if(abs(epsGp/betx+1.).le.1e-5)print *,'pro',i,epsGp,x
        gamV=1.
        x=zzt
        x2=zzt0
        betx=max(0.,betDpp(i,iclpro,icltar)-alppar-betxx+min(0.,gamb))
        if(ii.eq.-1)then
          epsGt=abs(epscrh)*x
          if(epscrh.le.0.)then
            epsG=epsG-epsGt
            epsGt=-min(betx,epsGt)
          endif
        elseif((i.gt.idxD0.and.i.le.idxD).or.iomega.ge.2)then
          epsGt=x!*abs(betDpp(i,iclpro,icltar))
          if(epscrh.le.0)then
            if(x.gt.0.)epsG=epsG-epsGt*x2/x
            epsGt=-min(betx,epsGt)
          endif
        elseif(i.gt.idxD)then
          epsGt=abs(epscrd)*x
          if(epscrd.le.0)then
            epsG=epsG-epsGt
            epsGt=0. !-min(betx,epsGt)
          endif
        else
          epsGt=abs(epscrd)*x
          if(epscrd.lt.0)then
            epsG=epsG-epsGt
            epsGt=0. !-min(betx,epsGt)
          endif
        endif
c        if(abs(epsGt/betx+1.).le.1e-5)print *,'tar',i,epsGt,x

        epsG=0.5*epsG

c        gamV=gamV
      endif

      end   
c----------------------------------------------------------------------
      double precision function om51k(n,k,iqq)   !---MC---
c----------------------------------------------------------------------
c     om1 * 0.5 (real with screening corrections)
c                       iqq=0 .... soft
c n  - pomeron index    iqq=1 .... gg
c k  - pair index       iqq=2 .... qg
c                       iqq=3 .... gq
c                       iqq=4 .... qq
c----------------------------------------------------------------------
#include "aaa.h"
#include "ems.h"
      double precision xp,xm,xh,yp,om51p,epsG,epsGp,epsGt,Fgfactork
      double precision plc,s
      common/cems5/plc,s

      xh=xpr(n,k)
      yp=ypr(n,k)
      sy=sngl(s*xh)

      om51k=0.d0
      if(xh.le.0.d0.or.xh.gt.1.d0)return
      xp=xppr(n,k)
      xm=xmpr(n,k)
      b=bk(k) !bhpr(n,k)

      i=min(iqq,1)
      epsG=epsilongs(i,k)
      epsGp=epsilongp(i,k)
      epsGt=epsilongt(i,k)
      if(iqq.eq.0)then
        om51k=om51p(sy,xh,yp,b,iqq)
     &       *xp**epsGp*xm**epsGt*s**epsG
      elseif(iqq.eq.1)then
        om51k=om51p(sy,xh,yp,b,iqq)
     &       *xp**epsGp*xm**epsGt*Fgfactork(n,k,xp,xm)*s**epsG
      elseif(iqq.eq.2)then
        om51k=om51p(sy,xh,yp,b,iqq)
     &       *xm**epsGt*Fgfactork(n,k,1d0,xm)*s**epsG
      elseif(iqq.eq.3)then
        om51k=om51p(sy,xh,yp,b,iqq)
     &       *xp**epsGp*Fgfactork(n,k,xp,1d0)*s**epsG
      elseif(iqq.eq.4)then
        om51k=om51p(sy,xh,yp,b,iqq)
      else
        call utstop('Wrong iqq in om51k&')
      endif
      end

     
c----------------------------------------------------------------------
      function fscra(x)
c----------------------------------------------------------------------
#include "par.h"
c      fscra=egyscr+x**epscrw !0.
      fscra=0.
      x2=(x/egyscr)!**0.5
      if(epscrw+x2.gt.1)fscra=log(epscrw+x2)!**epscrs!**2      

      end

c----------------------------------------------------------------------
      function fscro(x,rho)
c----------------------------------------------------------------------
#include "aaa.h"
#include "par.h"
      dummy=rho
      fscro=0.
      x2=x/egyscr
c      if(x2.gt.1.)fscro=sqrt(log(x2)**2+fscro**2)
      if(epscrw+x2.gt.1.)fscro=log(epscrw+x2)
c      fscro=epscrw*fscro!+znurho*rho 
      end

