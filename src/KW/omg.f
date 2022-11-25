C
C  This file is part of EPOS4
C  Copyright (C) 2022 research institutions and authors (See CREDITS file)
C  This file is distributed under the terms of the GNU General Public License version 3 or later
C  (See COPYING file for the text of the licence)
C

c         ------------------------------------------------
c         since EPOS3086 we have om5 = om1
c                  (in earlier versions   om5 = 0.5 * om1)
c         change in : om51nf,om5Jk,om5J,om1f,om1fJ,omIgamint,omGam,omGamint
c         ------------------------------------------------
c----------------------------------------------------------------------
c
c            parameterized "G" -- in the program called "om..."
c
c----------------------------------------------------------------------
c The two elementary fuctions of our approach, the profile fuction G
c and the Phi exponent G~, are here referred to as Gpro and Gtilde.
c Both functions can be written as
c
c               G = \sum_type  alp * xp**bet * xm**betp
c
c The parameters alp, bet, betp are calculated in GfunParK (k-mode,
c b is takento the one of pair k) or GfunPar (b-mode: arbitrary b) as
c
c  Gpro:   bet  = betD'  + epsGp + gamD*b**2 + epsG -alppar
c          bet' = betD'' + epsGt + gamD*b**2 + epsG -alppar
c
c          alp  = alpD*f * (s/s0)**(epsG) * exp(-b**2/delD)
c
c  Gtilde: bet~  = bet  + 1
c          bet~' = bet' + 1
c
c          alp~  = alp * gam(bet~)          * gam(bet~')
c                      * gam(1+alppro)      * gam(1+alptar)
c                      * gam(1+alppro+bet~) * gam(1+alptar+bet~')
c                      * (1+epsGt')         * (1+epsGt')
c
c The parameters depend implicitely on s.
c
c In the program we use om1 = Gpro
c  (they differ by a constant which is actually one)
c and om5 = om1
c All functions related to om1 are called om1... .
c
c The inclusive Pomeron distributions are
c
c      PomInc...(xp,xm) = Gpro(xp,xm) * (1-xp)**alppro * (1-xm)**alptar
c
c The functional 
c
c      Pominc(fct,...) 
c
c is used to parameterize the deformation due to multiple scattering,
c called with fct=PomInc...
c
c----------------------------------------------------------------------


c----------------------------------------------------------------------
      subroutine GfunParK(irea)   !---MC---
c----------------------------------------------------------------------
c  calculates parameters alp,bet,betp of the G functions (k-mode)
c  and the screening exponents epsilongp(i,k), epsilongt(i,k), epsilongs(i,k)
c----------------------------------------------------------------------
c  Gpro parameters written to /comtilde/atilde(,)btildep(,),btildepp(,)
c Gtilde parameters written to /cgtilde/atildg(,)btildgp(,),btildgpp(,)
c  two subscripts: first=type, second=collision k
c Certain pieces of this routine are only done if irea is <= or = zero.
c----------------------------------------------------------------------
#include "aaa.h"
#include "sem.h"
#include "ems.h"
#include "par.h"
      double precision atildg,btildgp,btildgpp
      common/cgtilde/atildg(idxD0:idxD1,kollmx)
     *,btildgp(idxD0:idxD1,kollmx),btildgpp(idxD0:idxD1,kollmx)
      double precision utgam2,coefgdp,coefgdt
      parameter(nbkbin=100)
      common /kfitd/ xkappafit(nbkbin,nclegy,nclha,nclha),xkappa,bkbin
      common /cgtilnu/ cfbetpnp,cfbetppnp,cfbetpnm,cfbetppnm,cfalpro
     &,cfaltar,cfbpap,cfbpam,cfbppap,cfbppam
      double precision cfbetpnp,cfbetppnp,cfbetpnm,cfbetppnm,cfalpro
     &,cfaltar,cfbpap,cfbpam,cfbppap,cfbppam,gamv,eps,sy,om1mean!,om1min!,SigGFFk,SigGFFc
      parameter (eps=1.d-25,betxx=-0.9)
      parameter(nxeps=22,nyeps=32)
      common/cxeps1/w(0:nxeps,nyeps),y1(nyeps),y2(nyeps)
      common/cxeps2/dbp,b1p,b2p
      common/geom/rmproj,rmtarg,bmax,bkmx
      common/nucl3/phi,bimp
      common /cncl/xproj(mamx),yproj(mamx),zproj(mamx)
     *            ,xtarg(mamx),ytarg(mamx),ztarg(mamx)
      common/ckopjtg/kopj(mamx),kotg(mamx)
      common/cnglndi/nglndik
      double precision om1intgk,dbetsf,dbetpsf,dalp,xmn,salp
      logical test
      real segy
      data segy/0./
      save 
c      save SigGFFk,SigGFFc

      b1p=0.03
      b2p=bkmx*1.2
      dbp=(b2p-b1p)/nyeps
      call utprj('GfunParK ',ish,ishini,10)

      cfbetpnp=0.d0
      cfbetppnp=0.d0
      cfbetpnm=0.d0
      cfbetppnm=0.d0
      cfalpro=dble(ucfpro)
      cfaltar=dble(ucftar)
      cfbpap=1.d0
      cfbppap=1.d0
      cfbpam=1.d0
      cfbppam=1.d0
      om1mean=0d0
c      om1min=1d30
      epsG=0.
      eGpsf=0.
      eGtsf=0.
      eGsf=0.
      gbsf=0.
      sy=dble(engy*engy)
      xmn=1d0/sy
      nglndik=0
      test=.false.

      do k=1,koll
        do i=ntymi,ntymx
          atilde(i,k)=0.d0
          btildep(i,k)=0.d0
          btildepp(i,k)=0.d0
        enddo
        do i=idxD0,idxD1
          atildg(i,k)=0.d0
          btildgp(i,k)=0.d0
          btildgpp(i,k)=0.d0
        enddo
      enddo

      do k=1,koll
       do n=1,nprmx(k)
        do i=ntymi,ntymx
          atildp(i,n,k)=0.d0
          btildpp(i,n,k)=0.d0
          btildppp(i,n,k)=0.d0
        enddo
       enddo
      enddo

      ntymin=ntymi           !minimum indice for diagram type
      ntymax=ntymx           !maximum indice for diagram type
      if(iomega.ge.2)then    !no diffraction
        ntymin=ntynd
        ntymax=ntynd
      endif
      
c initialize screening corrections
      zpav=0.
      ztav=0.
      call getZparEpsiK(zpav,ztav)

      !########
      !########   Important
      !########   ---------
      !########
      !########  don't change the following loop structure 
      !########  (do kkk ... ; do k ... ; if(kkk... ; if(kkk.eq.2 ...  )
      !########   needed for an efficient saturation scale  interpolation
      !########   via "call ipo..."                                     
      !########    
      !########             the saturation scale code 
      !########                      is marked as      !### kw ###
      !########

ctp Remark on "important" : true if only small part affected by Q2s
ctp                         but in fact Q2s used most of the time
ctp                         so we don't lose much time setting it always
ctp                         and it avoids some mistake and the dependence
ctp                         on glauber defined pairs which are not
ctp                         good enough for Q2s
      
c      do kkk=1,2                                         !### kw ###
      do k=1,koll                                        !### kw ###
      ip=iproj(k)                                        !### kw ###
      it=itarg(k)                                        !### kw ###
      !print*,'xxxxxxxxx',ip,kopj(ip),it,kotg(it)        !### kw ###
c      if(  kkk.eq.1.and.kopj(ip)*kotg(it).eq.0           !### kw ###
c     ..or. kkk.eq.2.and.kopj(ip)*kotg(it).ne.0)then      !### kw ###

        !------------------------------------------------------
        !       alpha beta betap for   Gtilde    (->PhiExpo) 
        !------------------------------------------------------
 
        b=bk(k)
        b2=bk(k)*bk(k)
        b2a=bk(k)*bk(k)
        ip=iproj(k)
        it=itarg(k)

        if(b.lt.(nbkbin-1)*bkbin)then
          ibk=int(bk(k)/bkbin)+1
          if(isetcs.gt.1.and.iclegy.lt.iclegy2)then
            egy0=egylow*egyfac**float(iclegy-1)
            xkappa1=xkappafit(ibk,iclegy,iclpro,icltar)
     *         +(bk(k)-bkbin*float(ibk-1))/bkbin
     *         *(xkappafit(ibk+1,iclegy,iclpro,icltar)
     *         -xkappafit(ibk,iclegy,iclpro,icltar))
            xkappa2=xkappafit(ibk,iclegy+1,iclpro,icltar)
     *         +(bk(k)-bkbin*float(ibk-1))/bkbin
     *         *(xkappafit(ibk+1,iclegy+1,iclpro,icltar)
     *         -xkappafit(ibk,iclegy+1,iclpro,icltar))
            xkappa=xkappa1+(xkappa2-xkappa1)/log(egyfac)
     *         *(log(engy)-log(egy0))
          else
            xkappa=xkappafit(ibk,iclegy,iclpro,icltar)
     *         +(bk(k)-bkbin*float(ibk-1))/bkbin
     *         *(xkappafit(ibk+1,iclegy,iclpro,icltar)
     *         -xkappafit(ibk,iclegy,iclpro,icltar))
          endif
        else
          xkappa=1.
        endif

        gamV=gammaV(k)
        do i=idxDmin(iomega),idxDmax(iomega)
c            epsG=0.
c            epsGp=0.
c            epsGt=0.
          epsG =epsilongs(i,k)
          epsGp=epsilongp(i,k)
          epsGt=epsilongt(i,k)
          gamb=gamD(i,iclpro,icltar)*b2
          if(i.eq.0.and.alpDp(0,iclpro,icltar).gt.0)then
            eGsf =epsilongs(1,k)
            eGpsf=epsilongp(1,k)
            eGtsf=epsilongt(1,k)
            gbsf=alpDp(2,iclpro,icltar)*b2
c for soft+diff component, exponent should not be higher than soft+screening one
            dbetsf=dble(max(betxx,alpDp(1,iclpro,icltar)
     *         +eGpsf
     *         +gbsf-alppar))
            dbetpsf=dble(max(betxx,alpDp(1,iclpro,icltar)
     *         +eGtsf
     *         +gbsf-alppar))
          else
            dbetsf=1d35
            dbetpsf=1d35
          endif
          atildg(i,k)=cfalpro*cfaltar
     *            *gamv
          atildg(i,k)=atildg(i,k)
     *            *dble(xkappa*xkappa)
c          if(i.le.idxD)then
c            epsGp=epsilongp(i,k)
c            epsGt=epsilongt(i,k)
            btildgp(i,k)=min(dbetsf,dble(max(betxx,
     *                    betDp(i,iclpro,icltar)+epsGp
     *                    +gamb-alppar)))
            btildgpp(i,k)=min(dbetpsf,dble(max(betxx,
     *                    betDpp(i,iclpro,icltar)+epsGt
     *                    +gamb-alppar)))
            if(i.eq.0.and.alpDp(0,iclpro,icltar).gt.0)then
c save proper eps
            if(abs(btildgp(i,k)-dbetsf).lt.1d-10)epsGp=eGpsf
            if(abs(btildgpp(i,k)-dbetpsf).lt.1d-10)epsGt=eGtsf
c for soft+diff component, factor should not be lower than soft+screening one
              salp=dble(alpDp(0,iclpro,icltar))
     *       *exp(-dble(b2)/dble(alpDp(3,iclpro,icltar)))
     *              /xminDf**dble(eGsf+gbsf)
              salp=salp
     *      *xmn**(0.5d0*(dbetsf+dbetpsf-btildgp(i,k)-btildgpp(i,k)))
            else
              salp=0d0
            endif
            dalp=dble(alpD(i,iclpro,icltar))
     *         *exp(-dble(b2)/dble(delD(i,iclpro,icltar)))
     *            /xminDf**dble(epsG+gamb) !*dble(sy)**(betD(i,iclpro,icltar)+gamb+epsG)!ctp20160906 s-independant fit now
            if(dalp.gt.salp)then
              atildg(i,k)=atildg(i,k)*dalp
            else
              atildg(i,k)=atildg(i,k)*salp
              epsG=eGsf
            endif
            btildgp(i,k)=btildgp(i,k)+1d0
            btildgpp(i,k)=btildgpp(i,k)+1d0
c          else
c            absb=abs(bk(k))-bmxdif(iclpro,icltar)
c            b2a=absb*absb
c            atildg(i,k)=atildg(i,k)
c     *                 *sy**dble(betD(i,iclpro,icltar)+epsG)
c     *                 *dble(exp(-b2a/delD(i,iclpro,icltar)))
c            btildgp(i,k)=                                     !WWWdiffr
c     *      dble(betDp(i,iclpro,icltar)-alppar+epsGp)+1.d0    !WWWdiffr
c            btildgpp(i,k)=                                    !WWWdiffr
c     *      dble(betDpp(i,iclpro,icltar)-alppar+epsGt)+1.d0   !WWWdiffr
c          endif
          coefgdp=utgam2(1.d0+dble(alplea(iclpro))+btildgp(i,k))
          coefgdt=utgam2(1.d0+dble(alplea(icltar))+btildgpp(i,k))
          atildg(i,k)=atildg(i,k)
     *            *utgam2(btildgp(i,k))*utgam2(btildgpp(i,k))
     *            /coefgdp/coefgdt
          call xEpsilonPrep(irea,epsG,epsGp,epsGt,gamb,zpav,ztav)
        enddo

        !------------------------------------------------------
        !         alpha beta betap for   Gpro
        !------------------------------------------------------
        if(irea.le.0)then
        ip=iproj(k)
        it=itarg(k)
        b2=bk(k)*bk(k)
        do i=ntymin,ntymax
c            epsG=0.
c            epsGp=0.
c            epsGt=0.
          epsG =epsilongs(i,k)
          epsGp=epsilongp(i,k)
          epsGt=epsilongt(i,k)
          gamb=gamD(i,iclpro,icltar)*b2
          if(i.eq.0.and.alpDp(0,iclpro,icltar).gt.0)then
            eGsf =epsilongs(1,k)
            eGpsf=epsilongp(1,k)
            eGtsf=epsilongt(1,k)
            gbsf=alpDp(2,iclpro,icltar)*b2
c for soft+diff component, exponent should not be higher than soft+screening one
            dbetsf=dble(max(betxx,alpDp(1,iclpro,icltar)
     *         +eGpsf
     *         +gbsf-alppar))
            dbetpsf=dble(max(betxx,alpDp(1,iclpro,icltar)
     *         +eGtsf
     *         +gbsf-alppar))
          else
            dbetsf=1d35
            dbetpsf=1d35
          endif
          atilde(i,k)=dble(alpD(i,iclpro,icltar)) *gamv  !ckw
c          if(i.le.ntynd)then
c            epsGp=epsilongp(i,k)
c            epsGt=epsilongt(i,k)
            btildep(i,k)=min(dbetsf,dble(max(betxx,
     *                   betDp(i,iclpro,icltar)+epsGp + gamb - alppar)))
            btildepp(i,k)=min(dbetsf,dble(max(betxx,
     *                  betDpp(i,iclpro,icltar)+epsGt + gamb - alppar)))
            if(i.eq.0.and.alpDp(0,iclpro,icltar).gt.0)then
c for soft+diff component, factor should not be lower than soft+screening one
              salp=dble(alpDp(0,iclpro,icltar)) *gamv  !ckw
     *       *exp(-dble(b2)/dble(alpDp(3,iclpro,icltar)))
     *              /xminDf**dble(eGsf+gbsf)
              salp=salp
     *      *xmn**(0.5d0*(dbetsf+dbetpsf-btildep(i,k)-btildepp(i,k)))
            else
              salp=0d0
            endif
            atilde(i,k)=max(salp,atilde(i,k)
     *              *exp(-dble(b2)/dble(delD(i,iclpro,icltar)))
     *              /xminDf**dble(epsG+gamb)) !*dble(sy)**(betD(i,iclpro,icltar)+gamb+epsG))!ctp20160906 s-independant fit now
c          else
c            absb=abs(bk(k))-bmxdif(iclpro,icltar)
c            b2a=absb*absb
c            atilde(i,k)=atilde(i,k)
c     *                            *sy**dble(betD(i,iclpro,icltar)+epsG)
c     *                            *dble(exp(-b2a/delD(i,iclpro,icltar)))
c
c            btildep(i,k)=dble(betDp(i,iclpro,icltar)-alppar+epsGp)     !WWWdiffr
c            btildepp(i,k)=dble(betDpp(i,iclpro,icltar)-alppar+epsGt)   !WWWdiffr
c          endif
          test=btildep(i,k)+1d0.lt.-eps.or.btildepp(i,k)+1d0.lt.-eps
          if(test)goto 999
        enddo
       endif

        !---------------------------------------------------------------
        !    former integom1
        !----------------------------------------------------------------

        om1intc(k)=om1intgk(k,1d0,1d0)

        !---------------------------------------------------------------
        !    mean number of collision
        ! 
        ! In pA or AA It is possible to use the <Npom> for each given
        ! pair with SigGFFk
        !----------------------------------------------------------------

c        om1MeanCol(1,k)=dble(xnpomcol)
c       if(k.eq.1)call fsigmak(k)  !set sigmak and om1MeanCol
        if(iomega.eq.2)then
          if(abs(segy-engy).gt.1e-5.or.(maproj.gt.1.or.matarg.gt.1))then
            segy=engy
c          if(bsatur.gt.0.)then
cc        if(nrevt.eq.0.and.k.eq.1)then
c            call NhpomQ2s(k,SigGFFk,SigGFFc)
cc        endif
c            om1MeanCol(1,k)=SigGFFk
c            om1MeanCol(2,k)=SigGFFc
c          else
c            om1MeanCol(1,k)=dble(xnhardcol)
c           om1MeanCol(2,k)=dble(xnsatcol)xnpomcol
            call fsigmak(k)  !set sigmak and om1MeanCol(1)
c            om1MeanCol(1,k)=dble(xnpomcol)   !overwrite to have fixed pp value
            om1MeanCol(1,k)=dble(xnhardcol)   !overwrite to have fixed pp value
c           om1MeanCol(1,k)=0.5d0*dble(xnpomcol+xnhardcol)   !overwrite to have fixed pp value
            if(k.eq.1)then
              om1MeanCol(2,k)=om1MeanCol(1,k)**dble(1.+vparam(3))
     .                  *(1d0+om1MeanCol(1,k)*log(dble(maproj+matarg-1))
     .             *dble(1.+vparam(1)))
            else
              om1MeanCol(2,k)=om1MeanCol(2,1)
            endif
c          endif
c            om1min=min(om1min,om1MeanCol(1,k))
c            om1mean=om1mean+om1MeanCol(1,k)
c            om1MeanCol(1,k)=om1MeanCol(1,1)
c            om1MeanCol(2,k)=om1MeanCol(1,k)   !for koll=1
c            sigmak(1,k)=sigmak(1,1)
c            sigmak(2,k)=sigmak(2,1)
            bglndif=max(0.,sngl(sigmak(2,k))/10./pi) !10= fm^2 -> mb
            bglndif=sqrt(bglndif)
            if(bk(k).le.bglndif)nglndik=nglndik+1
c          print *,dble(xnpomcol),om1MeanCol(1,k),om1MeanCol(2,k),k
          endif
        endif

        !---------------------------------------------------------------
        !    formerly in emsigr
        !----------------------------------------------------------------

        if(irea.eq.0)then
        !determine length of k-th line of grid

        o=max(1.e-5,min(sngl(om1intc(k)),float(npommx)))!if GFF used for propo
         if(ish.ge.7)write(ifch,*)'emsigr:k,o',k,o
        n=0
        if(o.le.50)then
          p=1./(exp(o)-1)
        else
          p=0.
        endif
10      n=n+1
        p=p*o/n
        if(ish.ge.7)write(ifch,*)'emsigr:n,p',n,p
        if((p.gt.1e-5.or.n.lt.int(o)).and.n.lt.npommx
     * .and.n.lt.nprmax)goto 10
        !n=min(npommx,n*log((r2had(iclpro)+r2had(icltar))/r2part))  !if local b
        n=min(npommx,n*2)  !to allow backup of hard pom
        if(ish.ge.5)write(ifch,*)'k-th line: nmax,b',n,bk(k)
        npr(0,k)=n
        nprmx(k)=n
        nprt(k)=0
        do i=1,4
          npr(i,k)=0
        enddo
        endif

        !---------------------------------------------------------------
        !    former GfunParP (used to define parameters which depend on
        ! position n in a pair (for b fluctuations for instance))
        ! Now just a copy of the pair value.
        !----------------------------------------------------------------

        do n=1,nprmx(k)
        do i=ntymi,ntymx
          atildp(i,n,k)=atilde(i,k)
          btildpp(i,n,k)=btildep(i,k)
          btildppp(i,n,k)=btildepp(i,k)
        enddo
        enddo
       
       
c      endif
      enddo         !koll
c     enddo
c      if(koll.gt.1)then
cc        om1mean=om1mean/dble(koll)
c        do k=1,koll
cc          om1MeanCol(1,k)=om1min !fix to average between the lowest value and the average
cc         om1MeanCol(1,k)=0.5d0*(om1mean+om1min) !fix to average between the lowest value and the average
c          ip=iproj(k)
c          it=itarg(k)
c          om1mean=0.d0
c          do kxy=1,koll         !ckw21****consider the whole event, makes sens for pA, not for AA***
c            if(iproj(kxy).eq.ip.or.itarg(kxy).eq.it)then
c              om1mean=om1mean+om1MeanCol(1,kxy)
c            endif
c          enddo 
c          om1MeanCol(2,k)=om1mean  !average number of Pom for all connection of proj and targ og this pair.
cc        print *,k,om1mean
c        enddo                   !koll
c      endif
        
 999  if(test.or.ish.ge.10)then
        ifout=ifch
        if(test)then
          ifout=ifmt
          kmin=k
          kmax=k
        else
          kmin=1
          kmax=koll
        endif
        do k=kmin,kmax
          ip=iproj(k)
          it=itarg(k)
          write(ifout,*)' k,b,ip,it,gamb,alppar',k,bk(k),ip,it
     *                                           ,alppar
          write(ifout,*)'*******************************************'
          write(ifout,*)" atildg,btilgp,btildgpp "
          do ii=ntymin,ntymax
            write(ifout,*)ii,atildg(ii,k),btildgp(ii,k),btildgpp(ii,k)
          enddo
          write(ifout,*)'*******************************************'
          write(ifout,*)" atilde,btildep,btildepp "
          do ii=ntymin,ntymax
            write(ifout,*)ii,atilde(ii,k),btildep(ii,k),btildepp(ii,k)
          enddo
        enddo
        call printGfunParK
        if(test)call utstop('Error in epos-omg in GfunPark&')
      endif

      call utprjx('GfunParK ',ish,ishini,10)

      return
      end

c----------------------------------------------------------------------
      subroutine GfunParP
c----------------------------------------------------------------------
c for the moment metropolis should not depend on parton pair position
c (if fluctuations too large, doesn't work properly) ctp131022
c----------------------------------------------------------------------
#include "aaa.h"
#include "sem.h"
#include "ems.h"
#include "par.h"

      return

c      do k=1,koll
c        do n=1,nprmx(k)
c
c        if(abs(facpos).gt.100.)then
c          bhpr(n,k)=bk(k) 
c        endif
c
c        b2=bhpr(n,k)*bhpr(n,k)
c        do i=ntymin,ntymax

c         epsG =epsilongs(i,k)*exp(-(b2-b2a)/b2xscr)
c         epsGp=epsilongp(i,k)*exp(-(b2-b2a)/b2xscr)
c         epsGt=epsilongt(i,k)*exp(-(b2-b2a)/b2xscr)
c         gamb=gamD(i,iclpro,icltar)*b2
c     *       *(r2had(iclpro)+r2had(icltar))/(2.*r2part)
c
c         atildp(i,n,k)=dble(alpD(i,iclpro,icltar))
c         if(i.lt.0)then                     !use bk for metropolis
c           atildp(i,n,k)=atildp(i,n,k)*dble(
c     *             exp(-b2/delD(i,iclpro,icltar)
c     *                    *(r2had(iclpro)+r2had(icltar))/(2.*r2part))
c     *             *(1.+(r2had(iclpro)+r2had(icltar))/(slopom*log(sy))
c     *             )*sy**(betD(i,iclpro,icltar)
c     *                   +gamb+epsG))
c           btildpp(i,n,k)=dble(betDp(i,iclpro,icltar)
c     *                   +epsGp
c     *                   +gamb-alppar)
c           btildppp(i,n,k)=dble(betDpp(i,iclpro,icltar)
c     *                   +epsGt
c     *                   +gamb-alppar)
c         else           !keep average b for diffractive Pomeron
c           atildp(i,n,k)=atilde(i,k)
c           btildpp(i,n,k)=btildep(i,k)
c           btildppp(i,n,k)=btildepp(i,k)
c         endif
c
c         if(btildpp(i,n,k)+1d0.lt.-eps
c     *     .or.btildppp(i,n,k)+1d0.lt.-eps)then
c           write(ifmt,*)' k,n,b,ip,it,gamb,alppar',k,n,bhpr(n,k),ip,it
c     *                                            ,gamb,alppar
c           write(ifmt,*)' gammaV,epsGP0-2,epsGT0-2,epsGS0-2'
c     *          ,gammaV(k),epsilongp(0,k),epsilongp(1,k),epsilongp(2,k)
c     *                    ,epsilongt(0,k),epsilongt(1,k),epsilongt(2,k)
c     *                    ,epsilongs(0,k),epsilongs(1,k),epsilongs(2,k)
c           write(ifmt,*)'*******************************************'
c           write(ifmt,*)" atildp,btildpp,btildppp "
c           do ii=ntymin,ntymax
c             write(ifmt,*)ii,atildp(ii,n,k),btildpp(ii,n,k)
c     *                      ,btildppp(ii,n,k)
c           enddo
c           call utstop('Error in epos-omg in GfunPark&')
c         endif
       
c        enddo
c        enddo
c      enddo
c      return
      end

c----------------------------------------------------------------------
      subroutine GfunPar(zzip,zzit,eflp,eflt,mm,ii,b,spp,alp,bet,betp
     &                  ,epsp,epst,epss,gamvv)
c----------------------------------------------------------------------
c  calculates parameters alp,bet,betp of the G functions for pp (b-mode)
c  and the screening exponents epsp,epst,epss,gamvv
c----------------------------------------------------------------------
c  zzip:additional z component (nuclear effect projectile side)
c  zzit:additional z component (nuclear effect target side)
c  m=1: profile function Gpro,  ii=0: soft+diff  ii=2: Y+   ii=10 (hard only)
c  m=2: Gtilde,                 ii=1: semi       ii=3: Y-         (for m=1)
c  b: impact param, spp: pp energy squared (should correspond to engy**2)
c  m<0 to avoid asymetry nuclear effect    
c----------------------------------------------------------------------

#include "aaa.h"
#include "sem.h"
#include "par.h"
#include "ems.h"
      parameter(nbkbin=100)
      common /kfitd/ xkappafit(nbkbin,nclegy,nclha,nclha),xkappa,bkbin
      common /kwrite/ xkapZ
      double precision utgam2,coefgdp,coefgdt,dalp,dbet,dbetp,dbetsf
     .,dbetpsf,eps,xmn,xmn0,salp
      parameter(eps=1.d-20,betxx=-0.9)

      call utprj('GfunPar ',ish,ishini,10)

      if(isetcs.lt.-1)stop 'isetcs too low for Gfunpar, wrong call !'

      m=abs(mm)
      ee=sqrt(spp)
c      bglaub2=FbGlaub2(ee)
      rs=r2had(iclpro)+r2had(icltar)+slopom*log(spp)
      bglaub2=4.*.0389*rs
      if(ish.ge.10)write(ifch,*)'Gf in:',m,ii,b,bglaub2,spp
     &                                    ,iclpro,icltar
      i=abs(ii)
      if(ii.eq.10)i=1
      b2=b*b
      cfalpro=ucfpro
      cfaltar=ucftar
      gamb=gamD(i,iclpro,icltar)*b2
      xmn=1.d0/dble(spp)      !minimum mass 
      xmn0=xmn
      if(iregge.ne.0.and.s0min.lt.0.5d0)then
        xmn0=dble(s0min)*xmn0 !minimum limit=s0min GeV**2 for Reggeon
      elseif(iomega.lt.2)then
        xmn0=dble(s0min)*xmn0 !minimum limit=s0min GeV**2 for Reggeon
      endif

      zzp=0.
      zzt=0.
      zzpUni=0d0
      zztUni=0d0
      epsGp=0.
      epsGt=0.
      epsG=0.
      eGpsf=0.
      eGtsf=0.
      eGsf=0.
      gbsf=0.
      rhosf=0.
      gamV=1.
      
      if(iscreen.ne.0)then
       if(zzip.lt.0.and.zzit.lt.0.and.i.eq.1)then
c special definition to compute parameters without calling getZparEpsi
        betx=max(0.,betDp(1,iclpro,icltar)-alppar-betxx+min(0.,gamb))
        epsGp=abs(epscrx*zzip*betDp(1,iclpro,icltar))
        if(epscrh.lt.0)then
          epsG=epsG-epsGp
          epsGp=-min(betx,epsGp)
        endif
        betx=max(0.,betDpp(1,iclpro,icltar)-alppar-betxx+min(0.,gamb))
        epsGt=abs(epscrx*zzit*betDpp(1,iclpro,icltar))
        if(epscrh.lt.0)then
          epsG=epsG-epsGt
          epsGt=-min(betx,epsGt)
        endif
        epsG=0.5*epsG          
       else
c full definition of screening parameters (depends on some definitions which are energy dependent)
        call getZparEpsi(mm,ii,b,zzip,zzit,eflp,eflt,epsG,epsGp,epsGt)
        if(i.eq.0.and.alpDp(0,iclpro,icltar).gt.0)then
          call getZparEpsi(mm,1,b,zzip,zzit,eflp,eflt,eGsf,eGpsf,eGtsf)
          gbsf=alpDp(2,iclpro,icltar)*b2    
          rhosf=eGsf+gbsf
        elseif(ii.eq.10)then
          gbsf=alpDpp(2,iclpro,icltar)*b2    
          rhosf=gbsf
        endif
       endif
      endif

      rho=epsG+gamb!+betD(i,iclpro,icltar)+gamb   !ctp20160906 s-independant fit now 

      if(m.eq.1)then

          if(i.eq.0.and.alpDp(0,iclpro,icltar).gt.0)then
c for soft+diff component, exponent should not be higher than soft+screening one
            dbetsf=dble(max(betxx,alpDp(1,iclpro,icltar)
     *         +eGpsf
     *         +gbsf-alppar))
            dbetpsf=dble(max(betxx,alpDp(1,iclpro,icltar)
     *         +eGtsf
     *         +gbsf-alppar))
          elseif(ii.eq.10)then
            dbetsf=dble(max(betxx,alpDpp(1,iclpro,icltar)
     *         +gbsf-alppar))
            dbetpsf=dble(max(betxx,alpDpp(1,iclpro,icltar)
     *         +gbsf-alppar))
          else
            dbetsf=1d35
            dbetpsf=1d35
          endif
c        if(i.le.idxD)then
          dbet=min(dbetsf,dble(max(betxx,betDp(i,iclpro,icltar)
     *         +epsGp
     *         +gamb-alppar)))
          dbetp=min(dbetpsf,dble(max(betxx,betDpp(i,iclpro,icltar)
     *         +epsGt
     *         +gamb-alppar)))
          if(i.eq.0.and.alpDp(0,iclpro,icltar).gt.0)then
c save proper eps
            if(abs(dbet-dbetsf).lt.1d-10)epsGp=eGpsf
            if(abs(dbetp-dbetpsf).lt.1d-10)epsGt=eGtsf
c for soft+diff component, factor should not be lower than soft+screening one
            salp=dble(alpDp(0,iclpro,icltar))
     *       *exp(min(50d0
     *          ,dble(-rhosf*log(xmn0)-b2/alpDp(3,iclpro,icltar))))
c recalculate alpha for proper exponent
            salp=salp*xmn**(0.5d0*(dbetsf+dbetpsf-dbet-dbetp))
          elseif(ii.eq.10)then
            salp=dble(alpDpp(0,iclpro,icltar))
     *       *exp(min(50d0
     *          ,dble(-rhosf*log(xmn0)-b2/alpDpp(3,iclpro,icltar))))
          else
            salp=0d0
          endif
          dalp=dble(alpD(i,iclpro,icltar))
     *     *exp(min(50d0,
     *          dble(-rho*log(xmn0)-b2/delD(i,iclpro,icltar))))
          if(ii.eq.10)then
            if(dalp.ge.salp)then
c if hard component below fit, then use fit of hard component without screening
              dalp=salp
              dbet=dbetsf
              dbetp=dbetpsf
              epsG=0.
              epsGp=0.
              epsGt=0.
            endif
          elseif(dalp.le.salp)then
            dalp=salp
            epsG=eGsf
          endif
c       else
c          absb=abs(b)-bmxdif(iclpro,icltar)
c          b2a=absb*absb
c          dalp=dalp
c     *       *exp(min(50d0,dble((betD(i,iclpro,icltar)+epsG)*log(spp)
c     *            -b2a/delD(i,iclpro,icltar))))
c          dbet=dble(betDp(i,iclpro,icltar)-alppar+epsGp)    !WWWdiffr
c          dbetp=dble(betDpp(i,iclpro,icltar)-alppar+epsGt)  !WWWdiffr
c        endif

        if((dbet+1.d0).lt.-eps.or.(dbetp+1.d0).lt.-eps)then
        write(ifch,*)'m,i,b,spp,alp,bet,betp',m,i,b,spp,dalp,dbet,dbetp
        write(ifch,*)'betDp(i),betDpp(i),epsGp,epsGt,gamb'
     * ,betDp(i,iclpro,icltar),betDpp(i,iclpro,icltar),epsGp,epsGt,gamb
          call utstop('Error : beta < -1 in Gfunpar in epos-omg&')
        endif



      elseif(m.eq.2)then
      
        xkappa=1.
c        if(i.eq.0.and.b.lt.(nbkbin-1)*bkbin)then
        if(b.lt.(nbkbin-1)*bkbin.and.isetcs.ge.1)then
          ibk=int(b/bkbin)+1
          if(isetcs.gt.1.and.iclegy.lt.iclegy2)then
            egy0=egylow*egyfac**float(iclegy-1)
            xkappa1=xkappafit(ibk,iclegy,iclpro,icltar)
     *         +(b-bkbin*float(ibk-1))/bkbin
     *         *(xkappafit(ibk+1,iclegy,iclpro,icltar)
     *         -xkappafit(ibk,iclegy,iclpro,icltar))
            xkappa2=xkappafit(ibk,iclegy+1,iclpro,icltar)
     *         +(b-bkbin*float(ibk-1))/bkbin
     *         *(xkappafit(ibk+1,iclegy+1,iclpro,icltar)
     *         -xkappafit(ibk,iclegy+1,iclpro,icltar))
            xkappa=xkappa1+(xkappa2-xkappa1)/log(egyfac)
     *         *(log(ee)-log(egy0))
c            xkappa=(1.+facmc*exp(-b2/bglaub2))*xkappa
          else
            xkappa=xkappafit(ibk,iclegy,iclpro,icltar)
     *         +(b-bkbin*float(ibk-1))/bkbin
     *         *(xkappafit(ibk+1,iclegy,iclpro,icltar)
     *         -xkappafit(ibk,iclegy,iclpro,icltar))
c            xkappa=(1.+facmc*exp(-b2/bglaub2))*xkappa
          endif
        endif
        xkapZ=xkappa

        falp=cfalpro*cfaltar *gamV*xkappa*xkappa

        if(i.eq.0.and.alpDp(0,iclpro,icltar).gt.0)then
c for soft+diff component, exponent should not be higher than soft+screening one
          dbetsf=dble(max(betxx,alpDp(1,iclpro,icltar)
     *         +eGpsf
     *         +gbsf-alppar))
          dbetpsf=dble(max(betxx,alpDp(1,iclpro,icltar)
     *         +eGtsf
     *         +gbsf-alppar))
        else
          dbetsf=1d35
          dbetpsf=1d35
        endif
c        if(i.le.idxD)then
          dbet=min(dbetsf,dble(max(betxx,betDp(i,iclpro,icltar)
     *        +epsGp
     *        +gamb-alppar)))
          dbetp=min(dbetpsf,dble(max(betxx,betDpp(i,iclpro,icltar)
     *        +epsGt
     *        +gamb-alppar)))

          if(i.eq.0.and.alpDp(0,iclpro,icltar).gt.0)then

c save proper eps
            if(abs(dbet-dbetsf).lt.1d-10)epsGp=eGpsf
            if(abs(dbetp-dbetpsf).lt.1d-10)epsGt=eGtsf
c for soft+diff component, factor should not be lower than soft+screening one
            salp=dble(alpDp(0,iclpro,icltar))
     *       *exp(min(50d0
     *          ,dble(-rhosf*log(xmn0)-b2/alpDp(3,iclpro,icltar))))
c recalculate alpha for proper exponent
            salp=salp*xmn**(0.5d0*(dbetsf+dbetpsf-dbet-dbetp))
          else
            salp=0d0
          endif
          dalp=dble(alpD(i,iclpro,icltar))
     *     *exp(min(50d0,
     *          dble(-rho*log(xmn0)-b2/delD(i,iclpro,icltar))))
          if(dalp.gt.salp)then
            dalp=dble(falp)*dalp
          else
            dalp=dble(falp)*salp
            epsG=eGsf
          endif
          dbet=dbet+1d0
          dbetp=dbetp+1d0
c        else
c          absb=abs(b)-bmxdif(iclpro,icltar)
c          b2a=absb*absb
c          dalp=dalp
c     *     *exp(min(50d0,dble((betD(i,iclpro,icltar)+epsG)*log(spp)
c     *          -b2a/delD(i,iclpro,icltar))))
c          dbet=dble(betDp(i,iclpro,icltar)-alppar+1.+epsGp)     !WWWdiffr
c          dbetp=dble(betDpp(i,iclpro,icltar)-alppar+1.+epsGt)   !WWWdiffr
c        endif
        coefgdp=utgam2(1.d0+dble(alplea(iclpro))+dbet)
        coefgdt=utgam2(1.d0+dble(alplea(icltar))+dbetp)
        dalp=dalp*utgam2(dbet)*utgam2(dbetp)/coefgdp/coefgdt
      
      else

        stop'GproPar: wrong m value.              '

      endif


      alp=sngl(dalp)
      bet=sngl(dbet)
      betp=sngl(dbetp)
      if(abs(bet+1.).lt.eps.or.abs(betp+1.).lt.eps)
     .call utstop("bet=-1&")
      epsp=epsGp
      epst=epsGt
      epss=epsG
      gamvv=gamV
c update value of Z used (should not for SigGFFk)
c      zzip=zzp
c      zzit=zzt

      if(ii.ne.10)then
        alpUni(i,max(m,1))=dalp
        betUni(i,max(m,1))=dbet
        betpUni(i,max(m,1))=dbetp
        epssUni(i)=dble(epss)
        epspUni(i)=dble(epsp)
        epstUni(i)=dble(epst)
      endif
c      zzpUni=dble(zzp)
c      zztUni=dble(zzt)


      if(ish.ge.10)write(ifch,*)'   GfunPar :',alp,bet,betp,epsp,epst
     &                                    ,epss,gamvv,zzip,zzit

      call utprjx('GfunPar ',ish,ishini,10)
      end

c----------------------------------------------------------------------
      function FbGlaub2(x)
c----------------------------------------------------------------------
c  calculates (glauber radius)^2 from pp cross section (data fit)
c(estimated if not already calculated --> not any more to have smoother xs)
c----------------------------------------------------------------------
c  x: pp energy
c----------------------------------------------------------------------

#include "aaa.h"

c      if(sigine.eq.0.)then
        if(iclpro+icltar.eq.3)then !pi+p
          siginex=20.+0.08*log(x)**3.-0.004*log(x)**4.
        elseif(iclpro+icltar.eq.5)then !K+p
          siginex=16.+0.08*log(x)**3.-0.004*log(x)**4.
        else
          siginex=30.+0.095*log(x)**3.-0.004*log(x)**4.
        endif
c      else
c       siginex=sigine
c      endif
      FbGlaub2=siginex/10./pi

      return
      end


c----------------------------------------------------------------------
      double precision function om1(xh,yp,b)   !---test---
c----------------------------------------------------------------------
c om1 = G * C * gamd    (C and gamd usually 1)
c xh - fraction of the energy squared s for the pomeron;
c b - impact parameter between the pomeron ends;
c yp - rapidity for the pomeron;
c----------------------------------------------------------------------
#include "aaa.h"
#include "sem.h"
#include "par.h"

      double precision Gf,xp,xm,xh,yp

      Gf=0.d0
      xp=sqrt(xh)*exp(yp)
      xm=xh/xp
      spp=engy**2
      do i=idxDmin(iomega),idxDmax(iomega)
        zp=0.
        zt=0.
        call Gfunpar(zp,zt,0.,0.,1,i,b,spp,alp,bet,betp,epsp,epst,epss
     .              ,gamv)
        Gf=Gf+dble(alp)*xp**dble(bet)*xm**dble(betp)
      enddo
      om1=Gf
      end

c----------------------------------------------------------------------
      double precision function om1intb(b)   !---test---
c----------------------------------------------------------------------
c  om1 integrated over xp and xm for given b
c  Calculation by analytical integration
c----------------------------------------------------------------------
#include "aaa.h"
#include "sem.h"
#include "par.h"
      double precision cint,cint2,eps
      parameter(eps=1.d-20)

      spp=engy*engy
      cint=0.d0
      do i=idxDmin(iomega),idxDmax(iomega)
        zp=0.
        zt=0.
      call Gfunpar(zp,zt,0.,0.,1,i,b,spp,alp,bet,betp,epsp,epst,epss
     .              ,gamv)
        cint2=dble(gamv*alp)
        if((bet+1.0).gt.eps)then
          cint2=cint2/dble(bet+1.0)
        else
          cint2=-cint2*log(xminDf)
        endif
        if((betp+1.0).gt.eps)then
          cint2=cint2/dble(betp+1.0)
        else
          cint2=-cint2*log(xminDf)
        endif
        cint=cint+cint2
      enddo

      if(cint.lt.-eps)then
        write(*,*) 'WARNING ! om1intb in epos-omg is <0 !!!!!'
        write(*,*) 'WARNING ! => om1intb set to 1e-3 !!!!!'
        write(*,*) 'WARNING ! => bmax=3.5 fm !!!!!'
        cint=1.d-3
      endif

      om1intb=cint

      return
      end

c----------------------------------------------------------------------
      double precision function om1intbi(b,iqq)   !---test---
c----------------------------------------------------------------------
c  om1 integrated over xp and xm for given b
c  Calculation by analytical integration of contribution iqq
c  iqq = 0 - soft
c  iqq = 1 - sh
c  iqq = 2 - diff
c----------------------------------------------------------------------
#include "aaa.h"
#include "sem.h"
#include "par.h"
      double precision eps,cint
      parameter(eps=1.d-20)

      om1intbi=0.D0
      zp=0.
      zt=0.
      spp=engy*engy
      if(iqq.eq.1)then
        imin=1
        call Gfunpar(zp,zt,0.,0.,1,1,b,spp,alp,bet,betp,epsp,epst,epss
     .              ,gamv)
        if(alp.le.0.)imin=idxD
        imax=idxD
      elseif(iqq.eq.2)then
        if(iomega.ge.2)return
        imin=idxD1-1
        imax=idxD1
      else
        imin=idxDmin(iomega)
        imax=idxDmin(iomega)
        call Gfunpar(zp,zt,0.,0.,1,0,b,spp,alp,bet,betp,epsp,epst,epss
     .              ,gamv)
        if(alp.lt.0)return
      endif
      do i=imin,imax
      call Gfunpar(zp,zt,0.,0.,1,i,b,spp,alp,bet,betp,epsp,epst,epss
     .              ,gamv)
        cint=dble(gamv*alp)
        if(dble(bet+1.0).gt.eps)then
          cint=cint/dble(bet+1.0)
        else
          cint=-cint*log(xminDf)
        endif
        if(dble(betp+1.0).gt.eps)then
          cint=cint/dble(betp+1.0)
        else
          cint=-cint*log(xminDf)
        endif
        om1intbi=om1intbi+cint
      enddo
      if(om1intbi.lt.-eps)then
        write(*,*) 'WARNING ! om1intbi in epos-omg is <0 !!!!!'
        write(*,*) 'WARNING ! => om1intbi set to 1e-6 !!!!!'
        write(*,*) 'WARNING ! => bmax=3.5 fm !!!!!'
        om1intbi=1.d-6
      endif


      return
      end


c----------------------------------------------------------------------
      double precision function om1intby(spp,b,iqq,imod)   !---test---
c----------------------------------------------------------------------
c  om1=int dxp dxm Gdiff(x)*Phi 
c  integrated over xp and xm for given b
c  Calculation by numerical integration of diff contribution for iqq contri.
c                       iqq=5 .... Reggeon
c                       iqq=-5 ... diff(soft (6 to 9)) for fit
c                       iqq=-50 .. diff(soft (6 to 9)) for MC
c                       iqq=6 .... SD- (pro non excited, mass to targ) real
c                       iqq=-6 ... SD- (pro non excited, mass to targ) fit
c                       iqq=-60 .. SD- (pro non excited, mass to targ) MC
c                       iqq=7 .... SD+ (targ non excited, mass to pro) real
c                       iqq=-7 ... SD+ (targ non excited, mass to pro) fit
c                       iqq=-70 .. SD+ (targ non excited, mass to pro) MC
c                       iqq=8 .... CD  (soft : pro and targ non excited
c                                       mass to resonance or string)     real
c                       iqq=-8 ... CD  (soft : pro and targ non excited) fit
c                       iqq=-80 .. CD  (soft : pro and targ non excited) ~MC
c                       iqq=9 .... DD  (mass to pro an targ) real
c                       iqq=-9 ... DD  (mass to pro an targ) fit
c                       iqq=-90 .. DD  (mass to pro an targ) MC
c  imod - calculation without (0-diagram) or with Phi (1-cross section)
c     (-1 is for test with analytical integration)
c  !!! Gfunpar should be called before for proper z and b (if imod.ne.-1) !!!!
c----------------------------------------------------------------------
#include "aaa.h"
#include "sem.h"
#include "par.h"
      double precision cint,w,xmmax,xmmin,xpmin,xpmax,xp,xm,om5J,ww
     .,PhiUnit,xvpr,om1intbi
      common /ar3/  x1(7),a1(7)
      common/cems13/xvpr(0:3)
      double precision dgx1,dga1
      common /dga20/ dgx1(10),dga1(10)

      om1intby=0d0

      if(imod.eq.-1)then
        if(iqq.eq.-50.or.iqq.eq.-80) om1intby=om1intby+om1intbi(b,0)
     .            *0.5d0*dble(gampar**2)/(1d0+0.5d0*dble(gampar**0))
        if(iqq.eq.-50.or.iqq.eq.-90) om1intby=om1intby+om1intbi(b,2)
     .                                  /(1d0+0.5d0*dble(gampar**2))
        if(iqq.eq.-50.or.iqq.eq.-60) om1intby=om1intby+om1intbi(b,2)
        if(iqq.eq.-50.or.iqq.eq.-70) om1intby=om1intby+om1intbi(b,3)
        return
      elseif(imod.eq.2)then
        do i=idxDmin(iomega),idxDmax(iomega)
          call Gfunpar(0.,0.,0.,0.,1,i,b,spp,xa,xb,xc,xd,xe,xf,xg)
          call Gfunpar(0.,0.,0.,0.,2,i,b,spp,xa,xb,xc,xd,xe,xf,xg)
        enddo
      endif
      
c xp and xm integration with Pom
      iqqp=abs(iqq)
      if(iqqp.gt.10)iqqp=iqqp/10
      cint=0d0
      xpmin=s0min/dble(spp)
      xpmax=1d0
      if(iqqp.eq.8.or.iqqp.eq.6)xpmax=dble(xmxrem(1))
      if(iqq.gt.0)xpmin=max(xpmin,xpmax*1d-7)          !ctp20180419 ??????????? whitout this, the diff cross-section increase too fast at very high energy ????????
      xmmin=s0min/dble(spp)
      xmmax=1d0
      if(iqqp.eq.8.or.iqqp.eq.7)xmmax=dble(xmxrem(2))
      if(iqq.gt.0)xmmin=max(xmmin,xmmax*1d-7)          !ctp20180419 ??????????? whitout this, the diff cross-section increase too fast at very high energy ????????
c      if(iqq.gt.0)xmmax=min(xmmax,xmmin*1d10)
c     cross section for 1 diffractive Pomeron exchange (multiple scattering negligible)
      if(xmmax.gt.xmmin.and.xpmax.gt.xpmin)then

        do n=1,2
          do i=1,10 !
c            xp=xpmin*(xpmax/xpmin)**dble(.5+x1(i)*(n-1.5))
            xp=xpmin*(xpmax/xpmin)**(.5d0+dgx1(i)*(n-1.5))
            w=0d0
            do m=1,2
              do j=1,10 !7
c                xm=xmmin*(xmmax/xmmin)**dble(.5+x1(j)*(m-1.5))
                xm=xmmin*(xmmax/xmmin)**(.5d0+dgx1(j)*(m-1.5))
c                w=w+dble(a1(j))*om5J(xp,xm,b,iqq,iqq)*xp*xm
                ww=om5J(xp,xm,b,iqq,iqq)
                if(imod.ge.1)ww=ww*PhiUnit(1d0-xp,1d0-xm)
                w=w+dga1(j)*ww*xp*xm
              enddo
            enddo
c            cint=cint+dble(a1(i))*0.5d0*log(xmmax/xmmin)*w
            cint=cint+dga1(i)*0.5d0*log(xmmax/xmmin)*w
          enddo
        enddo
        om1intby=cint*0.5d0*log(xpmax/xpmin)
      endif

c      print *,'om1intby',iqq,b,xmxrem,om1intby
      return
      end

c----------------------------------------------------------------------
      double precision function om1intbyk(k)   !---test---
c----------------------------------------------------------------------
c  om1=int dxp dxm Gdiff(x)*Phi 
c  integrated over xp and xm for given pair k (including nuclear effects)
c  Calculation by numerical integration for total diff contribution
c  !!! Gfunpark should be called before for proper z and b !!!!
c----------------------------------------------------------------------
#include "aaa.h"
#include "ems.h"
#include "sem.h"
#include "par.h"
      double precision cint,w,xmmax,xmmin,xpmin,xpmax,xp,xm,omGamk,ww
     .,PhiExpoK,xpsup,xmsup
      common /ar3/  x1(7),a1(7)
      double precision dgx1,dga1
      common /dga20/ dgx1(10),dga1(10)

      om1intbyk=0d0

      if(iomega.ge.2)return
     
c xp and xm integration with Pom
      cint=0d0
      xpmin=XminDf
      xpmax=1d0
      xmmin=XminDf
      xmmax=1d0
c cross section for 1 diffractive Pomeron exchange (multiple scattering negligible)
      if(xmmax.gt.xmmin.and.xpmax.gt.xpmin)then

        do n=1,2
          do i=1,10 !
c            xp=xpmin*(xpmax/xpmin)**dble(.5+x1(i)*(n-1.5))
            xp=xpmin*(xpmax/xpmin)**(.5d0+dgx1(i)*(n-1.5))
            xpsup=max(0d0,1d0-xp/dble(xmxrem(1))) !large xp suppressed
            w=0d0
            do m=1,2
              do j=1,10 !7
c                xm=xmmin*(xmmax/xmmin)**dble(.5+x1(j)*(m-1.5))
                xm=xmmin*(xmmax/xmmin)**(.5d0+dgx1(j)*(m-1.5))
                xmsup=max(0d0,1d0-xm/dble(xmxrem(2))) !large xm suppressed
                ww=0d0
                ww=ww+omGamk(1,k,xp,xm,2,2)
     .                      *xpsup
                ww=ww+omGamk(1,k,xp,xm,3,3)
     .                      *xmsup
                ww=ww+omGamk(1,k,xp,xm,0,0)
     .                      *0.5d0*dble(gampar**2)*xpsup*xmsup
     .                      /(1d0+0.5d0*dble(gampar**2)*xpsup*xmsup)
                ww=ww+omGamk(1,k,xp,xm,0,0)
     .                      /(1d0+0.5d0*dble(gampar**2)*xpsup*xmsup)
                ww=ww*PhiExpok(k,1d0-xp,1d0-xm)
     &               *(1d0-xp)**dble(alplea(iclpro))
     &               *(1d0-xm)**dble(alplea(icltar))
                w=w+dga1(j)*ww*xp*xm
              enddo
            enddo
c            cint=cint+dble(a1(i))*0.5d0*log(xmmax/xmmin)*w
            cint=cint+dga1(i)*0.5d0*log(xmmax/xmmin)*w
          enddo
        enddo
        om1intbyk=cint*0.5d0*log(xpmax/xpmin)
      endif
      return
      end

c----------------------------------------------------------------------
      double precision function om1intbhj(b,kk)   !---MC---
c----------------------------------------------------------------------
c  om1=int dxp dxm Ghard(x)*F*F 
c  integrated over xp and xm for given b (and pair kk if !=0)
c  (if kk<0 return om1=int dxp dxm xp*xm*Ghard(x)*F*F)
c----------------------------------------------------------------------
#include "aaa.h"
#include "ems.h"
#include "sem.h"
#include "par.h"
      double precision cint,w,xmmax,xmmin,xpmin,xpmax,xp,xm,om5J,xvpr,ww
c     .,om1intbh
      common /ar3/  x1(7),a1(7)
      common/cems13/xvpr(0:3)

      om1intbhj=0d0
      spp=engy*engy
      if(kk.eq.0)then
c     define variables for for om5J and PhiUnit
        zp=0.
        zt=0.
        i=-1
c      do i=idxDmin(iomega),idxDmax(iomega)
       call Gfunpar(zp,zt,0.,0.,1,i,b,spp,alpx,betx,betpx,epsp,epst,epss
     .,gamv)
        i=1
       call Gfunpar(zp,zt,0.,0.,1,i,b,spp,alpx,betx,betpx,epsp,epst,epss
     .,gamv)
c      enddo
      elseif(kk.lt.0)then
        k=abs(kk)
        alpUni(1,1)=atilde(1,k)
        betUni(1,1)=btildep(1,k)
        betpUni(1,1)=btildepp(1,k)
        epssUni(-1)=epsilongs(-1,k)
        epspUni(-1)=epsilongp(-1,k)
        epstUni(-1)=epsilongt(-1,k)
      endif

c xp and xm integration with Pom
      cint=0d0
      xpmin=dble(4.*q2nmin)/dble(spp)
      xpmax=1d0
      xmmin=dble(4.*q2nmin)/dble(spp)
      xmmax=1d0
c cross section for 1 diffractive Pomeron exchange (multiple scattering negligible)
      if(xmmax.gt.xmmin.and.xpmax.gt.xpmin)then

        do n=1,2
          do i=1,7
            xp=xpmin*(xpmax/xpmin)**dble(.5+x1(i)*(n-1.5))
            w=0d0
            do m=1,2
              do j=1,7
                xm=xmmin*(xmmax/xmmin)**dble(.5+x1(j)*(m-1.5))
                if(xp*xm*dble(spp).gt.1.d0)then
c                  ww=min(om5J(xp,xm,b,1,1),om5J(xp,xm,b,-3,-3)) !same definition as in Womty
                  ww=om5J(xp,xm,b,101,104) !+om5J(xp,xm,b,11,11) !soft+gg=0 at q2nmin and should not be counted to get hard scaling (do not contribute to high pt)
c                  ww=0.5d0*(om5J(xp,xm,b,-10,-10)
c     .                      *ww/(om5J(xp,xm,b,0,0)+ww)
c     .                     +min(ww,om5J(xp,xm,b,-3,-3)))
                  ww=min(ww,om5J(xp,xm,b,-3,-3))
                  if(kk.lt.0)ww=ww*xp*xm
                  w=w+dble(a1(j))*xp*xm*max(0d0,ww)
     .               *(1d0-xp)**alplea(iclpro)
     .               *(1d0-xm)**alplea(icltar)
                endif
              enddo
            enddo
            cint=cint+dble(a1(i))*0.5d0*log(xmmax/xmmin)*w
          enddo
        enddo
        om1intbhj=cint*0.5d0*log(xpmax/xpmin)
      endif


c      print *,'om1intbhj',iqq,b,om1intbhj!,om1intbh(b)

      return
      end


cc----------------------------------------------------------------------
c      double precision function om1intbh(b)   !---MC---
cc----------------------------------------------------------------------
cc  om1*F*F integrated over xp and xm for given b
cc  Calculation by analytical integration for hard Pom (fit)
cc  2018.03.29 Not used any more.
cc----------------------------------------------------------------------
c#include "aaa.h"
c#include "ems.h"
c#include "sem.h"
c      double precision cint,gamom,deltap,deltam
c      double precision utgam2,Fp,Fm
c
c      spp=engy**2
c      om1intbh=0.d0
c      Fp=dble(ucfpro)   !gamma(1+alplea)
c      Fm=dble(ucftar)
c
c      cint=0.d0
c      i=10
c      
c        zp=0.
c        zt=0.
c      call Gfunpar(zp,zt,0.,0.,1,i,b,spp,alp,bet,betp,epsp,epst,epss,gamv)
c        gamom=dble(alp*gamv)
c        deltap=dble(bet)
c        deltam=dble(betp)
c        cint=cint+gamom*utgam2(deltap+1.d0)*utgam2(deltam+1.d0)
c     &            /utgam2(2.d0+deltap+dble(alplea(iclpro)))
c     &            /utgam2(2.d0+deltam+dble(alplea(icltar)))
c
c      om1intbh=max(0d0,cint*Fp*Fm)
c
cc      if(om1intbh.lt.-1.d-10)then
cc        write(*,*) 'WARNING ! om1intbh in omg is <0 !!!!!'
cc        write(*,*) 'WARNING ! => om1intbh set to 0. !!!!!'
cc        om1intbh=0.d0
cc      endif
c
c      return
c      end
c

c----------------------------------------------------------------------
      double precision function om1intbc(b)   !---MC---
c----------------------------------------------------------------------
c  om1*F*F integrated over xp and xm for given b
c  Calculation by analytical integration
c----------------------------------------------------------------------
#include "aaa.h"
#include "ems.h"
#include "sem.h"
      double precision cint,gamom,deltap,deltam
      double precision utgam2,Fp,Fm

      spp=engy**2
      om1intbc=0.d0
      Fp=dble(ucfpro)   !gamma(1+alplea)
      Fm=dble(ucftar)

      cint=0.d0
      do i=idxDmin(iomega),idxDmax(iomega)
      
        zp=0.
        zt=0.
      call Gfunpar(zp,zt,0.,0.,1,i,b,spp,alp,bet,betp,epsp,epst,epss
     .              ,gamv)
        gamom=dble(alp*gamv)
        deltap=dble(bet)
        deltam=dble(betp)
        cint=cint+gamom*utgam2(deltap+1.d0)*utgam2(deltam+1.d0)
     &            /utgam2(2.d0+deltap+dble(alplea(iclpro)))
     &            /utgam2(2.d0+deltam+dble(alplea(icltar)))
      enddo

      om1intbc=max(0d0,cint*Fp*Fm)

c      if(om1intbc.lt.-1.d-10)then
c        write(*,*) 'WARNING ! om1intbc in omg is <0 !!!!!'
c        write(*,*) 'WARNING ! => om1intbc set to 0. !!!!!'
c        om1intbc=0.d0
c      endif

      return
      end

c----------------------------------------------------------------------
      double precision function om1intbci(iqq)   !---MC---
c----------------------------------------------------------------------
c     om1*F*F integrated over xp and xm with analytical integration
c     with parameters from Gfunpar called before
c     iqq=1 - non diff only
c----------------------------------------------------------------------
#include "aaa.h"
#include "ems.h"
#include "sem.h"
      double precision cint,gamom,deltap,deltam
      double precision utgam2,Fp,Fm

      om1intbci=0.d0
      Fp=dble(ucfpro)   !gamma(1+alplea)
      Fm=dble(ucftar)

      imin=ntymin
      imax=ntymax
      if(iqq.eq.1.or.iomega.eq.2)then
        imax=idxD
        imin=idxD
      endif

      cint=0.d0
      do i=imin,imax
      
        gamom=alpUni(i,1)
        deltap=betUni(i,1)
        deltam=betpUni(i,1)
        cint=cint+gamom*utgam2(deltap+1.d0)*utgam2(deltam+1.d0)
     &            /utgam2(2.d0+deltap+dble(alplea(iclpro)))
     &            /utgam2(2.d0+deltam+dble(alplea(icltar)))
      enddo

      om1intbci=max(0d0,cint*Fp*Fm)


      return
      end


c----------------------------------------------------------------------
      double precision function om1intgk(ki,xprem,xmrem)   !---MC---
c----------------------------------------------------------------------
c  om1*F(xprem-xp)*F(xmrem-xm) integrated over xp and xm for given k
c  Calculation by analytical integration
c  if ki<0, use om1 with saturation
c----------------------------------------------------------------------
#include "aaa.h"
#include "ems.h"
#include "sem.h"
#include "par.h"
      double precision cint,gamom,deltap,deltam,xprem,xmrem,utgam2

      k=abs(ki)
      corr=0.
c      if(ki.lt.0)then
c        bc=bcutz(fegypp,b2xscr)
c        b2b=max(min(0.,bc),abs(bk(k))-abs(bc))
c        b2a=max(0.,b2b)**2
c        b2b=b2b**2
c        corr=exp(-(b2a-b2b)/b2xscr)-1.
c      endif

      om1intgk=0.d0
      Fp=dble(ucfpro)   !gamma(1+alplea)
      Fm=dble(ucftar)

      cint=0.d0
      do i=ntymin,ntymax
        gamom=dble(atilde(i,k))
        deltap=dble(btildep(i,k))
        deltam=dble(btildepp(i,k))
        if(ki.lt.0.and.i.eq.1)then
c use om1 with saturated Z (not decreasing for b->0)
          epss=epsilongs(1,k)
          if(epss.lt.0.)gamom=gamom/xminDf**dble(corr*epss)
          deltap=deltap+dble(corr*epsilongp(1,k))
          deltam=deltam+dble(corr*epsilongt(1,k))
        endif
c Alternative om1*(xprem-xp)**alpea*(xmrem-xm)**alplea and integration from 0 to 1
        cint=cint+gamom*utgam2(deltap+1.d0)*utgam2(deltam+1.d0)
     &            /utgam2(2.d0+deltap+dble(alplea(iclpro)))*Fp
     &            /utgam2(2.d0+deltam+dble(alplea(icltar)))*fm
     &            *xprem**(1d0+deltap+dble(alplea(iclpro)))
     &            *xmrem**(1d0+deltam+dble(alplea(icltar)))
c Alternative om1*(xpr-xp)*(xmr-xm) and integration from 0 to xr
c        cint=cint+gamom/(deltap+1.d0)/(deltam+1.d0)
c     &            /(2.d0+deltap) /(2.d0+deltam)
c     &            *xprem**(deltap+2.d0)
c     &            *xmrem**(deltam+2.d0)
c Alternative om1 only and integration from 0 to xr
c        cint=cint+gamom/(deltap+1.d0)/(1.d0+deltam)
c     &            *xprem**(deltap+1.d0)
c     &            *xmrem**(deltam+1.d0)
      enddo
      om1intgk=cint

      return
      end

c----------------------------------------------------------------------
      double precision function om1intgck(n,k,xprem,xmrem)   !---MC---
c----------------------------------------------------------------------
c  om1*(xprem-xp)*(xmrem-xm) integrated over xp and xm for given 
c  position (n,k)
c  Calculation by analytical integration
c----------------------------------------------------------------------
#include "aaa.h"
#include "ems.h"
#include "sem.h"
      double precision cint,gamom,deltap,deltam,xprem,xmrem,utgam2

      om1intgck=0.d0
      Fp=dble(ucfpro)   !gamma(1+alplea)
      Fm=dble(ucftar)

      cint=0.d0
      do i=ntymin,ntymax
        gamom=atildp(i,n,k)
        deltap=btildpp(i,n,k)
        deltam=btildppp(i,n,k)
c Alternative om1*(xpr-xp)**alpea*(xmr-xm)**alplea and integration from 0 to 1
        cint=cint+gamom*utgam2(deltap+1.d0)*utgam2(deltam+1.d0)
     &            /utgam2(2.d0+deltap+dble(alplea(iclpro)))*Fp
     &            /utgam2(2.d0+deltam+dble(alplea(icltar)))*fm
     &            *xprem**(1d0+deltap+dble(alplea(iclpro)))
     &            *xmrem**(1d0+deltam+dble(alplea(icltar)))
c Alternative om1*(xpr-xp)*(xmr-xm) and integration from 0 to xr
c        cint=cint+gamom/(deltap+1.d0)/(deltam+1.d0)
c     &            /(2.d0+deltap) /(2.d0+deltam)
c     &            *xprem**(deltap+2.d0)
c     &            *xmrem**(deltam+2.d0)
c Alternative om1 only and integration from 0 to xr
c        cint=cint+gamom/(deltap+1.d0)/(1.d0+deltam)
c     &            *xprem**(deltap+1.d0)
c     &            *xmrem**(deltam+1.d0)
      enddo
      om1intgck=cint

      return
      end

c----------------------------------------------------------------------
      double precision function om1intgc(b)   !---xs---
c----------------------------------------------------------------------
c  om1*(1-xp)*(1-xm) integrated over xp and xm for given b
c  Calculation by analytical integration
c----------------------------------------------------------------------
#include "aaa.h"
#include "par.h"
#include "sem.h"
      double precision cint,gamom,deltap,deltam,eps
      parameter(eps=1.d-20)

      spp=engy**2
      om1intgc=0.d0


      cint=0.d0

      do i=idxDmin(iomega),idxDmax(iomega)
        zp=0.
        zt=0.
      call Gfunpar(zp,zt,0.,0.,1,i,b,spp,alp,bet,betp,epsp,epst,epss
     .              ,gamv)
        gamom=dble(alp*gamv)
        deltap=dble(bet)
        deltam=dble(betp)
c        if((deltap+1.d0).gt.eps)then
          gamom=gamom/(deltap+1.d0)
c        else
c          gamom=-gamom*log(xminDf)
c        endif
c        if((deltam+1.d0).gt.eps)then
          gamom=gamom/(deltam+1.d0)
c        else
c          gamom=-gamom*log(xminDf)
c        endif
        cint=cint+gamom /(2.d0+deltap) /(2.d0+deltam)
      enddo
      om1intgc=cint

      if(om1intgc.lt.eps)then
        write(*,*) b,deltap,deltam,gamom,cint
        write(*,*) 'WARNING ! om1intgc in epos-omg is <0 !!!!!'
        write(*,*) 'WARNING ! => om1intgc set to 0. !!!!!'
        om1intgc=0.d0
      endif

      return
      end


c----------------------------------------------------------------------
      double precision function om1x(x,gam,del,alp,xr)   !---MC---
c----------------------------------------------------------------------
c function for om1xk with x change to exp(x)
c----------------------------------------------------------------------
      double precision x,xr,gam,del,alp

      om1x=gam*exp(x*del)*(exp(xr)-exp(x))**alp

      return
      end

c----------------------------------------------------------------------
      double precision function om1xI(x,gam,del,alp,xr)   !---MC---
c----------------------------------------------------------------------
c integrated om1x between with x change to exp(x)
c----------------------------------------------------------------------
      double precision x,xr,gam,del,alp,betai,dgammln

      if(abs(alp-1d0).lt.1d-6)then         !alp=1
        om1xI=gam*(exp(xr)*exp(x)**(1.d0+del)/(1.d0+del)
     &            -exp(x)**(2.d0+del)/(2.d0+del))
      else
        om1xI=gam*betai(del+1d0,alp+1d0,exp(x-xr))
     .    *exp(xr*(1d0+del+alp)
     .        +dgammln(1d0+del)+dgammln(1d0+alp)-dgammln(2d0+del+alp))
      endif
      return
      end


c----------------------------------------------------------------------
      double precision function om1xpk(fct,fctI,xp,xpremi,xmremi,k)   !---test---
c----------------------------------------------------------------------
c \int dxm fct(xp,xm)   (normalised)
c k - pair indice;
c----------------------------------------------------------------------
#include "aaa.h"
#include "par.h"
#include "ems.h"

      double precision xp,cint,deltap,deltam,gamomx,x0
     &,gamom,xpremi,xmremi,lxprem,lxmrem,lxp,fct,fctI
      external fct,fctI


      om1xpk=0.d0
      if(xp.gt.xpremi)return
      lxprem=log(xpremi)
      lxmrem=log(xmremi)
      x0=log(xminDf)
      lxp=log(xp)

      imin=ntymin
      imax=ntymax

      cint=0.d0
      do i=imin,imax
        if(imin-imax.eq.0.or.i.ne.ntynd)then
        deltap=btildep(i,k)
        deltam=btildepp(i,k)
        gamom=atilde(i,k)
c integration over xm
        gamom=gamom*(fctI(lxmrem,1d0,deltam,dble(alplea(icltar)),lxmrem)
     .              -fctI(x0,1d0,deltam,dble(alplea(icltar)),lxmrem))
c integration over xp
       gamomx=gamom*(fctI(lxprem,1d0,deltap,dble(alplea(iclpro)),lxprem)
     .              -fctI(x0,1d0,deltap,dble(alplea(iclpro)),lxprem))
        cint=cint+gamomx
        om1xpk=om1xpk+fct(lxp,gamom,deltap,dble(alplea(iclpro)),lxprem)
        endif
      enddo

      om1xpk=om1xpk/cint

      return
      end
c----------------------------------------------------------------------
      double precision function om1xmk(fct,fctI,xp,xm,xpremi,xmremi,k)!---test---
c----------------------------------------------------------------------
c fct(xp,xm) for xp fixed  (normalised)
c k - pair indice;
c----------------------------------------------------------------------
#include "aaa.h"
#include "par.h"
#include "ems.h"

      double precision xp,xm,cint,deltap,deltam,eps,gamomx
     &,gamom,xpremi,xmremi,lxprem,lxmrem,lxm,fct,fctI,lxp,x0
      external fct,fctI

      parameter(eps=1.d-20)


      om1xmk=0.d0
      if(xp.gt.xpremi)return
      if(xm.gt.xmremi)return
      lxprem=log(xpremi)
      lxmrem=log(xmremi)
      x0=log(xminDf)
      lxp=log(xp)
      lxm=log(xm)

      imin=ntymin
      imax=ntymax

      cint=0.d0
      do i=imin,imax
        if(imin-imax.eq.0.or.i.ne.ntynd)then
        deltap=btildep(i,k)
        deltam=btildepp(i,k)
        gamom=fct(lxp,atilde(i,k),deltap,dble(alplea(iclpro)),lxprem)
c integration over xm
        gamomx=fctI(lxmrem,gamom,deltam,dble(alplea(icltar)),lxmrem)
     .        -fctI(x0,gamom,deltam,dble(alplea(icltar)),lxmrem)
        cint=cint+gamomx
        om1xmk=om1xmk+fct(lxm,gamom,deltam,dble(alplea(icltar)),lxmrem)
        endif
      enddo

      om1xmk=om1xmk/cint

      return
      end

c----------------------------------------------------------------------
      double precision function om1xpr(atil,btilp,btilpp
     &                                   ,xpremi,xmremi,ir)   !---MC---
c----------------------------------------------------------------------
c Random number generated from the function om1x. We solve the equation
c which give om1xprk by Newton-Raphson + secant method.
c k - pair indice;
c ir - 1 to get xp, -1 to get xm (be carrefull to inverse xpremi et xmremi
c when calling with ir=-1)
c----------------------------------------------------------------------
#include "aaa.h"
#include "par.h"
#include "ems.h"
      double precision x,x0,x1,gamomx(ntymi:ntymx),lxmrem
      double precision deltap(ntymi:ntymx),lxprem,deltam,gamom
      double precision xpremi,xmremi,xNewtsecant,om1x,om1xI
      double precision atil(ntymi:ntymx),alp
     &                ,btilp(ntymi:ntymx),btilpp(ntymi:ntymx)
      external om1x,om1xI


      om1xpr=0.d0
      lxprem=log(xpremi)
      lxmrem=log(xmremi)
      x0=log(xminDf)
      x1=lxprem
      imin=ntymin
      imax=ntymax

      do i=imin,imax
        if(ir.gt.0)then
          deltap(i)=btilp(i)
          deltam=btilpp(i)
          alp=dble(alplea(iclpro))
          icl=icltar
        else
          deltap(i)=btilpp(i)
          deltam=btilp(i)
          alp=dble(alplea(icltar))
          icl=iclpro
        endif
        gamom=atil(i)
c integration over xm
        gamomx(i)=om1xI(lxmrem,gamom,deltam,dble(alplea(icl)),lxmrem)
     .           -om1xI(x0,gamom,deltam,dble(alplea(icl)),lxmrem)
      enddo

      x=xNewtsecant(om1x,om1xI,x0,x1,gamomx,deltap,alp,imin,imax)

      om1xpr=exp(x)

      return
      end

c----------------------------------------------------------------------
      double precision function om1xprk(fct,fctI,n,k,xpremi,xmremi,iri)!---MC---
c----------------------------------------------------------------------
c Random number generated from the function fct. We solve the equation
c which give om1xpork by Newton-Raphson + secant method.
c k - pair indice; n - position indice;
c ir - 1 to get xp, -1 to get xm (be carrefull to exchange xpremi and xmremi
c when calling with ir=-1)
c----------------------------------------------------------------------
#include "aaa.h"
#include "par.h"
#include "ems.h"
      double precision x,x0,x1,gamomx(ntymi:ntymx),lxmrem,alp
      double precision deltap(ntymi:ntymx),lxprem,deltam,gamom
      double precision xpremi,xmremi,xNewtsecant,fct,fctI,drangen
      external fct,fctI


      om1xprk=0.d0
      if(abs(iri).eq.10)then
        ir=iri/10
      else
        ir=iri
      endif
      lxprem=log(xpremi)
      lxmrem=log(xmremi)
      imin=ntymin
      imax=ntymax
      x0=log(xminDf)
      x1=lxprem

      do i=imin,imax
        if(ir.gt.0)then
          deltap(i)=btildpp(i,n,k)
          deltam=btildppp(i,n,k)
          alp=dble(alplea(iclpro))
          icl=icltar
        else
          deltap(i)=btildppp(i,n,k)
          deltam=btildpp(i,n,k)
          alp=dble(alplea(icltar))
          icl=iclpro
        endif
        gamom=atildp(i,n,k)
c integration over xm
        if(imin-imax.eq.0.or.i.ne.ntynd)then
        gamomx(i)=fctI(lxmrem,gamom,deltam,dble(alplea(icl)),lxmrem)
     .         -fctI(x0,gamom,deltam,dble(alplea(icl)),lxmrem)
      else
        gamomx(i)=0d0
      endif
      enddo
 10   continue
      x=xNewtsecant(fct,fctI,x0,x1,gamomx,deltap,alp,imin,imax)

      if(abs(iri).eq.10.and.
     .   drangen(x).gt.(alp+1d0)*(1d0-exp(x)/xpremi)**alp/xpremi)then
        x0=log(xminDf)
        x1=lxprem
        goto 10
      endif
      om1xprk=exp(x)

      return
      end

c----------------------------------------------------------------------
      double precision function om1xmrk(fct,fctI,n,k,xp
     .                                  ,xpremi,xmremi,iri)   !---MC---
c----------------------------------------------------------------------
c Random number generated from the function fct. We solve the equation
c which give om1xmork by Newton-Raphson + secant method.
c k - pair indice; n - position indice;
c ir - 1 to get xp, -1 to get xm (be carrefull to exchange xpremi and xmremi
c when calling with ir=-1)
c----------------------------------------------------------------------
#include "aaa.h"
#include "par.h"
#include "ems.h"
      double precision x,xp,x0,x1,gamomx(ntymi:ntymx),lxmrem,lxprem
      double precision deltam(ntymi:ntymx),xpremi,lxp,drangen
      double precision xmremi,xNewtsecant,fct,fctI,alp
      external fct,fctI


      om1xmrk=0.d0
      if(xp.gt.xpremi)return
      if(abs(iri).eq.10)then
        ir=iri/10
      else
        ir=iri
      endif
      lxprem=log(xpremi)
      lxmrem=log(xmremi)
      lxp=log(xp)
      x0=log(xminDf)
      x1=lxmrem
      imin=ntymin
      imax=ntymax

      if(ir.gt.0)then
        do i=imin,imax
        if(imin-imax.eq.0.or.i.ne.ntynd)then
          gamomx(i)=fct(lxp,atildp(i,n,k),btildpp(i,n,k)
     .                     ,dble(alplea(iclpro)),lxprem)
        else
          gamomx(i)=0d0
        endif
          deltam(i)=btildppp(i,n,k)
        enddo
        alp=dble(alplea(icltar))
      else
        do i=imin,imax
        if(imin-imax.eq.0.or.i.ne.ntynd)then
          gamomx(i)=fct(lxp,atildp(i,n,k),btildppp(i,n,k)
     .                     ,dble(alplea(icltar)),lxprem)
        else
          gamomx(i)=0d0
        endif
          deltam(i)=btildpp(i,n,k)
        enddo
        alp=dble(alplea(iclpro))
      endif
 10   continue
      x=xNewtsecant(fct,fctI,x0,x1,gamomx,deltam,alp,imin,imax)
      if(abs(iri).eq.10.and.
     .   drangen(x).gt.(alp+1d0)*(1d0-exp(x)/xmremi)**alp/xmremi)then
        x0=log(xminDf)
        x1=lxmrem
        goto 10
      endif

      om1xmrk=exp(x)

      return
      end

c----------------------------------------------------------------------
      double precision function om1xk(xh,xpremi,xmremi,k)   !---MC---
c----------------------------------------------------------------------
c \int dxp om1   (normalised)
c k - pair indice;
c----------------------------------------------------------------------
#include "aaa.h"
#include "par.h"
#include "ems.h"

      double precision xh,gamomx(ntymi:ntymx),cint,alpp(ntymi:ntymx)
     &,delta(ntymi:ntymx),deltap,deltam,eps
     &,gamom,xpremi,xmremi,xprem,xmrem,eymx,eymn
      parameter(eps=1.d-20)


      om1xk=0.d0
      xprem=xpremi !1.d0
      xmrem=xmremi !1.d0
      if(xh.gt.xprem*xmrem)return
      eymx=xprem/sqrt(xh)
      eymn=sqrt(xh)/xmrem

      imin=ntymin
      imax=ntymax

      cint=0.d0
      do i=imin,imax
        gamomx(i)=atilde(i,k)
        deltap=btildep(i,k)
c        if(abs(deltap+1.d0).lt.eps)return
        deltam=btildepp(i,k)
c        if(abs(deltam+1.d0).lt.eps)return

        delta(i)=(deltap+deltam)*0.5d0
        alpp(i)=deltap-deltam
        gamom=gamomx(i)
c integration over xp
        gamom=gamom*xprem**(deltap+1.d0)/(deltap+1.d0)
c integration over xm
        gamom=gamom*xmrem**(deltam+1.d0)/(deltam+1.d0)
        cint=cint+gamom
c integration over y
        if(abs(alpp(i)).gt.eps)then
          om1xk=om1xk+gamomx(i)/alpp(i)*xh**delta(i)
     .                                 *(eymx**alpp(i)-eymn**alpp(i))
        else
          om1xk=om1xk+gamomx(i)*xh**delta(i)*log(eymx/eymn)
        endif
      enddo

      om1xk=om1xk/cint

      return
      end

c----------------------------------------------------------------------
      double precision function om1yk(xh,yp,xpremi,xmremi,k)   !---test---
c----------------------------------------------------------------------
c om1 normalized for fixed xp
c xh - fraction of the energy squared s for the pomeron;
c k - pair indice;
c----------------------------------------------------------------------
#include "aaa.h"
#include "ems.h"
      double precision xh,yp,gamomy(ntymi:ntymx),alpp(ntymi:ntymx),cint
      double precision deltap,deltam,eps,xpremi,xmremi,xprem,xmrem
     .                ,eymx,eymn
      parameter(eps=1.d-20)

      om1yk=0.d0

      xprem=xpremi !1.d0
      xmrem=xmremi !1.d0
      if(xh.gt.xprem*xmrem)return
      eymx=xprem/sqrt(xh)
      eymn=sqrt(xh)/xmrem

      imin=ntymin
      imax=ntymax

      cint=0.d0
      do i=imin,imax
        gamomy(i)=atilde(i,k)
        deltap=btildep(i,k)
        deltam=btildepp(i,k)
        alpp(i)=deltap-deltam
        gamomy(i)=gamomy(i)*xh**((deltap+deltam)*0.5d0)
        if(abs(alpp(i)).gt.eps)then
c integration over y
          cint=cint+gamomy(i)/alpp(i)*(eymx**alpp(i)-eymn**alpp(i))
c not integrated function
          om1yk=om1yk+gamomy(i)*exp(alpp(i)*yp)
        else
          cint=cint+gamomy(i)*log(eymx/eymn)
          om1yk=om1yk+gamomy(i)
        endif
      enddo

      om1yk=om1yk/cint


      return
      end

c----------------------------------------------------------------------
      double precision function om1xrk(n,k,xpremi,xmremi)   !---test---
c----------------------------------------------------------------------
c Random number generated from the function om1xk. We solve the equation
c which give om1xrk by Newton-Raphson + secant method.
c k - pair indice;
c----------------------------------------------------------------------
#include "aaa.h"
#include "ems.h"
#include "par.h"

      double precision x,x0,x1,gamomx(ntymi:ntymx),eps,prec,drangen
      double precision xt,fx,fpx,r,f1,f0,cint,deltx,alpp(ntymi:ntymx)
     &,delta(ntymi:ntymx),deltap(ntymi:ntymx),deltam(ntymi:ntymx)
     &,f0t,f1t,xpremi,xmremi,xprem,xmrem,lxmax,f00!,gamom,cint2
      parameter (eps=1.d-20)


      om1xrk=0.d0

      imin=ntymin
      imax=ntymax
      xprem=xpremi !1.d0
      xmrem=xmremi !1.d0

c      cint2=0.d0
      do i=imin,imax
        gamomx(i)=atildp(i,n,k)
        deltap(i)=btildpp(i,n,k)
        deltam(i)=btildppp(i,n,k)
        delta(i)=(deltap(i)+deltam(i))*0.5d0
        alpp(i)=deltap(i)-deltam(i)
c        gamom=gamomx(i)
c integration over xp
c        gamom=gamom/(deltap(i)+1.d0)
c     .       *(xprem**(deltap(i)+1.d0)-xminDf**(deltap(i)+1.d0))
c integration over xm
c        gamom=gamom/(deltam(i)+1.d0)
c     .       *(xmrem**(deltam(i)+1.d0)-xminDf**(deltam(i)+1.d0))
c        cint2=cint2+gamom
      enddo

      x0=xminDf
      x1=xprem*xmrem
      lxmax=log(x1)
      f0=0.d0
      f1=0.d0
      do i=imin,imax
        if(abs(alpp(i)).lt.eps)then
          f0=f0+gamomx(i)/(delta(i)+1.d0)*x0**(delta(i)+1.d0)
     &        *(lxmax-log(x0)+1.d0/(delta(i)+1.d0))
          f1=f1+gamomx(i)/(delta(i)+1.d0)*x1**(delta(i)+1.d0)
     &        *(lxmax-log(x1)+1.d0/(delta(i)+1.d0))
        else
          f0=f0+gamomx(i)/alpp(i)
     &         *(x0**(deltam(i)+1.d0)*xprem**alpp(i)/(deltam(i)+1.d0)
     &          -x0**(deltap(i)+1.d0)/xmrem**alpp(i)/(deltap(i)+1.d0))
          f1=f1+gamomx(i)/alpp(i)
     &         *(x1**(deltam(i)+1.d0)*xprem**alpp(i)/(deltam(i)+1.d0)
     &          -x1**(deltap(i)+1.d0)/xmrem**alpp(i)/(deltap(i)+1.d0))
        endif
      enddo
      f00=f0
      cint=f1-f00
c      print *,cint,cint2
      f0=-(f0-f00)/cint
      f1=-(f1-f00)/cint
      ntry=0
 11   ntry=ntry+1
      r=drangen(dble(ntry))
      f0t=f0+r
      f1t=f1+r
      if(f1t*f0t.ge.eps.and.ntry.lt.100)goto 11
      if(f1t*f0t.ge.eps)then
        do i=imin,imax
         write(ifmt,*)i,gamomx(i),deltap(i),deltam(i),alpp(i),delta(i)
        enddo
        write(ifmt,*)x0,f0,f0t,x1,f1,f1t,r,cint,ntry,bk(k),k,f0-f1
        call utstop('om1xrk (1)&')
      endif
      f0=f0t
      f1=f1t
c      if(f1*f0.gt.eps)then
c        call utmsg('om1xrk')
c        write(ifch,*)'Poblem with x0, no root ! --> om1xrk=xminDf'
c        write(ifmt,*)'Poblem with x0, no root ! --> om1xrk=xminDf'
c        write(ifmt,*)f0,f1,cint,r
c        call utmsgf
c        om1xrk=x0
c        return
c      endif
      if(abs(f0).lt.eps) then
        om1xrk=x0
        return
      endif
      if(abs(f1).lt.eps) then
        om1xrk=x1
        return
      endif
c      x=(x1+x0)*0.5D0
      x=sqrt(x1*x0)
      deltx=abs(x1-x0)

      ntry=0
      fx=0.d0
      fpx=0.d0
      xt=x
 111  continue

      if(ntry.le.1000)then
      fx=0.d0
      fpx=0.d0
      do i=imin,imax
        if(abs(alpp(i)).lt.eps)then
          fx=fx+gamomx(i)/(delta(i)+1.d0)*x**(delta(i)+1.d0)
     &         *(lxmax-log(x)+1.d0/(delta(i)+1.d0))
          fpx=fpx+gamomx(i)*x**delta(i)*(lxmax-log(x))
        else
          fx=fx+gamomx(i)/alpp(i)
     &         *(x**(deltam(i)+1.d0)*xprem**alpp(i)/(deltam(i)+1.d0)
     &          -x**(deltap(i)+1.d0)/xmrem**alpp(i)/(deltap(i)+1.d0))
          fpx=fpx+gamomx(i)/alpp(i)*(x**deltam(i)*xprem**alpp(i)
     &                              -x**deltap(i)/xmrem**alpp(i))
        endif
      enddo
      fx=-(fx-f00)/cint+r
      fpx=fpx/cint
      xt=x-fx/fpx

      if (f0*fx.lt.-eps) then
        f1=fx
        x1=x
      else
        f0=fx
        x0=x
      endif
      if ((xt.lt.x0-eps.or.xt.gt.x1+eps).and.abs(f1-f0).gt.eps) then
        xt=x1-f1*(x1-x0)/(f1-f0)
      endif

      else

        xt=0d0
        write(ifmt,*)'Warning in om1xrk, to much try !'

      endif

      if(abs(x-xt).gt.deltx*0.5d0) then
        xt=sqrt(x1*x0)
      endif
      deltx=abs(x-xt)
      if(abs(x).gt.eps)then
        prec=deltx/x
      else
        prec=0d0
        call utstop('Problem in om1xrk&')
      endif

      if (prec.gt.1.d-3.and.abs(f1-f0).gt.eps.and.ntry.le.1000)then
         x=xt
         ntry=ntry+1
         goto 111
      endif

      om1xrk=x

      return
      end

c----------------------------------------------------------------------
      double precision function om1yrk(n,k,xh,xpremi,xmremi)   !---test---
c----------------------------------------------------------------------
c Random number generated from the function om1yk(xh). We solve the
c equation which give om1yrk by Newton-Raphson + secant method.
c xh - fraction of the energy squared s for the pomeron;
c k - pair indice;
c----------------------------------------------------------------------
#include "aaa.h"
#include "ems.h"

      double precision x,x0,x1,gamomy(ntymi:ntymx),eps,prec,drangen
      double precision xt,fx,fpx,f1,f0,cint,deltx,alpp(ntymi:ntymx)
     &,deltap,deltam,f00
     &,f0t,f1t!,gamom,cint2
      parameter (eps=1.d-20)
      double precision xh,r,xpremi,xmremi,xprem,xmrem!,lnx,lnxpr,lxmr

c      lnxpr=log(xpremi) !0.d0
c      lnxmr=log(xmremi) !0.d0
c      lnx=log(xh)

c      r=dble(rangen())

cc      om1yrk=(0.5d0-r)*log(xh)
c      om1yrk=0.5d0*lnx-lnxmr+r*(lnxpr+lnxmr-lnx)


      om1yrk=0d0

      xprem=xpremi !1.d0
      xmrem=xmremi !1.d0
      if(xh.gt.xprem*xmrem)return
      eymx=xprem/sqrt(xh)
      eymn=sqrt(xh)/xmrem

      imin=ntymin
      imax=ntymax

c      cint2=0.d0
      do i=imin,imax
        gamomy(i)=atildp(i,n,k)
        deltap=btildpp(i,n,k)
        deltam=btildppp(i,n,k)
        alpp(i)=deltap-deltam
        gamomy(i)=gamomy(i)*xh**((deltap+deltam)*0.5d0)
c        if(abs(alpp(i)).ge.eps)then
cc integration over y
c          cint2=cint+gamomy(i)/alpp(i)*(eymx**alpp(i)-eymn**alpp(i))
c        else
c          cint2=cint+gamomy(i)*log(eymx/eymn)
c        endif
      enddo

      x0=log(eymn)
      x1=log(eymx)
      f0=0.d0
      f1=0.d0
      do i=imin,imax
        if(abs(alpp(i)).lt.eps)then
          f0=f0+gamomy(i)*x0
          f1=f1+gamomy(i)*x1
        else
          f0=f0+gamomy(i)/alpp(i)*exp(x0*alpp(i))
          f1=f1+gamomy(i)/alpp(i)*exp(x1*alpp(i))
        endif
      enddo
      f00=f0
      cint=f1-f00
c      print *,cint,cint2
      f0=-(f0-f00)/cint
      f1=-(f1-f00)/cint
      ntry=0
 11   ntry=ntry+1
      r=drangen(dble(ntry))
      f0t=f0+r
      f1t=f1+r
      if(f1t*f0t.ge.eps.and.ntry.lt.100)goto 11
      if(f1t*f0t.ge.eps)then
        do i=imin,imax
         write(ifmt,*)i,gamomy(i),alpp(i)
        enddo
        write(ifmt,*)x0,f0,f0t,x1,f1,f1t,r,cint,ntry,bk(k),k
        call utstop('om1yrk (1)&')
      endif
      f0=f0t
      f1=f1t
c      if(f1*f0.gt.eps)then
c        call utmsg('om1xrk')
c        write(ifch,*)'Poblem with x0, no root ! --> om1yrk=ymin'
c        write(ifmt,*)'Poblem with x0, no root ! --> om1yrk=ymin'
c        write(ifmt,*)f0,f1,cint,r
c        call utmsgf
c        om1yrk=x0
c        return
c      endif
      if(abs(f0).lt.eps) then
        om1yrk=x0
        return
      endif
      if(abs(f1).lt.eps) then
        om1yrk=x1
        return
      endif
      x=(x1+x0)*0.5D0
c      x=sqrt(x1*x0)
      deltx=abs(x1-x0)

      ntry=0
      fx=0.d0
      fpx=0.d0
      xt=x
 111  continue

      if(ntry.le.1000)then
      fx=0.d0
      fpx=0.d0
      do i=imin,imax
        if(abs(alpp(i)).lt.eps)then
          fx=fx+gamomy(i)*x
          fpx=fpx+gamomy(i)
        else
          fx=fx+gamomy(i)/alpp(i)*exp(x*alpp(i))
          fpx=fpx+gamomy(i)*exp(x*alpp(i))
        endif
      enddo
      fx=-(fx-f00)/cint+r
      fpx=fpx/cint
      xt=x-fx/fpx

      if (f0*fx.lt.-eps) then
        f1=fx
        x1=x
      else
        f0=fx
        x0=x
      endif
      if ((xt.lt.x0-eps.or.xt.gt.x1+eps).and.abs(f1-f0).gt.eps) then
        xt=x1-f1*(x1-x0)/(f1-f0)
      endif

      else

        xt=0d0
        write(ifmt,*)'Warning in om1xrk, to much try !'

      endif

      if(abs(x-xt).gt.deltx*0.5d0) then
c        xt=sqrt(x1*x0)
        xt=(x1+x0)*0.5D0
      endif
      deltx=abs(x-xt)
      if(abs(x).gt.eps)then
        prec=deltx/x
      else
        prec=deltx
      endif

      if (prec.gt.1.d-3.and.abs(f1-f0).gt.eps.and.ntry.le.1000)then
         x=xt
         ntry=ntry+1
         goto 111
      endif

      om1yrk=x

      return
      end

c----------------------------------------------------------------------
      double precision function om1xo(x,gam,del,alp,xr)   !---MC---
c----------------------------------------------------------------------
c function for om1xok with x change to exp(x)
c----------------------------------------------------------------------
      double precision x,xr,gam,del,dum

      dum=xr
      dum=alp
      om1xo=gam*exp(x*del)

      return
      end

c----------------------------------------------------------------------
      double precision function om1xoI(x,gam,del,alp,xr)   !---MC---
c----------------------------------------------------------------------
c integrated om1xo between with x change to exp(x)
c----------------------------------------------------------------------
      double precision x,xr,gam,del,dum

      dum=xr
      dum=alp
      om1xoI=gam/(del+1d0)*exp(x*(del+1d0))

      return
      end

c----------------------------------------------------------------------
      double precision function xNewtsecant(fct,fctI,x0,x1
     .                                     ,gam,del,alp,imin,imax)   !---MC---
c----------------------------------------------------------------------
c Random number generated from the function fct. We solve the equation
c by Newton-Raphson + secant method.
c fct - function to be followed by random numbers 
c       (form fct=sum on i from imin to imax of 
c       (gam(i)*x**del(i)*g(alp,x,x0,x1))
c        where g(alp,x,x0,x1) can be any regular function of alp, x, x0 and x1)
c fctI - integral of fct
c x0 - minimum value for random number
c x1 - maximum value for random number
c----------------------------------------------------------------------
#include "aaa.h"
#include "ems.h"
      double precision x,x0,x1,x2,gam(ntymi:ntymx),eps,f0t,f1t,f00
      double precision xt,fx,fpx,r,f1,f0,cint,deltx,prec,drangen
      double precision del(ntymi:ntymx),alp
      parameter (eps=1.d-20)
      double precision fct,fctI
      external fct,fctI

      xNewtsecant=0d0
      x2=x1
      f0=0.d0
      f1=0.d0
      do i=imin,imax
        f0=f0+fctI(x0,gam(i),del(i),alp,x2)
        f1=f1+fctI(x1,gam(i),del(i),alp,x2)
      enddo
      f00=f0
      cint=f1-f00
      f0=-(f0-f00)/cint
      f1=-(f1-f00)/cint
      ntry=0
 11   ntry=ntry+1
      r=drangen(dble(ntry))
      f0t=f0+r
      f1t=f1+r
      if(f1t*f0t.ge.eps.and.ntry.lt.100)goto 11
      if(f1t*f0t.ge.eps)then
        do i=imin,imax
         write(ifmt,*)i,gam(i),del(i)
        enddo
        write(ifmt,*)x0,f0,f0t,x1,f1,f1t,r,cint,x2,ntry
        call utstop('xNewtsecant not possible !&')
      endif
      f0=f0t
      f1=f1t
      if(abs(f0).le.eps) then
        xNewtsecant=x0
        return
      endif
      if(abs(f1).le.eps) then
        xNewtsecant=x1
        return
      endif
      x=0.5d0*(x1+x0)
      deltx=abs(x1-x0)


      ntry=0

 111  continue
      if(ntry.le.1000)then
      fx=0.d0
      fpx=0.d0
      do i=imin,imax
        fx=fx+fctI(x,gam(i),del(i),alp,x2)
        fpx=fpx+fct(x,gam(i),del(i),alp,x2)
      enddo
      fx=-(fx-f00)/cint+r
      fpx=fpx/cint
      xt=x-fx/fpx

      if (f0*fx.lt.0.D0) then
        f1=fx
        x1=x
      else
        f0=fx
        x0=x
      endif
      if ((xt.lt.x0-eps.or.xt.gt.x1+eps).and.abs(f1-f0).gt.eps) then
        xt=x1-f1*(x1-x0)/(f1-f0)
      endif

       else

        xt=0.
        write(ifmt,*)'Warning in xNewtsecant, to much try !'

      endif


      if(abs(x-xt).gt.deltx*0.5d0) then
        xt=(x1+x0)*0.5D0
      endif
      deltx=abs(x-xt)
      if(abs(x).gt.eps)then
        prec=abs(deltx/x)
      else
        prec=0d0
        call utstop('Problem in xNewtsecant...&')
      endif

      if (prec.gt.1.d-3.and.abs(f1-f0).gt.eps.and.ntry.le.1000) then
         x=xt
         ntry=ntry+1
         goto 111
      endif

      xNewtsecant=x

      return
      end


c----------------------------------------------------------------------
      double precision function om51(xh,yp,b,iqq1,iqq2)   !---test---
c----------------------------------------------------------------------
c xh - xplus*xminus     iqq=-1 ... fit        (om1)
c yp - rapidity         iqq=0 .... soft
c b - impact param      iqq=1 .... gg
c iq1 - min iq          iqq=2 .... qg
c iq2 - max iq          iqq=3 .... gq
c                       iqq=4 .... qq
c                       iqq=5 .... R   (reggeon exchange : pion exchange
c                                       pro and/or targ excited)
c                       iqq=6 .... Y+  (no-connection possible on proj)
c                       iqq=7 .... Y-  (no-connection possible on targ)
c                       iqq=8 .... X   (no-connection possible on pro and tar)
c                       iqq=9 .... DD  (excitation of pro and tar)
c                       for fit (with Chad*x**alpdif factor)
c                       iqq=-5 ... all diff
c                       iqq=-6 ... SD- (pro non excited, mass to targ)
c                       iqq=-7 ... SD+ (targ non excited, mass to pro)
c                       iqq=-8 ... CD
c                       iqq=-9 ... DD
c                       iqq=11 ... soft + gg born
c                       iqq=11 ... soft + gg born (bottom)
c----------------------------------------------------------------------
#include "aaa.h"
#include "sem.h"
#include "par.h"
      double precision xp,xm,xh,yp,om51p,om1,om1di,fom51d,om0
     .                ,zb,rp,zbd,rd

      om51=0.d0
      sy=engy*engy
      if(xh.gt.1.d0)return

c      if(sy*xh.le.4.*q2zmin)
c     . .and.iqq1.gt.0.and.iqq2.lt.5)return
      xp=sqrt(xh)*exp(yp)
      xm=xh/xp

      iqq=min(iqq1,iqq2)
      iqqx=max(iqq1,iqq2)
      if(iqq.eq.-1)then
        om51=om1(xh,yp,b)
c        om51=0.5d0*om1(xh,yp,b)
        return
      elseif((iqqx.ge.0.and.abs(iqq).le.5).or.iqq.eq.11
     .                                    .or.iqq.eq.12)then
        i2=max(abs(iqq1),abs(iqq2))
        i1=min(abs(iqq1),abs(iqq2))
        do i=i1,i2
          om51=om51+om51p(sy*sngl(xh),xh,yp,b,i)
        enddo
      endif
      if(iqq.lt.12.and.iomega.lt.2.and.iqqx.gt.5.or.iqq.le.-5)then
        rp=dble(r2had(iclpro)+r2had(icltar)
     .       +slopom*log(max(1.,sngl(dble(sy)*xh))))
        zb=exp(-dble(b)**2/(4.d0*.0389d0*rp))
        om0=om51p(sngl(dble(sy)*xh),xh,yp,b,-5) !soft Pomeron for diff
        if(iqq.lt.0)then
          if(iqq.eq.-5)then
            imin=-9
            imax=-6
          else
            imax=min(-6,iqqx)
            imin=max(-9,iqq)
          endif
        else
          imin=max(6,iqq)
          imax=min(9,iqqx)
        endif
c        write(ifch,*)'om51',sy,xh,yp,b,rp,iqq1,iqq2,imin,imax,om0
        do i=imin,imax
          om1di=om0*fom51d(dble(sy),xp,xm,zb,rp,zbd,rd,i)
          om51=om51+om1di
        enddo
      endif
c      print *,xh,yp,b,iqq1,iqq2,om51

      end

c----------------------------------------------------------------------
      double precision function fom51d(sy,xp,xm,zb,rp,zbd,rd,iqq)   !---MC---
c----------------------------------------------------------------------
c Factor to apply to G_soft to get diffractive amplitude for given
c square energy sy, light cone momenta, x+ and x- and impact parameter
c b (zb=exp(-b**2/(4.*.0389*rp)) where rp is given for soft int.).
c iqq=6 .... Y+  (no-connection possible on proj) (SD-)
c iqq=7 .... Y-  (no-connection possible on targ) (SD+)
c iqq=8 .... X   (no-connection possible on pro and tar) (CD)
c iqq=9 .... loop  (excitation of pro and tar) (DD)
c----------------------------------------------------------------------
#include "aaa.h"
#include "sem.h"
#include "ems.h"
#include "par.h"
      double precision xp,xm,rd,zb,zbd,rp,delp,sy,rom51r

      fom51d=0.5d0*r3pomti*rp/zb
      delp=dble(1.-2.*alppom+alppar)
      if(abs(iqq).eq.6)then            !Y+
        rd=rp-dble(2.*slopom)*log(xp)+r2hads(1)
        zbd=zb**(rp/rd)
        fom51d=fom51d*dble(gampar*wdiff(iclpro)/gamhad(iclpro))
     .               *xp**delp*(xm*dble(s0min))**dble(1.-alppom)
     .               *max(0d0,1d0-xp/dble(xmxrem(1))) !large xp suppressed
     .               *zbd/rd
      elseif(abs(iqq).eq.7)then        !Y-
        rd=rp-dble(2.*slopom)*log(xm)+r2hads(1)
        zbd=zb**(rp/rd)
        fom51d=fom51d*dble(gampar*wdiff(icltar)/gamhad(icltar))
     .               *xm**delp*(xp*dble(s0min))**dble(1.-alppom)
     .               *max(0d0,1d0-xm/dble(xmxrem(2))) !large xm suppressed
     .               *zbd/rd
      else                       !X or loop
        rd=rp-dble(2.*slopom)*log(xp*xm)+r2hads(2)
        zbd=zb**(rp/rd)
        fom51d=fom51d*r3pomti*(xp*xm)**delp
     .  *dble(s0min**(1-alppom))*zbd/rd

        if(abs(iqq).eq.8 .and. xmxrem(1)/=0. .and. xmxrem(2)/=0.)then !X !VV avoid division by zero 
        !if(abs(iqq).eq.8)then         !X
          fom51d=0.5d0*fom51d*dble(gampar**2)
     .  *dble(wdiff(iclpro)/gamhad(iclpro)*wdiff(icltar)/gamhad(icltar))
     .             *max(0d0,1d0-xp/dble(xmxrem(1)))
     .             *max(0d0,1d0-xm/dble(xmxrem(2)))
c        else                     !loop 
ctp2016 With this factor it doesn't work: missing integral cancelling s**dels????
c          fom51d=fom51d*dble((sy/s0min)**(alppom-1.))
        endif
      endif
      fpomd=1d0
c better to keep change in both xp and xm component to keep correct rapidity shape independently of alpdif
      if(iqq.lt.0)fpomd=fpomd*(xp*xm)**alpdif
     .                 *dble(chadr(abs(iqq)-5,iclpro,icltar))

c      write(ifch,*)'fom51d',iqq,xp,xm,zb,iqq,fom51d,fpomd

      rom51r=sy
c      if(iregge.ne.0)then
c        fom51d=fom51d*(fpomd+rom51r(sy,xp,xm,zbd,rd,iqq)) !Reggeon fraction
c      else
        fom51d=fom51d*fpomd
c      endif


      end

c----------------------------------------------------------------------
      double precision function rom51r(sy,xp,xm,zbd,rd,iqq)   !---MC---
c----------------------------------------------------------------------
c Factor to apply to G_diff to get Reggeon amplitude for given
c square energy sy, light cone momenta, x+ and x- and impact parameter
c b (zbd=exp(-b**2/(4.*.0389*rd)) where rd is given for diff int.).
c iqq=6 .... Y+  (no-connection possible on proj) (SD-)
c iqq=7 .... Y-  (no-connection possible on targ) (SD+)
c iqq=8 .... X   (no-connection possible on pro and tar) (CD)
c iqq=9 .... loop  (excitation of pro and tar) (DD)
c----------------------------------------------------------------------
#include "aaa.h"
#include "sem.h"
#include "ems.h"
#include "par.h"
      double precision xp,xm,rd,zbd,rr,delp,sy

      rom51r=(sy/dble(s0min))**dble(alpreg(iclreg)-alppom)*rd
     .      /dble(r3pom)
      delp=dble(alpreg(iclreg)+alppom-1.)
      rr=dble(sloreg(iclreg))
      if(abs(iqq).eq.6)then            !Y+
        rr=dble(r2reg(iclpro,ichproj)+r2had(icltar))
     .    +rr*log(max(1d0,xp*sy/dble(s0min)))
        rom51r=rom51r*dble(gamreg(iclpro,ichproj))
     .        *zbd**(rd/rr)/rr*xp**delp
      elseif(abs(iqq).eq.7)then        !Y-
        rr=dble(r2reg(icltar,ichtarg)+r2had(iclpro))
     .    +rr*log(max(1d0,xm*sy/dble(s0min)))
        rom51r=rom51r*dble(gamreg(icltar,ichtarg))
     .        *zbd**(rd/rr)/rr*xm**delp
      else                       !X or loop
        rr=dble(r2reg(iclpro,ichproj)+r2reg(icltar,ichtarg))
     .    +rr*log(max(1d0,xp*xm*sy/dble(s0min)))
        rom51r=rom51r*dble(gamreg(iclpro,ichproj)
     .                    *gamreg(icltar,ichtarg)/r3pom)
     .        *zbd**(rd/rr)/rr*(xp*xm)**delp
      endif

c      write(ifch,*)'rom51r',iqq,xp,xm,zbd,rr,rd/rr,rom51r

      end


c----------------------------------------------------------------------
      double precision function om5Jk(n,k,iqq)   !---MC---
c----------------------------------------------------------------------
c partial om5
c n     - Pomeron index
c k     - Pair index
c iqq=-4- fit SD
c iqq=-3- fit non-diff
c iqq=-2- fit DD
c iqq=-1- fit
c iqq=0 - soft
c iqq=1 - gg
c iqq=2 - qg
c iqq=3 - gq
c iqq=4 - qq
c iqq=5 - Reggeon
c iqq=-5 - all diff for fit
c iqq=-50- all diff for MC
c iqq=6 .... SD- (Y+ pro non excited, mass to targ, real)
c iqq=-6 ... SD- (Y+ pro non excited, mass to targ, fit)
c iqq=-60 .. SD- (Y+ pro non excited, mass to targ, MC)
c iqq=7 .... SD+ (Y- targ non excited, mass to pro, real)
c iqq=-7 ... SD+ (Y- targ non excited, mass to pro, fit)
c iqq=-70 .. SD+ (Y- targ non excited, mass to pro, MC)
c iqq=8 .... CD  (X soft : pro and targ non excited
c                 mass to resonance or string, real)
c iqq=-8 ... CD  (X soft : pro and targ non excited
c                 mass to resonance or string, fit)
c iqq=-80 .. CD  (X soft : pro and targ non excited
c                 mass to resonance or string, MC)
c iqq=9 .... DD  (loop mass to pro and targ, real)
c iqq=-9 ... DD  (loop mass to pro and targ, fit)
c iqq=-90 .. DD  (loop mass to pro and targ, MC)
c iqq=11 ... - soft with gg born
c iqq=12 ... - soft with gg born (bottom)
c iqq=100 .. real soft with screening soft
c iqq=101 .. real gg with screening hard
c iqq=102 .. real qg with screening hard
c iqq=103 .. real gq with screening hard
c iqq=104 .. real qq with screening hard
c Screening correction applied to fit and diff diagrams (but not to others except for 100)
c----------------------------------------------------------------------
#include "aaa.h"
#include "ems.h"
#include "sem.h"
#include "par.h"

      double precision xh,yp,om51,omGamk,xp,xm,xpsup,xmsup

      xp=xppr(n,k)
      xm=xmpr(n,k)
      xh=xpr(n,k)
      yp=ypr(n,k)
      b=bk(k)

      if(iqq.le.-50)then
        xpsup=max(0d0,1d0-xp/dble(xmxrem(1))) !large xp suppressed
        xmsup=max(0d0,1d0-xm/dble(xmxrem(2))) !large xm suppressed
        om5Jk=0d0
        if(iqq.eq.-60.or.iqq.eq.-50)om5Jk=om5Jk+omGamk(n,k,xp,xm,2,2)
     .                      *xpsup
        if(iqq.eq.-70.or.iqq.eq.-50)om5Jk=om5Jk+omGamk(n,k,xp,xm,3,3)
     .                      *xmsup
        if(iqq.eq.-80.or.iqq.eq.-50)om5Jk=om5Jk+omGamk(n,k,xp,xm,0,0)
     .                      *0.5d0*dble(gampar**2)*xpsup*xmsup
     .                      /(1d0+0.5d0*dble(gampar**2)*xpsup*xmsup)
        if(iqq.eq.-90.or.iqq.eq.-50)om5Jk=om5Jk+omGamk(n,k,xp,xm,0,0)
     .                      /(1d0+0.5d0*dble(gampar**2)*xpsup*xmsup)
      elseif(iqq.eq.-1)then
        om5Jk=omGamk(n,k,xp,xm,ntymin,ntymax)
      elseif(iqq.eq.-2)then
        om5Jk=omGamk(n,k,xp,xm,0,0)
      elseif(iqq.eq.-3)then
        om5Jk=omGamk(n,k,xp,xm,1,1)
      elseif(iqq.eq.-4)then
        om5Jk=omGamk(n,k,xp,xm,2,3)
      else
        if(iqq.lt.100)then
          om5Jk=om51(xh,yp,b,iqq,iqq)
        else
          om5Jk=om51(xh,yp,b,iqq-100,iqq-100)
        endif
c screening in Diffraction
        if(iqq.ge.100.or.
     .    (abs(iqq).ge.5.and.iqq.ne.5.and.iqq.ne.11.and.iqq.ne.12))then
c          if(iqq.ne.0.and.iqq.lt.100)om5Jk=om5Jk*dble(rexddf)
          if(abs(iqq).eq.6)then
            i=2
          elseif(abs(iqq).eq.7)then
            i=3
          elseif(iqq.eq.100)then
            i=1
          elseif(iqq.gt.100)then    !special screening for hard
            i=-1
          else
            i=0
          endif
          om5Jk=om5Jk*xp**(epsilongp(i,k))*xm**(epsilongt(i,k))
     .           /XminDf**(epsilongs(i,k))
        endif
      endif

      return
      end

c----------------------------------------------------------------------
      double precision function om5J(xp,xm,b,iqq1,iqq2) !---xs---
c----------------------------------------------------------------------
c                       iqq=-4 ... fit SD        (call of Gfunpar needed)
c                       iqq=-3 ... fit non-diff
c                       iqq=-2 ... fit DD
c                       iqq=-1 ... fit all
c xp - xplus            iqq=-10 .. fit all        (om1)
c xm - xminus           iqq=0 .... soft
c b - impact param      iqq=1 .... gg
c                       iqq=2 .... qg
c                       iqq=3 .... gq
c                       iqq=4 .... qq
c                       iqq=5 .... Reggeon
c                       iqq=-5 ... all diff(soft (6 to 9) for fit
c                       iqq=-50 .. all diff for MC
c                       iqq=6 .... SD- (pro non excited, mass to targ, real)
c                       iqq=-6 ... SD- (pro non excited, mass to targ, fit)
c                       iqq=-60 .. Y+ diff MC
c                       iqq=7 .... SD+ (targ non excited, mass to pro, real)
c                       iqq=-7 ... SD+ (targ non excited, mass to pro, fit)
c                       iqq=-70 .. Y- diff MC
c                       iqq=8 .... CD  (soft : pro and targ non excited
c                                       mass to resonance or string, real)
c                       iqq=-8 ... CD  (soft : pro and targ non excited
c                                       mass to resonance or string, fit)
c                       iqq=-80 .. X diff MC (approximation)
c                       iqq=9 .... DD  (mass to pro and targ, real)
c                       iqq=-9 ... DD  (mass to pro and targ, fit)
c                       iqq=-90 .. loop diff MC (approximation)
c                       iqq=100 .. real soft with screening soft
c                       iqq=101 .. real gg with screening hard
c                       iqq=102 .. real qg with screening hard
c                       iqq=103 .. real gq with screening hard
c                       iqq=104 .. real qq with screening hard
c Screening correction applied to diff diagrams and 100
c r2hads used to compensate approximation on calculation of diff xs (no MPI)
c Warning, Gfunpar has to be called for all diagrams before
c----------------------------------------------------------------------
#include "aaa.h"
#include "ems.h"
#include "sem.h"
#include "par.h"
      double precision xp,xm,xh,yp,om51,om1f,xpsup,xmsup
      om5J=0.d0
      i=0
      xh=xp*xm
      if(xp.gt.0.d0.and.xm.gt.0.d0)then
        yp=0.5d0*log(xp/xm)
      else
        yp=0.d0
      endif
      iqq=min(iqq1,iqq2)
      if(iqq.eq.-10)then
        om5J=om51(xh,yp,b,-1,-1)
      elseif(iqq.eq.-1)then
        om5J=om1f(xp,xm,ntymin,ntymax)
      elseif(iqq.eq.-2)then
        om5J=om1f(xp,xm,0,0)
      elseif(iqq.eq.-3)then
        om5J=om1f(xp,xm,1,1)
      elseif(iqq.eq.-4)then
        om5J=om1f(xp,xm,2,3)
      elseif(iqq.le.-50)then
        xpsup=max(0d0,1d0-xp/dble(xmxrem(1))) !large xp suppressed
        xmsup=max(0d0,1d0-xm/dble(xmxrem(2))) !large xm suppressed
        om5J=0d0
        if(iqq.eq.-60.or.iqq.eq.-50)om5J=om5J+om1f(xp,xm,2,2)*xpsup
        if(iqq.eq.-70.or.iqq.eq.-50)om5J=om5J+om1f(xp,xm,3,3)*xmsup
        if(iqq.eq.-80.or.iqq.eq.-50)om5J=om5J+om1f(xp,xm,0,0)
     .                      *0.5d0*dble(gampar**2)*xpsup*xmsup
     .                     /(1d0+0.5d0*dble(gampar**2)*xpsup*xmsup)
        if(iqq.eq.-90.or.iqq.eq.-50)om5J=om5J+om1f(xp,xm,0,0)
     .                     /(1d0+0.5d0*dble(gampar**2)*xpsup*xmsup)
c        om5J=om5J*dble(rexddf)
      else
        if(iqq.lt.100)om5J=om51(xh,yp,b,iqq1,iqq2)
c screening in Diffraction
        if(iqq.ge.100.or.
     .    (abs(iqq).ge.5.and.iqq.ne.11.and.iqq.ne.12))then
c          if(iqq.ne.0.and.iqq.lt.100)om5J=om5J*dble(rexddf)
          if(abs(iqq).eq.6)then
            i=2
          elseif(abs(iqq).eq.7)then
            i=3
          elseif(iqq.eq.100)then
            i=1
            om5J=om51(xh,yp,b,0,0)
          elseif(iqq.gt.100)then    !special screening for hard
            i=-1
            om5J=om51(xh,yp,b,iqq1-100,iqq1-100)
          else
            i=0
          endif
          om5J=om5J*xp**epspUni(i)*xm**epstUni(i)
     .         /XminDf**epssUni(i)
        endif
      endif


      end

c----------------------------------------------------------------------
      double precision function om1f(xp,xm,iq1,iq2)   !---MC---
c----------------------------------------------------------------------
c om1f = G from fit
c xp - fraction of the energy squared s for the pomeron on projectile side
c xm - fraction of the energy squared s for the pomeron on target side
c iq1 and iq2 to select one component (0-soft, 1-idxD-sh, idxD1-diff)
c----------------------------------------------------------------------
#include "aaa.h"
#include "sem.h"
#include "par.h"
      double precision xp,xm

      om1f=0.d0
      imin=max(idxDmin(iomega),iq1)
      imax=min(idxDmax(iomega),iq2)
      if(alpUni(0,1).le.0.)then
        if(imin.eq.0.and.imax.eq.0)then
          return
        elseif(imin.eq.1)then
          imin=0
        endif
      endif
      do i=imin,imax
        om1f=om1f+alpUni(i,1)*xp**betUni(i,1)*xm**betpUni(i,1)
      enddo
c      om1f=om1f*0.5d0
      end

c----------------------------------------------------------------------
      double precision function om1fJ(xp,xm,iq1,iq2)   !---test---
c----------------------------------------------------------------------
c om1fJ = G from fit
c xp - fraction of the energy squared s for the pomeron on projectile side
c xm - fraction of the energy squared s for the pomeron on target side
c iq1 and iq2 to select one component (0-soft, 1-idxD-sh, idxD1-diff)
c----------------------------------------------------------------------
#include "aaa.h"
#include "sem.h"
#include "par.h"
      double precision xp,xm,alp,bet,betp!,om1f,ww,om1ft

      om1fJ=0.d0
c      if(idxD0.ne.0.or.idxD1.ne.2) stop "Problem in om1fJ"
      imin=max(idxDmin(iomega),iq1)
      imax=min(idxDmax(iomega),iq2)
c      if(iscreen.ne.0)then
c        ww=2d0*om1f(xp,xm,0,1)
c        alp=alpUni(1,1)/engy**(2*epssUni(1))
c        bet=betUni(1,1)-epspUni(1)
c        betp=betpUni(1,1)-epstUni(1)
c        om1ft=alp*xp**bet*xm**betp     !hard without screening
c      endif
      if(alpUni(0,1).le.0.)then
        if(imin.eq.0.and.imax.eq.0)then
          return
        elseif(imin.eq.1)then
          imin=0
        endif
      endif
      do i=imin,imax
        alp=alpUni(i,1)
        bet=betUni(i,1)
        betp=betpUni(i,1)
c screening correction for Pomeron type
c        if(iscreen.ne.0.and.i.le.1)then
c          if(i.eq.0)then
c            om1fJ=om1fJ+max(0d0,ww-om1ft)
c          else                  !no screening for hard component
c            om1fJ=om1fJ+min(ww,om1ft)
c          endif
c        else
          om1fJ=om1fJ+alp*xp**bet*xm**betp
c        endif
      enddo
c      om1fJ=om1fJ*0.5d0
      end


ctp-------------------------------------------------------------
        double precision function Fgfactor(xp,xm,iqq)
ctp-------------------------------------------------------------
      double precision xp,xm
      idummy=iqq
      dummy=xp
      dummy=xm
      Fgfactor=1d0
      end

ctp-------------------------------------------------------------
        double precision function Fgfactork(n,k,xp,xm)   !---MC---
ctp-------------------------------------------------------------
      double precision xp,xm
      dummy=xp
      dummy=xm
      dummy=n*k
      Fgfactork=1d0
      end

c----------------------------------------------------------------------
      double precision function omIgamint(b,iqq)   !---test---
c----------------------------------------------------------------------
c - integrated GFF
c b - impact parameter between the pomeron ends;
c iqq=0 effective DD + CD
c iqq=1 hard
c iqq=2 Y+
c iqq=3 Y-
c----------------------------------------------------------------------
#include "aaa.h"
#include "ems.h"
#include "sem.h"
#include "par.h"

      double precision Df

      Df=0.d0
      sy=engy*engy
      omIgamint=0.d0
      imin=idxDmin(iomega)
      imax=idxDmax(iomega)

      if(iqq.ge.idxD0.and.iqq.le.idxD1)then
        imin=iqq
        imax=iqq
      endif
      coefp=1.+alplea(iclpro)
      coeft=1.+alplea(icltar)

      do i=imin,imax
        zp=0.
        zt=0.
        call Gfunpar(zp,zt,0.,0.,1,i,b,sy,alpx,betx,betpx,epsp,epst,epss
     .  ,gamv)
        betp=1.+betx
        betpp=1.+betpx
        Df=alpx*dble(utgam1(betp)*utgam1(betpp)*ucfpro
     *         *ucftar/utgam1(betp+coefp)/utgam1(betpp+coeft))
        omIgamint=omIgamint+Df
      enddo

c      omIgamint=omIgamint

c      omIgamint=omIgamint*0.5d0

      return
      end

cc-----------------------------------------------------------------------
c      subroutine WomTyDif(w,n,k)
cc-----------------------------------------------------------------------
cc - w(ity) for diff or all Pomeron
cc the probability of the type of the same final state.
cc k - pair indice; n - position indice;
cc    ity = 0   - CD+DD = fit - (soft+hard) - Reggeon - SD
cc        = 1   - diff from MC = fit (diff+non-diff) - (soft+hard) - Reggeon
cc        = 2   - SD- (Y+) fraction for soft+hard
cc        = 3   - SD+ (Y-) fraction for soft+hard
cc        = 4   - DPE (X)  fraction for soft+hard
cc        = 5   - Reggeon
cc        = 6   - soft SD- (mass to targ) from fit
cc        = 7   - soft SD+ (mass to proj) from fit
cc        = 8   - soft CD  (no mass)
cc        = 9   - soft DD  (double mass)
cc By definition we have MC=soft+hard+w(5)+w(1) (soft+hard is defined in WomTy) 
cc so fit-(soft+hard)=w(1)+w(5)=w(0)+w(5)+w(6)+w(7) and w(0)=w(8)+w(9)
cc w(9) is used as remaining difference (can be 0)
cc-----------------------------------------------------------------------
c#include "aaa.h"
c#include "par.h"
c#include "ems.h"
c#include "sem.h"
c      double precision plc,s
c      common/cems5/plc,s
c      double precision w(0:9),xp,xm,omGamk,om5Jk,wdh,wsd,wh
c     .                ,rp,zb,fom51d,zbd,rd
c
c      if(iregge.ne.0)w(5)=om5Jk(n,k,5)
c      if(iomega.ge.2)return
c      xp=xppr(n,k)
c      xm=xmpr(n,k)
c      rp=dble(r2had(iclpro)+r2had(icltar)
c     .       +slopom*log(max(1.,sngl(s*xpr(n,k)))))
c      zb=exp(-dble(bk(k))**2/(4.d0*.0389d0*rp))
c      wh=dble(q2kmin(0,n,k))  !hard
c      w(1)=max(0d0,omGamk(n,k,xp,xm,ntymin,ntymax)
c     .            -wh-dble(q2kmin(-1,n,k))-w(5)) !remove soft contrinution
cc fraction of high mass diff in soft or hard contribution
c      wdh=1d0
c      do i=2,4
c        w(i)=fom51d(s,xp,xm,zb,rp,zbd,rd,i+4)
c        if(i.eq.2)then !correction factor for hard diffraction
c          w(i)=dble(wdiff(iclpro))*w(i)
c        elseif(i.eq.3)then
c          w(i)=dble(wdiff(icltar))*w(i)
c        else
c          w(i)=dble(wdiff(iclpro)*wdiff(icltar))*w(i)
c        endif
c        wdh=wdh+w(i)
c      enddo
cc low mass diff contributions (remove hard from total for SD)
c      w(6)=max(0d0,om5Jk(n,k,-60)-w(2)*wh/wdh)
c      w(7)=max(0d0,om5Jk(n,k,-70)-w(3)*wh/wdh)
c      wsd=w(6)+w(7)
cc If w(1) is not large enough to cover SD contributions, rescale SD
c      if(wsd.gt.w(1).and.w(1).gt.0d0)then
c        w(6)=w(6)/wsd*w(1)
c        w(7)=w(7)/wsd*w(1)
c      endif
cc central diff = what is left
c      w(0)=max(0d0,w(1)-w(6)-w(7))
c      if(w(0).gt.0d0)then
c        w(8)=max(0d0,om5Jk(n,k,-8)-w(4)*wh/wdh)
cc If w(0) is not large enough to cover CD contribution, rescale CD
c        w(8)=min(w(0),w(8))
cc complete w(1) with w(9) to have correct sum always
c        w(9)=max(0d0,w(0)-w(8))
c      endif
cc      write(ifch,*)"WomTyDif",wh,w
c
c      return
c      end
c

c-----------------------------------------------------------------------
      double precision function Womegak(xp,xm,xprem,xmrem,n,k)   !---MC---
c-----------------------------------------------------------------------
c - sum(omGam(xp,xm))*(1-xp)*(1-xm) for group of cut enhanced
c diagram giving the same final state (without nuclear effect).
c xp,xm - fraction of the light cone momenta of the pomeron;
c k - pair index n - position index
c-----------------------------------------------------------------------
#include "aaa.h"
#include "ems.h"
      double precision xp,xm,xprem,xmrem!,dum

      Womegak=0.d0

      do i=ntymin,ntymax
        Womegak=Womegak+atildp(i,n,k)*xp**btildpp(i,n,k)
     &                               *xm**btildppp(i,n,k)
      enddo

c Alternative omGam*(xprem-xp)**alplea*(xmrem-xm)**alplea
      Womegak=Womegak*(xprem-xp)**alplea(iclpro)
     &               *(xmrem-xm)**alplea(icltar)
c      dum=xprem
c      dum=xmrem
c Alternative omGam*(xpr-xp)*(xmr-xm)
c      Womegak=Womegak*(xprem-xp)*(xmrem-xm)
c Alternative omGam only
c      dum=xprem
c      dum=xmrem

      return
      end

c----------------------------------------------------------------------
      double precision function omGam(xp,xm,bh)   !---test---
c-----------------------------------------------------------------------
c Cut diagram part for calculation of probability distribution
c xp,xm impulsion fraction of remnant
c bh - impact parameter between the pomeron ends;
c-----------------------------------------------------------------------
#include "aaa.h"
#include "ems.h"
      double precision om51,xp,xm,xh,yp,eps!,omYgam
      parameter (eps=1.d-20)

      omGam=0.d0
      if(xp.lt.eps.or.xm.lt.eps)return
      xh=xp*xm
      if(abs(xh).gt.1.d-10)then
        yp=0.5d0*log(xp/xm)
      else
        yp=0.d0
      endif

c call ipoOm5Tables if om51 is used with something else than (-1,-1)
      omGam=om51(xh,yp,bh,-1,-1)

c      omGam=2.d0*omGam

      return
      end

c----------------------------------------------------------------------
      double precision function omGamk(n,k,xp,xm,iq1,iq2)   !---MC---
c-----------------------------------------------------------------------
c Cut diagram part for calculation of probability distribution (for omega)
c xp,xm - light cone momentum fraction
c iq1- minimum type index
c iq2- maximum type index
c bh - impact parameter between the pomeron ends;
c-----------------------------------------------------------------------
#include "aaa.h"
#include "ems.h"
      double precision xp,xm,a,bp,bpp
      omGamk=0.d0
      imin=max(idxD0,iq1)
      imax=min(idxD1,iq2)
      if(iomega.ge.2)then
        if(.not.(imax.lt.ntynd.or.imin.gt.ntynd))then
          imin=max(imin,ntynd)
          imax=min(imax,ntynd)
        else
          return
        endif
      endif
c      if(atildp(0,n,k).le.0.d0)then
c        if(imin.eq.0.and.imax.eq.0)then
c          return
c        elseif(imin.eq.1)then
c          imin=0
c        endif
c      endif
      
      do i=imin,imax
      
        if(imin-imax.eq.0.or.i.ne.ntynd)then
          
          a=atildp(i,n,k)
          bp=btildpp(i,n,k)
          bpp=btildppp(i,n,k)
          
          omGamk=omGamk + a * xp**bp * xm**bpp
        
        endif
      
      enddo

      return
      end

c----------------------------------------------------------------------
      double precision function omGamint(bh)   !---test---
c-----------------------------------------------------------------------
c Integrated cut diagram part for calculation of probability distribution
c bh - impact parameter between the pomeron ends;
c-----------------------------------------------------------------------
#include "aaa.h"
      double precision omIgamint!,omYgamint

      omGamint=omIgamint(bh,0)
c      omGamint=2.d0*omIgamint(bh,0)

      return
      end





c----------------------------------------------------------------------
      block data dgdata
c----------------------------------------------------------------------
c constants for numerical integration (gaussian weights)
c----------------------------------------------------------------------
      double precision dgx1,dga1
      common /dga20/ dgx1(10),dga1(10)


      data dgx1/
     &   .765265211334973D-01,
     &   .227785851141645D+00,
     &   .373706088715420D+00,
     &   .510867001950827D+00,
     &   .636053680726515D+00,
     &   .746331906460151D+00,
     &   .839116971822219D+00,
     &   .912234428251326D+00,
     &   .963971927277914D+00,
     &   .993128599185095D+00/
      data dga1/
     &   .152753387130726D+00,
     &   .149172986472604D+00,
     &   .142096109318382D+00,
     &   .131688638449177D+00,
     &   .118194531961518D+00,
     &   .101930119817233D+00,
     &   .832767415767047D-01,
     &   .626720483341090D-01,
     &   .406014298003871D-01,
     &   .176140071391506D-01/

      end


c----------------------------------------------------------------------
      double precision function Phiexact(zzip,zzit,fj,xp,xm,s,b) !---test---
c----------------------------------------------------------------------
c    Exact expression of the Phi function for pp collision
c    zzip : additionnal component for Z (nuclear effect projectile side)
c    zzit : additionnal component for Z (nuclear effect target side)
c    fj   : overall factor for cross section (elastic or inelastic)
c    xp,xm: momentum fraction
c    s    : energy square
c    b    : impact parameter
c----------------------------------------------------------------------
#include "aaa.h"
#include "sem.h"
#include "ems.h"
      double precision al(idxD0:idxD1),betp(idxD0:idxD1)
     *,z,xIrst!,ffacto
      double precision zp(idxD0:idxD1),Phitmp,betpp(idxD0:idxD1)
     *,yp,ym,xm,xp
      double precision eps
      parameter(eps=1.d-20)
      dimension ipr(idxD0:idxD1),imax(idxD0:idxD1)

      if(idxD0.ne.0.or.idxD1.ne.3) stop "Problem in PhiExact"
      PhiExact=0.d0
      Phitmp=0.d0

      if(xp.gt.eps.and.xm.gt.eps.and.xp.le.1.d0+eps
     &   .and.xm.le.1.d0+eps)then


       do i=idxD0,idxD1
        imax(i)=0
        ipr(i)=0
        zp(i)=1.d0
        al(i)=0.d0
        betp(i)=0.d0
        betpp(i)=0.d0
       enddo

       do i=idxDmin(iomega),idxDmax(iomega)
        imax(i)=10+max(5,int(log10(s)))
        if(b.ge.1.)imax(i)=4+max(3,int(log10(sqrt(s))))
        imax(i)=min(30,imax(i))
       enddo
       do i=idxDmin(iomega),idxDmax(iomega)
         zpp=zzip
         zpt=zzit
         call Gfunpar(zpp,zpt,0.,0.,1,i,b,s,alpx,betx,betpx,epsp,epst
     .   ,epss,gamv)
         betp(i)=dble(betx)+1.d0
         betpp(i)=dble(betpx)+1.d0
         al(i)=dble(alpx*gamv)
       enddo

       do ipr0=0,imax(0)
          ipr(0)=ipr0
          zp(0)=1.d0
        if (ipr(0).ne.0) zp(0)=(-dble(fj)*al(0))**ipr(0)*facto(ipr(0))
        do ipr1=0,imax(1)
           ipr(1)=ipr1
           zp(1)=1.d0
        if (ipr(1).ne.0) zp(1)=(-dble(fj)*al(1))**ipr(1)*facto(ipr(1))
        do ipr2=0,imax(2)
           ipr(2)=ipr2
           zp(2)=1.d0
        if (ipr(2).ne.0) zp(2)=(-dble(fj)*al(2))**ipr(2)*facto(ipr(2))
        do ipr3=0,imax(3)
           ipr(3)=ipr3
           zp(3)=1.d0
        if (ipr(3).ne.0) zp(3)=(-dble(fj)*al(3))**ipr(3)*facto(ipr(3))
          yp=0.d0
          ym=0.d0
          z=1.d0
          isum=0
          do i=idxDmin(iomega),idxDmax(iomega)
            yp=yp+dble(ipr(i))*betp(i)
            ym=ym+dble(ipr(i))*betpp(i)
            isum=isum+ipr(i)
            z=z*zp(i)
          enddo

          z=z*xIrst(1,xp,yp,betp,ipr)
          z=z*xIrst(2,xm,ym,betpp,ipr)

          Phitmp=Phitmp+z

         enddo
         enddo
         enddo
        enddo



      endif

      PhiExact=Phitmp


      return
      end


c----------------------------------------------------------------------
      double precision function PhiExpoK(k,xp,xm)   !---MC---
c----------------------------------------------------------------------
c    Exponential expression of the Phi function for pp collision
c    for given k
c----------------------------------------------------------------------
#include "aaa.h"
#include "sem.h"
#include "ems.h"

      double precision xp,xm,Phitmp,Gt1
      double precision atildg,btildgp,btildgpp
      common/cgtilde/atildg(idxD0:idxD1,kollmx)
     *,btildgp(idxD0:idxD1,kollmx),btildgpp(idxD0:idxD1,kollmx)


      Phitmp=0.d0

      Phitmp=0.d0
      Gt1=0.d0
      do i=idxDmin(iomega),idxDmax(iomega)
       Gt1=Gt1+atildg(i,k)*xp**btildgp(i,k)*xm**btildgpp(i,k)
      enddo

      Phitmp=exp(-Gt1)

      PhiExpoK=Phitmp

      return
      end

c----------------------------------------------------------------------
      double precision function PhiExpoK2(k,xp,xm)   !---xs---
c----------------------------------------------------------------------
c    Exponential expression of the Phi function for pp collision
c    for given k without diffractive part
c----------------------------------------------------------------------
#include "aaa.h"
#include "sem.h"
#include "ems.h"

      double precision xp,xm,Phitmp,Gt1
      double precision atildg,btildgp,btildgpp
      common/cgtilde/atildg(idxD0:idxD1,kollmx)
     *,btildgp(idxD0:idxD1,kollmx),btildgpp(idxD0:idxD1,kollmx)


      Phitmp=0.d0

      Phitmp=0.d0
      Gt1=0.d0
      do i=idxD,idxD
       Gt1=Gt1+atildg(i,k)*xp**btildgp(i,k)*xm**btildgpp(i,k)
      enddo

      Phitmp=exp(-Gt1)

      PhiExpoK2=Phitmp

      return
      end

c----------------------------------------------------------------------
      double precision function Phiexpo(zzip,zzit,fj,xp,xm,s,b)   !---MC---
c----------------------------------------------------------------------
c    Exponential expression of the Phi function for pp collision
c    for given b
c input :
c    zzip : additionnal component for Z (nuclear effect projectile side)
c    zzit : additionnal component for Z (nuclear effect target side)
c    fj   : overall factor for cross section (elastic or inelastic)
c    xp,xm: momentum fraction
c    s    : energy square
c    b    : impact parameter
c----------------------------------------------------------------------
#include "aaa.h"
#include "sem.h"
#include "ems.h"
#include "par.h"

      parameter(nbkbin=100)
      common /kfitd/ xkappafit(nbkbin,nclegy,nclha,nclha),xkappa,bkbin
      double precision AlTi
      double precision BeTip,BeTipp
      double precision xp,xm,Phitmp,Gt1

      Gt1=0.d0
      do i=idxDmin(iomega),idxDmax(iomega)
        zp=zzip
        zt=zzit
        call Gfunpar(zp,zt,0.,0.,2,i,b,s,alpx,betx,betpx,epsp,epst,epss
     &              ,gamv)
        BeTip =dble(betx)
        BeTipp=dble(betpx)
        AlTi  =dble(alpx)
        Gt1=Gt1+AlTi*xp**BeTip*xm**BeTipp*dble(fj*xkappa**(fj-1.))

      enddo

      Phitmp=exp(-Gt1)

      PhiExpo=Phitmp
     &     *xp**dble(alplea(iclpro))
     &     *xm**dble(alplea(icltar))

c      if(xp.ge.1d0.and.xm.ge.1d0)print *,'phi',s,b,PhiExpo,imax

      return
      end

c----------------------------------------------------------------------
      double precision function PhiUnit(xp,xm)   !---test---
c----------------------------------------------------------------------
c    Exponential expression of the Phi function for pp collision
c    for given b
c----------------------------------------------------------------------
#include "aaa.h"
#include "sem.h"
#include "ems.h"
#include "par.h"

      double precision AlTi
      double precision BeTip,BeTipp
      double precision xp,xm,Phitmp,Gt1


      Gt1=0.d0
      do i=idxDmin(iomega),idxDmax(iomega)
        BeTip =betUni(i,2)
        BeTipp=betpUni(i,2)
        AlTi  =alpUni(i,2)
        Gt1=Gt1+AlTi*xp**BeTip*xm**BeTipp
c        write(ifch,*)'Phiunit',i,xp,xm,Gt1,AlTi,BeTip,BeTipp
      enddo

      Phitmp=exp(-Gt1)
c        write(ifch,*)'Phiunit',xp,xm,Phitmp

      PhiUnit=Phitmp
     &     *xp**dble(alplea(iclpro))
     &     *xm**dble(alplea(icltar))


      return
      end


cc----------------------------------------------------------------------
c      double precision function PhiUnit(xp,xm,s,b)   !---inu---
cc----------------------------------------------------------------------
c#include "aaa.h"
c      double precision xp,xm,PhiExpo,Znorm
c
c      PhiUnit=Phiexpo(0.,0.,1.,xp,xm,s,b)
c     &          /Znorm(s,b)
c
c      return
c      end
c
c     To have Znorm=1 with PhiExact, imax() should be large enough  
c----------------------------------------------------------------------
      double precision function Hrst(s,b,xp,xm)   !test
c----------------------------------------------------------------------
#include "aaa.h"
#include "ems.h"
#include "sem.h"
#include "par.h"
      parameter(idxD2=idxD1)
      double precision GbetUni,GbetpUni,HbetUni,HbetpUni,HalpUni
      common/DGamUni/GbetUni(  idxD0:idxD2),HbetUni(  idxD0:idxD2),
     &               GbetpUni(idxD0:idxD2),HbetpUni(idxD0:idxD2),
     &               HalpUni(idxD0:idxD2)
      double precision al(idxD0:idxD2),betp(idxD0:idxD2)
     *,z,xJrst!,ffacto
      double precision zp(idxD0:idxD2),Htmp,betpp(idxD0:idxD2)
     *,yp,ym,xp,xm
      dimension ipr(idxD0:idxD2),imax(idxD0:idxD2)

      if(idxD0.ne.0.or.idxD1.ne.3) stop "Problem in Hrst"
      Hrst=0.d0
      do i=idxD0,idxD2
        imax(i)=0
        ipr(i)=0
        zp(i)=1.d0
        al(i)=0.d0
      enddo

      if(xp.ge.0.d0.and.xm.ge.0.d0.and.xp.lt.1.d0.and.xm.le.1.d0)then

      imax0=idxD1
      if(iomega.ge.2)imax0=idxD
      imax1=idxD2
      if(iomega.ge.2)imax1=idxD

      do i=idxDmin(iomega),imax1
        imax(i)=6+int(log10(100.*s)/3.)
c        if(i.gt.idxD)imax(i)=imax(i)*2
        if(b.ge.2.)imax(i)=max(3,imax(i)/2)
        imax(i)=min(30,imax(i))
      enddo

      Htmp=0.d0
      do i=idxDmin(iomega),imax1
        betp(i)=HbetUni(i)
        betpp(i)=HbetpUni(i)
        al(i)=HalpUni(i)
      enddo

      do ipr0=0,imax(0)
c       write(ifch,*)'Hrst ipr0,xp,xm :',ipr0,xp,xm,Htmp
        ipr(0)=ipr0
        zp(0)=1.d0
        if (ipr(0).ne.0) zp(0)=al(0)**ipr(0)*facto(ipr(0))
        do ipr1=0,imax(1)
          ipr(1)=ipr1
          zp(1)=1.d0
          if (ipr(1).ne.0) zp(1)=al(1)**ipr(1)*facto(ipr(1))
          do ipr2=0,imax(2)
            ipr(2)=ipr2
            zp(2)=1.d0
            if (ipr(2).ne.0) zp(2)=al(2)**ipr(2)*facto(ipr(2))
           do ipr3=0,imax(3)
            ipr(3)=ipr3
            zp(3)=1.d0
            if (ipr(3).ne.0) zp(3)=al(3)**ipr(3)*facto(ipr(3))
            if (ipr(0)+ipr(1)+ipr(2)+ipr(3).ne.0) then
c            if (ipr(0)+ipr(1)+ipr(2).ne.0) then
             yp=0.d0
             ym=0.d0
             z=1.d0
             do i=idxDmin(iomega),imax1
               yp=yp+dble(ipr(i))*betp(i)
               ym=ym+dble(ipr(i))*betpp(i)
               z=z*zp(i)
             enddo
             z=z*xJrst(xp,yp,GbetUni,ipr)
             z=z*xJrst(xm,ym,GbetpUni,ipr)
             Htmp=Htmp+z
            endif
           enddo
          enddo
        enddo
      enddo
      Hrst=Htmp
          
      endif


      return
      end

c----------------------------------------------------------------------
      double precision function HrstI(s,b,xp,xm)   !test
c----------------------------------------------------------------------
#include "aaa.h"
#include "ems.h"
#include "sem.h"
#include "par.h"
      parameter(idxD2=idxD1)
      double precision GbetUni,GbetpUni,HbetUni,HbetpUni,HalpUni
      common/DGamUni/GbetUni(  idxD0:idxD2),HbetUni(  idxD0:idxD2),
     &               GbetpUni(idxD0:idxD2),HbetpUni(idxD0:idxD2),
     &               HalpUni(idxD0:idxD2)
      double precision al(idxD0:idxD2),betp(idxD0:idxD2)
     *,z,xJrstI!,ffacto
      double precision zp(idxD0:idxD2),Htmp,betpp(idxD0:idxD2)
     *,yp,ym,xp,xm
      dimension ipr(idxD0:idxD2),imax(idxD0:idxD2)

      if(idxD0.ne.0.or.idxD1.ne.3) stop "Problem in HrstI"
      HrstI=0.d0
      do i=idxD0,idxD2
        imax(i)=0
        ipr(i)=0
        zp(i)=1.d0
        al(i)=0.d0
      enddo

      Htmp=0.d0

      if(xp.ge.0.d0.and.xm.ge.0.d0.and.xp.lt.1.d0.and.xm.le.1.d0)then


      imax0=idxD1
      if(iomega.ge.2)imax0=idxD
      imax1=idxD2
      if(iomega.ge.2)imax1=idxD

      do i=idxDmin(iomega),imax1
        imax(i)=max(3,int(log10(s)/2.))
c        if(i.gt.idxD)imax(i)=imax(i)*2
        if(b.ge.2.)imax(i)=max(2,imax(i)/2)
        imax(i)=min(30,imax(i))
      enddo

      do i=idxDmin(iomega),imax1
        betp(i)=HbetUni(i)
        betpp(i)=HbetpUni(i)
        al(i)=HalpUni(i)
      enddo
      do ipr0=0,imax(0)
        ipr(0)=ipr0
        zp(0)=1.d0
        if (ipr(0).ne.0) zp(0)=al(0)**ipr(0)*facto(ipr(0))
        do ipr1=0,imax(1)
          ipr(1)=ipr1
          zp(1)=1.d0
          if (ipr(1).ne.0) zp(1)=al(1)**ipr(1)*facto(ipr(1))
          do ipr2=0,imax(2)
            ipr(2)=ipr2
            zp(2)=1.d0
            if (ipr(2).ne.0) zp(2)=al(2)**ipr(2)*facto(ipr(2))
           do ipr3=0,imax(3)
            ipr(3)=ipr3
            zp(3)=1.d0
            if (ipr(3).ne.0) zp(3)=al(3)**ipr(3)*facto(ipr(3))
            if (ipr(0)+ipr(1)+ipr(2)+ipr(3).ne.0) then
c            if (ipr(0)+ipr(1)+ipr(2).ne.0) then
             yp=0.d0
             ym=0.d0
             z=1.d0
             do i=idxDmin(iomega),imax1
               yp=yp+dble(ipr(i))*betp(i)
               ym=ym+dble(ipr(i))*betpp(i)
               z=z*zp(i)
             enddo
             z=z*xJrstI(xp,yp,GbetUni,ipr)
             z=z*xJrstI(xm,ym,GbetpUni,ipr)
             Htmp=Htmp+z
            endif
           enddo
          enddo
        enddo
      enddo

      endif

      HrstI=Htmp

      return
      end



cc----------------------------------------------------------------------
c      double precision function HrstI(s,xp,xm)   !---inu---
cc----------------------------------------------------------------------
c#include "aaa.h"
c#include "ems.h"
c#include "sem.h"
c#include "par.h"
c      double precision al(idxD0:idxD1),betp(idxD0:idxD1)
c     *,z,xJrstI!,ffacto
c      double precision zp(idxD0:idxD1),Htmp,betpp(idxD0:idxD1)
c     *,yp,ym,xp,xm
c      dimension ipr(idxD0:idxD1),imax(idxD0:idxD1)
c
c      if(idxD0.ne.0.or.idxD1.ne.3) stop "Problem in HrstI"
c      HrstI=0.d0
c      do i=idxD0,idxD1
c        imax(i)=0
c        ipr(i)=0
c        zp(i)=1.d0
c        al(i)=0.d0
c      enddo
c
c      if(xp.ge.0.d0.and.xm.ge.0.d0.and.xp.lt.1.d0.and.xm.lt.1.d0)then
c
c      HrstI=0.d0
c
c
c      do i=idxDmin(iomega),idxDmax(iomega)
c        imax(i)=max(5,int(log10(s)))
cc        if(i.ge.2)imax(i)=imax(i)*2
c        imax(i)=min(30,imax(i))
c      enddo
c      Htmp=0.d0
c        do i=idxDmin(iomega),idxDmax(iomega)
c          betp(i)=betUni(i,1)+1.d0
c          betpp(i)=betpUni(i,1)+1.d0
c          al(i)=alpUni(i,1)
c        enddo
c
c        do ipr0=0,imax(0)
c           ipr(0)=ipr0
c           zp(0)=1.d0
c          if (ipr(0).ne.0) zp(0)=al(0)**ipr(0)*facto(ipr(0))
c         do ipr1=0,imax(1)
c            ipr(1)=ipr1
c            zp(1)=1.d0
c          if (ipr(1).ne.0) zp(1)=al(1)**ipr(1)*facto(ipr(1))
c         do ipr2=0,imax(2)
c            ipr(2)=ipr2
c            zp(2)=1.d0
c          if (ipr(2).ne.0) zp(2)=al(2)**ipr(2)*facto(ipr(2))
c         do ipr3=0,imax(3)
c            ipr(3)=ipr3
c            zp(3)=1.d0
c          if (ipr(3).ne.0) zp(3)=al(3)**ipr(3)*facto(ipr(3))
c             if (ipr(0)+ipr(1)+ipr(2)+ipr(3).ne.0) then
c             yp=0.d0
c             ym=0.d0
c             z=1.d0
c             do i=idxDmin(iomega),idxDmax(iomega)
c               yp=yp+dble(ipr(i))*betp(i)
c               ym=ym+dble(ipr(i))*betpp(i)
c               z=z*zp(i)
c             enddo
c             z=z*xJrstI(xp,yp,betp,ipr)
c             z=z*xJrstI(xm,ym,betpp,ipr)
c             Htmp=Htmp+z
c           endif
c          enddo
c         enddo
c       enddo
c
c       HrstI=Htmp
c
c      endif
c
c      return
c      end
c

cc----------------------------------------------------------------------
c        double precision function ffacto(n)   !---test---
cc----------------------------------------------------------------------
c
c        ffacto=1.D0
c        do i=1,n
c          ffacto=ffacto*dble(i)
c        enddo
c        return
c        end
c

c----------------------------------------------------------------------
      double precision function xIrst(id,x,y,bet,ipr)   !---test---
c----------------------------------------------------------------------
#include "aaa.h"
      double precision y,gammag,utgam2,x,bet(idxD0:idxD1)
      dimension ipr(idxD0:idxD1)

      if(id.eq.1)iclrem=iclpro
      if(id.eq.2)iclrem=icltar
      if(y.le.160.)then
       xIrst=gammag(iclrem,y)*x**dble(alplea(iclrem))
      else
       xIrst=0
      endif
      if(xIrst.gt.0.d0)then
        do i=idxDmin(iomega),idxDmax(iomega)
          if(ipr(i).ne.0.and.bet(i).gt.1.d-10)
     &         xIrst=xIrst*utgam2(bet(i))**dble(ipr(i))
        enddo
        if (abs(y).gt.1.d-10) xIrst=xIrst*x**y
      endif
      return
      end


c----------------------------------------------------------------------
      double precision function xJrst(x,y,Gbeta,ipr)   !---inu---
c----------------------------------------------------------------------
#include "aaa.h"
      parameter(idxD2=idxD1)
      double precision y,utgam2,x,Gbeta(idxD0:idxD2),eps,gam
      dimension ipr(idxD0:idxD2)

      eps=1.d-10


      gam=utgam2(y)

      if(gam.lt.1.d99)then

      if ((x-1.d0).gt.eps.or.(y-1.d0).gt.eps) then
                        xJrst=(1.d0-x)**(y-1.d0)/gam
                        do i=idxDmin(iomega),idxDmax(iomega)
      if (ipr(i).ne.0)   xJrst=xJrst*Gbeta(i)**dble(ipr(i))
            
                        enddo
          else
c            write (*,*) 'Warning in xJrst, infinite value !'
                        xJrst=(1.d0-x+eps)**(y-1.d0)/gam
                        do i=idxDmin(iomega),idxDmax(iomega)
      if (ipr(i).ne.0)   xJrst=xJrst*Gbeta(i)**dble(ipr(i))
                        enddo
      endif
      else
        xJrst=0.d0
      endif

      return
      end


c----------------------------------------------------------------------
      double precision function xJrstI(x,y,Gbeta,ipr)   !---inu---
c----------------------------------------------------------------------
c Function used for the integration of H*Phi. We do the changement of
c variable (1-x)=z**alpha. The power alpha can be change if necessary.
c----------------------------------------------------------------------
#include "aaa.h"
      parameter(idxD2=idxD1)
      double precision y,utgam2,x,Gbeta(idxD0:idxD2),alpha,w,gam
      dimension ipr(idxD0:idxD2)

      alpha=4.d0
      w=alpha*(y-1.d0)+alpha-1.d0
      imax=idxD2
      if(iomega.ge.2)imax=idxD

      gam=utgam2(y)

      if(gam.lt.1.d99)then

      if (w.ge.0)then

                        xJrstI=alpha*x**w/gam
                        do i=idxDmin(iomega),imax
      if (ipr(i).ne.0)   xJrstI=xJrstI*Gbeta(i)**dble(ipr(i))
                        enddo

         else
           write(*,*) 'x,y,bet,ipr,w',x,y,Gbeta,ipr,w
          stop 'Error in xJrstI in epos-omg, integration not possible'
       endif

      else
        xJrstI=0.d0
      endif


      return
      end

c----------------------------------------------------------------------
      double precision function HPhiInt(s,b)   !---inu---
c----------------------------------------------------------------------
c  Set integrated over xp and xm (x and y) H(x,y)*Phi(x,y) for a
c  given b by gauss method
c  PhiExact should be used to check that Znorm=1
c----------------------------------------------------------------------
#include "aaa.h"
      parameter(idxD2=idxD1)
      double precision GbetUni,GbetpUni,HbetUni,HbetpUni,HalpUni
      common/DGamUni/GbetUni(  idxD0:idxD2),HbetUni(  idxD0:idxD2),
     &               GbetpUni(idxD0:idxD2),HbetpUni(idxD0:idxD2),
     &               HalpUni(idxD0:idxD2)
      double precision xhm,x,y,yhm,w,Hrst,utgam2,PhiUnit!,PhiExact
c      double precision zp2,zm2,HrstI,eps
c      common /ar3/  x1(7),a1(7)
      common /ar9/    x9(3),a9(3)

      eps=0d0 !1.d-5

      imax0=idxD1
      imax1=idxD2
      if(iomega.ge.2)then
        imax0=idxD
        imax1=idxD
      endif
      do i=idxDmin(iomega),imax0
        HbetUni(i)=betUni(i,1)+1.d0
        HbetpUni(i)=betpUni(i,1)+1.d0
        GbetUni(i)=utgam2(HbetUni(i))
        GbetpUni(i)=utgam2(HbetpUni(i))
        HalpUni(i)=alpUni(i,1)
      enddo

      w=0.d0
      xhm=.5d0*(1d0-eps)
      yhm=.5d0*(1d0-eps)
      do m=1,2
        do i=1,3
c        do i=1,7
          x=xhm+dble((2*m-3)*x9(i))*xhm
c          write(ifmt,*)'HPhiInt, xp int :',x
          do n=1,2
            do j=1,3
c            do j=1,7
              y=yhm+dble((2*n-3)*x9(j))*yhm
             w=w+dble(a9(i)*a9(j))*Hrst(s,b,x,y)
     &          *PhiUnit(x,y)
c     &          *Phiexact(0.,0.,1.,x,y,s,b)
            enddo
          enddo
        enddo
      enddo

      HPhiInt=w*xhm*yhm


c      w=0.d0
c      xhm=.5d0*eps
c      yhm=.5d0*eps
c      do m=1,2
c        do i=1,7
c          x=1d0-eps+xhm+dble((2*m-3)*x1(i))*xhm
c          do n=1,2
c            do j=1,7
c              y=1d0-epsyhm+dble((2*n-3)*x1(j))*yhm
c              zp2=1.d0-x**4
c              zm2=1.d0-y**4
c              w=w+dble(a1(i)*a1(j))*HrstI(s,x,y)
cc     &             *PhiUnit(zp2,zm2)
c     &             *Phiexact(0.,0.,1.,zp2,zm2,s,b)
c            enddo
c          enddo
c        enddo
c      enddo
c
c      HPhiInt=HPhiInt+w*xhm*yhm

      return
      end



c----------------------------------------------------------------------
      subroutine Kfit(iiclegy)
c----------------------------------------------------------------------
c q2pmin defined in mkParamTable for tabulation (b vary with fixed Q2s)
c----------------------------------------------------------------------
#include "aaa.h"
#include "sem.h"
      double precision Znorm
      parameter(nbkbin=100)
      common /kfitd/ xkappafit(nbkbin,nclegy,nclha,nclha),xkappa,bkbin
      parameter (nmax=30)
      logical lnoch

      if(iiclegy.le.0.or.iiclegy.gt.iclegy2)then

        if(iiclegy.eq.0.or.abs(iiclegy).gt.iclegy2)then
          do iiitar=1,nclha
            do iiipro=1,nclha
              do iiiegy=1,nclegy
                do iiibk=1,nbkbin
                  xkappafit(iiibk,iiiegy,iiipro,iiitar)=1.
                enddo
              enddo
            enddo
          enddo          
        else
          iiiegy=abs(iiclegy)
          do iiibk=1,nbkbin
            xkappafit(iiibk,iiiegy,iclpro,icltar)=1.
          enddo
        endif

      else

      s=engy*engy
      if(isetcs.le.1)then
c      if(isetcs.le.2)then       !??????????to speed up
        eps=0.01
      else
        eps=0.001
      endif

      !write(ifmt,*)"Fit xkappa ..."
      if(ish.ge.5)then
        write(ifch,*)"Kfit s,bkbin,iclegy,ipro,itar"
     *       ,s,bkbin,iiclegy,iclpro,icltar
      endif


      b=0.
      xkfs=0.
      deltas=0.
      if(isetcs.le.1.or.iiclegy.eq.iclegy2)then
        xkf=1.
      else
        xkf=xkappafit(1,iiclegy+1,iclpro,icltar)
      endif
      delta=0.


      do 5 ib=1,nbkbin-1
        if(ib.gt.1.and.ish.ge.3)write(ifch,*)"    End",b,delta,xkf
        b=float(ib-1)*bkbin
        xkappafit(ib,iiclegy,iclpro,icltar)=1.
        if(b.gt.3.+0.05*log(s))then
          xkf=1.
          goto 5
        endif
        delta=1.-sngl(Znorm(s,b))
        if(delta.le.0d0)then            !to comment to correct z<1
          if(xkf.gt.1.)then
            xkappafit(ib,iiclegy,iclpro,icltar)=xkf
            delta=1.-sngl(Znorm(s,b))
          else
            xkf=1.
          endif
        else       !to comment to correct z<1
          xkf=1.       !to comment to correct z<1
          goto 5         !accept Znorm<1 if Znorm calculated with PhiExact is < 1 !to comment to correct z<1
        endif       !to comment to correct z<1
        if(abs(delta).lt.eps)then
            xkfs=xkf-delta
            deltas=delta
          goto 5
        elseif(ib.le.nbkbin-1)then

          if(delta.gt.0.d0)then
            xkf0=1.
            xkf1=xkf
            delta0=delta
            xkf2=xkf-delta0
            xkappafit(ib,iiclegy,iclpro,icltar)=xkf2
            delta1=1.-sngl(Znorm(s,b))
            if(delta1.le.0.d0)then
              xkf0=xkf2
              xkf1=xkf
              delta=delta1
              xkf=xkf0
            else
              xkf1=max(delta0,xkf2)
              xkf0=max(0.,xkf1-5.*delta1)
              xkf=xkf1
            endif
          else
            xkf0=xkf
            xkf1=1.-delta
            xkf2=xkf
            delta1=delta
          endif

          if(ib.eq.1)then
            deltas=delta
            xkfs=max(0.00001,1.-delta)
          endif

          if(delta.le.deltas)xkf=xkfs
          if(ish.ge.3)write(ifch,*)"    Start",ib,b,delta,xkf,xkf0,xkf1
          if(xkf.eq.xkf2)delta=delta1

          n=0
          delta0=delta
          lnoch=.true.
 10       continue
          n=n+1
          if(n.le.nmax.and.abs(xkf1-xkf0).gt.1.e-3)then
            if(abs(xkf-xkf2).gt.1e-6.or.abs(delta).gt.abs(deltas))then
              xkappafit(ib,iiclegy,iclpro,icltar)=xkf
              delta=1.-sngl(Znorm(s,b))
            endif
            if(ish.ge.5)write(ifch,*)"    step",ib,n,delta,xkf,delta0
            if(delta*delta0.ge.0.)then
              if(lnoch.and.abs(delta).gt.abs(delta0))goto 5
            else
              lnoch=.false.
            endif
            if(abs(delta).gt.eps)then
              if(delta.gt.0.)then
                xkf1=xkf
                if(xkf.gt.100.)then
                  xkf=sqrt(xkf1)
                else
                  xkf=(xkf1+xkf0)*0.5
                endif
                delta0=delta
              else
                xkf0=xkf
                xkf=(xkf1+xkf0)*0.5
                delta0=delta
              endif
              goto 10
            endif
          else
            if(ish.ge.2.and.delta.gt.0.1)
     *      write(ifmt,*)"Warning in Kfit, nmax reached : xkappafit=1."
     *                  ,delta
            xkappafit(ib,iiclegy,iclpro,icltar)=xkf
          endif
        endif

 5    continue

      if(ish.ge.3)write(ifch,*)"    End",b,delta,xkf
      if(xkf.gt.1.+eps)write(ifmt,*)
     *     "Warning in Kfit, xkappafit not yet 1"
      xkappafit(nbkbin,iiclegy,iclpro,icltar)=1.
      
      endif

      return
      end

c----------------------------------------------------------------------
      double precision function Znorm(s,b)   !---inu---
c----------------------------------------------------------------------
c Normalization function : should be 1 if Phiexact is used here and in
c HPhiInt.
c----------------------------------------------------------------------
#include "aaa.h"
      common /kwrite/ xkapZ
      double precision HPhiInt,PhiUnit!,PhiExact

      iomegasave=iomega
      iomega=2
      call DefXminDf(dble(s))

c      write(ifmt,*)'Z calculation for (s,b) :',s,b
      do i=idxDmin(iomega),idxDmax(iomega)
        zp=0.
        zt=0.
        call Gfunpar(zp,zt,0.,0.,1,i,b,s,alpx,betx,betpx,epsp,epst,epss
     &              ,gamv)
        zp=0.
        zt=0.
        call Gfunpar(zp,zt,0.,0.,2,i,b,s,alpx,betx,betpx,epsp,epst,epss
     &              ,gamv)
      enddo
      Znorm=HPhiInt(s,b)
c      write(ifch,*)'int',Znorm,' phi',Phiexact(0.,0.,1.,1.d0,1.d0,s,b)
      Znorm=Znorm
     &       +PhiUnit(1.d0,1.d0)
c     &       +Phiexact(0.,0.,1.,1.d0,1.d0,s,b)

      !write(ifmt,*)'Z=',Znorm,xkapZ,b
      iomega=iomegasave
      return
      end


c------------------------------------------------------------
      double precision function gammag(iclrem,x)   !---test---
c----------------------------------------------------------------------
#include "aaa.h"
#include "sem.h"
      double precision x,utgam2

      gammag=utgam2(dble(alplea(iclrem))+1.D0)
     &       /utgam2(dble(alplea(iclrem))+1.D0+x)

      return
      end


c
cc----------------------------------------------------------------------
c      double precision function GCorr(xp,xm,b,k,rho) !---MC---
cc----------------------------------------------------------------------
cc Fonction to transform om1->om1_Sat (Z saturating for b->0)
cc----------------------------------------------------------------------
c#include "aaa.h"
c#include "ems.h"
c#include "par.h"
c      double precision xp,xm,xdum
c
c      xdum=xp
c      xdum=xm
c      if(k.eq.0)then
c        if(abs(epspUni(1)-1e-5).lt.0.)then
c          spp=engy*engy
c          zp=0.
c          zt=0.
c          call Gfunpar(zp,zt,rho,1,1,b,spp,alp,bet,betp,epsp,epst
c     &           ,epss,gamv)
c        else
c          epsp=epspUni(1)
c          epst=epstUni(1)
c          epss=epssUni(1)
c        endif
c      else
c        epsp=epsilongp(1,k)
c        epst=epsilongt(1,k)
c        epss=epsilongs(1,k)
c      endif
c
c      fegy=fegypp+rho
c      bc=bcutz(fegy,b2xscr)
c      b2b=max(min(0.,bc),abs(b)-abs(bc))
c      b2a=max(0.,b2b)**2
c      b2b=b2b**2
c      fmn=0.
c      if(abs(b).lt.abs(bc))fmn=epscrb
c      corr=exp(-b2a/b2xscr)/exp(-fmn*b2b/b2xscr)-1.
c      betc=epsp*corr
c      betpc=epst*corr
c      GCorr=1.d0 !xp**betc*xm**betpc
c      if(epss.lt.0.)GCorr=GCorr/xminDf**(corr*epss)
cc      print *,b,xp,xm,corr,epsp,epst,epss,betc,betpc,GCorr
c
c      return
c      end
c
cc----------------------------------------------------------------------
c      function SatCorr(xp,xm,b,k,rho) !---MC---
cc----------------------------------------------------------------------
cc correction for saturation from MC to Q2s (decreasing Z to have more
cc Pom but increasing Q2s)
cc----------------------------------------------------------------------
c#include "aaa.h"
c#include "ems.h"
c#include "par.h"
c      double precision xp,xm,GCorr
c
c      SatCorr=sngl(1d0/GCorr(xp,xm,b,k,rho))   !G_MC->Geff*SatCorr=G(Q2s)
cc      print *,'satcorr',b,xp,xm,rho,SatCorr
c
c
c      return
c      end

c----------------------------------------------------------------------
      function SigGFF()   !---MC---
c----------------------------------------------------------------------
c integral d2b of number of Pom including screening
c----------------------------------------------------------------------
#include "aaa.h"
#include "par.h"
      double precision om1intbc!,om1intbh,PhiExpo!,omIgamint
      common/geom/rmproj,rmtarg,bmax,bkmx
      common /ar3/  x1(7),a1(7)

      spp=engy**2
      q2pmin(1)=q2nmin
      q2pmin(2)=q2nmin
      call ipoOm5Tables(1)
      bmid=bkmx/2.
      SigGFF=0.
      do i=1,7
        do m=1,2
          bb=bmid*(1.+(2.*m-3)*x1(i))
          SigGFF=SigGFF+bb*a1(i)*sngl(om1intbc(bb))
        enddo
      enddo

      SigGFF=SigGFF*2.*pi*bmid*10d0   !*10 because of fm**2->mb


      return
      end


cc----------------------------------------------------------------------
c      subroutine NhpomQ2s(k,SigGFFk,SigGFFc)   !---MC---
cc----------------------------------------------------------------------
cc integral d2b of number of hard Pom including screening for a given pair
cc SigGFFk : mean number of hard Pom
cc SigGFFc : mean number of hard Pom taking into account SatCorr
cc----------------------------------------------------------------------
c#include "aaa.h"
c#include "ems.h"
c#include "sem.h"
c#include "par.h"
c      common/cems5/plc,s
c      double precision s,plc
c      double precision cint,Sig,cintx,PhiUnit,omdiff,om5J,wwc,ww0,cintc!,utgam2
c     .,xpmin,xpmax,xmmin,xmmax,xp,xm,w,ww,www,om1f,SigGFFk,SigGFFc
c     .,xpc,xmc
c      common/geom/rmproj,rmtarg,bmax,bkmx
c      common /ar3/  x1(7),a1(7)
c
c      Sig=0.d0
c      SigGFFk=0.d0
c      SigGFFc=0.d0
c
c      return     !not used for the moment (and time consumming if used)
c
c      Fp=dble(ucfpro)           !gamma(1+alplea)
c      Fm=dble(ucftar)
c
c      spp=sngl(s)
c      rho=zparpnx(2,k)+zpartnx(2,k)
cc xp and xm integration with Pom
c      xpmin=max(dble(4.*q2nmin)/s,8.d-4*log(s))
cc      xpmin=max(1.d0/s,5.d-5)
c      xpmax=1d0
c      xmmin=max(dble(4.*q2nmin)/s,4.d-4*log(s)) !s0min/s   !use xmin=1e-3*log(s) to get proper result with the (faster) single precision integration. Comparison to MC should give ~1 at all energies  
cc      xpmin=xmmin
c      xmmax=1d0
c      bmid=bkmx/2.
c      q2pmin(1)=q2nmin          !start at minimum Q2s
c      q2pmin(2)=q2pmin(1)
c      call ipoOm5Tables(1)
c      binscal=1d0!/dble(xnpomcol)!*dble(sigcut/sigine)  !use <Npom> instead of <Nhard> here
c      zp=0.
c      zt=0.
c      do i=1,7
c        do m=1,2
c          bb=bmid*(1.+(2.*m-3)*x1(i))
c
c          cint=0.d0
c          cintc=0.d0
c          do ii=idxDmin(iomega),idxDmax(iomega)
c      
c            if(ii.eq.1)then
c              zp=zparpnx(1,k)
c              zt=zpartnx(1,k)
c            else
c              zp=zparpnx(0,k)
c              zt=zpartnx(0,k)
c            endif
c
c            call Gfunpar(zp,zt,rho,1,ii,bb,spp,alp,bet,betp,epsp,epst
c     &           ,epss,gamv)
cc            gamom=dble(alp*gamv)
cc            deltap=dble(bet)
cc            deltam=dble(betp)
cc            cint=cint+gamom*utgam2(deltap+1.d0)*utgam2(deltam+1.d0)
cc     &            /utgam2(2.d0+deltap+dble(alplea(iclpro)))
cc     &            /utgam2(2.d0+deltam+dble(alplea(icltar)))
c            call Gfunpar(zp,zt,rho,2,ii,bb,spp,alp,bet,betp,epsp,epst
c     &            ,epss,gamv)
c          enddo
cc          cint=cint*Fp*Fm
cc cross section for 1 diffractive Pomeron exchange (multiple scattering negligible)
c          cintx=0d0
c          omdiff=0d0
c          if(xmmax.gt.xmmin.and.xpmax.gt.xpmin)then
c
c            do n=1,2
c              do jj=1,7
c                xp=xpmin*(xpmax/xpmin)**dble(.5+x1(jj)*(n-1.5))
c                xpc=xmmin*(xmmax/xmmin)**dble(.5+x1(jj)*(n-1.5))
c                w=0d0
c                www=0d0
c                wwc=0d0
c                do l=1,2
c                  do j=1,7
c                    xm=xpmin*(xpmax/xpmin)**dble(.5+x1(j)*(l-1.5))
c                    xmc=xmmin*(xmmax/xmmin)**dble(.5+x1(j)*(l-1.5))
c                    ww0=dble(a1(j))
c                    w=w+ww0*om5J(xp,xm,bb,-50,-50)*xp*xm
c     .               *PhiUnit(1d0-xp,1d0-xm)
c                    ww=om1f(xp,xm,0,10)       !fit
c                    omt=0d0
c                    do ii=1,4
c                      omt=omt+binscal*om5J(xp,xm,bb,ii,ii)
c                      if(omt.ge.ww)then
c                        omt=ww
c                        goto 10
c                      endif
c                    enddo
cc                    omt=omt+min(ww-omt,om5J(xp,xm,bb,11,11)) !can take a long time (not tabulated) while not very important
c 10                 continue
cc                    print *,xp,xm,omt,www,ww
c                    wwc=wwc+ww0*dble(SatCorr(xp,xm,bb,0,rho))*xpc*xmc
c                    www=www+ww0*omt*xp*xm
c     .               *(1d0-xp)**alplea(iclpro)*(1d0-xm)**alplea(icltar)
c
c                enddo
c                enddo
c                cintx=cintx+dble(a1(jj))*0.5d0*log(xpmax/xpmin)*w
c                cint=cint+dble(a1(jj))*0.5d0*log(xpmax/xpmin)*www
c                cintc=cintc+dble(a1(jj))*0.5d0*log(xmmax/xmmin)*wwc
c              enddo
c            enddo
c            omdiff=cintx*0.5d0*log(xpmax/xpmin)
c            cint=cint*0.5d0*log(xpmax/xpmin)
c            cintc=cintc*0.5d0*log(xmmax/xmmin)
c          endif
c          Sig=Sig+dble(bb*a1(i))*max(0d0,1d0-PhiUnit(1d0,1d0)-omdiff)  !cut cross-section
cc         Sig=Sig+dble(bb*a1(i))*max(0d0,1d0-PhiUnit(1d0,1d0))  !inelastic cross-section
c          SigGFFk=SigGFFk+dble(bb*a1(i))*max(0d0,cint) !cross-section of Pom production
c          SigGFFc=SigGFFc+dble(bb*a1(i))*max(0d0,cintc)
c     .                   *max(0d0,1d0-PhiUnit(1d0,1d0)-omdiff) !cross-section of Pom production
c        enddo
c      enddo
c
cc      print *,'siggffk',(SigGFFk*2.*pi*bmid*10d0),(Sig*2.*pi*bmid*10d0)
cc     .     ,zp,zt,rho,k
c      SigGFFc=SigGFFc/max(0.0001d0,Sig)       !normalization cancel
c      SigGFFk=SigGFFk/max(0.0001d0,Sig)       !normalization cancel
c
c      return
c      end





c####################################################################################
c#############   former chk #########################################################
c####################################################################################



c----------------------------------------------------------------------
      double precision function xmpi(x,a)   !---MC---
c-----------------------------------------------------------------------
c     Transformation function to fit PomInc distribution with MPI
c     iqq=0, use x directly
c     iqq=1, increase x by col first
c     a and d define "coef" which fixes the slop toward x->1 (if iqq=0 or
c     1/(1+b*log(col)) if iqq=1). If coef is large, looks like a threshold
c     effect, if small, we get a smoother transition from max to min in
c     Pominc.
c     c will fix the minimum in Pominc (x closer and closer to 1 when c
c     increase (needs large value)
c-----------------------------------------------------------------------
      implicit none
#include "par.h"
      double precision x,x1,x2
      real a(4),coef
      x1=x
      x1=(x-xminDf)*a(1)
      x1=xminDf+x1
      x1=min(1.d0,x1)
      coef=abs(a(2))/(1.+a(3))
      if(a(2).lt.0.)coef=1./coef
      x2=x1
      xmpi=x1**(1.d0/(1.d0+dble(a(4))*x2**coef))

      end

c----------------------------------------------------------------------
      double precision function PomInc(fct,x,bb,col1,col2,col3,iqq)   !---MC---
c-----------------------------------------------------------------------
c     Transformation function to fit PomInc distribution with MPI
c     col1, col2 counts Poms connected to other nucleons (nuclear effect) 
c     iqq=1, for x
c     iqq=2, for x+ or x-
c     iqq>0 bb dependence
c     iqq<0 bb independent      
c-----------------------------------------------------------------------
      implicit none
#include "aaa.h"
      double precision fct,x,x1,x2,xmpi
      real bb,col1,col2,col3,a(5),b(4),c(4),col!,frac
      integer iqq
      external fct
      

      if(abs(iqq).eq.1)then
        a(1)=1.
        a(2)=2.2
        a(3)=1.77 !2.
        a(4)=3. !7.
        a(5)=1.
        col=dsatur*log(max(1.,col1)*max(1.,col2))
c     .       *(1.+vparam(1)*abs(col1-col2))) !col1+col2

        b(1)=1.
        b(2)=2. !0.25
        b(3)=0.1*log(1.+col)
        b(4)=0.0*log(1.+col)
        x1=xmpi(x,b)

        c(1)=max(1.,col3)**0.8 !1.
        c(2)=-(2.+1.*log(1.+col))
        c(3)=0.25*log(1.+max(0.,col3-1.)+10.*col)
c        c(4)=0.25*(300.*col**2+max(0.,col3-1.)**2.)
        c(4)=0.25*(80.*col**2+max(0.,col3-1.)**2.)
        x2=xmpi(x,c)

      elseif(abs(iqq).eq.2)then

        a(1)=0.75
        a(2)=10.
        a(3)=1
        a(4)=10. !!25.
        a(5)=1.
        col=dsatur*log(1.+max(col1,0.01*col2))
c     .            *(1.+vparam(1)*0.03*abs(col1-col2)) !col1+col2


        b(1)=1.
        b(2)=1.
        b(3)=0.1*log(1.+max(0.,col3-1)+0.*col)
        b(4)=0.5*max(0.,col3+0.0*col)**0.5
        x1=xmpi(x,b)

c        frac=min(1.,max(1.,col1)/max(1.,col2))
        c(1)=max(1.,col3)**0.5
        c(2)=-2.
        c(3)=1.*log(1.+max(0.,col3-1.)+3.5*log(1.+col))
        c(4)=2.*(100.*(col)**2+max(0.,col3-1.))
c        c(4)=2.*(1000.*(col)**2+max(0.,col3-1.))
        x2=xmpi(x,c)

      else
        stop 'error 20200429 - should not be in PomInc'
      endif

      if(iqq.gt.0)then
        PomInc=fct(x1,bb)/dble(max(1.,col3)
     .        +a(2)*min(1.,col)**2*log(dble(maproj*matarg)))**a(1)
c     .                   /dble(1.+a(2)*col**a(5))
        PomInc=Pominc
     .        +fct(x2,bb)*dble(col*a(4)+max(0.,col3-1))**a(3)
c     .                   *dble(max(1.,col))**a(4) !dble(1.+a(4)*col)
        !PomInc=fct(x,bb) !Test (no deformation)
      else
        PomInc=fct(x1)/dble(max(1.,col3)
     .        +a(2)*min(1.,col)**2*log(dble(maproj*matarg)))**a(1)
c     .                /dble(1.+a(2)*col**a(5))
        PomInc=Pominc
     .        +fct(x2)*dble(col*a(4)+max(0.,col3-1.))**a(3)
c     .                *dble(max(1.,col))**a(4) !*dble(1.+a(4)*col)
        !write(*,'(a,f10.7,5f10.3)'),'PomInc()',x,fct(x),PomInc
        !.,fct(x1)/dble( 
        !.   max(1.,col3) + a(2)*min(1.,col)**2*log(dble(maproj*matarg))
        !.              )**a(1)
        !.,col, a(2)*min(1.,col)**2*log(dble(maproj*matarg))
        !PomInc=fct(x) !Test (no deformation)
      endif
      end


c----------------------------------------------------------------------
      double precision function PomInc2(fct,xp,xm,bb,col1,col2,col3)   !---MC---
c-----------------------------------------------------------------------
c     Transformation function to fit PomInc distribution with MPI
c     for xp with col1 collision, xm with col2 collisions and bb
c-----------------------------------------------------------------------
      implicit none
#include "aaa.h"
      double precision fct,xp,xm,xp1,xp2,xm1,xm2,xmpi
      real bb,col1,col2,col3,a(5),b(4),c(4),col
      external fct
      
      a(1)=1.
      a(2)=2.2
      a(3)=2.
      a(4)=3.
      a(5)=1.
      col=dsatur*log(max(1.,col1)*max(1.,col2))
c     .                 *(1.+vparam(1)*abs(col1-col2))) !col1+col2
      
      b(1)=1.
      b(2)=1.
      b(3)=0.1*log(1.+col)
      b(4)=0.0*log(1.+col)
      xp1=xmpi(xp,b)
      xm1=xmpi(xm,b)

      c(1)=max(1.,col3)**0.5
      c(2)=-0.5*(2.+1.*log(1.+col))
      c(3)=0.25*log(1.+max(0.,col3-1.)+10.*col)
      c(4)=0.25*(80.*col**2+max(0.,col3-1.)**2.)
      xp2=xmpi(xp,c)
      xm2=xmpi(xm,c)

c      b(1)=1.
c      b(2)=1.
c      b(3)=0.1*log(1.+max(0.,col3-1)+0.1*col1)
c      b(4)=0.5*max(0.,col3+0.025*col1)**0.5
c      xp1=xmpi(xp,b)
c      b(3)=0.1*log(1.+max(0.,col3-1)+0.1*col2)
c      b(4)=0.5*max(0.,col3+0.025*col2)**0.5
c      xm1=xmpi(xm,b)
c
c      c(1)=max(1.,col3)**0.5
c      c(2)=-2.
c      c(3)=1.*log(1.+max(0.,col3-1.)+0.05*col1)
c      c(4)=2.*(0.15*(col1)**2+max(0.,col3-1.))
c      xp2=xmpi(xp,c)
c      c(3)=1.*log(1.+max(0.,col3-1.)+0.05*col2)
c      c(4)=2.*(0.15*(col2)**2+max(0.,col3-1.))
c      xm2=xmpi(xm,c)

      PomInc2=fct(xp1,xm1,bb)/dble(max(1.,col3)
     .         +a(2)*min(1.,col)**2*log(dble(maproj*matarg)))**a(1)
c     .                       /dble(1.+a(2)*col**a(5))!/dble(1.+a(2)*log(1.+col))
      PomInc2=Pominc2
     .       +fct(xp2,xm2,bb)*dble(col*a(4)+max(0.,col3-1.))**a(3)
c     .                       *dble(max(1.,col))**a(4) !*dble(1.+a(4)*col)

      !PomInc2=fct(xp,xm,bb) !Test (no deformation)
      end






cc----------------------------------------------------------------------
c      double precision function PomIncII(b)   !---check---
cc----------------------------------------------------------------------
cc  integral_dx_dy om1*F_remn*F_remn for a given b   !---check---
cc----------------------------------------------------------------------
c#include "aaa.h"
c#include "ems.h"
c#include "sem.h"
c#include "par.h"
c       double precision cint,gamom(idxD0:idxD1),deltap(idxD0:idxD1)
c     &,deltapp(idxD0:idxD1),utgam2
c
cc Calculation by analytical integration (perfect but it changes
cc if om1 change):
c
c      s=engy**2
c      do i=idxDmin(iomega),idxDmax(iomega)
c        zp=0.
c        zt=0.
c        call Gfunpar(zp,zt,0.,0.,1,i,b,s,alp,bet,betp,epsp,epst,epss,gamv)
c        gamom(i)=dble(alp*gamv)
c        deltap(i)=dble(bet)
c        deltapp(i)=dble(betp)
c
cc Integration possible only if delta(i)>-1
c
c       if(deltap(i).le.-1.d0.or.deltapp(i).le.-1.d0)
c     &       stop 'Error in epos-par-300 in PomIncII'
c      enddo
c
c      cint=0.d0
c      do i=idxDmin(iomega),idxDmax(iomega)
c       cint=cint+gamom(i)*utgam2(deltap(i)+1.d0)*utgam2(deltapp(i)+1.d0)
c     &            *dble(ucfpro*ucftar)
c     &            /utgam2(dble(alplea(iclpro))+deltap(i)+2.d0)
c     &            /utgam2(dble(alplea(icltar))+deltapp(i)+2.d0)
c      enddo
c
c      PomIncII=cint
c
c      return
c      end
c

c----------------------------------------------------------------------
        double precision function PomIncXIExact(x)   !---check---
c----------------------------------------------------------------------
c integral d2b PomIncXExact
c----------------------------------------------------------------------
#include "aaa.h"
#include "par.h"
      double precision x,PomIncXExact
      common /ar3/  x1(7),a1(7)
      common/geom/rmproj,rmtarg,bmax,bkmx

      bmid=bkmx/2.
      PomIncXIExact=0.d0
      do i=1,7
        do m=1,2
          bb=bmid*(1.+(2.*m-3)*x1(i))
        PomIncXIExact=PomIncXIExact+dble(bb*a1(i))*PomIncXExact(x,bb)
        enddo
      enddo

      PomIncXIExact=PomIncXIExact*dble(2.*pi*bmid)

      return
      end

c----------------------------------------------------------------------
        double precision function PomIncXIUnit(x)   !---check---
c----------------------------------------------------------------------
c integral d2b PomIncXUnit
c----------------------------------------------------------------------
#include "aaa.h"
#include "par.h"
      double precision x,PomIncXUnit
      common /ar3/  x1(7),a1(7)
      common/geom/rmproj,rmtarg,bmax,bkmx

      bmid=bkmx/2.
      PomIncXIUnit=0.d0
      do i=1,7
        do m=1,2
          bb=bmid*(1.+(2.*m-3)*x1(i))
       PomIncXIUnit=PomIncXIUnit+dble(bb*a1(i))*PomIncXUnit(x,bb)
        enddo
      enddo

      PomIncXIUnit=PomIncXIUnit*dble(2.*pi*bmid)

      return
      end

c----------------------------------------------------------------------
        double precision function PomIncPIExact(x)   !---check---
c----------------------------------------------------------------------
c integral d2b PomIncPExact
c----------------------------------------------------------------------
#include "aaa.h"
#include "par.h"
      double precision x,PomIncPExact
      common/geom/rmproj,rmtarg,bmax,bkmx
      common /ar3/  x1(7),a1(7)

      bmid=bkmx/2.
      PomIncPIExact=0.d0
      do i=1,7
        do m=1,2
          bb=bmid*(1.+(2.*m-3)*x1(i))
       PomIncPIExact=PomIncPIExact+dble(bb*a1(i))*PomIncPExact(x,bb)
        enddo
      enddo

      PomIncPIExact=PomIncPIExact*dble(2.*pi*bmid)

      return
      end

c----------------------------------------------------------------------
        double precision function PomIncPIUnit(x)   !---check---
c----------------------------------------------------------------------
c integral d2b PomIncPUnit
c----------------------------------------------------------------------
#include "aaa.h"
#include "par.h"
      double precision x,PomIncPUnit
      common/geom/rmproj,rmtarg,bmax,bkmx
      common /ar3/  x1(7),a1(7)

      bmid=bkmx/2.
      PomIncPIUnit=0.d0
      do i=1,7
        do m=1,2
          bb=bmid*(1.+(2.*m-3)*x1(i))
          PomIncPIUnit=PomIncPIUnit+dble(bb*a1(i))*PomIncPUnit(x,bb)
        enddo
      enddo

      PomIncPIUnit=PomIncPIUnit*dble(2.*pi*bmid)

      return
      end

c----------------------------------------------------------------------
        double precision function PomIncMIExact(x)   !---check---
c----------------------------------------------------------------------
c integral d2b PomIncMExact
c----------------------------------------------------------------------
#include "aaa.h"
#include "par.h"
      double precision x,PomIncMExact
      common/geom/rmproj,rmtarg,bmax,bkmx
      common /ar3/  x1(7),a1(7)

      bmid=bkmx/2.
      PomIncMIExact=0.d0
      do i=1,7
        do m=1,2
          bb=bmid*(1.+(2.*m-3)*x1(i))
          PomIncMIExact=PomIncMIExact+dble(bb*a1(i))*PomIncMExact(x,bb)
        enddo
      enddo

      PomIncMIExact=PomIncMIExact*dble(2.*pi*bmid)

      return
      end

c----------------------------------------------------------------------
        double precision function PomIncMIUnit(x)   !---check---
c----------------------------------------------------------------------
c integral d2b PomIncMUnit
c----------------------------------------------------------------------
#include "aaa.h"
#include "par.h"
      double precision x,PomIncMUnit
      common/geom/rmproj,rmtarg,bmax,bkmx
      common /ar3/  x1(7),a1(7)

      bmid=bkmx/2.
      PomIncMIUnit=0.d0
      do i=1,7
        do m=1,2
          bb=bmid*(1.+(2.*m-3)*x1(i))
        PomIncMIUnit=PomIncMIUnit+dble(bb*a1(i))*PomIncMUnit(x,bb)
        enddo
      enddo

      PomIncMIUnit=PomIncMIUnit*dble(2.*pi*bmid)

      return
      end

c----------------------------------------------------------------------
        double precision function PomIncMExact(xm,b)   !---check---
c----------------------------------------------------------------------
c incluse Pomeron distribution \int dx+ { G F_remn F_remn }
c----------------------------------------------------------------------
#include "aaa.h"
#include "sem.h"
#include "ems.h"
      double precision AlTiP,BeTip,al,bep,bepp,xpInt,utgam2,xm

      s=engy**2
      PomIncMExact=0.d0
      do i=idxDmin(iomega),idxDmax(iomega)
        zp=0.
        zt=0.
        call Gfunpar(zp,zt,0.,0.,1,i,b,s,alp,bet,betp,epsp,epst,epss
     *              ,gamv)
        bep =dble(bet)
        bepp=dble(betp)
        al  =dble(alp*gamv)

        BeTip=bep+1.d0
        xpInt=utgam2(BeTip)*dble(ucfpro)
     *                    /utgam2(1.d0+dble(alplea(iclpro))+BeTip)
        AlTiP=al*xpInt
        PomIncMExact=PomIncMExact+AlTiP*xm**bepp
     *                            *(1.d0-xm)**dble(alplea(icltar))
      enddo

      return
      end

c----------------------------------------------------------------------
        double precision function PomIncMExactk(xm,bkk)   !---MC---
c----------------------------------------------------------------------
c incluse Pomeron distribution \int dx+ { G F_remn F_remn }
c----------------------------------------------------------------------
#include "aaa.h"
#include "sem.h"
#include "ems.h"
      double precision AlTiP,BeTip,al,bep,bepp,xpInt,utgam2,xm

      k=nint(bkk)
      PomIncMExactk=0.d0
      do i=ntymin,ntymax
        bep  = btildep(i,k)
        bepp = btildepp(i,k)
        al= atilde(i,k)

        BeTip=bep+1.d0
        xpInt=utgam2(BeTip)*dble(ucfpro)
     *                    /utgam2(1.d0+dble(alplea(iclpro))+BeTip)
        AlTiP=al*xpInt
        PomIncMExactk=PomIncMExactk+AlTiP*xm**bepp
     *                            *(1.d0-xm)**dble(alplea(icltar))
      enddo

      return
      end

c----------------------------------------------------------------------
      double precision function PomIncMUnit(xm,b)   !---check---
c----------------------------------------------------------------------
c incluse  Unitarized Pomeron distribution  \int dx+
c----------------------------------------------------------------------
#include "aaa.h"
#include "sem.h"
#include "ems.h"
      double precision Df,xp,xm,G2,w,xpm
      double precision PoInU!,Znorm
      common /ar3/  x1(7),a1(7)

      s=engy**2

c Calculation by numeric integration :
      w=0.d0
      xpm=.5d0
      do m=1,2
        do j=1,7
          xp=xpm*(1.d0+dble((2.*m-3.)*x1(j)))
          Df=0.d0
          do i=idxDmin(iomega),idxDmax(iomega)
            zp=0.
            zt=0.
            call Gfunpar(zp,zt,0.,0.,1,i,b,s,alp,bet,betp,epsp,epst,epss
     .                  ,gamv)
            Df=Df+dble(alp)*xp**dble(bet)*xm**dble(betp)
          enddo
          G2=Df
          w=w+dble(a1(j))*PoInU(xp,xm,s,b)*G2
        enddo
      enddo
      w=w*xpm


      PomIncMUnit=w!/Znorm(s,b)

      return
      end


c----------------------------------------------------------------------
      double precision function PomIncPExact(xp,b)   !---check---
c----------------------------------------------------------------------
c incluse Pomeron distribution  \int dx- { G F_remn F_remn }
c----------------------------------------------------------------------
#include "aaa.h"
#include "sem.h"
#include "ems.h"
      double precision AlTiP,BeTipp,al,bep,bepp,xmInt,utgam2,xp

      s=engy**2
      PomIncPExact=0.d0
      do i=idxDmin(iomega),idxDmax(iomega)
        zp=0.
        zt=0.
        call Gfunpar(zp,zt,0.,0.,1,i,b,s,alp,bet,betp,epsp,epst,epss
     *              ,gamv)
        bep=dble(bet)
        bepp=dble(betp)
        al=dble(alp*gamv)
        BeTipp=bepp+1.d0
        xmInt=utgam2(BeTipp)*dble(ucftar)
     *                    /utgam2(1.d0+dble(alplea(icltar))+BeTipp)
        AlTiP=al*xmInt
        PomIncPExact=PomIncPExact+AlTiP*xp**bep
     *                            *(1.d0-xp)**dble(alplea(iclpro))
      enddo

      return
      end

c----------------------------------------------------------------------
      double precision function PomIncPExactk(xp,bkk)   !---MC---
c----------------------------------------------------------------------
c incluse Pomeron distribution  \int dx- { G F_remn F_remn }
c----------------------------------------------------------------------
#include "aaa.h"
#include "sem.h"
#include "ems.h"
      double precision AlTiP,BeTipp,al,bep,bepp,xmInt,utgam2,xp

      k=nint(bkk)
      PomIncPExactk=0.d0
      do i=ntymin,ntymax
        bep  = btildep(i,k)
        bepp = btildepp(i,k)
        al   = atilde(i,k)
        BeTipp=bepp+1.d0
        xmInt=utgam2(BeTipp)*dble(ucftar)
     *                    /utgam2(1.d0+dble(alplea(icltar))+BeTipp)
        AlTiP=al*xmInt
        PomIncPExactk=PomIncPExactk+AlTiP*xp**bep
     *                            *(1.d0-xp)**dble(alplea(iclpro))
      enddo

      return
      end

c----------------------------------------------------------------------
      double precision function PomIncPUnit(xp,b)   !---check---
c----------------------------------------------------------------------
c incluse  Unitarized Pomeron distribution  \int dx-
c----------------------------------------------------------------------
#include "aaa.h"
#include "sem.h"
      double precision Df,xp,xm,G2,w,xmm
      double precision PoInU!,Znorm
      common /ar3/  x1(7),a1(7)

      s=engy**2

c Calculation by numeric integration :
      w=0.d0
      xmm=.5d0
      do m=1,2
        do j=1,7
          xm=xmm*(1.d0+dble((2.*m-3.)*x1(j)))
          Df=0.d0
          do i=idxDmin(iomega),idxDmax(iomega)
            zp=0.
            zt=0.
            call Gfunpar(zp,zt,0.,0.,1,i,b,s,alp,bet,betp,epsp,epst,epss
     .                  ,gamv)
            Df=Df+alp*xp**dble(bet)*xm**dble(betp)
          enddo
          G2=Df
          w=w+dble(a1(j))*PoInU(xp,xm,s,b)*G2
        enddo
      enddo
      w=w*xmm


      PomIncPUnit=w!/Znorm(s,b)

      return
      end


c----------------------------------------------------------------------
        double precision function PomIncJExact(b)   !---check---
c----------------------------------------------------------------------
c integral of Pomeron distribution  \int dy dx { G F_remn F_remn }
c----------------------------------------------------------------------
#include "aaa.h"
      double precision allea,PomIncXExact,xh
      common /ar3/  x1(7),a1(7)

      allea=2.d0+dble(alplea(iclpro)+alplea(icltar))
      PomIncJExact=0.d0
      do i=1,7
        do m=1,2
          xh=1.d0-(.5d0+dble(x1(i)*(float(m)-1.5)))**(1.d0/allea)
          PomIncJExact=PomIncJExact+dble(a1(i))
     &       *PomIncXExact(xh,b)/(1.d0-xh)**(allea-1.d0)
        enddo
      enddo
      PomIncJExact=PomIncJExact/allea/2.d0

      return
      end


c----------------------------------------------------------------------
        double precision function PomIncJUnit(b)   !---check---
c----------------------------------------------------------------------
c integral of Pomeron distribution  \int dy dx { G F_remn F_remn }
c----------------------------------------------------------------------
#include "aaa.h"
      double precision PomIncXUnit,xh,xhm
      common /ar3/  x1(7),a1(7)

      PomIncJUnit=0.d0
      xhm=.5d0
      do i=1,7
        do m=1,2
          xh=xhm*(1.d0+dble(x1(i)*(2.*float(m)-3.)))
          PomIncJUnit=PomIncJUnit+dble(a1(i))
     &                                *PomIncXUnit(xh,b)
        enddo
      enddo
      PomIncJUnit=PomIncJUnit*xhm

      return
      end


c----------------------------------------------------------------------
        double precision function PomIncExacth(b)   !---MC---
c----------------------------------------------------------------------
c integral of Pomeron distribution  \int dy dx { G F_remn F_remn }
c----------------------------------------------------------------------
#include "aaa.h"
#include "ems.h"
      double precision allea,PomIncXExacth,xh,xmin,xmax
      common /ar3/  x1(7),a1(7)

c define variables for for om5J and PhiUnit
      spp=engy**2
      xmin=1d0/dble(spp)
      xmax=1d0
      imax=idxD1
      if(iomega.ge.2)imax=idxD
      do i=idxD0,imax
      zp=0.
      zt=0.
      call Gfunpar(zp,zt,0.,0.,1,i,b,spp,alpx,betx,betpx,epsp,epst,epss
     .            ,gamv)
      enddo
      allea=2.d0+dble(alplea(iclpro)+alplea(icltar))
      PomIncExacth=0.d0
      do i=1,7
        do m=1,2
          xh=1.d0-(.5d0+dble(x1(i)*(float(m)-1.5)))**(1.d0/allea)
          PomIncExacth=PomIncExacth+dble(a1(i))
     &       *PomIncXExacth(b,xh)/(1.d0-xh)**(allea-1.d0)
        enddo
      enddo
      PomIncExacth=PomIncExacth/allea/2.d0

      return
      end


c----------------------------------------------------------------------
      double precision function PomIncXExacth(b,xh)   !---MC---
c----------------------------------------------------------------------
c incluse hard Pomeron distribution  \int dy { G F_remn F_remn }
c----------------------------------------------------------------------
#include "aaa.h"
#include "sem.h"
#include "ems.h"
      double precision xh,Df,xp,xm,w,ymax,om5J
      common /ar3/  x1(7),a1(7)

c Calculation by numeric integration :
      w=0.d0
      ymax=-.5d0*log(xh)
      do m=1,2
        do j=1,7
          xp=sqrt(xh)*exp(dble((2.*m-3.)*x1(j))*ymax)
          xm=xh/xp
          Df=om5J(xp,xm,b,-10,-10)
     *      *(1.d0-xp)**dble(alplea(iclpro))
     *      *(1.d0-xm)**dble(alplea(icltar))
          w=w+dble(a1(j))*Df
        enddo
      enddo
      w=w*ymax


      PomIncXExacth=w

      return
      end

c----------------------------------------------------------------------
        double precision function PomIncXExact(xh,b)   !---check---
c----------------------------------------------------------------------
c incluse Pomeron distribution  \int dy { G F_remn F_remn }
c (optimized integration but with no y dependance)
c----------------------------------------------------------------------
#include "aaa.h"
#include "sem.h"
      double precision AlTiP,bep,bepp,factor,factor1
      double precision xpmin,xh,xp,xm,ymax,y
      common /ar3/  x1(7),a1(7)

      s=engy**2
      PomIncXExact=0.d0
      do i=idxDmin(iomega),idxDmax(iomega)
        zp=0.
        zt=0.
        call Gfunpar(zp,zt,0.,0.,1,i,b,s,alpx,betx,betpx,epsp,epst,epss
     *              ,gamv)
        bep  =betx
        bepp =betpx
        AlTiP=alpx*gamv
        PomIncXExact=PomIncXExact+AlTiP*xh**((bep+bepp)/2.d0)
      enddo

      factor=0.d0
      allea=min(alplea(iclpro),alplea(icltar))+1.
      xpmin=max(sqrt(xh),exp(-1.d0))
      do i=1,7
        do m=1,2
          xp=1.d0-(1.d0-xpmin)*(.5d0+dble(x1(i)*(float(m)-1.5)))
     *                                         **(1.d0/dble(allea))
          xm=xh/xp
          factor=factor+dble(a1(i))
     *        *((1.d0-xp)**dble(alplea(iclpro)-allea+1.)
     *        *(1.d0-xm)**dble(alplea(icltar))
     *        +(1.d0-xp)**dble(alplea(icltar)-allea+1.)
     *        *(1.d0-xm)**dble(alplea(iclpro)))/xp
        enddo
      enddo
      factor=factor*(1.d0-xpmin)**dble(allea)/dble(allea)


      if(xpmin.gt.1.00001d0*sqrt(xh))then
        ymax=-log(xh)-2.d0
        factor1=0.d0
        do i=1,7
          do m=1,2
            y=ymax*dble(x1(i)*(2*m-3))
            xp=sqrt(xh*exp(y))
            xm=xh/xp
            factor1=factor1+dble(a1(i))*(1.d0-xp)**dble(alplea(iclpro))
     *                                 *(1.d0-xm)**dble(alplea(icltar))
          enddo
        enddo
        factor=factor+factor1*ymax
      endif

      factor=factor/2.d0

      PomIncXExact=PomIncXExact*factor

      return
      end

c----------------------------------------------------------------------
        double precision function PomIncXExactk(xh,bkk)   !---MC---
c----------------------------------------------------------------------
c incluse Pomeron distribution  \int dy { G F_remn F_remn }
c (optimized integration but with no y dependance)
c----------------------------------------------------------------------
#include "aaa.h"
#include "sem.h"
#include "ems.h"
      double precision AlTiP,bep,bepp,factor,factor1
      double precision xpmin,xh,xp,xm,ymax,y
      common /ar3/  x1(7),a1(7)

      k=nint(bkk)
      PomIncXExactk=0.d0
      do i=ntymin,ntymax
        bep  = btildep(i,k)
        bepp = btildepp(i,k)
        AlTiP= atilde(i,k)
        PomIncXExactk=PomIncXExactk+AlTiP*xh**((bep+bepp)/2.d0)
      enddo

      factor=0.d0
      allea=min(alplea(iclpro),alplea(icltar))+1.
      xpmin=max(sqrt(xh),exp(-1.d0))
      do i=1,7
        do m=1,2
          xp=1.d0-(1.d0-xpmin)*(.5d0+dble(x1(i)*(float(m)-1.5)))
     *                                         **(1.d0/dble(allea))
          xm=xh/xp
          factor=factor+dble(a1(i))
     *        *((1.d0-xp)**dble(alplea(iclpro)-allea+1.)
     *        *(1.d0-xm)**dble(alplea(icltar))
     *        +(1.d0-xp)**dble(alplea(icltar)-allea+1.)
     *        *(1.d0-xm)**dble(alplea(iclpro)))/xp
        enddo
      enddo
      factor=factor*(1.d0-xpmin)**dble(allea)/dble(allea)


      if(xpmin.gt.1.00001d0*sqrt(xh))then
        ymax=-log(xh)-2.d0
        factor1=0.d0
        do i=1,7
          do m=1,2
            y=ymax*dble(x1(i)*(2*m-3))
            xp=sqrt(xh*exp(y))
            xm=xh/xp
            factor1=factor1+dble(a1(i))*(1.d0-xp)**dble(alplea(iclpro))
     *                                 *(1.d0-xm)**dble(alplea(icltar))
          enddo
        enddo
        factor=factor+factor1*ymax
      endif

      factor=factor/2.d0

      PomIncXExactk=PomIncXExactk*factor

      return
      end


c----------------------------------------------------------------------
      double precision function PomIncXUnit(xh,b)   !---check---
c----------------------------------------------------------------------
c incluse  Unitarized Pomeron distribution  \int dy
c----------------------------------------------------------------------
#include "aaa.h"
#include "sem.h"
      double precision xh,Df,xp,xm,w
      double precision PoInU,ymax!,Znorm
      common /ar3/  x1(7),a1(7)

      s=engy**2
ctp060829      sy=s*sngl(xh)
c Calculation by numeric integration :
      w=0.d0
      ymax=-.5d0*log(xh)
      do m=1,2
        do j=1,7
          xp=sqrt(xh)*exp(dble((2.*m-3.)*x1(j))*ymax)
          xm=xh/xp
          Df=0.d0
          do i=idxDmin(iomega),idxDmax(iomega)
            zp=0.
            zt=0.
            call Gfunpar(zp,zt,0.,0.,1,i,b,s,alp,bet,betp,epsp,epst,epss
     .                  ,gamv)
            Df=Df+alp*xp**dble(bet)*xm**dble(betp)
          enddo
          w=w+dble(a1(j))*Df*PoInU(xp,xm,s,b)
        enddo
      enddo
      w=w*ymax


      PomIncXUnit=w!/Znorm(s,b)

      return
      end



c----------------------------------------------------------------------
      double precision function PoInU(xp,xm,s,b)   !---check---
c----------------------------------------------------------------------
c Function : PhiU(1-xp,1-xm) + /int(H(z+ + x+,z- + x-)PhiU(z+,z-)dz+dz-)
c----------------------------------------------------------------------
#include "aaa.h"
      double precision xp,xm,zp,zm,zp2,zm2,zpmin,zmmin,deltzp,deltzm
      double precision zpm,zmm,w,HrstI,PhiUnit,Hrst,eps!,PhiExact
      common /ar3/  x1(7),a1(7)
      parameter(idxD2=idxD1)
      double precision GbetUni,GbetpUni,HbetUni,HbetpUni,HalpUni
      common/DGamUni/GbetUni(  idxD0:idxD2),HbetUni(  idxD0:idxD2),
     &               GbetpUni(idxD0:idxD2),HbetpUni(idxD0:idxD2),
     &               HalpUni(idxD0:idxD2)
      double precision utgam2

      eps=1.d-25

      do i=idxDmin(iomega),idxDmax(iomega)
       zpp=0.
       zpt=0.
       call Gfunpar(zpp,zpt,0.,0.,1,i,b,s,alpx,betx,betpx,epsp,epst,epss
     .             ,gamv)
       zpp=0.
       zpt=0.
       call Gfunpar(zpp,zpt,0.,0.,2,i,b,s,alpx,betx,betpx,epsp,epst,epss
     .             ,gamv)
       HbetUni(i)=betUni(i,1)+1.d0
       HbetpUni(i)=betpUni(i,1)+1.d0
       GbetUni(i)=utgam2(HbetUni(i))
       GbetpUni(i)=utgam2(HbetpUni(i))
       HalpUni(i)=alpUni(i,1)
      enddo

      if (1.d0-xp-eps.gt.0.d0.and.1.d0-xm-eps.gt.0.d0) then
        w=0.d0
        zpmin=1.d0-xp-eps
        zmmin=1.d0-xm-eps
        zpm=.5d0*zpmin
        zmm=.5d0*zmmin
c        zpm=1d0/dble(s)
c        zmm=1d0/dble(s)
        do m=1,2
          do i=1,7
            zp=zpm+dble((2*m-3)*x1(i))*zpm
c            zp=zpm*(zpmin/zpm)**dble(.5+x1(i)*(m-1.5))
            do n=1,2
              do j=1,7
                zm=zmm+dble((2*n-3)*x1(j))*zmm
c                zm=zmm*(zmmin/zmm)**dble(.5+x1(j)*(n-1.5))
            w=w+dble(a1(i)*a1(j))*Hrst(s,b,zp+xp,zm+xm)!*zp*zm
     &                              *PhiUnit(zp,zm)
c     &                             *Phiexact(0.,0.,1.,zp,zm,s,b)
             enddo
           enddo
         enddo
       enddo
       PoInU=w*zpm*zmm!*log(zpm/zpmin)*log(zmm/zmmin)*0.25d0!
       deltzp=eps
       deltzm=eps
      else
        PoInU=0.d0
        zpmin=0.d0
        zmmin=0.d0
        deltzp=1.d0-xp
        deltzm=1.d0-xm
      endif

      w=0.d0
      zpm=0.d0
      zmm=0.d0
      if(abs(deltzp).gt.1.d-10.and.abs(deltzm).gt.1.d-10)then
      zpm=.5d0*deltzp
      zmm=.5d0*deltzm
      do m=1,2
        do i=1,7
          zp=zpmin+zpm+dble((2*m-3)*x1(i))*zpm
          do n=1,2
            do j=1,7
              zm=zmmin+zmm+dble((2*n-3)*x1(j))*zmm
              zp2=1.d0-xp-zp**2
              zm2=1.d0-xm-zm**2
              w=w+dble(a1(i)*a1(j))*HrstI(s,b,zp,zm)
     &             *PhiUnit(zp2,zm2)
c     &             *Phiexact(0.,0.,1.,zp2,zm2,s,b)
            enddo
          enddo
        enddo
      enddo
      endif

      PoInU=PoInU+w*zpm*zmm
     &           +PhiUnit(1.d0-xp,1.d0-xm)
c     &           +Phiexact(0.,0.,1.,1.d0-xp,1.d0-xm,s,b)

      return
      end

c----------------------------------------------------------------------
        double precision function PomIncExact(xp,xm,b)   !---check---
c----------------------------------------------------------------------
c inclusive Pomeron distribution  { G F_remn F_remn }
c----------------------------------------------------------------------
#include "aaa.h"
#include "sem.h"
      double precision Df,xp,xm

      Df=0.d0
      s=engy**2
      do i=idxDmin(iomega),idxDmax(iomega)
        zp=0.
        zt=0.
        call Gfunpar(zp,zt,0.,0.,1,i,b,s,alp,bet,betp,epsp,epst,epss
     *              ,gamv)
        Df=Df+alp*gamv*xp**dble(bet)*xm**dble(betp)
      enddo

      PomIncExact=Df
     *            *(1.d0-xp)**dble(alplea(iclpro))
     *            *(1.d0-xm)**dble(alplea(icltar))

      return
      end

c----------------------------------------------------------------------
        double precision function PomIncExactk(xp,xm,bkk)   !---MC---
c----------------------------------------------------------------------
c inclusive Pomeron distribution  { G F_remn F_remn }
c----------------------------------------------------------------------
#include "aaa.h"
#include "sem.h"
#include "ems.h"
      double precision Df,xp,xm

      k=nint(bkk)
      Df=0.d0
      do i=ntymin,ntymax
        Df=Df+atilde(i,k)*xp**btildep(i,k)*xm**btildepp(i,k)
      enddo

      PomIncExactk=Df
     *            *(1.d0-xp)**dble(alplea(iclpro))
     *            *(1.d0-xm)**dble(alplea(icltar))

      return
      end

c----------------------------------------------------------------------
        double precision function PomIncUnit(xp,xm,b)   !---check---
c----------------------------------------------------------------------
c inclusive Pomeron distribution  { Sum{int{G*Phi} }
c----------------------------------------------------------------------
#include "aaa.h"
#include "par.h"
#include "sem.h"
      double precision PoInU,xp,xm,om1,xh,yp


      xh=xp*xm
      yp=0.d0
      if(xm.ne.0.d0)yp=0.5d0*log(xp/xm)
      PomIncUnit=om1(xh,yp,b)*PoInU(xp,xm,engy*engy,b)


      return
      end

c----------------------------------------------------------------------
        double precision function PomIncDiffExact(xp,xm,b)   !---check---
c----------------------------------------------------------------------
c inclusive Diffractive Pomeron distribution  { G F_remn F_remn }
c----------------------------------------------------------------------
#include "aaa.h"
#include "sem.h"
      double precision xp,xm

      PomIncDiffExact=0.d0
      s=engy**2
      if(iomega.ge.2)return
      i=idxD1
      zp=0.
      zt=0.
      call Gfunpar(zp,zt,0.,0.,1,i,b,s,alp,bet,betp,epsp,epst,epss,gamv)
      PomIncDiffExact=alp*gamv*xp**dble(bet)*xm**dble(betp)
     *            *(1.d0-xp)**dble(alplea(iclpro))
     *            *(1.d0-xm)**dble(alplea(icltar))

      return
      end

c----------------------------------------------------------------------
        double precision function PomIncDiffUnit(xp,xm,b)   !---check---
c----------------------------------------------------------------------
c inclusive Pomeron distribution  { Sum{int{G*Phi} }
c----------------------------------------------------------------------
#include "aaa.h"
#include "par.h"
#include "sem.h"
      double precision PoInDiff,xp,xm,om51,xh,yp,utgam2
      parameter(idxD2=idxD1)
      double precision GbetUni,GbetpUni,HbetUni,HbetpUni,HalpUni
      common/DGamUni/GbetUni(  idxD0:idxD2),HbetUni(  idxD0:idxD2),
     &               GbetpUni(idxD0:idxD2),HbetpUni(idxD0:idxD2),
     &               HalpUni(idxD0:idxD2)
 
      PomIncDiffUnit=0d0
      if(iomega.ge.2)return
      s=engy*engy
      do i=idxDmin(iomega),idxDmax(iomega)
        zpp=0.
        zpt=0.
        call Gfunpar(zpp,zpt,0.,0.,1,i,b,s,alpx,betx,betpx,epsp,epst
     .              ,epss,gamv)
        zpp=0.
        zpt=0.
       call Gfunpar(zpp,zpt,0.,0.,2,i,b,s,alpx,betx,betpx,epsp,epst,epss
     .             ,gamv)
        HbetUni(i)=betUni(i,1)+1.d0
        HbetpUni(i)=betpUni(i,1)+1.d0
        GbetUni(i)=utgam2(HbetUni(i))
        GbetpUni(i)=utgam2(HbetpUni(i))
        HalpUni(i)=alpUni(i,1)
      enddo
      xh=xp*xm
      yp=0.d0
      if(xm.ne.0.d0)yp=0.5d0*log(xp/xm)
      PomIncDiffUnit=om51(xh,yp,b,-5,-5)*PoInDiff(xp,xm,s,b)

 
      return
      end

c----------------------------------------------------------------------
        double precision function PomIncDiffk(k)   !---xs---
c----------------------------------------------------------------------
c integral of Diff. Pomeron distribution  \int dy dx { Gdiff F_remn F_remn }
c for pair k
c----------------------------------------------------------------------
#include "aaa.h"
#include "ems.h"
      parameter(idxD2=idxD1)
      double precision GbetUni,GbetpUni,HbetUni,HbetpUni,HalpUni
      common/DGamUni/GbetUni(  idxD0:idxD2),HbetUni(  idxD0:idxD2),
     &               GbetpUni(idxD0:idxD2),HbetpUni(idxD0:idxD2),
     &               HalpUni(idxD0:idxD2)
      double precision PomIncXDiff,xh,xhm,utgam2
      common /ar3/  x1(7),a1(7)

      PomIncDiffk=0.d0
      if(iomega.ge.2)return
c prepare some variables for next calculations
      s=engy**2
      b=bk(k)
      do i=idxDmin(iomega),idxDmax(iomega)
       HbetUni(i)=btildep(i,k)+1.d0
       HbetpUni(i)=btildepp(i,k)+1.d0
       GbetUni(i)=utgam2(HbetUni(i))
       GbetpUni(i)=utgam2(HbetpUni(i))
       HalpUni(i)=atilde(i,k)
      enddo
      xhm=.5d0
      do i=1,7
        do m=1,2
          xh=xhm*(1.d0+dble(x1(i)*(2.*float(m)-3.)))
          PomIncDiffk=PomIncDiffk+dble(a1(i))
     &                                *PomIncXDiff(xh,s,b)
        enddo
      enddo
      PomIncDiffk=PomIncDiffk*xhm

      return
      end

c----------------------------------------------------------------------
        double precision function PomIncDiff(s,b)   !---xs---
c----------------------------------------------------------------------
c integral of Diff. Pomeron distribution  \int dy dx { Gdiff F_remn F_remn }
c give sigma*<number of Pom diff> and not sigma_diff 
c----------------------------------------------------------------------
#include "aaa.h"
      parameter(idxD2=idxD1)
      double precision GbetUni,GbetpUni,HbetUni,HbetpUni,HalpUni
      common/DGamUni/GbetUni(  idxD0:idxD2),HbetUni(  idxD0:idxD2),
     &               GbetpUni(idxD0:idxD2),HbetpUni(idxD0:idxD2),
     &               HalpUni(idxD0:idxD2)
      double precision PomIncXDiff,xh,xhm,utgam2
      common /ar3/  x1(7),a1(7)

      PomIncDiff=0.d0
      if(iomega.ge.2)return
c prepare some variables for next calculations
      do i=idxDmin(iomega),idxDmax(iomega)
       zpp=0.
       zpt=0.
       call Gfunpar(zpp,zpt,0.,0.,1,i,b,s,alpx,betx,betpx,epsp,epst,epss
     .             ,gamv)
       zpp=0.
       zpt=0.
       call Gfunpar(zpp,zpt,0.,0.,2,i,b,s,alpx,betx,betpx,epsp,epst,epss
     .             ,gamv)
       HbetUni(i)=betUni(i,1)+1.d0
       HbetpUni(i)=betpUni(i,1)+1.d0
       GbetUni(i)=utgam2(HbetUni(i))
       GbetpUni(i)=utgam2(HbetpUni(i))
       HalpUni(i)=alpUni(i,1)
      enddo
      xhm=.5d0
      do i=1,7
        do m=1,2
          xh=xhm*(1.d0+dble(x1(i)*(2.*float(m)-3.)))
          PomIncDiff=PomIncDiff+dble(a1(i))
     &                                *PomIncXDiff(xh,s,b)
        enddo
      enddo
      PomIncDiff=PomIncDiff*xhm

      return
      end

c----------------------------------------------------------------------
      double precision function PomIncXDiff(xh,s,b)   !---check---
c----------------------------------------------------------------------
c incluse  Unitarized Diffractive Pomeron distribution  \int dy
c be carreful : initialization of DGamUni needed before calling this function
c----------------------------------------------------------------------
#include "aaa.h"
#include "sem.h"
      double precision xh,Df,xp,xm,w
      double precision PoInDiff,ymax!,Znorm
      common /ar3/  x1(7),a1(7)

      PomIncXDiff=0d0
      if(iomega.ge.2)return
ctp060829      sy=s*sngl(xh)
c Calculation by numeric integration :
      w=0.d0
      ymax=-.5d0*log(xh)
      i=idxD1
      do m=1,2
        do j=1,7
          xp=sqrt(xh)*exp(dble((2.*m-3.)*x1(j))*ymax)
          xm=xh/xp
          Df=alpUni(i,1)*xp**betUni(i,1)*xm**betpUni(i,1)
          w=w+dble(a1(j))*Df*PoInDiff(xp,xm,s,b)
        enddo
      enddo
      w=w*ymax


      PomIncXDiff=w!/Znorm(s,b)

      return
      end

c----------------------------------------------------------------------
        double precision function PomIncXIDiff(x)   !---check---
c----------------------------------------------------------------------
c integral d2b PomIncXDiff
c----------------------------------------------------------------------
#include "aaa.h"
#include "par.h"
      double precision x,PomIncXDiff,utgam2
      common /ar3/  x1(7),a1(7)
      common/geom/rmproj,rmtarg,bmax,bkmx
      parameter(idxD2=idxD1)
      double precision GbetUni,GbetpUni,HbetUni,HbetpUni,HalpUni
      common/DGamUni/GbetUni(  idxD0:idxD2),HbetUni(  idxD0:idxD2),
     &               GbetpUni(idxD0:idxD2),HbetpUni(idxD0:idxD2),
     &               HalpUni(idxD0:idxD2)

      bmid=bkmx/2.
      PomIncXIDiff=0.d0
      s=engy**2
      if(iomega.ge.2)return
      do j=1,7
        do m=1,2
          bb=bmid*(1.+(2.*m-3)*x1(j))

c prepare some variables for next calculations
          do i=idxDmin(iomega),idxDmax(iomega)
            zpp=0.
            zpt=0.
            call Gfunpar(zpp,zpt,0.,0.,1,i,bb,s,alpx,betx,betpx,epsp
     .                  ,epst,epss,gamv)
            zpp=0.
            zpt=0.
            call Gfunpar(zpp,zpt,0.,0.,2,i,bb,s,alpx,betx,betpx,epsp
     .                  ,epst,epss,gamv)
            HbetUni(i)=betUni(i,1)+1.d0
            HbetpUni(i)=betpUni(i,1)+1.d0
            GbetUni(i)=utgam2(HbetUni(i))
            GbetpUni(i)=utgam2(HbetpUni(i))
            HalpUni(i)=alpUni(i,1)
          enddo
          PomIncXIDiff=PomIncXIDiff
     &                +dble(bb*a1(j))*PomIncXDiff(x,s,bb)
        enddo
      enddo

      PomIncXIDiff=PomIncXIDiff*dble(2.*pi*bmid)

      return
      end



c----------------------------------------------------------------------
      double precision function PoInDiff(xp,xm,s,b)   !---check---
c----------------------------------------------------------------------
c Function : PhiU(1-xp,1-xm) + /int(H(z+ + x+,z- + x-)PhiU(z+,z-)dz+dz-)
c but for diff. Pomeron only (in H)
c----------------------------------------------------------------------
#include "aaa.h"
      double precision xp,xm,zp,zm,zpmin,zmmin
      double precision zpm,zmm,w,PhiUnit,Hdiff,eps!,PhiExact
      common /ar3/  x1(7),a1(7)

      eps=1.d-15

      PoInDiff=0d0
      if(iomega.ge.2)return

      if (1.d0-xp-eps.gt.0.d0.and.1.d0-xm-eps.gt.0.d0) then
        w=0.d0
        zpmin=1.d0-xp
        zmmin=1.d0-xm
c        zpm=.5d0*zpmin
c        zmm=.5d0*zmmin
        zpm=1d0/dble(s)
        zmm=1d0/dble(s)
        do m=1,2
          do i=1,7
c            zp=zpm+dble((2*m-3)*x1(i))*zpm
            zp=zpm*(zpmin/zpm)**dble(.5+x1(i)*(m-1.5))
            do n=1,2
              do j=1,7
c                zm=zmm+dble((2*n-3)*x1(j))*zmm
                zm=zmm*(zmmin/zmm)**dble(.5+x1(j)*(n-1.5))
            w=w+dble(a1(i)*a1(j))*Hdiff(s,b,zp+xp,zm+xm)*zp*zm
     &                              *PhiUnit(zp,zm)
c     &                             *Phiexact(0.,0.,1.,zp,zm,s,b)
             enddo
           enddo
         enddo
       enddo
       PoInDiff=w*log(zpm/zpmin)*log(zmm/zmmin)*0.25d0!*zpm*zmm!
      else
       PoInDiff=0.d0
      endif

c      write(ifch,*)'Poindiff',xp,poindiff,PhiUnit(1.d0-xp,1.d0-xm)
      PoInDiff=PoInDiff+PhiUnit(1.d0-xp,1.d0-xm)
c     &           +Phiexact(0.,0.,1.,1.d0-xp,1.d0-xm,s,b)

      return
      end

c----------------------------------------------------------------------
      double precision function Hdiff(s,b,xp,xm)   !test
c----------------------------------------------------------------------
#include "aaa.h"
#include "ems.h"
#include "sem.h"
#include "par.h"
      parameter(idxD2=idxD1)
      double precision GbetUni,GbetpUni,HbetUni,HbetpUni,HalpUni
      common/DGamUni/GbetUni(  idxD0:idxD2),HbetUni(  idxD0:idxD2),
     &               GbetpUni(idxD0:idxD2),HbetpUni(idxD0:idxD2),
     &               HalpUni(idxD0:idxD2)
      double precision al(idxD0:idxD2),betp(idxD0:idxD2)
     *,z,xJrst!,ffacto
      double precision zp(idxD0:idxD2),Htmp,betpp(idxD0:idxD2)
     *,yp,ym,xp,xm
      dimension ipr(idxD0:idxD2),imax(idxD0:idxD2)

      if(idxD0.ne.0.or.idxD1.ne.3) stop "Problem in Hdiff"
      Hdiff=0.d0
      if(iomega.ge.2)return
      do i=idxD0,idxD2
        imax(i)=0
        ipr(i)=0
        zp(i)=1.d0
        al(i)=0.d0
      enddo

      if(xp.ge.0.d0.and.xm.ge.0.d0.and.xp.lt.1.d0.and.xm.le.1.d0)then

      imax1=idxD2

      do i=idxDmin(iomega),imax1
        imax(i)=6+int(log10(100.*s)/3.)
c        if(i.gt.idxD)imax(i)=imax(i)*2
        if(b.ge.2.)imax(i)=max(3,imax(i)/2)
        imax(i)=min(30,imax(i))
      enddo

      Htmp=0.d0
      do i=idxDmin(iomega),imax1
        betp(i)=HbetUni(i)
        betpp(i)=HbetpUni(i)
        al(i)=HalpUni(i)
      enddo

      do ipr0=0,imax(0)
c       write(ifch,*)'Hdiff ipr0,xp,xm :',ipr0,xp,xm,Htmp
        ipr(0)=ipr0
        zp(0)=1.d0
        if (ipr(0).ne.0) zp(0)=al(0)**ipr(0)*facto(ipr(0))
c        do ipr1=0,imax(1)
c          ipr(1)=ipr1
c          zp(1)=1.d0
c          if (ipr(1).ne.0) zp(1)=al(1)**ipr(1)*facto(ipr(1))
          do ipr2=0,imax(2)
            ipr(2)=ipr2
            zp(2)=1.d0
            if (ipr(2).ne.0) zp(2)=al(2)**ipr(2)*facto(ipr(2))
           do ipr3=0,imax(3)
            ipr(3)=ipr3
            zp(3)=1.d0
            if (ipr(3).ne.0) zp(3)=al(3)**ipr(3)*facto(ipr(3))
            if (ipr(0)+ipr(1)+ipr(2)+ipr(3).ne.0) then
c            if (ipr(0)+ipr(1)+ipr(2).ne.0) then
             yp=0.d0
             ym=0.d0
             z=1.d0
             do i=idxDmin(iomega),imax1
               yp=yp+dble(ipr(i))*betp(i)
               ym=ym+dble(ipr(i))*betpp(i)
               z=z*zp(i)
             enddo
             z=z*xJrst(xp,yp,GbetUni,ipr)
             z=z*xJrst(xm,ym,GbetpUni,ipr)
             Htmp=Htmp+z
            endif
           enddo
c          enddo
        enddo
      enddo
      Hdiff=Htmp
          
      endif


      return
      end

c----------------------------------------------------------------------
      double precision function HPhiDiff(s,b)   !---inu---
c----------------------------------------------------------------------
c  Set integrated over xp and xm (x and y) Hdiff(x,y)*Phi(x,y) for a
c  given b by gauss method
c  Doesn't seem to give proper results for xs (missing factor 2 ???) !!!
c----------------------------------------------------------------------
#include "aaa.h"
      parameter(idxD2=idxD1)
      double precision GbetUni,GbetpUni,HbetUni,HbetpUni,HalpUni
      common/DGamUni/GbetUni(  idxD0:idxD2),HbetUni(  idxD0:idxD2),
     &               GbetpUni(idxD0:idxD2),HbetpUni(idxD0:idxD2),
     &               HalpUni(idxD0:idxD2)
      double precision xmm,xp,xm,xpm,w,utgam2,PhiUnit!,PhiExact
      double precision Hdiff,eps
      common /ar3/  x1(7),a1(7)
c      common /ar9/    x9(3),a9(3)

      eps=0d0 !1.d-5
      HPhiDiff=0d0

      if(iomega.ge.2)return
c prepare some variables for next calculations
      do i=idxDmin(iomega),idxDmax(iomega)
       zpp=0.
       zpt=0.
       call Gfunpar(zpp,zpt,0.,0.,1,i,b,s,alpx,betx,betpx,epsp,epst,epss
     .             ,gamv)
       zpp=0.
       zpt=0.
       call Gfunpar(zpp,zpt,0.,0.,2,i,b,s,alpx,betx,betpx,epsp,epst,epss
     .             ,gamv)
       HbetUni(i)=betUni(i,1)+1.d0
       HbetpUni(i)=betpUni(i,1)+1.d0
       GbetUni(i)=utgam2(HbetUni(i))
       GbetpUni(i)=utgam2(HbetpUni(i))
       HalpUni(i)=alpUni(i,1)
      enddo

      w=0.d0
      xpm=1d0/dble(s)
      xmm=1d0/dble(s)
      do m=1,2
c        do i=1,3
        do i=1,7
          xp=xpm*(1d0/xpm)**dble(.5+x1(i)*(m-1.5))
c          x=xhm+dble((2*m-3)*x9(i))*xhm
c          write(ifmt,*)'HPhiInt, xp int :',x
          do n=1,2
c            do j=1,3
            do j=1,7
              xm=xmm*(1d0/xmm)**dble(.5+x1(j)*(n-1.5))
c              y=yhm+dble((2*n-3)*x9(j))*yhm
             w=w+dble(a1(i)*a1(j))*HDiff(s,b,xp,xm)
     &          *PhiUnit(xp,xm)*xp*xm
c     &          *Phiexact(0.,0.,1.,x,y,s,b)
            enddo
          enddo
        enddo
      enddo

      HPhiDiff=w*log(xmm)*log(xpm)*0.25d0

      return
      end


c----------------------------------------------------------------------
        double precision function Gammapp(sy,b,mtmp)   !---check---
c----------------------------------------------------------------------
#include "aaa.h"
#include "par.h"
#include "sem.h"
#include "ems.h"
      parameter(mmax=20)
      double precision Gtab,xxpu(mmax),xxmu(mmax),PhiExpo
      double precision Zm,xp,xm,pGtab,om1,sxp,sxm,xh,yp

      Gammapp=0.d0

      xp=1.d0
      xm=1.d0
      nmcint=20000
      nmax=nmcint
      do i=1,mmax
        xxpu(i)=0.d0
        xxmu(i)=0.d0
      enddo
      Zm=0.d0

        n=0

 10     continue
          n=n+1
          if(mod(n,10000).eq.0)write(*,*)
     &              "Calculation of Gammapp(",mtmp,")->",n
          sxp=0.d0
          sxm=0.d0
          pGtab=1.d0
          do i=1,mtmp
            rnau=rangen()!*sngl(xp-sxp)
            xxpu(i)=dble(rnau)
            sxp=sxp+xxpu(i)
            rnau=rangen()!*sngl(xm-sxm)
            xxmu(i)=dble(rnau)
            sxm=sxm+xxmu(i)
          enddo
          if(sxp.lt.xp.and.sxm.lt.xm)then
            do i=1,mtmp
              xh=xxpu(i)*xxmu(i)
              if(abs(xh).gt.1.d-30)then
                yp=0.5d0*log(xxpu(i)/xxmu(i))
              else
                yp=0.d0
              endif
              Gtab=om1(xh,yp,b)
              pGtab=pGtab*Gtab
            enddo
            Zm=Zm+pGtab*Phiexpo(0.,0.,1.,xp-sxp,xm-sxm,sy,b)
          endif
        if(n.lt.nmax)goto 10
        Zm=Zm/dble(nmax)!**2.d0*(xp*xm)**dble(mtmp)

      Gammapp=Zm

      return
      end

c----------------------------------------------------------------------
        double precision function GammaGauss(sy,b,mtmp)   !---check---
c----------------------------------------------------------------------
#include "aaa.h"
#include "par.h"
#include "sem.h"
#include "ems.h"
      parameter(mmax=3)
      common /psar7/ delx,alam3p,gam3p
      double precision xpmin,xmmin,Gst,zm(mmax),zp(mmax)
     *,xpmax,xmmax,zpmin(mmax),zmmin(mmax),zpmax(mmax)
      double precision xp,xm,pGtab,omGam,dzp(mmax),Gp1,Gm1,xmin,eps
     *,sxp,sxm,PhiExpo,zmmax(mmax),dzm(mmax),Gp2,Gm2,Gp3,Gm3,G0
c     *,PhiExact
      common /ar3/  x1(7),a1(7)
c      double precision dgx1,dga1
c      common /dga20/ dgx1(10),dga1(10)

      GammaGauss=0.d0
      xp=1.d0
      xm=1.d0
      xmin=1.d-13
      eps=1.d-15

      if(mtmp.eq.0)then
        nmax1=0
        jmax1=0
        nmax2=0
        jmax2=0
        nmax3=0
        jmax3=0
      elseif(mtmp.eq.1)then
        nmax1=2
        jmax1=7
        nmax2=0
        jmax2=0
        nmax3=0
        jmax3=0
      elseif(mtmp.eq.2)then
        nmax1=2
        jmax1=7
        nmax2=2
        jmax2=7
        nmax3=0
        jmax3=0
      elseif(mtmp.eq.3)then
        nmax1=2
        jmax1=7
        nmax2=2
        jmax2=7
        nmax3=2
        jmax3=7
      else
        write(*,*)"m not between 0 and 3, return ..."
        return
      endif

        xpmin=xmin
        xmmin=xmin
        xpmax=1.d0
        xmmax=1.d0
      do i=1,mmax
        zp(i)=0.d0
        zm(i)=0.d0
        dzp(i)=0.d0
        dzm(i)=0.d0
        zmmin(i)=0.d0
        zpmax(i)=0.d0
        zpmin(i)=0.d0
        zmmax(i)=0.d0
      enddo

        G0=1.d0

        if(mtmp.eq.0)then
          sxp=xp
          sxm=xm
          G0=Phiexpo(0.,0.,1.,sxp,sxm,sy,b)
        endif


c        write(*,*)'x+/-',xmmin,xmmax,xpmin,xpmax


        dzm(1)=0.d0
        if(abs(xmmin-xmmax).ge.eps.and.mtmp.ge.1)then
        zmmax(1)=-log(xmmin)
        zmmin(1)=-log(xmmax)
        if(abs(xmmin-xmin).lt.eps)then
          zmmin(1)=-log(min(xmmax,1.d0-xmmin-xmmin))
          zmmax(1)=-log(max(xmmin,1.d0-xmmax-xmmax))
        endif
        dzm(1)=(zmmax(1)-zmmin(1))/2.d0
        endif

        dzp(1)=0.d0
        if(abs(xpmin-xpmax).ge.eps.and.mtmp.ge.1)then
        zpmax(1)=-log(xpmin)
        zpmin(1)=-log(xpmax)
        if(abs(xpmin-xmin).lt.eps)then
          zpmin(1)=-log(min(xpmax,1.d0-xpmin-xpmin))
          zpmax(1)=-log(max(xpmin,1.d0-xpmax-xpmax))
        endif
        dzp(1)=(zpmax(1)-zpmin(1))/2.d0
        endif

c        write(*,*)'bornes1=',exp(-zpmax(1)),exp(-zpmin(1))
c     &,exp(-zmmax(1)),exp(-zmmin(1))

        Gp1=0.d0
        do np1=1,nmax1
        do jp1=1,jmax1
          zp(1)=zpmin(1)+dzp(1)*(1.d0+dble(2.*np1-3.)*dble(x1(jp1)))
          Gm1=0.d0
          if(dzm(1).eq.0.d0)then
            nmax1=1
            jmax1=1
          endif
          do nm1=1,nmax1
            do jm1=1,jmax1
              if(dzm(1).ne.0.d0)then
              zm(1)=zmmin(1)+dzm(1)*(1.d0+dble(2.*nm1-3.)*dble(x1(jm1)))
              else
              zm(1)=zp(1)
              endif

              if(mtmp.eq.1)then
              sxp=xp
              sxm=xm
              do i=1,mtmp
                sxp=sxp-exp(-zp(i))
                sxm=sxm-exp(-zm(i))
              enddo
              pGtab=1.d0
              k=0
              do l=1,mtmp
                k=k+1
                if(dzp(k).ge.0.d0.and.dzm(k).ge.0.d0)then
                  Gst=omGam(exp(-zp(k)),exp(-zm(k)),b)
                  pGtab=pGtab*Gst
                  if(Gst.eq.0.d0)
     &                write(*,*)'j1=',k,exp(-zp(k)),exp(-zm(k))
     &     ,exp(-zpmin(k)),exp(-zpmax(k)),dzp(k),dzm(k),jp1
                else
                  pGtab=0.d0
                  write(*,*)'error1 ?',dzp(k),dzm(k)
                endif
              enddo
              if(sxp.gt.0.d0.and.sxm.gt.0.d0)then
                if(dzm(1).ne.0.d0)then
                  Gm1=Gm1+pGtab*
     &dble(a1(jm1))*Phiexpo(0.,0.,1.,sxp,sxm,sy,b)
     &*exp(-zm(1))
c     &dble(a1(jm1))*Phiexact(0.,0.,1.,sxp,sxm,sy,b)
c     &*exp(-zm(1))
                  else
                  Gp1=Gp1+pGtab*
     &dble(a1(jp1))*Phiexpo(0.,0.,1.,sxp,sxm,sy,b)
     &*exp(-zp(1))
c     &dble(a1(jp1))*Phiexact(0.,0.,1.,sxp,sxm,sy,b)
c     &*exp(-zp(1))
                endif
c          write(*,*)'m=1',mtmp,Gm1,Gp1,pGtab,sxp,sxm
              endif
              endif

        dzp(2)=0.d0
        if(abs(xpmin-xpmax).ge.eps.and.mtmp.ge.2)then
        zpmin(2)=-log(min(min(xpmax,1.d0-exp(-zp(1))),
     &                     1.d0-xpmin-exp(-zp(1))))
        zpmax(2)=-log(max(xpmin,1.d0-xpmax-exp(-zp(1))))
      if(abs(xpmax+xpmax+xpmax-3.d0*dble(1./delx)).lt.eps)then
          zpmin(2)=-log(xpmax)
          zpmax(2)=-log(xpmin)
        endif
          dzp(2)=(zpmax(2)-zpmin(2))/2.d0
        endif

        dzm(2)=0.d0
        if(abs(xmmin-xmmax).ge.eps.and.mtmp.ge.2)then
            zmmin(2)=-log(min(min(xmmax,1.d0-exp(-zm(1))),
     &                     1.d0-xmmin-exp(-zm(1))))
            zmmax(2)=-log(max(xmmin,1.d0-xmmax-exp(-zm(1))))
      if(abs(xmmax+xmmax+xmmax-3.d0*dble(1./delx)).lt.eps)then
            zmmin(2)=-log(xmmax)
            zmmax(2)=-log(xmmin)
          endif
          dzm(2)=(zmmax(2)-zmmin(2))/2.d0
        endif
c        write(*,*)'bornes2=',exp(-zpmax(2)),exp(-zpmin(2))
c     &,exp(-zmmax(2)),exp(-zmmin(2)),xpmax(2),1.d0-exp(-zp(1))
c     &,1.d0-xpmin(3)-exp(-zp(1)),xpmin(2),1.d0-xpmax(3)-exp(-zp(1))

        Gp2=0.d0
        do np2=1,nmax2
        do jp2=1,jmax2
          zp(2)=zpmin(2)+dzp(2)*(1.d0+dble(2.*np2-3.)*dble(x1(jp2)))
          Gm2=0.d0
          if(dzm(2).eq.0.d0)then
            nmax2=1
            jmax2=1
          endif
          do nm2=1,nmax2
            do jm2=1,jmax2
              if(dzm(2).ne.0.d0)then
              zm(2)=zmmin(2)+dzm(2)*(1.d0+dble(2.*nm2-3.)*dble(x1(jm2)))
              else
              zm(2)=zp(2)
              endif

              if(mtmp.eq.2)then
              sxp=xp
              sxm=xm
              do i=1,mtmp
                sxp=sxp-exp(-zp(i))
                sxm=sxm-exp(-zm(i))
              enddo
              pGtab=1.d0
              k=0
              do l=1,mtmp
                k=k+1
                if(dzp(k).ge.0.d0.and.dzm(k).ge.0.d0)then
               Gst=omGam(exp(-zp(k)),exp(-zm(k)),b)
                  pGtab=pGtab*Gst
                  if(Gst.eq.0.d0)
     &                write(*,*)'j2=',k,exp(-zp(k)),exp(-zm(k))
     &     ,exp(-zpmin(k)),exp(-zpmax(k)),dzp(k),dzm(k),jp1,jp2
                else
                  pGtab=0.d0
                  write(*,*)'error2 ?',dzp(k),dzm(k)
                endif
              enddo
              if(sxp.gt.0.d0.and.sxm.gt.0.d0)then
                if(dzm(2).ne.0.d0)then
                  Gm2=Gm2+pGtab*
     &dble(a1(jm2))*Phiexpo(0.,0.,1.,sxp,sxm,sy,b)
     &*exp(-zm(2))
c     &dble(a1(jm2))*Phiexact(0.,0.,1.,sxp,sxm,sy,b,mk)
c     &*exp(-zm(2))
                  else
                  Gp2=Gp2+pGtab*
     &dble(a1(jp2))*Phiexpo(0.,0.,1.,sxp,sxm,sy,b)
     &*exp(-zp(2))
c     &dble(a1(jp2))*Phiexact(0.,0.,1.,sxp,sxm,sy,b,mk)
c     &*exp(-zp(2))
                endif
c          write(*,*)'m=2',mtmp,Gm2,Gp2,pGtab,sxp,sxm
              endif
              endif

        dzp(3)=0.d0
        if(abs(xpmin-xpmax).ge.eps.and.mtmp.ge.3)then
        zpmin(3)=-log(min(xpmax,1.d0-exp(-zp(1))-exp(-zp(2))))
        zpmax(3)=-log(xpmin)
        dzp(3)=(zpmax(3)-zpmin(3))/2.d0
        endif

        dzm(3)=0.d0
        if(abs(xmmin-xmmax).ge.eps.and.mtmp.ge.3)then
        zmmin(3)=-log(min(xmmax,1.d0-exp(-zm(1))-exp(-zm(2))))
        zmmax(3)=-log(xmmin)
        dzm(3)=(zmmax(3)-zmmin(3))/2.d0
        endif

c        write(*,*)'bornes3=',exp(-zpmax(3)),exp(-zpmin(3))
c     &,exp(-zmmax(3)),exp(-zmmin(3))

        Gp3=0.d0
        do np3=1,nmax3
        do jp3=1,jmax3
          zp(3)=zpmin(3)+dzp(3)*(1.d0+dble(2.*np3-3.)*dble(x1(jp3)))
          Gm3=0.d0
          if(dzm(3).eq.0.d0)then
            nmax3=1
            jmax3=1
          endif
          do nm3=1,nmax3
            do jm3=1,jmax3
              if(dzm(3).ne.0.d0)then
              zm(3)=zmmin(3)+dzm(3)*(1.d0+dble(2.*nm3-3.)*dble(x1(jm3)))
              else
              zm(3)=zp(3)
              endif

              sxp=xp
              sxm=xm
              do i=1,mtmp
                sxp=sxp-exp(-zp(i))
                sxm=sxm-exp(-zm(i))
              enddo
              pGtab=1.d0
              k=0
              do l=1,mtmp
                k=k+1
                if(dzp(k).ge.0.d0.and.dzm(k).ge.0.d0)then
               Gst=omGam(exp(-zp(k)),exp(-zm(k)),b)
                  pGtab=pGtab*Gst
                  if(Gst.eq.0.d0)
     &                write(*,*)'j3=',k,exp(-zp(k)),exp(-zm(k))
     &   ,exp(-zpmin(k)),exp(-zpmax(k)),dzp(k),dzm(k),jp1,jp2,jp3
                else
                  pGtab=0.d0
                  write(*,*)'error3 ?',k,dzp(k),dzm(k)
                endif
              enddo
              if(sxp.gt.0.d0.and.sxm.gt.0.d0)then
                if(dzm(3).ne.0.d0)then
                  Gm3=Gm3+pGtab
     &*dble(a1(jm3))*Phiexpo(0.,0.,1.,sxp,sxm,sy,b)
     &*exp(-zm(3))
                  else
                  Gp3=Gp3+pGtab
     &*dble(a1(jp3))*Phiexpo(0.,0.,1.,sxp,sxm,sy,b)
     &*exp(-zp(3))
                endif
              endif
            enddo
          enddo
         if(dzm(3).ne.0.d0)Gp3=Gp3+Gm3*dble(a1(jp3))*exp(-zp(3))*dzm(3)
         nmax3=2
         jmax3=7
        enddo
      enddo
              if(mtmp.gt.2.and.dzm(2).ne.0.d0)then
                Gm2=Gm2+Gp3*dble(a1(jm2))*exp(-zm(2))*dzp(3)
              elseif(mtmp.gt.2)then
                Gp2=Gp2+Gp3*dble(a1(jp2))*exp(-zp(2))*dzp(3)
              endif
            enddo
          enddo
         if(dzm(2).ne.0.d0)Gp2=Gp2+Gm2*dble(a1(jp2))*exp(-zp(2))*dzm(2)
         nmax2=2
         jmax2=7
        enddo
      enddo
              if(mtmp.gt.1.and.dzm(1).ne.0.d0)then
                Gm1=Gm1+Gp2*dble(a1(jm1))*exp(-zm(1))*dzp(2)
              elseif(mtmp.gt.1)then
                Gp1=Gp1+Gp2*dble(a1(jp1))*exp(-zp(1))*dzp(2)
              endif
            enddo
          enddo
         if(dzm(1).ne.0.d0)Gp1=Gp1+Gm1*dble(a1(jp1))*exp(-zp(1))*dzm(1)
         nmax1=2
         jmax1=7
        enddo
      enddo

      if(mtmp.gt.0)G0=Gp1*dzp(1)
      write(*,*)"int:",G0

      GammaGauss=GammaGauss+G0

      return

      end

c-----------------------------------------------------------------------
      double precision function omWi(sy,b)   !---check---
c-----------------------------------------------------------------------
c cut enhanced diagram integrated over xp, xm, xpr,xmr
c (with ideal G)
c b - impact parameter between the pomeron ends;
c sy- total energy
c-----------------------------------------------------------------------
#include "aaa.h"
#include "par.h"
#include "sem.h"
#include "ems.h"
      common /psar7/ delx,alam3p,gam3p
      double precision xpmin,xmmin,zp,zm,alpp,alpm,xjacp,xjacm
     *,xpmax,xmmax,zpmin,zmmin,zpmax,chg
      double precision xp,xm,pGtab,omGam,dzp,Gp1,Gm1,xmin,eps
     *,sxp,sxm,PhiExpo,zmmax,dzm!,gamp,gamm,gampp,gammp
c     *,PhiExact
      common /ar3/  x1(7),a1(7)
c      double precision dgx1,dga1
c      common /dga20/ dgx1(10),dga1(10)

      omWi=0.d0

        xmin=1.d-30
        eps=1.d-15
        chg=1.d0/dble(delx)
        b2=b*b
        gamb=gamD(1,iclpro,icltar)
ctp060829        gamp=dble(gamb*b2/2.-alppar)
ctp060829        gamm=dble(gamb*b2/2.-alppar)
ctp060829        gampp=1.d0+2.d0*gamp
ctp060829        gammp=1.d0+2.d0*gamm

        nmax=2
        jmax=7

        xpmin=xmin
        xmmin=xmin
        xpmax=1.d0
        xmmax=1.d0
        zpmin=0.d0
        zmmin=0.d0
        zpmax=0.d0
        zmmax=0.d0
        zp=0.d0
        zm=0.d0
        dzp=0.d0
        dzm=0.d0



        do intp=1,2
        do intm=1,2

          if(intp.eq.1)then

          xpmin=xmin
          xpmax=chg
          alpp=(1.d0+2.d0*dble(gamb*b2/2.))

          else

          xpmin=chg
          xpmax=1.d0
          alpp=1.d0!(1.d0+2.d0*dble(gamb*b2/2.))

          endif

          if(intm.eq.1)then

          xmmin=xmin
          xmmax=chg
          alpm=(1.d0+2.d0*dble(gamb*b2/2.))

          else

          xmmin=chg
          xmmax=1.d0
          alpm=1.d0!(1.d0+2.d0*dble(gamb*b2/2.))

          endif
c        write(*,*)'x+/-',intp,intm,xmmin,xmmax,xpmin,xpmax,alpp,alpm


        dzm=0.d0
        if(abs(xmmin-xmmax).ge.eps)then
          if(alpm.eq.0.d0)then
            zmmax=-log(xmmin)
            zmmin=-log(xmmax)
          else
            zmmin=xmmin**alpm
            zmmax=xmmax**alpm
          endif
          dzm=(zmmax-zmmin)/2.d0
        endif

        dzp=0.d0
        if(abs(xpmin-xpmax).ge.eps)then
          if(alpp.eq.0.d0)then
            zpmax=-log(xpmin)
            zpmin=-log(xpmax)
          else
            zpmin=xpmin**alpp
            zpmax=xpmax**alpp
          endif
          dzp=(zpmax-zpmin)/2.d0
        endif


        Gp1=0.d0

        if(abs(dzp).gt.eps.and.abs(dzm).gt.eps)then
c        write(*,*)'Ca passe ...'

        do np1=1,nmax
        do jp1=1,jmax
          zp=zpmin+dzp*(1.d0+dble(2.*np1-3.)*dble(x1(jp1)))
c          zp=zpmin+dzp*(1.d0+dble(2.*np1-3.)*dgx1(jp1))
          if(alpp.eq.0.d0)then
            xp=exp(-zp)
            xjacp=xp
          else
            xp=zp**(1.d0/alpp)
            xjacp=zp**(1.d0/alpp-1.d0)/alpp
          endif

          Gm1=0.d0
          do nm1=1,nmax
            do jm1=1,jmax
                zm=zmmin+dzm*(1.d0+dble(2.*nm1-3.)*dble(x1(jm1)))
c                zm=zmmin+dzm*(1.d0+dble(2.*nm1-3.)*dgx1(jm1))
                if(alpm.eq.0.d0)then
                  xm=exp(-zm)
                  xjacm=xm
                else
                  xm=zm**(1.d0/alpm)
                  xjacm=zm**(1.d0/alpm-1.d0)/alpm
                endif

              sxp=1.d0-xp
              sxm=1.d0-xm
              pGtab=1.d0
              if(dzp.ge.0.d0.and.dzm.ge.0.d0)then
                pGtab=omGam(xp,xm,b)
                if(pGtab.eq.0.d0)
     &          write(*,*)'j1=',xp,xm,xmmin,xmmax,dzp,dzm,jp1
              else
                pGtab=0.d0
                write(*,*)'error ?',dzp,dzm
              endif
              if(sxp.gt.0.d0.and.sxm.gt.0.d0)then
                if(dzm.ne.0.d0)then
                  Gm1=Gm1+pGtab*
     &dble(a1(jm1))*Phiexpo(0.,0.,1.,sxp,sxm,sy,b)*xjacm
c     &dga1(jm1)*Phiexpo(0.,0.,1.,sxp,sxm,sy,b)*xjacm
c     &dble(a1(jm1))*Phiexact(0.,0.,1.,sxp,sxm,sy,b)*xjacm
                  else
                  Gp1=Gp1+pGtab*
     &dble(a1(jp1))*Phiexpo(0.,0.,1.,sxp,sxm,sy,b)*xjacp
c     &dga1(jp1)*Phiexpo(0.,0.,1.,sxp,sxm,sy,b)*xjacp
c     &dble(a1(jp1))*Phiexact(0.,0.,1.,sxp,sxm,sy,b)*xjacp
                endif
c          write(*,*)'m=1',mtmp,Gm1,Gp1,pGtab,sxp,sxm
              endif

            enddo
          enddo
          if(dzm.ne.0.d0)Gp1=Gp1+Gm1*dble(a1(jp1))*dzm*xjacp
c          if(dzm.ne.0.d0)Gp1=Gp1+Gm1*dga1(jp1)*dzm*xjacp
        enddo
        enddo
        endif

        omWi=omWi+Gp1*dzp

        enddo
        enddo


      return
      end

c-----------------------------------------------------------------------
      double precision function Womint(sy,bh)   !---check---
c-----------------------------------------------------------------------
c - chi~(xp,xm) for group of cut enhanced diagram giving
c the same final state integrated over xpr and xmr (with ideal G)
c bh - impact parameter between the pomeron ends;
c xh - fraction of the energy squared s for the pomeron;
c yp - rapidity for the pomeron;
c-----------------------------------------------------------------------
#include "aaa.h"
      double precision omWi

      Womint=omWi(sy,bh)

      return
      end



c-----------------------------------------------------------------------
      double precision function WomGamint(bh)   !---check---
c-----------------------------------------------------------------------
c - chi~(xp,xm) for group of integrated cut enhanced diagram giving
c the same final for proposal.
c bh - impact parameter between the pomeron ends;
c xh - fraction of the energy squared s for the pomeron;
c yp - rapidity for the pomeron;
c-----------------------------------------------------------------------
#include "aaa.h"
      double precision omGamint

      WomGamint=omGamint(bh)

      return
      end

