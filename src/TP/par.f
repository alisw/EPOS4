C
C  This file is part of EPOS4
C  Copyright (C) 2022 research institutions and authors (See CREDITS file)
C  This file is distributed under the terms of the GNU General Public License version 3 or later
C  (See COPYING file for the text of the licence)
C


c----------------------------------------------------------------------
      subroutine paramini(imod)
c----------------------------------------------------------------------
c-----------------------------------------------------------------------
c common initialization procedure
c if isetcs = 0, alpD, betD, etc ... in inirj are not used and xkappa=1
c if isetcs = 1, alpD, betD, etc ... in inirj are not used but xkappa.ne.1
c if isetcs = 2, alpD, betD, xkappa, etc ... in inirj are used and
c                cross section from calculation in inics are read.
c    if inics doesn't exist, it produces only the calculated part of it.
c if isetcs = 3, alpD, betD, xkappa, etc ... in inirj are used and
c                cross section from simulation in inics are read.
c    if inics doesn't exist, it produces the calculated AND the
c    simulated part of it both for ionudi=1 and 3. Only the values for
c    ionudi=1 (elastic for diffraction counted in xs) are always correct.
c    AA xs with ionudi=3 do not always correspond to MC simulations.
c-----------------------------------------------------------------------
c  Set  parameters of the parametrisation of the eikonals.
c
c xDfit=Sum(i=0,idx1)(alpD(i)*xp**betDp(i)*xm**betDpp(i)*s**betD(i)
c                            *xs**(gamD(i)*b2)*exp(-b2/delD(i))
c
c Parameters stored in /Dparam/ (aaa.h)
c     if imod<=0, do settings only for iclpro (imod<0 save results in alDs)
c     if imod>0, do settings for iclpro=2 and iclpro 
c    (if imod=2 force to fit only soft part)
c----------------------------------------------------------------------

#include "aaa.h"
#include "ems.h"
#include "sem.h"
#include "par.h"
      parameter(nbkbin=100)
      common /kfitd/ xkappafit(nbkbin,nclegy,nclha,nclha),xkappa,bkbin
      double precision PhiExact,y,Dsoftshval,gzdp(0:4),gzdm(0:4)!,om51!,omdif,PhiExpo
      logical soft

      call utpri('parini',ish,ishini,3)

c Initialisation of the variables

      call Class('paramini  ')

c Variables used for xparg (only graphical application)


      emaxDf=engy            !energy for fit subroutines
      smaxDf=emaxDf**2       !energy squared for fit subroutines
c      sfshlim=100.
      nptf=10                !number of point for the fit
      xmaxDf=1.d0            !maximum limit for the fit
      call DefXminDf(dble(smaxDf))
c Remark : spmin growth like x**esatur(=0.3) while sy growth like x.
c          As a consequence if at some point xmin, spmin>(s*xmin) then
c          om=0 for any x<xmin. No fit possible below minimum spmin
c          defined by minimum x=1/E (if xp=1/s, xm=1. giving larger spmin)
      spmin=4.*q2nmin    !Definition of spmin in psvin (sem)
      sfshlim=10.*spmin          !transition energy for soft->hard
      xggfit=sfshlim/smaxDf
      xshmin=max(1d0/dble(smaxDf),1.d-1)           !minimum limit for sh variance fit
      xmxrem(1)=1.   !temporary entries
      xmxrem(2)=1.
c      soft=smaxDf.lt.sfshlim
      soft=Dsoftshval(1.d0,0.d0,0.,2)
     .     .le.0.001D0*Dsoftshval(1.d0,0.d0,0.,3)
c      soft=.true.
c      if(soft) xshmin=1.d0
      xfitmin=0.01d0   !sh fitted between f(1) and f(1)*xfitmin
      if(soft) xfitmin=1.d0

      bmaxDf=2.
     . *sqrt(4.*0.039*(r2had(iclpro)+r2had(icltar)+slopom*log(smaxDf)))!maximum b for variance fit

c security parameters
c      bxx=7.1**2                 !change minimum value of gam but should be large enough to take into account all possibilities
      betxx=-0.9              !minimum value for beta to have final beta not less than -1 even with a negative gam (happen with diffraction)

      if(iomega.gt.2)stop "iomega should be <=2 !"
      idxDmin(0)=idxD0       !minimum indice for parameters
      idxDmin(1)=idxD0       !minimum indice for parameters
      idxDmin(2)=idxD        !minimum indice for parameters no diffraction
      idxDmax(0)=idxD1       !maximum indice for parameters
      idxDmax(1)=idxD1       !maximum indice for parameters
      idxDmax(2)=idxD        !maximum indice for parameters no diffraction
      ntymin=ntymi           !minimum indice for diagram type
      ntymax=ntymx           !maximum indice for diagram type
      if(iomega.ge.2)then    !no diffraction
        ntymin=ntynd
        ntymax=ntynd
      endif

      ucfpro=utgam1(1.+alplea(iclpro))
      ucftar=utgam1(1.+alplea(icltar))

      b2xscr=1.
      fegypp=0.
      epscrx=0.
      zbcutx=0.
c initialize some variable for screening (before bmxdif amd Kfit)
      if(iscreen.eq.1)then
        fegypp=fscra(engy)
        b2xscr=4.*.0389
     &       *(r2had(iclpro)+r2had(icltar)+slopom*log(max(1.,engy**2)))
        epscrx=epscrxi
      endif

c      satrad=0.
c      if(iscreen.eq.1)then
c
cc caculate the radius where Z is saturated at epscrx to define the bases
cc of nuclear shadowing
c        rho=1.
c        fegyAA=fegypp+znurho !fscro(engy/egyscr,rho)
c        if(fegyAA.gt.epscrx)satrad=-b2xscr*log(epscrx/fegyAA)
cc        bglaubx=min(zbrmax,zbrads*sqrt(satrad))
c        bglaubx=sqrt(satrad)
c        zbcutx=0.
cc        print *,'par ---->',sqrt(satrad),bglaubx,zbcutx,epscrx
c      endif

      if(isetcs.le.1.or.imod.lt.0)call Kfit(-iclegy)             !xkappafit not used or not known : if arg<0, set xkappafit to 1

c Initialisation of the parameters

      do i=1,nbpf
        do j=1,4
          parDf(i,j)=0.
        enddo
      enddo

      if(isetcs.le.-2)return

      if(imod.eq.2.or.ishpom.eq.0)soft=.true.

      iiiclegy=iclegy

c for pi or K - p crosse section calculation, we need alpD, ... for
c iclpro=1 or 3 and iclpro=2
        iiiclpro1=iclpro
        iiidirec=1
        if(imod.le.0)then
          iiiclpro2=iclpro
        else
          iiiclpro2=2
          if(iiiclpro1.lt.iiiclpro2)iiidirec=-1
        endif
        iclprosave=iclpro

        do iiipro=iiiclpro2,iiiclpro1,iiidirec

          iclpro=iiipro

      if(ish.ge.4)write(ifch,*)'gamini'
     &    ,iscreen,iclpro,icltar,iclegy,smaxDf,sfshlim,spmin,soft,dels



c First try for fit parameters

        iscreen0=iscreen
        iomega0=iomega
        if(iomega0.lt.2)then
c maximum momentum fraction for mass (m*<<sqrt(s/(2.m.R/hbar/c)) limit due to coherence lenght in diffraction (Feinberg and Pomeranchuk)) (temporary value)
          bmx=3.*sqrt(4.*.0389*(r2had(iclpro)+r2had(icltar)
     .               +slopom*log(max(1.,smaxDf)))) !approximative size at low energy
          xmxrem(1)=min(1.,     !max(1.01*amhadr(iclpro)/e,
     .                   xmxmas/(10.135*amproj*bmx)) !) !projectile (10.135(GeV-1.fm-1)=2./hbar/c=2/0.197(GeV.fm))
          xmxrem(2)=min(1.,     !max(1.01*amhadr(icltar)/e,
     .                   xmxmas/(10.135*amtarg*bmx)) !)
          if(s0min/smaxDf.ge.min(xmxrem(1),xmxrem(2)))iomega0=2
c          print *,'--->',xmxrem
        endif

        if(isetcs.le.1.or.imod.lt.0)then !if set mode, do fit

c linear fit of the 4 components of G as a function of x
        do i=1,4
          chadr(i,iclpro,icltar)=1.   !make sure diff. renormalization is set to 1. first
        enddo
        alpsf=0.
        betsf=1.
        gamsf=0.
        delsf=1.
        alpsh=0.
        betsh=1.
        gamsh=0.
        delsh=1.
        alpsp=0.      
        betsp=0.     
        gamsp=0.
        delsp=1.
        alpst=0.      
        betst=0.     
        gamst=0.
        delst=1.
        alpsmn=0.
        betsmn=0.
        gamsmn=0.
        delsmn=1.
        alphmx=0.
        bethmx=0.
        gamhmx=0.
        delhmx=1.



c momentum fraction x dependence 
        if(iomega0.lt.2)then
          call pompar(alpst,betst,6) ! SD- (Y+) contribution
          betx=betst-alppar
          if(betx.lt.betxx)then
            if(ish.ge.1)write(ifmt,*)'Warning, bet(2) too negative !'
            alpst=alpst*xminDf**(0.5*(betx-betxx))
            betst=betxx+alppar
          endif
          if(iclpro.ne.icltar)then
            call pompar(alpsp,betsp,7) ! SD+ (Y-) contribution
            betx=betsp-alppar
            if(betx.lt.betxx)then
              if(ish.ge.1)write(ifmt,*)'Warning, bet(3) too negative !'
              alpsp=alpsp*xminDf**(0.5*(betx-betxx))
              betsp=betxx+alppar
            endif
          else
            alpsp=alpst
            betsp=betst
          endif
          parDf(2,3)=(alpsp+alpst)
          parDf(5,3)=0.5*max(betsp,betst)
        else
          parDf(2,3)=0.
          parDf(5,3)=0.
        endif
c        if(soft)then    !soft first
c          call pompar(alpsf,betsf,3) !soft (taking into R and diff)
c          if(betsf-alppar.lt.betxx)then
c            if(ish.ge.1)write(ifmt,*)'Warning, bet(0) too negative !'
c            alpsf=alpsf/xminDf**(betxx-betsf+alppar)
c            betsf=betxx+alppar
c          endif
c          parDf(1,3)=alpsf
c          parDf(4,3)=betsf
c          call pompar(alpsh,betsh,-4) !soft (minus previous fit)
c          if(betsh-alppar.lt.betxx)then
c            if(ish.ge.1)write(ifmt,*)'Warning, bet(1) too negative !'
c            alpsh=alpsh/xminDf**(betxx-betsh+alppar)
c            betsh=betxx+alppar
c          endif
c          parDf(3,3)=alpsh
c          parDf(6,3)=betsh
c        else           !hard first
          call pompar(alpsmn,betsmn,10) !soft only
          call pompar(alphmx,bethmx,2) !hard only
          call pompar(alpsh,betsh,4) !all
          betx=betsh-alppar
          if(betx.lt.betxx)then
            if(ish.ge.1)write(ifmt,*)'Warning, bet(1) too negative !'
            alpsh=alpsh*xminDf**(0.5*(betx-betxx))
            betsh=betxx+alppar
          endif
          parDf(3,3)=alpsh
          parDf(6,3)=betsh
          call pompar(alpsf,betsf,-3) !central Diff + R + soft
          betx=betsf-alppar
          if(betx.lt.betxx)then
            if(ish.ge.1)write(ifmt,*)'Warning, bet(0) too negative !'
            alpsf=alpsf*xminDf**(0.5*(betx-betxx))
            betsf=betxx+alppar
          endif
          parDf(1,3)=alpsf
          parDf(4,3)=betsf
c        endif
c          if(betsh.ge.2.
c     .   .or.parDf(3,3).lt.Dsoftshval(1.d0,0d0,0.,3))then
c            soft=.true.
c            call pompar(alpsa,betsa,10) !sa used to make better fit of soft (hard negligible)
c          else
c correction for saturation
c            call pompar(alpsa,betsa,-10) !sa
c            if(betsa.lt.2.)then
c              alpsa=sngl(Dsoftshval(1.d0,0.d0,0.,0)
c     .                  -dble(alpsh)*smaxDf**betsh
c     .                  -dble(alpsf)*smaxDf**betsf)
c     .             *smaxDf**(-betsa)
c            else
c              alpsa=0.
c              betsa=betsh*2.5          
c            endif
c          endif
c        else
cc          call pompar(alpsa,betsa,10) !sa used to make better fit of soft (hard negligible)
c          alpsh=1.e-6*alpsf     !other parameters fitted on soft part (=OK)
c          betsh=betsf     !other parameters fitted on soft part (=OK)
c          parDf(3,3)=alpsh
c          parDf(6,3)=betsh          
c        endif

c        betsh=0.31
c        alpsh0=sngl(Dsoftshval(1d0,0.d0,0.,2))*smaxDf**(-betsh)
c        alpsh1=sngl(Dsoftshval(1d0,0.d0,0.,0))*smaxDf**(-betsh)
c        if(alpsh0.lt.alpsf)alpsh1=alpsh0
c        if(alpsh0*smaxDf**betsh.lt.alpsf*smaxDf**betsf)alpsh1=alpsh0
c        if(smaxDf.gt.100.*sfshlim)then
c          xfmin=1.e-2
c          alpsh2=sngl(Dsoftshval(dble(xfmin),0.d0,0.,0))
c     &             *(xfmin*smaxDf)**(-betsh)
c        else
c          alpsh2=-1e10
c        endif
c        alpsh=max(alpsh1,alpsh2)
c        alpsh=max(alpsh0,alpsh2)
c Gaussian fit of the 2 components of G as a function of x and b

c Impact parameter dependence
c It is more stable to change gam to avoid bet to be to close to -1 than 
c changing bet according to gam (bet will be too limited if gam is too negative)

c diffractive contribution
        if(iomega0.lt.2)then
          call variance(delst,gamst,6) ! SD- (Y+) contribution
c          if(gamst.lt.0)then
c            delst=sigma2(sqrt(XminDf),6)
c            gamst=max(gamst,(alppar-0.99-betst)/bxx)
c          endif
          if(iclpro.ne.icltar)then
            call variance(delsp,gamsp,7) ! SD+ (Y-) contribution
c            if(gamsp.lt.0)then
c              delsp=sigma2(sqrt(XminDf),7)
c              gamsp=max(gamsp,(alppar-0.99-betsp)/bxx)
c            endif
          else
            gamsp=gamst
            delsp=delst
          endif
        endif        

        if(isopom.ne.0.or.(.not.soft.and.ishpom.ne.0))then
          call variance(delsmn,gamsmn,10) !soft only
          call variance(delhmx,gamhmx,2)  !hard only
          call variance(delsh,gamsh,4)
          gamsf=gamsh
          delsf=delsh
          if(iomega0.lt.2)then      !use initial values from diff. first
            delsf=delst
          if(r2hads(1)*r2hads(2).gt.0.)delsf=delsf*(r2hads(2)/r2hads(1))
            gamsf=gamst
          elseif(iregge.ne.0.or.iomega.lt.2)then
            call variance(delsf,gamsf,3) !use iomega and not iomega0 because DD is not affected by xmxrem
          endif
        endif
c        if(gamsf.lt.0)then
c          delsf=sigma2(sqrt(XminDf),3)
c          gamsf=max(gamsf,(alppar-0.99-betsf)/bxx)
c        endif
c        if(gamsh.lt.0)then
c          delsh=sigma2(1d0,4)
c          gamsh=max(gamsh,(alppar-0.99-betsh)/bxx)
c        endif

c        betsh=0.3
        
c Fit GFF


c       fit parameters
        numminDf=4             !minimum indice
        numparDf=6             !maximum indice
        betac=1.               !temperature for chi2
        betae=100.               !temperature for error
        fparDf=1.             !variation amplitude for range

        nmcxDf=5000          !number of try for x fit


        if(.not.soft)then
cc          call paramx           !alpD and betD
c          alpini=alpsf
c          call pompar(alpsf,betsf,-3) !soft (taking into account hard)
c          parDf(1,3)=max(parDf(1,3)/3.,alpsf)    !do not use 0 as minimum because it is unstable
c          frac=parDf(1,3)/alpini
c          if(isopom.ne.0.and.frac.lt.1.)then
c            call variance(delnew,gamnew,-1)
c            gamsf=max(0.,frac*gamsf+(1.-frac)*gamnew)
c            delsf=frac*delsf+(1.-frac)*delnew
c          endif
ctp          parDf(1,3)=max(parDf(1,3)/3.
ctp     .                  ,sngl(Dsoftshval(xmaxDf,0.d0,0.,3))-alpsh)
c          if(betsf.lt.2.)parDf(4,3)=betsf
c          parDf(4,3)=betsf

cc        if(smaxDf.ge.sfshlim)then
cc          call pompar(alpsf,betsf,-1) !soft (taking into account hard)
cc          parDf(1,3)=max(1e-6,alpsf-alpsh)

        endif

c      if(soft.or.
c     .dble(alpsmn)*xminDf**betsmn.ge.dble(alpsf)*xminDf**betsf)then !low diff
        alpsmn=0.
        betsmn=0.
        gamsmn=0.
        delsmn=1.
c      endif

c once the fit is done without diff renormalization, we can use the parameters to compute diff cross-section and fix renormalization
        if(iomega.lt.2)then

c do fine adjustement in 2 steps for rescaling factor of diff xs (first without
c screening, second with screening)
         alpst0=alpst
         alpsp0=alpsp
         do nn=0,1
          alpD(idxD0,iclpro,icltar)=alpsf
          alpDp(idxD0,iclpro,icltar)=alpsmn
          alpDpp(idxD0,iclpro,icltar)=alphmx
          betD(idxD0,iclpro,icltar)=0.
          betDp(idxD0,iclpro,icltar)=betsf
          betDpp(idxD0,iclpro,icltar)=betsf
          gamD(idxD0,iclpro,icltar)=gamsf
          delD(idxD0,iclpro,icltar)=delsf
          
          alpD(1,iclpro,icltar)=alpsh
          alpDp(1,iclpro,icltar)=betsmn
          alpDpp(1,iclpro,icltar)=bethmx
          betD(1,iclpro,icltar)=0.
          betDp(1,iclpro,icltar)=betsh
          betDpp(1,iclpro,icltar)=betsh
          gamD(1,iclpro,icltar)=gamsh
          delD(1,iclpro,icltar)=delsh
          
          alpD(2,iclpro,icltar)=alpst
          alpDp(2,iclpro,icltar)=gamsmn
          alpDpp(2,iclpro,icltar)=gamhmx
          betD(2,iclpro,icltar)=0.
          betDp(2,iclpro,icltar)=betst
          betDpp(2,iclpro,icltar)=0. !alpdif
          gamD(2,iclpro,icltar)=gamst
          delD(2,iclpro,icltar)=delst
          
          alpD(3,iclpro,icltar)=alpsp
          alpDp(3,iclpro,icltar)=delsmn
          alpDpp(3,iclpro,icltar)=delhmx
          betD(3,iclpro,icltar)=0.
          betDp(3,iclpro,icltar)=0. !alpdif
          betDpp(3,iclpro,icltar)=betsp
          gamD(3,iclpro,icltar)=gamsp
          delD(3,iclpro,icltar)=delsp

c compute true diff cross-sections without screening effect to have more stable results
          if(nn.eq.0)iscreen=0
          call sigmadiff(gzdp,1,1)
c compute diff cross-sections with factor x**alpdif
          call sigmadiff(gzdm,-1,1)
          if(nn.eq.0)iscreen=iscreen0
c compute renormalization factor
          do i=1,4
            if(gzdm(i).gt.0d0)chadr(i,iclpro,icltar)=gzdp(i)/gzdm(i)
     .                                       *chadr(i,iclpro,icltar)
c            print *,'ini',nn,i,gzdp(i)/gzdm(i)
c     .   ,(chadr(i,iclpro,icltar))**(-1./alpdif)
          enddo

c rescale SD components
          alpst=chadr(1,iclpro,icltar)*alpst0
          alpsp=chadr(2,iclpro,icltar)*alpsp0
          parDf(2,3)=(alpsp+alpst)
c refit hard-SD component
c          if(Dsoftshval(1d0,0d0,0.,-4).gt.0d0)then
c            call pompar(alpsh,betsh,-4) !all
c            betx=betsh-alppar
c            if(betx.lt.betxx)then
c              if(ish.ge.1)write(ifmt,*)'Warning, bet(1) too negative !'
c              alpsh=alpsh*xminDf**(0.5*(betx-betxx))
c              betsh=betxx+alppar
c            endif
c            parDf(3,3)=alpsh
c            parDf(6,3)=betsh
c          endif
c refit DD+CD+soft-hard component
          call pompar(alpsf,betsf,-3) !central Diff + R + soft - hard
          betx=betsf-alppar
          if(betx.lt.betxx)then
            if(ish.ge.1)write(ifmt,*)'Warning, bet(0) too negative !'
            alpsf=alpsf*xminDf**(0.5*(betx-betxx))
            betsf=betxx+alppar
          endif
c          call variance(delsf,gamsf,3)
c          if(gamsf.lt.0)then
c            delsf=sigma2(sqrt(XminDf),3)
c            gamsf=max(gamsf,(alppar-0.99-betsf)/bxx)
c          endif
          parDf(1,3)=alpsf
          parDf(4,3)=betsf
c          call variance(delsf,gamsf,8)
          call variance(delsf,dum,8)
         enddo
c final value of alpD, ... saved later in this subroutine

        endif

      else                      !else parameters from table (cs.i)
        nbpsf=iDxD0
        if(iclegy2.gt.1)then
          al=1.+(log(engy)-log(egylow))/log(egyfac) !energy class
          i2=min(iiiclegy+1,iclegy2)
          i1=i2-1
        else
          i1=iclegy
          i2=iclegy
          al=float(iclegy)
        endif
        dl=al-i1
        dl1=max(0.,1.-dl)
                                !linear interpolation
        alpsf=alpDs(nbpsf,i2,iclpro,icltar)*dl
     &       +alpDs(nbpsf,i1,iclpro,icltar)*dl1
        alpsh=alpDs(1,i2,iclpro,icltar)*dl
     &       +alpDs(1,i1,iclpro,icltar)*dl1
        alpsmn=alpDps(nbpsf,i2,iclpro,icltar)*dl
     &        +alpDps(nbpsf,i1,iclpro,icltar)*dl1
        alphmx=alpDpps(nbpsf,i2,iclpro,icltar)*dl
     &        +alpDpps(nbpsf,i1,iclpro,icltar)*dl1
        betsf=betDps(nbpsf,i2,iclpro,icltar)*dl
     &       +betDps(nbpsf,i1,iclpro,icltar)*dl1
        betsh=betDps(1,i2,iclpro,icltar)*dl
     &       +betDps(1,i1,iclpro,icltar)*dl1
        betsmn=alpDps(1,i2,iclpro,icltar)*dl
     &        +alpDps(1,i1,iclpro,icltar)*dl1
        bethmx=alpDpps(1,i2,iclpro,icltar)*dl
     &        +alpDpps(1,i1,iclpro,icltar)*dl1
        gamsf=gamDs(nbpsf,i2,iclpro,icltar)*dl
     &       +gamDs(nbpsf,i1,iclpro,icltar)*dl1
        gamsh=gamDs(1,i2,iclpro,icltar)*dl
     &       +gamDs(1,i1,iclpro,icltar)*dl1
        gamsmn=alpDps(2,i2,iclpro,icltar)*dl
     &       +alpDps(2,i1,iclpro,icltar)*dl1
        gamhmx=alpDpps(2,i2,iclpro,icltar)*dl
     &       +alpDpps(2,i1,iclpro,icltar)*dl1
        delsf=delDs(nbpsf,i2,iclpro,icltar)*dl
     &       +delDs(nbpsf,i1,iclpro,icltar)*dl1
        delsh=delDs(1,i2,iclpro,icltar)*dl
     &       +delDs(1,i1,iclpro,icltar)*dl1
        delsmn=alpDps(3,i2,iclpro,icltar)*dl
     &        +alpDps(3,i1,iclpro,icltar)*dl1
        delhmx=alpDpps(3,i2,iclpro,icltar)*dl
     &        +alpDpps(3,i1,iclpro,icltar)*dl1

        alpst=alpDs(2,i2,iclpro,icltar)*dl
     &       +alpDs(2,i1,iclpro,icltar)*dl1

        betst=betDps(2,i2,iclpro,icltar)*dl
     &       +betDps(2,i1,iclpro,icltar)*dl1
        gamst=gamDs(2,i2,iclpro,icltar)*dl
     &       +gamDs(2,i1,iclpro,icltar)*dl1
        delst=delDs(2,i2,iclpro,icltar)*dl
     &       +delDs(2,i1,iclpro,icltar)*dl1
        alpsp=alpDs(3,i2,iclpro,icltar)*dl
     &       +alpDs(3,i1,iclpro,icltar)*dl1
        betsp=betDpps(3,i2,iclpro,icltar)*dl
     &       +betDpps(3,i1,iclpro,icltar)*dl1
        gamsp=gamDs(3,i2,iclpro,icltar)*dl
     &       +gamDs(3,i1,iclpro,icltar)*dl1
        delsp=delDs(3,i2,iclpro,icltar)*dl
     &       +delDs(3,i1,iclpro,icltar)*dl1

        do i=1,4
        chadr(i,iclpro,icltar)=chadrs(i,i2,iclpro,icltar)*dl
     &                        +chadrs(i,i1,iclpro,icltar)*dl1
        enddo

      endif

c     if energy too small to have semi-hard interaction -> only soft diagram

c      if(soft.and.idxD0.eq.0)then !no hard
c        alpsh=1.e-6*alpsf   !other parameters fitted on soft part (=OK)
c        gamsh=gamsf
c        delsh=delsf
c      endif


c Record parameters


      alpD(idxD0,iclpro,icltar)=alpsf
      alpDp(idxD0,iclpro,icltar)=alpsmn
      alpDpp(idxD0,iclpro,icltar)=alphmx
      betD(idxD0,iclpro,icltar)=0.
      betDp(idxD0,iclpro,icltar)=betsf
      betDpp(idxD0,iclpro,icltar)=betsf
      gamD(idxD0,iclpro,icltar)=gamsf
      delD(idxD0,iclpro,icltar)=delsf

      alpD(1,iclpro,icltar)=alpsh
      alpDp(1,iclpro,icltar)=betsmn
      alpDpp(1,iclpro,icltar)=bethmx
      betD(1,iclpro,icltar)=0.
      betDp(1,iclpro,icltar)=betsh
      betDpp(1,iclpro,icltar)=betsh
      gamD(1,iclpro,icltar)=gamsh
      delD(1,iclpro,icltar)=delsh

      alpD(2,iclpro,icltar)=alpst
      alpDp(2,iclpro,icltar)=gamsmn
      alpDpp(2,iclpro,icltar)=gamhmx
      betD(2,iclpro,icltar)=0.
      betDp(2,iclpro,icltar)=betst
      betDpp(2,iclpro,icltar)=alpdif
      gamD(2,iclpro,icltar)=gamst
      delD(2,iclpro,icltar)=delst

      alpD(3,iclpro,icltar)=alpsp
      alpDp(3,iclpro,icltar)=delsmn
      alpDpp(3,iclpro,icltar)=delhmx
      betD(3,iclpro,icltar)=0.
      betDp(3,iclpro,icltar)=alpdif
      betDpp(3,iclpro,icltar)=betsp
      gamD(3,iclpro,icltar)=gamsp
      delD(3,iclpro,icltar)=delsp

c For the Plots
      parDf(1,3)=alpsf
      parDf(4,3)=betsf
      parDf(6,3)=betsh
      parDf(3,3)=alpsh
      parDf(5,3)=0.5*max(betsp,betst)
      parDf(2,3)=alpst+alpsp

c compute true diff cross-sections
c          call sigmadiff(gzdp,1)
c compute diff cross-sections with factor x**alppar
c          call sigmadiff(gzdm,-10)
c comput renormalization factor
c          do i=1,4
c            print *,'fin',i,gzdp(i)/gzdm(i)
c          enddo

c      idqp=iclpro
c      if(idqp.eq.3)iqdp=1
c      idqt=icltar
c      if(idqt.eq.3)iqdt=1
cc      alphigh(1)=alppom !2.*(1.-alppom+alppar)
cc      alphigh(2)=alppom !2.*(1.-alppom+alppar)
c      rexndii(iclpro)=rexpdif(iclpro)+rexdif(iclpro)+rexndi(iclpro)
c      rexndii(icltar)=rexpdif(icltar)+rexdif(icltar)+rexndi(icltar)
c      rexdip=0.5*(alphigh(2)+alphigh(1))
cc      rexdit=alphigh(2)
cc      rexdip=max(-0.99,min(rexdip,rexdit))
c      rexdit=rexdip
c      reg2diff=r3pom**2*wdiff(iclpro)*wdiff(icltar)*smaxDf**alpdif
c     .        *sngl(xminDf)**(-alppar)/gamreg
c      if(reg2diff.lt.1.)then
c        rexdip=-alpreg+reg2diff**edifac*(rexdip+alpreg)
cc        print *,'ici',reg2diff,rexdip,rexdit
c        rexdit=rexdip
c      endif
c      if(iomega0.lt.2)then
cc initiale values for chadr (needed for the fit)
c        if(rexdip.le.-1..or.rexdit.le.-1.)call utstop
c     .  ('alpreg or alphigh>1 not allowed !&')
c        gamD(idxD1,iclpro,icltar)=0.
c        delD(idxD1,iclpro,icltar)=4.*.0389*(gwidth*
c     &    (rexres(iclpro)+rexres(icltar))+slopoms*log(smaxDf))
cc        delD(idxD1,iclpro,icltar)=4.*.0389*(gwidth*
cc     &    (rexres(iclpro)+rexres(icltar)))
cc        gamD(idxD1,iclpro,icltar)=4.*.0389*slopoms
cc     &                           /delD(idxD1,iclpro,icltar)**2
c        alpDp(idxD1,iclpro,icltar)=0.
c        alpDpp(idxD1,iclpro,icltar)=0.
c        betD(idxD1,iclpro,icltar)=0.
c        betDp(idxD1,iclpro,icltar)=rexdit
c        betDpp(idxD1,iclpro,icltar)=rexdip
c        alpdifs=0.5*(betDpp(idxD1,iclpro,icltar)
c     .              +betDp(idxD1,iclpro,icltar))
c        betDp(idxD1,iclpro,icltar)=betDp(idxD1,iclpro,icltar)+alppar
c        betDpp(idxD1,iclpro,icltar)=betDpp(idxD1,iclpro,icltar)+alppar
cc maximum momentum fraction for mass (m*<<sqrt(s/(2.m.R/hbar/c)) limit due to coherence lenght in diffraction (Feinberg and Pomeranchuk)) (temporary value)
c        bmx=3.*sqrt(delD(idxD1,iclpro,icltar)) !approximative size at low energy
c        xmxrem(1)=min(1.,       !max(1.01*amhadr(iclpro)/e,
c     .                   xmxmas/(10.135*amproj*bmx)) !) !projectile (10.135(GeV-1.fm-1)=2./hbar/c=2/0.197(GeV.fm))
c        xmxrem(2)=min(1.,       !max(1.01*amhadr(icltar)/e,
c     .                   xmxmas/(10.135*amtarg*bmx)) !)
cc to simplify calculation take the value from the direct calculation of all diff contributions (it includes a factor s**adif which is included in alpdifs in om51 for diff contribution but not in betD because of the reggeon calculation)
cc here it allows us to suppress SD at low energy using using xmxrem (not used
cc otherwise to avoid too large contribution of diff diagram in cross section
cc at high energy)
c        omdif=om51(xminDf,0d0,0.,6,9)
c        alpD(idxD1,iclpro,icltar)=sngl(omdif*xminDf**(-alpdifs))
cc        print *,'---->',omdif,alpdifs,xminDf**(-alpdifs),xmxrem
cc     .  ,alpD(idxD1,iclpro,icltar)*xminDf**(-alppar)
c
c      else
c        alpD(idxD1,iclpro,icltar)=0.
c        alpDp(idxD1,iclpro,icltar)=0.
c        alpDpp(idxD1,iclpro,icltar)=0.
c        betD(idxD1,iclpro,icltar)=0.
c        betDp(idxD1,iclpro,icltar)=0.
c        betDpp(idxD1,iclpro,icltar)=0.
c        gamD(idxD1,iclpro,icltar)=0.
c        delD(idxD1,iclpro,icltar)=1.
c      endif

c bcut for diffractive pomeron
c      bmxdif(iclpro,icltar)=conbmxdif()     !important to do it before kfit,  b

c define b where there is saturation
      if(iscreen.eq.1)zbcutx=bcutz(smaxDf,zbrads,abs(epscrx))

      if(isetcs.le.1.or.imod.lt.0)then
        if(isetcs.gt.0)call Kfit(iclegy)
c for plots or table record alpDs, betDs, etc ...
        do i=idxD0,idxD1
         alpDs(i,iclegy,iclpro,icltar)=alpD(i,iclpro,icltar)
         alpDps(i,iclegy,iclpro,icltar)=alpDp(i,iclpro,icltar)
         alpDpps(i,iclegy,iclpro,icltar)=alpDpp(i,iclpro,icltar)
         betDs(i,iclegy,iclpro,icltar)=betD(i,iclpro,icltar)
         betDps(i,iclegy,iclpro,icltar)=betDp(i,iclpro,icltar)
         betDpps(i,iclegy,iclpro,icltar)=betDpp(i,iclpro,icltar)
         gamDs(i,iclegy,iclpro,icltar)=gamD(i,iclpro,icltar)
         delDs(i,iclegy,iclpro,icltar)=delD(i,iclpro,icltar)
        enddo
        do i=1,4
          chadrs(i,iclegy,iclpro,icltar)=chadr(i,iclpro,icltar)
        enddo

      endif


c Print results

      if(ish.ge.4)then
        write(ifch,*)"parameters for iclpro:",iclpro
        write(ifch,*)"alp,bet,gam,del sf:",alpsf,betsf,gamsf,delsf
        write(ifch,*)"alp,bet,gam,del sh:",alpsh,betsh,gamsh,delsh
        write(ifch,*)"alp,bet,gam,del st:",alpst,betst,gamst,delst
        write(ifch,*)"alp,bet,gam,del sp:",alpsp,betsp,gamsp,delsp
        write(ifch,*)"alp,bet,gam,del mn:",alpsmn,betsmn,gamsmn,delsmn
        write(ifch,*)"alp,bet,gam,del mx:",alphmx,bethmx,gamhmx,delhmx
        write(ifch,*)'diff renorm.:',(chadr(i,iclpro,icltar),i=1,4)
        write(ifch,*)'xmxrem:',xmxrem
        write(ifch,*)'screening saturation:',epscrx,b2xscr,zbcutx
c        write(ifch,*)'nuclear shadowing:',sqrt(satrad),bglaubx,zbcutx
        y=Phiexact(0.,0.,1.,1.d0,1.d0,smaxDf,0.)
        write(ifch,*)'PhiExact=',y
        if(ish.ge.4)write(ifch,'(a,12f8.4)')
     *   'xkappafit:',(xkappafit(m,iclegy,iclpro,icltar),m=1,12)
      endif

      enddo

      
      if(iclpro.ne.iclprosave)stop'problem in parini with iclpro'


      call utprix('parini',ish,ishini,3)

      return
      end

c----------------------------------------------------------------------
      subroutine Class(text)
c----------------------------------------------------------------------
#include "aaa.h"
#include "par.h"
      parameter (eps=1.e-5)    !to correct for precision problem)
      character*10 text
      if(iclegy1.eq.iclegy2)then
        iclegy=iclegy1
      else
      iclegy=1+int( (log(engy)-log(egylow))/log(egyfac) + eps ) !energy class
      if(iclegy.gt.iclegy2)then
         write(ifch,*)'***********************************************'
         write(ifch,*)'Warning in ',text
         write(ifch,*)'Energy above the range used for the fit of D:'
         write(ifch,*)egylow*egyfac**(iclegy1-1),egylow*egyfac**iclegy2
         write(ifch,*)'***********************************************'
         iclegy=iclegy2
      endif
      if(iclegy.lt.iclegy1)then
         write(ifch,*)'***********************************************'
         write(ifch,*)'Warning in ',text
         write(ifch,*)'Energy below the range used for the fit of D:'
         write(ifch,*)egylow*egyfac**(iclegy1-1),egylow*egyfac**iclegy2
         write(ifch,*)'***********************************************'
         iclegy=iclegy1
      endif
      endif
      end


c----------------------------------------------------------------------

      subroutine DefXminDf(sy)

c----------------------------------------------------------------------
c  Return the power beta and the factor alpha of the fit of the eikonal
c of a pomeron of type iqq : D(X)=alpha*(X)**beta
c----------------------------------------------------------------------
      implicit none
#include "aaa.h"
#include "par.h"
      double precision sy

      xminDf=1.d0/sy            !minimum mass to get excitation 
      if(iregge.ne.0.and.s0min.lt.0.5d0)then
        xminDf=dble(s0min)*xminDf !minimum limit=s0min GeV**2 for Reggeon
      elseif(iomega.lt.2)then
        xminDf=dble(s0min)*xminDf !minimum limit=s0min GeV**2 for diff
c        xminDf=0.5d0*xminDf     !minimum limit=0.5 GeV**2 for diff
      endif
        
      end
      
c----------------------------------------------------------------------

        subroutine pompar(alpha,beta,iqq)

c----------------------------------------------------------------------
c  Return the power beta and the factor alpha of the fit of the eikonal
c of a pomeron of type iqq : D(X)=alpha*(X)**beta
c----------------------------------------------------------------------

#include "aaa.h"
#include "sem.h"
#include "par.h"
      double precision X,D1,D0,X0,D,droot,Y
      double precision Dsoftshval,xmax
      double precision xlnXs(maxdataDf),xlnD(maxdataDf),sigma(maxdataDf)

      do i=1,nptf
        sigma(i)=1.e-2
      enddo

      iscr=iqq
      if(abs(iqq).eq.3) then       !small x (soft, Reggeon, CD, DD)

        X0=0.95d0/dble(smaxDf)!sqrt(dble(s0min))/dble(smaxDf)
        X=1.05d0/dble(smaxDf)
c        if(Dsoftshval(X,0.d0,0.,3)
c     *    -Dsoftshval(X0,0.d0,0.,3).gt.0d0.or.X.ge.X0)then
c        if(Dsoftshval(X0,0.d0,0.,iscr).gt.0d0.or.X0.le.XminDf)then
        if(X0.le.XminDf)then
          X0=X
          xmax=min(1d0,3d0*X0)
        else       !at low energy do not substract hard contribution to fit R
          iscr=abs(iscr)
          X0=xminDf
          xmax=min(0.9d0/dble(smaxDf),3d0*X0,1d0)
        endif
        if(ish.ge.4)write (ifch,'(1x,a,i2,a,1p,2e13.6)')
     *                    'pompar (',iscr,') x0,xmax=',X0,xmax

        npt=0
        do i=0,nptf-1
          X=X0
          if (i.ne.0) X=X*(xmax/X0)**(dble(i)/dble(nptf-1))
          D=max(1d-20,Dsoftshval(X,0.d0,0.,iscr))
          if(D.le.1.d-20)then
            if(ish.ge.1.and.iscr.ge.0)write(ifmt,*)
     &    "Warning in pompar ! Dsoftshval(0) could be negative"
c            sigma(i+1)=1.d5
          else
            npt=npt+1
            xlnXs(npt)=log(X)
            xlnD(npt)=log(D)
          endif
        enddo

c Fit of D(X) between X0 and xmax

        if(npt.gt.3)then
          call fit(xlnXs,xlnD,npt,sigma,1,a,beta)
          if(beta.gt.10.)beta=10.

          alpha=exp(min(50.,max(-50.,a)))
c         alpha=sngl(Dsoftshval(X0,0.d0,0.,iscr))
c     &       *(sngl(X0)*smaxDf)**(-beta)

        else   !nothing to fit
          beta=10.
          alpha=0.
        endif
        

      elseif(abs(iqq).eq.2.and.xfitmin.ne.1.d0) then    !high x (sh)

c        xmax=max(0.01d0,min(1d0,dble(sfshlim*100./smaxDf)))
        dec=1000.
        xmax=1d0 !max(0.05d0,min(1d0,dble(xggfit*dec)))  !fit dec above cut

c Definition of D0=D(X0)

        D1=Dsoftshval(xmax,0.d0,0.,iscr)

        if(D1.gt.1d-6)then
          D0=xfitmin*D1

c Calculation of X0 and D(X)

          X0=droot(D0,D1,xmax,iscr)
c          if(D0.gt.Dsoftshval(xminDf,0.d0,0.,3))then
c            X0=min(max(xmax*0.01d0,xggfit),X0)
c          endif
          
        else
          X0=xmax
        endif
        if(ish.ge.4)write (ifch,'(1x,a,i2,a,1p,2e13.6)')
     *                    'pompar (',iqq,') x0,xmax=',X0,xmax
        
        npt=0
        do i=0,nptf-1
          X=X0
          if (i.ne.0) X=X*(xmax/X0)**(dble(i)/dble(nptf-1))
          D=max(1.d-20,Dsoftshval(X,0.d0,0.,iscr))
          if(D.le.1.d-20)then
            if(ish.ge.1.and.iscr.ge.0)write(ifmt,*)
     &    "Warning in pompar ! Dsoftshval(1) could be negative"
c            sigma(i+1)=1.d5
          else
            npt=npt+1
            xlnXs(npt)=log(X)
            xlnD(npt)=log(D)
          endif
        enddo


c Fit of D(X) between X0 and xmax

        if(npt.gt.3)then
          call fit(xlnXs,xlnD,npt,sigma,1,a,beta)
          if(beta.gt.10.)beta=10.
          
          alpha=exp(min(50.,max(-50.,a)))
c        alpha=sngl(Dsoftshval(xmax,0.d0,0.,iscr))
c     &       *(sngl(xmax)*smaxDf)**(-beta)
        else  !very short range
          beta=10.
          alpha=0.
        endif


      elseif(iqq.eq.6.or.iqq.eq.7) then

        xmax=min(1d0,dble(xmxrem(1))*dble(xmxrem(2)))
c        print *,'pompar',xmax,sqrt(xminDf)
        X0=sqrt(dble(s0min))/dble(smaxDf)
c        if(X0.gt.xmax)xmax=1d0

c Calculation of X0 and D(X)

        if(ish.ge.4)write (ifch,'(1x,a,i2,a,1p,2e13.6)')
     *                    'pompar (',iqq,') x0,xmax=',X0,xmax

        npt=0
        do i=0,nptf-1
          X=X0
          if (i.ne.0) X=X*(xmax/X0)**(dble(i)/dble(nptf-1))
          if(iqq.eq.6)then    !fit x+
            Y=0.5d0*log(X)
          else                !fit x-
            Y=-0.5d0*log(X)
          endif
          D=max(1.d-20,Dsoftshval(X,Y,0.,iscr))
          if(D.le.1.d-20)then
c            write(ifch,*)
c     &    "Warning in pompar ! Dsoftshval(10) could be negative"
c            sigma(i+1)=1.d5
          else
            npt=npt+1
            xlnXs(npt)=log(X)
            xlnD(npt)=log(D)
          endif
cqq          write(ifch,*)i,X,D
        enddo


c Fit of D(X) between X0 and xmax
        if(npt.gt.3)then
          call fit(xlnXs,xlnD,npt,sigma,1,a,beta)
          if(beta.gt.10.)beta=10.

          alpha=exp(min(50.,max(-50.,a)))
c          alpha=sngl(Dsoftshval(xmax,0.d0,0.,iscr))
c     &       *(sngl(xmax)*smaxDf)**(-beta)
        else   !very short range
          beta=10.
          alpha=0.
        endif


      else                      !iqq=-4 or |iqq|=2 and xfitmin=1

c Calculation of X0 and D(X)

        xmax=1d0
        X0=max(2d0/smaxDf,0.1d0)
c        xmax=max(2.d0/dble(smaxDf),
c     &       min(max(0.03d0,sqrt(xminDf)),0.1d0))

        if(ish.ge.4)write (ifch,'(1x,a,i2,a,1p,2e13.6)')
     *                    'pompar (',iqq,') x0,xmax=',X0,xmax

        npt=0
        do i=0,nptf-1
          X=X0
          if (i.ne.0) X=X*(xmax/X0)**(dble(i)/dble(nptf-1))
          D=max(1.d-20,Dsoftshval(X,0.d0,0.,iscr))
          if(D.le.1.d-20)then
            if(ish.ge.1.and.iscr.ge.0)write(ifmt,*)
     &    "Warning in pompar ! Dsoftshval(-1) could be negative"
c            sigma(i+1)=1.d5
          else
            npt=npt+1
            xlnXs(npt)=log(X)
            xlnD(npt)=log(D)
          endif
        enddo

c Fit of D(X) between X0 and xmax

        if(npt.gt.3)then
          call fit(xlnXs,xlnD,npt,sigma,1,a,beta)
          if(beta.gt.10.)beta=10.


          alpha=exp(min(50.,max(-50.,a)))
c        alpha=sngl(Dsoftshval(xmax,0.d0,0.,iscr))
c     &           *(sngl(xmax)*smaxDf)**(-beta)
        else   !very short range
          beta=10.
          alpha=0.
        endif


      endif

      if(abs(alpha).lt.1e-10.or.abs(beta).gt.10)then       !bad fit
        alpha=0.
        beta=2.
c      else      !fix alpha at X=1 and not xmax
c        alpha=alpha/xmax**beta
      endif

        if(ish.ge.4)write(ifch,*) '%%%%%%%%%%%%% pompar %%%%%%%%%%%%%'
        if(ish.ge.4)write(ifch,*) 'alpD ini =',alpha,' betD ini=',beta

      return
      end

c----------------------------------------------------------------------

      double precision function droot(d0,d1,xmax,iscr)

c----------------------------------------------------------------------
c Find x0 which gives f(x0)=D(x0*S)-d0=0
c iqq=0 soft pomeron
c iqq=1 semi-hard pomeron
c----------------------------------------------------------------------
#include "aaa.h"
#include "sem.h"
#include "par.h"
      double precision Dsoftshval,d0,d1,d2,x0,x1,x2,f0,f1,f2,xmax
      parameter (kmax=1000)


      k=0
      x0=xggfit!max(xggfit,1d-5)
      x1=xmax
 5    d2=dabs(Dsoftshval(x0,0.d0,0.,iscr))
      if(d2.lt.1.d-10.and.x0.lt.x1)then
        x0=x0*1.1d0
c        print *,"droot 0",x0,x1,d0,d1,d2
        goto 5
      elseif(d2.gt.d0)then
        droot=x0
c        print *,"droot 1",x0,x1,d0,d1,d2
        return
      endif
      f0=d2-d0
      f1=d1-d0
      if(f0*f1.lt.0.d0)then


 10   x2=dsqrt(x0*x1)
      d2=dabs(Dsoftshval(x2,0.d0,0.,iscr))
      f2=d2-d0
      k=k+1
c        write (ifch,*) '******************* droot **************'
c        write (ifch,*) x0,x1,x2,f0,f1,f2

      if (f0*f2.lt.0.D0) then
        x1=x2
        f1=f2
      else
        x0=x2
        f0=f2
      endif

      if (dabs((x1-x0)/x1).gt.(1.D-5).and.k.le.kmax.and.x1.ne.x0) then
        goto 10
      else
        if (k.gt.kmax) then
          write(ifch,*)'??? Warning in Droot: Delta=',dabs((x1-x0)/x1)
c.........stop 'Error in Droot, too many steps'
        endif
        droot=dsqrt(x1*x0)
      endif

      else
        droot=dsqrt(x1*x0)
      endif

      return
      end

c----------------------------------------------------------------------

      double precision function drootom(d0,dmax,bmax,eps)

c----------------------------------------------------------------------
c Find b0 which gives f(b0)=(1-exp(-om(b0,iqq)))/dmax-d0=0
#include "aaa.h"
#include "sem.h"
#include "par.h"
      double precision om1intbc,d0,d1,d2,f0,f1,f2,dmax
      parameter (kmax=1000)


      k=0
      b0=0.
      b1=bmax
      d2=(1.d0-exp(-om1intbc(b1)))/dmax
      if(d2.gt.d0)then
        drootom=b1
c        write(*,*)"drootom exit (1)",b0,b1,d0,d1,d2
        return
      endif
      d1=(1.d0-exp(-om1intbc(b0)))/dmax
      f0=d1-d0
      f1=d2-d0
      if(f0*f1.lt.0.d0)then


 10   b2=0.5*(b0+b1)
      d2=(1.d0-dexp(-om1intbc(b2)))/dmax
      f2=d2-d0
      k=k+1
c      write (*,*) '******************* drootom **************'
c      write (*,*) b0,b1,b2,f0,f1,f2

      if (f1*f2.lt.0.D0) then
        b0=b2
        f0=f2
      else
        b1=b2
        f1=f2
      endif

      if (abs(f2).gt.eps.and.k.le.kmax.and.b1.ne.b0) then
        goto 10
      else
        if (k.gt.kmax) then
          write(ifch,*)'??? Warning in Drootom: Delta=',abs((b1-b0)/b1)
c.........stop 'Error in Droot, too many steps'
        endif
        drootom=0.5*(b1+b0)
      endif

      else
c        write(*,*)"drootom exit (2)",b0,b1,d0,d1,d2
        drootom=0.5*(b1+b0)
      endif

      return
      end

c----------------------------------------------------------------------
      subroutine variance(r2,alp,iqi)
c----------------------------------------------------------------------
c fit sigma2 into : 1/sigma2(x)=1/r2-alp*log(x*s)
c  iqi=0 -> soft pomeron
c  iqi=1 -> semi-hard pomeron
c  iqi=2 -> sum
c----------------------------------------------------------------------
#include "aaa.h"
#include "sem.h"
#include "par.h"
      double precision Xs(maxdataDf),vari(maxdataDf),sigma(maxdataDf)
      double precision X,X0,xmax

      do i=1,nptf
        sigma(i)=1.d-2
      enddo
      r2=0.

      iqq=iqi
      if(iqq.ne.4.or.xshmin.gt.0.95d0)then
        X0=dble(1d0/smaxDf)  
c        if(iomega.lt.2)X0=1d0/dble(smaxDf)
        xmax=sqrt(X0) !fit the entire range to avoid problem at large x
c        if(iomega.ge.2.and.iregge.ne.0)xmax=0.9d0/smaxDf
c        if(iqq.lt.0)then
c          X0=1d0/max(2d0,dble(smaxDf)/100.d0)
c          xmax=1d0
c       endif
        if(iqq.gt.5.and.iqq.lt.10)then
          xmax=X0
          X0=xminDf
          if(iqq.eq.8)iqq=3
        endif
      elseif(iqq.eq.4)then
        X0=xshmin
        xmax=xmaxDf
      else
        X0=1d0/log(max(exp(2.d0),dble(smaxDf)/1.d3))
        xmax=xmaxDf
      endif

      if(iqq.ge.0)then

        npt=0
        do i=0,nptf-1
          X=X0
          if (i.ne.0) X=X*(xmax/X0)**(dble(i)/dble(nptf-1))
          Xs(i+1)=log(X/XminDf)       !to get delD=sigma2 at XminDf
c          sig2=sigma1i(X)
          sig2=sigma2(X,iqq)
c          write(ifch,*)'sigma',iqq,X,sig2
c          if(sig2.le.0)then
c            sigma(i+1)=1d5
c            vari(i+1)=0d0
c          else
          if(sig2.gt.0)then
            npt=npt+1
            vari(npt)=1.d0/dble(sig2)
          endif
        enddo
 

        if(npt.gt.3)then
c Fit of the variance of D(X,b) between X0 and xmaxDf

        call fit(Xs,vari,npt,sigma,1,tr2,talp)
c talp should not be positive to avoid problem at large b : take value at x=1
c        talp=min(0d0,talp)

        r2=1./tr2
        alp=-talp
c in principle, the formula to convert 1/(del+eps*log(sy)) into
c  1/del*(1-eps/del*log(sy)) is valid only if eps/del*log(sy)=alp*r2*log(sy)
c is small. In practice, since the fit of G(x) being an approximation, each
c component of the fit should not be taken separatly but we should consider
c G as one function. Then it works even with larger alp (gamD). But if alp 
c is too large then even a subdominant contribution can become too large so 
c fit has to be done on the full range to stay close to the original value 
c even if the diagram is not a priori important there.
c        ttt=alp*r2*log(smaxDf)
c        if(ttt.gt.0.5)
c     &    write(ifmt,*)'Warning, G(b) parametrization not optimal : ',
c     &          'gamD too large compared to delD !',ttt

      elseif(npt.gt.0)then
        r2=0.
        xcnt=0.
        do i=1,npt
c          if(vari(i+1).gt.0.)then
            r2=r2+sngl(1d0/vari(i+1))
            xcnt=xcnt+1.
c          endif
        enddo
        r2=r2/xcnt
        alp=0.
      else
        r2=1.
        alp=0.
c        call utstop('In variance, initial(1) sigma2 not def !&')
      endif

      else
        if(iqq.eq.-3)r2=sigma2(xminDf,3)
        if(iqq.eq.-4)r2=sigma2(xmaxDf,3)
        if(r2.le.0.) call utstop
     &('In variance, initial(2) sigma2 not def!&')
        alp=0.
      endif

      if(ish.ge.4)then
        write(ifch,*) '%%%%%%%%%% variance ini %%%%%%%%%%%%'
        write(ifch,*) 'X0=',X0
        write(ifch,*) 'delD ini=',r2
        write(ifch,*) 'gamD ini=',alp
      endif

      return
      end



c----------------------------------------------------------------------

        function sigma2(x,iqq)

c----------------------------------------------------------------------
c Return the variance for a given x of :
c For G :
c iqq=0 the soft pomeron
c iqq=1 the semi-hard and valence quark pomeron
c iqq=2 the sum
c----------------------------------------------------------------------

#include "aaa.h"
#include "par.h"
      double precision x,y,Dsoftshval,sfsh,eps,range,sig2!,omNpuncut
      external varifit
      double precision varifit,Db(maxdataDf),bf(maxdataDf)

      bmax=bmaxDf
      sig2=bmax
      bmin=0
      eps=1.d-30
      ierro=0

      if(iqq.le.0)then
        iscr=0
        range=sig2
        sfsh=Dsoftshval(x,0.d0,0.,iscr)
c        if(iqq.lt.0)then
c         bmin=0.5*bmax
c         sfsh=0.5d0*(sfsh+dble(parDf(1,3))*x**parDf(4,3))
c        endif
        if (dabs(sfsh).gt.eps) then
        do i=0,nptf-1
          bf(i+1)=dble(bmin+float(i)*(bmax-bmin)/float(nptf-1))
          Db(i+1)=Dsoftshval(x,0.d0,sngl(bf(i+1)),iscr)/sfsh
        enddo
        else
          ierro=1
        endif
      elseif(iqq.eq.1.and.xshmin.lt..95d0)then
        iscr=4
        range=sig2
        sfsh=Dsoftshval(x,0.d0,0.,iscr)
        if (dabs(sfsh).gt.eps) then
        do i=0,nptf-1
          bf(i+1)=dble(bmin+float(i)*(bmax-bmin)/float(nptf-1))
          Db(i+1)=Dsoftshval(x,0.d0,sngl(bf(i+1)),iscr)
          Db(i+1)=Db(i+1)/sfsh
        enddo
        else
          ierro=1
        endif
      else
        sig2=2.d0*sig2
        range=sig2
        iscr=iqq
        if(iqq.eq.6)then        !fit x+
          y=0.5d0*log(x)
        elseif(iqq.eq.7)then                    !fit x-
          y=-0.5d0*log(x)
        else
          y=0d0
        endif
        sfsh=Dsoftshval(x,y,0.,iscr)
        if (dabs(sfsh).gt.eps) then
        do i=0,nptf-1
          bf(i+1)=dble(bmin+float(i)*(bmax-bmin)/float(nptf-1))
          Db(i+1)=Dsoftshval(x,y,sngl(bf(i+1)),iscr)
     &              /sfsh
        enddo
        else
          ierro=1
        endif
      endif

c Fit of D(X,b) between -bmaxDf and bmaxDf

      if(ierro.ne.1)then
        nptft=nptf
        call minfit(varifit,bf,Db,nptft,sig2,range)
        sigma2=sngl(sig2)
        if(sigma2.gt.bmax)sigma2=0.
      else
        sigma2=0.
      endif
c      write(ifch,*)'sigma2',iqq,iscr,sfsh,x,sigma2
          
      return
      end

c----------------------------------------------------------------------

        subroutine paramx

c----------------------------------------------------------------------
c updates the 4 parameters alpsf,betsf,alpsh,betsh by fitting GFF
c  parDf(1,3) parDf(4,3) ... alp, bet soft
c  parDf(2,3) parDf(5,3) ... alp, bet sat
c  parDf(3,3) parDf(6,3) ... alp, bet semihard
c----------------------------------------------------------------------

#include "aaa.h"
#include "par.h"

      double precision Dsoftshpar,chi2,error

      external Dsoftshpar

      dimension range(nbpf)

      call givedatax

      !determine parameter range
      do i=numminDf,numparDf
        range(i)=fparDf*parDf(i,3)
        parDf(i,1)=parDf(i,3)-range(i)
        parDf(i,2)=parDf(i,3)+range(i)
        if(i.eq.4)parDf(i,1)=-0.99
      enddo


   !   write(ifch,*) '%%%%%%%%%%%%%%%%%%% fitx %%%%%%%%%%%%%%%%%%%%%%%'

      call fitx(Dsoftshpar,nmcxDf,chi2,error)

   !   write(ifch,*) 'chi2=',chi2
   !   write(ifch,*) 'err=',error
   !   write(ifch,*) 'alpD(1)=',parDf(1,3),' betD(1)=',parDf(4,3)
   !   write(ifch,*) 'alpD(2)=',parDf(2,3),' betD(2)=',parDf(5,3)
   !   write(ifch,*) 'alpD(3)=',parDf(3,3),' betD(3)=',parDf(4,3)

      return
      end


c----------------------------------------------------------------------
      subroutine givedatax
c----------------------------------------------------------------------

#include "aaa.h"
#include "sem.h"
#include "par.h"
      double precision X,X0,X1,Dsoftshval,Xseuil

      numdataDf=nptf

      X0=xminDf
      X1=xmaxDf !min(100d0*X0,xmaxDf)
      Xseuil=1d0 !min(1.d0,xggfit*1d4)
c      print *,'--------------->',Xseuil

c Fit of G(X) between X0 and X1
      do i=0,nptf-1
        X=X0
        if (i.ne.0) X=X*(X1/X0)**(dble(i)/dble(nptf-1))
        datafitD(i+1,2)=max(1.d-10,
     &                  Dsoftshval(X,0.d0,0.,0))
        datafitD(i+1,1)=X
        datafitD(i+1,3)=1.
        if (X.gt.Xseuil)
     &  datafitD(i+1,3)=exp(-min((Xseuil/X-1.),50.d0))
      enddo

      return
      end





c----------------------------------------------------------------------

      function sigma1i(x)

c----------------------------------------------------------------------
c Return the variance of the sum of the soft pomeron and the semi-hard
c pomeron for a given x.
c----------------------------------------------------------------------

#include "aaa.h"
#include "par.h"
      double precision x,Dsoftshval,Dint


      iscr=0
      Dint=Dsoftshval(x,0.d0,0.,iscr)

      sigma1i=0.
      if(Dint.ne.0.)
     &sigma1i=sngl(-1.d0/log(Dsoftshval(x,0.d0,1.,iscr)
     &   /Dint))


      return
      end


c----------------------------------------------------------------------

      SUBROUTINE minfit(func,x,y,ndata,a,range)

c----------------------------------------------------------------------
c Given a set of data points x(1:ndata),y(1:ndata), and the range of
c the parameter a, fit it on function func by minimizing chi2.
c In input a define the expected value of a, and on output they
c correspond to  the fited value.
c ---------------------------------------------------------------------
#include "aaa.h"
      double precision x(ndata),y(ndata),func,a,range,Smin,Som,a1,a2,eps
     *,amin,rr,yp,drangen
      parameter (eps=1.d-5)
      external func


      Smin=1.d20
      amin=a



      a1=a-range
      a2=a+range

      do j=1,2000
        rr=drangen(amin)
        a=a1+(a2-a1)*rr
        k=0

 10     if(a.lt.0.d0.and.k.lt.100) then
          rr=dble(rangen())
          a=a1+(a2-a1)*rr
          k=k+1
          goto 10
        endif
        if(k.ge.100) call utstop
     &('Always negative variance in minfit ...&')

        Som=0.d0
        do k=1,ndata
             yp=min(1.d10,func(x(k),a))  !proposal function
              Som=Som+(yp-y(k))**2.d0
        enddo
        if(Som.lt.Smin)then

          if(Smin.lt.1.)then
            if(a.gt.amin)then
              a1=amin
            else
              a2=amin
            endif
          endif
          amin=a
          Smin=Som
        endif
        if(Smin.lt.eps)goto 20
      enddo

 20   continue
      a=amin

      return
      end



c----------------------------------------------------------------------
      subroutine fitx(func,nmc,chi2,err)
c----------------------------------------------------------------------
c  Determines parameters of the funcion func
c  representing the best fit of the data.
c  At the end of the run, the "best" parameters are stored in parDf(n,3).
c  The function func has to be defined via "function" using the parameters
c  parDf(n,3), n=1,numparDf .
c  Parameters as well as data are stored on /fitpar/:
c    numparDf: number of parameters  (input)
c    parDf: array containing parameters:
c         parDf(n,1): lower limit    (input)
c         parDf(n,2): upper limit    (input)
c         parDf(n,3): current parameter (internal and output = final result)
c         parDf(n,4): previous parameter (internal)
c    numdataDf: number of data points  (input)
c    datafitD: array containing data:
c         datafitD(i,1): x value       (input)
c         datafitD(i,2): y value       (input)
c         datafitD(i,3): error         (input)
c----------------------------------------------------------------------

#include "aaa.h"
#include "par.h"
      double precision func,x,chi2,err,drangen
      external func

 !     write (ifch,*) 'numparDf,numminDf',numparDf,numminDf


c initial configuration (better if one start directly with initial one)

c      do n=numminDf,numparDf
c              parDf(n,3)=parDf(n,1)+rangen()*(parDf(n,2)-parDf(n,1))
c      enddo

      chi2=0.
      err=0.
      do i=1,numdataDf
        x=datafitD(i,1)
        fx=func(x)
        chi2=chi2+(log(fx)-log(datafitD(i,2)))**2/datafitD(i,3)**2
        err=err+(log(fx)-log(datafitD(i,2)))/datafitD(i,3)**2
      enddo
      err=abs(err)/dble(numdataDf)

c metropolis iteration

      do i=1,nmc
c        if(mod(i,int(float(nmc)/1000.)).eq.0)then
          betac=betac*(1.+1./float(nmc))!1.05
          betae=betae*(1.+1./float(nmc))!1.05
c        endif
c        if(mod(i,int(float(nmc)/20.)).eq.0)write(ifch,*)i,chi2,err

        do n=numminDf,numparDf
          parDf(n,4)=parDf(n,3)
        enddo
        chi2x=chi2
        errx=err

        n=numminDf+int(rangen()*(numparDf-numminDf+1))
        n=max(n,numminDf)
        n=min(n,numparDf)
c              if(mod(i,int(float(nmc)/20.)).eq.0)write(ifch,*)n

        parDf(n,3)=parDf(n,1)+rangen()*(parDf(n,2)-parDf(n,1))

        chi2=0
        err=0
        do j=1,numdataDf
          x=datafitD(j,1)
          fx=func(x)
          chi2=chi2+(log(fx)-log(datafitD(j,2)))**2/datafitD(j,3)**2
          err=err+(log(fx)-log(datafitD(j,2)))/datafitD(j,3)**2
        enddo
        err=abs(err)/dble(numdataDf)

        if(chi2.gt.chi2x.and.drangen(chi2)
     $     .gt.exp(-min(50.d0,max(-50.d0,dble(betac)*(chi2-chi2x))))
     &     .or.err.gt.errx.and.drangen(err)
     $     .gt.exp(-min(50.d0,max(-50.d0,dble(betae)*(err-errx))))
     &                                                     ) then
          do n=numminDf,numparDf
            parDf(n,3)=parDf(n,4)
          enddo
          chi2=chi2x
          err=errx
        endif

      enddo

      return
      end


c----------------------------------------------------------------------

        SUBROUTINE fit(x,y,ndata,sig,mwt,a,b)

c----------------------------------------------------------------------
c Given a set of data points x(1:ndata),y(1:ndata) with individual standard
c deviations sig(1:ndata), fit them to a straight line y = a + bx by
c minimizing chi2 .
c Returned are a,b and their respective probable uncertainties siga and sigb,
c the chi-square chi2, and the goodness-of-fit probability q (that the fit
c would have chi2 this large or larger). If mwt=0 on input, then the standard
c deviations are assumed to be unavailable: q is returned as 1.0 and the
c normalization of chi2 is to unit standard deviation on all points.
c ---------------------------------------------------------------------

        implicit none
        INTEGER mwt,ndata
        double precision sig(ndata),x(ndata),y(ndata),chi2
        REAL a,b,siga,sigb !,q
        INTEGER i
        double precision sigdat,ss,st2,sx,sxoss,sy,t,wt


        sx=0.d0                 !Initialize sums to zero.
        sy=0.d0
        st2=0.d0
        b=0.
        a=0.
        if(mwt.ne.0) then ! Accumulate sums ...
          ss=0.d0
          do i=1,ndata          !...with weights
            wt=1.d0/(sig(i)**2)
            ss=ss+wt
            sx=sx+x(i)*wt
            sy=sy+y(i)*wt
          enddo
        else
          do i=1,ndata          !...or without weights.
            sx=sx+x(i)
            sy=sy+y(i)
          enddo
          ss=dble(ndata)
        endif
        sxoss=sx/ss
        if(mwt.ne.0) then
          do i=1,ndata
            t=(x(i)-sxoss)/sig(i)
            st2=st2+t*t
            b=b+t*y(i)/sig(i)
          enddo
        else
          do i=1,ndata
            t=x(i)-sxoss
            st2=st2+t*t
            b=b+t*y(i)
          enddo
        endif
        b=b/sngl(st2)                 !Solve for a, b, oe a , and oe b .
        a=sngl((sy-sx*b)/ss)
        siga=sngl(sqrt((1.+sx*sx/(ss*st2))/ss))
        sigb=sngl(sqrt(1./st2))
        chi2=0.d0                 !Calculate chi2 .
c        q=1.
        if(mwt.eq.0) then
          do i=1,ndata
            chi2=chi2+(y(i)-dble(a)-dble(b)*x(i))**2
          enddo

c For unweighted data evaluate typical sig using chi2, and adjust
c the standard deviations.

          sigdat=sqrt(chi2/(ndata-2))
          siga=siga*sngl(sigdat)
          sigb=sigb*sngl(sigdat)
        else
          do i=1,ndata
            chi2=chi2+((y(i)-dble(a)-dble(b)*x(i))/sig(i))**2
          enddo
        endif

c        if(chi2.ge.0.2)then
c          b=(y(ndata)-y(1))/(x(ndata)-x(1))
c          a=y(ndata)-b*x(ndata)
c        endif


c        write(31,*) x,y
c        write(31,*) '$$$$$$$$ fit : a,b,siga,sigb,chi2,q $$$$$$$$$'
c        write(31,*) a,b,siga,sigb,chi2!???????????????


        return
        END



c----------------------------------------------------------------------

      double precision function varifit(x,var)

c----------------------------------------------------------------------
      double precision x,var

      varifit=dexp(-min(50.d0,x**2.d0/var))

      return
      end



c----------------------------------------------------------------------

      double precision function Dsoftshval(x,y,b,iscr)

c----------------------------------------------------------------------
c iscr=0 soft
c iscr=1 sum of om5p (i), i from 0 to 4 * F * F
c iscr=2 sum of om5p (i), i from 1 to 4 (semihard + valence quark)
c iscr=3 sum of om5p (i), i for soft, Reggeon, DD and DPE only
c iscr=4 soft + hard
c iscr>4 Reggeon or individual diffractive contribution rescaled to take
c        into account the factor x**alpdif introduced to have bet<-1
c iscr=-1 sum of om5p (i), i from 0 to 9 (with regge and Diff) minus fit of SD
c iscr=-3 same as iscr=3 minus fit of hard
c iscr=-4 same as iscr=4 minus fit of SD and hard
c----------------------------------------------------------------------
      double precision x,om51,y,xp,xm
#include "aaa.h"
#include "sem.h"
#include "par.h"

      Dsoftshval=0.d0
        
      if(abs(iscr).eq.3)then    !G0 ((soft+)DD+DPE+R)

        Dsoftshval=om51(x,y,b,0,0)
        if(iregge.ne.0)Dsoftshval=Dsoftshval+om51(x,y,b,5,5)       !Reggeon
        if(iomega.lt.2)then                           !DD+DPE
            Dsoftshval=Dsoftshval
     .                +om51(x,y,b,-8,-9)
        endif
c        if(Dsoftshval.le.0d0)Dsoftshval=om51(x,y,b,0,0)

      elseif(iscr.eq.1)then
        xp=dsqrt(x)*dexp(y)
        if(dabs(xp).ge.1.d-15)then
          xm=x/xp
        else
          xm=1.d0
          write(ifch,*)'Warning in Dsoftshval in epos-par'
        endif
        Dsoftshval=om51(x,y,b,0,4)+om51(x,y,b,11,11)
        Dsoftshval=Dsoftshval*(1.d0-xm)**dble(alplea(icltar))
     &                       *(1.d0-xp)**dble(alplea(iclpro))
      elseif(iscr.eq.2)then
c hard component only (used to test if purely soft)
        Dsoftshval=om51(x,y,b,1,4)
      elseif(abs(iscr).eq.4)then
c soft + hard component (used for G1)
        Dsoftshval=om51(x,y,b,0,4)+om51(x,y,b,11,11)
            
      elseif(iscr.eq.10)then     !soft
        
        Dsoftshval=om51(x,y,b,0,0)+om51(x,y,b,11,11)
      
      elseif(iscr.le.0)then     !all
        
        Dsoftshval=om51(x,y,b,0,-5)+om51(x,y,b,11,11)
      
      else
        xp=dsqrt(x)*dexp(y)
        if(dabs(xp).ge.1.d-15)then
          xm=x/xp
        else
          xm=1.d0
          write(ifch,*)'Warning in Dsoftshval(2) in epos-par'
        endif
c individual diff components
        if(iscr.ge.6.and.iscr.le.9)then   !diff for MC
          iqq=-iscr
        else
          iqq=iscr
        endif
        Dsoftshval=om51(x,y,b,iqq,iqq)
      endif


c      Dsoftshval=2.d0*Dsoftshval
      Dsoftshval=Dsoftshval
     &     /(x**dble(-alppar))

      if(iscr.lt.0)then
        Dsoftshval=Dsoftshval
     &             -dble(parDf(2,3))*x**parDf(5,3)
        if(iscr.eq.-3)Dsoftshval=Dsoftshval
     &             -dble(parDf(3,3))*x**parDf(6,3)
      endif

c      write(ifch,*)'Dsoftshval',iscr,chadr,b,y,x,Dsoftshval

      return
      end

c----------------------------------------------------------------------

      double precision function Dsoftshpar(x)

c----------------------------------------------------------------------
      double precision x
#include "aaa.h"
#include "par.h"

      Dsoftshpar=
     &        dble(parDf(1,3))*x**parDf(4,3)
     &       +dble(parDf(2,3))*x**parDf(5,3)
     &       +dble(parDf(3,3))*x**parDf(6,3)
      Dsoftshpar=min(max(1.d-15,Dsoftshpar),1.d15)

      return
      end
