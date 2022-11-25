C
C  This file is part of EPOS4
C  Copyright (C) 2022 research institutions and authors (See CREDITS file)
C  This file is distributed under the terms of the GNU General Public License version 3 or later
C  (See COPYING file for the text of the licence)
C

!###############################################################################
!###############################################################################
!           cross section and CR related stuff
!###############################################################################
!###############################################################################

!
!
!       xsectionpar()  should be updated, when necessary !!!!!!
!
!


c-------------------------------------------------------------------------------
      subroutine epocrossc(niter,gtot,gprod,gabs,gcoh,gqel,gdd)
c-------------------------------------------------------------------------------
c epocrossc - nucleus-nucleus (nucleus-hydrogen) interaction cross sections
c by calculation with real nuclear profiles and eikonal (simplified simulations)
c gtot  - total cross section
c gprod - production cross section (all diffraction included)
c gabs  - cut Pomerons cross section (no diffraction at all)
c gdd   - proj (ionudi=2) or proj or targ (ionudi=0/3) excited diffraction
c         cross section
c gcoh  - coherent (elastic with respect to the projectile) cross section
c      (non excited diff proj if ionudi=2, non excited proj+targ if ionudi=0/3)
c
c Be careful : this function is not symmetric for gdd and gqel (only projectile
c diffraction) in case of ionudi=2.
c (target diffraction is not treated explicitely and contributes to
c gprod, gdd, gcoh and gtot).
c
c WARNING : results are sure only in case of ionudi=1 (no substraction from
c           diffractive part) in particular for AA with A > 10 (nuclear diff
c           not well described). For pA seems to be OK with ionudi 2 and 3.
c
c code from QGSJET programs by S.Ostapchenko
c-------------------------------------------------------------------------------
#include "aaa.h"
#include "ems.h"
      common /cncl/xproj(mamx),yproj(mamx),zproj(mamx)
     *            ,xtarg(mamx),ytarg(mamx),ztarg(mamx)
      common/geom/rmproj,rmtarg,bmax,bkmx
      dimension wabs(28),wdd(28),wcoh(28),wprod(28),wqel(28)
     &         ,b0(28),ai(28)
      common /ar3/ x1(7),a1(7)
      double precision xgabs,xgdd,xgprod,xgcoh,xgqel

      call utpri('epocrs',ish,ishini,2)
      if(ish.ge.2)write(ifch,201)niter,bmax
      kollini=koll        !koll modified
      do i=1,7
       b0(15-i)=bmax*sqrt((1.+x1(i))/2.)
       b0(i)=bmax*sqrt((1.-x1(i))/2.)
       ai(i)=a1(i)*bmax**2*pi*5.05        !factor change cs
       ai(15-i)=ai(i)
      enddo
      if(maproj.gt.1.or.matarg.gt.1)then
        difn=max(difnuc(maproj),difnuc(matarg))
      else
        difn=1.
      endif
      do i=1,7
        tp=(1.+x1(i))/2.
        tm=(1.-x1(i))/2.
        b0(14+i)=bmax-log(tp)*difn
        b0(29-i)=bmax-log(tm)*difn
        ai(14+i)=a1(i)*b0(14+i)/tp*10.*difn*pi
        ai(29-i)=a1(i)*b0(29-i)/tm*10.*difn*pi
      enddo
      do i=1,28
       wabs(i)=0.
       wdd(i)=0.
       wprod(i)=0.
       wcoh(i)=0.
       wqel(i)=0.
      enddo
      do nc=1,niter
        if(maproj.eq.1)then
          xproj(1)=0.
          yproj(1)=0.
          zproj(1)=0.
        else
          call conxyz('p',mamx,xproj,yproj,zproj,ypjtl-yhaha)
        endif
        if(matarg.eq.1)then
          xtarg(1)=0.
          ytarg(1)=0.
          ztarg(1)=0.
        else
          call conxyz('t',mamx,xtarg,ytarg,ztarg,yhaha)
        endif

        do i=1,28
          call epogcr(b0(i),xgabs,xgdd,xgprod,xgcoh,xgqel)
          wabs(i)=wabs(i)+sngl(xgabs)
          wdd(i)=wdd(i)+sngl(xgdd)
          wprod(i)=wprod(i)+sngl(xgprod)
          wcoh(i)=wcoh(i)+sngl(xgcoh)
          wqel(i)=wqel(i)+sngl(xgqel)
        enddo
      enddo

      gabs=0.
      gdd=0.
      gcoh=0.
      gprod=0.
      gqel=0.
      do i=1,28
       wabs(i)=wabs(i)/niter
       wdd(i)=wdd(i)/niter
       wcoh(i)=wcoh(i)/niter
       wprod(i)=wprod(i)/niter
       wqel(i)=wqel(i)/niter
       gabs=gabs+ai(i)*wabs(i)
       gdd=gdd+ai(i)*wdd(i)
       gcoh=gcoh+ai(i)*wcoh(i)
       gqel=gqel+ai(i)*wqel(i)
       gprod=gprod+ai(i)*wprod(i)
      enddo


      gtot=gprod+gcoh            !total=all cut (with diff) + all uncut
      if(ish.ge.2)then
        if(ionudi.eq.2)then
          write (ifch,202)gtot,gprod,gabs,gdd,gcoh,gqel
        else
          write (ifch,203)gtot,gprod,gabs,gdd,gcoh,gqel
        endif
      endif
        
201   format(2x,'epocrossc - A-B interaction cross sections,'
     *,' N of iter.:',i5,' bmax:',f5.2)
202   format(2x,'epocrossc: gtot=',e10.3,2x,'gprod=',e10.3,2x
     *,'gabs=',e10.3/4x,'gpe=',e10.3,2x,'gcoh=',e10.3,'gqel=',e10.3)
 203  format(2x,'epocrossc: gtot=',e10.3,2x,'gprod=',e10.3,2x
     *,'gabs=',e10.3/4x,'gdd=',e10.3,2x,'gcoh=',e10.3,'gsd=',e10.3)


      koll=kollini
      call utprix('epocrs',ish,ishini,2)

      return
      end

c-------------------------------------------------------------------------------
      subroutine epogcr(b,gabs,gdd,gprod,gcoh,gqel)
c-------------------------------------------------------------------------------
c epogcr - integrands (b-profiles) for nucleus-nucleus cross sections
c b - impact parameter
c code from QGSJET programs by S.Ostapchenko
c-------------------------------------------------------------------------------
#include "aaa.h"
#include "ems.h"
#include "par.h"
      common /cncl/xproj(mamx),yproj(mamx),zproj(mamx)
     *            ,xtarg(mamx),ytarg(mamx),ztarg(mamx)
      common/geom/rmproj,rmtarg,bmax,bkmx
      double precision vin,gabs,gdd,gprod,gcoh,fdd,gqel,fdt,vdt,vcu,vdp

      if(ish.ge.9)write (ifch,201)b
      gprod=1d0
      gabs=1d0
      gdd=1d0
      fdd=1d0
      fdt=1d0
      bx=0
      by=0

      if(maproj.eq.1.and.matarg.eq.1)then
        if(b.gt.bkmx)then
          koll=0
        else
          koll=1
          bk(1)=b
          iproj(1)=1
          itarg(1)=1
          lproj(1)=1
          ltarg(1)=1
          lproj3(1)=1
          ltarg3(1)=1
          kproj3(1,1)=1
          ktarg3(1,1)=1
          kproj(1,1)=1
          ktarg(1,1)=1
        endif
      else
        bx=b
        by=0.
        koll=0
        do i=1,maproj
          lproj(i)=0
          lproj3(i)=0
        enddo
        do j=1,matarg
          ltarg(j)=0
          ltarg3(j)=0
        enddo
        do 12 i=1,maproj
        do 11 j=1,matarg
          bij=sqrt((xproj(i)+bx-xtarg(j))**2+(yproj(i)+by-ytarg(j))**2)
          if(bij.gt.bkmx)goto 11

          koll=koll+1
          if(koll.gt.kollmx)call utstop('epogcr: kollmx too small&')
          bk(koll)=bij
          bkx(koll)=xproj(i)+bx-xtarg(j)
          bky(koll)=yproj(i)+by-ytarg(j)
          iproj(koll)=i
          itarg(koll)=j
          lproj(i)=lproj(i)+1
          ltarg(j)=ltarg(j)+1
          kproj(i,lproj(i))=koll
          ktarg(j,ltarg(j))=koll

 11     continue
 12     continue
      endif
      if(koll.eq.0)then
        gabs=0d0
        gdd=0d0
        gprod=0d0
        gcoh=0d0
        gqel=0d0
        goto 1000
      endif
c      if(iscreen.ne.0)call CalcScrPair(b)

      irea=-1
      call GfunParK(irea)
      if(ionudi.eq.0
     &  .and.(maproj.ne.1.or.matarg.ne.1).and.nglevt.eq.0)then
        gabs=0d0
        gdd=0d0
        gprod=0d0
        gcoh=0d0
        gqel=0d0
        goto 1000
      endif
      do n=1,maproj
       call epov(n,vin,vcu,vdp,vdt)
       gprod=gprod*vin
       gdd=gdd*vcu
       fdd=fdd*vdp
       fdt=fdt*vdt
      enddo
      gprod=min(gprod,1.d0)
      gcoh=1d0-2d0*sqrt(gprod)+gprod
      gprod=1d0-gprod
      gdd=max(0d0,gdd)       !diffractive part
      gabs=max(0d0,gprod-gdd)          !cut (no diffraction)
      if(ionudi.eq.2.and.maproj+matarg.gt.2)then
        gqel=fdd      !quasielastic = diffractive without excited proj.
        gdd=max(0d0,gdd-fdd)             !only excited projectile diffraction
      else
        gqel=fdd+fdt  !SD = diffractive without excited proj. or targ (if a proj is not excited then the target has to be. Can not compute low mass case where excited nucleon simply break the nucleus but don't produce new particles).
        gdd=max(0d0,gdd-gqel)  !DD = only diffraction but on both side
      endif
 1000 continue
      if(ish.ge.9)write (ifch,202)gabs,gdd,gprod,gcoh,gqel,fdd,fdt

201   format(2x,'epogcr-integrands for nucleus-nucleus cross sections,'
     *,' b=',e10.3)
202   format(2x,'epogcr: gabs=',e10.3,2x,'gdd=',e10.3,2x,'gprod=',e10.3
     *,2x,'gcoh=',e10.3,2x,'gqel=',e10.3,2x,'fdd=',e10.3,' fdt=',e10.3)
      return
      end

c=============================================================================
      subroutine epov(n,vin,vdi,vdp,vdt)
c epov - eikonal factors for nucleus-nucleus interaction
c (used for cross-section calculation)
c n - projectile nucleon indice
c vin - all uncut pomerons
c vdi - only diff pomerons
c vdp - single diffractive without excitation for projectile (target is excited)
c vdt - single diffractive without excitation for target (projectile is excited)
c nkol - number of collision
c code from QGSJET programs by S.Ostapchenko
c----------------------------------------------------------------------------
#include "aaa.h"
#include "ems.h"
      common /cncl/xproj(mamx),yproj(mamx),zproj(mamx)
     *            ,xtarg(mamx),ytarg(mamx),ztarg(mamx)
      double precision vvv2,vvv1,dvp,dvt,vin,vdi,vdt,vdp,PhiExpoK
     *                ,om1intby

      if(ish.ge.9)write (ifch,201)xproj(n),yproj(n)

      imin=idxDmin(iomega)
      imax=idxDmax(iomega)
      spp=engy**2
      vin=0.d0
      vdi=0.d0
      vvv1=1.d0
      vvv2=1.d0
      dvp=1.d0
      dvt=1.d0
      do m=1,lproj(n)
        k=kproj(n,m)
        vvv1=vvv1*max(0.d0,PhiExpoK(k,1.d0,1.d0))  !no inelastic diagram
        zzp=zparpnx(1,k)
        zzt=zpartnx(1,k)
        do i=imin,imax          !initialization of parameters for b1
          call Gfunpar(zzp,zzt,0.,0.,1,i,bk(k),sy,a,b,c,d,e,f,g)
          call Gfunpar(zzp,zzt,0.,0.,2,i,bk(k),sy,a,b,c,d,e,f,g)
        enddo
        vvv2=vvv2*max(0.d0,om1intby(spp,bk(k),-50,1))         !only diff diagram (mass)
        dvp=dvp*om1intby(spp,bk(k),-60,1)     !probability for projectile not to be excited (exchange of only one SD diagram with target excitation)
        dvt=dvt*om1intby(spp,bk(k),-70,1)     !target is not excited
      enddo
      vdi=vvv2
      vin=vvv1                    !exp(-2 * chi)
      vdp=dvp
      vdt=dvt

      if(ish.ge.9)write (ifch,202)vin,vdi,vdp,vdt
      if(ish.ge.9)write (ifch,203)

201   format(2x,'epov - eikonal factor: nucleon coordinates x=',
     *e10.3,2x,'y=',e10.3)
202   format(2x,'vin=',e10.3,2x,'vdi=',e10.3,2x,'vdp=',e10.3
     *,'vdt=',e10.3)
203   format(2x,'epov - end')
      return
      end

c------------------------------------------------------------------------
      subroutine psfz(iqq,gz2,b)
c-----------------------------------------------------------------------
c hadron-nucleus cross sections calculation
c b - impact parameter squared
C iqq - 1 = elastic cross section
C       2 = inelastic cross section
c  Not yet updated for Q2smin !!!!!!!!!!!!!!!!!!!!!!????????????????????
c-----------------------------------------------------------------------
      double precision PhiExpo
#include "aaa.h"
#include "ems.h"
#include "par.h"
      common /ar3/ x1(7),a1(7)
      external pttcs,pprcs

      gz2=0.
      e1=exp(-1.)
      rs=r2had(iclpro)+r2had(icltar)+slopom*log(engy**2)
      rpom=4.*.0389*rs
      rpom=rpom*facmc
      kollini=koll        !koll modified in zzfz
      koll=1
      
      if(iscreen.ne.0.and.(maproj.gt.1.or.matarg.gt.1))then
        call zzfz(zzp,zzt,kollth,b)
        koll=kollth
      else
        zzp=0.
        zzt=0.
      endif
      sy=engy**2

      do i1=1,7
      do m=1,2
        z=.5+x1(i1)*(m-1.5)
        zv1=exp(-z)
        zv2=(e1*z)
        b1=sqrt(-rpom*log(zv1))
        b2=sqrt(-rpom*log(zv2))

        if(maproj.eq.1.and.matarg.eq.1)then
          cg1=1.
          cg2=1.
        elseif(matarg.eq.1)then
          cg1=ptrot(pprcs,b,b1)
          cg2=ptrot(pprcs,b,b2)
        else
          cg1=ptrot(pttcs,b,b1)
          cg2=ptrot(pttcs,b,b2)
        endif

        vv21=sngl(Phiexpo(zzp,zzt,1.,1.d0,1.d0,sy,b1))
        vv22=sngl(Phiexpo(zzp,zzt,1.,1.d0,1.d0,sy,b2))
        if(iqq.ne.1)then
          gz2=gz2+a1(i1)*(cg1*(1.-vv21)+cg2*(1.-vv22)/z)
        else
          vv11=sngl(Phiexpo(zzp,zzt,0.5,1.d0,1.d0,sy,b1))
          vv12=sngl(Phiexpo(zzp,zzt,0.5,1.d0,1.d0,sy,b2))
          gz2=gz2+a1(i1)*(cg1*(vv21-2.*vv11+1.)
     &                   +cg2*(vv22-2.*vv12+1.)/z)
        endif
      enddo
      enddo
      gz2=gz2*rpom*0.5

      koll=kollini

      return
      end

c------------------------------------------------------------------------
      subroutine zzfz(zzp,zzt,kollth,b)
c-----------------------------------------------------------------------
c hadron-nucleus cross sections calculation
c b - impact parameter squared
C xsfct - 0.5 = total cross section
C         1.0 = inelastic cross section
c-----------------------------------------------------------------------
      common /psar50/ zznuc,b2xnuc
#include "aaa.h"
#include "ems.h"
#include "par.h"
      common /ar3/ x1(7),a1(7)
      external  pttcs,pprcs,pttzz,pprzz

      zzp=0.
      zzt=0.
      kollth=1
      if(iscreen.eq.0.or.(maproj.eq.1.and.matarg.eq.1))return

      rs=r2had(iclpro)+r2had(icltar)+slopom*log(engy**2)
      rpom=4.*.0389*rs
      bgl2=2.*rpom*epscrp
      zzpp=min(abs(epscrx),fscra(engy))
c caculate the radius where Z is saturated at epscrx to define the bases
c of nuclear shadowing
      satrad=0.
      if(zzpp.gt.0.)satrad=-bgl2*log(abs(epscrx)/zzpp)
      bglx=sqrt(max(0.1,satrad))
      fzbrmax=0.
      if(zbrmax.gt.0)fzbrmax=zbrmax
      fzbcut=1.
      fzbrads=1.
      if(bglx.gt.0)fzbrads=bglx
      fnuc=1.2*fzbcut/fzbrads
      b2xnuc=bgl2+4.*fzbrmax*sqrt(float(maproj*matarg))*fnuc
      rpom=rpom*facmc

      e1=exp(-1.)

      colp=0.
      colt=0.
      do i1=1,7
      do m=1,2
        z=.5+x1(i1)*(m-1.5)
        zv1=exp(-z)
        zv2=(e1*z)
        b1=sqrt(-rpom*log(zv1))
        b2=sqrt(-rpom*log(zv2))


        if(maproj.gt.1)then
          cg1=ptrot(pprcs,b,b1)
          cg2=ptrot(pprcs,b,b2)
          colnuc=a1(i1)*(cg1+cg2/z)
          colp=colp+colnuc
          rho=0.05
          zznuc=fscro(engy/egyscr,rho)
          zp1=ptrot(pprzz,b,b1)
          zp2=ptrot(pprzz,b,b2)
          zzp=zzp+a1(i1)*(zp1+zp2/z)
        endif
        if(matarg.gt.1)then
          cg1=ptrot(pttcs,b,b1)
          cg2=ptrot(pttcs,b,b2)
          colnuc=a1(i1)*(cg1+cg2/z)
          colt=colt+colnuc
          rho=0.05
          zznuc=fscro(engy/egyscr,rho)
          zt1=ptrot(pttzz,b,b1)
          zt2=ptrot(pttzz,b,b2)
          zzt=zzt+a1(i1)*(zt1+zt2/z)
        endif

      enddo
      enddo
      colp=sqrt(colp)
      colt=sqrt(colt)
      if(colp.gt.1.)then
        kollth=nint(max(1.,colp))
        colp=fnuc*log(colp)
        zzp=sqrt(zzp)
        zzp=0.01*zzp*colp*bgl2
c saturation
        zzp=min(zzp,colp*abs(epscrx))
      else
        zzp=0.
      endif
      if(colt.gt.1.)then
        kollth=nint(max(1.,kollth+colt))
        colt=fnuc*log(colt)
        zzt=sqrt(zzt)
        zzt=0.01*zzt*colt*bgl2
c saturation
        zzt=min(zzt,colt*abs(epscrx))
      else
        zzt=0.
      endif
c      zzp=zzp*2.   !correction to have formula=MC
c      zzt=zzt*2.

c      print *,'ici',b,zzp,zzt,kollth,b2xnuc,colp,colt

      return
      end

cc------------------------------------------------------------------------
c      subroutine zzfz(zzp,zzt,kollth,b)
cc-----------------------------------------------------------------------
cc Z theoretical calculation of Z function
cc-----------------------------------------------------------------------
c      common /psar50/ zznuc,b2xnuc
c#include "aaa.h"
c#include "ems.h"
c#include "par.h"
c      common /ar3/ x1(7),a1(7)
c      external  pttcs,pprcs,pttzz,pprzz
c
c      sy=engy*engy
c      kollth=1
c
c      rs=r2had(iclpro)+r2had(icltar)+slopom*log(sy)
c      rpom=4.*.0389*rs
c
c      e1=exp(-1.)
c
c      colp=0.
c      colt=0.
c      zzt=0.
c      zzp=0.
c      if(maproj.eq.1.and.matarg.eq.1)return
c      do i1=1,7
c      do m=1,2
c        z=.5+x1(i1)*(m-1.5)
c        zv1=exp(-z)
c        zv2=(e1*z)
c        b1=sqrt(-rpom*log(zv1))
c        b2=sqrt(-rpom*log(zv2))
c
c
c        if(maproj.gt.1)then
c          cg1=ptrot(pprcs,b,b1)
c          cg2=ptrot(pprcs,b,b2)
c          colnuc=a1(i1)*(cg1+cg2/z)
c          colp=colp+colnuc
c          zt1=cg1*Zpair??(sy,b1)
c          zt2=cg2*Zpair??(sy,b2)
c          zzt=zzt+a1(i1)*(zt1+zt2/z)
c        endif
c        if(matarg.gt.1)then
c          cg1=ptrot(pttcs,b,b1)
c          cg2=ptrot(pttcs,b,b2)
c          colnuc=a1(i1)*(cg1+cg2/z)
c          colt=colt+colnuc
c          zp1=cg1*Zpair??(sy,b1)
c          zp2=cg2*Zpair??(sy,b2)
c          zzp=zzp+a1(i1)*(zp1+zp2/z)
c        endif
c
c      enddo
c      enddo
c      if(colt.eq.0.)then
c        zzt=sqrt(zzt)
c        kollth=nint(max(1.,sqrt(colp)))
c      elseif(colp.eq.0.)then
c        zzp=sqrt(zzp)
c        kollth=nint(max(1.,sqrt(colt)))
c      else
c        zzp=sqrt(zzp)
c        zzt=sqrt(zzt)
c        kollth=nint(max(1.,sqrt(colt*colp)))
c      endif
c
c
c      return
c      end
c

c------------------------------------------------------------------------
      function ptgau(func,bm,ipt,iqq)
c-----------------------------------------------------------------------
c impact parameter integration for impact parameters <bm -
c for nucleus-nucleus and hadron-nucleus cross-sections calculation
c ipt=1 : projectile, ipt=2 : target
c iqq=1 : elastic xsection, iqq=2 : inelastic cross section
c-----------------------------------------------------------------------
#include "aaa.h"
      common /ar3/ x1(7),a1(7)
      external func

      ptgau=0.
      do i=1,7
      do m=1,2
        b=bm*sqrt(.5+x1(i)*(m-1.5))
        ptgau=ptgau+func(b,ipt,iqq)*a1(i)
      enddo
      enddo
      ptgau=ptgau*bm**2*pi*.5
      return
      end

c------------------------------------------------------------------------
      function ptgau1(bm,ipt,iqq)
c-----------------------------------------------------------------------
c impact parameter integration for impact parameters >bm -
c for hadron-nucleus cross-sections calculation
c ipt=1 : projectile, ipt=2 : target
c iqq=1 : elastic xsection, iqq=2 : inelastic cross section
c-----------------------------------------------------------------------
#include "aaa.h"
      common /ar5/    x5(2),a5(2)

      ptgau1=0.
      if(ipt.eq.1)then
        difn=difnuc(maproj)
      else
        difn=difnuc(matarg)
      endif
      do i=1,2
        b=bm+x5(i)*difn
        ptgau1=ptgau1+ptfau(b,ipt,iqq)*a5(i)*exp(x5(i))*b*2.*pi*difn
      enddo
      return
      end
c------------------------------------------------------------------------
      function ptgau2(bm,iqq)
c-----------------------------------------------------------------------
c impact parameter integration for impact parameters >bm -
c for nucleus-nucleus cross-sections calculation
c iqq=1 : elastic xsection, iqq=2 : inelastic cross section
c-----------------------------------------------------------------------
#include "aaa.h"
      common /ar5/    x5(2),a5(2)

      ptgau2=0.
      difn=difnuc(maproj)+difnuc(matarg)
      do i=1,2
        b=bm+x5(i)*difn
        ptgau2=ptgau2+ptfauAA(b,iqq)*a5(i)*exp(x5(i))*b*2.*pi*difn
      enddo
      return
      end


c------------------------------------------------------------------------
      function ptfau(b,ipt,iqq)
c-----------------------------------------------------------------------
c ptfau - integrands for hadron-nucleus cross-sections calculation
c ipt=1 : projectile, ipt=2 : target
c iqq=1 : elastic xsection, iqq=2 : inelastic cross section
c-----------------------------------------------------------------------
#include "aaa.h"
      common /psar35/ anorm,anormp

      call psfz(iqq,gz2,b)

      if(ipt.eq.1)then
        ptfau=1.-max(0.,(1.-anormp*gz2))**maproj
      else
        ptfau=1.-max(0.,(1.-anorm*gz2))**matarg
      endif

      return
      end

c------------------------------------------------------------------------
      function ptfauAA(b,iqq)
c-----------------------------------------------------------------------
c ptfau - integrands for hadron-nucleus cross-sections calculation
c iqq=1 : elastic xsection, iqq=2 : inelastic cross section
c-----------------------------------------------------------------------
#include "aaa.h"
      common /ar3/    x1(7),a1(7)
      common /psar35/ anorm,anormp
      external pprcs

      ptfauAA=0.
      e1=exp(-1.)
      rs=r2had(iclpro)+r2had(icltar)+slopom*log(engy**2)
      rpom=4.*.0389*rs
      rpom=rpom*facmc
      do i1=1,7
      do m=1,2
        z=.5+x1(i1)*(m-1.5)
        zv1=exp(-z)
        zv2=(e1*z)
        b1=sqrt(-rpom*log(zv1))
        b2=sqrt(-rpom*log(zv2))
        call psfz(iqq,gz21,b1)
        call psfz(iqq,gz22,b2)
        ptfau1=max(0.,(1.-anorm*gz21))**matarg
        ptfau2=max(0.,(1.-anorm*gz22))**matarg
        cg1=ptrot(pprcs,b,b1)
        cg2=ptrot(pprcs,b,b2)
        ptfauAA=ptfauAA+a1(i1)*(cg1*(1.-ptfau1)+cg2*(1.-ptfau2)/z)
      enddo
      enddo
      ptfauAA=ptfauAA*rpom/2.
      ptfauAA=1.-max(0.,(1.-anormp*ptfauAA))**maproj

      return
      end

c------------------------------------------------------------------------
      function ptrot(func,s,b)
c-----------------------------------------------------------------------
c convolution of nuclear profile functions (axial angle integration)
c-----------------------------------------------------------------------
      common /ar8/ x2(4),a2
      external func

      ptrot=0.
      do i=1,4
        sb1=b**2+s**2-2.*b*s*(2.*x2(i)-1.)
        sb2=b**2+s**2-2.*b*s*(1.-2.*x2(i))
       ptrot=ptrot+(func(sb1)+func(sb2))
      enddo
      ptrot=ptrot*a2
      return
      end

c------------------------------------------------------------------------
      function pttcs(b0)
c-----------------------------------------------------------------------
c ptt - nuclear profile function value at imp param squared b*difnuc**2
c-----------------------------------------------------------------------
#include "aaa.h"
      common /psar34/ rrr,rrrm
      common /ar5/    x5(2),a5(2)
      common /ar9/    x9(3),a9(3)

      b=b0/difnuc(matarg)**2
      pttcs=0.
      zm=rrrm**2-b
      if(zm.gt.4.*b)then
        zm=sqrt(zm)
      else
        zm=2.*sqrt(b)
      endif

      do i=1,3
        z1=zm*(1.+x9(i))*0.5
        z2=zm*(1.-x9(i))*0.5
        quq=sqrt(b+z1**2)-rrr
        if (quq.lt.85.)pttcs=pttcs+a9(i)/(1.+exp(quq))
        quq=sqrt(b+z2**2)-rrr
        if (quq.lt.85.)pttcs=pttcs+a9(i)/(1.+exp(quq))
      enddo
      pttcs=pttcs*zm*0.5

      dt=0.
      do i=1,2
        z1=x5(i)+zm
        quq=sqrt(b+z1**2)-rrr-x5(i)
        if (quq.lt.85.)dt=dt+a5(i)/(exp(-x5(i))+exp(quq))
      enddo

      pttcs=pttcs+dt
      return
      end


c------------------------------------------------------------------------
      function pttzz(b0)
c-----------------------------------------------------------------------
c ptt - nuclear Z function value at imp param squared b*difnuc**2
c-----------------------------------------------------------------------
#include "aaa.h"
#include "par.h"
      common /psar34/ rrr,rrrm
      common /psar50/ zznuc,b2xnuc
      common /ar5/    x5(2),a5(2)
      common /ar9/    x9(3),a9(3)

      pttzz=0.
      b=b0/difnuc(matarg)**2
      absb=max(1.e-9,sqrt(b0))
      bsq=absb*absb
      zm=rrrm**2-b
      if(zm.gt.4.*b)then
        zm=sqrt(zm)
      else
        zm=2.*sqrt(b)
      endif

      do i=1,3
        z1=zm*(1.+x9(i))*0.5
        z2=zm*(1.-x9(i))*0.5
        quq=sqrt(b+z1**2)-rrr
        if (quq.lt.85.)pttzz=pttzz+a9(i)/(1.+exp(quq))
        quq=sqrt(b+z2**2)-rrr
        if (quq.lt.85.)pttzz=pttzz+a9(i)/(1.+exp(quq))
      enddo
      pttzz=pttzz*zm*0.5

      dt=0.
      do i=1,2
        z1=x5(i)+zm
        quq=sqrt(b+z1**2)-rrr-x5(i)
        if (quq.lt.85.)dt=dt+a5(i)/(exp(-x5(i))+exp(quq))
      enddo

      pttzz=max(0.,(pttzz+dt)-1.)*zznuc*exp(-bsq/2./b2xnuc)

      return
      end

c------------------------------------------------------------------------
      function pprcs(b0)
c-----------------------------------------------------------------------
c ppr - nuclear profile function value at imp param squared b*difnuc**2
c-----------------------------------------------------------------------
#include "aaa.h"
      common /psar41/ rrrp,rrrmp
      common /ar5/    x5(2),a5(2)
      common /ar9/    x9(3),a9(3)

      b=b0/difnuc(maproj)**2
      pprcs=0.
      zm=rrrmp**2-b
      if(zm.gt.4.*b)then
        zm=sqrt(zm)
      else
        zm=2.*sqrt(b)
      endif

      do i=1,3
        z1=zm*(1.+x9(i))*0.5
        z2=zm*(1.-x9(i))*0.5
        quq=sqrt(b+z1**2)-rrrp
        if (quq.lt.85.)pprcs=pprcs+a9(i)/(1.+exp(quq))
        quq=sqrt(b+z2**2)-rrrp
        if (quq.lt.85.)pprcs=pprcs+a9(i)/(1.+exp(quq))
      enddo
      pprcs=pprcs*zm*0.5

      dt=0.
      do i=1,2
        z1=x5(i)+zm
        quq=sqrt(b+z1**2)-rrrp-x5(i)
        if (quq.lt.85.)dt=dt+a5(i)/(exp(-x5(i))+exp(quq))
      enddo

      pprcs=pprcs+dt
      return
      end

c------------------------------------------------------------------------
      function pprzz(b0)
c-----------------------------------------------------------------------
c ppr - Z nuclear function value at imp param squared b*difnuc**2
c-----------------------------------------------------------------------
#include "aaa.h"
#include "par.h"
      common /psar41/ rrrp,rrrmp
      common /psar50/ zznuc,b2xnuc
      common /ar5/    x5(2),a5(2)
      common /ar9/    x9(3),a9(3)

      pprzz=0.
      b=b0/difnuc(maproj)**2
      absb=max(1.e-9,sqrt(b0))
      bsq=absb*absb
      zm=rrrmp**2-b
      if(zm.gt.4.*b)then
        zm=sqrt(zm)
      else
        zm=2.*sqrt(b)
      endif

      do i=1,3
        z1=zm*(1.+x9(i))*0.5
        z2=zm*(1.-x9(i))*0.5
        quq=sqrt(b+z1**2)-rrrp
        if (quq.lt.85.)pprzz=pprzz+a9(i)/(1.+exp(quq))
        quq=sqrt(b+z2**2)-rrrp
        if (quq.lt.85.)pprzz=pprzz+a9(i)/(1.+exp(quq))
      enddo
      pprzz=pprzz*zm*0.5

      dt=0.
      do i=1,2
        z1=x5(i)+zm
        quq=sqrt(b+z1**2)-rrrp-x5(i)
        if (quq.lt.85.)dt=dt+a5(i)/(exp(-x5(i))+exp(quq))
      enddo

      pprzz=max(0.,(pprzz+dt)-1.)*zznuc*exp(-bsq/2./b2xnuc)

      return
      end

c------------------------------------------------------------------------------
      function pscrse(ek,mapr,matg,iqq)
c------------------------------------------------------------------------------
c hadron-nucleus (hadron-proton) and nucl-nucl particle production cross section
c ek     - lab kinetic energy for the interaction
c maproj - projec mass number
c matarg - target mass number
c iqq=1    - ela cross section
c     >2   - ine cross section (2 used for cut (changing iomega), 3 uses table,
c                               4 used for ine without table)
c------------------------------------------------------------------------------
      dimension wk(3),wa(3),wb(3)
#include "aaa.h"
      common /psar33/ asect(7,4,7),asectn(7,7,7)
      common /psar34/ rrr,rrrm
      common /psar35/ anorm,anormp
      common /psar41/ rrrp,rrrmp
      external ptfau,ptfauAA

      pscrse=0.
      call idmass(1120,amt1)
      call idmass(1220,amt2)
      amtar=0.5*(amt1+amt2)
      if(matg.eq.1)amtar=amt1
      if(mapr.eq.1)then
        call idmass(idproj,ampro)
      else
        ampro=amtar
      endif
      egy=ek+ampro
c      p=sqrt(max(0.,egy**2-ampro**2))
      egy=sqrt( 2*egy*amtar+amtar**2+ampro**2 )
      call DefXminDf(dble(egy*egy))

      if(isetcs.le.1.or.iqq.ne.3)then
        maprojsave=maproj
        matargsave=matarg
        engysave=engy
        maproj=mapr
        matarg=matg
        engy=egy
        if(matg.eq.1.and.mapr.eq.1)then
          if(iqq.eq.1)then !sig ela
            call psfz(1,gz2,0.)
          else             !sig ine
            call psfz(2,gz2,0.)
          endif
          gin=gz2*pi*10.
        elseif(mapr.eq.1)then
          rad=radnuc(matg)
          bm=rad+2.
          rrr=rad/difnuc(matg)
          rrrm=rrr+log(9.)
          anorm=1.5/pi/rrr**3/(1.+(pi/rrr)**2)/difnuc(matg)**2
          if(iqq.ne.1)then
            gin=(ptgau(ptfau,bm,2,2)+ptgau1(bm,2,2))*10. !sig ine
          else
            gin=(ptgau(ptfau,bm,2,1)+ptgau1(bm,2,1))*10. !sig ela
          endif
        elseif(matg.eq.1)then
          rad=radnuc(mapr)
          bm=rad+2.
          rrrp=rad/difnuc(mapr)
          rrrmp=rrrp+log(9.)
          anormp=1.5/pi/rrrp**3/(1.+(pi/rrrp)**2)/difnuc(mapr)**2
          if(iqq.ne.1)then
            gin=(ptgau(ptfau,bm,1,2)+ptgau1(bm,1,2))*10. !sig ine
          else
            gin=(ptgau(ptfau,bm,1,1)+ptgau1(bm,1,1))*10. !sig ela
          endif
         else
          rad=radnuc(matg)+1.
          radp=radnuc(mapr)+1.
          bm=rad+radp+2.
          rrr=rad/difnuc(matg)
          rrrm=rrr+log(9.)
          rrrp=radp/difnuc(mapr)
          rrrmp=rrrp+log(9.)
          anorm=1.5/pi/rrr**3/(1.+(pi/rrr)**2)/difnuc(matg)**2
          anormp=1.5/pi/rrrp**3/(1.+(pi/rrrp)**2)/difnuc(mapr)**2
          if(iqq.ne.1)then
            gin=(ptgau(ptfauAA,bm,2,2)+ptgau2(bm,2))*10. !sig ine
          else
            gin=(ptgau(ptfauAA,bm,2,1)+ptgau2(bm,1))*10. !sig ela
          endif
        endif
        pscrse=gin
        maproj=maprojsave
        matarg=matargsave
        engy=engysave
      else
        ye=log10(max(1.,egy/1.5))+1.
        je=min(5,int(ye))

        wk(2)=ye-je
        wk(3)=wk(2)*(wk(2)-1.)*.5
        wk(1)=1.-wk(2)+wk(3)
        wk(2)=wk(2)-2.*wk(3)

        ya=matg
        ya=log(ya)/.69315+1.
        ja=min(int(ya),4)
        wa(2)=ya-ja
        wa(3)=wa(2)*(wa(2)-1.)*.5
        wa(1)=1.-wa(2)+wa(3)
        wa(2)=wa(2)-2.*wa(3)

        if(mapr.eq.1)then

          do i=1,3
            do m=1,3
              pscrse=pscrse+asect(je+i-1,iclpro,ja+m-1)*wk(i)*wa(m)
            enddo
          enddo

        else

          yb=mapr
          yb=log(yb)/.69315+1.
          jb=min(int(yb),4)
          wb(2)=yb-jb
          wb(3)=wb(2)*(wb(2)-1.)*.5
          wb(1)=1.-wb(2)+wb(3)
          wb(2)=wb(2)-2.*wb(3)

          do i=1,3
            do m=1,3
              do n=1,3
            pscrse=pscrse+asectn(je+i-1,jb+n-1,ja+m-1)*wk(i)*wa(m)*wb(n)
              enddo
            enddo
          enddo

        endif

        pscrse=exp(pscrse)
      endif
      return
      end

c------------------------------------------------------------------------------
      function eposcrse(ek,mapro,matar,id)
c------------------------------------------------------------------------------
c inelastic cross section of epos
c (id=0 corresponds to air)
c ek     - kinetic energy for the interaction
c maproj - projec mass number     (1<maproj<64)
c matarg - target mass number     (1<matarg<64)
c------------------------------------------------------------------------------
#include "aaa.h"

      eposcrse=0.
      if(id.eq.0)then
        do k=1,3
          mt=int(airanxs(k))
          eposcrse=eposcrse+airwnxs(k)*pscrse(ek,mapro,mt,3)
        enddo
      else
        eposcrse=pscrse(ek,mapro,matar,3)
      endif

      return
      end

c------------------------------------------------------------------------------
      function eposinecrse(ek,mapro,matar,id)
c------------------------------------------------------------------------------
c inelastic cross section of epos not using tabulated xs
c (id=0 corresponds to air)
c ek     - kinetic energy for the interaction
c maproj - projec mass number     (1<maproj<64)
c matarg - target mass number     (1<matarg<64)
c------------------------------------------------------------------------------
#include "aaa.h"

      eposinecrse=0.
      if(id.eq.0)then
        do k=1,3
          mt=int(airanxs(k))
          eposinecrse=eposinecrse+airwnxs(k)*pscrse(ek,mapro,mt,4)
        enddo
      else
        eposinecrse=pscrse(ek,mapro,matar,4)
      endif

      return
      end

c------------------------------------------------------------------------------
      function eposelacrse(ek,mapro,matar,id)
c------------------------------------------------------------------------------
c elastic cross section of epos
c (id=0 corresponds to air)
c ek     - kinetic energy for the interaction
c maproj - projec mass number     (1<maproj<64)
c matarg - target mass number     (1<matarg<64)
c------------------------------------------------------------------------------
#include "aaa.h"

      eposelacrse=0.
      if(id.eq.0)then
        do k=1,3
          mt=int(airanxs(k))
          eposelacrse=eposelacrse+airwnxs(k)*pscrse(ek,mapro,mt,1)
        enddo
      else
        eposelacrse=pscrse(ek,mapro,matar,1)
      endif

      return
      end


c------------------------------------------------------------------------------
      function eposcutcrse(ek,mapro,matar,id)
c------------------------------------------------------------------------------
c total cross section of epos
c (id=0 corresponds to air)
c ek     - kinetic energy for the interaction
c maproj - projec mass number     (1<maproj<64)
c matarg - target mass number     (1<matarg<64)
c------------------------------------------------------------------------------
#include "aaa.h"

      eposcutcrse=0.
      iomegasave=iomega
      iomega=2
      if(id.eq.0)then
        do k=1,3
          mt=int(airanxs(k))
          eposcutcrse=eposcutcrse+airwnxs(k)*pscrse(ek,mapro,mt,2)
        enddo
      else
        eposcutcrse=pscrse(ek,mapro,matar,2)
      endif
      iomega=iomegasave
      call DefXminDf(dble(engy*engy))

      return
      end

c------------------------------------------------------------------------------
      subroutine crseaaEpos(sigt,sigi,sigc,sige)
c------------------------------------------------------------------------------
c nucleus-nucleus (hadron) cross section of epos from simplified (realistic)
c simulations
c (id=0 corresponds to air)
c  sigt = sig tot
c  sigi = sig inelastic (cut + projectile diffraction)
c  sigc = sig cut
c  sige = sig elastic (includes target diffraction)
c------------------------------------------------------------------------------
#include "aaa.h"
      niter=1000
      if(idtarg.eq.0)then
        sigt=0.
        sigc=0.
        sigi=0.
        sige=0.
        sigd=0.
        sigql=0.
        do k=1,3
          matarg=int(airanxs(k))
          call epocrossc(niter,xsigt,xsigi,xsigc,xsige,xsigql,xsigd)
          sigt=sigt+airwnxs(k)*xsigt
          sigi=sigi+airwnxs(k)*xsigi
          sigc=sigc+airwnxs(k)*xsigc
          sige=sige+airwnxs(k)*xsige
          sigd=sigd+airwnxs(k)*xsigd
          sigql=sigql+airwnxs(k)*xsigql
        enddo
      else
        call epocrossc(niter,sigt,sigi,sigc,sige,sigql,sigd)
      endif
      if(ionudi.eq.2)then
        sige=sige+sigql      !add non-excited diffractive projectile to elastic
        sigi=sigi-sigql      !do not count non-excited diffractive projectile in inelastic
        if(maproj+matarg.gt.2)then
          sigc=sigc+sigd*0.95   !for absorbtion cross section remove 5% of the
                                !excited projectile diffractive cross section
                                !which "looks like" non excited (approximation)
        endif
      endif
      end



c-------------------------------------------------------------
      subroutine iniCR
c-------------------------------------------------------------
#include "aaa.h"
      if(idtarg.eq.0)then   !in case of Air target, initialize with Argon nucleus
        if(model.eq.6)then     !no Argon in Sibyll
          latarg=7
          matarg=14
        else
          latarg=20
          matarg=40
        endif
      endif
      if((idproj.ne.1120.and.(laproj.ne.-1.or.maproj.ne.1))
     &  .or.maproj.le.0)
     &call utstop('Invalid projectile setup !&')
c      if((idtarg.ne.1120.and.(latarg.ne.-1.or.matarg.ne.1))
c        .or.matarg.le.0)
c      call utstop('Invalid target setup !&')

      if(iabs(idtarg).ne.1120.and.iabs(idtarg).ne.1220.and.idtarg.ne.0)
     &  call utstop('Invalid target !&')
      if((((idtarg.eq.-1120.or.iabs(idtarg).eq.1220)
     &    .and.(latarg.ne.-1.or.matarg.ne.1))
     &    .and.(idtarg.ne.1120.or.latarg.lt.0))
     &    .or.matarg.le.0)
     &  call utstop('Invalid target setup !&') 
      end

c-------------------------------------------------------------
      subroutine iniXsection
c-------------------------------------------------------------
#include "aaa.h"
        if(isigma.eq.1.and.ionudi.ne.1)then
          write(ifmt,'(a)')
     &  '!---------------------------------------'
          write(ifmt,'(a)')
     &  '! Use isigma=2 for MC X section calc. '
          write(ifmt,'(a,i2,a,i2)')
     &  '!             isigma=',isigma,', ionudi=',ionudi
          write(ifmt,'(a)')
     &  '!---------------------------------------'
        endif
      end

c-----------------------------------------------------------------------
      function xsectionpar()
c-----------------------------------------------------------------------
#include "aaa.h"
      data ncntxs/0/
      save ncntxs
      ncntxs=ncntxs+1
      x=engy
      stot=21*x**0.173+44.*x**(-0.8) 
      sela=30.*(x-1)**(-3)+17*x**(-0.47)+0.3*alog(x)**2
      xsectionpar=max(30.,stot-sela)
      if(ncntxs.eq.1)then
      write(ifmt,*)'********* sig_inel =',xsectionpar
     .,'   parameterized!!! ' 
      endif
      end



c-----------------------------------------------------------------------
      subroutine sigmaint(g0,gz)
c-----------------------------------------------------------------------
c hadron-hadron cross sections integration
c-----------------------------------------------------------------------
      common /ar3/  x1(7),a1(7)
#include "aaa.h"
#include "par.h"
#include "sem.h"
#include "ems.h"
      double precision PhiExpo,vvv11,vvv12,vvv21!,HPhiDiff
     *,vvv22,ww01,ww02,ww11,ww12,ww21,ww22,gz(0:3)
     *,PhiExact!,vvv11e,vvv12e,vvv21e,vvv22e


      kollini=koll
      koll=1
      sy=engy**2
      rs=r2had(iclpro)+r2had(icltar)+slopom*log(sy)
        rpom=4.*.0389*rs  !((hbar*c)**2=0.0389 GeV^2.fm^2, and r2had in GeV^-2)
        rpom=rpom*facmc
        e1=exp(-1.)
c compute r2hads from proper calculation
c        r2hads(iclpro)=1.
c        r2hads(icltar)=1.
c        rdiff=0.!25*sqrt(rpom)
c        r2hads(iclpro)=sqrt(max(1.,sngl(HPhiDiff(sy,rdiff)
c     .                                 /om1intby(sy,rdiff,5,1))))
c        r2hads(icltar)=r2hads(iclpro)
c        if(ish.ge.2)write(ifch,*)
c        print *,
c     .              'r2hads:',r2hads(iclpro),r2hads(icltar)
        gz(0)=0.d0
        gz(1)=0.d0
        gz(2)=0.d0
        gz(3)=0.d0

        do i1=1,7
        do m=1,2

          z=.5+x1(i1)*(m-1.5)
          zv1=exp(-z)
          zv2=(e1*z)
          b1=sqrt(-rpom*log(zv1))
          b2=sqrt(-rpom*log(zv2))
          zz=0.!znurho

c          write(ifch,*)'xs b loop ->',i1,m,b1,b2

          if(isetcs.lt.0)then
            vvv11=max(0.d0,PhiExact(zz,zz,.5,1.d0,1.d0,sy,b1))
            vvv21=max(0.d0,PhiExact(zz,zz,1.,1.d0,1.d0,sy,b1))
          else
            vvv11=max(0.d0,PhiExpo(zz,zz,.5,1.d0,1.d0,sy,b1))
            vvv21=max(0.d0,PhiExpo(zz,zz,1.,1.d0,1.d0,sy,b1))
          endif

          if(isetcs.lt.0)then
            vvv12=max(0.d0,PhiExact(zz,zz,.5,1.d0,1.d0,sy,b2))
            vvv22=max(0.d0,PhiExact(zz,zz,1.,1.d0,1.d0,sy,b2))
          else
            vvv12=max(0.d0,PhiExpo(zz,zz,.5,1.d0,1.d0,sy,b2))
            vvv22=max(0.d0,PhiExpo(zz,zz,1.,1.d0,1.d0,sy,b2))
          endif

          ww11=1.d0-vvv11
          ww12=1.d0-vvv12
          ww21=1.d0-vvv21
          ww22=1.d0-vvv22
          ww01=vvv21-2d0*vvv11+1d0
          ww02=vvv22-2d0*vvv12+1d0

          gz(0)=gz(0)+a1(i1)*(ww01+ww02/z)
          gz(1)=gz(1)+a1(i1)*(ww11+ww12/z)
          gz(2)=gz(2)+a1(i1)*(ww21+ww22/z)
          gz(3)=gz(3)+a1(i1)*(ww11*z+ww12/z*(1.-log(z)))

c          write(ifch,*)'<- xs b loop',gz

        enddo
        enddo
        g0=pi*rpom*10./2.                 !common factor (pi*rpom because of b->z, 10 to have mbarn and /2. because z is from -1 to 1 but we need 0 to 1.

        koll=kollini

      return
      end

c-----------------------------------------------------------------------
      subroutine fsigmak(k)
c-----------------------------------------------------------------------
c     hadron-hadron cut cross section for pair k
c     Warning : it gives slightly different value even for p-p because of
c     the updated value of xmxrem (system dependent)
c-----------------------------------------------------------------------
      common /ar3/  x1(7),a1(7)
#include "aaa.h"
#include "par.h"
#include "sem.h"
#include "ems.h"
      common/cnparticip/jproj(2,mamx),jtarg(2,mamx),efluct(6,mamx)
      double precision PhiUnit,vvv21,vvv22,ww21,ww22,gz2,gzd1,gzd2,gzd0
     *                ,om1intby,gzh0,gzh1,gzh2,om1intbci!,om1intbhj


      sy=engy**2
      rs=r2had(iclpro)+r2had(icltar)+slopom*log(sy)
      rpom=4.*.0389*rs          !((hbar*c)**2=0.0389 GeV^2.fm^2, and r2had in GeV^-2)
      rpom=rpom*facmc
      e1=exp(-1.)
      gz2=0.d0
      gzd0=0.d0
      gzh0=0.d0
      zzp=zparpnx(1,k)     !centrality dep of binary scaling is not good
      zzt=zpartnx(1,k)     !if nuclear xs used
      ip=iproj(k)
      it=itarg(k)
      ef1=efluct(1,ip)
      ef2=efluct(2,it)
      ef3=efluct(3,ip)
      ef4=efluct(4,it)
      eflup=ef1
      if(eflup.gt.0.)eflup=eflup+ef3
      eflut=ef2
      if(eflut.gt.0.)eflut=eflut+ef4
      imin=ntymin
      imax=ntymax
c      zrhoincsave=zrhoinc
c      epscrssave=epscrs
c      epscrgsave=epscrg
c      zrhoinc=0.
c      epscrs=0.
c      epscrg=0.

        do i1=1,7
        do m=1,2

          z=.5+x1(i1)*(m-1.5)
          zv1=exp(-z)
          zv2=(e1*z)
          b1=sqrt(-rpom*log(zv1))
          call Gfunpar(zzp,zzt,eflup,eflut,1,-1,b1,sy,a,b,c,d,e,f,g)
          do i=imin,imax          !initialization of parameters for b1
            call Gfunpar(zzp,zzt,eflup,eflut,1,i,b1,sy,a,b,c,d,e,f,g)
            call Gfunpar(zzp,zzt,eflup,eflut,2,i,b1,sy,a,b,c,d,e,f,g)
          enddo
          vvv21=max(0.d0,PhiUnit(1.d0,1.d0))
          ww21=1.d0-vvv21
c sum of all diffractive diagram
          gzd1=max(0d0,om1intby(sy,b1,-50,1))
c for Npom
c          gzh1=max(0d0,om1intbhj(b1,k)) !max(0d0,om1intbci(0))
          gzh1=max(0d0,om1intbci(1))


          b2=sqrt(-rpom*log(zv2))
          call Gfunpar(zzp,zzt,eflup,eflut,1,-1,b2,sy,a,b,c,d,e,f,g)
          do i=imin,imax          !initialization of parameters for b1
            call Gfunpar(zzp,zzt,eflup,eflut,1,i,b2,sy,a,b,c,d,e,f,g)
            call Gfunpar(zzp,zzt,eflup,eflut,2,i,b2,sy,a,b,c,d,e,f,g)
          enddo
          vvv22=max(0.d0,PhiUnit(1.d0,1.d0))
          ww22=1.d0-vvv22
c sum of all diffractive diagram
          gzd2=max(0d0,om1intby(sy,b2,-50,1))
c for Npom
c          gzh2=max(0d0,om1intbhj(b2,k)) !max(0d0,om1intbci(0))
          gzh2=max(0d0,om1intbci(1))

c integral
          gz2=gz2+a1(i1)*(ww21+ww22/z)
          gzd0=gzd0+a1(i1)*(gzd1+gzd2/z)
          gzh0=gzh0+a1(i1)*(gzh1+gzh2/z)

        enddo
      enddo

!inelastic
        sigmak(1,k)=gz2*dble(pi*rpom*10./2.) !common factor (pi*rpom because of b->z, 10 to have mbarn and /2. because z is from -1 to 1 but we need 0 to 1.
!cut
        sigmak(2,k)=(gz2-gzd0)*dble(pi*rpom*10./2.) !common factor (pi*rpom because of b->z, 10 to have mbarn and /2. because z is from -1 to 1 but we need 0 to 1.
c        om1MeanCol(1,k)=gzh0/max(0.0001d0,gz2)!*sigmak(2,k)/sigmak(1,k)
c     .                     *xfrachard
        om1MeanCol(1,k)=gzh0/max(0.0001d0,gz2-gzd0)
c        om1MeanCol(1,k)=xfrachard*gzh0/max(0.0001d0,gz2-gzd0)
c      om1MeanCol(1,k)=max(om1MeanCol(1,k)**(1d0/dble(koll)),
c     .                      om1MeanCol(2,k))
c        om1MeanCol(1,k)=gzh0/max(0.0001d0,gz2)      !mean hard npom from AGK (correspond to npr(4,k) if esature properly adjusted) for pair k in nucleus
c        om1MeanCol(2,k)=max(0d0,om1intbhj(bk(k),-k))      !hard npom from AGK (correspond to npr(4,k) if esature properly adjusted) for given bk(k)

c        print *,'fsigmak',k,sigmak(1,k),sigmak(2,k),om1MeanCol(1,k)
c     .                   ,zzp,zzt,eflup,eflut,imin,imax
c     .                                             ,om1MeanCol(2,k)
c        stop
c      zrhoinc=zrhoincsave
c      epscrs=epscrssave
c      epscrg=epscrgsave
      return
      end

c-----------------------------------------------------------------------
      subroutine sigmadiff(gzd,iqq,imod)
c-----------------------------------------------------------------------
c hadron-hadron diff cross sections integration
c iqq=1 : real diagram
c iqq=-1 : renormalized diagrams for MC
c imod - calculation without (0-diagram) or with Phi (1-cross section) 
c-----------------------------------------------------------------------
      common /ar3/  x1(7),a1(7)
#include "aaa.h"
#include "par.h"
#include "sem.h"
#include "ems.h"
      double precision om1intby,gzd(0:4),gzd1,gzd2

      gzd(0)=0.d0
      gzd(1)=0.d0
      gzd(2)=0.d0
      gzd(3)=0.d0
      gzd(4)=0.d0

      if(.not.((abs(iqq).eq.1.or.iqq.eq.-10).and.iomega.lt.2))return

      imin=idxD0
      imax=idxD1
      sy=engy**2
      rs=r2had(iclpro)+r2had(icltar)+slopom*log(sy)
        rpom=4.*.0389*rs  !((hbar*c)**2=0.0389 GeV^2.fm^2, and r2had in GeV^-2)
        rpom=rpom*facmc
        e1=exp(-1.)

        do i1=1,7
        do m=1,2

          z=.5+x1(i1)*(m-1.5)
          zv1=exp(-z)
          zv2=(e1*z)
          b1=sqrt(-rpom*log(zv1))
          do i=imin,imax          !initialization of parameters for b1
            call Gfunpar(zzp,zzt,0.,0.,1,i,b1,sy,a,b,c,d,e,f,g)
            call Gfunpar(zzp,zzt,0.,0.,2,i,b1,sy,a,b,c,d,e,f,g)
          enddo

c sum of all diffractive diagram
          gzd1=a1(i1)*max(0d0,om1intby(sy,b1,-50,imod))
          gzd(0)=gzd(0)+gzd1
c soft SD+ diffractive diagram
          gzd1=a1(i1)*max(0d0,om1intby(sy,b1,iqq*6,imod))
          gzd(1)=gzd(1)+gzd1
c soft SD- diffractive diagram
          gzd1=a1(i1)*max(0d0,om1intby(sy,b1,iqq*7,imod))
          gzd(2)=gzd(2)+gzd1
c both side diffractive diagram
          gzd1=a1(i1)*max(0d0,om1intby(sy,b1,iqq*8,imod))
          gzd(3)=gzd(3)+gzd1
c central exclusive diffractive diagram
          gzd1=a1(i1)*max(0d0,om1intby(sy,b1,iqq*9,imod))
          gzd(4)=gzd(4)+gzd1

          b2=sqrt(-rpom*log(zv2))
          do i=imin,imax          !initialization of parameters for b1
            call Gfunpar(zzp,zzt,0.,0.,1,i,b2,sy,a,b,c,d,e,f,g)
            call Gfunpar(zzp,zzt,0.,0.,2,i,b2,sy,a,b,c,d,e,f,g)
          enddo
c sum of all diffractive diagram
          gzd2=a1(i1)*max(0d0,om1intby(sy,b2,-50,imod))/z
          gzd(0)=gzd(0)+gzd2
c soft SD+ diffractive diagram
          gzd2=a1(i1)*max(0d0,om1intby(sy,b2,iqq*6,imod))/z
          gzd(1)=gzd(1)+gzd2
c soft SD- diffractive diagram
          gzd2=a1(i1)*max(0d0,om1intby(sy,b2,iqq*7,imod))/z
          gzd(2)=gzd(2)+gzd2
c both side diffractive diagram
          gzd2=a1(i1)*max(0d0,om1intby(sy,b2,iqq*8,imod))/z
          gzd(3)=gzd(3)+gzd2
c central diffractive diagram
          gzd2=a1(i1)*max(0d0,om1intby(sy,b2,iqq*9,imod))/z
          gzd(4)=gzd(4)+gzd2


        enddo
        enddo

      return
      end

c-----------------------------------------------------------------------
      subroutine sigmahard(gzh)
c-----------------------------------------------------------------------
c hadron-hadron hard cross sections integration
c-----------------------------------------------------------------------
      common /ar3/  x1(7),a1(7)
#include "aaa.h"
#include "par.h"
#include "sem.h"
#include "ems.h"
      double precision om1intbci,om1intbhj,gzh(0:2),gzh1,gzh2,xnh1 !,exnh1

      gzh(0)=0.d0
      gzh(1)=0.d0
      gzh(2)=0.d0
      zz=0.
      imin=ntymin
      imax=ntymax
      sy=engy**2
      iomegasave=iomega
      iomega=2
      call DefXminDf(dble(sy))
      imin=ntynd
      imax=ntynd

      rs=r2had(iclpro)+r2had(icltar)+slopom*log(sy)
        rpom=4.*.0389*rs  !((hbar*c)**2=0.0389 GeV^2.fm^2, and r2had in GeV^-2)
        rpom=rpom*facmc
        e1=exp(-1.)

        do i1=1,7
        do m=1,2

          z=.5+x1(i1)*(m-1.5)
          zv1=exp(-z)
          zv2=(e1*z)
          b1=sqrt(-rpom*log(zv1))
          b2=sqrt(-rpom*log(zv2))
          do i=imin,imax          !initialization of parameters for b1
            call Gfunpar(0.,0.,0.,0.,1,i,b1,sy,a,b,c,d,e,f,g)
          enddo
          vvv21=1d0 !max(0.d0,1d0-PhiExpo(zz,zz,1.,1.d0,1.d0,sy,b1))
c     .                      -om1intby(sy,b1,-50))
          vvv22=1d0 !max(0.d0,1d0-PhiExpo(zz,zz,1.,1.d0,1.d0,sy,b2))
c     .                      -om1intby(sy,b2,-50))

c all diagrams
          gzh1=a1(i1)*max(0d0,om1intbci(0))
          gzh(0)=gzh(0)+gzh1
c hard diagram
          gzh1=a1(i1)*max(0d0,om1intbhj(b1,1))
          gzh(1)=gzh(1)+gzh1
c hard diagram with bsatur
          gzh1=1d0
          gzh2=0d0
          n=0
          xnh1=max(0d0,om1intbci(1))
c          if(bsatur.gt.0.)then
c            exnh1=exp(-xnh1)
c            do while(n.lt.npommx.and.(dble(n).le.xnh1.or.gzh1.gt.1d-7))
c              n=n+1
c              gzh1=xnh1**n*facto(n)*exnh1
c              gzh2=gzh2+gzh1*n**bsatur
c            enddo
c          else
            gzh2=xnh1
c          endif
          gzh1=a1(i1)*gzh2
          gzh(2)=gzh(2)+gzh1

          do i=imin,imax          !initialization of parameters for b1
            call Gfunpar(0.,0.,0.,0.,1,i,b2,sy,a,b,c,d,e,f,g)
          enddo
c all diagrams
          gzh2=a1(i1)*max(0d0,om1intbci(0))/z
          gzh(0)=gzh(0)+gzh2
c hard diagram
          gzh2=a1(i1)*max(0d0,om1intbhj(b2,1))/z
          gzh(1)=gzh(1)+gzh2
c hard diagram with bsatur
          gzh1=1d0
          gzh2=0d0
          n=0
          xnh1=max(0d0,om1intbci(1))
c          if(bsatur.gt.0.)then
c            exnh1=exp(-xnh1)
c            do while(n.lt.npommx.and.(dble(n).le.xnh1.or.gzh1.gt.1d-7))
c              n=n+1
c              gzh1=xnh1**n*facto(n)*exnh1
c              gzh2=gzh2+gzh1*n**bsatur
c            enddo
c          else
            gzh2=xnh1
c          endif
          gzh2=a1(i1)*gzh2/z
          gzh(2)=gzh(2)+gzh2

        enddo
        enddo

      iomega=iomegasave
      call DefXminDf(dble(sy))

      return
      end

c-----------------------------------------------------------------------
      subroutine xsigma
c-----------------------------------------------------------------------
c hadron-hadron and hadron-nucleus cross sections calculation
c b - impact parameter squared (in case of hadron-nucleus interaction);
c-----------------------------------------------------------------------
#include "aaa.h"
#include "sem.h"
      double precision gz(0:3),GZ0(2),gzp(0:4),gzh(0:2)
c Model 2 Common
      COMMON /Q_AREA1/  IA(2),ICZ,ICP
      COMMON /Q_AREA6/  PIQGS,BM,AM
      COMMON /Q_AREA15/ FP(5),RQ(5),CD(5)
      COMMON /Q_AREA7/  RP1
      COMMON /Q_AREA16/ CC(5)
      double precision RP1,FP,RQ,CD,PIQGS,BM,AM,CC,GDP,GDT,GDD

C...Total cross sections in Pythia
      double precision SIGT
      COMMON/PYINT7/SIGT(0:6,0:6,0:5)

c Model 5 Common
      COMMON/HIPARNT/HIPR1(100), IHPR2(50), HINT1(100), IHNT2(50)

c theoretical cross sections
      sigcut=0.
      sigtot=0.
      sigela=0.
      sigine=0.
      sigtotold=0.
      sigtotf=0.
      sigelaold=0.
      sigelaf=0.
      sigineold=0.
      sigineaa=0.
      sigtotaa=0.
      sigelaaa=0.
      sigcutaa=0.
      sigdif=0.
      sloela=0.
      sigsd=0.
      sigdd=0.
      sigcd=0.
      xnpomcol=0.
      xnhardcol=0.
      xnsatcol=0.
c simulated cross sections
      sigineex=0.          !calculated in ems if isigma>0
      sigdifex=0.
      sigsdex=0.


      call utpri('xsigma',ish,ishini,4)

      xs=xsectionpar()      !ctp just as a reference not used here anymore)
c      rexndi(iclpro)=rexndii(iclpro)
c      rexndi(icltar)=rexndii(icltar)
        if(ish.ge.2)write(ifch,*)'Xsigma : rexdif/ndi=',rexdif(iclpro)
     &                                                 ,rexdif(icltar)
     &                                                 ,rexndi(iclpro)
     &                                                 ,rexndi(icltar)

      if(model.eq.1)then                        !epos

        if(icltar.ne.2)stop'xsigma: only icltar=2 allowed.'

        call sigmaint(g0,gz)
        call sigmadiff(gzp,-10,1)
        call sigmahard(gzh)
        
c       theoretical number of Pomerons
        if(iokoll.eq.0)then
          xnpomcol=gzh(0)/max(0.0001d0,gz(2)-gzp(0))
          xnhardcol=gzh(1)/max(0.0001d0,gz(2)-gzp(0))
c        xnsatcol=gzh(2)/max(0.0001d0,gz(2)-gzp(0))
          xfrachard=min(1.,xnhardcol/xnpomcol)
          xnsatcol=xfrachard*gzh(0)/max(0.0001d0,gz(2)-gzp(0))
        else
          xnpomcol=1d0
          xnhardcol=1d0
          xfrachard=1d0
          xnsatcol=1d0
        endif
          
        print *,'col',xnpomcol,xnhardcol,xnsatcol,xfrachard 

        sigdifold=gzp(0) * g0   !xs in mb

        sigelaold=g0*gz(0)               !elastic cross section
        rs=g0*0.4091    !=g0/pi/10.*2./4./.0389
        if(gz(1).gt.0d0)sloela=2.*rs*gz(3)/gz(1)

        sigineold=g0*gz(2)               !inelastic pomerons cross-section
        sigtotold=2.*g0*gz(1)                  !tot cross-section
        sigcut=sigineold-sigdifold             !cut cross section
        x=engy
c fit to data
        !sigtotf=14.5*x**0.21+20.*x**(-0.2)+19.*(x-1.)**(-1)
        !sigelaf=35.*(x-1)**(-2.8)+17.*x**(-0.47)+0.31*log(x)**2
        sigtotf=21*x**0.17+44.*x**(-0.8) 
        sigelaf=30.*(x-1)**(-3)+17*x**(-0.47)+0.3*alog(x)**2
c        sigtotfp=sigtotf
c        sigelafp=sigelaf
        if(iclpro.eq.1)then        !pi+p
          sigtotf=10.*(x-1)**(-3)+16.*x**0.13+40.*x**(-1.2)
          sigelaf=25.*(x-1)**(-3)+6.*x**(-0.4)+0.15*log(x)**2.
        elseif(iclpro.eq.3)then    !K+p
          sigtotf=13.*x**0.15+35.*x**(-1.5)
          sigelaf=15.*(x-1)**(-3)+5.*x**(-0.4)+0.1*log(x)**2
        elseif(iclpro.eq.4)then    !D+p
          sigtotf=0.!12.5*x**0.15+35.*x**(-1.5)
          sigelaf=0.!15.*(x-1)**(-3)+3.*x**(-0.4)+0.2*alog(x)**2
        endif

c        if(engy.lt.20.)then
c          sigcoul=max(0.,sigtotf-sigtotold)
c        else
c          sigcoul=0.
c        endif

        sigdif=sigdifold
c        sigdelaf=max(0.,(sigelaf+sigineold-sigtotf))
c        sigdelaf=max(0.,(sigelaf-sigelaold))

        sigsd=(gzp(1)+gzp(2)) * g0 !xs in mb
        sigcd=gzp(3) * g0 !xs in mb
        sigdd=gzp(4) * g0       !xs in mb
c        print *,(gzp(1)+gzp(2)) * g0,gzp(3)*g0,gzp(4)*g0,sigdif
c        if(r2hads(iclpro)*r2hads(icltar).gt.0.)then
c          dtot=gzp(1)+gzp(2)+gzp(3)+gzp(4)
c          sigsd=min(sigsd,sigdif*sngl(gzp(1)+gzp(2))/dtot) !if SD contribution is larger than diff diagram contribution
c          sigdd=max(sigdd,sigdif*sngl(gzp(4))/dtot)
c        endif
        sigdela=0.

c        if(engy.lt.10.)sigela=max(sigelaf,sigela)
        sigine=sigineold
        sigela=sigelaold
        if(ionudi.ne.1)then
          sigela=sigela+sigdela
          sigine=sigine-sigdela
        endif
c        sigine=xs    !??????????
        if(engy.lt.30.)then
          sigcoul=max(0.,sigelaf-sigela)
        else
          sigcoul=0.
        endif
        sigela=sigela+sigcoul
        sigtot=sigine+sigela
        sigineaa=eposcrse(ekin,maproj,matarg,idtarg)     !TP????? to update !


      elseif(model.eq.2)then

        g0=real(PIQGS*RP1/CD(ICZ)*AM**2*10.D0)
        CALL m2XXFZ(0.D0,GZ0)
        gz(1)=GZ0(1)
        gz(2)=GZ0(2)
        gz(3)=0d0
        rs=g0*0.4091    !=g0/pi/10.*2./4./.0389
        sigcut=g0*gz(2)/2.               !cut pomerons cross-section
        sigtot=g0*gz(1)                  !tot cross-section
        gz(0)=sigtot-sigcut
        sigela=gz(0)*CC(ICZ)*CC(2)       !elastic cross section
c GDP - projectile diffraction cross section
        GDP=(1.D0-CC(ICZ))*CC(2)*gz(0)
c GDT - target diffraction cross section
        GDT=(1.D0-CC(2))*CC(ICZ)*gz(0)
c  GDD - double diffractive cross section
        GDD=(1.D0-CC(ICZ))*(1.D0-CC(2))*gz(0)
        sigsd=GDT+GDP
        sigdd=GDD
        sigdif=sigsd+sigdd
        sigine=sigcut+sigdif
        rs=g0*0.4091    !=g0/pi/10.*2./4./.0389
        if(gz(1).gt.0.)sloela=2.*rs*gz(3)/gz(1)
        sigdifold=sigtot-sigcut-sigela       !diffractive cross section
        sigineaa=qgsincs

      elseif(model.eq.3)then

        call m3SIGMA(ekin,idproj,1120,1,1,sigi,sige)
        sigine=sigi
        sigela=sige
        sigcut=sigine
        sigtot=sigine+sigela
        sigdif=sigtot-sigcut-sigela       !diffractive cross section
        sigineaa=gheincs

      elseif(model.eq.4)then         !PYTHIA

        sigsd=sngl(SIGT(0,0,2)+SIGT(0,0,3))
        sigela=sngl(SIGT(0,0,1))
        sigcut=sngl(SIGT(0,0,5))
        sigtot=sngl(SIGT(0,0,0))
        sigine=sigtot-sigela
        sigdif=sigtot-sigcut-sigela       !diffractive cross section
        sigineaa=pytincs

      elseif(model.eq.5)then           !HIJING

        sigsd=HIPR1(33)*HINT1(12)
        sigdif=0.
        sigcut=0.
        sigtot=HINT1(13)
        sigine=HINT1(12)
        sigela=sigtot-sigine
        sigineaa=hijincs

      elseif(model.eq.6)then                  !for Sibyll

        call m6SIGMA(iclpro,engy,stot,sela,sine,sdifr,slela,Rho)
        sigtot=stot
        sigela=sela
        sigine=sine
        sigdif=sdifr
        sloela=slela
        sigcut=sigtot-sigdif-sigela     ! cut cross section
        sigsd=sigdif/2.
        sigineaa=sibincs

      elseif(model.eq.7.or.model.eq.11)then                  !for QGSJET-II

        call m7SIGMA(stot,scut,sine,slela)
        sigtot=stot
        sigcut=scut
        sigine=sine
        sloela=slela
        sigela=sigtot-sigine     ! elastic cross section
        sigdif=sigine-sigcut
        sigsd=sigdif
        sigineaa=qgsIIincs

      elseif(model.eq.8)then                  !for PHOJET

        call m8SIGMA(stot,scut,sine,sela,slela,ssd)
        sigtot=stot
        sigcut=scut
        sigine=sine
        sloela=slela
        sigela=sela 
        sigdif=sigine-sigcut
        sigsd=ssd
        sigineaa=phoincs

      elseif(model.eq.9)then                  !for Fluka

        call m9SIGMA(stot,sine,sela)
        sigtot=stot
        sigine=sine
        sigcut=sigine
        sigela=sela 
        sigineaa=fluincs

      elseif(model.eq.10)then                  !for Urqmd

        sigtot=urqincs
        sigineaa=urqincs

      endif


             if(isigma.ge.1)then  !===============!

      if(ish.ge.1.and.noebin.ge.0)
     *write (ifmt,225)engy,ekin,sigtot,sigtotf,sigtotold
     *,sigine,sigtotf-sigelaf,sigineold
     *,sigela,sigelaf,sigelaold,sigcut,sloela,sigdif,sigsd
     *,sigineaa
      if(ish.ge.1.and.ifch.ne.ifmt)
     *write (ifch,225)engy,ekin,sigtot,sigtotf,sigtotold
     *,sigine,sigtotf-sigelaf,sigineold
     *,sigela,sigelaf,sigelaold,sigcut,sloela,sigdif,sigsd
     *,sigineaa

c (from tabulation) for pA/AA
      if((isigma.eq.2.and.noebin.ge.0)
     &    .or..not.(maproj.eq.1.and.matarg.eq.1.or.iokoll.ne.0))then

        sigtotaa=0.
        sigelaaa=0.
        if(model.eq.1)then
          if(isigma.ne.2)then
            if(maproj.gt.5.and.matarg.gt.5)then
              if(isetcs.eq.3)then
                write(ifmt,*)
     &        'Total and elastic cross-sections may be wrong,',
     &        'use isigma=2 instead ! ',
     &        '(if you care about it ...)'
              else
                write(ifmt,*)
     &        'Cross-sections may be wrong, use isigma=2 instead ! ',
     &        '(if you care about it ...)'
                sigineaa=eposinecrse(ekin,maproj,matarg,idtarg)
              endif
            endif
c eposcrse depends of ionudi while eposinecrse corresponds to ionudi=1 always
            sigelaaa=eposelacrse(ekin,maproj,matarg,idtarg)
            sigtotaa=sigelaaa+sigineaa
            sigcutaa=eposcutcrse(ekin,maproj,matarg,idtarg)
            if(ionudi.gt.1)then
c First order approximation. Better to use isigma=2 for that
              difpart=max(0.,sigineaa-sigcutaa)
c  excited diffractive part
              sigqela=0d0
              sdpart=1.d0-sigqela
              sigqela=sigqela*difpart
              sigineaa=sigineaa-sigqela
              sigelaaa=sigelaaa+sigqela
c here cut is absorbtion xs : cut + 95 % of excited diff.
              sigcutaa=sigcutaa+0.95*difpart*sdpart
            elseif(ionudi.eq.0)then
              write(ifmt,*)
     &        'Cross-section can not be calculated with ionudi=0'
            endif
          else
              write(ifmt,*)
     &        'Compute EPOS Cross-section (can take a while...)'
            call crseaaEpos(sigtotaa,sigineaa,sigcutaa,sigelaaa)
          endif
        elseif(isigma.eq.2.and.matarg.gt.1)then
         write(ifmt,*)
     &        'Compute model Cross-section (can take a while...)'
         call crseaaModel(sigtotaa,sigineaa,sigcutaa,sigelaaa)
        endif
        if(ish.ge.1.and.noebin.ge.0)
     &  write (ifmt,226)sigtotaa,sigineaa,sigcutaa,sigelaaa
        if(ish.ge.1.and.ifch.ne.ifmt)
     &  write (ifch,226)sigtotaa,sigineaa,sigcutaa,sigelaaa

      endif


             endif  !================!





225   format(' hadron-proton cross sections for ',f10.2,' GeV',
     *'  (ekin:',g13.5,' GeV)'/
     *4x,'total cross section:           ',f8.2,3x,f8.2,3x,f8.2/
     *4x,'inelastic cross section:       ',f8.2,3x,f8.2,3x,f8.2/
     *4x,'elastic cross section:         ',f8.2,3x,f8.2,3x,f8.2/
     *4x,'cut cross section:             ',f8.2/
     *4x,'elastic slope parameter:       ',f8.2/
     *4x,'diffr. cross section:          ',f8.2,14x,f8.2/
     *4x,'inelastic (tab) cross section: ',f8.2)
 226  format(' hadron/nucleus-hadron/nucleus cross sections'/
     *4x,'total pA/AA cross section:     ',f8.2/
     *4x,'inelastic pA/AA cross section: ',f8.2/
     *4x,'cut pA/AA cross section:       ',f8.2/
     *4x,'elastic pA/AA cross section:   ',f8.2)

      call utprix('xsigma',ish,ishini,4)

      return
      end




