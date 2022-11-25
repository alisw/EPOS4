C
C  This file is part of EPOS4
C  Copyright (C) 2022 research institutions and authors (See CREDITS file)
C  This file is distributed under the terms of the GNU General Public License version 3 or later
C  (See COPYING file for the text of the licence)
C

csp*********************************************************
csp         Test Routines
csp*********************************************************
c-----------------------------------------------------------------------
      subroutine psaevp
c-----------------------------------------------------------------------
c Warning : in external PDF, PDF for q include val q + sea q
c factor for each quark type is q**2 ((2/3)**2=1/2.25 and (-1/3)**2=1/9
c-----------------------------------------------------------------------
#include "aaa.h"
#include "sem.h"
#include "par.h"
      double precision pifpartone !,fzeroVal
c      dimension q2imin(2)
      qq=xpar1
      jmod=nint(xpar2)
      iologb=1
      engysave=engy
      iclprosave=iclpro
      icltarsave=icltar
      iclpro=2
      icltar=2

      q2cmin(1)=max(q2nmin,xpar11)
      q2cmin(2)=q2cmin(1)
      call ipoBasics(q2cmin)
      call idmass(1120,amtarg)

      do i=1,nrbins
        if(iologb.eq.0)then
          xx=xminim+(xmaxim-xminim)*(i-.5)/nrbins
        else
          xx=xminim*(xmaxim/xminim)**((i-.5)/nrbins)
        endif

        s=qq*(1.-xx)/xx+amtarg**2 
        engy=sqrt(s)
c        alpff(iclpro)=gamhad(iclpro) *SoftAtt(dble(xx),q2sft,1)
c        alpff(icltar)=gamhad(icltar) *SoftAtt(dble(xx),q2sft,2)
c        betff(1)=   -alppar   !done in ipoBasics
c        betff(2)=   -alppar
c        call fzeroQ2min(s,dble(xx),1,q2cmin(1))   !set q2pmin(1)
        q2cmin(2)=q2cmin(1)
        q2pmin(1)=q2cmin(1)
        q2pmin(2)=q2pmin(1)
c        if(jmod.le.1.or.jmod/10.eq.1.or.jmod.eq.5.or.jmod.eq.56
c     & .or.jmod.eq.4.or.jmod/10.eq.4)then
c        call ipoBasics(q2cmin)
c        print *,qq,xx,q2cmin(1),dels(1),glu2sea(1),betpom
         if(isetcs.gt.-2.and.(
     &     jmod.le.1.or.jmod.eq.16.or.jmod.eq.4.or.jmod.eq.46
     &              .or.jmod.eq.5.or.jmod.eq.56))then
c initialization of some variable needed for our stuff
           inicnt=inicnt+1
           call psaini
           call mkParamTables
           call psdsig(qq,s)      !calculate sigma_qq_p use for color dipole approach
         endif
c        endif

        ar(i,1)=xx
        ar(i,2)=0.
        ar(i,3)=0.
        if(jmod.eq.0)then            !evolution+matrix element (ours) without screening
          ww=qq/xx
          ar(i,3)=max(psdsh1(ww,qq,1,dqsh,6,0)+psdsh1(ww,qq,1,dqsh,6,1),
     *    psdh(ww,qq,1,0,0)+psdh(ww,qq,1,0,1)+
     *    psdsh1(ww,qq,1,dqsh,0,0)+psdsh1(ww,qq,1,dqsh,0,1)
     *    )/(4.*pi**2*alfe)*qq
c          ar(i,3)=fzeroVal(dble(xx),qq,2,2)
        elseif(jmod.eq.1)then        !evolution+matrix element (ours) with screening
          ww=qq/xx
          ar(i,3)=max(psdsh(ww,qq,1,dqsh,6,0)+psdsh(ww,qq,1,dqsh,6,1),
     *     psdh(ww,qq,1,0,0)+psdh(ww,qq,1,0,1)
     *    +psdsh(ww,qq,1,dqsh,0,0)+psdsh(ww,qq,1,dqsh,0,1)
     *    )/(4.*pi**2*alfe)*qq
        elseif(jmod.eq.10)then        !direct no emission
          ww=qq/xx
          ar(i,3)=(psdh(ww,qq,1,-1,0)+psdh(ww,qq,1,-1,1)
     *    +psdsh(ww,qq,1,dqsh,-1,0)+psdsh(ww,qq,1,dqsh,-1,1)
     *    )/(4.*pi**2*alfe)*qq
        elseif(jmod.eq.11)then        !singlet (light)
          ww=qq/xx
          ar(i,3)=(psdh(ww,qq,1,1,0)+psdh(ww,qq,1,1,1)
     *    +psdsh(ww,qq,1,dqsh,1,0)+psdsh(ww,qq,1,dqsh,1,1)
     *    )/(4.*pi**2*alfe)*qq
        elseif(jmod.eq.12)then        !non-singlet
          ww=qq/xx
          ar(i,3)=(psdh(ww,qq,1,2,0)+psdh(ww,qq,1,2,1)
     *    +psdsh(ww,qq,1,dqsh,2,0)+psdsh(ww,qq,1,dqsh,2,1)
     *    )/(4.*pi**2*alfe)*qq
        elseif(jmod.eq.13)then        !charm 
          ww=qq/xx
          ar(i,3)=(psdh(ww,qq,1,3,0)+psdh(ww,qq,1,3,1)
     *    +psdsh(ww,qq,1,dqsh,3,0)+psdsh(ww,qq,1,dqsh,3,1)
     *    )/(4.*pi**2*alfe)*qq
        elseif(jmod.eq.14)then        !point-like
          ww=qq/xx
          ar(i,3)=(psdh(ww,qq,1,4,0)+psdh(ww,qq,1,4,1)
     *    +psdsh(ww,qq,1,dqsh,4,0)+psdsh(ww,qq,1,dqsh,4,1)
     *    )/(4.*pi**2*alfe)*qq
        elseif(jmod.eq.15)then        !resolved (VDM EPOS)
          ww=qq/xx
          ar(i,3)=(psdh(ww,qq,1,5,0)+psdh(ww,qq,1,5,1)
     *    +psdsh(ww,qq,1,dqsh,5,0)+psdsh(ww,qq,1,dqsh,5,1)
     *    )/(4.*pi**2*alfe)*qq
        elseif(jmod.eq.16)then        !resolved (VDM Engel)
          ww=qq/xx
          ar(i,3)=(psdh(ww,qq,1,6,0)+psdh(ww,qq,1,6,1)
     *    +psdsh(ww,qq,1,dqsh,6,0)+psdsh(ww,qq,1,dqsh,6,1)
     *    )/(4.*pi**2*alfe)*qq

        elseif(jmod.eq.2)then    !all input EPOS
          ar(i,3)=(
     *(pspdfg(xx,q2cmin(1),qq,2,1)+2.*pspdfg(xx,q2cmin(1),qq,2,-1))/2.25 !u sea+val
     *+(pspdfg(xx,q2cmin(1),qq,2,2)+2.*pspdfg(xx,q2cmin(1),qq,2,-2))/9.   !d sea+val
     *+pspdfg(xx,q2cmin(1),qq,2,-3)*2./9.)                            !s sea
c         if(naflav.ge.4)ar(i,3)=ar(i,3)+pspdfg(xx,q2cmin(1),qq,2,-4)     !c sea
c     *    *2./2.25
        elseif(jmod.eq.21)then    !u val input EPOS
          ar(i,3)=pspdfg(xx,q2cmin(1),qq,2,1)/2.25    
c          ar(i,3)=pspdfg(xx,q2cmin(1),q2cmin(1),2,1)/2.25    
        elseif(jmod.eq.22)then    !d val input EPOS
          ar(i,3)=pspdfg(xx,q2cmin(1),qq,2,2)/9.
        elseif(jmod.eq.23)then    !sea input EPOS
          ar(i,3)=pspdfg(xx,q2cmin(1),qq,2,-1)*2./2.25
     *+(pspdfg(xx,q2cmin(1),qq,2,-2)+pspdfg(xx,q2cmin(1),qq,2,-3))*2./9.
        elseif(jmod.eq.20)then    !glu input EPOS
          ar(i,3)=pspdfg(xx,q2cmin(1),qq,2,0)

        elseif(jmod.eq.3)then    !grv94
          ar(i,3)=(psdfh4(xx,qq,0.,2,1)+2.*psdfh4(xx,qq,0.,2,-1))/2.25
     *    +(psdfh4(xx,qq,0.,2,2)+2.*psdfh4(xx,qq,0.,2,-2))/9.
     *    +2.*psdfh4(xx,qq,0.,2,-3)/9.  !
        elseif(jmod.eq.31)then    !u val grv94
          ar(i,3)=psdfh4(xx,qq,0.,2,1)/2.25
        elseif(jmod.eq.32)then    !d val grv94
          ar(i,3)=psdfh4(xx,qq,0.,2,2)/9.
        elseif(jmod.eq.33)then    !sea grv94
          ar(i,3)=(2.*psdfh4(xx,qq,0.,2,-1))/2.25
     *    +(2.*psdfh4(xx,qq,0.,2,-2))/9.
     *    +2.*psdfh4(xx,qq,0.,2,-3)/9.  !
        elseif(jmod.eq.30)then    !gluon grv94
          ar(i,3)=psdfh4(xx,qq,0.,2,0)

        elseif(jmod.eq.4)then         !ours new (missing charm contribution !)
          ii=2
          ar(i,3)=(pifpartone(1,dble(xx),qq,1,ii,0)        !u val+sea
     *         +   pifpartone(1,dble(xx),qq,-1,ii,0))/2.25 !ub sea
     *           +(pifpartone(1,dble(xx),qq,2,ii,0)        !d val+sea
     *         +   pifpartone(1,dble(xx),qq,-2,ii,0))/9.   !db sea
     *           +(pifpartone(1,dble(xx),qq,3,ii,0)        !s sea
     *         +   pifpartone(1,dble(xx),qq,-3,ii,0))/9.   !sb sea
          ar(i,3)=ar(i,3)+(pifpartone(1,dble(xx),qq,4,ii,0)
     *                    +pifpartone(1,dble(xx),qq,-4,ii,0))/2.25
          ar(i,3)=ar(i,3)+(pifpartone(1,dble(xx),qq,5,ii,0)
     *                    +pifpartone(1,dble(xx),qq,-5,ii,0))/9.
c          ww=qq/xx
c          ar(i,3)=max(ar(i,3),(psdgh(ww,qq,6,0)+psdgh(ww,qq,6,1)   !add resolved contrib.
c     *    )/(4.*pi**2*alfe)*qq)
        elseif(jmod.eq.444)then         !ours new (missing charm contribution !)
          ii=2
          ar(i,3)=(pifpartone(1,dble(xx),qq,4,ii,0)
     *            +pifpartone(1,dble(xx),qq,-4,ii,0))/2.25
        elseif(jmod.eq.41)then         !u val
          ar(i,3)=(pifpartone(1,dble(xx),qq,1,2,2))/2.25  !u val
c          ar(i,3)=(fzeroVal(dble(xx),qq,2,1))/2.25  !u val
        elseif(jmod.eq.42)then         !d val
          ar(i,3)=(pifpartone(1,dble(xx),qq,2,2,2))/9     !d val
c          ar(i,3)=(fzeroVal(dble(xx),qq,2,2))/9  !d val
        elseif(jmod.eq.43)then         !sea
          ar(i,3)=(pifpartone(1,dble(xx),qq,1,ii,1)        !u val+sea
     *         +   pifpartone(1,dble(xx),qq,-1,ii,1))/2.25 !ub sea
     *           +(pifpartone(1,dble(xx),qq,2,ii,1)        !d val+sea
     *         +   pifpartone(1,dble(xx),qq,-2,ii,1))/9.   !db sea
     *           +(pifpartone(1,dble(xx),qq,3,ii,1)        !s sea
     *         +   pifpartone(1,dble(xx),qq,-3,ii,1))/9.   !sb sea
          ar(i,3)=ar(i,3)+(pifpartone(1,dble(xx),qq,4,ii,1)
     *                    +pifpartone(1,dble(xx),qq,-4,ii,1))/2.25
          ar(i,3)=ar(i,3)+(pifpartone(1,dble(xx),qq,5,ii,1)
     *                    +pifpartone(1,dble(xx),qq,-5,ii,1))/9.
        elseif(jmod.eq.44)then         !c sea
          ar(i,3)=(pifpartone(1,dble(xx),qq, 4,2,1)
     *            +pifpartone(1,dble(xx),qq,-4,2,1))/2.25
        elseif(jmod.eq.45)then         !b sea
          ar(i,3)=(pifpartone(1,dble(xx),qq,5,2,1)
     *            +pifpartone(1,dble(xx),qq,-5,2,1))/9
        elseif(jmod.eq.46)then        !resolved (VDM Engel) (input from fit)
          ww=qq/xx
          ar(i,3)=(psdgh(ww,qq,6,0)+psdgh(ww,qq,6,1)
     *    )/(4.*pi**2*alfe)*qq
        elseif(jmod.eq.40)then         !sea
          ar(i,3)=pifpartone(1,dble(xx),qq,0,2,0)  !g

        elseif(jmod.eq.5)then        !evolution+matrix element (input from fit)
          ww=qq/xx
          ar(i,3)=max(psdgh(ww,qq,6,0)+psdgh(ww,qq,6,1),
     *    psdgh(ww,qq,0,0)+psdgh(ww,qq,0,1)
     *    )/(4.*pi**2*alfe)*qq
        elseif(jmod.eq.50)then        !direct no emission (input from fit)
          ww=qq/xx
          ar(i,3)=(psdgh(ww,qq,-1,0)+psdgh(ww,qq,-1,1)
     *    )/(4.*pi**2*alfe)*qq
        elseif(jmod.eq.51)then        !singlet light (input from fit)
          ww=qq/xx
          ar(i,3)=(psdgh(ww,qq,1,0)+psdgh(ww,qq,1,1)
     *    )/(4.*pi**2*alfe)*qq
        elseif(jmod.eq.52)then        !non-singlet (input from fit)
          ww=qq/xx
          ar(i,3)=(psdgh(ww,qq,2,0)+psdgh(ww,qq,2,1)
     *    )/(4.*pi**2*alfe)*qq
        elseif(jmod.eq.53)then        !charm (input from fit)
          ww=qq/xx
          ar(i,3)=(psdgh(ww,qq,3,0)+psdgh(ww,qq,3,1)
     *    )/(4.*pi**2*alfe)*qq
        elseif(jmod.eq.54)then        !point-like (input from fit)
          ww=qq/xx
          ar(i,3)=(psdgh(ww,qq,4,0)+psdgh(ww,qq,4,1)
     *    )/(4.*pi**2*alfe)*qq
        elseif(jmod.eq.55)then        !resolved (VDM EPOS) (input from fit)
          ww=qq/xx
          ar(i,3)=(psdgh(ww,qq,5,0)+psdgh(ww,qq,5,1)
     *    )/(4.*pi**2*alfe)*qq
        elseif(jmod.eq.56)then        !resolved (VDM Engel) (input from fit)
          ww=qq/xx
          ar(i,3)=(psdgh(ww,qq,6,0)+psdgh(ww,qq,6,1)
     *    )/(4.*pi**2*alfe)*qq

        elseif(jmod.eq.6)then    !cteq6
          ar(i,3)=(cteq6(xx,qq,0.,2,1)+2.*cteq6(xx,qq,0.,2,-1))/2.25
     *    +(cteq6(xx,qq,0.,2,2)+2.*cteq6(xx,qq,0.,2,-2))/9.
     *    +2.*cteq6(xx,qq,0.,2,-3)/9.  !
c     *    +cteq6(xx,qq,0.,2,-4)*2./2.25+cteq6(xx,qq,0.,2,-5)*2./9.
        elseif(jmod.eq.61)then    !u val cteq6
          ar(i,3)=cteq6(xx,qq,0.,2,1)/2.25
        elseif(jmod.eq.62)then    !d val cteq6
          ar(i,3)=cteq6(xx,qq,0.,2,2)/9.
        elseif(jmod.eq.63)then    !sea cteq6
          ar(i,3)=(2.*cteq6(xx,qq,0.,2,-1))/2.25
     *    +(2.*cteq6(xx,qq,0.,2,-2))/9.
     *    +2.*cteq6(xx,qq,0.,2,-3)/9.  !
     *    +cteq6(xx,qq,0.,2,-4)*2./2.25+cteq6(xx,qq,0.,2,-5)*2./9.
        elseif(jmod.eq.64)then    !sea cteq6
          ar(i,3)=(2.*cteq6(xx,qq,0.,2,-4))/2.25
        elseif(jmod.eq.65)then    !sea cteq6
          ar(i,3)=(2.*cteq6(xx,qq,0.,2,-5))/9.
        elseif(jmod.eq.60)then    !gluon cteq6
          ar(i,3)=cteq6(xx,qq,0.,2,0)
        elseif(jmod.eq.6143)then    !cteq14
          ar(i,3)=  ( xCTEQ(xx,qq, 1) + xCTEQ(xx,qq,-1) ) /2.25
     *            + ( xCTEQ(xx,qq, 2) + xCTEQ(xx,qq,-2) ) /9.
     *            +         2.*xCTEQ(xx,qq,-3)             /9.  
        elseif(jmod.eq.6145)then    !cteq14
          ar(i,3)=  ( xCTEQ(xx,qq, 1) + xCTEQ(xx,qq,-1) ) /2.25
     *            + ( xCTEQ(xx,qq, 2) + xCTEQ(xx,qq,-2) ) /9.
     *            +         2.*xCTEQ(xx,qq,-3)             /9.  
     *            +         2.*xCTEQ(xx,qq,-4)             /2.25
     *            +         2.*xCTEQ(xx,qq,-5)             /9.  
        elseif(jmod.eq.61444)then    !cteq14
          ar(i,3)=          2.*xCTEQ(xx,qq,-4)             /2.25
        elseif(jmod.eq.7)then    !psdpdf
          ar(i,3)=(psdpdf(xx,qq,0.,2,1)+2.*psdpdf(xx,qq,0.,2,-1))/2.25
     *    +(psdpdf(xx,qq,0.,2,2)+2.*psdpdf(xx,qq,0.,2,-2))/9.
     *    +2.*psdpdf(xx,qq,0.,2,-3)/9.  !
        elseif(jmod.eq.71)then    !u val psdpdf
          ar(i,3)=psdpdf(xx,qq,0.,2,1)/2.25
        elseif(jmod.eq.72)then    !d val psdpdf
          ar(i,3)=psdpdf(xx,qq,0.,2,2)/9.
        elseif(jmod.eq.73)then    !sea psdpdf
          ar(i,3)=(2.*psdpdf(xx,qq,0.,2,-1))/2.25
     *    +(2.*psdpdf(xx,qq,0.,2,-2))/9.
     *    +2.*psdpdf(xx,qq,0.,2,-3)/9.  !
        elseif(jmod.eq.70)then    !gluon psdpdf
          ar(i,3)=psdpdf(xx,qq,0.,2,0)

        elseif(jmod.eq.8)then    !grv08
          ar(i,3)=(grv08(xx,qq,0.,2,1)+2.*grv08(xx,qq,0.,2,-1))/2.25
     *    +(grv08(xx,qq,0.,2,2)+2.*grv08(xx,qq,0.,2,-2))/9.
     *    +2.*grv08(xx,qq,0.,2,-3)/9.  !
        elseif(jmod.eq.81)then    !u val grv08
          ar(i,3)=grv08(xx,qq,0.,2,1)/2.25
        elseif(jmod.eq.82)then    !d val grv08
          ar(i,3)=grv08(xx,qq,0.,2,2)/9.
        elseif(jmod.eq.83)then    !sea grv08
          ar(i,3)=(2.*grv08(xx,qq,0.,2,-1))/2.25
     *    +(2.*grv08(xx,qq,0.,2,-2))/9.
     *    +2.*grv08(xx,qq,0.,2,-3)/9.  !
        elseif(jmod.eq.80)then    !gluon grv08
          ar(i,3)=grv08(xx,qq,0.,2,0)

        elseif(jmod.eq.9)then    !F2 fractal fit
          ar(i,3)=fractalpdf(xx,qq,0.,2,0)
        endif
        ar(i,4)=0.
      enddo
      engy=engysave
      iclpro=iclprosave
      icltar=icltarsave
      return
      end

c------------------------------------------------------------------------
      function pspdfg(xx,qqs,qq,iclpro0,j)
c-----------------------------------------------------------------------
c pspdfg - parton distribution function - x*f(x,Q2)
c input from EPOS + QCD evolution
c qq  - virtuality scale
c qqs - initial virtuality for the input distributions
c iclpro0 - hadron class
c j   - parton type
c-----------------------------------------------------------------------
      double precision z,psuds,fzeroGlu,fzeroVal,FzeroQua
      common/ar3/    x1(7),a1(7)
#include "aaa.h"
#include "sem.h"

      dummy=iclpro0
       
      pspdfg=0.
c      if(qq.lt.q2ini)return

c      pspdfg=psdpdf(xx,qqs,0.,iclpro0,j)
cc      if(j.gt.0)pspdfg=pspdfg+psdpdf(xx,qqs,0.,iclpro0,-j)  !+sea contr.
      ii=1        !projectile side
      if(j.eq.0)then
        pspdfg=fzeroGlu(dble(xx),min(qq,qqs),ii) 
      elseif(j.gt.0)then
        pspdfg=fzeroVal(dble(xx),min(qq,qqs),ii,j)
      else
        pspdfg=FzeroQua(dble(xx),min(qq,qqs),ii,999)/naflav/2.
      endif
      if(qq.lt.q2sft)then      !no perturbative calculation
        goto 999
      elseif(qq.lt.qqs)then     !q2sft<qq<q2cmin (no evolution and more gluons)
        if(j.le.0)pspdfg=pspdfg*psuds(qq,j)/psuds(q2sft,j)
        goto 999
      endif
      pspdfg=pspdfg*psuds(qq,j)/psuds(qqs,j)

      xmin=xx/(1.-q2ini/qq)
      if(xmin.ge.1.)return

      dpd1=0.
      dpd2=0.
      xm=max(xmin,.3)
      do i=1,7         !numerical integration over zx
      do m=1,2
        zx=1.-(1.-xm)*(.5+(m-1.5)*x1(i))**.25
        z=xx/zx

        gl=fzeroGlu(dble(zx),qqs,ii)
        uv=fzeroVal(dble(zx),qqs,ii,1)
        dv=fzeroVal(dble(zx),qqs,ii,2)
        sea=FzeroQua(dble(zx),qqs,ii,999)/naflav/2.      !includes naflav and q + qb (*naflav*2)
        if(j.eq.0)then
          pj=gl 
        elseif(j.eq.1)then
          pj=uv
        elseif(j.eq.2)then
          pj=dv
        else
          pj=sea
        endif

        if(j.eq.0)then
          aks=psevi(qqs,qq,z,2,1)                  !quark contribution
          akg=psevi(qqs,qq,z,1,1)                  !gluon contribution
          akns=0.
        elseif(j.gt.0)then
          aks=0.
          akg=0.
          akns=psevi(qqs,qq,z,3,2)            !nonsinglet contribution
        else
          akg=psevi(qqs,qq,z,1,2)/naflav/2.         !gluon contribution
          akns=psevi(qqs,qq,z,3,2)            !nonsinglet contribution
          aks=(psevi(qqs,qq,z,2,2)-akns)/naflav/2.  !quark contribution
        endif

        fz=akg*gl+akns*pj+aks*(uv+dv+naflav*2.*sea)
c        if(j.gt.0)fz=fz+akns*psdpdf(zx,qqs,0.,iclpro0,-j)

        dpd1=dpd1+a1(i)*fz/zx**2/(1.-zx)**3
      enddo
      enddo
      dpd1=dpd1*(1.-xm)**4/8.*xx

      if(xm.gt.xmin)then
        do i=1,7         !numerical integration
        do m=1,2
          zx=xx+(xm-xx)*((xmin-xx)/(xm-xx))**(.5-(m-1.5)*x1(i))
          z=xx/zx

          gl=fzeroGlu(dble(zx),qqs,ii)
          uv=fzeroVal(dble(zx),qqs,ii,1)
          dv=fzeroVal(dble(zx),qqs,ii,2)
          sea=FzeroQua(dble(zx),qqs,ii,999)/naflav/2.      !includes naflav and q + qb (*naflav*2)
          if(j.eq.0)then
            pj=gl 
          elseif(j.eq.1)then
            pj=uv
          elseif(j.eq.2)then
            pj=dv
          else
            pj=sea
          endif

          if(j.eq.0)then
            aks=psevi(qqs,qq,z,2,1)                  !quark contribution
            akg=psevi(qqs,qq,z,1,1)                  !gluon contribution
            akns=0.
          elseif(j.gt.0)then
            aks=0.
            akg=0.
            akns=psevi(qqs,qq,z,3,2)             !nonsinglet contribution
          else
            akg=psevi(qqs,qq,z,1,2)/naflav/2.         !gluon contribution
            akns=psevi(qqs,qq,z,3,2)            !nonsinglet contribution
            aks=(psevi(qqs,qq,z,2,2)-akns)/naflav/2.  !quark contribution
          endif

          fz=akg*gl+akns*pj+aks*(uv+dv+naflav*2.*sea)
c          if(j.gt.0)fz=fz+akns*psdpdf(zx,qqs,0.,iclpro0,-j)

          dpd2=dpd2+a1(i)*fz*(1.-xx/zx)/zx
        enddo
        enddo
        dpd2=dpd2*log((xm-xx)/(xmin-xx))*.5*xx
      endif
      pspdfg=pspdfg+dpd2+dpd1
 999  return
      end

c------------------------------------------------------------------------
      function pspdfgfit(xx,qqs,qq,iclpro0,j)
c-----------------------------------------------------------------------
c pspdfg - parton distribution function - x*f(x,Q2)
c input from fit + QCD evolution
c qq  - virtuality scale
c qqs - initial virtuality for the input distributions
c iclpro0 - hadron class
c j   - parton type
c-----------------------------------------------------------------------
      double precision z,psuds
      common/ar3/    x1(7),a1(7)
#include "sem.h"

      pspdfgfit=0.
c      if(qq.lt.q2ini)return

      pspdfg=psdpdf(xx,qqs,0.,iclpro0,j)
c      if(j.gt.0)pspdfg=pspdfg+psdpdf(xx,qqs,0.,iclpro0,-j)  !+sea contr.
      pspdfg=pspdfg*psuds(qq,j)/psuds(qqs,j)

      xmin=xx/(1.-q2ini/qq)
      if(xmin.ge.1.)return

      dpd1=0.
      dpd2=0.
      xm=max(xmin,.3)
      do i=1,7         !numerical integration over zx
      do m=1,2
        zx=1.-(1.-xm)*(.5+(m-1.5)*x1(i))**.25
        z=xx/zx

        if(j.eq.0)then
          aks=psevi(qqs,qq,z,2,1)                  !quark contribution
          akg=psevi(qqs,qq,z,1,1)                  !gluon contribution
          akns=0.
        else
          akg=psevi(qqs,qq,z,1,2)/naflav/2.         !gluon contribution
          akns=psevi(qqs,qq,z,3,2)            !nonsinglet contribution
          aks=(psevi(qqs,qq,z,2,2)-akns)/naflav/2.  !quark contribution
        endif

        fz=akg*psdpdf(zx,qqs,0.,iclpro0,0)
     *  +akns*psdpdf(zx,qqs,0.,iclpro0,j)
     *  +aks*(psdpdf(zx,qqs,0.,iclpro0,1)+
     *  2.*psdpdf(zx,qqs,0.,iclpro0,-1)
     *  +psdpdf(zx,qqs,0.,iclpro0,2)+2.*psdpdf(zx,qqs,0.,iclpro0,-2)
     *  +2.*psdpdf(zx,qqs,0.,iclpro0,-3))
c        if(j.gt.0)fz=fz+akns*psdpdf(zx,qqs,0.,iclpro0,-j)

        dpd1=dpd1+a1(i)*fz/zx**2/(1.-zx)**3
      enddo
      enddo
      dpd1=dpd1*(1.-xm)**4/8.*xx

      if(xm.gt.xmin)then
        do i=1,7         !numerical integration
        do m=1,2
          zx=xx+(xm-xx)*((xmin-xx)/(xm-xx))**(.5-(m-1.5)*x1(i))
          z=xx/zx

          if(j.eq.0)then
            aks=psevi(qqs,qq,z,2,1)                  !quark contribution
            akg=psevi(qqs,qq,z,1,1)                  !gluon contribution
            akns=0.
          else
            akg=psevi(qqs,qq,z,1,2)/naflav/2.         !gluon contribution
            akns=psevi(qqs,qq,z,3,2)            !nonsinglet contribution
            aks=(psevi(qqs,qq,z,2,2)-akns)/naflav/2.  !quark contribution
          endif

          fz=akg*psdpdf(zx,qqs,0.,iclpro0,0)
     *    +akns*psdpdf(zx,qqs,0.,iclpro0,j)
     *    +aks*(psdpdf(zx,qqs,0.,iclpro0,1)
     *    +2.*psdpdf(zx,qqs,0.,iclpro0,-1)
     *    +psdpdf(zx,qqs,0.,iclpro0,2)+2.*psdpdf(zx,qqs,0.,iclpro0,-2)
     *    +2.*psdpdf(zx,qqs,0.,iclpro0,-3))
c          if(j.gt.0)fz=fz+akns*psdpdf(zx,qqs,0.,iclpro0,-j)

          dpd2=dpd2+a1(i)*fz*(1.-xx/zx)/zx
        enddo
        enddo
        dpd2=dpd2*log((xm-xx)/(xmin-xx))*.5*xx
      endif
      pspdfg=pspdfg+dpd2+dpd1
      return
      end

c------------------------------------------------------------------------
      function psdgh(s,qq,iqq,long)
c-----------------------------------------------------------------------
c psdgh - sigma_gamma*_p (PDF from fit)
c all contributions eq 10.59 and 10.60 Pys Rep with phi from PDF fit at q2dmin
c s - energy squared for the interaction (hadron-hadron),
c iqq - contribution type : 
c  -1 - direct no emission, 0 - all, 1 - direct singlet, 2 - direct non-singlet
c                           3 - direct charm, 4 - point like, 5 - resolved (VDM)
c                           6 - resolved (VDM Engel)
c long - 0 (Transverse (10.59) ) - 1 (Longitudinal (10.60) )
c-----------------------------------------------------------------------
      common/ar3/    x1(7),a1(7)
      common /cnsta/ pi,pii,hquer,prom,piom,ainfin
#include "sem.h"
      double precision psuds

      psdgh=0.

      if(iqq.eq.6)then
        psdgh=psdvdm(qq,s,long)
        return
      endif


      xd=qq/s
      qqs=min(qq,q2dmin)
      if(long.eq.0.and.iqq.le.0)then    !direct no emission
        psdgh=2.*(psdpdf(xd,qqs,0.,2,-1)+psdpdf(xd,qqs,0.,2,-2)+
     *  psdpdf(xd,qqs,0.,2,-3))/4.5
        if(qq.lt.q2ini)psdgh=psdgh*qq/q2ini !correction below q2ini
        psdgh=(psdpdf(xd,qqs,0.,2,1)/2.25+psdpdf(xd,qqs,0.,2,2)/9.
     *  +psdpdf(xd,qqs,0.,2,3)/9.+psdgh)
     *  *4.*pi**2*alfe/qq
        if(qq.ge.q2dmin)psdgh=psdgh*psuds(qq,1)/psuds(q2dmin,1)
      else
        psdgh=0.
      endif

      if(qq.le.q2ini)return
      dgh=0.
      if(long.eq.0)then
        s2min=qq/(1.-q2ini/qq)
      else
        s2min=4.*max(q2dmin,qcmass**2)+qq
        s2min=s2min/(1.-4.*q2ini/(s2min-qq))
      endif
      if(s.lt.s2min)return
      xmin=s2min/s

      if(xmin.lt.1..and.iqq.ge.0)then
        do i=1,7          !numerical integration over z1
        do m=1,2
          if(long.eq.0)then
            z1=qq/s+(xmin-qq/s)*((1.-qq/s)/(xmin-qq/s))
     *      **(.5+(m-1.5)*x1(i))
          else
            z1=.5*(1.+xmin+(2*m-3)*x1(i)*(1.-xmin))
          endif
          call psdint(z1*s,qq,sds,sdn,sdb,sdt,sdr,1,long)
          call psdint(z1*s,qq,sdsg,sdng,sdbg,sdtg,sdrg,0,long)
          if(iqq.eq.1)then              !direct charm (singlet=with emission + born=no emission)
            sdt=sds
            sdtg=sdsg+sdbg
            sdn=sdb
          elseif(iqq.eq.2)then          !direct non-singlet
            sdt=0.
            sdtg=0.
          elseif(iqq.eq.3)then          !born charm (no emission)
            sdt=0.
            sdtg=sdbg
            sdn=sdb
          elseif(iqq.eq.4)then          !direct light (singlet + point-light)
            sdt=max(0.,sdt-sds-sdr)
            sdtg=max(0.,sdtg-sdsg-sdrg)
            sdn=0.
          elseif(iqq.eq.5)then          !resolved (VDM)
            sdt=sdr
            sdtg=sdrg
            sdn=0.
          else                  !all
            sdt=max(0.,sdt-sdr) !don't use VDM from psdint but psdvdm
            sdtg=max(0.,sdtg-sdrg)     !don't use VDM from psdint but psdvdm
            sdtg=sdtg+sdbg !max(sdtg,sdbg) !charm
            sdn=sdn+sdb !max(sdn,sdb)    !charm
          endif
          tu=psdpdf(z1,q2dmin,0.,2,1)
          td=psdpdf(z1,q2dmin,0.,2,2)
          ts=psdpdf(z1,q2dmin,0.,2,3)
          tg=psdpdf(z1,q2dmin,0.,2,0)
          tsea=2.*(psdpdf(z1,q2dmin,0.,2,-1)+psdpdf(z1,q2dmin,0.,2,-2)
     *    +psdpdf(z1,q2dmin,0.,2,-3))
          gy=sdn*(tu/2.25+td/9.+ts/9.+tsea/4.5)+sdtg*tg/4.5
     *    +sdt*(tu+td+ts+tsea)/4.5
          dgh=dgh+a1(i)*gy*(1.-qq/s/z1)
        enddo
        enddo
        dgh=dgh*log((1.-qq/s)/(xmin-qq/s))*.5
      endif
      psdgh=psdgh+dgh
c      if(iqq.eq.5)print *,qq,sqrt(s),sdr,psdgh*0.0389*10000.
      return
      end
c------------------------------------------------------------------------
      function psdsh1(s,qq,ipt,dqsh,iqq,long)
c-----------------------------------------------------------------------
c psdsh - semihard gamma-p interaction  - sigma_gamma*_p (sea) without screening
c all contributions eq 10.59 and 10.60 Pys Rep with phi(sea)
c s - energy squared for the interaction
c qq - Q2
c ipt - proj/Tar,
c long - 0 (Transverse (10.59) ) - 1 (Longitudinal (10.60) )
c-----------------------------------------------------------------------
      common /ar3/    x1(7),a1(7)
#include "aaa.h"
#include "sem.h"
      double precision psuds,FzeroQuaZZ,fzeroGluZZ


      psdsh1=0.
      dqsh=0.
      dvdm=0.
      iscreensave=iscreen

      if(iqq.eq.6)then
        psdsh1=psdvdm(qq,s,long)
        return
      endif

      xd=qq/s
      if(long.eq.0.and.(idisco.eq.0.or.idisco.eq.1)
     *               .and.iqq.le.0)then
        dqsh=FzeroQuaZZ(dble(xd),ipt)/xd**dels(ipt)/4.5
     *  *4.*pi**2*alfe/qq
c     *  *gamhad(iclpro0)*ffrr*4.*pi   !now in FzeroQuaZZ
        if(qq.lt.q2ini)dqsh=dqsh*qq/q2ini !correction below q2ini
        if(qq.ge.q2cmin(ipt))dqsh=dqsh*psuds(qq,1)/psuds(q2cmin(ipt),1)
      else
        dqsh=0.
      endif
      psdsh1=dqsh
      if(qq.le.q2ini)return
      if(long.eq.0)then
        s2min=qq/(1.-q2ini/qq)
      else
        s2min=qq+4.*max(q2cmin(ipt),qcmass**2)
      endif
      if(s.lt.s2min)return
      xmin=s2min/s
      xmin=xmin**(delh-dels(ipt))
      dsh=0.
      if(xmin.lt.1..and.iqq.ge.0)then
c numerical integration over z1
        do i=1,7
        do m=1,2
          z1=(.5*(1.+xmin-(2*m-3)*x1(i)*(1.-xmin)))**(1./
     *    (delh-dels(ipt)))
          call psdint(z1*s,qq,sdsg,sdng,sdbg,sdtg,sdrg,0,long)
          call psdint(z1*s,qq,sdsq,sdnq,sdbq,sdtq,sdrq,1,long)
          if(iqq.eq.1)then              !direct charm (singlet=with emission + born=no emission)
            sdtq=sdsq
            sdtg=sdsg+sdbg
            sdnq=sdbq
          elseif(iqq.eq.2)then          !direct non-singlet
            sdtq=0.
            sdtg=0.
          elseif(iqq.eq.3)then          !born charm (no emission)
            sdtq=0.
            sdtg=sdbg
            sdnq=sdbq
          elseif(iqq.eq.4)then          !direct light (singlet + point-light)
            sdtq=max(0.,sdtq-sdsq-sdrq)
            sdtg=max(0.,sdtg-sdsg-sdrg)
            sdnq=0.
          elseif(iqq.eq.5)then          !resolved (VDM)
            sdtq=sdrq
            sdtg=sdrg
            sdnq=0.
          else                          !all
            sdtq=max(0.,sdtq-sdrq)     !don't use VDM from psdint but psdvdm
            sdtg=max(0.,sdtg-sdrg)     !don't use VDM from psdint but psdvdm
            sdtg=sdtg+sdbg !max(sdtg,sdbg)        !charm
            sdnq=sdnq+sdbq !max(sdnq,sdbq)        !charm
          endif
          dsh=dsh+a1(i)/z1**delh*(sdtg*fzeroGluZZ(dble(z1),ipt)
     *    +(sdtq+sdnq)*FzeroQuaZZ(dble(z1),ipt))
c simplified in equation: 
c                                 *z1**(-1-dels(ipt))  !in front of FzeroQuaZZ and fzeroGluZZ
c                                 *z1/z1**dble(delh-dels(ipt))  !Jacobian
        enddo
        enddo
        dsh=dsh*(1.-xmin)/(delh-dels(ipt))/2.
      endif
      psdsh1=dqsh+dsh/4.5!*ffrr*4.*pi*gamhad(iclpro0)  !now in FzeroQuaZZ
c      psdsh1=dsh/4.5!*ffrr*4.*pi*gamhad(iclpro0)  !now in FzeroQuaZZ

      iscreen=iscreensave
      return
      end

c------------------------------------------------------------------------
      function psdsh(s,qq,ipt,dqsh,iqq,long)
c-----------------------------------------------------------------------
c psdsh - semihard gamma-p interaction  - sigma_gamma*_p (sea) with screening
c all contributions eq 10.59 and 10.60 Pys Rep with phi(sea)
c s - energy squared for the interaction
c qq - Q2
c ipt - proj/tar,
c long - 0 (Transverse (10.59) ) - 1 (Longitudinal (10.60) )
c-----------------------------------------------------------------------
      common /ar3/    x1(7),a1(7)
#include "aaa.h"
#include "sem.h"
      double precision psuds,FzeroQua,fzeroGlu


      psdsh=0.
      dqsh=0.
      dvdm=0.

      if(iqq.eq.6)then
        psdsh=0. !psdvdm(qq,s,long)
        return
      endif


      xd=qq/s
      if(long.eq.0.and.(idisco.eq.0.or.idisco.eq.1)
     *               .and.iqq.le.0)then
        dqsh=FzeroQua(dble(xd),qq,1,999)/4.5
     *  *4.*pi**2*alfe/qq
c     *  *gamhad(iclpro0)*ffrr*4.*pi   !now in FzeroQuaZZ
        if(qq.lt.q2ini)dqsh=dqsh*qq/q2ini !correction below q2ini
        if(qq.ge.q2cmin(ipt))dqsh=dqsh*psuds(qq,1)/psuds(q2cmin(ipt),1)
      else
        dqsh=0.
      endif
      psdsh=dqsh
      if(qq.le.q2ini)return
      if(long.eq.0)then
        s2min=qq/(1.-q2ini/qq)
      else
        s2min=qq+4.*max(q2cmin(ipt),qcmass**2)
      endif
      if(s.lt.s2min)return

      xmin=s2min/s
      xmin=xmin**(delh-dels(ipt))
      dsh=0.
      if(xmin.lt.1..and.iqq.ge.0)then
c numerical integration over z1
        do i=1,7
        do m=1,2
          z1=(.5*(1.+xmin-(2*m-3)*x1(i)*(1.-xmin)))**(1./
     *    (delh-dels(ipt)))
          call psdint(z1*s,qq,sdsg,sdng,sdbg,sdtg,sdrg,0,long)
          call psdint(z1*s,qq,sdsq,sdnq,sdbq,sdtq,sdrq,1,long)
          if(iqq.eq.1)then              !direct charm (singlet=with emission + born=no emission)
            sdtq=sdsq
            sdtg=sdsg+sdbg
            sdnq=sdbq
          elseif(iqq.eq.2)then          !direct non-singlet
            sdtq=0.
            sdtg=0.
          elseif(iqq.eq.3)then          !born charm (no emission)
            sdtq=0.
            sdtg=sdbg
            sdnq=sdbq
          elseif(iqq.eq.4)then          !direct light (singlet + point-light)
            sdtq=max(0.,sdtq-sdsq-sdrq)
            sdtg=max(0.,sdtg-sdsg-sdrg)
            sdnq=0.
          elseif(iqq.eq.5)then          !resolved (VDM)
            sdtq=sdrq
            sdtg=sdrg
            sdnq=0.
          else                          !all
            sdtq=max(0.,sdtq-sdrq)     !don't use VDM from psdint but psdvdm
            sdtg=max(0.,sdtg-sdrg)     !don't use VDM from psdint but psdvdm
            sdtg=sdtg+sdbg !max(sdtg,sdbg)        !charm
            sdnq=sdnq+sdbq !max(sdnq,sdbq)        !charm
          endif
          dsh=dsh+a1(i)/z1**(delh-dels(ipt))*(sdtg
     *    *fzeroGlu(dble(z1),qq,1)+(sdtq+sdnq)
     .    *FzeroQua(dble(z1),qq,1,999))
c simplified in equation: 
c                                 *z1**(-1)  !in front of FzeroQua and fzeroGlu
c                                 *z1/z1**dble(delh-dels)  !Jacobian
        enddo
        enddo
        dsh=dsh*(1.-xmin)/(delh-dels(ipt))/2.
      endif
      psdsh=dqsh+dsh/4.5!*ffrr*4.*pi*gamhad(iclpro0)  !now in FzeroQua
c      psdsh=dsh/4.5
      return
      end


c------------------------------------------------------------------------
      function psdh(s,qq,ipt,iqq,long)
c-----------------------------------------------------------------------
c pshard - hard gam-quark interaction cross-section - sigma_gamma*_p (val)
c all contributions eq 10.59 and 10.60 Pys Rep with phi(val)
c ipt - proj/tar
c long - 0 (Transverse (10.59) ) - 1 (Longitudinal (10.60) )
c-----------------------------------------------------------------------
      common /ar3/   x1(7),a1(7)
#include "sem.h"
#include "aaa.h"
      double precision psuds

      psdh=0.
      if(iqq.eq.6)return

      xd=qq/s
      qqs=min(qq,q2cmin(ipt))
      if(long.eq.0.and.(idisco.eq.0.or.idisco.eq.1)
     *                .and.iqq.le.0)then
        psdh=(psdpdf(xd,qqs,0.,2,1)/2.25+
     *  psdpdf(xd,qqs,0.,2,2)/9.)
     *  *4.*pi**2*alfe/qq
        if(qq.gt.qqs)psdh=psdh*psuds(qq,1)/psuds(qqs,1)
      else
        psdh=0.
      endif
      if(qq.le.q2ini)return
      dh=0.
      if(long.eq.0)then
        s2min=qq/(1.-q2ini/qq)
      else
        s2min=4.*max(q2cmin(ipt),qcmass**2)+qq
        s2min=s2min/(1.-4.*q2ini/(s2min-qq))
      endif
      if(s.lt.s2min)return
      xmin=s2min/s
      if(xmin.lt.1..and.iqq.ge.0)then
        do i=1,7          !numerical integration over z1
        do m=1,2
          if(long.eq.0)then
            z1=qq/s+(xmin-qq/s)*((1.-qq/s)/(xmin-qq/s))
     *      **(.5+(m-1.5)*x1(i))
          else
            z1=.5*(1.+xmin+(2*m-3)*x1(i)*(1.-xmin))
          endif
          call psdint(z1*s,qq,sds,sdn,sdb,sdt,sdr,1,long)
          if(iqq.eq.1)then              !direct charm (singlet=with emission)
            sdt=sds
            sdn=0.
          elseif(iqq.eq.2)then          !direct non-singlet
            sdt=0.
          elseif(iqq.eq.3)then          !born charm (no emission)
            sdt=0.
            sdn=0.                      !no charm production from valence quark
          elseif(iqq.eq.4)then          !direct light (singlet + point-light)
            sdt=max(0.,sdt-sds-sdr)
            sdn=0.
          elseif(iqq.eq.5)then          !resolved (VDM)
            sdt=sdr
            sdn=0.
          else                          !all
            sdt=max(0.,sdt-sdr)     !don't use VDM from psdint but psdvdm
          endif
          tu=psdpdf(z1,qqs,0.,2,1)
          td=psdpdf(z1,qqs,0.,2,2)
          gy=sdt*(tu+td)/4.5+sdn*(tu/2.25+td/9.)
          if(long.eq.0)then
            gy=gy*(1.-qq/s/z1)
          else
            gy=gy/z1
          endif
          dh=dh+a1(i)*gy
        enddo
        enddo
        if(long.eq.0)then
          dh=dh*log((1.-qq/s)/(xmin-qq/s))*.5
        else
          dh=dh*(1.-xmin)*.5
        endif
      endif
      psdh=psdh+dh
      return
      end


