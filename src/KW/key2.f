c
c  This file is part of EPOS4
c  Copyright (C) 2022 research institutions and authors (See CREDITS file)
c  This file is distributed under the terms of the GNU General Public License version 3 or later
c  (See COPYING file for the text of the licence)
c

c-----------------------------------------------------------------------
      subroutine setParamEdep(q2imin)  !parameter E-dependence
c-----------------------------------------------------------------------
c     the soft contribution should not be much lower than the hard one,
c     even for large Q2 because it will increase the low mass hard jet
c     too much (too many low pt particles).
c     so alpq should be set such that the ratio/hard soft doesn't change
c     much with Q2s (otherwise it distort scaling) -> now fixed in WomTy
c     alpq should be large enough not to change this but should be present
c     to have some amplitude at the threshold increasing with Q2s
c-----------------------------------------------------------------------
#include "aaa.h"
#include "sem.h"
#include "par.h"
      dimension q2imin(2)
      dummy=q2imin(1)
      gamzz=1.
c      return
      if(iscreen.eq.0)return
c      dummy=q2imin(2)
c      x=log(engy/200.) !20
c      !~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
c      ! gamzz : 3077 method revived
c      !          also ProPoTy changed to get rid of soft if gamzz=0
c      !~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
c      w= 1 - x
      alpq=alpsft  !0.45 !0.75 !(0.33+alppom)/2. !(alppar+alppom)/2. !(alppar+dels(1)+1.)/2.
      w=min(1.000,q2nmin/q2imin(1))**alpq
c      alpq=(alppar+dels(2)+1.)/2.
      w=w*min(1.000,q2nmin/q2imin(2))**alpq
c      w=min(1.000,SoftSat(xminDf,q2sft,1)*SoftSat(xminDf,q2sft,2))
c     .            *q2sft/q2cmin(1)*q2sft/q2cmin(2))
c      w=xminDf**(0.5*(dels(1)+dels(2))-alppom+1.)
c      w=w*exp(-max(0.,q2imin(2)-q2sft))
c      print *,w,SoftAtt(xminDf,q2sft,1),SoftSat(xminDf,q2sft,1)
      w=max(1e-27,w)
      gamzz=w       !soft suppression
      end

c-----------------------------------------------------------------------
      subroutine iniq2mn  !Define nodes for Q2s tabulation
c-----------------------------------------------------------------------
#include "aaa.h"
#include "sem.h"
#include "tab.h"
      qzzz=q2zmin*1.5 !KW1908: cross section calculations are problematic for  values < 1.5
      delq=1.5
      do n=1,maxq2mn
        q2mnval(n) = qzzz*delq**(n-1)
        write(q2mnch(n),'(i3)')n
        if(n.lt.10)q2mnch(n)(2:2)='0'
        if(n.lt.100)q2mnch(n)(1:1)='0'
      enddo
      end
c-----------------------------------------------------------------------
      subroutine ipoli(q,nq,frac)  !prepares interpolation  for above table
c-----------------------------------------------------------------------
#include "aaa.h"
#include "tab.h"
      qmi=log(q2mnval(1))
      delq=log(q2mnval(2))-qmi
      xx=(q-qmi)/delq +1
      nq=int(xx)
      nq=max(nq,1)
      nq=min(nq,maxq2mn-1)
      frac=min(1., xx-nq ) !more stable than free extrapolation
      !w(0)=1.-frac w(1)=frac
      end
c-----------------------------------------------------------------------
      subroutine ipoli3(q,nq,frac)  
c-----------------------------------------------------------------------
c prepares interpolation  for above table
c-----------------------------------------------------------------------
#include "aaa.h"
#include "tab.h"
      qmi=log(q2mnval(1))
      delq=log(q2mnval(2))-qmi
      xx=(q-qmi)/delq +1
      nq=int(xx)
      nq=max(nq,1)
      nq=min(nq,maxq2mn-2)
      frac=xx-nq 
      end

c-----------------------------------------------------------------------
      double precision function SoftSat(xx,qq,iflav,ipt)
c-----------------------------------------------------------------------
c Formula set to reproduce F2 at very low Q2
c-----------------------------------------------------------------------
#include "aaa.h"
#include "sem.h"
      double precision xx,fglu

      SoftSat=1d0
c      return
      
c      alpq=(dels(ipt)+1.-betff(ipt))/3.
c      SoftAtt=1.+(min(1.,qq/q2cmin(ipt))**alpq-1.)*exp(-500.0d0*xx)
      ss=q2cmin(ipt)*(1.-xx)/xx+1.
      ssmx=qq*(1.-xx)/xx+1.
      SoftSat=(1.d0+dble(min(1.,ssmx/ss)-1.))**zoeinc !*exp(-100.0d0*xx)
c     .       /xx**dble((betpom(ipt)/betpomi)**zzsoft-1.)
c     .       *(1.d0+dble(min(1.,400./sqrt(ss))-1.))  !to add saturation at very small x (but doesn't really change the inclusive cross-section (400. set not see the effect in measured F2)
c      print *,xx,qq,ss,ssmx,q2cmin(ipt)
      if(zzsoft.ge.0.)then
        fglu=dble(glusea/max(1.,q2cmin(ipt)/qq)**zzsoft)
      else
        fglu=dble(1.+(glusea-1.)/max(1.,q2cmin(ipt)/qq)**zzsoft)
      endif
      if(iflav.eq.0)then
        SoftSat=SoftSat/dble(1.-glusea)*(1d0-fglu)
      else
        SoftSat=SoftSat/dble(glusea)*fglu
      endif
c      print *,xx,qq,ss,ssmx,SoftSat
      return
      end

c-----------------------------------------------------------------------
      function XsectionSat(xx,qq,ipt)  !cross-section qq dependence
c-----------------------------------------------------------------------
c Formula modify born xsection in saturated mode.
c-----------------------------------------------------------------------
#include "aaa.h"
#include "sem.h"
      double precision xx
      dmy=xx 

c increase below Q2s only
      
c      XsectionSat=(qq/q2cmin(ipt)+max(0.,q2cmin(ipt)-qq)
c     .           *max(0d0,xx-1d0))**factjpsi
c      XsectionSat=(1.+max(0.,q2cmin(ipt)-qq))**factjpsi
      XsectionSat=(q2cmin(ipt)/qq)**factjpsi
c      XsectionSat=max(1.e-5,(q2cmin(ipt)/qq-1.))**factjpsi
      return
      end

c-----------------------------------------------------------------------
      function fglusea(xx,qq,ipt)  !glusea qq dependence
c-----------------------------------------------------------------------
c Formula set to reproduce F2 at very low Q2
c-----------------------------------------------------------------------
#include "aaa.h"
#include "sem.h"
      double precision xx
      dmy=ipt
      dmy=xx 
      dmy=qq

c increase below Q2s only (less gluons (saturation)
      
      fglusea=glusea!*SoftSat(xx,qq,ipt)       
c      alpq=(0.33+dels(ipt)+1.)/2.         !alppar set on F2 for Q2<Q20
c      fglusea=1.+(glusea-1.)/(q2cmin(ipt)/qq)**zzsoft
c      fglusea=1.+(glusea-1.)/(1.+max(0.,q2cmin(ipt)-qq))**0.14
c     .       *(1.+(min(1.,qq/q2cmin(ipt))**alpq-1.)*exp(-30.0d0*xx))

      return
      end

c-----------------------------------------------------------------------
      subroutine fieloss(enp,dens,r12,ppj,fiel,jcone)
c-----------------------------------------------------------------------
      ! energy loss per unit lenth, for initial estimate
      ! enp = pt of segment
      ! dens = number of segments per cell (in jintop)
      ! ppj = momentum of jet in case of cone
      ! r12 = distance to the jet
      !--------------------------------------------------
      ! fully equibrated plasma
      ! aE/dL ~ sqrt(T**3*E) (BDMPS2008)
      ! -> dE/dL ~ dens^3/8 * sqrt(E)
      !--------------------------------------------------
#include "aaa.h"
      double precision enp
      dummy=ppj
      jcone=0
      call felossget6(fconepp,fconeaa,fbulkpp,fbulkaa,rxcone,Zx)
      Z=max(0.,ng1evt-2.)/Zx   !ng1evt=Npart
      call  getSystemType(isys,amassMax,amassAsy)
      if(isys.le.2)call getRnpom(0.5,Z) !pp,pA,dA
      Z=min(1.,Z)
      zex=Z
      if(r12.lt.rxcone)then 
        jcone=1    
      endif              
      elex=1.0
      gel=1
      if(jcone.eq.0)then
        q1=fbulkpp
        q2=fbulkaa
        fel=q1+zex*(q2-q1)
      else !-----> jet
        elex=0.375
        q1=fconepp
        q2=fconeaa
        fel=q1+zex*(q2-q1)
        gel=enp**0.5
      endif
      fiel = 0.25 * gel * fel * dens**elex
c      print *,jcone,enp,dens,fiel
      end


c-----------------------------------------------------------------------
      subroutine setfacpos(k)
c-----------------------------------------------------------------------
      ! Defines facpos which determines spatial parton distribution.
      ! Bigger value gives smoother energy density distribution. 
      ! Subroutine conigr needs probably improvement !!!
      !-------------------------------------------------------
#include "aaa.h"
#include "ems.h"
#include "par.h"
      common/geom/rmproj,rmtarg,bmax,bkmx
      dummy=k 
      Z=max(0.,ng1evt-2.)/400.   !ng1evt=Npart  
      call  getSystemType(isys,amassMax,amassAsy)
      if(isys.le.2)call getRnpom1(0.5,Z) !pp,pA,dA
      Z=min(1.,Z)
      if(iscreen.eq.1)then
        facpos(1)=facposf+Z*(facposz-facposf)
        facpos(2)=facposf+Z*(facposz-facposf)
      endif
      end

c-------------------------------------------------------------
      subroutine getZpar(zpav,ztav) !called in TP/xsc
c------------------------------------------------------------- 
c called in TP/xcs.f:getZparEpsiK
c output: 
c        zparpro(i,k)  with i = -1, 0 , 1, 2         
c        zpartar(i,k)  with i = -1, 0 , 1, 2     
c        zparpnx(j,k)  with j =     0 , 1, 2    \ Meaning? 
c        zpartnx(j,k)  with j =     0 , 1, 2    /  Not really used!!
c--------------------------------------------------------------  
ctp 20180307 : there is 3 "flavors" of Zpart 
      !0=diff (no staturation (large)): use for diffraction G contributions
      !1=non-diff for x (saturation which could depend on b value and nuclear
      !2=non-diff for s  density effect): use for non-diff G contributions
ckw gammaV stuff removed
c-------------------------------------------------------------      
#include "aaa.h"
#include "sem.h"
#include "ems.h"
#include "par.h"
      common /cncl/xproj(mamx),yproj(mamx),zproj(mamx)
     *            ,xtarg(mamx),ytarg(mamx),ztarg(mamx)
      dimension ZpartPro(0:2),ZpartTar(0:2)
      common/cnparticip/jproj(2,mamx),jtarg(2,mamx),efluct(6,mamx)
      real eflu(2),efly(2),efli(2)

      call screen1paramget6(efluPP,eflyPP,efliPP,efluPA,eflyPA,efliPA)
      sy=engy*engy

      b2x=b2xscr
      alpfom=0.                 !alpfomi*fegypp
      bcut=0.                   !zbcut*b2x
      ztav=0.
      zpav=0.
      rho=0.
      rhop=0.
      rhot=0.

      do k=1,koll

        gammaV(k)=1.d0
        zparpro(-1,k)=0.         
        zparpro(0,k)=0.         
        zparpro(1,k)=0.         
        zparpro(2,k)=0.         
        zpartar(-1,k)=0.         
        zpartar(0,k)=0.         
        zpartar(1,k)=0.         
        zpartar(2,k)=0.         
        zparpnx(0,k)=0.
        zparpnx(1,k)=0.
        zparpnx(2,k)=0.
        zpartnx(0,k)=0.
        zpartnx(1,k)=0.
        zpartnx(2,k)=0.

        if(iscreen.eq.1)then

          fga1=0
          fga2=0
          call Zparticip(k,ZpartPro,ZpartTar)  !ip-it not counted
          ip=iproj(k)
          it=itarg(k)
          isgp=jproj(2,ip)
          isgt=jtarg(2,it)
          do i=1,2
            eflu(i)=0
            efly(i)=0
            efli(i)=0
          enddo
          call  getSystemType(isys,amassMax,amassAsy)
          if(iokoll.ne.0)then
            eflu(1)=efluPA
            efly(1)=eflyPA 
            efli(1)=efliPA
          elseif(maproj.eq.1.and.matarg.eq.1)then !---pp----
            do i=1,2
            eflu(i)=efluPP
            efly(i)=eflyPP
            efli(i)=efliPP
            enddo
          elseif(maproj.eq.1)then     !---pA---
            eflu(1)=efluPA 
            efly(1)=eflyPA 
            efli(1)=efliPA
          elseif(matarg.eq.1)then     !---Ap---
            eflu(2)=efluPA  
            efly(2)=eflyPA  
            efli(2)=efliPA
          endif
          zgl=ng1evt

        !-------targ partons seen by proj changing targ exponent ------

c Z for Q2s
          zkp=ZpartPro(1)
          zkt=ZpartTar(1)
          zpartar(-1,k)
     .    =Zpair(zkt,zkp,zgl,eflu(1),efly(1),efli(1),bk(k),-2)
c Z for diff
          zkp=ZpartPro(0)
          zkt=ZpartTar(0)
          zpartnx(0,k)=zkt             !Z without current pair for diff
          zpartar(0,k)
     .    =Zpair(zkt,zkp,zgl,eflu(1),efly(1),efli(1),bk(k),0)    
c Z for xp and xm dependence of G
          zkp=ZpartPro(1)
          zkt=ZpartTar(1)
          zpartnx(1,k)=zkt              !Z without current pair
          zpartar(1,k)
     .    =Zpair(zkt,zkp,zgl,eflu(1),efly(1),efli(1),bk(k),1)     
c Z for normalization of G 
          !zkt=ZpartTar(2)
          zpartnx(2,k)=ZpartTar(2)              !maximum b of nuclear Z
          zpartar(2,k)
     .    =Zpair(zkt,zkp,zgl,eflu(1),efly(1),efli(1),bk(k),2)  
c          print *,bk(k),zpartar(2,k),zkt,zkp

          ztav=ztav+zpartar(2,k)

         !------proj partons seen by targ changing proj exponent------

c          rho=rhop
c Z for Q2s
          zkp=ZpartPro(1)
          zkt=ZpartTar(1)
          zparpro(-1,k)
     .    =Zpair(zkp,zkt,zgl,eflu(2),efly(2),efli(2),bk(k),-2)
c Z for diff
          zkp=ZpartPro(0)
          zkt=ZpartTar(0)
          zparpnx(0,k)=zkp                  !Z without current pair
          zparpro(0,k)
     .    =Zpair(zkp,zkt,zgl,eflu(2),efly(2),efli(2),bk(k),0)
c Z for xp and xm dependence of G
          zkp=ZpartPro(1)
          zkt=ZpartTar(1)
          zparpnx(1,k)=zkp                 !Z without current pair
          zparpro(1,k)
     .    =Zpair(zkp,zkt,zgl,eflu(2),efly(2),efli(2),bk(k),1)     
c Z for normalization of G 
          !zkp=ZpartPro(2)
          zparpnx(2,k)=ZpartPro(2)        !maximum b of nuclear Z
          zparpro(2,k)
     .    =Zpair(zkp,zkt,zgl,eflu(2),efly(2),efli(2),bk(k),2)   

          zpav=zpav+zparpro(2,k)

          if(lproj(ip).gt.1.or.ltarg(it).gt.1)then
            Zpart=ZpartPro(0)+ZpartTar(0)
          endif

        endif                   ! no screening

      enddo

      end

c-----------------------------------------------------------------------
      subroutine  Zparticip(k,ZpartPro,ZpartTar)
c-----------------------------------------------------------------------
      !  returns the z's which determine the saturation scales
      !      attention: ip-it not counted
      !----------------------------------------------------------
#include "aaa.h"
#include "ems.h"
      double precision plc,s,drangen,PhiExpo
      dimension ZpartPro(0:2),ZpartTar(0:2)
      common/cems5/plc,s
      common/geom/rmproj,rmtarg,bmax,bkmx
      common /cncl/xproj(mamx),yproj(mamx),zproj(mamx)
     *            ,xtarg(mamx),ytarg(mamx),ztarg(mamx)

      ip=iproj(k)
      it=itarg(k)
      ZpartPro(0)=0
      ZpartTar(0)=0
      ZpartPro(1)=0
      ZpartTar(1)=0
      ZpartPro(2)=bkmx !*0.5/facposz   !affects facpos, to change it for pp, change 0.5 here 
      ZpartTar(2)=bkmx !*0.5/facposz
c      jproj(3,ip)=1
c      jtarg(3,it)=1
      iomegasave=iomega
      iomega=2  !use only non-diff amplitude to define enhanced diagram profile function
      call DefXminDf(s)
      !....... targ nucleons hit by proj -> ZpartTar
      if(lproj(ip).gt.0)then  
        do li=1,lproj(ip)
          kk=kproj(ip,li)
          itt=itarg(kk)
          if(itt.ne.it)then  !attention: ip-it not counted
            bij=
     .      sqrt((xtarg(it)-xtarg(itt))**2+(ytarg(it)-ytarg(itt))**2)
            rho=0. !rhop+conrho(2,rtarg)
c Diff
            if(bij.le.zdfbmx*bkmx)then
              w=Zpair(rho,0.,0.,0.,0.,0.,bk(kk),0)
              ZpartTar(0)=ZpartTar(0)+w
            endif
c Saturated (standard)
            if(drangen(dble(bij)).lt.
     .         1d0-Phiexpo(rho,rho,zbrmax,1d0,1d0,sngl(s),bij))then 
              w=Zpair(rho,0.,0.,0.,0.,0.,bk(kk),-1)
              ZpartTar(1)=ZpartTar(1)+w
c              jtarg(3,itt)=1
              ZpartTar(2)=max(ZpartTar(2),bk(kk))
            endif
          endif
        enddo
      endif
      !...........proj nucleons hit by targ -> ZpartPro
      if(ltarg(it).gt.0)then
        do li=1,ltarg(it)
          kk=ktarg(it,li)
          ipp=iproj(kk)
          if(ipp.ne.ip)then !attention: ip-it not counted
            bij=
     .      sqrt((xproj(ip)-xproj(ipp))**2+(yproj(ip)-yproj(ipp))**2)
            rho=0. !rhot+conrho(1,rproj)
c Diff
            if(bij.le.zdfbmx*bkmx)then
              w=Zpair(rho,0.,0.,0.,0.,0.,bk(kk),0)
              ZpartPro(0)=ZpartPro(0)+w
            endif
c Saturated (standard)
            if(drangen(dble(bij)).lt.
     .         1d0-Phiexpo(rho,rho,zbrmax,1d0,1d0,sngl(s),bij))then 
              w=Zpair(rho,0.,0.,0.,0.,0.,bk(kk),-1)
              ZpartPro(1)=ZpartPro(1)+w
c              jproj(3,ipp)=1
              ZpartPro(2)=max(ZpartPro(2),bk(kk))
            endif
          endif
        enddo
      endif
      iomega=iomegasave
      call DefXminDf(s)
      end

c-----------------------------------------------------------------------
      function  Zpair(rho1,rho2,zgl,eflux,eflyx,eflix,bi,iqq)
c-----------------------------------------------------------------------
c KW: has nothing in common with the original Zpair.
c
c Now Zpair computes essentially zparpro(1/2,k) and zpartar(1/2,k) in getZpar.
c Use Zp1,Zp2,Zt1,Zt2 for  these.
c  epsG = -Zp2 -Zt2
c  epsGp=-min(...,Zp1)
c  epsGt=-min(...,Zt1)
c     all being negative in case of positive Z
c  atildp   = ... / xminDf**epsG  (in GfunParK)
c  btildpp  = ... +epsGp
c  btildppp = ... +epsGt
c  Womegak = ... + atildp * xp**btildpp * xm**btildppp
c So finally:                               1|*        suppression at large xp,xm
c                  xp^epsGp  *  xm^epsGt     |   *     for positive zpar
c             G ~ -----------------------    |     *   (-> negatice eps
c                     xminDF^epsG            |       *
c                                            +-----------> xp,xm
c So positive Zpair makes screening.
c-----------------------------------------------------------------------
c  returns the z which determine the saturation scales for a given pair
c  via cross-section calculation
c  iqq=-2 for non-diff no sat for Q2s
c  iqq=-1 for non-diff no sat
c  iqq=0 for diff no sat
c  iqq=1 for non-diff with sat and x screening
c  iqq=2 for non-diff with sat and s screening
c  DO NOT CHANGE WITHOUT CHECKING CROSS-SECTIONS AND MULTIPLICITIES AND TP/xsc.f
c-----------------------------------------------------------------------
#include "aaa.h"
#include "par.h"
      parameter(eps=0.01)
      Zpair=0.
      b=abs(bi)
c      sig=sign(1.,bi)
      eflu=eflux
      efly=eflyx
      efli=eflix
      if(maproj.eq.1.and.matarg.eq.1)then
        efli1=eflu
        efli2=efly
        eflu=0 !eflu, efly cannot be used for pp
        efly=0 !would change xsections
      else
        efli1=0.5
        efli2=1.5
      endif
      if(iscreen.gt.0)then
        dummy=zgl
        zrho=rho1+rho2
        drho=0.
        if(zrho.gt.0.)drho=abs(rho1-rho2)
        b2x=b2xscr            !use amplitude width to define bcut (to avoid double bump profiles)
        fegy=fegypp
        fmn=0.
        bx=0.
        bc=0.
        if(iqq.eq.0)then
          b2x=epscrp**2         !max(0.5*b2x,b2x+epscrp*fegy)
          epsmx=0.                !abs(epscrx)
          b2y=exp(-epscrb/fegy)  !/(slopom*log(max(1.,engy*engy))))
        else
          b2y=0.
          epsmx=epscrx
          if(iqq.ge.-1)then
            b2x=2.*b2x          !max(0.5*b2x,b2x+epscrp*fegy)  !behavior of Q2s(b) depends on that
            bx=sqrt(b2x)
            bc=sign(1.,epsmx)
            if(iqq.lt.0)bc=abs(bc)
          else                 !iqq=-2 for Q2s
            b2x=2.*b2x        !bsatur
          endif
          epsmx=abs(epsmx)
          call screen2paramget5(f1,f2,f3,f4,f5)
          epsmx=epsmx+ znurho * zrho 
          if(abs(f1).gt.1e-5)then
            Z=max(0.,ng1evt-2.)/f5
            Z=min(1.,Z)
            epsmx=epsmx * (1+f1*Z**f2)   !f3,f4 not used
          endif
          fegy=epsmx  
          !max(epsmx,epscrw) !+zrhoinc*zrho)
        endif
        if(abs(bc).gt.0.)bc=bc*zbcutx !bcutz(fegy,b2x,epsmx)
        b2a=max(0.,abs(b)-bx)**2/b2x!*4.
        b2c=b**2/b2xscr!/epscrs
        Z0=0.
        if(iqq.ne.0)then
          Zpair=fegy        !*exp(-b2c)!+epscrw*exp(-b2b)
          if(bc.lt.0.)then
            b2b=(abs(b)/abs(bc))**2
            if(iqq.eq.1)then
              Z0=epsmx*epscrg+zrhoinc*drho -efly -eflu !*max(1.,zgl-10)  !*sig
            elseif(iqq.eq.2)then
              Z0=epsmx*epscrs+zrhoinc*drho -efly -eflu !*max(1.,zgl-10)  !*sig
            endif
            Z0=Z0*exp(-b2b/(bc*bc)*4.)
          endif
          if(abs(b).ge.bx)Zpair=Zpair*exp(-b2a)
        else
          b2b=max(0.,abs(b)+0.5*sqrt(b2xscr))**2/b2x*4.
          if(b2y.gt.0.)then
            b2b=min(b2b,b2b**b2y)
          else
            b2b=0.
          endif
          Zpair=fegy*exp(-b2b)
          Zpair=Zpair+zrho
        endif
        Zpair=Zpair+Z0
        Zpair=max(0.,Zpair)
        Zpair=Zpair*(1-efli*exp(-(b/efli1)**efli2))
      endif
      end

c-----------------------------------------------------------------------
      function  bcutz(fegy,b2x,epsmx)
c-----------------------------------------------------------------------
c  returns the impact parameter at which screening saturates
c-----------------------------------------------------------------------
#include "aaa.h"
#include "par.h"
      double precision PhiExpo
      dummy=epsmx
      
      bcutz=0.
c      if(fegy.gt.epsmx)then
c        bcutz=-sqrt(b2x)*log(epsmx/fegy)
c      endif
      if(b2x.gt.0.)then
        s=fegy
        zp=0.
        zt=0.
        bmn=0.
        b=bmn
        Phi=sngl(Phiexpo(zp,zt,1.,1d0,1d0,s,b))
        if(Phi.lt.b2x)then
          bmx=7.
          b=bmx
          Phi=sngl(Phiexpo(zp,zt,1.,1d0,1d0,s,b))
          if(Phi.gt.b2x)then
            do while (Phi.lt.0.9*b2x.or.Phi.gt.1.1*b2x)
              b=bmn+0.5*(bmx-bmn)
              Phi=sngl(Phiexpo(zp,zt,1.,1d0,1d0,s,b))
c              print *,b,phi,bmn,bmx,b2x
              if(Phi.lt.0.9*b2x)then
                bmn=b
              elseif(Phi.gt.1.1*b2x)then
                bmx=b
              else
                Phi=b2x
              endif
            enddo
          endif
          bcutz=b
        endif
c      if(iomega.ne.2)bcutz=bcutz*0.5
      endif
      end

c-----------------------------------------------------------------------
      subroutine iniCore(iii)
c-----------------------------------------------------------------------
#include "aaa.h"
#include "ho.h"
      integer nxicomax,nyicomax,nzicomax
      parameter (nxicomax=161,nyicomax=161,nzicomax=45) !nb of bins (uneven!!!)

      common/cnfrx/nfrx 
      common/cdelcore/delcore,egycore
      logical lprint 
      integer nxico,nyico,nzico
      real xminico,xmaxico,yminico,ymaxico,zminico,zmaxico

      call getIcoDim(nxico,nyico,nzico)
      call getIcoBound(
     .     xminico,xmaxico,
     .     yminico,ymaxico,
     .     zminico,zmaxico
     .     )

      if(nfrx.gt.0)return
      lprint=iorsdf.ne.3.and.iorsdf.ne.6  
     .             .and..not.(iorsdf.eq.5.and.ispherio+ihlle.eq.0)
     .             .or.ish.ge.1
      ia=0
      if(maproj.eq.1.or.matarg.eq.1)then
         rr=2.5              !if this is changed please change in hnb.f and ind.f for the effective flow definition
      elseif(maproj.eq.2.or.matarg.eq.2)then
         rr=5.0
      else
        rpj=1.19*maproj**(1./3.)-1.61*maproj**(-1./3.) + 4.
        rtg=1.19*matarg**(1./3.)-1.61*matarg**(-1./3.) + 4.
        rmx=max(rpj,rtg)
        rmn=min(rpj,rtg)
        !We compute overlap in y direction of two spheres with radii r=rmx and r=rmn
        !with rmx >= rmn for an impact parameter b=bimevt in x direction
        y=0 
        if(bimevt.lt.rpj+rtg-0.5)then
         !We compute the intersection of two circles (x+-b/2)^2+y^2=r^2. Subtracting (to get x) and then 
         !adding the two circle equations gives: y^2=-(b^2-(rmx-rmn)^2)*(b^2-(rmx+rmn)^2)/(4*b^2)
         !with extremum at b^2 = sqrt((rmx-rmn)^2*(rmx+rmn)^2) = (rmx-rmn)*(rmx+rmn) = rmx^2 - rmn^2 
         if(bimevt**2.gt.rmx**2-rmn**2)then  !b^2 > b^2(extremum)
           x=(rmx**2-rmn**2)/(2*bimevt) !From subtracting the two cirxcle equations
           y=sqrt(rmx**2-(x+bimevt/2)**2) !from inverting one of the circle equations
         else !small circle is inside the big one
           y=rmn
         endif
        endif
        rr=max( 5. , y )
        ia=3
      endif
      rr1=rr
      dmx=max(1.19*maproj**(1./3.)-1.61*maproj**(-1./3.)
     .       ,1.19*matarg**(1./3.)-1.61*matarg**(-1./3.))
      dmx=dmx/(engy/2/0.94)
      dmx=nint(dmx*100)/100.
      xmaxico =  rr  !used in jintpo
      Z=max(0.,ng1evt-2.)/400.
      tauzer = tauzer1 + Z * ( tauzer2 - tauzer1 ) !initial tau for hydro   
      !tauzer  = max( tauzer , dmx )  
      fee=yhaha/7.9866929    !reference = 2.76 TeV
      sgzico=min( 1.0 , max( 0.5 , fee ) )
      if(ihyskip.eq.1)dsegce=100000.

      !-----------------------------------------------------------------
      if(iii.eq.1) then
         call setIcoDim(nxico,nyico,nzico)
         call setIcoBound(
     .        xminico,xmaxico,
     .        yminico,ymaxico,
     .        zminico,zmaxico
     .        )
         return 
      end if
      !-----------------------------------------------------------------

      rr2=4*delcore
      rr=(rr1+rr2)/2.
      rr=max(rr,2.5)
      gg=egycore /(16*fee) /2.5**2 *2    /50
      if(rr.lt.2.5001.and.gg.gt.1)then
        rr=2.5+0.5*log(gg)
      endif
      floini  = 0  !0.06  !initial flow parameter
      call getRng1(Z) !Z=max(0.,ng1evt-2.)/400.  
      call  getSystemType(isys,amassMax,amassAsy)
      if(isys.le.2)call getRnpom(0.5,Z) !pp,pA,dA
      Z=min(1.,Z)
      radeft=radeft1+Z*(radeft2-radeft1)
      !~~~~~~~~~~ico binning~~~~~~~~
      if(lprint)
     .write(ifmt,'(a,f6.2,a,$)')'xmaxico changed from',xmaxico,' to '
      xminico = -rr
      xmaxico =  rr
      if(lprint)
     .write(ifmt,'(f6.2)')xmaxico
      yminico = -rr
      ymaxico =  rr
      zminico = -11.25 *fee
      zmaxico =  11.25 *fee
      nxico = 3+2*int(rr/radeft)  !1 + 2 * nint(14+0.55*rr)
      if(nxico.gt.nxicomax)then
        write(ifmt,'(a,i3,a,$)')'nxico=',nxico,' reduced to '
        nxico=nxicomax
        write(ifmt,'(i3)')nxico
        write(ifmt,'(a,f5.2,a,$)')'radeft=',radeft,' increased to '
        radeft=rr/((nxico-3)/2)
        write(ifmt,'(f5.2,a)')radeft,' (to avoid: increase nxicomax)' 
      endif
      nyico = nxico
      nzico = 45

      !~~~~~~~hydro grid~~~~~~~~~~~~~~
      rrr = rr * 0.94
      rrr = nint(rrr*10)/10.
      zzz = zmaxico * 0.94
      zzz = nint(zzz*10)/10.
      xminhy = -rrr
      xmaxhy =  rrr
      yminhy = -rrr
      ymaxhy =  rrr
      zminhy = -zzz
      zmaxhy = zzz
      nxhy = 51
      nyhy = 51
      nzhy = 27
      !~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
      !                        tau step
      !~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
      !hydro 3+1D in Bjorken coordinates faces same problem:
      !If we define
      !  dtau = c_t * dx  (*),  dtau / tau = c_z * deta  (**)
      !then the conditions of stability are:
      !  c_t <= 0.5,  c_z <= 0.5
      !If c_t or c_z stays about those values for a small amount of time,
      !then instability may emerge but not have enough time to develop
      !(and its onset is suppressed with numerical diffusion later).
      !If c_t or c_z is much smaller than 1, then numerical diffusion
      !increases in a given direc_tion - transverse or longitudinal.
      !Optimal values would be (found experimentally):
      !---------------------------
      c_t = 0.25
      c_z = 0.5
      !---------------------------
      !In addition, for better integration of geometrical source terms
      !in hydro equations (which also minimizes energy violation) there
      !should be
      !  dtau << tau
      !In practice this condition is satisfied once c_z = 0.5 in (**)
      !and deta is small enough
      !~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~!!~~~~~~~~~~~~~
      dtauhy = c_t * (xmaxhy-xminhy)/(nxhy-1)   !tau step
      !~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
      !   tau step changes
      !      tauone, tautwo, tauthree: dilute grid and increase tau step 
      !~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
      !   scale factors
      !      ifaahlle : change size of grid in all dimensions
      !      ifathlle : change transverse grid density
      !                   changes as well tau step (/ifathlle)
      !      ifazhlle : change longitudinal grid density
      !~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
      ifaahlle = 1
      ifazhlle = 2
      if(rr.le.3)then  !small systems
        ifathlle = 1
        tauone = 0.8      
        tautwo = 1e30
        tauthree = 1e30
      elseif(rr.le.7)then  !medium systems
        ifathlle = 2
        tauone = 0.4      
        tautwo = 1.0
        tauthree = 1e30
      else             !big systems  
        ifathlle = 4
        tauone = 0.2001
        tautwo = 1.0
        tauthree = 1e30 !tauzer + 3.0
      endif
      !~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
      del_tau = dtauhy/ifathlle !used in hydro calculations
      del_tau_z  = c_z * (zmaxhy-zminhy)/(nzhy-1) * tauzer
      !~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
      if(lprint)then
        write(ifmt,'(a,f7.3)')'tau check: del_tau/del_tau_z ='
     .,del_tau/del_tau_z
        write(ifmt,'(a)')'ico binning and hydro grid:'
        write(ifmt,'(i5,2f7.2,4x,i5,2f7.2)')
     .     nxico, xminico, xmaxico, nxhy, xminhy, xmaxhy
        write(ifmt,'(i5,2f7.2,4x,i5,2f7.2)')
     .     nyico, yminico, ymaxico,nyhy, yminhy, ymaxhy
        write(ifmt,'(i5,2f7.2,4x,i5,2f7.2,4x,f9.4)')
     .     nzico, zminico, zmaxico,nzhy, zminhy, zmaxhy , dtauhy
        write(ifmt,'(a,3i4,f7.3)')'ifaahlle,ifathlle,ifazhlle,epsfin:'
     .   ,ifaahlle,ifathlle,ifazhlle,epsfin
        write(ifmt,'(a,2f6.2)')'fofac fofai:',fofac,fofai
      endif

      call setIcoDim(nxico,nyico,nzico)
      call setIcoBound(
     .     xminico,xmaxico,
     .     yminico,ymaxico,
     .     zminico,zmaxico
     .     )

      end

