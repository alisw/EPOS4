C
C  This file is part of EPOS4
C  Copyright (C) 2022 research institutions and authors (See CREDITS file)
C  This file is distributed under the terms of the GNU General Public License version 3 or later
C  (See COPYING file for the text of the licence)
C

c         ------------------------------------------------
c         since EPOS3086 we have om5 = om1
c                  (in earlier versions   om5 = 0.5 * om1)
c         change in : fnorm ???????????????????
c         ------------------------------------------------
c-----------------------------------------------------------------------
      subroutine psahot(kcol,ncolp,iret)
c-----------------------------------------------------------------------
c psahot - showering (semihard parton-parton interaction)
c-----------------------------------------------------------------------
#include "aaa.h"
#include "sem.h"
#include "ems.h"
#include "tab.h"
#include "par.h"
      double precision plc,s,qt2,pt2Born,sBorn,plprt  !,dda,ddb,ddc,dddd,ddaa,ddcc,ddrat
      common/cems5/plc,s
      real ept1(4),ept2(4)
      double precision ept(4),wpt(2),pLadder(4),wplus,wminus,wp0i,wm0i
     *,wp0,wm0,si,smin,xmin,xp1,wpi,wmi,xp,xm,wp1,wm1,wp2,wm2,z1
     *,wp1i,wm1i,pth,xp2,xg,smax,xmax,xpi,xmi
     *,xmax1,xmin1,zp0,psutz,xpmax0,xpmin0,gb0,tmax0,tmin0,zpm,gb
     *,gbyj,tmax,tmin,t,x1min,x1max,t1min,t1max
     *,t1x,xq1,qq,xpmin,xpmax,psuds,ss,xpsat,rapsat,rapmaxsat
     *,uv1,dv1,uv2,dv2,drangen,xomin,xumin
     *,gbymin,gbymax,rrr,dq2mass,pxh1,pxh2,pyh1,pyh2
     *,xxp1pom,xyp1pom,xxp2pom,xyp2pom,xxm1pom,xym1pom,xxm2pom,xym2pom
     *,xxp1poi,xyp1poi,xxp2poi,xyp2poi,xxm1poi,xym1poi,xxm2poi,xym2poi
     *,wpq(2),wmq(2)
     *,wudt,wwuu,wwdd,wwud,wwdu
     *,zzmax,zzmin,xxmax,xxmin,xx1,xx2,zzmn,zzmx,zzz,uuu
     *,xxmax1,xxmin1,xx1x,xx2x,xpBorn,xmBorn,xzero,xnew
      double precision ee44(0:3,0:3)
      double precision EsatValTil,EsatSeaTil
      double precision EsaturValTil,val,val0
      double precision checkWeightSums(2,-klasmax:klasmax)
      common /ccheckWeightSums/ checkWeightSums
      dimension ep3(4),bx(6),qqs(2),q2mass(2),ncr(2)
     *,qmin(2),iqc(2),ncc(2,2),amqt(4),amqti(4)
      dimension ey(3)
      parameter (mjstr=20000)
      common /psar7/  delx,alam3p,gam3p
      common /psar29/eqj(4,mjstr),iqj(mjstr),ncj(2,mjstr),ioj(mjstr),nj
      common /psar30/ iorj(mjstr),ityj(mjstr),bxj(6,mjstr),q2j(mjstr)
      common /cgauss7/ xgauss7(7),wgauss7(7)
      common /testj/  ajeth(4),ajete(5),ajet0(7)
      parameter (ntim=1000)
      common/cprt/pprt(5,ntim),q2prt(ntim),idaprt(2,ntim),idprt(ntim)
     &,iorprt(ntim),jorprt(ntim),nprtj
      common/emsptl/nppr(npommx,kollmx),npproj(mamx),nptarg(mamx)
      integer icp(2),ict(2),iEOL(2),ihit(4)
      real sTYP
      common/fiveSTYP/sTYP(0:4,mxpsar3)
      common /cnemis/nemis(2)
      integer icp1(2),icp2(2),icm1(2),icm2(2)
      integer jcp(nflav,2),jct(nflav,2),jcpr(nflav,2),jctr(nflav,2)
      common/shat_common/shat
      real shat
      character*3 cside 
      real sTPex(4)
      common/cprtx/nprtjx,pprtx(5,2)/ciptl/iptl
      common/prtdur/idprtx(2)
      common/iortimsh2/iordur(ntim),iodur(ntim)        !ala
      common/cptldur/nptldur(2)                        !ala
      common/clcsh/ntrylcsh
      common/cidpomr/idpomr
      character*12 strg
      character *4 chj(10)
      real vSO(4)
      real xkk(11),ykk(0:3,11)
      real fkk(11)
      logical loback
      character fat*30
      double precision epscutSud2
      common/cchtest/ichtest
      data test/0/ntest/0/njstr/mxstr/
      save test,ntest,njstr
      data kerr1 /0/ kerr2 /0/ kerr3 /0/ kerr4 /0/ kerr5 /0/ 
      data kerr6 /0/ kerr7 /0/ 
      save kerr1,kerr2,kerr3,kerr4,kerr5,kerr6,kerr7
      data kktwo /0/ sstwo /0./ ssone /0./ kkone /0/
      save kktwo,sstwo,ssone,kkone
      data ncount/0/ ncount2/0/ ncount3/0/ ncount4/0/
      save ncount, ncount2, ncount3, ncount4
      data ncount10/0/ ncount11/0/ 
      save ncount10, ncount11
      common /cnquud12/ nquu1,nqud1,nquu2,nqud2
      integer ienvi
      common /cienvi/ ienvi
      integer ioTestFact
      common /ctestfact/ ioTestFact
      integer laddTestFact
      common /claddtestfact/ laddTestFact
      integer modeDeltaFzero
      common /cmodeDeltaFzero/ modeDeltaFzero
      integer iSwitchSides
      common /cSwitchSides/iSwitchSides
      integer noptTestFact
      common /cnopttestfact/ noptTestFact
      integer ihqTestFact
      common /chqtestfact/ ihqTestFact
      parameter (jdimtt=6,kdimtt=19)
      real wartt(jdimtt,0:kdimtt)
      real vartt(jdimtt,0:kdimtt)
      data vartt /120*0/
      save vartt
      double precision iutime(5),tiu3,tiu4,tidi
      logical ltime  
      real eytst(3)
      integer ncount5,lcount5
      data ncount5/0/
      save ncount5
      integer ncount6,lcount6
      data ncount6/0/
      save ncount6
      common/cpsutz/dxpsutz
      parameter (mxrapsata=8)
      real rapsata(mxrapsata)
      parameter (mxsplt=8)
      integer kdir(mxsplt)
      real wimi(mxsplt),wipi(mxsplt)
      call utpri('psahot',ish,ishini,3)

      call maxsize_get(1, mxsplit )
      if(mxsplt.lt.mxsplit)stop'ERROR mxsplt too small'

      !--- in ---

      call defParamsSat(kcol
     .  ,ixsplit,wtsplit,wtsplit2,qqsmall,qqfac,qqsfac,qqmax
     .  ,rapsata,asymmsat)

      call  getSystemType(isys,amassMax,amassAsy)
      iptl=nppr(ncolp,kcol)
      ip=iproj(kcol)
      it=itarg(kcol)
      idpomr=idhpr(ncolp,kcol)
      bpomr=bk(kcol) !bhpr(ncolp,kcol)
      idfpomr=idfpr(ncolp,kcol) !link of pomeron, 1=pro&tar

      if(npr(3,kcol).le.1)then
        ptiprj=0. ! if Npom<=1
        ptitrg=0. ! if Npom<=1
        ptisat=0. ! if Npom<=1
        ptiposii=0.
      else
        ptiprj=ptipom + ptipomi * sqrt(npr(3,kcol)-1.) !randow walk !ckw21 
        ptitrg=ptipom + ptipomi * sqrt(npr(3,kcol)-1.) 
        ptiposii=       ptiposi * sqrt(npr(3,kcol)-1.)  !log(max(q2cmin(1),q2cmin(2))/q2sft)
        ptisat=ptipos + ptiposii
      endif
 
      idsprj=idsppr(ncolp,kcol)
      idstrg=idstpr(ncolp,kcol)

      if(ish.ge.3)write(ifch,*)'Start psahot:',ip,it
     *,ncolp,kcol,iptl,idpomr,idfpomr,bpomr

      !--- restart ---

      igo=0 
      nptl0=nptl
      nptl1=nptl
      ntryhot=0
      ntxt80save=ntxt80
  1   ntryhot=ntryhot+1
      ntxt80=ntxt80save
      !if(ntryhot.gt.1)print*,'Redo ntryhot igo:',ntryhot,igo

      !--- in / out ---

      wp0i=sqrt(xpr(ncolp,kcol))*exp(ypr(ncolp,kcol))*plc      
      wm0i=sqrt(xpr(ncolp,kcol))*exp(-ypr(ncolp,kcol))*plc     

      xxp1pomi=xxp1pr(ncolp,kcol) !SE pt
      xyp1pomi=xyp1pr(ncolp,kcol)
      xxp2pomi=xxp2pr(ncolp,kcol)
      xyp2pomi=xyp2pr(ncolp,kcol)
      xxm1pomi=xxm1pr(ncolp,kcol)
      xym1pomi=xym1pr(ncolp,kcol)
      xxm2pomi=xxm2pr(ncolp,kcol)
      xym2pomi=xym2pr(ncolp,kcol)

      idp1pom=idp1pr(ncolp,kcol)   !String end type, from ProSeTy(k,n)
      idp2pom=idp2pr(ncolp,kcol)   ! will be modified, but only for val
      idm1pom=idm1pr(ncolp,kcol)   
      idm2pom=idm2pr(ncolp,kcol)

      do jsplit=1,mxsplit
      do ii=1,2
      call ptprboo_set(ii,jsplit,ncolp,kcol, 0. )
      call rapprboo_set(ii,jsplit,ncolp,kcol, 0. )
      enddo  
      call xpprbor_set(jsplit,ncolp,kcol, 0. )
      call xmprbor_set(jsplit,ncolp,kcol, 0. )
      call idbor_set(1,jsplit,ncolp,kcol, 0 )
      call idbor_set(2,jsplit,ncolp,kcol, 0 )
      call gbpom_set(1,jsplit,ncolp,kcol, 0. )
      call gbpom_set(2,jsplit,ncolp,kcol, 0. )
      enddo
      nemispr(1,ncolp,kcol)=0 
      nemispr(2,ncolp,kcol)=0 
      q2bor(ncolp,kcol)=0 
      shatpr(ncolp,kcol)=0 
      jsplit=1

      !--- init ---

      xxp1sum=0
      xyp1sum=0
      xxp2sum=0
      xyp2sum=0
      xxm1sum=0
      xym1sum=0
      xxm2sum=0
      xym2sum=0

      !--- exit ---

      iret=1
      nptl=nptl1
      if(ntryhot.gt.100)print*,'ntryhot too big ==> exit',idpomr
      if(ntryhot.gt.100)then
        nptl=nptl0
        if(ish.ge.1)
     .  write(ifmt,'(a)')'WARNING pasahot: ntryhot too big -> exit'
        goto 16
      endif

      kxsplit=1 ! Finally best dn_c/dy for kxsplit=1

      !--- REMOVED split into ixsplit sub poms (C) and kxsplit sub poms (A,B)

      !- - - - - - - - - - - - - - - - - - - - - - - - - - - - - -   
      !SPLIT-METHOD T: Only for tests !!!
      ! - - - - - - - - - - - - - - - - - - - - - - - - - - - - - 
      if(idpomr.lt.0)then 
        fxsplit=1./ixsplit   !method C 
      else
        fxsplit=1.
      endif
      if(0.25*wp0i/kxsplit*wm0i/kxsplit.le.100.0)kxsplit=1 !method A/B

      !--- split iteration ---

      jsplit=0
      sump=0
      summ=0 
 666  continue
      jsplit=jsplit+1
      if(jsplit.gt.mxsplit)stop'ERROR Increase mxsplit'

      ifav=-1
      ntxt80save=ntxt80
      nptl2=nptl
      ntrysplit=0
  2   ntrysplit=ntrysplit+1
      !if(ntrysplit.gt.1)print*,'Redo 2',igo,ntrysplit

      if(ntrysplit.gt.1)then
        ntxt80=ntxt80save
        nptl=nptl2
        do ii=1,2
        call ptprboo_set(ii,jsplit,ncolp,kcol, 0. )
        call rapprboo_set(ii,jsplit,ncolp,kcol, 0. )
        enddo  
        call xpprbor_set(jsplit,ncolp,kcol, 0. )
        call xmprbor_set(jsplit,ncolp,kcol, 0. )
        call idbor_set(1,jsplit,ncolp,kcol, 0 )
        call idbor_set(2,jsplit,ncolp,kcol, 0 )
        call gbpom_set(1,jsplit,ncolp,kcol, 0. )
        call gbpom_set(2,jsplit,ncolp,kcol, 0. )
      endif

      if(jsplit.eq.1.and.ntrysplit.eq.1)kdir(1)=0
      if(idpomr.lt.0.and.jsplit.eq.1.and.ntrysplit.eq.1)then ! sat pom
        wimi(1)=wm0i
        wipi(1)=wp0i
      endif 

      if(idpomr.lt.0)then ! sat pom
        wp0=wipi(jsplit)
        wm0=wimi(jsplit)
      else
        wp0=wp0i
        wm0=wm0i
      endif
     
      if(idpomr.lt.0)then ! sat pom  
        fsplit=1./kxsplit !simple sharing for SE pt
      else
        fsplit=1. 
      endif

      xxp1pom=xxp1pomi*fsplit
      xyp1pom=xyp1pomi*fsplit
      xxp2pom=xxp2pomi*fsplit
      xyp2pom=xyp2pomi*fsplit
      xxm1pom=xxm1pomi*fsplit
      xym1pom=xym1pomi*fsplit
      xxm2pom=xxm2pomi*fsplit
      xym2pom=xym2pomi*fsplit

      !--- initializations ---

      fat='                              '
      gbymin=1d30
      gbymax=-1d30
      nrejtot=0
      ktry41=0
      ktry10=0
      ktry11=0
      delta=0.01 !bg
      nj=0
      q2finsave=q2fin

      call timer(iutime)
      tiu3=iutime(3)
      tiu4=iutime(4)
      ltime=.false.

      !kkkkkkkkkkk!KW1909
      if(ioTestFact.ne.0)then
        q2cmin(1)=q2sft !q2nmin !actually already 
        q2cmin(2)=q2sft !q2nmin !defined in ems
        if(ioTestFact.ge.3)then
          idpomr=0 !sea-sea
        elseif(ioTestFact.le.-3)then
          idpomr=1 !val-sea
          idfpomr=1
        elseif(ioTestFact.eq.-2)then
          idpomr=2 !sea-val
          idfpomr=1
        elseif(ioTestFact.eq.-1)then
          idpomr=3 !val-val
          idfpomr=1
        endif
      endif
      !kkkkkkkkkkk

      if(idfpomr.eq.0)stop'idfpomr=0'
      if(ish.ge.3)then
        write(ifch,'(1x,i4,3x,4(e11.5,1x),i2,1x,i3)')iptl
     *  ,sqrt(pptl(1,iptl)**2+pptl(2,iptl)**2),pptl(3,iptl)
     *  ,pptl(4,iptl),pptl(5,iptl)
     *  ,istptl(iptl),ityptl(iptl)
      endif
      istptl(iptl)=31

      do i=1,2
        icp(i)=icproj(i,ip)
        ict(i)=ictarg(i,it)
      enddo

      !--- reinitialize all intermediate stack values ---

      call initPartonList(1,njstr)
      q2fin=q2finsave
      iret=0
      nj=0
      if(iremn.ge.2)then
        call iddeco(icp,jcp)
        call iddeco(ict,jct)
      endif
      factbqq=1.
      factbgg=1.

      do njx=nj+1,nj+10
        if(ncj(1,njx).ne.0.or.ncj(2,njx).ne.0)then
          write(ifmt,'(a)')
     .    'ERROR psahot (1) ncj non-zero, wrong initialization'
        endif
      enddo

      !kkkkkkkkkkk!KW1908
      if(ioTestFact.ne.0)then
        wp0=plc
        wm0=plc
      endif
      !kkkkkkkkkkk

      wpi=0d0
      wmi=0d0
      wp1=0d0
      wm1=0d0
      wp2=0d0
      wm2=0d0
      wpq(1)=0d0
      wmq(1)=0d0
      wpq(2)=0d0
      wmq(2)=0d0
      
      xp=wp0/plc
      xm=wm0/plc
      ss=wp0*wm0                !lc+*lc- for the semihard interaction
      sPomeron=ss

      !--- initial pt from string ends ---

      amqt(1)=sqrt(sngl(xxp1pom**2+xyp1pom**2))
      amqt(2)=sqrt(sngl(xxp2pom**2+xyp2pom**2))
      amqt(3)=sqrt(sngl(xxm2pom**2+xym2pom**2))
      amqt(4)=sqrt(sngl(xxm1pom**2+xym1pom**2))
      amqpt=amqt(1)+amqt(2)+amqt(3)+amqt(4)

      qqcut=max(q2cmin(1),q2cmin(2)) !virtuality cutoff
      !~~~~~~
      !An extra variable qqcut is not needed but we keep it, it could take the 
      !role of some imposed high pt cutoff. Used presently for sat poms !!!!!!! change it <=======
      !~~~~~~ 

      !-----------------------------------------------------------
      ! Pseusosoft Pom PDF (PPPDF): ifav, rapsatrel
      !-----------------------------------------------------------

      if(idpomr.lt.0)then
        probot=0.020!0.035!0.05!0.1!0.02
        if(ifav.lt.0)then 
         if(rangen().gt.probot)then !...charm pseudo pom
          ifav=4
         else!..........................bottom pseudo pom
          ifav=5
         endif
        endif
        if(ifav.ge.4)then
          c=2*rangen()-1
          rapsatrel=c
          ptisatrap=0
          if(ptisat.gt.0.)ptisatrap= min( 1.0 , ptipos) + ptiposii !avoid big values for HF
        else
          a=0 !0.20       !afrapsat
          b=3*(0.5-a)     !frapsat(x)=a+b*x^2  \int = 1
          c=2*rangen()-1  !use superposition method w*f1+(1-w)*f2, w=2a
          if(rangen().gt.2*a) c = sign(1.,c) * (abs(c))**(1./3)
          rapsatrel=c
          ptisatrap= ptisat  * (2.1-1.9*abs(rapsatrel))
        endif
      endif

      !------------------------------------------------------
      ! Pseusosoft Pom PDF (PPPDF): s2min, s2max  and  qqsat 
      !------------------------------------------------------

      facq2tim=1
      if(idpomr.le.-1)then ! sat pom    
        qqsma=qqsmall
        probot=0.020!0.035!0.05!0.1!0.02
        if(ifav.eq.5)
     .   qqsma=qqsmall+(qbmass-qcmass)*2.5!2!3!4! Threshold effect: 2-2.5 jump, then no change
        qqupp=min(0.25*sngl(ss), max(q2cmin(1),q2cmin(2)) )
        qqlow=0.5 !qqupp!q2sft
        qqff=qqupp*qqfac
        hh=qqsma**2/qqff
        if(rangen().lt.1./(1.+hh))then
          qqsat=-alog(rangen())+hh
        else
          qqsat=hh !*rangen()
        endif
        qqsat=qqsat*qqff 
        coelafSave=coelaf
        if(ifav.eq.4)then !........charm pseudo pom
          coelaf=9
          call HardScale(4,qqsat,pt2xx, 1 ) 
          if(pt2xx.lt.0.001)pt2xx=-log(rangen())-log(rangen()) !=(2-2*sqrt(rangen())) !law ~ \ = 1-0.5x{0-2}
          qqssfac=2.0!4.0
          qqmass=qcmass**2
          facq2tim=1
          ptemi=1.0*(rangen()+rangen()) !law ~ /\ = x{0-1};2-x{1-2}
        elseif(ifav.eq.5)then!.....bottom pseudo pom
          coelaf=9
          call HardScale(5,qqsat,pt2xx, 1 ) 
          if(pt2xx.lt.0.001)pt2xx=10*(2-2*sqrt(rangen())) !law ~ \ = 1-0.5x{0-2}
          qqssfac=5!10.0
          qqmass=qbmass**2
          facq2tim=5
          ptemi=3.5*(rangen()+rangen()) !law ~ /\ = x{0-1};2-x{1-2}
        else
          coelaf=coelafSave
          call HardScale(1,qqsat,pt2xx, 1 ) 
          if(pt2xx.lt.0.001)stop'ERROR 230815'
          qqssfac=2.0!4.0
          qqmass=0.
          facq2tim=1
          ptemi=1.0*(rangen()+rangen()) !law ~ /\ = x{0-1};2-x{1-2}
        endif
        if(ifav.ge.4.and.qqsma**2.ge.coelaf*qqmass)then !otherwise spike in pt distr
          print*,'qqsma,coelaf,qqmass:',qqsma,coelaf,qqmass
          stop'ERROR 23052023'
        endif
        coelaf=coelafSave
        qqs(1)=qqsat
        qqs(2)=qqsat
        qqcut=qqsat
        qqxxx=max(qqsat,pt2xx)
        s2min=4.*(qqxxx+qqmass) * qqssfac  
        s2min1_save=s2min
        s2max=4.*(qqxxx+qqmass) * qqssfac 
        !===> maybe update needed for "frej" and "fscale1sat"  
        s2min=min(s2min,ss)
        s2max=min(s2max,ss)
        s2min=(s2min+rangen()*(s2max-s2min))*0.999
        s2max=s2min/0.999*1.001
      else
        qqs(1)=q2cmin(1)
        qqs(2)=q2cmin(2)
        s2min=4.*qqcut  
        s2max=ss
      endif
      s2min2_save=s2min

      if(sngl(ss).le.(sqrt(s2min)+amqpt)**2)then 
        if(kxsplit.gt.1)then
          if(jsplit.lt.kxsplit)goto 666
        endif
        if(ish.ge.1)
     .  write(ifmt,'(a)')
     .  'WARNING psahot: Insufficient pomeron mass -> exit'
        iret=1
        imark=1
        goto 16
      endif

      !----------------------------------------------------
      !              generate intrinsic pt
      !----------------------------------------------------

      if(idpomr.ne.3
     .   .and..not.(ioTestFact.ne.0.or.noptTestFact.eq.1))then
c       alpha=-1.    !generation of pt as pt**(alpha-1)
c       a little bit of pt is necessary in case of sea quarks
        ptmx=0.2*(sngl(sqrt(ss))-sqrt(s2min)-amqpt)
        if(idpomr.ge.0)then
          pt=min(ptmx,ptiprj)
          pt=max(0.001,pt)
          pt=pt*ranptjcut(ptmx/pt)
        else
          pt=min(ptisatrap,ptmx)
          pt=max(0.001,pt)
          pt=pt*ranptjcut(ptmx/pt)
        endif
        if(idpomr.le.0.or.idpomr.eq.2.and.pt.gt.0.)then
          phi=2.*pi*rangen()
          pxh1=dble(pt*cos(phi))
          pyh1=dble(pt*sin(phi))
        else
          pxh1=0d0
          pyh1=0d0
        endif
        if(idpomr.ge.0)then
          pt=min(ptmx,ptitrg)
          pt=max(0.001,pt)
          pt=pt*ranptjcut(ptmx/pt)
        else
          pt=min(ptisatrap,ptmx)
          pt=max(0.001,pt)
          pt=pt*ranptjcut(ptmx/pt)
        endif
        if(idpomr.le.0.or.idpomr.eq.1.and.pt.gt.0.)then
          phi=2.*pi*rangen()
          pxh2=dble(pt*cos(phi))
          pyh2=dble(pt*sin(phi))
        else
          pxh2=0d0
          pyh2=0d0
        endif
      else
        pxh1=0d0
        pyh1=0d0
        pxh2=0d0
        pyh2=0d0
      endif
      amqpt=amqpt+sqrt((pxh1+pxh2)**2+(pyh1+pyh2)**2)   !intrinsic pt compensated in soft part
      if(idpomr.le.0)then        !sea-sea and sat
        continue
      elseif(idpomr.eq.1)then    !val-sea 
        pxh1=pxh1+xxp1pom
        pyh1=pyh1+xyp1pom
        amqpt=amqpt-amqt(1)
      elseif(idpomr.eq.2)then    !sea-val 
        pxh2=pxh2+xxm2pom
        pyh2=pyh2+xym2pom
        amqpt=amqpt-amqt(3)
      elseif(idpomr.eq.3)then    !val-val 
        pxh1=xxp1pom
        pyh1=xyp1pom
        pxh2=xxm2pom
        pyh2=xym2pom
        amqpt=amqpt-amqt(1)-amqt(3)
      else
        stop'unknown pomeron'
      endif
      !kkkkkkkkkkk!KW1909
      if(ioTestFact.ne.0.or.noptTestFact.eq.1)then
        pxh1=0d0
        pyh1=0d0
        pxh2=0d0
        pyh2=0d0
      endif
      !kkkkkkkkkkk
      pth=(pxh1+pxh2)**2+(pyh1+pyh2)**2
c      print *,idpomr,pth

      !-----------------------------------------------------------
      ! smin,smax,xmin,xmax
      !-----------------------------------------------------------

      if(idpomr.ge.0)then
        smin = ( sqrt(dble(s2min))+sqrt(pth) )**2    !hard pomeron s_min
        smax = ( sqrt(dble(s2max))-dble(amqpt) )**2  !hard pomeron s_max
      else
        smin = s2min+pth   
        smax = s2max+pth   
      endif
      pth1_save=pth
      s2minPlusPth_save=smin
      xmin = smin/ss
      xmax = smax/ss
      if(idpomr.ge.0)xmin = min(xmin,(sqrt(s2max)-ammsqd)**2/ss)
      call setmirror(0)

      !-----------------------------------------------------------
      ! Pseusosoft Pom PDF (PPPDF): xpsat
      !-----------------------------------------------------------

      if(idpomr.lt.0)then
        mirror=0
        !if(maproj.eq.matarg.and.rangen().lt.0.5)mirror=1
        !call setmirror(mirror)
        rappomi=0.5*log(wp0i/wm0i)  !rap of initial Pom in Lab frame
        rappom=0.5*log(wp0/wm0)   !rap of split Pom in Lab frame
        zpmsat=(xmin+xmax)/2.
        xpmaxsat=psutz(ss,zpmsat*ss,dble(amqpt**2)) 
        rapmaxsat=0.5*log(xpmaxsat**2/zpmsat)
        !rapsat=-rapmaxsat+rangen()*2*rapmaxsat !fine tuning sat 
        if(maproj.gt.matarg)then
          rap1=max(0.,rapmaxsat-rapsata(1))
          rap2=max(0.,rapmaxsat-rapsata(2))
        else
          rap1=max(0.,rapmaxsat-rapsata(2))
          rap2=max(0.,rapmaxsat-rapsata(1))
        endif
        if(mirror.eq.1)call swap(rap1,rap2)
        n555=0
  555   continue
        n555=n555+1
        !rapsat=-rap1
        !if(rangen().lt.0.5)rapsat=rap2
        rapsat=(-rap1+rap2)/2 + rapsatrel * (rap2+rap1)/2
        if(abs(rapsat).gt.rapsata(3).and.n555.lt.20)goto 555
        wtrapsat=1. ! Finally best dn_c/dy for wtrapsat=1.
        !- - - - - - - - - - - - - - - - - 
        !Redo if ifav and finally created flavor do not agree 
        !  see 'ifav redo' below
        !- - - - - - - - - - - - - - - - - 
        xpsat=sqrt(zpmsat)*exp(rapsat)         
      endif
   
      !----------------------------------------------------
      ! definitions, checks
      !----------------------------------------------------

      if(xmax.le.xmin.or.xmax.gt.1d0)then
        if(noptTestFact.eq.0.and.idpomr.ne.-1)then
          igo=1
          goto 1                !to redefine pth if too large
        else
          iret=1
          imark=2
          goto 16                !no pth, reject
        endif
      endif
      if(iremn.ge.2.and.ioTestFact.eq.0)then
        if(iclpro.eq.2)then
          if(iabs(idptl(ip)).eq.1120)then !proj=proton
            nquu1=jcp(1,1)+jcp(1,2)
            nqud1=jcp(2,1)+jcp(2,2)
          elseif(iabs(idptl(ip)).eq.1220)then !proj=neutron
            nquu1=jcp(2,1)+jcp(2,2)
            nqud1=jcp(1,1)+jcp(1,2)
          else    !to avoid flavor problem with exotic projectile (but should not happen (only gg)
            nquu1=0
            nqud1=0
          endif
        elseif(iclpro.eq.1)then
          if(iabs(idptl(ip)).eq.120)then
            nquu1=jcp(1,1)+jcp(1,2)
            nqud1=jcp(2,1)+jcp(2,2)
          else    !to avoid flavor problem with exotic projectile (but should not happen (only gg)
            nquu1=0
            nqud1=0
          endif
        elseif(iclpro.eq.3)then
          if(iabs(idptl(ip)).eq.130)then !proj=Kch
            nquu1=jcp(1,1)+jcp(1,2)
            nqud1=jcp(3,1)+jcp(3,2)
          elseif(iabs(idptl(ip)).eq.230)then  !proj=K0
            nquu1=jcp(2,1)+jcp(2,2)
            nqud1=jcp(3,1)+jcp(3,2)
          else    !to avoid flavor problem with exotic projectile (but should not happen (only gg)
            nquu1=0
            nqud1=0
          endif
        else                    !charm
          if(iabs(idptl(ip)).eq.140)then
            nquu1=jcp(1,1)+jcp(1,2)
            nqud1=jcp(4,1)+jcp(4,2)
          elseif(iabs(idptl(ip)).eq.240)then
            nquu1=jcp(2,1)+jcp(2,2)
            nqud1=jcp(4,1)+jcp(4,2)
          elseif(iabs(idptl(ip)).eq.340)then
            nquu1=jcp(3,1)+jcp(3,2)
            nqud1=jcp(4,1)+jcp(4,2)
          else
            nquu1=jcp(4,1)
            nqud1=jcp(4,2)
          endif
        endif
        if(iabs(idptl(maproj+it)).eq.1220)then !targ=neutron
          nquu2=jct(2,1)+jct(2,2)
          nqud2=jct(1,1)+jct(1,2)
        else
          nquu2=jct(1,1)+jct(1,2)
          nqud2=jct(2,1)+jct(2,2)
        endif
      else
        nquu1=2
        nqud1=1
        nquu2=2
        nqud2=1
      endif

      if(idpomr.le.0)then        !sea-sea 
        iqq=max(-1,idpomr)   ! -1 - sat, 0 - sea-sea, 1 - val-sea, 2 - sea-val, 3 - val-val
      elseif(idpomr.eq.1)then    !val-sea 
        iqq=1
        if(nquu1+nqud1.le.0)then
          write(ifmt,'(a)')'WARNING No more valence quark (1)'
          iret=1
          imark=3
          goto 16
        endif
      elseif(idpomr.eq.2)then    !sea-val 
        iqq=2
        if(nquu2+nqud2.le.0)then
          write(ifmt,'(a)')'WARNING No more valence quark (2)'
          iret=1
          imark=4
          goto 16
        endif
      elseif(idpomr.eq.3)then    !val-val 
        iqq=3
        if(nquu1+nqud1+nquu2+nqud2.le.0)then
          if(ish.ge.1)
     .    write(ifmt,'(a,6i5)')
     .    'WARNING psahot: No more val quark (3) -> exit '
     .         ,ip,it,nquu1,nqud1,nquu2,nqud2
          iret=1
          imark=5
          goto 16
        endif
      else
        stop'unknown pomeron'
      endif

      ipomtype=iqq

      nj0=nj
      ih=iproj(kcol)
      jh=itarg(kcol)
      do l=1,4
        bx(l)=xorptl(l,iptl)
      enddo
      bx(5)=tivptl(1,iptl)
      bx(6)=tivptl(2,iptl)
      ity=ityptl(iptl)

      xomin=0.d0
      xumin=0.d0
      coef_gb0=30.
      if(iscreen.eq.0)coef_gb0=50.

      if(idpomr.ge.0)then
        xpmax0=psutz(smax,smin,dble(amqpt**2)) !max x+/- for the hard P
      else
        xpmax0=psutz(ss,smin,dble(amqpt**2)) !max x+/- for the hard P
      endif 
      if(xpmax0.le.0)write(*,*)'negative xpmax0:',idpomr,xpmax0
     .,sngl(smax),sngl(smin),amqpt**2
      !if(idpomr.lt.0)print*,'TESTpsahot',xpmax0,dxpsutz,'    '
      !.,sngl(smax),sngl(smin),amqpt**2
      xpmin0=min(1d0,xmin/xpmax0)              !min x+ for the hard pomeron
      xpmin0=max(xpmin0,xomin)
      if(xumin.gt.0.)xpmax0=min(xpmax0,xmax/xumin)
      if(iqq.eq.1.or.iqq.eq.2)xpmax0=min(xpmax0,.99998d0) !to have psdfh4>0
      rp=r2had(iclpro)+r2had(icltar)
     .-slopom*log(xmin) !*log(1d0/s)
      z=exp(-bpomr**2/(4.*.0389*rp)) !coef for rejection
      if(z.eq.0)then
       write(ifch,*)'psahot : z,ih,jh ! -> ',z,ih,jh
       call gakli2(ih,ih)
       call gakli2(jh,jh)
       call gakli2(iptl,iptl)
       stop
      endif

      if(smax.le.smin)then
        if(ish.ge.1)
     .  write (ifmt,*)'WARNING psahot: smax < smin -> exit ',smax,smin
        iret=1
        imark=6
        goto 16
      endif

      xp1=wp0/plc              !lc+ share for the semihard interaction
      xp2=wm0/plc              !lc- share for the semihard interaction

      do njx=nj+1,nj+10
        if(ncj(1,njx).ne.0.or.ncj(2,njx).ne.0)then
          write(ifmt,'(a)')
     .    'ERROR psahot (2) ncj non-zero, wrong initialization'
        endif
      enddo



c###########################################################################################
c###########################################################################################
c        determine LC momenta wpi,wmi for hard Pomeron
c###########################################################################################
c###########################################################################################


      if(ish.ge.4)write(ifch,*)
     & 'determine LC momenta wpi,wmi for hard Pomeron'


      iq1=0
      iq2=0
      wwdd=0.d0
      do i=0,3
      do j=0,3
        ee44(i,j)=0
      enddo
      enddo

      if(ioTestFact.eq.0.and.iqq.ge.0)then
        call tabuSTYPZ(iqq,smin-pth,smax,qqs,qqcut) !tabulates sTYPz -> /ctabSTYPZ/
      endif

      !-----------------------------------------------------------
      if(iqq.eq.3)then          !     val-val  iqq=3
      !-----------------------------------------------------------

        if(ish.ge.4)write(ifch,*)'val-val'
        xmin1=xmin**dble(delh+.4)
        xmax1=xmax**dble(delh+.4)
        zp0=dsqrt(xmin)
        if(zp0.ge.1.d0)call utstop('zp0 in sem&')
        !........ kinematical bounds
        tmax0=dlog((1.d0+dsqrt(1.d0-zp0))/(1.d0-dsqrt(1.d0-zp0)))
        tmin0=dlog((1.d0+dsqrt(1.d0-xpmax0))
     *                   /(1.d0-dsqrt(1.d0-xpmax0)))

        !-----------------------
        !rejection normalization
        !-----------------------

        sss=sngl(smax-pth)
        if(ioTestFact.eq.0)then
          call pipolSTYPZ(sss)                       !interpolates sTYPz
        else
          call setQvalues(qqs(1),qqs(2),qqcut) 
          call calcSTYP(iqq,sss)                 !computes sTYPz -> /zeroSTYP/
        endif
        call calcWW4(iqq,ee44,sgg,sgq,sqg,sqq)  !ee44 here not used!!
        sqqp = sqg !special meaning (sqg is zero)
        sqaq = sgq !special meaning (sgq is zero)
        if(iclpro.eq.4)stop'add if...endif in calcWW4'

        if(iremn.ge.2)then
          if(nquu1.gt.nqud1.or.iclpro.ne.2)then
            uv1=zp0**(-dels(1))*EsaturValTil(1,zp0,xp1,qqs(1),1,iclpro)
            dv1=zp0**(-dels(1))*EsaturValTil(1,zp0,xp1,qqs(1),2,iclpro)
          else                  !if nquu<nqud => no u or no d
            uv1=zp0**(-dels(1))*EsaturValTil(1,zp0,xp1,qqs(1),2,iclpro)
            dv1=uv1
          endif
          if(nquu1.eq.0)uv1=0d0
          if(nqud1.eq.0)dv1=0d0
          if(nquu2.gt.nqud2)then
            uv2=zp0**(-dels(2))*EsaturValTil(2,zp0,xp2,qqs(2),1,icltar)
            dv2=zp0**(-dels(2))*EsaturValTil(2,zp0,xp2,qqs(2),2,icltar)
          else                  !if nquu<nqud => no u or no d
            uv2=zp0**(-dels(2))*EsaturValTil(2,zp0,xp2,qqs(2),2,icltar)
            dv2=uv2
          endif
          if(nquu2.eq.0)uv2=0d0
          if(nqud2.eq.0)dv2=0d0
        else
          uv1=zp0**(-dels(1))*EsaturValTil(1,zp0,xp1,qqs(1),1,iclpro)
          dv1=zp0**(-dels(1))*EsaturValTil(1,zp0,xp1,qqs(1),2,iclpro)
          uv2=zp0**(-dels(2))*EsaturValTil(2,zp0,xp2,qqs(2),1,icltar)
          dv2=zp0**(-dels(2))*EsaturValTil(2,zp0,xp2,qqs(2),2,icltar)
        endif
        wwuu=uv1*uv2*sqq
        if(iclpro.eq.2)then
          wwdd=dv1*dv2*sqq
        elseif(iclpro.eq.1)then
          wwdd=dv1*dv2*sqaq
        elseif(iclpro.eq.3)then
          wwdd=dv1*dv2*sqqp
        elseif(iclpro.eq.4)then
          wwuu=uv1*uv2*sqqp
          wwdd=0.
        endif
        wwud=uv1*dv2*sqqp
        wwdu=dv1*uv2*sqqp
        wudt=wwuu+wwdd+wwud+wwdu
        frej=1.     !coef for rejection
        gb0=dble(wudt)/xmax**dble(delh)/xmin**0.4
     *   *(tmax0-tmin0)*z*coef_gb0/frej
     *   *(1.d0-zp0)**dble(.5)
           ! /(4.*pi*(r2had(iclpro)+r2had(icltar)))       !cancel in gb

        !-------------------------
        ! xp,xm proposal/rejection 
        !-------------------------
   
        ntryrej=0
 3      zpm=(xmin1+drangen(xmin1)*(xmax1-xmin1))**dble(1./(delh+.4)) !zpm proposition
        ntryrej=ntryrej+1
        if(ntryrej.gt.999999)then
          write(ifmt,'(a,a)')
     .    'WARNING psahot too many rejections for val-val',' -> exit'
          imark=7
          goto 16
        endif

        sss=sngl(zpm*ss-pth)
        if(ioTestFact.eq.0)then
          call pipolSTYPZ(sss)                 !interpolates sTYPz -> /zeroSTYP/ 
        else
          call setQvalues(qqs(1),qqs(2),qqcut) 
          call calcSTYP(iqq,sss)                    !computes sTYPz -> /zeroSTYP/
        endif
        call calcWW4(iqq,ee44,sgg,sgq,sqg,sqq) !ee44 here not used!!
        sqqp = sqg !special meaning (sqg is zero)
        sqaq = sgq !special meaning (sgq is zero)
        if(iclpro.eq.4)stop'add if...endif in calcWW4'

        xpmax=psutz(ss,zpm*ss,dble(amqpt**2))  !max x+ for sh pomeron
        if(xpmax.lt.0)write(*,*)'xpmax',xpmax
        tmax=dlog((1.d0+dsqrt(1.d0-dsqrt(zpm)))
     *             /(1.d0-dsqrt(1.d0-dsqrt(zpm))))
        tmin=dlog((1.d0+dsqrt(1.d0-xpmax))/(1.d0-dsqrt(1.d0-xpmax)))
        t=(tmin+drangen(tmin)*(tmax-tmin))
        xp=1.d0-((1.d0-dexp(-t))/(1.d0+dexp(-t)))**2  !x+_v proposition
        xm=zpm/xp                             !x-_v
        if(xm.gt.xp.and.ish.ge.1)write(ifmt,*)'xm,xp',xm,xp
        gb=(tmax-tmin)
     *     *(1.d0-xp)**(.5)
        if(rangen().lt..5)then
          xp=xm
          xm=zpm/xp
        endif

        if(iremn.ge.2)then
          if(nquu1.gt.nqud1.or.iclpro.ne.2)then
            uv1=xp**(-dels(1))*EsaturValTil(1,xp,xp1,qqs(1),1,iclpro)
            dv1=xp**(-dels(1))*EsaturValTil(1,xp,xp1,qqs(1),2,iclpro)
          else                  !if nquu<nqud => no u or no d
            uv1=xp**(-dels(1))*EsaturValTil(1,xp,xp1,qqs(1),2,iclpro)
            dv1=uv1
          endif
          if(nquu1.eq.0)uv1=0d0
          if(nqud1.eq.0)dv1=0d0
          if(nquu2.gt.nqud2)then
            uv2=xm**(-dels(2))*EsaturValTil(2,xm,xp2,qqs(2),1,icltar)
            dv2=xm**(-dels(2))*EsaturValTil(2,xm,xp2,qqs(2),2,icltar)
          else                  !if nquu<nqud => no u or no d
            uv2=xm**(-dels(2))*EsaturValTil(2,xm,xp2,qqs(2),2,icltar)
            dv2=uv2
          endif
          if(nquu2.eq.0)uv2=0d0
          if(nqud2.eq.0)dv2=0d0
        else
          uv1=xp**(-dels(1))*EsaturValTil(1,xp,xp1,qqs(1),1,iclpro)
          dv1=xp**(-dels(1))*EsaturValTil(1,xp,xp1,qqs(1),2,iclpro)
          uv2=xm**(-dels(2))*EsaturValTil(2,xm,xp2,qqs(2),1,icltar)
          dv2=xm**(-dels(2))*EsaturValTil(2,xm,xp2,qqs(2),2,icltar)
        endif
        wwuu=uv1*uv2*sqq
        if(iclpro.eq.2)then
          wwdd=dv1*dv2*sqq
        elseif(iclpro.eq.1)then
          wwdd=dv1*dv2*sqaq
        elseif(iclpro.eq.3)then
          wwdd=dv1*dv2*sqqp
        elseif(iclpro.eq.4)then
          wwuu=uv1*uv2*sqqp
          wwdd=0.
        endif
        wwud=uv1*dv2*sqqp
        wwdu=dv1*uv2*sqqp
        if(idfpomr.eq.1)then
          rh=r2had(iclpro)+r2had(icltar)!+slopom*real(dlog(1./zpm))
          z1=exp(-bpomr**2/(4.*0.0389*rh))!/(4.*pi*rh)        !4*pi*rh factor cancel in gb0
        else
          stop'idfpomr_NE_1'
        endif
        wudt=wwuu+wwdd+wwud+wwdu

        gb=gb*dble(wudt)/zpm**dble(delh+0.4)/(delh+0.4)
     *       *z1*(xmax1-xmin1)/gb0
        if((gb.gt.1.d0.and.ish.ge.2).or.ish.ge.9)write (ifch,*)
     *      'gb-qq,iclpro,zpm,xp,tmax,tmin,xpmax,coef',
     *      gb,iclpro,zpm,xp,tmax,tmin,xpmax,coef_gb0
        gbymax=max(gbymax,gb)
        gbymin=min(gbymin,gb)
        if(gb.gt.1)then   !do not accept these events otherwise jet distribution false !!!!
          if(ish.ge.0)write(ifmt,*)iqq,' val-val gb=',gb
          goto 3
        endif
        if(drangen(gb).gt.gb)goto 3 !rejection

c        print *,'val',ntryrej
        !j=aint((real(xp)+delta/2)/delta) !bg
        !cpt(j+1)=cpt(j+1)+1    !bg

         aks=rangen()*wudt
         if(aks.le.wwuu)then
           if(iclpro.le.2)then
             iq1=1
           elseif(iclpro.eq.3)then
             if(iabs(idptl(ip)).eq.130)then !proj=Kch
               iq1=1
             else !proj=K0
               iq1=2
             endif
           else   !charm
             if(iabs(idptl(ip)).eq.140)then
               iq1=1
             elseif(iabs(idptl(ip)).eq.240)then
               iq1=2
             elseif(iabs(idptl(ip)).eq.340)then
               iq1=3
             else
               iq1=4
             endif
           endif
           iq2=1
         elseif(aks.le.wwuu+wwdd)then
           if(iclpro.eq.2)then
             iq1=2
           elseif(iclpro.eq.1)then
             iq1=-2
           elseif(iclpro.eq.3)then
             iq1=-3
           else
             iq1=-4
           endif
           iq2=2
         elseif(aks.le.wwuu+wwdd+wwud)then
           if(iclpro.le.2)then
             iq1=1
           elseif(iclpro.eq.3)then
             if(iabs(idptl(ip)).eq.130)then !proj=Kch
               iq1=1
             else !proj=K0
               iq1=2
             endif
           else   !charm
             if(iabs(idptl(ip)).eq.140)then
               iq1=1
             elseif(iabs(idptl(ip)).eq.240)then
               iq1=2
             elseif(iabs(idptl(ip)).eq.340)then
               iq1=3
             else
              iq1=4
            endif
          endif
          iq2=2
        else
          if(iclpro.eq.2)then
            iq1=2
          elseif(iclpro.eq.1)then
            iq1=-2
          elseif(iclpro.eq.3)then
            iq1=-3
          else
            iq1=-4
          endif
          iq2=1
        endif

        wpi=xp*wp0       !lc+ for the semihard interaction
        wmi=xm*wm0       !lc- for the semihard interaction
        wp1=(wp0-wpi)
        wm1=(wm0-wmi)
        wp1=wp1*psutz(wp1*wm1,dble(amqt(2)**2),dble(amqt(4)**2))
        wm1=wm1-amqt(2)**2/wp1

      !---------------------------------------------------------
      else        ! sea-sea  val-sea  sea-val   iqq=-1,0,1,2
      !---------------------------------------------------------

        if(ish.ge.4)write(ifch,*)'sea-sea  val-sea  sea-val'
     *                           ,iqq,xmin,xmax
        delss=0.5*(dels(1)+dels(2))  !ctp20170105 to be checked if it is the best value to replace unique dels before
        xmin1=xmin**(delh-delss)
        xmax1=xmax**(delh-delss)
        zp0=sqrt(xmin) !xpmin0 !tp

        !-----------------------
        !rejection normalization
        !-----------------------

        if(iqq.lt.0)then   
          do i=0,3
          do j=0,3
            ee44(i,j)=EsatSeaTil(ipomtype,1,xpmin0 ,qqs(1),i)
     .               *EsatSeaTil(ipomtype,2,xpmin0 ,qqs(2),j)
            if(i.ge.2)ee44(i,j)=ee44(i,j)*wtsplit2*wtrapsat
            if(j.ge.2)ee44(i,j)=ee44(i,j)*wtsplit2*wtrapsat
          enddo
          enddo
        elseif(iqq.eq.0)then   
          do i=0,3
          do j=0,3
            ee44(i,j)=EsatSeaTil(ipomtype,1,xpmin0 ,qqs(1),i)
     .               *EsatSeaTil(ipomtype,2,xpmin0 ,qqs(2),j)
          enddo
          enddo
        elseif(iqq.eq.1)then
          val=EsatValTil(1,zp0,xp1,qqs(1),iclpro,0,uv1,dv1)
          do j=0,3
            ee44(0,j)=0
            ee44(1,j)=val*EsatSeaTil(ipomtype,2,xpmin0 ,qqs(2),j)
            ee44(2,j)=0
            ee44(3,j)=0
          enddo
        elseif(iqq.eq.2)then
          val=EsatValTil(2,zp0,xp2,qqs(2),icltar,0,uv1,dv1)
          do i=0,3
            ee44(i,0)=0
            ee44(i,1)=EsatSeaTil(ipomtype,1,xpmin0 ,qqs(1),i)*val
            ee44(i,2)=0
            ee44(i,3)=0
          enddo
        endif 
        if(iclpro.eq.4)stop'treat iclpro=4 in calcSTYP'
        sss=sngl(smax-pth)
        if(ioTestFact.eq.0.and.iqq.ge.0)then
          call pipolSTYPZ(sss)                 !interpolates sTYPz -> /zeroSTYP/ 
        else
          call setQvalues(qqs(1),qqs(2),qqcut) 
          call calcSTYP(iqq,sss)                !computes sTYPz ->  /zeroSTYP/
        endif
        ! - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
        call calcWW4(iqq,ee44,wwgg,wwgq,wwqg,wwqq) !Complete weights containing CS & Esat
        ! - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -

        if(ioTestFact.ne.0)then !extreme case of xp1=xp2=1, val is very small
          if(iqq.eq.1.or.iqq.eq.2) then               
            frej=0.0005   
          elseif(iqq.eq.0)then
            frej=0.001
          else  !iqq.lt.0
            frej=0.01 
          endif                                    
        else !normal case
          if(iqq.ge.1.and.val.gt.0.d0) then               
            frej=4./(30.*(engy)**0.66/val)
          elseif(iqq.eq.0)then
            frej=0.333
          else !sat pom
            frej=1e3 ! 1.0e7             
            if(qqsat.lt.0.10)frej=frej*0.02
          endif                                    
        endif
 
        wwall0=wwgg+wwgq+wwqg+wwqq
        gb0=wwall0/xmax**delh*coef_gb0
     *     /(r2had(iclpro)+r2had(icltar))*z    !4*pi factor cancel in gb
     *     /frej     
c        if(iqq.ne.0)print*,'gb0',iqq,gb0,xmax,val,wwall0,xmax**delh
c     *                    ,coef_gb0
c       *     ,(r2had(iclpro)+r2had(icltar))*z,frej

        if(ioTestFact.lt.8.and.gb0.le.0d0)then
          write(ifmt,'(a,i5,f8.2,2f8.3)')
     .    'WARNING gb0 LE 0',iqq,sss,wwall0,qqs(1)
          if(iqq.ge.0)then
            call utstop("Problem with gb0 in psahot !&")
          else
            iret=1
            imark=8
            goto 16
          endif
        endif              
        z10=z/(r2had(iclpro)+r2had(icltar))
        val0=val

        !-------------------------
        ! xp,xm proposal/rejection 
        !-------------------------

        ntryrej=0
 41     continue ! <======== 41 =========

        !we use dxp * dxm = dzpm * dxp/xp, so we have a factor at rejection:
        !      Jacob1= 1/xp
        zpm=(xmin1+drangen(xmin1)*(xmax1-xmin1))**dble(1./(delh-delss))
        !factor: Jacob2=(xmax1-xmin1)/(delh-delss)*zpm**(1-delh+delss)

        !kkkkkkkkkkkk!KW1908
        if(ioTestFact.ge.6) zpm = 0.9**2
        if(ioTestFact.ge.8)  gb0 = 1
        !kkkkkkkkkkkkkkkkkkk 

        if(iclpro.eq.4)stop'treat iclpro=4 in calcSTYP'
        if(ioTestFact.eq.0.and.iqq.ge.0)then
          sss=sngl(zpm*ss-pth)
          call pipolSTYPZ(sss)                   !interpolates sTYPz -> /zeroSTYP/
        else
          sss=sngl(zpm*ss-pth)
          call setQvalues(qqs(1),qqs(2),qqcut) 
          call calcSTYP(iqq,sss)                   !computes sTYPz ->  /zeroSTYP/
        endif
        xpmax=psutz(ss,zpm*ss,dble(amqpt**2)) !max x+ for sh pomeron
        xpmin=zpm/xpmax
        xpmin=max(xpmin,xomin)
        if(xumin.gt.0.)xpmax=min(xpmax,zpm/xumin)
        if(iqq.lt.0)then
          !write(*,'(a,i5,5f12.2)')'TEST-xpmax',iqq
          !.  ,xpmin,xpmax,ss,zpm*ss,dble(amqpt**2)
          xpmax=xpsat*1.001    !xp is fixed
          xpmin=xpsat          ! to xpsat
        endif
        !kkkkkkkkkkkkkkkkkkk!KW1909
        if(ioTestFact.ne.0.and.ioTestFact.le.5)then 
          xpmax=1
          xpmin=zpm/xpmax
        endif
        !kkkkkkkkkkkkkkkkkkkkkkkkkk 

        gbyj=1d0
        if(iqq.le.2)then   ! better to use sea quark distribution for all cases
          xp=xpmin*(xpmax/xpmin)**drangen(xpmin) !lc+ share for the hard interaction
          !factor: Jacob3= xp*log(xpmax/xpmin)
          xm=zpm/xp             !lc- share for the hard interaction
          !-------------------------------
          !xp-xm swap 
          !-------------------------------
          !if(...)then
          !  xp=xm
          !  xm=zpm/xp
          !endif
          !-------------------------------
        else
          tmax=dlog((1.d0+dsqrt(1.d0-dsqrt(zpm)))
     *             /(1.d0-dsqrt(1.d0-dsqrt(zpm))))
          tmin=dlog((1.d0+dsqrt(1.d0-xpmax))/(1.d0-dsqrt(1.d0-xpmax)))
          t=(tmin+drangen(tmin)*(tmax-tmin))
          xp=1.d0-((1.d0-dexp(-t))/(1.d0+dexp(-t)))**2 !x+_v proposition
          xm=zpm/xp
          !factor: Jacob3= (tmax-tmin)*sqrt(1.d0-xp)
          gbyj=(tmax-tmin)/xp/log(xpmax/xpmin)
     *     *(1.d0-xp)**(.5)
          !if(iqq.eq.2)then !never happens !!!
          !  xp=xm
          !  xm=zpm/xp
          !endif
        endif
        !if(xm.gt.1.)print*,'TEST',iqq,xm,xp/zpm
        
        !kkkkkkkkkkkkkkkkkkk!KW1908
        if(ioTestFact.ge.6) xp = 0.9
        if(ioTestFact.ge.6) xm = 0.9
        !kkkkkkkkkkkkkkkkkkkkkkkkkk 

        if(iqq.le.0)then
          rh=r2had(iclpro)+r2had(icltar)-slopom*sngl(log(zpm))
          z1=z**(rp/rh)/rh    !4*pi factor cancel in gb0
        elseif(iqq.eq.1)then
          if(idfpomr.eq.1.or.idfpomr.eq.2)then
            rh=r2had(iclpro)+r2had(icltar)-slopom*sngl(log(xm))
            z1=z**(rp/rh)/rh    !4*pi factor cancel in gb0
          else
            stop 'idfpomr_NE_1_2'
          endif
        elseif(iqq.eq.2)then
          if(idfpomr.eq.1.or.idfpomr.eq.3)then
            rh=r2had(iclpro)+r2had(icltar)-slopom*sngl(log(xp))
            z1=z**(rp/rh)/rh    !4*pi factor cancel in gb0
          else
            stop'idfpomr_NE_1_3'
          endif
        endif

        if(iqq.lt.0)then
          do i=0,3
          do j=0,3
            ee44(i,j)=EsatSeaTil(ipomtype,1,xp ,qqs(1),i)
     .               *EsatSeaTil(ipomtype,2,xm ,qqs(2),j)
            if(i.ge.2)ee44(i,j)=ee44(i,j)*wtsplit2*wtrapsat
            if(j.ge.2)ee44(i,j)=ee44(i,j)*wtsplit2*wtrapsat
          enddo
          enddo
        elseif(iqq.eq.0)then
          do i=0,3
          do j=0,3
            ee44(i,j)=EsatSeaTil(ipomtype,1,xp ,qqs(1),i)
     .               *EsatSeaTil(ipomtype,2,xm ,qqs(2),j)
          enddo
          enddo
        elseif(iqq.eq.1)then
          val=EsatValTil(1,xp,xp1,qqs(1),iclpro,1,uv1,dv1)
          do j=0,3
            ee44(0,j)=0
            ee44(1,j)=val*EsatSeaTil(ipomtype,2,xm,qqs(2),j)
            ee44(2,j)=0
            ee44(3,j)=0
          enddo
        elseif(iqq.eq.2)then
          val=EsatValTil(2,xm,xp2,qqs(2),icltar,1,uv1,dv1)
          do i=0,3
            ee44(i,0)=0
            ee44(i,1)=EsatSeaTil(ipomtype,1,xp ,qqs(1),i)*val
            ee44(i,2)=0
            ee44(i,3)=0
          enddo
        endif
        ! - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
        call calcWW4(iqq,ee44,wwgg,wwgq,wwqg,wwqq) ! Complete weights containing CS & Esat
        ! - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -

        wwall=wwgg+wwgq+wwqg+wwqq
        gbyj=gbyj*(wwgg+wwgq+wwqg+wwqq)*z1
        !rejection function factors from trafos and EsaturTil:
        ! Jacob1 * Jacob2 * zpm**(-1-delss) 
        ! =1/xp * (xmax1-xmin1)/(delh-delss)*zpm**(1-delh+delss) * zpm**(-1-delss) 
        gbyj=gbyj/zpm**dble(delh)/(delh-delss)*(xmax1-xmin1)
     .  *log(xpmax/xpmin)
     .  /gb0
        gbymax=max(gbymax,gbyj)
        gbymin=min(gbymin,gbyj)
        if(.not.(gbyj.gt.0..or.gbyj.le.0.))stop'ERROR 08092019 NaN'
        if(gbyj.gt.1d0.and.ioTestFact.le.5)then
        write(ifmt,'(a,f8.2,i5,6f9.2)')'WARNING gbyj',gbyj,iqq !,xp1,xp2
     .                       ,qqs(1),zpm*ss,s2min,s2max !,rapsat
        endif

        if(iqq.ge.0.and.ioTestFact.le.5.and.drangen(gbyj).gt.gbyj
     .    .or.iqq.lt.0.and.gbyj.le.0.d0             !sat accepted when proba non-zero
     .    .or.ioTestFact/2.eq.3.and.gbyj.le.0.d0)then ! 6,7 accepted when proba non-zero (8,9,10 always accepted)
          ntryrej=ntryrej+1
          ktry41=ktry41+1
          nrejtot=nrejtot+1
          if(ntryrej.gt.999999)then
            write(ifmt,'(a,i3,9e10.2)')'WARNING ntryrej for iqq ='
     .     ,iqq,gbyj,gb0,rapsat !see 555
            igo=2
            goto 1 !redo psahot
          endif
          goto 41   !rejection ==============> 41
        endif
        !if(iqq.eq.-1)print*,'GBY',ntryrej,gbymin,gbymax,iqq,rapsat

        !j=aint((real(xp)+delta/2)/delta)            !bg
        !cpt(j+1)=cpt(j+1)+1                      !bg

        wpi=wp0*xp               !lc+ for the hard interaction
        wmi=wm0*xm               !lc+-for the hard interaction
        wp1=wp0-wpi         !remaining lc+ for the soft part
        wm1=wm0-wmi         !remaining lc- for the soft part

        !-------------------------
        !  iq1, iq2
        !-------------------------

        if(iqq.eq.1)then  
          aks=rangen()
          if(aks.le.uv1/(uv1+dv1))then
            if(iclpro.le.2)then
              iq1=1
            elseif(iclpro.eq.3)then
              if(iabs(idptl(ip)).eq.130)then !proj=Kch
                iq1=1
              else !proj=K0
                iq1=2
              endif
            else              !charm
              if(iabs(idptl(ip)).eq.140)then
                iq1=1
              elseif(iabs(idptl(ip)).eq.240)then
                iq1=2
              elseif(iabs(idptl(ip)).eq.340)then
                iq1=3
              else
                iq1=4
              endif
            endif
          else
            if(iclpro.eq.2)then
              iq1=2
            elseif(iclpro.eq.1)then
              iq1=-2
            elseif(iclpro.eq.3)then
              iq1=-3
            else
              iq1=-4
            endif
          endif
        elseif(iqq.eq.2)then
          aks=rangen()
          if(aks.le.uv1/(uv1+dv1))then
            iq2=1
          else
            iq2=2
          endif
        endif         

      !-------------------------------
      endif                     ! iqq (0,1,2,3)
      !-------------------------------

      call setQvalues(qqs(1),qqs(2),qqcut) 
      call calcSTYP(iqq,sss)          !computes sTYP ->  /fiveSTYP/
      call completeSTYP(iqq,ee44)     !computes sTYP = sTYP * ee44 ->  /fiveSTYP/
      if(ish.ge.6)write (ifch,*)'Initial Pt :'
     .  ,sqrt(pth),q2cmin,qqs,ptiprj,ptitrg

      if(ish.ge.8)write (ifch,*)'Momentum check (0)'
     .  ,wp0-wpi-wp1-wp2-wpq(1)-wpq(2),wm0-wmi-wm1-wm2-wmq(1)-wmq(2)

      do njx=nj+1,nj+10
        if(ncj(1,njx).ne.0.or.ncj(2,njx).ne.0)then
          write(ifmt,'(a)')
     .    'ERROR psahot (3) ncj non-zero, wrong initialization'
        endif
      enddo

c      print *,'hard pom',iqq,ntryrej,ktry41,gb0,gbyj
c      if(iqq.eq.0)stop



c##############################################################################################
c##############################################################################################
c        flavor and momenta of end partons of the hard Pomeron
c##############################################################################################
c##############################################################################################

      !~~~~~~~~~~~~~~~~~~
      !if(ioTestFact.ne.0)write(ifmt,'(a,3i7,4f10.3)')
      !.'TestFact iqq ntry ntryhot wp0 wm0 wpi wmi ='
      !.,iqq,ntryrej,ntryhot,wp0, wm0, wpi, wmi
      !~~~~~~~~~~~~~~~~~~

      call timer(iutime)
      tidi=iutime(3)-tiu3 + (iutime(4)-tiu4)*1d-3
      if(ltime)write(ifmt,'(a,i5)')'cpu time psahot 1/5'
     .,nint(tidi*1000)
      tiu3=iutime(3)
      tiu4=iutime(4)


      if(ish.ge.4)write(ifch,*)
     &  'flavor and momenta of end partons of the hard Pomeron'
      ntryx=0
      amqti(1)=amqt(1)
      xxp1poi=xxp1pom
      xyp1poi=xyp1pom
      amqti(2)=amqt(2)
      xxp2poi=xxp2pom
      xyp2poi=xyp2pom
      amqti(3)=amqt(3)
      xxm2poi=xxm2pom
      xym2poi=xym2pom
      amqti(4)=amqt(4)
      xxm1poi=xxm1pom
      xym1poi=xym1pom
      amqpti=amqpt
      xpi=xp
      xmi=xm

      amqt(1)=amqti(1)
      xxp1pom=xxp1poi
      xyp1pom=xyp1poi
      amqt(2)=amqti(2)
      xxp2pom=xxp2poi
      xyp2pom=xyp2poi
      amqt(3)=amqti(3)
      xxm2pom=xxm2poi
      xym2pom=xym2poi
      amqt(4)=amqti(4)
      xxm1pom=xxm1poi
      xym1pom=xym1poi
      amqpt=amqpti
      
      dq2mass=0d0

      wp1i=wp1
      wm1i=wm1
      xp=xpi
      xm=xmi

      call initPartonList(nj0+1,njstr)
      nj=nj0          !initialization for the number of jets
      if(ish.ge.8)write (ifch,*)'6-ww,smin',wpi*wmi,smin,xp,xm,wp1,wm1

      rrr=drangen(wpi)
      jqq=0
      iqc(1)=0
      iqc(2)=0

      call iddeco(icp,jcp)     !save valence quarks in jcp and jct
      call iddeco(ict,jct)

      if(iremn.ge.2)then
        do j=1,2
          do n=1,nrflav
            jcpr(n,j)=jcpref(n,j,ip)   !remnant sea quarks in jcpr and jctr
            jctr(n,j)=jctref(n,j,it)
          enddo
          do n=nrflav+1,nflav
            jcpr(n,j)=0
            jctr(n,j)=0
          enddo
        enddo
      endif
      do i=1,2
        icp1(i)=0
        icp2(i)=0
        icm1(i)=0
        icm2(i)=0
      enddo

      !---------------------------------------------------------------------------
      ! if one or 2 of the 2 initial pomeron ends is sea, choose whether the first
      ! emitted parton is a gluon or a quark to define the hard interaction
      !-------------------------------------------------
      ! iqq=0 sea-sea
      !        jqq=0 : gg interaction after soft emission
      !        jqq=1 : qg interaction after soft emission
      !        jqq=2 : gq interaction after soft emission
      !        jqq=3 : qq interaction after soft emission
      !-------------------------------------------------
      ! iqq=1 val-sea
      !        jqq=0 : qg interaction after soft emission
      !        jqq=1 : qq interaction after soft emission
      !-------------------------------------------------
      ! iqq=2 sea-val
      !        jqq=0 : gq interaction after soft emission
      !        jqq=1 : qq interaction after soft emission
      !-------------------------------------------------
      ! iqq=3 val-val
      !        jqq=3
      !-------------------------------------------------
      ! If the initial parton is a quark, the flavor has to be defined and the
      ! antiflavor remains as string end.
      !----------------------------------------------------------------------------

      if(iqq.eq.1)then
        if(rrr.lt.wwqq/(wwqg+wwqq))jqq=1 
      elseif(iqq.eq.2)then
        if(rrr.lt.wwqq/(wwgq+wwqq))jqq=1 
      elseif(iqq.le.0)then
        rrr=rrr*(wwgg+wwqg+wwgq+wwqq)
        if(rrr.lt.wwqg)then
          jqq=1
        elseif(rrr.lt.wwqg+wwgq)then
          jqq=2
        elseif(rrr.lt.wwqg+wwgq+wwqq)then
          jqq=3
        endif
c      elseif(iqq.lt.0)then
c        jqq=0
      else 
        jqq=3
      endif

      if(ish.ge.6)write(ifch,*)
c      print*,
     .       "Initial parton type:",iqq,jqq,wwgg,wwqq,wwgq,qqs

      if((iqq-1)*(iqq-3).eq.0)then !~~~~~~~ val on projectile side
c valence to start SL cascade
        iqc(1)=iq1                                !proj=particle
        if(iabs(idptl(ip)).eq.1220)iqc(1)=3-iq1      !proj=neutron
        if(idptl(ip).lt.0)iqc(1)=-iqc(1)               !proj=antiparticle
        ifl1=iabs(iqc(1))
c check string ends for mesons :
c switch q<->aq if needed (done randomly in ProSeTyp)
        if(iabs(idptl(ip)).lt.1000)then
          if(iqc(1).gt.0.and.idp1pom.ne.2)then
            idp2pom=idp1pom
            idp1pom=2
          elseif(iqc(1).lt.0.and.idp2pom.ne.2)then
            idp1pom=idp2pom
            idp2pom=2
          endif
        endif
c store flavor of used valence quark in icp1(1) (quark) or icp2(2) (antiquark)
        idsp=100
        if(iqc(1).gt.0)then
          icp1(1)=ifl1
          if(idp1pom.ne.2.and.ioTestFact.eq.0)
     &    call utstop("psahot: Problem with SE (1)!&")
        else
          icp2(2)=ifl1
          if(idp2pom.ne.2.and.ioTestFact.eq.0)
     &    call utstop("psahot: Problem with SE (2)!&")
        endif
      else !~~~~~~~~no valence quark on projectile side
        idsp=idsprj
      endif

      if((iqq-2)*(iqq-3).eq.0)then !~~~~~~~~ val on target side
c valence to start SL cascade
        iqc(2)=iq2             !tar=particle (can not be pion or kaon !)
        if(iabs(idptl(maproj+it)).eq.1220)iqc(2)=3-iq2     !targ=neutron
        if(idptl(maproj+it).lt.0)iqc(2)=-iqc(2)  !targ=antiparticle
        ifl2=iabs(iqc(2))

c store flavor of used valence quark in icm1(1) (quark) or icm2(2) (antiquark)
        idst=100
        if(iqc(2).gt.0)then
          icm1(1)=ifl2
          if(idm1pom.ne.2.and.ioTestFact.eq.0)
     &    call utstop("psahot: Problem with SE (3)!&")
        else
          icm2(2)=ifl2
          if(idm2pom.ne.2.and.ioTestFact.eq.0)
     &    call utstop("psahot: Problem with SE (4)!&")
        endif
      else !~~~~~~~~~~no valence quark on target side
        idst=idstrg
      endif

c determine soft string end flavors icp1,icp2,icm1,icm2

      if(jsplit.eq.1)then
        iret=0
        call fstrfl(jcpr,jctr,jcp,jct,icp1,icp2,icm1,icm2
     *  ,idp1pom,idp2pom
     *  ,idm1pom,idm2pom ,idsp,idst,iret)
        if(iret.ne.0)then
          if(ish.ge.1)
     .    write(ifmt,'(a,i5)')'WARNING psahot: SE flav pb-> exit '
          iret=1
          imark=9
          goto 16
        endif
      else
        iflp=1+3*rangen()
        icp1(1)=10**(6-iflp)
        icp1(2)=0
        icp2(1)=0
        icp2(2)=10**(6-iflp)
        iflt=1+3*rangen()
        icm1(1)=10**(6-iflt)
        icm1(2)=0
        icm2(1)=0
        icm2(2)=10**(6-iflt)
      endif

      !randomize order proj/targ
      ipt1=1
      ipt2=2
      ipto=1
      if(rangen().le.0.5)then
        ipt1=2
        ipt2=1
        ipto=-1
      endif
      
      !-----------------------------------------------------------------
      ! Determine flavor of quarks simultaneously on both sides  KW1808c 
      !-----------------------------------------------------------------

      !detailed values, q in the following means light quark
      laddtype=0
      i=0
      call weightsFlavorPairs(iqq,jqq,xkk,ykk) !xkk,ykk from sTYP

      !-------------------------------
      if(iqq.ne.3)then    !not val-val
      !-------------------------------

      if(iqq.le.0)then  !sea-sea     
        if(jqq.eq.0)then            !gg 
          iqc(1)=0
          iqc(2)=0
          i=1
        elseif(jqq.eq.1)then        !qg
          i=igetRNDindex(3,xkk)
          if(i.eq.1)then    !light
            j=min(max(1,nint(noflit))
     .       ,max(1,int(rangen()*noflit+1))) *(2.*int(.5+rangen())-1.) 
          elseif(i.eq.2)then !charm
            j=4 *(2.*int(.5+rangen())-1.) 
          elseif(i.eq.3)then  !bottom
            j=5 *(2.*int(.5+rangen())-1.) 
          endif
          iqc(1)=j
          iqc(2)=0
        elseif(jqq.eq.2)then        !gq
          i=igetRNDindex(3,xkk)
          if(i.eq.1)then    !light
            j=min(max(1,nint(noflit))
     .       ,max(1,int(rangen()*noflit+1))) *(2.*int(.5+rangen())-1.) 
          elseif(i.eq.2)then !charm
            j=4 *(2.*int(.5+rangen())-1.) 
          elseif(i.eq.3)then  !bottom
            j=5 *(2.*int(.5+rangen())-1.) 
          endif
          iqc(1)=0
          iqc(2)=j
        elseif(jqq.eq.3)then        !qq
          i=igetRNDindex(11,xkk)
          if(i.eq.1)then     ! l l , l al, l l'    (l means light)
            iqc(1) = min(max(1,nint(noflit))
     .       ,max(1,int(rangen()*noflit+1))) *(2.*int(.5+rangen())-1.)
            iqc(2) = min(max(1,nint(noflit))
     .       ,max(1,int(rangen()*noflit+1))) *(2.*int(.5+rangen())-1.)
          elseif(i.eq.2)then  ! l c 
            iqc(1) = min(max(1,nint(noflit))
     .       ,max(1,int(rangen()*noflit+1))) *(2.*int(.5+rangen())-1.)
            iqc(2) =        4          * (2.*int(.5+rangen())-1.) 
          elseif(i.eq.3)then  ! c l  !KW1808e case 2 and 3 need to be treated individually
            iqc(1) =        4          * (2.*int(.5+rangen())-1.) 
            iqc(2) = min(max(1,nint(noflit))
     .       ,max(1,int(rangen()*noflit+1))) *(2.*int(.5+rangen())-1.)
          elseif(i.eq.4)then  ! l b 
            iqc(1) = min(max(1,nint(noflit))
     .       ,max(1,int(rangen()*noflit+1))) *(2.*int(.5+rangen())-1.)
            iqc(2) =        5          * (2.*int(.5+rangen())-1.) 
          elseif(i.eq.5)then  ! b l  !KW1808e case 4 and 5 need to be treated individually
            iqc(1) =        5          * (2.*int(.5+rangen())-1.) 
            iqc(2) = min(max(1,nint(noflit))
     .       ,max(1,int(rangen()*noflit+1))) *(2.*int(.5+rangen())-1.)
          elseif(i.eq.6)then  !  c c
            iqc(1) =        4          * (2.*int(.5+rangen())-1.) 
            iqc(2) = iqc(1) 
          elseif(i.eq.7)then  !  c ac
            iqc(1) =        4          * (2.*int(.5+rangen())-1.) 
            iqc(2) = -iqc(1) 
          elseif(i.eq.8)then  !  b b
            iqc(1) =        5          * (2.*int(.5+rangen())-1.) 
            iqc(2) = iqc(1) 
          elseif(i.eq.9)then  !  b ab
            iqc(1) =        5          * (2.*int(.5+rangen())-1.) 
            iqc(2) = -iqc(1) 
          elseif(i.eq.10)then  !  c b
            iqc(1) =        4          * (2.*int(.5+rangen())-1.) 
            iqc(2) =        5          * (2.*int(.5+rangen())-1.)
          elseif(i.eq.11)then  !  b c
            iqc(1) =        5          * (2.*int(.5+rangen())-1.) 
            iqc(2) =        4          * (2.*int(.5+rangen())-1.)
          else
            stop'ERROR 20082018a'
          endif
        else
          stop'ERROR 19082018a'
        endif
      elseif(iqq.eq.1.or.iqq.eq.2)then !val-sea sea-val
        if(iqq.eq.1)then
          ivvv=iqc(1) !already determined
        elseif(iqq.eq.2)then
          ivvv=iqc(2) !already determined
        endif
        if(jqq.eq.0)then             !vg
          isea=0
          i=1
        elseif(jqq.eq.1)then         !vq
          i=igetRNDindex(5,xkk)
          if(i.eq.1)then      ! l l
            isea = ivvv
          elseif(i.eq.2)then  ! l al
            isea = -ivvv
          elseif(i.eq.3)then  ! l l'
  17        isea = min(max(1,nint(noflit))
     .       ,max(1,int(rangen()*noflit+1))) *(2.*int(.5+rangen())-1.)
            if(isea.eq.ivvv.or.isea.eq.-ivvv)goto 17
          elseif(i.eq.4)then  ! l c  or  c l
            isea =       4          * (2.*int(.5+rangen())-1.) 
          else                ! l b  or  b l
            isea =       5          * (2.*int(.5+rangen())-1.) 
          endif
        else
          stop'ERROR 19082018b'
        endif
        if(iqq.eq.1)then
          iqc(2)=isea
        elseif(iqq.eq.2)then
          iqc(1)=isea
        endif
      else
        stop'ERROR 19082018d'
      endif
      if(i.eq.0)stop'ERROR 28062019b'
      igroup=i
      laddtype=igetRNDindex4(ykk(1,igroup),ykk(2,igroup),ykk(3,igroup))
      if(iqq.lt.0.and.laddtype.ne.1)
     .call utstop("laddtype.ne.1 with iqq<0 !&")
      !if(abs(iqc(1)).eq.4.or.abs(iqc(2)).eq.4)
      !.write(ifmt,*)'ladder end flavors iqq jqq',iqc(1),iqc(2),iqq,jqq
      if(ish.ge.4)then !~~~~~~~~~~~~~~print~~~~~~~~~~~~~~~~~~~~~~~~~
        write(ifmt,'(80a1/80a1)')('#',idx=1,160)
        write(ifmt,'(a,2i3,3x,a,i3,3f8.3,3x,a,2i3)')'Pom Type',iqq,jqq
     .,'LaddType',laddtype,(ykk(k,i),k=1,3),'Flav',iqc(1),iqc(2)
      endif !~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

c in case of gluon emission, compensate intr pt in string end

      if(iqc(1).eq.0)then        
        rrr=drangen(xp)
        xxp1pom=xxp1pom-pxh1*rrr !amqt(1)
        xxp2pom=xxp2pom-pxh1*(1d0-rrr) !amqt(2)
        rrr=drangen(xm)
        xyp1pom=xyp1pom-pyh1*rrr !amqt(1)
        xyp2pom=xyp2pom-pyh1*(1d0-rrr) !amqt(2)
        amqt(1)=sqrt(sngl(xxp1pom**2+xyp1pom**2))
        amqt(2)=sqrt(sngl(xxp2pom**2+xyp2pom**2))
      endif
      if(iqc(2).eq.0)then    
        rrr=drangen(xm)
        xxm1pom=xxm1pom-pxh2*rrr !amqt(4)
        xxm2pom=xxm2pom-pxh2*(1d0-rrr) !amqt(3)
        rrr=drangen(xp)
        xym1pom=xym1pom-pyh2*rrr !amqt(4)
        xym2pom=xym2pom-pyh2*(1d0-rrr) !amqt(3)
        amqt(3)=sqrt(sngl(xxm2pom**2+xym2pom**2))
        amqt(4)=sqrt(sngl(xxm1pom**2+xym1pom**2))
      endif
      amqpt=amqt(1)+amqt(2)+amqt(3)+amqt(4)


c in case of sea quark emission, fix momentum here to take it from the soft part

      !----------------------------------------------------------
      ! wp0        lc+ for the semihard interaction
      ! wm0        lc- for the semihard interaction
      ! wpi=wp0*xp       lc+ for the hard interaction
      ! wmi=wm0*xm       lc+-for the hard interaction
      !
      ! wp1i=wp1=wp0-wpi          remaining lc+ for the soft part
      ! wm1i=wm1=wm0-wmi          remaining lc- for the soft part
      !
      ! wpq(1)                 emitted parton
      ! wmq(1)
      !   wm1i=wm1i-wmq(1)            remove from
      !   wp1i=wp1i-wpq(1)            soft part
      !
      ! wmq(2)                 emitted parton
      ! wpq(2)
      !   wp1i=wp1i-wpq(2)            remove from
      !   wm1i=wm1i-wmq(2             soft part
      !         _____    _____ 
      !              | g |     soft         ioj = 7 8
      !              | l |      
      !              v u ^      
      !              |   |      
      !              |...|_____ emission     ioj = -7 -8
      !         hard |         
      !---------------------------------------------------------      

      xg=5./sqrt(wp0*wm0)

      do ipt=ipt1,ipt2,ipto !~~~~~~~~~~~~~~loop over proj / targ ~~~~~~~

      if(ipt.eq.1)then      !~~~~~projectile

        if(jqq.eq.3.or.(jqq.eq.1.and.iqq.ne.1))then !~~~~~quark initiate space like cascade

          q2mass(1)=0.
          if(abs(iqc(1)).eq.4)then
            q2mass(1)=qcmass**2
          elseif(abs(iqc(1)).eq.5)then
            q2mass(1)=qbmass**2
          endif
          dq2mass=dq2mass+dble(q2mass(1))
          q2mass(1)=pxh1**2+pyh1**2+q2mass(1)
            
          wpq(1) = exp(rangen()*log(xg)) * wp1i   ! 1/x from xg to 1
          wmq(1) = exp(rangen()*log(xg)) * wm1i   !
          wm1i=wm1i-wmq(1)   !remove momentum- from soft part of other side
          wp1i=wp1i-wpq(1)   !remove momentum+ from soft part of other side

        endif                   !~~~~~quark initiate space like cascade    

      else                      !~~~~~target
 
        if(jqq.ge.2.or.(jqq.eq.1.and.iqq.eq.1))then !~~~~~quark initiate space like cascade

          q2mass(2)=0.
          if(abs(iqc(2)).eq.4)then
            q2mass(2)=qcmass**2
          elseif(abs(iqc(2)).eq.5)then
            q2mass(2)=qbmass**2
          endif
          dq2mass=dq2mass+dble(q2mass(2))
          q2mass(2)=pxh2**2+pyh2**2+q2mass(2)
          
          wmq(2) = exp(rangen()*log(xg)) * wm1i  
          wpq(2) = exp(rangen()*log(xg)) * wp1i 
          wp1i=wp1i-wpq(2)
          wm1i=wm1i-wmq(2)

        endif                                        !~~~~~quark initiate space like cascade

      endif                  !~~~~~target
      enddo                  !randomized order

      !momenta now fixed
      wp1=wp1i !remaining lc+ for the soft part
      wm1=wm1i !remaining lc- for the soft part

      if(ish.ge.8)write (ifch,*)'Momentum check (1)'
     .  ,wp0-wpi-wp1-wp2-wpq(1)-wpq(2),wm0-wmi-wm1-wm2-wmq(1)-wmq(2)
         
      if(iqq.le.0)then          ! ----------- sea-sea ------------- iqq=0

        if(ioTestFact.le.7)then
        !share momentum of soft part wp1 and wm1 in 2 strings with wp1,wm1 and wp2,wm2
        call pslcsh(wp1,wm1,wp2,wm2,amqt,iqq)
        if(ntrylcsh.gt.99)then
          if(ish.ge.1)
     .      write(ifmt,'(a,2e10.3,i5)')
     .      'WARNING psahot: return to start, after pslcsh',wp2,wm2,iqq
          igo=7
          goto 1 !redo psahot
        endif
        endif
        
      else                      ! -------------- val-sea  sea-val --- iqq=1,2

        if(ish.ge.8)write (ifch,*)'q_sea mass check',wp1*wm1,amqpt**2
        if(iqq.eq.1)then        !valence quark-gluon hard interaction
          amq1=amqt(3)**2
          s24=(amqt(2)+amqt(4))**2
        else
          amq1=amqt(1)**2
          s24=(amqt(2)+amqt(4))**2
          xp=xm
          xm=zpm/xp
        endif
        x1max=psutz(wp1*wm1,dble(amq1),dble(s24))
        x1min=dble(amq1)/x1max/wp1/wm1
        t1min=(1.d0/x1max-1.d0)
        t1max=(1.d0/x1min-1.d0)
        igo=8
        if(t1max.le.0.and.ioTestFact.le.7)goto 1 !redo psahot
        if(ioTestFact.le.7)then
 5        t1x=t1min*(t1max/t1min)**drangen(t1min)
          if(ish.ge.8)write (ifch,*)'t1x,t1min,t1max',t1x,t1min,t1max
          xq1=1.d0/(1.d0+t1x)
          if(drangen(xq1).gt.(xq1*(1.d0-xq1))**(1.+(-alpqua(3-iqq))))
     .    goto 5
        else
          xq1=0
        endif
        if(iqq.eq.1)then         !valence quark-gluon hard interacti
          wm2=wm1*(1.d0-xq1)   !share wm1
          wm1=wm1*xq1
          wp1=wp1-dble(amq1)/wm1   !remove momentum given to the string
          if(ish.ge.8)write (ifch,*)'q_sea+ mass check',
     *        wp1*wm2,s24
          wp1=wp1*psutz(wp1*wm2,dble(amqt(2)**2),dble(amqt(4)**2))
          wm2=wm2-dble(amqt(2))**2/wp1
        else
          wp2=wp1*(1.d0-xq1)
          wp1=wp1*xq1
          wm1=wm1-dble(amq1)/wp1
          if(ish.ge.8)write (ifch,*)'q_sea- mass check',
     *        wp2*wm1,s24
          wm1=wm1*psutz(wm1*wp2,dble(amqt(4)**2),dble(amqt(2)**2))
          wp2=wp2-amqt(4)**2/wm1
        endif

      endif !iqq

      !-------------------------------
      else  !val-val
      !-------------------------------

      !iqc(1),iqc(2) already determinded
      if(iqc(1).eq.iqc(2))then
        i=1 
      elseif(iqc(1).eq.-iqc(2))then
        i=2
      else
        i=3
      endif 
      igroup=i
      laddtype=igetRNDindex4(ykk(1,igroup),ykk(2,igroup),ykk(3,igroup))

      !-------------------------------
      endif  !not val-val / val-val
      !-------------------------------

      if(ish.ge.6)then
        write (ifch,*)'Momentum check (2)'
     .  ,wp0-wpi-wp1-wp2-wpq(1)-wpq(2)
     .  ,wm0-wmi-wm1-wm2-wmq(1)-wmq(2),amqt
      endif
      
      !~~~~~~~~~~~~~~~~~~~~~~~~
      !if(ioTestFact.ne.0)write(ifmt,'(a,3i7,2f10.3,2(/10x,6f10.3))')
      !.'TestFact1 iqq ntry ntryhot w+ w-',iqq,ntryrej,ntryhot
      !.,wp0-wpi-wp1-wp2-wpq(1)-wpq(2), wm0-wmi-wm1-wm2-wmq(1)-wmq(2)
      !.,wp0,wpi,wp1,wp2,wpq(1),wpq(2), wm0,wmi,wm1,wm2,wmq(1),wmq(2)
      !~~~~~~~~~~~~~~~~~~~~~~~~

      !kkkkkkkkkkkkkkkkkkkkkkkkkkkkkkk
      if(ioTestFact.ge.8)then
      wp1=0.001
      wp2=0.001
      wm1=0.001
      wm2=0.001
      endif
      !kkkkkkkkkkkkkkkkkkkkkkkkkkkkkkk


c####################################################################################
c####################################################################################
c    store TIME-LIKE partons related to Pomeron end ( into /psar29/ ) 
c               concerns momenta and  color connections
c####################################################################################
c####################################################################################

      !------------------------------------------------------
      ! parton properties, for final (TL) partons, index nj, in  /psar29/
      !   eqj(1:4,nj) - 4-momentum (qgs) ncj(1,nj)
      !   bxj(1:4,nj) - coordinates of formation point;\
      !   iqj(k) - parton id 
      ! jq = "color orientation"                                
      ! ncj(1:2,nj) color connections for TL to TL partons
      ! ncc(1:2,1) color connection to TL for current SL (proj) 
      ! ncc(1:2,2) color connection to TL for current SL (targ)
      !-------------------------------------------------------

      do njx=nj+1,nj+10
        if(ncj(1,njx).ne.0.or.ncj(2,njx).ne.0)then
          write(ifmt,'(a)')
     .    'ERROR psahot (4) ncj non-zero, wrong initialization'
        endif
      enddo


      do ipt=ipt1,ipt2,ipto                                 !proj/targ
      if(ipt.eq.1)then                                      !~~~~~~~~~~~~~~~projectile~~~~~~~~~~~~~~~~

       njii=nj

       if((iqq-1)*(iqq-3).eq.0)then                         !~~~~~valence projectile side


c define the second string end from sea
        if(iqc(1).gt.0)then
          ifl=idtra(icp2,1,0,2)
        elseif(iqc(1).lt.0)then
          ifl=idtra(icp1,1,0,2)
        else
          call utstop('No quark for hard Pomeron+ in psahot!&')
        endif

        if(ish.ge.5)write(ifch,*)'flavor vq+,sqb+',iqc(1)
     &                                            ,ifl

        nj=nj+1
        iodur(nj)=0                             !ala
        chj(nj)='proj'
c put the sea quark (or diquark) into the parton list as string end
        if(abs(ifl).eq.3)then
          iqj(nj)=ifl*4/3
        elseif(abs(ifl).ge.4)then
          iqj(nj)=ifl*10
        else
          iqj(nj)=ifl
        endif
        ioj(nj)=7
        ityj(nj)=30

        eqj(1,nj)=.5*sngl(wp1+dble(amqt(2))**2/wp1)
        eqj(2,nj)=wp1-eqj(1,nj)
        eqj(3,nj)=xxp2pom
        eqj(4,nj)=xyp2pom
        if(ish.ge.8)write (ifch,*)'q_v+ mass check',(eqj(1,nj)-
     *    eqj(2,nj))*(eqj(1,nj)+eqj(2,nj))-eqj(3,nj)**2-eqj(4,nj)**2
        eqj(1,nj)=sqrt(eqj(2,nj)**2+eqj(3,nj)**2+eqj(4,nj)**2)
        ncc(1,1)=nj
        ncc(2,1)=0
        call eqjcheck(1,nj)

       else      ! (iqq <= 0 or 2)                            !~~~~~sea projectile side

c define the 2 string ends from sea

        iflq=idtra(icp1,1,0,2)
        iflqb=idtra(icp2,1,0,2)

        if(ish.ge.5)write(ifch,*)'flavor sq+,sqb+',iflq,iflqb

        nj=nj+1
        iodur(nj)=0                             !ala
        chj(nj)='proj'
c put the sea quarks (or diquarks) into the parton list as string end
        if(abs(iflqb).eq.3)then
          iqj(nj)=iflqb*4/3
        elseif(abs(iflqb).ge.4)then
          iqj(nj)=iflqb*10
        else
          iqj(nj)=iflqb
        endif
        ioj(nj)=7
        ityj(nj)=30
        if(abs(iflq).eq.3)then
          iqj(nj+1)=iflq*4/3
        elseif(abs(iflq).ge.4)then
          iqj(nj+1)=iflq*10
        else
          iqj(nj+1)=iflq
        endif
        ioj(nj+1)=7
        ityj(nj+1)=30
        !print*,'ENDPARTONS',nptl,ncolp,ntrysplit, nj,iqj(nj),iqj(nj+1)

        eqj(1,nj)=.5*sngl(wp1+dble(amqt(1))**2/wp1)
        eqj(2,nj)=wp1-eqj(1,nj)
        eqj(3,nj)=xxp1pom
        eqj(4,nj)=xyp1pom
        if(ish.ge.8)write (ifch,*)'q_s1+ mass check',(eqj(1,nj)-
     *    eqj(2,nj))*(eqj(1,nj)+eqj(2,nj))-eqj(3,nj)**2-eqj(4,nj)**2
        eqj(1,nj)=sqrt(eqj(2,nj)**2+eqj(3,nj)**2+eqj(4,nj)**2)

        eqj(1,nj+1)=.5*sngl(wp2+dble(amqt(2))**2/wp2)
        eqj(2,nj+1)=wp2-eqj(1,nj+1)
        eqj(3,nj+1)=xxp2pom
        eqj(4,nj+1)=xyp2pom
        nj=nj+1
        iodur(nj)=0                             !ala
        chj(nj)='proj'
        if(ish.ge.8)write (ifch,*)'q_s2+ mass check',(eqj(1,nj)-
     *    eqj(2,nj))*(eqj(1,nj)+eqj(2,nj))-eqj(3,nj)**2-eqj(4,nj)**2
        eqj(1,nj)=sqrt(eqj(2,nj)**2+eqj(3,nj)**2+eqj(4,nj)**2)

        call eqjcheck(2,nj)
        call eqjcheck(3,nj+1)

        if(jqq.eq.0.or.(iqq.le.0.and.jqq.eq.2))then          !~~~~~gluon ini SL cascade projectile

          ncc(1,1)=nj-1
          ncc(2,1)=nj

        else                                                 !~~~~~quark ini SL cascade projectile

          ! add the corresponding anti-parton in 
          !the list to compensate the evolving one
          !-------------------------------!
          !    ubar _____  g  _____ d     !
          !              | l |            !  color 
          !              | u |            !  flow
          !              v o ^            !
          !              | n |            !  NOT
          !       s _____|...|            !  flavor
          !                  |            !  flow
          !                  |            !
          !                 sbar          !
          !-------------------------------!
          nj=nj+1
          iodur(nj)=0                             !ala
          chj(nj)='proj'
          call idtrafo_KW_SO( -iqc(1) , iqj(nj) , 1 ) !compute iqj(nj)
          ioj(nj)=-7
          ityj(nj)=30

          eqj(1,nj)=.5*sngl(wpq(1)+dble(q2mass(1))/wpq(1))
          eqj(2,nj)=-eqj(1,nj)+sngl(wpq(1))
          !compensate intr pt ==> already commented ==> REMOVED after 4.0.2.q6
          if(idpomr.lt.0)then !sat pom !to avoid strange correlation
           phi=2*pi*rangen() 
           eqj(3,nj)=ptemi*cos(phi)
           eqj(4,nj)=ptemi*sin(phi)
          elseif(idpomr.ge.0)then !nor pom
           eqj(3,nj)=-pxh1
           eqj(4,nj)=-pyh1
          endif
          call eqjcheck(4,nj)
          if(ish.ge.8)write (ifch,*)'q_s3+ mass check',(eqj(1,nj)-
     *      eqj(2,nj))*(eqj(1,nj)+eqj(2,nj))-eqj(3,nj)**2-eqj(4,nj)**2
     *     ,eqj(1,nj),wpq(1),wmq(1),q2mass(1),iqj(nj),nj
          eqj(1,nj)=sqrt(eqj(2,nj)**2+eqj(3,nj)**2+eqj(4,nj)**2
     *                   +q2mass(1))

        if(iqc(1).gt.0)then
            ncj(1,nj)=nj-1
            ncj(1,nj-1)=nj
            ncj(2,nj)=0
            ncj(2,nj-1)=0
            ncc(1,1)=nj-2
            ncc(2,1)=0
          else
            ncj(1,nj)=nj-2
            ncj(1,nj-2)=nj
            ncj(2,nj)=0
            ncj(2,nj-2)=0
            ncc(1,1)=nj-1
            ncc(2,1)=0
          endif
        endif                         !~~~~~  

       endif                          !~~~~~
       !do ij=njii+1,nj
       !write(ifmt,*)'ncj P',iqj(ij),ij,(ncj(j,ij),j=1,2),iqq,jqq
       !enddo

      else                                                   !~~~~~~~~~~~~~~~~~target~~~~~~~~~~~~~~~
 
       njii=nj

       if((iqq-2)*(iqq-3).eq.0)then                          !~~~~~valence target side

c define the second string end from sea
        if(iqc(2).gt.0)then
          ifl=idtra(icm2,1,0,2)
        elseif(iqc(2).lt.0)then
          ifl=idtra(icm1,1,0,2)
        else
          call utstop('No quark for hard Pomeron- in psahot!&')
        endif

        if(ish.ge.5)write(ifch,*)'flavor vq-,sqb-',iqc(2)
     &                                            ,ifl

        nj=nj+1
        iodur(nj)=0                             !ala
        chj(nj)='targ'
c put the sea quark (or diquark) into the parton list as string end
        if(abs(ifl).eq.3)then
          iqj(nj)=ifl*4/3
        elseif(abs(ifl).ge.4)then
          iqj(nj)=ifl*10
        else
          iqj(nj)=ifl
        endif
        ioj(nj)=8
        ityj(nj)=30

        eqj(1,nj)=.5*sngl(wm1+dble(amqt(4))**2/wm1)
        eqj(2,nj)=eqj(1,nj)-sngl(wm1)
        eqj(3,nj)=xxm1pom
        eqj(4,nj)=xym1pom
        if(ish.ge.8)write (ifch,*)'q_v- mass check',(eqj(1,nj)
     *    +eqj(2,nj))*(eqj(1,nj)-eqj(2,nj))-eqj(3,nj)**2-eqj(4,nj)**2
        eqj(1,nj)=sqrt(eqj(2,nj)**2+eqj(3,nj)**2+eqj(4,nj)**2)
        ncc(1,2)=nj
        ncc(2,2)=0
        call eqjcheck(5,nj)

       else        ! (iqq = 0 or 1)                         !~~~~~~ sea target side

c define the 2 string ends from sea

        iflq=idtra(icm1,1,0,2)
        iflqb=idtra(icm2,1,0,2)


        if(ish.ge.5)write(ifch,*)'flavor sq-,sqb-',iflq,iflqb

        nj=nj+1
        iodur(nj)=0                             !ala
        chj(nj)='targ'
c put the sea quarks (or diquarks) into the parton list as string end
        if(abs(iflqb).eq.3)then
          iqj(nj)=iflqb*4/3
        elseif(abs(iflqb).ge.4)then
          iqj(nj)=iflqb*10
        else
          iqj(nj)=iflqb
        endif
        ioj(nj)=8
        ityj(nj)=30
        if(abs(iflq).eq.3)then
          iqj(nj+1)=iflq*4/3
        elseif(abs(iflq).ge.4)then
          iqj(nj+1)=iflq*10
        else
          iqj(nj+1)=iflq
        endif
        ioj(nj+1)=8
        ityj(nj+1)=30
        !print*,'ENDPARTONS',nptl,ncolp,ntrysplit, nj,iqj(nj),iqj(nj+1)

        eqj(1,nj)=.5*sngl(wm1+dble(amqt(3))**2/wm1)
        eqj(2,nj)=eqj(1,nj)-sngl(wm1)
        eqj(3,nj)=xxm2pom
        eqj(4,nj)=xym2pom
        if(ish.ge.8)write (ifch,*)'q_s1- mass check',(eqj(1,nj)-
     *    eqj(2,nj))*(eqj(1,nj)+eqj(2,nj))-eqj(3,nj)**2-eqj(4,nj)**2
        eqj(1,nj)=sqrt(eqj(2,nj)**2+eqj(3,nj)**2+eqj(4,nj)**2)

        eqj(1,nj+1)=.5*sngl(wm2+dble(amqt(4))**2/wm2)
        eqj(2,nj+1)=eqj(1,nj+1)-sngl(wm2)
        eqj(3,nj+1)=xxm1pom
        eqj(4,nj+1)=xym1pom
        nj=nj+1
        iodur(nj)=0                             !ala
        chj(nj)='targ'
        if(ish.ge.8)write (ifch,*)'q_s2- mass check',(eqj(1,nj)-
     *    eqj(2,nj))*(eqj(1,nj)+eqj(2,nj))-eqj(3,nj)**2-eqj(4,nj)**2
        eqj(1,nj)=sqrt(eqj(2,nj)**2+eqj(3,nj)**2+eqj(4,nj)**2)
        call eqjcheck(6,nj)
        call eqjcheck(7,nj+1)

        if(jqq.eq.0.or.(iqq.le.0.and.jqq.eq.1))then         !~~~~~gluon ini SL cascade target

          ncc(1,2)=nj-1
          ncc(2,2)=nj

        else                                                !~~~~~quark ini SL cascade target

          ! add the corresponding anti-parton in 
          !the list to compensate the emitted one
          !-------------------------------!
          !                  |            !  color
          !                  |            !  flow
          !                  |s           !
          !    sbar _____ ...|            !  NOT
          !              | g |            !   
          !              | l |            !  flavor
          !              v u ^            !  color
          !              | o |            !  flow
          !       u _____| n |_____ dbar  !
          !-------------------------------!
          nj=nj+1
          iodur(nj)=0                            
          chj(nj)='targ'
          call idtrafo_KW_SO( -iqc(2) , iqj(nj) , 1 ) !compute iqj(nj)
          ioj(nj)=-8
          ityj(nj)=30

          eqj(1,nj)=.5*sngl(wmq(2)+dble(q2mass(2))/wmq(2))
          eqj(2,nj)=eqj(1,nj)-sngl(wmq(2))
          !compensate intr pt ==> already commented ==> REMOVED after 4.0.2.q6
          if(idpomr.lt.0)then !sat pom !to avoid strange correlation
           phi=2*pi*rangen() 
           eqj(3,nj)=ptemi*cos(phi)
           eqj(4,nj)=ptemi*sin(phi)
          elseif(idpomr.ge.0)then !nor pom
           eqj(3,nj)=-pxh2
           eqj(4,nj)=-pyh2
          endif
          call eqjcheck(8,nj)
          if(ish.ge.8)write (ifch,*)'q_s3- mass check',(eqj(1,nj)-
     *       eqj(2,nj))*(eqj(1,nj)+eqj(2,nj))-eqj(3,nj)**2-eqj(4,nj)**2
     *      ,eqj(1,nj),wpq(2),wmq(2),q2mass(2),iqj(nj),nj
          eqj(1,nj)=sqrt(eqj(2,nj)**2+eqj(3,nj)**2+eqj(4,nj)**2
     *                   +q2mass(2))
          
          if(iqc(2).gt.0)then
            ncj(1,nj)=nj-1
            ncj(1,nj-1)=nj
            ncj(2,nj)=0
            ncj(2,nj-1)=0
            ncc(1,2)=nj-2
            ncc(2,2)=0
          else
            ncj(1,nj)=nj-2
            ncj(1,nj-2)=nj
            ncj(2,nj)=0
            ncj(2,nj-2)=0
            ncc(1,2)=nj-1
            ncc(2,2)=0
          endif

        endif                         !~~~~~

       endif                          !~~~~~~

      !do ij=njii+1,nj
      !write(ifmt,*)'ncj 2',iqj(ij),ij,(ncj(j,ij),j=1,2)
      !enddo

      endif                                                !~~~~~~~~~~~~~~~~~target~~~~~~~~~~~~~~~~
      enddo                                                !proj/targ

      do ij=1,nj
        call psmirror1(eqj(2,ij))
      enddo

      if(ish.ge.4)then !~~~~~~~~~~~~~~print~~~~~~~~~~~~~~~~~~~~~~~~~
        write(ifmt,'(80a1)')('-',i=1,80)
        ncount=ncount+1
        write(ifmt,'(a,i4)')'TL list before evolution   count',ncount
        do ij=1,nj
        call idtrafo_KW_SO( iKW , iqj(ij) , -1 ) 
        write(ifmt,'(i4,3x,3i4,f10.2,4x,a4)')iKW
     .  ,ncj(1,ij),ij,ncj(2,ij),eqj(2,ij),chj(ij)
        enddo
        write(ifmt,'(80a1)')('=',i=1,80)
      endif !~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~



c#######################################################################################
c#######################################################################################
c            HARD       Born first, then backward evolution
c#######################################################################################
c#######################################################################################

      call timer(iutime)
      tidi=iutime(3)-tiu3 + (iutime(4)-tiu4)*1d-3
      if(ltime)write(ifmt,'(a,i10)')'cpu time psahot 2/5'
     .,nint(tidi*1000)
      tiu3=iutime(3)
      tiu4=iutime(4)

      jj=1
 
      !-----------------------------------------
      !  iqc(1) , iqc(2)    flavor of end partons 
      !-----------------------------------------

      nj00=nj
      si=wpi*wmi-pth     !total energy squared for the hard
      wpiwmi_save=wpi*wmi
      pth4_save=pth

      !~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~! 
      ! si remains EVERYWHERE the Mandelstam s !!!! !
      !       NOT si=si-dble(q2mass)                !
      !~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~!
 
      if(ish.ge.7)write (ifch,*)'si,wpi,wmi',si,wpi,wmi

      !---------------------------------------------
      !   ept(1) = energy of ladder  (ept will be updated during emissons)
      !   ept(2) = pz of ladder       (4-vector 
      !   ept(3) = px of ladder       SO convention)
      !   ept(4) = py of ladder
      !---------------------------------------------

      ept(1)=.5d0*(wpi+wmi)   
      ept(2)=.5d0*(wpi-wmi)   
      ept(3)=(pxh1+pxh2)      
      ept(4)=(pyh1+pyh2)  

      do i=1,4
        pLadder(i)=ept(i) !pladder unchanged  
      enddo

      if(ish.ge.6)then
        esum=ept(1)
        psum1=ept(3)
        psum2=ept(4)
        psum3=ept(2)
        write(ifch,*)iqc,ept
        do inj=1,nj
          write(ifch,*)inj,iqj(inj),eqj(1,inj),eqj(2,inj)
     .                    ,eqj(3,inj),eqj(4,inj)
          esum=esum+eqj(1,inj)
          psum1=psum1+eqj(3,inj)
          psum2=psum2+eqj(4,inj)
          psum3=psum3+eqj(2,inj)
        enddo
        write(ifch,*)'4-Momentum check (1)',esum,pptl(4,iptl),
     .  esum/pptl(4,iptl),psum3,pptl(3,iptl),
     .  (psum1**2+psum**2+(esum+psum3)*(esum-psum3))/ss,iqq,jqq
      endif

      nprod=0 !bg ga

      do njx=nj+1,nj+10
        if(ncj(1,njx).ne.0.or.ncj(2,njx).ne.0)then
          write(ifmt,'(a)')
     .    'ERROR psahot (5) ncj non-zero, wrong initialization'
        endif
      enddo

      !---check color connection table
      do njx=nj+1,nj+10
        if(ncj(1,njx).ne.0.or.ncj(2,njx).ne.0)then
          write(ifmt,'(a)')
     .    'ERROR psahot (6) ncj non-zero, wrong initialization'
        endif
      enddo

      pt2=ept(3)**2+ept(4)**2        !pt squared  of ladder
      pt=sqrt(sngl(pt2))             !pt of ladder

      !---------------------------------------------------------
      ! wpt(1) ... lc+ of ladder
      ! wpt(2) ... lc- of ladder
      !---------------------------------------------------------

      wpt(1)=ept(1)+ept(2)           
      wpt(2)=ept(1)-ept(2)       

      sis=sngl(si)              !energy of the ladder
      sispth=sis+pth
        
      !do i=1,50
      !   x=1./50*i 
      !   print*,'psevi', x, psevi(10.0,100.00,dble(x),1,1)
      !enddo
      !stop'CHECK psevi'   
      !do i=1,20   
      !  x=i*0.05
      !  a=EsaturGluonTil(dble(x),1.5,1)/x**(dels(1))
      !  b=EsaturQuarkTil(dble(x),1.5,1,1)/x**(dels(1))
      !  c=EsaturValTil(1,dble(x),1d0,1.5,1,2)/x**(dels(1))
      !  write(ifmt,'(a,f6.2,3x,3f8.3)')'Esatur',x,a,b,c
      !enddo
      !stop'CHECK Esatur' 



      !~~~~~~~~~~~~~~~~~~~
      !if(ioTestFact.ne.0)write(ifmt,'(a,3i7,5f10.3)')
      !.'TestFact1 iqq ntry ntryhot sqrt(s) pLadder'
      !.,iqq,ntryrej,ntryhot,sqrt(sis),ept
      !~~~~~~~~~~~~~~~~~~~
     


      !KKKKKKKKKKKKKKKKKKKKKKKKKKKKKKKKKKKKKKKKKKKKKKKKKKKKKKKKKKK
      !KKKKKKKKKKKKKKKKKKKKKKKKKKKKKKKKKKKKKKKKKKKKKKKKKKKKKKKKKKK
      !KK                                                       KK 
      !KK  KW1908                                               KK 
      !KK     Completely new procedure with backward evolution  KK
      !KK                                                       KK
      !KKKKKKKKKKKKKKKKKKKKKKKKKKKKKKKKKKKKKKKKKKKKKKKKKKKKKKKKKKK
      !KKKKKKKKKKKKKKKKKKKKKKKKKKKKKKKKKKKKKKKKKKKKKKKKKKKKKKKKKKK

      !MC settings 

      ntryBORN=1000000 
      fscale4=0.12*qqcut**2    !scale for rejection for Born generation laddtype=4
      fscale23=0.04*qqcut**2     !scale for rejection for Born generation laddtype=2,3
      fscale1=1.5 !0.3       !scale for rejection for Born generation laddtype=1
      ntryBWE=100000
      if(iqq.le.0)then
        fscaleBWE=0.005         !scale for rejection for backward evolution
      else
        fscaleBWE=0.002         !scale for rejection for backward evolution
      endif

      ipm=min(idpomr,0) ! -1 (SAT) or 0 (NOR)

      !KKKKKKKKKKKKKKKKKKKKKKKKKKKKKKKKKKKKKKKKKKKKKKKKKKKKKKKKKKK
      !KK         Generating BORN                               KK  
      !KKKKKKKKKKKKKKKKKKKKKKKKKKKKKKKKKKKKKKKKKKKKKKKKKKKKKKKKKKK
      !KK    for given end-of-ladder flavors mEOL, nEOL         KK 
      !KK   and given ladder four momentum ept (SO convention)  KK
      !KK       we generate tBorn, xmBorn, xpBorn,              KK 
      !KK    and flavors mIB,nIB (Born in) ,mOB,nOB (Born out)  KK 
      !KKKKKKKKKKKKKKKKKKKKKKKKKKKKKKKKKKKKKKKKKKKKKKKKKKKKKKKKKKK
       
      iiiter=1 
      do iiiTest=1,iiiter ! >1 for testing

      j2=iqc(jj) 
      l2=iqc(3-jj)

      !kkkkkkkkkkkkkkkkkkkkkkkkkkkkkKW1909
      !kkkkkkkkkkkkkkkkkkkkkkkkkkkkk
      if(ioTestFact.ne.0)then
        if(ioTestFact.le.5)then
c          modeDeltaFzero=0          !tp in psaini
          if(laddTestFact.gt.0)laddtype=laddTestFact
          pLadder(3)=0 
          pLadder(4)=0
        endif 
        if(ioTestFact.ge.6)then
c          modeDeltaFzero=1 !fzero ~ delta(1-x)           !tp in psaini
          q2cmin(1)=q2nmin
          q2cmin(2)=q2nmin
          qqs(1)=q2nmin
          qqs(2)=q2nmin
          laddtype=laddTestFact
          pLadder(3)=0 
          pLadder(4)=0
        endif
        if(ioTestFact.ge.8)then
          sis=engy**2 *0.9**2
          pLadder(1)=engy *0.9
          do i=2,4
            pLadder(i)=0
          enddo
        endif
        if(ioTestFact.ge.8.and.ioTestFact.le.9)then
          j2=0
          l2=0
        endif
        if(ioTestFact.eq.10)then
          j2=1
          l2=0
        endif
        if(ioTestFact.eq.11)then
          j2=1
          l2=1
        endif
      endif
      !kkkkkkkkkkkkkkkkkkkkkkkkkkkkk
      !kkkkkkkkkkkkkkkkkkkkkkkkkkkkk

      qmin(1)=qqs(1)            
      qmin(2)=qqs(2)
      q1=qmin(jj)
      q2=qmin(3-jj)
      qqcut=max(qmin(1),qmin(2)) !virtuality cutoff
      call eolpib( j2 , l2 , jl ) 
      call eolpi( jl , j1 , l1 )

      yLadder=0.5*log((pLadder(1)+pLadder(2))/(pLadder(1)-pLadder(2)))

      call eol44(jl,ifgr1,ifgr2) 
      ee44x=ee44(ifgr1,ifgr2)    !divide sTYP by ee44x
      if(ee44x.le.1e-10)then
        ncount6=ncount6+1 
        if(ncount6.eq.1.or.ncount6.eq.10.or.ncount6.eq.100
     .    .or.ncount6.eq.1000)then
          lcount6=log10(1.*ncount6)+1
          write(ifmt,'(a,i1,a,e10.3)')'WARNING ',lcount6
     .    ,' psahot ee44x =',ee44x
        endif 
        ee44x=max(ee44x,1.e-10)
      endif

      !call checkCS(q1,q2,qqcut,sis,j1,l1 !checks interpolation quality (slows program down)
      !.,sTYP(1,jl)/ee44x,sTYP(2,jl)/ee44x
      !.,sTYP(3,jl)/ee44x,sTYP(4,jl)/ee44x
      !.,sTPex,ihit,0) !1 to force print
 
  12  continue
      if(iiiter.gt.1)then !possibility to artificially changing values
        !laddtype=4 
        !jl=1+rangen()*20
        !call eolpi( jl , j1 , l1 )
        !call eolpic( j1 , l1 , j2 , l2 ) 
      endif

      ml=jl
      mEOL=j2 !end of ladder parton 1
      nEOL=l2 !end of ladder parton 2
      if(idpomr.ge.0)then
        call getLimitsEoL(jl,sis,qqcut , smn,tmn,tmx,iret)
      else !sat pom -> compute tmn for klas 4, to allow gg->ccar <-------------
        call getBornKin(ifav,sis,pt2xx,smn,tmn,tmx,iret) 
      endif
      if(iret.eq.1)then
        write(ifmt,
     .  '(a,2i5/a,3f8.2/a,f8.2/a,3f8.2/a,f8.2/a,f8.2/a,f8.2/a,3f8.2)')
     .  'ERROR 02082019',idpomr,jl
     .  ,' qqcut,pt2xx,qcmass**2',qqcut,pt2xx,qcmass**2
     .  ,' s2min1         ',s2min1_save
     .  ,' s2min2,pth     ',s2min2_save,pth1_save,pth4_save
     .  ,' s2minPlusPth   ',s2minPlusPth_save
     .  ,' wpiwmi_save    ',wpiwmi_save
     .  ,' sis            ',sis
     .  ,' smn,tmn,tmx    ',smn,tmn,tmx
        stop
      endif

      smin1=smn
      call getQmaxEoL(jl,sis,scamax)   !KW1908

      if(iiiter.gt.1)then
        if(iabs(j2).eq.4.and.q1.lt.qcmass**2) goto 12 !should not happen in a regular run
        if(iabs(j2).eq.5.and.q1.lt.qbmass**2) goto 12
        if(iabs(l2).eq.4.and.q2.lt.qcmass**2) goto 12
        if(iabs(l2).eq.5.and.q2.lt.qbmass**2) goto 12
      endif

      if(iabs(j2).eq.4.and.q1.lt.qcmass**2) stop'ERROR 02072019'
      if(iabs(j2).eq.5.and.q1.lt.qbmass**2) stop'ERROR 02072019'
      if(iabs(l2).eq.4.and.q2.lt.qcmass**2) stop'ERROR 02072019'
      if(iabs(l2).eq.5.and.q2.lt.qbmass**2) stop'ERROR 02072019'

      !-------------------------------------------------
      !   Born generation laddtype=4    both sided emissions
      !-------------------------------------------------

      if(laddtype.eq.4)then ! both sided emissions
        !inspired from Gauss integration in psjet, see discussion before "function psjet"

      call setLadderindex(laddtype)
      ienvi=5
      fmax=0
      fsum=0
      ntry=0
  10  ntry=ntry+1
      if(ntry.eq.ntryBORN)then
        write(ifmt,'(a,$)')'WARNING Many reject ladd4'
        !call checkCS(q1,q2,qqcut,sis,j1,l1 
        !.,sTYP(1,jl)/ee44x,sTYP(2,jl)/ee44x
        !.,sTYP(3,jl)/ee44x,sTYP(4,jl)/ee44x
        !.,sTPex,ihit,0)  
        !if(ihit(1)+ihit(2)+ihit(3)+ihit(4).eq.0)print*,' '
        write(ifmt,'(f7.1,2f15.11,3f12.8,i3,a)')sis,fmax,fsum/ntry
     .  ,psjet(q1,q2,qqcut,sis,j1,l1,0),sTYP(4,jl)/ee44x,qqcut,iqq
     .  ,'SKIP POM'
        iret=1
        imark=10
        goto 16 !skip Pomeron
      endif  
      !---------proposal---------
      !changing variables:
      !dtdx1*dx2=dtdz*dx1/x1  z=x1*x2 We use t,z,x1 (dont forget the 1/x1 at the end) 
      zzmin=dble((smin1+pth)/sispth)   
      zzmax=(1-epscutSud2(dble(scamax)))*(1-epscutSud2(dble(scamax)))
      zzmx=zzmax**(-delh)
      zzmn=zzmin**(-delh)
      uuu=rangen()*2-1
      zzz=((zzmx+zzmn-uuu*(zzmn-zzmx))/2)**(-1/delh)  !<------
      zPropi=1/delh/2*(zzmn-zzmx)*zzz**(1+delh) !Jacobien = 1 / proposal (Propi)
      sh=zzz*sispth-pth
      call getLimitsEoL(jl,sh,qqcut , smn,tmn,tmx,iret)
      u=rangen()*2-1
      t=2*tmn/(1+tmn/tmx-u*(1-tmn/tmx))  !<------
      tPropi=t**2*0.5*(1/tmn-1/tmx)      !Jacobien = 1 / proposal (Propi)
      xxmax=-1.    !to have xxmin.ge.xxmax if no klas is allowed (CTP : is that correct ??????)
      do klas=1,klasmax 
        call HardScale(klas,qqcut,pt2min,1) !get pt2min 
        call getBornKin(klas,sh,pt2min, sminDmy,tmin1Dmy,tmax1Dmy,iret) !get tmax1
        if(iret.eq.0)then
          call getBornPtr(klas,sh,sngl(t) , pt2) !get pt2
          call HardScale(klas,scale,pt2,2) !get scale
          xxmax=max(xxmax,1.d0-epscutSud2(dble(scale)))
          xxmax=min(xxmax,1.0d0-1.0d-9) !pt2 and scale might ne negative
        endif
      enddo
      xxmin=max(sqrt(zzz),zzz/xxmax)  !we consider lower triangle of x1-x2 plane, x1>x2
      if(xxmin.ge.xxmax)goto 10
      if(xxmax.le..8d0)then
        p1=1
      elseif(xxmin.ge..8d0)then
        p1=0
      else
        p1=0.5
      endif
      v=rangen()
      uuu=rangen()*2-1
      if(v.lt.p1)then
        xxmax1=min(xxmax,.8d0)
        xx1=xxmin*(xxmax1/xxmin)**dble(.5*(1+uuu))  !<------
        x1Propi=xx1*log(xxmax1/xxmin)/2          !Jacobien = 1 / proposal (Propi)
      else
        xxmin1=max(xxmin,.8d0)
        xx1=1-(1-xxmax)*((1-xxmin1)/(1-xxmax))**((1-uuu)/2) !<------
        x1Propi=xx1*(1/xx1-1)/2*log((1-xxmin1)/(1-xxmax)) !Jacobien = 1 / proposal (Propi)
      endif 
      if(.not.(xx1.gt.0..or.xx1.le.0.))then
        print*,zzz,xx1,v.lt.p1, zzz,xxmin,xxmax
        stop'ERROR xx1=NaN'
      endif
      allPropi=zPropi*tPropi*x1Propi
      xx2=zzz/xx1
      if(rangen()<0.5)then
        xx1x=xx2
        xx2x=xx1
      else
        xx1x=xx1
        xx2x=xx2
      endif
      !---------rejection procedure (isi=1)------------ 
      !-----or final detailed calculation (isi=-1)-----
      isi=1  
      jx=j1 !parton pair
      lx=l1 !  group
      ienvi=6
  11  f=0      
      do klas=1,klasmax 
        call HardScale(klas,qqcut,pt2min,1)  !get pt2min 
        call getBornKin(klas,sh,pt2min, sminDmy,tmin1,tmax1,iret) !get tmax1
        if(iret.eq.0)then
         call getBornPtr(klas,sh,sngl(t) , pt2) !compute pt2
         if(pt2.ge.pt2min)then
          call HardScale(klas,scale,pt2,2)  !compute scale
          !if(klas.eq.4.and.nrevt.eq.1.and.ntry.eq.2)ish=9
          a=psjeti(ipm,isi*klas,dble(tmax1),q1,q2,scale,t   !double tmax,t
     .       ,xx1x,xx2x,dble(sh),jx,lx,0)  !double x1,x2,s
          !~~~~~~~~~~~~~~~~~TEST~~~~~~~~~~~~~ to check if psjeti = ffsig / zzz  (OK)
          !if(ioTestFact.ge.6)then !makes no sense for < 6
          !b=ffsigSNGL(sis,klas,pt2,xx1x,xx2x,0) / zzz
          !idif=999
          !if(b.gt.0)idif=nint((a-b)/b*100.)
          !if(a.ne.0..or.b.ne.0.)then
          !write(ifmt,'(a,3i3,3e11.3,$)')'BTH-EMI',j2,l2,klas,t,xx1x,xx2x
          !write(ifmt,'(2e12.3,3i5)')a,b,idif !,nrevt,ntry
          !endif 
          !if(ish.ge.9)stop'BTH-EMI'
          !endif
          !~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
          a = a * (2*pi*pssalf(scale/qcdlam))**2 
          f = f + a      
         else
          iret=1
         endif
        endif
        if(iret.ne.0.and.isi.lt.0.and.klas.eq.klasmax)then !call psjeti needed even for iret.ne.0
          dmy=psjeti(ipm,-9999,dble(tmax1),q1,q2,scale,t   
     .       ,xx1x,xx2x,dble(sh),jx,lx,0)  
        endif
        call getKlasString(klas,strg)
        !write(ifmt,'(5f7.2,f11.0,2i3,e12.3,2x,a)')
        !.  scale,t,pt2,xx1,xx2,sh,j1,l1,a,strg
      enddo
      f = f * pi / sh**2  
      f = f * allPropi / xx1
      f = f *fscale4 !arbitrary scale for rejection
      if(isi.eq.1)then
        nn=0
        if(ncount10.gt.0)nn=ncount11/ncount10
        if(f.gt.1)write(ifmt,'(a,3f8.2,i8)')'WARNING rejection ratio 4'
     .  ,f,sqrt(t),sqrt(sis),nn
        fmax=max(fmax,f)
        fsum=fsum+f
        fx=f
        if(rangen().gt.f)goto 10 !reject
        ncount10=ncount10+1
        ncount11=ncount11+ntry
        if(ish.ge.4)then!~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
          write(ifmt,'(a,2i3,2f11.4,i8)')'accepted',j2,l2
     .    ,f,sqrt(t),ntry  !,sTYP(4,jl)  
        endif!~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
        jx=j2 !specific
        lx=l2 !parton pair
        call setBornInWeightsZero !initialize probability table for Born-in flavors 
        isi=-1  !this activates computation and storage of Born-in flavors in psjeti
        goto 11 
      endif
      call getBornIn(mIB,nIB) !get stored Born-in flavors 
      call getBornOutRdm(mOB,nOB)!Born-out (generated based on /cweightb/)
      klas=igetKlas( mIB,nIB,  mOB,nOB )
      call getBornPtr(klas,sh,sngl(t) , pt2) !compute pt2
      sBorn=sh
      tBorn=t
      pt2Born=pt2 
      xpBorn=xx1x
      xmBorn=xx2x
      !y calculation approximate, ignoring pth
      a=sqrt(pt2/sis)
      det=1.-4.*a**2/zzz
      if(det.lt.-0.000001)then
      write(ifmt,*)"psahot: det<0 (1) ",det,a,zzz,xx1x,xx2x,pt2,sis,sh,t
        call utstop("Negative determinant in psahot (1) !&")
      elseif(det.le.0.)then    !only precision problem
        det=0.
      endif
      if(rangen().lt.0.5)then
        yBorn=-log( 2*a/xx1x / ( 1-sqrt( det ) ) )
      else
        yBorn=-log( 2*a/xx1x / ( 1+sqrt( det ) ) )
      endif
      vparam(84)=vparam(84)+1 
      if(mIB.eq.0.and.nIB.eq.0)then
        vparam(85)=vparam(85)+1 
        if(abs(mOB).eq.4.and.abs(nOB).eq.4)vparam(86)=vparam(86)+1
      endif    
      if(mIB*nIB*mOB*nOB.lt.0.or.mIB*nIB.gt.0.and.mOB*nOB.eq.0)then
        write(ifmt,*)'ERROR invalid flavors (4)',mIB,nIB,mOB,nOB
        stop
      endif

      !----------------------------------------------------------
      ! Born generation laddtype=2,3  one sided emissions (2=pr,3=tg)
      !----------------------------------------------------------

      elseif(laddtype.eq.2.or.laddtype.eq.3)then ! one sided emissions

      call setLadderindex(laddtype)
      ienvi=7
      fmax=0
      fsum=0
      ntry=0
  20  ntry=ntry+1
      if(ntry.eq.ntryBORN)then
        write(ifmt,'(a,$)')'WARNING Many reject ladd23'
        !call checkCS(q1,q2,qqcut,sis,j1,l1 
        !.,sTYP(1,jl)/ee44x,sTYP(2,jl)/ee44x,sTYP(3,jl)/ee44x
        !.,sTYP(4,jl)/ee44x,sTPex,ihit,0)  
        !if(ihit(1)+ihit(2)+ihit(3)+ihit(4).eq.0)print*,' '
        if(laddtype.eq.2)sexa=psjet1(q1,q2,qqcut,sis,j1,l1,0)
        if(laddtype.eq.3)sexa=psjet1(q2,q1,qqcut,sis,l1,j1,0)
        write(ifmt,'(f7.1,f9.5,f11.7,2f10.6,a)')sis,fmax,fsum/ntry
     .  ,sexa,sTYP(laddtype,jl)/ee44x,' SKIP POM'
        ! if(laddtype.eq.3)write(ifmt,'(a,4f10.3,3i4,2f10.3)')
        !.'EXACT.NE.IPOL: q1,q2,qqcut,sis,j,l,jl,psjet1,sTYP(3,)  ='
        !.,q2,q1,qqcut,sis,l1,j1,jl,sexa,sTYP(laddtype,jl)/ee44x
        iret=1
        imark=11
        goto 16 !skip Pomeron
      endif  
      xmin=dble((smin1+pth)/sispth)
      xmax=1-epscutSud2(dble(scamax))     
      if(xmax.le..8d0)then
        p1=1
      elseif(xmin.ge..8d0)then
        p1=0
      else
        p1=0.5
      endif
      v=rangen()
      if(v.lt.p1)then
        zmx=min(xmax,.8d0)**(-delh)
        zmn=xmin**(-delh)
        u=rangen()*2-1
        z=((zmx+zmn-u*(zmn-zmx))/2)**(-1/delh)   !<------
        zPropi=1/delh/2*(zmn-zmx)*z**(1+delh)        !<------
        sh=z*sispth-pth
        call getLimitsEoL(jl,sh,qqcut , smn,tmn,tmx,iret)
        u=rangen()*2-1
        t=2*tmn/(1+tmn/tmx-u*(1-tmn/tmx))  !<------
        tPropi=t**2*0.5*(1/tmn-1/tmx)         !<------
      else
        zmn=max(xmin,.8d0)
        zmx=xmax
        u=rangen()*2-1
        zu=u                                   !<------
        z=1-(1-zmx)*((1-zmn)/(1-zmx))**((1-u)/2) !<-----
        zPropi=(1-z)/2*log((1-zmn)/(1-zmx))         !<------
        sh=z*sispth-pth
        call getLimitsEoL(jl,sh,qqcut , smn,tmn,tmx,iret)
        u=rangen()*2-1
        tu=u                             !<------
        t=2*tmn/(1+tmn/tmx-u*(1-tmn/tmx))  !<------
        tPropi=t**2*0.5*(1/tmn-1/tmx)         !<------
      endif
      allPropi=zPropi*tPropi
      fextra=1
      sh=z*sispth-pth
      !--------------rejection procedure (isi=1)-------------
      !---------or final detailed calculation (isi=-1)-------
      isi=1  
      jx=j1 !parton pair
      lx=l1 !  group
  21  f=0      
      do klas=1,klasmax
        checkWeightSums(1,isi*klas)=0
        checkWeightSums(2,isi*klas)=0
        call HardScale(klas,qqcut,pt2min,1)  !get pt2min 
        pt2min=max(0.,pt2min)
        call getBornKin(klas,sh,pt2min, sminDmy,tmin1,tmax1,iret) !get tmax1
        if(iret.eq.0)then
         !if(rangen().lt.0.5)t=2*tmax1-t
         call getBornPtr(klas,sh,sngl(t) , pt2) !compute pt2
         if(pt2.ge.pt2min)then
          call HardScale(klas,scale,pt2,2)  !compute scale
          if(laddtype.eq.2)then
            jxx=jx
            lxx=lx
            q1xx=q1
            q2xx=q2
          else !switch sides
            iSwitchSides=1
            jxx=lx
            lxx=jx
            q1xx=q2
            q2xx=q1
          endif
          !if(klas.eq.1.and.nrevt.eq.1.and.ntry.eq.2)ish=9 
          a=psjetj(ipm,isi*klas,dble(tmax1),q1xx,q2xx,scale,t  !double tmax,t
     .       ,dble(z),dble(sh),jxx,lxx)           !double z,sh
          a=a*psuds(scale,lxx)/psuds(q2,lxx)
          a= a * fextra
          !~~~~~~~~~~~~~~~~~TEST~~~~~~~~~~~~~ to check if psjetj = ffsig / z  (OK)
          !if(ioTestFact.ge.6)then !makes no sense for < 6
          !b=ffsigSNGL(sis,klas,pt2,dble(z),1.d0,0) / z
          !idif=999
          !if(b.gt.0)idif=nint((a-b)/b*100.)
          !if(a.ne.0..or.b.ne.0.)then
          !write(ifmt,'(a,3i3,4e12.3,$)')'ONE-EMI',j2,l2,klas,t,z,a,b
          !write(ifmt,'(i8,8x,2i8)')idif !,nrevt,ntry
          !endif 
          !if(ish.ge.9)stop'ONE-EMI' 
          !endif
          !~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
          a=a*(2*pi*pssalf(scale/qcdlam))**2
          f = f + a 
         else
          iret=1
         endif
        endif
        if(iret.ne.0.and.isi.lt.0.and.klas.eq.klasmax)then !call psjetj needed even for iret.ne.0
          dmy=psjetj(ipm,-9999,dble(tmax1),q1xx,q2xx,scale,t  
     .       ,dble(z),dble(sh),jxx,lxx)           
        endif
        call getKlasString(klas,strg)
        !write(ifmt,'(4f7.2,f11.0,2i3,e12.3,2x,a)')
        !.  scale,t,pt2,z,sh,j1,l1,a,strg
      enddo 
      f = f * pi / sh**2  
      !~~~~~~~~~~~~~~~~~TEST~~~~~~~~~~~~~use only klas=1 
      !fori=f
      !fbbb= ffsigSNGL(sis,1,pt2,dble(z),1.d0,0) / z
      !.*(2*pi*pssalf(scale/qcdlam))**2 * pi / sh**2 
      !f= 4  / z   * 0.5*(1-(z-0.01)*40)/pt2**2.8 !use same parameterization in ffsigii   
      !print*,'TEST fparam'  ,z,sqrt(pt2),fori,fbbb,f
      !~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
      f = f * allPropi
      f = f *fscale23 !arbitrary scale for rejection
      if(isi.eq.1)then
        nn=0
        if(ncount10.gt.0)nn=ncount11/ncount10
        if(f.gt.1)write(ifmt,'(a,2f8.2,i8)')'WARNING rejection ratio 23'
     .  ,f,sqrt(t),nn
        fmax=max(fmax,f)
        fsum=fsum+f
        fx=f
        if(rangen().gt.f)goto 20 !reject
        ncount10=ncount10+1
        ncount11=ncount11+ntry
        if(ish.ge.4)then!~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
          write(ifmt,'(a,2i3,2f11.4,i8)')'accepted',j2,l2
     .    , f,sqrt(t),ntry  !,sTYP(2,jl)+sTYP(3,jl)  
        endif!~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
        jx=j2 !specific
        lx=l2 !parton pair
        call setBornInWeightsZero !initialize probability table for Born-in flavors 
        isi=-1  !this activates computation and storage of Born-in flavors in psjetj
        goto 21 
      endif
      if(laddtype.eq.2)then
        call getBornIn(mIB,nIB) !get stored Born-in flavors 
        call getBornOutRdm(mOB,nOB)!Born-out (generated based on /cweightb/)
      else !switch sides back
        call getBornIn(nIB,mIB) !get stored Born-in flavors 
        call getBornOutRdm(nOB,mOB)!Born-out (generated based on /cweightb/)
        iSwitchSides=0
      endif
      klas=igetKlas( mIB,nIB,  mOB,nOB )
      call getBornPtr(klas,sh,sngl(t) , pt2) !compute pt2
      sBorn=sh
      tBorn=t
      pt2Born=pt2 
      !y calculation approximate, pth ignored
      a=sqrt(pt2/sis)
      det=1.-4.*a**2/z
      if(det.lt.-0.000001)then
        write(ifmt,*)"psahot: det<0 (2) ",det,a,z,pt2,sis,sh,t
        call utstop("Negative determinant in psahot (2) !&")
      elseif(det.le.0.)then    !only precision problem
        det=0.
      endif
      vparam(87)=vparam(87)+1 
      if(mIB.eq.0.and.nIB.eq.0)then
        vparam(88)=vparam(88)+1 
        if(abs(mOB).eq.4.and.abs(nOB).eq.4)vparam(89)=vparam(89)+1
      endif    
      if(laddtype.eq.2)then
        xpBorn=z
        xmBorn=1.0d0 
        yBorn=log( 2*a/( 1-sqrt( det ) ) ) !bigger solution 
        !yBorn=log( 2*a/( 1-sqrt( 1+4*a**2/z ) ) ) 
        if(nIB.ne.nEOL)then
          write(ifmt,*)'EOL to IB',mEOL,nEOL,'-->',mIB,nIB
          stop'ERROR invalid flavor change (2)'
        endif
      else
        xpBorn=1.0d0
        xmBorn=z 
        yBorn=-log( 2*a/( 1-sqrt( det ) ) ) !smaller solution 
        !yBorn=-log( 2*a/( 1+sqrt( 1-4*a**2/z ) ) )  
        if(mIB.ne.mEOL)then
          write(ifmt,*)'EOL to IB',mEOL,nEOL,'-->',mIB,nIB
          stop'ERROR invalid flavor change (3)'
        endif
      endif
      !check: inverse formua:
      ! zp= a*exp(-yBorn) / ( 1-a*exp(yBorn) )  !should be = z 
      if(mIB*nIB*mOB*nOB.lt.0.or.mIB*nIB.gt.0.and.mOB*nOB.eq.0)then
        write(ifmt,*)'ERROR invalid flavors (2)',mIB,nIB,mOB,nOB
        stop
      endif

      !-------------------------------------------------
      !   Born generation laddtype=1    no emissions
      !-------------------------------------------------

      elseif(laddtype.eq.1)then ! no emissions
        !inspired from Gauss integration in psjetk, see discussion before "function psjet"

      ienvi=9
      fmax=0
      fsum=0
      tmx_save=tmx
      fscale1sat=1
      if(iqq.lt.0)then
        call setLadderindex(5)
        jdis=-1
        tmx=tmn+0.2
        scale=qqsat
        fscale1sat=50  !<===== depends strongly on "fine tuning sat"
        if(tmn.lt.4.00)fscale1sat=fscale1sat*0.5
        if(tmn.lt.3.00)fscale1sat=fscale1sat*0.5
        if(tmn.lt.2.00)fscale1sat=fscale1sat*0.5
        if(tmn.lt.1.50)fscale1sat=fscale1sat*0.1
        if(tmn.lt.1.00)fscale1sat=fscale1sat*0.3
        if(tmn.lt.0.60)fscale1sat=fscale1sat*0.3
        if(tmn.lt.0.30)fscale1sat=fscale1sat*0.1
        if(tmn.lt.0.10)fscale1sat=fscale1sat*0.2
        if(engy.lt.1000)then
        if(tmn.lt.2.00)fscale1sat=fscale1sat*0.2
        if(tmn.lt.1.50)fscale1sat=fscale1sat*5.0
        endif
        !+++++++++++++++++++++++++
        !cKW2108 
        !To be consistent with psborn (no integration) we should use t=tmn.
        !To keep the normal procedure, we generate t bewteen tmn and tmn+dt 
        !+++++++++++++++++++++++++
      else
        call setLadderindex(1)
        jdis=0
      endif
      ntry=0
  30  ntry=ntry+1
      if(ntry.eq.ntryBORN)then
        if(iqq.ge.0)then
          write(ifmt,'(a,$)')'WARNING Many reject ladd1; SKIP POM '
          call checkCS(q1,q2,qqcut,sis,j1,l1 
     .    ,sTYP(1,jl)/ee44x,sTYP(2,jl)/ee44x,sTYP(3,jl)/ee44x
     .    ,sTYP(4,jl)/ee44x,sTPex,ihit,0)  
          if(ihit(1)+ihit(2)+ihit(3)+ihit(4).eq.0)print*,' '
          iret=1
          imark=12
          goto 16 !skip Pomeron
        else !sat pom
          write(ifmt,'(a,f8.2)')'WARNING Many reject PomSat; tmn =',tmn 
        endif
      endif  
      u=rangen()*2-1
      tu=u                             !<------
      t=2*tmn/(1+tmn/tmx-u*(1-tmn/tmx))  !<------
      allPropi=t**2*0.5*(1/tmn-1/tmx)         !<------
      z=1.
      sh=sis
      !-------------rejection procedure (isi=1)-------------
      !------- or final detailed calculation (isi=-1)-------
      isi=1  
      jx=j1 !parton pair
      lx=l1 !  group
  31  f=0   
      do klas=1,klasmax
        if(iqq.ge.0)then ! not PomSat
          call HardScale(klas,qqcut,pt2min,1)  !get pt2min 
          call getBornKin(klas,sh,pt2min, sminDmy,tmni,tmxi,iret) !get tmxi
        else
          call getBornPtr(klas,sh,sngl(t) , pt2) !compute pt2
          !if(isi*klas.eq.-4)print*,'TEST',sh,sngl(t) , pt2
          pt2min=0.001
          iret=1
          if(pt2.ge.pt2min)then
            call getBornKin(klas,sh,pt2, sminDmy,tmni,tmxi,iret) !get tmxi and iret
          endif
          !if(isi*klas.eq.-4)print*,'TEST',sh,pt2, sminDmy,tmni,tmxi
        endif 
        if(iret.eq.0)then
         !if(rangen().lt.0.5)t=2*tmax1-t
         if(iqq.ge.0)call getBornPtr(klas,sh,sngl(t) , pt2) !compute pt2
         if(pt2.ge.pt2min)then 
          if(iqq.ge.0)call HardScale(klas,scale,pt2,2)  !compute scale
          !if(klas.eq.1.and.nrevt.eq.1.and.ntry.eq.2)ish=9 
          a=psjetk(ipm,isi*klas,dble(tmxi),q1,q2,scale
     .      ,t,dble(sh),jx,lx,jdis)  !double tmax,t,s
          if(iqq.ge.0)then ! not PomSat
            a=a*psuds(scale,j1)*psuds(scale,l1)
            a=a/psuds(q1,j1)/psuds(q2,l1) 
          endif
          !~~~~~~~~~~~~~~~~~TEST~~~~~~~~~~~~~ to check if psjetk = ffsig   (OK)
          !if(ioTestFact.ge.6)then !makes no sense for < 6
          !b=ffsigSNGL(sis,klas,pt2,1.d0,1.d0,0)
          !idif=999
          !if(b.gt.0)idif=nint((a-b)/b*100.)
          !if(a.ne.0..or.b.ne.0.)then
          !write(ifmt,'(a,3i3,2e12.3,i8)')'NO-EMI',j2,l2,klas,a,b,idif
          !endif
          !if(ish.ge.9)stop'NO-EMI'
          !endif
          !~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
          a=a*(2*pi*pssalf(scale/qcdlam))**2
          f = f + a 
         else
          iret=1
         endif
        endif
        if(iret.ne.0.and.isi.lt.0.and.klas.eq.klasmax)then !call psjetk needed even for iret.ne.0
          dmy=psjetk(ipm,-9999,dble(tmxi),q1,q2,scale,t,dble(sh),jx,lx
     .              ,jdis) 
        endif
        call getKlasString(klas,strg)
        !write(ifmt,'(4f7.2,f11.0,2i3,e12.3,2x,a)')
        !.  scale,t,pt2,z,sh,j1,l1,a,strg
      enddo 
      f = f * pi / sh**2  
      f = f * allPropi
      f = f *fscale1*fscale1sat !arbitrary scale for rejection
      if(isi.eq.1)then
        if(f.gt.1)write(ifmt,'(a,i3,4f8.2)')'WARNING rejection ratio 1:'
     .  ,iqq,f,qqsat,tmn,pt2xx
        fmax=max(fmax,f)
        fsum=fsum+f
        fx=f
        if(rangen().gt.f)goto 30 !reject
        if(ish.ge.4)then!~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
          write(ifmt,'(a,2i3,2f11.4,i8)')'accepted',j2,l2
     .    , f,sqrt(t),ntry  !,sTYP(1,jl)  
        endif!~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~  
        jx=j2 !specific
        lx=l2 !parton pair
        call setBornInWeightsZero !initialize probability table for Born-in flavors 
        isi=-1   !this activates computation and storage of Born-in flavors in psjetk (selection of input quark flavor AND klas of born process) 
        goto 31 
      endif
      call getBornIn(mIB,nIB) !get stored Born-in flavors 
      call getBornOutRdm(mOB,nOB)!Born-out (generated based on /cweightb/)
      klas=igetKlas( mIB,nIB,  mOB,nOB )
      call getBornPtr(klas,sh,sngl(t) , pt2) !compute pt2
      xt=sqrt(pt2/(sis/4)) !ladder cms
      sBorn=sh
      tBorn=t
      pt2Born=pt2 
      if(iqq.lt.0)then
        vparam(81)=vparam(81)+1 
        if(mIB.eq.0.and.nIB.eq.0)then
          !print*,'TEST-sat-kl4',jsplit,pt2,qqcut,sngl(t),sh
          !.  ,nint(getWeighta(4,0,0)/getWeighta(1,0,0)*1000), mOB,nOB 
          vparam(82)=vparam(82)+1 
          if(abs(mOB).eq.4.and.abs(nOB).eq.4)then
            vparam(83)=vparam(83)+1
            ntxt80=ntxt80+1
            !write(txt80(ntxt80),*)
            !.'TEST-sat44',npr(3,kcol),igetNpom(),estpom
          endif 
        endif    
        !write(*,'(a,5f10.1)')'TEST-psahot',tmn,tmx,t,sh,sqrt(pt2)
        !REMOVE commented ifav redo after 4.0.2.r4
      else
        vparam(90)=vparam(90)+1 
        if(mIB.eq.0.and.nIB.eq.0)then
          vparam(91)=vparam(91)+1 
          if(abs(mOB).eq.4.and.abs(nOB).eq.4)vparam(92)=vparam(92)+1
        endif    
      endif 
      a=sqrt(pt2/sis)
      det=1.-4.*a**2
      if(det.lt.-0.000001)then
        write(ifmt,*)"psahot: det<0 (3) ",det,a,pt2,sis,sh,t
        call utstop("Negative determinant in psahot (3) !&")
      elseif(det.le.0.)then    !only precision problem
        det=0.
      endif
      yBorn=-log( 2*a/( 1-sqrt( det ) ) )
      !yp=acosh(1./xt) !from s+t+u=0 (see pub/19epos4)  gives same result
      xpBorn=1.0d0
      xmBorn=1.0d0
      if(mIB.ne.mEOL)then
        write(ifmt,*)'ERROR invalid flavors (1a)',mEOL,nEOL,mIB,nIB
        stop
      endif
      if(nIB.ne.nEOL)then
        write(ifmt,*)'ERROR invalid flavors (1b)',mEOL,nEOL,mIB,nIB
        stop
      endif
      if(mIB*nIB*mOB*nOB.lt.0.or.mIB*nIB.gt.0.and.mOB*nOB.eq.0)then
        write(ifmt,*)'ERROR invalid flavors (1c)',mIB,nIB,mOB,nOB
        stop
      endif
      tmx=tmx_save

      endif

      if(ish.ge.3)then
        write(ifmt,'(a,3(2i3,2x),2x,a,i4,4x,a,i3,4x,a,f7.4,$)')
     .  'HARD  Flav=',j2,l2,mIB,nIB,mOB,nOB,'Ntry=',ntry
     .  ,'LaddType=',laddtype,'Fmax=',fmax   
        write(ifmt,'(4x,a,f9.1,3x,a,4f9.3)')'s=',sBorn
        !.,'ykk=',ykk(1,igroup),ykk(2,igroup),ykk(3,igroup)
        !.,1-ykk(1,igroup)-ykk(2,igroup)-ykk(3,igroup)
      endif

      !-------------------------------------------------------
      !at this point, we know 
      !     sBorn, pt2Born, tBorn, xpBorn, xmBorn, yBorn
      ! and the flavors mEOL,nEOL,  mIB,nIB,  mOB,nOB 
      ! and the ladder momentum pLadder (4 vector SO convention)
      !------------------------------------------------------- 

      call timer(iutime)
      tidi=iutime(3)-tiu3 + (iutime(4)-tiu4)*1d-3
      if(ltime)write(ifmt,'(a,i15)')'cpu time psahot 3/5'
     .,nint(tidi*1000)
      tiu3=iutime(3)
      tiu4=iutime(4)

      ienvi=10
      iret=0

      !KKKKKKKKKKKKKKKKKKKKKKKKKKKKKKKKKKKKKKKKKKKKKKKKKKKKKKKKKKKK
      !KK                                                        KK
      !KK             TEST9  TEST10 TEST11                       KK
      !KK                                                        KK
      !KKKKKKKKKKKKKKKKKKKKKKKKKKKKKKKKKKKKKKKKKKKKKKKKKKKKKKKKKKKK

      if(ioTestFact.ge.9.and.ioTestFact.le.11)then
        pt=sqrt(pt2Born)        
        yPairCMS=0.5*log(xpBorn/xmBorn)
        yStar=yBorn-yPairCMS
        yBorn2=yPairCMS-yStar
        rap1=yLadder+yBorn   !ignore rotations for the moment
        rap2=yLadder+yBorn2
        ni=kdimtt
        !~~~~~~~~~~~~~~~~~~~~~~
        ! Check ffsigii
        !~~~~~~~~~~~~~~~~~~~~~~
        !useful in case of prints in ffsigii
        !dmy1=ffsigi(2.25 , 1000. , 0 , 1 )
        !ish=11
        !dmy2=ffsigi(2.25 , 1006. , 0 , 1 )
        !stop'TEST9'
        !~~~~~~~~~~~~~~~~~~~~~~~~~~~~
        ! Check rapidity distribution
        !~~~~~~~~~~~~~~~~~~~~~~~~~~~~
        ! Case of small pt: Broad U shape
        ! Getting narrower and less 
        !  pronounced with increasing pt 
        !~~~~~~~~~~~~~~~~~~~~~~~~~~~~
        if(iopcnt.le.0)then !no batch
        irap1= rap1 + ni/2+0.5 +1
        do j=1,jdimtt
         if(pt.gt.j+0.0.and.pt.le.j+1.0)then
          if(irap1.ge.1.and.irap1.le.ni)vartt(j,irap1)=vartt(j,irap1)+1
         endif   
        enddo
        if(nrevt.ge.nevent-10)then ! +99 -10
         ncount3=ncount3+1
         if(ncount3.eq.1)then   
          varttsum=0
          do j=1,jdimtt !pt-loop
           vartt(j,0)=0
           do i=1,ni
            vartt(j,0)=vartt(j,0)+vartt(j,i)
           enddo
           varttsum=varttsum+vartt(j,0)
          enddo
          q2cmin(1)=q2nmin
          q2cmin(2)=q2nmin
          warttsum=0
          do j=1,jdimtt !pt-loop
           pt2j=(0.5+j)**2 
           wartt(j,0)=0
           do i=1,ni !y-loop
            y0=sign(1000,i-ni/2-1)+(i-ni/2-1) 
            wartt(j,i)=ffsigi(pt2j , y0 , 0 , 1 )
            wartt(j,0)=wartt(j,0)+wartt(j,i)
           enddo
           warttsum=warttsum+wartt(j,0)
          enddo
          write(ifmt,'(3a)')'----y-----dn/dydpt(pt=1.5)'
     .    ,'---dn/dydpt(pt=2.5)','---dn/dydpt(pt=3.5)---'
          write(ifmt,'(8x,3(4x,a,4x,a,3x))')'Simu','Theo'
     .    ,'Simu','Theo','Simu','Theo'
          do i=1,ni
            iy=i-ni/2-1
            write(ifmt,'(i5,3x,2f8.3,3x,2f8.3,3x,2f8.3)') iy
     .      ,vartt(1,i)/vartt(1,0),wartt(1,i)/wartt(1,0)
     .      ,vartt(2,i)/vartt(2,0),wartt(2,i)/wartt(2,0)
     .      ,vartt(3,i)/vartt(3,0),wartt(3,i)/wartt(3,0)
          enddo
          write(ifmt,'(2a)')'----pt-----dn/dptSimu--dn/dptTheo'
     .                     ,'--(y integrated)----'
          do j=1,jdimtt
            write(ifmt,'(f6.2,2f12.4)') 0.5+j
     .      ,vartt(j,0)/varttsum,wartt(j,0)/warttsum
          enddo
         endif
        endif
        endif!end no batch
        !~~~~~~~~~~~~~~~~~~~~ 
        call ptprboo_set(2,jsplit,ncolp,kcol, 0. )
        do ii=1 ,1 !,2
          if(ii.eq.1)rap=rap1
          if(ii.eq.2)rap=rap2
          call ptprboo_set(ii,jsplit,ncolp,kcol, pt )
          call rapprboo_set(ii,jsplit,ncolp,kcol, rap )
        enddo
c        print *,'here',rap1,rap2,yBorn,yBorn2,Yladder,yPairCMS
        iret=55 !special code, for analysing pt spectra directly
        imark=13
        goto 16 
      endif

      !KKKKKKKKKKKKKKKKKKKKKKKKKKKKKKKKKKKKKKKKKKKKKKKKKKKKKKKKKKKK
      !KK                                                        KK
      !KK                        TEST8                           KK
      !KK                                                        KK
      !KKKKKKKKKKKKKKKKKKKKKKKKKKKKKKKKKKKKKKKKKKKKKKKKKKKKKKKKKKKK

      if(ioTestFact.eq.8)then
        stop'ioTestFact=8 needs update; pt not defined ...'
        !computing pt and rap according to procedure in p s a b o r   (in LAB frame)
        ept(1)=(wplus+wminus)/2
        ept(2)=(wplus-wminus)/2  !ept is 4 momentum of IB parton pair
        ept(3)=pLadder(3)
        ept(4)=pLadder(4)
        sbo= ept(1)**2-ept(2)**2-ept(3)**2-ept(4)**2
        !if(sbo.gt.100)then
        !  write(ifmt,'(a,2f10.2,3x,4f10.2,$)')'TEST',sbo,sBorn,ept
        !endif
        call psdeftr(dble(sbo),ept,eytst)   !for boost to CMS, defines eytst 
        !if(sbo.gt.100)write(ifmt,'(3x,3f10.3)')eytst
        call pscs(bcos,bsin)          !cos and sin of the polar angle
        call getMassSquared( mOB , q2m1 ) 
        call getMassSquared( nOB , q2m2 ) 
        z=psutz(dble(sbo),dble(pt2+q2m1),dble(pt2+q2m2)) 
        wp3=z*sqrt(sbo)
        wm3=(pt2+q2m1)/wp3 
        if(.not.(wp3.le.0..or.wp3.ge.0.)     !NaN catch
     .    .or.wp3.gt.1e35 .or. wp3.lt.-1e35    !Infinity catch
     .    .or..not.(wm3.le.0..or.wm3.ge.0.)     !NaN catch
     .    .or.wm3.gt.1e35 .or. wm3.lt.-1e35 )then !Infinity catch
          write(ifmt,'(a,$)')'ERROR psahot TEST NaN / Inf catch ; '
          write(ifmt,*)'wp3 wm3 z sbo =',wp3,wm3,z,sbo,sBorn
     .    ,wplus,wminus,ept
          ep3(1)=pt
          ep3(2)=pt 
        else
          ep3(1)=.5*sngl(wp3+wm3)        !4-momentum for out born parton 1
          ep3(2)=.5*sngl(wp3-wm3)        !   in Born CMS   
        endif
        ep3(3)=pt*bcos
        ep3(4)=pt*bsin
        call psdefrot(ep3,s0xh,c0xh,s0h,c0h)   !rotation to z-axis
        if(q2m1.gt.0.5d0.and.abs(mOB).ne.abs(nOB)) then
          en1=(sbo+q2m1)/(2.d0*sqrt(sbo))  
          en2=(sbo-q2m1)/(2.d0*sqrt(sbo))
        else
          en1=sqrt(sbo)/2
          en2=en1
        endif
        do ii=1,2
          if(ii.eq.1)rap=rap1
          if(ii.eq.2)rap=rap2
          !en1, en2 transferred to tim (pprtx) then to psreti 
          !to do rotation, boost, which amounts to:
          if(ii.eq.1)en=en1
          if(ii.eq.2)en=en2
          vSO(1)=en
          vSO(2)=en*(3-2*ii)  !ignore mt
          vSO(3)=0
          vSO(4)=0
          !here we do not consider pth
          call psrotat(vSO,s0xh,c0xh,s0h,c0h)
          call pstrans(vSO,eytst,1)    ! <-------- trafo to LAB  -- makes difference  compared to f_PDF
          ptxx=sqrt(vSO(3)**2+vSO(4)**2) !     for f_PDF we need to assume pt_BornCMS = pt_LAB
          rapxx=sign(1.,vSO(2))*
     *        log((sqrt(vSO(2)**2+ptxx**2)+abs(vSO(2)))/ptxx)
          !if(sis.gt.100.and.pt.gt.6.)then
          !  write(ifmt,'(4x,2i4,4x,2f8.3,$)')mOB,nOB,pt,ptxx
          !  write(ifmt,'(4x,2f8.3,3x,10f10.3)')rap,rapxx
          !  !.,sBorn, tBorn, xpBorn, xmBorn 
          !  !.,pLadder , bcos,bsin
          !endif 
          call ptprboo_set(ii,jsplit,ncolp,kcol, pt )
          call rapprboo_set(ii,jsplit,ncolp,kcol, rap )
        enddo
        iret=55 !special code, for analysing pt spectra directly
        imark=14
        goto 16 
      endif


      !KKKKKKKKKKKKKKKKKKKKKKKKKKKKKKKKKKKKKKKKKKKKKKKKKKKKKKKKKKKK
      !KK                                                        KK
      !KK  backward evolution (BWE)                              KK
      !KK       generate t and z for each iteration              KK
      !KK       z is the momentum fraction of the branching      KK
      !KK                                                        KK
      !KKKKKKKKKKKKKKKKKKKKKKKKKKKKKKKKKKKKKKKKKKKKKKKKKKKKKKKKKKKK

      !ienvi still 10

      klas=igetKlas( mIB,nIB,  mOB,nOB )
      if(klas.eq.0)then
        write(ifmt,'(a,4i4)')'ERROR Wrong klas from igetKlas'
     .  ,mIB,nIB,  mOB,nOB,  laddtype
        stop
      endif
      pt2=pt2Born
      !call getBornPtr(klas,sBorn,tBorn , pt2) !compute pt2
      call HardScale(klas,scale,pt2,2)  !compute scale
      call makeSudx

      do iside=1,2 !proj/targ !------------iside-------------->

      if(iside.eq.1)then
  
        !ladder end flavor and virtuality
        iend=mEOL
        tend=q1  
        xend=1
        !start values for evolution:
        kzero=mIB 
        tzero=scale 
        xzero=xpBorn 
        cside='+++'

      else
  
        !ladder end flavor and virtuality
        iend=nEOL
        tend=q2  
        xend=1
        !start values for evolution:
        kzero=nIB 
        tzero=scale 
        xzero=xmBorn 
        cside='---'

      endif

      kold=99
      told=99.999
      xold=99.999
      nemi=0
      ntry=0 
      n40=0
  40  continue
      n40=n40+1
      if(ish.ge.3.and.n40.gt.1)
     .write(ifmt,'(2a,3i3,2x,3f8.2,2x,3f7.3,i7)')cside
     .,' k t x',iend,kzero,kold,tend,tzero,told,xend,xzero,xold,n40
      if(xzero.gt.0.9999999999d0) goto 45 !no emission, nothing to do
  
      !generate t (named tnew)

      r=rangen()
      a=tend
      b=tzero
      psuds_a = psuds(  a  ,kzero)
      psevj_t = psevj(tend,tzero,xzero,iend,kzero)
      if(psuds_a.eq.0.0.or.psevj_t.eq.0.0)then
        write(ifmt,'(a)')'WARNING rsh FPE skip Pom'
        iret=1
        imark=16
        goto 16 !skip Pomeron
      endif
      fa=psuds(tzero,kzero)
     .  /psuds_a
     .  *psevj(tend,  a  ,xzero,iend,kzero)
     .  /psevj_t  -  r
      fb=1-r
      fc=100 
      nt=0
      do while(nt.lt.100.and.abs(fc).gt.0.001)
        nt=nt+1
        c=(a+b)/2 
        fc=psuds(tzero,kzero)
     .    /psuds(  c  ,kzero)
     .    *psevj(tend,  c  ,xzero,iend,kzero)
     .    /psevj(tend,tzero,xzero,iend,kzero)  -  r
        !write(ifmt,'(a,3f8.2,3x,3f10.6)')'BW find root',a,b,c,fa,fb,fc
        if(fa*fc.le.0.)then
          b=c
          fb=fc
        else
          a=c
          fa=fc
        endif
      enddo
      tnew=c 

      !determine probability W1 for z=x 
      ! (which amounts to momentum fraction unity before the branching)
      
      loback=iend.eq.0
      zmax=1-epscutSud2(dble(tnew)) !upper limit for z
      zmin= xzero/(1-epscutSud2(dble(tnew))) !upper limit for xzero/z from psevj
      zmin=max(zmin,xzero)
      zmed=0.1 
      !print*,'TEST',tnew,xzero,zmax,psevj(tend,tnew,xzero,iend,kzero)
      rat1=0
      rat2=0
      f=0
      if(zmin.gt.zmax)then
        !print*,'PSAHOT',zmax,zmin,xzero
        !stop'back evol: zmin > zmax'   
        goto 43
      endif
      f1min=1e30
      f1max=0 
      f1=0 
      if(zmin.lt.zmed)then
        zmax1=min(zmed,zmax)
        do m=1,2
        do i=1,7
          u=(2*m-3)*xgauss7(i)
          if(loback)then 
            z=zmin*(zmax1/zmin)**((1+u)/2)
            !Jacob=z*log(zmax1/zmin)/2
          else
            z=(zmax1+zmin+(zmax1-zmin)*u)/2
            !Jacob=(zmax1-zmin)/2
          endif
          fia=0
          do ia=-5,5
            fia=fia+psevj(tend,tnew,xzero/z,iend,ia)
     .       *pssalf(tnew/qcdlam)*psfapj(dble(z),ia,kzero) !<-- w/o noflav factor
          enddo 
          if(loback)then
            f1min=min(f1min,fia*log(zmax1/zmin)/2)
            f1max=max(f1max,fia*log(zmax1/zmin)/2)
            !write(ifmt,'(a,2f8.3,3x,f8.3,3x,2f8.3)')
            !.'W1 A',u,z,fia*log(zmax1/zmin)/2,f1min,f1max
          else
            f1min=min(f1min,fia/z*(zmax1-zmin)/2)
            f1max=max(f1max,fia/z*(zmax1-zmin)/2)
            !write(ifmt,'(a,2f8.3,3x,f8.3,3x,2f8.3)')
            !.'W1 B',u,z,fia/z*(zmax1-zmin)/2,f1min,f1max
          endif
          f1=f1+wgauss7(i)*fia
        enddo
        enddo
        if(loback)then
          f1=f1*log(zmax1/zmin)/2 ! /z*z cancels out
        else
          f1=f1/z*(zmax1-zmin)/2
        endif 
      endif
      f2min=1e30
      f2max=0 
      f2=0 
      if(zmax.gt.zmed)then
        zmin1=max(zmed,zmin)
        do m=1,2
        do i=1,7
          u=(2*m-3)*xgauss7(i)
          z=1-(1-zmax)*((1-zmin1)/(1-zmax))**((1-u)/2)
          !Jacob=(1-z)*log((1-zmin1)/(1-zmax))/2
          fia=0
          do ia=-5,5
            fia=fia+psevj(tend,tnew,xzero/z,iend,ia)
     .       *pssalf(tnew/qcdlam)*psfapj(dble(z),ia,kzero) !<-- w/o noflav factor
          enddo 
          fia=fia/z*(1-z)  !(1-z) from Jacob
          f2min=min(f2min,fia*log((1-zmin1)/(1-zmax))/2)
          f2max=max(f2max,fia*log((1-zmin1)/(1-zmax))/2)
          !write(ifmt,'(a,2f8.3,3x,f8.3,3x,2f8.3)')
          !.'W1 C',u,z,fia*log((1-zmin1)/(1-zmax))/2,f2min,f2max
          f2=f2+wgauss7(i)*fia
        enddo
        enddo
        f2=f2*log((1-zmin1)/(1-zmax))/2
      endif
      if(f1min.ne.0.)rat1=f1max/f1min
      if(f2min.ne.0.)rat2=f2max/f2min
      f=f1+f2
  43  fint=f
      g=pssalf(sngl(tnew/qcdlam))*psfapj(xzero,iend,kzero)
     .   * psuds(tnew,iend) / psuds(tend,iend)
      W1=g/(f+g)
      icase=2
      if(rangen().lt.W1)icase=1

      ! Treat the two cases, z=x and z>x

      if(icase.eq.1)then
        xnew=1 !momentum fraction before splitting
        if(ish.ge.5)write(ifmt,'(15x,a,2f8.3)')cside//' xnew/old'
     .  ,xnew,xzero
        !emitted TL parton
        nemi=nemi+1
        call storeTL(iside,nemi,tnew,tzero,xnew,xzero,iend,kzero)
        !next iteration for SL parton(will be the last)
        kold=kzero
        told=tzero
        xold=xzero
        kzero=iend
        tzero=tnew 
        xzero=xnew 
        goto 40
      else !generate z
        fmax=0
        ntry=0
  44    ntry=ntry+1
        if(ntry.eq.ntryBWE)then
          write(ifmt,'(a)')'WARNING too many rejections (BWE); SKIP POM'
          iret=1
          imark=15
          goto 16 !skip Pomeron
        endif  
        if(zmax.le.zmed)then
          p1=1
        elseif(zmin.ge.zmed)then
          p1=0
        else
          p1=0.5
        endif
        v=rangen()
        if(v.lt.p1)then !small x
          zmax1=min(zmed,zmax)
          u=rangen()*2-1
          if(loback)then
            z=zmin*(zmax1/zmin)**((1+u)/2)
            zPropi=z*log(zmax1/zmin)/2
          else
            z=(zmax1+zmin+(zmax1-zmin)*u)/2
            zPropi=(zmax1-zmin)/2
          endif
        else !large x
          zmin1=max(zmed,zmin)
          u=rangen()*2-1
          z=1-(1-zmax)*((1-zmin1)/(1-zmax))**((1-u)/2)
          zPropi=(1-z)*log((1-zmin1)/(1-zmax))/2
        endif
        !rejection procedure 
        f=0
        do ia=-5,5
          f=f+psevj(tend,tnew,xzero/z,iend,ia)
     .     *pssalf(tnew/qcdlam)*psfapj(dble(z),ia,kzero) !<-- w/o noflav factor
        enddo 
        f=f/z*zPropi
        if(fint.ne.0.)f = f / fint
        f = f * fscaleBWE !arbitrary scale for rejection
        if(f.gt.1)write(ifmt,'(a,f8.2)')'WARNING rejection ratio BWE',f
        fmax=max(fmax,f)
        if(rangen().gt.f)goto 44 !reject
        xnew=xzero/z
        if(ish.ge.5)write(ifmt,'(15x,a,2f8.3,3x,a,2f8.3,i8,3x,a,2f8.3)')
     .  cside//' xnew/old'
     .  ,xnew,xzero,'fmax fint ntry',fmax,fint,ntry,'rat1,2',rat1,rat2
        !generate flavor
        sum=0
        do ia=-5,5
          fkk(ia+6)=psevj(tend,tnew,xzero/z,iend,ia)
     .     *psfapj(dble(z),ia,kzero) !<-- w/o noflav factor
          sum=sum+fkk(ia+6)
        enddo 
        !~~~~~~~~~~check~~~~~~~~~
        if(sum.le.0.)then
          write(ifmt,'(a,$)')'ERROR psahot BWE Splitting weights zero ;'
          write(ifmt,*)(fkk(ia+6),ia=-5,5)
          stop
        endif
        !~~~~~~~~~~~~~~~~~~~~~~~~
        ia=igetRNDindex(11,fkk)-6
        !emitted TL parton
        nemi=nemi+1
        call storeTL(iside,nemi,tnew,tzero,xnew,xzero,ia,kzero)
        !next iteration
        kold=kzero
        told=tzero
        xold=xzero
        kzero=ia
        tzero=tnew 
        xzero=xnew 
        goto 40
      endif

  45  continue

      nemis(iside)=nemi

      enddo ! <--------------iside---------------

      enddo !iiiTest 

      !KKKKKKKKKKKKKKKKKKKKKKKKKKKKKKKKKKKKKKKKKKKKKKKKKKKKKKKKKKK
      !KK                                                       KK
      !KK      realize emissions (now forward evolution)        KK
      !KK                                                       KK
      !KKKKKKKKKKKKKKKKKKKKKKKKKKKKKKKKKKKKKKKKKKKKKKKKKKKKKKKKKKK

      ienvi=15

      iEOL(1)=mEOL
      iEOL(2)=nEOL

      do ii=1,2  !both sides
      do n=nemis(ii),1,-1 !forward evolution
      ! ----------------------emissions---------------------->

      !process TL parton, including TL cascade

      call getTL(ii,n,t1,t2,xx1,xx2,i1,i2) !get t,x,       1 |
      qq=t2                                !flavors i        |--emitted
      x=xx2/xx1                            !(1=in/2=out)   2 | 
      if(i1.ne.0.and.i2.ne.0.and.i1.ne.i2)then
        write(ifmt,'(a,$)')'ERROR psahot Emissions Wrong flavors;  '
     .  ,'i1 i2 =',i1,i2
        stop
      endif
      iq1=i1-i2 !flavor of emitted parton 
      !~~~~~~~~~~check~~~~~~~~~
      if(abs(iq1).eq.6)then
        write(ifmt,'(a,$)')'ERROR psahot Evol Top detected ;  '
        write(ifmt,*)'i1 i2 iq1 =',i1,i2,iq1
        stop
      endif
      !~~~~~~~~~~~~~~~~~~~~~~~~
      ncr(2)=0
      if(i1.eq.0.and.i2.eq.0)then ! g-g 
        jt=1
        jq=int(1.5+rangen())
        ncr(1)=ncc(jq,ii)
      elseif(i1.ne.0.and.i2.eq.0)then ! q-g
        jt=2                      
        if(i1.gt.0)then !quark     
          jq=1                    
        else !antiquark           
          jq=2                    
        endif                     
        ncr(1)=0
      elseif(i1.ne.0.and.i2.ne.0)then ! q-q
        jt=3
        if(i1.gt.0)then
          jq=1
        else
          jq=2
        endif
        ncr(1)=ncc(1,ii)
      elseif(i1.eq.0.and.i2.ne.0)then ! g-q
        jt=4
        if(i2.gt.0)then !SL quark
          jq=2  
        elseif(i2.lt.0)then !SL antiquark
          jq=1 
        endif
        ncr(1)=ncc(jq,ii)
      endif
      am2=0
      if(abs(iq1).eq.4)am2=qcmass**2
      if(abs(iq1).eq.5)am2=qbmass**2
      !~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
      !be k the momentum of the new SL parton and q
      !the momentum of the emitted parton. Thesis BG 
      !eq.IV.15 : q^2 = (1-x)*|k^2|/x - kt^2/x. We have 
      !as well qt^2=kt^2. Considering q^2=m^2 we get
      !qt^2 = Q^2 * (1-x) - x * m^2, so this would give
      !qt2=qq*(1-x)-x*am2, which is however incompatible 
      !with DGLAP (possibly negative for perfectly valid x and qq) 
      !So we use max(0.,...)
      !~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ 
      qt2=max(0.,qq*(1-x)-x*am2)
      qt=dsqrt(qt2)
      call pscs(bcos,bsin)
      ep3(3)=sngl(qt)*bcos  
      ep3(4)=sngl(qt)*bsin   
      !~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
      ! wpt(1) ... lc+ of ladder = lc+ upper EoL parton
      ! wpt(2) ... lc- of ladder = lc+ lower EoL parton
      !~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
      wplc=wpt(ii)
      qt2x=qt2+am2
      eprt=max( sqrt(qt2x) ,
     .  .5d0* ( (1.d0-x)*wplc + qt2x/(1.d0-x)/wplc ) )  !PR eq (6.16) 
      pl=((1.d0-x)*wplc-eprt)*dble(3-2*ii)
      plprt=max(1.d-3,eprt**2-qt2x) !plprt=0 gives NaN in psdefrot
      ep3(1)=sngl(eprt)
      ep3(2)=sngl(dsqrt(plprt))*sign(1.,pl)
      wpt(ii)=wpt(ii)*x !update wpt
      qp2max=max(qcdlambda(noflav(sngl(qt2))),qt2x) 
      if(iq1.eq.10) qp2max=0. 
      zeta=sqrt(qp2max/si)/sqrt(x*(1.-x))
      if(iq1.eq.0)then 
        iq2ini=9  !gluon is 9 in timsh1
        jo=iq2ini
        if(zeta.gt.zetacut)jo=-jo
      else
        iq2ini=iq1
        jo=iq2ini
      endif
      ey(1)=1.
      ey(2)=1.
      ey(3)=1.
      do i=1,4
        a=ep3(i) 
        if(.not.(a.gt.0..or.a.le.0.))stop'NaN catch 03092019'
        if(a.gt.1e35 .or. a.lt.-1e35 )stop'Infinity catch 01092019'
      enddo
      if(iq1.ne.10)then
        call psdefrot(ep3,s0xh,c0xh,s0h,c0h)
        !--------------------------------------------------------
        call timsh1(qp2max,sngl(eprt),iq2ini,jo,ey,0.,1.,0.,1.
     .                                  ,s0xh,c0xh,s0h,c0h) 
        !--------------------------------------------------------
        amprt=pprt(5,1)**2 !the mass computed in timsh1
        !~~~~~~~~~~~~checks~~~~~~~~~~~
        a=pprt(5,1)
        if(.not.(a.gt.0..or.a.le.0.))stop'NaN catch 04092019'
        do i=1,4
          a=pprt(i,1)
          if(.not.(a.gt.0..or.a.le.0.))then
            ncount5=ncount5+1 
            if(mod(log10(float(ncount5)),1.0).eq.0.0)then
              lcount5=log10(float(ncount5))+1
              write(*,'(a,i1,2a,i2)')
     .        'WARNING ',lcount5,' psahot; NaN after timsh1 corrected'
            endif
            pprt(1,1)=ep3(3)
            pprt(2,1)=ep3(4)
            pprt(3,1)=ep3(2)
            pprt(4,1)=ep3(1)
          endif  
        enddo
        !~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
      endif                                   
      iEOL(ii)=i2            !update EOL flavor

      !parton registration

      if(iq1.eq.10) then         
        nprod=nprod+1
        nptl=nptl+1
        desptl(nptl)=factgam !???????????????  !weight for individual photon
        pptl(1,nptl)=ep3(3)
        pptl(2,nptl)=ep3(4)
        pptl(3,nptl)=ep3(2)
        pptl(4,nptl)=ep3(1)
        pptl(5,nptl)=0
        idptl(nptl)=10
        iorptl(nptl)=iptl
        istptl(nptl)=0
        jorptl(nptl)=0
        do i=1,4
          xorptl(i,nptl)=bx(i)
        enddo
        tivptl(1,nptl)=bx(5) !bg ga    bx=coordo
        tivptl(2,nptl)=bx(6)
        ityptl(nptl)=73       !bg ga   73 for fragmentation photons
        ifrptl(1,nptl)=0      !bg ga
        ifrptl(2,nptl)=0        
        qsqptl(nptl)=0.        
        zpaptl(1,nptl)=zpaptl(1,iptl)
        zpaptl(2,nptl)=zpaptl(2,iptl)
      else
         !-------------------------------------------------------------
         call putInfoPsreti(ipto,ii,39
     .   ,ncc(1,1),ncc(2,1),ncc(1,2),ncc(2,2))
         call psreti(ncr,jq,1,ey,0.,1.,0.,1.,s0xh,c0xh,s0h,c0h,iret) !updates color conn list ncj, updates ncr
         !-------------------------------------------------------------
         if(iret.ne.0)stop'ERROR 24072019'
      endif 

      !color connection

      !----------------------------!
      !      ncc(1)    ncc(2)      ! moves to left
      !                            ! in case of   
      !      ncc(4)    ncc(3)      ! vanishing elements 
      !----------------------------!

      if(jt.eq.1)then
        ncc(jq,ii)=ncr(2)             
      elseif(jt.eq.2)then            
        ncc(jq,ii)=ncc(1,ii)         
        ncc(3-jq,ii)=ncr(1)          
      elseif(jt.eq.3)then            
        ncc(1,ii)=ncr(2)
      elseif(jt.eq.4)then
        ncc(1,ii)=ncc(3-jq,ii)
        ncc(2,ii)=0  !KW1904: was missing
      endif

      ! <------------------end  emissions-----------------------
      enddo
      enddo

      call timer(iutime)
      tidi=iutime(3)-tiu3 + (iutime(4)-tiu4)*1d-3
      if(ltime)write(ifmt,'(a,i20)')'cpu time psahot 4/5'
     .,nint(tidi*1000)
      tiu3=iutime(3)
      tiu4=iutime(4)

      !KKKKKKKKKKKKKKKKKKKKKKKKKKKKKKKKKKKKKKKKKKKKKKKKKKKKKKKKKKK
      !KK                                                       KK
      !KK                realize Born                           KK
      !KK                                                       KK
      !KKKKKKKKKKKKKKKKKKKKKKKKKKKKKKKKKKKKKKKKKKKKKKKKKKKKKKKKKKK

      ienvi=20

      !-------------------------------------------------------
      !at this point, we know  
      !kinematic vars double    sBorn, pt2Born, 
      !kinematic vars sngl      tBorn, xpBorn, xmBorn, yBorn
      !the flavors              mIB,nIB,  mOB,nOB (in-Born,out-Born)
      !the ladder 4-momentum    pLadder  (SO convention) 
      !Pom index (from nppr)    iptl
      !Pom type (0,1,2,3)       iqq 
      !Pom formation point      bx
      !------------------------------------------------------- 

      !~~~~~~~~~~check~~~~~~~~~
        if(abs(mOB).eq.6.or.abs(nOB).eq.6)then
          write(ifmt,'(a,$)')'ERROR psahot Born Top detected ;  '
          write(ifmt,*)'mOB nOB =',mOB,nOB
          stop
        endif
      !~~~~~~~~~~~~~~~~~~~~~~~~

      klas=igetKlas( mIB,nIB,  mOB,nOB )
      call getBornPtr2(klas,sBorn,dble(tBorn) , pt2Born) !compute pt2Born

      !ept is 4 momentum of IB parton pair 
      !ept(1)=(pLadder(1)*(xpBorn+xmBorn) + pLadder(2)*(xpBorn-xmBorn))/2     !(wplus+wminus)/2
      !ept(2)=(pLadder(1)*(xpBorn-xmBorn) + pLadder(2)*(xpBorn+xmBorn))/2     !(wplus-wminus)/2
      wplus = (pLadder(1)+pLadder(2)) * xpBorn        
      wminus= (pLadder(1)-pLadder(2)) * xmBorn  
      ept(1)=(wplus+wminus)/2d0
      ept(2)=(wplus-wminus)/2d0
      ept(3)=pLadder(3)
      ept(4)=pLadder(4)
      ept1(1)=.5d0*wplus  
      ept1(2)=.5d0*wplus 
      ept1(3)=pxh1            
      ept1(4)=pyh1
      ept2(1)=.5d0*wminus
      ept2(2)=-.5d0*wminus  
      ept2(3)=pxh2            
      ept2(4)=pyh2
      sBornX=sBorn
      sBorn=s1234(ept)
      if(pt2Born.gt.sBorn/4)then
        write(ifmt,'(a,f10.5)')'WARNING pt2Born correction',
     .  abs(pt2Born-(sBorn/4))/(sBorn/4)
        pt2Born=sBorn/4
      endif
      !print*,'Born energy squared',sBorn,sBorn/sBornX
      !print*,'Born LAB',
      !.sngl(ept(1)),sngl(ept(2)),sngl(ept(3)),sngl(ept(4))
      !print*,'  P2 LAB',ept2
      !print*,'  P1 LAB',ept1
      call psdeftr(sBorn,ept,ey)                              ! ey - boost to Born cms 
      call pstrans(ept1,ey,-1)
      !call pstrans(ept2,ey,-1)
      !print*,'Born CMS',
      !.sngl(ept(1)),sngl(ept(2)),sngl(ept(3)),sngl(ept(4)),' ey',ey
      !print*,'  P2 CMS',ept2
      !print*,'  P1 CMS',ept1
      call psdefrot(ept1,s0xi,c0xi,s0i,c0i)                ! s0xi,c0xi,s0i,c0i - rotation to z-axis
      !print*,'  P1 ROT',ept1,' euler',s0xi,c0xi,s0i,c0i
      qq0=max(qqs(1),qqs(2))    !initial PE virtualities, before emissions
      call xpprbor_set(jsplit,ncolp,kcol, sngl((ept(1)+ept(2))/plc) )
      call xmprbor_set(jsplit,ncolp,kcol, sngl((ept(1)-ept(2))/plc) )
      call idbor_set(1,jsplit,ncolp,kcol, mIB )
      call idbor_set(2,jsplit,ncolp,kcol, nIB )
      if(jsplit.eq.1)nemispr(1,ncolp,kcol)=nemis(1)
      if(jsplit.eq.1)nemispr(2,ncolp,kcol)=nemis(2)
      if(jsplit.eq.1)q2bor(ncolp,kcol)=qq 
      shatpr(ncolp,kcol)=max(shatpr(ncolp,kcol),sBorn)
      shat=sBorn
      !if(ipomtype.lt.0)then
      !  print*,'TEST bf Born',max(qqs(1),qqs(2)),pt2Born
      !endif

      !---------------------------------------------------------
      call psabor(sBorn,qq0,pt2Born,mIB,nIB,mOB,nOB,ncc,min(iqq,0),iptl
     .       ,bx,jtp,ey,s0xi,c0xi,s0i,c0i,ipomtype,qqupp)
      !---------------------------------------------------------
c      if(iqq.eq.0)write(40,*)qq,sss,zpm,sBorn,gbyj,gb0,gbymax
 
      call psphwr(iqq,mIB,nIB,mOB,nOB)   !out Born partons (pprtx) to /cptl/

      do ii=1,2
        ptprboo_=sqrt(pprtx(1,ii)**2+pprtx(2,ii)**2)
        rapprboo_=0.
        if(abs(pprtx(3,ii)).gt.0..and.ptprboo_.gt.0.)then
          rapprboo_=sign(1.,pprtx(3,ii))*
     *       log((sqrt(pprtx(3,ii)**2+ptprboo_**2)
     *       +abs(pprtx(3,ii)))/ptprboo_)
        endif
        if(ioTestFact.eq.-6.and.jqq.ne. 0      !not qg
     . .or.ioTestFact.eq.-7.and.jqq.ne. 1)then  !not qq 
           ptprboo_=0.
           rapprboo_=100000.
        endif
        if(ihqTestFact.eq.1)then !outBorn at least one HQ
          if(  abs(idprtx(ii)).ne.4.and.abs(idprtx(3-ii)).ne.4
     .    .and.abs(idprtx(ii)).ne.5.and.abs(idprtx(3-ii)).ne.5)then
           ptprboo_=0.
           rapprboo_=100000.          
           !else
           !if(ii.eq.1)print*,'HQ Born',mEOL,nEOL,mIB,nIB,mOB,nOB
          endif
        elseif(ihqTestFact.eq.2)then !outBorn at least one b
          if(  abs(idprtx(ii)).ne.5.and.abs(idprtx(3-ii)).ne.5)then
           ptprboo_=0.
           rapprboo_=100000.          
          endif
        elseif(ihqTestFact.eq.3)then !outBorn at least one c no b
          if( (abs(idprtx(ii)).ne.4.and.abs(idprtx(3-ii)).ne.4)
     .    .or.(abs(idprtx(ii)).eq.5.or.abs(idprtx(3-ii)).eq.5)  )then
           ptprboo_=0.
           rapprboo_=100000.          
          endif
        elseif(ihqTestFact.eq.4)then !outBorn one c one b
          if(.not.(
     .        (abs(idprtx(ii)).eq.4.and.abs(idprtx(3-ii)).eq.5)
     .    .or.(abs(idprtx(ii)).eq.5.and.abs(idprtx(3-ii)).eq.4) ))then
           ptprboo_=0.
           rapprboo_=100000.          
          endif
        endif
        call ptprboo_set(ii,jsplit,ncolp,kcol, ptprboo_ )
        call rapprboo_set(ii,jsplit,ncolp,kcol, rapprboo_ )
      enddo

      !KKKKKKKKKKKKKKKKKKKKKKKKKKKKKKKKKKKKKKKKKKKKKKKKKKKKKKKKKKK
      !KK                                                       KK
      !KK              TEST7                                    KK
      !KK                                                       KK
      !KKKKKKKKKKKKKKKKKKKKKKKKKKKKKKKKKKKKKKKKKKKKKKKKKKKKKKKKKKK

      if(ioTestFact.eq.7)then 
        jtest=0
        if(jtest.eq.1)then
        pt=sqrt(pt2Born)        
        call ptprboo_get(1,jsplit,ncolp,kcol, pt1 )
        call ptprboo_get(2,jsplit,ncolp,kcol, pt2 )
        yPairCMS=0.5*log(xpBorn/xmBorn)
        yStar=yBorn-yPairCMS
        yBorn2=yPairCMS-yStar
        rap1=yLadder+yBorn   !ignore rotations for the moment
        rap2=yLadder+yBorn2
        call rapprboo_get(1,jsplit,ncolp,kcol, y1 )
        call rapprboo_get(2,jsplit,ncolp,kcol, y2 )
        ymi1=min(rap1,rap2)
        ymi2=min(y1,y2)
        yma1=max(rap1,rap2)
        yma2=max(y1,y2)
        write(ifmt,'(a,4f7.3,2x,3f7.3,1x,4i3,4x,4i4)')'psahot TEST1  '
     .  ,ymi1,ymi2,yma1,yma2,pt,pt1,pt2
     .  ,nint((ymi1-ymi2)/ymi2*100)
     .  ,nint((yma1-yma2)/yma2*100)
     .  ,nint((pt-pt1)/pt1*100)
     .  ,nint((pt-pt2)/pt2*100),mIB,nIB,mOB,nOB
        endif
      endif

      !KKKKKKKKKKKKKKKKKKKKKKKKKKKKKKKKKKKKKKKKKKKKKKKKKKKKKKKKKKK
      !KK                                                       KK
      !KK   construct parton chains following color flow        KK
      !KK                                                       KK
      !KKKKKKKKKKKKKKKKKKKKKKKKKKKKKKKKKKKKKKKKKKKKKKKKKKKKKKKKKKK


      if(nj.ne.0.)then
        do i=1,nj
          do l=1,6
            bxj(l,i)=bx(l)
          enddo
          iorj(i)=iptl
        enddo
      endif

      !===================
      call psjarr(jtp,jfl)
      !===================

      !~~~~~ifav redo~~~~~
      if(idpomr.lt.0.and.ntrysplit.lt.20)then
        call getistptl(nptl2+1,ist1)
        call getistptl(nptl2+2,ist2)
        if(ist1.ne.25.or.ist2.ne.25)stop'ERROR 230816'
        nheavy=0
        do i=nptl2+3,nptl
          call getidptl(i,idi)
          call getistptl(i,isti)
          if(abs(idi).eq.4.or.abs(idi).eq.5)nheavy=nheavy+1
        enddo
        if((ifav.ge.4.and.nheavy.eq.0)
     .  .or.(ifav.lt.4.and.nheavy.gt.0))then
          ifav=1
          njstr=nj
          goto 2
        endif
        !if(nheavy.gt.0)then !.or.ntrysplit.gt.2)then
        !  print*,'----------------',ifav,mOB,nOB,ntrysplit,nheavy
        !  !do i=nptl2+3,nptl
        !  !  call getidptl(i,idi)
        !  !  call getistptl(i,isti)
        !  !  print*,i,idi,isti
        !  !enddo
        !endif
      endif
      !~~~~~ifav redo END~~~~~


      !~~~~~~~~~~~~~~~~~~~TEST~~~~~~~~~~~~~~~~~~~~~~~
      !sum=0
      !do i=nptl1+1,nptl
      !  if(istptl(i).eq.20)then
      !    sum=sum+pptl(4,i)
      !    print*,i,idptl(i),pptl(4,i),sum,pptl(4,iptl),sBorn
      !  endif
      !enddo
      !stop'TEST'
      !~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
      
      if(ish.ge.6)then
        esum=0.
        psum1=0.
        psum2=0.
        psum3=0.
        do i=nptl1+1,nptl
c        print *,i,idptl(i),istptl(i)
          if(istptl(i).eq.20)then
            esum=esum+pptl(4,i)
            psum1=psum1+pptl(1,i)
            psum2=psum2+pptl(2,i)
            psum3=psum3+pptl(3,i)
          endif
        enddo
        write(ifch,*)'4-Momentum check (2)',esum,pptl(4,iptl)
     .       ,esum/pptl(4,iptl),psum3,pptl(3,iptl)
     .       ,(psum1**2+psum**2+(esum+psum3)*(esum-psum3))/ss
      endif

      if(jfl.gt.0)then
        !if(jfl.ne.3)stop'ERROR 07042019'
        !jfl=3: chain too light -> redo complete psahot
        igo=9
        njstr=nj
        if(idpomr.lt.0.and.ntrysplit.lt.20)then
          if(abs(mOB).eq.4.and.abs(nOB).eq.4)vparam(93)=vparam(93)+1
          goto 2 !only for sat poms
        endif
        igo=12
        goto 1 !redo psahot
      endif

      !KKKKKKKKKKKKKKKKKKKKKKKKKKKKKKKKKKKKKKKKKKKKKKKKKKKKKKKK
      !KK                  finish                            KK
      !KKKKKKKKKKKKKKKKKKKKKKKKKKKKKKKKKKKKKKKKKKKKKKKKKKKKKKKK

      call gbpom_set(1,jsplit,ncolp,kcol, sngl(gbymin) )
      call gbpom_set(2,jsplit,ncolp,kcol, sngl(gbymax) )

      xxp1sum=xxp1sum+xxp1pom
      xyp1sum=xyp1sum+xyp1pom
      xxp2sum=xxp2sum+xxp2pom
      xyp2sum=xyp2sum+xyp2pom
      xxm1sum=xxm1sum+xxm1pom
      xym1sum=xym1sum+xym1pom
      xxm2sum=xxm2sum+xxm2pom
      xym2sum=xym2sum+xym2pom

      !update idhpr
      if(idpomr.eq.0)idhpr(ncolp,kcol)=10*jqq

      njstr=nj
      if(idpomr.lt.0.and.jsplit.lt.kxsplit)goto 666 !only for sat poms      

      xxp1pr(ncolp,kcol)=xxp1sum
      xyp1pr(ncolp,kcol)=xyp1sum
      xxp2pr(ncolp,kcol)=xxp2sum
      xyp2pr(ncolp,kcol)=xyp2sum
      xxm1pr(ncolp,kcol)=xxm1sum
      xym1pr(ncolp,kcol)=xym1sum
      xxm2pr(ncolp,kcol)=xxm2sum
      xym2pr(ncolp,kcol)=xym2sum

      idp1pr(ncolp,kcol)=idp1pom   !String end type, from ProSeTy(k,n)
      idp2pr(ncolp,kcol)=idp2pom   ! will be modified, but only for val
      idm1pr(ncolp,kcol)=idm1pom   
      idm2pr(ncolp,kcol)=idm2pom

      goto 18

  16  continue
      !if(idpomr.lt.0)print*,'Sat Pom SKIP',jsplit,imark

  18  call utprix('psahot',ish,ishini,3)

      jsplitpom(ncolp,kcol)=jsplit    
      fxsplit=1.

      call timer(iutime)
      tidi=iutime(3)-tiu3 + (iutime(4)-tiu4)*1d-3
      if(ltime)write(ifmt,'(a,i25,2i9)')'cpu time psahot 5/5'
     .,nint(tidi*1000),idpomr,laddtype

      njstr=nj
      q2fin=q2finsave
      factbgg=1.     !for safety
      factbqq=1.
      !print*,'rejection counts 1 41 10 11',ntryhot,ktry41,ktry10,ktry11

      ienvi=99 !psahot exit
      return

      end

c-----------------------------------------------------------------------
      subroutine checkCS(q1,q2,qqcut,sis,j1,l1
     .,s1,s2,s3,s4,sTPex,ihit,iforce)
c-----------------------------------------------------------------------
      integer ihit(4)
      real sar(4)
      real sTPex(4)
      sar(1)=s1  
      sar(2)=s2 
      sar(3)=s3
      sar(4)=s4
      sTPex(1) = psborn(q1,q2,qqcut,sis,j1,l1,0,0)
      sTPex(2) = psjet1(q1,q2,qqcut,sis,j1,l1,0)
      sTPex(3) = psjet1(q2,q1,qqcut,sis,l1,j1,0)
      sTPex(4) =  psjet(q1,q2,qqcut,sis,j1,l1,0)
      do i=1,4
        ihit(i)=0
        dif=abs(sar(i)-sTPex(i))/max(1e-10,sTPex(i))
        if(dif.gt.0.1)ihit(i)=1
      enddo
      call getMonitorFileIndex(ifmtx)
      if(ihit(1)+ihit(2)+ihit(3)+ihit(4).gt.0.or.iforce.eq.1)then
        if(iforce.eq.0)write(ifmtx,'(a,$)')
     .  'WARNING psahot HARD sTYP unprecise'
        if(iforce.eq.1)write(ifmtx,'(a,$)')'  sTYP & exact'
        do i=1,4
          if(ihit(i).ne.0)
     .    write(ifmtx,'(2x,i2,2f10.6,$)')i,sar(i),sTPex(i)
        enddo
        write(ifmtx,*)' '
      endif
      end

c-----------------------------------------------------------------------
      subroutine storeTL(ii,n,t1,t2,x1,x2,i1,i2)
c-----------------------------------------------------------------------
      double precision x1,x2,x1arr,x2arr
      parameter (maxemiss=100)
      common /ctimelike/ t1arr(2,maxemiss), t2arr(2,maxemiss)
     .                 , x1arr(2,maxemiss), x2arr(2,maxemiss)
     .                 , i1arr(2,maxemiss), i2arr(2,maxemiss)
      t1arr(ii,n)=t1
      t2arr(ii,n)=t2
      x1arr(ii,n)=x1
      x2arr(ii,n)=x2
      i1arr(ii,n)=i1
      i2arr(ii,n)=i2
      end

c-----------------------------------------------------------------------
      subroutine getTL(ii,n,t1,t2,x1,x2,i1,i2)
c-----------------------------------------------------------------------
      double precision x1,x2,x1arr,x2arr
      parameter (maxemiss=100)
      common /ctimelike/ t1arr(2,maxemiss), t2arr(2,maxemiss)
     .                 , x1arr(2,maxemiss), x2arr(2,maxemiss)
     .                 , i1arr(2,maxemiss), i2arr(2,maxemiss)
      t1=t1arr(ii,n)
      t2=t2arr(ii,n)
      x1=x1arr(ii,n)
      x2=x2arr(ii,n)
      i1=i1arr(ii,n)
      i2=i2arr(ii,n)
      end


c------------------------------------------------------------------------
      subroutine psjarr(jt,jfl)
c-----------------------------------------------------------------------
c
c   rearranging partons into chains F-g-g-g-g-...-g-g-L
c
c      where F is the first is a "non-gluon" (q,qbar,qq,qqbar)
c      and L as well a "non-gluon", but of opposite color compared to F
c      g is a gluon
c
c   and writing the chains into /cptl/ (one entry per parton)
c
c   jt  - process index
c   jfl - flag for rejection ( jfl > 0 )
c
c-----------------------------------------------------------------------
c
c in:
      parameter (mjstr=20000)
      common /psar29/ eqj(4,mjstr),iqj(mjstr),ncj(2,mjstr),ioj(mjstr),nj
c
c parton properties, for parton index k, in  /psar29/
c   eqj(1:4,k) - 4-momentum (qgs) 
c   bxj(1:4,k) - coordinates of formation point
c   iqj(k) - parton id 
c   ncj(1:2,k) - indices of colour connected partons for TL partons 
c nj - number of partons
c
c out:
c      parton chains written into /cptl/  in aaa.h
c
c-----------------------------------------------------------------------
      dimension mark(mjstr)
      double precision  ept(4)
#include "aaa.h"
#include "sem.h"
      data ncount/0/
      save ncount 
      data ncount2/0/
      save ncount2 
      real ksequ(100)
      nsequ=0

      if(nj.eq.0)then
        jfl=1
        write(ifmt,*)'ERROR in psjarr: nj=0'
        !stop'ERROR 05042019c'
        goto 999
      endif
      njpar=0

      if(nj.le.4.and.jt.eq.14)then
        if(nj.eq.4)then
          ncj(1,1)=4
          ncj(1,4)=1
          ncj(1,2)=3
          ncj(1,3)=2
          goto 15
         elseif(nj.lt.4)then
          write(*,*)'PROBLEME : qq->gamma+gamma and gluon string end'
        endif
      endif

      !------ complete color connection array ncj ----------

      do k=1,nj
        if(iqj(k).eq.0)then   !gluon must have two neighbours
          if(ncj(1,k).eq.0)then !first neigbour missing
            do kk=1,nj          !look to which parton it is connected
              if(ncj(1,kk).eq.k)then
                if(ncj(2,k).ne.kk)ncj(1,k)=kk !if not already connected : add connection
              elseif(ncj(2,kk).eq.k)then
                if(ncj(1,k).ne.kk)ncj(1,k)=kk
              endif
            enddo
           endif
           if(ncj(2,k).eq.0)then !second neigbour missing
            do kk=1,nj
              if(ncj(2,kk).eq.k)then
                if(ncj(1,k).ne.kk)ncj(2,k)=kk
              elseif(ncj(1,kk).eq.k)then
                if(ncj(2,k).ne.kk)ncj(2,k)=kk
              endif
            enddo
           endif
         else !quark must have at least one neigbhour
          if(ncj(1,k).eq.0)then !first neigbour missing
           do kk=1,nj          !look to which parton it is connected
              if(ncj(1,kk).eq.k)then
                if(ncj(2,k).ne.kk)ncj(1,k)=kk !if not already connected : add connection
              elseif(ncj(2,kk).eq.k)then
                if(ncj(1,k).ne.kk)ncj(1,k)=kk
              endif
            enddo
           endif
        endif
      enddo

15    continue ! <<<<<<------------------

      if(ish.ge.3)then
        write (ifch,*)'psjarr: nj',nj
        do k=1,nj
          if(iqj(k).ne.0)ncj(2,k)=0
          write(ifch,'(a,i4)')' parton',k
          write(ifch,'(i6,2x,4e10.3,2x,2i3)')iqj(k)
     *    ,(eqj(j,k),j=1,4),(ncj(j,k),j=1,2)
        enddo
      endif

      jfl=2
      do i=1,nj
        mark(i)=1 !available parton
      enddo

      nptl0=nptl  !nptl will be updated for each final jet, for the last one: nptl=nptl0+nj
                  !   and the last has to be a quark state

1     continue ! <<<<<--------------------------- 

      do ij=1,nj
        if(mark(ij).ne.0.and.iqj(ij).ne.0)goto 2 !bg so we begin by a quark state
      enddo

2     continue ! <<<<<--------------------------- 

      jfirst=1
      if(iabs(iqj(ij)).le.2)then       !u,d -> pion mass
        am1=amhadr(1)
      elseif(iabs(iqj(ij)).eq.4)then   !s -> kaon mass
        am1=amhadr(3)
      elseif(iabs(iqj(ij)).eq.40)then !bg charm
        am1=qcmass
      elseif(iabs(iqj(ij)).eq.50)then !bg bottom
        am1=qbmass                    !bg
      else                            !diquark -> nucleon mass
        am1=amhadr(2)
      endif
      do i=1,4
        ept(i)=0.d0
      enddo

3     continue

      if(ij.le.0)then
        ncount2=ncount2+1
        if(ncount2.le.5)then
          write(ifmt,'(a)')
     .    'ERROR in psjarr: ij .le. 0'
        endif
        jfl=6
        goto 999 
        !stop'ERROR 05042019d' 
      endif
      !===============
       call pspawr(ij) !nptl=nptl+1 in this function
      !===============
      mark(ij)=0  !parton no more available
      nsequ=nsequ+1
      if(nsequ.le.100)ksequ(nsequ)=ij
      do i=1,4
        ept(i)=ept(i)+dble(eqj(i,ij))
      enddo

      if(iqj(ij).ne.0)then !---------- no gluon

        if(jfirst.ne.1)then  !---------- L parton in chain F-g-g-g-g-g-L 

          if(iabs(iqj(ij)).le.2)then       !u,d -> pion mass
            am2=amhadr(1)
          elseif(iabs(iqj(ij)).eq.4)then   !s -> kaon mass
            am2=amhadr(3)
          elseif(iabs(iqj(ij)).eq.40)then !bg charm
            am2=qcmass
          elseif(iabs(iqj(ij)).eq.50)then !bg bottom
            am2=qbmass
          else                            !diquark -> nucleon mass
            am2=amhadr(2)
          endif
          amj=(am1+am2+stmass)**2
          sm=psnorm(ept)
          if(sm.lt.amj)then
            jfl=3
            if(ish.ge.6)write(ifch,'(a,i3,4f8.4,i5)')'chain rejected'
     .      ,jfl,sm,amj,am1,am2,iqj(ij)
            nptl=nptl0
            goto 999
          endif

          if(nptl-nptl0.lt.nj)then !bg if nptl-nptl0=nj-> all jets registered-> end
            igo=13
            goto 1    !bg else one looks for a new string and so for a new quark 
          else
            if(ish.ge.3)then
              write (ifch,*)'psjarr: nptl',nptl
              do k=nptl0+1,nptl
                write(ifch,'(a,i4)')' particle',k
                write(ifch,'(i5,2x,6e10.3)')idptl(k)
     *          ,(pptl(j,k),j=1,5),qsqptl(k)
              enddo
            endif
            jfl=0
            goto 999
          endif

        else  !  ------- L parton in chain L-g-g-g-g-g-R 

          jfirst=0
          njpar=ij
          ij=ncj(1,ij)
          goto 3   !bg know looking for gluons

        endif

      else ! ------- gluon

        if(ncj(1,ij).eq.njpar)then
          njdau=ncj(2,ij)
        elseif(ncj(2,ij).eq.njpar)then !KW1902
          njdau=ncj(1,ij)
        else !KW1902
c          do i=1,min(nsequ,100)
c            ijk=ksequ(i)
c            write(ifmt,*)' ------> ',iqj(ijk),ijk,(ncj(j,ijk),j=1,2)
c          enddo
          ncount=ncount+1
          if(ncount.le.5)then
            write(ifmt,'(a)')
     .      'ERROR in psjarr: Impossible color connections'
          endif
          jfl=5
          goto 999
          !stop'ERROR 05042019e' 
        endif
        njpar=ij
        ij=njdau
        goto 3     !bg now looking for gluons

      endif

  999 continue
      !jfl=0: OK 
      !jfl=3: chain too light, this may happen. 
      !all other jfl values represent errors, to be fixed! 

      !if(jfl.ne.3.and.jfl.ne.0)stop'ERROR 05042019f' 

      return
      end

c------------------------------------------------------------------------
      subroutine eqjcheck(kk,kj)
c------------------------------------------------------------------------
#include "aaa.h"
      parameter (mjstr=20000)
      common /psar29/ eqj(4,mjstr),iqj(mjstr),ncj(2,mjstr),ioj(mjstr),nj
      common /psar30/ iorj(mjstr),ityj(mjstr),bxj(6,mjstr),q2j(mjstr)

      if(kj.le.0)then
        write(ifmt,'(a,i3,a,i3)')'WARNING eqj  kk =',kk,'  kj =',kj
      endif
      do i=1,4
        xx=eqj(i,kj)
        if(.not.(xx.le.0..or.xx.ge.0.))then
          write(ifmt,'(a,i3,a,$)')'ERROR eqj NaN catch'
     .    ,kk,'   i kj eqj(,):'
          write(ifmt,*)i,kj,eqj(i,kj)
          stop
        endif
      enddo
      end

c------------------------------------------------------------------------
      subroutine pspawr(kj)
c-----------------------------------------------------------------------
c pspawr - writing final parton kj into particle list
c------------------------------------------------------------------------
c Input: chain of partons
      parameter (mjstr=20000)
      common /psar29/ eqj(4,mjstr),iqj(mjstr),ncj(2,mjstr),ioj(mjstr),nj
      common /psar30/ iorj(mjstr),ityj(mjstr),bxj(6,mjstr),q2j(mjstr)
      common /psar31/ ktyj(mjstr)
c------------------------------------------------------------------------
c nj - number of partons
c eqj(1:4,k) - 4-momentum (qgs) for k-th parton;
c bxj(1:4,k) - coordinates for k-th parton formation point;
c iqj(k) - ID (qgs) for k-th parton;
c ncj(1:2,k) - colour connected partons indexes for k-th parton;
c ioj(j) - flavor of mother parton
c ktyj(j) - parton type (from SL = 1,2, from Born = 11,12,13)
c------------------------------------------------------------------------
#include "aaa.h"
#include "sem.h"
      common/ciptl/iptl
      parameter (ntim=1000)                           !ala
      common/iortimsh2/iordur(ntim),iodur(ntim)       !ala
      common/cptldur/nptldur(2)                       !ala

      call eqjcheck(10,kj)
      nptl=nptl+1
      if(ish.ge.9)write (ifch,*)'nptl,kj (sto)',nptl,kj
      if(nptl.ge.mxptl.or.kj.le.0)then
       write (ifmt,*)'nptl,kj',nptl,kj
       call alist('Error in pspawr: nptl or kj out of bounds &',1,nptl)
       call utstop('nptl or kj out of bounds&')
      endif

      if(ifrptl(1,iptl).eq.0)ifrptl(1,iptl)=nptl !bg ifrptl= index of first child
      ifrptl(2,iptl)=nptl

      pptl(1,nptl)=eqj(3,kj)
      pptl(2,nptl)=eqj(4,kj)
      pptl(3,nptl)=eqj(2,kj)
      pptl(4,nptl)=eqj(1,kj)
      pptl(5,nptl)=0.
      idptl(nptl)=psidd(iqj(kj))
      if(iabs(idptl(nptl)).eq.4 .or.iabs(idptl(nptl)).eq.5)then    !bg
        pmass2=eqj(1,kj)**2-eqj(3,kj)**2-eqj(2,kj)**2-eqj(4,kj)**2 !bg
        if(pmass2.ge.0.)then                                       !bg
          pptl(5,nptl)=sqrt(pmass2)                                !bg
        else                                                       !bg
          pptl(5,nptl)=0                                           !bg
        endif                                                      !bg
      endif                                                        !bg
      iorptl(nptl)=iorj(kj)
      if(iodur(kj).gt.0)then                              !ala
       radptl(nptl)=nptldur(iodur(kj))                    !ala
      else                                                !ala
       radptl(nptl)=0                                     !ala
      endif                                               !ala
      !print*,'hard origin = iaaptl for partons: ',kj,nptl,iaaptl(nptl)
      jorptl(nptl)=ioj(kj)
      ifrptl(1,nptl)=ktyj(kj)
      istptl(nptl)=20
      do i=1,4
        xorptl(i,nptl)=bxj(i,kj)
      enddo
      tivptl(1,nptl)=bxj(5,kj)
      tivptl(2,nptl)=bxj(6,kj)
      ityptl(nptl)=ityj(kj)
      zpaptl(1,nptl)=zpaptl(1,iorj(kj))
      zpaptl(2,nptl)=zpaptl(2,iorj(kj))
c register to which big string the particle belongs to
      qsqptl(nptl)=q2j(kj)
       !    write(*,'(a,2i4,i6,f8.3)')'.... ',kj,nptl,idptl(nptl)
       !*     ,sqrt(pptl(1,nptl)**2+pptl(2,nptl)**2)
      return
      end

c-----------------------------------------------------------------------
      subroutine psreti(ncr,jort,nfprt
     .          ,ey,s0xi,c0xi,s0i,c0i,s0xh,c0xh,s0h,c0h,iret)
c-----------------------------------------------------------------------
c Updates the chain of partons (= kinks) following color flow
c starting from the time-like parton nfprt being the start of a TL cascade
c (either from Born or the SP cascade).
c It uses information from /cprt/, already filled from calls to tim... (TL cascade)
c----------------------------------------------------------------------- 
c Notation: jet = final state time-like partons
c----------------------------------------------------------------------- 
c ncr(i) - colour connections for the parton
c jort - color orientation for gluons (=1 if +color goes first, =-1 otherwise)

      parameter (ntim=1000)
      common/cprt/pprt(5,ntim),q2prt(ntim),idaprt(2,ntim),idprt(ntim)
     &,iorprt(ntim),jorprt(ntim),nprtj
      common /ep3g/ngam,ep3gam(5,ntim)       !bg ga
      common/ciptl/iptl                      !bg ga
      dimension ey(3)
      common/cidpomr/idpomr

c nprtj - number of partons in the jet (including virtual ones)
c pprt - 5-momenta for the partons
c idprt - parton id
c iorprt - parent parton position in the list
c idaprt - daughter partons positions in the list
c ey,s0xi,c0xi,s0i,c0i,s0xh,c0xh,s0h,c0h - boost, rotations Born

c----------------------------------------------
c output: chain of partons, updated
c----------------------------------------------
      parameter (mjstr=20000)
      common /psar29/ eqj(4,mjstr),iqj(mjstr),ncj(2,mjstr),ioj(mjstr),nj
      common /psar30/ iorj(mjstr),ityj(mjstr),bxj(6,mjstr),q2j(mjstr)
      common /psar31/ ktyj(mjstr)
c-----------------------------------------------
c nj - number of final partons in jet
c eqj(i,j) - 4-momentum for the final parton j
c iqj(j) - flavour for the final parton j
c ncj(m,j) - colour connections for the final parton j
c ioj(j) - flavor of mother parton
c ktyj(j) - parton type (from SL = 1,2, from Born = 11,12,13)
c-----------------------------------------------------------------------
      dimension ep3(4),ncr(2),nci(2,ntim)
#include "aaa.h"
#include "sem.h"
      common/cprtx/nprtjx,pprtx(5,2)
      common/prtdur/idprtx(2)                !ala
      common/iortimsh2/iordur(ntim),iodur(ntim)         !ala
      common/cptldur/nptldur(2)                         !ala
      data ncount/0/
      save ncount 
      data ncount2/0/
      save ncount2 

      call utpri('psreti',ish,ishini,4)

      call getInfoPsreti(ipto,jj,ibo,nl1,nr1,nl2,nr2)

      iret=0

      if(ish.ge.4)then !~~~~~~~~~~~~~~print~~~~~~~~~~~~~~~~~~~~~~~~~
        write(ifmt,'(80a1)')('=',i=1,80)
        write(ifmt,'(a,2i4,2(4x,2i4))')'psreti initial: ipto jj ncc =',
     .  ipto,jj,nl1,nr1,nl2,nr2
        write(ifmt,'(42x,a,2i4,4x,a,i4)')'ncr =',ncr,'jort =',jort
        do ij=1,nj
        call idtrafo_KW_SO( iKW , iqj(ij) , -1 ) 
        write(ifmt,'(i4,3x,3i4)')iKW,ncj(1,ij),ij,ncj(2,ij)
        enddo
      endif !~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~


      iprt=nfprt
      if(ish.ge.6)then
        write (ifch,*)'nprtj',nprtj
        do i=1,nprtj
          write (ifch,*)'i,ic,np,ndd',i,idprt(i),iorprt(i),
     *    idaprt(1,i),idaprt(2,i)
        enddo
      endif

      !-----------------------------------------------------------------
      !  color connection information:
      !-----------------------------------------------------------------
      !        ncr(1:2) -   color connection of one parton, 
      !                    the one creating the time-like cascade
      !                     (input to this routine)
      !-----------------------------------------------------------------
      !   ncj(1:2,k) - indices of color connected partons for TL partons 
      !              in /psar29/    nj - number of partons (max value of k)
      !               (will be updated in this routine)
      !-----------------------------------------------------------------
      !   nci(2,ntim) - color connection referring to TL cascade
      !-----------------------------------------------------------------

      !-----------------------------------------------------------------
      !  flavor information: the array iqj in /psar29/ uses 
      !              (compared to "normal id")
      !-----------------------------------------------------------------
      !  id=9     :  iqj(nj) = 0
      !  |id| < 3 :  iqj(nj) = id
      !  |ad| = 3 :  iqj(nj) = id * 4 / 3
      !  else     :  iqj(nj) = id * 10
      !-----------------------------------------------------------------

      nci(1,nfprt)=ncr(1)
      if(idprt(nfprt).eq.9)then
        nci(2,nfprt)=ncr(2)
      else
        nci(2,nfprt)=0
      endif

      if(nprtjx.eq.2)then !out Born before timelike cascade
       ep3(1)=pprtx(4,iprt)
       ep3(2)=pprtx(3,iprt)
       ep3(3)=pprtx(1,iprt)
       ep3(4)=pprtx(2,iprt)
       !print*,'PTL TIM  ',ep3
       call psrotat(ep3,s0xh,c0xh,s0h,c0h)     !rotation Born pt -- 25 partons
       !print*,'PTL ROT  ',ep3,'   ',s0xh,c0xh,s0h,c0h
       call psrotat(ep3,s0xi,c0xi,s0i,c0i)     !rotation intr pt -- 25 partons
       !print*,'PTL ROT  ',ep3,'   ',s0xi,c0xi,s0i,c0i
       call pstrans(ep3,ey,1)                 !back to LAB    -- 25 partons
       !print*,'PTL BOOST',ep3,'   ',ey
       call psmirror(ep3)
       pprtx(4,iprt)=ep3(1)
       pprtx(3,iprt)=ep3(2)
       pprtx(1,iprt)=ep3(3)
       pprtx(2,iprt)=ep3(4)
       idprtx(iprt)=idprt(iprt)
       !~~~~check~~~~
       !if(iprt.eq.2)then
       ! pt=min(pprtx(1,1)**2+pprtx(2,1)**2,pprtx(1,2)**2+pprtx(2,2)**2)
       ! if(sqrt(pt).gt.3)then
       !  phi=polar(pprtx(1,2),pprtx(2,2))-polar(pprtx(1,1),pprtx(2,1))
       !  phi=abs(phi)
       !  if(phi.gt.pi)phi=2*pi-phi
       !  write(ifmt,'(2(a,i10,4f10.2,f14.2,i7/))')
       !.   'psreti: p_OB1 =',idprtx(1),(pprtx(k,1),k=1,4),phi/pi,jj,
       !.   '        p_OB2 =',idprtx(2),(pprtx(k,2),k=1,4),phi/pi,jj
       ! endif
       !endif 
      endif

      !cbg---photon--registration --->
      if(ngam.gt.0)then
        do j=1,ngam
          nptl=nptl+1
          ep3(1)=ep3gam(1,j)
          ep3(2)=ep3gam(2,j)
          ep3(3)=ep3gam(3,j)
          ep3(4)=ep3gam(4,j)
          ityptl(nptl)=74       !bg ga 74 for final fragmentation photons.  73 for initial
          if(ep3gam(5,j).lt.0.) then !bg In tim.f, parton with Q**2>FacS2>0 -->  Q**2= -Q**2
            ep3gam(5,j)=-ep3gam(5,j)
            ityptl(nptl)=71  !bg registered as direct photon because q2> (MF)**2. MF = factorization scale
          endif
          desptl(nptl)=ep3gam(5,j)               !individual weight for photons
          call psrotat(ep3,s0xh,c0xh,s0h,c0h)
          call psrotat(ep3,s0xi,c0xi,s0i,c0i)     !rotation intr pt
          call pstrans(ep3,ey,1)
          call psmirror(ep3)
          pptl(1,nptl)=ep3(3)
          pptl(2,nptl)=ep3(4)
          pptl(3,nptl)=ep3(2)
          pptl(4,nptl)=ep3(1)
          if(ep3(1).lt.0.)then                             !bg ga sometimes when p**2 very small
            pptl(4,nptl)=sqrt(ep3(3)**2+ep3(4)**2+ep3(2)**2)
          endif
          pptl(5,nptl)=0
          idptl(nptl)=10
          iorptl(nptl)=iptl
          istptl(nptl)=0
          jorptl(nptl)=0
          do i=1,4
            xorptl(i,nptl)=xorptl(i,iptl) !bg ga voir avec Klaus
          enddo
          tivptl(1,nptl)=tivptl(1,iptl) !bg ga voir avec Klaus
          tivptl(2,nptl)=tivptl(2,iptl)
          ityptl(nptl)=74       !bg ga 74 frag photon in FSR.  73 for ISR
          ifrptl(1,nptl)=0
          ifrptl(2,nptl)=0
          qsqptl(nptl)=0.
          zpaptl(1,nptl)=zpaptl(1,iptl)
          zpaptl(2,nptl)=zpaptl(2,iptl)
        enddo
      endif
      ngam=0
      !cbg <-----

1     continue ! <<<<<-----------------------------------

      idau1=idaprt(1,iprt)
      idau2=idaprt(2,iprt)
      icp=idprt(iprt)

      if(idau1.ne.0.)then    ! parton having daughters

        icd1=idprt(idau1)
        if(icp.eq.9)then
          if(icd1.ne.9)then      !g -> qq~
            nci(1,idau1)=nci(jort,iprt)
            nci(1,idau2)=nci(3-jort,iprt)
          else                    !g -> gg
            nci(1,idau1)=nci(1,iprt)
            nci(2,idau1)=0
            nci(2,idau2)=nci(2,iprt)
            nci(1,idau2)=0
          endif
        else                      !q -> qg
          nci(1,idau1)=0
          if(icp*(3-2*jort).gt.0)then
            nci(1,idau2)=nci(1,iprt)
            nci(2,idau2)=0
          else
            nci(1,idau2)=0
            nci(2,idau2)=nci(1,iprt)
          endif
        endif
        if(ish.ge.4)then !~~~~~~~~~~~~~~print~~~~~~~~~~~~~~~~~~~~~~~~
        write (ifmt,'(a,4(2i4,3x))')'nr id  idau1/2  nci  ncj+'
     .  ,iprt,icp,idau1,idau2,(nci(k,iprt),k=1,2),(ncj(k,nj+1),k=1,2)
        endif !~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
        iprt=idau1
        goto 1

      else ! final parton, no daughters

        nj=nj+1

        if(ncj(1,nj).ne.0.or.ncj(2,nj).ne.0)then
          ncount2=ncount2+1
          if(ncount2.le.5)then
            write(ifmt,'(a)')
     .      'ERROR psreti : ncj non-zero, wrong initialization'
          endif
          iret=1
          return
        endif

        ep3(1)=pprt(4,iprt)
        ep3(2)=pprt(3,iprt)
        ep3(3)=pprt(1,iprt)
        ep3(4)=pprt(2,iprt)
        !~~~~~~~~~~check~~~~~~~~
        do i=1,4
          eqj(i,nj)=ep3(i)
        enddo
        call eqjcheck(11,nj)
        !print*,'PSRETI',nj,ep3
        !~~~~~~~~~~~~~~~~~~~~~~~
        call psrotat(ep3,s0xh,c0xh,s0h,c0h)  ! -- 21 partons
        call psrotat(ep3,s0xi,c0xi,s0i,c0i)  ! -- 21 partons
        call pstrans(ep3,ey,1)               ! -- 21 partons
        call psmirror(ep3)
        !~~~~~~~~~~check~~~~~~~~
        !print*,'      ',nj,ep3
        x=ep3(1)
        y=ep3(2)
        if(.not.(x.le.0..or.x.ge.0.).or..not.(y.le.0..or.y.ge.0.))then
          write(ifmt,'(a,$)')'ERROR psahot NaN 12a ;  '
          write(ifmt,*)'jj s0 ey sc sci ='
     .                ,jj,ey,s0xi,c0xi,s0i,c0i,s0xh,c0xh,s0h,c0h
          stop
        endif
        !~~~~~~~~~~~~~~~~~~~~~~~
        do i=1,4
          eqj(i,nj)=ep3(i)
        enddo
        call eqjcheck(12,nj)

        if(icp.eq.9)then
          iqj(nj)=0
        elseif(iabs(icp).lt.3)then
          iqj(nj)=icp
        elseif(iabs(icp).eq.3)then
          iqj(nj)=icp*4/3
        else
          iqj(nj)=icp*10
        endif

        ityj(nj)=ibo

        ioj(nj)=jorprt(iprt) !flavor of mother parton
        q2j(nj)=q2prt(iprt)
        ktyj(nj)=jj
        
        if(nprtjx.eq.2)then !out Born before timelike cascade  !ala
          iodur(nj)=iordur(iprt)                          !ala
        else                                              !ala
          iodur(nj)=0                                     !ala
        endif                                             !ala

        if(iqj(nj).ne.0)then
          njc=nci(1,iprt)
          if(njc.ne.0)then
            ncj(1,nj)=njc
            iqc=iqj(njc)
            if(iqc.ne.0)then
              ncj(1,njc)=nj
            else
              if(iqj(nj).gt.0)then
                ncj(2,njc)=nj
              else
                ncj(1,njc)=nj
              endif
            endif
          else
            nci(1,iprt)=nj
          endif
        else !----gluon----
          do m=1,2
            if(jort.eq.1)then
              m1=m
            else
              m1=3-m
            endif
            njc=nci(m1,iprt)
            if(njc.ne.0)then
              ncj(m,nj)=njc
              iqc=iqj(njc)
              if(iqc.ne.0)then
                ncj(1,njc)=nj
              else
                ncj(3-m,njc)=nj
              endif
            else
              nci(m1,iprt)=nj
            endif
          enddo
        endif

      endif

      if(ish.ge.4)then !~~~~~~~~~~~~~~~~print~~~~~~~~~~~~~~~~~~~~~~
        write (ifmt,'(a,5(2i4,3x))')'nr id  idau1/2  nci  ncj+'
     .  ,iprt,icp,idau1,idau2,(nci(k,iprt),k=1,2),(ncj(k,nj),k=1,2)
      endif !~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

2     continue
      if(iprt.ne.nfprt)then
        icp=idprt(iprt)
        ipar=iorprt(iprt)
        idau1=idaprt(1,ipar)
        idau2=idaprt(2,ipar)
        if(ish.ge.6)then
          write (ifch,*)'2-iprt,icp,idau1,idau2,ipar',
     *    iprt,icp,idau1,idau2,ipar,nci(1,iprt)
          if(icp.eq.9)write (ifch,*)nci(2,iprt)
        endif

        if(idau1.eq.iprt)then
          if(iprt.eq.0.or.idau2.eq.0.or.ipar.eq.0)then
            write (ifmt,*)'Error 1-iprt,icp,idau1,idau2,ipar',
     *           iprt,icp,idau1,idau2,ipar,nci(1,iprt)
            call utstop("Null array index in psreti (1) !&") 
          endif
          if(icp.eq.9)then                   !g -> gg
            nci(1,ipar)=nci(1,iprt)
            nci(1,idau2)=nci(2,iprt)
          else
            icpar=idprt(ipar)
            if(icpar.eq.9)then               !g -> qq~
              nci(jort,ipar)=nci(1,iprt)
            else                             !q -> qg
              if(icp*(3-2*jort).gt.0)then
                nci(2,idau2)=nci(1,iprt)
              else
                nci(1,idau2)=nci(1,iprt)
              endif
            endif
          endif
          iprt=idau2
          goto 1

        else
          if(iprt.eq.0.or.idau1.eq.0.or.ipar.eq.0)then
            write (ifmt,*)'Error 2-iprt,icp,idau1,idau2,ipar',
     *           iprt,icp,idau1,idau2,ipar,nci(1,iprt)
            call utstop("Null array index in psreti (2) !&") 
          endif
          if(icp.eq.9)then
            icpar=idprt(ipar)
            if(icpar.eq.9)then                !g -> gg
              nci(2,ipar)=nci(2,iprt)
              nci(2,idau1)=nci(1,iprt)
            else                              !q -> qg
              if(icpar*(3-2*jort).gt.0)then
                nci(1,ipar)=nci(1,iprt)
                nci(1,idau1)=nci(2,iprt)
              else
                nci(1,ipar)=nci(2,iprt)
                nci(1,idau1)=nci(1,iprt)
              endif
            endif
          else
            nci(3-jort,ipar)=nci(1,iprt)
          endif
          iprt=ipar
          goto 2
        endif
      else
        if(ish.ge.6)write (ifch,*)'3-iprt,nci',iprt,nci(1,iprt)
c     *                             ,pprtx(1,iprt),pprtx(2,iprt)
        ncr(1)=nci(1,nfprt)
        if(idprt(nfprt).eq.9)ncr(2)=nci(2,nfprt)
      endif

      if(ish.ge.4)then !~~~~~~~~~~~~~~print~~~~~~~~~~~
        write(ifmt,'(a,2i4,i8)')'psreti final' 
        do ij=1,nj
        call idtrafo_KW_SO( iKW , iqj(ij) , -1 ) 
        write(ifmt,'(i4,3x,3i4,5f10.2)')iKW
     .  ,ncj(1,ij),ij,ncj(2,ij),(eqj(j,ij),j=1,4)
     .  ,eqj(1,ij)**2-eqj(2,ij)**2-eqj(3,ij)**2-eqj(4,ij)**2
        call eqjcheck(9,ij)
        enddo
      endif !~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

      do ij=1,nj
        i1=ncj(1,ij) 
        if(i1.ne.0)then !ij connected to i1
          if(ncj(1,i1).ne.ij.and.ncj(2,i1).ne.ij)then
            write(ifmt,'(a,i3,a,i3,a,a,i3,a)')'ERROR psjeti;'
     .      ,ij,' connected to',i1,' not inverse;'
     .      ,' last flavor:',icp,'  (see check)'
            write(ifch,'(a,2i4,i8)')'psreti: flavor connection mismatch' 
            do ijx=1,nj
              call idtrafo_KW_SO( iKW , iqj(ijx) , -1 ) 
              write(ifch,'(i4,3x,3i4,5f10.2)')iKW
     .        ,ncj(1,ijx),ijx,ncj(2,ijx),(eqj(j,ijx),j=1,4)
     .        ,eqj(1,ijx)**2-eqj(2,ijx)**2-eqj(3,ijx)**2-eqj(4,ijx)**2
            enddo
            stop
            !iret=1
            !return
          endif
        endif  
        i2=ncj(2,ij) 
        if(i2.ne.0)then !ij connected to i2
          if(ncj(1,i2).ne.ij.and.ncj(2,i2).ne.ij)then
            write(ifmt,'(a,i3,a,i3,a,a,i3,a)')'ERROR psjeti'
     .      ,ij,' connected to',i2,' not inverse'
     .      ,' last flavor:',icp,'  (see check)'
            write(ifch,'(a,2i4,i8)')'psreti: flavor connection mismatch' 
            do ijx=1,nj
              call idtrafo_KW_SO( iKW , iqj(ijx) , -1 ) 
              write(ifch,'(i4,3x,3i4,5f10.2)')iKW
     .        ,ncj(1,ijx),ijx,ncj(2,ijx),(eqj(j,ijx),j=1,4)
     .        ,eqj(1,ijx)**2-eqj(2,ijx)**2-eqj(3,ijx)**2-eqj(4,ijx)**2
            enddo
            stop
            !iret=1
            !return
          endif
        endif  
      enddo

      call utprix('psreti',ish,ishini,4)

      return
      end

      subroutine putInfoPsreti(ipto,jj,ibo,nl1,nr1,nl2,nr2)
      common/cpsreti/iptox,jjx,ibox,nl1x,nr1x,nl2x,nr2x
      iptox=ipto
      jjx=jj
      ibox=ibo
      nl1x=nl1
      nr1x=nr1
      nl2x=nl2
      nr2x=nr2
      end  

      subroutine getInfoPsreti(ipto,jj,ibo,nl1,nr1,nl2,nr2)
      common/cpsreti/iptox,jjx,ibox,nl1x,nr1x,nl2x,nr2x
      ipto=iptox
      jj=jjx
      ibo=ibox
      nl1=nl1x
      nr1=nr1x
      nl2=nl2x
      nr2=nr2x
      end

c------------------------------------------------------------------------
      subroutine redefineMassPartons
c-----------------------------------------------------------------------
c     This subroutine has been created to solve a problem with heavy quark mass. 
c     The subroutine gakfra (fra.f) needs pptl(5,)=0. This subroutine use the real kinematic, i.e. E and P as:
c     m**2=E**2-p**2, but we put by hand pptl(5,)=0 instead of pptl(5,)=0.
c     Later in bas.f, if pptl(5,i) diff of m_i then the kinematic is redefined in order to get the relation :
c     pptl(5,i)**2 = m_i**2 = E_i**2-p_i**2 . We don't want that since the real kinematics has been provided.
c     
c     In this subroutine we simply do pptl(5,i)=m_i
c-----------------------------------------------------------------------
#include "aaa.h"

      call idmass(4,amc)
      call idmass(5,amb)
      do i=1,nptl
        if(iabs(idptl(i)).eq.4.and.istptl(i).eq.21)then
          pptl(5,i)=amc
          p5=pptl(4,i)**2-pptl(3,i)**2-pptl(2,i)**2-pptl(1,i)**2
          diff=abs(p5-amc**2)
          if(diff.gt.0.5)write(*,*)'wrong kinematic',i,amc,p5,
     $   pptl(4,i),pptl(3,i),pptl(2,i),pptl(1,i),istptl(i),idptl(i),nptl
          if(diff.gt.0.5)stop 'need to coment lines for heavy quark mass
     $                                                        in pspawr'
        elseif(iabs(idptl(i)).eq.5.and.istptl(i).eq.21)then
          pptl(5,i)=amb
          p5=pptl(4,i)**2-pptl(3,i)**2-pptl(2,i)**2-pptl(1,i)**2
          diff=abs(p5-amb**2)
          if(diff.gt.0.5)write(*,*)'wrong kinematic',i,amb,p5,
     $    pptl(4,i),pptl(3,i),pptl(2,i),pptl(1,i),istptl(i),idptl(i)
        endif
      enddo

      end

c-----------------------------------------------------------------------
      integer function igetRNDindex4(w1,w2,w3) 
c-----------------------------------------------------------------------
      real xkk(4)
      !if(w1+w2+w3.gt.1)then
      !  write(ifmt,'(a,4f8.3)')'WARNING w1 w2 w3 sum ='
      !.  ,w1,w2,w3,w1+w2+w3
      !  if(w1+w2+w3.gt.1.05)stop'ERROR 27062019' 
      !endif
      xkk(1)=w1
      xkk(2)=w2
      xkk(3)=w3
      xkk(4)=1-w1-w2-w3
      r=rangen()
      w=0.
      i=0
      do while (w.lt.r)
        i=i+1
        if(i.gt.4)then
          call getMonitorFileIndex(ifmt)
          write(ifmt,'(a)')'WARNING igetRNDindex4 i > 4'
          igetRNDindex4=4
          return 
        endif
        w=w+xkk(i)
      enddo          
      igetRNDindex4=i 
      end

c-----------------------------------------------------------------------
      integer function igetRNDindex(m,xkk) !KW1808e
c-----------------------------------------------------------------------
      real xkk(m)
      integer ienvi
      common /cienvi/ ienvi
      sum=0
      do i=1,m
        sum=sum+xkk(i)
      enddo
      !~~~~~~~~~~check~~~~~~~~~
      if(sum.le.0.)then
        write(ifmt,'(a,$)')'ERROR igetRNDindex Empty array;  '
        write(ifmt,*)' ienvi =',ienvi
        stop
      endif
      !~~~~~~~~~~~~~~~~~~~~~~~~
      r=rangen()*sum
      w=0.
      i=0
      do while (w.lt.r)
        i=i+1
        w=w+xkk(i)
      enddo          
      igetRNDindex=i 
      end

c-----------------------------------------------------------------------
      subroutine initPartonList(nj1,nj2)
c-----------------------------------------------------------------------
      parameter (mjstr=20000)
      common /psar29/eqj(4,mjstr),iqj(mjstr),ncj(2,mjstr),ioj(mjstr),nj
      common /psar30/ iorj(mjstr),ityj(mjstr),bxj(6,mjstr),q2j(mjstr)
      common /psar31/ ktyj(mjstr)
      if(nj2.lt.nj1)return
      do i=nj1,nj2
      ncj(1,i)=0
      ncj(2,i)=0
      eqj(1,i)=0
      eqj(2,i)=0
      eqj(3,i)=0
      eqj(4,i)=0
      iqj(i)=0
      ioj(i)=0
      iorj(i)=0
      ityj(i)=0
      bxj(1,i)=0
      bxj(2,i)=0
      bxj(3,i)=0
      bxj(4,i)=0
      bxj(5,i)=0
      bxj(6,i)=0
      q2j(i)=0
      ktyj(i)=0
      enddo
      end

c------------------------------------------------------------------------
      subroutine psabor(si,qq0,qt2,mIB,nIB,mOB,nOB,ncc,jdis,iptl
     .            ,coordo,jt,ey,s0xi,c0xi,s0i,c0i,ipomtype,qqupp)
c------------------------------------------------------------------------
c
c realizes Born process
c
c   si ......... Mandestam s
c   qq0 ........ maximum Q2 at Pomeron ends
c   qt2 ......... pt^2
c   mIB,nIB,mOB,nOB ..... incoming/outgoing parton ids (0=g, -5, ...,5 for quarks
c   ncc(2,k) ... indices of color conn partons of incoming parton k (1 and 2) (see epos4 manual)
c   jdis ....... type of process
c                  -1 ... saturated Pom, calculate born AT qq0
c                   0 ... hadronic process in pp
c                   1 .... resolved photon process
c   iptl ....... index of Pomeron
c   coordi(6) .. x,y,z,t,tiv1,tiv2 (not used)
c   jt .......... process number
c   ey,s0xi,c0xi,s0i,c0i ... boost and rotation parameters
c------------------------------------------------------------------------
      parameter (lsub=5)
      double precision psutz,si,qt2,qq,z,wp3,wm3
     &,q2m1,q2m2,en1,en2,qqborn
      dimension ep3(4),iIB(2),iOB(2),ncc(2,2),ncr(2),coordo(6)
      dimension ey(3)
      parameter (mjstr=20000)
      common /psar29/eqj(4,mjstr),iqj(mjstr),ncj(2,mjstr),ioj(mjstr),nj
      parameter (ntim=1000)
      common/cprt/pprt(5,ntim),q2prt(ntim),idaprt(2,ntim),idprt(ntim)
     &,iorprt(ntim),jorprt(ntim),nprtj
      common/cprtx/nprtjx,pprtx(5,2)
      common/prtdur/idprtx(2) 
      common/cqqmx/qqmx1,qqmx2
      common/cchtest/ichtest
#include "aaa.h"
#include "sem.h"
      data kkone /0/ kktwo /0/ kkthr/0/ kkzer /0/ 
      save kkone,kktwo,kkthr,kkzer
      common/chkbor/ichkbor
      integer ienvi
      common /cienvi/ ienvi
      data kerr3 /0/ kerr4 /0/
      save kerr3, kerr4 

      ienvi=101  !psabor

      !-----------------------------------------------------------------
      !  color connection information:
      !-----------------------------------------------------------------
      !   ncj(1:2,k) - indices of color connected partons for TL partons 
      !              in /psar29/    nj - number of partons (max value of k)
      !   ncc(1:2,j) indices of color connected partons for incoming 
      !         space-like partons (j=1 and 2) referring to the /psar29/ list
      !-----------------------------------------------------------------
      ! ncj in /psar29/ (concerning the time-like partons) 
      ! + ncc (concerning the space-like incoming partons)  
      ! give complete information about color connections
      !-----------------------------------------------------------------
      ! The picture will be completed after the timelike emissions, using the 
      !  array  ncr(1:2) - color connection of one parton, the one creating 
      !  the time-like cascade
      ! We call psreti twice, for two time-like partons (jets), with input "ncr".
      !  Connections already decided after chosing the Born process, will be 
      !  put into "ncj". Connections still "pending" will be put into "ncr"  
      !-----------------------------------------------------------------

      call utpri('psabor',ish,ishini,5)

c inBorn types

      if(mIB.eq.0.and.nIB.eq.0)then !gg
        ibo   = 31
      elseif(mIB.eq.0.and.abs(nIB).le.3
     .   .or.nIB.eq.0.and.abs(mIB).le.3)then !gq
        ibo   = 32
      elseif(abs(mIB).le.3.and.abs(nIB).le.3)then !qq
        ibo   = 33
      elseif(mIB.eq.0.and.abs(nIB).gt.3
     .   .or.nIB.eq.0.and.abs(mIB).gt.3)then !gQ
        ibo   = 34
      elseif(abs(mIB).le.3.and.abs(nIB).gt.3
     .   .or.abs(nIB).le.3.and.abs(mIB).gt.3)then !qQ
        ibo   = 35
      elseif(abs(mIB).gt.3.and.abs(nIB).gt.3)then !QQ
        ibo   = 36
      else 
        stop'ERROR 08062020'
      endif 

c initial stuff

      !write(ifmt,*)'color connections (start psabor)'
      !do ij=1,nj
      !write(ifmt,*)'  ',iqj(ij),ij,(ncj(j,ij),j=1,2)
      !enddo
      !print*,'    ',mIB,   ncc(1,1),ncc(2,1)
      !print*,'    ',nIB,   ncc(1,2),ncc(2,2)

      j2=mIB  ! id parton 1  (0,+-1,+-2,+-3,+-4,+-5)
      l2=nIB  ! id parton 2  (0,+-1,+-2,+-3,+-4,+-5)
      klas=igetKlas( mIB,nIB,  mOB,nOB )
      call HardScale2(klas,qq,qt2,2) !get qq
      qqborn=qq
      if(jdis.lt.0)qqborn=dble(qq0) 

      if(rangen().lt..5.and.jdis.gt.0)then  
        ii=2    !targ side              
      else                    
        ii=1    !proj side              
      endif
      iIB(1)=mIB
      iIB(2)=nIB
      iOB(1)=mOB
      iOB(2)=nOB

c determine transv momentum

      call pscs(bcos,bsin)          !cos and sin of the polar angle
      qt=sngl(sqrt(qt2))            !p_t

c determination of qt2 in the lab frame for photon

      if(iIB(ii).eq.0 .and.iIB(3-ii).ne.0)then ! ---- g q -> q +gamma ----
        z=psutz(si,qt2,qt2) !bg
        if(ii.eq.2)z=1d0-z
        wp3=z*sqrt(si)
        wm3=(qt2)/wp3
        ep3(1)=.5*sngl(wp3+wm3)               !4-momentum for out born parton 1
        ep3(2)=.5*sngl(wp3-wm3)
        ep3(3)=-qt*bcos
        ep3(4)=-qt*bsin         !bg minus because gamma is the second output
        call psrotat(ep3,s0xi,c0xi,s0i,c0i)     !rotation in LF 
        call pstrans(ep3,ey,1)                !back to LAB
        call psmirror(ep3)
        ptgam=sqrt(ep3(3)**2+ep3(4)**2)
      elseif(iIB(3-ii).eq.0 .and.iIB(ii).ne.0)then ! ---- q g -> q + gamma ----
        z=psutz(si,qt2,qt2) !bg
        if(ii.eq.1)z=1d0-z
        wp3=z*sqrt(si)
        wm3=(qt2)/wp3
        ep3(1)=.5*sngl(wp3+wm3)               !4-momentum for out born parton 1
        ep3(2)=.5*sngl(wp3-wm3)
        ep3(3)=-qt*bcos
        ep3(4)=-qt*bsin         !bg minus because gamma is the second output
        call psrotat(ep3,s0xi,c0xi,s0i,c0i)     !rotation in LF 
        call pstrans(ep3,ey,1)                !back to LAB
        call psmirror(ep3)
        ptgam=sqrt(ep3(3)**2+ep3(4)**2)
      elseif(iIB(ii).eq.-iIB(3-ii).and.iIB(ii).ne.0)then ! ---- q qbar -> g + gamma ----
        z=psutz(si,qt2,qt2) !bg
        if(ii.eq.1)z=1d0-z
        wp3=z*sqrt(si)
        wm3=(qt2)/wp3
        ep3(1)=.5*sngl(wp3+wm3)               !4-momentum for out born parton 1
        ep3(2)=.5*sngl(wp3-wm3)
        ep3(3)=-qt*bcos
        ep3(4)=-qt*bsin         !bg minus because gamma is the second output
        call psrotat(ep3,s0xi,c0xi,s0i,c0i)     !rotation in LF 
        call pstrans(ep3,ey,1)                !back to LAB
        call psmirror(ep3)
        ptgam=sqrt(ep3(3)**2+ep3(4)**2)
      endif

c get process number jt and update color configuration array ncj and determine jq, jq2

             !----------------------!  *** ncc(A,B) array *** 
             !         A=1          !     
             ! B=1__        ___B=2  !    A refers to projectile   
             !      |      |        !         or target
             !      v      ^        !      
             !      |      |        !   B refers to color or
             !      --------        !      anticolor side 
             !      | Born |        !   
             !      --------        !   ncc(A,B) refers to the
             !      |      |        !  TL parton index at the end of
             !      v      ^        !     the corresponding leg  
             !      |      |        ! 
             ! B=2__|      |__B=1   !  *** ncr(1) -> ncr(2) ***
             !         A=2          !  defines flow of emissions
             !----------------------!  in case of gq starting from q side

      ncr(2)=0

      if(iIB(ii).eq.0.and.iIB(3-ii).eq.0)then  !--------------------------g+g------------------------------

        jq=int(1.5+rangen()) !emission on color or anticolor side
        ncr(1)=ncc(jq,ii)
        if(iOB(ii).eq.0.and.iOB(3-ii).eq.0)then             !------------g+g->g+g ---------jt=1--------
          if(rangen().lt..5)then
            jt=1       !anticolor-color annihilation 
            ncr(2)=0   !connect color lines opposte to emission
            nj1=ncc(3-jq,ii)
            nj2=ncc(jq,3-ii)      
            if(nj1*nj2.eq.0)       !-------------!    !-------------!    
     .      stop'ERROR 24072019a'  !   ___| |    !    !   | |___    !
            if(iqj(nj1).ne.0)then  !   ___  |    !    !   |  ___    !
              ncj(1,nj1)=nj2       !      | |    !    !   | |       !
            else                   !   ___| |    ! or !   | |___    !
              ncj(jq,nj1)=nj2      !   ___  |    !    !   |  ___    !
            endif                  !      | |    !    !   | |       !
            if(iqj(nj2).ne.0)then  !-------------!    !-------------!
              ncj(1,nj2)=nj1
            else
              ncj(3-jq,nj2)=nj1
            endif
          else                     !-------------!   !----------g+g->g+g --------jt=2--------
                                   !   __| |__   !
            jt=2                   !   __   __   !
            ncr(2)=ncc(3-jq,3-ii)  !     | |     !
                                   !-------------!  
          endif                   
          !print*,'BORN NR',jt,jq,nj1,nj2
          !.,ncj(1,nj1),ncj(2,nj1),ncj(1,nj2),ncj(2,nj2)
        else                                       !-----------g+g->q+q~ ---------jt=9--------  
          jt=9                        
          if(iOB(ii).ne.-iOB(3-ii))stop'ERROR 24072019b'
          nj1=ncc(3-jq,ii)
          nj2=ncc(jq,3-ii)
          if(nj1*nj2.eq.0)stop'ERROR 24072019c'
          if(iqj(nj1).ne.0)then
            ncj(1,nj1)=nj2
          else
            ncj(jq,nj1)=nj2
          endif
          if(iqj(nj2).ne.0)then
            ncj(1,nj2)=nj1
          else
            ncj(3-jq,nj2)=nj1
          endif
        endif
        if(iOB(ii)*(3-2*jq).lt.0)call switchComponents(iOB)

      elseif(iIB(ii)*iIB(3-ii).eq.0)then       !--------------------------- q + g -----------------------------

        if(iIB(ii)+iIB(3-ii).gt.0)then
          jq=1                             !q
        else
          jq=2                             !q~
        endif
        if(iOB(ii)*iOB(3-ii).eq.0.and.abs(iIB(ii)+iIB(3-ii)).le.5)then  !q+g->q+g
          if(rangen().lt..5)then !anticolor-color annihilation
            if(iIB(ii).eq.0)then                !-------------g+q->g+q  ----------jt=3----------
              jt=3                            
              ncr(1)=ncc(jq,ii)
              nj1=ncc(3-jq,ii)
              nj2=ncc(1,3-ii)
              if(nj1*nj2.eq.0)stop'ERROR 24072019d'
              if(iqj(nj1).ne.0)then
                ncj(1,nj1)=nj2
              else
                ncj(jq,nj1)=nj2
              endif
              if(iqj(nj2).ne.0)then
                ncj(1,nj2)=nj1
              else
                ncj(3-jq,nj2)=nj1
              endif
              if(iOB(ii).ne.0)call switchComponents(iOB)
            else
              jt=4                              !--------------q+g->q+g ---------jt=4----------
              ncr(1)=0
              nj1=ncc(1,ii)          !-------------!
              nj2=ncc(3-jq,3-ii)     !        |q   !
              if(nj1*nj2.eq.0)       !  q___  |    !
     .        stop'ERROR 24072019e'  !      | |    !  create first q w/o connections 
              if(iqj(nj1).ne.0)then  !   ___| |    !  do connections in next step after g emission
                ncj(1,nj1)=nj2       !   ___  |    !
              else                   !      |g|    !
                ncj(3-jq,nj1)=nj2    !-------------!
              endif
              if(iqj(nj2).ne.0)then
                ncj(1,nj2)=nj1
              else
                ncj(jq,nj2)=nj1
              endif
              if(iOB(ii).eq.0)call switchComponents(iOB)
            endif
          else   !color transfer     !-------------!
            if(iIB(ii).eq.0)then     !     |g|     ! !---------g+q->g+q ----------jt=5----------
              jt=5                   !   __| |__q  !
              ncr(2)=ncc(3-jq,ii)    !   __        !
              ncr(1)=ncc(1,3-ii)     !     |       !
              !start from below      !     |q      ! 
              if(iOB(ii).ne.0)call   !-------------!
     .        switchComponents(iOB)
            else                                     !---------q+g->q+g ---------jt=6----------   
               
              jt=6  
              ncr(1)=ncc(jq,3-ii)
              if(iOB(ii).eq.0)call switchComponents(iOB)

            endif
          endif
        else        ! q+g->q+gamma (+-color annihilation)
          if(iIB(ii).eq.0)then
            jt=11                            !-------------g+q->q+gamma ---------jt=11---------  
            ncr(1)=ncc(jq,ii)
            nj1=ncc(3-jq,ii)
            nj2=ncc(1,3-ii)
            if(nj1*nj2.eq.0)stop'ERROR 24072019f'
            if(iqj(nj1).ne.0)then
              ncj(1,nj1)=nj2
            else
              ncj(jq,nj1)=nj2
            endif
            if(iqj(nj2).ne.0)then
              ncj(1,nj2)=nj1
            else
              ncj(3-jq,nj2)=nj1
            endif
            if(iOB(ii).eq.10)call switchComponents(iOB) !photon second
          else
            jt=12                            !-------------q+g->q+gamma ----------jt=12---------     
            ncr(1)=ncc(jq,3-ii)            
            nj1=ncc(1,ii)
            nj2=ncc(3-jq,3-ii)
            if(nj1*nj2.eq.0)stop'ERROR 24072019g'
            if(iqj(nj1).ne.0)then
              ncj(1,nj1)=nj2
            else
              ncj(3-jq,nj1)=nj2
            endif
            if(iqj(nj2).ne.0)then
              ncj(1,nj2)=nj1
            else
              ncj(jq,nj2)=nj1
            endif
            if(iOB(ii).eq.10)call switchComponents(iOB) !photon second
          endif
        endif
                                                   
      elseif(iIB(ii)*iIB(3-ii).gt.0)then      ! q + q' both quarks (q) or both antiquarks (a)

                                      !-----------qq'->qq' (both q or a) -------jt=7--------   
        jt=7                    !-------------!      
        if(iIB(ii).gt.0)then    !     q       !                    
          jq=1                  !     |       !                  
        else                    !     |  __q  !                 
          jq=2                  ! q'__|g|     !
        endif                   !       |     !
        ncr(1)=ncc(1,3-ii)      !       |     !
                                !       q'    !
                                !-------------!

      else    ! q + q'  one is quark the other antiquark 
              !    (not necessarily same flavor)

        if(iIB(ii).gt.0)then    
          jq=1
        else
          jq=2
        endif

        if(iOB(ii).ne.0.and.iOB(ii).ne.10)then    !----------- qo->q'o' ---------jt=8-------- 
          jt=8                 
          ncr(1)=0                      !--------------!    q means some flavor (>0,<0)
          nj1=ncc(1,ii)                 !  q'    q .   ! o means some flavor of opposite sign
          nj2=ncc(1,3-ii)               !    \   | .   ! 
          if(nj1*nj2.eq.0)              !     \__| .   !
     .     stop'ERROR 24072019g'        !      __g .   !  g is s or t 
          if(iqj(nj1).ne.0)then         !     /  | .   !     channel gluon
            ncj(1,nj1)=nj2              !  o'/   | .   !
          else                          !        o     ! o may be -q or not
            ncj(3-jq,nj1)=nj2           !--------------!
          endif
          if(iqj(nj2).ne.0)then
            ncj(1,nj2)=nj1
          else
            ncj(jq,nj2)=nj1
          endif
          if(iOB(ii)*iIB(ii).lt.0)call switchComponents(iOB) !photon second
        elseif(iOB(ii).eq.0.and.iOB(3-ii).eq.0)then
          jt=10                                   !------------qq~->gg  ---------jt=10-------- 
          ncr(1)=ncc(1,ii)                   !(color transfer)
          ncr(2)=0
        elseif(iOB(ii)*iOB(3-ii).eq.0.and.iOB(ii)+iOB(3-ii).eq.10)then
          jt=13                             !--------------qq~->g+gamma ---------jt=13-------- 
          ncr(1)=ncc(1,ii)
          ncr(2)=ncc(1,3-ii)
          if(iOB(ii).eq.10)call switchComponents(iOB) !photon second
        elseif(iOB(ii).eq.10.and.iOB(3-ii).eq.10)then
          jt=14                          !-------------qq~->gamma+gamma ---------jt=14-------- 
          nj1=ncc(1,ii)
          nj2=ncc(1,3-ii)
          if(nj1*nj2.eq.0)stop'ERROR 24072019h'
          if(iqj(nj1).ne.0)then
            ncj(jq,nj1)=nj2
          else
            ncj(3-jq,nj1)=nj2
          endif
          if(iqj(nj2).ne.0)then
            ncj(3-jq,nj2)=nj1
          else
            ncj(jq,nj2)=nj1
          endif
        else
          write(ifmt,*)'ERROR PSABOR'
     .    ,iIB(ii),iIB(3-ii),iOB(ii),iOB(3-ii)
          stop'no process index jt found'
        endif

      endif

      if(ish.ge.4)write(ifmt,'(80a1)')('=',ite=1,80)
      if(ish.ge.3)write(ifmt,'(a,i3)')'BORN NR',jt

      if(jt.ne.8.and.jt.ne.9)then
        jq2=jq
      else
        jq2=3-jq
      endif

c compute squared masses

      q2m1=0.d0      
      q2m2=0.d0                 
      if(iabs(iOB(ii)).eq.4)then   
        q2m1=dble(qcmass**2)        
      elseif(iabs(iOB(ii)).eq.5)then       
        q2m1=dble(qbmass**2)          
      endif                          
      if(iabs(iOB(3-ii)).eq.4)then      
        q2m2=dble(qcmass**2)         
      elseif(iabs(iOB(3-ii)).eq.5)then      
        q2m2=dble(qbmass**2)         
      endif                           

      !write(ifmt,*)'color connections; process',jt, ncr
      !do ij=1,nj
      !write(ifmt,*)'      ',iqj(ij),ij,(ncj(j,ij),j=1,2)
      !enddo

c momentum of the first parton, compute ROTATION parameters

      z=psutz(si,qt2+q2m1,qt2+q2m2) !bg
      if((jt.eq.11.and.ii.eq.1).or.(jt.eq.12.and.ii.eq.2)
     $ .or.(jt.eq.13.and.ii.eq.2))z=1d0-z
      wp3=z*sqrt(si)
      wm3=(qt2+q2m1)/wp3 !bg only q2m1 because wm3 is for iOB(ii)

      ep3(1)=.5*sngl(wp3+wm3)               !4-momentum for out born parton 1
      ep3(2)=.5*sngl(wp3-wm3)
      ep3(3)=qt*bcos
      ep3(4)=qt*bsin
      !if(ep3(1).ge.3.)write(ifmt,'(/a,i1,a,4f10.2)')
      !.'psabor: p_OB',1,' =',ep3(3),ep3(4),ep3(2),ep3(1)
      call psdefrot(ep3,s0xh,c0xh,s0h,c0h)   !compute parameters for spatial rotation to z-axis

c time-like cascade
      zeta=2.                        !2=back-to-back emission (angle=pi)
      if(iOB(ii).eq.0)then
        iq2ini1=9
        jo1=iq2ini1
        if(zeta.gt.zetacut)jo1=-jo1
      else
        iq2ini1=iOB(ii)
        jo1=iq2ini1
      endif
      if(iOB(3-ii).eq.0)then
        iq2ini2=9
        jo2=iq2ini2
        if(zeta.gt.zetacut)jo2=-jo2
      else
        iq2ini2=iOB(3-ii)
        jo2=iq2ini2
      endif

      if(q2m1.gt.0.5d0.and.abs(iOB(ii)).ne.abs(iOB(3-ii))) then
        en1=(si+q2m1)/(2.d0*sqrt(si))  
        en2=(si-q2m1)/(2.d0*sqrt(si))
      else
        en1=sqrt(si)/2
        en2=en1
      endif

      !Q2_max for TL cascade
      !---nexus method +mass terms---
      if(jt.le.10)then
        qq1=qt2+q2m1
        qq2=qt2+q2m2
      elseif(jt.le.13)then 
        qq1=qt2+q2m1
        qq2=0
      else
        qq1=0
        qq2=0
      endif
      qq1=qq1*facq2tim
      qq2=qq2*facq2tim
      dummy=qqupp
      !---modif qq1,2 sat pom---
      if(ipomtype.lt.0)then !fine tuning sat
      !  call idmass(4,am4)
      !  qq1=8 !4.*q2fin+2.*am4**2+0.1 !min value allowing ccbar split
      !  qq2=qq1                    
      !  if(mIB.eq.0.and.nIB.eq.0.and.mOB.eq.0.and.nOB.eq.0)then !ggggsat
      !    qq1=9. !qqupp
      !    qq2=9. !qqupp
      !  endif
      endif
      !---benjamin method---  
      !if(jt.le.13)then
      !  !qq1,qq2 = min(pt^2,E^2-pt^2) to avoid kinematics incompatibilities for pt^2 > E^2/2
      !  !since we have Q^2 < E^2 - pt^2, so with Q^2 = pt^2: pt^2 < E^2/2 
      !  qq1=sngl((en1-sqrt(qt2))*(en1+sqrt(qt2)))
      !  if(qq1.gt.sngl(qt2+q2m1))qq1=sngl(qt2+q2m1)
      !  if(jt.le.10)then
      !    if(q2m1.gt.0.5d0.and.abs(iOB(ii)).ne.abs(iOB(3-ii))) then
      !      qq2=sngl((en2-sqrt(qt2))*(en2+sqrt(qt2)))
      !    else
      !      qq2=qq1
      !    endif
      !    if(qq2.gt.sngl(qt2+q2m2))qq2=sngl(qt2+q2m2)
      !  elseif(jt.le.13)then
      !    qq2=0.
      !  endif
      !else
      !  qq1=0.
      !  qq2=0.
      !endif

      !KW2006  Things are not symmetric, I leave it for the moment ... ??????????????????????   

      !bg pz**2=E**2-pt**2-m**2 > 0 if parton on-shell    
      if(en1**2.lt.(qt2+q2m1)-0.0001
     . .or.qq1.lt.sngl(q2m1).or.qq2.lt.sngl(q2m2) )then
        if(en1**2.lt.(qt2+q2m1)-0.0001
     .    .or.qq1.lt.sngl(q2m1)-0.001
     .    .or.qq2.lt.sngl(q2m2)-0.001) then
          write(ifmt,'(a,i4,3x,4i3,1x,2f10.4,4f9.4,i4)')
     .    'WARNING psabor bef timsh2',ipomtype,j2,l2,iOB(ii),iOB(3-ii) 
     .    ,en1**2-qt2-q2m1,en1**2,qt2,sngl(q2m1),qq1-q2m1,qq1,jt
        endif
        if(qq1.le.sngl(q2m1))qq1=sngl(q2m1)+1e-5
        if(qq2.le.sngl(q2m2))qq2=sngl(q2m2)+1e-5
        if(qq1.le.0.)qq1=1e-5
        if(qq2.le.0.)qq2=1e-5
      endif

      q2finsave=q2fin
      if(sngl(en1).gt.getJetEmin())then
        !determine pt
        pt=utJetPt(sngl(en1),ey,s0xi,c0xi,s0i,c0i,s0xh,c0xh,s0h,c0h)
        if(pt.gt.ptrmin)then
          q2fin=1e30
          call incrementJetNdijet()
          !print*,'TEST=Jet',sngl(en1),getJetEmin(),pt,ptrmin
        endif
      endif
      
      !===========================================================
      call timsh2(qq1,qq2,sngl(en1),sngl(en2),iq2ini1,iq2ini2,jo1,jo2
     . ,ey,s0xi,c0xi,s0i,c0i,s0xh,c0xh,s0h,c0h,ipomtype,mIB,nIB,mOB,nOB)  !final state cascade
      !============================================================

      q2fin=q2finsave

      qqmx1=qq1
      qqmx2=qq2

c color connections

      if(jt.le.10)then    !NO gammas involved

        if(ish.ge.6)write(ifch,*)'jt,jq,ncr:',jt,jq,ncr
     .                                       ,ey,s0xh,c0xh,s0h,c0h
       !===============================================
        call putInfoPsreti(0,11,ibo,ncc(1,1),ncc(2,1),ncc(1,2),ncc(2,2))
        call psreti(ncr,jq,1
     .          ,ey,s0xi,c0xi,s0i,c0i,s0xh,c0xh,s0h,c0h,iret) !color conn. reconstruction
       !===============================================

        if(iret.ne.0)then
          write(ifmt,'(a)')'ERROR psahot Nonzero psreti-iret ;'
          stop
        endif

        !write(ifmt,*)'color connections (after psreti)', ncr
        !do ij=1,nj
        !write(ifmt,*)'      ',iqj(ij),ij,(ncj(j,ij),j=1,2)
        !enddo

        if(jt.eq.1)then
          ncr(1)=ncr(2)
          ncr(2)=ncc(3-jq,3-ii)
        elseif(jt.eq.2)then
          ncr(2)=ncc(3-jq,ii)
          ncr(1)=ncc(jq,3-ii)
        elseif(jt.eq.3)then
          ncr(1)=ncr(2)
        elseif(jt.eq.4)then
          ncr(2)=ncr(1)
          ncr(1)=ncc(jq,3-ii)
        elseif(jt.eq.5)then
          ncr(1)=ncc(jq,ii)
          ncr(2)=0
        elseif(jt.eq.6)then
          ncr(2)=ncc(3-jq,3-ii)
          ncr(1)=ncc(1,ii)
        elseif(jt.eq.7)then
          ncr(1)=ncc(1,ii)
        elseif(jt.eq.9)then
          ncr(1)=ncc(3-jq,3-ii)
        elseif(jt.eq.10)then
          ncr(1)=ncr(2)
          ncr(2)=ncc(1,3-ii)
        endif
        if(ish.ge.6)write(ifch,*)'jt,jq2,ncr:',jt,jq2,ncr

       !================================================
        call putInfoPsreti(0,12,ibo,ncc(1,1),ncc(2,1),ncc(1,2),ncc(2,2))
        call psreti(ncr,jq2,2
     .          ,ey,s0xi,c0xi,s0i,c0i,s0xh,c0xh,s0h,c0h,iret) !color conn. reconstr.
       !================================================

        if(iret.ne.0)then
          jt=996
          goto 1111 !error, redo psahot
        endif

      elseif(jt.le.13)then !one gamma involved

        if(ish.ge.6)write(ifch,*)'jt,jq,ncr:',jt,jq,ncr

       !================================================
        call putInfoPsreti(0,13,ibo,ncc(1,1),ncc(2,1),ncc(1,2),ncc(2,2))
        call psreti(ncr,jq,1
     .          ,ey,s0xi,c0xi,s0i,c0i,s0xh,c0xh,s0h,c0h,iret) !color conn. reconstruction
       !================================================

        ep3(1)=pprt(4,2)
        ep3(2)=pprt(3,2)
        ep3(3)=pprt(1,2)
        ep3(4)=pprt(2,2)
        call psrotat(ep3,s0xh,c0xh,s0h,c0h)  !special rotation for photon.
        call pstrans(ep3,ey,1)
        nptl=nptl+1
        if(jt.eq.11 .or.jt.eq.12)then
          desptl(nptl)= 0 ! ??? factgam ???   !individual weight for photons
        else                             
          desptl(nptl)=0 ! ???
        endif              
        pptl(1,nptl)=ep3(3)
        pptl(2,nptl)=ep3(4)
        pptl(3,nptl)=ep3(2)
        pptl(4,nptl)=ep3(1)
        pptl(5,nptl)=0
        idptl(nptl)=10
        iorptl(nptl)=iptl
        istptl(nptl)=0
        jorptl(nptl)=0
        do i=1,4
          xorptl(i,nptl)=coordo(i)
        enddo
        tivptl(1,nptl)=coordo(5)
        tivptl(2,nptl)=coordo(6)
        ityptl(nptl)=71
        ifrptl(1,nptl)=0
        ifrptl(2,nptl)=0
        qsqptl(nptl)=qt2
        zpaptl(1,nptl)=zpaptl(1,iptl)
        zpaptl(2,nptl)=zpaptl(2,iptl)
        pprtx(4,2)=pptl(4,nptl) !bg Out born photon
        pprtx(3,2)=pptl(3,nptl) !bg pprtx used in psphwr for out born partons
        pprtx(1,2)=pptl(1,nptl) !bg Here, photons are always in the "2" side
        pprtx(2,2)=pptl(2,nptl)
        idprtx(2)=10

      else !two gammas involed

        if(ish.ge.6)write(ifch,*)'jt,iqc:',jt,iqc
        do j=1,2
          ep3(1)=pprt(4,j)
          ep3(2)=pprt(3,j)
          ep3(3)=pprt(1,j)
          ep3(4)=pprt(2,j)
          call psrotat(ep3,s0xh,c0xh,s0h,c0h)  !special rotation for photon.
          call pstrans(ep3,ey,1)
          nptl=nptl+1
          desptl(nptl)= 0 !individual weight for photons
          pptl(1,nptl)=ep3(3)
          pptl(2,nptl)=ep3(4)
          pptl(3,nptl)=ep3(2)
          pptl(4,nptl)=ep3(1)
          pprtx(4,j)=pptl(4,nptl)
          pprtx(3,j)=pptl(3,nptl) !bg pprtx used in psphwr for out born partons
          pprtx(1,j)=pptl(1,nptl) !bg Here, two photons
          pprtx(2,j)=pptl(2,nptl)
          idprtx(j)=10
          pptl(5,nptl)=0
          idptl(nptl)=10
          iorptl(nptl)=iptl
          istptl(nptl)=0
          jorptl(nptl)=0
          do i=1,4
            xorptl(i,nptl)=coordo(i)
          enddo
          tivptl(1,nptl)=coordo(5)
          tivptl(2,nptl)=coordo(6)
          ityptl(nptl)=72
          ifrptl(1,nptl)=0
          ifrptl(2,nptl)=0
          qsqptl(nptl)=qt2
          zpaptl(1,nptl)=zpaptl(1,iptl)
          zpaptl(2,nptl)=zpaptl(2,iptl)
        enddo

      endif

      call utprix('psabor',ish,ishini,5)

c      print *,'psabor: jt      iqc ',jt,qt,iqc

      ienvi=110  !psabor exit
      return
      
 1111 continue

      if(ish.ge.6)write(ifch,*)'Escape psabor : ',jt
      !write(ifmt,*)'color connections(final)', ncr
      !do ij=1,nj
      !write(ifmt,*)'      ',iqj(ij),ij,(ncj(j,ij),j=1,2),'   ',jt
      !enddo

      ienvi=110  !psabor exit
      return
      end

c-----------------------------------------------------------------------
      subroutine switchComponents(iOB)
c-----------------------------------------------------------------------
      integer iOB(2)
      isave=iOB(1)
      iOB(1)=iOB(2)
      iOB(2)=isave
      end

c-----------------------------------------------------------------------
      subroutine idtrafo_KW_SO(iKW,iSO,ii)
c-----------------------------------------------------------------------
      if(ii.eq.1)then ! iKW -> iSO
        if(iabs(iKW).eq.3)then
          iSO=iKW*4/3
        elseif(iabs(iKW).ge.4)then
          iSO=iKW*10
        else
          iSO=iKW
        endif
      elseif(ii.eq.-1)then ! iKW <- iSO
        if(iabs(iSO).eq.4)then
          iKW=iSO/4*3
        elseif(iabs(iSO).eq.40.or.iabs(iSO).eq.50)then
          iKW=iSO/10
        else
          iKW=iSO
        endif
      else
         stop'ERROR idtrafo_KW_SO'
      endif
      end 
 


