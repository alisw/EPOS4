C
C  This file is part of EPOS4
C  Copyright (C) 2022 research institutions and authors (See CREDITS file)
C  This file is distributed under the terms of the GNU General Public License version 3 or later
C  (See COPYING file for the text of the licence)
C

c         ------------------------------------------------
c         since EPOS3086 we have om5 = om1
c                  (in earlier versions   om5 = 0.5 * om1)
c         change in : om52pp, om51p, om51, ffom12, ffom11
c         ------------------------------------------------

c###########################################################################
c###########################################################################
c###########################################################################
c###########################################################################
c
c          ffsig stuff  (extremely important for consistency checks !!!)
c
c            based on EPOS PDFs  (fpartone)
c
c###########################################################################
c###########################################################################
c###########################################################################
c###########################################################################

c-----------------------------------------------------------------------
      integer function noflav(qq) !KW1811 
c-----------------------------------------------------------------------
#include "sem.h"
#include "aaa.h"
      facfl4=2. 
      facfl5=2.5 
      nf=5                    
      if(qq.lt.(facfl4*qcmass)**2) then    
        nf=3                              
      elseif(qq.lt.(facfl5*qbmass)**2)then 
        nf=4                               
      endif           
      noflav=nf           
      end

c-----------------------------------------------------------------------
      integer function noflav0(qq) !TP1912
c-----------------------------------------------------------------------
#include "sem.h"
#include "aaa.h"
      nf=5                    
      if(qq.lt.qcmass**2) then    
        nf=3                              
      elseif(qq.lt.qbmass**2)then 
        nf=4                               
      endif           
      noflav0=nf           
      end

c-----------------------------------------------------------------------
      double precision function ffsig(fpdf,sis,klas,qt,x1,x2,ihq)     !KW1811 charm added
c-----------------------------------------------------------------------
c
c     fpdf(x1) * fpdf(x2) *  Born
c
c          Born = sum |M^2| / g^4 as tabulated in psbori  
c
c          fpdf should be f(x)*x, with f being the usual PDF 
c
c              acually fpartone(x) = f(x)*x
c
c-----------------------------------------------------------------------
c fpdf ... external function used as PDF
c sis .... ladder s or total (pp) s (Mandelstam s)
c klas ... Born kinematics class
c qt ..... squared trasverse momentum
c x1,x2 .. light cone momentum fractions
c-----------------------------------------------------------------------
c possible to select contributions from born with (p=g,q ; q=u,d,s ; Q=c,b) : 
c ihq = 0 - all
c     = 1 - pp' -> pp' 
c     = 2 - pp~ -> p'p'~
c     = 3 - qq~ -> gg
c     = 4 - charm prod
c     = 5 - bottom prod
c     <0 same but below Q2s
c-----------------------------------------------------------------------
#include "sem.h"
#include "aaa.h"
      double precision ffborn,fpdf,x1,x2,z,t,s,  smin,ttt,tmax
     .,f(-5:5),h(-5:5)
      external fpdf
      integer ioTestFact
      common /ctestfact/ ioTestFact
      integer laddTestFact
      common /claddtestfact/ laddTestFact
      integer ihqTestFact
      common /chqtestfact/ ihqTestFact
c      dimension q2cminsave(2)

      s=sis*x1*x2
      ffsig=0d0
      z=x1*x2

      if(ihqTestFact.eq.1)then !outBorn at least one HQ
        if(klas.eq.1.or.klas.eq.6.or.klas.eq.7)return
      elseif(ihqTestFact.eq.2)then !outBorn at least one b
        if(klas.eq.1.or.klas.eq.6.or.klas.eq.7
     . .or.klas.eq.2.or.klas.eq.4.or.klas.eq.8.or.klas.eq.10)return
      elseif(ihqTestFact.eq.3)then !outBorn at least one c no b
        if(klas.eq.1.or.klas.eq.6.or.klas.eq.7
     . .or.klas.eq.3.or.klas.eq.5.or.klas.eq.9.or.klas.eq.11
     . .or.klas.eq.12)return
      elseif(ihqTestFact.eq.4)then !outBorn one c and one b
        if(klas.ne.12)return
      endif

      call HardScale(klas,sc,qt,2)      !get scale sc
      call getBornKin2(klas,s,dble(qt)  ,  smin,ttt,tmax,iret) !get t and tmax
      t=ttt

      if(iret.ne.0)return
      
      n1=1
      n2=3
      factff=1.
      if(abs(ihq).le.3)then
        if(abs(ihq).gt.0)then
          n1=abs(ihq)
          n2=abs(ihq)
        endif
      else
        if(abs(ihq).eq.4)then
          if(klas.eq.2.or.klas.eq.12)then
            factff=0.5 !only one charm
          elseif(klas.eq.4.or.klas.eq.8.or.klas.eq.10)then
            factff=1.  !pair of charm
          else
            return
          endif
        elseif(abs(ihq).eq.5)then
          if(klas.eq.3.or.klas.eq.12)then
            factff=0.5 !only one bottom
          elseif(klas.eq.5.or.klas.eq.9.or.klas.eq.11)then
            factff=1.  !pair of bottom
          else
            return
          endif
        else
          return
        endif
      endif

c      if(qt.ge.max(q2cmin(1),q2cmin(2)))then

        is=1
        factbqq=1.
        factbgg=1.

        do i=-5,5
          f(i)=0
          h(i)=0
        enddo
        nf=noflav(sc)
        do i=-nf,nf
          f(i)=fpdf(1,x1,sc,i,2,0)
          h(i)=fpdf(2,x2,sc,i,2,0)
        enddo

        gg =   f(0) * h(0)

        gq =         f(0)     *    (h(1)+h(2)+h(3)+h(-1)+h(-2)+h(-3))
     .     +  (f(1)+f(2)+f(3)+f(-1)+f(-2)+f(-3))     *     h(0)     

        qq =    f( 1)*h( 1)  +  f( 2)*h( 2)  +  f( 3)*h( 3)
     .     +    f(-1)*h(-1)  +  f(-2)*h(-2)  +  f(-3)*h(-3)

        qa =  f( 1)*h(-1)  +  f( 2)*h(-2)  +  f( 3)*h(-3)
     .     +  f(-1)*h( 1)  +  f(-2)*h( 2)  +  f(-3)*h( 3)

        qqp=  f(-3)  *  (h(-2)+h(-1)+h(1)+h(2)) 
     .     +  f(-2)  *  (h(-3)+h(-1)+h(1)+h(3)) 
     .     +  f(-1)  *  (h(-3)+h(-2)+h(2)+h(3)) 
     .     +  f( 3)  *  (h(-2)+h(-1)+h(1)+h(2)) 
     .     +  f( 2)  *  (h(-3)+h(-1)+h(1)+h(3)) 
     .     +  f( 1)  *  (h(-3)+h(-2)+h(2)+h(3)) 

        gc =  f(0) * (h(4)+h(-4))  +  (f(4)+f(-4)) * h(0)

        gb =  f(0) * (h(5)+h(-5))  +  (f(5)+f(-5)) * h(0) 

        qc =  (f(1)+f(2)+f(3)+f(-1)+f(-2)+f(-3))   *   (h(4)+h(-4))  
     .     +    (f(4)+f(-4))    *    (h(1)+h(2)+h(3)+h(-1)+h(-2)+h(-3))

        qb =  (f(1)+f(2)+f(3)+f(-1)+f(-2)+f(-3))   *   (h(5)+h(-5))  
     .     +    (f(5)+f(-5))    *    (h(1)+h(2)+h(3)+h(-1)+h(-2)+h(-3))

        cc =  f(4)*h(4)   +  f(-4)*h(-4)  

        cac=  f(4)*h(-4)  +  f(-4)*h(4)  

        bb =  f(5)*h(5)   +  f(-5)*h(-5)  

        bab=  f(5)*h(-5)  +  f(-5)*h(5)  

        cb =  (f(4)+f(-4)) *  (h(5)+h(-5)) 
     .     +  (f(5)+f(-5)) *  (h(4)+h(-4)) 

c      else
c        return
c      endif
      
      ffsig= ffsig + 
     *         ffborn(0,klas,tmax,s,t,gg,gq,qq,qa,qqp ,n1,n2, is*3 ) +
     *         ffborn(0,klas,tmax,s,t,0.,gc,0.,0. ,qc ,n1,n2, is*4 ) +
     *         ffborn(0,klas,tmax,s,t,0.,gb,0.,0. ,qb ,n1,n2, is*5 ) +
     *         ffborn(0,klas,tmax,s,t,0.,0.,cc,cac,0. ,n1,n2, is*6 ) +
     *         ffborn(0,klas,tmax,s,t,0.,0.,bb,bab,cb ,n1,n2, is*7 )
      if(ish.ge.9)then
        write(ifmt,'(a,6x,3e12.3,3x,$)')'ffsig  = ',ffsig/z,s,t !Sudakov not included
        write(ifmt,'(5e12.3)')gg/z,gq/z,qq/z,qa/z,qqp/z  
        !write(ifmtx,'(3x,6e12.3,$)')gc/z,qc/z,gb/z,qb/z,cc/z,cac/z
        !write(ifmtx,'(3e12.3)')bb/z,bab/z,cb/z
      endif

      ffsig=ffsig*dble(factff)

c     print *,'sig',s,t,sc,x1,x2,g1,uv1,dv1,sea1,g2,uv2,dv2,sea2,ffsig
c      q2cmin(1)=q2cminsave(1)
c      q2cmin(2)=q2cminsave(2)
      end

c------------------------------------------------------------------------
      double precision function ffsigSumKlas(fpdf,sis,qq,x1,x2,ihq) 
c------------------------------------------------------------------------
c      \sum_klas ffsig * pi * alphas^2
c------------------------------------------------------------------------
#include "sem.h"
#include "aaa.h"
      double precision fpdf,x1,x2,ffsig,ft
      external fpdf
      ft=0
      do klas=1,klasmax
        call HardScale(klas,sc,qq,2)      !get scale sc
        ft=ft+ffsig(fpdf,sis,klas,qq,x1,x2,ihq) 
     .        * pi * (2*pi*pssalf(sc/qcdlam))**2  !*pi*alphas^2
      enddo
      ffsigSumKlas=ft
      end


c=======================================================================
      real function ffsigi(qq,y0,ihqx,jj)    !*** called from ems ***
c=======================================================================
      implicit none
      real qq,y0,ffsigii
      double precision pifpartone,fpartone,xCTEQ2
      external pifpartone,fpartone,xCTEQ2
      integer ioTestFact,ihqx,jj
      common /ctestfact/ ioTestFact
      if(jj.eq.6)then 
        ffsigi=ffsigii(xCTEQ2,qq,y0,ihqx)!CTEQ14
      elseif(jj.eq.1)then 
        if(ioTestFact.ne.0.and.ioTestFact.ne.12)then
          ffsigi=ffsigii(fpartone,qq,y0,ihqx) !ffsigiixxx(qq,y0)!(test)
        else
          ffsigi=ffsigii(pifpartone,qq,y0,ihqx)
        endif
      else
        stop'ERROR 14042022'
      endif
      end
c=======================================================================
 

c------------------------------------------------------------------------
      function ffsigii(fpdf,qq,y0,ihqx)             !former psjx1  (sto)
c------------------------------------------------------------------------
c
c  \int dy dx ffsigSumKlas / s**2 
c                          / 2   ! y interval  2 * Delta_y
c                          / 2   ! condition t < sqrt(s)/2, t > sqrt(s)/2 included in psbori
c
c      ffsigSumKlas = {\sum_klas ffsig} * pi * alphas**2  
c
c        ffsig =  fpdf(x1) * fpdf(x2) *  Born (with fpdf(x)=x*f(x))
c
c        Born = sum |M^2| / g^4 as tabulated in psbori  
c
c          attention: fpartone(x) = f(x)*x, with f being the usual PDF 
c
c------------------------------------------------------------------------
c   computation of jet cross section dsigma_jet/dpt 
c      in rapidity range -y0,y0] divided by 2
c      qq = pt**2,  
c------------------------------------------------------------------------ 
c ihq = 0 - all
c     = 1 - pp' -> pp' 
c     = 2 - pp~ -> p'p'~
c     = 3 - qq~ -> gg
c     = 4 - gg -> QQ~
c     = 5 - qq~ -> QQ~
c     = 11- sat (but integration might not be correct (q2sft<qq<qqcut) ...)
c if ihq<0, cross-section without emission
c-----------------------------------------------------------------------
#include "sem.h"
#include "aaa.h"
      double precision ffsigSumKlas,ft,ffsigid,fx,s,z,xx,fpdf,ylim
     &                 ,sh,xx1,xx2,xt,ymax,ymin,y,xmin,xmax,xmax0
      external fpdf
      integer ioTestFact
      common /ctestfact/ ioTestFact
      integer laddTestFact
      common /claddtestfact/ laddTestFact
      integer modeDeltaFzero
      common /cmodeDeltaFzero/ modeDeltaFzero

      ig=7
      ig1=7
      s=dble(engy)**2
      !kkkkkkkkkkkkkkkkkkkkk
      if(ioTestFact.ne.0.and.ioTestFact.le.5)s=dble(engy)**2    
      if(ioTestFact.ge.6)s=dble(engy)**2 *0.9**2   
      !kkkkkkkkkkkkkkkkkkkkk
      sis=s
      ihq=ihqx

      ffsigii=0.

      if(s.le.4d0*dble(qq))return
      if(qq.lt.q2sft)return
      xmax0=1.d0

      qqcut=max(q2cmin(1),q2cmin(2))     
      xt=2d0*sqrt(dble(qq)/s)
      ymax=min(dble(y0),log(1d0/xt+sqrt((1d0/xt-1d0)*(1d0/xt+1d0))))
      ymin=-ymax

      if(abs(y0).ge.1000)then
        yi=mod(abs(y0),1000.)
        yi=sign(yi,y0)
        ymin=yi-0.5
        ymax=yi+0.5
        ylim=log(1d0/xt+sqrt((1d0/xt-1d0)*(1d0/xt+1d0)))
        if(ymin.gt.ylim)return 
        if(ymax.lt.-ylim)return 
        ymax=min(ymax,ylim)
        ymin=max(ymin,-ylim)
      endif                     

      !--------------------------------------------------
      !    y integration
      !--------------------------------------------------

      !kkkkkkkkkkkkkkkkkkkkkkkkkkkkkkkkk!KW1908
      !kk  for modeDeltaFzero = 0  (normal cas)
      !kk  the following  contains "emission + no emission"  
      !kk         on both sides, 
      !kk  for modeDeltaFzero = 1 in includes
      !kk  only "emission" on both sides 
      !kk  (emission means at east one emission)
      !kkkkkkkkkkkkkkkkkkkkkkkkkkkkkkkkk

      ffsigid=0d0
      do i=1,ig
      do m=1,2
        y=.5d0*(ymax+ymin+(ymin-ymax)*dble((2*m-3)*tgss(ig,i)))
        xmin=xt**2/2.d0/(2.d0-xt*exp(-y))   !condition x2<1
        xmax=min(xmax0,1.d0-xt*exp(y)/2.d0) !condition x1<1
        fx=0.d0
        if(xmax.gt.xmin)then
          do i1=1,ig1
          do m1=1,2
            xx=xmin*(xmax/xmin)**dble(.5+tgss(ig1,i1)*(m1-1.5)) !Jacob=0.5*xx*log(xmax/xmin)
            xx1= xt*exp(y)/2 + xx
            xx2= xt*exp(-y)*xx1 / 2 / xx !from s+t+u=0
            z=xx1*xx2
            sh=z*s
            ft=ffsigSumKlas(fpdf,sis,qq,xx1,xx2,ihq)  
            fx=fx+dble(wgss(ig1,i1))*ft/sh**2
            !~~~~~~~~~~~~~TEST~~~~~~~~~~~~~~~~~~~~~~~~
            !if(igetIsh().eq.11)then 
            ! if(i.eq.1.and.m.eq.2.and.m1.eq.1)then
            ! call HardScale(1,sc,qq,2)      !get scale sc
            ! write(ifmt,'(a,2f9.3,2e10.3,$)')'FSIGII TEST',y,sc,xx1,xx2
            ! write(ifmt,'(f16.0,f12.6,$)')ft,fpdf(1,xx1,sc,0,2,0)
            ! write(ifmt,'(f12.6)')fpdf(2,xx2,sc,0,2,0)
            ! endif
            !endif
            !~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
          enddo
          enddo
          fx=fx*0.5d0*log(xmax/xmin)    ! xx-Jacob / xx 
        endif
        ffsigid=ffsigid+dble(wgss(ig,i))*fx
      enddo
      enddo
      ffsigii = ffsigii +  ffsigid * (ymax-ymin)/2  ! y-Jacob

      if(ioTestFact.ne.0.and.modeDeltaFzero.eq.0)then
        write(ifmt,'(a,f9.2,e12.3)')'ffsigii',sqrt(qq),ffsigii
      endif

      !kkkkkkkkkkkkkkkkkkkkkkkkkkkkkkkkkkkkk!KW1908
      !kk  for modeDeltaFzero = 1  ( fsoft = delta(1-x) )
      !kk  add explicitely no-emission contributions, 
      !kk not contibuting  to the above integral
      !kk due to delta(1-x) function in fpartone
      !kkkkkkkkkkkkkkkkkkkkkkkkkkkkkkkkkkkk

      if(modeDeltaFzero.eq.1)then  

      !only lower side emissions 

      ffsigid=0d0
      do i=1,ig
      do m=1,2
        y=.5d0*(ymax+ymin+(ymin-ymax)*dble((2*m-3)*tgss(ig,i)))
        fx=0.d0
        xx1=1d0
        xx=xx1-xt*exp(y)/2 
        xx2= xt/2*exp(-y)*xx1 / xx !from s+t+u=0
        if(xx2.lt.xmax0)then
          z=xx1*xx2
          sh=z*s
          ft=ffsigSumKlas(fpdf,sis,qq,xx1,xx2,ihq)  
          fx=fx+ft/sh**2
        endif
        ffsigid=ffsigid+dble(wgss(ig,i))*fx
      enddo
      enddo
      ffsigii = ffsigii +  ffsigid * (ymax-ymin)/2  ! y-Jacob

      !only upper side emissions 

      ffsigid=0d0
      do i=1,ig
      do m=1,2
        y=.5d0*(ymax+ymin+(ymin-ymax)*dble((2*m-3)*tgss(ig,i)))
        fx=0.d0
        xx2=1d0
        xx=xx2-xt*exp(-y)/2 
        xx1= xt/2*exp(y)*xx2 / xx !from s+t+u=0
        if(xx1.lt.xmax0)then
          z=xx1*xx2
          sh=z*s
          ft=ffsigSumKlas(fpdf,sis,qq,xx1,xx2,ihq) !qq is pt^2, transformed to t 
          fx=fx+ft/sh**2
        endif
        ffsigid=ffsigid+dble(wgss(ig,i))*fx
      enddo
      enddo
      ffsigii = ffsigii +  ffsigid * (ymax-ymin)/2  ! y-Jacob

      !no emissions

      ffsigid=0d0
      fx=0.d0
      xx2=1d0
      xx1=1d0
      y=acosh(1./xt) !from s+t+u=0 (see pub/19epos4) !KW2305 ???????????? actually two solutions +-acosh(..) -> factor 2
      if(y.ge.ymin-0.001.and.y.le.ymax+0.001)then
        z=xx1*xx2
        sh=z*s
        ft=ffsigSumKlas(fpdf,sis,qq,xx1,xx2,ihq) 
        fx=fx+ft/sh**2
      endif
      ffsigid=ffsigid+fx
      ffsigii = ffsigii + ffsigid 
 
      endif !modeDeltaFzero.eq.1

      !-----------------------------------------------
      !   common factors
      !-----------------------------------------------

      !factor pi * alphas^2 already included in ffsigSumKlas
      !factor 2 for jets (rather than dijets) included in ems
 
      ffsigii=ffsigii   
     .   /2   ! y interval  2 * Delta_y
     .   /2   ! two contributions, t < tmax+, t > tmax, included in psbori
              !  which is useful for integrated xsections, to integrate t<tmax+
              !  but for inclusive xsection we integrate all

      ffsigii=ffsigii  *2*sqrt(qq)    ! dsig/dpt = 2*pt*dsig/dpt^2

      !nexus: factor  pt (in addition to pi*alpha^2)
      !        ours is 2 times smaller, counting dijets ?  

      return
      end

c------------------------------------------------------------------------
      double precision function fpartone(iii,xx,qq,j,je,ji)
c   former pspdf0 (sha)
c-----------------------------------------------------------------------
c
c  parton distribution function for proton  ( actually x*f(x) !!!!!!! )
c
c iii = proj or targ
c xx = light cone momentum fraction
c qq = virtuality scale
c j = flavor 
c          <0 ... antiquarks
c          0 .... gluon
c          >0 ... quarks  
c je = emission type
c          0 ... no emissions
c          1 ... emissions
c          2 ... all
c ji = quark type
c          0 sum of sea + val
c          1 initial sea + evolution from sea and val
c          2 initial val + evolution from val only      
c - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
c   integral over fzeroGlu, fzeroVal,FzeroQua
c-----------------------------------------------------------------------
      common/ar3/    x1(7),a1(7)
#include "aaa.h"
#include "sem.h"
      double precision fz,zx,z,akns,dpd1,dpd2,epscutSud2
     *,psuds,xx,xmin,xm,dnoflav 
      double precision fzeroVal,fzeroGlu,FzeroQua
      integer ioTestFact
      common /ctestfact/ ioTestFact
      integer laddTestFact
      common /claddtestfact/ laddTestFact
      integer iSwitchSides
      common /cSwitchSides/iSwitchSides
      double precision sef(-5:5),aksf(2),qua,sef0(-5:5)
      double precision xdelta
      common/cdelta/xdelta
      integer modeDeltaFzero
      common /cmodeDeltaFzero/ modeDeltaFzero

      fpartone=0d0

      if(qq.lt.2.*qcdlam)goto 999

      ii=iii
      if(iSwitchSides.eq.1)ii=3-iii

      xdelta=0.99999999999999d0 !to define delta function
      !So for q2cmin>q2sft, nf will increase and new contributions will be used from Esatur
      nf=noflav(min(qq,q2cmin(ii)))

      nfqq=noflav(qq)
      dnoflav=dble(nfqq)

      nfqq0=noflav0(qq)
      !some charm is produced between qcmass**2 (noflav0) and 4*qcmass**2 (noflav) from pQCD evolution from non charm (results compatible with data).
      ! but to avoid HQ below mass threshold this is necessary (to avoid some contribution from pQCD evolution)
      if(abs(j).gt.nfqq0)return
c      if(dnoflav.lt.5d0.and.j.eq.-5)return
c      if(j.eq.0.and.ji.eq.1)then
c        fpartone=psdpdf(sngl(xx),qq,0.,iclpro,j)
c      elseif((j.eq.1.or.j.eq.2).and.ji.eq.2)then
c        fpartone=psdpdf(sngl(xx),qq,0.,iclpro,j)
c      elseif(j.le.-1.and.ji.eq.1)then
c        fpartone=psdpdf(sngl(xx),qq,0.,iclpro,j)
c      endif
c      return

c      print *,'fpartone',iii,xx,qq,j,je,ji, ioTestFact,laddTestFact
      
      if(laddTestFact.eq.2.and.ii.eq.1)goto 888 !emiss proj
      if(laddTestFact.eq.3.and.ii.eq.2)goto 888 !emiss targ
      if(laddTestFact.eq.4)goto 888  !emiss both sides

      if(je.eq.1)goto 888

c ...... f_0 * sudakov.........

      !kkkkkkkkkkkkkkkkkkkkkkkkkkkkkkkkkkkkkkkkkkkkkkk!KW1908
      !kk  modeDeltaFzero = 1
      !kk  fzero functions are considered to be delta(xx-1)
      !kk     So here we set fpartone=Sudakov for xx=1
      !kk                and zero elsewhere
      !kk     BUT in addition, special care is needed when 
      !kk     integrating over fpartone(..xx..), see ffsigi
      !kkkkkkkkkkkkkkkkkkkkkkkkkkkkkkkkkkkkkkkkkkkkkkk

      if(modeDeltaFzero.eq.1.and.xx.le.xdelta)then
        fpartone=0
      else !modeDeltaFzero=0 or xx>xdelta
        i=j
        sef(i)=0d0
        if(abs(i).le.nf)then
          if(i.eq.0)then
            if(ji.ne.2)sef(i)=fzeroGlu(xx,min(qq,q2cmin(ii)),ii)
          elseif(i.eq.1.or.i.eq.2)then
           sef(i)=0d0
           if(ji.eq.0.or.ji.eq.2)sef(i)=sef(i)
     .                          +fzeroVal(xx,min(qq,q2cmin(ii)),ii,i) ! do not add valence quark if we want sea only
           if(ji.le.1)sef(i)=sef(i)+FzeroQua(xx,min(qq,q2cmin(ii)),ii,i) !per flavor (and not if valence quark contribution only)
          else
            sef(i)=FzeroQua(xx,min(qq,q2cmin(ii)),ii,i) !per flavor
          endif
        endif
        fpartone=sef(i)
c       no pertubative calculation below q2sft and non linear pertubative calculation nelow q2cmin, in both cases there no evolution and then no need for sudakov factor
        if(qq.lt.q2cmin(ii))goto 999
        fpartone=fpartone*psuds(qq,j)/psuds(q2cmin(ii),j)
      endif

      if(laddTestFact.eq.1)goto 999  !no emiss
      if(laddTestFact.eq.2.and.ii.eq.2)goto 999 !no emiss targ
      if(laddTestFact.eq.3.and.ii.eq.1)goto 999 !no emiss proj

      if(je.eq.0)goto 999

c......... integral f_0 E_qcd............

 888  continue
      
      if(qq.lt.q2cmin(ii))goto 999

      xmin=xx/(1.d0-epscutSud2(dble(qq)))

      !kkkkkkkkkkkkkkkkkkkkkkkkkkkkkkkkkkkkkkkkkkkkkkk!KW1908
      !kk  modeDeltaFzero=1:  
      !kk  fzero functions are considered to be delta(1-zx)
      !kk      -->  no integration
      !kkkkkkkkkkkkkkkkkkkkkkkkkkkkkkkkkkkkkkkkkkkkkkk

      if(modeDeltaFzero.eq.1)then
       fpartone=0d0
       if(xx.le.xdelta)then
        zx=xdelta
        z=xx/zx
        do n=-5,5
          sef(n)=0
          sef0(n)=0
        enddo
        if(ji.ne.2)then
          sef(0)=fzeroGlu(zx,q2cmin(ii),ii)
          sef(1)=fzeroVal(zx,q2cmin(ii),ii,1)
          sef(2)=fzeroVal(zx,q2cmin(ii),ii,2)
          if(ji.eq.1)then
            sef0(1)=sef(1)
            sef0(2)=sef(2)
          endif
          do n=-nf,nf
            if(n.ne.0)then      !no sea contribution for gluons or if we want valence only
              qua=FzeroQua(zx,q2cmin(ii),ii,n) !allow n dependence
              sef(n)=sef(n)+qua
            endif
          enddo
        else !in case we want valence quark contribution only, only valence from defined flavor j is used
          sef(j)=fzeroVal(zx,q2cmin(ii),ii,j)

        endif
        fz=0.d0
        if(j.eq.0)then
          if(ji.ne.2)then
            do n=-nf,nf
              k=min(abs(n),1)+1
              fz=fz+sef(n)*psevi(q2cmin(ii),qq,z,k,1)
            enddo
          endif
          !~~~~~~~~~~~~~~
          !if(xx.gt.0.99d0)then
          !ish=9 
          !print*,z,psevi(q2cmin(ii),qq,z,1,1)
          !stop'FPARTONE TEST'
          !endif
          !~~~~~~~~~~~~~~~
        else
          akns=psevi(q2cmin(ii),qq,z,3,2) !nonsinglet contribution
          fz=max(0d0,sef(j)-sef0(j))*akns 
          if(ji.ne.2)then  !no singlet contribution for val only
            aksf(2)=(psevi(q2cmin(ii),qq,z,2,2)-akns) /dnoflav/2.d0 !singlet contribution
            aksf(1)=psevi(q2cmin(ii),qq,z,1,2) /dnoflav/2.d0
            do n=-nf,nf
              k=min(abs(n),1)+1
              fz=fz+sef(n)*aksf(k)
            enddo
          endif
        endif
        fpartone  = fz   *  xx    ! <-------------------------
       endif
       goto 999
      endif 

 !integrations \int d(zx) { 1/zx**2 * zx*f(zx) * E(xx/zx) } * xx   with  zx*f(zx) = fzero(zx)
 
      if(xmin.lt.1.d0)then

        dpd1=0.d0
        dpd2=0.d0
        xm=max(xmin,0.3d0)
 
 !numerical integration xm -> 1

        do i=1,7
        do m=1,2
          zx=1.d0-(1.d0-xm)*(.5d0+(dble(m)-1.5d0)*dble(x1(i)))**.25d0
          z=xx/zx
          do n=-5,5
            sef(n)=0
            sef0(n)=0
          enddo
          if(ji.ne.2)then
            sef(0)=fzeroGlu(zx,q2cmin(ii),ii)
            sef(1)=fzeroVal(zx,q2cmin(ii),ii,1)
            sef(2)=fzeroVal(zx,q2cmin(ii),ii,2)
            if(ji.eq.1)then
              sef0(1)=sef(1)
              sef0(2)=sef(2)
            endif
            do n=-nf,nf
              if(n.ne.0)then !no sea contribution for gluons or if we want valence only
                qua=FzeroQua(zx,q2cmin(ii),ii,n) !allow n dependence
                sef(n)=sef(n)+qua
              endif
            enddo
          else                  !in case we want valence quark contribution only, only valence from defined flavor j is used
            sef(j)=fzeroVal(zx,q2cmin(ii),ii,j)
          endif
          fz=0.

          if(j.eq.0)then
            if(ji.ne.2)then
              do n=-nf,nf 
                k=min(abs(n),1)+1
                fz=fz+sef(n)*psevi(q2cmin(ii),qq,z,k,1)
              enddo
            endif
          else
            akns=psevi(q2cmin(ii),qq,z,3,2) !nonsinglet contribution
            fz=max(0d0,sef(j)-sef0(j))*akns 
            if(ji.ne.2)then  !no singlet contribution for val only
              aksf(2)=(psevi(q2cmin(ii),qq,z,2,2)-akns) /dnoflav/2.d0 !singlet contribution
              aksf(1)=psevi(q2cmin(ii),qq,z,1,2) /dnoflav/2.d0
              do n=-nf,nf
                k=min(abs(n),1)+1
                fz=fz+sef(n)*aksf(k)
              enddo
            endif
          endif
          dpd1=dpd1 + dble(a1(i)) * fz
     *                          / zx**2  
     *                   /(1.d0-zx)**3 !Jacob
        enddo
        enddo
        dpd1=dpd1 * (1.d0-xm)**4 / 4.d0 * 0.5d0 !Jacob
     *                      * xx   ! <-------------------------

 !numerical integration  xmin -> xm

        if(xm.gt.xmin)then

          do i=1,7
          do m=1,2
            zx=xx+(xm-xx)
     &         *((xmin-xx)/(xm-xx))**(.5d0-(dble(m)-1.5d0)*dble(x1(i)))
            z=xx/zx
            do n=-5,5
              sef(n)=0
              sef0(n)=0
            enddo
            if(ji.ne.2)then
              sef(0)=fzeroGlu(zx,q2cmin(ii),ii)
              sef(1)=fzeroVal(zx,q2cmin(ii),ii,1)
              sef(2)=fzeroVal(zx,q2cmin(ii),ii,2)
              if(ji.eq.1)then
                sef0(1)=sef(1)
                sef0(2)=sef(2)
              endif
              do n=-nf,nf
                if(n.ne.0)then !no sea contribution for gluons or if we want valence only
                  qua=FzeroQua(zx,q2cmin(ii),ii,n) !allow n dependence
                  sef(n)=sef(n)+qua
                endif
              enddo
            else                !in case we want valence quark contribution only, only valence from defined flavor j is used
              sef(j)=fzeroVal(zx,q2cmin(ii),ii,j)
            endif
            fz=0.d0
            if(j.eq.0)then
              if(ji.ne.2)then
                do n=-nf,nf 
                  k=min(abs(n),1)+1
                  fz=fz+sef(n)*psevi(q2cmin(ii),qq,z,k,1)
                enddo
              endif
            else
              akns=psevi(q2cmin(ii),qq,z,3,2) !nonsinglet contribution
              fz=max(0d0,sef(j)-sef0(j))*akns 
              if(ji.ne.2)then  !no singlet contribution for val only
                aksf(2)=(psevi(q2cmin(ii),qq,z,2,2)-akns) /dnoflav/2.d0 !singlet contribution
                aksf(1)=psevi(q2cmin(ii),qq,z,1,2) /dnoflav/2.d0
                do n=-nf,nf
                  k=min(abs(n),1)+1
                  fz=fz+sef(n)*aksf(k)
                enddo
              endif
            endif
            dpd2=dpd2+dble(a1(i))*fz
     *                         / zx**2   ! <-------------------------
     *                      * (zx-xx) !Jacob
          enddo
          enddo
          dpd2=dpd2 * log((xm-xx)/(xmin-xx)) * .5d0 !Jabob
     *                       *  xx    ! <-------------------------
        endif
        fpartone=fpartone+dpd2+dpd1

      endif


  999 continue

      !if(j.lt.0)then

      !kkkkkkkkkkkkkkkkk !KW1908
      !kk  I completely removed this part (see version 3275 if still needed)
      !kk     all these special trials should be  clearly separated from the main part
      !kk     and should be off by default
      !kkkkkkkkkkkkkkkk
         
      !endif

      fpa=fpartone     ! NaN catch
      if(.not.(fpa.le.0..or.fpa.ge.0.))then
        print *,fpartone,iii,xx,qq,j,je,ji,dpd2,dpd1
        stop'ERROR 25102012b'
      endif

      return
      end

c-----------------------------------------------------------------------
      double precision function pifpartone(ii,xx,qq,jj,je,ji)
c polynomial interpolation of partone
c qq is scale
c-----------------------------------------------------------------------
#include "sem.h"
#include "aaa.h"
      double precision fptn,x2tmin,x2tmax,epscutSud2
      parameter(npifmax=40,kxxmax=20,kqqmax=10)
      common/tabfptn/fptn(kxxmax,kqqmax,2,24,npifmax)
     *              ,tabfptn0(-5:2,0:2,0:0,npifmax)
      common/cq2tmin2/x2tmin(npifmax),x2tmax,q2tmin(npifmax)
     *               ,q2tmax(npifmax),q2tmx
      data npifpartone /0/
      save npifpartone
      double precision wi(3),wj(3),xxk,qqk,xx,xxmin,xxmax
     *,fpartone
c      common/q2tabmin/q2tab1min(10),q2tab2min(10)
      integer t!,ixx
c      common/testinter/testinterfpartone(60,20,42)
c     *,testq(60,20,42)
      integer ioTestFact
      common /ctestfact/ ioTestFact
      integer modeDeltaFzero
      common /cmodeDeltaFzero/ modeDeltaFzero
      double precision xdelta
      common/cdelta/xdelta

      !kkkkkkkkkkkkkkkkkkkkkkkkkkkkkkkkkkkkkkkkkkkkkkk
      !kk  modeDeltaFzero=1:  
      !kk  fzero functions are considered to be delta(1-xx)
      !kk     Here we set pifpartone=fpartone for xx=1
      !kk                and zero elsewhere
      !kk     BUT in addition, special care is needed when 
      !kk     integrating over fpartone(..xx..), see ffsigi
      !kkkkkkkkkkkkkkkkkkkkkkkkkkkkkkkkkkkkkkkkkkkkkkk
      if(modeDeltaFzero.eq.1)then
        if(xx.gt.xdelta)then 
          pifpartone=fpartone(ii,xx,qq,jj,je,ji)
          return
        else
          pifpartone=0
          return
        endif
      endif

      if(qq.le.q2cmin(ii).or.ji.ne.0.or.qq.gt.1e12)then
        pifpartone=fpartone(ii,xx,qq,jj,je,ji)
        return
      endif

      nnn=0
      do  n=1,min(npifmax,npifpartone)
c ???? 0.25 can be decreased for better precision
        if(abs(q2cmin(ii)-q2tmin(n)).lt.0.25
     *               .and.qq.le.q2tmax(n))then
          nnn=n
          goto 1
        endif
      enddo

 1    if(nnn.eq.0)then
        npifpartone=npifpartone+1
        if(npifpartone.lt.npifmax)then
          call MakeFpartonTable(nt,ii,qq)
          nnn=nt
        else
          nnn=npifmax
          x2tmin(nnn)=1d0
        endif
      endif

      if(xx.lt.x2tmin(nnn))then       !if xx is too low or limit of table reached, fpartone is used
        pifpartone=fpartone(ii,xx,qq,jj,je,ji)
        return
      endif

      j=jj
      if(j.ge.3)j=-j
      t=nint(tabfptn0(j,je,ji,nnn))

      if(t.eq.0)then  !return 0 for case where t=0

        pifpartone=0.d0
        return

      elseif(xx.eq.1.d0)then

        pifpartone=0.d0
        return

c      elseif(((t.eq.2).or.(t.eq.3).or.(t.eq.5).or.(t.eq.7)
c     *.or.(t.eq.8)
c     *.or.(t.eq.10).or.(t.eq.12).or.(t.eq.15)))then
c
c        if(qq.eq.q2cmin(ii))then
c          pifpartone=0.d0
c          return
c        endif
c        if(xx.ge.0.9d0)then
c          ixx=int(xx*100d0-90d0)+1
c          if(ixx.gt.10)ixx=10
c          if((t.eq.3).or.(t.eq.5))then
c            if(qq.le.q2tab1min(ixx))then
c              pifpartone=0.d0
c              return
c            endif
c          else
c            if(qq.le.q2tab2min(ixx))then
c              pifpartone=0.d0
c              return
c            endif
c
c          endif
c
c        endif

      elseif((je.eq.1).and.(xx/(1d0-epscutSud2(dble(qq))).gt.1))then

        pifpartone=0.d0
        return

      endif


      !if f>0, interpolate f

      t=nint(tabfptn0(j,je,ji,nnn))
      qqmin=q2tmin(nnn)
      qqmax=q2tmax(nnn)
      qqk=1.d0+dble(log(qq/qqmin)/log(qqmax/qqmin)*(kqqmax-1))
      kqq=int(qqk)
      if(kqq.lt.1)kqq=1
      if(kqq.gt.(kqqmax-2))kqq=kqqmax-2
      wj(2)=qqk-dble(kqq)
      wj(3)=wj(2)*(wj(2)-1.)*.5d0
      wj(1)=1.d0-wj(2)+wj(3)
      wj(2)=wj(2)-2.d0*wj(3)

      pifpartone=0.d0
      xxmin=x2tmin(nnn)

      xxmax=x2tmax
      if(xx.lt.xxmax.and.xxmin.lt.xxmax)then

        xxk=1.d0+log(xx/xxmin)/log(xxmax/xxmin)*(kxxmax-1)
        kxx=int(xxk)
        if(kxx.lt.1)kxx=1
        if(kxx.gt.(kxxmax-2))kxx=kxxmax-2

        wi(2)=xxk-dble(kxx)
        wi(3)=wi(2)*(wi(2)-1.)*.5d0
        wi(1)=1.d0-wi(2)+wi(3)
        wi(2)=wi(2)-2.d0*wi(3)

        do kq=1,3
        do kx=1,3
          pifpartone=pifpartone+fptn(kxx+kx-1,kqq+kq-1,1,t,nnn)
     *              *wi(kx)*wj(kq)
        enddo
      enddo

c        pifpartone=exp(pifpartone)

      else

        xxmax=1.d0


        xxk=1.d0+(xx-xxmin)/(xxmax-xxmin)*dble(kxxmax-1)
        kxx=int(xxk)
        if(kxx.lt.1)kxx=1
        if(kxx.gt.(kxxmax-2))kxx=kxxmax-2

        wi(2)=xxk-dble(kxx)
        wi(3)=wi(2)*(wi(2)-1.)*.5d0
        wi(1)=1.d0-wi(2)+wi(3)
        wi(2)=wi(2)-2.d0*wi(3)

        do kq=1,3
        do kx=1,3
          pifpartone=pifpartone+fptn(kxx+kx-1,kqq+kq-1,2,t,nnn)
     *              *wi(kx)*wj(kq)
        enddo
        enddo


      endif

c      print*,'pif',xx,kxx,qq,kqq,j,jj,je,ji,pifpartone

      pif=pifpartone     ! NaN catch
      if(.not.(pif.le.0..or.pif.ge.0.))then
        print*,'pif',xx,qq,j,je,ji,pifpartone
        stop'ERROR 25102012'
      endif

      pifpartone=max(0d0,pifpartone)

      end

c-----------------------------------------------------------------------
      subroutine MakeFpartonTable(nt,ii,qqi)
c-----------------------------------------------------------------------

#include "sem.h"
#include "aaa.h"
      double precision fpartone
      double precision fptn,x2tmin,x2tmax,xvalueEvTab
      parameter(npifmax=40,kxxmax=20,kqqmax=10)
      common/tabfptn/fptn(kxxmax,kqqmax,2,24,npifmax)
     *              ,tabfptn0(-5:2,0:2,0:0,npifmax)
      common/cq2tmin2/x2tmin(npifmax),x2tmax,q2tmin(npifmax)
     *               ,q2tmax(npifmax),q2tmx
      real t
      double precision xx,xxmin,xxmax
      data ncstable2/0/
      save ncstable2
      ncstable2=ncstable2+1
      ncstable2=min(npifmax,ncstable2)
      nt=ncstable2

      if(nt.eq.1)then
        x2tmax=0.1d0
        if(abs(noebin).gt.1)then !special case to plot F2 for instance
          if(engmax.le.0.)
     *    stop '/n/n Please define engmax to run pifpartone .../n/n'
          q2tmx=engmax**2/4.
        else
          q2tmx=engy**2/4.
        endif
      endif
      
      q2tmin(nt)=q2cmin(ii)
      q2tmax(nt)=min(1.e4*qqi,1e12)
      x2tmin(nt)=xvalueEvTab(1,dble(q2tmax(nt)))     !dble(q2tmin(nt)/q2tmx)

      write (ifmt,'(a,i3,1x,2(e8.2,1x),$)')'(FpartonTable',nt
     *                             ,q2tmin(nt),q2tmax(nt)

        t=0.
        do j=-5,2
         do je=0,2
          do ji=0,0
c          if(   ( (j.le.-1).and.(je.eq.0).and.(ji.eq.2) )  ! sea-qu val  no-emiss
c     *     .or. ( (j.eq.0).and.(je.eq.0).and.(ji.eq.2)  )  ! g      val  no-emiss
c     *     .or. ( ((j.eq.1).or.(j.eq.2)).and.(ji.eq.1)  )  ! u,d    sea 
c     *     )then  
c           tabfptn0(j,je,ji,ncstable2)=0.
c           else
           t=t+1.
           tabfptn0(j,je,ji,ncstable2)=t
c          endif
          enddo
         enddo
        enddo

      !first part: tab for extrapolate x<0.1, fill between 0 and 1
      !xx,qq logarithmic scale

      call ipoBasics(q2cmin)

      qqmin=q2tmin(nt)
      qqmax=q2tmax(nt)
      xxmin=x2tmin(nt)
      xxmax=x2tmax
      if(xxmin.lt.xxmax)then
      do j=-5,2
       do je=0,2
        do ji=0,0

          t=tabfptn0(j,je,ji,ncstable2)
          if(t.gt.0)then
            write(ifmt,'(a,$)')'.'
            do kqq=1,kqqmax
              qq=qqmin*(qqmax/qqmin)**((kqq-1.)/(kqqmax-1.))
              do kxx=1,kxxmax
                xx=xxmin*(xxmax/xxmin)**(dble(kxx-1)/dble(kxxmax-1))
       fptn(kxx,kqq,1,nint(t),ncstable2)=fpartone(ii,xx,qq,j,je,ji)
c          print *,ii,xx,qq,j,je,ji,fpartone(ii,xx,qq,j,je,ji),t
          enddo
         enddo

         endif

        enddo
       enddo
      enddo
      endif

      !2part tab between 0.1 and 1.
      !xx linear scale, qq logarithmic scale,
      !tabfptn0(j,je,ji)>0,determine case where f>0

      xxmax=1.d0
      do j=-5,2
       do je=0,2
         do ji=0,0
           t= tabfptn0(j,je,ji,ncstable2)
           if(t.gt.0.)then

             do kqq=1,kqqmax
               qq=qqmin*(qqmax/qqmin)**((kqq-1.)/(kqqmax-1.))
               do kxx=1,kxxmax
                 xx=xxmin+(xxmax-xxmin)*((kxx-1.)/(kxxmax-1.))
           fptn(kxx,kqq,2,nint(t),ncstable2)=fpartone(ii,xx,qq,j,je,ji)
               enddo
             enddo
             
           endif
         enddo
       enddo
      enddo

      write (ifmt,'(a,$)')'done)'

      end


c------------------------------------------------------------------------
      double precision function EsaturTil(iii,xx,xpe,qq,j,kk,ji)
c   former pspdf0 (sha)
c-----------------------------------------------------------------------
c
c     parton density below Q2s - Esatur(xx,xpe)=xx**(-1-dels)*EsaturTil(xx,xpe)
c     convolution of Esoft(q2sft) with E_QCD(q2st,qq)
c
c iii = proj or targ
c xx = light cone momentum fraction
c qq = virtuality scale
c xpe=energy fraction taken from remnant (by q_val+aq pair) (for val only)
c k= iclpro or icltar
c j = flavor 
c          <0 ... antiquarks
c          0 .... gluon
c          >0 ... quarks  
c ji = parton type
c         -1 sea : initial sea + evolution from sea and val without SoftSat
c          1 sea : initial sea + evolution from sea and val
c          2 val : initial val + evolution from val only      
c-----------------------------------------------------------------------
      common/ar3/    x1(7),a1(7)
      common /ar14/ x14(14),a14(14)
      double precision dgx1,dga1
      common /dga20/ dgx1(10),dga1(10)
#include "aaa.h"
#include "sem.h"
      double precision fz,zx,z,akns,dpd1,dpd2,xpe,SoftSat
     *,psuds,xx,xmin,xm,dnoflav,epscutSud2
      double precision EsoftGluonTil,EsoftQuarkTil,EsoftValTil
      integer ioTestFact
      common /ctestfact/ ioTestFact
      integer laddTestFact
      common /claddtestfact/ laddTestFact
      integer iSwitchSides
      common /cSwitchSides/iSwitchSides
      double precision sef(-5:5),aksf(2),qua,sef0(-5:5)
      integer modeDeltaFzero
      common /cmodeDeltaFzero/ modeDeltaFzero

      ii=iii
      if(iSwitchSides.eq.1)ii=3-iii

      nf=noflav(min(qq,q2sft))

      EsaturTil=0d0

      if(xx.ge.1d0)return

      nfqq=noflav(qq)
      dnoflav=dble(nfqq)

c      if(abs(j).eq.0)print *,'esatur',iii,xx,qq,j,kk,ji,nf,nfqq

c ...... f_0 * sudakov.........

      i=j 
      if(abs(i).le.nf)then
        if(i.eq.0)then
          if(ji.ne.2)sef(i)=EsoftGluonTil(xx,min(qq,q2sft),ii)
        elseif((i.eq.1.or.i.eq.2).and.ji.eq.2)then
          sef(i)=EsoftValTil(ii,xx,xpe,min(qq,q2sft),i,kk)
        elseif(abs(ji).eq.1)then
          sef(i)=EsoftQuarkTil(xx,min(qq,q2sft),ii,i)
        else
          print *,iii,xx,xpe,qq,j,kk,ji
          stop "wrong EsaturTil arguments (1) !"
        endif
        EsaturTil=sef(i)
        if(qq.le.q2sft)goto 999
        EsaturTil=EsaturTil*psuds(qq,j)/psuds(q2sft,j)
      elseif(qq.le.q2sft)then
        goto 999
      endif

c......... integral f_0 E_qcd............

      xmin=xx/(1.d0-epscutSud2(dble(qq)))

 !integrations \int d(zx) { 1/zx * zx**(-1-dels)*Esofttil(zx) * E_QCD(xx/zx) } * xx**(1+dels)   
 
      if(xmin.lt.1.d0)then

        dpd1=0.d0
        dpd2=0.d0
        xm=1d0 !max(xmin,0.1d0)
 
 !numerical integration xm -> 1

c        do i=1,7
c        do m=1,2
c          zx=1.d0-(1.d0-xm)*(.5d0+(dble(m)-1.5d0)*dble(x1(i)))!**.25d0
c          z=xx/zx
c          do n=-5,5
c            sef(n)=0
c            sef0(n)=0
c          enddo
c          if(ji.ne.2)then
c            sef(0)=zx**dble(-dels(ii))*EsoftGluonTil(zx,q2sft,ii)
c            sef(1)=zx**dble(-dels(ii))*EsoftValTil(ii,zx,xpe,q2sft,1,kk)
c            sef(2)=zx**dble(-dels(ii))*EsoftValTil(ii,zx,xpe,q2sft,2,kk)
c            sef0(1)=sef(1)
c            sef0(2)=sef(2)
c            do n=-nf,nf
c              if(n.ne.0)then !no sea contribution for gluons or if we want valence only
c                qua=zx**dble(-dels(ii))*EsoftQuarkTil(zx,q2sft,ii,n)
c                sef(n)=sef(n)+qua
c              endif
c            enddo
c          else                  !in case we want valence quark contribution only, only valence from defined flavor j is used
c            sef(j)=zx**dble(-dels(ii))*EsoftValTil(ii,zx,xpe,q2sft,j,kk)
c          endif
c          fz=0.d0
c          if(j.eq.0)then
c            if(ji.ne.2)then
c              do n=-nf,nf  
c                k=min(abs(n),1)+1
c                fz=fz+sef(n)*psevi(q2sft,qq,z,k,1)
c              enddo
c            endif
c          else
c            akns=psevi(q2sft,qq,z,3,2) !nonsinglet contribution
c            fz=max(0d0,sef(j)-sef0(j))*akns 
c            if(ji.ne.2)then  !no singlet contribution for val only
c              aksf(2)=(psevi(q2sft,qq,z,2,2)-akns) /dnoflav/2.d0 !singlet contribution
c              aksf(1)=psevi(q2sft,qq,z,1,2) /dnoflav/2.d0
c              do n=-nf,nf
c                k=min(abs(n),1)+1
c                fz=fz+sef(n)*aksf(k)
c              enddo
c            endif
c          endif
c          dpd1=dpd1 + dble(a1(i)) * fz
c     *                          / zx**2  
cc     *                   /(1.d0-zx)**3 !Jacob
c        enddo
c        enddo
c        dpd1=dpd1 * 0.5d0 * (1.d0-xm)!**4 / 4.d0  !Jacob
c     *                      * xx**dble(1d0+dels(ii))   ! <---------------

 !numerical integration  xmin -> xm

        if(xm.gt.xmin)then

          do i=1,10
          do m=1,2
            zx=xx+(xm-xx)
     &         *((xmin-xx)/(xm-xx))**(.5d0-(dble(m)-1.5d0)*dgx1(i))
            z=xx/zx
            do n=-5,5
              sef(n)=0
              sef0(n)=0
            enddo
            if(ji.ne.2)then
              sef(0)=EsoftGluonTil(zx,q2sft,ii)
              sef(1)=EsoftValTil(ii,zx,xpe,q2sft,1,kk)
              sef(2)=EsoftValTil(ii,zx,xpe,q2sft,2,kk)
              sef0(1)=sef(1)
              sef0(2)=sef(2)
              do n=-nf,nf
                if(n.ne.0)then  !no sea contribution for gluons or if we want valence only
                  qua=EsoftQuarkTil(zx,q2sft,ii,n)
                  sef(n)=sef(n)+qua
                endif
              enddo
            else                !in case we want valence quark contribution only, only valence from defined flavor j is used
            sef(j)=EsoftValTil(ii,zx,xpe,q2sft,j,kk)
            endif
            fz=0.d0
            if(j.eq.0)then
              if(ji.ne.2)then
                do n=-nf,nf
                  k=min(abs(n),1)+1
                  fz=fz+sef(n)*psevi(q2sft,qq,z,k,1)
                enddo
              endif
            else
              akns=psevi(q2sft,qq,z,3,2) !nonsinglet contribution
              fz=max(0d0,sef(j)-sef0(j))*akns 
              if(ji.ne.2)then  !no singlet contribution for val only
                aksf(2)=(psevi(q2sft,qq,z,2,2)-akns) /dnoflav/2.d0 !singlet contribution
                aksf(1)=psevi(q2sft,qq,z,1,2) /dnoflav/2.d0
                do n=-nf,nf  
                  k=min(abs(n),1)+1
                  fz=fz+sef(n)*aksf(k)
                enddo
              endif
            endif
            dpd2=dpd2+dga1(i)*fz
     *                         / zx**(2+dels(ii))   ! <-------------------------
     *                      * (zx-xx) !Jacob
          enddo
          enddo
          dpd2=dpd2 * log((xm-xx)/(xmin-xx)) * .5d0 !Jabob
     *                      * xx**dble(1.+dels(ii))   ! <---------------------
        endif
        EsaturTil=EsaturTil+dpd2+dpd1

      endif

  999 continue
      !correction below q2cmin
      if(abs(ji).eq.1)then
        if(qq.gt.q2sft)EsaturTil=EsaturTil
     .               *(1.d0+0.15d0/(1d0+0.08d0*dble(log(qq/q2sft))**2)
     .               *(1d0-dble(exp(0.22*(q2sft-qq)))))  !correction factor ???
        if(ji.gt.0)EsaturTil=EsaturTil*SoftSat(xx,qq,j,ii)  !SoftSat contains q2cmin
      endif

c     if(abs(j).eq.0)print *,'esaturtil',EsaturTil,SoftSat(xx,qq,j,ii)
c     .               (1.d0+0.15d0/(1d0+0.08d0*dble(log(qq/q2sft))**2)
c     .               *(1d0-dble(exp(0.22*(q2sft-qq)))))  !correction factor ???

      Esat=EsaturTil    ! NaN catch
      if(.not.(Esat.le.0..or.Esat.ge.0.))then
        print *,EsaturTil,iii,xx,xpe,qq,j,kk,ji
        stop'ERROR 09122019a'
      endif
      
      return
      end

c-----------------------------------------------------------------------
      double precision function EsaturTilpi(iii,xx,qq,j)
c polynomial interpolation of EsaturTil
c-----------------------------------------------------------------------
#include "sem.h"
#include "aaa.h"
#include "tab.h"
      double precision fstr,x2tmin,x2tmax,SoftSat
      parameter(kxxmax=40,kqqmax=40)  !kqqmax=80 for nip=2
      common/tabfstr/fstr(kxxmax,kqqmax,-5:5,2,2),x2tmin,x2tmax,q2tmin
     *               ,q2tmax,q2tmx
      double precision wi(3),wj(3),xxk,qqk,xx,xxmin,xxmax,EsaturTil
      data ncstable2/0/
      save ncstable2

      nip=3 !2 !3  !two or three point interpolation

      if(.false..or.qq.le.q2sft.or.qq.gt.q2mnval(maxq2mn))then

        if(iii.eq.1)then
          kk=iclpro
        else
          kk=icltar
        endif
        EsaturTilpi=EsaturTil(iii,xx,0.1d0,qq,j,kk,-1)
        
      else
        
      
      if(ncstable2.eq.0)then
        ncstable2=ncstable2+1
        call MakeEsaturTilTable
      endif
      
      qqmin=q2tmin
      qqmax=q2tmax
      qqk=1.d0+dble(log(qq/qqmin)/log(qqmax/qqmin)*(kqqmax-1))
      kqq=int(qqk)
      kqq=max(kqq,1)
      kqq=min(kqq,kqqmax-nip+1)
      wj(2)=qqk-dble(kqq)
      if(nip.eq.2)then
        wj(1)=1d0-wj(2)
      else
        wj(3)=wj(2)*(wj(2)-1.d0)*.5d0
        wj(1)=1.d0-wj(2)+wj(3)
        wj(2)=wj(2)-2.d0*wj(3)
      endif

      EsaturTilpi=0.d0
      xxmin=x2tmin

      xxmax=x2tmax
      if(xx.lt.xxmax.and.xxmin.lt.xxmax)then

        ix=1
        xxk=1.d0+log(xx/xxmin)/log(xxmax/xxmin)*(kxxmax-1)

      else

        xxmax=1.d0

        ix=2
        xxk=1.d0+(xx-xxmin)/(xxmax-xxmin)*dble(kxxmax-1)

      endif

      kxx=int(xxk)
      kxx=max(kxx,1)
      kxx=min(kxx,kxxmax-nip+1)
      
      wi(2)=xxk-dble(kxx)
      if(nip.eq.2)then
        wi(1)=1d0-wi(2)
      else
        wi(3)=wi(2)*(wi(2)-1.d0)*.5d0
        wi(1)=1.d0-wi(2)+wi(3)
        wi(2)=wi(2)-2.d0*wi(3)
      endif
      
      do kq=1,nip
        do kx=1,nip
          EsaturTilpi=EsaturTilpi+fstr(kxx+kx-1,kqq+kq-1,j,iii,ix)
     *              *wi(kx)*wj(kq)
        enddo
      enddo

      EsaturTilpi=max(0d0,(1.d0+0.2d0/dble(nip**2))*EsaturTilpi) !2% correction factor for fine tuning



      endif
      EsaturTilpi=EsaturTilpi*SoftSat(xx,qq,j,iii) !SoftSat contains q2cmin so it should not be included in the tabulation
c     print*,'Esatpi',xx,kxx,qq,kqq,j,EsaturTilpi

      Estr=EsaturTilpi     ! NaN catch
      if(.not.(Estr.le.0..or.Estr.ge.0.))then
        print*,'Esatpi',xx,qq,j,EsaturTilpi
        stop'\n\n ERROR 10012020\n\n'
      endif

      end

c-----------------------------------------------------------------------
      subroutine MakeEsaturTilTable
c-----------------------------------------------------------------------

#include "sem.h"
#include "aaa.h"
#include "tab.h"
      double precision EsaturTil
      double precision fstr,x2tmin,x2tmax
      parameter(kxxmax=40,kqqmax=40)
      common/tabfstr/fstr(kxxmax,kqqmax,-5:5,2,2),x2tmin,x2tmax,q2tmin
     *               ,q2tmax,q2tmx
      double precision xx,xxmin,xxmax

      x2tmax=0.1d0
      if(abs(noebin).gt.1)then  !special case to plot F2 for instance
        if(engmax.le.0.)
     *    stop '/n/n Please define engmax to run EsaturTilpi .../n/n'
          x2tmin=1.d0/dble(engmax**2)
      else
        x2tmin=1.d0/dble(engy**2)
      endif
      
      q2tmin=q2sft
      q2tmax=q2mnval(maxq2mn)

      write (ifmt,'(a,1x,2g8.2,1x,$)')'EsaturTable ',q2tmin,q2tmax


      !first part: tab for extrapolate x<0.1, fill between 0 and 1
      !xx,qq logarithmic scale

      qqmin=q2tmin
      qqmax=q2tmax
      xxmin=x2tmin
      xxmax=x2tmax
      if(xxmin.lt.xxmax)then
        do ii=1,2
          if(ii.eq.1)then
            kk=iclpro
          else
            kk=icltar
          endif
          if(iclpro.eq.icltar.and.ii.eq.2)then
            do j=-5,5
              do kqq=1,kqqmax
                do kxx=1,kxxmax
                  fstr(kxx,kqq,j,ii,1)=fstr(kxx,kqq,j,1,1)
                enddo
              enddo
            enddo
          else
            do j=-5,5
              write(ifmt,'(a,$)')'.'
              do kqq=1,kqqmax
                qq=qqmin*(qqmax/qqmin)**((kqq-1.)/(kqqmax-1.))
                do kxx=1,kxxmax
                  xx=xxmin*(xxmax/xxmin)**(dble(kxx-1)/dble(kxxmax-1))
       fstr(kxx,kqq,j,ii,1)=EsaturTil(ii,xx,0.05d0,qq,j,kk,-1)
c          print *,ii,xx,qq,j,EsaturTil(iii,xx,0d0,qq,j,kk,-1)
                enddo
              enddo
            enddo

         endif

       enddo
      endif

      !2part tab between 0.1 and 1.
      !xx linear scale, qq logarithmic scale,

      xxmax=1.d0
      do ii=1,2
        if(ii.eq.1)then
          kk=iclpro
        else
          kk=icltar
        endif
        if(iclpro.eq.icltar.and.ii.eq.2)then
          do j=-5,5
            do kqq=1,kqqmax
              do kxx=1,kxxmax
                fstr(kxx,kqq,j,ii,2)=fstr(kxx,kqq,j,1,2)
              enddo
            enddo
          enddo
        else
          do j=-5,5
            write(ifmt,'(a,$)')'.'
            do kqq=1,kqqmax
              qq=qqmin*(qqmax/qqmin)**((kqq-1.)/(kqqmax-1.))
              do kxx=1,kxxmax
                xx=xxmin+(xxmax-xxmin)*((kxx-1.)/(kxxmax-1.))
       fstr(kxx,kqq,j,ii,2)=EsaturTil(ii,xx,0.5d0,qq,j,kk,-1)
c          print *,ii,xx,qq,j,EsaturTil(iii,xx,0d0,qq,j,kk,-1)
              enddo
            enddo
          enddo

        endif

      enddo

      write (ifmt,'(a)')'done'

      end



c########################################################################
c########################################################################
c########################################################################
c
c                fzero   
c
c########################################################################
c########################################################################
c########################################################################


c------------------------------------------------------------------------
      function fzerotest(z)          !temporary test
c------------------------------------------------------------------------
#include "aaa.h"
#include "sem.h"
      double precision Fsea,xp,fzeroval

      qq=q2nmin
      ipt=1           !always projectile
      xp=min(1d0,dble(z))**(1d0/dble(1.+betff(ipt)))
      fzerotest=(xp*Fsea(xp,ipt)+fzeroVal(xp,qq,ipt,1)
     .         +fzeroVal(xp,qq,ipt,2))*xp**(-betff(ipt))         !jacobian
      end

c------------------------------------------------------------------------
      subroutine fzeroQ2min(s,z,ipt,q2mean)
c-----------------------------------------------------------------------
c
c     average Q2sft for PDF calculation = \int Q2*x.f(x) / \int x.f(x)
c
c s - center of mass energy square
c z - light cone x
c ipt - 1=proj 2=targ
c-----------------------------------------------------------------------
      double precision z!,xpmin,xp,zz,Fsea,EsoftGluonTil,Fval,EsoftValTil
c     .,EsoftQuarkTil,EsoftTil,fnorm,fzeroQ2
#include "aaa.h"
#include "sem.h"
#include "par.h"
      common /ar3/   x1(7),a1(7)


      q2pmin(ipt)=q2nmin
      q2mean=q2pmin(ipt)
      dum=s
      dum=z

      return

c      if(iscreen.ne.0)then
c old procedure removed 3.247
c  ...

      end

c------------------------------------------------------------------------
      double precision function fzeroVal(z,qq,ipt,j)
c-----------------------------------------------------------------------
c
c        fzeroVal = z*f0(z)
c
c   f0 = Fval & Esatur      &=convolution
c   f0=\int_z^1 dxp/xp Fval(xp) Esatur(z/xp,xp,qq)
c
c   Fval(xp) = (1-xp)**alplea(k)
c
c   Esatur(zz,xp) = EsaturValTil(zz,xp,qq) * zz**(-1-dels)
c
c z - light cone momentum fraction of parton
c ipt - 1=proj 2=targ
c Integration only OK for 0 < alppar < 0.5, otherwise numerical integration
c is not good enough (gauss with singularity around z). If alppar is larger
c than 0.5, psdpdf (=fzeroval by construction) is used.
c (betff=-alppar)      
c-----------------------------------------------------------------------
      double precision xpmin,z,xp,zz,Fval,EsaturValTil
#include "aaa.h"
#include "sem.h"
      double precision dgx1,dga1
      common /dga20/   dgx1(10),dga1(10)
      integer ioTestFact
      common /ctestfact/ ioTestFact
      integer idhprTestFact
      common /cidhprtestfact/ idhprTestFact

      ii = ipt
      if(ipt.eq.1)then
        k=iclpro
      else
        k=icltar
      endif
      fzeroVal=0

      !kkkkkkkkkkkkkkkkkkkkk!KW1908
      if(ioTestFact.ne.0)then
        if(ioTestFact.eq.-1)then !val-val
          goto 999
        elseif(ioTestFact.eq.-2)then !sea-val
          if(ipt.eq.1)return
          goto 999
        elseif(ioTestFact.le.-3)then !val-sea
          if(ipt.eq.2)return
          goto 999
        elseif(ioTestFact.ge.3)then !sea-sea
          return
        endif
        stop'wrong ioTestFact (1)'
      endif
      if(idhprTestFact.ne.999)then
        if(idhprTestFact.eq.0)return !sea-sea
        if(idhprTestFact.eq.1.and.ipt.eq.2)return !val-sea     
        if(idhprTestFact.eq.2.and.ipt.eq.1)return !sea-val
      endif
      !kkkkkkkkkkkkkkkkkkkkk


      !test
      !fzeroVal=psdpdf(sngl(z),min(qq,q2cmin(ii)),0.,k,j)
      !return
      !end test

      fzeroval=0.d0

      !if(-betff(ii).le.0.51)then  !numerical integration doesn't work for large alppar
      xp=0d0
      alpq=alpqua(ipt)

      xpmin=z**dble(-alpq)      !not need for dels here since it's cancelled in Esaturvaltil
      do i=1,10
      do m=1,2
       xp=(.5d0*(1d0+xpmin+dble(2*m-3)*dgx1(i)*(1d0-xpmin)))
     *         **(1d0/dble(-alpq))
       zz=min(0.9999999999d0,z/xp)
       fzeroval=fzeroval+dga1(i)*Fval(xp,ii)
     *   *xp**dble(alpq)
     *      *zz**dble(-1.-dels(ipt))*EsaturValTil(ii,zz,xp,qq,j,k)
      !---------------------------------------
      !*   /xp
      !*   *xp/xp**dble(-alpq)    !Jacobian
      !*   *  *zz**dble(-1.-dels(ipt))*EsaturValTil(ii,zz,xp,qq,j,k)
      !------------------------------------------
      enddo
      enddo
      fzeroval=fzeroval*.5d0*(1d0-xpmin)/dble(-alpq)
     .     *z                   !z*f(z)


c      print *,z,qq,fzeroval,fzeroval/psdpdf(sngl(z),qq,0.,k,j)
      
      return

      !kkkkkkkkkkkkkkkkkkkkk!KW1910
 999  continue !using modeDeltaFzero=0, but F=delta(x-1)
      fzeroVal=EsaturValTil(ipt,z,1d0,qq,j,k)*z**dble(-dels(ipt))
      return
      !kkkkkkkkkkkkkkkkkkkkk

      end

c------------------------------------------------------------------------
      double precision function FzeroQua(z,qq,ipt,iflav) !per flavor
c-----------------------------------------------------------------------
c
c   FzeroQua = z*f(z)
c
c   f = F & EsaturQuark         &=convolution
c
c   F(x) = function Fsea(...)*Fgfactor(...)
c
c   EsaturQuark(zz) = zz**(-1.-dels) * EsaturQuarkTil (zz)
c
c z - light cone x of the quark,
c ipt - 1=proj 2=targ
c-----------------------------------------------------------------------
#include "aaa.h"
#include "sem.h"
      double precision xpmin,xp,z,zz,Fsea,Fgfactor,EsaturQuarkTil
     .,zzk,valueK
      common /ar3/   x1(7),a1(7)
      integer ioTestFact
      common /ctestfact/ ioTestFact
      integer idhprTestFact
      common /cidhprtestfact/ idhprTestFact
      data zzK/0d0/ 
      data qqK /0./ 
      data iptK /0/ 
      save zzK,qqK,iptK,valueK

      !----------tp test---------------
      !if(ipt.eq.1)then
      ! k=iclpro
      !else
      ! k=icltar
      !endif
      !FzeroQua=psdpdf(sngl(z),qq,0.,k,-1)*dble(naflav)*2.d0
      !return
      !---------------------------------

      if(.not.(iflav.ge.1.and.iflav.le.5)
     ..and..not.(-iflav.ge.1.and.-iflav.le.5))stop'ERROR 06092019d'

      FzeroQua=0.d0

      nf=noflav0(qq)
      if(abs(iflav).gt.nf)return     !like in EsaturQuark

      !kkkkkkkkkkkkkkkkkkkkk!KW1908
      if(ioTestFact.ne.0)then
        if(ioTestFact.eq.-1)then !val-val
          return
        elseif(ioTestFact.eq.-2)then !sea-val
          if(ipt.eq.2)return
          goto 999
        elseif(ioTestFact.eq.-5.or. !val-sea,qq 
     .    ioTestFact.eq.-3 .or.
     .    ioTestFact.eq.-6 .or.ioTestFact.eq.-7)then !val-sea
          if(ipt.eq.1)return
          goto 999
        elseif(ioTestFact.eq.-4)then !val-sea,qg
          return
        elseif(ioTestFact.ge.3.and.ioTestFact.le.5)then
          goto 999
        elseif(ioTestFact.ge.6.and.ioTestFact.le.11)then
          FzeroQua=EsaturQuarkTil(z,qq,ipt,iflav) !being 1 or 0, using modeDeltaFzero=1 (fzero...=delta(x-1))
          return
        elseif(ioTestFact.eq.12)then
          goto 999
c          FzeroQua=EsaturQuarkTil(z,qq,ipt,iflav) !being 1 or 0, using modeDeltaFzero=1 (fzero...=delta(x-1))
c          return
        endif
        stop'wrong ioTestFact (3)'
      endif
      if(idhprTestFact.ne.999)then
        if(idhprTestFact.eq.1.and.ipt.eq.1)return !val-sea     
        if(idhprTestFact.eq.2.and.ipt.eq.2)return !sea-val
        if(idhprTestFact.eq.3)return !val-val
      endif
      !kkkkkkkkkkkkkkkkkkkkk

      if(z.eq.zzK.and.qq.eq.qqK.and.ipt.eq.iptK.and.abs(iflav).le.3)then
        FzeroQua=valueK
c        print*,'FzeroQua',iflav,z,zzK,qq,qqK,ipt,iptK,FzeroQua
        return
      endif

      xpmin=z
      xpmax=1
      xpmin=xpmin**dble(1+betff(ipt)+dels(ipt))
      xpmax=xpmax**dble(1+betff(ipt)+dels(ipt))
      if(xpmin.lt.xpmax)then
      do i=1,7
      do m=1,2
        xp=(.5d0*(xpmax+xpmin+dble(2*m-3)*dble(x1(i))*(xpmax-xpmin)))
     *   **(1d0/dble(1+betff(ipt)+dels(ipt)))
        zz=min(0.9999999999d0,z/xp)
        !-----------------------------------------------------------
        ! Integrand:  (several terms cancel)
        !   F(xp) / xp * zz**dble(-1.-dels) * EsaturQuarkTil (zz)
        !                / xp**dble(betff(ipt)+dels)               !Jacob(1)
        !-----------------------------------------------------------
        FzeroQua=FzeroQua+dble(a1(i))*Fsea(xp,ipt)
     *  *Fgfactor(1d0,xp,1)*EsaturQuarkTil(zz,qq,ipt,iflav)
     *  *xp**dble(-betff(ipt))
      enddo
      enddo
      endif
      FzeroQua=FzeroQua
     *   *.5d0*(xpmax-xpmin)/dble(1.+betff(ipt)+dels(ipt))        !Jacob(2)
     *   *z**dble(-dels(ipt))                            ! z**(-1-dels(ipt)) * z
      fz=FzeroQua     ! NaN catch
      if(.not.(fz.le.0..or.fz.ge.0.))then
        print *,'FzeroQua',fz,z,qq,ipt,iflav
        stop'ERROR 25102012d'
      endif
      
      if(FzeroQua.gt.0d0.and.abs(iflav).le.3)then   !different for HQ because of Esatur and not Esoft !
        zzK = z
        qqK = qq 
        iptK = ipt
        valueK = FzeroQua
      else
        zzK = 0d0
        qqK = 0.
        iptK = 0
      endif

      return

      !kkkkkkkkkkkkkkkkkkkkk!KW1910
 999  continue !using modeDeltaFzero=0, but F=delta(x-1)
      FzeroQua=EsaturQuarkTil(z,qq,ipt,iflav)*z**dble(-dels(ipt)) 
      return
      !kkkkkkkkkkkkkkkkkkkkk

      end

c------------------------------------------------------------------------
      double precision function fzeroGlu(z,qq,ipt)
c-----------------------------------------------------------------------
c
c   fzeroGlu = z*f(z)
c
c   f = F & EsaturGluon         &=convolution
c
c   F(x) = Fsea(x) * Fgfactor(x)
c
c   EsaturGluon(x) = x**(-1.-dels) * EsaturGluonTil(x)
c
c z - light cone x
c ipt - 1=proj 2=targ
c-----------------------------------------------------------------------
#include "aaa.h"
#include "sem.h"
      double precision xpmin,xp,z,zz,Fsea,Fgfactor,EsaturGluonTil
      common /ar3/   x1(7),a1(7)
      external fzerotest
      integer ioTestFact
      common /ctestfact/ ioTestFact
      integer idhprTestFact
      common /cidhprtestfact/ idhprTestFact

      fzeroGlu=0

      !kkkkkkkkkkkkkkkkkkkkk!KW1908
      if(ioTestFact.ne.0)then
        if(ioTestFact.eq.-1)then !val-val
          return
        elseif(ioTestFact.eq.-2)then !sea-val
          if(ipt.eq.2)return
          goto 999 
        elseif(ioTestFact.eq.-5)then !val-sea,qq
          return
        elseif(ioTestFact.eq.-4 .or. !val-sea,qg
     .    ioTestFact.eq.-3 .or.
     .    ioTestFact.eq.-6 .or.ioTestFact.eq.-7)then !val-sea
          if(ipt.eq.1)return
          goto 999 
        elseif(ioTestFact.ge.3.and.ioTestFact.le.5)then !sea-sea ; sea-sea gq ; sea-sea gg 
          goto 999 
        elseif(ioTestFact.ge.6.and.ioTestFact.le.11)then
          fzeroGlu=EsaturGluonTil(z,qq,ipt)  !being 1 or 0, using modeDeltaFzero=1 (fzero...=delta(x-1))
          return
        elseif(ioTestFact.eq.12)then
          goto 999 
c          fzeroGlu=EsaturGluonTil(z,qq,ipt)  
c          return
        endif
        stop'wrong ioTestFact (5)'
      endif
      if(idhprTestFact.ne.999)then
        if(idhprTestFact.eq.1.and.ipt.eq.1)return !val-sea     
        if(idhprTestFact.eq.2.and.ipt.eq.2)return !sea-val
        if(idhprTestFact.eq.3)return !val-val
      endif
      !kkkkkkkkkkkkkkkkkkkkk
   

      !tp test----------
      !if(ipt.eq.1)then
      ! k=iclpro
      !else
      ! k=icltar
      !endif
      !fzeroGlu=psdpdf(sngl(z),q2cmin(ipt),0.,k,0)
      !return
      !-----------------

      fzeroGlu=0.d0
      xpmin=z
      xpmax=1
      xpmin=xpmin**dble(1+betff(ipt)+dels(ipt))
      xpmax=xpmax**dble(1+betff(ipt)+dels(ipt))
      if(xpmin.lt.xpmax)then
      do i=1,7
      do m=1,2
        xp=(.5d0*(xpmax+xpmin+dble(2*m-3)*dble(x1(i))*(xpmax-xpmin)))
     *    **(1d0/dble(1+betff(ipt)+dels(ipt)))
        zz=min(0.9999999999d0,z/xp)
        !-----------------------------------------------------------
        ! Integrand:  (several terms cancel)
        !  F(xp)/xp * zz**dble(-1.-dels) * EsaturGluonTil(zz)
        !               /xp**dble(betff(ipt)+dels)            !Jacob(1)
        !     
        !-----------------------------------------------------------
        fzeroGlu=fzeroGlu+dble(a1(i))*Fsea(xp,ipt)
     *  *Fgfactor(1d0,xp,1)*EsaturGluonTil(zz,qq,ipt)
     *  *xp**dble(-betff(ipt))
      enddo
      enddo
      endif
      fzeroGlu=fzeroGlu
     *       *.5d0*(xpmax-xpmin)/dble(dels(ipt)+betff(ipt)+1.)  !Jacob(2)
     *       *z**dble(-dels(ipt))                    ! z**(-1-dels(ipt)) * z

      ! test normalization (fz=1 if gamhad correctly defined) (eq. C.31 in Phys Rept)
c            xpmin=1./engy*2
c            xpmin=xpmin**dble(1.+betff(ipt))
c            call uttraq(fzerotest,sngl(xpmin),1.,fz)
c            fz=fz/(1.+betff(ipt))
c            print *,q2cmin(1),gamhad(2),engy,z,qq,betff(ipt),fz

      fz=fzeroGlu     ! NaN catch
      if(.not.(fz.le.0..or.fz.ge.0.))then
        print*, fzeroGlu,xpmin,z
        stop'ERROR 25102012c NaN catch'
      endif
      return

      !kkkkkkkkkkkkkkkkkkkkk!KW1910
 999  continue !using modeDeltaFzero=0, but F=delta(x-1)
      fzeroGlu=EsaturGluonTil(z,qq,ipt)*z**dble(-dels(ipt)) 
      return
      !kkkkkkkkkkkkkkkkkkkkk

      end



c########################################################################
c########################################################################
c########################################################################
c
c                Esoft   
c
c########################################################################
c########################################################################
c########################################################################


c-----------------------------------------------------------------------
      double precision function EsoftValTil(ii,x,xpe,qq,j,k)
c-----------------------------------------------------------------------
c  EsoftVal = x**(-1-dels)*EsoftValTil
c  x=fraction of xpe going into valence quark
c  xpe=energy fraction taken from remnant (by q_val+aq pair)
c  k= iclpro or icltar
c Exact results if direct comparison between fzeroval and psdpdf
c----------------------------------------------------------------------
#include "aaa.h"
#include "sem.h"
      double precision x,xpe,xq,alpq,delsd,coefxu1!,utgam2
      common /psar40/ coefxu1(2,nclha)

      alpq=dble(alpqua(ii))  !=0.5*(1+appar)
      delsd=dble(dels(ii))
      xq=x*xpe          !part of the valence quark
      EsoftValTil=dble(psdpdf(sngl(xq),qq,0.,k,j))  !xq*F(xq) -> xq cancelled with x**(1+delsd)
c      print *,ii,k,EsoftValTil,coefxu1(ii,k),delsd,alplea(k),alpq,xpe,x
      EsoftValTil=EsoftValTil
     *  *max(1d-20,(1d0-x)*xpe)**(-alpq)   !part for the anti-q: (1-x)*xpe
     *  *max(1d-20,1.d0-xq)**(-1d0+alpq-dble(alplea(k))) !cf PSHARF
     *  *coefxu1(ii,k)
     *  *x**(delsd)          !no soft preevolution (to compensate factor added everywhere for consistency with EsoftGluon or EsoftQuark)  = (xq/xpe)**(1+delsd)/xq*xpe
c     *  *xpe -> cancelled with x**(1+delsd)
     *  *xpe**(-betff(ii))/alpff(k)   !because multiplied by Fsea in EsaturTil instead of Fval
c     *  *utgam2(2.d0+dble(alplea(k))-alpq)  !normalization to get psdpdf
c     *  /utgam2(1.d0+dble(alplea(k)))/utgam2(1.d0-alpq) !replaced by coefxu1
     *  *sqrt(2.*alpq)       !correction factor to get correct integral
      end

c------------------------------------------------------------------------
      double precision function EsoftQuarkTil(zz,qq,ipt,iflav) !per flavor
c-----------------------------------------------------------------------
c   EsoftQuark = zz^(-1-dels) * EsoftQuarkTil
c   EsoftQuarkTil = gamsoft*glusea/(1-glusea)*int_zz^1 dz Pqg(z) EsoftGluon
c   Pqg=Nf*(x**2+(1-x)**2)   Altarelli Parisi splitting function
c-----------------------------------------------------------------------
c   iflav ... flavor considered (+-i, i=1,2,3,4,5) !not used presently
c             999 means any flavor (all considered equal)
c             but in any case: we consider quark distribution per flavor 
c             (we use Pqg/Nf=(x**2+(1-x)**2) directly)
c-----------------------------------------------------------------------
#include "aaa.h"
#include "sem.h"
      double precision zmin,z,zz,zzK,valueK!,SoftSat
      common /ar3/   x1(7),a1(7)
      integer ioTestFact
      common /ctestfact/ ioTestFact
      data zzK/0d0/ 
      data qqK /0./ 
      data iptK /0/ 
      save zzK,qqK,iptK,valueK

      if(.not.(iflav.eq.999)
     ..and..not.(iflav.ge.1.and.iflav.le.5)
     ..and..not.(-iflav.ge.1.and.-iflav.le.5))stop'ERROR 06092019b'


      EsoftQuarkTil=0.d0
      j=abs(iflav)
      if(j.eq.999)j=3
      nf=noflav(qq)
      if(j.gt.nf)return
     
      is=iSpecialQuark(ipt,iflav)
      if(is.ge.0)then
        EsoftQuarkTil=dble(is)
        return 
      endif 

      if(zz.eq.zzK.and.qq.eq.qqK.and.ipt.eq.iptK)then
        EsoftQuarkTil=valueK
c        print*,'FzeroQua',iflav,zz,zzK,qq,qqK,ipt,iptK,EsoftQuarkTil
        return
      endif

      zmin=zz
      zmin=zmin**dble(1.+dels(ipt))
      do i=1,7
      do m=1,2
        z=(.5d0*(1.d0+zmin+dble(2*m-3)*dble(x1(i))*(1.d0-zmin)))
     *  **(1.d0/dble(1.+dels(ipt)))
        !-----------------------------------------------------------
        ! Integrand:  (several terms cancel)
        ! (1-zz/z))**dble(betpom) *(zz/z)**dble(-1.-dels) !from EsoftGluon
        !           *(z**2+(1.d0-z)**2)/z                 !Pqg
        !           *z/z**dble(1.+dels)                   !Jacobian
        !------------------------------------------------------------
        EsoftQuarkTil=EsoftQuarkTil+dble(a1(i))
     *   *max(1.d-15,(1.d0-zz/z))**dble(betpom(ipt)) !from EsoftGluon
c     *   *SoftSat(zz/z,qq,iflav,ipt)
     *                           *0.5*(z**2+(1.d0-z)**2)
      enddo
      enddo
      EsoftQuarkTil=EsoftQuarkTil*(1.d0-zmin)/dble(1.+dels(ipt))

      !factor zz**dble(1.+dels)cancels against term in integrand

      EsoftQuarkTil=dble(gamsoft(ipt)*fglusea(zz,qq,ipt))*EsoftQuarkTil

      if(EsoftQuarkTil.gt.0d0)then
        zzK = z
        qqK = qq 
        iptK = ipt
        valueK = EsoftQuarkTil
      else
        zzK = 0d0
        qqK = 0.
        iptK = 0
      endif

      end

c-----------------------------------------------------------------------
      double precision function EsoftGluonTil(zz,qq,ipt)
c-----------------------------------------------------------------------
c   EsoftGluon = zz^(-1.-dels) * EsoftGluonTil
c-----------------------------------------------------------------------
#include "aaa.h"
#include "sem.h"
      double precision zz!,SoftSat
      integer ioTestFact
      common /ctestfact/ ioTestFact

      EsoftGluonTil=0

      is=iSpecialGluon(ipt)
      if(is.ge.0)then
        EsoftGluonTil=dble(is)
        return
      endif

      EsoftGluonTil=dble(gamsoft(ipt)*(1.-fglusea(zz,qq,ipt)))
     .                     * (1.d0-zz)**dble(betpom(ipt))
c     .                     * SoftSat(zz,qq,0,ipt)
       end

c-----------------------------------------------------------------------
      integer function iSpecialQuark(ipt,iflav) !KW20
c-----------------------------------------------------------------------
c Checks for special test cases
c Returns: 0 or 1 for special cases
c           -1    for normal case, real calculation required
c-----------------------------------------------------------------------
      integer ioTestFact
      common /ctestfact/ ioTestFact
      iSpecialQuark=0
      !kkkkkkkkkkkkkkkkkkkkk!KW1908
      if(ioTestFact.ne.0)then
      if(ioTestFact.eq.-1)then !val-val
        return
      elseif(ioTestFact.eq.-2)then !sea-val
        if(ipt.eq.2)return  
        goto 999 
      elseif(ioTestFact.eq.-5.or. !val-sea,qq
     .  ioTestFact.eq.-3 .or.ioTestFact.eq.-6 .or.ioTestFact.eq.-7)then !val-sea
        if(ipt.eq.1)return  
        goto 999 
      elseif(ioTestFact.eq.-4)then !val-sea,qg
        return  
      elseif(ioTestFact.eq.3)then
        goto 999 
      elseif(ioTestFact.eq.4)then
        if(ipt.eq.1)return
        if(iflav.eq.1.or.iflav.eq.999)goto 999 !ipt=2
        return 
      elseif(ioTestFact.eq.5)then
        return  
      elseif(ioTestFact.ge.6.and.ioTestFact.le.12)then
        if(ioTestFact.eq.6)then
          if(ipt.eq.2.and.(iflav.eq.1.or.iflav.eq.999))iSpecialQuark=1
        elseif(ioTestFact.eq.10)then
          if(ipt.eq.1.and.(iflav.eq.1.or.iflav.eq.999))iSpecialQuark=1
        elseif(ioTestFact.eq.11)then
          if((iflav.eq.1.or.iflav.eq.999))iSpecialQuark=1
        elseif(ioTestFact.eq.12)then
c          return  
          goto 999
        endif
        !more detailed specificaton : calcSTYP
        return  
      endif
      stop'wrong ioTestFact (4)'
      endif
      !kkkkkkkkkkkkkkkkkkkkk
 999  continue !normal case
      iSpecialQuark=-1
      end

c-----------------------------------------------------------------------
      integer function iSpecialGluon(ipt) !KW20
c-----------------------------------------------------------------------
c Checks for special test cases
c Returns: 0 or 1 for special cases
c           -1    for normal case, real calculation required
c-----------------------------------------------------------------------
      integer ioTestFact
      common /ctestfact/ ioTestFact
      iSpecialGluon=0
      !kkkkkkkkkkkkkkkkkkkkk!KW1908
      if(ioTestFact.ne.0)then
      if(ioTestFact.eq.-1)then !val-val
        return  
      elseif(ioTestFact.eq.-2)then !sea-val
        if(ipt.eq.2)return  
        goto 999 
      elseif(ioTestFact.eq.-5)then !val-sea,qq
        return
      elseif(ioTestFact.eq.-4 .or. !val-sea,qg
     .  ioTestFact.eq.-3 .or.ioTestFact.eq.-6  .or.ioTestFact.eq.-7)then !val-sea
        if(ipt.eq.1)return  
        goto 999 
      elseif(ioTestFact.eq.3)then !sea-sea 
        goto 999  
      elseif(ioTestFact.eq.4)then !sea-sea gq
        if(ipt.eq.1)goto 999  
        if(ipt.eq.2)return
      elseif(ioTestFact.eq.5)then !sea-sea gg
        goto 999  
      elseif(ioTestFact.ge.6.and.ioTestFact.le.12)then
        if(ioTestFact.eq.6)then
          if(ipt.eq.1)iSpecialGluon=1
        elseif(ioTestFact.ge.7.and.ioTestFact.le.9)then
          iSpecialGluon=1
        elseif(ioTestFact.eq.10)then
          if(ipt.eq.2)iSpecialGluon=1
        elseif(ioTestFact.eq.12)then
          goto 999  
c          iSpecialGluon=1
        endif
        return
      endif
      stop'wrong ioTestFact (6)'
      endif
      !kkkkkkkkkkkkkkkkkkkkk
  999 iSpecialGluon=-1
      end



c########################################################################
c########################################################################
c########################################################################
c
c                Esat   
c
c########################################################################
c########################################################################
c########################################################################

c-----------------------------------------------------------------------
      double precision function EsatValTil(ii,x,xpe,qq,k,ivo,uv1,dv1)
c-----------------------------------------------------------------------
c ii : 1 = proj   2 = targ 
c-----------------------------------------------------------------------
c iremn via aaa.h
#include "aaa.h"
      common /cnquud12/ nquu1,nqud1,nquu2,nqud2
      double precision x,xpe,uv1,dv1
      double precision EsaturValTil
      if(ii.eq.1)then
        nquu=nquu1
        nqud=nqud1
      else
        nquu=nquu2
        nqud=nqud2
      endif 
      if(iremn.ge.2.and.ivo.eq.1)then
        if(nquu.gt.nqud)then
          uv1=EsaturValTil(ii,x,xpe,qq,1,k)
          dv1=EsaturValTil(ii,x,xpe,qq,2,k)
        else              !if nquu<nqud => no u or no d
          uv1=EsaturValTil(ii,x,xpe,qq,1,k)
          dv1=uv1
        endif
        if(nquu.eq.0)uv1=0d0
        if(nqud.eq.0)dv1=0d0
      else
        uv1=EsaturValTil(ii,x,xpe,qq,1,k)
        dv1=EsaturValTil(ii,x,xpe,qq,2,k)
      endif
      EsatValTil=uv1+dv1
      end

c-----------------------------------------------------------------------
      double precision function EsatSeaTil(ipomtype,ipt,zz,qq,J)
c-----------------------------------------------------------------------
c ipt : 1 = proj   2 = targ 
c J : flavor type : 0=glu, 1=light, 2=charm, 3=bottom
c-----------------------------------------------------------------------
      double precision EsaturQuarkTil,EsaturGluonTil
      double precision zz
      dummy=ipomtype
      if(J.eq.0)then
        EsatSeaTil=EsaturGluonTil(zz,qq,ipt)
      elseif(J.eq.1)then 
        EsatSeaTil=EsaturQuarkTil(zz,qq,ipt,1)
      elseif(J.eq.2)then 
          EsatSeaTil=EsaturQuarkTil(zz,qq,ipt,4)
      elseif(J.eq.3)then
          EsatSeaTil= EsaturQuarkTil(zz,qq,ipt,5)
      else
        stop'ERROR 15062020a'
      endif
      end

c-----------------------------------------------------------------------
      double precision function EsaturValTil(ii,x,xpe,qq,j,k)
c-----------------------------------------------------------------------
c  EsaturVal = x**(-1-dels)*EsaturValTil
c  x=fraction of xpe going into valence quark
c  xpe=energy fraction taken from remnant (by q_val+aq pair)
c  k= iclpro or icltar
c Exact results if direct comparison between fzeroval and psdpdf
c-----------------------------------------------------------------------
#include "aaa.h"
      double precision x,xpe,EsoftValTil!,EsaturTil
      integer ioTestFact
      common /ctestfact/ ioTestFact
      
      EsaturValTil=0d0

      !kkkkkkkkkkkkkkkkkkkkk!KW1910
      if(ioTestFact.ne.0)then
      if(ioTestFact.eq.-1)then !val-val
        goto 999 
      elseif(ioTestFact.eq.-2)then !sea-val
        if(ii.eq.1)return  
        goto 999 
      elseif(ioTestFact.le.-3)then !val-sea
        if(ii.eq.2)return  
        goto 999 
      elseif(ioTestFact.ge.3)then !sea-sea
        return
      endif
      stop'wrong ioTestFact (2)'
      endif
      !kkkkkkkkkkkkkkkkkkkkk

 999  continue !compute EsaturValTil normally
c      EsaturValTil=EsaturTil(ii,x,xpe,qq,j,k,2)
c     *  *xpe**(betff(ii))*alpff(k)   !not included in Fval
      EsaturValTil=EsoftValTil(ii,x,xpe,qq,j,k)
     *  *xpe**(betff(ii))*alpff(k)   !not included in Fval
c      print *,'---->',EsaturValTil,xpe**(betff(ii))*alpff(k)
      return
      end
      
c-----------------------------------------------------------------------
      double precision function EsaturQuarkTil(zz,qq,ipt,iflav)
c-----------------------------------------------------------------------
c   EsaturQuark = zz^(-1-dels) * EsaturQuarkTil
c   iflav ... flavor considered (+-i, i=1,2,3,4,5) 
c             we consider quark distribution per flavor 
c-----------------------------------------------------------------------
#include "aaa.h"
      double precision zz,EsaturTilpi
      integer ioTestFact
      common /ctestfact/ ioTestFact
      
      if(.not.(iflav.eq.999)
     ..and..not.(iflav.ge.1.and.iflav.le.5)
     .     .and..not.(-iflav.ge.1.and.-iflav.le.5))stop'ERROR 10122019a'
      if(iflav.eq.999)then
        j=-3
      else
        j=iflav
      endif
      EsaturQuarkTil=0d0

      !allow charm down to qcmass because some contribution coming from evolution 
      !(and in good agreement with data for charm production below 2*qcmass)
      nf=noflav0(qq)
      if(abs(j).gt.nf)return

      is=iSpecialQuark(ipt,iflav)
      if(is.ge.0)then
        EsaturQuarkTil=dble(is)
        return 
      endif 


      if(ipt.eq.1)then
        k=iclpro
      else
        k=icltar
      endif
      EsaturQuarkTil=EsaturTilpi(ipt,zz,qq,j)
      

      return
      end
      
c-----------------------------------------------------------------------
      double precision function EsaturGluonTil(zz,qq,ipt)
c-----------------------------------------------------------------------
c   EsaturGluon = x^(-1.-dels) * EsaturGluonTil
c-----------------------------------------------------------------------
#include "aaa.h"
      double precision zz,EsaturTilpi!,EsaturTil
      integer ioTestFact
      common /ctestfact/ ioTestFact

      EsaturGluonTil=0

      is=iSpecialGluon(ipt)
      if(is.ge.0)then
        EsaturGluonTil=dble(is)
        return
      endif
      
      EsaturGluonTil=EsaturTilpi(ipt,zz,qq,0)

c      if(ipt.eq.1)then
c        k=iclpro
c      else
c        k=icltar
c      endif
c      EsaturGluonTil=EsaturTil(ipt,zz,0d0,qq,0,k,1)
      
      return
      end



c########################################################################
c########################################################################
c########################################################################
c
c               Vertex functions
c
c########################################################################
c########################################################################
c########################################################################

c-------------------------------------------------------------
      double precision function Fval(xpe,ipt)
c-------------------------------------------------------------
c
c    vertex function Fval
c-------------------------------------------------------------
#include "sem.h"
#include "aaa.h"
      double precision xpe

      if(ipt.eq.1)then
        k=iclpro
      else
        k=icltar
      endif

c      Fval=
c     .  dble(alpff(k))*xpe**dble(betff(ipt))*(1d0-xpe)**dble(alplea(k))
      Fval=(1d0-xpe)**dble(alplea(k))

      end


c-------------------------------------------------------------
      double precision function Fsea(xpe,ipt)
c-------------------------------------------------------------
c
c vertex function Fsea
c
c-------------------------------------------------------------
#include "sem.h"
#include "aaa.h"
      double precision xpe

      if(ipt.eq.1)then
        k=iclpro
      else
        k=icltar
      endif

c      if(iscreen.ne.0)then

c        stop'\n\n  ERROR 14112012 \n\n'

c      else

        Fsea=
     .  dble(alpff(k))*xpe**dble(betff(ipt))*(1d0-xpe)**dble(alplea(k))
c       .  dble(gamhad(k))*xpe**dble(-alppar)*(1d0-xpe)**dble(alplea(k))

c      endif

      end



c###########################################################################
c###########################################################################
c###########################################################################
c###########################################################################

c  KW1908:   some more ffsig stuff (not updated)   

c###########################################################################
c###########################################################################
c###########################################################################
c###########################################################################


c------------------------------------------------------------------------
      function ffsigiut(xx1,xx2,jpp,je1,je2)
c------------------------------------------------------------------------
c
c   \int(dt) \int(du)  ffsig *s/sh**3 *2*pi*alpha**2 *delta(uh+th+sh)
c
c   Here q2cmin defined as the minimum Q2s given by xx (x_pom > xx of the
c   parton). Real value of q2cmin averaged in fpartone. q2cmin should be
c   redefined before each ffsigiut call to avoid change by another function
c-----------------------------------------------------------------------
      common /ar3/   x1(7),a1(7)
#include "sem.h"
#include "aaa.h"
      double precision t,ffsigj,ft,ffsigiutd
      integer ienvi
      common /cienvi/ ienvi
      ienvi=101

      ig=7
      s=engy**2
      sh=s*xx1*xx2
      ffsigiut=0.
      qqcut=max(q2cmin(1),q2cmin(2))
      qq=qqcut

      ffsigiut=0
  
      do klas=1,klasmax
       ffsigiutd=0.d0
       call HardScale(klas,qq,pt2min,1)      !get pt2min
       call getBornKin(klas,sh,pt2min, smin, tmin, tmax,iret) !get tmin,tmax
       if(iret.ne.0)then
        do i=1,ig
        do m=1,2
          t=2d0*tmin/(1d0+tmin/tmax-dble(tgss(ig,i)*(2*m-3))
     &         *(1d0-tmin/tmax))
          call getBornPtr(klas,sh,sngl(t) , pt2)  !get pt2
          call HardScale(klas, scale , pt2 ,2) !get scale
          ft=ffsigj(klas,pt2,dble(xx1),dble(xx2),jpp,je1,je2)
     *        /sh**3
     *        *dble(2*pi*pssalf(scale/qcdlam))**2
          ffsigiutd=ffsigiutd+dble(wgss(ig,i))*ft*t**2
        enddo
        enddo
       ffsigiutd=ffsigiutd*0.5d0*(1d0/tmin-1d0/tmax)
       endif
       ffsigiut=ffsigiut + ffsigiutd
      enddo

      ffsigiut=ffsigiut
     *     *2*pi
     *     *s      !  /2   !CS for parton pair   factor is WRONG

      return
      end

c-----------------------------------------------------------------------
      double precision function ffsigj(klas,qt,x1,x2,jpp,je1,je2)  !KW1811 charm added
c-----------------------------------------------------------------------
c
c      \sum  x1*f_i(x1,qt) * x2*f_k(x2,qt) * B_ik
c
c        B_ik = psbori = contribution to Born xsection:
c                         dsigmaBorn/d2pt/dy
c                          = s/pi * delta(s+t+u) * 2*pi*alpha**2 /s**2 * B_ik
c  qt = pt2
c
c  x1, x2 = light cone momentum fractions
c
c  x*f_j(x,qt) = function fparton(x,qt,j)
c
c-----------------------------------------------------------------------
c jpp: type of Pomeron
c          1 ... sea-sea
c          2 ... val-sea
c          3 ... sea-val
c          4 ... val-val
c          5 ... all
c je = emission type
c          0 ... no emissions
c          1 ... emissions
c          2 ... all
c-----------------------------------------------------------------------
#include "sem.h"
#include "aaa.h"
      double precision ffborn,pifpartone,x1,x2,s,t,smin,ttt,tmax

      s=dble(engy*engy)*x1*x2
      ffsigj=0d0
      
      call HardScale(klas,sc,qt,2)      !get scale sc
      call getBornKin2(klas,s,dble(qt)  ,  smin,ttt,tmax,iret) !get t and tmax
      if(iret.ne.0)return
      t=ttt

c parton cross-section between q2sft and q2cmin

      if(factsat.ge.0..and.(jpp.ge.5.or.jpp.eq.1).and.je1.ne.1
     .        .and.je2.ne.1.and.qt.le.max(q2cmin(1),q2cmin(2)))then
        is=-1
        factbgg=1./XsectionSat(x1,qt,1)/XsectionSat(x2,qt,2)
        if(qt.gt.qbmass**2)then
          factbqq=factbgg
        else
          factbqq=1.
        endif
        fact=factsat
        g1=  pifpartone(1,x1,sc, 0,0,1)
        uv1= 0d0
        dv1= 0d0
        sea1=pifpartone(1,x1,sc,-3,0,1)
        ch1=pifpartone(1,x1,sc,-4,0,1)
        bt1=pifpartone(1,x1,sc,-5,0,1)
        g2=  pifpartone(2,x2,sc, 0,0,1)
        uv2= 0d0
        dv2= 0d0
        sea2=pifpartone(2,x2,sc,-3,0,1)
        ch2=pifpartone(2,x2,sc,-4,0,1)
        bt2=pifpartone(2,x2,sc,-5,0,1)
      
      elseif(qt.gt.max(q2cmin(1),q2cmin(2)).and.jpp.lt.11)then

        fact=1.
        factbqq=1.
        factbgg=1.
        is=1
        
      if(jpp.ne.5)then
        ji1=mod(jpp+1,2)+1
        ji2=(jpp+1)/2
        sea1=pifpartone(1,x1,sc,-3,je1,ji1)
        ch1=pifpartone(1,x1,sc,-4,je1,ji1)
        bt1=pifpartone(1,x1,sc,-5,je1,ji1)
        g1=  pifpartone(1,x1,sc, 0,je1,ji1)
        uv1= pifpartone(1,x1,sc, 1,je1,ji1)
        dv1= pifpartone(1,x1,sc, 2,je1,ji1)
        sea2=pifpartone(2,x2,sc,-3,je2,ji2)
        ch2=pifpartone(2,x2,sc,-4,je2,ji2)
        bt2=pifpartone(2,x2,sc,-5,je2,ji2)
        g2=  pifpartone(2,x2,sc, 0,je2,ji2)
        uv2= pifpartone(2,x2,sc, 1,je2,ji2)
        dv2= pifpartone(2,x2,sc, 2,je2,ji2)
      else
        sea1=pifpartone(1,x1,sc,-3,je1,1)+pifpartone(1,x1,sc,-3,je1,2)
        ch1=pifpartone(1,x1,sc,-4,je1,1)+pifpartone(1,x1,sc,-4,je1,2)
        bt1=pifpartone(1,x1,sc,-5,je1,1)+pifpartone(1,x1,sc,-5,je1,2)
        g1=  pifpartone(1,x1,sc, 0,je1,1)+pifpartone(1,x1,sc, 0,je1,2)
        uv1= pifpartone(1,x1,sc, 1,je1,1)+pifpartone(1,x1,sc, 1,je1,2)
        dv1= pifpartone(1,x1,sc, 2,je1,1)+pifpartone(1,x1,sc, 2,je1,2)
        sea2=pifpartone(2,x2,sc,-3,je2,1)+pifpartone(2,x2,sc,-3,je2,2)
        ch2=pifpartone(2,x2,sc,-4,je2,1)+pifpartone(2,x2,sc,-4,je2,2)
        bt2=pifpartone(2,x2,sc,-5,je2,1)+pifpartone(2,x2,sc,-5,je2,2)
        g2=  pifpartone(2,x2,sc, 0,je2,1)+pifpartone(2,x2,sc, 0,je2,2)
        uv2= pifpartone(2,x2,sc, 1,je2,1)+pifpartone(2,x2,sc, 1,je2,2)
        dv2= pifpartone(2,x2,sc, 2,je2,1)+pifpartone(2,x2,sc, 2,je2,2)
      endif

      else
        
        return

      endif

      call ffvalues()

      ffsigj=ffsigj+
     *         ffborn(0,klas,tmax,s,t,gg,gq,qq,qa,qqp ,1,3, is*3 ) +
     *         ffborn(0,klas,tmax,s,t,0.,gc,0.,0. ,qc ,1,3, is*4 ) +
     *         ffborn(0,klas,tmax,s,t,0.,gb,0.,0. ,qb ,1,3, is*5 ) +
     *         ffborn(0,klas,tmax,s,t,0.,0.,cc,cac,0. ,1,3, is*6 ) +
     *         ffborn(0,klas,tmax,s,t,0.,0.,bb,bab,cb ,1,3, is*7 )

      ffsigj=ffsigj*dble(fact)

      end

c-----------------------------------------------------------------------
      double precision function ffsigNo(klas,qt,x1,x2,ihq)  !former psjy - No emissions
c-----------------------------------------------------------------------
#include "sem.h"
#include "aaa.h"
      double precision ffborn,pifpartone,x1,x2,t,s,  smin,ttt,tmax

      s=dble(engy*engy)*x1*x2
      ffsigNo=0d0

      call HardScale(klas,sc,qt,2)      !get scale sc
      call getBornKin2(klas,s,dble(qt)  ,  smin,ttt,tmax,iret) !get t and tmax
      if(iret.ne.0)return
      t=ttt

      if(abs(ihq).le.3)then
        n1=1
        n2=3
        if(abs(ihq).gt.0)then
          n1=abs(ihq)
          n2=abs(ihq)
        endif
      else
        if(abs(ihq).eq.4)then
          n1=11
          n2=12
        elseif(abs(ihq).eq.5)then
          n1=13
          n2=14
        else
          return
        endif
      endif

c parton cross-section between q2sft and q2cmin

      if(factsat.ge.0.and.qt.le.max(q2cmin(1),q2cmin(2)))then

        is=-1
        factbgg=1./XsectionSat(x1,qt,1)/XsectionSat(x2,qt,2)
        if(qt.gt.qbmass**2)then
          factbqq=factbgg
        else
          factbqq=1.
        endif
        fact=factsat

        g1=  pifpartone(1,x1,sc, 0,0,1)
        uv1= 0d0
        dv1= 0d0
        sea1=pifpartone(1,x1,sc, -3,0,1)
        ch1=pifpartone(1,x1,sc, -4,0,1)
        bt1=pifpartone(1,x1,sc, -5,0,1)
        g2=  pifpartone(2,x2,sc, 0,0,1)
        uv2= 0d0
        dv2= 0d0
        sea2=pifpartone(2,x2,sc, -3,0,1)
        ch2=pifpartone(2,x2,sc, -4,0,1)
        bt2=pifpartone(2,x2,sc, -5,0,1)
        
      elseif(qt.gt.max(q2cmin(1),q2cmin(2)))then

c parton cross-section above q2cmin

        is=1
        fact=1.
        factbqq=1.
        factbgg=1.
  
       g1=  pifpartone(1,x1,sc, 0,0,1)+pifpartone(1,x1,sc, 0,0,2)
       uv1= pifpartone(1,x1,sc, 1,0,1)+pifpartone(1,x1,sc, 1,0,2)
       dv1= pifpartone(1,x1,sc, 2,0,1)+pifpartone(1,x1,sc, 2,0,2)
       sea1=pifpartone(1,x1,sc,-3,0,1)+pifpartone(1,x1,sc,-3,0,2)
       ch1=pifpartone(1,x1,sc,-4,0,1)+pifpartone(1,x1,sc,-4,0,2)
       bt1=pifpartone(1,x1,sc,-5,0,1)+pifpartone(1,x1,sc,-5,0,2)
       g2=  pifpartone(2,x2,sc, 0,0,1)+pifpartone(2,x2,sc, 0,0,2)
       uv2= pifpartone(2,x2,sc, 1,0,1)+pifpartone(2,x2,sc, 1,0,2)
       dv2= pifpartone(2,x2,sc, 2,0,1)+pifpartone(2,x2,sc, 2,0,2)
       sea2=pifpartone(2,x2,sc,-3,0,1)+pifpartone(2,x2,sc,-3,0,2)
       ch2=pifpartone(2,x2,sc,-4,0,1)+pifpartone(2,x2,sc,-4,0,2)
       bt2=pifpartone(2,x2,sc,-5,0,1)+pifpartone(2,x2,sc,-5,0,2)

      else

        return
        
      endif

      call ffvalues()
      ffsigNo = ffsigNo +
     *         ffborn(0,klas,tmax,s,t,gg,gq,qq,qa,qqp ,n1,n2, is*3 ) +
     *         ffborn(0,klas,tmax,s,t,0.,gc,0.,0. ,qc ,n1,n2, is*4 ) +
     *         ffborn(0,klas,tmax,s,t,0.,gb,0.,0. ,qb ,n1,n2, is*5 ) +
     *         ffborn(0,klas,tmax,s,t,0.,0.,cc,cac,0. ,n1,n2, is*6 ) +
     *         ffborn(0,klas,tmax,s,t,0.,0.,bb,bab,cb ,n1,n2, is*7 )

      ffsigNo=ffsigNo*dble(fact)
      

      end

c-----------------------------------------------------------------------
      double precision function ffsigWi(klas,qt,x1,x2,ihq)    !former psjy - with emissions
c-----------------------------------------------------------------------
#include "sem.h"
#include "aaa.h"
      double precision ffborn,pifpartone,x1,x2,t,s,  smin,ttt,tmax

      s=dble(engy*engy)*x1*x2
      ffsigWi=0d0

      call HardScale(klas,sc,qt,2)      !get scale sc
      call getBornKin2(klas,s,dble(qt)  ,  smin,ttt,tmax,iret) !get t and tmax
      t=ttt

      if(iret.ne.0)return

      if(ihq.le.3)then
        n1=1
        n2=3
        if(ihq.gt.0)then
          n1=ihq
          n2=ihq
        endif
      else
        if(ihq.eq.4)then
          n1=11
          n2=12
        elseif(ihq.eq.5)then
          n1=13
          n2=14
        else
          return
        endif
      endif

      g1=  pifpartone(1,x1,sc, 0,1,1)+pifpartone(1,x1,sc, 0,1,2)
      uv1= pifpartone(1,x1,sc, 1,1,1)+pifpartone(1,x1,sc, 1,1,2)
      dv1= pifpartone(1,x1,sc, 2,1,1)+pifpartone(1,x1,sc, 2,1,2)
      sea1=pifpartone(1,x1,sc,-3,1,1)+pifpartone(1,x1,sc,-3,1,2)
      ch1=pifpartone(1,x1,sc,-4,1,1)+pifpartone(1,x1,sc,-4,1,2)
      bt1=pifpartone(1,x1,sc,-5,1,1)+pifpartone(1,x1,sc,-5,1,2)
      g2=  pifpartone(2,x2,sc, 0,1,1)+pifpartone(2,x2,sc, 0,1,2)
      uv2= pifpartone(2,x2,sc, 1,1,1)+pifpartone(2,x2,sc, 1,1,2)
      dv2= pifpartone(2,x2,sc, 2,1,1)+pifpartone(2,x2,sc, 2,1,2)
      sea2=pifpartone(2,x2,sc,-3,1,1)+pifpartone(2,x2,sc,-3,1,2)
      ch2=pifpartone(2,x2,sc,-4,1,1)+pifpartone(2,x2,sc,-4,1,2)
      bt2=pifpartone(2,x2,sc,-5,1,1)+pifpartone(2,x2,sc,-5,1,2)
      call ffvalues()
      ffsigWi= ffborn(0,klas,tmax,s,t,gg,gq,qq,qa,qqp ,n1,n2, 3 ) +
     *         ffborn(0,klas,tmax,s,t,0.,gc,0.,0. ,qc ,n1,n2, 4 ) +
     *         ffborn(0,klas,tmax,s,t,0.,gb,0.,0. ,qb ,n1,n2, 5 ) +
     *         ffborn(0,klas,tmax,s,t,0.,0.,cc,cac,0. ,n1,n2, 6 ) +
     *         ffborn(0,klas,tmax,s,t,0.,0.,bb,bab,cb ,n1,n2, 7 )
      end

c------------------------------------------------------------------------
      function ffsigiWi(qq,y0,ihq)                   !former psjx1  (sto)
c------------------------------------------------------------------------
c   integral over  ffsig....
c------------------------------------------------------------------------
c
c    dsigma/dpt_jet =  \int dy \int dx1  ffsig(fpartone,sis,x1,x2(x1))
c
c x1=xplus, x2=xminus
c x2=x2(x1) due to u+t+s=0
c ( s=x1*x2*spp, t/spp=-x1*xt*exp(-y)/2, u/spp=-x2*xt*exp(y)/2 )
c
c qq = pt**2,  xt=2.*sqrt(qq/s)
c rapidity range: 0 to y0
c
c ihq = 0 - all
c     = 1 - pp' -> pp' 
c     = 2 - pp~ -> p'p'~
c     = 3 - qq~ -> gg
c     = 4 - gg -> QQ~
c     = 5 - qq~ -> QQ~
c
c-----------------------------------------------------------------------
#include "sem.h"
#include "aaa.h"
      double precision z,sh,ffsigWi,ft
     &                ,xx1,xx2,xt,ymax,ymin,y,xmin,xmax
      ig=7
      ig1=7
      s=engy**2
      ffsigiWi=0.d0
      if(s.le.4.*qq)return
      if(qq.lt.max(q2cmin(1),q2cmin(2)))return
      ffsigi=0.
      xt=2d0*sqrt(dble(qq)/dble(s))
      ymax=min(dble(y0),log(1d0/xt+sqrt((1d0/xt-1d0)*(1d0/xt+1d0))))
      ymin=-ymax                          !final result must be divided by 2
      do i=1,ig
      do m=1,2
        y=.5d0*(ymax+ymin+(ymin-ymax)*dble((2*m-3)*tgss(ig,i)))
       !for xx1-integration, use variable x=xx1-xt*exp(y)/2.,with xmin<x<xmax
        xmin=xt**2/2.d0/(2.d0-xt*exp(-y))                    !condition x2<1
        xmax=1.d0-xt*exp(y)/2.d0                             !condition x1<1
        fx=0.
        do i1=1,ig1
        do m1=1,2
          xx1=xt*exp(y)/2.d0+xmin*(xmax/xmin)**dble(.5
     &                                           +tgss(ig1,i1)*(m1-1.5))
          xx2=xt*exp(-y)*xx1/(2.d0*xx1-xt*exp(y))
          z=xx1*xx2
          sh=z*s
          ft=0
          do klas=1,klasmax
            ft=ft+ffsigWi(klas,qq,xx1,xx2,ihq)
          enddo
          fx=fx+wgss(ig1,i1)*ft/sh**2*sh
        enddo
        enddo
        fx=fx*0.5*log(xmax/xmin)       !dx/x=0.5*log(xmax/xmin)dt (gauss)
        ffsigi=ffsigi+wgss(ig,i)*fx
      enddo
      enddo
      ffsigiWi=ffsigi*0.5*(ymax-ymin)    !dy=0.5*(ymax-ymin)dt (gauss)
     *  *(2*pi*pssalf(qq/qcdlam))**2  !alpha = 2*pi*pssalf
csp     /s**2
     *   *2*pi*sqrt(qq)                 !d2pt=2*pi*pt*dpt
     *   *2                            !2 jets are produced per hard process
     *   /2   ! y interval  2 * Delta_y
     *   /2   ! condition t < sqrt(s)/2,
              !     since t > sqrt(s)/2 is automatically included,
              !      see psbori
      return
      !kkkkkkkkkkkkkkkkkkkkk

      end

c------------------------------------------------------------------------
      function ffsigiNo(qq,y0,ihq)      !former psjx1  (sto)
c------------------------------------------------------------------------
c   integral over  ffsig....
c------------------------------------------------------------------------
c
c    dsigma/dpt_jet =  \int dy \int dx1  ffsig(fpartone,sis,x1,x2(x1))
c
c x1=xplus, x2=xminus
c x2=x2(x1) due to u+t+s=0
c ( s=x1*x2*spp, t/spp=-x1*xt*exp(-y)/2, u/spp=-x2*xt*exp(y)/2 )
c
c qq = pt**2,  xt=2.*sqrt(qq/s)
c rapidity range: 0 to y0
c
c ihq = 0 - all
c     = 1 - pp' -> pp' 
c     = 2 - pp~ -> p'p'~
c     = 3 - qq~ -> gg
c     = 4 - gg -> QQ~
c     = 5 - qq~ -> QQ~
c
c-----------------------------------------------------------------------
#include "sem.h"
#include "aaa.h"
      double precision z,sh,ffsigNo,ft
     &                ,xx1,xx2,xt,ymax,ymin,y,xmin,xmax

      ig=7
      ig1=7
      s=engy**2
      ffsigiNo=0.
      if(s.le.4.*qq)return
      if(qq.lt.max(q2cmin(1),q2cmin(2)))return
      ffsigi=0.
      xt=2d0*sqrt(dble(qq)/dble(s))
      ymax=min(dble(y0),log(1d0/xt+sqrt((1d0/xt-1d0)*(1d0/xt+1d0))))
      ymin=-ymax                          !final result must be divided by 2
      do i=1,ig
      do m=1,2
        y=.5d0*(ymax+ymin+(ymin-ymax)*dble((2*m-3)*tgss(ig,i)))
       !for xx1-integration, use variable x=xx1-xt*exp(y)/2.,with xmin<x<xmax
        xmin=xt**2/2.d0/(2.d0-xt*exp(-y))                    !condition x2<1
        xmax=1.d0-xt*exp(y)/2.d0                             !condition x1<1
        fx=0.
        do i1=1,ig1
        do m1=1,2
          xx1=xt*exp(y)/2.d0+xmin*(xmax/xmin)**dble(.5
     &                                           +tgss(ig1,i1)*(m1-1.5))
          xx2=xt*exp(-y)*xx1/(2.d0*xx1-xt*exp(y))
          z=xx1*xx2
          sh=z*s
          ft=0
          do klas=1,klasmax
            ft=ft+ffsigNo(klas,qq,xx1,xx2,ihq)
          enddo
          fx=fx+wgss(ig1,i1)*ft/sh**2*sh
        enddo
        enddo
        fx=fx*0.5*log(xmax/xmin)       !dx/x=0.5*log(xmax/xmin)dt (gauss)
        ffsigi=ffsigi+wgss(ig,i)*fx
      enddo
      enddo
      ffsigiNo=ffsigi*0.5*(ymax-ymin)    !dy=0.5*(ymax-ymin)dt (gauss)
     *  *(2*pi*pssalf(qq/qcdlam))**2  !alpha = 2*pi*pssalf
csp     /s**2
     *   *2*pi*sqrt(qq)                 !d2pt=2*pi*pt*dpt
     *   *2                            !2 jets are produced per hard process
     *   /2   ! y interval  2 * Delta_y
     *   /2   ! condition t < sqrt(s)/2,
              !     since t > sqrt(s)/2 is automatically included,
              !      see psbori
      return
      end




cc------------------------------------------------------------------------
c      function EsoftGluonInt(zz)  !removed 3.247
cc-----------------------------------------------------------------------

cc------------------------------------------------------------------------
c      function EsoftQuarkInt(zz)  !removed 3.247
cc-----------------------------------------------------------------------

c------------------------------------------------------------------------
c      function EsoftQZero(zz)      !removed 3.247
c-----------------------------------------------------------------------


c------------------------------------------------------------------------!bg
      subroutine gbxpxm(xp1,xp2,xmax1,xmin1,pth)
c------------------------------------------------------------------------
#include "aaa.h"
#include "sem.h"
#include "ems.h"
      double precision xp1,xp2,pth
      double precision xmax1,xmin1,dum

      dum=xp1
      dum=xp2
      dum=xmax1
      dum=xmin1
      dum=pth

      end

c########################################################################
c########################################################################
c########################################################################
      function om52pp(sy,xpp,xpm,b,iqq,je1,je2)   !modified om11pp
c########################################################################
c########################################################################
c########################################################################
c
c  om52pp    =  G_semi   / Fpart
c
c            =  G_semi_basic * factk / sigine * 10
c
c There is a factor .0389 because the elementary cross sections are in GeV-2.
c There a factor 1/.0389 in the normalization of the b dependence
c g(b)=exp(-b**2/(4.*.0389*R2))/(4*pi*.0389*R2).
c These two factors always cancel and are therefore not considered.
c
c  The factor  10 to get the cross section in fm**2 units
c  The factor factk is the K-factor for the elementary cross sections
c
c  G_semi_basic is the G without all these factors, formula see below
c
c-----------------------------------------------------------------------
c
c      sy  - energy squared for the hard interaction
c
c      b - impact parameter
c
c         in case of iqq=3: no b dependence (b is dummy)  !!!!!
c
c      iqq = 0  - sea-sea,
c      iqq = 1  - val-sea,    ( val(+) sea(-) )
c      iqq = 2  - sea-val,    ( sea(+) val(-) )
c      iqq = 3  - val-val,
c      iqq =10  - sat,
c
c      je = emission type
c               0 ... no emissions
c               1 ... emissions
c            else ... all
c
c----------------------------------------------------------------------
c
c  in the following: EsaturSea = EsaturQuark + EsaturGluon
c
c                    "\Sum" means sum over all combinations
c
c  for the different iqq:
c
c  0: G_semi_basic = \Sum \int dxp dxm EsaturSea(xp)*EsaturSea(xm)*sigmaHard(xp*xm*sy)
c  1: G_semi_basic = \Sum \int dxp dxm EsaturVal(xp)*EsaturSea(xm)*sigmaHard(xp*xm*sy)
c  2: G_semi_basic = \Sum \int dxp dxm EsaturSea(xp)*EsaturVal(xm)*sigmaHard(xp*xm*sy)
c  3: G_semi_basic = \Sum \int dxp dxm EsaturVal(xp)*EsaturVal(xm)*sigmaHard(xp*xm*sy)

c-----------------------------------------------------------------------
      common /ar3/    x1(7),a1(7)
      common /ar9/ x9(3),a9(3)
      common /psar7/  delx,alam3p,gam3p
#include "aaa.h"
#include "sem.h"
      double precision EsatSeaTil,EsatValTil
      double precision xp,xm,xx,psuds,xppii,uv1,dv1
      double precision ee44(0:3,0:3),sud01(0:1),sud02(0:1)
      dimension qq(2)
      logical q2int

      if(iqq.lt.0.or.(iqq.gt.3.and..not.iqq.eq.10))
     .stop'om52pp: unvalid  iqq'

      om52pp=0.
ccccccccccccccccccccccccccccc

      !  EsaturQuarkTil(xp,qq(1),1,999  <----- 999 is WRONG unless qq(1) small
      !
      !                 needs consideration of flavors, 
      !               look for  ee44(i,j)=  and  calcSTYP  and  calcWWXX

      return

ccccccccccccccccccccccccccccc

      ef1=0
      ef2=0
      ef3=0
      ef4=0
      if(iqq.lt.10)then
      if( je1.ge.1             .and. je2.ge.1)             ef1=1
      if( je1.ge.1             .and.(je2.eq.0.or.je2.eq.2))ef2=1
      if((je1.eq.0.or.je1.eq.2).and. je2.ge.1)             ef3=1
      if((je1.eq.0.or.je1.eq.2).and.(je2.eq.0.or.je2.eq.2))ef4=1
      else
      ef4=1  
      endif

      qqcut=max(q2cmin(1),q2cmin(2))
      if(iqq.lt.10)then
        ipomtype=iqq 
        spmin=4.*qqcut
        spmax=0.             !not used
        xmax=1.
        qq(1)=q2cmin(1)
        qq(2)=q2cmin(2)
      elseif(factsat.ge.0.)then
        ipomtype=-1
        qqcut=min(0.25*sy,qqcut)
        spmax=4.*qqcut
        q2mass=0. !qcmass**2
        qq(1)=max(q2sft,2.*q2mass)
        qq(2)=max(q2sft,2.*q2mass)
        spmin=4.*max(qq(1),qq(2))
c       spmax=spmax+4.*q2mass+0.5*_zos_inc_*qqcut
        spmax=min(sy,spmax)
c        qqcut=max(qqcut,2.*q2mass)
c        qqcut=min(qqcut,qbmass**2)
c        if(sy.lt.spmax)goto 999
        xmax=1. !spmax/sy
        if(max(qq(1),qq(2)).ge.qqcut)goto 999
      else
        return
      endif
      if(sy.le.spmin)goto 999

      !----------------------------------------------------------
      !the following is used used to compute (in the loops)
      !the factor g(b) for the b dependence, without 1/.0389:
      !g'(b)= 1/(4*pi*r2hhx) * exp(-b**2/(4.*.0389*r2hhx))
      !where r2hhx depends on zh or xm .
      !----------------------------------------------------------
      irh=0
      if(b.lt.0.)irh=1
      r2hh=r2had(iclpro)+r2had(icltar)
      r2hhs=r2hh+slopom*log(max(1.,sy))
c      r2hh=2.*r2part+(r2had(iclpro)-r2part)*(q2nmin/q2cmin(1))**1.
c     .              +(r2had(icltar)-r2part)*(q2nmin/q2cmin(2))**1.  !tp?????????
      zbb=exp(-b**2/(4.*.0389*r2hhs))
      fopi=4.*pi

      ipt=0
      icl=0
      iclv=2   !PDF fitted on proton, so deconvolution done for nucleons
      xppii=0d0
      if(iqq.eq.1)then
        xppii=xpp
        ipt=1  !KW2006 2        !soft on targ side
        icl=iclpro
      elseif(iqq.eq.2)then !kw????? missing in older versions
        xppii=xpm
        ipt=2 !KW2006 1        !soft on proj side
        icl=icltar
      endif

      xx=dble(sy)
      if(iqq.eq.1)then
        qqc=q2cmin(2)
      elseif(iqq.eq.2)then
        qqc=q2cmin(1)
      else
        qqc=max(q2cmin(1),q2cmin(2))
      endif
      
      delss=0.5*(dels(1)+dels(2)) !just to optimize integration
      if(iqq.eq.3)delss=-0.5    !better value for hard (no soft evolution)

      
c     double integration over qq(1) and qq(2) from q2sft to q2cmin
c     since psborn depends on max(qq(1),qq(2)) only, split the integration
c     in 2 pieces to avoid unnecessary psborn call
c     change variable qq to qi=1./qq for better precision

c do integration only for iqq>10
      if(iqq.ge.10)then
        iptqmx=2
        iqmx=3
        mqmx=2
        q2int=.true.
      else
        iptqmx=1
        iqmx=1
        mqmx=1
        q2int=.false.
      endif
      
      do iptq=1,iptqmx

        om52pp1=0.
        qi1max=1./q2sft
        qi1min=1./min(q2cmin(iptq),sy/4.)
        do iq1=1,iqmx
          do 100 mq1=1,mqmx
           if(q2int)then
            qi1=.5*(qi1max+qi1min + (2*mq1-3)*x9(iq1) * (qi1max-qi1min))
            qq(iptq)=1./qi1
            qq(3-iptq)=qq(iptq)       !temporary for SoftAtt

            spmin=4.*qq(iptq)
            if(spmin.gt.sy)goto 100
            
            xmin=spmin/sy
            !comments removed 3443o
           else
             qi1max=2.
             qi1min=0.
             q1i=1.
             xmin=spmin/sy
           endif

      !zzzzzzzzzzzzzzzzzzzzzzzzzzzzzzzzzzzzzzzzzzzzzzzzzzzzzzzzzzzz
      ximin=0.
      !zzzzzzzzzzzzzzzzzzzzzzzzzzzzzzzzzzzzzzzzzzzzzzzzzzzzzzzzzzzzz


      xmi1=xmin**(delh-delss)
c      alpq=0.5*(alpqua(1)+alpqua(2))
      xma1=xmax**(delh-delss)     

      !-----------------------------------------------------------------
      ! in all integrals t = (2*m-3)*x1(i)  (or similar)
      ! is the Gaussian basic integration variable, with limits -1 and 1
      !
      ! we compute \int ... dxp dxp (x_plus, x_minus)
      ! this is transformed into \int ... dzh dw,
      ! with zh=xp*xm, and where the second variable depends on iqq
      !-----------------------------------------------------------------

      om52pp2=0.

      !zh integration:
      do i=1,7
        do m=1,2
          zh=(.5*(xma1+xmi1-(2*m-3)*x1(i)*(xma1-xmi1)))
     *              **(1./(delh-delss))
          !-------------------------------------------
          !   zh-t Jacobian for this transformation:
          !  0.5*(xma1-xmi1)/(delh-delss)*zh*(1-delh+delss)
          !-------------------------------------------
          sss=zh*sy
          if(q2int)then
            qqcut=max(qq(1),qq(2))
            call setQvalues(qq(iptq),qq(3-iptq),qqcut) !needed for call calcSTYPa
          else
            call setQvalues(qq(1),qq(2),qqcut) 
          endif
          call calcSTYPa(iqq,ef1,ef2,ef3,ef4,sss) !computes sTYPz -> /zeroSTYP/
          !call calcSTYP(iqq,sss) 

          if(iqq.eq.0.or.iqq.eq.10)then

            !-----------------------------
            !the second variable is here uu
            !defined via xp=zh**uu
            !xp,xm-zh,uu Jacobian: -log(zh)
            !umin,umax: z**uu>ximin, z**(1-uu)>ximin
            !------------------------------
            stg=0.
            umax=1
            if(ximin.gt.zh)umax=log(ximin)/log(zh)
            umin=1-umax
            if(umin.lt.umax)then

              qi2max=1./q2sft
              qi2min=1./qq(iptq) !second integral only up to qq(iptq)
              do iq2=1,iqmx
                do mq2=1,mqmx

             if(q2int)then
                    
            qi2=.5*(qi2max+qi2min + (2*mq2-3)*x9(iq2) * (qi2max-qi2min))
            qq(3-iptq)=1./qi2
            
            sud01(0)=psuds(qq(1),0)   !from psborn
            sud01(1)=psuds(qq(1),1)   !from psborn
            sud02(0)=psuds(qq(2),0)   !from psborn
            sud02(1)=psuds(qq(2),1)   !from psborn

             else

            sud01(0)=1.d0
            sud01(1)=1.d0
            sud02(0)=1.d0
            sud02(1)=1.d0
            qi2max=2.
            qi2min=0.
            qi2=1.

             endif
            
             stg2=0.
             do i1=1,7
              do m1=1,2
                uu=.5*( umax+umin + (2*m1-3)*x1(i1) * (umax-umin) )
                !-----------------------------
                !   u-t Jacobian:  .5*(umax-umin)
                !-----------------------------
                xp=dble(zh**uu)
                xm=dble(zh)/xp
                do i4=0,3
                do j4=0,3
                 ee44(i4,j4)
     .            =EsatSeaTil(ipomtype,1,xp ,qq(1),i4)*sud01(min(i4,1))
     .            *EsatSeaTil(ipomtype,2,xm ,qq(2),j4)*sud02(min(j4,1))
                enddo
                enddo
                call calcWW4(0,ee44,wwgg,wwgq,wwqg,wwqq) ! Complete weights containing CS & Esat
                r2hhx=r2hh-slopom*log(zh)
                !if(iseahot.eq.0)then
                !  dstg= glu1*glu2*sgg*sud0gg
                !elseif(iseahot.eq.1)then
                !  dstg= sea1*glu2*sqg*sud0qg
                !elseif(iseahot.eq.2)then
                !  dstg= glu1*sea2*sgq*sud0gq
                !elseif(iseahot.eq.3)then
                !  dstg= sea1*sea2*sqq*sud0qq
                !else
                !  dstg= glu1*glu2*sgg*sud0gg
     *          !     +glu1*sea2*sgq*sud0gq+sea1*glu2*sqg*sud0qg
     *          !       +sea1*sea2*sqq*sud0qq
                !endif
                dstg=wwgg+wwgq+wwqg+wwqq 
                if(irh.eq.0)
     *          dstg=dstg * zbb**(r2hhs/r2hhx) / r2hhx  !g(b)
                stg2=stg2+a1(i1)*dstg
              enddo
             enddo
            
             if(q2int)stg2=a9(iq2)*stg2/qi2**2    !qi2 integration with Jacobian
             
             stg=stg+stg2*.5*(umax-umin)      !u-t Jacobian
            !factor zh**(-1-dels) from E=x**(-1-dels)*Etil cancels against zh-t Jacobian
              enddo                !qq(3-iptq) integration
             enddo
             stg=stg*.5*(qi2max-qi2min)
            endif
            om52pp2=om52pp2-log(zh)*a1(i)*stg/zh**(delh) !zh-t and xp/xm-zh/u Jacobian

          elseif(iqq.eq.3)then

            call calcWW4(3,ee44,sgg,sgq,sqg,sqq) ! here does not use ee44
            !-----------------------------
            !the second variable is here uu
            !defined via xp=1.d0-dble( ((1.d0-ww)/(1.d0+ww))**2 )
            !             with   ww=exp(-1.d0*uu)
            !xm-zh Jacobian:1/xp
            !xp-uu Jacobian sqrt(1-xp)*xp*w
            !------------------------------
            stq=0.
            xpmax=1.
            xpmin=zh
            umax=sqrt(1.-xpmin) !u=ln((1+sqrt(1-xp))/(1-sqrt(1-xp)))
            umax=log((1.+umax)/(1.-umax))
            !kw????? xpmin definition wrong in earlier versions ???????
            umin=0
            if(umax.gt.1.e-20)then
              do i1=1,7
                do m1=1,2
                  uu=.5*( umax+umin + (2*m1-3)*x1(i1) * (umax-umin) )
                  !-----------------------------
                  !   uu-t Jacobian:  .5*(umax-umin)
                  !-----------------------------
                  ww=exp(-1.d0*uu)
                  xp=1.d0-dble( ((1.d0-ww)/(1.d0+ww))**2 )
                  xm=dble(zh)/xp
                  if(xp*dble(xpp).le..9999d0.and.xm*dble(xpm).le..9999d0
     *              .or.xm*dble(xpp).le..9999d0
     *              .and.xp*dble(xpm).le..9999d0)then
                    sqqp=sqg !KW1808e special meaning, defined like this in calcSTYPa
                    sqaq=sgq !KW1808e special meaning, defined like this in calcSTYPa
                    stq=stq+a1(i1)
     *              *(psharf(xp,dble(xpp),xm,dble(xpm),sqq,sqqp,sqaq)
     *               +psharf(xm,dble(xpp),xp,dble(xpm),sqq,sqqp,sqaq))
     *               *sngl(sqrt(max(1d-12,1d0-xp)))   !xm-zh and xp-uu Jacobian
                  endif
                enddo
              enddo
              stq=stq*0.5*(umax-umin)                          !u-t Jacobian
     *               *0.5            !average psharf for 2 combinations of x+ and x-
            endif
            om52pp2=om52pp2+a1(i)*stq/zh**(delh-delss-1.)        !zh-t Jacobian

          elseif(iqq.eq.1.or.iqq.eq.2)then

            !-----------------------------
            !the second variable is here uu
            !defined via xp=cos(uu)**2
            !xm-zh Jacobian:1/xp
            !xp-u Jacobian:  2*sqrt(1d0-xp)*sqrt(xp)
            !------------------------------
            stq=0.
            xpmax=1
            if(ximin.gt.zh)xpmax=zh/ximin
            xpmin=zh
            umin=acos(sqrt(xpmax))
            umax=acos(sqrt(xpmin))
            if(umin.lt.umax)then
             do i1=1,7
              do m1=1,2
                uu=.5*( umax+umin + (2*m1-3)*x1(i1) * (umax-umin) )
                !-----------------------------
                !   uu-t Jacobian:  .5*(umax-umin)
                !-----------------------------
                xp=dble(cos(uu)**2)
                xm=dble(zh)/xp
                dstq=0
                if(xp*xppii.lt..99999d0.and.xp.ne.1.)then !bg
                  if(iqq.eq.1)then
                    val=EsatValTil(1,xp,xppii,qq(1),iclpro,0,uv1,dv1)
                    do j4=0,3
                      ee44(0,j4)=0
                      ee44(1,j4)=val*EsatSeaTil(ipomtype,2,xm,qq(2),j4)
                      ee44(2,j4)=0
                      ee44(3,j4)=0
                    enddo
                  elseif(iqq.eq.2)then
                    val=EsatValTil(2,xp,xppii,qq(2),icltar,0,uv1,dv1)
                    do i4=0,3
                      ee44(i4,0)=0
                      ee44(i4,1)=EsatSeaTil(ipomtype,1,xm,qq(1),i4)*val
                      ee44(i4,2)=0
                      ee44(i4,3)=0
                    enddo
                  endif
                  call calcWW4(iqq,ee44,wwgg,wwgq,wwqg,wwqq) ! Complete weights containing CS & Esat
                  r2hhx=r2hh-slopom*log(xm)
                  dstq=(wwgg+wwgq+wwqg+wwqq)
     *             *2*sngl(sqrt(max(1d-12,1d0-xp))*xp**(-0.5)) !xm-zh,xp-u Jacobian
                  if(irh.eq.0)
     *              dstq=dstq * zbb**(r2hhs/r2hhx) / r2hhx !g(b)
                  stq=stq+a1(i1)*dstq
                endif
              enddo
             enddo
            endif
            stq=stq*0.5*(umax-umin)                          !u-t Jacobian

            !zh**(-1-dels) from E=xpm**(-1-dels)*Etil cancels against zh-t Jacobian term
             om52pp2=om52pp2+a1(i)*stq/zh**(delh)              !zh-t Jacobian

          else
            stop'om52pp: unvalid  iqq (2).                 '
          endif

        enddo
      enddo

      if(q2int)om52pp2=a9(iq1)*om52pp2/qi1**2    !qi1 integration with Jacobian
      
      om52pp1=om52pp1+om52pp2*0.5*(xma1-xmi1)/(delh-delss) !zh-t Jacobian
c     *                       /sy**delh

 100  continue                          !qq(iptq) integration
      enddo

      om52pp=om52pp+om52pp1*0.5*(qi1max-qi1min)

      enddo                     !iptq loop
      
      if(iqq.eq.3.and.irh.eq.0)om52pp=om52pp*zbb**(r2hhs/r2hh)/r2hh !g(b)
      if(iqq.ge.10)om52pp=factsat*om52pp
      if(irh.eq.0)om52pp=om52pp/fopi                                !g(b)

      om52pp=om52pp * factk / sigine * 10     !*2. !/ 2.
 999  continue

      !test
      !gbasic = om52pp
      !. / ( factk / sigine * 10 / 2. )
      !print*,' '
      !print*,'22222 gbasic' ,iqq,sy/engy/engy, b,  gbasic

      end

c########################################################################
c########################################################################
c########################################################################
      function om11pp(sy,xpp,zb,iqq)   !former psfsh
c########################################################################
c########################################################################
c########################################################################
c
c      attention: G for val-val is NOT here => function om11pp3(...)
c
c-----------------------------------------------------------------------
c
c   om11pp =  G_semi / Fpart  / ( factk / sigine * 10 * sy**delh )
c
c          =  G_semi_basic / sy**delh
c
c expression for  G_semi_basic see function om52pp
c
c  A factor .0389  (because the elementary cross sections are in GeV-2)
c  is always cancelled aginst a factor 1/.0389 in the g(b) term
c  and therefore not considered.
c-----------------------------------------------------------------------
c om11pp - semihard interaction eikonal (om)
c sy  - energy squared of om
c xpp - momentum fraction of the valence quark
c z   - impact parameter factor, z=exp(-b**2/rp),
c iqq - type of the hard interaction:
c   0 - gg, 1 - qg, 2 - gq, 3 - gg(int), 4 - gg(proj), 5 - qg(proj),
c   6 - gg(int)|b=0, 7 - <b^2*gg(int)>, 8 - gg(proj)|b=0,
c   9 - <b^2*gg(proj)>, 10 - qg(proj)|b=0, 11 - <b^2*qg(proj)>
c
c Remark : spmin growth like x**esatur(=0.3) while sy growth like x.
c          As a consequence if at some point xmin, spmin>(s*xmin) then
c          om11p=0 for any x<xmin.
c-----------------------------------------------------------------------
      common /ar3/    x1(7),a1(7)
      common /psar7/  delx,alam3p,gam3p
#include "aaa.h"
#include "sem.h"
      double precision EsatSeaTil,EsatValTil
      double precision xp,xm,uv1,dv1
     .,xmin,xmax,xma1,xmi1,zh,xx,xpp
      double precision ee44(0:3,0:3)
      common/cpriom/npriom
      dimension qq(2)
      ipomtype=iqq

      z=zb
      irh=0
      if(zb.lt.0.)then
        z=1
        irh=1
      endif

      om11pp=0.
      if(iqq.eq.0.or.iqq.eq.3.or.iqq.eq.4
     *.or.iqq.eq.6.or.iqq.eq.7.or.iqq.eq.8.or.iqq.eq.9
     *.or.iclpro.lt.4.and.(iqq.eq.1.or.iqq.eq.5
     *.or.iqq.eq.10.or.iqq.eq.11)
     *.or.icltar.lt.4.and.iqq.eq.2)then
        spmin=4.*max(q2cmin(1),q2cmin(2)) !??????????????????????
      elseif(iclpro.eq.4 .or.icltar.eq.4) then !bg charm
        spmin=4.*max(q2cmin(1),q2cmin(2))+2.*qcmass**2
      else !bg bottom
        spmin=4.*max(q2cmin(1),q2cmin(2))+2.*qbmass**2
      endif
      if(sy.le.spmin)goto 999

      qq(1)=q2cmin(1)
      qq(2)=q2cmin(2)

      iclv=2   !PDF fitted on proton, so deconvolution done for nucleons
      if(iqq.eq.2)then
        ipt=2 !kw2006  1       !soft on proj side
        icl=icltar
      else
        ipt=1 !KW2006  2       !soft on targ side
        icl=iclpro
      endif
      delss=0.5*(dels(1)+dels(2))

      xmax=1.d0 
      xmin=dble(spmin/sy)
c limit energy ranged allowed for soft part according to Q2s
      xx=dble(sy)
      if(iqq.eq.1)then
        qqc=q2cmin(2)
      elseif(iqq.eq.2)then
        qqc=q2cmin(1)
      else
        qqc=max(q2cmin(1),q2cmin(2))
      endif
c      if(iqq.le.2)xmin=xmin*dble(SoftAtt(xx,qqc,2))
c      if(iqq.eq.0.or.iqq.eq.1)xmin=xmin*dble(SoftAtt(xx,q2cmin(2),2))
c      if(iqq.eq.0.or.iqq.eq.2)xmin=xmin*dble(SoftAtt(xx,q2cmin(1),1))
c      if(xmin.gt.xmax-dble(spmin/sy))return   !no calculation if not enough energy

      !zzzzzzzzzzzzzzzzzzzzzzzzzzzzzzzzzzzzzzzzzzzzzzzzzzzzzzzzzzzz
      ximin=0.
      !zzzzzzzzzzzzzzzzzzzzzzzzzzzzzzzzzzzzzzzzzzzzzzzzzzzzzzzzzzzzz

      xmi1=xmin**(delh-delss)
      xma1=xmax**(delh-delss)

      rp=r2had(iclpro)+r2had(icltar)+slopom*log(max(1.,sy))
c      r2hh=2.*r2part+(r2had(iclpro)-r2part)*(q2nmin/q2cmin(1))**1.
c     .              +(r2had(icltar)-r2part)*(q2nmin/q2cmin(2))**1.   !tp?????????
c      alpq=0.5*(alpqua(1)+alpqua(2))

c numerical integration over zh
      do i=1,7
      do m=1,2
        zh=(.5d0*(xma1+xmi1-dble(2*m-3)*dble(x1(i))*(xma1-xmi1)))
     *                                   **(1./(delh-delss))
        if(iqq.eq.0.or.iqq.eq.3.or.iqq.eq.4
     *  .or.iqq.eq.6.or.iqq.eq.7.or.iqq.eq.8.or.iqq.eq.9
     *  .or.iclpro.lt.4.and.(iqq.eq.1.or.iqq.eq.5
     *  .or.iqq.eq.10.or.iqq.eq.11)
     *  .or.icltar.lt.4.and.iqq.eq.2)then   
          sss=sngl(zh*sy)
          !~~~~~~~~~~~~~~~~~  
          call setQvalues(qq(1),qq(2),qqc) !needed for call calcSTYP
          call calcSTYP(iqq,sss)             !computes sTYPz -> /zeroSTYP/
          !~~~~~~~~~~~~~~~~~ 
        elseif(iclpro.eq.4 .or.icltar.eq.4) then ! KW1811??????? check charm proj case !!!
          call psjti0(sngl(zh*sy),sgq,sgqb,4,0)
          call psjti0(sngl(zh*sy),sqq,sqqb,4,1)
        else                                     !KW1811??????? check 
          call psjti0(sngl(zh*sy),sgq,sgqb,5,0)
          call psjti0(sngl(zh*sy),sqq,sqqb,5,1)
        endif

        if(iqq.eq.0.or.iqq.eq.3.or.iqq.eq.4
     *  .or.iqq.eq.6.or.iqq.eq.7.or.iqq.eq.8.or.iqq.eq.9)then

          stg=0.
          umax=1
          if(ximin.gt.zh)umax=log(ximin)/log(zh)
          umin=1-umax
          if(umin.lt.umax)then
          do i1=1,7
          do m1=1,2
            uu=.5*( umax+umin + (2*m1-3)*x1(i1) * (umax-umin) )
            xp=zh**uu
            xm=zh/xp
            do i4=0,3
            do j4=0,3
              ee44(i4,j4)=EsatSeaTil(ipomtype,1,xp ,qq(1),i4)
     .                   *EsatSeaTil(ipomtype,2,xm ,qq(2),j4)
            enddo
            enddo
            call calcWW4(iqq,ee44,wwgg,wwgq,wwqg,wwqq) ! Complete weights containing CS & Esat
            if(iqq.eq.0)then
              rh=r2had(iclpro)+r2had(icltar)-slopom*log(zh)
            elseif(iqq.eq.3.or.iqq.eq.4)then
              rh=1.
            elseif(iqq.eq.6.or.iqq.eq.7)then
              rh=alam3p-slopom*log(zh)
            elseif(iqq.eq.8.or.iqq.eq.9)then
              rh=r2had(iclpro)+.5*alam3p-slopom*log(zh)
            endif
            if(irh.eq.1)rh=1
            dstg=wwgg+wwgq+wwqg+wwqq 
            dstg=dstg*z**(rp/rh)/rh
            if(iqq.eq.7.or.iqq.eq.9)dstg=dstg*rh**2
            stg=stg+a1(i1)*dstg
                           !~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
                           if(nancheck(stg).eq.1
     .                     .or.npriom.eq.1)then
                            val=0-a1(i)*log(zh)*a1(i1)*dstg*(umax-umin)
     .                      /zh**(delh) *(xma1-xmi1)/(delh-delss)
     .                      /sy**delh/(4.*pi)/4.
                            print*,'om11pp:  ',val,xp,xm,xmin,xmax
                            if(nancheck(stg).eq.1)stop
                           endif
                           !~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
          enddo
          enddo
          endif
          stg=stg*(umax-umin)                 !u-t Jacobian          /2 later
          !zh**(-1-dels) from E=xpm**(-1-dels)*Etil cancels against zh-t Jacobian term
          om11pp=om11pp-a1(i)*log(zh)*stg/zh**(delh) !zh-t and xp/xm-zh/u Jacobian

        elseif(iqq.eq.1.or.iqq.eq.2)then

          stq=0.
          xpmax=1
          if(ximin.gt.zh)xpmax=zh/ximin
          xpmin=zh
          umin=acos(sqrt(xpmax))
          umax=acos(sqrt(xpmin))
          if(umin.lt.umax)then
          do i1=1,7
          do m1=1,2
            uu=.5*( umax+umin + (2*m1-3)*x1(i1) * (umax-umin) )
            xp=dble(cos(uu)**2)
            xm=zh/xp
            if(xp*xpp.lt..99999d0)then
              if(iqq.eq.1)then
                val=EsatValTil(1,xp,xpp,qq(1),iclpro,0,uv1,dv1)
                do j4=0,3
                  ee44(0,j4)=0
                  ee44(1,j4)=val*EsatSeaTil(ipomtype,2,xm,qq(2),j4)
                  ee44(2,j4)=0
                  ee44(3,j4)=0
                enddo
              elseif(iqq.eq.2)then
                val=EsatValTil(2,xp,xpp,qq(2),icltar,0,uv1,dv1)
                do i4=0,3
                  ee44(i4,0)=0
                  ee44(i4,1)=EsatSeaTil(ipomtype,1,xm ,qq(1),i4)*val
                  ee44(i4,2)=0
                  ee44(i4,3)=0
                enddo
              endif
              call calcWW4(iqq,ee44,wwgg,wwgq,wwqg,wwqq) ! Complete weights containing CS & Esat
              if(iqq.le.2)then
                rh=r2had(iclpro)+r2had(icltar)-slopom*sngl(log(xm))
              elseif(iqq.eq.5)then
                rh=1.
              elseif(iqq.le.10.or.iqq.le.11)then
                rh=r2had(iclpro)+.5*alam3p-slopom*sngl(log(xm))
              endif
              if(irh.eq.1)rh=1
              dstq=0
              if(xp.ne.1.d0)
     *        dstq=(wwgg+wwgq+wwqg+wwqq)
     *        *z**(rp/rh)/rh
     *        *sngl(xp**(-0.5d0)*sqrt(1d0-xp)) !xm-x and xp-u Jacobian   *2 later
              if(iqq.eq.11)dstq=dstq*rh**2
              stq=stq+a1(i1)*dstq
            endif
          enddo
          enddo
          endif
          stq=stq*(umax-umin)                      !u-t Jacobian       /2 later
          !zh**(-1-dels) from E=xpm**(-1-dels)*Etil cancels against zh-t Jacobian term
          om11pp=om11pp+a1(i)*stq/zh**(delh)
        else
          stop'ERROR 210907'
        endif
      enddo
      enddo
ctp      if(om11pp.gt.0..and.iqq.eq.0)
ctp     .write(ifch,*)'ici',sy,q2cmin,dels,betpom,xtest/om11pp

      om11pp=om11pp*(xma1-xmi1)/(delh-delss)    !zh-t Jacobian         /2 later
     *                       /sy**delh


      if(iqq.eq.0)then
        om11pp=om11pp/(4.*pi)/4. !from 1/(4*pi*lambda_NN) and 1/4 from Jacobian
      elseif(iqq.le.2)then
        om11pp=om11pp/(4.*pi)/2. !from 1/(4*pi*lambda_NN) and 1/2 from Jacobian
      elseif(iqq.eq.3)then
        om11pp=om11pp/4./pi/4.
     *        *4.*.0389
      elseif(iqq.eq.6)then
        om11pp=om11pp/4./pi/4.
      elseif(iqq.eq.7)then
        om11pp=om11pp/4./pi/4.
     *        *(4.*.0389)**2
      elseif(iqq.eq.4.or.iqq.eq.8.or.iqq.eq.9)then
        om11pp=om11pp/4./pi/4.
        if(iqq.eq.4)om11pp=om11pp*4.*.0389
        if(iqq.eq.9)om11pp=om11pp*(4.*.0389)**2
      elseif(iqq.eq.5.or.iqq.eq.10.or.iqq.eq.11)then
        om11pp=om11pp/4./pi/2.
        if(iqq.eq.5)om11pp=om11pp*4.*.0389
        if(iqq.eq.11)om11pp=om11pp*(4.*.0389)**2
      endif

      !For the test option "irh.eq.1", the 1/4./pi factor is not wanted:
      if(irh.eq.1)om11pp=om11pp * 4.*pi

 999  continue

      !test
      !gbasic = om11pp   * sy**delh
      !bzb=sign(1.,zb)*sqrt(-log(abs(zb))*4.*.0389*rp)
      !print*,'11111 gbasic',iqq,sy/engy/engy,bzb,gbasic
      !should be identical to gbasic in om52pp

      om11pp=max(1e-35,om11pp)

                       !~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
                       if(nancheck(om11pp).eq.1
     .                   .or.nancheck(log(om11pp)).eq.1)then
                         write(ifmt,*)'ERROR om11pp = ',om11pp,iqq
                         stop
                       endif
                       !~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

      end

c########################################################################
c########################################################################
c########################################################################
      function om11pp3(sy,xpp,xpm)
c########################################################################
c########################################################################
c########################################################################
c
c   om11pp3 = G_semi_qq /sy**delh
c
c        without g(b)  !!!!
c
c       and without .0389 factor because later the g(b) will be used
c       without 1/.0389 factor
c
c   Fpart does not matter here (=1 anyhow)
c
c-----------------------------------------------------------------------
c sy - energy squared for the pomeron,
c xpp - lc+ for the pomeron (valence quark),
c xpm - lc- for the pomeron (valence quark)
c-----------------------------------------------------------------------
      common /ar3/   x1(7),a1(7)
#include "aaa.h"
#include "sem.h"
      common /ckinirj/kinirj
      double precision xpp,xpm,xp,xm,zh,uu,umax,umin,xpmin,xpmax,ww
      double precision ee44(0:3,0:3)
      dimension qq(2)

      qq(1)=q2cmin(1)
      qq(2)=q2cmin(2)

      om11pp3=0.
      if(iclpro.lt.4.and.icltar.lt.4)then              !bg
        spmin=4.*max(q2cmin(1),q2cmin(2)) !??????????????????????
      elseif(iclpro.eq.4.or.icltar.eq.4) then          !bg charm
        spmin=4.*max(q2cmin(1),q2cmin(2))+2.*qcmass**2
      else                                             !bg bottom
        spmin=4.*max(q2cmin(1),q2cmin(2))+2.*qbmass**2
      endif
      if(sy.le.spmin)goto 999

      xmin=spmin/sy             !min hard pomeron mass share

      xmi1=xmin**(delh+.5)
      xma1=1.

      do i=1,7
      do m=1,2
        zh=dble(.5*(xma1+xmi1-(2*m-3)*x1(i)*(xma1-xmi1)))
     *     **(1./(delh+.5))
        if(iclpro.lt.4.and.icltar.lt.4)then              !bg
          call setQvalues(qq(1),qq(2),max(qq(1),qq(2))) 
          call calcSTYP(3,sngl(zh*sy))               !computes sTYPz -> /zeroSTYP/
          call calcWW4(3,ee44,sgg,sgq,sqg,sqq)   !ee44 here not used!!
          sqqp = sqg !special meaning (sqg is zero)
          sqaq = sgq !special meaning (sgq is zero)
        elseif(iclpro.eq.4.or.icltar.eq.4) then          !bg charm
          stop'22062020a update this'
          call psjti0(sngl(zh*sy),sqq,sqqb,4,1)
          sqq=0.
          sqaq=0.
        else                                             !bg bottom
          stop'22062020b update this'
          call psjti0(sngl(zh*sy),sqq,sqqb,5,1)
          sqq=0.
          sqaq=0.
        endif
        stq=0.  !int^1_(sqrt(z)) dx_p / x_p / sqrt(1-x_p) =int^(tmax)_(0) dt
        xpmax=1.d0
        xpmin=zh
        umax=sqrt(1d0-xpmin) !u=ln((1+sqrt(1-xp))/(1-sqrt(1-xp)))
        umax=log((1.d0+umax)/(1.d0-umax))
        !kw????? umax definition wrong in earlier versions ???????
        umin=0.d0
        if(umax.gt.1.d-20)then
        do i1=1,7
        do m1=1,2
          uu=.5*( umax+umin + (2*m1-3)*x1(i1) * (umax-umin) )
          ww=exp(-1.d0*uu)
          xp=1.d0-dble( ((1.d0-ww)/(1.d0+ww))**2 )
          xm=zh/xp
          if(xp*xpp.le..9999d0.and.xm*xpm.le..9999d0.or.
     *    xm*xpp.le..9999d0.and.xp*xpm.le..9999d0)then
          stq=stq+a1(i1)*(
     *    psharf(xp,xpp,xm,xpm,sqq,sqqp,sqaq)+
     *    psharf(xm,xpp,xp,xpm,sqq,sqqp,sqaq))
     *             *(1d0-xp)**0.5        !xm-x and xp-u Jacobian
          endif
        enddo
        enddo
        stq=stq*0.5*(umax-umin)                          !u-t Jacobian
     *          *0.5         !average psharf for 2 combinations of x+ and x-
        endif
        om11pp3=om11pp3+a1(i)*stq/zh**(delh-0.5)     !zh-t Jacobian
      enddo
      enddo
      om11pp3=om11pp3*(xma1-xmi1)/(delh+.5)/2.     !Jacobian factors
     *  /sy**delh

      !test
      !gbasic = om11pp3    * sy**delh
      !btest=0.5
      !if(kinirj.eq.0)
      !. print*,'11111 gbasic',3, sy/engy/engy,btest
      !.   ,gbasic,gbasic*gimpdep(btest)
      ! gbasic should be identical to gbasic in om52pp for negative b,z
      ! gbasic*gimpdep(btest) should be identical to gbasic in om52pp
      !          for positive b,z, and btest=b

 999  continue

      om11pp3=max(1e-35,om11pp3)

                         !~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
                         if(nancheck(om11pp3).eq.1
     .                     .or.nancheck(log(om11pp3)).eq.1)then
                           write(ifmt,*)'ERROR om11pp3 = ',om11pp3
                           stop
                         endif
                         !~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

      return
      end

      function gimpdep(b)
      ! g(b) without 1/.0389 factor
#include "aaa.h"
#include "sem.h"
      r2hh=r2had(iclpro)+r2had(icltar)
      gimpdep=exp(-b**2/(4.*.0389*r2hh)) / 4. /pi / r2hh
      end

c########################################################################
c########################################################################
c########################################################################
      function om11pp4(sy,zb,iqq)   !---MC--- 
c########################################################################
c########################################################################
c########################################################################
c wrong routine, but anyway irrelevant. Only the q2s dependence counts, 
c and we multiply anyway with a factor q2s*n
c----------------------------------------------------------------------
c sat pom
c sy - Pomeron mass     
c zb - impact param
c iqq - 10 only 
c----------------------------------------------------------------------
#include "aaa.h"
#include "sem.h"
#include "par.h"
      common /ar3/ x1(7),a1(7)
      common /ar4/ x4(2),a4(2)
      common /ar9/ x9(3),a9(3)
      double precision EsatSeaTil,xp,xm,dstg,stg
     .     ,xmin,xmax,xma1,xmi1,zh,xx,psuds,fxmin
     .     ,sud01(0:1),sud02(0:1)
      common/cpriom/npriom
      double precision ee44(0:3,0:3)
      dimension qq(2)

c      write(ifmt,*)'enter',sy,zb,iqq,q2cmin
      om11pp4=0.
c      return
      if(q2cmin(1).le.q2sft.or.q2cmin(2).le.q2sft)return
      if(iqq.ne.10)return
      ipomtype=0  !we take the wrong value (instread of -1) in order not to change this routine 

      z=zb
c      qqcut=min(0.25*sy,min(q2cmin(1),q2cmin(2)))
c take into account intrinsic pt
c      spmax=4.*(qqcut+min(_zos_inc_**2,qqcut))
c      spmax=4.*qqcut
      q2mass=0. !qcmass**2
c      if(iqq.eq.12)q2mass=qbmass**2
c      qq(1)=max(q2sft,sqrt(q2cmin(1)))
c      qq(2)=max(q2sft,sqrt(q2cmin(2)))
c      qq(1)=qq(1)+q2mass
c      qq(2)=qq(2)+q2mass
      qq(1)=max(q2sft,2.*q2mass)
      qq(2)=max(q2sft,2.*q2mass)
      spmin=4.*max(qq(1),qq(2))
c      spmax=spmax+4.*q2cmass**2+0.5*_zos_inc_*qqcut
c      qqcut=max(qqcut,2.*q2mass)
c      qqcut=min(qqcut,qbmass**2)
c      spmax=min(sy,spmax)
      if(sy.lt.spmin)return

c      write(ifmt,*)'pass'

      delss=0.5*(dels(1)+dels(2))
      rp=r2had(iclpro)+r2had(icltar)+slopom*log(max(1.,sy))
      xx=dble(sy)
      xmax=1d0 !dble(spmax/sy)   !can be limited by spmax, but not needed

      fxmin=1d0
c limit energy ranged allowed for soft part according to Q2s
c      fxmin=fxmin*dble(SoftAtt(xx,q2cmin(1),1))    
c      fxmin=fxmin*dble(SoftAtt(xx,q2cmin(2),2))
      
c     double integration over qq(1) and qq(2) from q2sft to q2cmin
c     since psborn depends on max(qq(1),qq(2)) only, split the integration
c     in 2 pieces to avoid unnecessary psborn call
c     change variable qq to qi=1./qq for better precision
      
      do ipt=1,2

        om11pp1=0d0
        qi1max=1./q2sft                   !squ
        qi1min=1./min(q2cmin(ipt),sy/4.)  !squ
c        qi1min=q2sft                       !lin
c        qi1max=min(q2cmin(ipt),sy/4.)      !lin
        do iq1=1,3 !2    !lin/squ
          do 100 mq1=1,2
c            qi1=.5*(qi1max+qi1min + (2*mq1-3)*x4(iq1) * (qi1max-qi1min)) !lin
            qi1=.5*(qi1max+qi1min + (2*mq1-3)*x9(iq1) * (qi1max-qi1min))  !squ
c            qq(ipt)=qi1         !lin
            qq(ipt)=1./qi1       !squ
            qq(3-ipt)=qq(ipt)       !temporary for SoftAtt

            spmin=4.*qq(ipt)

            if(spmin.gt.sy)goto 100
            
            xmin=dble(spmin)/dble(sy)
            !comments removed 3443o

      xmi1=xmin**(delh-delss)
      xma1=xmax**(delh-delss)

c      r2hh=2.*r2part+(r2had(iclpro)-r2part)*(q2nmin/q2cmin(1))**1.
c     .              +(r2had(icltar)-r2part)*(q2nmin/q2cmin(2))**1.   !tp?????????

c numerical integration over zh
      om11pp2=0d0
      do i=1,7
      do m=1,2
        zh=(.5d0*(xma1+xmi1-dble(2*m-3)*dble(x1(i))*(xma1-xmi1)))
     *                                   **(1./(delh-delss))

        sss=sngl(zh*sy)
        qqcut=max(qq(1),qq(2))
        call setQvalues(qq(ipt),qq(3-ipt),qqcut) 
        call calcSTYP(iqq,sss)   !KW2108   ,sgg,sgq,sqg,sqq)   !Compiler error!!!!!

        stg=0.
        umax=1
        umin=1-umax
        if(umin.lt.umax)then
          qi2max=1./q2sft           !squ
          qi2min=1./qq(ipt)        !second integral only up to qq(ipt)  !squ
c          qi2min=q2sft         !lin
c          qi2max=qq(ipt)        !second integral only up to qq(ipt)   !lin
          do iq2=1,3 !2    !lin/squ
          do mq2=1,2
c            qi2=.5*(qi2max+qi2min + (2*mq2-3)*x4(iq2) * (qi2max-qi2min)) !lin
            qi2=.5*(qi2max+qi2min + (2*mq2-3)*x9(iq2) * (qi2max-qi2min))  !squ
            qq(3-ipt)=1./qi2    !squ
c            qq(3-ipt)=qi2  !lin
            
            sud01(0)=psuds(qq(1),0)   !from psborn
            sud01(1)=psuds(qq(1),1)   !from psborn
            sud02(0)=psuds(qq(2),0)   !from psborn
            sud02(1)=psuds(qq(2),1)   !from psborn

          stg2=0.
          do i1=1,7
          do m1=1,2
            uu=.5*( umax+umin + (2*m1-3)*x1(i1) * (umax-umin) )
            xp=zh**uu
            xm=zh/xp
            do i4=0,3
              do j4=0,3
                ee44(i4,j4)
     .          =EsatSeaTil(ipomtype,1,xp ,qq(1),i4)*sud01(min(i4,1))
     .          *EsatSeaTil(ipomtype,2,xm ,qq(2),j4)*sud02(min(j4,1))
              enddo
            enddo
            call calcWW4(0,ee44,wwgg,wwgq,wwqg,wwqq) ! Complete weights containing CS & Esat

c            glu1=EsaturGluonTil(xp,qq(1),1)
c            sea1=EsaturQuarkTil(xp,qq(1),1,999) !per flavor
c            glu2=EsaturGluonTil(xm,qq(2),2)
c            sea2=EsaturQuarkTil(xm,qq(2),2,999) !per flavor
            rh=r2had(iclpro)+r2had(icltar)-slopom*log(zh)
c              rh=r2hh
c     *          -slopom*(log(xp)*(q2nmin/q2cmin(1))**1.
c     *                  +log(xm)*(q2nmin/q2cmin(2))**1.) !tp???????
c            dstg= glu1*glu2*dble(sgg)*sud0gg
c     *           +glu1*sea2*dble(sgq)*sud0gq+sea1*glu2*dble(sqg)*sud0qg
c     *           +sea1*sea2*dble(sqq)*sud0qq
            dstg=wwgg+wwgq+wwqg+wwqq 
            stg2=stg2+a1(i1)*dstg*z**(rp/rh)/rh

c            write(ifmt,*)'la',zh,sss,xp,xm,glu1,glu2,dstg
c                           !~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
c                           if(nancheck(stg).eq.1
c     .                     .or.npriom.eq.1)then
c                            val=0-a1(i)*log(zh)*a1(i1)*dstg*(umax-umin)
c     .                      /zh**(delh) *(xma1-xmi1)/(delh-delss)
c     .                        /sy**delh/(4.*pi)/4.
c                            write(ifmt,'(a,f7.4,2f7.2,2f9.2,5f7.2)')
c     .                         'om11pp4:  ',val,xp,xm,sgg,glu1,glu2
c                              if(nancheck(stg).eq.1)stop
c                           endif
c                           !~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
          enddo
          enddo
        stg=stg+a9(iq2)*stg2*(umax-umin)/qi2**2     !u-t and qi2 Jacobian !squ         /2 later
c        stg=stg+a4(iq2)*stg2*(umax-umin)     !u-t and qi2 Jacobian  /2 later  !lin
!zh**(-1-dels) from E=xpm**(-1-dels)*Etil cancels against zh-t Jacobian term
      enddo   !qq(3-ipt) integration
      enddo
      stg=stg*(qi2max-qi2min)

        else

      stg=0d0    
          
        endif
        om11pp2=om11pp2-a1(i)*log(zh)*stg/zh**(delh) !zh-t and xp/xm-zh/u Jacobian
c      write(ifmt,*)'ici',sy,q2cmin,qq,om11pp2,sgg,glu1,glu2,stg

      enddo
      enddo

      om11pp1=om11pp1+a9(iq1)*om11pp2/(delh-delss)*(xma1-xmi1)/qi1**2 !zh-t and qi1 Jacobian (without 1/2)   !squ
c      om11pp1=om11pp1+a4(iq1)*om11pp2/(delh-delss)*(xma1-xmi1) !zh-t and qi1 Jacobian (without 1/2)   !lin

 100  continue                          !qq(ipt) integration
      enddo

      om11pp4=om11pp4+om11pp1*(qi1max-qi1min)
      
      enddo   !ipt loop
      
      om11pp4=om11pp4/dble(4.*pi)/16d0  !from 1/(4*pi*lambda_NN) and 1/16 from 4 Jacobians (u, zh-t, qq1, qq2)
     *                       /sy**delh

      om11pp4=max(1e-35,om11pp4)
c     .     *(1d0-(dble(sy)*xx)**dble(spmin/spmax))**2
c      print *,'om11pp4',sy,om11pp4,spmin,spmax,xmin,xmax,q2cmin

      end

c------------------------------------------------------------------------
      function psharf(zp,zpp,zm,zpm,sqq,sqqp,sqaq)
c-----------------------------------------------------------------------
#include "sem.h"
#include "aaa.h"
      double precision EsaturValTil,zp,zpp,zm,zpm,fff,uv1,dv1,uv2,dv2

      if(zp.le..9999d0.and.zm.le..9999d0)then
        uv1=zp**dble(-1.-dels(1))
     *     *EsaturValTil(1,zp,zpp,q2cmin(1),1,iclpro)
        dv1=zp**dble(-1.-dels(1))
     *     *EsaturValTil(1,zp,zpp,q2cmin(1),2,iclpro)
        uv2=zm**dble(-1.-dels(2))
     *     *EsaturValTil(2,zm,zpm,q2cmin(2),1,icltar)
        dv2=zm**dble(-1.-dels(2))
     *     *EsaturValTil(2,zm,zpm,q2cmin(2),2,icltar)
        fff=0.d0
        if(iclpro.eq.2.and.icltar.eq.2)then       !proton
          fff=dble(sqq)*(uv1*uv2+dv1*dv2)+dble(sqqp)*(uv1*dv2+dv1*uv2)
        elseif(iclpro.eq.1.or.icltar.eq.1)then   !pion
          fff=dble(sqq)*uv1*uv2+dble(sqaq)*dv1*dv2
     *       +dble(sqqp)*(uv1*dv2+dv1*uv2)
        elseif(iclpro.eq.3.or.icltar.eq.3)then   !kaon
          fff=dble(sqq)*uv1*uv2+dble(sqqp)*(uv1*dv2+dv1*uv2+dv1*dv2)
        elseif(iclpro.eq.4.or.icltar.eq.4)then   !J/psi
          fff=dble(sqq)*uv1*(uv2+dv2)
        elseif(iclpro.eq.5.or.icltar.eq.5)then   !bg
          fff=dble(sqq)*uv1*(uv2+dv2)
        endif
        psharf=sngl(fff)
      else
        psharf=0.
      endif
      return
      end


c########################################################################
c########################################################################
c########################################################################
c
c                       ffom12   functions
c
c########################################################################
c########################################################################
c########################################################################

c----------------------------------------------------------------------
      function ffom12aii(iq1,iq2,je1,je2)   !---test---
c----------------------------------------------------------------------
#include "aaa.h"
      ig=5
      xmin=1./engy**2
      xmax=1
      r2=0
      do i2=1,ig
      do m2=1,2
c       xm=xmin+(xmax-xmin)*(.5+tgss(ig,i2)*(m2-1.5))
       xm=xmin*(xmax/xmin)**(.5+tgss(ig,i2)*(m2-1.5))
       r1=0
       do i1=1,ig
       do m1=1,2
c        xp=xmin+(xmax-xmin)*(.5+tgss(ig,i1)*(m1-1.5))
         xp=xmin*(xmax/xmin)**(.5+tgss(ig,i1)*(m1-1.5))
        f=ffom12a(xp,xm,iq1,iq2,je1,je2)
c        r1=r1+wgss(ig,i1)*f
        r1=r1+wgss(ig,i1)*f*xp
       enddo
       enddo
c       f=r1*0.5*(xmax-xmin)
       f=r1*0.5*log(xmax/xmin)*xm
       r2=r2+wgss(ig,i2)*f
      enddo
      enddo
      ffom12aii=r2*0.5*log(xmax/xmin)
c      ffom12aii=r2*0.5*(xmax-xmin)
      end

c----------------------------------------------------------------------
      function ffom12aj(x,iq1,iq2,je1,je2)   !---test---
c----------------------------------------------------------------------
#include "aaa.h"
      ig=5
      xmin=x
      xmax=1
      r2=0
      do i2=1,ig
      do m2=1,2
       xm=xmin*(xmax/xmin)**(.5+tgss(ig,i2)*(m2-1.5))
       xp=x/xm
       f=ffom12a(xp,xm,iq1,iq2,je1,je2)
       r2=r2+wgss(ig,i2)*f      !*xm /xm cancel
      enddo
      enddo
      ffom12aj=r2*0.5*log(xmax/xmin)
      end

c----------------------------------------------------------------------
      function ffom12ai(xp,xmini,iq1,iq2,je1,je2)   !---test---
c----------------------------------------------------------------------
#include "aaa.h"
      ig=5
      xmin=min(1.,max(xmini,1./(xp*engy**2)))
      xmax=1
      r2=0
      do i2=1,ig
      do m2=1,2
c       xm=xmin+(xmax-xmin)*(.5+tgss(ig,i2)*(m2-1.5))
       xm=xmin*(xmax/xmin)**(.5+tgss(ig,i2)*(m2-1.5))
       f=ffom12a(xp,xm,iq1,iq2,je1,je2)
c       r2=r2+wgss(ig,i2)*f
       r2=r2+wgss(ig,i2)*f*xm
       if(.not.(r2.le.0..or.r2.ge.0.))then
         print *,'ffom12ai',r2,q2cmin,xmini,f,xmin,xm,xmax
     .                     ,iq1,iq2,je1,je2
         write(ifch,*)'ffom12ai',r2,q2cmin,xmini,f,xmin,xm,xmax
     .                          ,iq1,iq2,je1,je2
         call utstop("\n\n ffom NaN catch \n\n&")
       endif
      enddo
      enddo
c      ffom12ai=r2*0.5*(xmax-xmin)
      ffom12ai=r2*0.5*log(xmax/xmin)
      end

c----------------------------------------------------------------------
      function ffom12a(xp,xm,iq1,iq2,je1,je2)   !---test---
c----------------------------------------------------------------------
#include "aaa.h"
      common/geom/rmproj,rmtarg,bmax,bkmx

      ig=5
      bmid=bkmx/2.
      r=0.d0
      do i=1,ig
        do m=1,2
          bb=bmid*(1.+(2.*m-3)*tgss(ig,i))
          f=ffom12(xp,xm,bb,iq1,iq2,je1,je2)
          r=r+bb*wgss(ig,i)*f
        enddo
      enddo
      ffom12a=r*2.*pi*bmid
      return
      end


c----------------------------------------------------------------------
      function ffom12(xp,xm,b,iq1,iq2,je1,je2)   !---test---
c----------------------------------------------------------------------
c
c      2*om52*F*F == PomInc
c
c  xp - xplus
c  xm - xminus
c
c  b - impact parameter
c                              iq=1 .... sea-sea
c  iq1 - min iq                iq=2 .... val-sea
c  iq2 - max iq                iq=3 .... sea-val
c                              iq=4 .... val-val
c  je = emission type (projectile and target side)
c          0 ... no emissions
c          1 ... emissions
c       else ... all
c
c----------------------------------------------------------------------
#include "aaa.h"
#include "sem.h"
      double precision Fsea,Fgfactor,Fval,Fvertex

      sy=engy*engy
      xh=xm*xp
      ffom12=0
      do i=iq1,iq2
       if(i.eq.1.or.i.eq.11)then       !sea-sea
        Fvertex=Fsea(dble(xp),1)
     *         *Fsea(dble(xm),2)
     *         *Fgfactor(dble(xp),dble(xm),1)
        ffom12=ffom12+om52pp(sy*xh,1.,1.,b,i-1,je1,je2)*sngl(Fvertex)
c        ffom12=ffom12+2.*om52pp(sy*xh,1.,1.,b,0,je1,je2)*sngl(Fvertex)
       elseif(i.eq.2)then   !val-sea
        Fvertex=Fval(dble(xp),1)
     *         *Fsea(dble(xm),2)
     *         *Fgfactor(1d0,dble(xm),1)
        ffom12=ffom12+om52pp(sy*xh,xp,1.,b,1,je1,je2)*sngl(Fvertex)
c        ffom12=ffom12+2.*om52pp(sy*xh,xp,1.,b,1,je1,je2)*sngl(Fvertex)
       elseif(i.eq.3)then   !sea-val
        Fvertex=Fsea(dble(xp),1)
     *         *Fval(dble(xm),2)
     *         *Fgfactor(dble(xp),1d0,1)
        ffom12=ffom12+om52pp(sy*xh,1.,xm,b,2,je1,je2)*sngl(Fvertex)
c        ffom12=ffom12+2.*om52pp(sy*xh,1.,xm,b,2,je1,je2)*sngl(Fvertex)
       elseif(i.eq.4)then   !val-val
        Fvertex=Fval(dble(xp),1)
     *         *Fval(dble(xm),2)
        ffom12=ffom12+om52pp(sy*xh,xp,xm,b,3,je1,je2)*sngl(Fvertex)
c        ffom12=ffom12+2.*om52pp(sy*xh,xp,xm,b,3,je1,je2)*sngl(Fvertex)
       endif
      enddo

      end


c########################################################################
c########################################################################
c########################################################################
c
c                       ffom11   functions
c
c########################################################################
c########################################################################
c########################################################################

c----------------------------------------------------------------------
      function ffom11aj(x,iq1,iq2)   !---test---
c----------------------------------------------------------------------
#include "aaa.h"
      ig=5
      xmin=x
      xmax=1
      r2=0
      do i2=1,ig
      do m2=1,2
       xm=xmin*(xmax/xmin)**(.5+tgss(ig,i2)*(m2-1.5))
       xp=x/xm
       f=ffom11a(xp,xm,iq1,iq2)
       r2=r2+wgss(ig,i2)*f      !*xm /xm cancel
      enddo
      enddo
      ffom11aj=r2*0.5*log(xmax/xmin)
      end

c----------------------------------------------------------------------
      function ffom11ai(xp,xmini,iq1,iq2)   !---test---
c----------------------------------------------------------------------
#include "aaa.h"
      ig=5
      xmin=min(1.,max(xmini,1./(xp*engy**2)))
      xmax=1
      r2=0
      do i2=1,ig
      do m2=1,2
c       xm=xmin+(xmax-xmin)*(.5+tgss(ig,i2)*(m2-1.5))
       xm=xmin*(xmax/xmin)**(.5+tgss(ig,i2)*(m2-1.5))
       f=ffom11a(xp,xm,iq1,iq2)
c       r2=r2+wgss(ig,i2)*f
       r2=r2+wgss(ig,i2)*f*xm
      enddo
      enddo
c      ffom11ai=r2*0.5*(xmax-xmin)
      ffom11ai=r2*0.5*log(xmax/xmin)
      end

c----------------------------------------------------------------------
      function ffom11a(xp,xm,iq1,iq2)   !---test---
c----------------------------------------------------------------------
c
c      int(db) om1ff /sigine*10
c
c  xp - xplus                  iq=-1 ... fit
c  xm - xminus                 iq=0 .... diff
c                              iq=1 .... gg
c  iq1 - min iq                iq=2 .... qg
c  iq2 - max iq                iq=3 .... gq
c                              iq=4 .... qq
c----------------------------------------------------------------------
#include "aaa.h"
      common/geom/rmproj,rmtarg,bmax,bkmx

      ig=5
      bmid=bkmx/2.
      r=0.d0
      do i=1,ig
        do m=1,2
          bb=bmid*(1.+(2.*m-3)*tgss(ig,i))
          f=ffom11(xp,xm,bb,iq1,iq2)
          r=r+bb*wgss(ig,i)*f
        enddo
      enddo
      ffom11a=r*2.*pi*bmid
      return
      end

c----------------------------------------------------------------------
      function ffom11(xp,xm,b,iq1,iq2)   !---test---
c----------------------------------------------------------------------
c
c       2*om5*F*F == PomInc
c
c         factor  1 / sigine * 10   included !!!!!!!!!!!!!!!
c
c  xp - xplus                  iq=-1 ... fit
c  xm - xminus                 iq=0 .... diff
c  b - impact parameter        iq=1 .... gg
c  iq1 - min iq                iq=2 .... qg
c  iq2 - max iq                iq=3 .... gq
c                              iq=4 .... qq
c----------------------------------------------------------------------
#include "aaa.h"
#include "sem.h"
      double precision om51

      ffom11=0.

      xmin=min(1.,1./(xp*engy**2))

      if(xm.ge.xmin)then

       xh=xm*xp
       yp=0.5*log(xp/xm)
c       ffom11=2.*sngl(om51(dble(xh),dble(yp),b,iq1,iq2))
       ffom11=sngl(om51(dble(xh),dble(yp),b,iq1,iq2))
     *     *(1-xm)**alplea(icltar)*(1-xp)**alplea(iclpro)

      elseif(xm.lt.0.)then   !xm integration

       ig=5
       xmax=1.
       r=0
       do i=1,ig
       do m=1,2
        xmm=xmin*(xmax/xmin)**(.5+tgss(ig,i)*(m-1.5))
        xh=xmm*xp
        yp=0.5*log(xp/xmm)
c        f=2.*sngl(om51(dble(xh),dble(yp),b,iq1,iq2))
        f=sngl(om51(dble(xh),dble(yp),b,iq1,iq2))
     *     *(1-xmm)**alplea(icltar)*(1-xp)**alplea(iclpro)
        r=r+wgss(ig,i)*f*xmm
       enddo
       enddo
       ffom11=r*0.5*log(xmax/xmin)

      endif

      ffom11= min(1.e10,ffom11 / sigine * 10)

      end



c########################################################################
c########################################################################
c########################################################################
c
c      parton-parton CS  psjet, psjet1, psborn
c
c      All functions have arguments
c
c        q1 = virtuality at current end of the ladder
c        q2 = virtuality at opposite end of the ladder
c        qqcut = sdditional virtuality cutoff 
c
c        Currently there is no need for an additiona variable (qqcut)
c        Howver, qqcut may be used connected to high pt triggered events
c
c-----------------------------------------------------------
c
c  Variable transforms, Jacobiens, limits   
c
c      the variable u is defined in [-1,1] to allow Gauss integration
c      the variable transform X->u should be such that f(X(u))*Jacob
c          is a smooth function in [-1,1]  (Jacob=dX/du)
c      This is not only useful for integration, but also for generation
c 
c  in all cases, for given q1, q2, s, j, l
c
c      j1=j                                  
c      l1=l                                 
c      call eolpic( j1 , l1 , j2 , l2 )
c      call eolpib( j2 , l2 , jl ) 
c      call getLimitsEoL(jl,s,qq , smin,tmn,tmx,iret)
c      qq=max(q1,q2) 
c
c  psjet
c 
c      dx1*dx2=dz*dx1/x1  z=x1*x2 
c      zmax=(1-epscutSud2(dble(q1)))*(1-epscutSud2(dble(q2)))      
c      zmin=max(smin/dble(s),1d0-zmax)   
c      zmx=zmax**(-delh)
c      zmn=zmin**(-delh)
c      z=((zmx+zmn-u*(zmn-zmx))/2)**(-1/delh)   !<------
c      Jacob=1/delh/2*(zmn-zmx)*z**(1+delh)        !<------
c      sh=z*dble(s)
c      t=2*tmn/(1+tmn/tmx-u*(1-tmn/tmx))  !<------
c      Jacob=t**2*0.5*(1/tmn-1/tmx)    !<------
c      xmax=0 
c      do klas=1,klasmax 
c       call getBornPtr(klas,sngl(sh),sngl(t) , pt2) !compute pt2
c       call HardScale(klas,scale,pt2,2)  !compute scale
c       xmax=max(xmax,1.d0-epscutSud2(dble(scale)))
c      enddo
c      xmin=max(dsqrt(z),z/xmax)  
c      if(xmax.gt..8d0)then
c        xmin1=max(xmin,.8d0)
c        x1=1-(1-xmax)*((1-xmin1)/(1-xmax))**((1-u)/2)  !<------
c        Jacob/x1=(1/x1-1)/2*log((1-xmin1)/(1-xmax))    !<------
c      endif 
c      if(xmin.lt..8d0)then
c        xmax1=min(xmax,.8d0)
c        x1=xmin*(xmax1/xmin)**dble(.5*(1+u))  !<------
c        Jacob/x1=log(xmax1/xmin)/2       !<------
c      endif
c      ---for rejection---
c      f=0      
c      do klas=1,klasmax 
c        call HardScale(klas,qq,pt2min,1) !get pt2min
c        pt2min=max(0.,pt2min)
c        call getBornKin(klas,s,pt2min, smn, tmn, tmx,iret) !get tmx
c        call getBornPtr(klas,s,t , pt2) !compute pt2
c        call HardScale(klas,scale,pt2,2)  !compute scale
c        f=f+psjeti(0,klas,tmx,q1,q2,scale,t,x1,x2,sh,j,l,jdis)
c      enddo
c
c  psjet1
c
c      xmax=1.d0-epscutSud2(dble(q1))       
c      xmin=max(smin/dble(s),1d0-xmax)    
c      if(xmax.gt..8d0)then
c        zmn=max(xmin,.8d0)
c        zmx=xmax
c        z=1-(1-zmx)*((1-zmn)/(1-zmx))**((1-u)/2) !<-----
c        Jacob=(1-z)/2*log((1-zmn)/(1-zmx))    !<------
c        sh=z*dble(s)
c        t=2*tmn/(1+tmn/tmx-u*(1-tmn/tmx))  !<------
c        Jacob=t**2*0.5*(1/tmn-1/tmx)    !<------
c      endif
c      if(xmin.lt..8d0)then
c        zmx=min(xmax,.8d0)**(-delh)
c        zmn=xmin**(-delh)
c        z=((zmx+zmn-u*(zmn-zmx))/2)**(-1/delh)   !<------
c        Jacob=1/delh/2*(zmn-zmx)*z**(1+delh)        !<------
c        sh=z*dble(s)
c        t=2*tmn/(1+tmn/tmx-u*(1-tmn/tmx))  !<------
c        Jacob=t**2*0.5*(1/tmn-1/tmx)    !<------
c      endif
c      ---for rejection---
c      f=0      
c      do klas=1,klasmax
c        call HardScale(klas,qq,pt2min,1) !get pt2min
c        pt2min=max(0.,pt2min)
c        call getBornKin(klas,s,pt2min, smn, tmn, tmx,iret) !get tmx
c        call getBornPtr(klas,s,t , pt2) !compute pt2
c        call HardScale(klas,scale,pt2,2)  !compute scale
c        f=f+psjetj(0,klas,tmx,q1,q2,scale,t,z,sh,j,l)*psuds(scale,l)/psuds(q2,l)
c      enddo 
c
c  psborn
c
c      t=2*tmn/(1+tmn/tmx-u*(1-tmn/tmx))  !<------
c      Jacob=t**2*0.5*(1/tmn-1/tmx)    !<------
c      ---for rejection---
c      f=0      
c      do klas=1,klasmax
c        call HardScale(klas,qq,pt2min,1) !get pt2min
c        pt2min=max(0.,pt2min)
c        call getBornKin(klas,s,dble(pt2min), smn, tmn, tmx,iret) !get tmx
c        call getBornPtr(klas,s,t , pt2) !compute pt2
c        call HardScale(klas,scale,pt2,2)  !compute scale
c        f=f+psjetk(0,klas,tmx,q1,q2,scale,t,s,j,l,jdis)
c        f=f*psuds(scale,j1)*psuds(scale,l1)
c        f=f/psuds(q1,j1)/psuds(q2,l1) 
c      enddo
c
c  Other factors:
c
c    A = pi * (2*pi*pssalf(scale/qcdlam))**2 / s**2 
c
c         because of dsigma/dt = A * psbori
c
c########################################################################
c########################################################################
c########################################################################


c-----------------------------------------------------------------------
      function psjet(q1,q2,qqcut,s,j,l,jdis)
c-----------------------------------------------------------------------
c                     Integral over psjeti  
c-----------------------------------------------------------------------
c            both sides parton ladder cross-section
c              (at least one emission on each side)
c-----------------------------------------------------------------------
c q1 - virtuality cutoff at current end of the ladder;
c q2 - virtuality cutoff at opposite end of the ladder;
c qqcut - p_t cutoff for the born process;
c s - c.m. energy squared for the scattering;
c j - parton type at current end of the ladder (0 - g, 1,2 etc. - q);
c l - parton type at opposite end of the ladder (0 - g, 1,2 etc. - q).
c-----------------------------------------------------------------------
      double precision xx1,xx2,qq,xmin,xmax,xmin1,xmax1,t
     *,sh,z,ft,fx1,fx2,ss, smn, tmn, tmx,epscutSud2
      common /ar3/   x1(7),a1(7)
#include "aaa.h"
#include "sem.h"
      common/ccctest/iiitest
      integer ienvi
      common /cienvi/ ienvi
      ienvi=102
      iiitest=0

      psjet=0.

      if(iabs(j).eq.3.and.  noflav(q1).lt.4 ) return ! 1 charm
      if(iabs(j).eq.4.and.  noflav(q1).lt.5 ) return ! 1 bottom
      if(iabs(l).eq.3.and.  noflav(q2).lt.4 ) return ! 2 charm
      if(iabs(l).eq.4.and.  noflav(q2).lt.5 ) return ! 2 bottom

      if(jdis.eq.0)then
        ik=1
      else
        ik=4
      endif
      qq=dble(max(q1/ik,q2))
      qq=max(qq,dble(qqcut))
      qmin=qq
      ss=dble(s)

      j1=j                                  
      l1=l                                 
      call eolpic( j1 , l1 , j2 , l2 )
      call eolpib( j2 , l2 , jl ) 
      call getLimitsEoL(jl,s,qmin , smin,tmin,tmax,iret)
      if(iret.ne.0)return

      call getQmaxEoL(jl,s,scamax)   !KW1908

      !sabmin=smin   ! min s after branching
      !compute ssx = minimum s before branching
      !ssx=sabmin/(1.d0-epscutSud2(dble(scamax)))
      !.  /(1.d0-epscutSud2(dble(scamax))) 
      !if(s.lt.ssx)return
      !phi=acos(1.-54.*q2ini/s)/3.
      !zmax=(1.+2.*cos(phi))**2/9.    
      !zmin=(1.-cos(phi)+sqrt(3.d0)*sin(phi))/3.
      !zmin=max(zmin**2,sngl(sabmin/dble(s)))
      !zmax=(1-epscutSud2(dble(scamax)))*(1-epscutSud2(dble(scamax)))       !KW1812
      !zmin=max(sabmin/dble(s),1d0-zmax)    !KW1812

      zmin=dble(smin)/dble(s)                                               !KW2004
      zmax=(1-epscutSud2(dble(scamax)))*(1-epscutSud2(dble(scamax)))       !KW1812
      if(zmin.ge.zmax)return               !KW1812

      zmin=zmin**(-delh)
      zmax=zmax**(-delh)
      do i=1,7
      do m=1,2
       z=dble(.5*(zmax+zmin+(zmin-zmax)*(2*m-3)*x1(i)))**(-1./delh)
       xmin=dsqrt(z)
       sh=z*dble(s)
       qqx=qq !max(qq,dble(q2ini)/(1.d0-dsqrt(z))) 
       do klas=1,klasmax
        call HardScale(klas,qqx,pt2min,1)  !get pt2min 
        ft=0.d0
        !NEEDED FOR abs(coelaf-1).lt.1e-5.and.abs(coekaf-1).lt.1e-5 :
        if(pt2min.gt.0)then 
        call getBornKin2(klas,sh,dble(pt2min), smn, tmn, tmx,iret) !get tmn,tmx
        if(iret.eq.0)then
         do i1=1,7
         do m1=1,2
          t=2.d0*tmn/(1.d0+tmn/tmx-dble(x1(i1)*(2*m1-3))
     &    *(1.d0-tmn/tmx))
          call getBornPtr(klas,sngl(sh),sngl(t) , pt2) !compute pt2
          call HardScale(klas,scale,pt2,2)  !compute scale
          xmax=1.d0-epscutSud2(dble(scale))
          xmin=max(dsqrt(z),z/xmax)   !integration over lower triangle of x1-x2 plane (x1<x2)
          fx1=0.d0                    !therfore we sum psjeti(...xx1,xx2...) and psjeti(...xx2,xx1...)
          fx2=0.d0
          if(xmin.lt.xmax)then
          if(xmax.gt..8d0)then
            xmin1=max(xmin,.8d0)
            do i2=1,7
            do m2=1,2
              xx1=1.d0-(1.d0-xmax)*((1.d0-xmin1)/(1.d0-xmax))**
     *        dble(.5+x1(i2)*(m2-1.5))
              xx2=z/xx1
              a=psjeti(0,klas,tmx,q1,q2,scale,t,xx1,xx2,sh,j,l,jdis) !compensate triangular 
              b=psjeti(0,klas,tmx,q1,q2,scale,t,xx2,xx1,sh,j,l,jdis) !integration domaine
              fb=(a+b)*pssalf(scale/qcdlam)**2
              !if(fb.lt.0.)then !negative cross section
              !  write(ifmt,*) 'neg CS :   fb = ',fb,a,b
              !  stop'ERROR 18122018'
              !endif
              fx1=fx1+dble(a1(i2)*fb)*(1.d0/xx1-1.d0)          
            enddo
            enddo
            fx1=fx1*dlog((1.d0-xmin1)/(1.d0-xmax))
          endif
          if(xmin.lt..8d0)then
            xmax1=min(xmax,.8d0)
            do i2=1,7
            do m2=1,2
              xx1=xmin*(xmax1/xmin)**dble(.5+x1(i2)*(m2-1.5))
              xx2=z/xx1
              a=psjeti(0,klas,tmx,q1,q2,scale,t,xx1,xx2,sh,j,l,jdis) !compensate triangular
              b=psjeti(0,klas,tmx,q1,q2,scale,t,xx2,xx1,sh,j,l,jdis) !integration domaine
              fb=(a+b)*pssalf(scale/qcdlam)**2
              fx2=fx2+dble(a1(i2)*fb)
            enddo
            enddo
            fx2=fx2*dlog(xmax1/xmin)
          endif
          endif
          ft=ft+dble(a1(i1))*(fx1+fx2)*t**2
         enddo
         enddo
         psjet=psjet + a1(i)*sngl(ft*z**(1.+delh)/sh**2)
     &                    *(1.d0/tmn-1.d0/tmx)!part of t-Jacobian, MUST be inside klas loop
        endif !iret
        endif !pt2min.gt.0
       enddo !klas
      enddo
      enddo
      tt1=psjet
      psjet=psjet*(zmin-zmax)/delh*4*pi**3 
     *         /8.    !t-z-x-Jacobians

      tt=psjet
      if(.not.(tt.le.0..or.tt.ge.0.) !NaN catch
     .   .or.tt.gt.1e35 .or. tt.lt.-1e35 )then !Infinity catch
         write(ifmt,*) 'NaN / Inf catch :   psjet = ',tt
         write(ifmt,*) 'tt1:   ',tt1
         stop'ERROR 26112018'
      endif

      return
      end

c------------------------------------------------------------------------
      function psjet1(q1,q2,qqcut,s,j,l,jdis)
c-----------------------------------------------------------------------
c                    Integral over psjetj
c-----------------------------------------------------------------------
c              ordered parton ladder cross-section
c  (at least one emission on current side, no emission on opposite side)
c-----------------------------------------------------------------------
c q1 - virtuality cutoff at current end of the ladder;
c q2 - virtuality cutoff at opposite end of the ladder;
c qqcut - p_t cutoff for the born process;
c s - c.m. energy squared for the scattering;
c j - parton type at current end of the ladder (0 - g, 1,2 etc. - q);
c l - parton type at opposite end of the ladder (0 - g, 1,2 etc. - q).
c-----------------------------------------------------------------------
      double precision xx,z,qq,xmax,xmin,xmin1
     *,sh,t,xmax1,fx1,fx2,psuds, smn, tmn, tmx,epscutSud2
      common /ar3/   x1(7),a1(7)
#include "aaa.h"
#include "sem.h"
      integer ienvi
      common /cienvi/ ienvi
      ienvi=103

      psjet1=0.

      if(iabs(j).eq.3.and.  noflav(q1).lt.4 ) return ! 1 charm
      if(iabs(j).eq.4.and.  noflav(q1).lt.5 ) return ! 1 bottom
      if(iabs(l).eq.3.and.  noflav(q2).lt.4 ) return ! 2 charm
      if(iabs(l).eq.4.and.  noflav(q2).lt.5 ) return ! 2 bottom

      if(jdis.eq.0.or.jdis.eq.2)then
        ik=1
      else
        ik=4
      endif
      if(jdis.eq.0)then
        qq=dble(max(q1,q2))
      elseif(jdis.eq.1)then
        qq=dble(max(q1/4.,q2))
      else
        qq=dble(max(q1,q2/4.))
      endif
      qq=max(qq,dble(qqcut))
      qmin=qq

      j1=j                                  
      l1=l                                 
      call eolpic( j1 , l1 , j2 , l2 )
      call eolpib( j2 , l2 , jl ) 
      call getLimitsEoL(jl,s,qmin , smin,tmin,tmax,iret)
      if(iret.ne.0)return

      call getQmaxEoL(jl,s,scamax)   !KW1908

      !sabmin=dble(smin)   ! min s after branching
      !compute smin = minimum s before branching
      !smin=sngl(sabmin/(1.d0-epscutSud2(dble(scamax)/ik)))  !KW1812  
      !if(s.le.smin)return
      !call getM2( j2 , l2 , q2mass ) 
      !q2ms=dble(q2mass)/dble(s)
      !q2inis=dble(q2ini)/dble(s)
      !xmax=.5d0*(1.d0+q2ms)+dsqrt(.25d0*(1.d0-q2ms)**2-4.d0*q2inis)
      !xmin=max(1.d0+q2ms-xmax,s2min/dble(s))
      !xmax=1.d0-epscutSud2(dble(scamax))       !KW1812
      !xmin=max(sabmin/dble(s),1d0-xmax)    !KW1812

      xmin=dble(smin)/dble(s)            !KW2004
      xmax=1.d0-epscutSud2(dble(scamax))       !KW1812
      if(xmin.ge.xmax)return               !KW1812

      fx1=0.d0
      fx2=0.d0
      if(xmax.gt..8d0)then
        xmin1=max(xmin,.8d0)
        do i=1,7
        do m=1,2
         z=1.d0-(1.d0-xmax)*((1.d0-xmin1)/(1.d0-xmax))**
     *   (.5d0+dble(x1(i)*(m-1.5)))
         sh=z*dble(s)
         xx=z
         qqx=qq !max(qq,dble(q2ini)/(1.d0-dsqrt(z))) 
         do klas=1,klasmax !here we must sum over all klas, since jl is not in-Born pair
          call HardScale(klas,qqx,pt2min,1)  !get pt2min 
          ft=0.
          !NEEDED FOR abs(coelaf-1).lt.1e-5.and.abs(coekaf-1).lt.1e-5 :
          if(pt2min.gt.0)then 
          call getBornKin2(klas,sh,dble(pt2min), smn, tmn, tmx,iret)
          if(iret.eq.0)then
           do i1=1,7
           do m1=1,2
             t=2.d0*tmn/(1.d0+tmn/tmx-dble(x1(i1)*(2*m1-3))
     &       *(1.d0-tmn/tmx))
             call getBornPtr(klas,sngl(sh),sngl(t) , pt2) !get pt2
             call HardScale(klas, scale , pt2 ,2)  !get scale
c      if(.not.(scale.le.0..or.scale.ge.0.))then !NaN catch
c            print*,scale
c            stop'NaN in loop 1 in psjet1'
c      endif
             if(jdis.eq.0)then
               scale1=scale
               scale2=scale
             elseif(jdis.eq.1)then
               scale1=scale*4.
               scale2=scale 
             elseif(jdis.eq.2)then
               scale1=scale 
               scale2=scale*4. 
             endif
             fb=
     .       psjetj(0,klas,tmx,q1,q2,scale1,t,xx,sh,j,l)
     .       *pssalf(scale/qcdlam)**2*psuds(scale2,l)
             ft=ft+a1(i1)*fb*sngl(t**2)
             if(.not.(ft.le.0..or.ft.ge.0.))then !NaN catch
             print*,ft, psjetj(0,klas,tmx,q1,q2,scale1,t,xx,sh,j,l)
             stop'NaN in psjet1 ft' 
             endif
           enddo
           enddo
           fx1=fx1 + dble(a1(i)*ft)/sh**2*(1.d0-z)
     &                *(1.d0/tmn-1.d0/tmx) !part of t-x Jacobien, MUST be inside klas loop
           !~~~~~~~~~~~~~~~~~~~~~~
           if(.not.(fx1.le.0..or.fx1.ge.0.))then !NaN catch
           print*,'klas,fx1,sh**2,tmn,tmx,qqx,pt2min:'
     .     ,klas,fx1,sh**2,tmn,tmx,qqx,pt2min
           stop'NaN in psjet1 fx1' 
           endif 
           !~~~~~~~~~~~~~~~~~~~~~~ 
          endif
          endif
         enddo !klas
        enddo
        enddo
        fx1=fx1*dlog((1.d0-xmin1)/(1.d0-xmax))
      endif
      if(xmin.lt..8d0)then
        xmax1=min(xmax,.8d0)**(-delh)
        xmin1=xmin**(-delh)
        do i=1,7
        do m=1,2
         z=(.5d0*(xmax1+xmin1+(xmin1-xmax1)*dble((2*m-3)*x1(i))))
     *   **(-1./delh)
         sh=z*dble(s)
         xx=z
         qqx=qq !max(qq,dble(q2ini)/(1.d0-z)/ik) 
         do klas=1,klasmax
          call HardScale(klas,qqx,pt2min,1)  !get pt2min 
          ft=0.d0
          !NEEDED FOR abs(coelaf-1).lt.1e-5.and.abs(coekaf-1).lt.1e-5 :
          if(pt2min.gt.0)then 
          call getBornKin2(klas,sh,dble(pt2min), smn, tmn, tmx,iret)
          if(iret.eq.0)then
           do i1=1,7
           do m1=1,2
             t=2.d0*tmn/(1.d0+tmn/tmx-dble(x1(i1)*(2*m1-3))
     &       *(1.d0-tmn/tmx))
             call getBornPtr(klas,sngl(sh),sngl(t) , pt2) !get pt2
             call HardScale(klas, scale , pt2 ,2)  !get scale
             if(jdis.eq.0)then
               scale1=scale
               scale2=scale
             elseif(jdis.eq.1)then
               scale1=scale*4.
               scale2=scale
             elseif(jdis.eq.2)then
               scale1=scale
               scale2=scale*4.
             endif
             fb=
     .       psjetj(0,klas,tmx,q1,q2,scale1,t,xx,sh,j,l)
     .       *pssalf(scale/qcdlam)**2*psuds(scale2,l)
             ft=ft+a1(i1)*fb*sngl(t**2)
             if(.not.(ft.le.0..or.ft.ge.0.))then !NaN catch
             print*,ft, psjetj(0,klas,tmx,q1,q2,scale1,t,xx,sh,j,l)
             stop'NaN in psjet1 fx2' 
             endif
           enddo
           enddo
           fx2=fx2 + dble(a1(i)*ft)/sh**2*z**(1.+delh)
     &          *(1.d0/tmn-1.d0/tmx) !part of t-Jacobien, MUST be inside klas loop
          endif
          endif
c      if(.not.(fx2.le.0..or.fx2.ge.0.))then !NaN catch
c         print*,'NaN in loop 2 in psjet1',fx2,t,tmn,tmx
c         stop'NaN in loop 2 in psjet1'
c      endif
         enddo !klas
        enddo
        enddo
        fx2=fx2*(xmin1-xmax1)/dble(delh)
      endif
      tt1=psuds(q2,l)
      psjet1=sngl((fx1+fx2)/psuds(q2,l))*4*pi**3 
     .             /4   !(t-z-Jacobian)
      tt=psjet1
      if(.not.(tt.le.0..or.tt.ge.0.) !NaN catch
     .   .or.tt.gt.1e35 .or. tt.lt.-1e35 )then !Infinity catch
         write(ifmt,*) 'NaN / Inf catch :   psjet1 = ',tt
         write(ifmt,*) 'psuds,fx1,fx2:',tt1,fx1,fx2
         stop'ERROR 13032019'
      endif
      return
      end

c-----------------------------------------------------------------------
      function psborn(q1,q2,qqcut,ss,j,l,jdis,ipm)
c-----------------------------------------------------------------------
c
c    hard 2->2 parton scattering born cross-section
c     integrated over t from tmin to tmax   
c     including Sudakov on both sides
c
c-----------------------------------------------------------------------
c
c q1 - virtuality cutoff at current end of the ladder;
c q2 - virtuality cutoff at opposite end of the ladder;
c qqcut - p_t cutoff for the born process;
c s - c.m. energy squared for the scattering;
c jdis : -1 - soft gluon interaction below qqcut (-> no integration on t)
c         0 - normal hard above qqcut
c         1 - for DIS
c j - parton type at current end of the ladder 
c l - parton type at opposite end of the ladder 
c def of j,l corresonds to l1,j1  of eolpi arguments (inverted order !!!)
c-----------------------------------------------------------------------
#include "aaa.h"
#include "sem.h"
      common /ar3/   x1(7),a1(7)
      double precision sud0,psbornd,psbornk,psuds,fb,s,qq,t  
      double precision smn,tmn,tmx
      integer ienvi
      common /cienvi/ ienvi
      !real crash(2)
      !icrash=3
      iret=0
      iret0=0
      iret1=0

      psborn=0.

      if(iabs(j).eq.3.and.  noflav(q1).lt.4 ) return ! 1 charm
      if(iabs(j).eq.4.and.  noflav(q1).lt.5 ) return ! 1 bottom
      if(iabs(l).eq.3.and.  noflav(q2).lt.4 ) return ! 2 charm
      if(iabs(l).eq.4.and.  noflav(q2).lt.5 ) return ! 2 bottom

      ienviSave=ienvi
      ienvi=104

      s=dble(ss)
      if(jdis.le.0)then
        qq=dble(max(q1,q2))
      else
        qq=dble(max(q1/4.,q2))
      endif
      if(jdis.ge.0)then
        qq=max(qq,dble(qqcut))
      elseif(qq.gt.qqcut)then
        goto 999
      endif
      qmin=qq 

      j1=j                                  
      l1=l                                 
      call eolpic( j1 , l1 , j2 , l2 )
      call eolpib( j2 , l2 , jl ) 
      call getLimitsEoL(jl,ss,qmin , smin,tmin,tmax,iret0)
      if(iret0.ne.0)goto 999

      if(jdis.ge.0)then
        sud0=psuds(q1,j1)*psuds(q2,l1) !KW1808d
      else
        sud0=1d0      !for jdis<0 put sudakov factor outside to avoid q1/q2 dependence (only max(q1,q2))
      endif
      iret1=1

      psbornd=0.d0
      do kklas=1,igetKlasMaxEoL(jl) !only jl compatible klas values (no evolution)
         klas=igetKlasEoL(jl,kklas)
         qqsi=qq
        if(jdis.ge.0)then ! not PomSat
           call HardScale(klas,qqsi,pt2min,1) !get pt2min
           pt2min=max(0.,pt2min)
           call getBornKin2(klas,s,dble(pt2min), smn, tmn, tmx,iret)
         else ! PomSat
           tmn=tmin
           call getBornPtr(klas,sngl(s),sngl(tmn) , pt2) !compute pt2
           pt2min=0.001
           iret=1
           if(pt2.ge.pt2min)then
             call getBornKin(klas,sngl(s),pt2, sminDmy,tmni,tmxi,iret) !get iret
           endif
         endif
         psbornk=0.d0
         if(iret.eq.0)then
           iret1=0
           imx=7
           mmx=2
           if(jdis.lt.0)then  !if below qqcut, compute for tmn only (not int.)
             imx=1
             mmx=1
             tmx=tmn
           endif
           do i=1,imx
           do m=1,mmx
             t=2.d0*tmn/(1.d0+tmn/tmx-dble(x1(i))*(2*m-3)
     &       *(1.d0-tmn/tmx))
             call getBornPtr(klas,ss,sngl(t) , pt2)  !get pt2
             call HardScale(klas, scale , pt2 ,2) !get scale
             if(jdis.le.0)then
               scale1=scale
             else
               scale1=scale*4.
             endif
             fb=psjetk(ipm,klas,tmx,q1,q2,scale,t,s,j,l,jdis)
             !if(m.eq.1)print*,'CHECK psborn',klas, s,t,fb
             fb=fb*dble(pssalf(scale/qcdlam))**2
             if(jdis.ge.0)then
               fb=fb*psuds(scale1,j1)*psuds(scale,l1)  !cKW2108 Born NOT for -1
               psbornk=psbornk+dble(a1(i))*fb*t**2
             else   !no integration
               psbornk=fb*2d0   !factor 2 to compensate 1/2 from Jacobian hidden somewhere ... (TP20190503)
             endif
           enddo
           enddo
           if(jdis.ge.0)psbornk=psbornk*(1./tmn-1./tmx) !part of t-Jacobien, must be inside klas loop
         endif
         psbornd=psbornd+psbornk
      enddo !klas
      psbornd=psbornd*dble(2.*pi**3)/s**2/sud0*2
     *    /2d0   !part of t-Jacobian

      psborn=sngl(psbornd)
 999  if(ish.ge.7)write(ifch,*)'psborn:jdis,qq,qqcut,iret'
     &            ,psborn,jdis,qq,qqcut,j,l,jl,smin,ss,iret0,iret1
c      if(jdis.lt.0)
c     .print*,'psborn',q1,q2,qqcut,ss,j,l,psborn,iret,iret0,iret1
      ienvi=ienviSave
      return
      end

c-----------------------------------------------------------------------
      function psjeti(ipm,klasxx,tmax,q1,q2,scale,t,xx1,xx2,s,jx,lx
     .               ,jdis)
c-----------------------------------------------------------------------
c
c klasxx > 0 : fast version, regoups individual xsections
c           jx,lx have to be understood as j1,l1 representing classes of xsections
c
c klasxx < 0 : detailed version, computes individual xsections  
c           jx,lx have to be understood as j2,l2 representing individial pair
c
c            (see subroutine initEoL for these definitions)
c
c klasxx = -9999 : goto 99
c           This is needed in case psjeti was not called for klas=klasmax 
c
c-----------------------------------------------------------------------
c           Integrand of psjet (both sides parton ladder cross-section)
c-----------------------------------------------------------------------
c
c      E~qcd_ji * E~qcd_lk * B_ik
c
c        B_ik = psbori = contribution to Born xsection:
c                         dsigmaBorn/d2pt/dy
c                         = s/pi * delta(s+t+u) * 2*pi*alpha**2 /s**2 * B_ik
c
c        E~qcd: at least one emission
c
c klasx - absolute value of klasx : Born kinematics class
c         sign of klas negative -> detailed calculation
c q1  - virtuality cutoff at current end of the ladder
c q2  - virtuality cutoff at opposite end of the ladder
c scale  - factorization scale
c xx1 - feynman x for the first parton for the born process
c xx2 - feynman x for the second parton for the born process
c s   - c.m. energy squared for the born scattering
c t   - invariant variable for the scattering |(p1-p3)**2|,
c-----------------------------------------------------------------------
c j  - parton type at current end of the ladder
c l  - parton type at opposite end of the ladder
c-----------------------------------------------------------------------
c  !!  This function gives non-zero results only for selected choices
c  !!  of j,l with each one representing a group of pairs having the same 
c  !!  cross section. But psjeti returns the cross section for one pair
c  !!  and not the sum over the group as in earlier versions.
c  !!  These "selected choices" of j,l are precicely those whose cross  
c  !!  sections will be tabulated
c-----------------------------------------------------------------------
      implicit none
      integer ipm,klas,klasx,klasxx,jx,lx,kl,jii,lii
      real psjeti,psevi,a
      double precision xx1,xx2,xx1ii,xx2ii,ffborn,s,t,tmax,u,pt2
      double precision xx1Save
#include "sem.h"
      real q1,q2,scale, q1Save,q1ii,q2ii
      real cha1,bot1,cha2,bot2
      integer j,l,jdis,nf,jSave,mp,n1,n2,jdummy,nof,noflav,nfjet,ifmtx
      parameter (nfjet=5)
      double precision ff(-nfjet:nfjet,-nfjet:nfjet) 
      double precision ffk(klasmax,-nfjet:nfjet,-nfjet:nfjet) 
      double precision ak1(-nfjet:nfjet,-nfjet:nfjet)
     .,ak2(-nfjet:nfjet,-nfjet:nfjet)
      integer kflav1(-nfjet:nfjet),kflav2(-nfjet:nfjet)
      real ak1NS,ak1S,ak2NS,ak2S
      real aks1,akns1,aks2,akns2,akg1,akg2,akq1s,akc1s,akb1s,akq1n
     .,akc1n,akb1n,akq2s,akc2s,akb2s,akq2n,akc2n,akb2n
      real gg,gq,qq,qa,qqp,gc,gb,qc,qb,cc,cac,bb,bab,cb
      double precision ggX,gqX,qqX,qaX,qqpX,gcX,gbX,qcX,qbX,ccX,cacX
     .,bbX,babX,cbX
      double precision ggY,gqY,qqY,qaY,qqpY,gcY,gbY,qcY,qbY,ccY,cacY
     .,bbY,babY,cbY
      integer nf0,nf1,nf2,nf3,nf4  ,ifchx,ish,igetIsh
      integer i,k,j1,j2,l1,l2,jl,m,n
      logical lqq,lqa,lqp,lqc,lqb,lcc,lcac,lbb,lbab,lcb
      logical lgq,lgc,lgb
      double precision weightx
      !integer ladderindex
      !integer ival
      !real val5,val6
      real test10(10)
      common /ctest10/ test10 
      ish=igetIsh()

      if(klasxx.eq.-9999)then
        klasx=-klasmax
        klas=klasmax
        goto 99
      else
        klasx=klasxx
      endif

      if(.not.(xx1.gt.0..or.xx1.le.0.))stop'ERROR psjeti xx1=NaN'
      if(.not.(xx2.gt.0..or.xx2.le.0.))stop'ERROR psjeti xx1=NaN'

      if(klasx.lt.0)then !-----------detailed-----------
        j2=jx !initial parton
        l2=lx !    pair
        call eolpib( j2 , l2 , jl )
        call eolpi( jl , j1 , l1 )  
        j=j1 !initial parton
        l=l1 !  pair group
      else !----------------------end detailed-----------
        j=jx !initial parton
        l=lx !  pair group
      endif  

      klas=abs(klasx)
      if(klas.lt.1.or.klas.gt.klasmax)stop'psjeti: wrong klas'

      psjeti=0
      q1ii=q1
      q2ii=q2
      xx1ii=xx1
      xx2ii=xx2
      jii=j
      lii=l

      if(iabs(j).eq.3.and.  noflav(q1).lt.4 ) goto 999 ! 1 charm
      if(iabs(j).eq.4.and.  noflav(q1).lt.5 ) goto 999 ! 1 bottom
      if(iabs(l).eq.3.and.  noflav(q2).lt.4 ) goto 999 ! 2 charm
      if(iabs(l).eq.4.and.  noflav(q2).lt.5 ) goto 999 ! 2 bottom

      gg=0.
      gq=0.
      qq=0.
      qa=0.
      qqp=0.
      gc=0. 
      gb=0. 
      qc=0. 
      qb=0. 
      cc=0. 
      cac=0.
      bb=0. 
      bab=0.
      cb=0. 

      nof=noflav(scale)
      n1=1   !subprocess
      n2=3   !number range

      if(klasx.lt.0)then !-----------detailed-----------

        !flavor restrictions

        do i=-nfjet,nfjet
          kflav1(i)=1
          kflav2(i)=1
          if(nof.lt.abs(i))then
            kflav1(i)=0
            kflav2(i)=0
          endif
        enddo 
        
        !ak matrices

        !-------------------------------------------------
        ! reminder
        !     psevi: 1 1 ... gluon -> gluon
        !            2 1 ... quark -> gluon
        !            1 2 ... gluon -> quark
        !            3 2 ... quark -> quark non singlet
        !            2 2 ... quark -> quark all
        !                      singlet = all - non singlet
        !-------------------------------------------------
        ak1(0,0)=psevi(q1,scale,xx1,1,1)                  
        ak2(0,0)=psevi(q2,scale,xx2,1,1)                  
        ak1(0,1)=psevi(q1,scale,xx1,1,2)/noflav(scale)/2. 
        ak2(0,1)=psevi(q2,scale,xx2,1,2)/noflav(scale)/2. 
        ak1NS=psevi(q1,scale,xx1,3,2)   !nonsinglet 
        ak1S=(psevi(q1,scale,xx1,2,2)-ak1NS)/noflav(scale)/2. !singlet 
        ak2NS=psevi(q2,scale,xx2,3,2)   !nonsinglet 
        ak2S=(psevi(q2,scale,xx2,2,2)-ak2NS)/noflav(scale)/2. !singlet 
        ak1(1,0)=psevi(q1,scale,xx1,2,1)  
        ak2(1,0)=psevi(q2,scale,xx2,2,1)  
        do i=-nfjet, nfjet
        do k=-nfjet, nfjet
          if(i.eq.0)then !g->g,q
            ak1(i,k)=ak1(0,min(1,abs(k)))*kflav1(k)
            ak2(i,k)=ak2(0,min(1,abs(k)))*kflav2(k)
          elseif(k.eq.0)then !q->g
            ak1(i,k)=ak1(min(1,abs(i)),0)
            ak2(i,k)=ak2(min(1,abs(i)),0)
          elseif(i.eq.k)then !q->q
            ak1(i,k)=(ak1NS+ak1S)*kflav1(k)
            ak2(i,k)=(ak2NS+ak2S)*kflav2(k)
          else !q->q'
            ak1(i,k)=ak1S*kflav1(k)
            ak2(i,k)=ak2S*kflav2(k)
          endif
        enddo
        enddo

        !ff matrices

        do i=-nfjet,nfjet
        do k=-nfjet,nfjet
          ff(i,k)=ak1(j2,i)*ak2(l2,k)
          ffk(klas,i,k)=ff(i,k)
        enddo
        enddo

      endif !-----------end detailed-----------

      nf=3                
      nf0= nf*2
      nf1=(nf*2-1)
      nf2=(nf*2-2)
      nf3=(nf*2-3)
      nf4=(nf*2-4)

      jdummy=jdis !not used

      cha1=0.
      bot1=0.
      cha2=0.
      bot2=0. 
      if(nof.ge.4)then
        cha1=1. !consider evolution towards charm on current side 
        cha2=1. !consider evolution towards charm on opposite side 
      endif
      if(nof.ge.5)then
        bot1=1. !consider evolution towards bottom on current side 
        bot2=1. !consider evolution towards bottom on opposite side 
      endif

      !psjeti est symmetrique: map j,l to l,j 
      mp=0
      if(j.eq.1.and.l.eq.0)mp=1        ! j,l = 1,0  
      if(j.eq.3.and.l.eq.0)mp=1        ! j,l = 3,0  
      if(j.eq.4.and.l.eq.0)mp=1        ! j,l = 4,0  
      if(j.eq.3.and.l.eq.1)mp=1        ! j,l = 3,1
      if(j.eq.4.and.l.eq.1)mp=1        ! j,l = 4,1
      if(j.eq.4.and.l.eq.3)mp=1        ! j,l = 4,3
      if(mp.eq.1)then 
        q1Save=q1 
        q1=q2
        q2=q1Save
        xx1Save=xx1
        xx1=xx2
        xx2=xx1Save
        jSave=j
        j=l
        l=jSave
      endif

      if(j.eq.0.and.l.eq.0)then !---------------initial gluon-gluon---------------        

        akg1=psevi(q1,scale,xx1,1,1)                  !gluon contribution
        akg2=psevi(q2,scale,xx2,1,1)                  !gluon contribution
        aks1=psevi(q1,scale,xx1,1,2)/noflav(scale)/2. !KW1811  !singlet contribution per quark
        aks2=psevi(q2,scale,xx2,1,2)/noflav(scale)/2. !KW1811  !singlet contribution per quark
        akq1s=aks1
        akq2s=aks2
        akc1s=aks1*cha1                               !bg1808 charm
        akc2s=aks2*cha2                               !bg1808 charm
        akb1s=aks1*bot1                               !bg1808 bottom
        akb2s=aks2*bot2                               !bg1808 bottom 

        gg=  akg1*akg2
        gq=  (akg1*akq2s+akq1s*akg2)*nf0 
        gc=  (akg2*akc1s+akc2s*akg1)*2.
        gb=  (akg2*akb1s+akb2s*akg1)*2.
        qq=  akq1s*akq2s*nf0
        qa=  akq1s*akq2s*nf0
        qqp= akq1s*akq2s*nf0*(nf-1)*2.
        qc=  (akc2s*akq1s+akc1s*akq2s)*4.*nf 
        qb=  (akb2s*akq1s+akb1s*akq2s)*4.*nf
        cc=  akc1s*akc2s*2.
        cac= akc1s*akc2s*2.
        bab= akb1s*akb2s*2.
        bb=  akb1s*akb2s*2.
        cb=  (akc1s*akb2s+akb1s*akc2s)*4.  

      elseif(j.eq.0)then !-----------------initial gluon-quark-----------------

        if(   j.eq.0.and.l.eq.1     ! j,l = 0,1   
     .  .or.  j.eq.0.and.l.eq.3     ! j,l = 0,3
     .  .or.  j.eq.0.and.l.eq.4     ! j,l = 0,4
     .  )then
          continue !ok
        else
          stop'ERROR 30082018'
        endif 

        aks1=psevi(q1,scale,xx1,1,2)/noflav(scale)/2.  !KW1811     !singlet contribution
        akns2=psevi(q2,scale,xx2,3,2)                           !nonsinglet contribution
        aks2=(psevi(q2,scale,xx2,2,2)-akns2)/noflav(scale)/2. !KW1811  !singlet contribution

        akg1=psevi(q1,scale,xx1,1,1)                !gluon 
        akq1s=aks1                                     !light
        akc1s=aks1*cha1                                !charm
        akb1s=aks1*bot1                                !bottom

        akg2=psevi(q2,scale,xx2,2,1)               !gluon 
        akq2s=aks2                                    !light s
        akq2n=akns2                                   !light ns
        akc2s=aks2*cha2                             !charm s
        akc2n=akns2*cha2                            !charm ns
        akb2s=aks2*bot2                               !bottom s 
        akb2n=akns2*bot2                              !bottom ns

        if(l.eq.1)then                 !  0,1
          akc2n=0
          akb2n=0
        elseif(l.eq.3)then             !  0,3
          akq2n=0
          akb2n=0
        elseif(l.eq.4)then             !  0,4
          akq2n=0
          akc2n=0
        else
          stop'ERROR 28082018'
        endif 

        lgq= l.eq.1  
        lgc= l.eq.3 
        lgb= l.eq.4

        gg=  akg1 *akg2
        gq=  akg1 *(akq2n+akq2s*nf0) + akq1s*nf0*akg2
        gc=  akg1 *(akc2n+akc2s*2)    +   akc1s*2*akg2
        gb=  akg1 *(akb2n+akb2s*2)    +   akb1s*2*akg2
        qq=  akq1s*(akq2n+akq2s*nf0)
        qa=  akq1s*(akq2n+akq2s*nf0)
        qqp= akq1s*(akq2n+akq2s*nf0)*(nf-1)*2.
        qc=  akq1s*2*nf*(akc2n+akc2s*2)+akc1s*2*(akq2n+akq2s*2*nf)  
        qb=  akq1s*2*nf*(akb2n+akb2s*2)+akb1s*2*(akq2n+akq2s*2*nf) 
        cc=  akc1s*(akc2n+akc2s*2)
        cac= akc1s*(akc2n+akc2s*2)
        bb=  akb1s*(akb2n+akb2s*2)
        bab= akb1s*(akb2n+akb2s*2)
        cb=  akc1s*2*(akb2n+akb2s*2)+akb1s*2*(akc2n+akc2s*2)

      else !-----------------initial quark-quark-----------------------
 
        if(   j.eq.1.and.l.eq.1     ! j,l = 1,1   
     .  .or.  j.eq.2.and.l.eq.-2    ! j,l = 2,-2
     .  .or.  j.eq.2.and.l.eq.1     ! j,l = 2,1
     .  .or.  j.eq.1.and.l.eq.3     ! j,l = 1,3
     .  .or.  j.eq.1.and.l.eq.4     ! j,l = 1,4
     .  .or.  j.eq.3.and.l.eq.3     ! j,l = 3,3
     .  .or.  j.eq.3.and.l.eq.-3    ! j,l = 3,-3
     .  .or.  j.eq.4.and.l.eq.4     ! j,l = 4,4
     .  .or.  j.eq.4.and.l.eq.-4    ! j,l = 4,-4
     .  .or.  j.eq.3.and.l.eq.4     ! j,l = 3,4
     .  )then
          continue !ok
        else
          stop'ERROR 28082018'
        endif 

        akns1=psevi(q1,scale,xx1,3,2)                 !nonsinglet contribution
        aks1=(psevi(q1,scale,xx1,2,2)-akns1)/noflav(scale)/2. !singlet contribution
        akns2=psevi(q2,scale,xx2,3,2)                 !nonsinglet contribution
        aks2=(psevi(q2,scale,xx2,2,2)-akns2)/noflav(scale)/2. !singlet contribution

        akg1=psevi(q1,scale,xx1,2,1)                !gluon 
        akq1s=aks1                                    !light s
        akq1n=akns1                                   !light ns
        akc1s=aks1*cha1                             !charm s
        akc1n=akns1*cha1                            !charm ns
        akb1s=aks1*bot1                               !bottom s 
        akb1n=akns1*bot1                              !bottom ns

        akg2=psevi(q2,scale,xx2,2,1)                !gluon 
        akq2s=aks2                                    !light s
        akq2n=akns2                                   !light ns
        akc2s=aks2*cha2                             !charm s
        akc2n=akns2*cha2                            !charm ns
        akb2s=aks2*bot2                               !bottom s 
        akb2n=akns2*bot2                              !bottom ns

        if(abs(j).lt.3)then     
          akc1n=0
          akb1n=0
        elseif(abs(j).eq.3)then
          akq1n=0
          akb1n=0
        elseif(abs(j).eq.4)then
          akq1n=0
          akc1n=0
        else
          stop'ERROR 28082018c'
        endif

        if(abs(l).lt.3)then     
          akc2n=0
          akb2n=0
        elseif(abs(l).eq.3)then
          akq2n=0
          akb2n=0
        elseif(abs(l).eq.4)then
          akq2n=0
          akc2n=0
        else
          stop'ERROR 28082018d'
        endif

        lqq= j.eq.1.and.l.eq.1  
        lqa= j.eq.2.and.l.eq.-2 
        lqp= j.eq.2.and.l.eq.1
        lqc= j.eq.1.and.l.eq.3 
        lqb= j.eq.1.and.l.eq.4
        lcc= j.eq.3.and.l.eq.3
        lcac=j.eq.3.and.l.eq.-3
        lbb= j.eq.4.and.l.eq.4
        lbab=j.eq.4.and.l.eq.-4
        lcb= j.eq.3.and.l.eq.4

        gg = akg1*akg2
        gq = akg1*(akq2n+akq2s*nf0) + (akq1n+akq1s*nf0)*akg2
        gc =  akg1*(akc2n+akc2s*2) + (akc1n+akc1s*2)*akg2
        gb =  akg1*(akb2n+akb2s*2)+(akb1n+akb1s*2)*akg2
        if(lqq)then
          qq      =(akq1n+akq1s)*(akq2n+akq2s) + akq1s*nf1 *akq2s
        elseif(lqa)then
          qq      =akq1s*(akq2n+akq2s) +  (akq1n+akq1s*nf1)*akq2s
        elseif(lqp)then
          qq      =akq1s*(akq2n+akq2s)+ (akq1n+akq1s*nf1)*akq2s
        elseif(lqc)then
          qq      =     (akq1n+akq1s*nf0  ) *akq2s
        elseif(lqb)then
          qq      =     (akq1n+akq1s*nf0  ) *akq2s
        else!all w/o light quarks
          qq      = akq1s*nf0 * akq2s 
        endif
        if(lqq)then
          qa      =akq1s*(akq2n+akq2s) +  (akq1n+akq1s*nf1)*akq2s
        elseif(lqa)then
          qa      =(akq1n+akq1s)*(akq2n+akq2s) + akq1s*nf1 *akq2s
        elseif(lqp)then
          qa      =akq1s*(akq2n+akq2s) +  (akq1n+akq1s*nf1)*akq2s
        elseif(lqc)then
          qa      =     (akq1n+akq1s* nf0   ) *akq2s
        elseif(lqb)then
          qa      =     (akq1n+akq1s* nf0   ) *akq2s
        else! all w/o light quarks
          qa    = akq1s*nf0 * akq2s
        endif 
        if(lqq.or.lqa)then
          qqp    =(akq1n+akq1s*2)*akq2s*nf2+akq1s*nf2*(akq2n+akq2s*nf2)
        elseif(lqp)then
          qqp = (akq1n+akq1s*nf2)*(akq2n+akq2s*nf2) + akq1s*2*akq2s*nf2
        elseif(lqc.or.lqb)then
          qqp =    (akq1n+akq1s*nf0) * akq2s*nf2
        else!for all w/o light quarks
          qqp     = akq1s*nf0 * akq2s*nf2  
        endif
        qc = (akq1n+akq1s*nf0)*(akc2n+akc2s*2)  !qc
     .        +(akc1n+akc1s*2)*(akq2n+akq2s*nf0)  !cq
        qb = (akq1n+akq1s*nf0)*(akb2n+akb2s*2)
     .      +(akb1n+akb1s*2)*(akq2n+akq2s*nf0) 
        if(lqc)then
          cc= akc1s*(akc2n+akc2s*2)
        elseif(lcc)then
          cc= (akc1n+akc1s)*(akc2n+akc2s)+akc1s*akc2s
        elseif(lcac)then
          cc= (akc1n+akc1s)*akc2s+akc1s*(akc2n+akc2s)
        elseif(lcb)then
          cc= (akc1n+akc1s*2)*akc2s
        else!all w/o charm
          cc = akc1s*2 * akc2s 
        endif 
        if(lqc)then
          cac= akc1s*(akc2n+akc2s*2)
        elseif(lcac)then
          cac= (akc1n+akc1s)*(akc2n+akc2s)+akc1s*akc2s
        elseif(lcc)then
          cac= (akc1n+akc1s)*akc2s+akc1s*(akc2n+akc2s)
        elseif(lcb)then
          cac= (akc1n+akc1s*2)*akc2s
        else!for all w/o charm
          cac = akc1s*2 * akc2s 
        endif 
        if(lqb)then
          bb= akb1s*(akb2n+akb2s*2)
        elseif(lbb)then
          bb= (akb1n+akb1s)*(akb2n+akb2s)+akb1s*akb2s
        elseif(lbab)then
          bb= (akb1n+akb1s)*akb2s+akb1s*(akb2n+akb2s)
        elseif(lcb)then
          bb= akb1s*(akb2n+akb2s*2)
        else!for all w/o bottom
          bb = akb1s*2 * akb2s 
        endif
        if(lqb)then
          bab= akb1s*(akb2n+akb2s*2)
        elseif(lbb)then
          bab= (akb1n+akb1s)*akb2s+akb1s*(akb2n+akb2s)
        elseif(lbab)then
          bab= (akb1n+akb1s)*(akb2n+akb2s)+akb1s*akb2s
        elseif(lcb)then
          bab= akb1s*(akb2n+akb2s*2)
        else !for all w/o bottom
          bab = akb1s*2 * akb2s 
        endif
        cb = (akc1n+akc1s*2)*(akb2n+akb2s*2)
     .        +(akb1n+akb1s*2)*(akc2n+akc2s*2)

      endif

      psjeti=ffborn(ipm,klas,tmax,s,t,gg,gq,qq,qa ,qqp,n1,n2,nf) +
     *       ffborn(ipm,klas,tmax,s,t,0.,gc,0.,0. ,qc ,n1,n2,4 ) +
     *       ffborn(ipm,klas,tmax,s,t,0.,gb,0.,0. ,qb ,n1,n2,5 ) +
     *       ffborn(ipm,klas,tmax,s,t,0.,0.,cc,cac,0. ,n1,n2,6 ) +
     *       ffborn(ipm,klas,tmax,s,t,0.,0.,bb,bab,cb ,n1,n2,7 )

      if(ish.ge.9)then
        call getMonitorFileIndex(ifmtx)
        write(ifmtx,'(a,2i3,3e12.3,3x,$)')'psjeti = ',j,l,psjeti,s,t 
        write(ifmtx,'(5e12.3)')gg,gq,qq,qa,qqp
        !write(ifmtx,'(3x,9e12.3)')gc,qc,gb,qb,cc,cac,bb,bab,cb
      endif

      call getBornPtr2(klas,s,t , pt2) 
      call getBornU2(klas,s,t,u)
      psjeti = psjeti ! / abs(u-t) !* u * t / pt2 / (u+t)

      call evTest1(gg,gq,gc,gb,qq,qa,qqp,qc,qb,cc,cac,bb,bab,cb)

      if(klasx.lt.0)then !-------------------detailed------------------
        call setFFtypeZero(ggX,gqX,gcX,gbX,qqX,qaX,qqpX,qcX,qbX
     .             ,ccX,cacX,bbX,babX,cbX) !set all values to zero
        do m=-nfjet,nfjet ! in Born parton
        do n=-nfjet,nfjet ! in Born parton
          call setFFtypeZero(                       !set all values to zero
     .      ggY,gqY,gcY,gbY,qqY,qaY,qqpY,qcY,qbY,ccY,cacY,bbY,babY,cbY) 
          call getFFtype(ff,m,n,                    !given ff,m,n: get gg,gq etc 
     .      ggY,gqY,gcY,gbY,qqY,qaY,qqpY,qcY,qbY,ccY,cacY,bbY,babY,cbY) 
          call getFFtypeSngl(                          !trafo dble to sngl
     .      ggY,gqY,gcY,gbY,qqY,qaY,qqpY,qcY,qbY,ccY,cacY,bbY,babY,cbY,
     .      gg,gq,gc,gb,qq,qa,qqp,qc,qb,cc,cac,bb,bab,cb)
          !call setBornIn(m,n) !store Born-in flavors, just for info, to be known in psabori
          !call setIshy(11)
          weightx =
     .       ffborn(ipm,klas,tmax,s,t,gg,gq,qq,qa ,qqp,n1,n2,nf) +
     .       ffborn(ipm,klas,tmax,s,t,0.,gc,0.,0. ,qc ,n1,n2,4 ) +
     .       ffborn(ipm,klas,tmax,s,t,0.,gb,0.,0. ,qb ,n1,n2,5 ) +
     .       ffborn(ipm,klas,tmax,s,t,0.,0.,cc,cac,0. ,n1,n2,6 ) +
     .       ffborn(ipm,klas,tmax,s,t,0.,0.,bb,bab,cb ,n1,n2,7 )
          !call setIshy(0)
          !+++++++++++++++
          !call getLadderindex(ladderindex)
          !if((klas.eq.1.or.klas.eq.4).and.n.eq.0.and.m.eq.0)
          !.print*,'CHECKlad',ladderindex,klas,weightx
          !+++++++++++++++
          call setBornInWeight(klas,m,n,weightx)    !fill weight table
          call sumFFtype(  !add Y-value to X-value
     .      ggX,gqX,gcX,gbX,qqX,qaX,qqpX,qcX,qbX,ccX,cacX,bbX,babX,cbX,
     .      ggY,gqY,gcY,gbY,qqY,qaY,qqpY,qcY,qbY,ccY,cacY,bbY,babY,cbY)
        enddo
        enddo
        call setTmax(klas,tmax) !store tmax for given klas   
        call evTest2(3,j,l,klas
     .  ,ggX,gqX,gcX,gbX,qqX,qaX,qqpX,qcX,qbX,ccX,cacX,bbX,babX,cbX)
      endif !----------------end detailed------------------

  99  continue

      if(klasx.lt.0)then !------------- detailed------------------
        if(klas.eq.klasmax)then
          test10(1)=4
          call getBornInRdm(kl,m,n)      !generate kl,m,n
          !++++++++++++++
          !call getLadderindex(ladderindex)
          !ival=0
          !if(ladderindex.eq.5.and.kl.eq.4)ival=5
          !if(ladderindex.lt.5.and.kl.eq.4)ival=6
          !if(ival.gt.0)call incVparam(ival,1.)
          !call getVparam(5,val5)
          !call getVparam(6,val6)
          !print*,'CHECKlad-----------------------------',ladderindex,kl
          ! .,nint(val5),nint(val6)
          !++++++++++++++
          call getTmax(kl,tmax)            !get stored tmax for kl  
          do i=-nfjet,nfjet
          do k=-nfjet,nfjet
            ff(i,k)=ffk(kl,i,k)            !get ff for kl
          enddo
          enddo
          call setBornIn(m,n)                !store m,n
          call getFFtype(ff,m,n,              !given ff,m,n: get gg,gq etc 
     .      ggY,gqY,gcY,gbY,qqY,qaY,qqpY,qcY,qbY,ccY,cacY,bbY,babY,cbY) 
          call getFFtypeSngl(                   !trafo dble to sngl
     .      ggY,gqY,gcY,gbY,qqY,qaY,qqpY,qcY,qbY,ccY,cacY,bbY,babY,cbY,
     .      gg,gq,gc,gb,qq,qa,qqp,qc,qb,cc,cac,bb,bab,cb)
          call setBornOutWeightsZero !initialize weight table
          call getCheckFileIndex(ifchx)
          if(ish.eq.-1007)rewind(ifchx)
          !compute Born-out weight table (using m,n via getBornIn)
          a=ffborn(ipm,kl,tmax,s,t,gg,gq,qq,qa ,qqp,n1,n2,nf) +
     .      ffborn(ipm,kl,tmax,s,t,0.,gc,0.,0. ,qc ,n1,n2,4 ) +
     .      ffborn(ipm,kl,tmax,s,t,0.,gb,0.,0. ,qb ,n1,n2,5 ) +
     .      ffborn(ipm,kl,tmax,s,t,0.,0.,cc,cac,0. ,n1,n2,6 ) +
     .      ffborn(ipm,kl,tmax,s,t,0.,0.,bb,bab,cb ,n1,n2,7 )        !psjeti
          !table will be used via getBornOutRdm to generate the Born-out flavors
          !~~~
          ! call checkBornOutWeights(4,kl,m,n
          !.,gg,gq,gc,gb,qq,qa,qqp,qc,qb,cc,cac,bb,bab,cb)
          !~~~
        endif
      endif !-------------------end detailed------------------


    !print*,'---psjeti---',s,t,gg,ffborn(ipm,klas,tmax,s,t,gg,gq,qq,qa ,qqp,n1,n2,nf)

cc      if(psjeti.lt.0.)then !negative cross section
cc        call getMonitorFileIndex(ifmtx)
cc        write(ifmtx,*) 'neg CS :   psjeti = ',psjeti
cc        write(ifmtx,*) 'q1,q2,qt : ',q1,q2,qt
cc        write(ifmtx,*) 't,xx1,xx2 : ',t,xx1,xx2
cc        write(ifmtx,*) 's :',s
cc        write(ifmtx,*) 'j,l,jdis : ',j,l,jdis
cc        write(ifmtx,*) 'qc :',qc
cc        stop'ERROR 18122018b'
cc      endif

 999  continue

      q1=q1ii
      q2=q2ii
      xx1=xx1ii
      xx2=xx2ii
      j=jii
      l=lii

      return
      end !psjeti

c-----------------------------------------------------------------------
      function psjetj(ipm,klasxx,tmax,q1,q2,scale,t,xx,s,jx,lx)
c-----------------------------------------------------------------------
c
c klasxx > 0 : fast version, regoups individual xsections
c           jx,lx have to be understood as j1,l1 representing classes of xsections
c
c klasxx < 0 : detailed version, computes individual xsections  
c           jx,lx have to be understood as j2,l2 representing individial pair
c
c            (see subroutine initEoL for these definitions)
c
c klasxx = -9999 : goto 99
c           This is needed in case psjetj was not called for klas=klasmax 
c
c-----------------------------------------------------------------------
c            Integrand of psjet1 (ordered ladder cross-section) 
c-----------------------------------------------------------------------
c q1 - virtuality cutoff at current end of the ladder,
c q2 - virtuality cutoff at opposite end of the ladder,
c scale - born process scale,
c t  - invariant variable for the scattering |(p1-p3)**2|,
c xx - feynman x for the first parton for the born process
c s  - c.m. energy squared for the born scattering,
c j  - parton type at current end of the ladder
c l  - parton type at opposite end of the ladder
c-----------------------------------------------------------------------
c  !!  This function gives non-zero results only for selected choices
c  !!  of j,l with each one representing a group of pairs having the same 
c  !!  cross section. But psjetj returns the cross section for one pair
c  !!  and not the sum over the group as in earlier versions.
c  !!  These "selected choices" of j,l are precicely those whose cross  
c  !!  sections will be tabulated
c-----------------------------------------------------------------------
      implicit none
      integer ipm,klas,klasx,klasxx,kl
      real psjetj,psevi,a
      double precision tmax,u,pt2
      real q1,q2,scale
      integer j,l,n1,n2,nf,nof,noflav,nfjet,i,k,jx,lx,ifmtx
#include "sem.h"
      double precision xx,ffborn,s,t
      logical cha,bot,cha2,bot2
      real gg,gq,qq,qa,qqp,gc,gb,qc,qb,cc,cac,bab,bb,cb
      double precision ggX,gqX,qqX,qaX,qqpX,gcX,gbX,qcX,qbX,ccX,cacX
     .,bbX,babX,cbX
      double precision ggY,gqY,qqY,qaY,qqpY,gcY,gbY,qcY,qbY,ccY,cacY
     .,bbY,babY,cbY
      parameter (nfjet=5)
      double precision ff(-nfjet:nfjet,-nfjet:nfjet) 
      double precision ffk(klasmax,-nfjet:nfjet,-nfjet:nfjet) 
      double precision ak1(-nfjet:nfjet,-nfjet:nfjet)
      integer kflav1(-nfjet:nfjet)
      real ak1NS,ak1S
      double precision checkWeightSums(2,-klasmax:klasmax)
      common /ccheckWeightSums/ checkWeightSums
      integer j1,j2,l1,l2,jl,m,n
      double precision weightx
      real ak11,ak21,ak12,akns,aks
      integer nf0,nf1,nf2,nf3,nf4  ,ifchx,ish,igetIsh
      !integer ladderindex
      !integer ival
      !real val5,val6
      real crash(2)
      integer icrash
      real test10(10)
      common /ctest10/ test10 
      icrash=3
      ish=igetIsh()

      if(klasxx.eq.-9999)then
        klasx=-klasmax
        klas=klasmax
        goto 99 
      else
        klasx=klasxx
      endif

      if(klasx.lt.0)then !-----------detailed-----------
        j2=jx !initial parton
        l2=lx !    pair
        call eolpib( j2 , l2 , jl )
        call eolpi( jl , j1 , l1 )  
        j=j1 !initial parton
        l=l1 !  pair group
      else !-----------------------end detailed-----------
        j=jx !initial parton
        l=lx !  pair group
      endif  

      klas=abs(klasx)
      n1=1   !subprocess
      n2=3   !number range
      nof=noflav(scale)

      if(.not.(xx.gt.0..or.xx.le.0.))crash(icrash)=1 !NaN catch

      psjetj=0

      if(iabs(j).eq.3.and.  noflav(q1).lt.4 ) goto 999 ! 1 charm
      if(iabs(j).eq.4.and.  noflav(q1).lt.5 ) goto 999 ! 1 bottom
      if(iabs(l).eq.3.and.  noflav(q2).lt.4 ) goto 999 ! 2 charm
      if(iabs(l).eq.4.and.  noflav(q2).lt.5 ) goto 999 ! 2 bottom

      cha=.false.
      bot=.false.
      if(nof.ge.4)then
        cha=.true. !consider evolution towards charm on current side 
      endif
      if(nof.ge.5)then
        bot=.true. !consider evolution towards bottom on current side 
      endif

      !KW2006 cha2, bot2 added
      cha2=.false.
      bot2=.false.
      if(nof.ge.4)cha2=.true. 
      if(nof.ge.5)bot2=.true.

      if(klasx.lt.0)then !-----------detailed-----------

        !flavor restrictions

        do i=-nfjet,nfjet
          kflav1(i)=1
          if(nof.lt.abs(i))then
            kflav1(i)=0
          endif
        enddo 
        
        !ak matrices

        !-------------------------------------------------
        ! reminder
        !     psevi: 1 1 ... gluon -> gluon
        !            2 1 ... quark -> gluon
        !            1 2 ... gluon -> quark
        !            3 2 ... quark -> quark non singlet
        !            2 2 ... quark -> quark all
        !                      singlet = all - non singlet
        !-------------------------------------------------
        ak1(0,0)=psevi(q1,scale,xx,1,1)                  
        ak1(0,1)=psevi(q1,scale,xx,1,2)/noflav(scale)/2. 
        ak1NS=psevi(q1,scale,xx,3,2)   !nonsinglet 
        ak1S=(psevi(q1,scale,xx,2,2)-ak1NS)/noflav(scale)/2. !singlet 
        ak1(1,0)=psevi(q1,scale,xx,2,1)  
        do i=-nfjet, nfjet
        do k=-nfjet, nfjet
          if(i.eq.0)then !g->g,q
            ak1(i,k)=ak1(0,min(1,abs(k)))*kflav1(k)
          elseif(k.eq.0)then !q->g
            ak1(i,k)=ak1(min(1,abs(i)),0)
          elseif(i.eq.k)then !q->q
            ak1(i,k)=(ak1NS+ak1S)*kflav1(k)
          else !q->q'
            ak1(i,k)=ak1S*kflav1(k)
          endif
        enddo
        enddo

        !ff matrices

        do i=-nfjet,nfjet
        do k=-nfjet,nfjet
          if(k.eq.l2)then
            ff(i,k)=ak1(j2,i)*kflav1(k)  !KW2006 *kflav1(k) added
          else
            ff(i,k)=0
          endif 
          ffk(klas,i,k)=ff(i,k)
        enddo
        enddo
        if(klas.eq.2)then
         test10(6)=jx
         test10(7)=lx
         test10(8)=scale
         test10(9)=0
         do i=-3,3 
         test10(9)=test10(9)+abs(ff(i,4))+abs(ff(i,-4))
         test10(9)=test10(9)+abs(ff(4,i))+abs(ff(-4,i))
         enddo
        endif

      endif !-----------end detailed-----------

      gg=0.
      gq=0.
      qq=0.
      qa=0.
      qqp=0.
      gc=0. 
      gb=0. 
      qc=0. 
      qb=0. 
      cc=0. 
      cac=0.
      bab=0.
      bb=0. 
      cb=0. 

      nf=3 !light quark (u,d,s)
      nf0= nf*2
      nf1=(nf*2-1)
      nf2=(nf*2-2)
      nf3=(nf*2-3)
      nf4=(nf*2-4)


      ! j = current      <=== emission
      ! l = opposite    

      ak11=psevi(q1,scale,xx,1,1)
      ak12=psevi(q1,scale,xx,1,2)/noflav(scale)/2.  !KW1811   !contribution per quark
      ak21=psevi(q1,scale,xx,2,1)                          !gluon contribution
      akns=psevi(q1,scale,xx,3,2)                          !nonsinglet contribution
      aks=(psevi(q1,scale,xx,2,2)-akns)/noflav(scale)/2. !KW1811  !singlet contribution per quark

      !Remark (KW2006)
      !In principle we use indices j,l as a code for pair groups, 
      !so 2,1 stands for two light (anti)quarks of different flavor
      !which covers many pairs like ud,us,du,ds,ubardbar,dbarubar...
      !wheras 1,2 means nothing, not defined
      !BUT since we exchange j and l in some cases, we accept here 
      !codes like 1,2,  whichs amount to the same as 2,1 

      if(l.eq.0)then

        if(j.eq.0)then               !KW1808c  j,l = 0,0

          !always
            gg=ak11
            gq=ak12*nf0 
          if(cha)then
            gc=ak12*2 
          endif 
          if(bot)then
            gb=ak12*2 
          endif 

        elseif(j.eq.1)then           !KW1808c  j,l = 1,0

          !always
            gg=ak21
            gq=akns+aks*nf0 
          if(cha)then
            gc=aks*2 
          endif
          if(bot)then
            gb=aks*2 
          endif 

        elseif(j.eq.3)then           !KW1808c  j,l = 3,0 

          !always
            gg=ak21
            gq=aks*nf0 
          if(cha)then
            gc=akns+aks*2
          endif 
          if(bot)then
            gb=aks*2
          endif 

        elseif(j.eq.4)then           !KW1808c  j,l = 4,0 

          !always
            gg=ak21
            gq=aks*nf0 
          if(cha)then
            gc=aks*2
          endif
          if(bot)then
            gb=akns+aks*2
          endif 

        else
          stop'ERROR 25082018'
        endif

      elseif(j.eq.0)then

        if(l.eq.1) then              !KW1808c  j,l = 0,1    (stands for g,q)

          !always
            gq=ak11
            qq=ak12
            qa=ak12
            qqp=ak12*nf2 
          if(cha)then
            qc=ak12*2 
          endif 
          if(bot)then
            qb=ak12*2 
          endif 

        elseif(l.eq.3.and.cha2) then    !KW1808c  j,l = 0,3  !KW2006 .and.cha2 added

          !always
            gc=ak11
            qc=ak12*nf0
          if(cha)then
            cc=ak12
            cac=cc
          endif 
          if(bot)then
            cb=ak12*2.
          endif 

        elseif(l.eq.4.and.bot2) then   !KW1808c  j,l = 0,4   !KW2006 .and.bot2 added

          !always
            gb=ak11
            qb=ak12*nf0
          if(bot)then
            bb=ak12
            bab=bb
          endif 
          if(cha)then
            cb=ak12*2.
          endif 

        endif

      elseif(j.eq.1.and.l.eq.1)then  !KW1808c  j,l = 1,1  

        !always
          gq=ak21
          qq=akns+aks     
          qa=aks         
          qqp=aks*(nf-1)*2.
        if(cha)then
          qc=aks*2
        endif 
        if(bot)then
          qb=aks*2
        endif 

      elseif((j.eq.2.and.l.eq.-2).or.(j.eq.-2.and.l.eq.2))then !KW1808c  j,l = 2,-2  

        !always
          gq=ak21
          qq=aks
          qa=(akns+aks)
          qqp=aks*nf2
        if(cha)then
          qc=aks*2
        endif 
        if(bot)then
          qb=aks*2
        endif 

      elseif((j.eq.2.and.l.eq.1).or.(j.eq.1.and.l.eq.2))then  !KW1808c  j,l = 2,1 

        !always
          gq=ak21
          qq=aks
          qa=aks
          qqp=akns+aks*nf2
        if(cha)then
          qc=aks*2
        endif 
        if(bot)then
          qb=aks*2
        endif 

      elseif(j.eq.1.and.l.eq.3.and.cha2)then  !KW1808c  j,l = 1,3  !KW2006 .and.cha2 added

        !always
          gc=ak21
          qc=akns+aks*nf0
        if(cha)then
          cc=aks
          cac=aks
        endif 
        if(bot)then
          cb=aks*2 
        endif 

      elseif(j.eq.3.and.l.eq.1)then  !KW1808c  j,l = 3,1 

        !always
        gq=ak21
        qq=aks
        qa=aks
        qqp=aks*nf2
        if(cha)then
          qc=akns+aks*2.
        endif
        if(bot)then
          qb=aks*2 
        endif 

      elseif(j.eq.1.and.l.eq.4.and.bot2)then  !KW1808c  j,l = 1,4  !KW2006 .and.bot2 added

        !always
          gb=ak21
          qb=akns+aks*nf0
        if(bot)then
          bb=aks
          bab=aks
        endif 
        if(cha)then
          cb=aks*2 
        endif 

      elseif(j.eq.4.and.l.eq.1)then  !KW1808c  j,l = 4,1 

        !always
          gq=ak21
          qq=aks
          qa=aks
          qqp=aks*nf2
        if(bot)then
          qb=akns+aks*2.
        endif 
        if(cha)then
          qc=aks*2 
        endif 

      elseif(j.eq.3.and.l.eq.3.and.cha2)then  !KW1808c  j,l = 3,3  !KW2006 .and.cha2 added

        !always
          gc=ak21
          qc=aks*nf0
        if(cha)then
          cc=akns+aks
          cac=aks
        endif 
        if(bot)then
          cb=aks*2 
        endif 

      elseif((j.eq.3.and.l.eq.-3.or.j.eq.-3.and.l.eq.3).and.cha2)then !KW1808c  j,l = 3,-3  !KW2006 .and.cha2 added

        !always
          gc=ak21
          qc=aks*nf0
        if(cha)then
          cac=akns+aks
          cc=aks
        endif 
        if(bot)then
          cb=aks*2 
        endif 
        if(klas.eq.2)then
         test10(2)=j
         test10(3)=l
         test10(4)=scale
         test10(5)=abs(gc)+abs(qc)
        endif

      elseif(j.eq.4.and.l.eq.4.and.bot2)then  !KW1808c  j,l = 4,4 !KW2006 .and.bot2 added

        !always
          gb=ak21
          qb=aks*nf0
        if(bot)then
          bb=akns+aks
          bab=aks
        endif 
        if(cha)then
          cb=aks*2 
        endif 

      elseif((j.eq.4.and.l.eq.-4.or.j.eq.-4.and.l.eq.4).and.bot2)then !KW1808c  j,l = 4,-4 !KW2006 .and.bot2 added

        !always
          gb=ak21
          qb=aks*nf0
        if(bot)then
          bab=akns+aks
          bb=aks
        endif 
        if(cha)then
          cb=aks*2 
        endif 

      elseif(j.eq.4.and.l.eq.3.and.cha2)then  !KW1808c  j,l = 4,3  !KW2006 .and.cha2 added

        !always
          gc=ak21
          qc=aks*nf0
        if(bot)then
          cb=akns+aks*2
        endif 
        if(cha)then
          cc=aks
          cac=aks 
        endif 

      elseif(j.eq.3.and.l.eq.4.and.bot2)then  !KW1808c  j,l = 3,4 !KW2006 .and.bot2 added

        !always
          gb=ak21
          qb=aks*nf0
        if(cha)then
          cb=akns+aks*2
        endif 
        if(bot)then
          bb=aks
          bab=aks 
        endif 

      endif

      psjetj=ffborn(ipm,klas,tmax,s,t,gg,gq,qq,qa ,qqp,n1,n2,nf) +
     *       ffborn(ipm,klas,tmax,s,t,0.,gc,0.,0. ,qc ,n1,n2,4 ) +
     *       ffborn(ipm,klas,tmax,s,t,0.,gb,0.,0. ,qb ,n1,n2,5 ) +
     *       ffborn(ipm,klas,tmax,s,t,0.,0.,cc,cac,0. ,n1,n2,6 ) +
     *       ffborn(ipm,klas,tmax,s,t,0.,0.,bb,bab,cb ,n1,n2,7 )
      if(klasx.gt.0)then
      !checkWeightSums(1,klasx)=checkWeightSums(1,klasx)+abs(psjetj)
      !checkWeightSums(2,klasx)=checkWeightSums(2,klasx)+1
      endif

      if(ish.ge.9)then
        call getMonitorFileIndex(ifmtx)
        write(ifmtx,'(a,2i3,3e12.3,3x,$)')'psjetj = ',j,l,psjetj,s,t !Sudakov not included
        !write(ifmtx,'(5e12.3)')gg,gq,qq,qa,qqp
        write(ifmtx,'(3x,9e12.3)')gc,qc,gb,qb,cc,cac,bb,bab,cb
      endif

      call getBornPtr2(klas,s,t , pt2) 
      call getBornU2(klas,s,t,u)
      psjetj = psjetj ! / abs(u-t) * u * t / pt2 / (u+t)

      call evTest1(gg,gq,gc,gb,qq,qa,qqp,qc,qb,cc,cac,bb,bab,cb)

      if(klasx.lt.0)then !------------- detailed------------------
        !the values gg,gq etc do not depend on klas, so in principle they could 
        !be computed and tested once (if klas.eq.1 ...), and then stored for (a bit) more efficiency
        call setFFtypeZero(ggX,gqX,gcX,gbX,qqX,qaX,qqpX,qcX,qbX
     .             ,ccX,cacX,bbX,babX,cbX) !set all values to zero
        do m=-nfjet,nfjet ! in Born parton
        do n=-nfjet,nfjet ! in Born parton
          call setFFtypeZero(                       !set all values to zero
     .      ggY,gqY,gcY,gbY,qqY,qaY,qqpY,qcY,qbY,ccY,cacY,bbY,babY,cbY) 
          call getFFtype(ff,m,n,                    !given ff,m,n: get gg,gq etc 
     .      ggY,gqY,gcY,gbY,qqY,qaY,qqpY,qcY,qbY,ccY,cacY,bbY,babY,cbY) 
          call getFFtypeSngl(                          !trafo dble to sngl
     .      ggY,gqY,gcY,gbY,qqY,qaY,qqpY,qcY,qbY,ccY,cacY,bbY,babY,cbY,
     .      gg,gq,gc,gb,qq,qa,qqp,qc,qb,cc,cac,bb,bab,cb)
          weightx =
     .       ffborn(ipm,klas,tmax,s,t,gg,gq,qq,qa ,qqp,n1,n2,nf) +
     .       ffborn(ipm,klas,tmax,s,t,0.,gc,0.,0. ,qc ,n1,n2,4 ) +
     .       ffborn(ipm,klas,tmax,s,t,0.,gb,0.,0. ,qb ,n1,n2,5 ) +
     .       ffborn(ipm,klas,tmax,s,t,0.,0.,cc,cac,0. ,n1,n2,6 ) +
     .       ffborn(ipm,klas,tmax,s,t,0.,0.,bb,bab,cb ,n1,n2,7 )
             !+++++++++++++++
             !call getLadderindex(ladderindex)
             !if((klas.eq.1.or.klas.eq.4).and.n.eq.0.and.m.eq.0)
             !.print*,'CHECKlad',ladderindex,klas,weightx
             !+++++++++++++++
          call setBornInWeight(klas,m,n,weightx)    !fill weight table
          !checkWeightSums(1,klasx)=checkWeightSums(1,klasx)+abs(weightx)
          call sumFFtype(  !add Y-value to X-value
     .      ggX,gqX,gcX,gbX,qqX,qaX,qqpX,qcX,qbX,ccX,cacX,bbX,babX,cbX,
     .      ggY,gqY,gcY,gbY,qqY,qaY,qqpY,qcY,qbY,ccY,cacY,bbY,babY,cbY)
        enddo
        enddo
        !checkWeightSums(2,klasx)=checkWeightSums(2,klasx) +1
        !. +abs(ff( 0,4))+abs(ff( 0,-4))+abs(ff(4, 0))+abs(ff(-4, 0))
        !. +abs(ff( 1,4))+abs(ff( 1,-4))+abs(ff(4, 1))+abs(ff(-4, 1))
        !. +abs(ff( 2,4))+abs(ff( 2,-4))+abs(ff(4, 2))+abs(ff(-4, 2))
        !. +abs(ff( 3,4))+abs(ff( 3,-4))+abs(ff(4, 3))+abs(ff(-4, 3))
        !. +abs(ff(-1,4))+abs(ff(-1,-4))+abs(ff(4,-1))+abs(ff(-4,-1))
        !. +abs(ff(-2,4))+abs(ff(-2,-4))+abs(ff(4,-2))+abs(ff(-4,-2))
        !. +abs(ff(-3,4))+abs(ff(-3,-4))+abs(ff(4,-3))+abs(ff(-4,-3))
        call setTmax(klas,tmax) !store tmax    
        call evTest2(2,j,l,klas
     .  ,ggX,gqX,gcX,gbX,qqX,qaX,qqpX,qcX,qbX,ccX,cacX,bbX,babX,cbX)
      endif !----------------end detailed------------------

  99  continue

      if(klasx.lt.0)then !------------- detailed------------------
        if(klas.eq.klasmax)then
          test10(1)=23
          call getBornInRdm(kl,m,n)         !generate kl,m,n
          !+++++++++++++++++
          !call getLadderindex(ladderindex)
          !ival=0
          !if(ladderindex.eq.5.and.kl.eq.4)ival=5
          !if(ladderindex.lt.5.and.kl.eq.4)ival=6
          !if(ival.gt.0)call incVparam(ival,1.)
          !call getVparam(5,val5)
          !call getVparam(6,val6)
          !print*,'CHECKlad-----------------------------',ladderindex,kl
          !.,nint(val5),nint(val6)
          !+++++++++++++++++
          call getTmax(kl,tmax)            !get stored tmax for kl    
          do i=-nfjet,nfjet
          do k=-nfjet,nfjet
            ff(i,k)=ffk(kl,i,k)            !get ff for kl
          enddo
          enddo
          call setBornIn(m,n)               !store Born-In flavors   
          call getFFtype(ff,m,n,              !given ff,m,n: get gg,gq etc 
     .      ggY,gqY,gcY,gbY,qqY,qaY,qqpY,qcY,qbY,ccY,cacY,bbY,babY,cbY) 
          call getFFtypeSngl(                   !trafo dble to sngl
     .      ggY,gqY,gcY,gbY,qqY,qaY,qqpY,qcY,qbY,ccY,cacY,bbY,babY,cbY,
     .      gg,gq,gc,gb,qq,qa,qqp,qc,qb,cc,cac,bb,bab,cb)
          call setBornOutWeightsZero      !initialize weight table
          call getCheckFileIndex(ifchx)
          if(ish.eq.-1007)rewind(ifchx)
          !compute Born-out weight table (using m,n via getBornIn)
          a=ffborn(ipm,kl,tmax,s,t,gg,gq,qq,qa ,qqp,n1,n2,nf) +
     .      ffborn(ipm,kl,tmax,s,t,0.,gc,0.,0. ,qc ,n1,n2,4 ) +
     .      ffborn(ipm,kl,tmax,s,t,0.,gb,0.,0. ,qb ,n1,n2,5 ) +
     .      ffborn(ipm,kl,tmax,s,t,0.,0.,cc,cac,0. ,n1,n2,6 ) +
     .      ffborn(ipm,kl,tmax,s,t,0.,0.,bb,bab,cb ,n1,n2,7 )       !psjetj
          !table will be used via getBornOutRdm to generate the Born-out flavors
          !~~~
      call checkBornOutWeights(23,kl,m,n
     .,gg,gq,gc,gb,qq,qa,qqp,qc,qb,cc,cac,bb,bab,cb)
          !~~~
        endif !klas=klasmax
      endif !----------------end detailed------------------

 999  continue 

      if(.not.(psjetj.le.0..or.psjetj.ge.0.))then !NaN catch
      print*,psjetj,j,l,ak11,scale
      stop'NaN in psjetj'
      endif
      end

c-----------------------------------------------------------------------
      function psjetk(ipm,klasxx,tmx,q1,q2,scale,t,s,jx,lx,jdis)
c-----------------------------------------------------------------------
c
c klasxx > 0 : fast version, regoups individual xsections
c           jx,lx have to be understood as j1,l1 representing classes of xsections
c
c klasxx < 0 : detailed version, computes individual xsections  
c           jx,lx have to be understood as j2,l2 representing individual pair
c
c            (see subroutine initEoL for these definitions)
c
c klasxx = -9999 : goto 99
c           This is needed in case psjetk was not called for klas=klasmax 
c
c q1 - virtuality cutoff at current end of the ladder;
c q2 - virtuality cutoff at opposite end of the ladder;
c-----------------------------------------------------------------------
c           Integrand of psborn (no emissions)
c-----------------------------------------------------------------------
      implicit none
#include "sem.h"
      real psjetk,a
      real q1,q2,scale
      double precision  ffborn
      integer ipm,klas,klasx,klasxx,j,l,jdis,is,i,k,n1,n2,nfjet,jx,lx,kl
      double precision s,t,tmx,u,pt2
      real gg,gq,qq,qa,qqp,gc,gb,qc,qb,cc,cac,bab,bb,cb
      double precision ggX,gqX,qqX,qaX,qqpX,gcX,gbX,qcX,qbX,ccX,cacX
     .,bbX,babX,cbX
      double precision ggY,gqY,qqY,qaY,qqpY,gcY,gbY,qcY,qbY,ccY,cacY
     .,bbY,babY,cbY
      parameter (nfjet=5)
      double precision ff(-nfjet:nfjet,-nfjet:nfjet) 
      double precision ffk(klasmax,-nfjet:nfjet,-nfjet:nfjet) 
      integer j1,j2,l1,l2,jl,m,n,ishini,ish ,ifchx,igetIsh,ifmtx
      double precision weightx
      integer kflav1(-nfjet:nfjet),nof,noflav
      !integer ladderindex
      !integer ival
      !real val5,val6
      logical cha,bot
      real test10(10)
      common /ctest10/ test10 

      ish=igetIsh()
      call utpri('psjetk',ish,ishini,9)

      if(klasxx.eq.-9999)then
        klasx=-klasmax
        klas=klasmax
        goto 99 
      else
        klasx=klasxx
      endif

      if(klasx.lt.0)then !-----------detailed-----------
        j2=jx !initial parton
        l2=lx !    pair
        call eolpib( j2 , l2 , jl )
        call eolpi( jl , j1 , l1 )  
        j=j1 !initial parton
        l=l1 !  pair group
      else !----------------------end detailed-----------
        j=jx !initial parton
        l=lx !  pair group
      endif  

      psjetk=0

      if(iabs(j).eq.3.and.  noflav(q1).lt.4 ) goto 999 ! 1 charm
      if(iabs(j).eq.4.and.  noflav(q1).lt.5 ) goto 999 ! 1 bottom
      if(iabs(l).eq.3.and.  noflav(q2).lt.4 ) goto 999 ! 2 charm
      if(iabs(l).eq.4.and.  noflav(q2).lt.5 ) goto 999 ! 2 bottom

      klas=abs(klasx)
      nof=noflav(scale)

      if(jdis.ge.0.or.klasx.gt.0)then   !use factbgg only for detailed calculations
        is=1
      else
        is=-1
      endif

      !KW2006 cha, bot added
      cha=.false.
      bot=.false.
      if(nof.ge.4)cha=.true. 
      if(nof.ge.5)bot=.true.

      if(klasx.lt.0)then !-----------detailed-----------

        !flavor restrictions

        do i=-nfjet,nfjet
          kflav1(i)=1
          if(nof.lt.abs(i))then
            kflav1(i)=0
          endif
        enddo 

        !ff matrices

        do i=-nfjet,nfjet
        do k=-nfjet,nfjet
          if(i.eq.j2.and.k.eq.l2)then
            ff(i,k)=kflav1(i)*kflav1(k)  !KW2006 *kflav1 added
          else
            ff(i,k)=0
          endif 
          ffk(klas,i,k)=ff(i,k)
        enddo
        enddo

      endif !-----------end detailed-----------


      gg=0.
      gq=0.
      gc=0. 
      gb=0. 
      qq=0.
      qa=0.
      qqp=0.
      qc=0. 
      qb=0. 
      cc=0. 
      cac=0.
      bab=0.
      bb=0. 
      cb=0. 

      if(l.eq.0)then
        if(j.eq.0)then               !j,l = 0,0
          gg=1
        elseif(j.eq.1)then           !j,l = 1,0
          gq=1
        elseif(j.eq.3)then           !j,l = 3,0 
          if(cha)gc=1
        elseif(j.eq.4)then           !j,l = 4,0 
          if(bot)gb=1
        else
          print *,'j,l',j,l
          stop'ERROR 30082018'
        endif
      elseif(j.eq.0)then
        if(l.eq.1) then              !j,l = 0,1    (stands for 0,q)
          gq=1
        elseif(l.eq.3) then          !j,l = 0,3
          if(cha)gc=1
        elseif(l.eq.4) then          !j,l = 0,4
          if(bot)gb=1
        else
          print *,'j,l',j,l
           stop'ERROR 30082018'
        endif
      elseif(j.eq.1.and.l.eq.1)then  !j,l = 1,1  
        qq=1
      elseif(j.eq.2.and.l.eq.-2)then !j,l = 2,-2  
        qa=1
      elseif(j.eq.2.and.l.eq.1)then  !j,l = 2,1 
        qqp=1
      elseif(j.eq.1.and.l.eq.3)then  !j,l = 1,3 
        if(cha)qc=1
      elseif(j.eq.3.and.l.eq.1)then  !j,l = 3,1 
        if(cha)qc=1
      elseif(j.eq.1.and.l.eq.4)then  !j,l = 1,4 
        if(bot)qb=1
      elseif(j.eq.4.and.l.eq.1)then  !j,l = 4,1 
        if(bot)qb=1
      elseif(j.eq.3.and.l.eq.3)then  !j,l = 3,3 
        if(cha)cc=1
      elseif(j.eq.3.and.l.eq.-3)then !j,l = 3,-3 
        if(cha)cac=1
      elseif(j.eq.4.and.l.eq.4)then  !j,l = 4,4
        if(bot)bb=1
      elseif(j.eq.4.and.l.eq.-4)then !j,l = 4,-4
        if(bot)bab=1
      elseif(j.eq.4.and.l.eq.3)then  !j,l = 4,3
        if(bot)cb=1
      elseif(j.eq.3.and.l.eq.4)then  !j,l = 3,4
        if(bot)cb=1
      else
        print *,'j,l',j,l
        stop'ERROR 30082018'
      endif
      n1=1
      n2=3

      psjetk=ffborn(ipm,klas,tmx,s,t,gg,gq,qq,qa,qqp,n1,n2,3 ) +
     *       ffborn(ipm,klas,tmx,s,t,0.,gc,0.,0. ,qc ,n1,n2,4 ) +
     *       ffborn(ipm,klas,tmx,s,t,0.,gb,0.,0. ,qb ,n1,n2,5 ) +
     *       ffborn(ipm,klas,tmx,s,t,0.,0.,cc,cac,0. ,n1,n2,6 ) +
     *       ffborn(ipm,klas,tmx,s,t,0.,0.,bb,bab,cb ,n1,n2,7 )

      if(ish.ge.9)then
        call getMonitorFileIndex(ifmtx)
        write(ifmtx,'(a,2i3,3e12.3,3x,$)')'psjetk = ',j,l,psjetk,s,t !Sudakov not included
        !write(ifmtx,'(5e12.3)')gg,gq,qq,qa,qqp
        write(ifmtx,'(3x,9e12.3)')gc,qc,gb,qb,cc,cac,bb,bab,cb
      endif

      call getBornPtr2(klas,s,t , pt2) 
      call getBornU2(klas,s,t,u)

      psjetk = psjetk ! / abs(u-t) !* u * t / pt2 / (u+t) 

      call evTest1(gg,gq,gc,gb,qq,qa,qqp,qc,qb,cc,cac,bb,bab,cb)

      if(klasx.lt.0)then        !------------- detailed------------------
        !the values gg,gq etc do not depend on klas, so in principle they could 
        !be computed and tested once (if klas.eq.1 ...), and then stored for (a bit) more efficiency
        call setFFtypeZero(ggX,gqX,gcX,gbX,qqX,qaX,qqpX,qcX,qbX
     .             ,ccX,cacX,bbX,babX,cbX) !set all values to zero
        do m=-nfjet,nfjet ! in Born parton
        do n=-nfjet,nfjet ! in Born parton
          call setFFtypeZero(                       !set all values to zero
     .      ggY,gqY,gcY,gbY,qqY,qaY,qqpY,qcY,qbY,ccY,cacY,bbY,babY,cbY) 
          call getFFtype(ff,m,n,                    !given ff,m,n: get gg,gq etc 
     .      ggY,gqY,gcY,gbY,qqY,qaY,qqpY,qcY,qbY,ccY,cacY,bbY,babY,cbY) 
          call getFFtypeSngl(                          !trafo dble to sngl
     .      ggY,gqY,gcY,gbY,qqY,qaY,qqpY,qcY,qbY,ccY,cacY,bbY,babY,cbY,
     .      gg,gq,gc,gb,qq,qa,qqp,qc,qb,cc,cac,bb,bab,cb)
          !call setBornIn(m,n)    !store Born-in flavors, just for info, to be known in psabori
          !call setIshy(11)
          weightx =
     .       ffborn(ipm,klas,tmx,s,t,gg,gq,qq,qa ,qqp,n1,n2,is*3) +
     .       ffborn(ipm,klas,tmx,s,t,0.,gc,0.,0. ,qc ,n1,n2,is*4 ) +
     .       ffborn(ipm,klas,tmx,s,t,0.,gb,0.,0. ,qb ,n1,n2,is*5 ) +
     .       ffborn(ipm,klas,tmx,s,t,0.,0.,cc,cac,0. ,n1,n2,is*6 ) +
     .       ffborn(ipm,klas,tmx,s,t,0.,0.,bb,bab,cb ,n1,n2,is*7 )
          !call setIshy(0)
          !++++++++++++
          !call getLadderindex(ladderindex)
          !if((klas.eq.1.or.klas.eq.4).and.n.eq.0.and.m.eq.0)
          !.print*,'CHECKlad',ladderindex,klas,weightx
          !++++++++++++
          call setBornInWeight(klas,m,n,weightx)  !fill weight table
          call sumFFtype(  !add Y-value to X-value
     .      ggX,gqX,gcX,gbX,qqX,qaX,qqpX,qcX,qbX,ccX,cacX,bbX,babX,cbX,
     .      ggY,gqY,gcY,gbY,qqY,qaY,qqpY,qcY,qbY,ccY,cacY,bbY,babY,cbY)
        enddo
        enddo
        call setTmax(klas,tmx) !store tmx    
        call evTest2(1,j,l,klas
     .  ,ggX,gqX,gcX,gbX,qqX,qaX,qqpX,qcX,qbX,ccX,cacX,bbX,babX,cbX)
      endif !----------------end detailed------------------

  99  continue

      if(klasx.lt.0)then !------------- detailed------------------
        if(klas.eq.klasmax)then
          test10(1)=1
          call getBornInRdm(kl,m,n)      !generate kl,m,n
          !++++++++
          !call getLadderindex(ladderindex)
          !ival=0
          !if(ladderindex.eq.5.and.kl.eq.4)ival=5
          !if(ladderindex.lt.5.and.kl.eq.4)ival=6
          !if(ival.gt.0)call incVparam(ival,1.)
          !call getVparam(5,val5)
          !call getVparam(6,val6)
          !print*,'CHECKlad-----------------------------',ladderindex,kl
          !.,nint(val5),nint(val6)
          !+++++++++
          call getTmax(kl,tmx)             !get stored tmx for kl    
          do i=-nfjet,nfjet
          do k=-nfjet,nfjet
            ff(i,k)=ffk(kl,i,k)            !get ff for kl
          enddo
          enddo
          call setBornIn(m,n)             !store Born-in flavors  
          call getFFtype(ff,m,n,              !given ff,m,n: get gg,gq etc 
     .      ggY,gqY,gcY,gbY,qqY,qaY,qqpY,qcY,qbY,ccY,cacY,bbY,babY,cbY) 
          call getFFtypeSngl(                   !trafo dble to sngl
     .      ggY,gqY,gcY,gbY,qqY,qaY,qqpY,qcY,qbY,ccY,cacY,bbY,babY,cbY,
     .      gg,gq,gc,gb,qq,qa,qqp,qc,qb,cc,cac,bb,bab,cb)
          call setBornOutWeightsZero !initialize weight table
          call getCheckFileIndex(ifchx)
          if(ish.eq.-1007)rewind(ifchx)
          !compute Born-out weight table (using m,n via getBornIn)
          a=ffborn(ipm,kl,tmx,s,t,gg,gq,qq,qa ,qqp,n1,n2,is*3) +
     .      ffborn(ipm,kl,tmx,s,t,0.,gc,0.,0. ,qc ,n1,n2,is*4 ) +
     .      ffborn(ipm,kl,tmx,s,t,0.,gb,0.,0. ,qb ,n1,n2,is*5 ) +
     .      ffborn(ipm,kl,tmx,s,t,0.,0.,cc,cac,0. ,n1,n2,is*6 ) +
     .      ffborn(ipm,kl,tmx,s,t,0.,0.,bb,bab,cb ,n1,n2,is*7 )      !psjetk
          !table will be used via getBornOutRdm to generate the Born-out flavors
          !~~~
          ! call checkBornOutWeights(1,kl,m,n
          !.,gg,gq,gc,gb,qq,qa,qqp,qc,qb,cc,cac,bb,bab,cb)
          !~~~
        endif
      endif !----------------end detailed------------------

 999  continue

      call utprix('psjetk',ish,ishini,9)
      end
 
c------------------------------------------------------------------------
      double precision function ffborn(ipm,klas,tmax,s,t
     .                        ,ggi,gqi,qqi,qai,qqpi,n1,n2,nfi)
c------------------------------------------------------------------------
#include "sem.h"
      double precision psbori,s,t,tmax
c      real crash(2)
c      icrash=3

      nf=abs(nfi)
      is=sign(1,nfi)

      gg = ggi
      gq = gqi
      qq = qqi
      qa = qai
      qqp= qqpi
      
      call pdfparamget(1,facpdf4) 
      call pdfparamget(3,facpdf5) 
      if(ipm.lt.0)then
      call pdfparamget(2,facpdf4)  !sat pom
      call pdfparamget(4,facpdf5)  !sat pom
      endif
      if(nf.eq.4)then
        gq = gqi *facpdf4
        qqp= qqpi*facpdf4
      elseif(nf.eq.5)then
        gq = gqi *facpdf5
        qqp= qqpi*facpdf5
      elseif(nf.eq.6)then
        qq = qqi *facpdf4*facpdf4
        qa = qai *facpdf4*facpdf4
      elseif(nf.eq.7)then
        qq = qqi *facpdf5*facpdf5
        qa = qai *facpdf5*facpdf5
        qqp= qqpi*facpdf4*facpdf5
      endif
 
      ffborn=0.d0               

      if(gg.le.0..and.gq.le.0..and.qq.le.0..and.qa.le.0..and.qqp.le.0.)
     .return

      if(nf.eq.3)then !bg light quarks

       do n=n1,n2
        factqqb=1.
        if(n.eq.3)factqqb=2. !bg for qq~ -> gg because of outgoing identical particles
        factgg=1.
        if(n.eq.1)factgg=2. !bg for gg -> gg because of outgoing identical particles
        call setFactBorn(gg/factgg)
        ffborn=ffborn+(psbori(klas,s,t,0, 0,is*n)
     *                +psbori(klas,s,2*tmax-t,0, 0,is*n))   *gg/factgg 
        call setFactBorn(gq)
        ffborn=ffborn+(psbori(klas,s,t,0, 1,is*n)
     *                +psbori(klas,s,2*tmax-t,0, 1,is*n))   *gq 
        call setFactBorn(qq)
        ffborn=ffborn+(psbori(klas,s,t,1, 1,is*n)
     *                +psbori(klas,s,2*tmax-t,1, 1,is*n))/2 *qq 
        call setFactBorn(qa/factqqb)
        ffborn=ffborn +(psbori(klas,s,t,1,-1,is*n)
     *                 +psbori(klas,s,2*tmax-t,1,-1,is*n))  *qa/factqqb
        call setFactBorn(qqp)
        ffborn=ffborn +(psbori(klas,s,t,1, 2,is*n)
     *                 +psbori(klas,s,2*tmax-t,1, 2,is*n))  *qqp
       enddo

      elseif(nf.eq.4)then !charm

       do n=n1,n2
        call setFactBorn(gq)
        ffborn=ffborn  +
     *  (psbori(klas,s,t,4,0,is*n)+psbori(klas,s,2*tmax-t,4,0,is*n))*gq
        call setFactBorn(qqp)
        ffborn=ffborn  +
     *  (psbori(klas,s,t,4,1,is*n)+psbori(klas,s,2*tmax-t,4,1,is*n))*qqp
       enddo

      elseif(nf.eq.5)then       !bottom

       do n=n1,n2              
        call setFactBorn(gq)
        ffborn=ffborn  +     
     *  (psbori(klas,s,t,5,0,is*n)+psbori(klas,s,2*tmax-t,5,0,is*n))*gq
        call setFactBorn(qqp)
        ffborn=ffborn  +     
     *  (psbori(klas,s,t,5,1,is*n)+psbori(klas,s,2*tmax-t,5,1,is*n))*qqp
       enddo

      elseif(nf.eq.6)then       !charm-charm 

       do n=n1,n2              
       call setFactBorn(qq/2)
       ffborn=ffborn  +     
     *(psbori(klas,s,t,4,4,is*n)+psbori(klas,s,2*tmax-t,4,4,is*n))*qq/2 
       call setFactBorn(qa)
       ffborn=ffborn  +     
     *(psbori(klas,s,t,4,-4,is*n)+psbori(klas,s,2*tmax-t,4,-4,is*n))*qa 
       enddo

      elseif(nf.eq.7)then       !bottom-bottom and bottom-charm

       do n=n1,n2              
       call setFactBorn(qq/2)
       ffborn=ffborn  +      
     * (psbori(klas,s,t,5,5,is*n)+psbori(klas,s,2*tmax-t,5,5,is*n))*qq/2
       call setFactBorn(qa)
       ffborn=ffborn  +      
     * (psbori(klas,s,t,5,-5,is*n)+psbori(klas,s,2*tmax-t,5,-5,is*n))*qa 
       call setFactBorn(qqp)
       ffborn=ffborn  +      
     * (psbori(klas,s,t,4,5,is*n)+psbori(klas,s,2*tmax-t,4,5,is*n)) *qqp 
       enddo

        else
         stop'\n\n ERROR 201012 \n\n'
      endif

cc      call getMonitorFileIndex(ifmtx)
cc      if(ffborn.lt.0.)then !negative cross section
cc        write(ifmtx,*) 'neg CS :   ffborn = ',ffborn
cc        write(ifmtx,*) 's,t,nf : ',s,t,nf
cc        write(ifmtx,*) 'gg,gq,qq,qa,qqp:',gg,gq,qq,qa,qqp
cc        stop'ERROR 18122018c'
cc      endif

c      if(.not.(ffborn.le.0..or.ffborn.ge.0.))then !NaN catch
c      print*,ffborn,nf,gg,gq,qq,qa,qqp
c      crash(icrash)=1. !to force crash
c      stop'NaN in ffborn'
c      endif

      end


c#########################################################################################
c#########################################################################################
c
c                      set - get functions     and evTest....
c
c#########################################################################################
c#########################################################################################

c
c    get overview with " grep -e 'ine\ set\|ine\ get' KW/sem.f "
c

c------------------------------------------------------------------------
      subroutine setTmax(klas,tmax)                    !given klas: store tmax
c------------------------------------------------------------------------
      implicit none
#include "sem.h"
      integer klas
      double precision tmax,tmaxar
      common /ctmaxklas/ tmaxar(klasmax)
      tmaxar(klas)=tmax
      end

c------------------------------------------------------------------------
      subroutine getTmax(klas,tmax)                       !get stored tmax
c------------------------------------------------------------------------
      implicit none
#include "sem.h"
      integer klas
      double precision tmax,tmaxar
      common /ctmaxklas/ tmaxar(klasmax)
      tmax=tmaxar(klas)
      end

c------------------------------------------------------------------------
      subroutine setFactBorn(fac)                       !store bornfactor f*f
c------------------------------------------------------------------------
      common /cbornfactor/bornfactor
      bornfactor=fac
      end

c------------------------------------------------------------------------
      subroutine getFactBorn(fac)                    !get stored bornfactor f*f
c------------------------------------------------------------------------
      common /cbornfactor/bornfactor
      fac=bornfactor
      end

c-----------------------------------------------------------------------
      subroutine setFFtypeZero(ggX,gqX,gcX,gbX          !set all values to zero
     .                 ,qqX,qaX,qqpX,qcX,qbX,ccX,cacX,bbX,babX,cbX) 
c-----------------------------------------------------------------------
      implicit none
      double precision ggX,gqX,gcX,gbX,qqX,qaX,qqpX,qcX,qbX
     .                  ,ccX,cacX,bbX,babX,cbX
      ggX=0.d0
      gqX=0.d0
      gcX=0.d0 
      gbX=0.d0 
      qqX=0.d0
      qaX=0.d0
      qqpX=0.d0
      qcX=0.d0 
      qbX=0.d0 
      ccX=0.d0 
      cacX=0.d0
      bbX=0.d0 
      babX=0.d0
      cbX=0.d0 
      end 
c-----------------------------------------------------------------------
      subroutine getFFtype(ff,m,n,ggY,gqY,gcY         !given m,n: get gg,gq etc
     .  ,gbY,qqY,qaY,qqpY,qcY,qbY,ccY,cacY,bbY,babY,cbY)
c-----------------------------------------------------------------------
      implicit none
      integer m,n,mab,nab
      double precision ggY,gqY,gcY,gbY,qqY,qaY,qqpY,qcY,qbY
     .                  ,ccY,cacY,bbY,babY,cbY
      integer nfjet
      parameter (nfjet=5)
      double precision ff(-nfjet:nfjet,-nfjet:nfjet) 
      mab=abs(m)
      nab=abs(n)  
      if(m.eq.0.and.n.eq.0)then
        ggY = ff(m,n) 
      elseif(m.eq.0.and.nab.le.3.or.n.eq.0.and.mab.le.3)then
        gqY = ff(m,n)
      elseif(m.eq.0.and.nab.eq.4.or.n.eq.0.and.mab.eq.4)then
        gcY = ff(m,n)
      elseif(m.eq.0.and.nab.eq.5.or.n.eq.0.and.mab.eq.5)then
        gbY = ff(m,n)
      elseif(m.eq.n.and.nab.le.3)then
        qqY = ff(m,n)
      elseif(m.eq.-n.and.nab.le.3)then
        qaY = ff(m,n)
      elseif(nab.le.3.and.mab.le.3)then
        qqpY= ff(m,n)
      elseif(mab.le.3.and.nab.eq.4.or.nab.le.3.and.mab.eq.4)then
        qcY = ff(m,n)
      elseif(mab.le.3.and.nab.eq.5.or.nab.le.3.and.mab.eq.5)then
        qbY = ff(m,n)
      elseif(m.eq.4.and.n.eq.4.or.m.eq.-4.and.n.eq.-4)then
        ccY = ff(m,n)
      elseif(m.eq.4.and.n.eq.-4.or.m.eq.-4.and.n.eq.4)then
        cacY= ff(m,n)
      elseif(m.eq.5.and.n.eq.5.or.m.eq.-5.and.n.eq.-5)then
        bbY = ff(m,n)
      elseif(m.eq.5.and.n.eq.-5.or.m.eq.-5.and.n.eq.5)then
        babY= ff(m,n)
      elseif(mab.le.4.and.nab.eq.5.or.nab.le.4.and.mab.eq.5)then
        cbY = ff(m,n)
      endif
      end
c-----------------------------------------------------------------------
      subroutine getFFtypeSngl(                       !trafo dble to sngl
     .    ggY,gqY,gcY,gbY,qqY,qaY,qqpY,qcY,qbY,ccY,cacY,bbY,babY,cbY,
     .    gg,gq,gc,gb,qq,qa,qqp,qc,qb,cc,cac,bb,bab,cb)
c-----------------------------------------------------------------------
      implicit none
      double precision ggY,gqY,gcY,gbY,qqY,qaY,qqpY,qcY,qbY
     .                  ,ccY,cacY,bbY,babY,cbY
      real gg,gq,gc,gb,qq,qa,qqp,qc,qb,cc,cac,bb,bab,cb
      gg = ggY
      gq = gqY
      gc = gcY
      gb = gbY
      qq = qqY
      qa = qaY
      qqp= qqpY
      qc = qcY
      qb = qbY
      cc = ccY
      cac= cacY
      bb = bbY
      bab= babY
      cb = cbY
      end  
c-----------------------------------------------------------------------
      subroutine sumFFtype(                          !add Y-value to X-value
     .   ggX,gqX,gcX,gbX,qqX,qaX,qqpX,qcX,qbX,ccX,cacX,bbX,babX,cbX,
     .   ggY,gqY,gcY,gbY,qqY,qaY,qqpY,qcY,qbY,ccY,cacY,bbY,babY,cbY)
c-----------------------------------------------------------------------
      implicit none
      double precision ggX,gqX,gcX,gbX,qqX,qaX,qqpX,qcX,qbX
     .                  ,ccX,cacX,bbX,babX,cbX
      double precision ggY,gqY,gcY,gbY,qqY,qaY,qqpY,qcY,qbY
     .                  ,ccY,cacY,bbY,babY,cbY
      ggX =ggX +ggY
      gqX =gqX +gqY
      gcX =gcX +gcY
      gbX =gbX +gbY
      qqX =qqX +qqY
      qaX =qaX +qaY
      qqpX=qqpX+qqpY
      qcX =qcX +qcY
      qbX =qbX +qbY
      ccX =ccX +ccY
      cacX=cacX+cacY
      bbX =bbX +bbY
      babX=babX+babY
      cbX =cbX +cbY
      end
c-----------------------------------------------------------------------
      subroutine setBornIn(i,j)                        !store Born-in flavors 
c-----------------------------------------------------------------------   
      implicit none
      integer mInB,nInB,i,j
      common /cinborn/mInB,nInB
      mInB = i
      nInB = j
      end
c-----------------------------------------------------------------------
      subroutine getBornIn(i,j)                    !get stored Born-in flavors 
c-----------------------------------------------------------------------
      implicit none
      integer mInB,nInB,i,j
      common /cinborn/mInB,nInB
      i = mInB
      j = nInB
      end
c-----------------------------------------------------------------------
      subroutine setBornInWeightsZero
c-----------------------------------------------------------------------
      implicit none
#include "sem.h"
      integer nfjet,kl,minB,ninB
      parameter (nfjet=5)
      double precision weighta
      common /cweighta/ weighta(klasmax,-nfjet:nfjet,-nfjet:nfjet)
      do kl=1,klasmax
      do minB=-nfjet,nfjet
      do ninB=-nfjet,nfjet
        weighta(kl,minB,ninB)=0
      enddo
      enddo 
      enddo
      end
c-----------------------------------------------------------------------
      subroutine setBornInWeight(kl,m,n,w)        !given kl,m,n: set weight w
c-----------------------------------------------------------------------
      ! kl = klas, m,n = Born-in flavors
      implicit none
#include "sem.h"
      integer nfjet,kl,m,n
      parameter (nfjet=5)
      double precision weighta,w
      common /cweighta/ weighta(klasmax,-nfjet:nfjet,-nfjet:nfjet)
        weighta(kl,m,n)=w
      end
c-----------------------------------------------------------------------
      subroutine getBornInRdm(kl,m,n)          !get kl,m,n according to weighta
c-----------------------------------------------------------------------
      implicit none
#include "sem.h"
      integer nfjet !,nn
      real test10(10)
      common /ctest10/ test10 
      parameter (nfjet=5)
      double precision weighta
      common /cweighta/ weighta(klasmax,-nfjet:nfjet,-nfjet:nfjet)
      integer kl,m,n
      double precision checkWeightSums(2,-klasmax:klasmax)
      common /ccheckWeightSums/ checkWeightSums
      double precision sum1,sum,r
      real rangen
      sum=0.d0
      do kl=1,klasmax
      do m=-nfjet,nfjet
      do n=-nfjet,nfjet
        sum=sum+weighta( kl , m, n ) 
      enddo
      enddo
      enddo
      r=dble(rangen())*sum
      sum1=sum
      sum=0.d0
      kl=1
      m=-nfjet
      n=-nfjet - 1
      do while (sum.lt.r)
        n=n+1
        if(n.gt.nfjet)then
          n=-nfjet
          m=m+1
        endif
        if(m.gt.nfjet)then
          m=-nfjet
          kl=kl+1
        endif
        sum=sum+weighta( kl , m, n ) 
      enddo    
      if(abs(m).gt.nfjet.or.abs(n).gt.nfjet)then
       print*,'getBornInRdm: r sum test10 = ',r,sum,test10
       !if(test10(1).eq.23.)then !ladd 
       !print*,'checkWeightSums =',(checkWeightSums(1,nn),nn=1,klasmax)
       !print*,'checkWeightSums =',(checkWeightSums(1,-nn),nn=1,klasmax)
       !print*,'checkWeightSums =',(checkWeightSums(2,nn),nn=1,klasmax)
       !print*,'checkWeightSums =',(checkWeightSums(2,-nn),nn=1,klasmax)
       !endif
       stop'in getBornInRdm; flavor >5'
      endif
      end
c-----------------------------------------------------------------------
          subroutine checkBornOutWeights(ladd,kl,m,n
     .    ,gg,gq,gc,gb,qq,qa,qqp,qc,qb,cc,cac,bb,bab,cb)
c-----------------------------------------------------------------------
      implicit none
#include "sem.h"
      integer ladd,kl,m,n,mm,nn,nfjet
      real gg,gq,gc,gb,qq,qa,qqp,qc,qb,cc,cac,bb,bab,cb
      parameter (nfjet=5)
      double precision weightb(-nfjet:nfjet,-nfjet:nfjet)
      common /cweightb/ weightb
      double precision weighta
      common /cweighta/ weighta(klasmax,-nfjet:nfjet,-nfjet:nfjet)
      double precision sum
      sum=0.d0
      do m=-nfjet,nfjet
      do n=-nfjet,nfjet
        sum=sum+weightb( m, n ) 
      enddo
      enddo
      if(sum.eq.0d0)then
        print*,'ERROR in checkBornOutWeights'
        print*,'ladd,kl,m,n =',ladd,kl,m,n
        print*,'gg,gq,gc,gb,qq,qa,qqp =',gg,gq,gc,gb,qq,qa,qqp
        print*,'qc,qb,cc,cac,bb,bab,cb =',qc,qb,cc,cac,bb,bab,cb
        do mm=-nfjet,nfjet
          print*,'weighta (kl,',mm,',...) ='
          print*,(weighta( kl , mm, nn ), nn= -nfjet,nfjet)
        enddo
        stop'checkBornOutWeights: zero weights'
      endif
      end
c-----------------------------------------------------------------------
      subroutine setBornOutWeightsZero
c-----------------------------------------------------------------------
      implicit none
      integer nfjet,moutB,noutB
      parameter (nfjet=5)
      double precision weightb(-nfjet:nfjet,-nfjet:nfjet)
      common /cweightb/ weightb
      do moutB=-nfjet,nfjet
      do noutB=-nfjet,nfjet
        weightb(moutB,noutB)=0
      enddo
      enddo 
      end
c-----------------------------------------------------------------------
      subroutine getBornOutRdm(m,n)              !get m,n according to weightb
c-----------------------------------------------------------------------
      implicit none
      double precision sum,r
      integer nfjet,m,n
      integer mInB,nInB
      parameter (nfjet=5)
      double precision weightb(-nfjet:nfjet,-nfjet:nfjet)
      real rangen
      common /cweightb/ weightb
      sum=0.d0
      do m=-nfjet,nfjet
      do n=-nfjet,nfjet
        sum=sum+weightb( m, n )
      enddo
      enddo
      r=dble(rangen())*sum
      sum=0.d0
      m=-nfjet
      n=-nfjet - 1
      do while (sum.lt.r)
        n=n+1
        if(n.gt.nfjet)then
          n=-nfjet
          m=m+1
        endif
        sum=sum+weightb( m, n ) 
      enddo          
      if(abs(m).gt.nfjet.or.abs(n).gt.nfjet)then
        call getBornIn(mInB,nInB) !get stored Born-in flavors 
        print*,'getBornOutRdm: mInB,nInB:',mInB,nInB
        print*,'getBornOutRdm: r sum = ',r,sum
        stop'in getBornOutRdm; flavor >5'
      endif
      end
c-----------------------------------------------------------------------
      subroutine evTest1(gg,gq,gc,gb,qq,qa,qqp,qc,qb,cc,cac,bb,bab,cb)
c-----------------------------------------------------------------------
      implicit none
#include "aaa.h"
      real gg,gq,gc,gb,qq,qa,qqp,qc,qb,cc,cac,bb,bab,cb 
      double precision ggX,gqX,gcX,gbX,qqX,qaX,qqpX,qcX,qbX,ccX,cacX
     .,bbX,babX,cbX 
      common/cevtest/
     .ggX,gqX,gcX,gbX,qqX,qaX,qqpX,qcX,qbX,ccX,cacX,bbX,babX,cbX 
      ggX =gg
      gqX =gq
      gcX =gc
      gbX =gb
      qqX =qq
      qaX =qa
      qqpX=qqp
      qcX =qc
      qbX =qb
      ccX =cc
      cacX=cac
      bbX =bb
      babX=bab
      cbX =cb
      end
c-----------------------------------------------------------------------
      subroutine evTest2(ii,j,l,klas
     .,gg,gq,gc,gb,qq,qa,qqp,qc,qb,cc,cac,bb,bab,cb)
c-----------------------------------------------------------------------
      implicit none
#include "aaa.h"
      double precision gg,gq,gc,gb,qq,qa,qqp,qc,qb,cc,cac,bb,bab,cb 
      double precision ggX,gqX,gcX,gbX,qqX,qaX,qqpX,qcX,qbX,ccX,cacX
     .,bbX,babX,cbX
      common/cevtest/
     .ggX,gqX,gcX,gbX,qqX,qaX,qqpX,qcX,qbX,ccX,cacX,bbX,babX,cbX 
      integer klas,ii,j,l
      if(ggX.ne.0) call evTest3(ii,j,l,klas,'gg ',ggX,gg  )
      if(gqX.ne.0) call evTest3(ii,j,l,klas,'gq ',gqX,gq  )
      if(gcX.ne.0) call evTest3(ii,j,l,klas,'gc ',gcX,gc  )
      if(gbX.ne.0) call evTest3(ii,j,l,klas,'gb ',gbX,gb  )
      if(qqX.ne.0) call evTest3(ii,j,l,klas,'qq ',qqX,qq  )
      if(qaX.ne.0) call evTest3(ii,j,l,klas,'qa ',qaX,qa  )
      if(qqpX.ne.0)call evTest3(ii,j,l,klas,'qqp',qqpX,qqp)
      if(qcX.ne.0) call evTest3(ii,j,l,klas,'qc ',qcX,qc  )
      if(qbX.ne.0) call evTest3(ii,j,l,klas,'qb ',qbX,qb  )
      if(ccX.ne.0) call evTest3(ii,j,l,klas,'cc ',ccX,cc  )
      if(cacX.ne.0)call evTest3(ii,j,l,klas,'cac',cacX,cac)
      if(bbX.ne.0) call evTest3(ii,j,l,klas,'bb ',bbX,bb  )
      if(babX.ne.0)call evTest3(ii,j,l,klas,'bab',babX,bab)
      if(cbX.ne.0) call evTest3(ii,j,l,klas,'cb ',cbX,cb  )
      if(gg.ne.0 .and.ggX.eq.0) call evTest4(ii,j,l,klas,'gg ',ggX,gg  )
      if(gq.ne.0 .and.gqX.eq.0) call evTest4(ii,j,l,klas,'gq ',gqX,gq  )
      if(gc.ne.0 .and.gcX.eq.0) call evTest4(ii,j,l,klas,'gc ',gcX,gc  )
      if(gb.ne.0 .and.gbX.eq.0) call evTest4(ii,j,l,klas,'gb ',gbX,gb  )
      if(qq.ne.0 .and.qqX.eq.0) call evTest4(ii,j,l,klas,'qq ',qqX,qq  )
      if(qa.ne.0 .and.qaX.eq.0) call evTest4(ii,j,l,klas,'qa ',qaX,qa  )
      if(qqp.ne.0.and.qqpX.eq.0)call evTest4(ii,j,l,klas,'qqp',qqpX,qqp)
      if(qc.ne.0 .and.qcX.eq.0) call evTest4(ii,j,l,klas,'qc ',qcX,qc  )
      if(qb.ne.0 .and.qbX.eq.0) call evTest4(ii,j,l,klas,'qb ',qbX,qb  )
      if(cc.ne.0 .and.ccX.eq.0) call evTest4(ii,j,l,klas,'cc ',ccX,cc  )
      if(cac.ne.0.and.cacX.eq.0)call evTest4(ii,j,l,klas,'cac',cacX,cac)
      if(bb.ne.0 .and.bbX.eq.0) call evTest4(ii,j,l,klas,'bb ',bbX,bb  )
      if(bab.ne.0.and.babX.eq.0)call evTest4(ii,j,l,klas,'bab',babX,bab)
      if(cb.ne.0 .and.cbX.eq.0) call evTest4(ii,j,l,klas,'cb ',cbX,cb  )
      end
      subroutine evTest3(ii,j,l,klas,text,val1,val2)
#include "aaa.h"
      double precision val1,val2
      character text*3
      if(abs(val2/val1-1.0).gt.0.001)then
        write(ifmt,'(a,4i3,2x,a,f8.3,2f14.9)')'evTest Values different'
     .  ,ii,j,l,klas,text,val2/val1,val1,val2
      endif
      end
      subroutine evTest4(ii,j,l,klas,text,val1,val2)
#include "aaa.h"
      double precision val1,val2
      character text*3
        write(ifmt,'(a,4i3,2x,a,2f8.3)')'evTest values zero nonzero'
     .  ,ii,j,l,klas,text,val1,val2
      end


c########################################################################
c########################################################################
c########################################################################
c
c                    Evolution  psev...
c
c        completely redone in 3317  (old version see 3314)
c
c########################################################################
c########################################################################
c########################################################################

c------------------------------------------------------------------------
      subroutine getGaussLegendre7(x1,a1)
c------------------------------------------------------------------------
      implicit none
      double precision x1(7),a1(7),xx1(7),aa1(7)
      integer i
      data xx1/.9862838d0,.9284349d0,.8272013d0,.6872929d0,.5152486d0,
     *.3191124d0,.1080549d0/
      data aa1/.03511946d0,.08015809d0,.1215186d0,.1572032d0,
     *.1855384d0,.2051985d0,.2152639d0/
      do i=1,7
        x1(i)=xx1(i)
        a1(i)=aa1(i)
      enddo
      end
 
c-----------------------------------------------------------------------
      double precision function psevz2(qi,qj,xx,m,k)
c-----------------------------------------------------------------------
c Evolution function for precisely one emission, see PR B32
c-----------------------------------------------------------------------
c Approximate formula: We use for gg and qq :
c   int dQ2/Q2 1/4.5*log(Q2/alm) * Sud * Sud * AP 
c     ~  2/9*log( log(qj/alm)/log(qi/alm) ) * Sud(qi,qj) * AP
c-----------------------------------------------------------------------
c The function is also used as factor for tabulation.
c-----------------------------------------------------------------------
      implicit double precision (a-h,o-z)
      real val,psfap
      psevz2=0d0
      if(xx.gt.1.d0-epscutSud2(qj))return
      call getQcdlam(val)
      alm=val
      if(m.eq.1.and.k.eq.1.or.m.ne.1.and.k.ne.1)then         !11,22,32
       psevz2=1.d0/4.5d0/psudi2(qi,m-1)*psudi2(qj,m-1)
     . *dlog(dlog(qj/alm)/dlog(qi/alm))    
      else                                                   !12,21
       psevz2=.3d0/(dlog(epscutSud2(qj))+.75d0)
     . *( psudi2(qj,0)/psudi2(qi,0)-psudi2(qj,1)/psudi2(qi,1) )     
      endif
      psevz2=psevz2*dble(psfap(sngl(qj),xx,m-1,k-1))
      end

c-----------------------------------------------------------------------
      double precision function psev2(qin,qfi,xx,j,l) !args double
c-----------------------------------------------------------------------
c Updates evolution function E(qi,qf,xx) (needs to be iterated)
c j,l = 1 : gluon;   j,l = 2,3 : quark
c-----------------------------------------------------------------------
      implicit double precision (a-h,o-z)
      real val,psfap
      double precision x1(7),a1(7)
      call getGaussLegendre7(x1,a1)
      call getQcdlam(val)
      alm=val
      q1=qin
      qm=qfi
      qj=qfi

      psev2=0.d0
      zmax=1.d0-epscutSud2(dble(qfi))
      zmin=xx/zmax
      if(zmin.ge.zmax)return

      do i1=1,7
      do m1=1,2
       qi=q1*(qj/q1)**(.5d0+x1(i1)*(m1-1.5d0))

       fz1=0.d0
       fz2=0.d0
       fz3=0.d0
       zmin1=max(.2d0,zmin)
       zmax1=min(.2d0,zmax)
       zmax1=min(5.d0*xx,zmax1)
       zmax2=min(zmin1,zmax)
       zmin2=max(zmax1,zmin)

       if(zmax1.gt.zmin)then
        do i=1,7
        do m=1,2
         z=xx+(zmin-xx)*((zmax1-xx)/(zmin-xx))**(.5d0+(m-1.5d0)*x1(i))
         do k=1,2
          if(j.ne.3.or.k.ne.1)then
           fz1=fz1+a1(i)*psevi2(q1,qi,xx/z,j,k)
     .              *dble(psfap(sngl(qi),z,k-1,l-1))*(1.d0-xx/z)
          endif
         enddo
        enddo
        enddo
        fz1=fz1*dlog((zmax1-xx)/(zmin-xx))
       endif
       if(zmin1.lt.zmax)then
        do i=1,7
        do m=1,2
         z=1.d0-(1.d0-zmax)*((1.d0-zmin1)/(1.d0-zmax))
     *   **(.5d0+x1(i)*(m-1.5d0))
         do k=1,2
          if(j.ne.3.or.k.ne.1)then
           fz2=fz2+a1(i)*psevi2(q1,qi,xx/z,j,k)
     *     *dble(psfap(sngl(qi),z,k-1,l-1))
     *     *(1.d0/z-1.d0)
          endif
         enddo
        enddo
        enddo
        fz2=fz2*dlog((1.d0-zmin1)/(1.d0-zmax))
       endif
       if(zmax2.gt.zmin2)then
        do i=1,7
        do m=1,2
         z=zmin2*(zmax2/zmin2)**(.5d0+x1(i)*(m-1.5d0))
         do k=1,2
          if(j.ne.3.or.k.ne.1)then
           fz3=fz3+a1(i)*psevi2(q1,qi,xx/z,j,k)
     *      *dble(psfap(sngl(qi),z,k-1,l-1))
          endif
         enddo
        enddo
        enddo
        fz3=fz3*dlog(zmax2/zmin2)
       endif
       psev2=psev2+a1(i1)*(fz1+fz2+fz3)/psudi2(qi,l-1)*pssalf2(qi/alm)
      enddo
      enddo

      psev2=psev2*dlog(qj/q1)/4.d0*psudi2(qj,l-1)

      end
 
c-----------------------------------------------------------------------
      double precision function psevev2(qin,qfm,qfi,xx,j,l) !args double
c-----------------------------------------------------------------------
c Computes \int dz/z E(qi,qe,xx/z)*E(qe,qf,z)
c j,l = 1 : gluon;   j,l = 2,3 : quark
c-----------------------------------------------------------------------
      implicit double precision (a-h,o-z)
c      real val
      double precision x1(7),a1(7)
      call getGaussLegendre7(x1,a1)
c      call getQcdlam(val)
c      alm=val
      q1=qin
      qm=qfm
      qj=qfi

      psevev2=0.d0
      zmax=1.d0-epscutSud2(dble(qfi))
      zmin=xx/zmax
      if(zmin.ge.zmax)return

      fz1=0.d0
      fz2=0.d0
      fz3=0.d0
      zmin1=max(.2d0,zmin)
      zmax1=min(.2d0,zmax)
      zmax1=min(5.d0*xx,zmax1)
      zmax2=min(zmin1,zmax)
      zmin2=max(zmax1,zmin)

      if(zmax1.gt.zmin)then
       do i=1,7
       do m=1,2
        z=xx+(zmin-xx)*((zmax1-xx)/(zmin-xx))**(.5d0+(m-1.5d0)*x1(i))
        do k=1,2
         if(j.ne.3)then
          fz1=fz1+a1(i)*psevi2(q1,qm,xx/z,j,k)*psevi2(qm,qj,z,k,l)
     *    *(1.d0-xx/z)
         elseif(k.ne.1)then
          fz1=fz1+a1(i)*psevi2(q1,qm,xx/z,3,2)*psevi2(qm,qj,z,3,2)
     *    *(1.d0-xx/z)
         endif
        enddo
       enddo
       enddo
       fz1=fz1*dlog((zmax1-xx)/(zmin-xx))
      endif
      if(zmin1.lt.zmax)then
       do i=1,7
       do m=1,2
        z=1.d0-(1.d0-zmax)*((1.d0-zmin1)/(1.d0-zmax))
     *  **(.5d0+x1(i)*(m-1.5d0))
        do k=1,2
         if(j.ne.3)then
          fz2=fz2+a1(i)*psevi2(q1,qm,xx/z,j,k)*psevi2(qm,qj,z,k,l)
     *    *(1.d0/z-1.d0)
         elseif(k.ne.1)then
          fz2=fz2+a1(i)*psevi2(q1,qm,xx/z,3,2)*psevi2(qm,qj,z,3,2)
     *    *(1.d0/z-1.d0)
         endif
        enddo
       enddo
       enddo
       fz2=fz2*dlog((1.d0-zmin1)/(1.d0-zmax))
      endif
      if(zmax2.gt.zmin2)then
       do i=1,7
       do m=1,2
        z=zmin2*(zmax2/zmin2)**(.5d0+x1(i)*(m-1.5d0))
        do k=1,2
         if(j.ne.3)then
          fz2=fz2+a1(i)*psevi2(q1,qm,xx/z,j,k)*psevi2(qm,qj,z,k,l)
         elseif(k.ne.1)then
          fz2=fz2+a1(i)*psevi2(q1,qm,xx/z,3,2)*psevi2(qm,qj,z,3,2)
         endif
        enddo
       enddo
       enddo
       fz3=fz3*dlog(zmax2/zmin2)
      endif
      psevev2=(fz1+fz2+fz3)/2.d0

      end

c-----------------------------------------------------------------------
      double precision function psevi2(q1,qq,xx,m,l) !args double
c-----------------------------------------------------------------------
c Interpolates evolution function
c m,l = 1 : gluon;   m,l = 2,3 : quark
c-----------------------------------------------------------------------
c       m l: 1 1 ... gluon -> gluon
c            2 1 ... quark -> gluon
c            1 2 ... gluon -> quark
c            3 2 ... quark -> quark non singlet
c            2 2 ... quark -> quark all
c                             singlet = all - non singlet
c-----------------------------------------------------------------------
      implicit double precision (a-h,o-z)
      real val
c      double precision x1(7),a1(7)
      dimension wi(3),wj(3),wk(3)
      common /evkarr/ evk(40,40,100,3,2)

c      call getGaussLegendre7(x1,a1)
c      call getQcdlam(val)
c      alm=val
      call getEpmax(val)
      spmax=val

      psevi2=0.d0
      fpsevz2=psevz2(q1,qq,xx,m,l)
      if(fpsevz2.eq.0.d0)return
      if(q1.ge..9999d0*spmax)goto 1

      if(xx.le..1d0)then
       yx=37.d0-dlog(.1d0/xx)/dlog(.1d0*spmax)*36.d0
       k=max(1,int(yx))
       k=min(k,35)
      elseif(xx.le..9d0)then
       yx=(xx-.1d0)*40.d0+37.d0
       k=max(37,int(yx))
       k=min(k,67)
      else
       yx=dlog(10.d0*(1.d0-xx))
     .       / log( 10.d0 * epscutSud2(dble(qq)) ) * 31.d0 + 69.d0
c       if(yx.gt.100d0)then
c         write(*,*)'Warning x~1 in psevi2 (out of range)',xx,yx
c         yx=100d0
c       endif
       k=max(69,int(yx))
       k=min(k,98)
      endif
      wk(2)=yx-k
      wk(3)=wk(2)*(wk(2)-1.d0)*.5d0
      wk(1)=1.d0-wk(2)+wk(3)
      wk(2)=wk(2)-2.d0*wk(3)
      !if(1.d0-xx.lt.epscutSud2(dble(qq)))print*,xx,yx,k

      qli=log(q1)/dlog(spmax)*39.d0+1.d0           ! q1 -> i
      qlj=log(qq/q1)/dlog(spmax/q1)*39.d0+1.d0     ! qq -> j
      i=max(1,int(1.0001d0*qli))
      i=min(i,38)
      wi(2)=qli-i
      wi(3)=wi(2)*(wi(2)-1.d0)*.5d0
      wi(1)=1.d0-wi(2)+wi(3)
      wi(2)=wi(2)-2.d0*wi(3)

      j=max(1,int(1.0001d0*qlj))
      j=min(j,38)
      wj(2)=qlj-j
      wj(3)=wj(2)*(wj(2)-1.d0)*.5d0
      wj(1)=1.d0-wj(2)+wj(3)
      wj(2)=wj(2)-2.d0*wj(3)

      if((wk(1).lt.-100 .or. wk(2).lt.-100 .or. wk(3).lt.-100)
     .  .and.evk(i,j,k,m,l).gt.1d-10)then
        call getMonitorFileIndex(ifmt)
        write(ifmt,*)'WARNING psevi2: x far outside range'
     .   ,xx,yx,(evk(i,j,k2,m,l),k2=k-2,k)
      endif

      do k1=1,3
       k2=k+k1-1
      do j1=1,3
      do i1=1,3
       psevi2=psevi2+evk(i+i1-1,j+j1-1,k2,m,l)*wi(i1)*wj(j1)*wk(k1)
       !if(xx.gt.0.99d0)then
       ! print*,'PSEVI',i+i1-1,j+j1-1,k2,evk(i+i1-1,j+j1-1,k2,m,l),psevi2
       !endif
      enddo
      enddo
      enddo
      tmp=psevi2

  1   psevi2 = exp(psevi2) * fpsevz2

      !~~~~~~~~~~~~~~~~~~~~~~~~
      ! if(igetIsh().ge.9)then
      ! if(xx.gt.0.99d0)write(*,'(a,5e12.3)')'PSEVI TEST',xx,tmp,q1,qq
      !.,psevz2(q1,qq,xx,m,l)
      ! do i=1,19
      !   x=0.980+0.001*i
      !   print*,m,l,x,psevz2(q1,qq,x,m,l)
      ! enddo
      ! endif
      !~~~~~~~~~~~~~~~~~~~~~~~~

      a=psevi2     ! NaN catch
      if(.not.(a.le.0..or.a.ge.0.))then
        print *,q1,qq,xx,m,l,a
        stop'ERROR 24022020'
      endif

      end

c-----------------------------------------------------------------------
      function psevi(q1,qq,xx,m,l) !q1,qq sngl
c-----------------------------------------------------------------------
c interpolates evolution function
c m,l = 1 : gluon;   m,l = 2,3 : quark
c-----------------------------------------------------------------------
      double precision xx,psevi2
      psevi=psevi2(dble(q1),dble(qq),xx,m,l)
      end

c-----------------------------------------------------------------------
      function psevj(q1,qq,xx,i,k)  !i,k: -5,...,5 (single partons)
c-----------------------------------------------------------------------
      double precision xx
      real ak1(0:1,0:1)
      nof=noflav(qq)
      kflav1=1
      if(nof.lt.abs(k))then
        kflav1=0
      endif
      ak1(0,0)=psevi(q1,qq,xx,1,1)                  
      ak1(0,1)=psevi(q1,qq,xx,1,2)/noflav(qq)/2. 
      ak1NS=psevi(q1,qq,xx,3,2)   !nonsinglet 
      ak1S=(psevi(q1,qq,xx,2,2)-ak1NS)/noflav(qq)/2. !singlet 
      ak1(1,0)=psevi(q1,qq,xx,2,1)  
      if(i.eq.0)then !g->g,q
        psevj=ak1(0,min(1,abs(k)))*kflav1
      elseif(k.eq.0)then !q->g
        psevj=ak1(min(1,abs(i)),0)
      elseif(i.eq.k)then !q->q
        psevj=(ak1NS+ak1S)*kflav1
      else !q->q'
        psevj=ak1S*kflav1
      endif
      end



c########################################################################
c########################################################################
c########################################################################
c
c                    Utilities
c
c########################################################################
c########################################################################
c########################################################################



c------------------------------------------------------------------------
      subroutine pscs(c,s)
c-----------------------------------------------------------------------
c pscs - cos (c) and sin (s) generation for uniformly distributed angle
c-----------------------------------------------------------------------
1     s1=2.*rangen()-1.
      s2=2.*rangen()-1.
      s3=s1*s1+s2*s2
      if(s3.gt.1.)goto 1
      s3=sqrt(s3)
      c=s1/s3
      s=s2/s3
      return
      end

c------------------------------------------------------------------------
      subroutine psdefrot(ep,s0x,c0x,s0,c0)
c-----------------------------------------------------------------------
c psdefrot - determination of cos & sin of Euler angles of ep(2:4)
c s0x, c0x - sin and cos for the xy-rotation
c s0, c0 - sin and cos for the zx'-rotation;
c ep is returned as rotated vector
c-----------------------------------------------------------------------
      integer ienvi
      common /cienvi/ ienvi
      dimension ep(4)
      pt2=ep(3)**2+ep(4)**2
      if(pt2.ne.0.)then
        pt=sqrt(pt2)
c xy plane 
        c0x=ep(3)/pt
        s0x=ep(4)/pt
c zx' plane
        pl=sqrt(pt2+ep(2)**2)
        s0=pt/pl
        c0=ep(2)/pl
      else
        c0x=1.
        s0x=0.
        pl=abs(ep(2))
        s0=0.
        c0=ep(2)/pl
      endif
      !~~~~~~~~~~check~~~~~~~~
      if(.not.(c0.le.0..or.c0.ge.0.))then
        call getMonitorFileIndex(ifmtx)
        write(ifmtx,'(a,$)')'ERROR psdefrot NaN (no stop);  '
        write(ifmtx,*)'ienvi ep pt2 c0 =',ienvi,ep,pt2,c0
      endif
      !~~~~~~~~~~check~~~~~~~~
      ep(2)=pl
      ep(3)=0.
      ep(4)=0.
      return
      end

c------------------------------------------------------------------------
      subroutine psdeftr(s,ep,ey)
c-----------------------------------------------------------------------
c Determines lorentz boosts parameters ey for cms defined by ep
c ey(1),ey(2),ey(3) refer to boost along z,x,y-axis respectively
c ey(i)=exp(rapidity(i))
c 4-vector ep is supposed to have mass squared s (not checked!!)
c ep is returned as boosted vector
c-----------------------------------------------------------------------
      dimension ey(3)
      double precision ep(4),s,ww,wp,wm
      do i=1,3
        if(ep(i+1).eq.0.d0)then
          ey(i)=1.
        else
          wp=ep(1)+ep(i+1)
          wm=ep(1)-ep(i+1)
          if(wp.gt.1.e-8.and.wm/wp.lt.1.e-8)then
            ww=s
            do l=1,3
              if(l.ne.i)ww=ww+ep(l+1)**2
            enddo
            wm=ww/wp
          elseif(wm.gt.1.d-8.and.wp/wm.lt.1.d-8)then
            ww=s
            do l=1,3
              if(l.ne.i)ww=ww+ep(l+1)**2
            enddo
            wp=ww/wm
          endif
          ey(i)=sngl(sqrt(wm/wp)) !wp'=wm' => e^y*wp=e^-y*wm
          ep(1)=wp*ey(i)
          ep(i+1)=0.
        endif
      enddo
      ep(1)=sqrt(s)
      return
      end

c------------------------------------------------------------------------
      function psidd(icc)
c-----------------------------------------------------------------------
c psidd - kink type decoder
c-----------------------------------------------------------------------
      if(icc.eq.0)then                    !g
        psidd=9
      elseif(iabs(icc).le.2)then          !u,u~,d,d~
        psidd=icc
      elseif(iabs(icc).eq.4)then          !s,s~
        psidd=icc/4*3
      elseif(iabs(icc).gt.10)then         !c,c~ etc.
        psidd=icc/10
      elseif(icc.eq.3)then                !ud
        psidd=1200
      elseif(icc.eq.-3)then               !u~d~
        psidd=-1200
      elseif(icc.eq.6)then                !uu
        psidd=1100
      elseif(icc.eq.-6)then               !u~u~
        psidd=-1100
      elseif(icc.eq.7)then                !dd
        psidd=2200
      elseif(icc.eq.-7)then               !d~d~
        psidd=-2200
      else
        psidd=0
        write (*,*)'psidd?????????',icc
      endif
      return
      end

cc------------------------------------------------------------------------
c       function pslam(s,a,b)
cc-----------------------------------------------------------------------
cc kinematical function for two particle decay - maximal pt-value
cc a - first particle mass squared,
cc b - second particle mass squared,
cc s - two particle invariant mass squared
cc-----------------------------------------------------------------------
c       pslam=.25/s*(s+a-b)**2-a
c       return
c       end
c


c------------------------------------------------------------------------
      subroutine pslcsh(wp1,wm1,wp2,wm2,samqt,iqq)
c-----------------------------------------------------------------------
c pslcsh - sh pomeron lc momentum sharing between two strings
c------------------------------------------------------------------------
      double precision amqt(4),yqm(4),yqm1(4),xlp(4),xlm(4),am23,sx,y2
     *,wp1,wp2,wm1,wm2,s,sq,psutz,yqmax,y,amjp,amjm,y1,s12,s34,x34
     *,psu(4),drangen,yqmax1
      dimension samqt(4)
      common/clcsh/ntrylcsh
#include "aaa.h"
      integer ncount,lcount,mcount,kcount
      data ncount/0/ mcount/0/
      save ncount,mcount
      call getMonitorFileIndex(ifmtx)

      s=wp1*wm1
      sq=dsqrt(s)
      ntrylcsh2=0
 2    ntrylcsh2=ntrylcsh2+1
      amqpt=samqt(1)+samqt(2)+samqt(3)+samqt(4)
      do i=1,4
        amqt(i)=dble(samqt(i))
        psu(i)=psutz(s,amqt(i)**2,(dble(amqpt)-amqt(i))**2)
        if(psu(i).le.0)then
          if(iqq.eq.-1)then
            do j=1,4
              samqt(j)=samqt(j)*0.9
            enddo
            if(ntrylcsh2.lt.100)goto 2 !to avoid redo whole psahot
          endif
          !~~~~~~~
          mcount=mcount+1 
          if(mcount.eq.1.or.mcount.eq.10.or.mcount.eq.100
     .      .or.mcount.eq.1000)then
            kcount=log10(1.*mcount)+1
            write(ifmtx,'(a,i1,a,5f8.3,i4)')'WARNING ',kcount,
     .     ' pslcsh A - Redo psahot;'
     .      ,sqrt(s),amqt(i),(dble(amqpt)-amqt(i)),wp1,wm1,iqq
          endif
          !~~~~~~~
          ntrylcsh=999
          return
        endif
        yqm(i)=dlog(sq/amqt(i)*psu(i))
      enddo
      yqmax=max(yqm(1),yqm(2))
c      write (ifmt,*)'pslcsh',sq,wp1,wm1,amqt,yqm,yqmax,psu

      ntrylcsh=0
 1    ntrylcsh=ntrylcsh+1

      if(ntrylcsh.gt.100)return
      y=yqmax*drangen(yqmax)
      j=int(1.5+rangen())
      if(y.gt.yqm(j))goto 1
      amjp=amqt(j)*dexp(y)
      amjm=amqt(j)*dexp(-y)
      do i=3,4
        am23=amqt(3-j)+amqt(7-i)
        sx=(am23+amjp)*(am23+amjm)
        psu(i)=psutz(s,amqt(i)**2,sx)
        if(psu(i).le.0)goto1
        yqm1(i)=dlog(sq/amqt(i)*psu(i))
      enddo
      yqmax1=max(yqm1(3),yqm1(4))
      if(max(yqm(3),yqm(4)).eq.0.d0)then
        if(iqq.eq.-1)then
          do j=1,4
            samqt(j)=samqt(j)*0.9
          enddo
          if(ntrylcsh2.lt.200)goto 2 !to avoid redo whole psahot
        endif
        !~~~~~~~
        ncount=ncount+1 
        if(ncount.eq.1.or.ncount.eq.10.or.ncount.eq.100
     .    .or.ncount.eq.1000)then
          lcount=log10(1.*ncount)+1
          write(ifmtx,'(a,i1,a,2e10.3,i4)')'WARNING ',lcount,
     .   ' pslcsh B - Redo psahot;',yqm(3),yqm(4),iqq
        endif
        !~~~~~~~
        ntrylcsh=999
        return    
      endif
      if(drangen(yqmax1).gt.yqmax1/max(yqm(3),yqm(4)))goto 1

      y1=yqmax1*drangen(yqmax1)
      j1=int(3.5+rangen())
      if(y1.gt.yqm1(j1))goto 1


      amjp1=amqt(j1)*exp(y1)
      amjm1=amqt(j1)*exp(-y1)
      s12=(amqt(3-j)+amjp)*(amqt(3-j)+amjm)
      s34=(amqt(7-j1)+amjp1)*(amqt(7-j1)+amjm1)
      y2=dlog(sq/(amqt(3-j)+amjp)*psutz(s,s12,s34))


      xlp(j)=amqt(j)/sq*dexp(y+y2)
      xlm(j)=amqt(j)/sq*dexp(-y-y2)
      xlp(3-j)=amqt(3-j)/sq*dexp(y2)
      xlm(3-j)=amqt(3-j)/sq*dexp(-y2)
      x34=1.-xlm(1)-xlm(2)
      xlm(7-j1)=x34/(1.+amjp1/amqt(7-j1))
      xlm(j1)=x34-xlm(7-j1)
c      write (ifch,*)'xlc',xlp(1),xlp(2),xlm(3),xlm(4)
      if(drangen(x34).gt.(xlp(1)*xlp(2))**(-alpqua(1))
     *                  *(xlm(3)*xlm(4))**(-alpqua(2))*
     *  (xlp(j)*(1.d0-xlp(j))*xlm(j1)*(1.d0-xlm(j1))))goto 1

      wp2=xlp(2)*wp1
      wp1=xlp(1)*wp1
      wm2=xlm(4)*wm1
      wm1=xlm(3)*wm1
c      write (ifch,*)'wp1,wm1,wp2,wm2',wp1,wm1,wp2,wm2
      if(.not.(wm1.le.0..or.wm1.ge.0.))then
      !print*,'******** NaN catch wm1 *********',wm1
      goto 1
      endif
      if(.not.(wp1.le.0..or.wp1.ge.0.))then
      !print*,'******** NaN catch wp1 *********',wp1
      goto 1
      endif
      return
      end

c------------------------------------------------------------------------
      function psnorm(ep)
c-----------------------------------------------------------------------
c 4-vector squared calculation
c-----------------------------------------------------------------------
      double precision sm2,ep(4)
      sm2=(ep(1)-ep(2))*(ep(1)+ep(2))
      do i=3,4
        sm2=sm2-ep(i)**2
      enddo
      psnorm=sngl(sm2)
      return
      end

c------------------------------------------------------------------------
      subroutine psrotat(ep,s0x,c0x,s0,c0)
c-----------------------------------------------------------------------
c psrotat - spacial rotation to the lab. system for 4-vector ep
c s0, c0 - sin and cos for the zx-rotation;
c s0x, c0x - sin and cos for the xy-rotation
c-----------------------------------------------------------------------
      dimension ep(4),ep1(3)

      ep1(3)=ep(4)               ! y
      ep1(2)=ep(2)*s0+ep(3)*c0   ! x
      ep1(1)=ep(2)*c0-ep(3)*s0   ! z

      ep(2)=ep1(1)                   ! z
      ep(4)=ep1(2)*s0x+ep1(3)*c0x    ! y
      ep(3)=ep1(2)*c0x-ep1(3)*s0x    ! x
      return
      end

c-----------------------------------------------------------------------
      double precision function pssalf2(qq)
c-----------------------------------------------------------------------
      double precision qq
      pssalf2=pssalf(sngl(qq))
      end

c-----------------------------------------------------------------------
      function pssalf(qq)
c-----------------------------------------------------------------------
c pssalf - effective qcd coupling (alpha_s/2/pi)
c              qq means Q2 / qcdlam !!!!!!!!!
c-----------------------------------------------------------------------
c saturate at 1/2pi (not reached for q2min=0.25) but not 0.5/2pi
c because this leads to a too small
c rise of F2 at small x a large Q2 and then a too low jet incl. xs
c may be because the saturation is not taken into account in psuds for
c instance (TP20170215)
c TP20200306 : introduce smoother transition using
c     alpha_s~exp(-Q2/(4.adskap**2)) below Q2=0.9 GeV**2
c     (According to Brodsky et al Phys. Lett. B 750, 528 (2015);
c      J. Phys. G 44, 105005 (2017))         
c-----------------------------------------------------------------------
#include "aaa.h"
#include "sem.h"

      qq0=qq*qcdlam
      if(qq0.lt.0.9) then
        pssalf=exp(-qq/adskap)/(2.*pi)
      else
        nf = noflav(qq0)       !KW1811
        xqq = qq0/qcdlambda(nf)  !TP2007
        !other option: 
        !nf = nofeff  !nofeff defined in bas.f
        !xqq = qq0/qcdlam
        pssalf=2./(11.-float(nf)/1.5)/log(xqq) 
        if(pssalf.lt.0.or.pssalf.gt.1./(2.*pi))then
          pssalf=1./(2.*pi)
        endif
      endif
      return
      end

c------------------------------------------------------------------------
      subroutine pstrans(ep,ey,jj)
c-----------------------------------------------------------------------
c pstrans - lorentz boosts according to the parameters ey ( determining
c shift along the z,x,y-axis respectively (ey(1),ey(2),ey(3)))
c jj=1 - inverse transformation to the lab. system;
c jj=-1 - direct transformation
c-----------------------------------------------------------------------
      dimension ey(3),ep(4)

      if(jj.eq.1)then
c lorentz transform to lab. system according to 1/ey(i) parameters
        do i=1,3
          if(ey(4-i).ne.1.)then
            wp=(ep(1)+ep(5-i))/ey(4-i)
            wm=(ep(1)-ep(5-i))*ey(4-i)
            ep(1)=.5*(wp+wm)
            ep(5-i)=.5*(wp-wm)
          endif
        enddo
      else
c lorentz transform to lab. system according to ey(i) parameters
        do i=1,3
          if(ey(i).ne.1.)then
            wp=(ep(1)+ep(i+1))*ey(i)
            wm=(ep(1)-ep(i+1))/ey(i)
            ep(1)=.5*(wp+wm)
            ep(i+1)=.5*(wp-wm)
          endif
        enddo
      endif
      return
      end

c------------------------------------------------------------------------
      double precision function psuds(q,m) ! SL Sudakov interpolation !KW1811  
c-----------------------------------------------------------------------
      double precision psudi2
      psuds=psudi2(dble(q),m)
      end

c------------------------------------------------------------------------
      double precision function psudi2(q,m) ! SL Sudakov interpolation !KW1811  
c-----------------------------------------------------------------------
c psudi2 - spacelike sudakov formfactor interpolation
c q - maximal value of the effective momentum,
c m - type of parton (0 - g, 1,2, etc. - q)
c-----------------------------------------------------------------------
#include "aaa.h"
      double precision q
      dimension wi(2)
      common /psar15/ sudx(40,2)
#include "sem.h"

      j=min(iabs(m),1)+1

      if(q.gt.q2zmin)then
        !i=1,40: qi=q2zmin*exp(.5*(i-1))
        xi=log(q/q2zmin)*2.+1.
        i=int(xi)
        if(i.lt.1)i=1
        if(i.gt.39)i=39
        frac=xi-i
        wi(1)=1-frac
        wi(2)=frac
        psudi2=0
        do i1=1,2
          psudi2=psudi2+dble(sudx(i+i1-1,j)*wi(i1))
        enddo
        psudi2=exp(-psudi2)
      else
        psudi2=1.d0
      endif
      return
      end

c-----------------------------------------------------------------------
      subroutine makeSudx
c-----------------------------------------------------------------------
#include "aaa.h"
      double precision sudspa
      common /psar15/ sudx(40,2)
      common /csudy/ sudy(40,2)
      data ncount /0/
      save ncount
      data ((sudy(i,j),j=1,2),i=1,40)   !copied from version 3.246
     .  /    0.40861 ,     0.14706
     .  ,    0.70676 ,     0.27443
     .  ,    1.08492 ,     0.44145
     .  ,    1.54306 ,     0.64634
     .  ,    2.07839 ,     0.88665
     .  ,    2.68651 ,     1.15969
     .  ,    3.36237 ,     1.46285
     .  ,    4.10093 ,     1.79370
     .  ,    4.89745 ,     2.15007
     .  ,    5.74765 ,     2.53003
     .  ,    6.64768 ,     2.93188
     .  ,    7.59412 ,     3.35413
     .  ,    8.58393 ,     3.79545
     .  ,    9.61440 ,     4.25467
     .  ,   10.68313 ,     4.73076
     .  ,   11.78800 ,     5.22278
     .  ,   12.92710 ,     5.72991
     .  ,   14.09870 ,     6.25140
     .  ,   15.30127 ,     6.78658
     .  ,   16.53340 ,     7.33483
     .  ,   17.79380 ,     7.89559
     .  ,   19.08132 ,     8.46836
     .  ,   20.39487 ,     9.05265
     .  ,   21.73346 ,     9.64802
     .  ,   23.09617 ,    10.25409
     .  ,   24.48214 ,    10.87047
     .  ,   25.89060 ,    11.49681
     .  ,   27.32078 ,    12.13278
     .  ,   28.77201 ,    12.77809
     .  ,   30.24363 ,    13.43244
     .  ,   31.73504 ,    14.09556
     .  ,   33.24565 ,    14.76720
     .  ,   34.77492 ,    15.44713
     .  ,   36.32235 ,    16.13511
     .  ,   37.88745 ,    16.83092
     .  ,   39.46976 ,    17.53438
     .  ,   41.06884 ,    18.24528
     .  ,   42.68428 ,    18.96345
     .  ,   44.31570 ,    19.68870
     .  ,   45.96270 ,    20.42087 
     .  /
      ncount=ncount+1
      if(ncount.ne.1)return
      iopt=1  !  <------------------
      if(iopt.eq.2)then
        write(ifmt,'(75a/3a/75a)')('#',i=1,75),'WARNING makeSudx '
     .             ,'Temporarily we take the table sudx '
     .             ,'as tabulated from 3.246. '
     .             , ('#',i=1,75)
      endif
      do i=1,40
        qi=q2zmin*exp(.5*(i-1))
        do j=1,2
          if(iopt.eq.1)sudx(i,j) = -log( sudspa(qi,j-1) )
          if(iopt.eq.2)sudx(i,j) = sudy(i,j)  !test
        enddo 
      enddo
      return
      end

c-----------------------------------------------------------------------
      subroutine printSudx !print Sudakov for different definitions
c-----------------------------------------------------------------------
#include "aaa.h"
      common /psar15/ sudx(40,2)
      common /csudy/ sudy(40,2)
      double precision qii,sudspaSO
      call makeSudx
      write(ifmt,'(70a)')('-',i=1,70)
      write(ifmt,'(3x,a/3x,a,a/3x,a/3x,a)')
     .'Sudakov Delta(q2) comparing new and old (Nexus) method,'
     .,'considering evolution from q2ini to q2 ; '
     .,'we print -log(Delta(q2))'  
     .,'        q      gluon   gluon   gluon      quark   quark   quark'
     .,'                old     new     serg       old     new     serg'
      do i=1,40
        qi=q2zmin*exp(.5*(i-1))
        qii=dble(qi)
        q=sqrt(qi) 
        if(mod(i-1,2).eq.0.and.q.lt.2000)then 
          write(ifmt,'(f12.2,3x,3f8.2,7x,3f8.2)')
     .    q , sudy(i,1) , sudx(i,1), -log( sudspaSO(qii,0) )
     .      , sudy(i,2) , sudx(i,2), -log( sudspaSO(qii,1) )
        endif
      enddo
      write(ifmt,'(70a)')('-',i=1,70)
      stop'printSudx'
      end


c-----------------------------------------------------------------------
      double precision function epscutSud2(qq)
c-----------------------------------------------------------------------
      double precision qq
      !-------------------------------------
      !if the following line is activated: 
      !change also qqx=... in sem.f
      !-------------------------------------
        !real val
        !call getQ2ini(val)
        !epscutSud2=  dble(val) / qq   !3312
        !Before uncommenting the above lines:
        !change explicit epscutSud2 dependence in  psevi2
      !-------------------------------------
        epscutSud2= .01d0             !SO
      !-------------------------------------
      dmy=qq
      end
      function epscutSud(qq)
      double precision epscutSud2
      epscutSud=epscutSud2(dble(qq))
      end
c-----------------------------------------------------------------------
      function qminEvol(q1,q2)
c-----------------------------------------------------------------------
      qminEvol=q1          !SO 
      !qminEvol=max(q1,q2) !3312
      dmy=q2
      end


c-----------------------------------------------------------------------
      double precision function sudspa(q,k) 
c-----------------------------------------------------------------------
c
c   In epos3312/KW/suda.f there is a version with pssalf(q2*(1-z)) 
c         but that is not the right way, since in pssalf 
c            we need q2 and not q2*(1-z)
c
c-----------------------------------------------------------------------
c spacelike sudakov formfactor from q2ini to q
c k = 0 : evolving parton is a gluon ; else it's a (anti) quark
c pij : integrated splitting functions, see 265
c zma, zmi : born max and min for z
c-----------------------------------------------------------------------
#include "aaa.h"
#include "sem.h"
      common /ar3/    x1(7),a1(7)
      double precision zma,zmi,delz
      double precision epscutSud2 !function
      q2minSud           =  q2zmin      ! 1.0 !(SO)   !
      qcdlamSud          =  qcdlam      ! .04 !(SO)   !
      !--------
      !remarks, concerning q2ini difference -log(sudspa) to -log(sudspaSO) (between Q=1 and Q=1800)
      ! using the same epscutSud2, q2minSud, qcdlamSud ... maxi +10% diff (integr more accurate in sudspa)
      ! using the same epscutSud2, q2minSud .............. again maxi +10% difference
      ! using the same epscutSud2 ........................ again maxi +10% difference
      ! default ......................................... maxi +100% difference
      !-------- 
      sudspa=0.d0
      if(q.gt.q2minSud)then
        do i=1,7
        do m=1,2
          ! the integrand is roughly 1/t so we define a variable v
          ! v=log(t)  dv = dt/t
          v1 = log(q2minSud) !integration
          v2 = log(q)     !limits
          x = x1(i)*(2.*m-3.)
          v = 0.5*(v1+v2-x*(v1-v2))  
          t=exp(v)
          if(ioangord.ne.0)stop'ERROR 14022019' !angle ordering: we don't do that any more
          zma = 1.d0 -  epscutSud2(dble(t))
          zmi = epscutSud2(dble(t))         
          alpha=pssalf(t/qcdlamSud) 
          delz=zma-zmi
          if(delz.le.0.) then
            pij=0
          elseif(k.eq.0) then 
            pijA=(primSFg2g(zma)-primSFg2g(zmi))*0.5     !g->g   factor 0.5, see pub/19epos4/serg2.pdf, page 14
            pijB=primSFg2q(zma)-primSFg2q(zmi) !g->q /Nf
            pij = pijA + pijB*noflav(t)                  !g->q  NO factor 2, see pub/19epos4/serg2.pdf, page 14
          else
            pijA=primSFq2q(zma)-primSFq2q(zmi) !q->q 
            pijB=0                                       !only q->q, see pub/19epos4/serg2.pdf, page 14
            pij = pijA + pijB
          endif
          pij = pij * alpha
          !if(m.eq.1)write(ifmt,'(i3,f14.2,9f9.3)')k,t,pij1,pijA,pij2,pijB
          sudspa=sudspa+dble(a1(i)*pij)   !1/t absent due to varible change.   alpha=alpha_s/2/pi
          if( (1.d0-zma)/(1.d0-zmi) .eq. 0.d0 )then
            stop'ERROR 17112018a' !One needs to improve precision (beyond double)
          endif
        enddo
        enddo
        sudspa=sudspa*dble(0.5*(v2-v1)) !jacobien
        sudspa=dexp(-sudspa)
        if(.not.(sudspa.le.0.d0.or.sudspa.ge.0.d0))then
          write(*,*)q,q2minSud,zma,zmi,k
          stop 'sudspa,sem.f'
        endif
      else
        sudspa=1.d0
      endif
      
      return
      end

c-----------------------------------------------------------------------
      double precision function sudspaSO(q,k) 
c-----------------------------------------------------------------------
c spacelike sudakov formfactor from q2minSud to q
c k = 0 : evolving parton is a gluon ; else it's a (anti) quark
c-----------------------------------------------------------------------
      implicit double precision (a-h,o-z)
      j=k+1         ! type of parton (1 - g, 2 - q)
      alm=.04d0          !lambda_qcd squared  
      epscut=.01d0            !pt-resolution scale (technical)
      q2minSud=1       !(value used directly)
c-----------------------------------------------------------------------
c Simple formula with nf=3 and using the limit eps (=epscut) -> 0 for int_eps^(1-eps) dx ...
c We have Sudakov=exp(-AB), with
c   A = int dQ2/Q2 1/4.5*log(Q2/L2) = 2/9*log( log(Q2/L2)/log(q2minSud/L2) )
c for gluon:
c   B = int_eps^(1-eps) dx ( 1/2 Pg->g + nf*Pg->q )
c     = int_eps^(1-eps) dx ( 3*(z/(1-z)+(1-z)/z+z*(1-z)) + nf/2*(z**2+(1-z)**2) )
c     = -6 * ( ln(eps) + 3/4 )
c  AB = log( log(Q2/L2)/log(q2minSud/L2) ) * (-4/3) * ( ln(eps) + 3/4 )
c  for quark: 
c   B = int_eps^(1-eps) dx ( Pq->q )
c     = int_eps^(1-eps) dx ( 4/3 * (1+z**2)/(1-z) )
c     = -8/3 * ( ln(eps) + 3/4  )
c  AB = log( log(Q2/L2)/log(q2minSud/L2) ) * (-16/27) * ( ln(eps) + 3/4 )
c-----------------------------------------------------------------------
      if(q.gt.1.d0)then
       qgsudx=dlog(dlog(q/alm)/dlog(1.d0/alm))*(.75d0+dlog(epscut))
       if(j.eq.1)then
        qgsudx=exp(qgsudx/.75d0)
       else
        qgsudx=exp(qgsudx*16.d0/27.d0)
       endif
      else
       qgsudx=1.d0
      endif
      sudspaSO=qgsudx
      end

c-----------------------------------------------------------------------
c  primitives of splitting functions  primSF !KW1811
c-----------------------------------------------------------------------
      function primSFg2g(z) !g->g
      double precision z
      !SF: 6*(z/(1-z)+(1-z)/z+z(1-z))
      primSFg2g=6*(log(z)-log(1-z)-2*z+z**2/2-z**3/3)
      end
c---------------------------------------------------------------------------------------
      function primSFg2q(z) !g->q  divided by Nf !!!!
      double precision z
      !SF: 0.5*(z**2+(1-z)**2)
      primSFg2q=(z**3-(1-z)**3)/6 
      end
c---------------------------------------------------------------------------------------
      function primSFq2q(z) !q->q  
      double precision z
      !SF: 4/3 * (1+z**2)/(1-z)  
      primSFq2q=-2./3.*(z**2+2*z+4*log(1-z))
      end
c---------------------------------------------------------------------------------------
      function primSFq2g(z) !q->g    avec z->1-z)
      double precision z
      !SF: 4/3 * (1+(1-z)**2)/z  
      primSFq2g=2./3.*((1-z)**2+2*(1-z)+4*log(z))
      end

c------------------------------------------------------------------------
      function psfapj(x,j,l)  !for j,l being given individual flavors
c-----------------------------------------------------------------------
      double precision x
#include "sem.h"

      if(j.eq.0)then
        if(l.eq.0)then
          psfapj=((1.d0-x)/x+x/(1.d0-x)+x*(1.d0-x))*6.d0
        else
          psfapj=(x**2+(1.d0-x)**2)*0.5 
        endif
      else
        if(l.eq.0)then
          psfapj=(1.d0+(1.d0-x)**2)/x/.75d0
        else
          psfapj=0
          if(j.eq.l)psfapj=(x**2+1.d0)/(1.d0-x)/.75d0
        endif
      endif
      return
      end

c------------------------------------------------------------------------
      function psfap(qq,x,j,l)  !for j,l being gluon or quark
c-----------------------------------------------------------------------
c psfap - altarelli-parisi splitting function (multiplied by x)
c qq virtuality squared
c x - light cone momentum share value,
c j - type of the parent parton (0-g;1,2,etc.-q)
c l - type of the daughter parton (0-g;1,2,etc.-q)
c-----------------------------------------------------------------------
      double precision x
#include "sem.h"

      if(j.eq.0)then
        if(l.eq.0)then
          psfap=((1.d0-x)/x+x/(1.d0-x)+x*(1.d0-x))*6.d0
        else
          psfap=(x*x+(1.d0-x)**2)*noflav(qq)  !actually 0.5 * 2*noflav
        endif
      else
        if(l.eq.0)then
          psfap=(1.d0+(1.d0-x)**2)/x/.75d0
        else
          psfap=(x*x+1.d0)/(1.d0-x)/.75d0
        endif
      endif
      return
      end

c------------------------------------------------------------------------!bg ga
      function psfapel(x,m)
c-----------------------------------------------------------------------
c psfap - altarelli-parisi splitting function (multiplied by x) for q->q+emission gamma
c x - light cone momentum share value,
c j - type of the parent parton (0-g;1,2,etc.-q)
c-----------------------------------------------------------------------
      double precision x
#include "sem.h"

      if(iabs(m).eq.1) then
        e2=4./9. !bg square electric charge
      elseif(iabs(m).eq.2) then
        e2=1./9.
      elseif(iabs(m).eq.3) then
        e2=1./9.
      elseif(iabs(m).eq.4) then
        e2=4./9.
      elseif(iabs(m).eq.5) then
        e2=1./9.
      else
        stop 'wrong m value'
      endif
      psfapel=e2*(1.d0+x**2)/(1.d0-x)
      return
      end

c------------------------------------------------------------------------
      double precision function psutz(s,a,b)
c-----------------------------------------------------------------------
c psutz - decay of mass^2 s into two particles having mT^2 of a and b.
c  P+ = sqrt(s) -> z*sqrt(s) + (1-z)*sqrt(s)           (1)
c  P- = sqrt(s) -> a/(z*sqrt(s)) + b/((1-z)*sqrt(s))   (2)
c psutz = z (LC momentum fraction) of first decay particle
c Equation (2) is a quadratic equation for z, z**2+B*z+C=0
c  with B = -(s+a-b)/s and C = a/s
c In the following x=-B/2 and dx=sqrt((B/2)**2-C), so z=x+sqrt(dx) 
c-----------------------------------------------------------------------
      common/cpsutz/dxpsutz
      double precision a1,b1,s1,x,dx,s,a,b
      a1=dsqrt(a)
      b1=dsqrt(b)
      s1=dsqrt(s)
      x=(1.d0+(a1-b1)*(a1+b1)/s)/2.d0

      dx=(x-a1/s1)*(x+a1/s1)
      dxpsutz=dx !just for print in psahot
c      x=.5*(1.+(a-b)/s)
c      dx=(x*x-a/s)
      if(dx.gt.0.d0)then
        x=x+dsqrt(dx)
      else !take p+=p- for particle 1
        x=a1/s1
      endif
      psutz=min(0.999999999d0,x)
      return
      end

c------------------------------------------------------------------------ !bg enlever?
c      double precision function psutzQ(s,a,alpha)
c-----------------------------------------------------------------------
c psutz - kinematical function for two particle decay - light cone moment
c a - pt^2
c s - two particle invariant mass
c alpha - take into account virtuality with Q^2_ini = alpha*pt^2
c Remark : For the moment, for simplification, if the second particule is a photon, alpha=0
c-----------------------------------------------------------------------
c      double precision x,s,a

c      if((s/a).lt.(2.*alpha+4.))then
c        write(*,*)'alpha too big',s/a,alpha !bg enlever
c        alpha=0.
c      endif

c      x=0.5d0*(1.d0+dsqrt(dble(1.-4.*(1.+alpha)/(2.*alpha+s/a))))
c      if(x.gt.1.d0.or.x.lt.0.d0)stop'wrong x, psutzQ in sem.f'
c      psutzQ=min(0.999999999d0,x)
c      return
c      end




c------------------------------------------------------------------------ !bg enlever?
      double precision function psutzQ(s,a,q)
c-----------------------------------------------------------------------
c psutz - kinematical function for two particle decay - light cone moment
c a - pt^2
c s - two particle invariant mass
c q - initial virtuality for the timelike parton
c-----------------------------------------------------------------------
      double precision x,dx,s,a

      dx=1.d0-4.d0*(dble(q)+a)/(s+dble(q))
      if(dx.lt.-0.0001d0) then
        write(*,*)'dx, s , pt^2 , Q^2 :',dx,s,a,q
        stop 'dx<0 psutzQ, sem.f'
      elseif(dx.lt.0.d0) then  !here pt^2~s/4 and then dx -> 0
        x=0.5d0
      else
        x=0.5d0*(1.d0+dsqrt(dx))
      endif

      if(x.gt.1.d0.or.x.lt.0.d0)stop'wrong x, psutzQ in sem.f'
      psutzQ=min(0.999999999d0,x)
      return
      end




c------------------------------------------------------------------------
      block data ptdata
c-----------------------------------------------------------------------
c constants for numerical integration (gaussian nodes and weights)
c   x1(7),a1(7)  gaussian legendre for n=14
c   xh14(14),wh14(14) gauss hermite for n=28 
c-----------------------------------------------------------------------
      common /ar3/ x1(7),a1(7)
      common /ar4/ x4(2),a4(2)
      common /ar5/ x5(2),a5(2)
      common /ar8/ x2(4),a2
      common /ar9/ x9(3),a9(3)
      common /ar14/ x14(14),a14(14)
      common /arh5/xh5(5),wh5(5)
      common /arh6/xh6(6),wh6(6)
      common /cgauss7/ xgauss7(7),wgauss7(7)

      data xgauss7/.9862838,.9284349,.8272013,.6872929,.5152486,
     *.3191124,.1080549/
      data wgauss7/.03511946,.08015809,.1215186,.1572032,
     *.1855384,.2051985,.2152639/
      data x1/.9862838,.9284349,.8272013,.6872929,.5152486,
     *.3191124,.1080549/
      data a1/.03511946,.08015809,.1215186,.1572032,
     *.1855384,.2051985,.2152639/
      data x2/.00960736,.0842652,.222215,.402455/
      data a2/.392699/
      data x4/ 0.339981,0.861136/
      data a4/ 0.652145,0.347855/
      data x5/.585786,3.41421/
      data a5/.853553,.146447/
      data x9/.93247,.661209,.238619/
      data a9/.171324,.360762,.467914/
      !from https://keisan.casio.com/exec/system/1280624821
      data x14/0.0550792899,0.1645692821,0.2720616276,0.3762515161
     .,0.475874225,0.5697204718,0.656651094,0.735610878,0.8056413709
     .,0.8658925226,0.9156330264,0.9542592806,0.9813031654,0.9964424976/
      data a14/0.110047013,0.1087111923,0.1060557659,0.1021129676
     .,0.096930658,0.0905717444,0.0831134172,0.0746462142,0.065272924
     .,0.0551073457,0.0442729348,0.0329014278,0.0211321126,0.0091242826/
      !from https://keisan.casio.com/exec/system/1329114617
      data xh5/0.342901,1.03661,1.75668,2.53273,3.4361591188/
      data wh5/0.610862,0.240138,0.0338743,0.00134364,7.64043E-006/
      data xh6/0.314240,0.947788,1.59768,2.27950,3.02063,3.88972/
      data wh6/0.570135,0.260492,0.0516079,0.00390539,8.57368E-5
     .,2.65855E-7/

      end

c------------------------------------------------------------------------
      block data gjdata
c-----------------------------------------------------------------------
c constants for numerical integration (gauss-Jacobi weights)
c for Integral ( A <= x <= B ) (B-x)^alpha (x-A)^beta f(x) dx
c is to be approximated by
c        Sum ( 1 <= i <= order ) w(i) * f(x(i))
c with A=-1, B=1, alpha=beta=-0.75
c-----------------------------------------------------------------------
      double precision dx10,da10
      common /ar10/ dx10(10),da10(10)
      data dx10/  0.7948197113078417D-01,0.2364390429616058D0
     .           ,0.3874261741092539D0,0.5286310257944824D0
     .           ,0.6564882417344481D0,0.7677694505306535D0
     .           ,0.8596647146888806D0,0.9298532055552277D0
     .           ,0.9765599789306546D0,0.9985637292444355D0/
      data da10/ 0.7489979559374379D0,0.3456591409996251D0
     .          ,0.2630010942298393D0,0.2228908167854909D0
     .          ,0.1989151597748371D0,0.1832757320590487D0
     .          ,0.1727459069916109D0,0.1657467734816896D0
     .          ,0.1614405906016221D0,0.1593843834309191D0/
      end


c------------------------------------------------------------------------
      double precision function fzeroGluZZ(z,ipt)   ! former psftild
c-----------------------------------------------------------------------
c    (no screening)
c
c    fzeroGluZZComplete = z * fzeroGluZZ * z^(-1.-dels)
c
c
c fzeroGluZZ = gamhad
c  * int(du)/u u^(betff+dels+1) (1-u)^alplea * EsaturGluonTil(z/u)
c (betff=-alppar)
c
c z - light cone x of the gluon,
c ipt - pro/tar
c
c No screening in vertex function !
c-----------------------------------------------------------------------
      double precision xpmin,xp,z,Fsea,EsaturGluonTil
#include "aaa.h"
      common /ar3/   x1(7),a1(7)
#include "sem.h"

      fzeroGluZZ=0.d0
      xpmin=z
      xpmin=xpmin**dble(betff(ipt)+dels(ipt))
      do i=1,7
      do m=1,2
        xp=(.5d0*(1d0+xpmin+dble(2*m-3)*dble(x1(i))*(1d0-xpmin)))**(1d0/
     *            dble(betff(ipt)+dels(ipt)))
c        fzeroGluZZ=fzeroGluZZ+a1(i)*Fsea(xp,ipt)*(z/xp)**dble(-1.-dels)
c     *                             *EsaturGluonTil(z/xp,ipt)
c     *                             /xp                            !convolution
c     *                             *xp/xp**dble(betff(ipt)+dels)  !Jacobian
        fzeroGluZZ=fzeroGluZZ+a1(i)*Fsea(xp,ipt)
     *                             *EsaturGluonTil(z/xp,q2cmin(ipt),ipt)
     *                             *xp**(1.-betff(ipt))
      enddo
      enddo
      fzeroGluZZ=fzeroGluZZ*.5d0*(1d0-xpmin)/dble(betff(ipt)+dels(ipt))
c simplified calculation
c     *                                     *z**dble(1.+dels)
      return
      end

c------------------------------------------------------------------------
      double precision function FzeroQuaZZ(z,ipt)     ! former psftile
c-----------------------------------------------------------------------
c    (no screening)
c
c    FzeroQuaZZComplete = z * FzeroQuaZZ * z^(-1-dels)
c
c  gamsoft = 8*pi*s0*gampar*gamtilde
c integration over semihard pomeron light cone momentum share xp==u
c
c FzeroQuaZZ = gamhad *
c   * int(du)/u u^(betff+dels+1) (1-u)^alplea * EsaturQuarkTil (z/u)
c (betff=-alppar)
c      
c z - light cone x of the quark,
c ipt - proj/tar
c
c No screening in vertex function !
c-----------------------------------------------------------------------
      double precision xpmin,xp,z,Fsea,EsaturQuarkTil
      common /ar3/   x1(7),a1(7)
#include "aaa.h"
#include "sem.h"

      if(ipt.eq.1)then
        k=iclpro
      else
        k=icltar
      endif

      FzeroQuaZZ=0.d0
      xpmin=z
      xpmin=xpmin**dble(betff(ipt)+dels(ipt))
      do i=1,7
      do m=1,2
        xp=(.5d0*(1d0+xpmin+dble(2*m-3)*dble(x1(i))*(1d0-xpmin)))**(1d0/
     *            dble(betff(ipt)+dels(ipt)))
c        FzeroQuaZZ=FzeroQuaZZ+a1(i)*Fsea(xp,ipt)*(z/xp)**dble(-1.-dels(ipt))
c     *                             * EsaturQuarkTil (z/xp,ipt)
c     *                             / xp                            !convolution
c     *                             * xp/xp**dble(betff(ipt)+dels(ipt))  !Jacobian
c simplified calculation
        FzeroQuaZZ=FzeroQuaZZ+a1(i)*Fsea(xp,ipt)
     *                        *EsaturQuarkTil(z/xp,q2cmin(ipt),ipt,999)
     *                             *xp**(1.-betff(ipt))
      enddo
      enddo
      FzeroQuaZZ=FzeroQuaZZ*.5d0*(1d0-xpmin)/dble(betff(ipt)+dels(ipt))
c simplified calculation
c     *                                     *z**dble(1.+dels(ipt))
      return
      end


c########################################################################
c########################################################################
      subroutine psaini
c########################################################################
c########################################################################

#include "aaa.h"
#include "par.h"
#include "sem.h"
#include "tab.h"
      integer ioTestFact
      common /ctestfact/ ioTestFact
      integer laddTestFact
      common /claddtestfact/ laddTestFact
      integer modeDeltaFzero
      common /cmodeDeltaFzero/ modeDeltaFzero
      common /psar2/  edmax,epmax


      call utpri('psaini',ish,ishini,4)


c number of flavors in fragmentation not less than active flavor in hard string

      if(inicnt.eq.1)then
        if(q2fin.le.4.*qcdlam)stop'ERROR q2fin =< 4.*qcdlam forbiden'
        adskap=0.9/qcdlam/(log(log(0.9/qcdlam))-log(4.*pi/9.)) !parameter for alpha_s~exp(-Q2/(4.adskap**2)) below Q2=0.9 GeV**2 (According to Brodsky et al Phys. Lett. B 750, 528 (2015); J. Phys. G 44, 105005 (2017))          
        qcdlambda(3)=qcdlam     !lambda_qcd squared for 3 flavors
        qcdlambda(4)=qcdlam*0.68 !lambda_qcd squared for 4 flavors
        qcdlambda(5)=qcdlam*0.38 !lambda_qcd squared for 5 flavors
        if(abs(qcmass-1.27).ge.0.1.or.abs(qbmass-4.18).ge.0.1)then
          write(ifmt,'(a,f5.2)')'Please check qcdlambda for alpha_s fit'
          print *,qcmass,qbmass
c          stop
        endif
        edmax=edmaxi            !1.e12     defined in epos-bas
        epmax=epmaxi            !1.e12     defined in epos-bas
      endif


      call setBasics

      if(inicnt.gt.1)return

      if(iappl.le.1.or.iappl.eq.8)then
        !call printSudx 
        call mkEvTable       !  ev.i  -- evolution tables
        !call testEvolution
        if(iappl.le.1)call mkCsOmTables
        call mkParamTables
        if(ioTestFact.ge.6.and.ioTestFact.le.11)then
          modeDeltaFzero=1      !fzero ~ delta(1-x)
        else
          modeDeltaFzero=0
        endif
        if(ioTestFact.ge.4.and.ioTestFact.le.11)noflit=0.5

      endif

      

cTEST compare psjti0 psjti pijet psjet
c      q2cmin(1)=1.5 ! 6.25
c      q2cmin(2)=1.5 ! 6.25
c      qq=max( q2cmin(1), q2cmin(2))
c      do i=1,20
c      ss=10*2**(i-1)
c      call ipoCSzeroTables(q2cmin)
c      call psjti0(ss,sj,sjb,0,0)    
c      sj2=  psjti(1,q2cmin(1),qq,ss  ,0,0,  0)
c      sj3=    pijet( 2 ,ss,0,0)  
c     .        +pijet( 1 ,ss,0,0)  
c     .        +pijet( -1,ss,0,0)  
c     .        +pijet( 0 ,ss,0,0)  
c      q1=q2cmin(1)
c      q2=q2cmin(2)
c      sj4=psjet(q1,q2,qq,ss,0,0,0)
c     .   +psjet1(q1,q2,qq,ss,0,0,0)
c     .   +psjet1(q2,q1,qq,ss,0,0,0)
c     .   +psborn(q1,q2,qq,ss,0,0,0,0)
c      if(i.eq.1)print*,' '
c      print*,'TEST',ss,sj,sj2,sj4
c      enddo 
c      stop'TEST'
cTEST

      end

c------------------------------------------------------------------
      subroutine testEvolution
c------------------------------------------------------------------
      common/ar3/    x1(7),a1(7)
#include "aaa.h"
#include "par.h"
#include "sem.h"
      double precision zx,z,xx,xmin,xm,epscutSud2
      q2cmin(1)=q2nmin
      q2cmin(2)=q2nmin
      call ipoBasics(q2cmin) 
      nrbins=15
      xminim=1e-5
      xmaxim=1
          qq=1200 !q2
      xx=xminim*(xmaxim/xminim)**(( 1 -.5)/nrbins)
      xmin=xx/(1.d0-epscutSud2(dble(qq)))
      xm=max(xmin,0.3d0)
      if(xm.gt.xmin)then  ! xm=0.3  ->  xmin=10^-5  : 1 -> 500000    
         do i=1,7
         do m=2,2
            zx=xx+(xm-xx)
     &         *((xmin-xx)/(xm-xx))**(.5d0-(dble(m)-1.5d0)*dble(x1(i)))
            z=xx/zx
            print*,'TEST',sngl(z)
     7 ,    psevi(q2cmin(1),qq,z,1,1)  
     7 ,    psevi(q2cmin(1),qq,z,2,1)  
         enddo
         enddo
      endif
      stop'in testEvolution'
      end

c------------------------------------------------------------------
      integer function nancheck(a)
c------------------------------------------------------------------
        if(.not.(a.gt.0..or.a.le.0.)
     .    .or.a.gt.1e35 .or. a.lt.-1e35    )then
          nancheck=1
        else
          nancheck=0
        endif
        end

c------------------------------------------------------------------
      subroutine checkpsar4
c------------------------------------------------------------------
#include "aaa.h"
#include "tab.h"
c      parameter (myom=7)
      common /psar4/  fhgg(11,100,10,8),fhqg(11,100,10,8)
      common /psar4b/ fhgq(11,100,10,8),fhqq(11,100,8)
      common /psar4c/ fhss(11,100,10,8)
      do i=1,11
      do j=1,100
      do l=mnclpt,mxclpt
        do k=1,10
        a=fhss(i,j,k,l)
        if(.not.(a.gt.0..or.a.le.0.)
     .     .or.a.gt.1e35 .or. a.lt.-1e35    )then
         write(ifmt,*)' checkpsar4    fhss    ',i,j,k,a
        endif
        a=fhgg(i,j,k,l)
        if(.not.(a.gt.0..or.a.le.0.)
     .     .or.a.gt.1e35 .or. a.lt.-1e35    )then
         write(ifmt,*)' checkpsar4    fhgg    ',i,j,k,a
        endif
        a=fhqg(i,j,k,l)
        if(.not.(a.gt.0..or.a.le.0.)
     .     .or.a.gt.1e35 .or. a.lt.-1e35    )then
           write(ifmt,*)' checkpsar4    fhqg    ',i,j,k,l,a
        endif
        a=fhgq(i,j,k,l)
        if(.not.(a.gt.0..or.a.le.0.)
     .     .or.a.gt.1e35 .or. a.lt.-1e35    )then
           write(ifmt,*)' checkpsar4     fhgq   ',i,j,k,l,a
        endif
        enddo
        a=fhqq(i,j,l)
        if(.not.(a.gt.0..or.a.le.0.)
     .     .or.a.gt.1e35 .or. a.lt.-1e35    )then
          write(ifmt,*)' checkpsar4     fhqq   ',i,j,l,a
        endif
      enddo
      enddo
      enddo
      end

c----------------------------------------------------------------------
      subroutine xAlphaSx !copied from xpr
c----------------------------------------------------------------------
#include "aaa.h"
#include "sem.h"
      parameter(nbib=1000)

      q1=sqrt(qcdlam)
      q2=100000

c************************full-red = Znorm *****************************

      write(ifhi,'(a)')       'openhisto name pssalf'
      write(ifhi,'(a)')       'htyp lru xmod log ymod lin'
      write(ifhi,*)       'xrange ',q1,q2,' yrange 0.05 0.7'
      write(ifhi,'(a)')    'text 0 0 "xaxis  Q (GeV)"'
      write(ifhi,'(a)')    'text 0 0 "yaxis  [a]?s!(Q) "'
      write(ifhi,'(a)')       'array 2'
      do k=0,nbib
        q=q1
        if (k.ne.0) q=q*(q2/q1)**(float(k)/float(nbib))
        qq=q*q
        y=pssalf(qq/qcdlam)*2.*pi
        write(ifhi,*)q,y
      enddo
      write(ifhi,'(a)')    '  endarray'
      write(ifhi,'(a)')    'closehisto plot 0-'
      end     

c-----------------------------------------------------------------------
      subroutine ffvalues()
c-----------------------------------------------------------------------
c parton distributions after evolution from q2soft to sc :
c     g1,g2 ...... gluons  
c     uv1,uv2 .... u valence quarks (nonsinglet evolution)
c     dv1,dv2 .... d valence quarks (nonsinglet evolution)
c     sea1,sea2 .. sea quarks (per flavor)
c "1" refers to the projectile and "2" to the target side
c-----------------------------------------------------------------------
#include "sem.h"
#include "aaa.h"
      stop'ERROR 06092019c dont use this routine any more' 
      end
      
c-----------------------------------------------------------------------
c      function ffsigSNGL(sis,klas,qt,x1,x2,ihq)  !removed in 3445n 
c                                                (recover if needed)
c-----------------------------------------------------------------------

