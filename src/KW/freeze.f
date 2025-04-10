C
C  This file is part of EPOS4
C  Copyright (C) 2022 research institutions and authors (See CREDITS file)
C  This file is distributed under the terms of the GNU General Public License version 3 or later
C  (See COPYING file for the text of the licence)
C

#if !__BS__ && !__TP__

c---------------------------------------------------------------------
c       subroutine sofo2x(nfr)  !spherio stuff removed in version 3244
c---------------------------------------------------------------------


c##############################################################################
      subroutine FroS(iflag)
c##############################################################################

c------------------------------------------------------------------------------
c     Freeze Out (FO) surface determination via call FoSur
c       and plots
c------------------------------------------------------------------------------
c iflag  1 ... call xxHo...    to plot fluid properties
c              call FoSur      to determine FO surface
c              call xxHoFo...  to plot FO properties         
c        2 ... call xxHo...Eta to plot eta dependence of FO properties
c        3 ... destroy hydro tables      
c------------------------------------------------------------------------------

#include "aaa.h"
#include "ho.h"
      common/cetaplot/etaplot(20),ietaplot(20),mxetaplot
      common/jcen/jcentr,jmxcentr
      data etaplot / 20*0. /
      fc=dtauhy/0.25*fdtau
      modu=5

      call getHo(idu,j,kdu)

      if(j.ne.2.and.ntauhy.eq.0)then
      if(iflag.eq.2)goto 3
      return
      endif

      GOTO (1,2,3,4,5),iflag

  1   continue

      if(ifmt.ne.6)open(ifmt,file=fnmt(1:nfnmt),access='append')
      mxetaplot=4  !8
      etaplot(1)=0
      etaplot(2)=2
      etaplot(3)=4
      etaplot(4)=5
      etaplot(5)=0
      etaplot(6)=-2
      etaplot(7)=-4
      etaplot(8)=-5
      call DefineEtaPlot

      if(iHyEpsilon.eq.1)call xxHoEpsilon(modu)
      if(iHyEntropy.eq.1)call xxHoEntropy(modu)
      if(iHyTemperature.eq.1)call xxHoTemperature(modu)
      if(iHyRadVelocity.eq.1)call xxHoRadVelocity(modu)
      if(iHyLongVelocity.eq.1)call xxHoLongVelocity(modu)
      if(iHyAverages.eq.1)call xxHoAverages
      if(iHyBaryon.eq.1)call xxHoBaryon(modu)

      if(j.ne.2)call FoSur

      do i=1,1   !2    !!!!!!!!!!!!!!!
      if(iHyFoVol.eq.1)call xxHoFoVol(i)
      if(iHyFoRadius.eq.1)call xxHoFoRadius(i)
      if(iHyFoEpsilon.eq.1)call xxHoFoEpsilon(i)
      if(iHyFoRadVelocity.eq.1)call xxHoFoRadVelocity(i)
      if(iHyFoTangVelocity.eq.1)call xxHoFoTangVelocity(i)
      if(iHyFoBarmu.eq.1)call xxHoFoBarmu(i)
      enddo

      do i=1,iHyEmpty
        call xEmpty
      enddo

      call clop(3)
      return

   2  continue

      if(ifmt.ne.6)open(ifmt,file=fnmt(1:nfnmt),access='append')
      if(iHyEpsilonEta.ge.1)call xxHoEpsilonEta(modu)
      if(iHyBaryonEta.ge.1)call xxHoBaryonEta(modu)
      if(iHyBarmuEta.eq.1)call xxHoBarmuEta(modu)
      if(iHyEntropyEta.eq.1)call xxHoEntropyEta(modu)
      call clop(3)

   3  continue

      if(ntauhy.eq.0)return
      call memo(1,'destroy hydro tables;')
      call destroyvelio()
      call destroybario()
      call destroyepsio()
      call destroyemuio()
      call memo(2,';')
      return

   4  continue

      if(ifmt.ne.6)open(ifmt,file=fnmt(1:nfnmt),access='append')
      ii=0
      if(iHyEpsilonEta.ge.1)then
      ii=ii+1
      do n=1,jmxcentr
       jcentr=n
       call xxHoEpsilonEta(modu)
      enddo
      endif
      if(iHyBaryonEta.ge.1)then
      ii=ii+1
      do n=1,jmxcentr
       jcentr=n
       call xxHoBaryonEta(modu)
      enddo
      endif
      if(iHyBarmuEta.eq.1)then
      ii=ii+1
      do n=1,jmxcentr
       jcentr=n
       call xxHoBarmuEta(modu)
      enddo
      endif
      if(iHyEntropyEta.eq.1)then
      ii=ii+1
      do n=1,jmxcentr
       jcentr=n
       call xxHoEntropyEta(modu)
      enddo
      endif
      if(ii.gt.0)call xxHoEtc
      call clop(3)
      return

   5  continue

      call histo2BeginFigure()
      if(iHoEpsilon.eq.1)call ooHoEpsilon(1,6, 2*fc)
      if(iHoRadVel.eq.1)call ooHoRadVel(1,6, 2*fc)
      if(iHoTanVel.eq.1)call ooHoTanVel(1,3, 2*fc)
      if(iHoEpsilonEtaY.eq.1)call ooHoEpsilonEtaY(1,2, 4*fc)
      !if(iHoEpsilonEtaY.eq.1)call ooHoVyEtaY(1,2, 4*fc)
      !if(iHoEpsilonEtaY.eq.1)call ooHoVzEtaY(1,2, 4*fc)
      if(iHoEpsilonEtas8.eq.1)call ooHoEpsilonEtas(8,1,2, 4*fc)
      if(iHoEpsilonEtas6.eq.1)call ooHoEpsilonEtas(6,1,2, 4*fc)
      if(iHoTemperatureEtas6.eq.1)
     .                 call ooHoTemperatureEtas(6, 15,15, 0.5)
      if(iHoEpsilonTau6.ge.1)then
      ii=iHoEpsilonTau6/10
      jj=mod(iHoEpsilonTau6,10)
      call ooHoEpsilonTau(6, 1,6*ii, fc*(10.+jj)/10.)
      endif
      if(iHoTangVelTau6.eq.1)call ooHoTangVelTau(6, 1,6, 2*fc)
      if(iHoRadVelTau6.ge.1)then
      ii=iHoRadVelTau6/10
      jj=mod(iHoRadVelTau6,10)
      call ooHoRadVelTau(6, 1,6*ii, fc*(10.+jj)/10.)
      endif

      end

c##############################################################################
      subroutine FoSur
c##############################################################################

c------------------------------------------------------------------------------
c  FO surface determination
c------------------------------------------------------------------------------

#include "aaa.h"
#include "ho.h"
#include "so.h"

      common/jcen/jcentr,jmxcentr
      real v1(3),v2(3)
      real dr2dtau(mxsurf, netahxx,ntauhxx,nphihxx)
     .    ,dr2dphi(mxsurf, netahxx,ntauhxx,nphihxx)
     .    ,dr2deta(mxsurf, netahxx,ntauhxx,nphihxx)
      data ncntfo/0/
      save ncntfo
      data ncntfo2/0/
      save ncntfo2

      integer getAccumJerr

      if(kfrout.eq.1)then
        write(ifmt,'(a,f7.4,a)')
     .  'computing freeze out surface for Tfo=',tfrout,' ...    '
      elseif(kfrout.eq.2.or.kfrout.eq.12.or.kfrout.eq.22)then
        write(ifmt,'(a,f7.4,a)')
     .  'computing freeze out surface for Efo=',efrout,' ...    '
      elseif(kfrout.eq.3)then
        write(ifmt,'(a,a)')
     .  'computing cone freeze out surface (Test 1)',' ...    '
      elseif(kfrout.eq.4)then
        write(ifmt,'(a,a)')
     .  'computing cone freeze out surface (Test 2)',' ...    '
      endif

      do neta=1,nzhy
      eta=etahy(neta)

      do nphi=1,nphihy
       phi=phihy(nphi)
       do ntau=1,ntauhy
        tau=getHydynTauhy(ntau)
        do nsurf=1,maxsurf
         radaa(nsurf,  neta,ntau,nphi)=0.
         velaa(nsurf,1,neta,ntau,nphi)=0.
         velaa(nsurf,2,neta,ntau,nphi)=0.
         velaa(nsurf,3,neta,ntau,nphi)=0.
         epsaa(nsurf,  neta,ntau,nphi)=0.
         baraa(nsurf,1,neta,ntau,nphi)=0.
         baraa(nsurf,2,neta,ntau,nphi)=0.
         baraa(nsurf,3,neta,ntau,nphi)=0.
         do i=1,4
         suraa(nsurf,i,neta,ntau,nphi)=0.
         enddo 
         do ivi=1,10
          emuaa(nsurf,ivi,neta,ntau,nphi)=0.
         enddo 
         dr2dphi(nsurf,neta,ntau,nphi)=0.
         dr2dtau(nsurf,neta,ntau,nphi)=0.
         dr2deta(nsurf,neta,ntau,nphi)=0.
        enddo
        nsurf=0
        dq12=1.
        do nrad=nradhy,2,-1
          r2=radhy(nrad)
          call epsioget(1,neta,ntau,nphi,nrad, e2 )
          call barioget(1,neta,ntau,nphi,nrad, a2 )
          call barioget(2,neta,ntau,nphi,nrad, b2 )
          call barioget(3,neta,ntau,nphi,nrad, c2 )
          t2=PiEos(5,e2,a2*iochem,b2*iochem,c2*iochem)
          r1=radhy(nrad-1)
          call epsioget(1,neta,ntau,nphi,nrad-1, e1 )
          call barioget(1,neta,ntau,nphi,nrad-1, a1 )
          call barioget(2,neta,ntau,nphi,nrad-1, b1 )
          call barioget(3,neta,ntau,nphi,nrad-1, c1 )
          t1=PiEos(5,e1,a1*iochem,b1*iochem,c1*iochem)
          dr=r2-r1
          !write(*,'(a,i4,5f7.3)')'+++++',nrad,e2,a2,b2,c2,t2
          if(kfrout.eq.1)then
            q1=t1
            q2=t2
            qfo=tfrout
          elseif(kfrout.eq.2.or.kfrout.eq.12.or.kfrout.eq.22)then
            q1=e1
            q2=e2
            qfo=efrout
          elseif(kfrout.eq.3)then !cone surface, for testing
            ntaumax=12 !arbitrary choice
            rtau=radhy(nradhy) *
     .      ( 1 - (tau-getHydynTauhy(1))/
     .           (getHydynTauhy(ntaumax)-getHydynTauhy(1)) )  
            rtau=max(0.,rtau)
            q1=-r1 
            q2=-r2
            qfo=-rtau
          elseif(kfrout.eq.4)then !cone surface, eta dependent, for testing
            ntaumax=12 !arbitrary choice
            etamax=6 !arbitrary choice
            etaa=abs(eta)
            rmaxx=radhy(nradhy) 
            rmax=rmaxx * ( 1 - etaa/etamax )  
            rmax=max(0.,rmax)
            rtau=rmax -
     .        rmaxx*(tau-getHydynTauhy(1))/
     .           (getHydynTauhy(ntaumax)-getHydynTauhy(1))   
            rtau=max(0.,rtau)
            q1=-r1 
            q2=-r2
            qfo=-rtau
          else
            stop'####### ERROR 22032018 #######'
          endif 
          ihit=0 
          if(((q1.lt.qfo.and.q2.gt.qfo.or.q1.gt.qfo.and.q2.lt.qfo)
     .    .or.q2.eq.qfo.and.(q1-q2)*dq12.gt.0.
     .    .or.q1.eq.qfo.and.nrad.eq.1
     .    .or.nsurf.eq.0.and.q2.gt.qfo)
     .    .and..not.(nsurf.eq.0.and.q1.le.qfo.and.q2.gt.qfo))ihit=1
          if(ihit.eq.1)then
            if(nsurf.ge.2)then
              if(abs(r2-radaa(nsurf,neta,ntau,nphi)).lt.r2-r1)then
                !print*,'remove 2 close layers at',neta,ntau,nphi 
                radaa(nsurf,neta,ntau,nphi)=0.
                nsurf=nsurf-1
                ihit=0
              endif
            endif
          endif
          if(ihit.eq.1)then
            if(nsurf+1.gt.mxsurf)then
               if(ish.ge.3)
     .         write(ifmt,*)'WARNING FoSur: mxsurf too small ',nsurf+1
               ihit=0
            endif
          endif
          if(ihit.eq.1)then !~~~~~ hit = 1
          nsurf=nsurf+1
          if( nsurf.le.maxsurf)then
                     !***** error print out,
                     !***** maybe removed when code is stable *****
                     if(mod(nsurf,2).eq.0.and.q1.gt.q2
     .                .or.(mod(nsurf,2).eq.1.and.q1.lt.q2
     .               .and..not.(nsurf.eq.1.and.q2.gt.qfo)) )then
                     do n=nrad-1,min(nrad+5,nradhy)
                     call epsioget(1,neta,ntau,nphi,n, e22 )
                     call barioget(1,neta,ntau,nphi,n, a22 )
                     call barioget(2,neta,ntau,nphi,n, b22 )
                     call barioget(3,neta,ntau,nphi,n, c22 )
                     t22=PiEos(5,e22,a22,b22,c22)
                     write(ifmt,*)n,e22,t22
     .                ,q1.lt.qfo.and.q2.gt.qfo  ,q1.gt.qfo.and.q2.lt.qfo
     .               ,q2.eq.qfo   ,(q1-q2)*dq12.gt.0.  ,(q1-q2)*dq12
     .                  ,q1-q2, dq12
                     enddo
                     endif
                     !***** end error output *****
           if(mod(nsurf,2).eq.1.and.q1.lt.q2
     .     .and..not.(nsurf.eq.1.and.q2.gt.qfo)  )
     .      stop'ERROR FoSur: q1 < q2  => check!!! '
           if(mod(nsurf,2).eq.0.and.q1.gt.q2)
     .      stop'ERROR FoSur: q1 > q2  => check!!! '
           call velioget(1,neta,ntau,nphi,nrad-1, v1(1) )
           call velioget(2,neta,ntau,nphi,nrad-1, v1(2) )
           call velioget(3,neta,ntau,nphi,nrad-1, v1(3) )
           call velioget(1,neta,ntau,nphi,nrad, v2(1) )
           call velioget(2,neta,ntau,nphi,nrad, v2(2) )
           call velioget(3,neta,ntau,nphi,nrad, v2(3) )
           vr1=v1(1)*cos(phi)+v1(2)*sin(phi)
           vr2=v2(1)*cos(phi)+v2(2)*sin(phi)
           call setAccumJerr(21,getAccumJerr(21)+1)
           if(nsurf.eq.1.and.q2.gt.qfo)then
             if(ish.ge.3)
     .       write(ifmt,*)'WARNNG FoSur:  nrad limit reached',r2
             call SetAccumJerr(22,getAccumJerr(22)+1)  !FoSur: nrad limit reached
             f=1
           elseif(q1.ne.q2)then
             f=(q1-qfo)/(q1-q2)
           else
             f=0.5
           endif
           r=r1+f*dr
           radaa(nsurf,neta,ntau,nphi)=r
           !if(nsurf.eq.4)then
           !  do i7=1,nradhy 
           !  call epsioget(1,neta,ntau,nphi,i7, er )
           !  write(*,'(3i5,4x,f7.3)')neta,ntau,nphi,er / qfo
           !  enddo
           !  print*,' '
           !endif
           w1=v1(1)
           w2=v2(1)
           v=w1+f*(w2-w1)
           if(v.gt.1.)stop'\n\n STOP in FoSur: vx > 1.  '
           vx=v
           velaa(nsurf,1,neta,ntau,nphi)=vx !<--------------------
           w1=v1(2)
           w2=v2(2)
           v=w1+f*(w2-w1)
           if(v.gt.1.)stop'\n\n STOP in FoSur: vy > 1.  '
           vy=v
           velaa(nsurf,2,neta,ntau,nphi)=vy !<--------------------
           w1=v1(3)
           w2=v2(3)
           v=w1+f*(w2-w1)
           if(v.gt.1.)stop'\n\n STOP in FoSur: vz > 1.  '
           vz=v
           velaa(nsurf,3,neta,ntau,nphi)=vz !<--------------------
           baraa(nsurf,1,neta,ntau,nphi)=a1+f*(a2-a1)
           baraa(nsurf,2,neta,ntau,nphi)=b1+f*(b2-b1)
           baraa(nsurf,3,neta,ntau,nphi)=c1+f*(c2-c1)
           epsaa(nsurf,neta,ntau,nphi)=e1+f*(e2-e1)
           do ivi=1,10
             call emuioget(ivi,neta,ntau,nphi,nrad-1, w1 )
             call emuioget(ivi,neta,ntau,nphi,nrad, w2 )
             w=w1+f*(w2-w1)
             emuaa(nsurf,ivi,neta,ntau,nphi)=w
           enddo
          endif
          endif !~~~~~ hit = 1
          dq12=q1-q2
          if(ntau.eq.1)then
            do ivi=1,10
              call emuioget(ivi,neta,ntau,nphi,nrad, w )
              emuzz(ivi,neta,nrad,nphi)=w
            enddo
            do kb=1,3
              call barioget(kb,neta,ntau,nphi,nrad, w )
              barzz(kb,neta,nrad,nphi)=w
            enddo
            do kv=1,3
              call velioget(kv,neta,ntau,nphi,nrad, w )
              velzz(kv,neta,nrad,nphi)=w
            enddo
          endif 
        enddo !nrad
       enddo !ntau
      enddo !nphi

      enddo !neta

      write(ifmt,'(a)')'freeze out surface computation done'
      write(ifmt,'(a,$)')'making aa tables ...'

      do nsurf=1,maxsurf

      do neta=1,nzhy
        do nphi=1,nphihy
            n=ntauhy
            do while(radaa(nsurf,neta,n,nphi).eq.0.0.and.n.gt.2)
              n=n-1
            enddo
            n=n+1
            n=min(n,ntauhy)
            ntauhec(neta,nphi)=n
        enddo
      enddo

      dphi=phihy(2)-phihy(1)
      dtau=getHydynTauhy(2)-getHydynTauhy(1)
      deta=etahy(2)-etahy(1)
      do neta=1,nzhy
        nem=neta-1
        nem=max(nem,1)
        nep=neta+1
        nep=min(nep,nzhy)
        do nphi=1,nphihy
          npp=nphi+1
          npp=min(npp,nphihy)
          npm=nphi-1
          npm=max(npm,1)
          do ntau=1,ntauhec(neta,nphi)
            ntm=ntau-1
            ntm=max(ntm,1)
            ntp=ntau+1
            ntp=min(ntp,ntauhec(neta,nphi))
            if( radaa(nsurf,neta,ntau,nphi) .gt. 0. )then
              dr2dphi(nsurf,neta,ntau,nphi)
     .           =(radaa(nsurf,neta,ntau,npp )
     .           -radaa(nsurf,neta,ntau,npm )) / ((npp-npm)*dphi)
              dr2dtau(nsurf,neta,ntau,nphi)
     .           =(radaa(nsurf,neta,ntp ,nphi)
     .           -radaa(nsurf,neta,ntm ,nphi)) / ((ntp-ntm)*dtau)
              dr2deta(nsurf,neta,ntau,nphi)
     .           =(radaa(nsurf,nep ,ntau,nphi)
     .           -radaa(nsurf,nem ,ntau,nphi)) / ((nep-nem)*deta)
            endif
          enddo
        enddo
      enddo

      do neta=1,nzhy
       do nphi=1,nphihy
        phi=phihy(nphi)
        do ntau=1,ntauhec(neta,nphi)
         if( radaa(nsurf,neta,ntau,nphi) .gt. 0. )then
          tau=getHydynTauhy(ntau)
          rad=radaa(nsurf,neta,ntau,nphi)
          vx=velaa(nsurf,1,neta,ntau,nphi)
          vy=velaa(nsurf,2,neta,ntau,nphi)
          vz=velaa(nsurf,3,neta,ntau,nphi)
          vv=sqrt(vx**2+vy**2+vz**2)
          gm=1./sqrt((1-vv)*(1+vv))
          suraa(nsurf,4,neta,ntau,nphi)                    !0
     .     = -dr2dtau(nsurf,neta,ntau,nphi)*rad*tau
          suraa(nsurf,1,neta,ntau,nphi)                    !1
     .     = rad*tau*cos(phi)+dr2dphi(nsurf,neta,ntau,nphi)*tau*sin(phi)
          suraa(nsurf,2,neta,ntau,nphi)                    !2
     .     = rad*tau*sin(phi)-dr2dphi(nsurf,neta,ntau,nphi)*tau*cos(phi)
          suraa(nsurf,3,neta,ntau,nphi)                    !3
     .     = -dr2deta(nsurf,neta,ntau,nphi)*rad  
          dVs=  suraa(nsurf,4,neta,ntau,nphi)  *gm               !0
     .        + suraa(nsurf,1,neta,ntau,nphi)  *gm*vx            !1 
     .        + suraa(nsurf,2,neta,ntau,nphi)  *gm*vy            !2
     .        + suraa(nsurf,3,neta,ntau,nphi)  *gm*vz            !3
          vlmaa(nsurf,neta,ntau,nphi)=abs(dVs)
         endif
        enddo !tau
       enddo !phi
      enddo !neta

      enddo ! nsurf
     
      write(ifmt,'(a)')'  done'

      end
  

c######################################################################################
      subroutine EpoMiF 
c######################################################################################

c--------------------------------------------------------------------------------------
c    Micro-canonical FO via FO hyper-surface  
c--------------------------------------------------------------------------------------

c--------------------------------------------------------------------------------------
c kfrout <= 2 ... Nothing (returns immediately)
c         21 .... Using T for FO surface
c         22 .... Using epsilon for FO surface
c--------------------------------------------------------------------------------------

#include "aaa.h"
#include "ho.h"
      common/cranphi/ranphi
c      common/cee1ico/ee1ico,eistico,ee1hll
      parameter (mxeflo=21)
      real eflosa(mxeflo),eflose(mxeflo)
      real eflosu(0:mxsurf),eflos(0:mxsurf),pf(4),s(4),t(10)
      real amsu(0:mxsurf),amsui(0:mxsurf),eflosum
      real dMeta(netahxx),dFsum(3),dFlav(3,netahxx),nspleta(netahxx)
      real FoEta(1000),FoTau(1000)  
      common/citsy/itsy(4,4)
      parameter (ncluxx=netahxx*10)
      integer nfuse(2,ncluxx)
      real dMfuse(ncluxx),dFfuse(3,ncluxx),dGfuse(3,ncluxx)
      integer idrfuse(ncluxx),iflav(ncluxx),jflav(ncluxx)
      real u(4),dFcms(3),dF(3),p(4),psu(4),psa(4),psi(4),psum(4)
      real PcluNF(4),psumNF(4),uu(4)
      parameter (netaclu=30)
      parameter (mradtau=max(ntauhxx,nradhxx))
      real weta(netaclu),wsu(netaclu,mxsurf+1)
      real wradtau(netaclu,mxsurf+1,mradtau)
      real wphi(netaclu,mxsurf+1,mradtau,nphihxx)
      integer jc(nflav,2),ic(2)
      integer ncountmic
      real v5(5),ve3(3)
      real dMetaSum(netahxx)
      data dMetaSum / netahxx*0 /
      data ncountmic/0/ ncount/0/
      save ncountmic,ncount

      call utpri('epomif',ish,ishini,4)

      if(kfrout.le.2)return

      if(ntauhy.eq.0)goto 9999

      call checkGrandCanon(0,0,0.,ve3) !initialize to zero

      !early hadronization at large eta, later one at small eta -> eta-tau correlation
      !therfore enough to consider dFlav(i,neta) and not in addition ntau dependence
      do i=1,3
      do neta=1,netahxx
        dFlav(i,neta)=0   
      enddo
      enddo

      phinull=phievt+ranphi

c create dmass_

      call memo(1,'create dmass objects ;')
      call dmass_c0create(mxsurf)
      call dmass_c1create(mxsurf,nzhy)
      call dmass_c2create(mxsurf,nzhy,ntauhy)
      call dmass_c3create(mxsurf,nzhy,ntauhy,nphihy)
      call dmass_i0create(  1   )
      call dmass_i1create(  1   ,nzhy)
      call dmass_i2create(  1   ,nzhy,nradhy)
      call dmass_i3create(  1   ,nzhy,nradhy,nphihy)
      call memo(2,';')
      
c compute dmass_ level 3 and dFlav

      dphi=phihy(2)-phihy(1)
      dtau=getHydynTauhy(2)-getHydynTauhy(1)
      deta=etahy(2)-etahy(1)
      drad=radhy(2)-radhy(1)  

      do i=1,netahxx
        dMeta(i)=0
      enddo
      eflosum=0
      do i=1,mxeflo
        eflose(i)=0
      enddo

      do nsurf=0,maxsurf !~~~~~sum over surface layers~~~~~~

      ! nsurf = 0 : 
      !
      !   surface "zero" : const tau  (tau initial) outside cylinder (for layer 1)

      ! nsurf > 0 : 
      !
      !   surface "nsurf" : cylinder r = r(phi,eta,tau), layer number nsurf

      eflosu(nsurf)=0
      amsu(nsurf)=0
      amsui(nsurf)=0

      do ntau=1,ntauhy !~~~~~~ tau sum ~~~~~~~~~~~
      tau=getHydynTauhy(ntau)

      eflos(nsurf)=0
      do i=1,mxeflo
        eflosa(i)=0
      enddo

      do nrad=2,nradhy !~~~~~~ rad sum ~~~~~~~~~~ (needed for surface "zero")
      rd=radhy(nrad)

      do neta=1,nzhy  !~~~~~~ eta sum ~~~~~~~~~~~~
      nzz=neta-int(nzhy/2.)+int(mxeflo/2.)
      eta=etahy(neta)

      do nphi=2,nphihy !~~~~~~ phi sum ~~~~~~~~~~
      phi=phihy(nphi)

      if(iSkipSurfForCheck(nsurf,neta,ntau,nphi).eq.1)goto 888

      if(nsurf.eq.0.and.ntau.eq.1 .or. nsurf.gt.0.and.nrad.eq.2)then!~~~~ surface type ~~~~

          !~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
          !  quantities in Minkowski space, Bjorken frame
          !    identical to tilde quantities in Milne
          !~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~  

          ! surface "zero"
          !
          !   emuzz(k,neta,nrad,nphi)         ! T^k (10-vector) 
          !   barzz(k,neta,nrad,nphi)         ! B-, Q-, C-density (3-vector) 
          !   velzz(k,neta,nrad,nphi)         ! x-, y-, z-velocity (3-vector) 
          !
          ! surface "nsurf"
          !
          !   velaa(nsurf,1,neta,ntau,nphi)   ! x velocity
          !   velaa(nsurf,2,neta,ntau,nphi)   ! y velocity
          !   velaa(nsurf,3,neta,ntau,nphi)   ! z velocity
          !   radaa(nsurf,  neta,ntau,nphi)   ! radius
          !   epsaa(nsurf,  neta,ntau,nphi)   ! energy density
          !   baraa(nsurf,1,neta,ntau,nphi)   ! B- 
          !   baraa(nsurf,2,neta,ntau,nphi)   ! Q- densities
          !   baraa(nsurf,3,neta,ntau,nphi)   ! S-   
          !   suraa(nsurf,k,neta,ntau,nphi)   ! dSigma_k (4-vector) 
          !   emuaa(nsurf,k,neta,ntau,nphi)   ! T^k (10-vector) 
          !
          ! mass distribution 
          !
          !  call dmass_i0get(  1  )                   !surface zero 
          !  call dmass_i1get(  1  ,neta)              ! 1-nzhy
          !  call dmass_i2get(  1  ,neta,nrad)         !  1-nradhy
          !  call dmass_i3get(  1  ,neta,nrad,nphi)    !   2-nphihy
          !  call dmass_c0get(nsurf)                   !1-maxsurf
          !  call dmass_c1get(nsurf,neta)              ! 1-nzhy
          !  call dmass_c2get(nsurf,neta,ntau)         !  1-ntauhy
          !  call dmass_c3get(nsurf,neta,ntau,nphi)    !   2-nphihy

          rad=radaa(1,neta,ntau,nphi)   

          if(nsurf.eq.0)then
            do k=1,4
              s(k)=0
            enddo
            if(rd.gt.rad)s(4)=tau*rd*dphi*drad*deta
            do k=1,10
              t(k)=emuzz(k,neta,nrad,nphi) 
            enddo 
            wB=barzz(1,neta,nrad,nphi)  
            wQ=barzz(2,neta,nrad,nphi)  
            wS=barzz(3,neta,nrad,nphi)  
          else!nsurf.gt.0
            isig=mod(nsurf,2)*2-1
            do k=1,4
              s(k)=suraa(nsurf,k,neta,ntau,nphi)*dphi*dtau*deta*isig
            enddo 
            do k=1,10
              t(k)=emuaa(nsurf,k,neta,ntau,nphi) 
            enddo 
            wB=baraa(nsurf,1,neta,ntau,nphi)   ! B- 
            wQ=baraa(nsurf,2,neta,ntau,nphi)   ! Q- densities
            wS=baraa(nsurf,3,neta,ntau,nphi)   ! S-   
          endif
          dFcms(1) = wB+wQ         !u
          dFcms(2) = 2*wB-wQ+wS    !d
          dFcms(3) =  - wS         !s

          ! energy momentum flow vector per surface element 

          do i=1,4
            pf(i)=0
            do j=1,4
              pf(i) = pf(i) + t ( itsy(i,j) ) * s(j)
            enddo
          enddo

          ! energy flow in Minkowski, Lab 

          eflo=(pf(3)*sinh(eta)+pf(4)*cosh(eta))
          !if(nsurf.eq.4.and.eflo.ne.0.)then
          !  print*,'TEST',ntau,neta,nphi
          !endif
          eflos(nsurf)=eflos(nsurf)+eflo
          if(nzz.ge.1.and.nzz.le.mxeflo)then
            eflosa(nzz)=eflosa(nzz)+eflo
            eflose(nzz)=eflose(nzz)+eflo
          endif
          eflosu(nsurf)=eflosu(nsurf)+eflo

          ! check masses

          dM2=pf(4)**2-pf(1)**2-pf(2)**2-pf(3)**2
          if(dM2.gt.0.)then
            dM=sqrt(dM2)
            if(nsurf.eq.0)then !nsurf = 0
              call dmass_i3set(  1  ,neta,nrad,nphi,dM)
            else             !nsurf = 1 , 2
              call dmass_c3set(nsurf,neta, ntau ,nphi,dM)
            endif  
            amsu(nsurf)=amsu(nsurf)+dM
            dMeta(neta)= dMeta(neta)+dM
            dMetaSum(neta)= dMetaSum(neta)+dM
          elseif(dM2.lt.0.)then
            dM=sqrt(-dM2)
            amsui(nsurf)=amsui(nsurf)+dM
          endif  
          if(ish.ge.8)then
            if(nphi.eq.2)write(ifch,'(4i5,$)')nsurf,ntau,nrad,neta 
            write(ifch,'(f5.0,$)')dM2*10000
            if(nphi.eq.nphihy)write(ifch,*)' '
          endif  

          ! four-velocity from pf

          vcutfreeze=0.85
          do i=1,4
            uu(i)=0
            if(i.eq.4)uu(i)=1
            if(dM.ne.0..and.pf(4).ne.0.)then
              v=sqrt(pf(1)**2+pf(2)**2+pf(3)**2)/pf(4)
              if(v.lt.vcutfreeze)then !to avoid rare cases with v->1 giving artificial tails in pt at low E 
                uu(i)=pf(i)/dM   !four velocity of mass element
              endif
            endif
          enddo  

          !~~~~~~~~~~
          !if(nsurf.eq.1)then
          !if(neta.eq.15.and.ntau.eq.2.and.nphi.eq.42)then
          !!if(dM2.gt.0..and.dM.lt.0.01)then
          !v=sqrt(uu(1)**2+uu(2)**2+uu(3)**2)/uu(4)
          !pfp=sqrt(pf(1)**2+pf(2)**2+pf(3)**2)
          !write(ifmt,*)'CHECKfreeze ',v,dM2,dM, pf(4),pfp
          !endif
          !endif
          !~~~~~~~~~~

          ! four-velocity from surface

          if(nsurf.eq.0)then
            do k=1,3
              u(k)=velzz(k,neta,nrad,nphi)   
            enddo
          else!nsurf.gt.0
            do k=1,3
              u(k)=velaa(nsurf,k,neta,ntau,nphi)
            enddo
          endif
          gmx=max(1e-8, 1.-u(1)**2-u(2)**2-u(3)**2)
          gm=1./sqrt(gmx)
          u(4)=1
          do k=1,4
            u(k)=u(k)*gm
          enddo
          if(ish.ge.8)then
            if(dM2.ne.0.)then
            write(ifch,*)nsurf,ntau,nrad,neta,nphi,'     '
     .     ,dM2,dM,uu(3),u(3)
            endif
          endif

          ! choice of four-velocity 

          iii=2  !1 = via surface,  2 = via pf (correct)
          if(iii.eq.2)then
            do i=1,4
              u(i)=uu(i)
            enddo  
            gm=u(4)
            if(nsurf.eq.0)then
              do k=1,3
                velzz(k,neta,nrad,nphi)=u(k)/gm 
                ve3(k)=u(k)/gm 
              enddo
            else!nsurf.gt.0
              do k=1,3
                velaa(nsurf,k,neta,ntau,nphi)=u(k)/gm
                ve3(k)=u(k)/gm 
              enddo
            endif
          endif

          ! flavor flow per surface element 

          do i=1,3 !u,d,s
            dF(i)=0
            do j=1,4
              dF(i) = dF(i) + dFcms(i)*u(j) * s(j)
            enddo
            dFlav(i,neta)=dFlav(i,neta)+dF(i)
          enddo

          ! check

          imess=0
          if(nphi.eq.nphihy.and.neta.eq.nzhy)imess=max(ntau,nrad) 
          if(dM2.gt.0.)then
            dM=sqrt(dM2)
            call checkGrandCanon(1,imess,dM,ve3) !fill
          else
            call checkGrandCanon(1,imess,0.,ve3) !fill
          endif

      endif !surface type

 888  continue

      enddo !nphi

      enddo !neta

      enddo !nrad

      if(nsurf.eq.0.and.ntau.eq.1)then
      nen=log10(getIcoEe1ico())
      fac=10.**(4.-nen)
      if(ish.ge.3)write(ifmt,'(a,2f10.2,1x,21i6)')'Eflow0'
     .,eflosu(0),eflosu(0),(nint(eflosa(i)),i=1,mxeflo)   !eflosa(i)*fac
      endif

      if(ish.ge.3)then
        if(nsurf.eq.1.and.eflos(1).gt.0.)
     .  write(ifmt,'(a,2f10.2,1x,21i6)')'Eflow1'
     .  ,eflos(1),eflosu(1)+eflosu(0),(nint(eflosa(i)),i=1,mxeflo) 
      endif

      enddo !ntau

      eflosum=eflosum+eflosu(nsurf)
      esu=0
      do i=1,mxeflo
        esu=esu+eflose(i)
      enddo
      if(ish.ge.3)then
        if(nsurf.eq.1)
     .  write(ifmt,'(a,4x,f10.2,$)')'Eflow0+1',esu
        if(nsurf.eq.1)
     .  write(ifmt,'(5x,21i6)')(nint(eflose(i)),i=1,mxeflo) 
      endif  
      if(ish.ge.0)then
        if(nsurf.eq.1)
     .  write(ifmt,'(a,4x,2f10.2)')'Eflow0+1'
     .  ,eflosu(0)+eflosu(1),eflosum
        if(nsurf.gt.1)
     .  write(ifmt,'(a,i1,2f12.2)')'Eflow',nsurf
     .  ,eflosu(nsurf),eflosum
      endif
 
      enddo !nsurf

c print  

      mo=int(mxeflo/2.)
      if(ish.ge.4)then
        write(ifmt,'(a,5f8.1,$)')'real masses:',amsu
        write(ifmt,'(3x,a,5f8.1)')'imag masses:',amsui
      endif
      if(ish.ge.3)then
        write(ifmt,'(a,i3,$)')'mass eta distribut >=',1+nzhy/2-mo
        write(ifmt,'(3x,21f6.1)')(dMeta(i+nzhy/2-mo),i=1,mxeflo) 
        write(ifmt,'(a,3x,a,$)')'   mass*cosh(eta)    ','   '
        esu=0
        do i=1,mxeflo
          neta=i+nzhy/2-mo
          eta=etahy(neta)
          esu=esu+dMeta(neta)*cosh(eta)
          write(ifmt,'(i6,$)')nint(dMeta(neta)*cosh(eta))
        enddo
        write(ifmt,'(a)')' '
      endif  
      if(ish.ge.4)then
        write(ifmt,'(a,i3,$)')'flav eta distribut >=',1+nzhy/2-mo
        write(ifmt,'(3x,21f6.1)')
     .  (dFlav(1,n)+dFlav(2,n)+dFlav(3,n),n=1+nzhy/2-mo,1+nzhy/2+mo)
      endif
      if(ish.ge.3)then
      write(ifmt,'(a,f10.2,3x,a,f10.2)')'EfluidIni',getIcoEe1hll(),
     .'Sum mass*cosh(eta)',esu  
      endif

c compute dmass_ levels 1,2,3,     cumulative for levels 2,3, NOT for 1,     0 not used

      do nsurf=0,maxsurf 
      do neta=1,nzhy  
      do ntau=1,ntauhy 
      do nrad=2,nradhy 
      if(nsurf.eq.0.and.ntau.eq.1 .or. nsurf.gt.0.and.nrad.eq.2)then!~~~~ surface type ~~~~
        if(nsurf.eq.0)nn=nrad
        if(nsurf.eq.0)nnmin=2
        if(nsurf.ne.0)nn=ntau
        if(nsurf.ne.0)nnmin=1
        do nphi=2,nphihy 
          if(nsurf.eq.0)then
            call dmass_i3get(  1  ,neta,nn,nphi,dM)
            if(nphi.gt.2)call dmass_i3acc(  1  ,neta,nn,nphi,dM)
          else !(nsurf.ne.0)
            call dmass_c3get(nsurf,neta,nn,nphi,dM)
            if(nphi.gt.2)call dmass_c3acc(nsurf,neta,nn,nphi,dM)
          endif 
        enddo !nphi
        dM=0
        if(nsurf.eq.0)call dmass_i3get(  1  ,neta,nn,nphihy,dM)
        if(nsurf.ne.0)call dmass_c3get(nsurf,neta,nn,nphihy,dM)
        if(nsurf.eq.0)call dmass_i2set(  1  ,neta,nn,dM)
        if(nsurf.ne.0)call dmass_c2set(nsurf,neta,nn,dM)
        if(nn.gt.nnmin)then
          if(nsurf.eq.0)call dmass_i2acc(  1  ,neta,nn,dM)
          if(nsurf.ne.0)call dmass_c2acc(nsurf,neta,nn,dM)
        endif
      endif!~~~~ surface type ~~~~
      enddo !nrad
      enddo !ntau
      dM=0
      if(nsurf.eq.0)call dmass_i2get(  1  ,neta,nradhy,dM)
      if(nsurf.ne.0)call dmass_c2get(nsurf,neta,ntauhy,dM)
      !~~~test~~~
      !if(neta.ge.nzhy/2-2.and.neta.le.nzhy/2+2.and.nsurf.eq.1)then 
      !  write(ifmt,'(a,2i4,f9.2)')'nsurf neta mass:',nsurf,neta,dM
      !  do i=1,min(12,ntauhy)
      !    call dmass_c2get(nsurf,neta,i,dMt)
      !    write(ifmt,'(f6.2,$)')dMt
      !  enddo
      !  write(ifmt,'(a)')' '
      !endif 
      !~~~~~~
      if(nsurf.eq.0)call dmass_i1set(  1  ,neta,dM)
      if(nsurf.ne.0)call dmass_c1set(nsurf,neta,dM)
      enddo !neta
      enddo !nsurf

c create clusters (fused eta slices): compute dMfuse and dFfuse     

      amaxi=300
      amin=2
      neta=0
      neta1=0
      dMsum=0       
      do k=1,3
        dFsum(k)=0 
      enddo
      nclu=0
  10  neta=neta+1
      call dmass_i1get(  1  ,neta,dM)
      dMsum=dMsum+dM
      do nsurf=1,maxsurf
        call dmass_c1get(nsurf,neta,dM)
        dMsum=dMsum+dM
      enddo
      do k=1,3
        dFsum(k)=dFsum(k)+dFlav(k,neta)
      enddo
      if(dMsum.lt.0.01)then
        if(neta.lt.nzhy/2)then
          neta1=neta
          goto 10
        endif
      endif
      if(neta1+1.lt.nzhy/2.and.neta.lt.nzhy*3/4)then
        if(dMsum.lt.amin)goto 10
        if(neta-neta1.lt.nfrout.and.dMsum.lt.amaxi)goto 10
      endif
      nclu=nclu+1
      nfuse(1,nclu)=neta1+1
      nfuse(2,nclu)=neta
      dMfuse(nclu)=dMsum
      do k=1,3
        dFfuse(k,nclu)=dFsum(k) 
      enddo
      neta1=neta
      dMsum=0 
      do k=1,3
        dFsum(k)=0 
      enddo
      if(neta.lt.nzhy)goto 10

      nclu=nclu+1
  20  nclu=nclu-1
      dM=dMfuse(nclu)
      if(nclu.gt.1.and.dM.lt.0.01)then
        dMfuse(nclu-1)=dMfuse(nclu-1)+dM
        do k=1,3
          dFfuse(k,nclu-1)=dFfuse(k,nclu-1)+dFfuse(k,nclu)
        enddo
        goto 20 
      endif
      if(nclu.gt.1)then
        if(dM.lt.amin
     .    .or.(nfuse(2,nclu)-nfuse(1,nclu)+1.lt.nfrout
     .           .and.dM.lt.amaxi))then
          dMfuse(nclu-1)=dMfuse(nclu-1)+dM
          do k=1,3
            dFfuse(k,nclu-1)=dFfuse(k,nclu-1)+dFfuse(k,nclu)
          enddo
          nfuse(2,nclu-1)=nfuse(2,nclu)
          goto 20 
        endif
      endif

      n=0
  30  n=n+1
      if(n.eq.nclu)goto 35
      dM=dMfuse(n)
      if(dM.lt.amin
     ..or.(nfuse(2,n)-nfuse(1,n)+1.lt.nfrout.and.dM.lt.amaxi))then
        dMfuse(n)=dMfuse(n)+dMfuse(n+1)
        do k=1,3
          dFfuse(k,n)=dFfuse(k,n)+dFfuse(k,n+1)
        enddo
        nfuse(2,n)=nfuse(2,n+1)
        do nn=n+2,nclu
          dMfuse(nn-1)=dMfuse(nn)
          do k=1,3
            dFfuse(k,nn-1)=dFfuse(k,nn)
          enddo
          nfuse(1,nn-1)=nfuse(1,nn)
          nfuse(2,nn-1)=nfuse(2,nn)
        enddo
        nclu=nclu-1
        n=n-1
      endif 
      if(n.lt.nclu)goto 30
  35  continue

c split clusters in case of high masses

      do neta=1,netahxx
        nspleta(neta)=1
      enddo
      n=0
  40  n=n+1
      if(dMfuse(n).gt.amaxi)then
        nspl=dMfuse(n)/amaxi+1
        if(nclu+nspl.gt.ncluxx)then
          write(ifmt,'(a,i5)')'PROBLEM for cluster',n
          write(ifmt,*)'nclu=',nclu,'   nspl=',nspl,'   ncluxx=',ncluxx
          do m=1,nclu
            write(ifmt,'(a,3i5,4f8.2)')'cluster',m,nfuse(1,m)
     .     ,nfuse(2,m),dMfuse(m),(dFfuse(k,m),k=1,3)
          enddo
          write(ifmt,'(a)')'ERROR 11052019 EpoMiF split clusters'
          stop
        endif
        !create nspl spaces to place sub-clusters
        do m=nclu+nspl,n+nspl+1,-1
          dMfuse(m)=dMfuse(m-nspl)
          do k=1,3
            dFfuse(k,m)=dFfuse(k,m-nspl)
          enddo
          nfuse(1,m)=nfuse(1,m-nspl)
          nfuse(2,m)=nfuse(2,m-nspl)
        enddo
        nclu=nclu+nspl
        !split n to n,n+nspl
        do m=nfuse(1,n),nfuse(2,n)
          nspleta(m)=nspl+1
        enddo
        dMfuse(n)=dMfuse(n)/(nspl+1)
        do k=1,3
          dFfuse(k,n)=dFfuse(k,n)/(nspl+1)
        enddo
        do m=n+1,n+nspl
          dMfuse(m)=dMfuse(n)
          do k=1,3
            dFfuse(k,m)=dFfuse(k,n)
          enddo
          nfuse(1,m)=nfuse(1,n)
          nfuse(2,m)=nfuse(2,n)
        enddo 
        n=n+nspl
      endif
      if(n.lt.nclu)goto 40
      do m=1,nclu
        a=dMfuse(m)
        if(.not.(a.gt.0..or.a.le.0.))then !NaN catch
          write(ifmt,'(a)')'PROBLEM concerning dMfuse'
          write(ifmt,*)'dMfuse=',dMfuse(m),m,nfuse(1,m),nfuse(2,m)
          write(ifmt,'(a)')'ERROR 12052019b EpoMiF: dMfuse is NaN'
          stop
        endif
      enddo

c fuse cluster from unsuccessful  decay

      ncluE=0
  50  continue
      if(ncluE.gt.0)then
        n=ncluE
        if(n.eq.nclu)n=n-1
        !fuse n and n+1
        dMfuse(n)=dMfuse(n)+dMfuse(n+1)
        do k=1,3
          dFfuse(k,n)=dFfuse(k,n)+dFfuse(k,n+1)
        enddo
        nfuse(2,n)=nfuse(2,n+1)
        do nn=n+2,nclu
          dMfuse(nn-1)=dMfuse(nn)
          do k=1,3
            dFfuse(k,nn-1)=dFfuse(k,nn)
          enddo
          nfuse(1,nn-1)=nfuse(1,nn)
          nfuse(2,nn-1)=nfuse(2,nn)
        enddo
        nclu=nclu-1
      endif

c determine integer baryon number for clusters

      do n=1,nclu
        isum=0
        do k=1,3
          b=dFfuse(k,n)
          is=sign(1.,b)
          b=abs(b) 
          ib=int(b)
          if(rangen().lt.b-ib)ib=ib+1 
          dGfuse(k,n)=ib*is
          isum=isum+dGfuse(k,n)
        enddo
        is=sign(1,isum)
        ik=mod(abs(isum),3)
        if(ik.eq.1)ik=-1
        if(ik.eq.2)ik=1
        ik=ik*is
        if(ik.ne.0)then 
          k=min(3,1+int(rangen()*3))
          dGfuse(k,n)=dGfuse(k,n)+ik
        endif
        iflav(n)=isum
        jflav(n)=ik
      enddo

c compute cluster sum

      dMsum=0 
      do n=1,nclu
        dMsum=dMsum+dMfuse(n)
      enddo

c add cluster sum to cptl
    
      nptl=nptl+1
      idptl(nptl)=700000000
      pptl(1,nptl)=0
      pptl(2,nptl)=0
      pptl(3,nptl)=0
      pptl(4,nptl)=dMsum
      pptl(5,nptl)=dMsum
      ityptl(nptl)=60
      istptl(nptl)=11
      tivptl(2,nptl)=1.
      iorptl(nptl)=0
      jorptl(nptl)=0
      xorptl(1,nptl)=0
      xorptl(2,nptl)=0
      xorptl(3,nptl)=0
      xorptl(4,nptl)=0
      ifrptl(1,nptl)=0
      ifrptl(2,nptl)=0
      nptla=nptl

      nptlsum=0
      npisum=0

c add clusters to cptl, decay them
    
      do i=1,4
        psum(i)=0
        psumNF(i)=0
      enddo
      Psum1=0
      Esum1=0
      scalesum=0
      nscalesum=0
      if(ish.ge.3)then
        write(ifmt,'(5a)')
     .  '    P Ptls      Mclu    M Ptls'
     .  ,'   eta  McluCosh Ptls(bNF)  McluSinh Ptls(bNF) ' !boostedNoFlow
     .  ,'     E & P Ptls (LF)','    E & P Ptls (BjF)'
     .  ,'   Energy   & Pz  & Py  & Px  Ptls (CMS)'   
      endif

      do ncl=1,nclu

        nptl=nptl+1
        call utrepl(nptl,nptla) !copy nptla to nptl

        ! determine idr of cluster for cptl
  
        do k=1,6
          jc(k,1)=0
          jc(k,2)=0
        enddo
        do k=1,3
          j=nint(dGfuse(k,ncl))
          if(j.gt.0)jc(k,1)=j
          if(j.lt.0)jc(k,2)=-j
        enddo
        idr=0
        do  k=1,3
        do  ij=1,2
          if(jc(k,ij).ge.10)idr=7*10**8
        enddo
        enddo
        if(idr/10**8.ne.7)then
          call idenco(jc,ic,ireten)
          if(ic(1).eq.0.and.ic(2).eq.0)then
            ic(1)=100000
            ic(2)=100000
          endif
          idr=8*10**8+ic(1)*100+ic(2)/100
        else 
          idr=idr
     *    +mod(jc(1,1)+jc(2,1)+jc(3,1)+jc(4,1),10**4)*10**4
     *    +mod(jc(1,2)+jc(2,2)+jc(3,2)+jc(4,2),10**4)
          call idtrbi(jc,ibxx1,ibxx2,ibxx3,ibxx4)
          call setibptl(nptl,ibxx1,ibxx2,ibxx3,ibxx4)
        endif
        idrfuse(ncl)=idr
        if(ish.ge.4)
     .  write(ifmt,'(a,2i5,4f8.2,2i5,i14)')'cluster',nfuse(1,ncl)
     .  ,nfuse(2,ncl),dMfuse(ncl),(dGfuse(k,ncl),k=1,3),iflav(ncl)
     .  ,jflav(ncl),idr

        ! decay

        idptl(nptl)=idrfuse(ncl)
        iorptl(nptl)=nptla 
        pptl(4,nptl)=dMfuse(ncl)
        pptl(5,nptl)=dMfuse(ncl)
        ifrptl(1,nptl)=0
        ifrptl(2,nptl)=0
        ip=nptl
        nptlc=nptl
        !write(ifmt,'(a,i6,e12.3,i12)')'MiC decay',ip,pptl(5,ip),idptl(ip)
        call hnbaaa(8,ip,iret)
        if(iret.ne.0)then
          if(nclu.gt.1)then
            write(ifmt,*)'NB decay not possible, try fusing',ncl
            ncluE= ncl
            goto50
          endif
          write(ifmt,*)'WARNING EpoMiF hnbaaa iret; Mass', dMfuse(ncl)
     .    ,'   eta', nfuse(1,ncl),nfuse(2,ncl),'   nclu',nclu
          cycle
        endif
        do i=1,4
          p(i)=0
        enddo
        do n=nptlc+1,nptl
          call getpptl(n,v5(1),v5(2),v5(3),v5(4),dmy)
          do i=1,4
            p(i)=p(i)+v5(i)
          enddo
        enddo
        amsup=sqrt(p(4)**2-p(1)**2-p(2)**2-p(3)**2)
        ppsup=sqrt(p(1)**2+p(2)**2+p(3)**2)
        amclu=pptl(5,ip)
        Pclu1=0
        Eclu1=0
        do neta=nfuse(1,ncl),nfuse(2,ncl)
          eta=etahy(neta)
          Pclu1=Pclu1+dMeta(neta)*sinh(eta)/nspleta(neta)
          Eclu1=Eclu1+dMeta(neta)*cosh(eta)/nspleta(neta)
        enddo
        if(.not.(ppsup.gt.0..or.ppsup.le.0.))then !NaN catch
          write(ifmt,'(a)')'PROBLEM concerning p of cluster'
          write(ifmt,*)'p=',p,pptl(5,ip)
          write(ifmt,'(a)')'ERROR 12052019 EpoMiF: p is NaN'
          stop
        endif

        ! cumulative eta distr    weta(nn=1:netx)    !neta=nfuse(1,ncl)-1+nn

        netx=nfuse(2,ncl)-nfuse(1,ncl)+1
        if(netx.gt.netaclu)stop'ERROR: netaclu too small'
        do nn=1,netx
          neta=nfuse(1,ncl)-1+nn
          weta(nn)=0
          call dmass_i1get(  1  ,neta,w)
          weta(nn)=weta(nn)+max(w,0.)
          do nsurf=1,maxsurf
            call dmass_c1get(nsurf,neta,w)
            weta(nn)=weta(nn)+max(0.,w)
          enddo
        enddo
        do nn=2,netx
          weta(nn)=weta(nn)+weta(nn-1)
        enddo
        if(ish.ge.4)then
          write(ifmt,'(a,2i4,$)')'eta distr',nfuse(1,ncl),nfuse(2,ncl)
          write(ifmt,'(10f6.2)')(weta(nn),nn=1,netx)
        endif

        ! cumulative nsurf distr    wsu(nn,nsu=1:nsux)    !nsurf=nsu-1

        nsux=1+maxsurf
        do nn=1,netx
          neta=nfuse(1,ncl)-1+nn
          eta=etahy(neta)
          call dmass_i1get(  1  ,neta,w)
          wsu(nn,1)=max(w,0.)
          do nsurf=1,maxsurf
            call dmass_c1get(nsurf,neta,w)
            wsu(nn,1+nsurf)=wsu(nn,nsurf)+max(0.,w)
          enddo
          if(ish.ge.4)then
            write(ifmt,'(a,i4,f7.2,3x,a,i3,$)')'eta',neta,eta,'surf',0
            write(ifmt,'(i3,7f6.2)')maxsurf,(wsu(nn,nsu),nsu=1,nsux)
          endif
        enddo

        ! cumulative radtau distr         ! nrad=ntr+1 / nval=nradhy-1 (nsurf 0)
        !     wradtau(nn,nsu,ntr=1:nval)  ! ntau=ntr   / nval=ntauhy   (nsurf>0)

        if(ish.ge.4)write(ifmt,'(a)')'rad/tau distr'
        do nn=1,netx
        do nsu=1,nsux
          neta=nfuse(1,ncl)-1+nn
          eta=etahy(neta)
          nsurf=nsu-1
          if(nsurf.eq.0)then
            nval=nradhy-1
            do nrad=2,nradhy 
              call dmass_i2get(  1  ,neta,nrad,dM) !already cumulative
              wradtau(nn,nsu,nrad-1)=max(0.,dM)
            enddo
          else !nsurf.ne.0
            nval=ntauhy
            do ntau=1,ntauhy 
              call dmass_c2get(nsurf,neta,ntau,dM) !already cumulative
              wradtau(nn,nsu,ntau)=max(0.,dM)
            enddo
          endif 
          if(ish.ge.4)then
            write(ifmt,'(a,i4,f7.2,3x,a,$)')'eta',neta,eta,'surf'
            write(ifmt,'(i3,3x,a,2i4,$)')nsurf,'radtau',1,nval
            do ntr=1,nval,max(1,nval/10)
            write(ifmt,'(f6.2,$)')wradtau(nn,nsu,ntr)
            enddo
            write(ifmt,'(a)')' '
          endif
        enddo
        enddo

        ! cumulative phi distr    wphi(nn,nsu,ntr,nph=1:nphihy-1)   !nphi=nph+1

        do nn=1,netx
        neta=nfuse(1,ncl)-1+nn
        do nsu=1,nsux
        nsurf=nsu-1
        if(nsurf.eq.0)nval=nradhy-1
        if(nsurf.ne.0)nval=ntauhy
        do ntr=1,nval
          do nph=1,nphihy-1 
            nphi=nph+1
            if(nsurf.eq.0)then
              nrad=ntr+1 
              call dmass_i3get(  1  ,neta,nrad,nphi,dM)
            else !nsurf.ne.0
              ntau=ntr
              call dmass_c3get(nsurf,neta,ntau,nphi,dM)
            endif
            wphi(nn,nsu,ntr,nph)=max(0.,dM)
          enddo
        enddo
        enddo
        enddo

        ! daughter indices

        call setifrptl(ip,nptlc+1,nptl)

        ! determine random neta, eta 

        nit=0
        xi=100000
 200    nit=nit+1

        do n=nptlc+1,nptl

          r=rangen()*weta(netx) 
          w=0.
          i=0
          do while (w.lt.r)
            i=i+1
            w=weta(i) 
          enddo          
          nn=i 
          neta=nfuse(1,ncl)-1+nn
          call set2ibptl(n,1,neta)
          eta=etahy(neta)

        enddo

        ! compare E,P ptls (only Bj flow) with Mcosh,Msinh (sum eta slices)

        nnfuse=nfuse(2,ncl)-nfuse(1,ncl)+1
        if(nnfuse.gt.1)then 
          Pclu2=0 
          Eclu2=0 
          do n=nptlc+1,nptl
            call get2ibptl(n,1,neta)
            eta=etahy(neta)
            call getpptl(n,v5(1),v5(2),v5(3),v5(4),v5(5))
            amt=sqrt(v5(1)**2+v5(2)**2+v5(5)**2)
            p(3)=v5(3)
            p(4)=v5(4)
            rap=sign(1.,p(3))*alog((p(4)+abs(p(3)))/amt) 
            Pclu2=Pclu2+amt*sinh( eta+rap )
            Eclu2=Eclu2+amt*cosh( eta+rap )
          enddo
          x = max( abs((Eclu2-Eclu1)/Eclu1) , abs((Pclu2-Pclu1)/Eclu1) ) 
          if(x.lt.xi)then
            xi=x
            do n=nptlc+1,nptl
              call get2ibptl(n,1,neta) 
              call set2ibptl(n,2,neta) ! keep best neta selection
            enddo
          endif 
          xb=x
          if(nit.le.500*nnfuse)goto 200
          do n=nptlc+1,nptl
            call get2ibptl(n,2,neta)
            call set2ibptl(n,1,neta) ! take best neta selection
          enddo
        endif

        ! check energy per cluster

        do i=3,4
          PcluNF(i)=0
        enddo
        do n=nptlc+1,nptl
          nptlsum=nptlsum+1
          call getidptl(n,idxxxx)
          if(abs(idxxxx).eq.120)npisum=npisum+1
          call getpptl(n,v5(1),v5(2),v5(3),v5(4),v5(5))
          amt=sqrt(v5(1)**2+v5(2)**2+v5(5)**2)
          call get2ibptl(n,1,neta)
          eta=etahy(neta)
          p(3)=v5(3)
          p(4)=v5(4)
          rap=0
          iok=0
          if(amt.ne.0.)iok=1
          if(iok.eq.1.and.(p(4)+abs(p(3)))/amt.gt.0.)then
            rap=sign(1.,p(3))*alog((p(4)+abs(p(3)))/amt) 
          else
            write(ifmt,'(a)')'WARNING EpoMiF FP issue'
          endif
          PcluNF(3)=PcluNF(3)+amt*sinh(eta+rap) 
          PcluNF(4)=PcluNF(4)+amt*cosh(eta+rap) 
        enddo
        do i=3,4
          psumNF(i)=psumNF(i)+PcluNF(i)
        enddo
        Psum1=Psum1+Pclu1
        Esum1=Esum1+Eclu1
        if(ish.ge.3)then
          write(ifmt,'(3f10.2,$)')ppsup,amclu,amsup  ! P Ptls   Mclu   M Ptls
          write(ifmt,'(2i3,$)')nfuse(1,ncl),nfuse(2,ncl) ! eta
          write(ifmt,'(4f10.2,$)')Eclu1,PcluNF(4),Pclu1,PcluNF(3) ! McluCosh Ptls(bNF)  McluSinh Ptls(bNF)
        endif

        ! consider flow
 
        do n=nptlc+1,nptl

          call setiorptl(n,ip)
          call setjorptl(n,0)
          call setistptl(n,0)
          call setifrptl(n,0,0)
          call get2ibptl(n,1,neta)
          eta=etahy(neta)
          nzz=neta-int(nzhy/2.)+int(mxeflo/2.)
          nn=neta-nfuse(1,ncl)+1

          ! determine random nsurf 

          r=rangen()*wsu(nn,nsux) 

          w=0.
          i=0
          do while (w.lt.r)
            i=i+1
            w=wsu(nn,i)
          enddo          
          nsu=i
          nsurf=nsu-1

          ! determine random nrad / ntau

          if(nsurf.eq.0)nval=nradhy-1 
          if(nsurf.ne.0)nval=ntauhy
          r=rangen()*wradtau(nn,nsu,nval)
          w=0.
          i=0
          do while (w.lt.r)
            i=i+1
            w=wradtau(nn,nsu,i)
          enddo          
          ntr=max(1,i) 
          if(nsurf.eq.0)nrad=ntr+1
          if(nsurf.ne.0)ntau=ntr

          ! print phi table

          if(ish.ge.4)then
            write(ifmt,'(a,i4,f7.2,$)')'========> eta',neta,eta
            write(ifmt,'(3x,a,i3,3x,a,i4,a,$)')'surf',nsurf
     .      ,'nradtau',ntr,'  '
            wx=0
            do nph=1,nphihy-1
              w=wphi(nn,nsu,ntr,nph)
              iw=100*(w-wx)/wphi(nn,nsu,ntr,nphihy-1)
              if(iw.lt.0)stop'ERROR 05052019' 
              if(iw.gt.99)then
                write(ifmt,'(i4,1x,$)')iw
              elseif(iw.gt.9)then
                write(ifmt,'(i3,1x,$)')iw
              else
                write(ifmt,'(i1,$)')iw
              endif
              wx=w
            enddo
            write(ifmt,'(a)')' '
          endif

          ! determine random nphi, phi

          r=rangen()*wphi(nn,nsu,ntr,nphihy-1) 
          w=0.
          i=0
          do while (w.lt.r)
            i=i+1
            w=wphi(nn,nsu,ntr,i)
          enddo          
          nph=i 
          nphi=nph+1
          phi=phihy(nphi)

          ! determine flow 4-velocity u

          if(nsurf.eq.0)then
            do k=1,3
              u(k)=velzz(k,neta,nrad,nphi)   
            enddo
          else!nsurf.gt.0
            do k=1,3
              u(k)=velaa(nsurf,k,neta,ntau,nphi)
            enddo
          endif
          gmx=max(1e-8, 1.-u(1)**2-u(2)**2-u(3)**2)
          !~~~~~~~
          !if(n.eq.702)then
          !!variable u is actually velocity v at this point
          !v=sqrt(u(1)**2+u(2)**2+u(3)**2)
          !write(ifmt,*)'CHECKfreeze',n,v,gmx,nsurf,neta,ntau,nphi
          !endif
          !~~~~~~~
          gm=1./sqrt(gmx)
          u(4)=1
          do k=1,4
            u(k)=u(k)*gm
          enddo
          u5=sqrt(u(4)**2-u(1)**2-u(2)**2-u(3)**2)
          if(ish.ge.4)then
            write(ifmt,'(a,i4,f7.2,3x,a,i3,3x,a,i4,3x,a,i4,f7.2,$)')
     .      'particle: eta',neta,eta,'surf',nsurf,'nradtau',ntr,
     .      'phi',nphi,phi
            write(ifmt,'(3x,a,4i4,f6.2)')'4vel',(int(100*u(k)),k=1,4),u5
          endif

          ! temporarily

          call setxorptl(n,u(1),u(2),u(3),u(4))
          call set2ibptl(n,2,nsurf)
          call set2ibptl(n,3,ntr)
          call set2ibptl(n,4,nphi)
          !if(nsurf.eq.0)nrad=ntr+1;if(nsurf.ne.0)ntau=ntr
          !phi=phihy(nphi)
          call setradptl(n,rangen()) !to avoid changing random number sequence later

        enddo !n

        ! Boost back into working frame (Bjorken frame, Minkowski space vectors)
        !  to account for the flow
        ! Then boost to CMS frame

        do i=1,4
          psu(i)=0 
          psa(i)=0
          psi(i)=0
        enddo
        Eclu3=0 
        Eclu4=0 
        do n=nptlc+1,nptl
          call getxorptl(n,u(1),u(2),u(3),u(4))
          call getpptl(n,p(1),p(2),p(3),p(4),ama)
          do i=1,4 
            psu(i)=psu(i)+p(i)
          enddo
          Eclu3=Eclu3+p(4)
          !boost back from LF to Bjorken frame
          call utlob3(-1,u(1),u(2),u(3),u(4),1e0,p(1),p(2),p(3),p(4))
          Eclu4=Eclu4+p(4)
          do i=1,4 
            psa(i)=psa(i)+p(i)
          enddo
          ! boost to CMS frame
          amt=sqrt(p(1)**2+p(2)**2+ama**2)
          ptr=sqrt(p(1)**2+p(2)**2)
          pha=sign(1.,p(2))*acos(p(1)/ptr)
          rap=sign(1.,p(3))*alog((p(4)+abs(p(3)))/amt) 
          call get2ibptl(n,1,neta)
          dleta=(etahy(2)-etahy(1))
          eta=etahy(neta)-dleta/2+rangen()*dleta
          call setpptl(n,ptr*cos(pha+phinull)
     .                  ,ptr*sin(pha+phinull)
     .                  ,amt*sinh(eta+rap)
     .                  ,amt*cosh(eta+rap)
     .                  ,ama) 
          call getpptl(n,v5(1),v5(2),v5(3),v5(4),v5(5))
          !~~~~~~~~~~~~ 
          !if(ptr.gt.5.)then
          !v=sqrt(u(1)**2+u(2)**2+u(3)**2)/u(4)
          !!if(v.gt.0.7)then
          !write(ifmt,'(a,i6,2f10.4)')'CHECKfreeze',n,ptr,v
          !!stop'CHECKfreeze'
          !endif
          !~~~~~~~~~~~~ 
          do i=1,4 
            psi(i)=psi(i)+v5(i)
            psum(i)=psum(i)+v5(i)
          enddo
          call get2ibptl(n,2,nsurf)
          call get2ibptl(n,3,ntr)
          call get2ibptl(n,4,nphi)
          !to improve rad determination via interpolation
          !for the different variables: 
          !nodes:   i-1 --*-- i --*-- i+1
          !          1        2        3
          ! AB :          A       B   weights  wA1 wA2 wB2 wB3  (0.5, at boarders 0,1)
          ! W:  A+r(B-A) = A(1-r)+Br  weights  wA = 1-r  wB = r
          ! weights for 1,2,3:   w1 = (1-r)*wA1, w2 = (1-r)*wA2+r*wB2, w3=r*wB3
          !then compute rad=.. correspondingly
          dphi=phihy(2)-phihy(1)
          phi=phihy(nphi)-dphi/2+rangen()*dphi
          if(nsurf.eq.0)then
            nrad=ntr+1
            drad=radhy(2)-radhy(1)
            rad=radhy(nrad)-drad/2+rangen()*drad
            ntau=1
            tau=getHydynTauhy(ntau)  
          else !nsurf.ne.0
            ntau=ntr
            taumn=getHydynTauhy(ntau)
            if(ntau.gt.1)taumn=(getHydynTauhy(ntau-1)
     .           +getHydynTauhy(ntau))/2
            taumx=getHydynTauhy(ntau)
            if(ntau.lt.ntauho)taummx=(getHydynTauhy(ntau)
     .           +getHydynTauhy(ntau+1))/2
            tau=taumn+rangen()*(taumx-taumn)
            rad=radaa(nsurf,neta,ntau,nphi)      ! ----------> should be improved -> interpolation 
          endif 
          ncount=ncount+1
          if(ncount.le.1000)then
            FoEta(ncount)=eta
            FoTau(ncount)=tau
          endif
          call setxorptl(n,rad*cos(phi+phinull)
     .                    ,rad*sin(phi+phinull)
     .                    ,tau*sinh(eta)
     .                    ,tau*cosh(eta) )
          call getxorptl(n,dmy1,dmy2,dmy3,xo4)
          call getidptl(n,idxx)
          call getpptl(n,dmy1,dmy2,dmy3,p4xx,p5xx)
          tiv1=xo4
          call idtau(idxx,p4xx,p5xx,taugm)
          call getradptl(n,r)
          tiv2= tiv1+taugm*(-alog(r)) 
          call settivptl(n,tiv1,tiv2)
          call setityptl(n,60)
          call setifrptl(n,0,0)
          call setiorptl(n,nptlc)
          call setjorptl(n,0)
          call setistptl(n,0)
        enddo!n
        Pclu3=sqrt(psu(1)**2+psu(2)**2+psu(3)**2)
        Pclu4=sqrt(psa(1)**2+psa(2)**2+psa(3)**2)

        if(ish.ge.3)
     .  write(ifmt,'(8f10.2)')
     .   Eclu3,Pclu3 ! E and P Ptls (LF) 
     .  ,Eclu4,Pclu4 ! E and P Ptls (BjF) 
     .  ,psi(4),psi(3),psi(2),psi(1) ! E and Pz and Py and Px Ptls (CMS) 

      enddo !ncl

      write(ifmt,'(a,f10.2,i6,i9,i9)')
     .'Clusters dMsum nclu nptlsum npisum'
     .,dMsum,nclu,nptlsum,npisum
      if(ish.ge.3)then
        write(ifmt,'(a,$)')'Summing all clusters         '  
        write(ifmt,'(6x,4f10.2,$)')Esum1,psumNF(4),Psum1,psumNF(3)   
        write(ifmt,'(40x,4f10.2)')(psum(i),i=4,1,-1)
      endif

      ncountmic=ncountmic+1     
      if(ncountmic.eq.1)write(ifmt,'(a)')
     .'WARNING rad determination to be improved' !see above "improve rad..." discussion

      ncountmicplot=1
      if(ncountmic.eq.ncountmicplot)then
        open(unit=101,file=fndt(1:nfndt-5)//'-mic.py',status='unknown')
        write(101,'(a)')'#!/usr/bin/python                             '
        write(101,'(a)')'import matplotlib.pyplot as plt               '
        write(101,'(a)')'##############################################'
        write(101,'(a)')'plt.figure(figsize=(7,5))                     '
        write(101,'(a)')'plt.rcParams["font.size"] = 18                '
        write(101,'(a)')'plt.xlabel(u"space-time rapidity \u03b7")     '
        write(101,'(a)')'plt.ylabel(u"dM / d\u03b7")             '
        write(101,'(a)')'plt.xlim(-6,6)       '
        !write(101,'(a)')'plt.ylim(0,9.5)     '
        write(101,'(100(a,f8.2))')'X=['
     .  ,(etahy(neta),',',neta=1,nzhy-1),etahy(nzhy),']'         
        write(101,'(100(a,f8.2))')'Y=['
     .  ,(dMetaSum(neta),',',neta=1,nzhy-1),dMetaSum(nzhy),']'       
        write(101,'(a)')'Y[:]=[x / 1 for x in Y]     '
        write(101,'(a)')'plt.plot(X,Y,"rs",markersize=10);     '   !"bs"
        !write(101,'(a)')'plt.legend();     '
        write(101,'(a)')'plt.tight_layout();;     '
        write(101,'(a)')'##############################################'
        write(101,'(a)')'plt.figure(figsize=(7,5))                     '
        write(101,'(a)')'plt.rcParams["font.size"] = 18                '
        write(101,'(a)')'plt.xlabel(u"space-time rapidity \u03b7")     '
        write(101,'(a)')'plt.ylabel(u"proper time \u03c4")             '
        write(101,'(a)')'plt.xlim(-7,7)       '
        !write(101,'(a)')'plt.ylim(0,9.5)     '
        ncount=min(1000,ncount)
        if(ncount.ge.2)then
        write(101,'(100(a,f8.2))')'X=['
     .  ,(FoEta(n),',',n=1,ncount-1),FoEta(ncount),']'         
        write(101,'(100(a,f8.2))')'Y=['
     .  ,(FoTau(n),',',n=1,ncount-1),FoTau(ncount),']'       
        endif
        write(101,'(a)')'Y[:]=[x / 1 for x in Y]     '
        write(101,'(a)')'plt.plot(X,Y,"rs",markersize=10);     '   !"bs"
        !write(101,'(a)')'plt.legend();     '
        write(101,'(a)')'plt.tight_layout();;     '
        !write(101,'(a)')'##############################################'
        !write(101,'(a)')'label = "data 1");     '
        !write(101,'(a)')'plt.plot([ 2.5,2.5],[7,2], "-r" );     '
        !write(101,'(a)')'plt.plot([ 3.5,3.5],[7,2], "-r" );     '
        !write(101,'(a)')'plt.text(2,2,"text");     '
        !tau:  u"proper time \u03c4"
        write(101,'(a)')'##############################################'
        write(101,'(a)')'plt.show();     '
        close(101)
      endif

      call checkGrandCanon(2,0,0.,ve3) !plot

 9999 continue

      call utprix('epomif',ish,ishini,4)
      if(kfrout.ge.3.and.kfrout.lt.20)stop' (Testrun)'
      end


c###################################################################################
      subroutine EpoF(iflag)
c###################################################################################

c-----------------------------------------------------------------------------------
c   Grand-canonical FO via FO hyper-surface  
c-----------------------------------------------------------------------------------

c-----------------------------------------------------------------------------------
c kfrout >=20 ... Nothing (returns immediately)
c         1 ..... Using T for FO surface
c         2 ..... Using epsilon for FO surface
c iflag  1 ... initialization, concerning particle list essentially
c        2 ... MC freeze out EbE
c-----------------------------------------------------------------------------------

#include "aaa.h"
#include "ho.h"
      real probdi(100),dprob
      common/crap/nraphi,jfac /crapi/delrapi
      common/ceefrac/eefrac
      common/citsy/itsy(4,4)
      real u(3),q(4),s(4)
      real wrii(2),wwii(2),wxii(2),wyii(2),wzii(2),weii(2)
      real wt0ii(10,nradhxx,2)
      real w1ii(2),w2ii(2),w3ii(2)
      real wsii(4,2),wtii(10,2)
      real wi(2),wj(2),wk(2),wq(2)
      integer mi(2),mj(2),mk(2),mq(2)
      common/jcen/jcentr,jmxcentr
c      common/cratioee/ratioee,ratioeex(ntauhxx)
      data ncntrffo/0/
      save ncntrffo
      data ncnt3fo/0/
      save ncnt3fo

      real getHydynRatioeex
      integer getAccumJerr

      if(kfrout.ge.20)return

      if(iflag.eq.1)eefrac=1

      iigg=0 !1 = Test activated (comparing pions and Omegas), 0 = not
      if(iigg.ne.0.and.iigg.ne.1)stop'####### ERROR 05122016 #######'
      inupi=0 !pions
      ggsum1=0
      ggsum2=0
      inuom=0 !Omegas
      ggsum3=0
      ggsum4=0

      if(ntauhy.eq.0)return

      if(nofreeze.eq.1.or.eefrac.le.0.)then
        if(iflag.eq.2.and.eefrac.le.0.)then
          write(ifch,'(a)')'EpoF: eefrac < 0 => STOP'
          write(ifmt,'(a)')'EpoF: eefrac < 0 => STOP'
          stop'\n\n ERROR  EpoF: eefrac < 0 \n\n'
        endif
        return
      endif

      GOTO (1,10),iflag

  1   continue !initialization, concerning particle list essentially
      call hynDefineParticles(1)
      call hynDefineRapScale
      return

  10  continue !here starts the MC freeze out EbE
      ncntrffo=ncntrffo+1
      if(ncntrffo.eq.1)write(ifmt,'(a)')'epos freeze out EpoF'
      call clop(3)

      !special case
      if(ireadhyt.gt.0.and.jcorona.le.0)then
        call GetImpactParam
      endif

      ier=1
      nsimu=1 + 999*iigg
      bimp=bimevt

      !write(ifmt,'(a,$)')' RestFrameFO ...    '
      ier=0
      pi=3.1415927
      hbar=0.197327

      suel1=0
      suel2=0

      !~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
      !several surface layers possible 
      ! so we sum over them (nsurf)
      !(but small contribution from nsurf > 1 
      !~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

      do nsurf=1,maxsurf !~~~~~sum over surface layers~~~~~~

      do i=1,100
       probdi(i)=0
      enddo
      dprob=0.1

      ntauho=ntauhy  !EpoF

      !write(ifmt,*)'++++++',ntauho,jfac,nsimu

      dphi=phihy(2)-phihy(1)
      dtau=getHydynTauhy(2)-getHydynTauhy(1)
      deta=etahy(2)-etahy(1)
      dleta=(etahy(2)-etahy(1))/jfac
      dall=dphi*dtau*dleta*0.125
      !write(ifmt,*)'++++++',dphi,dtau,deta,dleta,dall

      !~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~   
      !the surface is given as radius(phi,eta,tau)
      ! so we sum over phi,eta,tau
      !~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~   

      do ntau=1,ntauho !~~~~~~ tau sum ~~~~~~~~~~~

      tau=getHydynTauhy(ntau)
      ftau=2
      if(ntau.eq.1.or.ntau.eq.ntauho)ftau=1

      do neta=2,nzhy  !~~~~~~ eta sum ~~~~~~~~~~~~

      do nphi=1,nphihy !~~~~~~ phi sum ~~~~~~~~~~

        phi=phihy(nphi)
        fphi=2
        if(nphi.eq.1.or.nphi.eq.nphihy)fphi=1

        do ii=1,2
          if(ii.eq.1)mt=neta-1
          if(ii.eq.2)mt=neta
          wwii(ii)=vlmaa(nsurf,mt,ntau,nphi)
          wwii(ii)=max(0.,wwii(ii))
          wxii(  ii)=velaa(nsurf,1,mt,ntau,nphi)   ! x velocity
          wyii(  ii)=velaa(nsurf,2,mt,ntau,nphi)   ! y velocity
          wzii(  ii)=velaa(nsurf,3,mt,ntau,nphi)   ! z velocity
          wrii(  ii)=radaa(nsurf,  mt,ntau,nphi)   ! radius
          weii(  ii)=epsaa(nsurf,  mt,ntau,nphi)   ! energy density
          w1ii(  ii)=baraa(nsurf,1,mt,ntau,nphi)   ! B- 
          w2ii(  ii)=baraa(nsurf,2,mt,ntau,nphi)   ! Q- densities
          w3ii(  ii)=baraa(nsurf,3,mt,ntau,nphi)   ! S-   
          do k=1,4                                 !   dSigma_nu
            wsii(k,ii)=suraa(nsurf,k,mt,ntau,nphi) ! -------------- 
          enddo                                    ! dtau deta dphi
          do k=1,10   
            wtii(k,ii)=emuaa(nsurf,k,mt,ntau,nphi) ! T^mu,nu 
          enddo 
          if(ntau.eq.1)then
            do nrad=1,nradhy
              do k=1,10   
              wt0ii(k,nrad,ii)=emuzz(k,mt,nrad,nphi)  ! T^mu,nu at initial tau
              enddo
            enddo
          endif
c          if(neta.eq.2.and.ii.eq.1.or.ii.eq.2) 
c     .    write(ifch,'(a,5i4,4f8.3)')'EpoF  sur cen mt t phi  e b123'
c     .   ,nsurf,mt,ntau,nphi
c     .   ,epsaa(nsurf,mt,ntau,nphi)
c     .   ,(baraa(nsurf,i3,mt,ntau,nphi),i3=1,3) 
        enddo

        do j=0,jfac  !~~~~~~~j sum ~~~~~~~~~~

          !divide each eta interval into jfac subintervals
          !to get a finer eta grid.
          !values of variables on the surface are 
          !obtained via linear interpolation

          eta=etahy(neta-1)+j*(etahy(2)-etahy(1))/jfac

          !meaning of the following variables: 
          ! v1   x velocity
          ! v2   y velocity
          ! v3   z velocity
          ! rad  radius
          ! e    energy density
          ! b1   B-density
          ! b2   Q-density
          ! b3   S-density
          ! s    dSigma_nu
          ! t    T^mu,nu
          v1=  wxii(  1)+j/float(jfac)*(wxii(  2)-wxii(  1))
          v2=  wyii(  1)+j/float(jfac)*(wyii(  2)-wyii(  1))
          v3=  wzii(  1)+j/float(jfac)*(wzii(  2)-wzii(  1))
          rad= wrii(  1)+j/float(jfac)*(wrii(  2)-wrii(  1))
          e=   weii(  1)+j/float(jfac)*(weii(  2)-weii(  1))
          b1=  w1ii(  1)+j/float(jfac)*(w1ii(  2)-w1ii(  1))
          b2=  w2ii(  1)+j/float(jfac)*(w2ii(  2)-w2ii(  1))
          b3=  w3ii(  1)+j/float(jfac)*(w3ii(  2)-w3ii(  1))
          do k=1,4
          s(k)=wsii(k,1)+j/float(jfac)*(wsii(k,2)-wsii(k,1))
          enddo 
 
          f=2
          if(j.eq.0.or.j.eq.jfac)f=1
          fall=fphi*ftau*f
          dVs=wwii(1)+j/float(jfac)*(wwii(2)-wwii(1))
          dVs=dVs*dall*fall
          gmx=max(1e-8, 1.-v1*v1-v2*v2-v3*v3)
          gm=1./sqrt(gmx)
          ey=max(e,0.)
          if(ey.ge.oEeos-0.01)then
            ncnt3fo=ncnt3fo+1
            call setAccumJerr(13,getAccumJerr(13)+1)
            if(ncnt3fo.le.50)then
              write(ifmt,*)'ey=',ey,'  b1=',b1,'  b2=',b2,'  b3=',b3
     .        ,'  ',nsurf,neta,ntau,nphi
              write(ifmt,*)weii(1),weii(2),j/float(jfac)
              write(ifmt,*)epsaa(nsurf,neta-1,ntau,nphi)
     .        ,epsaa(nsurf,neta-1,ntau,nphi)
              write(ifmt,*)epsaa(nsurf,neta,ntau,nphi)
     .        ,epsaa(nsurf,neta,ntau,nphi)
            endif
          endif

          ! equation of state to get tmpfo, chmB, chmQ, chmS, ppp
          call PimEos(ey, b1, b2, b3,   tmpfo, chmB, chmQ, chmS, ppp)
          
          !for the obtained chmB, chmQ, chmS
          !determine coefficients for interpolation
          !using predefined tables of yields of particles
          !for given chemical potentials

c          write(ifch,'(a,4f8.3,4x,4f8.3)')
c     .      'EpoFEos  eps b123  T ch123:'
c     .      ,ey, b1, b2, b3,   tmpfo, chmB, chmQ, chmS

          cox=volex*ppp/tmpfo
          ! ~~~~~~~~~~~~~~~~ T
          mte=2
          do while(mte.lt.mtemp.and.temfo(mte).lt.tmpfo)
            mte=mte+1
          enddo
          mq(1)=mte-1
          mq(2)=mte
          frac=(tmpfo-temfo(mq(1)))/(temfo(mq(2))-temfo(mq(1)))
          frac=min(frac,1.0)
          frac=max(frac,0.0)
          wq(1)=1-frac
          wq(2)=frac
          chm=chmB *iochem ! ~~~~~~~~~~~~~~~~ B
          mch=2
          do while(mch.lt.mchem.and.chepoB(mch).lt.chm)
            mch=mch+1
          enddo
          mi(1)=mch-1
          mi(2)=mch
          frac=(chm-chepoB(mi(1)))/(chepoB(mi(2))-chepoB(mi(1)))
          frac=min(frac,1.0)
          frac=max(frac,0.0)
          wi(1)=1-frac
          wi(2)=frac
          chm=chmQ *iochem ! ~~~~~~~~~~~~~~~~ Q
          mch=2
          do while(mch.lt.mchem.and.chepoQ(mch).lt.chm)
            mch=mch+1
          enddo
          mj(1)=mch-1
          mj(2)=mch
          frac=(chm-chepoQ(mj(1)))/(chepoQ(mj(2))-chepoQ(mj(1)))
          frac=min(frac,1.0)
          frac=max(frac,0.0)
          wj(1)=1-frac
          wj(2)=frac
          chm=chmS *iochem ! ~~~~~~~~~~~~~~~~ S
          mch=2
          do while(mch.lt.mchem.and.chepoS(mch).lt.chm)
            mch=mch+1
          enddo
          mk(1)=mch-1
          mk(2)=mch
          frac=(chm-chepoS(mk(1)))/(chepoS(mk(2))-chepoS(mk(1)))
          frac=min(frac,1.0)
          frac=max(frac,0.0)
          wk(1)=1-frac
          wk(2)=frac
          ! ~~~~~~~~~~~~~~~~
          finc=10
          volu=abs(dVs)*finc
          io3=ioclude-3
          if(io3.lt.1.or.io3.gt.2)
     .     stop'in EposFreezeOut: wrong ioclude (150808) '

          do m=0,nspes+2 !sum over hadron species

            !do interpolation to get the values for 
            ! actual chemical potentials

            !ggstat(m) is the thermal yield summed aup for all
            !species up to number m

            !so ggstat(nspes) is the total yield
            ! ggstat(nspes+2) is simply unity here
            ! we use the option io3 = 1

            ggstat(m)=0
            do iq=1,2
            do i=1,2
            do jj=1,2
            do k=1,2
            ggstat(m)=ggstat(m)
     .       +wq(iq)*wi(i)*wj(jj)*wk(k)
     .        * ffstat(mq(iq),io3,mi(i),mj(jj),mk(k),m)
            enddo
            enddo
            enddo
            enddo

          enddo

c          write(ifch,'(a,2f8.3,2i4,2f8.3,i7,2f8.3)')
c     .      'EpoF eps T mq wq id ff1 ff2:'
c     .     ,ey,tmpfo,mq(1),mq(2),wq(1),wq(2),ispes(1)
c     .     ,ffstat(mq(1),1,11,11,11,1),ffstat(mq(2),1,11,11,11,1)

          if(iigg.eq.1)then !~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
             write(*,'(a,$)')' iigg '
             mte=5 
             iixx=1
             write(*,'(2f9.5,3x,$)')
     .       (ffstat(mte,1,11,11,11, iixx ) 
     .           - ffstat(mte,1,11,11,11, iixx-1 ))
     .              / ffstat(mte,1,11,11,11, nspes ) 
     .       , (ggstat( iixx ) - ggstat( iixx-1 )) / ggstat( nspes )  
             ggsum1=ggsum1+(ggstat(iixx)-ggstat(iixx-1))/ggstat(nspes)
             iixx=2
             write(*,'(3x,2f9.5,3x,$)')
     .        (ffstat(mte,1,11,11,11, iixx ) 
     .            - ffstat(mte,1,11,11,11, iixx-1 ))
     .              / ffstat(mte,1,11,11,11, nspes ) 
     .       , (ggstat( iixx ) - ggstat( iixx-1 )) / ggstat( nspes ) 
             ggsum2=ggsum2+(ggstat(iixx)-ggstat(iixx-1))/ggstat(nspes)
             iixx=348
             write(*,'(3x,2f9.5,3x,$)')
     .        (ffstat(mte,1,11,11,11, iixx ) 
     .            - ffstat(mte,1,11,11,11, iixx-1 ))
     .              / ffstat(mte,1,11,11,11, nspes ) 
     .       , (ggstat( iixx) - ggstat( iixx-1 )) / ggstat( nspes ) 
             ggsum3=ggsum3+(ggstat(iixx)-ggstat(iixx-1))/ggstat(nspes)
             iixx=349
             write(*,'(3x,2f9.5,3x)')
     .        (ffstat(mte,1,11,11,11, iixx ) 
     .            - ffstat(mte,1,11,11,11, iixx-1 ))
     .              / ffstat(mte,1,11,11,11, nspes ) 
     .       , (ggstat( iixx ) - ggstat( iixx-1 )) / ggstat( nspes ) 
             ggsum4=ggsum4+(ggstat(iixx)-ggstat(iixx-1))/ggstat(nspes)
          endif !~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

          if(getHydynRatioeex(ntau).eq.0.)
     .         stop'####### ERROR 19022016 #######'
          yie=volu*ggstat(nspes)*ggstat(nspes+2)*exp(-cox)
     .         *fofac*eefrac/getHydynRatioeex(ntau)
           
          !volu              volume
          !ggstat(nspes)     total yield
          !ggstat(nspes+2)   1
          !exp(-cox)         excluded volume term
          !fofac             correction factor
          !eefrac            small (and well controlled) correction factor
          !                  to correct for fluid jet interactions
          !                   (1 in the simplest case) 

          !we use the yield in the first place as proposal
          ! and later  accept with some formula to have finally Cooper Frye

          do nsim=1,nsimu ! nsimu=1 usually

            np=yie
            if(rangen().le.yie-np)np=np+1

            if(np.gt.0)then

              do n=1,np !the number of hadrons is determined from yie

                r=rangen()*yie
                m=0
                do while(yie*ggstat(m)/ggstat(nspes).lt.r)
                 m=m+1
                enddo
                id=ispes(m)
                am=aspes(1,m)
                iQ=iQspes(m)
                ch=0
                do iq=1,2
                do i=1,2
                do jj=1,2
                do k=1,2
                ch=ch+ wq(iq)*wi(i)*wj(jj)*wk(k)
     .                  *chot(mq(iq),mi(i),mj(jj),mk(k),m)
                enddo
                enddo
                enddo
                enddo
                cot=ch+cohi(m)-cox
                !print*,'+++++++',cot,volex*ppp/tmpfo
                js=jspes(m)
                call setAccumJerr(19,getAccumJerr(19)+1)
                if(js.eq.-1.and.cot.gt.aspes(io3,m)/tmpfo)then  !mesons
                 call setAccumJerr(20,getAccumJerr(20)+1) !mu > m => BE condensate
                 frac=0
                endif
                !generate thermal momenta in local frame
                x=hynRanBoseFermi(aspes(io3,m),cot,js,tmpfo,istat)
                p=x*tmpfo
                e=sqrt(p**2+aspes(1,m)**2)
                ex=sqrt(p**2+aspes(io3,m)**2)
                !write(ifmt,*)id,e
                u(3)=2.*rangen()-1.
                angle=2.*pi*rangen()
                u(1)=sqrt(1.-u(3)**2)*cos(angle)
                u(2)=sqrt(1.-u(3)**2)*sin(angle)
                q(1)=p*u(1)
                q(2)=p*u(2)
                q(3)=p*u(3)
                q(4)=ex
                !boost into working frame
                call utlob3(-1, v1*gm , v2*gm , v3*gm , gm ,1e0
     .                   , q(1), q(2), q(3), q(4))
                fof=hynFOFactor2(nsurf,s,q)
                fof=fof*dall*fall
                probab=fof/volu/e
                ij=1+probab/dprob
                if(ij.ge.1.and.ij.le.100)
     .           probdi(ij)=probdi(ij)+1
                !~~~~~ accept or reject ~~~~~~
                if(rangen().le.probab)then !accept
                  if(io3.eq.2)then
                   q(4)=sqrt(q(1)**2+q(2)**2+q(3)**2+aspes(1,m)**2)
                  endif
                  if((q(4)+q(3)).ne.0..and.(q(4)-q(3)).ne.0.)then
                   rap=eta +  0.5*log((q(4)+q(3))/(q(4)-q(3)))
                  else
                   rap=eta + 4*(rangen()-0.5)
                   write(ifmt,*)'warning in RestFrameFO: q4 +- q3 = '
     .             ,(q(4)+q(3)),(q(4)-q(3))
                  endif
                  pt2=q(1)**2+q(2)**2
                  ptr=sqrt(pt2)
                   !if(ptr.gt.10)print*,'+++++++',ptr,gm
                  pha=sign(1.,q(2))*acos(q(1)/ptr)
                  if(pha.lt.0.)pha=pha+2*pi
                  if(iigg.eq.1)then !Test
                    if(abs(id).eq.211)inupi=inupi+1  
                    if(abs(id).eq.3334)inuom=inuom+1 
                  endif
                  call hynFOStore(id,iQ,ptr,pha,am,rap,rad,phi,tau,eta)
                !~~~~~~~~~~~~~~~~~~~~~~~~~~~
                else
                endif

              enddo !np

            endif !np.gt.0

          enddo !nsim

        enddo !j

c  99    continue

      enddo !nphi

      enddo !neta

      enddo !ntau

      enddo !nsurf

      !write(ifmt,'(a)')'  done'
      if(ish.ge.3)call alist('list after FoPtl&',1,nptl)

      if(iigg.eq.1)then !Test
        if(inupi.gt.0)anurat=1000.*inuom/inupi 
        print*,' iigg ',inuom,inupi,'   ',anurat,'    '
     .         ,1000*(ggsum3+ggsum4)/(ggsum1+ggsum2)
        stop' iigg '
      endif

      return

      end


c---------------------------------------------------------------------
       subroutine rhytabx(nfr)  !read hydro table
c---------------------------------------------------------------------
#include "aaa.h"
#if __ROOT__
      if(ireadhyt.gt.0)then
        call setHo(1,0,1)
        if(nfr.eq.0)then
          call InitializeHyperbola !does not affect results but necessary
          do ird=1,ireadhyt
            call ReadHydroTable
            if(ispherio.ne.0.or.jspherio.ne.0.and.ird.eq.1)then
              write(ifmt,'(a)')'enter spherio; nfr=-1 ... '
              call clop(2)
              call spherio(nrevt,-1)
              call clop(1)
              write(ifmt,'(a)')'exit spherio'
            endif
            if(ihlle.ne.0.or.ioeos.ne.0.and.ird.eq.1)then
              write(ifmt,'(a)')'enter hlle; nfr=-1 ... '
              call hlle(-1)
              write(ifmt,'(a)')'exit hlle'
            endif
            call fros(1)
          enddo
         call fros(4)
         call fros(3)
         call epof(1)
        endif
        call EpoMiF
        call epof(2)
        call decayall(1,0)
        call xxEos
      endif
#else
        stop "ROOT needed for rhytabx"
#endif
      end

c---------------------------------------------------------------------
      subroutine xfro
c---------------------------------------------------------------------
#include "aaa.h"
      if(iSpaceTime.eq.7)call xFreezeOutTauEta
      if(iSpaceTime.eq.7)call xFreezeOutTauX
      end

c---------------------------------------------------------------------
      subroutine GetImpactParam
c----------------------------------------------------------------------
#include "aaa.h"
      data ncntbim/0/
      save ncntbim
      ncntbim=ncntbim+1
      if(ncntbim.eq.1)call getbcheck
      call getb(b1,b2)
      end
      !~~~~~~~~~
      subroutine getb(b1,b2)
#include "aaa.h"
#include "ems.h"
      common/cranphi/ranphi
      common /cncl/xproj(mamx),yproj(mamx),zproj(mamx)
     *            ,xtarg(mamx),ytarg(mamx),ztarg(mamx)
      common/geom/rmproj,rmtarg,bmax,bkmx
      rmproj=1.19*maproj**(1./3.)-1.61*maproj**(-1./3.)+fctrmx*.54
      rmtarg=1.19*matarg**(1./3.)-1.61*matarg**(-1./3.)+fctrmx*.54
      b1=bminim
      b2=min(rmproj+rmtarg,bmaxim)
 10   bimevt=sqrt(b1**2+(b2**2-b1**2)*rangen())
      phievt=0
      ranphi=0
      bx=bimevt
      by=0
      call conxyz('p',mamx,xproj,yproj,zproj,ypjtl-yhaha)
      call conxyz('t',mamx,xtarg,ytarg,ztarg,yhaha)
      bkmx=conbmd()
      do i=1,maproj
      lproj(i)=0
      enddo
      do j=1,matarg
      ltarg(j)=0
      enddo
      koll=0
      do 12 i=1,maproj
      do 11 j=1,matarg
      bij=sqrt((xproj(i)+bx-xtarg(j))**2+(yproj(i)+by-ytarg(j))**2)
      if(bij.gt.bkmx)goto 11
      koll=koll+1
      lproj(i)=lproj(i)+1
      ltarg(j)=ltarg(j)+1
 11   continue
 12   continue
      if(koll.eq.0)goto10
      nglevt=koll
      npjevt=0
      ntgevt=0
      do i=1,maproj
      if(lproj(i).gt.0)npjevt=npjevt+1
      enddo
      do j=1,matarg
      if(ltarg(j).gt.0)ntgevt=ntgevt+1
      enddo
      end
      !~~~~~~~~
      subroutine getbcheck
#include "aaa.h"
      parameter (kbin=1000)
      real zbim(kbin)
      if(abs(b2-b1).lt.1)return
      write(ifmt,'(a,$)')
     .' check GetImpactParam ...  '
      call getb(b1,b2)
      do k=1,kbin
        zbim(k)=0
      enddo
      nsim=2000
      do i=1,nsim
        call getb(b1,b2)
        k=1+(bimevt-b1)/(b2-b1)*kbin
        if(k.ge.1.and.k.le.kbin)zbim(k)=zbim(k)+1./nsim
      enddo
      k=1
      sum=zbim(1)
      do while(sum.le.0.8)
        k=k+1
        b=b1+k*(b2-b1)/kbin
        sum=sum+zbim(k)
        zbim(k)=sum
      enddo
      write(ifmt,'(2f6.2,$)')b,bref80
      write(ifmt,'(a)')'   done'
      end

c------------------------------------------------------------------------------
      subroutine DefineEtaPlot
c------------------------------------------------------------------------------
#include "aaa.h"
#include "ho.h"
      common/cetaplot/etaplot(20),ietaplot(20),mxetaplot
      do ne=1,mxetaplot
        ieta=1
        etasoll=etaplot(ne)
        difmin=9999
        do neta=1,nzhy
          eta=etahy(neta)
          if(abs(eta-etasoll).lt.difmin)then
            difmin=abs(eta-etasoll)
            ieta=neta
          endif
        enddo
        ietaplot(ne)=ieta
      enddo
      end

c-------------------------------------------------------------------------
      subroutine getxyrange(neta,ntau)
c-------------------------------------------------------------------------
c      use ArrayModule, only: getd
c      use hocoModule, only: epsc
      double precision getHydynEpsc
      
#include "aaa.h"
#include "ho.h"
      common/crange/nxrange(nyhxx,2),nyrange(nxhxx,2)
      tol=0.001
      do nx=1,nxhy
       !---min
       ny=0
       ww=0.
       do while(ny.lt.nyhy.and.ww.lt.tol)
        ny=ny+1
        wwx=ww
        ww=0.
        do i=1,3
        ww=max(ww,abs(getHydynEpsc(neta,ntau,nx,ny)))
        enddo
       enddo
       if(ny.eq.nyhy.and.ww.eq.0)ny=nyhy
       !if(ny.gt.1)write(*,'(5x,2i5,2f10.4,$)')nx,ny,ww,wwx
       !if(ny.eq.1)write(*,'(5x,2i5,f10.4,10x,$)')nx,ny,ww
       nyrange(nx,1)=ny
       !---max
       ny=nyhy+1
       ww=0.
       do while(ny.gt.1.and.ww.lt.tol)
        ny=ny-1
        wwx=ww
        ww=0.
        do i=1,3
        ww=max(ww,abs(getHydynEpsc(neta,ntau,nx,ny)))
        enddo
       enddo
       if(ny.eq.1.and.ww.eq.0)ny=0
       !if(ny.lt.nyhy)write(*,'(2i5,2f10.4,$)')nx,ny,ww,wwx
       !if(ny.eq.nyhy)write(*,'(2i5,f10.4,10x,$)')nx,ny,ww
       nyrange(nx,2)=ny
       !write(ifmt,*)'                 ',nyrange(nx,1),nyrange(nx,2)
      enddo
      do ny=1,nyhy
       !---min
       nx=0
       ww=0.
       do while(nx.lt.nxhy.and.ww.lt.tol)
        nx=nx+1
        wwx=ww
        ww=0.
        do i=1,3
        ww=max(ww,abs(getHydynEpsc(neta,ntau,nx,ny)))
        enddo
       enddo
       if(nx.eq.nxhy.and.ww.eq.0)nx=nxhy
       !if(nx.gt.1)write(*,'(5x,2i5,2f10.4,$)')ny,nx,ww,wwx
       !if(nx.eq.1)write(*,'(5x,2i5,f10.4,10x,$)')ny,nx,ww
       nxrange(ny,1)=nx
       !---max
       nx=nxhy+1
       ww=0.
       do while(nx.gt.1.and.ww.lt.tol)
        nx=nx-1
        wwx=ww
        ww=0.
        do i=1,3
        ww=max(ww,abs(getHydynEpsc(neta,ntau,nx,ny)))
        enddo
       enddo
       if(nx.eq.1.and.ww.eq.0)nx=0
       !if(nx.lt.nxhy)write(*,'(2i5,2f10.4,$)')ny,nx,ww,wwx
       !if(nx.eq.nxhy)write(*,'(2i5,f10.4,10x,$)')ny,nx,ww
       nxrange(ny,2)=nx
       !write(ifmt,*)'                 ',nxrange(ny,1),nxrange(ny,2)
      enddo
      !if(neta.ne.7.or.ntau.ne.1)return
      !ny=22
      !write(ifmt,*)'------',neta,ntau,ny
      !write(ifmt,*)(epsc(neta,ntau,nx,ny),nx=1,41)
      end

c--------------------------------------------------------------------------------------
      subroutine PimEos(e, b1, b2, b3,   Tmp, chmB, chmQ, chmS, p)
c--------------------------------------------------------------------------------------
#include "aaa.h"
      double precision ee, nb, nq, ns, tt, mub, muq, mus, pp
      if(ioeos.eq.0.or.ioeos.eq.1)then
        Tmp=PiEos0(5,e,b1)
        chmB=PiEos0(6,e,b1)
        chmQ=0
        chmS=PiEos0(7,e,b1)
        p=PiEos0(4,e,b1)
      elseif(ioeos/10.eq.2.or.ioeos.eq.3)then
        ee=e
        nb=b1
        nq=b2
        ns=b3
        call eoshlle(ee, nb, nq, ns, tt, mub, muq, mus, pp)
        Tmp=tt
        chmB=mub
        chmQ=muq
        chmS=mus
        p=pp
      else
        stop'\n\n STOP in  PimEos, wrong choice of ioeos\n\n'
      endif
      end

c--------------------------------------------------------------------------------------
      function PiEos(k,e,b1,b2,b3)
c--------------------------------------------------------------------------------------
      ! k=2  entropy density
      !   4  pressure
      !   5  temperature
      !   6  baryon chemical potential
      !   7  strangeness chemical potential
      !   8  charge chemical potential
      ! e=energy density
      ! b1=baryon density
      ! b2=charge density
      ! b3=strangeness density
      !-----------------------------
#include "aaa.h"
#include "ho.h"
      double precision ee, nb, nq, ns, tt, mub, muq, mus, pp
      if(ioeos.ne.0)then
        ee=e
        nb=b1
        nq=b2
        ns=b3
        call eoshlle(ee, nb, nq, ns, tt, mub, muq, mus, pp)
        if(k.eq.2)w=0
        if(k.eq.4)w=pp
        if(k.eq.5)w=tt
        if(k.eq.6)w=mub
        if(k.eq.7)w=mus
        if(k.eq.8)w=muq
        PiEos=w
      else
        stop'STOP in  PiEos, ioeos is zero'
      endif
      end

c--------------------------------------------------------------------------------------
      function PiEos0(k,e,b)
c--------------------------------------------------------------------------------------
#include "ho.h"
      real wi(2),wj(2)

      n=2
      do while(n.lt.mxBeos.and.Beos(n).lt.b)
       n=n+1
      enddo
      n1=n-1
      n2=n
      frac=(b-Beos(n1))/(Beos(n2)-Beos(n1))
      wi(1)=1-frac
      wi(2)=frac
      i1=n1

      n=2
      do while(n.lt.mxEeos.and.Eeos(n).lt.e)
       n=n+1
      enddo
      !if(e.gt.oEeos)write(ifmt,*)'PiEos0: e = ',e,' out of range'
      n1=n-1
      n2=n
      frac=(e-Eeos(n1))/(Eeos(n2)-Eeos(n1))
      wj(1)=1-frac
      wj(2)=frac
      j1=n1

      w=0
      do i=1,2
      do j=1,2
      w=w+wi(i)*wj(j)*eost(k,i1+i-1,j1+j-1)
      enddo
      enddo

      PiEos0=w

      end

c--------------------------------------------------------------------------------------
      subroutine hynFOStore(id,iQ,ptr,pha,am,rap,rad,phi,tau,eta)
c--------------------------------------------------------------------------------------
#include "aaa.h"
#include "ho.h"
      common/cranphi/ranphi
      common/crap/nraphi,jfac
      !double precision dtime,utdble
      data ncntfoa/0/
      save ncntfoa
      ncntfoa=ncntfoa+1
      if(ncntfoa.eq.1)nptlb=nptl
      phinull=phievt+ranphi
      nptl=nptl+1
      if(nptl.gt.mxptl)
     . call utstop('hynFOStore: mxptl too small&')
      idptl(nptl)=id !->overwritten
      ibptl(1,nptl)=iQ
      pptl(1,nptl)=ptr*cos(pha+phinull)
      pptl(2,nptl)=ptr*sin(pha+phinull)
      pptl(3,nptl)=sqrt(am**2+ptr**2)*sinh(rap)
      pptl(4,nptl)=sqrt(am**2+ptr**2)*cosh(rap)
      pptl(5,nptl)=am
      ityptl(nptl)=61 !->overwritten
      !from now on we use idEPOS !!! 
      idepos=idtrafo('pdg','nxs',id)
      idptl(nptl)=idepos
      ityptl(nptl)=60
            !~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
            !betatx=pptl(1,nptl)/pptl(4,nptl)
            !betaty=pptl(2,nptl)/pptl(4,nptl)
            !betatz=pptl(3,nptl)/pptl(4,nptl)
            !------------------------------------
            ! |\vec{beta}|>1  happens due to limited accuracy (single prec)
            ! but is not really a problem, since the proper p(3) and p(4)
            ! may be reconstructred via a double precision calculation.
            !       first p -> pp (double), and then:
            !       aamt=sqrt(pp(5)**2+pp(1)**2+pp(2)**2)
            !       rrap=sign(1d0,pp(3))*log((pp(4)+abs(pp(3)))/aamt) 
            !       pp(3)=aamt*sinh(rrap)
            !       pp(4)=aamt*cosh(rrap)
            !------------------------------------
            !if((betatx**2+betaty**2+betatz**2).gt.1.0d0)then
            !  write(*,*)'NaN ; beta>1 in hynFOStore'
            !  write(*,*)'NaN    ',ityptl(nptl),idptl(nptl)
            !  write(*,*)'NaN    ', betatx,betaty,betatz
            !  write(*,*)'NaN    ',pptl(1,nptl),pptl(2,nptl)
            !.                 ,pptl(3,nptl) ,pptl(4,nptl)
            !endif
            !~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
      dtau=getHydynTauhy(2)-getHydynTauhy(1)
      taurd=tau-dtau/2+rangen()*dtau
      dleta=(etahy(2)-etahy(1))/jfac
      etard=eta-dleta/2+rangen()*dleta
      dphi=phihy(2)-phihy(1)
      phird=phi-dphi/2+rangen()*dphi
      xorptl(1,nptl)=rad*cos(phird+phinull)
      xorptl(2,nptl)=rad*sin(phird+phinull)
      xorptl(3,nptl)=taurd*sinh(etard)
      xorptl(4,nptl)=taurd*cosh(etard)
      iorptl(nptl)=0
      jorptl(nptl)=0
      istptl(nptl)=0
      ifrptl(1,nptl)=0
      ifrptl(2,nptl)=0
      tivptl(1,nptl)=xorptl(4,nptl)
      call idtau(idptl(nptl),pptl(4,nptl),pptl(5,nptl),taugm) !uses idEPOS
      r=rangen()
      tivptl(2,nptl)=tivptl(1,nptl)+taugm*(-alog(r))
      radptl(nptl)=0.
      dezptl(nptl)=0.
      itsptl(nptl)=0
      rinptl(nptl)=0
      zpaptl(1,nptl)=0.      !to be updated ???
      zpaptl(2,nptl)=0.      !energy density
      !WWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWW
      !  write(ifch,'(a,2i8,i4,5f10.2,4x,4f10.2,i4)')
      !. ' WWWWWW    hynFOStore   ',nptl,idptl(nptl),ityptl(nptl)
      !. ,sqrt(pptl(1,nptl)**2+pptl(2,nptl)**2)
      !. , xorptl(4,nptl),xorptl(3,nptl),xorptl(1,nptl),xorptl(2,nptl)
      !. ,tivptl(2,nptl),r,-alog(r),taugm
      !WWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWW
      end

c--------------------------------------------------------------------------------------
      function hynFOFactor2(nsurf,s,q)
c--------------------------------------------------------------------------------------
      real q(4) , s(4)
      fof  =    s(4) *q(4)      !0
     .        + s(1) *q(1)      !1 
     .        + s(2) *q(2)      !2
     .        + s(3) *q(3)      !3
      ii=mod(nsurf,2)*2-1
      hynFOFactor2=fof  * ii
      end

c--------------------------------------------------------------------------------------
      subroutine hynDefineParticles(iii)
c--------------------------------------------------------------------------------------
      !  chepoB(mch),chepoQ(mch),chepoS(mch),chot(mte,i,j,k,m)
      !  ffstat(mte,ii,i,j,k,m)
      !-----------------------------------------
#include "aaa.h"
#include "ho.h"
      double precision  gammaS, B, mS, delta0, muc, volex0,
     . emax, e0, anmax, bnmax, an0, muBmax, muQmax, muSmax, Taccuracy,
     . muBmx, muQmx, muSmx,
     . aaa, bbb, temmax
      common /ceospar/  gammaS, B, mS, delta0, muc, volex0,
     . emax, e0, anmax, bnmax, an0, muBmax, muQmax, muSmax, Taccuracy,
     . muBmx, muQmx, muSmx,
     . aaa, bbb,  temmax,  nxe, nxn, maxIter
      character dname*30
      common /camaq/amaq(0:3)
      logical hirpce
      double precision das
      common/ciuskip/iuskip
      data istatgetp /0/
      save istatgetp
      common/cistat/istateos,istatptl

#if __ROOT__
      if(iuskip.eq.2)
     .stop'\n\n hynDefineParticles: iuskip.ne.2 required\n\n'
      kspes=1
      if(istatptl.eq.0)then
        write(ifmt,'(a)')'load ptl and dky table ... '
        call loadptldky(
     .   fn3p(1:nfn3p)//CHAR(0)
     .  ,fn3d(1:nfn3d)//CHAR(0))
        istatptl=1
      else
        write(ifmt,'(a)')'ptl and dky table already loaded'
      endif
      nspesxx=igetnparticles()
      istatgetp=istatgetp+1
      if(istatgetp.eq.1)then !~~~~~~~
        m=0
        last=0
        do
          call igetparticle(dname,das,iQ,id,igs,js,iB,iS,last,iflag)  !----->use getNextParticle(...idEPOS...) + trafo idPDG
          !iflag=1 : 54 "basic" hadrons
          !iflag=2 : charm=bottom=top=0 
          !iflag=3 : all
          !---------------------------- 
          as=das
          gs=igs
          if(abs(last).eq.1)then
            m=m+1
            !idepos=idtrafo('pdg','nxs',id)
            !print*,m,'  ',dname(1:22),'  ',idepos,'  iflag',iflag
            if(m.gt.mspes)then
              print*,'####### m mspes :',m,mspes 
              stop'####### ERROR in hynDefineParticles (2) #######'
            endif
            aspes(1,m)=as
            iQspes(m)=iQ
            ispes(m)=id
            jspes(m)=js
            gspes(m)=gs
            iBspes(m)=iB
            iSspes(m)=iS
          endif
          if(last.lt.0)goto1
        enddo
   1    continue
        nspes=m
      endif !~~~~~~
      write(ifmt,'(2a,i3,a,i4)')'load ptl and dky table done.'
     ., ' Use',kspes,' to ',nspes
      
      if(ioclude.eq.5)then
        stop'ERROR hynDefineParticles: iclude=5 not working currently'
        !do m=kspes,nspes
        !  aspes(2,m)=0
        !enddo
      endif

      if(iii.eq.2)return

      if(kfrout.ge.3)then
        write(ifmt,'(a)')'hynDefineParticles -> skip (Testrun)'
        return
      endif

      write(ifmt,'(a,2f7.4,a,$)')
     . 'hynDefineParticles FO T range',temin,temax,' ... '

      istat=1  !0=Boltzmann, 1=Bose/Fermi
      pi=3.1415927
      hbar=0.197327

      call gaulag(xlag,wlag,nlag,0.)
      do n=1,nlag
      wlag(n)=wlag(n)*exp(xlag(n))
      enddo

      do mte=1,mtemp
        temfo(mte)=temin+(mte-1.)/(mtemp-1.)*(temax-temin)
      enddo
      chemx=muBmx
      do mch=1,mchem
        chepoB(mch)=-chemx+(mch-1.)/(mchem-1.)*2*chemx
      enddo
      chemx=muQmx
      do mch=1,mchem
        chepoQ(mch)=-chemx+(mch-1.)/(mchem-1.)*2*chemx
      enddo
      chemx=muSmx
      do mch=1,mchem
        chepoS(mch)=-chemx+(mch-1.)/(mchem-1.)*2*chemx
      enddo

      do m=kspes,nspes
        id=ispes(m)
        cohi(m)=0
        if(hirpce)cohi(m)=cothi(id)
        do mte=1,mtemp
        tefo=temfo(mte)
        do i=1,mchem
        do j=1,mchem
        do k=1,mchem
          ch=iBspes(m)*chepoB(i)+iQspes(m)*chepoQ(j)+iSspes(m)*chepoS(k)
          chot(mte,i,j,k,m)=ch/tefo
        enddo
        enddo
        enddo
        enddo
      enddo

      facphase=1./(2*pi*hbar)**3
      mit=int(mchem/2.)+1
      do mte=1,mtemp
      tefo=temfo(mte)
      do i=1,mchem
      do j=1,mchem
      do k=1,mchem
      !write(*,'(/25x,a,$)')'yie:'
       do m=0,kspes-1
        ffstat(mte,1,i,j,k,m)=0
        ffstat(mte,2,i,j,k,m)=0
       enddo
       eesum=0
       hhsum=0
       do m=kspes,nspes
        id=ispes(m)
        am=aspes(1,m)
        amx=aspes(2,m)
        esum=0
        fsum=0
        gsum=0
        hsum=0
        do n=1,nlag
          x=xlag(n)  ! p/tefo
          e=sqrt(am**2+x**2*tefo**2)
          w=exp(-sqrt(am**2/tefo**2+x**2)) 
     .          * exp(chot(mte,i,j,k,m)+cohi(m))
          fsum=fsum+wlag(n)*x**2*w /(1+istat*jspes(m)*w)
          esum=esum+wlag(n)*x**2*w /(1+istat*jspes(m)*w) *e
          wx=exp(-sqrt(amx**2/tefo**2+x**2))
     .          * exp(chot(mte,i,j,k,m)+cohi(m))
          gsum=gsum+wlag(n)*x**2*wx /(1+istat*jspes(m)*wx)
          hsum=hsum+wlag(n)*x**2*wx /(1+istat*jspes(m)*wx) *e
        enddo
        esum=esum * facphase * gspes(m) * 4 * pi * tefo**2  * tefo
        fsum=fsum * facphase * gspes(m) * 4 * pi * tefo**2  * tefo
        gsum=gsum * facphase * gspes(m) * 4 * pi * tefo**2  * tefo
        hsum=hsum * facphase * gspes(m) * 4 * pi * tefo**2  * tefo
        ffstat(mte,1,i,j,k,m)=ffstat(mte,1,i,j,k,m-1)+fsum
         !if(i.eq.mit.and.j.eq.mit.and.k.eq.mit)print*
         !.,'statistical weight   ',m,ispes(m),fsum,ffstat(mte,1,i,j,k,m)
        ffstat(mte,2,i,j,k,m)=ffstat(mte,2,i,j,k,m-1)+gsum
        eesum=eesum+esum
        hhsum=hhsum+hsum
        !if(m.lt.6)write(*,'(i5,a,f8.5,$)')id,':',fsum
       enddo
       ffstat(mte,1,i,j,k,nspes+1)=eesum
       ffstat(mte,1,i,j,k,nspes+2)=1
       ffstat(mte,2,i,j,k,nspes+1)=hhsum
       ffstat(mte,2,i,j,k,nspes+2)=eesum/hhsum
       !write(ifmt,'(f5.2,$)')ffstat(mte,1,i,j,k,nspes)
      enddo
      enddo
      enddo
      enddo
      write(ifmt,'(a)')'  done'
#else
        stop "ROOT needed for hynDefineParticles"
#endif

      end

c--------------------------------------------------------------------------------------
      subroutine getQBM(idpdg,iQ,iB,aM)
c--------------------------------------------------------------------------------------
      !returns charge iQ and baryon number iB for particle idpdg
      !-----------------------------------------------------------
#include "aaa.h"
#include "ho.h"
      common/cistat/istateos,istatptl
      if(istatptl.eq.0)call hynDefineParticles(2)
      do m=1,mspes
        if(ispes(m).eq.idpdg)then
        iQ=iQspes(m)
        iB=iBspes(m)
        aM=aspes(1,m)
        return
        endif
      enddo
      write(ifmt,*)'******* idpdg not found '
      write(ifmt,*)'*******    idpdg = ',idpdg
      write(ifmt,*)'******* list of particles (ispes iQspes iBspes):'   
      do m=1,mspes
        write(ifmt,*)'******* ',m,'   ',ispes(m),iQspes(m),iBspes(m)
      enddo
      write(ifmt,*)'******* STOP  ERROR 18072014 **********'
      stop'###### ERROR 18072014 ######'
      end  

c---------------------------------------------------------------------------------
      function hynRanBoseFermi(am,cot,js,tfo,istat)
c---------------------------------------------------------------------------------
c generates    x**2*exp(b-sqrt(x**2+a**2))            x**2
c x=p/tfo acc  --------------------------- = ------------------------
cto Bose/Fermi  1+-exp(b-sqrt(x**2+a**2))    exp(sqrt(x**2+a**2)-b)+1
c-----------------------------------------------------------------------
#include "aaa.h"
      a=am/tfo
      b=cot
      if(istat.eq.0)then                             !Boltznann
        x=hynRanTherm(a)
      elseif(js.eq.1.and.b.gt.a)then                  !Fermi
        !proposal:
        !------------------------
        !     x**2
        !   ----------
        !   exp(x-b)+1
        !------------------------
        !generate t=x-b
        cr1=1.+3./b+6./b**2+6./b**3
        cr2=3./b
        cr3=3./b+6./b**2
   3    zuk=rangen()*cr1-1.
        if(zuk.le.0.)then
          tt=b*(rangen()**.3333-1.)
        elseif(zuk.le.cr2 )then
          tt=-log(rangen())
        elseif(zuk.lt.cr3 )then
          tt=-log(rangen())-log(rangen())
        else
          tt=-log(rangen())-log(rangen())-log(rangen())
        endif
        if(rangen().gt.1./(1.+exp(-abs(tt))))goto 3
        x=tt+b
        pacc=(exp(x-b)+1)/(exp(sqrt(x**2+a**2)-b)+1)
        if(rangen().gt.pacc)goto3
      elseif(js.eq.1)then                            !Fermi
  2     x=hynRanTherm(a)
        p=x*tfo
        e=sqrt(p**2+am**2)
        w=exp(b-e/tfo)
        pacc=1./(1.+w)
        if(rangen().gt.pacc)goto2
      elseif(js.eq.-1.and.b.lt.a)then                !Bose
  1     x=hynRanTherm(a)
        p=x*tfo
        e=sqrt(p**2+am**2)
        w=exp(b-e/tfo)
        w0=exp(b-am/tfo)
        pacc=(1.-w0)/(1.-w)
        if(rangen().gt.pacc)goto1
      elseif(js.eq.-1)then                           !Bose
        x=sqrt(b**2-a**2)
        !write(ifmt,*)'***** BE condensate *****'
      else
        stop'in hynRanBoseFermi: unknown statistics (080726) '
      endif
      hynRanBoseFermi=x
      end

c------------------------------------------------------------------------------
      subroutine xhynRanBoseFermi
c------------------------------------------------------------------------------
      parameter(kkmax=100)
      common/ccc/count(kkmax),dx
      call ranfini(0d0,0,-1) !initialize some parameters
      call ranfini(1234456d0,2,1)
      open(55,file='zz.histo')
      write(55,*)'newpage  zone 4 4 1 '
      dx=1
      xmax=40
      kmax=xmax/dx

      do js=-1,1,2

      do kk=1,16
      b=(kk-1)*2
      do k=1,kmax
      count(k)=0
      enddo
      nsim=10000
      do i=1,nsim
      x=hynRanBoseFermi(1.6,b,js,.158,1)
      ix=1+x/dx
      if(ix.ge.1.and.ix.le.kmax)count(ix)=count(ix)+1
      enddo
      do k=1,kmax
      count(k)=count(k)/nsim/dx
      enddo
      write(55,*)'openhisto '
      write(55,*)'xrange 0  ',xmax,' yrange 1e-10 1e-1 '
      write(55,*)'htyp lin xmod lin ymod log'
      write(55,*)'txt  "xaxis x "'
      write(55,*)'txt  "yaxis y "'
      write(55,'(a,f7.2,a)')'text 0.6 0.85 "b = ',b,' "'
      write(55,*)'array 2'
      do k=1,kmax
      x=dx*(k-0.5)
      write(55,*)x,count(k)/x**2
      enddo
      write(55,*)'endarray closehisto  plot 0-'
      qua=0
      dl=xmax/200
      do i=0,200
       fac=2
       if(i.eq.0)fac=1
       if(i.eq.200)fac=1
       x=i*dl
       ww=exp( sqrt(x**2+10.0**2)-b )
       if(ww+js.gt.0.)qua=qua+x**2/( ww+js ) *fac*dl*0.5
      enddo
      v=1/qua
      if(js.eq.1)then
      write(55,*)'openhisto '
      write(55,*)'xrange 0  ',xmax,' yrange 1e-10 1e-1  '
      write(55,*)'htyp lin xmod lin ymod log'
      write(55,*)'txt  "xaxis x "'
      write(55,*)'array 2'
      do k=1,kmax
      x=dx*(k-0.5)
      write(55,*)x,v/(exp(x-b)+1)
      enddo
      write(55,*)'endarray closehisto  plot 0-'
      endif
      write(55,*)'openhisto '
      write(55,*)'xrange 0  ',xmax,' yrange 1e-10 1e-1  '
      write(55,*)'htyp lin xmod lin ymod log'
      write(55,*)'txt  "xaxis x "'
      write(55,*)'txt  "yaxis y "'
      write(55,*)'array 2'
      do k=1,kmax
      x=dx*(k-0.5)
      ww=exp(-b+sqrt(x**2+10.0**2))
      if(ww+js.gt.0.)write(55,*)x,v/(ww+js)
      enddo
      write(55,*)'endarray closehisto  plot 0'
      enddo

      enddo

      stop
      end

c------------------------------------------------------------------------------
      function hynRanTherm(a)
c------------------------------------------------------------------------------
c   generates a random number according to f(x) ~ x**2*exp(-sqrt(x**2+a**2))
c   in the interval from zero to infinity
c------------------------------------------------------------------------------
      !ntry=0
  1   i=2
      if(rangen().le.a**3/(a**3+3*a**2+6*a+6))i=1
      !ntry=ntry+1
      if(i.eq.1)then
        x=a*rangen()**(1./3.)
      if(rangen().gt.exp(a-sqrt(x**2+a**2)))goto1
      elseif(i.eq.2)then
        !f(z)~a**2+2*a*z+z**2)*exp(-z)  from zero to infty
        r=rangen()
        if(r.lt.a**2/(a**2+2*a+2))then
           z=-log(rangen())
        elseif(r.lt.(a**2+2*a)/(a**2+2*a+2))then
           z=-log(rangen())-log(rangen())
        else
           z=-log(rangen())-log(rangen())-log(rangen())
        endif
        x=a+z
      if(rangen().gt.exp(x-sqrt(x**2+a**2)))goto1
      endif
      hynRanTherm=x
      end

c---------------------------------------------------------------------------------
      subroutine hynDefineRapScale
c---------------------------------------------------------------------------------
#include "aaa.h"
#include "ho.h"
      common/crap/nraphi,jfac /crapi/delrapi
      delrapi=0.2
      deta=etahy(2)-etahy(1)
      write(ifmt,'(a,$)')'hynDefineRapScale  ...'
      jfac=1
      delrap=deta
      do while(delrap.gt.delrapi)
        jfac=jfac+1
        delrap=deta/jfac
      enddo
      write(ifmt,'(f7.2,i3,f7.2,5x,$)')deta,jfac,delrap
      nraphi=5
      rapmn=-delrap*nraphi
      do while(rapmn.gt.-2.)
        nraphi=nraphi+1
        if(nraphi.gt.20)stop'(1539070608)   '
        rapmn=-delrap*nraphi
      enddo
      write(ifmt,'(i4,$)')nzhy
      write(ifmt,'(a)')'  done'
      end

c-----------------------------------------------------------------------
      subroutine makeTestInvTable(fun,e0,emax,nxe, n0,nmax,nxn)
c-----------------------------------------------------------------------
#include "aaa.h"
      character*80 fn
      external fun
      integer ixe,ixnb,ixnq,ixns,ii
      double precision xe_max,xn_max
      double precision  e0,emax,n0,nmax,  xe,xnb,xnq,xns
      double precision p,T,mub,muq,mus,  e, nb, nq, ns
      xe_max = log(emax/e0+1)
      xn_max = log(nmax/n0+1)
      ii=index(fnhi(1:nfnhi),".histo")-1
      fn=fnhi(1:ii)//".eos3f"
      write(ifmt,'(2a)')'write eos table into ',fn(1:ii+4)
      call clop(3)
      open(15,file=fn(1:ii+4),status='unknown')
      write(15,*)emax,e0,nmax,n0,nxe,nxn
      do  ixe=0,nxe-1
      xe = ixe*xe_max/(nxe-1)
      e = e0*(exp(xe)-1.)
      do ixnb=0,nxn-1
      do ixnq=0,nxn-1
      do ixns=0,nxn-1
        xnb = -xn_max + ixnb*2.*xn_max/(nxn-1)
        xnq = -xn_max + ixnq*2.*xn_max/(nxn-1)
        xns = -xn_max + ixns*2.*xn_max/(nxn-1)
        nb = xnb/abs(xnb+1e-20)*n0*(exp(abs(xnb))-1.)
        nq = xnq/abs(xnq+1e-20)*n0*(exp(abs(xnq))-1.)
        ns = xns/abs(xns+1e-20)*n0*(exp(abs(xns))-1.)
        call fun(e,nb,nq,ns,p,t,mub,muq,mus)
        write(15,*)p,T,mub,muq,mus
        call clop(3)
        if(ixnb.eq.nxn/2.and.ixnq.eq.nxn/2.and.ixns.eq.nxn/2)then
          write(ifmt,*)'e bin ',ixe+1,'  e = ',e
     .       ,'  p = ',p,'  t = ',t
        endif
      enddo
      enddo
      enddo
      enddo
      close(15)
      end


cooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooo
cooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooo
cooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooo
cooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooo
cooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooo
cooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooo


c-----------------------------------------------------------------------
      subroutine ooHoEpsilon(itau1,itau2,dtaupl)
c-----------------------------------------------------------------------
c      use ArrayModule, only: getd
c      use hocoModule, only: epsc, barc
c      use hocoModule, only: barc
      double precision getHydynEpsc, getHydynBarc
#include "aaa.h"
#include "ho.h"
      character ctau*5, cbim*3, czz*5, txt*32
      character*11  bbtxt(2)
      character*16  bbt16(2)
      common/jcen/jcentr,jmxcentr
      real bb(2,nxhy,nyhy),bbmax(2)
      modu=nint(abs(dtaupl)/(getHydynTauhy(2)-getHydynTauhy(1)))
      if( nint(fdtau) .eq. 0 )then   
        itau2x=itau1
      elseif( fdtau .gt . 0. )then 
        itau2x=itau2
      else
        itau2x=(ntauhy-1)/modu+1
      endif 

      do itau=itau1,itau2x !~~~~~~~~~~~~itau loop~~~~~~~~~~~~~~

      nxoff=4  !5
      nyoff=4  !5
      delx5=(xmaxhy-xminhy)/(nxhy-1)*0.5
      dely5=(ymaxhy-yminhy)/(nyhy-1)*0.5
      ntau=1+modu*(itau-1) 
      tau=getHydynTauhy(ntau)

      if(tau.le.tauzer+tauup)then !~~~~~~~~~~~~tauup~~~~~~~~~~~~

      write(ctau,'(f5.2)')tau
      write(cbim,'(a,i2)')'J',jcentr
      txt='energy density [GeV/fm^{3}]  '
            
      do j=1,1   !3   !~~~~~~~~~~~~j loop (z) ~~~~~~~~~~~~
               
      iz=nzhy/2+1 + (j-1)*3
      deleta5=(zmaxhy-zminhy)/(nzhy-1)*0.5
      zz=zminhy+(iz-1)*deleta5*2
      if(zz.ge.0.)write(czz,'(f4.1)')zz
      if(zz.lt.0.)write(czz,'(f5.1)')zz
      call RootCanvas(-2)
      call histo2BeginFigure()
      call histo2BeginPlot(1, 1)
      call RootHisto(2,
     * txt//'(#eta_{s}='//czz//', #tau ='//ctau//' fm/c)   '//cbim//';'
     *,nxhy-2*nxoff
     *,xminhy-delx5+nxoff*delx5*2
     *,xmaxhy+delx5-nxoff*delx5*2
     *,nyhy-2*nyoff
     *,yminhy-dely5+nyoff*dely5*2
     *,ymaxhy+dely5-nyoff*dely5*2)
      call histo2BeginSubPlot(
     * 'src/KW/freeze.f:ooHoEpsilon(itau1,itau2,dtaupl)'
     *,'ooHoEpsilon')
      call histo2AddHeader(
     * txt//'($\eta_{s}$ ='//czz//', $\tau$ ='//ctau//' fm/c) '//cbim
     *,xminhy-delx5+nxoff*delx5*2
     *,xmaxhy+delx5-nxoff*delx5*2
     *,delx5*2
     *,yminhy-dely5+nyoff*dely5*2
     *,ymaxhy+dely5-nyoff*dely5*2
     *,dely5*2
     *,'x [fm]'
     *,'y [fm]')
      bbmax(1)=0. 
      bbmax(2)=0. 
      ccmax=0 
      do ix=1,nxhy
         xx=xminhy+(ix-1)*delx5*2
         do iy=1,nyhy
            yy=yminhy+(iy-1)*dely5*2
            eps=getHydynEpsc(iz,ntau,ix,iy)
            b1=getHydynBarc(1,iz,ntau,ix,iy)
            b2=getHydynBarc(2,iz,ntau,ix,iy)
            b3=getHydynBarc(3,iz,ntau,ix,iy)
            tem=PiEos(5,eps,b1,b2,b3)
            !if(tem.gt.tfrout)
            call RootFill(2,xx,yy,eps)
            call histo2FillArray(2,xx,yy,eps)
            bb(1,ix,iy)=eps
            bb(2,ix,iy)=b2
           !bb(2,ix,iy)=tem*1000
            bbtxt(1)='epsilon    ' 
            bbtxt(2)='charge dens' 
           !bbtxt(2)='temperature' 
            bbmax(1)=max(bbmax(1),eps)
            bbmax(2)=max(bbmax(2),abs(b2))
         enddo
      enddo
      call RootDraw('colz;','x [fm];','y [fm];')  !surf1 cont4Z
      call histo2EndSubPlot()   !surf1 cont4Z
      call histo2EndPlot()

      enddo !~~~~~~~~~~~~j loop (z) ~~~~~~~~~~~~

      iprint2d=1
      if(iprint2d.eq.1)then
        bbt16(1)='                '
        bbt16(2)='(*** means -100)'  
        n1=nxhy/2-19  !9
        n2=n1+40   !20
        do m=1,2 !~~~~~~~~~~~~m loop~~~~~~~~~~~~
        bbmaxtot=0 
        xx=-9999
        yy=-9999
        zz=-9999
        do iz=1,nzhy
         do ix=1,nxhy
          do iy=1,nyhy
            if(m.eq.1)bbx=    getHydynEpsc(  iz,ntau,ix,iy)
            if(m.eq.2)bbx=abs(getHydynBarc(2,iz,ntau,ix,iy))
            if(bbx.gt.bbmaxtot)then
              bbmaxtot=bbx
              xx=xminhy+(ix-1)*delx5*2
              yy=yminhy+(iy-1)*dely5*2
              zz=zminhy+(iz-1)*deleta5*2
            endif
          enddo
         enddo
        enddo
        write(ifch,*)'=== print2d (f.f) === 100 * ',bbtxt(m),' / '
     .  ,bbmax(m),'(max at z=0)    === tau = ',tau,'       ',bbt16(m) 
        write(ifch,*)'Total maximum (3D):', bbmaxtot
     .,'  Position:',xx,yy,zz
        write(ifch,'(4x,21i6)')
     .  (nint((xminhy+(i-1)*delx5*2)*100),i=n1,n2,2)
        do k=n2,n1,-1
          y=yminhy+(k-1)*dely5*2
          z=0
          if(bbmax(m).ne.0.)z=bb(m,i,k)/bbmax(m)
          write(ifch,'(i5,2x,41i3)')
     .    nint(y*100),(nint(z*100),i=n1,n2)
        enddo        
        enddo !~~~~~~~~~~~~m loop~~~~~~~~~~~~
      endif

      endif !~~~~~~~~~~~~tauup~~~~~~~~~~~~

      enddo !~~~~~~~~~~~~itau loop~~~~~~~~~~~~~~

      end

c-----------------------------------------------------------------------
      subroutine ooHoRadVel(itau1,itau2,dtaupl)
c-----------------------------------------------------------------------
c      use ArrayModule, only: getd
c      use hocoModule, only: epsc, velc, barc
c      use hocoModule, only: velc, barc
      double precision getHydynEpsc, getHydynBarc
      double precision getHydynVelc
#include "aaa.h"
#include "ho.h"
      
      character ctau*5, cbim*3, czz*5, txt*27
      common/jcen/jcentr,jmxcentr
      do itau=itau1,itau2
      nxoff=4  !5
      nyoff=4  !5
      delx5=(xmaxhy-xminhy)/(nxhy-1)*0.5
      dely5=(ymaxhy-yminhy)/(nyhy-1)*0.5
      modu=nint(dtaupl/(getHydynTauhy(2)-getHydynTauhy(1)))
      ntau=1+modu*(itau-1)
      tau=getHydynTauhy(ntau)
      if(tau.le.tauzer+tauup)then
      write(ctau,'(f5.2)')tau
      write(cbim,'(a,i2)')'J',jcentr
      txt='rad velocity [% of c]      '
      call histo2BeginFigure()
      call histo2BeginPlot(1, 1)
      do j=1,1   !2
      iz=nzhy/2+1 + (j-1)*3
      deleta5=(zmaxhy-zminhy)/(nzhy-1)*0.5
      zz=zminhy+(iz-1)*deleta5*2
      if(zz.ge.0.)write(czz,'(f4.1)')zz
      if(zz.lt.0.)write(czz,'(f5.1)')zz
      call RootCanvas(-2)
      call RootHisto(2,
     * txt//'(#eta_{s}='//czz//', #tau ='//ctau//' fm/c)   '//cbim//';'
     *,nxhy-2*nxoff
     *,xminhy-delx5+nxoff*delx5*2
     *,xmaxhy+delx5-nxoff*delx5*2
     *,nyhy-2*nyoff
     *,yminhy-dely5+nyoff*dely5*2
     *,ymaxhy+dely5-nyoff*dely5*2)
      call histo2BeginSubPlot(
     * 'srd/KW/freeze.f: subroutine ooHoRadVel(itau1,itau2,dtaupl)'
     *,'ooHoRadVel')
      call histo2AddHeader(
     * txt//'($\eta_{s}$ ='//czz//', $\tau$ ='//ctau//' fm/c) '//cbim
     *,xminhy-delx5+nxoff*delx5*2
     *,xmaxhy+delx5-nxoff*delx5*2
     *,2*delx5
     *,yminhy-dely5+nyoff*dely5*2
     *,ymaxhy+dely5-nyoff*dely5*2
     *,2*dely5
     *,'x [fm];'
     *,'y [fm];')

      do ix=1,nxhy
         xx=xminhy+(ix-1)*delx5*2
         do iy=1,nyhy
            yy=yminhy+(iy-1)*dely5*2
            vx=getHydynVelc(1,iz,ntau,ix,iy)
            vy=getHydynVelc(2,iz,ntau,ix,iy)
            r=sqrt(xx**2+yy**2)
            phi=0
            if(r.gt.0.)phi=sign(1.,yy)*acos(xx/r)
            vtg=vx*sin(phi)-vy*cos(phi)
            vrd=vx*cos(phi)+vy*sin(phi)
            eps=getHydynEpsc(  iz,ntau,ix,iy)
            b1= getHydynBarc(1,iz,ntau,ix,iy)
            b2= getHydynBarc(2,iz,ntau,ix,iy)
            b3= getHydynBarc(3,iz,ntau,ix,iy)
            tem=PiEos(5,eps,b1,b2,b3)
            w=vrd*100
            !w=sqrt(vx**2+vy**2)*100
            if(tem.gt.tfrout.and.w.gt.1e-4) then
               call RootFill(2,xx,yy,w)
               call histo2FillArray(2,xx,yy,w)
            endif
         enddo
      enddo
      call RootDraw('cont4Z;','x [fm];','y [fm];')  !surf1
      call histo2EndSubPlot()   !surf1 cont4Z
      enddo
      call histo2EndPlot()
      endif
      enddo
      end

c-----------------------------------------------------------------------
      subroutine ooHoTanVel(itau1,itau2,dtaupl)
c-----------------------------------------------------------------------
c      use ArrayModule, only: getd
c      use hocoModule, only: epsc, velc, barc
c      use hocoModule, only: velc, barc
      double precision getHydynEpsc, getHydynBarc
      double precision getHydynVelc
#include "aaa.h"
#include "ho.h"
      character ctau*5, cbim*3, czz*5, txt*27
      common/jcen/jcentr,jmxcentr
      do itau=itau1,itau2
      nxoff=0  !5
      nyoff=0  !5
      delx5=(xmaxhy-xminhy)/(nxhy-1)*0.5
      dely5=(ymaxhy-yminhy)/(nyhy-1)*0.5
      modu=nint(dtaupl/(getHydynTauhy(2)-getHydynTauhy(1)))
      ntau=1+modu*(itau-1)
      tau=getHydynTauhy(ntau)
      if(tau.le.tauzer+tauup)then
      write(ctau,'(f5.2)')tau
      write(cbim,'(a,i2)')'J',jcentr
      txt='tang velocity [% of c]     '
      do j=1,2
      iz=nzhy/2+1 + (j-1)*3
      deleta5=(zmaxhy-zminhy)/(nzhy-1)*0.5
      zz=zminhy+(iz-1)*deleta5*2
      if(zz.ge.0.)write(czz,'(f4.1)')zz
      if(zz.lt.0.)write(czz,'(f5.1)')zz
      call RootCanvas(-2)
      call histo2BeginFigure()
      call histo2BeginPlot(1, 1)
      call RootHisto(2,
     * txt//'(#eta_{s}='//czz//', #tau ='//ctau//' fm/c)   '//cbim//';'
     *,nxhy-2*nxoff
     *,xminhy-delx5+nxoff*delx5*2
     *,xmaxhy+delx5-nxoff*delx5*2
     *,nyhy-2*nyoff
     *,yminhy-dely5+nyoff*dely5*2
     *,ymaxhy+dely5-nyoff*dely5*2)
      call histo2BeginSubPlot(
     * 'src/KW/freeze.f:ooHoTanVel(itau1,itau2,dtaupl)'
     *,'ooHoTanVel')
      call histo2AddHeader(
     * txt//'($\eta_{s}$ ='//czz//', $\tau$ ='//ctau//' fm/c) '//cbim
     *,xminhy-delx5+nxoff*delx5*2
     *,xmaxhy+delx5-nxoff*delx5*2
     *,2*delx5
     *,yminhy-dely5+nyoff*dely5*2
     *,ymaxhy+dely5-nyoff*dely5*2
     *,2*dely5
     *,'x [fm]'
     *,'y [fm]')
      do ix=1,nxhy
         xx=xminhy+(ix-1)*delx5*2
         do iy=1,nyhy
            yy=yminhy+(iy-1)*dely5*2
            vx=getHydynVelc(1,iz,ntau,ix,iy)
            vy=getHydynVelc(2,iz,ntau,ix,iy)
            r=sqrt(xx**2+yy**2)
            phi=0
            if(r.gt.0.)phi=sign(1.,yy)*acos(xx/r)
            vtg=vx*sin(phi)-vy*cos(phi)
            vrd=vx*cos(phi)+vy*sin(phi)
            eps=getHydynEpsc(  iz,ntau,ix,iy)
            b1= getHydynBarc(1,iz,ntau,ix,iy)
            b2= getHydynBarc(2,iz,ntau,ix,iy)
            b3= getHydynBarc(3,iz,ntau,ix,iy)
            tem=PiEos(5,eps,b1,b2,b3)
            if(tem.gt.tfrout) then
               call RootFill(2,xx,yy,vtg*100)
               call histo2FillArray(2,xx,yy,vtg*100)
            endif
         enddo
      enddo
      call RootDraw('cont4Z;','x [fm];','y [fm];')  !surf1
      call histo2EndSubPlot()   !surf1 cont4Z
      call histo2EndPlot()
      enddo
      endif
      enddo
      end

c----------------------------------------------------------------------------
      subroutine oo2(ica,nyjjj,nzjjj
     . ,aryz,nybin,ymin,ymax,nzbin,zmin,zmax,txt1,txt2,txty,txtz)
c----------------------------------------------------------------------------

#include "aaa.h"
#include "ho.h"
      real aryz(nyjjj,nzjjj)
      character txt*29,txty*6,txtz*8,txt1*4,txt2*4
      txt='                             '
      txt(1:6)=txty
      txt(13:20)=txtz
      delz=(zmax-zmin)/(nzbin-1)
      dely=(ymax-ymin)/(nybin-1)
      call RootCanvas(ica)
      call histo2BeginFigure()
      call histo2BeginPlot(1, 1)
      call RootHisto(2,
     * txt//';'
     *,nybin
     *,ymin-dely/2
     *,ymax+dely/2
     *,nzbin
     *,zmin-delz/2
     *,zmax+delz/2
     *)
      call histo2BeginSubPlot(
     * 'src/KW/freeze.f:oo2(ica,nyjjj,nzjjj
     * ,aryz,nybin,ymin,ymax,nzbin,zmin,zmax,txt1,txt2,txty,txtz)'
     *,'oo2')
      call histo2AddHeader(
     * txt 
     *,ymin-dely/2
     *,ymax+dely/2
     *,dely
     *,zmin-delz/2
     *,zmax+delz/2
     *,delz
     *,'x [fm];'
     *,'y [fm];')
      do iy=1,nybin
      yy=ymin+(iy-1)*dely
      do iz=1,nzbin
      zz=zmin+(iz-1)*delz
          call RootFill(2,yy,zz,aryz(iy,iz))
          !print*,'+++++',iy,iz,yy,zz,aryz(iy,iz)
          call histo2FillArray(2,yy,zz,aryz(iy,iz))
      enddo
      enddo
      call RootDraw('colz;',txt1//';',txt2//';')  !surf1
      call histo2EndSubPlot() !surf1 cont4Z
      call histo2EndPlot()
      end

c----------------------------------------------------------------------------
      subroutine oo2bin(ica,nyjjj,nzjjj
     . ,aryz,nybin,ymin,ymax,nzbin,zmin,zmax,txt1,txt2,txty,txtz)
c----------------------------------------------------------------------------

#include "aaa.h"
#include "ho.h"
      real aryz(nyjjj,nzjjj)
      character txt*45,txty*10,txtz*35,txt1*4,txt2*4
      txt(1:10)=txty
      txt(11:45)=txtz
      delz=(zmax-zmin)/nzbin
      dely=(ymax-ymin)/nybin
      call RootCanvas(ica)
      call histo2BeginFigure()
      call histo2BeginPlot(1, 1)
      call RootHisto(2,
     * txt//';'
     *,nybin
     *,ymin
     *,ymax
     *,nzbin
     *,zmin
     *,zmax
     *)
      call histo2BeginSubPlot(
     * 'src/KW/freeze.f:oo2bin(ica,nyjjj,nzjjj
     * ,aryz,nybin,ymin,ymax,nzbin,zmin,zmax,txt1,txt2,txty,txtz)'
     *,'oo2bin')
      call histo2AddHeader(
     * txt
     *,ymin
     *,ymax
     *,dely
     *,zmin
     *,zmax
     *,delz
     *,txt1//';'
     *,txt2//';')
      do iy=1,nybin
      yy=ymin+(iy-0.5)*dely
      do iz=1,nzbin
      zz=zmin+(iz-0.5)*delz
          call RootFill(2,yy,zz,aryz(iy,iz))
          !print*,'+++++',iy,iz,yy,zz,aryz(iy,iz)
          call histo2FillArray(2,yy,zz,aryz(iy,iz))
      enddo
      enddo
      call RootDraw('lego2;',txt1//';',txt2//';')  !surf1
      call histo2EndSubPlot()  !surf1
      call histo2EndPlot()
      end

c----------------------------------------------------------------------------
      subroutine ooHoIniYEta(nyjjj,nzjjj
     . ,aryz,nybin,ymin,ymax,nzbin,zmin,zmax,txtx,txty)
c----------------------------------------------------------------------------

#include "aaa.h"
#include "ho.h"
      real aryz(nyjjj,nzjjj)
      character cbim*3, txt*29,txtx*6,txty*8
      common/jcen/jcentr,jmxcentr
      write(cbim,'(a,i2)')'J',jcentr
      txt='                             '
      txt(1:6)=txtx
      txt(13:20)=txty
      call RootCanvas(-3)
      call histo2BeginFigure()
      call RootHisto(2,
     * txt//cbim//';'
     *,nzbin
     *,zmin
     *,zmax
     *,nybin
     *,ymin
     *,ymax)
      call histo2BeginPlot(1, 1)
      call histo2BeginSubPlot(
     * 'src/KW/freeze.f: ooHoIniYEta(nyjjj,nzjjj
     * ,aryz,nybin,ymin,ymax,nzbin,zmin,zmax,txtx,txty)'
     *,'ooHoIniYEta')
      call histo2AddHeader(
     * txt//cbim
     *,zmin
     *,zmax
     *,(zmax-zmin)/nzbin
     *,ymin
     *,ymax
     *,(ymax-ymin)/nybin
     *,'x [fm];'
     *,'y [fm];')
      do iz=1,nzbin
         zz=zmin+(iz-0.5)*(zmax-zmin)/nzbin
         do iy=1,nybin
            yy=ymin+(iy-0.5)*(ymax-ymin)/nybin
            call RootFill(2,zz,yy,aryz(iy,iz))
            call histo2FillArray(2,zz,yy,aryz(iy,iz))
            !print*,'+++++',iy,iz,yy,zz,aryz(iy,iz)
         enddo
      enddo
      call RootDraw('lego2;','#eta_{s};','y [fm];')  !surf1
      call histo2EndSubPlot()  !surf1
      call histo2EndPlot()
      end

c-----------------------------------------------------------------------
      subroutine ooHoEpsilonEtaY(itau1,itau2,dtaupl)
c-----------------------------------------------------------------------
c      use ArrayModule, only: getd
c      use hocoModule, only: epsc, barc
c      use hocoModule, only: barc
      double precision getHydynEpsc, getHydynBarc

#include "aaa.h"
#include "ho.h"
      character ctau*5, cbim*3, txt*32
      common/jcen/jcentr,jmxcentr
      do itau=itau1,itau2
      netaoff=1  !4
      nyoff=5
      dely5=(ymaxhy-yminhy)/(nyhy-1)*0.5
      deleta5=(zmaxhy-zminhy)/(nzhy-1)*0.5
      modu=nint(dtaupl/(getHydynTauhy(2)-getHydynTauhy(1)))
      ntau=1+modu*(itau-1)
      tau=getHydynTauhy(ntau)
      if(tau.le.tauzer+tauup)then
      write(ctau,'(f5.2)')tau
      write(cbim,'(a,i2)')'J',jcentr
      txt='energy density [GeV/fm^{3}] '
      call RootCanvas(-3)
      call histo2BeginFigure()
      call histo2BeginPlot(1, 1)
      call RootHisto(2,
     * txt//'(x=0, tau='//ctau//' fm/c)   '//cbim//';'
     *,nzhy-2*netaoff
     *,zminhy-deleta5+netaoff*deleta5*2
     *,zmaxhy+deleta5-netaoff*deleta5*2
     *,nyhy-2*nyoff
     *,yminhy-dely5+nyoff*dely5*2
     *,ymaxhy+dely5-nyoff*dely5*2)
      call histo2BeginSubPlot(
     * 'srd/KW/freeze.f:ooHoEpsilonEtaY(itau1,itau2,dtaupl)'
     *,'ooHoEpsilonEtaY')
      call histo2AddHeader(
     * txt//'(x=0, $\tau$ ='//ctau//' fm/c)   '//cbim//';'
     *,zminhy-deleta5+netaoff*deleta5*2
     *,zmaxhy+deleta5-netaoff*deleta5*2
     *,deleta5*2
     *,yminhy-dely5+nyoff*dely5*2
     *,ymaxhy+dely5-nyoff*dely5*2
     *,dely5*2
     *,'$\eta_{s}$'
     *,'y [fm]')

      ix=nxhy/2+1
      do iz=1,nzhy
         zz=zminhy+(iz-1)*deleta5*2
         do iy=1,nyhy
            yy=yminhy+(iy-1)*dely5*2
            eps=getHydynEpsc(  iz,ntau,ix,iy)
            b1= getHydynBarc(1,iz,ntau,ix,iy)
            b2= getHydynBarc(2,iz,ntau,ix,iy)
            b3= getHydynBarc(3,iz,ntau,ix,iy)
            tem=PiEos(5,eps,b1,b2,b3)
            if(tem.gt.tfrout) then
               call RootFill(2,zz,yy,eps)
               call histo2FillArray(2,zz,yy,eps)
            endif
         enddo
      enddo
      call RootDraw('lego2;','#eta_{s};','y [fm];')  !surf1
      call histo2EndSubPlot() !surf1
      call histo2EndPlot()
      endif
      enddo
      end

c-----------------------------------------------------------------------
      subroutine ooHoVyEtaY(itau1,itau2,dtaupl)
c-----------------------------------------------------------------------
c      use ArrayModule, only: getd
c      use hocoModule, only: epsc, velc, barc
c      use hocoModule, only: velc, barc
      double precision getHydynEpsc, getHydynBarc
      
#include "aaa.h"
#include "ho.h"
      character ctau*5, cbim*3, txt*29
      common/jcen/jcentr,jmxcentr
      do itau=itau1,itau2
      netaoff=1  !4
      nyoff=5
      dely5=(ymaxhy-yminhy)/(nyhy-1)*0.5
      deleta5=(zmaxhy-zminhy)/(nzhy-1)*0.5
      modu=nint(dtaupl/(getHydynTauhy(2)-getHydynTauhy(1)))
      ntau=1+modu*(itau-1)
      tau=getHydynTauhy(ntau)
      if(tau.le.tauzer+tauup)then
      write(ctau,'(f5.2)')tau
      write(cbim,'(a,i2)')'J',jcentr
      txt='v_{y}  '
      call RootCanvas(-3)
      call histo2BeginFigure()
      call histo2BeginPlot(1, 1)
      call RootHisto(2,
     * txt//'(x=0, tau='//ctau//' fm/c)   '//cbim//';'
     *,nzhy-2*netaoff
     *,zminhy-deleta5+netaoff*deleta5*2
     *,zmaxhy+deleta5-netaoff*deleta5*2
     *,nyhy-2*nyoff
     *,yminhy-dely5+nyoff*dely5*2
     *,ymaxhy+dely5-nyoff*dely5*2)
      call histo2BeginSubPlot(
     * 'srd/KW/freeze.f:ooHoVyEtaY(itau1,itau2,dtaupl)'
     *,'ooHoVyEtaY')
      call histo2AddHeader(
     * txt//'(x=0, $\tau$ ='//ctau//' fm/c)   '//cbim
     *,zminhy-deleta5+netaoff*deleta5*2
     *,zmaxhy+deleta5-netaoff*deleta5*2
     *,deleta5*2
     *,yminhy-dely5+nyoff*dely5*2
     *,ymaxhy+dely5-nyoff*dely5*2
     *,dely5*2
     *,'$\eta_{s}'
     *,'y [fm]')
      ix=nxhy/2+1
      do iz=1,nzhy
         zz=zminhy+(iz-1)*deleta5*2
         do iy=1,nyhy
            yy=yminhy+(iy-1)*dely5*2
            vy= getHydynVelc(2,iz,ntau,ix,iy)
            eps=getHydynEpsc(  iz,ntau,ix,iy)
            b1= getHydynBarc(1,iz,ntau,ix,iy)
            b2= getHydynBarc(2,iz,ntau,ix,iy)
            b3= getHydynBarc(3,iz,ntau,ix,iy)
            tem=PiEos(5,eps,b1,b2,b3)
            if(tem.gt.tfrout) then
               call RootFill(2,zz,yy,vy)
               call histo2FillArray(2,zz,yy,vy)
            endif
         enddo
      enddo
      call RootDraw('lego2;','#eta_{s};','y [fm];')  !surf1
      call histo2EndSubPlot() !surf1
      call histo2EndPlot()
      endif
      enddo
      end

c-----------------------------------------------------------------------
      subroutine ooHoVzEtaY(itau1,itau2,dtaupl)
c-----------------------------------------------------------------------
c      use ArrayModule, only: getd
c      use hocoModule, only: epsc, velc, barc
c      use hocoModule, only: velc, barc
      double precision getHydynEpsc, getHydynBarc
      double precision getHydynVelc            
#include "aaa.h"
#include "ho.h"
      character ctau*5, cbim*3, txt*29
      common/jcen/jcentr,jmxcentr
      do itau=itau1,itau2
      netaoff=1  !4
      nyoff=5
      dely5=(ymaxhy-yminhy)/(nyhy-1)*0.5
      deleta5=(zmaxhy-zminhy)/(nzhy-1)*0.5
      modu=nint(dtaupl/(getHydynTauhy(2)-getHydynTauhy(1)))
      ntau=1+modu*(itau-1)
      tau=getHydynTauhy(ntau)
      if(tau.le.tauzer+tauup)then
      write(ctau,'(f5.2)')tau
      write(cbim,'(a,i2)')'J',jcentr
      txt='v_{z}  '
      call RootCanvas(-3)
      call histo2BeginFigure()
      call histo2BeginPlot(1, 1)
      call RootHisto(2,
     * txt//'(x=0, tau='//ctau//' fm/c)   '//cbim//';'
     *,nzhy-2*netaoff
     *,zminhy-deleta5+netaoff*deleta5*2
     *,zmaxhy+deleta5-netaoff*deleta5*2
     *,nyhy-2*nyoff
     *,yminhy-dely5+nyoff*dely5*2
     *,ymaxhy+dely5-nyoff*dely5*2)
      call histo2BeginSubPlot(
     * 'src/KW/freeze.f:ooHoVzEtaY(itau1,itau2,dtaupl)'
     *,'ooHoVzEtaY')
      call histo2AddHeader(
     * txt//'(x=0, $\tau$ ='//ctau//' fm/c) '//cbim
     *,zminhy-deleta5+netaoff*deleta5*2
     *,zmaxhy+deleta5-netaoff*deleta5*2
     *,deleta5*2
     *,yminhy-dely5+nyoff*dely5*2
     *,ymaxhy+dely5-nyoff*dely5*2
     *,dely5*2
     *,'$\eta_{s}$'
     *,'y [fm]')

      ix=nxhy/2+1
      do iz=1,nzhy
         zz=zminhy+(iz-1)*deleta5*2
         do iy=1,nyhy
            yy=yminhy+(iy-1)*dely5*2
            vz= getHydynVelc(3,iz,ntau,ix,iy)
            eps=getHydynEpsc(  iz,ntau,ix,iy)
            b1= getHydynBarc(1,iz,ntau,ix,iy)
            b2= getHydynBarc(2,iz,ntau,ix,iy)
            b3= getHydynBarc(3,iz,ntau,ix,iy)
            tem=PiEos(5,eps,b1,b2,b3)
            if(tem.gt.tfrout) then
               call RootFill(2,zz,yy,vz)
               call histo2FillArray(2,zz,yy,vz)
            endif
         enddo
      enddo
      call RootDraw('lego2;','#eta_{s};','y [fm];')  !surf1
      call histo2EndSubPlot() !surf1
      call histo2EndPlot()
      endif
      enddo
      end

c-----------------------------------------------------------------------
      subroutine ooHoEpsilonEtas(npad,itau1,itau2,dtaupl)
c-----------------------------------------------------------------------
c      use ArrayModule, only: getd
c      use hocoModule, only: epsc, barc
c      use hocoModule, only: barc
      double precision getHydynEpsc, getHydynBarc
#include "aaa.h"
#include "ho.h"
      character ctau*5, cbim*3, czz*5, txt*33
      common/jcen/jcentr,jmxcentr
      do itau=itau1,itau2
      nxoff=0  !3
      nyoff=0  !3
      delx5=(xmaxhy-xminhy)/(nxhy-1)*0.5
      dely5=(ymaxhy-yminhy)/(nyhy-1)*0.5
      deleta5=(zmaxhy-zminhy)/(nzhy-1)*0.5
      modu=nint(dtaupl/(getHydynTauhy(2)-getHydynTauhy(1)))
      ntau=1+modu*(itau-1)
      tau=getHydynTauhy(ntau)
      if(tau.le.tauzer+tauup)then
      write(ctau,'(f5.2)')tau
      write(cbim,'(a,i2)')'J',jcentr
      call RootCanvas(npad)
      call histo2BeginFigure()
      txt='energy density [GeV/fm^{3}]   ('
      call RootPave(txt//'(tau='//ctau//' fm/c)   '//cbim//';')
      call histo2BeginPlot(npad/2, 2)
      call RootPadMult(npad)
      do i=1,2
      do j=1,npad/2
      call RootPad((i-1)*npad/2+j)
      jj=2
      if(npad.eq.6)jj=3
      iz=nzhy/2+1 + (3-2*i) * (1+(j-1)*jj)
      zz=zminhy+(iz-1)*deleta5*2
      if(zz.ge.0.)write(czz,'(f4.1)')zz
      if(zz.lt.0.)write(czz,'(f5.1)')zz
      call RootHisto(2,'#eta_{s} ='//czz//';'
     *,nxhy-2*nxoff
     *,xminhy-delx5+nxoff*delx5*2
     *,xmaxhy+delx5-nxoff*delx5*2
     *,nyhy-2*nyoff
     *,yminhy-dely5+nyoff*dely5*2
     *,ymaxhy+dely5-nyoff*dely5*2)
      call histo2BeginSubPlot(
     * 'src/KW/freeze.f:ooHoEpsilonEtas(npad,itau1,itau2,dtaupl)'
     *,'ooHoEpsilonEtas')
      call histo2AddHeader(
     * txt//'($\tau$ ='//ctau//' fm/c) '//cbim//' $\eta_{s}$ ='//czz
     *,xminhy-delx5+nxoff*delx5*2
     *,xmaxhy+delx5-nxoff*delx5*2
     *,delx5*2
     *,yminhy-dely5+nyoff*dely5*2
     *,ymaxhy+dely5-nyoff*dely5*2
     *,dely5*2
     *,'x [fm]'
     *,'y [fm]')
      
      do ix=1,nxhy
         xx=xminhy+(ix-1)*delx5*2
         do iy=1,nyhy
            yy=yminhy+(iy-1)*dely5*2
            eps=getHydynEpsc(  iz,ntau,ix,iy)
            b1= getHydynBarc(1,iz,ntau,ix,iy)
            b2= getHydynBarc(2,iz,ntau,ix,iy)
            b3= getHydynBarc(3,iz,ntau,ix,iy)
            tem=PiEos(5,eps,b1,b2,b3)
            if (tem.gt.tfrout) then
               call RootFill(2,xx,yy,eps)
               call histo2FillArray(2,xx,yy,eps)
            endif
         enddo
      enddo
      call RootDraw('cont4Z;','x [fm];','y [fm];')  !surf1
      call histo2EndSubPlot()
      enddo
      enddo
      call histo2EndPlot()
      endif
      enddo
      end

c-----------------------------------------------------------------------
      subroutine ooHoTemperatureEtas(npad,itau1,itau2,dtaupl)
c-----------------------------------------------------------------------
c      use ArrayModule, only: getd
c      use hocoModule, only: epsc, barc
c      use hocoModule, only: barc
      double precision getHydynEpsc, getHydynBarc
#include "aaa.h"
#include "ho.h"
      character ctau*5, cbim*3, czz*5, txt*21
      common/jcen/jcentr,jmxcentr
      do itau=itau1,itau2
      nxoff=0  !3
      nyoff=0  !3
      delx5=(xmaxhy-xminhy)/(nxhy-1)*0.5
      dely5=(ymaxhy-yminhy)/(nyhy-1)*0.5
      deleta5=(zmaxhy-zminhy)/(nzhy-1)*0.5
      modu=nint(dtaupl/(getHydynTauhy(2)-getHydynTauhy(1)))
      ntau=1+modu*(itau-1)
      tau=getHydynTauhy(ntau)
      if(tau.le.tauzer+tauup)then
      write(ctau,'(f5.2)')tau
      write(cbim,'(a,i2)')'J',jcentr
      call RootCanvas(npad)
      call histo2BeginFigure()
      txt='temp-temp_fo [MeV]  ('
      call RootPave(txt//'(tau='//ctau//' fm/c)   '//cbim//';')
      call histo2BeginPlot(npad/2, 2)
      call RootPadMult(npad)
      do i=1,2
      do j=1,npad/2
      call RootPad((i-1)*npad/2+j)
      jj=2
      if(npad.eq.6)jj=3
      iz=nzhy/2+1 + (3-2*i) * (1+(j-1)*jj)
      zz=zminhy+(iz-1)*deleta5*2
      if(zz.ge.0.)write(czz,'(f4.1)')zz
      if(zz.lt.0.)write(czz,'(f5.1)')zz
      call RootHisto(2,'#eta_{s} ='//czz//';'
     *,nxhy-2*nxoff
     *,xminhy-delx5+nxoff*delx5*2
     *,xmaxhy+delx5-nxoff*delx5*2
     *,nyhy-2*nyoff
     *,yminhy-dely5+nyoff*dely5*2
     *,ymaxhy+dely5-nyoff*dely5*2)
      call histo2BeginSubPlot(
     * 'src/KW/freeze.f:ooHoTemperatureEtas(npad,itau1,itau2,dtaupl)'
     *,'ooHoTemperatureEtas')
      call histo2AddHeader(
     * txt//'($\tau$ ='//ctau//' fm/c) '//cbim//' $\eta_{s}$ ='//czz
     *,xminhy-delx5+nxoff*delx5*2
     *,xmaxhy+delx5-nxoff*delx5*2
     *,delx5*2
     *,yminhy-dely5+nyoff*dely5*2
     *,ymaxhy+dely5-nyoff*dely5*2
     *,dely5*2
     *,'x [fm]'
     *,'y [fm]')

      do ix=1,nxhy
         xx=xminhy+(ix-1)*delx5*2
         do iy=1,nyhy
            yy=yminhy+(iy-1)*dely5*2
            eps=getHydynEpsc(  iz,ntau,ix,iy)
            b1= getHydynBarc(1,iz,ntau,ix,iy)
            b2= getHydynBarc(2,iz,ntau,ix,iy)
            b3= getHydynBarc(3,iz,ntau,ix,iy)
            tem=PiEos(5,eps,b1,b2,b3)
            val=0.
            if(tem.gt.tfrout)val=(tem-tfrout)*1000
            if(val.gt.0.) then 
               call RootFill(2,xx,yy,val)
               call histo2FillArray(2,xx,yy,val)
            endif
         enddo
      enddo
      call RootDraw('cont4Z;','x [fm];','y [fm];')  !surf1
      call histo2EndSubPlot() !surf1
      enddo
      enddo
      call histo2EndPlot()
      endif
      enddo
      end

c-----------------------------------------------------------------------
      subroutine ooHoEpsilonTau(npad,itau1,itau2,dtaupl)
c-----------------------------------------------------------------------
c      use ArrayModule, only: getd
c      use hocoModule, only: epsc, barc
c      use hocoModule, only: barc
      double precision getHydynEpsc, getHydynBarc

#include "aaa.h"
#include "ho.h"
      character ctau*5, cbim*3, txt*33
      common/jcen/jcentr,jmxcentr
      itau=itau1-1
      do jtau=1+(itau1-1)/npad,1+(itau2-1)/npad
      nxoff=0  !3
      nyoff=0  !3
      delx5=(xmaxhy-xminhy)/(nxhy-1)*0.5
      dely5=(ymaxhy-yminhy)/(nyhy-1)*0.5
      deleta5=(zmaxhy-zminhy)/(nzhy-1)*0.5
      modu=nint(dtaupl/(getHydynTauhy(2)-getHydynTauhy(1)))
      write(cbim,'(a,i2)')'J',jcentr
      call RootCanvas(npad)
      call histo2BeginFigure()
      txt='energy density [GeV/fm^{3}]   ('
      call RootPave(txt//'#eta_{s}=0)   '//cbim//';')
      call histo2BeginPlot(npad/2, 2)
      call RootPadMult(npad)
      do i=1,2
      do j=1,npad/2
      call RootPad((i-1)*npad/2+j)
      itau=itau+1
      if(itau.gt.itau2)return
      ntau=1+modu*(itau-1)
      tau=getHydynTauhy(ntau)
      write(ctau,'(f5.2)')tau
      iz=nzhy/2+1
      call RootHisto(2,'tau = '//ctau//';'
     *,nxhy-2*nxoff
     *,xminhy-delx5+nxoff*delx5*2
     *,xmaxhy+delx5-nxoff*delx5*2
     *,nyhy-2*nyoff
     *,yminhy-dely5+nyoff*dely5*2
     *,ymaxhy+dely5-nyoff*dely5*2)
      call histo2BeginSubPlot(
     * 'src/KW/freeze.f:ooHoEpsilonTau(npad,itau1,itau2,dtaupl)' 
     *,'ooHoEpsilonTau')
      call histo2AddHeader(
     * txt//'$\eta_{s}$ =0) '//cbim//' $\tau$ = '//ctau
     *,xminhy-delx5+nxoff*delx5*2
     *,xmaxhy+delx5-nxoff*delx5*2
     *,delx5*2
     *,yminhy-dely5+nyoff*dely5*2
     *,ymaxhy+dely5-nyoff*dely5*2
     *,dely5*2
     *,'x [fm]'
     *,'y [fm]')

      do ix=1,nxhy
         xx=xminhy+(ix-1)*delx5*2
         do iy=1,nyhy
            yy=yminhy+(iy-1)*dely5*2
            eps=getHydynEpsc(  iz,ntau,ix,iy)
            e = getHydynEpsc(  iz,ntau,ix,iy)
            b1= getHydynBarc(1,iz,ntau,ix,iy)
            b2= getHydynBarc(2,iz,ntau,ix,iy)
            b3= getHydynBarc(3,iz,ntau,ix,iy)
            tem=PiEos(5,e,b1,b2,b3)
            if(tem.gt.tfrout) then
               call RootFill(2,xx,yy,eps)
               call histo2FillArray(2,xx,yy,eps)
            endif
         enddo
      enddo
      call RootDraw('cont4Z;','x [fm];','y [fm];')
      call histo2EndSubPlot()
      enddo
      enddo
      call histo2EndPlot()
      enddo
      end

c-----------------------------------------------------------------------
      subroutine ooHoTangVelTau(npad,itau1,itau2,dtaupl)
c-----------------------------------------------------------------------
c      use ArrayModule, only: getd
c      use hocoModule, only: epsc, barc, velc
c      use hocoModule, only: barc, velc
      double precision getHydynEpsc, getHydynBarc
      double precision getHydynVelc

#include "aaa.h"
#include "ho.h"            
      character ctau*5, cbim*3, txt*26
      common/jcen/jcentr,jmxcentr
      itau=itau1-1
      do jtau=1+(itau1-1)/npad,1+(itau2-1)/npad
      nxoff=0  !3
      nyoff=0  !3
      delx5=(xmaxhy-xminhy)/(nxhy-1)*0.5
      dely5=(ymaxhy-yminhy)/(nyhy-1)*0.5
      deleta5=(zmaxhy-zminhy)/(nzhy-1)*0.5
      modu=nint(dtaupl/(getHydynTauhy(2)-getHydynTauhy(1)))
      write(cbim,'(a,i2)')'J',jcentr
      call RootCanvas(npad)
      call histo2BeginFigure()
      txt='tang velocity [% of c]   ('
      call RootPave(txt//'#eta_{s}=0)   '//cbim//';')
      call histo2BeginPlot(npad/2, 2)
      call RootPadMult(npad)
      do i=1,2
      do j=1,npad/2
      call RootPad((i-1)*npad/2+j)
      itau=itau+1
      if(itau.gt.itau2)return
      ntau=1+modu*(itau-1)
      tau=getHydynTauhy(ntau)
      write(ctau,'(f5.2)')tau
      iz=nzhy/2+1
      call RootHisto(2,'tau = '//ctau//';'
     *,nxhy-2*nxoff
     *,xminhy-delx5+nxoff*delx5*2
     *,xmaxhy+delx5-nxoff*delx5*2
     *,nyhy-2*nyoff
     *,yminhy-dely5+nyoff*dely5*2
     *,ymaxhy+dely5-nyoff*dely5*2)
      call histo2BeginSubPlot(
     * 'src/KW/freeze.f:ooHoTangVelTau(npad,itau1,itau2,dtaupl)'
     *,'ooHoTangVelTau')
      call histo2AddHeader(
     * txt//'$\eta_{s}$ = 0) '//cbim//' $\tau$ = '//ctau
     *,xminhy+nxoff*delx5*2
     *,xminhy+(nxhy-nxoff-1)*delx5*2
     *,delx5*2
     *,yminhy+nyoff*dely5*2
     *,yminhy+(nyhy-nyoff-1)*dely5*2
     *,dely5*2
     *,'x [fm];'
     *,'y [fm];')

      do ix=1+nxoff,nxhy-nxoff
         xx=xminhy+(ix-1)*delx5*2
         do iy=1+nyoff,nyhy-nyoff
            yy=yminhy+(iy-1)*dely5*2
            vx=getHydynVelc(1,iz,ntau,ix,iy)
            vy=getHydynVelc(2,iz,ntau,ix,iy)
            r=sqrt(xx**2+yy**2)
            phi=0
            if(r.gt.0.)phi=sign(1.,yy)*acos(xx/r)
            vtg=vx*sin(phi)-vy*cos(phi)
            vrd=vx*cos(phi)+vy*sin(phi)
            e =getHydynEpsc(  iz,ntau,ix,iy)
            b1=getHydynBarc(1,iz,ntau,ix,iy)
            b2=getHydynBarc(2,iz,ntau,ix,iy)
            b3=getHydynBarc(3,iz,ntau,ix,iy)
            tem=PiEos(5,e,b1,b2,b3)
            if(tem.gt.tfrout) then
               call RootFill(2,xx,yy,vtg*100)
               call histo2FillArray(2,xx,yy,vtg*100)
            endif
         enddo
      enddo
      call RootDraw('cont4Z;','x [fm];','y [fm];')
      call histo2EndSubPlot()
      enddo
      enddo
      call histo2EndPlot()
      enddo
      end

c-----------------------------------------------------------------------
      subroutine ooHoRadVelTau(npad,itau1,itau2,dtaupl)
c-----------------------------------------------------------------------
c      use ArrayModule, only: getd
c      use hocoModule, only: epsc, barc, velc
c      use hocoModule, only: barc, velc
      double precision getHydynEpsc, getHydynBarc
      double precision getHydynVelc

#include "aaa.h"
#include "ho.h"
      character ctau*5, cbim*3, txt*26
      common/jcen/jcentr,jmxcentr
      itau=itau1-1
      do jtau=1+(itau1-1)/npad,1+(itau2-1)/npad
      nxoff=0  !3
      nyoff=0  !3
      delx5=(xmaxhy-xminhy)/(nxhy-1)*0.5
      dely5=(ymaxhy-yminhy)/(nyhy-1)*0.5
      deleta5=(zmaxhy-zminhy)/(nzhy-1)*0.5
      modu=nint(dtaupl/(getHydynTauhy(2)-getHydynTauhy(1)))
      write(cbim,'(a,i2)')'J',jcentr
      call RootCanvas(npad)
      call histo2BeginFigure()
      txt='rad velocity [% of c]   ('
      call RootPave(txt//'#eta_{s}=0)   '//cbim//';')
      call histo2BeginPlot(npad/2, 2)
      call RootPadMult(npad)
      do i=1,2
      do j=1,npad/2
      call RootPad((i-1)*npad/2+j)
      itau=itau+1
      if(itau.gt.itau2)return
      ntau=1+modu*(itau-1)
      tau=getHydynTauhy(ntau)
      write(ctau,'(f5.2)')tau
      iz=nzhy/2+1
      call RootHisto(2,'tau = '//ctau//';'
     *,nxhy-2*nxoff
     *,xminhy-delx5+nxoff*delx5*2
     *,xmaxhy+delx5-nxoff*delx5*2
     *,nyhy-2*nyoff
     *,yminhy-dely5+nyoff*dely5*2
     *,ymaxhy+dely5-nyoff*dely5*2)
      call histo2BeginSubPlot(
     * 'src/KW/freeze.f:ooHoRadVelTau(npad,itau1,itau2,dtaupl)'  
     *,'ooHoRadVelTau')
      call histo2AddHeader(
     * txt//'$\eta_{s}$=0) '//cbim//' $\tau$ = '//ctau
     *,xminhy-delx5+nxoff*delx5*2
     *,xmaxhy+delx5-nxoff*delx5*2
     *,delx5*2
     *,yminhy-dely5+nyoff*dely5*2
     *,ymaxhy+dely5-nyoff*dely5*2
     *,dely5*2
     *,'x [fm]'
     *,'y [fm]')

      do ix=1+nxoff,nxhy-nxoff
         xx=xminhy+(ix-1)*delx5*2
         do iy=1+nyoff,nyhy-nyoff
            yy=yminhy+(iy-1)*dely5*2
            vx=getHydynVelc(1,iz,ntau,ix,iy)
            vy=getHydynVelc(2,iz,ntau,ix,iy)
            r=sqrt(xx**2+yy**2)
            phi=0
            if(r.gt.0.)phi=sign(1.,yy)*acos(xx/r)
            vtg=vx*sin(phi)-vy*cos(phi)
            vrd=vx*cos(phi)+vy*sin(phi)
            e =getHydynEpsc(  iz,ntau,ix,iy)
            b1=getHydynBarc(1,iz,ntau,ix,iy)
            b2=getHydynBarc(2,iz,ntau,ix,iy)
            b3=getHydynBarc(3,iz,ntau,ix,iy)
            tem=PiEos(5,e,b1,b2,b3)
            if(tem.gt.tfrout.and.vrd.gt.0.) then
               call RootFill(2,xx,yy,vrd*100)
               call histo2FillArray(2,xx,yy,vrd*100)
            endif
         enddo
      enddo
      call RootDraw('cont4Z;','x [fm];','y [fm];')
      call histo2EndSubPlot()
      enddo
      enddo
      call histo2EndPlot()
      enddo
      end




cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc




c------------------------------------------------------------------------------
      subroutine xxHoFoVol(nsurf)
c------------------------------------------------------------------------------
#include "aaa.h"
#include "ho.h"
      parameter (mgeni=20)
      common/cgenar/genar(mgeni)
      common/cetaplot/etaplot(20),ietaplot(20),mxetaplot
      common/jcen/jcentr,jmxcentr
      common/czhits/zhits(0:ncenthxx),zfhits(0:ncenthxx)
      logical d,f,g
      call getlogi(d,f,g)
      if(ntauhy.le.0)then
        tauhymax=0
      else
        tauhymax=getHydynTauhy(ntauhy)
      endif
      if(jcentr.gt.0)taumaxi=rclass(jtable(7),1,jcentr)
      do jeta=1,mxetaplot
      neta=ietaplot(jeta)
      if(g)write(ifhi,'(a)')    '!####################################'
      if(g)write(ifhi,'(a,i3)') '!        hydro freeze out vol    '
      if(g)write(ifhi,'(a)')    '!####################################'
      if(g)write(ifhi,'(a)') '!newpage'
      do ii=1,4
        if(ii.eq.1)nphi=2
        if(ii.ge.2)nphi=1+(ii-1)*nphihy/4
        if(g)write(ifhi,'(a,4i1)')'openhisto name vol-'
     .   ,jcentr,jeta,ii
        if(g)write(ifhi,'(a)')'htyp lin'
        if(ii.eq.1)then !----------------------
         if(g)write(ifhi,'(a,f4.1)')'xrange 0. '
     .    ,taumaxi
         if(g)write(ifhi,'(a)')'txt  "xaxis [t] (fm/c)"'
         if(g)write(ifhi,'(a)')
     .    'txt "title[f] = 0^o! 90^o! 180^o! 270^o!  "'
         if(g)write(ifhi,'(a)') 'xmod lin ymod log yrange 1 auto '
         if(g)write(ifhi,'(a,f4.1,a)') 'text 0.1 0.7 "[c]='
     .    ,etahy(neta),'"'
         if(g)write(ifhi,'(a)')'txt "yaxis dV^*! / d[t] d[c] d[f]  "'
         if(g)write(ifhi,'(a,i2,a)')'text 0.83 0.55 "J',jcentr,'"'

        endif       !-------------------------------
        if(g)write(ifhi,*)'histoweight ',zfhits(jcentr)
        if(g)write(ifhi,'(a)')'array 2'
        do ngeni=2,mgeni
         !formely:ntau loop & tau=getHydynTauhy(ntau)-(getHydynTauhy(2)-getHydynTauhy(1))/2
         tau=max(tauzer1,tauzer2)
     .  +(ngeni-1)*(taumaxi-max(tauzer1,tauzer2))/(mgeni-1)
         ntau=1
         do while(getHydynTauhy(ntau).lt.tau.and.ntau.lt.ntauhxx)
         ntau=ntau+1
         enddo
         if(ntau.gt.1)then
          if(abs(tau-getHydynTauhy(ntau-1)).
     .           lt.abs(tau-getHydynTauhy(ntau)))ntau=ntau-1
         endif
         if(d.and.tauhymax.lt.tau)then
           y=0
         elseif(d)then
           y=vlmaa(nsurf,neta,ntau,nphi)
         endif
         if(f)call addHo(0+nsurf,jeta,ngeni,ii,y)
         if(d)y=max(1.,y)
         if(g)write(ifhi,'(2e11.3)')tau,y
        enddo
        if(g)write(ifhi,'(a)') 'endarray closehisto '
        if(ii.lt.4.and.g)write(ifhi,'(a)') 'plot 0-'
        if(ii.eq.4.and.g)write(ifhi,'(a)') 'plot 0'
      enddo
      enddo
      end
c------------------------------------------------------------------------------
      subroutine xxHoFoEpsilon(nsurf)
c------------------------------------------------------------------------------
#include "aaa.h"
#include "ho.h"
      parameter (mgeni=20)
      common/cetaplot/etaplot(20),ietaplot(20),mxetaplot
      common/jcen/jcentr,jmxcentr
      common/czhits/zhits(0:ncenthxx),zfhits(0:ncenthxx)
      logical d,f,g
      call getlogi(d,f,g)
      if(ntauhy.le.0)then
        tauhymax=0
      else
        tauhymax=getHydynTauhy(ntauhy)
      endif
      if(jcentr.gt.0)taumaxi=rclass(jtable(7),1,jcentr)
      do jeta=1,mxetaplot
      neta=ietaplot(jeta)
      if(g)write(ifhi,'(a)')    '!####################################'
      if(g)write(ifhi,'(a,i3)') '!        hydro freeze out epsilon    '
      if(g)write(ifhi,'(a)')    '!####################################'
      if(g)write(ifhi,'(a)') '!newpage'
      do ii=1,4
        if(ii.eq.1)nphi=1
        if(ii.ge.2)nphi=1+(ii-1)*nphihy/4
        if(g)write(ifhi,'(a,4i1)')'openhisto name e-'
     .    ,jcentr,jeta,ii
        if(g)write(ifhi,'(a)')       'htyp lin '
        if(ii.eq.1)then !----------------------
         if(g)write(ifhi,'(a,f4.1)')'xmod lin xrange 0. '
     .    ,taumaxi
         if(g)write(ifhi,'(a)')'txt  "xaxis [t] (fm/c)"'
         if(g)write(ifhi,'(a)')
     .    'txt "title[f] = 0^o! 90^o! 180^o! 270^o! "'
         if(g)write(ifhi,'(a)') 'ymod lin yrange 0 auto '
         if(g)write(ifhi,'(a,f4.1,a)') 'text 0.45 0.1 "[c]='
     .    ,etahy(neta),'"'
         if(g)write(ifhi,'(a)')'txt "yaxis FO epsilon (GeV/fm3) "'
         if(g)write(ifhi,'(a,i2,a)')'text 0.83 0.55 "J',jcentr,'"'

        endif       !-------------------------------
        if(g)write(ifhi,*)'histoweight ',zfhits(jcentr)
        if(g)write(ifhi,'(a)')'array 2'
        do ngeni=2,mgeni
         !formely:ntau loop & tau=getHydynTauhy(ntau)-(getHydynTauhy(2)-getHydynTauhy(1))/2
         tau=max(tauzer1,tauzer2)
     .  +(ngeni-1)*(taumaxi-max(tauzer1,tauzer2))/(mgeni-1)
         ntau=1
         do while(getHydynTauhy(ntau).lt.tau.and.ntau.lt.ntauhxx)
         ntau=ntau+1
         enddo
         if(ntau.gt.1)then
          if(abs(tau-getHydynTauhy(ntau-1))
     .           .lt.abs(tau-getHydynTauhy(ntau)))ntau=ntau-1
         endif
         if(d.and.tauhymax.lt.tau)then
            y=0
          elseif(d)then
            y=epsaa(nsurf,neta,ntau,nphi)
          endif
          if(f)call addHo(42+nsurf,jeta,ngeni,ii,y)
          if(g)write(ifhi,'(2e11.3)')getHydynTauhy(ntau),y
        enddo
        if(g)write(ifhi,'(a)') 'endarray'
        if(g)write(ifhi,'(a)') 'closehisto '
        if(ii.lt.4.and.g)write(ifhi,'(a)') 'plot 0-'
        if(ii.eq.4.and.g)write(ifhi,'(a)') 'plot 0'
      enddo
      enddo
      end


c------------------------------------------------------------------------------
      subroutine xxHoFoRadius(nsurf)
c------------------------------------------------------------------------------
#include "aaa.h"
#include "ho.h"
      parameter (mgeni=20)
      common/cetaplot/etaplot(20),ietaplot(20),mxetaplot
      common/jcen/jcentr,jmxcentr
      common/czhits/zhits(0:ncenthxx),zfhits(0:ncenthxx)
      logical d,f,g
      call getlogi(d,f,g)
      if(ntauhy.le.0)then
        tauhymax=0
      else
        tauhymax=getHydynTauhy(ntauhy)
      endif
      if(jcentr.gt.0)taumaxi=rclass(jtable(7),1,jcentr)
      do jeta=1,mxetaplot
      neta=ietaplot(jeta)
      if(g)write(ifhi,'(a)')    '!####################################'
      if(g)write(ifhi,'(a,i3)') '!        hydro freeze out radius    '
      if(g)write(ifhi,'(a)')    '!####################################'
      if(g)write(ifhi,'(a)') '!newpage'
      do ii=1,4
        if(ii.eq.1)nphi=1
        if(ii.ge.2)nphi=1+(ii-1)*nphihy/4
        if(g)write(ifhi,'(a,4i1)')'openhisto name r-'
     .    ,jcentr,jeta,ii
        if(g)write(ifhi,'(a)')       'htyp lin '
        if(ii.eq.1)then !----------------------
         if(g)write(ifhi,'(a,f4.1)')'xmod lin xrange 0. '
     .    ,taumaxi
         if(g)write(ifhi,'(a)')'txt  "xaxis [t] (fm/c)"'
         if(g)write(ifhi,'(a)')
     .    'txt "title[f] = 0^o! 90^o! 180^o! 270^o! "'
         if(g)write(ifhi,'(a)') 'ymod lin yrange 0 auto '
         if(g)write(ifhi,'(a,f4.1,a)') 'text 0.45 0.1 "[c]='
     .    ,etahy(neta),'"'
         if(g)write(ifhi,'(a)')'txt "yaxis FO radius (fm) "'
         if(g)write(ifhi,'(a,i2,a)')'text 0.83 0.55 "J',jcentr,'"'

        endif       !-------------------------------
        if(g)write(ifhi,*)'histoweight ',zfhits(jcentr)
        if(g)write(ifhi,'(a)')'array 2'
        do ngeni=2,mgeni
         !formely:ntau loop & tau=getHydynTauhy(ntau)-(getHydynTauhy(2)-getHydynTauhy(1))/2
         tau=max(tauzer1,tauzer2)
     .  +(ngeni-1)*(taumaxi-max(tauzer1,tauzer2))/(mgeni-1)
         ntau=1
         do while(getHydynTauhy(ntau).lt.tau.and.ntau.lt.ntauhxx)
         ntau=ntau+1
         enddo
         if(ntau.gt.1)then
          if(abs(tau-getHydynTauhy(ntau-1))
     .           .lt.abs(tau-getHydynTauhy(ntau)))ntau=ntau-1
         endif
         if(d.and.tauhymax.lt.tau)then
            y=0
          elseif(d)then
            y=radaa(nsurf,neta,ntau,nphi)
          endif
          if(f)call addHo(2+nsurf,jeta,ngeni,ii,y)
          if(g)write(ifhi,'(2e11.3)')getHydynTauhy(ntau),y
        enddo
        if(g)write(ifhi,'(a)') 'endarray'
        if(g)write(ifhi,'(a)') 'closehisto '
        if(ii.lt.4.and.g)write(ifhi,'(a)') 'plot 0-'
        if(ii.eq.4.and.g)write(ifhi,'(a)') 'plot 0'
      enddo
      enddo
      end

c------------------------------------------------------------------------------
      subroutine xxHoFoRadVelocity(nsurf)
c------------------------------------------------------------------------------
#include "aaa.h"
#include "ho.h"
      parameter (mgeni=20)
      common/cetaplot/etaplot(20),ietaplot(20),mxetaplot
      common/jcen/jcentr,jmxcentr
      common/czhits/zhits(0:ncenthxx),zfhits(0:ncenthxx)
      logical d,f,g
      call getlogi(d,f,g)
      if(ntauhy.le.0)then
        tauhymax=0
      else
        tauhymax=getHydynTauhy(ntauhy)
      endif
      if(jcentr.gt.0)taumaxi=rclass(jtable(7),1,jcentr)
      do jeta=1,mxetaplot
      neta=ietaplot(jeta)
      if(g)write(ifhi,'(a)')    '!####################################'
      if(g)write(ifhi,'(a,i3)') '!  hydro freeze out rad velocity     '
      if(g)write(ifhi,'(a)')    '!####################################'
      if(g)write(ifhi,'(a)') '!newpage'
      do ii=1,4
        if(ii.eq.1)nphi=1
        if(ii.ge.2)nphi=1+(ii-1)*nphihy/4
        phi=phihy(nphi)
        if(g)write(ifhi,'(a,4i1)')'openhisto name vr-'
     .    ,jcentr,jeta,ii
        if(g)write(ifhi,'(a)')       'htyp lin '
        if(ii.eq.1)then !----------------------
         if(g)write(ifhi,'(a,f4.1)')'xmod lin xrange 0. '
     .    ,taumaxi
         if(g)write(ifhi,'(a)')'txt  "xaxis [t] (fm/c)"'
         if(g)write(ifhi,'(a)')
     .     'txt "title[f] = 0^o! 90^o! 180^o! 270^o!  "'
         if(g)write(ifhi,'(a)') 'ymod lin yrange auto auto '
         if(g)write(ifhi,'(a,f4.1,a)') 'text 0.45 0.1 "[c]='
     .    ,etahy(neta),'"'
         if(g)write(ifhi,'(a)')'txt "yaxis rad FO velocity "'
         if(g)write(ifhi,'(a,i2,a)')'text 0.83 0.55 "J',jcentr,'"'

        endif       !-------------------------------
        if(g)write(ifhi,*)'histoweight ',zfhits(jcentr)
        if(g)write(ifhi,'(a)')'array 2'
        do ngeni=2,mgeni
         !formely:ntau loop & tau=getHydynTauhy(ntau)-(getHydynTauhy(2)-getHydynTauhy(1))/2
         tau=max(tauzer1,tauzer2)
     .  +(ngeni-1)*(taumaxi-max(tauzer1,tauzer2))/(mgeni-1)
         ntau=1
         do while(getHydynTauhy(ntau).lt.tau.and.ntau.lt.ntauhxx)
         ntau=ntau+1
         enddo
         if(ntau.gt.1)then
          if(abs(tau-getHydynTauhy(ntau-1))
     .           .lt.abs(tau-getHydynTauhy(ntau)))ntau=ntau-1
         endif
         if(d.and.tauhymax.lt.tau)then
            y=0.
          elseif(d)then
            vx=velaa(nsurf,1,neta,ntau,nphi)
            vy=velaa(nsurf,2,neta,ntau,nphi)
            vrd=vx*cos(phi)+vy*sin(phi)
            vtg=vx*sin(phi)-vy*cos(phi)
            y=vrd
          endif
          if(f)call addHo(4+nsurf,jeta,ngeni,ii,y)
          if(g)write(ifhi,'(2e11.3)')tau,y
          !if(f.and.y.ne.0..and.jeta.eq.1.and.ii.eq.1)
          !&write(ifmt,*)'CHECK xxplot ',tau,y
        enddo
        if(g)write(ifhi,'(a)') 'endarray'
        if(g)write(ifhi,'(a)') 'closehisto '
        if(ii.lt.4.and.g)write(ifhi,'(a)') 'plot 0-'
        if(ii.eq.4.and.g)write(ifhi,'(a)') 'plot 0'
      enddo
      enddo
      end

c------------------------------------------------------------------------------
      subroutine xxHoFoTangVelocity(nsurf)
c------------------------------------------------------------------------------
#include "aaa.h"
#include "ho.h"
      parameter (mgeni=20)
      common/cetaplot/etaplot(20),ietaplot(20),mxetaplot
      common/jcen/jcentr,jmxcentr
      common/czhits/zhits(0:ncenthxx),zfhits(0:ncenthxx)
      logical d,f,g
      call getlogi(d,f,g)
      if(ntauhy.le.0)then
        tauhymax=0
      else
        tauhymax=getHydynTauhy(ntauhy)
      endif
      if(jcentr.gt.0)taumaxi=rclass(jtable(7),1,jcentr)
      do jeta=1,mxetaplot
      neta=ietaplot(jeta)
      if(g)write(ifhi,'(a)')    '!####################################'
      if(g)write(ifhi,'(a,i3)') '!    hydro freeze out tang velocity  '
      if(g)write(ifhi,'(a)')    '!####################################'
      if(g)write(ifhi,'(a)') '!newpage'
      do ii=1,5
        if(ii.eq.1)nphi=1
        if(ii.ge.2)nphi=1+(ii-2)*nphihy/4+nphihy/8
        phi=phihy(nphi)
        if(g)write(ifhi,'(a,4i1)')'openhisto name vt-'
     .    ,jcentr,jeta,ii
        if(g)write(ifhi,'(a)')       'htyp lin '
        if(ii.eq.1)then !----------------------
         if(g)write(ifhi,'(a,f4.1)')'xmod lin xrange 0. '
     .    ,taumaxi
         if(g)write(ifhi,'(a)')'txt  "xaxis [t] (fm/c)"'
         if(g)write(ifhi,'(a)')
     .    'txt "title[f] = 0^o! 45^o! 135^o! 225^o! 315^o! "'
         if(g)write(ifhi,'(a)') 'ymod lin yrange auto auto '
         if(g)write(ifhi,'(a,f4.1,a)') 'text 0.65 0.05 "   [c]='
     .    ,etahy(neta),'"'
         if(g)write(ifhi,'(a)')'txt "yaxis tang FO velocity "'
         if(g)write(ifhi,'(a,i2,a)')'text 0.83 0.55 "J',jcentr,'"'

        endif       !-------------------------------
        if(g)write(ifhi,*)'histoweight ',zfhits(jcentr)
        if(g)write(ifhi,'(a)')'array 2'
        do ngeni=2,mgeni
         !formely:ntau loop & tau=getHydynTauhy(ntau)-(getHydynTauhy(2)-getHydynTauhy(1))/2
         tau=max(tauzer1,tauzer2)
     .  +(ngeni-1)*(taumaxi-max(tauzer1,tauzer2))/(mgeni-1)
         ntau=1
         do while(getHydynTauhy(ntau).lt.tau.and.ntau.lt.ntauhxx)
         ntau=ntau+1
         enddo
         if(ntau.gt.1)then
          if(abs(tau-getHydynTauhy(ntau-1))
     .           .lt.abs(tau-getHydynTauhy(ntau)))ntau=ntau-1
         endif
         if(d.and.tauhymax.lt.tau)then
            y=0
          elseif(d)then
            vx=velaa(nsurf,1,neta,ntau,nphi)
            vy=velaa(nsurf,2,neta,ntau,nphi)
            vrd=vx*cos(phi)+vy*sin(phi)
            vtg=vx*sin(phi)-vy*cos(phi)
            y=vtg
          endif
          if(f)call addHo(6+nsurf,jeta,ngeni,ii,y)
          if(g)write(ifhi,'(2e11.3)')tau,y
        enddo
        if(g)write(ifhi,'(a)') 'endarray'
        if(g)write(ifhi,'(a)') 'closehisto '
        if(ii.lt.5.and.g)write(ifhi,'(a)') 'plot 0-'
        if(ii.eq.5.and.g)write(ifhi,'(a)') 'plot 0'
      enddo
      enddo
      end

c------------------------------------------------------------------------------
      subroutine xxHoFoBarmu(nsurf)
c------------------------------------------------------------------------------
#include "aaa.h"
#include "ho.h"
      parameter (mgeni=20)
      common/cetaplot/etaplot(20),ietaplot(20),mxetaplot
      common/jcen/jcentr,jmxcentr
      common/czhits/zhits(0:ncenthxx),zfhits(0:ncenthxx)
      logical d,f,g
      call getlogi(d,f,g)
      if(ntauhy.le.0)then
        tauhymax=0
      else
        tauhymax=getHydynTauhy(ntauhy)
      endif
      if(jcentr.gt.0)taumaxi=rclass(jtable(7),1,jcentr)
      do jeta=1,mxetaplot
      neta=ietaplot(jeta)
      if(g)write(ifhi,'(a)')    '!####################################'
      if(g)write(ifhi,'(a,i3)') '!        hydro freeze out bmu    '
      if(g)write(ifhi,'(a)')    '!####################################'
      if(g)write(ifhi,'(a)') '!newpage'
      do ii=1,4
        if(ii.eq.1)nphi=1
        if(ii.ge.2)nphi=1+(ii-1)*nphihy/4
        if(g)write(ifhi,'(a,4i1)')'openhisto name bmu-'
     .    ,jcentr,jeta,ii
        if(g)write(ifhi,'(a)')       'htyp lin '
        if(ii.eq.1)then !----------------------
         if(g)write(ifhi,'(a,f4.1)')'xmod lin xrange 0. '
     .    ,taumaxi
         if(g)write(ifhi,'(a)')'txt  "xaxis [t] (fm/c)"'
         if(g)write(ifhi,'(a)') 'txt  "title  [f]=0^o!   [f]=90^o!   "'
         if(g)write(ifhi,'(a)') 'ymod lin yrange 0 auto '
         if(g)write(ifhi,'(a,f4.1,a)') 'text 0.1 0.1 "   [c]='
     .    ,etahy(neta),'"'
         if(g)write(ifhi,'(a)')'txt "yaxis FO [m]?B! (GeV) "'
         if(g)write(ifhi,'(a,i2,a)')'text 0.83 0.55 "J',jcentr,'"'

        endif       !-------------------------------
        if(g)write(ifhi,*)'histoweight ',zfhits(jcentr)
        if(g)write(ifhi,'(a)')'array 2'
        do ngeni=2,mgeni
          !formely:ntau loop & tau=getHydynTauhy(ntau)-(getHydynTauhy(2)-getHydynTauhy(1))/2
          tau=max(tauzer1,tauzer2)
     .  +(ngeni-1)*(taumaxi-max(tauzer1,tauzer2))/(mgeni-1)
          ntau=1
          do while(getHydynTauhy(ntau).lt.tau.and.ntau.lt.ntauhxx)
          ntau=ntau+1
          enddo
          if(ntau.gt.1)then
          if(abs(tau-getHydynTauhy(ntau-1))
     .            .lt.abs(tau-getHydynTauhy(ntau)))ntau=ntau-1
          endif
          if(d.and.tauhymax.lt.tau)then
            y=0
          elseif(d)then
            e=epsaa(nsurf,neta,ntau,nphi)
            b1=baraa(nsurf,1,neta,ntau,nphi)
            b2=baraa(nsurf,2,neta,ntau,nphi)
            b3=baraa(nsurf,3,neta,ntau,nphi)
            barmu=PiEos(6,e,b1,b2,b3)
            y=barmu
          endif
          if(f)call addHo(8+nsurf,jeta,ngeni,ii,y)
          if(g)write(ifhi,'(2e11.3)')tau,y
        enddo
        if(g)write(ifhi,'(a)') 'endarray'
        if(g)write(ifhi,'(a)') 'closehisto '
        if(ii.lt.4.and.g)write(ifhi,'(a)') 'plot 0-'
        if(ii.eq.4.and.g)write(ifhi,'(a)') 'plot 0'
      enddo
      enddo
      end

c------------------------------------------------------------------------------
      subroutine xxHoEpsilon(modu)
c------------------------------------------------------------------------------
#include "aaa.h"
#include "ho.h"
      parameter (mgeni=20)
      common/cetaplot/etaplot(20),ietaplot(20),mxetaplot
      common/jcen/jcentr,jmxcentr
      common/czhits/zhits(0:ncenthxx),zfhits(0:ncenthxx)
      logical d,f,g
      call getlogi(d,f,g)
      if(ntauhy.le.0)then
        tauhymax=0
      else
        tauhymax=getHydynTauhy(ntauhy)
      endif
      if(jcentr.gt.0)taumaxi=rclass(jtable(7),1,jcentr)
      if(jcentr.gt.0)radmaxi=rclass(jtable(7),2,jcentr)
      !write(ifmt,*)'+++++',ikoevt,jcentr,d,f,g
      dlt=(taumaxi-tauzer)/(mgeni-1)*modu
      rdmx=radmaxi
      do jeta=1,mxetaplot
      neta=ietaplot(jeta)
      if(g)write(ifhi,'(a)')       '!##################################'
      if(g)write(ifhi,'(a,i3)')    '!   hydro epsilon vs r     '
      if(g)write(ifhi,'(a)')       '!##################################'
      if(g)write(ifhi,'(a)') '!newpage'
      iphi=0
      do nphi=1,1+nphihy/2,nphihy/4
       phi=phihy(nphi)
       iphi=iphi+1
       itau=0
       do ngeni=1,mgeni
         !formely:ntau loop & tau=getHydynTauhy(ntau)
         tau=max(tauzer1,tauzer2)
     .  +(ngeni-1)*(taumaxi-max(tauzer1,tauzer2))/(mgeni-1)
         ntau=1
         do while(getHydynTauhy(ntau).lt.tau.and.ntau.lt.ntauhxx)
         ntau=ntau+1
         enddo
         if(ntau.gt.1)then
          if(abs(tau-getHydynTauhy(ntau-1))
     .           .lt.abs(tau-getHydynTauhy(ntau)))ntau=ntau-1
         endif
         if(mod(ngeni-1,modu).eq.0.and.(ngeni-1)/modu.lt.5)then
         itau=itau+1
         if(g)write(ifhi,'(a,4i1)')'openhisto name eps-'
     .       ,jcentr,jeta,iphi,itau
         if(itau.eq.1.and.iphi.eq.1.and.g)write(ifhi,'(a)')'htyp lxx '
         if(itau.ne.1.and.iphi.eq.1.and.g)write(ifhi,'(a)')'htyp lin '
         if(itau.eq.1.and.iphi.eq.2.and.g)write(ifhi,'(a)')'htyp lxp '
         if(itau.ne.1.and.iphi.eq.2.and.g)write(ifhi,'(a)')'htyp lip '
         if(itau.eq.1.and.iphi.eq.3.and.g)write(ifhi,'(a)')'htyp lxa '
         if(itau.ne.1.and.iphi.eq.3.and.g)write(ifhi,'(a)')'htyp lia '
         if(iphi.eq.1.and.itau.eq.1)then !----------------------
          if(g)write(ifhi,'(a,f5.2)')'xmod lin xrange 0. ' ,rdmx
          if(g)write(ifhi,'(a)')'txt  "xaxis r (fm)"'
          if(g)write(ifhi,'(a)')
     .     'txt  "title[f]: __ 0^o!  ... 90^o!  - - 180^o! "'
          if(g)write(ifhi,'(a,f4.2,a)')'text 0.50 0.75 "[Dt]='
     .     ,dlt,'fm/c"'
          if(g)write(ifhi,'(a)') 'ymod log yrange 0.01 auto '
          if(g)write(ifhi,'(a,f4.1,a)') 'text 0.70 0.3 "[c]='
     .     ,etahy(neta),'"'
          if(g)write(ifhi,'(a)')'txt "yaxis [e] (GeV/fm^3!)"'
         if(g)write(ifhi,'(a,i2,a)')'text 0.83 0.55 "J',jcentr,'"'
         endif        !------------------------------
         if(g)write(ifhi,*)'histoweight ',zfhits(jcentr)
         if(g)write(ifhi,'(a)')'array 2'
         do kgeni=1,mgeni
           rad=kgeni*radmaxi/mgeni
           nrad=1
           do while(radhy(nrad).lt.rad.and.nrad.lt.nradhy)
           nrad=nrad+1
           enddo
           if(nrad.gt.1)then
           if(abs(rad-radhy(nrad-1)).lt.abs(rad-radhy(nrad)))nrad=nrad-1
           endif
           if(d.and.(radhy(nradhy).lt.rad.or.tauhymax.lt.tau))then
            y=0
           elseif(d)then
            call epsioget(1,neta,ntau,nphi,nrad, y )
           endif
           if(f)call addHo(10+itau,jeta,kgeni,iphi,y)
           if(g)write(ifhi,'(2e11.3)')rad,y
         enddo
         if(g)write(ifhi,'(a)') 'endarray closehisto plot 0-'
        endif
       enddo
      enddo
      if(g)write(ifhi,'(a)') 'openhisto htyp fkp xmod lin ymod log'
      if(g)write(ifhi,'(a,2e11.3,a)')'xrange',0.,rdmx
     . ,' yrange 0.01 auto '
      if(g)write(ifhi,'(a)') 'function  0.01 closehisto  plot 0'
      enddo
      end

c------------------------------------------------------------------------------
      subroutine xxHoEntropy(modu)
c------------------------------------------------------------------------------
#include "aaa.h"
#include "ho.h"
      parameter (mgeni=20)
      common/cetaplot/etaplot(20),ietaplot(20),mxetaplot
      common/jcen/jcentr,jmxcentr
      common/czhits/zhits(0:ncenthxx),zfhits(0:ncenthxx)
      logical d,f,g
      call getlogi(d,f,g)
      if(ntauhy.le.0)then
        tauhymax=0
      else
        tauhymax=getHydynTauhy(ntauhy)
      endif
      if(jcentr.gt.0)taumaxi=rclass(jtable(7),1,jcentr)
      if(jcentr.gt.0)radmaxi=rclass(jtable(7),2,jcentr)
      dlt=(taumaxi-tauzer)/(mgeni-1)*modu
      rdmx=radmaxi
      do jeta=1,mxetaplot
      neta=ietaplot(jeta)
      if(g)write(ifhi,'(a)')       '!##################################'
      if(g)write(ifhi,'(a,i3)')    '!   hydro sigma vs r     '
      if(g)write(ifhi,'(a)')       '!##################################'
      if(g)write(ifhi,'(a)') '!newpage'
      iphi=0
      do nphi=1,1+nphihy/4,nphihy/4
       phi=phihy(nphi)
       iphi=iphi+1
       itau=0
       do ngeni=1,mgeni
         !formely:ntau loop & tau=getHydynTauhy(ntau)
         tau=max(tauzer1,tauzer2)
     .  +(ngeni-1)*(taumaxi-max(tauzer1,tauzer2))/(mgeni-1)
         ntau=1
         do while(getHydynTauhy(ntau).lt.tau.and.ntau.lt.ntauhxx)
         ntau=ntau+1
         enddo
         if(ntau.gt.1)then
          if(abs(tau-getHydynTauhy(ntau-1))
     .           .lt.abs(tau-getHydynTauhy(ntau)))ntau=ntau-1
         endif
         if(mod(ngeni-1,modu).eq.0.and.(ngeni-1)/modu.lt.5)then
         itau=itau+1
         if(g)write(ifhi,'(a,5i1)')'openhisto name sig-'
     .       ,jcentr,jeta,iphi,itau
         if(itau.eq.1.and.iphi.eq.1.and.g)write(ifhi,'(a)')'htyp lxx '
         if(itau.ne.1.and.iphi.eq.1.and.g)write(ifhi,'(a)')'htyp lin '
         if(itau.eq.1.and.iphi.eq.2.and.g)write(ifhi,'(a)')'htyp lxp '
         if(itau.ne.1.and.iphi.eq.2.and.g)write(ifhi,'(a)')'htyp lip '
         if(iphi.eq.1.and.itau.eq.1)then !----------------------
          if(g)write(ifhi,'(a,f5.2)')'xmod lin xrange 0. ' ,rdmx
          if(g)write(ifhi,'(a)')'txt  "xaxis r (fm)"'
          if(g)write(ifhi,'(a)')'txt  "title__ [f]=0^o! ... [f]=90^o! "'
          if(g)write(ifhi,'(a,f4.2,a)')'text 0.50 0.75 "[Dt]='
     .     ,dlt,'fm/c"'
          if(g)write(ifhi,'(a)') 'ymod log yrange 0.1 auto '
          if(g)write(ifhi,'(a,f4.1,a)') 'text 0.65 0.9 "   [c]='
     .     ,etahy(neta),'"'
          if(g)write(ifhi,'(a)')'txt "yaxis [s] (1/fm^3!)"'
         if(g)write(ifhi,'(a,i2,a)')'text 0.83 0.55 "J',jcentr,'"'

         endif        !------------------------------
         if(g)write(ifhi,*)'histoweight ',zfhits(jcentr)
         if(g)write(ifhi,'(a)')'array 2'
         do kgeni=1,mgeni
           rad=kgeni*radmaxi/mgeni
           nrad=1
           do while(radhy(nrad).lt.rad.and.nrad.lt.nradhy)
           nrad=nrad+1
           enddo
           if(nrad.gt.1)then
           if(abs(rad-radhy(nrad-1)).lt.abs(rad-radhy(nrad)))nrad=nrad-1
           endif
           if(d.and.(radhy(nradhy).lt.rad.or.tauhymax.lt.tau))then
            y=0
           elseif(d)then
            call epsioget(2,neta,ntau,nphi,nrad, y )
           endif
           if(f)call addHo(15+itau,jeta,kgeni,iphi,y)
           if(g)write(ifhi,'(2e11.3)')rad,y
         enddo
         if(g)write(ifhi,'(a)') 'endarray closehisto plot 0-'
        endif
       enddo
      enddo
      if(g)write(ifhi,'(a)') 'openhisto htyp fkp xmod lin ymod log'
      if(g)write(ifhi,'(a,2e11.3,a)')'xrange',0.,rdmx
     . ,' yrange 0.01 auto '
      if(g)write(ifhi,'(a)') 'function  0.01 closehisto  plot 0'
      enddo
      end

c------------------------------------------------------------------------------
      subroutine xxHoTemperature(modu)
c------------------------------------------------------------------------------
#include "aaa.h"
#include "ho.h"
      parameter (mgeni=20)
#include "so.h"
      common/cetaplot/etaplot(20),ietaplot(20),mxetaplot
      common/jcen/jcentr,jmxcentr
      common/czhits/zhits(0:ncenthxx),zfhits(0:ncenthxx)
      logical d,f,g
      call getlogi(d,f,g)
      if(ntauhy.le.0)then
        tauhymax=0
      else
        tauhymax=getHydynTauhy(ntauhy)
      endif
      if(jcentr.gt.0)taumaxi=rclass(jtable(7),1,jcentr)
      if(jcentr.gt.0)radmaxi=rclass(jtable(7),2,jcentr)
      dlt=(taumaxi-tauzer)/(mgeni-1)*modu
      rdmx=radmaxi
      do jeta=1,mxetaplot
      neta=ietaplot(jeta)
      if(g)write(ifhi,'(a)')       '!##################################'
      if(g)write(ifhi,'(a,i3)')    '!   hydro temp vs r     '
      if(g)write(ifhi,'(a)')       '!##################################'
      if(g)write(ifhi,'(a)') '!newpage'
      iphi=0
      do nphi=1,1+nphihy/4,nphihy/4
       phi=phihy(nphi)
       iphi=iphi+1
       itau=0
       do ngeni=1,mgeni
         !formely:ntau loop & tau=getHydynTauhy(ntau)
         tau=max(tauzer1,tauzer2)
     .  +(ngeni-1)*(taumaxi-max(tauzer1,tauzer2))/(mgeni-1)
         ntau=1
         do while(getHydynTauhy(ntau).lt.tau.and.ntau.lt.ntauhxx)
         ntau=ntau+1
         enddo
         if(ntau.gt.1)then
          if(abs(tau-getHydynTauhy(ntau-1))
     .           .lt.abs(tau-getHydynTauhy(ntau)))ntau=ntau-1
         endif
         if(mod(ngeni-1,modu).eq.0.and.(ngeni-1)/modu.lt.5)then
         itau=itau+1
         if(g)write(ifhi,'(a,5i1)')'openhisto name temp-'
     .       ,jcentr,jeta,iphi,itau
         if(itau.eq.1.and.iphi.eq.1.and.g)write(ifhi,'(a)')'htyp lxx '
         if(itau.ne.1.and.iphi.eq.1.and.g)write(ifhi,'(a)')'htyp lin '
         if(itau.eq.1.and.iphi.eq.2.and.g)write(ifhi,'(a)')'htyp lxp '
         if(itau.ne.1.and.iphi.eq.2.and.g)write(ifhi,'(a)')'htyp lip '
         if(iphi.eq.1.and.itau.eq.1)then !----------------------
          if(g)write(ifhi,'(a,f5.2)')'xmod lin xrange 0. ' ,rdmx
          if(g)write(ifhi,'(a)')'txt  "xaxis r (fm)"'
          if(g)write(ifhi,'(a)')'txt  "title__ [f]=0^o! . . [f]=90^o! "'
          if(g)write(ifhi,'(a,f4.2,a)')'text 0.50 0.75 "[Dt]='
     .     ,dlt,'fm/c"'
          if(g)write(ifhi,'(a,f6.3,a)') 'ymod lin yrange ',tempfinal
     .     ,' 0.30'
          if(g)write(ifhi,'(a,f4.1,a)') 'text 0.65 0.9 "   [c]='
     .     ,etahy(neta),'"'
          if(g)write(ifhi,'(a)')'txt "yaxis T (GeV)"'
         if(g)write(ifhi,'(a,i2,a)')'text 0.83 0.55 "J',jcentr,'"'

         endif        !------------------------------
         if(g)write(ifhi,*)'histoweight ',zfhits(jcentr)
         if(g)write(ifhi,'(a)')'array 2'
         do kgeni=1,mgeni
           rad=kgeni*radmaxi/mgeni
           nrad=1
           do while(radhy(nrad).lt.rad.and.nrad.lt.nradhy)
           nrad=nrad+1
           enddo
           if(nrad.gt.1)then
           if(abs(rad-radhy(nrad-1)).lt.abs(rad-radhy(nrad)))nrad=nrad-1
           endif
           if(d.and.(radhy(nradhy).lt.rad.or.tauhymax.lt.tau))then
            y=0
           elseif(d)then
            call epsioget(1,neta,ntau,nphi,nrad, e )   !e.lt.oEeos???
            call barioget(1,neta,ntau,nphi,nrad, b1 )
            call barioget(2,neta,ntau,nphi,nrad, b2 )
            call barioget(3,neta,ntau,nphi,nrad, b3 )
            temp=PiEos(5,e,b1,b2,b3)
            y=temp
           endif
           if(f)call addHo(20+itau,jeta,kgeni,iphi,y)
           if(g)write(ifhi,'(2e11.3)')rad,y
         enddo
         if(g)write(ifhi,'(a)') 'endarray closehisto plot 0-'
        endif
       enddo
      enddo
      if(g)write(ifhi,'(a)') 'openhisto htyp fkp xmod lin'
      if(g)write(ifhi,'(a,2e11.3,a)')'xrange',0.,rdmx
      if(g)write(ifhi,'(a,f6.3,a)') 'ymod lin yrange ',tempfinal,' 0.30'
      if(g)write(ifhi,'(a)') 'function  0.01 closehisto  plot 0'
      enddo
      end

c------------------------------------------------------------------------------
      subroutine xxHoRadVelocity(modu)
c------------------------------------------------------------------------------
#include "aaa.h"
#include "ho.h"
      parameter (mgeni=20)
      common/cetaplot/etaplot(20),ietaplot(20),mxetaplot
      common/jcen/jcentr,jmxcentr
      common/czhits/zhits(0:ncenthxx),zfhits(0:ncenthxx)
      logical d,f,g
      call getlogi(d,f,g)
      if(ntauhy.le.0)then
        tauhymax=0
      else
        tauhymax=getHydynTauhy(ntauhy)
      endif
      if(jcentr.gt.0)taumaxi=rclass(jtable(7),1,jcentr)
      if(jcentr.gt.0)radmaxi=rclass(jtable(7),2,jcentr)
      dlt=(taumaxi-tauzer)/(mgeni-1)*modu
      rdmx=radmaxi
      do jeta=1,mxetaplot
      neta=ietaplot(jeta)
      if(g)write(ifhi,'(a)')       '!##################################'
      if(g)write(ifhi,'(a,i3)')    '!   hydro velocity vs r       '
      if(g)write(ifhi,'(a)')       '!##################################'
      if(g)write(ifhi,'(a)') '!newpage'
      iphi=0
      do nphi=1,1+nphihy/2,nphihy/4
       phi=phihy(nphi)
       iphi=iphi+1
       itau=0
       do ngeni=1,mgeni
         !formely:ntau loop & tau=getHydynTauhy(ntau)
         tau=max(tauzer1,tauzer2)
     .  +(ngeni-1)*(taumaxi-max(tauzer1,tauzer2))/(mgeni-1)
         ntau=1
         do while(getHydynTauhy(ntau).lt.tau.and.ntau.lt.ntauhxx)
         ntau=ntau+1
         enddo
         if(ntau.gt.1)then
          if(abs(tau-getHydynTauhy(ntau-1))
     .           .lt.abs(tau-getHydynTauhy(ntau)))ntau=ntau-1
         endif
         if(mod(ngeni-1,modu).eq.0.and.(ngeni-1)/modu.lt.5)then
         itau=itau+1
         if(g)write(ifhi,'(a,5i1)')'openhisto name radvel-'
     .       ,jcentr,jeta,iphi,itau
         if(itau.eq.1.and.iphi.eq.1.and.g)write(ifhi,'(a)')'htyp lxx '
         if(itau.ne.1.and.iphi.eq.1.and.g)write(ifhi,'(a)')'htyp lin '
         if(itau.eq.1.and.iphi.eq.2.and.g)write(ifhi,'(a)')'htyp lxp '
         if(itau.ne.1.and.iphi.eq.2.and.g)write(ifhi,'(a)')'htyp lip '
         if(itau.eq.1.and.iphi.eq.3.and.g)write(ifhi,'(a)')'htyp lxa '
         if(itau.ne.1.and.iphi.eq.3.and.g)write(ifhi,'(a)')'htyp lia '
         if(iphi.eq.1.and.itau.eq.1)then !----------------------
          if(g)write(ifhi,'(a,f5.2)')'xmod lin xrange 0. ' ,rdmx
          if(g)write(ifhi,'(a)')'txt  "xaxis r (fm)"'
          if(g)write(ifhi,'(a)')
     .     'txt  "title[f]: __  0^o!  ... 90^o!  - - 180^o!"'
          if(g)write(ifhi,'(a,f4.2,a)')'text 0.05 0.7 "[Dt]='
     .     ,dlt,'fm/c"'
          if(g)write(ifhi,'(a)') 'ymod lin yrange auto auto '
          if(g)write(ifhi,'(a,f4.1,a)') 'text 0.70 0.3 "[c]='
     .     ,etahy(neta),'"'
          if(g)write(ifhi,'(a)')
     .          'txt "yaxis rad velocity"'
         if(g)write(ifhi,'(a,i2,a)')'text 0.83 0.55 "J',jcentr,'"'

         endif       !-------------------------------
         if(g)write(ifhi,*)'histoweight ',zfhits(jcentr)
         if(g)write(ifhi,'(a)')'array 2'
         do kgeni=1,mgeni
           rad=kgeni*radmaxi/mgeni
           nrad=1
           do while(radhy(nrad).lt.rad.and.nrad.lt.nradhy)
           nrad=nrad+1
           enddo
           if(nrad.gt.1)then
           if(abs(rad-radhy(nrad-1)).lt.abs(rad-radhy(nrad)))nrad=nrad-1
           endif
           if(d.and.(radhy(nradhy).lt.rad.or.tauhymax.lt.tau))then
            y=0
           elseif(d)then
            call velioget(1,neta,ntau,nphi,nrad, y1 ) 
            call velioget(2,neta,ntau,nphi,nrad, y2 )
            y=y1*cos(phi)+y2*sin(phi)
           endif
           if(f)call addHo(25+itau,jeta,kgeni,iphi,y)
           if(g)write(ifhi,'(2e11.3)')rad,y
         enddo
         if(g)write(ifhi,'(a)') 'endarray'
         if(g)write(ifhi,'(a)') ' closehisto '
         if(g)write(ifhi,'(a)') 'plot 0-'
        endif
       enddo
      enddo
      if(g)write(ifhi,'(a)')'openhisto htyp fkp xmod lin ymod lin'
      if(g)write(ifhi,'(a)')'xrange 0 1 yrange 0 1 '
      if(g)write(ifhi,'(a)')'function 0 from 20 to 21 closehisto plot 0'
      enddo
      end

c------------------------------------------------------------------------------
      subroutine xxHoLongVelocity(modu)
c------------------------------------------------------------------------------
#include "aaa.h"
#include "ho.h"
      parameter (mgeni=20)
      common/cetaplot/etaplot(20),ietaplot(20),mxetaplot
      common/jcen/jcentr,jmxcentr
      common/czhits/zhits(0:ncenthxx),zfhits(0:ncenthxx)
      logical d,f,g
      call getlogi(d,f,g)
      if(ntauhy.le.0)then
        tauhymax=0
      else
        tauhymax=getHydynTauhy(ntauhy)
      endif
      if(jcentr.gt.0)taumaxi=rclass(jtable(7),1,jcentr)
      if(jcentr.gt.0)radmaxi=rclass(jtable(7),2,jcentr)
      dlt=(taumaxi-tauzer)/(mgeni-1)*modu
      rdmx=radmaxi
      do jeta=1,mxetaplot
      neta=ietaplot(jeta)
      if(g)write(ifhi,'(a)')       '!##################################'
      if(g)write(ifhi,'(a,i3)')    '!   hydro velocity vs r       '
      if(g)write(ifhi,'(a)')       '!##################################'
      if(g)write(ifhi,'(a)') '!newpage'
      iphi=0
      do nphi=1,1+nphihy/2,nphihy/4
       phi=phihy(nphi)
       iphi=iphi+1
       itau=0
       do ngeni=1,mgeni
         !formely:ntau loop & tau=getHydynTauhy(ntau)
         tau=max(tauzer1,tauzer2)
     .  +(ngeni-1)*(taumaxi-max(tauzer1,tauzer2))/(mgeni-1)
         ntau=1
         do while(getHydynTauhy(ntau).lt.tau.and.ntau.lt.ntauhxx)
         ntau=ntau+1
         enddo
         if(ntau.gt.1)then
          if(abs(tau-getHydynTauhy(ntau-1))
     .           .lt.abs(tau-getHydynTauhy(ntau)))ntau=ntau-1
         endif
         if(mod(ngeni-1,modu).eq.0.and.(ngeni-1)/modu.lt.5)then
         itau=itau+1
         if(g)write(ifhi,'(a,5i1)')'openhisto name longvel-'
     .       ,jcentr,jeta,iphi,itau
         if(itau.eq.1.and.iphi.eq.1.and.g)write(ifhi,'(a)')'htyp lxx '
         if(itau.ne.1.and.iphi.eq.1.and.g)write(ifhi,'(a)')'htyp lin '
         if(itau.eq.1.and.iphi.eq.2.and.g)write(ifhi,'(a)')'htyp lxp '
         if(itau.ne.1.and.iphi.eq.2.and.g)write(ifhi,'(a)')'htyp lip '
         if(itau.eq.1.and.iphi.eq.3.and.g)write(ifhi,'(a)')'htyp lxa '
         if(itau.ne.1.and.iphi.eq.3.and.g)write(ifhi,'(a)')'htyp lia '
         if(iphi.eq.1.and.itau.eq.1)then !----------------------
          if(g)write(ifhi,'(a,f5.2)')'xmod lin xrange 0. ' ,rdmx
          if(g)write(ifhi,'(a)')'txt  "xaxis r (fm)"'
          if(g)write(ifhi,'(a)')
     .     'txt  "title[f]: __ 0^o!  ... 90^o!  - - 180^o!"'
          if(g)write(ifhi,'(a,f3.1,a)')'text 0.05 0.7 "[Dt]='
     .     ,dlt,'fm/c"'
          if(g)write(ifhi,'(a)') 'ymod lin yrange auto auto '
          if(g)write(ifhi,'(a,f4.1,a)') 'text 0.70 0.75 "[c]='
     .     ,etahy(neta),'"'
          if(g)write(ifhi,'(a)')
     .          'txt "yaxis v?z!([t],[c],[f],r)/c   (comov)"'
         if(g)write(ifhi,'(a,i2,a)')'text 0.83 0.55 "J',jcentr,'"'

         endif       !-------------------------------
         if(g)write(ifhi,*)'histoweight ',zfhits(jcentr)
         if(g)write(ifhi,'(a)')'array 2'
         do kgeni=1,mgeni
           rad=kgeni*radmaxi/mgeni
           nrad=1
           do while(radhy(nrad).lt.rad.and.nrad.lt.nradhy)
           nrad=nrad+1
           enddo
           if(nrad.gt.1)then
           if(abs(rad-radhy(nrad-1)).lt.abs(rad-radhy(nrad)))nrad=nrad-1
           endif
           if(d.and.(radhy(nradhy).lt.rad.or.tauhymax.lt.tau))then
            y=0
           elseif(d)then
            call velioget(3,neta,ntau,nphi,nrad, y )
           endif
           if(f)call addHo(30+itau,jeta,kgeni,iphi,y)
           if(g)write(ifhi,'(2e11.3)')rad,y
         enddo
         if(g)write(ifhi,'(a)') 'endarray'
         if(g)write(ifhi,'(a)') ' closehisto '
         if(g)write(ifhi,'(a)') 'plot 0-'
        endif
       enddo
      enddo
      if(g)write(ifhi,'(a)')'openhisto htyp fkp xmod lin ymod lin'
      if(g)write(ifhi,'(a)')'xrange 0 1 yrange 0 1 '
      if(g)write(ifhi,'(a)')'function 0 from 20 to 21 closehisto plot 0'
      enddo
      end

c------------------------------------------------------------------------------
      subroutine xxHoBaryon(modu)
c------------------------------------------------------------------------------
#include "aaa.h"
#include "ho.h"
      parameter (mgeni=20)
      common/cetaplot/etaplot(20),ietaplot(20),mxetaplot
      common/jcen/jcentr,jmxcentr
      common/czhits/zhits(0:ncenthxx),zfhits(0:ncenthxx)
      logical d,f,g
      call getlogi(d,f,g)
      if(ntauhy.le.0)then
        tauhymax=0
      else
        tauhymax=getHydynTauhy(ntauhy)
      endif
      if(jcentr.gt.0)taumaxi=rclass(jtable(7),1,jcentr)
      if(jcentr.gt.0)radmaxi=rclass(jtable(7),2,jcentr)
      dlt=(taumaxi-tauzer)/(mgeni-1)*modu
      rdmx=radmaxi
      do jeta=1,mxetaplot
      neta=ietaplot(jeta)
      if(g)write(ifhi,'(a)')       '!##################################'
      if(g)write(ifhi,'(a,i3)')    '!   hydro baryon vs r     '
      if(g)write(ifhi,'(a)')       '!##################################'
      if(g)write(ifhi,'(a)') '!newpage'
      iphi=0
      do nphi=1,1+nphihy/4,nphihy/4
       phi=phihy(nphi)
       iphi=iphi+1
       itau=0
       do ngeni=1,mgeni
         !formely:ntau loop & tau=getHydynTauhy(ntau)
         tau=max(tauzer1,tauzer2)
     .  +(ngeni-1)*(taumaxi-max(tauzer1,tauzer2))/(mgeni-1)
         ntau=1
         do while(getHydynTauhy(ntau).lt.tau.and.ntau.lt.ntauhxx)
         ntau=ntau+1
         enddo
         if(ntau.gt.1)then
          if(abs(tau-getHydynTauhy(ntau-1))
     .           .lt.abs(tau-getHydynTauhy(ntau)))ntau=ntau-1
         endif
         if(mod(ngeni-1,modu).eq.0.and.(ngeni-1)/modu.lt.5)then
         itau=itau+1
         if(g)write(ifhi,'(a,5i1)')'openhisto name bar-'
     .       ,jcentr,jeta,iphi,itau
         if(itau.eq.1.and.iphi.eq.1.and.g)write(ifhi,'(a)')'htyp lxx '
         if(itau.ne.1.and.iphi.eq.1.and.g)write(ifhi,'(a)')'htyp lin '
         if(itau.eq.1.and.iphi.eq.2.and.g)write(ifhi,'(a)')'htyp lxp '
         if(itau.ne.1.and.iphi.eq.2.and.g)write(ifhi,'(a)')'htyp lip '
         if(iphi.eq.1.and.itau.eq.1)then !----------------------
          if(g)write(ifhi,'(a,f5.2)')'xmod lin xrange 0. ' ,rdmx
          if(g)write(ifhi,'(a)')'txt  "xaxis r (fm)"'
          if(g)write(ifhi,'(a)')'txt  "title__ [f]=0^o! . . [f]=90^o! "'
          if(g)write(ifhi,'(a,f4.2,a)')'text 0.50 0.75 "[Dt]='
     .     ,dlt,'fm/c"'
          if(g)write(ifhi,'(a)') 'ymod log yrange 0.0001 auto '
          if(g)write(ifhi,'(a,f4.1,a)') 'text 0.65 0.9 "   [c]='
     .     ,etahy(neta),'"'
          if(g)write(ifhi,'(a)')'txt "yaxis n?B! (1/fm^3!)"'
         if(g)write(ifhi,'(a,i2,a)')'text 0.83 0.55 "J',jcentr,'"'

         endif        !------------------------------
         if(g)write(ifhi,*)'histoweight ',zfhits(jcentr)
         if(g)write(ifhi,'(a)')'array 2'
         do kgeni=1,mgeni
           rad=kgeni*radmaxi/mgeni
           nrad=1
           do while(radhy(nrad).lt.rad.and.nrad.lt.nradhy)
           nrad=nrad+1
           enddo
           if(nrad.gt.1)then
           if(abs(rad-radhy(nrad-1)).lt.abs(rad-radhy(nrad)))nrad=nrad-1
           endif
           if(d.and.(radhy(nradhy).lt.rad.or.tauhymax.lt.tau))then
             y=0
           elseif(d)then
             call barioget(1,neta,ntau,nphi,nrad, y )
           endif
           if(f)call addHo(35+itau,jeta,kgeni,iphi,y)
           if(g)write(ifhi,'(2e11.3)')rad,y
         enddo
         if(g)write(ifhi,'(a)') 'endarray closehisto plot 0-'
        endif
       enddo
      enddo
      if(g)write(ifhi,'(a)') 'openhisto name hy-dummy1 htyp fkv'
      if(g)write(ifhi,'(a)') 'xmod lin ymod log'
      if(g)write(ifhi,'(a,e11.3,a)')'xrange 0 ',rdmx
     . ,' yrange 0.0001 auto '
      if(g)write(ifhi,'(a)') 'function 0.0001 closehisto  plot 0'
      enddo
      end

c------------------------------------------------------------------------------
      subroutine xxHoEpsilonEta(modu)
c------------------------------------------------------------------------------
#include "aaa.h"
#include "ho.h"
      parameter (mgeni=20)
      common/cetaplot/etaplot(20),ietaplot(20),mxetaplot
      common/jcen/jcentr,jmxcentr
      common/czhits/zhits(0:ncenthxx),zfhits(0:ncenthxx)
      logical d,f,g
      call getlogi(d,f,g)
      if(ntauhy.le.0)then
        tauhymax=0
      else
        tauhymax=getHydynTauhy(ntauhy)
      endif
      if(jcentr.gt.0)taumaxi=rclass(jtable(7),1,jcentr)
      if(jcentr.gt.0)radmaxi=rclass(jtable(7),2,jcentr)
      dlt=(taumaxi-tauzer)/(mgeni-1)*modu
      etamx=etahy(nzhy)+1
      if(g)write(ifhi,'(a)')       '!##################################'
      if(g)write(ifhi,'(a,i3)')    '!   hydro epsilon vs eta      '
      if(g)write(ifhi,'(a)')       '!##################################'
      if(g)write(ifhi,'(a)') '!newpage'
      itau=0
      do ngeni=1,mgeni
        !formely:ntau loop & tau=getHydynTauhy(ntau)
        tau=max(tauzer1,tauzer2)
     .  +(ngeni-1)*(taumaxi-max(tauzer1,tauzer2))/(mgeni-1)
        ntau=1
        do while(getHydynTauhy(ntau).lt.tau.and.ntau.lt.ntauhxx)
        ntau=ntau+1
        enddo
        if(ntau.gt.1)then
         if(abs(tau-getHydynTauhy(ntau-1))
     .          .lt.abs(tau-getHydynTauhy(ntau)))ntau=ntau-1
        endif
        if(mod(ngeni-1,modu).eq.0.and.(ngeni-1)/modu.lt.5)then
         itau=itau+1
         if(g)write(ifhi,'(a,3i1)')'openhisto htyp lin name epseta-'
     .    ,jcentr,itau
         if(itau.eq.1)then !----------------------
          if(g)write(ifhi,'(a,f5.2)')'xmod lin xrange 0. ' ,etamx
          if(g)write(ifhi,'(a)')'txt  "xaxis [c]"'
          if(g)write(ifhi,'(a)')'txt  "title "'
          if(g)write(ifhi,'(a,f3.1,a)')'text 0.55 0.85 "[Dt]='
     .     ,dlt,'fm/c"'
          if(g)write(ifhi,'(a)') 'yrange 0.01 50 '
          if(g)write(ifhi,'(a)')'txt "yaxis [e] (GeV/fm^3!)"'
          if(g)write(ifhi,'(a,i2,a)')'text 0.83 0.55 "J',jcentr,'"'
          if(g.and.iHyEpsilonEta.eq.1)write(ifhi,'(a)')'ymod log '
          if(g.and.iHyEpsilonEta.eq.2)write(ifhi,'(a)')'ymod lin '
         endif        !------------------------------
         if(g)write(ifhi,*)'histoweight ',zfhits(jcentr)
         if(g)write(ifhi,'(a)')'array 2'
         do neta=1,nzhy
           if(d.and.tauhymax.lt.tau)then
             y=0
           elseif(d)then
             y=epsij(neta,ntau)
           endif
           if(f)call addHo(41,1,neta,itau,y)
           if(g)write(ifhi,'(2e11.3)')etahy(neta),y
         enddo
         if(g)write(ifhi,'(a)') 'endarray closehisto plot 0-'
        endif
      enddo
      if(g)write(ifhi,'(a)') 'openhisto htyp fkp xmod lin ymod log'
      if(g)write(ifhi,'(a,2e11.3,a)')'xrange',0.,etamx
     . ,' yrange 0.01 auto '
      if(g)write(ifhi,'(a)') 'function  0.01 closehisto  plot 0'
      end

c------------------------------------------------------------------------------
      subroutine xxHoBaryonEta(modu)
c------------------------------------------------------------------------------
#include "aaa.h"
#include "ho.h"
      parameter (mgeni=20)
      common/cetaplot/etaplot(20),ietaplot(20),mxetaplot
      common/jcen/jcentr,jmxcentr
      common/czhits/zhits(0:ncenthxx),zfhits(0:ncenthxx)
      character*10 fafa
      logical d,f,g
      call getlogi(d,f,g)
      if(ntauhy.le.0)then
        tauhymax=0
      else
        tauhymax=getHydynTauhy(ntauhy)
      endif
      if(jcentr.gt.0)taumaxi=rclass(jtable(7),1,jcentr)
      if(jcentr.gt.0)radmaxi=rclass(jtable(7),2,jcentr)
      fafa='(a,f5.2,a)'
      dlt=(taumaxi-tauzer)/(mgeni-1)*modu
      etamx=etahy(nzhy)+1
      if(g)write(ifhi,'(a)')       '!##################################'
      if(g)write(ifhi,'(a,i3)')    '!   hydro baryon vs eta      '
      if(g)write(ifhi,'(a)')       '!##################################'
      if(g)write(ifhi,'(a)') '!newpage'
      itau=0
      do ngeni=1,mgeni
        !formely:ntau loop & tau=getHydynTauhy(ntau)
        tau=max(tauzer1,tauzer2)
     .  +(ngeni-1)*(taumaxi-max(tauzer1,tauzer2))/(mgeni-1)
        ntau=1
        do while(getHydynTauhy(ntau).lt.tau.and.ntau.lt.ntauhxx)
        ntau=ntau+1
        enddo
        if(ntau.gt.1)then
         if(abs(tau-getHydynTauhy(ntau-1))
     .          .lt.abs(tau-getHydynTauhy(ntau)))ntau=ntau-1
        endif
         if(mod(ngeni-1,modu).eq.0.and.(ngeni-1)/modu.lt.5)then
         itau=itau+1
         if(g)write(ifhi,'(a,3i1)')'openhisto htyp lin name bareta-'
     .    ,jcentr,itau
         if(itau.eq.1)then !----------------------
          if(g)write(ifhi, fafa)'xmod lin xrange 0. ' ,etamx
          if(g)write(ifhi,'(a)')'txt  "xaxis [c]"'
          if(g)write(ifhi,'(a)')'txt  "title "'
          if(g)write(ifhi, fafa)'text 0.55 0.85 "[Dt]=',dlt,'fm/c"'
          if(g)write(ifhi,'(a)')'yrange 0.0001 2 '
          if(g)write(ifhi,'(a)')'txt "yaxis n?B! (1/fm^3!)"'
          if(g)write(ifhi,'(a,i2,a)')'text 0.83 0.55 "J',jcentr,'"'
          if(g.and.iHyBaryonEta.eq.1)write(ifhi,'(a)')'ymod log '
          if(g.and.iHyBaryonEta.eq.2)write(ifhi,'(a)')'ymod lin '
         endif        !------------------------------
         if(g)write(ifhi,*)'histoweight ',zfhits(jcentr)
         if(g)write(ifhi,'(a)')'array 2'
         do neta=1,nzhy
           if(d.and.tauhymax.lt.tau)then
             y=0
           elseif(d)then
             y=barij(1,neta,ntau)
           endif
           if(f)call addHo(41,2,neta,itau,y)
           if(g)write(ifhi,'(2e11.3)')etahy(neta),y
         enddo
         if(g)write(ifhi,'(a)') 'endarray closehisto  plot 0-'
        endif
      enddo
      if(g)write(ifhi,'(a)') 'openhisto htyp fkp'
      if(g)write(ifhi,'(a)') 'xmod lin ymod log'
      if(g)write(ifhi,'(a,e11.3,a)')'xrange 0 ',etamx
     . ,' yrange 0.0001 auto'
      if(g)write(ifhi,'(a)') 'function 0.0001 closehisto  plot 0'
      end

c------------------------------------------------------------------------------
      subroutine xxHoBarmuEta(modu)
c------------------------------------------------------------------------------
#include "aaa.h"
#include "ho.h"
      parameter (mgeni=20)
      common/cetaplot/etaplot(20),ietaplot(20),mxetaplot
      common/jcen/jcentr,jmxcentr
      common/czhits/zhits(0:ncenthxx),zfhits(0:ncenthxx)
      logical d,f,g
      call getlogi(d,f,g)
      if(ntauhy.le.0)then
        tauhymax=0
      else
        tauhymax=getHydynTauhy(ntauhy)
      endif
      if(jcentr.gt.0)taumaxi=rclass(jtable(7),1,jcentr)
      if(jcentr.gt.0)radmaxi=rclass(jtable(7),2,jcentr)
      dlt=(taumaxi-tauzer)/(mgeni-1)*modu
      etamx=etahy(nzhy)+1
      if(g)write(ifhi,'(a)')       '!##################################'
      if(g)write(ifhi,'(a,i3)')    '!   hydro baryon vs eta      '
      if(g)write(ifhi,'(a)')       '!##################################'
      if(g)write(ifhi,'(a)') '!newpage'
      itau=0
      do ngeni=1,mgeni
        !formely:ntau loop & tau=getHydynTauhy(ntau)
        tau=max(tauzer1,tauzer2)
     .  +(ngeni-1)*(taumaxi-max(tauzer1,tauzer2))/(mgeni-1)
        ntau=1
        do while(getHydynTauhy(ntau).lt.tau.and.ntau.lt.ntauhxx)
        ntau=ntau+1
        enddo
        if(ntau.gt.1)then
         if(abs(tau-getHydynTauhy(ntau-1))
     .          .lt.abs(tau-getHydynTauhy(ntau)))ntau=ntau-1
        endif
        if(mod(ngeni-1,modu).eq.0.and.(ngeni-1)/modu.lt.5)then
         itau=itau+1
         if(g)write(ifhi,'(a,3i1)')'openhisto htyp lin name bareta-'
     .    ,jcentr,itau
         if(itau.eq.1)then !----------------------
          if(g)write(ifhi,'(a,f5.2)')'xmod lin xrange 0. ' ,etamx
          if(g)write(ifhi,'(a)')'txt  "xaxis [c]"'
          if(g)write(ifhi,'(a)')'txt  "title "'
          if(g)write(ifhi,'(a,f4.2,a)')'text 0.55 0.85 "[Dt]='
     .     ,dlt,'fm/c"'
          if(g)write(ifhi,'(a)') 'ymod log yrange 0.0001 2 '
          if(g)write(ifhi,'(a)')'txt "yaxis [m]?B! (1/fm^3!)"'
         if(g)write(ifhi,'(a,i2,a)')'text 0.83 0.55 "J',jcentr,'"'

         endif        !------------------------------
         if(g)write(ifhi,*)'histoweight ',zfhits(jcentr)
         if(g)write(ifhi,'(a)')'array 2'
         do neta=1,nzhy
           if(d.and.tauhymax.lt.tau)then
             y=0
           elseif(d)then
             e=epsij(neta,ntau)
             b1=barij(1,neta,ntau)
             b2=barij(2,neta,ntau)
             b3=barij(3,neta,ntau)
             barmu=PiEos(6,e,b1,b2,b3)
             y=barmu
           endif
           if(f)call addHo(41,3,neta,itau,y)
           if(g)write(ifhi,'(2e11.3)')etahy(neta),y
         enddo
         if(g)write(ifhi,'(a)') 'endarray closehisto  plot 0-'
        endif
      enddo
      if(g)write(ifhi,'(a)') 'openhisto htyp fkp'
      if(g)write(ifhi,'(a)') 'xmod lin ymod log'
      if(g)write(ifhi,'(a,e11.3,a)')'xrange 0 ',etamx
     . ,' yrange 0.0001 auto'
      if(g)write(ifhi,'(a)') 'function 0.0001 closehisto  plot 0'
      end

c------------------------------------------------------------------------------
      subroutine xxHoEntropyEta(modu)
c------------------------------------------------------------------------------
#include "aaa.h"
#include "ho.h"
      parameter (mgeni=20)
      common/cetaplot/etaplot(20),ietaplot(20),mxetaplot
      common/jcen/jcentr,jmxcentr
      common/czhits/zhits(0:ncenthxx),zfhits(0:ncenthxx)
      logical d,f,g
      call getlogi(d,f,g)
      if(ntauhy.le.0)then
        tauhymax=0
      else
        tauhymax=getHydynTauhy(ntauhy)
      endif
      if(jcentr.gt.0)taumaxi=rclass(jtable(7),1,jcentr)
      if(jcentr.gt.0)radmaxi=rclass(jtable(7),2,jcentr)
      dlt=(taumaxi-tauzer)/(mgeni-1)*modu
      etamx=etahy(nzhy)+1
      if(g)write(ifhi,'(a)')       '!##################################'
      if(g)write(ifhi,'(a,i3)')    '!   hydro entropy vs eta      '
      if(g)write(ifhi,'(a)')       '!##################################'
      if(g)write(ifhi,'(a)') '!newpage'
      itau=0
      do ngeni=1,mgeni
        !formely:ntau loop & tau=getHydynTauhy(ntau)
        tau=max(tauzer1,tauzer2)
     .  +(ngeni-1)*(taumaxi-max(tauzer1,tauzer2))/(mgeni-1)
        ntau=1
        do while(getHydynTauhy(ntau).lt.tau.and.ntau.lt.ntauhxx)
        ntau=ntau+1
        enddo
        if(ntau.gt.1)then
         if(abs(tau-getHydynTauhy(ntau-1))
     .          .lt.abs(tau-getHydynTauhy(ntau)))ntau=ntau-1
        endif 
         if(mod(ngeni-1,modu).eq.0.and.(ngeni-1)/modu.lt.5)then
         itau=itau+1
         if(g)write(ifhi,'(a,3i1)')'openhisto htyp lin name enteta-'
     .    ,jcentr,itau
         if(itau.eq.1)then !----------------------
          if(g)write(ifhi,'(a,f5.2)')'xmod lin xrange 0. ' ,etamx
          if(g)write(ifhi,'(a)')'txt  "xaxis [c]"'
          if(g)write(ifhi,'(a)')'txt  "title "'
          if(g)write(ifhi,'(a,f4.2,a)')'text 0.55 0.85 "[Dt]='
     .     ,dlt,'fm/c"'
          if(g)write(ifhi,'(a)') 'ymod log yrange 0.1 auto '
          if(g)write(ifhi,'(a)')'txt "yaxis [s] (1/fm^3!)"'
         if(g)write(ifhi,'(a,i2,a)')'text 0.83 0.55 "J',jcentr,'"'

         endif        !------------------------------
         if(g)write(ifhi,*)'histoweight ',zfhits(jcentr)
         if(g)write(ifhi,'(a)')'array 2'
         do neta=1,nzhy
           if(d.and.tauhymax.lt.tau)then
             y=0
           elseif(d)then
             y=sigij(neta,ntau)
           endif
           if(f)call addHo(41,4,neta,itau,y)
           if(g)write(ifhi,'(2e11.3)')etahy(neta),y
         enddo
         if(g)write(ifhi,'(a)') 'endarray closehisto plot 0-'
        endif
      enddo
      if(g)write(ifhi,'(a)') 'openhisto htyp fkp'
      if(g)write(ifhi,'(a)') 'xmod lin ymod log'
      if(g)write(ifhi,'(a,e11.3,a)')'xrange 0 ',etamx
     . ,' yrange 0.1 auto '
      if(g)write(ifhi,'(a)') 'function 0.1'
      if(g)write(ifhi,'(a)') 'closehisto  plot 0'
      end

c------------------------------------------------------------------------------
      subroutine xxHoEtc
c------------------------------------------------------------------------------
#include "aaa.h"
#include "ho.h"
      parameter (mgeni=20)
      common/czhits/zhits(0:ncenthxx),zfhits(0:ncenthxx)
      logical d,f,g
      call getlogi(d,f,g)
      if(g)write(ifhi,'(a)')   '!##################################'
      if(g)write(ifhi,'(a,i3)')'!   hydro etc      '
      if(g)write(ifhi,'(a)')   '!##################################'
      if(g)write(ifhi,'(a)') '!newpage'
      if(g)write(ifhi,'(a)')'openhisto htyp lin name etc'
      if(g)write(ifhi,'(a)')'xmod lin xrange 0.5 5.5'
      if(g)write(ifhi,'(a)')'txt  "xaxis centrality"'
      if(g)write(ifhi,'(a)')'txt  "title "'
      if(g)write(ifhi,'(a)')'ymod log yrange 0.001 1.2 '
      if(g)write(ifhi,'(a)')'txt "yaxis fraction of fluid events"'
      if(g)write(ifhi,*)'histoweight ',0d0
      if(g)write(ifhi,'(a)')'array 3'
      do i=1,jmxcentr
       y=0.
       if(zhits(i).ne.0.)y=zfhits(i)/zhits(i)
       if(g)write(ifhi,'(3e11.3)')1.0*i,y, zhits(i)
      enddo
      if(g)write(ifhi,'(a)') 'endarray closehisto plot 0'
      end

c------------------------------------------------------------------------------
      subroutine xxHoEtcTest
c------------------------------------------------------------------------------
#include "aaa.h"
#include "ho.h"
      common/jcen/jcentr,jmxcentr
      common/czhits/zhits(0:ncenthxx),zfhits(0:ncenthxx)
      logical d,f,g
      call getlogi(d,f,g)
      if(g)write(ifhi,'(a)')   '!##################################'
      if(g)write(ifhi,'(a,i3)')'!   hydro etc Test     '
      if(g)write(ifhi,'(a)')   '!##################################'
      if(g)write(ifhi,'(a)') '!newpage'
      if(g)write(ifhi,'(a)')'openhisto htyp pnt name etc'
      if(g)write(ifhi,'(a)')'xmod lin xrange 0.5 5.5'
      if(g)write(ifhi,'(a)')'txt  "xaxis centrality"'
      if(g)write(ifhi,'(a)')'txt  "title "'
      if(g)write(ifhi,'(a)')'ymod log yrange 0.001 1.2 '
      if(g)write(ifhi,'(a)')'txt "yaxis fraction of fluid events"'
      i=jcentr
      if(g)write(ifhi,*)'histoweight ',zhits(i)
      if(g)write(ifhi,'(a)')'array 3'
      y=0.
      if(zhits(i).ne.0.)y=zfhits(i)/zhits(i)
      if(g)write(ifhi,'(2e11.3)')1.0*i,y, zhits(i)
      if(g)write(ifhi,'(a)') 'endarray closehisto plot 0'
      end

c------------------------------------------------------------------------------
      subroutine xxHoAverages
c------------------------------------------------------------------------------
#include "aaa.h"
#include "ho.h"
      parameter (mgeni=20)
      common/cetaplot/etaplot(20),ietaplot(20),mxetaplot
      common/jcen/jcentr,jmxcentr
      common/czhits/zhits(0:ncenthxx),zfhits(0:ncenthxx)
      logical d,f,g
      call getlogi(d,f,g)
      if(ntauhy.le.0)then
        tauhymax=0
      else
        tauhymax=getHydynTauhy(ntauhy)
      endif
      if(jcentr.gt.0)taumaxi=rclass(jtable(7),1,jcentr)
      if(jcentr.gt.0)radmaxi=rclass(jtable(7),2,jcentr)
      do jeta=1,mxetaplot
      neta=ietaplot(jeta)
      if(g)write(ifhi,'(a)')    '!####################################'
      if(g)write(ifhi,'(a,i3)') '!        averages    '
      if(g)write(ifhi,'(a)')    '!####################################'
      if(g)write(ifhi,'(a)') '!newpage'
      if(g)write(ifhi,'(a,4i1)')'openhisto name av-'
     .,jcentr,jeta
      if(g)write(ifhi,'(a)')       'htyp lin '
      !----------------------
      if(g)write(ifhi,'(a,f4.1)')'xmod lin xrange 0. '
     . ,getHydynTauhy(ntauhxx)
      if(g)write(ifhi,'(a)')'txt  "xaxis [t] (fm/c)"'
      if(g)write(ifhi,'(a)') 'txt  "title   "'
      if(g)write(ifhi,'(a)') 'ymod lin yrange auto auto '
      if(g)write(ifhi,'(a,f4.1,a)') 'text 0.3 0.9 "   [c]='
     . ,etahy(neta),'"'
      if(g)write(ifhi,'(a)')'txt ""yaxis [e]?x!  [e]?p! "L#v?t!"G#  ""'
      if(g)write(ifhi,'(a,i2,a)')'text 0.83 0.55 "J',jcentr,'"'

      !-------------------------------
      if(g)write(ifhi,*)'histoweight ',zfhits(jcentr)
        if(g)write(ifhi,'(a)')'array 2'
      do ngeni=1,mgeni
        !formely:ntau loop & tau=getHydynTauhy(ntau)
        tau=max(tauzer1,tauzer2)
     .  +(ngeni-1)*(taumaxi-max(tauzer1,tauzer2))/(mgeni-1)
        ntau=1
        do while(getHydynTauhy(ntau).lt.tau.and.ntau.lt.ntauhxx)
        ntau=ntau+1
        enddo
        if(ntau.gt.1)then
         if(abs(tau-getHydynTauhy(ntau-1))
     .          .lt.abs(tau-getHydynTauhy(ntau)))ntau=ntau-1
        endif 
        if(d.and.tauhymax.lt.tau)then
          y=0
        elseif(d)then
          y=eccxav(neta,ntau)
        endif
        if(f)call addHo(42,jeta,ntau,1,y)
        if(g)write(ifhi,'(2e11.3)')getHydynTauhy(ntau),y
      enddo
      if(g)write(ifhi,'(a)') 'endarray closehisto plot 0-'
      if(g)write(ifhi,'(a,2i1)')'openhisto htyp lin'
      if(g)write(ifhi,*)'histoweight ',zfhits(jcentr)
        if(g)write(ifhi,'(a)')'array 2'
      do ngeni=1,mgeni
        !formely:ntau loop & tau=getHydynTauhy(ntau)
        tau=max(tauzer1,tauzer2)
     .  +(ngeni-1)*(taumaxi-max(tauzer1,tauzer2))/(mgeni-1)
        ntau=1
        do while(getHydynTauhy(ntau).lt.tau.and.ntau.lt.ntauhxx)
        ntau=ntau+1
        enddo
        if(ntau.gt.1)then
         if(abs(tau-getHydynTauhy(ntau-1))
     .          .lt.abs(tau-getHydynTauhy(ntau)))ntau=ntau-1
        endif
        if(d.and.tauhymax.lt.tau)then
          y=0
        elseif(d)then
          y=eccpav(neta,ntau)
        endif
        if(f)call addHo(42,jeta,ntau,2,y)
        if(g)write(ifhi,'(2e11.3)')getHydynTauhy(ntau),y
      enddo
      if(g)write(ifhi,'(a)') 'endarray closehisto plot 0-'
      if(g)write(ifhi,'(a,4i1)')'openhisto name avt-'
     .,jcentr,jeta
      if(g)write(ifhi,'(a)')       'htyp lin '
      if(g)write(ifhi,*)'histoweight ',zfhits(jcentr)
        if(g)write(ifhi,'(a)')'array 2'
      do ngeni=1,mgeni
        !formely:ntau loop & tau=getHydynTauhy(ntau)
        tau=max(tauzer1,tauzer2)
     .  +(ngeni-1)*(taumaxi-max(tauzer1,tauzer2))/(mgeni-1)

        ntau=1
        do while(getHydynTauhy(ntau).lt.tau.and.ntau.lt.ntauhxx)
        ntau=ntau+1
        enddo
        if(ntau.gt.1)then
         if(abs(tau-getHydynTauhy(ntau-1))
     .          .lt.abs(tau-getHydynTauhy(ntau)))ntau=ntau-1
        endif
        if(d.and.tauhymax.lt.tau)then
          y=0
        elseif(d)then
          y=vtraav(neta,ntau)
        endif
        if(f)call addHo(42,jeta,ntau,3,y)
        if(g)write(ifhi,'(2e11.3)')getHydynTauhy(ntau),y
      enddo
      if(g)write(ifhi,'(a)') 'endarray closehisto  plot 0'
      enddo
      end

c------------------------------------------------------------------------------
      subroutine xxxSource
c------------------------------------------------------------------------------
#include "aaa.h"
#include "ho.h"
       common/kcen/kcentr,kmxcentr
       if(iSource.gt.0)then
          call getKcentr
          call setHo(1,1,0)
          call xxSource
          if(nrevt.eq.nevent)then
            write(ifmt,'(a)')'plot Source'
            call setHo(0,2,1)
            do kcentr=1,kmxcentr
             call xxSource
            enddo
          endif
        endif
        end

c------------------------------------------------------------------------------
      subroutine xxSource
c------------------------------------------------------------------------------
#include "aaa.h"
#include "ho.h"
      common/kcen/kcentr,kmxcentr
      common/cti/xspamax,kspamax
      common/czkhits/zkhits(0:ncenthxx)
      character*3 b
      double precision p(5),x(4),v(5),k_t,q(4)
      parameter (kspamaxx=40)
      logical d,f,g
      integer isum(4)
      real fkt(0:4)
      data fkt /0.15, 0.25, 0.35, 0.45, 0.60/
      data isum/4*0/
      save isum
      xspamax=20
      if(maproj.eq.1.and.matarg.eq.1)xspamax=5
      kspamax=kspamaxx
      call getlogi(d,f,g)
      b='(a)'
      k2=iSource
      nKtBins=nrclass(k2)   !number of Kt bins (<=10)
      if(nKtBins.ne.4)stop'xxSource error\n\n'
      fkt(0)    =   rclass(k2,1,1)
      do i=1,nrclass(k2)
      if(abs(fkt(i-1)-rclass(k2,1,i)).gt.0.0001)
     . stop'xxSource error 2\n\n'
      fkt(i)    =   rclass(k2,2,i)
      print*,'+++++++',i,fkt(i-1),fkt(i)
      enddo


      do kkk=2,6,2
      ! 1 ... prim ptl distr
      ! 2 ... all ptl distr
      ! 3 ... prim pair distr
      ! 4 ... all pair distr
      ! 5 ... prim pair distr in LCMS
      ! 6 ... all pair distr in LCMS
      if(d)then
        do i=nptlpt+1,nptl
          if(istptl(i).le.1)then
          id=idptl(i)
          if(id.eq. 120.or.id.eq. 211
     .   .or.id.eq. 130.or.id.eq. 321
     .   .or.id.eq.1120.or.id.eq.2212)then
          rap=9999
          pt=sqrt(pptl(1,i)**2+pptl(2,i)**2)
          pti=pt
          amt=sqrt(pptl(5,i)**2+pt**2)
          pp=pptl(4,i)+abs(pptl(3,i))
          if(pp.gt.0.d0)rap=sign(1.,pptl(3,i))*log(pp/amt)
          if(abs(rap).lt.0.5.and.pt.gt.0.15.and.pt.lt.0.8)then
            ior=iorptl(i)
            if(kkk.eq.1.and.ior.le.0.and.ior.ne.-999.or.kkk.eq.2)then
              do k=1,4
                if(k.eq.1.and.(id.eq. 120.or.id.eq. 211)
     .         .or.k.eq.2.and.(id.eq. 120.or.id.eq. 211)
     .         .or.k.eq.3.and.(id.eq. 130.or.id.eq. 321)
     .         .or.k.eq.4.and.(id.eq.1120.or.id.eq.2212))then
                do m=1,4
                  xi=xorptl(m,i)
                  yi=1
                  no=0
                  if(k.eq.2)yi=pptl(m,i)
                  if(k.eq.2)no=(kkk-1)*4+m
                  call addSou((kkk-1)*16+(k-1)*4+m, 0, no, xi, yi)
                enddo
                endif
              enddo
            elseif(kkk.eq.3.and.ior.le.0.and.ior.ne.-999
     .       .or.kkk.eq.4.  .or.kkk.eq.6
     .       .or.kkk.eq.5.and.ior.le.0.and.ior.ne.-999)then
              if(id.eq.120.or.id.eq.211)then  !pions
              do j=i+1,nptl
                if(istptl(j).le.1)then
                id=idptl(j)
                if(id.eq.120.or.id.eq.211)then  !pions
                rap=9999
                pt=sqrt(pptl(1,j)**2+pptl(2,j)**2)
                ptj=pt
                amt=sqrt(pptl(5,j)**2+pt**2)
                pp=pptl(4,j)+abs(pptl(3,j))
                if(pp.gt.0.d0)rap=sign(1.,pptl(3,j))*log(pp/amt)
                if(abs(rap).lt.0.5.and.pt.gt.0.15.and.pt.lt.0.8)then
                  ior=iorptl(j)
                  if(kkk.eq.3.and.ior.le.0.and.ior.ne.-999
     .              .or.kkk.eq.4 .or.kkk.eq.6
     .              .or.kkk.eq.5.and.ior.le.0.and.ior.ne.-999)then
                    do l=1,4
                      x(l)=xorptl(l,i)-xorptl(l,j)
                    enddo
                    do l=1,4
                      p(l)=pptl(l,i)+pptl(l,j)
                      q(l)=pptl(l,i)-pptl(l,j)
                    enddo
                    p(5)=sqrt(p(4)**2-p(1)**2-p(2)**2-p(3)**2)
                    ptij=sqrt(p(1)**2+p(2)**2)
                    call utlob2(1,p(1),p(2),p(3),p(4),p(5)
     .              ,q(1),q(2),q(3),q(4),103)
                    qq=sqrt(q(1)**2+q(2)**2+q(3)**2)
                    if(qq.le.0.075)then
                      !k_t=0.5*( pti + ptj )
                      k_t=0.5*ptij
                      if(kkk.eq.5.or.kkk.eq.6)then
                        v(1)=0
                        v(2)=0
                        v(3)=p(3)/p(4)
                        v(4)=1
                        v(5)=sqrt(1-v(3)**2)
                        call utlob2(1,v(1),v(2),v(3),v(4),v(5)
     .                 ,x(1),x(2),x(3),x(4),101)
                        call utlob2(1,v(1),v(2),v(3),v(4),v(5)
     .                 ,p(1),p(2),p(3),p(4),102)
                        qt=sqrt(p(1)**2+p(2)**2)
                        cph = p(1)/qt
                        sph = p(2)/qt
                        Rout_diff =  x(1)*cph + x(2)*sph;
                        Rside_diff = x(2)*cph - x(1)*sph;
                        Rlong_diff = x(3)
                        x(1) = Rout_diff / sqrt(2.)
                        x(2) = Rside_diff/ sqrt(2.)
                        x(3) = Rlong_diff/ sqrt(2.)
                        x(4) = x(4)      / sqrt(2.)
                      endif
                      do kt=1,4
                      if(fkt(kt-1).lt.k_t.and.k_t.lt.fkt(kt))then
                      do m=1,4
                      del=x(m)
                      call addSou(32+(kkk-3)*16+(kt-1)*4+m,0,0,del,1.0)
                      enddo
                      endif
                      enddo
                    endif
                  endif
                endif
                endif
                endif
              enddo
              endif
            endif
          endif
          endif
          endif
        enddo
      endif
      if(g)then
        xspamx=int(xspamax)
        dx=2*xspamax/kspamax
        do kt=1,4
        do m=1,4
        write(ifhi,'(a)')   '!##################################'
        write(ifhi,'(a,i3)')'!   xxSource     '
        write(ifhi,'(a)')   '!##################################'
        write(ifhi,'(a)') '!newpage'
        write(ifhi,'(a,4i1)')'openhisto htyp lin name xxSource-'
     .   ,jcentr,kkk,kt,m
        write(ifhi,'(a)')'xmod lin ymod lin'
        write(ifhi,'(a,2f8.2)')'xrange ',-xspamx,xspamx
        if(kkk.eq.1)write(ifhi,b)'txt "title prim ptls"'
        if(kkk.eq.2)write(ifhi,b)'txt "title all ptls"'
        if(kkk.eq.3)write(ifhi,b)'txt "title prim pairs"'
        if(kkk.eq.4)write(ifhi,b)'txt "title all pairs"'
        if(kkk.eq.5)write(ifhi,b)'txt "title prim pairs in LCMS"'
        if(kkk.eq.6)write(ifhi,b)'txt "title all pairs in LCMS"'
        write(ifhi,'(a,i1,a)')'txt "xaxis x',m,' (fm)"'
        no=0
        if(kkk.lt.3.and.kt.eq.2)then
         no=(kkk-1)*4+m
         write(ifhi,'(a,i1,a)')'txt "yaxis mean p',m,' (GeV/c)"'
         write(ifhi,'(a)')'yrange auto auto'
        else
         write(ifhi,'(a,i1,a)')'txt "yaxis dn/dx',m,' (1/fm)"'
        endif
        write(ifhi,'(a,i2,a)')'text 0.83 0.75 "J',jcentr,'"'
        if(kkk.ge.3)then
         write(ifhi,'(a,i2,a)')'text 0.05 0.75 "iKt',kt,'"'
        else
         if(kt.eq.1)write(ifhi,'(a)')'text 0.05 0.75 "pi+"'
         if(kt.eq.2)write(ifhi,'(a)')'text 0.05 0.75 "pi+"'
         if(kt.eq.3)write(ifhi,'(a)')'text 0.05 0.75 "K+"'
         if(kt.eq.4)write(ifhi,'(a)')'text 0.05 0.75 "p"'
        endif
        write(ifhi,*)'histoweight ',zkhits(kcentr)
        write(ifhi,'(a)')'array 2'
        sum=0
        do k=1,kspamax
          call addSou((kkk-1)*16+(kt-1)*4+m,  k,no,xi,y)
          sum=sum+y*2*xspamax/kspamax
          write(ifhi,'(2e11.3)')xi,y
        enddo
        write(ifhi,'(a)') 'endarray '
        if(m.eq.1)write(ifhi,'(a,f8.2,a)')'text 0.65 0.50 "',sum,'"'
        write(ifhi,'(a)') 'closehisto plot 0'
        enddo
        enddo
      endif
      enddo
      end

c----------------------------------------------------------------------
      subroutine addSou(i,k,no,x,y)
c----------------------------------------------------------------------
#include "aaa.h"
#include "ems.h"
#include "ho.h"
      common/kcen/kcentr,kmxcentr
      common/czkhits/zkhits(0:ncenthxx)
      parameter (mxdiag=96,kspamaxx=40)
      common/cti/xspamax,kspamax
      common/caddSou/yspa(mxdiag,0:5,kspamaxx)
      common/clogi/icalcy,iaddmo,integ
      integer nrm(8,0:5,kspamaxx)
      data  iaddSou /0/
      save iaddSou,nrm
      iaddSou=iaddSou+1
      if(iaddSou.eq.1)then
        do ix=1,mxdiag
        do jx=0,5
        do kx=1,kspamax
        yspa(ix,jx,kx)=0
        enddo
        enddo
        enddo
        do ix=1,8
        do jx=0,5
        do kx=1,kspamax
        nrm(ix,jx,kx)=0
        enddo
        enddo
        enddo
      endif
      xspamin=-xspamax
      if(iaddmo.eq.1)then
       kk=1+(x-xspamin)/(xspamax-xspamin)*kspamaxx
       if(kk.ge.1.and.kk.le.kspamax)then
        if(no.gt.0)
     .  nrm(no,kcentr,kk)=nrm(no,kcentr,kk)+1
        yspa(i,kcentr,kk)=yspa(i,kcentr,kk)+y
       endif
      elseif(iaddmo.eq.2)then
       anrm=zkhits(kcentr)*nfreeze
       if(no.gt.0)anrm=nrm(no,kcentr,k)
       y=0
       x=xspamin+(k-0.5)*(xspamax-xspamin)/kspamaxx
       dx=2*xspamax/kspamax
       if(anrm.gt.0.)y=yspa(i,kcentr,k)/anrm/dx
      endif
      end

c----------------------------------------------------------------------
      subroutine addHo(i,k,l,m,y)
c----------------------------------------------------------------------
#include "aaa.h"
#include "ems.h"
#include "ho.h"
      common/jcen/jcentr,jmxcentr
      common/caddho/yho(44,0:5,8,55,5)
      common/clogi/icalcy,iaddmo,integ
      common/czhits/zhits(0:ncenthxx),zfhits(0:ncenthxx)
      data  iaddho /0/
      save iaddho
      iaddho=iaddho+1
      if(iaddho.eq.1)then
        do ix=1,42
        do jx=0,5
        do kx=1,8
        do lx=1,50
        do mx=1,5
        yho(ix,jx,kx,lx,mx)=0
        enddo
        enddo
        enddo
        enddo
        enddo
      endif
      if(iaddmo.eq.1)then
       yho(i,jcentr,k,l,m)=yho(i,jcentr,k,l,m)+y
      elseif(iaddmo.eq.2)then
       y=0
       if(zfhits(jcentr).gt.0.)y=yho(i,jcentr,k,l,m)/zfhits(jcentr)
      endif
      end

c----------------------------------------------------------------------
      subroutine setHo(i,j,k)
c----------------------------------------------------------------------
      ! i=1 ... calc y
      ! j=1 ... add in addHo
      !   2 ... normalise in addHo
      ! k=1 ... write into xxHo
      common/clogi/icalcy,iaddmo,integ
      icalcy=i
      iaddmo=j
      integ=k
      end

c----------------------------------------------------------------------
      subroutine getHo(i,j,k)
c----------------------------------------------------------------------
      common/clogi/icalcy,iaddmo,integ
      i=icalcy
      j=iaddmo
      k=integ
      end

c----------------------------------------------------------------------
      subroutine getlogi(d,f,g)
c----------------------------------------------------------------------
      common/clogi/icalcy,iaddmo,integ
      logical d,f,g
      common/jcen/jcentr,jmxcentr
#include "ho.h"

      if(icalcy.eq.0)then
        d=.false.
      else
        d=.true.
      endif

      if(iaddmo.eq.0)then
        f=.false.
      else
        f=.true.
      endif

      if(integ.eq.0)then
        g=.false.
      else
        g=.true.
      endif

      if(jcentr.eq.0)then
        g=.false.
      endif

      end

#endif

c----------------------------------------------------------------------
      subroutine getncenthy
c----------------------------------------------------------------------
#include "aaa.h"
#include "ho.h"
      imax=0
      do i=1,100
        if(zclass(3,i).gt.1e-5)imax=i
      enddo
      ncenthy=imax
      end

c----------------------------------------------------------------------
      subroutine getJcentr
c----------------------------------------------------------------------
      common/jcen/jcentr,jmxcentr
      call getcentr(1,jcentr,jmxcentr)
      end

c----------------------------------------------------------------------
      subroutine getKcentr
c----------------------------------------------------------------------
      common/kcen/kcentr,kmxcentr
      call getcentr(2,kcentr,kmxcentr)
      end

c----------------------------------------------------------------------
      subroutine setcentrVar(x)
c----------------------------------------------------------------------
#include "aaa.h"
#include "ho.h"
#include "ems.h"
      common/cirtfile/irtfile
      common/ccc20/icc20  /ciext4/iext4
      
      if(iext4.eq.0)then 
      
      if(izmode.eq.0)then
        im=mod(irtfile-1,20)+1
        if(icc20.ne.1)stop'ERROR 19042014 '
        bimevt=(zclass(1,im)+zclass(2,im))/2
      elseif(izmode.eq.1)then
        bimevt=x
      elseif(izmode.eq.2)then
        nprt(1)=x
      elseif(izmode.eq.3)then
        nglevt=x
      elseif(izmode.eq.4)then
        segevt=x
      else
        stop'in setcentrVar: wrong izmode'
      endif
      
      else !read from ordered tree
      
        bimevt=0
        nprt(1)=0
        nglevt=0
        segevt=0
      
      endif
      
      end

c----------------------------------------------------------------------
      subroutine getJKNcentr
c----------------------------------------------------------------------
#include "aaa.h"
#include "ho.h"
#include "ems.h"
      common/cen/ncentr
      common/jcen/jcentr,jmxcentr /kcen/kcentr,kmxcentr
      kmax=0
      do k=1,100
        if(zclass(3,k).gt.0.)kmax=kmax+1
      enddo
      ncenthy=kmax
      kk=0
      if(izmode.eq.0.or.izmode.eq.1)then
        do k=1,100
         if(zclass(3,k).gt.0.and.
     .    bimevt.ge.zclass(1,k).and.bimevt.le.zclass(2,k))then
          kk=k
         endif
        enddo
      elseif(izmode.eq.2)then
        zp=nprt(1)
        do k=1,100
         if(zclass(3,k).gt.0.and.
     .    zp.ge.zclass(1,k).and.zp.le.zclass(2,k))then
          kk=k
         endif
        enddo
      elseif(izmode.eq.3)then
        zp=nglevt
        do k=1,5
         if(zclass(3,k).gt.0.and.
     .    zp.ge.zclass(1,k).and.zp.le.zclass(2,k))then
          kk=k
         endif
        enddo
        write(ifmt,'(a,i4)')'centrality class',kk
      else
        stop'in getJKNcentr: wrong jzmode\n\n'
      endif
      ncentr=kk

      if(jzmode(1).gt.0)call getJcentr
      if(jzmode(2).gt.0)call getKcentr

      if(jcentr.gt.ncenthxx)stop'getJKNcentr\n\n'

      end

#if !__BS__ && !__TP__

c##################################################################################
c##################################################################################
c##################################################################################
c##################################################################################
c##################################################################################

c------------------------------------------------------------------------------
      subroutine xEmpty
c------------------------------------------------------------------------------
#include "aaa.h"
      logical d,f,g
      call getlogi(d,f,g)
      if(g)write(ifhi,'(a)')       '!##############################'
      if(g)write(ifhi,'(a)')       '!   Empty     '
      if(g)write(ifhi,'(a)')       '!##############################'
      if(g)write(ifhi,'(a)')'openhisto  xrange 0 1'
      if(g)write(ifhi,'(a)')'array 2 0. 0. 1. 0. endarray'
      if(g)write(ifhi,'(a)')'closehisto plot 0'
      end

c----------------------------------------------------------------------
      subroutine athermal(iret)
c----------------------------------------------------------------------
#include "aaa.h"
#include "ho.h"
      real u(3),q(4)
      data ncntmi/0/
      save ncntmi
      ncntmi=ncntmi+1

      if(ncntmi.eq.1)then
      phievt=0
      ranphi=0
      call hynDefineParticles(1)
      endif
      call InitializeHyperbola

      iret=0
      nevt=0
      nptl=0
      io3=ioclude-3
      mz=int(mchem/2.)+1
      mte=13 !T=130
      temp=0.130
      do m=0,nspes+2
        ggstat(m)=ffstat(mte,io3,mz,mz,mz,m)
      enddo
      yie= volu * ggstat(nspes) * ggstat(nspes+2) * kfrout
      np=yie
      if(rangen().le.yie-np)np=np+1
      if(np.gt.0)then
        do n=1,np
          r=rangen()*yie
          m=0
          do while(yie*ggstat(m)/ggstat(nspes).lt.r)
           m=m+1
          enddo
          id=ispes(m)
          am=aspes(1,m)
          iQ=iQspes(m)
          cot=0.
          js=jspes(m)
          x=hynRanBoseFermi(am,cot,js,temp,istat)
          p=x*temp
          e=sqrt(p**2+aspes(1,m)**2)
          u(3)=2.*rangen()-1.
          angle=2.*pi*rangen()
          u(1)=sqrt(1.-u(3)**2)*cos(angle)
          u(2)=sqrt(1.-u(3)**2)*sin(angle)
          q(1)=p*u(1)
          q(2)=p*u(2)
          q(3)=p*u(3)
          q(4)=e
          pt2=q(1)**2+q(2)**2
          ptr=sqrt(pt2)
          amt=sqrt(pt2+am**2)
          pha=sign(1.,q(2))*acos(q(1)/ptr)
          if(pha.lt.0.)pha=pha+2*pi
          rap=sign(1.,q(3))*alog((q(4)+abs(q(3)))/amt)
          call hynFOStore(id,iQ,ptr,pha,am,rap,0.,0.,1.,0.)
        enddo !np
      endif !np.gt.0
      end

c---------------------------------------------------------------------
      subroutine TestThermal
c---------------------------------------------------------------------
      !Testing the generators  hynRanTherm and hynRanBoseFermi
      !call should be placed in the beginning of bas.f
      !----------------------------------------------------
      real bol(10),bol2(10)
      nn=10000000
      do kk=1,2
      tt=.130 
      if(kk.eq.1)am=0.14
      if(kk.eq.2)am=0.94
      a=am/tt
      if(kk.eq.1)fac=.633
      if(kk.eq.2)fac=44.5
      print*,'Mass : ',am
      do i=1,10
      bol(i)=0
      bol2(i)=0
      enddo
      del=1
      do j=1,nn
      x=hynRanTherm(a)
      y=hynRanBoseFermi(am,0.,kk*2-3,tt,1)
      i=1+x/del
      if(i.gt.0.and.i.le.10)bol(i)=bol(i)+1
      i=1+y/del
      if(i.gt.0.and.i.le.10)bol2(i)=bol2(i)+1
      enddo
      do i=1,10
      x= (i-0.5)*del
      write(*,'(4f8.3)')x, bol(i)/nn/del
     .,fac*x**2*exp(-sqrt(x**2+a**2))
     .,bol2(i)/nn/del      
      enddo
      enddo
      stop'in TestThermal'
      end

c---------------------------------------------------------------------
      subroutine TestThermal1b
c---------------------------------------------------------------------
      real bol(10)
      nn=10000000
      do kk=1,2
      tt=.130
      if(kk.eq.1)am=0.14
      if(kk.eq.2)am=0.94
      a=am/tt
      if(kk.eq.1)fac=.633/tt
      if(kk.eq.2)fac=44.5/tt
      print*,'Mass : ',am
      do i=1,10
      bol(i)=0
      enddo
      del=0.2
      do j=1,nn
      x=tt*hynRanTherm(a)
      i=1+x/del
      if(i.gt.0.and.i.le.10)bol(i)=bol(i)+1
      enddo
      do i=1,10
      x= (i-0.5)*del
      xx=x/tt
      write(*,'(4f8.3)')x, bol(i)/nn/del
     .,fac*xx**2*exp(-sqrt(xx**2+a**2))
      enddo
      enddo
      stop'in TestThermal1b'
      end

c---------------------------------------------------------------------
      subroutine TestThermal2
c---------------------------------------------------------------------
      ! compute            dn  
      !              ----------------(y=0)    from hynRanTherm
      !                2 pi mt dmt dy
      ! 
      ! compare with  dn =  fac*exp(-E/T) d3p 
      !                   with d3p = 2*pi*mt*E* dmt*dy (here E = mt)
      !---------------------------------------------------------------
      real bol(10)
      damt=0.2
      dy=0.2
      nn=10000000
      pi=3.14159
      do kk=1,2
      tt=.130
      if(kk.eq.1)am=0.14
      if(kk.eq.2)am=0.94
      a=am/tt
      if(kk.eq.1)fac=.633/tt /tt**2/4/pi
      if(kk.eq.2)fac=44.5/tt /tt**2/4/pi
      print*,'Mass : ',am
      do i=1,10
      bol(i)=0
      enddo
      do j=1,nn
      p=tt*hynRanTherm(a)
      p3=p*(2.*rangen()-1.)
      pt2=p**2-p3**2
      amt=sqrt(am**2+pt2) 
      e=sqrt(p3**2+amt**2)
      y=log((e+p3)/amt)
      if(abs(y).lt.dy/2)then
       x=amt-am
       i=1+x/damt
       if(i.gt.0.and.i.le.10)bol(i)=bol(i)+1
      endif 
      enddo
      do i=1,10
      amtx= (i-0.5)*damt
      amt=amtx+am
      x=amtx 
      e=amt
      write(*,'(4f11.6)')x, bol(i)/nn  / damt / dy / amt / 2 / pi
     .,fac * exp(-e/tt) * amt
      enddo
      enddo
      stop'in TestThermal1b'
      end

      subroutine xEnergy
#include "aaa.h"
#include "ho.h"
      parameter (mspex=54)
      common/ccsum/eesum,hhsum
      mz=int(mchem/2.)+1
      mte=2 !T=130
      eesum=ffstat(mte,1,mz,mz,mz,mspes+2)
      write(ifhi,'(a)')       '!--------------'
      write(ifhi,'(a,i3)')    '!   energy     '
      write(ifhi,'(a)')       '!--------------'
      write(ifhi,'(a)') ' openhisto htyp pgs xmod lin ymod lin '
      write(ifhi,'(a)') ' array 2'
      write(ifhi,'(2e11.3)')eesum*volu,0.005
      write(ifhi,'(a)') ' endarray closehisto plot 0'
      end

#endif

c---------------------------------------------------------------------------------
      subroutine xTests(neta)
c---------------------------------------------------------------------------------
c      use ArrayModule, only: getd
c      use hocoModule, only: epsc
      double precision getHydynEpsc

#include "aaa.h"
#include "ho.h"
      common/jcen/jcentr,jmxcentr

      write(ifhi,'(a)') '!-----------------------------------------!'
      write(ifhi,'(a)') '!               xTests                    !'
      write(ifhi,'(a)') '!-----------------------------------------!'
      eta=etahy(neta)
      do ntau=1,8
      tau=getHydynTauhy(ntau)
      write(ifhi,'(a)') 'openhisto htyp lin xmod lin ymod log '
      write(ifhi,'(a,2f8.2,a)') 'xrange -12 12 yrange 0.001 auto'
      write(ifhi,'(a,f4.1,a)')'text 0.1 0.9  "[c]=',eta,'"'
      write(ifhi,'(a,f4.1,a)')'text 0.5 0.9  "[t]=',tau,'"'
      write(ifhi,'(a,i1,a)')'text .05 0.8  "jcent=',jcentr,'"'
      write(ifhi,'(a)') 'txt  "xaxis x"'
      write(ifhi,'(a)') 'txt  "yaxis [e]"'
      write(ifhi,'(a)') 'array 2'
      do nx=1,nxhy
        x=xminhy+(nx-1)*(xmaxhy-xminhy)/(nxhy-1)
        y=getHydynEpsc(neta,ntau,nx,nyhy/2+1)
        if(y.gt.1e-3)write(ifhi,'(2e11.3)')x,y
      enddo
      write(ifhi,'(a)') ' endarray closehisto plot 0-'
      write(ifhi,'(a)') 'openhisto htyp fkp'
      write(ifhi,'(a)') 'xmod lin ymod log'
      write(ifhi,'(a)') 'xrange -12 12 yrange 0.001 auto '
      write(ifhi,'(a)') 'function 0.001 closehisto plot 0'
      enddo
      end

c------------------------------------------------------------------------------
      subroutine xFreezeOutTauX
c------------------------------------------------------------------------------
#include "aaa.h"
      do neta=-3,3
      eta=neta
      nhis=1
      npl=0
      nplx=0              !was not defined ???????????? ctp20180425
      deleta=1
      eta1=eta-deleta/2
      eta2=eta+deleta/2
      taumax=10
      do n=1,nptl
        if(ityptl(n).ge.60)then
          tau=0
          tau2=xorptl(4,n)**2-xorptl(3,n)**2
          if(tau2.gt.0.)tau=sqrt(tau2)
          if(tau.lt.taumax)then
           rap=
     *     .5*alog((xorptl(4,n)+xorptl(3,n))/(xorptl(4,n)-xorptl(3,n)))
           if(rap.ge.eta1.and.rap.le.eta2)then
           if(abs(xorptl(2,n)).le.2)then
            npl=npl+1
            if(npl.eq.1)then
              if(nplx.gt.1)
     *        write(ifhi,'(a)')       '  endarray closehisto plot 0-'
              write(ifhi,'(a)')       '!-------------------------------'
              write(ifhi,'(a)')       '!   tau-x      '
              write(ifhi,'(a)')       '!-------------------------------'
              write(ifhi,'(a)')       '!newpage'
              write(ifhi,'(a,i1)')    'openhisto name t-r-0-',nhis
              write(ifhi,'(a)')       'htyp prl xmod lin ymod lin'
              write(ifhi,'(a)')       'xrange -10 10'
              write(ifhi,'(a,f5.1)')  'yrange 0 ',taumax
              write(ifhi,'(a)')       'txt  "xaxis x (fm)"'
              write(ifhi,'(a)')       'txt  "yaxis [t] (fm/c)"'
              write(ifhi,'(a,f4.1,a)')'text 0.05 0.75 "[c]=',eta,'"'
              write(ifhi,'(a)')    'text 0.05 0.87 ""-2fm"L#y"L#2fm""'
              write(ifhi,'(a,f4.1,a,f4.1,a)')'text 0.6 0.75 "b='
     *         ,bimevt,'fm"'
              write(ifhi,'(a)')       'array 2'
            endif
            write(ifhi,'(2e13.5)')xorptl(1,n),tau
            if(npl.eq.1000)then
              write(ifhi,'(a)')    '  endarray closehisto plot 0-'
              nhis=nhis+1
              npl=0
            endif
           endif
           endif
          endif
        endif
      enddo
      if(npl.gt.0)write(ifhi,'(a)')  '  endarray closehisto plot 0'
      enddo
      do neta=-3,3
      eta=neta
      nhis=1
      npl=0
      deleta=1
      eta1=eta-deleta/2
      eta2=eta+deleta/2
      taumax=10
      do n=1,nptl
        if(ityptl(n).ge.60)then
          tau=0
          tau2=xorptl(4,n)**2-xorptl(3,n)**2
          if(tau2.gt.0.)tau=sqrt(tau2)
          if(tau.lt.taumax)then
           rap=
     *     .5*alog((xorptl(4,n)+xorptl(3,n))/(xorptl(4,n)-xorptl(3,n)))
           if(rap.ge.eta1.and.rap.le.eta2)then
           rad=0.
           rad2=xorptl(1,n)**2+xorptl(2,n)**2
           if(rad2.gt.0.)rad=sqrt(rad2)
            npl=npl+1
            if(npl.eq.1)then
              if(nplx.gt.1)
     *        write(ifhi,'(a)')       '  endarray closehisto plot 0-'
              write(ifhi,'(a)')       '!-------------------------------'
              write(ifhi,'(a)')       '!   tau-rad      '
              write(ifhi,'(a)')       '!-------------------------------'
              write(ifhi,'(a)')       '!newpage'
              write(ifhi,'(a,i1)')    'openhisto name t-r-0-',nhis
              write(ifhi,'(a)')       'htyp prl xmod lin ymod lin'
              write(ifhi,'(a)')       'xrange 0 10'
              write(ifhi,'(a,f5.1)')  'yrange 0 ',taumax
              write(ifhi,'(a)')       'txt  "xaxis r (fm)"'
              write(ifhi,'(a)')       'txt  "yaxis [t] (fm/c)"'
              write(ifhi,'(a,f4.1,a)')'text 0.05 0.75 "[c]=',eta,'"'
              write(ifhi,'(a)')    'text 0.05 0.87 ""-2fm"L#y"L#2fm""'
              write(ifhi,'(a,f4.1,a,f4.1,a)')'text 0.6 0.75 "b='
     *         ,bimevt,'fm"'
              write(ifhi,'(a)')       'array 2'
            endif
            write(ifhi,'(2e13.5)')rad,tau
            if(npl.eq.1000)then
              write(ifhi,'(a)')    '  endarray closehisto plot 0-'
              nhis=nhis+1
              npl=0
            endif
          endif
          endif
        endif
      enddo
      if(npl.gt.0)write(ifhi,'(a)')  '  endarray closehisto plot 0'
      enddo
      end

c------------------------------------------------------------------------------
      subroutine xFreezeOutTauEta
c------------------------------------------------------------------------------
#include "aaa.h"
      taumax=10
      nhis=1
      npl=0
      do n=1,nptl
        if(ityptl(n).ge.60)then
          tau=0
          tau2=xorptl(4,n)**2-xorptl(3,n)**2
          if(tau2.gt.0.)tau=sqrt(tau2)
          if(tau.lt.taumax)then
            npl=npl+1
            if(npl.eq.1)then
              write(ifhi,'(a)')      '!-------------------------------'
              write(ifhi,'(a)')      '!   tau-eta      '
              write(ifhi,'(a)')      '!-------------------------------'
              write(ifhi,'(a)')      '!newpage'
              write(ifhi,'(a,i1)')   'openhisto name t-eta-',nhis
              write(ifhi,'(a)')      'htyp prl xmod lin ymod lin'
              write(ifhi,'(a)')      'xrange -4 4'
              write(ifhi,'(a,f5.1)')  'yrange 0 ',taumax
              write(ifhi,'(a)')    'txt  "xaxis [c] "'
              write(ifhi,'(a)')    'txt  "yaxis [t] (fm/c)"'
              write(ifhi,'(a,f4.1,a)')'text 0.6 0.75 "b=',bimevt,'fm"'
              write(ifhi,'(a)')       'array 2'
            endif
            eta=
     *     .5*alog((xorptl(4,n)+xorptl(3,n))/(xorptl(4,n)-xorptl(3,n)))
            write(ifhi,'(2e13.5)') eta,tau
            if(npl.eq.1000)then
              write(ifhi,'(a)')    '  endarray closehisto plot 0-'
              nhis=nhis+1
              npl=0
            endif
          endif
        endif
      enddo
      if(npl.ne.0)write(ifhi,'(a)')  '  endarray closehisto plot 0'
      if(npl.eq.0)stop'xFreezeOutTZ: no particles!!!!!            '
      end

c----------------------------------------------------------------------
      subroutine InitializeHyperbola
c----------------------------------------------------------------------
#include "aaa.h"
      double precision tpro,zpro,ttar,ztar,ttaus,detap,detat
      common /cttaus/  tpro,zpro,ttar,ztar,ttaus,detap,detat
      ttaus=tauzer
      etapro=(ypjtl-yhaha)*etafac
      etatar=-yhaha*etafac
      detap=dble(etapro)
      detat=dble(etatar)
      tpro=dcosh(detap)
      zpro=dsinh(detap)
      ttar=dcosh(detat)
      ztar=dsinh(detat)
      end

c----------------------------------------------------------------------
      subroutine amicro(iret)
c----------------------------------------------------------------------
c  microcanonical decay of cluster specified via keu...ket, tecm, volu
c----------------------------------------------------------------------
#include "aaa.h" 
#include "ho.h" 
      integer jc(nflav,2),ic(2)
      data ncntmi/0/
      save ncntmi
      ncntmi=ncntmi+1

      call utpri('amicro',ish,ishini,4)

      if(ncntmi.eq.1)then
      phievt=0
      ranphi=0
      !call hynDefineParticles(1)
      endif
      call InitializeHyperbola


      nevt=0
      nptl=0

      do i=1,6
      jc(i,1)=0
      jc(i,2)=0
      enddo

      if(keu.ge.0)jc(1,1)=keu
      if(ked.ge.0)jc(2,1)=ked
      if(kes.ge.0)jc(3,1)=kes
      if(kec.ge.0)jc(4,1)=kec
      if(keb.ge.0)jc(5,1)=keb
      if(ket.ge.0)jc(6,1)=ket
      if(keu.lt.0)jc(1,2)=-keu
      if(ked.lt.0)jc(2,2)=-ked
      if(kes.lt.0)jc(3,2)=-kes
      if(kec.lt.0)jc(4,2)=-kec
      if(keb.lt.0)jc(5,2)=-keb
      if(ket.lt.0)jc(6,2)=-ket
      idr=0
      do  nf=1,nflav
        do  ij=1,2
          if(jc(nf,ij).ge.10)idr=7*10**8
        enddo
      enddo
      if(idr/10**8.ne.7)then
        call idenco(jc,ic,ireten)
        if(ic(1).eq.0.and.ic(2).eq.0)then
          ic(1)=100000
          ic(2)=100000
        endif
        idr=8*10**8+ic(1)*100+ic(2)/100
        if(ish.ge.5)write(ifch,'(a,i9)')' id:',idr
      else 
        idr=idr
     *       +mod(jc(1,1)+jc(2,1)+jc(3,1)+jc(4,1),10**4)*10**4
     *       +mod(jc(1,2)+jc(2,2)+jc(3,2)+jc(4,2),10**4)
        call idtrbi(jc,ibptl(1,1),ibptl(2,1),ibptl(3,1),ibptl(4,1))
      endif

      nptl=nptl+1
      idptl(nptl)=idr
      pptl(1,nptl)=0
      pptl(2,nptl)=0
      pptl(3,nptl)=0
      pptl(4,nptl)=tecm
      pptl(5,nptl)=tecm
      istptl(nptl)=11
      tivptl(2,nptl)=1.
      nptl=nptl+1
      idptl(nptl)=idr
      pptl(1,nptl)=0
      pptl(2,nptl)=0
      pptl(3,nptl)=0
      pptl(4,nptl)=tecm
      pptl(5,nptl)=tecm
      istptl(nptl)=11
      tivptl(2,nptl)=1.
      iorptl(nptl)=nptl-1
     
      nptlb=nptl
      ip=nptl
      call hnbaaa(8,ip,iret)
      if(iret.ne.0)then
        print*,'STOP in amicro: hnbaaa iret = ',iret
        stop
      endif
      ifrptl(1,nptlb)=nptlb+1
      ifrptl(2,nptlb)=nptl
      x=xorptl(1,nptlb)
      y=xorptl(2,nptlb)
      z=xorptl(3,nptlb)
      t=xorptl(4,nptlb)
      do n=nptlb+1,nptl
        iorptl(n)=nptlb
        jorptl(n)=0
        istptl(n)=0
        ifrptl(1,n)=0
        ifrptl(2,n)=0
        xorptl(1,n)=x 
        xorptl(2,n)=y 
        xorptl(3,n)=z
        xorptl(4,n)=t
        tivptl(1,n)=t
        call idtau(idptl(n),pptl(4,n),pptl(5,n),taugm)
        r=rangen()
        tivptl(2,n)=t+taugm*(-alog(r))
        ityptl(n)=60
      enddo

      call utprix('amicro',ish,ishini,4)
      return
      end

c----------------------------------------------------------------------
      subroutine getMcentr  ! -> mcentrf() -> mmxcentrf()
c----------------------------------------------------------------------
#include "aaa.h"
      common/mcen/mcentr,mmxcentr
      common/nucl3/phi,bimp
      bimevtsave=bimevt
      bimevt=bimp
      call getcentr(3,mcentr,mmxcentr)  ! 3 -> M centrality
      bimevt=bimevtsave
      end
      integer function mcentrf()
      common/mcen/mcentr,mmxcentr
      mcentrf=mcentr
      end
      integer function mmxcentrf()
      common/mcen/mcentr,mmxcentr
      mmxcentrf=mmxcentr
      end

c----------------------------------------------------------------------
      subroutine getcentr(i,icentr,imxcentr) !i: 1 -> I, 2 -> J, 3 -> M
c----------------------------------------------------------------------
#include "aaa.h"
#include "ho.h"
#include "ems.h"
      kk=0
      kkx=0
      x=centrVar(jzmode(i))
      !write(ifmt,*)'++++',i,jzmode(i),nprt(1),ikoevt
      if(jtable(i).gt.0)then
        do k=1,nrclass(jtable(i))
        xmin=rclass(jtable(i),1,k)
        xmax=rclass(jtable(i),2,k)
        !write(ifmt,*)'++++',i,xmin,xmax,x
        if(x.ge.xmin.and.x.le.xmax)kk=k
        enddo
        kkx=nrclass(jtable(i))
      endif
      icentr=kk
      imxcentr=kkx
      !write(ifmt,*)'++++',i,icentr
      end

c----------------------------------------------------------------------
      function centrVar(j) ! C1=bim C2=Npom C8=Nhpom
c----------------------------------------------------------------------
#include "aaa.h"
#include "ho.h"
#include "ems.h"
      if(j.eq.0)then
        x=0
      elseif(j.eq.1)then !C1
        x=bimevt
      elseif(j.eq.2)then !C2
        !x=nprt(1)
        !x=ikoevt
        x=0
        do i=1,nptl
         call getistptl(i,istx)
         if(istx.eq.30.or.istx.eq.31)x=x+1
        enddo
      elseif(j.eq.3)then
        x=nglevt
      elseif(j.eq.4)then
        x=segevt
      elseif(mod(j,10).eq.5)then
        if(iposi.lt.2)stop'in centrVar: wrong iposi\n\n'
        k=j/10
        nbi=nrclass(k)
        if(nbi.ne.3)stop'in centrVar: wrong nbi\n\n'
        eta1=rclass(k,1,1)
        eta2=rclass(k,2,1)
        pt1=rclass(k,1,2)
        pt2=rclass(k,2,2)
        amult=rclass(k,1,3)
        adivi=rclass(k,2,3)
        mul=0
        do i=maproj+matarg+1,nptl
          if(istptl(i).eq.0)then
            pt=pptl(1,i)**2+pptl(2,i)**2
            pp=sqrt(pptl(1,i)**2+pptl(2,i)**2+pptl(3,i)**2)
            if(pt.gt.0.)then
              pt=sqrt(pt)
              eta=sign(1.,pptl(3,i))*alog((pp+abs(pptl(3,i)))/pt)
            else
              eta=1000.
            endif
            if(abs(idptl(i)).ge.100
     $     .and.abs(idptl(i)).lt.10000)then
              call idchrg( 3 ,idptl(i),ch)
              if(abs(ch).gt.0.1)then
                if( eta.ge.eta1.and.eta.le.eta2
     *         .and. pt.ge. pt1.and. pt.le. pt2) mul=mul+1
              endif
            endif
          endif
        enddo
        x=mul*amult/adivi
      elseif(mod(j,10).eq.6)then !number of Pomerons
        npom=0
        do n=1,nptl
        if(istptl(n).eq.30.or.istptl(n).eq.31)npom=npom+1
        enddo
        x=npom
      elseif(mod(j,10).eq.7)then !vzero
        x=fmux(2.8, 5.1, 0., 100., 1., 1., 1., -0.468,0.,0.,0.) 
        x=x*0.10  !rescaled!!! to fit into 0-40 range, as Pomerons
      elseif(j.eq.8)then !C8
        x=0
        do i=1,nptl
         call getistptl(i,istx)
         call getidptl(i,idx)
         if((istx.eq.30.or.istx.eq.31).and.int(idx/1000000).eq.3)x=x+1
        enddo
      else
        stop'in centrVar: wrong j    '
      endif
      centrVar =x
      end




