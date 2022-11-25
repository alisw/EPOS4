C
C  This file is part of EPOS4
C  Copyright (C) 2022 research institutions and authors (See CREDITS file)
C  This file is distributed under the terms of the GNU General Public License version 3 or later
C  (See COPYING file for the text of the licence)
C

c------------------------------------------------------------------------
      subroutine psphwr(iqq,mIB,nIB,mOB,nOB)                       !ala
c------------------------------------------------------------------------
c iqq - hard contribution type
c mIB,nIB,mOB,nOB ..... incoming/outgoing parton ids (0=g, -5, ...,5 for quarks
c------------------------------------------------------------------------
c copy of subroutine pspawr(kj)
#include "aaa.h"
#include "sem.h"
      common/ciptl/iptl
      common /cnemis/nemis(2) !???????????????? to be added in psahot
      common/cidpomr/idpomr
      common/prtdur/idprtx(2)           ! = pprt (defined in rsh.f)
      common/cprtx/nprtjx,pprtx(5,2)      ! defined in rsh.f
      common/cptldur/nptldur(2)
      common/cqqmx/qqmx1,qqmx2
      dummy=mOB+nOB+iqq
      !-----------------------
      !    ist=25 partons 
      !-----------------------
      do ipr=1,2
        nptl=nptl+1
        nptldur(ipr)=nptl
        pptl(1,nptl)=pprtx(1,ipr) 
        pptl(2,nptl)=pprtx(2,ipr)
        pptl(3,nptl)=pprtx(3,ipr)
        pptl(4,nptl)=pprtx(4,ipr)
        pptl(5,nptl)=0.
        idptl(nptl)=idprtx(ipr)
        iorptl(nptl)=0
        jorptl(nptl)=   1 + idpomr  !kw
        zpaptl(1,nptl)=zpaptl(1,iptl)
        zpaptl(2,nptl)=zpaptl(2,iptl)
        do i=1,2
          ifrptl(i,nptl)=0
        enddo
        do i=1,4
          xorptl(i,nptl)=xorptl(i,iptl)
        enddo
        istptl(nptl)=25
      enddo
      qsqptl(nptl-1) = qqmx1
      qsqptl(nptl)   = qqmx2
      ! similar ityptl definitions for 21 partons in psabor (via variable ibo)
      if(mIB.eq.0.and.nIB.eq.0)then !gg
        if(mOB.eq.0.and.nOB.eq.0)then !gggg
          ityptl(nptl-1) = 30
          ityptl(nptl)   = 30
        else                     !ggNOTgg
          ityptl(nptl-1) = 31
          ityptl(nptl)   = 31
        endif
      elseif(mIB.eq.0.and.abs(nIB).le.3
     .   .or.nIB.eq.0.and.abs(mIB).le.3)then !gq
        ityptl(nptl-1) = 32
        ityptl(nptl)   = 32
      elseif(abs(mIB).le.3.and.abs(nIB).le.3)then !qq
        ityptl(nptl-1) = 33
        ityptl(nptl)   = 33
      elseif(mIB.eq.0.and.abs(nIB).gt.3
     .   .or.nIB.eq.0.and.abs(mIB).gt.3)then !gQ
        ityptl(nptl-1) = 34
        ityptl(nptl)   = 34
      elseif(abs(mIB).le.3.and.abs(nIB).gt.3
     .   .or.abs(nIB).le.3.and.abs(mIB).gt.3)then !qQ
        ityptl(nptl-1) = 35
        ityptl(nptl)   = 35
      elseif(abs(mIB).gt.3.and.abs(nIB).gt.3)then !QQ
        ityptl(nptl-1) = 36
        ityptl(nptl)   = 36
      else 
        stop'ERROR 08062020'
      endif
      !if(abs(idptl(nptl)).eq.4.or.abs(idptl(nptl-1)).eq.4)
      !.print*,'25partonsCharm IB OB:',mIB,nIB,mOB,nOB,'    '
      !.,ityptl(nptl-1),ityptl(nptl)
      pt=min(pprtx(1,1)**2+pprtx(2,1)**2,pprtx(1,2)**2+pprtx(2,2)**2)
      ! phi=polar(pprtx(1,2),pprtx(2,2))-polar(pprtx(1,1),pprtx(2,1))
      ! phi=abs(phi)
      ! if(phi.gt.pi)phi=2*pi-phi
      ! write(ifmt,'(2(a,2i5,4f10.2,f14.2,i7/))')
      ! . 'psphwr: 25_1 = ',nptl-1,idprtx(1),(pprtx(k,1),k=1,4),phi/pi,jj,
      ! . '        25_2 = ',nptl  ,idprtx(2),(pprtx(k,2),k=1,4),phi/pi,jj
      call charmwrite(iqq)
      end

c------------------------------------------------------------------------
      subroutine findHeavyQuarkPairs(nmax)   ! Quarkonium-Method B
c------------------------------------------------------------------------
c Looks for ccbar pairs among children of 25 partons
c------------------------------------------------------------------------
      integer id(nmax),jor(nmax),ity(nmax)
      real p(4,nmax),x(4,nmax)
      do kk=1,2
      !- - - - - - - - - - - - - - - - - - - - - - - - - - - 
      if(kk.eq.1)then
       idkk=4
       ida=44 
      elseif(kk.eq.2)then
        idkk=5
        ida=55
      endif
      n=0
      call getnptl(np)
      do i=1,np
        call getistptl(i,ist)
        !if(ist.eq.25)then
        !  call getidptl(i,idx)
        !  if(abs(idx).eq.4)then
        !    print*,'findHeavyQuarkPairs 25 ----------- ',idx
        !  endif
        !endif
        if(ist.eq.21)then  
          call getidptl(i,idx)
          if(abs(idx).eq.idkk)then
            n=n+1
            if(n.gt.nmax)stop'ERROR findHeavyQuarkPairs Increase nmax'
            id(n)=idx
            call getpptl(i,p(1,n),p(2,n),p(3,n),p(4,n),p5dmy) 
            call getxorptl(i,x(1,n),x(2,n),x(3,n),x(4,n))
            !print*,'findHeavyQuarkPairs ', idx 
            !--get ityxx and jorxx  
            jorxx=-1
            ityxx=-1 
            ixx=i
            idxx=idx
            istxx=ist
            call getistptl(ixx-1,istxxx)
            !print*,'CHK    ',ixx,idxx,istxx,istxxx
            do while(istxx.ne.25.and.istxxx.ge.21.and.istxxx.le.28)
              ixx=ixx-1
              call getidptl(ixx,idxx)
              call getistptl(ixx,istxx)
              call getistptl(ixx-1,istxxx)
              !print*,'CHK    ',ixx,idxx,istxx,istxxx
            enddo
            if(istxx.eq.25)then
              call getjorptl(ixx,jorxx)
              call getityptl(ixx,ityxx)
              !print*,'CHK===>',ixx,idxx,istxx,jorxx,ityxx
            endif
            !write(*,'(a,20a,2i5)')'CHK',('-',ji=1,20),jorxx,ityxx
            jor(n)=jorxx
            ity(n)=ityxx
          endif   
        endif
      enddo
      do j=1,n
        do i=1,j-1
          if(id(i).eq.-id(j))then
            distHQ=sqrt((x(1,i)-x(1,j))**2+(x(2,i)-x(2,j))**2)
            if(distHQ.lt.0.01)then
              p1=p(1,i)+p(1,j)
              p2=p(2,i)+p(2,j)
              p3=p(3,i)+p(3,j)
              p4=p(4,i)+p(4,j)
              x1=(x(1,i)+x(1,j))/2
              x2=(x(2,i)+x(2,j))/2
              x3=(x(3,i)+x(3,j))/2
              x4=(x(4,i)+x(4,j))/2
              call writeCharmPairs(p1,p2,p3,p4,x1,x2,x3,x4
     .         ,ida,jor(j),ity(j)) 
            endif
          endif
        enddo
      enddo
      !- - - - - - - - - - - - - - - - - - - - - - - - - - - 
      enddo !kk

      end
      
c------------------------------------------------------------------------
      subroutine writeCharmPairs(p1,p2,p3,p4,x1,x2,x3,x4,ida,jor,ity)   ! Quarkonium-Method B           
c------------------------------------------------------------------------
#include "aaa.h"
      !return !use Method A
      call  idmass(4   , amCq)
      call  idmass(240 , amCd) !D  !lowest mass open charm resonance
      call  idmass(341 , amCs) !Ds*  !open charm resonance with strangeness
      call  idmass(5   , amBq)
      call  idmass(150 , amBm) !B+
      call  idmass(451 , amBd) !Bc*
      call  idmass(441 , amJp)  !J/psi
      call  idmass(448 , amps)  !psi(3770)
      call  idmass(551 , amUp) !Upsilon
      call  idmass(558 , amUs) !Upsilon(4s)
      si=(p4-p3)*(p4+p3)-p1**2-p2**2
      p5=sqrt(max(0.,si))
      si=p5
      pt= sqrt(p1**2+p2**2)
      amt=sqrt(p1**2+p2**2+p5**2)
      rap=sign(1.,p3)*log((p4+abs(p3))/amt)
      id=0
      amJ2=sqrt(4.*amCd**2+4.*qkonia1**2)
      amU2=sqrt(4.*amBd**2+4.*qkonia2**2)  
      amJ1=amJp
      amU1=amUp
      weightSatJpsi=1  
      weightSatUpsi=1  
      if(jor.eq.0)then
        amJ1=amJ1-1.7!2
        amJ2=amJ2+1.7!2
        amU1=0!amU1-2
        amU2=1000!amU2+2
        weightSatJpsi=0.6071!1.0!1.5
        weightSatUpsi=1.0!2.0
      endif
      !if(ida.eq.44)print*,'CHECK-writeCharmPairs',jor,si
      !.,amJ2,'  ',si.lt.amJ2
      if(ida.eq.44.and.si.gt.amJ1.and.si.lt.amJ2)then
        weightJpsi=0.028 !0.017  !0.015    !physics parameter
        factorJpsi=35.7143!58.8235!66.6666  !technical parameter to increase stat, MUST BE THE SAME IN OPTNS (also in HEP) <==========
        weightJpsi=weightJpsi*factorJpsi    
        if(abs(weightJpsi-1.).gt.1e-3)stop'ERROR 211008a' !weight should be 1 by construction
        id=441
        p5=amJp
      elseif(ida.eq.55.and.si.gt.amU1.and.si.lt.amU2)then
        weightUpsi=0.140  !0.120  !physics parameter
        factorUpsi=7.1428 !8.3333 !technical parameter to increase stat, MUST BE THE SAME IN OPTNS (also in HEP) <==========
        weightUpsi=weightUpsi*factorUpsi
        if(abs(weightUpsi-1.).gt.1e-3)stop'ERROR 211008b'  !weight should be 1 by construction
        id=551
        p5=amUp
      endif
      if(id.gt.0)then
        if(id.eq.441)weightSat=weightSatJpsi
        if(id.eq.551)weightSat=weightSatUpsi
        iweightSat=int(weightSat)
        dweightSat=weightSat-iweightSat
        if(dweightSat.lt.0.)write(ifmt,'(2a,f8.4)')'WARNING '
     .  ,'writeCharmPairs dweightSat =',dweightSat
        amt=sqrt(p1**2+p2**2+p5**2) !overwrire amt using new p5
        do kk=1,iweightSat+1
          if(kk.le.iweightSat.or.rangen().lt.dweightSat)then
            nptl=nptl+1
            idptl(nptl)=id
            if(kk.ge.2)then
              phi=2*pi*rangen()
              p1=pt*cos(phi)
              p2=pt*sin(phi)
            endif
            pptl(1,nptl)=p1
            pptl(2,nptl)=p2
            !rap method (p3 method see 3438u9 and older)
            pptl(3,nptl)=amt*cosh(rap)
            pptl(4,nptl)=amt*sinh(rap) 
            pptl(5,nptl)=p5
            iorptl(nptl)=0
            istptl(nptl)=26
            jorptl(nptl)=jor
            xorptl(1,nptl)=x1
            xorptl(2,nptl)=x2
            xorptl(3,nptl)=x3
            xorptl(4,nptl)=x4
            ityptl(nptl)=ity                  
            ifrptl(1,nptl)=0
            ifrptl(2,nptl)=0  
          endif
        enddo
      endif 
      return
      end

c----------------------------------------------------
      subroutine charmwrite(iqq)!Quarkonium-Method A !====> removed (3438v)     
      dummy=iqq
      end !need to be updated if revived
      subroutine charmwr(p1,p2,p3,p4,ida1,ida2,iprod)!Quarkonium-Method A   !====> removed (3438v)     
c      also called from tim
      dummy=p1+p2+p3+p4+ida1+ida2+iprod
      end !need to be updated if revived
c------------------------------------------------------------------------

c----------------------------------------------------------------------------------------
      subroutine cumu(m,n,x)
c----------------------------------------------------------------------------------------
#include "aaa.h"
      parameter (nfjet=5)
      double precision weightb(-nfjet:nfjet,-nfjet:nfjet)
      common /cweightb/ weightb
      double precision x 
      weightb(m,n)=weightb(m,n)+x
      call getBornIn(mInB,nInB) !get stored Born-in flavors 
      if(ish.eq.-1007)then
        if(x.gt.0.d0)write(ifch,*)'in out val',mInB,nInB,m,n,x 
      endif
      end

c----------------------------------------------------------------------------------------
      double precision function psbori(klas,s,t,j,l,ni)               
c----------------------------------------------------------------------------------------
c contribution to the born cross-section:
c
c   dsigmaBorn/d2pt/dy = s/pi * delta(s+t+u) * 2*pi*alpha**2 /s**2 *psbori
c
c klas - reaction class  !---->  klas defined in subroutine getKlasString
c s - c.m. energy squared for the born scattering,
c t - invariant variable for the born scattering |(p1-p3)**2|,
c j - parton type at current end of the ladder (0 - g, 1,-1,2,... - q)
c l - parton type at opposite end of the ladder (0 - g, 1,-1,2,... - q)
c ni - subprocess number
c----------------------------------------------------------------------------------------
#include "aaa.h"
#include "sem.h"
#include "tab.h"
#include "ems.h"
      double precision u,s,t,t1,u1,t2,u2,xm,q2m1,q2m2
      common/chkbor/ichkbor
      data kerr1 /0/ kerr2 /0/
      save kerr1,kerr2
      real crash(2)
      icrash=3
      ichkbor=0

      call getBornIn(mInB,nInB) !get stored Born-in flavors
      call getFactBorn(bofac)  !get stored bornfactor f*f
      !only relevant for final iteration with ff=bofac corresponding to mInB,nInB


      psbori=0d0
      n=abs(ni)
c      write(ifch,*)"psbori",klas,s,t,j,l,ni,factbqq,factbgg

      if(klas.eq.99)return

      if(klas.lt.1.or.klas.gt.klasmax) then
        write(*,*)'-----> wrong choice: klas = ',klas
        crash(icrash)=1. !to force crash
        stop'wrong klas'
      endif

      q2m1=0d0
      q2m2=0d0
      ja=iabs(j)
      la=iabs(l)

      select case (klas)   
      case(1)
        if(ja.gt.3.or.la.gt.3) return
      case(2)
        if(ja.eq.4.and.la.lt.4.or.ja.lt.4.and.la.eq.4)then
          q2m1=dble(qcmass**2)
        else
          return
        endif
      case(3)
        if(ja.eq.5.and.la.lt.4.or.ja.lt.4.and.la.eq.5)then
          q2m1=dble(qbmass**2)
        else
          return
        endif
      case(4:5)                         !for klas=4,5,6,7 not changed q2m1 and q2m2 because t1, t2 etc used  
        if(ja.gt.3.or.la.gt.3) return        !bg making the expression for the cross section more compact
      case(6)                 
        if(ja.ne.4.or.la.ne.4) return
      case(7)
        if(ja.ne.5.or.la.ne.5) return
      case(8)
        if(ja.ne.4.or.la.ne.4) return
        q2m1=dble(qcmass**2)
        q2m2=dble(qcmass**2)
      case(9)
        if(ja.ne.4.or.la.ne.4) return
        q2m1=dble(qcmass**2)
        q2m2=dble(qbmass**2)
      case(10)
        if(ja.ne.5.or.la.ne.5) return
        q2m1=dble(qbmass**2)
        q2m2=dble(qcmass**2)
      case(11)
        if(ja.ne.5.or.la.ne.5) return
        q2m1=dble(qbmass**2)
        q2m2=dble(qbmass**2)
      case(12)
        if(ja.eq.5.and.la.eq.4.or.ja.eq.4.and.la.eq.5) then
          q2m1=dble(qcmass**2)
          q2m2=dble(qcmass**2)
        else
          return
        endif
      end select 

      u=s-t-2d0*q2m1-2d0*q2m2

      !------------------------------------------------------------------
      !in case of heavy quark creation lights->QQbar and QQbar->lights the  
      !cross sections are expressed in terms of t1 and u1 (ccbar) or t2  
      !and u2 (bbbar), t and u are only used in all other cases
      !------------------------------------------------------------------

      if(klas.lt.4.or.klas.gt.7)then
      if(u.le.0.d0) then
        write(*,*)'-----> psbori: u<0 : ',klas,u,s
        write(*,*)'   ',t,2d0*q2m1,2d0*q2m2
        crash(icrash)=1. !to force crash
        stop
      endif
      endif

      if(t.le.0.d0) then
        write(ifmt,'(//a,i2,2(a,e9.3)//)')'ERROR psbori: t<0  klas='
     .  ,klas,'  t=',t,'  s=',s
        stop
      endif

      t1=t+dble(qcmass**2)        !bg1811
      t2=t+dble(qbmass**2)        !bg1811
      u1=s-t1
      u2=s-t2
      xm=q2m1/s/u                 !bg xm used for Ql --> Ql

      ncgg=0                      !bg Used to define thresholds
      ncqq=0
      nbgg=0
      nbqq=0
      if(nbflav.eq.5) then
        if(u2.ge.dble(qbmass)**2.and.
     &                  s.ge.4.d0*dble(qbmass)**2) then
          ncgg=1
          ncqq=1
          nbgg=1
          nbqq=1
        elseif(u1.ge.dble(qcmass)**2.and.
     &                  s.ge.4.d0*dble(qcmass)**2) then
          ncgg=1
          ncqq=1
        endif
      elseif(nbflav.eq.4) then
        if(u1.ge.dble(qcmass)**2.and.
     &                  s.ge.4.d0*dble(qcmass)**2) then
          ncgg=1
          ncqq=1
        endif
      endif

      !-----------------
      select case (klas)
      !-----------------   

      !-----------------------------------------------------------------
      case(1)  ! ll' --> ll'       
      !-----------------------------------------------------------------
        select case(n)
        case(1)
          if(j.eq.0.and.l.eq.0)then                   !gg->gg
            psbori=(3d0-t*u/s**2+s*u/t**2+s*t/u**2)*4.5d0
            if(ni.lt.0)psbori=psbori*fxsplit
c            if(ni.lt.0)psbori=factbgg*psbori
            call cumu( mInB,nInB , psbori*bofac ) 
            !call getIshy(ishy)
            !if(ishy.eq.11.and.bofac.gt.0.)then
            !print*,'BORN1',mInB,nInB,psbori*bofac
            !endif
          elseif(j*l.eq.0)then                        !gq->gq
            psbori=(s**2+u**2)/t**2+(s/u+u/s)/2.25d0
c            if(ni.lt.0)psbori=factbgg*psbori
            call cumu( mInB,nInB , psbori/2*bofac )
            call cumu( nInB,mInB , psbori/2*bofac )
          elseif(j.eq.l)then                          !qq->qq
            psbori=((s**2+u**2)/t**2+(s**2+t**2)/u**2)/2.25d0
     *      -s**2/t/u/3.375d0                                 !bg
c            if(ni.lt.0)psbori=factbgg*psbori
            call cumu( mInB,nInB , psbori*bofac ) 
          elseif(j.eq.-l)then                         !qq~->qq~
            psbori=((s**2+u**2)/t**2+(u**2+t**2)/s**2)/2.25
     *      +u**2/t/s/3.375d0
c            if(ni.lt.0)psbori=factbqq*psbori
            call cumu( mInB,nInB , psbori/2*bofac )
            call cumu( nInB,mInB , psbori/2*bofac )
          else                                        !qq'->qq'
            psbori=(s**2+u**2)/t**2/2.25d0
c            if(ni.lt.0)psbori=factbgg*psbori
            call cumu( mInB,nInB , psbori/2*bofac )
            call cumu( nInB,mInB , psbori/2*bofac )
          endif
        case(2)
          if(j.eq.0.and.l.eq.0)then                   !gg->qq~
            psbori=3d0*((t/u+u/t)/6.0d0-0.375d0*(t*t+u*u)/s**2)
c            if(ni.lt.0)psbori=factbqq*psbori
            do moutB=1,3
              call cumu( moutB,-moutB , psbori/6*bofac )
              call cumu( -moutB,moutB , psbori/6*bofac )
            enddo
          elseif(j.eq.-l)then                         !qq~->q'q'~
            psbori=2d0*(t*t+u*u)/s**2 *4.0d0/9.0d0
c            if(ni.lt.0)psbori=factbqq*psbori
            do moutB=1,3
              call cumu( moutB,-moutB , psbori/6*bofac )
              call cumu( -moutB,moutB , psbori/6*bofac )
            enddo
          else
            psbori=0d0
          endif
        case(3)
          if(j.ne.0.and.j.eq.-l)then                  !qq~->gg
            psbori=(32.d0/27.d0*(t/u+u/t)-(t*t+u*u)/s**2/.375d0)
c            if(ni.lt.0)psbori=factbgg*psbori
            call cumu( 0,0 , psbori*bofac )
          else
            psbori=0.d0
          endif
        case(4)
          !-------------------------------------------------------------
          !n=4 for photon product processes, make e_q**2 =2/9.,
          !   the average value of charge squared for all types of quarks.
          !-------------------------------------------------------------
          if(j.ne.0.and.j.eq.-l)then                   !qq~->g+gamma
            psbori=16d0*(u/t+t/u)/81.d0
          elseif (j*l.eq.0.and.j+l.ne.0) then          !q(q~)g->q(q~)+gamma
            psbori=2d0*(u/s+s/u)/27.d0
          else
            psbori=0d0
          endif
c          if(ni.lt.0)psbori=factbgg*psbori
        case(5)
          if(j.ne.0.and.j.eq.-l)then                   !qq~->gamma+gamma
            psbori=8d0*(t/u+u/t)/81.d0/3.d0 !bg ga
             psbori=0d0 !bg enlever , temporary to avoid hard gamma which produce fragmentation problem in psahot
          else
            psbori=0d0
          endif
c          if(ni.lt.0)psbori=factbgg*psbori
        end select 
      !-----------------------------------------------------------------
      case(2:3)                     ! Ql --> Ql   2:c, 3:b
      !-----------------------------------------------------------------
        if(n.eq.1) then                                   !bg 1811
          if(j*l.eq.0)then                                !Qg->Qg
            psbori=(s**2+u**2)/t**2+(s/u+u/s)/2.25d0
     *      -4d0*q2m1/t+xm*(xm*t**2-t)/.5625d0
     *      +4d0*q2m1*xm
c            if(ni.lt.0)psbori=factbgg*psbori
            call cumu( mInB,nInB , psbori/2*bofac )
            call cumu( nInB,mInB , psbori/2*bofac )
          elseif(ja.lt.4.or.la.lt.4)then                  !Qq->Qq
            psbori=((s-q2m1)**2+(u+q2m1)**2-2d0*q2m1*t)
     &      /t**2/2.25d0                                       !bg1811
c            if(ni.lt.0)psbori=factbgg*psbori
            call cumu( mInB,nInB , psbori/2*bofac )
            call cumu( nInB,mInB , psbori/2*bofac )
          else
            psbori=0d0                                        !bg1811
          endif
        elseif(n.eq.4) then                                   !bg1811 Q(Q~)g->Q(Q~)+gamma
          if (j*l.eq.0.and.j+l.ne.0) then                     !bg ga
            psbori=2d0*(u/s+s/u)/27.d0                        !bg ga  Mass is missing
          else                                                !bg ga
            psbori=0d0                                        !bg ga
          endif                                               !bg ga
        endif                                                 !bg ga
      !-----------------------------------------------------------------
      case(4)               ! llbar --> CCbar
      !-----------------------------------------------------------------
        if(n.eq.2) then
          if(ja.eq.0.and.la.eq.0) then                  !gg --> CCbar
            if(ncgg.gt.0)
     *      psbori=(4.d0/3.d0-3.d0*u1*t1/s**2)/8.d0 *(t1/u1+u1/t1+
     *      4.d0*dble(qcmass)**2*s/u1/t1*(1d0-dble(qcmass)**2*s/u1/t1))
c            if(ni.lt.0)psbori=factbqq*psbori
            call cumu( 4,-4 , psbori/2*bofac )
            call cumu( -4,4 , psbori/2*bofac )
          elseif(ja.gt.0.and.j.eq.-l) then             !qqbar --> CCbar
            if(ncqq.gt.0)
     *      psbori=4.d0/9.d0*((t1**2+u1**2)/s**2+2.d0*qcmass**2/s)
c            if(ni.lt.0)psbori=factbqq*psbori
            call cumu( 4,-4 , psbori/2*bofac )
            call cumu( -4,4 , psbori/2*bofac )
          endif
        endif
      !-----------------------------------------------------------------
      case(5)              ! llbar --> BBbar
      !-----------------------------------------------------------------
        if(n.eq.2) then
          if(ja.eq.0.and.la.eq.0) then                  !gg --> BBbar
            if(nbgg.gt.0)
     *      psbori=(4.d0/3.d0-3.d0*u2*t2/s**2)/8.d0 *(t2/u2+u2/t2+
     *      4.d0*dble(qbmass)**2*s/u2/t2*(1d0-dble(qbmass)**2*s/u2/t2))
c          if(ni.lt.0)psbori=factbqq*psbori   !to decrease Upsilon 
            call cumu( 5,-5 , psbori/2*bofac )
            call cumu( -5,5 , psbori/2*bofac )
          elseif(ja.gt.0.and.j.eq.-l) then             !qqbar --> BBbar
            if(nbqq.gt.0)
     *      psbori=4.d0/9.d0*((t2**2+u2**2)/s**2+2.d0*qbmass**2/s)
c          if(ni.lt.0)psbori=factbqq*psbori   !to decrease Upsilon 
            call cumu( 5,-5 , psbori/2*bofac )
            call cumu( -5,5 , psbori/2*bofac )
          endif
        endif
      !-----------------------------------------------------------------
      case(6)          ! CCbar --> qqbar , CCbar --> gg 
      !-----------------------------------------------------------------
        if(n.eq.2.and.j.eq.-l)then
          psbori=4.d0/9.d0*((t1**2+u1**2)/s**2+2.d0*qcmass**2/s)*3d0  !bg factor 3 because of the sum u,d,s outgoing light quarks flavor
c          if(ni.lt.0)psbori=factbqq*psbori
          do moutB=1,3
            call cumu( moutB,-moutB , psbori/6*bofac )
            call cumu( -moutB,moutB , psbori/6*bofac )
          enddo
        elseif(n.eq.3.and.j.eq.-l) then
          psbori=(4.d0/3.d0-3.d0*u1*t1/s**2)/8.d0 *(t1/u1+u1/t1+
     *      4.d0*dble(qcmass)**2*s/u1/t1*(1d0-dble(qcmass)**2*s/u1/t1))
          psbori=psbori*dble(8.**2)/dble(3.**2)                       !bg CCbar --> gg = 8^2/3^2 * gg --> CCbar
c          if(ni.lt.0)psbori=factbgg*psbori
          call cumu( 0,0 , psbori*bofac )
        endif
      !-----------------------------------------------------------------
      case(7)          ! BBbar --> qqbar , BBbar --> gg 
      !-----------------------------------------------------------------
        if(n.eq.2.and.j.eq.-l)then                                   !BBbar --> qqbar
          psbori=4.d0/9.d0*((t2**2+u2**2)/s**2+2.d0*qbmass**2/s)*3d0
c          if(ni.lt.0)psbori=factbqq*psbori
          do moutB=1,3
            call cumu( moutB,-moutB , psbori/6*bofac )
            call cumu( -moutB,moutB , psbori/6*bofac )
          enddo
        elseif(n.eq.3.and.j.eq.-l) then                              !BBbar --> gg
          psbori=(4.d0/3.d0-3.d0*u2*t2/s**2)/8.d0 *(t2/u2+u2/t2+
     *      4.d0*dble(qbmass)**2*s/u2/t2*(1d0-dble(qbmass)**2*s/u2/t2))
          psbori=psbori*dble(8.**2)/dble(3.**2)
c          if(ni.lt.0)psbori=factbgg*psbori
          call cumu( 0,0 , psbori*bofac )
        endif
      !-----------------------------------------------------------------
      case(8)       ! QQ --> QQ or QQbar --> QQbar   Q=c
      !-----------------------------------------------------------------
        if(n.eq.1) then
          if(j.eq.l) then                                      !bg QQ --> QQ
            psbori=((s-2d0*q2m1)**2+(u+2d0*q2m1)**2)/t**2/2.25d0
     &             -4d0*q2m1/t/2.25d0+
     &             ((s-2d0*q2m1)**2+(t+2d0*q2m1)**2)/u**2/2.25d0
     &             -4d0*q2m1/u/2.25d0-
     &             (s**2-8d0*q2m1*s+8d0*q2m1**2)/t/u/3.375d0
c            if(ni.lt.0)psbori=factbqq*psbori
            call cumu( mInB,nInB , psbori*bofac )
          elseif(j.eq.-l)then                                  !bg QQbar --> QQbar
            psbori=((s-2d0*q2m1)**2+(u+2d0*q2m1)**2)/t**2/2.25d0
     &             -4d0*q2m1/t/2.25d0+
     &             ((u+2d0*q2m1)**2+(t+2d0*q2m1)**2)/s**2/2.25d0
     &             +4d0*q2m1/s/2.25d0+
     &             (u**2+8d0*q2m1*u+8d0*q2m1**2)/t/s/3.375d0
c            if(ni.lt.0)psbori=factbqq*psbori
            call cumu( mInB,nInB , psbori/2*bofac )
            call cumu( nInB,mInB , psbori/2*bofac )
          endif
        endif
      !-----------------------------------------------------------------
      case(9:10)              ! 9: CCbar --> BBbar  10: BBbar --> CCbar
      !-----------------------------------------------------------------
        if(n.eq.2.and.j.eq.-l) then
          psbori=((t+q2m1+q2m2)**2+(u+q2m1+q2m2)**2)/s**2/2.25d0
     &    +2d0*(q2m1+q2m2)/s/2.25d0
c          if(ni.lt.0)psbori=factbqq*psbori
          call cumu( mInB,nInB , psbori/2*bofac )
          call cumu( nInB,mInB , psbori/2*bofac )
        endif
      !-----------------------------------------------------------------
      case(11)     ! QQ --> QQ or QQbar --> QQbar    Q=b
      !-----------------------------------------------------------------
        if(n.eq.1) then
          if(j.eq.l) then                                      !bg QQ --> QQ
            psbori=((s-2d0*q2m1)**2+(u+2d0*q2m1)**2)/t**2/2.25d0
     &             -4d0*q2m1/t/2.25d0+
     &             ((s-2d0*q2m1)**2+(t+2d0*q2m1)**2)/u**2/2.25d0
     &             -4d0*q2m1/u/2.25d0-
     &             (s**2-8d0*q2m1*s+8d0*q2m1**2)/t/u/3.375d0
c            if(ni.lt.0)psbori=factbgg*psbori
            call cumu( mInB,nInB , psbori*bofac )
          elseif(j.eq.-l)then                                  !bg QQbar --> QQbar
            psbori=((s-2d0*q2m1)**2+(u+2d0*q2m1)**2)/t**2/2.25d0
     &             -4d0*q2m1/t/2.25d0+
     &             ((u+2d0*q2m1)**2+(t+2d0*q2m1)**2)/s**2/2.25d0
     &             +4d0*q2m1/s/2.25d0+
     &             (u**2+8d0*q2m1*u+8d0*q2m1**2)/t/s/3.375d0
c            if(ni.lt.0)psbori=factbqq*psbori
            call cumu( mInB,nInB , psbori/2*bofac )
            call cumu( nInB,mInB , psbori/2*bofac )
          endif
        endif
      !-----------------------------------------------------------------
      case(12)              ! CB --> CB
      !-----------------------------------------------------------------
        if(n.eq.1) then
          psbori=((s-q2m1-q2m2)**2+(u+q2m1+q2m2)**2)/t**2/2.25d0
     &    -2d0*(q2m1+q2m2)/t/2.25d0
        endif
c        if(ni.lt.0)psbori=factbgg*psbori
        call cumu( mInB,nInB , psbori/2*bofac )
        call cumu( nInB,mInB , psbori/2*bofac )

      !---------
      end select
      !---------
     
c      write(ifch,*)"out",psbori
c      if(j.eq.0.and.l.eq.0)print *,"out",psbori,factbgg,ni,klas

      if(psbori.ne.psbori) then
        write(*,*)'klas, n, j, l :',klas,n,j,l
        stop'psbori=NAN, sem.f'
      endif
      if(psbori.lt.0d0) then
        write(*,*)'klas, n, j, l :',klas,n,j,l
        stop'psbori<0, sem.f'
      endif
      return
      end



