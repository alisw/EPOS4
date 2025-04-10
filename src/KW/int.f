C
C  This file is part of EPOS4
C  Copyright (C) 2022 research institutions and authors (See CREDITS file)
C  This file is distributed under the terms of the GNU General Public License version 3 or later
C  (See COPYING file for the text of the licence)
C

c-----------------------------------------------------------------------
      subroutine bjinta(ier)
c-----------------------------------------------------------------------
c  fin. state interactions and decays
c-----------------------------------------------------------------------
#include "aaa.h"
      double precision tpro,zpro,ttar,ztar,ttaus,detap,detat
      common/cttaus/tpro,zpro,ttar,ztar,ttaus,detap,detat /ctimel/ntc
      common/col3/ncol,kolpt,ncoli
      double precision ttaun,ttau0,rcproj,rctarg
      common/cttaun/ttaun /cttau0/ttau0 /geom1/rcproj,rctarg
      logical go

      call utpri('bjinta',ish,ishini,4)

      ier=0

      if(ncol.eq.0.and.iappl.eq.2)goto1000
      if(nevt.ne.1.or.ifrade.eq.0)goto1000

      if(iappl.eq.4.or.iappl.eq.9)then
        goto5000
      endif

      ttaus=tauzer
      ttau0=dsqrt(rcproj*rctarg)
      call jtauin   ! initialize hyperbola

      if(iappl.ne.1)goto 5000


c     no-secondary-interactions or parton-ladder-fusion
c     -------------------------------------------------
      if(iorsce.eq.0.and.iorsdf.eq.0.and.iorshh.eq.0
     &      .or.iorsdf.eq.3.or.iorsdf.eq.5.or.iorsdf.eq.6)then
        if(iorsdf.eq.3.or.iorsdf.eq.5.or.iorsdf.eq.6)then
          if(nclean.gt.0.and.nptl.gt.mxptl/5)then
            ! if nptl already very big, clean up useless particles in cptl list.
            !(do not use it when gakstr() is called (some information lost)
            nptli=maproj+matarg+1
            do iii=nptli,nptl
              go=.true.
              if(nclean.eq.1.and.istptl(iii).le.istmax)go=.false.
              if(go.and.mod(istptl(iii),10).ne.0)istptl(iii)=99
            enddo
            nptl0=nptl
            call utclea(nptli,nptl0)
          endif
          nptlbpo=nptl
          if(ish.ge.2)call alist('core and droplets&',nptlbpo+1,nptl)
          call jintpo(iret)   !core and clustering from overlapping Pomerons
          if(iret.eq.1)goto 1001
          if(iret.ne.1000)then
            call decaydroplets(nptlbpo,iret)
            if(iret.eq.1)goto 1001
          endif
        endif
        goto 5000
      else
        stop'bjinta: not supported any more (310305).     '
      endif

5000  continue

      !call hllexx !commented => keep ist=5 for iorsdf.eq.5.and.ispherio+ihlle.eq.0

      call  getSystemType(isys,amassMax,amassAsy)
      if(isys.le.2)call findHeavyQuarkPairs(100)

      nptlbd=nptl

        if(ispherio.eq.1)goto779  !skip decay
        if(iorsdf.eq.5.and.ispherio+ihlle.ne.0)goto779    !skip decay
        if(iorsdf.eq.5.and.ispherio+ihlle.eq.0)goto779    !skip decay  !may be comented
        if(iorsdf.eq.6)goto779    !skip decay
        if(ifrade.eq.0)goto779    !skip decay
        if(idecay.eq.0)goto779    !skip decay

      do i=1,nptl
        if(ityptl(i).ne.61.and.idptl(i).eq.331)call utphiBW(i)  !BW mass smearing
        if(ityptl(i).eq.61.and.idptl(i).eq.333)call utphiBW(i)  !BW mass smearing
      enddo 

      if(ish.ge.2)call aalist('list after decays&',0,0)
      if(iappl.eq.4.or.iappl.eq.7.or.iappl.eq.9)then
        nptli=1
      else
        nptli=maproj+matarg+1
      endif
      np1=nptli
41    np2=nptl
      nptli=np1
      ip=np1-1
      do while (ip.lt.np2)
      ip=ip+1
      call getistptl(ip,istx)
      if(istx.eq.0)then
      call hdecas(ip,iret)
      if(iret.eq.1)goto 1001
      if(iret.eq.-1)goto 42
c remove useless particles if not enough space
      if(nclean.gt.0.and.nptl.gt.mxptl/2)then
        nnnpt=0
        do iii=nptli,ip
          go=.true.
          call getistptl(iii,istx)
          if(nclean.eq.1.and.istx.le.istmax)go=.false.
          if(go.and.mod(istx,10).ne.0)then
            istx=99
            call setistptl(iii,istx)
            nnnpt=nnnpt+1
          endif
        enddo
        call maxsize_get(2, mx2ptl ) 
        if(nnnpt.gt.mx2ptl-nptl)then
          nptl0=nptl
          call utclea(nptli,nptl0)
          np2=np2-nnnpt
          ip=ip-nnnpt
          nptli=ip
        endif
      endif
      endif
42    continue
      enddo
      nptli=max(nptli,np1)
      np1=np2+1
      if(np1.le.nptl)then
      if(ish.ge.2)then
      if(ish.ge.3)call alist('partial list&',0,0)
      do 6 ip=np1,nptl
        call aalist('&',ip,ip)
6     continue
      endif
      goto 41
      endif
  779 continue

      if(ish.ge.3)call alist('complete list&',1,nptl)

      do j=1,nptl
        call getidptl(j,ja)
        call getistptl(j,istj)
        ja=abs(ja) 
        if( (ja.ge.1.and.ja.le.6.or.ja.eq.9) .and. istj.eq.0 )then
          write(ifmt,'(2a)')'WARNING ist=0 parton detected '
     .      ,'(from rare decay probably) ==> check it '
        endif 
      enddo

c     on shell check
c     --------------
c      if(iappl.eq.1)call jresc

1000  continue
      call utprix('bjinta',ish,ishini,4)
      return

1001  continue
      ier=1
      goto 1000

      end

c----------------------------------------------------------------------
      subroutine decaydroplets(nptla,iret)
c----------------------------------------------------------------------
c  cut and paste from jintpo version 3236
c-----------------------------------------------------------------------
c Decay droplets not included in clusters
c----------------------------------------------------------------------
#include "aaa.h"
#include "ems.h"
      common/cdelzet/delzet,delsce /cvocell/vocell,xlongcell
      iret=0
      delzet=0.
      xlongcell=0.
      do mm=1,nptla
        nptlb=nptl
        if(istptl(mm).eq.10)then
          if(ish.ge.5)write(ifch,*)'Decay remaining droplet :',mm
          if(nptlb.gt.mxptl-10)
     &    call utstop('decaydroplets: mxptl too small (2)&')
         !~~~~~~~~~~~~~~~~~~~~~~~
          call hnbaaa(9,mm,iret) ! 9 = remnant droplets
         !~~~~~~~~~~~~~~~~~~~~~~~
          if(iret.eq.0.and.nptl.ne.nptlb)then ! ---successful decay---
            istptl(mm)=istptl(mm)+1
            ifrptl(1,mm)=nptlb+1
            ifrptl(2,mm)=nptl
            t=tivptl(2,mm)
            x=xorptl(1,mm)+(t-xorptl(4,mm))*pptl(1,mm)/pptl(4,mm)
            y=xorptl(2,mm)+(t-xorptl(4,mm))*pptl(2,mm)/pptl(4,mm)
            z=xorptl(3,mm)+(t-xorptl(4,mm))*pptl(3,mm)/pptl(4,mm)
            do 21 n=nptlb+1,nptl
              iorptl(n)=mm
              jorptl(n)=0
              istptl(n)=0
              ifrptl(1,n)=0
              ifrptl(2,n)=0
              radius=0.8*sqrt(rangen())
              phi=2*pi*rangen()
              ti=t
              zi=z
              xorptl(1,n)=x + radius*cos(phi)
              xorptl(2,n)=y + radius*sin(phi)
              xorptl(3,n)=zi
              xorptl(4,n)=ti
              tivptl(1,n)=ti
              call idtau(idptl(n),pptl(4,n),pptl(5,n),taugm)
              r=rangen()
              tivptl(2,n)=ti+taugm*(-alog(r))
              radptl(n)=0.
              dezptl(n)=0.
              itsptl(n)=0
              rinptl(nptl)=-9999
   21       continue
          else                  ! Unsuccessful decay
            if(ish.ge.1)write(ifch,*)
     *         '***** Unsuccessful remnant cluster decay'
     *             ,' --> redo event.'
          endif
        endif
      enddo
      end

c----------------------------------------------------------------------
c      integer function index25second(i)     !removed in 3438e
c----------------------------------------------------------------------

c----------------------------------------------------------------------
      subroutine index2521(i,i25first,i21last)
c----------------------------------------------------------------------
      i25first=0
      i21last=0
      call getSystemABE(ia,ib,e)
      if(ia.eq.0.or.ib.eq.0)return
      call getistptl(i,isti)
      call getityptl(i,ityi)
      if(ityi/10.ne.3)return
      if(isti.eq.31.or.isti.eq.32)return
      if(isti.eq.29)return
      if(isti.eq.26)return
      call getidptl(i,idi)
      if(idi.eq.9995)return !non isol pi0
      if(idi.eq.9996)return !isol pi0
      if(idi.eq.9997)return !non isol photon
      if(idi.eq.9998)return !isol photon
      if(idi.eq.9999)return !jet
      if(mod(abs(idi),100).eq.88.and.abs(idi).gt.9999)return
      if(i.le. nptlpt)return
      call getiorptl(i,iori)
      call getjorptl(i,jori)
      !print*,iori,jori,i,idi,isti
      call getistptl(iori,isto)
  77  continue
      if(iori.le.0)then
        print*,'i id ist ity ior = ',i,idi,isti,ityi,iori
        stop'ERROR index2521  iori LE 0'
      endif
      if(isto.eq.2)return !ghost particles
      if(isto.eq.1)then
        ioriOld=iori
        call getjorptl(iori,joroOld)
        call getidptl(iori,idoOld)
        istoOld=isto
        call getiorptl(iori,iori)
        call getistptl(iori,isto)
        !print*,iori,joroOld,ioriOld,idoOld,istoOld
        goto 77
      endif
      if(isto.ne.29)then
        print*,'i id ist ity ist(iori) = ',i,idi,isti,ityi
     .  ,isto
        call alist('ERROR&',1,i+1)
        stop'ERROR index2521  ist(iori) NE 29'
      endif  
      call getiorptl(iori,ioro)
      call getjorptl(iori,joro)
      call getidptl(iori,ido)
      !print*,ioro,joro,iori,ido,isto
      i1parton=ioro
      i2parton=joro
      istring=iori
      i25=i1parton
      call getiorptl(i25,ior25)
      call getjorptl(i25,jor25)
      call getidptl(i25,id25)
      call getistptl(i25,ist25)
      !print*,ior25,jor25,i25,id25,ist25
  78  i25=i25-1
      call getiorptl(i25,ior25)
      call getjorptl(i25,jor25)
      call getidptl(i25,id25)
      call getistptl(i25,ist25)
      !print*,ior25,jor25,i25,id25,ist25
      if(i25.eq.1)stop'ERROR 18062012b'
      if(ist25.ne.25)goto 78
      call getistptl(i25-1,ist25m1)
      if(ist25m1.ne.25)
     . stop'ERROR 18062012c'
      call getiorptl(i25-1,ior25m1)
      call getiorptl(i1parton,ior1parton)
      if(ior25m1.ne.ior1parton.and.ior25m1.ne.0)
     . stop'ERROR 18062012d'
      !index25second=i25
      i25first=i25-1
      i21last=i2parton
      iii=0
      do k=i25first+2,i21last
        call getistptl(k,istk)
        if(istk.eq.21)then
          iii=iii+1
          if(iii.eq.1)kref=k  
          call getiorptl(k,iork)
          call getiorptl(kref,iorkref)
          if(iork.ne.iorkref)then
            call getiorptl(i25first,ior25first)
            print*,i,istring,i25first,k, iork,ior25first
            stop'####### ERROR 29042017 #######'
          endif
        endif
      enddo
      end

c----------------------------------------------------------------------
      subroutine jet2521(i,ibig,phi,rap,pp)
c----------------------------------------------------------------------
#include "aaa.h"
      ibig=0
      phi=0
      rap=0
      if(istptl(i).ne.21.and.istptl(i).ne.25)return
      pt=max(1e-6, sqrt(pptl(1,i)**2+pptl(2,i)**2))
      if(pt.ge.1.0)ibig=1
      pp=sqrt(pt**2+pptl(3,i)**2)
      phi=polar( pptl(1,i) , pptl(2,i) )
      rap=sign(1.,pptl(3,i))*log((pp+abs(pptl(3,i)))/pt)
      end

cc-----------------------------------------------------------------------
c      subroutine jrad(i,nq,na,jc,rad)
cc-----------------------------------------------------------------------
cc     return hadron radius (data taken from huefner and povh)
cc-----------------------------------------------------------------------
c#include "aaa.h"
c      integer jc(nflav,2),kc(nflav)
c
c      id=iabs(idptl(i))
c      am=pptl(5,i)
c      if(id.lt.10000)then
c       k=mod(id,10)
c      else
c       k=1
c      endif
c      do l=1,nflav
c       kc(l)=iabs(jc(l,1)-jc(l,2))
c      enddo
c
c      if(nq.eq.0)then   ! mesons
c       if(kc(1).eq.0.and.kc(2).eq.0.and.kc(3).eq.0.and.kc(4).eq.0)then
c        if(k.eq.0)then           ! flavor singlet pseudoscalar mesons
c         if(am.ge.0.000)then
c          rad=0.64                 ! pi0
c          if(am.ge.0.500)then
c           rad=0.60                ! eta
c           if(am.ge.0.900)then
c            rad=0.40               ! eta prime
c            if(am.ge.2.900)then
c             rad=0.17              ! eta charm
c            endif
c           endif
c          endif
c         else
c          write(ifch,*)
c     *    'i:',i,' id:',idptl(i),' k:',k,' m:',am
c          write(ifch,*)'jc:',(jc(l,1),l=1,6),(jc(l,2),l=1,6)
c          call utstop('jrad: meson radius not defined&')
c         endif
c        else                    ! flavor singlet vector mesons
c         if(am.ge.0.000)then
c          rad=0.72                ! rho,omega
c          if(am.ge.1.000)then
c           rad=0.46               ! phi
c           if(am.ge.3.000)then
c            rad=0.20              ! J/psi
c           endif
c          endif
c         else
c          write(ifch,*)
c     *    'i:',i,' id:',idptl(i),' k:',k,' m:',am
c          write(ifch,*)'jc:',(jc(l,1),l=1,6),(jc(l,2),l=1,6)
c          call utstop('jrad: meson radius not defined&')
c         endif
c        endif
c       elseif(kc(3).eq.0.and.kc(4).eq.0)then  ! nonstrange, noncharmed
c        if(k.eq.0)then
c         rad=0.64  ! pi
c        else
c         rad=0.72  ! resonances
c        endif
c       elseif(kc(3).ne.0.and.kc(4).eq.0)then  ! strange
c        if(k.eq.0)then
c         rad=0.59  ! kaons
c        else
c         rad=0.68  ! kaon resonances
c        endif
c       else                                   ! charmed
c        write(ifch,*)'i:',i,' id:',idptl(i)
c        call utstop('jrad: radius of meson not defined&')
c       endif
c      else   !baryons
c       if(kc(4).gt.0)then       ! charmed
c        write(ifch,*)
c     *  'i:',i,' id:',idptl(i),' k:',k,' m:',am
c        write(ifch,*)'i:',i,' id:',idptl(i)
c        call utstop('jrad: radius of charmed baryon not defined&')
c       elseif(kc(3).eq.0)then   ! nonstrange
c        if(k.eq.0)then
c         rad=0.82  !nucleons
c        else
c         rad=1.00  !resonances
c        endif
c       elseif(kc(3).eq.1)then   ! strange
c        if(k.eq.0)then
c         rad=0.76  !lambda, sigma
c        else
c         rad=0.93  !resonances
c        endif
c       elseif(kc(3).eq.2)then   ! double strange
c        if(k.eq.0)then
c         rad=0.71  !cascades
c        else
c         rad=0.87  !resonances
c        endif
c       elseif(kc(3).ge.3)then   ! triple strange
c        rad=0.79  !omega
c       else
c        write(ifch,*)
c     *  'i:',i,' id:',idptl(i),' k:',k,' m:',am
c        write(ifch,*)
c     *  'q:',(jc(l,1),l=1,6),' qbar:',(jc(l,2),l=1,6),
c     *  ' |q-qbar|:',(kc(l),l=1,6)
c        call utstop('jrad: should not happen&')
c       endif
cc string fragments with |#q|>3
c       if(na.gt.3)then
c        a=(na/3.)**(1./3.)
c        if(ish.ge.7)then
c         call utmsg('jrad ')
c         write(ifch,*)
c     *   'i:',i,' id:',idptl(i),' k:',k,' m:',am
c         write(ifch,*)
c     *   'q:',(jc(l,1),l=1,6),' qbar:',(jc(l,2),l=1,6),
c     *   ' |q-qbar|:',(kc(l),l=1,6)
c         write(ifch,*)'nq:',nq,' na:',na,' r:',rad,' ar:',a*rad
c         call utmsgf
c        endif
c        rad=rad*a
c       endif
c      endif
c
c      if(ish.ge.7)then
c       write(ifch,*)
c     * 'i:',i,' id:',idptl(i),' k:',k,' m:',am,' rad:',rad
c       write(ifch,*)'jc:',(jc(l,1),l=1,6),(jc(l,2),l=1,6)
c      endif
c
c      return
c      end
c
c-----------------------------------------------------------------------
      subroutine jresc
c-----------------------------------------------------------------------
#include "aaa.h"
      double precision pa(5),pj(5),utdble
      integer ipptl(mxptl)

      call utpri('jresc ',ish,ishini,4)

      iret=0
      nptlpt=maproj+matarg
      np=0
      do i=nptlpt+1,nptl
       if(istptl(i).eq.0
     * .and.idptl(i).lt.10000.and.pptl(5,i).gt.0.01)then
        np=np+1
        ipptl(np)=i
       endif
      enddo
      if(np.lt.2)goto1001
      do ii=1,np
       i=ipptl(ii)
       if(mod(iabs(idptl(i)),10).lt.2)then
        call idmass(idptl(i),ami)
        dm=abs(ami-pptl(5,i))
        if(dm.gt.0.001)then
         ntry=0
1        continue
         ntry=ntry+1
2        jj=1+int(rangen()*np)
         j=ipptl(jj)
         if(ish.ge.4)write(ifch,*)i,j,istptl(j)
         if(j.eq.i)goto2
         if(mod(iabs(idptl(j)),10).lt.2)then
          call idmass(idptl(j),amj)
         else
          amj=pptl(5,j)
         endif
         do l=1,5
          pa(l)=utdble(pptl(l,i))
          pj(l)=utdble(pptl(l,j))
         enddo
         if(ish.ge.4)write(ifch,'(70a1)')('-',l=1,70)
         if(ish.ge.4)write(ifch,11)i,idptl(i),'before:',pa,'want:',ami
         if(ish.ge.4)write(ifch,11)j,idptl(j),'before:',pj,'want:',amj
         call jrescl(pa,dble(ami),pj,dble(amj),iret)
         if(iret.eq.1)then
          if(ntry.le.50)then
           goto1
          else
           goto1001
          endif
         endif
         if(ish.ge.4)write(ifch,11)i,idptl(i),' after:',pa
         if(ish.ge.4)write(ifch,11)j,idptl(j),' after:',pj
         if(ish.ge.4)write(ifch,'(70a1)')('-',l=1,70)
         do l=1,5
          pptl(l,i)=sngl(pa(l))
          pptl(l,j)=sngl(pj(l))
         enddo
        endif
       endif
      enddo
11    format(i5,1x,i5,1x,a,1x,5(d8.2,1x),a,1x,e8.2)

1000  continue
      call utprix('jresc ',ish,ishini,4)
      return

1001  continue
      if(ish.ge.1)then
        write(ifmt,'(a)')'jresc: could not put on shell'
      endif
      goto1000

      end

c-----------------------------------------------------------------------
      subroutine jrescl(p1,am1,p2,am2,iret)
c-----------------------------------------------------------------------
c rescale momenta of two particles such that the masses assume given
c values.
c input:
c   p1, p2: momenta of the two particles
c   am1, am2: desired masses of the two particles
c output:
c   p1, p2: new momenta of the two particles
c-----------------------------------------------------------------------
#include "aaa.h"
      double precision p1(5),p2(5)
     *                ,p1n(5),p2n(5)
     *                ,a1,a2,a12,am1,am2
     *                ,b1,b2,c,d,e,f,g,p,q,r

      call utpri('jrescl',ish,ishini,7)

      iret=0
      a1=p1(5)**2
      a2=p2(5)**2
      a12=p1(4)*p2(4)-p1(3)*p2(3)-p1(2)*p2(2)-p1(1)*p2(1)
      if(a12.le.(a1+a2))then
       if(ish.ge.7)write(ifch,*)'a_12 < a_1 + a_2'
       if(ish.ge.7)write(ifch,*)a12,' < ',a1+a2
c      goto1001
      endif

11    format(5(d9.3,1x))
      if(ish.ge.7)write(ifch,11)p1,a1
      if(ish.ge.7)write(ifch,11)p2,a2
      if(ish.ge.7)write(ifch,*)a12

      c=(a1+a12)/(a2+a12)
      d=(a1-am1**2-a2+am2**2)/(a2+a12)*0.5d0

      e=a1-2d0*a12*c+a2*c**2
      f=2d0*(a1-a12*(c+d)+a2*c*d)
      g=a1-2d0*a12*d+a2*d**2-am1**2

      p=f/e
      q=g/e
      r=p**2-4d0*q

      if(ish.ge.7)write(ifch,*)'c:',c,' d:',d
      if(ish.ge.7)write(ifch,*)'e:',e,' f:',f,' g:',g
      if(ish.ge.7)write(ifch,*)'p:',p,' q:',q,' r:',r
      if(r.lt.0d0)goto1001

      b1=-0.5d0*(p-dsqrt(r))

      b2=b1*c+d

      if(ish.ge.7)write(ifch,*)'b_1:',b1,' b_2:',b2

      do i=1,4
       p1n(i)=(1d0+b1)*p1(i)-b2*p2(i)
       p2n(i)=(1d0+b2)*p2(i)-b1*p1(i)
      enddo

      a1=p1n(4)**2-p1n(3)**2-p1n(2)**2-p1n(1)**2
      a2=p2n(4)**2-p2n(3)**2-p2n(2)**2-p2n(1)**2
      if(a1.gt.0d0.and.a2.gt.0d0)then
       do i=1,4
        p1(i)=p1n(i)
        p2(i)=p2n(i)
       enddo
       p1(5)=dsqrt(a1)
       p2(5)=dsqrt(a2)
       if(ish.ge.7)write(ifch,11)p1,a1
       if(ish.ge.7)write(ifch,11)p2,a2
      else
       goto1001
      endif

      if(p1(4).lt.0..or.p2(4).lt.0.)goto1001

1000  continue
      call utprix('jrescl',ish,ishini,7)
      return

1001  continue
      iret=1
      goto1000
      end

c-----------------------------------------------------------------------
      subroutine jtain(i,x,y,z,t,n,iopt)
c-----------------------------------------------------------------------
c returns intersection (x,y,z,t) of ptl-i-trajectory with exact hyperbola
c input:
c   i: particle number
c   iopt: formation time considered (0) or not (1)
c output:
c   x,y,z,t: 4-vector of intersection point
c   n: exit code
c       n=0: ok
c       n=1: ptl lives later
c       n=2: ptl lives earlier
c       n=9: tiv1>tiv2
c-----------------------------------------------------------------------
#include "aaa.h"

      !call jtainSimple(i,x,y,z,t,n,iopt)
      call jtainCompli(i,x,y,z,t,n,iopt)

      end

c-----------------------------------------------------------------------
      subroutine jtainSimple(i,x,y,z,t,n,iopt)
c-----------------------------------------------------------------------
c returns intersection (x,y,z,t) of ptl-i-trajectory with taus-line.
c  using the approximations:
c     1) exact hyperbola
c     2) particle pointing back to origin
c               =====> ONLY FOR TESTING <======
c-----------------------------------------------------------------------
c input:
c   i: particle number
c   iopt: formation time considered (0) or not (1)
c output:
c   x,y,z,t: 4-vector of intersection point
c   n: exit code
c       n=0: ok
c       n=1: ptl lives later
c       n=2: ptl lives earlier
c       n=9: tiv1>tiv2
c-----------------------------------------------------------------------
#include "aaa.h"
      double precision tpro,zpro,ttar,ztar,ttaus,detap,detat
      common/cttaus/tpro,zpro,ttar,ztar,ttaus,detap,detat
      double precision vv,tt,ti1,ti2,gamz
      double precision spt2m2E,p4,utdble
      common/ctfi/tin,tfi
      vv=sign(min(1.d0,abs(utdble(pptl(3,i)))/utdble(pptl(4,i)))
     &   ,dble(pptl(3,i)))
      if(abs(vv).ge.1.d0)then
        spt2m2E=utdble(pptl(1,i))**2+utdble(pptl(2,i))**2
     &  +utdble(pptl(5,i))**2
        p4=sqrt(utdble(pptl(3,i))**2+spt2m2E)
       !to avoid precision problem, replace abs(p3)/p4 by sqrt(1-(pt2+m2)/E2)
        spt2m2E=min(1.d0,sqrt(spt2m2E)/p4)
        vv=sign(sqrt((1d0+spt2m2E)*(1d0-spt2m2E)),dble(pptl(3,i)))
      endif
      gamz=1d0/sqrt(1d0-vv)/sqrt(1d0+vv)
      tt=gamz*ttaus !assuming particle points back to origin !!!!!
      t=tt
      z=tt*vv
      x=xorptl(1,i)       !+(t-xorptl(4,i))*pptl(1,i)/pptl(4,i)
      y=xorptl(2,i)       !+(t-xorptl(4,i))*pptl(2,i)/pptl(4,i)
      n=0

      ti1=0d0
      if(iopt.eq.0)ti1=utdble(tivptl(1,i))
      if(iopt.eq.1)ti1=utdble(xorptl(4,i))
      ti2=utdble(tivptl(2,i))

      if(ti1.gt.ti2)then
        n=9
      elseif(t.ge.sngl(ti2))then
        n=2
      elseif(t.le.sngl(ti1))then
        n=1
      endif

      end

c-----------------------------------------------------------------------
      subroutine jtainCompli(i,x,y,z,t,n,iopt)
c-----------------------------------------------------------------------
c returns intersection (x,y,z,t) of ptl-i-trajectory with taus-line.
c input:
c   i: particle number
c   iopt: formation time considered (0) or not (1)
c output:
c   x,y,z,t: 4-vector of intersection point
c   n: exit code
c       n=0: ok
c       n=1: ptl lives later
c       n=2: ptl lives earlier
c       n=9: tiv1>tiv2
c-----------------------------------------------------------------------
#include "aaa.h"
      double precision tpro,zpro,ttar,ztar,ttaus,detap,detat
      common/cttaus/tpro,zpro,ttar,ztar,ttaus,detap,detat
      double precision vv,zza,zz,tt,xo3,xo4,ti1,ti2,derr,dd
      double precision ttp,zzp,ttt,zzt,vvt,vvp,spt2m2E,p4,utdble
      common/ctfi/tin,tfi
      double precision ttau0
      common/cttau0/ttau0

      n=0

      tin=0
      tfi=0

      derr=1d-2
      ttp=tpro*ttaus
      zzp=zpro*ttaus
      ttt=ttar*ttaus
      zzt=ztar*ttaus
      vv=sign(min(1.d0,abs(utdble(pptl(3,i)))/utdble(pptl(4,i)))
     &                ,dble(pptl(3,i)))


      if(abs(vv).ge.1.d0)then
        spt2m2E=utdble(pptl(1,i))**2+utdble(pptl(2,i))**2
     &              +utdble(pptl(5,i))**2
        !if(pptl(4,i).le.0.)then
        p4=sqrt(utdble(pptl(3,i))**2+spt2m2E)
        !else
        !p4=dble(pptl(4,i))
        !endif
        !ctp to avoid precision problem, replace abs(p3)/p4 by sqrt(1-(pt2+m2)/E2)
        spt2m2E=min(1.d0,sqrt(spt2m2E)/p4)
        vv=sign(sqrt((1d0+spt2m2E)*(1d0-spt2m2E)),dble(pptl(3,i)))
      endif
      xo3=utdble(xorptl(3,i))
      xo4=utdble(xorptl(4,i))
      zza=xo3-xo4*vv
      if(iopt.eq.0)then
        ti1=utdble(tivptl(1,i))
      elseif(iopt.eq.1)then
        ti1=xo4
      else
        ti1=0
        call utstop("Wrong iopt in jtain !&")
      endif
      ti2=utdble(tivptl(2,i))

      if(ti1.gt.ti2)then
        n=9
        goto1
      endif

      zfi=sngl(xo3+(ti2-xo4)*vv)
      call jtaus(zfi,tzfi,szfi)
      tfi=tzfi
      if(tfi.ge.sngl(ti2))then
        n=2
        goto1
      endif

      zin=sngl(xo3+(ti1-xo4)*vv)
      call jtaus(zin,tzin,szin)
      tin=tzin
      if(tin.le.sngl(ti1))then
        n=1
        goto1
      endif


    1 continue

      !at this point we know : ti1 = start time of trajectory
      !                        ti2 = end time of trajectory

      if(ttaus.le.ttau0)then  !inner part, usually ttau0 = 0

        tt=ttaus
        zz=xo3+(tt-xo4)*vv
        if(tt.lt.ti1.and.n.eq.0)n=1
        if(tt.ge.ti2.and.n.eq.0)n=2
        goto1000
        
      else ! this is the important piece
           ! We compoute the timm tt (and then x,y,z) which correspond
           ! to the intersection of the trajectory with the generalized
           ! hyperbola, see  Phys.Rept. 232 (1993) 87-299, chaper 13

        vvt=zzt/ttt
        vvp=zzp/ttp
        tt=(ttt+(zza-zzt)*vvt)/(1-vv*vvt)
        zz=xo3+(tt-xo4)*vv
        if(zz.le.zzt)then
          if(tt.lt.ti1.and.n.eq.0)n=1
          if(tt.ge.ti2.and.n.eq.0)n=2
          goto1000
        endif
        tt=(ttp+(zza-zzp)*vvp)/(1-vv*vvp)
        zz=xo3+(tt-xo4)*vv
        if(zz.ge.zzp)then
          if(tt.lt.ti1.and.n.eq.0)n=1
          if(tt.ge.ti2.and.n.eq.0)n=2
          goto1000
        endif
        dd=1-vv**2
        if(sngl(dd).eq.0..and.vv.gt.0.)then
          tt=-(ttaus**2+zza**2)/2d0/zza
        elseif(sngl(dd).eq.0..and.vv.lt.0.)then
          tt=(ttaus**2+zza**2)/2d0/zza
        else
          tt=(zza*vv+dsqrt(zza**2+ttaus**2*dd))/dd
        endif
        zz=xo3+(tt-xo4)*vv
        if(tt.lt.ti1.and.n.eq.0)n=1
        if(tt.ge.ti2.and.n.eq.0)n=2
        if(dabs(ttaus**2-(tt+zz)*(tt-zz)).gt.derr*ttaus**2.and.
     *    dabs(ttaus**2-(tt+zz)*(tt-zz)).gt.derr)then
          if(ish.ge.1)then
            call utmsg('jtain')
            write(ifch,*)'*****  ttaus**2 .ne. (tt+zz)*(tt-zz)'
            write(ifch,*)sngl(ttaus**2),sngl((tt+zz)*(tt-zz))
            call utmsgf
          endif
          goto1000
        endif

      endif

1000  continue

      t=sngl(tt)
      z=sngl(zz)
      x=xorptl(1,i)   +(t-xorptl(4,i))*pptl(1,i)/pptl(4,i)
      y=xorptl(2,i)   +(t-xorptl(4,i))*pptl(2,i)/pptl(4,i)
      return
      end

c-----------------------------------------------------------------------
      subroutine jtaix(i,tau,zor,tor,z,t)
c-----------------------------------------------------------------------
c     returns intersection z,t of ptl-i-trajectory with hyperbola h.
c        h: (t-tor)**2-(z-zor)**2=tau**2 .
c        zor, tor double precision.
c-----------------------------------------------------------------------
#include "aaa.h"
      double precision tor,zor,tors,zors,vv,cc,dd,ttau,derr,tt,zz
     .,utdble
      derr=1d-3
      ttau=utdble(tau)
      zors=utdble(xorptl(3,i))-zor
      tors=utdble(xorptl(4,i))-tor
      vv=utdble(pptl(3,i))/utdble(pptl(4,i))
      vv=dmin1(vv,1d0)
      vv=dmax1(vv,-1d0)
      cc=zors-tors*vv
      dd=1d0-vv**2
      dd=dmax1(dd,0d0)
           if(dd.eq.0d0.and.cc.eq.0d0)then
      if(tau.eq.0.)tt=0d0
      if(tau.ne.0.)tt=utdble(ainfin)
      zz=tt
      goto1000
           elseif(dd.eq.0d0)then
      tt=-(ttau**2+cc**2)/2d0/cc/vv
           elseif(dd.lt.1e-8)then
      tt=-(ttau**2+cc**2)/2d0/cc/vv
      call utmsg('jtaix')
      write(ifch,*)'*****  dd = ',dd,'    treated as zero'
      call utmsgf
           else
      tt=(cc*vv+dsqrt(cc**2+ttau**2*dd))
      tt=tt/dd
           endif
      zz=cc+tt*vv
      if(dabs(ttau**2-(tt+zz)*(tt-zz)).gt.derr*ttau**2
     *.and.dabs(ttau**2-(tt+zz)*(tt-zz)).gt.derr
     *.and.tors**2-zors**2.lt.1e6)then
      if(ish.ge.1)then
      call utmsg('jtaix')
      write(ifch,*)'*****  ttau**2 .ne. (tt+zz)*(tt-zz)'
      write(ifch,*)sngl(ttau**2),sngl((tt+zz)*(tt-zz))
      write(ifch,*)'tau,t,z:'
      write(ifch,*)tau,tt,zz
      write(ifch,*)'#,id(ptl):',i,idptl(i)
      write(ifch,*)'zor,tor(str):',zor,tor
      write(ifch,*)'zors,tors,p,e(ptl):'
      write(ifch,*)sngl(zors),sngl(tors),pptl(3,i),pptl(4,i)
      call utmsgf
      endif
      endif
1000  z=sngl(zz+zor)
      t=sngl(tt+tor)
      return
      end

c-----------------------------------------------------------------------
      subroutine jtaug(su,so,g,y)
c-----------------------------------------------------------------------
c  returns g factor and rapidity y for given su, so
c-----------------------------------------------------------------------
#include "aaa.h"
      double precision tpro,zpro,ttar,ztar,ttaus,detap,detat
      common/cttaus/tpro,zpro,ttar,ztar,ttaus,detap,detat
      double precision ttp,zzp,ttt,zzt,ssp,sst,ssu,sso,ss1,ss2,gg
     *,ssav,yyav,hh,utdble
      double precision ttau0
      common/cttau0/ttau0

      ssu=utdble(su)
      sso=utdble(so)

      if(ssu.ge.sso)then
      sso=(ssu+sso)*0.5d0 + utdble(abs(dezzer))*ttaus*0.5d0
      ssu=(ssu+sso)*0.5d0 - utdble(abs(dezzer))*ttaus*0.5d0
      so=real(sso)
      su=real(ssu)
      endif
      if(ssu.ge.sso)then
        write(ifmt,*)ssu,sso,dble(abs(dezzer))*ttaus*0.5d0
        stop'STOP: sr jtaug: ssu.ge.sso'
      endif

      g=1

      if(ttaus.le.ttau0)return

      ttp=tpro*ttaus
      zzp=zpro*ttaus
      ttt=ttar*ttaus
      zzt=ztar*ttaus
      ssp=ttaus*0.5d0*dlog((ttp+zzp)/(ttp-zzp))
      sst=ttaus*0.5d0*dlog((ttt+zzt)/(ttt-zzt))

      ssav=(ssu+sso)/2d0
      yyav=ssav/ttaus
      if(ssav.ge.ssp)yyav=detap
      if(ssav.le.sst)yyav=detat

      gg=0
      if(ssu.lt.sst)gg=gg + dcosh(detat-yyav) * (dmin1(sst,sso)-ssu)
      if(sso.gt.ssp)gg=gg + dcosh(detap-yyav) * (sso-dmax1(ssp,ssu))
      if(ssu.lt.ssp.and.sso.gt.sst)then
      ss1=dmax1(ssu,sst)
      ss2=dmin1(sso,ssp)
      gg=gg+ttaus*( dsinh(ss2/ttaus-yyav)-dsinh(ss1/ttaus-yyav) )
      endif
      gg=gg/(sso-ssu)

      hh=0
      if(ssu.lt.sst)hh=hh + dsinh(detat-yyav) * (dmin1(sst,sso)-ssu)
      if(sso.gt.ssp)hh=hh + dsinh(detap-yyav) * (sso-dmax1(ssp,ssu))
      if(ssu.lt.ssp.and.sso.gt.sst)then
      ss1=dmax1(ssu,sst)
      ss2=dmin1(sso,ssp)
      hh=hh+ttaus*( dcosh(ss2/ttaus-yyav)-dcosh(ss1/ttaus-yyav) )
      endif
      hh=hh/(sso-ssu)

      yyav=yyav+0.5d0*dlog((gg+hh)/(gg-hh))

      gg=0
      if(ssu.lt.sst)gg=gg + dcosh(detat-yyav) * (dmin1(sst,sso)-ssu)
      if(sso.gt.ssp)gg=gg + dcosh(detap-yyav) * (sso-dmax1(ssp,ssu))
      if(ssu.lt.ssp.and.sso.gt.sst)then
      ss1=dmax1(ssu,sst)
      ss2=dmin1(sso,ssp)
      gg=gg+ttaus*( dsinh(ss2/ttaus-yyav)-dsinh(ss1/ttaus-yyav) )
      endif
      gg=gg/(sso-ssu)

      g=sngl(gg)
      y=sngl(yyav)
      return
      end

c-----------------------------------------------------------------------
      subroutine jtaui(s,ts,zs)
c-----------------------------------------------------------------------
c  returns time ts and coord zs corresponding to ttaus and inv. length s
c-----------------------------------------------------------------------

      double precision tpro,zpro,ttar,ztar,ttaus,detap,detat
      common/cttaus/tpro,zpro,ttar,ztar,ttaus,detap,detat
      double precision ttau0
      common/cttau0/ttau0

      double precision ttp,zzp,ttt,zzt,ssp,sst,ss,zeta,utdble

      zs=s
      ts=sngl(ttaus)

      if(ttaus.le.ttau0)return

      ttp=tpro*ttaus
      zzp=zpro*ttaus
      ttt=ttar*ttaus
      zzt=ztar*ttaus
      ssp=ttaus*0.5d0*dlog((ttp+zzp)/(ttp-zzp))
      sst=ttaus*0.5d0*dlog((ttt+zzt)/(ttt-zzt))
      ss=utdble(s)

           if(ss.le.sst)then
      zs=sngl(zzt+ttar*(ss-sst))
      ts=sngl(ttt+(utdble(zs)-zzt)*zzt/ttt)
           elseif(ss.ge.ssp)then
      zs=sngl(zzp+tpro*(ss-ssp))
      ts=sngl(ttp+(utdble(zs)-zzp)*zzp/ttp)
           else
      zeta=ss/ttaus
      ts=sngl(ttaus*dcosh(zeta))
      zs=sngl(ttaus*dsinh(zeta))
           endif

      return
      end

c-----------------------------------------------------------------------
      subroutine jtauin
c-----------------------------------------------------------------------
c initializes equal time hyperbola at ttaus
c-----------------------------------------------------------------------
#include "aaa.h"
      double precision tpro,zpro,ttar,ztar,ttaus,detap,detat
      common/cttaus/tpro,zpro,ttar,ztar,ttaus,detap,detat
      double precision ttau0,rcproj,rctarg
      common/geom1/rcproj,rctarg
      common/cttau0/ttau0

      call utpri('jtauin',ish,ishini,6)

      if(ttaus.gt.ttau0)then
       if(rcproj.gt.1d-10)then
        detap=dmin1(dble((ypjtl-yhaha)*etafac),dlog(ttaus/rcproj))
       else
        detap=dble((ypjtl-yhaha)*etafac)
       endif
       if(rctarg.gt.1d-10)then
        detat=dmax1(dble(-yhaha*etafac),dlog(rctarg/ttaus))
       else
        detat=dble(-yhaha*etafac)
       endif
       tpro=dcosh(detap)
       zpro=dsinh(detap)
       ttar=dcosh(detat)
       ztar=dsinh(detat)
      else
       detap=0d0
       detat=0d0
       tpro=0d0
       zpro=0d0
       ttar=0d0
       ztar=0d0
      endif

      if(ish.ge.6)then
       write(ifch,*)'hyperbola at tau=',ttaus
       write(ifch,*)'r_p:',rcproj,' r_t:',rctarg
       write(ifch,*)'y_p:',detap,' y_t:',detat
       write(ifch,*)'t_p:',tpro,' z_p:',zpro
       write(ifch,*)'t_t:',ttar,' z_t:',ztar
      endif

      call utprix('jtauin',ish,ishini,6)
      return
      end

c-----------------------------------------------------------------------
      subroutine jtaus(z,tz,sz)
c-----------------------------------------------------------------------
c  returns time tz and inv length sz corresponding to ttaus and z
c-----------------------------------------------------------------------
#include "aaa.h"

      !call jtausSimple(z,tz,sz)
      call jtausCompli(z,tz,sz)

      end

c-----------------------------------------------------------------------
      subroutine jtausSimple(z,tz,sz)
c-----------------------------------------------------------------------
c  returns time tz and inv length sz corresponding to ttaus and z
c-----------------------------------------------------------------------
#include "aaa.h"
      double precision tpro,zpro,ttar,ztar,ttaus,detap,detat
      common/cttaus/tpro,zpro,ttar,ztar,ttaus,detap,detat
      double precision zz,tzz,utdble
      zz=utdble(z)
      tzz=dsqrt(ttaus**2+zz**2)
      tz=sngl(tzz)
      sz=sngl(ttaus*0.5d0*dlog((tzz+zz)/(tzz-zz)))
      end

c-----------------------------------------------------------------------
      subroutine jtausCompli(z,tz,sz)
c-----------------------------------------------------------------------
c  returns time tz and inv length sz corresponding to ttaus and z
c-----------------------------------------------------------------------
#include "aaa.h"
      double precision tpro,zpro,ttar,ztar,ttaus,detap,detat
      common/cttaus/tpro,zpro,ttar,ztar,ttaus,detap,detat
      double precision ttau0
      common/cttau0/ttau0

      double precision ttp,zzp,ttt,zzt,zz,tzz,utdble

      tz=sngl(ttaus)
      sz=z

      if(ttaus.le.ttau0)return

      ttp=tpro*ttaus
      zzp=zpro*ttaus
      ttt=ttar*ttaus
      zzt=ztar*ttaus
      zz=utdble(z)

           if(zz.le.zzt)then
      tz=sngl(ttt+(zz-zzt)*zzt/ttt)
      sz=sngl(ttaus*detat+(zz-zzt)/ttar)
           elseif(zz.ge.zzp)then
      tz=sngl(ttp+(zz-zzp)*zzp/ttp)
      sz=sngl(ttaus*detap+(zz-zzp)/tpro)
           else
      if(sngl(ttaus).ge.ainfin)then
      tz=sngl(ttaus)
      sz=0.
      if(ish.ge.1)then
      call utmsg('jtaus ')
      write(ifch,*)'*****  large ttaus; set tz=ttaus, sz=0'
      write(ifch,*)'ttaus=',ttaus,'zz=',zz
      call utmsgf
      endif
      else
      tzz=dsqrt(ttaus**2+zz**2)
      tz=sngl(tzz)
      sz=sngl(ttaus*0.5d0*dlog((tzz+zz)/(tzz-zz)))
      endif
           endif
      return
      end

c-----------------------------------------------------------------------
      subroutine jtaux(t,z,ttaux)
c-----------------------------------------------------------------------
c  returns ttaux (-> tau-line) for given t and z.
c  ttaux: double precision.
c-----------------------------------------------------------------------
      double precision ttaux,tt,zz,rcproj,rctarg,zt1,zp1,zt2,zp2,ttau0
     .,utdble
      common/geom1/rcproj,rctarg
      common/cttau0/ttau0
      double precision tpro,zpro,ttar,ztar,ttaus,detap,detat
      common/cttaus/   tpro,zpro,ttar,ztar,ttaus,detap,detat

      tt=utdble(t)
      zz=utdble(z)

      if(tt.gt.ttau0)then
       zt1=rctarg-tt
       zp1=tt-rcproj
       zt2=ztar/ttar*tt
       zp2=zpro/tpro*tt
       if(zz.le.dmax1(zt1,zt2))then
        if(zt1.gt.zt2)then
         ttaux=rctarg*dsqrt((tt-zz)/(2d0*rctarg-tt-zz))
        else
         ttaux=(ttar*tt-ztar*zz)/(ttar**2-ztar**2)
        endif
       elseif(zz.ge.dmin1(zp1,zp2))then
        if(zp1.lt.zp2)then
         ttaux=rcproj*dsqrt((tt+zz)/(2d0*rcproj-tt+zz))
        else
         ttaux=(tpro*tt-zpro*zz)/(tpro**2-zpro**2)
        endif
       else
        ttaux=dsqrt(tt**2-zz**2)
       endif
      else
       ttaux=tt
      endif

      return
      end

c-----------------------------------------------------------------------
      subroutine xjden1(ii,itau,x,y,rad,o,u)
c-----------------------------------------------------------------------
c ii=0: initialization
c ii=1: determining dense regions in space of individual events
c       x,y,rad: tranverse coordinates and radius of particle i
c       o,u: specifies long range: u < s < o (s: long coordinate)
c ii=2: plot of individual event
c       xpar1: itau ; valid: 1,..,10
c       xpar2: z-range: -xpar2 < z < xpar2
c       xpar3, x-range: -xpar3 < x < xpar3
c       xpar4, y-range: -xpar4 < y < xpar4
c-----------------------------------------------------------------------
#include "aaa.h"
      double precision tpro,zpro,ttar,ztar,ttaus,detap,detat
      common/cttaus/   tpro,zpro,ttar,ztar,ttaus,detap,detat

      if(idensi.ne.1)stop'STOP in xjden1: idensi must be set 1'

      dlcoox=0.5
      dlcooy=0.5

           if(ii.eq.0)then

      do i=1,nzeta
      do j=1,mxcoox
      do k=1,mxcooy
      kdensi(itau,i,j,k)=0
      enddo
      enddo
      enddo

           elseif(ii.eq.1)then

      if(itau.lt.1.or.itau.gt.mxtau)return

      tau=sngl(ttaus)
      zu=u/tau
      zo=o/tau

            do 1 i=1,nzeta
      zi=-nzeta*dlzeta/2-dlzeta/2+i*dlzeta
      if(zu.gt.zi.or.zo.lt.zi)goto1
      do 2 j=1,mxcoox
      xj=-mxcoox*dlcoox/2-dlcoox/2+j*dlcoox
      do 3 k=1,mxcooy
      yk=-mxcooy*dlcooy/2-dlcooy/2+k*dlcooy
      if((x-xj)**2+(y-yk)**2.gt.rad**2)goto3
      kdensi(itau,i,j,k)=1
    3 continue
    2 continue
    1 continue

           elseif(ii.eq.2)then

      itaux=nint(xpar1)
      if(itaux.gt.mxtau)stop'STOP in xjden1: itaux too large'

      iz=nint(xpar2/dlzeta)
      ix=nint(xpar3/dlcoox)
      iy=nint(xpar4/dlcooy)
      if(iz.gt.nzeta/2)stop'STOP in xjden1: zeta-range too large'
      if(ix.gt.mxcoox/2)stop'STOP in xjden1: x-range too large'
      if(iy.gt.mxcooy/2)stop'STOP in xjden1: y-range too large'

      do k=mxcooy/2+1-iy,mxcooy/2+iy
      write(ifhi,'(a,f7.2)')      '! tau: ',tauv(itaux)
      write(ifhi,'(a)')      'openhisto'
      write(ifhi,'(a,2f7.2)')'xrange',-iz*dlzeta,iz*dlzeta
      write(ifhi,'(a,2f7.2)')'yrange',-ix*dlcoox,ix*dlcoox
      write(ifhi,'(a)')      'set ityp2d 3'
      write(ifhi,'(a)') 'txt  "xaxis space-time rapidity [z]"'
      write(ifhi,'(a)') 'txt  "yaxis transverse coordinate x (fm)"'
      write(ifhi,'(a,i4)')   'array2d',2*iz
      do j=mxcoox/2+1-ix,mxcoox/2+ix
      write(ifhi,'(40i2)')    (kdensi(itaux,i,j,k),
     *                        i=nzeta/2+1-iz,nzeta/2+iz)
      enddo
      write(ifhi,'(a)')       '  endarray'
      write(ifhi,'(a)')       'closehisto plot2d'
      enddo

           else

      stop'STOP in xjden1: wrong option'

           endif

      return
      end

c-----------------------------------------------------------------------
      subroutine xjden2(ii,itau,x,y,rad,s)
c-----------------------------------------------------------------------
c ii=0: initialization
c ii=1: determining dense regions in space of individual events
c       x,y,rad: tranverse coordinates and radius of particle i
c       s: long coordinate
c ii=2: plot of individual event
c       xpar1: itau ; valid: 1,..,10
c       xpar2: s-range: -xpar2 < s < xpar2
c       xpar3, x-range: -xpar3 < x < xpar3
c       xpar4, y-range: -xpar4 < y < xpar4
c-----------------------------------------------------------------------
#include "aaa.h"
      double precision tpro,zpro,ttar,ztar,ttaus,detap,detat
      common/cttaus/   tpro,zpro,ttar,ztar,ttaus,detap,detat
      parameter (mxcoos=60)
      common/cdensh/kdensh(matau,mxcoos,mxcoox,mxcooy),ktot(matau)
      character cy*3

      dlcoox=0.5
      dlcooy=0.5
      dlcoos=0.5

           if(ii.eq.0)then

      do i=1,mxcoos
      do j=1,mxcoox
      do k=1,mxcooy
      kdensh(itau,i,j,k)=0
      enddo
      enddo
      enddo
      ktot(itau)=0

           elseif(ii.eq.1)then

      if(itau.lt.1.or.itau.gt.mxtau)return

      tau=sngl(ttaus)
      z=s/tau

      do 1 i=1,mxcoos
      si=-mxcoos*dlcoos/2-dlcoos/2+i*dlcoos
      do 2 j=1,mxcoox
      xj=-mxcoox*dlcoox/2-dlcoox/2+j*dlcoox
      do 3 k=1,mxcooy
      yk=-mxcooy*dlcooy/2-dlcooy/2+k*dlcooy
      if(((x-xj)**2+(y-yk)**2+(z-si)**2).gt.rad**2)goto3
      kdensh(itau,i,j,k)=kdensh(itau,i,j,k)+1
      ktot(itau)=ktot(itau)+1
    3 continue
    2 continue
    1 continue

           elseif(ii.eq.2)then

      itaux=nint(xpar1)
      if(itaux.gt.mxtau)stop'STOP in xjden2: itaux too large'

      is=nint(xpar2/dlcoos)
      ix=nint(xpar3/dlcoox)
      iy=nint(xpar4/dlcooy)
      if(is.gt.mxcoos/2)stop'STOP in xjden2: s-range too large'
      if(ix.gt.mxcoox/2)stop'STOP in xjden2: x-range too large'
      if(iy.gt.mxcooy/2)stop'STOP in xjden2: y-range too large'

      do k=mxcooy/2+1-iy,mxcooy/2+iy
      write(cy,'(f3.1)')-mxcooy*dlcooy/2-dlcooy/2+k*dlcooy
      write(ifhi,'(a)')      'openhisto'
      write(ifhi,'(a,2f7.2)')'xrange',-is*dlcoos,is*dlcoos
      write(ifhi,'(a,2f7.2)')'yrange',-ix*dlcoox,ix*dlcoox
      write(ifhi,'(a)')      'set ityp2d 5'
      write(ifhi,'(a)') 'txt  "xaxis [z] "'
      write(ifhi,'(a)')
     *'txt  "yaxis  x (fm), y='//cy//' fm"'
      write(ifhi,'(a,i4)')   'array2d',2*is
      do j=mxcoox/2+1-ix,mxcoox/2+ix
      do i=mxcoos/2+1-is,mxcoos/2+is
      write(ifhi,'(e11.3)')
     *kdensh(itaux,i,j,k)/dlcooy/dlcoos/dlcoox/ktot(itaux)
      enddo
      enddo
      write(ifhi,'(a)')       '  endarray'
      write(ifhi,'(a)')       'closehisto plot2d'
      enddo

           else

      stop'STOP in xjden2: wrong option'

           endif

      return
      end

cc-----------------------------------------------------------------------
c      subroutine postscript(iii,ii,ic)
cc-----------------------------------------------------------------------
c#include "aaa.h"
c      character*10 color(10)
c      if(iii.eq.0)then
c      ifps=21
c      open(unit=ifps,file='zzz.ps',status='unknown')
c      WRITE(ifps,'(a)') '%!PS-Adobe-2.0'
c      WRITE(ifps,'(a)') '%%Title: tt2.fig'
c      WRITE(ifps,'(a)') '%%Orientation: Portrait'
c      WRITE(ifps,'(a)') '%%BeginSetup'
c      WRITE(ifps,'(a)') '%%IncludeFeature: *PageSize A4'
c      WRITE(ifps,'(a)') '%%EndSetup'
c      WRITE(ifps,'(a)') '%%EndComments'
c      WRITE(ifps,*) '/l {lineto} bind def'
c      WRITE(ifps,*) '/rl {rlineto} bind def'
c      WRITE(ifps,*) '/m {moveto} bind def'
c      WRITE(ifps,*) '/rm {rmoveto} bind def'
c      WRITE(ifps,*) '/s {stroke} bind def'
c      WRITE(ifps,*) '/gr {grestore} bind def'
c      WRITE(ifps,*) '/gs {gsave} bind def'
c      WRITE(ifps,*) '/cp {closepath} bind def'
c      WRITE(ifps,*) '/tr {translate} bind def'
c      WRITE(ifps,*) '/sc {scale} bind def'
c      WRITE(ifps,*) '/sd {setdash} bind def'
c      WRITE(ifps,*) '/sdo {[.01 .05] 0 sd} bind def'
c      WRITE(ifps,*) '/sdf {[1 .0] 0 sd} bind def'
c      WRITE(ifps,*) '/n {newpath} bind def'
c      WRITE(ifps,*) '/slw {setlinewidth } bind def'
c      write(ifps,*) '/srgb {setrgbcolor} bind def'
c      write(ifps,*) '/lgrey      { 0 0.95 0.95 srgb} bind def'
c      write(ifps,*) '/black      { 0 0 0 srgb} bind def'
c      write(ifps,*) '/red        { 1 0 0 srgb} bind def  '
c      write(ifps,*) '/green      { 0 1 0  srgb} bind def  '
c      write(ifps,*) '/blue       { 0 0 1  srgb} bind def  '
c      write(ifps,*) '/yellow     { 1 0.5 0  srgb} bind def  '
c      write(ifps,*) '/turquoise  { 0 1 1  srgb} bind def  '
c      write(ifps,*) '/purple     { 1 0 1  srgb} bind def  '
cc      write(ifps,*) '/  {   srgb} bind def  '
cc      write(ifps,*) '/  {   srgb} bind def  '
c      write(ifps,*) '/ef {eofill} bind def'
c      WRITE(ifps,'(a)') '%%EndProlog'
c      WRITE(ifps,*) 'gsave'
c      WRITE(ifps,*) '/Helvetica findfont 10 scalefont setfont'
c      color(9)='lgrey     '
c      color(1)='black     '
c      color(2)='red       '
c      color(3)='green     '
c      color(4)='blue      '
c      color(7)='yellow    '
c      color(5)='turquoise '
c      color(6)='purple    '
c      np=0
c         elseif(iii.eq.1)then
c      np=np+1
c      write(ifps,'(a,i4)') '%%Page: number ',np
c      write(ifps,'(a)') 'gsave'
c      WRITE(ifps,*) '100 700 tr'
c      scale=0.125
c      WRITE(ifps,*) 1./scale,1./scale,' sc'
c      WRITE(ifps,*) scale/2.,' slw'
c      WRITE(ifps,*) '/Helvetica findfont ',15.*scale
c     &     ,' scalefont setfont'
c      write(ifps,*) color(1),' n ',smin,xmin,' m ( tau:',tau,') show '
c
c      WRITE(ifps,*) '/Helvetica findfont ',5.*scale
c     &     ,' scalefont setfont'
c
c
c      yb=-2.
c      dy=4./12.
c      yb=yb-dy/2
c      do iyb=0,11
c        yb=yb+dy
c        WRITE(ifps,*) 'gsave'
c        WRITE(ifps,*) (xmax-xmin)*1.1*float(int(iyb/4))
c     &       ,-(xmax-xmin)*1.1*mod(iyb,4),' tr'
c        write(ifps,*) ' n ',smin,xmin,' m ',smax,xmin,' l '
c     &       ,smax,xmax,' l ',smin,xmax,' l cp s '
cc.......particles in layer iyb.............
c        do i=1,nptl
c          if(ii.gt.0)then
c              write(ifps,*)  color(mod(i,5)+2)
c     &             ,' n ',u,x-r,' m ',o,x-r,' l '
c     &             ,o,x+r,' l ',u,x+r,' l cp s '
c              write(ifps,*) ' n ',u,x-r,' m (',i,ior,') show '
c          else
c              write(ifps,*) ' n ',s,x,r,0,360,' arc ',color(ic),' s '
c              write(ifps,*) ' n ',s-r,x,' m (',i,io,') show '
c          endif
c 10     enddo
c        write(ifps,*) color(1),' n ',smin,xmin,' m (',yb,') show '
c        WRITE(ifps,*) 'grestore'
c      enddo                    !yb bin
c      write(ifps,'(a)') 'grestore'
c      write(ifps,*) 'showpage'
c          elseif(iii.eq.2)then
c       write(ifps,*) 'gr'
c
c       write(ifps,'(a)') '%%Trailer'
c       write(ifps,'(a,i4)') '%%Pages: ',np
c       write(ifps,'(a)') '%%EOF'
c       close(unit=ifps)
c          endif
c
c      return
c      end
c

c------------------------------------------------------------------------------
      subroutine xtauev(iii)
c------------------------------------------------------------------------------
      jdum=iii
      end
c------------------------------------------------------------------------------
      subroutine wimi
c------------------------------------------------------------------------------
      end
c------------------------------------------------------------------------------
      subroutine wimino
c------------------------------------------------------------------------------
      end
c------------------------------------------------------------------------------
      subroutine xspace(iii)
c------------------------------------------------------------------------------
      jdum=iii
      end
c------------------------------------------------------------------------------
      subroutine wclu
c------------------------------------------------------------------------------
      end
c------------------------------------------------------------------------------
      subroutine wclufi
c------------------------------------------------------------------------------
      end
c------------------------------------------------------------------------------
      subroutine wtime(iii)
c------------------------------------------------------------------------------
      jdum=iii
      end
