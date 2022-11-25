C
C  This file is part of EPOS4
C  Copyright (C) 2022 research institutions and authors (See CREDITS file)
C  This file is distributed under the terms of the GNU General Public License version 3 or later
C  (See COPYING file for the text of the licence)
C

c-----------------------------------------------------------------------
      subroutine conaa(nfr,iret)
c-----------------------------------------------------------------------
c  determines interaction configuration
c-----------------------------------------------------------------------

#include "aaa.h"
#include "ems.h"
#include "sem.h"
#include "par.h"

      common/geom/rmproj,rmtarg,bmax,bkmx
      common/cmean/xmean,ymean
      common/nucl3/phi,bimp
      common /cncl/xproj(mamx),yproj(mamx),zproj(mamx)
     *            ,xtarg(mamx),ytarg(mamx),ztarg(mamx)
      common /cnclxx/xprojx(mamx),yprojx(mamx),zprojx(mamx)
     *            ,xtargx(mamx),ytargx(mamx),ztargx(mamx)
      common/cfacmss/facmss
      common/cninx/ninx
      common/cnparticip/jproj(2,mamx),jtarg(2,mamx),efluct(6,mamx)

      call utpri('conaa ',ish,ishini,4)

      iret=0

c     initialisations
c     ---------------

      vel=tanh(ypjtl-yhaha)+tanh(yhaha)

c     determine phi, bimp, coll, iproj, itarg, x/y/zproj, x/y/ztarg
c     ---------------------------------------------------------------

           if(iokoll.ne.0)then

      koll=abs(iokoll)
      do k=1,koll
        do n=1,4
          coord(n,k)=0.
        enddo
      enddo
      bimp=0
      phi=0
      bij=bkmx*sqrt(rangen())     !koll=iokoll=matarg=maproj fix the number of Pomeron for one given impact parameter (one pair = one Pomeron, all pair with same b)
      do k=1,koll
        bk(k)=bij
        iproj(k)=k
        itarg(k)=k
        phi=2.*pi*rangen()
        xproj(k)=0d0
        yproj(k)=0d0
        zproj(k)=0d0
        lproj(k)=1
        xtarg(k)=bij*cos(phi)
        ytarg(k)=bij*sin(phi)
        bkx(k)=xtarg(k)
        bky(k)=ytarg(k)
        ztarg(k)=0d0
        ltarg(k)=1
        kproj(k,1)=k
        ktarg(k,1)=k
        ltarg3(k)=0
        lproj3(k)=0
        ktarg3(k,1)=0
        kproj3(k,1)=0
        jproj(1,k)=nint(sign(1.,rangen()-rexdif(iclpro)))
        jtarg(1,k)=nint(sign(1.,rangen()-rexdif(icltar)))
        jproj(2,k)=0 !nint(sign(1.,rangen()-zfluct))
        jtarg(2,k)=0 !nint(sign(1.,rangen()-zfluct))
      enddo

           elseif(maproj.eq.1.and.matarg.eq.1)then

      if(bimevt.lt.0.or.nfr.eq.0)then
        b1=bminim
        b2=amin1(bkmx,bmaxim)
        if(b1.gt.b2)call utstop('conaa: bmin > bmax&')
        bimp=sqrt(b1**2+(b2**2-b1**2)*rangen())
        if(nbarray.gt.0)bimp=barray(mod(nrevt,nbarray)+1)
        phi=phimin+rangen()*(phimax-phimin)
      else
        phi=phievt
        bimp=bimevt
      endif
      koll=1
      do n=1,4
        coord(n,1)=0.
      enddo
      bk(1)=bimp
      iproj(1)=1
      itarg(1)=1
      !!!!!!removed  a phi def here, made big problem !!!!!!!!!!!!!!!!!!!!!!!
      bkx(1)=bimp*cos(phi)
      bky(1)=bimp*sin(phi)
      xproj(1)=0.
      yproj(1)=0.
      zproj(1)=0.
      xtarg(1)=0.
      ytarg(1)=0.
      ztarg(1)=0.
      lproj(1)=1
      ltarg(1)=1
      lproj3(1)=1
      ltarg3(1)=1
      kproj3(1,1)=1
      ktarg3(1,1)=1
      kproj(1,1)=1
      ktarg(1,1)=1
      jproj(1,1)=nint(sign(1.,rangen()-rexdif(iclpro)))
      jtarg(1,1)=nint(sign(1.,rangen()-rexdif(icltar)))
      jproj(2,1)=0 !nint(sign(1.,rangen()-zfluct))
      jtarg(2,1)=0 !nint(sign(1.,rangen()-zfluct))

           else

      if(ninx.eq.1.or.iphsd.eq.0)then !generate random positions EPOS+PHSD
        call conxyz('p',mamx,xproj,yproj,zproj,ypjtl-yhaha)
        call conxyz('t',mamx,xtarg,ytarg,ztarg,yhaha)
        if(iphsd.ne.0)then !store positons 
          do i=1,maproj
            xprojx(i)=xproj(i)
            yprojx(i)=yproj(i)
            zprojx(i)=zproj(i)
          enddo
          do j=1,matarg
            xtargx(j)=xtarg(j)
            ytargx(j)=ytarg(j)
            ztargx(j)=ztarg(j)
          enddo
        endif
      else !take the stored positons from the first event
        do i=1,maproj
          xproj(i)=xprojx(i)
          yproj(i)=yprojx(i)
          zproj(i)=zprojx(i)
        enddo
        do j=1,matarg
          xtarg(j)=xtargx(j)
          ytarg(j)=ytargx(j)
          ztarg(j)=ztargx(j)
        enddo
      endif 
      bx=0
      by=0
      if(maproj.gt.0)then
      if(bimevt.lt.0)then
        b1=bminim
        b2=amin1(rmproj+rmtarg,bmaxim)
        if(b1.gt.b2)call utstop('conaa: bmin > bmax&')
        bimp=sqrt(b1**2+(b2**2-b1**2)*rangen())
        if(nbarray.gt.0)bimp=barray(mod(nrevt,nbarray)+1)
        phi=phimin+rangen()*(phimax-phimin)
      else
        phi=phievt
        bimp=bimevt
      endif
      bx=cos(phi)*bimp
      by=sin(phi)*bimp
      endif
      if(jpsi.lt.0)then !modify b
        bx=xtarg(1)
        by=ytarg(1)
      endif
      if(maproj.eq.0)goto1000
      koll=0
      do i=1,maproj
        lproj(i)=0
        lproj3(i)=0
        jproj(1,i)=nint(sign(1.,rangen()-rexdif(iclpro)))
        jproj(2,i)=0 !nint(sign(1.,rangen()-zfluct))
        !dipole:
        ang=2*pi*rangen()
        diproj(i,1)=disize*cos(ang)
        diproj(i,2)=disize*sin(ang)
      enddo
      do j=1,matarg
        ltarg(j)=0
        ltarg3(j)=0
        jtarg(1,j)=nint(sign(1.,rangen()-rexdif(icltar)))
        jtarg(2,j)=0 !nint(sign(1.,rangen()-zfluct))
        !dipole:
        ang=2*pi*rangen()
        ditarg(j,1)=disize*cos(ang)
        ditarg(j,2)=disize*sin(ang)
      enddo
      do 12 i=1,maproj
      do 11 j=1,matarg
        if(jpsi.lt.0.and.ztarg(j).le.ztarg(1))goto11
        bij=sqrt((xproj(i)+bx-xtarg(j))**2+(yproj(i)+by-ytarg(j))**2)
        if(ish.ge.7)write(ifch,*)'i_p:',i,' i_t:',j,' b_ij:',bij
        if(bij.gt.bkmx)goto 11
        koll=koll+1
        if(koll.gt.kollmx)call utstop('conaa: kollmx too small&')
        bk(koll)=bij
        bkx(koll)=xproj(i)+bx-xtarg(j)
        bky(koll)=yproj(i)+by-ytarg(j)
        iproj(koll)=i
        itarg(koll)=j
        lproj(i)=lproj(i)+1
        ltarg(j)=ltarg(j)+1
        kproj(i,lproj(i))=koll
        ktarg(j,ltarg(j))=koll
11    continue
12    continue
 
      do k=1,koll
        do n=1,4
          coord(n,k)=0.
        enddo
      enddo

           endif

      if(ish.ge.3)write(ifch,*)'koll=',koll
      if(koll.eq.0)goto 1001

      call geoglauber

      call setParamsEfluct()

c     determine coord
c     ---------------
      xmean=0
      ymean=0
      do kl=1,koll
      i=iproj(kl)
      j=itarg(kl)
      dist=ztarg(j)-zproj(i)
      coord(1,kl)=(xproj(i)+xtarg(j))*0.5
      coord(2,kl)=(yproj(i)+ytarg(j))*0.5
      xmean=xmean+coord(1,kl)
      ymean=ymean+coord(2,kl)
      coord(3,kl)=(zproj(i)+ztarg(j))*0.5
      coord(4,kl)=dist/vel
      enddo
      if(koll.gt.0)then
        xmean=xmean/koll
        ymean=ymean/koll
        do kl=1,koll
          coord(1,kl)=coord(1,kl)-xmean
          coord(2,kl)=coord(2,kl)-ymean
        enddo
      endif

c      if(iscreen.ne.0)call CalcScrPair(bimp)

c     exit
c     ----
1000  continue
      if(ish.ge.5)then
      write(ifch,*)'ia,x/y/zproj:'
      do mm=1,maproj
      write(ifch,*)mm,xproj(mm),yproj(mm),zproj(mm)
      enddo
      write(ifch,*)'ia,x/y/ztarg:'
      do mm=1,matarg
      write(ifch,*)mm,xtarg(mm),ytarg(mm),ztarg(mm)
      enddo
      write(ifch,*)'iret',iret
      endif
      call utprix('conaa ',ish,ishini,4)
      return

1001  continue !iret=1 causes redo of whole collision
      iret=1
      if(ish.ge.3)then
      write(ifch,*)
      write(ifch,*)'***** subroutine conaa:'
      write(ifch,*)'***** no nucleon pair found --> no interaction'
      write(ifch,*)
      endif
      goto 1000

      end

c-----------------------------------------------------------------------
      subroutine conigr(n,k)
c-----------------------------------------------------------------------
c  determines position of parton pair at position (n,k)
c-----------------------------------------------------------------------
#include "aaa.h"
#include "ems.h"      

      !KW210403
      !These are "effective" Pomeron positons, which may include initial random
      !walks between valence quarks and Pomeron end.
      !Parton positions are  very important for initial energy distribution asymmetry 
      !and therefore for everything related to flow, which is used to fix parameters.
      !Should be completely independent from the Metroplis stuff, and most importantly:
      !parameters should not be changed to improve the "non-hydro" part!!! 
      !It takes weeks to run full hydro to optimize these parameters.   

      !random:
      call ranpospair(xi,yi,ri,  xj,yj,rj)
      i=iproj(k)
      j=itarg(k)
      !dipole:
      isi=1
      if(rangen().le.0.5)isi=-1
      dix=isi*diproj(i,1)*0.5
      diy=isi*diproj(i,2)*0.5
      isi=1 
      if(rangen().le.0.5)isi=-1
      djx=isi*ditarg(j,1)*0.5
      djy=isi*ditarg(j,2)*0.5
      !add random and dipole
      coordpr(1,n,k)=(xi+xj+dix+djx)*0.5+coord(1,k)
      coordpr(2,n,k)=(yi+yj+diy+djy)*0.5+coord(2,k)
      end

c-----------------------------------------------------------------------
      subroutine ranpospair(xi,yi,ri,xj,yj,rj)
c----------------------------------------------------------------------- 
#include "aaa.h"
#include "ems.h"
      nko=0 !variable used in versions <= 2.14
      call ranpos(xi,yi,ri,1)
      call ranpos(xj,yj,rj,2)
      end

c-----------------------------------------------------------------------
      subroutine ranpos(xxx,yyy,rrr,ipt)
c-----------------------------------------------------------------------
c ipt - 1 for proj, 2 for targ
c-----------------------------------------------------------------------
#include "aaa.h"
      common/cbglaub/bglaub
ctp aaa=0.25 fm for a classical proton radius of 1fm ???
ckw these are not the positions of the constituent quarks of the proton, 
ckw  but of the string end partons, which should depend on Z
      if(facpos(ipt).gt.0.)then !~~~~~~exp distribution
        aaa=0.25 !in fm  (exp(-r/aaa) 
        rrr = ( -log(rangen()) -log(rangen()) -log(rangen()) ) * aaa
      elseif(facpos(ipt).lt.0.)then !~~~~~~~gauss distr 
        rrr=sqrt( -2*log(rangen())*(cos(2*pi*rangen()))**2
     .            -2*log(rangen())*(cos(2*pi*rangen()))**2
     .            -2*log(rangen())*(cos(2*pi*rangen()))**2  )
     .     / sqrt(3.)
      else
        print*,ipt,  facpos(ipt)
        stop'03122011' 
      endif !~~~~~~~~
      afacpos=mod(abs(facpos(ipt)),100.)
      rrr=rrr*afacpos
      costhet=2*rangen()-1
      sinthet=sqrt(1-costhet**2)
      phi = 2.*pi*rangen()
      xxx=rrr*sinthet*cos(phi)
      yyy=rrr*sinthet*sin(phi)
      end
      

c-----------------------------------------------------------------------
      subroutine xGeometry(iii)
c-----------------------------------------------------------------------
#include "aaa.h"
#include "ems.h"
#include "sem.h"
      common/xgeom/nnn,naa(kollmx),nbb(kollmx)
      character*5 fmt1,fmt2

      if(iii.eq.0)then

      do k=1,kollmx
        naa(k)=0
      enddo
      nnn=0


      elseif(iii.eq.1)then

      ngl=0
      do k=1,koll
        r=bk(k)
        if(r.le.sqrt(sigine/10./pi))ngl=ngl+1
      enddo
      if(ngl.ne.0)then
        nnn=nnn+1
        naa(ngl)=naa(ngl)+1
      endif

      elseif(iii.eq.2)then

      if(xpar1.eq.0..and.xpar2.eq.0.)then
       print*,'---------- xGeometry -----------------------'
       return
      endif
      x1=1-0.01*xpar2
      x2=1-0.01*xpar1
      kmx=0
      nbb(1)=naa(1)
      do k=2,kollmx
        if(naa(k).ne.0.)kmx=k
        nbb(k)=nbb(k-1)+naa(k)
      enddo
      k1=0
      k2=0
      do k=1,kmx
        x=nbb(k)/float(nnn)
        if(x.lt.x1)k1=k
        if(x.lt.x2)k2=k
      enddo
      k1=k1+1
      k2=k2+1
      x=0
      av=0
      su=0
      p1=0.
      p2=0.
      do k=1,kmx
        xb=x
        x=nbb(k)/float(nnn)
        y=naa(k)/float(nnn)
        dx=x-xb
        p=0.
        if(k.eq.k1)then
          p=(x-x1)/dx
          p1=p
        elseif(k.eq.k2)then
          p=(x2-xb)/dx
          p2=p
        elseif(k.gt.k1.and.k.lt.k2)then
          p=1
        endif
      av=av+y*p*k
      su=su+y*p
      enddo
      av=av/su
      n1=nint(100*p1)
      n2=nint(100*p2)
      if(n1.eq.0)then
        k1=k1+1
        n1=100
      endif
      if(k1.le.9)fmt1='i1,4x'
      if(k1.gt.9.and.k1.le.99)fmt1='i2,3x'
      if(k1.gt.99.and.k1.le.999)fmt1='i3,2x'
      if(k1.gt.999.and.k1.le.9999)fmt1='i4,1x'
      if(k2.le.9)fmt2='i1,4x'
      if(k2.gt.9.and.k2.le.99)fmt2='i2,3x'
      if(k2.gt.99.and.k2.le.999)fmt2='i3,2x'
      if(k2.gt.999.and.k2.le.9999)fmt2='i4,1x'
      write(6,'(i4,a,i5,a,'//fmt1//',i6,a,i5,a,'//fmt2//',5x,a,f8.2)')
     &       nint(xpar2),':MIN',n1,'%',k1
     &      ,nint(xpar1),':MAX',n2,'%',k2   ,'av:',av
      endif

      end

c-----------------------------------------------------------------------
      function conbmx()
c-----------------------------------------------------------------------
      double precision om1intbc,p,eps
#include "aaa.h"
#include "sem.h"

      iomegasave=iomega
      iomega=2
      call DefXminDf(dble(engy**2))

      conbmx=0.
      b1=0.
      b2=7.
c      return
      eps=5.0d-4
      p=1.d0-dexp(-om1intbc(b2))
      if(p.gt.2.d0*eps)goto 999

      ntry=0

10    ntry=ntry+1
      b=b1+.5*(b2-b1)

      p=(1.d0-dexp(-om1intbc(b)))

      if(p.gt.eps)then
       if(p.gt.2.d0*eps)then
        b1=b
       else
        conbmx=b
       goto 999
       endif
      else
       if(p.lt.eps/5.d0)then
        b2=b
       else
        conbmx=b
        goto 999
       endif
      endif

      if(ntry.le.1000)goto 10
      write(ifmt,*)'Too much try in conbmx ... bmax=',b
      conbmx=b

 999  continue
      iomega=iomegasave
      call DefXminDf(dble(engy**2))

      return

      end

c-----------------------------------------------------------------------
      function conbmd()
c-----------------------------------------------------------------------
      double precision om1intbc,p,eps
#include "aaa.h"
#include "sem.h"

      conbmd=0.
      b1=0.
      b2=7.
      eps=0.5

10    b=b1+.5*(b2-b1)
      p=(1.d0-dexp(-om1intbc(b)))

      if(abs(p-eps).lt.0.01)then
        conbmd=b
        return
      elseif(p.gt.eps)then
        b1=b
      elseif(p.lt.eps)then
        b2=b
      endif

      goto 10
      end

c-----------------------------------------------------------------------
      function conbmxndif()
c-----------------------------------------------------------------------
      double precision om1intbc,p,eps
#include "aaa.h"
#include "sem.h"


      iomegasave=iomega
      iomega=2
      b1=0.
      b2=7.
      conbmxndif=b2
      eps=5.d-4
      p=1.d0-dexp(-om1intbc(b2))
      if(p.gt.2.d0*eps)goto 100

      ntry=0

10    ntry=ntry+1
      b=b1+.5*(b2-b1)

      p=(1.d0-dexp(-om1intbc(b)))

      if(p.gt.eps)then
       if(p.gt.2.d0*eps)then
        b1=b
       else
        conbmxndif=b
       goto 100
       endif
      else
       if(p.lt.eps/5.d0)then
        b2=b
       else
        conbmxndif=b
        goto 100
       endif
      endif

      if(ntry.le.1000)goto 10
      write(ifmt,*)'Too much try in conbmxndif ... bkmxndif=',b
      conbmxndif=b
 100  iomega=iomegasave
c      print *,'conbmxdnif',conbmxndif
      return

      end

c-----------------------------------------------------------------------
      function conbmxdif()
c-----------------------------------------------------------------------
c find b to have (1-exp(-om))pmax=pdiff
c-----------------------------------------------------------------------
      double precision om1intbc,pmax,drootom,pdiff
#include "aaa.h"
#include "sem.h"

      conbmxdif=0.

      if(facdif.ge.1.)return
      if(facdif.gt.0)then
        pdiff=dble(facdif)
        corr=1.
      else
        pdiff=0.99d0
        corr=abs(facdif)
      endif

      b1=0.5      
      bmax=7.
      iomegasave=iomega
      iomega=2

      eps=1.e-5
      pmax=1.d0-dexp(-om1intbc(b1))
      if(pmax.gt.eps)then
        conbmxdif=corr*sngl(drootom(pdiff,pmax,bmax,eps))
      endif
c      print *,'conbmxdif',conbmxdif,corr,pmax
      iomega=iomegasave

      return

      end

c-----------------------------------------------------------------------
      subroutine conre
c-----------------------------------------------------------------------
c  initializes remnants
c-----------------------------------------------------------------------
#include "ems.h"
#include "aaa.h"

      call utpri('conre ',ish,ishini,6)

c     proj
c     ----
      la=laproj
      ma=iabs(maproj)
      las=0
      mas=0
      do l=1,ma
      if(la.lt.0)then
      ia=iabs(idproj/10)
      is=idproj/iabs(idproj)
      if(ia.ne.111.and.ia.ne.222.and.ia.ne.333)id=idproj/10*10
      if(ia.eq.111.or. ia.eq.222.or. ia.eq.333)id=idproj/10*10+is
      if(ia.eq.213)id=1230*is
      if(iabs(idproj).lt.20)id=idproj
      else
      id=1220
      if(rangen().le.(la-las)*1./(ma-mas))id=1120
      if(id.eq.1120)las=las+1
      mas=mas+1
      endif
      ic1=idtrai(1,id,1)
      ic2=idtrai(2,id,1)
      icproj(1,l)=ic1
      icproj(2,l)=ic2
      enddo

c     targ
c     ----
      la=latarg
      ma=iabs(matarg)
      las=0
      mas=0
      do l=1,ma
      if(la.lt.0)then
      ia=iabs(idtarg/10)
      is=idtarg/iabs(idtarg)
      if(ia.ne.111.and.ia.ne.222.and.ia.ne.333)id=idtarg/10*10
      if(ia.eq.111.or. ia.eq.222.or. ia.eq.333)id=idtarg/10*10+is
      if(ia.eq.213)id=1230*is
      if(iabs(idtarg).lt.20)id=idtarg
      else
      id=1220
      if(rangen().le.(la-las)*1./(ma-mas))id=1120
      if(id.eq.1120)las=las+1
      mas=mas+1
      endif
      ic1=idtrai(1,id,1)
      ic2=idtrai(2,id,1)
      ictarg(1,l)=ic1
      ictarg(2,l)=ic2
      enddo

      call utprix('conre ',ish,ishini,6)
      return
      end

c-----------------------------------------------------------------------
      subroutine conrl
c-----------------------------------------------------------------------
c  initializes target remnant in case of appl lepton
c-----------------------------------------------------------------------
#include "ems.h"
      common/nucl1/laproj,maproj,latarg,matarg,core,fctrmx
      common/hadr2/iomodl,idproj,idtarg,wexcit

c     targ
c     ----
      la=latarg
      ma=iabs(matarg)
      las=0
      mas=0
      do l=1,ma
      if(la.lt.0)then
      id=idtarg
      else
      id=1220
      if(rangen().le.(la-las)*1./(ma-mas))id=1120
      if(id.eq.1120)las=las+1
      mas=mas+1
      endif
      ic1=idtrai(1,id,1)
      ic2=idtrai(2,id,1)
      ictarg(1,l)=ic1
      ictarg(2,l)=ic2
      enddo

      return
      end

c-----------------------------------------------------------------------
      subroutine conwr
c-----------------------------------------------------------------------
c     writes /cptl/
c-----------------------------------------------------------------------
#include "aaa.h"
#include "ems.h"
      double precision XA(64,3),XB(64,3),BQGS,BMAXQGS,BMAXNEX,BMINNEX
      COMMON /Q_QGSNEX1/ XA,XB,BQGS,BMAXQGS,BMAXNEX,BMINNEX
      common/nucl3/phi,bimp
      common /cncl/xproj(mamx),yproj(mamx),zproj(mamx)
     *,xtarg(mamx),ytarg(mamx),ztarg(mamx)
      parameter(iapmax=207)
      double precision bqgs2,bmaxqgs2,bmaxnex2,bminnex2,xan,xbn
      common /qgsIInex1/xan(iapmax,3),xbn(iapmax,3)
     *,bqgs2,bmaxqgs2,bmaxnex2,bminnex2
      common/photrans/phoele(4),ebeam,noevt
      integer ic(2)

      call utpri('conwr ',ish,ishini,6)

      bx=cos(phi)*bimp
      by=sin(phi)*bimp

c     write /cptl/
c     ------------
      nptl=0

           if(iokoll.ne.0)then   ! precisely matarg collisions

      nptl=nptl+1
      do 3 i=1,4
3     xorptl(i,nptl)=0
      tivptl(1,nptl)=-ainfin
      tivptl(2,nptl)=0
      istptl(nptl)=1
      iorptl(nptl)=-1
      jorptl(nptl)=0
      do 1 k=1,koll
      nptl=nptl+1
      do 4 i=1,4
4     xorptl(i,nptl)=0
      tivptl(1,nptl)=-ainfin
      tivptl(2,nptl)=0
      istptl(nptl)=1
      iorptl(nptl)=-1
      jorptl(nptl)=0
1     continue

           elseif(iappl.ne.7)then

c             print *,'proj'

      do 6 i=1,maproj
      nptl=nptl+1
      istptl(nptl)=0
      iorptl(nptl)=0
      jorptl(nptl)=0
      if(model.eq.2)then       !QGSJet
      xproj(i)=XA(i,1)
      yproj(i)=XA(i,2)
      zproj(i)=XA(i,3)
      istptl(nptl)=1
      iorptl(nptl)=-1
      elseif(model.eq.7)then       !QGSJetII
      xproj(i)=xan(i,1)
      yproj(i)=xan(i,2)
      zproj(i)=xan(i,3)
      istptl(nptl)=1
      iorptl(nptl)=-1
      elseif(model.ge.3)then       !Gheisha, ...
      istptl(nptl)=1
      iorptl(nptl)=-1
      endif
      xorptl(1,nptl)=xproj(i)+bx/2
      xorptl(2,nptl)=yproj(i)+by/2
      xorptl(3,nptl)=zproj(i)
      xorptl(4,nptl)=0
      tivptl(1,nptl)=-ainfin
c for visualisation uncomment
c-c   tivptl(1,nptl)=-100
      tivptl(2,nptl)= ainfin
c             print *,i,xorptl(1,nptl),xorptl(2,nptl),xorptl(3,nptl)
6     continue
c             print *,'targ'
      do 7 i=1,matarg
      nptl=nptl+1
      istptl(nptl)=0
      iorptl(nptl)=0
      jorptl(nptl)=0
      if(model.eq.2)then       !QGSJet
      xtarg(i)=XB(i,1)
      ytarg(i)=XB(i,2)
      ztarg(i)=XB(i,3)
      istptl(nptl)=1
      iorptl(nptl)=-1
      elseif(model.eq.7)then       !QGSJetII
      xtarg(i)=xbn(i,1)
      ytarg(i)=xbn(i,2)
      ztarg(i)=xbn(i,3)
      istptl(nptl)=1
      iorptl(nptl)=-1
      elseif(model.ge.3)then       !Gheisha, ...
      istptl(nptl)=1
      iorptl(nptl)=-1
      endif
      xorptl(1,nptl)=xtarg(i)-bx/2
      xorptl(2,nptl)=ytarg(i)-by/2
      xorptl(3,nptl)=ztarg(i)
      xorptl(4,nptl)=0
      tivptl(1,nptl)=-ainfin
c for visualisation uncomment
c-c   tivptl(1,nptl)=-100
      tivptl(2,nptl)= ainfin
c             print *,i,xorptl(1,nptl),xorptl(2,nptl),xorptl(3,nptl)
7     continue
      if(abs(idprojin).eq.12)then   !electron for fake DIS
c electron projectile
        nptl=nptl+1
        istptl(nptl)=41
        iorptl(nptl)=-1
        jorptl(nptl)=-1
        iorptl(1)=nptl         !pi0 (porjectile) coming from lepton
        xorptl(1,nptl)=bx/2
        xorptl(2,nptl)=by/2
        xorptl(3,nptl)=0.
        xorptl(4,nptl)=0.
        tivptl(1,nptl)=-ainfin
        tivptl(2,nptl)=0.
c target nucleons (in lab frame)
        do i=1,matarg
          nptl=nptl+1
          istptl(nptl)=41
          iorptl(nptl)=-1
          jorptl(nptl)=-1
          xorptl(1,nptl)=xtarg(i)-bx/2
          xorptl(2,nptl)=ytarg(i)-by/2
          xorptl(3,nptl)=ztarg(i)
          xorptl(4,nptl)=0
          tivptl(1,nptl)=-ainfin
c         for visualisation uncomment
c         -c   tivptl(1,nptl)=-100
          tivptl(2,nptl)= ainfin
        enddo
c electron remnant
        nptl=nptl+1
        istptl(nptl)=0
        iorptl(nptl)=maproj+matarg+1
        jorptl(nptl)=-1
        xorptl(1,nptl)=bx/2
        xorptl(2,nptl)=by/2
        xorptl(3,nptl)=0.
        xorptl(4,nptl)=0.
        tivptl(1,nptl)=0.
        tivptl(2,nptl)= ainfin
      endif

          endif

      nptl=0
      if(iappl.le.2)then
      do i=1,maproj
      nptl=nptl+1
      ic(1)=icproj(1,i)
      ic(2)=icproj(2,i)
      id=idtra(ic,0,0,0)
c      id=idtra(ic,0,0,3)      !tp071107 imix=3 ??????????
      call idmass(id,ams)
      idptl(nptl)=id
      pptl(1,nptl)=0.
      pptl(2,nptl)=0.
      pptl(3,nptl)=pnullx
      pptl(4,nptl)=sqrt(pnullx**2+ams**2)
      pptl(5,nptl)=ams
      ifrptl(1,nptl)=0
      ifrptl(2,nptl)=0
      ityptl(nptl)=0
      zpaptl(1,nptl)=0.
      zpaptl(2,nptl)=0.
      enddo
      endif
      if(iappl.ne.7)then
      do i=1,matarg
      nptl=nptl+1
      ic(1)=ictarg(1,i)
      ic(2)=ictarg(2,i)
      id=idtra(ic,0,0,0)
c      id=idtra(ic,0,0,3)      !tp071107 imix=3 ??????????
      call idmass(id,ams)
      idptl(nptl)=id
      pptl(1,nptl)=0.
      pptl(2,nptl)=0.
      pptl(3,nptl)=-pnullx
      pptl(4,nptl)=sqrt(pnullx**2+ams**2)
      pptl(5,nptl)=ams
      ifrptl(1,nptl)=0
      ifrptl(2,nptl)=0
      ityptl(nptl)=0
      zpaptl(1,nptl)=0.
      zpaptl(2,nptl)=0.
      enddo
      if(abs(idprojin).eq.12)then   !electron for fake DIS
c electron projectile
        nptl=nptl+1
        id=idprojin
        call idmass(id,ams)
        idptl(nptl)=id
        pptl(1,nptl)=0.
        pptl(2,nptl)=0.
        pptl(3,nptl)=sqrt(max(0.,(elepti+ams)*(elepti-ams)))
        pptl(4,nptl)=elepti
        pptl(5,nptl)=ams
        ifrptl(1,nptl)=1
        ifrptl(2,nptl)=1
        ityptl(nptl)=40
        zpaptl(1,nptl)=0.
        zpaptl(2,nptl)=0.
c target nucleons (in lab frame)
        do i=1,matarg
          nptl=nptl+1
          idptl(nptl)=idptl(maproj+i)
          pptl(1,nptl)=0.
          pptl(2,nptl)=0.
          pptl(3,nptl)=-pnll
          pptl(4,nptl)=ebeam
          pptl(5,nptl)=pptl(5,maproj+i)
          ifrptl(1,nptl)=maproj+i
          ifrptl(2,nptl)=maproj+i
          ityptl(nptl)=50
          zpaptl(1,nptl)=0.
          zpaptl(2,nptl)=0.
        enddo
c electron remnant
        nptl=nptl+1
        idptl(nptl)=id
        pptl(1,nptl)=phoele(1)
        pptl(2,nptl)=phoele(2)
        pptl(3,nptl)=phoele(3)
        pptl(4,nptl)=phoele(4)
        pptl(5,nptl)=ams
        ifrptl(1,nptl)=0
        ifrptl(2,nptl)=0
        ityptl(nptl)=40
        zpaptl(1,nptl)=0.
        zpaptl(2,nptl)=0.
      endif

      else

      nptl=nptl+1
      id=idproj
      call idmass(id,ams)
      idptl(nptl)=id
      pptl(1,nptl)=0.
      pptl(2,nptl)=0.
      pptl(3,nptl)=pnullx
      pptl(4,nptl)=sqrt(pnullx**2+ams**2)
      pptl(5,nptl)=ams
      ifrptl(1,nptl)=0
      ifrptl(2,nptl)=0
      ityptl(nptl)=0
      iorptl(nptl)=-1
      jorptl(nptl)=0
      istptl(nptl)=0
      do 5 i=1,4
 5      xorptl(i,nptl)=0
      tivptl(1,nptl)=0
      tivptl(2,nptl)=0
      zpaptl(1,nptl)=0.
      zpaptl(2,nptl)=0.
      endif

c     exit
c     ----

      call utprix('conwr ',ish,ishini,6)
      return
      end

c------------------------------------------------------------------------
      subroutine conxyz(ch,n,x,y,z,ynuc)
c-----------------------------------------------------------------------
#include "aaa.h"

      real x(n),y(n),z(n)
      character ch*1

      if(ch.eq.'p')then
      massnr=maproj
      iii=1
      elseif(ch.eq.'t')then
      massnr=matarg
      iii=2
      else
      massnr=0
      iii=0
      call utstop('conxyz: nucleus neither proj nor targ&')
      endif

      if(massnr.eq.0)return
      if(massnr.gt.n)call utstop('conxyz: massnr.gt.n&')
      if(massnr.eq.1)then
      x(1)=0
      y(1)=0
      z(1)=0
      return
      endif

      rad=radnuc(massnr)
      !write(ifmt,'(a,2f7.3)')'rad,dif:',rad,difnuc(massnr)

      if(massnr.ge.10)then !---wood-saxon density---

        rad=rad/difnuc(massnr)
        cr1=1.+3./rad+6./rad**2+6./rad**3
        cr2=3./rad
        cr3=3./rad+6./rad**2
        do i=1,massnr
   1      zuk=rangen()*cr1-1.
          if(zuk.le.0.)then
            tt=rad*(rangen()**.3333-1.)
          elseif(zuk.le.cr2 )then
            tt=-log(rangen())
          elseif(zuk.lt.cr3 )then
            tt=-log(rangen())-log(rangen())
          else
            tt=-log(rangen())-log(rangen())-log(rangen())
          endif
          if(rangen().gt.1./(1.+exp(-abs(tt))))goto 1
          rim=tt+rad
          zz=rim*(2.*rangen()-1.)
          rim=sqrt(rim*rim-zz*zz)
          z(i)=zz*difnuc(massnr)
          call pscs(c,s)
          x(i)=rim*c*difnuc(massnr)
          y(i)=rim*s*difnuc(massnr)
          if(hacore.gt.0.)then
          do j=1,i-1
           di2= (x(i)-x(j))**2 + (y(i)-y(j))**2+ (z(i)-z(j))**2
           if(di2.lt.hacore**2)goto1
          enddo
          endif
c          rrr=radnuc(massnr)    !mj     hard sphere nuclear coordinates
c          zzz=(2.*rangen()-1.)   !mj
c          phi=2*3.14159*rangen()   !mj
c          x(i)=rrr*sqrt(1-zzz**2)*cos(phi)  !mj
c          y(i)=rrr*sqrt(1-zzz**2)*sin(phi)   !mj
c          z(i)=rrr*zzz     !mj
        enddo

      elseif(massnr.ge.3)then  ! ---gaussian density---

        rad=rad*sqrt(2.*massnr/(massnr-1.))   !van hove simulation
        do l=1,3
          summ=0.
          do i=1,massnr-1
            j=massnr-i
            aks=rad *(rangen()+rangen()+rangen()-1.5)
            k=j+1
            if(l.eq.1)x(k)=summ-aks*sqrt(float(j)/k)
            if(l.eq.2)y(k)=summ-aks*sqrt(float(j)/k)
            if(l.eq.3)z(k)=summ-aks*sqrt(float(j)/k)
            summ=summ+aks/sqrt(float(j*k))
          enddo
          if(l.eq.1)x(1)=summ
          if(l.eq.2)y(1)=summ
          if(l.eq.3)z(1)=summ
        enddo

      elseif(massnr.eq.2)then  ! ---deuteron---

        !.........r=t*difnuc(massnr), t~exp(-2*t)*(1-exp(-a*t))
        a=radnuc(massnr)
  2     t=-0.5*alog(rangen())  !~exp(-2*t)
        if(rangen().gt.(1-exp(-a*t))**2)goto2
        r=t*difnuc(massnr)
        zz=r*(2.*rangen()-1.)
        call pscs(c,s)
        rxy=sqrt(r*r-zz*zz)
        z(1)=0.5*zz
        x(1)=0.5*rxy*c
        y(1)=0.5*rxy*s
        z(2)=-z(1)
        x(2)=-x(1)
        y(2)=-y(1)

      else

        stop'conxyz: wrong massnr.     '

      endif

c...plot preparation

      rmax=(radnuc(massnr)+3)
      drnucl(iii)=rmax/mxnucl
      nrnucl(iii)=nrnucl(iii)+1
      do i=1,massnr
        r=sqrt(x(i)**2+y(i)**2+z(i)**2)
        k=1+int(r/drnucl(iii))
        if(k.le.mxnucl)rnucl(k,iii)=rnucl(k,iii)+1
      enddo

c...lorentz trafo

      fac=1./cosh(ynuc) 
      !fac=  0  !!!!???????????????????????  
      do i=1,massnr
      z(i)=  z(i)*fac
      enddo

      return
      end

c-----------------------------------------------------------------------
      subroutine conini
c-----------------------------------------------------------------------
#include "aaa.h"

      imax=max(maproj,matarg)
      if(idtargin.eq.0)imax=max(imax,40)
      do massnr=1,mamxx
        dif=0.54
        if(massnr.gt.imax.or.massnr.eq.1)then
          dif=0
          rad=0
        elseif(massnr.eq.197)then
          dif=0.562
          rad=6.5
        elseif(massnr.ge.10)then
          rad=1.12*massnr**0.33333-0.86*massnr**(-0.33333)
        elseif(massnr.ge.3)then
          rad=.9*float(massnr)**.3333
        elseif(massnr.eq.2)then
          dif=4.316
          rad=4.68
        endif
        difnuc(massnr)=dif
        radnuc(massnr)=rad
      enddo
      end

c-----------------------------------------------------------------------
      subroutine xConNuclDens(iii)
c-----------------------------------------------------------------------
c plots distribution of nucleons in nuclei
c  iii = 1 (proj) or 2 (targ)
c-----------------------------------------------------------------------
#include "aaa.h"
      if(model.ne.1)return
      if(iii.eq.1)then
        massnr=maproj
      else!if(iii.eq.2)then
        massnr=matarg
      endif
      if(massnr.eq.1)return
      a=1./4.316
      b=4.68
      write(ifhi,'(a)') '!-----------------------------------------'
      write(ifhi,'(a)') '!          nuclear density          '
      write(ifhi,'(a)') '!-----------------------------------------'
      write(ifhi,'(a)')       'openhisto'
      if(massnr.ge.10)write(ifhi,'(a)')'htyp lin xmod lin ymod lin'
      if(massnr.lt.10)write(ifhi,'(a)')'htyp lin xmod lin ymod log'
      write(ifhi,'(a)')       'text 0 0 "title nuclear density"'
      write(ifhi,'(a)')       'text 0.99 -0.15 " r (fm) "'
      write(ifhi,'(a)')       'text 0 0 "yaxis rho(r)"'
      write(ifhi,'(a,2e11.3)')'xrange',0.,mxnucl*drnucl(iii)
      write(ifhi,'(3a)')'yrange',' 0 ',' auto'
      write(ifhi,'(a)')       'array 2'
      do j=1,mxnucl
      r=(j-0.5)*drnucl(iii)
      d=0.5*drnucl(iii)
      write(ifhi,'(2e12.4)')  r,rnucl(j,iii)/nrnucl(iii)/
     *                     (4./3.*pi*((r+d)**3-(r-d)**3))
      enddo
      write(ifhi,'(a)')       '  endarray'
      write(ifhi,'(a)')       'closehisto plot 0-'
      write(ifhi,'(a)')       'openhisto'
      write(ifhi,'(a)')       'htyp lbo '
      write(ifhi,'(a)')       'array 2'
      do j=1,mxnucl
      r=(j-0.5)*drnucl(iii)
      rr=2*r
      rho=1
      if(massnr.eq.2)then
        rho=1.00*((1-exp(-b*a*rr))*exp(-a*rr)/rr)**2
      elseif(massnr.eq.197)then
        rho=0.16/(1+exp((r-6.5)/0.562))
      elseif(massnr.ge.10)then
        rho=0.16/(1+exp((r-radnuc(massnr))/difnuc(massnr)))
      endif
      write(ifhi,'(2e12.4)')  r,rho
      enddo
      write(ifhi,'(a)')       '  endarray'
      write(ifhi,'(a)')       'closehisto plot 0'
      end

c-----------------------------------------------------------------------
      subroutine xConThick(iii)
c-----------------------------------------------------------------------
      ! plots sigma_pp *T_A (b)  (=average number of collisions)
      ! T_A = thickness function
      !  iii = 1 (proj) or 2 (targ)
      !----------------------------------------------------------------
#include "aaa.h"
      parameter(iconimax=20,iconkmax=100)
      real thick(2,0:iconimax)
      imax=iconimax
      kmax=iconkmax
      if(model.ne.1)return
      ramx=mxnucl*drnucl(iii)
      do i=0,imax
        x=i/float(imax)*ramx
        sum=0
        rho0=conrho(iii,0.)
        h=ramx/kmax
        do k=1,kmax
          z=k*h
          r=sqrt(x*x+z*z)
          rho2=conrho(iii,r)
          z=(k-0.5)*h
          r=sqrt(x*x+z*z)
          rho1=conrho(iii,r)
          sum=sum+h/6.*(rho0+4*rho1+rho2)
          rho0=rho2
        enddo
        sum=sum*2 ! integral fro -infty to +infty
        thick(1,i)=x
        thick(2,i)=sum
      enddo
      write(ifhi,'(a)')       'openhisto'
      write(ifhi,'(a)')       'htyp lbo '
      write(ifhi,'(a)')       'txt "xaxis b (fm)" '
      write(ifhi,'(a)')       'txt "yaxis [s]?pp! T?A! (b) " '
      write(ifhi,'(a)')       'array 2'
      do i=0,imax
        write(ifhi,'(2e12.4)') thick(1,i),sigine/10.*thick(2,i)
      enddo
      write(ifhi,'(a)')       '  endarray'
      write(ifhi,'(a)')       'closehisto plot 0'

      end

c-----------------------------------------------------------------------
      function conrho(iii,r)
c-----------------------------------------------------------------------
      ! nuclear density
      !  iii = 1 (proj) or 2 (targ)
      !----------------------------------------------------------------
#include "aaa.h"
      conrho=1.
      if(model.ne.1)return
      if(iii.eq.1)then
        massnr=maproj
      else!if(iii.eq.2)then
        massnr=matarg
      endif
      if(massnr.eq.1)return
      a=1./4.316
      b=4.68
      rr=2*r
      rho=1
      if(massnr.eq.2.and.rr.gt.0.)then
        rho=1.00*((1-exp(-b*a*rr))*exp(-a*rr)/rr)**2
      elseif(massnr.eq.197)then
        rho=0.16/(1+exp((r-6.5)/0.562))
      elseif(massnr.ge.10)then
        rho=0.16/(1+exp((r-radnuc(massnr))/difnuc(massnr)))
      endif
      conrho=rho
      end




c----------------------------------------------------------------------
      subroutine geoglauber
c----------------------------------------------------------------------
#include "aaa.h"
#include "ems.h"
      common/geom/rmproj,rmtarg,bmax,bkmx
      common/nucl3/phi,bimp
      common /cncl/xproj(mamx),yproj(mamx),zproj(mamx)
     *            ,xtarg(mamx),ytarg(mamx),ztarg(mamx)
      common/ckopjtg/kopj(mamx),kotg(mamx)
      common/cbglaub/bglaub
      logical phi_zero

      phi_zero=.false.


! calculate collision number according to Glauber ------------------

c      rs=r2had(iclpro)+r2had(icltar)+slopom*log(sy)
c      bglaub2=4.*.0389*rs
      xs=sigine   !xsectionpar()   !?????????????????????????????????
      bglaub=max(0.,xs/10./pi)        !10= fm^2 -> mb
c      bglhard=bglaub*min(1.,xfrachard)
      bglaub=sqrt(bglaub)
      nglevt=0
      nglndif=0
      nglhard=0
      bglndif=max(0.,sigcut/10./pi) !10= fm^2 -> mb
      bglhard=bglndif*xfrachard
      bglndif=sqrt(bglndif)
      bglhard=sqrt(bglhard)
      do i=1,maproj
        kopj(i)=0
      enddo
      do i=1,matarg
        kotg(i)=0
      enddo
      if(abs(phimin).lt.1e-5.and.abs(phimax).lt.1e-5)phi_zero=.true.
      if(phi_zero)then
        x2av=0
        y2av=0
        xyav=0
        xav=0
        yav=0
        bx=cos(phi)*bimp
        by=sin(phi)*bimp
      else     !not used
        bx=0.
        by=0.
      endif  
      do ko=1,koll
        ip=iproj(ko)
        it=itarg(ko)
        r=bk(ko)
        if(phi_zero)then
          r2=(xproj(ip)+bx-xtarg(it))**2+(yproj(ip)+by-ytarg(it))**2
          if((maproj.gt.1.or.matarg.gt.1).and.
     .       abs(r**2-r2).gt.1e-3)then
            print*,'geoglauber: r**2,r2:',r**2,r2
            stop
          endif
        endif
        if(r.le.bglaub)then
          nglevt=nglevt+1
          kopj(ip)=kopj(ip)+1
          kotg(it)=kotg(it)+1
c pair dependent hard interactions
c          bglhard=max(0.,xfrachard*sngl(sigmak(1,ko))/10./pi) !10= fm^2 -> mb
c          bglhard=sqrt(bglhard)
          if(r.le.bglhard)nglhard=nglhard+1
        endif
        if(r.le.bglndif)nglndif=nglndif+1
      enddo
      ng1evt=0
      ng11evt=0
      ng12evt=0
      ng2evt=0
      ngspecp=laproj+latarg
      ngspecn=maproj+matarg-ngspecp
      do i=1,maproj
        if(kopj(i).ge.1)then
          ng1evt=ng1evt+1
          ng11evt=ng11evt+1
          if(idptl(i).eq.1120)ngspecp=ngspecp-1
          if(idptl(i).eq.1220)ngspecn=ngspecn-1
        endif 
        if(kopj(i).ge.2)ng2evt=ng2evt+1
      enddo   
      do i=1,matarg
        if(kotg(i).ge.1)then
          ng1evt=ng1evt+1
          ng12evt=ng12evt+1
          if(idptl(maproj+i).eq.1120)ngspecp=ngspecp-1
          if(idptl(maproj+i).eq.1220)ngspecn=ngspecn-1
        endif
        if(kotg(i).ge.2)ng2evt=ng2evt+1
      enddo   
      eglevt=0
      fglevt=0
      rglevt=0
      sglevt=0
      if(phi_zero)then
        do i=1,maproj
          if(kopj(i).ge.1)then
            x=xproj(i)+bx/2
            y=yproj(i)+by/2
            xav=xav+x
            yav=yav+y
            x2av=x2av+x*x
            y2av=y2av+y*y
            xyav=xyav+x*y
            !print*,x,y
          endif
        enddo
        do i=1,matarg
          if(kotg(i).ge.1)then
            x=xtarg(i)-bx/2
            y=ytarg(i)-by/2
            xav=xav+x
            yav=yav+y
            x2av=x2av+x*x
            y2av=y2av+y*y
            xyav=xyav+x*y
            !print*,x,y
          endif
        enddo
        if(ng1evt.gt.0)then
        x2av=x2av/float(ng1evt)
        y2av=y2av/float(ng1evt)
        xyav=xyav/float(ng1evt)
        xav= xav /float(ng1evt)
        yav= yav /float(ng1evt)
        eglevt=(y2av-x2av)/(y2av+x2av)
        fglevt=sqrt( (y2av-x2av)**2+4*(xyav-xav*yav)**2 )  /(y2av+x2av)
        rglevt=ng2evt/float(ng1evt)
        sglevt=pi*sqrt(x2av)*sqrt(y2av)
        !~~~~~~~~~~~~~~~~~~~~
        xx=x2av ! <x**2>
        yy=y2av ! <y**2>
        xy=xyav ! <x*y>
        dta=0.5*abs(xx-yy)
        ev1=(xx+yy)/2+sqrt(dta**2+xy**2)
        ev2=(xx+yy)/2-sqrt(dta**2+xy**2)
        yy=ev1
        xx=ev2
        fglevt=(yy-xx)/(yy+xx)
        !~~~~~~~~~~~~~~~~~~~~
        endif
        !if(ng1evt.gt.0)print*,bimp,eglevt,bx,by,phi
      endif  
      !print*,'TESTgeoglauber',bimp,bglaub,ng1evt
  
      end

c-----------------------------------------------------------------------
      function percentile()
c-----------------------------------------------------------------------
#include "aaa.h"  
      common/ccc20/icc20 
      common/nucl3/phi,bimp
      real wi(2)
      if(izmode.ne.1)stop'\n\n  ERROR 24022013c \n\n'
      if(icc20.ne.1)stop'\n\n  ERROR 24022013b \n\n'
      do i=1,20
        ii=i
        if(bimp.ge.zclass(1,i).and.bimp.le.zclass(2,i))goto 1
      enddo
  1   cmin=(ii-1)*0.05
      cmax= ii*0.05
      bmin=zclass(1,ii)
      bmax=zclass(2,ii)
      frac=(bimp-bmin)/(bmax-bmin)
      wi(1)=1-frac
      wi(2)=frac
      percentile   = wi(1)*cmin +wi(2)*cmax
      end
      
      



