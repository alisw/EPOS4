C
C  This file is part of EPOS4
C  Copyright (C) 2022 research institutions and authors (See CREDITS file)
C  This file is distributed under the terms of the GNU General Public License version 3 or later
C  (See COPYING file for the text of the licence)
C

c-----------------------------------------------------------------------
c  inicon table means storage of IcoE, IcoV, IcoF
c     IcoE ... energy density in the comoving frame
c     IcoV ... flow velocity in hyperbolic coordinates
c               v=(vx,vy,veta=vz/tau)
c                with the velocity in the frame moving with rapidity eta
c     IcoV ... net flavor density in the comoving frame
c  the indices of these fields ix,iy,iz also refer to hyperbolic coordinates
c     (x,y,eta) at given tau
c     the corresponding grid is defined in ico.h
c-----------------------------------------------------------------------

c-----------------------------------------------------------------------
      subroutine IniCon(nin,nev,ierrico)
c-----------------------------------------------------------------------
#include "aaa.h"
#include "ico.h"
#include "ems.h"
      double precision seedf
      character*80 fn
      logical table
      common/cranphi/ranphi
      integer ioicoplot
      common /cioicoplot/ioicoplot
      data ncntmico/0/
      save ncntmico
      inicop=0
      if(iorsdf.eq.3)inicop=1
      if(iorsdf.eq.5)inicop=2
      ierrico=0

      if(icotabr.eq.1)then  !read from file

        inquire(file=fnio(1:nfnio),exist=table)
        if(.not.table)then
          write(ifmt,'(//10x,3a//)')'STOP in IniCon: file "'
     *             ,fnio(1:nfnio),'" not found !!!'
          stop
        endif
        write(ifmt,'(3a)')'read from ',fnio(1:nfnio),' ...'
        open(97,file=fnio(1:nfnio),status='old')
        read(97,*) iversn
        if(iversn.lt.201)stop'\n\n redo table with version >= 201\n\n'
        read(97,*) laprojx,maprojx,latargx,matargx
        read(97,*) engyx
        read(97,*) bminimx,bmaximx,ikolmnx,ikolmxx
        bimevt=(bminimx+bmaximx)/2
        nprt(1)=(ikolmnx+ikolmxx)/2
        phievt=0
        ranphi=0
        read(97,*) tauzerx
        read(97,*) iabs_ninicon
        read(97,*) nxicox,nyicox,nzicox
        read(97,*)xminicox,xmaxicox,yminicox,ymaxicox,zminicox,zmaxicox
        if(laprojx.ne.laproj.or.maprojx.ne.maproj
     . .or.latargx.ne.latarg.or.matargx.ne.matarg)
     .  stop'IniCon: different nuclei\n\n'
        if(engyx.ne.engy)
     .  stop'IniCon: different energy\n\n  '
        if(abs(bmaximx-bmaxim).gt.0.001)
     .  stop'IniCon: different max impact parameter\n\n'
        if(abs(bminimx-bminim).gt.0.001)
     .  stop'IniCon: different min impact parameter\n\n'
        if(abs(ikolmxx-ikolmx).gt.0.001)
     .  stop'IniCon: different max kol\n\n'
        if(abs(ikolmnx-ikolmn).gt.0.001)
     .  stop'IniCon: different min kol\n\n'
        if(abs(tauzerx-tauzer).gt.0.001)
     .  stop'IniCon: different tauzer\n\n'
        if(nxicox.ne.nxico.or.nyicox.ne.nyico
     . .or.nzicox.ne.nzico)
     .  stop'IniCon:  different nxico/nyico/nzico\n\n'
        if( xminicox.ne.xminico.or.xmaxicox.ne.xmaxico
     . .or.yminicox.ne.yminico.or.ymaxicox.ne.ymaxico
     . .or.zminicox.ne.zminico.or.zmaxicox.ne.zmaxico)
     .  stop'IniCon: different x/y/zminico or x/y/zmaxico\n\n'
        read(97,*)
     .   (((IcoE(ix,iy,iz)         ,ix=1,nxico),iy=1,nyico),iz=1,nzico)
     . ,((((IcoV(i,ix,iy,iz),i=1,3),ix=1,nxico),iy=1,nyico),iz=1,nzico)
     . ,((((IcoF(i,ix,iy,iz),i=1,3),ix=1,nxico),iy=1,nyico),iz=1,nzico)
        close(97)
        if(nsmooth.gt.0)call SmoothIco
        call amaxIcoE(wmax,x,y,z)
        write(ifmt,'(a,f9.2,a,3f8.2)')'wmax',wmax,' x,y,eta  ',x,y,z

      elseif(icotabr.eq.0)then  ! calculate

        if(nin.eq.1)then
          write(ifmt,'(2a,i7,a)')'calculate inicon table ',
     &    'based on',iabs(ninicon),'  ico-events ...'
          if(inicop.eq.1)then
          stop'\n\n IcoStrSeg no more  supported in 212 and later \n\n'
          else
          call IcoStr(1,ierrico)
          endif
        endif

        call ranfgt(seedf)
        if((nin.le.iabs(ninicon).and.mod(nin,modsho).eq.0)
     &     .or.iabs(ninicon).eq.1   )
     &        write(ifmt,'(a,i7,a,f6.2,a,d27.16)')'ico-event ',nin
     &                   ,'   bimevt:',bimevt,'   seedf:',seedf
        if(nin.le.iabs(ninicon).and.jwseed.eq.1)then
          open(unit=1,file=fnch(1:nfnch-5)//'see',status='unknown')
          write(1,'(a,i7,a,f6.2,a,d27.16)')'ico-event ',nin
     &                   ,'   bimevt:',bimevt,'   seedf:',seedf
         close(1)
        endif
        !if(nin.le.iabs(ninicon).and.mod(nin,modsho).eq.0)
        !&         write(ifmt,*)'+++++ time: ',timefin-timeini
        seedc=seedf
        if(inicop.eq.1)then
         stop'\n\n IcoStrSeg no more  supported in 212 and later \n\n'
        else
         call IcoStr(2,ierrico)
         if(ierrico.eq.1)goto 1001
        endif

        if(nin.eq.iabs(ninicon))then
          iabs_ninicon=iabs(ninicon)
          write(ifmt,'(a)')'normalize engy-mom tensor'
          if(inicop.eq.1)then
           stop'\n\n IcoStrSeg no more  supported in 212 and later \n\n'
          else
           call IcoStr(3,ierrico)
          endif
          if(nsmooth.gt.0)call SmoothIco
          ncntmico=ncntmico+1
          !!!if(ncntmico.ne.8)goto99 ! for debugging
          if(icotabm.gt.0)then  ! write table
            ii=index(fnhi(1:nfnhi),".histo")-1
            if(jcentrality.ge.0)then
              fn=fnhi(1:ii)//".ico"
              iixx=ii+4
            else
              iix=ii
              do while(fnhi(ii:ii).ne.'-')
              ii=ii-1
              enddo
              fn(1:ii)=fnhi(1:ii)
              imax=0
              do i=1,100
                if(zclass(3,i).gt.1e-5)imax=i
              enddo
              icentrality=1+mod((iopcnt-1)/(-jcentrality),imax)
              if(icentrality.le.9)then
                write(fn(ii+1:ii+1),'(i1)')icentrality
                fn(ii+2:ii+2)='-'
                fn=fn(1:ii+2)//fnhi(ii+1:iix)//".ico"
                iixx=iix+6
              elseif(icentrality.le.99)then
                write(fn(ii+1:ii+2),'(i2)')icentrality
                fn(ii+3:ii+3)='-'
                fn=fn(1:ii+3)//fnhi(ii+1:iix)//".ico"
                iixx=iix+7
              else
                stop'in IniCon: icentrality.gt.99\n\n'
              endif
            endif
            write(ifmt,'(2a)')'write inicon table into ',fn(1:iixx)
            open(97,file=fn(1:iixx),status='unknown')
            write(97,*) iversn
            write(97,*) laproj,maproj,latarg,matarg
            write(97,*) engy
            write(97,*) bminim,bmaxim,ikolmn,ikolmx
            write(97,*) tauzer
            write(97,*) iabs_ninicon
            write(97,*) nxico,nyico,nzico
            write(97,*) xminico,xmaxico,yminico,ymaxico,zminico,zmaxico
            write(97,*)
     .   (((IcoE(ix,iy,iz)         ,ix=1,nxico),iy=1,nyico),iz=1,nzico)
     . ,((((IcoV(i,ix,iy,iz),i=1,3),ix=1,nxico),iy=1,nyico),iz=1,nzico)
     . ,((((IcoF(i,ix,iy,iz),i=1,3),ix=1,nxico),iy=1,nyico),iz=1,nzico)
            !if(ispherio.eq.0)call xIco2
            close(97)
          endif
          icotabt=0
          if(icotabt.eq.1)then
            ii=index(fnhi(1:nfnhi),".histo")-1
            fn=fnhi(1:ii)//".tmunu"
            iixx=ii+6
            open(97,file=fn(1:iixx),status='unknown')
            write(97,*) nxico,nyico,nzico
            write(97,*) xminico,xmaxico,yminico,ymaxico,zminico,zmaxico
            write(97,*)
     .     (((((IcoT(i1,i2,ix,iy,iz),i1=1,4),i2=1,4),
     ,      ix=1,nxico),iy=1,nyico),iz=1,nzico)
            close(97)
          endif
          !!  99         continue
          call xIco(nev)
          call xIco2
          if(ioicoplot.eq.1)call xIcoPlot(nev)
          call amaxIcoE(wmax,x,y,z)
          write(ifmt,'(a,f9.2,a,3f8.2)')'wmax',wmax,' x,y,eta  ',x,y,z
        endif

      else

        stop'IniCon: wrong choice'

      endif

 1001 continue
      call clop(3)
      end

c-----------------------------------------------------------------------
      subroutine amaxIcoE(wmax,x,y,z)
c-----------------------------------------------------------------------
#include "aaa.h"
#include "ico.h"
      common/cwmaxIcoE/wmaxIcoE
      wmax=0.
      x=0.
      y=0.
      z=0.
      i=0
      j=0
      k=0
      do ix=1,nxico
      do iy=1,nyico
      do iz=1,nzico
      w=IcoE(ix,iy,iz)
      if(w.gt.wmax)then
        wmax=w
        i=ix
        j=iy
        k=iz
      endif
      enddo
      enddo
      enddo
      wmaxIcoE=wmax
      if(wmax.gt.0.)then
        x=xminico+(i-0.5)*(xmaxico-xminico)/nxico
        y=yminico+(j-0.5)*(ymaxico-yminico)/nyico
        z=zminico+(k-0.5)*(zmaxico-zminico)/nzico
      endif

      end

c-----------------------------------------------------------------------
      subroutine SmoothIco
c-----------------------------------------------------------------------
#include "aaa.h"
#include "ico.h"
      double precision  IcoTx, IcoCx
      real IcoEx, IcoVx, IcoFx
      dimension
     . IcoTx(4,4,nxicomax,nyicomax,nzicomax)
     .,IcoCx(3,4,nxicomax,nyicomax,nzicomax)
     .,IcoEx(nxicomax,nyicomax,nzicomax)
     .,IcoVx(3,nxicomax,nyicomax,nzicomax)
     .,IcoFx(3,nxicomax,nyicomax,nzicomax)
      double precision w(-5:5),w4(-2:2),w6(-3:3)
      data w4 /1.d0, 4.d0, 6.d0, 4.d0, 1.d0/
      data w6 /1.d0, 6.d0, 15.d0, 20.d0, 15.d0, 6.d0, 1.d0/

      write(ifmt,*)'smoothing initial conditions, nsmooth =',nsmooth

      nn=nsmooth/2
      do n=-nn,nn
        if(nsmooth.eq.4)then
          w(n)=w4(n)/16d0
        elseif(nsmooth.eq.6)then
          w(n)=w6(n)/64d0
        else
          stop'in SmoothIco \n\n'
        endif
      enddo

      do ix=1,nxico
        do iy=1,nyico
          do iz=1,nzico
            IcoEx(ix,iy,iz)=IcoE(ix,iy,iz)
            do i=1,3
            IcoVx(i,ix,iy,iz)=IcoV(i,ix,iy,iz)
            IcoFx(i,ix,iy,iz)=IcoF(i,ix,iy,iz)
            enddo
            do j=1,4
            do i=1,4
            IcoTx(i,j,ix,iy,iz)=IcoT(i,j,ix,iy,iz)
            enddo
            do i=1,3
            IcoCx(i,j,ix,iy,iz)=IcoC(i,j,ix,iy,iz)
            enddo
            enddo
          enddo
        enddo
      enddo

      do ix=1,nxico
        do iy=1,nyico
          do iz=1,nzico
            IcoE(ix,iy,iz)=0.
            do i=1,3
            IcoV(i,ix,iy,iz)=0.
            IcoF(i,ix,iy,iz)=0.
            enddo
            do j=1,4
            do i=1,4
            IcoT(i,j,ix,iy,iz)=0.
            enddo
            do i=1,3
            IcoC(i,j,ix,iy,iz)=0.
            enddo
            enddo
          enddo
        enddo
      enddo

      do ix=1,nxico
      do iy=1,nyico
        do nx=max(1,ix-nn),min(nxico,ix+nn)
        do ny=max(1,iy-nn),min(nyico,iy+nn)
          do iz=1,nzico
            IcoE(ix,iy,iz)=IcoE(ix,iy,iz)
     .     +IcoEx(nx,ny,iz)*w(ix-nx)*w(iy-ny)
            do i=1,3
            IcoV(i,ix,iy,iz)=IcoV(i,ix,iy,iz)
     .     +IcoVx(i,nx,ny,iz)*w(ix-nx)*w(iy-ny)
            IcoF(i,ix,iy,iz)=IcoF(i,ix,iy,iz)
     .     +IcoFx(i,nx,ny,iz)*w(ix-nx)*w(iy-ny)
            enddo
            do j=1,4
            do i=1,4
            IcoT(i,j,ix,iy,iz)=IcoT(i,j,ix,iy,iz)
     .     +IcoTx(i,j,nx,ny,iz)*w(ix-nx)*w(iy-ny)
            enddo
            do i=1,3
            IcoC(i,j,ix,iy,iz)=IcoC(i,j,ix,iy,iz)
     .     +IcoCx(i,j,nx,ny,iz)*w(ix-nx)*w(iy-ny)
            enddo
            enddo
          enddo
        enddo
        enddo
      enddo
      enddo
      end

c-----------------------------------------------------------------------
      subroutine IcoStr(iflag,ierrico)
c-----------------------------------------------------------------------
c     energy momentum tensor and flavor current
c     and corresponding contractions from string method
c
c     iflag = 1 initialization
c           = 2 sum up density
c           = 3 normalization
c
c     output: common blocks /Ico1/ - /Ico5/
c-----------------------------------------------------------------------
#include "aaa.h"
#include "ico.h"
      double precision tpro,zpro,ttar,ztar,ttaus,detap,detat
      common/cttaus/tpro,zpro,ttar,ztar,ttaus,detap,detat
      common /Ico10/ntot,ish0,irs
      double precision v(4),u(4),eps,tt(4,4),gamma
      common/cranphi/ranphi
      common/cninx/ninx
      common/celoss/eloss
      common/cee1ico/ee1ico,eistico,ee1hll
      double precision  IcoTx, IcoCx,IcoTxx, IcoCxx
      real IcoEx, IcoVx, IcoFx,IcoExx, IcoVxx, IcoFxx
      real gint(3),g2int(3)
      common arryz(nyicomax,nzicomax),arryz2(nyicomax,nzicomax)
     .   ,arryz3(nyicomax,nzicomax),arryz4(nyicomax,nzicomax)
      common/cgint/gint,g2int,sigmad,sigmaz,lmbdad,lmbdaz,dnorma,znorma
      common/ciniflo/ecccoe(5),eccphi(5)
      dimension
     . IcoTx(4,4,nxicomax,nyicomax,nzicomax)
     .,IcoCx(3,4,nxicomax,nyicomax,nzicomax)
     .,IcoEx(nxicomax,nyicomax,nzicomax)
     .,IcoVx(3,nxicomax,nyicomax,nzicomax)
     .,IcoFx(3,nxicomax,nyicomax,nzicomax)
     .,IcoTxx(4,4,nxicomax,nyicomax,nzicomax)
     .,IcoCxx(3,4,nxicomax,nyicomax,nzicomax)
     .,IcoExx(nxicomax,nyicomax,nzicomax)
     .,IcoVxx(3,nxicomax,nyicomax,nzicomax)
     .,IcoFxx(3,nxicomax,nyicomax,nzicomax)
      common/credonoco/kredonoco
      real pskip(4),pcoro(4)
      common /cnoskip/pcore(4),qcore(4)
      save pskip,pcoro
      common/cchkengy/ichkengy/cchkengy2/esollxx,eistxx
      common /cnoskipadd/qcoreadd(4)
      real ar4(4),psum(4)
      parameter(ihhmax=5)
      real dhh1(ihhmax),dhh2(ihhmax)
      data dhh1 / 0.1, 0.3, -0.3,  0.3,  -0.3 /
      data dhh2 / 0.1, 0.3,  0.3, -0.3,  -0.3 /
      data ncntico/0/
      save ncntico

      nxicomx=nxicomax
      nyicomx=nyicomax
      nzicomx=nzicomax

      isegmin=0  !??????????????????????????????????????
      if(iabs(ninicon).ne.1)isegmin=0

      ish0=0                    ! set to zero if you don't want the output

      GOTO (1,10,100),iflag
      return

      !...............................................................
      !                       initialization
      !...............................................................

 1    do ix=1,nxico
        do iy=1,nyico
          do iz=1,nzico
            IcoE(ix,iy,iz)=0
            do i2=1,3
            IcoV(i2,ix,iy,iz)=0
            IcoF(i2,ix,iy,iz)=0
            enddo
            do i1=1,4
            do i2=1,3
            IcoC(i2,i1,ix,iy,iz)=0.
            enddo
            do i2=1,4
            IcoT(i1,i2,ix,iy,iz)=0.
            enddo
            enddo
          enddo
        enddo
      enddo
      ntot=0
      return

      !................................................................
      !                 fill arrays
      !................................................................

 10   continue
      call utpri('icostr',ish,ishini,4)
      ntot=ntot+1
      do ix=1,nxico
        do iy=1,nyico
          do iz=1,nzico
            do i1=1,4
            do i2=1,3
            IcoCx(i2,i1,ix,iy,iz)=0.
            enddo
            do i2=1,4
            IcoTx(i1,i2,ix,iy,iz)=0.
            enddo
            enddo
          enddo
        enddo
      enddo
      do iii=1,4
      pskip(iii)=0
      pcore(iii)=0
      qcore(iii)=0
      pcoro(iii)=0
      enddo
      do iy=1,nyico
      do iz=1,nzico
      arryz(iy,iz)=0
      arryz2(iy,iz)=0
      arryz3(iy,iz)=0
      enddo
      enddo
                   ! here s and t are just the sigma,tau
                   ! of the kinky-string model.
                   ! call them s,t to avoid confusion
                   ! with the ivariant tau
      gint(1)=0.68268
      gint(2)=0.95450
      gint(3)=0.99730
      g2int(1)=1-exp(-1/2.0)
      g2int(2)=1-exp(-4/2.0)
      g2int(3)=1-exp(-9/2.0)
      dxico=(xmaxico-xminico)/nxico
      dyico=(ymaxico-yminico)/nyico
      dzico=(zmaxico-zminico )/nzico
      !---------------------------
      sigmad=radeft !r-smearing
      if(sigmad.lt.0.999*dxico)stop'####### ERROR 06022016 #######'
      sigmaz=sgzico !eta-smearing
      !---------------------------
      lmbdad=nint(cutico)
      lmbdaz=3
      dnorma=(2.*pi)*sigmad**2           *g2int(lmbdad)
      znorma=(2.*pi)**0.5*sigmaz*tauzer  *gint(lmbdaz)
      nix=nint(lmbdad*sigmad/dxico)
      niy=nint(lmbdad*sigmad/dyico)
      niz=nint(lmbdaz*sigmaz/dzico)

      ttaus=tauzer
      zz=0.
      istrs=icostrg
      irmn=icoremn
      ispec=icospec
      nkmax=0
      nsmax=0
      do i=1,nptl
        if((istptl(i).eq.20.or.istptl(i).eq.21)
     &  .and.ityptl(i).ge.20.and.ityptl(i).le.39)nsmax=i
        if(istptl(i).eq.20.or.istptl(i).eq.21)nkmax=i
        if(istptl(i).eq.7)istptl(i)=5
      enddo


      !..........particles (spectators,remnants)..........................

      do i=1,nptl
        !~~~~~~~~~~spectators~~~~~~~~~~
        if(ispec.ge.1.and.istptl(i).eq.5.and.i.le.maproj+matarg)then
          irs=i
          call IcoAdd(1,i,i,1
     .           ,nxicomx,nyicomx,nzicomx
     .           ,IcoTx,IcoCx,IcoEx,IcoVx,IcoFx)
          if(ish0.ge.5)then
            write(ifmt,*)'Spectator'
            write(ifmt,'(i5,1x,i10,1x,19(1pg14.7,1x))')
     $       i,idptl(i),(pptl(k,i),k=1,5)
     $             ,(xorptl(k,i),k=1,4)
     $             ,sqrt(pptl(3,i)**2
     $             +pptl(2,i)**2+pptl(1,i)**2)
          endif
        endif
        !~~~~~~~~remnants~~~~~~~~~~~~~
        if(irmn.ge.1.and.
     &  istptl(i).eq.5.and.ityptl(i).ge.40)then
          if(iorptl(i).ne.0)then
          if(.not.(iorptl(i).ge.1.and.istptl(iorptl(i)).eq.29))then
             irs=i
             call IcoAdd(2,i,i,1
     .           ,nxicomx,nyicomx,nzicomx
     .           ,IcoTx,IcoCx,IcoEx,IcoVx,IcoFx)
             if(ish0.ge.3)then
               write(ifmt,*)'Remnant (resonance)'
               write(ifmt,'(i5,1x,i10,1x,19(1pg14.7,1x))')
     $          i,idptl(i),(pptl(k,i),k=1,5)
     $             ,(xorptl(k,i),k=1,4)
     $             ,sqrt(pptl(3,i)**2
     $             +pptl(2,i)**2+pptl(1,i)**2)
             endif
          endif
          endif
        endif
      enddo

      if(istrs.eq.0.and.irmn.eq.0)return   !strings from remnant

      nk1=1
      if(istrs.eq.0)nk1=nsmax+1
      if(irmn.eq.0)nkmax=nsmax

      if(nkmax.eq.0)return          !no string

      !...........strings.............................................

      nrseg=0
      do while(nk1.le.nkmax)
        if(istptl(nk1).eq.20.or.istptl(nk1).eq.21)then
          nk2=nk1+1
          do while(idptl(nk2).eq.9)
            nk2=nk2+1
          enddo
          ifr=ifrptl(1,nk1)
          do n=nk1+1,nk2
            if(ifrptl(1,n).ne.ifr)stop'in IcoStr: impossible'
          enddo
          irs=0
          n1=ifrptl(1,ifr)
          n2=ifrptl(2,ifr)
          nrseg=nrseg+n2-n1+1
          frac=0
          do i=n1,n2
            if(istptl(i).eq.5)frac=frac+1
          enddo
          ifrac=nint(frac)
          frac=frac/(n2-n1+1)
          if(frac.ge.1.1)then  !0.999   NEVER
            do i=n1,n2
              istptl(i)=7
              do iii=1,4
              pcore(iii)=pcore(iii)+pptl(iii,i)
              enddo
            enddo
            !write(ifch,*)'+++++',ifr,n1,n2,nint(frac*100),
            !&      '   ',(istptl(i),i=n1,n2)
          else
            !if(frac.gt.0.9999)
            !.      stop'\n\n STOP in IcoStr: should be treated as string \n\n'
            !print*,'++++ifrac++++',ifrac
            if(ifrac.lt.isegmin)then
              do i=n1,n2
                if(istptl(i).eq.5)then
                  do iii=1,4
                  pskip(iii)=pskip(iii)+pptl(iii,i)
                  enddo
                endif
              enddo
            else
              call IcoAdd(3,n1,n2,ifrac
     .               ,nxicomx,nyicomx,nzicomx
     .               ,IcoTx,IcoCx,IcoEx,IcoVx,IcoFx)
            endif
            goto54
          endif
          stop'\n\n30102011d\n\n'
 54       nk1=nk2+1
        else
          nk1=nk1+1
        endif
      enddo                     ! next string
      reno=1
      if(qcore(4).gt.0.)reno=pcore(4)/qcore(4)

      !~~~diagonalize T

      if(ish.ge.7)write(ifch,'(a)')
     .'----------------------- diagonalize T -----------------------'
      !call ooHoIniYEta(nyicomax,nzicomax,arryz
      !.,nyico,yminico, ymaxico,nzico,zminico,zmaxico,'Ico-1 ','segments')
      !call ooHoIniYEta(nyicomax,nzicomax,arryz2
      !.,nyico,yminico, ymaxico,nzico,zminico,zmaxico,'Ico-2 ','segments')
      !call ooHoIniYEta(nyicomax,nzicomax,arryz3
      !.,nyico,yminico, ymaxico,nzico,zminico,zmaxico,'Ico-3 ','Tmunu-00')
      do iy=1,nyico
      do iz=1,nzico
      arryz4(iy,iz)=IcoTx(4,4,nxIco/2+1,iy,iz)
      enddo
      enddo
      !call ooHoIniYEta(nyicomax,nzicomax,arryz4
      !.,nyico,yminico, ymaxico,nzico,zminico,zmaxico,'Ico-4 ','Tmunu-00')

      ee1=0
      ee2=0
      do ix=1,nxico
        do iy=1,nyico
          do iz=1,nzico
            do i=1,4
            do k=1,4
             tt(i,k)=IcoTx(i,k,ix,iy,iz)
            enddo
            enddo
            if(ish.ge.7)write(ifch,'(3i3,16f8.2)')ix,iy,iz
     .      ,((tt(i,k),k=1,4),i=1,4)
           !~~~~~~~~~~~~~~~~~~
            call DiagTmunu(tt,u,v,gamma,eps,nonz,ix,iy,iz)
           !~~~~~~~~~~~~~~~~~~
            if(nonz.eq.1)then
              eps=eps*reno
              if(pcore(4).ne.0.)then
                eps=eps*(pcore(4)+pskip(4)+eloss)/(pcore(4))
              endif
              forfac=1.0  !for testing /= 1.0, should be 1.0 for real runs
              eps=eps*forfac
              IcoEx(ix,iy,iz)=eps
              IcoVx(1,ix,iy,iz)=v(1)
              IcoVx(2,ix,iy,iz)=v(2)
              IcoVx(3,ix,iy,iz)=v(3) / tauzer
              do i=1,3
               IcoFx(i,ix,iy,iz)=IcoCx(i,4,ix,iy,iz)*u(4)
               do  k=1,3
                IcoFx(i,ix,iy,iz)=IcoFx(i,ix,iy,iz)
     .                  -IcoCx(i,k,ix,iy,iz)*u(k)
               enddo
              enddo
              if((    IcoFx(1,ix,iy,iz).gt.1e5
     .            .or.IcoFx(2,ix,iy,iz).gt.1e5
     .            .or.IcoFx(3,ix,iy,iz).gt.1e5)
     .         .and.50.gt.ncntico)then
               ncntico=ncntico+1
               write(ifmt,*)'ico: IcoFx_very_big: ',
     .          IcoFx(1,ix,iy,iz),IcoFx(2,ix,iy,iz),IcoFx(3,ix,iy,iz)
     .          ,'      ',u
               !write(ifmt,*)'      IcoFx(',ix,iy,iz,')=',
               !.          IcoFx(1,ix,iy,iz),IcoFx(2,ix,iy,iz),IcoFx(3,ix,iy,iz)
              endif
              xix=xminico+(ix-0.5)*(xmaxico-xminico)/nxico
              yiy=yminico+(iy-0.5)*(ymaxico-yminico)/nyico
              eta=zminico+(iz-0.5)*(zmaxico-zminico)/nzico
              u0=u(4)
              u3=u(3)
              press=eps/3.
              e=eps
              pi00=0.
              pi03=0.
              pi33=0.
              ee1=ee1+eflow(e,press,eta,u0,u3,pi00,pi03,pi33)
     .         *tauzer*dxico*dyico*dzico
              ee2=ee2+eps*cosh(eta)*tauzer*dxico*dyico*dzico
              !if(ix.eq.nxico/2+1.and.iy.eq.nyico/2+1)
              !.         print*,'+++++++',xix,yiy,eta,eps
            else
              IcoEx(ix,iy,iz)=0
              do i=1,3
               IcoVx(i,ix,iy,iz)=0
               IcoFx(i,ix,iy,iz)=0
              enddo
            endif
c            write(ifch,'(a,3i5,3f8.3)')'ico  x y z b123'
c     .     ,ix,iy,iz,(IcoFx(i3,ix,iy,iz),i3=1,3)
          enddo
        enddo
      enddo

      !~~~rescale~~~~~~~~

      call setParamsIco() 
      call core1paramget7(fico1,fico2,fico3,fico4,fico5,fico6,fico7) 
      Z=max(0.,ng1evt-2.)/fico7
      ficoscale = -1*
     . ( fico1 + min( fico2*Z**fico4 , fico5-fico6*(Z-1) ) )
      gfac=1+ficoscale  
      !factor to increase energy density compensate decrease of multiplicity due to hydro
      if(ee1.gt.0..and.abs(gfac-1).gt.0.00001)then
        isif= sign(1.,ficoscale)
        fc=2
        !first value (a):
        a=0 
        gscal=1-a*isif  
        call ScaleIcoX2XX(gscal,gfac
     .  ,nxicomx,nyicomx,nzicomx,IcoTx,IcoCx,IcoEx,IcoVx,IcoFx
     .  ,IcoTxx,IcoCxx,IcoExx,IcoVxx,IcoFxx, ee1xx , ee2xx )
        fa=ee1xx/ee1-1  
        !second value (b):
        b=0.7
        if(isif.lt.0)b=1.5
        gscal=1-b*isif       !narrower (broader) distribution vs eta 
        call ScaleIcoX2XX(gscal,gfac
     .  ,nxicomx,nyicomx,nzicomx,IcoTx,IcoCx,IcoEx,IcoVx,IcoFx
     .  ,IcoTxx,IcoCxx,IcoExx,IcoVxx,IcoFxx, ee1xx , ee2xx )
        fb=ee1xx/ee1-1   
        if(fb*isif.gt.0.)then !should be less than zero
          print*,'ERROR 18062017 rescale in ico: ee1xx ee1 isif:'
     .    , ee1xx,ee1,isif
          stop
        endif 
        iscal=0
        do while(iscal.lt.20.and.abs(fc).gt.0.001)
         iscal=iscal+1
         c=(a+b)/2 
         gscal=1-c*isif  ! < 1 makes narrower distributions vs eta
         call ScaleIcoX2XX(gscal,gfac
     .   ,nxicomx,nyicomx,nzicomx,IcoTx,IcoCx,IcoEx,IcoVx,IcoFx
     .   ,IcoTxx,IcoCxx,IcoExx,IcoVxx,IcoFxx, ee1xx , ee2xx )
         fc=ee1xx/ee1-1
         !print*,'ico rescale abc  ',a,b,c,'     ',fa,fb,fc
         if(fa*fc.le.0.)then
           b=c
           fb=fc
         else
           a=c
           fa=fc
         endif
        enddo
        if(abs(fc).gt.0.001)then
          write(ifch,*) 
     .    'WARNING rescale in ico: fc = ',fc,' -> too big!'
        else     
         do ix=1,nxico
          do iy=1,nyico
           do iz=1,nzico
                       IcoEx(ix,iy,iz)       = IcoExx(ix,iy,iz)
             do i2=1,3
                       IcoVx(i2,ix,iy,iz)    = IcoVxx(i2,ix,iy,iz)
                       IcoFx(i2,ix,iy,iz)    = IcoFxx(i2,ix,iy,iz)
             enddo
             do i1=1,4
             do i2=1,3
                       IcoCx(i2,i1,ix,iy,iz) = IcoCxx(i2,i1,ix,iy,iz)
             enddo
             do i2=1,4
                       IcoTx(i1,i2,ix,iy,iz) = IcoTxx(i1,i2,ix,iy,iz)
             enddo
             enddo
           enddo
          enddo
         enddo
        endif
        write(ifmt,'(a,i7,2f8.2)')'ico rescale ',iscal,gscal-1,gfac
      endif

      !~~~get random phi

      call RotateIcoX2XX(phievt
     .,nxicomx,nyicomx,nzicomx,IcoTx,IcoCx,IcoEx,IcoVx,IcoFx
     .,IcoTxx,IcoCxx,IcoExx,IcoVxx,IcoFxx)
      call IniFlowIni(nxicomx,nyicomx,nzicomx,IcoExx)
      if(iranphi.eq.1)then
       ranphi=eccphi(2)
      else
       ranphi=0
      endif
      if(ninx.eq.1)write(ifmt,*)'+++++ranphi+++++',ranphi,eccphi(2)
      phinll=phievt+ranphi
      call RotateIcoX2XX(phinll
     .,nxicomx,nyicomx,nzicomx,IcoTx,IcoCx,IcoEx,IcoVx,IcoFx
     .,IcoTxx,IcoCxx,IcoExx,IcoVxx,IcoFxx)

      !~~~~~~sum

      do ix=1,nxico
        do iy=1,nyico
          do iz=1,nzico
              IcoE(ix,iy,iz)
     .       =IcoE(ix,iy,iz)+IcoExx(ix,iy,iz)
              do i2=1,3
              IcoV(i2,ix,iy,iz)
     .       =IcoV(i2,ix,iy,iz)+IcoVxx(i2,ix,iy,iz)
              IcoF(i2,ix,iy,iz)
     .       =IcoF(i2,ix,iy,iz)+IcoFxx(i2,ix,iy,iz)
              enddo
              do i1=1,4
              do i2=1,3
              IcoC(i2,i1,ix,iy,iz)
     .       =IcoC(i2,i1,ix,iy,iz)+IcoCxx(i2,i1,ix,iy,iz)
              enddo
              do i2=1,4
              IcoT(i1,i2,ix,iy,iz)
     .       =IcoT(i1,i2,ix,iy,iz)+IcoTxx(i1,i2,ix,iy,iz)
              enddo
              enddo
          enddo
        enddo
      enddo

      if(ninx.eq.1)then
  111   format(1x,a,2f13.1,2f15.1)
        write(ifmt,111)'+++++pskip++++++',pskip
        write(ifmt,111)'+++++pcore++++++',pcore
        write(ifmt,111)'+++++qcore++++++',qcore
        write(ifmt,'(1x,a,41x,f15.1)')'+++++eloss++++++',eloss
        espec=0
        do i=1,nptlpt
          if(istptl(i).eq.0)espec=espec+pptl(4,i)
        enddo
        write(ifmt,'(1x,a,41x,f15.1)')'+++++espec++++++',espec
        do iii=1,4
         pcoro(iii)=0
         do i=1,nptl
          if(istptl(i).eq.0)pcoro(iii)=pcoro(iii)+pptl(iii,i)
         enddo
        enddo
        write(ifmt,111)'+++++pcoro++++++',pcoro
        einit=maproj*engy/2+matarg*engy/2
        if(pcore(4).ne.0.)then
        efina=pskip(4)+pcore(4)+pcoro(4)+eloss
        else
        efina=pskip(4)+pcoro(4)
        endif
        write(ifmt,111)'+++++pall+++++++'
     . ,pskip(1)+pcore(1)+pcoro(1)
     . ,pskip(2)+pcore(2)+pcoro(2)
     . ,pskip(3)+pcore(3)+pcoro(3)
     . ,efina
        write(ifmt,'(1x,a,41x,f15.1)')'+++++eini+++++++',einit
        do j=1,4
          psum(j)=0
        enddo
        do i=1,maproj+matarg
          call  getpptl(i,ar4(1),ar4(2),ar4(3),ar4(4),dum)
          do j=1,4
            psum(j)=psum(j)+ar4(j)
          enddo
        enddo
        write(ifmt,111)'+++++pini+++++++',psum
        ecoeloss=pcore(4)+eloss
        if(pcore(4).ne.0.)
     .  write(ifmt,'(1x,a,41x,f15.1,a)')'+++ecore+eloss++',ecoeloss,' *'
        write(ifmt,'(1x,a,41x,f15.1)')  '+++EfluidSimpl++',ee2
        write(ifmt,'(1x,a,41x,f15.1,a)')'+++++Efluid+++++',ee1,' *'
        ee1ico=ee1

        if(pcore(4).ne.0.)then
        egain=ee1-ecoeloss
        write(ifmt,'(1x,a,41x,f15.1)')'+++++egain++++++',egain
        esoll=pcoro(4)-egain
        eist=pcoro(4)
        esollxx=esoll
        eistxx=eist
        ichkengy=1
        call chkengy(ierrchk)
        if(ierrchk.eq.1)then
          ierrico=1
          goto 1001
        endif
        ichkengy=2
        else
        ichkengy=2
        esollxx=-1
        endif
        eistico=esollxx
        !call checkengy('after rescaling ')

        if(abs(einit-efina).gt.0.01*einit.and.irescl.ge.1)then
          write(ifmt,*)'REDO; ico: Eall /= Eini:',efina,einit
          ierrico=1
          if(pcore(4).eq.0.)kredonoco=1
          goto 1001
        endif

        if(abs(reno-1.).gt.0.2
     .   .and.abs(pcore(4)-qcore(4)).gt.3.)then
         if(pcore(4).eq.0.)then
          write(ifmt,*)'WARNING ico: pcore(4) = 0 '
          ierrico=1
          kredonoco=1
          goto 1001
         endif
        endif

        if(abs(reno-1.).gt.0.2
     .   .and.abs(pcore(4)-qcore(4)).gt.3.)then
          if(ish.ge.0)
     .    write(ifmt,*)'WARNING  reno = ',pcore(4),'/',qcore(4)
     .    ,'=',reno,'   => REDO'
          jerr(18)=jerr(18)+1  ! pcore(4) .NE. qcore(4)
              ierrico=1
              goto 1001
        endif

        if(pcore(4).ne.0.
     .    .and.(abs(ee1-ecoeloss).gt.0.10*ecoeloss
     .    .and.abs(ee1-ecoeloss).gt.0.10*einit))then
          write(ifmt,*)'REDO; ico: int{T00} =',ee1
     .      ,'  NE  ecore+eloss =',ecoeloss
          ierrico=1
          if(pcore(4).eq.0.)kredonoco=1
          goto 1001
        endif

      endif

 1001 continue
      call utprix('icostr',ish,ishini,4)
      return

      !......................................................................
c                 normalization
      !......................................................................

 100  ish2=0                    ! set to 1 if you want screen-output

      za=abs(ninicon)
      do ix=1,nxico
        do iy=1,nyico
          do iz=1,nzico
            IcoE(ix,iy,iz)
     .     =IcoE(ix,iy,iz)/za
            do i2=1,3
              IcoV(i2,ix,iy,iz)
     .       =IcoV(i2,ix,iy,iz)/za
            enddo
            do i2=1,3
              IcoF(i2,ix,iy,iz)
     .       =IcoF(i2,ix,iy,iz)/za
            enddo
            do i1=1,4
              do i2=1,3
                IcoC(i2,i1,ix,iy,iz)
     .         =IcoC(i2,i1,ix,iy,iz)/za
              enddo
              do i2=1,4
                IcoT(i1,i2,ix,iy,iz)
     .         =IcoT(i1,i2,ix,iy,iz)/za
              enddo
            enddo
          enddo
        enddo
      enddo
      return
      end

c--------------------------------------------------------------------------
      subroutine chkengy(ierrchk)
c--------------------------------------------------------------------------
#include "aaa.h"
      common/cchkengy/ichkengy/cchkengy2/esollxx,eistxx
      ierrchk=0
      if(ichkengy.eq.2)return
      if(ichkengy.eq.0)then
        if(esollxx.lt.0.)return
      endif
      esoll=esollxx
      if(ichkengy.eq.1)then
        eist=eistxx
      elseif(ichkengy.eq.0)then
         eist=0
         do i=1,nptl
          if(istptl(i).eq.0)eist=eist+pptl(4,i)
         enddo
      endif
      if(abs(esoll/eist-1.).lt.0.0001)return
      ntrygain=1
 7777 ntrygain=ntrygain+1
      if(ntrygain.gt.100)then
        write(ifmt,*)'REDO; ico: ntrygain to big'
        ierrchk=1
        return
      endif
      scalegain=esoll/eist
      !write(ifmt,*)'   ntry scale',ntrygain,scalegain
      eist=0
      do i=1,nptl
        if(istptl(i).eq.0)then
          amt2=pptl(5,i)**2+pptl(1,i)**2+pptl(2,i)**2
          pptl(3,i)= pptl(3,i)*scalegain
          pptl(4,i)=sqrt(pptl(3,i)**2+amt2)
          eist=eist+pptl(4,i)
        endif
      enddo
      if(abs(esoll/eist-1.).gt.0.01)goto 7777
      if(mod(nrevt+1,modsho).eq.0)
     .write(ifmt,'(a,f8.0,a,i3,f8.3)')
     . '   ecoro',eist,'   ntry scale',ntrygain,scalegain
      end

c--------------------------------------------------------------------------
      subroutine storepqcore
c--------------------------------------------------------------------------
      common /cnoskip/pcore(4),qcore(4)
      common /cnoskipxx/pcorexx(4),qcorexx(4)
      do ii=1,4
        qcorexx(ii)=qcore(ii)
        pcorexx(ii)=pcore(ii)
      enddo
      end

c--------------------------------------------------------------------------
      subroutine restorepqcore
c--------------------------------------------------------------------------
      common /cnoskip/pcore(4),qcore(4)
      common /cnoskipxx/pcorexx(4),qcorexx(4)
      do ii=1,4
        qcore(ii)=qcorexx(ii)
        pcore(ii)=pcorexx(ii)
      enddo
      end

c-----------------------------------------------------------------------
      subroutine IcoAdd(jtyp,n1,n2,ifrac
     .               ,nxicomx,nyicomx,nzicomx
     .               ,IcoTx,IcoCx,IcoEx,IcoVx,IcoFx)
c-----------------------------------------------------------------------
      ! jtyp: 1=spectator  2=remnant  3=string
      !-----------------------------------------------------------------
#include "aaa.h"
#include "ico.h"
      common/cxyzt/xptl(mxptl+29),yptl(mxptl+29),zptl(mxptl+29)
     * ,tptl(mxptl+29),optl(mxptl+29),uptl(mxptl+29),sptl(mxptl+29)
     *,rptl(mxptl+29,3)
      dimension xo(4),ps(5)
      double precision  IcoTx, IcoCx
      real IcoEx, IcoVx, IcoFx
      parameter(mxrapseg=1000)
      real rapseg(mxrapseg),rapsegxx(mxrapseg),etasegxx(mxrapseg)
      integer iord(mxrapseg),iflxx(3,mxrapseg)
      real facreno(2),facrenox(2)
      dimension
     . IcoTx(4,4,nxicomax,nyicomax,nzicomax)
     .,IcoCx(3,4,nxicomax,nyicomax,nzicomax)
     .,IcoEx(nxicomax,nyicomax,nzicomax)
     .,IcoVx(3,nxicomax,nyicomax,nzicomax)
     .,IcoFx(3,nxicomax,nyicomax,nzicomax)
      real gint(3),g2int(3)
      common/cgint/gint,g2int,sigmad,sigmaz,lmbdad,lmbdaz,dnorma,znorma
      common /cnoskip/pcore(4),qcore(4)
      common /cnoskipadd/qcoreadd(4)
      parameter(ihhmax=5)
      real dhh1(ihhmax),dhh2(ihhmax)
      data dhh1 / 0.1, 0.3, -0.3,  0.3,  -0.3 /
      data dhh2 / 0.1, 0.3,  0.3, -0.3,  -0.3 /

      if(ifrac.eq.0)return

      call utpri('icoadd',ish,ishini,4)

      etamin=100000
      etamax=-100000
      rapmin=100000
      rapmax=-100000
      ipl=0
      imi=0
      if(ish.ge.6)then
        write(ifch,*)' enter IcoAdd; particle list:',ifrac
        do n=n1,n2
          write(ifch,*)n,idptl(n),istptl(n)
        enddo
      endif
      jjmin=0
      do i=n1,n2
        ii=i-n1+1
        call idflav(idptl(i), iflxx(1,ii), iflxx(2,ii), iflxx(3,ii)
     .   ,idummy,idummy)
      enddo
      if(rangen().lt.0.5)then
        ii1=1
        ii2=n2-n1+1
        iii=1
      else
        ii1=n2-n1+1
        ii2=1
        iii=-1
      endif
      iannihilate=0
      if(iannihilate.eq.1)then!~~~~
       do ii=ii1,ii2,iii
        i=n1+ii-1
        if(istptl(i).eq.5)then
          if(iflxx(2,ii).eq.-iflxx(3,ii))then
            iflxx(2,ii)=0
            iflxx(3,ii)=0
          endif
          if(ii.ne.ii2)then
            do jj=ii+iii,ii2,iii
              j=n1+jj-1
              if(istptl(j).eq.5)then
                do ki=1,3
                do kj=1,3
                  if(iflxx(ki,ii).eq.-iflxx(kj,jj))then
                    iflxx(ki,ii)=0
                    iflxx(kj,jj)=0
                  endif
                enddo
                enddo
              endif
            enddo
          endif
        endif
       enddo
      endif!~~~~
      do ilo=1,2 !~~~~~~~~~~~~~~ilo~~~~~~~~~~~~~~~~~
        ietarap=1
        ifull=0
        if(ilo.eq.2)then
          ifull=1
          ietarap=0
          do j=1,n2-n1+1
           rapmn=2*1e10
           do jj=1,n2-n1+1
             if(rapseg(jj).lt.rapmn)then
               rapmn=rapseg(jj)
               jjmin=jj
             endif
           enddo
           if(jjmin.eq.0)then
             write(ifmt,*)'##### ERROR 12062015 : '
     .       ,n1,n2,(rapseg(jjtst),jjtst=1,n2-n1+1)
             stop
           endif
           iord(j)=jjmin
           rapseg(jjmin)=3*1e10
           !print*,'+++iord+++',j,iord(j)
          enddo
        endif
        ntryico=0
        facreno(1)=1.
        facreno(2)=1.
        facrenox(1)=0.
        facrenox(2)=0.
        ihh=0
        nnn=0
        eee1=0.
        eee2=0.
        eee0=0.
        eee00=0.
        iout=0
        dh1=0.
        dh2=0.
        call storepqcore
  11    continue
        ntryico=ntryico+1
        nord=0
        if(ilo.eq.2)then
          if(ntryico.gt.1.or.iout.gt.0)then
            do ii=n1,n2
              i=ii
              if(istptl(i).eq.7)istptl(i)=5
            enddo
            call restorepqcore
          endif
          dra=sigmaz
          rapmid=
     .    (rapmin+dra*facreno(1)+rapmax-dra*facreno(2))/2.
          rapminx=rapmin+dra*facreno(1)
          rapminx=min(rapminx,rapmid-dra/2.)
          rapmaxx=rapmax-dra*facreno(2)
          rapmaxx=max(rapmaxx,rapmid+dra/2.)
        endif
        !write(ifmt,*)'+++++++',ntryico
        ist57=0
        nn2=0
        do ii=n1,n2
          if(ilo.eq.1)i=ii
          if(ilo.eq.2)i=n1-1+iord(ii-n1+1)
          if(istptl(i).eq.5)n2x=ii
        enddo
        do ii=n1,n2 !~~~~~~~~~~~~ii~~~~~~~~~~~~~
          if(ilo.eq.1)i=ii
          if(ilo.eq.2)i=n1-1+iord(ii-n1+1)
          if(ilo.eq.1)rapseg(ii-n1+1)=10000+ii-n1+1
          if(istptl(i).eq.5)then
            ist57=1
            if(ilo.eq.2)istptl(i)=7
            !call idflav(idptl(i),ifl1,ifl2,ifl3,idummy,idummy)
            ifl1=iflxx(1,i-n1+1)
            ifl2=iflxx(2,i-n1+1)
            ifl3=iflxx(3,i-n1+1)
            if(mod(leadcore,10).eq.1)then !ifl123=0 for strings
             if(jtyp.eq.3)then !strings
             !~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
             !we overwrite the ifl, because they should essentially
             !annihilate (always qqbar from string breaks), but they
             !dont, because some string segments are moving out
             !of the fluid, others not. The latter ones should not take
             !their flavors with them, but rather pick up flavor at the
             !fluid surface.
             !~~~~~~~~~~~~~~~~~~~~~~~~~~~~
             ifl1=0
             ifl2=0
             ifl3=0
             endif
            endif
            !if(ilo.eq.2.and.ntryico.eq.1)then
            ! print*,' IcoAdd ifl ',i,ifl1,ifl2,ifl3
            !endif
            !call jtain(i,xo(1),xo(2),xo(3),xo(4),nn,1)
            xo(1)=xptl(i)
            xo(2)=yptl(i)
            xo(3)=zptl(i)
            xo(4)=tptl(i)
            irs=i
            do iii=1,5
              ps(iii)=pptl(iii,i)
            enddo
            if(ilo.eq.2)then
              do iii=1,4
              pcore(iii)=pcore(iii)+ps(iii)
              enddo
              nord=nord+1
              if(nord.eq.1)then
                esum1=0
                do iii=1,4
                  qcoreadd(iii)=0.
                enddo
               endif
              !eta=etamin+(nord-1)*(etamax-etamin)/(ifrac-1)
              if(rapmin.lt.100000.and.rapmax.gt.-100000)then
                if(ifrac.gt.1)then
                  eta=rapminx+(nord-1)*(rapmaxx-rapminx)/(ifrac-1)
                else
                  eta=(rapminx+rapmaxx)/2
                endif
              else !no valid rapidity found
                write(ifmt,'(a)')
     .          'WARNING icoadd Eta changed to avoid crash'
                eta=0 !avoid crash in cosh(eta)
              endif
              amt=ps(5)**2+ps(1)**2+ps(2)**2
              if(amt.gt.0.)then
                amt=sqrt(amt)
              else
                amt=0
              endif 
              e2=amt*cosh(eta)
              esum1=esum1+ps(4)
              rapxx=rapsegxx(i-n1+1)
              etaxx=etasegxx(i-n1+1)
              if(ish.ge.5.and.nnn.eq.1)
     .        write(ifch,'(a,3f8.2,i6,2(3x,4f8.2))')
     .         ' ++     eta etaxx rapxx ',eta,etaxx,rapxx,i
            !.         ,pptl(3,i)
            !.         ,amt*cosh(rap),e2
            !.         ,rapminx ,rapmaxx
     .        ,xptl(i),yptl(i),zptl(i),tptl(i)  ,xo
              !print*,'IcoAdd',jtyp,'    ',ifl1,ifl2,ifl3
            endif
            call IcoAddptl(ps,xo,ifl1,ifl2,ifl3,ietarap,ifull
     .       ,eta,rap,nxicomx,nyicomx,nzicomx
     .       ,IcoTx,IcoCx,IcoEx,IcoVx,IcoFx,nnn,ii-n2x,iout)
            if(ilo.eq.1)then
              if(i-n1+1.gt.mxrapseg)
     .         stop'\n\n ERROR 24062011\n\n'
              rapseg(i-n1+1)=rap
              rapsegxx(i-n1+1)=rap
              etasegxx(i-n1+1)=eta
              !print*,'+++rapseg+++',i-n1+1,rapseg(i-n1+1)
              if(eta.gt.etamax)ipl=ipl+1
              if(eta.lt.etamin)imi=imi+1
              etamin=min(etamin,eta)
              etamax=max(etamax,eta)
              if(rap.gt.-1e9.and.rap.lt.1e9)then
                rapmin=min(rapmin,rap)
                rapmax=max(rapmax,rap)
              endif
            endif
          endif
        enddo  !~~~~~~~~~~~ii~~~~~~~~~~~
        if(ilo.eq.2)then
          esum2=qcoreadd(4)
          if(iout.eq.1)then !ptl outside ico grid
            iout=2 
            write(ifmt,*)'redo'
            ntryico=0
            goto 11
          endif
          eee=esum1/esum2-1  !should be zero
          ddd=abs(esum1-esum2)
          if(ish.ge.6)
     .     write(ifch,'(a,i2,3f7.2,3x,2f8.2,2f10.2,3x
     .     ,2f7.2,i8)')
     .    ' IcoAddptl-',ntryico,facreno(1)
     .    ,facreno(2),eee,esum1,esum2,qcore(4),pcore(4)
     .    ,rapminx,rapmaxx,n2-n1+1
          limi=25
          if(ntryico.lt.limi
     .     .and.abs(eee).gt.0.01
     .     .and.(ddd.gt.1.0.or.ntryico<12))then
            if(ntryico.gt.98)then
              write(ifmt,*)'ERROR: too many REDO'
              stop'\n\n ERROR 16102011 \n\n'
            endif
            imo=mod(ntryico-1,3)+1
            if(imo.eq.1)then
              if(ntryico.eq.1)ihh=1
              if(ntryico.eq.1)eee00=eee
              dh1=dhh1(ihh)
              dh2=dhh2(ihh)
              if(ntryico.gt.1)then
                if(((abs(eee1-eee0).lt.0.05*abs(eee0))
     .            .and.(abs(eee2-eee0).lt.0.05*abs(eee0)))
     .           .or.(abs(eee).gt.abs(eee00)   ) )then
                  facreno(1)=facrenox(1)
                  facreno(2)=facrenox(2)
                  ihh=1+ihh
                  ihh=1+mod(ihh-1,ihhmax)
                  dh1=dhh1(ihh)
                  dh2=dhh2(ihh)
                endif
              endif
              facrenox(1)=facreno(1)
              facrenox(2)=facreno(2)
              eee0=eee
              facreno(1)=facrenox(1)+dh1
              facreno(2)=facrenox(2)
              goto 11
            elseif(imo.eq.2)then
              eee1=eee
              facreno(1)=facrenox(1)
              facreno(2)=facrenox(2)+dh2
              goto 11
            elseif(imo.eq.3)then
              eee2=eee
            endif
            !E0+grad{E}*\vec{h}=0
            !\vec{h}=-grad{E}*h
            !->h=E0/grad{E}^2
            !grad{E}=((E1-E0)/dh1,E2-E0/dh2)
            ge2=((eee0-eee1)/dh1)**2+((eee0-eee2)/dh2)**2
            if( (abs(eee1-eee0).lt.0.05*abs(eee0))
     .      .and.(abs(eee2-eee0).lt.0.05*abs(eee0)) )then
              facreno(1)=facrenox(1)
              facreno(2)=facrenox(2)
              goto 11
            endif
            h1=(eee0-eee1)/dh1*eee0/ge2
            h2=(eee0-eee2)/dh2*eee0/ge2
            facreno(1)=facrenox(1)+h1
            facreno(2)=facrenox(2)+h2
            facreno(1)=max(facreno(1),-4.)
            facreno(2)=max(facreno(2),-4.)
            facreno(1)=min(facreno(1),4.)
            facreno(2)=min(facreno(2),4.)
            goto 11
          endif
          if(nnn.eq.0)then
            if(ntryico.eq.limi)then
              facreno(1)=1
              facreno(2)=1
            endif
            nnn=1
            goto 11
          endif
          if(ish.ge.6)write(ifch,*)'  done'
        endif
      enddo !~~~~~~~~~~~~~~~ilo~~~~~~~~~~~~~~~~~~
      if(ish.ge.6)then
        write(ifch,*)' exit IcoAdd; particle list:'
        do n=n1,n2
          write(ifch,*)n,idptl(n),istptl(n)
        enddo
      endif
      call utprix('icoadd',ish,ishini,4)
      end

c--------------------------------------------------------------------------
      subroutine IcoAddptl(ps,xo,ifl1,ifl2,ifl3,ietarap,ifull,eta,rap
     .,nxicomx,nyicomx,nzicomx,IcoTx,IcoCx,IcoEx,IcoVx,IcoFx,nnn
     .,iin2,iout)
c--------------------------------------------------------------------------
      dimension ps(5),xo(4),ifl(3),pp(4),qq(4),nfl(3)
      common /Ico10/ntot,ish0,irs
#include "aaa.h"
#include "ico.h"
      double precision  IcoTx, IcoCx
      real IcoEx, IcoVx, IcoFx
      common /cnoskip/pcore(4),qcore(4)
      common /cnoskipxx/pcorexx(4),qcorexx(4)
      common /cnoskipadd/qcoreadd(4)
      logical ldum
      common arryz(nyicomax,nzicomax),arryz2(nyicomax,nzicomax)
     .   ,arryz3(nyicomax,nzicomax),arryz4(nyicomax,nzicomax)
      dimension
     . IcoTx(4,4,nxicomx,nyicomx,nzicomx)
     .,IcoCx(3,4,nxicomx,nyicomx,nzicomx)
     .,IcoEx(nxicomx,nyicomx,nzicomx)
     .,IcoVx(3,nxicomx,nyicomx,nzicomx)
     .,IcoFx(3,nxicomx,nyicomx,nzicomx)
      real gint(3),g2int(3)
      common/cgint/gint,g2int,sigmad,sigmaz,lmbdad,lmbdaz,dnorma,znorma
      !real qq(4)
      data ncnticoadd/0/
      save ncnticoadd

c Just to avoid warnings with gfortran
      ldum=.false.
      if(ldum)print *,IcoEx,IcoVx,IcoFx

c compute eta, rap

      if(ietarap.eq.1)then
        if(xo(3).ge.xo(4))then
          eta=1e10
        elseif(xo(3).le.-xo(4))then
          eta=-1e10
        else
          eta= 0.5d0*log((xo(4)+xo(3))/(xo(4)-xo(3)))
        endif
        amt=ps(5)**2+ps(1)**2+ps(2)**2
        if(amt.gt.0..and.ps(4)+abs(ps(3)).gt.0.d0)then
          amt=sqrt(amt)
          rap=sign(1.,ps(3))*alog((ps(4)+abs(ps(3)))/amt)
        else
         amt=0
         rap=0
        endif 
      endif

      if(ifull.eq.0)return

c tests
      if(nnn.eq.1)then
        iiz=  (eta-zminIco)/(zmaxIco-zminIco )*nzIco+1
        iiy=(xo(2)-yminIco)/(ymaxIco-yminIco)*nyIco+1
        iix=(xo(1)-xminIco)/(xmaxIco-xminIco)*nxIco+1
        if(iix.eq.nxIco/2+1.and.iiy.ge.1.and.iiy.le.nyIco
     .  .and.iiz.ge.1.and.iiz.le.nzIco)arryz(iiy,iiz)=arryz(iiy,iiz)+1
      endif

c take ps(1..4) and do contribution to grid-points
      ifl(1)=ifl1
      ifl(2)=ifl2
      ifl(3)=ifl3
      nfl(1)=0
      nfl(2)=0
      nfl(3)=0
      do i=1,3
        if(ifl(i).ne.0)then
          iafl=abs(ifl(i))
          if(iafl.le.3)nfl(iafl)=nfl(iafl)+sign(1,ifl(i))
        endif
      enddo
      amt=ps(5)**2+ps(1)**2+ps(2)**2
      if(amt.gt.0..and.ps(4)+abs(ps(3)).gt.0.d0)then
        amt=sqrt(amt)
        rap=sign(1.,ps(3))*alog((ps(4)+abs(ps(3)))/amt)
      else
       amt=0
       rap=0
      endif 
      delz=(zmaxIco-zminIco )/nzIco
      delx=(xmaxIco-xminIco)/nxIco
      dely=(ymaxIco-yminIco)/nyIco
      sum=0
      wwsum=0
      qqsum=0
      ncnticoadd=ncnticoadd+1
      if(iout.ge.1)then !redo due to ptl outside ico grid
        xo(1)=max( xo(1) , xminIco+(float(  1  )-0.5)*delx )
        xo(1)=min( xo(1) , xminIco+(float(nxIco)-0.5)*delx )
        xo(2)=max( xo(2) , yminIco+(float(  1  )-0.5)*dely )
        xo(2)=min( xo(2) , yminIco+(float(nyIco)-0.5)*dely )
        !write(ifmt,*)'xo(1/2) = ',xo(1),xo(2)
      endif
      do iz=1,nzIco
        zz=zminIco +(float(iz)-0.5)*delz
        z2=( eta - zz )**2
        !.......boost ps() into co-moving frame with rapidity=zz
        pp(1)=ps(1)
        pp(2)=ps(2)
        pp(3)=amt*sinh(zz)
        pp(4)=amt*cosh(zz)
        qq(1)=pp(1)
        qq(2)=pp(2)
        qq(3)=pp(3)
        qq(4)=pp(4)
        !WWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWW
        !if(ncnticoadd.lt.7)then
        !www=0
        !if(z2.le.lmbdaz**2*sigmaz**2)
        !.   www=   exp(-z2/2./sigmaz**2) / znorma
        ! wwsum=wwsum+www*delz*tauzer
        ! qqsum=qqsum+www*delz*tauzer*qq(4)
        !write(ifch,'(a,2f7.2,4x,f10.2,4x,3f8.2,4x,2f10.2,$)')
        !.   '       ',zz,eta - zz, qq(4)
        !.    ,www,www*delz*tauzer,wwsum
        !.   ,  www*delz*tauzer*qq(4),qqsum
        !endif
        !WWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWW
        !after boost:
        pp(3)=0
        pp(4)=amt
        do ix=1,nxIco
          xx=xminIco+(float(ix)-0.5)*delx
          do iy=1,nyIco
            yy=yminIco+(float(iy)-0.5)*dely
            d2=(xo(1)-xx)**2+(xo(2)-yy)**2
            if(d2.gt.lmbdad**2*sigmad**2.
     &            or.z2.gt.lmbdaz**2*sigmaz**2)then
              w=0.
            else
              w= exp(-d2/2./sigmad**2) / dnorma
     &         * exp(-z2/2./sigmaz**2) / znorma
              sum=sum+w*delx*dely*delz*tauzer
              if(pp(4).gt.0.)then
                do iii=1,4
                  dqc=qq(iii)*w*delx*dely*delz*tauzer
                  qcore(iii)=qcore(iii)+dqc
                  qcoreadd(iii)=qcoreadd(iii)+dqc
                enddo
                if(nnn.eq.1)then !~~~~~nnn~~~~~
                  !...particles...
                  if(ix.eq.nxIco/2+1)arryz2(iy,iz)=arryz2(iy,iz)+w
                  !....engy-momentum tensor...
                  if(ix.eq.nxIco/2+1)arryz3(iy,iz)=arryz3(iy,iz)+
     &                   w*pp(4)
                  do i1=1,4
                    do i2=1,4
                      IcoTx(i1,i2,ix,iy,iz)=IcoTx(i1,i2,ix,iy,iz)+
     $                   w*pp(i1)*pp(i2)/pp(4)
                    enddo
                  enddo
                  !....velocity field of particles....
                  do iafl=1,3
                    do i1=1,4
                      f=nfl(iafl)
                      IcoCx(iafl,i1,ix,iy,iz)=IcoCx(iafl,i1,ix,iy,iz)+
     $                f*w*pp(i1)/pp(4)
                    enddo
                  enddo
                endif !~~~~~nnn~~~~~~
              endif
            endif
          enddo
        enddo
        !if(ncnticoadd.lt.7)write(ifch,*)'    ',qcoreadd(4),sum
      enddo
      !print*,'+++sum+++++++',sum
      if(iout.eq.0)then
       if((xo(1).lt.xminIco+(float(  1  )-0.5)*delx
     .  .or.xo(1).gt.xminIco+(float(nxIco)-0.5)*delx
     .  .or.xo(2).lt.yminIco+(float(  1  )-0.5)*dely
     .  .or.xo(2).gt.yminIco+(float(nyIco)-0.5)*dely)
     .  .and.sum.eq.0..and.iin2.eq.0)then
        write(ifmt,'(a,6f6.2,$)')
     .  'WARNING ptl outside ico grid:'
     .  ,xminIco+(float(1)-0.5)*delx
     .  ,xminIco+(float(nxIco)-0.5)*delx
     .  ,yminIco+(float(1)-0.5)*dely
     .  ,yminIco+(float(nyIco)-0.5)*dely , xo(1),xo(2)
        iout=1
       endif
      endif

      end

c-----------------------------------------------------------------------------
c subroutine DiagTmunuRedu(Tmunu,eps)
c subroutine DiagTmunuFullRedu(Tmunu,eps) !
c subroutine DiagTmunu(Tmunu,u,v,gamma,eps,nonz,ix,iy,iz)
c subroutine DiagTmunuFull(Tmunu,u,v,gamma,eps,nonz,ix,iy,iz)
c-----------------------------------------------------------------------------
c T^mu^nu(x) = sum_{particles} p^mu*p^nu/p^0 * f_smooth(x-x_particle) 
c-----------------------------------------------------------------------------
c solve lambda * T * lambda = diag(eps,P,P,P), for symmetric T (T.Kodama and KW)
c 
c lambda^a_b = g      g*vx        g*vy        g*vz          =symmetric!
c              g*vx   a*vx*vx+1   a*vy*vx     a*vz*vx
c              g*vy   a*vx*vy     a*vy*vy+1   a*vz*vy
c              g*vz   a*vx*vz     a*vy*vz     a*vz*vz+1
c with g=gamma, a=g*g/(g+1)
c
c so T*lambda(v)=lambda(-v)*diag(eps,P,P,P)
c first column: four equations
c                    eps=T00+T0k*vk
c                -eps*vi=Ti0+Tik*vk
c solved iteratively to get eps, vi
c returns u,v,eps if nonz=1 (otherwise zero tensor)
c Method in PHSD "diagonalizing Tmunu" is equivalent, but only 
c  concerning ONE eigenvalue (epsilon), solving T*g*u=epsilon*u
c------------------------------------------------------------------------------
c DiagTmunu:
c ----------
c Consider only longitudinal boost (v1=v2=0) 
c => two equations eps=T00+T03*v3  and   -eps*v3=T30+T33*v3
c =>  -(T00+T03*v3)*v3=T30+T33*v3
c =>  T00*v3 + T03*v3**2 + T30 + T33*v3=0
c =>  T03*v3**2 + (T00+T33)*v3 + T30 = 0
c =>  v3**2 + (T00+T33)/T03*v3 + 1 = 0
c------------------------------------------------------------------------------

      !--------------------------------------------------
      !   Consider longitudinal boost, reduced arguments:
      !--------------------------------------------------
      subroutine DiagTmunuRedu(Tmunu,eps)
      double precision Tmunu(4,4),eps
      double precision u(4),v(4),gamma
      integer nonz,ix,iy,iz
      ix=0
      iy=0
      iz=0
      eps=0
      call DiagTmunu(Tmunu,u,v,gamma,eps,nonz,ix,iy,iz)
      end

      !-------------------------------------------
      !    Consider full boost, reduced arguments:
      !-------------------------------------------
      subroutine DiagTmunuFullRedu(Tmunu,eps)
      double precision Tmunu(4,4),eps
      double precision u(4),v(4),gamma
      integer nonz,ix,iy,iz
      ix=0
      iy=0
      iz=0
      eps=0
      call DiagTmunuFull(Tmunu,u,v,gamma,eps,nonz,ix,iy,iz)
      end

      !-----------------------------------
      !   Consider longitudinal boost:
      !-----------------------------------
      subroutine DiagTmunu(Tmunu,u,v,gamma,eps,nonz,ix,iy,iz)
#include "aaa.h"
#include "ico.h"
      double precision Tmunu(4,4),u(4),v(4),gamma,eps
      integer nonz,ix,iy,iz
      double precision beta,b

      nonz=0
      sum=0d0
      do i=1,4
       do k=1,4
       sum=sum+dabs(Tmunu(i,k))
       end do
      end do
      if(sum.eq.0.0d0)return

      nonz=1
      if(abs(Tmunu(4,3)).le.1e-20)then
        beta=0
      else
        b=(Tmunu(4,4)+Tmunu(3,3))/Tmunu(4,3)
        if(abs(b).lt.2)then
          write(ifmt,*)'DiagTmunu(',ix,iy,iz,'):   b < 2'
          nonz=0
          return
        endif
        if(b.ge.0)then
          beta=(b-sqrt(b**2-4)) / 2
        else
          beta=(b+sqrt(b**2-4)) / 2
        endif
      endif
      if(.not.(beta.le.0..or.beta.ge.0.))stop'DiagTmunu NaN catch'

      v(1)=0
      v(2)=0
      v(3)=beta
      v2=0.d0
      do i=1,3
        v2=v2+v(i)**2
      enddo
      if(v2.eq.1.0)then
        gamma=ainfin
        eps=0.
        u(4)=0.
        u(1)=0.
        u(2)=0.
        u(3)=0.
        v(1)=0.
        v(2)=0.
        v(3)=0.
      else
        gamma=1./sqrt(abs(1.-v2))
        eps=
     .  (Tmunu(4,4) - 2*beta*Tmunu(4,3) + beta**2*Tmunu(3,3))*gamma**2
        u(4)=gamma
        u(1)=v(1)*gamma
        u(2)=v(2)*gamma
        u(3)=v(3)*gamma
      endif

      end

      !---------------------------
      !    Consider full boost:
      !---------------------------
      subroutine DiagTmunuFull(Tmunu,u,v,gamma,eps,nonz,ix,iy,iz)
#include "aaa.h"
#include "ico.h"
      double precision Tmunu(4,4),u(4),v(4),gamma,eps
      integer nonz,ix,iy,iz
      double precision g,a, Lor(4,4),w(4),sum
      double precision vx(3),tt(4,4),err,sg(4)

      sg(4)=1d0
      v(4)=1d0
      do i=1,3
       sg(i)=-1d0
      enddo
      sum=0d0
      do i=1,4
       do k=1,4
       sum=sum+dabs(Tmunu(i,k))
       end do
      end do
      nonz=0
      if(sum.eq.0.0d0)return
      nonz=1

      do k=1,3
       v(k)=0.
      end do
      eps=0

      do lrep=1,100
       epsx=eps
       do i=1,3
        vx(i)=v(i)
       end do
       eps=Tmunu(4,4)
       do k=1,3
        eps=eps-Tmunu(4,k)*vx(k)
       end do
       if(eps.le.0d0)then
         if(.not.(eps.gt.-1e-2.or.eps.gt.-1e-2*sum))then
           am2=Tmunu(4,4)**2
     .       -Tmunu(4,3)**2-Tmunu(4,2)**2-Tmunu(4,1)**2
           write(ifmt,*)'DiagTmunu: ***** negative epsilon *****'
           write(ifmt,*)'  ix,iy,iz:',ix,iy,iz
           write(ifmt,*)'  sum(abs(Tmunu))=',sngl(sum)
     .      ,'   m**2=',am2,'   eps=',sngl(eps)
           write(ifmt,*)'  Tmunu(4,nu):',(sngl(Tmunu(4,nu)),nu=1,4)
         endif
         nonz=0
         return
       endif
       do i=1,3
        Tv=0d0
        do k=1,3
         Tv=Tv+Tmunu(i,k)*vx(k)
        end do
        v(i)=(Tmunu(i,4)-Tv)/eps
       end do
       if(lrep.gt.60)then
        do i=1,3
         v(i)=0.5d0*(vx(i)+v(i))
        enddo
       endif
       !write(ifmt,*)'Tmunu: ',lrep,abs(eps-epsx),(abs(v(i)-vx(i)),i=1,3)
       err=1d-6
       if(lrep.gt.50)err=1d-5
       if(lrep.gt.89)err=1d-4
       if(abs(eps-epsx).lt.err.and.abs(v(1)-vx(1)).lt.err
     . .and.abs(v(2)-vx(2)).lt.err.and.abs(v(3)-vx(3)).lt.err)goto1
        do i=1,4
          w(i)=0
          do k=1,4
          w(i)=w(i)+Tmunu(i,k)*v(k)*sg(k)
          enddo
          w(i)=w(i)-eps*v(i)
        enddo
        if(lrep.gt.95
     ..and.w(1)*w(1)+w(2)*w(2)+w(3)*w(3)+w(4)*w(4).lt.err)goto1
      end do

  1   v2=0.d0
      do i=1,3
        v2=v2+v(i)**2
      enddo
      if(v2.ge.1.)then
        if(eps.gt.1e-2)then
          write(ifmt,*)'DiagTmunu: ***** v2 ge 1 ***** '
          write(ifmt,*)'DiagTmunu:   v2=',v2,'    eps=',eps
          write(ifmt,*)'DiagTmunu:   v=',v
          write(ifmt,*)'DiagTmunu:   ix,iy,iz=',ix,iy,iz
        else
          nonz=0
          return
        endif
      endif
      gamma=1./sqrt(abs(1.-v2))
      u(4)=gamma
      u(1)=v(1)*gamma
      u(2)=v(2)*gamma
      u(3)=v(3)*gamma

      !~~~check
      g=gamma
      a=g*g/(g+1d0)
      Lor(4,4)=g
      do k=1,3
       Lor(4,k)=-g*v(k)
       Lor(k,4)=Lor(4,k)
      enddo
      do i=1,3
      do k=1,3
       Lor(i,k)=a*v(i)*v(k)
      enddo
      enddo
      do k=1,3
       Lor(k,k)=Lor(k,k)+1
      enddo
      do i=1,4
      do k=1,4
        tt(i,k)=0d0
        do m=1,4
        do n=1,4
          tt(i,k)=tt(i,k)+Lor(i,m)*Tmunu(m,n)*Lor(n,k)
        enddo
        enddo
      enddo
      enddo
      err=err*1d2
      if(tt(1,4).gt.err.or.tt(2,4).gt.err.or.tt(3,4).gt.err)then
        write(ifmt,'(2a,i5,a)')
     .   ' ************ nonzero T14 or T24 or T34 '
     .  ,'after ',lrep,' iterations'
        write(ifmt,*)'ix,iy,iz=',ix,iy,iz
        do i=1,4
        write(ifmt,'(4f9.5,2x,4f9.5)')
     .  (sngl(tmunu(i,k)),k=1,4),(sngl(tt(i,k)),k=1,4)
        enddo
      endif

      end

c----------------------------------------------------------------------------
      subroutine xstr(nk1,nk2,sigma,tau,x,gp,gm,ixy)
c----------------------------------------------------------------------------
c     calculates position if string as a function of sigma and tau
c     returns coordinate vector x(1..4) and velocity vector g(1..4)
c     string defined from nk1 to nk2 in pptl list
c----------------------------------------------------------------------------

#include "aaa.h"
      logical pr
      dimension x(4),gp(4),gm(4)
      mmod(i)=mod(mod(i,nob)+nob,nob)

      data pr/.false./
      spt=sigma+tau
      smt=sigma-tau
c      if(sigma.gt.19.975) pr=.true.
      if(pr)print *,'sigma-region:',sigma,tau,smt,spt
      a=0.
      i=0
      nob=(nk2-nk1+1)*2
      ii=i+nk1
      if(i.ge.nob/2) ii=nk1+(nob-1)-i
c.....positioning of index i
      do while( smt-a.gt.pptl(4,ii) .or. smt-a.lt.0. )
        if(smt-a.lt.0.)then
          i=mmod(i-1)
          ii=i+nk1
          if(i.ge.nob/2) ii=nk1+(nob-1)-i
          a=a-pptl(4,ii)
          if(pr)print *,'skip -1:',i,ii,nob,a,smt,pptl(4,ii)
        else
          a=a+pptl(4,ii)
          i=mmod(i+1)
          ii=i+nk1
          if(i.ge.nob/2) ii=nk1+(nob-1)-i
          if(pr)print *,'skip +1:',i,ii,nob,a,smt
        endif
      enddo
c.....smt,smp are in one band ?
      if(spt-a.lt.pptl(4,ii))then !   on band-case
        if(pr)print *,'one:',i,(spt-smt),sigma,tau
        i1=ii
        i2=ii
        do k=1,4
          x(k)=pptl(k,ii)*((spt-smt))/pptl(4,ii)/2
          gp(k)=pptl(k,ii)/pptl(4,ii)
          gm(k)=0.
        enddo
      else
c.....  loop to band of spt
        if(pr)print *,'first:',i,(pptl(4,ii)-(smt-a))/pptl(4,ii)
        i1=ii
        do k=1,4
          x(k)=pptl(k,ii)*(pptl(4,ii)-(smt-a))/pptl(4,ii)/2
          gp(k)= 0.5*pptl(k,ii)/pptl(4,ii)
          gm(k)=-0.5*pptl(k,ii)/pptl(4,ii)
        enddo
        a=a+pptl(4,ii)
        i=mmod(i+1)
        ii=i+nk1
        if(i.ge.nob/2) ii=nk1+(nob-1)-i
        do while( spt-a.ge.pptl(4,ii) )
          if(pr)print *,'med:',i,a,1.0
          do k=1,4
            x(k)=x(k)+pptl(k,ii)/2
          enddo
          a=a+pptl(4,ii)
          i=mmod(i+1)
          ii=i+nk1
          if(i.ge.nob/2) ii=nk1+(nob-1)-i
        enddo
        if(pr)print *,'last:',i,a,(spt-a)/pptl(4,ii)
        i2=ii
        do k=1,4
          x(k)=x(k)+pptl(k,ii)*(spt-a)/pptl(4,ii)/2
          gp(k)=gp(k)+0.5*pptl(k,ii)/pptl(4,ii)
          gm(k)=gm(k)+0.5*pptl(k,ii)/pptl(4,ii)
        enddo
      endif
      if(pr)print *,'x:',sigma,tau,x
      if(pr)write (*,'("gap:",12(g12.6,1x))') gp
      if(pr)write (*,'("gam:",12(g12.6,1x))') gm
      ixy=i1+nob*i2
      end

c---------------------------------------------------------------------
      subroutine RotateIcoX2XX(phi
     .,nxicomx,nyicomx,nzicomx,IcoTx,IcoCx,IcoEx,IcoVx,IcoFx
     .,IcoTxx,IcoCxx,IcoExx,IcoVxx,IcoFxx)
c---------------------------------------------------------------------
#include "aaa.h"
#include "ico.h"
      double precision aa(2,2)
      double precision  IcoTx, IcoCx,IcoTxx, IcoCxx
      real IcoEx, IcoVx, IcoFx,IcoExx, IcoVxx, IcoFxx
      dimension
     . IcoTx(4,4,nxicomx,nyicomx,nzicomx)
     .,IcoCx(3,4,nxicomx,nyicomx,nzicomx)
     .,IcoEx(nxicomx,nyicomx,nzicomx)
     .,IcoVx(3,nxicomx,nyicomx,nzicomx)
     .,IcoFx(3,nxicomx,nyicomx,nzicomx)
     .,IcoTxx(4,4,nxicomx,nyicomx,nzicomx)
     .,IcoCxx(3,4,nxicomx,nyicomx,nzicomx)
     .,IcoExx(nxicomx,nyicomx,nzicomx)
     .,IcoVxx(3,nxicomx,nyicomx,nzicomx)
     .,IcoFxx(3,nxicomx,nyicomx,nzicomx)
      double precision  gi(2),gj(2)
      delx=(xmaxico-xminico)/nxico
      dely=(ymaxico-yminico)/nyico
      do ix=1,nxico
        xp=xminico+(ix-0.5)*(xmaxico-xminico)/nxico
        do iy=1,nyico
          yp=yminico+(iy-0.5)*(ymaxico-yminico)/nyico
          do iz=1,nzico
            xv=  cos(phi)*xp - sin(phi)*yp
            yv=  sin(phi)*xp + cos(phi)*yp
            i=(xv-(xminico-0.5*delx))/delx
            i=max(1,i)
            i=min(nxico-1,i)
            xi=xminico+(i-0.5)*(xmaxico-xminico)/nxico
            frac=(xv-xi)/delx
            gi(1)=1-frac
            gi(2)=frac
            j=(yv-(yminico-0.5*dely))/dely
            j=max(1,j)
            j=min(nyico-1,j)
            yj=yminico+(j-0.5)*(ymaxico-yminico)/nyico
            frac=(yv-yj)/dely
            gj(1)=1-frac
            gj(2)=frac
            IcoExx(ix,iy,iz)=0
            do i2=1,3
            IcoVxx(i2,ix,iy,iz)=0
            IcoFxx(i2,ix,iy,iz)=0
            enddo
            do i1=1,4
            do i2=1,3
            IcoCxx(i2,i1,ix,iy,iz)=0
            enddo
            enddo
            do i1=1,4
            do i2=1,4
            IcoTxx(i1,i2,ix,iy,iz)=0
            enddo
            enddo
            dlta=-0.001
            if(gi(1).ge.dlta .and. gi(2).ge.dlta .and.
     .         gj(1).ge.dlta .and. gj(2).ge.dlta       )then
             do n=1,2
             do m=1,2
              IcoExx(ix,iy,iz)=IcoExx(ix,iy,iz)
     .                            +gi(n)*gj(m)*IcoEx(i+n-1,j+m-1,iz)
              do i2=1,3
              IcoVxx(i2,ix,iy,iz)=IcoVxx(i2,ix,iy,iz)
     .                         +gi(n)*gj(m)*IcoVx(i2,i+n-1,j+m-1,iz)
              IcoFxx(i2,ix,iy,iz)=IcoFxx(i2,ix,iy,iz)
     .                         +gi(n)*gj(m)*IcoFx(i2,i+n-1,j+m-1,iz)
              enddo
              do i1=1,4
              do i2=1,3
                IcoCxx(i2,i1,ix,iy,iz)=IcoCxx(i2,i1,ix,iy,iz)
     .                      +gi(n)*gj(m)*IcoCx(i2,i1,i+n-1,j+m-1,iz)
              enddo
              enddo
              do i1=1,4
              do i2=1,4
                IcoTxx(i1,i2,ix,iy,iz)=IcoTxx(i1,i2,ix,iy,iz)
     .                      +gi(n)*gj(m)*IcoTx(i1,i2,i+n-1,j+m-1,iz)
              enddo
              enddo
             enddo
             enddo
            endif
            aa(1,1)= cos(phi)
            aa(1,2)= sin(phi)
            aa(2,1)=-sin(phi)
            aa(2,2)= cos(phi)
            !~~~~~~~~~~~~
            a=0
            b=0
            do m=1,2
            a=a+aa(1,m)*IcoVxx(m,ix,iy,iz)
            b=b+aa(2,m)*IcoVxx(m,ix,iy,iz)
            enddo
            IcoVxx(1,ix,iy,iz)=a
            IcoVxx(2,ix,iy,iz)=b
            !~~~~~~~~~~~~
            do i2=1,3
            a=0
            b=0
            do m=1,2
            a=a+aa(1,m)*IcoCxx(i2,m,ix,iy,iz)
            b=b+aa(2,m)*IcoCxx(i2,m,ix,iy,iz)
            enddo
            IcoCxx(i2,1,ix,iy,iz)=a
            IcoCxx(i2,2,ix,iy,iz)=b
            enddo
            !~~~~~~~~~~~~
            a=0
            b=0
            c=0
            d=0
            do m=1,2
            do n=1,2
            a=a+aa(1,m)*aa(1,n)*IcoTxx(m,n,ix,iy,iz)
            b=b+aa(1,m)*aa(2,n)*IcoTxx(m,n,ix,iy,iz)
            c=c+aa(2,m)*aa(1,n)*IcoTxx(m,n,ix,iy,iz)
            d=d+aa(2,m)*aa(2,n)*IcoTxx(m,n,ix,iy,iz)
            enddo
            enddo
            IcoTxx(1,1,ix,iy,iz)=a
            IcoTxx(1,2,ix,iy,iz)=b
            IcoTxx(2,1,ix,iy,iz)=c
            IcoTxx(2,2,ix,iy,iz)=d
            !~~~~~~~~~~~~~~~
            vv=IcoVxx(1,ix,iy,iz)**2+IcoVxx(2,ix,iy,iz)**2
     .       +IcoVxx(3,ix,iy,iz)**2*tauzer**2
            if(vv.gt.0.)vv=sqrt(vv)
            if(vv.gt.1.)then
             if(IcoExx(ix,iy,iz).gt.0.01)then
             !write(ifmt,*)'WARNING: velocty bigger than 1 -> 1 ; '
             !.       ,' energy density: ',IcoExx(ix,iy,iz)
             endif
             do i2=1,3
             IcoVxx(i2,ix,iy,iz)=IcoVxx(i2,ix,iy,iz)/vv*0.99
             enddo
            endif
          enddo
        enddo
      enddo
      end

c---------------------------------------------------------------------
      subroutine ScaleIcoX2XX(gscal,gfac
     .,nxicomx,nyicomx,nzicomx,IcoTx,IcoCx,IcoEx,IcoVx,IcoFx
     .,IcoTxx,IcoCxx,IcoExx,IcoVxx,IcoFxx, ee1xx ,ee2xx)
c---------------------------------------------------------------------
#include "aaa.h"
#include "ico.h"
      double precision  IcoTx, IcoCx,IcoTxx, IcoCxx
      real IcoEx, IcoVx, IcoFx,IcoExx, IcoVxx, IcoFxx
      dimension
     . IcoTx(4,4,nxicomx,nyicomx,nzicomx)
     .,IcoCx(3,4,nxicomx,nyicomx,nzicomx)
     .,IcoEx(nxicomx,nyicomx,nzicomx)
     .,IcoVx(3,nxicomx,nyicomx,nzicomx)
     .,IcoFx(3,nxicomx,nyicomx,nzicomx)
     .,IcoTxx(4,4,nxicomx,nyicomx,nzicomx)
     .,IcoCxx(3,4,nxicomx,nyicomx,nzicomx)
     .,IcoExx(nxicomx,nyicomx,nzicomx)
     .,IcoVxx(3,nxicomx,nyicomx,nzicomx)
     .,IcoFxx(3,nxicomx,nyicomx,nzicomx)
      double precision  gi(2)
      real etadi(8:38)
      do i=8,38
        etadi(i)=0
      enddo
      !choose method (1 or 2) 
      methodico=1 !2
      method=methodico
      if(gscal.le.1.)method=1
      nloop=method
      weight=1./nloop
      !definitions
      dxico=(xmaxico-xminico)/nxico
      dyico=(ymaxico-yminico)/nyico
      dzico=(zmaxico-zminico )/nzico
      delz=(zmaxico-zminico)/nzico
      !do rescaling
      ee1xx=0
      ee2xx=0
      do ix=1,nxico
        do iy=1,nyico
          do iz=1,nzico
            IcoExx(ix,iy,iz)=0
            do i2=1,3
            IcoFxx(i2,ix,iy,iz)=0
            enddo
            do i1=1,4
            do i2=1,3
            IcoCxx(i2,i1,ix,iy,iz)=0
            enddo
            enddo
            do i1=1,4
            do i2=1,4
            IcoTxx(i1,i2,ix,iy,iz)=0
            enddo
            enddo
            do i2=1,3
              IcoVxx(i2,ix,iy,iz)=IcoVx(i2,ix,iy,iz) 
            enddo
            zxx=zminico+(iz-0.5)*delz
            do nlo=1,nloop 
              isi=3-2*nlo
              if(method.eq.1)then 
                z=zxx/gscal  !rescale grid
              elseif(method.eq.2)then 
                 z=zxx+isi*5*(gscal-1) !shift left and right
              endif
              i=(z-(zminico-0.5*delz))/delz
              i=max(1,i)
              i=min(nzico-1,i)
              zi=zminico+(i-0.5)*delz
              frac=(z-zi)/delz
              gi(1)=1-frac
              gi(2)=frac
              dlta=-0.001
              if(gi(1).ge.dlta .and. gi(2).ge.dlta)then
              do n=1,2
                IcoExx(ix,iy,iz)=IcoExx(ix,iy,iz)
     .                      +gi(n)*IcoEx(ix,iy,i+n-1) *gfac *weight
                do i2=1,3
                IcoFxx(i2,ix,iy,iz)=IcoFxx(i2,ix,iy,iz)
     .                         +gi(n)*IcoFx(i2,ix,iy,i+n-1) *weight
                enddo
                do i1=1,4
                do i2=1,3
                IcoCxx(i2,i1,ix,iy,iz)=IcoCxx(i2,i1,ix,iy,iz)
     .                      +gi(n)*IcoCx(i2,i1,ix,iy,i+n-1) *weight
                enddo
                enddo
                do i1=1,4
                do i2=1,4
                  IcoTxx(i1,i2,ix,iy,iz)=IcoTxx(i1,i2,ix,iy,iz)
     .                      +gi(n)*IcoTx(i1,i2,ix,iy,i+n-1) *weight
                enddo
                enddo
              enddo
              endif
            enddo !nlo 
            v2=IcoVxx(1,ix,iy,iz)**2+IcoVxx(2,ix,iy,iz)**2
     .                +IcoVxx(3,ix,iy,iz)**2
            if(v2.lt.1.)then
              gamma=1./sqrt(1.-v2)
              u0=gamma
              u3=IcoVxx(3,ix,iy,iz)*gamma
              eps=IcoExx(ix,iy,iz)
              if(iz.ge.8.and.iz.le.38)etadi(iz)=etadi(iz)+eps
              press=eps/3.
              e=eps
              pi00=0.
              pi03=0.
              pi33=0.
              !take original eta grid value:
              eta=zminico+(iz-0.5)*delz 
              !=>distribution of values rescaled
              !  gscal < 1 =>  narrower distribution 
              ee1xx=ee1xx+eflow(e,press,eta,u0,u3,pi00,pi03,pi33)
     .         *tauzer*dxico*dyico*dzico
              ee2xx=ee2xx+eps*cosh(eta)*tauzer*dxico*dyico*dzico
            endif
          enddo
        enddo
      enddo
      !write(ifmt,'(31i6)')(nint(etadi(i)),i=8,38)
      end



c-----------------------------------------------------------------------
      subroutine checkengy(text)
c-----------------------------------------------------------------------
#include "aaa.h"
      character text*16, txtwar*7
      common/cee1ico/ee1ico,eistico,ee1hll
      real sum(4)
      txtwar='       '
      einit=maproj*engy/2+matarg*engy/2
      do i=1,4
        sum(i)=0
      enddo
      do n=1,nptl
        call getistptl(n,ist)
        call getityptl(n,ity)
        if(ist.eq.0)then
          call getpptl(n,p1,p2,p3,p4,p5)
          sum(1)=sum(1)+p1
          sum(2)=sum(2)+p2
          sum(3)=sum(3)+p3
          sum(4)=sum(4)+p4
        endif
      enddo
      txtwar='       '
      if(sum(4)/einit.gt.1.25)txtwar='WARNING'
      !if(ish.ge.1.or.txtwar.ne.'       ')then
        write(ifmt,'(a,$)')'checkengy '
        write(ifmt,'(2a,2f5.1,2f9.0,f6.2,2a)')text,' ist=0 p='
     .  ,sum,sum(4)/einit,' ',txtwar
      !endif
      end

c-----------------------------------------------------------------------
      subroutine xIco3
c-----------------------------------------------------------------------
#include "aaa.h"
#include "ico.h"
      logical ok
      dxico=(xmaxico-xminico)/nxico
      dyico=(ymaxico-yminico)/nyico
      dzico=(zmaxico-zminico )/nzico
c      if(  zclass(3,1).lt.0.0001
c     ..and.zclass(3,3).lt.0.0001
c     ..and.zclass(3,5).lt.0.0001)return
      delz=(zmaxico-zminico)/nzico
      nmax=nzico
      vmin=zminico
      vmax=zmaxico
      delv=delz
c      do ii=1,3
      ok=.true. !.false.
c      if( ii.eq.1.and.bimevt.ge.zclass(1,1).and.bimevt.le.zclass(2,1)
c     ..or.ii.eq.2.and.bimevt.ge.zclass(1,3).and.bimevt.le.zclass(2,3)
c     ..or.ii.eq.3.and.bimevt.ge.zclass(1,5).and.bimevt.le.zclass(2,5))
c     .ok=.true.
      write(ifhi,'(a)')       '!##################################'
      write(ifhi,'(a,i3)')    '!   epsi       '
      write(ifhi,'(a)')       '!##################################'
      write(ifhi,'(a)')'openhisto name epsi'
      write(ifhi,'(a)')'xrange -9 11'
      write(ifhi,'(a)')'yrange 1 1000'
c      if(ii.eq.1)write(ifhi,'(a)')'txt  "title centr 1"'
c      if(ii.eq.2)write(ifhi,'(a)')'txt  "title centr 3"'
c      if(ii.eq.3)write(ifhi,'(a)')'txt  "title centr 5"'
      write(ifhi,'(a)') 'htyp lin xmod lin ymod log '
      write(ifhi,'(a)')'txt  "xaxis [c] "'
      write(ifhi,'(a)')'txt "yaxis sum [e] "'
      hiwe=abs(ninicon)
      if(.not.ok)hiwe=0d0
      write(ifhi,'(a,e22.14)')'histoweight ',hiwe
      write(ifhi,'(a)')'array 2'
      do n=1,nmax
        v=vmin+(n-0.5)*delv
        if(ok)then
        w=0
        do nx=1,nxico
        do ny=1,nyico
        w=w+IcoE(nx,ny,n)*tauzer*dxico*dyico*dzico
        enddo
        enddo
        else
        w=0
        endif
        write(ifhi,'(2e13.5)')v,w
      enddo
      write(ifhi,'(a)') 'endarray closehisto '
      write(ifhi,'(a)') 'plot 0'
c      enddo
      end

c-----------------------------------------------------------------------
      subroutine xIco(nev)
c-----------------------------------------------------------------------
#include "aaa.h"
#include "ico.h"
#include "ems.h"
      real IcE(20,nxicomax,nzicomax)
      integer nIcE(20)
      save IcE,nIcE
      if(abs(ninicon).gt.1)return
      ixmin=1      !1+nxico/8*2
      ixmax=nxico  !1+nxico/8*6
      iymin=1      !1+nyico/8*2
      iymax=nyico  !1+nyico/8*6
      izmin=1      !1+nzico/6*1
      izmax=nzico  !1+nzico/6*5
      nx=ixmax-ixmin+1
      ny=iymax-iymin+1
      nz=izmax-izmin+1
      xmin=xminico+(ixmin-1)*(xmaxico-xminico)/nxico
      ymin=yminico+(iymin-1)*(ymaxico-yminico)/nyico
      zmin=zminico+(izmin-1)*(zmaxico-zminico)/nzico
      xmax=xminico+(ixmax-0)*(xmaxico-xminico)/nxico
      ymax=yminico+(iymax-0)*(ymaxico-yminico)/nyico
      zmax=zminico+(izmax-0)*(zmaxico-zminico)/nzico
      if(nev.lt.0)goto777
      if(nev.eq.1)then
        do k=1,5
          nIcE(k)=0
          do ix=ixmin,ixmax
            do iz=izmin,izmax
              IcE(k,ix,iz)=0
            enddo
          enddo
        enddo
      endif
      kk=0
      if(izmode.eq.1)then
        do k=1,100
         if(zclass(3,k).gt.0.and.
     .    bimevt.ge.zclass(1,k).and.bimevt.le.zclass(2,k))then
          kk=k
         endif
        enddo
        write(ifmt,'(a,i4)')'centrality class',kk
      elseif(izmode.eq.2)then
        zp=nprt(1)
        do k=1,100
         if(zclass(3,k).gt.0.and.
     .    zp.ge.zclass(1,k).and.zp.le.zclass(2,k))then
          kk=k
         endif
        enddo
        write(ifmt,'(a,i4,a,i10,a)')'centrality class',kk
     .   ,' ; ', nprt(1),' elementary scattering(s)'
      elseif(izmode.eq.3)then
        write(ifmt,*)nprt(1),' NN scatterings'
        zp=nglevt
        do k=1,5
         if(zclass(3,k).gt.0.and.
     .    zp.ge.zclass(1,k).and.zp.le.zclass(2,k))then
          kk=k
         endif
        enddo
      elseif(izmode.eq.4)then
        zp=segevt
        do k=1,100
         if(zclass(3,k).gt.0.and.
     .    zp.ge.zclass(1,k).and.zp.le.zclass(2,k))then
          kk=k
         endif
        enddo
        write(ifmt,'(a,i4,a,f6.2)')'centrality class',kk
     .   ,' ; segment plateau height = ', zp
      else
        stop'in xIco: wrong izmode\n\n'
      endif
      if(kk.ge.1.and.kk.le.5)then
        nIcE(kk)=nIcE(kk)+1
        do ix=ixmin,ixmax
          do iz=izmin,izmax
            IcE(kk,ix,iz)=IcE(kk,ix,iz)+IcoE(ix,nyico/2+1,iz)
          enddo
        enddo
      endif
      return

 777  continue
      do k=1,5
        do ix=ixmin,ixmax
          do iz=izmin,izmax
            if(nIcE(k).ne.0)IcE(k,ix,iz)=IcE(k,ix,iz)/nIcE(k)
          enddo
        enddo
      enddo
      write(ifhi,'(a)')   '!----------------------------------'
      write(ifhi,'(a,i3)')'!             xIco       '
      write(ifhi,'(a)')   '!----------------------------------'
      write(ifhi,'(a)')   '!newpage'
      do k=1,5
        do iz=izmin,izmax,2
        z=zminico+(iz-0.5)*(zmaxico-zminico)/nzico
        if(iz.le.9)write(ifhi,'(a,i1,i1)')'openhisto name epsi',k,iz
        if(iz.gt.9)write(ifhi,'(a,i1,i2)')'openhisto name epsi',k,iz
        write(ifhi,'(a)')'htyp lin xmod lin ymod log'
        if(iz.eq.izmin)then
          write(ifhi,'(a,2e11.3)')'xrange',xmin,xmax
          write(ifhi,'(a,2e11.3)')'yrange 0.001 auto'
          write(ifhi,'(a)')'txt  "xaxis x (fm)"'
          write(ifhi,'(a)')'txt "yaxis [e] (GeV/fm^3!)"'
          write(ifhi,'(a,f5.2,a)')'text 0.30 0.88 "eta=',z,'"'
          write(ifhi,'(a)')'text 0.07 0.07 "y=0"'
        elseif(iz.gt.izmax-2)then
          write(ifhi,'(a,f5.2,a)')'text 0.70 0.88 "...',z,'"'
        endif
        write(ifhi,'(a,d22.14)')'histoweight ',float(nIcE(k))
        write(ifhi,'(a)')'array 2'
        do ix=ixmin,ixmax
          x=xminico+(ix-0.5)*(xmaxico-xminico)/nxico
          w=IcE(k,ix,iz)
          write(ifhi,'(2e13.5)')x,w
        enddo
        write(ifhi,'(a)') 'endarray closehisto '
        if(iz.gt.izmax-2)then
          write(ifhi,'(a)') 'plot 0'
        else
          write(ifhi,'(a)') 'plot 0-'
        endif
        enddo
      enddo
      end

c-----------------------------------------------------------------------
      subroutine xIco2
c-----------------------------------------------------------------------
#include "aaa.h"
#include "ico.h"
      common/ciprint2d/iprint2d

      if(iprint2d.eq.1)then

      dx=(xmaxico-xminico)/nxico
      dy=(ymaxico-yminico)/nyico
      dz=(zmaxico-zminico)/nzico
      di=8
      ixmin=1+nxico/2-8
      ixmax=1+nxico/2+8
      iymin=1+nyico/2-8
      iymax=1+nyico/2+8
      izmin=1+nzico/2
      izmax=1+nzico/2

      do iz=izmin,izmax

      z=zminico+(iz-0.5)*dz
      write(ifch,*)'=== print2d (ico.f) === IcoE -- eta = ',z
      write(ifch,'(7x,99i4)')
     .(nint( (xminico+(ix-0.5)*dx)*100 ),ix=ixmin,ixmax)
      do iy=iymax,iymin,-1
        iiy=nint( (yminico+(iy-0.5)*dy)*100 )
        write(ifch,'(i5,2x,99i4)')iiy
     .  ,(nint(IcoE(ix,iy,iz)),ix=ixmin,ixmax)
      enddo

      enddo

      endif
      end

c-----------------------------------------------------------------------
      subroutine xIniCon(word)
c-----------------------------------------------------------------------
c this unit creates root files for plots
c-----------------------------------------------------------------------
c the argument "word" is some string of the form:
c      X(N)
c where X is some text without spaces, and N some integer number.
c The text X is a keyword, used in the following as
c     if(word(1:length).eq.X)then
c            ...
c     endif
c where length is the length of the string X.  The number N is used as 
c parameter for creating plots (called ival in the program). In principle 
c we want to visualise three-dimensional arrays, f(-nx:nx,-ny:ny,-nz:nz). 
c To do so, we plot for example f(ix,iy,0) vs ix and iy. To have a 
c family of plots we use the parameter N, to plot f(ix,iy,N) for different
c values of N.
c-----------------------------------------------------------------------
c The section  if(word(1:length).eq.X)then ... endif is activated by 
c putting the following in the optns file:
c
c      root->X(N)
c
c with a valid keyword X (see below in the code)
c-----------------------------------------------------------------------

#include "aaa.h"
#include "ico.h"
        common /Ico10/ntot,ish0,irs
        character word*(*)
        character fmt*8

        length1=index(word,';')-1
        fmt='(a6,a  )'
        if(length1.le.9)then
          write(fmt(6:6),'(i1)')length1
        elseif(length1.le.99)then
          write(fmt(6:7),'(i2)')length1
        else
          stop'ERROR 29032021'
        endif

        nzoff=6
        nyoff=7

        zminoff=zminico+nzoff*(zmaxico-zminico)/nzico
        zmaxoff=zmaxico-nzoff*(zmaxico-zminico)/nzico
        yminoff=yminico+nyoff*(ymaxico-yminico)/nyico
        ymaxoff=ymaxico-nyoff*(ymaxico-yminico)/nyico
        kval=0
        ival=0
        kmax=1
        do k=2,length1
          if(word(k:k).eq.'(')kval=k+1
          if(word(k:k).eq.')')kmax=k-1
        enddo
        length=length1
        if(kval.gt.0)then
          length=kval-2
          read(word(kval:kmax),*)ival
        endif
        write(ifmt,fmt)'Graph ',word(1:length1)

                if(word(1:length).eq.'IniConEnergydensityEtaY')then

        call RootCanvas(1)
        call histo2BeginFigure(1)
        call histo2BeginPlot(1, 1)  
        call RootHisto(2,'energy density [GeV/fm^{3}]    (y=0);'
     *   ,nzico-2*nzoff,zminoff,zmaxoff
     *   ,nyico-2*nyoff,yminoff,ymaxoff)
        call histo2BeginSubPlot(word(1:length)
     *        ,'energy density [GeV/f$m^{3}$]    (y=0)' 
     *        ,'energy density [GeV/f$m^{3}$]    (y=0)'
     *        ,zminoff,zmaxoff,(zmaxico-zminico)/nzico
     *        ,yminoff,ymaxoff,(ymaxico-yminico)/nyico
     *        ,'$\eta_{s}$','y [fm]')
        ix=nxico/2+1+ival
        do iz=1+nzoff,nzico-nzoff
          zz=zminico+(float(iz)-0.5)*(zmaxico-zminico)/nzico
          do iy=1+nyoff,nyico-nyoff
            yy=yminico+(float(iy)-0.5)*(ymaxico-yminico)/nyico
            call RootFill(2,zz,yy,IcoE(ix,iy,iz))
            call histo2FillArray(2,zz,yy,IcoE(ix,iy,iz))
          enddo
        enddo
        call RootDraw('contz;','#eta_{s};','y [fm];')
        call histo2EndSubPlot()
        call histo2EndPlot()

                elseif(word(1:length).eq.'IniConEnergydensityXY')then

        iz=nzico/2+1+ival
        zz=zminico+(float(iz)-0.5)*(zmaxico-zminico)/nzico
        write(fmt,'(f7.2)')zz
        call RootCanvas(-3)
        call histo2BeginFigure(-3)
        call histo2BeginPlot(1, 1)
        call RootHisto(2
     *   ,'energy density [GeV/fm^{3}]    (#eta_{s}='//fmt//');'
     *   ,nxico,xminico,xmaxico
     *   ,nyico,yminico,ymaxico)
        call histo2BeginSubPlot(word(1:length) 
     *        ,'energy density [GeV/fm^{3}] ($\eta_{s}$ = '//fmt//')'
     *        ,'energy density [GeV/fm^{3}] ($\eta_{s}$ = '//fmt//')'
     *        ,xminico
     *        ,xmaxico
     *        ,(xmaxico-xminico)/nxico
     *        ,yminico
     *        ,ymaxico
     *        ,(ymaxico-yminico)/nyico
     *        ,'x [fm]','y [fm]')
        do ix=1,nxico
          xx=xminico+(float(ix)-0.5)*(xmaxico-xminico)/nxico
          do iy=1,nyico
            yy=yminico+(float(iy)-0.5)*(ymaxico-yminico)/nyico
            call RootFill(2,xx,yy,IcoE(ix,iy,iz))
            call histo2FillArray(2,xx,yy,IcoE(ix,iy,iz))
          enddo
        enddo
        call RootDraw('colz;','x [fm];','y [fm];')
        call histo2EndSubPlot()
        call histo2EndPlot()

                elseif(word(1:length).eq.'IniConXvelocityXY')then

        call RootCanvas(1)
        call histo2BeginFigure(-3)
        call histo2BeginPlot(1, 1)
         call histo2BeginSubPlot(word(1:length)
     *        ,'x velocity    ($\eta_{s}$ = 0)'
     *        ,'x velocity    ($\eta_{s}$ = 0)'
     *        ,yminico
     *        ,ymaxico
     *        ,(ymaxico-yminico)/nyico
     *        ,xminico
     *        ,xmaxico
     *        ,(xmaxico-xminico)/nxico
     *        ,'x [fm]','y [fm]')
        call RootHisto(2,
     *   'x velocity    (#eta_{s}=0);',nyico,yminico,ymaxico
     *                               ,nxico,xminico,xmaxico)
        iz=nzico/2+1+ival
        do iy=1,nyico
          yy=yminico+(float(iy)-0.5)*(ymaxico-yminico)/nyico
          do ix=1,nxico
            xx=xminico+(float(ix)-0.5)*(xmaxico-xminico)/nxico
            call RootFill(2,yy,xx,IcoV(1,ix,iy,iz))
            call histo2FillArray(2,yy,xx,IcoV(1,ix,iy,iz))
          enddo
        enddo
        call RootDraw('contz;','y [fm];','x [fm];')
        call histo2EndSubPlot()
        call histo2EndPlot()

                elseif(word(1:length).eq.'IniConXvelocityXEta')then

        call RootCanvas(1)
         call histo2BeginFigure(-3)
         call histo2BeginPlot(1, 1)
         call histo2BeginSubPlot(word(1:length)
     *        ,'x velocity    (y=0);'
     *        ,'x velocity    (y=0);'
     *        ,zminoff
     *        ,zmaxoff
     *        ,(zmaxico-zminico)/nzico
     *        ,xminico
     *        ,xmaxico
     *        ,(xmaxico-xminico)/nxico
     *        ,'x [fm]','y [fm]')
        call RootHisto(2,'x velocity    (y=0);',
     *   nzico-2*nzoff,zminoff,zmaxoff,nxico,xminico,xmaxico)
        iy=nyico/2+1+ival
        do iz=1+nzoff,nzico-nzoff
          zz=zminico+(float(iz)-0.5)*(zmaxico-zminico)/nzico
          do ix=1,nxico
            xx=xminico+(float(ix)-0.5)*(xmaxico-xminico)/nxico
            call RootFill(2,zz,xx,IcoV(1,ix,iy,iz))
            call histo2FillArray(2,zz,xx,IcoV(1,ix,iy,iz))
          enddo
        enddo
        call RootDraw('contz;','#eta_{s};','x [fm];')
        call histo2EndSubPlot()
        call histo2EndPlot()

                elseif(word(1:length).eq.'IniConRapidityXY')then

        call RootCanvas(1)
        call histo2BeginFigure(-3)
        call histo2BeginPlot(1, 1)
        call histo2BeginSubPlot(word(1:length)
     *        ,'rapidity    ($\eta_{s}$ = 0);'
     *        ,'rapidity    ($\eta_{s}$ = 0);'
     *        ,yminico
     *        ,ymaxico
     *        ,(ymaxico-yminico)/nyico
     *        ,xminico
     *        ,xmaxico
     *        ,(xmaxico-xminico)/nxico
     *        ,'x [fm]','y [fm]')
        call RootHisto(2,'rapidity    (#eta_{s}=0);',nyico,yminico
     *                               ,ymaxico,nxico,xminico,xmaxico)
        iz=nzico/2+1+ival
        zz=zminico+(float(iz)-0.5)*(zmaxico-zminico)/nzico
        do iy=1,nyico
          yy=yminico+(float(iy)-0.5)*(ymaxico-yminico)/nyico
          do ix=1,nxico
            xx=xminico+(float(ix)-0.5)*(xmaxico-xminico)/nxico
            vz=IcoV(3,ix,iy,iz)
            en=IcoE(ix,iy,iz)
            rap=0
            if(en.gt.0.)rap=0.5*alog((1+vz)/(1-vz))+zz
            call RootFill(2,yy,xx,rap)
            call histo2FillArray(2,yy,xx,rap)
          enddo
        enddo
        call RootDraw('contz;','y [fm];','x [fm];')
        call histo2EndSubPlot()
        call histo2EndPlot()

                elseif(word(1:length).eq.'IniConRapidityXEta')then

        call RootCanvas(1)
        call histo2BeginFigure(1)
        call histo2BeginPlot(1, 1)
        call histo2BeginSubPlot(word(1:length)
     *        ,'rapidity    (y=0)'
     *        ,'rapidity    (y=0)'
     *        ,zminoff
     *        ,zmaxoff
     *        ,(zmaxico-zminico)/nzico
     *        ,xminico
     *        ,xmaxico
     *        ,(xmaxico-xminico)/nxico
     *        ,'x [fm]','y [fm]')
        call RootHisto(2,'rapidity    (y=0);'
     *   ,nzico-2*nzoff,zminoff,zmaxoff,nxico,xminico,xmaxico)
        iy=nyico/2+1+ival
        do iz=1+nzoff,nzico-nzoff
          zz=zminico+(float(iz)-0.5)*(zmaxico-zminico)/nzico
          do ix=1,nxico
            xx=xminico+(float(ix)-0.5)*(xmaxico-xminico)/nxico
            vz=IcoV(3,ix,iy,iz)
            en=IcoE(ix,iy,iz)
            rap=0
            if(en.gt.0.)rap=0.5*alog((1+vz)/(1-vz))+zz
             call RootFill(2,zz,xx,rap)
             call histo2FillArray(2,zz,xx,rap)
          enddo
        enddo
        call RootDraw('contz;','#eta_{s};','x [fm];')
        call histo2EndSubPlot()
        call histo2EndPlot()

                elseif(word(1:length).eq.'IniConUDflavorXY')then

        call RootCanvas(1)
        call histo2BeginFigure(1)
        call histo2BeginPlot(1, 1)
        call histo2BeginSubPlot(word(1:length)
     *        ,'rapidity    ($\eta_{s}$ = 0)'
     *        ,'rapidity    ($\eta_{s}$ = 0)'
     *        ,yminico
     *        ,ymaxico
     *        ,(ymaxico-yminico)/nyico
     *        ,xminico
     *        ,xmaxico
     *        ,(xmaxico-xminico)/nxico
     *        ,'x [fm]','y [fm]')
        call RootHisto(2,'ud flavor    (#eta_{s}=0);',nyico,yminico
     *                               ,ymaxico,nxico,xminico,xmaxico)
        iz=nzico/2+1+ival
        do iy=1,nyico
          yy=yminico+(float(iy)-0.5)*(ymaxico-yminico)/nyico
          do ix=1,nxico
            xx=xminico+(float(ix)-0.5)*(xmaxico-xminico)/nxico
            call RootFill(2,yy,xx,IcoF(1,ix,iy,iz)+IcoF(2,ix,iy,iz))
            call histo2FillArray(2,
     *                      yy,xx,IcoF(1,ix,iy,iz)+IcoF(2,ix,iy,iz))
          enddo
        enddo
        call RootDraw('contz;','y [fm];','x [fm];')
        call histo2EndSubPlot()
        call histo2EndPlot()

                elseif(word(1:length).eq.'IniConUDflavorXEta')then

        call RootCanvas(1)
        call histo2BeginFigure(1)
        call histo2BeginPlot(1, 1)
        call histo2BeginSubPlot(word(1:length)
     *        ,'ud flavor    (y=0)'
     *        ,'ud flavor    (y=0)'
     *        ,zminoff
     *        ,zmaxoff
     *        ,(zmaxico-zminico)/nzico
     *        ,xminico
     *        ,xmaxico
     *        ,(xmaxico-xminico)/nxico
     *        ,'x [fm]','y [fm]')
        call RootHisto(2,'ud flavor    (y=0);'
     *   ,nzico-2*nzoff,zminoff,zmaxoff,nxico,xminico,xmaxico)
        iy=nyico/2+1+ival
        do iz=1+nzoff,nzico-nzoff
          zz=zminico+(float(iz)-0.5)*(zmaxico-zminico)/nzico
          do ix=1,nxico
            xx=xminico+(float(ix)-0.5)*(xmaxico-xminico)/nxico
            call RootFill(2,zz,xx,IcoF(1,ix,iy,iz)+IcoF(2,ix,iy,iz))
            call histo2FillArray(2,zz,xx,
     *                            IcoF(1,ix,iy,iz)+IcoF(2,ix,iy,iz))
          enddo
        enddo
        call RootDraw('contz;','#eta_{s};','x [fm];')
        call histo2EndSubPlot()
        call histo2EndPlot()

                elseif(word(1:length).eq.'IniConUflavorXY')then

        call RootCanvas(1)
        call histo2BeginFigure(1)
        call histo2BeginPlot(1, 1)
        call histo2BeginSubPlot(word(1:length)
     *        ,'u flavor    ($\eta_{s}$ = 0)'
     *        ,'u flavor    ($\eta_{s}$ = 0)'
     *        ,yminico
     *        ,ymaxico
     *        ,(ymaxico-yminico)/nyico
     *        ,xminico
     *        ,xmaxico
     *        ,(xmaxico-xminico)/nxico
     *        ,'x [fm]','y [fm]')

        call RootHisto(2,'u flavor    (#eta_{s}=0);',nyico,yminico
     *                               ,ymaxico,nxico,xminico,xmaxico)
        iz=nzico/2+1+ival
        do iy=1,nyico
          yy=yminico+(float(iy)-0.5)*(ymaxico-yminico)/nyico
          do ix=1,nxico
            xx=xminico+(float(ix)-0.5)*(xmaxico-xminico)/nxico
            call RootFill(2,yy,xx,IcoF(1,ix,iy,iz))
            call histo2FillArray(2,yy,xx,IcoF(1,ix,iy,iz))
          enddo
        enddo
        call RootDraw('contz;','y [fm];','x [fm];')
        call histo2EndSubPlot() 
        call histo2EndPlot() 

                elseif(word(1:length).eq.'IniConUflavorXEta')then

        call RootCanvas(1)
        call histo2BeginFigure(1)
        call histo2BeginPlot(1, 1)
        call histo2BeginSubPlot(word(1:length)
     *        ,'u flavor    (y=0)'
     *        ,'u flavor    (y=0)'
     *        ,zminoff
     *        ,zmaxoff
     *        ,(zmaxico-zminico)/nzico
     *        ,xminico
     *        ,xmaxico
     *        ,(xmaxico-xminico)/nxico
     *        ,'x [fm]','y [fm]')

        call RootHisto(2,'u flavor    (y=0);'
     *   ,nzico-2*nzoff,zminoff,zmaxoff,nxico,xminico,xmaxico)
        iy=nyico/2+1+ival
        do iz=1+nzoff,nzico-nzoff
          zz=zminico+(float(iz)-0.5)*(zmaxico-zminico)/nzico
          do ix=1,nxico
            xx=xminico+(float(ix)-0.5)*(xmaxico-xminico)/nxico
            call RootFill(2,zz,xx,IcoF(1,ix,iy,iz))
            call histo2FillArray(2,zz,xx,IcoF(1,ix,iy,iz))
          enddo
        enddo
        call RootDraw('contz;','#eta_{s};','x [fm];')
        call histo2EndSubPlot()
        call histo2EndPlot()

                elseif(word(1:length).eq.'IniConDflavorXY')then

        call RootCanvas(1)
        call histo2BeginFigure(1)
        call histo2BeginPlot(1, 1)
        call histo2BeginSubPlot(word(1:length)
     *        ,'d flavor    ($\eta_{s}$ = 0)'
     *        ,'d flavor    ($\eta_{s}$ = 0)'
     *        ,yminico
     *        ,ymaxico
     *        ,(ymaxico-yminico)/nyico
     *        ,xminico
     *        ,xmaxico
     *        ,(xmaxico-xminico)/nxico
     *        ,'x [fm]','y [fm]')

        call RootHisto(2,'d flavor    (#eta_{s}=0);',nyico,yminico
     *                               ,ymaxico,nxico,xminico,xmaxico)
        iz=nzico/2+1+ival
        do iy=1,nyico
          yy=yminico+(float(iy)-0.5)*(ymaxico-yminico)/nyico
          do ix=1,nxico
            xx=xminico+(float(ix)-0.5)*(xmaxico-xminico)/nxico
            call RootFill(2,yy,xx,IcoF(2,ix,iy,iz))
            call histo2FillArray(2,yy,xx,IcoF(2,ix,iy,iz))
          enddo
        enddo
        call RootDraw('contz;','y [fm];','x [fm];')
        call histo2EndSubPlot()
        call histo2EndPlot()

                elseif(word(1:length).eq.'IniConDflavorXEta')then

        call RootCanvas(1)
        call histo2BeginFigure(1)
        call histo2BeginPlot(1, 1)
        call histo2BeginSubPlot(word(1:length)
     *        ,'d flavor    (y=0)'
     *        ,'d flavor    (y=0)'
     *        ,zminoff
     *        ,zmaxoff
     *        ,(zmaxico-zminico)/nzico
     *        ,xminico
     *        ,xmaxico
     *        ,(xmaxico-xminico)/nxico
     *        ,'x [fm]','y [fm]')

        call RootHisto(2,'d flavor    (y=0);'
     *   ,nzico-2*nzoff,zminoff,zmaxoff,nxico,xminico,xmaxico)
        iy=nyico/2+1+ival
        do iz=1+nzoff,nzico-nzoff
          zz=zminico+(float(iz)-0.5)*(zmaxico-zminico)/nzico
          do ix=1,nxico
            xx=xminico+(float(ix)-0.5)*(xmaxico-xminico)/nxico
            call RootFill(2,zz,xx,IcoF(2,ix,iy,iz))
            call histo2FillArray(2,zz,xx,IcoF(2,ix,iy,iz))
          enddo
        enddo
        call RootDraw('contz;','#eta_{s};','x [fm];')
        call histo2EndSubPlot()
        call histo2EndPlot()

                elseif(word(1:length).eq.'IniConSflavorXY')then

        call RootCanvas(1)
        call histo2BeginFigure(1)
        call histo2BeginPlot(1, 1) 
        call histo2BeginSubPlot(word(1:length)
     *        ,'s flavor    ($\eta_{s}$ = 0)'
     *        ,'s flavor    ($\eta_{s}$ = 0)'
     *        ,yminico
     *        ,ymaxico
     *        ,(ymaxico-yminico)/nyico
     *        ,xminico
     *        ,xmaxico
     *        ,(xmaxico-xminico)/nxico
     *        ,'x [fm]','y [fm]')

        call RootHisto(2,'s flavor    (#eta_{s}=0);',nyico,yminico
     *                               ,ymaxico,nxico,xminico,xmaxico)
        iz=nzico/2+1+ival
        do iy=1,nyico
          yy=yminico+(float(iy)-0.5)*(ymaxico-yminico)/nyico
          do ix=1,nxico
            xx=xminico+(float(ix)-0.5)*(xmaxico-xminico)/nxico
            call RootFill(2,yy,xx,IcoF(3,ix,iy,iz))
            call histo2FillArray(2,yy,xx,IcoF(3,ix,iy,iz))
          enddo
        enddo
        call RootDraw('contz;','y [fm];','x [fm];')
        call histo2EndSubPlot()
        call histo2EndPlot()

                elseif(word(1:length).eq.'IniConSflavorXEta')then

        call RootCanvas(1)
        call histo2BeginFigure(1)
        call histo2BeginPlot(1, 1)
        call histo2BeginSubPlot(word(1:length)
     *        ,'s flavor    (y=0)'
     *        ,'s flavor    (y=0)'
     *        ,zminoff
     *        ,zmaxoff
     *        ,(zmaxico-zminico)/nzico
     *        ,xminico
     *        ,xmaxico
     *        ,(xmaxico-xminico)/nxico
     *        ,'x [fm]','y [fm]')

        call RootHisto(2,'s flavor    (y=0);'
     *   ,nzico-2*nzoff,zminoff,zmaxoff,nxico,xminico,xmaxico)
        iy=nyico/2+1+ival
        do iz=1+nzoff,nzico-nzoff
          zz=zminico+(float(iz)-0.5)*(zmaxico-zminico)/nzico
          do ix=1,nxico
            xx=xminico+(float(ix)-0.5)*(xmaxico-xminico)/nxico
            call RootFill(2,zz,xx,IcoF(3,ix,iy,iz))
            call histo2FillArray(2,zz,xx,IcoF(3,ix,iy,iz))
          enddo
        enddo
        call RootDraw('contz;','#eta_{s};','x [fm];')
        call histo2EndSubPlot() 
        call histo2EndPlot()

      endif

      end

c#######################################################################
c#######################################################################
c#######################################################################
c#######################################################################
c#######################################################################
c#######################################################################
c#######################################################################

c-----------------------------------------------------------------------
      subroutine xSpaceTime(kkk)
c-----------------------------------------------------------------------
#include "aaa.h"
#include "ico.h"
#include "ems.h"
      common/cmean/xmean,ymean
      common/cshft/xshft,yshft
      common/ciprint2d/iprint2d
      call srand(86456)
      bx=cos(phievt)*bimevt
      by=sin(phievt)*bimevt
      xsh=bx/2
      ysh=by/2
      phinll=phievt+ranphi
      cc=cos(phinll)
      ss=sin(phinll)
      xshft=0
      yshft=0
      if(maproj*matarg.eq.1.and.nprt(1).lt.15)then
      !print*,'number of Pomerons:',nprt(1)
      endif
      if(kkk.eq.1)then
         call xCoreCorona(1,0,0)
      elseif(kkk.eq.2)then
         call xCoreCorona(2,0,0)
      elseif(kkk.eq.3)then
         call xCoreCorona(3,0,0)
      elseif(kkk.eq.31)then
        iprint2d=1
        x= xmean + xsh   !   undo         target center
        y= ymean + ysh   ! mean shift      -> (0,0)
        xshft= cc*x+ss*y
        yshft=-ss*x+cc*y 
        call xCoreCorona(3,0,0)
        !call xBinary(2.)
      elseif(kkk.eq.4)then
         call xCoreCorona(3,1,0)
      elseif(kkk.eq.5)then
         call xCoreCorona(5,1,0)
      elseif(kkk.eq.6)then
         call xCoreCorona(6,0,0)
      elseif(kkk.eq.8)then
         stop'ERROR 31122013 ######################################'
           !call xSpaceTime(8) obsolete, use call xSpaceTime(31)
      endif
      end

c-----------------------------------------------------------------------
      subroutine xBinary(xmax)
c-----------------------------------------------------------------------
#include "aaa.h"
#include "ico.h"
#include "ems.h"
      write(ifhi,'(a)')       '!----------------------------------'
      write(ifhi,'(a)')       '!   coord            '
      write(ifhi,'(a)')       '!----------------------------------'
      write(ifhi,'(a)')      'openhisto name coord'
      write(ifhi,'(a,2e11.3)')'xrange',-xmax*1.35,xmax*1.35
      write(ifhi,'(a,2e11.3)')'yrange',-xmax,xmax
      write(ifhi,'(a)')       'htyp pbc xmod lin ymod lin'
      write(ifhi,'(a)')       'array 2'
      do k=1,koll
        write(ifhi,'(2e11.3)') coord(1,k) , coord(2,k)
      enddo
      write(ifhi,'(a)') '  endarray closehisto plot 0-'
      write(ifhi,'(a)')       '!----------------------------------'
      write(ifhi,'(a)')       '!   coordpr          '
      write(ifhi,'(a)')       '!----------------------------------'
      write(ifhi,'(a)')       'openhisto name binary'
      write(ifhi,'(a)')       'htyp pxk xmod lin ymod lin'
      write(ifhi,'(a,2e11.3)')'xrange',-xmax,xmax
      write(ifhi,'(a,2e11.3)')'yrange',-xmax,xmax
      write(ifhi,'(a)')       'array 2'
      do k=1,koll
        do n=1,nprmx(k)
           write(ifhi,'(2e11.3)')coordpr(1,n,k),coordpr(2,n,k)
        enddo
      enddo
      write(ifhi,'(a)')    '  endarray closehisto plot 0'
      end

c-----------------------------------------------------------------------
      subroutine xCoreCorona(kkk,iii,jjj)
c-----------------------------------------------------------------------
c  red circle   = core
c  green circle = corona
c  blue circle  = Pomerons
c  black star   = collison 
c-----------------------------------------------------------------------
c     space-time evolution of core and corona
c
c     cluster ............   ist=11  ity=60
c     core particles .....   ist=0   ity=60
c     corona particles ...   ist=0   ity/=60
c
c    iii=1: plot also possible Pomerons
c    jjj>0: multiplicity trigger (useful for pp)
c------------------------------------------------xSpace-----------------------
#include "aaa.h"
#include "ico.h"
#include "ems.h"
      common/cxyzt/xptl(mxptl+29),yptl(mxptl+29),zptl(mxptl+29)
     * ,tptl(mxptl+29),optl(mxptl+29),uptl(mxptl+29),sptl(mxptl+29)
     *,rptl(mxptl+29,3)
      common/cdelzet/delzet,delsce
      common/cdelcore/delcore,egycore
      parameter (myy=12,mrr=21)
      real yy(myy),yy1(myy),rr(mrr),rr1(mrr),rrr(mrr,mrr)
      common/cranphi/ranphi
      character ch1s*3,ch2s*3,ch1r*3,ch2r*3,ch1*3,ch2*3,chh*3
      common/cshft/xshft,yshft
      common/ciprint2d/iprint2d

      call utpri('xcorec',ish,ishini,4)

      bx=cos(phievt)*bimevt
      by=sin(phievt)*bimevt
      xsh=bx/2
      ysh=by/2
      phinll=phievt+ranphi
      cc=cos(phinll)
      ss=sin(phinll)
      if(ish.ge.5)then
        do i=1,nptl
          if(istptl(i).eq.5.or.istptl(i).eq.7)then
            amt=pptl(5,i)**2+pptl(1,i)**2+pptl(2,i)**2
            rap=1000
            if(amt.gt.0..and.pptl(4,i).gt.0.)then
              amt=sqrt(amt)
              rap=sign(1.,pptl(3,i))*
     .        alog((pptl(4,i)+abs(pptl(3,i)))/amt)
            endif
            write(ifch,*)'xSpaceTime:',i,dezptl(i),rap,pptl(3,i)
          endif
        enddo
      endif
      phinll=phievt+ranphi
      if(phinll.lt.-pi)phinll=phinll+2*pi
      if(phinll.gt.pi)phinll=phinll-2*pi
      rapmax=6
      radmax=10
      r1=0.00001
      if(maproj.gt.1)r1=radnuc(maproj)
      r2=0.00001
      if(matarg.gt.1)r2=radnuc(matarg)
      a=0
      ch1s='   '
      ch2s='   '
      ch1r='   '
      ch2r='   '
      fac=0.
      if(kkk.eq.1.or.kkk.eq.5)Then
       a=2
       ch1s='prk' !core from string
       ch2s='pgk' !coro from string
       ch1r='prk' !core from remnant
       ch2r='pgk' !coro from remnant
       fac=0.37
      elseif(kkk.eq.2)Then
       a=8.2
       ch1s='prk' !core from string
       ch2s='pgk' !coro from string
       ch1r='prk' !core from remnant
       ch2r='pgk' !coro from remnant
       fac=0.5
      elseif(kkk.eq.3)Then
       a=2.2
       ch1s='prk' !core from string
       ch2s='pgk' !coro from string
       ch1r='prk' !core from remnant
       ch2r='pgk' !coro from remnant
       fac=0.5
      endif
      radmax=a
      n1=koievt
      n2=0
      do k=1,koll
       if(itpr(k).gt.0)n2=n2+1
      enddo
      n3=nglevt
      npoms=0
      do i=1,nptl
       if(istptl(i).eq.30.or.istptl(i).eq.31)npoms=npoms+1
      enddo

      if(jjj.gt.0)then
      multy1=0
       do i=maproj+matarg+1,nptl
        if(istptl(i).eq.0)then
          amt=pptl(5,i)**2+pptl(1,i)**2+pptl(2,i)**2
          rap=1000
          if(amt.gt.0..and.pptl(4,i).gt.0.)then
            amt=sqrt(amt)
            rap=sign(1.,pptl(3,i))*alog((pptl(4,i)+abs(pptl(3,i)))/amt)
          endif
          ch=0
          if(abs(idptl(i)).ge.100.and.abs(idptl(i)).lt.10000)then
            call idchrg( 7 ,idptl(i),ch)
            if(abs(ch).gt.0.1.and.abs(rap).le.1.)multy1=multy1+1
          endif
         endif
        enddo
        ih1=jjj/100
        ih2=mod(jjj,100)
        if(0.5*multy1.lt.ih1.or.0.5*multy1.gt.ih2)goto 999
      endif

      xmax=a*1.35 !1+int(a*1.5)
      xunt=-xmax
      xob=+xmax
      yunt=-a
      yob=+a

      do loopeta=4,10 !+++++++++++++++++++++++++++++++++

      loo1=loopeta/10
      loo2=mod(loopeta,10)
      etamin=-6.5+(loopeta-1.)
      etamax=etamin+1
      eta=(etamin+etamax)/2

      do loor=1,2 !+++++++++++++++++++++++++++++++++

      if(loor.eq.1)then
        ch1=ch1s
        ch2=ch2s
      endif
      if(loor.eq.2)then
        ch1=ch1r
        ch2=ch2r
      endif

      if(kkk.ge.1.and.kkk.le.5)then !~~~~~~~~~~~~~~~~~~

      write(ifhi,'(a)')       '!---------------------------------'
      write(ifhi,'(a)')       '!   core particles                '
      write(ifhi,'(a)')       '!---------------------------------'
      write(ifhi,'(a)')       '!newpage'
      write(ifhi,'(a,3i1)')      'openhisto name ste',loo1,loo2,loor
      write(ifhi,'(a)')       'htyp '//ch1//' xmod lin ymod lin'
      write(ifhi,'(a,2e11.3)')'xrange',xunt,xob
      write(ifhi,'(a,2e11.3)')'yrange',yunt,yob
      write(ifhi,'(a)')    'txt "xaxis x (fm)"'
      write(ifhi,'(a)')    'txt "yaxis y (fm)"'

      if(kkk.le.3.or.kkk.eq.5)then !~~~~~~~~~~~~~~~~~~

       write(ifhi,'(a)')       'text 0.05 0.94 "core-"'
       write(ifhi,'(a)')       'text 0.02 0.87 "corona"'
       write(ifhi,'(a,f4.1,a)')'text 0.05 0.07 "',bimevt,'fm"'
       write(ifhi,'(a,i3,a)')  'text 0.65 0.07 "',npoms
     .  ,' Pomerons"'
       write(ifhi,'(a,f5.2,a)')'text 0.70 0.93 "[c] = ',eta,'"'
       write(ifhi,'(a)')       'array 2'
       ncore=0
       iorp=0
       do i=1,nptl
        if(dezptl(i).gt.etamin.and.dezptl(i).le.etamax
     .   .and.(  (loor.eq.1.and.ityptl(i).ge.20.and.ityptl(i).le.39)
     .         .or.(loor.eq.2.and.ityptl(i).ge.40)  )
     .   .and.(istptl(i).eq.5.or.istptl(i).eq.7)
     .    .and.i.ge.minfra.and.i.le.maxfra)then
         ior=iorptl(i)
         ! print*,'55555',i,idptl(i) , istptl(i),' '
         !.  , xptl(i),yptl(i),' ',dezptl(i),ior,iorp,delcore
         ra=0
         if(ior.eq.iorp)ra=0.05
         x=xptl(i)
         y=yptl(i)
         write(ifhi,'(2e13.5)')
     .           cc*x+ss*y +xshft  +ra*cos(2*3.14159*rand())
     .      ,   -ss*x+cc*y +yshft  +ra*sin(2*3.14159*rand())
         ncore=ncore+1
         iorp=ior
        endif
       enddo
      elseif(kkk.eq.4)then !~~~~~~~~~~~~~~~~~
       write(ifhi,'(a)')       'array 2'
       ncore=0
       iz=nzico/2+1
       do ix=1,nxico
        do iy=1,nyico
         x=xminico+(ix-0.5)*(xmaxico-xminico)/nxico
         y=yminico+(iy-0.5)*(ymaxico-yminico)/nyico
         !print*,x,y,IcoE(ix,iy,iz)
         if(IcoE(ix,iy,iz).gt.2)
     .    write(ifhi,'(2e11.3)')   x   +xshft
     .      ,                      y   +yshft
        enddo
       enddo

      endif !~~~~~~~~~~~~~~~~~

      write(ifhi,'(a)')    'endarray closehisto plot 0-'
      write(ifhi,'(a)')       '!----------------------------------'
      write(ifhi,'(a)')       '!   corona particles               '
      write(ifhi,'(a)')       '!----------------------------------'
      write(ifhi,'(a,3i1)')      'openhisto name sta',loo1,loo2,loor
      write(ifhi,'(a)')       'htyp '//ch2//' xmod lin ymod lin'
      write(ifhi,'(a)')       'array 2'
      ncorona=0
      iorp=0
      do i=1,nptl
       if(dezptl(i).gt.etamin.and.dezptl(i).le.etamax
     .   .and.(  (loor.eq.1.and.ityptl(i).ge.20.and.ityptl(i).le.39)
     .         .or.(loor.eq.2.and.ityptl(i).ge.40)  )
     .  .and.(istptl(i).eq.0.or.istptl(i).eq.1)
     .    .and.i.ge.minfra.and.i.le.maxfra)then
         ior=iorptl(i)
         ra=0
         if(ior.eq.iorp)ra=0.05
         x=xptl(i)
         y=yptl(i)
         write(ifhi,'(2e13.5)')
     .       cc*x+ss*y  +xshft +ra*cos(2*3.14159*rand())
     .    , -ss*x+cc*y  +yshft +ra*sin(2*3.14159*rand())
         ncorona=ncorona+1
         iorp=ior
       endif
      enddo
      write(ifhi,'(a)')    '  endarray'
      write(ifhi,'(a)')    'closehisto'
      !print*,'b=',bimevt,'   ncorona:ncore =  ',ncorona,':',ncore
      if(loor.eq.1)write(ifhi,'(a)')  'plot 0-'

      endif !~~~~~~~~~~

      enddo !+++++++++++++++++++++++++++++++++++++++++++

      if(iii.eq.1)then

        write(ifhi,'(a)')       'plot 0-'
        write(ifhi,'(a)')       '!----------------------------------'
        write(ifhi,'(a)')       '!   possible Pomerons            '
        write(ifhi,'(a)')       '!----------------------------------'
        write(ifhi,'(a)')      'openhisto name coord'
        write(ifhi,'(a,2e11.3)')'xrange',xunt,xob
        write(ifhi,'(a,2e11.3)')'yrange',yunt,yob
        write(ifhi,'(a)')       'htyp pxk xmod lin ymod lin'
        write(ifhi,'(a)')       'array 2'
        do k=1,koll
          do n=1,nprmx(k)
c            if(idpr(n,k).gt.0)then
              x=coordpr(1,n,k)
              y=coordpr(2,n,k)
              write(ifhi,'(2e11.3)')
     *           cc*x+ss*y  +xshft
     *         ,-ss*x+cc*y  +yshft
c            endif
          enddo
        enddo
        write(ifhi,'(a)') '  endarray closehisto'
        write(ifhi,'(a)')       'plot 0-'
        write(ifhi,'(a)')       '!----------------------------------'
        write(ifhi,'(a)')       '!   possible collisions            '
        write(ifhi,'(a)')       '!----------------------------------'
        write(ifhi,'(a)')      'openhisto name coord'
        write(ifhi,'(a,2e11.3)')'xrange',xunt,xob
        write(ifhi,'(a,2e11.3)')'yrange',yunt,yob
        write(ifhi,'(a)')       'htyp pxs xmod lin ymod lin'
        write(ifhi,'(a)')       'array 2'
        do k=1,koll
              x=coord(1,k)
              y=coord(2,k)
              write(ifhi,'(2e11.3)')
     *          cc*x+ss*y  +xshft
     *        ,-ss*x+cc*y  +yshft
        enddo
        write(ifhi,'(a)') '  endarray closehisto'

c        write(ifhi,'(a)') ' plot 0-'
c        write(ifhi,'(a)')       '!----------------------------------'
c        write(ifhi,'(a)')       '!   string positions             '
c        write(ifhi,'(a)')       '!----------------------------------'
c        write(ifhi,'(a)')   'openhisto name coord htyp pot'
c        write(ifhi,'(a)')   'array 2'
c        nk1=1
c        do while(nk1.le.nptl)
c          if(istptl(nk1).eq.20.or.istptl(nk1).eq.21)then
c            nk2=nk1+1
c            do while(idptl(nk2).eq.9)
c              nk2=nk2+1
c            enddo
c            x=xorptl(1,nk1)
c            y=xorptl(2,nk1)
c            write(ifhi,'(2e11.3)')
c     *     cc*x+ss*y  +xshft
c     *   ,-ss*x+cc*y  +yshft
c            nk1=nk2+1
c          else
c            nk1=nk1+1
c          endif
c        enddo
c        write(ifhi,'(a)') '  endarray closehisto'

      endif

      if(kkk.ge.1.and.kkk.le.5.and.koll.lt.100)then

        write(ifhi,'(a)')        'plot 0-'
        write(ifhi,'(a)')       '!----------------------------------'
        write(ifhi,'(a)')       '!   cut Pomerons             '
        write(ifhi,'(a)')       '!----------------------------------'
        write(ifhi,'(a)')      'openhisto name coo'
        write(ifhi,'(a)')       'htyp pbc xmod lin ymod lin'
        write(ifhi,'(a,2e11.3)')'xrange',xunt,xob
        write(ifhi,'(a,2e11.3)')'yrange',yunt,yob
        write(ifhi,'(a)')       'array 2'
        do k=1,koll
          do n=1,nprmx(k)
            if(idpr(n,k).gt.0)then
              x=coordpr(1,n,k)
              y=coordpr(2,n,k)
              write(ifhi,'(2e11.3)')
     *        cc*x+ss*y  +xshft
     * ,     -ss*x+cc*y  +yshft
            endif
          enddo
        enddo
        write(ifhi,'(a)')     'endarray closehisto'
        write(ifhi,'(a)')     'plot 0-'
        write(ifhi,'(a)')       '!----------------------------------'
        write(ifhi,'(a)')       '!   collisions            '
        write(ifhi,'(a)')       '!----------------------------------'
        write(ifhi,'(a)')      'openhisto name coord'
        write(ifhi,'(a,2e11.3)')'xrange',xunt,xob
        write(ifhi,'(a,2e11.3)')'yrange',yunt,yob
        write(ifhi,'(a)')       'htyp pks xmod lin ymod lin'
        write(ifhi,'(a)')       'array 2'
        do k=1,koll
          if(nprt(k).gt.0)then
            if(loopeta.eq.7)then 
              i=iproj(k)
              j=itarg(k)
              write(ifmt,'(a,2i4,$)')'Dipole',i,j
              x=diproj(i,1)*0.5
              y=diproj(i,2)*0.5 
              write(ifmt,'(2f8.2,$)')           
     *          cc*x+ss*y+xshft
     *        ,-ss*x+cc*y+yshft
              x=ditarg(j,1)*0.5
              y=ditarg(j,2)*0.5
              write(ifmt,'(3x,2f8.2)')           
     *          cc*x+ss*y+xshft
     *        ,-ss*x+cc*y+yshft
            endif
            x=coord(1,k)
            y=coord(2,k)
            write(ifhi,'(2e11.3)')
     *          cc*x+ss*y+xshft
     *        ,-ss*x+cc*y+yshft
          endif
        enddo
        write(ifhi,'(a)') '  endarray closehisto'

      endif

      if(kkk.ge.2.and.r1.ne.0.0)then
        ncirc=max(1.,300*r1)
        chh='pyc'
        if(ncirc.gt.300)chh='pyl'
        write(ifhi,'(a)')    'plot 0-'
        write(ifhi,'(a)')       '!----------------------------------'
        write(ifhi,'(a)')       '!   hard spheres             '
        write(ifhi,'(a)')       '!----------------------------------'
        write(ifhi,'(a,2i1)')   'openhisto name stc',loo1,loo2
        write(ifhi,'(a)')   'htyp '//chh
        write(ifhi,'(a)')   'array 2'
        do j=1,ncirc
         ang=2*pi*j/ncirc
         x=bx/2
         y=by/2
         a= cc*x+ss*y
         b=-ss*x+cc*y
         x=r1*cos(ang) +a +xshft
         y=r1*sin(ang) +b +yshft
         if(x.gt.xunt.and.x.lt.xob
     .   .and.y.gt.yunt.and.y.lt.yob)
     .   write(ifhi,'(2e13.5)')x,y
        enddo
        write(ifhi,'(a)')    '  endarray'
        write(ifhi,'(a)')    'closehisto'
      endif
      if(kkk.ge.2.and.r2.ne.0.0)then
        ncirc=max(1.,300*r2)
        chh='pyc'
        if(ncirc.gt.300)chh='pyl'
        write(ifhi,'(a)')    'plot 0-'
        write(ifhi,'(a,2i1)')   'openhisto name std' ,loo1,loo2
        write(ifhi,'(a)') 'htyp '//chh
        write(ifhi,'(a)')   'array 2'
        do j=1,ncirc
         ang=2*pi*j/ncirc
         x=-bx/2
         y=-by/2
         a= cc*x+ss*y
         b=-ss*x+cc*y
         x=r2*cos(ang) +a  +xshft
         y=r2*sin(ang) +b  +yshft
         if(x.gt.xunt.and.x.lt.xob
     .   .and.y.gt.yunt.and.y.lt.yob)
     .    write(ifhi,'(2e13.5)')x,y
        enddo
        write(ifhi,'(a)')    '  endarray'
        write(ifhi,'(a)')    'closehisto'
      endif

      write(ifhi,'(a)')    'plot 0'

      enddo !+++++++++++++++++++++++++++++++++++++++++++

   !....................................................................

      if(kkk.eq.3)goto 999

      delrad=2*radmax/float(mrr)
      delrap=2*rapmax/float(myy)

   !....................................................................
   !....................................................................

      do m=1,mrr
      do n=1,mrr
        rrr(m,n)=0
      enddo
      enddo
      do n=1,nptl
       if(n.ge.minfra.and.n.le.maxfra
     . .and.(istptl(n).eq.5.or.istptl(n).eq.7))then
        routp=-sin(phinll)*(xptl(n)+xshft)+cos(phinll)*(yptl(n)+yshft)
        rinp = cos(phinll)*(xptl(n)+xshft)+sin(phinll)*(yptl(n)+yshft)
        rapx=dezptl(n)
        if(abs(rapx).le.1)then
          eco=0
          amt=sqrt(pptl(5,n)**2+pptl(1,n)**2+pptl(2,n)**2)
          i=(rinp+radmax)/delrad+1
          j=(routp+radmax)/delrad+1
          if(i.ge.1.and.i.le.mrr.and.j.ge.1.and.j.le.mrr)
     .      rrr(i,j)=rrr(i,j)+amt
        endif
       endif
      enddo
      if(iprint2d.eq.1)then
        write(ifch,*)'=== print2d (ico.f) === dE* / dx / dy / deta'
        write(ifch,*)'Rotation angle: ',phinll/2/pi*360,' deg'
     . ,'       Shift: ',nint(xshft*100),nint(yshft*100),'  [fm/100]'
        write(ifch,'(7x,21i4)')
     .  (nint((-radmax+(i-0.5)*delrad)*100),i=1,mrr)
        do j=mrr,1,-1
          y=-radmax+(j-0.5)*delrad
          write(ifch,'(i5,2x,21i4)')
     .    nint(y*100),(nint(rrr(i,j)/2/delrad**2),i=1,mrr)
        enddo
      endif

   !....................................................................
   !....................................................................

      do m=1,mrr
        rr(m)=0
      enddo
      do n=1,nptl
        if(n.ge.minfra.and.n.le.maxfra
     .  .and.(istptl(n).eq.5.or.istptl(n).eq.7))then
        routp=-sin(phinll)*(xptl(n)+xshft)+cos(phinll)*(yptl(n)+yshft)
        rinp = cos(phinll)*(xptl(n)+xshft)+sin(phinll)*(yptl(n)+yshft)
        rapx=dezptl(n)
        if(abs(rapx).le.1.and.abs(routp).le.delrad/2)then
          eco=0
          amt=pptl(5,n)**2+pptl(1,n)**2+pptl(2,n)**2
          if(amt.gt.0..and.pptl(4,n)+abs(pptl(3,n)).gt.0.d0)then
            amt=sqrt(amt)
            rap=sign(1.,pptl(3,n))*alog((pptl(4,n)+abs(pptl(3,n)))/amt)
            eco=amt   !*cosh(rap-rapx)
          endif
          m=(rinp+radmax)/delrad+1
          if(m.ge.1.and.m.le.mrr)rr(m)=rr(m)+eco
        endif
        endif
      enddo
      write(ifhi,'(a)')'!---------------------------------------------'
      write(ifhi,'(a)')'!    core segment energy per d[c]dxdy vs x    '
      write(ifhi,'(a)')'!   (same as histogram rinp eco... in optns)  '
      write(ifhi,'(a)')'!---------------------------------------------'
      write(ifhi,'(a)')       '!newpage'
      write(ifhi,'(a)')    'openhisto name rapx'
      write(ifhi,'(a)')       'htyp lin xmod lin ymod lin'
      write(ifhi,'(a,2f7.3)') 'xrange ',-radmax,radmax
      write(ifhi,'(a)')        'yrange 0 auto '
      write(ifhi,'(a)')  'txt "title initial energy          "'
      write(ifhi,'(a,f4.1,a)')'text 0.05 0.70 "[c]=0"'
      write(ifhi,'(a,f4.1,a)')'text 0.05 0.60 "y=0"'
      write(ifhi,'(a,f4.1,a)')'text 0.65 0.9 "',bimevt,'fm"'
      write(ifhi,'(a)')  'txt  "xaxis x (fm)"'
      write(ifhi,'(a)')  'txt  "yaxis dE/d[c]dxdy "'
      write(ifhi,'(a)')       'array 2'
      do m=1,mrr
        write(ifhi,'(2e11.3)')-radmax+(m-0.5)*delrad, rr(m)/2/delrad**2
        rr1(m)=rr(m)
      enddo
      write(ifhi,'(a)')  '  endarray closehisto plot 0-'
   !....................................................................
      do m=1,mrr
        rr(m)=0
      enddo
      do n=1,nptl
        if(n.ge.minfra.and.n.le.maxfra.and.istptl(n)/2.eq.0)then
        routp=-sin(phinll)*(xptl(n)+xshft)+cos(phinll)*(yptl(n)+yshft)
        rinp = cos(phinll)*(xptl(n)+xshft)+sin(phinll)*(yptl(n)+yshft)
        rapx=dezptl(n)
        if(abs(rapx).le.1.and.abs(routp).le.delrad/2)then
          eco=0
          amt=pptl(5,n)**2+pptl(1,n)**2+pptl(2,n)**2
          if(amt.gt.0..and.pptl(4,n)+abs(pptl(3,n)).gt.0.d0)then
            amt=sqrt(amt)
            rap=sign(1.,pptl(3,n))*alog((pptl(4,n)+abs(pptl(3,n)))/amt)
            eco=amt   !*cosh(rap-rapx)
          endif
          m=(rinp+radmax)/delrad+1
          if(m.ge.1.and.m.le.mrr)rr(m)=rr(m)+eco
        endif
        endif
      enddo
      write(ifhi,'(a)')'!---------------------------------------------'
      write(ifhi,'(a)')'!    corona segment energy per d[c]dxdy vs x  '
      write(ifhi,'(a)')'!   (same as histogram rinp eco... in optns)  '
      write(ifhi,'(a)')'!---------------------------------------------'
      write(ifhi,'(a)')    'openhisto name rapx'
      write(ifhi,'(a)')       'htyp lin xmod lin ymod lin'
      write(ifhi,'(a)')       'array 2'
      do m=1,mrr
        write(ifhi,'(2e11.3)')-radmax+(m-0.5)*delrad
     .    , (rr1(m)+rr(m))/2/delrad**2
      enddo
      write(ifhi,'(a)')  '  endarray closehisto plot 0-'
      write(ifhi,'(a)')  'plot 0'

   !....................................................................
   !....................................................................

      do m=1,mrr
        rr(m)=0
      enddo
      do n=1,nptl
        if(n.ge.minfra.and.n.le.maxfra
     .   .and.(istptl(n).eq.5.or.istptl(n).eq.7))then
        routp=-sin(phinll)*(xptl(n)+xshft)+cos(phinll)*(yptl(n)+yshft)
        rinp = cos(phinll)*(xptl(n)+xshft)+sin(phinll)*(yptl(n)+yshft)
        rapx=dezptl(n)
        if(abs(rapx).le.1.and.abs(rinp).le.delrad/2)then
          eco=0
          amt=pptl(5,n)**2+pptl(1,n)**2+pptl(2,n)**2
          if(amt.gt.0..and.pptl(4,n)+abs(pptl(3,n)).gt.0.d0)then
            amt=sqrt(amt)
            rap=sign(1.,pptl(3,n))*alog((pptl(4,n)+abs(pptl(3,n)))/amt)
            eco=amt   !*cosh(rap-rapx)
          endif
          m=(routp+radmax)/delrad+1
          if(m.ge.1.and.m.le.mrr)rr(m)=rr(m)+eco
        endif
        endif
      enddo
      write(ifhi,'(a)')'!---------------------------------------------'
      write(ifhi,'(a)')'!    core segment energy per d[c]dxdy vs y    '
      write(ifhi,'(a)')'!   (same as histogram routp eco... in optns)  '
      write(ifhi,'(a)')'!---------------------------------------------'
      write(ifhi,'(a)')       '!newpage'
      write(ifhi,'(a)')    'openhisto name rout'
      write(ifhi,'(a)')       'htyp lin xmod lin ymod lin'
      write(ifhi,'(a,2f7.3)') 'xrange ',-radmax,radmax
      write(ifhi,'(a)')        'yrange 0 auto '
      write(ifhi,'(a)')  'txt "title initial energy          "'
      write(ifhi,'(a,f4.1,a)')'text 0.05 0.70 "[c]=0"'
      write(ifhi,'(a,f4.1,a)')'text 0.05 0.60 "x=0"'
      write(ifhi,'(a,f4.1,a)')'text 0.65 0.9 "',bimevt,'fm"'
      write(ifhi,'(a)')  'txt  "xaxis y (fm)"'
      write(ifhi,'(a)')  'txt  "yaxis dE/d[c]dxdy "'
      write(ifhi,'(a)')       'array 2'
      do m=1,mrr
        write(ifhi,'(2e11.3)')-radmax+(m-0.5)*delrad, rr(m)/2/delrad**2
        rr1(m)=rr(m)
      enddo
      write(ifhi,'(a)')  '  endarray closehisto plot 0-'
   !....................................................................
      do m=1,mrr
        rr(m)=0
      enddo
      do n=1,nptl
        if(n.ge.minfra.and.n.le.maxfra.and.istptl(n)/2.eq.0)then
        routp=-sin(phinll)*(xptl(n)+xshft)+cos(phinll)*(yptl(n)+yshft)
        rinp = cos(phinll)*(xptl(n)+xshft)+sin(phinll)*(yptl(n)+yshft)
        rapx=dezptl(n)
        if(abs(rapx).le.1.and.abs(rinp).le.delrad/2)then
          eco=0
          amt=pptl(5,n)**2+pptl(1,n)**2+pptl(2,n)**2
          if(amt.gt.0..and.pptl(4,n)+abs(pptl(3,n)).gt.0.d0)then
            amt=sqrt(amt)
            rap=sign(1.,pptl(3,n))*alog((pptl(4,n)+abs(pptl(3,n)))/amt)
            eco=amt   !*cosh(rap-rapx)
          endif
          m=(routp+radmax)/delrad+1
          if(m.ge.1.and.m.le.mrr)rr(m)=rr(m)+eco
        endif
        endif
      enddo
      write(ifhi,'(a)')'!---------------------------------------------'
      write(ifhi,'(a)')'!    corona segment energy per d[c]dxdy vs y  '
      write(ifhi,'(a)')'!   (same as histogram routp eco... in optns)  '
      write(ifhi,'(a)')'!---------------------------------------------'
      write(ifhi,'(a)')    'openhisto name rout'
      write(ifhi,'(a)')       'htyp lin xmod lin ymod lin'
      write(ifhi,'(a)')       'array 2'
      do m=1,mrr
        write(ifhi,'(2e11.3)')-radmax+(m-0.5)*delrad
     .   , (rr1(m)+rr(m))/2/delrad**2
      enddo
      write(ifhi,'(a)')  '  endarray closehisto plot 0'

   !....................................................................
   !....................................................................

      do m=1,myy
        yy(m)=0
      enddo
      do n=1,nptl
        if(n.ge.minfra.and.n.le.maxfra
     .   .and.(istptl(n).eq.5.or.istptl(n).eq.7))then
        routp=-sin(phinll)*(xptl(n)+xshft)+cos(phinll)*(yptl(n)+yshft)
        rinp = cos(phinll)*(xptl(n)+xshft)+sin(phinll)*(yptl(n)+yshft)
        if(abs(rinp).le.delrad/2.and.abs(routp).le.delrad/2)then
          rapx=dezptl(n)
          eco=0
          amt=pptl(5,n)**2+pptl(1,n)**2+pptl(2,n)**2
          if(amt.gt.0..and.pptl(4,n)+abs(pptl(3,n)).gt.0.d0)then
            amt=sqrt(amt)
            rap=sign(1.,pptl(3,n))*alog((pptl(4,n)+abs(pptl(3,n)))/amt)
            eco=amt   !*cosh(rap-rapx)
          endif
          m=(rapx+rapmax)/delrap+1
          if(m.gt.myy)m=myy
          if(m.ge.1.and.m.le.mrr)yy(m)=yy(m)+eco
        endif
        endif
      enddo
      write(ifhi,'(a)')'!---------------------------------------------'
      write(ifhi,'(a)')'!    core segment energy per d[c]dxdy         '
      write(ifhi,'(a)')'!           vs space-time rapidity rapx       '
      write(ifhi,'(a)')'!   (same as histogram rapx eco... in optns)  '
      write(ifhi,'(a)')'!---------------------------------------------'
      write(ifhi,'(a)')       '!newpage'
      write(ifhi,'(a)')    'openhisto name rapx'
      write(ifhi,'(a)')       'htyp lin xmod lin ymod lin'
      write(ifhi,'(a,2f7.3)') 'xrange ',-rapmax,rapmax
      write(ifhi,'(a)')        'yrange 0 auto '
      write(ifhi,'(a)')  'txt "title initial energy          "'
      write(ifhi,'(a,f4.1,a)')'text 0.05 0.70 "x=0"'
      write(ifhi,'(a,f4.1,a)')'text 0.05 0.60 "y=0"'
      write(ifhi,'(a,f4.1,a)')'text 0.65 0.9 "',bimevt,'fm"'
      write(ifhi,'(a)')  'txt  "xaxis space-time rapidity [c] "'
      write(ifhi,'(a)')  'txt  "yaxis dE/d[c]dxdy "'
      write(ifhi,'(a)')       'array 2'
      do m=1,myy
        write(ifhi,'(2e11.3)')-rapmax+(m-0.5)*delrap
     .             , yy(m)/delrad**2/delrap
        yy1(m)=yy(m)
      enddo
      write(ifhi,'(a)')  '  endarray closehisto plot 0-'
   !....................................................................
      delrap=2*rapmax/float(myy)
      do m=1,myy
        yy(m)=0
      enddo
      do n=1,nptl
        if(n.ge.minfra.and.n.le.maxfra.and.istptl(n)/2.eq.0)then
        routp=-sin(phinll)*(xptl(n)+xshft)+cos(phinll)*(yptl(n)+yshft)
        rinp = cos(phinll)*(xptl(n)+xshft)+sin(phinll)*(yptl(n)+yshft)
        if(abs(rinp).le.delrad/2.and.abs(routp).le.delrad/2)then
          rapx=dezptl(n)
          eco=0
          amt=pptl(5,n)**2+pptl(1,n)**2+pptl(2,n)**2
          if(amt.gt.0..and.pptl(4,n)+abs(pptl(3,n)).gt.0.d0)then
            amt=sqrt(amt)
            rap=sign(1.,pptl(3,n))*alog((pptl(4,n)+abs(pptl(3,n)))/amt)
            eco=amt   !*cosh(rap-rapx)
          endif
          m=0
          if(abs(rapx).lt.1e10)m=(rapx+rapmax)/delrap+1
          if(m.ge.1.and.m.le.myy)yy(m)=yy(m)+eco
        endif
        endif
      enddo
      write(ifhi,'(a)')'!---------------------------------------------'
      write(ifhi,'(a)')'!    corona segment energy per d[c]dxdy       '
      write(ifhi,'(a)')'!           vs space-time rapidity rapx       '
      write(ifhi,'(a)')'!   (same as histogram rapx eco... in optns)  '
      write(ifhi,'(a)')'!---------------------------------------------'
      write(ifhi,'(a)')    'openhisto name rapx'
      write(ifhi,'(a)')       'htyp lin xmod lin ymod lin'
      write(ifhi,'(a,2f7.3)') 'xrange ',-rapmax,rapmax
      write(ifhi,'(a)')        'yrange 0 auto '
      write(ifhi,'(a)')  'txt "title initial energy          "'
      write(ifhi,'(a)')  'txt  "xaxis space-time rapidity [c] "'
      write(ifhi,'(a)')  'txt  "yaxis dE/d[c]dxdy "'
      write(ifhi,'(a)')       'array 2'
      do m=1,myy
        write(ifhi,'(2e11.3)')-rapmax+(m-0.5)*delrap
     .    , (yy1(m)+yy(m))/delrad**2/delrap
      enddo
      write(ifhi,'(a)')  '  endarray closehisto plot 0-'
      write(ifhi,'(a)')  'plot 0'
   !....................................................................
   !....................................................................

      write(ifhi,'(a)')
     .'openhisto xrange 0 1  array 2  endarray closehisto plot 0'
      write(ifhi,'(a)')
     .'openhisto xrange 0 1  array 2  endarray closehisto plot 0'
      write(ifhi,'(a)')
     .'openhisto xrange 0 1  array 2  endarray closehisto plot 0'

  999 continue

      call utprix('xcorec',ish,ishini,4)
      return
      end

c-----------------------------------------------------------------------
      subroutine xIcoPlot(nev)
c-----------------------------------------------------------------------
      common/ceve/neve
      neve=nev

      call xxIcoPlot('epsI','z','lin',  0., 1)
      call xxIcoPlot('flaB','z','lin',  0., 1)
      call xxIcoPlot('epsi','x','log',  0., 0)
      call xxIcoPlot('epsi','y','log',  0., 1)

c      call xxIcoPlot('epsi','x','log',  4., 0)
c      call xxIcoPlot('epsi','x','log',  8., 0)
c      call xxIcoPlot('epsi','x','log', -4., 0)
c      call xxIcoPlot('epsi','x','log', -8., 1)
c      call xxIcoPlot('epsi','y','log',  4., 0)
c      call xxIcoPlot('epsi','y','log',  8., 0)
c      call xxIcoPlot('epsi','y','log', -4., 0)
c      call xxIcoPlot('epsi','y','log', -8., 1)
c      call xxIcoPlot('epsi','z','lin',  0., 1)
c      call xxIcoPlot('epsi','x','lin',  0., 1)

      end

c-----------------------------------------------------------------------
      subroutine xxIcoPlot(cipl,cxi,cmod,zeta,kpl)
c-----------------------------------------------------------------------
#include "aaa.h"
#include "ico.h"
      common/ceve/neve
      common/xyz/delx,dely,delz,nxzero,nyzero,nzzero
      character cx*1,cxi*1,cipl*4,name*15,cfl*1,cmod*3
      common/cname/ii,name
      real asm(10,3,5,100),iasm(10,3,5)
      logical g
      data nax/0/ nicnt/0/
      save nax,asm,iasm,nicnt
      nax=nax+1
      nicnt=nicnt+1
      g=.false.
      if(neve.gt.nevent-nfreeze)g=.true.

      if(nicnt.eq.1)then
      do i=1,10
      do j=1,3
      do k=1,5
      iasm(i,j,k)=0
      do m=1,50
      asm(i,j,k,m)=0
      enddo
      enddo
      enddo
      enddo
      endif
      ii=0
      jj=0
      kk=0
      if(cipl.eq.'epsi')ii=1
      if(cipl.eq.'velx')ii=2
      if(cipl.eq.'vely')ii=3
      if(cipl.eq.'velz')ii=4
      if(cipl.eq.'flau')ii=5
      if(cipl.eq.'flad')ii=6
      if(cipl.eq.'flas')ii=7
      if(cipl.eq.'flab')ii=8
      if(cipl.eq.'epsI')ii=9
      if(cipl.eq.'flaB')ii=10
      if(cxi.eq.'x')jj=1
      if(cxi.eq.'y')jj=2
      if(cxi.eq.'z')jj=3
      if(nint(zeta).eq.-8)kk=1
      if(nint(zeta).eq.-4)kk=2
      if(nint(zeta).eq. 0)kk=3
      if(nint(zeta).eq. 4)kk=4
      if(nint(zeta).eq. 8)kk=5
      if(neve.eq.1)then
      iasm(ii,jj,kk)=iasm(ii,jj,kk)+1
      if(iasm(ii,jj,kk).gt.1)stop'variable occurs twice'
      endif

      delx=(xmaxico-xminico)/nxico
      dely=(ymaxico-yminico)/nyico
      delz=(zmaxico-zminico)/nzico
      nxzero=nxico/2+1
      nyzero=nyico/2+1
      nzzero=nzico/2+1

      cx=cxi

      nmax=0
      vmin=0.
      vmax=0.
      delv=0.
      if(cx.eq.'x')then
        nmax=100
        vmin=-10
        vmax=10
      elseif(cx.eq.'y')then
        nmax=100
        vmin=-10
        vmax=10
      elseif(cx.eq.'z')then
        nmax=70
        vmin=-7
        vmax=7
      endif
      delv=(vmax-vmin)/nmax
      if(nmax.gt.100)stop'nmax too small'

       cfl=cipl(4:4)
       if(cfl.eq.'B'.or.cfl.eq.'I')then
         fcr=delx*dely
       else
         fcr=1
       endif

  99  format(2e13.5)

      if(cipl(1:3).eq.'eps')then

       if(g)write(ifhi,'(a)')    '!##################################'
       if(g)write(ifhi,'(a,i3)') '!   epsi       '
       if(g)write(ifhi,'(a)')    '!##################################'
       if(g)write(ifhi,'(a)') '!newpage'
       if(g)write(ifhi,'(a,i1,a,i3,a)')
     . 'openhisto name epsi-'//'-',500+nint(zeta),'-'//cx
       if(g)write(ifhi,'(a,2f9.3)')'xrange',vmin,vmax
       if(g)write(ifhi,'(a)')'yrange auto auto'
       if(g)write(ifhi,'(a)') 'htyp lin xmod lin ymod '//cmod
       if(nax.eq.1)then
       if(g)write(ifhi,'(a)')'txt  "xaxis '//cx//' (fm)"'
       if(g.and.cfl.eq.'i')
     .         write(ifhi,'(a)')'txt "yaxis [e] (GeV/fm^3!)"'
       if(g.and.cfl.eq.'I')write(ifhi,'(a)')'txt "yaxis dE/dz (GeV/fm)"'
       endif
       if(g)write(ifhi,'(a)')'array 2'
       do n=1,nmax
         v=vmin+(n-0.5)*delv
         w=0
         if(cx.eq.'x')w=w+wIcoPI3(v,0.,zeta,0)
         if(cx.eq.'y')w=w+wIcoPI3(0.,v,zeta,0)
         if(cx.eq.'z')w=w+wIcoPI1(v,0)*fcr
         asm(ii,jj,kk,n)=asm(ii,jj,kk,n)+w
         if(g)w=asm(ii,jj,kk,n)/nevent
         if(g)write(ifhi,99)v,w
       enddo
       if(g)write(ifhi,'(a)') 'endarray closehisto '
       if(kpl.eq.1)then
         if(g)write(ifhi,'(a)') 'plot 0'
       else
         if(g)write(ifhi,'(a)') 'plot 0-'
       endif

      elseif(cipl.eq.'velz')then

       if(g)write(ifhi,'(a)')    '!##################################'
       if(g)write(ifhi,'(a,i3)') '!   velz       '
       if(g)write(ifhi,'(a)')    '!##################################'
       if(g)write(ifhi,'(a)') '!newpage'
       if(g)write(ifhi,'(a,i1,a,i3,a)')
     . 'openhisto name velz-'//'-',500+nint(zeta),'-'//cx
       if(g)write(ifhi,'(a,2f9.3)')'xrange',vmin,vmax
       if(g)write(ifhi,'(a)')'yrange auto auto'
       if(g)write(ifhi,'(a)') 'htyp lin xmod lin ymod lin '
       if(nax.eq.1)then
         if(g)write(ifhi,'(a)')'txt  "xaxis '//cx//' (fm)"'
         if(g)write(ifhi,'(a)')'txt "yaxis v?z!  "'
       endif
       if(g)write(ifhi,'(a)')'array 2'
       do n=1,nmax
         v=vmin+(n-0.5)*delv
         if(cx.eq.'x')w=w+wIcoPI3(v,0.,zeta,3)
         if(cx.eq.'y')w=w+wIcoPI3(0.,v,zeta,3)
         if(cx.eq.'z')w=w+wIcoPI3(0.,0., v ,3)
         asm(ii,jj,kk,n)=asm(ii,jj,kk,n)+w
         if(g)write(ifhi,99)v,asm(ii,jj,kk,n)/nevent
       enddo
       if(g)write(ifhi,'(a)') 'endarray closehisto '
       if(kpl.eq.1)then
         if(g)write(ifhi,'(a)') 'plot 0'
       else
         if(g)write(ifhi,'(a)') 'plot 0-'
       endif

      elseif(cipl.eq.'velx')then

       if(g)write(ifhi,'(a)')    '!##################################'
       if(g)write(ifhi,'(a,i3)') '!   velx       '
       if(g)write(ifhi,'(a)')    '!##################################'
       if(g)write(ifhi,'(a)') '!newpage'
       if(g)write(ifhi,'(a,i1,a,i3,a)')
     . 'openhisto name velx-'//'-',500+nint(zeta),'-'//cx
       if(g)write(ifhi,'(a,2f9.3)')'xrange',vmin,vmax
       if(g)write(ifhi,'(a)')'yrange -1 1'
       if(g)write(ifhi,'(a)') 'htyp lin xmod lin ymod lin '
       if(nax.eq.1)then
         if(g)write(ifhi,'(a)')'txt  "xaxis '//cx//' (fm)"'
         if(g)write(ifhi,'(a)')'txt "yaxis v?x! "'
       endif
       if(g)write(ifhi,'(a)')'array 2'
       do n=1,nmax
         v=vmin+(n-0.5)*delv
         if(cx.eq.'x')w=w+wIcoPI3(v,0.,zeta,1)
         if(cx.eq.'y')w=w+wIcoPI3(0.,v,zeta,1)
         if(cx.eq.'z')w=w+wIcoPI3(0.,0., v ,1)
         asm(ii,jj,kk,n)=asm(ii,jj,kk,n)+w
         if(g)write(ifhi,99)v,asm(ii,jj,kk,n)/nevent
       enddo
       if(g)write(ifhi,'(a)') 'endarray closehisto '
       if(kpl.eq.1)then
         if(g)write(ifhi,'(a)') 'plot 0'
       else
         if(g)write(ifhi,'(a)') 'plot 0-'
       endif


      elseif(cipl.eq.'vely')then

       if(g)write(ifhi,'(a)')    '!##################################'
       if(g)write(ifhi,'(a,i3)') '!   vely       '
       if(g)write(ifhi,'(a)')    '!##################################'
       if(g)write(ifhi,'(a)') '!newpage'
       if(g)write(ifhi,'(a,i1,a,i3,a)')
     . 'openhisto name velx-'//'-',500+nint(zeta),'-'//cx
       if(g)write(ifhi,'(a,2f9.3)')'xrange',vmin,vmax
       if(g)write(ifhi,'(a)')'yrange -1 1'
       if(g)write(ifhi,'(a)') 'htyp lin xmod lin ymod lin '
       if(nax.eq.1)then
         if(g)write(ifhi,'(a)')'txt  "xaxis '//cx//' (fm)"'
         if(g)write(ifhi,'(a)')'txt "yaxis v?y! "'
      endif
       if(g)write(ifhi,'(a)')'array 2'
       do n=1,nmax
         v=vmin+(n-0.5)*delv
         if(cx.eq.'x')w=w+wIcoPI3(v,0.,zeta,2)
         if(cx.eq.'y')w=w+wIcoPI3(0.,v,zeta,2)
         if(cx.eq.'z')w=w+wIcoPI3(0.,0., v ,2)
         asm(ii,jj,kk,n)=asm(ii,jj,kk,n)+w
         if(g)write(ifhi,99)v,asm(ii,jj,kk,n)/nevent
       enddo
       if(g)write(ifhi,'(a)') 'endarray closehisto '
       if(kpl.eq.1)then
         if(g)write(ifhi,'(a)') 'plot 0'
       else
         if(g)write(ifhi,'(a)') 'plot 0-'
       endif

      elseif(cipl(1:3).eq.'fla')then

       if(g)write(ifhi,'(a)')    '!##################################'
       if(g)write(ifhi,'(a,i3)') '!   flavor       '
       if(g)write(ifhi,'(a)')    '!##################################'
       fc=1
       if(cfl.eq.'u')then
       ii1=1
       ii2=1
       fc=1
       elseif(cfl.eq.'d')then
       ii1=2
       ii2=2
       fc=1
       elseif(cfl.eq.'s')then
       ii1=3
       ii2=3
       fc=1
       elseif(cfl.eq.'b'.or.cfl.eq.'B')then
       ii1=1
       ii2=3
       fc=1./3.
       endif
       if(g)write(ifhi,'(a)') '!newpage'
       if(g)write(ifhi,'(a,2i1,a,i1,a,i3,a)')
     . 'openhisto name falvor-',ii1,ii2,'-'//'-',500+nint(zeta),'-'//cx
       if(g)write(ifhi,'(a)')'openhisto'
       if(g)write(ifhi,'(a,2f9.3)')'xrange',vmin,vmax
       if(g)write(ifhi,'(a)')'yrange auto auto '
       if(g)write(ifhi,'(a)') 'htyp lin xmod lin ymod lin '
       if(nax.eq.1)then
         if(g)write(ifhi,'(a)')'txt  "xaxis '//cx//' (fm)"'
         if(g)write(ifhi,'(a)')
     .    'txt "yaxis f?'//cfl//'!'//' "'
       endif
       if(g)write(ifhi,'(a)')'array 2'
       do n=1,nmax
         v=vmin+(n-0.5)*delv
         w=0
         do ii=ii1,ii2
           if(cx.eq.'x')w=w+wIcoPI3(v,0.,zeta,ii+3)*fc
           if(cx.eq.'y')w=w+wIcoPI3(0.,v,zeta,ii+3)*fc
           if(cx.eq.'z')w=w+wIcoPI1(v,ii+3)*fc*fcr
         enddo
         asm(ii,jj,kk,n)=asm(ii,jj,kk,n)+w
         if(g)write(ifhi,99)v,asm(ii,jj,kk,n)/nevent
       enddo
       if(g)write(ifhi,'(a)') 'endarray closehisto '
       if(kpl.eq.1)then
         if(g)write(ifhi,'(a)') 'plot 0'
       else
         if(g)write(ifhi,'(a)') 'plot 0-'
       endif

      endif

      if(kpl.eq.1)nax=0
      end

c-----------------------------------------------------------------------
      function wIcoPI3(x,y,z,m) !polynomial interpolation
c-----------------------------------------------------------------------
#include "aaa.h"
#include "ico.h"
      real wx(2),wy(2),wz(2)
      tau0 = tauzer
      delxico=(xmaxico-xminico)/nxico
      delyico=(ymaxico-yminico)/nyico
      delzico=(zmaxico-zminico)/nzico
      w=0
      jx=1+(x-xminico-delxico/2)/delxico
      jy=1+(y-yminico-delyico/2)/delyico
      jz=1+(z-zminico-delzico/2)/delzico
      if(    jx.gt.0.and.jx.lt.nxico
     .  .and.jy.gt.0.and.jy.lt.nyico
     .  .and.jz.gt.0.and.jz.lt.nzico)then
        xj=xminico+(jx-0.5)*delxico
        yj=yminico+(jy-0.5)*delyico
        zj=zminico+(jz-0.5)*delzico
        frac=(x-xj)/delxico
        wx(1)=1-frac
        wx(2)=frac
        frac=(y-yj)/delyico
        wy(1)=1-frac
        wy(2)=frac
        frac=(z-zj)/delzico
        wz(1)=1-frac
        wz(2)=frac
        do i=1,2
        do j=1,2
        do k=1,2
          select case (m)   
          case(0)
        w=w+wx(i)*wy(j)*wz(k)*IcoE(  jx+i-1,jy+j-1,jz+k-1)
          case(1)
        w=w+wx(i)*wy(j)*wz(k)*IcoV(1,jx+i-1,jy+j-1,jz+k-1)
          case(2)
        w=w+wx(i)*wy(j)*wz(k)*IcoV(2,jx+i-1,jy+j-1,jz+k-1)
          case(3)
        w=w+wx(i)*wy(j)*wz(k)*IcoV(3,jx+i-1,jy+j-1,jz+k-1)*tau0
          case(4)
        w=w+wx(i)*wy(j)*wz(k)*IcoF(1,jx+i-1,jy+j-1,jz+k-1)
          case(5)
        w=w+wx(i)*wy(j)*wz(k)*IcoF(2,jx+i-1,jy+j-1,jz+k-1)
          case(6)
        w=w+wx(i)*wy(j)*wz(k)*IcoF(3,jx+i-1,jy+j-1,jz+k-1)
          end select 
        enddo
        enddo
        enddo
      else
        stop'xIcoPI Out of range'
      endif
      wIcoPI3=w
      end


c-----------------------------------------------------------------------
      function wIcoPI1(z,m) !polynomial interpolation
c-----------------------------------------------------------------------
#include "aaa.h"
#include "ico.h"
      real wz(2)
      tau0 = tauzer
      delzico=(zmaxico-zminico)/nzico
      w=0
      jz=1+(z-zminico-delzico/2)/delzico
      if(jz.gt.0.and.jz.lt.nzico)then
        zj=zminico+(jz-0.5)*delzico
        frac=(z-zj)/delzico
        wz(1)=1-frac
        wz(2)=frac
        do i=1,nxico
        do j=1,nyico
        do k=1,2
          select case (m)   
          case(0)
        w=w+wz(k)*IcoE(  i,j,jz+k-1)
          case(1)
        w=w+wz(k)*IcoV(1,i,j,jz+k-1)
          case(2)
        w=w+wz(k)*IcoV(2,i,j,jz+k-1)
          case(3)
        w=w+wz(k)*IcoV(3,i,j,jz+k-1)*tau0
          case(4)
        w=w+wz(k)*IcoF(1,i,j,jz+k-1)
          case(5)
        w=w+wz(k)*IcoF(2,i,j,jz+k-1)
          case(6)
        w=w+wz(k)*IcoF(3,i,j,jz+k-1)
          end select 
        enddo
        enddo
        enddo
      else
        stop'xIcoPI Out of range'
      endif
      wIcoPI1=w
      end



