C
C  This file is part of EPOS4
C  Copyright (C) 2022 research institutions and authors (See CREDITS file)
C  This file is distributed under the terms of the GNU General Public License version 3 or later
C  (See COPYING file for the text of the licence)
C

c-----------------------------------------------------------------------
      subroutine WriteHydroTable
c-----------------------------------------------------------------------
c      use ArrayModule, only: getd
c      use hocoModule, only: epsc, velc, barc, sigc
c      use hocoModule, only: velc, barc, sigc
      double precision getHydynEpsc, getHydynVelc 
      double precision getHydynBarc, getHydynSigc

#include "aaa.h"
#include "ho.h"
      
      character*80 fn
      real val
      data ncnthyt/0/
      save ncnthyt
      ncnthyt=ncnthyt+1
      !!!if(ncnthyt.ne.8)return ! for debugging
      if(igethyt.lt.2)return
      tauminhy=getHydynTauhy(1)
      taumaxhy=getHydynTauhy(ntauhy)
      ii=index(fnhi(1:nfnhi),".histo")-1
      fn=fnhi(1:ii)//".out"
      iixx=ii+4
      if(ifmt.ne.6)open(ifmt,file=fnmt(1:nfnmt),access='append')
      write(ifmt,'(2a)')'write hydro table into ',fn(1:iixx)
      call clop(3)
      open(11,file=fn(1:iixx),status='unknown')
      write(11,*)iversn
      write(11,*)laproj,maproj,latarg,matarg
      write(11,*)engy
      write(11,*)bminim,bmaxim,ikolmn,ikolmx
      write(11,*)tauzer,ifathlle,dtauhy
      write(11,*)0
      write(11,*)nzhy,ntauhy,nxhy,nyhy
      write(11,*)zminhy,zmaxhy,tauminhy,taumaxhy
      write(11,*)xminhy,xmaxhy,yminhy,ymaxhy

c      write(11,*)
c     *  (((((getd(velc,[i,neta,ntau,nx,ny]),i=1,3),neta=1,nzhy)
c     *          ,ntau=1,ntauhy),nx=1,nxhy),ny=1,nyhy)
c     * ,((((getd(epsc,[neta,ntau,nx,ny]),neta=1,nzhy)
c     *          ,ntau=1,ntauhy),nx=1,nxhy),ny=1,nyhy)
c     * ,((((getd(sigc,[neta,ntau,nx,ny]),neta=1,nzhy)
c     *          ,ntau=1,ntauhy),nx=1,nxhy),ny=1,nyhy)
c     * ,(((((getd(barc,[i,neta,ntau,nx,ny]),i=1,3),neta=1,nzhy)
c     *          ,ntau=1,ntauhy),nx=1,nxhy),ny=1,nyhy)


      do ny=1,nyhy
         do nx=1,nxhy
            do ntau=1,ntauhy
               do neta=1,nzhy
                  do i=1,3
                     write(11,*)getHydynVelc(i,neta,ntau,nx,ny)
                  end do
               end do 
            end do
         end do
      end do

      do ny=1,nyhy
         do nx=1,nxhy
            do ntau=1,ntauhy
               do neta=1,nzhy
                     write(11,*)getHydynEpsc(neta,ntau,nx,ny)
               end do 
            end do
         end do
      end do

      do ny=1,nyhy
         do nx=1,nxhy
            do ntau=1,ntauhy
               do neta=1,nzhy
                     write(11,*)getHydynSigc(neta,ntau,nx,ny)                  
               end do 
            end do
         end do
      end do

      do ny=1,nyhy
         do nx=1,nxhy
            do ntau=1,ntauhy
               do neta=1,nzhy
                  do i=1,3
                     write(11,*)getHydynBarc(i,neta,ntau,nx,ny)
                  end do
               end do 
            end do
         end do
      end do

      do ny=1,nyhy
       do nx=1,nxhy
        do ntau=1,ntauhy
         do nz=1,nzhy
          do ivi=1,10
            call emucget(ivi,nz,ntau,nx,ny,val)
            write(11,*)val
          enddo !ivi
         enddo !nz
        enddo !ntau
       enddo !nx
      enddo !ny

      close(11)
      end

c-----------------------------------------------------------------------
      subroutine ReadHydroTable
c-----------------------------------------------------------------------
c      use ArrayModule, only: setDouble, getd
c      use hocoModule, only: epsc, velc, barc, sigc
c      use hocoModule, only: velc, barc, sigc
      double precision getHydynEpsc, getHydynVelc
      double precision getHydynBarc, getHydynSigc
      
#include "aaa.h"
#include "ho.h"
#include "ems.h"
      logical table
      common/jcen/jcentr,jmxcentr/kcen/kcentr,kmxcentr
      common/czhits/zhits(0:ncenthxx),zfhits(0:ncenthxx) 
      common/czkhits/zkhits(0:ncenthxx)
      common/cntauhyori/ntauhyori
      real pixx(10)
      integer itsy
      real ww
      double precision readData
      real val
      common/citsy/itsy(4,4)
      real wi(3),wj(3)
      data nnn/0/
      data ncntrdhy/0/
      save ncntrdhy
      
      if(ireadhyt.eq.0)return

      if(ifmt.ne.6)open(ifmt,file=fnmt(1:nfnmt),access='append')

      nnn=nnn+1

      inquire(file=fnho(nnn)(1:nfnho(nnn)),exist=table)
      if(.not.table)then
        write(ifmt,'(//10x,3a//)')'STOP in ReadHydroTable: file "'
     *             ,fnho(nnn)(1:nfnho(nnn)),'" not found !!!'
        stop
      endif
      write(ifmt,'(3a,$)')'read from '
     * ,fnho(nnn)(1:nfnho(nnn)),' ... '
      open(11,file=fnho(nnn)(1:nfnho(nnn)),status='old')
      read(11,*)iversn
      if(iversn.lt.201)stop'ERROR redo table with version >= 201'
      read(11,*)laprojx,maprojx,latargx,matargx
      read(11,*)engyx
      read(11,*)bminimx,bmaximx,ikolmnx,ikolmxx
      bevt=bminimx
      if(bmaximx**2-bminimx**2.ne.0.)then
        bevt=2./3.*(bmaximx**3-bminimx**3)/(bmaximx**2-bminimx**2)
      endif
      bimevt=bevt
      nprt1=(ikolmnx+ikolmxx)/2
      nprt(1)=nprt1
      phievt=0
      ranphi=0
      read(11,*)tauzer,ifathlle,dtauhy
      read(11,*)iabs_ninicon
      read(11,*)nzhy,ntauhy,nxhy,nyhy
      read(11,*)zminhy,zmaxhy,tauminhy,taumaxhy
      read(11,*)xminhy,xmaxhy,yminhy,ymaxhy
      if(laprojx.ne.laproj.or.maprojx.ne.maproj
     . .or.latargx.ne.latarg.or.matargx.ne.matarg)
     .stop'***** ERROR in ReadHydroTable: different nuclei *****'
      if(engyx.ne.engy)
     .stop'***** ERROR in ReadHydroTable: different energy *****'

c      read(11,*)
c     *  (((((velc(i,neta,ntau,nx,ny),i=1,3),neta=1,nzhy)
c     *          ,ntau=1,ntauhy),nx=1,nxhy),ny=1,nyhy)
c     * ,((((epsc(neta,ntau,nx,ny),neta=1,nzhy)
c     *          ,ntau=1,ntauhy),nx=1,nxhy),ny=1,nyhy)
c     * ,((((sigc(neta,ntau,nx,ny),neta=1,nzhy)
c     *          ,ntau=1,ntauhy),nx=1,nxhy),ny=1,nyhy)
c     * ,(((((barc(i,neta,ntau,nx,ny),i=1,3),neta=1,nzhy)
c     *          ,ntau=1,ntauhy),nx=1,nxhy),ny=1,nyhy)

      iemuc=1
      call memo(1,'create emuc;')
      call createemuc(10,nzhy,ntauhxx,nxhy,nyhy)
      call memo(2,';')

      do ny=1,nyhy
         do nx=1,nxhy
            do ntau=1,ntauhy
               do neta=1,nzhy
                  do i=1,3
                     read(11,*) readData
                     call setHydynVelc(i,neta,ntau,nx,ny,readData)
                  end do
               end do 
            end do
         end do
      end do

      do ny=1,nyhy
         do nx=1,nxhy
            do ntau=1,ntauhy
               do neta=1,nzhy
                  read(11,*) readData
                  call setHydynEpsc(neta,ntau,nx,ny,readData)
               end do 
            end do
         end do
      end do

      do ny=1,nyhy
         do nx=1,nxhy
            do ntau=1,ntauhy
               do neta=1,nzhy
                  read(11,*) readData
                  call setHydynSigc(neta,ntau,nx,ny,readData)                  
               end do 
            end do
         end do
      end do

      do ny=1,nyhy
         do nx=1,nxhy
            do ntau=1,ntauhy
               do neta=1,nzhy
                  do i=1,3
                     read(11,*) readData
                     call setHydynBarc(i,neta,ntau,nx,ny,readData)
                  end do
               end do 
            end do
         end do
      end do

      do ny=1,nyhy
       do nx=1,nxhy
        do ntau=1,ntauhy
         do nz=1,nzhy
          do ivi=1,10
            read(11,*,end=88)val
            call emucset(ivi,nz,ntau,nx,ny, val )
          enddo !ivi
         enddo !nz
        enddo !ntau
       enddo !nx
      enddo !ny

      close(11)
      write(ifmt,'(a)')'  done'
      write(ifmt,'(a,f7.2)')'impact parameter =',bevt
      write(ifmt,'(a,i5,f7.2)')'nzhy,zmaxhy =',nzhy,zmaxhy
      call clop(3)

      goto77

  88  continue 
      write(ifmt,'(5i5)')ivi,nz,ntau,nx,ny
      stop'ERROR 230805 End of file'  

      !-----------------------!
       entry ManiHydroTable
      !-----------------------!

      nnn=nnn+1
      if(ifmt.ne.6)open(ifmt,file=fnmt(1:nfnmt),access='append')
      ntauhyori=ntauhy
        !if(ntauhy.le.2)ntauhy=0
      tauminhy=tauzer
      if(ntauhy.gt.0)then
      taumaxhy=getHydynTauhy(ntauhy)
      else
      taumaxhy=tauminhy
      endif
      bevt=bimevt
      nprt1=nprt(1)
        !if(ntauhy.eq.0.and.tauup.gt.0.05)then
        !do n=1,nptl
        !if(istptl(n).eq.7)istptl(n)=0
        !enddo
        !endif
      
  77  continue

      do ntau=1,ntauhxx
       mmsteps=4*ifathlle        !must identical to def in h.f
       del_tau=dtauhy/ifathlle   !must identical to def in h.f
       delta_tau_out = mmsteps * del_tau
       tau=tauminhy+(ntau-1)*delta_tau_out !for output hydro table
       call setHydynTauhy(ntau,tau)
      enddo

      if(nnn.eq.1)then
        do n=0,jmxcentr
        zhits(n)=0
        zkhits(n)=0
        zfhits(n)=0
        enddo
      endif

      do neta=1,nzhy
       eta=zminhy+(neta-1)*(zmaxhy-zminhy)/(nzhy-1)
       etahy(neta)=eta
      enddo

      rmaxhy=max(xmaxhy,ymaxhy)
      do nrad=1,nradhy
       rad=(nrad-1.)/(nradhy-1.)*rmaxhy
       radhy(nrad)=rad
      enddo

      call getJcentr
      call getKcentr
      zhits(jcentr)=zhits(jcentr)+1
      zkhits(kcentr)=zkhits(kcentr)+1
      if(ntauhy.ne.0)zfhits(jcentr)=zfhits(jcentr)+1
      write(ifmt,'(a,2i5,f9.1)')
     .'hyt.f ManiHydroTable: ntauhy jcentr zfhits = '
     .,ntauhy,jcentr,zfhits(jcentr)
  
      if(ntauhy.eq.0)goto999

      call memo(1,'create hydro tables;')
      call createvelio(3 ,nzhy,ntauhy,nphihy,nradhy)
      call createbario(3 ,nzhy,ntauhy,nphihy,nradhy)
      call createepsio(2 ,nzhy,ntauhy,nphihy,nradhy)
      call createemuio(10,nzhy,ntauhy,nphihy,nradhy)
      call memo(2,';')

      !neta=1
      !call xTests(neta)

      write(ifmt,'(a,$)')'calculating averages ...'
      do neta=1,nzhy
        do ntau=1,ntauhy
          esum=0
          xx=0
          yy=0
          vxx=0
          vyy=0
          vt=0
          do nx=1,nxhy
          do ny=1,nyhy
             e=max(0.d0,getHydynEpsc(neta,ntau,nx,ny))
             esum=esum+e
             x=xminhy+(nx-1)*(xmaxhy-xminhy)/(nxhy-1)
             y=yminhy+(ny-1)*(ymaxhy-yminhy)/(nyhy-1)
             xx=xx+x**2*e
             yy=yy+y**2*e
             vx=getHydynVelc(1,neta,ntau,nx,ny)
             vy=getHydynVelc(2,neta,ntau,nx,ny)
             vxx=vxx+vx**2*e
             vyy=vyy+vy**2*e
             if(vx**2+vy**2.ne.0.)
     .            vt=vt+sqrt(vx**2+vy**2)*e
          enddo
          enddo
          if(xx+yy.ne.0.)
     .     eccxav(neta,ntau)=(yy-xx)/(xx+yy)
          if(esum.ne.0.)
     .     vt=vt/esum
          eccpav(neta,ntau)=0
          if(vt*vt.gt.0.001)
     .      eccpav(neta,ntau)=(vxx-vyy)/(vxx+vyy)
          vtraav(neta,ntau)=vt
        enddo
      enddo
      write(ifmt,'(a)')'  done'

      write(ifmt,'(a,$)')'transforming to polar coordinates ...'
      delx=(xmaxhy-xminhy)/(nxhy-1)
      dely=(ymaxhy-yminhy)/(nyhy-1)
      do nphi=1,nphihy
       phi=(nphi-1.)/(nphihy-1.)*2*pi
       phihy(nphi)=phi
      enddo
      do neta=1,nzhy
       eta=etahy(neta)
       gam=cosh(eta)
       do ntau=1,ntauhy
        call getxyrange(neta,ntau)
        do nphi=1,nphihy
         phi=phihy(nphi)
         do nrad=1,nradhy
           !~~~epsio,velio,bario initialization not needed any more
           rad=radhy(nrad)
           x=rad*cos(phi)
           y=rad*sin(phi)
           vx=0
           vy=0
           vz=0
           ep=0
           ba=0
           b2=0
           b3=0
           sg=0
           do ivi=1,10
             pixx(ivi)=0
           enddo
           xx=(x-xminhy)/delx +1
           yy=(y-yminhy)/dely +1
           nx=int(xx)
           ny=int(yy)
           ipoint=2
           if(ipoint.eq.3)then !~~~~3point interpolation~~~~
           if(nx.lt.1)nx=1
           if(ny.lt.1)ny=1
           if(nx.gt.(nxhy-2))nx=nxhy-2
           if(ny.gt.(nyhy-2))ny=nyhy-2
           f=xx-nx
           fii=f
           wi(3)=f*(f-1.)*.5
           wi(1)=1.-f+wi(3)
           wi(2)=f-2.*wi(3)
           f=yy-ny
           fjj=f
           wj(3)=f*(f-1.)*.5
           wj(1)=1.-f+wj(3)
           wj(2)=f-2.*wj(3)
           if(  fii.ge.0.0.and.fii.le.2.0
     .     .and.fjj.ge.0.0.and.fjj.le.2.0)then
            do i=1,3
             do j=1,3
              do ivi=1,10
               call emucget(ivi,neta,ntau,nx+i-1,ny+j-1, ww )
               pixx(ivi)=pixx(ivi)+ ww * wi(i)*wj(j)
              enddo
              vx=vx+getHydynVelc(1,neta,ntau,nx+i-1,ny+j-1)*wi(i)*wj(j)
              vy=vy+getHydynVelc(2,neta,ntau,nx+i-1,ny+j-1)*wi(i)*wj(j)
              vz=vz+getHydynVelc(3,neta,ntau,nx+i-1,ny+j-1)*wi(i)*wj(j)
              ep=ep+getHydynEpsc(  neta,ntau,nx+i-1,ny+j-1)*wi(i)*wj(j)
              ba=ba+getHydynBarc(1,neta,ntau,nx+i-1,ny+j-1)*wi(i)*wj(j)
              b2=b2+getHydynBarc(2,neta,ntau,nx+i-1,ny+j-1)*wi(i)*wj(j)
              b3=b3+getHydynBarc(3,neta,ntau,nx+i-1,ny+j-1)*wi(i)*wj(j)
              sg=sg+getHydynSigc(  neta,ntau,nx+i-1,ny+j-1)*wi(i)*wj(j)
             enddo
            enddo
           endif
           elseif(ipoint.eq.2)then !~~~~2point interpolation~~~~
           if(nx.ge.1.and.ny.ge.1
     .    .and.nx.le.nxhy-1.and.ny.le.nyhy-1)then
           frac=xx-nx
           wi(1)=1-frac
           wi(2)=frac
           frac=yy-ny
           wj(1)=1-frac
           wj(2)=frac
           do i=1,2
            do j=1,2
              do ivi=1,10
               call emucget(ivi,neta,ntau,nx+i-1,ny+j-1, ww )
               pixx(ivi)=pixx(ivi)+ ww * wi(i)*wj(j)
              enddo
              vx=vx+getHydynVelc(1,neta,ntau,nx+i-1,ny+j-1)* wi(i)*wj(j)
              vy=vy+getHydynVelc(2,neta,ntau,nx+i-1,ny+j-1)* wi(i)*wj(j)
              vz=vz+getHydynVelc(3,neta,ntau,nx+i-1,ny+j-1)* wi(i)*wj(j)
              ep=ep+getHydynEpsc(  neta,ntau,nx+i-1,ny+j-1)* wi(i)*wj(j)
              ba=ba+getHydynBarc(1,neta,ntau,nx+i-1,ny+j-1)* wi(i)*wj(j)
              b2=b2+getHydynBarc(2,neta,ntau,nx+i-1,ny+j-1)* wi(i)*wj(j)
              b3=b3+getHydynBarc(3,neta,ntau,nx+i-1,ny+j-1)* wi(i)*wj(j)
              sg=sg+getHydynSigc(  neta,ntau,nx+i-1,ny+j-1)* wi(i)*wj(j)
            enddo
           enddo
           endif
           endif !~~~~~~~~~~~~~~~~~~
           if(vz.gt.1.)stop'***** ERROR in ManiHydroTable: vz > 1.'
           if(ep.ge.20*oEeos.and.50.gt.ncntrdhy.and.ireadhyt.ne.1)then
            ncntrdhy=ncntrdhy+1
            write(ifmt,*)'ManiHydroTable: energy density very big: ',ep
            print*,oEeos,ireadhyt
            !write(ifmt,*)'         e(',neta,ntau,nphi,nrad,')=',ep
            !write(ifmt,*)nx,ny
            !.       ,((epsc(neta,ntau,nx+i-1,ny+j-1),j=1,2),i=1,2)
           endif
c           write(ifch,'(a,4i5,3f8.3)')'hyt  c t f r  b123'
c     .     ,neta,ntau,nphi,nrad,ba,b1,b2
           call velioset(1,neta,ntau,nphi,nrad,  vx)
           call velioset(2,neta,ntau,nphi,nrad,  vy)
           call velioset(3,neta,ntau,nphi,nrad,  vz)
           call epsioset(1,neta,ntau,nphi,nrad,  ep)
           call barioset(1,neta,ntau,nphi,nrad,  ba)
           call barioset(2,neta,ntau,nphi,nrad,  b2)
           call barioset(3,neta,ntau,nphi,nrad,  b3)
           call epsioset(2,neta,ntau,nphi,nrad,  sg)
           do ivi=1,10
           call emuioset(ivi,neta,ntau,nphi,nrad, pixx(ivi) )
           enddo
         enddo !nrad
        enddo !nphi
       enddo !ntau
      enddo !neta

      !~~~ for plots~~~~
      do neta=1,nzhy
       do ntau=1,ntauhy
        call epsioget(1,neta,ntau,1,1, epsij(neta,ntau) )
        do kij=1,3
         call barioget(kij,neta,ntau,1,1, barij(kij,neta,ntau) )
        enddo
        call epsioget(2,neta,ntau,1,1, sigij(neta,ntau) )
       enddo
      enddo
      write(ifmt,'(a)')'  done'

      call eflowcheck(3)

 999  continue 
      call clop(3)
      end
c-----------------------------------------------------------------------
      subroutine eflowcheck(ii)
c-----------------------------------------------------------------------
#include "aaa.h"
#include "ho.h"
      real tmunu(4,4),t(10),s(4),pf(4)
      common/citsy/itsy(4,4)

      if(ii.eq.0)return

      if(ii.eq.1)then ! fixed time surface, x-y (agrees perfectly with tests in h.f)

      d1= (xmaxhy-xminhy) /(nxhy-1)
      d2= (ymaxhy-yminhy) /(nyhy-1)
      d3= (zmaxhy-zminhy) /(nzhy-1)
      do ntau=1,ntauhy
       tau=getHydynTauhy(ntau)
       eetau=0
       do nz=1,nzhy
        do nx=1,nxhy
         do ny=1,nyhy
           z=zminhy  +(nz-1)  *d3
           eta=z
           do i=1,4
             do j=1,4
               !\tilde{T}^ij
               call emucget( itsy(i,j) ,nz,ntau,nx,ny, tmunu(i,j) )
             enddo
           enddo
           ! \tilde{T}^ij*\tilde{dSigma}_j
           ! Then trafo to Minkowski, lab, 
           eetau=eetau
     .     + ( tmunu(3,4)*sinh(eta)+tmunu(4,4)*cosh(eta) )*tau*d1*d2*d3
         enddo
        enddo
       enddo
       write(ifmt,'(a,f12.3)')'Eflow through t=const-surface:',eetau
      enddo
      return

      elseif(ii.eq.2)then ! fixed time surface, r-phi  (r range too small to check energy conservation for large t)

      d1=radhy(2)-radhy(1)
      d2=phihy(2)-phihy(1)
      d3= (zmaxhy-zminhy) /(nzhy-1)
      !etausu=0
      do ntau=1,ntauhy
       tau=getHydynTauhy(ntau)
       eetau=0
       do nz=1,nzhy
        do nphi=2,nphihy
         phi=phihy(nphi)
         do nrad=2,nradhy
         rad=radhy(nrad)
           eta=zminhy  +(nz-1)  *d3
           do i=1,4
             do j=1,4
               !\tilde{T}^ij
               call emuioget( itsy(i,j) ,nz,ntau,nphi,nrad, tmunu(i,j) )
             enddo
           enddo
           ! \tilde{T}^ij*\tilde{dSigma}_j
           ! Then trafo to Minkowski, lab
           eetau=eetau
     .    +(tmunu(3,4)*sinh(eta)+tmunu(4,4)*cosh(eta))*tau*rad*d1*d2*d3
         enddo
        enddo
       enddo
       write(ifmt,'(a,f12.3,a,f7.3)')'Eflow through t=const-surface-2:'
     . ,eetau, '   radmax = ', radhy(nradhy)
      enddo
      return

      elseif(ii.eq.3)then! fixed radius for all times + fixed final time (agrees perfectly with tests in h.f)

      d1=getHydynTauhy(2)-getHydynTauhy(1)
      d2=phihy(2)-phihy(1)
      d3= (zmaxhy-zminhy) /(nzhy-1)
      etausu=0
      ! fixed radius for all times 
      do ntau=2,ntauhy
       tau=(getHydynTauhy(ntau-1)+getHydynTauhy(ntau))/2
       eetau=0
       do nz=1,nzhy
        do nphi=2,nphihy
         phi=phihy(nphi)
         nrad=nradhy
         rad=radhy(nrad)
           eta=zminhy  +(nz-1)  *d3
           do i=1,10
             call emuioget( i ,nz,ntau-1,nphi,nrad, w1 )
             call emuioget( i ,nz,ntau  ,nphi,nrad, w2 )
             t(i)=(w1+w2)/2
           enddo
           ! \tilde{T}^ij*\tilde{dSigma}_j
           ! Then trafo to Minkowski, lab
           s(1)=tau*rad*cos(phi)
           s(2)=tau*rad*sin(phi) 
           s(3)=0
           s(4)=0
           do i=1,4
            pf(i)=0
            do n=1,4
              pf(i) = pf(i) + t ( itsy(i,n) ) * s(n)
            enddo
           enddo
           eetau=eetau+(pf(3)*sinh(eta)+pf(4)*cosh(eta))*d1*d2*d3
        enddo
       enddo
       etausu=etausu+eetau
       if(ish.ge.4)then
         if(ntauhy-ntau.lt.5)
     .   write(ifmt,'(a,3f12.3)')'Eflow through r=const-surface:',
     .   eetau, etausu,tau
       endif
      enddo
      ! fixed final time
      d1=radhy(2)-radhy(1)
      d2=phihy(2)-phihy(1)
      d3= (zmaxhy-zminhy) /(nzhy-1)
      ntau=ntauhy
       tau=getHydynTauhy(ntau)
       eetau=0
       do nz=1,nzhy
        do nphi=2,nphihy
         phi=phihy(nphi)
         do nrad=2,nradhy
         rad=radhy(nrad)
           eta=zminhy  +(nz-1)  *d3
           do i=1,4
             do j=1,4
               !\tilde{T}^ij
               call emuioget( itsy(i,j) ,nz,ntau,nphi,nrad, tmunu(i,j) )
             enddo
           enddo
           ! \tilde{T}^ij*\tilde{dSigma}_j
           ! Then trafo to Minkowski, lab
           eetau=eetau
     .    +(tmunu(3,4)*sinh(eta)+tmunu(4,4)*cosh(eta))*tau*rad*d1*d2*d3
         enddo
        enddo
       enddo
       etausu=etausu+eetau
       if(ish.ge.4)then
         write(ifmt,'(a,2f12.3)')'Eflow through t=const-surface:',
     .   eetau, etausu
       endif
      return
 
      elseif(ii.eq.4)then ! rad = rmax to zero  (cone) (agrees within 5% with tests in h.f, for the default grid)
                                         !   agrees perfectly with test kfrout = 3  in f.f  , for the same grid     

      d1=getHydynTauhy(2)-getHydynTauhy(1)
      d2=phihy(2)-phihy(1)
      d3= (zmaxhy-zminhy) /(nzhy-1)
      etausu=0
      ntaumax=12 ! arbitrary choice
      do ntau=1,ntaumax  !in principle the first step should only count 0.5, but contribution is anyway very small
       tau=getHydynTauhy(ntau)
       eetau=0
       do nz=1,nzhy
        eta=zminhy  +(nz-1)  *d3
        do nphi=2,nphihy
          phi=phihy(nphi)
          rad= radhy(nradhy) *
     .         ( 1 - (tau-getHydynTauhy(1))/
     .         (getHydynTauhy(ntaumax)-getHydynTauhy(1)) )  
          rad=max(0.,rad)
          nrad=2
          do while (radhy(nrad).lt.rad)
            nrad=nrad+1
          enddo
          r1=radhy(nrad-1)
          r2=radhy(nrad)
          f=(rad-r1)/(r2-r1)
          do i=1,10
            call emuioget( i ,nz,ntau,nphi,nrad-1, w1 )
            call emuioget( i ,nz,ntau,nphi,nrad,   w2 )
            t(i)=w1+f*(w2-w1)
          enddo
          ! \tilde{T}^ij*\tilde{dSigma}_j
          ! Then trafo to Minkowski, lab
          s(1)=tau*rad*cos(phi)
          s(2)=tau*rad*sin(phi) 
          s(3)=0
          s(4)=tau*rad*radhy(nradhy)/
     .         (getHydynTauhy(ntaumax)-getHydynTauhy(1))
          do i=1,4
           pf(i)=0
           do n=1,4
             pf(i) = pf(i) + t ( itsy(i,n) ) * s(n)
           enddo
          enddo
          eetau=eetau+(pf(3)*sinh(eta)+pf(4)*cosh(eta))*d1*d2*d3
        enddo
       enddo
       etausu=etausu+eetau
       write(ifmt,'(a,2f12.3)')'Eflow through cone-surface:'
     . ,eetau, etausu
      enddo  

      endif

      end

c-----------------------------------------------------------------------
      subroutine jetfluid
c-----------------------------------------------------------------------
c      use ArrayModule, only: getd
c      use hocoModule, only: epsc, velc, barc
c      use hocoModule, only: velc, barc
      double precision getHydynEpsc, getHydynVelc, getHydynBarc
      
#include "aaa.h"
#include "ho.h"
      common/cranphi/ranphi 
      common/cijetfluid/ijetfluid
      character txt*8
      logical lplot
      real temcxx(netahxx,ntauhxx,nxhxx,nyhxx)
      !real arxy(nxhxx,nyhxx)
      real arxy2(nxhxx,nyhxx),arxy3(nxhxx,nyhxx)
      real velo(3),velox(3),qq(3),u(3),p(5),px(5)
      integer ndim(4)
      double precision xo(4),v(3),xmin(4),del(4),x,y,z,t,tau,xoo(4)
      double precision za,ta,vz,a,b,c,dt,tauxx
      double precision gmbj,etaxx,v1til,v2til,v3til,gmtil,v1xx,v2xx,v3xx
      double precision rapxx
      real q(4)
      real velcxx(3,netahxx,ntauhxx,nxhxx,nyhxx)
      dimension velx(3,nxhxx,nyhxx),temx(nxhxx,nyhxx)
      dimension velxx(3,nxhxx,nyhxx),temxx(nxhxx,nyhxx)
c      common/cee1ico/ee1ico,eistico,ee1hll
      common/ceefrac/eefrac
      real wk(2),wi(2),wj(2)
      parameter (nquark=3)  
      integer ic(2),jc(nflav,2)
      common /camaq/amaq(0:3) 
      real quamass(3)
      !data quamass/0.,0.,0.120/   !from Yuri, he took from Hama
      ! little effect since m < T and m < pt, and not well justified
      data quamass/0.330,0.330,0.500/   !much more effective and better justified
      call utpri('jetflu',ish,ishini,4)
      
      if(ihlle+ioeos.le.0)return
      if(ntauhy.eq.0)return
      if(ijetfluid.eq.1)then
        return  ! normal case
      endif
      
      !write(ifmt,*)'    +-----------------------+'
      !write(ifmt,*)'    |   jetfluid            |'
      !write(ifmt,*)'    +-----------------------+'
      
      tfo=tfrout
      
      eecobf=0
      eecobfx=0
      eegain=0
      
      nxhxxx=nxhxx
      nyhxxx=nyhxx
      phinull=phievt+ranphi 
      if(ish.ge.7)
     .    write(ifmt,'(a,i3)')'phinull=',nint(phinull/3.14159*180)
      do iz=1,nzhy
        do ntau=1,ntauhy
          lplot=.false.
          if(ish.ge.7.and.iz.eq.nzhy/2+1
     .    .and.(ntau.eq.8.or.ntau.eq.10))lplot=.true.
          do ix=1,nxhy
          do iy=1,nyhy
            eps=getHydynEpsc(  iz,ntau,ix,iy)
            b1 =getHydynBarc(1,iz,ntau,ix,iy)
            b2 =getHydynBarc(2,iz,ntau,ix,iy)
            b3 =getHydynBarc(3,iz,ntau,ix,iy)
            temx(ix,iy)=PiEos(5,eps,b1,b2,b3)
            velx(1,ix,iy)=getHydynVelc(1,iz,ntau,ix,iy)
            velx(2,ix,iy)=getHydynVelc(2,iz,ntau,ix,iy)
            velx(3,ix,iy)=getHydynVelc(3,iz,ntau,ix,iy)
            !tm=temx(ix,iy)
            !arxy(ix,iy)=0
            !if(lplot.and.tm.gt.tfo)arxy(ix,iy)=tm*1000
          enddo
          enddo
          call RotateHydroX2XX(-phinull
     .    ,nxhxxx,nyhxxx,temx,velx,temxx,velxx)
          do ix=1,nxhy
          do iy=1,nyhy
            temcxx(iz,ntau,ix,iy)=temxx(ix,iy)
            velcxx(1,iz,ntau,ix,iy)=velxx(1,ix,iy)
            velcxx(2,iz,ntau,ix,iy)=velxx(2,ix,iy)
            velcxx(3,iz,ntau,ix,iy)=velxx(3,ix,iy)
            tm=temxx(ix,iy)
            arxy2(ix,iy)=0
            arxy3(ix,iy)=0
            if(lplot.and.tm.gt.tfo)arxy2(ix,iy)=tm*1000
            if(lplot.and.tm.gt.tfo)arxy3(ix,iy)=
     .        sqrt(velxx(1,ix,iy)**2+velxx(2,ix,iy)**2)      
          enddo
          enddo
          !~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
          if(lplot)then
          tau=tauminhy+(ntau-1)*(taumaxhy-tauminhy)/(ntauhy-1)
          if(ish.ge.7)write(txt,'(f8.3)')tau
          call oo2(-4, nxhxx,nyhxx,arxy2,nxhy,xminhy,xmaxhy
     .     ,nyhy,yminhy,ymaxhy,'x   ','y   ','tem xy',txt)
          call oo2(-4, nxhxx,nyhxx,arxy3,nxhy,xminhy,xmaxhy
     .     ,nyhy,yminhy,ymaxhy,'x   ','y   ','vel xy',txt)
          endif
          !~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
        enddo
      enddo
      xmin(1)=xminhy
      xmin(2)=yminhy
      xmin(3)=zminhy
      xmin(4)=tauminhy
      ndim(1)=nxhy
      ndim(2)=nyhy
      ndim(3)=nzhy
      ndim(4)=ntauhy
      del(1)=(xmaxhy-xminhy)/(nxhy-1)
      del(2)=(ymaxhy-yminhy)/(nyhy-1)
      del(3)=(zmaxhy-zminhy)/(nzhy-1)
      del(4)=(taumaxhy-tauminhy)/(ntauhy-1)
      if(ish.ge.7)write(ifmt,'(a,2i5)')'min,max:',(nptlpt+1),nptlbd

      do n=1,nptl
         if(ityptl(n).ge.20.and.ityptl(n).le.29)ityptl(n)=20
         if(ityptl(n).ge.30.and.ityptl(n).le.39)ityptl(n)=30
      enddo

      do n=(nptlpt+1),nptlbd
        iaaptl(n)=n
        if(istptl(n).eq.0.and.ityptl(n).ge.20.and.ityptl(n).le.39)then
          eecobf=eecobf+pptl(4,n)
        endif
      enddo

      nfree=nptlbd-(nptlpt+1)+1
      do nxxn=(nptlpt+1),nptlbd
        nrd=1+rangen()*nfree
        nrd=max(1,nrd)
        nrd=min(nfree,nrd)
        nrd=(nptlpt+1)-1+nrd
        n=iaaptl(nrd)
        do nso=nrd,(nptlpt+1)-1+nfree-1
          iaaptl(nso)=iaaptl(nso+1)
        enddo
        iaaptl((nptlpt+1)-1+nfree)=n
        nfree=nfree-1
      enddo

      eecobfxx=0
      do nxxn=nptlbd,(nptlpt+1),-1
        n=iaaptl(nxxn)
        if(istptl(n).eq.0)then
        eecobfxx=eecobfxx+pptl(4,n)
        if(istptl(n).eq.0.and.ityptl(n).ge.20.and.ityptl(n).le.39)then
        eecobfx=eecobfx+pptl(4,n)
        endif
        endif
      enddo
      if(abs(eecobf-eecobfx).gt.1e-4.and.
     . abs(eecobf-eecobfx).gt.1e-4*eecobf)stop'ERROR 02012012'
      !write(ifmt,*)'+++++++++++++++ eecobf,x,xx',eecobf,eecobfx,eecobfxx

      !write(ifmt,'(6a)')'        n'
      !.,'   eegain','   getIcoEe1hll()','    eerem','   eefrac'


      i888888=0
      g1=0
      g2=0
      n2=0
      eefrac=1
      do nxxn=nptlbd,(nptlpt+1),-1
      n=iaaptl(nxxn)
      nohea=1
      if(istptl(n).eq.0)then
      call idquac(n,idum1,idum2,idum3,jc)
      if(jc(4,1).ne.0.or.jc(4,2).ne.0)nohea=0
      if(jc(5,1).ne.0.or.jc(5,2).ne.0)nohea=0
      endif   
      if(istptl(n).eq.0.and.ityptl(n).ge.20.and.ityptl(n).le.39
     .  .and.nohea.eq.1)then
        if(eefrac.le.0.65)then
           goto 199
        endif
        xo(1)=xorptl(1,n)
        xo(2)=xorptl(2,n)
        xo(3)=xorptl(3,n)
        xo(4)=xorptl(4,n)
        io=iorptl(n)
        if(istptl(io).ne.29)stop'*****  ERROR 31122011  *****'
        ioo=iorptl(io)
        xoo(1)=xorptl(1,ioo)
        xoo(2)=xorptl(2,ioo)
        xoo(3)=xorptl(3,ioo)
        xoo(4)=xorptl(4,ioo)
        do im=1,5
         p(im)=pptl(im,n)
         px(im)=pptl(im,n)
        enddo
        etapxx=getetap(p(1),p(2),p(3))
        rap=99999
        ppl=pptl(4,n)+pptl(3,n)
        pmi=pptl(4,n)-pptl(3,n)
        if(ppl.gt.0.0.and.pmi.gt.0.0)rap= 0.5d0*log(ppl/pmi)
        ptr=p(1)**2+p(2)**2
        if(ptr.gt.0.)ptr=sqrt(ptr)
        v(1)=pptl(1,n)/pptl(4,n)
        v(2)=pptl(2,n)/pptl(4,n)
        v(3)=pptl(3,n)/pptl(4,n)
        gmm=99999
        gmm=1-v(1)**2-v(2)**2-v(3)**2
        if(gmm.gt.0.)gmm=sqrt(gmm)
        if(gmm.gt.0.)gmm=1/gmm
        if(v(1).ne.0.0.or.v(2).ne.0.0)then
          if(ish.ge.6)then
            write(ifch,'(a,i5,f9.2,3x,3f7.2,i7)')
     .       ' ++++++++++++++++++v++',n,ptr,v(1),v(2),v(3),idptl(n)
            write(ifch,'(a,2(4x,4f7.2))')' +++++ori ' ,xoo,xo
          endif
          ihit=0
          nhit=0
          temp1=999.999
          temp2=999.999
          temp3=999.999
          i=9999
          j=9999
          k=9999
          za=xo(3)
          ta=xo(4)
          vz=v(3)  
          z=xo(3)  
          t=xo(4)
          if(z.ge.t)goto 99
          if(z.le.-t)goto 99
          eta= 0.5d0*log((t+z)/(t-z))
          tau=0d0
          tau=(t-z)*(t+z)
          if(tau.gt.0d0)tau=sqrt(tau)
          l=1+(tau-xmin(4))/del(4)
          l=max(l,1)
  77      continue  
          i=999
          j=999
          k=999
          if(l.ge.ndim(4))goto 99
          tau=xmin(4)+(l-1)*del(4)
          !solve (X=dt) :  (ta+X)**2-(za+X*vz)**2=tau**2
          a=(1-vz)*(1+vz)
          if(a.le.0.)goto 99
          if(ta.le.za)goto 99
          b=2*ta-2*za*vz
          c=ta**2-za**2-tau**2
          dt=sqrt((b/2d0/a)**2-c/a)-b/2d0/a
          x=xo(1)+dt*v(1)
          y=xo(2)+dt*v(2)
          z=xo(3)+dt*v(3)
          t=xo(4)+dt
          tauxx=0
          tauxx=(t-z)*(t+z)
          if(tauxx.gt.0.)tauxx=sqrt(tauxx)
          if(abs(tau-tauxx).gt.0.001d0.and.
     .    abs(tau-tauxx).gt.0.001d0*tauxx)
     .    write(ifmt,*)' WARNING jetfluid: tau tauxx = ',tau, tauxx
          i=1+(x-xmin(1))/del(1)
          j=1+(y-xmin(2))/del(2)
          if(i.lt.1.or.i.ge.ndim(1))goto 99
          if(j.lt.1.or.j.ge.ndim(2))goto 99
          if(z.ge.t)goto 99
          if(z.le.-t)goto 99
          eta= 0.5d0*log((t+z)/(t-z))
          k=1+(eta-xmin(3))/del(3)
          if(k.lt.1.or.k.ge.ndim(3))goto 99
          if(p(3).ge.p(4))goto 99
          if(p(3).le.-p(4))goto 99
          rap= 0.5d0*log((p(4)+p(3))/(p(4)-p(3)))
          !~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
          if(ish.ge.6)write(ifch,'(a,4f7.2, 4i3,$)')
     .     ' +++++ xyzt:', x,y,z,t,i-ndim(1)/2-1,j-ndim(2)/2-1,k,l
          !~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
          frac=(eta-xmin(3))/del(3)-(k-1)
          wk(1)=1-frac
          wk(2)=frac
          frac=(x-xmin(1))/del(1)-(i-1)
          wi(1)=1-frac
          wi(2)=frac
          frac=(y-xmin(2))/del(2)-(j-1)
          wj(1)=1-frac
          wj(2)=frac
          temp=0
          velo(1)=0
          velo(2)=0
          velo(3)=0
          do kk=1,2
          do ii=1,2
          do jj=1,2
          temp=temp+wk(kk)*wi(ii)*wj(jj)
     .     *temcxx(k+kk-1,l,i+ii-1,j+jj-1)
          do mm=1,3
          velo(mm)=velo(mm)+wk(kk)*wi(ii)*wj(jj)
     .     *velcxx(mm,k+kk-1,l,i+ii-1,j+jj-1)
          enddo
          enddo
          enddo
          enddo
          temp3=temp2
          temp2=temp1
          temp1=temp
          xsurf=0.
          ysurf=0.
          zsurf=0.
          tsurf=0.
          if(temp.gt.tfo)then
            ihit=1
            nhit=nhit+1
            ix=i
            jx=j
            kx=k
            lx=l
            tempx=temp
            velox(1)=velo(1)
            velox(2)=velo(2)
            velox(3)=velo(3)
            etax=eta
            xsurf=x
            ysurf=y
            zsurf=z
            tsurf=t
          endif   
          if(ish.ge.6)write(ifch,'(3x,3f7.3,i3)')temp3,temp2,temp1,ihit
          l=l+1
          if(temp1.gt.tfo.or.temp2.gt.tfo.or.temp3.gt.tfo)goto 77
          if(ihit.eq.0)goto 98
          !if(nhit*del(4).gt.1.0)goto 95
          
          dxy= (xsurf-xoo(1))**2 + (ysurf-xoo(2))**2
          if(dxy.gt.0.)then
            dxy=sqrt(dxy)
          else
            dxy=0
          endif  
          if(ish.ge.6)write(ifch,'(a,f8.3,a,2f8.3)')' +++++ dxy:',dxy
     .    ,'    velox:',velox(1),velox(2)
          if(dxy.lt.cutdxy)goto 98
          
          etaxx=etax
          gmbj=cosh(etaxx)
          v1xx=velox(1)
          v2xx=velox(2)
          v3xx=velox(3)
          rapxx=0
          if(1d0-v3xx.gt.0d0.and.1d0+v3xx.gt.0d0)
     .     rapxx=0.5d0*log((1d0+v3xx)/(1d0-v3xx))
          v1til=v1xx/gmbj  !in pp-cms = calc frame
          v2til=v2xx/gmbj
          v3til=tanh(etaxx+rapxx)   !=0
          !if(1d0-v1til**2-v2til**2-v3til**2.le.0d0)then
          v3til=tanh(etaxx)  ! ***** ALWAYS, more consistent *****
          !endif
          if(1d0-v1til**2-v2til**2-v3til**2.le.0d0)then
            ngm=0
            do while(1d0-v1til**2-v2til**2-v3til**2.le.0d0
     .      .and.ngm.lt.100)
            ngm=ngm+1
            v1til=v1til*0.9
            v2til=v2til*0.9
            enddo
          endif
          if(v(1)**2+v(2)**2.lt.v1til**2+v2til**2)goto 95
          
          if(ijetfluid.eq.2)goto 99
          
          !~~~~~~~~~~~~~~~~add flow to segments~~~~~~~~~~~~~
          qq(1)=0
          qq(2)=0
          qq(3)=0
          !idxxxx=idptl(n)  !WWWW31
          !call idqua6(n,nq,ns,nc,nb,nt,na)
          call idquaRan(nq,ns,nc,nb,nt,na)
          surat=flratio(quamass(3),quamass(1),tfo)
          if(ish.ge.7)write(ifch,*)'etax=',etax,'   quamass=',quamass
          if(ish.ge.7)write(ifch,*)
     .     '+++++id,nq,na,ns ',idptl(n),nq,na,'     ',ns,nc,nb,nt
          if(nc.ne.0)stop'***** ERROR 04112011b *****'
          if(nb.ne.0)stop'***** ERROR 04112011c *****'
          if(nt.ne.0)stop'***** ERROR 04112011d *****'
          if(na.ne.2.and.na.ne.3)stop'***** ERROR 04112011e *****'
          ic(1)=0
          ic(2)=0
          nsqu=0
          do j=1,na
            rg=rangen()
            suratxx=surat
            do isu=1,nsqu
            suratxx=suratxx*0.5
            enddo
            ifl=1+rg*(2+suratxx)
            ifl=max(1,ifl)
            ifl=min(3,ifl)
            if(ifl.eq.3)nsqu=nsqu+1
            if(ish.ge.7)write(ifch,*)'+++++ifl',ifl,'   ',rg
            isg=1
            if(na.eq.3.and.nq.lt.0)isg=-1
            if(na.eq.2.and.j.eq.2)isg=-1
            if(isg.eq. 1)ic(1)=ic(1)+10**(nflav-ifl)
            if(isg.eq.-1)ic(2)=ic(2)+10**(nflav-ifl)
          enddo
          id=idtra(ic,0,0,0)
          if(ish.ge.7)write(ifch,*)'+++++ic,id',ic,'   ',id
          idxxx=id
          ia=abs(id)
          if(mod(id,10).eq.0.and.rangen().lt.0.5)id=id/ia*(ia+1)
          if(ia/10.eq.123.and.rangen().lt.0.33333)id=id/ia*2130
          call idmass(id,am)
          p(5)=am
          if(ish.ge.7)write(ifch,*)'+++++id, am ',id, am
          !~~~~~~~~
          call idqua(n,nq,ns,na)
          if(ish.ge.7)
     .    write(ifch,*)'+++++id,nq,na,ns ',id,nq,na,'     ',ns
          !if(abs(id).eq.3331)print*,'WWWW31',idxxxx,id
          nqx=nq
          nsx=ns
          nax=na
          nn=0
  76      continue 
          nn=nn+1 
          if(ns.ne.0)then
           am=quamass(3)
           nq=nq-ns/iabs(ns)
           ns=ns-ns/iabs(ns)
           na=na-1
          elseif(nq.ne.0)then
           am=quamass(1)
           nq=nq-nq/iabs(nq)
           na=na-1
          elseif(na.ne.0)then
           am=quamass(1)
           na=na-1
          endif
          cot=0
          x=hynRanBoseFermi(am,cot,1,tfo,1) 
          pp=x*tfo
          e=sqrt(pp**2+am**2)
          u(3)=2.*rangen()-1.
          angle=2.*pi*rangen()
          u(1)=sqrt(1.-u(3)**2)*cos(angle)
          u(2)=sqrt(1.-u(3)**2)*sin(angle)
          q(1)=pp*u(1)
          q(2)=pp*u(2)
          q(3)=pp*u(3)
          q(4)=e
          if(ish.ge.7)
     .     write(ifch,'(f7.3,3x,4f7.3,3x,3f7.3,4x,3i3,2x,3i3)')
     .     am,q,velox,nqx,nsx,nax ,nq,ns,na
          gm=1./sqrt(1-velox(1)**2-velox(2)**2)
          gmtil=1d0/sqrt(1d0-v1til**2-v2til**2-v3til**2)
          u1til=gmtil*v1til
          u2til=gmtil*v2til
          u3til=gmtil*v3til
          u4til=gmtil
          call utlob3(-1, u1til , u2til , u3til , u4til ,1e0
     .             , q(1), q(2), q(3), q(4))
          qq(1)=qq(1)+q(1)
          qq(2)=qq(2)+q(2)
          qq(3)=qq(3)+q(3)
          !~~~~~~~~~~~~~~~~~~~
          if(.not.(q(3).le.0..or.q(3).ge.0.))then
           write(ifch,*)'**** NaN catch q3 ****',n,q(3),gmtil,v3til**2
     .     ,v1til**2+v2til**2+v3til**2    
           stop'*****  ERROR 14092012b  *****'
          endif
          !~~~~~~~~~~~~~~~~~~~
          if(ish.ge.7)write(ifch,'(10x,4f7.3)')q
          if(na.ne.0)                                 goto 76
          !----------------------------------------------------------->
          if(nq.ne.0)stop'***** ERROR 06072011 *****'
          if(ns.ne.0)stop'***** ERROR 06072011c *****'
          if(nn.ne.2.and.nn.ne.3)stop'***** ERROR 06072011b *****'
          if(ish.ge.7)write(ifch,'(a,5f7.3)')'  qq:',qq
          if(ish.ge.7)write(ifch,'(a,5f7.3)')'  p: ',px
          !~~~~~
          fcs=(1-ptfra*sqrt(2.)/ptr)
          if(fcs.gt.0)then
            p(1)=p(1)*fcs
            p(2)=p(2)*fcs
          endif
          !~~~~~
          p(1)=p(1)+qq(1)
          p(2)=p(2)+qq(2)
          if(ish.ge.7)write(ifch,'(a,5f7.3)')'  p: ',p
          amt=sqrt(p(1)**2+p(2)**2+p(5)**2)
          p(3)=amt*sinh(rap)
          p(4)=amt*cosh(rap)
          !p(3)=p(3)+qq(3)
          !p(4)=sqrt(amt**2+p(3)**2)
          !~~~~~~~~~~~~~~~~~~~
          if(.not.(p(4).le.0..or.p(4).ge.0.))then
           write(ifch,*)'**** NaN catch p4 ****',n,p(4),amt,p(3)
           stop'*****  ERROR 14092012  *****'
          endif
          !~~~~~~~~~~~~~~~~~~~
          etap=getetap(p(1),p(2),p(3))
          if(ish.ge.7)write(ifch,'(a,5f7.3)')'  p: ',p
          !????????????????????????????????????????????????????????????
          !pt99=sqrt(p(1)**2+p(2)**2)
          !if(pt99.gt.1..and.pt99.lt.2.)then
          !write(ifch,'(4f8.3,2f9.3)')px(1),px(2),px(3),px(4),etapxx,v3xx
          !write(ifch,'(4f8.3,2f9.3)')  p(1),p(2),p(3),p(4),etap,rapxx
          !write(ifch,*)'jetpeak',n,amt,amt**2+p(3)**2
          !endif
          !????????????????????????????????????????????????????????????
          !~~~~~~~~~~~~~~~~~~~~~
          eegainxx=eegain+p(4)-px(4)
          eeremxx=getIcoEe1hll()-eegainxx
          eefracxx=1
          if(getIcoEe1hll().gt.0.)eefracxx=eeremxx/getIcoEe1hll()          
          if(eefracxx.le.0.65)goto 199
          !~~~~~~~~~~~~~~~~~~~~~
          if(ityptl(n).ge.20.and.ityptl(n).le.29)ityptl(n)=29
          if(ityptl(n).ge.30.and.ityptl(n).le.39)ityptl(n)=39
          idptl(n)=id
          do im=1,5
           pptl(im,n)= p(im)
          enddo
          xorptl(1,n)=xsurf
          xorptl(2,n)=ysurf
          xorptl(3,n)=zsurf
          xorptl(4,n)=tsurf
          iorptl(n)=0
          if(ish.ge.5) !.and.abs(rap).lt.0.5.and.abs(px(5)).gt.0.2)
     .     write(ifch,'(1x,a,i6,4f6.3,4x,4f5.1,4x,2f7.2)')
     .    ' OK   ',n,temp1,temp2,temp3,tfo,xo,ptr,gmm
          !~~~~~~~~~~~~~~~~~~~
          !if(nax.eq.3.and.nsqu.eq.1)then 
          !write(ifch,*)'++squ+++',n,'  ',ic,'   ',idxxx,id
          !call alist('&',n,n)
          !endif 
          !~~~~~~~~~~~~~~~~~~~
        else
          if(ish.ge.5) !.and.abs(rap).lt.0.5.and.abs(px(5)).gt.0.2)
     .    write(ifch,*)('-',m=1,80),' skip v = ',v(1),v(2)
        endif
        g1=g1+p(4)-px(4)

        goto 97
                
 99     continue   
        if(ish.ge.5) !.and.abs(rap).lt.0.5.and.abs(px(5)).gt.0.2)
     .  write(ifch,'(1x,a,i6,3x,8i3,3x,4f5.1,3x,2f7.2,a,5f7.1)')
     .   ' SKIP ',n
     .  ,i,ndim(1),j,ndim(2),k,ndim(3),l,ndim(4),xo,ptr,gmm
     .  ,' ',px(1),px(2),px(3),px(4),px(5)
        goto 97
        
 98     continue
        if(ish.ge.5) !.and.abs(rap).lt.0.5.and.abs(px(5)).gt.0.2)
     .  write(ifch,'(1x,a,i6,4f6.3,4x,4f5.1,4x,2f7.2,a,5f7.1)')
     .  ' NOHIT',n,temp1,temp2,temp3,tfo,xo,ptr,gmm
     .  ,' ',px(1),px(2),px(3),px(4),px(5)
        goto 97

 95     continue   
        if(ish.ge.5) !.and.abs(rap).lt.0.5.and.abs(px(5)).gt.0.2)
     .  write(ifch,'(1x,a,i6,3x,8i3,3x,4f5.1,3x,2f7.2,a,5f7.1)')
     .   ' TO FLUID ',n
     .  ,i,ndim(1),j,ndim(2),k,ndim(3),l,ndim(4),xo,ptr,gmm
     .  ,' ',px(1),px(2),px(3),px(4),px(5)
        eegain=eegain-px(4)
        istptl(n)=7
        g2=g2-px(4)
        n2=n2+1
        goto 97

 97     continue
  
        eegain=eegain+p(4)-px(4)
        eerem=getIcoEe1hll()-eegain
        eefrac=1
        if(getIcoEe1hll().gt.0.)eefrac=eerem/getIcoEe1hll()
        !write(ifmt,'(i9,3f9.0,f9.5)')n,eegain,getIcoEe1hll(),eerem,eefrac
        
      endif
      enddo

 199  continue
      eespec=0
      do n=1,nptlpt
        if(istptl(n).eq.0)eespec=eespec+pptl(4,n)
      enddo
      eecotot=0
      eecoafxx=0
      do n=(nptlpt+1),nptlbd
        if(istptl(n).eq.0)then
        eecotot=eecotot+pptl(4,n)
        if(.not.(eecotot.le.0..or.eecotot.ge.0.))then
         write(ifch,*)'**** NaN catch eecotot ****',n,eecotot,pptl(4,n)
         stop'*****  ERROR 11072012  *****'
        endif
        if(istptl(n).eq.0.and.ityptl(n).ge.20.and.ityptl(n).le.39)then
        eecoafxx=eecoafxx+pptl(4,n)
        endif
        endif
      enddo
      !write(ifmt,*)'+++++++ eecoafxx,eegain',eecoafxx,eegain
      eetot=eespec+eecotot+eefrac*getIcoEe1hll()
      !if(mod(nrevt+1,modsho).eq.0) 
      write(ifmt,'(a,f8.0,a,f6.3,a,f8.0,a,f8.0,a,2f8.0,i4,a,i3)')
     . 'jetfluid Eflu',getIcoEe1hll(),' f',eefrac
     . ,' Ecoro',eespec+eecotot
     . ,' Etot',eetot,' g',g1,g2,n2   !,' ctr',ikoevt (ONLY pp)

      call utprix('jetflu',ish,ishini,4)
      end
      
          !before adding qq
          !call utlob3(1, 0.,0., sinh(etax) ,cosh(etax),1.
          !.             , p(1), p(2), p(3), p(4))
          !if(ish.ge.7)write(ifch,'(a,5f7.3)')'  p:',p
          !after adding qq
          !call utlob3(-1, 0.,0., sinh(etax) ,cosh(etax),1.
          !.             , p(1), p(2), p(3), p(4))
          !if(ish.ge.7)write(ifch,'(a,5f7.3)')'  p:',p

c----------------------------------------------------------------------------
      function getetap(p1,p2,p3)
      pt=sqrt(p2**2+p1**2)
      x=0
      if(p3.ne.0..and.pt.ne.0.)x=sign(1.,p3)*
     *    alog((sqrt(p3**2+pt**2)+abs(p3))/pt)
      getetap=x
      end

c----------------------------------------------------------------------------
      subroutine idquaRan(nq,ns,nc,nb,nt,na)
c----------------------------------------------------------------------------
#include "aaa.h"
      nq=0
      ns=0
      nc=0
      nb=0
      nt=0
      na=2
      if(rangen().lt.fludiq)then
        na=3
        nq=3
        if(rangen().lt.0.5)nq=-3
      endif
      end

c----------------------------------------------------------------------------
      function flratio(am1,am2,tem)
c----------------------------------------------------------------------------
#include "ho.h"
      real ar(2)
      real qua(2)
      ar(1)=am1
      ar(2)=am2
      call gaulag(xlag,wlag,nlag,0.)
      do n=1,nlag
      wlag(n)=wlag(n)*exp(xlag(n))
      enddo
      do i=1,2
        a=ar(i)/tem
        fsum=0
        do n=1,nlag
          x=xlag(n)  
          fsum=fsum+wlag(n)*x**2 /(exp(sqrt(x**2+a**2))+1)
        enddo
        qua(i)=fsum
      enddo
      flratio=qua(1)/qua(2)
      end

c----------------------------------------------------------------------------
          subroutine RotateHydroX2XX(phi
     .     ,nxhxxx,nyhxxx,temx,velx,temxx,velxx)
c----------------------------------------------------------------------------
#include "aaa.h"
#include "ho.h"
      dimension velx(3,nxhxxx,nyhxxx),temx(nxhxxx,nyhxxx)
      dimension velxx(3,nxhxxx,nyhxxx),temxx(nxhxxx,nyhxxx)
      real gi(2),gj(2),aa(2,2)
      delx=(xmaxhy-xminhy)/(nxhy-1)
      dely=(ymaxhy-yminhy)/(nyhy-1)
      do ix=1,nxhy
      do iy=1,nyhy
        xp= xminhy  +(ix-1)  *delx
        yp= yminhy  +(iy-1)  *dely
        xv=  cos(phi)*xp - sin(phi)*yp
        yv=  sin(phi)*xp + cos(phi)*yp
        i=1+(xv-xminhy)/delx
        i=max(1,i)
        i=min(nxhy-1,i)
        xi=xminhy+(i-1)*delx
        frac=(xv-xi)/delx
        gi(1)=1-frac
        gi(2)=frac
        j=1+(yv-yminhy)/dely
        j=max(1,j)
        j=min(nyhy-1,j)
        yj=yminhy+(j-1)*dely
        frac=(yv-yj)/dely
        gj(1)=1-frac
        gj(2)=frac
        temxx(ix,iy)=0
        do i2=1,3
        velxx(i2,ix,iy)=0
        enddo
        dlta=-0.001
        if(gi(1).ge.dlta .and. gi(2).ge.dlta .and.
     .     gj(1).ge.dlta .and. gj(2).ge.dlta       )then
         do n=1,2
         do m=1,2
          temxx(ix,iy)=temxx(ix,iy)
     .                        +gi(n)*gj(m)*temx(i+n-1,j+m-1)
          do i2=1,3
          velxx(i2,ix,iy)=velxx(i2,ix,iy)
     .                     +gi(n)*gj(m)*velx(i2,i+n-1,j+m-1)
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
        a=a+aa(1,m)*velxx(m,ix,iy)
        b=b+aa(2,m)*velxx(m,ix,iy)
        enddo
        velxx(1,ix,iy)=a
        velxx(2,ix,iy)=b
      enddo
      enddo
      end
      
      !  write(ifch,*)'CHECK',i+n-1,j+m-1
      !.  ,gi(n)*gj(m),temx(i+n-1,j+m-1)
      !.            ,velx(1,i+n-1,j+m-1),velx(2,i+n-1,j+m-1)
      
      !  write(ifch,*)'CHECK',ix,iy,temxx(ix,iy)
      !. ,(velxx(i2,ix,iy),i2=1,2)
       
c-----------------------------------------------------------------------
      subroutine Pi2Hyt(xxx,eps,vv,iout)
c-----------------------------------------------------------------------
c      use ArrayModule, only: getd
c      use hocoModule, only: epsc, velc
c      use hocoModule, only: velc
      double precision getHydynEpsc, getHydynVelc
#include "aaa.h"
#include "ho.h"
      real xxx(4),vv(3)
      real wi(3),wj(3),wn(3),wm(3)
      x=xxx(1)
      y=xxx(2)
      z=xxx(3)
      t=xxx(4)
      iout=1
      eps=0
      vv(1)=0
      vv(2)=0
      vv(3)=0
      if(t+z.eq.0.)return
      if(t-z.eq.0.)return
      et=0.5*log((t+z)/(t-z))
      ta=sqrt((t+z)*(t-z))
      delx=(xmaxhy-xminhy)/(nxhy-1)
      dely=(ymaxhy-yminhy)/(nyhy-1)
      deleta=(zmaxhy-zminhy)/(nzhy-1)
      dltau=(taumaxhy-tauminhy)/(ntauhy-1)
      xx=(x-xminhy)/delx +1      
      yy=(y-yminhy)/dely +1
      ee=(et-zminhy)/deleta +1
      tt=(ta-tauminhy)/dltau +1
      nx=int(xx)
      ny=int(yy)
      ne=int(ee)
      nt=int(tt)
      if(nx.ge.1.and.ny.ge.1.and.ne.ge.1.and.nt.ge.1
     ..and.nx.le.nxhy-1.and.ny.le.nyhy-1
     ..and.ne.le.nzhy-1.and.nt.le.ntauhy-1)then
        frac=xx-nx
        wi(1)=1-frac
        wi(2)=frac
        frac=yy-ny
        wj(1)=1-frac
        wj(2)=frac
        frac=ee-ne
        wn(1)=1-frac
        wn(2)=frac
        frac=tt-nt
        wm(1)=1-frac
        wm(2)=frac
        do i=1,2
         do j=1,2
          do n=1,2
           do m=1,2
            w=wi(i)*wj(j)*wn(n)*wm(m)
            eps  =eps  + getHydynEpsc(  ne+n-1,nt+m-1,nx+i-1,ny+j-1) * w
            vv(1)=vv(1)+ getHydynVelc(1,ne+n-1,nt+m-1,nx+i-1,ny+j-1) * w
            vv(2)=vv(2)+ getHydynVelc(2,ne+n-1,nt+m-1,nx+i-1,ny+j-1) * w
            vv(3)=vv(3)+ getHydynVelc(3,ne+n-1,nt+m-1,nx+i-1,ny+j-1) * w
           enddo
          enddo
         enddo
        enddo
        iout=0
      endif
      end

