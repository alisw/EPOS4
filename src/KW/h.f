C
C  This file is part of EPOS4
C  Copyright (C) 2022 research institutions and authors (See CREDITS file)
C  This file is distributed under the terms of the GNU General Public License version 3 or later
C  (See COPYING file for the text of the licence)
C

c-----------------------------------------------------------------------
      subroutine itsyminit
c-----------------------------------------------------------------------
      !symmetric tensor to vector mapping
      !       index-0-1-2-3-              1-2-3--4-
      !   map   0-  1,2,4,7,                        -0        
      ! tensor  1-  2,3,5,8,     switch   3,5,8, 2  -1
      !   to    2-  4,5,6,9,    indices   5,6,9, 4  -2
      ! vector  3-  7,8,9,10    0 and 4   8,9,10,7  -3
      !                                   2,4,7, 1  -4
      !---------------------------------------------
      integer itsym(4,4) 
      data itsym /
     .              3,5,8,  2 
     .,             5,6,9,  4
     .,             8,9,10, 7
     .,             2,4,7,  1 /
      common/citsy/itsy(4,4)
      !---------------------------------------------
      ! itsy set to equal itsym below
      !   usage : tensor^{mu,nu} =  vec( itsy(mu,nu) )
      !--------------------------------------------- 
      do i=1,4
        do j=1,4
          itsy(i,j)=itsym(i,j)
        enddo
      enddo
      end


#ifndef __BS__

c#######################################################################
!-----------------------------------------------------------------------
!> @brief
!> removes final TLC partons originating from 25 (=Born) partons 
!> for 25 partons with pt > ptrmin and E > getJetEmin()
!>
!> To add in optns:
!>   !---------------------------------------------------------------------------
!>   print * 2      !print info to check file (see "list before fragmentation")
!>   set jet 1      !activate jet stuff
!>   set jetcheck 1 !print additional (jet related) info to checkfile
!>                  !find in the check file:
!>                  !"removeTLCfromBorn (1)" contains relevant piece of particle list
!>                  !"removeTLCfromBorn (2)" piece of list after removing TLC partons
!>                  !"list before fragmentation (2)"  --> new complete list
!>                  !Important:
!>                  !Our 21 partons (TLC) have at this stage still ist=20; 
!>                  !  among these ist=20 partons one identifies those coming from 
!>                  !  Born by the ifr1 (4th column) values bigger than 10
!>                  !If the Jet procedure is triggered, there are only 2 ifr1>10
!>                  !  partons,  since the TL cascade was suppressed
!>                  !These two ifr1>10 partons are actually replaced by zero momentum
!>                  !  partons (same flavor) which is needed to make color connections
!>   set jetemin 50 !min jet energy (Born CMS) to trigger Jet procedure
!>   set ptrmin 30  !min jet pt (Lab sys) to trigger Jet procedure
!>   !---------------------------------------------------------------------------
!
!> @param[in] idry to do a dry run (1) or the real one (2) 
!> @param[in] icheck to produce output to check file (1) or not (0) 
!-----------------------------------------------------------------------
      subroutine removeTLCfromBorn(idry,icheck)
c#######################################################################
#include "aaa.h"
      integer :: idry,icheck
      real p1,p2,p3,p4,p5
      if(idry.ne.1.and.idry.ne.2)stop '####### ERROR 231026d #######'
      if(idry.eq.1.and.icheck.eq.0)return
      if(idry.eq.2)call iniParticleHighpt25() !initialize Highpt25 list
c header
      if(icheck.eq.1)
     .write(ifch,'(/S1x,72a1/1x,a,i1,a,35a1/1x,72a1/)')
     .('#',k=1,72),'############  removeTLCfromBorn (',idry,')  ',
     .('#',k=1,35),('#',k=1,72)
c find 25 partons in the particle list
      nk1=maproj+matarg+1 
      do while(nk1.le.nptl)
        call getistptl(nk1,ist1)
        if(ist1.eq.25)then
          nk2=nk1+1
          call getistptl(nk2,ist2)
          if(ist2.ne.25)stop '####### ERROR 231026a #######'
          nk1x25=nk1 !25 parton
          nk2x25=nk2 !25 parton
          nk1=nk2+1 !first 21 parton
          call getidptl(nk1,id1) !should not be 9 (gluon)
          if(id1.eq.9)stop '####### ERROR 231026b #######'
          nk2=nk1+1
          call getistptl(nk2,ist2)
          do while(ist2.eq.20)
            nk2=nk2+1
            call getistptl(nk2,ist2)
          enddo
          nk2=nk2-1
          call getidptl(nk2,id2) !should not be 9 (gluon)
          if(id2.eq.9)stop '####### ERROR 231026c #######'
          !we have 21 partons from nk1 to nk2 (Born partons + SLC partons)
          nfromBorn=0
          do nk=nk1,nk2
            call getifrptl(nk,ifr1,ifr2) 
            if(ifr1.gt.10)then !Born partons (emitted from 25)
              nfromBorn=nfromBorn+1
            endif
          enddo 
          iOK=0
          if(nfromBorn.eq.2)then
            call getidptl(nk1x25,id1x25)
            call getidptl(nk2x25,id2x25)
            call getpptl(nk1x25,p1,p2,p3,p4,p5)
            pt1=p1**2+p2**2
            call getpptl(nk2x25,p1,p2,p3,p4,p5)
            pt2=p1**2+p2**2
            if(pt1.gt.ptrmin**2.or.pt2.gt.ptrmin**2)iOK=1
          endif
          if(iOK.eq.1)then
            !=====================================================================
            if(idry.eq.2)call addParticleHighpt25(nk1x25)
            if(idry.eq.2)call addParticleHighpt25(nk2x25)
            ! print list of partons
            if(icheck.eq.1)
     .      write(ifch,'(1x,97a1/1x,4x,a,9x,a/1x,97a1)')('-',k=1,97),
     .      'ior    jor      i   ifr1   ifr2       id ist ity',
     .      'p1         p2         p3         p4',('-',k=1,97)
            do k=1,4
            enddo
            irem=0
            i=nk1-2
            do while(i.le.nk2)
              call getifrptl(i,ifr1,ifr2) 
              idelete=0
              if(i.ge.nk1)then !21 partons
                if(ifr1.gt.10)then !only Born partons (emitted from 25)
                  if(idry.eq.2)idelete=1
                endif
              endif
              call getiorptl(i,ior)    
              call getjorptl(i,jor)
              call getidptl (i,id)
              call getistptl(i,ist)
              call getityptl(i,ity)
              call getpptl(i,p1,p2,p3,p4,p5)
              if(idelete.eq.1)then
                write(ifmt,'(a,5i7,i10,2i3,2x,e10.3)')
     .          'removeTLCfromBorn ',
     .          ior,jor,i,ifr1,ifr2,id,ist,ity,sqrt(p1**2+p2**2)
                call utreplaZero(i)  !replaced by zero momentum parton 
                call setidptl (i,id)   
                call setistptl(i,20) 
                call setiorptl(i,ior  ) 
                call getiorptl(i,ior)    
                call getjorptl(i,jor)
                call getifrptl(i,ifr1,ifr2) 
                call getidptl (i,id)
                call getistptl(i,ist)
                call getityptl(i,ity)
                call getpptl(i,p1,p2,p3,p4,p5)
              endif
              if(icheck.eq.1)
     .        write(ifch,'(1x,5i7,i10,2i3,2x,4(e10.3,1x))')
     .        ior,jor,i,ifr1,ifr2,id,ist,ity,p1,p2,p3,p4
              i=i+1
            enddo ! i loop 
            !=====================================================================
          endif
          nk1=nk2+1
        else
          nk1=nk1+1
        endif
      enddo
      end

c-----------------------------------------------------------------------
      subroutine hllex(nfr)
c-----------------------------------------------------------------------
c      use arrayModule, only: newDoubleArray, deleteDoubleArray
c      use hocoModule, only: epsc, sigc, velc, barc,
c      use hocoModule, only: sigc, velc, barc,
c     &     sigcShape, velcShape, barcShape 
#include "aaa.h"
#include "ho.h"
      common /jcen/jcentr,jmxcentr
      common/cwmaxIcoE/wmaxIcoE
      common/ciemuc/iemuc

      real getHydynTauhy
      integer getAccumNrevt

      if(ihlle.ne.0.and.ireadhyt.eq.0)then
        if(nfr.eq.0)then
          call InitializeHyperbola
          if(wmaxIcoE.gt.0.)then
            if(iemuc.eq.0)then
             call memo(1,'create emuc object;')
             call createemuc(10,nzhy,ntauhxx,nxhy,nyhy)
             call memo(2,';')
             iemuc=1
            endif
            write(ifmt,'(a,i3,a)')'enter hlle; nfr=',nfr,' ... '
            call clop(3)
            call hlle(1)
            call useHydroTable
            write(ifmt,'(a,f7.2)')'exit hlle, taumax=',
     .           getHydynTauhy(ntauhy)
            !call skipRandomNumbers
            call clop(3)
            call WriteHydroTable
          else
            write(ifmt,'(a,i3,a)')'skip hlle; eps=0'
            ntauhy=0
          endif
          call ManiHydroTable
          call clop(3)
          call setHo(1,1,0)
          call getJcentr
          call fros(1)
          call fros(2)
          jcentr_save=jcentr
          if(getAccumNrevt().gt.nevent-nfreeze)then
            call setHo(0,2,1)
            do jcentr=1,jmxcentr
              call fros(1)
            enddo
            call fros(4)
          endif
          jcentr=jcentr_save
          call fros(5)
          call epof(1)
        endif
          !if(ntauhy.eq.0.and.tauup.gt.0.05)then
          !  do nn=nptlpt+1,nptlbd
          !    if(istptl(nn).eq.5.or.istptl(nn).eq.7)istptl(nn)=0
          !  enddo
          !endif
        !call checkengy('before jetfluid ')
        call jetfluid
        !call checkengy('after jetfluid  ')
        nptlBeforeFO=nptl
        call EpoMiF
        call epof(2)
        if(ish.ge.2)
     .  call aalist('list after core freeze out&',nptlBeforeFO+1,nptl)
        !if(mod(nrevt+1,modsho).eq.0)
        call checkengy('after freeze out')
        call utrescxx(iret,0,104)  !after freeze out, not any more before call afinal
        call checkengy('after utrescxx  ')
        do i=1,nptl
        if(ityptl(i).ne.61.and.idptl(i).eq.331)call utphiBW(i)  !BW mass smearing
        if(ityptl(i).eq.61.and.idptl(i).eq.333)call utphiBW(i)  !BW mass smearing
        enddo 
        call decayall(1,0)
      endif
      if(idecay.eq.2)then
        do i=1,nptl
          if(istptl(i).eq.7)then
            istptl(i)=0
            ityptl(i)=60
          endif 
        enddo 
        call decayall(1,0)
      endif
      if(ntauhy.ne.0)then
        call memo(1,'destroy dmass objects ;')
        call dmass_c0destroy()
        call dmass_c1destroy()
        call dmass_c2destroy()
        call dmass_c3destroy()
        call dmass_i0destroy()
        call dmass_i1destroy()
        call dmass_i2destroy()
        call dmass_i3destroy()
        call memo(2,';')
      endif
!call memo(0,'exit hllex;')

      end subroutine hllex


c-----------------------------------------------------------------------
      subroutine hllexx
c-----------------------------------------------------------------------
#include "aaa.h"
      if(iorsdf.eq.5.and.ispherio+ihlle.eq.0)then
        !write(ifmt,'(a)')'core but no hydro, change ist=5 to 0'
        do i=1,nptl
          call getistptl(i,istx)
          if(istx.gt.3.and.istx.le.7)then
            if(istx.ne.5)stop'ERROR 19032022'
            call setistptl(i,0)
            !print*,'TESThllexx',idptl(i),istptl(i),istx
          endif 
        enddo
      endif
      end subroutine hllexx


c-----------------------------------------------------------------------
      subroutine hlle(mode)
c-----------------------------------------------------------------------
      implicit none
      integer,parameter::mxjerr=30
      integer,parameter::mxho=10
      integer icotabr
      integer nfnio
      integer       iopcnt
      integer ifmt
      real         tauzer,tauone,tautwo,tauthree,tauup
      integer ioeos, iozerof
      real etaos, zetaos
      integer ntaumx
      character*500 fnio
      real xminico,xmaxico,yminico,ymaxico,zminico,zmaxico
      integer,parameter::nxicomax=161
      integer,parameter::nyicomax=161
      integer,parameter::nzicomax=45
      integer nxico,nyico,nzico
      real getIcoIcoE, getIcoIcoV, getIcoIcoF
      integer,parameter::ntauhxx=120

      integer nxhy,nyhy,nzhy,ntauhy
      real oEeos
      real dtauhy
      real epsfin
      integer ienvar
      integer ifaahlle,ifathlle,ifazhlle
      real tfo
      real xminhy,xmaxhy,yminhy,ymaxhy,zminhy,zmaxhy
      real tfrout
      integer i, ierrchk, istep, istepx, itest, ix, iy, iz, 
     . j, jx, jy, jz, k, memory, mmsteps, mode,
     . ncnthl, ncnthlf, ntau, ntsteps, nxihy, nyihy, nzihy
      real cputime, del_tau, delx, dely, delz, eeitau, 
     . eejtau, eektau, eetau2, eetauaa, eetauxx, egain,  
     . eistxx, emx, eps, esollxx, fd, fs, fu, gamma, ggi, pi00, pi03, 
     . pi33, ratioee, tau, tau2, tau3, tidisu, u0, u3, zero, zz
      real eflow, eechk, gettimehlle
      real dummy1,dummy2,enfo

      integer getAccumJerr
      integer getCoreIcotabr
      real    getEospOeEos
      integer getEventIopcnt
      integer getFilesIfmt
      integer getFilesNfnio
      real    getHydynDtauhy
      real    getHydynEpsfin
      integer getHydynIenvar
      integer getHydynIfaahlle
      integer getHydynIfathlle
      integer getHydynIfazhlle
      integer getHydynNtaumx
      real    getHydynTfo
      real    getHydynTfrout
      real    getHydynTauhy
      real    getCcotiCoti
      real    getRescEtaos
      real    getRescZetaos
      real    getIcoEe1ico
      real    getIcoEistico
      real    getHydynEetau
      real    getHydynEetau2
      real    getHydynEmx
      real    getHydynRatioeex
c ------------------------------------------------------------

      double precision tau0,tau_max,dtau,frac,nb,nq,ns,tau1
      double precision e,vx,vy,vz,x,y,z,xj,yj,zj,efin,v2
      double precision zzz
      double precision ecc,vxcc,vycc,vzcc,nbcc,nqcc,nscc,xcc,ycc,zcc
      double precision getxhlle, getyhlle, getzhlle,wx(2),wy(2),wz(2)
      double precision delxico,delyico,delzico
      double precision xmnhy,xmxhy,ymnhy,ymxhy,zmnhy,zmxhy
      double precision etaS,zetaS
      double precision exxx, efull,etotsurf,nbsurf
      double precision iutime(5),tiu3,tiu4,tidi
      data ncnthl/0/
      save ncnthl
      data ncnthlf/0/
      save ncnthlf
      integer ncnthlle
      data ncnthlle/0/
      save ncnthlle

      real IcoEArray
      dimension IcoEArray(nxicomax,nyicomax,nzicomax)
      
      character onechar

c     read Fortran common variable values from Cpp

      call createCore()
      call createEosp()
      call createEvent()      
      call createResc()

      call getRescTau(tauzer, tauone, tautwo, tauthree, tauup)
      etaos    = getRescEtaos()
      zetaos   = getRescZetaos()
      icotabr  = getCoreIcotabr()
      oEeos    = getEospOeEos()
      iopcnt   = getEventIopcnt()
      ifmt     = getFilesIfmt()
      nfnio    = getFilesNfnio()
      ntaumx   = getHydynNtaumx()
      dtauhy   = getHydynDtauhy()
      tfrout   = getHydynTfrout()
      epsfin   = getHydynEpsfin()
      ienvar   = getHydynIenvar()
      ifaahlle = getHydynIfaahlle()
      ifathlle = getHydynIfathlle()
      ifazhlle = getHydynIfazhlle()
      tfo      = getHydynTfo()

      call getHY(dummy1,dummy2,enfo)
      if(epsfin.gt.0.5*enfo)stop'##### ERROR 03112024 #####'
      
      call getHydynIo(ioeos, iozerof)
      call getIcoBound(
     .     xminico, xmaxico,
     .     yminico, ymaxico,
     .     zminico, zmaxico
     . )
      call getIcoDim(nxico,nyico,nzico)
      call getHydynDim(nxhy,nyhy,nzhy,ntauhy)
      call getHydynBound(xminhy,xmaxhy,yminhy,ymaxhy,zminhy,zmaxhy)

c copy the C++ Array<char> *fnio to the Fortran array fnio
      do i=1,getFilesNfnio()
         call getFilesFnio(i, onechar)
         fnio(i:i)=onechar
      end do
 
      ncnthlle=ncnthlle+1
      call checkmemory(memory)
      write(ifmt,'(a,i8)')'start hlle;  memory',memory
      !call redirecterrorshlle('hydro_messages.txt'//CHAR(0))

      tfo=tfrout
      call setHydynTfo(tfo)

      if(mode.eq.-1)then
      call getEosHlle
      oEeos   = getEospOeEos()

      call destroyhlle()
      return
      endif

      call checktime('start hlle;') 

      itest=0            !testing if 1
      call actienvar(ienvar)
      efin=epsfin         !min edensity
      tau0 = tauzer
      tau1 = tauone
      tau2 = tautwo
      tau3 = tauthree
      tau_max = tauzer+tauup
      dtau = dtauhy / ifathlle
      xmnhy= xminhy  *ifaahlle   !hydro boundaries
      xmxhy= xmaxhy  *ifaahlle
      ymnhy= yminhy  *ifaahlle
      ymxhy= ymaxhy  *ifaahlle
      zmnhy= zminhy  *ifaahlle
      zmxhy= zmaxhy  *ifaahlle
      nxihy= (nxhy-1)  *ifaahlle*ifathlle +1  !hydro grid
      nyihy= (nyhy-1)  *ifaahlle*ifathlle +1
      nzihy= (nzhy-1)  *ifaahlle*ifazhlle +1
      delx=(xmxhy-xmnhy)/(nxihy-1)       !hydro steps
      dely=(ymxhy-ymnhy)/(nyihy-1)
      delz=(zmxhy-zmnhy)/(nzihy-1)
      write(ifmt,*)'cell size:',delx,' * ',dely,' * ',delz
      write(ifmt,*)'initial tau:',tau0,'   tau step:',dtau
      call setHydynMmx(1+nint((xminhy-xmnhy)/delx))     !for transfer
      call setHydynMmy(1+nint((yminhy-ymnhy)/dely))
      call setHydynMmz(1+nint((zminhy-zmnhy)/delz))

      call getEosHlle
      oEeos   = getEospOeEos()

      etaS=etaos
      zetaS=zetaos
      write(ifmt,'(a,2f10.5)')'set viscosity coeffs etaS and zetaS:'
     .,etaS,zetaS
      call inittrcoeff(etaS,zetaS)

      call memo(1,'create fluid;')
      call clop(3)
      call initfluidhlle(nxihy+4,nyihy+4,nzihy+4
     .                   ,xmnhy-2*delx,xmxhy+2*delx
     .                   ,ymnhy-2*dely,ymxhy+2*dely
     .                   ,zmnhy-2*delz,zmxhy+2*delz,dtau)
      call memo(2,';')
      write(ifmt,'(3(a,i3))')
     .'fluid size ',nxihy+4,' x ',nyihy+4,' x ',nzihy+4
      call clop(3)

      write(ifmt,'(a)')'set fluid initial conditions'

      call clop(3)
      if(ncnthlle.eq.1.and.itest.eq.1)then !*** test ******
          write(ifmt,'(1x,3a,$)')
     .   'TEST: read ico from ',fnio(1:nfnio),' ...'
          call initichlle(fnio(1:nfnio)//CHAR(0))
          write(ifmt,'(a)')'  done'
      endif

c fill IcoEArray array parameter from IcoE values
      do ix = 1, nxicomax
         do iy = 1, nyicomax
            do iz = 1, nzicomax
               IcoEArray(ix, iy, iz)=getIcoIcoE(ix, iy, iz)
            end do
         end do
      end do
c call iniflowini : IcoEArray is a fortran array
      call IniFlowIni(nxicomax,nyicomax,nzicomax,IcoEArray)
c fill IcoE from IcoEArray array parameter values
      do ix = 1, nxicomax
         do iy = 1, nyicomax
            do iz = 1, nzicomax
               call setIcoIcoE(ix, iy, iz, IcoEArray(ix, iy, iz))
            end do
         end do
      end do


      eeitau=0
      eejtau=0
      delxico=(xmaxico-xminico)/nxico
      delyico=(ymaxico-yminico)/nyico
      delzico=(zmaxico-zminico)/nzico
      do ix = 1, nxihy
      do iy = 1, nyihy
      do iz = 1, nzihy
        e=0
        vx=0
        vy=0
        vz=0
        fu=0
        fd=0
        fs=0
        x   = xmnhy  +(ix-1)  *(xmxhy-xmnhy)  /(nxihy-1)
        y   = ymnhy  +(iy-1)  *(ymxhy-ymnhy)  /(nyihy-1)
        z   = zmnhy  +(iz-1)  *(zmxhy-zmnhy)  /(nzihy-1)
        jx=1+(x-xminico-delxico/2)/delxico
        jy=1+(y-yminico-delyico/2)/delyico
        jz=1+(z-zminico-delzico/2)/delzico
        if(    jx.gt.0.and.jx.lt.nxico
     .    .and.jy.gt.0.and.jy.lt.nyico
     .    .and.jz.gt.0.and.jz.lt.nzico)then
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
          e=e  +wx(i)*wy(j)*wz(k)*getIcoIcoE(  jx+i-1,jy+j-1,jz+k-1)
          vx=vx+wx(i)*wy(j)*wz(k)*getIcoIcoV(1,jx+i-1,jy+j-1,jz+k-1)
          vy=vy+wx(i)*wy(j)*wz(k)*getIcoIcoV(2,jx+i-1,jy+j-1,jz+k-1)
          vz=vz+wx(i)*wy(j)*wz(k)*
     .         getIcoIcoV(3,jx+i-1,jy+j-1,jz+k-1)*tau0
          fu=fu+wx(i)*wy(j)*wz(k)*getIcoIcoF(1,jx+i-1,jy+j-1,jz+k-1)
          fd=fd+wx(i)*wy(j)*wz(k)*getIcoIcoF(2,jx+i-1,jy+j-1,jz+k-1)
          fs=fs+wx(i)*wy(j)*wz(k)*getIcoIcoF(3,jx+i-1,jy+j-1,jz+k-1)
          enddo
          enddo
          enddo
          e=max(e,0d0)
          eps=e
          zz=z
          !~~~~~~~~~~~~~~~~
          gamma=1
          ggi=1-vx**2-vy**2-vz**2
          if(ggi.gt.0.)gamma=1./sqrt(ggi)
          u0=gamma
          u3=gamma*vz
          pi00=0.
          pi03=0.
          pi33=0.
          eeitau=eeitau+eflow(eps,eps/3,zz,u0,u3,pi00,pi03,pi33)
     .                    *tauzer*delx*dely*delz
          !~~~~~~~~~~~~~~~~~~
          call addIniFlow(x,y,vx,vy)
          !~~~~~~~~~~~~~~~~
          gamma=1
          ggi=1-vx**2-vy**2-vz**2
          if(ggi.gt.0.)gamma=1./sqrt(ggi)
          u0=gamma
          u3=gamma*vz
          eejtau=eejtau+eflow(eps,eps/3,zz,u0,u3,pi00,pi03,pi33)
     .                    *tauzer*delx*dely*delz
          !~~~~~~~~~~~~~~~~~~
        endif
        nb=(fu+fd+fs)/3.
        nq=(2*fu-fd-fs)/3.
        ns=-fs

        if(e.ge.20*oEeos)then
        ncnthl=ncnthl+1
        call setAccumJerr(14,getAccumJerr(14)+1)
        if(50.gt.ncnthl)then
         write(ifmt,*)'hlle: energy density very big: ',e
         !write(ifmt,*)'         e(',ix,iy,iz,')=',e
         !write(ifmt,*)(((IcoE(jx+i-1,jy+j-1,jz+k-1),k=1,2),j=1,2),i=1,2)
        endif
        endif
        if((nb.gt.1e5.or.nq.gt.1e5.or.ns.gt.1e5)
     .   .and.50.gt.ncnthlf)then
         ncnthlf=ncnthlf+1
         write(ifmt,*)'hlle: nb/nq/ns density very big: ',nb,nq,ns
         !write(ifmt,*)'      nb/nq/ns(',ix,iy,iz,')=',nb,nq,ns
         !do m=1,3
         !write(ifmt,*)'         '
         !.   ,(((IcoF(m,jx+i-1,jy+j-1,jz+k-1),k=1,2),j=1,2),i=1,2)
         !enddo
        endif

        if(ncnthlle.eq.1.and.itest.eq.1)then !*** test ******
            xcc = getxhlle(ix)
            ycc = getyhlle(iy)
            zcc = getzhlle(iz)
            if(ioeos.eq.1)then
            call icgethlle(  xcc, ycc, zcc, ecc, nbcc, nqcc, nscc,
     .       vxcc, vycc, vzcc)
            elseif(ioeos/10.eq.2)then
            call icgethlle3f(xcc, ycc, zcc, ecc, nbcc, nqcc, nscc,
     .       vxcc, vycc, vzcc)
            else
            stop'\n\n STOP in  hlle, wrong choice of ioeos\n\n'
            endif
            if(jx.ne.nxico.and.jy.ne.nyico.and.jz.ne.nzico
     .      .and.jx.ne.0.and.jy.ne.0.and.jz.ne.0)then
            if(abs(xcc-x).gt.1d-6)stop'1\n\n'
            if(abs(ycc-y).gt.1d-6)stop'2\n\n'
            if(abs(zcc-z).gt.1d-6)stop'3\n\n'
            if(abs(ecc-e).gt.0.01)stop'4\n\n'
            if(ioeos/10.eq.2)then
            if(abs(nbcc-nb).gt.0.001)stop'5\n\n'
            if(abs(nqcc-nq).gt.0.001)stop'6\n\n'
            if(abs(nscc-ns).gt.0.001)stop'7\n\n'
            endif
            if(abs(vxcc-vx).gt.0.001)stop'8\n\n'
            if(abs(vycc-vy).gt.0.001)stop'9\n\n'
            if(abs(vzcc-vz).gt.0.001)stop'10\n\n'
            endif
        endif

        e=max(0d0,e)
        v2=vx**2+vy**2+vz**2
        if(v2.ge.1d0)then
         v2=v2*1.001
         vx=vx/sqrt(v2)
         vy=vy/sqrt(v2)
         vz=vz/sqrt(v2)
        endif
        if(iozerof.eq.1)then
          nb=0
          nq=0
          ns=0
        endif
        if(e.lt.1e-10)then
          vx=0
          vy=0
          vz=0
        endif
        call checkNaN(e, nb, nq, ns, vx, vy, vz)
        call icsethlle(ix+2,iy+2,iz+2, tau0, e, nb, nq, ns, vx, vy, vz)
c        write(ifch,'(a,3i5,3f8.3)')'hlle ini  x y z b123'
c     .     ,ix,iy,iz,nb, nq, ns
      enddo
      enddo
      enddo

      !ghosts:
      zzz=0d0

      do ix = 1, nxihy+4
      do iy = 1, nyihy+4
      do iz = 1,2
        call icsethlle(ix, iy, iz, tau0, zzz,zzz,zzz,zzz,zzz,zzz,zzz)
      enddo
      do iz = nzihy+3,nzihy+4
        call icsethlle(ix, iy, iz, tau0, zzz,zzz,zzz,zzz,zzz,zzz,zzz)
      enddo
      enddo
      enddo

      do ix = 1, nxihy+4
      do iz = 1, nzihy+4
      do iy = 1,2
        call icsethlle(ix, iy, iz, tau0, zzz,zzz,zzz,zzz,zzz,zzz,zzz)
      enddo
      do iy = nyihy+3,nyihy+4
        call icsethlle(ix, iy, iz, tau0, zzz,zzz,zzz,zzz,zzz,zzz,zzz)
      enddo
      enddo
      enddo

      do iy = 1, nyihy+4
      do iz = 1, nzihy+4
      do ix = 1,2
        call icsethlle(ix, iy, iz, tau0, zzz,zzz,zzz,zzz,zzz,zzz,zzz)
      enddo
      do ix = nxihy+3,nxihy+4
        call icsethlle(ix, iy, iz, tau0, zzz,zzz,zzz,zzz,zzz,zzz,zzz)
      enddo
      enddo
      enddo

      write(ifmt,'(a)')'initialize hydro'
      call clop(3)
      call inithydrohlle(tau0, tau_max, dtau)

      eektau=eejtau
      eejtau=eechk(tauzer,nxihy,nyihy,nzihy
     .    ,xmnhy,xmxhy,ymnhy,ymxhy,zmnhy,zmxhy)

      !call checkengy('before hyd loop ') !counts only ist=0, gives useless warning
      if(icotabr.eq.1)call setIcoEe1ico(eejtau)
      egain=eejtau-getIcoEe1ico()
      write(ifmt,'(a,4f8.0,3x,a,f8.0)')
     . '  ECHK4',getIcoEe1ico(),eeitau,eektau,eejtau,'  gain',egain
         !--------------------------------------------------------------
         !   ECHK4
         !--------------------------------------------------------------
         !ee1ico = E from Ico
         !eeitau = E from initial fluid table
         !eektau = E from initial fluid table + IniFlow
         !eejtau = E from initial fluid table + IniFlow + visc (eechk)
         !--------------------------------------------------------------
      if(icotabr.eq.0)then
        eistxx=getIcoEistico()
        esollxx=eistxx-egain
        call setIcoIchkengy(1)
        write(ifmt,'(a,2f8.0)')
     . '  eist esoll ',eistxx,esollxx
        call setIcoEistxx(eistxx)
        call setIcoEsollxx(esollxx)
        call chkengy(ierrchk)
      endif
      call setIcoEe1hll(eejtau)

      if(ncnthlle.eq.1.and.itest.eq.1)
     . write(ifmt,'(1x,a)')'TEST: passed succesfully.'

      ntau=1
      istep=1
      ntsteps=0
      mmsteps=4*ifathlle        !must identical to def in hyt.f
      del_tau=dtauhy/ifathlle   !must identical to def in hyt.f
      tau=tauzer
      call setHydynTauhy(ntau,tau)
      write(ifmt,'(i3,3x,a,f6.2,$)')ntau,'  tau',tau
      call TransferHlle(ntau)
      call checkmemory(memory)
      call envarget(exxx, efull,etotsurf,nbsurf)
      eetauxx=eejtau !from eechk
      ratioee=1
      call initializeHydynRatioeex(real(1.0))
      zero=0
      emx=getHydynEmx() 
c      xmx=getHydynXmx()
c      ymx=getHydynYmx()
c      zmx=getHydynZmx()
      write(ifmt,'(a,f8.2,2x,a,2f8.0,f6.3,3f8.0)')
     . '  emx',emx
     . ,'ECHK5',getHydynEetau(),eejtau,getHydynRatioeex(ntau),zero
                   !.      ,efull,etotsurf+exxx
                   !.      ,eeinn,xmx,ymx,zmx
      call timer(iutime)
      tiu3=iutime(3)
      tiu4=iutime(4)
      tidisu=0

      write(ifmt,'(a)')'hydro loop'
      call clop(3)

      do i = 0, ntaumx-1
        if(ntau+1.ge.ntauhxx)then
          write(ifmt,*)'ERROR ntauhxx too small'
          stop 'ERROR ntauhxx too small'
        endif
        istepx=istep
        if(tau.gt.tau1.and.ntsteps.eq.0)istep=2
        if(tau.gt.tau2.and.ntsteps.eq.0)istep=4
        if(tau.gt.tau3.and.ntsteps.eq.0)istep=8
        dtau = istep * del_tau
        if(istep.gt.1.and.istep.ne.istepx)then
          call memo(1,'reduce fluid;')
          call reducegridhlle(dtau,istep) !coarse-grain istep times
          call memo(2,';')                !compared to initial grid
          write(ifmt,'(3(a,i3))')'new fluid size ',(nxihy+3)/istep+1
     .    ,' x ',(nyihy+3)/istep+1,' x ',nzihy+4
        endif
        call dtauhlle(dtau)           
        call makestephlle(i)
        ntsteps=ntsteps+istep
        tau=tau+dtau
        if(ntsteps.gt.mmsteps)stop'####### ERROR 02112015 #######'
        if(ntsteps.eq.mmsteps)then !any change must be  consistent with
          ntau=ntau+1              !delta_tau_out def in hyt.f
          call setHydynTauhy(ntau,tau)
          if(abs(getHydynTauhy(ntau)-getHydynTauhy(ntau-1)-4*dtauhy)
     .         .gt.0.001)
     .    stop'####### ERROR 27122011 #######'
          call timer(iutime)
          tidi=iutime(3)-tiu3 + (iutime(4)-tiu4)*1d-3
          tidisu=tidisu+tidi
          tiu3=iutime(3)
          tiu4=iutime(4)
          write(ifmt,'(2i3,a,f6.2,$)')ntau,mmsteps/istep,'  tau',tau
          eetauaa=eechk(tau,nxihy,nyihy,nzihy
     .    ,xmnhy,xmxhy,ymnhy,ymxhy,zmnhy,zmxhy)
          ratioee=max(ratioee,eetauaa/eetauxx)
          call setHydynRatioeex(ntau,ratioee)
          call TransferHlle(ntau)
          call checkmemory(memory)
          call envarget(exxx, efull,etotsurf,nbsurf)
          !-------------------------------------------------------------
          ! A surface is defined by energy density e = eCrit, it is a 3D 
          ! surface (in 3+1D spacetime), and the precise value is set at:
          ! YK/fld.cpp:17:  const double eCrit = 0.5 ;
          ! Arguments of envarget: 
          ! exxx    : energy inside the surface (hydro cells with e>eCrit)
          ! efull   : full energy (all hydro cells)
          ! etotsurf: Eflow through the surface
          ! nbsurf  : flow of baryon charge through the surface, 
          !              \int n_B * d\sigma_mu u^\mu
          !--------------------------------------------------------------        
          emx=getHydynEmx()
          eetau2=getHydynEetau2()
          write(ifmt,'(a,f8.2,2x,a,f8.0,$)')
     .   '  emx',emx,'ECHK5',eetau2
c         xmx=getHydynXmx()
c         ymx=getHydynYmx()
c         zmx=getHydynZmx()
          write(ifmt,'(f8.0,f6.3,3f8.0)')eetauaa,getHydynRatioeex(ntau),
     .         tidi 
                   !.    ,efull,etotsurf+exxx
                   !.    ,eeinn,xmx,ymx,zmx
         !--------------------------------------------------------------
         !   ECHK5 (E means energy including visc)
         !--------------------------------------------------------------
         ! eetau2  = E from TransferHlle  (sparce grid)
         ! eetauaa = E from eechk    (full grid)
         ! ratioee = eetauaa / eetauxx (initial E from eechk)    
         ! efull   = E all hydro cells from envarget
         ! etotsurf+exxx
         !    etotsurf  = E flow through cylinder from envarget
         !    exxx      = E inside cylinder  from envarget 
         !--------------------------------------------------------------
         ! eeinn     = E from eechk  (full grid, inner range)
         !--------------------------------------------------------------
          call clop(3)
          if(emx.lt.efin)goto 2
          if(tau.gt.tau_max)goto 2
          ntsteps=0
        endif
      end do
  2   continue

      ntauhy=ntau
      if(iopcnt.gt.0)then
        j= 1+mod(iopcnt-1,20)
        call setCcotiCoti(1,j,getCcotiCoti(1,j)+1.0)
        call setCcotiCoti(2,j,getCcotiCoti(2,j)+tidisu)
      endif
      cputime=gettimehlle()
      call memo(1,'destroy fluid;')
      call destroyhlle()
      call memo(2,';')

      call setHydynNtauhy(ntauhy)

      call destroyResc()
      call destroyEvent()
      call destroyEosp()
      call destroyCore()

      !write(ifmt,*)'finish hlle;   tau =',tau,'   CPU time ='
      !. ,nint(cputime),' sec'
      call checktime('exit hlle;') 
      end subroutine hlle

c-----------------------------------------------------------------------
      subroutine TransferHlle(ntau)
c-----------------------------------------------------------------------
c      use ArrayModule, only: set
c      use hocoModule, only: epsc, velc, barc, sigc
c      use hocoModule, only: velc, barc, sigc
      
#include "aaa.h"
#include "ho.h"
      double precision e, p, nb, nq, ns, vx, vy, vz, viscflag
      double precision temc
      double precision pishear(10), pibulk
      real emx

      real u(4),tmunu(4,4)
      common/citsy/itsy(4,4)
      common/ctemc/temc(netahxx,ntauhxx,nxhxx,nyhxx)
      common/cpssc/pssc(netahxx,ntauhxx,nxhxx,nyhxx)
      common/cvifl/vifl(netahxx,ntauhxx,nxhxx,nyhxx)
c      common/cmaxi/emx,xmx,ymx,zmx
      common/ceetau/eetau,eetau2
c      common/cmm/mmx,mmy,mmz
      data ncnttrhl/0/
      save ncnttrhl

      integer getHydynMmx
      integer getHydynMmy
      integer getHydynMmz
      integer getAccumJerr

      emx=0
      call setHydynEetau(0)
      d1= (xmaxhy-xminhy) /(nxhy-1)
      d2= (ymaxhy-yminhy) /(nyhy-1)
      d3= (zmaxhy-zminhy) /(nzhy-1)
      do nz=1,nzhy
        do nx=1,nxhy
         do ny=1,nyhy
           nxhlle = getHydynMmx()+(nx-1)*ifathlle
           nyhlle = getHydynMmy()+(ny-1)*ifathlle
           nzhlle = getHydynMmz()+(nz-1)*ifazhlle
           !-------------------------------------------
           !  get values from hlle :
           !-------------------------------------------
           call getvalueshlle(nxhlle+2, nyhlle+2, nzhlle+2,
     .      e, p, nb, nq, ns, vx, vy, vz, viscflag)
           call getvisc10hlle(nxhlle+2,nyhlle+2,nzhlle+2,pishear,pibulk)
           !-------------------------------------------
           ! pishear is a 10-dim vector, the tensor is :
           !  pishear ( itsy(i,j) ), with i,j=1,4  (4 = tau)
           !-------------------------------------------
           ! shear stress tensor pishear from hlle is the
           !  tilde version of the tensor stricly in Milne coordinates:
           ! pishear = \tilde(\pi}^\mu,nu with \pi^\mu,nu in Milne
           ! vz = \tilde{u}^eta / \tilde{u}^tau
           !------------------------------------------- 
           ! Be "Bjorken frame" a frame moving with  y=eta with respect
           !   to the lab frame.
           ! It can be shown: 
           !   1) vz = velocity in Minkowski space (ordinary velocity), 
           !              in Bjorken frame 
           !   2) the tilde shear stress tensor in Milne is identical
           !     to the tensor in Minkowski space, in Bjorken frame
           !-------------------------------------------
           ! Be P^\mu = T^mu,nu * dSigma_nu     mu,nu = 1,4  (4=0)
           ! One can show:  
           !   P^\mu(Minkowski,lab frame) 
           !      = tilde K^mu_nu  tilde{T}^nu,rho * \tilde{dSigma}_rho
           !                 (  1                              )   <-- 1
           !  {K^mu_nu} =   |        1                          |  <-- 2
           !                |            cosh(eta)   sinh(eta)  |  <-- 3
           !                 (           sinh(eta)   cosh(eta) )   <-- 4 
           !-------------------------------------------
           if(e.ge.20*oEeos)then
             ncnttrhl=ncnttrhl+1
             call setAccumJerr(15,getAccumJerr(15)+1)
             if(50.gt.ncnttrhl)then
               write(ifmt,*)
     .         'TransferHlle: energy density very big: ',e
               !write(ifmt,*)'        e(',nz,ntau,nx,ny,')=',e
             endif
           endif
           if(e.gt.dble(emx))then
             emx=e
             call setHydynEmx(emx)
             i=nx
             j=ny
             k=nz
           endif
           call setHydynEpsc(nz,ntau,nx,ny,e)
           call setHydynSigc(nz,ntau,nx,ny,dble(0))
           call setHydynBarc(1,nz,ntau,nx,ny,nb)
           call setHydynBarc(2,nz,ntau,nx,ny,nq)
           call setHydynBarc(3,nz,ntau,nx,ny,ns)
           call setHydynVelc(1,nz,ntau,nx,ny,vx)
           call setHydynVelc(2,nz,ntau,nx,ny,vy)
           call setHydynVelc(3,nz,ntau,nx,ny,vz)
           vifl(nz,ntau,nx,ny)=viscflag
           !~~~~~~ compute T^mu^nu ~~~~~~~
           gamma=1
           ggi=1-vx**2-vy**2-vz**2
           if(ggi.gt.0.)gamma=1./sqrt(ggi)
           eps=e
           b1=nb
           b2=nq
           b3=ns
           p=PiEos(4,eps,b1,b2,b3)
           press=p+pibulk
           !~~~ \tilde(u}^i :
           u(1)=gamma*vx
           u(2)=gamma*vy
           u(3)=gamma*vz
           u(4)=gamma
           !~~~~~ \tilde{T}^ij
           do i=1,4
           do j=1,4
            tmunu(i,j) = (eps+press)*u(i)*u(j) + pishear(  itsy(i,j)  )
           enddo
           enddo
           tmunu(4,4)=tmunu(4,4) - press
           do i=1,3
             tmunu(i,i) = tmunu(i,i) + press
           enddo
           !~~~~~~~ fill emuc and some arrays ~~~~~
           do i=1,4
             do j=1,i
               call emucset( itsy(i,j) ,nz,ntau,nx,ny, tmunu(i,j) )
             enddo
           enddo
           temc(nz,ntau,nx,ny)=PiEos(5,eps,b1,b2,b3)  !*** needed in gethyval
           pssc(nz,ntau,nx,ny)=PiEos(4,eps,b1,b2,b3)  !*** needed in gethyval
           tau=getHydynTauhy(ntau)
           !~~~~~~~~~~~~~~~~~~~~~~
         enddo
        enddo
      enddo
      call setHydynXmx(xminhy  +(i-1)  *d1)
      call setHydynYmx(yminhy  +(j-1)  *d2)
      call setHydynZmx(zminhy  +(k-1)  *d3)
      !~~~~ test ~~~~~
      eetau2=0
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
           !    = tmunu(i,4) * tau * d1 * d2 * d3  ,i=1,4 
           ! Then trafo to Minkowski, lab, using K tensor, see above
           eetau2=eetau2
     .     + ( tmunu(3,4)*sinh(eta)+tmunu(4,4)*cosh(eta) )*tau*d1*d2*d3
         enddo
        enddo
      enddo

      call setHydynEetau2(eetau2)
      !..........
      !print*,'ntau, nxhy, nyhy',ntau, nxhy, nyhy
      !do i=1,16
      ! write(*,'(16i5)')(nint(100*epsc(12,ntau,i,j)),j=16,31)
      !enddo
      !if(ntau.eq.3)stop
      !..........
      end subroutine TransferHlle

c-----------------------------------------------------------------------
      real function eechk(tau,nxihy,nyihy,nzihy
     . ,xmnhy,xmxhy,ymnhy,ymxhy,zmnhy,zmxhy)
c-----------------------------------------------------------------------
#include "aaa.h"
#include "ho.h"
      double precision xmnhy,xmxhy,ymnhy,ymxhy
     . ,zmnhy,zmxhy, viscflag
      double precision pishear(10), pibulk
      double precision e, p, nb, nq, ns, vx, vy, vz
      real u(4),tmunu(4,4)
c      integer mmx,mmy,mmz
      common/citsy/itsy(4,4)
      common/ceflowvi/pivi,eflowvi
c      common/ceevtau/eevtau
c      common/ceeinn/eeinn
c      common/cmm/mmx,mmy,mmz

      eetau3=0
c      eevtau=0
c      eeinn=0
      d1= (xmxhy-xmnhy) /(nxihy-1)
      d2= (ymxhy-ymnhy) /(nyihy-1)
      d3= (zmxhy-zmnhy) /(nzihy-1)
      do nz=1,nzihy
        do nx=1,nxihy
         do ny=1,nyihy
           call getvalueshlle(nx+2, ny+2, nz+2,
     .      e, p, nb, nq, ns, vx, vy, vz, viscflag)
           call getvisc10hlle(nx+2, ny+2, nz+2, pishear, pibulk)
           !~~~~~~ compute T^mu^nu ~~~~~~~
           gamma=1
           ggi=1-vx**2-vy**2-vz**2
           if(ggi.gt.0.)gamma=1./sqrt(ggi)
           eps=e
           b1=nb
           b2=nq
           b3=ns
           p=PiEos(4,eps,b1,b2,b3)
           u(1)=gamma*vx
           u(2)=gamma*vy
           u(3)=gamma*vz
           u(4)=gamma
           press=p+pibulk
           do i=1,4
           do j=1,4
            tmunu(i,j) = (eps+press)*u(i)*u(j) + pishear(  itsy(i,j)  )
           enddo
           enddo
           tmunu(4,4)=tmunu(4,4) - press
           do i=1,3
             tmunu(i,i) = tmunu(i,i) + press
           enddo
           !~~~~ test ~~~~~
           z=zminhy  +(nz-1)  *d3
           eta=z
           !T^mu^nu*dSigma_nu=T^mu^0*dSigma_0 in Bjorken frame
           !  then boost of the resulting vector into cms :
           delee=
     .       ( tmunu(4,3)*sinh(eta)+tmunu(4,4)*cosh(eta) )*tau*d1*d2*d3
           eetau3=eetau3 + delee
c           if(  nx.ge.mmx.and.nx.le.nxihy+1-mmx
c     .     .and.ny.ge.mmy.and.ny.le.nyihy+1-mmy
c     .     .and.nz.ge.mmz.and.nz.le.nzihy+1-mmz  )then
c             eeinn=eeinn+delee
c           endif
           !~~~endtest~~~~
         enddo
        enddo
      enddo
      eechk=eetau3
      end

      subroutine FillZero(ntau)
c      use ArrayModule, only: set
c      use hocoModule, only: epsc, sigc, velc, barc
c      use hocoModule, only: sigc, velc, barc
      
#include "ho.h"

      do nz=1,nzhy
        do nx=1,nxhy
         do ny=1,nyhy
           call setHydynEpsc(nz,ntau,nx,ny,dble(0))
           call setHydynSigc(nz,ntau,nx,ny,dble(0))
           call setHydynBarc(1,nz,ntau,nx,ny,dble(0))
           call setHydynBarc(2,nz,ntau,nx,ny,dble(0))
           call setHydynBarc(3,nz,ntau,nx,ny,dble(0))
           call setHydynVelc(1,nz,ntau,nx,ny,dble(0))
           call setHydynVelc(2,nz,ntau,nx,ny,dble(0))
           call setHydynVelc(3,nz,ntau,nx,ny,dble(0))
         enddo
        enddo
      enddo
      end 

c-----------------------------------------------------------------------
        subroutine checkNaN(e, nb, nq, ns, vx, vy, vz)
c-----------------------------------------------------------------------
      double precision e, nb, nq, ns, vx, vy, vz, a(7),b
      a(1)=e
      a(2)=nb
      a(3)=nq
      a(4)=ns
      a(5)=vx
      a(6)=vy
      a(7)=vz
      do i=1,7
      b=a(i)
      if(.not.(b.le.0..or.b.ge.0.))then !NaN catch
        print*,'ERROR in checkNaN'
        print*,'   e, nb, nq, ns, vx, vy, vz :'
     .  ,e, nb, nq, ns, vx, vy, vz
        stop'\n\n ERROR in checkNaN \n\n'
      endif
      enddo
      end

c-----------------------------------------------------------------------
      subroutine getEosHlle
c-----------------------------------------------------------------------

      !------------------------------------------------
      !eos(1)= energy density (GeV/fm^3!)
      !eos(2)= entropy density (1/fm^3!)
      !eos(3)= baryon density (1/fm^3!)
      !eos(4)= pressure (GeV/fm^3!)"'
      !eos(5)= temperature (GeV)"'
      !eos(6)= baryon chemical potential (GeV)
      !eos(7)=  strangeness chemicalpotential (GeV)
      !-----------------------------------------------

#include "aaa.h"
#include "ho.h"
#if __ROOT__
      !external fhir
      double precision emax,e0,nmax,n0
      if(ioeos.eq.22)then
        call eosparams(1)
        call eostabs(1,1,1)
      elseif(ioeos.ne.0)then
        call eostabs(0,0,1)
      else
        stop'ERROR in getEosHlle, ioeos is zero'
      endif
      call eostest
      call eosrangeshlle(emax,e0,nmax,n0,ne,nn)
      uEeos=sngl(e0)
      oEeos=sngl(emax)
      mxEeos=ne
      mxBeos=nn
      write(ifmt,*)'uEeos=',uEeos,'  oEeos=',oEeos
      write(ifmt,*)'mxEeos=',mxEeos,'  mxBeos=',mxBeos
      call clop(3)
#else
        stop "ROOT needed for getEosHlle"
#endif
      end

      function cothi(id)  !chem/T of PCE a la Hirano
#include "aaa.h"
#include "ho.h"
      parameter (mxeoshi=50001 )
      common/ceoshi/eoshi(17,mxeoshi),meoshi
      real cothiar(17)
      character line*400
      data ncnthi /0/
      save ncnthi
      save cothiar
      ncnthi=ncnthi+1
      if(ncnthi.eq.1)then
        open(98,file=
     .  fnnx(1:nfnnx)//'../eos3f/Tbo/Hirano_mu.dat',status='old')
        meoshi=0
        do
         read(98,'(a)',end=778)line
         meoshi=meoshi+1
         read(line,*)(eoshi(i,meoshi),i=1,17)
         do i=1,17
         eoshi(i,meoshi)=eoshi(i,meoshi)*0.001
         enddo
         !print*,meoshi,eoshi(1,meoshi),eoshi(2,meoshi)
        enddo
  778   continue
        close(98)
        k=1
        do while(eoshi(1,k).gt.tfrout)
          !print*,'k,eoshi,tfrout',k,eoshi(1,k),tfrout
          k=k+1
        enddo
        f=(tfrout-eoshi(1,k))/(eoshi(1,k-1)-eoshi(1,k))
        do ihi=1,17
          chpot=eoshi(ihi,k)*(1-f) + eoshi(ihi,k-1)*f
          cothiar(ihi)=chpot/tfrout
          !print*,'tfrout,ihi,chpot,cothiar',tfrout,ihi,chpot,cothiar(ihi)
        enddo
      endif
      cothi=0
      if(meoshi.gt.0)then
       ihi=idxxHiranoTable(id)
       if(ihi.gt.0)cothi=cothiar(ihi)
      endif
      end
      integer function idxxHiranoTable(id)
      ida=abs(id)
      ihi=0
      if(ida.eq.120.or.id.eq.110)   ihi=2       !pion
      if(id.eq.220)                 ihi=3       !eta
      if(ida.eq.121.or.id.eq.111)   ihi=4       !rho
      if(id.eq.221)                 ihi=5       !omega
      !if(id.eq.)                        ihi=6         !sigma
      if(id.eq.330)                 ihi=7       !eta prime
      !if(id.eq.)                        ihi=8         !f_0
      !if(id.eq.)                        ihi=9         !a_0
      if(id.eq.331)                 ihi=10      !phi
      !if(id.eq.)                        ihi=11          !h_1
      if(ida.eq.130.or.ida.eq.230)  ihi=12      !K
      if(ida.eq.131.or.ida.eq.231)  ihi=13      !K star
      if(ida.eq.1120.or.ida.eq.1220)ihi=14      !N
      if(ida.eq.1111.or.ida.eq.1121)ihi=15      !Delta
      if(ida.eq.1221.or.ida.eq.2221)ihi=15      !Delta
      if(ida.eq.2130)               ihi=16      !Lambda
      if(ida.eq.1130.or.ida.eq.1230)ihi=17      !Sigma
      if(ida.eq.2230)               ihi=17      !Sigma
      idxxHiranoTable=ihi
      end

c-----------------------------------------------------------------------
      subroutine fhir(e,nb,nq,ns,p,t,mub,muq,mus)
c-----------------------------------------------------------------------
      double precision p,T,mub,muq,mus,  e, nb, nq, ns
      es=e
      b1=nb
      b2=nq
      b3=ns
      if(es.le.3.)then
        p   =  PiEos(4,es,b1,b2,b3)    !pressure
        T   =  PiEos(5,es,b1,b2,b3)    !temperature
      else
        call fhir3(es,ps,ts)
        p=ps
        T=ts
      endif
      mub =  0  !baryon chemical potential
      muq =  0  !charge chemical potential
      mus =  0  !strangeness chemical potential
      end
      subroutine fhir3(e,p,t)
      !..............
      ! Data exist up to e = 3. Above this analytic formula such as
      ! T = (30*(e-B)/pi^2/deg)^0.25
      ! P = (e-4*B)/3
      ! s = (e+P)/T
      ! where B^{1/4} = 0.24719 GeV  and deg=47.5
      !..............
#include "aaa.h"
      bag=0.24719**4
      deg=47.5
      fcr=0.197327**3
      t=(30*(e*fcr-bag)/pi**2/deg)**0.25
      p=(e-4*bag/fcr)/3.
      end


#endif


c-------------------------------------------------------------------------------
      subroutine IniFlowIni(nxicomx,nyicomx,nzicomx,IcoEx)
c-------------------------------------------------------------------------------
#include "aaa.h"
      common/cninx/ninx
      common/cecc/ecct,rrrt
      common/ciniflo/ecccoe(5),eccphi(5) /ciniflo2/rrr,vvv,xav0,yav0
      real avecos(5),avesin(5)
      integer, intent(in) :: nxicomx
      integer, intent(in) :: nyicomx
      integer, intent(in) :: nzicomx
      real, intent(in) :: IcoEx
      dimension IcoEx(nxicomx,nyicomx,nzicomx)
      integer nxico,nyico,nzico
      real xminico,xmaxico,yminico,ymaxico,zminico,zmaxico
      
      call getIcoDim(nxico,nyico,nzico)
      call getIcoBound(
     .     xminico, xmaxico,
     .     yminico, ymaxico,
     .     zminico, zmaxico
     . )


      pi=3.1415927
      do m=2,5
      eccphi(m)=0
      ecccoe(m)=0
      enddo
      rrr=0
      vvv=0
      xav0=0
      yav0=0
      x1=0
      y1=0
      aa=0
      do ix=1,nxico
        xi=xminico+(ix-0.5)*(xmaxico-xminico)/nxico
        do iy=1,nyico
          yi=yminico+(iy-0.5)*(ymaxico-yminico)/nyico
          do iz=1,nzico
            eps=IcoEx(ix,iy,iz)
            aa=aa+eps
            x1=x1+eps*xi
            y1=y1+eps*yi
          enddo
        enddo
      enddo
      if(aa.ne.0.)then
       xav0=x1/aa
       yav0=y1/aa
       rr=0
       rrxx=0
       do m=2,5
         avecos(m)=0
         avesin(m)=0
       enddo
       do ix=1,nxico
         xi=xminico+(ix-0.5)*(xmaxico-xminico)/nxico   - xav0
         do iy=1,nyico
           yi=yminico+(iy-0.5)*(ymaxico-yminico)/nyico - yav0
           rixri=xi**2+yi**2
           rixrixx=xi**2+yi**2
           !rixrixx=1.0  !   ******* modified weight, to avoid pb with isolated peaks
           phi=polar(xi,yi)
           do iz=1,nzico
             eps=IcoEx(ix,iy,iz)
             rr=rr+eps*rixri
             rrxx=rrxx+eps*rixrixx
             do m=2,5
               avecos(m)=avecos(m)+eps*rixrixx*cos(m*phi)
               avesin(m)=avesin(m)+eps*rixrixx*sin(m*phi)
             enddo
           enddo
         enddo
       enddo
       do m=2,5
       avecos(m)=avecos(m)/rrxx
       avesin(m)=avesin(m)/rrxx
       eccphi(m)=polar(avecos(m),avesin(m)) / m  - pi/m
       ecccoe(m)=sqrt(avecos(m)**2+avesin(m)**2)
       enddo
       rrr= sqrt(rr/aa)
       vvv=floini
       if(ninx.eq.1)
     .  write(ifmt,'(a,4f6.2,3x,4f7.4)')'ecc='
     .  ,(eccphi(m),m=2,5),(ecccoe(m),m=2,5)
       e2=ecccoe(2)
       e3=ecccoe(3)
      else
       e2=0
       e3=0
      endif
      if(ioecc.eq.1)then
        if(.not.(e2.gt.valecc.and.e3.gt.valecc))itrigg=0
      elseif(ioecc.eq.2)then
        if(.not.(e2.lt.valecc.and.e3.gt.valecc))itrigg=0
      elseif(ioecc.eq.3)then
        if(.not.(e2.lt.valecc.and.e3.lt.valecc))itrigg=0
      elseif(ioecc.eq.4)then
        if(.not.(e2.gt.valecc.and.e3.lt.valecc))itrigg=0
      endif
      if(ninx.eq.1)write(ifmt,'(a,i1,i5,2f7.4)')
     . 'itrigg=',itrigg,ioecc,e2,e3
      end



      subroutine addIniFlow(xo,yo,vx,vy)
      double precision xo,yo,vx,vy,b1,b2,v1,v2,gm,r,vvm,rnll,vv,v,a,bv
      common/ciniflo/ecccoe(5),eccphi(5) /ciniflo2/rrr,vvv,xav0,yav0
      save kkkk
      kkkk=kkkk+1
      vxbef=vx
      x=xo-xav0 !x-position relative to average
      y=yo-yav0 !y-position relative to average
      phi=polar(x,y)
      r=sqrt(x**2+y**2)
      do m=2,3
        rnll=rrr*sqrt(1-ecccoe(m)*cos(m*(phi-eccphi(m))))
        r0=rnll/rrr
        r0p=0.5/r0*m*ecccoe(m)*sin(m*(phi-eccphi(m)))
        vvm=vvv
        if(m.eq.2.and.ecccoe(2).lt.ecccoe(3))vvm=0
        if(m.eq.3.and.ecccoe(3).le.ecccoe(2))vvm=0
        vv=min(0.99999999999999d0,vvm*r/rnll)
        v=vv*sqrt(1+0.25*ecccoe(m)*cos(m*(phi-eccphi(m))))
        !vector perpendicular to r0(phi)*e_r: r0*e_r - r0p*e_phi
        a=sqrt(r0**2+r0p**2)
        b1=v/a*(r0*cos(phi)+r0p*sin(phi))
        b2=v/a*(r0*sin(phi)-r0p*cos(phi))
        gm=1d0
        if(1d0-b1**2-b2**2.gt.0d0)gm=1d0/sqrt(1d0-b1**2-b2**2)
        v1=vx
        v2=vy
        !vx=b1+v1
        !vy=b2+v2
        !formula 7.38, Relativite, C.Semay, Dunod
        bv=b1*v1+b2*v2
        vx=(v1+gm*b1*(1+gm/(gm+1)*bv)) / (gm*(1+bv))
        vy=(v2+gm*b2*(1+gm/(gm+1)*bv)) / (gm*(1+bv))
      enddo
      !print*,x,y,phi-eccphi(2),r,rnll,r/rnll
      !if(mod(kkkk,100).eq.0)
      !. print*,'+++++',nint(phi/pi*180),rnll/r*sqrt(vx**2+vy**2)
      !.,vv,rnll/r,vx,vxbef
      end


c----------------------------------------------------------------------
      real function eflow(epsuu,pressuu,eta,u0uu,u3uu,
     . pi00uu,pi03uu,pi33uu)
c----------------------------------------------------------------------
      common/ceflowvi/pivi,eflowvi
      double precision sh,ch, eps,press,u0,u3,pi00,pi03,pi33
      double precision u0x,u3x,pi00x,pi03x
      sh=sinh(eta)
      ch=cosh(eta)
      eps=epsuu
      press=pressuu
      u0=u0uu
      u3=u3uu
      pi00=pi00uu
      pi03=pi03uu
      pi33=pi33uu
      !--- transform from Bjorken to cms frame ---
      u0x=u0*ch+u3*sh
      u3x=u0*sh+u3*ch
      pi00x=ch**2*pi00+ch*sh*2*pi03+sh**2*pi33
      pi03x=ch**2*pi03+sh*ch*(pi00+pi33)+sh**2*pi03
      !--- compute T^00 and T^03 in the cms frame ---
      t00=(eps+press)*u0x*u0x-press
      t03=(eps+press)*u0x*u3x
      !eflow=(t00*ch-t03*sh)
      eflow=(eps+press)*u0x*u0-press*ch + pi00x*ch - pi03x*sh
      eflowvi=  pivi   *u0x*u0-pivi *ch + pi00x*ch - pi03x*sh
      if(.not.(eflow.le.0..or.eflow.ge.0.))then !NaN catch
        print*,'ERROR NaN:  eps,press,eta,u0,u3,pi00,pi03,pi33 :'
     .  ,eps,press,eta,u0,u3,pi00,pi03,pi33
        stop'\n\n ERROR  eflow NaN catch \n\n'
      endif
      end

c----------------------------------------------------------------------
      function polar(x,y)
c----------------------------------------------------------------------
      pi=3.1415927
      if(abs(x).gt.1.e-6)then
        phi=atan(y/x)
        if(x.lt.0.)phi=pi+phi
        if(phi.lt.0)phi=2*pi+phi
      else
        phi=0.5*pi
        if(y.lt.0)phi=phi+pi
      endif
      polar=phi
      end



