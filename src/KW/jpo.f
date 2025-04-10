C
C  This file is part of EPOS4
C  Copyright (C) 2022 research institutions and authors (See CREDITS file)
C  This file is distributed under the terms of the GNU General Public License version 3 or later
C  (See COPYING file for the text of the licence)
C

c----------------------------------------------------------------------
      subroutine jintpo(iret)
c----------------------------------------------------------------------
c  cut and paste from jintpo version 3236
c----------------------------------------------------------------------
c  "joining interacting Pomerons" 
c      => Plasma droplet formation via clustering  
c
c  1) Core-corona procedure based on cell centered 3D grid via
c  counting string segments (representing the Pomerons) and a energy
c  loss procedure of the segments moving out of the "plasma"
c
c  2) Cluster formation based on the grid and cluster decay
c      (in SR cluster, called at the end of jintpo)
c
c----------------------------------------------------------------------
c Depends on scenario (in particular the decay and its flow settings) :
c iorsdf = 3 ... parametrized flow approach 
c          5 ... core for hydro 
c          6 ... cluster approach for initial conditions 
c iorsdf transmitted via aaa.h
c----------------------------------------------------------------------
#include "aaa.h"
#include "ems.h"
      common/cxyzt/xptl(mxptl+29),yptl(mxptl+29),zptl(mxptl+29)
     * ,tptl(mxptl+29),optl(mxptl+29),uptl(mxptl+29),sptl(mxptl+29)
     *,rptl(mxptl+29,3)
      parameter(mxcl=4000,mxcli=1000)
      parameter(nxjjj=65,nyjjj=65,nzjjj=65)
      integer idropgrid(nxjjj,nyjjj,nzjjj)
     &       ,jdropgrid(nxjjj,nyjjj,nzjjj)
     &       ,kdropgrid(nxjjj,nyjjj,nzjjj)
     &       ,iorgrid(nxjjj,nyjjj,nzjjj)
     &       ,jclu(nzjjj)
     &       ,irep(mxcl),kclu(mxcl)
     &       ,nseg(mxcl),jc(nflav,2)
     &       ,nsegsuj
      common/jintpoc1/idropgrid,jdropgrid,jclu,kclu,nseg,nsegsuj
     &,jc
      real xcell,scell,zcell,delxce
      integer jjj,m1cell,m3cell,nptlb,nptla
      common/jintpoc2/xcell,scell,zcell,delxce
     &,jjj,m1cell,m3cell,nptlb,nptla
      double precision tpro,zpro,ttar,ztar,ttaus,detap,detat
      double precision ena,pza,taa2,pta2,utdble
     .,p5a2,ptff,ptaa
      common/celoss/eloss
      common/cmean/xmean,ymean
      common/cnfrx/nfrx
      real aryz(nyjjj,nzjjj)
      common   /cttaus/tpro,zpro,ttar,ztar,ttaus,detap,detat
      common/cdelzet/delzet,delsce /cvocell/vocell,xlongcell
      common/cranphi/ranphi  /credonoco/kredonoco
      common/cdelcore/delcore,egycore
      integer leading(500),iormn(nzjjj),iormx(nzjjj)
      logical first,two_remnant_sources,one_remnant_source
      common/cwmaxIcoE/wmaxIcoE
      logical lprint 
      save ncntpo
      data ncntpo /0/
      real xmaxico, getIcoXmax
      
      xmaxico=getIcoXmax()

      ncntpo=ncntpo+1
      call utpri('jintpo',ish,ishini,4)
      lprint=nfrx.eq.0
     .       .and.(iorsdf.ne.3.and.iorsdf.ne.6  
     .             .and..not.(iorsdf.eq.5.and.ispherio+ihlle.eq.0)
     .             .or.ish.ge.1)

      if(ish.ge.4)write(ifch,*)'event number',nrevt+1

      iret=0

      fee=yhaha/5.36                 !rapidity range/rap range at RHIC

      m1cell=xmaxico*fxcell
      xcell=xmaxico

      m3cell=21*fee*fzcell
      zcell=10*fee

      if(m1cell.gt.nxjjj )stop'\n\n m1cell too big\n\n'
      if(m3cell.gt.nzjjj )stop'\n\n m3cell too big\n\n'
      scell=zcell*ttaus
      delxce=2*xcell/m1cell
      delsce=2*scell/m3cell
      vocell= delxce * delxce * delsce
      delzet= delsce/ttaus
      xlongcell=1./delsce                  !zcell !delsce/ttaus
      nptla=nptl
      nptlsave=mxptl+1
      nsegsuj=max(3,nsegsu)


      if(lprint)then
      write(ifmt,'(a,$)')'jintpo - ttaus fxcell fzcell nsegce:'
      write(ifmt,'(4f7.2)')ttaus,fxcell,fzcell,dsegce*vocell
      write(ifmt,'(a,$)')'delxce delsce vocell:'
      write(ifmt,*)delxce,delsce,vocell
      write(ifmt,'(a,f6.2,i4,f6.2,i4)')'xcell,m1cell,zcell,m3cell:'
     . ,xcell,m1cell,zcell,m3cell
      esu=0
      do i=1,nptl
      if(istptl(i).eq.0)esu=esu+pptl(4,i)
      enddo
      write(ifmt,'(a,10x,i6,13x,f12.1)')
     . ' +++++Eall+++++(jintpo start)',ng1evt,esu
      endif

c...leading particles

      jj=0
      jjmax=maproj+matarg
      do n=1,nptla
        iaaptl(n)=0
        if(istptl(n).eq.0)then
          iaaptl(n)=1
          do j=1,jj
            if(pptl(4,n).gt.pptl(4,leading(j)))then
              do i=jj+1,j+1,-1
                leading(i)=leading(i-1)
              enddo
              leading(j)=n
              goto 2
            endif
          enddo
          leading(jj+1)=n
    2     continue
          jj=min(jj+1,jjmax)
        endif
      enddo

       jy=0
       do j=1,jjmax
         if(leading(j).le.maproj+matarg)then
           iaaptl(leading(j))=0 !exclude spectator leadings 
         else
           jy=jy+1           !count non spectator leadings (NSL)
         endif 
       enddo

      wleadcore1=(leadcore/10/10000)/1000.
      wleadcore2=(mod(leadcore/10,10000))/1000.
      if(wleadcore1.gt.1.0)stop'ERROR 14072022'
      if(wleadcore2.gt.1.0)stop'ERROR 15072022'
      call remn1paramget(3,f3)
      Z=max(0.,ng1evt-2.)/f3
      wleadcore=wleadcore1+Z*(wleadcore2-wleadcore1) !fraction of NSL to consider

      if(wleadcore.lt.1e-6)then !........exclude all leadings in core
       do j=1,jjmax
        iaaptl(leading(j))=0
        !print*,'TESTjpo EXCLUDE',j,leading(j),pptl(4,leading(j))
       enddo
      else
       !print*,'TESTjpo leadcore,wleadcore=',leadcore,wleadcore,Z,ng1evt
       jy=jy*(1.-wleadcore) !number of candidats to exclude
       jz=0
       do j=1,jjmax
         if(iaaptl(leading(j)).eq.1)then !possible candidates
           jz=jz+1
           if(jz.le.jy)then
             iaaptl(leading(j))=0 !exclude jy candidates
             !print*,'TESTjpo EXCLUDE',jz,jy,nint((1.-wleadcore)*100),'%'
           !else
           !  print*,'TESTjpo KEEP',jz,jy,nint((wleadcore)*100),'%'
           endif 
         !else
         !  print*,'TESTjpo SPECTATOR', j, leading(j).le.maproj+matarg
         endif 
       enddo
      endif

c...compute x,y,z

      if(ish.ge.6)write(ifch,*)'compute x,y,z'
      do n=1,nptla
        !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
        ! the following check is important for  ioclude= 4 or 5
        !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
        if(ioclude.ne.3.and.istptl(n).eq.10)then
        stop'\n\n remnant cluster detected in jintpo\n\n'
        endif
        if(istptl(n).eq.0.or.istptl(n).eq.10)then
          call jtain(n,xptl(n),yptl(n),zptl(n),tptl(n),nnn,0)
          !-----------------------------------------------------------
          !if(abs(taurem).gt.1e-6)then
          !if(ityptl(n).ge.40.and.ityptl(n).le.59.and.nnn.eq.1)then
          !if(ityptl(n).ge.40.and.ityptl(n).le.59)then
          !  if(pptl(4,n).gt.engy/5)then
          !    iaaptl(n)=0
          !  endif
          !endif
          !-----------------------------------------------------------
          call jtaus(zptl(n),tz,sptl(n))
          strap=1e10
          if(zptl(n).lt.0)strap=-1e10
          xpl=tptl(n)+zptl(n)
          xmi=tptl(n)-zptl(n)
          if(xmi.gt.0.0.and.xpl.gt.0.0)then
             strap=0.5*log(xpl/xmi)
          else
             iaaptl(n)=0
          endif
          dezptl(n)=strap       !space-time-rapidity
c           if(iorptl(n).eq.4142)write(ifch,*)'ini',n
c     & ,pptl(4,n),pptl(3,n),idptl(n),xptl(n),yptl(n),zptl(n),sptl(n)
        else
          xptl(n)=0.
          yptl(n)=0.
          zptl(n)=0.
          sptl(n)=0.
          dezptl(n)=1e10
        endif
      enddo
      ntry=0

c...valid particles

      do n=1,nptla
       !random number calls before if's
       !to keep random number sequence for testing
       rdm=rangen()
       if(istptl(n).eq.0)then
         call idtr7( idptl(n) , jc )
       else
         jc(4,1)=0
         jc(4,2)=0
       endif
       if(iaaptl(n).ne.0)then
        if(istptl(n).ne.10)then
          pt2=pptl(1,n)**2+pptl(2,n)**2
          am2tmp=(pptl(4,n)+pptl(3,n))*(pptl(4,n)-pptl(3,n))-pt2
          !~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
          call pthardparent(n,i1,i2,ptmax)
          !if(ptmax.gt.0.)print*,'+++++',n,ptmax,pt2
          !~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
          !if(istptl(n).ne.0.or.ityptl(n).ge.60)then
          if(istptl(n).ne.0)then  !include cluster hadrons
            iaaptl(n)=0
          elseif(abs(idptl(n)).le.100)then
            iaaptl(n)=0         !no fundamental particle (electron...)
         !elseif(ityptl(n).eq.41.or.ityptl(n).eq.51)then
         !  iaaptl(n)=0     !avoid particles coming from remnant droplet decay
         !                       !(already droplet decay products !)
          elseif(abs(am2tmp-pptl(5,n)*pptl(5,n)).gt.1.)then
c            print *,pptl(5,n),pptl(4,n),pptl(3,n),idptl(n)
            iaaptl(n)=0         !to discard off shell particles
          elseif(jc(4,1).ne.0.or.jc(4,2).ne.0)then
          !  iaaptl(n)=-2         !count only energy loss for charmed particles
             iaaptl(n)=0          !charmed particles not counted
          elseif(jc(5,1).ne.0.or.jc(5,2).ne.0)then
          !  iaaptl(n)=-2         !count only energy loss for bottom particles
             iaaptl(n)=0          !bottom particles not counted
          !elseif(pt2.lt.1.e-3)then
          !  iaaptl(n)=0           !to discard slow proton (spectators)
         !~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
         ! elseif(pt2.gt.(1.5*ptoldclu)**2)then
         !   if(maproj.lt.20.or.matarg.lt.20.or.ioquen.eq.0)iaaptl(n)=0
         ! elseif(pt2.gt.(0.5*ptoldclu)**2)then
         !   if(rdm.lt.(sqrt(pt2)-0.5*ptoldclu)/ptoldclu)then
         !     if(maproj.lt.20.or.matarg.lt.20.or.ioquen.eq.0)iaaptl(n)=0
         !   endif
         !~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
          elseif(pt2.gt.ptupp**2)then
            iaaptl(n)=0
         !~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
          elseif(pt2.gt.ptlow**2)then
            iaaptl(n)=-1
         !~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
         ! elseif(ityptl(n).ge.40.and.ityptl(n).le.59
         !.      .and.abs(idptl(n)).gt.1000)then
         !   iaaptl(n)=0   !????????????????????????????????????????????????
         !---------------------------------------------------------------
          endif
        endif
       endif
       if(kredonoco.eq.1)iaaptl(n)=0
      enddo
      if(kredonoco.eq.1)
     . write(ifmt,*)'jintpo: redoing no-core event -> all iaa=0'

c...Start cluster formation

c 8888 continue
c      ntry=ntry+1
c      if(ntry.gt.90)        !do not put more than 100 or change limit for p4mean
c     &call utstop('jintpo: cluster formation impossible ! &')
c      nptl=nptla

      do 1 k=1,m3cell
      iormn(k)=mxptl
      iormx(k)=0
      do 1 j=1,m1cell
      do 1 i=1,m1cell
      idropgrid(i,j,k)=0
      jdropgrid(i,j,k)=0
      kdropgrid(i,j,k)=0
      iorgrid(i,j,k)=0
  1   continue

c...count string segments in cell

      if(ish.ge.6)write(ifch,*)'count string segments in cell'
      do n=1,nptla
        if(iaaptl(n).ne.0)then
          i=1+(xptl(n)+xcell)/delxce
          j=1+(yptl(n)+xcell)/delxce
          k=1+(sptl(n)+scell)/delsce
          if(  i.ge.1.and.i.le.m1cell
     &    .and.j.ge.1.and.j.le.m1cell
     &    .and.k.ge.1.and.k.le.m3cell)then
            idropgrid(i,j,k)=idropgrid(i,j,k)+1
ctp20161220-----------------------------------------------------------
c Core produced in events with few particles because of many particles from the same string fall into the same grid bin
c not to count particles coming from the same source for the density calculation
            ior=iorptl(n)
c            if(ior.gt.0)ior=iorptl(ior)
c            if(ior.gt.0)ior=iorptl(ior)     !get source Pomeron id
            if(iorgrid(i,j,k).ne.ior)then
              kdropgrid(i,j,k)=kdropgrid(i,j,k)+1
              iorgrid(i,j,k)=ior
              iormn(k)=min(iormn(k),ior)
              iormx(k)=max(iormx(k),ior)
            endif
ctp20161220 ----------------------------------------------------------
          endif
        endif
      enddo
c check clusters origins : if all origins are unique in each bin k, not a real core 
      do k=1,m3cell
        if(iormn(k).eq.iormx(k))then
          do j=1,m1cell
            do i=1,m1cell
              kdropgrid(i,j,k)=0
            enddo
          enddo
        endif
      enddo

c...print

      if(ish.ge.5)then
      write(ifch,*)'-------idropgrid--------'
      write(ifch,*)' m3cell,m1cell:',m3cell,m1cell
      j1=0
      j2=0
      ksusu=0 
      do k=1,m3cell
        ksu=0
        j2=0
        do j=1,m1cell
          jsu=0
          do i=1,m1cell
          ksu=ksu+kdropgrid(i,j,k)
          jsu=jsu+idropgrid(i,j,k)
          enddo
          if(jsu.gt.0)j2=j
          if(jsu.eq.0.and.j2.eq.0)j1=j
        enddo
        if(ksu.gt.0)then
          ksusu=ksusu+ksu 
          write(ifch,'(i2,f7.2,20x,2i5)')k,(-scell+(k-0.5)*delsce)/ttaus
     .   ,ksu,ksusu
          write(ifch,'(9x,50i2)')((i),i=1,m1cell)
          do j=j2,j1+1,-1
          y=-xcell+(j-0.5)*delxce
          write(ifch,'(f5.2,2x,50i2)')y,j,(kdropgrid(i,j,k),i=1,m1cell)
          enddo
          write(ifch,'(9x,50i2)')((i),i=1,m1cell)
        endif
      enddo
      endif

      iprgr=0
      kmin=m3cell/2-4
      kmax=m3cell/2+6
      if(iprgr.ge.1)then
      ymin=-xcell
      ymax= xcell
      nybin=m1cell
      zmin=(-scell+(kmin-1)*delsce)/ttaus
      zmax=(-scell+(kmax-0)*delsce)/ttaus
      nzbin=kmax-kmin+1
      do j=1,m1cell
      do k=kmin,kmax
      aryz(j,k-kmin+1)=idropgrid(m1cell/2+1,j,k)
      enddo
      enddo
#if !__BS__ && !__TP__
      call ooHoIniYEta(nyjjj,nzjjj,aryz
     .  ,nybin,ymin,ymax,nzbin,zmin,zmax,'Int   ','segments')
#endif
      endif

c...check high pt segments

      !to use this part one has to define:
      !...  1 = valid
      !... -1 = valid but high pt
      !...  0 = not valid
      
      !TEST
      !n=igetNpom()
      !call getRnpom(1.0,Z)!Z=n/max_n/f
      !Z2=max(0.,ng1evt-2.)/35
      !write(*,'(a,2f10.2)')'TEST-Z',Z,Z2
      !ENDTEST 

      esu=0
      do i=1,nptla
      if(istptl(i).eq.0.or.istptl(i).eq.3)esu=esu+pptl(4,i)
      if(istptl(i).eq.25)rinptl(i)=0.
      if(istptl(i).eq.25)qsqptl(i)=0.
      if(istptl(i).eq.25)itsptl(i)=0
      enddo
      if(lprint)
     . write(ifmt,'(a,16x,f12.1)')
     . ' +++++Eall+++++(jintpo start check highpt)',esu

      if(ish.ge.6)write(ifch,*)'check high pt segments'
      if(qufac.gt.0.0)then
       ein=0
       elo=0
       ncore=0
       nesc=0
       nesc0=0
       encore=0
       enloss=0
       do n=1,nptla
        !pt25=sqrt(pptl(1,n)**2+pptl(2,n)**2)
        !if(istptl(n).eq.25)write(ifmt,*)
        !.   '25parton',n,istptl(n),pt25
        if(iaaptl(n).le.-1)then
          i=1+(xptl(n)+xcell)/delxce
          j=1+(yptl(n)+xcell)/delxce
          k=1+(sptl(n)+scell)/delsce
          if(   i.ge.1.and.i.le.m1cell
     &     .and.j.ge.1.and.j.le.m1cell
     &     .and.k.ge.1.and.k.le.m3cell)then
           iescape=0  
           timu=tptl(n)
           timo=xorptl(4,n)
           xa=xptl(n)
           ya=yptl(n)
           eta=dezptl(n)
           pz=pptl(3,n)
           en=pptl(4,n)
           taa2=pptl(1,n)**2+pptl(2,n)**2+pptl(5,n)**2
           taa=sqrt(taa2)
           pza=pz
           ena=sqrt(taa2+pza**2)
           enaxx=en
           p5a2=utdble(pptl(5,n)**2)
           ppa2=pptl(1,n)**2+pptl(2,n)**2
           pta2=utdble(ppa2)
           ppa2=ppa2+pptl(3,n)**2
           ptaa=sqrt(max(1d-6,pta2))
           phia=polar( pptl(1,n) , pptl(2,n) )
           r12=1e30
           p12=0.
           rap12=0
           phi12=0
           call index2521(n,i25first,i21last)
           if(i25first.gt.0)then
            do ij=i25first,i21last
             call jet2521(ij,ibig,phij,rapj,ppj)
             if(ibig.eq.1.or.ij-i25first.le.1)then  
              phia1=phia-phij
              if(phia1.lt.-3.14159)phia1=phia1+2*3.14159
              if(phia1.gt. 3.14159)phia1=phia1-2*3.14159
              amt=sqrt(max(1d-6, p5a2+pta2 ))
              ppa=sqrt(max(1e-6, ppa2 ))
              !pseudorapidity:
              rap1=sign(1.,pz)*log((ppa+abs(pz))/ptaa) - rapj
              r1=sqrt(phia1**2+rap1**2)
              if(r1.le.r12)then
                r12=min(r1,r12)
                p12=ppj
                rap12=rap1
                phi12=phia1
              endif
              !if(ptaa.gt.4)
              !.write(ifmt,'(a,3i7,f8.2,2i7,2f8.2)')
              !.'jetcone',n,iorptl(n),istptl(iorptl(n)),ptaa,ij
              !.,istptl(ij),r1,r12
             endif 
            enddo
           endif
           vz=pza/ena
           vx=pptl(1,n)/ena  !pp cms working frame
           vy=pptl(2,n)/ena  !pp cms working frame
           !------------------------------------------
           !we consider a purely transverse motion of the particle 
           !  in the frame where its vz=0, assuming a "fluid" at rest
           !  in this frame (Bjorken scenario, whis is roughly realized)
           !it is not necessary to boost into this frame, so we work with
           !  vx,vy in the pp cms working frame
           !-------------------------------------------
           !print*,'+++++++',eta,vz,pz/en
           delen=0
           delll=0
           jcone=0
           if(vx.ne.0.0.or.vy.ne.0.0)then
             if(abs(vx).ge.abs(vy))then
               ica=1
               rat=vy/vx
               is=sign(1.,vx)
               l=i
               va=xa
               wa=ya
             else
               ica=2
               rat=vx/vy
               is=sign(1.,vy)
               l=j
               va=ya
               wa=xa
             endif
             if(is.eq.-1)then
               imax=1
             else !if(is.eq. 1)
               imax=m1cell
             endif
             delen=0
             nsegpa=-1   !not to count the particle it-self at the starting point
             vr=sqrt(vx**2+vy**2)
             dll=delxce*sqrt(1+rat**2)  
             dlt=dll/vr
             do lu=l,imax,is
              vu=va+(lu-l)*delxce
              wu=wa+rat*(lu-l)*delxce
              mu=1+(wu+xcell)/delxce
              if(mu.ge.1.and.mu.le.m1cell)then
                if(ica.eq.1)then
                  ix=lu
                  jx=mu
                else
                  ix=mu
                  jx=lu
                endif
                !ctp nsegpa=nsegpa+idropgrid(ix,jx,k)
                !ctp dens=float(idropgrid(ix,jx,k))
                nsegpa=nsegpa+kdropgrid(ix,jx,k) !ctp20161220
                dens=float(kdropgrid(ix,jx,k))   !ctp20161220
                n3=3 !triple counting
                !if(iorsdf.eq.3)n3=2        !to compensate some missing absorption in effective core
                do ixx=ix-1,ix+1
                do jxx=jx-1,jx+1
                  if((ixx.ne.ix.or.jxx.ne.jx)
     .            .and. ixx.ge.1 .and. jxx.ge.1 
     .            .and. ixx.le.m1cell .and. jxx .le.m1cell)then 
                    !ctp dens=dens+float(idropgrid(ixx,jxx,k))
                    !ctp nsegpa=nsegpa+idropgrid(ixx,jxx,k)
                    dens=dens+float(kdropgrid(ixx,jxx,k)) !ctp20161220
                    nsegpa=nsegpa+kdropgrid(ixx,jxx,k)    !ctp20161220
                  endif
                enddo
                enddo
                if(idropgrid(ix,jx,k).gt.0)delll=delll+dll
                timu=timu+dlt
                tidi=timu-timo
                call fieloss(ptaa,dens,r12,p12,fiel,jcone) !,tidi)
                if(ish.ge.6)then
                !~~~~~~~~~~~~~~~~~~~~~~~~~~
                if(ptaa.gt.7.0)then !.and.dens.gt.1e-5.and.fiel.gt.1e-5 
                write(ifmt,
     .          '(a,f8.2,3x,2f8.2,3x,3f6.2,i3,3x,f11.5,i7,3x,6i5)')
     .          'HIGHPT',ptaa,timo,timu,rap12,phi12,r12,jcone,fiel
     .          ,ikoevt,n,idptl(n),istptl(n),ityptl(n),i25first,i21last
                endif
                !~~~~~~~~~~~~~~~~~~~~~~~~~~
                endif            
                delen=delen+fiel*dll
                if(ish.ge.5)then
                  write(ifch,'(i6,3x,2i3,4x,2i3,4x,2i3,4x,i3,2f7.3
     .           ,4x,i3,e12.3,i7,2e12.3)')
     .            n,k,ica,i,j,ix,jx,nint(dens),delen
     .         ,sqrt(pta2)-delen,nsegpa,dsegce*n3*vocell,jcone,r12,tidi
                endif  
                ! write(ifmt,'(7x,a,3f6.2,5x,2i3,2x,2(2f6.2,2x),f8.2)')'TRAJ   '
                !.          ,-xcell+(ix-0.5)*delxce
                !.          ,-xcell+(jx-0.5)*delxce
                !.          ,-scell+(k -0.5)*delsce
                !.          ,idropgrid(i,j,k),idropgrid(ix,jx,k)
                !.          ,dens,dens**(3./8.),del*escale,qu,dliz
              endif
             enddo
             if(float(nsegpa).le.dsegce*n3*vocell)delen=0  ! <=================== !!!
             if(iaaptl(n).eq.-2.and.jcone.ne.1)then
               !for HQ, count only energy loss in a jet cone
               delen=0.
               iaaptl(n)=0
             endif
           endif
           ptold=sngl(sqrt(pta2))
           ptnew=ptold-delen
           if(iaaptl(n).eq.-2)then
             if(ptnew.le.0.)then
               ptnew=rangen()*0.5
             endif
           endif
           if(ptnew.gt.0)then
             iescape=1
             p1i=pptl(1,n)
             p2i=pptl(2,n)
             p3i=pptl(3,n)
             p4i=pptl(4,n)
             p1=ptnew*pptl(1,n)/ptold
             p2=ptnew*pptl(2,n)/ptold
             taf=sqrt(p1**2+p2**2+pptl(5,n)**2)
             p3=taf/taa*pza
             p4=taf/taa*ena
             !check if momentum is large enough taking into account the particle mass
             !   if(iorsdf.eq.3.and.delen/p4i.gt.1.e-5.and.
             !  . ptnew**2.lt.yradpi*p5a2)iescape=0
             
             !with effective core low momentum particles can't be absorbed by medium, 
             !do it here (important to get same number of segment in real hydro effective core (empirical formula ...))
             if(iorsdf.eq.3.and.delen.gt.0.
     .       .and.ptnew.lt.max(0.,pptl(5,n)-0.5*(ptnew+ptold)))iescape=0
           endif
           if(iescape.eq.1)then  
              nptlsave=nptlsave-1 
              iaaptl(nptlsave)=n
              idptl(nptlsave)=idptl(n)
              do l=1,4
                pptl(l,nptlsave)=pptl(l,n)
              enddo
              iaaptl(n)=0
              idropgrid(i,j,k)=idropgrid(i,j,k)-1
              if(idropgrid(i,j,k).lt.0)stop'jintpo: not possible.    '
              pptl(1,n)=p1
              pptl(2,n)=p2
              pptl(3,n)=p3
              pptl(4,n)=p4
              call updateXor(n,n)
              elo= elo+ena-sqrt(ena**2-pta2+ptnew**2)  
              !write(ifch,*)n,'      escape'
              i25=i25first+1
              if(i25first.ne.0)then
                ptlim=50
                if(pptl(1,i25)**2+pptl(2,i25)**2.gt.ptlim**2
     .          .or.pptl(1,i25-1)**2+pptl(2,i25-1)**2.gt.ptlim**2)then
                  ptff=sqrt(pptl(1,n)**2+pptl(2,n)**2)
                  do ih=i25-1,i25
                    pth=sqrt(pptl(1,ih)**2+pptl(2,ih)**2)
                    phi=polar( pptl(1,ih) , pptl(2,ih) )
                    dphi=abs(phia-phi)
                    if(dphi.gt.3.14159)dphi=2*3.14159-dphi
                    if(dphi.lt.0)stop'\n\n  19062012 \n\n'
                    if(dphi.lt.1.)rinptl(ih)=rinptl(ih)+ptaa-ptff
                    if(dphi.lt.1.)qsqptl(ih)=qsqptl(ih)+delll
                    if(dphi.lt.1.)itsptl(ih)=itsptl(ih)+1
                    write(ifmt,'(a,2i6,2i3,7f7.2,i3)')
     .               ' ++++25parton++++'
     .               ,n,ih,idptl(ih),istptl(ih)
     .               ,ptff,ptaa,pth,dphi,rinptl(ih),delll
     .               ,qsqptl(ih),itsptl(ih)
                  enddo
                  write(ifmt,'(a)')' ++++25parton++++'
                endif
              endif             !~~~~ i25 ~~~~
              if(iorsdf.ne.5.and.delen.gt.1e-3)call fpartloss(
     .             i,j,k,n,p1i-p1,p2i-p2,p3i-p3,p4i-p4)
              nesc=nesc+1
              enloss=enloss+ena-pptl(4,n)
              if(ish.ge.4)
     .        write(ifmt,'(a,2i6,2f10.1,2f10.2,i5)')'******CORONA******'
     .        ,n,nesc,pptl(4,n),ena,ptnew,ptold,jcone
           else
              ncore=ncore+1
              encore=encore+pptl(4,n)
              if(ish.ge.4)
     .        write(ifmt,'(a,2i6,2f10.1,2f10.2,i5)')'*******CORE*******'
     .        ,n,ncore,pptl(4,n),ena,ptnew,ptold,jcone
                 !~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
                 ijetconecheck=0  ! 0 1 2
                 if(ijetconecheck.ge.1)then 
                 if(ptold.gt.6)then                     !0
                 r12=1e30
                 p12=0.
                 i777=0
                 if(i25first.gt.0)then                  !1
                 do ij=i25first,i21last
                 call jet2521(ij,ibig,phij,rapj,ppj)
                 if(ibig.eq.1.or.ij-i25first.le.1)then  !2
                 phia1=phia-phij
                 if(phia1.lt.-3.14159)phia1=phia1+2*3.14159
                 if(phia1.gt. 3.14159)phia1=phia1-2*3.14159
                 amt=sqrt(max(1d-6, p5a2+pta2 ))
                 ppa=sqrt(max(1e-6, ppa2 ))
                 !pseudorapidity:
                 rap1=sign(1.,pz)*log((ppa+abs(pz))/ptaa) - rapj
                 r1=sqrt(phia1**2+rap1**2)
                 if(r1.lt.r12)then
                   r12=min(r1,r12)
                   p12=ppj
                 endif
                 i777=i777+1
                 j777=111
                 if(ijetconecheck.ge.2)then 
                 write(ifmt,'(a,3i7,f8.2,2i7,2f8.2)')
     .           'jetcone 111111',n,iorptl(n)
     .           ,istptl(iorptl(n)),ptaa,ij
     .           ,istptl(ij),r1,r12
                 endif
                 endif                                  !2
                 enddo
                 endif                                  !1
                 if(i777.eq.0)then
                 j777=222
                 if(ijetconecheck.ge.2)then 
                 write(ifmt,'(a,2i7,3x,2i7,3x,2i7)')
     .           'jetcone 222222',n,iorptl(n)
     .           ,istptl(iorptl(n)),ityptl(iorptl(n))
     .           ,i25first,i21last
                 endif
                 endif  
                 write(ifmt,'(a,2i7,2f8.2,i7,i3,i4,a,i7)')
     .           'jetcone',jcone,nsegpa,ptold,delen
     .           ,istptl(iorptl(n)),ityptl(iorptl(n)),j777
     .           ,'      ID',idptl(n)
                 endif                                   !0
                 endif
                 !~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
              iaaptl(n)=1
              ein=ein+pptl(4,n)
              !~~~~~~~~~~~~~~~~~~~~~~~~~~
              !if(sqrt(pta2).gt.7.)then
              !write(*,'(a,f10.4,i5,$)')'HIGHPT CORE ',ptaa,jcone
              !print*,r12,i25,idptl(n),istptl(n),ityptl(n)
              !endif
              !~~~~~~~~~~~~~~~~~~~~~~~~~~
           endif
          else
            if(ish.ge.4)
     .      write(ifmt,'(a,2i6,2f10.1,2f10.2,i5)')'********OUT********'
     .      ,n,0,pptl(4,n)
            iaaptl(n)=0
          endif
        elseif(istptl(n).eq.0.and.iaaptl(n).eq.0)then
          if(ish.ge.4)
     .    write(ifmt,'(a,2i6,2f10.1,2f10.2,i5)')'***FORCED CORONA**'
     .    ,n,0,pptl(4,n)
          nesc0=nesc0+1
        elseif(istptl(n).eq.0.and.iaaptl(n).eq.1)then
          if(ish.ge.4)
     .    write(ifmt,'(a,2i6,2f10.1,2f10.2,i5)')'****FORCED CORE***'
     .    ,n,0,pptl(4,n)
          ncore0=ncore0+1
          encore=encore+pptl(4,n)
        endif
       enddo !~~~~~~ end n loop ~~~~~~~~~
       if(lprint)
     . write(ifmt,'(a,i6,3x,a,i6,4x,2(a,f10.2))')'Ncore=',ncore+ncore0
     .  ,'   Ncorona=',nesc+nesc0,'   Ecore=',encore,'   Eloss=',enloss
      else
       do n=1,nptla
        if(iaaptl(n).eq.-1)then
          iaaptl(n)=0
        endif
       enddo
      endif

      !ncore=0
      !nesc=0
      !do n=1,nptla
      !  if(istptl(n).eq.0.and.iaaptl(n).eq.1)ncore=ncore+1
      !  if(istptl(n).eq.0.and.iaaptl(n).eq.0)nesc=nesc+1
      !enddo
      !write(ifmt,'(a,i6,a,i6)')'Ncore=',ncore,'   Nescape=',nesc

c...identify particular fluids

      !KW 2010 still relevant !!
      !small fluids may have eloss >> ecore
      !which gives completely artificial results (after correction in ico)
      npomincor=3
      nmidmincor= -1 !40
      rref=0.30
      npom=0
      nmid=0
      nremn1=0
      nremn2=0
      nremn3=0
      npomer=0
      xm=0
      ym=0
      x2m=0
      y2m=0
      npat=npjevt+ntgevt
      rcox=0      
      one_remnant_source = .false.
      two_remnant_sources = .false.
      k1=0
      k2=0
      iremn1=0
      do n=1,nptla
       if(istptl(n).eq.30.or.istptl(n).eq.31)npom=npom+1
       k=1+(sptl(n)+scell)/delsce
       if(k.ge.m3cell/2-10 .and. k.le.m3cell/2+10
     .  .and.iaaptl(n).eq.1 )then
         nmid=nmid+1
         !~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
         ! identify core 
         !  from one remnant source
         !    one_remnant_source = .true.
         !  from two remnant sources with large rap difference:
         !    two_remnant_sources = .true.
         !~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
         if(ityptl(n).le.39)then
           npomer=npomer+1
         endif
         if(nremn2.gt.0.and.iorptl(n).ne.iremn2
     .   .and.ityptl(n).ge.40.and.ityptl(n).le.59)then
           nremn3=nremn3+1
         endif
         if(nremn3.eq.0.and.nremn1.gt.0.and.iorptl(n).ne.iremn1
     .   .and.ityptl(n).ge.40.and.ityptl(n).le.59)then
           nremn2=nremn2+1
           iremn2=iorptl(n)
           k2=k
         endif
         if(nremn2.eq.0.and.ityptl(n).ge.40.and.ityptl(n).le.59)then
           nremn1=nremn1+1
           iremn1=iorptl(n)
           k1=k
         endif
         one_remnant_source =
     .   nremn1.gt.0.and.nremn2.eq.0.and.nremn3.eq.0.and.npomer.eq.0
         two_remnant_sources =
     .   nremn1.gt.0.and.nremn2.gt.0.and.nremn3.eq.0.and.abs(k2-k1).gt.3
     .   .and.npomer.eq.0
         !print*,'nremn1,...',nremn1,nremn2,nremn3,npomer,iorptl(n)
         !.   ,k2-k1,one_remnant_source,two_remnant_sources
         !~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
         xm=xm+xptl(n)
         ym=ym+yptl(n)
         x2m=x2m+xptl(n)**2
         y2m=y2m+yptl(n)**2
       endif
      enddo
      if(nmid.gt.0)then
        xm=xm/nmid
        ym=ym/nmid
        x2m=x2m/nmid
        y2m=y2m/nmid
        rcox=(x2m-xm**2)+(y2m-ym**2)
        if(rcox.ge.0.)then
          rcox=sqrt(rcox)
        else
          write(ifmt,*)' rcox**2<0',rcox
          rcox=0.
        endif
      endif
      rcox=max(rcox,rref)
      icox=rcox*100
      rax=(rcox/rref)**2 !ratio of areas
      !kcox=ncox/rax
      !kcox=nmid/rax
      kcox=nmid/npat*2
      iundo=0
      if(ncore+ncore0.lt.nsegmincore)iundo=1
      if(enloss.gt.encore*ratiomaxcore)iundo=1
      !if(npom.le.npomincor)iundo=1
      !if(ncore.le.ncomincor)iundo=1
      !if(kcox.le.ncomincor)iundo=1
      !if(nmid.le.nmidmincor
      !. .or. nmid.le.nmidmincor*4.and.kcox.le.nmidmincor )iundo=1
      !if(nmid.le.nmidmincor)iundo=1
      if(one_remnant_source)iundo=1
      if(two_remnant_sources)iundo=1
      if(nfrx.gt.0.and.wmaxIcoE.le.1e-5)iundo=1 !if wmax zero for nfrx=0 
      if(iphsd.ne.0)iundo=0

c...undo core (transform core to corona) 

      if(iundo.eq.1)then !KW2010
        !write(ifch,'(a,5x,5i4)')'small flu',npat,npom,icox,nmid,kcox
        do k=1,m3cell
        do j=1,m1cell
        do i=1,m1cell
          idropgrid(i,j,k)=0
        enddo
        enddo
        enddo
        do i=1,nptla
          iaaptl(i)=0
          if(iorsdf.ne.5)call fpartdel(i,idum)
        enddo
        do nsa=mxptl,nptlsave,-1
          n=iaaptl(nsa)
          if(idptl(n).ne.idptl(nsa))
     .    stop'#######  ERROR 03062015########'
          do l=1,4
            pptl(l,n)=pptl(l,nsa) 
          enddo
        enddo
        nptl=nptla
        eloss=0.
        if(iorsdf.eq.5)call iniCore(2) !needed nevertheless, to define tables
        goto 1001
      !else
      !  write(ifch,'(a,55x,5i4)')'small flu',npat,npom,icox,nmid,kcox
      endif

c...compute eloss

      eloss=esu
      esu=0
      do i=1,nptla
        if(istptl(i).eq.0.or.istptl(i).eq.3)esu=esu+pptl(4,i)
        if(i.gt.1)then
        if(istptl(i).eq.25.and.(rinptl(i).gt.0.001
     .    .or.(istptl(i-1).eq.25.and.rinptl(i-1).gt.0.001)
     .    .or.(istptl(i+1).eq.25.and.rinptl(i+1).gt.0.001)))then
          pth=sqrt(pptl(1,i)**2+pptl(2,i)**2)
          if(itsptl(i).ne.0)then
          qsqptl(i)=qsqptl(i)/itsptl(i)
          else
          qsqptl(i)=0
          endif
          write(ifmt,'(a,i6,2i3,3f7.2)')' ++++25parton++++'
     .   ,i,idptl(i),istptl(i)
     .   ,pth,rinptl(i),qsqptl(i)
        endif
        endif
      enddo
      if(lprint)
     . write(ifmt,'(a,23x,f12.1)')
     . ' +++++Eall+++++(jintpo after eloss)',esu
      eloss=eloss-esu
      !write(ifch,*)'+++++ein,elo,eloss ',ein,elo,eloss

      nptla=nptl  !to take into account the "fake" particles in effective core
 8888 continue
      nsegsuj=max(3,nint(float(nsegsu)/1.02**min(90,ntry)))
      ntry=ntry+1
      if(ntry.gt.999)
     &call utstop('jintpo: cluster formation impossible ! &')
      nptl=nptla

c...identify clusters

      if(ish.ge.6)write(ifch,*)'identify clusters'
      do k=1,m3cell !~~~~~~k-loop
        jjj=0
        !here one loops over all  cells in the transverse plane
        !(i,j correspond to x,y) and one associates numbers jjj=1,2,3...
        !to the cells, such that at the end different numbers refer
        !to different unconnected clusters. A cluster is defined to be
        !a sequence of neighboring cells. 
        do j=1,m1cell
          ifirst=0
          first=.true.
          do i=1,m1cell
            if(idropgrid(i,j,k).ge.kigrid)then
              if(first)then
                ifirst=i
                jjj=jjj+1
                if(jjj.gt.mxcl)stop'jintpo: mxcl too small.   '
                irep(jjj)=0
                jj=jjj
                first=.false.
              endif
              jdropgrid(i,j,k)=jj
              if(j.gt.1)then
               if(jdropgrid(i,j-1,k).ne.0)then
                jjo=jdropgrid(i,j-1,k)
                if(jjo.lt.jj)then
                  if(jj.eq.jjj)jjj=jjj-1
                  jjx=jj
                  jj=jjo
                  do ii=ifirst,i
                    jdropgrid(ii,j,k)=jj
                    if(jdropgrid(ii,j-1,k).eq.jjx)then
                      if(jjx.gt.jjj)jjj=jjj+1
                        jja=jjx
                      jjb=jj
  90                  continue
                      if(irep(jja).eq.0.or.irep(jja).eq.jjb)then
                        irep(jja)=jjb
                      else
                        mn=min(irep(jja),jjb)
                        mx=max(irep(jja),jjb)
                        irep(jja)=mn
                        jja=mx
                        jjb=mn
                        goto90
                      endif
                    endif
                  enddo
                elseif(jdropgrid(i,j-1,k).gt.jj)then
                  irep(jjo)=jj
                endif
               endif
              endif
            else
              jdropgrid(i,j,k)=0
              first=.true.
            endif
          enddo
        enddo
        !~~~~cluster jj ---> cluster irep(jj)
        do jj=jjj,1,-1
         if(irep(jj).ne.0)then
           do j=1,m1cell
             do i=1,m1cell
               if(jdropgrid(i,j,k).eq.jj)jdropgrid(i,j,k)=irep(jj)
             enddo
           enddo
         endif
        enddo
        !~~~~~remove empty cluster indices
        jjjx=jjj
        jjj=0
        jj=0
        do jjx=1,jjjx
         if(irep(jjx).eq.0)then
           jj=jj+1
           jjj=jjj+1
         else
           do j=1,m1cell
             do i=1,m1cell
               if(jdropgrid(i,j,k).gt.jj)
     &           jdropgrid(i,j,k)=jdropgrid(i,j,k)-1
             enddo
           enddo
         endif
        enddo
        !~~~~~
        jclu(k)=jjj
      enddo !~~~~~~~~~~~~~~~~~ END k-loop

c...absolute clusters numbering (for all k)

      if(ish.ge.6)write(ifch,*)'absolute clusters numbering'
      jjj=jclu(1)
      do k=2,m3cell
        ncellmax=0
        do j=1,m1cell
          do i=1,m1cell
            if(jdropgrid(i,j,k).gt.0)then
              jdropgrid(i,j,k)=jdropgrid(i,j,k)+jjj
            endif
          enddo
        enddo
        jjj=jjj+jclu(k)
      enddo

c...set effective flow

      if(iorsdf.eq.3.or.iorsdf.eq.6)call SetPFE() !Parametrized Fluid Expansion

      nptlb=nptl+jjj

c...print
      
      if(ish.ge.5)then
      write(ifch,*)'-------jdropgrid--------'
      write(ifch,*)' m3cell,m1cell:',m3cell,m1cell
      j1=0
      j2=0
      do k=1,m3cell
        ksu=0
        do j=1,m1cell
        jsu=0
        do i=1,m1cell
        ksu=ksu+idropgrid(i,j,k)
        jsu=jsu+idropgrid(i,j,k)
        enddo
        if(jsu.gt.0)j2=j
        if(jsu.eq.0.and.j2.eq.0)j1=j
        enddo
        if(ksu.gt.0)then
          write(ifch,'(i2,f7.2)')k,(-scell+(k-0.5)*delsce)/ttaus
          write(ifch,'(9x,50i2)')((i),i=1,m1cell)
          do j=j2,j1+1,-1
          y=-xcell+(j-0.5)*delxce
          write(ifch,'(f5.2,2x,50i2)')y,j,(jdropgrid(i,j,k),i=1,m1cell)
          enddo
          write(ifch,'(9x,50i2)')((i),i=1,m1cell)
        endif
      enddo
      endif

      esu=0
      do i=1,nptla
        if(istptl(i).eq.0.or.istptl(i).eq.3)esu=esu+pptl(4,i)
      enddo
      if(lprint)
     . write(ifmt,'(a,11x,f12.1)')
     .' +++++Eall+++++(jintpo before marking segments)',esu

c...mark segments going into in clusters, count them

      if(ish.ge.6)write(ifch,*)'mark segments going into in clusters'

      !nlost=0
      !nlow=0
      !do n=1,nptla
      !  if(istptl(n).eq.0.and.iaaptl(n).eq.1)ncore=ncore+1
      !  if(istptl(n).eq.0.and.iaaptl(n).eq.0)nesc=nesc+1
      !enddo
      !write(ifmt,'(a,i6,a,i6)')'Ncore=',ncore,'   Nescape=',nesc

      do 96 n=1,nptla
        if(iaaptl(n).eq.0)goto96
        i=1+(xptl(n)+xcell)/delxce
        j=1+(yptl(n)+xcell)/delxce
        k=1+(sptl(n)+scell)/delsce
        if(    i.ge.1.and.i.le.m1cell
     &          .and.j.ge.1.and.j.le.m1cell
     &          .and.k.ge.1.and.k.le.m3cell)then
          jj=jdropgrid(i,j,k)
          if(jj.gt.0)then
              istptl(n) = 3
              nseg(jj)=nseg(jj)+1
          else
            if(iorsdf.ne.5)call fpartdel(n,nlow)
            nlow=nlow+1
            if(ish.ge.5)
     .      write(ifch,*)'******CORE but NOCLUSTER******',n
            write(ifch,*)i,j,k,jdropgrid(i,j,k)
            istptl(n) = 0   !restore istptl in case of more than 1 try
          endif
        else
          nlost=nlost+1
          if(ish.ge.5)
     .    write(ifch,*)'******CORE but OUT******',n
          istptl(n) = 0   !restore istptl in case of more than 1 try
        endif
  96  continue

      !ncore=0
      !nesc=0
      !do n=1,nptla
      !  if(istptl(n).eq.3)ncore=ncore+1
      !  if(istptl(n).eq.0)nesc=nesc+1
      !enddo
      !write(ifmt,'(a,i6,a,i6,a,i6,a,i6)')'Ncore=',ncore,'   Nescape='
      !. ,nesc,'   Nlost=',nlost,'   Nlow=',nlow

c...add segments moving towards clusters

      if(ish.ge.6)write(ifch,*)'add segments moving towards clusters'
      if(iocluin.eq.1)then
      do 93 n=1,nptla
        if(iaaptl(n).eq.0)goto93
        ihit=0
        i=1+(xptl(n)+xcell)/delxce
        j=1+(yptl(n)+xcell)/delxce
        k=1+(sptl(n)+scell)/delsce
        if(    i.ge.1.and.i.le.m1cell
     &          .and.j.ge.1.and.j.le.m1cell
     &          .and.k.ge.1.and.k.le.m3cell)then
          jgr=jdropgrid(i,j,k)
c no cluster at particle position
          if(jgr.eq.0)then
c look for custers on the way of the particle (toward center)
           if(i.ge.m1cell/2)then
            isi=-1
           else
            isi=1
           endif
           if(j.ge.m1cell/2)then
            jsi=-1
           else
            jsi=1
           endif
            do ii=i,i+2*isi,isi
             do jj=j,j+2*jsi,jsi
               if(.not.(ii.eq.i.and.jj.eq.j))then
                if(ii.ge.1.and.ii.le.m1cell
     .          .and.jj.ge.1.and.jj.le.m1cell)then
                 jg=jdropgrid(ii,jj,k)
c cluster jg is between particle and (x,y)=(0,0)
                 if(jg.gt.0)then
                  !if(nseg(jg).gt.50)then  !??????????????????????????????????????
                   x=xptl(n)
                   y=yptl(n)
                   vrad=( x*pptl(1,n)/pptl(4,n)+y*pptl(2,n)/pptl(4,n))
                    if(vrad.lt.0.)then
c vrad<0 means particles is going in the direction  of the (0,0) axis
c so cluster is on the way of the particle
                     ihit=1
                     goto94
                    endif
                  !endif
                 endif
                endif
               endif
             enddo
            enddo
          endif
        endif
   94   continue
        if(ihit.eq.1)then
         delx=delxce*(ii-i)
         dely=delxce*(jj-j)
c particle position at impact point (cell dimension)
         xn=xptl(n)+delx
         yn=yptl(n)+dely
         ix=1+(xn+xcell)/delxce
         jx=1+(yn+xcell)/delxce
         jgrx=jdropgrid(ix,jx,k)
c ix,jx,jgrx should be the same as ii,jj,jg
         if(jgrx.gt.0)then
c new position of the segment
          xptl(n)=xn
          yptl(n)=yn
          istptl(n) = 3
c no new cluster formed but jgrx is updated
          nseg(jgrx)=nseg(jgrx)+1
c move segment
          idropgrid(i,j,k)=idropgrid(i,j,k)-1        !remove segment from original position
          idropgrid(ix,jx,k)=idropgrid(ix,jx,k)+1    !add segment to new position
         endif
        endif
  93  continue
      endif

      esu=0
      do i=1,nptla
        if(istptl(i).eq.0.or.istptl(i).eq.3)esu=esu+pptl(4,i)
      enddo
      if(lprint)
     . write(ifmt,'(a,5x,f12.1)')
     .' +++++Eall+++++(jintpo after moving towards clusters)',esu

      !ncore=0
      !nesc=0
      !do n=1,nptla
      !  if(istptl(n).eq.3)ncore=ncore+1
      !  if(istptl(n).eq.0)nesc=nesc+1
      !enddo
      !write(ifmt,'(a,i6,a,i6)')'Ncore=',ncore,'   Nescape=',nesc

c...add segments to avoid holes

      iohole=0
      if(ish.ge.6)write(ifch,*)'add segments from holes'
      if(iohole.eq.1)then
      do 83 n=1,nptla
        if(iaaptl(n).eq.0)goto83
        ihit=0
        i=1+(xptl(n)+xcell)/delxce
        j=1+(yptl(n)+xcell)/delxce
        k=1+(sptl(n)+scell)/delsce
        if(    i.ge.1.and.i.le.m1cell
     &          .and.j.ge.1.and.j.le.m1cell
     &          .and.k.ge.1.and.k.le.m3cell)then
          jgr=jdropgrid(i,j,k)
          if(jgr.eq.0)then
            isgi=1
            if(rangen().gt.0.5)isgi=-1
            isgj=1
            if(rangen().gt.0.5)isgj=-1
            isgk=1
            if(rangen().gt.0.5)isgk=-1
            do ii=i-isgi,i+isgi,isgi
            do jj=j-isgj,j+isgj,isgj
            do kk=k-isgk,k+isgk,isgk
              if(  ii.ge.1.and.ii.le.m1cell
     .        .and.jj.ge.1.and.jj.le.m1cell
     .        .and.kk.ge.1.and.kk.le.m3cell)then
                if(jdropgrid(ii,jj,kk).gt.0)then
                  nplus=idropgrid(ii,jj,kk)+1
                  if(nplus.gt.kigrid)then
                    ihit=1
                    goto84
                  endif
                endif
              endif
            enddo
            enddo
            enddo
          endif
        endif
   84   continue
        if(ihit.eq.1)then
          istptl(n) = 3
          idropgrid(i,j,k)=idropgrid(i,j,k)+1
          jdropgrid(i,j,k)=jjj+1
        endif
  83  continue
      if(iorsdf.ne.5)stop'\n\n  update avoid hole proc\n\n'
      !the added particle should be put into a cluster
      endif

    !~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    !~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    !~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    !~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    !~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

      if(iorsdf.eq.5)then
        esu=0
        do i=1,nptla
          if(istptl(i).eq.0.or.istptl(i).eq.3)esu=esu+pptl(4,i)
        enddo
        if(lprint)
     .  write(ifmt,'(a,14x,f12.1)')
     .  ' +++++Eall+++++(jintpo after avoiding holes)',esu
      endif

      ncore=0
      nesc=0
      xcoremean=0
      ycoremean=0
      xcore2mean=0
      ycore2mean=0
      delxcore=0
      delycore=0
      do n=1,nptla
        if(istptl(n).eq.3)then
          ncore=ncore+1
          xcoremean=xcoremean+xptl(n)
          ycoremean=ycoremean+yptl(n)
          xcore2mean=xcore2mean+xptl(n)**2
          ycore2mean=ycore2mean+yptl(n)**2
          if(iorsdf.eq.5)then
            if(idptl(n).lt.1e4)then
              istptl(n) = 5
            else
              istptl(n) = 15
            endif
          endif
        else
          if(istptl(n).eq.0)nesc=nesc+1
        endif
      enddo
      if(ncore.gt.0)then
        xcoremean=xcoremean/ncore
        ycoremean=ycoremean/ncore
        xcore2mean=xcore2mean/ncore
        ycore2mean=ycore2mean/ncore
        delxcore=xcore2mean-xcoremean**2
        if(delxcore.ge.0.)then
          delxcore=sqrt(delxcore)
        else
          write(ifmt,*)' delxcore**2<0',delxcore
          delxcore=0.
        endif
        delycore=(ycore2mean-ycoremean**2)
        if(delycore.ge.0.)then
          delycore=sqrt(delycore)
        else
          write(ifmt,*)' delycore**2<0',delycore
          delycore=0.
        endif
      endif  
      delcore=max(delxcore,delycore)
        !write(ifmt,'(a,i6,a,i6)')'Ncore=',ncore,'   Nescape=',nesc
        !print*,xcoremean,ycoremean,delxcore,delycore
        !print*,xmean,ymean
      xmean=xmean+xcoremean
      ymean=ymean+ycoremean
        !print*,xmean,ymean

      if(iorsdf.eq.5)then

        do n=1,nptla
          xptl(n)=xptl(n)-xcoremean
          yptl(n)=yptl(n)-ycoremean
          !     if( istptl(n).eq.5)
          !.    print*,'5555555', n, idptl(n) , istptl(n),'      '
          !.    , xptl(n),yptl(n), '   ',dezptl(n)
        enddo
        do kl=1,koll
          coord(1,kl)=coord(1,kl)-xcoremean
          coord(2,kl)=coord(2,kl)-ycoremean
        enddo
        do k=1,koll
          do n=1,nprmx(k)
            coordpr(1,n,k)=coordpr(1,n,k)-xcoremean
            coordpr(2,n,k)=coordpr(2,n,k)-ycoremean
          enddo
        enddo
        esu=0
        eco=0
        do i=1,nptla
        if(istptl(i).eq.0.or.istptl(i).eq.5)esu=esu+pptl(4,i)
        if(istptl(i).eq.5)eco=eco+pptl(4,i)
        enddo
        if(lprint)
     .  write(ifmt,'(a,5x,2i7,4x,f12.1,f9.1)')
     .  ' +++++Eall+++++(jintpo before exit)',ncore,nesc,esu,eco
        egycore=eco
        call iniCore(2)
        goto 1001
      endif

    !~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    !~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    !~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    !~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    !~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

      call  cluster(ntry,iret)
      if(iret.eq.8888)goto 8888
      if(iret.eq.1000)goto 1000
 
 !~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
 1001 continue
 !~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

c droplet decays moved to int (subroutine decaydroplets) 

 1000 continue
      call utprix('jintpo',ish,ishini,4)
      end

c----------------------------------------------------------------------
      subroutine cluster(ntry,iret)
c----------------------------------------------------------------------
c cut and paste from jinpo_drop version 3236
c cut and paste from jinto version 3127
c----------------------------------------------------------------------
c Cluster formation (second part) and cluster decay 
c----------------------------------------------------------------------
c Depends on scenario (in particular the decay and its flow settings) :
c iorsdf = 3 ... parametrized flow approach 
c          5 ... (will not reach this part)
c          6 ... cluster approach for initial conditions 
c iorsdf transmitted via aaa.h
c----------------------------------------------------------------------
      implicit none
#include "aaa.h"
      integer iret
      real xptl,yptl,zptl,tptl,optl,uptl,sptl,rptl
      common/cxyzt/xptl(mxptl+29),yptl(mxptl+29),zptl(mxptl+29)
     * ,tptl(mxptl+29),optl(mxptl+29),uptl(mxptl+29),sptl(mxptl+29)
     *,rptl(mxptl+29,3)
      integer mxcl,mxcli,nxjjj,nyjjj,nzjjj
      parameter(mxcl=4000,mxcli=1000)
      parameter(nxjjj=65,nyjjj=65,nzjjj=65)
      integer idropgrid(nxjjj,nyjjj,nzjjj)
     &       ,jdropgrid(nxjjj,nyjjj,nzjjj)
     &       ,jclu(nzjjj)
     &       ,mmji(mxcl,mxcli)
     &       ,jccl(mxcl,nflav,2),nseg(mxcl),mseg(mxcl),kclu(mxcl)
     &       ,naseg(0:mxcl),nfseg(mxcl),nsegmx(mxcl)
     &       ,jc(nflav,2),ke(6),jcjj(nflav,2)
     &       ,nsegmt(mxptl)   !,ic(2),nclk(nzjjj)
     &       ,nsegsuj
      common/jintpoc1/idropgrid,jdropgrid,jclu,kclu,nseg,nsegsuj
     &,jc
      real xcell,scell,zcell,delxce
      integer jjj,m1cell,m3cell,nptlb,nptla
      common/jintpoc2/xcell,scell,zcell,delxce
     &,jjj,m1cell,m3cell,nptlb,nptla
      double precision tpro,zpro,ttar,ztar,ttaus,detap,detat
     &                ,pptld(5,mxcl),p4tmp
      common   /cttaus/tpro,zpro,ttar,ztar,ttaus,detap,detat
      real eloss,xmean,ymean
      common/celoss/eloss
      common/cmean/xmean,ymean
      integer nfrx,kredonoco
      common/cnfrx/nfrx
      real delzet,delsce,vocell,ranphi,xlongcell
      common/cdelzet/delzet,delsce /cvocell/vocell,xlongcell
      common/cranphi/ranphi  /credonoco/kredonoco
      double precision ptest(5),ttest,p52,xmxms,bp,am2tmp
     &,am2tmpmx,amcmin,amcmi0,utdble
      parameter(xmxms=200d0)      !maximum mass for a subcluster
      real aa,bb,cc,dd,p4max,pjj52,r,rini,sg,taugm,tecm0,tm,xrot,xx,xy
     &,rangen,y,yrot,yy,utamnu,AngleEllipsoid
      integer i,ier,ii,iimx,iini,j,ji,jj,k,l
     &,mjjseg,mjjsegsum,mm,nclu,nh,ni,njj,np,nptl0,nptlbcf,nptlij,ns
     &,nst,ntmp,ntry,n,idum,mbinsum,ior,nmn,nmx
      real amt,rap

      integer getAccumJerr
c...prepare /cptl/ for clusters

      if(ish.ge.6)write(ifch,*)'prepare /cptl/ for clusters'
      do jj=1,jjj
         nptl=nptl+1
          istptl(nptl)=12
          do l=1,4
            pptl(l,nptl)=0.
            xorptl(l,nptl)=0
          enddo
          sptl(nptl)=0
          uptl(nptl)=0
          optl(nptl)=0
          desptl(nptl)=0
c limit the maximum number of subcluster to half the number of particle
c (not to have empty subclusters)
          nsegmx(jj)=max(1,min(mxcli,
     .                     nint(float(nseg(jj))/float(nsegsuj))))
      enddo

c...prepare /cptl/ for subclusters

      if(ish.ge.6)write(ifch,*)'prepare /cptl/ for subclusters'
      mm=0
      do jj=1,jjj
        do ii=1,nsegmx(jj)
          mm=mm+1
          if(mm.gt.mxcl)stop'cluster: mxcl too small.        '
          mmji(jj,ii)=mm
          mseg(mm)=0
          nptl=nptl+1
          istptl(nptl)=10
          do l=1,4
            pptld(l,mm)=0d0
            pptl(l,nptl)=0.
            xorptl(l,nptl)=0
          enddo
          sptl(nptl)=0
          uptl(nptl)=0
          optl(nptl)=0
          desptl(nptl)=0
          do l=1,nflav
            jccl(mm,l,1)=0
            jccl(mm,l,2)=0
          enddo
          iorptl(nptl)=nptla+jj
          jorptl(nptl)=0
          if(ii.eq.1)ifrptl(1,nptla+jj)=nptlb+mm
          ifrptl(2,nptla+jj)=nptlb+mm
        enddo
      enddo

c...separate string segments, add dense area segments to clusters

      if(ish.ge.6)write(ifch,*)'separate string segments'
      do 98 n=1,nptla
        if(istptl(n).ne.3)goto 98
        i=1+(xptl(n)+xcell)/delxce
        j=1+(yptl(n)+xcell)/delxce
        k=1+(sptl(n)+scell)/delsce
        if(    i.ge.1.and.i.le.m1cell
     &          .and.j.ge.1.and.j.le.m1cell
     &          .and.k.ge.1.and.k.le.m3cell)then
          jj=jdropgrid(i,j,k)
c share string segments into nsegmx clusters with low enough mass
          if(jj.gt.0)then
            iimx=nsegmx(jj)
c check if a cluster is empty
            do ii=1,iimx
              mm=mmji(jj,ii)
              if(mseg(mm).eq.0)goto 10          !not to have an empty cluster
            enddo
c select one randomly
            ii=1+rangen()*iimx
            ii=min(ii,iimx)
c save starting point
            iini=ii
            am2tmpmx=1e30
 9          ntmp=mmji(jj,ii)
            am2tmp=(pptld(4,ntmp)+pptl(4,n)+pptld(3,ntmp)+pptl(3,n))
     &            *(pptld(4,ntmp)+pptl(4,n)-pptld(3,ntmp)-pptl(3,n))
c check cluster future mass
            if(am2tmp.gt.0.05d0*xmxms*xmxms)then       !start looking for a better place when energy is larger than 20% of the maximum
c if larger than the limit
              if(am2tmp.lt.am2tmpmx)then
c and all are larger than the limit, then use later the cluster with minimum mass saved here
                mm=mmji(jj,ii)
                am2tmpmx=am2tmp
              endif
c go to next cluster and try again
              ii=ii+1
              if(ii.gt.iimx)ii=1
              if(ii.ne.iini)then
                goto 9
              else
c loop is complete, use the cluster with the minimum mass
                goto 10
              endif
            endif
c else continue with this cluster
            mm=mmji(jj,ii)
 10         continue
            mseg(mm)=mseg(mm)+1
            ifrptl(1,n)=mm      !local use of ifrptl
c           write(ifch,*)'mseg',mm,mseg(mm),n,istptl(n)
c     & ,pptl(4,n),pptl(3,n),idptl(n),i,j,k,sptl(n)
            p4tmp=0d0
            do l=1,3
             pptld(l,mm)=  pptld(l,mm)  + utdble(pptl(l,n))
             p4tmp=p4tmp+utdble(pptl(l,n))**2
            enddo
            p4tmp=sqrt(p4tmp+utdble(pptl(5,n))**2)
            pptld(4,mm)=  pptld(4,mm)  + p4tmp
c           if(mm.eq.86)write(ifch,*)'other',n
c     & ,pptl(4,n),pptl(3,n),pptl(5,n),idptl(n),k,p4tmp
c     & ,pptld(4,mm),pptld(3,mm),pptld(2,mm),pptld(1,mm)
c     & ,(pptld(4,mm)+pptld(3,mm))*(pptld(4,mm)-pptld(3,mm))
            ! particle position at starting point
            ior=iorptl(n)
            if(mod(abs(idptl(n)),100).eq.88)ior=iorptl(ior)
            do l=1,4
            xorptl(l,nptlb+mm)=xorptl(l,nptlb+mm)+xorptl(l,ior)
            enddo
            ! particle position at ttaus
            xptl(nptlb+mm)=xptl(nptlb+mm)+xptl(n)
            yptl(nptlb+mm)=yptl(nptlb+mm)+yptl(n)
            zptl(nptlb+mm)=zptl(nptlb+mm)+zptl(n)
            tptl(nptlb+mm)=tptl(nptlb+mm)+tptl(n)
            sptl(nptlb+mm)=sptl(nptlb+mm)+sptl(n)
            aa=cos(phievt)
            bb=sin(phievt)
            cc=-sin(phievt)
            dd=cos(phievt)
            xrot=xptl(n)*aa+yptl(n)*bb
            yrot=xptl(n)*cc+yptl(n)*dd
            uptl(nptlb+mm)=uptl(nptlb+mm)+xrot**2
            optl(nptlb+mm)=optl(nptlb+mm)+yrot**2
            desptl(nptlb+mm)=desptl(nptlb+mm)+xrot*yrot
            call idtr7( idptl(n) , jc )
c            id=idptl(n)
c            ida=iabs(id/10)
c            ids=id/iabs(id)
c            if(ida.ne.111.and.ida.ne.222.and.ida.ne.333)id=id/10*10
c            if(ida.eq.111.or. ida.eq.222.or. ida.eq.333)id=id/10*10+ids
c            if(ida.eq.213)id=1230*ids
c            ic(1)=idtrai(1,id,1)
c            ic(2)=idtrai(2,id,1)
c            call iddeco(ic,jc)
            do l=1,nflav
              jccl(mm,l,1)=jccl(mm,l,1)+jc(l,1)
              jccl(mm,l,2)=jccl(mm,l,2)+jc(l,2)
            enddo
          else
            idropgrid(i,j,k)=0
          endif
        endif
  98  continue

c...associate segments to clusters

      if(ish.ge.6)write(ifch,*)'associate segments to clusters'
      naseg(0)=0
      do jj=1,jjj
        do ii=1,nsegmx(jj)
          mm=mmji(jj,ii)
          naseg(mm)=naseg(mm-1)+mseg(mm)
          nfseg(mm)=0
        enddo
      enddo
      do n=1,nptla
        if(istptl(n).ne.3)goto 97
        istptl(n) = 7
        mm=ifrptl(1,n)
        nfseg(mm)=nfseg(mm)+1
        nsegmt(naseg(mm-1)+nfseg(mm))=n
 97   continue
      enddo
      do jj=1,jjj
        nst=0
        do ii=1,nsegmx(jj)
          mm=mmji(jj,ii)
          if(mseg(mm).ne.nfseg(mm))stop'cluster: mseg.ne.nfseg        '
          nst=nst+mseg(mm)
        enddo
        if(nst.ne.nseg(jj))stop'sum(mseg(mm)).ne.nseg(jj)'
      enddo

c...finish cluster storage to /cptl/

      if(ish.ge.6)write(ifch,*)'finish cluster storage to /cptl/'
      xx=0.
      yy=0.
      xy=0.
      mjjsegsum=0
      do jj=1,jjj
       njj=nptla+jj
       mjjseg=0
       do l=1,nflav
         jcjj(l,1)=0
         jcjj(l,2)=0
       enddo
       do ii=1,nsegmx(jj)
        mm=mmji(jj,ii)
        n=nptlb+mm

        do l=1,nflav
          jc(l,1)=jccl(mm,l,1)
          jc(l,2)=jccl(mm,l,2)
          ke(l)=jc(l,1)-jc(l,2)
          jcjj(l,1)=jcjj(l,1)+jc(l,1)
          jcjj(l,2)=jcjj(l,2)+jc(l,2)
        enddo
        call idenct(jc,idptl(n)
     *  ,ibptl(1,n),ibptl(2,n),ibptl(3,n),ibptl(4,n))
        ttest=0d0
        do ji=1,4
         ptest(ji)=0d0
          do ns=naseg(mm-1)+1,naseg(mm)
            ni=nsegmt(ns)
            ptest(ji)=ptest(ji)+utdble(pptl(ji,ni))
          enddo
         ptest(ji)=abs(ptest(ji)-pptld(ji,mm))
         ttest=ttest+ptest(ji)
        enddo
        amcmin=utdble(utamnu(ke(1),ke(2),ke(3),ke(4),ke(5),ke(6),4))
        p52=(pptld(4,mm)+pptld(3,mm))*(pptld(4,mm)-pptld(3,mm))
     &      -pptld(1,mm)**2-pptld(2,mm)**2
        pptld(5,mm)=sqrt(max(0d0,p52))
        if(pptld(5,mm).lt.amcmin)then
        amcmi0=utdble(1.1*utamnu(ke(1),ke(2),ke(3),ke(4),ke(5),ke(6),5))
          if(idptl(n).ne.800000000)
     &       call setAccumJerr(2,getAccumJerr(2)+1)  !tp20140605 error here
          if(pptld(4,mm).gt.amcmi0)then
c give enough mass to the cluster by rescaling momentum
            if(pptld(4,mm).ge.amcmin)then
              pptld(5,mm)=amcmin
            else
              pptld(5,mm)=amcmi0
            endif
            bp=sqrt((pptld(4,mm)+pptld(5,mm))*(pptld(4,mm)-pptld(5,mm))
     &            /(pptld(3,mm)*pptld(3,mm)+pptld(2,mm)*pptld(2,mm)
     &                                      +pptld(1,mm)*pptld(1,mm)))
            pptld(1,mm)=bp*pptld(1,mm)
            pptld(2,mm)=bp*pptld(2,mm)
            pptld(3,mm)=bp*pptld(3,mm)
c           write(ifch,*)"ici ",n,sqrt(p52),pptld(5,mm),bp,pptld(4,mm),
c     &  sqrt(pptld(3,mm)*pptld(3,mm)
c     &                   +pptld(2,mm)*pptld(2,mm)
c     &                   +pptld(1,mm)*pptld(1,mm)
c     &                   +pptld(5,mm)*pptld(5,mm))
c      write(ifch,*)'droplet uds=',ke(1),ke(2),ke(3),'   E=',pptld(5,mm)
          elseif(idptl(n).eq.800000000)then   !pure energy loss : give it back
            do ns=naseg(mm-1)+1,naseg(mm)
              ni=nsegmt(ns)
              if(mod(abs(idptl(ni)),100).eq.88)then !restore lost energy
                call fpartdel(ni,idum)
              else
                print *,'cluster',ni,idptl(ni)
                call utstop('Should not happen in ind.f !&')
              endif
            enddo
          else
            call setAccumJerr(3,getAccumJerr(3)+1)
c            print *,pptld(5,mm),pptld(4,mm),sqrt(pptld(3,mm)*pptld(3,mm)
c     &                   +pptld(2,mm)*pptld(2,mm)
c     &                   +pptld(1,mm)*pptld(1,mm)
c     &                   +amcmi0**2),amcmi0,jc,idptl(n)
            pptld(5,mm)=amcmi0
            pptld(4,mm)=sqrt(pptld(3,mm)*pptld(3,mm)
     &                   +pptld(2,mm)*pptld(2,mm)
     &                   +pptld(1,mm)*pptld(1,mm)
     &                   +pptld(5,mm)*pptld(5,mm))
          endif
        endif
        if(ish.ge.2.and.(abs(ttest).gt.1.d0
     &    .or.pptld(5,mm).gt.xmxms))then
          call utmsg('cluster&')
          write(ifmt,*)'***** Warning in cluster !',ntry
          write(ifch,*)'***** cluster: momenta messed up (ttest > 0)'
          write(ifch,*)'*****',mm,n,mseg(mm),p52,ttest
          write(ifch,*)'*****',jj,nsegmx(jj),amcmin
          write(ifch,'(a,16x,5f15.4)')' *****',(pptld(ji,mm),ji=1,5)
          do ns=naseg(mm-1)+1,naseg(mm)
            ni=nsegmt(ns)
            write(ifch,'(a,i7,i9,5f15.4,f12.4)')' *****',ni,idptl(ni)
     *     ,(pptl(ji,ni),ji=1,4),pptl(5,ni)**2
     *     ,(pptl(4,ni)+pptl(3,ni))*(pptl(4,ni)-pptl(3,ni))
     *       -pptl(1,ni)**2-pptl(2,ni)**2
          enddo
        endif
        if(pptld(5,mm).gt.xmxms)then
c          nnn=naseg(mm)-naseg(mm-1)-1
c          if(nnn.lt.10)then     !if few particles only, skip all (some have very large momenta)
c            nh=nnn
c            nmn=naseg(mm-1)+1
c            nmx=naseg(mm)
c          else
            p4max=0.
            nh=0
            do ns=naseg(mm-1)+1,naseg(mm)
              ni=nsegmt(ns)
              if(pptl(4,ni).ge.p4max)then
                nh=ns
                p4max=pptl(4,ni)
              endif
            enddo
            nmn=nh
            nmx=nh
c          endif
          if(nh.le.0)then
            stop'Cannot be in cluster ...'
          else   !put back nh as normal particle
c           do not use this particle any more
            do ns=nmn,nmx
              nh=nsegmt(ns)
              iaaptl(nh)=0
              if(mod(abs(idptl(nh)),100).eq.88)then !restore lost energy
                call fpartdel(nh,idum)
              elseif(idptl(nh).lt.1e4)then
                istptl(nh) = 0
                ifrptl(1,nh) = 0
                ifrptl(2,nh) = 0
              else
                istptl(nh) = 10
              endif
            enddo
c but don't change the grid configuration (define once at the beginning)
          endif
         if(ish.ge.1)
     &    write(ifch,*)'***** Redo cluster without heavy particle :'
     &                 ,nsegmt(nmn),nsegmt(nmx),ntry
          goto 8888
        endif
        do l=1,5
         pptl(l,n)=sngl(pptld(l,mm))
        enddo
        mjjseg=mjjseg+mseg(mm)
        do l=1,4
          pptl(l,njj)=pptl(l,njj)+pptl(l,n)
          xorptl(l,njj)=xorptl(l,njj)+xorptl(l,n)
          xorptl(l,n)=xorptl(l,n)/float(mseg(mm))
        enddo
        !kw18  
        amt=pptl(5,n)**2+pptl(1,n)**2+pptl(2,n)**2
        if(amt.gt.0..and.pptl(4,n).gt.0.)then
          amt=sqrt(amt)
          rap=sign(1.,pptl(3,n))*alog((pptl(4,n)+abs(pptl(3,n)))/amt)
        else
          rap=100.
        endif
        xorptl(3,n)=tauzer*sinh(rap)
        xorptl(4,n)=tauzer*cosh(rap)
        !kw18 
        xptl(njj)=xptl(njj)+xptl(n)
        yptl(njj)=yptl(njj)+yptl(n)
        zptl(njj)=zptl(njj)+zptl(n)
        tptl(njj)=tptl(njj)+tptl(n)
        sptl(njj)=sptl(njj)+sptl(n)
        uptl(njj)=uptl(njj)+uptl(n)
        optl(njj)=optl(njj)+optl(n)
        desptl(njj)=desptl(njj)+desptl(n)
        xptl(n)=xptl(n)/float(mseg(mm))
        yptl(n)=yptl(n)/float(mseg(mm))
        zptl(n)=zptl(n)/float(mseg(mm))
        tptl(n)=tptl(n)/float(mseg(mm))
        sptl(n)=sptl(n)/float(mseg(mm))
        uptl(n)=uptl(n)/float(mseg(mm))
        optl(n)=optl(n)/float(mseg(mm))
        desptl(n)=desptl(n)/float(mseg(mm))
        istptl(n)=10
        ifrptl(1,n)=0
        ifrptl(2,n)=0
        tivptl(1,n)=xorptl(4,n)
        tivptl(2,n)=ainfin
        ityptl(n)=60
        radptl(n)=0
        dezptl(n)=0.
        zpaptl(1,n)=0.
        zpaptl(2,n)=0.
       enddo !ii
       do l=1,4
        xorptl(l,njj)=xorptl(l,njj)/float(mjjseg)
       enddo
       !kw18  
       amt=pptl(5,njj)**2+pptl(1,njj)**2+pptl(2,njj)**2
       if(amt.gt.0..and.pptl(4,njj).gt.0.)then
         amt=sqrt(amt)
         rap=sign(1.,pptl(3,njj))
     .    *alog((pptl(4,njj)+abs(pptl(3,njj)))/amt)
       else
         rap=100.
       endif
       xorptl(3,njj)=tauzer*sinh(rap)
       xorptl(4,njj)=tauzer*cosh(rap)
       !kw18 
       mjjsegsum=mjjsegsum+mjjseg
       xx=xx+uptl(njj)
       yy=yy+optl(njj)
       xy=xy+desptl(njj)
c       radptl(njj)=max(0.0001,sqrt(5./3.*(uptl(njj)+optl(njj))))
       xptl(njj)=xptl(njj)/float(mjjseg)
       yptl(njj)=yptl(njj)/float(mjjseg)
       zptl(njj)=zptl(njj)/float(mjjseg)
       tptl(njj)=tptl(njj)/float(mjjseg)
       sptl(njj)=sptl(njj)/float(mjjseg)
       uptl(njj)=uptl(njj)/float(mjjseg)
       optl(njj)=optl(njj)/float(mjjseg)
       desptl(njj)=desptl(njj)/float(mjjseg)
c       radptl(njj)=max(0.0001,sqrt(5./3.*(uptl(njj)+optl(njj))))
       pjj52=(pptl(4,njj)+pptl(3,njj))*(pptl(4,njj)-pptl(3,njj))
     &    -pptl(1,njj)**2-pptl(2,njj)**2
       pptl(5,njj)=0
       if(pjj52.gt.0)pptl(5,njj)=sqrt(pjj52)
       ityptl(njj)=60
       call idenct(jc,idptl(njj)
     *  ,ibptl(1,njj),ibptl(2,njj),ibptl(3,njj),ibptl(4,njj))
       zpaptl(1,njj)=0.
       zpaptl(2,njj)=0.
      enddo !jj

c store the number of bin of each cluster in radptl
      mbinsum=0
      do k=1,m3cell
        do j=1,m1cell
          do i=1,m1cell
            if(jdropgrid(i,j,k).gt.0)then
              njj=nptla+jdropgrid(i,j,k)
              mbinsum=mbinsum+1
              dezptl(njj)=dezptl(njj)+1.
            endif
          enddo
        enddo
      enddo

c...ranphi

      ranphi=0
      rini=0.
      if(mjjsegsum.gt.0)then
c        rini=max(0.01,sqrt(5./3.*(xx+yy))) !<r**2>=3/5*R**2 for sphere of radius R
        xx=xx/float(mjjsegsum)
        yy=yy/float(mjjsegsum)
        xy=xy/float(mjjsegsum)
        ranphi=AngleEllipsoid(xx,yy,xy)
        rini=max(0.01,sqrt(5./3.*(xx+yy))) !<r**2>=3/5*R**2 for sphere of radius R
      endif

c...print

      if(ish.ge.5)then
        write(ifch,*)'print'
        do k=1,m3cell
        write(ifch,*)'k=',k,'  jclu=',jclu(k)
        do j=m1cell,1,-1
        write(ifch,'(10i4,3x,10i4)')(idropgrid(i,j,k),i=1,m1cell)
     &                ,(jdropgrid(i,j,k),i=1,m1cell)
        enddo
        enddo
        write(ifch,'(a,a)')
     &    '    k   jj  nseg      mm  mseg     n      mass'
     &   ,'         s         y         z         t '
        do jj=1,jjj
         do ii=1,max(1,nint(1.*nseg(jj)/nsegsuj))
           mm=mmji(jj,ii)
           n=nptlb+mm
           sg=pptl(3,n)/abs(pptl(3,n))
           tm=sqrt(pptl(5,n)**2+pptl(1,n)**2+pptl(2,n)**2)
           y=sg*alog((pptl(4,n)+sg*pptl(3,n))/tm)
c           if(kclu(jj).eq.44)print *,tm,pptl(4,n),pptl(3,n),iorptl(n)
           write(ifch,'(2i5,i6,i8,2i6,5f10.3)')
     &       kclu(jj),jj,nseg(jj),mm,mseg(mm),n,pptl(5,n)
     &      ,sptl(n),y,xorptl(3,n),xorptl(4,n)
         enddo
        enddo
      endif

c...decay
      tecm0=0.
      iret=0
      if(jjj.gt.0)then     !decay only if some cluster produced
      if(ish.ge.5)write(ifch,*)'decay ...'
      if(ifrade.eq.0.or.ispherio.gt.0)goto 1000
      if(jdecay.eq.0)goto 1000
      nptlbcf=nptl
      nptl0=nptl
      ier=0
      hydt='---'
      if(hydt.eq.'---'.or.ier.ne.0)then  !ier>0 possible for very peripheral
        nclu=0
        ptest(1)=0d0
        ptest(2)=0d0
        ptest(3)=0d0
        ptest(4)=0d0
        ptest(5)=0d0
        do jj=1,jjj
c          radptl(nptla+jj)=float(jjj)
c          radptl(nptla+jj)=float(mbinsum)
          do ii=1,max(1,nint(1.*nseg(jj)/nsegsuj))
           mm=mmji(jj,ii)
           np=nptlb+mm
           nptlij=nptl
          !---------------------------
           call hnbaaa(iorsdf,np,iret)
          !---------------------------
           if(iret.ne.0)then           !decay failed, restore particles
             istptl(np)=istptl(np)+2
             do ns=naseg(mm-1)+1,naseg(mm)
               n=nsegmt(ns)
               if(mod(abs(idptl(n)),100).eq.88)then   !restore lost energy
                 call fpartdel(n,idum)
               elseif(idptl(n).lt.1e4)then
                 istptl(n) = 0
                 ifrptl(1,n) = 0
                 ifrptl(2,n) = 0
               else
                 istptl(n) = 10
               endif
             enddo
           else !decay ok
             !print*,'cluster decay',np,nptlij+1,nptl
             ifrptl(1,np)=nptlij+1
             ifrptl(2,np)=nptl
             do i=nptlij+1,nptl
               iorptl(i)=np
               jorptl(i)=0
             enddo 
             if(iorsdf.eq.3)then
               tecm0=tecm0+pptl(4,np)
               do i=nptl0+1,nptl
                 if(ityptl(i).eq.60)then
                   nclu=nclu+1
                   ptest(1)=ptest(1)+utdble(pptl(1,i))
                   ptest(2)=ptest(2)+utdble(pptl(2,i))
                   ptest(3)=ptest(3)+utdble(pptl(3,i))
                   ptest(4)=ptest(4)+utdble(pptl(4,i))
                 endif
               enddo
               nptl0=nptl
             endif
             !~~~~~~~~~~~~~~~~~~~~~~
             do ns=naseg(mm-1)+1,naseg(mm)
               n=nsegmt(ns) !particle index n of particle melted into cluster with index np
               !write(ifch,*)'cluster parents',n,np
               ifrptl(2,n)=-np
             enddo
             !~~~~~~~~~~~~~~~~~~~~~~  
           endif
          enddo
        enddo
      endif
      do jj=1,jjj
       do ii=1,max(1,nint(1.*nseg(jj)/nsegsuj))
        mm=mmji(jj,ii)
        np=nptlb+mm
        istptl(np)=istptl(np)+1
       enddo
      enddo

c
c a commented piece following  ### flowpp ### has removed
c       to be found (if needed again) in version 3236 
c

      do n=nptlbcf+1,nptl
       if(mod(abs(idptl(n)),100).ne.88)then   !should not happen ...
        istptl(n)=0
        ifrptl(1,n)=0
        ifrptl(2,n)=0
        tivptl(1,n)=xorptl(4,n)
        call idtau(idptl(n),pptl(4,n),pptl(5,n),taugm)
        r=rangen()
        tivptl(2,n)=tivptl(1,n)+taugm*(-alog(r))
        radptl(n)=0.
        dezptl(n)=0.
        itsptl(n)=0
        !variable needed for other purpose!rinptl(n)=kclu(jj)-m3cell/2
       else
         write(ifmt,*)'Really ???'
       endif
      enddo

      endif
      
      iret=0
      return

 8888 iret=8888
      return

 1000 iret=1000 !skip also the following droplet decay
      return

      end


c----------------------------------------------------------------------
      subroutine fpartloss(i,j,k,n,dp1,dp2,dp3,dp4)
c----------------------------------------------------------------------
c cut and paste from jinto_drop version 3236
c cut and paste from jinto version 3200
c----------------------------------------------------------------------
c Creates virtual core particles (with ist=3) carrying the lost energy 
c    of an escaped segment
c Counter idropgrid(i,j,k) increase by 1  (counting this particle)
c----------------------------------------------------------------------
      implicit none
#include "aaa.h"
      real dp1,dp2,dp3,dp4
      integer i,j,k,n,nn
      real xptl,yptl,zptl,tptl,optl,uptl,sptl,rptl
      common/cxyzt/xptl(mxptl+29),yptl(mxptl+29),zptl(mxptl+29)
     * ,tptl(mxptl+29),optl(mxptl+29),uptl(mxptl+29),sptl(mxptl+29)
     *,rptl(mxptl+29,3)
      integer mxcl,mxcli,nxjjj,nyjjj,nzjjj
      parameter(mxcl=4000,mxcli=1000)
      parameter(nxjjj=65,nyjjj=65,nzjjj=65)
      integer idropgrid(nxjjj,nyjjj,nzjjj)
     &       ,jdropgrid(nxjjj,nyjjj,nzjjj)
     &       ,jclu(nzjjj)
     &       ,nseg(mxcl),kclu(mxcl)
     &       ,jc(nflav,2),nsegsuj
      common/jintpoc1/idropgrid,jdropgrid,jclu,kclu,nseg,nsegsuj
     &,jc
      real xcell,scell,zcell,delxce
      integer jjj,m1cell,m3cell,nptlb,nptla
      common/jintpoc2/xcell,scell,zcell,delxce
     &,jjj,m1cell,m3cell,nptlb,nptla
      integer kk,l
      real taugm, rangen
c create fake particle to count energy lost in cluster
      idropgrid(i,j,k)=idropgrid(i,j,k)+1 !keep part of the segment in grid for cluster formation
      nptl=nptl+1
      nptla=nptla+1
      istptl(nptl)=3            !core particle
      iaaptl(nptl)=1
      pptl(1,nptl)=dp1
      pptl(2,nptl)=dp2
      pptl(3,nptl)=dp3
      pptl(4,nptl)=dp4
      do l=1,3
        xorptl(l,nptl)=xorptl(l,n)
      enddo
! mother
      tivptl(1,n)=xorptl(4,n)
      call idtau(idptl(n),pptl(4,n),pptl(5,n),taugm)
      tivptl(2,n)=tivptl(1,n)+taugm*(-alog(rangen()))
      ifrptl(1,n)=nptl
! daughter
      pptl(5,nptl)=(pptl(4,nptl)-pptl(3,nptl))
     &            *(pptl(4,nptl)+pptl(3,nptl))
     &             -pptl(1,nptl)**2-pptl(2,nptl)**2
      if(pptl(5,nptl).gt.0.)then
        pptl(5,nptl)=sqrt(pptl(5,nptl))
      else
        pptl(5,nptl)=1e-30
      endif
      nn=nptl
      call updateXor(n,nn)
      iorptl(nptl)=n
      jorptl(nptl)=0
      idptl(nptl)=idptl(n)*100+sign(88,idptl(n))
      ityptl(nptl)=ityptl(n)+2
      xptl(nptl)=xptl(n)
      yptl(nptl)=yptl(n)
      zptl(nptl)=zptl(n)
      sptl(nptl)=sptl(n)
      dezptl(nptl)=dezptl(n)
      zpaptl(1,nptl)=zpaptl(1,n)
      zpaptl(2,nptl)=zpaptl(2,n)
      if(ish.ge.6)write(ifch,*)'--> Virtual particle : ',nptl
     & ,idptl(nptl),iorptl(nptl),ityptl(nptl),(pptl(kk,nptl),kk=1,5)

      return
      end

c----------------------------------------------------------------------
      subroutine fpartdel(n,ncpt)
c----------------------------------------------------------------------
c cut and paste from jinto_drop version 3236
c cut and paste from jinto version 3200
c----------------------------------------------------------------------
c Deletes virtual core particles, energy back to "parent" particle 
c----------------------------------------------------------------------
      implicit none
#include "aaa.h"
      integer n,ncpt
      integer k,ior
      real taugm, rangen

      if(mod(abs(idptl(n)),100).eq.88)then !virtual particle is lost (no cluster here)
        ncpt=ncpt-1             !to get proper counting
        iaaptl(n) = 0           !don't use this particle later
        istptl(n) = 2   !invalid
        ior=iorptl(n)
        ifrptl(1,ior) = 0
        ifrptl(2,ior) = 0
        do k=1,4                !restore energy to mother particle
          pptl(k,ior)=pptl(k,ior)+pptl(k,n)
        enddo
c        pptl(4,ior)=sqrt(pptl(1,ior)**2+pptl(2,ior)**2
c     &       +pptl(3,ior)**2+pptl(5,ior)**2)
        call idtau(idptl(ior),pptl(4,ior),pptl(5,ior),taugm)
        tivptl(2,ior)=tivptl(1,ior)+taugm*(-alog(rangen()))
      endif
      return
      end

c----------------------------------------------------------------------
      subroutine SetPFE() !Parametrized Fluid Expansion
c----------------------------------------------------------------------
c cut and paste from jinto version 3200
c----------------------------------------------------------------------
c Needed for the parametrized flow approach
c (In  SetEffectiveFlow several k planes are mixed, which should not be 
c in the general framework)
c----------------------------------------------------------------------
c set the effective flow such that the different behavior observed in
c pp, pA and AA is reproduced :
c pp : strong radial flow and decrease of high multiplicity (long flow)
c pA : same radial flow than in pp but less multiplicity reduction
c AA : slow increase of radial flow and decrease of multiplicity
c This can be explained by the density of particles. Strong radial flow
c coming from very hot and small point due to MPI in a pair and less 
c reduction of multiplicity if the system is extended.
c----------------------------------------------------------------------
      implicit none
#include "aaa.h"
      real bglaub,f1,f2
      common/cbglaub/bglaub
      integer ncol,kolpt,ncoli
      common/col3/ncol,kolpt,ncoli
      real ycor,fradflo,yco,yrmax
      real xptl,yptl,zptl,tptl,optl,uptl,sptl,rptl
      common/cxyzt/xptl(mxptl+29),yptl(mxptl+29),zptl(mxptl+29)
     * ,tptl(mxptl+29),optl(mxptl+29),uptl(mxptl+29),sptl(mxptl+29)
     *,rptl(mxptl+29,3)
      integer mxcl,mxcli,nxjjj,nyjjj,nzjjj
      parameter(mxcl=4000,mxcli=1000)
      parameter(nxjjj=65,nyjjj=65,nzjjj=65)
      integer idropgrid(nxjjj,nyjjj,nzjjj)
     &       ,jdropgrid(nxjjj,nyjjj,nzjjj)
     &       ,jclu(nzjjj),nsegp4(mxcl)
     &       ,nseg(mxcl),kclu(mxcl)
     &       ,jc(nflav,2),nsegsuj
      common/jintpoc1/idropgrid,jdropgrid,jclu,kclu,nseg,nsegsuj
     &,jc
      double precision p4mean(5,mxcl),xymean(2,mxcl),amamax,amamin
     &,amacur,amc2,utdble,ectot,amctot
     &,xyomean(2,mxcl)!,xwcrmean,ywcrmean,xwcr2mean,ywcr2mean,amweight,wpow
      real xcell,scell,zcell,delxce
      integer jjj,m1cell,m3cell,nptlb,nptla
      common/jintpoc2/xcell,scell,zcell,delxce
     &,jjj,m1cell,m3cell,nptlb,nptla
      real delzet,delsce,vocell,xlongcell
      common/cdelzet/delzet,delsce /cvocell/vocell,xlongcell
      integer jj,nsegtot,i,j,k,n,kk,idum,jjs,jjf,jjc,jj0,k1,k2,iofuse
     &,ncelltot,ncellrad,ncellong,ncellmax,isgn,jjl,ior
      real rangen,voll,fgeo!,yyrmax
      logical lnma

      if(ish.ge.6)write(ifch,*)'count cells'
      jjs=0
      ncelltot=0
      ncellrad=0
      ncellong=0
      do k=1,m3cell
        ncellmax=0
        do j=1,m1cell
          do i=1,m1cell
            if(jdropgrid(i,j,k).gt.0)then
              ncelltot=ncelltot+1 ! count number of active cells
              ncellmax=ncellmax+1
            endif
          enddo
        enddo
        if(ncellmax.gt.0)ncellong=ncellong+1
        ncellrad=max(ncellrad,ncellmax)
        do jj=1,jclu(k)
          jjs=jjs+1
          kclu(jjs)=k
          nseg(jjs)=0
          nsegp4(jjs)=0
          xymean(1,jjs)=0d0
          xymean(2,jjs)=0d0
          xyomean(1,jjs)=0d0
          xyomean(2,jjs)=0d0
          p4mean(1,jjs)=0d0
          p4mean(2,jjs)=0d0
          p4mean(3,jjs)=0d0
          p4mean(4,jjs)=0d0
          p4mean(5,jjs)=0d0
         enddo
      enddo

c...calculate mean energy of segments going into in clusters

      if(ish.ge.6)write(ifch,*)
     &'calculate mean energy of segments going into in clusters'
      do 95 n=1,nptla
        if(iaaptl(n).eq.0)goto 95
        i=1+(xptl(n)+xcell)/delxce
        j=1+(yptl(n)+xcell)/delxce
        k=1+(sptl(n)+scell)/delsce
        if(    i.ge.1.and.i.le.m1cell
     &          .and.j.ge.1.and.j.le.m1cell
     &          .and.k.ge.1.and.k.le.m3cell)then
          jj=jdropgrid(i,j,k)
          if(jj.gt.0)then
            nsegp4(jj)=nsegp4(jj)+1
            do kk=1,4
              p4mean(kk,jj)=p4mean(kk,jj)+utdble(pptl(kk,n))
            enddo
            xymean(1,jj)=xymean(1,jj)+xptl(n)
            xymean(2,jj)=xymean(2,jj)+yptl(n)
            ior=iorptl(n)
            if(mod(abs(idptl(n)),100).eq.88)ior=iorptl(ior)
            xyomean(1,jj)=xyomean(1,jj)+xorptl(1,ior)
            xyomean(2,jj)=xyomean(2,jj)+xorptl(2,ior)
c           if(jj.eq.9)write(ifch,*)'after',n
c     & ,pptl(4,n),pptl(3,n),idptl(n),jj,nseg(jj),i,j,k
          else
            call fpartdel(n,idum)
         endif
        endif
 95   continue

c define flow according to total mass and volume of clusters

      ectot=0.d0
      amctot=0.d0
      nsegtot=0
      do jj=1,jjj
c check if cluster in not empty
        if(p4mean(4,jj).gt.0.d0)then
          nsegtot=nsegtot+nsegp4(jj)
          ectot=ectot+p4mean(4,jj)
          amc2=(p4mean(4,jj)+p4mean(3,jj))
     &     *(p4mean(4,jj)-p4mean(3,jj))-p4mean(1,jj)**2-p4mean(2,jj)**2
c          print *,p4mean(5,jj),sqrt(max(0.d0,amc2))
          if(amc2.gt.0d0)then
            p4mean(5,jj)=sqrt(amc2)
            amctot=amctot+p4mean(5,jj)
          endif
        else
        endif
      enddo

      if(iorsdf.eq.3)then
        call pfe3paramget2( f1 , f2 ) !yrmax,yco (we dont consider ycoj here)
      elseif (iorsdf.eq.6)then
        call pfe6paramget2( f1 , f2 ) !yrmax,yco
      else
        stop'ERROR 23082022'
      endif 
      yrmax=dble(f1)
      yco=f2
      if(yrmax.gt.1d-5)then
        fradflo=sngl(1d0/
     &  ((sinh(yrmax)*yrmax-cosh(yrmax)+1d0)/(yrmax**2/2d0)))
      else
        fradflo=1.
        yrmax=0d0
      endif
      if(yco.gt.1.e-2)then
        ycor=sinh(yco)/yco/fradflo
      else
        ycor=1.
      endif


      if(ish.ge.3)write(ifch,*)'yrmax,ycor,voll='
     .                         ,yrmax,ycor,voll,fgeo

      iofuse=1   !1=fuse in on k-plane only, 2=fuse as long as the mass is too low (no good results)
      jjs=jjj
      if(iofuse.gt.0.and.jjs.gt.0.and.yrmax.lt.ainfin)then
        amamin=utdble(ycor)
        amamax=min(50d0,utdble(ycor)*amctot/dble(jjs))
        amacur=0d0
        if(rangen().gt.0.5.or.iofuse.eq.1)then
          isgn=1
          jj=1
          jjl=jjs
        else
          isgn=-1
          jj=jjs
          jjl=1
        endif
        jjj=1
        jj0=jj
        jjf=0
        do while(isgn*jj.le.isgn*jjl)
          jjc=jjj
          k=kclu(jj0)
          amacur=amacur+p4mean(5,jj)
c          print *,k,jj,jjl,jclu(k),amacur,jjf,jjc
          if(isgn*jj.lt.isgn*jjl)then
            kk=kclu(jj+isgn)
            if(iofuse.eq.2)then
              lnma=amacur.lt.amamin
     .        .or.(amacur+p4mean(5,jj+isgn)).lt.amamax
              if(jj.eq.jjl-isgn)lnma=amacur.lt.amamin
     .                           .or.p4mean(5,jj+isgn).lt.amamax
            else
              lnma=.false.
            endif
          else
            if(iofuse.eq.2)then
              kk=k
              lnma=amacur.lt.amamin
            else
              kk=0
              lnma=.false.
            endif
          endif
          if(isgn*jj.lt.isgn*jjl.and.(lnma
     .                 .or.(amacur.lt.amamax.and.kk.eq.k)))then
c if mass too low and future mass not too large, fuse cluster jjf into jjc
            if(isgn*jjf.lt.isgn*jj0)jj0=jj
            jjf=jj+isgn
          else
c otherwise save cluster jjc and create a new one
            jjf=jj0
            amacur=0d0
            k2=kk
            kk=kclu(jj0)
            k1=kk
            do j=jj0,jj
              kclu(j)=0         !reset kclu
            enddo
            kclu(jjj)=nint(float(k2+k1)/2.) !define new kclu as the average one 
            jjj=jjj+1
            jj0=jj+isgn
          endif
          do j=1,m1cell
            do i=1,m1cell
              if(jdropgrid(i,j,kk).eq.jjf)then
                jdropgrid(i,j,kk)=jjc
              endif
            enddo
          enddo
          jj=jj+isgn
        enddo
        jjj=jjc       !final number of cluster
      endif


c redefine cluster masses after fusion to match the final list of clusters
      if(amctot.gt.0.d0)then

        if(jjj.ne.jjs)then
          jjc=0
          do k=1,m3cell
            do jj=1,jclu(k)
              jjc=jjc+1
              p4mean(1,jjc)=0d0
              p4mean(2,jjc)=0d0
              p4mean(3,jjc)=0d0
              p4mean(4,jjc)=0d0
              p4mean(5,jjc)=0d0
            enddo
          enddo

          do 195 n=1,nptla
            if(iaaptl(n).eq.0)goto 195
            i=1+(xptl(n)+xcell)/delxce
            j=1+(yptl(n)+xcell)/delxce
            k=1+(sptl(n)+scell)/delsce
            if(    i.ge.1.and.i.le.m1cell
     &           .and.j.ge.1.and.j.le.m1cell
     &           .and.k.ge.1.and.k.le.m3cell)then
              jj=jdropgrid(i,j,k)
              if(jj.gt.0)then
                do kk=1,4
                  p4mean(kk,jj)=p4mean(kk,jj)+utdble(pptl(kk,n))
                enddo
              endif
            endif
 195      continue
        endif
        do jj=1,jjj
          if(p4mean(4,jj).gt.0.d0.and.jjs.ne.jjj)then
            amc2=(p4mean(4,jj)+p4mean(3,jj))
     &     *(p4mean(4,jj)-p4mean(3,jj))-p4mean(1,jj)**2-p4mean(2,jj)**2
            if(amc2.gt.0d0)then
              p4mean(5,jj)=sqrt(amc2)
            else
              p4mean(5,jj)=0d0
            endif
          elseif(jjs.ne.jjj)then
            p4mean(5,jj)=0d0
          endif
          radptl(nptla+jj)=p4mean(5,jj)/amctot
        enddo
      endif


      return
      end

