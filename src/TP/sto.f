C
C  This file is part of EPOS4
C  Copyright (C) 2022 research institutions and authors (See CREDITS file)
C  This file is distributed under the terms of the GNU General Public License version 3 or later
C  (See COPYING file for the text of the licence)
C

c-----------------------------------------------------------------------
      subroutine estore
c-----------------------------------------------------------------------
c     writes the results of a simulation into the file with unit ifdt
c     contains a description of the stored variables.
c     modifiable by the user
c-----------------------------------------------------------------------
#include "aaa.h"

c  count the number of particles to be stored (--> nptevt)

      nptevt=0
      do i=1,nptl
        if(istptl(i).le.istmax)nptevt=nptevt+1
      enddo


      write(ifdt,*)nptevt,bimevt,phievt,kolevt,pmxevt,egyevt
     *           ,npjevt,ntgevt,qsqevt,typevt
      do n=1,nptl

       if(istptl(n).le.istmax)then !store events with istptl < istmax
        
        write(ifdt,*)n,idptl(n),pptl(1,n),pptl(2,n),pptl(3,n),pptl(4,n)
     *            ,pptl(5,n),iorptl(n),jorptl(n),istptl(n),ityptl(n)
     *            ,xorptl(1,n),xorptl(2,n),xorptl(3,n),xorptl(4,n)
     *            ,ifrptl(1,n),ifrptl(2,n),dezptl(n)

       endif

      enddo
      write(ifdt,'(A)') ' '         !to finish the file

      return
      end

c-----------------------------------------------------------------------
      subroutine hepmcstore(iout)
c-----------------------------------------------------------------------
c     writes the results of a simulation into the file with unit ifdt
c     contains a description of the stored variables.
c     iout < 0 : no nuclei as beam (not recognize by ROOT ???)
c-----------------------------------------------------------------------
#include "aaa.h"
      double precision phep2,vhep2
      integer       nhep2,isthep2,idhep2,jmohep2,jdahep2

      dimension isthep2(nmxhep),idhep2(nmxhep),jmohep2(2,nmxhep)
     &,jdahep2(2,nmxhep),phep2(5,nmxhep),vhep2(4,nmxhep)
      integer iepo2hep(mxptl),ihep2epo(nmxhep),ihep2epo2(nmxhep)

      double precision pprojin,ptargin

      logical lrcore,lrcor0,lclean

c  count the number of particles to be stored

      do i=1,nptl
        iepo2hep(i)=0    !initialize epos index to hep index
      enddo
      do j=1,nmxhep
        ihep2epo(j)=0    !initialize hep index to epos index
        ihep2epo2(j)=0   !initialize 2nd hep index to epos index
      enddo

c  store event variables in HEP common :


c information available :
c     nrevt.......... event number
      nevhep=nrevt
c     nptevt ........ number of (stored!) particles per event
c     bimevt ........ absolute value of impact parameter
c     phievt ........ angle of impact parameter
c     kolevt ........ number of collisions
c     pmxevt ........ reference momentum
c     egyevt ........ pp cm energy (hadron) or string energy (lepton)
c     npjevt ........ number of primary projectile participants
c     ntgevt ........ number of primary target participants
c     npnevt ........ number of primary projectile neutron spectators
c     nppevt ........ number of primary projectile proton spectators
c     ntnevt ........ number of primary target neutron spectators
c     ntpevt ........ number of primary target proton spectators
c     jpnevt ........ number of absolute projectile neutron spectators
c     jppevt ........ number of absolute projectile proton spectators
c     jtnevt ........ number of absolute target neutron spectators
c     jtpevt ........ number of absolute target proton spectators

c final particles (all mother/daughter + the one absorbed in core)
      if(istmax.gt.0.and.ioclude.ne.0)then
        istmaxhep=2
      else
        istmaxhep=1
      endif


c fill first hep stack before ordering
      nhep2=0

      ipro=0
      itar=0
      pprojin=0d0
      ptargin=0d0
      imolim=0
      if(istmax.eq.0
     &   .or.(model.ne.1.and.model.ne.5.and.mod(model,4).ne.0))  !model will full list of particles
     &imolim=maproj+matarg

      do i=1,nptl

      if(iout.ge.0)then
c initialize beam momenta
        if(i.le.maproj)then
          pprojin=pprojin+dble(pptl(3,i))
        elseif(i.le.maproj+matarg)then
          ptargin=ptargin+dble(pptl(3,i))
        endif
      endif

      if(istptl(i).le.istmaxhep.or.i.le.maproj+matarg)then !store events with istptl < istmax


c  store particle variables:

c     i ............. particle number
      io=iorptl(i)
      lclean=.false.
      if(istptl(i).le.1.and.io.eq.0.and.i.gt.maproj+matarg)then
        lclean=.true.                       !happens only after cleaning
        io=i
      endif
      iadd=1
      jm1io=0
      jm2io=0
      jd1io=0
      jd2io=0
      jm1hep=0
      jm2hep=0
      jd1hep=0
      jd2hep=0
      idio=0
      if(io.gt.0)then
        if(istptl(io).le.1.and.i.gt.maproj+matarg.and..not.lclean)then!mother is normal particle (incl. spectators and fragments)
          iadd=1
          if(iepo2hep(io).gt.imolim)then
            jm1hep=iepo2hep(io)
            if(jorptl(i).gt.0)jm2hep=iepo2hep(jorptl(i))
          endif
        elseif(istmax.gt.0.and.(iorptl(io).gt.0.or.lclean))then
c     create special father/mother to have the complete chain from beam to final particle
          if(lclean)then
            istptlio=99       !if cleaning defined: use all remnants
            iorptlio=io
          else
            do while(iorptl(iorptl(io)).gt.0)
              io=iorptl(io)
            enddo
            istptlio=istptl(io)
            iorptlio=iorptl(io)
          endif
          if(istptlio.eq.41)then !remnant
            if(istptl(i).eq.2.and.jdahep2(2,iorptl(io)).eq.0)then  !beam remnant used in core
              ifrptl(1,iorptl(io))=io !no to be used again as core mother
              jmohep2(2,iorptl(io))=1
c              print *,'rcore',i,iorptl(io)
            elseif(istptl(i).le.1)then
              iadd=2
              idio=93
              jm1io=iorptl(io)
              jd1io=nhep2+2
              jm1hep=nhep2+1
              jmohep2(2,jm1io)=-1
              jdahep2(2,jm1io)=i
c              print *,'remn',i,iorptl(io)
           endif
          elseif(istptl(i).le.1.and.istptlio.eq.31)then !string
            iadd=2
            idio=92
            if(pptl(3,i).ge.0.)then
              jm1io=iorptl(io)
            else
              jm1io=jorptl(io)
            endif
            jd1io=nhep2+2
            jm1hep=nhep2+1
            jdahep2(2,jm1io)=i
            jmohep2(2,jm1io)=-1
c              print *,'string',i,iorptl(io)
          elseif(istptl(i).le.1.and.istptlio.eq.11)then !core
            iadd=2
            idio=91
            jm1io=1
            dddmn=1e33
            if(jorptl(iorptlio).eq.0)then
              jm1io=1
              dddmn=1e33
              do k=1,maproj     !look for closest projectile nucleon
                ddd=1e34
                if(iorptl(k).lt.0.and.ifrptl(1,k).eq.0)
     &               ddd=(xorptl(1,k)-xorptl(1,io))**2
     &               +(xorptl(2,k)-xorptl(2,io))**2
                if(ddd.lt.dddmn)then
                  jm1io=k
                  dddmn=ddd
                endif
              enddo
              jm2io=1
              dddmn=1e33
              do k=maproj+1,maproj+matarg !look for closest target nucleon
                ddd=1e34
                if(iorptl(k).lt.0.and.ifrptl(1,k).eq.0)
     &               ddd=(xorptl(1,k)-xorptl(1,io))**2
     &               +(xorptl(2,k)-xorptl(2,io))**2
                if(ddd.lt.dddmn)then
                  jm2io=k
                  dddmn=ddd
                endif
              enddo
              iorptl(io)=jm1io
              jorptl(io)=jm2io
              jorptl(iorptl(io))=io
              jm1io=0
              jm2io=0
            endif
            if(pptl(3,i).ge.0.)then
              jm1io=iorptl(io)
            else
              jm1io=jorptl(io)
            endif

            jd1io=nhep2+2
            jm1hep=nhep2+1
            jdahep2(2,jm1io)=i
            jmohep2(2,jm1io)=-1
c              print *,'core',i,iorptl(io)
          elseif(istptl(i).le.1.and.istptlio.eq.51)then !nuclear fragment
            idio=90
            if(jorptl(io).gt.0)then
              iadd=2
              jm1io=iorptl(io)
c              jm2io=jorptl(io)
              jm1hep=nhep2+1
              jorptl(io)=0
            else
              iadd=1
              jm1hep=iepo2hep(io)
            endif
          elseif(istptlio.eq.99)then !particle after cleaning
            iadd=2
            idio=91
            if(pptl(3,i).ge.0.)then
              ipro=ipro+1
              if(ipro.gt.maproj)ipro=1
              jm1io=ipro
            else
              itar=itar+1
              if(itar.gt.matarg)itar=1
              jm1io=maproj+itar
            endif

            jd1io=nhep2+2
            jm1hep=nhep2+1
            jdahep2(2,jm1io)=i
            jmohep2(2,jm1io)=-1
c            print *,'clean',i,iorptl(io),pptl(3,i),jm1io
          endif
        endif
      else
        jm1hep=-1
        jm2hep=-1
        if(istptl(i).eq.2)then !beam remnant used in core
          jm2hep=1
        else
          jd1hep=ifrptl(1,i)
          jd2hep=ifrptl(2,i)
        endif
      endif
      
      if(istptl(i).gt.1.and.i.gt.maproj+matarg)goto 100    !skip non final particles (allowed before to define correctly some spectators going to core)

      idpdg=idtrafo('nxs','pdg',idptl(i))
      if(idpdg.eq.99)then
        print *,'Skip particle',i,idptl(i)
        goto 100
      endif

      if(nhep2+iadd.gt.nmxhep-2)then
        print *,'Warning : produced number of particles is too high'
        print *,'          Particle list is truncated at, ', nmxhep
        print *,'          Skip event            !'
        print *,'          CHANGE HEPEVT_EntriesAllocation in your ',
     &           'HepMC library (HEPEVT_Wrapper.h) to fix this !!!!'
        goto 10000
      endif
      do j=1,iadd
        nhep2=nhep2+1
        if(j.eq.1.and.iadd.eq.2)then
          ii=io
          ix=jm1io  !for position we use father position
          id=idio
          jm1=iepo2hep(jm1io)
          jm2=0 
          if(jm2io.gt.0)jm2=iepo2hep(jm2io)
          jd1=jd1io
          jd2=jd2io
          if(id.eq.90)ihep2epo2(nhep2)=ii
        else
          ii=i
          ix=i
          id=idpdg
          jm1=jm1hep
          jm2=jm2hep
          jd1=jd1hep
          jd2=jd2hep
          if(ii.gt.maproj+matarg.or.jd1.gt.0)ihep2epo2(nhep2)=ii
        endif
        iepo2hep(ii)=nhep2
c       idptl(i) ...... particle id
        idhep2(nhep2)=id
c       pptl(1,i) ..... x-component of particle momentum (GeV/c)
        phep2(1,nhep2)=dble(pptl(1,ii))
c       pptl(2,i) ..... y-component of particle momentum (GeV/c)
        phep2(2,nhep2)=dble(pptl(2,ii))
c       pptl(3,i) ..... z-component of particle momentum (GeV/c)
        phep2(3,nhep2)=dble(pptl(3,ii))
c       pptl(4,i) ..... particle energy  (GeV)
        phep2(4,nhep2)=dble(pptl(4,ii))
c       pptl(5,i) ..... particle mass    (GeV/c2)
        phep2(5,nhep2)=dble(pptl(5,ii))
c       istptl(i) ..... generation flag: last gen. (0) or not (1)
        isthep2(nhep2)=min(2,istptl(ii)+1) !in hep:1=final, 2=decayed
        if(i.le.maproj+matarg)isthep2(nhep2)=4 !beam particles
        if(id.ge.90.and.id.le.99)isthep2(nhep2)=10+isthep2(nhep2) !intermediate state
c       ityptl(i) ..... particle type (string, remnant ...)
c       xorptl(1,i) ... x-component of formation point (fm)
        vhep2(1,nhep2)=xorptl(1,ix)*1e-12 !conversion to mm
c       xorptl(2,i) ... y-component of formation point (fm)
        vhep2(2,nhep2)=xorptl(2,ix)*1e-12 !conversion to mm
c       xorptl(3,i) ... z-component of formation point (fm)
        vhep2(3,nhep2)=xorptl(3,ix)*1e-12 !conversion to mm
c     xorptl(4,i) ... formation time (fm/c)
        vhep2(4,nhep2)=xorptl(4,ix)*1E-12 !conversion to mm/c
c       tivptl(1,i) ... formation time (always in the pp-cms!)
c       tivptl(2,i) ... destruction time (always in the pp-cms!)
c       ifrptl(1,i) ..... particle number of first daughter (no daughter=0)
        jdahep2(1,nhep2)=jd1      !need a second loop to calculated proper indice
c       ifrptl(2,i) ..... particle number of last daughter (no daughter=0)
        jdahep2(2,nhep2)=jd2      !need a second loop to calculated proper indice
c       iorptl(i) ..... particle number of father (if .le. 0 : no father)
        jmohep2(1,nhep2)=jm1
c       jorptl(i) ..... particle number of mother (if .le. 0 : no mother)
        jmohep2(2,nhep2)=jm2

c      write(ifch,130)jmohep2(1,nhep2),jmohep2(2,nhep2),nhep2
c     &,jdahep2(1,nhep2),jdahep2(2,nhep2),idhep2(nhep2)
c     &,isthep2(nhep2),(phep2(k,nhep2),k=1,5),(vhep2(k,nhep2),k=1,4)

      enddo

 100  continue

      endif
      enddo

c copy first list in final list to define daughters of beam particles 

      nhep=0
      nhepio=0
      lrcor0=.true.   !link spectator remnants to core only once 
      lrcore=.false. 

c start with beam particles (except spectators producing fragments)

c define only 2 beam particles (projectile and target)      
      if(iout.ge.0)then
c projectile
c store initial target
        if(maproj.gt.1)then
          idprin=1000000000+maproj*10+laproj*10000
        else
          idprin=idprojin
        endif
        nhep=1
        nhepio=nhepio+1
        ii=1
        id=idtrafo('nxs','pdg',idprin)
        call idmass(idprin,amass)
        idhep(nhep)=id
        phep(1,nhep)=0d0
        phep(2,nhep)=0d0
        phep(3,nhep)=pprojin
        phep(4,nhep)=sqrt(pprojin**2+dble(amass)**2)
        phep(5,nhep)=dble(amass)
        isthep(nhep)=4 !in hep:beam particle
        vhep(1,nhep)=0d0  !Main vertex at 0.
        vhep(2,nhep)=0d0
        vhep(3,nhep)=0d0
        vhep(4,nhep)=0d0
        jdahep(1,nhep)=3
        jdahep(2,nhep)=maproj+matarg+2  !updated later if needed
        jmohep(1,nhep)=-1
        jmohep(2,nhep)=-1
c Target
        if(matarg.gt.1)then
          idtgin=1000000000+matarg*10+latarg*10000
        else
          idtgin=idtargin
        endif
        nhep=2
        nhepio=nhepio+1
        ii=maproj+1
        id=idtrafo('nxs','pdg',idtgin)
        call idmass(idtgin,amass)
        idhep(nhep)=id
        phep(1,nhep)=0d0
        phep(2,nhep)=0d0
        phep(3,nhep)=ptargin
        phep(4,nhep)=sqrt(ptargin**2+dble(amass)**2)
        phep(5,nhep)=dble(amass)
        isthep(nhep)=4 !in hep:beam particle
        vhep(1,nhep)=0d0  !Main vertex at 0.
        vhep(2,nhep)=0d0
        vhep(3,nhep)=0d0
        vhep(4,nhep)=0d0
        jdahep(1,nhep)=3
        jdahep(2,nhep)=maproj+matarg+2  !updated later if needed
        jmohep(1,nhep)=-1
        jmohep(2,nhep)=-1
      endif

      if(imolim.eq.0)then
        nhep=maproj+matarg   !beam particles will be copied later
        if(iout.ge.0)nhep=nhep+2   !add new beam particles
      endif

      nhep0=0
      nskip=0
      if(iout.ge.0)nskip=-2     !add 2 beam particles and may be skip nucleons
          
c reorder individual beam particles to make list of mothers for a given daughter
      do j=1,maproj+matarg

        if(imolim.ne.0)then

          if(iout.lt.0)then

c when no daughter/mother informations, simply copy beam particles
          nhep=nhep+1
          nhepio=nhepio+1
          idhep(nhep)=idhep2(j)
          phep(1,nhep)=phep2(1,j)
          phep(2,nhep)=phep2(2,j)
          phep(3,nhep)=phep2(3,j)
          phep(4,nhep)=phep2(4,j)
          phep(5,nhep)=phep2(5,j)
          isthep(nhep)=4
          vhep(1,nhep)=vhep2(1,j)
          vhep(2,nhep)=vhep2(2,j)
          vhep(3,nhep)=vhep2(3,j)
          vhep(4,nhep)=vhep2(4,j)
          jmohep(1,nhep)=-1
          jmohep(2,nhep)=-1
          jdahep(1,nhep)=0
          jdahep(2,nhep)=0
          iepo2hep(j)=nhep
          isthep2(j)=-isthep2(j) 
          ihep2epo2(j)=-nhep

          else

c skip individual beam particles in case of short list
          nskip=nskip+1

          endif

        else

        nhep0=nhep       !index of daughter list of current beam particle
        nhepi0=nhepio+1
        isthep(nhepi0)=0
        nio=0

c copy all daughters after the mother
        do k=maproj+matarg+1,nhep2

          if(jmohep2(1,k).eq.j.and.idhep2(k).ne.90)then
            if(isthep(nhepi0).eq.0)then       !first save mother beam particle
              nhepio=nhepio+1
              idhep(nhepio)=idhep2(j)
              phep(1,nhepio)=phep2(1,j)
              phep(2,nhepio)=phep2(2,j)
              phep(3,nhepio)=phep2(3,j)
              phep(4,nhepio)=phep2(4,j)
              phep(5,nhepio)=phep2(5,j)
              if(iout.lt.0)then
                isthep(nhepio)=4   !beam
                jmohep(1,nhepio)=-1
                jmohep(2,nhepio)=-1
              else
                isthep(nhepio)=14   !daughter of beam
                jmohep(1,nhepio)=1
                jmohep(2,nhepio)=2
              endif
              vhep(1,nhepio)=vhep2(1,j)
              vhep(2,nhepio)=vhep2(2,j)
              vhep(3,nhepio)=vhep2(3,j)
              vhep(4,nhepio)=vhep2(4,j)
              iepo2hep(j)=nhepio
              isthep2(j)=-isthep2(j) 
              ihep2epo2(j)=-nhepio
              lrcore=.false.     !link spectator remnants to core only once 
              if(lrcor0)then
               kk=k
               do while (kk.le.nhep2.and..not.lrcore)
                if(idhep2(kk).eq.91.and.jmohep2(1,kk).eq.j)lrcore=.true.
                kk=kk+1
               enddo
              endif
            endif
            if(lrcore)then           !save other mothers for same core
              lrcor0=.false.
              do i=1,maproj+matarg
                if(isthep2(i).gt.0.and.jmohep2(2,i).gt.0)then
                  nhepio=nhepio+1
                  nio=nio+1
                  idhep(nhepio)=idhep2(i)
                  phep(1,nhepio)=phep2(1,i)
                  phep(2,nhepio)=phep2(2,i)
                  phep(3,nhepio)=phep2(3,i)
                  phep(4,nhepio)=phep2(4,i)
                  phep(5,nhepio)=phep2(5,i)
                  if(iout.lt.0)then
                    isthep(nhepio)=4 !beam
                    jmohep(1,nhepio)=-1
                    jmohep(2,nhepio)=-1
                  else
                    isthep(nhepio)=14 !daughter of beam
                    jmohep(1,nhepio)=1
                    jmohep(2,nhepio)=2
                  endif
                  vhep(1,nhepio)=vhep2(1,i)
                  vhep(2,nhepio)=vhep2(2,i)
                  vhep(3,nhepio)=vhep2(3,i)
                  vhep(4,nhepio)=vhep2(4,i)
                  isthep2(i)=-isthep2(i) 
                  iepo2hep(i)=nhepio
                  ihep2epo2(i)=-nhepio
                endif
              enddo
            endif
            nhep=nhep+1
            idhep(nhep)=idhep2(k)
            phep(1,nhep)=phep2(1,k)
            phep(2,nhep)=phep2(2,k)
            phep(3,nhep)=phep2(3,k)
            phep(4,nhep)=phep2(4,k)
            phep(5,nhep)=phep2(5,k)
            isthep(nhep)=isthep2(k)
            vhep(1,nhep)=vhep2(1,k)
            vhep(2,nhep)=vhep2(2,k)
            vhep(3,nhep)=vhep2(3,k)
            vhep(4,nhep)=vhep2(4,k)
            jdahep(1,nhep)=0
            if(jdahep2(1,k).gt.0)jdahep(1,nhep)=-ihep2epo2(jdahep2(1,k))
            jdahep(2,nhep)=0
            if(jdahep2(2,k).gt.0)jdahep(2,nhep)=-ihep2epo2(jdahep2(2,k))
            jmohep(1,nhep)=nhepi0
            jmohep(2,nhep)=nhepi0+nio
            isthep2(k)=-isthep2(k)
            if(ihep2epo2(k).le.0)then
              ihep2epo2(k)=-nhep
            else   !for nuclear fragments
              ihep2epo(nhep)=ihep2epo2(k)
              if(ihep2epo(nhep).gt.0)iepo2hep(ihep2epo(nhep))=nhep
            endif
              
          endif

        enddo

        if(nhepio.ge.nhepi0)then
          do i=nhepi0,nhepi0+nio
            jdahep(1,i)=nhep0+1
            jdahep(2,i)=nhep
          enddo
        endif

        endif

      enddo

      iround=1
      if(iout.ge.0.and.nhep0.eq.0)iround=0    !copy first all mothers and then daughters if not full list from EPOS (in that case nhep0.ne.0)

      do iii=iround,1

c copy all other particles (secondary particles and spectators)

      do k=maproj+matarg+1,nhep2

          if(isthep2(k).gt.0.and.(iii.eq.1
     &                       .or.(iii.eq.0.and.jmohep2(1,k).le.0)))then

c look for mother of fragments
            nhepi0=nhepio+1
            if(nhep0.ne.0.and.jmohep2(1,k).gt.0
     &                   .and.jmohep2(1,k).le.maproj+matarg)then
c copy all mothers before the daughter to complete the beam remnant list
              do j=1,maproj+matarg

                if(jdahep2(1,j).eq.ihep2epo2(k).and.isthep2(j).gt.0)then
                  nhepio=nhepio+1
                  idhep(nhepio)=idhep2(j)
                  phep(1,nhepio)=phep2(1,j)
                  phep(2,nhepio)=phep2(2,j)
                  phep(3,nhepio)=phep2(3,j)
                  phep(4,nhepio)=phep2(4,j)
                  phep(5,nhepio)=phep2(5,j)
                  if(iout.lt.0)then
                    isthep(nhepio)=4 !beam
                    jmohep(1,nhepio)=-1
                    jmohep(2,nhepio)=-1
                  else
                    isthep(nhepio)=3 !daughter of beam
                    jmohep(1,nhepio)=1
                    jmohep(2,nhepio)=2
                  endif
                  vhep(1,nhepio)=vhep2(1,j)
                  vhep(2,nhepio)=vhep2(2,j)
                  vhep(3,nhepio)=vhep2(3,j)
                  vhep(4,nhepio)=vhep2(4,j)
                  jdahep(1,nhepio)=0
                  if(jdahep2(1,j).gt.0)jdahep(1,nhepio)=-jdahep2(1,j)
                  jdahep(2,nhepio)=0
                  if(jdahep2(2,j).gt.0)jdahep(2,nhepio)=-jdahep2(2,j)
                  isthep2(j)=-isthep2(j) 
                  iepo2hep(j)=nhepio
                  ihep2epo(nhepio)=j
                endif

              enddo

            endif

            nhep=nhep+1
            idhep(nhep)=idhep2(k)
            phep(1,nhep)=phep2(1,k)
            phep(2,nhep)=phep2(2,k)
            phep(3,nhep)=phep2(3,k)
            phep(4,nhep)=phep2(4,k)
            phep(5,nhep)=phep2(5,k)
            isthep(nhep)=isthep2(k)
            vhep(1,nhep)=vhep2(1,k)
            vhep(2,nhep)=vhep2(2,k)
            vhep(3,nhep)=vhep2(3,k)
            vhep(4,nhep)=vhep2(4,k)
            jdahep(1,nhep)=0
            if(jdahep2(1,k).gt.0)jdahep(1,nhep)=-ihep2epo2(jdahep2(1,k))
            jdahep(2,nhep)=0
            if(jdahep2(2,k).gt.0)jdahep(2,nhep)=-ihep2epo2(jdahep2(2,k))
            ihep2epo(nhep)=ihep2epo2(k)
            if(ihep2epo(nhep).gt.0)iepo2hep(ihep2epo(nhep))=nhep
            if(iii.eq.0)then     !this particle is a first mother directly link to beam
              jmohep(1,nhep)=1
              jmohep(2,nhep)=2
              jdahep(2,1)=nhep
              jdahep(2,2)=nhep
            else
              if(nhepio.lt.nhepi0)then
                if(jmohep2(1,k).gt.0)then
                  if(ihep2epo2(jmohep2(1,k)).le.0)then
                    jmohep(1,nhep)=-ihep2epo2(jmohep2(1,k))
                  else
                    jmohep(1,nhep)=iepo2hep(ihep2epo2(jmohep2(1,k)))
                  endif
                else
                  jmohep(1,nhep)=0
                endif
                if(jmohep2(2,k).gt.0)then
                  if(ihep2epo2(jmohep2(2,k)).le.0)then
                    jmohep(2,nhep)=-ihep2epo2(jmohep2(2,k))
                  else
                    jmohep(2,nhep)=iepo2hep(ihep2epo2(jmohep2(2,k)))
                  endif
                else
                  jmohep(2,nhep)=0
                endif
              else              !for nuclear fragments
                jmohep(1,nhep)=nhepi0
                jmohep(2,nhep)=nhepio             
              endif
            endif
            isthep2(k)=-isthep2(k)   !mark particle as copied to final list

          endif

      enddo

      enddo  !iround


      if(nhep+nskip.ne.nhep2.or.nhepio+nskip.ne.maproj+matarg)then
        print *,'Warning : number of particles changed after copy'
     &         ,nhepio+nskip,maproj+matarg,nskip
        nrem1=0
        do k=1,nhep2
        if(isthep2(k).eq.-4)then
          nrem1=nrem1+1
        endif
        if(isthep2(k).gt.0)
     &  print *,'         ',k,idhep2(k),jmohep2(1,k),isthep2(k)

     &         ,'from',ihep2epo2(k)
        enddo
        print *,'         ',nhep2-nskip,'->',nhep
        nrem2=0
        do k=1,nhep
          if((iout.ge.0.and.isthep(k).eq.3)
     &   .or.(iout.lt.0.and.isthep(k).eq.4))then
            nrem2=nrem2+1
          endif
        print *,'         ',k,idhep(k),isthep(k),'from',ihep2epo(k)
        enddo
        print *,'          Particle list not consistent, skip event !'
        print *,'         ',nrem1,'->',nrem2
c        stop
        goto 10000
      endif


c update daughter list with correct index

      do j=1,nhep


        i=ihep2epo(j)

        if(jdahep(1,j).lt.0)then

          jdahep(1,j)=iepo2hep(-jdahep(1,j))
          if(jdahep(2,j).lt.0)jdahep(2,j)=iepo2hep(-jdahep(2,j))
          

        elseif(i.gt.0.and.jdahep(1,j).eq.0)then

          ifr1=ifrptl(1,i)
          ifr2=ifrptl(2,i)
c         ifrptl(1,i) ..... particle number of first daughter (no daughter=0)
          if(ifr1.gt.0)then
            if(iepo2hep(ifr1).gt.0)then
              jdahep(1,j)=iepo2hep(ifr1)
            elseif(ifr2.gt.0)then  
c if first daughter not in the final list look for first finally saved daughter
              do while(iepo2hep(ifr1).eq.0.and.ifr1.lt.ifr2)
                ifr1=ifr1+1
              enddo
              jdahep(1,j)=iepo2hep(ifr1)              
            endif
          else
            jdahep(1,j)=0
          endif

          if(jdahep(2,j).eq.0)then
c         ifrptl(2,i) ..... particle number of last daughter (no daughter=0)
            if(ifr2.gt.0)then
              if(iepo2hep(ifr2).gt.0)then
                jdahep(2,j)=iepo2hep(ifr2)
              elseif(ifr1.gt.0)then  
c if last daughter not in the final list look for last finally saved daughter
                do while(iepo2hep(ifr2).eq.0.and.ifr1.lt.ifr2)
                  ifr2=ifr2-1
                enddo
                jdahep(2,j)=iepo2hep(ifr2)              
              endif
            else
              jdahep(2,j)=0
            endif
          endif

        endif

c      write(ifch,130)jmohep(1,j),jmohep(2,j),j,jdahep(1,j),jdahep(2,j)
c     &,idhep(j),isthep(j),(phep(k,j),k=1,5),(vhep(k,j),k=1,4)
c 130  format (1x,i6,i6,3x,i6,3x,i6,i6,i12,i4,8x,5(e8.2,1x)
c     *,4x,4(e8.2,1x))
        
      enddo


 9999 return
10000 nhep=0
      goto 9999
      end

c-----------------------------------------------------------------------
      subroutine lhestore(n)
c-----------------------------------------------------------------------
c     writes the results of a simulation into the file with unit ifdt
c     contains a description of the stored variables.
c     use Les Houches Event File as defined in hep-ph/0109068 for the
c     common block and hep-ph/0609017 for the XML output.
c     some code taken from example from Torbjrn Sjstrand
c     in http://www.thep.lu.se/~torbjorn/lhef
c-----------------------------------------------------------------------
#include "aaa.h"
 
C...User process event common block.
      INTEGER MAXNUP
      PARAMETER (MAXNUP=50000)  !extend array for file production
c      PARAMETER (MAXNUP=500)
      INTEGER NUP,IDPRUP,IDUP,ISTUP,MOTHUP,ICOLUP
      DOUBLE PRECISION XWGTUP,SCALUP,AQEDUP,AQCDUP,PUP,VTIMUP,SPINUP
      COMMON/HEPEUP/NUP,IDPRUP,XWGTUP,SCALUP,AQEDUP,AQCDUP,IDUP(MAXNUP),
     &ISTUP(MAXNUP),MOTHUP(2,MAXNUP),ICOLUP(2,MAXNUP),PUP(5,MAXNUP),
     &VTIMUP(MAXNUP),SPINUP(MAXNUP)
      SAVE /HEPEUP/


      integer iepo2hep(mxptl)

c  count the number of particles to be stored (--> nptevt)

      nhep=0
      do i=1,nptl
        if(istptl(i).le.istmax.and.abs(idptl(i)).le.10000)nhep=nhep+1
      enddo
      if(nhep.gt.MAXNUP)then
        print *,'Warning : produced number of particles is too high'
        print *,'          event is not stored'
        goto 1000
      endif

C...set event info and get number of particles.
      NUP=nhep             !number of particles
      IDPRUP=nint(abs(typevt))  !type of event (ND,DD,CD,SD)
      XWGTUP=1d0           !weight of event
      SCALUP=-1d0          !scale for PDF (not used)
      AQEDUP=-1d0          !alpha QED (not relevant)
      AQCDUP=-1d0          !alpha QCD (not relevant)

C...Copy event lines, omitting trailing blanks. 
C...Embed in <event> ... </event> block.
      write(ifdt,'(A)') '<event>' 
      write(ifdt,*)NUP,IDPRUP,XWGTUP,SCALUP,AQEDUP,AQCDUP
      nhep=0
      DO 220 i=1,nptl

        if(istptl(i).le.istmax.and.abs(idptl(i)).le.10000)then !store events with istptl < istmax

          nhep=nhep+1
c     i ............. particle number
c     idptl(i) ...... particle id
          idpdg=idtrafo('nxs','pdg',idptl(i))
          if(idpdg.eq.99)idpdg=0   !unknown particle
          iepo2hep(i)=nhep
c  store particle variables:
          IDUP(nhep)=idpdg
          if(iorptl(i).lt.0)then
          ISTUP(nhep)=-9      !incoming particle
          else
          ISTUP(nhep)=min(3,istptl(i)+1) !in LHEF:1=final, 2=decayed, 3=intermediate state
          endif
          if(iorptl(i).gt.0)then
            MOTHUP(1,nhep)=iepo2hep(iorptl(i))
          else
            MOTHUP(1,nhep)=0
          endif
c     jorptl(i) ..... particle number of mother (if .le. 0 : no mother)
          if(jorptl(i).gt.0)then
            MOTHUP(2,nhep)=iepo2hep(jorptl(i))
          else
            MOTHUP(2,nhep)=-1
          endif
          ICOLUP(1,nhep)=0        !color flow
          ICOLUP(2,nhep)=0        !color flow
          do J=1,5                !particle momentum (GeV/c)
            PUP(J,nhep)=dble(pptl(J,i))
          enddo
          VTIMUP(nhep)=(dble(tivptl(2,i))-dble(tivptl(1,i)))*1d-12 !life time c*tau in mm
          if(VTIMUP(nhep).gt.dble(ainfin)
     &   .or.VTIMUP(nhep).ne.VTIMUP(nhep))then
            write(ifch,*)'ici',VTIMUP(nhep),tivptl(2,i),tivptl(1,i)
     &                        ,i,nptl
            VTIMUP(nhep)=ainfin
            call utstop("aie&")
          endif
          SPINUP(nhep)=9           !polarization (not known)
          write(ifdt,*)IDUP(nhep),ISTUP(nhep),
     &      MOTHUP(1,nhep),MOTHUP(2,nhep),ICOLUP(1,nhep),ICOLUP(2,nhep),
     &      (PUP(J,nhep),J=1,5),VTIMUP(nhep),SPINUP(nhep)
        endif
  220 CONTINUE

c optional informations
      write(ifdt,*)'#geometry',bimevt,phievt

      write(ifdt,'(A)') '</event>' 

      if(n.eq.nevent)then
C...Successfully reached end of event loop: write closing tag
        write(ifdt,'(A)') '</LesHouchesEvents>' 
        write(ifdt,'(A)') ' ' 
      endif

 1000 continue

      return
      end


