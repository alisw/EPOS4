      subroutine oneHQevent(icol,totcol,timelim,ifail)

      use PBGhqdmesons
      use PBGhquarks
      use PBGmainpro1
      use PBGmainpro2
      use PBGnpartmax
      use PBGsystem
      use PBGsystem2
      use PBGhqprod
      use PBGinitnb
      use PBGgenvar
      use PBGsceplas
      use PBGhydro
      use PBGbookpsi
      use denysvar
      use genuineRemler
      use analcol
      use onerateHQparam
      implicit none
      integer ifail,ifailhyd
      logical first,ifsing,check
      data first/.true./
      save first
      integer i,i2,itime,ntime,k,ipart,ipart2,ires,iresnew
      integer icol,totcol,imed,nbclose
      integer partcount1,partcount2,nbsing
      double precision time,timelim
      integer nbtypeires(nres)
      double precision p_part(0:3),x_part(0:3),temp,densener,
     & betaplasma(3),pT,rapid,ran2
      real cphi,sphi
      call reseteventlist
      call resetparticlelist()
      if(isce==7) then
         if(.not.(ifhqfromhq)) then
            call def_part_init_epos3
         else
            call def_part_init_epos_512(ifcronin)
         endif
      else
         call def_part_init(tauf,avnc0ncoll,avnpsi0ncoll,avnb0ncoll,
     &          avnups0ncoll,A,idistspatial,ifcronin)
      endif
      taustep=timestep
      nbhqc=0
      nbhqb=0
      ipart=firstpart
      do while (ipart.NE.-1)
        if((partinfo_ires(ipart).eq.1).or.(partinfo_ires(ipart).eq.2)) 
     &        nbhqc=nbhqc+1
        if((partinfo_ires(ipart).eq.14).or.(partinfo_ires(ipart).eq.15)) 
     &        nbhqb=nbhqb+1
        ipart=partinfo_next(ipart)
      enddo
      nbhq=nbhqc+nbhqb
      allocate(typquarkdyn(nbhq))
      allocate(hqmother(nbhq))      
      allocate(timelimdyn(nbhq))
      allocate(ifcorona(nbhq))
      allocate(iresdyn(nbhq,0:nbtausteps))
      allocate(xdyn(nbhq,0:nbtausteps,0:3))
      allocate(pdyn(nbhq,0:nbtausteps,0:3))
      allocate(tempatHQpos(nbhq,0:nbtausteps))  
      allocate(imedatHQpos(nbhq,0:nbtausteps))  
      allocate(ifsinglet(nbhq/2,nbhq/2))
      allocate(ifsingletfull(nbhq,nbhq))
      do i=1,nbhq/2
         do i2=1,nbhq/2
            ifsinglet(i,i2)=.false.
         enddo
      enddo
      do i=1,nbhq
         do i2=1,nbhq
            ifsingletfull(i,i2)=.false.
         enddo
      enddo
      ipart=firstpart
      i=0
      partcount1=0
      nbsing=0
      iselec=0
      do while (ipart.NE.-1)
         partcount1=partcount1+1
         if((partinfo_ires(ipart).eq.1).or.(partinfo_ires(ipart).eq.14))
     &        i=i+1
         ipart2=firstpart
         i2=0
         partcount2=0
         do while (ipart2.NE.-1)
            partcount2=partcount2+1
            if((partinfo_ires(ipart2).eq.2).or.
     &           (partinfo_ires(ipart2).eq.15)) i2=i2+1
            if(((partinfo_ires(ipart).eq.1).or.
     &           (partinfo_ires(ipart).eq.14)).and.
     &        ((partinfo_ires(ipart2).eq.2).or.
     &           (partinfo_ires(ipart2).eq.15))) then
               ifsing=(ran2().lt.0.1111111111111)
               ifsinglet(i,i2)=ifsing
               ifsingletfull(partcount1,partcount2)=ifsing
               ifsingletfull(partcount2,partcount1)=ifsing
               if(ifsing) nbsing=nbsing+1
            endif
            ipart2=partinfo_next(ipart2)
         enddo
         pT=partinfo_p(ipart,1)**2+partinfo_p(ipart,2)**2
         if((pT.lt.25d0).and.partinfo_ires(ipart).le.2) then
            if(partinfo_p(ipart,3).gt.0.d0) then
               rapid=log((partinfo_p(ipart,0)+partinfo_p(ipart,3))/
     &              sqrt(pT+partinfo_mass(ipart)**2))
            else
               rapid=-log((partinfo_p(ipart,0)-partinfo_p(ipart,3))/
     &              sqrt(pT+partinfo_mass(ipart)**2))
            endif
            pT=sqrt(pT)
            if((abs(rapid).lt.1.5d0).and.iselec.lt.10) then
               iselec=iselec+1
               selectedparticles(iselec)=partcount1
            endif
         endif
         ipart=partinfo_next(ipart)
      enddo
      i=0
      ipart=firstpart
      do while (ipart.NE.-1)
         if(resinfo_ifquark(partinfo_ires(ipart))) then
            i=i+1
            if(partinfo_ires(ipart).eq.1) typquarkdyn(i)=1
            if(partinfo_ires(ipart).eq.2) typquarkdyn(i)=-1
            if(partinfo_ires(ipart).eq.14) typquarkdyn(i)=2
            if(partinfo_ires(ipart).eq.15) typquarkdyn(i)=-2
            hqmother(i)=partinfo_mother(ipart)
            tempatHQpos(i,0)=0.d0
            imedatHQpos(i,0)=-1
            iresdyn(i,0)=partinfo_ires(ipart)
            do k=0,3
               xdyn(i,0,k)=partinfo_r(ipart,k)
               pdyn(i,0,k)=partinfo_p(ipart,k)
            enddo
         endif
         ipart=partinfo_next(ipart)
      enddo
      time=tauf
      ntime=1
      call freefly(1.001*tau0)

      time=1.001*tau0
      itime=1
      i=0
      ipart=firstpart
      do while (ipart.NE.-1)
         ires=partinfo_ires(ipart)
         if((ires.ge.1).and.(ires.le.26).and.(ires.ne.3).and.
     &        (ires.ne.16)) then
            i=i+1
            do k=0,3
               x_part(k)=partinfo_r(ipart,k)
            enddo
            call givelocalhydro(x_part,temp,betaplasma,densener,
     &           ifailhyd,imed)
            tempatHQpos(i,itime)=temp
            imedatHQpos(i,itime)=imed
            iresdyn(i,itime)=partinfo_ires(ipart)
            do k=0,3
               xdyn(i,itime,k)=partinfo_r(ipart,k)
               pdyn(i,itime,k)=partinfo_p(ipart,k)
            enddo
            if(.not.(resinfo_ifquark(partinfo_ires(ipart)))) then
               timelimdyn(i)=time
            else
               timelimdyn(i)=1.D8
            endif
         endif
         ipart=partinfo_next(ipart)
      enddo
      nbhardcol=0
      allocate(ipartcol(nbhardcolmax))
      allocate(tempcol(nbhardcolmax))
      allocate(timecol(nbhardcolmax))
      allocate(dpcol(nbhardcolmax,0:3))
      if(.not.(ifinterpolation)) then
         allocate(pcolafter(nbhardcolmax,0:3))
      endif
      ifail=0
      do while(((time+timestep).le.timelim).and.(ifail.eq.0))
         if(infotrack.ge.1) then
            write(8,*) '--------------------------------------' 
            write(8,*) 'time Bjorken = ',time
         endif   
         time=time+timestep
         ntime=ntime+1
         itime=itime+1
 
         call update_part_chain(time,nbclose,ifail)
         if(ifail.ne.0) then
            write(6,*) 'skipping event'
            return
         endif
         if(itime.le.nbtausteps0) then
            avnbclose(itime)=avnbclose(itime)+nbclose
         endif
         i=0
         ipart=firstpart
         do while (ipart.NE.-1)
            ires=partinfo_ires(ipart)
         if((ires.ge.1).and.(ires.le.26).and.(ires.ne.3).and.
     &        (ires.ne.16)) then
               i=i+1
               do k=0,3
                  x_part(k)=partinfo_r(ipart,k)
               enddo
               call givelocalhydro(x_part,temp,betaplasma,densener,
     &              ifailhyd,imed)
               tempatHQpos(i,itime)=temp
               imedatHQpos(i,itime)=imed
               iresdyn(i,itime)=partinfo_ires(ipart)
               do k=0,3
                  xdyn(i,itime,k)=partinfo_r(ipart,k)
                  pdyn(i,itime,k)=partinfo_p(ipart,k)
               enddo
               if(resinfo_ifquark(ires).and.
     &              (.not.resinfo_ifquark(iresdyn(i,itime-1))))
     &              then
                  timelimdyn(i)=1.D8
               endif
               if(.not.resinfo_ifquark(ires)) then
                  if(xdyn(i,itime,0).lt.timelimdyn(i)) then
                     timelimdyn(i)=xdyn(i,itime-1,0)
                  endif
               endif
            endif
            ipart=partinfo_next(ipart)
         enddo
           call count_part_type(nbtypeires)
           call reseteventlist
        enddo
      i=0
      nbhqincorona=0
      ipart=firstpart
      do while (ipart.NE.-1)
         ires=partinfo_ires(ipart)
         i=i+1
         if(partinfo_corona(ipart)) then
            ifcorona(i)=.true.
            nbhqincorona=nbhqincorona+1
         else
            ifcorona(i)=.false.
         endif
         if((ires.eq.4).or.(ires.eq.5).or.(ires.eq.6).or.
     & (ires.eq.7).or.(ires.eq.8).or.(ires.eq.9).or.(ires.eq.10).or.
     & (ires.eq.11).or.(ires.eq.12).or.(ires.eq.13).or.(ires.eq.17).or.
     & (ires.eq.18).or.(ires.eq.19).or.
     & (ires.eq.20).or.(ires.eq.21).or.(ires.eq.22).or.(ires.eq.23).or.
     & (ires.eq.24).or.(ires.eq.25).or.(ires.eq.26)) then   
         else
            write(6,*) 'unrecognized ires:',ires
            write(6,*) 'for ipart:',ipart
            ifail=1
            return
         endif         
         ipart=partinfo_next(ipart)                      
      enddo
      ifail=0
      ipart=firstpart
      i=0
      do while (ipart.NE.-1)
         check=.false.
         time=0.d0
         itime=0
         i=i+1
         ires=iresdyn(i,itime)
         do while(time+timestep.le.timelim)
            time=time+timestep
            itime=itime+1            
            iresnew=iresdyn(i,itime)
            if(iresnew.eq.ires-3) then
               check=.true.
            endif
            if((iresnew.eq.ires+3).and.check) then
            endif
            ires=iresnew
         enddo
         ipart=partinfo_next(ipart)                      
      enddo
      if(icol.le.nbcolsaveEPOS) then
        if(nHQ+npart.gt.npartmax) then
           ifail=-1
           return
        endif
        ifail=0
        ipart=firstpart
        cphi=cos(phi+phir)
        sphi=sin(phi+phir)
        do while (ipart.NE.-1)
           ires=partinfo_ires(ipart)
           if((ires.eq.4).or.(ires.eq.5).or.(ires.eq.6).or.(ires.eq.7)
     & .or.(ires.eq.8).or.(ires.eq.9).or.(ires.eq.10).or.(ires.eq.11)
     & .or.(ires.eq.12).or.(ires.eq.13).or.(ires.eq.17).or.(ires.eq.18)
     & .or.(ires.eq.19).or.(ires.eq.20).or.(ires.eq.21).or.(ires.eq.22)
     & .or.(ires.eq.23).or.(ires.eq.24).or.(ires.eq.25)
     & .or.(ires.eq.26)) then
              nHQ=nHQ+1
              nd=nd+1   
              do k=0,3
                 p_part(k)=partinfo_p(ipart,k)
                 x_part(k)=partinfo_rghost(ipart,k)
              enddo
            if((ires.eq.4).or.(ires.eq.6).or.(ires.eq.8).or.
     & (ires.eq.10).or.(ires.eq.12))then
              idhquarks(nHQ)=1
            elseif((ires.eq.5).or.(ires.eq.7).or.(ires.eq.9).or.
     & (ires.eq.11).or.(ires.eq.13))then
              idhquarks(nHQ)=2
            elseif((ires.eq.17).or.(ires.eq.19).or.(ires.eq.21).or.
     & (ires.eq.23).or.(ires.eq.25))then
              idhquarks(nHQ)=14
            elseif((ires.eq.18).or.(ires.eq.20).or.(ires.eq.22).or.
     & (ires.eq.24).or.(ires.eq.26))then
              idhquarks(nHQ)=15
            endif

              phquarks(nHQ,0)=partinfo_pghost(ipart,0)
              phquarks(nHQ,3)=partinfo_pghost(ipart,3)
              phquarks(nHQ,1)=cphi*partinfo_pghost(ipart,1)-
     &             sphi*partinfo_pghost(ipart,2)
              phquarks(nHQ,2)=sphi*partinfo_pghost(ipart,1)+
     &             cphi*partinfo_pghost(ipart,2)
C              phquarks(nHQ,2)=sphi*partinfo_pghost(ipart,1)-
C     &             cphi*partinfo_pghost(ipart,2)
C BIG BUG IN v0
         if((ires.eq.4).or.(ires.eq.5).or.(ires.eq.6).or.
     & (ires.eq.7).or.(ires.eq.8).or.(ires.eq.9).or.(ires.eq.10).or.
     & (ires.eq.11).or.(ires.eq.12).or.(ires.eq.13)) then
                 mhquarks(nHQ)=mc
         elseif((ires.eq.17).or.(ires.eq.18).or.(ires.eq.19).or.
     &(ires.eq.20).or.(ires.eq.21).or.(ires.eq.22).or.(ires.eq.23).or.
     &(ires.eq.24).or.(ires.eq.25).or.(ires.eq.26)) then
                 mhquarks(nHQ)=mb
              else
                 call utstop('unrecognized ires in oneHQevent&')
              endif
              hqmeson(nd,1)=ires
              hqmeson(nd,2)=p_part(0)
              HQmeson(nd,3)=cphi*p_part(1)-sphi*p_part(2)
              HQmeson(nd,4)=sphi*p_part(1)+cphi*p_part(2)
              HQmeson(nd,5)=p_part(3)
              HQmeson(nd,6)=x_part(0)
              HQmeson(nd,7)=cphi*x_part(1)-sphi*x_part(2)
              HQmeson(nd,8)=sphi*x_part(1)+cphi*x_part(2)
              HQmeson(nd,9)=x_part(3)
           endif
           ipart=partinfo_next(ipart) 
        enddo
      endif
      call pack_it
      if((totcol/1000)*1000.eq.totcol) then
         write(6,*) 'event:',totcol
      endif
      end subroutine oneHQevent


      subroutine analysePSdensity(time,timestep)
      use denysvar
      use PSdens
      use PBGgenvar
      implicit none
      integer ipart,ipart2,nquark,nquarkbar,nbsing,j,k,it,irjpsi,ipjpsi,
     & partcount1,partcount2
      double precision time,timestep,pQ(0:3),pQbar(0:3),xQ(0:3),
     & xQbar(0:3),timeround,yQ,yQbar,mQ,mQbar,sigma,wignerrelat,pT,phi,
     & rapid,rjpsi,pjpsi,gamma,ergcm
      logical anal,first
      parameter(sigma=0.35d0)
      data first/.true./
      save first
      if(first) then
         allocate(PSdensnative(10,100,50))
         allocate(PSdensnativehighpT(10,100,50))
         allocate(avnbclosewigner(10))
         do it=1,10
            avnbclosewigner(it)=0.d0
            do j=1,100
               do k=1,50
                  PSdensnative(it,j,k)=0.d0
                  PSdensnativehighpT(it,j,k)=0.d0
               enddo
            enddo
         enddo
         first=.false.
      endif
      timeround=dble(int(time))
      anal=(time.ge.timeround).and.(time-timestep.lt.timeround).and.
     &  (timeround.ge.1.d0).and.(timeround.le.10.d0)
      if(.not.(anal)) return
      it=nint(timeround)
      ipart=firstpart
      nquark=0
      partcount1=0
      nbsing=0
      do while (ipart.NE.-1)
         partcount1=partcount1+1
         if(partinfo_ires(ipart).eq.1) nquark=nquark+1
         ipart2=firstpart
         nquarkbar=0
         partcount2=0
         do while (ipart2.NE.-1)
            partcount2=partcount2+1
            if(partinfo_ires(ipart2).eq.2) nquarkbar=nquarkbar+1
            if(ifsingletfull(partcount1,partcount2)) then
               nbsing=nbsing+1
               do k=0,3
                  pQ(k)=partinfo_p(ipart,k)
                  pQbar(k)=partinfo_p(ipart2,k)
                  xQ(k)=partinfo_r(ipart,k)
                  xQbar(k)=partinfo_r(ipart2,k)
               enddo   
               mQ=partinfo_mass(ipart)
               mQbar=partinfo_mass(ipart2)
               if(pQ(3).gt.0.d0) then
                  yQ=log((pQ(0)+pQ(3))/sqrt(mQ**2+pQ(1)**2+pQ(2)**2))
               else
                  yQ=-log((pQ(0)-pQ(3))/sqrt(mQ**2+pQ(1)**2+pQ(2)**2))
               endif
               if(pQbar(3).gt.0.d0) then
                  yQbar=log((pQbar(0)+pQbar(3))/
     &                 sqrt(mQbar**2+pQbar(1)**2+pQbar(2)**2))
               else
                  yQbar=-log((pQbar(0)-pQbar(3))/
     &                 sqrt(mQbar**2+pQbar(1)**2+pQbar(2)**2))
               endif
               if((abs(yQ).le.1.5d0).and.(abs(yQbar).le.1.5d0)) then
                  call gimmewignerrelatsinglet(pQ,pQbar,xQ,xQbar,sigma,
     &                 wignerrelat,pT,phi,rapid,rjpsi,pjpsi,gamma,ergcm)  
                  if(rjpsi.le.1.d0) then
                     avnbclosewigner(it)=avnbclosewigner(it)+1.d0
                  endif
                  if(rjpsi.lt.10.d0) then                      
                     irjpsi=int(rjpsi/0.1)+1
                     if(pjpsi.lt.5.d0) then
                        ipjpsi=int(pjpsi/0.1)+1                
                        if(pT.le.5.d0) then
                           PSdensnative(it,irjpsi,ipjpsi)=
     &                          PSdensnative(it,irjpsi,ipjpsi)+1.d0
                        else
                           PSdensnativehighpT(it,irjpsi,ipjpsi)=
     &                         PSdensnativehighpT(it,irjpsi,ipjpsi)+1.d0
                        endif
                     endif
                  endif
               endif
            endif
            ipart2=partinfo_next(ipart2)
         enddo
         ipart=partinfo_next(ipart)
      enddo
      end


      subroutine process_events(tinf,tsup)
      use PBGgenvar
      implicit none
      logical transf
      double precision tinf,tsup,x_event(0:3),tb_event,timeievent,
     &  delt,deltr2max,tb_newevent
      integer nnewevents,firstnewevent,i,j,test,scatcode,ii,nevents,
     &  nextevent,ievent,kf(0:4),oldfirstevent,ipart,prevpart,nextpart,
     &  kk,kkk,ipartcheck,newscatcode,inewevent,lastevent,
     &  oldfirstemptyevent,eventsprocessed,ipart2,veto,allveto(0:15)
      character*30 vetoname(0:15)
      data vetoname/'unknown','one part not in plasma',
     & 'no react for these part types','delta y too large',
     & 'delta r too large','in the past of one part','Bjork time^2 <0',
     & 'Bjorken time smaller tB inf','Bjorken time larger tB sup',
     & 'too far apart at closest','too hot around','too cold around',
     & 'pi b^2>sigma(s,T,...)',' ',' ',' '/
      integer ifmtx
      call getMonitorFileIndex(ifmtx)
      deltr2max=4*(tsup-tinf)**2
      do i=0,15
         allveto(i)=0
      enddo
      nnewevents=0
      firstnewevent=-1
      ipart=firstpart
      do while (ipart.NE.-1)  
         ipart2=ipart
         do while (ipart2.NE.-1)
            if(ipart.eq.ipart2) then
               call dissocfind(ipart,tsup,test,scatcode,
     &                     x_event,tb_event,delt)
            else
               call colfindb(ipart,ipart2,tinf,tsup,test,
     &                scatcode,x_event,tb_event,delt,veto)
               if(test.eq.0) allveto(veto)=allveto(veto)+1
            endif
            if(test.eq.1) then
               call addevent(firstnewevent,x_event,tb_event,delt,
     &                 scatcode,ipart,ipart2)
               if(infotrack.ge.3) write(6,*) 'Event added: ',scatcode
               nnewevents=nnewevents+1
               n_newevent(nnewevents)=firstnewevent
               if(itypeorder.eq.1) then
                  time_newevent(nnewevents)=x_event(0)          
               else
                  time_newevent(nnewevents)=tb_event          
               endif   
            endif 
            ipart2=partinfo_next(ipart2)
         enddo
         ipart=partinfo_next(ipart)
      enddo
      if(infotrack.ge.2) write(6,*) nnewevents,'original events created'
      if(nnewevents.eq.0) then 
        return
      endif 
      if(nnewevents.gt.1) then 
        call shell_sort(nnewevents,time_newevent,n_newevent)
      endif
      ii=1
      firstevent=n_newevent(ii)
      ievent=firstevent
      nextevent=n_newevent(ii)
      do ii=2,nnewevents
       nextevent=n_newevent(ii)
       eventinfo_next(ievent)=nextevent
       ievent=nextevent
      enddo
      eventinfo_next(nextevent)=-1
      nevents=nnewevents
      if(nevents.gt.neventmax) then
       write(ifmtx,*) 'increase neventmax parameter'
       close(ifmtx)
       stop
      endif
      if(infotrack.ge.3) write(6,*)'Initialization completed with',
     &  nnewevents,' future events'
      eventsprocessed=0
 10   if(nevents.ne.0) then
         if(infotrack.ge.4) write(6,*) 'still ',nevents,' to process'
         if(infotrack.ge.4) call event_chain_test()
         nevents=nevents-1
         i=eventinfo_i(firstevent)
         j=eventinfo_j(firstevent)
         if(infotrack.ge.3) then
           write(6,*) 'processing event at ordering time:',
     &       eventinfo_tord(firstevent),i,j
         endif
         tb_event=eventinfo_tb(firstevent)
         scatcode=eventinfo_scatcode(firstevent)
         if(scatcode.gt.0) then
            call resonancemaker(firstevent,kf)
            if(infotrack.ge.3) then
             write(6,*) 'One more J/psi:',resinfo_number(3)
            endif
         endif
         if(scatcode.lt.0) then
          call resonancedissoc(firstevent,kf)
         endif
         if(scatcode.eq.0) then
            write(ifmtx,*) 'scatcode not defined'
            close(ifmtx)
            stop
         endif
         oldfirstevent=firstevent
         firstevent=eventinfo_next(firstevent)
         firstnewevent=-1
         if(firstemptyevent.eq.neventmax) then
            write(ifmtx,*) 'neventmax is too small'
            close(ifmtx)
            stop
         endif
         eventinfo_next(oldfirstevent)=firstemptyevent
         firstemptyevent=oldfirstevent

         nnewevents=0
         ipart=firstpart
         prevpart=-1
 20      If(ipart.eq.-1) goto 30
         nextpart=partinfo_next(ipart)
         if((scatcode.NE.0).AND.((i.EQ.ipart).OR.(j.EQ.ipart))) then
            if(prevpart.ne.-1) then
               partinfo_next(prevpart)=partinfo_next(ipart)
            else
               firstpart=partinfo_next(ipart)
            endif
            if(firstemptypart.eq.npartmax) then
               write(ifmtx,*) 'npartmax is too small'
               close(ifmtx)
               stop
            endif
            partinfo_next(ipart)=firstemptypart
            firstemptypart=ipart
            npart=npart-1
         else
            do kk=1,kf(0)
               test=0
               ipartcheck=1
               do kkk=1,kf(0)
                  if(kf(kkk).eq.ipart) ipartcheck=0
               enddo
               if(ipartcheck.eq.1) then
                  call colfindb(ipart,kf(kk),tb_event,tsup,
     &              test,newscatcode,x_event,tb_newevent,delt,veto)
               endif
               if(ipart.eq.kf(kk)) then
                  call dissocfind(kf(kk),tsup,test,newscatcode,
     &             x_event,tb_newevent,delt)
               endif
               if(test.eq.1) then
                  if(infotrack.ge.4) write(6,*) 'new EV:',ipart,
     &               kf(kk),newscatcode,x_event(0)
                  nnewevents=nnewevents+1
                  if(tb_newevent.le.tb_event) then
                     write(ifmtx,*) 'BIG PROBLEM IN PROCESS EVENT'
                     write(ifmtx,*) 'WE ARE GOING BACKWARDS'
                     write(ifmtx,*) tinf,tsup,tb_event,tb_newevent,
     &                          x_event(0)
                     write(ifmtx,*) 'part:',ipart,kf(kk)
                     write(ifmtx,*) 'event type:',newscatcode
                     write(ifmtx,*) partinfo_r(ipart,0),
     &                    partinfo_r(ipart,1),partinfo_r(ipart,2),
     &                    partinfo_r(ipart,3)
                    write(ifmtx,*)  partinfo_p(ipart,0),
     &                    partinfo_p(ipart,1),partinfo_p(ipart,2),
     &                    partinfo_p(ipart,3)
                     close(ifmtx)
                     stop
                  endif   
                  call addevent(firstnewevent,x_event,tb_newevent,delt,
     &                 newscatcode,kf(kk),ipart)
                  n_newevent(nnewevents)=firstnewevent
                  if(itypeorder.eq.1) then
                     time_newevent(nnewevents)=x_event(0)          
                  else
                     time_newevent(nnewevents)=tb_newevent          
                  endif   
               endif
            enddo
            prevpart=ipart
         endif
         ipart=nextpart
         goto 20   
 30      continue
         if((nnewevents.gt.0).and.(infotrack.ge.3)) then
            write(6,*) nnewevents,' NEW EVENTS FOUND:',firstnewevent 
         endif   
         if(nnewevents.gt.0)then
            if(nnewevents.gt.1) then
               call shell_sort(nnewevents,time_newevent,n_newevent)
            endif
            ii=1
            firstnewevent=n_newevent(ii)
            inewevent=firstnewevent
            nextevent=n_newevent(ii)
            do ii=2,nnewevents
               nextevent=n_newevent(ii)
               eventinfo_next(inewevent)=nextevent
               inewevent=nextevent
            enddo
            eventinfo_next(nextevent)=-1
         endif
         nevents=nevents+nnewevents
         if(nevents.ge.neventmax) then
            write(ifmtx,*) 'Make neventmax larger!'
            close(ifmtx)
            stop
         endif
         ievent=firstevent
         lastevent=-1
         ii=1
         inewevent=n_newevent(ii)
         do while((ievent.NE.-1).or.(nnewevents.GT.0))         
 35         if((ievent.NE.-1).AND.((eventinfo_i(ievent).EQ.i).OR.
     &           (eventinfo_j(ievent).EQ.j)
     &           .OR.(eventinfo_i(ievent).EQ.j).OR.
     &           (eventinfo_j(ievent).EQ.i))) then
               if(firstemptyevent.EQ.neventmax) then
                  write(ifmtx,*) 'neventmax is too small'
                  close(ifmtx)
                  stop
               endif
               oldfirstemptyevent=firstemptyevent
               firstemptyevent=ievent
               if(lastevent.NE.-1) then
                  eventinfo_next(lastevent)=eventinfo_next(ievent)
               else
                  firstevent=eventinfo_next(ievent)
               endif
               nevents=nevents-1
               ievent=eventinfo_next(ievent)
               eventinfo_next(firstemptyevent)=oldfirstemptyevent
               goto 35
            endif
 40         transf=.false.
            if(nnewevents.GT.0) then
               if(ievent.EQ.-1) then
                  transf=.true.
               else
                  if(itypeorder.eq.1) then
                     timeievent=eventinfo_r(ievent,0)
                  else
                     timeievent=eventinfo_tb(ievent)
                  endif
                  transf=(time_newevent(ii).LT.timeievent)
               endif
            endif   
            if(transf) then
               if(lastevent.GE.0) then
                  eventinfo_next(lastevent)=inewevent
               else
                  firstevent=inewevent
               endif
               lastevent=inewevent
               eventinfo_next(inewevent)=ievent
               firstnewevent=n_newevent(ii+1)
               ii=ii+1
               nnewevents=nnewevents-1
               inewevent=n_newevent(ii)
               goto 40
            endif
            if(ievent.NE.-1) then
               lastevent=ievent
               ievent=eventinfo_next(ievent)
            endif
         enddo 
         eventsprocessed=eventsprocessed+1
         goto 10
      ENDIF
      return
      end


      subroutine dissocfind(ipart,tordfin,test,scatcode,x_event,
     & tb_event,delt)
      use PBGgenvar
      use PBGsceplas
      use PBGevolQ
      use PBGpawc
      implicit none
      integer test,ires,ipart,i,scatcode,ifail,imed,veto
      double precision tordfin,x0(0:3),p(0:3),tfin,gtimefin2,dtau,
     & xfin(0:3),tresol,tineffect,tfineffect,survivalprob,surv,
     & x_event(0:3),tb_event,ran2,delt,temp,beta(3),densener,effdeg
      integer ifmtx
      call getMonitorFileIndex(ifmtx)
      test=0
      veto=0      
      ires=partinfo_ires(ipart)
      if(resinfo_dissoc(ires).ne.1) then
         return
      endif   
      do i=0,3
         x0(i)=partinfo_r(ipart,i)
         p(i)=partinfo_p(ipart,i)
      enddo
      if(itypeorder.eq.1) then
         tfin=tordfin
      else
         tfin=gtimefin2(ipart,tordfin)
      endif   
      if(tfin.lt.x0(0)) then
         write(ifmtx,*) 'in dissocfind; ',tfin,' should be >',x0(0)
         write(ifmtx,*) ipart,p(0),p(1),p(2),p(3),resinfo_mass(ires)
         write(ifmtx,*) p(0)**2-p(1)**2-p(2)**2-p(3)**2-
     &        resinfo_mass(ires)**2
         write(ifmtx,*) x0(0),x0(3),tfin
         close(ifmtx)
         stop
      endif
      dtau=(tfin-x0(0))/p(0)
      do i=1,3
         xfin(i)=x0(i)+p(i)*dtau
      enddo
      xfin(0)=tfin
      tresol=0.05*(tfin-x0(0))
      call givewindow(x0,xfin,temptrans,1,tresol,ifail,
     &                tineffect,tfineffect)
      if(ifail.ne.0) return
      if(tineffect.lt.x0(0)) then
         write(ifmtx,*) 'stupidity in dissoc_find; tineffect<tin'
         write(ifmtx,*) tineffect,x0(0),tresol
         close(ifmtx)
         stop
      endif
      if(tfineffect.gt.tfin) then
         write(ifmtx,*) 'stupidity in dissoc_find; tfineffect>tfin'
         write(ifmtx,*) tfineffect,tfin,tresol
         close(ifmtx)
         stop
      endif
      if(tfineffect.lt.tineffect) then
         write(ifmtx,*) 'tfineffect<tineffect'
         write(ifmtx,*) x0(0),tfin,tineffect,tfineffect,tresol
         close(ifmtx)
         stop
      endif
      if((ires.eq.4).or.(ires.eq.5)) then
         write(ifmtx,*) 'In dissocfind: D or Dbar will reinteract w QGP'
         write(ifmtx,*) 'during [',tineffect,',',tfineffect,'] of 
     &    [',x0(0),',',tfin,']'
         write(ifmtx,*) 'case not considered up to now; do something'
         close(ifmtx)
         stop
      endif
      surv=survivalprob(x0,p,tineffect,tfineffect)         
      call givelocalhydro(x0,temp,beta,densener,ifail,imed)
      if(isce<7) then 
      if((imed.eq.1).and.ifdecrease) then
         effdeg=min(1.0d0,densener/entrans_max)
      else
         effdeg=1.d0
      endif
      end if
      write(6,*) 'part',ipart,' of type',ires,
     &           ' survives with a prob',surv
      if(ran2().lt.effdeg*(1-surv)) then
         test=1
         scatcode=-ires
         dtau=(ran2()*(tfineffect-tineffect)+tineffect-x0(0))/p(0)
         if(dtau.lt.0.) then
            write(ifmtx,*) 'tfineffect<tineffect'
            write(ifmtx,*) x0(0),tfin,tineffect,tfineffect,tresol,dtau
            close(ifmtx)
            stop
         endif
         do i=0,3
            x_event(i)=x0(i)+dtau*p(i)
         enddo   
         tb_event=sqrt(abs(x_event(0)*x_event(0)-x_event(3)*x_event(3)))
         delt=0.0
      endif   
      return
      end


      subroutine colfindb(ip1,ip2,tordin,tordfin,test,scatcode,x_event,
     &  tb_event,delta_t,veto)
      use PBGblcktransi
      use PBGgenvar
      use PBGpartinmedium
      use PBGbookpsi
      use PBGpawc
      implicit none
      integer ip1,ip2,k,ires1,ires2,test,scatcode,ifail,imed,veto
      double precision tordin,tordfin,tfin2,tb_event,x_event(0:3),m1,
     & m2,a,b,srap_event,pi
      double precision p1(0:3),p2(0:3),x1(0:3),x2(0:3),s,
     & xrel(0:3),m1inv,m2inv,u1(0:3),u2(0:3),u1u2,u1dx,u2dx,dtau1,dtau2,
     & den,xrelcoll(0:3),sigfus,sigmafus,delta_t,drelcm2,    
     & temp,betaplasma(0:3),densener,x1ip,x1im,x2ip,x2im
      parameter (pi=3.14159265358979323844d0)
      scatcode=0
      test=0
      veto=0
      if((partinfo_whichmedium(ip1).ne.0).or.
     &   (partinfo_whichmedium(ip2).ne.0)) then
         veto=1
         goto 123
      endif   
      ires1=partinfo_ires(ip1)
      ires2=partinfo_ires(ip2)
      if(ifreac(ires1,ires2).eq.0) then
         veto=2 
         goto123
      endif   
      if (infotrack.eq.3) then
        write(6,*) 'entering colfind ',ip1,ip2
        write(6,*) partinfo_p(ip1,0),partinfo_p(ip1,1),
     &   partinfo_p(ip1,2),partinfo_p(ip1,3)
        write(6,*) partinfo_p(ip2,0),partinfo_p(ip2,1),
     &   partinfo_p(ip2,2),partinfo_p(ip2,3)
      endif
      do k=0,3
        x1(k)=partinfo_r(ip1,k)
        x2(k)=partinfo_r(ip2,k)
      enddo
      if(itypeorder.eq.1) then
         if((x1(1)-x2(1))**2+(x1(2)-x2(2))**2+(x1(3)-x2(3))**2.gt.
     &       (2*tordfin-x1(0)-x2(0))**2+biggestb2(ires1,ires2)) then
            veto=3
            goto 123
         endif
      endif
      if(itypeorder.eq.2) then
         tfin2=tordfin*tordfin
         x1im=x1(0)-x1(3)
         x2ip=x2(0)+x2(3)
         a=x1im*x2ip/tfin2
         if(a.ge.1) then
            if ((a-1)**2/a.gt.biggestb2(ires1,ires2)/tfin2) then
                veto=3
               goto 123
            endif   
         else
            x2im=x2(0)-x2(3)
            x1ip=x1(0)+x1(3)
            b=x2im*x1ip/tfin2
            if(b.ge.1) then
               if ((b-1)**2/b.gt.biggestb2(ires1,ires2)/tfin2) then
                  veto=3
                  goto 123
               endif   
            endif   
         endif   
         if((x1(1)-x2(1))**2+(x1(2)-x2(2))**2.gt.(2*tordfin
     &        -partinfo_tb(ip1)-
     &        partinfo_tb(ip2))**2+biggestb2(ires1,ires2)) then
            veto=4
            goto 123
         endif   
      endif   
      do k=0,3
         p1(k)=partinfo_p(ip1,k)
         p2(k)=partinfo_p(ip2,k)
         xrel(k)=x2(k)-x1(k)
      enddo
      m1=resinfo_mass(ires1)
      m2=resinfo_mass(ires2)
      m1inv=1./m1
      m2inv=1./m2
      u1(0)=p1(0)*m1inv
      u2(0)=p2(0)*m2inv
      u1u2=u1(0)*u2(0)
      u1dx=u1(0)*xrel(0)
      u2dx=u2(0)*xrel(0)
      do k=1,3
         u1(k)=p1(k)*m1inv
         u2(k)=p2(k)*m2inv
         u1dx=u1dx-u1(k)*xrel(k)
         u2dx=u2dx-u2(k)*xrel(k)
         u1u2=u1u2-u1(k)*u2(k)
      enddo
      den=1./(1-u1u2**2)
      dtau1=den*(u1dx-u1u2*u2dx)
      if(dtau1.lt.-0.5) then
         veto=5
         goto 123
      endif   
      dtau2=den*(u1u2*u1dx-u2dx)
      if(dtau2.lt.-0.5) then
         veto=5
         goto 123
      endif
      x_event(0)=0.5*(x1(0)+x2(0)+u1(0)*dtau1+u2(0)*dtau2)
      x_event(3)=0.5*(x1(3)+x2(3)+u1(3)*dtau1+u2(3)*dtau2)
      tb_event=x_event(0)**2-x_event(3)**2
      if((tb_event.lt.0.d0).or.(x_event(0).le.0)) then
         veto=6
         goto 123
      endif   
      tb_event=sqrt(tb_event)
      if(itypeorder.eq.1) then
         if(x_event(0).lt.tordin) then
            veto=7
            goto 123
          else if(x_event(0).gt.tordfin) then
             veto=8
             goto 123
          endif   
      endif
      if(itypeorder.eq.2) then
         if(tb_event.lt.tordin) then
            veto=7
            goto 123
         else if(tb_event.gt.tordfin) then
            veto=8
            goto123
         endif   
      endif   
      do k=0,3
         xrelcoll(k)=xrel(k)+u2(k)*dtau2-u1(k)*dtau1
      enddo   
      drelcm2=xrelcoll(1)**2+xrelcoll(2)**2+xrelcoll(3)**2
     &       -xrelcoll(0)**2
      if(drelcm2.gt.biggestb2(ires1,ires2)) then
         veto=9
         goto 123
      endif   
      s=m1**2+m2**2+2*(p1(0)*p2(0)-p1(1)*p2(1)-p1(2)*p2(2)-p1(3)*p2(3))
      x_event(1)=0.5*(x1(1)+x2(1)+u1(1)*dtau1+u2(1)*dtau2)
      x_event(2)=0.5*(x1(2)+x2(2)+u1(2)*dtau1+u2(2)*dtau2)
      scatcode=3
      call givelocalhydro(x_event,temp,betaplasma,densener,ifail,imed)
      if((temp.gt.tempdissoc(scatcode)).or.(ifail.eq.-1)) then
         veto=10
         goto 123
      else
         if(densener.lt.entrans_min) then
            veto=11
            goto 123
         endif   
      endif
      sigfus=sigmafus(ires1,ires2,s)
      if((drelcm2.lt.0).or.(drelcm2.gt.(sigfus/pi))) then
         veto=12
         goto 123
      endif   
      if(partinfo_mother(ip1).eq.partinfo_mother(ip2)) then
         nbdiag=nbdiag+1
      else
         nboffdiag=nboffdiag+1
      endif
      delta_t=xrelcoll(0)
      test=1
      srap_event=0.5*log((x_event(0)+x_event(3))/(x_event(0)
     &         -x_event(3)))
 123  if (infotrack.eq.2) write(6,*) 'exiting colfindb'
      return
      end


      subroutine resonancemaker(ievent,kf)
      use PBGgenvar
      implicit none
      integer ievent,kf(0:3),alpha,scatcode,i,j,ires1,ires2
      double precision delt,x_event(0:3),pfin1(0:3),
     &  pfin2(0:3),ptot(0:3),prel(0:3),min1,min2,mfin1,mfin2
      integer ifmtx
      call getMonitorFileIndex(ifmtx)
      i=eventinfo_i(ievent)
      j=eventinfo_j(ievent)
      do alpha=0,3
        x_event(alpha)=eventinfo_r(ievent,alpha)
      enddo
      scatcode=eventinfo_scatcode(ievent)
      delt=eventinfo_delt(ievent)
      kf(0)=1
      kf(1)=firstemptypart
      firstemptypart=partinfo_next(firstemptypart)
      partinfo_next(kf(1))=firstpart
      partinfo_ires(kf(1))=scatcode
      partinfo_mother(kf(1))=0
      firstpart=kf(1)
      partinfo_mass(kf(1))=resinfo_mass(partinfo_ires(kf(1)))
      partinfo_whichmedium(kf(1))=0
      do alpha=0,3
         partinfo_r(kf(1),alpha)=x_event(alpha)
         ptot(alpha)=partinfo_p(i,alpha)+partinfo_p(j,alpha)
         prel(alpha)=0.5*(partinfo_p(j,alpha)-partinfo_p(i,alpha))
      enddo
      partinfo_tb(kf(1))=eventinfo_tb(ievent)
      min1=partinfo_mass(i)
      min2=partinfo_mass(j)
      mfin1=partinfo_mass(kf(1))
      mfin2=0.0
      ires1=partinfo_ires(i)
      ires2=partinfo_ires(j)
      if(ires1*ires2.eq.2) then
         call fusion_diff(ires1,ires2,ptot,prel,min1,min2,mfin1,mfin2,
     &                    pfin1,pfin2)
      else
         write(ifmtx,*) 'confirm that you still want to use isotopic'
         close(ifmtx)
         stop
      endif   
      if((abs(pfin1(0)**2-pfin1(1)**2-pfin1(2)**2-pfin1(3)**2
     & -mfin1**2).gt.1.D-7).or.
     &   (abs(pfin2(0)**2-pfin2(1)**2-pfin2(2)**2-pfin2(3)**2
     & -mfin2**2).gt.1.D-7)) then
         write(ifmtx,*) 'big prob in resonance maker'
         write(ifmtx,*) pfin1(0),pfin1(1),pfin1(2),pfin1(3),mfin1
         write(ifmtx,*) pfin2(0),pfin2(1),pfin2(2),pfin2(3),mfin2
         close(ifmtx)
         stop
      endif   
      do alpha=0,3
         partinfo_p(kf(1),alpha)=pfin1(alpha)
      enddo
      npart=npart+1
      resinfo_number(scatcode)=resinfo_number(scatcode)+1
      return
      end


      subroutine resonancedissoc(ievent,kf)
      use PBGgenvar
      implicit none
      integer ievent,kf(0:3),ipart,ires(0:3),alpha,ifail,imed,k,
     &    nbodies,ibody
      parameter (nbodies=2)
      double precision pres(0:3),x(0:3),tempevent,
     & betahydro(3),densener,mass,pgluon(0:3),s,ss,ptot(0:3),
     & prelin(0:3),
     & ucm(0:3),prelincm(0:3),sr2,cthcm,costhetafus,vprelincm(3),
     & vprelfincm(3),pf,
     & prelfincm(0:3),prelfin(0:3),pfin(3,0:3),chk(3)
      integer ifmtx
      call getMonitorFileIndex(ifmtx)
      ipart=eventinfo_i(ievent)
      ires(0)=-eventinfo_scatcode(ievent)
      if(ires(0).ne.3) then
         write(ifmtx,*) 'resonance decay is just for J/Psi'
         close(ifmtx)
         stop
      endif
      mass=partinfo_mass(ipart)
      do alpha=0,3
         pres(alpha)=partinfo_p(ipart,alpha)
         x(alpha)=eventinfo_r(ievent,alpha)   
      enddo
      call givelocalhydro(x,tempevent,betahydro,densener,ifail,imed)
      if(ifail.ne.0) then
         write(ifmtx,*) 'dissoc resonance: temp should be > 0 !!!'
         write(ifmtx,*) 'ifail is :',ifail
         write(ifmtx,*) 'temperature is :', tempevent
         close(ifmtx)
         stop
      endif   
      call flashpsi(pres,tempevent,betahydro,pgluon,s)
      ss=sqrt(s)
      do alpha=0,3
         prelin(alpha)=pres(alpha)-pgluon(alpha)
         ptot(alpha)=pres(alpha)+pgluon(alpha)
         ucm(alpha)=ptot(alpha)/ss
      enddo   
      call lorentz2(ucm,prelin,prelincm)
      sr2=(s-4*mc2)/(4*mc2-mpsi**2)
      if(sr2.lt.0.) then
         write(ifmtx,*) 'error 1 in resonance dissoc'
         close(ifmtx)
         stop
      endif         
      cthcm=costhetafus(sr2)
      vprelincm(1)=prelincm(1)
      vprelincm(2)=prelincm(2)
      vprelincm(3)=prelincm(3)
      pf=sqrt(s-4*mc2)/2
      call rotate3d(vprelincm,cthcm,pf,vprelfincm)
      prelfincm(0)=0.
      do k=1,3
         prelfincm(k)=vprelfincm(k)
         ucm(k)=-ucm(k)
      enddo
      call lorentz2(ucm,prelfincm,prelfin)
      do alpha=0,3
         pfin(1,alpha)=ptot(alpha)/2+prelfin(alpha)
         pfin(2,alpha)=ptot(alpha)/2-prelfin(alpha)
      enddo   
      chk(1)=pfin(1,0)**2-(pfin(1,1)**2+pfin(1,2)**2+pfin(1,3)**2+mc2)
      chk(2)=pfin(2,0)**2-(pfin(2,1)**2+pfin(2,2)**2+pfin(2,3)**2+mc2)
      if((abs(chk(1)).gt.1d-7).or.(abs(chk(2)).gt.1d-7)) then
         write(ifmtx,*) 'in resonance_dissoc; c or cbar not on shell'
         write(ifmtx,*) 'c:',chk(1),'; cbar:',chk(2) 
         close(ifmtx)
         stop
      endif   
      ires(1)=1
      ires(2)=2
      npart=npart+nbodies
      kf(0)=nbodies
      nbccbaremission=nbccbaremission+1
      do ibody=1,nbodies
        kf(ibody)=firstemptypart
        firstemptypart=partinfo_next(firstemptypart)
        partinfo_next(kf(ibody))=firstpart
        partinfo_ires(kf(ibody))=ires(ibody)
        partinfo_mother(kf(ibody))=nbccbaremission
        firstpart=kf(ibody)
        do alpha=0,3
           partinfo_p(kf(ibody),alpha)=pfin(ibody,alpha)
           partinfo_r(kf(ibody),alpha)=eventinfo_r(ievent,alpha)
        enddo   
        partinfo_tb(kf(ibody))=eventinfo_tb(ievent)
        partinfo_mass(kf(ibody))=resinfo_mass(ires(ibody))
        partinfo_whichmedium(kf(ibody))=0
      enddo
      resinfo_number(ires(0))=resinfo_number(ires(0))-1
      return
      end


      subroutine flashpsi(ppsi,temp,betahydro,pgluon,s)
      implicit none
      integer k
      double precision ppsi(0:3),temp,betahydro(3),pgluon(0:3),s,chck,
     & uhydro(0:4),ppsih(0:3),pg,ctheta,vpsi(3),vpg(3),pgluonh(0:3)
      integer ifmtx
      call getMonitorFileIndex(ifmtx)
      uhydro(0)=1/sqrt(1-(betahydro(1)**2+betahydro(2)**2
     &         +betahydro(3)**2))
      uhydro(1)=uhydro(0)*betahydro(1)
      uhydro(2)=uhydro(0)*betahydro(2)
      uhydro(3)=uhydro(0)*betahydro(3)
      call lorentz2(uhydro,ppsi,ppsih)
      call givesdiss(temp,ppsih(0),s,pg,ctheta)
      vpsi(1)=ppsih(1)
      vpsi(2)=ppsih(2)
      vpsi(3)=ppsih(3)
      call rotate3d(vpsi,ctheta,pg,vpg)
      pgluonh(0)=pg
      do k=1,3
         pgluonh(k)=vpg(k)
         uhydro(k)=-uhydro(k)
      enddo
      call lorentz2(uhydro,pgluonh,pgluon)
      chck=(pgluon(0)+ppsi(0))**2      
      do k=1,3
         chck=chck-(pgluon(k)+ppsi(k))**2
      enddo   
      if(abs(chck-s).gt.1d-7) then
         write(ifmtx,*) 'inv mass error in flashpsi;',chck,' vs. ',s
         close(ifmtx)
         stop
      endif
      return
      end


      subroutine testflashpsi(nbtest)
      use PBGpsiinfo
      implicit none
      integer k,nbtest,itest
      double precision betan,beta(3),pi,cbeta,sbeta,ran2,phibeta,
     & ppsin,cthpsi,sthpsi,phipsi,ppsi(0:3),temp,pgluon(0:3),s,
     & ptot(0:3)
      parameter (pi=3.14159265358979323844d0)      
      write(6,*) 'entering testflashpsi'
      do itest=1,nbtest
         betan=ran2()
         cbeta=1-2*ran2()
         sbeta=sqrt(1-cbeta**2)
         phibeta=2*ran2()*pi
         beta(3)=betan*cbeta
         beta(1)=betan*sbeta*cos(phibeta)
         beta(2)=betan*sbeta*sin(phibeta) 
         ppsin=10*ran2()
         ppsi(0)=sqrt(ppsin**2+mpsi**2)
         cthpsi=1-2*ran2()
         sthpsi=sqrt(1-cthpsi**2)
         phipsi=2*ran2()*pi
         ppsi(3)=ppsin*cthpsi
         ppsi(1)=ppsin*sthpsi*cos(phipsi)
         ppsi(2)=ppsin*sthpsi*sin(phipsi)
         temp=0.2+0.8*ran2()
         call flashpsi(ppsi,temp,beta,pgluon,s)
         do k=0,3
            ptot(k)=ppsi(k)+pgluon(k)
         enddo
      enddo   
      write(6,*) 'testflashpsi performed with succes'
      return
      end
      

      subroutine givesdiss(temp,epsi,s,pg,ctheta)
      use PBGpsiinfo
      implicit none
      integer nit,i
      parameter(nit=30)
      double precision temp,epsi,ppsi,thresh,beta,sig2,mpsi2,ds,s,stry,
     & w,wtry,ran2,cos2rand,integrantdissocbose,expinf,expsup,pg,ctheta
      integer ifmtx
      call getMonitorFileIndex(ifmtx)
      mpsi2=mpsi**2
      ppsi=sqrt(epsi**2-mpsi2)
      thresh=4*mc2
      if(ppsi.lt.temp) then
         beta=3
      else
         beta=4
      endif
      sig2=mpsi2/(mpsi2/(temp*(epsi+ppsi))+2*beta*(thresh/mpsi2-1.))
      ds=-2*sig2*(log(ran2())+log(ran2())+log(ran2())*cos2rand())
      s=ds+thresh
      w=integrantdissocbose(temp,ppsi,mpsi,s)/(exp(-ds/(2*sig2))
     &  *sqrt(ds)*ds)
      do 10 i=1,nit
         ds=-2*sig2*(log(ran2())+log(ran2())+log(ran2())*cos2rand())
         stry=ds+thresh
         wtry=integrantdissocbose(temp,ppsi,mpsi,stry)/
     &            (exp(-ds/(2*sig2))*sqrt(ds)*ds)
         if(wtry.lt.w)then
            if(w*ran2().gt.wtry) goto 10
         endif
         s=stry
         w=wtry
 10   continue
      ds=0.5*(s/mpsi2-1.)/temp
      expinf=1-exp(-ds*(epsi+ppsi))
      expsup=1-exp(-ds*mpsi2/(epsi+ppsi))
      pg=-temp*log(1-expinf*(expsup/expinf)**ran2())
      ctheta=(epsi-(s-mpsi2)/(2*pg))/ppsi
      if(abs(ctheta).gt.1) then
         write(ifmtx,*) 'problem in givesdiss: cos(theta)=',ctheta,' !'
         close(ifmtx)
         stop
      endif
      return
      end

