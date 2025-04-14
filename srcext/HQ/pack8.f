
      subroutine count_part_type(nbtypeires)
      use PBGgenvar
      implicit none 
      integer nbtypeires(nres),ires,ipart
      do ires=1,nres
        nbtypeires(ires)=0
      enddo
      ipart=firstpart
      do while (ipart.NE.-1)
         nbtypeires(partinfo_ires(ipart))=
     &     nbtypeires(partinfo_ires(ipart))+1
        ipart=partinfo_next(ipart)
      enddo
      if(infotrack.ge.1) write(6,*) (nbtypeires(ires),ires=1,nres)
      return
      end


      subroutine initpacking(nbcol,timestep)
      use PBGsystem
      use PBGgenvar
      use PBGfordisplay
      implicit none
      integer nbcol,j,k,nhist,ndiv2d,ndiv1d,ires,iclass,m
      double precision boundpxinf,boundpxsup,pi,
     & steppx,boundzinf,boundzsup,boundrtsup,timestep
      integer plasmascen
      integer ifmtx

      parameter (pi=3.14159265358979323844d0)
      call getMonitorFileIndex(ifmtx)
      if(ifinput) then
         read(515,*) rapsupobs,ycentname
      else
         rapsupobs=1.d0
         ycentname='y0'
      endif
      do k=0,ndivpsimax
         packnbpsi(k)=0.d0
      enddo
      do k=1,ndivpt1dmax
        packpty0charm(k)=0.d0
        packpty0dmeson(k)=0.d0
        packvacpty0charm(k)=0.d0
        packvacpty0dmeson(k)=0.d0
      enddo
      do ires=0,nres+4
         do k=1,ndivrap1dmax
            packrap(ires,k)=0.d0
            packmeanptsq(ires,k)=0.d0
            packvacrap(ires,k)=0.d0
            packvacmeanptsq(ires,k)=0.d0
            packvacrap_zero(ires,k)=0.d0
            packvacmeanptsq_zero(ires,k)=0.d0
         enddo
         do k=1,ndivpt1dmax
            packpty0(ires,k)=0.d0
            packpty0_0_20(ires,k)=0.d0
            packpty0_20_40(ires,k)=0.d0
            packpty0_40_60(ires,k)=0.d0
            packpty0_60_100(ires,k)=0.d0
            packpty0alongb(ires,k)=0.d0
            packpty0bis(ires,k)=0.d0
            packpty0bis_0_20(ires,k)=0.d0
            packpty0bis_20_40(ires,k)=0.d0
            packpty0bis_40_60(ires,k)=0.d0
            packpty0bis_60_100(ires,k)=0.d0
            packvacpty0(ires,k)=0.d0
            packvacpty0bis(ires,k)=0.d0
            packvacpty0_zero(ires,k)=0.d0
            packvacpty0bis_zero(ires,k)=0.d0
            packpty2(ires,k)=0.d0
            packpty2bis(ires,k)=0.d0
            packvacpty2(ires,k)=0.d0
            packvacpty2bis(ires,k)=0.d0
            packvacpty2_zero(ires,k)=0.d0
            packvacpty2bis_zero(ires,k)=0.d0
         enddo
         do k=1,ndivphi1dmax
            packphi(ires,k)=0.d0
            do iclass=0,6
               packcorrelphi(ires,iclass,k)=0.d0
            enddo                
         enddo
         do j=1,ndivpt1dmax
            do k=1,ndivphi1dmax
               packptphiy0(ires,j,k)=0.d0
               packptphiy0bis(ires,j,k)=0.d0
            enddo
         enddo
         do j=1,ndivpt1dmax
            do k=1,ndivpt1dmax
               packcorrelpty0(ires,j,k,1)=0.D0
               packcorrelpty0(ires,j,k,2)=0.D0
               packcorrelpty0(ires,j,k,3)=0.D0
            enddo
         enddo
         do j=1,2*ndivpt1dmax
           do k=1,ndivphi1dmax
           do iclass=1,5
           enddo
           enddo
         enddo
      enddo      
      do j=1,ndivpt1dmax
         packvacpty0nonpromptJpsi(j)=0.d0
         packpty0nonpromptJpsi(j)=0.d0
         do k=1,ndivphi1dmax
            packptphiy0nonpromptJpsi(j,k)=0.d0
         enddo
      enddo     
      incbas=1./nbcol
      ndiv2d=40
      ndiv1d=80
      hbkname='distr_fin.hbook'
      boundzinf=-60.
      boundzsup=60.
      if(abs(sqrts-200d0).lt.0.1) then
         boundrapinf=-5.
         boundrapsup=5.
      else if((abs(sqrts-5500d0).lt.0.1).or.(abs(sqrts-2760d0).lt.0.1)
     &        .or.(abs(sqrts-5020d0).lt.0.1)
     &        .or.(abs(sqrts-7000d0).lt.0.1)) then  
         boundrapinf=-8.
         boundrapsup=8.
      else
         write(ifmtx,*) 'unrecognized sqrt(s) in initpacking'
         close(ifmtx)
         stop
      endif
      ndivrap1d=40
      if(ndivrap1d.gt.ndivrap1dmax) then
         write(ifmtx,*) 'too large # of rapidity steps in initpacking'
         close(ifmtx)
         stop
      endif
      steprap1d=(boundrapsup-boundrapinf)/ndivrap1d
      incrap1d=1./steprap1d
      boundpxinf=-3.
      boundpxsup=3.
      steppx=(boundpxsup-boundpxinf)/ndiv2d

      if(abs(sqrts-200).lt.0.1) then
         boundptsup=15
         ndivpt1d=75
      else if((abs(sqrts-5500).lt.0.1).or.(abs(sqrts-2760d0).lt.0.1)
     &        .or.(abs(sqrts-5020d0).lt.0.1)
     &        .or.(abs(sqrts-7000d0).lt.0.1)) then  
         boundptsup=50
         ndivpt1d=100
      else
         write(ifmtx,*) 'unrecognized sqrt(s) in initpacking'
         close(ifmtx)
         stop
      endif
      if(ndivpt1d.gt.ndivpt1dmax) then
         write(ifmtx,*) 'too large # of pt steps in initpacking'
         close(ifmtx)
         stop
      endif
      steppt1d=boundptsup/ndivpt1d
      incpt1dy0=1/(2*steppt1d*rapsupobs)
      incpt1dy2=1./steppt1d
      boundphisup=2*pi
      ndivphi1d=40
      stepphi1d=boundphisup/ndivphi1d
      incphi1d=0.5/stepphi1d
      incphipt=incpt1dy0/stepphi1d
      boundrtsup=10.
      if(plasmascen()<7) then
        radstep=0.4d0
      else if(plasmascen()==7) then
        radstep=0.8d0
      endif
      if(ndivphi1d.gt.ndivphi1dmax) then
         write(ifmtx,*) 'too large # of phi steps in initpacking'
         close(ifmtx)
         stop
      endif
      nhist=0

         do k=1,ndivpt1dmax
            do m=1,ndivrap1dmax
               dncdptdyvac(k,m)=0.d0
            enddo
         enddo
      timestepforsave=max(timestep,0.1)
      do j=0,ndivtmax
         timesave(j)=j*timestepforsave
         do k=-5,5
            avdnbpsidyoftime(j,k)=0.d0
         enddo
         do k=1,ndivpt1dmax
            dncdptoftimey0(j,k)=0.0d0
            dncdptoftime(j,k)=0.0d0
            do m=1,ndivrap1dmax
              dncdptdyoftime(j,k,m)=0.d0
              dncdptdyvac(k,m)=0.d0
            enddo
         enddo
          do k=1,ndivr1dmax
             dncdroftimey0(j,k)=0.0d0
             dncdroftime(j,k)=0.0d0
          enddo
           do k=1,ndivr1dmax
             dptdroftimey0(j,k)=0.0d0
             dptdroftime(j,k)=0.0d0
          enddo
       enddo
      return
      end


      double precision function phi_of_p(p_x,p_y)
      implicit none
      double precision pi,p_x,p_y,psq,phi
      parameter (pi=3.14159265358979323844d0)
      psq=sqrt(p_x**2+p_y**2)
      if(psq.gt.0) then
         phi=acos(p_x/psq)
         if(p_y.lt.0) phi=2*pi-phi
      else
         phi=0.
      endif
      phi_of_p=phi
      return
      end

      integer function iphi_of_phi(phi)
      use PBGgenvar
      use PBGfordisplay
      implicit none

      integer iphi
      double precision phi 
      integer ifmtx
      call getMonitorFileIndex(ifmtx)
      iphi=int(phi/sngl(stepphi1d))+1
      if(iphi.gt.ndivphi1d) then
         if(phi.gt.(ndivphi1d*stepphi1d+0.001)) then
            write(ifmtx,*) 'phi beyond 2xpi in pack iphi of phi !!!'
            write(ifmtx,*) phi,iphi,ndivphi1d,stepphi1d
            close(ifmtx)
            stop
         else
            iphi=ndivphi1d 
         endif
      endif
      iphi_of_phi=iphi
      return
      end

      subroutine add_one_electron_vac(itypl,pl,ov)
      use PBGgenvar
      use PBGfordisplay
      implicit none
      integer irap,ipt,iphi,iphi_of_phi,itypl
      double precision rap,pt,phi,phi_of_p,pl(0:3),ov,branching(3),incr
      data branching/0.108,0.106,0.108/ 
      rap=0.5*log((pl(0)+pl(3))/(pl(0)-pl(3)))
      irap=1+int((rap-boundrapinf)/steprap1d)
      pt=sqrt(pl(1)**2+pl(2)**2)
      ipt=int(pt/steppt1d)+1
      phi=phi_of_p(pl(1),pl(2))
      iphi=iphi_of_phi(phi)
      if((irap.ge.1).and.(irap.le.ndivrap1d)) then
         incr=incrap1d*ov*branching(itypl)
         packvacrap(nres+itypl,irap)=packvacrap(nres+itypl,irap)+incr
         packvacrap(nres+4,irap)=packvacrap(nres+4,irap)+incr
         incr=incrap1d*ov*branching(itypl)*pt**2
         packvacmeanptsq(nres+itypl,irap)=
     &        packvacmeanptsq(nres+itypl,irap)+incr
         packvacmeanptsq(nres+4,irap)=packvacmeanptsq(nres+4,irap)+incr
      endif
      if(abs(rap).lt.rapsupobs) then
         if(ipt.le.ndivpt1d) then
            incr=ov*incpt1dy0*branching(itypl)
            packvacpty0(nres+itypl,ipt)=packvacpty0(nres+itypl,ipt)+incr
            packvacpty0(nres+4,ipt)=packvacpty0(nres+4,ipt)+incr
            incr=ov*incpt1dy0*branching(itypl)/pt
            packvacpty0bis(nres+itypl,ipt)=
     &           packvacpty0bis(nres+itypl,ipt)+incr                
            packvacpty0bis(nres+4,ipt)=packvacpty0bis(nres+4,ipt)+incr  
         endif   
      endif
      if((abs(rap).lt.2).and.(abs(rap).gt.1.5)) then
         if(ipt.le.ndivpt1d) then
            incr=ov*incpt1dy2*branching(itypl)
            packvacpty2(nres+itypl,ipt)=packvacpty2(nres+itypl,ipt)+incr
            packvacpty2(nres+4,ipt)=packvacpty2(nres+4,ipt)+incr
            incr=ov*incpt1dy2*branching(itypl)/pt
            packvacpty2bis(nres+itypl,ipt)=
     &           packvacpty2bis(nres+itypl,ipt)+incr                
            packvacpty2bis(nres+4,ipt)=packvacpty2bis(nres+4,ipt)+incr  
         endif   
      endif
      end

      subroutine add_one_JPsi_from_Bdecay(pPsi,ifvac)
      use PBGgenvar
      use PBGfordisplay
      implicit none
      logical ifvac
      integer irap,ipt,iphi_of_phi,iphi
      double precision pPsi(0:3),rap,pt,incr,phi,phi_of_p
      rap=0.5*log((pPsi(0)+pPsi(3))/(pPsi(0)-pPsi(3)))
      irap=1+int((rap-boundrapinf)/steprap1d)
      pt=sqrt(pPsi(1)**2+pPsi(2)**2)
      ipt=int(pt/steppt1d)+1
      phi=phi_of_p(pPsi(1),pPsi(2)) 
      iphi=iphi_of_phi(phi)
      if((irap.ge.1).and.(irap.le.ndivrap1d)) then
      endif
      if(abs(rap).lt.rapsupobs) then
         if(ipt.le.ndivpt1d) then
            incr=incpt1dy0*0.0115d0
            if(ifvac) then
             packvacpty0nonpromptJpsi(ipt)=packvacpty0nonpromptJpsi(ipt)
     &         +incr
            else
               packpty0nonpromptJpsi(ipt)=packpty0nonpromptJpsi(ipt)
     &         +incr
            endif
            incr=incphipt*0.0115d0
            packptphiy0nonpromptJpsi(ipt,iphi)=
     &                      packptphiy0nonpromptJpsi(ipt,iphi)+incr
         endif
      endif
      end

      subroutine pack_it_vacuum(ires,p)
      use PBGgenvar
      use PBGfordisplay
      use PBGcharge
      implicit none
      integer ires,irap,ipt,iphi,k,ilept,iphi_of_phi
      double precision p(0:3),rap,pt,phi,phi_of_p,pdb(0:3),pl(0:3)
     & ,p4JPsi(0:3)
      integer ifmtx
      call getMonitorFileIndex(ifmtx)
      if(p(0).lt.abs(p(3))) then
         write(ifmtx,*) 'p0<|p3| in pack_it_vacuum',p(0),p(3),ires
         close(ifmtx)
         stop
      endif
      rap=log((p(0)+p(3))/(p(0)-p(3)))/2
      irap=1+int((rap-boundrapinf)/steprap1d)
      pt=sqrt(p(1)**2+p(2)**2)
      ipt=int(pt/steppt1d)+1
      phi=phi_of_p(p(1),p(2))
      iphi=iphi_of_phi(phi)
      if((irap.ge.1).and.(irap.le.ndivrap1d)) then
         packvacrap(ires,irap)=packvacrap(ires,irap)
     &                        +incrap1d
         packvacmeanptsq(ires,irap)=packvacmeanptsq(ires,irap)
     &                             +pt**2
      endif   
      if(abs(rap).lt.rapsupobs) then
         if(ipt.le.ndivpt1d) then
            packvacpty0(ires,ipt)=packvacpty0(ires,ipt)
     &                           +incpt1dy0
            packvacpty0bis(ires,ipt)=packvacpty0bis(ires,ipt)
     &                              +incpt1dy0/pt
            if(ires.eq.1.or.ires.eq.2) then
              packvacpty0charm(ipt)=packvacpty0charm(ipt)
     &                             +incpt1dy0
            endif
            if(ires.eq.4.or.ires.eq.5) then
              packvacpty0dmeson(ipt)=packvacpty0dmeson(ipt)
     &                             +incpt1dy0
            endif            
            if(ncharged.gt.28) then
              packvacpty0bis_0_20(ires,ipt)=
     &          packvacpty0bis_0_20(ires,ipt)+incpt1dy0/pt
            endif
            if(ncharged.le.28.and.ncharged.ge.19) then
              packvacpty0bis_20_40(ires,ipt)=
     &          packvacpty0bis_20_40(ires,ipt)+incpt1dy0/pt    
            endif
            if(ncharged.lt.19.and.ncharged.ge.12) then
              packvacpty0bis_40_60(ires,ipt)=
     &          packvacpty0bis_40_60(ires,ipt)+incpt1dy0/pt
            endif
            if(ncharged.lt.12) then
              packvacpty0bis_60_100(ires,ipt)=
     &          packvacpty0bis_60_100(ires,ipt)+incpt1dy0/pt 
            endif          
         endif
      endif
      if((abs(rap).lt.2).and.(abs(rap).gt.1.5)) then
         if(ipt.le.ndivpt1d) then
            packvacpty2(ires,ipt)=packvacpty2(ires,ipt)
     &                           +incpt1dy2
            packvacpty2bis(ires,ipt)=packvacpty2bis(ires,ipt)
     &                              +incpt1dy2/pt
         endif
      endif
      if((ipt.le.ndivpt1d).and.(irap.le.ndivrap1d).and.(irap.ge.1)) then
            dncdptdyvac(ipt,irap)=dncdptdyvac(ipt,irap)
     &                        +incpt1dy0*2.0d0*incrap1d
      endif
      if((ires.eq.17).or.(ires.eq.18)) then
         call BdecayJPsi(p,p4JPsi)
         call add_one_JPsi_from_Bdecay(p4JPsi,.true.)
      endif
      if((ires.eq.4).or.(ires.eq.5).or.(ires.eq.9).or.(ires.eq.10)) then
         do k=0,3
            pdb(k)=p(k)
         enddo
         do ilept=1,10
            if(ires.le.nresc) then
               call weakdecayD(pdb,pl)
               call add_one_electron_vac(1,pl,0.1d0)
            else
               call weakdecayB(pdb,pl)
               call add_one_electron_vac(2,pl,0.1d0)
               call weakdecayBsecond(pdb,pl)
               call add_one_electron_vac(3,pl,0.1d0)
            endif
         enddo
      endif
      return
      end

      subroutine pack_it_vacuum_cronin(ires,p)
      use PBGgenvar
      use PBGfordisplay
      implicit none
      integer ires,irap,ipt,iphi,iphi_of_phi
      double precision p(0:3),rap,pt,phi,phi_of_p
      rap=log((p(0)+p(3))/(p(0)-p(3)))/2
      irap=1+int((rap-boundrapinf)/steprap1d)
      pt=sqrt(p(1)**2+p(2)**2)
      ipt=int(pt/steppt1d)+1
      phi=phi_of_p(p(1),p(2))
      iphi=iphi_of_phi(phi)
      if((irap.ge.1).and.(irap.le.ndivrap1d)) then
         packvacrap_cronin(ires,irap)=packvacrap_cronin(ires,irap)
     &                                +incrap1d
         packvacmeanptsq_cronin(ires,irap)=packvacmeanptsq_cronin(ires
     &                                    ,irap)
     &                                    +pt**2
      endif   
      if(abs(rap).lt.rapsupobs) then
         if(ipt.le.ndivpt1d) then
            packvacpty0_cronin(ires,ipt)=packvacpty0_cronin(ires,ipt)
     &                                 +incpt1dy0
            packvacpty0bis_cronin(ires,ipt)=
     &              packvacpty0bis_cronin(ires,ipt)+incpt1dy0/pt
         endif
      endif
      if((abs(rap).lt.2).and.(abs(rap).gt.1.5)) then
         if(ipt.le.ndivpt1d) then
            packvacpty2_cronin(ires,ipt)=packvacpty2_cronin(ires,ipt)
     &                                  +incpt1dy0
            packvacpty2bis_cronin(ires,ipt)=
     &   packvacpty2bis_cronin(ires,ipt)+incpt1dy2/pt
         endif
      endif
      return
      end


      subroutine pack_it_vacuum_zero(ires,p)
      use PBGgenvar
      use PBGfordisplay
      implicit none
      integer ires,irap,ipt,iphi,iphi_of_phi
      double precision p(0:3),rap,pt,phi,phi_of_p
      rap=log((p(0)+p(3))/(p(0)-p(3)))/2
      irap=1+int((rap-boundrapinf)/steprap1d)
      pt=sqrt(p(1)**2+p(2)**2)
      ipt=int(pt/steppt1d)+1
      phi=phi_of_p(p(1),p(2))
      iphi=iphi_of_phi(phi)
      if((irap.ge.1).and.(irap.le.ndivrap1d)) then
         packvacrap_zero(ires,irap)=packvacrap_zero(ires,irap)+incrap1d
         packvacmeanptsq_zero(ires,irap)=packvacmeanptsq_zero(ires,irap)
     &                                    +pt**2
      endif   
      if(abs(rap).lt.rapsupobs) then
         if(ipt.le.ndivpt1d) then
            packvacpty0_zero(ires,ipt)=packvacpty0_zero(ires,ipt)
     &                                 +incpt1dy0
            packvacpty0bis_zero(ires,ipt)=packvacpty0bis_zero(ires,ipt)
     &                                    +incpt1dy0/pt
         endif
      endif
      if((abs(rap).lt.2).and.(abs(rap).gt.1.5)) then
         if(ipt.le.ndivpt1d) then
            packvacpty2_zero(ires,ipt)=packvacpty2_zero(ires,ipt)
     &                                  +incpt1dy0
            packvacpty2bis_zero(ires,ipt)=packvacpty2bis_zero(ires,ipt)
     &                                     +incpt1dy2/pt
         endif
      endif
      return
      end


      subroutine packonetimeradius(time,timenextsave)
      use PBGgenvar
      use PBGfordisplay
      implicit none
      integer itimesave,ipart,ipt,k,irr
      integer inumberc,inumbercy0
      double precision time,timenextsave,rap,pt,p(0:3),x(0:3)
      double precision rradius
      integer incr1d
      integer ifmtx
      call getMonitorFileIndex(ifmtx)

      itimesave=int(time/timestepforsave)
      if(itimesave.gt.ndivtmax) then 
         write(ifmtx,*) 'itimesave>max foreseen in packonetime'
         close(ifmtx)
         stop
      endif
      incr1d=1.0d0/radstep/2.0
      inumbercy0=0
      inumberc=0
      ipart=firstpart
      do while (ipart.NE.-1)
         if(partinfo_ires(ipart).eq.1) then
            do k=0,3
               p(k)=partinfo_p(ipart,k)
               x(k)=partinfo_r(ipart,k)
            enddo
            rap=log((partinfo_p(ipart,0)+partinfo_p(ipart,3))/
     &           (partinfo_p(ipart,0)-partinfo_p(ipart,3)))/2
            pt=sqrt(p(1)**2+p(2)**2)
            ipt=int(pt/steppt1d)+1
            rradius=sqrt(x(1)**2+x(2)**2)
            irr=int(rradius/radstep)
            if(abs(rap).lt.rapsupobs) then
               inumbercy0=inumbercy0+1
               if(irr.le.ndivr1dmax) then
                dncdroftimey0(itimesave,ipt)=
     &           dncdroftimey0(itimesave,ipt)+incr1d
                dptdroftimey0(itimesave,ipt)=
     &           dptdroftimey0(itimesave,ipt)
     &           +pt/2.0/radstep
               endif
            endif
            if(abs(rap).lt.5.0) then
            inumberc=inumberc+1
            if(irr.le.ndivr1dmax) then
                dncdroftime(itimesave,ipt)=
     &           dncdroftime(itimesave,ipt)
     &            +incr1d/5.0
                dptdroftime(itimesave,ipt)=
     &           dptdroftime(itimesave,ipt)+pt/5.0/radstep
            endif
            endif
         endif
         ipart=partinfo_next(ipart)
      enddo
      timenextsave=timesave(itimesave+1)
      return
      end

      subroutine packonetimednc(time,timenextsave)
      use PBGgenvar
      use PBGfordisplay
      implicit none
      integer itimesave,ipart,ipt,k
      integer inumberc,inumbercy0
      double precision time,timenextsave,rap,pt,p(0:3)
      integer ifmtx
      call getMonitorFileIndex(ifmtx)

      itimesave=int(time/timestepforsave)
      if(itimesave.gt.ndivtmax) then 
         write(ifmtx,*) 'itimesave > max foreseen in packonetime'
         close(ifmtx)
         stop
      endif
      inumbercy0=0
      inumberc=0
      ipart=firstpart
      do while (ipart.NE.-1)
         if(partinfo_ires(ipart).eq.1) then
            do k=0,3
               p(k)=partinfo_p(ipart,k)
            enddo
            rap=log((partinfo_p(ipart,0)+partinfo_p(ipart,3))/
     &           (partinfo_p(ipart,0)-partinfo_p(ipart,3)))/2.0d0
            pt=sqrt(p(1)**2+p(2)**2)
            ipt=int(pt/steppt1d)+1
            if(abs(rap).lt.rapsupobs) then
               inumbercy0=inumbercy0+1
               if(ipt.le.ndivpt1d) then
                dncdptoftimey0(itimesave,ipt)=
     &           dncdptoftimey0(itimesave,ipt)+incpt1dy0
               endif
            endif
            if(abs(rap).lt.5.0) then
            inumberc=inumberc+1
            if(ipt.le.ndivpt1d) then
                dncdptoftime(itimesave,ipt)=
     &           dncdptoftime(itimesave,ipt)
     &            +incpt1dy0/5.0
            endif
            endif
         endif
         if(partinfo_ires(ipart).eq.4) then
            do k=0,3
               p(k)=partinfo_pghost(ipart,k)
            enddo
            rap=log((p(0)+p(3))/(p(0)-p(3)))/2.0d0
            pt=sqrt(p(1)**2+p(2)**2)
            ipt=int(pt/steppt1d)+1
            if(abs(rap).lt.rapsupobs) then
               if(ipt.le.ndivpt1d) then
                dncdptoftimey0(itimesave,ipt)=
     &           dncdptoftimey0(itimesave,ipt)+incpt1dy0
               endif
            endif
            if(abs(rap).lt.5.0) then
            if(ipt.le.ndivpt1d) then
                dncdptoftime(itimesave,ipt)=
     &           dncdptoftime(itimesave,ipt)
     &            +incpt1dy0/5.0
            endif
            endif
         endif
          ipart=partinfo_next(ipart)
      enddo
      timenextsave=timesave(itimesave+1)
      return
      end

      subroutine packonetimedncdptdy(time,timenextsave)
      use PBGgenvar
      use PBGfordisplay
      implicit none
      integer itimesave,ipart,irap,ipt,k
      double precision incpty2d
      double precision time,timenextsave,rap,pt,p(0:3)
      integer ifmtx
      call getMonitorFileIndex(ifmtx)

      itimesave=int(time/timestepforsave)
      if(itimesave.gt.ndivtmax) then 
         write(ifmtx,*) 'itimesave > max foreseen in packonetime'
         close(ifmtx)
         stop
      endif
      incpty2d=incpt1dy0*2.0d0*incrap1d
      ipart=firstpart
      do while (ipart.NE.-1)
         if(partinfo_ires(ipart).eq.1) then
            do k=0,3
               p(k)=partinfo_p(ipart,k)
            enddo
            rap=log((p(0)+p(3))/(p(0)-p(3)))/2.0d0
            irap=1+int((rap-boundrapinf)/steprap1d)
            pt=sqrt(p(1)**2+p(2)**2)
            ipt=int(pt/steppt1d)+1
            if(ipt.le.ndivpt1d.and.irap.le.ndivrap1d) then
            dncdptdyoftime(itimesave,ipt,irap)=
     &           dncdptdyoftime(itimesave,ipt,irap)+incpty2d
            endif
         endif
         if(partinfo_ires(ipart).eq.4) then
            do k=0,3
               p(k)=partinfo_pghost(ipart,k)
            enddo
            rap=log((p(0)+p(3))/(p(0)-p(3)))/2.0d0
            irap=1+int((rap-boundrapinf)/steprap1d)
            pt=sqrt(p(1)**2+p(2)**2)
            ipt=int(pt/steppt1d)+1
            if(ipt.le.ndivpt1d.and.irap.le.ndivrap1d) then
            dncdptdyoftime(itimesave,ipt,irap)=
     &           dncdptdyoftime(itimesave,ipt,irap)+incpt1dy0
            endif
         endif
         ipart=partinfo_next(ipart)
      enddo
      timenextsave=timesave(itimesave+1)
      return
      end


      subroutine packonetime(time,timenextsave)
      use PBGgenvar
      use PBGfordisplay
      implicit none
      integer itimesave,ipart,irap
      double precision time,timenextsave,rap
      integer ifmtx
      call getMonitorFileIndex(ifmtx)
      itimesave=int(time/timestepforsave)
      if(itimesave.gt.ndivtmax) then 
         write(ifmtx,*) 'itimesave > max foreseen in packonetime'
         close(ifmtx)
         stop
      endif
      ipart=firstpart
      do while (ipart.NE.-1)
         if(partinfo_ires(ipart).eq.3) then
            rap=log((partinfo_p(ipart,0)+partinfo_p(ipart,3))/
     &           (partinfo_p(ipart,0)-partinfo_p(ipart,3)))/2
            if(abs(rap).lt.5.5) then
               irap=nint(rap)
               avdnbpsidyoftime(itimesave,irap)=
     &              avdnbpsidyoftime(itimesave,irap)+1.d0
            endif
         endif
         ipart=partinfo_next(ipart)
      enddo
      timenextsave=timesave(itimesave+1)
      return
      end

      subroutine timeevolveNB(nntime,mm)
      use PBGmainpro1
      use PBGgenvar
      use PBGfordisplay
      implicit none
      integer nntime,ipart,nn
      integer mm,mi
      integer vNBc, vNBcbar,vNBD,vNBDbar
      integer timevNBc(1:200),timevNBcbar(1:200),timevNBD(1:200),
     &      timevNBDbar(1:200)
      save timevNBc,timevNBcbar,timevNBD,timevNBDbar
      integer call_og
      save call_og
      if(call_og==0) then
        do mi=1,200
         timevNBc(mi)=0 
         timevNBcbar(mi)=0 
         timevNBD(mi)=0 
         timevNBDbar(mi)=0 
        enddo
        call_og=1
      end if
      vNBc=0
      vNBcbar=0
      vNBD=0
      vNBDbar=0
      ipart=firstpart
      if(mm==0) then
      do while (ipart.NE.-1)
         if(partinfo_ires(ipart).eq.1) then
           vNBc=vNBc+1
         endif
         if(partinfo_ires(ipart).eq.2) then
           vNBcbar=vNBcbar+1
         endif
         if(partinfo_ires(ipart).eq.4) then
           vNBD=vNBD+1
         endif
         if(partinfo_ires(ipart).eq.5) then
           vNBDbar=vNBDbar+1
         endif
           ipart=partinfo_next(ipart)
      enddo
      timevNBc(nntime)=timevNBc(nntime)+vNBc
      timevNBcbar(nntime)=timevNBcbar(nntime)+vNBcbar
      timevNBD(nntime)=timevNBD(nntime)+vNBD
      timevNBDbar(nntime)=timevNBDbar(nntime)+vNBDbar
      else
        ipart=firstpart
        open(98,file='timeevolution_epos1.dat')
        write(98,*) 'particle conservation test timeevolveNB'
        write(98,*) nbcol
        do nn=1,200
        write(98,*) nn,timevNBc(nn)/dble(nbcol)
     &                   ,timevNBcbar(nn)/dble(nbcol)
     &                   ,timevNBD(nn)/dble(nbcol)
     &                   ,timevNBDbar(nn)/dble(nbcol)
        enddo 
        close(98)
      endif
      return
      end 


      subroutine inidistr_hq(mm)
      use PBGmainpro1
      use PBGgenvar
      use PBGhydro
      use PBGfordisplay
      implicit none
      integer mm,ipart

      integer initdistr_c(nbxmax,nbymax),initdistr_cbar(nbxmax,nbymax)
      save initdistr_c,initdistr_cbar
      double precision xx,yy
      integer xbin,ybin
      integer nn, ll
      integer ifmtx
      call getMonitorFileIndex(ifmtx)
      if((nxhyx.gt.nbxmax).or.(nyhyx.gt.nbymax)) then
         write(ifmtx,*) 'out of bound in inidistr_hq'
         write(ifmtx,*) nxhyx,nbxmax,nyhyx,nbxmax
         close(ifmtx)
         stop
      endif
      if(mm==0) then
      ipart=firstpart
      do while (ipart.NE.-1)
         if(partinfo_ires(ipart)==1.or.partinfo_ires(ipart)==2)then
         xx=partinfo_r(ipart,1)
         yy=partinfo_r(ipart,2)
         if((xx.ge.xmin).and.(xx.le.xmax).and.(yy.ge.ymin)
     &       .and.(yy.le.ymax))then
            xbin=int((xx-xmin)/xstep)+1
            ybin=int((yy-ymin)/ystep)+1
            initdistr_c(xbin,ybin)=initdistr_c(xbin,ybin)+1
            initdistr_cbar(xbin,ybin)=initdistr_cbar(xbin,ybin)+1
         end if
         end if
         ipart=partinfo_next(ipart)
      enddo
      else
        open(98,file='inidistr_hq100.dat')
        write(98,*) 'HQ initialisation test inidistr_hq'
        write(98,*) nbcol,xmin,xstep, ymin, ystep
        do nn=1,nxhyx
           do ll=1,nyhyx
              write(98,*) nn,ll,initdistr_c(nn,ll)/dble(nbcol)
     &                   ,initdistr_cbar(nn,ll)/dble(nbcol)
           enddo 
           write(98,*)
        enddo
        close(98)
      endif
      return
      end subroutine inidistr_hq


      subroutine pack_it_1part(ires,p,x,ov)
      use PBGgenvar
      use PBGfordisplay
      use PBGcharge
      implicit none
      logical sav
      integer nhstb,ipt,irap,iphi,ires,iphi_of_phi
      double precision x(0:3),p(0:3),rap,pt,phi,phi_of_p,srap,ov
     &  ,p4JPsi(0:3)
      integer ifmtx
      call getMonitorFileIndex(ifmtx)
      if(ires.gt.nres) then
         write(ifmtx,*) 'lepton seen in pack_it_1part:',ires
         close(ifmtx)
         stop
      endif
      if(ires.le.nresc) then
         nhstb=ires*10
         sav=.true.
      else
         sav=.false.
      endif
      srap=log((x(0)+x(3))/(x(0)-x(3)))/2
      rap=log((p(0)+p(3))/(p(0)-p(3)))/2
      irap=1+int((rap-boundrapinf)/steprap1d)
      pt=sqrt(p(1)**2+p(2)**2)
      ipt=int(pt/steppt1d)+1
      phi=phi_of_p(p(1),p(2)) 
      iphi=iphi_of_phi(phi)
      if((irap.ge.1).and.(irap.le.ndivrap1d)) then
         packrap(ires,irap)=packrap(ires,irap)+ov*incrap1d
         packmeanptsq(ires,irap)=packmeanptsq(ires,irap)+ov*pt**2
      endif   
      if(abs(rap).lt.rapsupobs) then
         if((pt.ne.0).and.sav) then
            if(pt.lt.1.D0) then
            else
               if(pt.lt.4.d0) then
               endif
            endif
         endif               
         if(pt.lt.boundptsup) then
            packpty0(ires,ipt)=packpty0(ires,ipt)+ov*incpt1dy0
            if(ires.eq.1.or.ires.eq.2) then
              packpty0charm(ipt)=packpty0charm(ipt)
     &                             +incpt1dy0
            endif
            if(ires.eq.4.or.ires.eq.5) then
              packpty0dmeson(ipt)=packpty0dmeson(ipt)
     &                             +incpt1dy0
            endif
            if(ncharged.gt.24) then
              packpty0_0_20(ires,ipt)=packpty0_0_20(ires,ipt)
     &                               +ov*incpt1dy0
              packpty0bis_0_20(ires,ipt)=packpty0bis_0_20(ires,ipt)
     &                               +ov*incpt1dy0/pt
            endif
            if(ncharged.le.24.and.ncharged.ge.16) then
              packpty0_20_40(ires,ipt)=packpty0_20_40(ires,ipt)
     &                               +ov*incpt1dy0
              packpty0bis_20_40(ires,ipt)=packpty0bis_20_40(ires,ipt)
     &                               +ov*incpt1dy0/pt    
            endif
            if(ncharged.lt.16.and.ncharged.ge.10) then
              packpty0_40_60(ires,ipt)=packpty0_40_60(ires,ipt)
     &                               +ov*incpt1dy0
              packpty0bis_40_60(ires,ipt)=packpty0bis_40_60(ires,ipt)
     &                               +ov*incpt1dy0/pt
            endif
            if(ncharged.lt.10) then
              packpty0_60_100(ires,ipt)=packpty0_60_100(ires,ipt)
     &                               +ov*incpt1dy0
              packpty0bis_60_100(ires,ipt)=packpty0bis_60_100(ires,ipt)
     &                               +ov*incpt1dy0/pt 
            endif          
            if(abs(p(2)).lt.abs(p(1))) then
              packpty0alongb(ires,ipt)=packpty0alongb(ires,ipt)
     &                              +ov*incpt1dy0
           endif
            packpty0bis(ires,ipt)=packpty0bis(ires,ipt)
     &                           +ov*incpt1dy0/pt
         endif
         packphi(ires,iphi)=packphi(ires,iphi)+ov*incphi1d
         if(pt.lt.boundptsup) then
            packptphiy0(ires,ipt,iphi)=packptphiy0(ires,ipt,iphi)
     &                              +ov*incphipt
            packptphiy0bis(ires,ipt,iphi)=
     &      packptphiy0bis(ires,ipt,iphi)+ov*incphipt/pt
         endif
      endif
      if((abs(rap).lt.2.).and.(abs(rap).gt.1.5)) then
         if(pt.lt.boundptsup) then
            packpty2(ires,ipt)=packpty2(ires,ipt)+ov*incpt1dy2
            packpty2bis(ires,ipt)=packpty2bis(ires,ipt)
     &         +ov*incpt1dy2/pt
         endif
      endif
      if((ires.eq.17).or.(ires.eq.18)) then
         call BdecayJPsi(p,p4JPsi)
         call add_one_JPsi_from_Bdecay(p4JPsi,.false.)
      endif
      return
      end


      subroutine pack_it_1lepton(itypl,p,ov)
      use PBGgenvar
      use PBGfordisplay
      implicit none
      integer itypl,ires,irest,nhstb,ipt,irap,iphi,iphi_of_phi
      double precision p(0:3),ov,rap,pt,phi,phi_of_p,incr
      integer ifmtx
      call getMonitorFileIndex(ifmtx)
      rap=log((p(0)+p(3))/(p(0)-p(3)))/2
      irap=1+int((rap-boundrapinf)/steprap1d)
      pt=sqrt(p(1)**2+p(2)**2)
      ipt=int(pt/steppt1d)+1
      phi=phi_of_p(p(1),p(2)) 
      iphi=iphi_of_phi(phi)
      nhstb=(nres+1)*10
      ires=nres+itypl
      irest=nres+4
      incr=ov*incbas
      if((irap.ge.1).and.(irap.le.ndivrap1d)) then
         packrap(ires,irap)=packrap(ires,irap)+incr
         packrap(irest,irap)=packrap(irest,irap)+incr
         incr=ov*incrap1d*pt**2
         packmeanptsq(ires,irap)=packmeanptsq(ires,irap)+incr
         packmeanptsq(irest,irap)=packmeanptsq(irest,irap)+incr
      endif   
      if(pt.lt.boundptsup) then
         incr=ov*incpt1dy0
         packpty0(ires,ipt)=packpty0(ires,ipt)+incr
         packpty0(irest,ipt)=packpty0(irest,ipt)+incr
         if(abs(p(2)).lt.abs(p(1))) then
            packpty0alongb(ires,ipt)=packpty0alongb(ires,ipt)+incr
            packpty0alongb(irest,ipt)=packpty0alongb(irest,ipt)+incr
         endif
         incr=ov*incpt1dy0/pt
         packpty0bis(ires,ipt)=packpty0bis(ires,ipt)+incr
         packpty0bis(irest,ipt)=packpty0bis(irest,ipt)+incr
      endif
      incr=ov*incphi1d
      packphi(ires,iphi)=packphi(ires,iphi)+incr
      packphi(irest,iphi)=packphi(irest,iphi)+incr
      if(pt.lt.boundptsup) then
         incr=ov*incphipt
         packptphiy0(ires,ipt,iphi)=packptphiy0(ires,ipt,iphi)+incr
         packptphiy0(irest,ipt,iphi)=packptphiy0(irest,ipt,iphi)+incr
         incr=ov*incphipt/pt
         packptphiy0bis(ires,ipt,iphi)=packptphiy0bis(ires,ipt,iphi)
     &                  +incr
         packptphiy0bis(irest,ipt,iphi)=
     &                      packptphiy0bis(irest,ipt,iphi)+incr
      endif
      if((abs(rap).lt.2.).and.(abs(rap).gt.1.5)) then
         if(pt.lt.boundptsup) then
            incr=ov*incpt1dy2
            packpty2(ires,ipt)=packpty2(ires,ipt)+incr
            packpty2(irest,ipt)=packpty2(irest,ipt)+incr
            incr=ov*incpt1dy2/pt
            packpty2bis(ires,ipt)=packpty2bis(ires,ipt)+incr
            packpty2bis(irest,ipt)=packpty2bis(irest,ipt)+incr
         endif
      endif
      return
      end


      subroutine pack_it_2part(ires1,ires2,p1,p2,ov,samemother)
      use PBGgenvar
      use PBGfordisplay
      implicit none
      integer ires1,ires2,iclass,iphi,iphi_of_phi
      logical samemother
      integer ipt1,ipt2
      double precision p1(0:3),p2(0:3),rap1,rap2,pt1,pt2,phi_of_p,phi1,
     &  phi2,deltaphi,ov,pi
      integer ifmtx
      parameter (pi=3.14159265358979323844d0)

      call getMonitorFileIndex(ifmtx)
      if((ires1.le.nres).and.(ires2.le.nres)) then
         if (ires2.ne.resinfo_antipart(ires1)) return
      else
         if(ires1.ne.ires2) then 
            write(ifmtx,*) 'unforseen ires1,ires2 in pack_it_2part'
            write(ifmtx,*) ires1,ires2
            close(ifmtx)
            stop
         endif
      endif
      rap1=log((p1(0)+p1(3))/(p1(0)-p1(3)))/2
      rap2=log((p2(0)+p2(3))/(p2(0)-p2(3)))/2
      pt1=sqrt(p1(1)**2+p1(2)**2)
      pt2=sqrt(p2(1)**2+p2(2)**2)
      phi1=phi_of_p(p1(1),p1(2))
      phi2=phi_of_p(p2(1),p2(2))
      deltaphi=phi2-phi1
	if(deltaphi.lt.0) deltaphi=deltaphi+2*pi
      iphi=iphi_of_phi(deltaphi)
      if((abs(rap1).lt.rapsupobs).and.(abs(rap2).lt.rapsupobs)) then
         ipt1=int(pt1/steppt1d)+1
         ipt2=int(pt2/steppt1d)+1
         if((ipt1.le.ndivpt1d).and.(ipt2.le.ndivpt1d)) then
            packcorrelpty0(ires1,ipt1,ipt2,1)=packcorrelpty0(ires1,ipt1,
     &      ipt2,1)+incpt1dy0**2
            if((abs(deltaphi).ge.0.75*pi).and.(abs(deltaphi).le.
     &           1.25*pi)) packcorrelpty0(ires1,ipt1,ipt2,2)=
     &           packcorrelpty0(ires1,ipt1,ipt2,2)+incpt1dy0**2
            if(samemother) packcorrelpty0(ires1,ipt1,ipt2,3)=
     &           packcorrelpty0(ires1,ipt1,ipt2,3)+incpt1dy0**2
         endif
         iclass=0
         if((pt1.lt.1).and.(pt2.lt.1)) then
            iclass=1
         endif
         if((pt1.ge.1).and.(pt1.lt.4).and.(pt2.ge.1).and.(pt2.lt.4)) 
     &   then
            iclass=2
         endif
         if((pt1.ge.4).and.(pt1.lt.10).and.(pt2.ge.4).and.(pt2.lt.10)) 
     &   then
            iclass=3
         endif
         if((pt1.ge.10).and.(pt1.lt.20).and.(pt2.ge.10).and.(pt2.lt.20))
     &   then
            iclass=4
         endif
         if((pt1.ge.20).and.(pt2.ge.20))then
            iclass=5
         endif
         if(iclass.gt.0) then 
            packcorrelphi(ires1,iclass,iphi)=packcorrelphi(ires1,iclass,
     &   iphi)+ov*incphi1d
         endif
      endif
      return
      end

      subroutine gen_all_leptons
      use PBGgenvar
      use PBGforleptons
      implicit none
      integer ipart,k,ires,ifail,imed
      double precision pparent(0:3),plept(0:3),x(0:3),temp,betaplasma(3)
     & ,densener,ran2 
      integer ifmtx
      call getMonitorFileIndex(ifmtx)
      nblepton=0
      ipart=firstpart
      do while (ipart.NE.-1)
         ires=partinfo_ires(ipart)
         if((ires.eq.1).or.(ires.eq.2).or.(ires.eq.14)
     &            .or.(ires.eq.15)) then
            write(ifmtx,*) 'not a meson in gen_all_leptons',ires,ipart
            write(ifmtx,*) (partinfo_p(ipart,k),k=0,3)
            do k=0,3
               x(k)=partinfo_p(ipart,k)
            enddo
            call givelocalhydro(x,temp,betaplasma,densener,ifail,imed)
            write(ifmtx,*) temp,ifail,imed
            close(ifmtx)
            stop
         endif
         if((ires.eq.4).or.(ires.eq.5)) then
            if(ran2().lt.0.108) then
               if(nblepton.eq.nbmaxlept) then
                  write(ifmtx,*) 'too many leptons in gen_all_leptons'
                  close(ifmtx)
                  stop
               else
                  nblepton=nblepton+1
               endif
               typlepton(nblepton)=ires-3
               origlepton(nblepton)=1
               do k=0,3
                  pparent(k)=partinfo_p(ipart,k)
               enddo
               call weakdecayD(pparent,plept)
               do k=0,3
                  tbplepton(nblepton,k)=plept(k)
               enddo
            endif
         else if((ires.eq.17).or.(ires.eq.18)) then
            if(ran2().lt.0.106) then
               if(nblepton.eq.nbmaxlept) then
                  write(ifmtx,*) 'too many leptons in gen_all_leptons'
                  close(ifmtx)
                  stop
               else
                  nblepton=nblepton+1
               endif
               typlepton(nblepton)=ires-6
               origlepton(nblepton)=2
               do k=0,3
                  pparent(k)=partinfo_p(ipart,k)
               enddo
               call weakdecayB(pparent,plept)
               do k=0,3
                  tbplepton(nblepton,k)=plept(k)
               enddo
            endif
            if(ran2().lt.0.108) then
               if(nblepton.eq.nbmaxlept) then
                  write(ifmtx,*) 'too many leptons in gen_all_leptons'
                  close(ifmtx)
                  stop
               else
                  nblepton=nblepton+1
               endif
               typlepton(nblepton)=ires-6
               origlepton(nblepton)=3
               do k=0,3
                  pparent(k)=partinfo_p(ipart,k)
               enddo
               call weakdecayBsecond(pparent,plept)
               do k=0,3
                  tbplepton(nblepton,k)=plept(k)
               enddo
            endif
         endif   
         ipart=partinfo_next(ipart)
      enddo
      end


      subroutine pack_it
      use PBGgenvar
      use PBGfordisplay
      use PBGforleptons
      implicit none
      integer ipart,ipart2,k,nbjpsi,ires,ires2
      double precision x(0:3),p(0:3),x2(0:3),p2(0:3),ov1,ov2

      ov1=1.d0
      ov2=0.01d0
      ipart=firstpart
      do while (ipart.NE.-1)
         do k=0,3
            x(k)=partinfo_r(ipart,k)
            p(k)=partinfo_p(ipart,k)
         enddo
         ipart2=firstpart                
         do while (ipart2.NE.-1)
             do k=0,3
               x2(k)=partinfo_r(ipart2,k)
               p2(k)=partinfo_p(ipart2,k)
            enddo
            if(ipart.ne.ipart2) then
               call pack_it_2part(partinfo_ires(ipart),
     &                            partinfo_ires(ipart2),
     &                            p,p2,ov1,.false.)
            endif
            ipart2=partinfo_next(ipart2)
         enddo
         call pack_it_1part(partinfo_ires(ipart),p,x,ov1)
         ipart=partinfo_next(ipart)
      enddo
      ipart=firstpart
      do while (ipart.NE.-1)
         ires=partinfo_ires(ipart)
         if((ires.eq.4).or.(ires.eq.5).or.(ires.eq.17)
     &        .or.(ires.eq.18)) then
            do k=0,3
               p(k)=partinfo_pghost(ipart,k)
               x(k)=partinfo_rghost(ipart,k)
            enddo
            ipart2=firstpart                
            do while (ipart2.NE.-1)
               ires2=partinfo_ires(ipart2)
               if(abs(ires-ires2).eq.1) then
                  do k=0,3
                     p2(k)=partinfo_pghost(ipart2,k)
                     x2(k)=partinfo_rghost(ipart2,k)
                  enddo
                  call pack_it_2part(ires-3,ires2-3,p,p2,ov1,.false.)
               endif
               ipart2=partinfo_next(ipart2)
            enddo
            call pack_it_1part(ires-3,p,x,ov1)
         endif
         ipart=partinfo_next(ipart)
      enddo


      nbjpsi=resinfo_number(3)      
      if(nbjpsi.lt.ndivpsimax) packnbpsi(nbjpsi)=packnbpsi(nbjpsi)+1 
      return
      end


      subroutine close_pack(nbsuccevent,glauber)
      use PBGmainpro1
      use PBGfordisplay
      use PBGievent
      use PBGglaubercentrality
      implicit none
#include "aaa.h" 
      integer ibase,i,j,k,nbsuccevent,ires,il,itime,itypq
      double precision pt,rap,phi,ratio1,ratio2,ratio3,ratio4,ratio5,
     & meanptsq(0:nres+4),meanvacptsq(0:nres+4),v2(0:nres+4,ndivpt1d),
     & Snum(0:nres+4,ndivpt1d),Sden(0:nres+4,ndivpt1d),
     & v2nonpromptpsi(ndivpt1d),Snumnonpromptpsi(ndivpt1d),
     & Sdennonpromptpsi(ndivpt1d),
     & ratio6, ratio7,ratiol(4)
      integer glauber
      character*1 namequark(2)
      data namequark/'c','b'/
     
 151  format(F10.3,' ',16(E12.5,' '))

      open(unit=ifnus2,file=fnus2(1:nfnus2),status='unknown')      
      ibase=0.

       write(ifnus2,*)nbsuccevent,glauber*nbcol/float(nbsuccevent)

      write(ifnus2,*) nbsuccevent,iev_0_20*nbcol,iev_20_40*nbcol,
     &                iev_40_60*nbcol,iev_60_100*nbcol

      if(iev_0_20.ne.0) then
        glauber_0_20=glauber_0_20/float(iev_0_20)
       else
        glauber_0_20=0
      endif
      if(iev_20_40.ne.0) then
        glauber_20_40=glauber_20_40/float(iev_20_40)
       else
        glauber_20_40=0
      endif

      if(iev_40_60.ne.0) then
        glauber_40_60=glauber_40_60/float(iev_40_60)
       else
        glauber_40_60=0
      endif 
      if(iev_60_100.ne.0) then
        glauber_60_100=glauber_60_100/float(iev_60_100)
       else
        glauber_60_100=0
      endif  


      write(ifnus2,*) glauber*nbcol/float(nbsuccevent),glauber_0_20,
     &                glauber_20_40,glauber_40_60,glauber_60_100

      write(ifnus2,*) '  '


      write(ifnus2,*) '******* pt *******  '
      write(ifnus2,*) '  '
      write(ifnus2,*) '1 mean pt**2 of charm (c,J/psi,D(itypq=1))and 
     &bottom (b,Y,B,itypq=2)' 
      write(ifnus2,*) 'meanptcb.res'
      write(ifnus2,*) ndivrap1d
      do itypq=1,2
         ibase=(itypq-1)*nresc
         do i=1,ndivrap1d
            rap=(i-0.5)*steprap1d+boundrapinf
            do j=ibase+1,ibase+1+nresc
               if(packrap(j,i).gt.0.)then   
                  meanptsq(j)=packmeanptsq(j,i)/packrap(j,i)
               else
                  meanptsq(j)=0.d0
               endif
               if(packvacrap(j,i).gt.0.)then   
                  meanvacptsq(j)=packvacmeanptsq(j,i)/packvacrap(j,i)
               else
                  meanvacptsq(j)=0.d0
               endif
            enddo
            write(ifnus2,151) rap,sqrt(meanptsq(ibase+4)),
     &       sqrt(meanptsq(ibase+3)),
     &       sqrt(meanvacptsq(ibase+4)),sqrt(meanvacptsq(ibase+3)),
     &       sqrt(meanptsq(ibase+1)),sqrt(meanvacptsq(ibase+1))
         enddo
	 write(ifnus2,*)
      enddo



      Write(ifnus2,*) '  '
      write(ifnus2,*) '2 mean pt**2 for all nonphotonic electrons'
      write(ifnus2,*) 'rmspt_elec.res'
      write(ifnus2,*) ndivrap1d



      do i=1,ndivrap1d
         rap=(i-0.5)*steprap1d+boundrapinf
         do j=nres+1,nres+4
            if(packrap(j,i).gt.0.)then   
               meanptsq(j)=packmeanptsq(j,i)/packrap(j,i)
            else
               meanptsq(j)=0.d0
            endif
            if(packvacrap(j,i).gt.0.)then   
               meanvacptsq(j)=packvacmeanptsq(j,i)/packvacrap(j,i)
            else
               meanvacptsq(j)=0.d0
            endif
         enddo


         write(ifnus2,151) rap,(sqrt(meanptsq(j)),j=nres+1,nres+4),
     &        (sqrt(meanvacptsq(j)),j=nres+1,nres+4)
      enddo

      write(ifnus2,*) ' '
      write(ifnus2,*) '3 ptdpt distribution c,D,J/psi(itypq=0) and b,Y, 
     &B(itypq=1) for |y|<1'
      write(ifnus2,*) 'distr_ptdpt_y0.res'
      write(ifnus2,*) ndivpt1d
      do itypq=1,2
         ibase=(itypq-1)*nresc         
         do i=1,ndivpt1d
            pt=(i-0.5)*steppt1d
            if(packpty0(ibase+4,i).gt.0) then 
               ratio1=packpty0(ibase+3,i)/packpty0(ibase+4,i)
            else
               ratio1=0.
            endif
            if(packvacpty0(ibase+4,i).gt.0) then 
               ratio2=packpty0(ibase+4,i)/packvacpty0(ibase+4,i)
               ratio6=2.d0*packpty0alongb(ibase+4,i)/
     &              packvacpty0(ibase+4,i)
               ratio7=2*ratio2-ratio6
            else
               ratio2=0.
               ratio6=0.d0
               ratio7=0.d0
            endif
            if(packvacpty0(ibase+3,i).gt.0) then 
               ratio3=packpty0(ibase+3,i)/packvacpty0(ibase+3,i)
            else
               ratio3=0.
            endif
            if(packvacpty0(ibase+1,i).gt.0) then 
               ratio4=packpty0(ibase+1,i)/packvacpty0(ibase+1,i)
            else
               ratio4=0.
            endif
            if(packvacpty0(ibase+3,i).gt.0) then 
               ratio5=packpty0(ibase+3,i)/packvacpty0(ibase+3,i)
            else
               ratio5=0.
            endif
            write(ifnus2,151) pt,packpty0(ibase+4,i)/nbsuccevent,
     &        packpty0(ibase+3,i)/nbsuccevent,packpty0bis(ibase+4,i)/
     &        nbsuccevent,packpty0bis(ibase+3,i)/nbsuccevent,ratio1,
     &        packvacpty0(ibase+4,i)/nbsuccevent,
     &        packvacpty0bis(ibase+4,i)/
     &        nbsuccevent,ratio2,packpty0bis(ibase+1,i)/nbsuccevent,
     &        packvacpty0bis(ibase+1,i)/nbsuccevent,ratio4,
     &        packvacpty0bis(ibase+3,i)/nbsuccevent,ratio5,ratio6,
     &        ratio7
            enddo
        write(ifnus2,*)'  '
      enddo
      write(ifnus2,*)'  '
      write(ifnus2,*)'4 p_t distributions |y|<1'
      write(ifnus2,*) 'distr_ptdpt_y0_bis.res'
      write(ifnus2,*) ndivpt1d,steppt1d
      do i=1,ndivpt1d  
         pt=(i-0.5)*steppt1d
         write(ifnus2,151) pt,
     &        packpty0dmeson(i)/nbsuccevent/2,
     &                        packvacpty0dmeson(i)/nbsuccevent/2,
     &                        packpty0charm(i)/nbsuccevent/2,
     &                        packvacpty0charm(i)/nbsuccevent/2,
     &                        packpty0(ibase+4,i)/nbsuccevent,
     &                        packvacpty0(ibase+4,i)/nbsuccevent,
     &                        packpty0(ibase+1,i)/nbsuccevent,
     &                        packvacpty0(ibase+1,i)/nbsuccevent
      enddo 
      if(iev_0_20*nbcol.eq.0)iev_0_20=1
      write(ifnus2,*) ' '
      write(ifnus2,*) '5 distr_ptdpt_y0_0_20,  |y|<1 c,D(itypq=0) b,
     &   Y(itypq=1)'
      write(ifnus2,*)'distr_ptc_y0_0_20.res' 
      write(ifnus2,*) iev_0_20*nbcol,glauber_0_20/iev_0_20*nbcol
      write(ifnus2,*) ndivpt1d,itypq
      do i=1,ndivpt1d
         pt=(i-0.5)*steppt1d
         if(iev_0_20.eq.0) then
            packpty0_0_20(ibase+4,i)=0
            packpty0bis_0_20(ibase+1,i)=0
            packvacpty0bis_0_20(ibase+1,i)=0 
         endif 
         write(ifnus2,*)pt,
     &           packpty0_0_20(ibase+4,i)/float(iev_0_20*nbcol),
     &           packpty0bis_0_20(ibase+1,i)/float(iev_0_20*nbcol),
     &           packvacpty0bis_0_20(ibase+1,i)/float(iev_0_20*nbcol)
      enddo
      if(iev_20_40.eq.0)iev_20_40=1
      write(ifnus2,*) ' '
      write(ifnus2,*) 
     &   '6 pt distribution 20_0_40 |y|<1 c,D(itypq=0) b,Y(itypq=1)'
      write(ifnus2,*)'distr_ptc_y0_20_40.res' 
      write(ifnus2,*) iev_20_40*nbcol,glauber_20_40/float(iev_20_40)
      write(ifnus2,*) ndivpt1d
      do i=1,ndivpt1d
         pt=(i-0.5)*steppt1d
         write(ifnus2,*)pt,
     &           packpty0_20_40(ibase+4,i)/float(iev_20_40*nbcol),
     &           packpty0bis_20_40(ibase+1,i)/float(iev_20_40*nbcol),
     &           packvacpty0bis_20_40(ibase+1,i)/float(iev_20_40*nbcol)
      enddo
      if(iev_40_60.eq.0)iev_40_60=1
      if(iev_40_60.ne.0) then
         write(ifnus2,*) ' '
         write(ifnus2,*) 
     &     '7 pt distribution 40_0_60 |y|<1 c,D(itypq=0),b,Y(itypq=1)'
         write(ifnus2,*)'distr_ptc_y0_40_60.res' 
         write(ifnus2,*) iev_40_60*nbcol,glauber_40_60/float
     &     (iev_40_60)
         write(ifnus2,*) ndivpt1d,itypq
         do i=1,ndivpt1d
            pt=(i-0.5)*steppt1d
            write(ifnus2,*)pt,
     &             packpty0_40_60(ibase+4,i)/float(iev_40_60*nbcol),
     &             packpty0bis_40_60(ibase+1,i)/float(iev_40_60*nbcol),
     &             packvacpty0bis_40_60(ibase+1,i)/
     &             float(iev_40_60*nbcol)
         enddo
      endif
      if(iev_60_100.eq.0)iev_60_100=1
      if(iev_60_100.ne.0) then
         write(ifnus2,*) ' '
         write(ifnus2,*) 
     &     '8 pt distribution 60_0_100 |y|<1 c,D(itypq=0) b,Y(itypq=1)'
         write(ifnus2,*)'distr_ptc_y0_60_100.res' 
         write(ifnus2,*) iev_60_100*nbcol,glauber_60_100/
     &     float(iev_60_100)
         write(ifnus2,*) ndivpt1d
         do i=1,ndivpt1d
            pt=(i-0.5)*steppt1d
            write(ifnus2,*)pt,
     &            packpty0_60_100(ibase+4,i)/float(iev_60_100*nbcol),
     &            packpty0bis_60_100(ibase+1,i)/float(iev_60_100*nbcol),
     &            packvacpty0bis_60_100(ibase+1,i)/float(iev_60_100*
     &            nbcol)
         enddo
      endif
      write(ifnus2,*) ' '
      write(ifnus2,*) '9 pt distribution c,D,J/psi(itypq=0)'
     &   ,' and b,Y,B(itypq=1) for <|y|<'
      write(ifnus2,*) 'distr_ptc_y2.res'
      write(ifnus2,*) ndivpt1d,itypq
      do i=1,ndivpt1d
         pt=(i-0.5)*steppt1d
         if(packpty2(ibase+4,i).gt.0) then 
            ratio1=packpty2(ibase+3,i)/packpty2(ibase+4,i)
         else
            ratio1=0.
         endif
         if(packvacpty2(ibase+4,i).gt.0) then 
            ratio2=packpty2(ibase+4,i)/packvacpty2(ibase+4,i)
         else
            ratio2=0.
         endif
         if(packvacpty2(ibase+3,i).gt.0) then 
            ratio3=packpty2(ibase+3,i)/packvacpty2(ibase+3,i)
         else
            ratio3=0.
         endif
         if(packvacpty2(ibase+1,i).gt.0) then 
            ratio4=packpty2(ibase+1,i)/packvacpty2(ibase+1,i)
         else
            ratio4=0.
         endif
         if(packvacpty2(ibase+3,i).gt.0) then 
            ratio5=packpty2(ibase+3,i)/packvacpty2(ibase+3,i)
         else
            ratio5=0.
         endif
         write(ifnus2,151) pt,packpty2(ibase+4,i)/nbsuccevent,
     &        packpty2(ibase+3,i)/nbsuccevent,packpty2bis(ibase+4,i)/
     &        nbsuccevent,packpty2bis(ibase+3,i)/nbsuccevent,ratio1,
     &        packvacpty2(ibase+4,i)/nbsuccevent,
     &        packvacpty2bis(ibase+4,i)/
     &        nbsuccevent,ratio2,packpty2bis(ibase+1,i)/nbsuccevent,
     &        packvacpty2bis(ibase+1,i)/nbsuccevent,ratio4,
     &        packvacpty2bis(ibase+3,i)/nbsuccevent,ratio5
      enddo
      Write(ifnus2,*) '  '
      write(ifnus2,*) '10 pt distribution prompy J/psi for |y|<1'
      Write(ifnus2,*) 'distr_pt_y0_promptJ/psi.res'
      write(ifnus2,*)  ndivpt1d
      do i=1,ndivpt1d
         pt=(i-0.5)*steppt1d
         if(packvacpty0nonpromptJpsi(i).gt.0) then 
            ratio1=packpty0nonpromptJpsi(i)/packvacpty0nonpromptJpsi(i)
         else
            ratio1=0.
         endif
         write(ifnus2,151) pt,packpty0nonpromptJpsi(i)/nbsuccevent,
     &         packvacpty0nonpromptJpsi(i)/nbsuccevent,ratio1
      enddo
      Write(ifnus2,*) '  '
      write(ifnus2,*) '13 the tranverse pt spectra of non-photonic '
     &   ,'leptons (at mid-rapidity)'
      Write(ifnus2,*) 'distr_ptc_elec_y0.res'
      write(ifnus2,*) ndivpt1d
      do i=1,ndivpt1d
         pt=(i-0.5)*steppt1d
         do il=1,4
            if(packvacpty0(nres+il,i).ne.0) then
               ratiol(il)=packpty0bis(nres+il,i)/
     &              packvacpty0bis(nres+il,i)
            else
               ratiol(il)=0.
            endif
         enddo
         write(ifnus2,151) pt,
     &    (packpty0bis(nres+il,i)/nbsuccevent,il=1,4),
     &    (packvacpty0bis(nres+il,i)/nbsuccevent,il=1,4),
     &    (ratiol(il),il=1,4)
      enddo
      Write(ifnus2,*) '  '
      write(ifnus2,*) '14 the tranverse pt spectra of non-photonic '
     &   ,'leptons (at all rapidity)'
      Write(ifnus2,*) 'distr_ptc_elec_y2.res'
      write(ifnus2,*) ndivpt1d
      do i=1,ndivpt1d
         pt=(i-0.5)*steppt1d
         do il=1,4
            if(packvacpty2(nres+il,i).ne.0) then
               ratiol(il)=packpty2bis(nres+il,i)/
     &              packvacpty2bis(nres+il,i)
            else
               ratiol(il)=0.
            endif
         enddo
         write(ifnus2,151) pt,
     &           (packpty2bis(nres+il,i)/nbsuccevent,il=1,4),
     &    (packvacpty2bis(nres+il,i)/nbsuccevent,il=1,4),
     &           (ratiol(il),il=1,4)
      enddo
      write(ifnus2,*) '******* rap *******  '
      Write(ifnus2,*) '  '
      write(ifnus2,*) '15 rap spectra c ccbar b bbar'
      Write(ifnus2,*) 'y_dist_hq.res'
      write(ifnus2,*) ndivrap1d
      do itypq=1,2
         ibase=(itypq-1)*nresc
         do i=1,ndivrap1d
            rap=(i-0.5)*steprap1d+boundrapinf
            if(packrap(ibase+4,i).gt.0) then 
               ratio1=packrap(ibase+3,i)/packrap(ibase+4,i)
               ratio2
     .            =packrap(ibase+3,i)/packrap(ibase+4,i)**2*nbsuccevent
            else
               ratio1=0.
               ratio2=0.
            endif
            if(packvacrap_zero(ibase+4,i).gt.0) then 
               ratio6=packrap(ibase+4,i)/packvacrap_zero(ibase+4,i)
            else
               ratio6=0.
            endif
             if(packvacrap(ibase+4,i).gt.0) then 
               ratio3=packrap(ibase+4,i)/packvacrap(ibase+4,i)
            else
               ratio3=0.
            endif
            if(packvacrap_zero(ibase+1,i).gt.0) then 
               ratio7=packrap(ibase+1,i)/packvacrap_zero(ibase+1,i)
            else
               ratio7=0.
            endif
             if(packvacrap(ibase+1,i).gt.0) then 
               ratio4=packrap(ibase+1,i)/packvacrap(ibase+1,i)
            else
               ratio4=0.
            endif
            if(packvacrap(ibase+3,i).gt.0) then 
               ratio5=packrap(ibase+3,i)/packvacrap(ibase+3,i)
            else
               ratio5=0.
            endif
            write(ifnus2,151) rap,packrap(ibase+4,i)/nbsuccevent,
     &        packrap(ibase+3,i)/nbsuccevent,ratio1,ratio2,
     &        packvacrap(ibase+4,i)/nbsuccevent,
     &         ratio3,packrap(ibase+1,i)/
     &        nbsuccevent,packvacrap(ibase+1,i)/nbsuccevent,ratio4,
     &        packvacrap(ibase+3,i)/nbsuccevent,ratio5,
     &        packvacrap_zero(ibase+4,i)/nbsuccevent,ratio6,
     &        packvacrap_zero(ibase+1,i)/nbsuccevent,ratio7
         enddo
	 Write(ifnus2,*) '  '
      enddo
      Write(ifnus2,*) '  '
      write(ifnus2,*) '16 rap spectra non-photonic e'
      Write(ifnus2,*) 'y_dist_elec.res'
      write(ifnus2,*) ndivrap1d
      do i=1,ndivrap1d
         rap=(i-0.5)*steprap1d+boundrapinf
         do il=1,4
            if(packvacrap(nres+il,i).ne.0) then
               ratiol(il)=packrap(nres+il,i)/packvacrap(nres+il,i)
            else
               ratiol(il)=0.
            endif
         enddo
         write(ifnus2,151) rap,(packrap(nres+il,i)/nbsuccevent,il=1,4),
     &     (packvacrap(nres+il,i)/nbsuccevent,il=1,4),
     &     (ratiol(il),il=1,4)
      enddo

      write(ifnus2,*) '******* phi *******  '
      Write(ifnus2,*) '  '
      write(ifnus2,*) '17 phi-spectra of const quarks'
     &    ,' at mid rapidity ??, ??(for various pt)'
      Write(ifnus2,*) 'dn_dphi_hq_y0.res'
      write(ifnus2,*) ndivphi1d,stepphi1d
      do j=1,ndivphi1d
         phi=(j-0.5)*stepphi1d
         write(ifnus2,151) phi,
     &    packptphiy0bis(0,5,j)/nbsuccevent,
     &    packptphiy0bis(0,10,j)/nbsuccevent,
     &    packptphiy0bis(0,15,j)/nbsuccevent,
     &    packptphiy0bis(0,20,j)/nbsuccevent,
     &    packptphiy0bis(0,30,j)/nbsuccevent,
     &    packptphiy0bis(0,40,j)/nbsuccevent
      enddo
      Write(ifnus2,*) '  '
      write(ifnus2,*) '18 phi-spectra of D mesons at mid rapidity??,'
     &   ,'??(for various pt)'
      Write(ifnus2,*) 'distr_phiD_y0.res'
      write(ifnus2,*) ndivphi1d
      do j=1,ndivphi1d
         phi=(j-0.5)*stepphi1d
         write(ifnus2,151) phi,
     &    packptphiy0bis(4,5,j)/nbsuccevent,
     &    packptphiy0bis(4,10,j)/nbsuccevent,
     &    packptphiy0bis(4,15,j)/nbsuccevent,
     &    packptphiy0bis(4,20,j)/nbsuccevent,
     &    packptphiy0bis(4,30,j)/nbsuccevent,
     &    packptphiy0bis(4,40,j)/nbsuccevent
      enddo
      Write(ifnus2,*) '  '
      write(ifnus2,*) '19 phi-spectra of leptons (all)'
     &   ,'at mid rapidity ??'
     &   ,'??for various pt'
      Write(ifnus2,*) 'distr_philept_y0.res'
      write(ifnus2,*) ndivphi1d
      do j=1,ndivphi1d
         phi=(j-0.5)*stepphi1d
         write(ifnus2,151) phi,
     &    packptphiy0bis(nres+4,5,j)/nbsuccevent,
     &    packptphiy0bis(nres+4,10,j)/nbsuccevent,
     &    packptphiy0bis(nres+4,15,j)/nbsuccevent,
     &    packptphiy0bis(nres+4,20,j)/nbsuccevent,
     &    packptphiy0bis(nres+4,30,j)/nbsuccevent,
     &    packptphiy0bis(nres+4,40,j)/nbsuccevent
      enddo
      Write(ifnus2,*) '  '
      write(ifnus2,*) ' 20 distr_phirel  -cbar correlations'
      Write(ifnus2,*) 'distr_phirel c cbar_y0.res'
      write(ifnus2,*) ndivphi1d
      do j=1,ndivphi1d
         phi=(j-0.5)*stepphi1d
         write(ifnus2,151) phi,packcorrelphi(1,1,j)/nbsuccevent,
     &    packcorrelphi(1,2,j)/nbsuccevent,packcorrelphi(1,3,j)
     &    /nbsuccevent,
     &    packcorrelphi(1,4,j)/nbsuccevent,packcorrelphi(1,5,j)
     &    /nbsuccevent
      enddo
C      open(9,file='distr_phirelbbbar_y0.res')
      Write(ifnus2,*) '  '
      write(ifnus2,*) '21 b-bbar correlations'
      Write(ifnus2,*) 'distr_phirelbbbar_y0.res'
      write(ifnus2,*) ndivphi1d
      do j=1,ndivphi1d
         phi=(j-0.5)*stepphi1d
         write(ifnus2,151) phi,packcorrelphi(6,1,j)/nbsuccevent,
     &    packcorrelphi(6,2,j)/nbsuccevent,packcorrelphi(6,3,j)
     &    /nbsuccevent,
     &    packcorrelphi(6,4,j)/nbsuccevent,packcorrelphi(6,5,j)
     &    /nbsuccevent
      enddo
C      open(9,file='distr_phirelDDbar_y0.res')
      Write(ifnus2,*) '  '
      write(ifnus2,*) '22 D-Dbar correlations'
      Write(ifnus2,*) 'distr_phirelDDbar_y0.res'
      write(ifnus2,*) ndivphi1d
      do j=1,ndivphi1d
         phi=(j-0.5)*stepphi1d
         write(ifnus2,151) phi,
     &    packcorrelphi(4,1,j)/nbsuccevent,
     &    packcorrelphi(4,2,j)/nbsuccevent,
     &    packcorrelphi(4,3,j)/nbsuccevent,
     &    packcorrelphi(4,4,j)/nbsuccevent,
     &    packcorrelphi(4,5,j)/nbsuccevent
      enddo
      Write(ifnus2,*) '  '
      write(ifnus2,*) '23 B-Bbar correlations'
      Write(ifnus2,*) 'distr_phirelBmesBbarmes_y0.res'
      write(ifnus2,*) ndivphi1d
      do j=1,ndivphi1d
         phi=(j-0.5)*stepphi1d
         write(ifnus2,151) phi,
     &    packcorrelphi(9,1,j)/nbsuccevent,
     &    packcorrelphi(9,2,j)/nbsuccevent,
     &    packcorrelphi(9,3,j)/nbsuccevent,
     &    packcorrelphi(9,4,j)/nbsuccevent,
     &    packcorrelphi(9,5,j)/nbsuccevent
      enddo
      Write(ifnus2,*) '  '
      write(ifnus2,*) '24 ep_em_phi correlations'
      Write(ifnus2,*) 'distr_phirelee_y0.res'
      write(ifnus2,*) ndivphi1d
      do j=1,ndivphi1d
         phi=(j-0.5)*stepphi1d
         write(ifnus2,151) phi,
     &    packcorrelphi(nres+4,1,j)/nbsuccevent,
     &    packcorrelphi(nres+4,2,j)/nbsuccevent,
     &    packcorrelphi(nres+4,3,j)/nbsuccevent,
     &    packcorrelphi(nres+4,4,j)/nbsuccevent,
     &    packcorrelphi(nres+4,5,j)/nbsuccevent
      enddo
      do i=1,ndivpt1d
         do j=1,ndivpt1d
         enddo
      enddo
      do i=1,ndivpt1d
         do j=1,ndivpt1d
         enddo
      enddo


      do ires=0,nres+4
         do j=1,ndivpt1d
            Snum(ires,j)=0.d0
            Sden(ires,j)=0.d0
            do k=1,ndivphi1d
               phi=(k-0.5)*stepphi1d
               Snum(ires,j)=Snum(ires,j)+dcos(2*phi)
     &                     *packptphiy0(ires,j,k)
               Sden(ires,j)=Sden(ires,j)+packptphiy0(ires,j,k)
            enddo
            Sden(ires,j)=Sden(ires,j)*stepphi1d*steppt1d
            Snum(ires,j)=Snum(ires,j)*stepphi1d*steppt1d
            if(Sden(ires,j).gt.0) then
               v2(ires,j)=Snum(ires,j)/Sden(ires,j)
            else
               v2(ires,j)=0.
            endif
         enddo
      enddo
      do j=1,ndivpt1d
         Snumnonpromptpsi(j)=0.d0
         Sdennonpromptpsi(j)=0.d0
         do k=1,ndivphi1d
            phi=(k-0.5)*stepphi1d
            Snumnonpromptpsi(j)=Snumnonpromptpsi(j)+dcos(2*phi)*
     &           packptphiy0nonpromptJpsi(j,k)
            Sdennonpromptpsi(j)=Sdennonpromptpsi(j)+
     &           packptphiy0nonpromptJpsi(j,k)
            enddo
            Sdennonpromptpsi(j)=Sdennonpromptpsi(j)*stepphi1d*steppt1d
            Snumnonpromptpsi(j)=Snumnonpromptpsi(j)*stepphi1d*steppt1d
            if(Sdennonpromptpsi(j).gt.0) then
               v2nonpromptpsi(j)=Snumnonpromptpsi(j)/Sdennonpromptpsi(j)
            else
               v2nonpromptpsi(j)=0.d0
            endif
         enddo

      write(ifnus2,*) '******* v2 *******  '
      write(ifnus2,*) ' '
      write(ifnus2,*) '27 elliptic flow of c at midrapidity  '
      write(ifnus2,*) 'elliptic_flowc_y0.res'
      write(ifnus2,*) ndivpt1d
      do j=1,ndivpt1d
         pt=(j-0.5)*steppt1d
         write(ifnus2,151) pt,v2(1,j),Snum(1,j),Sden(1,j)
      enddo
      write(ifnus2,*) ' '
      write(ifnus2,*) '28 elliptic flow of D at midrapidity  '
      write(ifnus2,*) 'elliptic_flowD_y0.res'
      write(ifnus2,*) ndivpt1d
      do j=1,ndivpt1d
         pt=(j-0.5)*steppt1d
         write(ifnus2,151) pt,v2(4,j),Snum(4,j),Sden(4,j)
      enddo
      do ires=0,nres+4
         do j=1,ndivpt1d
            Snum(ires,j)=0.d0
            Sden(ires,j)=0.d0
            do k=1,ndivphi1d
               phi=(k-0.5)*stepphi1d
               Snum(ires,j)
     &        =Snum(ires,j)+dcos(2*phi)*packptphiy0(ires,j,k)
               Sden(ires,j)=Sden(ires,j)+packptphiy0(ires,j,k)
            enddo
            Sden(ires,j)=Sden(ires,j)*stepphi1d*steppt1d
            Snum(ires,j)=Snum(ires,j)*stepphi1d*steppt1d
            if(Sden(ires,j).gt.0) then
               v2(ires,j)=Snum(ires,j)/Sden(ires,j)
            else
               v2(ires,j)=0.
            endif
         enddo
      enddo
      do j=1,ndivpt1d
         Snumnonpromptpsi(j)=0.d0
         Sdennonpromptpsi(j)=0.d0
         do k=1,ndivphi1d
            phi=(k-0.5)*stepphi1d
            Snumnonpromptpsi(j)=Snumnonpromptpsi(j)+dcos(2*phi)*
     &           packptphiy0nonpromptJpsi(j,k)
            Sdennonpromptpsi(j)=Sdennonpromptpsi(j)+
     &           packptphiy0nonpromptJpsi(j,k)
         enddo
         Sdennonpromptpsi(j)=Sdennonpromptpsi(j)*stepphi1d*steppt1d
         Snumnonpromptpsi(j)=Snumnonpromptpsi(j)*stepphi1d*steppt1d
         if(Sdennonpromptpsi(j).gt.0) then
            v2nonpromptpsi(j)=Snumnonpromptpsi(j)/Sdennonpromptpsi(j)
         else
            v2nonpromptpsi(j)=0.d0
         endif
      enddo
      Write(ifnus2,*) '  '
      write(ifnus2,*) '29 elliptic flow of c at midrapidity    '
      write(ifnus2,*) 'elliptic_flow_const_y0.res'
      write(ifnus2,*) ndivpt1d
      do j=1,ndivpt1d
         pt=(j-0.5)*steppt1d
         write(ifnus2,151) pt,v2(0,j),Snum(0,j),Sden(0,j)
      enddo
      Write(ifnus2,*) '  '
      write(ifnus2,*) '30 '
      write(ifnus2,*) 'elliptic_flowc_y0_bis.res'
      write(ifnus2,*) ndivpt1d
      do j=1,ndivpt1d
         pt=(j-0.5)*steppt1d
         write(ifnus2,151) pt,v2(1,j),Snum(1,j),Sden(1,j)
      enddo
      Write(ifnus2,*) '  '      
      write(ifnus2,*) '31 '
      write(ifnus2,*) 'elliptic_flowD_y0bis.res'
      write(ifnus2,*) ndivpt1d
      do j=1,ndivpt1d
         pt=(j-0.5)*steppt1d
         write(ifnus2,151) pt,v2(4,j),Snum(4,j),Sden(4,j)
      enddo
      Write(ifnus2,*) '  '
      write(ifnus2,*) '32 '
      write(ifnus2,*) 'elliptic_flowb_y0.res'
      write(ifnus2,*) ndivpt1d
      do j=1,ndivpt1d
         pt=(j-0.5)*steppt1d
         write(ifnus2,151) pt,v2(6,j),Snum(6,j),Sden(6,j)
      enddo
      Write(ifnus2,*) '  '
      write(ifnus2,*) '33 '
      write(ifnus2,*) 'elliptic_flowBmes_y0.res'
      write(ifnus2,*) ndivpt1d
      do j=1,ndivpt1d
         pt=(j-0.5)*steppt1d
         write(ifnus2,151) pt,v2(9,j),Snum(9,j),Sden(9,j)
      enddo
      Write(ifnus2,*) '  '
      write(ifnus2,*) ' 34'
      write(ifnus2,*) 'elliptic_flowlept_y0.res'
      write(ifnus2,*) ndivpt1d
      do j=1,ndivpt1d
         pt=(j-0.5)*steppt1d
         write(ifnus2,151) pt,(v2(nres+il,j),il=1,4)
     &     ,(Snum(nres+il,j),il=1,4),
     &      (Sden(nres+il,j),il=1,4)
      enddo
      Write(ifnus2,*) '  '
      write(ifnus2,*) ' 35'
      write(ifnus2,*) 'elliptic_flowPsi_y0.res'
      write(ifnus2,*) ndivpt1d
      do j=1,ndivpt1d
         pt=(j-0.5)*steppt1d
         write(ifnus2,151) pt,v2(3,j),Snum(3,j),Sden(3,j)
      enddo
      Write(ifnus2,*) '  '
      write(ifnus2,*) '36 '
      write(ifnus2,*) 'elliptic_flow_nonpromptpsi_y0.res'
      write(ifnus2,*) ndivpt1d
      do j=1,ndivpt1d
         pt=(j-0.5)*steppt1d
         write(ifnus2,151) pt,v2nonpromptpsi(j),Snumnonpromptpsi(j),
     &                Sdennonpromptpsi(j)
      enddo
      write(ifnus2,*) '******* time evol *******  '
      Write(ifnus2,*) '  '
      write(ifnus2,*) '37 '
      write(ifnus2,*) 'nbpsi'
      write(ifnus2,*)ndivpsimax 
      do i=0,ndivpsimax
        write(ifnus2,*) i,packnbpsi(i)/nbsuccevent 
      enddo

      Write(ifnus2,*) '  '
      write(ifnus2,*) '38 dpsidy_of_time'
      Write(ifnus2,*) '  '
      write(ifnus2,*)ndivtmax 
      do itime=0,ndivtmax
         write(ifnus2,151) timesave(itime),
     &      (avdnbpsidyoftime(itime,j)/nbsuccevent,j=-5,5)
      enddo

      Write(ifnus2,*) '  '
      write(ifnus2,*) '39 pt-and radial distribution of cquarks '
     &   ,'and D-mesons as afunction of time'
      Write(ifnus2,*) '  '
      write(ifnus2,*)ndivpt1d,ndivrap1dmax,ndivtmax 
      do itime=0,ndivtmax
         do j=1,ndivpt1d
            if(packvacpty0(1,j).gt.0) then 
               ratio4=dncdptoftimey0(itime,j)/packvacpty0(1,j)
            else
               ratio4=0.
            endif
            pt=(j-0.5)*steppt1d

         enddo

      enddo
      do itime=0,ndivtmax
      goto 999 
	 Write(ifnus2,*) '  '
         write(ifnus2,*) '40 pt-y-distribution of cquarks as a ??,
     &    ??function of time'
	 Write(ifnus2,*) '  '
         write(ifnus2,*)ndivpt1d,ndivrap1dmax 
         do j=1,ndivpt1d
            do k=1,ndivrap1dmax
             if(dncdptdyvac(j,k).gt.0) then 
               ratio4=dncdptdyoftime(itime,j,k)/dncdptdyvac(j,k)
             else
               ratio4=0.
             endif
             pt=(j-0.5)*steppt1d
             rap=(k-0.5)*steprap1d+boundrapinf

            enddo
         enddo
      enddo
 999  close(ifnus2)
      return
      end


      subroutine displaydist
      use PBGgenvar
      implicit none
      integer partcount,ipart,nextpart,k,nhist,ndivz,ndivx,ndivt,ires
      real boundzinf,boundzsup,boundxinf,boundxsup,stepz,stepx,x(0:3),
     & boundpzinf,boundpzsup,boundpxinf,boundpxsup,steppz,steppx,
     & p(0:3),boundtinf,boundtsup,stept 
      integer ifmtx
      call getMonitorFileIndex(ifmtx)
      ndivx=20
      ndivz=20
      ndivt=20
      boundzinf=-6.
      boundzsup=6.
      boundxinf=-7.
      boundxsup=7.
      boundtinf=0.
      boundtsup=6.
      stepz=(boundzsup-boundzinf)/ndivz
      stepx=(boundxsup-boundxinf)/ndivx
      stept=(boundtsup-boundtinf)/ndivt
      boundpzinf=-6.
      boundpzsup=6.
      boundpxinf=-4.
      boundpxsup=4.
      steppz=(boundpzsup-boundpzinf)/ndivz
      steppx=(boundpxsup-boundpxinf)/ndivx
      nhist=1
      do ires=1,3
        ipart=firstpart
        partcount=0
         do while (ipart.NE.-1)
          partcount=partcount+1
          if(partinfo_ires(ipart).eq.ires) then
            do k=0,3
              x(k)=partinfo_r(ipart,k)
              p(k)=partinfo_p(ipart,k)
            enddo 
          endif
          nextpart=partinfo_next(ipart)
          if(partcount.GT.npart) then 
            write(ifmtx,*) 'part_chain_test failed'
            close(ifmtx)
            stop
          endif
          ipart=nextpart
         enddo
         if(npart.ne.partcount) then
            write(6,*) 'FAILURE in displaydist: npart=',npart,
     &      'partcount=',partcount
            stop
         endif
         nhist=nhist+5
      enddo
      return
      end

