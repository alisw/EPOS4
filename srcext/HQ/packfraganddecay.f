      subroutine initrepartBraat
      use PBGforquantDBraat
      implicit none
      integer nbquantmax,nbquant,iq,ipb
      parameter (nbquantmax=1000)
      double precision dz,tbquantDP(0:nbquantmax),quantileDBraatPc,
     & quantileDBraatVc,tbquantDV(0:nbquantmax),pb
      integer ifmtx
      call getMonitorFileIndex(ifmtx)
      nbquant=1000
      if(nbquant.gt.nbquantmax) then
         write(ifmtx,*) 'nbquant > ',nbquantmax,'in initrepartBraat'
         close(ifmtx)
         stop
      endif
      dz=1.d0/nbquant
      tbquantDP(0)=0.d0
      tbquantDV(0)=0.d0
      do iq=1,nbquant
         tbquantDP(iq)=tbquantDP(iq-1)+quantileDBraatPc((iq-1)*dz,iq*dz)
         tbquantDV(iq)=tbquantDV(iq-1)+quantileDBraatVc((iq-1)*dz,iq*dz)
      enddo
      do iq=1,nbquant
         tbquantDP(iq)=tbquantDP(iq)/tbquantDP(nbquant)
         tbquantDV(iq)=tbquantDV(iq)/tbquantDV(nbquant)
      enddo
      nbpb=100
      if(nbpb.gt.nbpbmax) then
         write(ifmtx,*) 'nbpb > ',nbpbmax,'in initrepartBraat'
         close(ifmtx)
         stop
      endif
      dpb=1.d0/nbpb
      zpbDP(0)=0.d0
      iq=0
      do ipb=1,nbpb-1
         pb=ipb*dpb
         do while(tbquantDP(iq).lt.pb)
            iq=iq+1
         enddo
         zpbDP(ipb)=(iq-1+(pb-tbquantDP(iq-1))/(tbquantDP(iq)
     &               -tbquantDP(iq-1)))
     &                *dz
      enddo         
      zpbDP(nbpb)=1
      zpbDV(0)=0.d0
      iq=0
      do ipb=1,nbpb-1
         pb=ipb*dpb
         do while(tbquantDV(iq).lt.pb)
            iq=iq+1
         enddo
         zpbDV(ipb)=(iq-1+(pb-tbquantDV(iq-1))/(tbquantDV(iq)
     &                   -tbquantDV(iq-1)))
     &                *dz
      enddo         
      zpbDV(nbpb)=1
      end

      subroutine initrepartKart
      use PBGforquantDkart
      implicit none
      integer nbquantmax,nbquant,iq,ipb
      parameter (nbquantmax=1000)
      double precision dz,tbquantDK(0:nbquantmax),Dkartrepart,alpha,pb
      integer ifmtx
      call getMonitorFileIndex(ifmtx)
      alpha=34.d0
      nbquant=1000
      if(nbquant.gt.nbquantmax) then
         write(ifmtx,*) 'nbquant > ',nbquantmax,'in initrepartKart'
         close(ifmtx)
         stop
      endif
      dz=1.d0/nbquant
      do iq=0,nbquant
         tbquantDK(iq)= Dkartrepart(alpha,iq*dz)
      enddo
      nbpb=200
      if(nbpb.gt.nbpbmax) then
         write(ifmtx,*) 'nbpb > ',nbpbmax,'in initrepartKart'
         close(ifmtx)
         stop
      endif
      dpb=1.d0/nbpb
      zpbDK(0)=0.d0
      iq=0
      do ipb=1,nbpb-1
         pb=ipb*dpb
         do while(tbquantDK(iq).lt.pb)
            iq=iq+1
         enddo
         zpbDK(ipb)=(iq-1+(pb-tbquantDK(iq-1))/(tbquantDK(iq)
     &               -tbquantDK(iq-1)))
     &                *dz
      enddo         
      zpbDK(nbpb)=1
      zpbDK(0)=zpbDK(1)
      end

      subroutine peterson(zmin,itypQ,z)
      use PBGfragmentc
      implicit none
      integer itypQ,nit,i
      parameter (nit=50)
      double precision zmin,z,ran2,dtry,maxval
      integer ifmtx
      call getMonitorFileIndex(ifmtx)
      if(zmin.lt.peterzmax(itypQ)) then
         maxval=petermaxval(itypQ)
      else
         maxval=zmin*(1-zmin)**2/(peterepsh(itypQ)+zmin**2-
     &     (1+peterepsh(itypQ)-peterepsq(itypQ))*zmin)**2
      endif
      do i=1,nit
         z=zmin+(1-zmin)*ran2()
         dtry=z*(1-z)**2/(peterepsh(itypQ)+z**2-(1+peterepsh(itypQ)-
     &           peterepsq(itypQ))*z)**2
         if(dtry.gt.maxval) then 
            write(ifmtx,*) 'unexpected high dtry in peterson;d(',
     &           z,')=',dtry
            write(ifmtx,*) itypQ,peterepsq(itypQ),peterepsh(itypQ)
            close(ifmtx)
            stop
         endif
         if(maxval*ran2().lt.dtry) goto 10
      enddo
      write(ifmtx,*) 'no success in peterson after ',nit,' iterations.',
     & zmin,
     & petermaxval(itypQ),maxval            
      do i=1,nit
         z=zmin+(1-zmin)*ran2()
         dtry=z*(1-z)**2/(peterepsh(itypQ)+z**2-(1+peterepsh(itypQ)-
     &           peterepsq(itypQ))*z)**2
         write(ifmtx,*) z,dtry
      enddo
      close(ifmtx)
      stop
 10   return
      end


      double precision function DBraatPc(z)
      implicit none
      double precision z
      DBraatPc=0.1*(1-z)**2*(6+(-14.4+(14.28+(-7.704+1.9926*z)*z)*z)*z)
     &   *z/
     &  (1-0.9*z)**6
      end

      double precision function DBraatVc(z)
      implicit none
      double precision z
      DBraatVc=0.3*(1-z)**2*(2+(-5.6+(8.52+(-7.056+ 2.2842*z)*z)*z)*z)
     &  *z/
     &  (1-0.9*z)**6
      end

      double precision function quantileDBraatPc(zinf,zsup)
      implicit none
      double precision zinf,zsup,DBraatPc,eps,ss
      parameter(eps=1D-8)
      external DBraatPc
      call qromb(DBraatPc,zinf,zsup,ss,eps)
      quantileDBraatPc=ss
      return
      end

      double precision function quantileDBraatVc(zinf,zsup)
      implicit none
      double precision zinf,zsup,eps,ss
      parameter(eps=1D-8)
      external DBraatVc
      call qromb(DBraatVc,zinf,zsup,ss,eps)
      quantileDBraatVc=ss
      return
      end

      double precision function zofDBraatP(pb)
      use PBGforquantDBraat
      implicit none
      integer ipb
      double precision pb,pbinf,pbsup
      ipb=int(pb/dpb)
      pbinf=ipb*dpb
      pbsup=pbinf+dpb
      zofDBraatP=((pb-pbinf)*zpbDP(ipb+1)+(pbsup-pb)*zpbDP(ipb))/dpb
      end

      double precision function zofDBraatV(pb)
      use PBGforquantDBraat
      implicit none
      integer ipb
      double precision pb,pbinf,pbsup
      ipb=int(pb/dpb)
      pbinf=ipb*dpb
      pbsup=pbinf+dpb
      zofDBraatV=((pb-pbinf)*zpbDV(ipb+1)+(pbsup-pb)*zpbDV(ipb))/dpb
      end

      double precision function DKart(z,alpha)
      implicit none
      double precision z,alpha
      DKart=(1+alpha)*(2+alpha)*(1-z)*z**alpha
      end

      double precision function Dkartrepart(alpha,zsup)
      implicit none
      double precision zsup,alpha
      Dkartrepart=(1+(alpha+1)*(1-zsup))*zsup**(1+alpha)
      end

      double precision function zofDKart(pb)
      use PBGforquantDkart
      implicit none
      integer ipb
      double precision pb,pbinf,pbsup
      ipb=int(pb/dpb)
      pbinf=ipb*dpb
      pbsup=pbinf+dpb
      zofDKart=((pb-pbinf)*zpbDK(ipb+1)+(pbsup-pb)*zpbDK(ipb))/dpb
      end

      double precision function zforqfragCacciari(itypq)
      implicit none
      integer itypq
      double precision ran2,fracDP,zofDBraatP,zofDBraatV,
     &        zofDKart
      parameter(fracDP=0.41772233)
      integer ifmtx
      call getMonitorFileIndex(ifmtx)
      if(itypq.eq.1) then
         if(ran2().lt.fracDP) then
            zforqfragCacciari=zofDBraatP(ran2())
         else
            zforqfragCacciari=zofDBraatV(ran2())*0.93
         endif
      else if(itypq.eq.2) then
         zforqfragCacciari=zofDKart(ran2())
      else
         write(ifmtx,*) 'itypq unrecog in zforqfrag:',itypq
         close(ifmtx)
         stop
      endif
      end
         
      subroutine fragment(pq,iresq,ph,iresh)
      use PBGfragandcoal
      use PBGgenvar
      implicit none
      integer iresq,iresh,i
      double precision pq(0:3),ph(0:3),pqnorm2,pqnorm,pqplus,z,rap,pnew,
     &   pnew2,
     &  zforqfragCacciari
      integer ifmtx
      call getMonitorFileIndex(ifmtx)
      if((itypfrag.eq.1).or.(itypfrag.eq.2)) then
         pqnorm2=pq(1)**2+pq(2)**2+pq(3)**2
         pqnorm=sqrt(pqnorm2)
         pqplus=sqrt(pqnorm2+resinfo_mass(iresq)**2)+pqnorm
         if(itypfrag.eq.1) then
            z=1
         endif
         if(itypfrag.eq.2) then
            call peterson(resinfo_mass(iresh)/pqplus,itypquark(iresq),z)
         endif
         pnew=(z*pqplus-resinfo_mass(iresh)**2/(z*pqplus))/2
         if(pnew.lt.0) then
            write(ifmtx,*) '|p|^2<0 in fragment:',pnew
            write(ifmtx,*) pq(0),pq(1),pq(2),pq(3)
            close(ifmtx)
            stop
         endif
         pnew2=pnew**2
         rap=pnew/pqnorm
      else if(itypfrag.eq.3) then
         rap=zforqfragCacciari(itypquark(iresq))
         if(rap>1.0d0) then
           write(ifmtx,*) 'rap larger than 1 in fragment,rap= ',rap
           write(ifmtx,*) pq(0),pq(1),pq(2),pq(3)
           close(ifmtx)
           stop
         endif
       else
         write(ifmtx,*) 'itypfrag unrecognized in routine fragment'
         close(ifmtx)
         stop
      endif
      do i=1,3
         ph(i)=rap*pq(i)
         pnew2=ph(1)**2+ph(2)**2+ph(3)**2
      enddo
      ph(0)=sqrt(resinfo_mass(iresh)**2+pnew2)
      return
      end

      subroutine decayBtoDsemilepto(p4B,p4D)
      use PBGopenflavor
      implicit none
      logical accept,first
      integer i,it,nitmax
      parameter(nitmax=1000)
      double precision p4B(0:3),p4D(0:3),dw,dwsup,plsupd,pdcm,
     &  pd4cm(0:3),
     & u4B(0:3),cth,pt,phi,ran2,pi,rej
      parameter (pi=3.14159265358979323844)
      data first/.true./
      save first,plsupd,dwsup
      integer ifmtx
      call getMonitorFileIndex(ifmtx)
      if(first) then
         plsupd=(mbmes**2-md**2)/(2*mbmes)
         dwsup=sqrt((plsupd/md)**2+1)-1
         first=.false.
      endif
      accept=.false.
      it=0
      do while((.not.accept).and.(it.lt.nitmax))
         it=it+1
         if(it.ge.nitmax) then
            write(ifmtx,*) 'decay out of loop in decayBtoDsemilepto'
            write(ifmtx,*) (p4B(i),i=0,3)
            close(ifmtx)
            stop
         endif
         dw=ran2()*dwsup
         rej=0.537675d0*(1-1.22015*dw)/(1+1.36075*dw+1.25932*dw**2)
         pdcm=sqrt(dw*(2+dw))
         rej=rej*pdcm*(2+dw)**2
         if(rej-1.gt.1D-5) then
            write(ifmtx,*) 'rej > 1 in decayBtoD'
            close(ifmtx)
            stop
         endif
         accept=ran2().lt.rej
      enddo
      pdcm=pdcm*md
      cth=2*ran2()-1
      pt=sqrt(1-cth**2)*pdcm
      phi=2*pi*ran2() 
      pd4cm(0)=(dw+1)*md
      pd4cm(1)=pt*cos(phi) 
      pd4cm(2)=pt*sin(phi)
      pd4cm(3)=pdcm*cth
      do i=0,3
         u4B(i)=p4B(i)/mBmes
      enddo
      call boost(u4B,pd4cm,p4d)
      return
      end

      subroutine weakdecayD(pD,pl)
      use PBGopenflavor
      implicit none
      logical accept
      integer i,it,nitmax
      parameter(nitmax=1000)
      double precision pD(0:3),pl(0:3),uD(0:3),plcm,pl4cm(0:3),
     &  cth,pt,phi,ran2,pi,plsupe,rej
      parameter (pi=3.14159265358979323844,plsupe=1.)
      integer ifmtx
      call getMonitorFileIndex(ifmtx)
      accept=.false.
      it=0
      do while((.not.accept).and.(it.lt.nitmax))
         it=it+1
         if(it.ge.nitmax) then
            write(ifmtx,*) 'decay out of loop in weakdecayD'
            write(ifmtx,*) (pD(i),i=0,3)
            close(ifmtx)
            stop
         endif
         plcm=ran2()*plsupe
         rej=44.2365d0*plcm**2.5d0*(plsupe-plcm)**3
         accept=ran2().lt.rej
      enddo
      cth=2*ran2()-1
      pt=sqrt(1-cth**2)*plcm
      phi=2*pi*ran2() 
      pl4cm(0)=plcm 
      pl4cm(1)=pt*cos(phi) 
      pl4cm(2)=pt*sin(phi)
      pl4cm(3)=plcm*cth
      do i=0,3
         uD(i)=pD(i)/mD
      enddo
      call boost(uD,pl4cm,pl)
      return
      end

      subroutine weakdecayB(pB,pl)
      use PBGopenflavor
      implicit none
      logical accept
      integer i,it,nitmax
      parameter(nitmax=1000)
      double precision pB(0:3),pl(0:3),uB(0:3),plcm,pl4cm(0:3),
     &  cth,pt,phi,ran2,pi,plsupe,rej
      parameter (pi=3.14159265358979323844)
      integer ifmtx
      call getMonitorFileIndex(ifmtx)
      plsupe=0.5d0*mBmes
      accept=.false.
      it=0
      do while((.not.accept).and.(it.lt.nitmax))
         it=it+1
         if(it.ge.nitmax) then
            write(ifmtx,*) 'decay out of loop in weakdecayB'
            write(ifmtx,*) mD,mBmes,plcm
            write(ifmtx,*) (pB(i),i=0,3)
            close(ifmtx)
            stop
         endif
         plcm=ran2()*plsupe
         rej=1.01146*(1-plcm/plsupe)*plcm**2/(1+exp(9.0984d0
     &     *(plcm-2.0237d0))) 
         if(rej-1.gt.1D-5) then
            write(ifmtx,*) 'rej > 1 in weakdecayB'
            close(ifmtx)
            stop
         endif
         accept=ran2().lt.rej
      enddo
      cth=2*ran2()-1
      pt=sqrt(1-cth**2)*plcm
      phi=2*pi*ran2() 
      pl4cm(0)=plcm 
      pl4cm(1)=pt*cos(phi) 
      pl4cm(2)=pt*sin(phi)
      pl4cm(3)=plcm*cth
      do i=0,3
         uB(i)=pB(i)/mBmes
      enddo
      call boost(uB,pl4cm,pl)
      return
      end

      subroutine weakdecayBsecond(p4B,p4l)
      implicit none
      double precision p4B(0:3),p4D(0:3),p4l(0:3)
      call decayBtoDsemilepto(p4B,p4D)
      call weakdecayD(p4D,p4l)
      end

      double precision function pstarjpsifrombdecay()
      implicit none      
      double precision x,ran2
      double precision xmin,xmax,xrange
      parameter (xmin=0.049220087321580966d0,xmax=0.9752654483640086d0)
      parameter (xrange=xmax-xmin)
      integer ix,nbinx
      parameter (nbinx=200)
      double precision dBdx(nbinx)
      logical first
      data first/.true./
      save first,dBdx
      double precision dx
      parameter (dx=xrange/nbinx)
      double precision top
      integer ifmtx
      call getMonitorFileIndex(ifmtx)
      if(first) then
         top=0.d0
         do ix=1,nbinx
            x=(ix-0.5d0)*dx+xmin
            dBdx(ix)=0.000050517d0-0.00379d0*x+0.08878d0*x**2
     &           -0.85901d0*x**3
     &           +4.62988d0*x**4 - 14.35995*x**5 + 26.05286*x**6
     &           -27.20092*x**7 + 15.10314*x**8 - 3.45121*x**9
            if(dBdx(ix).lt.0) then
               write(ifmtx,*) 'dBdx smaller then 0 found'
               write(ifmtx,*) x
               close(ifmtx)
               stop
            endif
            if(dBdx(ix).gt.top) top=dBdx(ix)
         enddo
         do ix=1,nbinx
            dBdx(ix)=dBdx(ix)/top
         enddo
         first=.false.
      endif
 10   continue   
      x=xrange*ran2()
      ix=int(x/dx)+1
      if(ix.eq.0.) goto 10
      if(ix.gt.nbinx) goto 10
      if(ran2().gt.dBdx(ix)) goto 10   
      pstarjpsifrombdecay=(x+xmin)*2
      end


      subroutine BdecayJPsi(pB,p4JPsi)
      use PBGopenflavor
      implicit none
      double precision mpsivacuum
      parameter (mpsivacuum=3.097d0)
      integer i
      double precision pB(0:3),p4JPsi(0:3),uB(0:3),p4jpsicm(0:3),
     &  cth,pt,phi,ran2,pi,pstarjpsifrombdecay,pstar
      parameter (pi=3.14159265358979323844)
      pstar=pstarjpsifrombdecay()
      cth=2*ran2()-1
      pt=sqrt(1-cth**2)*pstar
      phi=2*pi*ran2() 
      p4jpsicm(0)=sqrt(pstar**2+mpsivacuum**2) 
      p4jpsicm(1)=pt*cos(phi) 
      p4jpsicm(2)=pt*sin(phi)
      p4jpsicm(3)=pstar*cth
      do i=0,3
         uB(i)=pB(i)/mBmes
      enddo
      call boost(uB,p4jpsicm,p4JPsi)
      return
      end
         
