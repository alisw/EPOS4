
      double precision function zlongmaxoner(t)
      use PBGsceplas
      use PBGhydro
      use PBGpartinmedium
      implicit none

      integer ix,iy,imed,ifail
      double precision x(0:3),temp,beta(3),densener,densenermax
      double precision t,r,a,bb
      double precision zza,zzb,z_dist,distzz,distz
      logical cont

      if((isce.eq.1).or.(isce.eq.2)) then
        a=ri0+vr0*t+a*r/2.*t**2
        b=li0/2+vl0*t
        zlongmaxoner=bb*sqrt(1-(r/a)**2)
        return
      endif
      if((isce.eq.3).or.(isce.eq.4).or.(isce.eq.5).or.(isce.eq.6).or.
     & (isce.eq.7).or.(isce.eq.8).or.(isce.eq.9).or.(isce.eq.10).or.
     & (isce.eq.11).or.(isce.eq.12).or.(isce.eq.13))then
         zlongmaxoner=li0/2+vl0*t
         return
      endif
      if(isce.eq.7)then
         cont=.true.
         x(0)=t
         zza=zmin
         densenermax=0.0d0
         do while(cont)
            x(3)=t*dsinh(zza)/dcosh(zza)
            do ix=1,nxhyx
               x(1)=(ix-13)*xstep
               do iy=0,nyhyx
                  x(2)=(iy-13)*ystep
                  call givelocalhydro(x,temp,beta,densener,ifail,imed)
                  if((ifail.eq.0).and.(densener.gt.densenermax)) then
                  densenermax=densener
                  endif
               enddo
            enddo
            if((densenermax.ge.entrans_max).or.(zza.gt.zmax))then
               cont=.false. 
            else
               zza=zza+zstep
            end if
         enddo
         distz=x(3)
         cont=.true.
         zzb=zmax
         do while(cont)
            x(3)=t*dsinh(zzb)/dcosh(zzb)
            do ix=1,nxhyx
               x(1)=(ix-13)*xstep
               do iy=0,nyhyx
                  x(2)=(iy-13)*ystep
                  call givelocalhydro(x,temp,beta,densener,ifail,imed)
                  if((ifail.eq.0).and.(densener.gt.densenermax)) then
                  densenermax=densener
                  endif
               enddo
            enddo
            if((densenermax.ge.entrans_max).or.(zzb.lt.zmin))then
               cont=.false. 
            else
               zzb=zzb-zstep
            end if
         enddo
         distzz=x(3)
         if((zza.lt.zmax).and.(zzb.gt.zmin)) then
            z_dist=distzz-distz
         else
            z_dist=0.0d0
         end if
         zlongmaxoner=z_dist
         return
      endif
      write(6,*) 'scenario not recognized in zlongmaxoner'
      stop
      end


      double precision function radius(A)
      integer A
      double precision onethird
      onethird=0.3333333333333333d0
      radius=1.12d0*dble(A)**(onethird)-0.86d0*dble(A)**(-onethird)
      return
      end

      double precision function woodsaxondistr(x,y,z,A)
      integer A
      double precision x,y,z,radius,rho0,d
      rho0=0.1691096d0
      d=0.54d0
      woodsaxondistr=rho0/(1+dexp((dsqrt(x**2+y**2+z**2)-radius(A))/d))
      return
      end

      double precision function nuclthick(x,y,A)
      implicit none
      integer A,n
      double precision x,y,c,d,s,woodsaxondistr,r2,r
      external woodsaxondistr
      if(A.eq.197) then
         r2=x*x+y*y
         if(r2.gt.81) then
            nuclthick=0
         else
            nuclthick=(r2*(r2*0.000182074-0.0267047)+1)/
     &        (r2*(r2*(r2*0.00000474052-0.0000709273)-0.00564324)
     &         +0.463348)
         endif
      else if(A.eq.208) then
         r2=x*x+y*y
         r=sqrt(r2)
         nuclthick=2.2204*exp(((-0.000524949*r+0.00525874)*r-
     &      0.0172032)*r2/
     &     (((-0.000544998*r+0.0246878)*r-0.268878)*r+1))
      else
         write(*,*) 'aaaaah!'
         c=-20d0
         d=20d0
         n=200
         call intrapz(woodsaxondistr,x,y,A,c,d,s,n)
         nuclthick=s
      endif
      return
      end

      double precision function nuclthick_x(x,A)
      integer A,n
      double precision x,c,d,s,nuclthick
      external nuclthick
      c=-20d0
      d=20d0
      n=200
      call intrapy(nuclthick,x,A,c,d,s,n)
      nuclthick_x=s
      return
      end


      double precision function transvprobcoll(x,y,b,A)
      implicit none
      integer A
      real b
      double precision x,y,nuclthick,bb

      bb=dble(b)
      transvprobcoll=nuclthick(x+bb/2.0d0,y,A)*nuclthick(x-bb/2.0d0,y,A)
      return
      end


      double precision function denspart(x,y,b,A)
      integer A
      real b
      double precision x,y,nuclthick,bb
      bb=dble(b)
      denspart=nuclthick(x+bb/2.0d0,y,A)*
     &         (1.0d0-exp(-5.0d0*nuclthick(x-bb/2.0d0,y,A)))+
     &         nuclthick(x-bb/2.0d0,y,A)*
     &         (1.0d0-exp(-5.0d0*nuclthick(x+bb/2.0d0,y,A)))
      return
      end





      subroutine azget1x(x,b,A,tauf)
      use PBGvariancesig
      implicit none
      integer A
      real b
      double precision ran2,x(0:3),tauf,zlongmaxoner,x1,x2,
     &     transvprobcoll,twoDgaussian
      x1=0
      x2=0
      x(0)=tauf
      call twoDnormaldeviates(sigmax,sigmay,sigma,x1,x2)
      do while(ran2().gt.(transvprobcoll(x1,x2,b,A)/twoDgaussian(sigmax,
     &     sigmay,sigma,x1,x2)))
         call twoDnormaldeviates(sigmax,sigmay,sigma,x1,x2)
      enddo
      x(1)=x1
      x(2)=x2
      x(3)=(2.0d0*ran2()-1.0d0)*zlongmaxoner(x(0))      
      return
      end  


      subroutine azget1xnpart(x,b,A,tauf)
      implicit none
      integer A,i,imax
      parameter(imax=50)
      real b
      double precision ran2,x(0:3),tauf,zlongmaxoner,radius,r,xmax,ymax,
     &  xtry,ytry,denspart,densmax,bb
      bb=dble(b)
      r=radius(A)+3*0.54
      xmax=r-bb/2.0d0
      ymax=sqrt(r*r-bb*bb/4.0d0)
      densmax=denspart(0.d0,0.d0,b,A)
      do i=1,imax
         xtry=(2.0d0*ran2()-1.0d0)*xmax
         ytry=(2.0d0*ran2()-1.0d0)*ymax
         if(ran2().lt.denspart(xtry,ytry,b,A)/densmax) goto 10
      enddo
      write(6,*) 'problem in azget1xpart:'
      write(6,*) 'was not able to generate y within ',imax,'tries'
      write(6,*) b
      do i=1,imax
         xtry=(2.0d0*ran2()-1.0d0)*xmax
         ytry=(2.0d0*ran2()-1.0d0)*ymax
         write(6,*) xtry,ytry,denspart(xtry,ytry,b,A)/densmax
      enddo
      stop
 10   x(0)=tauf
      x(1)=xtry
      x(2)=ytry
      x(3)=(2.0d0*ran2()-1.0d0)*zlongmaxoner(x(0))      
      return
      end  


      double precision function RdAuJpsiGdC(rap,b)
      implicit none
      double precision rap,b
      if(b.gt.7) then
         RdAuJpsiGdC=1
      else
         RdAuJpsiGdC=(1-b/7.)*(0.4+0.6/(1+exp(rap+0.6)/0.4))+b/7.
      endif
      return
      end
      

      double precision function RAuAunucJPsiGdC(x,y,b,rap)
      use PBGsystem
      implicit none
      double precision x,y,rap,RdAuJpsiGdC
      real b
      if(A.ne.197) then
         write(8,*) 'unappropriate A in RAuAunucJPsiGdC'
         close(8)
         stop
      endif
      RAuAunucJPsiGdC=RdAuJpsiGdC(rap,sqrt((x+b/2)**2+y**2))*
     &                RdAuJpsiGdC(-rap,sqrt((x-b/2)**2+y**2))
      return
      end


      double precision function shadowpsiLHC_PBPB_ofy(y)
      implicit none
      double precision y,y2
      y2=y**2
      if(abs(y).lt.4.937) then
         shadowpsiLHC_PBPB_ofy=exp(y2*(y2*(y2*(y2*4.86752E-6-
     &                          0.0000227718)
     &                         -0.00110302)+0.0204211)-0.538909)
      else
         shadowpsiLHC_PBPB_ofy=2
      endif
      return
      end

      double precision function shadowpsiLHC_PBPB_ofpt(pt)
      implicit none
      double precision pt
      if(pt.gt.30) then
          shadowpsiLHC_PBPB_ofpt=1.
       else
          shadowpsiLHC_PBPB_ofpt=(0.65+0.35*tanh(pt/10))
       endif
      return
      end

      double precision function shadowpsiLHC_PBPB_ofp(y)
      implicit none
      double precision y,shadowpsiLHC_PBPB_ofy
      shadowpsiLHC_PBPB_ofp=shadowpsiLHC_PBPB_ofy(y)
      return
      end

      double precision function maxantishadowingPsi()
      use PBGsystem
      implicit none
      double precision shadowpsiLHC_PBPB_ofp
      maxantishadowingPsi=0.
      if(abs(sqrtS-200.).lt.0.1) then
         maxantishadowingPsi=1.
      else if(abs(sqrtS-5250).lt.260) then
         if(A.eq.208) then
            maxantishadowingPsi=shadowpsiLHC_PBPB_ofp(10.d0)
         else
            call utstop('maxantishadowingPsi not defined for A<>208&')
         endif
       else
         call utstop('maxantishadowingPsi not defined for this E&')
      endif
      end


      logical function rejectnuclear(x,b,p)
      use PBGsystem
      implicit none
      double precision x(0:3),p(0:3),rap,RAuAunucJPsiGdC,ran2,
     & shadowpsiLHC_PBPB_ofp,maxantishadowingPsi
      real b
      rap=0.5*log((p(0)+p(3))/(p(0)-p(3)))
      if(abs(sqrtS-200).lt.0.1) then
         rejectnuclear=ran2().gt.RAuAunucJPsiGdC(x(1),x(2),b,rap)
      else if(abs(sqrtS-5250).lt.260) then
         rejectnuclear=ran2().gt.(shadowpsiLHC_PBPB_ofp(rap)
     &   /maxantishadowingPsi())
      else
         write(8,*) 'system not foreseen in rejectnuclear', A,sqrtS 
         close(8)
         stop
      endif
      return
      end
         

      subroutine azget2x(itypQ,xcm,x1,x2,b,A,tauf,iprofile)
      use PBGvariancesig
      implicit none
      integer A,iprofile,itypQ
      real b
      double precision pi,ran2,x1(0:3),x2(0:3),tauf,xcm(0:3),rrel,
     & thetarel, xrel,yrel,zlongmaxoner
      parameter(pi=3.14159265358979323844d0)
      if(iprofile.eq.1) then
         call azget1x(xcm,b,A,tauf)
      else
         call azget1xnpart(xcm,b,A,tauf)
      endif
      x1(0)=xcm(0) 
      x2(0)=xcm(0)
      if(itypQ.eq.1) then
         rrel=0.15d0*sqrt(ran2())
      else
         rrel=0.05d0*sqrt(ran2())
      endif
      thetarel=2.0d0*pi*ran2()
      xrel=rrel*dcos(thetarel)
      yrel=rrel*dsin(thetarel)
      x1(1)=xcm(1)+xrel
      x2(1)=xcm(1)-xrel
      x1(2)=xcm(2)+yrel
      x2(2)=xcm(2)-yrel
      x1(3)=(2.0d0*ran2()-1d0)*zlongmaxoner(x1(0))  
      x2(3)=(2.0d0*ran2()-1d0)*zlongmaxoner(x2(0))  
      return
      end


      subroutine getonep0(p,mass)
      use PBGsceplas
      use PBGpartinmedium
      implicit none
      integer k
      double precision pi,meansqrtp3
      parameter (pi=3.14159265358979323844d0)
      parameter (meansqrtp3=1.)
      double precision p(0:3),mass,bid
       call gauss(p(1),p(2))
       call gauss(p(3),bid)
       p(0)=mass**2 
       do k=1,3
          p(k)=p(k)*meansqrtp3
          p(0)=p(0)+p(k)**2
       enddo
       p(0)=sqrt(p(0))
      return
      end  

      double precision function distrib_c_rhic_pt(pt,ifact)
      implicit none
      integer ifact
      double precision pt,fact,bas
      bas=1./(1+(pt/2.20429)**2)**3.92985
      fact=1.
      if(ifact.eq.1) then
         fact=0.03551441
      endif
      if(ifact.eq.2) then
         fact=5.40451
      endif
      distrib_c_rhic_pt=fact*bas
      return
      end

      double precision function distrib_D_rhic_pt(pt,ifact)
      implicit none
      integer ifact
      double precision pt,fact,bas
      bas=1./(1+(pt/1.646)**2)**3.7656
      fact=1.
      if(ifact.eq.1) then
         fact=0.071
      endif
      if(ifact.eq.2) then
         fact=4.57386
      endif
      distrib_D_rhic_pt=fact*bas
      return
      end

      double precision function yrand_c(distry,ymax)
      implicit none
      integer imax,i,distry
      double precision pi,y,ytry2,ran2,bid,ymax
      parameter (imax=20)
      parameter (pi=3.14159265358979323844d0)
      if(distry.eq.1) then
         y=0.
      else if(distry.eq.2) then
         y=ymax*(2*ran2()-1)
      else if(distry.eq.3) then
         do i=1,imax
            bid=cos(pi*ran2())
            ytry2=-8.8628*Log(ran2())*bid*bid
            if(ran2().lt.exp(-ytry2*ytry2*(0.00574647+
     &           0.000240332*ytry2))) goto 10
         enddo
         write(6,*) 'problem in yrand_c:'
         write(6,*) 'was not able to generate y within ',imax,'tries'
         stop
 10      continue
         y=sign(sqrt(ytry2),bid)
      endif
      yrand_c=y
      return
      end


      subroutine getonept_c(mass,ptrans,mt2)
      use PBGsystem
      use PBGspectra
      implicit none
      double precision ptrans(2),mass,pi,ran2,pt2,pt,phi,mt2,temp,
     & get1ptfromspec
      parameter (temp=0.25)
      parameter (pi=3.14159265358979323844d0)
      integer ifmtx
      call getMonitorFileIndex(ifmtx)
      if(ifextspectra(1)) then
         pt=get1ptfromspec(1,ifshadowing(1))
         pt2=pt*pt
      else
         if(abs(sqrtS-200).lt.0.1) then
            pt2=4.85892*(ran2()**(-0.3413)-1)
         else if(abs(sqrtS-5250).lt.260) then
            pt2=10.583*(ran2()**(-0.4842)-1)
         else if(abs(sqrtS-2760).lt.0.1) then
            pt2=10.06*(ran2()**(-0.44677)-1)
         else
            write(ifmtx,*) 'HQ:unrecognized sqrt(s) in getonept_c'
            close(ifmtx)
            stop
         endif
         pt=sqrt(pt2)
      endif
      mt2=pt2+mass*mass
      phi=2*pi*ran2()
      ptrans(1)=pt*cos(phi)
      ptrans(2)=pt*sin(phi)
      return
      end

      double precision function yrand_b(distry,ymax)
      implicit none
      integer imax,i,distry
      double precision pi,y,ytry2,ran2,bid,ymax
      parameter (imax=20)
      parameter (pi=3.14159265358979323844d0)
      if(distry.eq.1) then
         y=0.
      else if(distry.eq.2) then
         y=ymax*(2*ran2()-1)
      else if(distry.eq.3) then
         do i=1,imax
            bid=cos(pi*ran2())
            ytry2=-2.98191*Log(ran2())*bid*bid
            if(ran2().lt.exp(-ytry2*ytry2*(0.00137026+
     &           0.00368154*ytry2))) goto 10
         enddo
         write(6,*) 'problem in yrand_b:'
         write(6,*) 'was not able to generate y within ',imax,'tries'
         stop
 10      continue
         y=sign(sqrt(ytry2),bid)
      endif
      yrand_b=y
      return
      end

      subroutine getonept_b(mass,ptrans,mt2)
      use PBGsystem
      use PBGspectra
      implicit none
      double precision ptrans(2),mass,pi,ran2,pt2,pt,phi,mt2,
     &  get1ptfromspec
      parameter (pi=3.14159265358979323844d0)
      integer ifmtx
      call getMonitorFileIndex(ifmtx)
      if(ifextspectra(2)) then 
         pt=get1ptfromspec(2,ifshadowing(2))
         pt2=pt*pt
      else
         if(abs(sqrtS-200).lt.0.1) then
            pt2=35.4774*(ran2()**(-0.30813)-1)
         else if(abs(sqrtS-5250).lt.260) then
            pt2=39.694*(ran2()**(-0.5572)-1) 
         else if(abs(sqrtS-2760).lt.0.1) then
            pt2=45.8*(ran2()**(-0.47229)-1)
         else
            write(ifmtx,*) 'unrecognized sqrt(s) in getonept_b'
            close(ifmtx)
            stop
         endif
         pt=sqrt(pt2)
      endif
      mt2=pt2+mass*mass
      phi=2*pi*ran2()
      ptrans(1)=pt*cos(phi)
      ptrans(2)=pt*sin(phi)
      return
      end

      subroutine getonep_Q(itypq,mass,distry,ymax,p,mt2)
      implicit none
      integer distry,itypq
      double precision ymax,p(0:3),mass,y,ptrans(2),mt2,yrand_c,yrand_b
      integer ifmtx
      call getMonitorFileIndex(ifmtx)
      if(itypq.eq.1) then
         y=yrand_c(distry,ymax)
         call getonept_c(mass,ptrans,mt2)
      else if(itypq.eq.2) then
         y=yrand_b(distry,ymax)
         call getonept_b(mass,ptrans,mt2)
      else
         write(ifmtx,*) 'HQ:itypq unrecognized in getonep_Q:',itypq
         close(ifmtx)
         stop
      endif
      p(1)=ptrans(1)
      p(2)=ptrans(2)
      p(3)=sqrt(mt2)*sinh(y)
      p(0)=sqrt(mt2+p(3)**2)         
      return
      end


      subroutine get2p_q(itypq,p1,p2,mass)
      use PBGsystem
      use PBGdistccbar
      use forydistr
      implicit none
      double precision p1(0:3),p2(0:3),mass,ran2,pi
      parameter(pi=3.14159265358979323844d0)
      double precision ptc,ptcbar,dphi,mtq2,mtc2,mtqbar2,mtcbar2,phi1
     & ,phi2,yqbar,yc,ycbar,yrand_c,yrand_b,ymax,p1T,p2T
      integer itypq,irec,distry
      integer ifmtx
      call getMonitorFileIndex(ifmtx)
 12   format(e14.5,e14.5,e14.5,e14.5,e14.5)
      if(abs(sqrts-200).lt.0.1) then 
         if(ifyat0) then
            distry=1
         else
            distry=3
         endif
         ymax=4
      else if((abs(sqrts-5250).lt.260).or.(abs(sqrts-2760).lt.0.1)) then
         if(ifyat0) then
            distry=1
         else
            distry=2
         endif
         ymax=4
      else
         write(ifmtx,*) 'HQ:unknown system found in get2p_Q'
         close(ifmtx)
         stop
      endif 
      if(idist.eq.1) then
         call getonep_Q(itypq,mass,distry,ymax,p1,mtq2)
         p2(1)=-p1(1)
         p2(2)=-p1(2)
         if(itypq.eq.1) then
            yqbar=yrand_c(distry,ymax)
         else if(itypq.eq.2) then
            yqbar=yrand_b(distry,ymax)
         else
            write(ifmtx,*) 'HQ:itypq not recognized in get2p_Q:',itypq
            close(ifmtx)
            stop
         endif
         p2(3)=sqrt(mtq2)*sinh(yqbar)
         p2(0)=sqrt(mtq2+p2(3)**2)      
      endif
      if(idist.eq.2) then
         call getonep_Q(itypq,mass,distry,ymax,p1,mtq2)
         call getonep_Q(itypq,mass,distry,ymax,p2,mtqbar2)
      endif
      if(idist.eq.3) then
         if(itypq.eq.1) then
            irec=1+int(ran2()*nbccbarindatabase)
            read(11,12,rec=irec) yc,ycbar,ptc,ptcbar,dphi
            mtc2=mass**2+ptc**2
            mtcbar2=mass**2+ptcbar**2
            p1(3)=sqrt(mtc2)*sinh(yc)
            p1(0)=sqrt(mtc2+p1(3)**2)
            p2(3)=sqrt(mtcbar2)*sinh(ycbar)
            p2(0)=sqrt(mtcbar2+p2(3)**2)
            phi2=2*pi*ran2()
            phi1=phi2+dphi
            p1(1)=ptc*cos(phi1)
            p1(2)=ptc*sin(phi1)
            p2(1)=ptcbar*cos(phi2)
            p2(2)=ptcbar*sin(phi2)
         else
            write(ifmtx,*) 'HQ:no database for b in get2p_q'
            close(ifmtx)
            stop
         endif
      endif
      if(idist.eq.4) then
         call getonep_Q(itypq,mass,distry,ymax,p1,mtq2)
         call getonep_Q(itypq,mass,distry,ymax,p2,mtqbar2)
         p1T=sqrt(p1(1)**2+p1(2)**2)
         p2T=sqrt(p2(1)**2+p2(2)**2)
         if(p1T.gt.0.d0) then 
            p2(1)=p1(1)*p2T/p1T
            p2(2)=p1(2)*p2T/p1T
         endif
      endif
      return
      end

      double precision function distrib_psi_rhic_pt(pt,ifact)
      implicit none
      integer ifact
      double precision pt,fact,bas
      bas=1./(1+(pt/3.94)**2)**6
      fact=1.
      if(ifact.eq.1) then
         fact=6.8D-5
      endif
      if(ifact.eq.2) then
         fact=0.102525
      endif
      distrib_psi_rhic_pt=fact*bas
      return
      end


      subroutine getonep_psi(p,mass)
      use PBGsystem
      implicit none
      integer imaxRHIC,imaxLHC,i
      parameter(imaxRHIC=20,imaxLHC=30)
      double precision p(0:3),mass,pt2,mt2,pt,ran2,pi,phi,y,ytry2,bid
      parameter (pi=3.14159265358979323844d0)      
      integer ifmtx
      call getMonitorFileIndex(ifmtx)
      if(abs(sqrtS-200).lt.0.1) then
         do i=1,imaxRHIC
            bid=cos(pi*ran2())
            ytry2=-Log(ran2())*bid*bid/0.168
            if(ran2().lt.exp(-ytry2*ytry2*(0.0007248
     &         +0.0008496*ytry2))) goto 10
         enddo
         write(ifmtx,*) 'HQ:problem in getonep_psi:'
         write(ifmtx,*) 'unable to generate ywithin ',imaxRHIC,'tries'
         close(ifmtx)
         stop
 10      continue
         pt2=15.5236*(ran2()**(-0.2)-1)
      else if(abs(sqrts-5250).lt.260) then
         do i=1,imaxLHC
            bid=cos(pi*ran2())
            ytry2=-log(ran2())*bid*bid/0.0107174
            if(ran2().lt.exp(-ytry2*ytry2*(0.0030485981
     &                          -3.476542592E-4*ytry2+ 
     &                          1.454073E-5*ytry2**2))) goto 15
         enddo
         write(ifmtx,*) 'problem in getonep_psi:'
         write(ifmtx,*) 'unable to generate y within ',imaxLHC,'tries'
         close(ifmtx)
         stop
 15      continue
         pt2=(ran2()**(-1./(2.55095-1.))-1)/0.3629186
      else
         write(8,*) 'only able to generate Psi for RHIC (at 200 GeV)'
         write(8,*) 'and for LHC at 5.5 TeV'
         close(8)
         stop
      endif 
      y=sign(sqrt(ytry2),bid)
      pt=sqrt(pt2)
      mt2=mass**2+pt2
      phi=2*pi*ran2()
      p(1)=pt*cos(phi)
      p(2)=pt*sin(phi)
      p(3)=sqrt(mt2)*sinh(y)
      p(0)=sqrt(mt2+p(3)**2)      
      return
      end


      subroutine get_ktkick_pA(A,b,kt)
      implicit none
      integer A
      double precision b,kt(0:3),avnbcoll,b2,delta0,ran2,kt2,theta,pi
      parameter(delta0=0.2,pi=3.14159265358979323844d0)      

      b2=b*b
      if(A.eq.197) then
         avnbcoll=(1+(-0.0364147+(0.000450308-1.8665e-6*b2)*b2)*b2)/
     &        (0.219715+(-0.00501456+(0.00002550504+1.116294e-6*b2)
     &        *b2)*b2)
      else if(A.eq.208) then
         if(b.lt.9.9) then
            avnbcoll=(1+(-0.0359115+(0.0004380525-1.79274d-6*b2)*b2)
     & *b2)/(0.215603+(-0.00493175+(0.0000263330+9.19488d-7*b2)
     & *b2)*b2)
         else
            avnbcoll=0
         endif
      else
         write(8,*) 'A not defined in get_ktkick_pA:',A
         close(8)
         stop
      endif
      if(avnbcoll.lt.0) avnbcoll=0
      kt2=sqrt(-avnbcoll*delta0*log(ran2()))
      theta=2*pi*ran2()
      kt(0)=0.D0
      kt(1)=kt2*cos(theta)
      kt(2)=kt2*sin(theta)
      kt(3)=0.D0
      return
      end


      subroutine get_ktkick_AA(xcm,A,b,kt)
      implicit none
      integer A
      double precision b,b1,b2,kt(0:3),kt1(0:3),kt2(0:3),xcm(0:3),ycm2
      ycm2=xcm(2)**2
      b1=sqrt((xcm(1)-b/2.d0)**2+ycm2)
      call get_ktkick_pA(A,b1,kt1)
      b2=sqrt((xcm(1)+b/2.d0)**2+ycm2)
      call get_ktkick_pA(A,b2,kt2)
      kt(0)=0.D0
      kt(1)=kt1(1)+kt2(1)
      kt(2)=kt1(2)+kt2(2)
      kt(3)=0.D0
      end

      subroutine def_part_init(tauf,avnc0ncoll,avnpsi0ncoll,avnb0ncoll,
     &   avnups0ncoll,A,idistspatial,ifcronin)
      use PBGgenvar
      use PBGsceplas
      use PBGpartinmedium
      use PBGimparam
      implicit none
      logical ifcronin,rejectnuclear,ifcoal
      integer A,idistspatial,poisson,nc0,nb0,npsi0ncoll,npsi0,newpart,
     & k,l,m,ifail,imed
      double precision avnc0ncoll,avnpsi0ncoll,avnb0ncoll,avnups0ncoll,
     & xcm(0:3),x1(0:3),x2(0:3),p1(0:3),p2(0:3),p(0:3),kt(0:3),tauf
     &,rgold,ran2,
     & timefly,x_fly(0:3),temp,beta(3),densener,betavac(3),
     & maxantishadowingPsi,ktsum,rapid,mt,pb_coal
      parameter(rgold=6.4d0)
      integer ifmtx
      call getMonitorFileIndex(ifmtx)
      betavac(1)=0.d0
      betavac(2)=0.d0
      kt(0)=0.d0
      kt(1)=0.d0
      kt(2)=0.d0
      kt(3)=0.d0

      if(isce<7) then

      nc0=poisson(avnc0ncoll)
      do k=1,nc0
         do l=1,2
            if(firstemptypart.gt.npartmax) then
               write(ifmtx,*) 'def_part_init: npartmax too small'
               close(ifmtx) 
               stop
            endif
            newpart=firstemptypart
            partinfo_ires(newpart)=l
            partinfo_mother(newpart)=k
            partinfo_mass(newpart)=resinfo_mass(l)
            if(l.eq.1) then 
               call azget2x(1,xcm,x1,x2,b,A,tauf,idistspatial)
               if(ifcronin) call get_ktkick_AA(xcm,A,dble(b),kt)
                    ktsum=ktsum+sqrt(kt(1)*kt(1)+kt(2)*kt(2))
               call get2p_q(1,p1,p2,
     &                      partinfo_mass(newpart))               
               call pack_it_vacuum(1,p1)
                    call pack_it_vacuum_cronin(1,p1+kt)
                    betavac(3)=p1(3)/p1(0)
                    call hadronizebasicnew(x1,p1,1,p,4,1,betavac,0,
     &                   pb_coal,ifcoal,.false.)
                    call pack_it_vacuum(4,p)
                    call hadronizebasicnew(x1,p1+kt,1,p,4,1,betavac,0,
     &                   pb_coal,ifcoal,.false.)
                    call pack_it_vacuum_cronin(4,p+kt)
                    call pack_it_vacuum(2,p2)
                    call pack_it_vacuum_cronin(2,p2+kt)
                    betavac(3)=p2(3)/p2(0)               
                    call hadronizebasicnew(x2,p2,2,p,5,1,betavac,0,
     &                   pb_coal,ifcoal,.false.)
                    call pack_it_vacuum(5,p)
                    call hadronizebasicnew(x2,p2+kt,2,p,5,1,betavac,0,
     &                   pb_coal,ifcoal,.false.)
                    call pack_it_vacuum_cronin(5,p+kt)
                rapid=0.5d0*log((p1(0)+p1(3))/(p1(0)-p1(3)))
                p1(1)=p1(1)+kt(1)
                p1(2)=p1(2)+kt(2)
                mt=sqrt(partinfo_mass(newpart)**2+p1(1)**2+p1(2)**2)
                p1(3)=mt*sinh(rapid)
                do m=1,3
                   partinfo_r(newpart,m)=x1(m)
                   partinfo_p(newpart,m)=p1(m)
                enddo
                partinfo_r(newpart,0)=x1(0)
                partinfo_tb(newpart)=sqrt(x1(0)**2-x1(3)**2)
             else
                rapid=0.5d0*log((p2(0)+p2(3))/(p2(0)-p2(3)))
                p2(1)=p2(1)+kt(1)
                p2(2)=p2(2)+kt(2)
                mt=sqrt(partinfo_mass(newpart)**2+p2(1)**2+p2(2)**2)
                p2(3)=mt*sinh(rapid)
                do m=1,3
                   partinfo_r(newpart,m)=x2(m)
                   partinfo_p(newpart,m)=p2(m)
                enddo
                partinfo_r(newpart,0)=x2(0)
                partinfo_tb(newpart)=sqrt(x2(0)**2-x2(3)**2)
             endif
            partinfo_p(newpart,0)=sqrt(partinfo_p(newpart,1)**2+
     &           partinfo_p(newpart,2)**2+partinfo_p(newpart,3)**2+
     &           partinfo_mass(newpart)**2)             
            partinfo_whichmedium(newpart)=-2
 

            firstemptypart=partinfo_next(firstemptypart)
            partinfo_next(newpart)=firstpart
            firstpart=newpart
            npart=npart+1
         enddo
      enddo
      resinfo_number(1)=nc0
      resinfo_number(2)=nc0
      nbccbaremission=nc0
      endif

      timefly=1.001*tau0
      if(avnpsi0ncoll.gt.0.d0) then
         npsi0ncoll=poisson(avnpsi0ncoll*maxantishadowingPsi())
      else
         npsi0ncoll=0
      endif
      npsi0=0
      do k=1,npsi0ncoll
         if(firstemptypart.gt.npartmax) then
            write(ifmtx,*) 'def_part_init: npartmax too small'
            close(ifmtx)
            stop
         endif
         newpart=firstemptypart
         if(idistspatial.eq.1) then
            call azget1x(x1,b,A,tauf)
         else
            call azget1xnpart(x1,b,A,tauf)
         endif
         if(ifcronin) call get_ktkick_AA(x1,A,dble(b),kt)
         call getonep_psi(p,resinfo_mass(3))
         if(ran2().lt.(1./maxantishadowingPsi())) then
            call pack_it_vacuum(3,p) 
         endif  
         if(rejectnuclear(x1,b,p)) goto 10
          rapid=0.5d0*log((p(0)+p(3))/(p(0)-p(3)))
          p(1)=p(1)+kt(1)
          p(2)=p(2)+kt(2)
          mt=sqrt(resinfo_mass(3)**2+p(1)**2+p(2)**2)
          p(3)=mt*sinh(rapid)
         p(0)=sqrt(resinfo_mass(3)**2+p(1)**2+p(2)**2+p(3)**2)
         call freeflyone(x1,p,mpsi,timefly,itypeorder,x_fly)
         call givelocalhydro(x_fly,temp,beta,densener,ifail,imed)
         if(imed.le.-1) then
            write(ifmtx,*) 'def_part_init: imed<0 after free fly ?!'
            close(ifmtx)
            stop
         endif
         if((imed.gt.0).or.(temp.lt.tempdissoc(3))) then
            partinfo_ires(newpart)=3
            partinfo_mother(newpart)=0
            partinfo_mass(newpart)=resinfo_mass(3)
            do m=0,3
               partinfo_r(newpart,m)=x1(m)
               partinfo_p(newpart,m)=p(m)
            enddo
            partinfo_tb(newpart)=sqrt(x1(0)**2-x1(3)**2)
            partinfo_whichmedium(newpart)=imed
            firstemptypart=partinfo_next(firstemptypart)
            partinfo_next(newpart)=firstpart
            firstpart=newpart
            npart=npart+1
            npsi0=npsi0+1
         else
         endif
 10      continue
      enddo
      resinfo_number(3)=npsi0
      resinfo_number(4)=0
      resinfo_number(5)=0
      if(isce<7) then

      nb0=poisson(avnb0ncoll)
      do k=1,nb0
         do l=1,2
            if(firstemptypart.gt.npartmax) then
               write(ifmtx,*) 'def_part_init: npartmax too small'
               close(ifmtx)
               stop
            endif
            newpart=firstemptypart
            partinfo_ires(newpart)=nresc+l
            partinfo_mother(newpart)=k
            partinfo_mass(newpart)=resinfo_mass(6)
            if(l.eq.1) then 
               call azget2x(2,xcm,x1,x2,b,A,tauf,idistspatial)
               if(ifcronin) call get_ktkick_AA(xcm,A,dble(b),kt)
               call get2p_q(2,p1,p2,partinfo_mass(newpart))
               call pack_it_vacuum(6,p1)
               betavac(3)=p1(3)/p1(0)
               call hadronizebasicnew(x1,p1,6,p,9,1,betavac,0,
     &                   pb_coal,ifcoal,.false.)
               call pack_it_vacuum(9,p)
               call pack_it_vacuum(7,p2)
               betavac(3)=p2(3)/p2(0)               
               call hadronizebasicnew(x2,p2,7,p,10,1,betavac,0,
     &                   pb_coal,ifcoal,.false.)
               call pack_it_vacuum(10,p)
                rapid=0.5d0*log((p1(0)+p1(3))/(p1(0)-p1(3)))
                p1(1)=p1(1)+kt(1)
                p1(2)=p1(2)+kt(2)
                mt=sqrt(partinfo_mass(newpart)**2+p1(1)**2+p1(2)**2)
                p1(3)=mt*sinh(rapid)
                do m=1,3
                   partinfo_r(newpart,m)=x1(m)
                   partinfo_p(newpart,m)=p1(m)
                enddo
                partinfo_r(newpart,0)=x1(0)
                partinfo_tb(newpart)=sqrt(x1(0)**2-x1(3)**2)
             else
                rapid=0.5d0*log((p2(0)+p2(3))/(p2(0)-p2(3)))
                p2(1)=p2(1)+kt(1)
                p2(2)=p2(2)+kt(2)
                mt=sqrt(partinfo_mass(newpart)**2+p2(1)**2+p2(2)**2)
                p2(3)=mt*sinh(rapid)
                do m=1,3
                   partinfo_r(newpart,m)=x2(m)
                   partinfo_p(newpart,m)=p2(m)
                enddo
                partinfo_r(newpart,0)=x2(0)
                partinfo_tb(newpart)=sqrt(x2(0)**2-x2(3)**2)
             endif

            partinfo_p(newpart,0)=sqrt(partinfo_p(newpart,1)**2+
     &           partinfo_p(newpart,2)**2+partinfo_p(newpart,3)**2+
     &           partinfo_mass(newpart)**2)
            partinfo_whichmedium(newpart)=-2
            firstemptypart=partinfo_next(firstemptypart)
            partinfo_next(newpart)=firstpart
            firstpart=newpart
            npart=npart+1
         enddo
      enddo
      resinfo_number(6)=nb0
      resinfo_number(10)=nb0
      nbbbbaremission=nb0

      endif

      if(avnups0ncoll.ne.0) then
         write(ifmtx,*) 'no upsilon included up to now in def_part'
         close(ifmtx)
         stop 
      endif
      if(infotrack.ge.1) then
       write(6,*) 'Initial Particles Created, npart=',npart
      endif
      return
      end




      subroutine def_part_init_epos_512(ifcronin)
      use PBGiniphasespace
      use PBGinixypos
      use PBGsystem
      use PBGblckinitptspect3
      use PBGsigmahqprod
      use PBGvariancesig
      use PBGgenvar
      use PBGsceplas
      use PBGpartinmedium
      use PBGhydro
      use PBGspectra
      use PBGdistccbarbis
      use PBGinitnbmult
      implicit none
      logical ifcronin,ifcoal
      integer newpart,k,l,m,i,numberprodpoints,ires1,ires2
      double precision xcm(0:3),x1(0:3),x2(0:3),p1(0:3),p2(0:3),p(0:3),
     & kt(0:3),ran2,betavac(3),rapid,mt,pb_coal,rrel,thetarel,xrel,
     & yrel,ii4,rtotal,rbbbar,rccbar,pi

      parameter (pi=3.14159265358979323844d0)
      integer sumcs,sumbs
      save sumcs,sumbs
      integer ifmtx
      call getMonitorFileIndex(ifmtx)

      if(abs(sqrtS-200)<0.001) then
         rccbar=0.013904d0
         rbbbar=0.00007d0
      elseif(abs(sqrtS-5250)<260) then
         rccbar=0.1591
         rbbbar=0.004534d0
      elseif(abs(sqrtS-2755)<10.) then
         rccbar=0.0819d0*4/3.d0*0.66
         rbbbar=0.0027d0
      else
         write(ifmtx,*) 'HQ: unforeseen sqrtS:',sqrtS
         close(ifmtx)
         stop
      endif
      rccbar=rccbar*avnc0ncollmult      
      rbbbar=rbbbar*avnb0ncollmult      
      rtotal=rccbar+rbbbar
      betavac(1)=0.d0
      betavac(2)=0.d0
      kt(0)=0.d0
      kt(1)=0.d0
      kt(2)=0.d0
      kt(3)=0.d0
      numberprodpoints=numberofNNcolls
      sumcs=0
      sumbs=0
      do i=1,numberofHQ/2
         k=int(ran2()*numberprodpoints)+1
         ii4=ran2()*rtotal
         if(ii4<rbbbar) then 
            if(.not.(ifshadowing(2)).or.(ifshadowing(2).and.
     &           (ran2().lt.globshadow(2)))) then
               sumbs=sumbs+1
               do l=1,2
                  if(firstemptypart.gt.npartmax) then
                     write(ifmtx,*) 'error in def_part_init_epos:' 
                     write(ifmtx,*) 'npartmax too small'
                     stop
                  endif
                  newpart=firstemptypart
                  partinfo_ires(newpart)=nresc+l
                  partinfo_mother(newpart)=k
                  partinfo_mass(newpart)=resinfo_mass(nresc+l)
                  if(l.eq.1) then 
                     ires1=nresc+1
                     ires2=nresc+2
                     xcm(0)=0.d0
                     xcm(3)=0.d0
                     xcm(1)=xyarray(k,1)
                     xcm(2)=xyarray(k,2)
                     x1(0)=xcm(0) 
                     x2(0)=xcm(0)
                     rrel=0.05d0*sqrt(ran2())
                     thetarel=2.0d0*pi*ran2()
                     xrel=rrel*dcos(thetarel)
                     yrel=rrel*dsin(thetarel)
                     x1(1)=xcm(1)+xrel
                     x2(1)=xcm(1)-xrel
                     x1(2)=xcm(2)+yrel
                     x2(2)=xcm(2)-yrel
                     x1(3)=xcm(3)
                     x2(3)=xcm(3)
                     if(ifcronin) then
                        call get_ktkick_AA(xcm,A,dble(b),kt)
                     end if
                     call get2p_q(2,p1,p2,partinfo_mass(newpart))
                     call pack_it_vacuum(ires1,p1)
                     betavac(3)=p1(3)/p1(0)
                           call hadronizebasicnew(x1,p1,ires1,p,
     &                          ires1+3,1,betavac,0,pb_coal,ifcoal,
     &                          .false.)
                     call pack_it_vacuum(ires1+3,p)
                     call pack_it_vacuum(ires2,p2)
                     betavac(3)=p2(3)/p2(0)               
                           call hadronizebasicnew(x2,p2,ires2,p,
     &                          ires2+3,1,betavac,0,
     &                          pb_coal,ifcoal,.false.)
                     call pack_it_vacuum(ires2+3,p)
                     if(ifcronin) then
                        mt=sqrt(partinfo_mass(newpart)**2+p1(1)**2
     &                       +p1(2)**2)
                        rapid=sign(0.5d0*log((p1(0)+abs(p1(3)))/
     &                       mt),p1(3))
                        p1(1)=p1(1)+kt(1)
                        p1(2)=p1(2)+kt(2)
                        mt=sqrt(partinfo_mass(newpart)**2+p1(1)**2
     &                       +p1(2)**2)
                        p1(3)=mt*sinh(rapid)
                     endif
                     do m=1,3
                        partinfo_r(newpart,m)=x1(m)
                        partinfo_p(newpart,m)=p1(m)
                     enddo
                     partinfo_r(newpart,0)=x1(0)
                     partinfo_tb(newpart)=sqrt(x1(0)**2-x1(3)**2)
                  else
                     if(ifcronin) then
                        mt=sqrt(partinfo_mass(newpart)**2+p2(1)**2
     &                       +p2(2)**2)
                        rapid=sign(0.5d0*log((p2(0)+abs(p2(3)))/ 
     &                       mt),p2(3))
                        p2(1)=p2(1)+kt(1)
                        p2(2)=p2(2)+kt(2)
                        mt=sqrt(partinfo_mass(newpart)**2+p2(1)**2
     &                       +p2(2)**2)
                        p2(3)=mt*sinh(rapid)
                     endif
                     do m=1,3
                        partinfo_r(newpart,m)=x2(m)
                        partinfo_p(newpart,m)=p2(m)
                     enddo
                     partinfo_r(newpart,0)=x2(0)
                     partinfo_tb(newpart)=sqrt(x2(0)**2-x2(3)**2)
                  endif                
                  partinfo_p(newpart,0)=
     &                 sqrt(partinfo_p(newpart,1)**2+
     &                 partinfo_p(newpart,2)**2+
     &                 partinfo_p(newpart,3)**2+
     &                 partinfo_mass(newpart)**2)
                  partinfo_whichmedium(newpart)=-2
                  partinfo_corona(newpart)=.true.
                  firstemptypart=partinfo_next(firstemptypart)
                  partinfo_next(newpart)=firstpart
                  firstpart=newpart
                  npart=npart+1
               enddo
            endif
         else
            if(.not.(ifshadowing(1)).or.(ifshadowing(1).and.
     &           (ran2().lt.globshadow(1)))) then
               sumcs=sumcs+1
               do l=1,2
                  if(firstemptypart.gt.npartmax) then
                     write(ifmtx,*) 'error in def_part_init:' 
                     write(ifmtx,*) 'npartmax too small'
                     stop
                  endif
                  newpart=firstemptypart
                  partinfo_ires(newpart)=l
                  partinfo_mother(newpart)=k
                  partinfo_mass(newpart)=resinfo_mass(l)
                  if(l.eq.1) then 
                     ires1=1
                     ires2=2
                     xcm(0)=0.d0
                     xcm(3)=0.d0
                     xcm(1)=xyarray(k,1)
                     xcm(2)=xyarray(k,2)
                     x1(0)=xcm(0) 
                     x2(0)=xcm(0)
                     if(typedistrel.eq.0) then
                        rrel=paramdistrel(1)/2*sqrt(ran2())
                        thetarel=2.0d0*pi*ran2()
                        xrel=rrel*dcos(thetarel)
                        yrel=rrel*dsin(thetarel)
                     elseif(typedistrel.eq.1) then
                        rrel=paramdistrel(1)/2
                        thetarel=2.0d0*pi*ran2()
                        xrel=rrel*dcos(thetarel)
                        yrel=rrel*dsin(thetarel)
                     else
                        call gauss(xrel,yrel)
                        xrel=xrel*paramdistrel(1)/sqrt(2.d0)/2.d0
                        yrel=yrel*paramdistrel(1)/sqrt(2.d0)/2.d0            
                     endif
                     x1(1)=xcm(1)+xrel
                     x2(1)=xcm(1)-xrel
                     x1(2)=xcm(2)+yrel
                     x2(2)=xcm(2)-yrel
                     x1(3)=xcm(3)
                     x2(3)=xcm(3)
                     if(ifcronin) then
                        call get_ktkick_AA(xcm,A,dble(b),kt)
                     end if
                     call get2p_q(1,p1,p2,partinfo_mass(newpart))
                     call pack_it_vacuum(ires1,p1)
                     betavac(3)=p1(3)/p1(0)
                     call hadronizebasicnew(x1,p1,ires1,p,
     &                          ires1+3,1,betavac,0,pb_coal,ifcoal,
     &                          .false.)
                     call pack_it_vacuum(ires1+3,p)
                     call pack_it_vacuum(ires2,p2)
                     betavac(3)=p2(3)/p2(0)               
                     call hadronizebasicnew(x2,p2,ires2,p,
     &                          ires2+3,1,betavac,0,
     &                          pb_coal,ifcoal,.false.)
                     call pack_it_vacuum(ires2+3,p)
                     if(ifcronin) then
                        mt=sqrt(partinfo_mass(newpart)**2+p1(1)**2
     &                       +p1(2)**2)
                        rapid=sign(0.5d0*log((p1(0)+abs(p1(3)))/
     &                       mt),p1(3))
                        p1(1)=p1(1)+kt(1)
                        p1(2)=p1(2)+kt(2)
                        mt=sqrt(partinfo_mass(newpart)**2+p1(1)**2
     &                       +p1(2)**2)
                        p1(3)=mt*sinh(rapid)
                     endif
                     do m=1,3
                        partinfo_r(newpart,m)=x1(m)
                        partinfo_p(newpart,m)=p1(m)
                     enddo
                     partinfo_r(newpart,0)=x1(0)
                     partinfo_tb(newpart)=sqrt(x1(0)**2-x1(3)**2)
                  else
                     if(ifcronin) then
                        mt=sqrt(partinfo_mass(newpart)**2+p2(1)**2
     &                       +p2(2)**2)
                        rapid=sign(0.5d0*log((p2(0)+abs(p2(3)))/
     &                       mt),p2(3))
                        p2(1)=p2(1)+kt(1)
                        p2(2)=p2(2)+kt(2)
                        mt=sqrt(partinfo_mass(newpart)**2+p2(1)**2
     &                       +p2(2)**2)
                        p2(3)=mt*sinh(rapid)
                     endif
                     do m=1,3
                        partinfo_r(newpart,m)=x2(m)
                        partinfo_p(newpart,m)=p2(m)
                     enddo
                     partinfo_r(newpart,0)=x2(0)
                     partinfo_tb(newpart)=sqrt(x2(0)**2-x2(3)**2)
                  endif                
                  partinfo_p(newpart,0)=
     &                 sqrt(partinfo_p(newpart,1)**2+
     &                 partinfo_p(newpart,2)**2+
     &                 partinfo_p(newpart,3)**2+
     &                 partinfo_mass(newpart)**2)
                  partinfo_whichmedium(newpart)=-2
                  firstemptypart=partinfo_next(firstemptypart)
                  partinfo_next(newpart)=firstpart
                  firstpart=newpart
                  npart=npart+1
               enddo
            endif
         endif
      enddo
      resinfo_number(1)=sumcs
      resinfo_number(2)=sumcs
      nbccbaremission=sumcs
      resinfo_number(6)=sumbs
      resinfo_number(10)=sumbs
      nbbbbaremission=sumbs
      return
      end subroutine def_part_init_epos_512

      subroutine def_part_init_epos3
      use PBGiniphasespace
      use PBGgenvar
      use PBGsceplas
      use PBGpartinmedium
      use PBGhydro
      use PBGfordisplay
      use PBGdistccbarbis
      implicit none
      logical ifunmatch,ifcoal
      integer k,kmother,l,newpart,m,id,ires,idprev,iresnewini
      double precision x1(0:3),p1(0:3),betavac(3),p(0:3)
      double precision ran2,rrel,thetarel,xrel,yrel,rapid,mt,pb_coal,pi
      integer ifmtx

      parameter (pi=3.14159265358979323844d0)
      call getMonitorFileIndex(ifmtx)
      betavac(1)=0.d0
      betavac(2)=0.d0            
      ifunmatch=.false.      
      do kmother=1,numberofHQ/2
         do l=1,2
            if(firstemptypart.gt.npartmax) then
               write(ifmtx,*) 'in def_part_init: npartmax too small'
               close(ifmtx)
               stop
            endif
            newpart=firstemptypart
            k=2*(kmother-1)+l
            id=HQarray(k,1)
            ires=0
            if(id.eq.4)ires=1
            if(id.eq.-4)ires=2
            if(id.eq.5)ires=14
            if(id.eq.-5)ires=15
            if(l.eq.1) then
               idprev=id
               if((ires.eq.14).or.(ires.eq.15)) then
                  rrel=0.05d0*sqrt(ran2())
                  thetarel=2.0d0*pi*ran2()
                  xrel=rrel*dcos(thetarel)
                  yrel=rrel*dsin(thetarel)
               elseif((ires.eq.1).or.(ires.eq.2)) then
                  if(typedistrel.eq.0) then
                     rrel=paramdistrel(2)/2*sqrt(ran2())
                     thetarel=2.0d0*pi*ran2()
                     xrel=rrel*dcos(thetarel)
                     yrel=rrel*dsin(thetarel)
                  elseif(typedistrel.eq.1) then
                     rrel=paramdistrel(2)/2
                     thetarel=2.0d0*pi*ran2()
                     xrel=rrel*dcos(thetarel)
                     yrel=rrel*dsin(thetarel)
                  else
                     call gauss(xrel,yrel)
                     xrel=xrel*paramdistrel(2)/sqrt(2.d0)/2.d0
                     yrel=yrel*paramdistrel(2)/sqrt(2.d0)/2.d0            
                  endif
               else
                  write(ifmtx,*) 'unknown id found in def_part_init:',id
                  close(ifmtx)
                  stop
               endif
            else
               if(idprev+id.ne.0) then
                  ifunmatch=.true.
               endif
            endif
            partinfo_ires(newpart)=ires
            partinfo_mother(newpart)=HQmother(k)
            if(l.eq.2) then
               if(HQmother(k).ne.HQmother(k-1)) then
                  write(ifmtx,*) '!! mothers:',HQmother(k),HQmother(k-1)
                  write(*,*) '!! mothers:',HQmother(k),HQmother(k-1)
                  close(ifmtx)
                  stop
               endif
            endif
            partinfo_mass(newpart)=resinfo_mass(ires)
            x1(1)=HQarray(k,7)-(-1)**l*xrel
            x1(2)=HQarray(k,8)-(-1)**l*yrel
            x1(3)=HQarray(k,9)
            partinfo_tb(newpart)=0.001
            x1(0)=sqrt(x1(3)**2+partinfo_tb(newpart)**2)
            do m=0,3
               partinfo_r(newpart,m)=x1(m)
            enddo
            p1(1)=HQarray(k,2)
            p1(2)=HQarray(k,3)
            p1(3)=HQarray(k,4)
            p1(0)=HQarray(k,5)
            mt=sqrt(p1(1)**2+p1(2)**2+HQarray(k,6)**2)
            rapid=dsign(1.d0,p1(3))*log((p1(0)+abs(p1(3)))/mt)
            mt=sqrt(p1(1)**2+p1(2)**2+partinfo_mass(newpart)**2)
            p1(3)=mt*sinh(rapid)
            p1(0)=mt*cosh(rapid)
            do m=0,3
               partinfo_p(newpart,m)=p1(m)
            enddo
            call pack_it_vacuum(ires,p1)
            betavac(3)=p1(3)/p1(0)

            call hadronizebasicnewall(x1,p1,ires,p,iresnewini,1,betavac,
     &      0,pb_coal,ifcoal,.false.)

            call pack_it_vacuum(iresnewini,p)
            partinfo_whichmedium(newpart)=-2
            partinfo_corona(newpart)=.true.
            firstemptypart=partinfo_next(firstemptypart)
            partinfo_next(newpart)=firstpart
            firstpart=newpart
            npart=npart+1
         enddo
      enddo
      resinfo_number(1)=numberofHQ/2
      resinfo_number(2)=numberofHQ/2
      nbccbaremission=numberofHQ/2  

C     &               partinfo_r(k,1),partinfo_r(k,2),
C     &                  partinfo_r(k,3) 
      return
      end
     

      subroutine freeflyone(x_in,p,mass,time,itorder,x_fin)
      implicit none
      integer k,itorder
      double precision x_in(0:3),x_fin(0:3),p(0:3),time,timefin,mass,
     &  gtimefin,mt2,dt
      integer ifmtx
      call getMonitorFileIndex(ifmtx)
         if (itorder.eq.1) then
            if(time.lt.x_in(0)) then
               write(ifmtx,*) 'error in freeflyone: part. time ',
     &            x_in(0),' > ',time
               close(ifmtx)
               stop
            else
               timefin=time
            endif
         else
            mt2=mass**2+p(1)**2+p(2)**2
            timefin=gtimefin(x_in(0),x_in(3),p(0),p(3),mt2,time)
         endif
         dt=(timefin-x_in(0))/p(0)
         do k=1,3
            x_fin(k)=x_in(k)+p(k)*dt
         enddo
         x_fin(0)=timefin
         return
         end

      subroutine freefly(time)
      use PBGgenvar
      use PBGsceplas
      use PBGpartinmedium
      use PBGevolQ
      implicit none
      integer ipart,k,ifail,imed
      double precision time,timefin,gtimefin2,mass,dt,x(0:3),temp,
     & beta(3),densener
      integer ifmtx
      call getMonitorFileIndex(ifmtx)
      ipart=firstpart
      do while (ipart.NE.-1)
         if (itypeorder.eq.1) then
            if(time.lt.partinfo_r(ipart,0)) then
               write(ifmtx,*) 'error in freefly: part. time ',
     &              partinfo_r(ipart,0),' > ',time
               close(ifmtx)
               stop
            else
               timefin=time
            endif
         else
            mass=resinfo_mass(partinfo_ires(ipart))
            timefin=gtimefin2(ipart,time)
         endif
         dt=(timefin-partinfo_r(ipart,0))/partinfo_p(ipart,0)
         
         do k=1,3
            partinfo_r(ipart,k)=partinfo_r(ipart,k)
     &                         +partinfo_p(ipart,k)*dt
         enddo
         partinfo_r(ipart,0)=timefin
         if(itypeorder.eq.1) then            
            partinfo_tb(ipart)=sqrt(partinfo_r(ipart,0)**2-
     &                              partinfo_r(ipart,3)**2)
         else
            partinfo_tb(ipart)=time
         endif   
         do k=0,3
            x(k)=partinfo_r(ipart,k)
         enddo
         call givelocalhydro(x,temp,beta,densener,ifail,imed)
         if(imed.le.-1) then
            write(ifmtx,*) 'imed<0 after free fly ?!'
            write(ifmtx,*) x(0),x(1),x(2),x(3),ifail
            close(ifmtx)
            stop
         endif
         partinfo_whichmedium(ipart)=imed 
         if(partinfo_whichmedium(ipart).gt.imedfinal) then
            if(resinfo_ifquark(partinfo_ires(ipart))) then
               call hadronize(ipart,1,beta)
            endif
         else
            if(temp.gt.tempdissoc(partinfo_ires(ipart))) then
               write(ifmtx,*) 'J/Psi found above diss Temp in free fly'
               write(ifmtx,*) partinfo_ires(ipart)
               close(ifmtx)
               stop
            endif   
         endif   
         ipart=partinfo_next(ipart)
      enddo
      return
      end


      subroutine update_part_chain_old(tord0)
      use JA
      use PBGblckhadronize
      use PBGgenvar
      use PBGsceplas
      use PBGpartinmedium
      use PBGfordisplay
      use PBGevolQ
      use PBGreductiondofs
      implicit none 
      integer partcount,ipart,k,ifail_med,ifaillang,imed,
     &  ifail_medbid,imedbid,ires
      double precision tord0,time,x_part(0:3),x_part_old(0:3),
     & p_part(0:3),deltap(0:3),gtimefin2,dt,dtau,temp,beta(3),densener,
     & effdeg
      double precision lambdaepos,reduction_dof
      external lambdaepos,reduction_dof
      ipart=firstpart
      partcount=0
      if(itypeorder.eq.1) time=tord0
      do while (ipart.NE.-1)
         jaipart=ipart
         partcount=partcount+1         
         ires=partinfo_ires(ipart)   
         if(itypeorder.eq.2) time=gtimefin2(ipart,tord0)
         do k=0,3
            p_part(k)=partinfo_p(ipart,k)
         enddo   
         dt=time-partinfo_r(ipart,0)
         dtau=dt/p_part(0)
         do k=1,3
            x_part_old(k)=partinfo_r(ipart,k)
            x_part(k)=x_part_old(k)+p_part(k)*dtau
         enddo
         x_part_old(0)=partinfo_r(ipart,0)
         x_part(0)=time
         call givelocalhydro(x_part,temp,beta,densener,ifail_med,imed)
         if(imed.lt.0) then
            write(8,*)'in update_part:part should be at least in QGP'
            write(120,*)'in update_part:part should be at least in QGP'
            close(8)
         endif   
         if((imed.ge.0).and.(imed.le.imedfinal)) then
            if(partinfo_whichmedium(ipart).gt.imedfinal) then
               if(resinfo_ifquark(ires)) then
                write(8,*) 'stupidity in update: no quark can reenter'
                  write(120,*) 'stupidity: no quark can reenter'
                  close(8)
                  stop
               endif
            endif
            if(.not.resinfo_ifquark(ires)) then
               if(temp.gt.tempdissoc(ires)) call melt(ipart)     
            else



                if(reddof<3) then
                effdeg=reduction_dof(1,10.0d0,temp)
                else
                effdeg=1.0
                endif

                if(isce<7) then
                if((imed.eq.0).or.(.not.(ifdecrease))) then
                   effdeg=1.d0
                else
                   effdeg=densener/entrans_max
                endif
                endif
               if(type_evol_Q.eq.1) then
                  call kickme_lang(p_part,partinfo_mass(ipart),dt,
     &                 resinfo_ifquark(ires),ires,temp,effdeg,beta,
     &                 ifaillang)  
                else

                  call kickme_bolt(p_part,dt,
     &                 resinfo_ifquark(ires),ires,temp,effdeg,beta,
     &                 deltap,ifaillang,.false.)
               endif   
               do k=0,3
                  partinfo_p(ipart,k)=p_part(k)
               enddo
            endif                
         endif  
	    if(ipart.le.5) then
         endif
         if((imed.gt.imedfinal).and.
     &       (partinfo_whichmedium(ipart).le.imedfinal)) then
            if(resinfo_ifquark(ires)) then
               call givelocalhydro(x_part_old,temp,beta,densener,
     &              ifail_medbid,imedbid)
               call hadronize(ipart,itypDproduct,beta)
            endif
         endif  
         partinfo_whichmedium(ipart)=imed
         do k=0,3
            partinfo_r(ipart,k)=x_part(k)
         enddo
         partinfo_tb(ipart)=sqrt(x_part(0)**2-x_part(3)**2)
         ipart=partinfo_next(ipart)
      enddo
      return
      end

      subroutine update_part_chain(tord0,nbclose,ifail)
      use JA
      use PBGblckhadronize
      use PBGgenvar
      use PBGsceplas
      use PBGpartinmedium
      use PBGfordisplay
      use PBGevolQ
      use PBGreductiondofs
      use denysvar
      use for2body
      use onerateHQparam
      use genuineRemler
      implicit none 
      logical iftimecheck,ifdisplay,accept
      integer partcount,ipart,k,ifail,ifail_med,ifaillang,imed,
     &  ifail_medbid,imedbid,ires,ipart2,partcount2,nbclose,
     &  ifailup,nbwarn,nbhqbis
      save nbwarn
      data nbwarn/0/
      double precision tord0,time,x_part(0:3),x_part_old(0:3),
     &  p_part(0:3),p_part_old(0:3),
     &  gtimefin2,dt,dtau,temp,beta(3),densener,
     &  effdeg,x_part2(0:3),
     &  mass1,mass2,p_part2(0:3),dx_part(0:3),dp_part(0:3),rrelinv,
     &  xMRU(0:3),tbinst,entotcm,x01fin(0:3),p01fin(0:3),deltp(0:3),
     &  x02fin(0:3),p02fin(0:3),dx_part2(0:3),dp_part2(0:3),dp_cm(3),
     &  coherefact,ran2
      integer nbcontribgt2
      logical, dimension(:), allocatable :: ifquarkQGP
      integer, dimension(:), allocatable :: imedsave,nbcontrib
      real, dimension(:), allocatable :: tempsave
      double precision, dimension(:,:), allocatable :: drstore,dpstore
      integer, dimension(:), allocatable :: minrrelpartner
      double precision, dimension(:,:), allocatable :: minrrel
      double precision, dimension(:), allocatable :: minrrelinv
      integer npart2
        double precision  absdeltap
      double precision lambdaepos,reduction_dof
      external lambdaepos,reduction_dof
      logical iffirst,ifinfo
      data iffirst/.true./
      save iffirst
      double precision wignerrelat,pTpsi,phipsi,rapidpsi,rjpsi,
     & pjpsi,gamma,ergcm
      integer ifmtx,counter
      call getMonitorFileIndex(ifmtx)
      ifinfo=.false.
      ifail=0
      counter=0
      nbclose=0
      npart2=0
      ipart=firstpart
      do while (ipart.NE.-1)
         npart2=npart2+1
         ipart=partinfo_next(ipart)
      enddo
      allocate(ifquarkQGP(npart2))
      allocate(imedsave(npart2))
      allocate(nbcontrib(npart2))
      nbcontribgt2=0
      allocate(tempsave(npart2))
      allocate(drstore(npart2,0:3))
      allocate(dpstore(npart2,3))
      allocate(minrrel(npart2,0:3))  
      allocate(minrrelpartner(npart2))
      allocate(minrrelinv(npart2))
      if(itypeorder.eq.1) then 
         time=tord0 
         iftimecheck=.false.
      elseif(itypeorder.eq.2) then
         iftimecheck=.true.
      else
         write(*,*) 'unknown itypeorder:',itypeorder
         write(ifmtx,*) 'unknown itypeorder:',itypeorder
         close(ifmtx)
         stop
      endif
      ipart=firstpart
      partcount=0
      do while (ipart.NE.-1)
         partcount=partcount+1         
         ires=partinfo_ires(ipart)   
         if(itypeorder.eq.2) time=gtimefin2(ipart,tord0)
         do k=0,3
            p_part(k)=partinfo_p(ipart,k)
         enddo   
         dt=time-partinfo_r(ipart,0)
         dtau=dt/p_part(0)
         do k=1,3
            x_part(k)=partinfo_r(ipart,k)+p_part(k)*dtau
         enddo
         x_part(0)=time
         call givelocalhydro(x_part,temp,beta,densener,ifail_med,
     &        imed)
         if(imed.lt.0) then
            write(ifmtx,*) 'update_par:pcle should be @ least in QGP'
            close(ifmtx)
         endif   
         tempsave(partcount)=temp
         imedsave(partcount)=imed
         ifquarkQGP(partcount)=.false.
         nbcontrib(partcount)=0
         if((imed.ge.0).and.(imed.le.imedfinal)) then
            if(resinfo_ifquark(ires)) then
               ifquarkQGP(partcount)=.true.
            endif
         endif  
         drstore(partcount,0)=0.d0
         minrrel(ipart,0)=0.d0
         minrrelpartner(ipart)=0
         do k=1,3
            drstore(partcount,k)=0.d0
            dpstore(partcount,k)=0.d0
            minrrel(ipart,k)=0.d0
         enddo
         minrrelinv(partcount)=100.d0
         ipart=partinfo_next(ipart)
      enddo
      ipart=firstpart
      partcount=0
      if(itypeorder.eq.1) time=tord0
      do while (ipart.NE.-1)
         partcount=partcount+1         
         if(ifquarkQGP(partcount)) then
            if(itypeorder.eq.2) time=gtimefin2(ipart,tord0)
            dt=time-partinfo_r(ipart,0)
            do k=0,3
               x_part(k)=partinfo_r(ipart,k)
               p_part(k)=partinfo_p(ipart,k)
            enddo  
            mass1=partinfo_mass(ipart)
            ipart2=ipart
            partcount2=partcount-1
            do while (ipart2.NE.-1)
               partcount2=partcount2+1
               if(ifquarkQGP(partcount2)) then
                  if(.not.((imedsave(partcount2).ge.0).and.
     &                 (imedsave(partcount2).le.imedfinal))) then
                     write(ifmtx,*) 'inconsistency evolve'
                     close(ifmtx)
                     stop
                  endif
                  mass2=partinfo_mass(ipart2)
                  if((ifsingletfull(partcount,partcount2)).and.
     &                 (abs(mass2-mass1).lt.0.01))  then
                     do k=0,3
                        p_part2(k)=partinfo_p(ipart2,k)
                        x_part2(k)=partinfo_r(ipart2,k)
                     enddo
                     temp=(tempsave(partcount)+tempsave(partcount2))/
     &                    2.d0
                     ifdisplay=.false.
                     call oneevolQQbarbis(mass1,x_part,x_part2,
     &                    p_part,
     &                    p_part2,if2bodySC,ifdisplay,
     &                    x01fin,p01fin,
     &                    entotcm,rrelinv,dx_part,dp_part,xMRU,ifailup)
                     if(rrelinv.lt.minrrelinv(partcount)) then
                        minrrelinv(ipart)=rrelinv
                        minrrelpartner(ipart)=ipart2
                     endif
                     if((ifailup.gt.0).and.(ifailup.le.6)) then 
                        if((ifailup.eq.1).or.(ifailup.eq.2)) then
                           write(ifmtx,*) 'error ',ifailup,'onevolQQ' 
                           write(ifmtx,*) 'part:',ipart
                        endif
                        if((ifailup.eq.3).or.(ifailup.eq.4)) then
                           write(ifmtx,*) 'error ',ifailup,'onevolQQ' 
                           write(ifmtx,*) 'part:',ipart2
                        endif
                        if((ifailup.eq.5).or.(ifailup.eq.6)) then
                           write(ifmtx,*) 'error ',ifailup,'onevolQQ' 
                           write(ifmtx,*) 'part1:',ipart
                           write(ifmtx,*) 'part2:',ipart2
                        endif
                        if(ifallowhqskip) then
                           ifail=ifailup
                           return
                        else
                           close(ifmtx)
                           stop
                        endif
                     endif
                     if(rrelinv.le.1.d0) then
                        call gimmewignerrelatsinglet(p_part,p_part2,
     &                       x_part,x_part2,0.25d0,wignerrelat,pTpsi,
     &                       phipsi,rapidpsi,rjpsi,pjpsi,gamma,ergcm) 
                        if(abs(rapidpsi).lt.1.d0) nbclose=nbclose+1
                     endif
                     if(sqrt(dx_part(1)**2+dx_part(2)**2+
     &                    dx_part(3)**2+dx_part(0)**2).gt.1D-3) then
                        nbcontrib(partcount)=nbcontrib(partcount)+1
                     endif
                     call oneevolQQbarbis(mass2,x_part2,x_part,
     &                    p_part2,
     &                    p_part,if2bodySC,ifdisplay,
     &                    x02fin,p02fin,entotcm,rrelinv,dx_part2,
     &                    dp_part2,xMRU,ifailup)
                     if((ifailup.gt.0).and.(ifailup.le.6)) then 
                        if((ifailup.eq.1).or.(ifailup.eq.2)) then
                           write(ifmtx,*) 'error ',ifailup,'onevolQQ' 
                           write(ifmtx,*) 'part:',ipart
                        endif
                        if((ifailup.eq.3).or.(ifailup.eq.4)) then
                           write(ifmtx,*) 'error ',ifailup,'onevolQQ' 
                           write(ifmtx,*) 'part:',ipart2
                        endif
                        if((ifailup.eq.5).or.(ifailup.eq.6)) then
                           write(ifmtx,*) 'error ',ifailup,'onevolQQ' 
                           write(ifmtx,*) 'part1:',ipart
                           write(ifmtx,*) 'part2:',ipart2
                        endif
                        if(ifallowhqskip) then
                           ifail=ifailup
                           return
                        else
                           close(ifmtx)
                           stop
                        endif
                     endif
                     if(rrelinv.le.1.d0) then
                        call gimmewignerrelatsinglet(p_part2,p_part,
     &                       x_part2,x_part,0.25d0,wignerrelat,pTpsi,
     &                       phipsi,rapidpsi,rjpsi,pjpsi,gamma,ergcm) 
                        if(abs(rapidpsi).lt.1.d0) nbclose=nbclose+1
                     endif
                     if(sqrt(dx_part2(1)**2+dx_part2(2)**2+
     &                    dx_part2(3)**2+dx_part2(0)**2).gt.1D-3)then
                        nbcontrib(partcount2)=nbcontrib(partcount2)+1
                     endif
                     if(correctshift) then
                        do k=1,3
                           dp_cm(k)=dp_part(k)+dp_part2(k)
                           dp_part(k)=dp_part(k)-dp_cm(k)/2.d0
                           dp_part2(k)=dp_part2(k)-dp_cm(k)/2.d0
                        enddo
                     endif
                     drstore(partcount,0)=drstore(partcount,0)+
     &                    dx_part(0)
                     do k=1,3
                        dpstore(partcount,k)=dpstore(partcount,k)+
     &                       dp_part(k)
                        drstore(partcount,k)=drstore(partcount,k)+
     &                       dx_part(k)
                     enddo
                     drstore(partcount2,0)=drstore(partcount2,0)+
     &                    dx_part2(0)
                     do k=1,3
                        dpstore(partcount2,k)=dpstore(partcount2,k)+
     &                       dp_part2(k)
                        drstore(partcount2,k)=drstore(partcount2,k)+
     &                       dx_part2(k)
                     enddo
                  endif
               endif
               ipart2=partinfo_next(ipart2)
            enddo
         endif
         if(nbcontrib(partcount).gt.1) then
            nbcontribgt2=nbcontribgt2+1
            tbinst=sqrt((xMRU(0)+drstore(partcount,0))**2-
     &           (xMRU(3)+drstore(partcount,3))**2)
            if((tbinst.gt.tord0).and.(nbwarn.le.100)) then
               nbwarn=nbwarn+1
               write(ifmtx,*) 'part:',ipart,':',nbcontrib(partcount),
     &              'contrib'
               write(ifmtx,*) tbinst,'>',tord0
            endif
         endif
         if(minrrelpartner(ipart).ne.0) then
            ipart2=minrrelpartner(ipart)
            do k=0,3
               minrrel(ipart,k)=partinfo_r(ipart2,k)-partinfo_r(ipart,
     &              k)
            enddo
         endif
         ipart=partinfo_next(ipart)
      enddo
      ipart=firstpart
      partcount=0
      nbhqbis=0
      if(itypeorder.eq.1) time=tord0
      do while (ipart.NE.-1)
         partcount=partcount+1         
         ires=partinfo_ires(ipart)   
         if((ires.ge.1).and.(ires.le.26).and.(ires.ne.3).and.
     &        (ires.ne.16)) nbhqbis=nbhqbis+1
         if(itypeorder.eq.2) time=gtimefin2(ipart,tord0)
         do k=0,3
            p_part_old(k)=partinfo_p(ipart,k)
            p_part(k)=p_part_old(k)
         enddo   
         dt=time-partinfo_r(ipart,0)
         dtau=dt/p_part(0)
         do k=1,3
            x_part_old(k)=partinfo_r(ipart,k)
            x_part(k)=x_part_old(k)+p_part(k)*dtau
         enddo
         x_part_old(0)=partinfo_r(ipart,0)
         x_part(0)=time
         if(if2bodySC.eq.1) then
            do k=1,3
               if(.not.isnan(dpstore(partcount,k))) then
                  p_part(k)=p_part(k)+dpstore(partcount,k)
               else
                  write(ifmtx,*) 'isnan in 1coll evol'
                  close(ifmtx)
                  stop
               endif
               if(.not.isnan(drstore(partcount,k))) then
                  x_part(k)=x_part(k)+drstore(partcount,k)
               else
                  write(ifmtx,*) 'isnan in 1coll evol'
                  close(ifmtx)
                  stop
               endif
            enddo
            x_part(0)=x_part(0)+drstore(partcount,0)
            p_part(0)=sqrt(partinfo_mass(ipart)**2+p_part(1)**2+
     &           p_part(2)**2+p_part(3)**2)
         endif
         call givelocalhydro(x_part,temp,beta,densener,ifail_med,imed)
         if(imed.lt.0) then
            write(ifmtx,*) 'update_part:part. should be at least in QGP'
            close(ifmtx)
            stop
         endif   
         if((imed.ge.0).and.(imed.le.imedfinal)) then
            if(partinfo_whichmedium(ipart).gt.imedfinal) then
               if(resinfo_ifquark(ires)) then
                  write(ifmtx,*) 'stupid in update:no quark can reenter'
                  close(ifmtx)
                  stop
               endif
            endif
            if(.not.resinfo_ifquark(ires)) then
               if(temp.gt.tempdissoc(ires)) then 
                  call melt(ipart)
               endif
            else
                if(reddof<3) then
                   effdeg=reduction_dof(1,10.0d0,temp)
                else
                   effdeg=1.0
                endif
                if(isce<7) then
                   if((imed.eq.0).or.(.not.(ifdecrease))) then
                      effdeg=1.d0
                   else
                      effdeg=densener/entrans_max
                   endif
                endif
               if(type_evol_Q.eq.1) then
                  call kickme_lang(p_part,partinfo_mass(ipart),dt,
     &            resinfo_ifquark(ires),ires,temp,effdeg,beta,ifaillang)  
               else
                  call kickme_bolt(p_part,dt,
     &             resinfo_ifquark(ires),ires,temp,effdeg,beta,deltp,
     &                 ifaillang,ifinfo)
                  accept=.false.
                  absdeltap=sqrt(deltp(1)**2+deltp(2)**2+deltp(3)**2)
                  if(absdeltap.gt.1D-9) then
                     if(ifsingletred.and.(minrrelpartner(ipart).gt.0)) 
     &                    then
                        coherefact=deltp(0)*minrrel(ipart,0)
                        do k=1,3
                           coherefact=coherefact-deltp(k)*
     &                          minrrel(ipart,k)
                        enddo       
                        coherefact=coherefact/0.197
                        if(ran2().lt.tanh(coherefact)**2) then
                           accept=.true.
                        else
                        endif
                     else
                        accept=.true.
                     endif
                  endif
                  if(accept) then
                     do k=0,3
                        p_part(k)=p_part(k)+deltp(k)
                     enddo
                  endif
               endif
               p_part(0)=sqrt(partinfo_mass(ipart)**2+p_part(1)**2+
     &           p_part(2)**2+p_part(3)**2)
               do k=0,3
                  partinfo_p(ipart,k)=p_part(k)
               enddo
            endif                
         endif  
         if((imed.gt.imedfinal).and.
     &       (partinfo_whichmedium(ipart).le.imedfinal)) then
            if(resinfo_ifquark(ires)) then
               call givelocalhydro(x_part_old,temp,beta,densener,
     &              ifail_medbid,imedbid)
               call hadronize(ipart,itypDproduct,beta)
            endif
         endif  
         partinfo_whichmedium(ipart)=imed
         do k=0,3
            partinfo_r(ipart,k)=x_part(k)
         enddo
         partinfo_tb(ipart)=sqrt(x_part(0)**2-x_part(3)**2)
         ipart=partinfo_next(ipart)
      enddo
         deallocate(ifquarkQGP)
         deallocate(imedsave)
         deallocate(nbcontrib)
         deallocate(tempsave)
         deallocate(drstore)
         deallocate(dpstore)
         deallocate(minrrel)
         deallocate(minrrelinv)
         deallocate(minrrelpartner)
      return
      end

      subroutine kickme_bolt(p,dt,ifquark,ires,temp,effdeg,beta,
     &                       deltap,ifail,ifinfo)
      use JA
      use PBGgenvar
      use PBGforratedat1
      use PBGforprocesstype
      use PBGnumbereventHQ
      implicit none
      logical ifquark,ifinfo
      integer ires,i,ifail,itypq,nbcol
      double precision p(0:3),dt,temp,beta(3),u(0:3),pp(0:3),
     & dtp,ppvec(3),pnorm,effdeg,deltap(0:3),pin(0:3),pfin(0:3)
CJA


      ifail=0
      
      if((.not.ifquark).or.(ratemult.le.0.)) then
         do i=0,3
            deltap(i)=0.d0
         enddo
         return
      endif
      do i=0,3
         pin(i)=p(i)
      enddo
CJA
CJA
      itypq=itypquark(ires)
      u(0)=1./sqrt(1-(beta(1)**2+beta(2)**2+beta(3)**2))
      u(1)=u(0)*beta(1)
      u(2)=u(0)*beta(2)
      u(3)=u(0)*beta(3)
CJA
CJA
      call lorentz2(u,p,pp)
      dtp=pp(0)/p(0)*dt
      if(ifinfo) then
         write(6,*) 'u0:',u(0)
         write(6,*) 'p0:',p(0)
         write(6,*) 'pp0:',pp(0)
         write(6,*) 'dtp:',dtp
      endif
      ppvec(1)=pp(1)
      ppvec(2)=pp(2)
      ppvec(3)=pp(3)
      pnorm=sqrt(pp(1)**2+pp(2)**2+pp(3)**2)
      if(pnorm.gt.pmax(itypq)) then
         ifail=1
      endif   
      call evolveboltzmann(itypq,pp,temp,effdeg,dtp,ifcol,ifrad,nbcol,
     & ifinfo)
      do i=1,3
         u(i)=-u(i)
      enddo
      call lorentz2(u,pp,pfin)
      do i=0,3
         deltap(i)=pfin(i)-pin(i)
      enddo
CJA
      if(abs(pin(0)-pfin(0)).gt.0.001) then
      jancoll=jancoll+1
      if(jahqcoll.eq.1)
     &write(73,555)icol,jancoll,jaipart,ires,jaityp,pfin(0),pfin(1),
     &pfin(2),pfin(3),pin(0),
     &pin(1),pin(2),pin(3),partinfo_r(jaipart,0),
     &partinfo_r(jaipart,1),
     &partinfo_r(jaipart,2),partinfo_r(jaipart,3),temp,beta(1),beta(2),
     &beta(3)
555   Format(1x,2i5,3i4,18f8.3)
      endif
CJA
      return
      end


      subroutine kickme_lang(p,mass,dt,ifquark,ires,temp,effdeg,beta,
     &                       ifail)
      use PBGgenvar
      use PBGforfpdat1
      implicit none
      logical ifquark
      integer ires,i,ifail,itypq
      double precision p(0:3),mass,dt,temp,beta(3),u(0:3),pp(0:3),dtp,
     & ppvec(3),deltap(3),pnorm,effdeg
      ifail=0
      if((.not.ifquark).or.(fpmult.le.0.)) return
      itypq=itypquark(ires)
      u(0)=1./sqrt(1-(beta(1)**2+beta(2)**2+beta(3)**2))
      u(1)=u(0)*beta(1)
      u(2)=u(0)*beta(2)
      u(3)=u(0)*beta(3)
      call lorentz2(u,p,pp)
      dtp=pp(0)/p(0)*dt
      ppvec(1)=pp(1)
      ppvec(2)=pp(2)
      ppvec(3)=pp(3)
      pnorm=sqrt(pp(1)**2+pp(2)**2+pp(3)**2)
      if(pnorm.gt.pfmax(itypq)) then
         ifail=1
      endif   
      call gimmedplangevin(itypq,ppvec,temp,effdeg,dtp,deltap)
      pp(0)=mass**2
      do i=1,3
         pp(i)=pp(i)+deltap(i)
         pp(0)=pp(0)+pp(i)**2
      enddo
      pp(0)=sqrt(pp(0))
      do i=1,3
         u(i)=-u(i)
      enddo
      call lorentz2(u,pp,p)
      return
      end

      subroutine gimmedplangevin(itypq,p,temp,effdeg,dt,dp)
      implicit none
      integer i,itypq
      double precision p(3),temp,effdeg,dt,dp(3),normp,ran2,adrag,blong,
     & btrans,dpl,dpt,pi,pperp(3),pperp2(3),phi,phi2,arad,blrad,btrad,
     & dplrad,dptrad,rap
      parameter (pi=3.14159265358979323844d0)      
      normp=sqrt(p(1)**2+p(2)**2+p(3)**2)
      call gimmefpcoeff(itypq,normp,temp,adrag,blong,btrans)
      adrag=effdeg*adrag 
      blong=effdeg*blong
      btrans=effdeg*btrans
      arad=effdeg*arad 
      blrad=effdeg*blrad
      btrad=effdeg*btrad
      dpl=-adrag*dt+sqrt(-4*dt*blong*log(ran2()))*cos(pi*ran2())
      dplrad=0.d0
      dpt=sqrt(-4*dt*btrans*log(ran2()))
      dptrad=0.d0
      if(normp.lt.1d-16) then
         dp(3)=dpl+dplrad
         phi=2*pi*ran2()
         phi2=2*pi*ran2()
         dp(1)=dpt*cos(phi)+dptrad*cos(phi2)
         dp(2)=dpt*sin(phi)+dptrad*sin(phi2)
      else
         call giveperpvec(p,pperp)
         call giveperpvec(p,pperp2)
         rap=(dpl+dplrad)/normp
         do i=1,3
            dp(i)=rap*p(i)+dpt*pperp(i)+dptrad*pperp2(i)
         enddo   
      endif
      return
      end

      subroutine gimmefpcoeff(itypq,p0,temp0,adrag,blong,btrans)
      use PBGforfpdat1
      use PBGforfpdat2
      implicit none
      double precision p0,temp0,p,temp,adrag,blong,btrans
      integer itypq,ipinf,ipsup,itinf,itsup
      double precision prap,dp,trap,dtemp,aii,asi,ass,ais
      adrag=0.
      blong=0.
      btrans=0.
      if(fpmult.le.0) return
      p=p0
      temp=temp0
      if(p.gt.pfmax(itypq)) p=0.95*pfmax(itypq)
      if(temp.gt.tempfmax(itypq)) temp=tempfmax(itypq)
      if(temp.lt.tempfmin(itypq)) return
      if(temp.eq.tempfmax(itypq)) then
         if(p.eq.pfmax(itypq))then
            adrag=fpa(itypq,nbp(itypq),nbtemp(itypq))
            blong=fpbl(itypq,nbp(itypq),nbtemp(itypq))
            btrans=fpbt(itypq,nbp(itypq),nbtemp(itypq))
         else
            prap=p/deltp(itypq)
            ipinf=int(prap)
            ipsup=ipinf+1
            dp=prap-ipinf
            adrag=dp*fpa(itypq,ipsup,nbtemp(itypq))+
     &            (1-dp)*fpa(itypq,ipinf,nbtemp(itypq))
            blong=dp*fpbl(itypq,ipsup,nbtemp(itypq))+
     &            (1-dp)*fpbl(itypq,ipinf,nbtemp(itypq))
            btrans=dp*fpbt(itypq,ipsup,nbtemp(itypq))+
     &            (1-dp)*fpbt(itypq,ipinf,nbtemp(itypq))
         endif
      else 
         trap=(temp-tempfmin(itypq))/delttemp(itypq)
         itinf=int(trap)
         itsup=itinf+1
         dtemp=trap-itinf
         if(p.eq.pfmax(itypq)) then
            adrag=dtemp*fpa(itypq,nbp(itypq),itsup)+
     &            (1-dtemp)*fpa(itypq,nbp(itypq),itinf)
            blong=dtemp*fpbl(itypq,nbp(itypq),itsup)+
     &            (1-dtemp)*fpbl(itypq,nbp(itypq),itinf)
            btrans=dtemp*fpbt(itypq,nbp(itypq),itsup)+
     &            (1-dtemp)*fpbt(itypq,nbp(itypq),itinf)
         else
            prap=p/deltp(itypq)
            ipinf=int(prap)
            ipsup=ipinf+1
            dp=prap-ipinf
            aii=(1-dp)*(1-dtemp)
            ais=(1-dp)*dtemp
            asi=dp*(1-dtemp)
            ass=dp*dtemp
            adrag=aii*fpa(itypq,ipinf,itinf)
     &            +ais*fpa(itypq,ipinf,itsup)+
     &            asi*fpa(itypq,ipsup,itinf)+ass*fpa(itypq,ipsup,itsup)
            blong=aii*fpbl(itypq,ipinf,itinf)
     &            +ais*fpbl(itypq,ipinf,itsup)+
     &          asi*fpbl(itypq,ipsup,itinf)+ass*fpbl(itypq,ipsup,itsup)
            btrans=aii*fpbt(itypq,ipinf,itinf)
     &             +ais*fpbt(itypq,ipinf,itsup)+
     &          asi*fpbt(itypq,ipsup,itinf)+ass*fpbt(itypq,ipsup,itsup)
         endif
      endif
      adrag=adrag*fpmult
      blong=blong*fpmult
      btrans=btrans*fpmult
      return
      end



      subroutine gimmewignerrelatsinglet(pQ,pQbar,xQ,xQbar,sigma,
     & wignerrelat,pT,phi,rapid,xrelnormcm,prelnormcm,gamma,ergcm)
      use PBGpsiinfo
      implicit none 
      double precision wignerrelat,wignerGxsinglet
      double precision pT,rapid,mT
      double precision pQ(0:3),pQbar(0:3),xQ(0:3),xQbar(0:3)
      double precision s,ptot(0:3),u4(0:3),prel(0:3),mQQbar,
     & gamma,prelnormcm,xrel(0:3),xrel2,xrelnormcm,xrelu,
     & prel2,prelu,massjsi,uT,phi,pi,sigma,ergcm
      parameter (massjsi=3.09,pi=3.141592653589793)
      integer i
      integer ifmtx
      call getMonitorFileIndex(ifmtx)
      if((pQ(0).gt.1D8).or.(pQbar(0).gt.1D8)) then  
         wignerrelat=0.d0
         pT=1D8 
         phi=0.d0
         rapid=0.d0
         xrelnormcm=1D5
         prelnormcm=1D5
         gamma=1D5
         ergcm=1D5
         return
      endif
      do i=0,3
         ptot(i)=pQ(i)+pQbar(i)
         prel(i)=(pQ(i)-pQbar(i))/2
         xrel(i)=xQ(i)-xQbar(i)
      enddo
      s=ptot(0)**2-(ptot(1)**2+ptot(2)**2+ptot(3)**2)
      if(s.le.0.d0) then
         wignerrelat=0.d0
         pT=1D8 
         phi=0.d0
         rapid=0.d0
         xrelnormcm=1D5
         prelnormcm=1D5
         gamma=1D5
         ergcm=1D5         
         return
      endif
      if(isnan(mQQbar)) then
         write(ifmtx,*) 'NAN invariant in mass gimmewignerrelatsinglet'
         close(ifmtx)
         stop
      endif
      mQQbar=sqrt(s)
      mT=sqrt(ptot(0)**2-ptot(3)**2)
      if(ptot(3).gt.0) then
         rapid=log((ptot(0)+ptot(3))/mT)
      else
         rapid=-log((ptot(0)-ptot(3))/mT)
      endif
      do i=0,3
         u4(i)=ptot(i)/mQQbar
      enddo
      gamma=u4(0)
      uT=sqrt(u4(1)**2+u4(2)**2)
      pT=uT*massjsi
      ergcm=mQQbar
      phi=atan2(u4(2),u4(1))+pi
      xrel2=xrel(0)**2-(xrel(1)**2+xrel(2)**2+xrel(3)**2)
      prel2=prel(0)**2-(prel(1)**2+prel(2)**2+prel(3)**2)
      xrelu=xrel(0)*u4(0)-xrel(1)*u4(1)-xrel(2)*u4(2)-xrel(3)*u4(3)
      prelu=prel(0)*u4(0)-prel(1)*u4(1)-prel(2)*u4(2)-prel(3)*u4(3)
      xrelnormcm=sqrt(xrelu**2-xrel2)
      prelnormcm=sqrt(prelu**2-prel2)
      wignerrelat=wignerGxsinglet(xrelnormcm,prelnormcm,sigma)
      return
      end

      double precision function wignerGxsinglet(r,p,sigma)
      use PBGpsiinfo
      implicit none
      double precision r,p,miu,hbc,g,sigma,arg
      parameter (g=3.d0/4.d0,hbc=0.1973)
      miu=hbc/sigma
      arg=(r/sigma)**2+(p/miu)**2
      if(arg.lt.30.d0) then
         wignerGxsinglet=8*g*exp(-arg)

      else
         wignerGxsinglet=0.d0
      endif
      return      
      end function 


      subroutine oneevolQQbarbis(mQ,x01,x02,p01,p02,
     & ifevol,ifdisplay,x01fin,p01fin,etotcmfin,rrelinv,
     & delta_x4Q,delta_p4Q,xMRU,ifail)
      use PBGmainpro2
      implicit none
      logical ifdisplay
      integer i,ifail,ifevol,nbwarns
      save nbwarns
      data nbwarns/0/
      double precision mq,x01(0:3),x02(0:3),p01(0:3),p02(0:3),
     & x01fin(0:3),p01fin(0:3),delta_x4Q(0:3),delta_p4Q(0:3),
     & rQcmsup,xMRU(0:3),etotcmfin,rrelinv
      parameter(rQcmsup=1.d0)
      integer ifmtx
      call getMonitorFileIndex(ifmtx)
      if(ifdisplay) write(ifmtx,*) 'entering oneQQbarevolvebis' 
      ifail=0
      do i=0,3
         if(isnan(x01(i))) then
            write(ifmtx,*) 'x01 ISNAN',x01
            ifail=1
         endif
      enddo
      do i=0,3
         if(isnan(x02(i))) then
            write(ifmtx,*) 'x02 ISNAN',x02
            ifail=2
         endif
      enddo
      do i=0,3
         if(isnan(p01(i))) then
            write(ifmtx,*) 'p01 ISNAN',p01
            ifail=3
         endif
      enddo
      do i=0,3
         if(isnan(p02(i))) then
            write(ifmtx,*) 'p02 ISNAN',p02
            ifail=4
         endif
      enddo
      if(ifail.gt.0) then
         write(ifmtx,*) 'existing'
         return
      endif
      etotcmfin=2*mQ
      rrelinv=1/mQ
      if(ifevol.eq.1) then
         write(ifmtx,*) 'NO 2-body implemented in this version'
         nbwarns=nbwarns+1
         if(nbwarns.ge.20) then
            close(ifmtx)
            stop
         endif
      endif
      do i=0,3
         x01fin(i)=xMRU(i)
         p01fin(i)=p01(i)
         delta_x4Q(i)=0.d0
         delta_p4Q(i)=0.d0
      enddo
      return
      end 
      
