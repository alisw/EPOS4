      double precision function cos2rand()
      implicit none
      double precision x,x2,y,r2,ran2
 401  x=ran2()
      y=ran2()
      x2=x**2
      r2=x2+y**2
      if(r2.gt.1.) goto 401
      cos2rand=x2/r2
      return
      end


      subroutine gauss(g1,g2)
      implicit none
      double precision g1,g2,x,y,r2,r,ran2
 401  x=1.d0-2*ran2()
      y=1.d0-2*ran2()
      r2=x*x+y*y
      if(r2.gt.1.0) goto 401
      r=sqrt(r2)
      g1=(x/r)*sqrt(-2*log(r2))
      g2=(y/x)*g1
      return
      end


      subroutine twoDnormaldeviates(sigmax,sigmay,sigma,fx,fy)
      implicit none
      real sigmax,sigmay,sigma
      double precision r,rmax,rho,rho2,theta,fx,fy,ran2,pi
      parameter (pi=3.14159265358979323844d0)  
      rmax=dble(sigmax*sigmay/(4.0d0*pi*sigma))
 402  r=rmax*ran2()
      if(r.eq.0) goto 402
      rho2=-2.d0*dlog(4.0d0*pi*dble(sigma)*r/dble(sigmax*sigmay))
      r=ran2()
      if(r.gt.0.5d0)then
         rho=dsqrt(rho2)
      else
         rho=-dsqrt(rho2)
      endif
      theta=pi*ran2()
      fx=dble(sigmax*rho*dcos(theta))
      fy=dble(sigmay*rho*dsin(theta))
      return
      end


      subroutine twopibdeviates(b1,b2,b)
      implicit none
      double precision ran2
      real b1,b2,b
      b=nint(sqrt(ran2()*(b2**2-b1**2)+b1**2))
      return
      end


      subroutine classify_impact(binf,bsup,nbcol,nbsubrun,nbimpact,
     &                          nbevent)
      implicit none 
      integer nbcol,nbsubrun,nbimpact,nbsubrunmax,nbimpactmax
      parameter (nbimpactmax=15,nbsubrunmax=100)
      integer nbevent(nbsubrunmax,0:nbimpactmax),
     & nbeventall(0:nbimpactmax),
     & ibinf,ibsup,ib,checksum,i,isr,npersub,remain 
      real binf,bsup,propor_event(0:nbimpactmax)
      double precision ran2
        write(*,*) 'nbimpact: ',nbimpact
      ibinf=nint(binf)
      ibsup=nint(bsup)
      if((nbsubrun.gt.nbsubrunmax).or.(nbimpact.gt.nbimpactmax).or.
     &    (ibsup.gt.nbimpact)) then
         write(8,*) 'problem in classify_impact'
         close(8)
         stop
      endif
      do ib=0,nbimpact
         propor_event(ib)=0.
      enddo
      if(ibinf.eq.ibsup) then
         propor_event(ibinf)=1.
      else
         propor_event(ibinf)=(ibinf+0.5)**2-binf**2
         do ib=ibinf+1,ibsup-1
            propor_event(ib)=(ib+0.5)**2-(ib-0.5)**2
         enddo
         propor_event(ibsup)=bsup**2-(ibsup-0.5)**2
         do ib=0,nbimpact
            propor_event(ib)=propor_event(ib)/(bsup**2-binf**2)
         enddo
      endif
      checksum=0
      do ib=0,nbimpact
         nbeventall(ib)=int(nbcol*propor_event(ib))
         checksum=checksum+nbeventall(ib)
      enddo
      do i=checksum+1,nbcol
         ib=nint(sqrt(ran2()*(bsup**2-binf**2)+binf**2))
         nbeventall(ib)=nbeventall(ib)+1
      enddo
      do ib=0,nbimpact
         write(8,*) 'b=',ib,'->',nbeventall(ib),' events'
      enddo
      do ib=0,nbimpact
         npersub=nbeventall(ib)/nbsubrun
         remain=nbeventall(ib)-nbsubrun*npersub
         do isr=1,nbsubrun
            nbevent(isr,ib)=npersub
         enddo
         nbevent(nbsubrun,ib)=nbevent(nbsubrun,ib)+remain
      enddo
      return
      end


      double precision function twoDgaussian(sigmax,sigmay,sigma,x,y)
      implicit none
      real sigmax,sigmay,sigma
      double precision x,y,pi
      parameter (pi=3.14159265358979323844d0) 
      twoDgaussian=dexp(-0.5d0*(x**2/dble(sigmax)**2+y**2/
     &             dble(sigmay)**2))/
     &             (2.0d0*pi*dble(sigma))
      return
      end


      subroutine giveperpvec(p,porth)
      implicit none
      integer i
      double precision p(3),porth(3),pt,normp2,normp,pi,proj,phi,ran2,
     & normporth
      parameter (pi=3.14159265358979323844d0)      
      normp2=p(1)**2+p(2)**2+p(3)**2
      normp=sqrt(normp2)
 10   porth(3)=2*ran2()-1
      pt=sqrt(1-porth(3)**2)
      phi=2*pi*ran2()
      porth(1)=pt*cos(phi)
      porth(2)=pt*sin(phi)
      if(normp.eq.0.d0) return
      proj=porth(1)*p(1)+porth(2)*p(2)+porth(3)*p(3)
      if(1-abs(proj).le.1d-4*normp) goto 10
      proj=proj/normp2
      normporth=0.
      do i=1,3
         porth(i)=porth(i)-proj*p(i)
         normporth=normporth+porth(i)**2
      enddo
      normporth=1./sqrt(normporth)
      do i=1,3
         porth(i)=normporth*porth(i)
      enddo
      proj=porth(1)*p(1)+porth(2)*p(2)+porth(3)*p(3)
      if(abs(proj).gt.1d-8) then
         write(8,*) 'problem in giveperpvec: not orthogonal'
         write(8,*) p(1),p(2),p(3)
         write(8,*) porth(1),porth(2),porth(3)
         close(8)
         stop
      endif
      return
      end
      

      double precision function gtimefin(t0,z0,en,pz,mt2,tbfin)
      implicit none
      double precision t0,z0,tbfin,tbin2,tbin,en,mt2,pz,a
      tbin2=t0**2-z0**2
      if(tbin2.lt.0) then 
         write(8,*) 'error in gtimefin: tbin2<0'
         close(8)
         stop
         return
      endif
      tbin=sqrt(tbin2)
      if(tbfin.lt.tbin) then
         write(8,*) 'error in gtimefin: part. Bkorken time ',
     &               tbin,' > ',tbfin
         write(8,*) t0,z0,en,pz,mt2
         close(8)
         stop
         return
      else
         a=(en*z0-pz*t0)
         gtimefin=(pz*a+en*sqrt(a*a+tbfin*tbfin*mt2))/mt2
      endif
      return
      end


      double precision function gtimefin2(ipart,tbfin)
      use PBGgenvar
      implicit none
      integer ipart
      double precision tbfin,tbin2,tbin,mt2,a
      tbin2=partinfo_r(ipart,0)**2-partinfo_r(ipart,3)**2
      if(tbin2.lt.0) then 
         write(6,*) 'error in gtimefin2: tbin2<0'
         stop
      endif
      tbin=sqrt(tbin2)
      if(tbfin.lt.tbin) then
         write(6,*) 'In gtimefin2: part.',ipart,
     &    ' has already Bkorken time ',tbin,' > asked ',tbfin
         write(6,*) partinfo_r(ipart,0),partinfo_r(ipart,3),
     &      partinfo_tb(ipart)
         stop
      else
         a=partinfo_p(ipart,0)*partinfo_r(ipart,3)-
     &     partinfo_p(ipart,3)*partinfo_r(ipart,0)
         mt2=partinfo_mass(ipart)**2+partinfo_p(ipart,1)**2+
     &       partinfo_p(ipart,2)**2
         gtimefin2=(partinfo_p(ipart,3)*a+partinfo_p(ipart,0)*
     &              sqrt(a*a+tbfin*tbfin*mt2))/mt2
      endif
      return
      end


      subroutine lorentz2(u,p,pp)
      implicit none
      double precision proj,u(0:3),p(0:3),pp(0:3)
      pp(0)=p(0)*u(0)-p(1)*u(1)-p(2)*u(2)-p(3)*u(3)      
      proj=(u(1)*p(1)+u(2)*p(2)+u(3)*p(3))/(1+u(0))-p(0) 
      pp(1)=p(1)+proj*u(1)
      pp(2)=p(2)+proj*u(2)
      pp(3)=p(3)+proj*u(3)
      return
      end


      subroutine boost(u,p,pp)
      implicit none
      double precision proj,u(0:3),p(0:3),pp(0:3)
      proj=p(1)*u(1)+p(2)*u(2)+p(3)*u(3)
      pp(0)=p(0)*u(0)+proj      
      proj=proj/(1+u(0))+p(0) 
      pp(1)=p(1)+proj*u(1)
      pp(2)=p(2)+proj*u(2)
      pp(3)=p(3)+proj*u(3)
      return
      end


      SUBROUTINE midexp(funk,aa,s,n)
      implicit none
      integer n
      double precision aa,s,funk
      external funk
      integer it,j
      double precision ddel,del,sum,tnm,x,func,a,b
      func(x)=funk(-log(x))/x
      b=exp(-aa)
      a=0
      if (n.eq.1) then
         s=(b-a)*func(0.5d0*(a+b))
      else
         it=3**(n-2)
         tnm=it
         del=(b-a)/(3.d0*tnm) 
         ddel=del+del
         x=a+0.5d0*del
         sum=0.d0
         do j=1,it
            sum=sum+func(x)
            x=x+ddel
            sum=sum+func(x)
            x=x+del
         enddo
         s=(s+(b-a)*sum/tnm)/3.d0
      endif
      return
      end


      subroutine dpolint(xa,ya,n,x,y,dy)
      implicit none
      integer n,nmax,i,m,ns
      parameter (nmax=10)
      double precision dy,x,y,xa(n),ya(n),den,dif,dift,ho,hp,w,
     & c(nmax),d(nmax)
      ns=1
      dif=abs(x-xa(1))
      do i=1,n
        dift=abs(x-xa(i))
        if(dift.lt.dif)then
          ns=i
          dif=dift
        endif
        c(i)=ya(i)
        d(i)=ya(i)
      enddo
      y=ya(ns)
      ns=ns-1 
      do m=1,n-1
        do i=1,n-m
          ho=xa(i)-x
          hp=xa(i+m)-x
          w=c(i+1)-d(i)
          den=ho-hp
          if(den.eq.0.) then
             write(*,*) 'failure in polint'
             read(*,'()')
          endif
          den=w/den
          d(i)=hp*den
          c(i)=ho*den
        enddo
        if (2*ns.lt.n-m) then
          dy=c(ns+1)
        else
          dy=d(ns)
          ns=ns-1
        endif
        y=y+dy
      enddo
      return
      end  


      subroutine qromb(func,a,b,ss,eps)
      implicit none
      integer jmax,jmaxp,k,km
      double precision a,b,func,ss,eps
      external func
      parameter (jmax=16, jmaxp=jmax+1, k=5, km=k-1)
      integer j
      double precision dss,h(jmaxp),s(jmaxp)
      h(1)=1.
      do j=1,jmax
        call dtrapzd(func,a,b,s(j),j)
        if(j.ge.k) then
          call dpolint(h(j-km),s(j-km),k,0.D0,ss,dss)
          if (abs(dss).le.eps*abs(ss)) return
        endif
        s(j+1)=s(j)
        h(j+1)=0.25d0*h(j)
      enddo
      write(*,*) 'too many steps in qromb'
      read(*,'()')
      end       


      subroutine dtrapzd(func,a,b,s,n)
      implicit none
      integer n
      real*8 a,b,s,func
      external func
      integer it,j
      real*8 del,sum,tnm,x
      if(n.eq.1) then
        s=0.5d0*(b-a)*(func(a)+func(b))
      else
        it=2**(n-2)
        tnm=it
        del=(b-a)/tnm
        x=a+0.5d0*del
        sum=0.d0
        do j=1,it
          sum=sum+func(x)
          x=x+del
        enddo
        s=0.5d0*(s+(b-a)*sum/tnm)
      endif
      return
      end
      
      subroutine intrapx(func,b,A,c,d,s,n)
      implicit none
      integer n,k,A
      real b
      double precision c,d,s,h,func
      external func
      h=(d-c)/n
      s=0d0
      do k=0,n-1
         s=s+func(c+k*h+b/2,A)*func(c+k*h-b/2,A)+
     &       func(c+(k+1)*h+b/2,A)*func(c+(k+1)*h-b/2,A)
      enddo
      s=0.5d0*h*s
      return
      end

      subroutine intrapy(func,x,A,c,d,s,n)
      implicit none
      integer n,k,A
      double precision c,d,s,h,func,x
      external func
      h=(d-c)/n
      s=0d0
      do k=0,n-1
         s=s+func(x,c+k*h,A)+func(x,c+(k+1)*h,A)
      enddo
      s=0.5d0*h*s
      return
      end

      subroutine intrapz(func,x,y,A,c,d,s,n)
      implicit none
      integer n,k,A
      double precision c,d,s,h,func,x,y
      external func
      h=(d-c)/n
      s=0d0
      do k=0,n-1
         s=s+func(x,y,c+k*h,A)+func(x,y,c+(k+1)*h,A)
      enddo
      s=0.5d0*h*s
      return
      end
      subroutine qtrapimproper(func,a,s,prec)
      implicit none
      integer jmax
      double precision a,func,s,prec
      external func
      parameter(jmax=20)
      integer j
      double precision olds
      olds =-1e30
      do j=1,jmax
        call midexp(func,a,s,j)
        if(abs(s-olds).lt.prec*abs(olds)) return
        olds=s
      enddo
      write(*,*) 'too many steps in qtrap'
      read(*,'()')
      end
 

      double precision function ran2()
      use PBGran2junk
      implicit none
      integer IM1,IM2,IMM1,IA1,IA2,IQ1,IQ2,IR1,IR2,NDIV
      double precision EPS,RNMX,AM
      parameter (IM1=2147483563,IM2=2147483399,AM=(1.0d0/IM1),
     &           IMM1=IM1-1)
      parameter (IA1=40014,IA2=40692,IQ1=53668,IQ2=52774,IR1=12211)
      parameter (IR2=3791,NDIV=(1+IMM1/NTAB))
      parameter (EPS=1.2d-7,RNMX=(1.0d0-EPS))
      integer J,K
      if(IDUM.le.0) then
         IDUM=max(-IDUM,1)
         IDUM2=IDUM
         do J=NTAB+8,1,-1
            K=(IDUM)/IQ1
            IDUM=IA1*(IDUM-K*IQ1)-K*IR1
            if(IDUM.lt.0) IDUM =IDUM+IM1
            if(J.le.NTAB) IV(J)=IDUM
         enddo
         IY=IV(1)
      endif
      K=(IDUM)/IQ1
      IDUM=IA1*(IDUM-K*IQ1)-K*IR1
      if(IDUM.lt.0) IDUM=IDUM+IM1
      K=IDUM2/IQ2
      IDUM2=IA2*(IDUM2-K*IQ2)-K*IR2
      if(IDUM2.lt.0) IDUM2=IDUM2+IM2
      J=1+IY/NDIV
      IY=IV(J)-IDUM2
      IV(J)=IDUM
      if(IY.lt.1) IY=IY+IMM1
      ran2=DBLE(min(AM*IY,RNMX))
      return
      end


      integer function poisson(av)
      implicit none
      integer i
      double precision av,x,ran2
      if(av.eq.0) then
         poisson=0
      else
         x=0.0d0
         i=-1
         do while (x.lt.1.0d0)
            x=x-log(ran2())/av
            i=i+1
         enddo
         poisson=i
      endif
      return
      end


      subroutine rotate3d(pref,cth,pfin,vpfin)
      implicit none
      integer imindir,k
      double precision pref(3),dir(3),vpfin(3),pfin,cth,eps1(3),
     & eps2(3),sth,phi,ran2,pi,cphi,sphi,invnorm,mindir,chck
      parameter(pi=3.14159265358979323844d0)
      invnorm=1./sqrt(pref(1)**2+pref(2)**2+pref(3)**2)
      imindir=0
      mindir=1.d0
      do k=1,3
         dir(k)=pref(k)*invnorm         
         if(abs(dir(k)).lt.mindir) then
            imindir=k
            mindir=abs(dir(k))
         endif
         vpfin(k)=0.d0
      enddo
      do k=1,3
         eps1(k)=-dir(imindir)*dir(k)
      enddo 
      eps1(imindir)=eps1(imindir)+1.
      invnorm=1./sqrt(eps1(1)**2+eps1(2)**2+eps1(3)**2)
      eps1(1)=eps1(1)*invnorm
      eps1(2)=eps1(2)*invnorm
      eps1(3)=eps1(3)*invnorm
      chck=dir(1)*eps1(1)+dir(2)*eps1(2)+dir(3)*eps1(3)
      if(abs(chck).gt.1d-8) then
         write(6,*) 'error 1 in rotate3d'
         stop
      endif
      eps2(1)=dir(2)*eps1(3)-dir(3)*eps1(2)
      eps2(2)=dir(3)*eps1(1)-dir(1)*eps1(3)
      eps2(3)=dir(1)*eps1(2)-dir(2)*eps1(1)
      chck=eps2(1)**2+eps2(2)**2+eps2(3)**2-1 
      if(abs(chck).gt.1d-8) then
         write(6,*) 'error 2 in rotate3d'
         stop
      endif
      sth=sqrt(1-cth**2)
      phi=2*pi*ran2() 
      cphi=cos(phi)
      sphi=sin(phi)
      do k=1,3
         vpfin(k)=pfin*(cth*dir(k)+sth*(cphi*eps1(k)+sphi*eps2(k)))
      enddo   
      chck=vpfin(1)*dir(1)+vpfin(2)*dir(2)+vpfin(3)*dir(3) 
      if(abs(chck/pfin-cth).gt.1d-7) then
         write(6,*) 'error in rotate3d: rotation angle not realized'
         stop
       endif  
      return
      end

      
      subroutine shell_sort(n,a,b)
      implicit none
      integer n,b(n)
      double precision a(n)
      integer i,j,inc,vb
      double precision va
      inc=1
 801  inc=3*inc+1
      if(inc.le.n) goto 801
 802  continue
      inc =inc/3
      do i=inc+1,n
         va=a(i)
         vb=b(i)
         j=i
 803     if (a(j-inc).gt.va) then
            a(j)=a(j-inc)
            b(j)=b(j-inc)
            j=j-inc
            if(j.le.inc) goto 804
            goto 803
         endif
 804     a(j)=va
         b(j)=vb
      enddo
      if (inc.gt.1) goto 802
      return
      end

      subroutine wedge(a,b,awedgeb)
      implicit none
      double precision a(3),b(3),awedgeb(3)
      awedgeb(1)=a(2)*b(3)-a(3)*b(2) 
      awedgeb(2)=a(3)*b(1)-a(1)*b(3)
      awedgeb(3)=a(1)*b(2)-a(2)*b(1) 
      end

