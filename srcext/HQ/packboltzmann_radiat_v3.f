
 
      double precision function alphas4(q2)
      implicit none
      double precision pi,sqlambda3, q2inf3,q2sup3, invbeta0,thresh3,q2
      parameter(pi=3.14159265358979323844d0, sqlambda3=0.04d0,
     &  q2inf3=-0.0005299712436888622d0, q2sup3=0.1396153933725991d0,
     &  invbeta0=4*pi/(11-2), thresh3=1.117d0)
      if(q2.gt.q2sup3) then
         alphas4=invbeta0/dlog(q2/sqlambda3)
      else
         if(q2.lt.q2inf3) then
            alphas4=invbeta0*(0.5-atan(log(-q2/sqlambda3)/pi)/pi)
         else
            alphas4=thresh3
         endif
      endif
      return
      end 

      double precision function alphas(q2)
      implicit none
      double precision q2,alphas4
      alphas=alphas4(q2)
      end

      double precision function alphasmax()
      implicit none
      double precision thresh3
      parameter(thresh3=1.117)
      alphasmax=thresh3
      end

      double precision function muoftemp(temp)
      implicit none
      integer nf
      parameter (nf=3)
      double precision temp,alphas,lambda3,pi
      parameter (lambda3=0.2d0,pi=3.14159265358979323844d0)
      alphas=4*pi/((11-2)*2*log(2*pi*temp/lambda3))
      muoftemp=sqrt((1+nf/6.d0)*4*pi*alphas)*temp
      return
      end


      double precision function mdebeff(temp,mdebinit)
      implicit none
      integer nf
      double precision temp,mdebinit,pi,mdeb2,mdeb2prev,tol,alphas4
      parameter(pi=3.14159265358979323844d0,nf=3,tol=1d-8)
      logical accept
      accept=.false.
      mdeb2prev=mdebinit**2
      do while(.not.accept)
         mdeb2=4*pi*(1+nf/6.)*alphas4(mdeb2prev)*temp**2
         accept=(abs(mdeb2-mdeb2prev).lt.tol)
         mdeb2prev=mdeb2
      enddo
      mdebeff=sqrt(mdeb2)
      return
      end

      subroutine initdebyemass
      use PBGfordebyeeff
      implicit none
      integer itemp
      double precision temp,mdebinit,mdebeff
      integer ifmtx
      call getMonitorFileIndex(ifmtx)
      tempmin=0.1d0
      tempmax=2.0d0
      nbtemp=191
      if(nbtemp.gt.nbtempmax) then
         write(ifmtx,*) 'nbtemp > nbtempmax in initdebyemass'
         close(ifmtx)
         stop
      endif
      deltemp=(tempmax-tempmin)/(nbtemp-1)
      mdebinit=0.1
      do itemp=1,nbtemp
         temp=tempmin+(itemp-1)*deltemp
         tbdebye(itemp)=mdebeff(temp,mdebinit)
         mdebinit=tbdebye(itemp)
      enddo
      return
      end

      double precision function mdebeffinterp(temp)
      use PBGfordebyeeff
      implicit none
      integer itemp
      double precision temp,dtemp
      if((temp.gt.tempmax)) then
         write(8,*) 'temp out of bounds in mdebeffinterp:',temp
         write(8,*) 'tempmax und min in mdebeffinterp:',tempmax,tempmin
         write(8,*) 'temp set to',tempmax
          temp=tempmax         
      else if(temp.lt.tempmin) then
         write(8,*) 'temp out of bounds in mdebeffinterp:',temp
         write(8,*) 'tempmax und min in mdebeffinterp:',tempmax,tempmin
         write(8,*) 'temp set to',tempmin
         temp=tempmin
      endif
      if(temp.eq.tempmax) then
         mdebeffinterp=tbdebye(nbtemp)
      else
         dtemp=(temp-tempmin)/deltemp
         itemp=int(dtemp)
         dtemp=dtemp-itemp
         itemp=itemp+1
         mdebeffinterp=(1-dtemp)*tbdebye(itemp)+dtemp*tbdebye(itemp+1)
      endif
      return
      end

      double precision function msqferm(m,s,t,mu)
      implicit none
      integer nf
      double precision s,m2,m,t,mu,smin,pi2
      parameter (nf=3,pi2=9.869604401089359d0)
      m2=m**2      
      smin=s-m2
      msqferm=pi2*512*nf*(smin**2+2*m2*t+(smin + t)**2)/(t-mu**2)**2
      return
      end

      double precision function msqglu1(m,s,t,u,mu)
      implicit none
      double precision s,m2,m,t,u,mu,smin,umin,pi2
      parameter (pi2=9.869604401089359d0)
      m2=m**2
      smin=s-m2
      umin=u-m2
      msqglu1=-pi2*96*32*smin*umin/(t-mu**2)**2
      return
      end

      double precision function msqglu2(m,s,u)
      implicit none
      double precision s,m2,m,u,smin,umin,pi2
      parameter (pi2=9.869604401089359d0)
      m2=m**2
      smin=s-m2
      umin=u-m2
      msqglu2=pi2*96*64*(2*m2*(m2+s)-smin*umin)/(9*smin**2)
      return
      end

      double precision function msqglu3(m,s,u)
      implicit none
      double precision s,m2,m,u,smin,umin,pi2
      parameter (pi2=9.869604401089359d0)
      m2=m**2
      smin=s-m2
      umin=u-m2
      msqglu3=pi2*96*64*(2*m2*(m2+u)-smin*umin)/(9*umin**2)
      return
      end

      double precision function msqglu4(m,s,t,u)
      implicit none
      double precision s,m2,m,t,u,smin,umin,pi2
      parameter (pi2=9.869604401089359d0)
      m2=m**2
      smin=s-m2
      umin=u-m2
      msqglu4=-pi2*96*16*m2*(4*m2-t)/(9*smin*umin)
      return
      end

        double precision function msqglu5(m,s,t,u,mu)
      implicit none
      double precision s,m2,m,t,u,mu,smin,umin,pi2
      parameter (pi2=9.869604401089359d0)
      m2=m**2
      smin=s-m2
      umin=u-m2
      msqglu5=pi2*96*16*(s*u-2*m2*s+m2**2)/(smin*(mu**2-t))
      return
      end

      double precision function msqglu6(m,s,t,u,mu)
      implicit none
      double precision s,m2,m,t,u,mu,smin,umin,pi2
      parameter (pi2=9.869604401089359d0)
      m2=m**2
      smin=s-m2
      umin=u-m2
      msqglu6=pi2*96*16*(smin*umin+ m2*(s-u))/(umin*(mu**2-t))
      return
      end

      double precision function msqgluall(m,s,t,mu,msqtab)
      implicit none
      integer i
      double precision m,s,t,mu,u,msqglu1,msqglu2,msqglu3,msqglu4,
     &  msqglu5,msqglu6,msqtab(6)
      u=2*m**2-s-t
      msqtab(1)=msqglu1(m,s,t,u,mu)
      msqtab(2)=msqglu2(m,s,u)
      msqtab(3)=msqglu3(m,s,u)
      msqtab(4)=msqglu4(m,s,t,u)
      msqtab(5)=msqglu5(m,s,t,u,mu)
      msqtab(6)=msqglu6(m,s,t,u,mu)
      msqgluall=0
      do i=1,6
         msqgluall=msqgluall+msqtab(i)
      enddo
      return
      end

      double precision function msqglu6pos(m,s,t,u,mu)
      implicit none
      double precision s,m2,m,t,u,mu,smin,umin,pi2
      parameter (pi2=9.869604401089359d0)
      m2=m**2
      smin=s-m2
      umin=u-m2
      msqglu6pos=pi2*96*16*smin/(mu**2-t)
      return
      end

      double precision function msqglutest(m,s,t,mu)
      implicit none
      double precision m,s,t,mu,u,msqglu1,msqglu2,msqglu3,msqglu4,
     &  msqglu6pos
      u=2*m**2-s-t
      msqglutest=msqglu1(m,s,t,u,mu)+msqglu2(m,s,u)+
     &                 msqglu3(m,s,u)+msqglu4(m,s,t,u)+
     &                 msqglu6pos(m,s,t,u,mu)
      return
      end


      double precision function msqfermrunning0(m,s,t,mu,alphast)
      implicit none
      double precision m,s,t,mu,msqferm,alphast
      msqfermrunning0=msqferm(m,s,t,mu)*alphast**2
      end

      double precision function msqfermrunning1(m,s,t,mu)
      implicit none
      double precision m,s,t,mu,msqfermrunning0,alphas
      msqfermrunning1=msqfermrunning0(m,s,t,mu,alphas(-t))
      end

      double precision function msqfermrunning1b(t)
      use PBGtomsqfermrunning1b
      implicit none
      double precision t,msqfermrunning1
      msqfermrunning1b=msqfermrunning1(m,s,t,mu)
      end

      double precision function msqfermrunning2(m,s,v,mu)
      implicit none
      double precision m,s,v,t,mu,msqfermrunning1
      t=(s-m**2)**2*v/(2*s)
      msqfermrunning2=msqfermrunning1(m,s,t,mu)
      end

      double precision function msqfermrunning2b(v)
      use PBGtomsqfermrunning2b
      implicit none
      double precision v,msqfermrunning2
      msqfermrunning2b=msqfermrunning2(m,s,v,mu)
      end


      double precision function msqglu1running0(m,s,t,u,mu,alphast)
      implicit none
      double precision s,m,t,u,mu,msqglu1,alphast
      msqglu1running0=msqglu1(m,s,t,u,mu)*alphast**2
      end

      double precision function msqglu1running1(m,s,t,u,mu)
      implicit none
      double precision m,s,t,u,mu,msqglu1running0,alphas
      msqglu1running1=msqglu1running0(m,s,t,u,mu,alphas(-t))
      end

      double precision function msqglu1running1b(t)
      use PBGtomsqglurunning1b
      implicit none
      double precision t,u,msqglu1running1
      u=2*m**2-s-t
      msqglu1running1b=msqglu1running1(m,s,t,u,mu)
      end

      double precision function msqglu2running0(m,s,u,alphass)
      implicit none
      double precision s,m,u,msqglu2,alphass
      msqglu2running0=msqglu2(m,s,u)*alphass**2
      end

      double precision function msqglu2running1(m,s,u)
      implicit none
      double precision m,s,u,msqglu2running0,alphas
      msqglu2running1=msqglu2running0(m,s,u,alphas(m**2-s))
      end

      double precision function msqglu2running1b(t)
      use PBGtomsqglurunning1b
      implicit none
      double precision t,u,msqglu2running1
      u=2*m**2-s-t
      msqglu2running1b=msqglu2running1(m,s,u)
      end

      double precision function msqglu3running0(m,s,u,alphasu)
      implicit none
      double precision s,m,u,msqglu3,alphasu
      msqglu3running0=msqglu3(m,s,u)*alphasu**2
      end

      double precision function msqglu3running1(m,s,t,u)
      implicit none
      double precision m,s,t,u,msqglu3running0,alphas
      msqglu3running1=msqglu3running0(m,s,u,alphas(s-m**2+t))
      end

      double precision function msqglu3running1b(t)
      use PBGtomsqglurunning1b
      implicit none
      double precision t,u,msqglu3running1
      u=2*m**2-s-t
      msqglu3running1b=msqglu3running1(m,s,t,u)
      end

      double precision function msqglu4running0(m,s,t,u,alphass,
     &                                alphasu)
      implicit none
      double precision s,m,t,u,msqglu4,alphass,alphasu
      msqglu4running0=msqglu4(m,s,t,u)*alphasu*alphass
      end

      double precision function msqglu4running1(m,s,t,u)
      implicit none
      double precision m,m2,s,t,u,msqglu4running0,alphas
      m2=m**2
      msqglu4running1=msqglu4running0(m,s,t,u,alphas(m2-s),
     &         alphas(m2-u))
      end

      double precision function msqglu4running1b(t)
      use PBGtomsqglurunning1b
      implicit none
      double precision t,u,msqglu4running1
      u=2*m**2-s-t
      msqglu4running1b=msqglu4running1(m,s,t,u)
      end

      double precision function msqglu5running0(m,s,t,u,mu,alphass
     &                             ,alphast)
      implicit none
      double precision s,m,t,u,mu,msqglu5,alphass,alphast
      msqglu5running0=msqglu5(m,s,t,u,mu)*alphast*alphass
      end

      double precision function msqglu5running1(m,s,t,u,mu)
      implicit none
      double precision m,smin,s,t,u,mu,msqglu5running0,alphas
      smin=s-m**2
      msqglu5running1=msqglu5running0(m,s,t,u,mu,alphas(-smin),
     &                     alphas(-t))
      end

      double precision function msqglu6running0(m,s,t,u,mu,alphast,
     &                   alphasu)
      implicit none
      double precision s,m,t,u,mu,msqglu6,alphasu,alphast
      msqglu6running0=msqglu6(m,s,t,u,mu)*alphasu*alphast
      end

      double precision function msqglu6running1(m,s,t,u,mu)
      implicit none
      double precision m,m2,s,t,u,mu,msqglu6running0,alphas
      m2=m**2
      msqglu6running1=msqglu6running0(m,s,t,u,mu,alphas(-t),
     &                    alphas(m2-u))
      end

      double precision function msqgluallrunning1(m,s,t,mu,msqtab)
      implicit none
      integer i
      double precision m,m2,s,smin,t,mu,mumin,u,alphas,alphass,alphast,
     & alphasu,
     & msqglu1running0,msqglu2running0,msqglu3running0,msqglu4running0,
     & msqglu5running0,msqglu6running0,msqtab(6)
      m2=m**2
      smin=s-m2
      mumin=smin+t
      u=m2-mumin
      alphast=alphas(-t)
      alphass=alphas(-smin)
      alphasu=alphas(mumin)
      msqtab(1)=msqglu1running0(m,s,t,u,mu,alphast)
      msqtab(2)=msqglu2running0(m,s,u,alphass)
      msqtab(3)=msqglu3running0(m,s,u,alphasu)
      msqtab(4)=msqglu4running0(m,s,t,u,alphass,alphasu)
      msqtab(5)=msqglu5running0(m,s,t,u,mu,alphass,alphast)
      msqtab(6)=msqglu6running0(m,s,t,u,mu,alphast,alphasu)
      msqgluallrunning1=0
      do i=1,6
         msqgluallrunning1=msqgluallrunning1+msqtab(i)
      enddo
      return
      end

      double precision function msqgluallrunning1b(t)
      use PBGtomsqglurunning1b
      implicit none
      double precision t,msqtab(6),msqgluallrunning1
      msqgluallrunning1b=msqgluallrunning1(m,s,t,mu,msqtab)
      end

      double precision function msqgluallrunning2(m,s,v,mu)
      implicit none
      double precision m,s,v,t,mu,msqtab(6),msqgluallrunning1
      t=(s-m**2)**2*v/(2*s)
      msqgluallrunning2=msqgluallrunning1(m,s,t,mu,msqtab)
      end

      double precision function msqgluallrunning2b(v)
      use PBGtomsqglurunning2b
      implicit none
      double precision v,msqgluallrunning2
      msqgluallrunning2b=msqgluallrunning2(m,s,v,mu)
      end

      double precision function msqglu6posrunning0(m,s,t,u,mu,alphast,
     &          alphasu)
      implicit none
      double precision s,m,t,u,mu,msqglu6pos,alphasu,alphast
      msqglu6posrunning0=msqglu6pos(m,s,t,u,mu)*alphasu*alphast
      end

      double precision function msqglu6posrunning1(m,s,t,u,mu)
      implicit none
      double precision m,m2,s,t,u,mu,msqglu6posrunning0,alphas
      m2=m**2
      msqglu6posrunning1=msqglu6posrunning0(m,s,t,u,mu,alphas(-t),
     &            alphas(m2-u))
      end

      double precision function msqglu6posrunning1b(t)
      use PBGtomsqglurunning1b
      implicit none
      double precision t,u,msqglu6posrunning1
      u=2*m**2-s-t
      msqglu6posrunning1b=msqglu6posrunning1(m,s,t,u,mu)
      end

      double precision function msqglutestrunning1(m,s,t,mu)
      implicit none
      double precision m,m2,s,smin,t,mu,mumin,u,alphas,alphass,alphast,
     & alphasu,
     & msqglu1running0,msqglu2running0,msqglu3running0,msqglu4running0,
     & msqglu6posrunning0
      m2=m**2
      smin=s-m2
      mumin=smin+t
      u=m2-mumin
      alphast=alphas(-t)
      alphass=alphas(-smin)
      alphasu=alphas(mumin)
      msqglutestrunning1=msqglu1running0(m,s,t,u,mu,alphast)+
     & msqglu2running0(m,s,u,alphass)+
     & msqglu3running0(m,s,u,alphasu)+
     & msqglu4running0(m,s,t,u,alphass,alphasu)+
     & msqglu6posrunning0(m,s,t,u,mu,alphast,alphasu)
      return
      end

      double precision function rejglurunning(m,s,t,mu)
      implicit none
      double precision m,m2,s,smin,t,mu,mumin,u,alphas,alphass,alphast,
     & alphasu,
     & msqglu1running0,msqglu2running0,msqglu3running0,msqglu4running0,
     & msqglu5running0,msqglu6running0,msqglu6posrunning0,sum0,rej
      integer ifmtx
      call getMonitorFileIndex(ifmtx)
      m2=m**2
      smin=s-m2
      mumin=smin+t
      u=m2-mumin
      alphast=alphas(-t)
      alphass=alphas(-smin)
      alphasu=alphas(mumin)
      sum0=msqglu1running0(m,s,t,u,mu,alphast)+
     &     msqglu2running0(m,s,u,alphass)+
     &     msqglu3running0(m,s,u,alphasu)+
     &     msqglu4running0(m,s,t,u,alphass,alphasu)
      rej=(sum0+msqglu5running0(m,s,t,u,mu,alphass,alphast)+
     &          msqglu6running0(m,s,t,u,mu,alphast,alphasu))/
     &    (sum0+msqglu6posrunning0(m,s,t,u,mu,alphast,alphasu))
      if(rej-1.gt.1d-6) then
         write(ifmtx,*) 'rej > 1 in rejglurunning:',rej
         write(ifmtx,*) m,s,t,mu
         close(ifmtx)
         stop
      endif
      rejglurunning=rej
      return
      end

      double precision function momferm0base(m,s,mu)
      implicit none
      integer nf
      double precision s,m,m2,smin,mu,mu2,pi2
      parameter (nf=3,pi2=9.869604401089359d0)
      m2=m**2
      mu2=mu**2
      smin=s-m2
      momferm0base=512*nf*pi2*(2 + 
     &   2*s*(mu2**2+2*smin**2+2*mu2*s)/(mu2*smin**2+mu2**2*s)- 
     &   4*s*(mu2+s)/smin**2*log(1+smin**2/(s*mu2)))
      return
      end
           
      double precision function momferm0approx(m,s,mu)
      implicit none
      integer nf
      double precision s,m,m2,smin,mu,mu2,pi2
      parameter (nf=3,pi2=9.869604401089359d0)
      m2=m**2
      mu2=mu**2
      smin=s-m2
      momferm0approx=1024*nf*pi2*(smin**2/mu2**2-
     &     (2*m2-mu2)*smin**4/(3.*m2**2*mu2**3))
      return
      end

      double precision function momferm0(m,s,mu)
      implicit none
      double precision s,m,mu,momferm0approx,momferm0base
      if(s.le.(m**2+mu*m/10.)) then
        momferm0=momferm0approx(m,s,mu)
      else
        momferm0=momferm0base(m,s,mu)
      endif
      return
      end

      double precision function momferm0asympt(s,mu)
      implicit none
      integer nf
      double precision s,mu,pi2
      parameter (nf=3,pi2=9.869604401089359d0)
      momferm0asympt=2048*nf*pi2*s/mu**2
      return
      end
        
      double precision function ratiomomferm0(m,s,mu)
      implicit none
      double precision m,s,mu,momferm0asympt,momferm0
      ratiomomferm0= momferm0(m,s,mu)/momferm0asympt(s,mu)
      return
      end                

      double precision function momglu0part1(m,s,mu)
      implicit none
      double precision s,m,m2,smin,mu,mu2,pi2
      parameter (pi2=9.869604401089359d0)
      m2=m**2
      mu2=mu**2
      smin=s-m2
      if(smin.eq.0) then
         momglu0part1=0
      else
         momglu0part1=pi2*96*64*s*((1/mu2-m2/(smin**2+mu2*s))- 
     &        log(1+smin**2/(s*mu2))/smin) 
      endif
      return
      end

      double precision function momglu0part1approx(m,s,mu)
      implicit none
      double precision m,s,mu,m2,smin,mu2,pi2
      parameter (pi2=9.869604401089359d0)
      m2=m**2
      mu2=mu**2
      smin=s-m2
      momglu0part1approx=pi2*96*32*(smin/(m*mu2))**2*(2*m2-smin)
      return
      end

      double precision function momglu0part2(m,s)
      implicit none
      double precision s,m,m2,pi2
      parameter (pi2=9.869604401089359d0)
      integer ifmtx
      call getMonitorFileIndex(ifmtx)
      m2=m**2
      if(s.eq.m2) then
         write(ifmtx,*) 'error in momglu0part2'
         close(ifmtx)
         stop
      else
         momglu0part2=pi2*96*64*(s+m2)**3/(9.*(s-m2)**2*s)
      endif
      return
      end

      double precision function momglu0part2approx(m,s)
      implicit none
      double precision m,s,m2,pi2
      parameter (pi2=9.869604401089359d0)
      m2=m**2
      momglu0part2approx=pi2*96*256*m2*(s+m2)/(9.*(s-m2)**2)
      return
      end

      double precision function momglu0part3(m,s)
      implicit none
      double precision s,m,m2,pi2
      parameter (pi2=9.869604401089359d0)
      integer ifmtx
      call getMonitorFileIndex(ifmtx)
      m2=m**2
      if(s.eq.m2) then
         write(ifmtx,*) 'error in momglu0part3'
         close(ifmtx)
         stop
      else
         momglu0part3=pi2*96*128*s*(4*m2+(s-3*m2)*log(s/m2))/
     &        (9.*(s-m2)**2)
      endif
      return
      end

      double precision function momglu0part3approx(m,s)
      implicit none
      double precision m,s,m2,pi2
      parameter (pi2=9.869604401089359d0)
      m2=m**2
      momglu0part3approx=pi2*96*256*m2*(s+m2)/(9.*(s-m2)**2)
      return
      end

      double precision function momglu0part4(m,s)
      implicit none
      double precision s,m,m2,smin,pi2
      parameter (pi2=9.869604401089359d0)
      integer ifmtx
      call getMonitorFileIndex(ifmtx)
      m2=m**2
      smin=s-m2
      if(smin.eq.0) then
         write(ifmtx,*) 'error in momglu0part4'
         close(ifmtx)
         stop
      else
         momglu0part4=pi2*96*32*m2/(9.*smin)*(s*(s+3*m2)/smin**2*
     &        log(s/m2)-1)
      endif
      return
      end

      double precision function momglu0part4approx(m,s)
      implicit none
      double precision m,s,m2,pi2
      parameter (pi2=9.869604401089359d0)
      m2=m**2
      momglu0part4approx=pi2*96*64*m2*(s+m2)/(9.*(s-m2)**2)
      return
      end

      double precision function momglu0part5(m,s,mu)
      implicit none
      double precision s,m,m2,mu,mu2,smin,pi2
      parameter (pi2=9.869604401089359d0)
      integer ifmtx
      call getMonitorFileIndex(ifmtx)
      m2=m**2
      smin=s-m2
      mu2=mu**2
      if(smin.eq.0) then
         write(ifmtx,*) 'error in momglu0part5'
         close(ifmtx)
         stop
      else
         momglu0part5=-pi2*96*32*s/smin*
     &        ((s*(s+mu2)-m2**2)*log(1+smin**2/(s*mu2))/smin**2-1)
      endif
      return
      end

      double precision function momglu0part6(m,s,mu)
      implicit none
      double precision s,m,m2,mu,mu2,smin,pi2
      parameter (pi2=9.869604401089359d0)
      integer ifmtx
      call getMonitorFileIndex(ifmtx)
      m2=m**2
      smin=s-m2
      mu2=mu**2
      if(smin.eq.0) then
         write(ifmtx,*) 'error in momglupart6'
         close(ifmtx)
         stop
      else
         momglu0part6=pi2*96*32*s/(smin*(smin+mu2))*(
     &   ((smin-m2+mu2)*(smin-m2)-m2**2)*log(1+smin**2/(s*mu2))/smin-
     &     m2*Log(s/m2))
      endif
      return
      end

      double precision function momglu0base(m,s,mu,momglu0tab)
      implicit none
      integer i
      double precision s,m,mu,momglu0part1,momglu0part2,momglu0part3,
     & momglu0part4,momglu0part5,momglu0part6,momglu0tab(6)
      momglu0tab(1)=momglu0part1(m,s,mu)
      momglu0tab(2)=momglu0part2(m,s)
      momglu0tab(3)=momglu0part3(m,s)
      momglu0tab(4)=momglu0part4(m,s)
      momglu0tab(5)=momglu0part5(m,s,mu)
      momglu0tab(6)=momglu0part6(m,s,mu)
      momglu0base=0
      do i=1,6
         momglu0base=momglu0base+momglu0tab(i)
      enddo
      return
      end
           
      double precision function momglu0approx(m,s)
      implicit none
      double precision s,m,m2,smin,pi2
      parameter (pi2=9.869604401089359d0)
      m2=m**2
      smin=s-m2
      momglu0approx=6144*pi2*m2/smin*(1+2*m2/smin)
      return
      end

      double precision function momglu0(m,s,mu,ifsmalls,momglu0tab)
      implicit none
      logical ifsmalls
      double precision s,m,mu,momglu0approx,momglu0base,momglu0tab(6)
      if(s.le.(m**2+mu*m/100.)) then
         momglu0=momglu0approx(m,s)
         ifsmalls=.true.
      else
         momglu0=momglu0base(m,s,mu,momglu0tab)
         ifsmalls=.false.
      endif
      return
      end

      double precision function momglu0approx3(m,s,mu)
      implicit none
      double precision m,s,mu,m2,pi2
      parameter (pi2=9.869604401089359d0)
      m2=m**2
      momglu0approx3=pi2*6144*s*(2*m2/(s-m2)**2+1/mu**2)
      return
      end

      double precision function ratiomomglu0(m,s,mu,ifsmalls,momglu0tab)
      implicit none
      logical ifsmalls
      double precision m,s,mu,momglu0,momglu0approx3,
     & momglu0tab(6)
      ratiomomglu0=momglu0(m,s,mu,ifsmalls,momglu0tab)/
     &  momglu0approx3(m,s,mu)
      return
      end

      double precision function momglu0part6pos(m,s,mu)
      implicit none
      double precision s,m,m2,mu,mu2,smin,pi2
      parameter (pi2=9.869604401089359d0)
      integer ifmtx
      call getMonitorFileIndex(ifmtx)
      m2=m**2
      smin=s-m2
      mu2=mu**2
      if(smin.eq.0) then
         write(ifmtx,*) 'error in momglupart6pos'
         close(ifmtx)
         stop
      else
         momglu0part6pos=pi2*96*32*s/smin*log(1+smin**2/(s*mu2))
      endif
      return
      end

      double precision function momglu0part6posapprox(m,s,mu)
      implicit none
      double precision s,m,mu,mu2,smin,pi2
      parameter (pi2=9.869604401089359d0)
      integer ifmtx
      call getMonitorFileIndex(ifmtx)
      smin=s-m**2
      mu2=mu**2
      if(smin.eq.0) then
         write(ifmtx,*) 'error in momglupart6posapprox'
         close(ifmtx)
         stop
      else
         momglu0part6posapprox=pi2*96*32*smin/mu2*(1-smin**2/(2*s*mu2))
      endif
      return
      end

      double precision function momferm0running1(m0,s0,mu0)
      use PBGtomsqfermrunning1b
      implicit none
      double precision m0,s0,mu0,tinf,tsup,msqfermrunning1b,eps,
     & ss1,ss2,ss
      parameter(tsup=-1D-16,eps=1D-8)
      external msqfermrunning1b
      m=m0
      s=s0
      mu=mu0
      tinf=(s-m**2)**2/s
      if(abs(tinf).gt.4) then
         call qromb(msqfermrunning1b,-2.d0,tsup,ss1,eps)
         call qromb(msqfermrunning1b,-tinf,-2.d0,ss2,eps)
         ss=ss1+ss2
      else
         call qromb(msqfermrunning1b,-tinf,tsup,ss,eps)
      endif
      momferm0running1=ss*2/tinf
      return
      end

      double precision function momferm0running2(m0,s0,mu0)
      use PBGtomsqfermrunning2b
      implicit none
      double precision m0,s0,mu0,vinf,vsup,msqfermrunning2b,ss,
     &  eps
      parameter(vinf=-2,vsup=-1D-16,eps=1D-7)
      external msqfermrunning2b
      m=m0
      s=s0
      mu=mu0
      call qromb(msqfermrunning2b,vinf,vsup,ss,eps)
      momferm0running2=ss
      return
      end

      double precision function momglu0part1running1(m0,s0,mu0)
      use PBGtomsqglurunning1b
      implicit none
      double precision m0,s0,mu0,tinf,tsup,msqglu1running1b,eps,
     & ss1,ss2,ss3,ss
      parameter(tsup=-1D-16,eps=1D-8)
      external msqglu1running1b
      m=m0
      s=s0
      mu=mu0
      tinf=(s-m**2)**2/s
      if(abs(tinf).gt.4) then
         call qromb(msqglu1running1b,-2.d0,tsup,ss1,eps)
         call qromb(msqglu1running1b,-tinf,-tinf+2.d0,ss2,eps)
         call qromb(msqglu1running1b,-tinf+2.d0,-2.d0,ss3,eps)
         ss=ss1+ss2+ss3
      else
         call qromb(msqglu1running1b,-tinf,tsup,ss,eps)
      endif
      momglu0part1running1=ss*2/tinf
      return
      end

      double precision function momglu0part2running1(m0,s0,mu0)
      use PBGtomsqglurunning1b
      implicit none
      double precision m0,s0,mu0,tinf,tsup,msqglu2running1b,eps,
     & ss1,ss2,ss3,ss
      parameter(tsup=-1D-16,eps=1D-8)
      external msqglu2running1b
      m=m0
      s=s0
      mu=mu0
      tinf=(s-m**2)**2/s
      if(abs(tinf).gt.4) then
         call qromb(msqglu2running1b,-2.d0,tsup,ss1,eps)
         call qromb(msqglu2running1b,-tinf,-tinf+2.d0,ss2,eps)
         call qromb(msqglu2running1b,-tinf+2.d0,-2.d0,ss3,eps)
         ss=ss1+ss2+ss3
      else
         call qromb(msqglu2running1b,-tinf,tsup,ss,eps)
      endif
      momglu0part2running1=ss*2/tinf
      return
      end

      double precision function momglu0part3running1(m0,s0,mu0)
      use PBGtomsqglurunning1b
      implicit none
      double precision m0,s0,mu0,tinf,tsup,msqglu3running1b,eps,
     & ss1,ss2,ss3,ss
      parameter(tsup=-1D-16,eps=1D-8)
      external msqglu3running1b
      m=m0
      s=s0
      mu=mu0
      tinf=(s-m**2)**2/s
      if(abs(tinf).gt.4) then
         call qromb(msqglu3running1b,-2.d0,tsup,ss1,eps)
         call qromb(msqglu3running1b,-tinf,-tinf+2.d0,ss2,eps)
         call qromb(msqglu3running1b,-tinf+2.d0,-2.d0,ss3,eps)
         ss=ss1+ss2+ss3
      else
         call qromb(msqglu3running1b,-tinf,tsup,ss,eps)
      endif
      momglu0part3running1=ss*2/tinf
      return
      end

      double precision function momglu0part4running1(m0,s0,mu0)
      use PBGtomsqglurunning1b
      implicit none
      double precision m0,s0,mu0,tinf,tsup,msqglu4running1b,eps,
     & ss1,ss2,ss3,ss
      parameter(tsup=-1D-16,eps=1D-8)
      external msqglu4running1b
      m=m0
      s=s0
      mu=mu0
      tinf=(s-m**2)**2/s
      if(abs(tinf).gt.4) then
         call qromb(msqglu4running1b,-2.d0,tsup,ss1,eps)
         call qromb(msqglu4running1b,-tinf,-tinf+2.d0,ss2,eps)
         call qromb(msqglu4running1b,-tinf+2.d0,-2.d0,ss3,eps)
         ss=ss1+ss2+ss3
      else
         call qromb(msqglu4running1b,-tinf,tsup,ss,eps)
      endif
      momglu0part4running1=ss*2/tinf
      return
      end

      double precision function momglu0part6posrunning1(m0,s0,mu0)
      use PBGtomsqglurunning1b
      implicit none
      double precision m0,s0,mu0,tinf,tsup,eps,
     & ss1,ss2,ss3,ss
      parameter(tsup=-1D-16,eps=1D-8)
      external msqglu6posrunning1b
      m=m0
      s=s0
      mu=mu0
      tinf=(s-m**2)**2/s
      if(abs(tinf).gt.4) then
         call qromb(msqglu6posrunning1b,-2.d0,tsup,ss1,eps)
         call qromb(msqglu6posrunning1b,-tinf,-tinf+2.d0,ss2,eps)
         call qromb(msqglu6posrunning1b,-tinf+2.d0,-2.d0,ss3,eps)
         ss=ss1+ss2+ss3
      else
         call qromb(msqglu6posrunning1b,-tinf,tsup,ss,eps)
      endif
      momglu0part6posrunning1=ss*2/tinf
      return
      end

      double precision function momglu0partrunning1(ipart,m,s,mu)
      implicit none
      integer ipart
      double precision m,s,mu,momglu0part1running1,momglu0part2running1,
     &  momglu0part3running1,momglu0part4running1,
     &  momglu0part6posrunning1
      integer ifmtx
      call getMonitorFileIndex(ifmtx)
      if(ipart.eq.1) then
         momglu0partrunning1=momglu0part1running1(m,s,mu)
      else if(ipart.eq.2) then
         momglu0partrunning1=momglu0part2running1(m,s,mu)
      else if(ipart.eq.3) then
         momglu0partrunning1=momglu0part3running1(m,s,mu)
      else if(ipart.eq.4) then
         momglu0partrunning1=momglu0part4running1(m,s,mu)
      else if(ipart.eq.5) then
         momglu0partrunning1=momglu0part6posrunning1(m,s,mu)
      else
         write(ifmtx,*) 'unrecognized ipart in momglu0partrunning1'
         close(ifmtx)
         stop
      endif
      end

      double precision function momglu0running1(m0,s0,mu0)
      use PBGtomsqglurunning1b
      implicit none
      double precision m0,s0,mu0,tinf,tsup,msqgluallrunning1b,
     & eps,
     & ss1,ss2,ss3,ss
      parameter(tsup=-1D-16,eps=1D-8)
      external msqgluallrunning1b
      m=m0
      s=s0
      mu=mu0
      tinf=(s-m**2)**2/s
      if(abs(tinf).gt.4) then
         call qromb(msqgluallrunning1b,-2.d0,tsup,ss1,eps)
         call qromb(msqgluallrunning1b,-tinf,-tinf+2.d0,ss2,eps)
         call qromb(msqgluallrunning1b,-tinf+2.d0,-2.d0,ss3,eps)
         ss=ss1+ss2+ss3
      else
         call qromb(msqgluallrunning1b,-tinf,tsup,ss,eps)
      endif
      momglu0running1=ss*2/tinf
      return
      end

      double precision function momglu0running2(m0,s0,mu0)
      use PBGtomsqglurunning2b
      implicit none
      double precision m0,s0,mu0,vinf,vsup,ss,eps
      parameter(vinf=-2,vsup=-1D-16,eps=1D-7)
      external msqgluallrunning2b
      m=m0
      s=s0
      mu=mu0
      call qromb(msqgluallrunning2b,vinf,vsup,ss,eps)
      momglu0running2=ss
      return
      end


      double precision function ratiomomferm0running(m,s,mu)
      implicit none
      double precision m,s,mu,momferm0asympt,momferm0running1,alphasmax
      ratiomomferm0running=momferm0running1(m,s,mu)/
     &         (alphasmax()**2*momferm0asympt(s,mu))
      return
      end                

      double precision function ratiomomglu0running(m,s,mu)
      implicit none
      double precision m,s,mu,momglu0approx3,momglu0running1,alphasmax
      ratiomomglu0running=momglu0running1(m,s,mu)/
     &         (alphasmax()**2*momglu0approx3(m,s,mu))
      return
      end                

      subroutine storeratioferm
      use PBGpsiinfo2
      use PBGforratioprbferm
      implicit none
      integer imu,ilns,iflavor
      double precision mu,mass(2),s,ratiomomferm0running
      character*40 dataprobnamec(2),dataprobnameb,dataprobname(2)
      data dataprobnamec/'../data.dir/ratioprobfermc_1500.dat',
     & '../data.dir/ratioprobfermc_1940.dat'/
      data dataprobnameb/'../data.dir/ratioprobfermb.dat'/
      integer ifmtx
      call getMonitorFileIndex(ifmtx)
      mumin=0.05d0
      mumax=2.d0
      deltmu=0.05d0
      nbmu=nint((mumax-mumin)/deltmu)
      if(nbmu.gt.nbmumax) then
         write(ifmtx,*) 'nbmu too large in storeratioferm'
         close(ifmtx)
         stop
      endif
      lnsmin=-3.d0
      lnsmax=3.d0
      deltlns=0.05d0
      nblns=nint((lnsmax-lnsmin)/deltlns)
      if(mcchoice.eq.1) then
         mass(1)=1.5d0
      else
         mass(1)=1.94d0
      endif
      mass(2)=5.1d0
      dataprobname(1)=dataprobnamec(mcchoice)
      dataprobname(2)=dataprobnameb
      do iflavor=1,2
         open(10,file=dataprobname(iflavor),status='new')
         write(10,*) mumin
         write(10,*) mumax
         write(10,*) nbmu
         write(10,*) lnsmin
         write(10,*) lnsmax
         write(10,*) nblns
         do imu=0,nbmu
            mu=mumin+imu*deltmu
            write(6,*) mu
            write(10,*) mu 
            do ilns=0,nblns
               s=mass(iflavor)**2*(1.d0+10.d0**(lnsmin+ilns*deltlns))
               ratioprobferm(iflavor,ilns,imu)=ratiomomferm0running(
     &           mass(iflavor),s,mu)
               write(10,*) ratioprobferm(iflavor,ilns,imu)
            enddo
         enddo
         close(10)
      enddo
      write(6,*) 'storage performed with apparent succes for fermions' 
      return
      end

      subroutine storeratioglu
      use PBGpsiinfo2
      use PBGforratioprbglu
      implicit none
      integer imu,ilns,iflavor,ib
      double precision mu,m,mass(2),s,ratiomomglu0running,
     & momglu0partrunning1,branch(5)
      character*40 dataprobnamec(2),dataprobnameb,dataprobname(2)
      data dataprobnamec/'../data.dir/ratioprobgluc_1500.dat',
     & '../data.dir/ratioprobgluc_1940.dat'/
      data dataprobnameb/'../data.dir/ratioprobglub.dat'/
      integer ifmtx
      call getMonitorFileIndex(ifmtx)
      mumin=0.05d0
      mumax=2.d0
      deltmu=0.05d0
      nbmu=nint((mumax-mumin)/deltmu)
      if(nbmu.gt.nbmumax) then
         write(ifmtx,*) 'nbmu too large in storeratioglu'
         close(ifmtx)
         stop
      endif
      lnsmin=-5.d0
      lnsmax=3.d0
      deltlns=0.02d0
      nblns=nint((lnsmax-lnsmin)/deltlns)
      if(nblns.gt.nblnsmax) then
         write(ifmtx,*) 'nblns too large in storeratioglu'
         close(ifmtx)
         stop
      endif
      if(mcchoice.eq.1) then
         mass(1)=1.5d0
      else
         mass(1)=1.94d0
      endif
      mass(2)=5.1d0
      dataprobname(1)=dataprobnamec(mcchoice)
      dataprobname(2)=dataprobnameb
      do iflavor=1,2
         m=mass(iflavor)
         open(10,file=dataprobname(iflavor),status='new')
         write(10,*) mumin
         write(10,*) mumax
         write(10,*) nbmu
         write(10,*) lnsmin
         write(10,*) lnsmax
         write(10,*) nblns
         do imu=0,nbmu
            mu=mumin+imu*deltmu
            write(6,*) mu
            write(10,*) mu 
            do ilns=0,nblns
               s=mass(iflavor)**2*(1.d0+10.d0**(lnsmin+ilns*deltlns))
            ratioprobglu(iflavor,ilns,imu,1)=ratiomomglu0running(m,s,mu)
               branch(1)=momglu0partrunning1(1,m,s,mu)
               do ib=2,5
                  branch(ib)=branch(ib-1)+momglu0partrunning1(ib,m,s,mu)
               enddo
               if(branch(5).eq.0) then 
                  write(ifmtx,*) 'tot branching=0;store ratio gluons'
                  write(ifmtx,*) iflavor,mu,s
                  close(ifmtx)
                  stop
               endif
               do ib=1,4
                ratioprobglu(iflavor,ilns,imu,ib+1)=branch(ib)/branch(5)
               enddo
               write(10,*) (ratioprobglu(iflavor,ilns,imu,ib),ib=1,5)
            enddo
         enddo
         close(10)
      enddo
      write(6,*) 'storage performed with apparent succes for gluons' 
      return
      end

      subroutine loadratioferm
      use PBGpsiinfo2
      use PBGforratioprbferm
      implicit none
#include "aaa.h"
      integer imu,ilns,lentrue,iflavor
      double precision mu,mutry
      character*40 dataprobnamec(2),dataprobnameb,dataprobname(2)
      data dataprobnamec/'RATES/ratioprobfermc_1500.dat',
     & 'RATES/ratioprobfermc_1940.dat'/
      data dataprobnameb/'RATES/ratioprobfermb.dat'/
      integer ifmtx
      call getMonitorFileIndex(ifmtx)

      dataprobname(1)=dataprobnamec(mcchoice)
      dataprobname(2)=dataprobnameb
      do iflavor=1,2
         open(152,file=fnus1(1:index(fnus1,' ')-1)
     &  //dataprobname(iflavor)(1:lentrue(dataprobname(iflavor),40))
     &  ,status='old')
         read(152,*) mumin
         read(152,*) mumax
         read(152,*) nbmu
         if(nbmu.gt.nbmumax) then
            write(ifmtx,*) 'nbmu too large in loadratioferm'
            close(ifmtx)
            stop
         endif
         deltmu=(mumax-mumin)/nbmu
         read(152,*) lnsmin
         read(152,*) lnsmax
         read(152,*) nblns
         if(nblns.gt.nblnsmax) then
            write(ifmtx,*) 'nblns too large in loadratioferm'
            close(ifmtx)
            stop
         endif
         deltlns=(lnsmax-lnsmin)/nblns
         do imu=0,nbmu
            mutry=mumin+imu*deltmu
            read(152,*) mu 
            if(abs(mu-mutry).gt.1d-6) then
               write(ifmtx,*) 'disagree in loadratioferm:',mu,mutry
               close(ifmtx)
               stop
            endif
            do ilns=0,nblns
               read(152,*) ratioprobferm(iflavor,ilns,imu)
            enddo
         enddo
         close(152)
      enddo
      return
      end

      double precision function gimmeratioprobferm(itypq,s,mu)
      use PBGpsiinfo
      use PBGupsinfo
      use PBGforratioprbferm
      implicit none
      integer itypq,imuinf,imusup,ilnsinf,ilnssup
      double precision mu,lns,s,lnsrap,murap,dlns,dmu,
     & aii,ais,asi,ass,log10
      parameter(log10=2.302585092994046)
      integer ifmtx
      call getMonitorFileIndex(ifmtx)
      gimmeratioprobferm=0.d0
      if((mu.gt.mumax).or.(mu.lt.mumin)) then
         write(6,*) 'mu out of bounds in ratioprobferm:',mu,mumin,mumax
                  if(mu.gt.mumax) then
         mu=mumax*0.999999d0
         else 
         mu=1.0000001d0*mumin 
         endif
      endif
      if(itypq.eq.1) then
         lns=log(s/mc2-1)/log10
      else
         if(itypq.eq.2) then
            lns=log(s/mb2-1)/log10
         else
            write(ifmtx,*) 'unrecognized itypquark in ratioprobferm:', 
     &           itypq
            close(ifmtx)   
            stop
         endif
      endif
      if(lns.lt.lnsmin) lns=lnsmin
      if(lns.gt.lnsmax) lns=lnsmax
      if(mu.eq.mumax) then
         if(lns.eq.lnsmax)then
            gimmeratioprobferm=ratioprobferm(itypq,nblns,nbmu)
         else
            lnsrap=(lns-lnsmin)/deltlns
            ilnsinf=int(lnsrap)
            ilnssup=ilnsinf+1
            dlns=lnsrap-ilnsinf
            gimmeratioprobferm=dlns*ratioprobferm(itypq,ilnssup,nbmu)+
     &          (1-dlns)*ratioprobferm(itypq,ilnsinf,nbmu)
         endif
      else 
         murap=(mu-mumin)/deltmu
         imuinf=int(murap)
         imusup=imuinf+1
         dmu=murap-imuinf
         if(lns.eq.lnsmax) then
            gimmeratioprobferm=dmu*ratioprobferm(itypq,nblns,imusup)+
     &         (1-dmu)*ratioprobferm(itypq,nblns,imuinf)
         else
            lnsrap=(lns-lnsmin)/deltlns
            ilnsinf=int(lnsrap)
            ilnssup=ilnsinf+1
            dlns=lnsrap-ilnsinf
            aii=(1-dlns)*(1-dmu)
            ais=(1-dlns)*dmu
            asi=dlns*(1-dmu)
            ass=dlns*dmu
            gimmeratioprobferm=
     &           aii*ratioprobferm(itypq,ilnsinf,imuinf)+
     &           ais*ratioprobferm(itypq,ilnsinf,imusup)+
     &           asi*ratioprobferm(itypq,ilnssup,imuinf)+
     &           ass*ratioprobferm(itypq,ilnssup,imusup)
         endif
      endif
      return
      end

      double precision function gimmeratioprobferm2(itypq,s,temp,kappa)
      implicit none 
      integer itypq
      double precision s,temp,kappa,mu,mdebeffinterp,gimmeratioprobferm
      mu=sqrt(kappa)*mdebeffinterp(temp)
      gimmeratioprobferm2=gimmeratioprobferm(itypq,s,mu)
      return
      end

      subroutine testratioprobferm
      implicit none
      integer i
      double precision s,mu,gimmeratioprobferm,ran2,rat,
     &  ratiomomferm0running
      do i=1,20
         mu=0.05+1.45*ran2()
         s=2.25d0*(1+10**(-3+ran2()*6))
         rat=gimmeratioprobferm(1,s,mu)
         write(6,*) mu,s,rat,ratiomomferm0running(1.5d0,s,mu)
      enddo
      return
      end

      subroutine loadratioglu
      use PBGpsiinfo2
      use PBGforratioprbglu
      implicit none
#include "aaa.h"
      integer imu,ilns,ib,iflavor,lentrue
      double precision mu,mutry
      character*40 dataprobnamec(2),dataprobnameb,dataprobname(2)
      data dataprobnamec/'RATES/ratioprobgluc_1500.dat',
     & 'RATES//ratioprobgluc_1940.dat'/
      data dataprobnameb/'RATES//ratioprobglub.dat'/
      integer ifmtx
      call getMonitorFileIndex(ifmtx)

      dataprobname(1)=dataprobnamec(mcchoice)
      dataprobname(2)=dataprobnameb
      do iflavor=1,2
         open(152,file=fnus1(1:index(fnus1,' ')-1)
     &  //dataprobname(iflavor)(1:lentrue(dataprobname(iflavor),40))
     &  ,status='old')
         read(152,*) mumin
         read(152,*) mumax
         read(152,*) nbmu
         if(nbmu.gt.nbmumax) then
            write(ifmtx,*) 'nbmu too large in loadratioglu'
            close(ifmtx)
            stop
         endif
         deltmu=(mumax-mumin)/nbmu
         read(152,*) lnsmin
         read(152,*) lnsmax
         read(152,*) nblns
         if(nblns.gt.nblnsmax) then
            write(ifmtx,*) 'nblns too large in loadratioglu'
            close(ifmtx)
            stop
         endif
         deltlns=(lnsmax-lnsmin)/nblns
         do imu=0,nbmu
            mutry=mumin+imu*deltmu
            read(152,*) mu 
            if(abs(mu-mutry).gt.1d-6) then
               write(ifmtx,*) 'disagree in loadratioglu:',mu,mutry
               close(ifmtx)
               stop
            endif
            do ilns=0,nblns
               read(152,*) (ratioprobglu(iflavor,ilns,imu,ib),ib=1,5)
            enddo
         enddo
         close(152)
      enddo
      return
      end

      double precision function gimmeratioprobglu(itypq,s,mu)
      use PBGpsiinfo
      use PBGupsinfo
      use PBGforratioprbglu
      implicit none
      integer itypq,imuinf,imusup,ilnsinf,ilnssup
      double precision mu,lns,s,lnsrap,murap,dlns,dmu,
     & aii,ais,asi,ass,log10
      parameter(log10=2.302585092994046)
      integer ifmtx
      call getMonitorFileIndex(ifmtx)
      gimmeratioprobglu=0.d0
      if((mu.gt.mumax).or.(mu.lt.mumin)) then
         write(6,*) 'mu out of bounds in ratioprobglu:',mu,mumin,mumax
         if(mu.gt.mumax) then
         mu=mumax*0.999999d0
         else 
         mu=1.0000001d0*mumin 
         endif
      endif
      if(itypq.eq.1) then
         lns=log(s/mc2-1)/log10
      else
         if(itypq.eq.2) then
            lns=log(s/mb2-1)/log10
         else
            write(ifmtx,*) 'unrecognized itypquark in ratioprobglu:',
     &           itypq
            close(ifmtx)
            stop
         endif
      endif
      if(lns.lt.lnsmin) lns=lnsmin
      if(lns.gt.lnsmax) lns=lnsmax
      if(mu.eq.mumax) then
         if(lns.eq.lnsmax)then
            gimmeratioprobglu=ratioprobglu(itypq,nblns,nbmu,1)
         else
            lnsrap=(lns-lnsmin)/deltlns
            ilnsinf=int(lnsrap)
            ilnssup=ilnsinf+1
            dlns=lnsrap-ilnsinf
            gimmeratioprobglu=dlns*ratioprobglu(itypq,ilnssup,nbmu,1)+
     &          (1-dlns)*ratioprobglu(itypq,ilnsinf,nbmu,1)
         endif
      else 
         murap=(mu-mumin)/deltmu
         imuinf=int(murap)
         imusup=imuinf+1
         dmu=murap-imuinf
         if(lns.eq.lnsmax) then
            gimmeratioprobglu=dmu*ratioprobglu(itypq,nblns,imusup,1)+
     &         (1-dmu)*ratioprobglu(itypq,nblns,imuinf,1)
         else
            lnsrap=(lns-lnsmin)/deltlns
            ilnsinf=int(lnsrap)
            ilnssup=ilnsinf+1
            dlns=lnsrap-ilnsinf
            aii=(1-dlns)*(1-dmu)
            ais=(1-dlns)*dmu
            asi=dlns*(1-dmu)
            ass=dlns*dmu
            gimmeratioprobglu=
     &           aii*ratioprobglu(itypq,ilnsinf,imuinf,1)+
     &           ais*ratioprobglu(itypq,ilnsinf,imusup,1)+
     &           asi*ratioprobglu(itypq,ilnssup,imuinf,1)+
     &           ass*ratioprobglu(itypq,ilnssup,imusup,1)
         endif
      endif
      return
      end

      double precision function gimmeratioprobglu2(itypq,s,temp,kappa)
      implicit none 
      integer itypq
      double precision s,temp,kappa,mu,mdebeffinterp,gimmeratioprobglu
      mu=sqrt(kappa)*mdebeffinterp(temp)
      gimmeratioprobglu2=gimmeratioprobglu(itypq,s,mu)
      return
      end

      subroutine testratioprobglu
      implicit none
      integer i
      double precision s,mu,gimmeratioprobglu,ran2,rat,
     &  ratiomomglu0running
      write(6,*) 'testratioprobglu'
      do i=1,20
         mu=0.05+1.45*ran2()
         s=2.25d0*(1+10**(-5+ran2()*8))
         rat=gimmeratioprobglu(1,s,mu)
         write(6,*) mu,s,rat,ratiomomglu0running(1.5d0,s,mu)
      enddo
      return
      end

      subroutine gimmebranchglu(itypq,s,mu,branch)
      use PBGpsiinfo
      use PBGupsinfo
      use PBGforratioprbglu
      implicit none
      integer itypq,imuinf,imusup,ilnsinf,ilnssup,ib
      double precision mu,lns,s,lnsrap,murap,dlns,dmu,
     & aii,ais,asi,ass,log10,branch(4)
      parameter(log10=2.302585092994046)
      integer ifmtx
      call getMonitorFileIndex(ifmtx)
      if((mu.gt.mumax).or.(mu.lt.mumin)) then
         write(6,*) 'mu out of bounds in gimmebranchglu:',mu,mumin,mumax
         if(mu.gt.mumax) then
         mu=mumax*0.999999d0
         else 
         mu=1.0000001d0*mumin 
         endif
      endif
      if(itypq.eq.1) then
         lns=log(s/mc2-1)/log10
      else
         if(itypq.eq.2) then
            lns=log(s/mb2-1)/log10
         else
            write(ifmtx,*) 'unrecognized itypquark in gimmebranchglu:',
     &           itypq
            close(ifmtx)
            stop
         endif
      endif
      if(lns.lt.lnsmin) lns=lnsmin
      if(lns.gt.lnsmax) lns=lnsmax
      if(mu.eq.mumax) then
         if(lns.eq.lnsmax)then
            do ib=2,5
               branch(ib-1)=ratioprobglu(itypq,nblns,nbmu,ib)
            enddo
         else
            lnsrap=(lns-lnsmin)/deltlns
            ilnsinf=int(lnsrap)
            ilnssup=ilnsinf+1
            dlns=lnsrap-ilnsinf
            do ib=2,5
               branch(ib-1)=dlns*ratioprobglu(itypq,ilnssup,nbmu,ib)+
     &          (1-dlns)*ratioprobglu(itypq,ilnsinf,nbmu,ib)
            enddo
         endif
      else 
         murap=(mu-mumin)/deltmu
         imuinf=int(murap)
         imusup=imuinf+1
         dmu=murap-imuinf
         if(lns.eq.lnsmax) then
            do ib=2,5
               branch(ib-1)=dmu*ratioprobglu(itypq,nblns,imusup,ib)+
     &         (1-dmu)*ratioprobglu(itypq,nblns,imuinf,ib)
            enddo
         else
            lnsrap=(lns-lnsmin)/deltlns
            ilnsinf=int(lnsrap)
            ilnssup=ilnsinf+1
            dlns=lnsrap-ilnsinf
            aii=(1-dlns)*(1-dmu)
            ais=(1-dlns)*dmu
            asi=dlns*(1-dmu)
            ass=dlns*dmu
            do ib=2,5
               branch(ib-1)=
     &              aii*ratioprobglu(itypq,ilnsinf,imuinf,ib)+
     &              ais*ratioprobglu(itypq,ilnsinf,imusup,ib)+
     &              asi*ratioprobglu(itypq,ilnssup,imuinf,ib)+
     &              ass*ratioprobglu(itypq,ilnssup,imusup,ib)
            enddo
         endif
      endif
      return
      end

      double precision function qrandlargeu(u,u0,temp)
      implicit none
      logical accept
      double precision u,u0,temp,q,ran2
      accept=.false.
      do while(.not.accept)
        q=-(u+u0)*log(ran2()*ran2())
        accept=(ran2().gt.exp(-2*q*u))
        enddo
        qrandlargeu=temp*q
      return
      end

      double precision function qrandsmallu(u,u0,temp)
      implicit none
      logical accept
      double precision u,u0,temp,q,ran2,thresh,alpha
      accept=.false.
      thresh=1/(1+1/(u0+u)**6)
      do while(.not.accept)
         q=-log(ran2()*ran2()*ran2())
         if(ran2().lt.thresh) then
            q=q*(u+u0)
         else
            q=q*(u0-u)
         endif
         alpha=q*u
         if(alpha.eq.0) then
            accept=(ran2().lt.1.d0)
         else
            accept=(ran2().lt.tanh(alpha)/alpha)
         endif
      enddo
      qrandsmallu=temp*q
      return
      end

      double precision function qrandfermtest(u,u0,temp)
        implicit none
        double precision u,u0,temp,qrandlargeu,qrandsmallu
           if(u.lt.0.317) then
          qrandfermtest=qrandsmallu(u,u0,temp)
        else
          qrandfermtest=qrandlargeu(u,u0,temp)
      endif                        
      return
      end        

      double precision function meanqfermtest(u,temp,nsample,mom2)
      implicit none
      integer nsample,i
      double precision qrandfermtest,u,u0,meanq,temp,oneq,mom2
      u0=sqrt(u**2+1)
      meanq=0.d0
      do i=1,nsample
         oneq=qrandfermtest(u,u0,temp)
         meanq=meanq+oneq
         mom2=mom2+oneq**2
      enddo
      meanqfermtest=meanq/nsample
      mom2=mom2/nsample
      return
      end

      subroutine testqfermtest
      implicit none
      double precision meanqfermtest,u,temp,mq,mom2
      integer nsample,i
      nsample=1000000
      temp=0.4d0
      do i=0,10
         u=i
         mq=meanqfermtest(u,temp,nsample,mom2)
         write(6,*) u,mq,mom2
      enddo
      end

      double precision function qrandferm(m,u,temp,mu)
      implicit none
      double precision m,u,u0,temp,mu,q,qrandfermtest,ratiomomferm0,
     & rej,ran2
      logical accept
      integer ifmtx
      call getMonitorFileIndex(ifmtx)
      u0=sqrt(1+u**2)
      accept=.false.
      do while(.not.accept)
        q=qrandfermtest(u,u0,temp)
        rej=ratiomomferm0(m,m*(m+2*q),mu)
        if(rej.gt.1) then
           write(ifmtx,*) 'reject function >1 in qrandferm'
           write(ifmtx,*) m,mu,temp,q
           close(ifmtx)
           stop
        endif
        accept=(ran2().lt.rej)
      enddo
      qrandferm=q
      return
      end        
        
      double precision function meanqferm(m,u,temp,mu,nsample,mom2)
      implicit none
      integer nsample,i
      double precision qrandferm,m,u,temp,mu,oneq,meanq,mom2
      meanq=0.d0
      mom2=0.d0
      do i=1,nsample
         oneq=qrandferm(m,u,temp,mu)
         meanq=meanq+oneq
         mom2=mom2+oneq**2
      enddo
      meanqferm=meanq/nsample
      mom2=mom2/nsample
      return
      end

      subroutine testqferm
      implicit none
      double precision meanqferm,m,u,temp,mu,mq,mom2
      integer nsample,i
      m=1.5
      nsample=1000000
      temp=0.4d0
      mu=1.1d0
      do i=0,30
         u=i*0.1
         mq=meanqferm(m,u,temp,mu,nsample,mom2)
         write(6,*) u,mq,mom2
      enddo
      end

      double precision function cthetarand(u,temp,q)
      implicit none
      double precision u,temp,q,alpha,em2alpha,ran2
      alpha=q*u/temp
      em2alpha=exp(-2*alpha)
      cthetarand=1+log(1-ran2()*(1-em2alpha))/alpha
      return
      end        

      double precision function sampleinvqsmallu(u,u0,temp)
      implicit none
      logical accept
      double precision u,u0,temp,u0plus,u0minus,intplus,intminus,
     & q,alpha,ran2
      u0plus=u0+u
      u0minus=u0-u
      intplus=u/u0plus
      intminus=u/u0minus
      accept=.false.
      do while(.not.accept) 
         if(ran2()*(intplus+intminus).lt.intplus) then
            q=-temp/u0plus*log(ran2())
         else
            q=-temp/u0minus*log(ran2())
         endif
         alpha=q*u/temp
         if(alpha.eq.0) then
            accept=(ran2().lt.1.d0)
         else
            accept=(ran2().lt.tanh(alpha)/alpha)
         endif
      enddo
      sampleinvqsmallu=q
      return
      end

      subroutine gimmeexpintegral(x,expintegral,ifail)
      implicit none
      integer iofx,ifail
      double precision x,eulergamma,dataexpint(39),xinf,xsup,
     & expintegral
      data dataexpint/0.201021,-0.0990729,-0.353281,-0.580223,-0.788823, 
     &               -0.984118,-1.16926,  -1.34637, -1.51693 ,-1.68206, 
     &               -1.84258, -1.99915,  -2.15228, -2.30239, -2.44983, 
     &               -2.59488, -2.73779,  -2.87876, -3.01797, -3.15556, 
     &               -3.29168, -3.42645,  -3.55995, -3.69229, -3.82354, 
     &               -3.95379, -4.08309,  -4.21151, -4.33909, -4.46589, 
     &               -4.59196, -4.71733,  -4.84204, -4.96612, -5.08961,
     &               -5.21254, -5.33493, -5.45681, -5.5782/
      parameter(eulergamma= 0.5772156649)
      integer ifmtx
      call getMonitorFileIndex(ifmtx)
      ifail=0
      if(isnan(x)) then
         ifail=-1
         return
      endif
      if(x.le.0.d0) then
         write(ifmtx,*) 'x<0 in expintegral'
         close(ifmtx)
         stop
      else 
         if(x.le.0.2) then
            expintegral=-log(x)-eulergamma+x-x**2/4
         else 
            if(x.ge.4) then
               expintegral=exp(-x)*(1/(1+x)+1/(1+x)**3)
            else
               iofx=int(10*(x-0.2))
               xinf=0.2+iofx*0.1
               xsup=xinf+0.1
               iofx=iofx+1
               if((iofx.lt.1).or.(iofx.ge.38)) then
                  write(ifmtx,*) 'error iox inexpintegral'
                  write(ifmtx,*) x,iofx
                  close(ifmtx)
                  stop
               endif
               expintegral=exp(dataexpint(iofx)*(xsup-x)/0.1+
     &                         dataexpint(iofx+1)*(x-xinf)/0.1)
            endif
         endif
      endif      
      return
      end

      double precision function invexpint(ei)
      implicit none
      integer iofei
      double precision ei,z0,zprev,z,eulergamma,invz,eiinf,eisup,
     & datainvexpint(100)
      data datainvexpint /3.21051, 2.66785, 2.3599, 2.14638, 1.98395,
     &  1.85348, 1.74487, 1.65213, 1.57142, 1.50013, 1.43642, 1.37892, 
     &  1.32662, 1.27871, 1.23457, 1.1937,  1.15568, 1.12018, 1.08692, 
     &  1.05565,1.02618, 0.998328,0.971948,0.946909,0.923095,0.900407, 
     &  0.878756,0.858063,0.838258,0.819277,0.801063,0.783567,0.76674, 
     &  0.750541,0.734933,0.719879,0.705347,0.691309,0.677737,0.664606,
     &  0.651893,0.639577,0.627638,0.616058,0.604819,0.593905,0.583302,
     &  0.572996,0.562973,0.553222,0.54373, 0.534487,0.525484,0.51671,
     &  0.508156,0.499813,0.491675,0.483733,0.475979,0.468408,0.461012,
     &  0.453785,0.446723,0.439818,0.433065,0.42646, 0.419998,0.413674,
     &  0.407484,0.401424,0.395489,0.389676,0.383981,0.3784, 0.372931, 
     &  0.367569,0.362313,0.357158,0.352102,0.347143,0.342277,0.337502,
     &  0.332816,0.328217,0.323701,0.319267,0.314913,0.310637,0.306437,
     &  0.30231, 0.298256,0.294272,0.290356,0.286508,0.282724,0.279005,
     &  0.275348,0.271752,0.268215,0.264737/
      parameter(eulergamma= 0.5772156649)
      integer ifmtx
      call getMonitorFileIndex(ifmtx)
      if(ei.le.0) then
         write(ifmtx,*) 'error in invexpint'
         write(ifmtx,*) ei
         close(ifmtx)
         stop
      else
         if(ei.le.0.01) then
            zprev=-1.
            z0=-log(ei)
            z=z0
            do while(abs(z-zprev).ge.1e-5)
               zprev=z
               invz=1/(1+zprev)
               z=z0+log(invz*(1+invz**2))
            enddo
         else
            if(ei.ge.1)then
               zprev=-1.
               z0=exp(-(ei+eulergamma))
               z=z0
               do while(abs(z-zprev).ge.1e-5)
                  zprev=z
                  z=z0*exp(zprev*(1-zprev/4.))
               enddo
            else
               iofei=int(100*(ei-0.01))
               eiinf=0.01+iofei*0.01
               eisup=eiinf+0.01
               iofei=iofei+1
               z=datainvexpint(iofei)*(eisup-ei)/0.01+
     &            datainvexpint(iofei+1)*(ei-eiinf)/0.01
            endif
         endif
      endif
      invexpint=z
      return
      end

      subroutine testinvexpint
      implicit none
      double precision invexpint
      write(6,*) invexpint(0.0023d0)
      write(6,*) invexpint(0.023d0)
      write(6,*) invexpint(0.1d0)
      write(6,*) invexpint(0.632d0)
      write(6,*) invexpint(1.03d0)
      write(6,*) invexpint(2.57d0)
      write(6,*) invexpint(4.17d0)
      write(6,*) invexpint(6.17d0)
      return
      end

      double precision function sampleinvqallu(u,u0,temp)
      implicit none
      logical accept
      double precision u,u0,temp,u0plus,u0minus,intplus,intminus,
     & expplus,expminus,eiplus,eiminus,q,alpha,ran2,
     & invexpint,ratio
      integer ifail1
      integer ifmtx
      call getMonitorFileIndex(ifmtx)
      u0plus=u0+u
      u0minus=u0-u
      if(u.eq.0) then
         write(ifmtx,*) 'u=0 in sampleinvqallu'
         close(ifmtx)
         stop
      else
         ratio=u0/u
      endif
      accept=.false.
      if(isnan(ratio)) then
         write(ifmtx,*) 'isnan ratio in sampleinvqallu'
         write(ifmtx,*) u, u0
         close(ifmtx)
         stop
      endif
      call gimmeexpintegral(ratio+1,eiplus,ifail1)
      call gimmeexpintegral(ratio-1,eiminus,ifail1)
      if (ran2()*log(u0plus/u0minus).gt.(eiminus-eiplus))then
         expplus=exp(-u0plus/u0)  
         expminus=exp(-u0minus/u0)  
         intplus=u/u0plus*(1-expplus)
         intminus=u/u0minus*(1-expminus)
         do while(.not.accept)
            if(ran2()*(intplus+intminus).le.intplus) then
               q=-temp/u0plus*log(1-ran2()*(1-expplus))
            else
               q=-temp/u0minus*log(1-ran2()*(1-expminus))
            endif
            alpha=q*u/temp
            accept=(ran2().le.tanh(alpha)/alpha)
         enddo
      else
         do while(.not.accept)
            q=temp/u0minus*invexpint(ran2()*eiminus)
            alpha=q*u/temp
            accept=(ran2().ge.exp(-2*alpha))
         enddo
      endif
      sampleinvqallu=q
      return
      end

      double precision function qrandglutest(u,u0,temp,mu)
      implicit none
      double precision u,u0,temp,mu,intq,intinvq,ran2,
     & qrandfermtest,sampleinvqsmallu,sampleinvqallu
      integer ifmtx
      call getMonitorFileIndex(ifmtx)
      intq=8*u0*temp**3/mu**2
      if(isnan(u0/u)) then
         write(ifmtx,*) 'isnan u in qrandglutest'
         write(ifmtx,*) u, u0
         close(ifmtx)
         stop
      endif
      if(u.eq.0) then
         intinvq=2*temp
      else
         intinvq=2*log(u0+u)*temp/u
      endif
      if(ran2()*(intq+intinvq).le.intq) then
         qrandglutest=qrandfermtest(u,u0,temp)
      else
         if(u.le.4) then
            qrandglutest=sampleinvqsmallu(u,u0,temp)
         else
            qrandglutest=sampleinvqallu(u,u0,temp)
         endif
      endif
      return
      end

      double precision function meanqglutest(u,temp,mu,nsample,mom2)
      implicit none
      integer nsample,i
      double precision q,qrandglutest,u,u0,meanq,temp,mu,mom2
      u0=sqrt(u**2+1)
      meanq=0.d0
      mom2=0.d0
      do i=1,nsample
         q=qrandglutest(u,u0,temp,mu)
         meanq=meanq+q
         mom2=mom2+q**2
      enddo
      meanqglutest=meanq/nsample
      mom2=mom2/nsample
      return
      end

      subroutine testqglutest
      implicit none
      double precision meanqglutest,u,temp,mu,mq,mom2
      integer nsample,i
      temp=0.4d0
      nsample=40000000
      mu=1.1d0
      do i=0,10
         u=i
         mq=meanqglutest(u,temp,mu,nsample,mom2)
         write(6,*) u,mq,mom2
      enddo
      end

      double precision function qrandglu(m,u,temp,mu,ifsmalls,
     &     momglu0tab)     
      implicit none
      logical ifsmalls
      double precision m,u,u0,temp,mu,q,qrandglutest,ratiomomglu0,
     & rej,ran2,momglu0tab(6),momglu0part1approx,momglu0part2approx
      logical accept
      integer ifmtx
      call getMonitorFileIndex(ifmtx)
      u0=sqrt(1+u**2)
      accept=.false.
      do while(.not.accept)
        q=qrandglutest(u,u0,temp,mu)
        rej=ratiomomglu0(m,m*(m+2*q),mu,ifsmalls,momglu0tab)
        if(rej.gt.1) then
           write(ifmtx,*) 'reject function >1 in qrandglu'
           write(ifmtx,*) m,mu,temp,q
           close(ifmtx)
           stop
        endif
        accept=(ran2().lt.rej)
      enddo
      if(ifsmalls) then
         momglu0tab(1)=momglu0part1approx(m,m*(m+2*q),mu)
         momglu0tab(2)=momglu0part2approx(m,m*(m+2*q))
         momglu0tab(3)=momglu0tab(2)
         momglu0tab(4)=momglu0tab(2)/4.
      endif
      qrandglu=q
      return
      end        
        
      double precision function meanqglu(m,u,temp,mu,nsample,mom2)
      implicit none
      integer nsample,i
      logical ifsmalls
      double precision q,qrandglu,m,u,temp,mu,meanq,momglu0tab(6),
     &  mom2
      meanq=0.d0
      mom2=0.d0
      do i=1,nsample
         q=qrandglu(m,u,temp,mu,ifsmalls,momglu0tab)
         meanq=meanq+q
         mom2=mom2+q**2
      enddo
      mom2=mom2/nsample
      meanqglu=meanq/nsample
      return
      end

      subroutine testqglu
      implicit none
      double precision meanqglu,m,u,temp,mu,mq,mom2
      integer nsample,i
      m=1.5d0
      temp=0.4d0
      nsample=10000000
      mu=1.1d0
      do i=5,10
         u=i
         mq=meanqglu(m,u,temp,mu,nsample,mom2)
         write(6,*) 'u=',u,mq,mom2
      enddo
      end

      double precision function qrandfermrunning(itypq,u,temp,mu)
      use PBGpsiinfo
      use PBGupsinfo
      implicit none
      integer itypq
      double precision m,u,u0,temp,mu,q,qrandfermtest,
     & gimmeratioprobferm,
     & rej,ran2
      logical accept
      integer ifmtx
      call getMonitorFileIndex(ifmtx)
      if(itypq.eq.1) then
         m=mc
      else
         if(itypq.eq.2) then
            m=mb
         else
           write(ifmtx,*)'unrecognized itypquark in qrandfermrunning:',
     &           itypq
            close(ifmtx)
            stop
         endif
      endif
      u0=sqrt(1+u**2)
      accept=.false.
      do while(.not.accept)
        q=qrandfermtest(u,u0,temp)
        rej=gimmeratioprobferm(itypq,m*(m+2*q),mu)
        if(rej.gt.1) then
           write(ifmtx,*) 'reject function >1 in qrandfermrunning'
           write(ifmtx,*) itypq,m,mu,temp,q
           close(ifmtx)
           stop
        endif
        accept=(ran2().lt.rej)
      enddo
      qrandfermrunning=q
      return
      end        

      double precision function meanqfermrunning(u,temp,mu,nsample,mom2)
      implicit none
      integer nsample,i
      double precision qrandfermrunning,u,temp,mu,meanq,oneq,mom2
      meanq=0.d0
      mom2=0.d0
      do i=1,nsample
         oneq=qrandfermrunning(1,u,temp,mu)
         meanq=meanq+oneq
         mom2=mom2+oneq**2
      enddo
      meanqfermrunning=meanq/nsample
      mom2=mom2/nsample
      return
      end
      
      subroutine testqfermrunning
      implicit none
      double precision meanqfermrunning,u,temp,mu,mq,mom2
      integer nsample,i
      temp=0.4d0
      nsample=1000000
      mu=1.1d0
      do i=0,10
         u=i
         mq=meanqfermrunning(u,temp,mu,nsample,mom2)
         write(6,*) u,mq,mom2
      enddo
      end
        

      double precision function qrandglurunning(itypq,u,temp,mu)
      use PBGpsiinfo
      use PBGupsinfo
      implicit none
      integer itypq
      double precision m,u,u0,temp,mu,q,qrandglutest,
     & gimmeratioprobglu,
     & rej,ran2
      logical accept
      integer ifmtx
      call getMonitorFileIndex(ifmtx)
      if(itypq.eq.1) then
         m=mc
      else
         if(itypq.eq.2) then
            m=mb
         else
            write(ifmtx,*)'unrecognized itypquark in qrandglurunning:',
     &           itypq
            close(ifmtx)
            stop
         endif
      endif
      u0=sqrt(1+u**2)
      if(isnan(u0/u)) then
         write(ifmtx,*) 'isnan u0/u in qrandglurunning'
         write(ifmtx,*) u, u0
         close(ifmtx)
         stop
      endif
      accept=.false.
      do while(.not.accept)
        q=qrandglutest(u,u0,temp,mu)
        rej=gimmeratioprobglu(itypq,m*(m+2*q),mu)
        if(rej.gt.1) then
           write(ifmtx,*) 'reject function >1 in qrandglurunning'
           write(ifmtx,*) itypq,m,mu,temp,q
           close(ifmtx)
           stop
        endif
        accept=(ran2().lt.rej)
      enddo
      qrandglurunning=q
      return
      end        

      double precision function meanqglurunning(u,temp,mu,nsample,mom2)
      implicit none
      integer nsample,i,itypq
      double precision qrandglurunning,u,temp,mu,meanq,oneq,mom2
      itypq=1
      meanq=0.d0
      mom2=0.d0
      do i=1,nsample
         oneq=qrandglurunning(itypq,u,temp,mu)
         meanq=meanq+oneq
         mom2=mom2+oneq**2
      enddo
      meanqglurunning=meanq/nsample
      mom2=mom2/nsample
      return
      end

      subroutine testqglurunning
      implicit none
      double precision meanqglurunning,u,temp,mu,mq,mom2
      integer nsample,i
      temp=0.4d0
      nsample=10000000
      mu=1.10877d0
      do i=4,10
         u=i
         mq=meanqglurunning(u,temp,mu,nsample,mom2)
         write(6,*) 'u=',u,mq,mom2
      enddo
      end

      double precision function trandferm(m,s,mu)
      implicit none
      logical accept
      double precision m,mu,s,ran2,m2,mu2,smin,iprel2,tsample,r
      m2=m**2
      mu2=mu**2
      smin=s-m2
      iprel2=s/smin**2
      accept=.false.
      do while(.not.accept)
         r=ran2()        
         tsample=-(1-r)/(r/mu2+iprel2)
         accept=(ran2().lt.0.5d0+((smin+tsample)**2/2
     &        +m2*tsample)/smin**2)
      enddo
      trandferm=tsample
      return
      end

      double precision function trandfermforward(m,s,mu)
      implicit none
      logical accept
      double precision m,mu,s,ran2,m2,mu2,smin,iprel2,tsample,r
      m2=m**2
      mu2=mu**2
      smin=(s-m2)
      iprel2=s/(0.25*smin**2)
      accept=.false.
      do while(.not.accept)
         r=ran2()        
         tsample=-(1-r)/(r/mu2+iprel2)
         accept=(ran2().lt.0.5d0+((smin+tsample)**2/2
     &       +m2*tsample)/smin**2)
      enddo
      trandfermforward=tsample
      return
      end

      double precision function trandglu1(m,s,mu)
      implicit none
      logical accept
      double precision m,mu,s,ran2,m2,mu2,smin,iprel2,tsample,r
      m2=m**2
      mu2=mu**2
      smin=s-m2
      iprel2=s/smin**2
      accept=.false.
      do while(.not.accept)
         r=ran2()        
         tsample=-(1-r)/(r/mu2+iprel2)
         accept=(ran2().lt.(1+tsample/smin))
      enddo
      trandglu1=tsample
      return
      end

      double precision function trandglu2(m,s)
      implicit none
      double precision m,s,ran2,m2,smin,splus
      m2=m**2
      smin=s-m2
      splus=s+m2
      if(ran2()*splus**2.lt.4*m2*s) then
         trandglu2=-ran2()*smin**2/s
      else
         trandglu2=(sqrt(ran2()*smin*splus+m2**2)/s-1)*smin
      endif
      return
      end

      double precision function trandglu3(m,s,mu)
      implicit none
      logical accept
      double precision m,mu,s,ran2,m2,mu2,smin,uminsample,rej
      integer ifmtx
      call getMonitorFileIndex(ifmtx)
      m2=m**2
      mu2=mu**2
      smin=s-m2
      accept=.false.
      do while(.not.accept)
         uminsample=-1/(ran2()/m2+1/smin)
         if(smin.gt.2*m2) then
            rej=(-uminsample*(smin-2*m2)+4*m2**2)/
     &          (smin*(smin-2*m2)+4*m2**2)
         else
            rej=(-uminsample*(smin-2*m2)+4*m2**2)*s/
     &          ((s**2+3*m2**2)*m2)
         endif
         if(rej.gt.1) then
            write(ifmtx,*) 'rej > 1 in trandglu3; m:',m,'/mu:',mu,
     &                 '/s:',s,'/u:',uminsample+m2
            close(ifmtx)
            stop
         endif
         accept=(ran2().lt.rej)
      enddo
      trandglu3=-smin-uminsample
      return
      end

      double precision function trandglu4(m,s,mu)
      implicit none
      logical accept
      double precision m,mu,s,ran2,m2,smin,uminsample,rej
      integer ifmtx
      call getMonitorFileIndex(ifmtx)
      m2=m**2
      smin=s-m2
      accept=.false.
      do while(.not.accept)
         uminsample=-smin*(m2/s)**ran2()
         rej=s*(s+3*m2+uminsample)/(s+m2)**2
         if(rej.gt.1) then
            write(ifmtx,*) 'rej > 1 in trandglu4; mu',m,'/mu=',mu,
     &                 '/s=',s,'u=',uminsample+m2
            close(ifmtx)
            stop
         endif
         accept=(ran2().lt.rej)
      enddo
      trandglu4=-smin-uminsample
      return
      end

      double precision function trandglu6pos(m,s,mu)
      implicit none
      double precision m,mu,s,ran2,mu2,smin
      mu2=mu**2
      smin=s-m**2
      trandglu6pos=mu2*(1-(1+smin**2/(mu2*s))**ran2())
      return
      end

      double precision function trandglu(m,s,mu,iftab,momglu0tab,
     &  ifsmalls)
      implicit none
      logical iftab,ifsmalls,accept
      double precision m,mu,s,ran2,r,integ(4),integall,
     & momglu0part1approx,momglu0part2approx,momglu0part1,momglu0part2,
     & momglu0part3,momglu0part4,momglu0part6pos,momglu0part6posapprox,
     & momglu0tab(6),trandglu1,trandglu2,trandglu3,trandglu4,
     & trandglu6pos,t,msqgluall,msqtab(6),msqall,msqtest,msqglu6pos,rej
      integer ifmtx
      call getMonitorFileIndex(ifmtx)
      if(.not.iftab)then
         if(ifsmalls) then
            momglu0tab(1)=momglu0part1approx(m,s,mu)
            momglu0tab(2)=momglu0part2approx(m,s)
            momglu0tab(3)=momglu0tab(2)
            momglu0tab(4)=momglu0tab(2)/4.
         else
            momglu0tab(1)=momglu0part1(m,s,mu)
            momglu0tab(2)=momglu0part2(m,s)
            momglu0tab(3)=momglu0part3(m,s)
            momglu0tab(4)=momglu0part4(m,s)
         endif
      endif            
      integ(1)=momglu0tab(1)
      integ(2)=integ(1)+momglu0tab(2)
      integ(3)=integ(2)+momglu0tab(3)
      integ(4)=integ(3)+momglu0tab(4)
      if(ifsmalls)then
         integall=integ(4)+momglu0part6posapprox(m,s,mu)    
      else
         integall=integ(4)+momglu0part6pos(m,s,mu)
      endif            
      accept=.false.
      do while(.not.accept)
         r=ran2()*integall
         if(r.lt.integ(1)) then
            t=trandglu1(m,s,mu)
         else if(r.lt.integ(2)) then
            t=trandglu2(m,s)
         else if(r.lt.integ(3)) then
            t=trandglu3(m,s,mu)
         else if(r.lt.integ(4)) then
            t=trandglu4(m,s,mu)
         else
            t=trandglu6pos(m,s,mu)
         endif
         msqall=msqgluall(m,s,t,mu,msqtab)
         msqtest=msqglu6pos(m,s,t,2*m**2-s-t,mu)+msqtab(1)+
     &           msqtab(2)+msqtab(3)+msqtab(4)
         rej=msqall/msqtest
         if(rej.gt.1) then
            write(ifmtx,*) 'rej>1 in trandglu:',rej,'/s:',s,'/t:',t
            write(ifmtx,*) msqall,msqtest
            write(ifmtx,*) msqtab
            write(ifmtx,*) trandglu6pos(m,s,mu)
            close(ifmtx)
            stop
         endif
         accept=(ran2().lt.rej)
      enddo
      trandglu=t
      return
      end

      double precision function meantrandglu(m,s,mu)
      implicit none
      logical ifsmalls
      integer nsample,i
      double precision m,s,mu,momglu0tab(6),trandglu,meant
      ifsmalls=(s.lt. m*(m+mu/100.))
      meant=0
      nsample=10000000
      do i=1,nsample
         meant=meant+trandglu(m,s,mu,.false.,momglu0tab,
     &  ifsmalls)
      enddo
      meantrandglu=meant/nsample
      return
      end

      double precision function trandfermrunning(m,s,mu)
      implicit none
      logical accept
      double precision m,mu,s,ran2,m2,mu2,smin,iprel2,tspl,r,alphas,
     & alphasmax,
     & rej
      integer ifmtx
      call getMonitorFileIndex(ifmtx)
      m2=m**2
      mu2=mu**2
      smin=s-m2
      iprel2=s/smin**2
      accept=.false.
      do while(.not.accept)
         r=ran2()        
         tspl=-(1-r)/(r/mu2+iprel2)
         rej=(0.5d0+((smin+tspl)**2/2+m2*tspl)/smin**2)*
     &        (alphas(-tspl)/alphasmax())**2
         if(rej-1.gt.1D-6) then
            write(ifmtx,*) 'rej>1 in trandfermrunning:',rej
            write(ifmtx,*) m,s,mu,tspl
            close(ifmtx)
            stop
         endif
         accept=(ran2().lt.rej)
      enddo
      trandfermrunning=tspl
      return
      end

      double precision function trandfermrunningfwd(m,s,mu)
      implicit none
      logical accept
      double precision m,mu,s,ran2,m2,mu2,smin,iprel2,tspl,r,alphas,
     & alphasmax,
     & rej
      integer ifmtx
      call getMonitorFileIndex(ifmtx)
      m2=m**2
      mu2=mu**2
      smin=s-m2
      iprel2=s/(0.25*smin**2)
      accept=.false.
      do while(.not.accept)
         r=ran2()        
         tspl=-(1-r)/(r/mu2+iprel2)
         rej=(0.5d0+((smin+tspl)**2/2+m2*tspl)/smin**2)*
     &        (alphas(-tspl)/alphasmax())**2
         if(rej-1.gt.1D-6) then
            write(ifmtx,*) 'rej>1 in trandfermrunning:',rej
            write(ifmtx,*) m,s,mu,tspl
            close(ifmtx)
            stop
         endif
         accept=(ran2().lt.rej)
      enddo
      trandfermrunningfwd=tspl
      return
      end

      double precision function trandglu1running(m,s,mu)
      implicit none
      logical accept
      double precision m,mu,s,ran2,mu2,smin,iprel2,tspl,r,alphas,
     & alphasmax,rej
      mu2=mu**2
      smin=s-m**2
      iprel2=s/smin**2
      accept=.false.
      do while(.not.accept)
         r=ran2()        
         tspl=-(1-r)/(r/mu2+iprel2)
         rej=(1+tspl/smin)*(alphas(-tspl)/alphasmax())**2
         accept=(ran2().lt.rej)
      enddo
      trandglu1running=tspl
      return
      end

      double precision function trandglu2running(m,s)
      implicit none
      double precision m,s,ran2,m2,m4,smin,tinf,splus,sp2,m2s,tspl,
     & spsm
      m2=m**2
      m4=m2**2
      smin=s-m2
      tinf=smin**2/s
      splus=s+m2
      spsm=splus*smin
      sp2=splus**2
      m2s=4*m2*s
         if(ran2()*sp2.lt.m2s) then
            tspl=-ran2()*tinf
         else
            tspl=(sqrt(ran2()*spsm+m4)/s-1)*smin
         endif
      trandglu2running=tspl
      return
      end

      double precision function trandglu3running(m,s)
      implicit none
      logical accept
      double precision m,s,ran2,m2,m4,smin,smin3,den,uminsample,
     & alphas,
     &  alphasmax,rej
      integer ifmtx
      call getMonitorFileIndex(ifmtx)
      m2=m**2
      m4=m2**2
      smin=s-m2
      smin3=s-3*m2
      if(smin3.gt.0) then
         den=(smin*smin3+4*m4)
      else
         den=s/((s**2+3*m4)*m2)
      endif
      accept=.false.
      do while(.not.accept)
         uminsample=-1/(ran2()/m2+1/smin)
         if(smin3.gt.0) then
            rej=(-uminsample*smin3+4*m4)/den*
     &             (alphas(-uminsample)/alphasmax())**2
         else
            rej=(-uminsample*smin3+4*m4)*den*
     &             (alphas(-uminsample)/alphasmax())**2
         endif
         if(rej-1.gt.1d-6) then
            write(ifmtx,*) 'rej > 1 in trandglu3running; m:',m,
     &                 '/s:',s,'/u:',uminsample+m2
            close(ifmtx)
            stop
         endif
         accept=(ran2().lt.rej)
      enddo
      trandglu3running=-smin-uminsample
      return
      end

      double precision function trandglu4running(m,s,mu)
      implicit none
      logical accept
      double precision m,mu,s,ran2,m2,smin,prefact,s3m2,uminsample,rej,
     &  alphas,alphasmax
      integer ifmtx
      call getMonitorFileIndex(ifmtx)
      m2=m**2
      smin=s-m2
      prefact=s/(s+m2)**2/alphasmax()
      s3m2=s+3*m2  
      accept=.false.
      do while(.not.accept)
         uminsample=-smin*(m2/s)**ran2()
         rej=prefact*(s3m2+uminsample)*alphas(-uminsample)
         if(rej-1.gt.1d-6) then
            write(ifmtx,*) 'rej > 1 in trandglu4; mm:',m,'/mu=',mu,
     &        '/s=',s,'u=',uminsample+m2,'rej:',rej
            close(ifmtx)
            stop
         endif
         accept=(ran2().lt.rej)
      enddo
      trandglu4running=-(smin+uminsample)
      return
      end

      double precision function trandglu6posrunning(m,s,mu)
      implicit none
      logical accept
      double precision m,mu,s,ran2,mu2,smin,smin2ovmu2s,tspl,alphas,
     & alphasmax,
     & rej
      smin=s-m**2
      mu2=mu**2
      smin2ovmu2s=1+smin**2/(mu2*s)
      accept=.false.
      do while(.not.accept)
         tspl=mu2*(1-smin2ovmu2s**ran2())
         rej=alphas(-tspl)*alphas(smin+tspl)/alphasmax()**2
         accept=(ran2().lt.rej)
      enddo
      trandglu6posrunning=tspl
      return
      end

      double precision function trandglurunning(itypq,m,s,mu,ipart)
      implicit none
      logical accept
      integer itypq,ipart
      double precision m,mu,s,ran2,branch(4),r,trandglu1running,
     & trandglu2running,trandglu3running,trandglu4running,
     & trandglu6posrunning,
     & t,rej,rejglurunning
      trandglurunning=0d0
      if(ipart.ne.0) then
         if(ipart.eq.1) then
            trandglurunning=trandglu1running(m,s,mu)
         else if(ipart.eq.2) then            
            trandglurunning=trandglu2running(m,s)
         else if(ipart.eq.3) then            
            trandglurunning=trandglu3running(m,s)
         else if(ipart.eq.4) then            
            trandglurunning=trandglu4running(m,s,mu)
         else if(ipart.eq.5) then            
            trandglurunning=trandglu6posrunning(m,s,mu)
         endif
         return
      endif
      call gimmebranchglu(itypq,s,mu,branch)
      t=0d0
      accept=.false.
      do while(.not.accept)
         r=ran2()
         if(r.lt.branch(1)) then
            t=trandglu1running(m,s,mu)
         else if(r.lt.branch(2)) then
            t=trandglu2running(m,s)
         else if(r.lt.branch(3)) then
            t=trandglu3running(m,s)
         else if(r.lt.branch(4)) then
            t=trandglu4running(m,s,mu)
         else
            t=trandglu6posrunning(m,s,mu)
         endif
         rej=rejglurunning(m,s,t,mu)
         accept=(ran2().lt.rej)
      enddo
      trandglurunning=t
      return
      end

      double precision function meantrandglurunning(itypq,m,s,mu,
     &   ipart,mom2)
      implicit none
      integer nsample,i,itypq,ipart
      double precision m,s,mu,trandglurunning,meant,mom2,onet
      meant=0
      mom2=0
      nsample=10000000
      do i=1,nsample
         onet=trandglurunning(itypq,m,s,mu,ipart)
         meant=meant+onet
         mom2=mom2+onet**2
      enddo
      meantrandglurunning=meant/nsample
      mom2=mom2/nsample
      return
      end

      subroutine testtrandrunning
      implicit none
      double precision s,mu,mt,mom2,meantrandglurunning
      write(6,*) 'entering testtrandrunning'
      mu=1.10877d0
      s=5.d0
      mt=meantrandglurunning(1,1.5d0,s,mu,0,mom2)
      write(6,*) s,mt,mom2
      s=10.d0
      mt=meantrandglurunning(1,1.5d0,s,mu,0,mom2)
      write(6,*) s,mt,mom2
      s=30.d0
      mt=meantrandglurunning(1,1.5d0,s,mu,0,mom2)
      write(6,*) s,mt,mom2
      s=60.d0
      mt=meantrandglurunning(1,1.5d0,s,mu,0,mom2)
      write(6,*) s,mt,mom2
      write(6,*) 'leaving testtrandrunning'
      end

      subroutine onecollrest(ifferm,m,uvec,temp,mu,pp,iftest)
      implicit none
      integer i
      logical ifferm,iftest,ifsmalls
      double precision m,uvec(1:3),temp,mu,pp(0:3),m2,q,u,ctheta,
     & qrandferm,qrandglu,cthetarand,qvec(1:3),porth(1:3),bi1,bi2,s,t,
     & trandferm,trandglu,prelcm2,cthetacm,ppt2,ppt,ep,ppl,ppl2,qpl,
     & qp,error,momglu0tab(6)
      integer ifmtx
      call getMonitorFileIndex(ifmtx)
      m2=m**2
      u=sqrt(uvec(1)**2+uvec(2)**2+uvec(3)**2)
      if(isnan(u))then
         write(ifmtx,*) 'isnan u in onecollrest'
         write(ifmtx,*) uvec(1),uvec(2),uvec(3)
         close(ifmtx)
         stop
      endif
      if(ifferm) then
         q=qrandferm(m,u,temp,mu)
      else
         q=qrandglu(m,u,temp,mu,ifsmalls,momglu0tab)
      endif   
      ctheta=cthetarand(u,temp,q)
      call giveperpvec(uvec,porth)
      bi1=ctheta/u
      bi2=sqrt(1-ctheta**2)
      do i=1,3
        qvec(i)=q*(bi1*uvec(i)+bi2*porth(i))
      enddo
      s=m*(m+2*q)
      prelcm2=(s-m2)**2/(4*s)
      if(ifferm) then
         t=trandferm(m,s,mu)
      else
         t=trandglu(m,s,mu,.true.,momglu0tab,ifsmalls)
      endif
      cthetacm=1+t/(2*prelcm2)
      ppt2=prelcm2*(1-cthetacm**2)
      ppt=sqrt(ppt2)
      ep=(2*m2-t)/(2*m)
      ppl2=max(0.D0,ep**2-m2-ppt2)
      ppl=sqrt(ppl2)
      ep=sqrt(ppl2+m2+ppt2)
      if(iftest) then
         qpl=q-ppl
         qp=sqrt(qpl**2+ppt2)
         error=abs(s-m2-2*(qp*ep-ppl*qpl+ppt2))
         if(error.gt.1e-8) then
            write(8,*) 'on shellness error in onecollrest:',error
         endif
      endif
      call giveperpvec(qvec,porth)
      bi1=ppl/q
      do i=1,3
         pp(i)=bi1*qvec(i)+ppt*porth(i)
      enddo
      pp(0)=ep
      return
      end

      subroutine getoneqperprandcst(itypq,ifferm,temp,m,s,mu,t,ifsmalls,
     & momglu0tab)
      implicit none
      logical ifferm,ifsmalls
      integer itypq
      double precision temp,m,s,mu,trandfermforward,trandglu,ran2,
     & t,ratio,qperp,
     & momglu0tab(6),qperpmax,ratioratemax
      integer ifmtx
      call getMonitorFileIndex(ifmtx)
      qperpmax=(s-m**2)/(2*sqrt(s))
      call gimmeratioratesoneqperp(itypq,qperpmax,temp,ratioratemax)
       if(s.le.(m**2+mu*m/100.)) then
         ifsmalls=.true.
      else
         ifsmalls=.false.
      endif
 10   continue
      if(ifferm) then
         t=trandfermforward(m,s,mu)
      else
         t=trandglu(m,s,mu,.true.,momglu0tab,ifsmalls)
      endif
      qperp=sqrt(-t)       
      if(qperp.gt.qperpmax) then
         goto 10
      endif
      call gimmeratioratesoneqperp(itypq,qperp,temp,ratio)
      if(ratio.gt.ratioratemax) then
         write(ifmtx,*) 'unexpected ratio in getoneqperprandcst'
         write(ifmtx,*) s,temp
         write(ifmtx,*) qperpmax, ratioratemax
         write(ifmtx,*) qperp,ratio
         close(ifmtx)
         stop
         goto 10
      endif
      if(ran2()*ratioratemax.gt.ratio) then
         goto 10
      endif
      return
      end


      subroutine genoneqperpandinvariantexitcst(itypq,ifferm,m,pl,q,
     &  temp,
     &  mu,ifsmalls,momglu0tab,ifailglob,s,qperp2,x,kperp,pptx,ppty,
     &  pplnew)
      use PBGgluprops
      implicit none
      logical ifferm,ifsmalls
      integer itypq,ifailglob,ifailLPM,ifailphase,nbcol
      double precision m,pl,q,s,temp,mu,ratio_ferm,ratio_glu,
     &  ratio,
     &  t,qperp2,x,kperp,cosphi,sinphi,pplus,qperp,pptx,ppty,ppt2,
     &  mpt2,
     &  pplnew,ran2,momglu0tab(6),mg
      integer ifmtx
      call getMonitorFileIndex(ifmtx)
      ifailglob=0
      call gimmeratiorates(itypq,q,temp,ratio_ferm,ratio_glu)
      if(ifferm) then
         ratio=ratio_ferm
      else
         ratio=ratio_glu
      endif
      if(ratio.gt.1.D0) then
         write(ifmtx,*) 'ratio=',ratio
         write(ifmtx,*) itypq,q,temp
         close(ifmtx)
         stop
      endif
      if(ran2().lt.ratio) then
         nbcol=1
      else
         nbcol=0
         ifailglob=1
         return
      endif
      s=m*(m+2*q)
      call getoneqperprandcst(itypq,ifferm,temp,m,s,mu,t,ifsmalls,
     & momglu0tab)

      qperp2=-t
      call get1x_veck_qQ_radiat(itypq,pl,temp,m,qperp2,x,kperp,cosphi,
     &   sinphi,ifailLPM)
      if(ifailLPM.ne.0) then
         ifailglob=2
         return
      else         
         pplus=m
         qperp=sqrt(qperp2)
         pptx=qperp-kperp*cosphi
         ppty=-kperp*sinphi
         ppt2=pptx**2+ppty**2
         mpt2=m**2+ppt2
         mg=c1*temp
         call gennewPradiat(x,s,pplus,mpt2,mg,kperp,qperp,pplnew,
     &      ifailphase)
         if(ifailphase.le.1) then
            ifailglob=0
         else
            ifailglob=3
         endif
      endif
      return
      end

      subroutine testgenoneqperpandinvariantexitcst
      implicit none
      integer i,itypq,ifailglob,nbsamp,nbsuccess,isamp
      logical ifferm,ifsmalls
      double precision mQ,pl,q,temp,mu,s,qperp2,x,kperp,
     & pptx,
     & ppty,pplnew,xav,muoftemp,plav,pplnewshift,branching,
     & momglu0tab(6)
      temp=0.3
      nbsamp=10000000
      ifferm=.true.
      mu=sqrt(0.15)*muoftemp(temp)
      itypq=2
      mQ=5.1D0
      pl=20.D0
      open(9,file='branching_of_s_SQCD.res',status='new')
      write(9,*) 'all rad quantities rescaled by simgam_col'      
      write(9,*) mQ,temp,mu
      write(9,*) 'q s sig_rad/sig_col <x> I1/sig_col'
      do i=1,60
         s=mQ**2*exp(0.1*i)
         q=(s-mQ**2)/(2*mQ)
         xav=0.D0
         plav=0.D0
         nbsuccess=0
         do isamp=1,nbsamp
            call genoneqperpandinvariantexitcst(itypq,ifferm,mQ,pl,q,
     & temp,
     &         mu,ifsmalls,momglu0tab,ifailglob,s,qperp2,x,kperp,pptx,
     & ppty,
     &          pplnew)
            if(ifailglob.eq.0) then
               nbsuccess=nbsuccess+1
               xav=xav+x
               plav=plav+pplnewshift
            endif
         enddo
         write(6,*) s,xav/nbsuccess
         branching=real(nbsuccess)/nbsamp
         write(9,*) q,s,branching,xav/nbsuccess,xav/nbsamp
      enddo
      close(9)
      end


      subroutine collrestrad(itypq,ifferm,m,pl,uvec,temp,mu,nbcol,pp,s,
     & qperp2,
     &  x)
      implicit none
      logical ifferm,ifsmalls
      integer itypq,nbcol,ifailglob,j      
      double precision m,m2,pl,uvec(3),u2,u,temp,mu,q,qrandferm,
     & qrandglu,s,
     & qperp2,x,kperp,pptx,ppty,pplnew,pp(0:3),ctheta,
     & cthetarand,
     & bi1,bi2,qvec(3),porth(3),vecet(3),momglu0tab(6)
      integer ifmtx
      call getMonitorFileIndex(ifmtx)
      pp(0)=m
      pp(1)=0.D0
      pp(2)=0.D0
      pp(3)=0.D0
      m2=m**2
      u2=uvec(1)**2+uvec(2)**2+uvec(3)**2
      u=sqrt(u2)
      if(ifferm) then
         q=qrandferm(m,u,temp,mu)
      else
         q=qrandglu(m,u,temp,mu,ifsmalls,momglu0tab)
      endif
      call genoneqperpandinvariantexitcst(itypq,ifferm,m,pl,q,temp,mu,
     &  ifsmalls,momglu0tab,ifailglob,s,qperp2,x,kperp,pptx,ppty,pplnew)
      if(ifailglob.ne.0) then
         nbcol=0
         return
      else
         nbcol=1
         ctheta=cthetarand(u,temp,q)
         call giveperpvec(uvec,porth)
         bi1=ctheta/u
         bi2=sqrt(1-ctheta**2)
         do j=1,3
            qvec(j)=bi1*uvec(j)+bi2*porth(j)
         enddo
         if(abs(qvec(1)**2+qvec(2)**2+qvec(3)**2-1).gt.1E-6) then
           write(ifmtx,*)'consist. test 1 failed;collresrunningradnew'
            close(ifmtx)
            stop
         endif
         call giveperpvec(qvec,porth)
         vecet(1)=qvec(2)*porth(3)-porth(2)*qvec(3)
         vecet(2)=qvec(3)*porth(1)-porth(3)*qvec(1)
         vecet(3)=qvec(1)*porth(2)-porth(1)*qvec(2)
         if(abs(vecet(1)**2+vecet(2)**2+vecet(3)**2-1).gt.1E-6) then
           write(ifmtx,*)'consist. test 2 failed;collresrunningradnew'
            close(ifmtx)
            stop
         endif
         pp(0)=m**2
         do j=1,3
            pp(j)=-pplnew*qvec(j)+pptx*porth(j)+ppty*vecet(j)
            pp(0)=pp(0)+pp(j)**2
         enddo
         pp(0)=sqrt(pp(0))
      endif
      return
      end


      subroutine onecoll(itypq,ifferm,m,p,temp,mu,ifrad,nbcollsucc,pp,
     &    x,s)
      implicit none
      logical ifferm,ifrad
      integer i,nbcol,itypq,nbcollsucc
      double precision m,p(0:3),temp,mu,pp(0:3),u4vec(0:3),uvec(3),
     & pprest(0:3),s,qperp2,pl,x 
      integer ifmtx
      call getMonitorFileIndex(ifmtx)
      nbcollsucc=0
      u4vec(0)=p(0)/m
      do i=1,3
         u4vec(i)=-p(i)/m
         uvec(i)=u4vec(i)
      enddo
      if(ifrad) then
         pl=sqrt(p(1)**2+p(2)**2+p(3)**2)
         call collrestrad(itypq,ifferm,m,pl,uvec,temp,mu,nbcol,pprest,s,
     &     qperp2,x)
         if(nbcol.eq.0) then
            return
         else 
            if(nbcol.ge.2) then
               write(ifmtx,*) 'need several successive coll;onecoll !'
               close(ifmtx)
               stop
            else
               call lorentz2(u4vec,pprest,pp)
               nbcollsucc=1
            endif
         endif
      else
         nbcollsucc=1
         call onecollrest(ifferm,m,uvec,temp,mu,pprest,.false.)
         call lorentz2(u4vec,pprest,pp)
      endif
      return
      end


      subroutine testonecoll
      implicit none
      logical ifferm,ifrad
      integer i,j,nsample,nbcol,itypq
      double precision pi,m,s,tbpin(0:30),pin(0:3),mu,temp,pp(0:3),
     & meandp,
     & momlong,momtrans,muoftemp,nbsucces,ratef,rateg,pv,
     & meande,x,avx,sav 
      parameter (pi=3.14159265358979323844d0)
 10   format(F10.5,' ', 12 (F12.7,' '))
      write(6,*) 'entering test onecol'
      m=1.5d0
      itypq=1
      temp=0.3d0
      mu=0.387299*muoftemp(temp)
      tbpin(0)=0.1D0
      do i=1,4
         tbpin(i)=i
      enddo
      do i=5,20
         tbpin(i)=4*(i-4)+4
      enddo
      do i=21,30
         tbpin(i)=10*(i-20)+68
      enddo
      
      write(6,*) 'generating transport coefficients'
      ifrad=.true.
      write(6,*) 'for radiative processes'      
      ifferm=.true.
      nsample=10000000
      write(6,*) 'for fermions:'
      open(9,
     & file='av_moment_modelC_radSQCD_ferm_c_T300_1minx_boundnog.res',
     &    status='new')
      write(9,*) m,temp,mu,nsample
      write(9,*)'p |rate|<s>|dp_z/coll|<x>|dp_z^2/coll|dp_t^2/coll| 
     &  dp_z/dt | dp_t^2/dt | dE/dt'  
      do i=0,30
         pin(1)=0
         pin(2)=0
         pin(3)=tbpin(i)
         pv=sqrt(pin(1)**2+pin(2)**2+pin(3)**2)
         pin(0)=sqrt(pv**2+m**2)
         meandp=0.D0
         meande=0.D0
         momlong=0.
         momtrans=0.
         avx=0.D0
         sav=0.D0
         nbsucces=0
         call gimmerates(itypq,pv,temp,ratef,rateg)
         do j=1,nsample
            call onecoll(itypq,ifferm,m,pin,temp,mu,ifrad,nbcol,pp,x,s) 
            if(nbcol.eq.1) then
               nbsucces=nbsucces+1
               meandp=meandp+pp(3)
               momlong=momlong+(pin(3)-pp(3))**2
               momtrans=momtrans+(pp(1)**2+pp(2)**2)
               meande=meande+pin(0)-pp(0)
               avx=avx+x
               sav=sav+s
            endif
         enddo
         meandp=pin(3)-meandp/nbsucces
         momlong=momlong/nbsucces
         momtrans=momtrans/nbsucces
         meande=meande/nbsucces
         avx=avx/nbsucces
         sav=sav/nbsucces
         nbsucces=nbsucces/nsample
         ratef=nbsucces*ratef
         write(6,10) pin(3),ratef,sav,meandp,avx,momlong,momtrans,
     &      ratef*meandp,ratef*momtrans,ratef*meande
         write(9,10) pin(3),ratef,sav,meandp,avx,momlong,momtrans,
     &      ratef*meandp,ratef*momtrans,ratef*meande
      enddo
      close(9)
      ifferm=.false.
      nsample=10000000
      write(6,*) 'for gluons:'
      open(9,file='av_moment_fixedalpha_rad_gluo.res',status='new')
      write(9,*) m,temp,mu,nsample
      write(9,*) 'p |rate|dp_z/coll|dp_z^2/coll|dp_t^2/coll|dp_z/dt| 
     & dp_t^2/dt'    
      do i=0,20
         pin(1)=0
         pin(2)=0
         pin(3)=tbpin(i)
         pv=sqrt(pin(1)**2+pin(2)**2+pin(3)**2)
         pin(0)=sqrt(pv**2+m**2)
         meandp=0.D0
         meande=0.D0
         momlong=0.
         momtrans=0.
         nbsucces=0
         call gimmerates(itypq,pv,temp,ratef,rateg)
         do j=1,nsample
           call onecoll(itypq,ifferm,m,pin,temp,mu,ifrad,nbcol,pp,x,s)  
            if(nbcol.eq.1) then
               nbsucces=nbsucces+1
               meandp=meandp+pp(3)
               momlong=momlong+(pin(3)-pp(3))**2
               momtrans=momtrans+(pp(1)**2+pp(2)**2)
               meande=meande+pin(0)-pp(0)
               
            endif
         enddo
         meandp=pin(3)-meandp/nbsucces
         momlong=momlong/nbsucces
         momtrans=momtrans/nbsucces
         meande=meande/nbsucces
         nbsucces=nbsucces/nsample
         rateg=nbsucces*rateg
         write(6,10) pin(3),rateg,meandp,momlong,momtrans,
     &      rateg*meandp,rateg*momtrans,rateg*meande
         write(9,10) pin(3),rateg,meandp,momlong,momtrans,
     &      rateg*meandp,rateg*momtrans,rateg*meande
      enddo
      close(9)     
      return
      end


      subroutine onecollrestrunning(itypq,ifferm,m,uvec,temp,mu,pp,
     & iftest)
      implicit none
      integer i,itypq
      logical ifferm,iftest
      double precision m,uvec(1:3),temp,mu,pp(0:3),m2,q,u,ctheta,
     & qrandfermrunning,qrandglurunning,cthetarand,qvec(1:3),porth(1:3),
     & bi1,bi2,s,t,trandfermrunning,trandglurunning,prelcm2,cthetacm,
     & ppt2,
     & ppt,ep,ppl,ppl2,qpl,qp,error
      integer ifmtx
      call getMonitorFileIndex(ifmtx)
      m2=m**2
      u=sqrt(uvec(1)**2+uvec(2)**2+uvec(3)**2)
      if(isnan(u/u))then
         write(ifmtx,*) 'isnan u in onecollrest running'
         write(ifmtx,*) uvec(1),uvec(2),uvec(3)
         write(ifmtx,*) u
         close(ifmtx)
         stop
      endif
      u=sqrt(uvec(1)**2+uvec(2)**2+uvec(3)**2)
      if(ifferm) then
         q=qrandfermrunning(itypq,u,temp,mu)
      else
         q=qrandglurunning(itypq,u,temp,mu)
      endif   
      ctheta=cthetarand(u,temp,q)
      call giveperpvec(uvec,porth)
      bi1=ctheta/u
      bi2=sqrt(1-ctheta**2)
      do i=1,3
        qvec(i)=q*(bi1*uvec(i)+bi2*porth(i))
      enddo
      s=m*(m+2*q)
      prelcm2=(s-m2)**2/(4*s)
      if(ifferm) then
         t=trandfermrunning(m,s,mu)
      else
         t=trandglurunning(itypq,m,s,mu,0)
      endif
      cthetacm=1+t/(2*prelcm2)
      ppt2=prelcm2*(1-cthetacm**2)
      ppt=sqrt(ppt2)
      ep=(2*m2-t)/(2*m)
      ppl2=max(0.D0,ep**2-m2-ppt2)
      ppl=sqrt(ppl2)
      ep=sqrt(ppl2+m2+ppt2)
      if(iftest) then
         qpl=q-ppl
         qp=sqrt(qpl**2+ppt2)
         error=abs(s-m2-2*(qp*ep-ppl*qpl+ppt2))
         if(error.gt.1e-8) then
            write(8,*) 'on shellness error in onecollrestrunning:',error
         endif
      endif
      call giveperpvec(qvec,porth)
      bi1=ppl/q
      do i=1,3
         pp(i)=bi1*qvec(i)+ppt*porth(i)
      enddo
      pp(0)=ep
      return
      end

    
      subroutine getoneqperprandbasic(itypq,ifferm,temp,m,s,mu,nit,t,
     &  qperp)
      implicit none
      logical ifferm
      integer itypq,nit,i
      double precision temp,m,s,mu,trandfermrunning,trandglurunning,
     & ran2,
     & ttry,t,ratiotry,ratio,qperp,qperptry
      if(ifferm) then
         t=trandfermrunning(m,s,mu)
      else
         t=trandglurunning(itypq,m,s,mu,0)
      endif
      qperp=sqrt(-t)
      call gimmeratioratesoneqperp(itypq,qperp,temp,ratio)
      do 10 i=1,nit
         if(ifferm) then
            ttry=trandfermrunning(m,s,mu)
         else
            ttry=trandglurunning(itypq,m,s,mu,0)
         endif
         qperptry=sqrt(-ttry)
         call gimmeratioratesoneqperp(itypq,qperptry,temp,ratiotry)
         if(ratiotry.lt.ratio)then
            if(ratio*ran2().gt.ratiotry) goto 10
         endif
         t=ttry
         qperp=qperptry
         ratio=ratiotry
 10   continue
      return
      end

      subroutine collrestrunningradbasic(itypq,ifferm,m,uvec,temp,mu,
     & nit,
     &  nbcol,qp)
      implicit none
      integer nit,j,itypq,nbcol,nbcolmax
      parameter (nbcolmax=10)
      logical ifferm
      double precision m,uvec(1:3),temp,mu,qp(nbcolmax,3),m2,q,u,
     & qrandfermrunning,qrandglurunning,ratio_ferm,ratio_glu,ran2,ratio,
     & porth(1:3),qperp,s,t
      integer ifmtx
      call getMonitorFileIndex(ifmtx)
      m2=m**2
      u=sqrt(uvec(1)**2+uvec(2)**2+uvec(3)**2)
      if(ifferm) then
         q=qrandfermrunning(itypq,u,temp,mu)
      else
         q=qrandglurunning(itypq,u,temp,mu)
      endif 
      call gimmeratiorates(itypq,q,temp,ratio_ferm,ratio_glu)
      if(ifferm) then
         ratio=ratio_ferm
      else
         ratio=ratio_glu
      endif
      nbcol=int(ratio)
      if(ran2().lt.ratio-nbcol) then
         nbcol=nbcol+1
      endif
      if(nbcol.eq.0) return
      if(nbcol.gt.nbcolmax) then
         write(ifmtx,*) 'too large nbcol generated;collrestrunningrad:',
     &  nbcol
         close(ifmtx)
         stop
      endif
      s=m*(m+2*q)
      do j=1,nbcol
         call getoneqperprandbasic(itypq,ifferm,temp,m,s,mu,nit,t,qperp)
         call giveperpvec(uvec,porth)
         qp(j,1)=qperp*porth(1)
         qp(j,2)=qperp*porth(2)
         qp(j,3)=qperp*porth(3)
      enddo
      return
      end


      subroutine genonekinfinrad(m,uvec,u2,q,qvec,t,prelcm2,deltapperp)
      implicit none 
      integer i
      double precision t,m,uvec(3),u2,q,qvec(3),prelcm2,ppt,pplovq,
     & porth(3),
     & deltap(3),deltapperp(3),proj 
      ppt=sqrt(-t-t**2/(4*prelcm2))
      pplovq=-t*(m+q)/(2*m*q**2)
      call giveperpvec(qvec,porth)
      proj=0
      do i=1,3
         deltap(i)=pplovq*qvec(i)+ppt*porth(i)
         proj=proj+deltap(i)*uvec(i)
      enddo
      proj=proj/u2
      do i=1,3
         deltapperp(i)=deltap(i)-proj*uvec(i)
      enddo
      return
      end

      subroutine getoneqperprand(faithfull,ifrunning,itypq,ifferm,temp,
     & m,uvec,
     & u2,q,s,prelcm2,qvec,mu,nit,t,qperpvec)
      implicit none
      logical ifrunning,faithfull,ifferm,ifsmalls
      integer itypq,nit,i
      double precision temp,m,q,s,prelcm2,uvec(3),u2,qvec(3),mu,
     & trandferm,trandglu,trandfermrunning,trandglurunning,ran2,
     & ttry,t,ratiotry,ratio,qperp,qperptry,qperpvec(3),deltapperp(3),
     & deltapperptry(3),momglu0tab(6)
      if(s.le.(m**2+mu*m/100.)) then
         ifsmalls=.true.
      else
         ifsmalls=.false.
      endif
      if(ifrunning)then
         if(ifferm) then
            t=trandfermrunning(m,s,mu)
         else
            t=trandglurunning(itypq,m,s,mu,0)
         endif
      else
         if(ifferm) then
            t=trandferm(m,s,mu)
         else
            t=trandglu(m,s,mu,.false.,momglu0tab,ifsmalls)
         endif
      endif
      if(faithfull) then 
         qperp=sqrt(-t)       
      else
         call genonekinfinrad(m,uvec,u2,q,qvec,t,prelcm2,deltapperp)
         qperp=sqrt(deltapperp(1)**2+deltapperp(2)**2+deltapperp(3)**2)
      endif
      call gimmeratioratesoneqperp(itypq,qperp,temp,ratio)
      do 10 i=1,4*nit
         if(ifrunning) then
            if(ifferm) then
               ttry=trandfermrunning(m,s,mu)
            else
               ttry=trandglurunning(itypq,m,s,mu,0)
            endif
         else
            if(ifferm) then
               ttry=trandferm(m,s,mu)
            else
               ttry=trandglu(m,s,mu,.true.,momglu0tab,ifsmalls)
            endif
         endif
         if(faithfull) then 
            qperptry=sqrt(-ttry)       
         else
            call genonekinfinrad(m,uvec,u2,q,qvec,ttry,prelcm2,
     &           deltapperptry)
            qperptry=sqrt(deltapperptry(1)**2+deltapperptry(2)**2+
     &                 deltapperptry(3)**2)
         endif
         call gimmeratioratesoneqperp(itypq,qperptry,temp,ratiotry)
         if(ratiotry.lt.ratio)then
            if(ratio*ran2().gt.ratiotry) goto 10
         endif
         t=ttry
         deltapperp(1)=deltapperptry(1)
         deltapperp(2)=deltapperptry(2)
         deltapperp(3)=deltapperptry(3)
         ratio=ratiotry
 10   continue
      if(faithfull) then
         call genonekinfinrad(m,uvec,u2,q,qvec,t,prelcm2,qperpvec)
      else
         qperpvec(1)=deltapperp(1)
         qperpvec(2)=deltapperp(2)
         qperpvec(3)=deltapperp(3)
      endif
      return
      end


      subroutine getoneqperprandsimp(itypq,ifferm,temp,m,s,mu,t)
      implicit none
      logical ifferm
      integer itypq
      double precision temp,m,s,mu,ran2,
     & trandfermrunningfwd,trandglurunning,t,ratio,qperp,qperpmax,
     & ratioratemax
      integer ifmtx
      call getMonitorFileIndex(ifmtx)
      qperpmax=(s-m**2)/(2*sqrt(s))
      call gimmeratioratesoneqperp(itypq,qperpmax,temp,ratioratemax)
 10   continue
      if(ifferm) then
         t=trandfermrunningfwd(m,s,mu)
      else
         t=trandglurunning(itypq,m,s,mu,0)
      endif
      qperp=sqrt(-t)       
      if(qperp.gt.qperpmax) then
         goto 10
      endif
      call gimmeratioratesoneqperp(itypq,qperp,temp,ratio)
      if(ratio.gt.ratioratemax) then
         write(ifmtx,*) 'unexpected ratio in getoneqperprandsimp'
         write(ifmtx,*) s,temp
         write(ifmtx,*) qperpmax, ratioratemax 
         write(ifmtx,*) qperp, ratio
         close(ifmtx)
         stop
         goto 10
      endif
      if(ran2()*ratioratemax.gt.ratio) then
         goto 10
      endif
      return
      end


      subroutine collrestrunningrad(itypq,ifferm,m,uvec,temp,mu,nit,
     & nbcol,
     &  qpvec,s,t)
      implicit none
      integer i,j,itypq,nbcol,nbcolmax,jp,nit
      parameter (nbcolmax=10)
      logical ifferm
      double precision m,uvec(1:3),temp,mu,m2,q,u,
     & qrandfermrunning,
     & qrandglurunning,ratio_ferm,ratio_glu,ran2,ratio,ctheta,
     & cthetarand,
     & qvec(1:3),porth(1:3),qpvec(nbcolmax,3),s,t,u2,prelcm2,bi1,
     & bi2,
     & qperpvec(3)
      integer ifmtx
      call getMonitorFileIndex(ifmtx)
      m2=m**2
      u2=uvec(1)**2+uvec(2)**2+uvec(3)**2
      u=sqrt(u2)
      if(ifferm) then
         q=qrandfermrunning(itypq,u,temp,mu)
      else
         q=qrandglurunning(itypq,u,temp,mu)
      endif 
      call gimmeratiorates(itypq,q,temp,ratio_ferm,ratio_glu)
      if(ifferm) then
         ratio=ratio_ferm
      else
         ratio=ratio_glu
      endif
      nbcol=int(ratio)
      if(ran2().lt.ratio-nbcol) then
         nbcol=nbcol+1
      endif
      if(nbcol.eq.0) return
      s=m*(m+2*q)
      prelcm2=(s-m2)**2/(4*s)
      if(nbcol.gt.nbcolmax) then
         write(ifmtx,*) 'too large nbcol generated;collrestrunningrad:',
     &      nbcol
         close(ifmtx)
         stop
      endif
      jp=0
      do j=1,nbcol
         ctheta=cthetarand(u,temp,q)
         call giveperpvec(uvec,porth)
         bi1=ctheta/u
         bi2=sqrt(1-ctheta**2)
         do i=1,3
            qvec(i)=q*(bi1*uvec(i)+bi2*porth(i))
         enddo
         call getoneqperprand(.false.,.true.,itypq,ifferm,temp,m,uvec,
     &        u2,q,s,
     &        prelcm2,qvec,mu,nit,t,qperpvec)
         if(-t/s.gt. 0.25) then
         else
            jp=jp+1
            do i=1,3
               qpvec(jp,i)=qperpvec(i)
            enddo
         endif
      enddo
      nbcol=jp
      return
      end


      subroutine genoneqperpandinvariantexit(itypq,ifferm,m,pl,q,temp,
     &  mu,
     &  ifailglob,s,qperp2,x,kperp,pptx,ppty,pplnew)
      use PBGgluprops
      implicit none
      logical ifferm
      integer itypq,ifailglob,ifailLPM,ifailphase,ifailphase2,nbcol
      double precision m,pl,q,s,temp,mu,ratio_ferm,ratio_glu,
     &  ratio,
     &  t,qperp2,x,kperp,cosphi,sinphi,pplus,qperp,pptx,ppty,ppt2,
     &  mpt2,
     &  pplnew,pplnewshift,ran2,pplnewbis,yorig,yshift
      double precision mg
      integer ifmtx
      call getMonitorFileIndex(ifmtx)
      ifailglob=0
      call gimmeratiorates(itypq,q,temp,ratio_ferm,ratio_glu)
      if(ifferm) then
         ratio=ratio_ferm
      else
         ratio=ratio_glu
      endif
      if(ratio.gt.1.D0) then
         write(ifmtx,*) 'ratio=',ratio
         write(ifmtx,*) itypq,q,temp
         close(ifmtx)
         stop
      endif
      if(ran2().lt.ratio) then
         nbcol=1
      else
         nbcol=0
         ifailglob=1
         return
      endif
      s=m*(m+2*q)
      call getoneqperprandsimp(itypq,ifferm,temp,m,s,mu,t)
      qperp2=-t
      call get1x_veck_qQ_radiat(itypq,pl,temp,m,qperp2,x,kperp,cosphi,
     &   sinphi,ifailLPM)
      if(ifailLPM.ne.0) then
         ifailglob=2
         return
      else         
         pplus=m
         qperp=sqrt(qperp2)
         pptx=qperp-kperp*cosphi
         ppty=-kperp*sinphi
         ppt2=pptx**2+ppty**2
         mpt2=m**2+ppt2
         mg=c1*temp
         call gennewPradiat(x,s,pplus,mpt2,mg,kperp,qperp,pplnew,
     &       ifailphase)
         if(ifailphase.le.1) then
            ifailglob=0
            pplus=pl+sqrt(m**2+pl**2)
            call gennewPradiat(x,s,pplus,mpt2,mg,kperp,qperp,
     &   pplnewbis,
     &       ifailphase2)
            if(ifailphase2.eq.2) then
               write(ifmtx,*) 'genoneqperpandinvariantexit:inconsistent'
               close(ifmtx)
               stop
            endif
            yorig=log(sqrt(1+(pl/m)**2)+pl/m)
            yshift=log(sqrt(1+pplnew**2/mpt2)+pplnew/sqrt(mpt2))
            yorig=yorig+yshift
            pplnewshift=sqrt(mpt2)*(exp(yorig)-exp(-yorig))/2.D0
            if(abs(pplnewbis/pplnewshift-1).gt.1E-8) then
            endif
         else
            ifailglob=3
         endif
      endif
      return
      end


      subroutine collrestrunningradnew(itypq,ifferm,m,pl,uvec,
     &  temp,mu,nbcol,pp,s,qperp2)
      use JA
      implicit none
      logical ifferm
      integer itypq,nbcol,ifailglob,j      
      double precision m,m2,pl,uvec(3),u2,u,temp,mu,q,qrandfermrunning,
     & qrandglurunning,s,qperp2,x,kperp,pptx,ppty,pplnew,
     & pp(0:3),ctheta,cthetarand,bi1,bi2,qvec(3),porth(3),vecet(3)
      integer ifmtx
      call getMonitorFileIndex(ifmtx)
      pp(0)=m
      pp(1)=0.D0
      pp(2)=0.D0
      pp(3)=0.D0
      m2=m**2
      u2=uvec(1)**2+uvec(2)**2+uvec(3)**2
      u=sqrt(u2)
      if(ifferm) then
         q=qrandfermrunning(itypq,u,temp,mu)
      else
         q=qrandglurunning(itypq,u,temp,mu)
      endif
      call genoneqperpandinvariantexit(itypq,ifferm,m,pl,q,temp,mu,
     &  ifailglob,s,qperp2,x,kperp,pptx,ppty,pplnew)
      if(ifailglob.ne.0) then
         nbcol=0
         return
      else
         if(ifferm)then
         jaityp=3
         else
         jaityp=4
         endif 
         
         nbcol=1
         ctheta=cthetarand(u,temp,q)
         call giveperpvec(uvec,porth)
         bi1=ctheta/u
         bi2=sqrt(1-ctheta**2)
         do j=1,3
            qvec(j)=bi1*uvec(j)+bi2*porth(j)
         enddo
         if(abs(qvec(1)**2+qvec(2)**2+qvec(3)**2-1).gt.1E-6) then
           write(ifmtx,*)'consist. test 1 failed;collresrunningradnew'
            close(ifmtx)
            stop
         endif
         call giveperpvec(qvec,porth)
         vecet(1)=qvec(2)*porth(3)-porth(2)*qvec(3)
         vecet(2)=qvec(3)*porth(1)-porth(3)*qvec(1)
         vecet(3)=qvec(1)*porth(2)-porth(1)*qvec(2)
         if(abs(vecet(1)**2+vecet(2)**2+vecet(3)**2-1).gt.1E-6) then
           write(ifmtx,*)'consist. test 2 failed;collresrunningradnew'
            close(ifmtx)
            stop
         endif
         pp(0)=m**2
         do j=1,3
            pp(j)=-pplnew*qvec(j)+pptx*porth(j)+ppty*vecet(j)
            pp(0)=pp(0)+pp(j)**2
         enddo
         pp(0)=sqrt(pp(0))
      endif
      return
      end


      subroutine testonecollrest
      use PBGforradiatGB
      implicit none
      logical ifferm
      integer i,j,k,nsample,nsamplesucc,nbcol,nbcolmax,tbnit(9),ip
      parameter(nbcolmax=10)
      double precision pi,m,pin(0:3),uvec(3),mu,temp,
     & qpvec(nbcolmax,3),
     & momtrans(9),mdebeffinterp,distqp(9,50)
      parameter (pi=3.14159265358979323844d0)
 10   format(F10.5,' ', 9 (F12.7,' '))
      temp=0.3d0
      m=1.5d0
      do k=1,9
         tbnit(k)=(k-1)*5
      enddo
      mu=mdebeffinterp(temp)
      open(9,file='qperp_av_ferm.res',status='new')
      write(9,*) m
      write(9,*) temp
      ifferm=.true.
      nsample=100000
      write(6,*) 'for fermions:'
      do i=5,30,5
         write(6,*) 'pin=',i,' GeV/c)'
         pin(1)=0
         pin(2)=0
         pin(3)=i
         pin(0)=sqrt(pin(3)**2+m**2)
         do j=1,3
            uvec(j)=-pin(j)/m
         enddo
         do k=1,9
            momtrans(k)=0.
            nsamplesucc=0
            do j=1,nsample
               call collrestrunningradbasic(1,ifferm,m,uvec,
     &              temp,mu,tbnit(k),
     &              nbcol,qpvec)
               if(nbcol.ge.1) then
                  nsamplesucc=nsamplesucc+1
                  momtrans(k)=momtrans(k)+sqrt(qpvec(1,1)**2
     &            +qpvec(1,2)**2)
               endif
            enddo
            write(6,*) '# of successes:',nsamplesucc
            if(nsamplesucc.ge.1) then               
               momtrans(k)=momtrans(k)/nsamplesucc
            else
               momtrans(k)=0.d0
            endif
         enddo
         write(9,10) pin(3),(momtrans(k),k=1,9)
      enddo
      close(9)
      open(9,file='qperp_av_gluo.res',status='new')
      write(9,*) m
      write(9,*) temp
      ifferm=.false.
      nsample=100000
      write(6,*) 'for gluons:'
      do i=5,30,5
         write(6,*) 'pin=',i,' GeV/c)'
         pin(1)=0
         pin(2)=0
         pin(3)=i
         pin(0)=sqrt(pin(3)**2+m**2)
         do j=1,3
            uvec(j)=-pin(j)/m
         enddo
         do k=1,9
            momtrans(k)=0.
            nsamplesucc=0
            do j=1,nsample
               call collrestrunningradbasic(1,ifferm,m,uvec,temp,mu,
     &              tbnit(k),
     &              nbcol,qpvec)
               if(nbcol.ge.1) then
                  nsamplesucc=nsamplesucc+1
                  momtrans(k)=momtrans(k)+sqrt(qpvec(1,1)**2
     &              +qpvec(1,2)**2)
               endif
            enddo
            write(6,*) '# of successes:',nsamplesucc
            if(nsamplesucc.ge.1) then               
               momtrans(k)=momtrans(k)/nsamplesucc
            else
               momtrans(k)=0.d0
            endif
         enddo
         write(9,10) pin(3),(momtrans(k),k=1,9)
      enddo
      close(9)
      ifferm=.true.
      nsample=5000000
      write(6,*) 'for fermions:'
      pin(1)=0
      pin(2)=0
      pin(3)=10
      do j=1,3
         uvec(j)=-pin(j)/m
      enddo
      do k=1,9
         nsamplesucc=0
         do i=1,50
            distqp(k,i)=0
         enddo
         do j=1,nsample
            call collrestrunningradbasic(1,ifferm,m,uvec,temp,mu,
     &             tbnit(k),
     &              nbcol,qpvec)
            do i=1,nbcol
               nsamplesucc=nsamplesucc+1
               ip=int(sqrt(qpvec(i,1)**2+qpvec(i,2)**2)*10.)+1
               if(ip.le.50) then
                  distqp(k,ip)=distqp(k,ip)+1
               endif
            enddo
         enddo
         write(6,*) '# of successes:',nsamplesucc
         do i=1,50
            distqp(k,i)=distqp(k,i)/((i-0.5)*0.1*nsamplesucc)
         enddo
      enddo
      open(9,file='dist_qperp_ferm.res',status='new')
      write(9,*) m
      write(9,*) temp
      do i=1,50
         write(9,10) (i-0.5)*0.1,(distqp(k,i),k=1,9)
      enddo
      close(9)
      ifferm=.false.
      nsample=5000000
      write(6,*) 'for gluons:'
      pin(1)=0
      pin(2)=0
      pin(3)=10
      do j=1,3
         uvec(j)=-pin(j)/m
      enddo
      do k=1,9
         nsamplesucc=0
         do i=1,50
            distqp(k,i)=0
         enddo
         do j=1,nsample
            call collrestrunningradbasic(1,ifferm,m,uvec,temp,mu,
     &           tbnit(k),
     &              nbcol,qpvec)
            do i=1,nbcol
               nsamplesucc=nsamplesucc+1
               ip=int(sqrt(qpvec(i,1)**2+qpvec(i,2)**2)*10.)+1
               if(ip.le.50) then
                  distqp(k,ip)=distqp(k,ip)+1
               endif
            enddo
         enddo
         write(6,*) '# of successes:',nsamplesucc
         do i=1,50
            distqp(k,i)=distqp(k,i)/((i-0.5)*0.1*nsamplesucc)
         enddo
      enddo
      open(9,file='dist_qperp_gluo.res',status='new')
      write(9,*) m
      write(9,*) temp
      do i=1,50
         write(9,10) (i-0.5)*0.1,(distqp(k,i),k=1,9)
      enddo
      close(9)
      return
      end


      subroutine gennewPradiat(x,s,pplus,mpt2,mg,kperp,qperp,
     &      ppl,ifail)
      implicit none
      integer ifail
      double precision pplus,x,b,bplus,bminus,s,mpt2,mpt,mgt2,qperp,
     & qperp2,yminus,qprimeminus,ppplus,ppminus,ppl,mg,kperp
      logical ifgluonmassinboundary
      parameter(ifgluonmassinboundary=.false.)
      if(ifgluonmassinboundary) then
         mgt2=mg**2+kperp**2
      else
         mgt2=kperp**2
      endif
      mpt=sqrt(mpt2)
      b=(1-x)*(s-mgt2/x)
      ifail=0
      bplus=b-(mpt+qperp)**2
      if(bplus.lt.0) then
         ifail=2
         return
      endif
      bminus=b-(mpt-qperp)**2
      yminus=(bplus+bminus+2*sqrt(bplus*bminus))/4
      qperp2=qperp**2
      qprimeminus=(yminus+qperp2)/((1-x)*pplus)
      ppplus=(1-x)*pplus-qperp2/qprimeminus
      ppminus=mpt2/ppplus
      ppl=(ppplus-ppminus)/2
      if(ppl.lt.0d0) then
         ifail=1
      endif
      return
      end
      

      subroutine onecollrunningV2(itypq,ifferm,m,p,temp,mu,ifrad,
     & nbcollsucc,
     &  pp,ppt2)
      implicit none
      logical ifferm,ifrad
      integer i,j,itypq,nbcol,nbcolmax,nit,nbcollsucc,ifail,ifail2
      parameter (nit=20)
      parameter(nbcolmax=10)
      double precision m,p(0:3),temp,mu,pp(0:3),u4vec(0:3),uvec(3),
     & pprest(0:3),
     &  qp(nbcolmax,3),qperp2,qperp,pptx,ppty,ppt2,x,kperp,cosphi,
     & sinphi,pplus,mg,
     & ppl,pplnew,al,atx,aty,vecet(3),t,mpt2,s,pl
      integer ifmtx
      call getMonitorFileIndex(ifmtx)
      if(itypq.ge.3) then
         write(ifmtx,*) 'unrecognized itypquark in onecollrunning:', 
     &        itypq
         close(ifmtx)
         stop
      endif
      nbcollsucc=0
      u4vec(0)=p(0)/m
      do i=1,3
         u4vec(i)=-p(i)/m
         uvec(i)=u4vec(i)
      enddo
      pl=sqrt(p(1)**2+p(2)**2+p(3)**2)
      if(.not.(ifrad)) then
         nbcollsucc=1
         call onecollrestrunning(itypq,ifferm,m,uvec,temp,mu,
     &           pprest,.false.)
         call lorentz2(u4vec,pprest,pp)
      else
         call collrestrunningrad(itypq,ifferm,m,uvec,temp,mu,nit,nbcol,
     &      qp,s,t)
         do i=0,3
            pp(i)=p(i)
         enddo
         if(nbcol.eq.0) then
            return
         else
            if(nbcol.ge.2) then
              write(ifmtx,*)'need sever. successive coll;onecollrunning'
               close(ifmtx)
               stop
            endif
            do i=1,nbcol
               qperp2=qp(i,1)**2+qp(i,2)**2+qp(i,3)**2
               call get1x_veck_qQ_radiat(itypq,pl,temp,m,qperp2,x,kperp,
     &           cosphi,sinphi,ifail)
               if(ifail.eq.0) then
                  ppl=sqrt(pp(1)**2+pp(2)**2+pp(3)**2)
                  pplus=pp(0)+ppl
                  pptx=sqrt(qperp2)-kperp*cosphi
                  ppty=kperp*sinphi
                  ppt2=pptx**2+ppty**2
                  mpt2=m**2+ppt2
                  qperp=sqrt(qperp2)
                  mg=0d0
                  call gennewPradiat(x,s,pplus,mpt2,mg,kperp,
     &                     qperp,pplnew,ifail2)
                  if(ifail2.le.1) then
                     nbcollsucc=nbcollsucc+1
                     al=pplnew/ppl
                     atx=pptx/qperp
                     aty=ppty/(qperp*ppl)
                     vecet(1)=pp(2)*qp(i,3)-qp(i,2)*pp(3)
                     vecet(2)=pp(3)*qp(i,1)-qp(i,3)*pp(1)
                     vecet(3)=pp(1)*qp(i,2)-qp(i,1)*pp(2)
                     pp(0)=m**2
                     do j=1,3
C BUG found in nov 2013 !!!!!!
                        pp(j)=al*pp(j)+atx*qp(i,j)+aty*vecet(j)
                        pp(0)=pp(0)+pp(j)**2
                     enddo
                     pp(0)=sqrt(pp(0))
                  endif
               endif
            enddo
         endif
      endif
      if(pp(0).gt.(p(0)+2*temp)) then
         pp(0)=p(0)
         pp(1)=p(1)
         pp(2)=p(2)
         pp(3)=p(3)
      endif
      return
      end


      subroutine onecollrunning(itypq,ifferm,m,p,temp,mu,ifrad,
     & nbcollsucc,pp)
      use JA
      implicit none
      logical ifferm,ifrad
      integer i,itypq,nbcol,nbcolmax,nit,nbcollsucc
      parameter (nit=20)
      parameter(nbcolmax=10)
      double precision m,p(0:3),temp,mu,pp(0:3),u4vec(0:3),uvec(3),
     & pprest(0:3),qperp2,s,pl,u
      integer ifmtx
      call getMonitorFileIndex(ifmtx)
      if(itypq.ge.3) then
         write(ifmtx,*) 'unrecognized itypquark in onecollrunning:', 
     &        itypq
         close(ifmtx)
         stop
      endif
         if(ifferm)then
         jaityp=1
         else
         jaityp=2
         endif

      nbcollsucc=0
      do i=0,3
         pp(i)=p(i)
      enddo            
      u4vec(0)=p(0)/m
      do i=1,3
         u4vec(i)=-p(i)/m
         uvec(i)=u4vec(i)
      enddo
      pl=sqrt(p(1)**2+p(2)**2+p(3)**2)
      u=pl/m
      if(isnan(u)) then
         write(ifmtx,*) 'isnan u in onecollrunning'
         close(ifmtx)
         stop
      endif
      if(isnan(u/u)) then
         return
      endif
      if(.not.(ifrad)) then
         nbcollsucc=1
         call onecollrestrunning(itypq,ifferm,m,uvec,temp,mu,
     &     pprest,.false.)
         call lorentz2(u4vec,pprest,pp)
      else
         call collrestrunningradnew(itypq,ifferm,m,pl,uvec,temp,mu,
     &     nbcol,
     &     pprest,s,qperp2)
         if(nbcol.eq.0) then
            return
         else 
            if(nbcol.ge.2) then
             write(ifmtx,*)'need several successive coll;onecollrunning'
               close(ifmtx)
               stop
            else
               nbcollsucc=1
               call lorentz2(u4vec,pprest,pp)
            endif
         endif
      endif
      return
      end


      subroutine oneaveragecollrunning(itypq,m,p,temp,mu,ifrad,
     &  nbcollsucc,pp)
      implicit none
      integer itypq,nbcollsucc
      logical ifrad,ifferm
      double precision m,p(0:3),temp,mu,pnorm,pp(0:3),ratef,rateg,
     & ratetot,ran2
      integer ifmtx
      call getMonitorFileIndex(ifmtx)
      pnorm=sqrt(p(1)**2+p(2)**2+p(3)**2)
      call gimmerates(itypq,pnorm,temp,ratef,rateg)
      ratetot=ratef+rateg
      if(ratetot.eq.0) then
         write(ifmtx,*) '0 rate in oneaveragecollrunning'
         close(ifmtx)
         stop
      endif
      ifferm=(ran2()*ratetot.lt.ratef)
      call onecollrunning(itypq,ifferm,m,p,temp,mu,ifrad,nbcollsucc,pp)
      return
      end


      subroutine testonecollrunning
      implicit none
      logical ifferm,ifrad
      integer i,j,nsample,nbcol,nbcoltot,itypq,ipt,ipz
      double precision pi,m,tbpin(0:30),pin(0:3),mu,temp,pp(0:3),
     & meandp,
     & momlong,momtrans,mdebeffinterp,nbsucces,ratef,rateg,ratetot,ppt2 
      parameter (pi=3.14159265358979323844d0)
      double precision distribpt(1:100,2)
      double precision distribpz(-100:100,2)
 10   format(F10.5,' ', 12 (F12.7,' '))
      write(6,*) 'entering test onecol'
      m=1.5d0
      itypq=1
      temp=0.3d0
      mu=0.4472136d0*mdebeffinterp(temp)      
      tbpin(0)=0.1D0
      do i=1,4
         tbpin(i)=i
      enddo
      do i=5,16
         tbpin(i)=5*(i-4)
      enddo
      do i=17,30
         tbpin(i)=10*(i-16)+50
      enddo

      goto 6
      write(6,*) 'generating dpt and dpl distributions'
      pin(1)=0
      pin(2)=0
      pin(3)=25
      ifrad=.false.
      write(6,*) 'for collisional'


      write(6,*) 'for collisional with fermions and gluons'
      m=1.5d0
      itypq=1
      pin(0)=sqrt(pin(3)**2+m**2)
      do i=1,100
         distribpt(i,1)=0
         distribpt(i,2)=0
      enddo
      do i=-100,100
         distribpz(i,1)=0
         distribpz(i,2)=0
      enddo
      nbcoltot=0
      nsample=10000000
      do j=1,nsample
         call oneaveragecollrunning(itypq,m,pin,temp,mu,ifrad,nbcol,
     &        pp)         
         if(nbcol.gt.0) then
            if(nbcol.gt.1) then
               write(6,*) 'problem'
               stop
            endif
            nbcoltot=nbcoltot+1
             ppt2=sqrt(pp(1)**2+pp(2)**2)
            if((ppt2.le.10).and.(ppt2.gt.0)) then
               ipt=int(ppt2*10)+1
               distribpt(ipt,1)=distribpt(ipt,1)+1.
            endif
            ipz=nint(pp(3)/(pin(3)+10)*100)
            if(abs(ipz).le.100) then
                distribpz(ipz,1)=distribpz(ipz,1)+1.
            endif
         endif
      enddo
      do i=1,100
         distribpt(i,1)=distribpt(i,1)/nbcoltot
      enddo
      do i=-100,100
         distribpz(i,1)=distribpz(i,1)/nbcoltot
      enddo
      itypq=2
      m=5.1d0
      pin(0)=sqrt(pin(3)**2+m**2)
      nbcoltot=0
      do j=1,nsample
         call oneaveragecollrunning(itypq,m,pin,temp,mu,ifrad,nbcol,
     &        pp)         
         if(nbcol.gt.0) then
            nbcoltot=nbcoltot+1
            ppt2=sqrt(pp(1)**2+pp(2)**2)
            if((ppt2.le.10).and.(ppt2.gt.0)) then
               ipt=int(ppt2*10)+1
               distribpt(ipt,2)=distribpt(ipt,2)+1.
            endif
            ipz=nint(pp(3)/(pin(3)+10)*100)
            if(abs(ipz).le.100) then
                distribpz(ipz,2)=distribpz(ipz,2)+1.
            endif
         endif
      enddo
      do i=1,100
         distribpt(i,2)=distribpt(i,2)/nbcoltot
      enddo
      do i=-100,100
         distribpz(i,2)=distribpz(i,2)/nbcoltot
      enddo
      open(9,file='dist_Pprime_perp_av_elast.res',status='new')
      write(9,*) temp,mu,pin(3),nsample,nbcoltot
      write(9,*) 'efficiency:', nbcoltot/(1.*nsample)
      write(9,*) 'Pperp | dN_c | dN_b'
      do i=1,100
         write(9,10) (i-0.5)/10,distribpt(i,1)*10,distribpt(i,2)*10
      enddo
      close(9)
      open(9,file='dist_Pprime_long_av_elast.res',status='new')
      write(9,*) temp,mu,pin(3),nsample,nbcoltot
      write(9,*) 'efficiency:', nbcoltot/(1.*nsample)
      write(9,*) 'Pz | dN_c | dN_b'
      do i=-100,100
         write(9,10) i*(pin(3)+10.)/100,distribpz(i,1)*100/pin(3),
     &      distribpz(i,2)*100/pin(3)
      enddo
      close(9)

      ifrad=.true.
      itypq=1
      m=1.5d0
      pin(0)=sqrt(pin(3)**2+m**2)
      write(6,*) 'for radiative'


      write(6,*) 'for radiative with fermions and gluons'
      m=1.5d0
      itypq=1
      pin(0)=sqrt(pin(3)**2+m**2)
      do i=1,100
         distribpt(i,1)=0
         distribpt(i,2)=0
      enddo
      do i=-100,100
         distribpz(i,1)=0
         distribpz(i,2)=0
      enddo
      nbcoltot=0
      nsample=10000000
      do j=1,nsample
         call oneaveragecollrunning(itypq,m,pin,temp,mu,ifrad,nbcol,
     &        pp)         
         if(nbcol.gt.0) then
            if(nbcol.gt.1) then
               write(6,*) 'problem'
               stop
            endif
            nbcoltot=nbcoltot+1
             ppt2=sqrt(pp(1)**2+pp(2)**2)
            if((ppt2.le.10).and.(ppt2.gt.0)) then
               ipt=int(ppt2*10)+1
               distribpt(ipt,1)=distribpt(ipt,1)+1.
            endif
            ipz=nint(pp(3)/(pin(3)+10)*100)
            if(abs(ipz).le.100) then
                distribpz(ipz,1)=distribpz(ipz,1)+1.
            endif
         endif
      enddo
      do i=1,100
         distribpt(i,1)=distribpt(i,1)/nbcoltot
      enddo
      do i=-100,100
         distribpz(i,1)=distribpz(i,1)/nbcoltot
      enddo
      itypq=2
      m=5.1d0
      pin(0)=sqrt(pin(3)**2+m**2)
      nbcoltot=0
      do j=1,nsample
         call oneaveragecollrunning(itypq,m,pin,temp,mu,ifrad,nbcol,
     &        pp)         
         if(nbcol.gt.0) then
            nbcoltot=nbcoltot+1
            ppt2=sqrt(pp(1)**2+pp(2)**2)
            if((ppt2.le.10).and.(ppt2.gt.0)) then
               ipt=int(ppt2*10)+1
               distribpt(ipt,2)=distribpt(ipt,2)+1.
            endif
            ipz=nint(pp(3)/(pin(3)+10)*100)
            if(abs(ipz).le.100) then
                distribpz(ipz,2)=distribpz(ipz,2)+1.
            endif
         endif
      enddo
      do i=1,100
         distribpt(i,2)=distribpt(i,2)/nbcoltot
      enddo
      do i=-100,100
         distribpz(i,2)=distribpz(i,2)/nbcoltot
      enddo
      open(9,file='dist_Pprime_perp_av_rad.res',status='new')
      write(9,*) temp,mu,pin(3),nsample,nbcoltot
      write(9,*) 'efficiency:', nbcoltot/(1.*nsample)
      write(9,*) 'Pperp | dN_c | dN_b'
      do i=1,100
         write(9,10) (i-0.5)/10,distribpt(i,1)*10,distribpt(i,2)*10
      enddo
      close(9)
      open(9,file='dist_Pprime_long_av_rad.res',status='new')
      write(9,*) temp,mu,pin(3),nsample,nbcoltot
      write(9,*) 'efficiency:', nbcoltot/(1.*nsample)
      write(9,*) 'Pz | dN_c | dN_b'
      do i=-100,100
         write(9,10) i*(pin(3)+10.)/100,distribpz(i,1)*100/pin(3),
     &      distribpz(i,2)*100/pin(3)
      enddo
      close(9)

      write(6,*) 'for ell + radiative with fermions + gluons'
      m=1.5d0
      itypq=1
      pin(0)=sqrt(pin(3)**2+m**2)
      do i=1,100
         distribpt(i,1)=0
         distribpt(i,2)=0
      enddo
      do i=-100,100
         distribpz(i,1)=0
         distribpz(i,2)=0
      enddo
      nbcoltot=0
      nsample=10000000
      do j=1,nsample
         call oneaveragecollrunning(itypq,m,pin,temp,mu,.false.,nbcol,
     &        pp)         
         if(nbcol.gt.0) then
            if(nbcol.gt.1) then
               write(6,*) 'problem'
               stop
            endif
            nbcoltot=nbcoltot+1
             ppt2=sqrt(pp(1)**2+pp(2)**2)
            if((ppt2.le.10).and.(ppt2.gt.0)) then
               ipt=int(ppt2*10)+1
               distribpt(ipt,1)=distribpt(ipt,1)+1.
            endif
            ipz=nint(pp(3)/(pin(3)+10)*100)
            if(abs(ipz).le.100) then
                distribpz(ipz,1)=distribpz(ipz,1)+1.
            endif
         endif
         call oneaveragecollrunning(itypq,m,pin,temp,mu,.true.,nbcol,
     &        pp)         
         if(nbcol.gt.0) then
            if(nbcol.gt.1) then
               write(6,*) 'problem'
               stop
            endif
            nbcoltot=nbcoltot+1
             ppt2=sqrt(pp(1)**2+pp(2)**2)
            if((ppt2.le.10).and.(ppt2.gt.0)) then
               ipt=int(ppt2*10)+1
               distribpt(ipt,1)=distribpt(ipt,1)+1.
            endif
            ipz=nint(pp(3)/(pin(3)+10)*100)
            if(abs(ipz).le.100) then
                distribpz(ipz,1)=distribpz(ipz,1)+1.
            endif
         endif
      enddo
      do i=1,100
         distribpt(i,1)=distribpt(i,1)/nbcoltot
      enddo
      do i=-100,100
         distribpz(i,1)=distribpz(i,1)/nbcoltot
      enddo
      itypq=2
      m=5.1d0
      pin(0)=sqrt(pin(3)**2+m**2)
      nbcoltot=0
      do j=1,nsample
         call oneaveragecollrunning(itypq,m,pin,temp,mu,.false.,nbcol,
     &        pp)         
         if(nbcol.gt.0) then
            nbcoltot=nbcoltot+1
            ppt2=sqrt(pp(1)**2+pp(2)**2)
            if((ppt2.le.10).and.(ppt2.gt.0)) then
               ipt=int(ppt2*10)+1
               distribpt(ipt,2)=distribpt(ipt,2)+1.
            endif
            ipz=nint(pp(3)/(pin(3)+10)*100)
            if(abs(ipz).le.100) then
                distribpz(ipz,2)=distribpz(ipz,2)+1.
            endif
         endif
         call oneaveragecollrunning(itypq,m,pin,temp,mu,.true.,nbcol,
     &        pp)         
         if(nbcol.gt.0) then
            nbcoltot=nbcoltot+1
            ppt2=sqrt(pp(1)**2+pp(2)**2)
            if((ppt2.le.10).and.(ppt2.gt.0)) then
               ipt=int(ppt2*10)+1
               distribpt(ipt,2)=distribpt(ipt,2)+1.
            endif
            ipz=nint(pp(3)/(pin(3)+10)*100)
            if(abs(ipz).le.100) then
                distribpz(ipz,2)=distribpz(ipz,2)+1.
            endif
         endif
      enddo
      do i=1,100
         distribpt(i,2)=distribpt(i,2)/nbcoltot
      enddo
      do i=-100,100
         distribpz(i,2)=distribpz(i,2)/nbcoltot
      enddo
      open(9,file='dist_Pprime_perp_av_mixed.res',status='new')
      write(9,*) temp,mu,pin(3),nsample,nbcoltot
      write(9,*) 'efficiency:', nbcoltot/(1.*nsample)
      write(9,*) 'Pperp | dN_c | dN_b'
      do i=1,100
         write(9,10) (i-0.5)/10,distribpt(i,1)*10,distribpt(i,2)*10
      enddo
      close(9)
      open(9,file='dist_Pprime_long_av_mixed.res',status='new')
      write(9,*) temp,mu,pin(3),nsample,nbcoltot
      write(9,*) 'efficiency:', nbcoltot/(1.*nsample)
      write(9,*) 'Pz | dN_c | dN_b'
      do i=-100,100
         write(9,10) i*(pin(3)+10.)/100,distribpz(i,1)*100/pin(3),
     &      distribpz(i,2)*100/pin(3)
      enddo
      close(9)


      write(6,*) 'generating transport coefficients'
      itypq=1
      m=1.5d0
      ifrad=.false.
      write(6,*) 'for elastic processes'      
      goto 5
      ifferm=.true.
      nsample=10000000
      write(6,*) 'for fermions:'
      open(9,file='av_moment_varia_el_ferm.res',status='new')
      write(9,*) m,temp,mu,nsample
      write(9,*) 'p |rate|dp_z/coll|dp_z^2/coll|dp_t^2/coll|dp_z/dt|
     &  dp_t^2/dt'      
      do i=0,30
         pin(1)=0
         pin(2)=0
         pin(3)=tbpin(i)
         pin(0)=sqrt(pin(3)**2+m**2)
         meandp=0
         momlong=0.
         momtrans=0.
         call gimmerates(itypq,pin(3),temp,ratef,rateg)
         do j=1,nsample
            call onecollrunning(itypq,ifferm,m,pin,temp,mu,ifrad,nbcol,
     &           pp) 
            meandp=meandp+pp(3)
            momlong=momlong+(pin(3)-pp(3))**2
            momtrans=momtrans+(pp(1)**2+pp(2)**2)
         enddo
         meandp=pin(3)-meandp/nsample
         momlong=momlong/nsample
         momtrans=momtrans/nsample
         write(6,10) pin(3),ratef,meandp,momlong,momtrans,ratef*meandp,
     &                ratef*momtrans
         write(9,10) pin(3),ratef,meandp,momlong,momtrans,ratef*meandp,
     &                ratef*momtrans
      enddo
      close(9)
      ifferm=.false.
      nsample=10000000
      write(6,*) 'for gluons:'
      open(9,file='av_moment_varia_el_gluo.res',status='new')
      write(9,*) m,temp,mu,nsample
      write(9,*) 'p |rate|dp_z/coll|dp_z^2/coll|dp_t^2/coll|dp_z/dt| 
     & dp_t^2/dt'        
      do i=0,30
         pin(1)=0
         pin(2)=0
         pin(3)=tbpin(i)
         pin(0)=sqrt(pin(3)**2+m**2)
         meandp=0
         momlong=0.
         momtrans=0.
         call gimmerates(itypq,pin(3),temp,ratef,rateg)
         do j=1,nsample
            call onecollrunning(itypq,ifferm,m,pin,temp,mu,ifrad,nbcol,
     &         pp)  
            meandp=meandp+pp(3)
            momlong=momlong+(pin(3)-pp(3))**2
            momtrans=momtrans+(pp(1)**2+pp(2)**2)
         enddo
         meandp=pin(3)-meandp/nsample
         momlong=momlong/nsample
         momtrans=momtrans/nsample
         write(6,10) pin(3),rateg,meandp,momlong,momtrans,rateg*meandp,
     &                rateg*momtrans
         write(9,10) pin(3),rateg,meandp,momlong,momtrans,rateg*meandp,
     &                rateg*momtrans
      enddo
      close(9)

 5    continue
      write(6,*) 'for both scattering on quarks and gluons'
      nsample=10000000
      open(9,file='av_moment_varia_el.res',status='new')
      write(9,*) m,temp,mu,nsample
      write(9,*) 'p|rate|dp_z/coll|dp_z^2/coll|dp_t^2/coll|dp_z/dt|
     &   dp_z^2/dt | dp_t^2/dt'  
      do i=0,30
         pin(1)=0
         pin(2)=0
         pin(3)=tbpin(i)
         pin(0)=sqrt(pin(3)**2+m**2)
         meandp=0
         momlong=0.
         momtrans=0.
         call gimmerates(itypq,pin(3),temp,ratef,rateg)
         do j=1,nsample
            call oneaveragecollrunning(itypq,m,pin,temp,mu,ifrad,nbcol,
     &        pp)         
            meandp=meandp+pp(3)
            momlong=momlong+(pin(3)-pp(3))**2
            momtrans=momtrans+(pp(1)**2+pp(2)**2)
         enddo
         meandp=pin(3)-meandp/nsample
         momlong=momlong/nsample
         momtrans=momtrans/nsample
         ratetot=ratef+rateg
         write(6,10) pin(3),ratetot,meandp,momlong,momtrans,
     &       meandp*ratetot,momlong*ratetot,momtrans*ratetot
         write(9,10) pin(3),ratetot,meandp,momlong,momtrans,
     &       meandp*ratetot,momlong*ratetot,momtrans*ratetot
      enddo
      close(9)

 6    continue
      ifrad=.true.
      write(6,*) 'for radiative processes'      
      ifferm=.true.
      nsample=10000000
      write(6,*) 'for fermions:'
      open(9,file='av_moment_varia_rad_ferm.res',status='new')
      write(9,*) m,temp,mu,nsample
      write(9,*) 'p |rate|dp_z/coll|dp_z^2/coll|dp_t^2/coll|dp_z/dt|
     &  dp_t^2/dt'  
      do i=0,30
         pin(1)=0
         pin(2)=0
         pin(3)=tbpin(i)
         pin(0)=sqrt(pin(3)**2+m**2)
         meandp=0
         momlong=0.
         momtrans=0.
         nbsucces=0
         call gimmerates(itypq,pin(3),temp,ratef,rateg)
         do j=1,nsample
            call onecollrunning(itypq,ifferm,m,pin,temp,mu,ifrad,nbcol,
     &         pp) 
            if(nbcol.eq.1) then
               nbsucces=nbsucces+1
               meandp=meandp+pp(3)
               momlong=momlong+(pin(3)-pp(3))**2
               momtrans=momtrans+(pp(1)**2+pp(2)**2)
            endif
         enddo
         meandp=pin(3)-meandp/nbsucces
         momlong=momlong/nbsucces
         momtrans=momtrans/nbsucces
         nbsucces=nbsucces/nsample
         ratef=nbsucces*ratef
         write(6,10) pin(3),ratef,meandp,momlong,momtrans,
     &      ratef*meandp,ratef*momtrans
         write(9,10) pin(3),ratef,meandp,momlong,momtrans,
     &      ratef*meandp,ratef*momtrans
      enddo
      close(9)
      ifferm=.false.
      nsample=10000000
      write(6,*) 'for gluons:'
      open(9,file='av_moment_varia_rad_gluo.res',status='new')
      write(9,*) m,temp,mu,nsample
      write(9,*) 'p |rate|dp_z/coll|dp_z^2/coll|dp_t^2/coll|dp_z/dt| 
     & dp_t^2/dt'    
      do i=0,30
         pin(1)=0
         pin(2)=0
         pin(3)=tbpin(i)
         pin(0)=sqrt(pin(3)**2+m**2)
         meandp=0
         momlong=0.
         momtrans=0.
         nbsucces=0
         call gimmerates(itypq,pin(3),temp,ratef,rateg)
         do j=1,nsample
            call onecollrunning(itypq,ifferm,m,pin,temp,mu,ifrad,nbcol,
     &         pp)  
            if(nbcol.eq.1) then
               nbsucces=nbsucces+1
               meandp=meandp+pp(3)
               momlong=momlong+(pin(3)-pp(3))**2
               momtrans=momtrans+(pp(1)**2+pp(2)**2)
            endif
         enddo
         meandp=pin(3)-meandp/nbsucces
         momlong=momlong/nbsucces
         momtrans=momtrans/nbsucces
         nbsucces=nbsucces/nsample
         rateg=nbsucces*rateg
         write(6,10) pin(3),rateg,meandp,momlong,momtrans,
     &      rateg*meandp,rateg*momtrans
         write(9,10) pin(3),rateg,meandp,momlong,momtrans,
     &      rateg*meandp,rateg*momtrans
      enddo
      close(9)     
  
      nsample=10000000
      write(6,*) 'for both scattering on quarks and gluons'
      open(9,file='av_moment_varia_rad.res',status='new')
      write(9,*) m,temp,mu,nsample
      write(9,*) 'p |rate|dp_z/coll|dp_z^2/coll|dp_t^2/coll|dp_z/dt|
     &   dp_z^2/dt | dp_t^2/dt'  
      do i=0,30
         pin(1)=0
         pin(2)=0
         pin(3)=tbpin(i)
         pin(0)=sqrt(pin(3)**2+m**2)
         meandp=0
         momlong=0.
         momtrans=0.
         nbsucces=0
         call gimmerates(itypq,pin(3),temp,ratef,rateg)
         do j=1,nsample
            call oneaveragecollrunning(itypq,m,pin,temp,mu,ifrad,nbcol,
     &        pp)         
            if(nbcol.eq.1) then
               nbsucces=nbsucces+1
               meandp=meandp+pp(3)
               momlong=momlong+(pin(3)-pp(3))**2
               momtrans=momtrans+(pp(1)**2+pp(2)**2)
            endif
         enddo
         meandp=pin(3)-meandp/nbsucces
         momlong=momlong/nbsucces
         momtrans=momtrans/nbsucces
         nbsucces=nbsucces/nsample
         ratetot=nbsucces*(rateg+ratef)
         write(6,10) pin(3),ratetot,meandp,momlong,momtrans,
     &   ratetot*meandp,
     &   ratetot*momlong,ratetot*momtrans
         write(9,10) pin(3),ratetot,meandp,momlong,momtrans,
     &   ratetot*meandp,
     &   ratetot*momlong,ratetot*momtrans
      enddo
      close(9)     

      write(6,*) 'for el + radiative processes'      
      nsample=10000000
      write(6,*) 'for both scattering on quarks and gluons'
      open(9,file='av_moment_varia_mixed.res',status='new')
      write(9,*) m,temp,mu,nsample
      write(9,*) 'p |rate|dp_z/coll|dp_z^2/coll|dp_t^2/coll|dp_z/dt|
     &   dp_z^2/dt | dp_t^2/dt'  
      do i=0,30
         pin(1)=0
         pin(2)=0
         pin(3)=tbpin(i)
         pin(0)=sqrt(pin(3)**2+m**2)
         meandp=0
         momlong=0.
         momtrans=0.
         nbsucces=0
         call gimmerates(itypq,pin(3),temp,ratef,rateg)
         do j=1,nsample
            call oneaveragecollrunning(itypq,m,pin,temp,mu,.false.,
     &        nbcol,
     &        pp)         
            if(nbcol.eq.1) then
               nbsucces=nbsucces+1
               meandp=meandp+pp(3)
               momlong=momlong+(pin(3)-pp(3))**2
               momtrans=momtrans+(pp(1)**2+pp(2)**2)
            else
               write(6,*) 'should always be success for coll'
               stop
            endif
            call oneaveragecollrunning(itypq,m,pin,temp,mu,.true.,
     &       nbcol,
     &        pp)         
            if(nbcol.eq.1) then
               nbsucces=nbsucces+1
               meandp=meandp+pp(3)
               momlong=momlong+(pin(3)-pp(3))**2
               momtrans=momtrans+(pp(1)**2+pp(2)**2)
            endif
         enddo
         meandp=pin(3)-meandp/nbsucces
         momlong=momlong/nbsucces
         momtrans=momtrans/nbsucces
         nbsucces=nbsucces/nsample
         ratetot=nbsucces*(rateg+ratef)
         write(6,10) pin(3),ratetot,meandp,momlong,momtrans,
     &   ratetot*meandp,
     &   ratetot*momlong,ratetot*momtrans
         write(9,10) pin(3),ratetot,meandp,momlong,momtrans,
     &   ratetot*meandp,
     &   ratetot*momlong,ratetot*momtrans
      enddo
      close(9)     

      return
      end

      subroutine evolveboltzmann(itypq,pin,temp,effdeg,tevol,ifcol,
     & ifrad,nbcoltot,ifinfo)

      use JA
      use PBGpsiinfo
      use PBGupsinfo
      use PBGforboltzmann
      use PBGreductiondofs
      implicit none
      logical first,ifferm,ifcol,ifrad,ifinfo
      integer i,itypq,nbcol,nbcoltot
      double precision m,temp,effdeg,tevol,tref,pin(0:3),mu,muoftemp,
     &  pint(0:3),pnorm,ratef,rateg,ratetot,ran2,mdebeffinterp,x,s
      double precision reduction_dof
      external reduction_dof
      integer ifmtx
      call getMonitorFileIndex(ifmtx)
      first=.true.
      nbcoltot=0
      if(itypq.eq.1) then
         m=mc
      else if(itypq.eq.2) then
         m=mb
      else
         write(ifmtx,*) 'unrecognized itypquark in evolveboltzmann:', 
     &        itypq
         close(ifmtx)
         stop
      endif
      tref=0
      if(ifinfo) write(6,*) 'entering evolveboltz with tevol:',tevol
      pnorm=sqrt(pin(1)**2+pin(2)**2+pin(3)**2)
      do while((tref.lt.tevol).or.(pnorm.gt.1D8))
         call gimmerates(itypq,pnorm,temp,ratef,rateg)
         if(reddof==3)then
            effdeg=reduction_dof(1,pnorm,temp)
         endif
         ratef=ratef*effdeg
         if(reddof==3)then
            effdeg=reduction_dof(0,pnorm,temp)
         endif
         rateg=rateg*effdeg 
         ratetot=ratef+rateg
         if(ifinfo) write(6,*) 'tref:',tref,',pnorm:',pin(0),pnorm,
     &   'rate:',ratetot
         if(ratetot.eq.0) then
            tref=2*tevol
         else
            tref=tref-log(ran2())/ratetot
         endif
         if(tref.le.tevol) then
            if(first) then
               if(boltz_model.eq.1) then
                  mu=muoftemp(temp)
               else if(boltz_model.eq.2) then
                  mu=mdebeffinterp(temp)
               else if(boltz_model.eq.3) then
                  mu=0.387299*muoftemp(temp)
               else if(boltz_model.eq.4) then
                  mu=0.4472136d0*mdebeffinterp(temp)
               else
                  write(ifmtx,*) 'unrecogn. model in evolveboltzmann:'
     &                    ,boltz_model
                  close(ifmtx)
                  stop
               endif                 
               first=.false.
            endif
            ifferm=(ran2()*ratetot.lt.ratef)
            if(ifcol) then
               if((boltz_model.eq.1).or.(boltz_model.eq.3)) then
                  call onecoll(itypq,ifferm,m,pin,temp,mu,.false.,nbcol,
     &                pint,
     &               x,s)
               else if((boltz_model.eq.2).or.(boltz_model.eq.4)) then
                  call onecollrunning(itypq,ifferm,m,pin,temp,
     &                 mu,.false.,
     &                 nbcol,pint)
               endif
               do i=0,3
                 pin(i)=pint(i)
               enddo
               nbcoltot=nbcoltot+1
            endif
            if(ifrad) then
               if((boltz_model.eq.1).or.(boltz_model.eq.3)) then
                  call onecoll(itypq,ifferm,m,pin,temp,mu,.true.,
     &                nbcol,pint,
     &                x,s)
               else if((boltz_model.eq.2).or.(boltz_model.eq.4)) then
                  call onecollrunning(itypq,ifferm,m,pin,temp,mu,.true.,
     &                 nbcol,pint)

               endif
               do i=0,3
                 pin(i)=pint(i)
               enddo
               nbcoltot=nbcoltot+nbcol
            endif
            pnorm=sqrt(pin(1)**2+pin(2)**2+pin(3)**2)
            pin(0)=sqrt(pnorm**2+m**2)
         endif
      enddo
      return
      end


      subroutine testboltzmann
      implicit none
      logical ifcol,ifrad
      integer i,j,k,nsample,nbcol,itypq,ipt,ipz
      double precision pi,m,tbpin(20),pin(0:3),mu,temp,pp(0:3),meandp,
     & momlong,momtrans,mdebeffinterp,ratetot,ppt2,
     & tevol,effdeg,nbcoltot(2) 
      parameter (pi=3.14159265358979323844d0)
      double precision distribpt(1:100,2)
      double precision distribpz(-100:100,2)
 10   format(F10.5,' ', 12 (F12.7,' '))
      write(6,*) 'entering test boltzman'
      itypq=1
      temp=0.4d0
      mu=0.4472136d0*mdebeffinterp(temp)      
      effdeg=1.d0
      tevol=0.4
      do i=1,4
         tbpin(i)=i
      enddo
      do i=5,16
         tbpin(i)=5*(i-4)
      enddo

      write(6,*) 'generating dpt and dpl distributions'
      pin(1)=0
      pin(2)=0
      pin(3)=25
      write(6,*) 'for collisional'
      ifcol=.true.
      ifrad=.false.
      m=1.5d0
      itypq=1
      pin(0)=sqrt(pin(3)**2+m**2)
      do i=1,100
         distribpt(i,1)=0
         distribpt(i,2)=0
      enddo
      do i=-100,100
         distribpz(i,1)=0
         distribpz(i,2)=0
      enddo
      nbcoltot(1)=0
      nbcoltot(2)=0
      nsample=1000000
      do j=1,nsample
         do k=0,3
            pp(k)=pin(k)
         enddo
         call evolveboltzmann(itypq,pp,temp,effdeg,tevol,ifcol,ifrad,
     &        nbcol,.false.)
         nbcoltot(1)=nbcoltot(1)+nbcol
         ppt2=sqrt(pp(1)**2+pp(2)**2)
         if((ppt2.le.10).and.(ppt2.ge.0)) then
            ipt=int(ppt2*10)+1
            distribpt(ipt,1)=distribpt(ipt,1)+1.
         endif
         ipz=nint(pp(3)/(pin(3)+10)*100)
         if(abs(ipz).le.100) then
            distribpz(ipz,1)=distribpz(ipz,1)+1.
         endif
      enddo
      do i=1,100
         distribpt(i,1)=distribpt(i,1)/(nsample*tevol)
      enddo
      do i=-100,100
         distribpz(i,1)=distribpz(i,1)/(nsample*tevol)
      enddo
      itypq=2
      m=5.1d0
      pin(0)=sqrt(pin(3)**2+m**2)
      nbcoltot(2)=nbcoltot(2)+nbcol
      do j=1,nsample
         do k=0,3
            pp(k)=pin(k)
         enddo
         call evolveboltzmann(itypq,pp,temp,effdeg,tevol,ifcol,ifrad,
     &                nbcol,.false.)
         nbcoltot(2)=nbcoltot(2)+nbcol
         ppt2=sqrt(pp(1)**2+pp(2)**2)
         if((ppt2.le.10).and.(ppt2.ge.0)) then
            ipt=int(ppt2*10)+1
            distribpt(ipt,2)=distribpt(ipt,2)+1.
         endif
         ipz=nint(pp(3)/(pin(3)+10)*100)
         if(abs(ipz).le.100) then
            distribpz(ipz,2)=distribpz(ipz,2)+1.
         endif
      enddo
      do i=1,100
         distribpt(i,2)=distribpt(i,2)/(nsample*tevol)
      enddo
      do i=-100,100
         distribpz(i,2)=distribpz(i,2)/(nsample*tevol)
      enddo
      open(9,file='dist_Pprime_perp_av_elast_evolveB.res',
     &           status='new')
      write(9,*) temp,mu,pin(3),nsample,nbcoltot(1)/(nsample*tevol),
     & nbcoltot(2)/(nsample*tevol)
      write(9,*) 'Pperp | dN_c | dN_b'
      do i=1,100
         write(9,10) (i-0.5)/10,distribpt(i,1)*10,distribpt(i,2)*10
      enddo
      close(9)
      open(9,file='dist_Pprime_long_av_elast_evolveB.res',status='new')
      write(9,*) temp,mu,pin(3),nsample,nbcoltot(1)/(nsample*tevol),
     & nbcoltot(2)/(nsample*tevol)
      write(9,*) 'Pz | dN_c | dN_b'
      do i=-100,100
         write(9,10) i*(pin(3)+10.)/100,distribpz(i,1)*100/pin(3),
     &      distribpz(i,2)*100/pin(3)
      enddo
      close(9)

      write(6,*) 'for radiative'
      ifcol=.false.
      ifrad=.true.
      m=1.5d0
      itypq=1
      pin(0)=sqrt(pin(3)**2+m**2)
      do i=1,100
         distribpt(i,1)=0
         distribpt(i,2)=0
      enddo
      do i=-100,100
         distribpz(i,1)=0
         distribpz(i,2)=0
      enddo
      nbcoltot(1)=0
      nbcoltot(2)=0
      nsample=1000000
      do j=1,nsample
         do k=0,3
            pp(k)=pin(k)
         enddo
         call evolveboltzmann(itypq,pp,temp,effdeg,tevol,ifcol,ifrad
     &          ,nbcol,.false.)
         nbcoltot(1)=nbcoltot(1)+nbcol
         ppt2=sqrt(pp(1)**2+pp(2)**2)
         if((ppt2.le.10).and.(ppt2.ge.0)) then
            ipt=int(ppt2*10)+1
            distribpt(ipt,1)=distribpt(ipt,1)+1.
         endif
         ipz=nint(pp(3)/(pin(3)+10)*100)
         if(abs(ipz).le.100) then
            distribpz(ipz,1)=distribpz(ipz,1)+1.
         endif
      enddo
      do i=1,100
         distribpt(i,1)=distribpt(i,1)/(nsample*tevol)
      enddo
      do i=-100,100
         distribpz(i,1)=distribpz(i,1)/(nsample*tevol)
      enddo
      itypq=2
      m=5.1d0
      pin(0)=sqrt(pin(3)**2+m**2)
      do j=1,nsample
         do k=0,3
            pp(k)=pin(k)
         enddo
         call evolveboltzmann(itypq,pp,temp,effdeg,tevol,ifcol,ifrad
     &           ,nbcol,.false.)
         nbcoltot(2)=nbcoltot(2)+nbcol
         ppt2=sqrt(pp(1)**2+pp(2)**2)
         if((ppt2.le.10).and.(ppt2.ge.0)) then
            ipt=int(ppt2*10)+1
            distribpt(ipt,2)=distribpt(ipt,2)+1.
         endif
         ipz=nint(pp(3)/(pin(3)+10)*100)
         if(abs(ipz).le.100) then
            distribpz(ipz,2)=distribpz(ipz,2)+1.
         endif
      enddo
      do i=1,100
         distribpt(i,2)=distribpt(i,2)/(nsample*tevol)
      enddo
      do i=-100,100
         distribpz(i,2)=distribpz(i,2)/(nsample*tevol)
      enddo
      open(9,file='dist_Pprime_perp_av_rad_evolveB.res',status='new')
      write(9,*) temp,mu,pin(3),nsample,nbcoltot(1)/(nsample*tevol),
     & nbcoltot(2)/(nsample*tevol)
      write(9,*) 'Pperp | dN_c | dN_b'
      do i=1,100
         write(9,10) (i-0.5)/10,distribpt(i,1)*10,distribpt(i,2)*10
      enddo
      close(9)
      open(9,file='dist_Pprime_long_av_rad_evolveB.res',status='new')
      write(9,*) temp,mu,pin(3),nsample,nbcoltot(1)/(nsample*tevol),
     & nbcoltot(2)/(nsample*tevol)
      write(9,*) 'Pz | dN_c | dN_b'
      do i=-100,100
         write(9,10) i*(pin(3)+10.)/100,distribpz(i,1)*100/pin(3),
     &      distribpz(i,2)*100/pin(3)
      enddo
      close(9)

      write(6,*) 'for collisional + radiative'
      ifcol=.true.
      ifrad=.true.
      m=1.5d0
      itypq=1
      pin(0)=sqrt(pin(3)**2+m**2)
      do i=1,100
         distribpt(i,1)=0
         distribpt(i,2)=0
      enddo
      do i=-100,100
         distribpz(i,1)=0
         distribpz(i,2)=0
      enddo
      nsample=1000000
      nbcoltot(1)=0
      nbcoltot(2)=0
      do j=1,nsample
         do k=0,3
            pp(k)=pin(k)
         enddo
         call evolveboltzmann(itypq,pp,temp,effdeg,tevol,ifcol,ifrad
     &          ,nbcol,.false.)
         nbcoltot(1)=nbcoltot(1)+nbcol
         ppt2=sqrt(pp(1)**2+pp(2)**2)
         if((ppt2.le.10).and.(ppt2.ge.0)) then
            ipt=int(ppt2*10)+1
            distribpt(ipt,1)=distribpt(ipt,1)+1.
         endif
         ipz=nint(pp(3)/(pin(3)+10)*100)
         if(abs(ipz).le.100) then
            distribpz(ipz,1)=distribpz(ipz,1)+1.
         endif
      enddo
      do i=1,100
         distribpt(i,1)=distribpt(i,1)/(nsample*tevol)
      enddo
      do i=-100,100
         distribpz(i,1)=distribpz(i,1)/(nsample*tevol)
      enddo
      itypq=2
      m=5.1d0
      pin(0)=sqrt(pin(3)**2+m**2)
      do j=1,nsample
         do k=0,3
            pp(k)=pin(k)
         enddo
         call evolveboltzmann(itypq,pp,temp,effdeg,tevol,ifcol,ifrad
     &            ,nbcol,.false.)
         nbcoltot(2)=nbcoltot(2)+nbcol
         ppt2=sqrt(pp(1)**2+pp(2)**2)
         if((ppt2.le.10).and.(ppt2.ge.0)) then
            ipt=int(ppt2*10)+1
            distribpt(ipt,2)=distribpt(ipt,2)+1.
         endif
         ipz=nint(pp(3)/(pin(3)+10)*100)
         if(abs(ipz).le.100) then
            distribpz(ipz,2)=distribpz(ipz,2)+1.
         endif
      enddo
      do i=1,100
         distribpt(i,2)=distribpt(i,2)/(nsample*tevol)
      enddo
      do i=-100,100
         distribpz(i,2)=distribpz(i,2)/(nsample*tevol)
      enddo
      open(9,file='dist_Pprime_perp_av_mixed_evolveB.res',status='new')
      write(9,*) temp,mu,pin(3),nsample,nbcoltot(1)/(nsample*tevol),
     & nbcoltot(2)/(nsample*tevol)
      write(9,*) 'Pperp | dN_c | dN_b'
      do i=1,100
         write(9,10) (i-0.5)/10,distribpt(i,1)*10,distribpt(i,2)*10
      enddo
      close(9)
      open(9,file='dist_Pprime_long_av_mixed_evolveB.res',status='new')
      write(9,*) temp,mu,pin(3),nsample,nbcoltot(1)/(nsample*tevol),
     & nbcoltot(2)/(nsample*tevol)
      write(9,*) 'Pz | dN_c | dN_b'
      do i=-100,100
         write(9,10) i*(pin(3)+10.)/100,distribpz(i,1)*100/pin(3),
     &      distribpz(i,2)*100/pin(3)
      enddo
      close(9)



      write(6,*) 'generating transport coefficients'
      itypq=1
      m=1.5d0
      ifcol=.true.
      ifrad=.false.
      write(6,*) 'for elastic processes'      
      nsample=1000000
      open(9,file='av_moment_varia_el_evolveB.res',status='new')
      write(9,*) m,temp,mu,nsample
      write(9,*) 'p | rate | dp_z/dt | dp_z^2/dt | dp_t^2/dt'  
      do i=1,14
         pin(1)=0
         pin(2)=0
         pin(3)=tbpin(i)
         pin(0)=sqrt(pin(3)**2+m**2)
         meandp=0
         momlong=0.
         momtrans=0.
         nbcoltot(1)=0.
         do j=1,nsample
            do k=0,3
               pp(k)=pin(k)
            enddo
            call evolveboltzmann(itypq,pp,temp,effdeg,tevol,ifcol,ifrad
     &         ,nbcol,.false.)
            nbcoltot(1)=nbcoltot(1)+nbcol
            meandp=meandp+pp(3)
            momlong=momlong+(pin(3)-pp(3))**2
            momtrans=momtrans+(pp(1)**2+pp(2)**2)
         enddo
         ratetot=nbcoltot(1)/(tevol*nsample)
         meandp=pin(3)-meandp/nsample
         momlong=momlong/nsample
         momtrans=momtrans/nsample
         write(6,10) pin(3),ratetot,meandp/tevol,momlong/tevol,
     &     momtrans/tevol
         write(9,10) pin(3),ratetot,meandp/tevol,momlong/tevol,
     &     momtrans/tevol
      enddo
      close(9)

      ifcol=.true.
      ifrad=.true.
      write(6,*) 'for elastic AND radiative processes'      
      nsample=1000000
      open(9,file='av_moment_varia_mixed_evolveB.res',status='new')
      write(9,*) m,temp,mu,nsample
      write(9,*) 'p | rate | dp_z/dt | dp_z^2/dt | dp_t^2/dt'  
      do i=1,14
         pin(1)=0
         pin(2)=0
         pin(3)=tbpin(i)
         pin(0)=sqrt(pin(3)**2+m**2)
         meandp=0
         momlong=0.
         momtrans=0.
         nbcoltot(1)=0.
         do j=1,nsample
            do k=0,3
               pp(k)=pin(k)
            enddo
            call evolveboltzmann(itypq,pp,temp,effdeg,tevol,ifcol,ifrad
     &                   ,nbcol,.false.)
            nbcoltot(1)=nbcoltot(1)+nbcol
            meandp=meandp+pp(3)
            momlong=momlong+(pin(3)-pp(3))**2
            momtrans=momtrans+(pp(1)**2+pp(2)**2)
         enddo
         ratetot=nbcoltot(1)/(tevol*nsample)
         meandp=pin(3)-meandp/nsample
         momlong=momlong/nsample
         momtrans=momtrans/nsample
         write(6,10) pin(3),ratetot,meandp/tevol,momlong/tevol,
     &     momtrans/tevol
         write(9,10) pin(3),ratetot,meandp/tevol,momlong/tevol,
     &     momtrans/tevol
      enddo
      close(9)

      ifcol=.false.
      ifrad=.true.
      write(6,*) 'for pure radiative processes'      
      nsample=1000000
      open(9,file='av_moment_varia_rad_evolveB.res',status='new')
      write(9,*) m,temp,mu,nsample
      write(9,*) 'p | rate | dp_z/dt | dp_z^2/dt | dp_t^2/dt'  
      do i=1,14
         pin(1)=0
         pin(2)=0
         pin(3)=tbpin(i)
         pin(0)=sqrt(pin(3)**2+m**2)
         meandp=0
         momlong=0.
         momtrans=0.
         nbcoltot(1)=0.
         do j=1,nsample
            do k=0,3
               pp(k)=pin(k)
            enddo
            call evolveboltzmann(itypq,pp,temp,effdeg,tevol,ifcol,ifrad
     &            ,nbcol,.false.)
            nbcoltot(1)=nbcoltot(1)+nbcol
            meandp=meandp+pp(3)
            momlong=momlong+(pin(3)-pp(3))**2
            momtrans=momtrans+(pp(1)**2+pp(2)**2)
         enddo
         ratetot=nbcoltot(1)/(tevol*nsample)
         meandp=pin(3)-meandp/nsample
         momlong=momlong/nsample
         momtrans=momtrans/nsample
         write(6,10) pin(3),ratetot,meandp/tevol,momlong/tevol,
     &     momtrans/tevol
         write(9,10) pin(3),ratetot,meandp/tevol,momlong/tevol,
     &     momtrans/tevol
      enddo
      close(9)

      return
      end


      subroutine gimmerates(ifa,p0,temp0,rate_ferm,rate_glu)
      use PBGforratedat1
      use PBGforratedat2
      implicit none
      double precision p0,temp0,p,temp,rate_ferm,rate_glu
      integer ifa,ipinf,ipsup,itinf,itsup
      double precision prap,dp,trap,dtemp,aii,asi,ass,ais
      rate_ferm=0.
      rate_glu=0.
      if(ratemult.le.0) return
      p=p0
      temp=temp0
      if(p.gt.pmax(ifa)) p=0.99*pmax(ifa)
      if(temp.lt.tempmin(ifa)) return
      if(temp.gt.tempmax(ifa)) temp=tempmax(ifa)
      if(temp.eq.tempmax(ifa)) then
         if(p.eq.pmax(ifa))then
            rate_ferm=rf(ifa,nbp(ifa),nbtemp(ifa))
            rate_glu=rg(ifa,nbp(ifa),nbtemp(ifa))
         else
            prap=p/deltp(ifa)
            ipinf=int(prap)
            ipsup=ipinf+1
            dp=prap-ipinf
            rate_ferm=dp*rf(ifa,ipsup,nbtemp(ifa))+
     &                (1-dp)*rf(ifa,ipinf,nbtemp(ifa))
            rate_glu=dp*rg(ifa,ipsup,nbtemp(ifa))+
     &                (1-dp)*rg(ifa,ipinf,nbtemp(ifa))
         endif
      else 
         trap=(temp-tempmin(ifa))/delttemp(ifa)
         itinf=int(trap)
         itsup=itinf+1
         dtemp=trap-itinf
         if(p.eq.pmax(ifa)) then
            rate_ferm=dtemp*rf(ifa,nbp(ifa),itsup)+
     &                (1-dtemp)*rf(ifa,nbp(ifa),itinf)
            rate_glu=dtemp*rg(ifa,nbp(ifa),itsup)+
     &                (1-dtemp)*rg(ifa,nbp(ifa),itinf)
         else
            prap=p/deltp(ifa)
            ipinf=int(prap)
            ipsup=ipinf+1
            dp=prap-ipinf
            aii=(1-dp)*(1-dtemp)
            ais=(1-dp)*dtemp
            asi=dp*(1-dtemp)
            ass=dp*dtemp
            rate_ferm=aii*rf(ifa,ipinf,itinf)+ais*rf(ifa,ipinf,itsup)+
     &            asi*rf(ifa,ipsup,itinf)+ass*rf(ifa,ipsup,itsup)
            rate_glu=aii*rg(ifa,ipinf,itinf)+ais*rg(ifa,ipinf,itsup)+
     &            asi*rg(ifa,ipsup,itinf)+ass*rg(ifa,ipsup,itsup)
         endif
      endif
      rate_ferm=rate_ferm*ratemult
      rate_glu=rate_glu*ratemult
      return
      end

      subroutine testrate
      implicit none
      integer i
      double precision p0,temp0,rate_ferm,rate_glu,ran2
      do i=1,10
         temp0=0.15+0.35*ran2()
         p0=ran2()*60
         call gimmerates(1,p0,temp0,rate_ferm,rate_glu)
         write(6,*) temp0,p0,rate_ferm,rate_glu
      enddo
      return
      end

      subroutine gimmeratiorates(ifa,p0,temp0,ratio_ferm,ratio_glu)
      use PBGforrateratiologdat1
      use PBGforrateratiologdat2
      implicit none
      double precision p0,temp0,logp,temp,ratio_ferm,ratio_glu
      integer ifa,ipinf,ipsup,itinf,itsup
      double precision prap,dp,trap,dtemp,aii,asi,ass,ais
      ratio_ferm=0.
      ratio_glu=0.
      temp=temp0
      logp=log(p0)
      if(logp.gt.logpmaxrad(ifa)) then
         logp=logpmaxrad(ifa) 
      endif
      if(temp.lt.tempminrad(ifa)) temp=tempminrad(ifa)
      if(temp.gt.tempmaxrad(ifa)) temp=tempmaxrad(ifa)
      if(temp.eq.tempmaxrad(ifa)) then
         if(logp.eq.logpmaxrad(ifa))then
            ratio_ferm=rtildef(ifa,nbprad(ifa),nbtemprad(ifa))
            ratio_glu=rtildeg(ifa,nbprad(ifa),nbtemprad(ifa))
         else
            prap=(logp-logpminrad(ifa))/deltprad(ifa)
            ipinf=max(int(prap),0)
            ipsup=ipinf+1
            dp=prap-ipinf
            ratio_ferm=dp*rtildef(ifa,ipsup,nbtemprad(ifa))+
     &                (1-dp)*rtildef(ifa,ipinf,nbtemprad(ifa))
            ratio_glu=dp*rtildeg(ifa,ipsup,nbtemprad(ifa))+
     &                (1-dp)*rtildeg(ifa,ipinf,nbtemprad(ifa))
         endif
      else 
         trap=(temp-tempminrad(ifa))/delttemprad(ifa)
         itinf=int(trap)
         itsup=itinf+1
         dtemp=trap-itinf
         if(logp.ge.logpmaxrad(ifa)) then
            ratio_ferm=dtemp*rtildef(ifa,nbprad(ifa),itsup)+
     &                (1-dtemp)*rtildef(ifa,nbprad(ifa),itinf)
            ratio_glu=dtemp*rtildeg(ifa,nbprad(ifa),itsup)+
     &                (1-dtemp)*rtildeg(ifa,nbprad(ifa),itinf)
         else
            prap=(logp-logpminrad(ifa))/deltprad(ifa)
            ipinf=max(int(prap),0)
            ipsup=ipinf+1
            dp=prap-ipinf
            aii=(1-dp)*(1-dtemp)
            ais=(1-dp)*dtemp
            asi=dp*(1-dtemp)
            ass=dp*dtemp
            ratio_ferm=aii*rtildef(ifa,ipinf,itinf)+
     &                 ais*rtildef(ifa,ipinf,itsup)+
     &                 asi*rtildef(ifa,ipsup,itinf)+
     &                 ass*rtildef(ifa,ipsup,itsup)
            ratio_glu=aii*rtildeg(ifa,ipinf,itinf)+
     &                ais*rtildeg(ifa,ipinf,itsup)+
     &                asi*rtildeg(ifa,ipsup,itinf)+
     &                ass*rtildeg(ifa,ipsup,itsup)
         endif
      endif
      ratio_ferm=exp(ratio_ferm)
      ratio_glu=exp(ratio_glu)
      return
      end


      subroutine testratiorate
      implicit none
      integer i,ifa
      parameter(ifa=1)
      double precision p0,temp0,rate_ferm,rate_glu,ran2
      do i=1,10
         temp0=0.15+0.35*ran2()
         p0=ran2()*60
         call gimmeratiorates(ifa,p0,temp0,rate_ferm,rate_glu)
         write(6,*) temp0,p0,rate_ferm,rate_glu
      enddo
      return
      end

      subroutine initratiorateoneqperplog
      use PBGpsiinfo2
      use PBGforrateratiologdat1
      use PBGforrateratiooneqlogdat2
      implicit none
#include "aaa.h"
      integer itemp,ip
      double precision temp
      character*4 addnamec(2), addnameb
      data addnamec/'1500','1940'/
      data addnameb/'5100'/
      integer ifmtx
      call getMonitorFileIndex(ifmtx)
      open(161,file=fnus1(1:index(fnus1,' ')-1)
     &  //'ratioradel_oneqperp_loglog_1500.dat'
     &  ,status='old')
      read(161,*) tempminrad(1)
      read(161,*) tempmaxrad(1)
      read(161,*) nbtemprad(1)
      if(nbtemprad(1).gt.nbtempmax) then
         write(ifmtx,*) 'error in initrateratio: nbtemp > nbtempmax'
         close(ifmtx)
         close(161)
         stop
      endif
      delttemprad(1)=(tempmaxrad(1)-tempminrad(1))/nbtemprad(1)
      read(161,*) logpminrad(1)
      read(161,*) logpmaxrad(1)
      read(161,*) nbprad(1)
      if(nbprad(1).gt.nbpmax) then
        write(ifmtx,*) 'error in initrate: nbprad > nbpmax'
        close(ifmtx)
        close(161)
        stop
      endif
      deltprad(1)=(logpmaxrad(1)-logpminrad(1))/nbprad(1)
      do itemp=0,nbtemprad(1)
        read(161,*) temp
        do ip=0,nbprad(1)
          read(161,*) rtQCD(1,ip,itemp),rtSQCD(1,ip,itemp)
        enddo
      enddo
      close(161)

      open(162,file=fnus1(1:index(fnus1,' ')-1)
     &  //'ratioradel_oneqperp_loglog_5100.dat'
     &  ,status='old')
      read(162,*) tempminrad(2)
      read(162,*) tempmaxrad(2)
      read(162,*) nbtemprad(2)
      if(nbtemprad(2).gt.nbtempmax) then
         write(ifmtx,*) 'error in initrateratio: nbtemp > nbtempmax'
         close(ifmtx)
         close(162)
         stop
      endif
      delttemprad(2)=(tempmaxrad(2)-tempminrad(2))/nbtemprad(2)
      read(162,*) logpminrad(2)
      read(162,*) logpmaxrad(2)
      read(162,*) nbprad(2)
      if(nbprad(2).gt.nbpmax) then
        write(ifmtx,*) 'error in initrate: nbprad > nbpmax'
        close(ifmtx)
        close(162)
        stop
      endif
      deltprad(2)=(logpmaxrad(2)-logpminrad(2))/nbprad(2)
      do itemp=0,nbtemprad(2)
        read(162,*) temp
        do ip=0,nbprad(2)
          read(162,*) rtQCD(2,ip,itemp),rtSQCD(2,ip,itemp)
        enddo
      enddo
      close(162)


      write(ifmtx,*)'initratiorateoneqperplog performed with success' 
      return
      end


      subroutine gimmeratioratesoneqperp(ifa,p0,temp0,ratio)
      use PBGforradiatGB
      use PBGforrateratiologdat1
      use PBGforrateratiooneqlogdat2
      implicit none
      logical iflargep
      double precision p0,temp0,logp,temp,ratio
      integer ifa,ipinf,ipsup,itinf,itsup
      double precision prap,dp,trap,dtemp,aii,asi,ass,ais
      double precision xmin,complement,pi
      parameter(xmin=0.05d0,pi=3.14159265358979323844d0)
      ratio=0.
      logp=log(p0)
      temp=temp0
      iflargep=.false.
      if(logp.gt.logpmaxrad(ifa)) then
         iflargep=.true.
         if(iffullQCD) then
            complement=-3/pi*(logp-logpmaxrad(ifa))*(3-4*xmin+xmin**2+
     &         4*log(xmin))
         else
            complement=-12/pi*(logp-logpmaxrad(ifa))*(1-xmin+log(xmin))
         endif
         logp=logpmaxrad(ifa)
      endif
      if(temp.lt.tempminrad(ifa)) temp=tempminrad(ifa)
      if(temp.gt.tempmaxrad(ifa)) temp=tempmaxrad(ifa)
      if(temp.eq.tempmaxrad(ifa)) then
         if(logp.eq.logpmaxrad(ifa))then
            if(iffullQCD) then
               ratio=rtQCD(ifa,nbprad(ifa),nbtemprad(ifa))
            else
               ratio=rtQCD(ifa,nbprad(ifa),nbtemprad(ifa))
            endif
         else
            prap=(logp-logpminrad(ifa))/deltprad(ifa)
            ipinf=max(int(prap),0)
            ipsup=ipinf+1
            dp=prap-ipinf
            if(iffullQCD) then
               ratio=dp*rtQCD(ifa,ipsup,nbtemprad(ifa))+
     &            (1-dp)*rtQCD(ifa,ipinf,nbtemprad(ifa))
            else
               ratio=dp*rtSQCD(ifa,ipsup,nbtemprad(ifa))+
     &            (1-dp)*rtSQCD(ifa,ipinf,nbtemprad(ifa))
            endif
         endif
      else 
         trap=(temp-tempminrad(ifa))/delttemprad(ifa)
         itinf=int(trap)
         itsup=itinf+1
         dtemp=trap-itinf
         if(logp.ge.logpmaxrad(ifa)) then
            if(iffullQCD) then
               ratio=dtemp*rtQCD(ifa,nbprad(ifa),itsup)+
     &                (1-dtemp)*rtQCD(ifa,nbprad(ifa),itinf)
            else
               ratio=dtemp*rtSQCD(ifa,nbprad(ifa),itsup)+
     &                (1-dtemp)*rtSQCD(ifa,nbprad(ifa),itinf)
            endif
         else
            prap=(logp-logpminrad(ifa))/deltprad(ifa)
            ipinf=max(int(prap),0)
            ipsup=ipinf+1
            dp=prap-ipinf
            aii=(1-dp)*(1-dtemp)
            ais=(1-dp)*dtemp
            asi=dp*(1-dtemp)
            ass=dp*dtemp
            if(iffullQCD) then
               ratio=aii*rtQCD(ifa,ipinf,itinf)
     &            +ais*rtQCD(ifa,ipinf,itsup)+
     &            asi*rtQCD(ifa,ipsup,itinf)+ass*rtQCD(ifa,ipsup,itsup)
            else
               ratio=aii*rtSQCD(ifa,ipinf,itinf)
     &               +ais*rtSQCD(ifa,ipinf,itsup)+
     &         asi*rtSQCD(ifa,ipsup,itinf)+ass*rtSQCD(ifa,ipsup,itsup)
            endif
         endif
      endif
      ratio=exp(ratio)
      if(iflargep) ratio=ratio+complement
      return
      end

      subroutine testratiorateoneqperp
      implicit none
      integer i,ifa
      parameter(ifa=2)
      double precision p0,temp0,ratio
      temp0=0.293
      do i=1,100
         p0=i*1.d0
         call gimmeratioratesoneqperp(ifa,p0,temp0,ratio)
         write(6,*) temp0,p0,ratio
      enddo
      return
      end

      double precision function omegadsigmadomegaQCDbug(x,mQ,mg,qperp2,
     &  a,c0,lna)
      implicit none
      double precision x,mQ,mg,qperp2,a,c0,sqta,lna
      a=((1-x)*mg**2+(x*mQ)**2)/qperp2
      sqta=sqrt(4*a+1)
      c0=3+sqta+(sqta+1)/a
      lna=log(c0/(sqta-1))
      omegadsigmadomegaQCDbug=
     &         ((1-x)*((a+0.5d0)*lna-1.d0)+0.25d0*x**2*lna)/sqta
      return 
      end


      double precision function omegadsigmadomegaQCD(x,mQ,mg,qperp2,a,
     &                                    c0,lna)
      implicit none
      double precision x,mQ,mg,qperp2,a,c0,sqta,lna
      a=((1-x)*mg**2+(x*mQ)**2)/qperp2
      sqta=sqrt(4*a+1)
      c0=3+sqta+(sqta+1)/a
      lna=log(c0/(sqta-1))/sqta
      omegadsigmadomegaQCD=(1-x)*((a+0.5d0)*lna-1.d0)+0.25d0*x**2*lna
      return 
      end


      double precision function omegadsigmadomegaSQCDbug(x,mQ,mg,qperp2,
     &  a,c0,lna)
      implicit none
      double precision x,mQ,mg,qperp2,a,c0,sqta,lna
      a=((1-x)*mg**2+(x*mQ)**2)/qperp2
      sqta=sqrt(4*a+1)
      c0=3+sqta+(sqta+1)/a
      lna=log(c0/(sqta-1))
      omegadsigmadomegaSQCDbug=(1-x)*((a+0.5d0)*lna-1.d0)/sqta
      return 
      end

      double precision function omegadsigmadomegaSQCD(x,mQ,mg,qperp2,a,
     & c0,lna)
      implicit none
      double precision x,mQ,mg,qperp2,a,c0,sqta,lna
      a=((1-x)*mg**2+(x*mQ)**2)/qperp2
      sqta=sqrt(4*a+1)
      c0=3+sqta+(sqta+1)/a
      lna=log(c0/(sqta-1))
      omegadsigmadomegaSQCD=(1-x)*((a+0.5d0)*lna/sqta-1.d0)
      return 
      end

      subroutine get1x_qQ_radiat(mQ,mg,qperp2,x,propSQCD,a,c0,lna)
      use PBGforradiatGB
      implicit none
      integer it
      double precision mQ,mg,qperp2,x,dsigdx0,dsigdx,
     & omegadsigmadomegaQCD,
     & omegadsigmadomegaSQCD,ran2,xmin,lnxmin,ratio,propSQCD,a,c0,lna
      parameter(xmin=0.05d0)
      logical accept
      integer ifmtx
      call getMonitorFileIndex(ifmtx)
      lnxmin=log(xmin)
      if(iffullQCD) then
         dsigdx0=max(omegadsigmadomegaQCD(0.D0,mQ,mg,qperp2,a,c0,lna),
     &           omegadsigmadomegaQCD(1.D0,mQ,mg,qperp2,a,c0,lna))
      else
         dsigdx0=omegadsigmadomegaSQCD(0.D0,mQ,mg,qperp2,a,c0,lna)
      endif
      accept=.false.
      it=0
      do while(.not.(accept))
         it=it+1
         x=exp(ran2()*lnxmin)
         if(iffullQCD) then
            dsigdx=omegadsigmadomegaQCD(x,mQ,mg,qperp2,a,c0,lna)
         else
            dsigdx=omegadsigmadomegaSQCD(x,mQ,mg,qperp2,a,c0,lna)
         endif
         ratio=dsigdx/dsigdx0
         if(ratio.gt.1) then
           write(ifmtx,*)'ratio >1 found in get1x_qQ_radiat for x=',
     &        x,', mg=',mg,', mQ=',mQ,', qperp^2=',qperp2
            write(ifmtx,*) 'ratio=',dsigdx/dsigdx0
            close(ifmtx)
            stop
         endif
         accept=(ran2().lt.ratio)
         if(it.gt.1000) then
           write(ifmtx,*)'1000 iteration in get1x_qQ_radiat wo success'
            write(ifmtx,*) mQ,mg,qperp2
            close(ifmtx)
            stop
         endif
      enddo
      if(iffullQCD) then
         propSQCD=omegadsigmadomegaSQCD(x,mQ,mg,qperp2,a,c0,lna)/dsigdx
         if(propSQCD.gt.1) then
            write(ifmtx,*) 'prop SQCD >1 in get1x_qQ_radiat'
            close(ifmtx)
            stop
         endif
      else
         propSQCD=1.d0
      endif
      return
      end

      subroutine testget1x_qQ_radiat
      implicit none
      integer i,ifa
      parameter(ifa=2)
      double precision mQ,mg,x,qperp2,ran2,xtot,probsqcd,a,c0,lna
      mQ=1.5
      do i=1,10
         mg=2*(0.15+0.35*ran2())
         qperp2=(ran2()*5)**2
         call get1x_qQ_radiat(mQ,mg,qperp2,x,probsqcd,a,c0,lna)
         write(6,*) mg,sqrt(qperp2),x
      enddo
      mg=0.6
      qperp2=1.44
      xtot=0
      do i=1,1000000
         call get1x_qQ_radiat(mQ,mg,qperp2,x,probsqcd,a,c0,lna)
         xtot=xtot+x
      enddo
      write(6,*) '<x> for mg=',mg,' and qperp^2=',qperp2,' is ',
     &          xtot/1000000
      return
      end

      subroutine get1_veck_radiat_Bdist(qperp2,a,c0,lna,kperp,cosphi,
     &           sinphi)
      implicit none
      double precision pi,qperp2,qperp,a,c0,lna,kperp2,kperp,cosphi,
     & sinphi,
     & ran2,c1,c12,tildeb,bid
      parameter (pi=3.14159265358979323844d0)      
      c1=c0*exp(-ran2()*lna)
      c12=c1*c1
      kperp2=(4*a**2-c12*a+6*c1*a+2*c1)/(c12+2*c1-4*a)*qperp2
      kperp=sqrt(kperp2)
      qperp=sqrt(qperp2)
      tildeb=2*kperp*qperp/((1+a)*qperp2+kperp2)
      bid=(1-tildeb)/((1+tildeb)*tan(ran2()*pi/2)**2)
      cosphi=(1-bid)/(1+bid)
      sinphi=sqrt(1-cosphi**2)
      if(ran2().gt.0.5d0) then
         sinphi=-sinphi
      endif
      return
      end

      subroutine get1_veck_radiat_Adist(mgt2,qperp2,a,c0,lna,kperp,
     &  cosphi,
     &  sinphi)
      implicit none
      logical accept
      integer it
      double precision qperp2,qperp,a,c0,lna,kperp2,kperp,cosphi,sinphi,
     & ran2,mgt2,rap,rej
      integer ifmtx
      call getMonitorFileIndex(ifmtx)
      accept=.false.
      it=0
      do while(.not.(accept))
         it=it+1
         call get1_veck_radiat_Bdist(qperp2,a,c0,lna,kperp,cosphi,
     & sinphi)
         qperp=sqrt(qperp2)
         kperp2=kperp**2
         rap=1+(qperp2-2*qperp*kperp*cosphi)/(kperp2+mgt2)
         rej=1-a*(rap+1/rap-2)
         if(rej.gt.1) then
            write(ifmtx,*) 'rejection>1 in get1_veck_radiat_Adist'
            write(ifmtx,*) qperp,kperp,cosphi
            close(ifmtx)
            stop
         endif
         accept=(ran2().lt.rej)
         if(it.gt.1000) then
            write(ifmtx,*) '>1000 iterations in get1_veck_radiat_Adist'
            write(ifmtx,*) mgt2,qperp2,a
            close(ifmtx)
            stop
         endif
      enddo
      return
      end

      subroutine test_get1_veck_radiat_dist
      implicit none
      integer i
      double precision mg,mQ,x,mgt2,qperp2,kperp,cosphi,sinphi,ktx,kty,
     & a,sqta,c0,lna,meankalongqperp,meanktransqperp,meankalongqperp2,
     & meanktransqperp2,meank 
      mg=0.6
      mQ=1.5
      x=0.1
      qperp2=4
      mgt2=(1-x)*mg**2+(x*mQ)**2
      a=mgt2/qperp2
      sqta=sqrt(4*a+1)
      c0=3+sqta+(sqta+1)/a
      lna=log(c0/(sqta-1))
      meankalongqperp=0
      meanktransqperp=0
      meankalongqperp2=0
      meanktransqperp2=0
      meank=0
      do i=1,2000000
         call get1_veck_radiat_Adist(mgt2,qperp2,a,c0,lna,kperp,cosphi,
     &   sinphi)
         ktx=kperp*cosphi
         kty=kperp*sinphi
         meankalongqperp=meankalongqperp+ktx
         meankalongqperp2=meankalongqperp2+abs(ktx)**1.5
         meanktransqperp=meanktransqperp+kty
         meanktransqperp2=meanktransqperp2+abs(kty)**1.5
         meank=meank+kperp
      enddo
      meankalongqperp=meankalongqperp/2000000
      meanktransqperp=meanktransqperp/2000000
      meankalongqperp2=meankalongqperp2/2000000
      meanktransqperp2=meanktransqperp2/2000000
      meank=meank/2000000
      write(6,*) '<kt>:',meank
      write(6,*) '<kt_x>:',meankalongqperp,' & <kt_y>:',meanktransqperp
      write(6,*) '<kt^1.5_x>:',meankalongqperp2,
     & ' & <kt^1.5_y>:',meanktransqperp2
      return
      end


      subroutine get1x_veck_qQ_radiat(itypq,p,temp,mQ,qperp2,x,kperp,
     &  cosphi,sinphi,ifail)
      use PBGforradiatGB
      use PBGforradiatLPM
      use PBGgluprops
      implicit none
      integer ifail,itypq
      double precision mQ,mg,qperp2,x,kperp,cosphi,sinphi,propsqcd,a,
     & c0,lna,
     &  mgt2,ran2,supprLPM,p,temp,supprDamp1,supprDamp3,suppr
      logical corrx2
      integer ifmtx
      call getMonitorFileIndex(ifmtx)
      ifail=0
      mg=c1*temp
      call get1x_qQ_radiat(mQ,mg,qperp2,x,propSQCD,a,c0,lna)
      suppr=1
      if(typelpmdamp.eq.1) then
         if(iflpm) then
            if(ifdamp) then
               suppr=supprDamp1(itypq,p,mQ,x,temp)
            else
               suppr=supprLPM(itypq,p,mQ,x,temp)
            endif
         endif
      else if(typelpmdamp.eq.3) then
         suppr=supprDamp3(itypq,p,mQ,x,temp,iflpm,ifdamp)
      else
         write(ifmtx,*) 'unrecog. typelpmdamp in get1x_veck_qQ_radiat'
         close(ifmtx)
         stop
      endif
      if(ran2().gt.suppr) then
         ifail=1
         return
      endif
      corrx2=.false.
      if(iffullQCD) then
         if(ran2().gt.propSQCD) corrx2=.true.
      endif

      if(corrx2) then
         call get1_veck_radiat_Bdist(qperp2,a,c0,lna,kperp,cosphi,
     &                sinphi)
      else
         mgt2=(1-x)*mg**2+(x*mQ)**2
         call get1_veck_radiat_Adist(mgt2,qperp2,a,c0,lna,kperp,
     &           cosphi,sinphi)
      endif
      return
      end


      subroutine test_get1x_veck_radiat_dist
      implicit none
      integer i,nb,itypq,ifail,nbperp,iperp,ifail2,nbtry,nbsuccess,ix
      double precision mQ,x,qperp2,kperp,cosphi,sinphi,ktx,
     & meankalongqperp,meanx,meankalongqperp2,meanx2,
     & temp,p,
     & distribpperp(100),ppt2,ppt,pperpmax,steppperp,s,pptx,ppty,qperp,
     & pplnew,pplus,mpt2,distribx(100),mg

      temp=0.3
      itypq=2
      mQ=5.1D0
      qperp2=0.1d0
      s=50.D0
      p=(s-mQ**2)/(2*temp)
      nbperp=100
      pperpmax=5.D0
      steppperp=pperpmax/nbperp
      do iperp=1,nbperp
         distribpperp(iperp)=0.D0
      enddo
      do ix=1,nbperp
         distribx(ix)=0.D0
      enddo
      nb=0
      nbtry=40000000
      nbsuccess=0      
      meankalongqperp=0.D0
      meankalongqperp2=0.D0
      meanx=0.D0
      meanx2=0.D0
      do i=1,nbtry
         call get1x_veck_qQ_radiat(itypq,p,temp,mQ,qperp2,x,kperp,
     &            cosphi,
     &   sinphi,ifail)
         if(ifail.eq.0) then
            pplus=sqrt(p**2+mQ**2)+p
            pptx=sqrt(qperp2)-kperp*cosphi
            ppty=kperp*sinphi
            ppt2=pptx**2+ppty**2
            mpt2=mQ**2+ppt2
            qperp=sqrt(qperp2)
            mg=0d0
            call gennewPradiat(x,s,pplus,mpt2,mg,kperp,
     &                     qperp,pplnew,ifail2)
            if(ifail2.le.1) then
               nbsuccess=nbsuccess+1
               ktx=kperp*cosphi
               if((x.gt.0).and.(x.lt.1).and.(abs(ktx).lt.5)) then 
                  nb=nb+1
                  meankalongqperp=meankalongqperp+ktx
                  meankalongqperp2=meankalongqperp2+abs(ktx)**1.5
                  meanx=meanx+x
                  meanx2=meanx2+x**2
                  ix=int(x/0.01)+1
                  distribx(ix)=distribx(ix)+1.D0/0.01
               endif
               ppt=sqrt(ppt2)
               if(ppt.lt.pperpmax) then
                  iperp=int(ppt/steppperp)+1
                  distribpperp(iperp)=distribpperp(iperp)+1.D0/steppperp
               endif
            endif
         endif    
      enddo
      meankalongqperp=meankalongqperp/nb
      meankalongqperp2=meankalongqperp2/nb
      meanx=meanx/nb
      meanx2=meanx2/nb
      write(6,*) nb
      write(6,*) '<kt_x>:',meankalongqperp,' & <kt_x^1.5>:',
     &                meankalongqperp2
      write(6,*) '<x>:',meanx,' & <x^2>:',meanx2
      open(9,file='distrib_pperp_SQCD_fixed_qt.res',status='new')
      write(9,*) mQ,temp,s,qperp2,nbtry,nbsuccess
      write(9,*) 'pt | dPg/dpt'  
      do iperp=1,nbperp
         write(9,*) (iperp-0.5)*steppperp,distribpperp(iperp)/nbtry
      enddo
      close(9)
      open(9,file='distrib_x_SQCD_fixed_qt.res',status='new')
      write(9,*) mQ,temp,s,qperp2,nbtry,nbsuccess
      write(9,*) 'x | dPg/dx'  
      do ix=1,nbperp
         write(9,*) (ix-0.5)*0.01,distribx(ix)/nbtry
      enddo
      close(9)
      return
      end


      subroutine test_qperp_and_exit_channel
      implicit none
      logical ifrunning,ifferm,iffaith
      integer itypq,nbsamp,nbsuccess,isamp,nit,nbdist,ipt,iz,i,ifail,
     & ifail2
      double precision s,mQ,uvec(3),q,qvec(3),u2,temp,mu,t,prelcm2,
     & steppt,
     & ptmax,distribpt(200),qperp2,x,kperp,cosphi,sinphi,pplus,pptx,
     & ppty,
     & ppt2,mpt2,qperp,pplnew,ppt,pvec,p(0:3),al,
     & mdebeffinterp,muoftemp,distribz(200),distribx(200),stepz,tav,
     & tavsucc,
     & xav,xav2,ppt2av,zav,mg
      parameter (nit=20)
      s=10
      temp=0.3
      nbsamp=10000000
      nbsuccess=0
      ifferm=.true.
      ifrunning=.false.
      iffaith=.true.
      if(ifrunning) then
         mu=0.4472136d0*mdebeffinterp(temp)      
      else
         mu=sqrt(0.15)*muoftemp(temp)
      endif
      itypq=1
      mQ=1.5D0
      p(1)=0.D0
      p(2)=0.D0
      p(3)=20.D0
      pvec=p(1)**2+p(2)**2+p(3)**2
      p(0)=sqrt(mQ**2+pvec)
      pvec=sqrt(pvec)
      do i=1,3
         uvec(i)=-p(i)/mQ
      enddo
      u2=uvec(1)**2+uvec(2)**2+uvec(3)**2      
      q=(s-mQ**2)/(2*mQ)
      qvec(1)=0.D0
      qvec(2)=0.D0
      qvec(3)=-q
      prelcm2=(s-mQ**2)**2/(4*s)
      ptmax=10.D0
      nbdist=100
      steppt=ptmax/nbdist
      stepz=1.D0/nbdist
      do ipt=1,nbdist
         distribpt(ipt)=0.D0
         distribz(ipt)=0.D0
         distribx(ipt)=0.D0
      enddo
      tav=0.D0
      xav=0.D0
      xav2=0.D0
      tavsucc=0.D0
      ppt2av=0.D0
      zav=0.D0
      do isamp=1,nbsamp
         call getoneqperprandsimp(itypq,ifferm,temp,mQ,s,mu,t)
         qperp2=-t
         tav=tav+qperp2
         call get1x_veck_qQ_radiat(itypq,pvec,temp,mQ,qperp2,x,kperp,
     & cosphi,
     &   sinphi,ifail)
         if(ifail.eq.0) then
            pplus=p(0)+pvec
            qperp=sqrt(qperp2)
            pptx=qperp-kperp*cosphi
            ppty=-kperp*sinphi
            ppt2=pptx**2+ppty**2
            mpt2=mQ**2+ppt2
            mg=0d0
            call gennewPradiat(x,s,pplus,mpt2,mg,kperp,
     &                     qperp,pplnew,ifail2)
            if(ifail2.le.1) then
               nbsuccess=nbsuccess+1
               tavsucc=tavsucc+qperp2
               xav=xav+x
               xav2=xav2+x**2
               ppt2av=ppt2av+ppt2
               if((x.lt.1.D0).and.(x.gt.0.D0)) then
                  iz=int(x/stepz)
                  distribx(iz)=distribx(iz)+1.D0/stepz
               endif
               al=pplnew/pvec
               zav=zav+al
               if((al.lt.2.D0).and.(al.gt.0.D0)) then
                  iz=int(al/stepz)
                  distribz(iz)=distribz(iz)+1.D0/stepz
               endif
               ppt=sqrt(ppt2)
               if(ppt.lt.ptmax) then
                  ipt=int(ppt/steppt)+1
                  distribpt(ipt)=distribpt(ipt)+1.D0/steppt
               endif
            endif
         endif
      enddo
      open(9,file='distrib_pperp_SQCD_fixed_s.res',status='new')
      write(9,*) mQ,s,pvec,temp,mu,nbsamp,nbsuccess
      write(9,*) ifferm,ifrunning,iffaith
      write(9,*) 'pt | dN/dpt'  
      do ipt=1,nbdist
         write(9,*) (ipt-0.5)*steppt,distribpt(ipt)/nbsamp
      enddo
      close(9)
      open(9,file='distrib_plong_SQCD_fixed_s.res',status='new')
      write(9,*) mQ,s,pvec,temp,mu,nbsamp,nbsuccess
      write(9,*) ifferm,ifrunning,iffaith
      write(9,*) 'z=pz/pzinit | dN/dz'  
      do iz=1,nbdist
         write(9,*) (iz-0.5)*stepz,distribz(iz)/nbsamp
      enddo
      close(9)
      open(9,file='distrib_x_SQCD_fixed_s.res',status='new')
      write(9,*) mQ,s,pvec,temp,mu,nbsamp,nbsuccess
      write(9,*) ifferm,ifrunning,iffaith
      write(9,*) 'x | dN/dx'  
      do iz=1,nbdist
         write(9,*) (iz-0.5)*stepz,distribx(iz)/nbsamp
      enddo
      close(9)
      write(6,*) '<t> no cut:',tav/nbsamp
      write(6,*) '<t> :',tavsucc/nbsuccess
      write(6,*) '<pt^2> :',ppt2av/nbsuccess
      write(6,*) '<ppl/pl>:',zav/nbsuccess
      write(6,*) '<x> :',xav/nbsuccess
      write(6,*) '<x^2> :',xav2/nbsuccess
      return
      end


      subroutine testgenoneqperpandinvariantexit
      implicit none
      integer itypq,ifailglob,nbsamp,nbsuccess,isamp
      logical ifferm,ifrunning
      double precision mQ,pl,q,temp,mu,s,qperp2,x,kperp,
     & pptx,
     & ppty,pplnew,pplnewshift,xav,mdebeffinterp,muoftemp,plav
      temp=0.3
      nbsamp=10000000
      nbsuccess=0
      ifferm=.true.
      ifrunning=.false.
      if(ifrunning) then
         mu=0.4472136d0*mdebeffinterp(temp)      
      else
         mu=sqrt(0.15)*muoftemp(temp)
      endif
      itypq=1
      mQ=1.5D0
      pl=20.D0
      s=200
      q=(s-mQ**2)/(2*mQ)
      xav=0.D0
      plav=0.D0
      do isamp=1,nbsamp
         call genoneqperpandinvariantexit(itypq,ifferm,mQ,pl,q,temp,mu,
     &  ifailglob,s,qperp2,x,kperp,pptx,ppty,pplnew)
         if(ifailglob.eq.0) then
            nbsuccess=nbsuccess+1
            xav=xav+x
            plav=plav+pplnewshift
         endif
      enddo
      write(6,*) 'prob succ:', real(nbsuccess)/nbsamp
      write(6,*) 'xav:',xav/nbsuccess
      write(6,*) '<pl>',plav/nbsuccess/pl
      end



      double precision function qhattilde(itypq,p0,temp0)
      use PBGpsiinfo
      use PBGupsinfo
      implicit none
      integer itypq
      double precision p0,temp0,adrag,blong,btrans,m
      call gimmefpcoeff(itypq,p0,temp0,adrag,blong,btrans)
      if(itypq.eq.1) then
         m=mc
      else if(itypq.eq.2) then
         m=mb
      endif
      if(p0.eq.0) then
         write(8,*) 'p0=0 in qhattilde' 
         qhattilde=4.5d0*btrans*m/0.001
      else
         qhattilde=4.5d0*btrans*sqrt(1+(m/p0)**2)
      endif
      return
      end

      double precision function lfmult(x,mg,mQ,en,qhat,phidec)
      implicit none
      double precision x,mg,mQ,en,qhat,phidec,bid,w
      bid=0.5d0*(mg**2+(x*mQ)**2)
      w=x*en
      lfmult=2*w*phidec/(bid+sqrt(qhat*w*phidec+bid**2))
      return
      end

      double precision function lfmultunit(x,mg,mQ,en,qhat,phidec)
      implicit none
      double precision x,mg,mQ,en,qhat,phidec,lfmult
      lfmultunit=lfmult(x,mg,mQ,en,qhat*0.1973d0,phidec)*0.1973d0
      return
      end

      double precision function lfmultmodelE(itypq,p,x,temp)
      use PBGpsiinfo
      use PBGupsinfo
      use PBGgluprops
      implicit none
      integer itypq
      double precision:: p,temp,x,mQ,en,qhat,qhattilde,lfmultunit,mg
      integer ifmtx
      call getMonitorFileIndex(ifmtx)
      if(itypq.eq.1) then
         mQ=mc
      else if(itypq.eq.2) then
         mQ=mb
      else
         write(ifmtx,*) 'no type of mass for ityp=',itypq
         close(ifmtx)
         stop
      endif
      en=sqrt(mQ**2+p**2)
      mg=c1*temp
      qhat=qhattilde(itypq,p,temp)
      lfmultmodelE=lfmultunit(x,mg,mQ,en,0.5*qhat,2.d0)
      return
      end function lfmultmodelE


      double precision function tildemeanfree(itypq,p,m,temp)
      implicit none
      integer itypq
      double precision p,m,temp,ratef,rateg,ratetot 
      call gimmerates(itypq,p,temp,ratef,rateg)
      ratetot=ratef+rateg
      tildemeanfree=8*p/(9*sqrt(p**2+m**2)*ratetot)
      return
      end

      double precision function nbcoh(itypq,p,m,x,temp)
      implicit none
      integer itypq
      double precision p,m,temp,x,tildemeanfree,lfmultmodelE
      nbcoh=max(lfmultmodelE(itypq,p,x,temp)
     &        /tildemeanfree(itypq,p,m,temp),1.D0)
      return
      end

      double precision function d2IdxdzGB(itypq,p,m,x,temp,alphas)
      implicit none
      integer itypq,nc
      parameter (nc=3)
      double precision pi,p,m,x,temp,alphas,tildemeanfree,qhattilde,mfp
      parameter (pi=3.14159265358979323844d0)
      mfp=tildemeanfree(itypq,p,m,temp)
      d2IdxdzGB=2*nc*alphas/(pi*mfp)*log(1+mfp*qhattilde(itypq,p,temp)/
     &   (3*(4*temp**2+(x*m)**2)))*(1-x)
      return
      end

      double precision function d2IdxdzLPM(itypq,p,m,x,temp,alphas)
      implicit none
      integer itypq,nc,ncohcr
      parameter (nc=3,ncohcr=10)
      double precision pi,p,m,x,temp,alphas,tildemeanfree,lfmult,
     & lfmultmodelE,
     & qhattilde,mfp,nbcoh,lf,dec,w
      parameter (pi=3.14159265358979323844d0)
      mfp=tildemeanfree(itypq,p,m,temp)
      lfmult=lfmultmodelE(itypq,p,x,temp)
      lf=max(lfmult,mfp)
      nbcoh=lf/mfp
      dec=1.d0-exp(-(nbcoh-1)/ncohcr)
      w=x*sqrt(m**2+p**2)
      d2IdxdzLPM=(2*nc+(40.d0/3-2*nc)*dec)*alphas/(pi*lf)*log(1+
     &   lf*qhattilde(itypq,p,temp)/(3*(4*temp**2+(x*m)**2+
     &   2*0.1973d0*w/lf*dec)))*(1-x)
      return
      end

      double precision function supprLPM(itypq,p,m,x,temp)
      implicit none
      integer itypq
      double precision p,m,x,temp,d2IdxdzLPM,d2IdxdzGB
      supprLPM=min(d2IdxdzLPM(itypq,p,m,x,temp,1.d0)/
     &  d2IdxdzGB(itypq,p,m,x,temp,1.d0),1.d0)
      return
      end

      double precision function supprDamp1(itypq,p,m,x,temp)
      use PBGgluprops
      implicit none
      integer itypq
      double precision:: p,m,x,temp,hbarc,en,mg
      double precision:: lfsing,lfmult,lf
      double precision:: GamD,ldamp
      double precision:: lftilde,ratiolf,supprLPM
      parameter(hbarc=0.197327d0)

      en=sqrt(p**2+m**2)
      mg=c1*temp
      GamD=c2*temp

      lfsing=2.d0*x*(1.d0-x)*en*hbarc/(mg**2+x**2*m**2)
      lfmult=lfsing*supprLPM(itypq,p,m,x,temp)
      lf=min(lfsing,lfmult)
      ldamp=hbarc/GamD
      lftilde=min(lf,ldamp)
      ratiolf=lftilde/lfsing
      supprDamp1=ratiolf
      return
      end function supprDamp1


      double precision function supprDamp3(itypq,p,m,x,temp,iflpm,
     &   ifdamp)
      use PBGgluprops
      implicit none
      integer itypq
      logical iflpm,ifdamp
      double precision:: p,m,x,temp,en,mg,hbarc
      double precision:: qhat,qhattilde
      double precision:: tftilde,tfGB,tfsing,tfmult
      double precision:: GamD,ratio,tdamp

      parameter(hbarc=0.197327d0)
      en=sqrt(p**2+m**2)
      mg=c1*temp
      GamD=c2*temp
      qhat=qhattilde(itypq,p,temp)
      tfGB=2.d0*x*en*hbarc/mg**2
      tfsing=2.d0*x*(1.d0-x)*en*hbarc/(x**2*m**2+mg**2*(1.d0-x))
      if(iflpm) then
         tfmult=sqrt(2.d0*x*en*hbarc/((1.d0-x)*qhat))
         if(ifdamp) then
            tdamp=hbarc/GamD
            tftilde=min(tfsing,tfmult,tdamp)
         else
            tftilde=min(tfsing,tfmult)
         endif
      else
         if(ifdamp) then
            tdamp=hbarc/GamD
            tftilde=min(tfsing,tdamp)
         else
            tftilde=tfsing
         endif 
      endif
      ratio=tftilde/tfsing
      supprDamp3=ratio
      return
      end function supprDamp3



