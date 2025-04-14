      subroutine initptspectrum_of_b(ib,itypeq,ispect,sigbid,ifshadow,
     & ishadow)
      use PBGblckinitptspect
      use PBGblckinitptspect2
      use PBGblckinitptspect3
      use PBGsystem
      implicit none
#include "aaa.h"       
      logical ifshadow
      integer ib,lentrue
      integer nbspect,itypeq,ishadow,ipt,nptbin2,ispect,j
      double precision pt,dndpt(4),shadow,ptstep2,
     & sigbid,sigload(4),shadowload(3),weight
      character*60 spectrname,spectrnameshadow,
     & spectrnamebas(2),spectrnamebasshadow(2)
      character*2 spectrnameshadow2,addname,tochain
      character*4 spectrnameshadow3
      data spectrnamebasshadow/
     &   'FONLLEMMI/shadowc','FONLLEMMI/shadowb'/
      data spectrnamebas/
     &   'FONLLEMMI/dncdpt','FONLLEMMI/dnbdpt'/
      integer ifmtx
      call getMonitorFileIndex(ifmtx)
      write(ifmtx,*) 'entering iniptspectrum'
      write(ifmtx,*) sqrtS,A,ib
      if(A.eq.129)then
         spectrnameshadow2='Xe'
      elseif(A.eq.197) then
         spectrnameshadow2='Au'
      elseif(A.eq.208) then
         spectrnameshadow2='Pb'
      else
         write(ifmtx,*) 'unrecognized A in initptspectrumofb'
         write(ifmtx,*) A
         close(ifmtx)
         stop
      endif
      if(abs(sqrts-200).lt.0.1)then
         spectrnameshadow3='0200'
      elseif(abs(sqrts-2755).lt.10.) then
         spectrnameshadow3='2760'
      elseif(abs(sqrts-5250).lt.260.) then
         spectrnameshadow3='5020'
      else
         write(ifmtx,*) 'unrecognized energy in initptspectrumofb'
         write(ifmtx,*) sqrts
         close(ifmtx)
         stop
      endif
      addname=tochain(ib)
      spectrname=spectrnamebas(itypeq)
      spectrname=spectrname(1:lentrue(spectrname,60))//spectrnameshadow3
     & //'.dat'
      spectrnameshadow=spectrnamebasshadow(itypeq)
      spectrnameshadow=spectrnameshadow(1:lentrue(spectrnameshadow,60))
     & //spectrnameshadow2//'_'//spectrnameshadow3//'_'//addname//'.dat'
      nofpt(itypeq,0)=0.d0
      nofptAA(itypeq,0)=0.d0
      write(6,*) 'opening file ',spectrname
      open(157,file=fnus1(1:index(fnus1,' ')-1)
     &  //spectrname(1:lentrue(spectrname,60)),status='old')
      read(157,*)
      read(157,*)
      read(157,*) nbspect
      read(157,*) (sigload(j),j=1,nbspect)
      sigbid=sigload(ispect)
      read(157,*)
      read(157,*)
      read(157,*) nptbin(itypeq),ptstep(itypeq)
      do ipt=1,nptbin(itypeq)
         read(157,*) pt,(dndpt(j),j=1,nbspect)
         nofpt(itypeq,ipt)=dndpt(ispect)
      enddo
      close(157)
      ptmax(itypeq)=pt+ptstep(itypeq)
      if(ifshadow) then
         write(6,*) 'opening file ',spectrnameshadow
         open(157,file=fnus1(1:index(fnus1,' ')-1)
     &  //spectrnameshadow(1:lentrue(spectrnameshadow,60)),status='old')
         read(157,*)
         read(157,*) nptbin2,ptstep2
         if((nptbin2.ne.nptbin(itypeq)).or.
     &        (ptstep2.ne.ptstep(itypeq))) then
            write(ifmtx,*) 'unmatched grid in initptspectrum'
            close(ifmtx)
            stop
         endif
      endif
      do ipt=1,nptbin(itypeq)
         if(ifshadow) then
            read(157,*) pt,(shadowload(j),j=1,3)
            shadow=shadowload(ishadow)
         else
            shadow=1.d0
         endif
            weight=1.d0
         nofptAA(itypeq,ipt)= nofptAA(itypeq,ipt-1)+nofpt(itypeq,ipt)*
     &       shadow*weight
         nofpt(itypeq,ipt)= nofpt(itypeq,ipt-1)+nofpt(itypeq,ipt)*
     &       weight
      enddo
      if(ifshadow) close(157)
      globshadow(itypeq)=nofptAA(itypeq,nptbin(itypeq))/
     &     nofpt(itypeq,nptbin(itypeq))
      if(globshadow(itypeq).gt.1) then
         write(ifmtx,*) 'global shadow > 1; change strategy'
         write(ifmtx,*) globshadow(itypeq)
         close(ifmtx)
         stop
      endif
      do ipt=1,nptbin(itypeq)-1
         nofpt(itypeq,ipt)= nofpt(itypeq,ipt)/
     &       nofpt(itypeq,nptbin(itypeq))
         nofptAA(itypeq,ipt)= nofptAA(itypeq,ipt)/
     &       nofptAA(itypeq,nptbin(itypeq))
      enddo
      nofpt(itypeq,nptbin(itypeq))=1.0d0
      nofptAA(itypeq,nptbin(itypeq))=1.0d0
      call buildspline(itypeq,ifshadow)
      return
      end


      double precision function interpnofpt(itypeq,pt,ifshadow)
      implicit none
      integer method,itypeq
      parameter (method=3)
      logical ifshadow
      double precision pt,interpnofpt_lin,interpnofpt_quad,
     &  interpnofpt_splin
      integer ifmtx
      call getMonitorFileIndex(ifmtx)
      if(method.eq.1) then
         interpnofpt=interpnofpt_lin(itypeq,pt,ifshadow)
      else if(method.eq.2) then
         interpnofpt=interpnofpt_quad(itypeq,pt,ifshadow)
      else if(method.eq.3) then
         interpnofpt=interpnofpt_splin(itypeq,pt,ifshadow)
      else
         write(ifmtx,*) 'interpnofpt: unknown method'
         close(ifmtx)
         stop
      endif
      end

      double precision function interpnofpt_splin(itypeq,pt,ifshadow)
      use PBGblckinitptspect
      implicit none
      integer itypeq,iinf,isup
      double precision pt,dptinf,dptsup,base,dpt
      logical first(2),ifshadow
      save first
      data first/.true.,.true./
      dpt=ptstep(itypeq)
      iinf=int(pt/ptstep(itypeq))
      isup=iinf+1
      dptinf=pt-iinf*ptstep(itypeq)
      dptsup=isup*ptstep(itypeq)-pt
      base=(tbcoeffM(itypeq,iinf)*(dptsup**3/dpt-dpt*dptsup)+
     &      tbcoeffM(itypeq,isup)*(dptinf**3/dpt-dpt*dptinf))/6.d0
      if(ifshadow) then
         interpnofpt_splin=base+(nofptAA(itypeq,iinf)*dptsup+ 
     &                       nofptAA(itypeq,isup)*dptinf)/dpt
      else
         interpnofpt_splin=base+(nofpt(itypeq,iinf)*dptsup+ 
     &                       nofpt(itypeq,isup)*dptinf)/dpt
      endif
      end

      double precision function interplagrange(n,tbx,tbval,x)
      implicit none
      integer n,i,j
      double precision tbx(5),tbval(5),x,tbweight(5)
      integer ifmtx
      call getMonitorFileIndex(ifmtx)
      if(n.gt.5) then
         write(ifmtx,*) 'error interplagrange'
         close(ifmtx)
         stop
      endif
      interplagrange=0.d0
      do i=1,n
         tbweight(i)=1.d0
         do j=1,i-1
            tbweight(i)=tbweight(i)*(x-tbx(j))/(tbx(i)-tbx(j))
         enddo
         do j=i+1,n
            tbweight(i)=tbweight(i)*(x-tbx(j))/(tbx(i)-tbx(j))
         enddo
         interplagrange=interplagrange+tbweight(i)*tbval(i)
      enddo
      end
      

      subroutine buildspline(itypeq,ifshadow)
      use PBGblckinitptspect
      use PBGblckinitptspect3
      implicit none
      logical ifshadow
      integer i,itypeq,n
      double precision tba(0:nptbinmax),tbb(0:nptbinmax),
     &  tbc(0:nptbinmax),tb2nd(0:nptbinmax),tbalpha(0:nptbinmax),
     &  tbbeta(0:nptbinmax),tby(0:nptbinmax)
      n=nptbin(itypeq)-1
      tba(0)=2.d0
      tbc(0)=1.d0
      if(ifshadow) then
         tb2nd(0)=6.d0*nofptAA(itypeq,1)/ptstep(itypeq)**2
      else
         tb2nd(0)=6.d0*nofpt(itypeq,1)/ptstep(itypeq)**2
      endif         
      do i=1,nptbin(itypeq)-1 
         tba(i)=4.d0
         tbc(i)=1.d0
         tbb(i)=1.d0
         if(ifshadow) then
            tb2nd(i)=6.d0*(nofptAA(itypeq,i+1)-2*nofptAA(itypeq,i)+
     &           nofptAA(itypeq,i-1))/ptstep(itypeq)**2
         else
            tb2nd(i)=6.d0*(nofpt(itypeq,i+1)-2*nofpt(itypeq,i)+
     &           nofpt(itypeq,i-1))/ptstep(itypeq)**2
         endif                     
      enddo
      tbalpha(0)=tba(0)
      tby(0)=tb2nd(0)
      do i=1,n
         tbbeta(i)=tbb(i)/tbalpha(i-1)
         tbalpha(i)=tba(i)-tbbeta(i)*tbc(i-1)
         tby(i)=tb2nd(i)-tbbeta(i)*tby(i-1)
      enddo
      tbcoeffM(itypeq,n+1)=0.d0
      tbcoeffM(itypeq,n)=tby(n)/tbalpha(n)
      do i=n-1,0,-1
         tbcoeffM(itypeq,i)=(tby(i)-tbc(i)*tbcoeffM(itypeq,i+1))/
     &        tbalpha(i)
      enddo
      return
      end

      double precision function interpnofpt_quad(itypeq,pt,ifshadow)
      use PBGblckinitptspect
      use PBGblckinitptspect2
      use PBGblckinitptspect3
      implicit none
      integer itypeq,iinf,isup,i
      logical ifshadow
      double precision pt,tbx(5),tbval(5),interplagrange
      integer ifmtx
      call getMonitorFileIndex(ifmtx)
      if(pt.le.0) then
         interpnofpt_quad=0.D0
      else if(pt.ge.ptmax(itypeq)) then  
         interpnofpt_quad=1.D0
      else
         iinf=int(pt/ptstep(itypeq))
         isup=iinf+2
         if(isup.gt.nptbin(itypeq)) then
            isup=nptbin(itypeq)
            iinf=isup-2
            if(iinf.lt.0) then
               write(ifmtx,*) 'problem 1 interpnofpt'
               close(ifmtx)
               stop
            endif
         endif
         do i=iinf,isup
            tbx(i-iinf+1)=i*ptstep(itypeq)
            if(ifshadow) then
               tbval(i-iinf+1)=nofptAA(itypeq,i)
            else
               tbval(i-iinf+1)=nofpt(itypeq,i)
            endif
         enddo
         interpnofpt_quad=interplagrange(3,tbx,tbval,pt)
      endif
      return
      end

      double precision function interpnofpt_lin(itypeq,pt,ifshadow)
      use PBGblckinitptspect
      use PBGblckinitptspect2
      use PBGblckinitptspect3
      implicit none
      integer itypeq,iinf
      logical ifshadow
      double precision pt,dmpt,dspt
      integer ifmtx
      call getMonitorFileIndex(ifmtx)
      if(pt.le.0) then
         interpnofpt_lin=0.D0
      else if(pt.ge.ptmax(itypeq)) then  
         interpnofpt_lin=1.D0
      else if(pt.le.ptstep(itypeq)) then
         if(ifshadow) then
           interpnofpt_lin=nofptAA(itypeq,1)*(pt/ptstep(itypeq))**2
         else
            interpnofpt_lin=nofpt(itypeq,1)*(pt/ptstep(itypeq))**2
         endif
      else
         dmpt=pt/ptstep(itypeq)
         iinf=int(dmpt)
         dmpt=dmpt-iinf
         dspt=1.d0-dmpt
         if((dmpt.gt.1).or.(dspt.gt.1).or.(dspt.lt.0)) then
            write(ifmtx,*) 'problem in interpnofpt'
            write(ifmtx,*) pt,iinf
            close(ifmtx)
            stop
         endif
         if(ifshadow) then
            interpnofpt_lin=dspt*nofptAA(itypeq,iinf)+
     &           dmpt*nofptAA(itypeq,iinf+1)
         else
            interpnofpt_lin=dspt*nofpt(itypeq,iinf)+
     &           dmpt*nofpt(itypeq,iinf+1)
         endif
      endif
      return
      end

      subroutine findrootptspec1(itypeq,ifshadow,val,pt0,f0,pt1,f1,eps)
      implicit none
      logical ifshadow
      integer itypeq
      double precision val,pt0,f0,pt1,f1,eps,ptnext,fnext,interpnofpt
      integer ifmtx
      call getMonitorFileIndex(ifmtx)
      if (abs(pt1-pt0).le.eps) return
      do while (abs(pt1-pt0).gt.eps)
         ptnext=(pt0+pt1)/2.d0
         fnext=interpnofpt(itypeq,ptnext,ifshadow)-val
         if(f1*fnext.lt.0.D0) then
            pt0=ptnext
            f0=fnext
         else if(f0*fnext.lt.0.D0) then
            pt1=ptnext
            f1=fnext
         else if(fnext.eq.0.d0) then
            pt0=ptnext
            pt1=ptnext
         else
            write(ifmtx,*) 'no obvious sign change in findrootpt'
            write(ifmtx,*)  pt0,f0,pt1,f1
            close(ifmtx)
            stop            
         endif
      enddo
      return
      end

      subroutine findrootptspec2(itypeq,ifshadow,val,pt0in,f0in,
     &  pt1in,f1in,ptsol,ifail)
      implicit none
      logical ifshadow
      integer ifail,itypeq,nbit,nbitmax
      double precision epsilon,diff,val,pt0,pt1,ptsol,f0,f1,ptnext,
     &  interpnofpt,pt0in,pt1in,f0in,f1in
      parameter (epsilon=1.D-10,nbitmax=50)
      pt0=pt0in
      pt1=pt1in
      f0=f0in
      f1=f1in
      diff=1.d0
      ifail=1
      nbit=0
      do while ((abs(diff).gt.epsilon).and.(nbit.lt.nbitmax))
         nbit=nbit+1
         ptnext=pt1-f1*(pt1-pt0)/(f1-f0)
         if((ptnext.ge.pt1in).or.(ptnext.le.pt0in)) then
            return
         endif
         diff=ptnext-pt1
         f0=f1
         f1=interpnofpt(itypeq,ptnext,ifshadow)-val
         pt0=pt1
         pt1=ptnext
      enddo
      if(abs(diff).lt.epsilon) then
         ifail=0
         ptsol=ptnext
      else
         write(6,*) 'did not converge in ',nbitmax,' iteratsions'
      endif
      return
      end

      double precision function findrootptspec(itypeq,ifshadow,val,
     &  ifail)
      use PBGblckinitptspect2
      implicit none
      logical ifshadow
      integer itypeq,ifail,nbitmax,nbit
      double precision val,ptsol,pt0,f0,pt1,f1,eps,epsstart
      parameter (nbitmax=20,epsstart=1.D-1)
      ifail=1
      pt0=0.d0
      f0=-val
      pt1=ptmax(itypeq)
      f1=1.d0-val
      eps=epsstart
      call findrootptspec1(itypeq,ifshadow,val,pt0,f0,pt1,f1,eps)
      nbit=0
      do while((ifail.ne.0).and.(nbit.lt.nbitmax))
         nbit=nbit+1
         call findrootptspec2(itypeq,ifshadow,val,pt0,f0,pt1,f1,ptsol,
     &        ifail)
         if(ifail.ne.0) then
            eps=eps/8
            call findrootptspec1(itypeq,ifshadow,val,pt0,f0,pt1,f1,eps)
         endif
      enddo
      if(ifail.eq.0) then 
         findrootptspec=ptsol
      else
         write(8,*) 'no solution found in findrootptspec' 
         write(8,*) val,nbit,pt0,pt1,ifail
         findrootptspec=0.D0
      endif
      return
      end

      double precision function get1ptfromspec(itypeq,ifshadow)
      implicit none
      logical ifshadow
      integer itypeq,ifail,nbit,nbitmax
      parameter(nbitmax=5)
      double precision ran2,findrootptspec
      integer ifmtx
      call getMonitorFileIndex(ifmtx)
      ifail=1
      nbit=0
      do while((ifail.eq.1).and.(nbit.lt.nbitmax))
         nbit=nbit+1
         get1ptfromspec=findrootptspec(itypeq,ifshadow,ran2(),ifail)
      enddo
      if(ifail.eq.1) then
         write(ifmtx,*) 'no solution found in get1ptfromspec' 
         close(ifmtx)
         stop
      endif
      return
      end

      subroutine testspectrum
      implicit none
      integer i
      double precision pt,interpnofpt
      open(unit=10,file='testspectrum.res')
      do i=1,200  
         pt=(i-0.5)*0.1d0
         write(10,*) pt,interpnofpt(1,pt,.false.)
      enddo
      close(10)
      end

      subroutine testdndpt
      use PBGblckinitptspect3
      implicit none
      logical ifshadow
      integer tb(100,2),nbit,i,ipt
      parameter(nbit=100000000)
      double precision pt,pmax,stepp,get1ptfromspec,ran2
      real raa(100) 
      parameter (pmax=50,stepp=0.5)
      do i=1,100
         tb(i,1)=0
         tb(i,2)=0
      enddo
      ifshadow=.false.
      do i=1,nbit
         pt= get1ptfromspec(1,ifshadow)
         ipt=int(pt/stepp)+1
         if(ipt.le.100) tb(ipt,1)=tb(ipt,1)+1
      enddo
      ifshadow=.true.
      do i=1,nbit
         if(ran2().lt.globshadow(1)) then
            pt= get1ptfromspec(1,ifshadow)
            ipt=int(pt/stepp)+1
            if(ipt.le.100) tb(ipt,2)=tb(ipt,2)+1
         endif
      enddo
      open(unit=10,file='distribpt.res')
      do ipt=1,100
         if(tb(ipt,1).gt.0) raa(ipt)=real(tb(ipt,2))/tb(ipt,1) 
         write(10,*) (ipt-0.5)*stepp,(tb(ipt,i),i=1,2),raa(ipt)
      enddo
      close(10)
      return
      end

