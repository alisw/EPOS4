
      double precision function sigmadissoc(ires,s)
      use PBGpsiinfo
      use PBGXsection
      implicit none
      integer ires
      double precision s,pi,hbarc
      parameter (pi=3.14159265358979323844d0)
      parameter (hbarc=0.197323)
      if(ires.ne.3) then
        write(6,*) 'not a correct state' 
        sigmadissoc=0.       
        return
      endif
      if(s.ge.4*mc2) then
        sigmadissoc=16**3*Pi/(6*9)*a0c**3*eps0c/hbarc*(s-4*mc2)**(1.5)/
     & (s-mpsi**2)**5*(4*mc2-mpsi**2)**(3.5)
      else
        sigmadissoc=0
      endif
      sigmadissoc=sigmadissoc*Xsectioncranck
      return
      end


      double precision function sigmafus(ires1,ires2,s)
      use PBGpsiinfo
      use PBGXsection
      implicit none
      integer ires1,ires2
      double precision s,pi,hbarc
      parameter (pi=3.14159265358979323844d0)
      parameter (hbarc=0.197323)
      sigmafus=0.d0
      if(((ires1.eq.1).and.(ires2.eq.2)).or.
     &   ((ires1.eq.2).and.(ires2.eq.1))) then
         if(s.ge.4*mc2) then
            sigmafus=4./3.*16**3*pi/(6*9)*a0c**3*eps0c/hbarc*
     &        sqrt(s-4*mc2)/(s*(s-mpsi**2)**3)*(4*mc2-mpsi**2)**3.5
         else
            sigmafus=0.d0
         endif
         sigmafus=sigmafus*Xsectioncranck
         return
      endif 
      return
      end 


      double precision function maxsigmafus(ires1,ires2,threshold)
      implicit none
      integer ires1,ires2
      double precision sigmafus,threshold,sstep,s0,s1,s2,sig0,sig1,sig2,
     &  sintinf,sintsup,sigintinf,sigintsup
      sstep=1.
      s1=threshold**2
      sig1=0.d0
      s0=s1-sstep
      sig0=0.d0
      s2=s1+sstep
      sig2=sigmafus(ires1,ires2,s2)
      do while(sig2.gt.sig1)
        s0=s1
        s1=s2
        sig0=sig1
        sig1=sig2
        s2=s1+sstep
        sig2=sigmafus(ires1,ires2,s2)
      enddo
      if((sig2.gt.sig1).or.(sig0.gt.sig1)) then
        write(8,*) 'problem type 1 in maxsigmafus'
        close(8)
        stop
      endif      
      do while((s2-s1).gt.(sstep/1000.))
         sintinf=(s0+s1)/2.
         sigintinf=sigmafus(ires1,ires2,sintinf)
         if(sigintinf.gt.sig1) then
            s2=s1
            sig2=sig1
            s1=sintinf
            sig1=sigintinf
         else
            sintsup=(s1+s2)/2.
            sigintsup=sigmafus(ires1,ires2,sintsup)
            if(sigintsup.gt.sig1) then
               s0=s1
               sig0=sig1
               s1=sintsup
               sig1=sigintsup
            else
               s0=sintinf
               sig0=sigintinf
               s2=sintsup
               sig2=sigintsup
            endif
         endif
      enddo
      maxsigmafus=sig1
      return
      end
      
      

      subroutine fusion_diff(ires1,ires2,ptot,prelin,min1,min2,mfin1,
     & mfin2,pfin1,pfin2)
      implicit none
      integer k,ires1,ires2
      double precision ptot(0:3),prelin(0:3),min1,min2,mfin1,mfin2,
     & pfin1(0:3),pfin2(0:3),s,ss,ucm(0:3),prelincm(0:3),sr2,
     & costhetafus,
     & costhetacm,vprelincm(3),vprelfincm(3),mfin12,mfin22,pcmfin,
     & pcmfinsq,
     & prelfincm(0:3),prelfin(0:3),bid,chck1,chck2
      if(ires1*ires2.ne.2) then
      write(8,*) 'subroutine fusion_diff just defined for c-cbar fusion'
         close(8)
         stop
      endif         
      s=ptot(0)**2-(ptot(1)**2+ptot(2)**2+ptot(3)**2)
      if(s.lt.(mfin1+mfin2)**2) then
         write(8,*) 'error in fusion_diff: particles cannot react'
         close(8)
         stop
      endif
      ss=sqrt(s)
      do k=0,3
        ucm(k)=ptot(k)/ss
      enddo
      call lorentz2(ucm,prelin,prelincm)
      sr2=(s-(min1+min2)**2)/abs((min1+min2)**2-(mfin1+mfin2)**2)
      costhetacm=costhetafus(sr2)
      vprelincm(1)=prelincm(1)
      vprelincm(2)=prelincm(2)
      vprelincm(3)=prelincm(3)
      mfin12=mfin1**2
      mfin22=mfin2**2
      pcmfinsq=((s-mfin12-mfin22)**2-4.d0*mfin12*mfin22)/(4.d0*s)
      pcmfin=sqrt(pcmfinsq)
      call rotate3d(vprelincm,costhetacm,pcmfin,vprelfincm)
      do k=1,3
         prelfincm(k)=vprelfincm(k)
         ucm(k)=-ucm(k)
      enddo
      prelfincm(0)=0.5*(sqrt(mfin12+pcmfinsq)-sqrt(mfin22+pcmfinsq))
      call lorentz2(ucm,prelfincm,prelfin)
      do k=0,3
         bid=0.5*ptot(k)
         pfin1(k)=bid+prelfin(k)
         pfin2(k)=bid-prelfin(k)
      enddo   
      chck1=pfin1(0)**2-(pfin1(1)**2+pfin1(2)**2+pfin1(3)**2+mfin1**2)
      chck2=pfin2(0)**2-(pfin2(1)**2+pfin2(2)**2+pfin2(3)**2+mfin2**2)
      if((abs(chck1).gt.1.d-8).or.(abs(chck2).gt.1.d-8)) then 
         write(8,*) 'error in fusion_diff: par fin off shell!'
         write(8,*) pfin1(0),pfin1(1),pfin1(2),pfin1(3),mfin1,chck1
         write(8,*) pfin2(0),pfin2(1),pfin2(2),pfin2(3),mfin2,chck2
         close(8)
         stop
      endif
      return
      end


      double precision function costhetafus(sr2)
      implicit none
      integer nit,i
      parameter(nit=15)
      double precision sr2,plow,phigh,ptot,plowrel,ran2,u,w,futry,utry,
     &   wtry,probu,probutry
      plow=128/(3*(1+sr2)**4)
      phigh=(1-1/(1+2*sr2)**3)/6.
      ptot=plow+phigh
      plowrel=plow/ptot
      u=futry(sr2,plowrel,phigh)
      w=probu(sr2,u)/probutry(sr2,u)
      do 10 i=1,nit
         utry=futry(sr2,plowrel,phigh)
         wtry=probu(sr2,utry)/probutry(sr2,utry)
         if(wtry.lt.w)then
            if(w*ran2().gt.wtry) goto 10
         endif
         u=utry
         w=wtry
 10   continue
      costhetafus=sign(u,(1-2*ran2()))
      return
      end
      
      double precision function probutry(sr2,u)
      implicit none
      double precision sr2,u
      probutry=64*(1-u**2)/(1+sr2)**4+sr2/(1+2*sr2*(1-u))**4
      return
      end

      double precision function futry(sr2,plowrel,phigh)
      implicit none
      double precision sr2,plowrel,phigh,ran2,pi
      parameter (pi=3.14159265358979323844d0)
      if(ran2().gt.plowrel) then
         futry=1-((1-6*ran2()*phigh)**(-0.3333333333333333)-1)/(2.*sr2)
      else
         futry=abs(2*cos((4*pi+acos(1-2*ran2()))/3.))
      endif
      return
      end

      double precision function probu(sr2,u)
      implicit none
      double precision sr2,u,v
      v=1-u**2
      probu=(sr2+64*v)/(1+sr2*v)**4
      return
      end

      double precision function probtheta(sr2,th)
      implicit none
      double precision sr2,th,v
      v=sin(th)**2
      probtheta=sin(th)*(sr2+64*v)/(1+sr2*v)**4
      return
      end


