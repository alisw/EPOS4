       double precision function kerneldissocbose(temp,p,mpsi,s)
      implicit none
      double precision hbarc,pi
      parameter (pi=3.14159265358979323844d0)
      parameter (hbarc=0.197323d0)
      double precision temp,p,mpsi,mpsi2,s,epsi
      if (temp.eq.0) then  
        kerneldissocbose=0.d0     
      else
        mpsi2=mpsi**2
        if(s.le.mpsi2) then
          kerneldissocbose=0.d0   
        else
          if(p.eq.0) then
            kerneldissocbose=2*(s-mpsi2)/((mpsi*hbarc)**3*Pi**2*
     &        (exp((s-mpsi2)/(2*mpsi*temp))-1))
          else
            epsi=sqrt(mpsi2+p**2)
            kerneldissocbose=2*temp/(hbarc**3*p*epsi*Pi**2)*
     &       (log(1-exp(-(p+epsi)*(s-mpsi2)/(2*mpsi2*temp)))-
     &        log(1-exp(-(s-mpsi2)/(2*(p+epsi)*temp))))
          endif
        endif
      endif  
      return
      end

       double precision function integrantdissocbose(temp,p,mpsi,s)
       implicit none
      double precision s,temp,p,mpsi,kerneldissocbose,sigmadissoc
      integrantdissocbose=kerneldissocbose(temp,p,mpsi,s)*
     &  (s-mpsi**2)/2*sigmadissoc(3,s)
      return
      end
       
      double precision function integrantdissocbose2(s)
      use PBGforkerneldissoc
      implicit none
      double precision s,kerneldissocbose,sigmadissoc
      integrantdissocbose2=kerneldissocbose(temp2,p2,mass,s)*
     &       (s-mass**2)/2*sigmadissoc(ires,s)
      return
      end


      double precision function decayratebose(temp,p,mpsi,mc)
      use PBGforkerneldissoc
      implicit none
      double precision temp,p,mpsi,mc,decay,liminf,eps,
     &  integrantdissocbose2
      external integrantdissocbose2
      parameter(eps=1.D-2)
      ires=3
      mass=mpsi
      p2=p
      temp2=temp 
      liminf=4*mc**2
      call qtrapimproper(integrantdissocbose2,liminf,decay,eps)
      decayratebose=decay
      return
      end


      double precision function decayratebosemc(temp,ppsi,mpsi,mc)
      implicit none
      double precision temp,ppsi,mpsi,mc,integrantdissocbose
      integer nitmax,nit
      parameter(nitmax=200000)
      double precision mpsi2,mean,sigma,thresh,epsi,beta,sig2,
     & precsiseek,
     & precis,ds,s,m1,m2,w,ran2,cos2rand
      parameter(precsiseek=1.D-2)
      thresh=4*mc**2
      mpsi2=mpsi**2
      epsi=sqrt(ppsi**2+mpsi2)
      if(ppsi.lt.temp) then
         beta=3
      else
         beta=4
      endif
      sig2=mpsi2/(mpsi2/(temp*(epsi+ppsi))+2*beta*(thresh/mpsi2-1.))
      nit=0
      m1=0
      m2=0
      do while((nit.lt.10).or.((nit.lt.nitmax)
     &           .and.(precis.gt.precsiseek)))
         nit=nit+1
         ds=-2*sig2*(log(ran2())+log(ran2())+log(ran2())*cos2rand())
         s=ds+thresh
         w=integrantdissocbose(temp,ppsi,mpsi,s)/
     &          (exp(-ds/(2*sig2))*sqrt(ds)*ds)
         m1=m1+w
         mean=m1/nit
         m2=m2+w**2
         sigma=sqrt(abs(m2/nit-mean**2))
         precis=sigma/(mean*sqrt(1.*nit))
      enddo   
      decayratebosemc=mean*3*sqrt(2*3.14159265358979323844*sig2**5)
      return
      end


      double precision function decayratejpsi(temp,densener,p)
      use PBGblcktransi
      use PBGpsiinfo
      use PBGXsection
      use PBGfordecaydat1
      use PBGfordecaydat2
      implicit none
      integer ipinf,ipsup,itinf,itsup
      double precision temp,densener,p,decayratebosemc,dp,dtemp,prap,
     & trap
      if(densener.lt.entrans_min) then
         decayratejpsi=0.
         return
      endif
      if((p.le.pmax).and.(temp.ge.tempmin).and.(temp.le.tempmax)) then
         if(temp.eq.tempmax) then
            if(p.eq.pmax)then
               decayratejpsi=decaygrid(nbp,nbtemp)
               return
            else
               prap=p/deltp
               ipinf=int(prap)
               ipsup=ipinf+1
               dp=prap-ipinf
               decayratejpsi=dp*decaygrid(ipsup,nbtemp)+
     &            (1-dp)*decaygrid(ipinf,nbtemp)
               return
            endif
         else 
            if(p.eq.pmax) then
               trap=(temp-tempmin)/delttemp
               itinf=int(trap)
               itsup=itinf+1
               dtemp=trap-itinf
               decayratejpsi=dtemp*decaygrid(nbp,itsup)+
     &            (1-dtemp)*decaygrid(nbp,itinf)            
               return
            else
               trap=(temp-tempmin)/delttemp
               itinf=int(trap)
               itsup=itinf+1
               dtemp=trap-itinf
               prap=p/deltp
               ipinf=int(prap)
               ipsup=ipinf+1
               dp=prap-ipinf
               decayratejpsi=
     &           dp*(dtemp*decaygrid(ipsup,itsup)+
     &               (1-dtemp)*decaygrid(ipsup,itinf))+
     &           (1-dp)*(dtemp*decaygrid(ipinf,itsup)+
     &               (1-dtemp)*decaygrid(ipinf,itinf))
               return  
            endif
         endif
      else
         decayratejpsi=Xsectioncranck*decayratebosemc(temp,p,mpsi,mc)
      endif
      return
      end


      double precision function covdecayratejpsi2(t)
      use PBGpsiinfo
      use PBGfordecaydat1
      use PBGforcovdecay2
      use PBGinfotr
      implicit none
      integer stathyd,imed
      double precision t,x(0:3),u(0:3),dtau,beta(3),epsiincellframe,
     & ppsiincellframe2,decayratejpsi,ppsiincellframe,densener
      dtau=(t-x0(0))/p(0)
      x(1)=x0(1)+p(1)*dtau
      x(2)=x0(2)+p(2)*dtau
      x(3)=x0(3)+p(3)*dtau
      x(0)=t
      covdecayratejpsi2=0.d0
      call givelocalhydro(x,temp,beta,densener,stathyd,imed)
      if(stathyd.ne.0) return
      u(0)=1./sqrt(1-(beta(1)**2+beta(2)**2+beta(3)**2))
      u(1)=u(0)*beta(1)
      u(2)=u(0)*beta(2)
      u(3)=u(0)*beta(3)
      epsiincellframe=u(0)*p(0)-(u(1)*p(1)+u(2)*p(2)+u(3)*p(3))
      ppsiincellframe2=epsiincellframe**2-mpsi**2
      if(ppsiincellframe2.lt.-1.d-10) then
         write(8,*) 'error in covdecayratejpsi2:',ppsiincellframe2
         close(8)
         stop
      endif   
      ppsiincellframe=sqrt(ppsiincellframe2)
      if(ppsiincellframe.gt.pmax) then
         if((ifail.eq.0).and.(infotrack.ge.5)) then
            write(8,*)'large u_hyd.p_psi encountered in covdecayratejp',
     &         epsiincellframe
            write(8,*) 'at point ',x(0),x(1),x(2),x(3),' one has:'
            write(8,*) 'temp= ',temp
            write(8,*) 'u_cell= ',u(0),u(1),u(2),u(3)
            write(8,*) 'and p_psi= ',p(0),p(1),p(2),p(3)
         endif
         ifail=ifail+1   
      endif
      covdecayratejpsi2=decayratejpsi(temp,densener,ppsiincellframe)*
     &        epsiincellframe/p(0)
      return
      end


      double precision function survivalprob(x0e,pe,tin,tfin)
      use PBGforcovdecay2
      use PBGinfotr
      implicit none
      integer i
      double precision x0e(0:3),pe(0:3),tin,tfin,decay,precis
      parameter(precis=1.D-3)
      external covdecayratejpsi2
      if(tfin.lt.tin) then
         write(8,*) 'tfin<tin in dissocprob'
         close(8)
         stop
      else
         if(tfin.eq.tin) then
            survivalprob=1.
            return
         else
            do i=0,3 
              x0(i)=x0e(i)
              p(i)=pe(i)
            enddo  
            ifail=0
            call qromb(covdecayratejpsi2,tin,tfin,decay,precis)
            survivalprob=exp(-decay)
            if((ifail.ne.0).and.(infotrack.ge.4)) then 
               write(8,*) 'needed ',ifail,
     &           ' sup. eval of decayrate to eval this surv. prob.' 
            endif
         endif
      endif   
      return
      end


      subroutine melt(ipart)
      use PBGgenvar
      implicit none
      integer ipart,ires,i,iresnew
      double precision rap
      ires=partinfo_ires(ipart)
      if((ires.eq.1).or.(ires.eq.2).or.(ires.eq.14).or.(ires.eq.15))then
         write(8,*) 'do not know to melt quarks ;-)'
         close(8)
         stop
      endif
      if((ires.eq.3).or.(ires.eq.16)) then
         return
      endif
       if((ires.eq.4).or.(ires.eq.6).or.(ires.eq.8).or.
     & (ires.eq.10).or.(ires.eq.12))then
         iresnew=1
         partinfo_ires(ipart)=iresnew
         rap=resinfo_mass(iresnew)/resinfo_mass(ires)
         do i=0,3
            partinfo_p(ipart,i)=rap*partinfo_p(ipart,i)
         enddo
         partinfo_mass(ipart)=resinfo_mass(iresnew)
        elseif((ires.eq.5).or.(ires.eq.7).or.(ires.eq.9).or.
     & (ires.eq.11).or.(ires.eq.13))then
         iresnew=2
         partinfo_ires(ipart)=iresnew
         rap=resinfo_mass(iresnew)/resinfo_mass(ires)
         do i=0,3
            partinfo_p(ipart,i)=rap*partinfo_p(ipart,i)
         enddo
         partinfo_mass(ipart)=resinfo_mass(iresnew)

         elseif((ires.eq.17).or.(ires.eq.19).or.(ires.eq.21).or.
     & (ires.eq.23).or.(ires.eq.25))then
         iresnew=14
         partinfo_ires(ipart)=iresnew
         rap=resinfo_mass(iresnew)/resinfo_mass(ires)
         do i=0,3
            partinfo_p(ipart,i)=rap*partinfo_p(ipart,i)
         enddo
         partinfo_mass(ipart)=resinfo_mass(iresnew)

         elseif((ires.eq.18).or.(ires.eq.20).or.(ires.eq.22).or.
     & (ires.eq.24).or.(ires.eq.26))then
         iresnew=15
         partinfo_ires(ipart)=iresnew
         rap=resinfo_mass(iresnew)/resinfo_mass(ires)
         do i=0,3
            partinfo_p(ipart,i)=rap*partinfo_p(ipart,i)
         enddo
         partinfo_mass(ipart)=resinfo_mass(iresnew)
         else 
           write(8,*) 'wrong ires in subroutine melt'
           close(8)
           stop        
       endif

      if(ires.gt.26) then
         write(8,*) 'undefined ires in melt'
         close(8)
         stop
      endif
      return
      end


      subroutine hadronize(ipart,itypH,betafluid)
      use PBGfragandcoal
      use PBGgenvar
      implicit none
      logical ifcoal
      integer ipart,ires,iresh,itypH,i
      double precision xQ(0:3),pQ(0:3),pmes(0:3),betafluid(3),
     & pb_coal
      ires=partinfo_ires(ipart)
      if(.not.(resinfo_ifquark(ires))) then
         write(8,*) 'in Hadronize: should be just for quarks'
         close(8)
         stop
      endif
      do i=0,3
         xQ(i)=partinfo_r(ipart,i)
         pQ(i)=partinfo_p(ipart,i)
      enddo
      call hadronizebasicnewall(xQ,pQ,ires,pmes,iresh,itypH,betafluid,
     & 10, pb_coal,ifcoal,.true.)

      partinfo_ires(ipart)=iresh
      do i=0,3
         partinfo_pghost(ipart,i)=pQ(i)
         partinfo_rghost(ipart,i)=xQ(i)
         partinfo_p(ipart,i)=pmes(i)
      enddo
       partinfo_mass(ipart)=resinfo_mass(iresh)
      return
      end



      subroutine hadronizebasicnewall(xq,pq,iresq,ph,iresh,itypH,betafl,
     &  nbcogen,pb_coal,ifcoal,ifsavehadroloc)
      use PBGfragandcoal
      use PBGgenvar
      use PBGfordisplay
      implicit none

      logical ifcoal,ifsavehadroloc,ifsave
      integer iresq,iresh,itypH,nbcogen,icogen,ires_const,tl,i
      integer ikineticfragtype,ipt,ifail
      parameter (ikineticfragtype=2)
      double precision xq(0:3),pq(0:3),ph(0:3),betafl(3),pqnorminfl,
     & pb_coal,ufl(0:3),pqinfl(0:3),phinfl(0:3),ov,pi,
     & pconstnorm2,pconstnorm,costh,pconsttr,phi,etherm,pconstinfl(0:3),
     & pconst(0:3),dsigunit(0:3),uQdsig,uQucell,ucelldsig,d,
     & dsigcov(0:3),dsigcovinfl(0:3),pqnorm,pqt(0:3),mtquark,mtmeson,
     & uflback(0:3),plightqinfl(0:3),plight(0:3),mass2Qq,ymeson,
     & pQq(0:3),ayHQ,probtotp,pb_coalH(5),pcoalD0,pcoalDs,pcoalLc,
     & pcoalXc,randhfH,randprobfrag,phinflh(0:3)
      parameter(pi=3.14159265358979323844d0)
      double precision ran2
 10   format(50(E12.5,' '))

      ifail=0
      ifsave=ifsavehadroloc.and.ifsavehadro
      ufl(0)=1/sqrt(1-(betafl(1)**2+betafl(2)**2+betafl(3)**2))
      ufl(1)=ufl(0)*betafl(1)
      ufl(2)=ufl(0)*betafl(2)
      ufl(3)=ufl(0)*betafl(3)
      call lorentz2(ufl,pq,pqinfl)

      pqnorminfl=sqrt(pqinfl(1)**2+pqinfl(2)**2+pqinfl(3)**2)
      if((itypcoal.ge.2).and.((itypH.eq.2).or.(itypH.eq.3))) then
         uQucell=pqinfl(0)/resinfo_mass(iresq)
         call getnormal(xq,dsigunit,tl,.false.)
         uQdsig=dsigunit(0)*pq(0)+dsigunit(1)*pq(1)+dsigunit(2)*pq(2)+
     &             dsigunit(3)*pq(3)
         uQdsig=uQdsig/resinfo_mass(iresq)
         if(uQdsig.lt.0) then
            call getnormal(xq,dsigunit,tl,.true.)
            ifail=1
         endif
         ucelldsig=ufl(0)*dsigunit(0)+ufl(1)*dsigunit(1)+ufl(2)*
     &        dsigunit(2)+ufl(3)*dsigunit(3)
      endif

      ifcoal=.false.
      pb_coal=0.
      if(itypH.eq.1) then
         pqnorm=sqrt(pq(1)**2+pq(2)**2+pq(3)**2)
         if((itypfrag.lt.3).and.(pqnorm.le.pcrithadr(itypquark(iresq)))) 
     &   then
            ifcoal=.true.
            pb_coal=1.
         endif
      endif

      if(itypH.eq.2) then 
         ifcoal=.true.
      elseif(itypH.eq.3) then
         pqnorm=sqrt(pq(1)**2+pq(2)**2+pq(3)**2)
         if(((itypfrag.eq.1).or.(itypfrag.eq.2)).and.
     &       (pqnorm.lt.pcrithadr(itypquark(iresq)))) then
            ifcoal=.true.
         else
            if(itypcoal.eq.1.or.itypcoal.eq.5) then
               call gimmeprobcoal(itypquark(iresq),pqnorminfl,pb_coal)
            else
               call gimmeprobcoalbis(itypquark(iresq),pq,
     &              tl,pqnorminfl,uQucell,uQdsig,ucelldsig,d,pb_coal)
               if(ifsave.and.(nbhadro.le.1000)) then
                  write(12,10) xq,pq,sqrt(xq(0)**2-xq(3)**2),
     &              sqrt(pq(1)**2+pq(2)**2),ufl,pqnorminfl,dsigunit,
     &              real(tl),
     &              uQucell,uQdsig,ucelldsig,d,pb_coal
               endif

               if(isnan(pb_coal)) then
                  write(8,*) 'NAN pb_coal found in hadronizebasicnew'
                  write(8,10) xq,pq,sqrt(xq(0)**2-xq(3)**2),
     &              sqrt(pq(1)**2+pq(2)**2),ufl,pqnorminfl,dsigunit,
     &              real(tl),
     &              uQucell,uQdsig,ucelldsig,d,pb_coal
                  close(8)
                  stop
               endif
            endif
            if(pb_coal.ge.1) then
               ifcoal=.true.
            else
               ifcoal=(ran2().lt.pb_coal)
            endif
         endif
      endif

      ayHQ=log((pq(0)+abs(pq(3)))/
     &          sqrt(pq(1)**2+pq(2)**2+resinfo_mass(iresq)**2))
      if((ifsave).and.((iresq.eq.1).or.(iresq.eq.2)).and.
     &     (ayHQ.le.rapsupobs)) then
         nbhadro=nbhadro+1
         ipt=int(sqrt(pq(1)**2+pq(2)**2)/0.2)+1
         if(ipt.le.200) then
            nbhadroofpt(ipt,0)=nbhadroofpt(ipt,0)+1.
            if(itypH.eq.2) then
               nbhadroofpt(ipt,tl)=nbhadroofpt(ipt,tl)+1.
            else if(itypH.eq.3) then
               nbhadroofpt(ipt,tl)=nbhadroofpt(ipt,tl)+1.
               avprobcoalofpt(ipt,tl)=avprobcoalofpt(ipt,tl)+pb_coal
               varprobcoalofpt(ipt,tl)=varprobcoalofpt(ipt,tl)+
     &              pb_coal**2
               if((pb_coal.lt.5.).and.(ipt.le.5)) then
                  ipt=int(pb_coal/0.1)+1
                  distribproba(ipt,tl)=distribproba(ipt,tl)+1.
               endif
            endif
         endif
      endif
      if(ifcoal) then
         uflback(0)=ufl(0)
         uflback(1)=-ufl(1)
         uflback(2)=-ufl(2)
         uflback(3)=-ufl(3)
         if(itypcoal.eq.1) then
            iresh=iresq+3
            call coalesce(pqinfl,iresq,phinfl,iresh)
            call lorentz2(uflback,phinfl,ph)
         elseif(itypcoal.eq.5) then
            call gimmeprobcoalHad(itypquark(iresq),pqnorminfl,pb_coalH)

            probtotp=pb_coalH(1)+pb_coalH(2)+pb_coalH(3)+
     & pb_coalH(4)+pb_coalH(5)
            pcoalD0=(pb_coalH(1))/probtotp
            pcoalDs=(pb_coalH(1)+pb_coalH(2))/probtotp
            pcoalLc=(pb_coalH(1)+pb_coalH(2)+pb_coalH(3))/probtotp
            pcoalXc=(pb_coalH(1)+pb_coalH(2)+pb_coalH(3)+
     & pb_coalH(4))/probtotp

            randhfH = ran2()
            if(randhfH.le.pcoalD0)then
            iresh=iresq+3
            elseif(randhfH.gt.pcoalD0.and.randhfH.le.pcoalDs)then
            iresh=iresq+5
            elseif(randhfH.gt.pcoalDs.and.randhfH.le.pcoalLc)then
            iresh=iresq+7
            elseif(randhfH.gt.pcoalLc.and.randhfH.le.pcoalXc)then
            iresh=iresq+9
            elseif(randhfH.gt.pcoalXc)then
            iresh=iresq+11
            endif
           call coalesceall(pqinfl,phinflh,iresh)            
            phinflh(0)=sqrt(resinfo_mass(iresh)**2+phinflh(1)**2+
     &   phinflh(2)**2+phinflh(3)**2)
            call lorentz2(uflback,phinflh,ph)
         else
            iresh=iresq+3
            dsigcov(0)=dsigunit(0)
            dsigcov(1)=-dsigunit(1)
            dsigcov(2)=-dsigunit(2)
            dsigcov(3)=-dsigunit(3)
            call lorentz2(ufl,dsigcov,dsigcovinfl)
            call coalescebis(pqinfl,dsigcovinfl,iresq,tl,uQucell,uQdsig,
     &              ucelldsig,plightqinfl,mass2Qq,phinfl,iresh)
            if(ifcoalinflcell) then
               call lorentz2(uflback,phinfl,ph)
            else
               call lorentz2(uflback,plightqinfl,plight)
               do i=0,3
                  pQq(i)=pq(i)+plight(i)
               enddo
               mtmeson=sqrt(mass2Qq+pQq(1)**2+pQq(2)**2)
               if(pQq(3).gt.0) then
                  ymeson=log((pQq(0)+pQq(3))/mtmeson)
               else
                  ymeson=-log((pQq(0)-pQq(3))/mtmeson)                     
               endif
               ph(1)=pQq(1)
               ph(2)=pQq(2)
               mtmeson=sqrt(resinfo_mass(iresh)**2+ph(1)**2+ph(2)**2)
               ph(3)=mtmeson*sinh(ymeson)
               ph(0)=mtmeson*cosh(ymeson)
            endif
         endif
         if((ifsave).and.((iresq.eq.1).or.(iresq.eq.2)).and.
     &     (ayHQ.le.rapsupobs)) then
            ipt=int(sqrt(ph(1)**2+ph(2)**2)/0.2)+1
            if(ipt.le.200) then
               nbmesonsofpthadrocoal(ipt,0)=nbmesonsofpthadrocoal(ipt,0)
     &              +1.
               nbmesonsofpthadrocoal(ipt,tl)=nbmesonsofpthadrocoal(ipt,
     &              tl)+1.
            endif
         endif
      else
         randprobfrag = ran2()
         if(randprobfrag.le.0.848)then
           iresh=iresq+3
         elseif(randprobfrag.gt.0.848.and.randprobfrag.le.0.928)then
           iresh=iresq+5
         elseif(randprobfrag.gt.0.928.and.randprobfrag.le.0.988)then
           iresh=iresq+7
         elseif(randprobfrag.gt.0.988.and.randprobfrag.le.0.998)then
           iresh=iresq+9
         elseif(randprobfrag.gt.0.998)then
           iresh=iresq+11
         endif
         if(ikineticfragtype.eq.1) then
            call fragment(pq,iresq,ph,iresh)
         else
            mtquark=sqrt(pq(1)**2+pq(2)**2+resinfo_mass(iresq)**2)
            pqt(0)=mtquark
            pqt(1)=pq(1)
            pqt(2)=pq(2)
            pqt(3)=0.d0
            call fragment(pqt,iresq,ph,iresh)
            mtmeson=ph(0)
            ph(3)=mtmeson*pq(3)/mtquark
            ph(0)=mtmeson*sqrt((pq(3)/mtquark)**2+1)
         endif
         if((ph(1)**2+ph(2)**2.gt.pq(1)**2+pq(2)**2)) then
            write(6,*) 'trans acceleration'
            write(6,*) ph(1)**2+ph(2)**2,pq(1)**2+pq(2)**2
         endif
         if((ifsave).and.((iresq.eq.1).or.(iresq.eq.2)).and.
     &     (ayHQ.le.rapsupobs)) then
            ipt=int(sqrt(ph(1)**2+ph(2)**2)/0.2)+1
            if(ipt.le.200) then
               nbmesonsofpthadrofrag(ipt,0)=nbmesonsofpthadrofrag(ipt,0)
     &              +1.
               if(itypH.eq.3) then
                  if(abs(tl).ne.1) then
                     write(8,*) '|tl|=1 expected'
                     write(8,*) xq
                     close(8)
                  endif
                  nbmesonsofpthadrofrag(ipt,tl)=
     &                 nbmesonsofpthadrofrag(ipt,tl)+1.
               endif
            endif
         endif
      endif
      if(nbcogen.ge.1) ov=1.d0/nbcogen
      ires_const=0
      do icogen=1,nbcogen
         call getetherm(mqconstit,tempcoal,etherm)
         pconstnorm2=etherm**2-mqconstit2
         pconstnorm=sqrt(pconstnorm2)
         costh=2*ran2()-1.d0
         pconsttr=sqrt(1-costh**2)*pconstnorm
         phi=2*pi*ran2()
         pconstinfl(0)=etherm
         pconstinfl(1)=pconsttr*cos(phi)
         pconstinfl(2)=pconsttr*sin(phi)
         pconstinfl(3)=pconstnorm*costh
         call lorentz2(ufl,pconstinfl,pconst)
         call pack_it_1part(ires_const,pconst,xq,ov)
      enddo

      return
      end


      subroutine hadronizebasicnew(xq,pq,iresq,ph,iresh,itypH,betafl,
     &  nbcogen,pb_coal,ifcoal,ifsavehadroloc)
      use PBGfragandcoal
      use PBGgenvar
      use PBGfordisplay
      implicit none
      logical ifcoal,ifsavehadroloc,ifsave
      integer iresq,iresh,itypH,nbcogen,icogen,ires_const,tl,i
      integer ikineticfragtype,ipt,ifail
      parameter (ikineticfragtype=2)
      double precision xq(0:3),pq(0:3),ph(0:3),betafl(3),pqnorminfl,
     & pb_coal,pb_coal0,ufl(0:3),pqinfl(0:3),phinfl(0:3),ov,
     & pconstnorm2,pconstnorm,costh,pconsttr,phi,etherm,pconstinfl(0:3),
     & pconst(0:3),dsigunit(0:3),uQdsig,uQucell,ucelldsig,d,
     & dsigcov(0:3),dsigcovinfl(0:3),pqnorm,pqt(0:3),mtquark,mtmeson,
     & uflback(0:3),plightqinfl(0:3),plight(0:3),mass2Qq,ymeson,
     & pQq(0:3),ayHQ,pi
      double precision ran2
      integer ifmtx

      parameter (pi=3.14159265358979323844d0)
      call getMonitorFileIndex(ifmtx)
 10   format(50(E12.5,' '))
      ifail=0
      ifsave=ifsavehadroloc.and.ifsavehadro
      ufl(0)=1/sqrt(1-(betafl(1)**2+betafl(2)**2+betafl(3)**2))
      ufl(1)=ufl(0)*betafl(1)
      ufl(2)=ufl(0)*betafl(2)
      ufl(3)=ufl(0)*betafl(3)
      call lorentz2(ufl,pq,pqinfl)
      pqnorminfl=sqrt(pqinfl(1)**2+pqinfl(2)**2+pqinfl(3)**2)
      if((itypcoal.ge.2).and.((itypH.eq.2).or.(itypH.eq.3))) then
         uQucell=pqinfl(0)/resinfo_mass(iresq)
         call getnormal(xq,dsigunit,tl,.false.)
         uQdsig=dsigunit(0)*pq(0)+dsigunit(1)*pq(1)+dsigunit(2)*pq(2)+
     &             dsigunit(3)*pq(3)
         uQdsig=uQdsig/resinfo_mass(iresq)
         if(uQdsig.lt.0) then
            call getnormal(xq,dsigunit,tl,.true.)
            ifail=1
         endif
         ucelldsig=ufl(0)*dsigunit(0)+ufl(1)*dsigunit(1)+ufl(2)*
     &        dsigunit(2)+ufl(3)*dsigunit(3)
      endif
      ifcoal=.false.
      pb_coal=0.
      if(itypH.eq.1) then
         pqnorm=sqrt(pq(1)**2+pq(2)**2+pq(3)**2)
         if((itypfrag.lt.3).and.(pqnorm.le.pcrithadr(itypquark(iresq)))) 
     &   then
            ifcoal=.true.
            pb_coal=1.
            call getnormal(xq,dsigunit,tl,.false.)
         endif
      endif
      if(itypH.eq.2) ifcoal=.true.
      if(itypH.eq.3) then 
         pqnorm=sqrt(pq(1)**2+pq(2)**2+pq(3)**2)
         if(((itypfrag.eq.1).or.(itypfrag.eq.2)).and.
     &       (pqnorm.lt.pcrithadr(itypquark(iresq)))) then
            ifcoal=.true.
         else
            if(itypcoal.eq.1) then
               call gimmeprobcoal(itypquark(iresq),pqnorminfl,pb_coal)
               call gimmeprobcoal(itypquark(iresq),pcrithadr(itypquark(
     &              iresq)),pb_coal0)
               pb_coal=pb_coal/pb_coal0               
            else
               call gimmeprobcoalbis(itypquark(iresq),pq,
     &              tl,pqnorminfl,uQucell,uQdsig,ucelldsig,d,pb_coal)
               if(ifsave.and.(nbhadro.le.1000)) then
                  write(12,10) xq,pq,sqrt(xq(0)**2-xq(3)**2),
     &              sqrt(pq(1)**2+pq(2)**2),ufl,pqnorminfl,dsigunit,
     &              real(tl),
     &              uQucell,uQdsig,ucelldsig,d,pb_coal
               endif
               if(isnan(pb_coal)) then
                  write(ifmtx,*) 'NAN pb_coal in hadronizebasicnew'
                  write(ifmtx,10) xq,pq,sqrt(xq(0)**2-xq(3)**2),
     &              sqrt(pq(1)**2+pq(2)**2),ufl,pqnorminfl,dsigunit,
     &              real(tl),
     &              uQucell,uQdsig,ucelldsig,d,pb_coal
                  close(ifmtx)
                  stop
               endif
            endif
            if(pb_coal.ge.1) then
               ifcoal=.true.
            else
               ifcoal=(ran2().lt.pb_coal)
            endif
         endif
      endif
      ayHQ=log((pq(0)+abs(pq(3)))/
     &          sqrt(pq(1)**2+pq(2)**2+resinfo_mass(iresq)**2))
      if((ifsave).and.((iresq.eq.1).or.(iresq.eq.2)).and.
     &     (ayHQ.le.rapsupobs)) then
         nbhadro=nbhadro+1
         ipt=int(sqrt(pq(1)**2+pq(2)**2)/0.2)+1
         if(ipt.le.200) then
            nbhadroofpt(ipt,0)=nbhadroofpt(ipt,0)+1.
            if(itypH.eq.2) then
               nbhadroofpt(ipt,tl)=nbhadroofpt(ipt,tl)+1.
            else if(itypH.eq.3) then
               nbhadroofpt(ipt,tl)=nbhadroofpt(ipt,tl)+1.
               avprobcoalofpt(ipt,tl)=avprobcoalofpt(ipt,tl)+pb_coal
               varprobcoalofpt(ipt,tl)=varprobcoalofpt(ipt,tl)+
     &              pb_coal**2 
               if((pb_coal.lt.5.).and.(ipt.le.5)) then
                  ipt=int(pb_coal/0.1)+1
                  distribproba(ipt,tl)=distribproba(ipt,tl)+1.
               endif
            endif
         endif
      endif
      if(ifcoal) then
         uflback(0)=ufl(0)
         uflback(1)=-ufl(1)
         uflback(2)=-ufl(2)
         uflback(3)=-ufl(3)
         if(itypcoal.eq.1) then  
            call coalesce(pqinfl,iresq,phinfl,iresh)
            call lorentz2(uflback,phinfl,ph)
         else
            dsigcov(0)=dsigunit(0)
            dsigcov(1)=-dsigunit(1)
            dsigcov(2)=-dsigunit(2)
            dsigcov(3)=-dsigunit(3)
            call lorentz2(ufl,dsigcov,dsigcovinfl)
            call coalescebis(pqinfl,dsigcovinfl,iresq,tl,uQucell,uQdsig,
     &              ucelldsig,plightqinfl,mass2Qq,phinfl,iresh)
            if(ifcoalinflcell) then
               call lorentz2(uflback,phinfl,ph)
            else
               call lorentz2(uflback,plightqinfl,plight)
               do i=0,3
                  pQq(i)=pq(i)+plight(i)
               enddo
               mtmeson=sqrt(mass2Qq+pQq(1)**2+pQq(2)**2)
               if(pQq(3).gt.0) then
                  ymeson=log((pQq(0)+pQq(3))/mtmeson)
               else
                  ymeson=-log((pQq(0)-pQq(3))/mtmeson)                     
               endif
               ph(1)=pQq(1)
               ph(2)=pQq(2)
               mtmeson=sqrt(resinfo_mass(iresh)**2+ph(1)**2+ph(2)**2)
               ph(3)=mtmeson*sinh(ymeson)
               ph(0)=mtmeson*cosh(ymeson)
            endif
         endif
         if((ifsave).and.((iresq.eq.1).or.(iresq.eq.2)).and.
     &     (ayHQ.le.rapsupobs)) then
            ipt=int(sqrt(ph(1)**2+ph(2)**2)/0.2)+1
            if(ipt.le.200) then
               nbmesonsofpthadrocoal(ipt,0)=nbmesonsofpthadrocoal(ipt,0)
     &              +1.
               nbmesonsofpthadrocoal(ipt,tl)=nbmesonsofpthadrocoal(ipt,
     &              tl)+1.
            endif
         endif
      else
         if(ikineticfragtype.eq.1) then
            call fragment(pq,iresq,ph,iresh)
         else         
            mtquark=sqrt(pq(1)**2+pq(2)**2+resinfo_mass(iresq)**2)
            pqt(0)=mtquark
            pqt(1)=pq(1)
            pqt(2)=pq(2)
            pqt(3)=0.d0
            call fragment(pqt,iresq,ph,iresh)
            mtmeson=ph(0)
            ph(3)=mtmeson*pq(3)/mtquark
            ph(0)=mtmeson*sqrt((pq(3)/mtquark)**2+1)
         endif
         if((ph(1)**2+ph(2)**2.gt.pq(1)**2+pq(2)**2)) then
            write(6,*) 'trans acceleration'
            write(6,*) ph(1)**2+ph(2)**2,pq(1)**2+pq(2)**2
         endif
         if((ifsave).and.((iresq.eq.1).or.(iresq.eq.2)).and.
     &     (ayHQ.le.rapsupobs)) then
            ipt=int(sqrt(ph(1)**2+ph(2)**2)/0.2)+1
            if(ipt.le.200) then
               nbmesonsofpthadrofrag(ipt,0)=nbmesonsofpthadrofrag(ipt,0)
     &              +1.
               if(itypH.eq.3) then
                  if(abs(tl).ne.1) then
                     write(ifmtx,*) '|tl|=1 expected'
                     write(ifmtx,*) xq
                     close(ifmtx)
                     stop
                  endif
                  nbmesonsofpthadrofrag(ipt,tl)=
     &                 nbmesonsofpthadrofrag(ipt,tl)+1.
               endif
            endif
         endif
      endif
      if(nbcogen.ge.1) ov=1.d0/nbcogen
      ires_const=0
      do icogen=1,nbcogen
         call getetherm(mqconstit,tempcoal,etherm)
         pconstnorm2=etherm**2-mqconstit2
         pconstnorm=sqrt(pconstnorm2)
         costh=2*ran2()-1.d0
         pconsttr=sqrt(1-costh**2)*pconstnorm
         phi=2*pi*ran2()
         pconstinfl(0)=etherm
         pconstinfl(1)=pconsttr*cos(phi)
         pconstinfl(2)=pconsttr*sin(phi)
         pconstinfl(3)=pconstnorm*costh
         call lorentz2(ufl,pconstinfl,pconst)
         call pack_it_1part(ires_const,pconst,xq,ov)
      enddo
      return
      end



      subroutine getnormalnum(x0,dsigunit,tl,verbose)
      use PBGhydro
      implicit none
      logical verbose
      integer k,ifail,tl,imedp,imedm
      double precision x0(0:3),derdense(0:3),xprime(0:3),tb,temp,
     & beta(3),densep,densem,dx,dtb,ch,sh,dt,dz,derdensedtb,normgrad,
     & dsigunit(0:3)
      parameter(dx=0.05)
      if(verbose) write(6,*) 'entering getnormalnum'
      tb=sqrt(x0(0)**2-x0(3)**2)
      ch=x0(0)/tb
      sh=x0(3)/tb
      do k=1,2
         x0(k)=x0(k)+dx
         call givelocalhydro(x0,temp,beta,densep,ifail,imedp)
         x0(k)=x0(k)-2*dx
         call givelocalhydro(x0,temp,beta,densem,ifail,imedm)
         derdense(k)=(densep-densem)/(2*dx)
         x0(k)=x0(k)+dx
         if(verbose) write(6,*) x0(k),densep,imedp,densem,imedm,
     &        derdense(k)          
         xprime(k)=x0(k)
      enddo
      if(tb.gt.2*tbmin) then
         dtb=dx
      else
         dtb=(tb-tbmin)*dx
      endif
      dt=dtb*ch
      dz=dtb*sh
      xprime(0)=x0(0)+dt
      xprime(3)=x0(3)+dz
      call givelocalhydro(xprime,temp,beta,densep,ifail,imedp)
      xprime(0)=x0(0)-dt
      xprime(3)=x0(3)-dz
      call givelocalhydro(xprime,temp,beta,densem,ifail,imedm)
      derdensedtb=(densep-densem)/(2*dtb)      
      if(verbose) write(6,*) dtb,ch,densep,imedp,densem,imedm,
     &     derdensedtb          
      derdense(0)=derdensedtb*ch
      derdense(3)=-derdensedtb*sh
      normgrad=derdensedtb**2-(derdense(1)**2+derdense(2)**2)
      if(normgrad.gt.0) then
         tl=+1
      else
         tl=-1
      endif
      normgrad=sqrt(abs(normgrad))
      do k=0,3
         dsigunit(k)=-derdense(k)/normgrad
      enddo
      return
      end


      subroutine getnormal(x0,dsigunit,tl,verbose)
      use PBGsceplas
      implicit none
      logical verbose
      integer tl
      double precision x0(0:3),dsigunit(0:3)
      if(isce.eq.7) then
         call getnormalnotboostinv(x0,dsigunit,tl)
      else
         call getnormalnum(x0,dsigunit,tl,verbose)
      endif
      end
      

      subroutine getnormalnotboostinv(x0,dsigunit,tl)
      use PBGhydro
      implicit none
      integer k,ifail,tl,imed
      double precision x0(0:3),derdense(0:3),xprime(0:3),tb,temp,
     & beta(3),densep,densem,dx,dtb,ch,sh,dt,dz,derdensedtb,normgrad,
     & dsigunit(0:3),eta,deta,derdensedeta,tb2
      parameter(dx=0.05,deta=0.01)
      tb2=x0(0)**2-x0(3)**2
      tb=sqrt(tb2)
      eta=0.5d0*log((x0(0)+x0(3))/(x0(0)-x0(3)))
      ch=x0(0)/tb
      sh=x0(3)/tb
      do k=1,2
         x0(k)=x0(k)+dx
         call givelocalhydro(x0,temp,beta,densep,ifail,imed)
         x0(k)=x0(k)-2*dx
         call givelocalhydro(x0,temp,beta,densem,ifail,imed)
         derdense(k)=(densep-densem)/(2*dx)
         x0(k)=x0(k)+dx
         xprime(k)=x0(k)
      enddo
      if(tb.gt.2*tbmin) then
         dtb=dx
      else
         dtb=(tb-tbmin)*dx
      endif
      dt=dtb*ch
      dz=dtb*sh
      xprime(0)=x0(0)+dt
      xprime(3)=x0(3)+dz
      call givelocalhydro(xprime,temp,beta,densep,ifail,imed)
      xprime(0)=x0(0)-dt
      xprime(3)=x0(3)-dz
      call givelocalhydro(xprime,temp,beta,densem,ifail,imed)
      derdensedtb=(densep-densem)/(2*dtb)      
      derdense(0)=derdensedtb*ch
      derdense(3)=-derdensedtb*sh
      xprime(0)=tb*cosh(eta+deta)
      xprime(3)=tb*sinh(eta+deta)
      call givelocalhydro(xprime,temp,beta,densep,ifail,imed)
      xprime(0)=tb*cosh(eta-deta)
      xprime(3)=tb*sinh(eta-deta)
      call givelocalhydro(xprime,temp,beta,densem,ifail,imed)
      derdensedeta=(densep-densem)/(2*deta)
      derdense(0)=derdense(0)-x0(3)/tb2*derdensedeta
      derdense(3)=derdense(3)+x0(0)/tb2*derdensedeta
      normgrad=derdensedtb**2-(derdense(1)**2+derdense(2)**2)
      if(normgrad.gt.0) then
         tl=+1
      else
         tl=-1
      endif
      normgrad=sqrt(abs(normgrad))
      do k=0,3
         dsigunit(k)=-derdense(k)/normgrad
      enddo
      return
      end




      subroutine getetherm(mass,temp,etherm)
      implicit none
      integer nit,i
      parameter(nit=30)
      double precision mass,temp,etherm,alpha,beta,ran2,f,ftry,g1,g2,w,
     &  wtry,u
      alpha=mass/temp
      beta=alpha/(1.5+alpha)
      f=0
      w=sqrt(2*alpha) 
      do 10 i=1,nit
         u=-2.0*log(ran2())
         call gauss(g1,g2) 
         u=u+g2**2
         if(ran2().gt.beta) u=u-2.0*log(ran2())
         ftry=u/2
         wtry=sqrt(2*alpha+ftry)
         if(wtry.lt.w)then
            if(w*ran2().gt.wtry) goto 10
         endif
         f=ftry
         w=wtry
 10   continue 
      etherm=mass+f*temp
      return
      end

      subroutine getonepcooper(mass,temp,n4,p4)
      implicit none
      integer nit,nit2,i,j
      parameter(nit=50,nit2=1000)
      double precision mass,temp,alpha,beta,ran2,f,ftry,g1,g2,w,wtry,u,
     & ctheta,cthetatry,n4(0:3),nv,p,nvec(3),porth(3),p4(0:3),stheta
      integer ifmtx
      call getMonitorFileIndex(ifmtx)
      alpha=mass/temp
      beta=alpha/(1.5+alpha)
      nv=sqrt(n4(1)**2+n4(2)**2+n4(3)**2)
      f=0
      w=sqrt(2*alpha) 
      ctheta=0
      do 10 i=1,nit
         wtry=-1.
         j=0
         do while((wtry.lt.0).and.(j.lt.nit2))
            j=j+1
            if(j.eq.nit2) then
               write(ifmtx,*) 'out of loop in getonepcooper'
               write(ifmtx,*) 'alpha=',alpha
               close(ifmtx)
               stop
            endif
            u=-2.0*log(ran2())
            call gauss(g1,g2) 
            u=u+g2**2
            if(ran2().gt.beta) u=u-2.0*log(ran2())
            ftry=u/2
            cthetatry=2*ran2()-1
            wtry=n4(0)+nv*sqrt(ftry*(2*alpha+ftry))/(ftry+alpha)*
     &           cthetatry
         enddo
         wtry=wtry*sqrt(2*alpha+ftry)
         if(wtry.lt.w)then
            if(w*ran2().gt.wtry) goto 10
         endif
         f=ftry
         w=wtry
         ctheta=cthetatry
 10   continue 
      p4(0)=mass+f*temp
      p=sqrt(p4(0)**2-mass**2)
      do i=1,3
         nvec(i)=n4(i)/nv
      enddo
      call giveperpvec(nvec,porth)
      stheta=sqrt(1-ctheta**2)
      do i=1,3
         p4(i)=p*(ctheta*nvec(i)+stheta*porth(i))
      enddo
      if(abs(p4(1)**2+p4(2)**2+p4(3)**2-p**2).gt.1e-4) then
         write(ifmtx,*) 'checksum error in getonecooper'
         close(ifmtx)
         stop
      endif
      return
      end

      subroutine getonepcooperbis(mass,temp,n4,u4,p4)
      implicit none
      integer i
      double precision mass,temp,n4(0:3),p4(0:3),u4(0:3),n4influid(0:3),
     & p4influid(0:3),ubackinfluid(0:3),n4cov(0:3),n4covinfluid(0:3)
      n4cov(0)=n4(0)
      n4cov(1)=-n4(1)
      n4cov(2)=-n4(2)
      n4cov(3)=-n4(3)
      call lorentz2(u4,n4cov,n4covinfluid)
      n4influid(0)=n4covinfluid(0)
      n4influid(1)=-n4covinfluid(1)
      n4influid(2)=-n4covinfluid(2)
      n4influid(3)=-n4covinfluid(3)
      call getonepcooper(mass,temp,n4influid,p4influid)
      ubackinfluid(0)=u4(0)
      do i=1,3
         ubackinfluid(i)=-u4(i)
      enddo
      call lorentz2(ubackinfluid,p4influid,p4)
      return
      end


      subroutine get1qforcoal(u0,u1,plight,costh,lambda2)
      use PBGblcktransi
      use PBGfragandcoal
      implicit none
      double precision u0,u1,g1,g2,g3,ran2,plight2,plight,plighttry,
     &  elight,alph,alphbb,w,wtry,vmin,vmax,dv,y,ytry,costh,lambda2
      integer nit,i
      parameter(nit=30)
      plight=0
      w=exp(-u0*mqconstit/temptrans)
      do 10 i=1,nit
         g1=-2.0*log(ran2())
         call gauss(g2,g3) 
         plight2=Lambda2*(g1+g2**2)
         plighttry=sqrt(plight2)
         elight=sqrt(plight2+mqconstit2)
         alph=u0*elight/temptrans
         alphbb=u1*plighttry/temptrans
         wtry=((1+alph)*sinh(alphbb)/alphbb-cosh(alphbb))
     &        /alph*exp(-alph)
         if(wtry.lt.w)then
            if(w*ran2().gt.wtry) goto 10
         endif
         plight=plighttry
         w=wtry
 10   continue 
      vmin=exp(-2*alphbb)
      vmax=1
      dv=vmax-vmin
      y=0
      w=alph+alphbb
      do 20 i=1,nit
         ytry=-log(vmin+ran2()*dv)
         if(ytry.gt.2*alphbb) goto 20
         wtry=ytry+alph-alphbb
         if(wtry.lt.w)then
            if(w*ran2().gt.wtry) goto 20
         endif
         y=ytry
         w=wtry
 20   continue
      costh=y/alphbb-1      
      return
      end

      subroutine get1qforcoalnew(mquds,pql)
      use PBGgenvar
      use PBGblcktransi
      use PBGfragandcoal
      implicit none
      double precision mquds,pql(0:3),phiq,thetaq,pqtest,fvalue,
     & fdismax,ran2,pi
      parameter (pi=3.14159265358979323844d0)
      fvalue  = 0.0 
      fdismax = 0.1
      do while(ran2()*fdismax.gt.fvalue)
        phiq=2.0*pi*ran2()
        thetaq=pi*ran2()
        pqtest=6.0*ran2()
        pql(1)=pqtest*sin(thetaq)*cos(phiq)
        pql(2)=pqtest*sin(thetaq)*sin(phiq)
        pql(3)=pqtest*cos(thetaq)
        pql(0)=sqrt(mquds**2+pqtest**2)
        fvalue=6.0/(exp(pql(0)/tempcoal)+1.0)*
     & pqtest**2*sin(thetaq)
      enddo


      return
      end


      subroutine coalesce(pq,iresq,ph,iresh)
      use PBGfragandcoal
      use PBGgenvar
      implicit none
      integer iresq,iresh,i
      double precision pq(0:3),ph(0:3),u12,uq0,uqvec(3),uqvecnorm,
     & plight,
     & costh,pqnorm2,plonghrel,ptransh,phincframe(0:3),
     & pperp(3),
     & ubackinfluid(0:3),lambda2
      integer ifmtx
      call getMonitorFileIndex(ifmtx)
      pqnorm2=pq(1)**2+pq(2)**2+pq(3)**2
      u12=(pqnorm2)/resinfo_mass(iresq)**2
      uq0=sqrt(u12+1)
      do i=1,3
         uqvec(i)=pq(i)/resinfo_mass(iresq)
      enddo
      uqvecnorm=sqrt(u12)
      if(itypquark(iresq).eq.1) then
         lambda2=lambdaD2
      else if(itypquark(iresq).eq.2) then
         lambda2=lambdaB2
      else
         write(ifmtx,*) 'itypquark not recognized in coalesce'
         close(ifmtx)
         stop
      endif
      call get1qforcoal(uq0,uqvecnorm,plight,costh,lambda2)
      plonghrel=plight*costh/sqrt(pqnorm2)
      ptransh=plight*sqrt(1-costh**2)
      call giveperpvec(uqvec,pperp)
      phincframe(0)=resinfo_mass(iresh)**2
      do i=1,3
         phincframe(i)=plonghrel*pq(i)+ptransh*pperp(i)
         phincframe(0)=phincframe(0)+phincframe(i)**2
         ubackinfluid(i)=-uqvec(i)
      enddo
      phincframe(0)=sqrt(phincframe(0))
      ubackinfluid(0)=uq0
      call lorentz2(ubackinfluid,phincframe,ph)
      return
      end

      subroutine coalesceall(pq,ph,iresh)
      use PBGfragandcoal
      use PBGgenvar
      implicit none
      integer iresh,j,iquark,icount
      logical accept,ifmeson
      double precision pq(0:3),ph(0:3),plq1(0:3),plq2(0:3),massinv,
     & ufl(0:3),pcomf(0:3),plcomf(0:3),wigner,wignerall,wignerdq,
     & ran2,deltap(0:3),pl1comf(0:3),pl2comf(0:3),pdq(0:3),ufldq(0:3),
     & deltaqabs,sigma,sigmadq,massinvdq,pdqcomf(0:3),deltaqdq,mq1,mq2,
     & mheavy,mdq,deltapl(0:3)
      parameter(sigmadq=3.7245)
      parameter(sigma=3.7245)
      if(iresh.le.13)then
        iquark=1
      else                   
        iquark=14
      endif
      mheavy=resinfo_mass(iquark) 
      if(iresh.eq.4.or.iresh.eq.5.or.iresh.eq.17.or.
     & iresh.eq.18)then
           mq1=mqconstit
           ifmeson=.true.
      elseif(iresh.eq.6.or.iresh.eq.7.or.iresh.eq.19.or.
     & iresh.eq.20)then
           mq1=msconstit
           ifmeson=.true.
      elseif(iresh.eq.8.or.iresh.eq.9.or.iresh.eq.21.or.
     & iresh.eq.22)then
          mq1=mqconstit
          mq2=mqconstit
          ifmeson=.false.
      elseif(iresh.eq.10.or.iresh.eq.11.or.iresh.eq.23.or.
     & iresh.eq.24)then
          mq1=mqconstit
          mq2=msconstit
          ifmeson=.false.
      elseif(iresh.eq.12.or.iresh.eq.13.or.iresh.eq.25.or.
     & iresh.eq.26)then
          mq1=msconstit
          mq2=msconstit
          ifmeson=.false.
      else 
           write(8,*) 'wrong iresh in subroutine coalesceall'
           close(8)
           stop
      endif
      if(ifmeson)then
           accept=.false.
           icount=0
           do while(.not.accept) 
             icount=icount+1
             call get1qforcoalnew(mq1,plq1)
             pq(0)=sqrt(mheavy**2+pq(1)**2+pq(2)**2+pq(3)**2)
             do j=0,3
               ph(j)=pq(j)+plq1(j)
             enddo
             massinv=sqrt(ph(0)**2-ph(1)**2-ph(2)**2-ph(3)**2)
             do j=0,3
              ufl(j)=ph(j)/massinv
             enddo
             call lorentz2(ufl,pq,pcomf)
             call lorentz2(ufl,plq1,plcomf)
             do j=0,3
               deltap(j)=(mq1*pcomf(j)-mheavy*plcomf(j))/(mheavy+mq1)
             enddo
             deltaqabs=sqrt(deltap(1)**2+deltap(2)**2+deltap(3)**2)

             wigner=exp(-deltaqabs*deltaqabs*sigma*sigma)
             if(ran2().lt.wigner.or.icount.gt.1000000)then
               accept=.true.
             endif
           enddo
      else
           accept=.false.
           icount=0
           do while(.not.accept) 
             icount=icount+1
             call get1qforcoalnew(mq1,plq1)
             call get1qforcoalnew(mq2,plq2)
             do j=0,3
               pdq(j)=plq1(j)+plq2(j)
             enddo
             massinvdq=sqrt(pdq(0)**2-pdq(1)**2-pdq(2)**2-pdq(3)**2)
             do j=0,3
               ufldq(j)=pdq(j)/massinvdq
             enddo  
             call lorentz2(ufldq,plq1,pl1comf)
             call lorentz2(ufldq,plq2,pl2comf)
             do j=0,3
                deltapl(j)=(mq2*pl1comf(j)-mq1*pl2comf(j))/(mq1+mq2)
             enddo
             deltaqdq=sqrt(deltapl(1)**2+deltapl(2)**2+deltapl(3)**2)
             wignerdq=exp(-deltaqdq*deltaqdq*sigmadq*sigmadq)

             mdq=mq1+mq2
             pdq(0)=sqrt(mdq**2+pdq(1)**2+pdq(2)**2+pdq(3)**2)
             do j=0,3
                ph(j)=pq(j)+pdq(j)
             enddo
             massinv=sqrt(ph(0)**2-ph(1)**2-ph(2)**2-ph(3)**2)
             do j=0,3
              ufl(j)=ph(j)/massinv
             enddo
             call lorentz2(ufl,pq,pcomf)
             call lorentz2(ufl,pdq,pdqcomf)

             do j=0,3
               deltap(j)=(mdq*pcomf(j)-mheavy*pdqcomf(j))/(mheavy+mdq)
             enddo
             deltaqabs=sqrt(deltap(1)**2+deltap(2)**2+deltap(3)**2)

               wigner=exp(-deltaqabs*deltaqabs*sigma*sigma)
             wignerall=wignerdq*wigner
               if(ran2().lt.wignerall.or.icount.gt.1000000)then
               accept=.true.
             endif
           enddo 
      endif
      return
      end


      subroutine get1qforcoalTL(d,dsig0,u0,ifail)
      implicit none
      logical accept
      integer ifail
      double precision d,dsig0,dsigv,v0,u0
      integer nit,i,nit2,j
      parameter(nit=30,nit2=10000)
      double precision beta,x,x2,g1,g2,ran2,xtry,w,wtry,u0bis,uv,v
      integer ifmtx
      call getMonitorFileIndex(ifmtx)
      ifail=0
      if(dsig0.lt.-1d0) then
         ifail=-1
      endif
      beta=d/(1.5+d)
      x=0
      w=sqrt(2*d+x) 
      dsigv=sqrt(dsig0**2+1)
      v0=dsigv/dsig0
      do 10 i=1,nit
         accept=.false.
         j=0
         do while(.not.accept)
            j=j+1
            if(j.gt.nit2) then
              write(ifmtx,*) 'max nit2 reached in get1qforcoalTL:',dsig0
              close(ifmtx)
              stop
            endif
            x2=-2.0*log(ran2())
            call gauss(g1,g2) 
            x2=x2+g2**2
            if(ran2().gt.beta) x2=x2-2.0*log(ran2())
            xtry=x2/2
            u0bis=(1+xtry/d)
            uv=sqrt(u0bis**2-1)
            accept=(uv.gt.-dsig0)
         enddo
         wtry=sqrt(2*d+xtry)
         if(u0bis.gt.dsigv) then
            v=v0*uv/u0bis
            wtry=wtry*(1.d0+v)**2/(4.d0*abs(v))
         endif
         if((wtry.lt.w).and.(i.ne.1)) then
            if(w*ran2().gt.wtry) goto 10
         endif
         x=xtry
         w=wtry
 10   continue 
      u0=1+x/d
      return
      end

      subroutine coalescebis(pqinfl,dsigcovfl,iresq,tl,uQucell,uQdsig,
     &   ucelldsig,plightqinfl,mass2Qq,phinfl,iresh)
      use PBGfragandcoal
      use PBGgenvar
      use PBGsceplas
      implicit none
      integer iresq,iresh,i,tl,ifail,nbfail
      save nbfail
      data nbfail/0/
      double precision pqinfl(0:3),dsigcovfl(0:3),uQucell,uQdsig,
     & ucelldsig,u0intilde,uqintilde(0:3),uvecintilde,utilde(0:3),
     & dsigcovut(0:3),dsigcovutvec(3),dsigutvecnorm,costh,bigv,pperp(3),
     & ulong,utrans,ubackinfluid(0:3),uqinfl(0:3),phinfl(0:3),mu,d,
     & check,plightqinfl(0:3),mass2Qq
      double precision ran2
      integer ifmtx
      call getMonitorFileIndex(ifmtx)
      mu=tempcoal/(mqconstit*alphadov**2)
      d=mqconstit/tempcoal*sqrt(1+2*mu*uQucell+mu**2)
      do i=1,3
         utilde(i)=(pqinfl(i)/(resinfo_mass(iresq)*alphadov**2))/d
      enddo
      utilde(0)=(pqinfl(0)/(resinfo_mass(iresq)*alphadov**2)+
     &          mqconstit/tempcoal)/d
      check=abs(utilde(0)**2-utilde(1)**2-utilde(2)**2-utilde(3)**2-1)
      if(check.gt.1E-5)then
         write(ifmtx,*) 'self consistency problem in coalesce bis'
         write(ifmtx,*) (utilde(i),i=0,3)
         write(ifmtx,*) check
         close(ifmtx)
         stop
      endif
      call lorentz2(utilde,dsigcovfl,dsigcovut)
      dsigutvecnorm=sqrt(dsigcovut(1)**2+dsigcovut(2)**2+
     &              dsigcovut(3)**2)
      if(tl.gt.0) then 
         call getetherm(1.d0,1/d,u0intilde)
      else
         call get1qforcoalTL(d,dsigcovut(0),u0intilde,ifail)
         if(ifail.ne.0) then
            if(nbfail.le.20) then
               nbfail=nbfail+1
               write(8,*) 'dsigtilde0<0' 
               write(8,*) pqinfl,dsigcovfl,tl,uQucell,uQdsig,ucelldsig
            endif
         endif
      endif
      uvecintilde=sqrt(u0intilde**2-1)
      bigv=dsigutvecnorm*uvecintilde/(u0intilde*dsigcovut(0))
      if((tl.gt.0).or.((bigv.gt.0).and.(bigv.le.1))) then 
         costh=(1-sqrt((1+bigv)**2-4*bigv*ran2()))/bigv
      else
         costh=(1-(1+bigv)*sqrt(ran2()))/bigv
      endif
      if(abs(costh).ge.1) then
         write(ifmtx,*) '|cos theta| >1 in coalescebis:',tl,bigv,costh
         write(ifmtx,*) dsigutvecnorm,uvecintilde,u0intilde,dsigcovut(0)
         close(ifmtx)
         stop
      endif
      dsigcovutvec(1)=dsigcovut(1)/dsigutvecnorm
      dsigcovutvec(2)=dsigcovut(2)/dsigutvecnorm
      dsigcovutvec(3)=dsigcovut(3)/dsigutvecnorm
      call giveperpvec(dsigcovutvec,pperp)
      uqintilde(0)=u0intilde
      ulong=costh
      utrans=sqrt(1-costh**2)
      do i=1,3
         uqintilde(i)=uvecintilde*(ulong*dsigcovutvec(i)+utrans*
     &        pperp(i))
      enddo
      ubackinfluid(0)=utilde(0)
      do i=1,3
         ubackinfluid(i)=-utilde(i)
      enddo
      call lorentz2(ubackinfluid,uqintilde,uqinfl)
      check=abs(uqinfl(0)**2-uqinfl(1)**2-uqinfl(2)**2-uqinfl(3)**2-1)
      if(check.gt.1E-5)then
         write(ifmtx,*) 'self consistency prob in coalescebis:uqinfl'
         write(ifmtx,*) (uqinfl(i),i=0,3)
         write(ifmtx,*) check
         close(ifmtx)
         stop
      endif
      plightqinfl(0)=mqconstit*uqinfl(0)
      mass2Qq=(pqinfl(0)+plightqinfl(0))**2
      phinfl(0)=resinfo_mass(iresh)**2
      do i=1,3
         plightqinfl(i)=mqconstit*uqinfl(i)
         phinfl(i)=pqinfl(i)+plightqinfl(i)
         mass2Qq=mass2Qq-phinfl(i)**2
         phinfl(0)=phinfl(0)+phinfl(i)**2
      enddo
      if(mass2Qq.le.0) then
         write(ifmtx,*) 'inv mass squared<=0 for Q-q system'
         write(ifmtx,*) mass2Qq
         write(ifmtx,*) pqinfl
         write(ifmtx,*) plightqinfl
         close(ifmtx)
         stop
      endif
      phinfl(0)=sqrt(phinfl(0))
      return
      end




      subroutine gimmeprobcoal(itypquark,p,pb_coal)
      use PBGforprobcoal
      implicit none
      integer itypquark
      integer ipinf,ipsup
      double precision p,pb_coal,prap,dp
      if((p.lt.0.d0).or.(p.gt.pmax(itypquark))) then 
         pb_coal=0.
         return
      endif
      prap=p/deltp(itypquark)
      ipinf=int(prap)
      ipsup=ipinf+1
      dp=prap-ipinf
      pb_coal=dp*probcoal(itypquark,ipsup)+(1-dp)
     &       *probcoal(itypquark,ipinf)
      return
      end

      subroutine gimmeprobcoalHad(itypquark,p,pb_coalH)
      use PBGforprobcoal
      implicit none
      integer i,itypquark
      integer ipinf,ipsup
      double precision p,prap,dp,pb_coalH(5)
      
      if((p.lt.0.d0).or.(p.gt.pmax(itypquark))) then
         do i=1,5
            pb_coalH(i)=0.
         enddo
         return
      endif
      prap=p/deltp(itypquark)
      ipinf=int(prap)
      ipsup=ipinf+1
      dp=prap-ipinf
      do i=1,5
      pb_coalH(i)=dp*probcoalHad(itypquark,i,ipsup)+(1-dp)
     &       *probcoalHad(itypquark,i,ipinf)
      enddo
      return
      end


      subroutine gimmeprobcoalTL(itq,d,yds,pb_coal,ifail)
      use PBGforprobcoalTL
      implicit none
      integer ifail,itq
      double precision d,yds,deff,ydseff,pb_coal
      integer iuinf,iusup,ivinf,ivsup
      real aii,asi,ass,ais,u,du,dup,v,dv,dvp
      ifail=0
      if(d.gt.dcoalmax(itq)) then
         pb_coal=0.d0
         return
      else
         if(d.lt.dcoalmin(itq)) then
            write(8,*) 'D donwto ',d,' encountered in gimmeprocoalTL'
            ifail=1
            deff=dcoalmin(itq)
         else
            deff=d
         endif
      endif
      if(yds.gt.ydsmax(itq)) then
         write(8,*) 'y_dS upto ',yds,' encountered in gimmeprocoalTL'
         ydseff=ydsmax(itq)
      else
         if(yds.lt.ydsmin(itq)) then
            write(8,*) 'y_dS donwto ',yds,' encountin gimmeprocoalTL'
            ydseff=ydsmin(itq)
         else
            ydseff=yds
         endif
      endif
      if(deff.eq.dcoalmax(itq)) then
         if(ydseff.eq.ydsmax(itq))then
            pb_coal=probcoalTL(itq,nbd(itq),nbyds(itq))
         else
            v=(ydseff-ydsmin(itq))/dyds(itq)
            ivinf=int(v)
            ivsup=ivinf+1
            dv=v-ivinf
            dvp=1-dv
            pb_coal=dv*probcoalTL(itq,nbd(itq),ivsup)+ 
     &              dvp*probcoalTL(itq,nbd(itq),ivinf)
         endif
      else 
         u=(deff-dcoalmin(itq))/deltdcoal(itq)
         iuinf=int(u)
         iusup=iuinf+1
         du=u-iuinf
         dup=1-du
         if(ydseff.eq.ydsmax(itq)) then
            pb_coal=du*probcoalTL(itq,iusup,nbyds(itq))+ 
     &              dup*probcoalTL(itq,iuinf,nbyds(itq))
         else
            v=(ydseff-ydsmin(itq))/dyds(itq)
            ivinf=int(v)
            ivsup=ivinf+1
            dv=v-ivinf
            dvp=1-dv
            aii=dup*dvp
            asi=du*dvp
            ass=du*dv
            ais=dup*dv
            pb_coal=aii*probcoalTL(itq,iuinf,ivinf)
     &             +asi*probcoalTL(itq,iusup,ivinf)
     &             +ass*probcoalTL(itq,iusup,ivsup)
     &             +ais*probcoalTL(itq,iuinf,ivsup)  
         endif
      endif
      return
      end

      subroutine gimmeprobcoalbis(itypquark,pq,
     &  tl,pqnorminfl,uQucell,uQdsig,ucelldsig,d,pb_coal)
      use PBGsceplas
      use PBGfragandcoal
      implicit none
      integer itypquark,tl,k,ifail
      double precision pq(0:3),
     & pqnorminfl,pb_coal,
     & pb_coal0,mu,uQucell,uQdsig,ucelldsig,d,dtilde0,yds,ff,sq
      integer ifmtx
      call getMonitorFileIndex(ifmtx)
      mu=tempcoal/(mqconstit*alphadov**2)
      sq=sqrt(1+2*mu*uQucell+mu**2)
      d=mqconstit/tempcoal*sq
      if(uQdsig.le.0.d0) then
         pb_coal=0.d0
         return
      endif
      if(tl.gt.0) then
         call gimmeprobcoal(itypquark,pqnorminfl,pb_coal0)
         ff=(mu+ucelldsig/uQdsig)/(mu+1.d0/uQucell)
         pb_coal=ff*pb_coal0
      else 
         if(itypcoal.eq.2) then
            dtilde0=(ucelldsig+mu*uQdsig)/d
         else
            dtilde0=(ucelldsig+mu*uQdsig)/sq
         endif
         yds=log(dtilde0+sqrt(1+dtilde0**2))
         call gimmeprobcoalTL(itypquark,d,yds,pb_coal0,ifail)
         if(ifail.ne.0) then
            write(ifmtx,*) 'gimmeprobcoalbis:',ifail
            write(ifmtx,*) itypquark,mu,uQucell,(pq(k),k=0,3)
            close(ifmtx)
            stop
         endif
         pb_coal=pb_coal0/uQdsig
      endif
      return
      end




      subroutine initsavehadro
      use PBGfragandcoal
      implicit none
      integer i,j
      nbhadro=0
      do i=1,200
         do j=-1,1
            nbhadroofpt(i,j)=0.
            avprobcoalofpt(i,j)=0.
            varprobcoalofpt(i,j)=0. 
            nbmesonsofpthadrocoal(i,j)=0.
            nbmesonsofpthadrofrag(i,j)=0.
            Haahadro(i,j)=0.
         enddo
         Haahadrocoal(i)=0.
      enddo
      do i=1,50
         do j=-1,1
            distribproba(i,j)=0.
         enddo
      enddo
      end

      subroutine closesavehadro(totcol)
      use PBGfragandcoal
      use PBGgenvar
      use PBGfordisplay
      implicit none
      integer i,j,jt,nbhqsmallpt(-1:1),totcol,lentrue
      double precision av(200,-1:1),sig(200,-1:1),ratio1,ratio2
      character*60 filename
      integer ifmtx
      call getMonitorFileIndex(ifmtx)
      do i=1,200
         if(nbmesonsofpthadrocoal(i,0).ne.(nbmesonsofpthadrocoal(i,-1)+
     &      nbmesonsofpthadrocoal(i,1))) then
            write(ifmtx,*) 'checksum prob frag:',i,
     &        nbmesonsofpthadrocoal(i,-1),
     &           nbmesonsofpthadrocoal(i,1),nbmesonsofpthadrocoal(i,0)
            close(ifmtx)
            stop
         endif
         do j=0,1
            jt=-1+2*j
            if(nbhadroofpt(i,jt).gt.0) then
               av(i,jt)=avprobcoalofpt(i,jt)/nbhadroofpt(i,jt)
               sig(i,jt)=varprobcoalofpt(i,jt)/nbhadroofpt(i,jt)-
     &              av(i,jt)**2
               if(sig(i,jt).lt.0) then
                  write(8,*) 'variance < 0 in closesavehadro'
                  write(8,*) i,jt,nbhadroofpt(i,jt),avprobcoalofpt(i,
     &                 jt),varprobcoalofpt(i,jt),av(i,jt),sig(i,jt)
               else
                  sig(i,jt)=sqrt(sig(i,jt))
               endif
               Haahadro(i,jt)=(nbmesonsofpthadrocoal(i,jt)+
     &       nbmesonsofpthadrofrag(i,jt))/nbhadroofpt(i,jt)   
            else
               av(i,jt)=0.d0
               sig(i,jt)=0.d0
            endif
         enddo
         if(nbhadroofpt(i,0).gt.0) then
            Haahadro(i,0)=(nbmesonsofpthadrocoal(i,0)+
     &       nbmesonsofpthadrofrag(i,0))/nbhadroofpt(i,0)
            Haahadrocoal(i)=nbmesonsofpthadrocoal(i,0)/nbhadroofpt(i,0)
         endif
      enddo
      do j=0,1
         jt=-1+2*j
         nbhqsmallpt(jt)=0.
         do i=1,5 
            nbhqsmallpt(jt)=nbhqsmallpt(jt)+nbhadroofpt(i,jt)
         enddo
      enddo
      filename='probdistribution_'//
     &      ycentname(1:lentrue(ycentname,10))//'.res'
      open(unit=9,file=filename)
      do i=1,50
         if(nbhqsmallpt(1).gt.0) then
            ratio1=distribproba(i,1)/nbhqsmallpt(1)
         else
            ratio1=0.d0
         endif
         if(nbhqsmallpt(-1).gt.0) then
            ratio2=distribproba(i,-1)/nbhqsmallpt(-1)
         else
            ratio2=0.d0
         endif
         write(9,*)(i-0.5)*0.1,ratio1,ratio2
      enddo
      close(9)
      filename='probcoalonSpaceLike_'//
     &      ycentname(1:lentrue(ycentname,10))//'.res'
      open(unit=9,file=filename)
      do i=1,200
         write(9,*) (i-0.5)*0.2,nbhadroofpt(i,1)/totcol,av(i,1),
     & sig(i,1),(nbmesonsofpthadrocoal(i,1)+nbmesonsofpthadrofrag(i,1))/
     &        totcol,nbmesonsofpthadrocoal(i,1)/totcol,
     &        nbmesonsofpthadrofrag(i,1)/totcol,Haahadro(i,1)
      enddo
      close(9)
      filename='probcoalonTimeLike_'//
     &      ycentname(1:lentrue(ycentname,10))//'.res'
      open(unit=9,file=filename)
      do i=1,200
         write(9,*) (i-0.5)*0.2,nbhadroofpt(i,-1)/totcol,av(i,-1),
     &        sig(i,-1),(nbmesonsofpthadrocoal(i,-1)+
     &  nbmesonsofpthadrofrag(i,-1))/totcol,nbmesonsofpthadrocoal(i,-1)/
     & totcol,nbmesonsofpthadrofrag(i,-1)/totcol,Haahadro(i,-1)
      enddo
      close(9)
      filename='spectrahadro_'//
     &      ycentname(1:lentrue(ycentname,10))//'.res'
      open(unit=9,file=filename)
      do i=1,200
         write(9,*) (i-0.5)*0.2,nbhadroofpt(i,0)/totcol,
     &    (nbmesonsofpthadrocoal(i,0)+nbmesonsofpthadrofrag(i,0))/
     &    totcol,nbmesonsofpthadrocoal(i,0)/totcol,
     & nbmesonsofpthadrofrag(i,0)/totcol,Haahadro(i,0),Haahadrocoal(i)
      enddo
      close(9) 
      return
      end



      double precision function reduction_dof(iglfr,ppp,temp)
      use PBGgenvar
      use PBGsceplas
      use PBGreductiondofs
      implicit none
 
      integer iglfr
      double precision temp,ppp
      double precision lambdaboundstates
      external lambdaboundstates
      double precision lambdaepos,effdof
      external lambdaepos
      double precision massreduction
      external massreduction
      integer ifmtx
      call getMonitorFileIndex(ifmtx)
      if(isce<7) then 
           effdof=1.0d0
      endif
         reduction_dof=effdof

      if(isce==7) then
      if(reddof==0) then
         if(temp>temptrans) then
           effdof=1.0d0
         else
           effdof=0.0d0
         endif
         reduction_dof=effdof
      elseif(reddof==1) then
         effdof=lambdaboundstates(temp)
         reduction_dof=effdof
      elseif(reddof==2) then
         effdof=lambdaepos(temp)
         reduction_dof=effdof
      elseif(reddof==3) then
         if(temp>temptrans) then
           effdof=massreduction(iglfr,temp,ppp)
         else
           effdof=0.0d0
         endif
         reduction_dof=effdof
       else
         write(ifmtx,*) 'reddof=',reddof,'invalid!!!'
         close(ifmtx)
         stop
      endif
      endif

      return
      end 

      double precision function massreduction(iglfr,temp,ppp)
      use PBGforratedat3
      implicit none
      double precision:: temp,ppp,reddo,deltappp,deltt
      double precision:: pmaxx,pminn,tempmaxx,tempminn
      double precision:: prap,trap,dp,dtemp
      double precision:: aii,ais,asi,ass
      integer:: j,iglfr
      integer::ipsup,ipinf,itsup,itinf
      pmaxx=80.0d0
      tempmaxx=0.934d0
      pminn=0.5d0
      tempminn=0.144d0
      deltt=0.01
      if(ppp>pmaxx) ppp=pmaxx
      if(temp>tempmaxx) temp=tempmaxx
      if(ppp<pminn) ppp=pminn
      if(temp<tempminn) temp=tempminn
      j=int((temp-tempminn)/deltt)+1
      if(j<0) j=0
       if(ppp<5.0d0) then
        deltappp=0.5
        if(temp==tempmaxx) then
           prap=ppp/deltappp
           ipinf=int(prap)
           ipsup=ipinf+1
           dp=prap-ipinf
           if(iglfr==1) then
           reddo=dp*redthmass(ipsup,80)+(1-dp)*redthmass(ipinf,80)
           elseif(iglfr==0) then
           reddo=dp*redthmassgl(ipsup,80)+(1-dp)*redthmassgl(ipinf,80)
           endif
        else
          trap=(temp-tempminn)/deltt+1
          itinf=int(trap)
          itsup=itinf+1
          dtemp=trap-itinf
          prap=ppp/deltappp
          ipinf=int(prap)
          ipsup=ipinf+1
          dp=prap-ipinf
          aii=(1-dp)*(1-dtemp)
          ais=(1-dp)*dtemp
          asi=dp*(1-dtemp)
          ass=dp*dtemp
          if(iglfr==1) then
          reddo=aii*redthmass(ipinf,itinf)+ais*redthmass(ipinf,itsup)
     &         +asi*redthmass(ipsup,itinf)+ass*redthmass(ipsup,itsup)
           elseif(iglfr==0) then
          reddo=aii*redthmassgl(ipinf,itinf)
     &          +ais*redthmassgl(ipinf,itsup)
     &   +asi*redthmassgl(ipsup,itinf)+ass*redthmassgl(ipsup,itsup)
           endif
          endif
  
      elseif(ppp>=5.0d0.and.ppp<10.0d0)then
        deltappp=1.0
        if(temp==tempmaxx) then
           prap=(ppp-5.0d0)/deltappp
           ipinf=int(prap)+10
           ipsup=ipinf+1
           dp=prap-int(prap)
           if(iglfr==1) then
             reddo=dp*redthmass(ipsup,80)+(1-dp)*redthmass(ipinf,80)
           elseif(iglfr==0) then
             reddo=dp*redthmassgl(ipsup,80)+(1-dp)*redthmassgl(ipinf,80)
           endif
         else
          trap=(temp-tempminn)/deltt+1
          itinf=int(trap)
          itsup=itinf+1
          dtemp=trap-itinf
          prap=(ppp-5.0d0)/deltappp
          ipinf=int(prap)+10
          ipsup=ipinf+1
          dp=prap-int(prap)
          aii=(1-dp)*(1-dtemp)
          ais=(1-dp)*dtemp
          asi=dp*(1-dtemp)
          ass=dp*dtemp
           if(iglfr==1) then
          reddo=aii*redthmass(ipinf,itinf)+ais*redthmass(ipinf,itsup)
     &         +asi*redthmass(ipsup,itinf)+ass*redthmass(ipsup,itsup)
           elseif(iglfr==0) then
         reddo=aii*redthmassgl(ipinf,itinf)+ais*redthmassgl(ipinf,itsup)
     &        +asi*redthmassgl(ipsup,itinf)+ass*redthmassgl(ipsup,itsup)
            endif
         endif
      elseif(ppp>=10.0d0.and.ppp<15.0d0)then
        if(temp==tempmaxx) then
           ipinf=15
           ipsup=ipinf+1
           dp=(ppp-10.0d0)/5.0d0
           if(iglfr==1) then
           reddo=dp*redthmass(ipsup,80)+(1-dp)*redthmass(ipinf,80)
           elseif(iglfr==0) then
           reddo=dp*redthmassgl(ipsup,80)+(1-dp)*redthmassgl(ipinf,80)
           endif
        else
          trap=(temp-tempminn)/deltt+1
          itinf=int(trap)
          itsup=itinf+1
          dtemp=trap-itinf
          ipinf=15
          ipsup=ipinf+1
          dp=(ppp-10.0d0)/5.0d0
          aii=(1-dp)*(1-dtemp)
          ais=(1-dp)*dtemp
          asi=dp*(1-dtemp)
          ass=dp*dtemp
           if(iglfr==1) then
          reddo=aii*redthmass(ipinf,itinf)+ais*redthmass(ipinf,itsup)
     &         +asi*redthmass(ipsup,itinf)+ass*redthmass(ipsup,itsup)
           elseif(iglfr==0) then
         reddo=aii*redthmassgl(ipinf,itinf)+ais*redthmassgl(ipinf,itsup)
     &        +asi*redthmassgl(ipsup,itinf)+ass*redthmassgl(ipsup,itsup)
 
         endif
         endif
         else
        deltappp=10.0
         if(temp==tempmaxx) then
           if(ppp==pmaxx)then
           if(iglfr==1) then
              reddo=redthmass(22,80)
           elseif(iglfr==0) then
              reddo=redthmassgl(22,80)
           endif
           else
             prap=(ppp-15.0d0)/deltappp
             ipinf=int(prap)+16
             ipsup=ipinf+1
             dp=prap-int(prap)
           if(iglfr==1) then
              reddo=dp*redthmass(ipsup,80)+(1-dp)*redthmass(ipinf,80)
           elseif(iglfr==0) then
             reddo=dp*redthmassgl(ipsup,80)+(1-dp)*redthmassgl(ipinf,80)
           endif
           endif
         else
          trap=(temp-tempminn)/deltt+1
          itinf=int(trap)
          itsup=itinf+1
          dtemp=trap-itinf
          if(ppp==pmaxx)then
           if(iglfr==1) then
          reddo=dtemp*redthmass(22,itinf)
     &             +(1.0d0-dtemp)*redthmass(22,itsup)
           elseif(iglfr==0) then
             reddo=dtemp*redthmassgl(22,itinf)
     &               +(1.0d0-dtemp)*redthmassgl(22,itsup)
           endif
          else
            prap=(ppp-15.0d0)/deltappp
            ipinf=int(prap)+16
            ipsup=ipinf+1
            dp=prap-int(prap)
            aii=(1-dp)*(1-dtemp)
            ais=(1-dp)*dtemp
            asi=dp*(1-dtemp)
            ass=dp*dtemp
           if(iglfr==1) then
             reddo=aii*redthmass(ipinf,itinf)+ais*redthmass(ipinf,itsup)
     &         +asi*redthmass(ipsup,itinf)+ass*redthmass(ipsup,itsup)
           elseif(iglfr==0) then
         reddo=aii*redthmassgl(ipinf,itinf)+ais*redthmassgl(ipinf,itsup)
     &        +asi*redthmassgl(ipsup,itinf)+ass*redthmassgl(ipsup,itsup)
           endif
         endif
 
       endif
       endif

      massreduction=reddo


      end


      double precision function lambdaboundstates(ttemp)
      use PBGreductiondofs
      implicit none
      double precision ttemp
      double precision lamb
      double precision Tc
      Tc=0.155d0

      if (ttemp.le.Tc) then
         lambdaboundstates=0.0d0
      elseif (ttemp.ge.cTc*Tc) then
         lambdaboundstates=1.0d0
      else
         lamb=1.0d0/exp((cTc*Tc-ttemp)/(ttemp-Tc)) 
         lambdaboundstates=lamb
      endif
      return
      end 

      double precision function lambdaepos(ttemp)
      implicit none
    
      double precision ttemp
      double precision lamb,xx,zz
      double precision Tc
      Tc=0.134737d0

      if (ttemp<Tc) then
         lambdaepos=0.0d0
      else
         xx=(ttemp-Tc)/0.24d0
         zz=xx/(1.0d0+xx/0.77d0)
         lamb=exp(-zz-3.0d0*zz*zz)
         lambdaepos=1.0d0-lamb
      endif
      return
      end
