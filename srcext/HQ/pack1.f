    
      subroutine givelocalhydro(x,temp,betaplasma,densener,ifail,imed)

      use PBGsceplas
      use PBGpartinmedium
      implicit none
      integer ifail,imed
      double precision x(0:3),temp,betaplasma(3),densener,edensoft
      imed=-2
         if(isce.eq.3) then
            call givelocalhydro3(x,temp,betaplasma,densener,ifail)
         endif
         if((isce.eq.5).or.(isce.eq.6)) then
            call givelocalhydrogrid(x,temp,betaplasma,densener,ifail)
         endif
         if((isce.eq.7)) then
            call givelocalhydroepos(x,temp,betaplasma,densener,ifail)
         endif
         if((isce.eq.1).or.(isce.eq.2).or.(isce.eq.4)) then
            call givelocalhydroanal(x,temp,betaplasma,ifail)
            if(ifail.eq.-1) then
               densener=1D10
            else
               densener=edensoft(temp)         
            endif
         endif
         if(ifail.eq.-1) then
            imed=-1
         else
            if(ifail.eq.1) then
               imed=3
            else
               if(ifail.eq.0) then
                 if(temp.ge.temptrans) then
                     imed=0
                  else
                     imed=1
                  endif
               else
                  write(8,*) 'unknown ifail code in givelocalhydro'
                  close(8)
                  stop
               endif
            endif
         endif  
      if(isce>7) then 
         write(8,*) 'unforseen plasma scenario index ingivelochydro'
         close(8)
         stop
      endif
      end
      

      subroutine givelocalhydro3(x,temp,betaplasma,densener,ifail)
      use PBGsceplas
      use PBGpartinmedium
      implicit none
      integer ifail
      logical first
      double precision x(0:3),temp,betaplasma(3),densener,tau,rt2,a,b,
     &  onethird,edensoft
      data first/.true./
      save first,onethird
      if(first) then
         onethird=1/3.
         first=.false.
      endif
      temp=0.
      betaplasma(1)=0.
      betaplasma(2)=0.
      betaplasma(3)=0.
      densener=0.
      ifail=1
      tau=x(0)**2-x(3)**2
      if (tau.le.0.) return
      tau=sqrt(tau)
      rt2=x(1)**2+x(2)**2
      a=rt2/ri0**2
      if(a.gt.1) return
      b=tau0/tau
      if (b.gt.1) then
         ifail=-1
         write(8,*)'unexpected in givelocalhydro3:plasma not formed yet'
      else
         temp=tempinit*(b*(1-a)**alphaprof)**onethird
         if(temp.gt.temptrans) then
            densener=edensoft(temp)            
            betaplasma(3)=x(3)/x(0)
            ifail=0
         endif
      endif  
      return
      end


      subroutine givelocalhydroanal(x,temp,betaplasma,ifail)
      use PBGsceplas
      use PBGpartinmedium
      implicit none
      integer ifail
      logical first
      double precision x(0:3),temp,betaplasma(3),tau,volpl,
     & rt,rl,onethird,vtransmax,factvtrans
      data first/.true./
      save first,volpl,onethird
      if(first) then
         volpl=(ri0+vr0*tau0+pbar/2.*tau0**2)**2*(li0/2+vl0*tau0)
         onethird=1/3.
         first=.false.
      endif
      tau=x(0)**2-x(3)**2
      if (tau.lt.0.) then 
         write(8,*) 'error in givelocalhydro: outside light cone'
         ifail=-1
         return
      endif
      tau=sqrt(tau)
      ifail=0
      temp=0.
      betaplasma(1)=0.
      betaplasma(2)=0.
      betaplasma(3)=0.
      if(tau.lt.tau0) then
         write(8,*) 'error in givelocalhydro: plasma is not formed yet'
         ifail=-1
         return
      endif
      if(typesce(isce).eq.1) then
        write(8,*) 'error in givelocalhydro: not defined for typesce=1'
        close(8)
        stop
      endif
      if(typesce(isce).eq.2) then
        rt=ri0+vr0*tau+pbar/2.*tau**2
        if((x(1)*x(1)+x(2)*x(2)).gt.rt*rt) then
           ifail=1
        else
           rl=li0/2+vl0*tau
           temp=tempinit*(volpl/(rl*rt*rt))**onethird
           vtransmax=vr0+pbar*tau
           if(vtransmax.gt.1) then
              write(8,*) 'beurk velocity in givelocalhydro: ',vtransmax
              close(8)
              stop
           endif
           factvtrans=vtransmax/rt*tau/x(0)
           betaplasma(1)=factvtrans*x(1)
           betaplasma(2)=factvtrans*x(2)
           betaplasma(3)=x(3)/x(0)
        endif
      endif
      return
      end


      double precision function edensoft(temp)
      use PBGblcktransi
      implicit none
      double precision temp
      if(temp.ge.temptrans) then
         edensoft=1600*temp**4+80*temp**3
      else
         edensoft=0.5*(temp/0.17)**3.5
      endif
      return
      end


      logical function inboundaryplas(x)
      use PBGsceplas
      use PBGpartinmedium
      implicit none
      logical inboundaryplasanal
      integer ifail
      double precision x(0:3),temp,beta(3),densener
      if((isce.eq.1).or.(isce.eq.2).or.(isce.eq.4)) then
         inboundaryplas=inboundaryplasanal(x)
         return
      endif
      if(isce.eq.3) then
         inboundaryplas=.false.
         call givelocalhydro3(x,temp,beta,densener,ifail)
         inboundaryplas=ifail.eq.0
         return
      endif
      if((isce.eq.5).or.(isce.eq.6)) then
         inboundaryplas=.false.
         call givelocalhydrogrid(x,temp,beta,densener,ifail)
         if(ifail.eq.-1) then 
            write(8,*) 'In function inboundaryplas:'
            write(8,*) 'ifail=-1 found in givelocalhydrogrid'
            write(8,*) 'dunno what to do with this for isce=5 or 6'
            write(8,*) x(0),x(1),x(2),x(3)
            close(8)
            stop
         else
            if(ifail.eq.1) return
            inboundaryplas=(densener.gt.entrans_max)
         endif
         return
      endif 
      write(8,*) 'Did not recognize the plasma scenario in inboundary'
      close(8)
      stop
      end
      

      logical function inboundaryplasanal(x)
      use PBGsceplas
      use PBGpartinmedium
      implicit none
      double precision x(0:3),a,b,tb2,tb
      inboundaryplasanal=.false.
      if(typesce(isce).eq.1) then
        a=ri0+vr0*x(0)+pbar/2.*x(0)**2
        b=li0/2+vl0*x(0)
        if((x(1)**2+x(2)**2)/a**2+x(3)**2/b**2.le.1) then
           inboundaryplasanal=.true.
        endif
      endif
      if(typesce(isce).eq.2) then
        tb2=x(0)**2-x(3)**2
        if(tb2.lt.0.) return 
        tb=sqrt(tb2)
        a=ri0+vr0*tb+pbar/2.*tb2
        b=li0/2+vl0*x(0)
        if(((x(1)**2+x(2)**2).lt.a**2).and.(abs(x(3)).le.b)) then
           inboundaryplasanal=.true.
        endif
      endif
      return
      end     


      subroutine givewindowhot(xin,xfin,tempcrit,tresol,ifail,tineffect,
     & tfineffect)
      implicit none
      integer ifail,imed,k,statin,statfin,statint
      logical ifinf,ifsup,monot
      double precision xin(0:3),xfin(0:3),tempcrit,tresol,tineffect,
     & tfineffect,betaplasma(3),xint(0:3),tempin,tempfin,tempint,ener,
     & tfloatsup,tfloatinf,invdtall,dt1,dt2,curv,tmaxest,tempmaxest
      ifail=0
      call givelocalhydro(xin,tempin,betaplasma,ener,statin,imed)
      ifinf=(((statin.eq.0).and.(tempin.gt.tempcrit)).or.(statin.eq.-1))
      call givelocalhydro(xfin,tempfin,betaplasma,ener,statfin,imed)
      ifsup=(((statfin.eq.0).and.(tempfin.gt.tempcrit))
     &  .or.(statfin.eq.-1)) 
      invdtall=1./(xfin(0)-xin(0))
      tineffect=xin(0)
      tfineffect=xfin(0)
      monot=.true.
      if(ifinf) then
         if(ifsup) then
            return
         else
            tfloatsup=xin(0)
            goto 20
         endif   
      else
         if(ifsup) then
            tfloatinf=xfin(0)
            goto 30
         else
            monot=.false.
         endif   
      endif
      do k=0,3
         xint(k)=0.5*(xin(k)+xfin(k))
      enddo
      call givelocalhydro(xint,tempint,betaplasma,ener,statint
     &    ,imed)      
      if((statint.eq.-1).or.((statint.eq.0)
     &  .and.(tempint.gt.tempcrit))) then
         tfloatsup=xint(0)
         tfloatinf=xint(0)
         goto 20
      endif
      if(statint.eq.1)then
         ifail=4
         return
      endif
      if(statfin.eq.1) then
         if(tempint.lt.tempin) then
            ifail=3
            return
         else   
            goto 123
         endif
      endif   
      if(statin.eq.1) then
         if(tempint.lt.tempfin) then
            ifail=3
            return
         else
            goto123
         endif
      endif  
      curv=(tempfin+tempin-2*tempint)
      if(curv.ge.0) then
         ifail=2
         return
      endif   
      tempmaxest=tempint-(tempin-tempfin)**2/(8*curv)
      tmaxest=xint(0)-(xfin(0)-xin(0))*(tempfin-tempin)/(4*curv)
      if((tempmaxest.lt.tempcrit).or.(tmaxest.lt.xin(0)).or.
     &     (tmaxest.gt.xfin(0)))then
         ifail=1
         return
      endif
      write(8,*) 'looking for maximum'
      xint(0)=tmaxest
      dt1=(xfin(0)-xint(0))*invdtall
      dt2=(xint(0)-xin(0))*invdtall
      do k=1,3
         xint(k)=dt2*xfin(k)+dt1*xin(k)
      enddo
      call givelocalhydro(xint,tempint,betaplasma,ener,statint,imed)
      if((statint.eq.1).or.(tempint.le.tempcrit)) then
         write(8,*) 'do not find promissed max in givewindowhot'
         ifail=1
         return
      endif   
      tfloatsup=xint(0)
      tfloatinf=xint(0)
 20   continue
      do while(tfineffect-tfloatsup.gt.tresol)
         xint(0)=0.5*(tfloatsup+tfineffect)
         dt1=(xfin(0)-xint(0))*invdtall
         dt2=(xint(0)-xin(0))*invdtall
         do k=1,3
            xint(k)=dt2*xfin(k)+dt1*xin(k)
         enddo
         call givelocalhydro(xint,tempint,betaplasma,ener,statint,imed)
         if(tempint.lt.tempcrit) then
            tfineffect=xint(0)
         else
            tfloatsup=xint(0)
         endif
      enddo       
      tfineffect=tfloatsup
      if(monot) return
 30   continue
      do while((tfloatinf-tineffect.gt.tresol).or.
     &          (.not.(tfloatinf.lt.tfineffect)))
         xint(0)=0.5*(tfloatinf+tineffect)
         dt1=(xfin(0)-xint(0))*invdtall
         dt2=(xint(0)-xin(0))*invdtall
         do k=1,3
            xint(k)=dt2*xfin(k)+dt1*xin(k)
         enddo
         call givelocalhydro(xint,tempint,betaplasma,ener,statint,imed)
         if(tempint.lt.tempcrit) then
            tineffect=xint(0)
         else
            tfloatinf=xint(0)
         endif
      enddo
      tineffect=tfloatinf
      return
 123  continue
      write(8,*) 'unforeseen case in givewindowhot'
      write(8,*) statin,statint,statfin
      write(8,*) tempin,tempint,tempfin
      ifail=1
      return
      end

      subroutine givewindowener(xin,xfin,epscrit,tresol,ifail,tineffect,
     & tfineffect)

      implicit none
      integer ifail,k,statin,statfin,statint,statintl,statintr,imed
      logical ifinf,ifsup,monot
      double precision xin(0:3),xfin(0:3),epscrit,tresol,tineffect,
     & tfineffect,betaplasma(3),xint(0:3),epsin,epsfin,epsint,temp,
     & tfloatsup,tfloatinf,invdtall,dt1,dt2,curv,tmaxest,epsmaxest,
     & xfloatsup(0:3),xfloatinf(0:3),xintl(0:3),xintr(0:3),epsintl,
     & epsintr
      logical goleft,goright
      ifail=0
      call givelocalhydro(xin,temp,betaplasma,epsin,statin,imed)
      ifinf=(((statin.eq.0).and.(epsin.gt.epscrit)).or.(statin.eq.-1))
      call givelocalhydro(xfin,temp,betaplasma,epsfin,statfin,imed)
      ifsup=(((statfin.eq.0).and.(epsfin.gt.epscrit))
     & .or.(statfin.eq.-1)) 
      invdtall=1./(xfin(0)-xin(0))
      tineffect=xin(0)
      tfineffect=xfin(0)
      monot=.true.
      if(ifinf) then
         if(ifsup) then
            return
         else
            tfloatsup=xin(0)
            goto 20
         endif   
      else
         if(ifsup) then
            tfloatinf=xfin(0)
            goto 30
         else
            monot=.false.
         endif   
      endif
      do k=0,3
         xint(k)=0.5*(xin(k)+xfin(k))
      enddo
      call givelocalhydro(xint,temp,betaplasma,epsint,statint
     &  ,imed)      
      if((statint.eq.-1).or.((statint.eq.0)
     & .and.(epsint.gt.epscrit))) then
         tfloatsup=xint(0)
         tfloatinf=xint(0)
         goto 20
      endif
      if(statint.eq.1)then
         ifail=5
         return
      endif
      if(statfin.eq.1) then
         if(epsint.lt.epsin) then
            ifail=4
            return
         else   
            goto 12
         endif
      endif   
      if(statin.eq.1) then
         if(epsint.lt.epsfin) then
            ifail=3
            return
         else
            goto12
         endif
      endif  
      goto 13
 12   continue
      do k=0,3
         xfloatsup(k)=xfin(k)
         xfloatinf(k)=xin(k)
      enddo
      do while((epsint.lt.epscrit).and.(.not.(statint.eq.-1)).and.
     &         ((xfloatsup(0)-xfloatinf(0)).gt.tresol))
         do k=0,3
            xintl(k)=0.5*(xfloatinf(k)+xint(k))
            xintr(k)=0.5*(xint(k)+xfloatsup(k))
         enddo
         call givelocalhydro(xintl,temp,betaplasma,epsintl,statintl
     &     ,imed)      
         call givelocalhydro(xintr,temp,betaplasma,epsintr,statintr
     &     ,imed)    
         goleft=((statintl.eq.-1).or.
     &           ((statintl.eq.0).and.(epsintl.gt.epsint))) 
         goright=((statintr.eq.-1).or.
     &           ((statintr.eq.0).and.(epsintr.gt.epsint))) 
         if(goleft)then
            if(goright)then            
             write(8,*) 'unforeseen case found in givewindowener: 2 max'
               ifail=6
               return
            else
               do k=0,3
                  xfloatsup(k)=xint(k)
                  xint(k)=xintl(k)
               enddo
               statint=statintl
               epsint=epsintl
            endif
         else
            if(goright)then            
               do k=0,3
                  xfloatinf(k)=xint(k)
                  xint(k)=xintr(k)
               enddo
               statint=statintr
               epsint=epsintr
            else
               do k=0,3
                  xfloatinf(k)=xintl(k)
                  xfloatsup(k)=xintr(k)
               enddo               
            endif   
         endif
      enddo   
      if((statint.eq.-1).or.(epsint.gt.epscrit))then
         tfloatsup=xint(0)
         tfloatinf=xint(0)
         write(8,*) '... some hot phase found indeed !'
         goto 20
      else
         ifail=1
         return
      endif
 13   continue
      curv=(epsfin+epsin-2*epsint)
      if(curv.ge.0) then
         ifail=2
         return
       endif   
      epsmaxest=epsint-(epsin-epsfin)**2/(8*curv)
      tmaxest=xint(0)-(xfin(0)-xin(0))*(epsfin-epsin)/(4*curv)
      if((epsmaxest.lt.epscrit).or.(tmaxest.lt.xin(0)).or.
     &     (tmaxest.gt.xfin(0)))then
         ifail=1
         return
      endif
      write(8,*) 'looking for maximum'
      xint(0)=tmaxest
      dt1=(xfin(0)-xint(0))*invdtall
      dt2=(xint(0)-xin(0))*invdtall
      do k=1,3
         xint(k)=dt2*xfin(k)+dt1*xin(k)
      enddo
      call givelocalhydro(xint,temp,betaplasma,epsint,statint,imed)
      if((statint.eq.1).or.(epsint.le.epscrit)) then
         write(8,*) 'do not find promissed max in givewindowener'
         ifail=1
         return
      endif   
      tfloatsup=xint(0)
      tfloatinf=xint(0)
 20   continue
      do while(tfineffect-tfloatsup.gt.tresol)
         xint(0)=0.5*(tfloatsup+tfineffect)
         dt1=(xfin(0)-xint(0))*invdtall
         dt2=(xint(0)-xin(0))*invdtall
         do k=1,3
            xint(k)=dt2*xfin(k)+dt1*xin(k)
         enddo
         call givelocalhydro(xint,temp,betaplasma,epsint,statint,imed)
         if(epsint.lt.epscrit) then
            tfineffect=xint(0)
         else
            tfloatsup=xint(0)
         endif
      enddo       
      tfineffect=tfloatsup
      if(monot) return
 30   continue
      do while((tfloatinf-tineffect.gt.tresol).or.
     &          (.not.(tfloatinf.lt.tfineffect)))
         xint(0)=0.5*(tfloatinf+tineffect)
         dt1=(xfin(0)-xint(0))*invdtall
         dt2=(xint(0)-xin(0))*invdtall
         do k=1,3
            xint(k)=dt2*xfin(k)+dt1*xin(k)
         enddo
         call givelocalhydro(xint,temp,betaplasma,epsint,statint,imed)
         if(epsint.lt.epscrit) then
            tineffect=xint(0)
         else
            tfloatinf=xint(0)
         endif
      enddo
      tineffect=tfloatinf
       write(8,*) 'unforeseen case found in givewindowener'
       write(8,*) xin(0),xin(1),xin(2),xin(3)
       write(8,*) xfin(0),xfin(1),xfin(2),xfin(3)
       write(8,*) statin,statint,statfin
       write(8,*) epscrit,epsin,epsint,epsfin
       ifail=1
      end


      subroutine givewindow(xin,xfin,paramcrit,whatparam,tresol,ifail,
     & tineffect,tfineffect)
      implicit none
      integer ifail,whatparam
      double precision xin(0:3),xfin(0:3),paramcrit,tresol,tineffect,
     &  tfineffect
      if(whatparam.eq.1) then
         call givewindowhot(xin,xfin,paramcrit,tresol,ifail,
     &        tineffect,tfineffect)
      else
         if(whatparam.eq.2) then
            call givewindowener(xin,xfin,paramcrit,tresol,ifail,
     &           tineffect,tfineffect)
         else
            write(8,*) 'givewindow : did not recognize', whatparam,
     &           'as a valid parameter indic'
            close(8)
            stop
         endif
      endif   
      return
      end


      double precision function lifetime_old(param,whatparam)
      use PBGsceplas
      use PBGpartinmedium
      implicit none
      integer k,ifail,whatparam
      double precision param,xin(0:3),xfin(0:3),tineffect,tfineffect
     &   ,tresol
      do k=1,3
         xin(k)=0.
         xfin(k)=0.
      enddo
      xin(0)=1.01*tau0
      xfin(0)=50.
      tresol=0.1 
      call givewindow(xin,xfin,param,whatparam,tresol,ifail,
     &     tineffect,tfineffect)
      if((abs(tineffect-xin(0)).gt.1d-7).or.(ifail.ne.0)) then
         write(8,*) 'problem in lifetime'
         write(8,*) abs(tineffect-xin(0)), ifail
         close(8)
         stop
      else
      endif   
      lifetime_old=tfineffect
      return
      end

      double precision function tempmaxf(tb)
      use PBGhydro
      implicit none
      integer ix,iy,imed,ifail
      double precision tb,x(0:3),temp,beta(3),densener,tempmax
      x(0)=tb
      x(3)=0
      tempmax=0
      do ix=0,nbx
         x(1)=ix*xstep
         do iy=0,nby
            x(2)=iy*ystep
            call givelocalhydro(x,temp,beta,densener,ifail,imed)
            if((ifail.eq.0).and.(temp.gt.tempmax)) then
               tempmax=temp
            endif
         enddo
      enddo
      tempmaxf=tempmax
      end

      double precision function densenermaxf(tb)
      use PBGhydro
      implicit none
      integer ix,iy,imed,ifail
      double precision tb,x(0:3),temp,beta(3),densener,densenermax

      x(0)=tb
      x(3)=0
      densenermax=0
      do ix=0,nbx
         x(1)=ix*xstep
         do iy=0,nby
            x(2)=iy*ystep
            call givelocalhydro(x,temp,beta,densener,ifail,imed)
            if((ifail.eq.0).and.(densener.gt.densenermax)) then
               densenermax=densener
            endif
         enddo
      enddo
      densenermaxf=densenermax
      end

      double precision function lifetime(param,whatparam)
      use PBGhydro
      implicit none
      logical trynext
      integer whatparam
      double precision param,tresol,tbase,timetry,densenermaxf,tempmaxf,
     &  paramnext,lifetime_old
      tbase=lifetime_old(param,whatparam)
      tresol=0.1 
      timetry=tbase+tresol
      trynext=.true.      
      do while(trynext.and.(timetry.lt.tbmax))
         if(whatparam.eq.1) then
            paramnext=tempmaxf(timetry)             
         else
            paramnext=densenermaxf(timetry)
         endif
         if(paramnext.ge.param) then
            tbase=timetry
         else
            trynext=.false.
         endif   
         timetry=timetry+tresol
      enddo
      if(trynext) then
         write(8,*) 'lifetime: ran off available times without success'
         close(8)
         stop
      endif
      lifetime=tbase
      return
      end

      double precision function tempmaxfepos(ntb)
      use PBGhydro
      implicit none
      integer ix,iy,iz,imed,ifail
      integer ntb
      double precision x(0:3),temp,beta(3),densener,tempmax
      double precision etaa,tauu

      tauu=tbmin+dble(ntb-1)*tbstep
      tempmax=0.d0
      do iz=1,nzhyx
        etaa=(iz-(nzhyx+1)/2)*zstep
        x(3)=tauu*dsinh(etaa)
        x(0)=tauu*dcosh(etaa)
      do ix=1,nxhyx
         x(1)=(ix-(nxhyx+1)/2)*xstep
         do iy=0,nyhyx
            x(2)=(iy-(nyhyx+1)/2)*ystep
            call givelocalhydro(x,temp,beta,densener,ifail,imed)
            if((ifail.eq.0).and.(temp.gt.tempmax)) then
                tempmax=temp
            endif
         enddo
      enddo
      enddo
      tempmaxfepos=tempmax
      end 


      subroutine checktemp
      use PBGhydro
      implicit none
      integer imed,ifail,ntb
      double precision x(0:3),temp,beta(3),densener
      double precision etaa,tauu
      open(unit=2220,file='temp_cent_of_time.dat')
      write(2220,*) ntauhyx
      do ntb =1,ntauhyx
         tauu=tbmin+dble(ntb-1)*tbstep
         etaa=0.d0
         x(3)=tauu*dsinh(etaa)
         x(0)=tauu*dcosh(etaa)
         x(1)=0.d0
         x(2)=0.d0
         call givelocalhydro(x,temp,beta,densener,ifail,imed)
         write(2220,*) x(0),x(3),temp,densener,ifail,imed
      enddo
      close(2220)
      end


      double precision function densenermaxfepos(ntb)
      use PBGhydro
      implicit none
      integer ix,iy,iz,imed,ifail
      integer ntb
      double precision x(0:3),temp,beta(3),densener,densenermax
      double precision etaa,tauu
 
      tauu=tbmin+dble(ntb-1)*tbstep
      densenermax=0
      do iz=1,nzhyx
      etaa=(iz-14)*zstep
      x(3)=tauu*dsinh(etaa)
      x(0)=tauu*dcosh(etaa)
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
      enddo
      densenermaxfepos=densenermax
      end


      subroutine givelifetimeepos(param,whatparam,lifetimeepos,ifail)
      use PBGhydro
      implicit none
      logical trynext
      integer ifail,whatparam
      integer ntimetry
      double precision param,
     &  densenermaxfepos
     &  ,tempmaxfepos, paramnext,lifetimeepos

      
      if(ntauhyx.eq.0) then
         lifetimeepos=0.d0
         ifail=1
         return
      endif
      if(whatparam.eq.1) then
      elseif(whatparam.eq.2) then
      endif
      ifail=0
      ntimetry=1
      trynext=.true.      
      do while(trynext.and.(ntimetry.le.ntauhyx))
         if(whatparam.eq.1) then
            paramnext=tempmaxfepos(ntimetry)             
         else
            paramnext=densenermaxfepos(ntimetry)
         endif
         if((ntimetry.eq.1).and.(paramnext.lt.param)) then
            write(8,*) 'lifetime: not even a plasma at 1st timestep'
            write(8,*) paramnext
            ifail=1
         endif
         if(paramnext.lt.param) then
            trynext=.false.
         endif   
         ntimetry=ntimetry+1
      enddo
      if(trynext) then
         write(8,*) 'lifetime: ran off times without finding end of QGP'
         write(8,*) 'scan stopped at',ntimetry
         close(8)
         stop
      endif
      lifetimeepos=tbmin+tbstep*dble(ntimetry-1)
      return
      end




      subroutine givelocalhydroepos(x,temp,beta,ener,ifail)


      implicit none
      integer ifail
      double precision x(0:3),tb,temp,beta(3),ener,etaa
      double precision velotrans
      real ftemp,fener,fvx,fvy,fvz
       
      tb=x(0)**2-x(3)**2
      if (tb.lt.0.) then 
         write(8,*) 'error in givelocalhydrogrid: outside light cone'
         write(8,*) x(0),x(1),x(2),x(3)
         ifail=-1
         return
      endif
      tb=sqrt(tb)
      etaa=0.5d0*log((x(0)+x(3))/(x(0)-x(3)))
      call azhydrogridepos(tb,x(1),x(2),etaa,ftemp,fvx,fvy,fvz,
     &                     fener,ifail)
      if(ifail.eq.-1) return
      temp=dble(ftemp)
      ener=dble(fener)
      if(temp<0.0001d0)then
         beta(1)=0d0
         beta(2)=0d0
         beta(3)=0d0
      else
         if((fvx.gt.1).or.(fvy.gt.1).or.(fvz.gt.1)) then 
            write(6,*) 'beurk velocity in givelocalhydrokolb (1)'
            write(*,*) fvx,fvy,fvz
            stop
         endif
         velotrans=tanh(etaa)
         beta(1)=dble(fvx)/cosh(etaa)/(1.0d0+velotrans*dble(fvz))
         beta(2)=dble(fvy)/cosh(etaa)/(1.0d0+velotrans*dble(fvz))
         beta(3)=(dble(fvz)+velotrans)/(1.0d0+velotrans*dble(fvz))
         if((beta(1).gt.1).or.(beta(2).gt.1).or.(beta(3).gt.1)) then 
            write(6,*) 'beurk velocity in givelocalhydrokolb (1)'
            write(*,*) beta(1),beta(2),beta(3)
            write(*,*) fvx,fvy,fvz
            write(*,*) 'velotrans = ',velotrans,", eta = ",etaa
            stop
         endif
       endif
      return
      end

      subroutine givelocalhydrogrid(x,temp,beta,ener,ifail)
      use PBGimparam
      implicit none
      integer ifail
      double precision x(0:3),temp,beta(3),ener,invgammalong
      real ftemp,vel,fener,fvx,fvy,tb,rt,bo
      tb=x(0)**2-x(3)**2
      if (tb.lt.0.) then 
         write(8,*) 'error in givelocalhydrogrid: outside light cone'
         write(8,*) x(0),x(1),x(2),x(3)
         ifail=-1
         return
      endif
      tb=sqrt(tb)
      if((b.eq.0.).and.ifaxial) then
         rt=sqrt(x(1)**2+x(2)**2)
         call hydrogrid(tb,rt,ftemp,vel,fener,ifail)
         if(ifail.eq.-1) return
         temp=ftemp
         ener=fener
         if((ifail.eq.1).or.(rt.eq.0)) then
            beta(1)=0.
            beta(2)=0.
            beta(3)=0.
         else   
            if(vel.gt.1) then
               write(8,*) 'beurk velocity in givelocalhydrogrid (1)'
               close(8)
               stop
            endif
            beta(3)=x(3)/x(0)
            invgammalong=dsqrt(1-beta(3)**2)
            bo=vel/rt*invgammalong
            beta(1)=bo*x(1) 
            beta(2)=bo*x(2) 
         endif   
      else
         call azhydrogrid(tb,x(1),x(2),ftemp,fvx,fvy,fener,ifail)
         if(ifail.eq.-1) return
         temp=ftemp
         ener=fener
         if((ifail.eq.1).or.((x(1).eq.0d0).and.(x(2).eq.0d0))) then
            beta(1)=0d0
            beta(2)=0d0
            beta(3)=0d0
         else
            if((fvx.gt.1).or.(fvy.gt.1)) then 
               write(6,*) 'beurk velocity in givelocalhydrokolb (1)'
               stop
            endif
            beta(3)=x(3)/x(0)
            invgammalong=dsqrt(1-beta(3)**2)
            if(x(1).lt.0) then
               beta(1)=-fvx*invgammalong
            else
               beta(1)=fvx*invgammalong
            endif
            if(x(2).lt.0) then
               beta(2)=-fvy*invgammalong
            else
               beta(2)=fvy*invgammalong
            endif
         endif
      endif
      return
      end


      subroutine hydrogrid(tb,rt,temp,vel,ener,ifail)
      use PBGhydro
      implicit none
      integer ifail,iuinf,iusup,ivinf,ivsup
      real tb,rt,temp,vel,ener,tempt,aii,asi,ass,ais
      real u,du,dup,v,dv,dvp
      ifail=-1
      if(tb.lt.tbmin) return
      temp=0.
      ener=0.
      vel=0.
      ifail=1
      if((tb.gt.tbmax).or.(rt.gt.rmax)) return
      if(tb.eq.tbmax) then
         if(rt.eq.rmax)then
            tempt=temphyd(nbtb,nbr)
            if(abs(tempt).lt.0.15) return
            vel=velhyd(nbtb,nbr)
            ener=enerhyd(nbtb,nbr)
         else
            v=rt/rstep
            ivinf=int(v)
            ivsup=ivinf+1
            dv=v-ivinf
            dvp=1-dv
            tempt=dv*temphyd(nbtb,ivsup)+ dvp*temphyd(nbtb,ivinf)
            if(abs(tempt).lt.0.15) return
            vel =dv*velhyd(nbtb,ivsup) + dvp*velhyd(nbtb,ivinf)
            ener=dv*enerhyd(nbtb,ivsup)
     &           +dvp*enerhyd(nbtb,ivinf)         
         endif
      else 
         u=(tb-tbmin)/tbstep
         iuinf=int(u)
         iusup=iuinf+1
         du=u-iuinf
         dup=1-du
         if(rt.eq.rmax) then
            tempt=du*temphyd(iusup,nbr)+ dup*temphyd(iuinf,nbr)
            if(abs(tempt).lt.0.15) return
            vel =du*velhyd(iusup,nbr) + dup*velhyd(iuinf,nbr)
            ener=du*enerhyd(iusup,nbr)+ dup*enerhyd(iuinf,nbr)
         else
            v=rt/rstep
            ivinf=int(v)
            ivsup=ivinf+1
            dv=v-ivinf
            dvp=1-dv
            aii=dup*dvp
            asi=du*dvp
            ass=du*dv
            ais=dup*dv
            tempt=aii*temphyd(iuinf,ivinf)+asi*temphyd(iusup,ivinf)
     &           +ass*temphyd(iusup,ivsup)+ais*temphyd(iuinf,ivsup)  
            if(abs(tempt).lt.0.15) return
            vel=aii*velhyd(iuinf,ivinf)+asi*velhyd(iusup,ivinf)
     &           +ass*velhyd(iusup,ivsup)+ais*velhyd(iuinf,ivsup)  
            ener=aii*enerhyd(iuinf,ivinf)+asi*enerhyd(iusup,ivinf)
     &           +ass*enerhyd(iusup,ivsup)+ais*enerhyd(iuinf,ivsup)  
         endif
      endif
      temp=tempt
      ifail=0
      return
      end

      subroutine azhydrogrideposold(tb,xx,yy,zz,temp,velx,vely,velz
     &                          ,ener,ifail)

      use PBGhydro
      implicit none
      integer ifail,iuinf,iusup,ivinf,ivsup,iwinf,iwsup,istinf,istsup
      real temp,velx,vely,velz,ener,enert
      real u,du,dup,v,dv,dvp,w,dw,dwp,st,dst,dstp
      double precision tb,xx,yy,zz
      real aiiii,aiisi,aisii,asiii,asisi,assii,aissi,asssi
      real aiiis,aiiss,aisis,asiis,asiss,assis,aisss,assss

      ifail=-1
      if(tb.lt.tbmin) return
      temp=0.0d0
      ener=0.0d0
      velx=0.0d0
      vely=0.0d0
      velz=0.0d0
      ifail=1
      if((tb.gt.tbmax).or.(dabs(xx).gt.xmax).or.(dabs(yy).gt.ymax)
     &     .or.(dabs(zz).gt.zmax))then
         return
      endif
******weights for all cases******
      u=(xx-xmin)/xstep
      iuinf=int(u)
      iusup=iuinf+1
      du=u-iuinf
      dup=1-du
      iuinf=iuinf+1
      iusup=iusup+1
      if(iuinf==nxhyx) iusup=nxhyx
      v=(yy-ymin)/ystep
      ivinf=int(v)
      ivsup=ivinf+1
      dv=v-ivinf
      dvp=1-dv
      ivinf=ivinf+1
      ivsup=ivsup+1
       if(ivinf==nyhyx) ivsup=nyhyx
      w=(tb-tbmin)/tbstep
      iwinf=int(w)
      iwsup=iwinf+1
      dw=w-iwinf
      dwp=1-dw
      iwinf=iwinf+1
      iwsup=iwsup+1
       if(iwinf==ntauhyx) iwsup=ntauhyx
      st=(zz-zmin)/zstep
      istinf=int(st)
      istsup=istinf+1
      dst=st-istinf
      dstp=1-dst
      istinf=istinf+1
      istsup=istsup+1
       if(istinf==nzhyx) istsup=nzhyx
      aiiii=dwp*dup*dvp*dstp
      aiiis=dwp*dup*dvp*dst
      aiisi=dwp*dup*dv*dstp
      aiiss=dwp*dup*dv*dst
      aisii=dwp*du*dvp*dstp
      aisis=dwp*du*dvp*dst
      asiii=dw*dup*dvp*dstp
      asiis=dw*dup*dvp*dst
      asisi=dw*dup*dv*dstp
      asiss=dw*dup*dv*dst
      aissi=dwp*du*dv*dstp
      aisss=dwp*du*dv*dst
      assii=dw*du*dvp*dstp
      assis=dw*du*dvp*dst
      asssi=dw*du*dv*dstp
      assss=dw*du*dv*dst

         if(velx>1.0d0)then
           write(*,*) tb,xx,yy,zz, 'velx= ', velx
         else if(vely>1.0d0)then
           write(*,*) tb,xx,yy,zz, 'vely= ', vely
         else if(velz>1.0d0)then
           write(*,*) tb,xx,yy,zz, 'velz= ', velz
          end if
      ener=enert
      ifail=0
      return
      end 

      subroutine azhydrogridepos(tb,xx,yy,zz,temp,velx,vely,velz
     &                          ,ener,ifail)
      use PBGhydro
      implicit none
      integer ifail,iuinf,iusup,ivinf,ivsup,iwinf,iwsup,istinf,istsup
      real temp,velx,vely,velz,ener,enert
      real u,du,dup,v,dv,dvp,w,dw,dwp,st,dst,dstp
      double precision tb,xx,yy,zz
      real aiiii,aiisi,aisii,asiii,asisi,assii,aissi,asssi
      real aiiis,aiiss,aisis,asiis,asiss,assis,aisss,assss
      real pss,epepsiiii,eptemiiii,v1iiii,v2iiii,v3iiii,
     &  epepsiisi,eptemiisi,v1iisi,v2iisi,v3iisi,
     &  epepsisii,eptemisii,v1isii,v2isii,v3isii,
     &  epepssiii,eptemsiii,v1siii,v2siii,v3siii,
     &  epepssisi,eptemsisi,v1sisi,v2sisi,v3sisi,
     &  epepsissi,eptemissi,v1issi,v2issi,v3issi,
     &  epepsssii,eptemssii,v1ssii,v2ssii,v3ssii,
     &  epepssssi,eptemsssi,v1sssi,v2sssi,v3sssi,
     &  epepsiiis,eptemiiis,v1iiis,v2iiis,v3iiis,
     &  epepsiiss,eptemiiss,v1iiss,v2iiss,v3iiss,
     &  epepsisis,eptemisis,v1isis,v2isis,v3isis,
     &  epepssiis,eptemsiis,v1siis,v2siis,v3siis,
     &  epepssiss,eptemsiss,v1siss,v2siss,v3siss,
     &  epepsisss,eptemisss,v1isss,v2isss,v3isss,
     &  epepsssis,eptemssis,v1ssis,v2ssis,v3ssis,
     &  epepsssss,eptemssss,v1ssss,v2ssss,v3ssss

      ifail=-1
      if(tb.lt.tbmin) return
      temp=0.0d0
      ener=0.0d0
      velx=0.0d0
      vely=0.0d0
      velz=0.0d0
      ifail=1
      if((tb.gt.tbmax).or.(dabs(xx).gt.xmax).or.(dabs(yy).gt.ymax)
     &     .or.(dabs(zz).gt.zmax))then
         return
      endif
******weights for all cases******
      u=(xx-xmin)/xstep
      iuinf=int(u)
      iusup=iuinf+1
      du=u-iuinf
      dup=1-du
      iuinf=iuinf+1
      iusup=iusup+1
      if(iuinf==nxhyx) iusup=nxhyx
      v=(yy-ymin)/ystep
      ivinf=int(v)
      ivsup=ivinf+1
      dv=v-ivinf
      dvp=1-dv
      ivinf=ivinf+1
      ivsup=ivsup+1
       if(ivinf==nyhyx) ivsup=nyhyx
      w=(tb-tbmin)/tbstep
      iwinf=int(w)
      iwsup=iwinf+1
      dw=w-iwinf
      dwp=1-dw
      iwinf=iwinf+1
      iwsup=iwsup+1
       if(iwinf==ntauhyx) iwsup=ntauhyx
      st=(zz-zmin)/zstep
      istinf=int(st)
      istsup=istinf+1
      dst=st-istinf
      dstp=1-dst
      istinf=istinf+1
      istsup=istsup+1
       if(istinf==nzhyx) istsup=nzhyx
      aiiii=dwp*dup*dvp*dstp
      aiiis=dwp*dup*dvp*dst
      aiisi=dwp*dup*dv*dstp
      aiiss=dwp*dup*dv*dst
      aisii=dwp*du*dvp*dstp
      aisis=dwp*du*dvp*dst
      asiii=dw*dup*dvp*dstp
      asiis=dw*dup*dvp*dst
      asisi=dw*dup*dv*dstp
      asiss=dw*dup*dv*dst
      aissi=dwp*du*dv*dstp
      aisss=dwp*du*dv*dst
      assii=dw*du*dvp*dstp
      assis=dw*du*dvp*dst
      asssi=dw*du*dv*dstp
      assss=dw*du*dv*dst
      if((istinf.lt.1).or.(iwinf.lt.1).or.(iuinf.lt.1).or.(ivinf.lt.1)
     &  .or.(iusup.gt.nxhyx).or.(ivsup.gt.nyhyx).or.(istsup.gt.nzhyx)
     &  .or.(iwsup.gt.ntauhyx)) then
         write(8,*) 'point out of mesh'
         write(8,*) tb,xx,yy,zz
         write(8,*) xmin,ymin,zmin,tbmin
         write(8,*) xmax,ymax,zmax,tbmax
         close(8)
         stop
      endif
      call gethyval(istinf,iwinf,iuinf,ivinf,epepsiiii,eptemiiii,pss,
     &     v1iiii,v2iiii,v3iiii)
      call gethyval(istinf,iwinf,iuinf,ivsup,epepsiisi,eptemiisi,pss,
     &     v1iisi,v2iisi,v3iisi)
      call gethyval(istinf,iwinf,iusup,ivinf,epepsisii,eptemisii,pss,
     &     v1isii,v2isii,v3isii)
      call gethyval(istinf,iwsup,iuinf,ivinf,epepssiii,eptemsiii,pss,
     &     v1siii,v2siii,v3siii)
      call gethyval(istinf,iwsup,iuinf,ivsup,epepssisi,eptemsisi,pss,
     &     v1sisi,v2sisi,v3sisi)
      call gethyval(istinf,iwinf,iusup,ivsup,epepsissi,eptemissi,pss,
     &     v1issi,v2issi,v3issi)
      call gethyval(istinf,iwsup,iusup,ivinf,epepsssii,eptemssii,pss,
     &     v1ssii,v2ssii,v3ssii)
      call gethyval(istinf,iwsup,iusup,ivsup,epepssssi,eptemsssi,pss,
     &     v1sssi,v2sssi,v3sssi)
      call gethyval(istsup,iwinf,iuinf,ivinf,epepsiiis,eptemiiis,pss,
     &     v1iiis,v2iiis,v3iiis)
      call gethyval(istsup,iwinf,iuinf,ivsup,epepsiiss,eptemiiss,pss,
     &     v1iiss,v2iiss,v3iiss)
      call gethyval(istsup,iwinf,iusup,ivinf,epepsisis,eptemisis,pss,
     &     v1isis,v2isis,v3isis)
      call gethyval(istsup,iwsup,iuinf,ivinf,epepssiis,eptemsiis,pss,
     &     v1siis,v2siis,v3siis)
      call gethyval(istsup,iwsup,iuinf,ivsup,epepssiss,eptemsiss,pss,
     &     v1siss,v2siss,v3siss)
      call gethyval(istsup,iwinf,iusup,ivsup,epepsisss,eptemisss,pss,
     &     v1isss,v2isss,v3isss)
      call gethyval(istsup,iwsup,iusup,ivinf,epepsssis,eptemssis,pss,
     &     v1ssis,v2ssis,v3ssis)
      call gethyval(istsup,iwsup,iusup,ivsup,epepsssss,eptemssss,pss,
     &     v1ssss,v2ssss,v3ssss)
        enert=aiiii*epepsiiii+aiisi*epepsiisi+aisii*epepsisii
     &       +asiii*epepssiii+asisi*epepssisi+aissi*epepsissi
     &       +assii*epepsssii+asssi*epepssssi+aiiis*epepsiiis
     &       +aiiss*epepsiiss+aisis*epepsisis+asiis*epepssiis
     &       +asiss*epepssiss+aisss*epepsisss+assis*epepsssis
     &       +assss*epepsssss
        temp=aiiii*eptemiiii+aiisi*eptemiisi+aisii*eptemisii
     &       +asiii*eptemsiii+asisi*eptemsisi+aissi*eptemissi
     &       +assii*eptemssii+asssi*eptemsssi+aiiis*eptemiiis
     &       +aiiss*eptemiiss+aisis*eptemisis+asiis*eptemsiis
     &       +asiss*eptemsiss+aisss*eptemisss+assis*eptemssis
     &       +assss*eptemssss
        velx=aiiii*v1iiii+aiisi*v1iisi+aisii*v1isii
     &       +asiii*v1siii+asisi*v1sisi+aissi*v1issi
     &       +assii*v1ssii+asssi*v1sssi+aiiis*v1iiis
     &       +aiiss*v1iiss+aisis*v1isis+asiis*v1siis
     &       +asiss*v1siss+aisss*v1isss+assis*v1ssis
     &       +assss*v1ssss
        vely=aiiii*v2iiii+aiisi*v2iisi+aisii*v2isii
     &       +asiii*v2siii+asisi*v2sisi+aissi*v2issi
     &       +assii*v2ssii+asssi*v2sssi+aiiis*v2iiis
     &       +aiiss*v2iiss+aisis*v2isis+asiis*v2siis
     &       +asiss*v2siss+aisss*v2isss+assis*v2ssis
     &       +assss*v2ssss
        velz=aiiii*v3iiii+aiisi*v3iisi+aisii*v3isii
     &       +asiii*v3siii+asisi*v3sisi+aissi*v3issi
     &       +assii*v3ssii+asssi*v3sssi+aiiis*v3iiis
     &       +aiiss*v3iiss+aisis*v3isis+asiis*v3siis
     &       +asiss*v3siss+aisss*v3isss+assis*v3ssis
     &       +assss*v3ssss
         if(abs(velx)>1.0d0)then
           write(*,*) tb,xx,yy,zz, 'velx= ', velx
         else if(abs(vely)>1.0d0)then
           write(*,*) tb,xx,yy,zz, 'vely= ', vely
         else if(abs(velz)>1.0d0)then
           write(*,*) tb,xx,yy,zz, 'velz= ', velz
          end if
      ener=enert
      ifail=0
      return
      end 


      subroutine azhydrogrid(tb,xx,yy,temp,velx,vely,ener,ifail)
      use PBGhydro
      implicit none
      integer ifail,iuinf,iusup,ivinf,ivsup,iwinf,iwsup
      real tb,temp,velx,vely,ener,enert,aii,asi,ass,ais,bii,bsi,bss,
     &  bis,cii,csi,css,cis
      real u,du,dup,v,dv,dvp,w,dw,dwp
      double precision xx,yy
      real aiii,aiis,aisi,asii,asis,assi,aiss,asss
      ifail=-1
      if(tb.lt.tbmin) return
      temp=0.
      ener=0.
      velx=0.
      vely=0.
      ifail=1
      if((tb.gt.tbmax).or.(dabs(xx).gt.xmax).or.(dabs(yy).gt.ymax))then
         return
      endif
******weights for all cases******
      u=dabs(xx)/xstep
      iuinf=int(u)
      iusup=iuinf+1
      du=u-iuinf
      dup=1-du
      v=dabs(yy)/ystep
      ivinf=int(v)
      ivsup=ivinf+1
      dv=v-ivinf
      dvp=1-dv
      w=(tb-tbmin)/tbstep
      iwinf=int(w)
      iwsup=iwinf+1
      dw=w-iwinf
      dwp=1-dw
      aii=dup*dvp
      asi=du*dvp
      ass=du*dv
      ais=dup*dv
      bii=dwp*dup
      bsi=dw*dup
      bss=dw*du
      bis=dwp*du
      cii=dwp*dvp
      csi=dw*dvp
      css=dw*dv
      cis=dwp*dv
      aiii=dwp*dup*dvp
      aiis=dwp*dup*dv
      aisi=dwp*du*dvp
      asii=dw*dup*dvp
      asis=dw*dup*dv
      aiss=dwp*du*dv
      assi=dw*du*dvp
      asss=dw*du*dv
********************************
      if(tb.eq.tbmax) then
         if(dabs(xx).eq.xmax)then
            if(dabs(yy).eq.ymax)then
               enert=eps(nbtb,nbx,nby)
               temp=tem(nbtb,nbx,nby)
               velx=vx(nbtb,nbx,nby)
               vely=vy(nbtb,nbx,nby)
            else
               enert=dv*eps(nbtb,nbx,ivsup)+dvp*eps(nbtb,nbx,ivinf)
               temp=dv*tem(nbtb,nbx,ivsup)+dvp*tem(nbtb,nbx,ivinf)
               velx=dv*vx(nbtb,nbx,ivsup)+dvp*vx(nbtb,nbx,ivinf)
               vely=dv*vy(nbtb,nbx,ivsup)+dvp*vy(nbtb,nbx,ivinf)
            endif
         else
            if(dabs(yy).eq.ymax) then
               enert=du*eps(nbtb,iusup,nby)+dup*eps(nbtb,iuinf,nby)
               temp=du*tem(nbtb,iusup,nby)+dup*tem(nbtb,iuinf,nby)
               velx=du*vx(nbtb,iusup,nby)+dup*vx(nbtb,iuinf,nby)
               vely=du*vy(nbtb,iusup,nby)+dup*vy(nbtb,iuinf,nby)
            else
               enert=aii*eps(nbtb,iuinf,ivinf)
     &               +asi*eps(nbtb,iusup,ivinf)+
     &              ais*eps(nbtb,iuinf,ivsup)
     &              +ass*eps(nbtb,iusup,ivsup)
               temp=aii*tem(nbtb,iuinf,ivinf)
     &              +asi*tem(nbtb,iusup,ivinf)+
     &              ais*tem(nbtb,iuinf,ivsup)
     &              +ass*tem(nbtb,iusup,ivsup)
               velx=aii*vx(nbtb,iuinf,ivinf)
     &              +asi*vx(nbtb,iusup,ivinf)+
     &              ais*vx(nbtb,iuinf,ivsup)
     &              +ass*vx(nbtb,iusup,ivsup)
               vely=aii*vy(nbtb,iuinf,ivinf)
     &              +asi*vy(nbtb,iusup,ivinf)+
     &              ais*vy(nbtb,iuinf,ivsup)
     &              +ass*vy(nbtb,iusup,ivsup)
            endif
         endif
      else
         if(dabs(xx).eq.xmax) then
            if(dabs(yy).eq.ymax) then
               enert=dw*eps(iwsup,nbx,nby)+dwp*eps(iwinf,nbx,nby)
               temp=dw*tem(iwsup,nbx,nby)+dwp*tem(iwinf,nbx,nby)
               velx=dw*vx(iwsup,nbx,nby)+dwp*vx(iwinf,nbx,nby)
               vely=dw*vy(iwsup,nbx,nby)+dwp*vy(iwinf,nbx,nby)
            else
               enert=cii*eps(iwinf,nbx,ivinf)+csi*eps(iwsup,nbx,ivinf)+
     &              cis*eps(iwinf,nbx,ivsup)+css*eps(iwsup,nbx,ivsup)
               temp=cii*tem(iwinf,nbx,ivinf)+csi*tem(iwsup,nbx,ivinf)+
     &              cis*tem(iwinf,nbx,ivsup)+css*tem(iwsup,nbx,ivsup)
               velx=cii*vx(iwinf,nbx,ivinf)+csi*vx(iwsup,nbx,ivinf)+
     &              cis*vx(iwinf,nbx,ivsup)+css*vx(iwsup,nbx,ivsup)
               vely=cii*vy(iwinf,nbx,ivinf)+csi*vy(iwsup,nbx,ivinf)+
     &              cis*vy(iwinf,nbx,ivsup)+css*vy(iwsup,nbx,ivsup)
            endif
         else
            if(dabs(yy).eq.ymax) then 
               enert=bii*eps(iwinf,iuinf,nby)+bsi*eps(iwsup,iuinf,nby)+
     &              bis*eps(iwinf,iusup,nby)+bss*eps(iwsup,iusup,nby)
               temp=bii*tem(iwinf,iuinf,nby)+bsi*tem(iwsup,iuinf,nby)+
     &              bis*tem(iwinf,iusup,nby)+bss*tem(iwsup,iusup,nby)
               velx=bii*vx(iwinf,iuinf,nby)+bsi*vx(iwsup,iuinf,nby)+
     &              bis*vx(iwinf,iusup,nby)+bss*vx(iwsup,iusup,nby)
               vely=bii*vy(iwinf,iuinf,nby)+bsi*vy(iwsup,iuinf,nby)+
     &              bis*vy(iwinf,iusup,nby)+bss*vy(iwsup,iusup,nby)
            else
           enert=aiii*eps(iwinf,iuinf,ivinf)+aiis*eps(iwinf,iuinf,ivsup)
     &          +aisi*eps(iwinf,iusup,ivinf)+asii*eps(iwsup,iuinf,ivinf)
     &          +asis*eps(iwsup,iuinf,ivsup)+aiss*eps(iwinf,iusup,ivsup)
     &          +assi*eps(iwsup,iusup,ivinf)+asss*eps(iwsup,iusup,ivsup)
           temp=aiii*tem(iwinf,iuinf,ivinf)+aiis*tem(iwinf,iuinf,ivsup)
     &          +aisi*tem(iwinf,iusup,ivinf)+asii*tem(iwsup,iuinf,ivinf)
     &          +asis*tem(iwsup,iuinf,ivsup)+aiss*tem(iwinf,iusup,ivsup)
     &          +assi*tem(iwsup,iusup,ivinf)+asss*tem(iwsup,iusup,ivsup)
           velx=aiii*vx(iwinf,iuinf,ivinf)+aiis*vx(iwinf,iuinf,ivsup)
     &          +aisi*vx(iwinf,iusup,ivinf)+asii*vx(iwsup,iuinf,ivinf)
     &          +asis*vx(iwsup,iuinf,ivsup)+aiss*vx(iwinf,iusup,ivsup)
     &          +assi*vx(iwsup,iusup,ivinf)+asss*vx(iwsup,iusup,ivsup)
            vely=aiii*vy(iwinf,iuinf,ivinf)+aiis*vy(iwinf,iuinf,ivsup)
     &          +aisi*vy(iwinf,iusup,ivinf)+asii*vy(iwsup,iuinf,ivinf)
     &          +asis*vy(iwsup,iuinf,ivsup)+aiss*vy(iwinf,iusup,ivsup)
     &          +assi*vy(iwsup,iusup,ivinf)+asss*vy(iwsup,iusup,ivsup)
            endif
         endif
      endif
      ener=enert
      ifail=0
      return
      end

      double precision function temprofil(tb,x,y)
      double precision x,y,pi
      real tb
      parameter(pi=3.14159265358979323844d0)
      temprofil=(1/(2*pi*(0.12*tb+0.4)))
     &     *dexp(-.5*(x**2/(5.31*tb+2.61)**2
     &     +y**2/(4.8*tb+5.1)**2))
      end
      
      double precision function enerprofil(tb,x,y)
      double precision x,y,pi
      real tb
      parameter(pi=3.14159265358979323844d0)
      enerprofil=(1/(2*pi*(0.0294*tb+0.0484)))*dexp(-.5*(x**2/(0.145*tb+
     &     2.81)**2+y**2/(-0.012*tb+3.95)**2))
      end
