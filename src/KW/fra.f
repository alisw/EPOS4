C
C  This file is part of EPOS4
C  Copyright (C) 2022 research institutions and authors (See CREDITS file)
C  This file is distributed under the terms of the GNU General Public License version 3 or later
C  (See COPYING file for the text of the licence)
C

c-----------------------------------------------------------------------
      subroutine gakfra(iret)
c-----------------------------------------------------------------------

#include "aaa.h"
#include "ems.h"
#include "sem.h"
      parameter (eta=0.4,etap=0.1)
      parameter (mxnbr=800,mxnba=10000)
      double precision pst,pf
      common /gag/pst(4,0:mxnba),nob,ipst(0:mxnba)
     $     ,nbr,ijb(2,0:mxnbr),xyb(2,0:mxnbr)
     &     ,ptb(4,0:mxnbr),iflb(0:mxnbr),ip(0:mxnbr),co,x,y
c      common /cpptr/ pptr(4,mxptl),ystr(mxptl)
      double precision p1,p12,p2
      dimension co(12),ind(0:mxnbr),p1(5),p2(5)
      dimension ic2(2),ic1(2),ic(2)!,fkap(3)
      dimension xs(3),xori(4)
      common/pbk/pbkl(3) /cnsbp/nsbp  /cn8ptl/n8ptl !bg
      common/czz/kky,krm,kdrm,kzrm,pui(3),pdi(3),psi(3),pci(3),pbi(3) !bg
     &,zzzz,itp,pduui(3),pdudi(3),pdusi(3),pduci(3),pdubi(3),pdddi(3) !bg
     &,pddsi(3),pddci(3),pddbi(3),pdssi(3),pdsci(3),pdsbi(3),pdcci(3) !bg
     &,pdcbi(3),pdbbi(3),fkap(3),kkc !bg
      double precision psg,utdble
      common/zpsg/psg(5)
      logical go
                   
      pf(i,j)=(pst(4,i)-pst(3,i))*(pst(4,j)+pst(3,j))
     &     -pst(2,i)*pst(2,j)-pst(1,i)*pst(1,j)

      nob=0
      call utpri('gakfra',ish,ishini,4)
      iret=0
c.....search string in pptl-list...................................
      nkmax=nptl
      nk1=maproj+matarg+1
 2    continue
      if(nclean.gt.0.and.nptl.gt.mxptl/2)then
        nnnpt=0
        do iii=maproj+matarg+1,nptl
          go=.true.
          if(nclean.eq.1.and.istptl(iii).le.istmax)go=.false.
          if(go.and.mod(istptl(iii),10).ne.0)then
            istptl(iii)=99
            nnnpt=nnnpt+1
          endif
        enddo
        if(nnnpt.gt.mxptl-nptl)then
          nptl0=nptl
          call utclea(maproj+matarg+1,nptl0)
          nkmax=nptl
          nk1=maproj+matarg+1
c          write(ifch,'(a)')' ---------after utclea------------'
c          call blist('&',nk1,nkmax)
        endif
      endif

      do while(nk1.le.nkmax)
        if(istptl(nk1).eq.20)then
          nk2=nk1+1
          do while(idptl(nk2).eq.9)
            if(istptl(nk2).ne.20)then
              print*,'gakfra parton chain wrong (ist=20 expected)'
              do nkjj=nk1,nk2
              print*,'n id ist:',nkjj, idptl(nkjj),istptl(nkjj)
              enddo
              stop 'ERROR231027c'
            endif
            nk2=nk2+1
          enddo
          goto 33                !ok, string from nk1 to nk2
        else
          if(istptl(nk1).eq.22)write(ifmt,'(a)')'Skip 22 parton'
          nk1=nk1+1
        endif
      enddo

      goto 9999                 !no more string


c.......................................................................
c                       decay string  nk1 - nk2
c.......................................................................


 33   continue

      if(iotst2.gt.0)then
        fpost = float(iotst2)/100.
        rrr=sqrt( -2*log(rangen()) ) /sqrt(2.00)  *fpost
         phi=2*pi*rangen()
         xorptl(1,nk1)=xorptl(1,nk1)+rrr*cos(phi)
         xorptl(2,nk1)=xorptl(2,nk1)+rrr*sin(phi)
      endif
      sora=0.035 !small spread for string origin
      xra=sqrt(-2*log(rangen()))*cos(2*pi*rangen())*sora
      yra=sqrt(-2*log(rangen()))*cos(2*pi*rangen())*sora
      xori(1)=xorptl(1,nk1)+xra
      xori(2)=xorptl(2,nk1)+yra
      xori(3)=iostro*xorptl(3,nk1)  
      xori(4)=iostro*xorptl(4,nk1)  
      !print*,'xori',nk1,xori(1),xori(2),ityptl(nk1)

      nsde=0 
      nptlini=nptl

  3   continue
      nsde=nsde+1
      !if(mod(nsde,100000).eq.0)then
        !write(ifmt,*)'Redo string decay',nk1,nk2,nsde
        !write(ifch,*)'Gakfra : redo event (2)'
        !iret=1
        !goto 9999
      !endif

      nob=-1
      if(ish.ge.3)then
        write(ifch,'(/1x,25a1/5x,a/1x,25a1/)')
     *  ('-',k=1,25),'string decay',('-',k=1,25)
        write(ifch,'(a)')' ---------string------------'
        call blist('&',nk1,nk2)
      endif
      
      if(abs(idptl(nk1)).ge.4.or.abs(idptl(nk2)).ge.4)then
        icharm=1
      else
        icharm=0
      endif

      itp=ityptl(nk1)
      kkc=0
c      kky=0
c      if(nk2-nk1.gt.1.or.iappl.eq.6)kky=1
      if(jorptl(nk1).ne.0.or.jorptl(nk2).ne.0)then
        kky=1
        if((abs(jorptl(nk1)).eq.7.or.abs(jorptl(nk1)).eq.8)
     *.and.(abs(jorptl(nk2)).eq.7.or.abs(jorptl(nk2)).eq.8)
     *       .and.nk2-nk1.gt.1)itp=itp+1
        if(abs(jorptl(nk1)).ge.4.or.abs(jorptl(nk2)).ge.4)kkc=1
      else
        kky=0
      endif
      krm=0
      if((itp.ge.50.and.itp.le.59).or.(itp.ge.40.and.itp.le.49))krm=1
      kdrm=0
      if(itp.eq.43.or.itp.eq.53)kdrm=1
c Excited remnant due to split
      kzrm=0
      if(itp.eq.44.or.itp.eq.54.or.itp.eq.46.or.itp.eq.56)kzrm=1

      do k=1,5
        p1(k)=0d0
      enddo
      x=q2sft
      iorrem=iorptl(nk1)
      if(iorrem.ne.0)pmior=pptl(5,iorrem)
      do i=nk1,nk2
        x=max(x,qsqptl(i)) !virtuality (Q2)
        do k=1,5
          !---------------
          a=pptl(k,i)
          if(.not.(a.gt.0..or.a.le.0.) !NaN catch
     .   .or.a.gt.1e35 .or. a.lt.-1e35 )then !Infinity catch
            write(ifmt,'(a,2i5,$)')
     .      'WARNING gakfra (NaN,Inf)',idptl(i),k
            write(ifmt,*)a
            iret=1
            return
          endif
          !---------------
          p2(k)=utdble(pptl(k,i))
        enddo
        p2(4)=sqrt(p2(3)**2+p2(2)**2+p2(1)**2+p2(5)**2)
        do k=1,4
          p1(k)=p1(k)+p2(k)
        enddo
      enddo
      xq2s=1+exp(-x/100.)
      p1(5)=(p1(4)-p1(3))*(p1(4)+p1(3))-p1(2)**2-p1(1)**2
      if(p1(5).gt.0.d0)then
        p1(5)=sqrt(p1(5))
        if(iorrem.eq.0)pmior=sngl(p1(5))
      else
        write(ifmt,'(a)')'WARNING gakfra String m2 < 0'
c        write(ifmt,'(a,5e10.3,$)')'5-mom str = '
c     .  ,p1(1),p1(2),p1(3),p1(4),p1(5)
c        do i=nk1,nk2
c          write(ifmt,'(a,5e10.3,$)')'  parton = ',(pptl(k,i),k=1,5)
c        enddo
        if(p1(5).lt.-1d-4)then
          !return only if it is not a problem due to limited precision (very high momentum)
          write(ifmt,'(a,e10.3,a)')
     .    'WARNING gakfra (negative m2)',p1(5),'  return'
          iret=1
          return
        endif 
        p1(5)=1.d0    !not to crash for the boost
        !mass known only for soft strings
        if(iorrem.ne.0.and.
     &    (ityptl(i).eq.42.or.ityptl(i).eq.52.or.ityptl(i)/10.eq.2))then
          if(ish.ge.1)
     &    write(ifch,*)'string mass : ',ityptl(iorrem),pptl(5,iorrem)
          p1(5)=utdble(pptl(5,iorrem))
        endif
        if(iorrem.eq.0)pmior=1e30
      endif
      do k=1,5
        psg(k)=p1(k)
      enddo

c.....pbreak ................................

c      if(iappl.eq.6)then
c        x=engy               !energy density = energy for (almost) point like particle renormalized to use energy (phenomenological radius "zipinc" used to define zpaptl for hadronic collisions)
c      elseif(zipinc.gt.0)then
c      else
c        x=p1(5)              !sub string mass (energy in cms)
c      endif
      pb0i=pbreak

c      if(pb0i.gt.0.d0.or.(pb0i.gt.-0.99999d0.and.krm.eq.1))then 
c      if(pb0i.gt.0.d0.or.(pb0i.gt.-0.99999d0.and.kky.eq.0))then
      if(pb0i.gt.0.d0)then
c around 0.4 for soft because of pi0 spectra for CR and pp multiplicity
c and pt spectra at RHIC
        pbk=abs(pb0i) !bg
      else
c fit is important not only in e+e- but for all particle production
        pb0=abs(pb0i)

c        pbk=min(pb0,0.065 + 3600./x**3 - 400./x**2 + 17./x)   !<--------
c       pbk=0.14 + 855./min(max(14.,x),91.2)**3 
c     &            - 111.1/min(max(14.,x),91.2)**2 
c     &            + 7.16/min(max(14.,x),91.2)   !<-------- hajo
ckw17        pbk=pb0+(1.-exp(-x/50.))*(pbreakg-pb0) !ckw17 old method
        pbk=pbreakg+(pb0-pbreakg)/max(1.,x/q2sft)**zipinc !ckw17 new method 
        fcharm=1.1
        if(icharm.eq.1)pbk=pbk*fcharm
c        pbk=pb0+(1.-exp(-sqrt(x)/zipinc))*(pbremxnbaakg-pb0)
c        pbk=pb0+(1.-exp(-sqrt(x)/10.))*(pbreakg-pb0) !bg

c       {{14,0.25},{22,0.22},{34,0.22},{91.2,0.15}} !new
      endif

      fbreakgg = 1. !1.5  !correction factor of pbreak in gg bands

      pbkl(1)=pbk                 
      pbkl(3)=pbk * fbreakgg           
      pbkl(2)=0.5*(pbkl(1)+pbkl(3))  

c      pbkl(1)=abs(pbreak)+abs(pbreakg)*x
c      pbkl(1)=abs(pbreak)
c      pbkl(3)=abs(pbreakg)
c      pbkl(2)=0.5*(pbkl(1)+pbkl(3))
c      if(kky.eq.1)pbkl(1)=pbreakg

      fkap(1)=fkappa
      fkap(3)=fkappag
      fkap(2)=0.5*(fkap(1)+fkap(3))

c      write(ifch,*)'pbk',kky,krm,kdrm,x,pbk,p1(5)
c     .                  ,zpaptl(1,nk1),zpaptl(2,nk1)

      zzzz=1.
      if(ish.ge.4)write(ifch,*)"String parameters:"
c      print *
     &                         ,itp,krm,kdrm,kzrm,fkap,zzzz,pbk,x !bg
     &                         ,zpaptl(1,nk1),zpaptl(2,nk1)
      if(fkap(1).le.0..or.zzzz.le.0.)
     &   call utstop("fkappa<=0. not allowed !&")

      pmqqx=pmqq+zopinc*log(max(1.,x/q2sft)) 
      pmqu2=pmqu**2
      difud=                       pmqd**2-pmqu2
      difus=                       pmqs**2-pmqu2
      difuuu=(          pmqqx+pmqu+pmqu)**2-pmqu2
      difuud=(pudd*     pmqqx+pmqd+pmqu)**2-pmqu2
      difuus=(puds*     pmqqx+pmqs+pmqu)**2-pmqu2
      difudd=(pudd*pudd*pmqqx+pmqd+pmqd)**2-pmqu2
      difuds=(pudd*puds*pmqqx+pmqs+pmqd)**2-pmqu2
      difuss=(puds*puds*pmqqx+pmqs+pmqs)**2-pmqu2
      if(nrflav.gt.3)then
       difuc=                       pmqc**2-pmqu2
       difuuc=(pudc*     pmqqx+pmqc+pmqu)**2-pmqu2
       difudc=(pudd*pudc*pmqqx+pmqc+pmqd)**2-pmqu2
       difusc=(puds*pudc*pmqqx+pmqc+pmqs)**2-pmqu2
       difucc=(pudc*pudc*pmqqx+pmqc+pmqc)**2-pmqu2 !bg
      endif

c  ig=1  qq bands
c  ig=2  qg bands
c  ig=3  gg bands

      fschwisgg=1.0 !0.3 !correction factor in Schwinger formula for s-quarks in gg bands
      xs(1)=1          !*xq2s
      xs(3)=fschwisgg  !*xq2s
      xs(2)=0.5*(xs(1)+xs(3)) 

      do ig=1,3
      pui(ig)=1.
      pdi(ig)=pui(ig)*exp(-pi*difud/fkap(ig)) !assymmetric u and d relevant because of assymmetry Kch / K0 in e+e-
      psi(ig)=pui(ig)*exp(-pi*difus*xs(ig)/fkap(ig))
      pduui(ig)=pui(ig)*exp(-pi*difuuu/fkap(ig))
      pdudi(ig)=pui(ig)*exp(-pi*difuud/fkap(ig))
      pdusi(ig)=pui(ig)*exp(-pi*difuus/fkap(ig))
      pdddi(ig)=pui(ig)*exp(-pi*difudd/fkap(ig))
      pddsi(ig)=pui(ig)*exp(-pi*difuds/fkap(ig))
      pdssi(ig)=pui(ig)*exp(-pi*difuss/fkap(ig))
      if(nrflav.gt.3)then
        pci(ig)=pui(ig)*exp(-pi*difuc/fkap(ig))
        pduci(ig)=pui(ig)*exp(-pi*difuuc/fkap(ig))
        pddci(ig)=pui(ig)*exp(-pi*difudc/fkap(ig))
        pdsci(ig)=pui(ig)*exp(-pi*difusc/fkap(ig))
        pdcci(ig)=pui(ig)*exp(-pi*difucc/fkap(ig))
        pdcbi(ig)=0. !pui(ig)*exp(-pi*difucb/fkap(ig)) !tp
        pbi(ig)=0.   !pui(ig)*exp(-pi*difub/fkap(ig))  !tp
        pdubi(ig)=0. !pui(ig)*exp(-pi*difuub/fkap(ig)) !tp
        pddbi(ig)=0. !pui(ig)*exp(-pi*difudb/fkap(ig)) !tp
        pdsbi(ig)=0. !pui(ig)*exp(-pi*difusb/fkap(ig)) !tp
        pdbbi(ig)=0. !pui(ig)*exp(-pi*difubb/fkap(ig)) !tp
      else
        pci(ig)=0.
        pduci(ig)=0.
        pddci(ig)=0.
        pdsci(ig)=0.
        pdcci(ig)=0.
        pdcbi(ig)=0. !bg
        pbi(ig)=0.   !bg
        pdubi(ig)=0. !bg
        pddbi(ig)=0. !bg
        pdsbi(ig)=0. !bg
        pdbbi(ig)=0. !bg
      endif
      enddo

c.....light like ...............................................

      if(ish.ge.6) call blist("before light like&",nk1,nk2)
      if(rangen().lt.0.5)then
        i1=nk1
        i2=nk1+1
      else
        i1=nk2-1
        i2=nk2
      endif
      ii=1
      do while(ii.ge.0)
        do i=1,4
          p2(i)=utdble(pptl(i,i1))+utdble(pptl(i,i2))
        enddo
c       p2(5)=sqrt(p2(4)**2-p2(3)**2-p2(2)**2-p2(1)**2)
        am1=pptl(5,i1)**2
        am2=pptl(5,i2)**2
        am12=max(0d0,
     &          (p2(4)+p2(3))*(p2(4)-p2(3))-p2(2)**2-p2(1)**2-am1-am2)/2
        if(am12**2.gt.am1*am2)then
          ak1=.5*((am2+am12)/sqrt(am12**2-am1*am2))-.5
          ak2=.5*((am1+am12)/sqrt(am12**2-am1*am2))-.5
        else
          ak1=0.
          ak2=0.
        endif
        if(ish.ge.6)write(ifch,'(a,2i4,9f12.6)') 'overlaps'
     $       ,i1,i2,ak1,ak2,sqrt(2.*am12)
     &       ,am1,am2
c test energy first
        j=4
        a1=(1.+ak1)*pptl(j,i1)-ak2*pptl(j,i2)
        a2=(1.+ak2)*pptl(j,i2)-ak1*pptl(j,i1)
c rescale momenta only if energy is positive
        if(a1.gt.0..and.a2.gt.0.)then
          pptl(j,i1)=a1
          pptl(j,i2)=a2
          do j=1,3
            a1=(1.+ak1)*pptl(j,i1)-ak2*pptl(j,i2)
            a2=(1.+ak2)*pptl(j,i2)-ak1*pptl(j,i1)
            pptl(j,i1)=a1
            pptl(j,i2)=a2
          enddo
        endif
        pptl(5,i1)=0.
        pptl(5,i2)=0.
        if(nk2-nk1.eq.1) ii=-1
        ii=ii-1
        if (nk1.eq.i1)then
          i1=nk2-1
          i2=nk2
        else
          i1=nk1
          i2=nk1+1
        endif
      enddo

c.....inverse  ...............................................

c      if(nk2.eq.nk1+1.and.nptl.lt.mxptl)then
c        call utrepla(nptl+1,nk2)
c        call utrepla(nk2,nk1)
c        call utrepla(nk1,nptl+1)
c      endif

      if(ish.ge.6) call blist("after  light like&",nk1,nk2)

c.....on-shell and copy ...............................................

      if(ish.ge.6) write (ifch,*) 'on-shell check'
      do i=nk1,nk2
        do k=1,5
          p2(k)=utdble(pptl(k,i))
        enddo
        p12=p2(4)
        p2(4)=sqrt(p2(3)**2+p2(2)**2+p2(1)**2+p2(5)**2)
        if(abs(p12-p2(4))/max(.1d0,p12,p2(4)).ge.1e-1.and.ish.ge.2)then
          write(ifch,*)'warning: on-shell problem:'
     &           ,i, idptl(i),(pptl(k,i),k=1,4),' (',sngl(p2(4)),') '
     &           ,pptl(5,i), istptl(i)
        endif
        if(ish.ge.6)  write (ifch,*)  p12 ,'->',p2(4)
        call utlob2(1,p1(1),p1(2),p1(3),p1(4),p1(5)
     $       ,p2(1),p2(2),p2(3),p2(4),60)
        f=1.0
        if(i.ne.nk1.and.i.ne.nk2)f=0.5
        nmax=1
        ff=1.
        aemax=1000.              ! here max band mass
c          write(ifch,*)"test",f,p2(4),aemax,f*p2(4).gt.aemax
        if(f*p2(4).gt.aemax) then
          nmax = int(f*p2(4)/aemax)+1
          ff=1./float(nmax)
c          print *,"nmax",nmax
c          write(ifch,*)"nmax",nmax
        endif
        nn=1
        do while(nn.le.nmax)
          f=.5
          if(i.eq.nk1.and.nn.eq. 1  ) f=1.
          if(i.eq.nk2.and.nn.eq.nmax) f=1.
          nob=nob+1
          if(nob.gt.mxnba-7)stop'gakfra: increase parameter mxnba '
c         if(nmax.ge.2) print *, nob,ff,f,i
c         if(nmax.ge.2) write(ifch,*)'nob',nob,ff,f,i
          do k=1,4
            pst(k,nob)=dble(ff*f)*p2(k)
            ipst(nob)=i   !parton index associated to band nob 
          enddo
          nn=nn+1
        enddo
      enddo                     !i

      do i1=nob-1,1,-1          ! copy again gluons
        nob=nob+1
        if(nob.gt.mxnba-7)stop'gakfra: increase parameter mxnba '
        do k=1,4
          pst(k,nob)=pst(k,i1)
          ipst(nob)=ipst(i1)
        enddo
      enddo
      nob=nob+1                 ! nob is number of kinks
      if(nob.gt.mxnba-7)stop'gakfra: increase parameter mxnba '
      if(ish.ge.6) then
        write (ifch,*) 'bands:'
        do i=0,nob-1
          write (ifch,'(i4,4g20.13)') i,(pst(k,i),k=1,4)
        enddo
      endif

c.....total momentum...............................................
      do k=1,4
        pst(k,nob)=0.d0
      enddo
      do i=0,nob-1
        do k=1,4
          pst(k,nob)=pst(k,nob)+pst(k,i)
        enddo
      enddo



c........... incement nptl and write total string on /cptl/ ............
c            ^^^^^^^^^^^^^      
      nptl=nptl+1
      if(nptl.gt.mxptl)call utstop('gakfra: mxptl too small ...&')
      call idtr5(idptl(nk1),ic1)
      call idtr5(idptl(nk2),ic2)
      ic(1)=ic1(1)+ic2(1)
      ic(2)=ic1(2)+ic2(2)
      idptl(nptl) = 8*10**8+ ic(1)*100 + ic(2)/100 !nbr+1
      istptl(nptl)=10
      qsqptl(nptl)=x
      zpaptl(1,nptl)=zpaptl(1,nk1)
      zpaptl(2,nptl)=zpaptl(2,nk1)
      n8ptl=nptl
      do i=1,5
        pptl(i,nptl)=p1(i)
      enddo
      if(ish.ge.3)then
        write(ifch,'(a)')' ---------total string------------'
        call blist('&',nptl,nptl)
      endif

c....................................................................
      ijb(1,0)=-1
      ijb(2,0)=-1
      ijb(1,1)=-1
      ijb(2,1)=-1
      iflb(0)=-idptl(nk1)
      iflb(1)= idptl(nk2)
      iflb(2)= 0               !to avoid misusage in HQ strings 
      do i=1,4
        ptb(i,0)=0.
        ptb(i,1)=0.
      enddo

c...................light string....................................
      amm=ammin(idptl(nk1),idptl(nk2))


      if(sqrt(max(0.d0,pf(nob,nob))).lt.amm)then
        id=idsp(  idptl(nk1),idptl(nk2) )
        if(ish.ge.1)then
          write (ifch,*)
     .       '-------string too light to decay--------'
          write (ifch,*) id,p1(5),amm,nk1,nk2
          write (ifch,*) id, sqrt(max(0.d0,pf(nob,nob))) ,amm,nptl+1
        endif
        if(id.ne.0)then
        am=sqrt(max(0.d0,pf(nob,nob)))
        call idress(id,am,idr,iadj)
        if(ish.ge.1)write (ifch,*) id,nptl+1,am,idr,iadj
        nptl=nptl+1
        if(nptl.gt.mxptl)
     &       call utstop('gakfra (light): mxptl too small ...&')
        do i=1,5
          pptl(i,nptl)=p1(i)
        enddo
        if(iadj.eq.1)pptl(5,nptl)=am   !otherwise it will be rejected at decay
        do i=nk1,nk2
          istptl(i)=21
          ifrptl(1,i)=n8ptl
          ifrptl(2,i)=0
        enddo
        istptl(n8ptl)=29        !code for fragmented string
        ityptl(n8ptl)=itp !ityptl(nk1)
        itsptl(n8ptl)=itsptl(nk1)
        rinptl(n8ptl)=-9999
        iorptl(n8ptl)=nk1
        jorptl(n8ptl)=nk2
        ifrptl(1,n8ptl)=n8ptl+1
        ifrptl(2,n8ptl)=nptl
        ityptl(nptl)=itp !ityptl(nk1)
        itsptl(nptl)=itsptl(nk1)
        rinptl(nptl)=-9999
        istptl(nptl)=0
        iorptl(nptl)=n8ptl
        jorptl(nptl)=0
        ifrptl(1,nptl)=0
        ifrptl(2,nptl)=0
        qsqptl(nptl)=qsqptl(n8ptl)
        zpaptl(1,nptl)=zpaptl(1,n8ptl)
        zpaptl(2,nptl)=zpaptl(2,n8ptl)
        if(idr.ne.0)then
          idptl(nptl)=idr
        else
          idptl(nptl)=id
        endif
        do j=1,4
          xorptl(j,nptl)=xori(j)
          xorptl(j,n8ptl)=xori(j)
        enddo
        tivptl(1,nptl)=xorptl(4,nptl)
        call idtau(idptl(nptl),pptl(4,nptl),pptl(5,nptl),taugm)
        tivptl(2,nptl)=tivptl(1,nptl)+taugm*(-alog(rangen()))
        if(abs(p1(5)-am).gt.0.01.and.ish.ge.1)then
          write (ifmt,'(2a,i6,2f8.2)') 'WARNING light string,',
     .    ' ptl off-shell, id mass minmass:',idr,p1(5),am
        endif
        goto 98
      endif
      endif

c................search breakpoints...................................

      n=0
      nsbp=1
      nbr=0
      iok=0

 11   continue
      if(nsbp.gt.10000)then
        if(ish.ge.1)then
        write(*,*)'Gakfra : string can not be fragmented, redo event !'
        write(*,*)nk1,idptl(nk1),nk2,idptl(nk2),istptl(nk1)
        write(ifch,*)'Gakfra : redo event !'
        endif
        iret=1
        goto 9999
      endif

      !----------------------!
      call gaksbp(0,1,iok)
      !----------------------!
      nbrtry=0
      goto 10

c..............delete breakpoint....................................
 9    if(ish.ge.4)write(ifch,*) 'delete breakpoint',n,' -> '
     &,nbr-1,' breakpoints left'
      call gakdel(n)              !hier loeschen

c..............no more breakpoints -> redo..........................
      if(nbr.eq.0)then
        nsbp=nsbp+1
        if(ish.ge.3)write (ifch,*)
     &    'no breakpoints left ... try again (',nsbp,')'
        if(ish.ge.3)write (ifch,*)' '
        goto 11
      endif

c............make index list of breakpoints still to be treated.........
c             (string-segment (ip=0) left or right of the BP)
 10   continue
      nbrtry=nbrtry+1
      nind=0
      do i=1,nbr
        if(ip(i).eq.0.or.ip(i+1).eq.0)then
          nind=nind+1
          ind(nind)=i
        endif
      enddo
      if(nbrtry.gt.1000000)then
        if(ish.ge.1)then
        write(*,*)'Gakfra : string can not be fragmented, redo event !'
        write(*,*)nk1,nk2,nbr,nind,pbk !bg
        endif
        iret=1
        goto 9999
      endif

c.....no more breakpoints to be treated...............................
      if(nind.eq.0) goto 100

c.....................................................................
      if(ish.ge.5)then
        write(ifch,*) 'breakpoints:'
        write(ifch,'(i3,i5)') 0,-iflb(0)
        do i=1,nbr
          write(ifch,'(i3,2i5,1x,4e12.4,2x,2i3,2f6.3)')
     &      i,iflb(i),-iflb(i),(ptb(j,i),j=1,4)
     &         ,ijb(1,i),ijb(2,i),xyb(1,i),xyb(2,i)
        enddo
        write(ifch,'(i3,i5)') nbr+1,iflb(nbr+1)
      endif

c.....choose breakpoint (index n), calculate masses, lambda.............
      r=rangen()
      nn=1+int(r*float(nind))
c      nn=1+(nind-1)*int(r+0.5)
      n=ind(nn)
      if(ish.ge.4)write(ifch,*) 'choose breakpoint',n
      nl=n-1
      nr=n+1
      do while (nl.gt.0.and.ip(nl).eq.0)
         nl=nl-1
      enddo
      do while (nr.lt.nbr+1.and.ip(nr+1).eq.0)
         nr=nr+1
      enddo
      if(ish.ge.6)write (ifch,*) 'nl,n,nr:',nl,n,nr
c      print *,'------------------------------------',1
      call gaksco(n-1,n,n+1,ijb(1,n),ijb(2,n),x1,x2,y1,y2)
      if(x2.le.x1.or.y2.le.y1)goto 9
cc      x=(xyb(1,n)-x1)/(x2-x1)
cc      y=(xyb(2,n)-y1)/(y2-y1)
      x=xyb(1,n)
      y=xyb(2,n)
      aml2=co(1)+x*co(2)+y*co(3)+x*y*co(4)
      amr2=co(5)+x*co(6)+y*co(7)+x*y*co(8)
      ala2=co(9)+x*co(10)+y*co(11)+x*y*co(12)

c.....determine id of both string-segments left and right of BP...............
      aml=sqrt(max(0.,aml2))
      idl=idsp( -iflb(n-1) , iflb(n)   )
      amr=sqrt(max(0.,amr2))
      idr=idsp( -iflb(n)   , iflb(n+1) )
      !icccl=0
      !icccr=0
      !if(abs(iflb(n-1)).eq.4.or.abs(iflb(n)).eq.4)icccl=1
      !if(abs(iflb(n)).eq.4.or.abs(iflb(n+1)).eq.4)icccr=1
      !if(icccl.eq.1.and.aml.lt.3.or.icccr.eq.1.and.amr.lt.3)then
      ! write(*,'(a,2i6,2f5.2,10i5)')
      !.' segments',idl,idr,aml,amr,n,nbr,(ip(nnn),nnn=1,nbr+1)
      !    goto 9 !redo string decay
      !endif

c.....if mass to small (because of spacelike pt)  reject..........
c      amin=0.0
c      if (aml2.le.amin.or.amr2.le.amin) goto 9

c.....if no left or right ptl id -> reject.........
      if(idl.eq.0.or.idr.eq.0)then
        if(ish.ge.5)write(ifch,*)'no left or right ptl id'
        goto 9
      endif

      iadjl=0
      iadjr=0
      idrl=0
      idrr=0

c.....if string-segment on left side, check for resonance...........
      rud=pmqu/(pmqu+pmqd)
      if(ip(n)  .eq.0) then
        aml=sqrt(max(0.,aml2+0.*min(0.,amr2))) !!!????
        amlxx=aml
        call idress(idl,aml,idrl,iadjl)
        r=rangen()
        if(idrl.eq.110.and.r.lt.rud)then
c          if (r.gt.eta+etap) goto 9      !rare numerical errors
          idl=220
c          aml=.549
c          if(r.lt.etap)aml=.958
c          amlxx=aml
          call idress(idl,aml,idrl,iadjl)
          iadjl=1
        endif
        if(ish.ge.5)write(ifch,'(a,i5,2f12.6,1i10,i5)')
     &    ' l:  ',idl,amlxx,aml,idrl,iadjl
      else
        if(ish.ge.5)write(ifch,'(a,i5,2f12.6,i10,a5)')
     &    ' l:  ',idl,aml,aml,ip(n),'ok'
      endif

c.....if string-segment on right side, check for resonance...........
      if(ip(n+1).eq.0) then
        amr=sqrt(max(0.,amr2+0.*min(0.,aml2))) !!!????
        amrxx=amr
        call idress(idr,amr,idrr,iadjr)
        r=rangen()
        if(idrr.eq.110.and.r.lt.rud)then
c          if (r.gt.eta+etap) goto 9    !rare numerical errors
          idr=220
c          amr=.549
c          if(r.lt.etap)amr=.958
c          amrxx=amr
          call idress(idr,amr,idrr,iadjr)
          iadjr=1
        endif
        if(ish.ge.5)write(ifch,'(a,i5,2f12.6,1i10,i5)')
     &    ' r:  ',idr,amrxx,amr,idrr,iadjr
      else
        if(ish.ge.5)write(ifch,'(a,i5,2f12.6,i10,a5)')
     &    ' r:  ',idr,amr,amr,ip(n+1),'ok'
      endif

c.....check segments left, right of BP.....do adjustments if needed.....
      iok=0
      if(ip(n+1).ne.0)then  !.........particle on right side
        if(idrl.eq.0)then
          call gaksbp(n-1,n,iok)
          if(iok.eq.1)goto 9
          goto 10
        endif
        if(iadjl.eq.1)then
           if(ish.ge.5)write(ifch,*)'mass adjustment 1'
           n2=n+1
           call gakanp(n-1,n,n2,aml**2,amr**2,0.,iok)
        endif
        if(iok.eq.1)goto 9
        ip(n)=idrl
      elseif(ip(n).ne.0)then !.........particle on left side
        if(idrr.eq.0)then
          call gaksbp(n,n+1,iok)
          if(iok.eq.1)goto 9
          goto 10
        endif
        if(iadjr.eq.1)then
           if(ish.ge.5)write(ifch,*)'mass adjustment 2'
           n2=n+1
           call gakanp(n-1,n,n2,aml**2,amr**2,0.,iok)
        endif
        if(iok.eq.1)goto 9
        ip(n+1)=idrr
      else       !.........adjusted mass neither left nor right
        if(idrr.eq.0.and.idrl.eq.0)then
          call gaksbp(n,n+1,iok)
          if(iok.eq.1)goto 9
          call gaksbp(n-1,n,iok)
          if(iok.eq.1)goto 9
          goto 10
        elseif(idrl.ne.0.and.idrr.eq.0)then
          if(iadjl.eq.1) then
           if(ish.ge.5)write(ifch,*)'mass adjustment 3'
           call gakanp(n-1,n,nr,aml**2,0.,ala2,iok)
          endif
        elseif(idrl.eq.0.and.idrr.ne.0)then
          if(iadjr.eq.1) then
           if(ish.ge.5)write(ifch,*)'mass adjustment 4'
           n2=n+1
           call gakanp(nl,n,n2,0.,amr**2,ala2,iok)
          endif
        else  !if(idrl.ne.0.and.idrr.ne.0)then
          if(iadjl.eq.1.or.iadjr.eq.1) then
           if(ish.ge.5)write(ifch,*)'mass adjustment 5'
           n2=n+1
           call gakanp(n-1,n,n2,aml**2,amr**2,0.,iok)
          endif
        endif
        if(iok.eq.1)goto 9
        ip(n)=idrl
        ip(n+1)=idrr
      endif
      if(ish.ge.4)then
        write(ifch,*) 'left/right particles:'
     &         ,ip(n),aml,'  ',ip(n+1),amr
      endif
      !if(abs(iflb(n-1)).eq.4.or.abs(iflb(n)).eq.4 .or.
      !.abs(iflb(n+1)).eq.4 )print*,'ip n n+1',ip(n),ip(n+1),aml,amr
      goto 10

c................................................................
c                         end of decay string  nk1 - nk2
c................................................................



 100  continue  ! <=========================================



c.....final list...............................................


c good for lambda-xf in pp but problem with K in pi interaction
c      if(kdrm.eq.1)then
c        if(nbr.eq.1.and.(!abs(iflb(2)).gt.1000.or.
c     *     abs(iflb(1)).gt.2))then
cc     *     abs(iflb(1)).gt.2.or.abs(iflb(0)).gt.1000))then
c          n=1
c          goto 9
cc        elseif(nbr.eq.2
cc     *     .and.(abs(iflb(0)).gt.1000.and.abs(iflb(2)).gt.1000))then
cc          n=2
cc          goto 9
c        endif
c      endif


      if(ish.ge.4)then
        write(ifch,*) '   ************ OK **************'
        write(ifch,*) 'final breakpoints:'
        write(ifch,'(i3,i5)') 0,-iflb(0)
        do i=1,nbr
          write(ifch,'(i3,2i5,1x,4e12.4,2x,2i3,2f6.3)')
     &      i,iflb(i),-iflb(i),(ptb(j,i),j=1,4)
     &         ,ijb(1,i),ijb(2,i),xyb(1,i),xyb(2,i)
        enddo
        write(ifch,'(i3,i5)') nbr+1,iflb(nbr+1)
      endif

c.....write particles in pptl-list................................
      if(ish.ge.3)then
        write(ifch,'(a)')' ---------produced particles---------'
      endif

      if(nptlini+nbr+1.gt.mxptl)
     &       call utstop('gakfra (end): mxptl too small ...&')
      do i=1,nbr+1
         nptl=nptl+1
         if(i.lt.nbr+1)then     !particle = left side of breakpoints
           call gaksco(i-1,i,i+1,ijb(1,i),ijb(2,i),x1,x2,y1,y2)
c           taubrr=pst(4,nob+7)+x*pst(4,nob+8)+y*pst(4,nob+9)
           x=xyb(1,i)
           y=xyb(2,i)
           aml2=co(1)+x*co(2)+y*co(3)+x*y*co(4)
           amr2=co(5)+x*co(6)+y*co(7)+x*y*co(8)

           ala2=co(9)+x*co(10)+y*co(11)+x*y*co(12)
           aml=sign(sqrt(abs(aml2)),aml2)
           amr=sign(sqrt(abs(amr2)),amr2)
           if(aml.le.0.d0.or.amr.le.0.d0)then
             if(ish.ge.4)write(ifch,*)
     & 'Negative mass, fragmentation not OK -> redo ...'
             n=i
             nptl=nptlini+1
             goto 9
           endif
           qsqptl(nptl)=ala2
           pptl(5,nptl)=aml
           do j=1,4
             pptl(j,nptl)=pst(j,nob+1)-
     &            x*pst(j,nob+2)+y*pst(j,nob+3)
c             pptr(j,nptl)=ptb(j,i)
           enddo
         else                   !last particle = right side of last breakpoint
           pptl(5,nptl)=amr
           qsqptl(nptl)=0.
           do j=1,4
             pptl(j,nptl)=pst(j,nob+4)+
     &            x*pst(j,nob+5)-y*pst(j,nob+6)
c             pptr(j,nptl)=ptb(j,i) !should be zero
           enddo
c           taubrr=0.
         endif
         idptl(nptl)=ip(i)
         if(ish.ge.7)call blist('&',nptl,nptl)
         if(pptl(4,nptl).le.0.d0)then
           if(ish.ge.4)write(ifch,*)
     & 'Negative energy, fragmentation not OK -> redo ...'
           n=i
           nptl=nptlini+1
           goto 9
         endif

      enddo

c..... fast charm ...............................

ctp this is needed to get the proper shape for lambda_c and D xp distribution in e+e-
ctp  be applied only in case of heavy quarks because otherwise it is
ctp less good for light particles. (2019/09)
      
c      if(nbr.le.2.and.p1(5).gt.2.*qmass(4).and.nsde.lt.10000)then
c      if(kkc.eq.1.and.nbr.le.2
c      if(nbr.eq.1
c     .   .and.(abs(iflb(0)).gt.1.or.abs(iflb(2)).gt.1)
c     .   .and.nsde.lt.10000)then
cc        print *,nbr,iflb(0),iflb(1),iflb(2)
cc        call blist('&',nptlini,nptl)
c        nptl=nptlini
c        goto 3
c       endif
c      if(abs(iflb(0)).gt.2.or.abs(iflb(nbr+1)).gt.2)then
c      if((nbr.eq.1.or.(nbr.eq.2
c     .   .and.(abs(iflb(1)).gt.2.or.abs(iflb(2)).gt.2)))
c     .   .and.nsde.lt.10000)then
c      if(nbr.le.2)then
      if(pmior.gt.0.)then
      if(nsde.lt.10000.and.max(pptl(4,nptl)
     .  ,pptl(4,nptl-nbr))/(pmior*(1.-0.66*rangen())).gt.0.5)then
c        print *,nbr,iflb(0),iflb(1),iflb(2)
        if(ish.ge.8)call blist('reject forward part.&',nptlini,nptl)
        nptl=nptlini
        goto 3
      endif
      endif
c     endif
c      npin=nptl
c      if(pptl(4,nptl-nbr).gt.pptl(4,nptl))npin=nptl-nbr
c      if(nsde.lt.10000.and.pptl(4,npin)
c     .   /(pmior*(1.-0.3*rangen()**(1.2/pptl(5,npin)))).gt.0.5)then
cc        print *,nbr,iflb(0),iflb(1),iflb(2)
c        if(ish.ge.8)call blist('reject forward part.&',nptlini,nptl)
c        nptl=nptlini
c        goto 3
c      endif

c      do k=nptlini+1,nptl
c        id=idptl(k)
c        !call idmass(id,am)
c        !zzz=pptl(4,k)**2-pptl(1,k)**2-pptl(2,k)**2-pptl(3,k)**2
c        !print*,'fra tst',k,id,pptl(4,k),pptl(5,k)
c        !.      ,sqrt(zzz),am
c        call idflav(id,ifl1,ifl2,ifl3,jspin,idu)
c        if(abs(ifl1).eq.4.or.abs(ifl2).eq.4.or.abs(ifl3).eq.4)then
c          ppl=pptl(4,k)+abs(pptl(3,k))
c          pmi=pptl(4,k)-abs(pptl(3,k)) 
c          x=engy/2-abs(pptl(3,k))
c          w=exp(-x**2/5)
c          if(x.lt.0.5.or.rangen().lt.w)then
c            !print*,'fast charm',id,ppl,pmi,0.5*log(ppl/pmi)
c            nptl=nptlini
c            goto 3
c          endif
c        endif
c      enddo

c      if(krm.eq.1)then
cc for remnant randomize direction if high mass forward
c        plong=0.
c        emax=0.
c        do n=nptlini+1,nptl
c          plong=max(plong,abs(pptl(3,n)))
c          emax=max(emax,pptl(4,n))
c        enddo
c        rlong=emax/pptl(5,nptlini)      !rlong < 0.5
c        smas1=0.
c        do i=1,1
c          if(abs(iflb(i)).gt.1000)then
c            pmas=qmass(0)
c          else
c            pmas=qmass(abs(iflb(i)))
c          endif
c          smas1=smas1+pmas
c        enddo
c        smas2=0.
c        do i=2,nbr+1
c          if(abs(iflb(i)).gt.1000)then
c            pmas=qmass(0)
c          else
c            pmas=qmass(abs(iflb(i)))
c          endif
c          smas2=smas2+pmas
c        enddo
c        the=0.25*pi*exp(-smas2/smas1)
c        the=min(atan(strcut/plong),the+rlong*pi*rangen())!     *strcut  ! *pptl(5,nptlini+1)
c        phi=2.*pi*rangen()
c        if(ish.ge.7)write(ifch,*)' rotation :',the,phi,fact,plong
c        do n=nptlini+1,nptl
c          if(ish.ge.7)call alistc('&',n,n)
c          if(ish.ge.7)write(ifch,'(a)')' ----->'
c          call utrotb(the,phi,pptl(1,n),pptl(2,n),pptl(3,n))
c          if(ish.ge.7)call alistc('&',n,n)
c        enddo
c      endif

      !print*,'TEST nk1,nk2,xq2s',nk1,nk2,xq2s

      nptl=nptlini+1  ! nptlini+1 is the string entry (with ist=29)

      if(ish.ge.7)then
        write(ifch,'(a)')' ---------produced particles---------'
      endif

      do i=1,nbr+1
         nptl=nptl+1

         call utlob4(-1,p1(1),p1(2),p1(3),p1(4),p1(5)
     $        ,pptl(1,nptl),pptl(2,nptl),pptl(3,nptl),pptl(4,nptl))
c         call utlob4(-1,p1(1),p1(2),p1(3),p1(4),p1(5)
c     $        ,pptr(1,nptl),pptr(2,nptl),pptr(3,nptl),pptr(4,nptl))


c........Origin..................................................
         istptl(nptl)=0
         iorptl(nptl)=n8ptl
         jorptl(nptl)=0
         ifrptl(1,nptl)=0
         ifrptl(2,nptl)=0

         r=rangen()
         amgener=1.0  !generic mass = 1 GeV
         tauran=-taustr*alog(r)*pptl(4,nptl)/pptl(5,nptl)    !amgener 
c         if(jpsi.ge.1)tauran=0.
c         tauran=max(taubrl,taubrr) !take formation time from string-theory
         xorptl(1,nptl)=xori(1)+pptl(1,nptl)/pptl(4,nptl)*tauran
         xorptl(2,nptl)=xori(2)+pptl(2,nptl)/pptl(4,nptl)*tauran
         xorptl(3,nptl)=xori(3)+pptl(3,nptl)/pptl(4,nptl)*tauran
         xorptl(4,nptl)=xori(4)+pptl(4,nptl)/pptl(4,nptl)*tauran
         !write(ifch,*),'*****fra*****',nptl
         !.    ,pptl(4,nptl)/pptl(5,nptl),tauran
         !.    ,xorptl(1,nptl),xorptl(2,nptl),xorptl(3,nptl),xorptl(4,nptl)
         tivptl(1,nptl)=xorptl(4,nptl)
         call idtau(idptl(nptl),pptl(4,nptl),pptl(5,nptl),taugm)
         tivptl(2,nptl)=tivptl(1,nptl)+taugm*(-alog(rangen()))
         ityptl(nptl)=itp !ityptl(nk1)
         itsptl(nptl)=itsptl(nk1)
         qsqptl(nptl)=qsqptl(n8ptl)
         zpaptl(1,nptl)=zpaptl(1,n8ptl)
         zpaptl(2,nptl)=zpaptl(2,n8ptl)
         rinptl(nptl)=-9999
         if(ish.ge.3)call blist('&',nptl,nptl)
ctp060829          taubrl=taubrr
      enddo
      iorptl(n8ptl)=nk1
      jorptl(n8ptl)=nk2
      ityptl(n8ptl)=itp !ityptl(nk1)
      do i=nk1,nk2
         istptl(i)=21
         ifrptl(1,i)=n8ptl
         ifrptl(2,i)=0
      enddo
      istptl(n8ptl)=29          !code for fragmented string
      ifrptl(1,n8ptl)=n8ptl+1
      ifrptl(2,n8ptl)=nptl
      xorptl(1,n8ptl)=xori(1)
      xorptl(2,n8ptl)=xori(2)
      xorptl(3,n8ptl)=xori(3)
      xorptl(4,n8ptl)=xori(4)
      if(ish.ge.5)then
        write(ifch,*)'string momentum sum:'
        do kk=1,5
          pptl(kk,nptl+1)=0
          do ii=n8ptl+1,nptl
            pptl(kk,nptl+1)=pptl(kk,nptl+1)+pptl(kk,ii)
          enddo
        enddo
        call alist2('&',n8ptl,n8ptl,nptl+1,nptl+1)
      endif



c.....another string?.........................................
 98   nk1=nk2+1
      goto 2                    !next string

c.....end of fragmentation.....................................
 9999 continue
      call utprix('gakfra',ish,ishini,4)
      return
      end

c---------------------------------------------------------------------
      subroutine gaksbp(ibr1,ibr2,iok)
c---------------------------------------------------------------------
      ! search break points
      !-----------------------------------------
      ! nbr ... number of break points
      !
      !
      !
      !--------------------------------------------------------------
#include "aaa.h"

      parameter (mxnbr=800,mxnba=10000)
      double precision pst,pf
      common /gag/pst(4,0:mxnba),nob,ipst(0:mxnba)
     $     ,nbr,ijb(2,0:mxnbr),xyb(2,0:mxnbr)
     &     ,ptb(4,0:mxnbr),iflb(0:mxnbr),ip(0:mxnbr),co,x,y
      common/pbk/pbkl(3) /cnsbp/nsbp /cn8ptl/n8ptl !bg
      double precision ax,ay,az,ae,am,bx,by,bz,be
      dimension co(12)
      logical weiter
      common/czz/kky,krm,kdrm,kzrm,pui(3),pdi(3),psi(3),pci(3),pbi(3) !bg
     &,zzzz,itp,pduui(3),pdudi(3),pdusi(3),pduci(3),pdubi(3),pdddi(3) !bg
     &,pddsi(3),pddci(3),pddbi(3),pdssi(3),pdsci(3),pdsbi(3),pdcci(3) !bg
     &,pdcbi(3),pdbbi(3),fkap(3),kkc !bg
      double precision psg
      common/zpsg/psg(5)
c      dimension pleft(5),pright(5)

      pf(i,j)=pst(4,i)*pst(4,j)-pst(3,i)*pst(3,j)
     &     -pst(2,i)*pst(2,j)-pst(1,i)*pst(1,j)
      mmod(i)=mod(mod(i,nob)+nob,nob)

      call utpri('gaksbp',ish,ishini,5)



      ib1=ibr1
      ib2=ibr2
c guess the initial band to set ig0 for pbk
      if(ib1.eq.0)then 
         jj=nob-1
      else
         jj=ijb(2,ib1)
      endif
      if(ib1.eq.0)then 
        ii=mmod(jj+1)
        if(jj.lt.nob/2)ii=mmod(nob-jj+1)
      else
        ii=ijb(1,ib1)
      endif
      ig0=1                                !quarks by default
      if(idptl(ipst(ii)).eq.-9)ig0=ig0+1    !gluon
      if(idptl(ipst(jj)).eq.-9)ig0=ig0+1    !gluon
      pbk=pbkl(ig0)  !bg
c      x=0.5*(qsqptl(ipst(ii))+qsqptl(ipst(jj)))
c      pbk=abs(pbreak)+(1.-exp(-x/2000.))*(pbreakg-abs(pbreak))
c continue
      iside=1
      if(rangen().lt.0.5)iside=2
      if(ish.ge.6)write(ifch,*)'iside:',iside
      Amxx=80./pbk !bg
      if(ish.ge.4)write(ifch,*)
     &'search brk between ib1=',ib1,' and ib2=',ib2
      ntry=0
      nbrold=nbr
      Amin=0.
      Amax=Amxx
                    !      <======================================
 26   ntry=ntry+1   !  <==========================================
                    !      <======================================
      ig=1
      if(jorptl(ipst(ii)).eq.-9)ig=ig+1
      if(jorptl(ipst(jj)).eq.-9)ig=ig+1
      pbk=pbkl(ig)  

      A0=-log(exp(-pbk*Amin)+rangen()                               !bg
     &  *(exp(-pbk*Amax)-exp(-pbk*Amin)))/pbk                       !bg
      if(ish.ge.6)write(ifch,*)'pbk, Amin, Amax, A0:',ig0,pbk, Amin !bg
     &                                              , Amax, A0
      A=A0
      Amaxn=0.
      if(iside.eq.2)goto 51
c.....................................................................
      if(ib1.eq.0)then          !startwert der j-schleife
         jj=nob-1
      else
         jj=ijb(2,ib1)
      endif
      do while (jj.gt.0)
         if(jj.eq.ijb(2,ib1) )then       !linker brk im Streifen jj?
            y1=xyb(2,ib1)
         else                   !nein
            y1=0.
         endif
         if(jj.eq.ijb(2,ib2))then       !rechter brk im Streifen jj?
            y2=xyb(2,ib2)
         else                   !nein
            y2=1.
         endif
         if(y1.eq.y2) goto 9999
         if(abs(y1-y2)/abs(y1+y2).le.1e-7) goto 9999
         if(ish.ge.6)write(ifch,*)'y1,y2,A',y1,y2,A
                                !startwert der i-schleife
         if(ib1.eq.0)then       !ohne linken brkpt
            ii=mmod(jj+1)
            if(jj.lt.nob/2)ii=mmod(nob-jj+1)
            if(ib2.le.nbr.and.mmod(ii+nob/2)
     &        .gt.mmod(ijb(1,ib2)+nob/2))then
              if(ish.ge.6)write(ifch,*) 'very special case',ii,jj
              goto12
            endif
         else
            ii=ijb(1,ib1)
            if(jj.lt.nob/2 .and.
     $           mmod(nob-jj+1+nob/2).gt. mmod(ii+nob/2)
     $           )ii=mmod(nob-jj+1)
         endif
         weiter=.true.
         aa=0.
         if(ii.eq.jj) then
            if(ish.ge.6)write(ifch,*) 'Rand erreicht'
            goto 15
         endif
         do while(weiter)
            if(ii.eq.ijb(1,ib1))then    !linker brk im Feld (ii,jj)
               x2=xyb(1,ib1)
            else
               x2=1.
            endif

            if(ii.eq.ijb(1,ib2))then    !rechter brk im Feld (ii,jj)
               x1=xyb(1,ib2)
            else
               x1=0.
            endif
            if(x1.eq.x2) goto 9999
            if(abs(x1-x2)/abs(x1+x2).le.1e-7) goto 9999
            f=1.0
            if(ipst(ii).ne.ipst(jj))aa=aa+2.*(x2-x1)*(y2-y1)*f*pf(ii,jj)
            if(ish.ge.6)write(ifch,*)'ii,jj,x1,x2,aa:',ii,jj
     $           ,ipst(ii).ne.ipst(jj),x1,x2,aa,pf(ii,jj)
            if(ii.eq.ijb(1,ib2))then
               weiter=.false.
            else
               ii=mmod(ii+1)
               if(ii.eq.jj.or.mmod(ii+jj).eq.0)weiter=.false.
            endif
         enddo
         Amaxn=Amaxn+aa
         if(aa.gt.A)goto 10
         A=A-aa
         if(jj.eq.ijb(2,ib2)) then
            if(ish.ge.6)write(ifch,*) 'brk erreicht'
            goto 15
         endif
 12      jj=mmod(jj-1)
      enddo
      goto 15

 10   continue
      yb=A/aa*(y2-y1)+y1
      if(ish.ge.6)write(ifch,*)'jj,yb ok:',jj,yb
      r=rangen()
      ra=aa*r
      if(ish.ge.6)write(ifch,*)'r,ra,aa',r,ra,aa
                                !nochmal die letzte ii-Schleife
      if(ib1.eq.0)then          !ohne linken brkpt
         ii=mmod(jj+1)
         if(jj.lt.nob/2)ii=mmod(nob-jj+1)
      else
        ii=ijb(1,ib1)
        if(jj.lt.nob/2 .and.
     $       mmod(nob-jj+1+nob/2).gt. mmod(ii+nob/2)
     $       )ii=mmod(nob-jj+1)
      endif
      do while(ra.gt.0.)
         if(ii.eq.ijb(1,ib1))then       !linker brk im Feld (ii,jj)
            x2=xyb(1,ib1)
         else
            x2=1.
         endif

         if(ii.eq.ijb(1,ib2))then       !rechter brk im Feld (ii,jj)
            x1=xyb(1,ib2)
         else
            x1=0.
         endif
         f=1.0
         ab=0.
         if(ipst(ii).ne.ipst(jj)) ab=2.*(x2-x1)*(y2-y1)*f*pf(ii,jj)
         if(ish.ge.6)write(ifch,*)'ii,jj,x1,x2,ab,ra:',ii,jj
     $        ,ipst(ii).ne.ipst(jj),x1,x2,ab,ra
         if(ab.gt.ra)then
            xb=ra/ab*(x2-x1)+x1
         else
            ii=mmod(ii+1)
         endif
         ra=ra-ab
      enddo
      if(ish.ge.5)write(ifch,*) 'breakpoint in field ',ii,jj
     & ,' at position ',xb,yb
      goto 95

c......................................................................
      !von rechts
 51   if(ib2.eq.nbr+1)then      !startwert der i-schleife
         ii=nob/2-1
      else
         ii=ijb(1,ib2)
      endif
      do while (ii.ne.nob/2)
         if(ii.eq.ijb(1,ib1))then !linker brk im Streifen (ii)
            x2=xyb(1,ib1)
         else
            x2=1.
         endif
         if(ii.eq.ijb(1,ib2))then !rechter brk im Streifen (ii)
            x1=xyb(1,ib2)
         else
            x1=0.
         endif
         if(x1.eq.x2) goto 9999
         if(abs(x1-x2)/abs(x1+x2).le.1e-7) goto 9999
         if(ish.ge.6)write(ifch,*)'x1,x2 A',x1,x2,A

         if(ib2.eq.nbr+1)then
            jj=mmod(ii+1)
            if(ii.gt.nob/2)jj=mmod(nob-ii+1)
            if(ib1.ge.1.and.jj.gt.ijb(2,ib1))then
               if(ish.ge.6)write(ifch,*) 'very special case',ii,jj
               goto 13
            endif
         else
            jj=ijb(2,ib2)
            if(ii.gt.nob/2 .and. mmod(nob-ii+1).gt.jj)jj=mmod(nob-ii+1)
         endif
         weiter=.true.
         aa=0.
         if(ii.eq.jj) then
            if(ish.ge.6)write(ifch,*) 'Rand erreicht'
            goto 15
         endif
         do while(weiter)
            if(jj.eq.ijb(2,ib1))then    !linker brk im Feld (ii,jj)
               y1=xyb(2,ib1)
            else
               y1=0.
            endif
            if(jj.eq.ijb(2,ib2))then !rechter brk im Feld (ii,jj)
               y2=xyb(2,ib2)
            else
               y2=1.
            endif
            if(y1.eq.y2) goto 9999
            if(abs(y1-y2)/abs(y1+y2).le.1e-7) goto 9999
            f=1.0
            if(ipst(ii).ne.ipst(jj))aa=aa+2.*(x2-x1)*(y2-y1)*f*pf(ii,jj)
            if(ish.ge.6)write(ifch,*)'ii,jj,x1,x2,aa:',ii,jj
     $           ,ipst(ii).ne.ipst(jj),x1,x2,aa
     &           ,pf(ii,jj)
            if(jj.eq.ijb(2,ib1))then
               weiter=.false.
            else
               jj=mmod(jj+1)
               if(jj.eq.ii.or.mmod(ii+jj).eq.0)weiter=.false.
            endif
         enddo
         Amaxn=Amaxn+aa
         if(aa.gt.A)goto 14
         A=A-aa
         if(ii.eq.ijb(1,ib1)) then
            if(ish.ge.6)write(ifch,*) 'brk erreicht'
            goto 15
         endif
 13      ii=mmod(ii-1)
      enddo
      goto 15



 14   continue
      xb=A/aa*(x2-x1)+x1
      if(ish.ge.6)write(ifch,*)'ii,xb ok:',ii,xb
      r=rangen()
      ra=aa*r
      if(ish.ge.6)write(ifch,*)'r,ra,aa',r,ra,aa
                                !nochmal die letzte jj-Schleife
      if(ib2.eq.nbr+1)then
         jj=mmod(ii+1)
         if(ii.gt.nob/2)jj=mmod(nob-ii+1)
      else
         jj=ijb(2,ib2)
         if(ii.gt.nob/2 .and. mmod(nob-ii+1).gt.jj)jj=mmod(nob-ii+1)
      endif

      do while(ra.gt.0.)
         if(jj.eq.ijb(2,ib1))then       !linker brk im Feld (ii,jj)
            y1=xyb(2,ib1)
         else
            y1=0.
         endif

         if(jj.eq.ijb(2,ib2))then       !rechter brk im Feld (ii,jj)
            y2=xyb(2,ib2)
         else
            y2=1.
         endif
         f=1.0
         ab=0.
         if(ipst(ii).ne.ipst(jj)) ab=2.*(x2-x1)*(y2-y1)*f*pf(ii,jj)
         if(ish.ge.6)write(ifch,*)'ii,jj,x1,x2,ab,ra:',ii,jj
     $        ,ipst(ii).ne.ipst(jj),x1,x2,ab,ra
         if(ab.gt.ra)then
            yb=ra/ab*(y2-y1)+y1
         else
            jj=mmod(jj+1)
         endif
         ra=ra-ab
      enddo
      if(ish.ge.5)write(ifch,*) 'breakpoint in field ',ii,jj
     & ,' at position ',xb,yb

 95   continue

c.....breakpoint accepted......................
      nbr=nbr+1
      if(nbr.gt.mxnbr-2) stop 'gaksbp: increase parameter mxnbr'
      do i=nbr+1,ib1+1,-1
         do j=1,2
            ijb(j,i)=ijb(j,i-1)
            xyb(j,i)=xyb(j,i-1)
         enddo
         do k=1,4
            ptb(k,i)=ptb(k,i-1)
         enddo
         iflb(i)=iflb(i-1)
         ip(i)=ip(i-1)
      enddo
      ip(ib1+1)=0
      ip(ib1+2)=0
      ijb(1,ib1+1)=ii
      ijb(2,ib1+1)=jj
      xyb(1,ib1+1)=xb
      xyb(2,ib1+1)=yb
      ig0=1                                !quarks by default
      if(idptl(ipst(ii)).eq.-9)ig0=ig0+1    !gluon
      if(idptl(ipst(jj)).eq.-9)ig0=ig0+1    !gluon

cc ..........diquark...............................................
c
c
c      pdiqu=pdiqua
c      if(kzrm.eq.1.or.krm.eq.0)then
cc        pdiqu=exp(log(pdiqua)/sqrt(1.+fkainc*(zzzz-1.)))
cc        pdiqu=exp(log(pdiqua)/(1.+fkainc*(zzzz-1.)))
c      endif
cc      if(iappl.eq.1.and.krm.eq.0)then
cc        pdiqu=pdiquak
cc      endif
c
c
c      if(nbr.le.2.and.abs(psg(3)/psg(4)).gt.diqcut)pdiqu=0.
c
c
c
c      if(rangen().lt.pdiqu.and.iabs(iflb(ib1)).lt.6
c     &  .and.iabs(iflb(ib2)).lt.6)then
c        jqu=2
c      else
c        jqu=1
c      endif
c

c ..........flavor...............................................

      ig=1
      if(jorptl(ipst(ii)).eq.-9)ig=ig+1
      if(jorptl(ipst(jj)).eq.-9)ig=ig+1
      !-----------------------------------------------------------------
      !subroutine timsho
      ! jorprt(nprtj) ...... special flavor information 
      !                     -9 = high Q2 gluon
      !                     gives value to all daughters except for 
      !                          high Q2 gloun-> quark
      !Later written to jorptl
      !-----------------------------------------------------------------
c      if(abs(iflb(ib1)).ge.4.or.abs(iflb(ib1+2)).ge.4)then !charm-X string
ccc      if(abs(jorptl(ipst(ii))).eq.4.or.abs(jorptl(ipst(jj))).eq.4)then !charm-X string
ccc      if(kkc.eq.1)then !HQ
ccc        print*,jorptl(ipst(ii)),jorptl(ipst(jj)),ib1+1,nbr
ccc     .        ,iflb(ib1),iflb(ib1+1),iflb(ib1+2)
c        fqq=exp(pi*pmqqx**2/fkap(ig)) !10
c      else
c        fqq=1
c      endif 

      fqq=1
      fqqq=1
      if(abs(iflb(ib1)).ge.4.or.abs(iflb(ib1+2)).ge.4)then!HQ-X string
        fqq=3       
        fqqq=3      
      endif
 
      !quark breakup
      pu=pui(ig)
      pd=pdi(ig)
      ps=psi(ig)
      pc=pci(ig)
      pb=pbi(ig) 
      !diquark breakup
      pduu=pduui(ig) *fqq
      pdud=pdudi(ig) *fqq
      pdus=pdusi(ig) *fqq *fqqq
      pduc=pduci(ig) *fqq *fqqq
      pdub=pdubi(ig) *fqq *fqqq
      pddd=pdddi(ig) *fqq
      pdds=pddsi(ig) *fqq *fqqq
      pddc=pddci(ig) *fqq *fqqq
      pddb=pddbi(ig) *fqq *fqqq
      pdss=pdssi(ig) *fqq *fqqq
      pdsc=pdsci(ig) *fqq *fqqq
      pdsb=pdsbi(ig) *fqq *fqqq
      pdcc=pdcci(ig) *fqq *fqqq
      pdcb=pdcbi(ig) *fqq *fqqq
      pdbb=pdbbi(ig) *fqq *fqqq


!      print *,pu,pd,ps,pc,difud,difus,difuc
c      print *,krm,kdrm,nbr,abs(psg(3)/psg(4))

c suppress forward (proj) or backward (targ) diquark production
c in non central strings (baryon spectra in pion collisions or NA49 data)
c ----> this should not be used for meson projectile because of kaon spectra
c       at 100 GeV (Barton, etc) <--------
c      if(kdrm.eq.1.and.ib1.le.1.and.abs(psg(3)/psg(4)).gt.strcut)then
c        ps=0.
c        pc=0.
c      endif

c suppress forward (proj) or backward (targ) diquark production
c in non central strings (baryon spectra in pion collisions or NA49 data)
c      if(kdrm.eq.1.and.ib1.le.1.and.abs(psg(3)/psg(4)).gt.diqcut)then
c     if(ib1.le.1.and.abs(psg(3)/psg(4)).gt.diqcut)then
c      if(diqcut.gt.0.)then
c        xy0=1.
c        if(ib1.le.1)xy0=xy0*xyb(1,ib1)*xyb(2,ib1)
c        if(ib2.ge.nbr-1)xy0=xy0*xyb(1,ib2)*xyb(2,ib2)
cc        write(ifch,*)'xy0',ib1,nbr,ib2,xy0
c        if((ib1.le.1.or.ib2.ge.nbr-1).and.
c     & xy0.le.0.)then
cc        if((ib1.le.1.or.ib2.ge.nbr-1).and.
cc     & rangen().lt.diqcut)then
cc     & abs(psg(3))/sqrt(psg(1)**2+psg(2)**2+psg(3)**2).gt.1.-diqcut)then
c          pduu=0.
c          pdud=0.
c          pdus=0.
c          pduc=0.
c          pdub=0.               !bg
c          pddd=0.
c          pdds=0.
c          pddc=0.
c          pddb=0.               !bg
c          pdss=0.
c          pdsc=0.
c          pdsb=0.               !bg
c          pdcc=0.
c          pdcb=0.               !bg
c          pdbb=0.               !bg
c        endif
c      endif

c suppress forward (proj) or backward (targ) strangeness (or heavier) production
c in non central strings (lambda spectra in pion collisions or NA49 data)
c      if(abs(psg(3)).gt.1e-5.and.strcut.gt.0.)then
c        if((ib1.le.1.or.ib2.ge.nbr-1).and.
c     & abs(psg(3))/sqrt(psg(1)**2+psg(2)**2+psg(3)**2).gt.1.-strcut)then
c          ps=0.
c          pc=0.
c          pb=0.
c          pdus=0.
c          pduc=0.
c          pdub=0.               !bg
c          pdds=0.
c          pddc=0.
c          pddb=0.               !bg
c          pdss=0.
c          pdsc=0.
c          pdsb=0.               !bg
c          pdcc=0.
c          pdcb=0.               !bg
c          pdbb=0.               !bg
c        endif
c      endif



      pdiqu=pduu+pdud+pdus+pduc+pdub+pddd+pdds+pddc+pddb+pdss+pdsc+pdsb
     * +pdcc+pdcb+pdbb !bg
      psum=pu+pd+ps+pc+pb+pdiqu !bg
c      print *,krm,kdrm,kky,zzzz,ps/(pu+pd+ps),pdiqu/psum

      r=rangen()*psum
      rdiqu=r-psum+pdiqu
      jqu=1
      ifl1=0
      ifl2=0
      if(rdiqu.gt.0.)then
        jqu=2
        r=rdiqu
        if(r.gt.pdiqu-pduu.and.pduu.gt.0.)then              ! --- LF --- 
          ifl1=1
          ifl2=1
        elseif(r.gt.pdiqu-pduu-pdud.and.pdud.gt.0.)then
          if(rangen().gt.0.5)then
            ifl1=1
            ifl2=2
          else
            ifl1=2
            ifl2=1
          endif
        elseif(r.gt.pdiqu-pduu-pdud-pddd.and.pddd.gt.0.)then
          ifl1=2
          ifl2=2
        elseif(r.gt.pdiqu-pduu-pdud-pddd-pdus.and.pdus.gt.0.)then
          if(rangen().gt.0.5)then
            ifl1=1
            ifl2=3
          else
            ifl1=3
            ifl2=1
          endif
        elseif(r.gt.pdiqu-pduu-pdud-pddd-pdus-pdds.and.pdds.gt.0.)then
          if(rangen().gt.0.5)then
            ifl1=2
            ifl2=3
          else
            ifl1=3
            ifl2=2
          endif
        elseif(r.gt.pdcc+pdsc+pddc+pduc+pdbb+pdcb+pdsb+pddb+pdub
     .         .and.pdss.gt.0.)then !bg
          ifl1=3
          ifl2=3
        elseif(r.gt.pdcc+pdsc+pddc+pdbb+pdcb+pdsb+pddb+pdub ! --- HF ---  
     .          .and.pduc.gt.0.)then !bg
          if(rangen().gt.0.5)then
            ifl1=1
            ifl2=4
          else
            ifl1=4
            ifl2=1
          endif
        elseif(r.gt.pdcc+pdsc+pdbb+pdcb+pdsb+pddb+pdub
     .         .and.pddc.gt.0.)then !bg
          if(rangen().gt.0.5)then
            ifl1=2
            ifl2=4
          else
            ifl1=4
            ifl2=2
          endif
        elseif(r.gt.pdcc+pdbb+pdcb+pdsb+pddb+pdub.and.pdsc.gt.0.)then !bg
          if(rangen().gt.0.5)then
            ifl1=3
            ifl2=4
          else
            ifl1=4
            ifl2=3
          endif
        elseif(r.gt.pdbb+pdcb+pdsb+pddb+pdub.and.pdcc.gt.0.)then !bg
          ifl1=4
          ifl2=4
        elseif(r.gt.pdbb+pdsb+pddb+pdub.and.pdcb.gt.0.)then !bg
          if(rangen().gt.0.5)then
            ifl1=4
            ifl2=5
          else
            ifl1=5
            ifl2=4
          endif
        elseif(r.gt.pdbb+pdsb+pddb.and.pdub.gt.0.)then !bg
          if(rangen().gt.0.5)then
            ifl1=1
            ifl2=5
          else
            ifl1=5
            ifl2=1
          endif
        elseif(r.gt.pdbb+pdsb.and.pddb.gt.0.)then !bg
          if(rangen().gt.0.5)then
            ifl1=2
            ifl2=5
          else
            ifl1=5
            ifl2=2
          endif
        elseif(r.gt.pdbb.and.pdsb.gt.0.)then !bg
          if(rangen().gt.0.5)then
            ifl1=3
            ifl2=5
          else
            ifl1=5
            ifl2=3
          endif
        elseif(pdbb.gt.0.)then !bg
          ifl1=5
          ifl2=5
        else
          ifl1=min(3,1+int(rangen()*3))
          ifl2=min(3,1+int(rangen()*3))
          write(ifmt,'(a,2i3,12e9.2)')'WARNING gakfra diquark flavor'
     .    ,ifl1,ifl2,r,pdiqua
     .    ,pdcc+pdsc+pddc+pduc+pdbb+pdcb+pdsb+pddb+pdub
     .    ,pdcc,pdsc,pddc,pduc,pdbb,pdcb,pdsb,pddb,pdub
        endif
      elseif(r.gt.pc+pd+ps+pb.and.pu.gt.0.)then !bg
        ifl1=1
      elseif(r.gt.pc+ps+pb.and.pd.gt.0.)then !bg
        ifl1=2
      elseif(r.gt.pc+pb.and.ps.gt.0.)then !bg
        ifl1=3
      elseif(r.gt.pb.and.pc.gt.0.)then !bg
        ifl1=4
      elseif(pb.gt.0.)then !bg
        ifl1=5
      else
        call utstop("Problem with flavor in gakfra !&")
      endif

      if(jqu.eq.2)then  !diquark ------
        ifl=-min(ifl1,ifl2)*1000-max(ifl1,ifl2)*100
      else              !quark ------
        ifl=ifl1
      endif

      if(iflb(ib1).lt.0.and.iflb(ib1).ge.-6)ifl=-ifl
      if(iflb(ib1).gt.1000)ifl=-ifl
      iflb(ib1+1)=ifl
      if(ish.ge.5)write(ifch,20) ig,ifl,pu,pd,ps,pc,pb,pdiqu/psum
     &,pduu,pdud,pdus,pduc,pdub,pddd,pdds,pddc,pddb,pdss,pdsc,pdsb,pdcc !bg
     &,pdcb,pdbb
 20   format('ig,flavor,pu,pd,ps,pc,pb,pdiqu:',i2,i6,6g13.5,/      !bg
     &,'  diquark u:',5g15.5,/,'  diquark d:',4g15.5,/
     &,'  diquark s:',3g15.5,10x,/,'diquark c:',2g15.5
     &,'  diquark b:',g15.5) !bg
c..............................pt.......................................
c      if(krm.ne.1.and.kdrm.ne.1)then
c        icub=ib1+1
c        !---------------------------------------------------------------------------
c        !  ib1+1 is the current break point index
c        !             (between 1 and nbr)
c        !  ijb(1,ib1) and ijb(2,ib1) are band indices
c        !         each index from 0 to nob-1 (nob= number of bands)
c        !         0 is left, then come the gluons, then antiquark, then agin gluons
c        !         like q - g - g - g - ~q - g - g - g
c        !--------------------------------------------------------------------------
c        call gaksco(icub-1,icub,icub+1
c     &    ,ijb(1,icub),ijb(2,icub),x1,x2,y1,y2)
c        x=xyb(1,icub)
c        y=xyb(2,icub)
c        aml2=co(1)+x*co(2)+y*co(3)+x*y*co(4)
c        amr2=co(5)+x*co(6)+y*co(7)+x*y*co(8)
c        aml=sign(sqrt(abs(aml2)),aml2)
c        amr=sign(sqrt(abs(amr2)),amr2)
c        !------segment left of current breakpoint icub -----
c        pleft(5)=aml
c        do j=1,4
c          pleft(j)=pst(j,nob+1)-x*pst(j,nob+2)+y*pst(j,nob+3)
c        enddo
c        call utlob4(-1,psg(1),psg(2),psg(3),psg(4),psg(5)
c     $     ,pleft(1),pleft(2),pleft(3),pleft(4))
c        !------segment right of current breakpoint icub-----
c        pright(5)=amr
c        do j=1,4
c          pright(j)=pst(j,nob+4)+x*pst(j,nob+5)-y*pst(j,nob+6)
c        enddo
c        call utlob4(-1,psg(1),psg(2),psg(3),psg(4),psg(5)
c     $     ,pright(1),pright(2),pright(3),pright(4))
c        !-------------------------
c        amt=pleft(5)**2+pleft(1)**2+pleft(2)**2
c        if(amt.gt.0..and.pleft(4)+abs(pleft(3)).gt.0.d0)then
c          amt=sqrt(amt)
c          yleft=sign(1.,pleft(3))*alog((pleft(4)+abs(pleft(3)))/amt)
c        else
c          yleft=0.                  !
c        endif
c        amt=pright(5)**2+pright(1)**2+pright(2)**2
c        if(amt.gt.0..and.pright(4)+abs(pright(3)).gt.0.d0)then
c          amt=sqrt(amt)
c          yright=sign(1.,pright(3))*alog((pright(4)+abs(pright(3)))/amt)
c        else
c          yright=0.                  !
c        endif
c        ybrk=(yleft+yright)/2.
c        yhax=2                  !0.5*yhaha
c        zzipx=1
c        if(ybrk.gt.yhax)then
c          zzipx=zzipx+(zzipp-1)
c        elseif(ybrk.gt.-yhax)then
c          zzipx=zzipx+(zzipp-1)*(ybrk+yhax)/(2*yhax)
c        endif
c        if(ybrk.lt.-yhax)then
c          zzipx=zzipx+(zzipt-1)
c        elseif(ybrk.lt.yhax)then
c          zzipx=zzipx+(zzipt-1)*(yhax-ybrk)/(2*yhax)
c        endif
c      endif
      delptfra=0.
      if(ifl2.ne.0)delptfra=delptfra+ptfrasr
c      if(ifl1.eq.3.and.ifl2.eq.0)delptfra=delptfra+ptfrasr
c      if(ifl1.eq.3.or.ifl2.eq.3)delptfra=delptfra+ptfrasr
      fnsbp=1
      if(nsbp.gt.9)fnsbp=0
c pt kink due to split elastic scattering (if zzzz.ne.1.)
c      if(kky.eq.1)then
      if(ig.eq.1)then       !quark and soft strings
c same pt for q-qb or qq-qbqb because of NA49 <pt> for ap + Barton
        pttra=(ptfraqq+delptfra)*(zzzz-1.)*fnsbp
        pt=ranpt()*(ptfraqq +  pttra)
      else
c same pt for q-qb or qq-qbqb because of pt at Tevatron ??? (less good at RHIC) 
        pttra=(ptfra+delptfra)*(zzzz-1.)*fnsbp
        pt=ranpt()*(ptfra   +  pttra)          !??????????????? tp100909
c        pt=ranpt()*(ptfra   +  pttra)/float(jqu)
      endif

      beta=2.*pi*rangen()

      if(ish.ge.5)then
        write(ifch,*) 'pt:',pt
      endif
      ptb(1,ib1+1)=pt*cos(beta)
      ptb(2,ib1+1)=pt*sin(beta)
      ptb(3,ib1+1)=0.
      ptb(4,ib1+1)=0.
      if(ish.ge.8)then
        write(ifch,*) 'the bands'
        write(ifch,'(4g12.6)') (pst(i,ii),i=1,4)
        write(ifch,'(4g12.6)') (pst(i,jj),i=1,4)
        write(ifch,*) 'pt before rotation'
        write(ifch,'(4f12.8)') (ptb(i,ib1+1),i=1,4)
      endif
      ax=pst(1,ii)+pst(1,jj)
      ay=pst(2,ii)+pst(2,jj)
      az=pst(3,ii)+pst(3,jj)
      ae=pst(4,ii)+pst(4,jj)
      am=sqrt(max(1d-10,(ae+az)*(ae-az)-ax**2-ay**2)) !?????????
      if(ish.ge.8)then
        write(ifch,*) 'boost vector ( region ) '
        write(ifch,'(5g12.6)') ax,ay,az,ae,am,pf(ii,jj)
      endif
      bx=pst(1,ii)
      by=pst(2,ii)
      bz=pst(3,ii)
      be=pst(4,ii)
      call utlob2(1,ax,ay,az,ae,am,bx,by,bz,be,60)
      if(ish.ge.8) then
        write (ifch,*) 'boost of b'
        write (ifch,*) 'bx,by,bz,be',bx,by,bz,be
      endif
      if(abs(bx)+abs(by)+abs(bz).gt.0.)then
        call utrot4(-1,bx,by,bz,ptb(1,ib1+1),ptb(2,ib1+1),ptb(3,ib1+1))
      else
        write(ifmt,*) 'null rot of pt',bx,by,bz
        write(ifmt,'(4f12.8)') (ptb(i,ib1+1),i=1,4)
      endif
      if(ish.ge.8) then
        write (ifch,*) 'rot of pt'
        write(ifch,'(4f12.8)') (ptb(i,ib1+1),i=1,4)
      endif

      call utlob4(-1,ax,ay,az,ae,am
     &  ,ptb(1,ib1+1),ptb(2,ib1+1),ptb(3,ib1+1),ptb(4,ib1+1))

      if(ish.ge.8) then
        write (ifch,*) 'backboost of pt'
        write(ifch,'(4f12.8)') (ptb(i,ib1+1),i=1,4)
      endif
c      if(az.eq.0..and.ay.eq.0.)az=1e-8    !not needed if following line commented
c      if(ish.ge.8)write(ifch,*)'rot vector:',ax,ay,az,ae,am

c.....call utrota(-1,sngl(ax),sngl(ay),sngl(az)  !already commented in nexus2
c....&              ,ptb(1,ib1+1),ptb(2,ib1+1),ptb(3,ib1+1))

      if(ish.ge.6)then
        write(ifch,*) 'pt'
        write(ifch,'(4g12.6)') (ptb(i,ib1+1),i=1,4)
        write (ifch,*) ijb(1,ib1+1),ijb(2,ib1+1),xyb(1,ib1+1)
     $       ,xyb(2,ib1+1)
      endif

      if(iside.eq.1)then
         ib1=ib1+1
         ib2=ib2+1
      endif


c      Amin=0.
c      if(Amax.lt.Amxx) goto 15
c      goto 25

 15   continue
      if(ish.ge.6)write(ifch,*) 'Amax:',Amax,Amaxn
      if (nbr.eq.nbrold) then
         Amax=Amaxn
         Amin=0.
         if(pbk*Amax.le.0..or.ntry.ge.1000)then  !bg
           goto 9999
         endif
         goto 26
      endif

      if(ish.ge.6)then
         write(ifch,*) 0,iflb(0)
         do i=1,nbr
           if(i.eq.ib2) write(ifch,*) '.................'
            write(ifch,'(i3,2x,2(i3),2(" ",g14.7),3x,i5,4(" ",g12.6))')
     &           i,ijb(1,i),ijb(2,i),xyb(1,i),xyb(2,i)
     &           ,iflb(i),(ptb(j,i),j=1,4)
           if(i.eq.ibr1) write(ifch,*) '.................'
         enddo
         write(ifch,*) nbr+1,iflb(nbr+1)
      endif

 9999 if(nbr.eq.nbrold)iok=1
      call utprix('gaksbp',ish,ishini,5)
      return
      end

c----------------------------------------------------------------------
      function ranptcut(xcut)
c----------------------------------------------------------------------
  12  xmx=50
      r=2.
      do while (r.gt.1.)
  11    x=sqrt(exp(rangen()*log(1+xmx**2))-1)
        if(x.eq.0.)goto 11
        r=rangen()  /  ( exp(-x)*(1+x**2) )
      enddo
      x=x/2.
      if(xcut.gt.0.)then
        if(rangen().lt.x/xcut)goto 12 !to avoid x>xcut, little change for x<<xcut
      endif
      ranptcut=x
      end

c----------------------------------------------------------------------
      function ranptdcut(xcut) !x*exp(-x**2)
c----------------------------------------------------------------------
 12   x=sqrt(-(log(rangen()))/(3.1415927/4.)) ! x * exp( -x**2/2 * (pi/2) )  which gives <x>=1
      if(xcut.gt.0.)then
        if(rangen().lt.x/xcut)goto 12 !to avoid x>xcut, little change for x<<xcut
      endif
      ranptdcut=x
      end

c----------------------------------------------------------------------
      function ranpticut(xcut) !x*exp(-sqrt(x))
c----------------------------------------------------------------------
 12   xmx=500
      r=2.
      do while (r.gt.1.)
        x=sqrt(exp(rangen()*log(1+xmx**2))-1)
        r=rangen()  /  ( exp(-sqrt(x))*(1+x**2)/5. )
      enddo
      x=x/20.
      if(rangen().gt.xcut/x)goto 12
      ranpticut=x
      end
      
c----------------------------------------------------------------------
      function ranptjcut(xcut) !x*exp(-x)
c----------------------------------------------------------------------
  12  continue
      x=-log(rangen())-log(rangen())
      x=x/2 !<x>=1
      if(x.gt.xcut)goto 12
      ranptjcut=x
      end

c----------------------------------------------------------------------
      function ranpt()
c----------------------------------------------------------------------

c .........exp(-x)
      xmx=50
      r=2.
      do while (r.gt.1.)
  11    x=sqrt(exp(rangen()*log(1+xmx**2))-1)
        if(x.eq.0.)goto11
        r=rangen()  /  ( exp(-x)*(1+x**2) )
      enddo
      ranpte=x/2.

c     -------------
      ranpt=ranpte
      return
c     -------------
cc .........exp(-x**2)
c      ranptg=sqrt(-log(rangen())/(3.1415927/4.)) !gauss
c
cc .........exp(-sqrt(x))
c      xmx=500
c      r=2.
c      do while (r.gt.1.)
c        x=sqrt(exp(rangen()*log(1+xmx**2))-1)
c        r=rangen()  /  ( exp(-sqrt(x))*(1+x**2)/5. )
c      enddo
c      ranpts=x/20.
c

      end

c----------------------------------------------------------------------
      function ranptd()
c----------------------------------------------------------------------

c .........exp(-x**2)
      ranptg=sqrt(-log(rangen())/(3.1415927/4.)) !gauss


c     -------------
      ranptd=ranptg
      return
c     -------------
c
cc .........exp(-x)
c      xmx=50
c      r=2.
c      do while (r.gt.1.)
c  11    x=sqrt(exp(rangen()*log(1+xmx**2))-1)
c        if(x.eq.0.)goto11
c        r=rangen()  /  ( exp(-x)*(1+x**2) )
c      enddo
c      ranpte=x/2.
c
cc .........exp(-sqrt(x))
c      xmx=500
c      r=2.
c      do while (r.gt.1.)
c        x=sqrt(exp(rangen()*log(1+xmx**2))-1)
c        r=rangen()  /  ( exp(-sqrt(x))*(1+x**2)/5. )
c      enddo
c      ranpts=x/20.





      end
c----------------------------------------------------------------------
      subroutine gakdel(ibr)
c----------------------------------------------------------------------
      parameter (mxnbr=800,mxnba=10000)
      double precision pst
      common /gag/pst(4,0:mxnba),nob,ipst(0:mxnba)
     $     ,nbr,ijb(2,0:mxnbr),xyb(2,0:mxnbr)
     &     ,ptb(4,0:mxnbr),iflb(0:mxnbr),ip(0:mxnbr),co,x,y
      dimension co(12)

      do i=ibr,nbr+1
         do j=1,2
            ijb(j,i)=ijb(j,i+1)
            xyb(j,i)=xyb(j,i+1)
         enddo
         do k=1,4
            ptb(k,i)=ptb(k,i+1)
         enddo
         iflb(i)=iflb(i+1)
         ip(i)=ip(i+1)
      enddo
      ip(ibr)=0
      nbr=nbr-1
      end

c----------------------------------------------------------------------
      subroutine gaksco(ibr1,ibr,ibr2,ib,jb,x1,x2,y1,y2)
c----------------------------------------------------------------------
#include "aaa.h"

      parameter (mxnbr=800,mxnba=10000)
      double precision pst,pf
      common /gag/pst(4,0:mxnba),nob,ipst(0:mxnba)
     $     ,nbr,ijb(2,0:mxnbr),xyb(2,0:mxnbr)
     &     ,ptb(4,0:mxnbr),iflb(0:mxnbr),ip(0:mxnbr),co,x,y
      dimension co(12)
      logical weiter

      pf(i,j)=pst(4,i)*pst(4,j)-pst(3,i)*pst(3,j)
     &     -pst(2,i)*pst(2,j)-pst(1,i)*pst(1,j)
      mmod(i)=mod(mod(i,nob)+nob,nob)

      call utpri('gaksco',ish,ishini,8)

      if(ish.ge.8)then
         write(ifch,*) 'zwischen brk:',ibr1,ibr,ibr2,'(',nob,')',nbr
         write(ifch,*) 'region:',ib,jb
      endif
      if(ib.eq.ijb(1,ibr1))then
         x2=xyb(1,ibr1)
      else
         x2=1.
      endif
      if(ib.eq.ijb(1,ibr2))then
         x1=xyb(1,ibr2)
      else
         x1=0.
      endif
      if(jb.eq.ijb(2,ibr1))then
         y1=xyb(2,ibr1)
      else
         y1=0.
      endif
      if(jb.eq.ijb(2,ibr2))then
         y2=xyb(2,ibr2)
      else
         y2=1.
      endif

c.....left side...................................................
      n=nob+1
      if(ish.ge.8)write(ifch,*)'x1,x2',x1,x2
      do i=1,4
cc        pst(i,n)=(x2-x1)*pst(i,ib)+ptb(i,ibr)-ptb(i,ibr1)
        pst(i,n)=     x2*pst(i,ib)+ptb(i,ibr)-ptb(i,ibr1)-y1*pst(i,jb)
      enddo
      if(ish.ge.8)write(ifch,*) 'add a  1 ii',1.,ib
      ii=mmod(ib-1)
      weiter=.true.
      if(ib.eq.ijb(1,ibr1))weiter=.false.
      do while(ii.ne.jb.and.weiter) !linker Rand??
         f1=1.
         if(ii.eq.ijb(1,ibr1))f1=xyb(1,ibr1)
         do i=1,4
            pst(i,n)=pst(i,n)+f1*pst(i,ii)
         enddo
         if(ish.ge.8)write(ifch,*) 'add a f1 ii',f1,ii
         if(ii.eq.ijb(1,ibr1))weiter=.false.
         ii=mmod(ii-1)
      enddo
      jj=mmod(jb+1)
      weiter=.not.weiter
      if(jb.eq.ijb(2,ibr1))weiter=.false.
      do while(weiter)
         f1=1.
         if(jj.eq.ijb(2,ibr1))f1=1.-xyb(2,ibr1)
         do i=1,4
            pst(i,n)=pst(i,n)+f1*pst(i,jj)
         enddo
         if(ish.ge.8)write(ifch,*) 'add b f1 ii',f1,jj
         if(jj.eq.ijb(2,ibr1))weiter=.false.
         jj=mmod(jj+1)
      enddo
      do i=1,4
cc        pst(i,n+1)=(x2-x1)*pst(i,ib)
        pst(i,n+1)=        pst(i,ib)
cc        pst(i,n+2)=(y2-y1)*pst(i,jb)
        pst(i,n+2)=        pst(i,jb)
      enddo
      co(1)= pf(n,n)
      co(2)=-2.*pf(n,n+1)
      co(3)= 2.*pf(n,n+2)
      co(4)=-2.*pf(n+1,n+2)
      if(ish.ge.8) then
        do i=n,n+2
          write (ifch,'(4g12.5)') (pst(j,i),j=1,4)
        enddo
      endif
      if(ish.ge.8)write(ifch,*) 'co left:',co(1),co(2),co(3),co(4)

c.....right side...................................................
      n=nob+4
      if(ish.ge.8)write(ifch,*)'y1,y2',y1,y2
      do i=1,4
cc         pst(i,n)=(y2-y1)*pst(i,jb)-ptb(i,ibr)+ptb(i,ibr2)
         pst(i,n)=    y2*pst(i,jb)-ptb(i,ibr)+ptb(i,ibr2)-x1*pst(i,ib)
      enddo
      if(ish.ge.8)write(ifch,*) 'add a  1 jj',1.,jb
      ii=mmod(ib+1)
      weiter=.true.
      if(ib.eq.ijb(1,ibr2))weiter=.false.
      do while(ii.ne.jb.and.weiter)
         f1=1.
         if(ii.eq.ijb(1,ibr2))f1=1.-xyb(1,ibr2)
         do i=1,4
            pst(i,n)=pst(i,n)+f1*pst(i,ii)
         enddo
         if(ish.ge.8)write(ifch,*) 'add a f1 ii',f1,ii
         if(ii.eq.ijb(1,ibr2))weiter=.false.
         ii=mmod(ii+1)
      enddo
      jj=mmod(jb-1)
      weiter=.not.weiter
      if(jb.eq.ijb(2,ibr2))weiter=.false.
      do while(weiter)
         f1=1.
         if(jj.eq.ijb(2,ibr2))f1=xyb(2,ibr2)
         do i=1,4
            pst(i,n)=pst(i,n)+f1*pst(i,jj)
         enddo
         if(ish.ge.8)write(ifch,*) 'add b f1 ii',f1,jj
         if(jj.eq.ijb(2,ibr2))weiter=.false.
         jj=mmod(jj-1)
      enddo
      do i=1,4
cc         pst(i,n+1)=(x2-x1)*pst(i,ib)
         pst(i,n+1)=        pst(i,ib)
cc         pst(i,n+2)=(y2-y1)*pst(i,jb)
         pst(i,n+2)=         pst(i,jb)
      enddo
      co(5)=pf(n,n)
      co(6)= 2.*pf(n,n+1)
      co(7)=-2.*pf(n,n+2)
      co(8)=-2.*pf(n+1,n+2)
      if(ish.ge.8) then
        do i=n,n+2
          write (ifch,'(4g12.5)') (pst(j,i),j=1,4)
        enddo
      endif
      if(ish.ge.8)write(ifch,*) 'co right:',co(5),co(6),co(7),co(8)

c.....lambda (absolute past).....................................
      n=nob+7
      do i=1,4
cc         pst(i,n)= x1*pst(i,ib)+y1*pst(i,jb)
         pst(i,n)= 0.d0
      enddo
      ii=mmod(ib+1)
      do while (mmod(ii+jb).ne.0)
         if(ish.ge.8)write(ifch,*) 'add lambda',ii
         do i=1,4
            pst(i,n)=pst(i,n)+pst(i,ii)
         enddo
         ii=mmod(ii+1)
      enddo
      do i=1,4
cc         pst(i,n+1)=(x2-x1)*pst(i,ib)
cc         pst(i,n+2)=(y2-y1)*pst(i,jb)
         pst(i,n+1)= pst(i,ib)
         pst(i,n+2)= pst(i,jb)
      enddo
      co(9)=     pf(n,n)
      co(10)= 2.*pf(n,n+1)
      co(11)= 2.*pf(n,n+2)
      co(12)= 2.*pf(n+1,n+2)
      if(ish.ge.8)write(ifch,*)'co lambda:',co(9),co(10),co(11),co(12)
      call utprix('gaksco',ish,ishini,8)
      end

c---------------------------------------------------------------------
      subroutine gakanp(ibr1,ibr,ibrr2,aml2,amr2,ala2,iok)
c---------------------------------------------------------------------
c   mass adjustment of fragments
c        ibr1-ibr   ibr-ibrr2
c   where ibr1,ibr,ibrr2 are break point indices
c   aml2,amr2 are the reqired squared masses (if zero -> any mass)
c   iok=0 (ok) or ok=1 (error)
c---------------------------------------------------------------------
#include "aaa.h"
      parameter (mxnbr=800,mxnba=10000,mxnin=2000)
      double precision pst!,pf
      common /gag/pst(4,0:mxnba),nob,ipst(0:mxnba)
     $     ,nbr,ijb(2,0:mxnbr),xyb(2,0:mxnbr)
     &     ,ptb(4,0:mxnbr),iflb(0:mxnbr),ip(0:mxnbr),co,x,y
      double precision ax,ay,az,ae,am,bx,by,bz,be,A,B,C
      dimension co(12),am2(0:2),nin(2,0:mxnin)
      logical weiter
c      pf(i,j)=pst(4,i)*pst(4,j)-pst(3,i)*pst(3,j)
c     &     -pst(2,i)*pst(2,j)-pst(1,i)*pst(1,j)
      mmod(i)=mod(mod(i,nob)+nob,nob)
      call utpri('gakanp',ish,ishini,6)

      ibr2=ibrr2
      if(ish.ge.6)write(ifch,*) ibr1,ibr,ibr2,aml2,amr2,ala2,iok
      iok=0
      ib=ijb(1,ibr)
      jb=ijb(2,ibr)
      ni=0
 10   do i=1,ni
        if((nin(1,i).eq.ib.and.nin(2,i).eq.jb)
     $       .or.(ipst(ib).eq.ipst(jb)))then
          iok=1
          if(ish.ge.4)then
            write(ifch,*) 'error ... endless loop'
            if(ib.eq.ipst(jb))  write(ifch,*) ' in zero mass region'
          endif
          goto 9999
        endif
      enddo
      ni=ni+1
      if(ni.gt.mxnin)stop'gakanp: increase parameter mxnin  '
      nin(1,ni)=ib
      nin(2,ni)=jb
      if(ish.ge.6)write(ifch,*)
      if(ish.ge.6)write(ifch,*) 'ib,jb=',ib,jb
      if(ni.ge.2)then
         if(ish.ge.6)write(ifch,*)'rotate pt to new band'
         if(ish.ge.6)write(ifch,*)'from',nin(1,ni-1),nin(2,ni-1)
         if(ish.ge.6)write(ifch,*)'  to',ib,jb
         if(ish.ge.6)write(ifch,'(4f12.8)') (ptb(i,ibr),i=1,4)
         ax=pst(1,nin(1,ni-1))+pst(1,nin(2,ni-1))
         ay=pst(2,nin(1,ni-1))+pst(2,nin(2,ni-1))
         az=pst(3,nin(1,ni-1))+pst(3,nin(2,ni-1))
         ae=pst(4,nin(1,ni-1))+pst(4,nin(2,ni-1))
         am=sngl(sqrt(max(1d-8,(ae+az)*(ae-az)-ax**2-ay**2))) !???????????????????????
         bx=pst(1,nin(1,ni-1))
         by=pst(2,nin(1,ni-1))
         bz=pst(3,nin(1,ni-1))
         be=pst(4,nin(1,ni-1))
         if(ish.ge.6)write (ifch,*) 'bx,by,bz,be',bx,by,bz,be
         call utlob2(1,ax,ay,az,ae,am,bx,by,bz,be,60)
         if(ish.ge.6)write (ifch,*) 'bx,by,bz,be',bx,by,bz,be
         call utlob4( 1,ax,ay,az,ae,am
     &        ,ptb(1,ibr),ptb(2,ibr),ptb(3,ibr),ptb(4,ibr))
         if(ish.ge.6)write(ifch,'(4f12.8)') (ptb(i,ibr),i=1,4)
         if(abs(bx)+abs(by)+abs(bz).gt.0.)then
           call utrot4( 1,bx,by,bz,ptb(1,ibr),ptb(2,ibr),ptb(3,ibr))
         else
           write(ifmt,*) 'null rot of pt (2)',bx,by,bz
           write(ifmt,'(4f12.8)') (ptb(i,ibr),i=1,4)
         endif
         if(ish.ge.6)write(ifch,'(4f12.8)') (ptb(i,ibr),i=1,4)
         ax=pst(1,ib)+pst(1,jb)
         ay=pst(2,ib)+pst(2,jb)
         az=pst(3,ib)+pst(3,jb)
         ae=pst(4,ib)+pst(4,jb)
         am=sngl(sqrt(max(1d-8,(ae+az)*(ae-az)-ax**2-ay**2))) !???????????????????????
         if(am.le.1.1e-4)then
           if(ish.ge.5)write(ifch,*)'error ... am<1.1e-4'
           iok=1
           goto 9999
         endif
         bx=pst(1,ib)
         by=pst(2,ib)
         bz=pst(3,ib)
         be=pst(4,ib)
         if(ish.ge.6)write (ifch,*) 'bx,by,bz,be',bx,by,bz,be
         if(ish.ge.6)write (ifch,*) 'ax,ay,az,ae',ax,ay,az,ae,am
         call utlob2(1,ax,ay,az,ae,am,bx,by,bz,be,60)
         if(ish.ge.6)write (ifch,*) 'bx,by,bz,be',bx,by,bz,be
         if(abs(bx)+abs(by)+abs(bz).gt.0.)then
           call utrot4(-1,bx,by,bz,ptb(1,ibr),ptb(2,ibr),ptb(3,ibr))
         else
           write(ifmt,*) 'null rot of pt (3)',bx,by,bz
           write(ifmt,'(4f12.8)') (ptb(i,ibr),i=1,4)
         endif
         if(ish.ge.6)write(ifch,'(4f12.8)') (ptb(i,ibr),i=1,4)
         call utlob4(-1,ax,ay,az,ae,am
     &        ,ptb(1,ibr),ptb(2,ibr),ptb(3,ibr),ptb(4,ibr))
         if(ish.ge.6)write(ifch,'(4f12.8)') (ptb(i,ibr),i=1,4)
      endif
      call gaksco(ibr1,ibr,ibr2,ib,jb,x1,x2,y1,y2)
c      if(ni.eq.1)print *,'------------------------------------',2
cc      x=(xyb(1,ibr)-x1)/(x2-x1)
cc      y=(xyb(2,ibr)-y1)/(y2-y1)
      x=xyb(1,ibr)
      y=xyb(2,ibr)
      am2(0)=aml2
      am2(1)=amr2
      am2(2)=ala2
      if(ish.ge.6)write(ifch,*) ibr1,ibr,ibr2,aml2,amr2,ala2,iok
      if(amr2.le.0.)then
         l1=2
         l2=0
      elseif(aml2.le.0.)then
         l1=1
         l2=2
      elseif(ala2.le.0.)then
         l1=1
         l2=0
      else
         stop' not like this , please...'
      endif
      if(ish.ge.6.and.amr2.le.0)write(ifch,*) 'adjust: 1',l1,l2
      if(ish.ge.6.and.aml2.le.0)write(ifch,*) 'adjust: 2' ,l1,l2
      if(ish.ge.6.and.ala2.le.0)write(ifch,*) 'adjust: 3',l1,l2
      i=4*l1
      j=4*l2
      A = dble(co(i+4))*dble(co(j+3)) - dble(co(j+4))*dble(co(i+3))
      B = dble(co(i+4))*dble(co(j+1)) - dble(co(i+3))*dble(co(j+2))
     &  + dble(co(i+2))*dble(co(j+3)) - dble(co(i+1))*dble(co(j+4))
     &  - dble(am2(l2))*dble(co(i+4)) + dble(am2(l1))*dble(co(j+4))
      C = dble(co(i+2))*dble(co(j+1)) - dble(co(i+1))*dble(co(j+2))
     &  + dble(am2(l1))*dble(co(j+2)) - dble(am2(l2))*dble(co(i+2))
      if (ish.ge.7) then
         write(ifch,*) 'ABC,q ',A,B,C,B**2-4.*A*C
         if(abs(A).gt.0.d0)then
            write(ifch,*) sqrt(max(0d0,B**2-4d0*A*C))/2d0/A-B/A/2d0
            write(ifch,*) -sqrt(max(0d0,B**2-4d0*A*C))/2d0/A-B/A/2d0
         endif
      endif
      x=0.
      y=0.
      xx=0.
      yy=0.
      if(abs(A).gt.1.d-20.and.B*B-4.*A*C.ge.0.d0)then
        y=sngl(sqrt(max(0.d0,B**2-4.*A*C))/2.d0/A-B/A/2.d0)
        if(abs(co(i+2)+y*co(i+4)).gt.0.)
     &       x=(am2(l1)-co(i+1)-y*co(i+3))/(co(i+2)+y*co(i+4))
      elseif(abs(A).le.1.d-20.and.abs(B).gt.0.d0)then
        y=-sngl(C/B)
        if(abs(co(i+2)+y*co(i+4)).gt.0.)
     &       x=(am2(l1)-co(i+1)-y*co(i+3))/(co(i+2)+y*co(i+4))
      else
        if(ish.ge.5)write(ifch,*)'error ... no solution of quadr equ'
        iok=1
        goto 9999
      endif
      if(abs(A).gt.1.d-20.and.B**2-4.*A*C.ge.0.d0)then
        yy=sngl(-sqrt(max(0.d0,B**2-4.*A*C))/2.d0/A-B/A/2.d0)
        if(abs(co(i+2)+yy*co(i+4)).gt.0.)
     &       xx=(am2(l1)-co(i+1)-yy*co(i+3))/(co(i+2)+yy*co(i+4))
      elseif(abs(A).le.1.d-20.and.abs(B).gt.0.d0)then
        yy=-sngl(C/B)
        if(abs(co(i+2)+yy*co(i+4)).gt.0.)
     &       xx=(am2(l1)-co(i+1)-yy*co(i+3))/(co(i+2)+yy*co(i+4))
      else
         if(ish.ge.5)write(ifch,*)'error ... no solution (2) '
         iok=1
         goto 9999
      endif
      if(ish.ge.6)then
         write(ifch,*) x ,y ,(co(i+2)+ y*co(i+4)),' OK '
         write(ifch,*) xx,yy,(co(i+2)+yy*co(i+4)),' OK '
      endif
      weiter=.true.
 50   if(x.gt.x1.and.x.lt.x2.and.y.gt.y1.and.y.lt.y2)then
cc         xyb(1,ibr)=x1+(x2-x1)*x
cc         xyb(2,ibr)=y1+(y2-y1)*y
         xyb(1,ibr)=x
         xyb(2,ibr)=y
         ijb(1,ibr)=ib
         ijb(2,ibr)=jb
         e1=pst(4,nob+1)-x*pst(4,nob+2)+y*pst(4,nob+3)
         e2=pst(4,nob+4)+x*pst(4,nob+5)-y*pst(4,nob+6)
         if( e1.lt.0. .or. e2.lt.0. ) then
           if(ish.ge.5)write(ifch,*)'error ... e1<0 or e2<0'
           iok=1
           goto 9999
         endif
         !amal2=co(1)+co(2)*x+co(3)*y+co(4)*x*y
         !amar2=co(5)+co(6)*x+co(7)*y+co(8)*x*y
         if(ish.ge.6)then
           amal2=co(1)+co(2)*x+co(3)*y+co(4)*x*y
           amar2=co(5)+co(6)*x+co(7)*y+co(8)*x*y
           write(ifch,*) 'brkshift:',xyb(1,ibr),xyb(2,ibr),ib,jb
           write (ifch,*)'E:',e1
           write (ifch,*)'E:',e2
           write(ifch,'(2(a6,1g12.6))') 'aml:'
     &          ,sqrt(max(0.,amal2)),'amr:',sqrt(max(0.,amar2))
         endif
         i=ibr1+1
         do while(i.le.ibr-1)
            if((mmod(ijb(1,i)+nob/2) .lt. mmod(ijb(1,ibr)+nob/2)
     &           .or.(mmod(ijb(1,i)+nob/2) .eq. mmod(ijb(1,ibr)+nob/2)
     &           .and.(xyb(1,i).gt.xyb(1,ibr))))
     &           .and.
     &           (ijb(2,i) .gt. ijb(2,ibr)
     &           .or.(ijb(2,i) .eq. ijb(2,ibr)
     &           .and.xyb(2,i).lt.xyb(2,ibr)))) goto 150
            if(ish.ge.6) then
               write(ifch,*) 'away:'
     &          ,i,xyb(1,i),xyb(2,i),ijb(1,i),ijb(2,i)
            endif
            call gakdel(i)
            i=i-1
            ibr=ibr-1
            ibr2=ibr2-1
 150        i=i+1
         enddo
         i=ibr+1
         do while (i.le.ibr2-1)
            if((mmod(ijb(1,i)+nob/2) .gt. mmod(ijb(1,ibr)+nob/2)
     &           .or.(mmod(ijb(1,i)+nob/2) .eq. mmod(ijb(1,ibr)+nob/2)
     &           .and.(xyb(1,i).lt.xyb(1,ibr))))
     &           .and.
     &           (ijb(2,i) .lt. ijb(2,ibr)
     &           .or.(ijb(2,i) .eq. ijb(2,ibr)
     &           .and.xyb(2,i).gt.xyb(2,ibr)))) goto 160
            if(ish.ge.6) then
               write(ifch,*) 'away:'
     &          ,i,xyb(1,i),xyb(2,i),ijb(1,i),ijb(2,i)
            endif
            call gakdel(i)
            ibr2=ibr2-1
            i=i-1
 160        i=i+1
         enddo
         goto 9999
      else
        if(x.gt.x2
     &    .and.ib.ne.ijb(1,ibr1) !brk-begrenzung
     &    .and.mmod(ib-1).ne.jb !linker oder rechter Rand
     &    .and.mmod(ib-1+jb).ne.0)then !oben oder unten
          ib=mmod(ib-1)
          goto 10
        endif
        if(x.lt.x1
     &    .and.ib.ne.ijb(1,ibr2) !brk-begrenzung
     &    .and.mmod(ib+1).ne.jb !linker oder rechter Rand
     &    .and.mmod(ib+1+jb).ne.0)then !oben oder unten
          ib=mmod(ib+1)
          goto 10
        endif
        if(y.gt.y2
     &    .and.jb.ne.ijb(2,ibr2) !brk-begrenzung
     &    .and.mmod(jb-1).ne.ib !linker oder rechter Rand
     &    .and.mmod(jb-1+ib).ne.0)then !oben oder unten
          jb=mmod(jb-1)
          goto 10
        endif
        if(y.lt.y1
     &    .and.jb.ne.ijb(2,ibr1) !brk-begrenzung
     &    .and.mmod(jb+1).ne.ib !linker oder rechter Rand
     &    .and.mmod(jb+1+ib).ne.0)then !oben oder unten
          jb=mmod(jb+1)
          goto 10
        endif
        if(weiter)then
          weiter=.false.
          x=xx
          y=yy
          goto 50
        endif
        if(ish.ge.5)write(ifch,*)'error ... x,y not in allowed range'
        iok=1
      endif
 9999 if(amr2.eq.0.) ibrr2=ibr2
      call utprix('gakanp',ish,ishini,6)
      end

cc----------------------------------------------------------------------
c      subroutine gakstr(ifl)
cc----------------------------------------------------------------------
cc
cc     calculates string-fragments by taking off pt of breakup
cc     do with ifl=1   undo with ifl=-1
cc
cc----------------------------------------------------------------------
c#include "aaa.h"
c      common /cpptr/ pptr(4,mxptl),ystr(mxptl)
c
c      do i=1,nptl
c        if(istptl(i).eq.29)then
c          nk1=ifrptl(1,i)
c          nk2=ifrptl(2,i)
c          do j=nk1,nk2
c            if ((istptl(j).eq.0.or.istptl(j-1).eq.0).and.j.ne.nk1) then
c              do k=1,4
c                pptl(k,j)=pptl(k,j)+pptr(k,j-1)*ifl
c              enddo
c              !write(ifch,*)"left side back  to ",j,(pptr(k,j-1),k=1,4)
c            endif
c            if ((istptl(j).eq.0.or.istptl(j+1).eq.0).and.j.ne.nk2) then
c              do k=1,4
c                pptl(k,j)=pptl(k,j)-pptr(k,j)*ifl
c              enddo
c              !write(ifch,*)"right side back to ",j,(-pptr(k,j),k=1,4)
c            endif
c            if(ifl.eq.-1.and.istptl(j).eq.0)then
c           e=pptl(1,j)**2+pptl(2,j)**2+pptl(3,j)**2
c     &           +pptl(5,j)**2
c              e=sqrt(e)
c           !dif=abs(e-pptl(4,j))
c           !if(dif.gt.0.01.and.dif/e.gt.0.1)print*,j,e,pptl(4,j)
c              pptl(4,j)=e
c            endif
c          enddo
c        endif
c        if(istptl(i).eq.0)then
c          if ( ifl.eq.1 ) then
c            ystr(i)=sign(1.,pptl(3,i))*alog((pptl(4,i)+abs(pptl(3,i)))
c     *           /sqrt(pptl(5,i)**2+pptl(1,i)**2+pptl(2,i)**2) )
c          endif
c        endif
c      enddo
c      end

c----------------------------------------------------------------------
      subroutine gakli2(nn1,nn2)
c----------------------------------------------------------------------

#include "aaa.h"
      double precision pgampr,rgampr
      common/cgampr/pgampr(5),rgampr(4)
c      double precision db1,db2,db3,db4,db5
      character label*8,idlabl*8
c      db1=0d0
c      db2=0d0
c      db3=rgampr(4)
c      db4=sqrt(rgampr(1)**2+rgampr(2)**2+rgampr(3)**2)
c      db5=sqrt( db4**2-rgampr(4)**2)
      n1=nn1
      n2=nn2
      if (n1.eq.0) n1=1
      if (n2.eq.0) n2=nptl
      write (ifch,'(1a4,5a12,4a4,a10,2a4)')
     &'no.','px','py','pz','E','m','ior','jor','if1','if2'
     &,'id','ist','ity'
      do i=n1,n2
         if (idptl(i).lt.10000)then
            label='        '
            label=idlabl(idptl(i))
         else
            label='        '
         endif
         chrg=0.
         if(iabs(idptl(i)).le.9999)call idchrg( 4 ,idptl(i),chrg)
         write (ifch,125) i,(pptl(j,i),j=1,5),iorptl(i),jorptl(i)
     &        ,ifrptl(1,i),ifrptl(2,i),idptl(i)
     $        ,chrg             !charge
     &        ,istptl(i),ityptl(i),label
      enddo
 125  format (1i4,5g18.10,4i6,1i10
     $     ,1f5.2               !charge
     $     ,2i4,'  ',A8
c     $     ,7g12.4,i5
     $     )
      return
      end

c----------------------------------------------------------------------
c      subroutine gakli4
cc----------------------------------------------------------------------
c
c#include "aaa.h"
c     parameter (mxnbr=800,mxnba=10000)
c      double precision pst
c      common /gag/pst(4,0:mxnba),nob,ipst(0:mxnba)
c     $     ,nbr,ijb(2,0:mxnbr),xyb(2,0:mxnbr)
c     &     ,ptb(4,0:mxnbr),iflb(0:mxnbr),ip(0:mxnbr),co,x,y
c      dimension co(12)
c      do i=0,nob-1
c        write (ifch,10) i,(pst(j,i),j=1,4)
c      enddo
c 10   format(1i4,5g18.10)
c      return
c      end
c
cc----------------------------------------------------------------------
c      subroutine gakli3
cc----------------------------------------------------------------------
c
c#include "aaa.h"
c      parameter (mxnbr=800,mxnba=10000)
c      double precision pst
c      common /gag/pst(4,0:mxnba),nob,ipst(0:mxnba)
c     $     ,nbr,ijb(2,0:mxnbr),xyb(2,0:mxnbr)
c     &     ,ptb(4,0:mxnbr),iflb(0:mxnbr),ip(0:mxnbr),co,x,y
c      dimension co(12),p1(5),p2(5)
c
c      write(ifch,*) 'particle list of string decay'
c      do i=1,5
c        p1(i)=0.
c      enddo
c      do i=1,nbr+1
c        if(i.lt.nbr+1)then
c          call gaksco(i-1,i,i+1,ijb(1,i),ijb(2,i),x1,x2,y1,y2)
c          if(x2.gt.x1)then
ccc            x=(xyb(1,i)-x1)/(x2-x1)
c            x=xyb(1,i)
c          else
c            x=0.
c          endif
c          if(y2.gt.y1)then
ccc            y=(xyb(2,i)-y1)/(y2-y1)
c            y=xyb(2,i)
c          else
c            y=0.
c          endif
c          aml2=co(1)+x*co(2)+y*co(3)+x*y*co(4)
c          amr2=co(5)+x*co(6)+y*co(7)+x*y*co(8)
c          aml=sign(sqrt(abs(aml2)),aml2)
c          amr=sign(sqrt(abs(amr2)),amr2)
c          do j=1,4
c            p2(j)=pst(j,nob+1)-x*pst(j,nob+2)+y*pst(j,nob+3)
c            p1(j)=p1(j)+p2(j)
c          enddo
c          p2(5)=aml
c        else
c          do j=1,4
c            p2(j)=pst(j,nob+4)+x*pst(j,nob+5)-y*pst(j,nob+6)
c            p1(j)=p1(j)+p2(j)
c          enddo
c          p2(5)=amr
c        endif
c        write(ifch,'(2i4,i6,a,i5,i10,5g14.6)') i-1,i
c     &    ,-iflb(i-1),'==',iflb(i),ip(i)
c     &    ,(p2(j),j=1,5)
c      enddo
c      am2=p1(4)**2-p1(3)**2-p1(2)**2-p1(1)**2
c      p1(5)=sign(sqrt(abs(am2)),am2)
c      write(ifch,'(12x,a60)')
c     &  '------------------------------------------------------------'
c      write(ifch,'(14x,5g14.6)') (p1(j),j=1,5)
c      write(ifch,*)
c
c      end
c

c---------------------------------------------------------------------
      subroutine idress(id,am,idr,iadj) !for mass adjustment
c---------------------------------------------------------------------
#include "aaa.h"
      amb=am
      call idres(id,am,idr,iadj,1)
      if(idr.eq.0)then
        return
      endif
      ids=max(mod(iabs(id)/100,10),mod(iabs(id)/10,10))
      ids=max(ids,mod(iabs(id)/1000,10))
      !if(ids.ge.4)write (ifch,'(2x,3i6,2f7.2,$)')id,idr,iadj,amb,am
      if(iabs(idr).le.999) then !~~~~~~~~~mesons
        if(ids.le.4)goto 777  !???? if the following is used, bad result
        if(ids.le.2)then
          idr=sign(iabs(id)+int(rangen()+0.5),id)
        elseif(ids.eq.3)then
          idr=sign(iabs(id)+int(rangen()+0.6),id)
        elseif(ids.eq.4)then
          idr=sign(iabs(id)+int(rangen()+0.75),id)
        else
          idr=sign(iabs(id)+int(rangen()+0.875),id)
        endif
        call idmass(idr,am)
      elseif(iabs(idr).le.9999)then !~~~~~~~~~baryons
        if(ids.le.3)goto 777
        if(mod(iabs(idr),10).gt.1)then
          if(iabs(id).ne.1111.and.iabs(id).ne.2221.and.iabs(id).ne.3331)
     $         then
            idr=sign(iabs(id)+1,id)
            call idmass(idr,am)
          else
            idr=id
            call idmass(idr,am)
          endif
        endif
      endif

  777 continue
      !if(ids.ge.4)write (ifch,'(2x,a,2i6)') '->',id,idr
       
      !if(ids.ge.4)iadj=0 !to avoid break point adjustments ! <========

      end


c---------------------------------------------------------------------
      SUBROUTINE gaksphe(sphe,r,mstu41)

C...Purpose: to perform sphericity tensor analysis to give sphericity,
C...aplanarity and the related event axes. stolen from jetset ;-)
#include "aaa.h"
      dimension sphe(4,3)
      DIMENSION SM(3,3),SV(3,3)

C...Calculate matrix to be diagonalized.
      NP=0
      JA=0
      JB=0
      JC=0
c      MSTU41=2
      DO 110 J1=1,3
         DO 100 J2=J1,3
            SM(J1,J2)=0.
 100     CONTINUE
 110  CONTINUE
      PS=0.
      DO 140 I=1,nptl
      IF(istptl(i).ne.0)  GOTO 140
      IF(MSTU41.GE.2) THEN
         ida=iabs(idptl(i))
         IF(ida.EQ.0.OR.ida.EQ.11.OR.ida.EQ.13.OR.ida.EQ.15) GOTO 140
         IF(MSTU41.GE.3) then
            call idchrg( 5 ,idptl(i),chrg)
            if (abs(chrg).le.0.1) goto 140
         endif
      ENDIF
      NP=NP+1
      PA=SQRT(pptl(1,i)**2+pptl(2,I)**2+pptl(3,i)**2)
      PWT=1.
      IF(ABS(r-2.).GT.0.001) PWT=MAX(1E-10,PA)**(r-2.)
      DO 130 J1=1,3
         DO 120 J2=J1,3
            SM(J1,J2)=SM(J1,J2)+PWT*pptl(j1,i)*pptl(j2,i)
 120     CONTINUE
 130  CONTINUE
      PS=PS+PWT*PA**2
 140  CONTINUE

C...Very low multiplicities (0 or 1) not considered.
      IF(NP.LE.1) THEN
        if(ish.ge.1)then
          CALL utmsg('sphe  ')
          write(ifch,*) 'too few particles for analysis'
          call utmsgf
        endif
        sphe(4,1)=-1.
        RETURN
      ENDIF
      DO 160 J1=1,3
         DO 150 J2=J1,3
            SM(J1,J2)=SM(J1,J2)/PS
 150     CONTINUE
 160  CONTINUE

C...Find eigenvalues to matrix (third degree equation).
      SQ=(SM(1,1)*SM(2,2)+SM(1,1)*SM(3,3)+SM(2,2)*SM(3,3)-SM(1,2)**2-
     &     SM(1,3)**2-SM(2,3)**2)/3.-1./9.
      SR=-0.5*(SQ+1./9.+SM(1,1)*SM(2,3)**2+SM(2,2)*SM(1,3)**2+SM(3,3)*
     &     SM(1,2)**2-SM(1,1)*SM(2,2)*SM(3,3))
     &     +SM(1,2)*SM(1,3)*SM(2,3)+1./27.
      SP=COS(ACOS(MAX(MIN(SR/SQRT(-SQ**3),1.),-1.))/3.)
      sphe(4,1)=1./3.+SQRT(-SQ)*MAX(2.*SP,SQRT(3.*(1.-SP**2))-SP)
      sphe(4,3)=1./3.+SQRT(-SQ)*MIN(2.*SP,-SQRT(3.*(1.-SP**2))-SP)
      sphe(4,2)=1.-sphe(4,1)-sphe(4,3)
      IF(sphe(4,2).LT.1E-5) THEN
        if(ish.ge.1)then
          CALL utmsg('gaksphe')
          write(ifch,*) 'all particles back-to-back'
          call utmsgf
        endif
        sphe(4,1)=-1.
        RETURN
      ENDIF

C...Find first and last eigenvector by solving equation system.
      DO 240 I=1,3,2
         DO 180 J1=1,3
            SV(J1,J1)=SM(J1,J1)-sphe(4,I)
            DO 170 J2=J1+1,3
               SV(J1,J2)=SM(J1,J2)
               SV(J2,J1)=SM(J1,J2)
 170        CONTINUE
 180     CONTINUE
         SMAX=0.
         DO 200 J1=1,3
            DO 190 J2=1,3
               IF(ABS(SV(J1,J2)).LE.SMAX) GOTO 190
               JA=J1
               JB=J2
               SMAX=ABS(SV(J1,J2))
 190        CONTINUE
 200     CONTINUE
         SMAX=0.
         DO 220 J3=JA+1,JA+2
            J1=J3-3*((J3-1)/3)
            RL=SV(J1,JB)/SV(JA,JB)
            DO 210 J2=1,3
               SV(J1,J2)=SV(J1,J2)-RL*SV(JA,J2)
               IF(ABS(SV(J1,J2)).LE.SMAX) GOTO 210
               JC=J1
               SMAX=ABS(SV(J1,J2))
 210        CONTINUE
 220     CONTINUE
         JB1=JB+1-3*(JB/3)
         JB2=JB+2-3*((JB+1)/3)
         sphe(JB1,I)=-SV(JC,JB2)
         sphe(jb2,I)=SV(JC,JB1)
         sphe(jb,I)=-(SV(JA,JB1)*sphe(jb1,I)+SV(JA,JB2)
     &        *sphe(jb2,I))/SV(JA,JB)
         PA=SQRT(sphe(1,I)**2+sphe(2,I)**2+sphe(3,I)**2)
         SGN=(-1.)**INT(rangen()+0.5)
         DO 230 J=1,3
            sphe(j,I)=SGN*sphe(j,I)/PA
 230     CONTINUE
 240  CONTINUE

C...Middle axis orthogonal to other two. Fill other codes.
      SGN=(-1.)**INT(rangen()+0.5)
      sphe(1,2)=SGN*(sphe(2,1)*sphe(3,3)-sphe(3,1)*sphe(2,3))
      sphe(2,2)=SGN*(sphe(3,1)*sphe(1,3)-sphe(1,1)*sphe(3,3))
      sphe(3,2)=SGN*(sphe(1,1)*sphe(2,3)-sphe(2,1)*sphe(1,3))

      do i=1,3
         do j=1,4
            pptl(j,nptl+i)=sphe(j,i)
         enddo
      enddo


C...Calculate sphericity and aplanarity. Select storing option.
ccc      SPH=1.5*(sphe(4,2)+sphe(4,3))
ccc      APL=1.5*sphe(4,3)

      RETURN
      END

C*********************************************************************

      SUBROUTINE gakthru(thru,mstu41)

C...Purpose: to perform thrust analysis to give thrust, oblateness
C...and the related event axes. stolen from jetset ;-)
#include "aaa.h"
      DIMENSION TDI(3),TPR(3)
      dimension thru(4,3)

C...Take copy of particles that are to be considered in thrust analysis.
      IAGR=0
      NP=0
      PS=0.
c      MSTU41=2
      MSTU44=4
      MSTU45=2
      PARU42=1.0
      PARU48=0.0000001
      DO 100 I=1,nptl
         IF(istptl(i).ne.0)goto 100
         ida=iabs(idptl(i))
         IF(ida.EQ.0.OR.ida.EQ.11.OR.ida.EQ.13.OR.ida.EQ.15)GOTO 100
         IF(MSTU41.GE.3) then
            call idchrg( 6 ,idptl(i),chrg)
            if (abs(chrg).le.0.1) goto 100
         endif

         IF(nptl+NP.GE.mxptl) THEN
            CALL utstop('gakthru: no more memory left in cptl&')
            thru(4,1)=-1.
            RETURN
         ENDIF
         NP=NP+1
c         K(nptl+NP,1)=23
         pptl(1,nptl+NP)=pptl(1,I)
         pptl(2,nptl+NP)=pptl(2,I)
         pptl(3,nptl+NP)=pptl(3,I)
         pptl(4,nptl+NP)=SQRT(pptl(1,I)**2+pptl(2,I)**2+pptl(3,I)**2)
         pptl(5,nptl+NP)=1.
         IF(ABS(PARU42-1.).GT.0.001)
     &        pptl(5,nptl+NP)=pptl(4,nptl+NP)**(PARU42-1.)
         PS=PS+pptl(4,nptl+NP)*pptl(5,nptl+NP)
 100  CONTINUE

C...Very low multiplicities (0 or 1) not considered.
      IF(NP.LE.1) THEN
         CALL utmsg('thru  ')
         write(ifch,*) 'too few particles for analysis'
         call utmsgf
         thru(4,1)=-1
         RETURN
      ENDIF

C...Loop over thrust and major. T axis along z direction in latter case.
      DO 320 ILD=1,2
      IF(ILD.EQ.2) THEN
c        PHI=GAKANG(pptl(1,nptl+NP+1),pptl(2,nptl+NP+1))
c        CALL lurot(nptl+1,nptl+NP+1,0.,-PHI)
c        THE=GAKANG(pptl(3,nptl+NP+1),pptl(1,nptl+NP+1))
c        CALL lurot(nptl+1,nptl+NP+1,-THE,0.)
         ax=pptl(1,nptl+NP+1)
         ay=pptl(2,nptl+NP+1)
         az=pptl(3,nptl+NP+1)
         do irot=nptl+1,nptl+NP+1
            call utrota(1,ax,ay,az
     &           ,pptl(1,irot),pptl(2,irot),pptl(3,irot))
         enddo
         if(np.eq.2)then
           pptl(1,nptl+NP+2)=1.
           pptl(2,nptl+NP+2)=0.
           pptl(3,nptl+NP+2)=0.
           pptl(4,nptl+NP+2)=0.
           goto 325
         endif
      ENDIF

C...Find and order particles with highest p (pT for major).
      DO 110 ILF=nptl+NP+4,nptl+NP+MSTU44+4
         pptl(4,ILF)=0.
 110  CONTINUE
      DO 160 I=nptl+1,nptl+NP
         IF(ILD.EQ.2) pptl(4,I)=SQRT(pptl(1,I)**2+pptl(2,I)**2)
         DO 130 ILF=nptl+NP+MSTU44+3,nptl+NP+4,-1
            IF(pptl(4,I).LE.pptl(4,ILF)) GOTO 140
            DO 120 J=1,5
               pptl(j,ILF+1)=pptl(j,ILF)
 120        CONTINUE
 130     CONTINUE
         ILF=nptl+NP+3
 140     DO 150 J=1,5
            pptl(j,ILF+1)=pptl(j,I)
 150     CONTINUE
 160  CONTINUE

C...Find and order initial axes with highest thrust (major).
      DO 170 ILG=nptl+NP+MSTU44+5,nptl+NP+MSTU44+15
         pptl(4,ILG)=0.
 170  CONTINUE
      NC=2**(MIN(MSTU44,NP)-1)
      DO 250 ILC=1,NC
         DO 180 J=1,3
            TDI(J)=0.
 180     CONTINUE
         DO 200 ILF=1,MIN(MSTU44,NP)
            SGN=pptl(5,nptl+NP+ILF+3)
            IF(2**ILF*((ILC+2**(ILF-1)-1)/2**ILF).GE.ILC) SGN=-SGN
            DO 190 J=1,4-ILD
               TDI(J)=TDI(J)+SGN*pptl(j,nptl+NP+ILF+3)
 190        CONTINUE
 200     CONTINUE
         TDS=TDI(1)**2+TDI(2)**2+TDI(3)**2
         DO 220 ILG=nptl+NP+MSTU44+MIN(ILC,10)+4,nptl+NP+MSTU44+5,-1
            IF(TDS.LE.pptl(4,ILG)) GOTO 230
            DO 210 J=1,4
               pptl(j,ILG+1)=pptl(j,ILG)
 210        CONTINUE
 220     CONTINUE
         ILG=nptl+NP+MSTU44+4
 230     DO 240 J=1,3
            pptl(j,ILG+1)=TDI(J)
 240     CONTINUE
         pptl(4,ILG+1)=TDS
 250  CONTINUE

C...  Iterate direction of axis until stable maximum.
      pptl(4,nptl+NP+ILD)=0.
      ILG=0
 260  ILG=ILG+1
      THP=0.
 270  THPS=THP
      DO 280 J=1,3
         IF(THP.LE.1E-10) TDI(J)=pptl(j,nptl+NP+MSTU44+4+ILG)
         IF(THP.GT.1E-10) TDI(J)=TPR(J)
         TPR(J)=0.
 280  CONTINUE
      DO 300 I=nptl+1,nptl+NP
         SGN=SIGN(pptl(5,I),TDI(1)*pptl(1,I)
     &        +TDI(2)*pptl(2,I)+TDI(3)*pptl(3,I))
         DO 290 J=1,4-ILD
            TPR(J)=TPR(J)+SGN*pptl(j,I)
 290     CONTINUE
 300  CONTINUE
      THP=SQRT(TPR(1)**2+TPR(2)**2+TPR(3)**2)/PS
      IF(THP.GE.THPS+PARU48) GOTO 270

C...  Save good axis. Try new initial axis until a number of tries agree.
      IF(THP.LT.pptl(4,nptl+NP+ILD)-PARU48.AND.ILG.LT.MIN(10,NC))
     &     GOTO 260
      IF(THP.GT.pptl(4,nptl+NP+ILD)+PARU48)
     $     THEN
         IAGR=0
         SGN=(-1.)**INT(rangen()+0.5)
         DO 310 J=1,3
            pptl(j,nptl+NP+ILD)=SGN*TPR(J)/(PS*THP)
 310     CONTINUE
         pptl(4,nptl+NP+ILD)=THP
         pptl(5,nptl+NP+ILD)=0.
      ENDIF
      IAGR=IAGR+1
      IF(IAGR.LT.MSTU45.AND.ILG.LT.MIN(10,NC)) GOTO 260
 320  CONTINUE

C...  Find minor axis and value by orthogonality.
 325  SGN=(-1.)**INT(rangen()+0.5)
      pptl(1,nptl+NP+3)=-SGN*pptl(2,nptl+NP+2)
      pptl(2,nptl+NP+3)=SGN*pptl(1,nptl+NP+2)
      pptl(3,nptl+NP+3)=0.
      THP=0.
      DO 330 I=nptl+1,nptl+NP
         THP=THP+pptl(5,I)*ABS(pptl(1,nptl+NP+3)*pptl(1,I)
     &        +pptl(2,nptl+NP+3)*pptl(2,I))
 330  CONTINUE
      pptl(4,nptl+NP+3)=THP/PS
      pptl(5,nptl+NP+3)=0.


C...  Fill axis information. Rotate back to original coordinate system.
      do irot=nptl+NP+1,nptl+NP+3
         call utrota(-1,ax,ay,az
     &        ,pptl(1,irot),pptl(2,irot),pptl(3,irot))
      enddo

      do ild=1,3
         do j=1,4
            thru(j,ild)=pptl(j,nptl+NP+ild)
         enddo
      enddo

C...Calculate thrust and oblateness. Select storing option.
ccc      THR=thru(4,1)
ccc      OBL=thru(4,2)-thru(4,3)

      RETURN
      END

      subroutine gakjet(ijadu)
#include "aaa.h"
      common/njjjj/njets(5,2,mxbins)
c      nmin=xpar1
c      nmax=xpar2
      if(nrevt.eq.1)then
        do i=1,5
          do j=1,mxbins
            njets(i,ijadu,j)=0
          enddo
        enddo
      endif
      do i=1,nrbins
        ycut=xminim*(xmaxim/xminim)**((float(i)-0.5)/nrbins)
c        if(iologe.ne.1)ycut=xminim+(xmaxim-xminim)/nrbins*(nrbins)
        nj=gaknjt(ycut,ijadu)
        if(nj.ge.1.and.nj.le.5)then
          njets(nj,ijadu,i)=njets(nj,ijadu,i)+1
        endif
      enddo
      end

      subroutine gakjto
#include "aaa.h"
      common/njjjj/njets(5,2,mxbins)
      n=xpar4
      ijadu=xpar3
      do i=1,nrbins
        ycut=xminim*(xmaxim/xminim)**((float(i)-0.5)/nrbins)
        write (ifhi,*) ycut,float(njets(n,ijadu,i))/nrevt
     $       ,sqrt(1.*njets(n,ijadu,i))/nrevt
      enddo
      end
C*********************************************************************
      function gaknjt(ycut,ijadu)
c
c     ijadu   1 =  JADE    2=DURHAM
c     ycut - max. distance
c
#include "aaa.h"

      a2j(i,j)=2.*pptl(4,i)*pptl(4,j)*(1.-(pptl(1,i)*pptl(1,j)
     &     +pptl(2,i)*pptl(2,j)+pptl(3,i)*pptl(3,j))
     &     /(sqrt(pptl(1,i)**2+pptl(2,i)**2+pptl(3,i)**2)
     &       *sqrt(pptl(1,j)**2+pptl(2,j)**2+pptl(3,j)**2)))/evis**2

      a2d(i,j)=2.*min(pptl(4,i)**2,pptl(4,j)**2)
     &     *(1.-(pptl(1,i)*pptl(1,j)
     &     +pptl(2,i)*pptl(2,j)+pptl(3,i)*pptl(3,j))
     &     /(sqrt(pptl(1,i)**2+pptl(2,i)**2+pptl(3,i)**2)
     &       *sqrt(pptl(1,j)**2+pptl(2,j)**2+pptl(3,j)**2)))/evis**2

      bet(i)=sqrt(pptl(1,i)**2+pptl(2,i)**2+pptl(3,i)**2)

      ska(i,j)=pptl(1,i)*pptl(1,j)+pptl(2,i)*pptl(2,j)
     &     +pptl(3,i)*pptl(3,j)

      a2c(i,j)= ((bet(i)*bet(j)-ska(i,j))
     &     *2.*bet(i)*bet(j)  /  (0.00001+bet(i)+bet(j))**2 )

      evis=0.
      nn=0
      do i=1,nptl
         if (istptl(i).eq.0) then
            nn=nn+1
            do j=1,5
               pptl(j,nptl+nn)=pptl(j,i)
            enddo
            evis=evis+pptl(4,i)
            jorptl(i)=nn
         endif
      enddo
      iflag=0
      i1=0
      j1=0
      do while (iflag.eq.0.and.nn.ge.2)
         a2min=ycut
         iflag=1
         do i=nptl+1,nptl+nn-1
            do j=i+1,nptl+nn
               if(ijadu.eq.1)then
                  a2=a2j(i,j)
               elseif(ijadu.eq.2) then
                  a2=a2d(i,j)
               else
                  a2=a2c(i,j)
               endif
               if (a2.lt.a2min) then
                  iflag=0
                  i1=i
                  j1=j
                  a2min=a2
               endif
            enddo
         enddo
         if (iflag.eq.0) then
            do i=1,4
               pptl(i,i1)=pptl(i,i1)+pptl(i,j1)
            enddo
            do i=1,nptl
               if(istptl(i).eq.0.and.jorptl(i).eq.j1-nptl)
     &              jorptl(i)=i1-nptl
               if(istptl(i).eq.0.and.jorptl(i)+nptl.gt.j1)
     &              jorptl(i)=jorptl(i)-1
            enddo
            do i=j1,nptl+nn
               do j=1,5
                  pptl(j,i)=pptl(j,i+1)
               enddo
            enddo
            nn=nn-1
         endif
      enddo
      do i=nptl+1,nptl+nn
         istptl(i)=-1
         jorptl(i)=i-nptl
         pptl(5,i)=sqrt(max(0.,(pptl(4,i)+pptl(3,i))
     &   *(pptl(4,i)-pptl(3,i))-pptl(2,i)**2-pptl(1,i)**2))
      enddo
      do i=nptl+1,nptl+nn-1
         do j=i+1,nptl+nn
            if(pptl(4,jorptl(i)+nptl).lt.pptl(4,jorptl(j)+nptl))then
               k=jorptl(i)
               jorptl(i)=jorptl(j)
               jorptl(j)=k
            endif
         enddo
      enddo
      do i=nptl+1,nptl+nn
         idptl(nptl+jorptl(i))=9910+i-nptl
      enddo
      do i=1,nptl
        jorptl(i)=0             !jorptl(nptl+jorptl(i))
      enddo
c      nptl=nptl+nn

      gaknjt=nn
      return
      end

      subroutine idtr5(id,ic)  ! for string ends (+-q, +-qq)
      integer ic(2)
      ic(1)=0
      ic(2)=0
      ii=1
      if(id.lt.0)ii=2
      i1=1
      if(iabs(id).gt.999)i1=3
      do i=i1,int(log(abs(float(id)))/log(10.))+1
        j=mod(iabs(id)/10**(i-1),10)
        if(j.gt.0)then
          ic(ii)=ic(ii)+10**(6-j)
        endif
      enddo
      return
      end

      function ammin(id1,id2)
      dimension ic1(2),ic2(2),jc2(6,2),jc1(6,2)
      call idtr5(id1,ic1)
      call idtr5(id2,ic2)
      call idcomk(ic1)
      call idcomk(ic2)
      call iddeco(ic1,jc1)
      call iddeco(ic2,jc2)
      ammin=utamnx(jc1,jc2)
      end

      function idsp(id1x,id2x)
      id1=id1x
      id2=id2x
      do i=1,2
        if(i.eq.1)id=id1  
        if(i.eq.2)id=id2  
        ia=iabs(id) 
        if(ia.gt.1000.and.mod(ia,100).eq.0)then
          ib=ia/1000
          ic=mod(ia/100,10)
          if(ib.gt.ic)then !isospin 1 diquark
            ie=(ic*1000+ib*100)*id/ia
            !print*,ia,ib,ic,ie
            if(i.eq.1)id1=ie  
            if(i.eq.2)id2=ie  
          endif
        endif
      enddo 
      ia1=iabs(id1)
      ia2=iabs(id2)
      if(ia1.ge.1000.and.ia2.ge.1000)then
        idsp=0
        isign=0
      elseif(ia1.le.1000.and.ia2.le.1000)then
        idsp=min(ia1,ia2)*100+max(ia1,ia2)*10
        isign=1
        if(max(ia1,ia2).ne.-min(id1,id2)) isign = -1
        if(idsp.eq.220)idsp=110
        if(idsp.eq.330)idsp=220
      else
        isign=1
        if(id1.lt.0.and.id2.lt.0)isign=-1
        idb=min(ia1,ia2)
        if(idb.eq.5)then
          idsp=0
          return
        endif
        ida=max(ia1,ia2)
        ida1=ida/1000
        ida2=mod(ida/100,10)
        if(idb.le.ida1)then
          idsp=idb*1000+ida/10
        elseif(idb.le.ida2)then
          idsp=ida1*1000+idb*100+ida2*10
        else
          idsp=ida+idb*10
        endif
        if(ida1.eq.ida2.and.ida2.eq.idb)idsp=idsp+1
      endif
      idsp=idsp*isign
      !if(id2.ne.id2x)print*,'TEST ',id2x,id2,idsp
      !if(id2.ne.id2x)stop
      return
      end

