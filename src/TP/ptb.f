C
C  This file is part of EPOS4
C  Copyright (C) 2022 research institutions and authors (See CREDITS file)
C  This file is distributed under the terms of the GNU General Public License version 3 or later
C  (See COPYING file for the text of the licence)
C



!   remark (KW 29.05.2013)
!
!   some basic paremeters depend on q2min
!    I tabulated them (without checking their use)
!
!    so please check !!!!!!!!!!!!!!!!



c-----------------------------------------------------------------------
      subroutine ipoBasics(q2imin)
c-----------------------------------------------------------------------
#include "aaa.h"
#include "sem.h"
#include "tab.h"
c      common /psar36/ alvc
c      real wi(0:1)
      real wi(0:2)
      parameter(mxbaspar=10+2*nclha)
      common/cbaspar/baspar(maxq2mn,mxbaspar)
      common /psar40/ coefxu1(2,nclha)
      double precision coefxu1
      real q2imin(2)

      if(q2imin(1).le.0.)then
        call setBasicq(q2nmin,1)
        call setBasicq(q2nmin,2)
        return 
      endif
      
c      call setBasicq(q2imin(1),1)
c      call setBasicq(q2imin(2),2)
c        print *,'basics',1,q2imin(1),dels(1),betpom(1),alpff(1)!,coefxu1(ir,2)
      return       !20200107 with Esatur, parameters do not depend on Q2s anymore

      do ir=1,2

        q=log(max(1.e-6,q2imin(ir)))
c        ipol=1
c        call ipoli(q,nq,frac)
c        wi(0)=1.-frac
c        wi(1)=frac
        ipol=2
        call ipoli3(q,nq,frac)
        wi(1)=frac
        wi(2)=wi(1)*(wi(1)-1.)*.5
        wi(0)=1.-wi(1)+wi(2)
        wi(1)=wi(1)-2.*wi(2)
     
        alpff(ir)  = 0.
c        betff(ir)  = 0.
c        betpom(ir) = 0.
c        dels(ir)   = 0.
c        gamsoft(ir)= 0.
c        do k=1,nclha
c          coefxu1(ir,k)=0d0
c        enddo
        do n=0,ipol
        alpff(ir)  = alpff(ir)+wi(n)*baspar(nq+n,ir)
c        betff(ir)  = betff(ir)+wi(n)*baspar(nq+n,ir+2)
c        betpom(ir) = betpom(ir)+wi(n)*baspar(nq+n,ir+4)
c        dels(ir)   = dels(ir)+wi(n)*baspar(nq+n,ir+6)
c        gamsoft(ir)= gamsoft(ir)+wi(n)*baspar(nq+n,ir+8)
c        do k=1,nclha
c          coefxu1(ir,k)=coefxu1(ir,k)
c     .                 +dble(wi(n)*baspar(nq+n,10+(ir-1)*nclha+k))
c        enddo
        enddo
c        alpqua(ir)=(1.-betff(ir))/2.
c        if(ir.eq.1)then
c          icl=iclpro
c          if(icl.eq.4)icl=1
c        else
c          icl=icltar
c        endif
c        gamhad(icl)=alpff(ir)
c        print *,'basics',ir,q2imin(ir),dels(ir),betpom(ir),alpff(ir)!,coefxu1(ir,2)
      enddo

      end

c----------------------------------------------------------------------
      subroutine setBasics
c----------------------------------------------------------------------
#include "aaa.h"
#include "ems.h"
#include "par.h"
#include "sem.h"
#include "tab.h"
      common /psar2/  edmax,epmax
c      common /psar7/  delx,alam3p,gam3p
      common /psar36/ alvc,qnorm,qnormp
      common /ar3/    x1(7),a1(7)
      double precision fzeroGluZZ,FzeroQuaZZ
      common /psar9/  alpr
      common /ckinirj/kinirj
      common /testj/  ajeth(4),ajete(5),ajet0(7)
      parameter(nbkbin=100)
      common /kfitd/ xkappafit(nbkbin,nclegy,nclha,nclha),xkappa,bkbin
      parameter(mxbaspar=10+2*nclha)
      common/cbaspar/baspar(maxq2mn,mxbaspar)
      common /psar40/ coefxu1(2,nclha)
      double precision coefxu1
      dimension q2imin(2)

      if(inicnt.eq.1)then
c        q2nmin=q2dmin !q2mnval(1)      !minimum for Q2s
        if(q2nmin.lt.q2zmin)
     .  call utstop('Minimum Q2 should be larger than q2zmin.&')
        difucb=0.               !bg
        difub=0.                !bg
        difuub=0.               !bg
        difudb=0.               !bg
        difusb=0.               !bg
        difubb=0.               !bg
        if(nrflav.le.3)then
          difuc=0.
          difuuc=0.
          difudc=0.
          difusc=0.
          difucc=0.
          rstrac(1)=0.
          rstrac(2)=0.
          rstrac(3)=0.
c          rstrac(4)=0.
        endif
c       set maximum value of Q2s to read in tables (to save time)
c       Warning at running time if value exceeded
        maxq2mx=min(maxq2mn
     .             ,5+nint(log(max(1,maproj*matarg)*max(engy,engmax))))
        maxq2mx=maxq2mn      !just use the maximum for the time being (2019.05.13)
      endif

      if(iappl.ne.6)then

      do i=1,4
      ajeth(i)=0.
      enddo
      do i=1,5
      ajete(i)=0.
      ajet0(i)=0.
      enddo
      ajet0(6)=0.
      ajet0(7)=0.


      if(isetcs.le.1)then              !for Kfit
c      if(isetcs.le.2)then              !for Kfit  ????????? to speed up
        bkbin=0.1
      else
        bkbin=0.05
      endif
      xkappa=1.

      do i=1,nclha
      do j=1,nclha
      do k=1,4
        chadr(k,i,j)=1. !make sure diff. renormalization is set to 1. first
      enddo
      enddo
      enddo


      
c fix enhanced diagrams at minimum energy = 2.5
c      delx=1.5 !sqrt(egymin*egymin/exp(1.))
c arbitrary value for alam3p (not good if too small (infinite loop in rsh))
c      alam3p=0.5*(r2had(1)+r2had(2)+r2had(3)) !0.6
c      gam3p=.1



      if(iappl.eq.1.or.iappl.eq.2.or.iappl.eq.8.or.iappl.eq.9)then

        !parton density normalization

        sq=log(log(q2nmin/.232**2)/log(.23/.232**2))
        du=2.997+.753*sq-.076*sq*sq
        qnorm=0.
        do i=1,7
        do m=1,2
          xx=.5+x1(i)*(m-1.5)
          xxq=1.-xx**(1./(1.+du))
          qnorm=qnorm+a1(i)*(psdpdf(xxq,q2nmin,0.,2,1)
     *                      +psdpdf(xxq,q2nmin,0.,2,2))/(1.-xxq)**du
        enddo
        enddo
        qnorm=qnorm*.5/(1.+du)
        qnormp=qnorm

        if(ish.ge.4) write (ifch,*)'qnorm tar',qnorm

        if(iclpro.ne.icltar)then

        iclp=iclpro
        if(iclp.eq.4)iclp=1
        sq=log(log(q2nmin/.232**2)/log(.25/.232**2))
        dpi=.367+.563*sq
        qnorm=0.
        do i=1,7
        do m=1,2
          xx=.5+x1(i)*(m-1.5)
          xxq=1.-xx**(1./(1.+dpi))
          qnorm=qnorm+a1(i)*(psdfh4(xxq,q2nmin,0.,iclp,1)
     *                      +psdfh4(xxq,q2nmin,0.,iclp,2))/(1.-xxq)**dpi
        enddo
        enddo
        qnorm=qnorm*.5/(1.+dpi)
        endif

        if(ish.ge.4) write (ifch,*)'qnorm pro',qnorm

      endif
    
      r3pomti=dble(r3pom)/dble(sqrt(8.*pi))

      stmass=.05               !string mass cutoff

      do ir=1,2       ! ======= proj/tar loop ========
        nq2mn1=1
        nq2mn2=maxq2mn
        if(iscreen.eq.0)then    !maximum q2min given by q2nmin
          q=sqrt(q2nmin)
          call ipoli(q,nq,frac)
          nq2mn2=nq+1
        endif
        do nq2mn=nq2mn1,nq2mn2  ! ========= q2jmin loop =========

          q2jmin=q2mnval(nq2mn)
          
          call setBasicq(q2jmin,ir)
          
          baspar(nq2mn,  ir)=alpff(ir)
          baspar(nq2mn,2+ir)=betff(ir)
          baspar(nq2mn,4+ir)=betpom(ir)
          baspar(nq2mn,6+ir)=dels(ir)
          baspar(nq2mn,8+ir)=gamsoft(ir)
          do k=1,nclha
          baspar(nq2mn,10+(ir-1)*nclha+k)=sngl(coefxu1(ir,k))
          enddo
      
        enddo                   ! ========= q2jmin loop =========
      enddo                     ! ======= proj/tar loop ========

      q2imin(1)=0.
      q2imin(2)=0.
      call ipoBasics(q2imin)
      q2cmin(1)=q2nmin
      q2cmin(2)=q2nmin

      if(abs(alpqua(1)).lt.1.e-6)call utstop('alpq should not be 0 !&')
      if(abs(alpqua(2)).lt.1.e-6)call utstop('alpq should not be 0 !&')
      alpr=-2.+alpqua(1)      !x-exponent for remnant mass (in DIS ... to be checked !)
      gamhad(iclpro)=alpff(1)
      gamhad(icltar)=alpff(2)

c for j/psi compute val pdf at q2nmin (fixed)      
      if(iclpro.eq.4)then
        gnorm=0.
        do i=1,7
          do m=1,2
            xx=.5+x1(i)*(m-1.5)
            xxg=xx**(1./(1.-dels(1)))
            gnorm=gnorm+a1(i)*(fzeroGluZZ(dble(xxg),4)
     &           +FzeroQuaZZ(dble(xxg),4))
          enddo
        enddo
        gnorm=gnorm/(1.-dels(1))*0.5 !*4.*pi*gamhad(4)*ffrr
        alvc=6./(1.-gnorm)-4.
      else
        alvc=0.
      endif
      if(ish.ge.4)write(ifch,*)'alvc',alvc

      endif
      
      end

c----------------------------------------------------------------------
      subroutine setBasicq(q2jmin,ir)
c----------------------------------------------------------------------
c ir=1 : projectile
c ir=2 : target
c----------------------------------------------------------------------
#include "aaa.h"
#include "par.h"
#include "sem.h"
      common /psar36/ alvc,qnorm,qnormp
      common /psar40/ coefxu1(2,nclha)
      double precision coefxu1,utgam2
c      common /psar36/ alvc
      common /ar3/    x1(7),a1(7)
c      double precision fzeroGluZZ,FzeroQuaZZ
 
      if(ir.eq.1)then
        icl=iclpro
        if(icl.eq.4)icl=1
c        qnorma=qnorm
      else
        icl=icltar
c        qnorma=qnormp
      endif


c      return

      dels(ir)=alppom-1.
      betpom(ir)=betpomi
      
      !correction only at very low energy (< q2jmin -> only for F2 at low Q2)
c      sat=max(0.,1.+zoeinc*max(-1.,log10(min(q2sft,q2jmin)/q2sft)))
c      dels(ir)=min(0.99,(alppom-1.)*sat**0.25)   !doesn't change pt spectra much
c
cctp      dels=(alppom-1.)*(1.-log(min(1.,max(q2ini,q2nmin)/q2jmin))
cctp     &                   /log(q2ini/q2jmin)) !alppom-1.
c
cc      sat=1.-min(1.,(q2nmin/q2jmin))
c      betpom(ir)=betpomi*sat**2    !larger betpom=softer pt spectra



      gamsoft(ir)=gamtil
     *       *utgam1(2.+betpom(ir)-dels(ir))/utgam1(1.-dels(ir))
     *                                  /utgam1(1.+betpom(ir))
c     *       *min(1.,q2nmin/q2jmin)
c     *       *exp(2.*(dels/(alppom-1.)-1.))      !correction for low energy

c      if(zoeinc.gt.0.)gamsoft(ir)=gamsoft(ir)/sat
c     *     **max(0.,1.-zoeinc*0.6) !change normalization of jet pt spectrum
      

      if(ish.ge.7)write(ifch,*)'ir,qq,betpom,gamsoft'
     *           ,ir,q2jmin,betpom(ir),gamsoft(ir)

c      betff(ir)=-(alpparh+(alppar-alpparh)/max(1.,sat)**0.5)
      betff(ir)=-alppar

      alpqua(ir)=(1.-betff(ir))/2.

c      print *,'ptb',engy,q2jmin,q2nmin,betpom(ir),dels(ir),gamsoft(ir)
c     &             ,betff(ir),sat

!     normalization to get psdpdf in EsoftValTil
      do k=1,nclha
        coefxu1(ir,k)=utgam2(2.d0+dble(alplea(k)-alpqua(ir)))  
     *  /utgam2(1.d0+dble(alplea(k)))/utgam2(1.d0-dble(alpqua(ir)))
c     *  *dble(2.*alpqua(1))
      enddo

      if(iappl.eq.1.or.iappl.eq.2.or.iappl.eq.8.or.iappl.eq.9)then

        !parton density normalization

c        if(q2jmin.gt.q2nmin)then
          
          if(icl.eq.2)then

            sq=log(log(q2jmin/.232**2)/log(.23/.232**2))
            du=2.997+.753*sq-.076*sq*sq
            qnormj=0.
            do i=1,7
              do m=1,2
                xx=.5+x1(i)*(m-1.5)
                xxq=1.-xx**(1./(1.+du))
                qnormj=qnormj+a1(i)*(psdpdf(xxq,q2jmin,0.,2,1)
     *               +psdpdf(xxq,q2jmin,0.,2,2))/(1.-xxq)**du
              enddo
            enddo
            qnormj=qnormj*.5/(1.+du)
            
            if(ish.ge.7)write (ifch,*)'qnormj baryon',qnormj
            
          else
            
            sq=log(log(q2jmin/.232**2)/log(.25/.232**2))
            dpi=.367+.563*sq
            qnormj=0.
            do i=1,7
              do m=1,2
                xx=.5+x1(i)*(m-1.5)
                xxq=1.-xx**(1./(1.+dpi))
                qnormj=qnormj+a1(i)*(psdfh4(xxq,q2jmin,0.,icl,1)
     *               +psdfh4(xxq,q2jmin,0.,icl,2))/(1.-xxq)**dpi
              enddo
            enddo
            qnormj=qnormj*.5/(1.+dpi)
            
            if(ish.ge.7) write (ifch,*)'qnorm meson',qnormj

          endif
          
c        else
c          qnormj=qnorma
c        endif

        alpff(ir)=(1.-qnormj)/utgam1(2.+betff(ir))
     *  /utgam1(alplea(icl)+1.)*utgam1(alplea(icl)+3.+betff(ir))
        
        !correction factor
        if(gamhadsi(icl).gt.0.)then
          alpff(ir)=alpff(ir)*gamhadsi(icl)
        endif

        gamhad(icl)=alpff(ir)

      endif
      
      if(ish.ge.7.and.betff(ir).lt.0.)
     .write(ifch,*)'alpff,betff',alpff(ir),betff(ir)
    

      
      end


c-----------------------------------------------------------------------
      subroutine mkParamTables
c-----------------------------------------------------------------------
c     make tables alpD,alpDp,alpDpp,betD,betDp,betDpp,gamD,delD
c     and inelastic cross-section
c-----------------------------------------------------------------------
#include "aaa.h"
#include "sem.h"
#include "tab.h"
#include "ems.h"
#include "par.h"
      common/producetab/ producetables              !used to link with CRMC
      logical producetables
      common/geom/rmproj,rmtarg,bmax,bkmx
      character textini*38
      external ptfau,ptfauAA
      parameter(nbkbin=100)
      common /kfitd/ xkappafit(nbkbin,nclegy,nclha,nclha),xkappa,bkbin
      dimension gamhad0(nclha),r2had0(nclha)
     *,alplea0(nclha),asect11(7,4,8),asect13(7,4,8),asect21(7,4,8)
     *,asect23(7,4,8),asect31(7,8,8),asect33(7,8,8)
     *,asect41(7,8,8),asect43(7,8,8)!,cgam(idxD)
      logical lcalc


      call utpri('mkPara',ish,ishini,4)


c--------------------------------------
c fit parameters and inelastic cross sections
c---------------------------------------

      if(isetcs.ge.2)then !--------------------

      if(ish.ge.4)write(ifch,*)'cross sections ...'
 6    continue
      inquire(file=fncs,exist=lcalc)
      if(lcalc)then
        write(ifmt,'(3a)')'read from ',fncs(1:nfncs),' ...'
        open(1,file=fncs(1:nfncs),status='old')
        read (1,*)alppar0,alplea0,alppom0,slopom0,
     *  gamhad0,r2had0,gampar0,
     *  qcdlam0,q2nmin0,q2ini0,betpom0,glusea0,naflav0,
     *  factk0,pt2cut0
        if(alppar0.ne.alppar)write(ifmt,'(a,2f8.4)')
     *  'cs: wrong alpqua',alppar0,alppar
        if(alppom0.ne.alppom)write(ifmt,'(a,2f8.4)')
     *  'cs: wrong alppom',alppom0,alppom
        if(slopom0.ne.slopom)write(ifmt,'(a,2f8.4)')
     *  'cs: wrong slopom',slopom0,slopom
        iii=2
c        if(gamhad0(iii).ne.gamhad(iii))write(ifmt,'(a,i1,a,2f8.4)')
c     *  'cs: wrong gamhad(',iii,')',gamhad0(iii),gamhad(iii)
        do iii=1,3
        if(r2had0(iii) .ne.r2had(iii) )write(ifmt,'(a,i1,a,2f8.4)')
     *  'cs: wrong r2had(',iii,')',r2had0(iii),r2had(iii)
        if(alplea0(iii).ne.alplea0(iii))write(ifmt,'(a,i1,a,2f8.4)')
     *  'cs: wrong alplea(',iii,')',alplea0(iii),alplea(iii)
        enddo
        if(gampar0.ne.gampar)write(ifmt,'(a,2f8.4)')
     *  'cs: wrong gampar',gampar0,gampar
        if(qcdlam0.ne.qcdlam)write(ifmt,'(a,2f8.4)')
     *  'cs: wrong qcdlam',qcdlam0,qcdlam
        if(q2nmin0 .ne.q2nmin )write(ifmt,'(a,2f8.4)')
     *  'cs: wrong q2nmin',q2nmin0,q2nmin
        if(q2ini0 .ne.q2ini )write(ifmt,'(a,2f8.4)')
     *  'cs: wrong q2ini',q2ini0,q2ini
        if(betpom0.ne.betpomi)write(ifmt,'(a,2f8.4)')
     *  'cs: wrong betpom',betpom0,betpomi
        if(glusea0.ne.glusea)write(ifmt,'(a,2f8.4)')
     *  'cs: wrong glusea',glusea0,glusea
        if(naflav0.ne.naflav)write(ifmt,'(a,2f8.4)')
     *  'cs: wrong naflav',naflav0,naflav
        if(factk0 .ne.factk )write(ifmt,'(a,2f8.4)')
     *  'cs: wrong factk', factk0,factk
        if(pt2cut0 .ne.pt2cut )write(ifmt,'(a,2f8.4)')
     *  'cs: wrong pt2cut', pt2cut0,pt2cut
        if(alppar0.ne.alppar.or.alppom0.ne.alppom
     *  .or.slopom0.ne.slopom!.or.gamhad0(2).ne.gamhad(2)
     *  .or.r2had0(1).ne.r2had(1).or.r2had0(2).ne.r2had(2)
     *  .or.r2had0(3).ne.r2had(3)
     *  .or.alplea0(1).ne.alplea(1).or.alplea0(2).ne.alplea(2)
     *  .or.alplea0(3).ne.alplea(3)
     *  .or.qcdlam0.ne.qcdlam.or.q2nmin0 .ne.q2nmin
     *  .or.q2ini0 .ne.q2ini
     *  .or.betpom0.ne.betpomi.or.glusea0.ne.glusea.or.naflav0.ne.naflav
     *  .or.factk0 .ne.factk .or.pt2cut0.ne.pt2cut)then
           write(ifmt,'(//a//)')'   cs has to be reinitialized!!!!'
           stop
        endif

        read(1,*)bkbin0,iclpro10,iclpro20,icltar10,icltar20,iclegy10
     *,iclegy20,egylow0,egymax0,iomega0,egyscr0,epscrw0,epscrp0,isetcs0

        if(bkbin0.ne.bkbin)write(ifmt,'(a,2f8.4)')
     *  'cs: wrong iclpro1',bkbin0,bkbin
        if(iclpro10.ne.iclpro1)write(ifmt,'(a,2i2)')
     *  'cs: wrong iclpro1',iclpro10,iclpro1
        if(iclpro20.ne.iclpro2)write(ifmt,'(a,2i2)')
     *  'cs: wrong iclpro2',iclpro20,iclpro2
        if(icltar10.ne.icltar1)write(ifmt,'(a,2i2)')
     *  'cs: wrong icltar1',icltar10,icltar1
        if(icltar20.ne.icltar2)write(ifmt,'(a,2i2)')
     *  'cs: wrong icltar2',icltar20,icltar2
        if(iclegy10.ne.iclegy1)write(ifmt,'(a,2i4)')
     *  'cs: wrong iclegy1',iclegy10,iclegy1
        if(iclegy20.ne.iclegy2)write(ifmt,'(a,2i4)')
     *  'cs: wrong iclegy2',iclegy20,iclegy2
        if(iomega0.ne.iomega)write(textini,'(a,2i8)')
     *  'cs: wrong iomega ',iomega0,iomega
        if(egylow0.ne.egylow)write(ifmt,'(a,2f8.4)')
     *  'cs: wrong egylow',egylow0,egylow
        if(egymax0.ne.egymax)write(ifmt,'(a,2f12.4)')
     *  'cs: wrong egymax',egymax0,egymax
        if(egyscr0.ne.egyscr)write(ifmt,'(a,2f8.4)')
     *  'cs: wrong egyscr ',egyscr0,egyscr
        if(epscrw0.ne.epscrw)write(ifmt,'(a,2f8.4)')
     *  'cs: wrong epscrw',epscrw0,epscrw
        if(epscrp0.ne.epscrp)write(ifmt,'(a,2f8.4)')
     *  'cs: wrong epscrp',epscrp0,epscrp
        if(isetcs0.lt.isetcs)write(ifmt,'(a,2f8.4)')
     *  'cs: wrong isetcs',isetcs0,isetcs
        if(iclpro10.ne.iclpro1.or.iclpro20.ne.iclpro2
     *   .or.icltar10.ne.icltar1.or.icltar20.ne.icltar2
     *   .or.iclegy10.ne.iclegy1.or.iclegy20.ne.iclegy2
     *   .or.egylow0.ne.egylow.or.egymax0.ne.egymax
     *   .or.egyscr0.ne.egyscr.or.epscrw0.ne.epscrw.or.isetcs0.lt.isetcs
     *   .or.epscrp0.ne.epscrp.or.bkbin0.ne.bkbin)then
           write(ifmt,'(//a//)')'   cs has to be reinitialized!!!!'
           stop
        endif
c        read(1,*)xkappafit,chadrs
c     *          ,alpDs,alpDps,alpDpps,betDs,betDps,betDpps,gamDs,delDs
        read(1,*)((((xkappafit(iiib,iiiegy,iiipro,iiitar)
     *       ,iiib=1,nbkbin),iiiegy=iclegy1,iclegy2)
     *       ,iiipro=iclpro1,iclpro2),iiitar=icltar1,icltar2)        
        read(1,*)((((chadrs(iii,iiiegy,iiipro,iiitar)
     *       ,iii=1,4),iiiegy=iclegy1,iclegy2),iiipro=iclpro1,iclpro2)
     *       ,iiitar=icltar1,icltar2)
        read(1,*)((((alpDs(iiidf,iiiegy,iiipro,iiitar),
     *   alpDps(iiidf,iiiegy,iiipro,iiitar),
     *   alpDpps(iiidf,iiiegy,iiipro,iiitar),
     *   betDs(iiidf,iiiegy,iiipro,iiitar),
     *   betDps(iiidf,iiiegy,iiipro,iiitar),
     *   betDpps(iiidf,iiiegy,iiipro,iiitar),
     *   gamDs(iiidf,iiiegy,iiipro,iiitar),
     *   delDs(iiidf,iiiegy,iiipro,iiitar)
     *  ,iiidf=idxD0,idxD1),iiiegy=iclegy1,iclegy2)
     *  ,iiipro=iclpro1,iclpro2),iiitar=icltar1,icltar2)
c        do iiipro=iclpro1,iclpro2
c        do iiitar=icltar1,icltar2
c        do iiiegy=iclegy1,iclegy2
c        do iiib=1,nbkbin
c          read(1,*)xkappafit(iiiegy,iiipro,iiitar,iiib)
c        enddo
cc        xkappafit(iiiegy,iiipro,iiitar,nbkbin)=1.
cc        do iiib=2,nbkbin-1
cc          if(xkappafit(iiiegy,iiipro,iiitar,iiib).lt.1.)then
cc            xkappafit(iiiegy,iiipro,iiitar,iiib)=max(1.,0.5*
cc     *        (xkappafit(iiiegy,iiipro,iiitar,iiib-1)
cc     *        +xkappafit(iiiegy,iiipro,iiitar,iiib+1)))
cc          endif
cc        enddo
c        do iiidf=idxD0,idxD
c         read(1,*)alpDs(iiidf,iiiegy,iiipro,iiitar),
c     *   alpDps(iiidf,iiiegy,iiipro,iiitar),
c     *   alpDpps(iiidf,iiiegy,iiipro,iiitar),
c     *   betDs(iiidf,iiiegy,iiipro,iiitar),
c     *   betDps(iiidf,iiiegy,iiipro,iiitar),
c     *   betDpps(iiidf,iiiegy,iiipro,iiitar),
c     *   gamDs(iiidf,iiiegy,iiipro,iiitar),
c     *   delDs(iiidf,iiiegy,iiipro,iiitar)
c        enddo
c        enddo
c        enddo
c        enddo
        if(isetcs.eq.2)then
          if(ionudi.eq.1)then
            read (1,*)asect,asect13,asect21,asect23,asectn
     *               ,asect33,asect41,asect43
          else  !ionudi=3
            read (1,*)asect11,asect,asect21,asect23,asect31
     *               ,asectn,asect41,asect43
          endif
        elseif(isetcs.eq.3)then
          if(ionudi.eq.1)then
            read (1,*)asect11,asect13,asect,asect23,asect31
     *               ,asect33,asectn,asect43
          else  !ionudi=3
            read (1,*)asect11,asect13,asect21,asect,asect31
     *               ,asect33,asect41,asectn
          endif
        else
           write(ifmt,'(//a//)')' Wrong isetcs in psaini !!!!'
        endif

        close(1)

        goto 7


      elseif(.not.producetables)then
        write(ifmt,*) "Missing epos cs file !"        
        write(ifmt,*) "Please correct the defined path ",
     &"or force production ..."
        stop

      elseif(maxq2mx.ne.maxq2mn)then

        write(ifmt,*) "Missing epos cs file !"        
        write(ifmt,*) "Energy too low to produce it properly !"        
        stop

      endif

      ifradesave=ifrade
      iremnsave=iremn
      idprojsave=idproj
      idprojinsave=idprojin
      idtargsave=idtarg
      idtarginsave=idtargin
      laprojsave=laproj
      latargsave=latarg
      maprojsave=maproj
      matargsave=matarg
      icltarsave=icltar
      iclprosave=iclpro
      engysave=engy
      pnllsave=pnll
      elabsave=elab
      ecmssave=ecms
      iclegysave=iclegy
      nrevtsave=nrevt
      neventsave=nevent
      ntevtsave=ntevt
      isetcssave=isetcs
      noebinsave=noebin
      isigmasave=isigma
      bminimsave=bminim
      bmaximsave=bmaxim
      bimevtsave=bimevt
      bkmxndifsave=bkmxndif
c      fctrmxsave=fctrmx
      ionudisave=ionudi


      isetcs=2
      isigma=1
      noebin=1
      idtarg=1120
      idtargin=1120
      bminim=0.
      bmaxim=10000.
      ifrade=0            !to save time, no fragmentation
      iremn=0             !to save time, simple remnants
      ionudi=3            !to have both ionudi=1 and 3 in tables

      write(ifmt,'(a,6i3)')'iclpro,icltar,iclegy :'
     .                ,iclpro1,iclpro2,icltar1,icltar2,iclegy1,iclegy2
      write(ifmt,'(a,$)')'cs does not exist -> tabulate fit param '
      
      q2cmin(1)=q2nmin
      q2cmin(2)=q2nmin
      q2pmin(1)=q2cmin(1)
      q2pmin(2)=q2cmin(2)
      call ipoOm5Tables(1)      !used for psvin and thus for om51p
      do iclpro=iclpro1,iclpro2 !hadron type (1 - pion, 2 - nucleon, 3 - kaon, 4 - charm)
      do icltar=icltar1,icltar2 !hadron type (2 - nucleon)
      do iclegy=iclegy2,iclegy1,-1
        write(ifmt,'(a,i2,a,$)')'.(',iclegy,')'
        engy=egyfac**(iclegy-1)*egylow
        call paramini(-1)
      enddo
      enddo
      enddo
      write(ifmt,'(a)')'.'

      write(ifmt,'(a)')'cs does not exist -> calculate tables  ...'

c initialize random numbers
      if(seedj.ne.0d0)then
        call ranfini(seedj,iseqsim,2)
      else
        stop 'seedi = 0 ... Please define it !'
      endif
      call aseed(2)

      laproj=-1
      maproj=1
      icltar=2
      do iclpro=1,4
       if(iclpro.lt.iclpro1.or.iclpro.gt.iclpro2)then
         do ie=1,7
           do iia=1,8
             asect11(ie,iclpro,iia)=0.
             asect21(ie,iclpro,iia)=0.
             asect13(ie,iclpro,iia)=0.
             asect23(ie,iclpro,iia)=0.
           enddo
         enddo
       else
         do ie=1,7
           engy=1.5*10.**(ie-1)
           call paramini(0)
c           bkmxndif=conbmxndif()
           if(ish.ge.1)
     &     write(ifch,*)'  calcul.   ',ie,'  (',iclpro,')',engy
           write(ifmt,*)'  calcul.   ',ie,'  (',iclpro,')',engy

           sigine=0.
           do iia=1,8
            matarg=2**(iia-1)
            if(matarg.eq.1)then !hadron-proton interaction
c ine=cut+diff
              call psfz(2,gz2,0.)
              gin=gz2*pi*10.
c cut
              iomegasave=iomega
              iomega=2
              call psfz(2,gz2,0.)
              iomega=iomegasave
              gcut=gz2*pi*10.
c diff
              difpart=gin-gcut
c  non excited projectile and target
              gqela=(1.-rexdif(iclpro))*(1.-rexdif(icltar))*difpart
              gin3=max(1.,gin-gqela)              
            else
              call conini
              rad=radnuc(matarg)
              bm=rad+2.
              rrr=rad/difnuc(matarg)
              rrrm=rrr+log(9.)
              anorm=1.5/pi/rrr**3/(1.+(pi/rrr)**2)/difnuc(matarg)**2
c             gela=(ptgau(ptfau,bm,2,1)+ptgau1(bm,2,1))*10. !sig_ela
c in=cut+diff
              gcut=(ptgau(ptfau,bm,2,2)+ptgau1(bm,2,2))*10. !sig_in
              gin=gcut
c cut
              iomegasave=iomega
              iomega=2
              gcut=(ptgau(ptfau,bm,2,2)+ptgau1(bm,2,2))*10. !sig_cut
              iomega=iomegasave
c diff
              difpart=gin-gcut
c  non excited projectile
              gqela=(1.-rexdif(iclpro))
     &             **(1.+rexres(iclpro)*float(matarg-1)**0.3)
c  non excited target
              gqela=gqela*(1.-rexdif(icltar))
              gqela=gqela*difpart
              gin3=max(1.,gin-gqela)
            endif
            if(ish.ge.1)write (ifch,226)matarg,gin,gin3
226         format(2x,'psaini: hadron-nucleus (',i3,') cross sections:'/
     *       4x,'gin,gin3=',2e10.3)
            write(ifmt,*)'  matarg,gin,gin3:',matarg,gin,gin3
            asect11(ie,iclpro,iia)=log(gin)
            asect13(ie,iclpro,iia)=log(gin3)
           enddo
         enddo

         if(isetcssave.ge.3)then

         if(iclpro.eq.1)then
          idprojin=120
         elseif(iclpro.eq.2)then
          idprojin=1120
         elseif(iclpro.eq.3)then
          idprojin=130
         endif
         do ie=1,7
          engy=1.5*10.**(ie-1)
           if(engy.le.egymin)engy=egymin
           if(engy.ge.egymax)engy=egymax
           write(ifmt,*)'  simul.   ',ie,'  (',iclpro,')',engy
           if(ish.ge.1)
     &     write(ifch,*)'  simul.   ',ie,'  (',iclpro,')',engy
           do iia=1,8
            matarg=2**(iia-1)
            latarg=min(1,matarg/2)
c            fctrmx=max(ftcrmxsave,float(matarg))          !to get stable pA and AA cross section, this number has to be large for large A
            ntevt=0
            nrevt=0
            pnll=-1.
            elab=-1.
            ecms=-1.
            ekin=-1.
            call conini
            call ainit
            nevent=50000
            if(iia.gt.40)nevent=nevent/10
c            if(matarg.eq.1)nevent=1
            call epocrossc(nevent,sigt,sigi,sigc,sige,sigql,sigd)
c do not count non-excited diffractive projectile in inelastic
            sigi3=sigi-sigql
            if(ish.ge.1)write (ifch,228)matarg,sigi,sigi3
 228        format(2x,'simul.: hadron-nucleus (',i3,') cross sections:'/
     *       4x,'gin,gin3=',2e10.3)
            write(ifmt,*)'  matarg,sigi,sigi3 :',matarg,sigi,sigi3
            asect21(ie,iclpro,iia)=log(sigi)
            asect23(ie,iclpro,iia)=log(sigi3)
c            do  n=1,nevent
c              ntry=0
c 222          ntevt=ntevt+1
c              iret=0
c              ntry=ntry+1
c              bimevt=-1.
c              if(ntry.lt.10000)then
cc if random sign for projectile, set it here
c                idproj=idprojin*(1-2*int(rangen()+0.5d0))
c                call emsaaa(iret)
c                if(iret.gt.0)goto 222
c              else
c                ntevt=ntry
c              endif
c            enddo
c            a=pi*bmax**2
c            if(a.gt.0..and.ntevt.gt.0.)then
c             xs=anintine/float(ntevt)*a*10.
c             write(ifmt,*)'  matarg,nevent,ntevt,bmax,xs :'
c     .       ,matarg,anintine,ntevt,bmax,xs
c             write(ifch,*)'  matarg,nevent,ntevt,bmax,xs :'
c     .       ,matarg,anintine,ntevt,bmax,xs
c             asect2(ie,iclpro,iia)=log(xs)
c            else
c             write(ifmt,*)' Problem ? ',iclpro,matarg,bmax,ntevt
c             asect2(ie,iclpro,iia)=0.
c            endif
          enddo
        enddo
        else
          do ie=1,7
            do iia=1,8
              asect21(ie,iclpro,iia)=0.
              asect23(ie,iclpro,iia)=0.
            enddo
          enddo
        endif
       endif
      enddo

      idprojin=1120
      iclpro=2
      icltar=2
      do ie=1,7
        engy=1.5*10.**(ie-1)
        call paramini(0)
c        bkmxndif=conbmxndif()
        if(ish.ge.1)
     &  write(ifch,*)'  calcul. AB  ',ie,engy
        write(ifmt,*)'  calcul. AB  ',ie,engy

        do iia=1,8
          maproj=2**(iia-1)
          laproj=max(1,maproj/2)
        do iib=1,8
          matarg=2**(iib-1)
          latarg=max(1,matarg/2)
          sigine=0.
          if(matarg.eq.1.and.maproj.eq.1)then !proton-proton interaction
c ine=cut+diff
            call psfz(2,gz2,0.)
            gin=gz2*pi*10.
c cut
            iomegasave=iomega
            iomega=2
            call psfz(2,gz2,0.)
            iomega=iomegasave
            gcut=gz2*pi*10.
c diff
            difpart=gin-gcut
c  non excited projectile and target
            gqela=(1.-rexdif(iclpro))*(1.-rexdif(icltar))*difpart
            gin3=max(1.,gin-gqela)              
          else
            call conini
            if(maproj.eq.1)then
              rad=radnuc(matarg)
              bm=rad+2.
              rrr=rad/difnuc(matarg)
              rrrm=rrr+log(9.)
              anorm=1.5/pi/rrr**3/(1.+(pi/rrr)**2)/difnuc(matarg)**2
c              gela=(ptgau(ptfau,bm,2,1)+ptgau1(bm,2,1))*10. !sig_ela
c in=cut+diff
              gcut=(ptgau(ptfau,bm,2,2)+ptgau1(bm,2,2))*10. !sig_in
              gin=gcut
c cut
              iomegasave=iomega
              iomega=2
              gcut=(ptgau(ptfau,bm,2,2)+ptgau1(bm,2,2))*10. !sig_cut
              iomega=iomegasave
c diff
              difpart=gin-gcut
c  non excited projectile
              gqela=(1.-rexdif(iclpro))
     &             **(1.+rexres(iclpro)*float(matarg-1)**0.3)
c  non excited target
              gqela=gqela*(1.-rexdif(icltar))**(1.+float(matarg)**0.3)
              gqela=gqela*difpart
              gin3=max(1.,gin-gqela)
            elseif(matarg.eq.1)then
              radp=radnuc(maproj)
              bm=radp+2.
              rrrp=radp/difnuc(maproj)
              rrrmp=rrrp+log(9.)
              anormp=1.5/pi/rrrp**3/(1.+(pi/rrrp)**2)/difnuc(maproj)**2
c              gtot=(ptgau(ptfau,bm,1,1)+ptgau1(bm,1,1))*10. !sig_in
c in=cut+diff
              gcut=(ptgau(ptfau,bm,1,2)+ptgau1(bm,1,2))*10. !sig_in
              gin=gcut     !in=cut+diff
c cut
              iomegasave=iomega
              iomega=2
              gcut=(ptgau(ptfau,bm,1,2)+ptgau1(bm,1,2))*10. !sig_cut
              iomega=iomegasave
c diff
              difpart=gin-gcut
c  non excited projectile
              gqela=(1.-rexdif(iclpro))**(1.+float(maproj)**0.3)
c  non excited target
              gqela=gqela*(1.-rexdif(icltar))
     &             **(1.+rexres(icltar)*float(maproj-1)**0.3)
              gqela=gqela*difpart
              gin3=max(1.,gin-gqela)
            else
              rad=radnuc(matarg)+1.
              radp=radnuc(maproj)+1.
              bm=rad+radp+2.
              rrr=rad/difnuc(matarg)
              rrrm=rrr+log(9.)
              rrrp=radp/difnuc(maproj)
              rrrmp=rrrp+log(9.)
              anorm=1.5/pi/rrr**3/(1.+(pi/rrr)**2)/difnuc(matarg)**2
              anormp=1.5/pi/rrrp**3/(1.+(pi/rrrp)**2)/difnuc(maproj)**2
c ine=cut+diff
c              gtot=(ptgau(ptfauAA,bm,2,1)+ptgau2(bm,1))*10.
              gcut=(ptgau(ptfauAA,bm,2,2)+ptgau2(bm,2))*10.
c              gin=gtot
              gin=gcut
c cut
              iomegasave=iomega
              iomega=2
              gcut=(ptgau(ptfauAA,bm,2,2)+ptgau2(bm,2))*10. !sig_cut
              iomega=iomegasave
c diff
              difpart=gin-gcut
c  non excited projectile
              gqelap=(1.-rexdif(iclpro))
     &             **(1.+rexres(iclpro)*float(matarg-1)**0.3)
              gqelap=gqelap**(1.+float(maproj)**0.3)
c  non excited target
              gqelat=(1.-rexdif(icltar))
     &             **(1.+rexres(icltar)*float(maproj-1)**0.3)
              gqelat=gqelat**(1.+float(maproj)**0.3)
              gqela=gqelap*gqelat*difpart
              gin3=gin-gqela
            endif
          endif
          if(ish.ge.1)write (ifch,227)maproj,matarg,gin,gin3
 227      format(2x,'psaini: nucleus-nucleus (',i3,'-',i3
     *       ,') cross sections:',/,4x,'gin,gin3=',2e10.3)
            write(ifmt,*)'  maproj,matarg,gin,gin3 :'
     *       ,maproj,matarg,gin,gin3
            asect31(ie,iia,iib)=log(gin)
            asect33(ie,iia,iib)=log(gin3)

          enddo
        enddo
      enddo

      if(isetcssave.ge.3)then

      do ie=1,7
        engy=1.5*10.**(ie-1)
        if(engy.le.egymin)engy=egymin
        if(engy.ge.egymax)engy=egymax
        write(ifmt,*)'  AB xs   ',ie,engy
        if(ish.ge.1)
     &  write(ifch,*)'  AB xs   ',ie,engy
        do iia=1,8
          maproj=2**(iia-1)
          laproj=max(1,maproj/2)
        do iib=1,8
          matarg=2**(iib-1)
          latarg=max(1,matarg/2)
c          fctrmx=max(ftcrmxsave,float(max(maproj,matarg))) !to get stable pA and AA cross section, this number has to be large for large A
          ntevt=0
          nrevt=0
          pnll=-1.
          elab=-1.
          ecms=-1.
          ekin=-1.
          call conini
          call ainit
          nevent=10000
          if(iia.gt.40)nevent=nevent/5
          if(iib.gt.40)nevent=nevent/5
c          if(maproj+matarg.eq.2)nevent=1
          call epocrossc(nevent,sigt,sigi,sigc,sige,sigql,sigd)
c do not count non-excited diffractive projectile in inelastic
          sigi3=sigi-sigql
          if(ish.ge.1)write (ifch,229)maproj,matarg,sigi,sigi3
 229      format(2x,'simul.: nucleus-nucleus (',i3,'-',i3
     *       ,') cross sections:',/,4x,'gin,gin3=',2e10.3)
         write(ifmt,*)'  maproj,matarg,sigi,sigi3 :',maproj,matarg
     &                                               ,sigi,sigi3
          asect41(ie,iia,iib)=log(sigi)
          asect43(ie,iia,iib)=log(sigi3)

c          do  n=1,nevent
c            ntry=0
c 223        ntevt=ntevt+1
c            iret=0
c            ntry=ntry+1
c            bimevt=-1.
c            if(ntry.lt.10000)then
c              call emsaaa(iret)
c              if(iret.gt.0)goto 223
c            else
c              ntevt=ntry
c            endif
c          enddo
c          a=pi*bmax**2
c          if(a.gt.0..and.ntevt.gt.0.)then
c            xs=anintine/float(ntevt)*a*10.
c          write(ifmt,*)'  maproj,matarg,nevent,ntevt,bmax,xs :'
c     &                         ,maproj,matarg,anintine,ntevt,bmax,xs
c          write(ifch,*)'  maproj,matarg,nevent,ntevt,bmax,xs :'
c     &                         ,maproj,matarg,anintine,ntevt,bmax,xs
c            asect4(ie,iia,iib)=log(xs)
c          else
c            write(ifmt,*)' Problem ? ',maproj,matarg,bmax,ntevt
c            asect4(ie,iia,iib)=0.
c          endif
        enddo
      enddo
      enddo
      else
        do ie=1,7
          do iia=1,8
            do iib=1,8
              asect41(ie,iia,iib)=0.
              asect43(ie,iia,iib)=0.
            enddo
          enddo
        enddo
      endif

      ifrade=ifradesave
      iremn=iremnsave
      idproj=idprojsave
      idprojin=idprojinsave
      idtarg=idtargsave
      idtargin=idtarginsave
      laproj=laprojsave
      latarg=latargsave
      maproj=maprojsave
      matarg=matargsave
      icltar=icltarsave
      iclpro=iclprosave
      engy=engysave
      pnll=pnllsave
      elab=elabsave
      ecms=ecmssave
      iclegy=iclegysave
      nrevt=nrevtsave
      nevent=neventsave
      ntevt=ntevtsave
      isetcs=isetcssave
      noebin=noebinsave
      isigma=isigmasave
      bminim=bminimsave
      bmaxim=bmaximsave
      bimevt=bimevtsave
      bkmxndif=bkmxndifsave
      ionudi=ionudisave
c      fctrmx=fctrmxsave

      write(ifmt,'(a)')'write to cs ...'
      open(1,file=fncs,status='unknown')
      write (1,*)alppar,alplea,alppom,slopom,gamhad,r2had,gampar,
     *qcdlam,q2nmin,q2ini,betpomi,glusea,naflav,factk,pt2cut
      write(1,*)bkbin,iclpro1,iclpro2,icltar1,icltar2,iclegy1,iclegy2
     *,egylow,egymax,iomega,egyscr,epscrw,epscrp,isetcs
      write(1,*)((((xkappafit(iiib,iiiegy,iiipro,iiitar)
     *       ,iiib=1,nbkbin),iiiegy=iclegy1,iclegy2)
     *       ,iiipro=iclpro1,iclpro2),iiitar=icltar1,icltar2)
      write(1,*)((((chadrs(iii,iiiegy,iiipro,iiitar)
     *       ,iii=1,4),iiiegy=iclegy1,iclegy2),iiipro=iclpro1,iclpro2)
     *       ,iiitar=icltar1,icltar2)
      write(1,*)((((alpDs(iiidf,iiiegy,iiipro,iiitar),
     *   alpDps(iiidf,iiiegy,iiipro,iiitar),
     *   alpDpps(iiidf,iiiegy,iiipro,iiitar),
     *   betDs(iiidf,iiiegy,iiipro,iiitar),
     *   betDps(iiidf,iiiegy,iiipro,iiitar),
     *   betDpps(iiidf,iiiegy,iiipro,iiitar),
     *   gamDs(iiidf,iiiegy,iiipro,iiitar),
     *   delDs(iiidf,iiiegy,iiipro,iiitar)
     *  ,iiidf=idxD0,idxD1),iiiegy=iclegy1,iclegy2)
     *  ,iiipro=iclpro1,iclpro2),iiitar=icltar1,icltar2)
c      write(1,*)xkappafit,chadrs
c     *         ,alpDs,alpDps,alpDpps,betDs,betDps,betDpps,gamDs,delDs
c      do iiipro=iclpro1,iclpro2
c       do iiitar=icltar1,icltar2
c        do iiiegy=iclegy1,iclegy2
c        do iiib=1,nbkbin
c          write(1,*)xkappafit(iiiegy,iiipro,iiitar,iiib)
c        enddo
c        do iiidf=idxD0,idxD
c         write(1,*)alpDs(iiidf,iiiegy,iiipro,iiitar),
c     *   alpDps(iiidf,iiiegy,iiipro,iiitar),
c     *   alpDpps(iiidf,iiiegy,iiipro,iiitar),
c     *   betDs(iiidf,iiiegy,iiipro,iiitar),
c     *   betDps(iiidf,iiiegy,iiipro,iiitar),
c     *   betDpps(iiidf,iiiegy,iiipro,iiitar),
c     *   gamDs(iiidf,iiiegy,iiipro,iiitar),
c     *   delDs(iiidf,iiiegy,iiipro,iiitar)
c        enddo
c        enddo
c       enddo
c      enddo
      write (1,*)asect11,asect13,asect21,asect23
     *          ,asect31,asect33,asect41,asect43

      close(1)


      goto 6

 7    continue

      endif !----------isetcs.ge.2-----------
 
c initialization of parameters
      
      call utprix('mkPara',ish,ishini,4)

      end

