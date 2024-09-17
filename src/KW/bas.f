C
C  This file is part of EPOS4
C  Copyright (C) 2022 research institutions and authors (See CREDITS file)
C  This file is distributed under the terms of the GNU General Public License version 3 or later
C  (See COPYING file for the text of the licence)
C

c---------------------------------------------------------------------
      subroutine setCounters
c---------------------------------------------------------------------
#include "aaa.h"
      if(jcorona.eq.0.and.(ireadhyt.gt.0.or.icotabr.gt.0)
     &   .or.jcorona.eq.-1 )then
        nevt=0
        nptl=0
        nptlpt=0
      elseif(iappl.eq.4.or.iappl.eq.9)then
        nptlpt=0
      else
        nptlpt=iabs(maproj)+iabs(matarg)
      endif
      end

c---------------------------------------------------------------------
      subroutine checktime(text)
c---------------------------------------------------------------------
#include "aaa.h"
      double precision iutime(5),tidi,tiu3ini,tiu4ini
      common/cchecktime/tiu3ini,tiu4ini
      character text*(*)
      n=index(text,';')
      call timer(iutime)
      if(text(1:n).eq.'start program;')then
        tiu3ini=iutime(3)
        tiu4ini=iutime(4)
      else
        tidi=iutime(3)-tiu3ini + (iutime(4)-tiu4ini)*1d-3
        write(ifmt,'(2a,f8.1,a)')text(1:n),'  cputime: ',tidi/60.,' min'
      endif
      end

c---------------------------------------------------------------------
      subroutine eposevent(n)
c---------------------------------------------------------------------
#include "aaa.h"
#include "ho.h"
#include "ems.h"
      common/photrans/phoele(4),ebeam,noevt /credonoco/kredonoco
      common/cchkengy/ichkengy  /cgbyjmax/gbyjmax
      common/cnfrx/nfrx
      common/cmxpom/mxpom
#if __EB__
      double precision pphsd(5),xphsd(4)
#else
      double precision iutime(5),tiu3,tiu4,tidi
#endif
      common/ifhq/ihq
      !call TestThermal2

      iutime=(/0.0d0, 0.0d0, 0.0d0, 0.0d0, 0.0d0/)
      tiu3=0.0d0
      tiu4=0.0d0
      tidi=0.0d0

      call eventStart(n,ihlle,npommx,kollmx)

      mxpom=0
      iposi=1
      ichkengy=0
      gbyjmax=0
      if(irewch.eq.1)rewind(ifch)
      if(nevent.eq.0)goto 11
      nfr=mod(n-1,nfreeze)
      nev=1+(n-1)/nfreeze
      nfrx=nfr
      if(nrevt.eq.nevent-nfreeze)call wrxx
      if(nrevt.eq.0)call wrxxx
      if(nrevt.eq.0.and.ifrade.gt.0)call hnbcreate 
      do nin=1,iabs(ninicon)
        kredonoco=0
        icou=0
   3    icou=icou+1
        if((ptrmin.gt.0..or.ihlle.eq.1.or.ispherio.eq.1.or.icocore.eq.1)
     .   .and.nin.eq.1.and.nfr.eq.0
     .   .and..not.(iorsdf.eq.5.and.ispherio+ihlle.eq.0))
     .    write(ifmt,'(80a)')('#',i=1,80)
        if((icotabr.eq.0.and.ireadhyt.eq.0.or.jcorona.ne.0)
     .    .and. (nfr.eq.0.or.nin.eq.1) )then !also needed for iphsd=9 for bim
          call aepos(isign(1,-ninicon)*nin,nfr,icou)
        endif
        iexit=0
        ireturn=0
        if(iphsd.eq.-99)stop'useful to check IC for phsd'
#if __EB__
        if(iphsd.eq.0)stop'ERROR: make without XBS needed' 
        if(iphsd.eq.1.or.iphsd.eq.2.or.iphsd.eq.9)then
          if(iphsd.eq.1.or.iphsd.eq.9)then
            num=iabs(ninicon)
            nu=nin
          elseif(iphsd.eq.2)then
            num=nfreeze
            nu=nfr+1
          endif
          nrevt=nrevt+1 
          write(ifmt,'(a,2i5,i9)')'Epos->Phsd ',nu,num,nptl
          if (iphsd.eq.1 .or. iphsd.eq.2) call eposphsd(nu,num)    !mj,  call EPOS+PHSD
          if (iphsd.eq.9 .and. nu.eq.num) call eposphsd9()  ! call pure PHSD, ignore EPOS - not optimal
          if(nu.eq.num)then
            nptlpt=0
            do i=1,num
              call getFromPhsdNptl(i,NrOfParticles)
              write(ifmt,'(a,2i5,i9)')'Epos<-Phsd ',i,NrOfParticles
              do k=1,NrOfParticles
                call  getFromPhsdPtl(i,k,ipdg,
     &          pphsd(1),pphsd(2),pphsd(3),pphsd(4),pphsd(5),
     &          xphsd(1),xphsd(2),xphsd(3),xphsd(4))
                do j=1,5
                  pptl(j,k)=pphsd(j)
                enddo
                do j=1,4
                  xorptl(j,k)=xphsd(j)
                enddo
                idptl(k)=idtrafo('pdg','nxs',ipdg)
                if(abs(idptl(k)).eq.230)then !K0
                  if(rangen().le.0.5)idptl(k)=20 !transform to Ks
                endif
                istptl(k)=0 
                ityptl(k)=60  
                iorptl(k)=0
                jorptl(k)=0
              enddo
              nptl=NrOfParticles
              call decayall(1,99) !to do all weak decays
              call xana  !Analysis for one event
#if           __ROOT__
              if(ifillTree.gt.0.or.ihepmc.eq.1.or.ihepmc.eq.2)call treestore(1)
#endif
            enddo
          endif
        endif
#else
        if(iphsd.ge.1)stop'ERROR: make XBS needed' 
#endif
        if(ihq.eq.1)then
          call exiteposevent(iexit,n,ifevent())
          if(ninicon.eq.1.and.iexit.eq.1)then
            write(ifmt,*)'No c quarks, event skipped'
            ireturn=1
            goto 11
          endif
        endif

        if(ireadhyt.eq.0.and.icocore.ne.0
     .  .and..not.(iorsdf.eq.5.and.ispherio+ihlle.eq.0))then
          if(nfr.eq.0)then
            itrigg=1
            call IniCon(nin,nev,ierrico)
            if(ierrico.eq.1)goto 3
            if(itrigg.ne.1)write(ifmt,'(a)')' REDO aepos'
            if(itrigg.ne.1)goto 3
          else
            do i=1,nptl
              if(istptl(i).eq.5)istptl(i)=7
            enddo
          endif
        endif
        call chkengy(ierrchk)
        if(ierrchk.eq.1)goto 3
        if(icotabm.eq.2)call xana
        call xSpaceTime(iSpaceTime)
      enddo

 11   continue

      call eventEndPrimary(n,ihlle)
      if(ireturn.eq.1)return

#ifndef __BS__
      if(nevent.eq.0)goto 10
      call setCounters
c      call getSystemType(isys)
c      if(isys.le.2)call findHeavyQuarkPairs(100)
      if(noevt.ne.1)nrevt=nrevt+1
      if(ihlle.ne.99)then
#if __ROOT__
        call rhytabx(nfr)
        if(noevt.eq.1)return    !fake DIS
        call hllex(nfr)
        call eosx
        if(nptlpt.eq.0)nptlbd=nptl
        call xfro
#else
        stop "ROOT needed for full hydro"
#endif
      elseif(noevt.eq.1)then                  !fake DIS
        return
      endif
      call aafinal
#if !__BS__ && !__TP__
      if(ihq.eq.1)then
        call timer(iutime)
        tiu3=iutime(3)
        tiu4=iutime(4)
        call aacharm
        call timer(iutime)
        tidi=iutime(3)-tiu3 + (iutime(4)-tiu4)*1d-3
        write(ifmt,'(a,f8.0)')'time for aacharm [sec] : ',tidi 
      endif
#endif
      if(ihacas.eq.1.and.nptl.gt.nptlpt.or.ihacas.ge.2)then
        call hacas(n)
      endif
#if !__BS__ && !__TP__
      if(ihq.eq.1)then
        call timer(iutime)
        tiu3=iutime(3)
        tiu4=iutime(4)
        call aacharm2
        call timer(iutime)
        tidi=iutime(3)-tiu3 + (iutime(4)-tiu4)*1d-3
        write(ifmt,'(a,f8.0)')'time for aacharm2 [sec] : ',tidi 
      endif
#endif
      !++++++++++++++++++++++++++++++++++++++++++++++++ 
      !extra set with ist=-1,-2 not any more used, since utrescxx done 
      !  after freeze out is reliable and robust, so now ist=0,1 represents
      !  the final set respecting energy-momentum conservation
      !But we stil leave the following which may be uncommented for test purposes
      !  while at the same time commenting "call  utrescxx" after freeze out 
      !++++++++++++++++++++++++++++++++++++++++++++++++ 
      !if(irescl.eq.2.or.irescl.eq.3.and.iappl.eq.1)then 
      ! call  utrescxx(iret,irescl-2) !before call afinal  
      ! if(iret.gt.0)call utmsg('##### ERROR 08092014 #####&')
      !endif
      !++++++++++++++++++++++++++++++++++++++++++++++++
      !do ntxt=1,ntxt80
      ! write(*,'(a)')txt80(ntxt)
      !enddo
      !ii1=igetYieldCharged() 
      !ii2= igetYieldCharm12() 
      !if(ntxt80.gt.0)print*,'TEST-Yields',ii1,ii2
      !++++++++++++++++++++++++++++++++++++++++++++++++
      call afinal
      iposi=2
      call xxxSource
      if(ifillTree.gt.0.or.ihepmc.eq.1.or.ihepmc.eq.2)call treestore(1)
      if(ifillH2.eq.1)call d2hstore(nev,nfr)
      if(istore.ge.1.and.istore.le.4) call bstore
      if(istore.eq.5) call ustore
      if(istore.eq.6) call hepmcstore(0)
      if(istore.eq.7) call lhestore(n)
      if(istore.eq.-1) call estore
      call or999
 10   continue
      if(icotabm.ne.2)call xana
      call aseed(1)
#if __ROOT__
      if(nevent.ne.0)call closeEos(nfrx)
#endif
      call closecccptl
      !call memo(0,'exit eposevent;')
#else
      if(noevt.ne.1)nrevt=nrevt+1
      call aafinal
      if(ihacas.eq.1.and.nptl.gt.nptlpt.or.ihacas.ge.2)then
        call hacas(n)
      endif
      call afinal
      iposi=2
      if(istore.ge.1.and.istore.le.4) call bstore
      if(istore.eq.5) call ustore
      if(istore.eq.6) call hepmcstore(0)
      if(istore.eq.7) call lhestore(n)
      if(istore.eq.-1) call estore
      if(icotabm.ne.2)call xana
      call aseed(1)
#endif
      end

c---------------------------------------------------------------------
      subroutine eposStart
c---------------------------------------------------------------------
#include "aaa.h"
      character*1000 line,fnjob
      common/jobfname/  fnjob,nfnjob
      double precision iutime(5)
c dimensions
      call maxsize_create(2)
      call maxsize_set(1, 6 ) !1=mxsplit
      call maxsize_set(2, 150000 ) !2=mx2ptl  ===> CHECK also mx3ptl
c creations
      call deformoptcreate(1)
      call foparamcreate(1) 
      call remn1paramcreate(6)
      call screen1paramcreate(6)
      call screen2paramcreate(5)
      call pfe3paramcreate(5)
      call pfe6paramcreate(2)
      call pdfparamcreate(4) !1=facpdf4,3=facpdf5 
      call saturparamcreate(2) 
      call core1paramcreate(7) !ficoscale
      call core2paramcreate(6) !feloss
      call eventvaricreate(19) 
c aaset
      call aaset(0)
      !call atitle(0)
      call xiniall
      call timer(iutime)
      cotid(1)=iutime(3)
      cotid(2)=iutime(4)
c     to read input from different file name given as argument
      j=-1
      call utword(line,i,j,0)
      fnjob=line(i:j)
      nfnjob=j-i+1
c added 2/9/2022 PBG
      call hqinitransfer
      end

c---------------------------------------------------------------------
      subroutine eventStart(n,ihlle,npommx,kollmx) !start of epos event
c---------------------------------------------------------------------
      call maxsize_get(1, mxsplit )
      if(n.eq.1.or.ihlle.eq.1)then
      call memo(1,'create bor objects ;')
      call xpprbor_create(mxsplit,npommx,kollmx)
      call xmprbor_create(mxsplit,npommx,kollmx)
      call ptprboo_create(2,mxsplit,npommx,kollmx)
      call rapprboo_create(2,mxsplit,npommx,kollmx)
      call gbpom_create(2,mxsplit,npommx,kollmx)
      call idbor_create(2,mxsplit,npommx,kollmx)
      call memo(2,';')
      endif
      end

c---------------------------------------------------------------------
      subroutine eventEndPrimary(n,ihlle) !end of primary of epos event
c---------------------------------------------------------------------
      if(n.eq.ifevent().or.ihlle.eq.1)then
      call memo(1,'destroy bor objects ;')
      call xpprbor_destroy()
      call xmprbor_destroy()
      call ptprboo_destroy()
      call rapprboo_destroy()
      call gbpom_destroy()
      call idbor_destroy()
      call memo(2,';')
      endif
      end

c---------------------------------------------------------------------
      subroutine eposEnd
c---------------------------------------------------------------------
#ifndef __BS__
#if __ROOT__
#include "aaa.h"
#endif
#endif
#ifndef __BS__
      call hnbdestroy
#if __ROOT__
      call closePtlDky
#endif
#endif
      call maxsize_destroy()
      call deformoptdestroy()
      call foparamdestroy() 
      call remn1paramdestroy()
      call screen1paramdestroy()
      call screen2paramdestroy()
      call pfe3paramdestroy()
      call pfe6paramdestroy()
      call pdfparamdestroy() 
      call saturparamdestroy() 
      call core1paramdestroy() 
      call core2paramdestroy() 
      call eventvaridestroy()
      end

c---------------------------------------------------------------------
      subroutine swopen
c---------------------------------------------------------------------
#include "aaa.h"
      if(istore.eq.3) write(ifdt,*) ' 0 0 '
      write(ifmt,'(a)')'rewind copy-file'
      rewind (ifcp)
      nopen=-1
      iecho=0
      end

c---------------------------------------------------------------------
      integer function icheckevtnbs()
c---------------------------------------------------------------------
#include "aaa.h"
      if(nevent.ge.0.or.nfull.gt.0.and.nfreeze.gt.0)then
        icheckevtnbs=1
      else
        icheckevtnbs=0
      endif
      end

c---------------------------------------------------------------------
      subroutine eventsini
c---------------------------------------------------------------------
#include "aaa.h"
      if(icheckevtnbs().eq.0)
     & stop'\n\n STOP in eventsini \n\n'
      if(icotabr.eq.1)nevent=1
      if(nfull.gt.0)then
        if(nfreeze.eq.0)then
          nofreeze=1
          nfreeze=1
        endif
        nevent=nfreeze*nfull
      endif
      if(mod(nevent,nfreeze).ne.0)
     &stop'nevent must be a multiple of nfreeze!!!!!!\n\n '
      if(noebin.ge.0)then
        write(ifmt,'(a,i10,a,f10.2,a)')'generate',nevent
     & ,'  events for engy =',engy,'  ...'
      endif
      call clop(3)
      end

c-----------------------------------------------------------------------
      subroutine aaset(iop)
c-----------------------------------------------------------------------
c     sets parameters and options, initializes counters ...
c-----------------------------------------------------------------------

#include "aaa.h"
#include "ems.h"
#include "ho.h"
#include "par.h"
#include "sem.h"
#include "ico.h"
#include "so.h"
#include "sf.h"
#include "tab.h"
      common/chnbcreate/ihnbcreate
      common/record/maxrec(2),irecty(30,2)
      common/cfacmss/facmss /cr3pomi/r3pomi/cifset/ifset
      common /ems12/iodiba,bidiba  ! defaut iodiba=0. if iodiba=1, study H-Dibaryon
      character*500 fndat,fnncs,fnIIdat,fnIIncs,fnII03dat,fnII03ncs,
     &fndpmjet,fndpmjetpho                  !qgs-II????????
      common/dpmjetfname/  fndpmjet,fndpmjetpho
      common/qgsfname/  fndat, fnncs, ifdat, ifncs
      common/qgsIIfname/fnIIdat, fnIIncs, ifIIdat, ifIIncs !qgs-II????????
      common/qgsII03fname/fnII03dat, fnII03ncs, ifII03dat, ifII03ncs !qgs-II????????
      common/ghecsquel/anquasiel,iquasiel /cjjtb/jjtb,jjeos
      common/cbincond/nozero,ibmin,ibmax  /crapcol/rapcol
      common/photrans/phoele(4),ebeam,noevt
      common/ciuelast/iuelast /ciuskip/iuskip
      common/ciuchaskip/iuchaskip/ciunostring/iunostring
      common/cistat/istateos,istatptl  /cieof/ieof
      common/cgefac/gefac
      common/camaq/amaq(0:3)
      common/cchkengy2/esollxx,eistxx
      common/cijetfluid/ijetfluid /ciotype/iotype
      common/mcen/mcentr,mmxcentr
      common/ccc20/icc20
      common/ciexhd/iexhd
      common/cnfifac/nfifac
      common/producetab/ producetables              !used to link CRMC
      logical producetables                         !with EPOS and QII
      character*30 hdtext(8)
      common/chdtext/hdtext
      common/cnparticip/jproj(2,mamx),jtarg(2,mamx),efluct(6,mamx)
      common/cmodshox/modshox
      common/ifhq/ihq
      common/civirtual/ivirtual
      parameter(maxit=500000)
      common/count/nacc,nrej,naccit(maxit),nptot,npit(maxit)
      integer ioTestFact
      common /ctestfact/ ioTestFact
      integer modeDeltaFzero
      common /cmodeDeltaFzero/ modeDeltaFzero
      integer laddTestFact
      common /claddtestfact/ laddTestFact
      integer iSwitchSides
      common /cSwitchSides/iSwitchSides
      integer noptTestFact
      common /cnopttestfact/ noptTestFact
      integer idhprTestFact
      common /cidhprtestfact/ idhprTestFact
      integer ihqTestFact
      common /chqtestfact/ ihqTestFact
      integer iffsigiEfficiency
      common /ciffsigiEfficiency/ iffsigiEfficiency
      integer ioicoplot
      common /cioicoplot/ioicoplot
      !gauss weights
      data (tgss(2,i),i=1,7)/ .3399810436,.8611363116    ,5*0.     /
      data (wgss(2,i),i=1,7)/ .6521451549,.3478548451    ,5*0.     /
      data (tgss(3,i),i=1,7)/ .2386192,.6612094,.9324700  ,4*0.    /
      data (wgss(3,i),i=1,7)/ .4679139,.3607616,.1713245  ,4*0.    /
      data (tgss(5,i),i=1,7)/ .1488743389,.4333953941,.6794095682
     *                       ,.8650633666,.9739065285    ,2*0.     /
      data (wgss(5,i),i=1,7)/ .2955242247,.2692667193,.2190863625
     *                       ,.1494513491,.0666713443    ,2*0.     /
      data (tgss(7,i),i=1,7)/ .9862838,.9284349,.8272013,.6872929
     *                       ,.5152486,.3191124,.1080549           /
      data (wgss(7,i),i=1,7)/ .03511946,.08015809,.1215186,.1572032
     *                       ,.1855384,.2051985,.2152639           /

      if(iop.eq.1)write(ifmt,'(a)')'default settings ...'
      if(iop.eq.1)goto 1001
      if(iop.eq.2)goto 1002

c  version

                  iversn=3451    !version number
                  iverso=0

c  application

      iappl=1          !choice for application (0,1,2,3,4)

c hq stuff

      ihq=0        !1 = HQ stuff activated, otherwise not
      ivirtual=0   !1 = virtual scattering is active, otherwise not

c  model

      model=1
      iquasiel=1       !allow (1) or not (0) quasi-elastic event in model 3
csp addition for jet-medium interaction
      medium=0

c  file names and units

      fnnx='path/epos '      !path epos name
      fngrv='pdf/ '          !path pdf tables
      fnch='z.check '         !check-file name
      fnhi='z.histo '         !histo-file name
      fndt='z.data '          !data-file name
      fnhm='z.hepmc '      !hepmc-file name
      fncp='z.copy '          !copy-file name
      fnii='aa '       !initl-file name
      fnid='di '       !inidi-file name
      fndr='dr '       !inidr-file name
      fnie='ev '       !iniev-file name
      fnrj='rj '       !inirj-file name
      fncs='cs '       !inics-file name
      fn3p='par '
      fn3d='dec '
      fn3g='f3g '
      fnhpf='zzz.hpf '         !hpf-file name
      nfnnx=index(fnnx,' ')-1   !length of path epos name
      nfnch=index(fnch,' ')-1   !length of check-file name
      nfnhi=index(fnhi,' ')-1   !length of histo-file name
      nfndt=index(fndt,' ')-1   !length of data-file name
      nfnhm=index(fnhm,' ')-1   !length of hepmc-file name
      nfncp=index(fncp,' ')-1   !length of copy-file name
      nfnii=index(fnii,' ')-1   !length of initl-file name
      nfnid=index(fnid,' ')-1   !length of inidi-file name
      nfndr=index(fndr,' ')-1   !length of inidr-file name
      nfnie=index(fnie,' ')-1   !length of iniev-file name
      nfnrj=index(fnrj,' ')-1   !length of inirj-file name
      nfncs=index(fncs,' ')-1   !length of inics-file name
      nfngrv=index(fngrv,' ')-1 !length of pdf path name
      nfn3p=index(fn3p,' ')-1
      nfn3d=index(fn3d,' ')-1
      nfn3g=index(fn3g,' ')-1
      nfnio=3
      fnio(1:3)='---'
      fnus1='none '
      fnus2='none '
      fnus3='none '
      fnus4='none '
      nfnus1=4
      nfnus2=4
      nfnus3=4
      nfnus4=4

      nfnhpf=index(fnhpf,' ')-1 !length of hpf-file name

      ifop=5     !optns-file unit
      ifmt=6     !std output unit
      ifcx=31    !check-file unit (for open)
      ifch=31    !check-file unit (for write)
      ifhi=35    !histo-file unit
      ifdt=51    !data-file unit
      ifcp=52    !copy-file unit
      ifin=53    !input-file unit
      ifhm=54    !hepmc-file unit
      ifnus1=41
      ifnus2=42
      ifnus3=43
      ifnus4=44

      hydt='---'

      producetables=.true.

c  initial seed

c following number should be less than kseq in RMMARD (=2 in EPOS)
      iseqini=2   !sequence number at start program
c seed for random number generator: at start program
      seedi=0d0   !.ne.0.
      iseqsim=1   !.ne.iseqini : sequence number at start program
c seed for random number generator: for first event
      seedj=0d0   !.ne.0.
c place to start for random number generator: for first event
      seedj2=0d0
      call ranfini(0d0,0,-1) !initialize some parameters


 1001 continue ! ----> entry iop=1 <-----

      call factoriel


c  some basic things

      nevent=1    !number of events
      nfull=-1          !number of full events (ico+hydro+f.o.)(if different from -1)
      nfreeze=1    !number of freeze out events per epos event
      ninicon=1    ! +-number of events per initial condition
                    ! if positive: keep same b, otherwise generate b each time
      engy=-1      !energy
      elab=-1      !energy lab
      ecms=-1      !energy centre of mass
      ekin=-1      !energy kinetic
      pnll=-1      !momentum
      rapcms=0     !rap of cms (used in xan for pPb)
      ebeam=-1     !beam energy for proton in fake DIS (pi0-proton collision)
      s0min=0.03  !absolute minimum energy square (reference s0) for reggeon (for Pomeron it is implicitely 1 GeV**2)
      egymin=1.    !minimum energy
      egymax=2.E+06 !maximum energy
      noebin=1     !number of energy bins
      engmin=0     !minimum energy
      engmax=0     !maximum energy
      iologe=-1     !0=linear bins   1=log bins (cms engy)
      iologl=-1     !0=linear bins   1=log bins (Kinetic engy)
      infragm=2    !nuclear fragmentation

c  printout options

      if(iop.eq.0)iprmpt=-2   !prompt (>0) or not (<0)
      ish=0      !1,2,3,4 ...: more and more detailed messages
      irandm=0   !prints all random numbers (1)
      irewch=0   !rewinds check file before each event (1)
      if(iop.eq.0)iecho=1    !verify option for input reading
      modsho=1   !message all modsho events
      modshox=modsho
      idensi=0   !must be 1 when subr xjden1 is used
      ishevt=0   !minimum event number for certain print options
      iwseed=1   !print out seed (1) or not
      jwseed=1   !print out seed in see file (1) or not
      jprint=0   !print out to ifmt
      ihepmc=0
      ihepframe=0 !determine in which frame the collision is defined for the .hepmc output (0:cms, 1:lab)

c  fragmentation and decay parameters

      fkappa=0.014  !String tension (GeV2) for quark string
      fkappag=0.028 !String tension (GeV2) for gluon string
      pud=0.433  !prob u d (from e+e- but used only in epos-dky/vedi)
      pudd=0.88  !d suppression in diquark break (e+e- data)
      puds=0.5   !s suppression in diquark break (important for Xi production in e+e- AND aXi in NA49)
      pudc=0.33  !c suppression in diquark break (??? no data) 0.33
      pudb=0.10  !b suppression in diquark break  !bg
      pmqu=0.003 !mass quark u for string fragm
      pmqd=0.004  !mass quark d for string fragm
      pmqs=0.076 !mass quark s for string fragm
      pmqc=0.9   !mass quark c for string fragm (real mass 1.15<m<1.35)
      pmqb=0.9   !mass quark b for string fragm                      !bg
      pmqq=0.114 !kw! 0.115 !antibaryon production (upper limit)
      strcut=1.   !cut factor for diffractive string fragmentation (1=no suppresion)
      diqcut=0.5  !baryon cut factor for diffractive string fragmentation (needed for pi+p->p/ap data (pz/E > diqcut : no diquark for first node)
      pdiqua= 0.1   !qq-qqbar probability in epos-dro/vedi (decazys only)
      ptfra=  0.35 !string break pt for gluon strings (pp semihard)
      ptfraqq=0.35 !string break pt for quark strings (e+e- and soft)
      ptfrasr=0.   !string break pt increase for diquark
      pbreak=-0.33 !break-parameter (~0.4 to match NA49 data and pi0 spectra fro CR)
                   !if -1<pbreak<0, take pbreak for soft
                   !and e+e- parameterization for hard strings
      pbreakg=0.08 !kw! 0.1 !assymptotic value of pbreak (OK from 0.08 to 0.1)
      zetacut=1.5  !g->ggq2 cut for special hadronization

c  fragmentation and decay options

      ndecay=0   !ndecay suppresses certain resonance decays
                 !0000001: all resonance decays suppressed
                 !0000010: k_short/long (+-20) decays suppressed
                 !0000100: lambda (+-2130) decays suppressed
                 !0001000: sigma (+-1130,+-2230) decays suppressed
                 !0010000: cascade (+-2330,+-1330) decays suppressed
                 !0100000: omega (+-3331) decays suppressed
                 !1000000: pi0 (110) decays suppressed
                 !also several "1"s may be used with obvious meaning
      maxres=99999 !maximum resonance spin
      aouni=0.   !technical parameter for development
      ibreit=0   !Breit-Wigner mass broadening (0=off, 1=on)
      ifoele=0   !forced electron decay of c and b hadrons 

c  lepton-nucleon and e+e- options

      iolept=1     !q2-x ditribution calculated (1) or taken from table (<0)
      ydmin=0      ! range of y
      ydmax=0
      qdmin=0      ! range of q**2
      qdmax=0
      themin=0     !minimum scattering angle
      themax=0     !maximum scattering angle
      elomin=0     !minimum energy of outgoing lepton
      elepti=0     !incoming lepton energy
      elepto=0     !outgoing lepton energy
      angmue=3.9645/360.*2*3.14159 !mue angle
      icinpu=0
      itflav=0     ! initial flavor for e+e-
      idisco=0     !deep inelastic contributions
                   !0=all, 1=direct-light, 2=direct-charm, 3=resolved

c  hadron-hadron options +++

      isetcs=3   !option to obtain pomeron parameters
                 !-2.....do not use fit at all (only good enough for some plots)
                 !-1.....use om11pp directly for the fit (no need for rj tables)
                 ! 0.....determine parameters but do not use Kfit
                 ! 1.....determine parameters and use Kfit
                 ! else..get from table
                 !         should be sufficiently detailed
                 !          say iclegy1=1,iclegy2=99
                 !         table is always done, more or less detailed!!!
                 !and option to use cross ection tables
                 ! 2....tabulation of formula
                 ! 3....tabulation of simulations
                 ! else...not
      iclegy1=1  !energy class lower limit ( 1 minimal    1 realistic    )
      iclegy2=99 !energy class upper limit ( 1  option   99 use of table )
      isigma=0   !option for xsection printing (always calculated now)
                 !  0=no, 1=yes : calculation (not good for ionudi=2)
                 !  2=AA pseudo simulations
      iremn=2    !0 = suppression of multiquark remnant (string end with
                 !    same flavor) -> reduce remnant excitation and
                 !     suppress droplet prod. in remnant
                 !1 = full multiquark remnant (no limitations)
                 !2 = multiquark remnant with valence quark conservation
                 !     and inelastic remnant as low mass droplet only
                 !3 = or suppression of multiquark remnant with
                 !    different string end flavors

c  testing factorization

      !all test should be done using "set  iokoll 1"

      ioTestFact=0 !factorization testing mode
                   !testing procedure should start with large ioTestFact, 
                   !in case of success successively go to smaller values
                   !12 to test fpartone with Esoft=1.
                   !6-11 should be done  using laddTestFact > 0 (for given ladder type)
                   !9,10,11: definitions in rsh, close to the generation procedures 
                   !  11 ... sea-sea, jl=4, fzero~delta(1-x), using shortcuts
                   !  10 ... sea-sea, jl=3, fzero~delta(1-x), using shortcuts   
                   !   9 ... sea-sea, jl=1, fzero~delta(1-x), using shortcuts
                   !   7 ... sea-sea, jl=1, fzero~delta(1-x)   
                   !   6 ... sea-sea, jl=2, fzero~delta(1-x)   
                   !5,4... should be done  using laddTestFact > 0, 
                   !in case of success also for laddTestFact = 0
                   !   5 ... sea-sea, jl=1
                   !   4 ... sea-sea, jl=2    
                   !-3 to 3 full psahot tests, essentially for laddTestFact = 0
                   !   3 ... sea-sea
                   !  -1 ... val-val 
                   !  -2 ... sea-val
                   !  -3 ... val-sea
                   !partial tests, needed in case of problems for -3 
                   !  -4 ... val-sea,qg
                   !  -5 ... val-sea,qq
                   !  -6 ... val-sea, but only count partons from qg 
                   !  -7 ... val-sea, but only count partons from qq 
                   !normal case
                   !   0 ...  full run

      laddTestFact=0 !forcing ladder type for testing (emission means "at least one emission"
                     !  0 ... no forcing
                     !  1 ... no emissions
                     !  2 ... proj side emissions
                     !  3 ... targ side emissons
                     !  4 ... both sided emissions 

      modeDeltaFzero=0 !interpretation of fzero, variable set automatically depending on iotestFact
                       !  0 ... normal case
                       !  1 ... fzero considered to be delta(1-x) 

      noptTestFact=0 !ladder end pt forced to zero and Q2s to Q2s_min
                     !useful for tests with ioTestFact=0
                     !  0 ... no
                     !  1 ... yes 

      idhprTestFact=999 !Pomeron type forced to be idhprTestFact (if NE 999)
                        !useful for tests with ioTestFact=0
                        !  999 ... no forcing
                        !  0 ..... sea-sea
                        !  1 ..... val-sea
                        !  2 ..... sea-val
                        !  3 ..... val-val 

      ihqTestFact=0 !factorization test only for outBorn HQ (1)

      iffsigiEfficiency=0 !efficient ffsigi calculation in batch (allows to avoid using parametrization)
                        !  0 ... no    (has to be changed via set 
                        !  1 ... yes    before nfull)

c  hadron-hadron options

      idprojin=1120 !projectile particle id
      idtargin=1120 !target particle id
      idproj=1120 !projectile particle id
      idtarg=1120 !target particle id
      iregge=1    !consider reggeons (1=yes 0=no)
      isopom=1    !consider soft pomerons (1=yes 0=no)
      ishpom=1    !consider semihard pomerons (1=yes 0=no)
      iscreen=1   !consider saturation effects (1=yes 0=no)
      iq2sat=1    !consider variable Q20 (0=fixed Q20=q2nmin, 1=no Npom dep and Nbin scal, 2=assymmetry with max first, 3=assymmetry with min first)
      irmdrop=0    !consider droplet for remnants (0=no,1=yes, 2=multiple strings)
      nprmax=10000 !maximum number of pomerons/reggeons
      intpol=3     !number of points for interpolation in psfli, psfaz
      ioems=2
      ipsahot=0   !call testpsahot (1)
      iseahot=-1   !which sea-sea component is plotted in xEmsP2 (-1 all 0=gg 1=qg 2=gq 3=qq)
      iomega=0    !option for omega calculation (if 1, no mass transfer, if 2, no mass diagram in G)
        !hadron excitation classes (used in psvin)
      icdp=2     !projectile hadron excitation
      icdt=2     !target hadron excitation
              !hadron classes: 1=pions,2=nucleons,3=kaons
      iclpro1= 2 !1   !projectile hadron class lower limit
      iclpro2= 2 !3   !projectile hadron class upper limit
      icltar1= 2   !target hadron class lower limit (should not be change (see epos-sem))
      icltar2= 2   !target hadron class upper limit (should not be change (see epos-sem))
      egylow=1.5  !lower limit of lowest energy bin
      delh=0.25   !effective overcriticity for the hard pom (techn param)
      factgam=1.         !enhancement factor for gammas
      factjpsi=1.0      !factor for jpsi production in saturated contribution
      factsat=1.0      !factor for saturated contribution
      fracsat=0.  !unused remove from sem.h
      factsft=0.  !unused remove from sem.h
      alpsft=0.45 !soft qq dependence exponent
      alpsat=1.0 !sat qq dependence exponent
      ptphot=9.999     !bg
      fInitialHF=0  !unused remove from sem.h
      disize=0.0 !dipole size
      call pdfparamset(1, 0.2 )!facpdf4
      call pdfparamset(2, 1.0 )!facpdf4sat
      call pdfparamset(3, 0.2 )!facpdf5
      call pdfparamset(4, 1.0 )!facpdf5sat

c  hadron-hadron parameters +++

      betpomi=   0.25   !gluon structure function parameter for the so
      glusea=    0.1    !sea quarks / (sea quarks + gluons) at x -> 0
      r2had(2)=  1.2    !r2 proton
      r2hads(2)= 1.25   !diff corr factor proton
      slopom=    0.06   !alpha prime pomeron
      slopoms=   0.     !not used
      gampar=    1.     !gamma parton soft ->  sig up, n up, softPom up
      gamtil=    0.08   !increase -> sig up, n up, hard Pom up
      alppar=    0.33   !alpha particip soft (not 1 !) increase -> sig up, n down, width y down
      alpparh=   0.     !alpha particip hard (not 1 !) increase -> sig up, n down, width y down
      alppom=    1.075  !alpha pomeron
      ptsend=    1.     !string end pt
      ptsendi=   0.3    !string end pt for non excited state and diquark
      ptsems=    0.2    !mass of string end (minimum pt)
      ptsecu=    1.     !cut factor for gaussian distribution
      ptdiff=    0.32   !pt for diffraction
      ioangord= 0  !angle ordering in space-like cascade (1) or not (0)
      coekaf=   1.0 !coefficient kappa_F for factorization scale
      q2nmin=    4.     !minimum used q2

      q2zmin=    1.   !absolute minimum for tables used for binning
             !also used as reference for tabulation
          ! SHOULD NOT BE CHANGED    

      q2dmin=    1.     !has to be properly defined in epos-qsh
                        !(meaning of q 2 m i n    changed, therefore q 2 m i n
                        !      in epos-qsh renamed q2dmin)
      q2cmin(1)=-99999  !current  q 2 m i n    proj (x dependent)
      q2cmin(2)=-99999  !current  q 2 m i n    targ
      q2pmin(1)=-99999  !maximum pair q 2 m i n    proj (x independent)
      q2pmin(2)=-99999  !maximum pair q 2 m i n    targ
      q2min= -99999     !should not be set here any more
      q2ini=     0.25   !resolution scale
      q2fin=     0.014  !q**2 cutoff timelike      decrease ->  high pt down
      amdrmax=   30.    !maximum mass leading droplet (<50 for stability)
      amdrmin=   10.    !minimum mass leading droplet
      facdif=    0.4    !factor for diffractive profile
      edifac=    0.333  !transition reggeon-diffraction
      ediflim=   15.    !limit in energy for increase of diffraction
      facmc=     1.0    !correction factor to match MC simulations (should be 1.)
      reminv=    0.13   !remnant inversion probability (inversion important for forward pi(0) spectra : consequences on Xmax)

      mxzzvex=1
      zzvexch(1)='01' !we keep (and have to define) the arrays though just in case a new variable would be need later

      rstrau(1)=1.   !pion !effective ratio of u sea over u sea basis
      rstrad(1)=1.         !effective ratio of d sea over u sea basis
      rstras(1)=0.2        !effective ratio of s sea over u sea basis (kaons in pipp250)
      rstrac(1)=1.         !effective ratio of c sea over u sea basis
      rstrau(2)=1.   !nucl !effective ratio of u sea over u sea basis
      rstrad(2)=1.         !effective ratio of d sea over u sea basis
      rstras(2)=0.7        !effective ratio of s sea over u sea basis
      rstrac(2)=1.        !effective ratio of c sea over u sea basis
      rstrau(3)=1.   !kaon !effective ratio of u sea over u sea basis
      rstrad(3)=1.         !effective ratio of d sea over u sea basis
      rstras(3)=0.2        !effective ratio of s sea over u sea basis (kaons in kpp250)
      rstrac(3)=1.        !effective ratio of c sea over u sea basis
c      rstrau(4)=1.   !j/psi!effective ratio of u sea over u sea basis
c      rstrad(4)=1.         !effective ratio of d sea over u sea basis
c      rstras(4)=0.8        !effective ratio of s sea over u sea basis
c      rstrac(4)=0.2        !effective ratio of c sea over u sea basis
      rstrasi=0.0      !effective ratio of strange sea over u sea increase
c wgtqqq (<1) define the probability that the diquark is taken
c (or introduced in  a meson) directly from the remnant
c it is for stopping not for baryon prod (put baryon in the center but do
c not change very large x of mesons). Value fixed with forward baryon
c in pion interactions. Active only with iremn>1
      wgtqqq(1)=0.22   !weight for val diq (as soft string ends for one pomeron) for pion
      wgtqqq(2)=0.22   !weight for val diq for nucleon
      wgtqqq(3)=0.22   !weight for val diq for kaon
c      wgtqqq(4)=0.22   !weight weight for val diq for J/Psi
c within 1-wgtqqq not to take a diquark, wgtdiq (<1) is
c the absolut probability to create a diquark as string end
      wgtdiq=0.15      !weight for seadiq - antidiq as soft string end
c in the 1-wgtqqq-wgtdiq probabilty not to have a q-aq string ends,
c wgtval is the probability to take the valence quark in soft interactions
c wgtsea is the probability to take a q-aq pair from the sea,
c if no valence quarks are available, then we use sea quarks
c these values can be arbitrary choosen since wgtval/wgtsea is used
      wgtval=0.15        !weight for valq - antiq as soft string ends
      wgtsea=0.85       !weight for seaq - antiq as soft string ends
      exmass=0.02         !excitation mass for remnant
      r3pom=0.1      !triple pomeron coupling (not used)
      wexcit=0.       !excitation in fremnu (for DIS)
      wproj=0.        !not used
      wtarg=0.        !not used
      q2sft=1.5      !minimum Q2 to have QCD perturbative calculation 
c     gamhad(1) defined in psaini
c     gamhad(2) defined in psaini
c     gamhad(3) defined in psaini
      gamhadsi(1)=0.55    !correction factor to gamma soft pion
      gamhadsi(2)=1.      !correction factor to gamma soft nucleon
      gamhadsi(3)=0.47    !correction factor to gamma soft kaon
c      gamhadsi(4)=-1.    !correction factor to gamma soft charm
      r2part=0.05     !r2 parton
      r2har(1)=0.05    !r2 hard pion
      r2har(2)=0.05    !r2 hard baryon
      r2har(3)=0.05    !r2 hard kaon
c      r2har(4)=0.05    !r2 hard charm
      r2had(1)=1.5    !r2 pion
      r2had(3)=0.8    !r2 kaon
c      r2had(4)=0.     !r2 charm
      r2hads(1)=1.1   !diff corr factor pion
      r2hads(3)=1.1   !diff corr factor kaon
c      r2hads(4)=1.25  !diff corr factor kaon
      wdiff(1)=0.5    !diffractive probability
      wdiff(2)=0.7    !diffractive probability
      wdiff(3)=0.53   !diffractive probability
c      wdiff(4)=0.1    !diffractive probability
      alplea(1)=0.7  !alpha leading pion
      alplea(2)=1.    !alpha leading proton
      alplea(3)=0.7    !alpha leading kaon
c      alplea(4)=1.    !alpha leading jpsi
      rexddf=25.  !high mass double diffraction increase factor
      xmxmas=0.5  !xmxrem(ipt)=xmxmas/(2*m(ipt)*R(itp)) : maximum momentum fraction allowed to be given 
                  !for remnant mass determination (very important for multiplicity and xf distri at low energy). 
                  !ipt = pro or tar and itp corresponding tar or pro.
      xmindiff=1.35  !factor for minimum energy of pion exchange in prorem 
                  !(change multiplicity and xf distribution of protons and neutrons)
      xminremn=1.25  !factor for minimum energy in prorem (change multiplicity of all remnants)
      rexres(1)=0.175 !pion remnant low mass excitation probability in nucleus
      rexres(2)=0.1   !nucleon low mass remnant excitation probability in nucleus
      rexres(3)=0.5   !kaon low mass remnant excitation probability in nucleus
c      rexres(4)=1.    !charm low mass remnant excitation probability in nucleus
c      rexres(1)=1.5  !relative width for diffractive Pomeron in pion
c      rexres(2)=1.   !relative width for diffractive Pomeron in nucleus
c      rexres(3)=1.   !relative width for diffractive Pomeron in kaon
c      rexres(4)=1.    !relative width for diffractive Pomeron in charm
      alpdif=0.32     !alpha mass diffractive energy dependance for cross section and metropolis
      alpdi(1)=1.     !alpha low mass diffractive for remnant without diquark
      alpdi(2)=1.     !alpha low mass diffractive for remnant with diquark
      alpndi(1)=1.    !alpha low mass non-diff for remnant without diquark
      alpndi(2)=1.    !alpha low mass non-diff for remnant with diquark
      alpsea=0.3      !alpha string end x for sea parton
      alpval=0.3      !alpha string end x for valence parton
      alpdiq=0.3      !alpha string end x for sea diquark
      ammsqq=0.28     !minimum mass string quark quark
      ammsqd=1.08    !minimum mass string quark diquark
      ammsdd=1.88    !minimum mass string diquark diquark
      delrex=0.5     !excitation mass to be added to the minimal remnant mass when remnant is not connected to string (nuclear splitting)
      alpdro(1)=2.   !factor for minimum mass of leading droplet (not less than 1.5 for kinematic reasons)
      alpdro(2)=0.3   !pt of leading droplet
      alpdro(3)=1.6  !alpha mass of leading droplet (iept=3)
      ptipom=1.5   !intrinsic pt of semi-hard Pomerons
      ptipomi= 1.  !increase of intrinsic pt with q2s
      ptipos=999.   !intrinsic pt of sat Pomerons                   ===> will be overwritten in key.f 
      ptiposi=999.  !increase of intrinsic pt with q2s for sat pom  ===>   unless both are zero
      facq2tim=1.  !factor for Q2 arguments of timsh2 

c  hadron-hadron parameters

      edmaxi=    1.e11  !defines edmax in epos-sem
      epmaxi=    1.e11  !maximum energy in semi-hard tables
      qcdlam=.04        !lambda_qcd squared
      nbflav=5          !maximum number of flavors in born MC and timelike 
      naflav=99999  !KW1811 should not appear any more
      noflit=3   !number of light flavors !KW1811 
      nofeff=5   !used in alpha_s (pssalf) !KW1811 
      nrflav=3       !number of active flavors in remnant and string fragm
      facnof=3.33    !factor to suppress initial charm at low q2 
      factk=2.         !k-factor value
      alfe=1./137.
      pt2cut=0.        !p_t-cutoff for the born process
      adskap=4.*(0.513)**2/qcdlam    !non-perturbative alpha_s parameter from AdS/QCD 
                      !(recalculated from our value of qcdlam and transition scale in psaini)
      cumpox=1.     !cutoff mass for pomerons
      xcupom=1.e-5  !x threshold for low mass virtual diffractive pomerons
      xcutsf=1.e-3   !1/xmin for semihard pom (not used)
      r3pomi=r3pom  !store

c dibaryon stuff

      iodiba=0.      ! if iodiba=1, study H-Dibaryon (not used (see ProRef in epos-ems ????)
      bidiba=0.030   ! epsilon of H-Dibaryon

c following parameters recalculated in xsigma ...

      alphigh(1)=0.5   !remnant without diquark mass index for MC (high mass)
      alphigh(2)=0.5   !remnant with diquark mass index for MC (high mass)
      rexdif(1)=0.5    !remnant low mass probability diffractive pion
      rexdif(2)=0.5  !remnant low mass probability diffractive proton
      rexdif(3)=0.5    !remnant low mass probability diffractive kaon
c      rexdif(4)=0.    !remnant low mass probability diffractive charmed
      rexpdif(1)=0.25    !remnant pomeron mass probability diffractive pion
      rexpdif(2)=0.25    !remnant pomeron mass probability diffractive proton
      rexpdif(3)=0.25    !remnant pomeron mass probability diffractive kaon
c      rexpdif(4)=0.    !remnant pomeron mass probability diffractive charmed
      rexndi(1)=0.4   !remnant excitation probability nondiffractive pion
      rexndi(2)=0.55  !remnant excitation probability nondiffractive proton
      rexndi(3)=0.6   !remnant excitation probability nondiffractive kaon
c      rexndi(4)=1.    !remnant excitation probability nondiffractive charmed
      !... up to here.

c screening splitting +++

      !Note : cross section/saturation value, change inelasticity and not
      !only cross section
      epscrw=0.1      !overall factor for Z (Zsame,Zother)-> pp xsect .... w_Z
      epscrp=3.       !b width param for Z     -> pp xsection ............ w_B
      egyscr=3.     !screening minimum energy -> pp xsection ........... s_M
      epscrd=0.1        !screening power for diffractive part
      epscrs=0.1    !screening power increase soft  -> pp xsctn ........ alp_S
      epscrh=2.5    !screening power increase hard  -> pp xsctn ........ alp_H
      epscrxi=0.42  !screening power maxi          -> pp xsctn
      zbcut=-1.     !factor to fix impact parameter at which Z saturates (in unit of bsat given by epscrx) if <0, all Z reduction when b->0
      epscrb=1.     !epsilon reduction factor when b->0 for xs     ->  multiplicity
      epscrg=0.     !factor to increase Z around b=0 only if zbcut<0
 
      !nuclear part of Z
      znurho=0.    !increase of Z due to nuclear density effect  -> pA xsctn (low E)
      zdfbmx=0.    !fraction of bkmx below which a pair is used to compute nuclear part of Z (diff Z)
      zbrmax=0.    !factor in Phiexpo to define probability to add a pair to nuclear screening
      zrhoinc=0.   !factor of nuclear Z for saturated non-diffractive part of Z
      zbrads=0.01  !limit on phi to define saturation 
      zfluct=0   !add fluctuations for pA

      !Z parameters
      zodinc=0.  !increase remnant excitation probability due to q2s
      zipinc= 0. !q2s dep pbk string fragmentation /mult/ ~0.24 e+e-
      zopinc= 1. !radius for energy density used for string fragmentation
      zoeinc= 0.6      !q2s factor in dels
      zdfinc=0.       !x dep. of alpha for the remnant mass
      zdrinc= 0.     !increase of droplet minimum mass (iremn>=2)
      zmsinc= 0.5    !?
      xzcut=0.   !factor for minimum x for a Pomeron to be used for nuclear splitting

c Reggeon parameters

      alpreg(1)=0.3   !alpha_reggeon
      alpreg(2)=0.5   !alpha_reggeon
      alpreg(3)=0.3   !alpha_reggeon
      sloreg(1)=0.74   !slope_reggeon
      sloreg(2)=0.36   !slope_reggeon
      sloreg(3)=0.57   !slope_reggeon
      gamreg(1,-1)=28.9  !gamma_reggeon pim+N
      gamreg(1,0)=25.5   !gamma_reggeon pi0+N
      gamreg(1,1)=23.5   !gamma_reggeon pip+N
      gamreg(2,-1)=40.7  !gamma_reggeon Nbar+N
      gamreg(2,0)=26.46   !gamma_reggeon strangeB+N
      gamreg(2,1)=16.46   !gamma_reggeon N+N
      gamreg(3,-1)=17.7  !gamma_reggeon kp+N
      gamreg(3,0)=12.56   !gamma_reggeon k0+N
      gamreg(3,1)=7.56   !gamma_reggeon km+N
      r2reg(1,-1)=0.55   !r^2_reggeon pim+N
      r2reg(1,0)=0.37    !r^2_reggeon pi0+N
      r2reg(1,1)=0.17    !r^2_reggeon pip+N
      r2reg(2,-1)=2.24   !r^2_reggeon Nbar+N
      r2reg(2,0)=1.4    !r^2_reggeon strangeB+N
      r2reg(2,1)=0.613    !r^2_reggeon N+N
      r2reg(3,-1)=0.79   !r^2_reggeon km+N
      r2reg(3,0)=0.66    !r^2_reggeon k0+N
      r2reg(3,1)=0.52    !r^2_reggeon kp+N

c  masses

      amhadr(1)=.14            !pion mass
      amhadr(2)=.939           !nucleon mass
      amhadr(3)=.496           !kaon mass
      amhadr(4)=1.868          !d-meson mass
      amhadr(6)=1.116          !lambda mass
      amhadr(5)=2.27           !lambda_c mass
      amhadr(7)=.548           !eta mass
      amhadr(8)=3.097          !J/psi mass
c      qcmass=1.6               !c-quark mass  (in idmass = 1.2 )
c      amhdibar=2.200           !h-dibaryon mass       !not used any more
      qmass(0)=pmqq           !diquark effective bounding energy (for pt distribtions)
      isospin(0)=0

      ! qumass, qdmass...,  qmass(1-6), isospin(1-6) defined in ainit

      amaq(0)= 0.
      amaq(1)= 0.337
      amaq(2)= 0.337
      amaq(3)= 0.486

c  nucleus-nucleus

      iokoll=0      !fix # of independent collisions (given by iokoll if not 0)
      laproj=0      !projectile charge number
      maproj=0      !projectile mass number
      latarg=0      !target charge number
      matarg=0      !target mass number
      core=0.34     !hard core distance(2*0.17)
c      ncolmx=100000 !maximum number of collisions
      fctrmx=10     !parameter determining range for density distribution
      bmaxim=10000  !maximum impact parameter
      bminim=0.     !minimum impact parameter
      bref80=0.     !reference impact parameter for 80%
      phimax=2*3.1415927 !maximum phi
      phimin=0      !minimum phi
      nq2mnfix=0  ! -1 (table making, batch mode)  1 (testing) 0 (else) 

      asatur=0.
      bsatur=0.
      csatur=1.
      dsatur=0.
      esatur=0.
      zzsoft=0.04   !soft attenuation for Q2s

      ipytune=350   !Pythia Perugia Tune (2011) for PYTHIA run
      ionudi=3      !nuclear diffraction included (>0) or not (0)
                    ! = 0 for RHIC nuclear data based on glauber
                    !     (no event with 0 collision "a la Glauber")
                    ! = 1 for cosmic ray simulations (diffraction without
                    !      excitation counted as inelastic
                    !      to be consistent with sigma_ine used for CR)
                    ! = 2 for fixed target accelerator data cross section
                    !     (diffractive without projectile excitation = elastic)
                    ! = 3 real inelastic (no trigger effect)
                    !     (diffractive without excitation = elastic)

      iotype=0   ! 0=inelastic  1=nondiffr

c multiple scattering

      ikolmn=0         !minimum number of inelastic collisions
      ikolmx=1000000   !maximum number of inelastic collisions
      nglmin=0         !minimum number of inelastic NN collisions
      nglmax=1000000   !maximum number of inelastic NN collisions
      segmin=0         !segm mult min
      segmax=1000000   !segm mult max

c trigger

      ptrmin=-1.0
      ioecc=0      !1,2,3,4 triggers on quadrant in ecc2-ecc3 plane
      valecc=0.0   !defines quadrants

c rescattering parameters +++

      iorsdf=3      !droplet formation
                    ! 0 ... turnedoff
                    ! 3 ... string segment ("effective" by default)
                    ! 5 ... string
      iorsce=0      !color exchange turned on(1) or off(0)
      iorshh=0      !other hadron-hadron int. turned on(1) or off(0)
      iocluin=1     !include inwards moving corona segments into clusters (1) or not (0)
      ioquen=0      !jet quenching option (0=no)
      tauzer=1.     !tau for core-corona procedure 
      taurem=0.0    !decay time of remnant elements
      taustr=1.    !formation time of ptls from string
      tauhac=0.0    !formation time in eposu
      iostro=1      !z-t of string origin set zero (0), x-y kept
      nsegsu=30     !number of segments per subcluster
      dsegce=11.68  !minimum number of hit segments for core per cell volume
      kigrid=1 ! 999 !now used in jintpo
      fxcell=6.0    !binning factor jintpo
      fzcell=1.0    !binning factor jintpo
      ptlow=1.0     !for core/corona separation
      ptupp=3.0     !for core/corona separation
      qufac=1.0    !energy loss factor ! .le.0 no Eloss
      quexpo=-9999    !not used presently
      ijetfluid=2   ! 1 = skip jetfluid
      cutdxy=1      !cutoff mass for dxy
      fludiq=0.25
      amimfs=1.0    !below this: elastic
      amimel=0.050  !below this: nothing
      delamf=1      !above this: color exch  !cutoffs for kinetic energy
      deuamf=1      !above this: nothing     !mass - minimum mass of pair
      epscri(1)=0.15!energy density for hnbaaa
      epscri(3)= -1 !read in from table
      rapcol=1.0    !rap maxi for coll in coload
      leadcore=1  !different options for leading particles in core

c core relevant

      call core2paramset6(0.,0.,0.,0.,0.,0.) !feloss
      radeft1=0.05 !radius of elementary flux tube (->max(radeft,cell size)
      radeft2=0.05 
      facposz=      0.33 !factor for Z dependence in setfacpos
      facposf=      1.0  !factor in setfacpos
      tauzer1=    0.70 !core formation ftime pp
      tauzer2=    1.00 !core formation ftime AA 
      nsegmincore=  1e7 !min number of segments for core
      ratiomaxcore= 0  !max ratio of eloss / core

c rescattering parameters

      amsiac=.8     !interaction mass
      amprif=0      !print option
      delvol=1      !print option
      deleps=1      !print option
      deltau=0.2    !initial delta_tau for space-time evolution
      numtau=80     !number of tau steps for space-time evolution
      dlzeta=0.5    !delta_zeta for longitudinal droplet binning
      etafac=1.0    !factor determining inner range
      facnuc=0.0    !factor for nuclear size to determine inner range
      hacore=-1.0    !hadron compressibility factor
      cepara=0.03   !parameter for excitation for color exchange
      dscale=1.    !scale parameter for hadron-hadron
      iceopt=1      !real color exchange (1) or just excitation (0)
      ihacas=0    ! call hacas (1) or not
      iuelast=0    !force elastic scattering (1) or not
      iuskip=0    !only decay in hacas, skip reaction (1)
      iuchaskip=0    !skip reaction certain channels (1)
      iunostring=0    !no strings in hacas (1)

c core (initial conditions for whatever, on hyperbola)

      cutico=3 !1    !cutoff parameter for string smoothing kernel in epos-ico
                    !  (as small as possible with stable results)
      sgzico=1.0    !sigmaz factor in ico
      dssico=0.2    !s step size for string integration in epos-ico
                    !  (as big as possible with stable results)
      nsmooth=0     !binomial smoothing parameter
      iranphi=0     !ranphi rotation (1) or not (0)
      jcorona=0     !force corona calc (1 = all, 2 = just geom)
                    !jcorona= -1 : NO corona in any case
      icocore=0     !consider core initial condition (1 or 2)
      icospec=1     !consider spectators
      icoremn=1     !consider hadrons from remnants
      icostrg=1     !consider strings
      icotabm=0     !make table for initial condition stuff
      icotabr=0     !read table for initial condition stuff
      nxico=0
      nyico=0
      nzico=0
      xminico=0.  !xrange for initial condition calculation
      xmaxico=8.     !.ne.0 to be able to run "core effective" (default)
      yminico=0.  !yrange for initial condition calculation
      ymaxico=0.
      zminico=0.  !eta range for initial condition calculation
      zmaxico=0.
      ioicoplot=0 !call xIcoPlot

c phsd 

      iphsd=0      !0 = no phsd, 9 = only phsd
                   !1 = epos+phsd, ninicon > 1
                   !2 = epos+phsd, nfreeze > 1

c eos making tables

      ieostabm=0    !make eos tables

c  spherio

      ispherio=0
      jspherio=0
      call aasetspherio
      igethyt=0
      ireadhyt=0

c  hlle

      ihlle=99       !"effective" by default
      tfrout=.166
      kfrout=22  !1 = GC(T), 2 = GC(eps), 22 = MiC(eps), 3+4+12 = Test  
      efrout=0.0    ! FO epsilon
      nfrout=10000 !number of eta bins for building subclusters for MiC
      fofac=1.0
      fofai=.0
      ntaumx=1e9
      ioeos=0
      iozerof=0
      etaos=0.08
      zetaos=0.0
      iochem=1 ! 1 = normal ; 0 = chem potential set zero

c hydro

      tauone=1e30   !hydro tau step / 2 for tau >  tauone
      tautwo=1e30   !hydro tau step / 4 for tau >  tautwo
      tauup=50      !upper limit for tau - tauzer for hydro
      floini=0      !initial flow parameter
      fdtau=1.0
      volex=0
      itabptl=3
      nradhy=75
      nphihy=61
      nxhy=0
      nyhy=0
      nzhy=0
      ntauhy=0
      xminhy=0.   !xrange for hydro output
      xmaxhy=0.
      yminhy=0.   !yrange for hydro output
      ymaxhy=0.
      zminhy=0. !etarange for hydro output
      zmaxhy=0.
      ifaahlle=0
      ifathlle=0
      ifazhlle=0
      epsfin=99999.
      dtauhy=0    !tau step for hydro
      maxsurf=4
      ienvar=0  !activate envar (1) or not (0) !1 = very slow! only test
      izmode=1
      do i=1,7
      jzmode(i)=0
      jtable(i)=0
      enddo
      irclass(1)=0
      irclass(2)=0
      ihdim(1)=0
      ihdim(2)=0
      ihdim(3)=0
      iSource=0
      iHyEpsilon=0
      iHyEpsilon2=0
      iHyEntropy=0
      iHyTemperature=0
      iHyRadVelocity=0
      iHyLongVelocity=0
      iHyAverages=0
      iHyBaryon=0
      iHyFoVol=0
      iHyFoRadius=0
      iHyFoEpsilon=0
      iHyFoRadVelocity=0
      iHyFoTangVelocity=0
      iHyFoBarmu=0
      iHyEpsilonEta=0
      iHyBaryonEta=0
      iHyBarmuEta=0
      iHyEntropyEta=0
      iHyEos=0
      iHyEmpty=0
      iHoEpsilon=0
      iHoEpsilonEtas8=0
      iHoEpsilonEtas6=0
      iHoTemperatureEtas8=0
      iHoTemperatureEtas6=0
      iHoEpsilonTau8=0
      iHoEpsilonTau6=0
      iHoTangVelTau6=0
      iHoRadVelTau6=0
      iHoRadVel=0
      iHoTanVel=0

c  cluster decay

      ioclude=3     !cluster decay option     ("effective" by default)
            ! ioclude=4  in CooperFrye hadron masses (normal case) 
            ! ioclude=5  in CooperFrye quark masses (mimic coalescence) 
      amuseg=3.0    !min mass for radial boost (limit=amuseg+yrmax(E))
      yradmx= 0.04  !max radial collective boost
      yradpi=2.5     !fix the eta shape at high pt
      yradpx=0.333     !exponent for pt random distribution
      yradpp=4.     !mass dependence of low mass flow (attenuation for heavy part)
      facecc=0.5    !eccentricity parameter
      rcoll=0.0     !radial collective flow param
      bag4rt=0.200  !bag constant ^1/4
      vrad=0.3
      facts=0.35     !gamma-s factor
      factb=0.35     !gamma-s factor baryons
      factq=1       !gamma-qqbar
      facmicro=1 !plot option

c  droplet decay initializations

         asuhax(1)=1.134  !two lowest multiplets
         asuhax(2)=1.301  !two lowest multiplets
         asuhax(3)=1.461  !two lowest multiplets
         asuhax(4)=1.673  !two lowest multiplets
         asuhax(5)=0.7700 !two lowest multiplets   rho
         asuhax(6)=0.8920 !two lowest multiplets   K*
         asuhax(7)=1.2320 !two lowest multiplets
         asuhay(1)=0.940  !lowest multiplet
         asuhay(2)=1.200  !lowest multiplet
         asuhay(3)=1.322  !lowest multiplet
         asuhay(4)=1.673  !lowest multiplet
         asuhay(5)=0.1400 !lowest multiplet
         asuhay(6)=0.4977 !lowest multiplet
         asuhay(7)=1.2320 !lowest multiplet

c  droplet specification

      keu=0     !u flavour of droplet
      ked=0     !d flavour of droplet
      kes=0     !s flavour of droplet
      kec=0     !c flavour of droplet
      keb=0     !b flavour of droplet
      ket=0     !t flavour of droplet
      tecm=10   !energy of droplet
      volu=70   !volume of droplet

c  metropolis and grand canonical

      ihnbcreate=0
      fitermet=2.0  !metropolis iterations factor
      felamet=0.50 !metropolis elastic factor
      iospec=24 !option for particle species (24 = full set, all others are tests)
      iocova=40 !LIPS for n <= iocova else NRPP 
                !  (to be choses such that both methods work for this value)  
      iopair=2  !double pair method (2), otherwise (3,4) test cases
      iozero=65 !relative weight of zeros (compared to hadrons)
                ! (-1) nspecs
                ! (-2) nspecs/sqrt(tecm/volu)
      ioflac=1  !test multipl distr without (1) or with (2) flavour conserv
                !  (2 only good for nspecs=3,7)
      iostat=1  !use boltzmann (1) or quantum (0) statistics in hgc-routines
      ioinco=1  !call hnbmin for initial configuration (0)
                !call hgcnbi for initial configuration to generate better
                !initial configuration (1)
                !call hgcnbi for initial configuration to generate optimal
                !initial configuration (2)
      iograc=1  !call hgcaaa in case of ioinco=0 (1)
      epsgc=2.  !required accuracy in hgcaaa 10**(-epsgc)
      iocite=0  !perform counting at metropolis iterations (1) or not (else)
      ioceau=0  !perform counting for exp. autocorrel. time (1) or not (else)
      iociau=0  !perform counting for int. autocorrel. time (1) or not (else)
      ioinct=0  !test grand canonical metropolis sampling (1)
                !to plot results call xhgccc, xhgcfl and xhgcam
      ioinfl=1  !conserve flavor in initial configuration in hgcnbi (1)
                !do not conserve flavor (0)
                !do not conserve flavor and energy (-1)
      iowidn=2  !width of total multiplicity distribution in hgcnbi
                ! sigma_tot -> sigma_tot/iowidn
                ! >0 unnormalized
                ! <0 normalized
      ionlat=2  !determine nlattc ,old method (0)
                !or determine nlattc in hgcnbi as:
                ! (1) max(1.3*<N>,<N>+2*sigma,6)
                ! (2) max(1.5*<N>,<N>+3*sigma,6)
                ! (3) max(2.0*<N>,<N>+4*sigma,6)
      iomom=1   !number of momenta to be changed in hnbodz
      ioobsv=0  !observable for autocorrelation time calculation
                !0: total multiplicity
                !else: particle id for particle species
      iosngl=0  !event # for which counting at metropolis iterations is done
      iorejz=0  !reject pair exchange with only zeros (1) or not (else)
      iompar=4  !parameter for windowing algorithm
      iozinc=0  !if iozevt>0: modifies iozero for testing (in sr hgcnbi)
      iozevt=0  !if >0: modifies iozero for testing (in sr hgcnbi)
      nadd=0    !number of pi0s added to minimum initial configuration
      iterma=-6 !>0: maximum number of iterations
                !<0: - number of iterations per corr time
      iterpr=10000 !iter-increment for printout
      iterpl=1  !iter-increment for plot
      iternc=50 !iterations not counted for analysis
      epsr=1e-4 !required accuracy in sr hnbraw
      keepr=1   !keep most random numbers rnoz in hnbodz (1)
                !  or update all (0)

c  strangelets

      iopenu=1      !option for minimum energy
                    !1: sum of hadron masses
                    !2: bag model curve with minimum at nonzero strangen
      themas=.51225 !parameter theta in berger/jaffe mass formula

c tests

      iotst1=0     !test
      iotst2=0     !test
      iotst3=0     !test
      iotst4=0     !test
      do i=1,mxparam
      vparam(i)=0
      enddo   

c  jpsi,qkonia (quarkonia)

      jpsi=0     !jpsi to be produced (1) or not (0)
      jpsifi=0   !jpsi final state interaction (1) or not (0)
      sigj=0.2   !jpsi nucleon cross section [fm**2]
      taumx=20   !max time for jpsi evolution
      nsttau=100 !time steps for jpsi evolution
      ijphis=0   !fill jpsi histograms (1) or not (0)
      qkonia1=0
      qkonia2=0
      qkonia3=0
      qkonia4=0

c  analysis: intermittency, space-time, droplets, formation time

      ymximi=2   !upper limit for rapidity interval for intermittency analysis
      imihis=0   !fill intermittency histograms (1) or not (0)
      isphis=0   !fill space-time histograms (1) or not (0)
      iologb=0   !0=linear bins   1=log bins
      ispall=1   !xspace: all ptls (1) or only interacting ptls (else)
      wtmini=-3  !tmin in xspace
      wtstep=1   !t-step in xspace
      iwcent=0   !only central point (1) or longitudinal distr (else) in xspace
      iclhis=0   !fill droplet histograms (1) or not (0)
      iwtime=0   !fill formation time histogram (1) or not (else)
      wtimet=100 !max time in wtime
      wtimei=0   !max mass in wtime
      wtimea=1000 !max mass in wtime


c  storing
      maxrec(1)=7
      irecty(1,1)=1
      irecty(2,1)=2
      irecty(3,1)=3
      irecty(4,1)=4
      irecty(5,1)=5
      irecty(6,1)=6
      irecty(7,1)=7
      maxrec(2)=14
      irecty(1,2)=1
      irecty(2,2)=2
      irecty(3,2)=3
      irecty(4,2)=4
      irecty(5,2)=5
      irecty(6,2)=6
      irecty(7,2)=7
      irecty(8,2)=8
      irecty(9,2)=9
      irecty(10,2)=10
      irecty(11,2)=11
      irecty(12,2)=12
      irecty(13,2)=13
      irecty(14,2)=14

c root

      ifillTree=0 !fill root tree
      iextree=1 !consider tree extensions up to iextree (0 for original format)
      ifemto=0 !analyse root tree (femto)
      ivmd=0   !analyse root tree (vmd)
      idih=0   !analyse root tree (dih)
      mixevt=5
      ifillH2=0 !fill 2d histo
      igrTree=0
      muTree=0

c  other

      gaumx=8    !range for gauss distribution
      nclean=1   !clean /cptl/ if nclean > 0
                 !(not for part with istptl<istmax if nclean=1 (do not change analysis)
                 ! for all part with ist.ne.0 if nclean > 1 (particle nb reduce at max))
      istore=0   !0: no storage to data-file
                 !-1: epos full info (fixed format)
                 !1: epos     standard
                 !2: OSC1997A standard
                 !3: OSC1999A standard
                 !4: pdg      standard
                 !5: calls routine ustore (modifiable by the user)
                 !6: HepMC standard
                 !7: Les Houches Event Format (LHEF) standard
      ioidch=1   !id choice for storage to data-file
      iframe=0   !frame specification production run
      jframe=0   !frame specification analysis
      kframe=0   !frame specification analysis (2nd frame)
                 ! 1:total
                 !11:nucleon-nucleon
                 !12:target
                 !21:gamma-nucleon
                 !22:lab
                 !32:sphericity
                 !33:thrust
      irescl=1   !momentum rescaling (1) or not (0)
      ifrade=1   !suppression of fragmentation and decay (0)
      idecay=1   !suppression of decay (0)
                 !idecay=2 : force decay in hllex even for no hydro
      ihdecay=1  !suppression of decay after hadr casc (0)
      jdecay=1   !suppression of cluster decay (0), concerns only ity=60 cluster
      ntrymx=10  !try-again parameter
      istmax=50   !analyse only istptl <= istmax
      istfor=-999 !analyse  forcing ist
      irdmpr=0   !random sign for projectile if 1
      ilprtg=1   !consider leading particle in projectile (1)
                 !or target  (-1) side
      iselect=0  !selection of hard processes
c  constants

      pi=3.1415927
      pii=1./(2*pi**2)
      hquer=0.197327
      prom=.94
      piom=.14
      ainfin=1e31

c air

      airanxs(1)=14.007
      airznxs(1)=7.
      airwnxs(1)=0.781
      airanxs(2)=15.999
      airznxs(2)=8.
      airwnxs(2)=.21
      airanxs(3)=39.948
      airznxs(3)=18.
      airwnxs(3)=0.009
      airavanxs=airanxs(1)*airwnxs(1)+airanxs(2)*airwnxs(2)
     &         +airanxs(3)*airwnxs(3)
      airavznxs=airznxs(1)*airwnxs(1)+airznxs(2)*airwnxs(2)
     &         +airznxs(3)*airwnxs(3)

c initializations

      jselect=0  
      fxsplit=1.
      idooptns=0
      ibeginwrite=0
      macrocurrent=0
      macrocurrentline=0
      macronr=0
      facpos(1)=1
      facpos(2)=1
      iSwitchSides=0
      nacc=0
      nrej=0
      ktnbod=-1
      do i=1,mamx
      jproj(1,i)=0
      jtarg(1,i)=0
      jproj(2,i)=0
      jtarg(2,i)=0
      enddo
      do i=1,8
      hdtext(i)='                              '
      enddo
      gefac=1
      istatom=0
      nbarray=0
      nfifac=1
      do npom=1,mxxpom
      nfillt(npom)=0
      enddo
      iexhd=0
      icc20=0
      mcentr=0
      mmxcentr=0
      esollxx=-1
      ranphi=0
      jjeos=0
      jjtb=1
      irootcproot=0
      iboein=0
      ieof=0
      isyst=0
      ihyskip=0
      istateos=0
      istatptl=0
      nofreeze=0
      noevt=0
      do i=1,100
      zclass(1,i)=0.
      zclass(2,i)=0.
      zclass(3,i)=0.
      zclass(4,i)=-1.
      zclass(5,i)=-1.
      enddo
      kexit=0
      do i=0,100
      zlimit(i)=0
      enddo
      do i=1,mrclass
      nrclass(i)=0
      enddo
      iopcnt=0
      ncenthy=0
      ixgeometry=0
      ixbDens=0
      ixtau=0
      iEmsB=0
      iEmsBg=0
      iEmsPm=0
      iEmsPx=0
      iEmsPDF=0
      iEmsSe=0
      iEmsDr=0
      iEmsRx=0
      iEmsI2=0
      iEmsI1=0
      iSpaceTime=0
      nemsi=0
      facmss=1.
      nstmax=0
      do 6 i=1,99
      prob(i)=0
      do 6 j=1,2
      icbac(i,j)=0
6     icfor(i,j)=0
      imsg=0
      do j=1,mxjerr
       jerr(j)=0
      enddo
      do j=1,mxcoti
      do i=1,mycoti
       coti(i,j)=0
      enddo
      enddo
      ntevt=0
      nrevt=0
      naevt=0
      nrstr=0
      nrptl=0
      nptlu=0
      do itau=1,mxtau
      volsum(itau)=0
      vo2sum(itau)=0
      nclsum(itau)=0
      do ivol=1,mxvol
      do ieps=1,mxeps
      clust(itau,ivol,ieps)=0
      enddo
      enddo
      enddo
      iutotc=0
      iutote=0
      nopen=0
      nopenr=0
      do npom=1,mxxpom
      nopent(npom)=0
      enddo
      nopend2h=0
      knxopen=0
      kchopen=0
      khiopen=0
      kdtopen=0
      khepmcopen=0
      klgopen=0
      ifdat=0
      ifncs=0
      xpar1=0.
      xpar2=0.
      xpar3=0.
      xpar4=0.
      xpar5=0.
      xpar6=0.
      xpar7=0.
      xpar8=0.
      xpar9=0.
      xpar10=0.
      xpar11=0.
      xpar12=0.
      xpar13=0.
      xpar14=0.
      xpar15=0.
      xpar16=0.
      xpar17=0.
      xpar98=0.
      xpar99=0.
      if(iop.eq.0)khisto=0
      nrclu=0
      nrnody=0
      do n=1,mxnody
      nody(n)=0
      enddo
      nrpri=0
      ctaumin=-1.
      do n=1,mxpri
      subpri(n)='      '
      ishpri(n)=0
      enddo
      nctcor=0
      ncttim=0
      do n=1,matau
      tauv(n)=0
      enddo
      ncnt=0
      nrnucl(1)=0
      nrnucl(2)=0
      do i=1,mxnucl
      rnucl(i,1)=0
      rnucl(i,2)=0
      rnuclo(i,1)=0
      rnuclo(i,2)=0
      bnucl(i,1)=0
      bnucl(i,2)=0
      bnucl(i,3)=0
      bnucl(i,4)=0
      enddo
      xbtot(1)=0.
      xbtot(2)=0.
      inicnt=0
      accept=0.
      reject=0.
      do n=1,matau
      tauv(n)=0.
      enddo
      anquasiel=0.
      nglacc=0
      ifset=1

 1002 continue ! ----> entry iop=2 <----


c  analysis

      xvaria='numptl'
      yvaria='ycmptl'
      normal=11
      xminim=-100
      xmaxim=100
      nrbins=100
      hisfac=1
      do nr=1,mxbins
      do l=1,5
      ar(nr,l)=0
      enddo
      enddo
      nozero=0
      ibmin=1
      ibmax=1e8

      return
      end

c-----------------------------------------------------------------------
      subroutine ustore
c-----------------------------------------------------------------------
c     writes the results of a simulation into the common hepevt
c     contains a description of the stored variables.
c     modifiable by the user
c-----------------------------------------------------------------------
#include "aaa.h"
      integer iepo2hep(mxptl)

c  count the number of particles to be stored (--> nptevt)

      nptevt=0
      do i=1,nptl
        iepo2hep(i)=-1    !initialize hep index to epos index
        if(istptl(i).le.istmax)nptevt=nptevt+1
      enddo

c  store event variables in HEP common :


c information available :
c     nrevt.......... event number
      nevhep=nrevt
c     nptevt ........ number of (stored!) particles per event
c     bimevt ........ absolute value of impact parameter
c     phievt ........ angle of impact parameter
c     kolevt ........ number of collisions
c     pmxevt ........ reference momentum
c     egyevt ........ pp cm energy (hadron) or string energy (lepton)
c     npjevt ........ number of primary projectile participants
c     ntgevt ........ number of primary target participants
c     npnevt ........ number of primary projectile neutron spectators
c     nppevt ........ number of primary projectile proton spectators
c     ntnevt ........ number of primary target neutron spectators
c     ntpevt ........ number of primary target proton spectators
c     jpnevt ........ number of absolute projectile neutron spectators
c     jppevt ........ number of absolute projectile proton spectators
c     jtnevt ........ number of absolute target neutron spectators
c     jtpevt ........ number of absolute target proton spectators

      nhep=0
      do i=1,nptl

      if(istptl(i).le.istmax.or.i.le.maproj+matarg)then !store events with istptl < istmax

        nhep=nhep+1
        if(nhep.gt.nmxhep)then
          print *,'Warning : produced number of particles is too high'
          print *,'          Particle list is truncated at ', nmxhep
          print *,'          CHANGE HEPEVT_EntriesAllocation in your ',
     &           'HepMC library (HEPEVT_Wrapper.h) to fix this !!!!'
          goto 1000
        endif

c  store particle variables:

c     i ............. particle number
c     idptl(i) ...... particle id
      idpdg=idtrafo('nxs','pdg',idptl(i))
      if(idpdg.ne.99)then
        idhep(nhep)=idpdg
        iepo2hep(i)=nhep
      else
        print *,'Skip particle',i,idptl(i)
        nhep=nhep-1
        goto 100
      endif
c     pptl(1,i) ..... x-component of particle momentum (GeV/c)
      phep(1,nhep)=dble(pptl(1,i))
c     pptl(2,i) ..... y-component of particle momentum (GeV/c)
      phep(2,nhep)=dble(pptl(2,i))
c     pptl(3,i) ..... z-component of particle momentum (GeV/c)
      phep(3,nhep)=dble(pptl(3,i))
c     pptl(4,i) ..... particle energy  (GeV)
      phep(4,nhep)=dble(pptl(4,i))
c     pptl(5,i) ..... particle mass    (GeV/c2)
      phep(5,nhep)=dble(pptl(5,i))
c     iorptl(i) ..... particle number of father (if .le. 0 : no father)
      if(iorptl(i).gt.0)then
        jmohep(1,nhep)=iepo2hep(iorptl(i))
      else
        jmohep(1,nhep)=-1
      endif
c     jorptl(i) ..... particle number of mother (if .le. 0 : no mother)
      if(jorptl(i).gt.0)then
        jmohep(2,nhep)=iepo2hep(jorptl(i))
      else
        jmohep(2,nhep)=-1
      endif
c     ifrptl(1,i) ..... particle number of first daughter (no daughter=0)
      jdahep(1,nhep)=0  !need a second loop to calculated proper indice
c     ifrptl(2,i) ..... particle number of last daughter (no daughter=0)
      jdahep(2,nhep)=0  !need a second loop to calculated proper indice
c     istptl(i) ..... generation flag: last gen. (0) or not (1)
      isthep(nhep)=min(2,istptl(i)+1)  !in hep:1=final, 2=decayed
      if(i.le.maproj+matarg)isthep(nhep)=4     !beam particles
c     ityptl(i) ..... particle type (string, remnant ...)
c     xorptl(1,i) ... x-component of formation point (fm)
      vhep(1,nhep)=xorptl(1,i)*1e-12 !conversion to mm
c     xorptl(2,i) ... y-component of formation point (fm)
      vhep(2,nhep)=xorptl(2,i)*1e-12 !conversion to mm
c     xorptl(3,i) ... z-component of formation point (fm)
      vhep(3,nhep)=xorptl(3,i)*1e-12 !conversion to mm
c     xorptl(4,i) ... formation time (fm/c)
      vhep(4,nhep)=xorptl(4,i)*1E-12 !conversion to mm/c
c     tivptl(1,i) ... formation time (always in the pp-cms!)
c     tivptl(2,i) ... destruction time (always in the pp-cms!)

 100   continue

      endif
      enddo

 1000 continue
c Second list to update daughter list (only if mothers are in list)
      if(istmax.ge.1)then
        nhep=0
        do i=1,nptl

          if(istptl(i).le.istmax)then !store events with istptl < istmax

            nhep=nhep+1
            if(nhep.gt.nmxhep)return

c           ifrptl(1,i) ..... particle number of first daughter (no daughter=0)
            if(ifrptl(1,i).gt.0)then
              jdahep(1,nhep)=iepo2hep(ifrptl(1,i))
            else
              jdahep(1,nhep)=0
            endif
c           ifrptl(2,i) ..... particle number of last daughter (no daughter=0)
            if(ifrptl(2,i).gt.0)then
              jdahep(2,nhep)=iepo2hep(ifrptl(2,i))
            else
              jdahep(2,nhep)=0
            endif

          endif
        enddo
      endif


      return
      end


c-----------------------------------------------------------------------
      subroutine bstora
c-----------------------------------------------------------------------
c     writes the results of a simulation into the file with unit ifdt
c     contains a description of the stored variables.
c-----------------------------------------------------------------------
#include "aaa.h"
C...User process initialization commonblock.
      INTEGER MAXPUP
      PARAMETER (MAXPUP=100)
      INTEGER IDBMUP,PDFGUP,PDFSUP,IDWTUP,NPRUP,LPRUP
      DOUBLE PRECISION EBMUP,XSECUP,XERRUP,XMAXUP
      COMMON/HEPRUP/IDBMUP(2),EBMUP(2),PDFGUP(2),PDFSUP(2),
     &IDWTUP,NPRUP,XSECUP(MAXPUP),XERRUP(MAXPUP),XMAXUP(MAXPUP),
     &LPRUP(MAXPUP)
      SAVE /HEPRUP/
      common/photrans/phoele(4),ebeam,noevt
      common/record/maxrec(2),irecty(30,2)
      character code*8,version*8,frame*4,ldum*888

      if(istore.eq.0)return

      code='EPOS   '
      if(model.eq.2)code='QGSJET01'
      if(model.eq.3)code='GHEISHA '
      if(model.eq.4)code='PYTHIA  '
      if(model.eq.5)code='HIJING  '
      if(model.eq.6)code='SIBYLL  '
      if(model.eq.7.or.model.eq.11)code='QGSJETII'
      if(model.eq.8)code='PHOJET  '
      if(model.eq.9)code='FLUKA   '
      if(model.eq.12)code='DPMJET  '
      write(version,'(f6.3,3x)')iversn/1000.+0.0001

      if(iframe.eq. 1)frame='ttcm'
      if(iframe.eq.11)frame='nncm'
      if(iframe.eq.12)frame='targ'
      if(iframe.eq.21)frame='gncm'
      if(iframe.eq.22)frame='lncm'
      ntest=1
      if (istore.eq.2) then     ! OSC1997A
        if(iappl.eq.3)then
          read(ifdt,'(A)')ldum
          read(ifdt,'(A)')ldum
          read(ifdt,'(A)')ldum
        else
        write (ifdt,'(a)') 'OSC1997A'
        write (ifdt,'(a)') 'final_id_p_x'
        write(ifdt,100) code,version
     *       ,maproj,laproj,matarg,latarg,frame,engy,ntest
 100    format(2(a8,'  '),'(',i3,',',i3,')+(',i3,',',i3,')',
     *       '  ',a4,'  ',e10.4,'  ',i8)
        maxrec(1)=4
        irecty(1,1)=1           !nevt
        irecty(2,1)=2
        irecty(3,1)=3
        irecty(4,1)=4
        maxrec(2)=11
        irecty(1,2)=1           !nr
        irecty(2,2)=2           !id
        irecty(3,2)=3           !px
        irecty(4,2)=4           !py
        irecty(5,2)=5           !pz
        irecty(6,2)=6           !E
        irecty(7,2)=7           !M
        irecty(8,2)=11          !x
        irecty(9,2)=12          !y
        irecty(10,2)=13         !z
        irecty(11,2)=14         !t
        endif
      elseif(istore.eq.3) then
        if(iappl.eq.3)then
          read(ifdt,'(A)')ldum
          read(ifdt,'(A)')ldum
          read(ifdt,'(A)')ldum
          read(ifdt,'(A)')ldum
        else
 201    format('# ',a)
        write (ifdt,201) 'OSC1999A'
        if(istmax.eq.0) then
          write (ifdt,201) 'final_id_p_x'
        elseif(istmax.ge.2) then
          write (ifdt,201) 'full_event_history'
        endif
 202    format('# ',a8,' ',a8)
        write(ifdt,202) code,version !3rd line
 203    format('# (',i3,',',i3,')+(',i3,',',i3,')',
     *       '  ',a4,'  ',e10.4,'  ',i8)
        write(ifdt,203) maproj,laproj,matarg,latarg,frame,engy,ntest
        endif
        maxrec(1)=5
        irecty(1,1)=2           !nevt
        irecty(2,1)=0           !zero
        irecty(3,1)=1           !additional information
        irecty(4,1)=3           !additional information
        irecty(5,1)=4           !additional information
        maxrec(2)=12
        irecty(1,2)=1           !nr
        irecty(2,2)=2           !id
        irecty(3,2)=10          !ist
        irecty(4,2)=3           !px
        irecty(5,2)=4           !py
        irecty(6,2)=5           !pz
        irecty(7,2)=6           !E
        irecty(8,2)=7           !M
        irecty(9,2)=11          !x
        irecty(10,2)=12         !y
        irecty(11,2)=13         !z
        irecty(12,2)=14         !t
                                ! nin nout [optional information]
                                ! ipart id ist px py pz p0 mass x y z t [optional information]
      elseif(istore.eq.7)then

C rename .data file .lhe file
      if(kdtopen.eq.1)close(ifdt)
      kdtopen=1
      fndt(nfndt-4:nfndt)=".lhe "
      nfndt=nfndt-1
      if(iappl.eq.3)then
        open(unit=ifdt,file=fndt(1:nfndt),status='old')
      else
        open(unit=ifdt,file=fndt(1:nfndt),status='unknown')
C...Write header info.
        write(ifdt,'(A)') '<LesHouchesEvents version="1.0">'
        write(ifdt,'(A)') '<!--'
        write(ifdt,'(A,A8,A,A8)') '# File generated with ',code,' '
     *                           ,version
        write(ifdt,'(A,I9)')'# Total number of min. bias events : '
     *                      ,nevent
        write(ifdt,'(A)') '# 4 types of subprocess are defined : '
        write(ifdt,'(A)')
     *  '#  ->  1 : Non Diffractive events AB-->X'
        write(ifdt,'(A)')
     *  '#  ->  2 : Double Diffractive events AB-->XX'
        write(ifdt,'(A)')
     *  '#  ->  3 : Central Diffractive events AB-->AXB'
        write(ifdt,'(A)')
     *  '#  ->  4 : Single Diffractive events AB-->XB or AB-->AX'
        write(ifdt,'(A)')
     *  '#geometry gives impact parameter (fm) and phi (rad) of events'
        write(ifdt,'(A)') '-->'

C...Set initialization info and get number of processes.
        IDBMUP(1)=idtrafo('nxs','pdg',idproj)  !projectile
        IDBMUP(2)=idtrafo('nxs','pdg',idtarg)  !target
        if(noebin.lt.0)then
        EBMUP(1)=dble(elepti)                 !energy beam proj
        EBMUP(2)=dble(ebeam)                  !energy beam targ
        PDFGUP(1)=-1d0            !PDFlib group code for proj PDF (lepton)
        PDFGUP(2)=1d0             !PDFlib group code for targ PDF (user defined)
        PDFSUP(1)=-1d0            !PDFlib set code for proj PDF (lepton)
        PDFSUP(2)=1d0             !PDFlib set code for targ PDF (user defined)
        else
        EBMUP(1)=dble(0.5*engy)                 !energy beam proj
        EBMUP(2)=dble(0.5*engy)                 !energy beam targ
        if(iappl.eq.6)then
        PDFGUP(1)=-1d0            !PDFlib group code for proj PDF (lepton)
        PDFGUP(2)=-1d0            !PDFlib group code for targ PDF (lepton)
        PDFSUP(1)=-1d0            !PDFlib set code for proj PDF (lepton)
        PDFSUP(2)=-1d0            !PDFlib set code for targ PDF (lepton)
        elseif(iappl.eq.7)then
        PDFGUP(1)=-1d0            !PDFlib group code for proj PDF (lepton)
        PDFGUP(2)=1d0             !PDFlib group code for targ PDF (user defined)
        PDFSUP(1)=-1d0            !PDFlib set code for proj PDF (lepton)
        PDFSUP(2)=1d0             !PDFlib set code for targ PDF (user defined)
        else
        PDFGUP(1)=1d0             !PDFlib group code for proj PDF (user defined)
        PDFGUP(2)=1d0             !PDFlib group code for targ PDF (user defined)
        PDFSUP(1)=1d0             !PDFlib set code for proj PDF (user defined)
        PDFSUP(2)=1d0             !PDFlib set code for targ PDF (user defined)
        endif
        endif
        IDWTUP=3                !weight=1 for all events
        NPRUP=4                 !number of subprocess (ND,DD,CD,SD)
        IPR=1                   !subprocesses (store non diffractive events)
        XSECUP(IPR)=dble(sigcut)*1d9 !cross section in pb
        XERRUP(IPR)=0d0         !statistical error
        XMAXUP(IPR)=1d0         !weight
        LPRUP(IPR)=1            !ND event (typevt=1)
        IPR=2                   !subprocesses (store double diffractive events)
        XSECUP(IPR)=dble(sigdd)*1d9 !cross section in pb
        XERRUP(IPR)=0d0         !statistical error
        XMAXUP(IPR)=1d0         !weight
        LPRUP(IPR)=2            !DD event (typevt=2)
        IPR=3                   !subprocesses (store single diffractive events)
        XSECUP(IPR)=dble(sigdif-sigdd-sigsd)*1d9 !cross section in pb
        XERRUP(IPR)=0d0         !statistical error
        XMAXUP(IPR)=1d0         !weight
        LPRUP(IPR)=3            !CD event (typevt=3)
        IPR=4                   !subprocesses (store single diffractive events)
        XSECUP(IPR)=dble(sigsd)*1d9 !cross section in pb
        XERRUP(IPR)=0d0         !statistical error
        XMAXUP(IPR)=1d0         !weight
        LPRUP(IPR)=4            !SD event (typevt=4)

C...Copy initialization lines, omitting trailing blanks.
C...Embed in <init> ... </init> block.
        write(ifdt,'(A)') '<init>'
        write(ifdt,*) IDBMUP(1),IDBMUP(2),EBMUP(1),EBMUP(2)
     &     ,PDFGUP(1),PDFGUP(2),PDFSUP(1),PDFSUP(2),IDWTUP,NPRUP
        DO 120 IPR=1,NPRUP
          write(ifdt,*) XSECUP(IPR),XERRUP(IPR),XMAXUP(IPR),LPRUP(IPR)
 120    CONTINUE
        write(ifdt,'(A)') '</init>'
       endif

      endif

      end

c-----------------------------------------------------------------------
      subroutine bstore
c-----------------------------------------------------------------------
c     writes the results of a simulation into the file with unit ifdt
c     contains a description of the stored variables.
c-----------------------------------------------------------------------

#include "aaa.h"
      common/record/maxrec(2),irecty(30,2)
      common/dimensi/k2(100)

      nptevt=0
      do n=1,nptl
        iok=1               !idcode simple
        if(istptl(n).gt.istmax
     &     .or.(ioidch.eq.2.and.idptl(n).gt.10000))then
          iok=0
        endif
      if (iok.eq.1) nptevt=nptevt+1
      enddo
 11   format (i6,' ',$)
 12   format (e12.6,' ',$)
 13   format (f3.0,' ',$)
      do i=1,maxrec(1)
        l=irecty(i,1)
        if(l.eq.0)write(ifdt,21) 0
        if(l.eq.1)write(ifdt,11)nrevt
        if(l.eq.2)write(ifdt,11)nptevt
        if(l.eq.3)write(ifdt,12)bimevt
        if(l.eq.4)write(ifdt,12)phievt
        if(l.eq.5)write(ifdt,11)kolevt
        if(l.eq.6)write(ifdt,12)pmxevt
        if(l.eq.7)write(ifdt,12)egyevt
        if(l.eq.8)write(ifdt,11)npjevt
        if(l.eq.9)write(ifdt,11)ntgevt
        if(l.eq.10)write(ifdt,11)npnevt
        if(l.eq.11)write(ifdt,11)nppevt
        if(l.eq.12)write(ifdt,11)ntnevt
        if(l.eq.13)write(ifdt,11)ntpevt
        if(l.eq.14)write(ifdt,11)jpnevt
        if(l.eq.15)write(ifdt,11)jppevt
        if(l.eq.16)write(ifdt,11)jtnevt
        if(l.eq.17)write(ifdt,11)jtpevt
        if(l.eq.20)write(ifdt,12)amproj
        if(l.eq.21)write(ifdt,12)amtarg
        if(l.eq.22)write(ifdt,12)qsqevt
        if(l.eq.23)write(ifdt,12)xbjevt
        if(l.eq.24)write(ifdt,13)typevt
      enddo
      write (ifdt,*)            !RETURN
 21   format (i6,' ',$)
 22   format (e12.6,' ',$)
 23   format (i10,' ',$)
      do n=1,nptl
        iok=1                   !idcode simple
        if(istptl(n).gt.istmax
     &     .or.(ioidch.eq.2.and.idptl(n).gt.10000))then
          iok=0
        endif
        if (iok.eq.1) then
          id=idptl(n)
          if(istore.eq.2.or.istore.eq.4.or.ioidch.eq.2)then
            id=idtrafo('nxs','pdg',idptl(n))
          endif
          do i=1,maxrec(2)
            l=irecty(i,2)
            if(l.eq.0)write(ifdt,21) 0
            if(l.eq.1)write(ifdt,21) n
            if(l.eq.2)write(ifdt,23) id
            if(l.eq.3.or.l.eq.17)write(ifdt,22) pptl(1,n)
            if(l.eq.4.or.l.eq.17)write(ifdt,22) pptl(2,n)
            if(l.eq.5.or.l.eq.17)write(ifdt,22) pptl(3,n)
            if(l.eq.6.or.l.eq.17)write(ifdt,22) pptl(4,n)
            if(l.eq.7.or.l.eq.17)write(ifdt,22) pptl(5,n)
            if(l.eq.8)write(ifdt,21) iorptl(n)
            if(l.eq.9)write(ifdt,21) jorptl(n)
            if(l.eq.10)write(ifdt,21) istptl(n)
            if(l.eq.11.or.l.eq.18)write(ifdt,22) xorptl(1,n)
            if(l.eq.12.or.l.eq.18)write(ifdt,22) xorptl(2,n)
            if(l.eq.13.or.l.eq.18)write(ifdt,22) xorptl(3,n)
            if(l.eq.14.or.l.eq.18)write(ifdt,22) xorptl(4,n)
            if(l.eq.19)write(ifdt,22) dezptl(n)
            if(l.eq.21)write(ifdt,21) ifrptl(1,n)
            if(l.eq.22)write(ifdt,21) ifrptl(2,n)
            if(l.eq.23)write(ifdt,21) ityptl(n)
            if(l.eq.15) then
              if(iorptl(n).gt.0)then
                write(ifdt,23) idptl(iorptl(n))
              else
                write(ifdt,23) 0
              endif
            endif
            if(l.eq.16) then
              if(jorptl(n).gt.0)then
                write(ifdt,23) idptl(jorptl(n))
              else
                write(ifdt,23) 0
              endif
            endif
          enddo
          write (ifdt,*)        !RETURN
        endif
      enddo
      return
      end


c-----------------------------------------------------------------------
      subroutine bread
c-----------------------------------------------------------------------
c     reads the results of a simulation into the file with unit ifdt
c     contains a description of the stored variables.
c-----------------------------------------------------------------------

#include "aaa.h"
      common/record/maxrec(2),irecty(30,2)
      character*255 line
      dimension inptl(mxptl)
      data ichkfile/0/
      save ichkfile
      logical info
C...User process event common block.
      INTEGER MAXNUP
      PARAMETER (MAXNUP=50000)  !extend array for file production
c      PARAMETER (MAXNUP=500)
      INTEGER NUP,IDPRUP,IDUP,ISTUP,MOTHUP,ICOLUP
      DOUBLE PRECISION XWGTUP,SCALUP,AQEDUP,AQCDUP,PUP,VTIMUP,SPINUP
      COMMON/HEPEUP/NUP,IDPRUP,XWGTUP,SCALUP,AQEDUP,AQCDUP,IDUP(MAXNUP),
     &ISTUP(MAXNUP),MOTHUP(2,MAXNUP),ICOLUP(2,MAXNUP),PUP(5,MAXNUP),
     &VTIMUP(MAXNUP),SPINUP(MAXNUP)
      SAVE /HEPEUP/

        if(istore.eq.-1.or.iappl.eq.-1)then

      ifinp=ifdt
      if(ichkfile.eq.0)then
        if(iappl.eq.-1)then
          inquire(file=fnin(1:nfnin),exist=info)
          if(info)then
            open(unit=ifin,file=fnin(1:nfnin),status='old')
            ifinp=ifin
          else
            call utstop('Cannot open file for conversion !&')
          endif
        endif
      endif

      read(ifinp,*,end=999)nptevt,bimevt,phievt,kolevt,pmxevt,egyevt
     *           ,npjevt,ntgevt,qsqevt,typevt
      if(nptevt.eq.0.or.nptevt.gt.mxptl)then
        print *,'sorry, '
        print *,'there is wrong particle number in the event record  '
        stop
      endif
      nptl=0
      do n=1,nptevt
        nptl=nptl+1
        read(ifinp,*)nidp,id,pp1,pp2,pp3,pp4,pp5,io,jo,is,it
     *                     ,xo1,xo2,xo3,xo4,ifr1,ifr2,dez
c keep the structure of the original event
        do while(nptl.lt.nidp)
          idptl(nptl)=0
          pptl(1,nptl)=0.
          pptl(2,nptl)=0.
          pptl(3,nptl)=0.
          pptl(4,nptl)=0.
          pptl(5,nptl)=0.
          iorptl(nptl)=0
          jorptl(nptl)=0
          istptl(nptl)=5
          ityptl(nptl)=0
          xorptl(1,nptl)=0.
          xorptl(2,nptl)=0.
          xorptl(3,nptl)=0.
          xorptl(4,nptl)=0.
          ifrptl(1,nptl)=0
          ifrptl(2,nptl)=0
          zpaptl(1,nptl)=0.
          zpaptl(2,nptl)=0.
          dezptl(nptl)=0
          nptl=nptl+1
        enddo
        idptl(nptl)=id
        pptl(1,nptl)=pp1
        pptl(2,nptl)=pp2
        pptl(3,nptl)=pp3
        pptl(4,nptl)=pp4
        pptl(5,nptl)=pp5
        iorptl(nptl)=io
        jorptl(nptl)=jo
        istptl(nptl)=is
        ityptl(nptl)=it
        xorptl(1,nptl)=xo1
        xorptl(2,nptl)=xo2
        xorptl(3,nptl)=xo3
        xorptl(4,nptl)=xo4
        ifrptl(1,nptl)=ifr1
        ifrptl(2,nptl)=ifr2
        dezptl(nptl)=dez
      enddo

        elseif(istore.eq.4)then

c skip intro
 10    read(ifdt,'(a)',end=999)line
       if(line(1:7).ne."<event>")goto 10

c      write(ifdt,'(A)') '<event>'
      read(ifdt,*)NUP,typevt,XWGTUP,SCALUP,AQEDUP,AQCDUP
      nhep=0
      nptl=0
      DO 220 i=1,nup


          nhep=nhep+1
          read(ifdt,*,end=999)IDUP(nhep),ISTUP(nhep),
     &      MOTHUP(1,nhep),MOTHUP(2,nhep),ICOLUP(1,nhep),ICOLUP(2,nhep),
     &      (PUP(J,nhep),J=1,5),VTIMUP(nhep),SPINUP(nhep)

          id=idtrafo('pdg','nxs',IDUP(nhep))
          if(id.eq.99)id=0   !unknown particle
          nptl=nptl+1
          idptl(nptl)=id
          if(ISTUP(nhep).eq.-9)then
            istptl(nptl)=1
            iorptl(nptl)=-1
            jorptl(nptl)=0
          else
            istptl(nptl)=ISTUP(nhep)-1
            iorptl(nptl)=MOTHUP(1,nhep)
            jorptl(nptl)=max(0,MOTHUP(2,nhep))
          endif
          do J=1,5                !particle momentum (GeV/c)
            pptl(J,nptl)=sngl(PUP(J,nhep))
          enddo

  220 CONTINUE

c optional informations
       read(ifdt,*,end=999)line,bimevt,phievt

       read(ifdt,*,end=999)line

c       read(ifdt,*,end=999)nptevt,bimevt,phievt,kolevt,pmxevt,egyevt
c     *           ,npjevt,ntgevt,qsqevt,typevt
c       if(nptevt.eq.0.or.nptevt.gt.mxptl)then
c         print *,'sorry, '
c         print *,'there is wrong particle number in the event record  '
c         stop
c       endif
c       do n=1,nptevt
c       read(ifdt,*)nidp,idptl(n),pptl(1,n),pptl(2,n),pptl(3,n),pptl(4,n)
c     *            ,pptl(5,n),iorptl(n),jorptl(n),istptl(n),ityptl(n)
c     *            ,xorptl(1,n),xorptl(2,n),xorptl(3,n),xorptl(4,n)
c       inptl(nidp)=n
c       if(iorptl(n).gt.0)iorptl(n)=inptl(iorptl(n))
c       if(jorptl(n).gt.0)jorptl(n)=inptl(jorptl(n))
c       enddo


        elseif(istore.eq.5)then

       read(ifdt,*,end=999)nptevt,bimevt,phievt,kolevt,pmxevt,egyevt
     *           ,npjevt,ntgevt,qsqevt,typevt
       if(nptevt.eq.0.or.nptevt.gt.mxptl)then
         print *,'sorry, '
         print *,'there is wrong particle number in the event record  '
         stop
       endif
       do n=1,nptevt
       read(ifdt,*)nidp,idptl(n),pptl(1,n),pptl(2,n),pptl(3,n),pptl(4,n)
     *            ,pptl(5,n),iorptl(n),jorptl(n),istptl(n),ityptl(n)
     *            ,xorptl(1,n),xorptl(2,n),xorptl(3,n),xorptl(4,n)
       inptl(nidp)=n
       if(iorptl(n).gt.0)iorptl(n)=inptl(iorptl(n))
       if(jorptl(n).gt.0)jorptl(n)=inptl(jorptl(n))
       enddo
       nptl=nptevt

        else

      info=.false.
      do n=1,mxptl
        inptl(n)=0
      enddo

      read(ifdt,'(a255)',end=1)line
 1    k=1
      nptevt=0
      do i=1,maxrec(1)
        l=irecty(i,1)
        if(l.eq.0)read(line(k:),'(i6)')ldummy !0
        if(l.eq.0) k=k+7
        if(l.eq.1)read(line(k:),'(i6)')ldummy !nrevt
        if(l.eq.1) k=k+7
        if(l.eq.2)read(line(k:),'(i6)')nptevt
        if(l.eq.2) k=k+7
        if(l.eq.3)read(line(k:),'(e12.6)')bimevt
        if(l.eq.3) k=k+13
        if(l.eq.4)read(line(k:),'(e12.6)')phievt
        if(l.eq.4) k=k+13
        if(l.eq.5)read(line(k:),'(i6)')kolevt
        if(l.eq.5) k=k+7
        if(l.eq.6)read(line(k:),'(e12.6)')pmxevt
        if(l.eq.6) k=k+13
        if(l.eq.7)read(line(k:),'(e12.6)')egyevt
        if(l.eq.7) k=k+13
        if(l.eq.8)read(line(k:),'(i6)')npjevt
        if(l.eq.8) k=k+7
        if(l.eq.9)read(line(k:),'(i6)')ntgevt
        if(l.eq.9) k=k+7
        if(l.eq.10)read(line(k:),'(i6)')npnevt
        if(l.eq.10) k=k+7
        if(l.eq.11)read(line(k:),'(i6)')nppevt
        if(l.eq.11) k=k+7
        if(l.eq.12)read(line(k:),'(i6)')ntnevt
        if(l.eq.12) k=k+7
        if(l.eq.13)read(line(k:),'(i6)')ntpevt
        if(l.eq.13) k=k+7
        if(l.eq.14)read(line(k:),'(i6)')jpnevt
        if(l.eq.14) k=k+7
        if(l.eq.15)read(line(k:),'(i6)')jppevt
        if(l.eq.15) k=k+7
        if(l.eq.16)read(line(k:),'(i6)')jtnevt
        if(l.eq.16) k=k+7
        if(l.eq.17)read(line(k:),'(i6)')jtpevt
        if(l.eq.17) k=k+7
        if(l.eq.20)read(line(k:),'(e12.6)')amproj
        if(l.eq.20) k=k+13
        if(l.eq.21)read(line(k:),'(e12.6)')amtarg
        if(l.eq.21) k=k+13
        if(l.eq.22)read(line(k:),'(e12.6)')qsqevt
        if(l.eq.22) k=k+13
        if(l.eq.23)read(line(k:),'(e12.6)')xbjevt
        if(l.eq.23) k=k+13
        if(l.eq.24)read(line(k:),'(f3.0)')typevt
        if(l.eq.24) k=k+4
      enddo
      if(nptevt.eq.0)then
        print *,'sorry, '
        print *,'there is no particle number in the event record  '
        stop
      endif
      do n=1,nptevt
        read(ifdt,'(a255)',end=2)line
 2      k=1
        do i=1,maxrec(2)
          l=irecty(i,2)
          if(l.eq.0)read(line(k:),'(i6)') ldummy
          if(l.eq.0) k=k+7
          if(l.eq.1)then
            read(line(k:),'(i6)') nidp
            if(nidp.gt.0.and.nidp.le.mxptl)then
              info=.true.
              inptl(nidp)=n
            endif
            k=k+7
          endif
          if(l.eq.2)read(line(k:),'(i10)') idptl(n)
          if(l.eq.2) k=k+11
          if(l.eq.3.or.l.eq.17)read(line(k:),'(e12.6)') pptl(1,n)
          if(l.eq.3.or.l.eq.17) k=k+13
          if(l.eq.4.or.l.eq.17)read(line(k:),'(e12.6)') pptl(2,n)
          if(l.eq.4.or.l.eq.17) k=k+13
          if(l.eq.5.or.l.eq.17)read(line(k:),'(e12.6)') pptl(3,n)
          if(l.eq.5.or.l.eq.17) k=k+13
          if(l.eq.6.or.l.eq.17)read(line(k:),'(e12.6)') pptl(4,n)
          if(l.eq.6.or.l.eq.17) k=k+13
          if(l.eq.7.or.l.eq.17)read(line(k:),'(e12.6)') pptl(5,n)
          if(l.eq.7.or.l.eq.17) k=k+13
          if(l.eq.8)then
            read(line(k:),'(i6)') iorptl(n)
            k=k+7
            if(info.and.iorptl(n).gt.0)iorptl(n)=inptl(iorptl(n))
          endif
          if(l.eq.9)then
            read(line(k:),'(i6)') jorptl(n)
            k=k+7
            if(info.and.jorptl(n).gt.0)jorptl(n)=inptl(jorptl(n))
          endif
          if(l.eq.10)read(line(k:),'(i6)') istptl(n)
          if(l.eq.10) k=k+7
          if(l.eq.11.or.l.eq.18)read(line(k:),'(e12.6)')xorptl(1,n)
          if(l.eq.11.or.l.eq.18) k=k+13
          if(l.eq.12.or.l.eq.18)read(line(k:),'(e12.6)')xorptl(2,n)
          if(l.eq.12.or.l.eq.18) k=k+13
          if(l.eq.13.or.l.eq.18)read(line(k:),'(e12.6)')xorptl(3,n)
          if(l.eq.13.or.l.eq.18) k=k+13
          if(l.eq.14.or.l.eq.18)read(line(k:),'(e12.6)')xorptl(4,n)
          if(l.eq.14.or.l.eq.18) k=k+13
c     if(i.eq.15)read(line(k:),'(i6)') idiptl(n)
          if(l.eq.15) k=k+7
c     if(i.eq.16)read(line(k:),'(i6)') idjptl(n)
          if(l.eq.16) k=k+7
          if(l.eq.19)read(line(k:),'(e12.6)') dezptl(n)
          if(l.eq.19) k=k+13
c          if(l.eq.21)read(line(k:),'(I6)') ifrptl(1,n)
          if(l.eq.21) k=k+7
c          if(l.eq.22)read(line(k:),'(I6)') ifrptl(2,n)
          if(l.eq.22) k=k+7
          if(l.eq.23)read(line(k:),'(I6)') ityptl(n)
          if(l.eq.23) k=k+7
        enddo
      enddo

      nptl=nptevt

        endif

      nevt=1
 999  continue
      end

c-----------------------------------------------------------------------
      subroutine aafinal
c-----------------------------------------------------------------------
c  * calculates xorptl(j,i), representing formation points.
c    (xorptl(j,i),j=1,4) is the 4-vector representing the space-time of
c    creation of particle i.
c-----------------------------------------------------------------------
#include "aaa.h"
      do iloo=1,nptl

        i=iloo
        if(iloo.gt.mxptl)then
          call restorecccptl(iloo,mxptl+2)
          i=mxptl+2
        endif

       ! projectile and target nucleons
       if(maproj.gt.0.and.matarg.gt.0)then
         if(i.le.maproj+matarg)then
           tivptl(1,i)=0
         endif
       endif

       if(idptl(i).ne.0.and.istptl(i).le.1)then
          if(    abs(tivptl(1,i)).le.ainfin
     .    .and.abs(xorptl(1,i)).le.ainfin
     .    .and.abs(xorptl(2,i)).le.ainfin
     .    .and.abs(xorptl(3,i)).le.ainfin
     .    .and.abs(xorptl(4,i)).le.ainfin
     .    .and.pptl(5,i).le.ainfin
     .    .and.pptl(4,i).gt.0.)then
c            if(ish.ge.4)call alistc('afinal&',i,i)
            t=tivptl(1,i)
            xorptl(1,i)=xorptl(1,i)+pptl(1,i)/pptl(4,i)*(t-xorptl(4,i))
            xorptl(2,i)=xorptl(2,i)+pptl(2,i)/pptl(4,i)*(t-xorptl(4,i))
            xorptl(3,i)=xorptl(3,i)+pptl(3,i)/pptl(4,i)*(t-xorptl(4,i))
            xorptl(4,i)=t
          else
            if(ish.ge.1)then
              if(iorptl(i).gt.0)idior=idptl(iorptl(i))
              !write(ifmt,'(a)')
              !.        '*** warning (afinal see check file): '
              write(ifch,'(a,i6,i10,i10,i3,1x,7(e7.1,1x))')
     .        '*** warning (afinal): ',
     .        i,idptl(i),idior,ityptl(i),tivptl(1,i), pptl(4,i)
     .        ,pptl(5,i),xorptl(1,i),xorptl(2,i),xorptl(3,i),xorptl(4,i)
            endif
            tivptl(1,i)=2*ainfin
            tivptl(2,i)=2*ainfin
            xorptl(1,i)=2*ainfin
            xorptl(2,i)=2*ainfin
            xorptl(3,i)=2*ainfin
            xorptl(4,i)=2*ainfin
          endif
        endif
        if(iloo.gt.mxptl)call dumpcccptl(mxptl+2,iloo)
      enddo
      end

c-----------------------------------------------------------------------
      subroutine afinal
c-----------------------------------------------------------------------
c  does some final calculations, to be called before call aasto.
c  * calculates nptlu, the maximum nptl for all events.
c  * in case of mod(iframe,10) .ne. 1, these vectors are transformed
c    (being originally in the "natural frame",
c    NB : boost of coordinates only if not a non-sense (otherwise put to inf)
c         always boost of momentum if possible (if not STOP !)
c  * calculates numbers of spectators:
c    npnevt (number of primary proj neutron spectators)
c    nppevt (number of primary proj proton spectators)
c    ntnevt (number of primary targ neutron spectators)
c    ntpevt (number of primary targ proton spectators)
c    jpnevt (number of absolute proj neutron spectators)
c    jppevt (number of absolute proj proton spectators)
c    jtnevt (number of absolute targ neutron spectators)
c    jtpevt (number of absolute targ proton spectators)
c-----------------------------------------------------------------------

#include "aaa.h"
      common/geom/rmproj,rmtarg,bmax,bkmx
      double precision pgampr,rgampr
      common/cgampr/pgampr(5),rgampr(4)
      common/cgbyjmax/gbyjmax /cicentrality/icentrality

      double precision pp1,pp2,pp3,pp4,pp5,om1intbc
      logical lclean
      call utpri('afinal',ish,ishini,4)

      nptlu=max0(nptl,nptlu)
      lclean=.false.
      if(mod(iframe,10).ne.1)then
        if(iframe.eq.12.or.iframe.eq.22)then    !targ
          pp1=0d0
          pp2=0d0
          pp3=dsinh(dble(yhaha))
          pp4=dcosh(dble(yhaha))
          pp5=1d0
        else
          stop'transformation not yet defined'
        endif
      endif


      do iloo=1,nptl


      i=iloo
      if(iloo.gt.mxptl)then
        call restorecccptl(iloo,mxptl+2)
        i=mxptl+2
      endif


      if(idptl(i).ne.0.and.istptl(i).le.1)then  !~~~~~~~~~~~

        if(pptl(5,i).le.ainfin
     .     .and.pptl(4,i).gt.0.)then

          if(    abs(tivptl(1,i)).le.ainfin
     .      .and.abs(xorptl(1,i)).le.ainfin
     .      .and.abs(xorptl(2,i)).le.ainfin
     .      .and.abs(xorptl(3,i)).le.ainfin
     .      .and.abs(xorptl(4,i)).le.ainfin)then

c Space-time boost
            if(mod(iframe,10).ne.1)then

              if(iframe.eq.12)then
                call utlob4(-1,pp1,pp2,pp3,pp4,pp5
     .               ,xorptl(1,i),xorptl(2,i),xorptl(3,i),xorptl(4,i))
              elseif(iframe.eq.22)then
c not the electron in lab frame in fake DIS
                if(.not.((abs(iappl).eq.1.or.iappl.eq.3)
     *             .and.i.eq.2*(maproj+matarg)+1))then
c put particle from cms to target frame
                  call utlob4(-1,pp1,pp2,pp3,pp4,pp5
     .                 ,xorptl(1,i),xorptl(2,i),xorptl(3,i),xorptl(4,i))
c do rotation of gamma in proton rest frame
                  call utrot4(-1,rgampr(1),rgampr(2),rgampr(3)
     .                 ,xorptl(1,i),xorptl(2,i),xorptl(3,i))
c boost in lab frame
                  call utlob4(-1,pgampr(1),pgampr(2),pgampr(3),pgampr(4)
     .            ,pgampr(5),xorptl(1,i),xorptl(2,i),xorptl(3,i)
     .            ,xorptl(4,i))
                endif
              else
                stop'transformation not yet defined'
              endif
            endif
          else
            tivptl(1,i)=ainfin
            xorptl(1,i)=ainfin
            xorptl(2,i)=ainfin
            xorptl(3,i)=ainfin
            xorptl(4,i)=ainfin
          endif

c Momentum boost
          if(mod(iframe,10).ne.1)then
            if(iframe.eq.12)then
              call utlob5(-yhaha
     .        , pptl(1,i), pptl(2,i), pptl(3,i), pptl(4,i), pptl(5,i))
            elseif(iframe.eq.22)then
c not the electron in lab frame in fake DIS
              if(.not.((abs(iappl).eq.1.or.iappl.eq.3)
     *           .and.i.eq.2*(maproj+matarg)+1))then
c put particle from cms to target frame
                call utlob5(-yhaha
     .          ,pptl(1,i), pptl(2,i), pptl(3,i), pptl(4,i), pptl(5,i))
c do rotation of gamma in proton rest frame
                call utrot4(-1,rgampr(1),rgampr(2),rgampr(3)
     .               , pptl(1,i), pptl(2,i), pptl(3,i))
c boost in lab frame
                call utlob4(-1,pgampr(1),pgampr(2),pgampr(3),pgampr(4)
     .         ,pgampr(5), pptl(1,i), pptl(2,i), pptl(3,i), pptl(4,i))
              endif
            endif
          endif
        elseif(model.eq.6)then
          lclean=.true.
          istptl(i)=99
        else
          call alist('list before stop in afinal&',1,nptl)
          write(ifmt,'(a,$)')'ERROR afinal: Negative energy, see check;'
          write(ifmt,*)'  i id p4 p5 = ',i,idptl(i),pptl(4,i),pptl(5,i)
          stop
        endif

      endif !~~~~~~~~~~~~~~~~~~


      if(iloo.gt.mxptl)call dumpcccptl(mxptl+2,iloo)

      enddo

      if(lclean)then
        nptl0=nptl
        call utclea(maproj+matarg+1,nptl0)
      endif


      if(ish.ge.2)then
        !if(model.eq.1)call alistf('EPOS&')
        if(model.eq.2)call alistf('QGSJET01&')
        if(model.eq.3)call alistf('GHEISHA&')
        if(model.eq.4)call alistf('PYTHIA&')
        if(model.eq.5)call alistf('HIJING&')
        if(model.eq.6)call alistf('SIBYLL 2.1&')
        if(model.eq.7.or.model.eq.11)call alistf('QGSJET II&')
        if(model.eq.8)call alistf('PHOJET&')
        if(model.eq.9)call alistf('FLUKA&')
        if(model.eq.10)call alistf('URQMD&')
        if(model.eq.12)call alistf('DPMJET&')
      endif

c      if(isto.eq.1)stop
c$$$      call testconex(2)

      npnevt=0
      nppevt=0
      ntnevt=0
      ntpevt=0
      jpnevt=0
      jppevt=0
      jtnevt=0
      jtpevt=0
      if(ish.ge.6)write(ifch,'(/31a1/a/31a1)')('-',l=1,31)
     *,'primary and absolute spectators',('-',l=1,31)
      if(ish.ge.6)write(ifch,'(/a//a/)')'projectile nucleons:'
     *,'     i    id   ior   ist'
      do i=1,maproj
      if(ish.ge.6)write(ifch,'(4i6)')i,idptl(i),iorptl(i),istptl(i)
      io=iorptl(i)
      id=idptl(i)
      is=istptl(i)
      if(io.eq.0.and.id.eq.1220)npnevt=npnevt+1
      if(io.eq.0.and.id.eq.1120)nppevt=nppevt+1
      if(io.eq.0.and.is.eq.0.and.id.eq.1220)jpnevt=jpnevt+1
      if(io.eq.0.and.is.eq.0.and.id.eq.1120)jppevt=jppevt+1
      enddo
      if(ish.ge.6)write(ifch,'(/a//a/)')'target nucleons:'
     *,'     i    id   ior   ist'
      do i=maproj+1,maproj+matarg
      if(ish.ge.6)write(ifch,'(4i6)')i,idptl(i),iorptl(i),istptl(i)
      io=iorptl(i)
      id=idptl(i)
      is=istptl(i)
      if(io.eq.0.and.id.eq.1220)ntnevt=ntnevt+1
      if(io.eq.0.and.id.eq.1120)ntpevt=ntpevt+1
      if(io.eq.0.and.is.eq.0.and.id.eq.1220)jtnevt=jtnevt+1
      if(io.eq.0.and.is.eq.0.and.id.eq.1120)jtpevt=jtpevt+1
      enddo
      if(ish.ge.6)then
      write(ifch,'(/a/)')'numbers of participants and spectators:'
      write(ifch,'(a,i4,a,i4)')'primary participants:   projectile:'
     *,npjevt,'   target:',ntgevt
      write(ifch,'(a,i4,a,i4)')'primary spectators:     projectile:'
     *,npnevt+nppevt,'   target:',ntnevt+ntpevt
      write(ifch,'(a,i4,a,i4)')
     *'primary spectator neutrons:   projectile:',npnevt
     *,'   target:',ntnevt
      write(ifch,'(a,i4,a,i4)')
     *'primary spectator protons:    projectile:',nppevt
     *,'   target:',ntpevt
      write(ifch,'(a,i4,a,i4)')'absolute spectators:    projectile:'
     *,jpnevt+jppevt,'   target:',jtnevt+jtpevt
      endif

c cross-section calculation (to be done before xana)
      if(ntevt.gt.0)then
        b1=bminim
        b2=min(bmax,bmaxim)
        a=pi*(b2**2-b1**2)
        if(iappl.eq.3.or.iappl.eq.-1)then      !read
          ntevt=nint(float(nevent)/sigine*a*10.)
          anintine=float(nevent)
          anintdiff=anintine*sigdif/sigine
          anintsdif=anintine*sigsd/sigine
        endif
        sigineex=anintine/float(ntevt)*a*10
        sigdifex=anintdiff/float(ntevt)*a*10
        sigsdex=anintsdif/float(ntevt)*a*10
        sigddex=anintddif/float(ntevt)*a*10
        if(iokoll.lt.0)then
          sigineex=sigineex*sngl(om1intbc(bb))/float(abs(iokoll))
        endif
      endif

      if(imihis.eq.1)call wimi
      if(imihis.eq.1.and.nrevt.eq.nevent)call wimino
      if(isphis.eq.1)call xspace(1)
      if(iclhis.eq.1)call wclu
      if(iclhis.eq.1.and.nrevt.eq.nevent)call wclufi
      if(iwtime.eq.1)call wtime(1)
      if(iwtime.eq.1.and.nrevt.eq.nevent)call wtime(2)

      if(ish.ge.8)call alistc('afinal&',1,nptl)

      call utprix('afinal',ish,ishini,4)
      return
      end

c-----------------------------------------------------------------------
      subroutine xjerr
c-----------------------------------------------------------------------
#include "aaa.h"
      common/ciexhd/iexhd
      character*56 cj(mxjerr)

      !-----------------------------------------------------------------
      cj(1) ='ems.f - sr ProSeF -  > 9 quarks per flavor attempted    '
      cj(3) ='ind.f - sr jintpo - finish cluster storage - p52.lt.0d0 '
      cj(5) ='lea.f - sr ProReF - decay droplet - Unsuccessful decay  '
      cj(6) ='lea.f - sr ProReF - remove hadrons - irmdropx.eq.irmdrop'
      cj(7) ='ems.f - sr emsaa - Treat hard Pomerons - ivi.lt.0       '
      cj(8) ='ems.f - sr emsaa - Diffractive Pt - iret.ne.0 (ProDiPt) '
      cj(9) ='dky.f - sr decayall - iret.eq.1 (hdecas)                '
      cj(10)='u.f   - sr uinitial - abs(ityptmp).gt.1000              '
      cj(11)='f.f   - sr FoSur - e.ge.oEeos                           '
      cj(12)='f.f   - sr FoSur - e.ge.oEeos                           '
      cj(13)='f.f   - sr EpoF - ey.ge.oEeos                           '
      cj(14)='h.f   - sr hlle - e.ge.20*oEeos                         '
      cj(15)='h.f   - sr TransferHlle - e.ge.20*oEeos                 '
      cj(16)='                                                        '
      cj(17)='                                                        '
      cj(18)='ico.f - sr IcoStr - pcore(4) .NE. qcore(4)              '
      cj(20)='f.f   - sr EpoF -  mu > m => BE condensate              '
      cj(22)='f.f   - sr FoSur - FoSur: nrad limit reached            '
      !-----------------------------------------------------------------
      if(iexhd.eq.1)then
      write(ifhi,'(a)')   '!--------------------------------'
      write(ifhi,'(a)')   '!          jerr                  '
      write(ifhi,'(a)')   '!--------------------------------'
      write(ifhi,'(a,i2)')  'openhisto name jerr  xrange 0 ',mxjerr
      write(ifhi,'(a)')  'htyp prs xmod lin ymod log yrange 0.5e-4 1.4'
      write(ifhi,'(a)')  'txt  "title  error summary "'
      write(ifhi,'(a)')  'txt  "xaxis  idx "'
      write(ifhi,'(a)')  'txt  "yaxis  jerr(idx) "'
      write(ifhi,'(a)')  'histoweight 1'
      write(ifhi,'(a)')  'array 2'
      do i=1,mxjerr
        yref=nevent
        if(i.eq.2.or.i.eq.4.or.i.eq.19.or.i.eq.21)yref=0.
        if(i.eq.3.or.i.eq.5.or.i.eq.20.or.i.eq.22)yref=jerr(i-1)
        if(i.eq.9.or.i.eq.10)yref=nptl
        if(i.eq.18)yref=nfull
        y=0
        if(yref.gt.0.)y=jerr(i)/yref
        write(ifhi,'(2e11.3)')1.*i, min(y,1.0)
        if(y.gt.0.01)write(ifmt,'(a,i2,3a,f5.2)')'WARNING(jerr:',i
     .       ,'): ',cj(i),' --> y = ',y
      enddo
      write(ifhi,'(a)') 'endarray'
      write(ifhi,'(a)') 'closehisto plot 0'

      do kk=1,3
      write(ifhi,'(a)')   '!--------------------------------'
      write(ifhi,'(a)')   '!  coti (computing time)         '
      write(ifhi,'(a)')   '!--------------------------------'
      write(ifhi,'(a,i1)')  'openhisto name coti',kk,' xrange 0 21'
      write(ifhi,'(a)')  'htyp lin xmod lin ymod log '
      write(ifhi,'(a)')  'txt  "xaxis  iopcnt "'
      write(ifhi,'(a)')  'txt  "yaxis  coti "'
      write(ifhi,'(a)')  'histoweight 0'
      write(ifhi,'(a)')  'array 3'
      yav=0
      jj= 1+mod(iopcnt-1,20)
      z=0.
      do j=1,20
        wgt=0
        if(j.eq.jj)wgt=1
        W=1
        if(zclass(3,j).gt.1e-5)W=zclass(3,j)
        if(kk.eq.1)z=coti(1,j) !hydro events
        if(kk.eq.2)z=1 !coti(4,j) !all events
        if(kk.eq.3)z=1 !coti(4,j) !all events
        y=0
        if(z.gt.0.)then
          if(kk.eq.1)y=coti(2,j)/z /60.
          if(kk.eq.2)y=coti(2,j)/z /60.
          if(kk.eq.3)y=coti(3,j)/z /60.
        endif
        yav=yav+y/20
        write(ifhi,'(3e11.3)')1.*j, y , wgt 
      enddo
      write(ifhi,'(a)') 'endarray'
      write(ifhi,'(a)')'txt  "title CPU (min) av= $sum2nbin "'
      if(kk.eq.1)write(ifhi,'(a)')'text 0.2 0.15 "CPU Hy / Hy Evt  "'
      if(kk.eq.2)write(ifhi,'(a)')'text 0.2 0.15 "CPU Hy / run   "'
      if(kk.eq.3)write(ifhi,'(a)')'text 0.2 0.15 "CPU All / run "'
      write(ifhi,'(a)') 'closehisto plot 0'
      enddo

      endif
      end

c-----------------------------------------------------------------------
      subroutine stopall(text)
c-----------------------------------------------------------------------
#include "aaa.h"
      character  text*(*) 
      call treeclose
      imax=index(text,'&')
      if(imax.gt.1)then
      write(ifch,'(/1x,72a1/1x,a,a/1x,72a1)')
     *('#',k=1,72),'############  STOP ', text(1:imax-1),('#',k=1,72)
      stop
      endif
      end

c-----------------------------------------------------------------------
      subroutine bfinal
c-----------------------------------------------------------------------
#include "aaa.h"
      double precision iutime(5)
      common/ifhq/ihq
      parameter(maxit=500000)
      common/count/nacc,nrej,naccit(maxit),nptot,npit(maxit)

      call d2hclose
      call treeclose
      call clop(3)
      if(noebin.lt.0)call phoGPHERAepo(2)  !statistic for fake DIS
      if(igrTree.gt.0.and.muTree.gt.0)then
       write(ifmt,'(a)')'igrTree,muTree>0 -> no second read -> stop'
       stop
      endif

      if(iopcnt.gt.0)then
        call timer(iutime)
        y=iutime(3)-cotid(1) + (iutime(4)-cotid(2))*1d-3
        j= 1+mod(iopcnt-1,20)
        coti(3,j)= coti(3,j) + y
        coti(4,j)= coti(4,j) + 1
      endif

#if !__BS__ && !__TP__
      if(ihq.eq.1)call aacharmclose
#endif
      if(nacc+nrej.gt.0)
     .write(ifmt,'(a,e9.3)')
     .'acceptance rate = ',float(nacc)/float(nacc+nrej)
      write(ifmt,'(a,5i8)')'nPomSoft nPomSat nPomHard All'
     .   ,nint(vparam(80))
     .,nint(vparam(81)),nint(vparam(84))
     .,nint(vparam(87)),nint(vparam(90))
      write(ifmt,'(a,8x,4i8)')'                          gg '
     .,nint(vparam(82)),nint(vparam(85))
     .,nint(vparam(88)),nint(vparam(91))
      write(ifmt,'(a,8x,5i8)')'                          ggc'
     .,nint(vparam(83)),nint(vparam(86))
     .,nint(vparam(89)),nint(vparam(92)),nint(vparam(93))
      write(ifmt,'(a,5x,4f8.3)')'                         ggc/gg '
     .,vparam(83)/max(1.,vparam(82))
     .,vparam(86)/max(1.,vparam(85))
     .,vparam(89)/max(1.,vparam(88))
     .,vparam(92)/max(1.,vparam(91))

      end

c-----------------------------------------------------------------------
      subroutine ainit
c-----------------------------------------------------------------------
#include "aaa.h"
#include "ems.h"
#include "ho.h"
#include "ico.h"
#include "sem.h"
#include "par.h"
#include "tab.h"
      parameter (nptj=129)
      common /cptj/xptj(nptj),qptj(nptj),wptj(nptj)
      common/geom/rmproj,rmtarg,bmax,bkmx
      double precision tpro,zpro,ttar,ztar,ttaus,detap,detat!,seedp
      common/cttaus/tpro,zpro,ttar,ztar,ttaus,detap,detat /ctain/mtain
      double precision rcproj,rctarg
      common/geom1/rcproj,rctarg
      common/photrans/phoele(4),ebeam,noevt
      common/cisk/iskmin,isk2min
      common/cicentrality/icentrality
      !external sptj

      call utpri('ainit ',ish,ishini,4)

#if __EB__
        if(iphsd.eq.0)stop'ERROR: make without XBS needed' 
#else
        if(igrTree.eq.0)then !No Tree read
          if(iphsd.ge.1)stop'ERROR: make XBS needed' 
        else
          nptlpt=0
        endif
#endif


      if(ifappl().eq.0)stop'\n\n STOP in ainit: iappl=0 \n\n'

      if(inicnt.eq.0.or.inicnt.eq.1)then
        inicnt=inicnt+1
      endif

      if(inicnt.eq.1)then
      call initEoL
      call idresi
      call itsyminit
      call idmass(1,qumass)
      qmass(1)=qumass !u quark effective mass (for pt distribtions)
      isospin(1)=1
      call idmass(2,qdmass)
      qmass(2)=qdmass !d quark effective mass (for pt distribtions)
      isospin(2)=-1
      call idmass(3,qsmass)
      qmass(3)=qsmass !s quark effective mass (for pt distribtions)
      isospin(3)=0
      call idmass(4,qcmass)
      !if(qcmass.ne.1.20)then !cKW22 This was done in order to keep the KWn tables
      !  qcmass=1.20
      !  write(ifmt,'(a,f5.2)')'WARNING qcmass reset to',qcmass
      !  nl=nlidtbl(4)       !idmass should give the same value than qcmass
      !  amtbl(nl)=qcmass
      !endif
      qmass(4)=qcmass !c quark effective mass (for pt distribtions)
      isospin(4)=0
      call idmass(5,qbmass)
      !if(qbmass.ne.4.60)then
      !  qbmass=4.60
      !  write(ifmt,'(a,f5.2)')'WARNING qbmass reset to',qbmass
      !  nl=nlidtbl(5)       !idmass should give the same value than qbmass
      !  amtbl(nl)=qbmass
      !endif
      qmass(5)=qbmass !b quark effective mass (for pt distribtions)
      isospin(5)=0
      call idmass(6,qtmass)
      qmass(6)=qtmass !t quark effective mass (for pt distribtions)
      isospin(6)=0
      endif

      if(igrTree.gt.0)then
      write(ifmt,'(a)')'skip ainit'
      goto 9999 !read from root
      endif

      if(inicnt.eq.1)then
        write(ifmt,'(a)')'initializations ...'
        call readcharm
        call iniq2mn
        call iniXsection
#if __ROOT__
        call eostabs(1,0,0)
#endif
      endif

#if __ROOT__
      call iniEos
#endif
      call iniRan
      if(model.ne.1.and.inicnt.eq.1)then
        if(model.eq.2)iversn=1000 !'QGSJET01'
        if(model.eq.3)iversn=1000 !'GHEISHA '
        if(model.eq.4)iversn=6110 !'PYTHIA  '
        if(model.eq.5)iversn=1380 !'HIJING   '
        if(model.eq.6)iversn=2100 !'SIBYLL  '
        if(model.eq.7)iversn=4000 !'QGSJETII-04'
        if(model.eq.8)iversn=1120 !'PHOJET  '
        if(model.eq.9)iversn=201125 !'FLUKA   '
        if(model.eq.11)iversn=3000 !'QGSJETII-03'
        if(model.eq.12)iversn=3060 !'DPMJET  '
        if(model.ne.1)iverso=iversn
        call IniModel(model)
      endif

      if(isphis.eq.1)iframe=11  !nncm
      if(icinpu.ge.1)elepti=engy
      !ctp060829  if(iopenu.eq.2)call smassi(themas)
      if(iopenu.eq.2.and.ish.eq.19)stop'change this?????????' !call smassp

      if(iappl.eq.5)then
      yhaha=0
      ypjtl=0
      endif

      if(iokoll.ne.0)then
        if(maproj.eq.1.and.matarg.eq.1)then
          maproj=abs(iokoll)
          matarg=abs(iokoll)
        else
          stop'invalid use of iokoll with maproj or matarg>1'
        endif
      endif

      if(ispherio.ne.0)jdecay=0
      if(ihacas.ne.0)ndecay=1
      if(ihacas.ne.0)idecay=0
      if(ifrade.eq.0)irescl=0
      idtarg=idtargin
      idproj=idprojin
      if(noebin.gt.1)then
        engy=-1
        ekin=-1
        if(iologe.eq.1)engy=
     *       engmin*(engmax/engmin)**((real(nrebin)-0.5)/noebin)
        if(iologe.eq.0.or.(iologe.lt.0.and.iologl.lt.0))engy=
     *       engmin+(engmax-engmin)/noebin*(nrebin-0.5)
        if(iologl.eq.1)ekin=
     *       engmin*(engmax/engmin)**((real(nrebin-0.5))/noebin)
        if(iologl.eq.0)ekin=
     *       engmin+(engmax-engmin)/noebin*(real(nrebin)-0.5)
        elab=-1
        ecms=-1
        pnll=-1
        if(jpsi.lt.0)then
  11      z=0.19*sqrt(-2*alog(rangen()))*cos(2*pi*rangen())
          engy=abs(z)*engmax
          if(engy.lt.egymin)goto11
        endif
      elseif(noebin.lt.0)then
        call iniFakeEA
      endif

      if(iappl.le.3)then  !------------------ < 3 ------------------

      call iniCR

      call idmass(idproj,amproj)
      call idmass(idtarg,amtarg)
      call idspin(idproj,ispin,jspin,istra)
      isoproj=sign(1,idproj)*ispin
      call idspin(idtarg,ispin,jspin,istra)
      isotarg=sign(1,idtarg)*ispin
      call idchrg( 1 ,idproj,chrg)
      ichproj=nint(chrg)
      call idchrg( 2 ,idtarg,chrg)
      ichtarg=nint(chrg)
      nre=0
      if(engy.ge.0.)nre=nre+1
      if(pnll.ge.0.)nre=nre+1
      if(elab.ge.0.)nre=nre+1
      if(ekin.ge.0.)nre=nre+1
      if(ecms.ge.0.)nre=nre+1
      if(nre.ne.1)stop'invalid energy definition'
      ifirstghe=0
 101  continue
      if(engy.gt.0.)then
        pnll=sqrt(amproj**2+amtarg**2)
        pnll=(engy-pnll)*(engy+pnll)*0.5/amtarg
        pnll=sqrt(max(0.,(pnll-amproj)*(pnll+amproj)))
        ! pnll=sqrt(max(0., ((engy**2-amproj**2-amtarg**2)/2/amtarg)**2
        !&                   -amproj**2) )
        elab=sqrt(pnll**2+amproj**2)
        ekin=elab-amproj
        ecms=engy
      elseif(ecms.gt.0.)then
        engy=ecms
        pnll=sqrt(amproj**2+amtarg**2)
        pnll=(engy-pnll)*(engy+pnll)*0.5/amtarg
        pnll=sqrt(max(0.,(pnll-amproj)*(pnll+amproj)))
        ! pnll=sqrt(max(0., ((engy**2-amproj**2-amtarg**2)/2/amtarg)**2
        !&                   -amproj**2) )
        elab=sqrt(pnll**2+amproj**2)
        ekin=elab-amproj
      elseif(elab.gt.0)then
        pnll=sqrt(max(0.,(elab-amproj)*(elab+amproj)))
        engy=sqrt( 2*elab*amtarg+amtarg**2+amproj**2 )
        ecms=engy
        ekin=elab-amproj
      elseif(pnll.gt.0)then
        elab=sqrt(pnll**2+amproj**2)
        engy=sqrt(2*sqrt(pnll**2+amproj**2)*amtarg+amtarg**2+amproj**2)
        ecms=engy
        ekin=elab-amproj
      elseif(ekin.gt.0.)then
        elab=ekin+amproj
        pnll=sqrt(max(0.,(elab-amproj)*(elab+amproj)))
        engy=sqrt( 2*elab*amtarg+amtarg**2+amproj**2 )
        ecms=engy
      endif
      if(model.eq.3.and.ifirstghe.eq.0)then    !det, trit and alp
        if(maproj.eq.2.and.laproj.eq.1)idproj=17
        if(maproj.eq.3.and.laproj.eq.1)idproj=18
        if(maproj.eq.4.and.laproj.eq.2)idproj=19
        if(idproj.ge.17.and.idproj.le.19)then
          elab=elab*maproj
          call idmass(idproj,amproj)
          maproj=1
          laproj=-1
          ifirstghe=1
          engy=-1
          ecms=-1
          pnll=-1
          ekin=-1
          goto 101
        endif
      endif

      if(pnll.le.0.001)call utstop('ainit: energy too low&')
      if(engy.gt.egymax)call utstop('ainit: energy too high&')
      s=engy**2
      pnullx=utpcm(engy,amproj,amtarg)
      yhaha=alog((sqrt(pnll**2+s)+pnll)/sqrt(s))
      ypjtl=alog((sqrt(pnll**2+amproj**2)+pnll)/amproj)
      if(noebin.lt.0)call iniFakeEA2

      elseif(iappl.eq.7)then !------------------ 7 ------------------

      call idmass(idproj,amproj)
      call idmass(idtarg,amtarg)
      if(elab.gt.0)then
        pnll=sqrt(max(0.,elab**2-amproj**2))
        engy=amproj
        ecms=engy
        ekin=elab-amproj
      elseif(pnll.gt.0)then
        elab=sqrt(pnll**2+amproj**2)
        engy=amproj
        ecms=engy
        ekin=elab-amproj
      elseif(ekin.gt.0.)then
        elab=ekin+amproj
        pnll=sqrt(max(0.,elab**2-amproj**2))
        engy=amproj
        ecms=engy
      else
        engy=amproj
        ecms=amproj
        elab=0.
        pnll=0.
        ekin=0.
      endif

      pnullx=0.
      ypjtl=alog((sqrt(pnll**2+amproj**2)+pnll)/amproj)
      yhaha=ypjtl

      elseif(iappl.eq.9)then !-----------------9------------------

        engy=tecm

      elseif(engy.gt.0.)then !------------------------------------

        ecms=engy

      endif                 !------------------------------------

      !call setParamEdep !moved to KW/tab.f
      call setParams()

      detap=(ypjtl-yhaha)*etafac
      detat=-yhaha*etafac
      tpro=dcosh(detap)
      zpro=dsinh(detap)
      ttar=dcosh(detat)
      ztar=dsinh(detat)

      egyevt=engy
      ekievt=ekin
      pmxevt=pnll

      if(iappl.gt.9)stop'update following statement'
      if(iappl.ge.5.and.iappl.le.8)then
      s=12.**2
      endif

      if(iappl.ge.5.and.iappl.le.8)then
      s=12.**2
      endif
      if(iappl.le.3)then
       if(maproj.gt.1)then
        rpj=1.19*maproj**(1./3.)-1.61*maproj**(-1./3.)
        rmproj=rpj+fctrmx*.54
        rcproj=dble(rpj/cosh(yhaha)*facnuc)
       else
        rmproj=0
        rcproj=dble(0.8/cosh(yhaha)*facnuc)
       endif
       if(matarg.gt.1)then
        rtg=1.19*matarg**(1./3.)-1.61*matarg**(-1./3.)
        rmtarg=rtg+fctrmx*.54
        rctarg=dble(rtg/cosh(yhaha)*facnuc)
       else
        rmtarg=0
        rctarg=dble(0.8/cosh(yhaha)*facnuc)
       endif

      endif

      call iclass(idproj,iclpro)
      call iclass(idtarg,icltar)
      if(icltar.eq.2)then
        iclreg=iclpro
        ichreg=ichproj
      elseif(iclpro.eq.2)then
        iclreg=icltar
        ichreg=ichtarg
      else
        call utstop("Both proj and targ .ne. nucleon not allowed !&")
      endif

      if(inicnt.eq.1)then
        call hdecin(.false.)
        if(iappl.eq.1.or.iappl.ge.5.and.iappl.le.8)then
          c=6
          !call utquaf(sptj,nptj,xptj,qptj,wptj,0.,.33*c,.66*c,c)
        endif
        if(model.eq.1)then
          if(iclegy2.gt.1)then
            egyfac=(egymax*1.0001/egylow)**(1./float(iclegy2-1))
          else
            egyfac=1.
          endif
          call conini
        else
          iorsce=0
          iorsdf=0
          iorshh=0
          iorsdf=0
        endif
        if(model.eq.1.or.model.eq.2)then
          call psaini
          call disini
        endif
      endif

      if(model.eq.1)then                   !only for epos
        koll=1      !because it's needed in Gfunpar
        if(iappl.le.3)then
          call emsini(engy,idproj,idtarg,0)      !cross-section (sig) set here
          if(ixtau.eq.1)call xtauev(0)
          if(iEmsB.eq.1)call xEmsB(0,0,0)
          if(iEmsBg.eq.1)call xEmsBg(0,0,0,0)
          if(iEmsPm.eq.1)call xEmsPm(0,0,0,0)
          if(iEmsPx.eq.1)call xEmsPx(0,0.,0.,0.,0.,0,0.,0.)
          if(iEmsPBx.eq.1)call xEmsP2(0,0,0,0.,0.,0.,0.,0.,0.,0.,0.)
          if(iEmsPDF.eq.1)call xEmsP3(0,0,0,0,0.,0)
          if(iEmsSe.eq.1)call xEmsSe(0,0.,0.,0,1)
          if(iEmsSe.eq.1)call xEmsSe(0,0.,0.,0,2)
          if(iEmsDr.eq.1)call xEmsDr(0,0.,0.,0)
          if(iEmsRx.eq.1)call xEmsRx(0,0,0.,0.)
        endif
      else
        call IniEvtModel
      endif



      if(idtarg.eq.0)idtarg=1120 !air = nucleus

      !call testXpXm

      !call testBW

      !call xhynRanBoseFermi

ccc      call MakeFpartonTable

c$$$      call testconex(1)

 9999 continue

      if(inicnt.eq.1.and.noebin.ge.0)then
        if(seedj2.ne.0d0)then
          call ranfcv(seedj2)
          write(ifmt,'(a)')
     &"Random number sequence does not start at 0 ... please wait !"
        endif
        if(seedj.ne.0d0)then
          call ranfini(seedj,iseqsim,2)
        else
          stop 'seedi = 0 ... Please define it !'
        endif
        call aseed(2)
      else                    !to use the proper random sequence
        call ranfini(seedc,iseqsim,0)
        if(noebin.ge.0.and.nevent.gt.0)call aseed(2)
      endif

      call utprix('ainit ',ish,ishini,4)
      call clop(3)
      return
      end

c---------------------------------------------------------------------
      subroutine wrxx
c---------------------------------------------------------------------
#include "aaa.h"
      character*80 twritexx
      common/cwritexx/nwritexx,twritexx(20)
      if(nwritexx.eq.0)return
      do n=1,nwritexx
      write(ifhi,'(a)')twritexx(n)
      enddo
      end

c---------------------------------------------------------------------
      subroutine wrxxx
c---------------------------------------------------------------------
#include "aaa.h"
      character*80 twritexxx
      common/cwritexxx/nwritexxx,twritexxx(50)
      if(nwritexxx.eq.0)return
      do n=1,nwritexxx
      write(ifhi,'(a)')twritexxx(n)
      enddo
      end

c---------------------------------------------------------------------
      subroutine aread
c---------------------------------------------------------------------
c  reads and interprets input commands
c---------------------------------------------------------------------

#include "aaa.h"
#include "ico.h"
#include "ho.h"
#include "par.h"
#include "ems.h"
#include "sem.h"
#include "so.h"
#include "sf.h"
#include "tab.h"
#include "xan.h"

      double precision histoweight
      common/chiswei/histoweight
      common/cyield/yield/cifset/ifset/caverg/averg
      common/csigma/sigma
      double precision key
      double precision val,val1,val2
      character*1000 line,linex,cline
      data nappl /0/
      common/record/maxrec(2),irecty(30,2)
      common/cfacmss/facmss /cr3pomi/r3pomi
      common /ems12/iodiba,bidiba  ! defaut iodiba=0. if iodiba=1, study H-Dibaryon
      character*1000 fnjob
      common/jobfname/  fnjob,nfnjob
      character*500 fndat,fnncs,fnIIdat,fnIIncs,fnII03dat,fnII03ncs,
     &fndpmjet,fndpmjetpho
      common/dpmjetfname/  fndpmjet,fndpmjetpho
      common/qgsfname/  fndat, fnncs, ifdat, ifncs
      common/qgsIIfname/fnIIdat, fnIIncs, ifIIdat, ifIIncs !qgs-II
      common/qgsII03fname/fnII03dat, fnII03ncs, ifII03dat, ifII03ncs !qgs-II03
      common/qgsnfname/ nfndat, nfnncs
      common/qgsIInfname/ nfnIIdat, nfnIIncs     !qgs-II
      common/qgsII03nfname/ nfnII03dat, nfnII03ncs     !qgs-II03
      common/ghecsquel/anquasiel,iquasiel
      common/cjjj/jjj,cline
      character*400 fnamein, cbasout
      common/croot2/fnamein,ix1 /croot6/cbasout,iba,nout
      character cmodel*21, fmt4*4
      common/cbincond/nozero,ibmin,ibmax
      common/photrans/phoele(4),ebeam,noevt /cjjtb/jjtb,jjeos
      common/cisk/iskmin,isk2min   /crapcol/rapcol
      common/ciuelast/iuelast /ciuskip/iuskip
      common/ciuchaskip/iuchaskip /ciunostring/iunostring
      character*80 twritexx,twritexxx
      character*4 cfmt
      common/cwritexx/nwritexx,twritexx(20) /cieof/ieof
      common/cwritexxx/nwritexxx,twritexxx(50)
      common/cicentrality/icentrality
      common/cihifcount/ihifcount(9),ifhix(9)
      integer ihifcounti(9)
      data itit/1/ ishxxx/1/  ihifcounti/9*0/  iskkey/2/
      common/cijetfluid/ijetfluid  /ciotype/iotype
      common /cnnnhis/nnnhis
      character cext1*10
      common/ccext1/cext1
      character cext3*10
      common/ccext3/iext3,cext3 /ciext4/iext4
      common/cgefac/gefac
      common/ciprotectinirj/iprotectinirj
      common/cigrpac/igrpac
      common/cnfifac/nfifac
      !-----should eventually move to aaa.h
      integer ioTestFact
      common /ctestfact/ ioTestFact
      integer laddTestFact
      common /claddtestfact/ laddTestFact
      integer noptTestFact
      common /cnopttestfact/ noptTestFact
      integer idhprTestFact
      common /cidhprtestfact/ idhprTestFact
      integer ihqTestFact
      common /chqtestfact/ ihqTestFact
      integer iffsigiEfficiency
      common /ciffsigiEfficiency/ iffsigiEfficiency
      character*30 hdtext(8)
      common/chdtext/hdtext
      common/cmodshox/modshox
      integer ioicoplot
      common /cioicoplot/ioicoplot
      real psum(4)
      !--------
      save itit,ishxxx,nskip,iskkey
      do i=1,9
      ihifcount(i)=ihifcounti(i)
      ifhix(i)=0
      enddo

      call setinp(fnjob,nfnjob,ifop)

      j=-1

      if(nopen.ne.-1)then       !only first read
        jcentrality=0
        nskip=0
        icentrality=0
      endif
      iaddplot=0
      do i=1,mxaddplot1
      do j=1,mxaddplot2
      caddplot(i,j)='                    '
      enddo
      enddo
      nhsto=0
      ndefine=0
      ncentrality=1
      nwritexx=0
      nwritexxx=0
      iskmin=1
      isk2min=1
      ifhiSave=0
      iprotectinirj=0
      jj=0
      ncontr=0
      bwidth=0.
      iPFE=0

      j=-1

    1 call utword(line,i,j,1)

          if(line(i:j).eq.'#define')then

      call setDefine(line,i,j,1)

          elseif(line(i:j).eq.'#define2')then

      call setDefine(line,i,j,2)

          elseif(line(i:j).eq.'not')then

       itit=0
       ishxxx=0
       iecho=0

          elseif(line(i:j).eq.'goto')then

      iechox=iecho
      iecho=0
      call utword(line,i,j,ne)
      ix=i
      jx=j
      linex=line
      call utword(line,i,j,ne)
      do while(line(i:j).ne.linex(ix:jx))
      call utword(line,i,j,ne)
      enddo
      iecho=iechox
      goto1

          elseif(line(i:j).eq.'#ifCentralityZero')then

      if(jcentrality.ne.0)then ! not min bias (centrality 0)
      iechox=iecho
      iecho=0
      call utword(line,i,j,ne)
      do while(line(i:j).ne.'#fiCentralityZero')
      call utword(line,i,j,ne)
      enddo
      iecho=iechox
      goto1
      endif

          elseif(line(i:j).eq.'#fiCentralityZero')then

      continue

          elseif(line(i:j).eq.'#ifNotReadingRoot')then

      if(igrTree.gt.0)then ! reading root
      iechox=iecho
      iecho=0
      call utword(line,i,j,ne)
      do while(line(i:j).ne.'#fi')
      call utword(line,i,j,ne)
      enddo
      iecho=iechox
      goto1
      endif

          elseif(line(i:j).eq.'#ifReadingRoot')then

      if(.not.(igrTree.gt.0))then ! not reading root
      iechox=iecho
      iecho=0
      call utword(line,i,j,ne)
      do while(line(i:j).ne.'#fi')
      call utword(line,i,j,ne)
      enddo
      iecho=iechox
      goto1
      endif

           elseif(line(i:j).eq.'application')then

      call utworn(line,j,ne)
      if(ne.eq.0.and.iprmpt.gt.0)write(ifmt,'(a)')'application?'
      call utword(line,i,j,0)
      if(nopen.ne.-1)then       !only first read
      if(line(i:j).eq.'conversion')iappl=-1
      if(line(i:j).eq.'analysis')  iappl=0
      if(line(i:j).eq.'hadron')    iappl=1
      if(line(i:j).eq.'geometry')  iappl=2
      if(line(i:j).eq.'read')      iappl=3
      if(line(i:j).eq.'thermal')     iappl=4
      if(line(i:j).eq.'kinky')     iappl=5
      if(line(i:j).eq.'ee')        iappl=6
      if(line(i:j).eq.'decay')     iappl=7
      if(line(i:j).eq.'lepton')    iappl=8
      if(line(i:j).eq.'micro')     iappl=9
      if(line(i:j).eq.'hydro')
     . stop'\n\n in aread; application hydro no longer supported \n\n'
      if(line(i:j).eq.'ee')    then
        naflav=5                ! number of flavors in ee
      endif
      nappl=nappl+1
      if(iappl.ne.0.and.nappl.gt.1)call aaset(1)
      if(iappl.eq.0.and.nappl.gt.1)call aaset(2)
      if(iappl.eq.0)jframe=iframe
      if(iappl.eq.0)kframe=iframe
      if(iappl.eq.1)iframe=0
      if(iappl.eq.2)iframe=0
      if(iappl.eq.4)iframe=1
      if(iappl.eq.5)iframe=1
      if(iappl.eq.6)iframe=1
      if(iappl.eq.7)iframe=0
      if(iappl.eq.8)iframe=21        !gncm
      if(iappl.eq.9)iframe=1
      endif

           elseif(line(i:j).eq.'call')then

      call utworn(line,j,ne)
      if(ne.eq.0.and.iprmpt.gt.0)write(ifmt,'(a)')'subroutine?'
      call utword(line,i,j,0)

      if(nopen.eq.-1)then       !-----only second run

      if(.not.(igrTree.gt.0))then !not reading root

      if(line(i:j).eq.'xjerr')call xjerr
      if(line(i:j).eq.'xParams')write(ifmt,*)'xParams obsolete'
      if(line(i:j).eq.'xIco3')call xIco3
      if(line(i:j).eq.'xIco')call xIco(-1)
c      if(line(i:j).eq.'xEnergy')call xEnergy
      if(line(i:j).eq.'xConThickProj')call xConThick(1)
      if(line(i:j).eq.'xConThickTarg')call xConThick(2)
      if(line(i:j).eq.'xConNuclDensProj')call xConNuclDens(1)
      if(line(i:j).eq.'xConNuclDensTarg')call xConNuclDens(2)
      if(line(i:j).eq.'xConNuclDensProjTarg')call xConNuclDens(1)
      if(line(i:j).eq.'xConNuclDensProjTarg')call xConNuclDens(2)
      if(line(i:j).eq.'xGeometry')call xGeometry(2)
      if(line(i:j).eq.'xbDens')call xbDens(2)
      if(line(i:j).eq.'xEpsilon')call xEpsilon(2)
      if(line(i:j).eq.'xZnucTheo')call xZnucTheo
      if(line(i:j).eq.'xRanPt')call xRanPt
      if(line(i:j).eq.'xParGam')call xParGam
      if(line(i:j).eq.'xParGampp')call xParGampp
      if(line(i:j).eq.'xParOmega1xy')call xParOmega1xy
c$$$      if(line(i:j).eq.'xParOmega3xyz')call xParOmega3xyz
      if(line(i:j).eq.'xParPro')call xParPro
      if(line(i:j).eq.'xParPro1')call xParPro1
      if(line(i:j).eq.'xParPomInc')call xParPomInc
      if(line(i:j).eq.'xParPomIncX')call xParPomIncX
      if(line(i:j).eq.'xParPomIncP')call xParPomIncP
      if(line(i:j).eq.'xParPomIncM')call xParPomIncM
      if(line(i:j).eq.'xParPomIncXI')call xParPomIncXI
      if(line(i:j).eq.'xParPomIncPI')call xParPomIncPI
      if(line(i:j).eq.'xParPomIncMI')call xParPomIncMI
      if(line(i:j).eq.'xParPomIncJ')call xParPomIncJ
      if(line(i:j).eq.'xParPomIncDiff')call xParPomIncDiff
      if(line(i:j).eq.'xParOmega1')call xParOmega1
c$$$      if(line(i:j).eq.'xParOmega3')call xParOmega3
c$$$      if(line(i:j).eq.'xParOmega5')call xParOmega5
      if(line(i:j).eq.'xParOmegaN')call xParOmegaN
      if(line(i:j).eq.'xParGauss')call xParGauss
      if(line(i:j).eq.'xParSigma')call xParSigma
c$$$      if(line(i:j).eq.'xParSigma2')call xParSigma2
c$$$      if(line(i:j).eq.'xScrD')call xScrD
      if(line(i:j).eq.'xFitD1')call xFitD1
c$$$      if(line(i:j).eq.'xExaD2')call xExaD2
      if(line(i:j).eq.'xbExaD')call xbExaD
c$$$      if(line(i:j).eq.'xbExaD2')call xbExaD2
      if(line(i:j).eq.'xbnExaD')call xbnExaD
c$$$      if(line(i:j).eq.'xbnExaD2')call xbnExaD2
      if(line(i:j).eq.'xFitD2')call xFitD2
      if(line(i:j).eq.'xFitD3')call xFitD3
      if(line(i:j).eq.'xbParD')call xbParD
c$$$      if(line(i:j).eq.'xParD2')call xParD2
      if(line(i:j).eq.'xGexaJ')call xGexaJ
      if(line(i:j).eq.'xGexaDiff')call xGexaDiff
      if(line(i:j).eq.'xyGexaDiff')call xyGexaDiff
      if(line(i:j).eq.'xbnParD')call xbnParD
      if(line(i:j).eq.'xsParD')call xsParD
      if(line(i:j).eq.'xfzero')call xfzero
c$$$      if(line(i:j).eq.'xmParD2')call xmParD2
      if(line(i:j).eq.'xyParD')call xyParD
c$$$      if(line(i:j).eq.'xyParD2')call xyParD2
      if(line(i:j).eq.'xParPhi1')call xParPhi1
      if(line(i:j).eq.'xParPhi')call xParPhi
      if(line(i:j).eq.'xParH')call xParH
      if(line(i:j).eq.'xParHPhiInt')call xParHPhiInt
      if(line(i:j).eq.'xParZ')call xParZ
      if(line(i:j).eq.'xAlphaS')call xAlphaS
      if(line(i:j).eq.'xAlphaSx')call xAlphaSx
      if(line(i:j).eq.'xtauev')call xtauev(2)
      if(line(i:j).eq.'xspace')call xspace(2)
      if(line(i:j).eq.'gakjto'   )call gakjto
      if(line(i:j).eq.'psaevp')call psaevp !stop'\n\n psaevp removed in epos218 \n\n'
c     if(line(i:j).eq.'pyarea')call pyarea
      if(line(i:j).eq.'xjden1')call xjden1(2,0,0.,0.,0.,0.,0.)
      if(line(i:j).eq.'xjden2')call xjden2(2,0,0.,0.,0.,0.)
c     if(line(i:j).eq.'xjdis' )call xjdis(2,0,0)

      if(model.eq.1)then

      if(line(i:j).eq.'xEmsB' )call xEmsB(2,0,0)
      if(line(i:j).eq.'xEmsBg')call xEmsBg(2,0,0,0)
      if(line(i:j).eq.'xEmsPm')call xEmsPm(2,0,0,0)
      if(line(i:j).eq.'xEmsPx')call xEmsPx(2,0.,0.,0.,0.,0,0.,0.)
      if(line(i:j).eq.'xEmsSe')call xEmsSe(2,0.,0.,0,1)
      if(line(i:j).eq.'xEmsSe')call xEmsSe(2,0.,0.,0,2)
      if(line(i:j).eq.'xEmsDr')call xEmsDr(2,0.,0.,0)
      if(line(i:j).eq.'xEmsRx')call xEmsRx(2,0,0.,0.)
      if(line(i:i+5).eq.'xEmsP2')then
        read(line(j-1:j-1),*)val
        idh=nint(val)
        read(line(j:j),*)val
        jex=nint(val)
        if(line(j-3:j-2).eq.'PE')
     &       call xEmsP2(2,idh,jex,0.,0.,0.,0.,0.,0.,0.,0.)
        if(line(j-4:j-2).eq.'PEx')
     &       call xEmsP2(5,idh,jex,0.,0.,0.,0.,0.,0.,0.,0.)
        if(line(j-3:j-2).eq.'IB')
     &       call xEmsP2(3,idh,jex,0.,0.,0.,0.,0.,0.,0.,0.)
        if(line(j-3:j-2).eq.'OB')
     &       call xEmsP2(4,idh,jex,0.,0.,0.,0.,0.,0.,0.,0.)
        if(line(j-4:j-2).eq.'GBi')
     &       call xEmsP2(6,idh,jex,0.,0.,0.,0.,0.,0.,0.,0.)
        if(line(j-4:j-2).eq.'GBa')
     &       call xEmsP2(7,idh,jex,0.,0.,0.,0.,0.,0.,0.,0.)
        if(line(j-4:j-2).eq.'PEB')
     &       call xEmsP2(8,idh,jex,0.,0.,0.,0.,0.,0.,0.,0.)
        if(line(j-4:j-2).eq.'IBP')
     &       call xEmsP2(9,idh,jex,0.,0.,0.,0.,0.,0.,0.,0.)
      endif
      if(line(i:i+6).eq.'xEmsPDF')then
        read(line(j-2:j-2),*)val
        idh=nint(val)
        read(line(j-1:j-1),*)val
        jex=nint(val)
        read(line(j:j),*)val
        idi=nint(val)
        if(line(j-3:j-3).eq.'p')
     &       call xEmsP3(2,idh,jex,idi,-1.,1)
        if(line(j-3:j-3).eq.'t')
     &       call xEmsP3(2,idh,jex,idi,-1.,2)
      endif
      if(line(i:j).eq.'xConxyzProj')
     &stop'xConxyzProj->xConNuclDensProj'
      if(line(i:j).eq.'xConxyzTarg')
     &stop'xConxyzTarg->xConNuclDensTarg'
      if(line(i:j).eq.'xConxyzProjTarg')
     &stop'xConxyzProjTarg->xConNuclDensProjTarg'

      endif  ! model 1

      else   ! reading root

      if(line(i:i+5).eq.'xEmsP2')then
      do ix=j+1,j+20
      if(line(ix:ix+4).eq.'write')then
      do ixx=ix+5,ix+20
      if(line(ixx:ixx+5).eq.'plot 0')then
      j=ixx+6
      endif
      if(line(ixx:ixx+6).eq.'plot 0-')then
      j=ixx+7
      endif
      enddo
      endif
      enddo
      endif

      endif  ! reading root

      elseif(model.eq.1)then  !first run and model 1
      if(.not.(igrTree.gt.0))then !not reading root

      if(line(i:j).eq.'xGeometry')then
       call xGeometry(0)
       ixgeometry=1
      elseif(line(i:j).eq.'xEpsilon')then
       call xEpsilon(0)
      elseif(line(i:j).eq.'xbDens')then
       ixbDens=1
      elseif(line(i:j).eq.'xtauev')then
       ixtau=1
      elseif(line(i:j).eq.'xEmsB')then
       iEmsB=1
      elseif(line(i:j).eq.'xEmsBg')then
       iEmsBg=1
      elseif(line(i:j).eq.'xEmsPm')then
       iEmsPm=1
      elseif(line(i:j).eq.'xEmsPx')then
       iEmsPx=1
      elseif(line(i:i+5).eq.'xEmsP2')then
       iEmsPBx=1
      elseif(line(i:i+6).eq.'xEmsPDF')then
       iEmsPDF=1
      elseif(line(i:j).eq.'xEmsSe')then
       iEmsSe=1
      elseif(line(i:j).eq.'xEmsDr')then
       iEmsDr=1
      elseif(line(i:j).eq.'xEmsRx')then
       iEmsRx=1
      elseif(line(i:j).eq.'xEmsI1')then
       iEmsI1=1
       if(iEmsI1+iEmsI2.eq.1)write(ifhi,'(a)')'newpage zone 3 4 1'
      elseif(line(i:j).eq.'xEmsI2')then
       iEmsI2=1
       if(iEmsI1+iEmsI2.eq.1)write(ifhi,'(a)')'newpage zone 3 4 1'
      elseif(line(i:j).eq.'xSpaceTime')   then
       iSpaceTime=1
      elseif(line(i:j).eq.'xSpaceTime(1)')then
       iSpaceTime=1
      elseif(line(i:j).eq.'xSpaceTime(2)')then
       iSpaceTime=2
      elseif(line(i:j).eq.'xSpaceTime(3)')then
       iSpaceTime=3
      elseif(line(i:j).eq.'xSpaceTime(4)')then
       iSpaceTime=4
      elseif(line(i:j).eq.'xSpaceTime(5)')then
       iSpaceTime=5
      elseif(line(i:j).eq.'xSpaceTime(6)')then
       iSpaceTime=6
      elseif(line(i:j).eq.'xSpaceTime(7)')then
       iSpaceTime=7
      elseif(line(i:j).eq.'xSpaceTime(8)')then
       iSpaceTime=8
      elseif(line(i:j).eq.'xSpaceTime(31)')then
       iSpaceTime=31
      elseif(line(i:j).eq.'xSpaceTime(32)')then
       iSpaceTime=32
      elseif(line(i:j).eq.'xSpaceTime(33)')then
       iSpaceTime=33
      elseif(line(i:j).eq.'xxSpaceTime')then
       stop'xxSpaceTime->xSpaceTime.'
      endif

      endif !not reading root
      endif !first run and model 1

           elseif(line(i:j).eq.'decayall')then

      nrnody=0

           elseif(line(i:j).eq.'echo')then

      call utworn(line,j,ne)
      if(ne.eq.0.and.iprmpt.gt.0)write(ifmt,'(a)')'on or off?'
      call utword(line,i,j,0)
      if(line(i:j).eq.'on')iecho=1
      if(line(i:j).eq.'off')iecho=0
      if(line(i:j).ne.'on'.and.line(i:j).ne.'off')stop'invalid option'
      if(nopen.eq.-1)iecho=0
      if(ishxxx.eq.0)iecho=0

           elseif(line(i:j).eq.'|')then

      iecho=0

           elseif(line(i:j).eq.'fdpmjet')then              !DPMJET

      call utworn(line,j,ne)
      if(ne.eq.0.and.iprmpt.gt.0)write(ifmt,'(a)')'file-type file-name?'
      call utword(line,i,j,0)
      linex=line
      ix=i
      jx=j
      call utworn(line,j,ne)
      if(ne.eq.0.and.iprmpt.gt.0)write(ifmt,'(a)')'file-name?'
      call utword(line,i,j,0)
      if(linex(ix:jx).eq.'dat')fndpmjet(1:j-i+1)=line(i:j)
      if(linex(ix:jx).eq.'dat')nfndat=j-i+1             !length of dpmjet.dat path
      if(nfndat.gt.1)ifdat=1

          elseif(line(i:j).eq.'fdpmjetpho')then              !DPMJET phojet fitpar file

      call utworn(line,j,ne)
      if(ne.eq.0.and.iprmpt.gt.0)write(ifmt,'(a)')'file-type file-name?'
      call utword(line,i,j,0)
      linex=line
      ix=i
      jx=j
      call utworn(line,j,ne)
      if(ne.eq.0.and.iprmpt.gt.0)write(ifmt,'(a)')'file-name?'
      call utword(line,i,j,0)
      if(linex(ix:jx).eq.'dat')fndpmjetpho(1:j-i+1)=line(i:j)
      if(linex(ix:jx).eq.'dat')nfndat=j-i+1             !length of dpmjet.dat path
      if(nfndat.gt.1)ifdat=1

           elseif(line(i:j).eq.'fqgsjet')then              !QGSJet

      call utworn(line,j,ne)
      if(ne.eq.0.and.iprmpt.gt.0)write(ifmt,'(a)')'file-type file-name?'
      call utword(line,i,j,0)
      linex=line
      ix=i
      jx=j
      call utworn(line,j,ne)
      if(ne.eq.0.and.iprmpt.gt.0)write(ifmt,'(a)')'file-name?'
      call utword(line,i,j,0)
      if(linex(ix:jx).eq.'dat')fndat(1:j-i+1)=line(i:j)
      if(linex(ix:jx).eq.'dat')nfndat=j-i+1             !length of qgsdat01-file name
      if(linex(ix:jx).eq.'ncs')fnncs(1:j-i+1)=line(i:j)
      if(linex(ix:jx).eq.'ncs')nfnncs=j-i+1             !length of sectnu-file name
      if(nfndat.gt.1)ifdat=1
      if(nfnncs.gt.1)ifncs=2

           elseif(line(i:j).eq.'fqgsjetII03')then              !QGSJET-II-03

      call utworn(line,j,ne)
      if(ne.eq.0.and.iprmpt.gt.0)write(ifmt,'(a)')'file-type file-name?'
      call utword(line,i,j,0)
      linex=line
      ix=i
      jx=j
      call utworn(line,j,ne)
      if(ne.eq.0.and.iprmpt.gt.0)write(ifmt,'(a)')'file-name?'
      call utword(line,i,j,0)
      if(linex(ix:jx).eq.'dat')fnII03dat(1:j-i+1)=line(i:j)
      if(linex(ix:jx).eq.'dat')nfnII03dat=j-i+1             !length of qgsjet-II.dat name
      if(linex(ix:jx).eq.'ncs')fnII03ncs(1:j-i+1)=line(i:j)
      if(linex(ix:jx).eq.'ncs')nfnII03ncs=j-i+1             !length of qgsjet-II.ncs name
      if(nfnII03dat.gt.1)ifII03dat=1
      if(nfnII03ncs.gt.1)ifII03ncs=2

           elseif(line(i:j).eq.'fqgsjetII')then              !QGSJET-II-04

      call utworn(line,j,ne)
      if(ne.eq.0.and.iprmpt.gt.0)write(ifmt,'(a)')'file-type file-name?'
      call utword(line,i,j,0)
      linex=line
      ix=i
      jx=j
      call utworn(line,j,ne)
      if(ne.eq.0.and.iprmpt.gt.0)write(ifmt,'(a)')'file-name?'
      call utword(line,i,j,0)
      if(linex(ix:jx).eq.'dat')fnIIdat(1:j-i+1)=line(i:j)
      if(linex(ix:jx).eq.'dat')nfnIIdat=j-i+1             !length of qgsjet-II.dat name
      if(linex(ix:jx).eq.'ncs')fnIIncs(1:j-i+1)=line(i:j)
      if(linex(ix:jx).eq.'ncs')nfnIIncs=j-i+1             !length of qgsjet-II.ncs name
      if(nfnIIdat.gt.1)ifIIdat=1
      if(nfnIIncs.gt.1)ifIIncs=2

           elseif(line(i:j).eq.'beginoptns')then

      call utword(line,i,j,0)
      if(line(i:j).eq.'spherio')call spherioptns(line,j)
      if(line(i:j).eq.'hq')call hqoptns(line,j)
      if(line(i:j).eq.'epos')continue

           elseif(line(i:j).eq.'endoptns')then

      call utword(line,i,j,0)
      linex=line
      ix=i
      jx=j
      if(linex(ix:jx).ne.'epos')stop 'beginoptns mismatch!'

           elseif(line(i:j).eq.'xfname'
     .        .or.line(i:j).eq.'RelativePathFileName'
     .        .or.line(i:j).eq.'ReFileName')then

      call utword(line,i,j,0)
      linex=line
      ix=i
      jx=j
      call utword(line,i,j,0)
      if(linex(ix:jx).eq.'ico'.or.linex(ix:jx).eq.'IcoTable')then
        nfnio=nfnnx+j-i+1
        fnio(1:nfnio)=fnnx(1:nfnnx)//line(i:j)
      elseif(linex(ix:jx).eq.'out'.or.linex(ix:jx).eq.'HydroTable')then
        nfnho(ireadhyt)=nfnnx+j-i+1
        fnho(ireadhyt)(1:nfnho(ireadhyt))=fnnx(1:nfnnx)//line(i:j)
      else
        stop'\n\n ERROR 29072011\n\n'
      endif

           elseif(line(i:j).eq.'writeeof')then

       ieof=1

           elseif(line(i:j).eq.'ext1')then

      cext1='          '
      call utword(line,i,j,0)
      if(line(i:j).ne.'-')then
        cext1=line(i:j)
      endif

           elseif(line(i:j).eq.'ext3')then

      cext3='          '
      call utword(line,i,j,0)
      iext3=j-i+1
      cext3(1:iext3)=line(i:j)

           elseif(line(i:j).eq.'ext4')then

      iext4=0
      call utword(line,i,j,0)
      if(line(i:j).eq.'0')iext4=1

           elseif(line(i:j).eq.'#if1')then

       call  setIf1(line,i,j)

           elseif(line(i:j).eq.'#else1')then

       call  setElse1(line,i,j)

           elseif(line(i:j).eq.'#if3')then

      call utword(line,i,j,0)
      i3skip=1
      n3=1
      if(line(i:i).eq.'#')then
        read(line(i+1:j),*)n3
        call utword(line,i,j,0)
      endif
      if(line(i:j).eq.cext3(1:iext3))i3skip=0
      do nuw=2,n3
        call utword(line,i,j,0)
        if(line(i:j).eq.cext3(1:iext3))i3skip=0
      enddo
      if(i3skip.eq.1)then
        do while(line(i:j).ne.'#fi')
          call utword(line,i,j,0)
        enddo
      else
        continue
      endif

           elseif(line(i:j).eq.'#if4')then

      call utword(line,i,j,0)
      i3skip=1
      n3=1
      if(line(i:i).eq.'#')then
        read(line(i+1:j),*)n3
        call utword(line,i,j,0)
      endif
      if(line(i:j).eq.cext3(1:iext3))i3skip=0
      do nuw=2,n3
        call utword(line,i,j,0)
        if(line(i:j).eq.cext3(1:iext3))i3skip=0
      enddo
      if(i3skip.eq.1)then
        do while(line(i:j).ne.'#fi4')
          call utword(line,i,j,0)
        enddo
      else
        continue
      endif

           elseif(line(i:j).eq.'#fi')then

      continue

           elseif(line(i:j).eq.'#fi1')then

      continue

           elseif(line(i:j).eq.'#fi4')then

      continue

           elseif(line(i:j).eq.'system')then

      call setSystem(line,i,j,isyst)

           elseif(line(i:j).eq.'rootcproot')then

      call setRootcproot(line,i,j,irootcproot,iboein)

           elseif(line(i:j).eq.'fname')then

      call utword(line,i,j,0)
      linex=line
      ix=i
      jx=j
      call utword(line,i,j,0)

      if(linex(ix:jx).eq.'pathnx')fnnx(1:j-i+1)=line(i:j)
      if(linex(ix:jx).eq.'pathep')fnnx(1:j-i+1)=line(i:j)
      if(linex(ix:jx).eq.'pathpdf')fngrv(1:j-i+1)=line(i:j)
      if(linex(ix:jx).eq.'check')fnch(1:j-i+1)=line(i:j)
      if(linex(ix:jx).eq.'mtr')fnmt(1:j-i+1)=line(i:j)
      if(linex(ix:jx).eq.'histo')fnhi(1:j-i+1)=line(i:j)
      if(linex(ix:jx).eq.'data') fndt(1:j-i+1)=line(i:j)
      if(linex(ix:jx).eq.'hepfile') fnhm(1:j-i+1)=line(i:j)
      if(linex(ix:jx).eq.'input')fnin(1:j-i+1)=line(i:j)
      if(linex(ix:jx).eq.'user1')fnus1(1:j-i+1+1)=line(i:j)//' '
      if(linex(ix:jx).eq.'user2')fnus2(1:j-i+1+1)=line(i:j)//' '
      if(linex(ix:jx).eq.'user3')fnus3(1:j-i+1+1)=line(i:j)//' '
      if(linex(ix:jx).eq.'user4')fnus4(1:j-i+1+1)=line(i:j)//' '
      if(linex(ix:jx).eq.'copy') fncp(1:j-i+1)=line(i:j)
      if(linex(ix:jx).eq.'initl') fnii(1:j-i+1)=line(i:j)
      if(linex(ix:jx).eq.'initl+') fnii(nfnii+1:nfnii+j-i+1)=line(i:j)
      if(linex(ix:jx).eq.'inidi') fnid(1:j-i+1)=line(i:j)
      if(linex(ix:jx).eq.'inidi+') fnid(nfnid+1:nfnid+j-i+1)=line(i:j)
      if(linex(ix:jx).eq.'inidr') fndr(1:j-i+1)=line(i:j)
      if(linex(ix:jx).eq.'inidr+') fndr(nfndr+1:nfndr+j-i+1)=line(i:j)
      if(linex(ix:jx).eq.'iniev') fnie(1:j-i+1)=line(i:j)
      if(linex(ix:jx).eq.'iniev+') fnie(nfnie+1:nfnie+j-i+1)=line(i:j)
      if(linex(ix:jx).eq.'inirj') fnrj(1:j-i+1)=line(i:j)
      if(linex(ix:jx).eq.'inirj+') fnrj(nfnrj+1:nfnrj+j-i+1)=line(i:j)
      if(linex(ix:jx).eq.'inics') fncs(1:j-i+1)=line(i:j)
      if(linex(ix:jx).eq.'inics+') fncs(nfncs+1:nfncs+j-i+1)=line(i:j)
      if(linex(ix:jx).eq.'partab')fn3p(1:j-i+1)=line(i:j)
      if(linex(ix:jx).eq.'dectab')fn3d(1:j-i+1)=line(i:j)
      if(linex(ix:jx).eq.'hpf') fnhpf(1:j-i+1)=line(i:j)
      if(linex(ix:jx).eq.'pathnx')nfnnx=j-i+1
      if(linex(ix:jx).eq.'pathep')nfnnx=j-i+1
      if(linex(ix:jx).eq.'pathpdf')nfngrv=j-i+1
      if(linex(ix:jx).eq.'check')nfnch=j-i+1
      if(linex(ix:jx).eq.'mtr')nfnmt=j-i+1
      if(linex(ix:jx).eq.'histo')nfnhi=j-i+1
      if(linex(ix:jx).eq.'data') nfndt=j-i+1
      if(linex(ix:jx).eq.'hepfile') nfnhm=j-i+1
      if(linex(ix:jx).eq.'input')nfnin=j-i+1
      if(linex(ix:jx).eq.'user1')nfnus1=j-i+1
      if(linex(ix:jx).eq.'user2')nfnus2=j-i+1
      if(linex(ix:jx).eq.'user3')nfnus3=j-i+1
      if(linex(ix:jx).eq.'user4')nfnus4=j-i+1
      if(linex(ix:jx).eq.'copy') nfncp=j-i+1
      if(linex(ix:jx).eq.'initl') nfnii=j-i+1
      if(linex(ix:jx).eq.'initl+')nfnii=nfnii+j-i+1
      if(linex(ix:jx).eq.'inidi') nfnid=j-i+1
      if(linex(ix:jx).eq.'inidi+')nfnid=nfnid+j-i+1
      if(linex(ix:jx).eq.'inidr') nfndr=j-i+1
      if(linex(ix:jx).eq.'inidr+')nfndr=nfndr+j-i+1
      if(linex(ix:jx).eq.'iniev') nfnie=j-i+1
      if(linex(ix:jx).eq.'iniev+')nfnie=nfnie+j-i+1
      if(linex(ix:jx).eq.'inirj') nfnrj=j-i+1
      if(linex(ix:jx).eq.'inirj+')nfnrj=nfnrj+j-i+1
      if(linex(ix:jx).eq.'inics') nfncs=j-i+1
      if(linex(ix:jx).eq.'inics+')nfncs=nfncs+j-i+1
      if(linex(ix:jx).eq.'partab')nfn3p=j-i+1
      if(linex(ix:jx).eq.'dectab')nfn3d=j-i+1
      if(linex(ix:jx).eq.'hpf')nfnhpf=j-i+1
      if(linex(ix:jx).eq.'ico'.or.linex(ix:jx).eq.'IcoTable')then
      nfnio=j-i+1
      fnio(1:nfnio)=line(i:j)
      elseif(linex(ix:jx).eq.'out'.or.linex(ix:jx).eq.'HydroTable')then
      nfnho(ireadhyt)=j-i+1
      fnho(ireadhyt)(1:nfnho(ireadhyt))=line(i:j)
      endif
      if(linex(ix:jx).eq.'check'.and.fnch(1:nfnch).ne.'none') then
        open(unit=ifcx,file=fnch(1:nfnch),status='unknown')
        kchopen=1
      elseif((linex(ix:jx).eq.'pathnx'.or.linex(ix:jx).eq.'pathep')
     .  .and.fnnx(1:nfnnx).ne.'none')then
        if(knxopen.eq.0)then
          knxopen=1
          call readidtable !should be called early
        endif
      elseif(linex(ix:jx).eq.'histo'.and.fnhi(1:nfnhi).ne.'none')then
        open(unit=ifhi,file=fnhi(1:nfnhi),status='unknown')
        khiopen=1
      elseif(linex(ix:jx).eq.'mtr'.and.fnmt(1:nfnmt).ne.'none')then
        ifmt=7
        open(unit=ifmt,file=fnmt(1:nfnmt),status='unknown')
        call atitle(0)
      elseif(linex(ix:jx).eq.'data'.and.fndt(1:nfndt).ne.'none')then
        open(unit=ifdt,file=fndt(1:nfndt),status='unknown')
        kdtopen=1
      elseif(linex(ix:jx).eq.'hepfile'.and.fnhm(1:nfnhm).ne.'none')then
        open(unit=ifhm,file=fnhm(1:nfnhm),status='unknown')
        khepmcopen=1  
      elseif(linex(ix:jx).eq.'copy'.and.fncp(1:nfncp).ne.'none')then
        open(unit=ifcp,file=fncp(1:nfncp),status='unknown')
        kcpopen=1
      elseif(linex(ix:jx).eq.'user1'.and.fnus1(1:nfnus1).ne.'none')then
        if(fnus1(nfnus1:nfnus1).ne.'/' )then
          open(unit=ifnus1,file=fnus1(1:nfnus1),status='unknown')
          write(ifmt,'(2a)')'open file ',fnus1(1:nfnus1)
        else
          write(ifmt,'(2a)')'def path fnus1 ',fnus1(1:nfnus1)
        endif
      elseif(linex(ix:jx).eq.'user2'.and.fnus2(1:nfnus2).ne.'none')then
        if(fnus2(nfnus2:nfnus2).ne.'/' )then
          open(unit=ifnus2,file=fnus2(1:nfnus2),status='unknown')
          write(ifmt,'(2a)')'open file ',fnus2(1:nfnus2)
        else
          write(ifmt,'(2a)')'def path fnus2 ',fnus2(1:nfnus2)
        endif
      elseif(linex(ix:jx).eq.'user3'.and.fnus3(1:nfnus3).ne.'none')then
        if(fnus3(nfnus3:nfnus3).ne.'/' )then
          open(unit=ifnus3,file=fnus3(1:nfnus3),status='unknown')
          write(ifmt,'(2a)')'open file ',fnus3(1:nfnus3)
        else
          write(ifmt,'(2a)')'def path fnus3 ',fnus3(1:nfnus3)
        endif
      elseif(linex(ix:jx).eq.'user4'.and.fnus4(1:nfnus4).ne.'none')then
        if(fnus4(nfnus4:nfnus4).ne.'/' )then
          open(unit=ifnus4,file=fnus4(1:nfnus4),status='unknown')
          write(ifmt,'(2a)')'open file ',fnus4(1:nfnus4)
        else
          write(ifmt,'(2a)')'def path fnus4 ',fnus4(1:nfnus4)
        endif
      endif

           elseif(line(i:j).eq.'frame')then

      call utworn(line,j,ne)
      if(ne.eq.0.and.iprmpt.gt.0)write(ifmt,'(a)')'frame?'
      call utword(line,i,j,0)
        if(nopen.ne.-1)then       ! event definition only, not analysis
        if(line(i:j).eq.'nucleon-nucleon')then
      if(iappl.eq.0)jframe=11
      if(abs(iappl).eq.1)iframe=11
      if(iappl.eq.3)iframe=11
      if(iappl.gt.3.or.iappl.eq.2)stop'invalid option nucleon-nucleon'
        elseif(line(i:j).eq.'target')then
      if(iappl.eq.0)jframe=12
      if(abs(iappl).eq.1)iframe=12
      if(iappl.eq.3)iframe=12
      if(iappl.gt.3.or.iappl.eq.2)stop'invalid option target'
        elseif(line(i:j).eq.'gamma-nucleon')then
      if(iappl.eq.0)then
        jframe=21
      elseif(iappl.le.3.and.iappl.ne.2)then
        iframe=21
      endif
      if((iappl+1)/2.eq.4)iframe=21
      if(iappl.ne.0.and.(iappl+1)/2.ne.4.and.iappl.ne.3
     *   .and.abs(iappl).ne.1)stop'invalid option gamma-nucleon'
        elseif(line(i:j).eq.'lab')then
      if(iappl.eq.0)then
        jframe=22
      elseif(iappl.le.3.and.iappl.ne.2)then
        iframe=22
      endif
      if((iappl+1)/2.eq.4)iframe=22
      if(iappl.ne.0.and.(iappl+1)/2.ne.4.and.iappl.ne.3
     *             .and.abs(iappl).ne.1)stop'invalid option lab'
        elseif(line(i:j).eq.'sphericity')then
      if(iappl.eq.0)jframe=32
      if(iappl.ne.0)stop'invalid option sphericity'
        elseif(line(i:j).eq.'thrust')then
      if(iappl.eq.0)jframe=33
      if(iappl.ne.0)stop'invalid option thrust'
        elseif(line(i:j).eq.'breit')then
          if(iappl.ne.0)stop'invalid option breit'
        else
      stop'frame not recognized'
        endif
        endif

           elseif(line(i:j).eq.'frame+')then

      call utword(line,i,j,0)

           elseif(line(i:j).eq.'epresolution')then

      call utword(line,i,j,0)
      read(line(i:j),*)kk
      do k=1,20 
        call utword(line,i,j,0)
        read(line(i:j),*)epreso(kk,k)
      enddo 

           elseif(line(i:j).eq.'binning')then

      call utworn(line,j,ne)
      if(ne.eq.0.and.iprmpt.gt.0)write(ifmt,'(a)')'log/lin?'
      call utword(line,i,j,0)
      if(line(i:j).eq.'lin')iologb=0
      if(line(i:j).eq.'log')iologb=1

           elseif(line(i:j).eq.'beginhisto'  !  | bh | beginanalysis
     .           .or.line(i:j).eq.'bh'
     .           .or.line(i:j).eq.'beginanalysis')then

      if(nopen.ne.-1)then !first read
        jjj=j
        cline=line
        call setCounters
        call xini
        j=jjj
        line=cline
      endif
      averg=-1e30

           elseif(line(i:j).eq.'endhisto' ! | eh | endanalysis
     .            .or.line(i:j).eq.'endanalysis'
     .            .or.line(i:j).eq.'eh')then

      if(nopen.eq.-1)then !second read
        nhsto=nhsto+1
        call xhis(nhsto)
      endif

           elseif(line(i:j).eq.'noweak')then

      if(nopen.ne.-1)stop'ERROR 08052021' !should not pass in first read 
      !only used in xini in first read

           elseif(line(i:j).eq.'weak')then

      if(nopen.ne.-1)stop'ERROR 08052021' !should not pass in first read 
      !only used in xini in first read

           elseif(line(i:j).eq.'noweak2')then

      if(nopen.ne.-1)stop'ERROR 08052021' !should not pass in first read 
      !only used in xini in first read

           elseif(line(i:j).eq.'noweak3')then

      if(nopen.ne.-1)stop'ERROR 08052021' !should not pass in first read 
      !only used in xini in first read

           elseif(line(i:j).eq.'histogram'
     .           .or.line(i:j).eq.'hi')then !-----------

      call utword(line,i,j,0)
      call utword(line,i,j,0)
      call utword(line,i,j,0)
      call utword(line,i,j,0)
      call utword(line,i,j,0)
      call utword(line,i,j,0)

           elseif(line(i:j).eq.'plot')then

      call utword(line,i,j,0)
      if(line(i:j).eq.'HydroEpsilon')       iHyEpsilon=1
      if(line(i:j).eq.'HydroEpsilon2')      iHyEpsilon2=1
      if(line(i:j).eq.'HydroEntropy')       iHyEntropy=1
      if(line(i:j).eq.'HydroTemperature')   iHyTemperature=1
      if(line(i:j).eq.'HydroRadVelocity')   iHyRadVelocity=1
      if(line(i:j).eq.'HydroLongVelocity')  iHyLongVelocity=1
      if(line(i:j).eq.'HydroAverages')      iHyAverages=1
      if(line(i:j).eq.'HydroBaryon')        iHyBaryon=1
      if(line(i:j).eq.'HydroFoVol')         iHyFoVol=1
      if(line(i:j).eq.'HydroFoRadius')      iHyFoRadius=1
      if(line(i:j).eq.'HydroFoEpsilon')     iHyFoEpsilon=1
      if(line(i:j).eq.'HydroFoRadVelocity') iHyFoRadVelocity=1
      if(line(i:j).eq.'HydroFoTangVelocity')iHyFoTangVelocity=1
      if(line(i:j).eq.'HydroFoBarmu')       iHyFoBarmu=1
      if(line(i:j).eq.'HydroEpsilonEta')    iHyEpsilonEta=1
      if(line(i:j).eq.'HydroEpsilonEtaLin') iHyEpsilonEta=2
      if(line(i:j).eq.'HydroBaryonEta')     iHyBaryonEta=1
      if(line(i:j).eq.'HydroBaryonEtaLin')  iHyBaryonEta=2
      if(line(i:j).eq.'HydroBarmuEta')      iHyBarmuEta=1
      if(line(i:j).eq.'HydroEntropyEta')    iHyEntropyEta=1
      if(line(i:j).eq.'HydroEos')           iHyEos=1
      if(line(i:j-3).eq.'HydroEmpty')read(line(j-1:j-1),*)iHyEmpty
      if(line(i:j-4).eq.'HydroEmpty')read(line(j-2:j-1),*)iHyEmpty

           elseif(line(i:i+5).eq.'root->')then

      if(nopen.eq.-1)then !second read
      if(line(i+6:i+11).eq.'IniCon')call xIniCon(line(i+6:j)//';')
      endif
      do k=1,30
        if(line(i+6+k:i+6+k+4).eq.'(x,y)')
     .  stop'ERROR: Use XY instead of (x,y) in optns file'
      enddo
      if(line(i+6:j).eq.'HydroEpsilonXY')iHoEpsilon=1
      if(line(i+6:j).eq.'HydroEpsilonEtas8XY')iHoEpsilonEtas8=1
      if(line(i+6:j).eq.'HydroEpsilonEtas6XY')iHoEpsilonEtas6=1
      if(line(i+6:j).eq.'HydroTemperatureEtas8XY')
     .                                       iHoTemperatureEtas8=1
      if(line(i+6:j).eq.'HydroTemperatureEtas6XY')
     .                                       iHoTemperatureEtas6=1
      if(line(i+6:j).eq.'HydroEpsilonTau8XY')iHoEpsilonTau8=1
      if(line(i+6:j).eq.'HydroTangVelTau6XY')iHoTangVelTau6=1

      if(line(i+6:j).eq.'HydroRadVelXY')iHoRadVel=1
      if(line(i+6:j).eq.'HydroTanVelXY')iHoTanVel=1
      if(line(i+6:i+20).eq.'HydroRadVelTau6')then
        if(line(i+21:i+22).eq.'XY')then
         iHoRadVelTau6=1
        elseif(line(i+22:i+23).eq.'XY')then
          read(line(i+21:i+21),'(i1)')iHoRadVelTau6
        elseif(line(i+23:i+24).eq.'XY')then
          read(line(i+21:i+22),'(i2)')iHoRadVelTau6
        endif
      endif
      if(line(i+6:i+21).eq.'HydroEpsilonTau6')then
        if(line(i+22:i+23).eq.'XY')then
          iHoEpsilonTau6=1
        elseif(line(i+23:i+24).eq.'XY')then
          read(line(i+22:i+22),'(i1)')iHoEpsilonTau6
        elseif(line(i+24:i+25).eq.'XY')then
          read(line(i+22:i+23),'(i2)')iHoEpsilonTau6
        endif
      endif
      if(line(i+6:j).eq.'HydroEpsilonEtaY')  iHoEpsilonEtaY=1
      if(line(i+6:i+16).eq.'MixedEvents')then
        read(line(i+18:j-1),*)val
        mixevt=nint(val)
      endif

          elseif(line(i:i+8).eq.'fillTree(')then

      write(ifmt,'(//70a/a/a/a/70a/)')('-',k=1,70)
     .,'Important info: New definitions of ior and jor, referring  '
     .,'                 to indices 0,1,2,... rather than 1,2,3,...' 
     .,'Use now fillTree4 rather than fillTree in optns'
     .,('-',k=1,70)
      stop'fillTree not known'

          elseif(line(i:i+8).eq.'fillTree4')then

       call setFillTree(line,i,j,irootcproot,ifillTree,izmode)

          elseif(line(i:i+3).eq.'fillH2')then

      ifillH2=1

          elseif(line(i:i+9).eq.'execSource')then

      if(line(i+11:i+11).ne.'B')stop'execSource error\n\n'
      read(line(i+12:i+12),*)iSource

          elseif(line(i+4:i+19).eq.'defineCentrality')then

      if(nopen.ne.-1)then       !only first read
      if(line(i+1:i+1).eq.'J')then
        jj=1
      elseif(line(i+1:i+1).eq.'K')then
        jj=2
      elseif(line(i+1:i+1).eq.'M')then
        jj=3
      else
      stop'defineCentrality error\n\n'
      endif
      read(line(i+22:i+22),*)jtable(jj)
      read(line(i+25:i+25),*)jzmode(jj)
      if(line(i+27:i+27).eq.'B')then
      read(line(i+28:i+28),*)i2
      jzmode(jj)=i2*10+jzmode(jj)
      endif
      !print*,'+++++++',jj,jtable(jj),jzmode(jj)
      endif
      if(jj.eq.3)call getMcentr
      if(zlimit(100).gt.0..and.jj.eq.3)then !Mcentrality 
        if(jzmode(3).eq.2)then !C2=Npom
          if(zlimit(100).lt.500)stop'CentralityLimit 100 too small'
        elseif(jzmode(3).eq.1)then !C1=bim
          if(zlimit(100).gt.500)stop'CentralityLimit 100 too big '
        endif
      endif 

          elseif(line(i+4:i+15).eq.'defineRanges')then

      if(nopen.ne.-1)then       !only first read
      if(line(i:i+1).ne.'ZJ')stop'ERROR 31122013b ###############'
      read(line(i+18:i+18),*)jtable(7)
      endif

          elseif(line(i:i+8).eq.'execFemto')then

      if(irootcproot.ne.10)then
      if(nopen.ne.-1)then       !only first read
      read(line(i+11:i+11),*)k1
      read(line(i+14:i+14),*)k2
      read(line(i+17:i+17),*)k3
      ifemto=k1*100+k2*10+k3
      ix=17
      if(i+ix+2.lt.j)then
      if(line(i+ix+2:i+ix+2).eq.'H')then
      read(line(i+ix+3:i+ix+3),*)kk
      if(kk.lt.1.or.kk.gt.3)stop'execFemto error 1\n\n'
      ihdim(kk)=1
      if(i+ix+5.lt.j)then
      if(line(i+ix+5:i+ix+5).eq.'H')then
      read(line(i+ix+6:i+ix+6),*)kk
      if(kk.lt.1.or.kk.gt.3)stop'execFemto error 2\n\n'
      ihdim(kk)=1
      if(i+ix+8.lt.j)then
      if(line(i+ix+8:i+ix+8).eq.'H')then
      read(line(i+ix+9:i+ix+9),*)kk
      if(kk.lt.1.or.kk.gt.3)stop'execFemto error 3\n\n'
      ihdim(kk)=1
      endif
      endif
      endif
      endif
      endif
      endif
      !print*,'+++++++',ifemto,ihdim
      endif
      endif

          elseif(line(i:i+8).eq.'orderTree')then

      if(irootcproot.eq.1)then
      if(nopen.ne.-1)then               !only first read
      call setOrderTree(line,i,j,irootcproot,maxpom,igrpac,nfifac)
      call orderTree(maxpom)
      endif
      endif

          elseif(line(i:i+6).eq.'execDih')then

      if(irootcproot.ne.10)then
      if(nopen.ne.-1)then               !only first read
      if(fndt(nfndt:nfndt).ne.'X')then
      read(line(i+9:i+9),*)k1
      read(line(i+12:i+12),*)k2
      idih=k1*10+k2
      ix=12
      if(i+ix+2.lt.j)then
      read(line(i+ix+2:i+ix+2),*)kk
      idih=idih+100*kk
      endif
      if(i+ix+4.lt.j)then
      read(line(i+ix+4:i+ix+4),*)kk
      idih=idih+1000*kk
      endif
      if(i+ix+6.lt.j)then
      read(line(i+ix+6:i+ix+6),*)kk
      idih=idih+10000*kk
      endif
      if(i+ix+8.lt.j)then
      read(line(i+ix+8:i+ix+8),*)kk
      idih=idih+100000*kk
      endif
      !call dih  !removed in 3444m4
      endif
      endif
      endif

          elseif(line(i:j).eq.'addplot')then

      iaddplot=iaddplot+1
      if(iaddplot.gt.mxaddplot1)stop'ERROR 29062013'
      call utword(line,i,j,0)
      jaddplot=0
      do while(line(i:j).ne.';')
      if(line(i:j).eq.'reset')then
        iaddplot=0
        do ii2=1,mxaddplot1
        do jj2=1,mxaddplot2
        caddplot(ii2,jj2)='                    '
        enddo
        enddo
      else
        jaddplot=jaddplot+1
        if(jaddplot.gt.mxaddplot2)stop'ERROR 29062013b'
        caddplot(iaddplot,jaddplot)(1:j-i+1)=line(i:j)
        maddplot(iaddplot)=jaddplot
      endif
      call utword(line,i,j,0)
      enddo

          elseif(line(i:i+6).eq.'execVmd')then

      if(nopen.ne.-1)then       !only first read
      read(line(i+9:i+9),*)k1
      ivmd=k1
      !print*,'+++++++',ivmd
      endif

          elseif(line(i:i+6).eq.'getTree')then

      call setGetTree(line,i,j,irootcproot,igrTree,muTree)

          elseif(line(i+4:i+13).eq.'defineBins')then

      if(nopen.ne.-1)then       !only first read
      ixin=i+15
      if(line(i+14:i+14).eq.'%')ixin=i+16
      read(line(i+1:i+1),*)ii
      if(ii.gt.mrclass)stop'mrclass too small\n\n'
      if(irclass(ii).ne.0)stop'Ki in use\n\n'
      irclass(ii)=1
      imo=0
      if(line(i+14:i+14).eq.'%')imo=1
      if(line(i+14:i+14).eq.'(')imo=2
      ix2=ixin-2
      val2=0
      nrclass(ii)=nrclass(ii)-1
      do
        if(ix2+2.gt.j)exit
        ix1=ix2+2
        if(line(ix1-1:ix1+1).eq.'___')ix1=ix1+2
        val1=val2
        ix2=ix1
        do while(line(ix2+1:ix2+1).ne.','.and.line(ix2+1:ix2+1).ne.')'
     .    .and.line(ix2+1:ix2+1).ne.';'.and.line(ix2+1:ix2+3).ne.'___')
        ix2=ix2+1
        enddo
        read(line(ix1:ix2),*)val2
        if(line(ix1-1:ix1-1).ne.';'.and.line(ix1-3:ix1-1).ne.'___')then
          nrclass(ii)=nrclass(ii)+1
          if(nrclass(ii).gt.0)then
          if(imo.eq.1)then
          rclass(ii,1,nrclass(ii))=zlimit(nint(val1))
          rclass(ii,2,nrclass(ii))=zlimit(nint(val2))
          kclass(ii,1,nrclass(ii))=nint(val1)
          kclass(ii,2,nrclass(ii))=nint(val2)
          else
          rclass(ii,1,nrclass(ii))=val1
          rclass(ii,2,nrclass(ii))=val2
          endif
          !print*,'+++++++++',ii,imo,nrclass(ii),
          !.rclass(ii,1,nrclass(ii))     ,rclass(ii,2,nrclass(ii))
          rclass(ii,3,nrclass(ii))=0
          endif
        endif
      enddo
      endif

           elseif(line(i+4:i+12).eq.'histoBins')then

      if(nopen.eq.-1)then       !only second read
      ixin=i+15
      if(line(i+13:i+13).ne.'%')stop'\n\n ERROR 30092012\n\n'
      read(line(i+1:i+1),*)ii
      read(line(ixin:j),'(a)')txt
      write(ifhi,'(a)')   '!--------------------------------'
      write(ifhi,'(a)')   '!          histoBins             '
      write(ifhi,'(a)')   '!--------------------------------'
      write(ifhi,'(a)')  'openhisto name histoBins  xrange 0 1'
      write(ifhi,'(a,i1,2x,3a)')'text 0.3 0.9 ""B',ii
     .,'(',line(ixin:j-1),')""'
      do i=1,nrclass(ii)
        if(i.le.6)then
          xi=0.05
          yi=0.8-i*0.1
        else
          xi=0.50
          yi=1.4-i*0.1
        endif
        if(kclass(ii,2,i).eq.100)then
        write(ifhi,'(a,2f4.1,a,i2,a,i4,a,i3,a)')
     . 'text',xi,yi,' ""',i,':',kclass(ii,1,i),'-',kclass(ii,2,i),'%""'
        else
        write(ifhi,'(a,2f4.1,a,i2,a,i4,a,i2,a)')
     . 'text',xi,yi,' ""',i,':',kclass(ii,1,i),'-',kclass(ii,2,i),'%""'
        endif
      enddo
      write(ifhi,'(a)')  'closehisto plot 0'

      endif

           elseif(line(i:i+6).eq.'histo->')then

      if(line(i+7:j).eq.'HydroEpsilon')       stop'histo-> obsolete'
      if(line(i+7:j).eq.'HydroEpsilon2')      stop'histo-> obsolete'
      if(line(i+7:j).eq.'HydroEntropy')       stop'histo-> obsolete'
      if(line(i+7:j).eq.'HydroTemperature')   stop'histo-> obsolete'
      if(line(i+7:j).eq.'HydroRadVelocity')   stop'histo-> obsolete'
      if(line(i+7:j).eq.'HydroLongVelocity')  stop'histo-> obsolete'
      if(line(i+7:j).eq.'HydroAverages')      stop'histo-> obsolete'
      if(line(i+7:j).eq.'HydroBaryon')        stop'histo-> obsolete'
      if(line(i+7:j).eq.'HydroFoVol')         stop'histo-> obsolete'
      if(line(i+7:j).eq.'HydroFoRadius')      stop'histo-> obsolete'
      if(line(i+7:j).eq.'HydroFoRadVelocity') stop'histo-> obsolete'
      if(line(i+7:j).eq.'HydroFoTangVelocity')stop'histo-> obsolete'
      if(line(i+7:j).eq.'HydroFoBarmu')       stop'histo-> obsolete'
      if(line(i+7:j).eq.'HydroEpsilonEta')    stop'histo-> obsolete'
      if(line(i+7:j).eq.'HydroBaryonEta')     stop'histo-> obsolete'
      if(line(i+7:j).eq.'HydroEntropyEta')    stop'histo-> obsolete'
      if(line(i+7:j).eq.'HydroEos')           stop'histo-> obsolete'

           elseif(line(i:j).eq.'idcode')then

      call utword(line,i,j,0)

           elseif(line(i:j).eq.'istuse')then

      call utword(line,i,j,0)

           elseif(line(i:j).eq.'xpara')then

      call utword(line,i,j,0)
      call utword(line,i,j,0)

           elseif(line(i:i+5).eq.'xparas'.or.line(i:i+2).eq.'xps')then

      if(line(i:j).eq.'xparas'.or.line(i:j).eq.'xps')then
        call utword(line,i,j,1)
        if(line(i:j).eq.'+10'.or.line(i:j).eq.'+15'
     .  .or.line(i:j).eq.'+20'.or.line(i:j).eq.'+25'
     .  .or.line(i:j).eq.'+30'.or.line(i:j).eq.'+35'
     .  .or.line(i:j).eq.'+40'.or.line(i:j).eq.'+45')then
          call utword(line,i,j,1)
        endif
      else
        call utword(line,i,j,1)
      endif 
      read(line(i:j),*)ipara
      ii=1
      do while (ii.le.ipara)
       call utword(line,i,j,1)
       if(line(i:i).eq.'>')then
        read(line(i+1:j),*)ii
        ii=ii+1
        call utword(line,i,j,1)
       endif 
       ii=ii+1
      enddo

           elseif(line(i:i+10).eq.'histoweight'
     .          .or.line(i:j).eq.'hw')then

      ncontr=0
      if(line(j-4:j).eq.'contr')then
       call utword(line,i,j,0)
       read(line(i:j),*)val
       ncontr=nint(val)
      endif
      if(nopen.eq.-1)then
      if(ncontr.eq.0)then
        write(ifhi,'(a,e22.14)')'histoweight ',histoweight
      else
        write(ifhi,'(a,e22.14)')'histoweight ',histowy(ncontr)
      endif
      endif

           elseif(line(i:j).eq.'histoweight-1'
     .          .or.line(i:j).eq.'hw-1')then

      if(nopen.eq.-1)then
      write(ifhi,'(a,e22.14)')'histoweight -1'
      endif

           elseif(line(i:j).eq.'input')then

      call setInput(line,i,j,nopen,iprmpt,iopcnt)

           elseif(line(i:j).eq.'xinput')then

      call setXinput(line,i,j,nopen,iprmpt,fnnx,nfnnx)

           elseif(line(i:j).eq.'nodecays')then

      call utword(line,i,j,0)
      do while (line(i:j).ne.'end')
       if(nrnody.ge.mxnody)then
        write(ifmt,'(a)')'too many nodecays; command ignored'
       else
        nrnody=nrnody+1
        read(line(i:j),*)val
        nody(nrnody)=nint(val)
       endif
      call utword(line,i,j,0)
      enddo

           elseif(line(i:j).eq.'nodecay')then

      call utword(line,i,j,0)
      if(nopen.ne.-1)then
      if(nrnody.ge.mxnody)then
      write(ifmt,'(a)')'too many nodecay commands; command ignored'
      j=1000
      i=j+1
      goto1
      endif
      nrnody=nrnody+1
      read(line(i:j),*)val
      nody(nrnody)=nint(val)
      endif

           elseif(line(i:j).eq.'dodecay')then

      call utword(line,i,j,0)
      if(nopen.ne.-1)then
      read(line(i:j),*)val
      idx=nint(val)
      nrn=0
      imv=0
      do while(nrn.lt.nrnody)
       nrn=nrn+1
       if(idx.eq.nody(nrn))then
         nrnody=nrnody-1
         imv=1
       endif
       if(imv.eq.1.and.nrn.le.nrnody)nody(nrn)=nody(nrn+1)
      enddo
      endif

           elseif(line(i:j).eq.'MinDecayLength')then

      call utword(line,i,j,0)
      if(nopen.ne.-1)then
      read(line(i:j),*)val
      ctaumin=val
      endif

           elseif(line(i:j).eq.'print')then

      call utworn(line,j,ne)
      call utword(line,i,j,0)
      if(line(i:j).eq.'monitor')then
        call utword(line,i,j,0)
        read(line(i:j),*)val
        jprint=nint(val)
      elseif(line(i:j).eq.'*')then
        nrpri=0
        call utword(line,i,j,0)
        read(line(i:j),*)val
        ish=nint(val)
      else
        nrpri=nrpri+1
        subpri(nrpri)='                    '
        subpri(nrpri)(1:j-i+1)=line(i:j)
        call utword(line,i,j,0)
        read(line(i:j),*)val
        ishpri(nrpri)=nint(val)
      endif

           elseif(line(i:j).eq.'printcheck')then

      call utworn(line,j,ne)
      if(ne.eq.0.and.iprmpt.gt.0)write(ifmt,'(a)')'screen or file?'
      call utword(line,i,j,0)
      if(line(i:j).eq.'screen')ifch=ifmt
      if(line(i:j).eq.'file')ifch=ifcx
      if(line(i:j).ne.'screen'.and.line(i:j).ne.'file')
     *write(ifmt,'(a)')'invalid option; command ignored'

           elseif(line(i:j).eq.'prompt')then

      call utworn(line,j,ne)
      if(ne.eq.0.and.iprmpt.gt.0)write(ifmt,'(a)')'on or off or auto?'
      call utword(line,i,j,0)
      if(line(i:j).eq.'on')iprmpt=2
      if(line(i:j).eq.'off')iprmpt=-2
      if(line(i:j).eq.'auto'.and.nopen.eq.0)iprmpt=1
      if(line(i:j).eq.'auto'.and.nopen.eq.1)iprmpt=-1
      if(line(i:j).ne.'on'.and.line(i:j).ne.'off'
     *.and.line(i:j).ne.'auto')stop'invalid option'

           elseif(line(i:j).eq.'return')then

      if(nopen.ne.-1)then
        close(ifop)
        nopen=nopen-1
        if(nopen.eq.0.and.iprmpt.eq.-1)iprmpt=1
ccc      close(20+nopen)
ccc      nopen=nopen-1
ccc      if(nopen.eq.0.and.iprmpt.eq.-1)iprmpt=1
      endif

           elseif(line(i:j).eq.'run')then

           elseif(line(i:j).eq.'rivet')then

      if(nopen.eq.-1)then       !only second read
       call setRivet(line,i,j)
      else
       call utword(line,i,j,0)
       call utword(line,i,j,0)
       if(line(j+2:j+7).eq.'calib:')call utword(line,i,j,0)
      endif

           elseif(line(i:j).eq.'runprogram')then

      if(nskip.lt.0)then
        jjtb=mod(-nskip,10)
        jjeos=-nskip/10
      endif
      if(fndt(nfndt:nfndt).eq.'X')then
        if(ishxxx.ne.0)stop'\n\n  ERROR 07072012a \n\n'
        if(iecho .ne.0)stop'\n\n  ERROR 07072012b \n\n'
        if(itit  .ne.0)stop'\n\n  ERROR 07072012c \n\n'
        ixx=1
        if(jcentrality.lt.0)then
        ixx=nint(zclass(3,icentrality))
        if(izmode.eq.2)ixx=1
        if(nskip.eq.1)ixx=1
        if(irootcproot.eq.1)ixx=1
        endif
        if(nskip.lt.0)ixx=jjtb
        write(ifdt,*)ixx
        open(17,file=fndt(1:nfndt-1)//'Y',status='unknown')
        iyy=ncentrality
        if(nskip.lt.0)iyy=1
        write(17,*)iyy
        close (17)
        open(17,file=fndt(1:nfndt-1)//'Z',status='unknown')
        write(17,*)igrTree
        close (17)
        open(17,file=fndt(1:nfndt-1)//'W',status='unknown')
        write(17,*)nfifac
        close (17)
        stop
      endif
      if(kchopen.eq.0.and.fnch(1:nfnch).ne.'none')then
        open(unit=ifcx,file=fnch(1:nfnch),status='unknown')
        kchopen=1
      endif
      if(khiopen.eq.0.and.fnhi(1:nfnhi).ne.'none')then
        open(unit=ifhi,file=fnhi(1:nfnhi),status='unknown')
        khiopen=1
      endif
      if(kdtopen.eq.0.and.fndt(1:nfndt).ne.'none')then
        open(unit=ifdt,file=fndt(1:nfndt),status='unknown')
        kdtopen=1
      endif
      if(kcpopen.eq.0.and.fncp(1:nfncp).ne.'none')then
        open(unit=ifcp,file=fncp(1:nfncp),status='unknown')
        kcpopen=1
      endif
      call clop(1)
      if(nopen.ne.-1)then
        return     !only for first read
      else
        if(ieof.eq.1)write(ifhi,'(a)')'!EOF'
        do i=1,9
        if(ieof.eq.1.and.ifhix(i).gt.0)write(ifhix(i),'(a)')'!EOF'
        enddo
        close(unit=ifcx)
        close(unit=ifhi)
        close(unit=ifdt)
        call RootClose()
        call histo2EndFigure()
        return
      endif

           elseif(line(i:j).eq.'if')then

      call utword(line,i,j,0)
      ix=i
      jx=j
      linex=line
      call utword(line,i,j,0)
      ifset=1
      read(line(i:j),*)val1
      call utword(line,i,j,0)
      read(line(i:j),*)val2
      if(linex(ix:jx).eq.'engy')then
        call idmass(idproj,amproj)
        call idmass(idtarg,amtarg)
        xxengy=0.
        if(engy.gt.0.)then
          xxengy=engy
        elseif(ecms.gt.0.)then
          xxengy=ecms
        elseif(elab.gt.0)then
          xxengy=sqrt( 2*elab*amtarg+amtarg**2+amproj**2 )
        elseif(pnll.gt.0)then
          xxengy=sqrt( 2*sqrt(pnll**2+amproj**2)*amtarg
     *                 +amtarg**2+amproj**2 )
        elseif(ekin.gt.0.)then
          xxelab=ekin+amproj
          xxengy=sqrt( 2*xxelab*amtarg+amtarg**2+amproj**2 )
        endif
        if(xxengy.lt.val1.or.xxengy.gt.val2)ifset=0
      elseif(linex(ix:jx).eq.'projtarg')then
        if(maproj.ne.nint(val1).or.matarg.ne.nint(val2))ifset=0
      elseif(linex(ix:jx).eq.'minmass')then
        if(min(maproj,matarg).lt.nint(val1)
     . .or.min(maproj,matarg).gt.nint(val2))ifset=0
      endif

           elseif(line(i:j).eq.'set'.and.ifset.eq.1)then

      call utword(line,i,j,0)
      linex=line
      ix=i
      jx=j
      call utword(line,i,j,0)
      if(linex(ix:jx).ne.'centrality'
     ..and.linex(ix:jx-1).ne.'hdtext')then
      read(line(i:j),*)val
      endif
c       general
      if(linex(ix:jx).eq.'model') model=nint(val)
      if(linex(ix:jx).eq.'iquasiel') iquasiel=nint(val)
      if(linex(ix:jx).eq.'iversn')iversn=nint(val)
      if(linex(ix:jx).eq.'iappl' )iappl=nint(val)
      if(linex(ix:jx).eq.'gefac')gefac=val
      if(linex(ix:jx).eq.'nevent')nevent=nint(val)
      if(linex(ix:jx).eq.'nfull') nfull=nint(val)
      if(linex(ix:jx).eq.'nfull'.and.gefac.ne.1.)then
        nfull=max(1,nint(nfull*gefac))
        if(nopen.ne.-1)
     .  write(ifmt,'(a,f5.2,a,i5)')'nfull multiplied by',gefac
     .    ,' -> ',nfull
      endif
      if(linex(ix:jx).eq.'nfull'.and.iffsigiEfficiency.eq.1)then
        if(iopcnt.eq.1)then
          nfull=0
          if(nopen.ne.-1)write(ifmt,'(a)')'nfull set to zero'
        endif
        !print*,'set nfull: iopcnt nfull',iopcnt,nfull
      endif
      if(linex(ix:jx).eq.'nfreeze')nfreeze=nint(val)
      if(linex(ix:jx).eq.'ninicon')ninicon=nint(val)
      if(linex(ix:jx).eq.'rapcms')rapcms=sngl(val)
      if(linex(ix:jx).eq.'s0min' )s0min=sngl(val)
      if(linex(ix:jx).eq.'egymin' )egymin=sngl(val)
      if(linex(ix:jx).eq.'egymax' )egymax=sngl(val)
csp addition for jet-medium interaction
      if(linex(ix:jx).eq.'medium' )medium=sngl(val)
c       constants
      if(linex(ix:jx).eq.'ainfin')ainfin=sngl(val)
c       printout options
      if(linex(ix:jx).eq.'iprmpt')iprmpt=nint(val)
      if(linex(ix:jx).eq.  'ish' )ish=nint(val)
      if(linex(ix:jx).eq.'ishsub')ishsub=nint(val)
      if(linex(ix:jx).eq.'irandm')irandm=nint(val)
      if(linex(ix:jx).eq.'irewch')irewch=nint(val)
      if(linex(ix:jx).eq.'iecho ')iecho =nint(val)
      if(linex(ix:jx).eq.'modsho')modsho=nint(val)
      if(linex(ix:jx).eq.'modsho')modshox=modsho
      if(linex(ix:jx).eq.'idensi')idensi=nint(val)
      if(linex(ix:jx).eq.'infragm')infragm=nint(val)
      if(linex(ix:jx).eq.'ibreit')ibreit=nint(val)
      if(linex(ix:jx).eq.'ifoele')ifoele=nint(val)
      if(linex(ix:jx).eq.'ishevt')ishevt=nint(val)
      if(linex(ix:jx).eq.'iwseed')iwseed=nint(val)
      if(linex(ix:jx).eq.'jwseed')jwseed=nint(val)
      if(linex(ix:jx).eq.'ihepmc')ihepmc=nint(val)
      if(linex(ix:jx).eq.'ihepframe')ihepframe=nint(val)
c       fragmentation and decay
      if(linex(ix:jx).eq.'pdiqua')pdiqua=sngl(val)
      if(linex(ix:jx).eq.'pud'   )pud   =sngl(val)
      if(linex(ix:jx).eq.'pmqu'  )pmqu  =sngl(val)
      if(linex(ix:jx).eq.'pmqd'  )pmqd  =sngl(val)
      if(linex(ix:jx).eq.'pmqs ' )pmqs  =sngl(val)
      if(linex(ix:jx).eq.'pmqc ' )pmqc  =sngl(val)
      if(linex(ix:jx).eq.'pmqb ' )pmqb  =sngl(val) !bg
      if(linex(ix:jx).eq.'pmqq ' )then
                                  pmqq  =sngl(val)
                                  qmass(0)=pmqq
      endif
      if(linex(ix:jx).eq.'fkappa' )fkappa   =sngl(val)
      if(linex(ix:jx).eq.'fkappag' )fkappag   =sngl(val)
      if(linex(ix:jx).eq.'pudd ' )pudd  =sngl(val)
      if(linex(ix:jx).eq.'puds ' )puds  =sngl(val)
      if(linex(ix:jx).eq.'pudc ' )pudc  =sngl(val)
      if(linex(ix:jx).eq.'pudb ' )pudb  =sngl(val)  !bg
      if(linex(ix:jx).eq.'strcut' )strcut   =sngl(val)
      if(linex(ix:jx).eq.'diqcut' )diqcut   =sngl(val)
      if(linex(ix:jx).eq.'ioptf ')ioptf =nint(val)
      if(linex(ix:jx).eq.'delrex')delrex=sngl(val)
      if(linex(ix:jx).eq.'ndecay')ndecay=nint(val)
      if(linex(ix:jx).eq.'maxres')maxres=nint(val)
      if(linex(ix:jx).eq.'pbreak')pbreak=sngl(val)
      if(linex(ix:jx).eq.'pbreakg')pbreakg=sngl(val)
      if(linex(ix:jx).eq.'zetacut')zetacut=sngl(val)
      if(linex(ix:jx).eq.'ptfra')ptfra=sngl(val)
      if(linex(ix:jx).eq.'ptfraqq')ptfraqq=sngl(val)
      if(linex(ix:jx).eq.'aouni ')aouni=sngl(val)
c       lepton-nucleon and e+e-
      if(linex(ix:jx).eq.'iolept')iolept=nint(val)
      if(linex(ix:jx).eq.'ydmin')ydmin=sngl(val)
      if(linex(ix:jx).eq.'ydmax')ydmax=sngl(val)
      if(linex(ix:jx).eq.'qdmin')qdmin=sngl(val)
      if(linex(ix:jx).eq.'qdmax')qdmax=sngl(val)
      if(linex(ix:jx).eq.'themin')themin=sngl(val)
      if(linex(ix:jx).eq.'themax')themax=sngl(val)
      if(linex(ix:jx).eq.'elomin')elomin=sngl(val)
      if(linex(ix:jx).eq.'engy' )engy=sngl(val)
      if(linex(ix:jx).eq.'elab' )elab=sngl(val)
      if(linex(ix:jx).eq.'ekin' )ekin=sngl(val)
      if(linex(ix:jx).eq.'ecms' )ecms=sngl(val)
      if(linex(ix:jx).eq.'ebeam' )ebeam=sngl(val)
      if(linex(ix:jx).eq.'elepti')elepti=sngl(val)
      if(linex(ix:jx).eq.'elepto')elepto=sngl(val)
      if(linex(ix:jx).eq.'angmue')angmue=sngl(val)
      if(linex(ix:jx).eq.'noebin')noebin=nint(val)
      if(linex(ix:jx).eq.'engmin')engmin=sngl(val)
      if(linex(ix:jx).eq.'engmax')engmax=sngl(val)
      if(linex(ix:jx).eq.'iologe')iologe=nint(val)
      if(linex(ix:jx).eq.'iologl')iologl=nint(val)
      if(linex(ix:jx).eq.'itflav')itflav=nint(val)
      if(linex(ix:jx).eq.'idisco')idisco=nint(val)
c       hadron-hadron
      if(linex(ix:jx).eq.'ioTestFact')ioTestFact=nint(val)
      if(linex(ix:jx).eq.'idhprTestFact')idhprTestFact=nint(val)
      if(linex(ix:jx).eq.'ihqTestFact')ihqTestFact=nint(val)
      if(linex(ix:jx).eq.'laddTestFact')laddTestFact=nint(val)
      if(linex(ix:jx).eq.'noptTestFact')noptTestFact=nint(val)
      if(linex(ix:jx).eq.'iffsigiEfficiency')iffsigiEfficiency=nint(val) 
      if(linex(ix:jx).eq. 'pnll' )pnll=sngl(val)
      if(linex(ix:jx).eq.'idproj')idprojin=nint(val)
      idproj=idprojin
      if(linex(ix:jx).eq.'idtarg')idtargin=nint(val)
      idtarg=idtargin
      if(idtarg.eq.0)idtarg=1120
      if(linex(ix:jx).eq.'ptq   ')ptq   =sngl(val)
      if(linex(ix:jx).eq.'rstrau(1)')rstrau(1)=sngl(val)
      if(linex(ix:jx).eq.'rstrad(1)')rstrad(1)=sngl(val)
      if(linex(ix:jx).eq.'rstras(1)')rstras(1)=sngl(val)
      if(linex(ix:jx).eq.'rstrac(1)')rstrac(1)=sngl(val)
      if(linex(ix:jx).eq.'rstrau(2)')rstrau(2)=sngl(val)
      if(linex(ix:jx).eq.'rstrad(2)')rstrad(2)=sngl(val)
      if(linex(ix:jx).eq.'rstras(2)')rstras(2)=sngl(val)
      if(linex(ix:jx).eq.'rstrac(2)')rstrac(2)=sngl(val)
      if(linex(ix:jx).eq.'rstrau(3)')rstrau(3)=sngl(val)
      if(linex(ix:jx).eq.'rstrad(3)')rstrad(3)=sngl(val)
      if(linex(ix:jx).eq.'rstras(3)')rstras(3)=sngl(val)
      if(linex(ix:jx).eq.'rstrac(3)')rstrac(3)=sngl(val)
c      if(linex(ix:jx).eq.'rstrau(4)')rstrau(4)=sngl(val)
c      if(linex(ix:jx).eq.'rstrad(4)')rstrad(4)=sngl(val)
c      if(linex(ix:jx).eq.'rstras(4)')rstras(4)=sngl(val)
c      if(linex(ix:jx).eq.'rstrac(4)')rstrac(4)=sngl(val)
      if(linex(ix:jx).eq.'rstrasi')rstrasi=sngl(val)
      if(linex(ix:jx).eq.'wgtval')wgtval=sngl(val)
      if(linex(ix:jx).eq.'wgtsea')wgtsea=sngl(val)
      if(linex(ix:jx).eq.'wgtdiq')wgtdiq=sngl(val)
      if(linex(ix:jx).eq.'wgtqqq(1)')wgtqqq(1)=sngl(val)
      if(linex(ix:jx).eq.'wgtqqq(2)')wgtqqq(2)=sngl(val)
      if(linex(ix:jx).eq.'wgtqqq(3)')wgtqqq(3)=sngl(val)
c      if(linex(ix:jx).eq.'wgtqqq(4)')wgtqqq(4)=sngl(val)
      if(linex(ix:jx).eq.'wproj ')wproj =sngl(val)
      if(linex(ix:jx).eq.'wtarg ')wtarg =sngl(val)
      if(linex(ix:jx).eq.'wexcit')wexcit=sngl(val)
c      if(linex(ix:jx).eq.'cutmsq')cutmsq=sngl(val)
      if(linex(ix:jx).eq.'cutmss')cutmss=sngl(val)
      if(linex(ix:jx).eq.'exmass')exmass=sngl(val)
      if(linex(ix:jx).eq.'iregge')iregge=nint(val)
      if(linex(ix:jx).eq.'isopom')isopom=nint(val)
      if(linex(ix:jx).eq.'ishpom')ishpom=nint(val)
      if(linex(ix:jx).eq.'iscreen')iscreen=nint(val)
      if(linex(ix:jx).eq.'iq2sat')iq2sat=nint(val)
      if(linex(ix:jx).eq.'iseahot')iseahot=nint(val)
      if(linex(ix:jx).eq.'irmdrop')irmdrop=nint(val)
      if(linex(ix:jx).eq.'nprmax')nprmax=nint(val)
      if(linex(ix:jx).eq.'intpol')intpol=nint(val)
      if(linex(ix:jx).eq.'isigma')isigma=nint(val)
      if(linex(ix:jx).eq.'iomega')iomega=nint(val)
      if(linex(ix:jx).eq.'isetcs')isetcs=nint(val)
      if(linex(ix:jx).eq.'iemsb' )iemsb= nint(val)
      if(linex(ix:jx).eq.'iemspm')iemspm=nint(val)
      if(linex(ix:jx).eq.'iemspx')iemspx=nint(val)
      if(linex(ix:jx).eq.'iemsse')iemsse=nint(val)
      if(linex(ix:jx).eq.'iemsdr')iemsdr=nint(val)
      if(linex(ix:jx).eq.'iemsrx')iemsrx=nint(val)
      if(linex(ix:jx).eq.'iemsi2')iemsi2=nint(val)
      if(linex(ix:jx).eq.'iemsi1')iemsi1=nint(val)
      if(linex(ix:jx).eq.'ioems' )ioems= nint(val)
      if(linex(ix:jx).eq.'ipsahot')ipsahot= nint(val)
      if(linex(ix:jx).eq.'ispacetime')ispacetime= nint(val)
c       unified parameters
      if(linex(ix:jx).eq.'iclpro1')iclpro1=nint(val)
      if(linex(ix:jx).eq.'iclpro2')iclpro2=nint(val)
      if(linex(ix:jx).eq.'icltar1')icltar1=nint(val)
      if(linex(ix:jx).eq.'icltar2')icltar2=nint(val)
      if(linex(ix:jx).eq.'iclegy1')iclegy1=nint(val)
      if(linex(ix:jx).eq.'iclegy2')iclegy2=nint(val)
      if(linex(ix:jx).eq.'egylow')egylow=sngl(val)
      if(linex(ix:jx).eq.'alppom')alppom=sngl(val)
      if(linex(ix:jx).eq.'gamhad(1)')stop'gamhad(1) not allowed'
      if(linex(ix:jx).eq.'gamhad(2)')stop'gamhad(2) not allowed'
      if(linex(ix:jx).eq.'gamhad(3)')stop'gamhad(3) not allowed'
      if(linex(ix:jx).eq.'gamhads(1)')gamhadsi(1)=sngl(val)
      if(linex(ix:jx).eq.'gamhads(2)')gamhadsi(2)=sngl(val)
      if(linex(ix:jx).eq.'gamhads(3)')gamhadsi(3)=sngl(val)
c      if(linex(ix:jx).eq.'gamhads(4)')gamhadsi(4)=sngl(val)
      if(linex(ix:jx).eq.'gampar')gampar=sngl(val)
      if(linex(ix:jx).eq.'gamtil')gamtil=sngl(val)
      if(linex(ix:jx).eq.'slopom')slopom=sngl(val)
      if(linex(ix:jx).eq.'slopoms')slopoms=sngl(val)
      if(linex(ix:jx).eq.'r2part' )r2part= sngl(val)
      if(linex(ix:jx).eq.'r2har(1)' )r2har(1)= sngl(val)
      if(linex(ix:jx).eq.'r2har(2)' )r2har(2)= sngl(val)
      if(linex(ix:jx).eq.'r2har(3)' )r2har(3)= sngl(val)
c      if(linex(ix:jx).eq.'r2har(4)' )r2har(4)= sngl(val)
      if(linex(ix:jx).eq.'r2had(1)' )r2had(1)= sngl(val)
      if(linex(ix:jx).eq.'r2had(2)' )r2had(2)= sngl(val)
      if(linex(ix:jx).eq.'r2had(3)' )r2had(3)= sngl(val)
c      if(linex(ix:jx).eq.'r2had(4)' )r2had(4)= sngl(val)
      if(linex(ix:jx).eq.'r2hads(1)' )r2hads(1)= sngl(val)
      if(linex(ix:jx).eq.'r2hads(2)' )r2hads(2)= sngl(val)
      if(linex(ix:jx).eq.'r2hads(3)' )r2hads(3)= sngl(val)
c      if(linex(ix:jx).eq.'r2hads(4)' )r2hads(4)= sngl(val)
      if(linex(ix:jx).eq.'r3pom'   )r3pom= sngl(val)
      if(linex(ix:jx).eq.'egyscr'  )egyscr= sngl(val)
      if(linex(ix:jx).eq.'epscrw'  )epscrw= sngl(val)
      if(linex(ix:jx).eq.'epscrx'  )epscrxi= sngl(val)
      if(linex(ix:jx).eq.'zbrads'  )zbrads= sngl(val)
      if(linex(ix:jx).eq.'zfluct'  )zfluct= sngl(val)
      if(linex(ix:jx).eq.'epscrs'  )epscrs= sngl(val)
      if(linex(ix:jx).eq.'epscrh'  )epscrh= sngl(val)
      if(linex(ix:jx).eq.'epscrp'  )epscrp= sngl(val)
      if(linex(ix:jx).eq.'epscrb'  )epscrb= sngl(val)
      if(linex(ix:jx).eq.'znurho'  )znurho= sngl(val)
      if(linex(ix:jx).eq.'epscrd'  )epscrd= sngl(val)
      if(linex(ix:jx).eq.'q2sft'   )q2sft= sngl(val)
      if(linex(ix:jx).eq.'zdfbmx'  )zdfbmx= sngl(val)
      if(linex(ix:jx).eq.'epscrg'  )epscrg= sngl(val)
      if(linex(ix:jx).eq.'wdiff(1)')wdiff(1)=sngl(val)
      if(linex(ix:jx).eq.'wdiff(2)')wdiff(2)=sngl(val)
      if(linex(ix:jx).eq.'wdiff(3)')wdiff(3)=sngl(val)
c      if(linex(ix:jx).eq.'wdiff(4)')wdiff(4)=sngl(val)
      if(linex(ix:jx).eq.'facdif')  facdif=sngl(val)
      if(linex(ix:jx).eq.'facmc')   facmc=sngl(val)
      if(linex(ix:jx).eq.'rexddf')  rexddf=sngl(val)
      if(linex(ix:jx).eq.'rexndi(1)')rexndi(1)=sngl(val)
      if(linex(ix:jx).eq.'rexndi(2)')rexndi(2)=sngl(val)
      if(linex(ix:jx).eq.'rexndi(3)')rexndi(3)=sngl(val)
c      if(linex(ix:jx).eq.'rexndi(4)')rexndi(4)=sngl(val)
      if(linex(ix:jx).eq.'alphigh(1)')alphigh(1)=sngl(val)
      if(linex(ix:jx).eq.'alphigh(2)')alphigh(2)=sngl(val)
      if(linex(ix:jx).eq.'rexpdif(1)')rexpdif(1)=sngl(val)
      if(linex(ix:jx).eq.'rexpdif(2)')rexpdif(2)=sngl(val)
      if(linex(ix:jx).eq.'rexpdif(3)')rexpdif(3)=sngl(val)
c      if(linex(ix:jx).eq.'rexpdif(4)')rexpdif(4)=sngl(val)
      if(linex(ix:jx).eq.'rexdif(1)')rexdif(1)=sngl(val)
      if(linex(ix:jx).eq.'rexdif(2)')rexdif(2)=sngl(val)
      if(linex(ix:jx).eq.'rexdif(3)')rexdif(3)=sngl(val)
c      if(linex(ix:jx).eq.'rexdif(4)')rexdif(4)=sngl(val)
      if(linex(ix:jx).eq.'rexres(1)')rexres(1)=sngl(val)
      if(linex(ix:jx).eq.'rexres(2)')rexres(2)=sngl(val)
      if(linex(ix:jx).eq.'rexres(3)')rexres(3)=sngl(val)
c      if(linex(ix:jx).eq.'rexres(4)')rexres(4)=sngl(val)
      if(linex(ix:jx).eq.'zrhoinc'  )zrhoinc=sngl(val)
      if(linex(ix:jx).eq.'alpreg(1)')alpreg(1)=sngl(val)
      if(linex(ix:jx).eq.'alpreg(2)')alpreg(2)=sngl(val)
      if(linex(ix:jx).eq.'alpreg(3)')alpreg(3)=sngl(val)
      if(linex(ix:jx).eq.'sloreg(1)')sloreg(1)=sngl(val)
      if(linex(ix:jx).eq.'sloreg(2)')sloreg(2)=sngl(val)
      if(linex(ix:jx).eq.'sloreg(3)')sloreg(3)=sngl(val)
      if(linex(ix:jx).eq.'gamreg(1,-1)')gamreg(1,-1)=sngl(val)
      if(linex(ix:jx).eq.'gamreg(1,0)')gamreg(1,0)=sngl(val)
      if(linex(ix:jx).eq.'gamreg(1,1)')gamreg(1,1)=sngl(val)
      if(linex(ix:jx).eq.'gamreg(2,-1)')gamreg(2,-1)=sngl(val)
      if(linex(ix:jx).eq.'gamreg(2,0)')gamreg(2,0)=sngl(val)
      if(linex(ix:jx).eq.'gamreg(2,1)')gamreg(2,1)=sngl(val)
      if(linex(ix:jx).eq.'gamreg(3,-1)')gamreg(3,-1)=sngl(val)
      if(linex(ix:jx).eq.'gamreg(3,0)')gamreg(3,0)=sngl(val)
      if(linex(ix:jx).eq.'gamreg(3,1)')gamreg(3,1)=sngl(val)
      if(linex(ix:jx).eq.'r2reg(1,-1)' )r2reg(1,-1)= sngl(val)
      if(linex(ix:jx).eq.'r2reg(1,0)' )r2reg(1,0)= sngl(val)
      if(linex(ix:jx).eq.'r2reg(1,1)' )r2reg(1,1)= sngl(val)
      if(linex(ix:jx).eq.'r2reg(2,-1)' )r2reg(2,-1)= sngl(val)
      if(linex(ix:jx).eq.'r2reg(2,0)' )r2reg(2,0)= sngl(val)
      if(linex(ix:jx).eq.'r2reg(2,1)' )r2reg(2,1)= sngl(val)
      if(linex(ix:jx).eq.'r2reg(3,-1)' )r2reg(3,-1)= sngl(val)
      if(linex(ix:jx).eq.'r2reg(3,0)' )r2reg(3,0)= sngl(val)
      if(linex(ix:jx).eq.'r2reg(3,1)' )r2reg(3,1)= sngl(val)
c      if(linex(ix:jx).eq.'amhdibar')amhdibar= sngl(val)
      if(linex(ix:jx).eq.'ptdiff')ptdiff=sngl(val)
      if(linex(ix:jx).eq.'alppar')alppar=sngl(val)
      if(linex(ix:jx).eq.'alpparh')alpparh=sngl(val)
      if(linex(ix:jx).eq.'alpsea')alpsea=sngl(val)
      if(linex(ix:jx).eq.'alpval')alpval=sngl(val)
      if(linex(ix:jx).eq.'alpdiq')alpdiq=sngl(val)
      if(linex(ix:jx).eq.'alplea(1)')alplea(1)=sngl(val)
      if(linex(ix:jx).eq.'alplea(2)')alplea(2)=sngl(val)
      if(linex(ix:jx).eq.'alplea(3)')alplea(3)=sngl(val)
c      if(linex(ix:jx).eq.'alplea(4)')alplea(4)=sngl(val)
      if(linex(ix:jx).eq.'alpdif')alpdif=sngl(val)
      if(linex(ix:jx).eq.'alpdi(1)')alpdi(1)=sngl(val)
      if(linex(ix:jx).eq.'alpdi(2)')alpdi(2)=sngl(val)
      if(linex(ix:jx).eq.'alpndi(1)')alpndi(1)=sngl(val)
      if(linex(ix:jx).eq.'alpndi(2)')alpndi(2)=sngl(val)
      if(linex(ix:jx).eq.'ammsqq')ammsqq=sngl(val)
      if(linex(ix:jx).eq.'ammsqd')ammsqd=sngl(val)
      if(linex(ix:jx).eq.'ammsdd')ammsdd=sngl(val)
      if(linex(ix:jx).eq.'cumpox')cumpox=sngl(val)
      if(linex(ix:jx).eq.'xcupom')xcupom=sngl(val)
      if(linex(ix:jx).eq.'xcutsf')xcutsf=sngl(val)
      if(linex(ix:jx).eq.'zzsoft')zzsoft=sngl(val)
      if(linex(ix:jx).eq.'ycupom')ycupom=sngl(val)
      if(linex(ix:jx).eq.'reminv')reminv=sngl(val)
      if(linex(ix:jx).eq.'ptsend')ptsend=sngl(val)
      if(linex(ix:jx).eq.'ptsendi')ptsendi=sngl(val)
      if(linex(ix:jx).eq.'ptsems')ptsems=sngl(val)
      if(linex(ix:jx).eq.'ptipom')ptipom=sngl(val)
      if(linex(ix:jx).eq.'ptipos')ptipos=sngl(val)
      if(linex(ix:jx).eq.'ptiposi')ptiposi=sngl(val)
      if(linex(ix:jx).eq.'facq2tim')facq2tim=sngl(val)
      if(linex(ix:jx).eq.'zdrinc')zdrinc=sngl(val)
      if(linex(ix:jx).eq.'zmsinc')zmsinc=sngl(val)
      if(linex(ix:jx).eq.'ediflim')ediflim=sngl(val)
      if(linex(ix:jx).eq.'edifac')edifac=sngl(val)
      if(linex(ix:jx).eq.'edmaxi')edmaxi=sngl(val)
      if(linex(ix:jx).eq.'epmaxi')stop'\n\n not supported \n\n'
      if(linex(ix:jx).eq.'zopinc')zopinc=sngl(val)
      if(linex(ix:jx).eq.'zipinc')zipinc=sngl(val)
      if(linex(ix:jx).eq.'fkainc')fkainc=sngl(val)
      if(linex(ix:jx).eq.'fkaadd')fkaadd=sngl(val)
      if(linex(ix:jx).eq.'zodinc')zodinc=sngl(val)
      if(linex(ix:jx).eq.'zbrmax')zbrmax=sngl(val)
      if(linex(ix:jx).eq.'zoeinc')zoeinc=sngl(val)
      if(linex(ix:jx).eq.'ptipomi')ptipomi=sngl(val)
      if(linex(ix:jx).eq.'xmxmas')xmxmas=sngl(val)
      if(linex(ix:jx).eq.'ptsecu')ptsecu=sngl(val)
      if(linex(ix:jx).eq.'zdfinc')zdfinc=sngl(val)
      if(linex(ix:jx).eq.'zbcut')zbcut=sngl(val)
      if(linex(ix:jx).eq.'xzcut') xzcut=sngl(val)
      if(linex(ix:jx).eq.'zxymin') zxymin=sngl(val)
      if(linex(ix:jx).eq.'zxyinc') zxyinc=sngl(val)
      if(linex(ix:jx).eq.'xminremn')xminremn=sngl(val)
      if(linex(ix:jx).eq.'xmindiff')xmindiff=sngl(val)
      if(linex(ix:jx).eq.'alpdro(1)')alpdro(1)=sngl(val)
      if(linex(ix:jx).eq.'alpdro(2)')alpdro(2)=sngl(val)
      if(linex(ix:jx).eq.'alpdro(3)')alpdro(3)=sngl(val)
      if(linex(ix:jx).eq.'amdrmax')amdrmax=sngl(val)
      if(linex(ix:jx).eq.'amdrmin')amdrmin=sngl(val)
      if(linex(ix:jx).eq.'iodiba')iodiba=nint(val)
      if(linex(ix:jx).eq.'bidiba')bidiba=sngl(val)
      if(linex(ix:jx).eq.'disize')disize=sngl(val)
      if(linex(ix:jx).eq.'facpdf4')   call pdfparamset(1, sngl(val) )
      if(linex(ix:jx).eq.'facpdf4sat')call pdfparamset(2, sngl(val) )
      if(linex(ix:jx).eq.'facpdf5')   call pdfparamset(3, sngl(val) )
      if(linex(ix:jx).eq.'facpdf5sat')call pdfparamset(4, sngl(val) )

c       hard pomeron parameters
      if(linex(ix:jx).eq.'ioangord')ioangord=nint(val)
      if(linex(ix:jx).eq.'coekaf')coekaf=sngl(val)
      if(linex(ix:jx).eq.'q2nmin')q2nmin=sngl(val)
      if(linex(ix:jx).eq.'q2ini' )q2ini=sngl(val)
      if(linex(ix:jx).eq.'q2fin' )q2fin=sngl(val)
      if(linex(ix:jx).eq.'qcdlam')qcdlam=sngl(val)
      if(linex(ix:jx).eq.'betpom')betpomi=sngl(val)
      if(linex(ix:jx).eq.'glusea')glusea=sngl(val)
      if(linex(ix:jx).eq.'factk' )factk=sngl(val)
      if(linex(ix:jx).eq.'naflav')naflav=nint(val)
      if(linex(ix:jx).eq.'nbflav')nbflav=nint(val)
      if(linex(ix:jx).eq.'nrflav')nrflav=nint(val)
      if(linex(ix:jx).eq.'facnof')facnof=sngl(val)
      if(linex(ix:jx).eq.'pt2cut')pt2cut=sngl(val)
      if(linex(ix:jx).eq.'factgam')factgam=sngl(val)
      if(linex(ix:jx).eq.'factjpsi')factjpsi=sngl(val)
      if(linex(ix:jx).eq.'factsat')factsat=sngl(val)
      if(linex(ix:jx).eq.'alpsft')alpsft=sngl(val)
      if(linex(ix:jx).eq.'alpsat')alpsat=sngl(val)
      if(linex(ix:jx).eq.'ptphot' )ptphot=sngl(val)
      if(linex(ix:jx).eq.'delh')delh=sngl(val)
c       nucleus-nucleus
      if(linex(ix:jx).eq.'iokoll')iokoll=nint(val)
      if(linex(ix:jx).eq.'laproj')laproj=nint(val)
      if(linex(ix:jx).eq.'maproj')maproj=nint(val)
      if(linex(ix:jx).eq.'latarg')latarg=nint(val)
      if(linex(ix:jx).eq.'matarg')matarg=nint(val)
      if(linex(ix:jx).eq.'core'  )core  =sngl(val)
c      if(linex(ix:jx).eq.'ncolmx')ncolmx=nint(val)
      if(linex(ix:jx).eq.'fctrmx')fctrmx=sngl(val)
      if(linex(ix:jx).eq.'bmaxim')bmaxim=sngl(val)
      if(linex(ix:jx).eq.'bminim')bminim=sngl(val)
      if(linex(ix:jx).eq.'bref80')bref80=sngl(val)
      if(linex(ix:jx).eq.'phimax')phimax=sngl(val)
      if(linex(ix:jx).eq.'phimin')phimin=sngl(val)
      if(linex(ix:jx).eq.'nq2mnfix')nq2mnfix=nint(val)
      if(linex(ix:jx).eq.'asatur')asatur=sngl(val)
      if(linex(ix:jx).eq.'bsatur')bsatur=sngl(val)
      if(linex(ix:jx).eq.'csatur')csatur=sngl(val)
      if(linex(ix:jx).eq.'dsatur')dsatur=sngl(val)
      if(linex(ix:jx).eq.'esatur')esatur=sngl(val)
c       rescattering parameters
      if(linex(ix:jx).eq.'iorsce')iorsce=nint(val)
      if(linex(ix:jx).eq.'iorsdf')iorsdf=nint(val)
      if(linex(ix:jx).eq.'iorshh')iorshh=nint(val)
      if(linex(ix:jx).eq.'iocluin')iocluin=nint(val)
      if(linex(ix:jx).eq.'ioquen')ioquen=nint(val)
      if(linex(ix:jx).eq.'hacore')hacore=sngl(val)
      if(linex(ix:jx).eq.'amimfs')amimfs=sngl(val)
      if(linex(ix:jx).eq.'amimel')amimel=sngl(val)
      if(linex(ix:jx).eq.'cepara')cepara=sngl(val)
      if(linex(ix:jx).eq.'dscale')dscale=sngl(val)
      if(linex(ix:jx).eq.'iceopt')iceopt=nint(val)
      if(linex(ix:jx).eq.'delamf')delamf=sngl(val)
      if(linex(ix:jx).eq.'deuamf')deuamf=sngl(val)
      if(linex(ix:jx).eq.'taustr')taustr=sngl(val)
      if(linex(ix:jx).eq.'tauhac')tauhac=sngl(val)
      if(linex(ix:jx).eq.'iostro')iostro=nint(val)
      if(linex(ix:jx).eq.'radeft1')radeft1=sngl(val)
      if(linex(ix:jx).eq.'radeft2')radeft2=sngl(val)
      if(linex(ix:jx).eq.'nsegsu')nsegsu=nint(val)
      if(linex(ix:jx).eq.'dsegce')dsegce=sngl(val)
      if(linex(ix:jx).eq.'kigrid')kigrid=nint(val)
      if(linex(ix:jx).eq.'fxcell')fxcell=sngl(val)
      if(linex(ix:jx).eq.'fzcell')fzcell=sngl(val)
      if(linex(ix:jx).eq.'ptlow') ptlow=sngl(val)
      if(linex(ix:jx).eq.'ptupp') ptupp=sngl(val)
      if(linex(ix:jx).eq.'qufac') qufac=sngl(val)
      if(linex(ix:jx).eq.'quexpo')quexpo=sngl(val)
      if(linex(ix:jx).eq.'ijetfluid')ijetfluid=nint(val)
      if(linex(ix:jx).eq.'cutdxy')cutdxy=sngl(val)
      if(linex(ix:jx).eq.'fludiq')fludiq=sngl(val)
      if(linex(ix:jx).eq.'epscri(1)')epscri(1)=sngl(val)
      if(linex(ix:jx).eq.'epscri(3)')epscri(3)=sngl(val)
      if(linex(ix:jx).eq.'rapcol')rapcol=sngl(val)
      if(linex(ix:jx).eq.'leadcore')leadcore=nint(val)
      if(linex(ix:jx).eq.'amsiac')amsiac=sngl(val)
      if(linex(ix:jx).eq.'amprif')amprif=sngl(val)
      if(linex(ix:jx).eq.'delvol')delvol=sngl(val)
      if(linex(ix:jx).eq.'deleps')deleps=sngl(val)
      if(linex(ix:jx).eq.'tauzer')tauzer=sngl(val)
      if(linex(ix:jx).eq.'tempoico')tauzer=sngl(val)
      if(linex(ix:jx).eq.'deltau')deltau=sngl(val)
      if(linex(ix:jx).eq.'numtau')numtau=nint(val)
      if(linex(ix:jx).eq.'ihacas')ihacas=int(val)
      if(linex(ix:jx).eq.'iuelast')iuelast=int(val)
      if(linex(ix:jx).eq.'iuskip')iuskip=int(val)
      if(linex(ix:jx).eq.'iuchaskip')stop'\n\n 13092012\n\n'
      if(linex(ix:jx).eq.'iunostring')iunostring=int(val)
      if(linex(ix:jx).eq.'dlzeta')dlzeta=sngl(val)
      if(linex(ix:jx).eq.'etafac')etafac=sngl(val)
      if(linex(ix:jx).eq.'facnuc')facnuc=sngl(val)
c       ico
      if(linex(ix:jx).eq.'cutico') cutico=sngl(val)
      if(linex(ix:jx).eq.'sgzico') sgzico=sngl(val)
      if(linex(ix:jx).eq.'dssico') dssico=sngl(val)
      if(linex(ix:jx).eq.'jcorona')jcorona=int(val)
      if(linex(ix:jx).eq.'icocore')icocore=int(val)
      if(linex(ix:jx).eq.'nsmooth')nsmooth=int(val)
      if(linex(ix:jx).eq.'iranphi')iranphi=int(val)
      if(linex(ix:jx).eq.'icospec')icospec=int(val)
      if(linex(ix:jx).eq.'icoremn')icoremn=int(val)
      if(linex(ix:jx).eq.'icostrg')icostrg=int(val)
      if(linex(ix:jx).eq.'icotabm')icotabm=int(val)
      if(linex(ix:jx).eq.'icotabr')icotabr=int(val)
      if(linex(ix:jx).eq.'nxico')nxico=nint(val)
      if(linex(ix:jx).eq.'nyico')nyico=nint(val)
      if(linex(ix:jx).eq.'nzico')nzico=nint(val)
      if(linex(ix:jx).eq.'xminico')xminico=sngl(val)
      if(linex(ix:jx).eq.'xmaxico')xmaxico=sngl(val)
      if(linex(ix:jx).eq.'yminico')yminico=sngl(val)
      if(linex(ix:jx).eq.'ymaxico')ymaxico=sngl(val)
      if(linex(ix:jx).eq.'zminico')zminico=sngl(val)
      if(linex(ix:jx).eq.'zmaxico')zmaxico=sngl(val)
      if(linex(ix:jx).eq.'ioicoplot')ioicoplot=nint(val)
c       phsd
      if(linex(ix:jx).eq.'iphsd')iphsd=nint(val)
c       eos
      if(linex(ix:jx).eq.'ieostabm')ieostabm=int(val)
c       spherio
      if(linex(ix:jx).eq.'ispherio')ispherio=nint(val)
      if(linex(ix:jx).eq.'jspherio')jspherio=nint(val)
c       hlle
      if(linex(ix:jx).eq.'ihlle')ihlle=nint(val)
      if(linex(ix:jx).eq.'tfrout')tfrout=sngl(val)
      if(linex(ix:jx).eq.'kfrout')kfrout=nint(val)
      if(linex(ix:jx).eq.'nfrout')nfrout=nint(val)
      if(linex(ix:jx).eq.'efrout')efrout=sngl(val)
      if(linex(ix:jx).eq.'fofac') fofac=sngl(val)
      if(linex(ix:jx).eq.'fofai') fofai=sngl(val)
      if(linex(ix:jx).eq.'ntaumx')ntaumx=nint(val)
      if(linex(ix:jx).eq.'ioeos')ioeos=nint(val)
      if(linex(ix:jx).eq.'iochem')iochem=nint(val)
      if(linex(ix:jx).eq.'iozerof')iozerof=nint(val)
      if(linex(ix:jx).eq.'etaos') etaos=sngl(val)
      if(linex(ix:jx).eq.'zetaos')zetaos=sngl(val)
c       hydro
      if(linex(ix:jx).eq.'tauone')tauone=sngl(val)
      if(linex(ix:jx).eq.'tautwo')tautwo=sngl(val)
      if(linex(ix:jx).eq.'tauup')tauup=sngl(val)
      if(linex(ix:jx).eq.'floini') floini=sngl(val)
      if(linex(ix:jx).eq.'nxhy')nxhy=nint(val)
      if(linex(ix:jx).eq.'nyhy')nyhy=nint(val)
      if(linex(ix:jx).eq.'nzhy')nzhy=nint(val)
      if(linex(ix:jx).eq.'xminhy')xminhy=sngl(val)
      if(linex(ix:jx).eq.'xmaxhy')xmaxhy=sngl(val)
      if(linex(ix:jx).eq.'yminhy')yminhy=sngl(val)
      if(linex(ix:jx).eq.'ymaxhy')ymaxhy=sngl(val)
      if(linex(ix:jx).eq.'zminhy')zminhy=sngl(val)
      if(linex(ix:jx).eq.'zmaxhy')zmaxhy=sngl(val)
      if(linex(ix:jx).eq.'ifaahlle')ifaahlle=nint(val)
      if(linex(ix:jx).eq.'ifathlle')ifathlle=nint(val)
      if(linex(ix:jx).eq.'ifazhlle')ifazhlle=nint(val)
      if(linex(ix:jx).eq.'epsfin')epsfin=sngl(val)
      if(linex(ix:jx).eq.'dtauhy')dtauhy=sngl(val)
      if(linex(ix:jx).eq.'volex') volex= sngl(val)
      if(linex(ix:jx).eq.'fdtau') fdtau= sngl(val)
      if(linex(ix:jx).eq.'itabptl')itabptl=nint(val)
      if(linex(ix:jx).eq.'maxsurf')maxsurf=nint(val)
      if(linex(ix:jx).eq.'ienvar')ienvar=nint(val)
      if(linex(ix:jx).eq.'izmode')izmode=nint(val)
c       droplet decay
      if(linex(ix:jx).eq.'dezzer')dezzer=sngl(val)
      if(linex(ix:jx).eq.'ioclude')ioclude=nint(val)
      if(linex(ix:jx).eq.'amuseg')amuseg=sngl(val)
      if(linex(ix:jx).eq.'yradmx')yradmx=sngl(val)
      if(linex(ix:jx).eq.'yradpp')yradpp=sngl(val)
      if(linex(ix:jx).eq.'yradpi')yradpi=sngl(val)
      if(linex(ix:jx).eq.'yradpx')yradpx=sngl(val)
      if(linex(ix:jx).eq.'facecc')facecc=sngl(val)
      if(linex(ix:jx).eq.'rcoll' )rcoll= sngl(val)
      if(linex(ix:jx).eq.'bag4rt')bag4rt=sngl(val)
      if(linex(ix:jx).eq.'taurem')taurem=sngl(val)
c       droplet specification
      if(linex(ix:jx).eq. 'keu'  )keu=nint(val)
      if(linex(ix:jx).eq. 'ked'  )ked=nint(val)
      if(linex(ix:jx).eq. 'kes'  )kes=nint(val)
      if(linex(ix:jx).eq. 'kec'  )kec=nint(val)
      if(linex(ix:jx).eq. 'keb'  )keb=nint(val)
      if(linex(ix:jx).eq. 'ket'  )ket=nint(val)
      if(linex(ix:jx).eq. 'tecm' )tecm=sngl(val)
      if(linex(ix:jx).eq. 'volu' )volu=sngl(val)
      if(linex(ix:jx).eq. 'vrad' )vrad=sngl(val)
      if(linex(ix:jx).eq. 'facts')facts=sngl(val)
      if(linex(ix:jx).eq. 'factb')factb=sngl(val)
      if(linex(ix:jx).eq. 'factq')factq=sngl(val)
      if(linex(ix:jx).eq. 'yrrope')call pfe6paramset(1, sngl(val) )
      if(linex(ix:jx).eq. 'ylrope')call pfe6paramset(2, sngl(val) )
      if(linex(ix:jx).eq. 'facmicro')facmicro=sngl(val)
      if(linex(ix:jx).eq.'inbxxx')inbxxx=sngl(val)
c       metropolis
      if(linex(ix:jx).eq.'fitermet')fitermet=sngl(val)
      if(linex(ix:jx).eq.'felamet')felamet=sngl(val)
      if(linex(ix:jx).eq.'iospec')iospec=nint(val)
      if(linex(ix:jx).eq.'iocova')iocova=nint(val)
      if(linex(ix:jx).eq.'iopair')iopair=nint(val)
      if(linex(ix:jx).eq.'iozero')iozero=nint(val)
      if(linex(ix:jx).eq.'ioflac')ioflac=nint(val)
      if(linex(ix:jx).eq.'iostat')iostat=nint(val)
      if(linex(ix:jx).eq.'ioinco')ioinco=nint(val)
      if(linex(ix:jx).eq.'iograc')iograc=nint(val)
      if(linex(ix:jx).eq.'epsgc' )epsgc=sngl(val)
      if(linex(ix:jx).eq.'iocite')iocite=nint(val)
      if(linex(ix:jx).eq.'ioceau')ioceau=nint(val)
      if(linex(ix:jx).eq.'iociau')iociau=nint(val)
      if(linex(ix:jx).eq.'ioinct')ioinct=nint(val)
      if(linex(ix:jx).eq.'ioinfl')ioinfl=nint(val)
      if(linex(ix:jx).eq.'iowidn')iowidn=nint(val)
      if(linex(ix:jx).eq.'ionlat')ionlat=nint(val)
      if(linex(ix:jx).eq.'iomom')iomom=nint(val)
      if(linex(ix:jx).eq.'ioobsv')ioobsv=nint(val)
      if(linex(ix:jx).eq.'iosngl')iosngl=nint(val)
      if(linex(ix:jx).eq.'iorejz')iorejz=nint(val)
      if(linex(ix:jx).eq.'iompar')iompar=nint(val)
      if(linex(ix:jx).eq.'iozinc')iozinc=nint(val)
      if(linex(ix:jx).eq.'iozevt')iozevt=nint(val)
      if(linex(ix:jx).eq. 'nadd' )nadd=nint(val)
      if(linex(ix:jx).eq.'iterma')iterma=nint(val)
      if(linex(ix:jx).eq.'itermx')stop'STOP: set iterma, not itermx'
      if(linex(ix:jx).eq.'iterpr')iterpr=nint(val)
      if(linex(ix:jx).eq.'iterpl')iterpl=nint(val)
      if(linex(ix:jx).eq.'iternc')iternc=nint(val)
      if(linex(ix:jx).eq.'epsr'  )epsr=sngl(val)
      if(linex(ix:jx).eq.'keepr' )keepr=nint(val)
c       strangelets
      if(linex(ix:jx).eq.'iopenu')iopenu=nint(val)
      if(linex(ix:jx).eq.'themas')themas=sngl(val)
c       tests
      if(linex(ix:jx).eq.'iotst1')iotst1=nint(val)
      if(linex(ix:jx).eq.'iotst2')iotst2=nint(val)
      if(linex(ix:jx).eq.'iotst3')iotst3=nint(val)
      if(linex(ix:jx).eq.'iotst4')iotst4=nint(val)
      if(linex(ix:jx).eq.'vparam(1)')vparam(1)=sngl(val)
      if(linex(ix:jx).eq.'vparam(2)')vparam(2)=sngl(val)
      if(linex(ix:jx).eq.'vparam(3)')vparam(3)=sngl(val)
      if(linex(ix:jx).eq.'vparam(4)')vparam(4)=sngl(val)
      if(linex(ix:jx).eq.'vparam(5)')vparam(5)=sngl(val)
      if(linex(ix:jx).eq.'vparam(6)')vparam(6)=sngl(val)
      if(linex(ix:jx).eq.'vparam(7)')vparam(7)=sngl(val)
      if(linex(ix:jx).eq.'vparam(8)')vparam(8)=sngl(val)
      if(linex(ix:jx).eq.'vparam(9)')vparam(9)=sngl(val)
      if(linex(ix:jx).eq.'vparam(10)')vparam(10)=sngl(val)
      if(linex(ix:jx).eq.'vparam(11)')vparam(11)=sngl(val)
      if(linex(ix:jx).eq.'vparam(12)')vparam(12)=sngl(val)
      if(linex(ix:jx).eq.'vparam(13)')vparam(13)=sngl(val)
      if(linex(ix:jx).eq.'vparam(14)')vparam(14)=sngl(val)
      if(linex(ix:jx).eq.'vparam(15)')vparam(15)=sngl(val)
      if(linex(ix:jx).eq.'vparam(16)')vparam(16)=sngl(val)
      if(linex(ix:jx).eq.'vparam(17)')vparam(17)=sngl(val)
      if(linex(ix:jx).eq.'vparam(18)')vparam(18)=sngl(val)
      if(linex(ix:jx).eq.'vparam(19)')vparam(19)=sngl(val)
      if(linex(ix:jx).eq.'vparam(20)')vparam(20)=sngl(val)
c       jpsi, qkonia (quarkonia)
      if(linex(ix:jx).eq.'jpsi  ')jpsi  =nint(val)
      if(linex(ix:jx).eq.'jpsifi')jpsifi=nint(val)
      if(linex(ix:jx).eq.'sigj  ')sigj  =sngl(val)
      if(linex(ix:jx).eq.'taumx ')taumx =sngl(val)
      if(linex(ix:jx).eq.'nsttau')nsttau=nint(val)
      if(linex(ix:jx).eq.'ijphis')ijphis=nint(val)
      if(linex(ix:jx).eq.'qkonia1')qkonia1=sngl(val)
      if(linex(ix:jx).eq.'qkonia2')qkonia2=sngl(val)
      if(linex(ix:jx).eq.'qkonia3')qkonia3=sngl(val)
      if(linex(ix:jx).eq.'qkonia4')qkonia4=sngl(val)

c       analysis: intermittency, space-time, droplets, formation time
      if(linex(ix:jx).eq.'ymximi')ymximi=sngl(val)
      if(linex(ix:jx).eq.'imihis')imihis=nint(val)
      if(linex(ix:jx).eq.'isphis')isphis=nint(val)
      if(linex(ix:jx).eq.'iologb')iologb=nint(val)
      if(linex(ix:jx).eq.'ispall')ispall=nint(val)
      if(linex(ix:jx).eq.'wtmini')wtmini=sngl(val)
      if(linex(ix:jx).eq.'wtstep')wtstep=sngl(val)
      if(linex(ix:jx).eq.'iwcent')iwcent=nint(val)
      if(linex(ix:jx).eq.'iclhis')iclhis=nint(val)
      if(linex(ix:jx).eq.'iwtime')iwtime=nint(val)
      if(linex(ix:jx).eq.'wtimet')wtimet=sngl(val)
      if(linex(ix:jx).eq.'wtimei')wtimei=sngl(val)
      if(linex(ix:jx).eq.'wtimea')wtimea=sngl(val)
c       core
      if(linex(ix:jx).eq.'feloss1')call core2paramset(1, sngl(val) )
      if(linex(ix:jx).eq.'feloss2')call core2paramset(2, sngl(val) )
      if(linex(ix:jx).eq.'feloss3')call core2paramset(3, sngl(val) )
      if(linex(ix:jx).eq.'feloss4')call core2paramset(4, sngl(val) )
      if(linex(ix:jx).eq.'feloss5')call core2paramset(5, sngl(val) )
      if(linex(ix:jx).eq.'feloss6')call core2paramset(6, sngl(val) )
      if(linex(ix:jx).eq.'facposz')facposz=sngl(val)
      if(linex(ix:jx).eq.'facposf')facposf=sngl(val)
      if(linex(ix:jx).eq.'tauzer1')tauzer1=sngl(val)
      if(linex(ix:jx).eq.'tauzer2')tauzer2=sngl(val)
      if(linex(ix:jx).eq.'nsegmincore')nsegmincore=nint(val)
      if(linex(ix:jx).eq.'ratiomaxcore')ratiomaxcore=sngl(val)
      if(linex(ix:jx).eq.'corcor(4)')
     .                stop'Use facposz,facposf instead of corcor(4)'
c       root
      if(linex(ix:jx).eq.'iextree')iextree=nint(val)
c       other
      if(linex(ix:jx).eq.'jjeos')jjeos =nint(val)
      if(linex(ix:jx).eq.'gaumx')gaumx =sngl(val)
      if(linex(ix:jx).eq.'nclean')nclean=nint(val)
      if(linex(ix:jx).eq.'istore')istore=nint(val)
      if(linex(ix:jx).eq.'ioidch')ioidch=nint(val)
      if(linex(ix:jx).eq.'iframe')iframe=nint(val)
      if(linex(ix:jx).eq.'jframe')jframe=nint(val)
      if(linex(ix:jx).eq.'labsys')stop'labsys no longer supported'
      if(linex(ix:jx).eq.'irescl')irescl=nint(val)
      if(linex(ix:jx).eq.'iremn')iremn=nint(val)
      if(linex(ix:jx).eq.'ifrade')ifrade=nint(val)
      if(linex(ix:jx).eq.'idecay')idecay=nint(val)
      if(linex(ix:jx).eq.'ihdecay')ihdecay=nint(val)
      if(linex(ix:jx).eq.'jdecay')jdecay=nint(val)
      if(linex(ix:jx).eq.'ntrymx')ntrymx=nint(val)
      if(linex(ix:jx).eq.'istmax')istmax=nint(val)
      if(linex(ix:jx).eq.'istfor')istfor=nint(val)
      if(linex(ix:jx).eq.'ionudi')ionudi=nint(val)
      if(linex(ix:jx).eq.'iotype')iotype=nint(val)
      if(linex(ix:jx).eq.'seedi') seedi =val
      if(linex(ix:jx).eq.'seedj') seedj =val
      if(linex(ix:jx).eq.'seedf') seedj2=val
      if(linex(ix:jx).eq.'ikolmn')ikolmn=nint(val)
      if(linex(ix:jx).eq.'ikolmx')ikolmx=nint(val)
      if(linex(ix:jx).eq.'nglmin')nglmin=nint(val)
      if(linex(ix:jx).eq.'nglmax')nglmax=nint(val)
      if(linex(ix:jx).eq.'ptrmin')ptrmin=val
      if(linex(ix:jx).eq.'ioecc')ioecc=nint(val)
      if(linex(ix:jx).eq.'valecc')valecc=sngl(val)
      if(linex(ix:jx).eq.'xvaria')xvaria='      '
      if(linex(ix:jx).eq.'xvaria')xvaria(1:j-i+1)=line(i:j)
      if(linex(ix:jx).eq.'yvaria')yvaria='      '
      if(linex(ix:jx).eq.'yvaria')yvaria(1:j-i+1)=line(i:j)
      if(linex(ix:jx).eq.'normal')normal=nint(val)
      if(linex(ix:jx).eq.'xminim')xminim=sngl(val)
      if(linex(ix:jx).eq.'xmaxim')xmaxim=sngl(val)
      if(linex(ix:jx).eq.'nrbins')nrbins=nint(val)
      if(linex(ix:jx).eq.'hisfac')hisfac=sngl(val)
      if(linex(ix:jx).eq.'xshift')xshift=sngl(val)
      if(linex(ix:jx).eq.'etacut')etacut=sngl(val)
      if(linex(ix:jx).eq.'xpar1' )xpar1=sngl(val)
      if(linex(ix:jx).eq.'xpar2' )xpar2=sngl(val)
      if(linex(ix:jx).eq.'xpar3' )xpar3=sngl(val)
      if(linex(ix:jx).eq.'xpar4' )xpar4=sngl(val)
      if(linex(ix:jx).eq.'xpar5' )xpar5=sngl(val)
      if(linex(ix:jx).eq.'xpar6' )xpar6=sngl(val)
      if(linex(ix:jx).eq.'xpar7' )xpar7=sngl(val)
      if(linex(ix:jx).eq.'xpar8' )xpar8=sngl(val)
      if(linex(ix:jx).eq.'xpar9' )xpar9=sngl(val)
      if(linex(ix:jx).eq.'xpar10')xpar10=sngl(val)
      if(linex(ix:jx).eq.'xpar11')xpar11=sngl(val)
      if(linex(ix:jx).eq.'xpar12')xpar12=sngl(val)
      if(linex(ix:jx).eq.'xpar13')xpar13=sngl(val)
      if(linex(ix:jx).eq.'xpar14')xpar14=sngl(val)
      if(linex(ix:jx).eq.'xpar15')xpar15=sngl(val)
      if(linex(ix:jx).eq.'xpar16')xpar16=sngl(val)
      if(linex(ix:jx).eq.'xpar17')xpar17=sngl(val)
      if(iscreen.ne.0)then
      if(linex(ix:jx).eq.'xpar98')xpar98=sngl(val)
      if(linex(ix:jx).eq.'xpar99')xpar99=sngl(val)
      endif
      if(linex(ix:jx).eq.'irdmpr')irdmpr=nint(val)
      if(linex(ix:jx).eq.'ilprtg')ilprtg=nint(val)
      if(linex(ix:jx).eq.'pytune')ipytune=nint(val)
c       hdtext
      if(linex(ix:jx-1).eq.'hdtext')then
      len=j-i+1
      len=min(len,30)
      if(linex(ix:jx).eq.'hdtext1')hdtext(1)(1:len)=line(i:i-1+len)
      if(linex(ix:jx).eq.'hdtext2')hdtext(2)(1:len)=line(i:i-1+len)
      if(linex(ix:jx).eq.'hdtext3')hdtext(3)(1:len)=line(i:i-1+len)
      if(linex(ix:jx).eq.'hdtext4')hdtext(4)(1:len)=line(i:i-1+len)
      if(linex(ix:jx).eq.'hdtext5')hdtext(5)(1:len)=line(i:i-1+len)
      if(linex(ix:jx).eq.'hdtext6')hdtext(6)(1:len)=line(i:i-1+len)
      if(linex(ix:jx).eq.'hdtext7')hdtext(7)(1:len)=line(i:i-1+len)
      if(linex(ix:jx).eq.'hdtext8')hdtext(8)(1:len)=line(i:i-1+len)
      endif
c       frame definitions
      if(linex(ix:jx).eq.'engy'.and.iframe.eq.0)iframe=11
      if(linex(ix:jx).eq.'ecms'.and.iframe.eq.0)iframe=11
      if(linex(ix:jx).eq.'elab'.and.iframe.eq.0)iframe=12
      if(linex(ix:jx).eq.'ekin'.and.iframe.eq.0)iframe=12
      if(linex(ix:jx).eq.'pnll'.and.iframe.eq.0)iframe=12
      if(linex(ix:jx).eq.'ebeam'.and.iframe.eq.0)iframe=21
      if(linex(ix:jx).eq.'noebin'.and.dabs(val-1.d0).gt.1.d-10)iframe=11
      if(linex(ix:jx).eq.'hydtab')hydt=line(i:j)
c       centrality
      if(linex(ix:jx).eq.'nskip')nskip=nint(val)
      if(linex(ix:jx).eq.'iopcnt')iopcnt=nint(val)
      if(linex(ix:jx).eq.'centrality')then
       if(nopen.ne.-1)then   !first read
        call setCentrality(line,i,j,iopcnt,ifmt
     .               ,icentrality,jcentrality,ffac,imax,ival,ishxxx)
        if(ffac.gt.1.000)then
          nfull=nfull*ffac
          modsho=modsho*ffac
          if(ishxxx.eq.1)
     .    write(ifmt,'(a,f6.2,a,i5)')'nfull multiplied by',ffac
     .    ,' -> ',nfull
          if(ishxxx.eq.1)
     .    write(ifmt,'(a,f6.2,a,i5)')'modsho multiplied by',ffac
     .    ,' -> ',modsho
        endif
        ncentrality=imax
        if(ival.ne.0)then
          if(izmode.eq.1)then
            bminim=zclass(1,icentrality)
            bmaxim=zclass(2,icentrality)
          elseif(izmode.eq.2)then
            ikolmn=zclass(1,icentrality)
            ikolmx=zclass(2,icentrality)
            if(zclass(4,icentrality).gt.-0.5)then
              bminim=zclass(4,icentrality)
              bmaxim=zclass(5,icentrality)
            endif
          elseif(izmode.eq.3)then
            nglmin=zclass(1,icentrality)
            nglmax=zclass(2,icentrality)
          elseif(izmode.eq.4)then
            segmin=zclass(1,icentrality)
            segmax=zclass(2,icentrality)
          endif
        endif
       else   !second read
        if(line(i:j).eq.'within')then
          call utword(line,i,j,0)
          if(line(i:j).ne.'{')stop'\n\n ERROR 19112011g\n\n'
          do while(line(i:j).ne.'}')
            call utword(line,i,j,0)
            if(line(i:j).eq.'2*')call utword(line,i,j,0)
          enddo
          if(line(i:j).ne.'}')stop'\n\n ERROR 19112011h\n\n'
        endif
       endif
      endif

           elseif(line(i:j).eq.'ifval')then

      isk=1
      call utword(line,i,j,0)
      if(line(i:j).eq.'centrality')then
        call utword(line,i,j,0)
        if(line(i:j).ne.'within')stop'\n\n ERROR 19112011c\n\n'
        call utword(line,i,j,0)
        if(line(i:j).ne.'{')stop'\n\n ERROR 19112011d\n\n'
        do while(line(i:j).ne.'}')
          call utword(line,i,j,0)
          if(line(i:j).ne.'}')then
          read(line(i:j),*)val1
          if(nint(val1).eq.icentrality)isk=0
          endif
        enddo
      elseif(line(i:j).eq.'rootcproot')then
        call utword(line,i,j,0)
        if(line(i:j).ne.'within')stop'\n\n ERROR 19112011e\n\n'
        call utword(line,i,j,0)
        if(line(i:j).ne.'{')stop'\n\n ERROR 19112011f\n\n'
        do while(line(i:j).ne.'}')
          call utword(line,i,j,0)
          if(line(i:j).ne.'}')then
          read(line(i:j),*)val1
          if(nint(val1).eq.irootcproot)isk=0
          endif
        enddo
      endif
      if(isk.eq.1)then
        if(nopen.ne.-1.and.iecho.ne.0)write(ifmt,'(a)')'SKIP'
        do while(line(i:j).ne.'endifval')
          call utword(line,i,j,0)
        enddo
        if(nopen.ne.-1.and.iecho.ne.0)write(ifmt,'(a)')'ENDSKIP'
      endif

           elseif(line(i:j).eq.'endifval')then

      continue

           elseif(line(i:j).eq.'set')then

      call utword(line,i,j,0)
      write(ifmt,'(2a)')line(i:j),' skipped'
      call utword(line,i,j,0)
      ifset=1

           elseif(line(i:j).eq.'kill')then

      call utword(line,i,j,0)
      write(ifmt,'(a)')'KILLED; ', line(i:j)
      stop'ERROR 29062011'
      ifset=1

           elseif(line(i:j).eq.'select')then !args: energyrange <valmin> <valmax>

      call utword(line,i,j,0)
      linex=line
      ix=i
      jx=j
      call utword(line,i,j,0)
      read(line(i:j),*)valmin
      call utword(line,i,j,0)
      read(line(i:j),*)valmax
      if(linex(ix:jx).eq.'energyrange')then
        val=engy
      else
        stop'ERROR wrong select argument'
      endif  
      jselect=1
      kselect=0
      if(val.ge.valmin.and.val.le.valmax)kselect=1
      if(kselect.eq.0)then!not selected
        call utword(line,i,j,0)
        do while(line(i+3:j).ne.'select')
          call utword(line,i,j,0)
        enddo
        if(line(i:j).eq.'endselect')jselect=0
      endif

           elseif(line(i:j).eq.'notselect')then

      if(jselect.eq.0)stop'ERROR select not active'
      if(kselect.eq.1)then!selected
        call utword(line,i,j,0)
        do while(line(i+3:j).ne.'select')
          call utword(line,i,j,0)
        enddo
        if(line(i:j).eq.'endselect')jselect=0
      endif

           elseif(line(i:j).eq.'endselect')then

        continue

           elseif(line(i:j).eq.'key')then

c key ... keyx in optns should be replaced 
c by s.th. like #ifBigSystem ... #fiBigSystem   <==============================
c where a "big system" is defined in the code,  <==============================
c based on maproj/matarg/engy

      call utword(line,i,j,0)
      linex=line
      ix=i
      jx=j
      call utword(line,i,j,0)
      if(ihlle.ne.99)then
        read(linex(ix:jx),*)val
        key=val
        read(line(i:j),*)val
        if(ecms.gt.0.)then
          ey=ecms
        elseif(elab.gt.0)then
          call idmass(idproj,amproj)
          call idmass(idtarg,amtarg)
          ey=sqrt( 2*elab*amtarg+amtarg**2+amproj**2 )
        else
          stop'\n\n ERROR 15062010B\n\n'
        endif
        isk=0
        w=mod(key,1d6)
        if(nint(ey).lt.nint(w)
     .  .or.nint(ey).gt.nint(val))isk=1
        key=(key-w)*1d-06
        w=mod(key,1d3)
        if(matarg.lt.w-20.or.matarg.gt.w+20)isk=1
        key=(key-w)*1d-3
        w=mod(key,1d2)
        if(latarg.lt.w-20.or.latarg.gt.w+20)isk=1
        key=(key-w)*1d-2
        w=mod(key,1d3)
        if(maproj.lt.w-20.or.maproj.gt.w+20)isk=1
        key=(key-w)*1d-3
        w=mod(key,1d3)
        if(laproj.lt.w-20.or.laproj.gt.w+20)isk=1
      else
        isk=1
      endif
      iskmin=min(iskmin,isk)
      iskkey=isk
      if(isk.eq.1)then
        do while(line(i:j).ne.'keyx')
          call utword(line,i,j,0)
        enddo
      endif

           elseif(line(i:j).eq.'keyi')then

      if(iskkey.ne.0.and.iskkey.ne.1)
     . stop'\n\n previous "key ... keyx" missing\n\n'
      if(iskkey.eq.0)then
        do while(line(i:j).ne.'keyx')
          call utword(line,i,j,0)
        enddo
      endif
      iskkey=2

           elseif(line(i:j).eq.'keyx')then

      continue

           elseif(line(i:j).eq.'stop')then  !same as return

      if(nopen.ne.-1)then
      close(20+nopen)
      nopen=nopen-1
      if(nopen.eq.0.and.iprmpt.eq.-1)iprmpt=1
      endif

           elseif(line(i:j).eq.'stopprogram')then

      !not used any more

           elseif(line(i:j).eq.'EndEposInput')then

      call clop(1)
      return

           elseif(line(i:j).eq.'string')then

      nstmax=nstmax+1
      ns=nstmax
      icinpu=0
      call utworn(line,j,ne)
      if(ne.eq.0.and.iprmpt.gt.0)
     *write(ifmt,'(a)')'string: prob icbac1 icbac2 icfor1 icfor2?'
      call utword(line,i,j,0)
      read(line(i:j),*)val
      prob(ns)=sngl(val)
      call utworn(line,j,ne)
      if(ne.eq.0.and.iprmpt.gt.0)
     *write(ifmt,'(a)')'string: icbac1 icbac2 icfor1 icfor2?'
      call utword(line,i,j,0)
      read(line(i:j),*)val
      icbac(ns,1)=nint(val)
      call utworn(line,j,ne)
      if(ne.eq.0.and.iprmpt.gt.0)
     *write(ifmt,'(a)')'string: icbac2 icfor1 icfor2?'
      call utword(line,i,j,0)
      read(line(i:j),*)val
      icbac(ns,2)=nint(val)
      call utworn(line,j,ne)
      if(ne.eq.0.and.iprmpt.gt.0)
     *write(ifmt,'(a)')'string: icfor1 icfor2?'
      call utword(line,i,j,0)
      read(line(i:j),*)val
      icfor(ns,1)=nint(val)
      call utworn(line,j,ne)
      if(ne.eq.0.and.iprmpt.gt.0)
     *write(ifmt,'(a)')'string: icfor2?'
      call utword(line,i,j,0)
      read(line(i:j),*)val
      icfor(ns,2)=nint(val)

           elseif(line(i:j).eq.'kinks')then

      nptl=0
      do k=1,4
        psum(k)=0    
      enddo
ctp290806      nrow=0
      nel=0
 10   continue
      call utword(line,i,j,0)
      if(line(i:j).eq.'endkinks')goto 12
      nel=nel+1
ctp290806      nrow=1+(nel-1)/4
      nc=mod(nel-1,4)+1
      read(line(i:j),*)a
      if(nc.eq.1)nptl=nptl+1
      if(nc.eq.1)idptl(nptl)=nint(a)
      if(nc.eq.2)pptl(1,nptl)=a
      if(nc.eq.3)pptl(2,nptl)=a
      if(nc.eq.4)then
        pptl(3,nptl)=a
        istptl(nptl)=20
        pptl(4,nptl)=sqrt(pptl(3,nptl)**2+pptl(2,nptl)**2
     $       +pptl(1,nptl)**2)
        do k=1,4
          psum(k)=psum(k)+pptl(k,nptl)
        enddo
      endif
      goto 10
 12   continue
      engy=sqrt(psum(4)**2-psum(1)**2-psum(2)**2-psum(3)**2)

           elseif(line(i:j).eq.'record')then

      call utworn(line,j,ne)
c      if(ne.eq.0.and.iprmpt.gt.0)
c     *     write(ifmt,'(a)')'kinks: icbac1 icbac2 icfor1 icfor2?'
      call utword(line,i,j,0)
      ir=0
      if(line(i:j).eq.'event')then
        ir=1
      elseif(line(i:j).eq.'particle')then
        ir=2
      else
        call utstop("Wrong definition for record!&")
      endif
      maxrec(ir)=0
 20   call utworn(line,j,ne)
c      if(ne.eq.0.and.iprmpt.gt.0)
c     *     write(6,'(a)')'<kinks-data (px-py-pz)>? (End=endkinks)'
      call utword(line,i,j,0)
      if(line(i:j).eq.'endrecord')then
         goto 22
      endif
      maxrec(ir)=maxrec(ir)+1
      irecty(maxrec(ir),ir)=-1
      if(ir.eq.1)then
        if(line(i:j).eq.'0') irecty(maxrec(ir),ir)=0
        if(line(i:j).eq.'nevt') irecty(maxrec(ir),ir)=1
        if(line(i:j).eq.'nptl') irecty(maxrec(ir),ir)=2
        if(line(i:j).eq.'b') irecty(maxrec(ir),ir)=3
        if(line(i:j).eq.'phi') irecty(maxrec(ir),ir)=4
        if(line(i:j).eq.'ncol') irecty(maxrec(ir),ir)=5
        if(line(i:j).eq.'pmx') irecty(maxrec(ir),ir)=6
        if(line(i:j).eq.'egy') irecty(maxrec(ir),ir)=7
        if(line(i:j).eq.'npj') irecty(maxrec(ir),ir)=8
        if(line(i:j).eq.'ntg') irecty(maxrec(ir),ir)=9
        if(line(i:j).eq.'npn') irecty(maxrec(ir),ir)=10
        if(line(i:j).eq.'npp') irecty(maxrec(ir),ir)=11
        if(line(i:j).eq.'ntn') irecty(maxrec(ir),ir)=12
        if(line(i:j).eq.'ntp') irecty(maxrec(ir),ir)=13
        if(line(i:j).eq.'jpn') irecty(maxrec(ir),ir)=14
        if(line(i:j).eq.'jpp') irecty(maxrec(ir),ir)=15
        if(line(i:j).eq.'jtn') irecty(maxrec(ir),ir)=16
        if(line(i:j).eq.'jtp') irecty(maxrec(ir),ir)=17
        if(line(i:j).eq.'amp') irecty(maxrec(ir),ir)=20
        if(line(i:j).eq.'amt') irecty(maxrec(ir),ir)=21
        if(line(i:j).eq.'qsq') irecty(maxrec(ir),ir)=22
        if(line(i:j).eq.'xbj') irecty(maxrec(ir),ir)=23
        if(line(i:j).eq.'typ') irecty(maxrec(ir),ir)=24
      else
        if(line(i:j).eq.'0') irecty(maxrec(ir),ir)=0
        if(line(i:j).eq.'i') irecty(maxrec(ir),ir)=1
        if(line(i:j).eq.'id') irecty(maxrec(ir),ir)=2
        if(line(i:j).eq.'p1') irecty(maxrec(ir),ir)=3
        if(line(i:j).eq.'p2') irecty(maxrec(ir),ir)=4
        if(line(i:j).eq.'p3') irecty(maxrec(ir),ir)=5
        if(line(i:j).eq.'p4') irecty(maxrec(ir),ir)=6
        if(line(i:j).eq.'p5') irecty(maxrec(ir),ir)=7
        if(line(i:j).eq.'fa') irecty(maxrec(ir),ir)=8
        if(line(i:j).eq.'mo') irecty(maxrec(ir),ir)=9
        if(line(i:j).eq.'st') irecty(maxrec(ir),ir)=10
        if(line(i:j).eq.'x1') irecty(maxrec(ir),ir)=11
        if(line(i:j).eq.'x2') irecty(maxrec(ir),ir)=12
        if(line(i:j).eq.'x3') irecty(maxrec(ir),ir)=13
        if(line(i:j).eq.'x4') irecty(maxrec(ir),ir)=14
        if(line(i:j).eq.'idfa') irecty(maxrec(ir),ir)=15
        if(line(i:j).eq.'idmo') irecty(maxrec(ir),ir)=16
        if(line(i:j).eq.'p') irecty(maxrec(ir),ir)=17
        if(line(i:j).eq.'x') irecty(maxrec(ir),ir)=18
        if(line(i:j).eq.'dez') irecty(maxrec(ir),ir)=19
        if(line(i:j).eq.'c1') irecty(maxrec(ir),ir)=21
        if(line(i:j).eq.'c2') irecty(maxrec(ir),ir)=22
        if(line(i:j).eq.'ty') irecty(maxrec(ir),ir)=23
      endif
      if(irecty(maxrec(ir),ir).eq.-1)then
        write(*,*) 'unknown variable ',line(i:j)
        stop
      endif
      goto 20
 22   continue

           elseif(line(i:j).eq.'CentralityClass')then

      call setCentralityClass(line,i,j,nopen)

           elseif(line(i:j).eq.'CentralityLimit')then

      call utword(line,i,j,0)
      read(line(i:j),*)val
      ival=nint(val)
      call utword(line,i,j,0)
      read(line(i:j),*)val
      zlimit(ival)=val

           elseif(line(i:j).eq.'eos')then

      !Choices of the equation of state (EoS)
      call utword(line,i,j,0)
      if(line(i:j).eq.'x3ff')then       !Phys. Rev. C 89, 064903 (2014)
        ioeos=22
      elseif(line(i:j).eq.'best')then   !Nucl.Phys.A 982 (2019) 183-185
        ioeos=6
      elseif(line(i:j).eq.'chiral')then !Yuirii’s one
        ioeos=30
      elseif(line(i:j).eq.'cem')then    !Phys. Rev. D 97, 114030 (2018)
        ioeos=35
      elseif(line(i:j).eq.'pnjl')then   !The „local one”
        ioeos=7
      elseif(line(i:j).eq.'off')then    !No EoS loaded
        ioeos=0
      endif

           elseif(line(i:j).eq.'hydro')then

      call utword(line,i,j,0)
      if(line(i:j).eq.'effective')then
        stop'hydro `effective` not used any more'
        ihlle=99
      elseif(line(i:j).eq.'x3ff')then
        stop'Use hydro hlle'
      elseif(line(i:j).eq.'hlle')then
        ihlle=1
      elseif(line(i:j).eq.'off')then
        ihlle=0
        ihyskip=0
      elseif(line(i:j).eq.'skip')then
        ihlle=1
        ihyskip=1
      else
        stop'Invalid hydro option'
      endif
      if(iPFE.eq.1.and.line(i:j).ne.'off')then
        write(ifmt,'(a)')'hydro turned off (replaced by PFE)' 
        ihlle=0
        ihyskip=0
      endif

           elseif(line(i:j).eq.'ftime')then

      call utword(line,i,j,0)
      if(line(i:j).eq.'on' ) continue
      if(line(i:j).eq.'off') taustr=0

           elseif(line(i:j).eq.'hacas')then

      call utword(line,i,j,0)
        if(line(i:j).eq.'full')then
      ihacas=1
      iuelast=0
      iuskip=0
      iuchaskip=0
        elseif(line(i:j).eq.'elastic')then
      ihacas=1
      iuelast=1
      iuskip=0
      iuchaskip=0
        elseif(line(i:j).eq.'skip')then
      ihacas=1
      iuelast=0
      iuskip=1
      iuchaskip=0
        elseif(line(i:j).eq.'chaskip')then
      ihacas=1
      iuelast=0
      iuskip=0
      iuchaskip=1
        elseif(line(i:j).eq.'off')then
      ihacas=0
      iuelast=0
      iuskip=0
      iuchaskip=0
        elseif(line(i:i).eq.'u')then
      ihacas=1
      iuelast=0
      iuskip=0
        if(line(i:j).eq.'u0')then !same as full
      iuchaskip=0
        elseif(line(i:j).eq.'u1')then
      iuchaskip=1
        elseif(line(i:j).eq.'u2')then
      iuchaskip=2
        elseif(line(i:j).eq.'u3')then
      iuchaskip=3
        elseif(line(i:j).eq.'u4')then
      iuchaskip=4
        else
          stop'\n\n STOP: invalid hacas u option \n\n   '
        endif
        else
          stop'\n\n STOP: invalid hacas option \n\n   '
        endif
        !write(ifmt,'(a,i1)')'hacas choice u',iuchaskip

           elseif(line(i:j).eq.'core')then

      call utword(line,i,j,0)
        if(line(i:j).eq.'off')then
      iorsdf=0
      icocore=0
      ioclude=3
        elseif(line(i:j).eq.'full')then
      iorsdf=5
      icocore=1
      ioclude=4
        elseif(line(i:j).eq.'rope')then
      iorsdf=6
        elseif(line(i:j).eq.'central')then
      iorsdf=5
      icocore=2
      ioclude=4
        elseif(line(i:j).eq.'PFE'.or.line(i:j).eq.'pfe')then !Parameterized Fluid Expansion
      iPFE=1 
      iorsdf=3
      icocore=0
      ioclude=3
      tauzer=1.
      taustr=1.
      ihlle=0
      ihyskip=0
        else
          stop'\n\n STOP: invalid core option \n\n   '
        endif

           elseif(line(i:j).eq.'satpom')then

      call utword(line,i,j,0)
      if(line(i:j).eq.'off')then
        factjpsi=0.0            !factor for jpsi production in saturated contribution
        factsat=0.0             !factor for saturated contribution
          
      elseif(line(i:j).eq.'full')then
        factjpsi=1.0            !factor for jpsi production in saturated contribution
        factsat=1.0             !factor for saturated contribution
      else
        stop'STOP: invalid satpom option   '
      endif

        
           elseif(line(i:j).eq.'idchoice')then

      call utword(line,i,j,0)
      if(line(i:j).eq.'nxs')then
        ioidch=1
      elseif(line(i:j).eq.'pdg')then
        ioidch=2
      else
        stop'invalid idchoice.     '
      endif

           elseif(line(i:j).eq.'make')then

      call utword(line,i,j,0)
      if(line(i:j).eq.'icotable')icotabm=1

           elseif(line(i:j).eq.'read')then

      call utword(line,i,j,0)
      if(line(i:j).eq.'ico'.or.line(i:j).eq.'IcoTable')then
        icotabr=1
        call utword(line,i,j,0)
        if(line(i:j).eq.'-nsmooth')then
          call utword(line,i,j,0)
          read(line(i:j),*)val
          nsmooth=nint(val)
          call utword(line,i,j,0)
        endif
        if(line(i:j).eq.'$EPO')then
          call utword(line,i,j,0)
          nfnio=nfnnx+j-i+1
          fnio(1:nfnio)=fnnx(1:nfnnx)//line(i:j)
        elseif(line(i:j).eq.'-')then
          continue
        else
          stop'\n\n in aread/read (1) \n\n'
        endif
        if(jcentrality.lt.0)then
         if(nopen.ne.-1)then   !only first read
          fmt4(1:4)='(i )'
          iin=1+log10(float(icentrality))
          write(fmt4(3:3),'(i1)')iin
          jjn=nfnio
          nfnio=jjn+iin
          write(fnio(jjn-3:jjn-4+iin),fmt4)icentrality
          fnio(jjn-3+iin:jjn+iin)='.ico'
          write(ifmt,'(a/2x,a)')'IcoTable set to ',fnio(1:jjn+iin)
         endif
        endif
      elseif(line(i:j).eq.'out'.or.line(i:j).eq.'HydroTable')then
        ireadhyt=ireadhyt+1
        if(ireadhyt.gt.mxho)stop'\n\n in aread/read (2) \n\n'
        call utword(line,i,j,0)
        if(line(i:j).eq.'$EPO')then
          call utword(line,i,j,0)
          nfnho(ireadhyt)=nfnnx+j-i+1
          fnho(ireadhyt)(1:nfnho(ireadhyt))=fnnx(1:nfnnx)//line(i:j)
        else
          stop'\n\n in aread/read (3) \n\n'
        endif
        if(jcentrality.lt.0)then
         if(nopen.ne.-1)then   !only first read
          fmt4(1:4)='(i )'
          iin=1+log10(float(icentrality))
          write(fmt4(3:3),'(i1)')iin
          jjn=nfnho(ireadhyt)
          nfnho(ireadhyt)=jjn+iin
          write(fnho(ireadhyt)(jjn-3:jjn-4+iin),fmt4)icentrality
          fnho(ireadhyt)(jjn-3+iin:jjn+iin)='.out'
          write(ifmt,'(a/2x,a)')'HydroTable set to '
     .   ,fnho(ireadhyt)(1:jjn+iin)
         endif
        endif
      endif

           elseif(line(i:j).eq.'output')then

      call utword(line,i,j,0)
      if(line(i:j).eq.'full' )      istore=-1
      if(line(i:j).eq.'epos' )      istore=1
      if(line(i:j).eq.'osc1997a' )  istore=2
      if(line(i:j).eq.'osc1999a' )  istore=3
      if(line(i:j).eq.'pdg' )       istore=4
      if(line(i:j).eq.'ustore' )    istore=5
      if(line(i:j).eq.'hepmc' )     istore=6
      if(line(i:j).eq.'lhef' )      istore=7

           elseif(line(i:j).eq.'model')then

      call utword(line,i,j,0)
      if(line(i:j).eq.'epos')then
        model=1
      else
        nij=j-i+1
        if(nij.gt.20)stop'cmodel too small'
        cmodel(1:nij)=line(i:j)
        cmodel(nij+1:nij+1)=' '
        call NumberModel(cmodel,model)
      endif
      if(abs(iappl).ne.1.and.iappl.ne.3.and.model.ne.1
     &.and..not.(model.eq.4.and.iappl.eq.7))
     &call utstop('Application not possible with this model&')

           elseif(line(i:j).eq.'trigger'.or.line(i:j).eq.'trg')then

      call utword(line,i,j,0)
      nsingle=0
      if(line(i:j).eq.'-single')then
        nsingle=1
        call utword(line,i,j,1)
      endif
      ntc=1
      if(line(i:j).eq.'or'.or.line(i:j).eq.'-or'
     ..or.line(i:j).eq.'contr')then
        imo=2
        if(line(i:j).eq.'-or')imo=-2
        if(line(i:j).eq.'contr')imo=3
        call utword(line,i,j,0)
        read(line(i:j),*)ztc
        ntc=nint(ztc)
        call utword(line,i,j,1)
      endif
      do n=1,ntc
        if(n.ne.1.and.imo.ne.-2)call utword(line,i,j,0)
        call utword(line,i,j,0)
        if(nsingle.ne.1)then
          call utword(line,i,j,0)
        endif
      enddo

           elseif(line(i:j).eq.'noerrorbut')then

      call utword(line,i,j,0)

           elseif(line(i:j).eq.'b')then

      nbarray=nbarray+1
      call utword(line,i,j,0)
      read(line(i:j),*)val
      barray(nbarray)=val

           elseif(line(i:j).eq.'message')then

      call utword(line,i,j,0)
      if(nopen.eq.-1)then      !only write in second read
      if(.not.(igrTree.gt.0))then !not reading root
      write(ifmt,'(a,$)')line(i:j)
      endif
      endif

           elseif(line(i:j).eq.'endmessage')then

      if(nopen.eq.-1)then      !only write in second read
      write(ifmt,'(a)')' '
      endif

           elseif(line(i:j).eq.'write'.or.line(i:j).eq.'writex'
     .   .or.line(i:j).eq.'writexx'.or.line(i:j).eq.'writexxx')then

      if(line(i:j).eq.'write')ii=1
      if(line(i:j).eq.'writex')ii=2
      if(line(i:j).eq.'writexx')ii=3
      if(line(i:j).eq.'writexxx')ii=4
      call utword(line,i,j,0)

      if(line(i:j).eq.'HydroTable')then

        igethyt=2

      elseif(line(i:j).eq.'IcoTable')then

        icotabm=1

      elseif(line(i:j).eq.'IcoTableAna')then

        icotabm=2

      else

      idol=0
      if(line(i:j).eq.'$')then
       idol=1
       call utword(line,i,j,0)
      endif
      divi=1.
      if(line(i:j).eq.'-divisor')then
       call utword(line,i,j,0)
       read(line(i:j),*)divi 
       call utword(line,i,j,0)
      endif
      yieldx=yield/divi
      !write: only write in second read; writex: only write in first read
      if(ii.eq.1.and.nopen.eq.-1.or.ii.eq.2.and.nopen.ne.-1)then
       call dollartext('$hdtext1;',line,i,j)
       call dollartext('$hdtext2;',line,i,j)
       call dollartext('$hdtext3;',line,i,j)
       call dollartext('$hdtext4;',line,i,j)
       call dollartext('$hdtext5;',line,i,j)
       call dollartext('$hdtext6;',line,i,j)
       call dollartext('$hdtext7;',line,i,j)
       call dollartext('$hdtext8;',line,i,j)
       call dollar('$iscreen;', float(iscreen),line,i,j)
       call dollar('$iq2sat;', float(iq2sat),line,i,j)
       if(j-4.ge.i)then
        do l=i,j-4
         if(line(l:l+4).eq.'$engy')then
           line(l:l+5)='      '
           iengy=alog10(engy)
           if(iengy.eq.0)write(line(l+3:l+3),'(i1)')int(engy)
           if(iengy.eq.1)write(line(l+3:l+4),'(i2)')int(engy)
           if(iengy.eq.2)write(line(l+2:l+4),'(i3)')int(engy)
           if(iengy.eq.3)write(line(l+1:l+4),'(i4)')int(engy)
           if(iengy.eq.4)write(line(l+0:l+4),'(i5)')int(engy)
           if(iengy.eq.5)write(line(l+0:l+5),'(i6)')int(engy)
           if(iengy.eq.6)write(line(l+0:l+6),'(i7)')int(engy)
         endif
        enddo
       endif
       if(j-4.ge.i)then
        do l=i,j-4
         if(line(l:l+4).eq.'$reac')then
           line(l:l)=' '
           if(maproj.eq.  1)line(l+1:l+2)=' p'
           if(maproj.eq.  2)line(l+1:l+2)=' d'
           if(maproj.eq.197)line(l+1:l+2)='Au'
           if(maproj.eq.208)line(l+1:l+2)='Pb'
           if(matarg.eq.  1)line(l+3:l+4)='p '
           if(matarg.eq. 12)line(l+3:l+4)='C '
           if(matarg.eq.197)line(l+3:l+4)='Au'
           if(matarg.eq.208)line(l+3:l+4)='Pb'
         endif
        enddo
       endif
       if(j-6.ge.i)then
        do l=i,j-6
         if(line(l:l+6).eq.'$iversn')then
           write(line(l:l+6),'(f5.3,2x)')
     *     iversn/1000.+0.0001
         endif
        enddo
       endif
       if(j-6.ge.i)then
        do l=i,j-6
         if(line(l:l+6).eq.'$iverso')then
           write(line(l:l+6),'(i5,2x)')
     *     iverso
         endif
        enddo
       endif
       if(j-8.ge.i)then
        do l=i,j-8
         if(line(l:l+8).eq.'$xx1yield')then
          write(line(l:l+8),'(f8.1,1x)')yieldx
         endif
         if(line(l:l+8).eq.'$xxxyield')then
          if(yieldx.lt.0.001)then
           write(line(l:l+8),'(f8.6,1x)')yieldx
          elseif(yieldx.lt.0.01)then
           write(line(l:l+8),'(1x,f7.5,1x)')yieldx
          elseif(yieldx.lt.0.1)then
           write(line(l:l+8),'(1x,f6.4,2x)')yieldx
          elseif(yieldx.lt.1.0)then
           write(line(l:l+8),'(1x,f5.3,3x)')yieldx
          elseif(yieldx.lt.10.)then
           write(line(l:l+8),'(1x,f6.3,2x)')yieldx
          elseif(yieldx.lt.100.)then
           write(line(l:l+8),'(1x,f7.3,1x)')yieldx
          else
           write(line(l:l+8),'(f8.1,1x)')yieldx
          endif
         endif
        enddo
       endif
       if(j-7.ge.i)then
        do l=i,j-7
         if(line(l:l+7).eq.'$xxyield')then
          if(yieldx.lt.1.0)then
           write(line(l:l+7),'(f7.5,1x)')yieldx
          elseif(yieldx.lt.10.)then
           write(line(l:l+7),'(f7.4,1x)')yieldx
          elseif(yieldx.lt.100.)then
           write(line(l:l+7),'(f7.3,1x)')yieldx
          else
           write(line(l:l+7),'(f7.1,1x)')yieldx
          endif
         endif
        enddo
       endif
       if(j-6.ge.i)then
        do l=i,j-6
         if(line(l:l+6).eq.'$xyield')then
          if(yieldx.lt.1.0)then
           write(line(l:l+6),'(f6.4,1x)')yieldx
          elseif(yieldx.lt.10.)then
           write(line(l:l+6),'(f6.3,1x)')yieldx
          elseif(yieldx.lt.100.)then
           write(line(l:l+6),'(f6.2,1x)')yieldx
          else
           write(line(l:l+6),'(f6.1,1x)')yieldx
          endif
         endif 
        enddo
       endif
       if(j-5.ge.i)then
        do l=i,j-5
         if(line(l:l+5).eq.'$yield')then
          if(yieldx.lt.1.0)then
           write(line(l:l+5),'(f5.3,1x)')yieldx
          elseif(yieldx.lt.100.)then
           write(line(l:l+5),'(f5.2,1x)')yieldx
          elseif(yieldx.lt.1000.)then
           write(line(l:l+5),'(f5.1,1x)')yieldx
          elseif(yieldx.lt.10000.)then
           write(line(l:l+5),'(f6.1)')yieldx
          else
           write(line(l:l+5),'(i6)')nint(yieldx)
          endif
         endif
        enddo
       endif
       if(j-8.ge.i)then
        do l=i,j-8
         if(line(l:l+8).eq.'$xxxaverg')then
          if(averg.lt.-1e29)stop'\n\n error averg\n\n'
          if(averg.lt.0.001)then
           write(line(l:l+8),'(f8.7,1x)')averg
          elseif(averg.lt.0.01)then
           write(line(l:l+8),'(f8.6,1x)')averg
          elseif(averg.lt.0.1)then
           write(line(l:l+8),'(f8.5,1x)')averg
          elseif(averg.lt.1.)then
           write(line(l:l+8),'(f8.4,1x)')averg
          elseif(averg.lt.10.)then
           write(line(l:l+8),'(f8.3)')averg
          else
           write(line(l:l+5),'(i6)')nint(averg)
          endif
         endif
        enddo
       endif
       if(j-5.ge.i)then
        do l=i,j-5
         if(line(l:l+5).eq.'$averg')then
          if(averg.lt.-1e29)stop'\n\n error averg\n\n'
          if(averg.lt.1.0)then
           write(line(l:l+5),'(f5.3,1x)')averg
          elseif(averg.lt.100.)then
           write(line(l:l+5),'(f5.2,1x)')averg
          elseif(averg.lt.1000.)then
           write(line(l:l+5),'(f5.1,1x)')averg
          elseif(averg.lt.10000.)then
           write(line(l:l+5),'(f6.1)')averg
          else
           write(line(l:l+5),'(i6)')nint(averg)
          endif
         endif
        enddo
        do l=i,j-5
         if(line(l:l+5).eq.'$sigma')then
          if(sigma.lt.1.0)then
           write(line(l:l+5),'(f5.3,1x)')sigma
          elseif(sigma.lt.100.)then
           write(line(l:l+5),'(f5.2,1x)')sigma
          elseif(sigma.lt.1000.)then
           write(line(l:l+5),'(f5.1,1x)')sigma
          elseif(sigma.lt.10000.)then
           write(line(l:l+5),'(f6.1)')sigma
          else
           write(line(l:l+5),'(i6)')nint(sigma)
          endif
         endif
        enddo
       endif
       if(idol.eq.0)then
        write(ifhi,'(a)')line(i:j)
       else
        write(ifhi,'(a,a,$)')line(i:j),' '
       endif
      elseif(ii.eq.3.and.nopen.ne.-1)then !writexx: only write in first read
       nwritexx=nwritexx+1
       if(nwritexx.gt.20)stop'\n\n ERROR 20062010\n\n'
       twritexx(nwritexx)=line(i:j)
      elseif(ii.eq.4.and.nopen.ne.-1)then !writexxx: only write in first read
       nwritexxx=nwritexxx+1
       if(nwritexxx.gt.50)stop'\n\n ERROR 14122013\n\n'
       twritexxx(nwritexxx)=line(i:j)
      endif

      endif

           elseif(line(i:j).eq.'nozero')then

      nozero=1

           elseif(line(i:j).eq.'ibmin')then

      call utword(line,i,j,0)
      read(line(i:j),*)val
      ibmin=nint(val)

           elseif(line(i:j).eq.'ibmax')then

      call utword(line,i,j,0)
      read(line(i:j),*)val
      ibmax=nint(val)

            elseif(line(i:j).eq.'openxhisto')then

      call utword(line,i,j,0)
      read(line(i:j),*)ihif
      if(nopen.eq.-1)then !second run
      ihifcount(ihif)=ihifcount(ihif)+1
      if(ifhi.eq.35+ihif)stop'\n\n 23012012 \n\n'
      ifhiSave=ifhi
      ifhi=35+ihif
      ifhix(ihif)=ifhi
      if(ihifcount(ihif).eq.1)then
       call getNames(3)
       cfmt='(i )'
       ndig=1
       if(nout.gt.0)ndig=1+log10(float(nout))
       write(cfmt(3:3),'(i1)')ndig
       cbasout(iba+1:iba+2)='xx'
       write(cbasout(iba+3:iba+3),'(i1)')ihif
       cbasout(iba+4:iba+4)='-'
       write(cbasout(iba+5:iba+4+ndig),cfmt)nout
       cbasout(iba+4+ndig+1:iba+4+ndig+6)='.histo'
       open(unit=ifhi,file=cbasout(1:iba+4+ndig+6),status='unknown')
      endif
      endif !second run

            elseif(line(i:j).eq.'closexhisto')then

      if(nopen.eq.-1)then !second run
      if(ifhiSave.eq.0)stop'\n\n 23012012b \n\n'
      ifhi=ifhiSave
      endif !second run


            elseif(line(i:j).eq.'writearray'.or.line(i:j).eq.'wa'
     $     .or.line(i:j).eq.'writehisto')then

      if(nopen.eq.-1)then !second run
       ih=0
       isinglebin=0
       if(line(i:j).eq.'writearray'.or.line(i:j).eq.'wa') ih=1
       call utword(line,i,j,0)
       if(line(i:j).eq.'s')then
        call utword(line,i,j,0)
        linex=line
        ix=i
        jx=j
        call utword(line,i,j,0)
        if(linex(ix:jx).eq.'inicon')stop'error 060307'
       else
        ioint=0
        iocontr=0
        if(line(i:j).eq.'int')then
         ioint=1
         call utword(line,i,j,0)
        endif
        biwi=1
        bishi=0
        bilog3=0
        if(line(i:j).eq.'-bilog')then
         call utword(line,i,j,0)
         read(line(i:j),*)bilog1
         call utword(line,i,j,0)
         read(line(i:j),*)bilog2
         call utword(line,i,j,0)
         read(line(i:j),*)bilog3
         call utword(line,i,j,0)
        endif
        if(line(i:j).eq.'-biwi')then
         call utword(line,i,j,0)
         read(line(i:j),*)val
         biwi=val
         bishi=-0.5
         call utword(line,i,j,0)
        endif
        if(line(i:j).eq.'-bishi')then
         call utword(line,i,j,0)
         read(line(i:j),*)val
         bishi=val
         call utword(line,i,j,0)
        endif
        if(line(i:j).eq.'contr')then
         iocontr=1
         call utword(line,i,j,0)
         read(line(i:j),*)val
         ncontr=nint(val)
         call utword(line,i,j,0)
        endif
        if(line(i:j).eq.'singlebin')then
         call utword(line,i,j,0)
         read(line(i:j),*)val
         bwidth=val
         call utword(line,i,j,0)
         isinglebin=1
         if(nrbins.ne.1)stop'\n\n single bin expected\n\n'
        endif
        read(line(i:j),*)val
        nco=nint(val)
        if(isinglebin.eq.1.and.nco.ne.2)
     .   stop'\n\n nco.ne.2 not possible for singlebin\n\n'
        if(ih.eq.1)write(ifhi,'(a,i3)')'array',nco
        if(ioint.eq.0)then
         sum=0
         averg=0
         do k=1,nrbins
          if(iocontr.eq.0.and.ionoerr.eq.0)then
            ar3=ar(k,3)
            ar4=ar(k,4)
          elseif(ionoerr.eq.1)then
            ar3=ar(k,3)
            ar4=ar(k,4)
          elseif(ionoerr.eq.2)then
            ar3=ar(k,3)
            ar4=ar(k,4)
            ar5=ar(k,5)
          else
            ar3=ary(k,ncontr)
            ar4=ardy(k,ncontr)
          endif
          iok=1
          if(k.lt.ibmin.or.k.gt.ibmax)iok=0
          if(nint(bilog3).gt.0)then
            del=log(bilog2/bilog1)/nint(bilog3)
            ar1=log(bilog1)+(ar(k,1)-0.5)*del
            ar1=exp(ar1)
          else
            ar1=(ar(k,1)+bishi)*biwi
          endif
          sum=sum+ar3
          averg=averg+ar1*ar3
          if(nco.eq.2)then
            if(nozero.eq.1.and.ar3.eq.0.)iok=0
            if(iok.eq.1)then
              if(isinglebin.eq.0)then
                write(ifhi,'(3e12.4)')ar1,ar3
              else
                write(ifhi,'(3e12.4)')ar1-bwidth/2,ar3
                write(ifhi,'(3e12.4)')ar1+bwidth/2,ar3
              endif
            endif
          elseif(nco.eq.3)then
            if(nozero.eq.1.and.ar3.eq.0..and.ar4.eq.0.)iok=0
            if(iok.eq.1)write(ifhi,'(3e12.4)')ar1,ar3,ar4
          elseif(nco.eq.4)then
          if(nozero.eq.1.and.ar3.eq.0..and.ar4.eq.0..and.ar5.eq.0.)iok=0
            if(iok.eq.1)write(ifhi,'(4e12.4)')ar1,ar3,ar4,ar5
          endif
         enddo
         if(sum.gt.0.)averg=averg/sum
        else
         if(isinglebin.eq.1)stop'\n\n not possible for singlebin\n\n'
         sum=0.
         sum2=0.
         sum3=0.
         err2=0.
         do k=1,nrbins
          if(iocontr.eq.0.and.ionoerr.eq.0)then
            ar3=ar(k,3)
            ar4=ar(k,4)
          elseif(ionoerr.eq.1)then
            ar3=ar(k,3)
            ar4=ar(k,4)
          elseif(ionoerr.eq.2)then
            ar3=ar(k,3)
            ar4=ar(k,4)
            ar5=ar(k,5)
          else
            ar3=ary(k,ncontr)
            ar4=ardy(k,ncontr)
          endif
          sum=sum+ar3*(ar(2,1)-ar(1,1))
          if(nco.eq.2)write(ifhi,'(3e12.4)')ar1,sum
          if(ionoerr.eq.0)then
            err2=err2+(ar4*(ar(2,1)-ar(1,1)))**2
            if(nco.eq.3)write(ifhi,'(3e12.4)')ar1,sum,sqrt(err2)
          elseif(ionoerr.eq.1)then
            sum2=sum2+(ar4*(ar(2,1)-ar(1,1)))
            if(nco.eq.3)write(ifhi,'(3e12.4)')ar1,sum,sum2
          elseif(ionoerr.eq.2)then
            sum2=sum2+(ar4*(ar(2,1)-ar(1,1)))
            sum3=sum3+(ar5*(ar(2,1)-ar(1,1)))
            if(nco.eq.3)write(ifhi,'(3e12.4)')ar1,sum,sum2
            if(nco.eq.4)write(ifhi,'(3e12.4)')ar1,sum,sum2,sum3
          endif
         enddo
        endif
        if(ih.eq.1)write(ifhi,'(a)')'endarray'
       endif
      else !nopen .ge. 0 -- first run
        call utword(line,i,j,0)
        if(line(i:j).eq.'s')then
          call utword(line,i,j,0)
          call utword(line,i,j,0)
        else
         if(line(i:j).eq.'int')then
          call utword(line,i,j,0)
         endif
        if(line(i:j).eq.'-bilog')then
         call utword(line,i,j,0)
         call utword(line,i,j,0)
         call utword(line,i,j,0)
         call utword(line,i,j,0)
        endif
         if(line(i:j).eq.'-biwi')then
          call utword(line,i,j,0)
          call utword(line,i,j,0)
         endif
         if(line(i:j).eq.'contr')then
          call utword(line,i,j,0)
          call utword(line,i,j,0)
         endif
         if(line(i:j).eq.'singlebin')then
          call utword(line,i,j,0)
          call utword(line,i,j,0)
         endif
         if(line(i:j).eq.'-bishi')then
          call utword(line,i,j,0)
          call utword(line,i,j,0)
         endif
        endif
      endif
      nozero=0
      ibmin=1
      ibmax=1e8

           elseif(line(i:j).eq.'END'.or.line(i:j).eq.'XXX')then

      continue

           else

      write(ifmt,'(72a1)')('-',k=1,72)
      write(ifmt,'(a)')line(1:80)
      write(ifmt,*)'i j :  ',i,j
      write(ifmt,'(a,a,a)')'ERROR: command "',line(i:j),'" not found'
      write(ifmt,'(72a1)')('-',k=1,72)
      j=1000
      stop

           endif

      if(itit.eq.1)then
        call atitle(6)
        itit=0
      endif

      i=j+1
      goto 1

      end

c-----------------------------------------------------------------------
      subroutine aseed(modus)
c-----------------------------------------------------------------------

#include "aaa.h"
      double precision seedf
      call utpri('aseed ',ish,ishini,3)

      call ranfgt(seedf)
      if(iwseed.eq.1)then
        if(nrevt.eq.0)then
          write(ifmt,'(a,i10,a,d27.16)')
     *     'seedj:',nint(seedj),'  seedf:',seedf
        elseif(mod(nrevt,modsho).eq.0)then
          if(modus.eq.1)
     *   write(ifmt,'(a,i10,5x,a,i10,a,d27.16)')
     *             'nrevt:',nrevt,'seedj:',nint(seedj),'  seedf:',seedf
          if(modus.eq.2)
     *   write(ifmt,'(a,i10,a,d27.16)')
     *         'seed:',nint(seedj),'  seedf:',seedf
        endif
        if(jwseed.eq.1)then
         open(unit=1,file=fnch(1:nfnch-5)//'see',status='unknown')
         write(1,'(a,i10,5x,a,i10,a,d27.16)')
     *           'nrevt:',nrevt,'seedj:',nint(seedj),' seedf:',seedf
         close(1)
        endif
      endif
      seedc=seedf

      call utprix('aseed ',ish,ishini,3)
      return
      end

c-----------------------------------------------------------------------
      subroutine aseedi
c-----------------------------------------------------------------------

#include "aaa.h"
      call utpri('aseedi',ish,ishini,3)

      if(ish.ge.1)write(ifmt,'(a,i10)')'seedi:',nint(seedi)

      call utprix('aseedi',ish,ishini,3)
      return
      end

c$$$c-----------------------------------------------------------------------
c$$$        subroutine aseed(modus)        !Flush ????
c$$$c-----------------------------------------------------------------------
c$$$
c$$$#include "aaa.h"
c$$$      double precision seedf
c$$$      call utpri('aseed',ish,ishini,4)
c$$$
c$$$      call ranfgt(seedf)
c$$$      if(modus.eq.2)then
c$$$        write(ifmt,'(a,d26.15)')'seed:',seedf
c$$$      elseif(modus.eq.1)then
c$$$        if(mod(nrevt,modsho).eq.0)then
c$$$          write(ifmt,100)'nrevt:',nrevt,'seedf:',seedf
c$$$          call flush(ifmt)
c$$$        endif
c$$$      endif
c$$$      seedc=seedf
c$$$
c$$$  100 format(a,i10,10x,a,d26.15)
c$$$      call utprix('aseed',ish,ishini,4)
c$$$      return
c$$$      end
c$$$
c-----------------------------------------------------------------------
      subroutine astati
c-----------------------------------------------------------------------

#include "aaa.h"
      common/geom/rmproj,rmtarg,bmax,bkmx
      common/ghecsquel/anquasiel,iquasiel

      call utpri('astati',ish,ishini,1)
      if(ish.ge.1.and.iappl.eq.1.)then
        if(abs(accept+reject).gt.1.e-5)write(ifch,'(a,f9.5)')
     *' EMS acc.rate:',accept/(accept+reject)
        if(antot.ne.0.)write(ifch,*)'initial soft,hard(%)'
     *                           ,ansf/antot*100.
     *                           ,ansh/antot*100.,' of' ,antot
        if(antotf.ne.0.)write(ifch,*)'final soft,hard(%)'
     *                           ,ansff/antotf*100.
     *                           ,anshf/antotf*100.,' of' ,antotf
        if(antotre.ne.0.)write(ifch,*)
     *                  'droplet,string(+d),reson(+d), (had)(%) '
        if(antotre.ne.0.)write(ifch,*)'     '
     *                           ,andropl/antotre*100.
     *                           ,anstrg0/antotre*100.
     *                           ,'(',anstrg1/antotre*100.,') '
     *                           ,anreso0/antotre*100.
     *                           ,'(',anreso1/antotre*100.,') '
        if(antotre.ne.0.)write(ifch,*)'     '
     *             ,' (',anghadr/antotre*100.,')',' of' ,antotre
       if(pp4ini.gt.0.)write(ifch,*)'Energy loss',(pp4ini-pp4max)/pp4ini
      write(ifch,*)'ine cross section:',sigineex
      write(ifch,*)'diffr cross section:',sigdifex
      write(ifch,*)'SD cross section:',sigsdex
c      if(model.eq.3)write(ifch,*)'quasi-elastic cross section:'
c     &,anquasiel/float(ntevt)*a*10
         endif

c$$$      call testconex(3)

      if(iprmpt.le.0)goto1000

      write(ifch,'(77a1)')('-',i=1,77)
      write(ifch,'(a)')'statistics'
      write(ifch,'(77a1)')('-',i=1,77)
      write(ifch,'(a,i6)')'nr of messages:',imsg
      write(ifch,'(a,i8)')'maximum nptl:',nptlu
      write(ifch,'(77a1)')('-',i=1,77)

1000  continue
      call utprix('astati',ish,ishini,1)

      call clop(3)
      return
      end


c-----------------------------------------------------------------------
!-----------------------------------------------------------------------
!> @brief
!> display a header at epos start
!
!> Displays "EPOS4." followed by the EPOS version from the file VERSION.txt
!> Displays the content of the file README.md until a horizontal rule '----'
! 
!> @param[in] ii writeUnit
!-----------------------------------------------------------------------
      subroutine atitle(ii)
c-----------------------------------------------------------------------
      implicit none
      integer, intent(in) :: ii
      integer :: ios, jj, writeUnit
      integer, parameter :: readUnit = 99
      integer, parameter :: maxLineSize = 80
      character(len=80) :: lines(30)
      character(len=15) :: version
      character(len=22) :: title
      character(len=20) :: formatString
      character(len=255) :: eposdir, readmeFile, versionFile
      integer :: n, i, nMax, lenMax, offset

#include "aaa.h"

      jj=ifmt
      if(ii.eq.6)jj=6
      writeUnit = jj

      call get_environment_variable("EPO4", eposdir)

c     Read the file VERSION.txt
      versionFile = trim(eposdir)//'VERSION.txt'
      open(unit=readUnit, file=versionFile, 
     .     iostat=ios, status='old')
      if ( ios /= 0 ) stop "Error opening file VERSION.txt"

      do
         read(unit=readUnit, fmt='(a)', iostat=ios) version
         if (ios /= 0) exit
      end do

      close(readUnit)

c     Read the file README.md
      readmeFile = trim(eposdir)//'README.md'
      open(unit=readUnit, file=readmeFile, 
     .     iostat=ios, status='old')
      if ( ios /= 0 ) stop "Error opening file README.md"

c     n : line number
      n = 0     
c     lenMax : maximal line length
      lenMax = 0
      do 
         n = n + 1
         read(readUnit, '(a)', iostat=ios) lines(n)
         lines(n) = trim(lines(n))
         if (ios /= 0) exit
c        read until a horizontal rule line
         if (lines(n).eq.'----') then
            exit
         endif 
         lenMax = max(lenMax, len(trim(lines(n))))
      end do
c     nMax : number of lines to print
      nMax = n-1    
      close(readUnit)

c     Print the header
c     print a line of maxLineSize characters "#"
      write(formatString, '("(a",i2,")")') maxLineSize
      write(unit=writeUnit, fmt=formatString) 
     .     repeat("#", maxLineSize)
c     print a line of maxLineSize-2 blank characters beginning and ending with "#"
      write(formatString, '("(a, a",i2,", a)")') maxLineSize-2
      write(unit=writeUnit, fmt=formatString)
     .     "#", repeat(" ", maxLineSize-2), "#"

c     print the EPOS version
      write(formatString, '("(a",i2,",a",i2,")")') 
     .         len("EPOS "), len(version)
      write(title,formatString) "EPOS ", version
c     print a centered line with title
c     compute the beginning offset
      offset = (maxLineSize - len(trim(title))) / 2      
      write(formatString, '("(a,a",i2,",a",i2,",a",i2,",a)")') 
     .         offset,
     .         len(trim(title)),
     .         maxLineSize-offset-len(trim(title))-2
      write(unit=writeUnit, fmt=formatString)
     .         "#", 
     .         repeat(" ", offset),
     .         title,
     .         repeat(" ", maxLineSize-offset-len(trim(title))-2),
     .         "#"

c     print the centered lines from README informations
c     compute the beginning offset
      offset = (maxLineSize - lenMax) / 2
      do i=2, nMax
         if (len(trim(lines(i))).gt.0) then
            write (formatString, '("(a, a",i2,",a",i2,",a",i2,",a)")') 
     .         offset,
     .         len(trim(lines(i))),
     .         maxLineSize-offset-len(trim(lines(i)))-2
            write(unit=writeUnit, fmt=formatString, advance='no')
     .         "#", 
     .         repeat(" ", offset),
     .         lines(i),
     .         repeat(" ", maxLineSize-offset-len(trim(lines(i)))-2),
     .         "#"
            write(unit=writeUnit,fmt='(a)') ''
         else
c           print a line of maxLineSize-2 blank characters beginning and ending with "#"
            write (formatString, '("(a, a",i2,", a)")') maxLineSize-2
            write(unit=writeUnit, fmt=formatString)
     .           "#", repeat(" ", maxLineSize-2), "#"
         end if
      end do
      write (formatString, '("(a",i2,")")') maxLineSize
      write(unit=writeUnit, fmt=formatString) 
     .     repeat("#", maxLineSize)

      return
      end

c-----------------------------------------------------------------------
      subroutine avehep
c-----------------------------------------------------------------------

#include "aaa.h"

      call utpri('avehep',ish,ishini,4)


      call utprix('avehep',ish,ishini,4)
      end


c-----------------------------------------------------------------------
      subroutine aepos(nin,nfr,icou)
c-----------------------------------------------------------------------

#include "aaa.h"
#include "ems.h"
      double precision eppass,etpass
      common/emnpass/eppass(mamx,4),etpass(mamx,4)
      common/photrans/phoele(4),ebeam,noevt
      common/cninx/ninx
      common/cikoevt/ikoevtxx /cnglevt/nglevtxx/ csegevt/segevtxx
      character ccc*23
      common/cnfrx/nfrx
      save ntry,ntevt0,andropl0,anstrg00,anstrg10,anreso00
      save anreso10,anghadr0,antotre0,anintdiff0,anintsdif0,anintine0
     .    ,anintddif0
      double precision iutime(5),tiu3,tiu4,tidi,tidisu,tidisu0
      data nptly/0/
      save nptly
      call utpri('aepos ',ish,ishini,4)

      call timer(iutime)
      tidisu=0
      tidisu0=0

      ninx=iabs(nin)
      nfrx=nfr

      if(icou.gt.1)goto 3        !error after aepos (outside)

      if(ish.ge.2.and.noevt.eq.0)then
          ccc='start event number     '
          write(ccc(19:23),'(i5)')nrevt+1
          call alist(ccc//'&',0,0)
      endif

c if random sign for projectile, set it here
      if(irdmpr.ne.0.and.maproj.eq.1)then
        idproj=idprojin*(1-2*int(rangen()+0.5d0))
        call emsini(engy,idproj,idtarg,1) !recall emsini to set initial valence quark properly
      endif

c for Air target, set the target nucleus
      if(idtargin.eq.0.and.model.ne.6)then
        call getairmol(latarg,matarg)
        if(ish.ge.2)write(ifch,*)'Air Target, select (Z,A) :'
     &                           ,latarg,matarg
      endif

      if(iappl.ne.5)then ! NOT kinky
        nptl=0
        do i=1,40
          idptl(i)=0
          istptl(i)=-999999
          ityptl(i)=-999999
        enddo
      endif

      if(iappl.eq.4)then !~~~~~~thermal~~~~~~
      nptlpt=0
      ntry=0
  1   ntry=ntry+1
      if(ntry.gt.100)stop'in aepos, to many athermal attempts.    '
      call athermal(iret)
      if(iret.ne.0)goto 1
      if(ish.ge.2)call alist('list before int/decays&',1,nptl)
      nevt=1
      nbdky=nptl
      call bjinta(ier)
      if(ier.ne.0)goto 1
      if(ish.ge.2)call alist('list after int/decays&',1,nptl)
      goto 1000
      else
      nptlpt=abs(maproj)+abs(matarg) !has to be defined here because it is used in utresc
      endif

      if(iappl.eq.9)then !~~~~~~amicro~~~~~~
      nptlpt=0
      ntry=0
  11  ntry=ntry+1
      if(ntry.gt.100)stop'in aepos, to many amicro attempts.    '
      call amicro(iret)
      if(iret.ne.0)goto 11
      if(ish.ge.2)call alist('list before int/decays&',1,nptl)
      nevt=1
      nbdky=nptl
      call bjinta(ier)
      if(ier.ne.0)goto 11
      if(ish.ge.2)call alist('list after int/decays&',1,nptl)
      goto 1000

      else
      nptlpt=abs(maproj)+abs(matarg) !has to be defined here because it is used in utresc
      endif

      ntry=0
      noevt=0
      nr3=0
      if(nin.le.1.and.(nfr.eq.0.or.ireadhyt.ne.0))bimevt=-1
c save statistic at last inelastic event
      ntevt0=ntevt
      andropl0=andropl
      anstrg00=anstrg0
      anstrg10=anstrg1
      anreso00=anreso0
      anreso10=anreso1
      anghadr0=anghadr
      antotre0=antotre
      anintdiff0=anintdiff
      anintsdif0=anintsdif
      anintddif0=anintddif
      anintine0=anintine
 3    continue !set value back to last inelastic event
      nr3=nr3+1
      if(nr3.gt.1)then
        tiu3=iutime(3)
        tiu4=iutime(4)
        call timer(iutime)
        tidi=iutime(3)-tiu3 + (iutime(4)-tiu4)*1d-3
        tidisu=tidisu+tidi
        tidisu0=tidisu0+tidi
        !print*,iutime(3)-tiu3,iutime(4)-tiu4,tidi
        if(tidisu0.gt.600d0)then
          write(ifmt,'(a,f12.4,2a,$)')
     .    ' redo;  accum time =',tidisu,' seconds;  ',
     .    'nr of attempted evts ='
          write(ifmt,*)nr3
          call clop(3)
          tidisu0=0
        endif
      endif
      ntevt=ntevt0
      andropl=andropl0
      anstrg0=anstrg00
      anstrg1=anstrg10
      anreso0=anreso00
      anreso1=anreso10
      anghadr=anghadr0
      antotre=antotre0
      anintdiff=anintdiff0
      anintsdif=anintsdif0
      anintddif=anintddif0
      anintine=anintine0
 2    if(nfr.eq.0)ntevt=ntevt+1
      iret=0
      ntry=ntry+1
      if(iappl.eq.1.or.iappl.eq.2)naevt=naevt+1
 5    nevt=0
      if(nrevt.eq.0)nptly=nptl
      if(iappl.ne.5)call cleanup
      if(irewch.eq.2)rewind(ifch)
      ntxt80=0
      icccptl=0
      nptl=0
      maxfra=0
      if(iappl.ne.5)then ! NOT kinky
        do i=1,40
          idptl(i)=0
          istptl(i)=0
          ityptl(i)=0
        enddo
      endif
    
      if(nfreeze.gt.1.and.nfr.gt.0.and.jcorona.le.0)goto 77

      minfra=mxptl

      if(iappl.eq.1.or.iappl.eq.2)then !---hadron---geometry---

      if(ntry.lt.1000
     ..and.engy.ge.egymin
     ..or.izmode.ge.2
     ..or.ikolmn.gt.10)then
        !if no inel scattering -> nothing !
        if(model.eq.1)then
          call emsaaa(nfr,ntry,iret)
          if(iret.eq.55)return 
          if(iret.eq.-2)goto 3
          call correctLowNpom(iskip)
          if(iskip.eq.1)then
            if(nin.le.1.and.nfr.eq.0)bimevt=-1
            goto 3
          endif
        else
          call emsaaaModel(model,idtargin,iret)
        endif
        if(iret.lt.0)then       !ncol=0
          if(noebin.lt.0)then
            if(ish.ge.1)then
              write(ifch,*)
              write(ifch,*)'no collision, again ...'
              write(ifch,*)
              write(ifch,*)('#',k=1,20)
              write(ifch,*)
            endif
            noevt=1             !fake DIS : if no collision,
            !return              ! try new Q2 and Y (in ainit)
            stop'ERROR fake DIS needs extra compile mode'
          endif
          goto 2
        elseif(iret.gt.0)then
          goto 3                !error
        endif
      else
        if(ish.ge.2)
     &  write(ifch,*)'Nothing done after ',ntry,' ntry ... continue'
        if(ish.ge.1)
     &  write(ifmt,*)'Nothing done after ',ntry,' ntry ... continue'
        ntevt=ntevt0+100
        iret=0
        nevt=1
        call conre     !define projectile and target (elastic scattering)
        call conwr     !when the MC is suppose to produce something but failed
      endif
      if(iappl.eq.2.or.jcorona.eq.2)then
        nevt=1
        goto 1000
      endif

      elseif(iappl.eq.3.or.iappl.eq.-1) then
        nevt=1
        call bread
        if(ish.ge.2)call alist('list after reading&',1,nptl)
        goto 1000

      elseif(iappl.eq.5) then !------kinky------

         nptl=nptly
         nevt=1
         do i=1,nptl
           istptl(i)=20
         enddo

      elseif(iappl.eq.6)then !--------ee-------

        call timann
        nevt=1

      elseif(iappl.eq.7)then !------decay------

        call conwr
        nevt=1

      elseif(iappl.eq.8)then !------lepton-----

        call psadis(iret)
        if(iret.gt.0)goto5
        nevt=1

      endif !----------------------------------------

      if(nevt.eq.0)stop'\n\n ERROR 12072011 \n\n'

      if(nfreeze.gt.1.and.nfr.eq.0.and.jcorona.le.0)call dumpList
 77   if(nfreeze.gt.1.and.nfr.gt.0.and.jcorona.le.0)call readList(nfr)

      !call etotcheck

      if(ish.ge.2)call alist('list before fragmentation&',1,nptl)
      nptlx=nptl+1
      if(iappl.ne.2.and.iappl.ne.7.and.nevt.eq.1.and.ifrade.ne.0)then
        call gakfra(iret)
        if(iret.gt.0)goto 3
        !call alist3('After gak:&',30) 
        !cbg   call redefineMassPartons
        if(iappl.eq.1)then
          call utrescxx(iret,0) !because of off-shell correction in gakfra
          if(iret.gt.0)goto 3
        endif
        maxfra=nptl   !after fragmentation
        if(ish.ge.2.and.model.eq.1)
     &              call alist('list after fragmentation&',nptlx,nptl)
        !~~~~~~~~~~~~~~~~~~~
        !do n=1,nptl
        !igo=0
        !call idflav(idptl(n),i1,i2,i3,jdu1,jdu2)
        !if(i1.gt.0.and.i2.gt.0.and.i3.gt.0)then
        !amt=sqrt(pptl(5,n)**2+pptl(1,n)**2+pptl(2,n)**2)
        !rap=sign(1.,pptl(3,n))*alog((pptl(4,n)+abs(pptl(3,n)))/amt)
        !if(ityptl(n)/10.eq.4.and.rap.lt.-2.)igo=1
        !if(ityptl(n)/10.eq.5.and.rap.gt. 2.)igo=1
        !if(igo.eq.1)then
        !print*,'TESTgak',n,idptl(n),rap,istptl(n),ityptl(n)
        !endif
        !endif
        !enddo
        !~~~~~~~~~~~~~~~~~~~
        !esu=0
        !do i=1,nptl
        !if(istptl(i).eq.0)esu=esu+pptl(4,i)
        !enddo
        !write(ifmt,'(a,41x,f15.1)')' +++++Eafrag+++++',esu
        if(irescl.eq.1)then
          call utghost(iret)
          if(iret.gt.0)goto 3
        endif
        !nptlx=nptl+1
      endif
      !call etotcheck

      if(ispherio.eq.1.and.irescl.eq.1)then
        call utrsph(iret)
        if(iret.gt.0)goto 3
      endif

      if(iappl.lt.5)call iniCore(1)

      if(iappl.ne.2.and.nevt.eq.1)then
        nbdky=nptl
        call bjinta(ier)
        if(ier.eq.1)goto 3
        call Segments(nfr,ntry,iret)
        if(iret.eq.-2)goto 3
        if(iappl.eq.1.and.irescl.eq.1)then
          call utrescxx(iret,0)  !after bjinta if irescl=1
          !call  utresc(iret,1.02)  !after bjinta if irescl=1
          if(iret.gt.0)goto 3
        endif
        if(ish.ge.1)then
          !if(ish.ge.2.and.ifrade.ne.0)
          !&    call alist('list after int/decays&',1,nptl)
          if(iappl.eq.1)then
            numbar=0
            pp4=0.
            do j=1,nptl
              if(istptl(j).eq.0)then
             if(idptl(j).gt. 1000.and.idptl(j).lt. 10000)numbar=numbar+1
             if(idptl(j).lt.-1000.and.idptl(j).gt.-10000)numbar=numbar-1
             if((((idptl(j).eq.1120.or.idptl(j).eq.1220)
     *           .and.idproj.gt.1000).or.(iabs(idptl(j)).gt.100
     *           .and.idproj.lt.1000)).and.pptl(4,j)
     *           .gt.pp4.and.pptl(3,j).gt.0.)pp4=pptl(4,j)
              endif
            enddo
            pp4max=pp4max+pp4
            pp4ini=pp4ini+pptl(4,1)
            nvio=isign(matarg,idtarg)-numbar
            if(iabs(idproj).gt.1000)then
              nvio=nvio+isign(maproj,idproj)
            elseif(iabs(idproj).eq.17)then
              nvio=nvio+isign(2,idproj)
            elseif(iabs(idproj).eq.18)then
              nvio=nvio+isign(3,idproj)
            elseif(iabs(idproj).eq.19)then
              nvio=nvio+isign(4,idproj)
            endif
            if(ish.ge.3)write (ifch,*)'- Baryon number conservation : '
     &                  ,nvio,' -'

          endif
          if(ish.ge.2.and.ifrade.ne.0)
     *    call alist('list after bjinta&',1,nptl)
        endif
      endif

      if((iappl.eq.1.or.iappl.eq.2).and.nevt.eq.0)then
        if(nin.le.1.and.(nfr.eq.0.or.ireadhyt.ne.0))bimevt=-1
        goto 2
      endif

      if(ifrade.ne.0.and.iappl.eq.2.and.
     $     idproj.eq.1120.and.idtarg.eq.1120)then
       numbar=0
       do j=1,nptl
        if(istptl(j).eq.0)then
         if(idptl(j).gt. 1000.and.idptl(j).lt. 10000)numbar=numbar+1
         if(idptl(j).lt.-1000.and.idptl(j).gt.-10000)numbar=numbar-1
        endif
       enddo
       nvio=maproj+matarg-numbar
       if(nvio.ne.0)then
        call alist('complete list&',1,nptl)
        write(6,'(//10x,a,i3//)')'ERROR: baryon number violation:',nvio
        write(6,'(10x,a//)')
     *        'a complete list has been printed into the check-file'
        stop
       endif
      endif


      !ifirst=0
      !if(nrevt+1.eq.1)ifirst=1
      !if(jpsi.gt.0)then
      !  npjpsi=0
      !  do i=1,jpsi
      !    call jpsifo(npjpsi)
      !    call jpsian(ifirst)
      !  enddo
      !  if(ish.ge.1)call jtauan(0,0)
      !  if(nrevt+1.eq.nevent)call jpsihi
      !endif

      if(ixtau.eq.1)call xtauev(1)

1000  continue
      nglacc=nglacc+nglevt
      if(ninx.eq.1.and.nfr.eq.0)ikoevtxx=ikoevt
      if(ninx.eq.1.and.nfr.eq.0)nglevtxx=nglevt
      if(ninx.eq.1.and.nfr.eq.0)segevtxx=segevt
      call utprix('aepos ',ish,ishini,4)
      call clop(3)
      return
      end

      subroutine etotcheck
#include "aaa.h"
      einit=maproj*engy/2+matarg*engy/2
      esu=0
      do i=1,nptl
        if(mod(istptl(i),10).eq.0)then
          esu=esu+pptl(4,i)
        endif
      enddo
      print*,'####### ETOT #######',esu,einit
      end

c-----------------------------------------------------------------------
      subroutine Segments(nfr,ntry,iret)
c-----------------------------------------------------------------------
#include "aaa.h"
      common/csegevt/segevtxx
      segevt=0
      do n=minfra,maxfra
      if(istptl(n).ge.5.and.istptl(n).le.7
     ..or.istptl(n).eq.0.or.istptl(n).eq.1)then
        !rapx=dezptl(n)
        !if(abs(rapx).le.2.5.)then
        segevt=segevt+1.
        !endif
      endif
      enddo
      !print*,'+++++',minfra,maxfra,segevt

      iret=0
      if(izmode.eq.4)then
        p1=segmin
        p2=segmax
        if(nfr.gt.0)then
          dseg=(p2-p1+1)/2.
          if(p2.gt.100000.and.zclass(3,1).gt.1e-5)then
            dseg=(zclass(2,1)-zclass(1,1))/2
          endif
          p1=max(p1,segevtxx-dseg)
          p2=min(p2,segevtxx+dseg)
        endif
        pp=segevt
        if(pp.lt.p1.or.pp.gt.p2)iret=-2
      endif

      !print*,'+++++',pp,p1,p2,nfr,iret
      if(izmode.eq.4.and.iret.eq.0.and.
     . (ihlle.eq.1.or.ispherio.eq.1.or.nfr.gt.0))
     .  write(ifmt,'(i7,a,f7.2,a,f7.2,a,f7.2,a)')
     .  ntry,' attempts to get segment multiplicity',pp
     .  ,' in',p1,' -',p2

      end

c-----------------------------------------------------------------------
      subroutine cleanup
c-----------------------------------------------------------------------
#include "aaa.h"
      do i=1,min(nptl,mxptl+29)
        do  k=1,5
          pptl(k,i)=0
        enddo
        iorptl(i)  =0
        jorptl(i)  =0
        idptl(i)   =0
        istptl(i)  =0
        tivptl(1,i)=0
        tivptl(2,i)=0
        ifrptl(1,i)=0
        ifrptl(2,i)=0
        ityptl(i)  =0
        iaaptl(i)  =0
        radptl(i)  =0
        dezptl(i)  =0
        itsptl(i)  =0
        rinptl(i)  =-9999
        do  k=1,4
          xorptl(k,i)=0
          ibptl(k,i) =0
        enddo
      enddo
      end

c-----------------------------------------------------------------------
      subroutine emsaaa(nfr,ntry,iret)
c-----------------------------------------------------------------------
c  basic EMS routine to determine Pomeron configuration
c-----------------------------------------------------------------------

#include "aaa.h"
#include "ems.h"
      common/nucl3/phi,bimp
      common/col3/ncol,kolpt,ncoli
      common/cninx/ninx  /ciotype/iotype
      logical glauber,only_geom
      common/cikoevt/ikoevtxx /cnglevt/nglevtxx
      common/czptav/zptav
      common/cnparticip/jproj(2,mamx),jtarg(2,mamx),efluct(6,mamx)
c      dimension ZpartPro(0:2),ZpartTar(0:2)

      call utpri('emsaaa',ish,ishini,4)
      if(ish.ge.3)call alist('Determine Pomeron Configuration&',0,0)

      iret=0
      if(iappl.eq.2.and.ionudi.gt.0)
     .stop'\n\n ERROR 22092012 (this option gives wrong results)\n\n'
      ! "geometry" only works for  ionudi = 0, otherwise use "hadron"
      only_geom=iappl.eq.2.or.jcorona.eq.2
      glauber=ionudi.eq.0.and.(maproj.ne.1.or.matarg.ne.1)
      if(glauber.and.only_geom)kexit=1

      nptl=0
      call conaa(nfr,iret)
      if(iret.gt.0.and..not.(glauber.and.only_geom))goto 1001
      call conre
      call conwr
      !call GfunParK(iret)
      !call setzzpos
      if(iappl.eq.2.and.ixgeometry.eq.1)call xGeometry(1)
             !zptav=0
             !do k=1,koll
             !call Zparticip(k,ZpartPro,ZpartTar)
             !zptav=zptav+ZpartPro(1)+ZpartTar(1)
             !enddo
             !npartic=0
             !do i=1,mamx
             !if(jproj(3,i).eq.1)npartic=npartic+1
             !if(jtarg(3,i).eq.1)npartic=npartic+1
             !enddo
             !zptpa=float(npartic)/(maproj+matarg)
             !zptgl=float(ng1evt)/(maproj+matarg)
             !zptav=zptav/maproj/matarg
             !print*,'ZPT',percentile(),zptav,zptpa,zptgl
      if(izmode.eq.3)then
        i1=nglmin
        i2=nglmax
        if(nfr.gt.0)then
        i1=nglevtxx
        i2=nglevtxx
        endif
        if(nglevt.lt.i1.or.nglevt.gt.i2)then
          iret=-2
          goto1000
        endif
      endif
      if(only_geom)then
        if(glauber.and.nglevt.eq.0)goto 1001
        bimevt=bimp
        goto 1000
      endif

      kk=0
      if(izmode.eq.1)then
       do k=1,100
        if(zclass(3,k).gt.0.and.
     .   bimp.ge.zclass(1,k).and.bimp.le.zclass(2,k))then
         kk=k
        endif
       enddo
      endif
      if(kk.gt.0)then
       del=(bimp-zclass(1,kk))/(zclass(2,kk)-zclass(1,kk))
       ctrx=(kk-1+del)*5
      else
       ctrx=-2.5
      endif
      ctrevt=kk*5-2.5

      if(iret.gt.0)goto 1000    !error
      if(glauber.and.nglevt.eq.0)goto 1001
      call emsaa(nfr,iret)
      if(iret.gt.0)goto 1000    !error or analysis mode (55)
      call emsaaDeb
      if(ncol.eq.0)goto 1001 !no interaction
      if(iotype.eq.1.and.abs(nint(typevt)).ne.1
     ..and.(maproj.gt.2.or.matarg.gt.2) )then
         bimevt=-1
         goto 1001 !no interaction
      endif

      !ikoevt=nprt(1)
      izp=0
      izh=0
      if(iappl.le.3)qsqevt=0.
      do i=1,nptl
       if(istptl(i).eq.30.or.istptl(i).eq.31)then
         izp=izp+1
         if(ityptl(i)/10.eq.3)then
           izh=izh+1
           if(iappl.le.3)qsqevt=qsqevt+qsqptl(i)
         endif
       endif
      enddo
      ikoevt=izp
      ikhevt=izh
      f_ikhevt=ikhevt
      call eventvariset(4,f_ikhevt) !z4zevt
      if(iappl.le.3.and.izh.gt.0)qsqevt=qsqevt/float(izh)
      if(izmode.eq.2
     .  .or.ikolmn.gt.0.or.ikolmx.lt.10000)then
        i1=ikolmn
        i2=ikolmx
        if(nfr.gt.0)then
        i1=ikoevtxx
        i2=ikoevtxx
        endif
        izp=ikoevt
        if(izp.lt.i1.or.izp.gt.i2)iret=-2
      endif
      ptrevt=0
      do i=maproj+matarg+1,minfra
      if(istptl(i).eq.25)then
        pt=sqrt(pptl(1,i)**2+pptl(2,i)**2)
        ptrevt=max(ptrevt,pt)
      endif
      enddo
      if(ptrevt.lt.ptrmin)then
        iret=-2
        goto1000
      endif
      !print*,'+++++++++++++++++',ikoevt,ptrevt

1000  continue
      if(iret.eq.0.or.iret.eq.55)call emsaaFin
      if(mod(ntry,20000).eq.0)write(ifmt,*)ntry,'th iteration, ',
     .izp,' collisions,  iret=',iret,'  b=',bimevt
      if(iret.eq.0.and.(ptrmin.gt.0..or.
     .  ihlle.eq.1.or.ispherio.eq.1
     .  .or.ikolmn.gt.10))then
        write(ifmt,'(a,i3,a,i7,a,i5,a,f6.2,a,f7.1)')
     . 'ini',ninx,':',
     .  ntry,' its to get',izp,' Poms,  b =',bimevt
     .,' ,  ptrevt =',ptrevt
      endif
      if(izmode.eq.3.and.iret.eq.0.and.(nglmin.gt.0.or.nglmax.le.99999))
     .  write(ifmt,*)'ninx =',ninx,': ',
     .  ntry,' attempt(s) to get ',nglevt,' NN collision(s) for b ='
     .  ,bimevt
      call utprix('emsaaa',ish,ishini,4)
      return
1001  iret=-1           !no interaction
      goto 1000

      end

c----------------------------------------------------------------------
      subroutine dumpList
c----------------------------------------------------------------------
#include "aaa.h"
      write(ifmt,'(a,$)')'dump cptl ...'
      open(17,file=fndt(1:nfndt-4)//'dump',status='unknown')
      write(17,*)nptl
      do n=1,nptl
      write(17,*)
     . (pptl(j,n),j=1,5)
     ., iorptl(n)
     ., idptl(n)
     ., istptl(n)
     .,(tivptl(j,n),j=1,2)
     .,(ifrptl(j,n),j=1,2)
     ., jorptl(n)
     .,(xorptl(j,n),j=1,4)
     .,(ibptl(j,n),j=1,4)
     ., ityptl(n)
     ., itsptl(n)
     ., iaaptl(n)
     ., radptl(n)
     ., desptl(n)
     ., dezptl(n)
     ., rinptl(n)
     ., qsqptl(n)
     .,(zpaptl(j,n),j=1,2)
      enddo
      close (17)
      write(ifmt,'(a)')'  done'
      end

c----------------------------------------------------------------------
      subroutine readList(nfr)
c----------------------------------------------------------------------
#include "aaa.h"
      if(modsho.eq.1)
     *write(ifmt,'(a,i3,a,$)')'read cptl for nfr =',nfr,' ...'
      open(17,file=fndt(1:nfndt-4)//'dump',status='unknown')
      nevt=1
      read(17,*)nptl
      do n=1,nptl
      read(17,*)
     . (pptl(j,n),j=1,5)
     ., iorptl(n)
     ., idptl(n)
     ., istptl(n)
     .,(tivptl(j,n),j=1,2)
     .,(ifrptl(j,n),j=1,2)
     ., jorptl(n)
     .,(xorptl(j,n),j=1,4)
     .,(ibptl(j,n),j=1,4)
     ., ityptl(n)
     ., itsptl(n)
     ., iaaptl(n)
     ., radptl(n)
     ., desptl(n)
     ., dezptl(n)
     ., rinptl(n)
     ., qsqptl(n)
     .,(zpaptl(j,n),j=1,2)
      enddo
      close (17)
      if(modsho.eq.1)
     .write(ifmt,'(a)')'  done'
      end

c----------------------------------------------------------------------
      subroutine alist(text,n1,n2)
c----------------------------------------------------------------------
c    ior  jor  i  ifr1  ifr2     id  ist  ity      pt  m  y
c----------------------------------------------------------------------
c       ist                                     ity
c                                  light cluster ........ 19
c   ptl ...  0 1                   soft pom ............. 20-23   25(reggeon)
c   clu ... 10 11                  hard pom low mass .... 30
c   ptn ... 20 21                  proj remnant ......... 40-49
c   str ... 29                     targ remnant ......... 50-59
c   pom ... 30 31 32(virtual)      cluster .............. 60
c   rem ... 40 41                  direct photon ........ 71,72
c----------------------------------------------------------------------
#include "aaa.h"
      common/cxyzt/xptl(mxptl+29),yptl(mxptl+29),zptl(mxptl+29)
     * ,tptl(mxptl+29),optl(mxptl+29),uptl(mxptl+29),sptl(mxptl+29)
     *,rptl(mxptl+29,3)
      character  text*(*)
      dimension pp(5)
      if(n1.gt.n2)return
      imax=index(text,'&')
      if(imax.gt.1)then
      write(ifch,'(/1x,72a1/1x,a,a,a,72a1)')
     *('#',k=1,72),'############  ',text(1:imax-1),'  '
     *,('#',k=1,57-imax)
      write(ifch,'(1x,72a1/)')('#',k=1,72)
      endif
      if(n1.eq.0.and.n2.eq.0)return
      if(imax.gt.1)then
      write(ifch,'(1x,a,a/1x,89a1)')
     *'    ior    jor      i   ifr1   ifr2       id ist ity',
     *'        pt         m         E         y'
     *,('-',k=1,89)
      endif

      do j=1,5
        pp(j)=0.
      enddo
      nqu=0
      nqd=0
      nqs=0

      do iloo=n1,n2

        i=iloo
        if(iloo.gt.mxptl)then
          call restorecccptl(iloo,mxptl+2)
          i=mxptl+2
        endif

        ptptl=pptl(1,i)**2.+pptl(2,i)**2.
        if(ptptl.le.0.)then
          ptptl=0.
        else
          ptptl=sqrt(ptptl)
        endif
        amtptl=pptl(1,i)**2.+pptl(2,i)**2.+pptl(5,i)**2.
        if(amtptl.le.0.)then
          amtptl=0.
          if(abs(idptl(i)).lt.10000)then
            call idmass(idptl(i),amtptl)
          endif
          amtptl=sqrt(amtptl*amtptl+pptl(1,i)**2.+pptl(2,i)**2.)
        else
          amtptl=sqrt(amtptl)
        endif
        jor=jorptl(i)
        ifr2=ifrptl(2,i)
ctp        if(istptl(i).eq.20.or.istptl(i).eq.21)then
ctp        jor=nint(radptl(i))
ctp        if(jor.ne.0)jor=jor*100+99
ctp        endif
        rap=0.
        if(amtptl.gt.0..and.pptl(4,i).gt.0.)
     &  rap=sign(1.,pptl(3,i))*alog((pptl(4,i)+abs(pptl(3,i)))/amtptl)
         write(ifch,'(1x,i7,i7,i7,i7,i7,i10,2i3,2x,4(e9.3,1x),$)')
     &          iorptl(i),jor,iloo,ifrptl(1,i),ifr2
     &    ,idptl(i),istptl(i),ityptl(i),ptptl,pptl(5,i),pptl(4,i),rap
c        write(ifch,*)' '
        write(ifch,'(3x,4(e9.3,1x),e9.3)')
     &   (xorptl(jj,i),jj=3,4),qsqptl(i)    !,tivptl(1,i)
c        if(istptl(i).ne.12)write(ifch,*)' '
c        if(istptl(i).eq.12)write(ifch,'(1x,3(e9.3,1x))')
c     &           sptl(i),sqrt(uptl(i)-xorptl(1,i)**2)
c     &            ,sqrt(optl(i)-xorptl(2,i)**2)
        if(mod(istptl(i),10).eq.0.and.n1.eq.1.and.n2.eq.nptl)then
          do j=1,4
            pp(j)=pp(j)+pptl(j,i)
          enddo
        endif

        if(iloo.gt.mxptl)call dumpcccptl(mxptl+2,iloo)

      enddo

      end

c----------------------------------------------------------------------
      subroutine alist3(text,i)
c---------------------------------------------------------------------- 
      character  text*(*)
      imax=index(text,'&')
      if(imax.gt.1)then
        call getMonitorFileIndex(ifmtx)
        call getidptl(i,id)
        call getpptl(i,p1,p2,p3,p4,p5)
        call getxorptl(i,x1,x2,x3,x4)
        write(ifmtx,'(a,a,i5,i10,2(4e11.3,2x))') text(1:imax-1)
     .  , '  i id  x p = ' , i , id
     .  , x1,x2,x3,x4 , p1,p2,p3,p4
      else
        stop'in alist3, use & in text'
      endif
      end

c----------------------------------------------------------------------
      subroutine blist(text,n1,n2)
c----------------------------------------------------------------------
#include "aaa.h"
      character  text*(*)
      dimension pp(5)
      if(n1.gt.n2)return
      imax=index(text,'&')
      if(imax.gt.1)then
      write(ifch,'(/1x,89a1/1x,a,a,a,90a1)')
     *('#',k=1,89),'#############  ',text(1:imax-1),'  '
     *,('#',k=1,74-imax)
      write(ifch,'(1x,89a1/)')('#',k=1,89)
      endif
      if(n1.eq.0.and.n2.eq.0)return
      if(imax.gt.1)then
      write(ifch,'(1x,a,a,a/1x,90a1)')
     *'   ior   jor     i  ifr1   ifr2      id ist ity',
     *'        pt      mass    energy','       rap'
     *,('-',k=1,90)
      endif

      do j=1,5
        pp(j)=0.
      enddo
      nqu=0
      nqd=0
      nqs=0
      nqc=0
      nqb=0
      do i=n1,n2
        amtptl=pptl(1,i)**2.+pptl(2,i)**2.+pptl(5,i)**2.
        if(amtptl.le.0.)then
          amtptl=0.
          if(abs(idptl(i)).lt.10000)then
            call idmass(idptl(i),amtptl)
          endif
          amtptl=sqrt(amtptl*amtptl+pptl(1,i)**2.+pptl(2,i)**2.)
        else
          amtptl=sqrt(amtptl)
        endif
        pt=pptl(1,i)**2.+pptl(2,i)**2.
        if(pt.gt.0.)pt=sqrt(pt)
        rap=0.
        if(amtptl.gt.0..and.pptl(4,i).gt.0.)
     &  rap=sign(1.,pptl(3,i))*alog((pptl(4,i)+abs(pptl(3,i)))/amtptl)
        write(ifch,125)iorptl(i),jorptl(i),i,ifrptl(1,i),ifrptl(2,i)
     &       ,idptl(i),istptl(i),ityptl(i)
     &       ,pt,pptl(5,i),pptl(4,i),rap
 125  format (1x,i6,i6,i6,i6,i6,i10,2i3,2x,5(e9.3,1x)
     *     ,f9.2,4x,5(e8.2,1x))
      if(mod(istptl(i),10).eq.0.and.idptl(i).gt.-100000
     *                      .and.n1.eq.1.and.n2.eq.nptl)then
          do j=1,4
            pp(j)=pp(j)+pptl(j,i)
          enddo
        endif
      enddo
      write(ifch,'(90a1)')('-',k=1,90)
      write(ifch,125)0,0,0,0,0,0,0,0
     & ,sqrt(pp(1)**2+pp(2)**2),pp(3),pp(4)
     & ,sqrt(max(0.,pp(4)-pp(3))*max(0.,pp(4)+pp(3))-pp(1)**2-pp(2)**2)
      write(ifch,*)' '
      end

c----------------------------------------------------------------------
      subroutine clist(text,n1,n2,ity1,ity2)
c----------------------------------------------------------------------
#include "aaa.h"
c      parameter(itext=40)
      character  text*(*)
      dimension pp(5)
      if(n1.gt.n2)return
      imax=index(text,'&')
      if(imax.gt.1)then
      write(ifch,'(/1x,a,a,a,90a1)')
     *'-------------  ',text(1:imax-1),'  ',('-',k=1,74-imax)
      endif
      if(n1.eq.0.and.n2.eq.0)return
      if(imax.gt.1)then
      write(ifch,'(1x,a,a/1x,90a1)')
     *'     i       id ist ity',
     *'        pt        pz        p0      mass'
     *,('-',k=1,90)
      endif

      do j=1,5
        pp(j)=0.
      enddo
      do i=n1,n2
        pt=sqrt(pptl(1,i)**2+pptl(2,i)**2)
        write(ifch,127)i,idptl(i),istptl(i),ityptl(i)
     &       ,pt,pptl(3,i),pptl(4,i),pptl(5,i)
 127    format (1x,i6,i10,2i3,2x,4(e9.3,1x))
        if(ityptl(i).ge.ity1.and.ityptl(i).le.ity2)then
          do j=1,4
            pp(j)=pp(j)+pptl(j,i)
          enddo
        endif
      enddo
      write(ifch,'(90a1)')('-',k=1,90)
      write(ifch,127)0,0,0,0
     & ,sqrt(pp(1)**2+pp(2)**2),pp(3),pp(4)
     & ,sqrt(max(0.,pp(4)-pp(3))*max(0.,pp(4)+pp(3))-pp(1)**2-pp(2)**2)
      write(ifch,*)' '
      end

c----------------------------------------------------------------------
      subroutine alistf(text)
c----------------------------------------------------------------------
#include "aaa.h"
      character  text*(*)
      dimension pp(5),erest(5),errp(4)
      n1=1
      if(iframe.eq.21.and.(abs(iappl).eq.1.or.iappl.eq.3))
     *n1=2*(maproj+matarg+1)
      n2=nptl
      imax=index(text,'&')
      if(imax.gt.1)then
      write(ifch,'(/1x,124a1/1x,a,a,a,108a1)')
     *('#',k=1,124),'#############  ',text(1:imax-1),'  '
     *,('#',k=1,108-imax)
      write(ifch,'(1x,124a1/)')('#',k=1,124)
      endif
      if(imax.gt.1)then
      write(ifch,'(1x,a,a,a/1x,124a1)')
     *'   ior   jor        i     ifr1   ifr2         id ist ity',
     *'            px         py         pz         p0       mass',
     *'       rap'
     *,('-',k=1,124)
      endif

      do j=1,4
        pp(j)=0.
        errp(j)=0.
      enddo
      pp(5)=0.
      do i=n1,n2
        if(mod(istptl(i),10).eq.0)then
        amtptl=pptl(1,i)**2.+pptl(2,i)**2.+pptl(5,i)**2.
        if(amtptl.le.0.)then
          amtptl=0.
          if(abs(idptl(i)).lt.10000)then
            call idmass(idptl(i),amtptl)
          endif
          amtptl=sqrt(amtptl*amtptl+pptl(1,i)**2.+pptl(2,i)**2.)
        else
          amtptl=sqrt(amtptl)
        endif
        rap=0.
        if(amtptl.gt.0..and.pptl(4,i).gt.0.)
     &  rap=sign(1.,pptl(3,i))*alog((pptl(4,i)+abs(pptl(3,i)))/amtptl)
        write(ifch,125)iorptl(i),jorptl(i),i,ifrptl(1,i),ifrptl(2,i)
     &       ,idptl(i),istptl(i),ityptl(i),(pptl(j,i),j=1,5),rap
c     &,(xorptl(j,i),j=1,4)
        do j=1,4
          pp(j)=pp(j)+pptl(j,i)
        enddo
        endif
      enddo
 125  format (1x,i6,i6,3x,i6,3x,i6,i6,i12,2i4,4x,5(e10.4,1x)
     *     ,f9.2,4x,4(e8.2,1x))
 126  format (51x,5(e10.4,1x))
 128  format (51x,65('-'))
      pp(5)=(pp(4)-pp(3))*(pp(4)+pp(3))-pp(2)**2-pp(1)**2
      if(pp(5).gt.0.)then
        pp(5)=sqrt(pp(5))
      else
        pp(5)=0.
      endif
      write (ifch,128)
      write (ifch,126) (pp(i),i=1,5)
      erest(1)=0.
      erest(2)=0.
      if(iframe.eq.22.and.(abs(iappl).eq.1.or.iappl.eq.3))then
        i=maproj+matarg+1
        erest(3)=pptl(3,i)+matarg*pptl(3,i+1)
        erest(4)=pptl(4,i)+matarg*pptl(4,i+1)
      else
        erest(3)=maproj*pptl(3,1)+matarg*pptl(3,maproj+1)
        erest(4)=maproj*pptl(4,1)
     &          +matarg*pptl(4,maproj+1)
      endif
      erest(5)=amproj
      write (ifch,129)  (erest(j),j=1,5)
 129  format (50x,'(',5(e10.4,1x),')')
      do j=1,4
      if(abs(pp(j)).gt.0.d0)errp(j)=100.*(pp(j)-erest(j))/pp(j)
      enddo
      write (ifch,130)  (errp(j),j=1,4)
 130  format (50x,'(',3x,4(f7.2,4x),2x,'err(%))')
      end

c----------------------------------------------------------------------
      subroutine alist2(text,n1,n2,n3,n4)
c----------------------------------------------------------------------
#include "aaa.h"
      character  text*(*)
      if(n1.gt.n2)return
      imax=index(text,'&')
      write(ifch,'(1x,a,a,a)')
     *'--------------- ',text(1:imax-1),' ---------------  '
      do i=n1,n2
      write(ifch,125)iorptl(i),jorptl(i),i,ifrptl(1,i),ifrptl(2,i)
     &,idptl(i),istptl(i),ityptl(i),(pptl(j,i),j=1,5)
c     &,(xorptl(j,i),j=1,4)
      enddo
      write(ifch,'(1x,a)')'----->'
      do i=n3,n4
      write(ifch,125)iorptl(i),jorptl(i),i,ifrptl(1,i),ifrptl(2,i)
     &,idptl(i),istptl(i),ityptl(i),(pptl(j,i),j=1,5)
c     &,(xorptl(j,i),j=1,4)
      enddo
 125  format (1x,i6,i6,3x,i6,3x,i6,i6,i12,2i4,4x,5(e8.2,1x))
c     *,4x,4(e8.2,1x))
      end

c----------------------------------------------------------------------
      subroutine alistc(text,n1,n2)
c----------------------------------------------------------------------
#include "aaa.h"
      character  text*(*)
      if(n1.gt.n2)return
      imax=index(text,'&')
      if(n1.ne.n2)write(ifch,'(1x,a,a,a)')
     *'--------------- ',text(1:imax-1),' ---------------  '
      do i=n1,n2
      write(ifch,130)iorptl(i),jorptl(i),i,ifrptl(1,i),ifrptl(2,i)
     &,idptl(i),istptl(i),ityptl(i),(pptl(j,i),j=1,5)
     &,(xorptl(j,i),j=1,4),tivptl(1,i),tivptl(2,i)
      enddo
 130  format (1x,i6,i6,3x,i6,3x,i6,i6,i12,2i4,4x,5(e8.2,1x)
     *,4x,6(e8.2,1x))
      end

c-----------------------------------------------------------------------
      subroutine clop(n)
c-----------------------------------------------------------------------
#include "aaa.h"
      if(ifmt.eq.6)return
      if(n.eq.1)then
        open(ifmt,file=fnmt(1:nfnmt),access='append')
      elseif(n.eq.2)then
        close(ifmt)
      elseif(n.eq.3)then
        close(ifmt)
        open(ifmt,file=fnmt(1:nfnmt),access='append')
      endif
      end

c-----------------------------------------------------------------------
      integer function ifebin()
c-----------------------------------------------------------------------
#include "aaa.h"
      ifebin=abs(noebin)
      end

c-----------------------------------------------------------------------
      integer function ifopen()
c-----------------------------------------------------------------------
#include "aaa.h"
      ifopen=nopen
      end

c-----------------------------------------------------------------------
      integer function ifappl()
c-----------------------------------------------------------------------
#include "aaa.h"
      ifappl=iappl
      end

c-----------------------------------------------------------------------
      integer function ifnoevt()
c-----------------------------------------------------------------------
      common/photrans/phoele(4),ebeam,noevt
      ifnoevt=noevt
      end

c-----------------------------------------------------------------------
      integer function ifevent()
c-----------------------------------------------------------------------
#include "aaa.h"
      common/cmodshox/modshox
      ifevent=max(1,nevent)
      if(igrTree.gt.0.and.muTree.gt.0)then
        write(ifmt,'(2a)')'igrTree,muTree>0 -> changed to zero events'
     .   ,' ; modsho reset'
        ifevent=0
        modsho=modshox
      endif
      end

c-----------------------------------------------------------------------
      subroutine setebin(nrebinx)
c-----------------------------------------------------------------------
#include "aaa.h"
      nrebin=nrebinx
      end

c-----------------------------------------------------------------------
      subroutine getReaction(iprojZ,iprojA,itargZ,itargA,fegyevt)
c-----------------------------------------------------------------------
#include "aaa.h"
      integer  iprojZ,iprojA,itargZ,itargA
      real fegyevt
      iprojZ=laproj !projectile Z
      iprojA=maproj !projectile A
      itargZ=latarg !target Z
      itargA=matarg !target A
      fegyevt=egyevt !pp cm energy (hadron) or string energy (lepton)
      end

c-----------------------------------------------------------------------
      subroutine getifch(ivalue)
c-----------------------------------------------------------------------
#include "aaa.h"
      ivalue=ifch
      end

c-------------------------------------------------------
      subroutine getecc(psi2,psi3,psi4,psi5,ecci2,ecci3,ecci4,ecci5)
c-----------------------------------------------------------------------
      ! psiM  = eccphi(M) = polar(<cos(M*phi)>,<sin(M*phi)>) / M  - pi/M
      ! ecciM = ecccoe(M) = sqrt(<cos(M*phi)>**2+<sin(M*phi)>**2)
      !   polar(x,y) = polar angle of complex number x+i*y
      !------------------------------------------------------------

      common/ciniflo/ecccoe(5),eccphi(5)
      psi2=eccphi(2)
      psi3=eccphi(3)
      psi4=eccphi(4)
      psi5=eccphi(5)
      ecci2=ecccoe(2)
      ecci3=ecccoe(3)
      ecci4=ecccoe(4)
      ecci5=ecccoe(5)
      end

c-----------------------------------------------------------------------
      subroutine getnuclei(laprojx,maprojx,latargx,matargx)
c-----------------------------------------------------------------------
#include "aaa.h"
      laprojx=laproj
      maprojx=maproj
      latargx=latarg
      matargx=matarg
      end

c-----------------------------------------------------------------------
c Interface routine getevt: get event info from EPOS
c Usage:
c   call getevt(nev,phi,phir,bim,egy,npt,ngl,kol)
c Variables:
c   nev ........ error code. 1=valid event, 0=invalid event
c   bim ........ absolute value of impact parameter
c   phi ........ angle of impact parameter
c   phir ....... angle of the string segment based event plane (n=2)
c                 with respect to the impact parameter vector
c                 if iranphi=1
c                phir = 0 if iranphi=0
c   npt ........ number of particles
c   ngl ........ number of Glauber collisions
c   kol ........ number of real EPOS collisions
c-----------------------------------------------------------------------
      subroutine getevt(nev,phi,phir,bim,egy,npl,ngl,kol)
#include "aaa.h"
#include "ems.h"
      common/cgetevt/kFkol(kollmx/4),npoms(kollmx/4)
     . ,nFkolni(kollmx/4,100)
      common/cranphi/ranphi
      nev=nevt
      phi=phievt
      phir=ranphi
      bim=bimevt
      egy=egyevt
      npl=nptl
      ngl=nglevt
      kol=0
      do k=1,koll
        npo=0
        do n=1,nprmx(k)
          if(idpr(n,k).gt.0)npo=npo+1
          !if(idpr(n,k).gt.0.and.shatpr(n,k).gt.400.)npo=npo+1
        enddo
        if(npo.gt.0)then
          kol=kol+1
          if(kol.gt.kollmx/4)stop'\n\n ERROR 28072011 \n\n'
          kFkol(kol)=k
          npoms(kol)=npo
          ni=0
          do n=1,nprmx(k)
            if(idpr(n,k).gt.0)then
              ni=ni+1
              nFkolni(kol,ni)=n
            endif
          enddo
        endif
      enddo
      end

c-----------------------------------------------------------------------
c Interface routine getptl: get particle list from EPOS
c Usage:
c   call getevt(nev,phi,phir,bim,egy,npt,ngl,kol) ! to get npt
c   do i=1,npt                                
c     call getptl(imo,ifa,i,ic1,ic2            
c  .  ,id,ist,ity,p1,p2,p3,p4,p5,x1,x2,x3,x4)  
c     ...
c   enddo
c Arguments:
c   imo ..... index of mother
c   ifa ..... index of father
c   i ....... index of particle
c   ic1 ..... index of first child
c   ic2 ..... index of last child
c   id ...... particle id
c   ist ..... status: 0 = last generation particle
c   ity ..... type of origin
c   p1 ..... x-component of particle momentum
c   p2 ..... y-component of particle momentum
c   p3 ..... z-component of particle momentum
c   p4 ..... particle energy
c   p5 ..... particle mass
c   x1 ... x-component of formation point
c   x2 ... y-component of formation point
c   x3 ... z-component of formation point
c   x4 ... formation time
c-----------------------------------------------------------------------
      subroutine getptl(imo,ifa,iloo,ic1,ic2
     . ,id,ist,ity,p1,p2,p3,p4,p5,x1,x2,x3,x4)
#include "aaa.h"
      j=iloo
      if(iloo.gt.mxptl)then
        call restorecccptl(iloo,mxptl+2)
        j=mxptl+2
      endif
      i=j
      imo=iorptl(i)
      ifa=jorptl(i)
      ic1=ifrptl(1,i)
      ic2=ifrptl(2,i)
      id=idptl(i)
      ist=istptl(i)
      ity=ityptl(i)
      p1=pptl(1,i)
      p2=pptl(2,i)
      p3=pptl(3,i)
      p4=pptl(4,i)
      p5=pptl(5,i)
      x1=xorptl(1,i)
      x2=xorptl(2,i)
      x3=xorptl(3,i)
      x4=xorptl(4,i)
      end

c-----------------------------------------------------------------------
c Interface routine getpos: get Pomeron positions from EPOS
c Usage:
c   call getevt(nev,phi,phir,bim,egy,npt,ngl,kol) ! to get kol
c   do ko=1,kol
c    do ni=1,ifnpom(ko)
c     call getpos(ko,ni,x,y,z,t)
c    enddo
c   enddo
c-----------------------------------------------------------------------
      subroutine getpos(kol,ni,x,y,z,t)
#include "aaa.h"
#include "ems.h"
      common/cgetevt/kFkol(kollmx/4),npoms(kollmx/4)
     . ,nFkolni(kollmx/4,100)
      common/cranphi/ranphi
      phi=   phievt+ranphi
      k=kFkol(kol)
      n=nFkolni(kol,ni)
      if(idpr(n,k).le.0)stop'\n\n ERROR 28072011b \n\n '
      x=  coordpr(1,n,k)*cos(phi)+coordpr(2,n,k)*sin(phi)
      y= -coordpr(1,n,k)*sin(phi)+coordpr(2,n,k)*cos(phi)
      z= coord(3,k)
      t= coord(4,k)
      end
      integer function ifnpom(kol)
#include "aaa.h"
#include "ems.h"
      common/cgetevt/kFkol(kollmx/4),npoms(kollmx/4)
     . ,nFkolni(kollmx/4,100)
      ifnpom=npoms(kol)
      end

c-----------------------------------------------------------------------
c Interface routine gethypar: get hydro table from EPOS
c Usage:
c   call gethypar(nzhyx,ntauhyx,nxhyx,nyhyx
c  . ,zminhyx,zmaxhyx,tauminhyx,taumaxhyx
c  . ,xminhyx,xmaxhyx,yminhyx,ymaxhyx)
c   do neta=1,nzhyx
c    do ntau=1,ntauhyx
c     do nx=1,nxhyx
c      do ny=1,nyhyx
c       call gethyval(neta,ntau,nx,ny,eps,tem,pss,v1,v2,v3)
c        ...
c      enddo
c     enddo
c    enddo
c   enddo
c-----------------------------------------------------------------------
      subroutine gethypar(nzhyx,ntauhyx,nxhyx,nyhyx
     .                ,zminhyx,zmaxhyx,tauminhyx,taumaxhyx
     .                ,xminhyx,xmaxhyx,yminhyx,ymaxhyx)
#include "aaa.h"
#include "ho.h"
      nzhyx    = nzhy
      ntauhyx  = ntauhy
      nxhyx    = nxhy
      nyhyx    = nyhy
      zminhyx  = zminhy
      zmaxhyx  = zmaxhy
      tauminhyx= tauminhy
      taumaxhyx= taumaxhy
      xminhyx  = xminhy
      xmaxhyx  = xmaxhy
      yminhyx  = yminhy
      ymaxhyx  = ymaxhy
      end
      subroutine gethyval(neta,ntau,nx,ny,eps,tem,pss,v1,v2,v3)
#include "aaa.h"
#include "ho.h"
      double precision temc
      common/ctemc/temc(netahxx,ntauhxx,nxhxx,nyhxx)
      common/cpssc/pssc(netahxx,ntauhxx,nxhxx,nyhxx)
      eps= epsc(neta,ntau,nx,ny)
      v1=velc(1,neta,ntau,nx,ny)
      v2=velc(2,neta,ntau,nx,ny)
      v3=velc(3,neta,ntau,nx,ny)
      b1=barc(1,neta,ntau,nx,ny)
      b2=barc(2,neta,ntau,nx,ny)
      b3=barc(3,neta,ntau,nx,ny)
      tem=temc(neta,ntau,nx,ny)
      pss=pssc(neta,ntau,nx,ny)
      end

      subroutine gettauzer(tau0)
      real         tauzer,tauone,tautwo,deltau,numtau,amsiac,amprif
     *,tauup,tauthree
      common/resc1/tauzer,tauone,tautwo,deltau,numtau,amsiac,amprif
     *,tauup,tauthree
      double precision tau0
      tau0=real(tauzer)
      end


c-----------------------------------------------------------------------
      subroutine iniRan
c-----------------------------------------------------------------------
c this subroutine should not be in gen if we want to link to another
c code with different random number generator such that gen.f is not
c used at all.
c-----------------------------------------------------------------------
#include "aaa.h"

      if(noebin.ge.0)then
        ntevt=0
        if(seedi.ne.0d0)then
          call ranfini(seedi,iseqini,1)
        else
          stop 'seedi = 0 ... Please define it !'
        endif
        seedc=seedi
        if(inicnt.eq.1)call aseedi
      elseif(inicnt.eq.1)then !fake DIS, initialization is part of the event
        if(seedj.ne.0d0)then
          if(seedj2.ne.0d0)then
            call ranfcv(seedj2)
            write(ifmt,'(a)')
     & "Random number sequence does not start at 0 ... please wait !"
          endif
          call ranfini(seedj,iseqsim,2)
        else
          stop 'seedj = 0 ... Please define it !'
        endif
        call aseed(2)
      endif

      end

c-----------------------------------------------------------------------
c get functions
c-----------------------------------------------------------------------

      subroutine getQcdlam(val)
      real val 
#include "sem.h"
      val=qcdlam
      end

      subroutine getQ2ini(val)
      real val 
#include "sem.h"
      val=q2ini
      end

      subroutine getEpmax(val)
      real val 
      common /psar2/  edmax,epmax
      val=epmax
      end


c-------------------------------BASIC----------------------------------------
#if __BS__ 
      subroutine athermal(n)
      m=n
      end
      subroutine dih
      end
      subroutine iniEos
      end
      subroutine eostabs(i,j,k)
      m=i
      m=j
      m=k
      end
c in version 3454d removed commented subroutines hdecin,hqini,exiteposevent,hdecas
      subroutine d2hclose
      end
#endif
c-------------------------------VERY BASIC----------------------------------
#if __BS__ && !__ROOT__
      subroutine getNames(n)
      m=n
      end
      subroutine pthardparent(i,j,k,x)
      m=i
      m=j
      m=k
      y=x
      end
      subroutine treeclose
      end
      subroutine ordertree(n)
      m=n
      end
      function fmux(q,w,e,r,t,y,u,o,p,a,s)
      a=q
      a=w
      a=e
      a=r
      a=t
      a=y
      a=u
      a=o
      a=p
      a=a
      a=s
      fmux=a
      end
#endif


