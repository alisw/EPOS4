C
C  This file is part of EPOS4
C  Copyright (C) 2022 research institutions and authors (See CREDITS file)
C  This file is distributed under the terms of the GNU General Public License version 3 or later
C  (See COPYING file for the text of the licence)
C

c------------------------------------------------------------------------
      subroutine eosInv(e, nb, t, mub, p)
c------------------------------------------------------------------------
c Universal inverted EoS function
c compute t, mub, p
c------------------------------------------------------------------------
#include "aaa.h"
      double precision e,nb,nq,ns,p,t,mub,muq,mus
      if(ioeos.eq.22)then     ! x3ff
        !muq=...
        !mus=... 
        call eoshlle (e, nb, nq, ns, t, mub, muq, mus, p) !interpolated
        !call eosohlle(e, nb, nq, ns, t, mub, muq, mus, p) !computed
      endif
      end

c------------------------------------------------------------------------
      subroutine eosOri(t, mub, e, nb, p)
c------------------------------------------------------------------------
c Universal original EoS function
c compute e, nb, p
c------------------------------------------------------------------------
#include "aaa.h"
      double precision e,nb,nq,ns,p,t,mub,muq,mus
      if(ioeos.eq.22)then     ! x3ff
        !nq=...
        !ns=...
        call eosihlle(t, mub, muq, mus, e, nb, nq, ns, p) 
      endif
      end

c---------------------------------------------------------------------
      subroutine iniEos
c---------------------------------------------------------------------
#include "aaa.h"
#include "ho.h"
      if(ioeos.eq.22)then     ! x3ff
        fn3f='z-eos4f.eos'
        fn3g=fnnx(1:nfnnx)//'src/KWt/eos1f.eos'
        if(ieostabm.eq.1)call tctabm !make table1
        if(ieostabm.eq.2)call eostabm !make table2
      endif
      !here nothing needed for other ioeos
      call seteostype(ioeos) !just to allow code in YK only for 22
      end

c---------------------------------------------------------------------
      subroutine closeEos(nfr)
c---------------------------------------------------------------------
#include "aaa.h"
#include "ho.h"
      common/cistat/istateos,istatptl
      common/cntauhyori/ntauhyori
      call clop(3)
      if(nfr+1.eq.nfreeze
     ..and.istateos.eq.1)then
        call memo(1,'destroy eos table;')
        if(ioeos.eq.22)then
          call destroyeoshlle()
        endif
        call memo(2,';')
        istateos=0
      endif
      end

c---------------------------------------------------------------------
      subroutine closePtlDky
c---------------------------------------------------------------------
#include "aaa.h"
#include "ho.h"
      common/cistat/istateos,istatptl
      call memo(1,'destroy ptldky table;')
      call destroyptldky()
      call memo(2,';')
      istatptl=0
      end

c---------------------------------------------------------------------
       subroutine eosx
c---------------------------------------------------------------------
#include "aaa.h"
      if(ioeos.gt.0)then
       if(ireadhyt.eq.0)then
         if(ihlle.eq.0.and.ioeos.gt.0)then
           write(ifmt,'(a)')'enter hlle in eosx; nfr=-1 ... '
           call hlle(-1)
           write(ifmt,'(a)')'exit hlle in eosx'
         endif
         call xxEos
       endif
      endif    
      end

c-----------------------------------------------------------------------
      subroutine eosparams(ipr)
c-----------------------------------------------------------------------
#include "aaa.h"
      double precision  gammaS, B, mS, delta0, muc, volex0,
     . emax, e0, anmax, bnmax, an0, muBmax, muQmax, muSmax, Taccuracy,
     . muBmx, muQmx, muSmx,
     . aaa, bbb, temmax
      integer nxe, nxn, maxIter
      common /ceospar/  gammaS, B, mS, delta0, muc, volex0,
     . emax, e0, anmax, bnmax, an0, muBmax, muQmax, muSmax, Taccuracy,
     . muBmx, muQmx, muSmx,
     . aaa, bbb,  temmax,  nxe, nxn, maxIter
       

      gammaS = 1.00       ! gamma_s, strangeness suppression factor
      B = 0.20            ! bag const [GeV/fm^3]
      mS = 0.12           ! strange quark mass [GeV]
      volex0 = 1.50 
      delta0 =0.24  ! delta_0, as in in Hama et al. paper
      aaa = 0.77
      bbb = 3.0
      muc = 0.7  !0.4    ! mu_c [GeV], as in draft (critical

      !-----------------------------------------------------------------
      ! e discretization: 
      !   In es.cpp, one computes del_xe = log(emax/e0+1)/(nxe-1)
      !   Parameters emax and nxe may be increased, but such that del_xe 
      !   remains unchanged, to keep same accuracy for a larger range 
      !   So: nxe  =  1 + log(emax/e0+1) / del_xe(=0.130)
      ! n discretization:
      !   In es.cpp, one computes del_xn = 2.*log(nmax/n0+1)/(nxn-1)
      !-----------------------------------------------------------------
      emax = 5000. !500!40  ! max epsilon
      e0 = 0.01   !scale on which the points for epsilon
                  !start to compress to give better accuracy
      anmax = 0.8   ! max density:
      bnmax = 0.6   ! nmax = anmax * E^bnmax
      an0 = 0.0001  ! n0 =  an0 * (E/2)^bnmax   
      nxe = 102 !84!65   
      nxn = 45 !33 
      
      muBmax = 1.0   ! maximal baryon chem. potential
      muQmax = 0.3   ! maximal electric charge chem. potential
      muSmax = 0.5   ! maximal strange chem. potential
      !In principle bigger values of muBmax,muQmax,muSmax are needed 
      !to invert EoS, but then EoS is not calculable,
      !because exp(mu/T) becomes inf in HG_exv2.
      !After fixing this in an approximate way, inversion does not converge. 
      !It concern T close to T_c
      muBmx = 1.0 !>> affects <<!  values
      muQmx = 0.3 !>>   eos1  <<!   used 
      muSmx = 0.5 !>>  table  <<!  for Tc 
      temmax = 1.500     ! maximal temperature  
      maxIter = 20       ! maximal iterations in  inversion procedure
      Taccuracy = 0.000005 ! maximal T accuracy in inversion procedure

      if(ipr.eq.1)write(ifmt,*)'B =   ',B  ,'  volex0 = ',volex0
      if(ipr.eq.1)write(ifmt,*)'muc = ',muc,'  delta0 = ',delta0
      if(ipr.eq.1)write(ifmt,*)'aaa = ',aaa,'  bbb = ',bbb
            
      call  setparameters(gammaS, B, mS, volex0, delta0, aaa, bbb, muc,
     .  emax, e0, anmax, bnmax, an0, nxe, nxn,  muBmax, muQmax, muSmax,
     .  muBmx, muQmx, muSmx,  temmax, maxIter, Taccuracy)
     
      end

c-----------------------------------------------------------------------
      subroutine tctabm
c-----------------------------------------------------------------------
#include "aaa.h"
      character*80 fn
      double precision  gammaS, B, mS, delta0, muc, volex0,
     . emax, e0, anmax, bnmax, an0, muBmax, muQmax, muSmax, Taccuracy,
     . muBmx, muQmx, muSmx,
     . aaa, bbb, temmax
      integer nxe, nxn, maxIter
      common /ceospar/  gammaS, B, mS, delta0, muc, volex0,
     . emax, e0, anmax, bnmax, an0, muBmax, muQmax, muSmax, Taccuracy,
     . muBmx, muQmx, muSmx,
     . aaa, bbb,  temmax,  nxe, nxn, maxIter
      write(ifmt,'(a)')'making Tc tab ...'
      call eosparams(1)
      call eostabs(1,0,0)
      ii=index(fnhi(1:nfnhi),".histo")-1
      fn=fnhi(1:ii)//".eos"
      if(iopcnt.eq.0)then
        istart=0 
        iend=nxn-1
      else
        istart=  (iopcnt-1)*2
        iend=min((iopcnt-1)*2+1,nxn-1)
      endif
      call maketctable(fn(1:ii+4)//CHAR(0),istart,iend)
      stop'tctabm finished successfully'
      end

c-----------------------------------------------------------------------
      subroutine eostabs(i,j,k)
c-----------------------------------------------------------------------
#include "aaa.h"
      double precision  gammaS, B, mS, delta0, muc, volex0,
     . emax, e0, anmax, bnmax, an0, muBmax, muQmax, muSmax, Taccuracy,
     . muBmx, muQmx, muSmx,
     . aaa, bbb, temmax
      integer nxe, nxn, maxIter
      common /ceospar/  gammaS, B, mS, delta0, muc, volex0,
     . emax, e0, anmax, bnmax, an0, muBmax, muQmax, muSmax, Taccuracy,
     . muBmx, muQmx, muSmx,
     . aaa, bbb,  temmax,  nxe, nxn, maxIter
      common/cistat/istateos,istatptl
      common/cestype/iestype
      logical ok

      if(i.eq.1)then
       write(ifmt,'(2a)')'load ',fn3p(1:nfn3p) !ptltab
       write(ifmt,'(2a)')'load ',fn3d(1:nfn3d) !dkytab
       if(istatptl.eq.0)then
        call loadptldky(
     .  fn3p(1:nfn3p)//CHAR(0)
     . ,fn3d(1:nfn3d)//CHAR(0)) 
        istatptl=1
        np=igetnparticles()
        write(ifmt,'(i4,a)')np,' particles'
       else
        write(ifmt,'(a)')'ptl/dky tabs already loaded'
       endif
      endif

      if(j.eq.1)then
       if(ioeos.eq.22)then
        nfn3g=index(fn3g,' ')-1
        write(ifmt,'(2a)')'load ',fn3g(1:nfn3g)
        call initeos3tc(fn3g(1:nfn3g)//CHAR(0))
       endif
      endif

      if(k.eq.1)then
       if(istateos.eq.0)then
        istateos=1
        if(ioeos.eq.22)then
          nfn3f=index(fn3f,' ')-1
          inquire(file=fn3f(1:nfn3f),exist=ok)
          if(.not.ok)then
            write(ifmt,'(3a)')'ERROR: File ',fn3f(1:nfn3f)
     .      ,' does not exist'
            stop'ERROR file not found'
          endif
          write(ifmt,'(2a)')'load ',fn3f(1:nfn3f)
          call memo(1,'create eos via initeoshlle3f;')
          call initeoshlle3f(fn3f(1:nfn3f)//CHAR(0)
     .      ,B,volex0,delta0,aaa, bbb) 
        elseif(ioeos.eq.6)then  !--> eos = new eoBEST(PathBEST) in ctrl.cpp
          write(*,*)'reading BEST table'
          call initbest(fnnx(1:nfnnx)//'src/MSt/'//CHAR(0))
        elseif(ioeos.eq.8)then  !--> eos = new eoBEST(PathBEST2) in ctrl.cpp
          write(*,*)'reading BEST2 table'
          call initbest(fnnx(1:nfnnx)//'srcext/EOS/BEST2/'//CHAR(0)) !BEST2 to be created
        elseif(ioeos.eq.30)then
          call memo(1,'create eos via initeoschiralhlle;')
          call initeoschiralhlle(fnnx(1:nfnnx)//
     .    'MSt/chiralsmall.dat'//CHAR(0),
     .    fnnx(1:nfnnx)//'MSt/chiraleos.dat'//CHAR(0))
        elseif(ioeos.eq.35)then
          call initcem(fnnx(1:nfnnx)//'MSt/CEMmini.dat'//CHAR(0),
     .    '..MSt/CEMsmall.dat'//CHAR(0),
     .    '..MSt/CEMbig.dat'//CHAR(0))
        elseif(ioeos.eq.7)then
          write(*,*)'reading PJNL table'
          call initpnjl(fnnx(1:nfnnx)//'MSt/'//CHAR(0))
        else
          stop'in eostabs: wrong ioeos choice'
        endif
        call memo(2,';')
       else
        write(ifmt,'(a)')'already loaded'
       endif 
      endif

      end

c-----------------------------------------------------------------------
      subroutine eostabm
c-----------------------------------------------------------------------
#include "aaa.h"
      double precision  gammaS, B, mS, delta0, muc, volex0,
     . emax, e0, anmax, bnmax, an0, muBmax, muQmax, muSmax, Taccuracy,
     . muBmx, muQmx, muSmx,
     . aaa, bbb, temmax
      integer nxe, nxn, maxIter
      common /ceospar/  gammaS, B, mS, delta0, muc, volex0,
     . emax, e0, anmax, bnmax, an0, muBmax, muQmax, muSmax, Taccuracy,
     . muBmx, muQmx, muSmx,
     . aaa, bbb,  temmax,  nxe, nxn, maxIter
      integer  istart,iend
      character*80 fnm,fn
      common/cestype/iestype /cjjtb/jjtb,jjeos
      
      write(ifmt,'(a)')'making eos tab ...'

      call eosparams(1)
      call eostabs(1,1,0)
                         
      !-----------------------------------
      ! in case of jjtb=1
      !   one needs nxn * (nxe-1) jobs 
      !-----------------------------------
      ! in case of jjtb>1
      ! one runs packages of jjtb units
      ! using an increment of jjtb in cepos  
      !  (only every jjtb-th job is executed)
      ! using nskip
      !-----------------------------------

      if(jjtb.gt.1.or.jjeos.gt.0)then
        k=iopcnt-1
        if(mod(k,jjtb).ne.0)stop'####### ERROR 21032011B #######'
        nn=(nxn-1)/jjtb+1
        nxxn=nn*jjtb 
        if(nxxn.gt.jjeos)stop'####### ERROR 09072015 #######'
                !optns file not compatible with dimensions     
        istart= k/jjeos+1
        iend=   k/jjeos+1 
        jstart= mod(k,jjeos)
        jend=   min( mod(k,jjeos)+jjtb-1 , nxn-1 )
      elseif(jjtb.eq.1)then
        k=iopcnt-1
        istart= k/nxn+1
        iend=   k/nxn+1
        jstart= mod(k,nxn)
        jend=   mod(k,nxn)
      else
        k=iopcnt-1
        istart=k
        iend=k
        jstart=0
        jend=nxn-1
      endif
      if(iend.gt.nxe-1)
     .stop'\n\n STOP in eostabm: iend too large\n\n'
      
      
      if(ifmt.ne.6)then
        fnm=fnmt(1:nfnmt)
        n=nfnmt
      else
        n=0
      endif
      ii=index(fnhi(1:nfnhi),".histo")-1
      fn=fnhi(1:ii)//".eos"
      iixx=ii+4
      call clop(2)
      if(iopcnt.eq.1)then
      !first do complete table for istart=0 iend=0 (fast)
      call makeinvtable(
     .fnm(1:n)//CHAR(0),fn(1:iixx)//CHAR(0)
     .,0,0,0,nxn-1)
      endif
      !then for given istart,iend, tables for jstart to jend
      call makeinvtable(
     .fnm(1:n)//CHAR(0),fn(1:iixx)//CHAR(0)
     .,istart,iend,jstart,jend) 
      !makeinvtable uses solveMix   (T is in GeV) 
      call clop(1)
      write(ifmt,'(2a)')'eos tab written to ',fn(1:iixx)
      stop
      end

c-----------------------------------------------------------------------
      subroutine xxEos
c-----------------------------------------------------------------------
#include "aaa.h"
#include "ho.h"
      if(iHyEos.eq.0)return
      if(ioeos.eq.0.or.nrevt.gt.1)return
      if(ioeos.ge.1)then
c       call xxEos1
c       call xxEos2
        call xxEos3(0.0, 0.5,0.0)
        call xxEos3(0.0, 0.5,0.3)
        call xxEos3(0.0, 0.5,0.45)
        stop'Regular stop in xxEos'
      endif
      end
      
c-----------------------------------------------------------------------
      subroutine eostest
c-----------------------------------------------------------------------
#include "ho.h"
      external eosohlle,eoshlle
      if(ioeos.eq.0.or.nrevt.gt.0)return
c      call eostest1(eosohlle)
c      call eostest2
c      call eostest4
c      call eostest1(eoshlle)
c      call eostest5 
c      call eostest6 
c      call eostest7
c      call eostest8
c       stop'Normal stop in eostest'
      end           

c-----------------------------------------------------------------------
      subroutine xxEos1 ! plot  EoS  eps,n ( T, mu)       (T is in GeV)
c-----------------------------------------------------------------------
      double precision  gammaS, B, mS, delta0, muc, volex0,
     . emax, e0, anmax, bnmax, an0, muBmax, muQmax, muSmax, Taccuracy,
     . muBmx, muQmx, muSmx,
     . aaa, bbb, temmax
      integer nxe, nxn, maxIter
      common /ceospar/  gammaS, B, mS, delta0, muc, volex0,
     . emax, e0, anmax, bnmax, an0, muBmax, muQmax, muSmax, Taccuracy,
     . muBmx, muQmx, muSmx,
     . aaa, bbb,  temmax,  nxe, nxn, maxIter
      double precision mubx,muqx,musx
#include "aaa.h"
#include "so.h"
#include "ho.h"
      integer n(3,4)
      data n /1,0,0,  0,1,0,  0,0,1, 1,0,0 /
      character*8 cmd
      character*4 ck(4,4),ci(4)
      data ck/'[m]B','[m]Q','[m]S',' T  ', '[m]Q','[m]S','[m]B',' T  '
     .       ,'[m]S','[m]B','[m]Q',' T  ', ' T  ','[m]Q','[m]S','[m]B'/
      data ci / ' E  ', ' nB ', ' nQ ', ' nS ' /
      double precision ee, n_b, n_q, n_s, pp,ttt,mub, muq, mus
      
      mubx=muBmx/4*3
      muqx=muQmx/4*3
      musx=muSmx/4*3
      temx=temmax/15*4

      do iii=1,4     ! y variable (1-4)
      do jjj=1,4     ! x variable + 3 fixed variables (1-4)
      do kkk=1,2         
      write(ifhi,'(a)')'!##############################################'
      write(ifhi,'(a)')'!        eos  eps, n ( T, mu)        '
      write(ifhi,'(a)')'!##############################################'
      cmd='lin'
      if(iii.eq.1)cmd='log'
      if(jjj.le.3)then
      x1=0
      x2=temx
      elseif(jjj.eq.4)then
      x1=0
      x2=mubx
      endif
      do k=1,3
        if(jjj.eq.1)ww=-(kkk*2-3)*mubx/2.*(k-1)
        if(jjj.eq.2)ww=-(kkk*2-3)*muqx/2.*(k-1)
        if(jjj.eq.3)ww=-(kkk*2-3)*musx/2.*(k-1)
        if(jjj.eq.4)then
           if(kkk.eq.1.and.k.eq.1)ww=.040
           if(kkk.eq.1.and.k.eq.2)ww=.060
           if(kkk.eq.1.and.k.eq.3)ww=.080
           if(kkk.eq.2.and.k.eq.1)ww=.100
           if(kkk.eq.2.and.k.eq.2)ww=.500
           if(kkk.eq.2.and.k.eq.3)ww=.900
        endif
        xtx=0.20+0.27*(k-1)
        write(ifhi,'(a/ a,2e15.6,a/ a/ a/ a/ a,f5.2,a,f6.3,a/ a/ a)') 
     . 'openhisto htyp lin xmod lin ymod '//cmd
     . ,'xrange ',x1,x2,'  yrange auto auto'
     . ,'txt  "yaxis '//ci(iii)//'"'
     . ,'txt  "xaxis '//ck(4,jjj)//'"'
     . ,'text 0.05 0.90 "'//ck(1,jjj)//': "'
     . ,'text ',xtx,' 0.90 "',ww,'"'
     . ,'text 0.05 1.03 "'//ck(2,jjj)//'='//ck(3,jjj)//'=0 "'
     . ,'array 2'
        if(jjj.le.3)then
        nnn=1
        nxx=200
        elseif(jjj.eq.4)then
        nnn=1
        nxx=200
        endif
        do ix =nnn,nxx-1
         if(jjj.le.3)then
           ttt = ix*x2/(nxx-1)
           mub=n(1,jjj)*ww
           muq=n(2,jjj)*ww
           mus=n(3,jjj)*ww
           x=ttt
         elseif(jjj.eq.4)then
          ttt=ww
          mub= ix*x2/(nxx-1)
          muq=0
          mus=0
          x=mub
         endif
         if(x.gt.x1.and.x.lt.x2)then
           call eosihlle(ttt, mub, muq, mus, ee, n_b, n_q, n_s, pp) 
           if(iii.eq.1)y=ee
           if(iii.eq.2)y=n_b
           if(iii.eq.3)y=n_q
           if(iii.eq.4)y=n_s
           if(ci(iii).ne.' E  '.or.ck(4,jjj).ne.' T '
     .      .or. abs(y).gt.1e-10)
     .      write(ifhi,'(2e15.6)')x,y
         endif
        enddo   
        write(ifhi,'(a)') ' endarray closehisto '
        if(k.eq.3)write(ifhi,'(a)')'plot 0'
        if(k.ne.3)write(ifhi,'(a)')'plot 0-'
      enddo
      enddo
      enddo
      enddo
      end

c-----------------------------------------------------------------------
      subroutine xxEos2    !   plot/check  iEoS  T,mu ( eps, n) 
c-----------------------------------------------------------------------
      ! iEoS =  inverse EoS   (T is in GeV)
      !-----------------------------------------------------------------
      ! blue points: interpol iEoS using eoshlle (ee..) for grid points 
      ! yellow pnts: interpol iEoS using eoshlle (ee..) for more pnts(*) 
      ! red stars:   cumputes iEoS using eosohlle(ee...) 
      !          (*)to check interpolation between grid points 
      !-----------------------------------------------------------------
      ! For the red star points one computes T(E...); then one uses 
      ! eosihlle to compute E( T(E...)...). It is checked that initial
      ! and final E are the same (same for the mu's).    
      !-----------------------------------------------------------------
      ! In case of 3-flavor-EoS:
      ! eoshlle(ee...)  interpol iEoS  using  eo3.cpp  EoS3f::eos(e...)  
      ! eosohlle(ee...) computes iEoS  using  es.cpp   solveMix(*e...) 
      ! eosihlle(T...)  computes  EoS  using  es.cpp   mix_3f(*T...) 
      !-----------------------------------------------------------------
      ! Info: makeinvtable makes table of inv EoS  using  solveMix
      !-----------------------------------------------------------------
      ! In YK the EoS is used via eos->p(...), most importantly in rmn.cpp
      ! transformPV(...) which transforms from Q to {e,nb,nq,ns,vx,vy,vz}
      ! Two more calls eos->eos(...) in  hdo.cpp, but they are called after
      ! getPrimVarHCenter, which by itself calls for transformPV(...).
      !-----------------------------------------------------------------

      double precision  nmax, n0
      double precision  gammaS, B, mS, delta0, muc, volex0,
     . emax, e0, anmax, bnmax, an0, muBmax, muQmax, muSmax, Taccuracy,
     . muBmx, muQmx, muSmx,
     . aaa, bbb, temmax
      integer nxe, nxn, maxIter
      common /ceospar/  gammaS, B, mS, delta0, muc, volex0,
     . emax, e0, anmax, bnmax, an0, muBmax, muQmax, muSmax, Taccuracy,
     . muBmx, muQmx, muSmx,
     . aaa, bbb,  temmax,  nxe, nxn, maxIter
#include "aaa.h"
#include "so.h"
#include "ho.h"
      integer n(3,4)
      data n /1,0,0,  0,1,0,  0,0,1, 1,0,0 /
      character*8 ctp(3)
      character*2 ck(4,4)
      character*4 ci(5)
      character*3 cy(5)
      character*4 fomt
      data ck/'nB','nQ','nS',' E',  'nQ','nS','nB',' E'
     .       ,'nS','nB','nQ',' E' , ' E','nQ','nS','nB'/
      data ci / '  T ', '[m]B', '[m]Q', '[m]S' , '  p '/
      data cy / 'lin' , 'lin' , 'lin' , 'lin' , 'log'/
      double precision ee, n_b, n_q, n_s, pp,ttt,mub, muq, mus
      double precision eee, nb, nq, ns, ppp
      real fk(3)
      xe_max = log(emax/e0+1) 
      nmax = anmax * emax**bnmax
      n0 = an0  !* (emax/2)**bnmax 
      xn_max = log(nmax/n0+1)
      ctp(2)='htyp pbm'
      fk(1)=3
      fk(2)=1
      fk(3)=1
      
      do iii=1,5     ! y variable (1-5)
      do jjj=1,4     ! x variable, fixed variables (1-4)
      do kkk=1,2         ! fixed variables set 1 and 2  of values
      write(ifhi,'(a)')'!##############################################'
      write(ifhi,'(a)')'!        eos  T,mu ( eps, n)        '
      write(ifhi,'(a)')'!##############################################'
      if(jjj.le.3)then ! E
        x1=0.001
        x2=emax
      elseif(jjj.eq.4)then ! nB
        x1=0.00001
        x2=nmax
      endif
      nxni=nxn-1
      do k=1,3 ! different values of fixed variables
      if(k.eq.1)ctp(1)='htyp prl'
      if(k.eq.2)ctp(1)='htyp pgl'
      if(k.eq.3)ctp(1)='htyp pkl'
      if(k.eq.1)ctp(3)='htyp prf'
      if(k.eq.2)ctp(3)='htyp pge'
      if(k.eq.3)ctp(3)='htyp pkd'
      do j=1,3,2 ! different methods (method 2 obsolete -> suppressed)
        if(jjj.le.3)then
          xnb = 14.5*(k-1)/32.  * xn_max
          ww=-(kkk*2-3)*n0*(exp(xnb)-1.) 
          xtx=0.16+0.27*(k-1)
          fomt='f7.2'
        elseif(jjj.eq.4)then
          ww=e0*
     .    ( exp( (10+30*(kkk-1)+10*(k-1))/65. * xe_max ) -1. )
          xtx=0.09+0.29*(k-1)
          fomt='f8.2'
        endif
        write(ifhi
     . ,'(a/ a,2e15.6,a/ a/ a/ a/ a,f5.2,a,'//fomt//',a/ a/ a)') 
     . 'openhisto '//ctp(j)//' xmod log ymod '//cy(iii)
     . ,'xrange ',x1,x2,'  yrange auto auto'
     . ,'txt  "yaxis '//ci(iii)//'"'
     . ,'txt  "xaxis '//ck(4,jjj)//'"'
     . ,'text 0.01 0.90 "'//ck(1,jjj)//': "'
     . ,'text ',xtx,' 0.90 "',ww,'"'
     . ,'text 0.05 1.03 "'//ck(2,jjj)//'='//ck(3,jjj)//'=0 "'
     . ,'array 2'
        if(jjj.le.3)then
          nnn=1
          nxx=nxe*fk(j)
        elseif(jjj.eq.4)then
         nxx=5*nxn*fk(j)
         nnn=nxx/2+1
        endif
        do ix =nnn,nxx-1
         if(jjj.le.3)then
           ee = e0*(exp( ix*xe_max/(nxx-1) )-1.)
           n_b=n(1,jjj)*ww
           n_q=n(2,jjj)*ww
           n_s=n(3,jjj)*ww
           x=ee
         elseif(jjj.eq.4)then
          ee=ww
          xnb = -xn_max + ix*2.*xn_max/(nxx-1)
          n_b= xnb/abs(xnb+1e-20)*n0*(exp(abs(xnb))-1.) 
          n_q=0
          n_s=0
          x=n_b
         endif
         if(x.gt.x1.and.x.lt.x2)then
           ioor=0
           if(j.eq.3.and.mod(ix-nnn,7).eq.0)then
             if(ix-nnn.eq.0)
     .       write(ifmt,'(a,4i4,3x,6a,f8.3)')'eosohlle',iii,jjj,kkk,k
     .         ,ci(iii),' vs ',ck(4,jjj),'   ',ck(1,jjj),' = ',ww
             call clop(3)
             call eosohlle(ee, n_b, n_q, n_s, ttt, mub, muq, mus, pp)
             if(mub.eq.999..or.muq.eq.999..or.mus.eq.999.)ioor=1
             !~~~~~~~~~~~~~~~~ 
             if(ioor.eq.0)then ! otherwise no solution
               call eosihlle(ttt, mub, muq, mus, eee, nb, nq, ns, ppp)
              if(ioeos.eq.21.or.ioeos.eq.22)then
                err=abs(ee-eee)+abs(n_b-nb)+abs(n_q-nq)+abs(n_s-ns)
              elseif(ioeos.eq.3)then
                err=abs(ee-eee)+abs(n_b-nb)
              else
              stop'in xxEos2'
              endif
              val=abs(ee)    +abs(n_b)    +abs(n_q)    +abs(n_s)
              if(err.gt.0.005.and.err/val.gt.0.01)
     .        write(ifmt,'(a,4f10.6,3x,2f10.3,a,3(2f10.6,a))')
     .        'UNPRECISE '
     .        ,ttt, mub, muq, mus
     .        ,ee,eee,'  '
     .        ,n_b,nb,'  '
     .        ,n_q,nq,'  '
     .        ,n_s,ns,'  '
             else
              !write(ifmt,'(a,4f7.3,a,5f7.1)')
              !.        'NO SOLUTION for e,nb,nq,ns = ',ee,n_b,n_q,n_s
              !.            ,' ==> ', ttt, mub, muq, mus, pp
             endif
             !~~~~~~~~~~~~~~~~           
           elseif(j.ne.3)then
             call eoshlle (ee, n_b, n_q, n_s, ttt, mub, muq, mus, pp) 
             if(mub.eq.999..or.muq.eq.999..or.mus.eq.999.)ioor=1
           endif
           if(ioor.eq.0)then  ! otherwise no solution
             if(iii.eq.1)y=ttt
             if(iii.eq.2)y=mub
             if(iii.eq.3)y=muq
             if(iii.eq.4)y=mus
             if(iii.eq.5)y=pp
             if( (j.le.2.and.mod(ix-nnn+1-k,3).eq.0)
     .       .or.(j.eq.3.and.mod(ix-nnn,7).eq.0) )then
               if(iii.eq.5)y=max(0.0001,y)
               write(ifhi,'(2e15.6)')x,y
             endif
           endif
         endif
        enddo   
        write(ifhi,'(a)') ' endarray closehisto '
        if(k.eq.3.and.j.eq.3)write(ifhi,'(a)')'plot 0'
        if(k.ne.3.or.j.ne.3)write(ifhi,'(a)')'plot 0-'
      enddo
      enddo
      enddo
      enddo
      enddo
      end

c-----------------------------------------------------------------------
      subroutine xxEos3(x1,x2,xmub)    !   eps,p (T) comp with lattice
c-----------------------------------------------------------------------
      double precision  gammaS, B, mS, delta0, muc, volex0,
     . emax, e0, anmax, bnmax, an0, muBmax, muQmax, muSmax, Taccuracy,
     . muBmx, muQmx, muSmx,
     . aaa, bbb, temmax
      integer nxe, nxn, maxIter
      common /ceospar/  gammaS, B, mS, delta0, muc, volex0,
     . emax, e0, anmax, bnmax, an0, muBmax, muQmax, muSmax, Taccuracy,
     . muBmx, muQmx, muSmx,
     . aaa, bbb,  temmax,  nxe, nxn, maxIter
#include "aaa.h"
#include "so.h"
#include "ho.h"
      double precision ee, n_b, n_q, n_s, pp,ttt,mub, muq, mus
      real tar(500),ear(500),par(500)
      parameter (hc=0.1973,Tc=0.170/hc)
      common/latt/temp(28),press(28),epsi(28)
      character yyy*11, xxx*6, ymd*3
      nn=2
      
      if(ioeos.eq.21.or.ioeos.eq.22.or.ioeos.eq.3)then
        imax=199
        do ix=1,imax
         ttt = x1 + ix*(x2-x1)/199
         tar(ix)=ttt
         mub=xmub
         muq=0
         mus=0
         call eosihlle(ttt, mub, muq, mus, ee, n_b, n_q, n_s, pp) 
         ear(ix)=ee
         par(ix)=pp
        enddo
      else
        imax=0
        do j=2,mxEeos,4
          imax=imax+1
          ear(imax)=Eeos(j)
          tar(imax)=eost(5,1,j)
          par(imax)=eost(4,1,j)
        enddo
      endif
        
      !call Lattice(npnts)
      call Lattice2010(npnts)

      do kk=1,2
      
      if(kk.eq.1)then
        y1=0.001
        y2=200
        fac=1
        xxx='(GeV) '
        yyy='(GeV/fm^3!)'
        ymd='log'
      else
        y1=0
        y2=30
        fac=1/hc
        xxx='(1/fm)'
        yyy='/ T^4!   '
        ymd='lin'
      endif
      
      write(ifhi,'(a)')'!##############################################'
      write(ifhi,'(a)')'!        eps(T) comp with lattice       '
      write(ifhi,'(a)')'!##############################################'
      write(ifhi,'(a,2i1,a/ a,2e15.6,a/a/ a/ a/ a)') 
     . 'openhisto name eos1',kk,nn,' htyp lin xmod lin ymod '//ymd
     .,'xrange ',x1*fac,x2*fac,'  yrange 0 auto ' ! ,y1,y2
     .,'txt  "title   " '
     .,'txt  "yaxis E '//yyy//'"'
     .,'txt  "xaxis T '//xxx//'"'
     .,'array 2'
      do ix =1,imax
         x=tar(ix)
         y=ear(ix)
         fac2=1.
         if(kk.eq.2)fac2=1./(x*fac)**4
         write(ifhi,'(2e15.8)')x*fac,y*fac*fac2
      enddo   
      write(ifhi,'(a)') ' endarray closehisto '
      write(ifhi,'(a)')'plot 0-'

      write(ifhi,'(a,2i1,a/ a,2e15.6,a,2e15.6/ a)') 
     . 'openhisto name eos2',kk,nn,' htyp pnt xmod lin ymod '//ymd
     .,'xrange ',x1*fac,x2*fac,'  yrange ',y1,y2
     .,'array 2'
      do ix =1,npnts
         x = temp(ix)
         y = epsi(ix)
         fac2=1.
         if(kk.eq.2)fac2=1./(x*fac)**4
         write(ifhi,'(2e15.6)')x*fac,y*fac*fac2
      enddo   
      write(ifhi,'(a)') ' endarray closehisto '
      write(ifhi,'(a)')'plot 0'

      write(ifhi,'(a)')'!##############################################'
      write(ifhi,'(a)')'!        p(T) comp with lattice       '
      write(ifhi,'(a)')'!##############################################'
      write(ifhi,'(a,2i1,a/ a,2e15.6,a/a/ a/ a/ a)') 
     . 'openhisto name eos3',kk,nn,' htyp lin xmod lin ymod '//ymd
     .,'xrange ',x1*fac,x2*fac,'  yrange 0 auto '
     .,'txt  "title   " '
     .,'txt  "yaxis p '//yyy//'"'
     .,'txt  "xaxis T '//xxx//'"'
     .,'array 2'
      do ix =1,imax
         x=tar(ix)
         y=par(ix)
         fac2=1.
         if(kk.eq.2)fac2=1./(x*fac)**4
         write(ifhi,'(2e15.8)')x*fac,y*fac*fac2
      enddo   
      write(ifhi,'(a)') ' endarray closehisto '
      write(ifhi,'(a)')'plot 0-'

      write(ifhi,'(a,2i1,a/ a,2e15.6,a/ a)') 
     . 'openhisto name eos4',kk,nn,' htyp pnt xmod lin ymod '//ymd
     .,'xrange ',x1*fac,x2*fac,'  yrange 0 auto'
     .,'array 2'
      do ix =1,npnts
         x = temp(ix)
         y = press(ix)
         fac2=1.
         if(kk.eq.2)fac2=1./(x*fac)**4
         write(ifhi,'(2e15.6)')x*fac,y*fac*fac2
      enddo   
      write(ifhi,'(a)') ' endarray closehisto '
      write(ifhi,'(a)')'plot 0'

      enddo

      write(ifhi,'(a)')'!##############################################'
      write(ifhi,'(a)')'!        dp/depsilon       '
      write(ifhi,'(a)')'!##############################################'
      write(ifhi,'(a,a/ a,2e15.6,a/ a/ a/ a/ a)') 
     . 'openhisto name eos5',' htyp lin xmod lin ymod '//ymd
     .,'xrange ',x1*fac,x2*fac,'  yrange 0 auto'
     .,'txt  "title " '
     .,'txt  "yaxis dp/d[e] "'
     .,'txt  "xaxis T (1/fm)"'
     .,'array 2'
      do ix =2,imax
         x=tar(ix)
         y=(par(ix)-par(ix-1))/(ear(ix)-ear(ix-1))
         fac=1/hc
         write(ifhi,'(2e17.8)')x*fac,y
      enddo   
      write(ifhi,'(a)') ' endarray closehisto '
      write(ifhi,'(a)')'plot 0'
 
      end

c------------------------------------------------------------------------------
      subroutine Lattice(npnts)
c------------------------------------------------------------------------------
      parameter (hc=0.1973,Tc=0.170/hc)
      common/latt/temp(28),press(28),epsi(28)
      real t2tc(14),p2T4(14),e2T4(14)
      ! T/Tc  no units
      data (t2tc(i),i=1,14) /0.80,0.87,0.96,1.02,1.07,1.14,1.20,1.28
     *          ,1.35,1.52,1.70,1.90,2.24,2.55/
      ! p/T^4,  no units
      data (p2T4(i),i=1,14) /0.05,0.15,0.39,0.60,0.86,1.12,1.40,1.66
     *          ,1.91,2.32,2.65,2.89,3.19,3.41/
      npnts=14
      do i=1,14
       temp(i)=t2tc(i)*Tc*hc                ! in GeV
       press(i)=p2T4(i)*(t2tc(i)*Tc)**4*hc  ! in GeV/fm^3
      enddo
      do i=2,13
        f1=p2T4(i-1)*(t2tc(i-1)*Tc)**4
        f2=p2T4(i  )*(t2tc(i  )*Tc)**4
        f3=p2T4(i+1)*(t2tc(i+1)*Tc)**4
        a=(t2tc(i  )*Tc - t2tc(i-1)*Tc)
        b=(t2tc(i+1)*Tc - t2tc(i  )*Tc)
        s=(f2-f1)/a*b/(a+b) + (f3-f2)/b*a/(a+b)
        s2T3=s / (t2tc(i)*Tc)**3
        e2T4(i)=  s2T3 - p2T4(i)
        epsi(i)=e2T4(i)*(t2tc(i)*Tc)**4*hc  ! in GeV/fm^3
      enddo
      end
      
c------------------------------------------------------------------------------
      subroutine Lattice2010(npnts)
c------------------------------------------------------------------------------
      !S. Borsanyi et al. arXiv:1007.2580
      !------------------------------------
      common/latt/temp(28),press(28),epsi(28)
      real tmp(28),p2T4(28),d2T4(28)
      parameter (hc=0.1973)
       !~~~~~~~~~~ T (MeV)
      data tmp/
     *100,115,129,134,139,143,147,152,158,162,166,170,175,185,200,215
     *,228,250,275,299,330,366,400,450,500,600,800,1000/
      !~~~~~~~~~~ pressure / T4 (no units)
      data p2T4/
     *0.16,0.23,0.30,0.34,0.40,0.45,0.51,0.58,0.68,0.76,0.85,0.94
     *,1.05,1.27,1.57,1.85,2.06,2.37,2.66,2.87,3.08,3.29,3.45
     *,3.62,3.73,3.90,4.09,4.19/      
      !~~~~~~~~~~~ (epsilon-3*pressure) / T4 (no units)
      data d2T4/
     *0.41,0.52,0.95,1.25,1.60,1.89,2.20,2.60,3.06,3.33,3.55,3.72
     *,3.87,3.98,3.92,3.75,3.55,3.18,2.78,2.46,2.13,1.82,1.59,1.31
     *,1.08,0.77,0.50,0.45/ 
      npnts=28
      do i=1,28
       temp(i)=  tmp(i)/1000.               ! in GeV
       press(i)=p2T4(i)*(temp(i)/hc)**4*hc  ! in GeV/fm^3
       epsi(i)= d2T4(i)*(temp(i)/hc)**4*hc  ! in GeV/fm^3
       epsi(i)=epsi(i)+3*press(i)
      enddo
      end
     
c-----------------------------------------------------------------------
      subroutine eostest2
c-----------------------------------------------------------------------
      double precision  nmax, n0
      double precision  gammaS, B, mS, delta0, muc, volex0,
     . emax, e0, anmax, bnmax, an0, muBmax, muQmax, muSmax, Taccuracy,
     . muBmx, muQmx, muSmx,
     . aaa, bbb, temmax
      integer nxe, nxn, maxIter
      common /ceospar/  gammaS, B, mS, delta0, muc, volex0,
     . emax, e0, anmax, bnmax, an0, muBmax, muQmax, muSmax, Taccuracy,
     . muBmx, muQmx, muSmx,
     . aaa, bbb,  temmax,  nxe, nxn, maxIter
      double precision xe,xnb,xnq,xns,e,nb,nq,ns  
      double precision xe_max,xn_max
      double precision ttt, mub, muq, mus, pp
      double precision tt, mu_b, mu_q, mu_s, p_p
      xe_max = log(emax/e0+1) 
      nmax = anmax * emax**bnmax
      n0 = an0  !* (emax/2)**bnmax 
      xn_max = log(nmax/n0+1)
      tol=0.006
      do ixe =0,nxe-1,5
       xe = ixe*xe_max/(nxe-1)
       e = e0*(exp(xe)-1.)
         write(ifmt,*)'+++++ xe e = ',xe,e
      do ixnb=0,nxn-1,5
      do ixnq=0,nxn-1,5
      do ixns=0,nxn-1,5
       xnb = -xn_max + ixnb*2.*xn_max/(nxn-1) 
       xnq = -xn_max + ixnq*2.*xn_max/(nxn-1) 
       xns = -xn_max + ixns*2.*xn_max/(nxn-1) 
       nb = xnb/abs(xnb+1e-20)*n0*(exp(abs(xnb))-1.) 
       nq = xnq/abs(xnq+1e-20)*n0*(exp(abs(xnq))-1.) 
       ns = xns/abs(xns+1e-20)*n0*(exp(abs(xns))-1.) 
       call eosohlle(e, nb, nq, ns, ttt, mub, muq, mus, pp) 
       call eoshlle (e, nb, nq, ns, tt, mu_b, mu_q, mu_s, p_p) 
       err=abs(tt-ttt)+abs(mu_b-mub)+abs(mu_q-muq)+abs(mu_s-mus)
       val=abs(tt)    +abs(mu_b)    +abs(mu_q)    +abs(mu_s)
       if(err.gt.tol.and.err/val.gt.tol)
     .     write(ifmt,'(a,4f9.5,a,4(2f9.5,a))')'ERROR ',e,nb,nq,ns,' ',
     .     tt,ttt,' '
     .    ,mu_b,mub,' '
     .    ,mu_q,muq,' '
     .    ,mu_s,mus,' '
       enddo  
       enddo  
       enddo  
       enddo  
       
       end

c-----------------------------------------------------------------------
      subroutine eostest4
c-----------------------------------------------------------------------
#include "aaa.h"
      double precision tt, mu_b, mu_q, mu_s, ee, n_b, n_q, n_s, pp


      do i=1,5
        tt=0.02+i*0.01
      print*,'------------------------------------------------------'
      do j=0,4
      mu_b=j
      mu_q=0
      mu_s=0
        call eosihlle(tt, mu_b, mu_q, mu_s, ee, n_b, n_q, n_s, pp)
        write(ifmt,'(a,3(2f9.5,3x))')' t e mu n:',tt,ee
     .   ,mu_b, n_b,    n_q, n_s
      enddo
      enddo
      end
            
c-----------------------------------------------------------------------
      subroutine eostest1(eos_hlle)
c-----------------------------------------------------------------------
#include "aaa.h"
      double precision tt, mu_b, mu_q, mu_s, ee, n_b, n_q, n_s, pp
     .  ,ttt,mub, muq, mus
      tol=0.01

      do kkk=1,4
        if(kkk.eq.1)then      
        mu_b=0.4
        mu_q=0.1
        mu_s=0.2
        elseif(kkk.eq.2)then
        mu_b=0.4
        mu_q=0
        mu_s=0
        elseif(kkk.eq.3)then
        mu_b=0
        mu_q=0.1
        mu_s=0
        elseif(kkk.eq.4)then
        mu_b=0
        mu_q=0
        mu_s=0.2
        endif
        tt=0.050+0.15/20.*i
        write(ifmt,'(a,3f9.5,a)')'test invtab',mu_b,mu_q,mu_s
        do i=1,20
         tt=0.050+0.15/20.*i
         call eosihlle(tt, mu_b, mu_q, mu_s, ee, n_b, n_q, n_s, pp)
         call eos_hlle(ee, n_b, n_q, n_s, ttt, mub, muq, mus, pp) 
         write(ifmt,'(2(4f7.3,3x))')
     .    tt, mu_b, mu_q, mu_s,ee, n_b, n_q, n_s
         err=abs(tt-ttt)+abs(mu_b-mub)+abs(mu_q-muq)+abs(mu_s-mus)
         val=abs(tt)    +abs(mu_b)    +abs(mu_q)    +abs(mu_s)
         if(err.gt.tol.and.err/val.gt.tol)
     .     write(ifmt,'(a,4f7.3,3x,4(2f7.3,a))')'ERROR '
     .    ,ee, n_b, n_q, n_s
     .    ,tt,ttt,'  '
     .    ,mu_b,mub,'  '
     .    ,mu_q,muq,'  '
     .    ,mu_s,mus,'  '
        enddo   
       enddo
       stop  
       end 
c-----------------------------------------------------------------------
      subroutine eostest5
c-----------------------------------------------------------------------
#include "aaa.h"
      double precision tt, mu_b, mu_q, mu_s, ee, n_b, n_q, n_s, p_p
     .  ,ttt,mub, muq, mus,pp    ,tttx,mubx, muqx, musx,ppx
        mu_b=0
        mu_q=0
        mu_s=0
        tt=0.5
        call eosihlle(tt, mu_b, mu_q, mu_s, ee, n_b, n_q, n_s, p_p)
        call eosohlle(ee, n_b, n_q, n_s, ttt, mub, muq, mus, pp) 
        call eoshlle(ee, n_b, n_q, n_s, tttx, mubx, muqx, musx, ppx) 
        write(ifmt,'(a,4f5.1,3x,4(3f7.3,a))')'TEST5 '
     .    ,ee, n_b, n_q, n_s
     .    ,tt,ttt,tttx,'  '
     .    ,mu_b,mub,mubx,'  '
     .    ,mu_q,muq,muqx,'  '
     .    ,mu_s,mus,musx,'  '
       end
c-----------------------------------------------------------------------
      subroutine eostest6
c-----------------------------------------------------------------------
#include "aaa.h"
      real a(4,2)
      data ((a(i,j),i=1,4),j=1,2)/
     . 100., 35.36 , 0, 0,
     .  100., 0 , 0, 0/
        double precision ee, n_b, n_q, n_s
     .  ,ttt,mub, muq, mus,pp    ,tttx,mubx, muqx, musx,ppx
        do j=2,2
        ee =a(1,j)
        n_b=a(2,j)          
        n_q=a(3,j)
        n_s=a(4,j)
        call eosohlle(ee, n_b, n_q, n_s, ttt, mub, muq, mus, pp) 
        call eoshlle( ee, n_b, n_q, n_s, tttx, mubx, muqx, musx, ppx) 
        write(ifmt,'(a,4f7.2,3x,5(2f9.3,a))')'TEST6 '
     .    ,ee, n_b, n_q, n_s
     .    ,ttt,tttx,'  '
     .    ,mub,mubx,'  '
     .    ,muq,muqx,'  '
     .    ,mus,musx,'  '
     .    ,pp,ppx,'  '
       enddo
       end

c-----------------------------------------------------------------------
      subroutine eostest7
c-----------------------------------------------------------------------
#include "aaa.h"
      double precision ee, n_b, n_q, n_s
     . ,ttt,mub, muq, mus,pp 
     .  ,eex, n_bx, n_qx, n_sx, p_p
      do i=1,1
       ee  = 0.164325E+00
       n_b = 35.36
       n_q = 0.
       n_s = 0.
       call eosohlle(ee, n_b, n_q, n_s, ttt, mub, muq, mus, pp) 
       call eosihlle(ttt, mub, muq, mus, eex, n_bx, n_qx, n_sx, p_p)
       write(ifmt,'(a,5(3f8.2,3x))')  'TEST7 '
     .   ,ee,eex,ttt,   n_b,n_bx,mub,    n_q, n_qx,muq,  n_s,n_sx,mus
     .   ,pp,p_p
      enddo
      end

c-----------------------------------------------------------------------
      subroutine eostest8
c-----------------------------------------------------------------------
#include "aaa.h"
      double precision
     .  ttt,mub, muq, mus
     .   ,eex, n_bx, n_qx, n_sx, p_p
       do i=1,20
       ttt = 0.0508301
       mub = 2*i
       muq =  -0.218604
       mus = 1.03501 
        call eosihlle(ttt, mub, muq, mus, eex, n_bx, n_qx, n_sx, p_p)
        write(ifmt,'(a,2(4f12.4,a))')  'TEST8 ',
     .     ttt, mub, muq, mus,'  ==>  ', eex, n_bx, n_qx, n_sx
       enddo
       end

