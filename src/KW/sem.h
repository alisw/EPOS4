C
C  This file is part of EPOS4
C  Copyright (C) 2022 research institutions and authors (See CREDITS file)
C  This file is distributed under the terms of the GNU General Public License version 3 or later
C  (See COPYING file for the text of the licence)
C

      integer mmrys
      parameter (mmrys=10)   !memory saving factor
      integer klasmax
      parameter(klasmax=12)
      real coekaf,coelaf
      common/ccoekaf/coekaf,coelaf
      integer ibin1mx,ibin2mx
      parameter (ibin1mx=20/mmrys,ibin2mx=40/mmrys)
      integer icdp,icdt
      common /psar1/icdp,icdt
      real  alfp,ffrr,delh,r3pom,adskap,q2sft
      common /psar3/alfp,ffrr,delh,r3pom,adskap,q2sft
      real qcdlam,q2min,q2ini,q2fin,pt2cut,
     *betpom,glusea,alfe, q2dmin,qcdlambda(3:5)
      real facnof,noflit,nofeff
      common/psar6/noflit,nofeff,facnof
      integer naflav,nbflav
      common /psar5/  qcdlam,q2min,q2ini,q2fin,pt2cut,
     *betpom(2),glusea,naflav,alfe, q2dmin,qcdlambda,nbflav
      integer  nuflav
      common /psar5b/ nuflav
      real  stmass ,amhadr,qcmass,qbmass
      common /psar8/  stmass ,amhadr(8),qcmass,qbmass
      real ydmin,ydmax,qdmin,qdmax,themin,themax,elomin
      common /psar12/ ydmin,ydmax,qdmin,qdmax,themin,themax,elomin
      real factk,factgam,ptphot,factjpsi,factsat,fracsat
      common /psar14/ factk,factgam,ptphot,factjpsi,factsat,fracsat
      real factbgg,factbqq,alpsft,alpsat,factsft
      common /psar14b/ factbgg,factbqq,alpsft,alpsat,factsft
      real fInitialHF
      common /psar14c/ fInitialHF
      double precision  cYfus,cYscr,cYdif,cXfus,cXdif,cXscr
      common /psar16/ cYfus,cYscr,cYdif,cXfus,cXdif,cXscr
      real dels,gamsoft
      common/cgamsoft/dels(2),gamsoft(2)
      double precision fptnxx
      integer kzz1max,kzz2max,kxp1max,kxp2max
      real  epsi,zzmaxa,zzmina,xpmina,xpmaxa,zz2min,zzmed1,zzmed2,ccc
      common/tabfptnxx/fptnxx(5,0:2,0:2,ibin1mx,ibin1mx,ibin2mx,ibin2mx)
     *,kzz1max,kzz2max,kxp1max,kxp2max,epsi,zzmaxa
     *,zzmina,xpmina,xpmaxa,zz2min(5,0:2,0:2,ibin1mx,ibin1mx,ibin2mx)
     *,zzmed1,zzmed2(5,0:2,0:2,ibin1mx,ibin1mx,ibin2mx)
     *,ccc(5,0:2,0:2,ibin1mx,ibin1mx,ibin2mx,ibin2mx)
      real testinterxp1,testinterxp2
      common /testinterpol/
     *testinterxp1(5,0:2,0:2,ibin1mx,ibin1mx,ibin2mx,ibin2mx),
     *testinterxp2(5,0:2,0:2,ibin1mx,ibin1mx,ibin2mx,ibin2mx)
       integer medium
       common /medium_modif/medium
       integer kzz1max_s,kzz2max_s,kxp1max_s,kxp2max_s
       real  epsi_s,fptnxx_s,zzmaxa_s,zzmina_s,xpmina_s,xpmaxa_s
     *,zz2min_s,zzmed1_s,zzmed2_s,ccc_s
       common/tabf_s/kzz1max_s,kzz2max_s,kxp1max_s,kxp2max_s,epsi_s,
     *fptnxx_s(5,0:2,0:2,ibin1mx,ibin1mx,ibin2mx,ibin2mx),zzmaxa_s
     *,zzmina_s,xpmina_s,xpmaxa_s
     *,zz2min_s(5,0:2,0:2,ibin1mx,ibin1mx,ibin2mx)
     *,zzmed1_s,zzmed2_s(5,0:2,0:2,ibin1mx,ibin1mx,ibin2mx)
     *,ccc_s(5,0:2,0:2,ibin1mx,ibin1mx,ibin2mx,ibin2mx)
      real testinterxp1_s,testinterxp2_s
      common /testinterpol_s/
     *testinterxp1_s(5,0:2,0:2,ibin1mx,ibin1mx,ibin2mx,ibin2mx),
     *testinterxp2_s(5,0:2,0:2,ibin1mx,ibin1mx,ibin2mx,ibin2mx)
