C
C  This file is part of EPOS4
C  Copyright (C) 2022 research institutions and authors (See CREDITS file)
C  This file is distributed under the terms of the GNU General Public License version 3 or later
C  (See COPYING file for the text of the licence)
C

c-------------------------------------------------------------------------
c               hydro output
c-------------------------------------------------------------------------

      integer iBspes, iHoEpsilon, iHoEpsilonEtaY, iHoEpsilonEtas6, 
     . iHoEpsilonEtas8, iHoEpsilonTau6, iHoEpsilonTau8, iHoRadVel, 
     . iHoRadVelTau6, iHoTanVel, iHoTangVelTau6, iHoTemperatureEtas6, 
     . iHoTemperatureEtas8, iHyAverages, iHyBarmuEta, iHyBaryon, 
     . iHyBaryonEta, iHyEmpty, iHyEntropy, iHyEntropyEta, iHyEos, 
     . iHyEpsilon, iHyEpsilon2, iHyEpsilonEta, iHyFoBarmu, iHyFoEpsilon, 
     . iHyFoRadVelocity, iHyFoRadius, iHyFoTangVelocity, iHyFoVol, 
     . iHyLongVelocity, iHyRadVelocity, iHyTemperature, iQspes, iSource, 
     . iSspes, ienvar, ifaahlle, ifathlle, ifazhlle, imakefzo, ireadfzo, 
     . ispes, istat, jspes, kfrout, klax, kspes, mchem, metahxx, metahy, 
     . mspes, mtemp, mxBeos, mxEeos, mxbimp, mxsurf, mxxBeos, mxxEeos, 
     . nbimp, nbimpx, ncenthxx, ncenthy, netahxx, nfrout, nlag, nphahy, 
     . nphihxx, nphihy, nptrhy, nradhxx, nradhy, nspes, ntauhec, 
     . ntauhxx, ntauhy, nxhxx, nxhy, nyhxx, nyhy, nzhy, nzhyxx

      real Beos, Eeos, aspes, baraa, barij, barzz, bimpar, bimparx, 
     . centhy, chepoB, chepoQ, chepoS, chot, cohi, dtauhy, eccpav, 
     . eccxav, efrout, emuaa, emuzz, eost, epsaa, epscrit, epsfin, 
     . epsij, etahy, fdtau, ffstat, fofac, fofai, ggstat, gspes, oBeos, 
     . oEeos, phihy, ptrmx, radaa, radhy, rmaxhy, sigij, suraa, 
     . taumax, taumaxhy, tauminhy, temax, temfo, temin, tfo, tfrout, 
     . uBeos, uEeos, velaa, velzz, vlmaa, vtraav, wlag, xlag, xmaxhy, 
     . xminhy, ymaxhy, yminhy, zetahy, zmaxhy, zminhy

#if __TP__
      parameter (nxhxx=1,nyhxx=1,netahxx=1,ntauhxx=1) !points, not bins!!
      parameter(metahxx=1)
      parameter (ncenthxx=1)
      parameter (nphihxx=1, nradhxx=1)
#else
      parameter (nxhxx=51,nyhxx=51,netahxx=54,ntauhxx=120) !points, not bins!!
      parameter(metahxx=netahxx/2)
      parameter (ncenthxx=20)
      parameter (nphihxx= 101, nradhxx=201)
#endif

      common/hocoxyeta/nxhy,nyhy, nzhy,ntauhy
      common/hocotau/tauminhy,taumaxhy
      common/hocobound/xminhy,xmaxhy,yminhy,ymaxhy,zminhy,zmaxhy
      common/hocotau2/dtauhy
      common/cifahlle/ifaahlle,ifathlle,ifazhlle
      common/cepsfin/epsfin

      ! tau = tauminhy+(ntau-1)*(taumaxhy-tauminhy)/(ntauhy-1)
      ! eta = zminhy  +(nz-1)  *(zmaxhy-zminhy)    /(nzhy-1)
      ! x   = xminhy  +(nx-1)  *(xmaxhy-xminhy)    /(nxhy-1)
      ! y   = yminhy  +(ny-1)  *(ymaxhy-yminhy)    /(nyhy-1)

      common/hoco10/vtraav(netahxx,ntauhxx)
     *          ,eccxav(netahxx,ntauhxx),eccpav(netahxx,ntauhxx)

      common/hoco1/ncenthy,nradhy,nphihy,nphahy,nptrhy,ptrmx
      
      common/hoco3/ntauhec(netahxx,nphihxx)
  
      common/hoco4/centhy(0:ncenthxx),etahy(netahxx)
     *             ,phihy(nphihxx),radhy(nradhxx)
      common/hoco5/zetahy(1-netahxx:netahxx-1)

      common/hoco7a/metahy     

      common/hoco8/taumax,nzhyxx(2),epscrit,rmaxhy

      common/hoco9/epsij(netahxx,ntauhxx)
     *            ,barij(3,netahxx,ntauhxx)
     *            ,sigij(netahxx,ntauhxx)

c------------------------------------------------------------------------

      parameter(mxsurf=4)
      common/hoco7/
     *   radaa(mxsurf,    netahxx,ntauhxx,nphihxx)
     *  ,velaa(mxsurf, 3 ,netahxx,ntauhxx,nphihxx)
     *  ,vlmaa(mxsurf,    netahxx,ntauhxx,nphihxx)
     *  ,epsaa(mxsurf,    netahxx,ntauhxx,nphihxx)
     *  ,baraa(mxsurf, 3 ,netahxx,ntauhxx,nphihxx)
     *  ,emuaa(mxsurf,10 ,netahxx,ntauhxx,nphihxx)
      common/hoco7b/emuzz(10 ,netahxx,nradhxx,nphihxx)
      common/hoco7c/barzz( 3 ,netahxx,nradhxx,nphihxx)
      common/hoco7d/velzz( 3 ,netahxx,nradhxx,nphihxx)
      common/csurfelem/
     *   suraa(mxsurf, 4 ,netahxx,ntauhxx,nphihxx)

c------------------------------------------------------------------------

      common/ctfo/tfo /cfzo/ireadfzo,imakefzo

      common/copt/istat,ienvar 
      parameter (mspes=500,mchem=21,mtemp=17,temin=0.010,temax=0.170)
      common/cspes11/nspes,ispes(mspes),aspes(2,mspes),gspes(mspes)
      common/cspes12/kspes,jspes(mspes),iBspes(mspes)
      common/cspes16/iQspes(mspes),iSspes(mspes)
      common/cspes13/chepoB(mchem),chot(mtemp,mchem,mchem,mchem,mspes)
      common/cspes17/cohi(mspes)
      common/cspes18/temfo(mtemp)
      common/cspes15/chepoQ(mchem),chepoS(mchem)
      common/cspes14/ffstat(mtemp,2,mchem,mchem,mchem,0:mspes+2)
      common/cspes14/ggstat(0:mspes+2)
      parameter (nlag=15)
      common /clag/xlag(nlag),wlag(nlag)
      parameter (klax=5)    

      parameter (mxbimp=500 )
      common/hoco12/bimpar(2,mxbimp),nbimp
     *             ,bimparx(2,mxbimp),nbimpx

c------------------------------------------------------------------------

      common/hoco13/tfrout,fofac,fofai,kfrout,efrout,nfrout
      
c------------------------------------------------------------------------

#if __TP__
      parameter (mxxBeos=1,mxxEeos=1)
#else
      parameter (mxxBeos=41,mxxEeos=2001)
#endif
      common/ceost/ eost(7,mxxBeos,mxxEeos),Beos(mxxBeos),Eeos(mxxEeos)
      common/ceost2/uEeos,oEeos,uBeos,oBeos,mxBeos,mxEeos

c------------------------------------------------------------------------
      common/xhoco/
     .  iHyEpsilon,    iHyEpsilon2,     iHyEntropy,    iHyTemperature
     . ,iHyRadVelocity,iHyLongVelocity, iHyAverages,   iHyBaryon
     . ,iHyFoVol,      iHyFoRadius,   iHyFoRadVelocity,iHyFoTangVelocity
     . ,iHyFoBarmu,    iHyEpsilonEta, iHyBaryonEta,    iHyEntropyEta
     . ,iHyEos,        iHyEmpty,      iHyBarmuEta,     iHoTanVel
     . ,iHoEpsilon,    iHoEpsilonEtaY,iHoRadVelTau6,   iHoRadVel
     . ,iHoEpsilonEtas8,iHoEpsilonEtas6,iHoEpsilonTau8,iHoEpsilonTau6
     . ,iHoTemperatureEtas8,iHoTemperatureEtas6,iHoTangVelTau6
     . ,iHyFoEpsilon
     
      common/xplots/
     . iSource, fdtau

