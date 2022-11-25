C
C  This file is part of EPOS4
C  Copyright (C) 2022 research institutions and authors (See CREDITS file)
C  This file is distributed under the terms of the GNU General Public License version 3 or later
C  (See COPYING file for the text of the licence)
C

      integer maxdataDf,nbpf
      parameter(maxdataDf=100,nbpf=10)
      double precision xminDf,xmaxDf,xfitmin,xshmin,xggfit
     &                ,datafitD
      integer         nmcxDf,nmcbDf,nptf
      real            smaxDf,bmaxDf,sfshlim
      common /parvar/ xminDf,xmaxDf,xfitmin,xshmin,xggfit
     &               ,nmcxDf,nmcbDf,nptf,smaxDf,bmaxDf,sfshlim
      real            parDf,fparDf,betac,betae
      integer         numparDf,numdataDf,numminDf
      common /fitpar/ datafitD(maxdataDF,3),parDf(nbpf,4),
     &                fparDf,betac,betae,numparDf,numdataDf,numminDf
      real            epscrw,epscrp,egyscr,iscreeni,epscrs,epscrx,epscrh
     &,epscrg,zbrads,epscrd,bglaubx,b2xscr,fegypp,zbcutx,epscrxi,epscrb
      common /epspar/epscrw,epscrp,egyscr,iscreeni,epscrs,epscrx,epscrh
     &,epscrg,zbrads,epscrd,bglaubx,b2xscr,fegypp,zbcutx,epscrxi,epscrb
