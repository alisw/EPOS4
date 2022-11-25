C
C  This file is part of EPOS4
C  Copyright (C) 2022 research institutions and authors (See CREDITS file)
C  This file is distributed under the terms of the GNU General Public License version 3 or later
C  (See COPYING file for the text of the licence)
C

      integer mxpsar3 , mxpsar1
      parameter (mxpsar3=20)
      integer jleol
      common /ceol/ jleol(mxpsar3,2)
      integer jlklas
      common /cjlklas/ jlklas(mxpsar3,3)
      integer maxq2mn
      parameter (maxq2mn=16)
      parameter (mxpsar1=20*mxpsar3)
      real cstot, csord, csbor, cstotzero,csborzer,cstotr,csordr
      real csqtotzero, csqborzer, csbor0
      common /psar19/  cstot(20,20,20,mxpsar3,2) ! ( q1, q2, e, flav, 1)
      integer kenwidth
      parameter (kenwidth=2)  ! could be 1,2,4,5,10,20 
      common /psar19b/ cstotr(20,20,2,kenwidth)
      common /psar20/  csord(         20,20,2*mxpsar1)
      common /psar20b/ csordr(  20,20,20,2 )
      common /psar21/csbor( 20, mxpsar1, 2),csbor0( 20, 20, mxpsar1)
      common /psar22/  cstotzero(     20,     mxpsar3)
     .                ,csborzer(      20,     mxpsar3)
      common /psarq22/ csqtotzero(maxq2mn,maxq2mn,20,mxpsar3)
     .                 ,csqborzer(maxq2mn,maxq2mn,20,mxpsar3)
      real csjet
      common/tabcsjt/ csjet(-1:2,2,20,1,1,mxpsar3,1000)
c      common/tabcsjt/ csjet(-1:2,2,20,1,1,  6    ,1000)

      integer    maxzzvex
      parameter (maxzzvex=1)
      integer     mxzzvex
      character*2 zzvexch
      common/cvex/mxzzvex,zzvexch(maxzzvex)

      integer nq2mnfix,maxq2mx
      real       q2mnval
      character*3 q2mnch
      common/csemev/ q2mnval(maxq2mn),q2mnch(maxq2mn)
      common/cnq2mnfix/ nq2mnfix,maxq2mx
      integer mnclpt,mxclpt
      parameter (mnclpt=6,mxclpt=6)
      real fhxss,fhxgg,fhxqg,fhxgq,fhxqq
      common /psarq4/
     *   fhxss(maxq2mn,maxq2mn,maxzzvex,11,100,10,mxclpt-mnclpt+1)
     *  ,fhxgg(maxq2mn,maxq2mn,maxzzvex,11,100,10,mxclpt-mnclpt+1)
     *  ,fhxqg(maxq2mn,maxq2mn,maxzzvex,11,100,10,mxclpt-mnclpt+1)
     *  ,fhxgq(maxq2mn,maxq2mn,maxzzvex,11,100,10,mxclpt-mnclpt+1)
     *  ,fhxqq(maxq2mn,maxq2mn,maxzzvex,11,100,mxclpt-mnclpt+1)

