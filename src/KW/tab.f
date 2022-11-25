C
C  This file is part of EPOS4
C  Copyright (C) 2022 research institutions and authors (See CREDITS file)
C  This file is distributed under the terms of the GNU General Public License version 3 or later
C  (See COPYING file for the text of the licence)
C

c=======================================================================
! Technical remark concerning the saturation scales (KW, 25.5.2013):
!
! We employ (with very exceptions for very local use) the following convention:
! the saturation scale is always named q2.min, where the point "." stands for
! a single character (so a Regex search with "q2.min" will find all variables)
!
! To accomodate two saturation scales (proj and targ side) we use
! two element arrays, like q2cmin(2)
!
! q2cmin(2) refer to "current" scales (incl x)
! q2pmin(2) refer to "pair" scales for tabulation (max Q2s for x=1)
! q2zmin refers to the minimum scale used for binning for tables
! q2nmin is the USED minimum scale (low x or large b) (tables do not depend on this)
!
! look for q2.min (Regex search) plus "??????" to locate questionable
! definitions, still to check  


c         ------------------------------------------------
c         since EPOS3086 we have om5 = om1
c                  (in earlier versions   om5 = 0.5 * om1)
c         change in : testOmExIpo
c         ------------------------------------------------
c=======================================================================




c####################################################################################
c####################################################################################
c####################################################################################
c####################################################################################
c#######
c#######                       flavor pair coding
c#######
c####################################################################################
c####################################################################################
c####################################################################################
c####################################################################################


c=======================================================================

c KW1808 August 2018 :
c KW1808    creation of tables csbor, csord, cstot, csborzer, cstotzero 
c KW1808    and the corresponding interpolation fuctions 
c KW1808    extended to 5 flavors

c-----------------------------------------------------------------------
c KW1808 Reminder: we use "basic functions" to make predefined tables, 
c           (then make interpolated tables (for given q2) in some cases
c             marked as ... -> ... -> ... in the following ) 
c           then  define "interpolation functions"
c ----------------------------------------------------------------------
c basic functions            tables              interpolation functions
c used to make tables                            based on the tables
c-----------------------------------------------------------------------
c   psborn,psjet            cstot                              psjti  
c   psborn,psjet1           csord                              psjti1 
c   psborn                  csbor                              psboi
c   psborn           csborzer -> csqborzer -> csborzer        sr psjti0
c   (uses psjti)     cstotzero -> csqtotzero -> cstotzero     sr psjti0

c All basic functions have two arguments referring to the parton type at 
c   current and opposite side, say j1 , l1 .
c All tables use one integer, say  jl to tabulate the corresponding 
c cross section values.
c All interpolation functions have two arguments referring to the parton 
c  type at  current and opposite side, say j2 , l2. 
c Subroutine initEoL defines the relations between the indices.
c=======================================================================


c-----------------------------------------------------------------------
      subroutine initEoL
c----------------------------------------------------------------------- 
c KW1808c end of ladder parton pair classes
c array jleol will be used as
c                        j1 = jleol(jl,1)
c                        l1 = jleol(jl,2)
c----------------------------------------------------------------------- 
#include "tab.h"
      integer ieol(mxpsar3,2)
      data ((ieol(jl,k),k=1,2),jl=1,mxpsar3) /
      !~~~~~~~~~~~~~!~~~~~!~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
      !    j1  l1   !  jl ! j2 l2
      !~~~~~~~~~~~~~!~~~~~!~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
                          !in the following, the symbols mean 
                          ! g ....... gluon (0)
                          ! q,q',a .. light quark or antiquark (+-1,2,3)
                          ! c ....... charm quark or antiquark   (+-4)
                          !  b ....... bottom quark or antiquark   (+-5)
      !~~~~~~~~~~~~~!~~~~~!~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
     .      0 , 0   !  1  !  g g    two gluons
     .   ,  0 , 1   !  2  !  g q    one gluon, one light (anti)quark
     .   ,  1 , 0   !  3  !  q g
     .   ,  1 , 1   !  4  !  q q    two light (anti)quarks of same flavor
     .   ,  2 ,-2   !  5  !  q a    light quark and antiquark of same flavor
     .   ,  2 , 1   !  6  !  q q'   two light (anti)quarks of different flavor
     .   ,  0 , 3   !  7  !  g c
     .   ,  3 , 0   !  8  !  c g
     .   ,  0 , 4   !  9  !  g b
     .   ,  4 , 0   ! 10  !  b g
     .   ,  1 , 3   ! 11  !  q c
     .   ,  3 , 1   ! 12  !  c q
     .   ,  1 , 4   ! 13  !  q b
     .   ,  4 , 1   ! 14  !  b q
     .   ,  3 , 3   ! 15  !  c c    (4,4 or -4,-4                  
     .   ,  3 ,-3   ! 16  !  c ac   (4,-4 or -4 4)                 
     .   ,  4 , 4   ! 17  !  b b    (5,5 or -5,-5)                 
     .   ,  4 ,-4   ! 18  !  b ab   (5,-5 or -5,5)                 
     .   ,  3 , 4   ! 19  !  c b    (+-4,+-5)                       
     .   ,  4 , 3   ! 20  !  b c    (+-5,+-4)                       
     .   /
      integer iklas(mxpsar3,3) !---->  klas defined in subroutine getKlasString
      data ((iklas(jl,k),k=1,3),jl=1,mxpsar3) /
     .    1,4,5,  1,0,0,  1,0,0,  1,0,0,  1,4,5,  
     .    1,0,0,  2,0,0,  2,0,0,  3,0,0,  3,0,0,
     .    2,0,0,  2,0,0,  3,0,0,  3,0,0,  8,0,0,
     .    6,8,9,  11,0,0, 7,10,11, 12,0,0,  12,0,0 /
      do jl=1,mxpsar3
        do k=1,2
          jleol(jl,k)=ieol(jl,k)
        enddo
      enddo
      do jl=1,mxpsar3
        do k=1,3
          jlklas(jl,k)=iklas(jl,k)
        enddo
      enddo
      end

c-----------------------------------------------------------------------
c KW1808 new subroutine to manage End Of Ladder (eol) Parton Indices
c Maps tabulation indices jl used to make tables csord etc
c   to eol parton indices j1,l1 used as arguments 
c        in psjetj, psjet1, psjeti, psjet 
c and maps eol parton indices j2,l2 used as arguments 
c        in interpolation functions psjti,psjti1,psboi,psjti0     
c   to tabulation indices jl.
c            j1,2 = eol parton type at current side
c            l1,2 = eol parton type at opposite side
c The indices j2,l2 take all possible values of j2: -5 to 5 and l2: -5 to 5, 
c For j1,l1 only the selection of paires in the table above are defined,
c    each one represents a GROUP of parton pairs 
c-----------------------------------------------------------------------
c The mapping jl --> j1 , l1 is done via
c        call eolpi( jl , j1 , l1 ) 
c The mapping  j2 , l2 --> jl is done via
c        call eolpib( j2 , l2 , jl ) 
c The mapping  j1 , l1 --> j2 , l2 is done via
c        call eolpic( j1 , l1 , j2 , l2 ) 
c (the functions are defined below)
c-----------------------------------------------------------------------

c-----------------------------------------------------------------------
      subroutine eolpi( jl , j1 , l1 )  !KW1808c  mapping  jl -> j1 l1  
c-----------------------------------------------------------------------
      implicit none
      integer jl,j1,l1
#include "tab.h"
      j1 = jleol(jl,1)
      l1 = jleol(jl,2)
      end

c-----------------------------------------------------------------------
      subroutine eolpib( j2 , l2 , jl )  !KW1808c mapping  j2 l2 -> jl    
c-----------------------------------------------------------------------
      implicit none
      integer j2,l2,jl,j2a,l2a
      j2a=abs(j2)
      l2a=abs(l2)
      if(j2a.le.3.and.l2a.le.3)then !light
          if(    j2.eq.0.and.l2.eq.0)then  ! g g
            jl=1
          elseif(j2.eq.0.and.l2.ne.0)then  ! g q
            jl=2
          elseif(j2.ne.0.and.l2.eq.0)then  ! q g
            jl=3
          elseif(j2.eq.l2)then             ! q q
            jl=4
          elseif(j2.eq.-l2)then            ! q a
            jl=5
          else                             ! q q' 
            jl=6
          endif
      elseif(j2.eq.0.or.l2.eq.0)then !gX Xg
          if(l2a.eq.4)then                 ! g c 
            jl=7
          elseif(j2a.eq.4)then             ! c g 
            jl=8
          elseif(l2a.eq.5)then             ! g b 
            jl=9
          elseif(j2a.eq.5)then             ! b g 
            jl=10
           endif
      elseif(j2a.le.3.or.l2a.le.3)then !qX Xq
          if(l2a.eq.4)then                 ! q c 
            jl=11
          elseif(j2a.eq.4)then             ! c q 
            jl=12
          elseif(l2a.eq.5)then             ! q b 
            jl=13
          elseif(j2a.eq.5)then             ! b q 
            jl=14
          endif
        elseif(j2a.eq.4.and.l2a.eq.4)then 
          if(j2.eq.l2)then                 ! c c 
            jl=15
          else                             ! c ac 
            jl=16
          endif
      elseif(j2a.eq.5.and.l2a.eq.5)then 
          if(j2.eq.l2)then                 ! b b 
            jl=17
          else                             ! b ab 
            jl=18
          endif
      elseif(j2a.eq.4.and.l2a.eq.5)then  ! c b
            jl=19
      elseif(j2a.eq.5.and.l2a.eq.4)then  ! b c
            jl=20
      else
          stop'ERROR 12082018' 
      endif
      end

c-----------------------------------------------------------------------
      subroutine eolpic( j1 , l1 , j2 , l2 ) !KW1808c mapping  j1 l1 -> j2 l2    
c-----------------------------------------------------------------------
      implicit none
      integer j1,l1,j2,l2
      j2=j1
      l2=l1
      if(abs(j2).ge.3)j2=sign(abs(j2)+1,j2)
      if(abs(l2).ge.3)l2=sign(abs(l2)+1,l2)
      end




c####################################################################################
c####################################################################################
c####################################################################################
c####################################################################################
c######
c######               reorganizing cross sections 
c######
c####################################################################################
c####################################################################################
c####################################################################################
c####################################################################################

c-----------------------------------------------------------------------
      subroutine setQvalues(q1,q2,q3)
c-----------------------------------------------------------------------
      real qqq1,qqq2,qqq3
      common/cqqq/qqq1,qqq2,qqq3
      qqq1=q1
      qqq2=q2
      qqq3=q3
      end

c-----------------------------------------------------------------------
      subroutine tabuSTYPZ(iii,smin,smax,qqs,qqcut) !KW2006 TP2003 
c-----------------------------------------------------------------------
c tabulates sTYPz -> /ctabSTYPZ/
c-----------------------------------------------------------------------
      implicit none
#include "tab.h"
      double precision smin,smax
      real qqs(2),qqcut,sss,scsmax,csbin
      integer iii,n,ncsbin,ncsbmx,jl
      parameter (ncsbmx=100,csbin=8.) !5 bins per decade
      real sTYPzTab,scsmin
      common /ctabSTYPZ/sTYPzTab(mxpsar3,ncsbmx),scsmin
      real sTYPz
      common/zeroSTYP/sTYPz(mxpsar3)
      call setQvalues(qqs(1),qqs(2),qqcut) !needed for call calcSTYP
      scsmin=sngl(smin)
      scsmax=sngl(smax)
      ncsbin=1+int(log10(scsmax/scsmin)*csbin)
      if(ncsbin.ge.ncsbmx-2)stop'increase ncsbmx in tabuSTYPZ'
      do n=1,ncsbin+2
        sss=sngl(scsmin*(scsmax/scsmin)**dble((n-1)/dble(ncsbin)))
        call calcSTYP(iii,sss) !computes sTYPz ->  /zeroSTYP/
        do jl=1,mxpsar3
          sTYPzTab(jl,n)=sTYPz(jl) 
        enddo
      enddo
      end
      
c-----------------------------------------------------------------------
      subroutine pipolSTYPZ(sss) !KW2006 TP2003 
c-----------------------------------------------------------------------
c  (polynomially) interpolates sTYPz -> /zeroSTYP/ 
c-----------------------------------------------------------------------
      implicit none
#include "tab.h"
      real sss,xn,fac,csbin
      integer n,ncsbmx,jl
      parameter (ncsbmx=100,csbin=8.)      !5 bins per decade
      real sTYPzTab,scsmin
      common /ctabSTYPZ/sTYPzTab(mxpsar3,ncsbmx),scsmin
      real sTYPz
      common/zeroSTYP/sTYPz(mxpsar3)
      xn=1.+log10(max(1.,sngl(sss)/scsmin))*csbin
      n=int(xn)
      fac=xn-float(n)
      do jl=1,mxpsar3
        sTYPz(jl)=fac*sTYPzTab(jl,n+1)+(1.-fac)*sTYPzTab(jl,n)
      enddo
      end
      
c-----------------------------------------------------------------------
      subroutine calcSTYP(iii,sss) !KW2006
c-----------------------------------------------------------------------
c computes sTYP -> /fiveSTYP/
c computes sTYPz -> /zeroSTYP/
c-----------------------------------------------------------------------
      implicit none
#include "aaa.h"
#include "sem.h"
#include "tab.h"
      integer ioTestFact
      common /ctestfact/ ioTestFact
      integer laddTestFact
      common /claddtestfact/ laddTestFact
      integer iii,jl,j1,l1,j2,l2,k
      real sss,psboi,psborn,psjti,psjti1,sssTra
      real qqq1,qqq2,qqq3
      common/cqqq/qqq1,qqq2,qqq3
      real sTYP
      common/fiveSTYP/sTYP(0:4,mxpsar3)
      real sTYPz
      common/zeroSTYP/sTYPz(mxpsar3)
      common/csssTra/sssTra

      !--------------------------------------------------------------
      !  sTYP(k,...) different ladder types 
      !  k=1 no emissions at all (only Born)
      !  k=2 at least one emission on current side, no emission on other side
      !  k=3 at least one emission on opposite side, no emission on other side
      !  k=4 at least one emisson on each side (not computed)
      !  k=0 sum of all 
      !  0 = 1 + 2 + 3 + 4
      !--------------------------------------------------------------

      sssTra=sss

      do jl=1,mxpsar3

        call eolpi( jl , j1 , l1 )      
        call eolpic( j1 , l1 , j2 , l2 ) 

        if(iii.ge.0.and.iii.le.3)then
          sTYP(0,jl)=psjti(1,qqq1,qqq3,sss, j2, l2, 0 ) 
          sTYP(1,jl)=psboi(qqq1,qqq2,qqq3,sss, j2, l2, 0 )  
          sTYP(2,jl) = psjti1(qqq1,qqq2,qqq3,sss,j2,l2,0) - sTYP(1,jl)
          sTYP(3,jl) = psjti1(qqq2,qqq1,qqq3,sss,l2,j2,0) - sTYP(1,jl)
          sTYP(4,jl) =  sTYP(0,jl)-sTYP(1,jl)-sTYP(2,jl)-sTYP(3,jl) 
        elseif(iii.eq.-1.or.iii.eq.10)then
          sTYP(0,jl)=psborn(qqq1,qqq2,qqq3,sss, j1, l1, -1, -1 )  
          sTYP(1,jl) = sTYP(0,jl)
          sTYP(2,jl) = 0.
          sTYP(3,jl) = 0.
          sTYP(4,jl) = 0.
        else 
          stop'ERROR 10032019d' 
        endif


      enddo  

      !kkkkkkkkkkk!KW1909
      if(abs(ioTestFact).ge.4)then
        do jl=1,mxpsar3
          if(ioTestFact.eq.4.and.jl.ne.2
     .   .or.ioTestFact.eq.5.and.jl.ne.1
     .   .or.ioTestFact.eq.6.and.jl.ne.2
     .   .or.ioTestFact.ge.7.and.ioTestFact.le.9.and.jl.ne.1
     .   .or.ioTestFact.eq.10.and.jl.ne.2
     .   .or.ioTestFact.eq.11.and.jl.ne.4
     .   )then
            do k=0,4 
              sTYP(k,jl)=0
            enddo
          endif
          if(laddTestFact.gt.0)then
            if(ioTestFact.eq. 4.or.ioTestFact.eq. 5
     .     .or.ioTestFact.eq.-4.or.ioTestFact.eq.-5
     .      )then 
              do k=1,4 
                if(k.ne.laddTestFact)sTYP(k,jl)=0
              enddo
              sTYP(0,jl)=sTYP(laddTestFact,jl)
            endif
          endif
        enddo
      endif
      !kkkkkkkkkkk

      do jl=1,mxpsar3
        sTYPz(jl) = sTYP(0,jl) 
      enddo

      end

c-----------------------------------------------------------------------
      subroutine completeSTYP(iii,ee44) !KW2006
c-----------------------------------------------------------------------
c computes sTYP = sTYP * ee44 ->  /fiveSTYP/
c-----------------------------------------------------------------------
      implicit none
#include "tab.h"
      integer iii,jl,j1,l1,j2,l2,k
      double precision ee44(0:3,0:3)
      real sTYP
      common/fiveSTYP/sTYP(0:4,mxpsar3)
      integer nfjet
      parameter (nfjet=5)
      integer ifgr(-nfjet:nfjet) !trafo to flavor groups 0,1,2,3
      ! flavors -5  -4  -3  -2  -1   0   1   2   3   4   5
      data ifgr/ 3 , 2 , 1 , 1 , 1 , 0 , 1 , 1 , 1 , 2 , 3 / 
      !cross section * Esat1 *Esat2 
      if(iii.ne.3)then
       do jl=1,mxpsar3
        call eolpi( jl , j1 , l1 )      
        call eolpic( j1 , l1 , j2 , l2 ) 
        do k=0,4 
          sTYP(k,jl) = sTYP(k,jl) * ee44( ifgr(j2) , ifgr(l2) )
        enddo
       enddo
      endif
      end

c-----------------------------------------------------------------------
      subroutine eol44(jl,ifgr1,ifgr2) 
c-----------------------------------------------------------------------
      integer nfjet
      parameter (nfjet=5)
      integer ifgr(-nfjet:nfjet) !trafo to flavor groups 0,1,2,3
      ! flavors -5  -4  -3  -2  -1   0   1   2   3   4   5
      data ifgr/ 3 , 2 , 1 , 1 , 1 , 0 , 1 , 1 , 1 , 2 , 3 / 
      call eolpi( jl , j1 , l1 )      
      call eolpic( j1 , l1 , j2 , l2 ) 
      ifgr1=ifgr(j2)
      ifgr2=ifgr(l2)
      end

c-----------------------------------------------------------------------
      subroutine calcWW4(iii,ee44,wwgg,wwgq,wwqg,wwqq)
c-----------------------------------------------------------------------
c taking sTYPz from /zeroSTYP/
c computes the 4 weights wwgg,wwgq,wwqg,wwqq
c containing factors ee44 = Esat1 * Esat2, but NOT for iii=3
c-----------------------------------------------------------------------
      implicit none
#include "sem.h"
#include "tab.h"
      integer iii, jl , j1 , l1, j2 , l2
      double precision ee44(0:3,0:3)
      real wwgg,wwgq,wwqg,wwqq
      real sTYPz
      common/zeroSTYP/sTYPz(mxpsar3)
      real sTYPzz(mxpsar3)
      integer nfjet
      parameter (nfjet=5)
      integer ifgr(-nfjet:nfjet) !trafo to flavor groups 0,1,2,3
      ! flavors -5  -4  -3  -2  -1   0   1   2   3   4   5
      data ifgr/ 3 , 2 , 1 , 1 , 1 , 0 , 1 , 1 , 1 , 2 , 3 / 

      if(iii.ne.3)then
       do jl=1,mxpsar3
        call eolpi( jl , j1 , l1 )      
        call eolpic( j1 , l1 , j2 , l2 ) 
          sTYPzz(jl) = sTYPz(jl) * ee44( ifgr(j2) , ifgr(l2) )
       enddo
      endif

      if(iii.le.0.or.iii.ge.4)then ! case iii = 0  and < 0 and > 3 
        wwgg = sTYPzz(1)                   ! gg  !1
        wwgq = sTYPzz(2) *noflit*2         ! gq  !6
     .       + sTYPzz(7) *2                ! gc  !2
     .       + sTYPzz(9) *2                ! gb  !2
        wwqg = sTYPzz(3) *noflit*2         ! qg  !6
     .       + sTYPzz(8) *2                ! cg  !2
     .       + sTYPzz(10)*2                ! bg  !2  
        wwqq = sTYPzz(4) *noflit*2         ! qq  !6 
     .       + sTYPzz(5) *noflit*2        ! qaq !6
     .       + sTYPzz(6) *noflit*2*(noflit*2-2) ! qqp  ! 24
     .       + sTYPzz(11)*noflit*2*2       ! qc  !12
     .       + sTYPzz(12)*noflit*2*2       ! cq  !12   
     .       + sTYPzz(13)*noflit*2*2       ! qb  !12
     .       + sTYPzz(14)*noflit*2*2       ! bq  !12 
     .       + sTYPzz(15)*2                ! cc  !2
     .       + sTYPzz(16)*2                ! cac !2
     .       + sTYPzz(17)*2                ! bb  !2
     .       + sTYPzz(18)*2                ! bab !2
     .       + sTYPzz(19)*4                ! cb  !4
     .       + sTYPzz(20)*4                ! bc  !4  
                                               ! 11*11 = 121  ! 121 OK
      elseif(iii.eq.1)then ! case iqq=1 val-sea   Here one quark is fixed (=valence)
        wwgg = 0                 
        wwgq = 0                 
        wwqg = sTYPzz(3)                  ! qg     !assuming no valence charm 
        wwqq = sTYPzz(4)                  ! qq
     .       + sTYPzz(5)                  ! qaq
     .       + sTYPzz(6) *(noflit*2-2)    ! qqp
     .       + sTYPzz(11)*2               ! qc
     .       + sTYPzz(13)*2               ! qb             
      elseif(iii.eq.2)then ! case iqq=2 sea-val   Here one quark is fixed (=valence)
        wwgg = 0                  
        wwqg = 0                  
        wwgq = sTYPzz(2)                  ! gq    !assuming no valence charm 
        wwqq = sTYPzz(4)                  ! qq
     .       + sTYPzz(5)                  ! qaq
     .       + sTYPzz(6) *(noflit*2-2)    ! qqp
     .       + sTYPzz(12)*2               ! cq
     .       + sTYPzz(14)*2               ! bq
      elseif(iii.eq.3)then ! case iqq=3 val-val   Here both quarks are fixed (=val-val)
        wwgg = 0
        wwqg = sTYPz(6)   ! qqp    !special meaning (wqg is zero)
        wwgq = sTYPz(5)   ! qaq    !special meaning (wgq is zero)
        wwqq = sTYPz(4)    ! qq
      else
        stop'ERROR 01092018' 
      endif
      end

c-----------------------------------------------------------------------
      subroutine weightsFlavorPairs(iii,jqq,xkk,ykk)
c-----------------------------------------------------------------------
c calculates weights of flavor pairs
c-----------------------------------------------------------------------
      implicit none
#include "sem.h"
#include "tab.h"
      integer iii,jqq,ifmtx,i,j,imax,jl,ibth,k,j1,l1,j2,l2
      integer iDebug
      real psjti1,psboi,psjti,psborn,psjet1,psjet
      real u(0:4),w(0:4)
      real sssTra
      real sTYP
      common/fiveSTYP/sTYP(0:4,mxpsar3)
      real qqq1,qqq2,qqq3
      common/cqqq/qqq1,qqq2,qqq3
      real xkk(11),ykk(0:3,11)
      common/csssTra/sssTra
      integer ncount,lcount
      data ncount/0/
      save ncount
      call getMonitorFileIndex(ifmtx)
      iDebug=0
      ibth=0
      do jl=1,mxpsar3
       !~~~~~~~~~~~~~~~~~~~~~~~~~~~
       if(sTYP(4,jl).lt.0.)then
        if(iDebug.eq.1)then
         !compute it directly
         call eolpi( jl , j1 , l1 )      
         call eolpic( j1 , l1 , j2 , l2 ) 
         !call psjti0(sssTra,u(0), u(1), j2, l2 )  
         u(0) = psjti(1,qqq1,qqq3,sssTra, j2, l2, 0 ) 
         u(1) = psboi(qqq1,qqq2,qqq3,sssTra, j2, l2, 0 )  
         u(2) = psjti1(qqq1,qqq2,qqq3,sssTra,j2,l2,0) - u(1)
         u(3) = psjti1(qqq2,qqq1,qqq3,sssTra,l2,j2,0) - u(1)
         u(4) =  u(0)-u(1)-u(2)-u(3) 
         w(1) = psborn(qqq1,qqq2,qqq3,sssTra,j1,l1,0,0)
         w(2) = psjet1(qqq1,qqq2,qqq3,sssTra,j1,l1,0)
         w(3) = psjet1(qqq2,qqq1,qqq3,sssTra,l1,j1,0)
         w(4) = psjet(qqq1,qqq2,qqq3,sssTra,j1,l1,0)
         w(0) = w(1) + w(2) + w(3) + w(4) 
         print*,'E q1 q2 q3',sssTra,qqq1,qqq2,qqq3
         write(ifmtx,'(a,i5,5f10.5)')'sTYP',jl,(u(k),k=0,4)
         write(ifmtx,'(a,i5,5f10.5)')'    ',jl,(w(k),k=0,4)
         !write(ifmtx,'(a,i5,5f10.4)')'    ',jl,(sTYPx(k,jl),k=0,4)
        endif  
        ibth=1
       endif
       !~~~~~~~~~~~~~~~~~~~~~~~~~~~
      enddo
      !~~~~~~~~~~~~~~~~~~~~~~~~~~~
      if(ibth.eq.1)then
        ncount=ncount+1 
        if(ncount.eq.1.or.ncount.eq.10.or.ncount.eq.100
     .    .or.ncount.eq.1000)then
          lcount=log10(1.*ncount)+1
          write(*,'(a,i1,2a,i2)')'WARNING ',lcount
     .    ,' tab Negative sTYP(4)'
     .    ,'  (tabulation?)  iDebug =',iDebug
        endif
      endif  
      !~~~~~~~~~~~~~~~~~~~~~~~~~~~
      imax=0
      if(iii.le.0)then !----------------------- sea-sea + iii<0 -----------------
       if(jqq.eq.0)then   !gg 
        imax=1
        do k=0,3 
          ykk(k,1)=sTYP(k,1) !gg
        enddo   
       elseif(jqq.eq.1)then !qg
        imax=3
        do k=0,3 
          ykk(k,1)=sTYP(k,3)*noflit*2 !qg
          ykk(k,2)=sTYP(k,8)*2   !cg
          ykk(k,3)=sTYP(k,10)*2  !bg
        enddo
       elseif(jqq.eq.2)then !gq
        imax=3
        do k=0,3 
          ykk(k,1)=sTYP(k,2)*noflit*2 !gq
          ykk(k,2)=sTYP(k,7)*2   !gc
          ykk(k,3)=sTYP(k,9)*2   !gb
        enddo
       elseif(jqq.eq.3)then  !qq
        imax=11
        do k=0,3 
          ykk(k,1) = sTYP(k,4)*noflit*2              
     .              +sTYP(k,5)*noflit*2              
     .              +sTYP(k,6)*noflit*2*(noflit*2-2) !qq
          ykk(k,2) = sTYP(k,11)*noflit*2*2           !qc
          ykk(k,3) = sTYP(k,12)*noflit*2*2           !cq
          ykk(k,4) = sTYP(k,13)*noflit*2*2           !qb
          ykk(k,5) = sTYP(k,14)*noflit*2*2           !bq
          ykk(k,6) = sTYP(k,15)*2                    !cc
          ykk(k,7) = sTYP(k,16)*2                    !cac
          ykk(k,8) = sTYP(k,17)*2                    !bb
          ykk(k,9) = sTYP(k,18)*2                    !bab
          ykk(k,10) = sTYP(k,19)*4                   !cb
          ykk(k,11) = sTYP(k,20)*4                   !bc
        enddo
       endif
      elseif(iii.eq.1)then !------------------------ val-sea ------------------------
       if(jqq.eq.0)then  !vg
        imax=1
        do k=0,3 
          ykk(k,1)=sTYP(k,3) !qg
        enddo
       elseif(jqq.eq.1)then !vq
        imax=5
        do k=0,3 
          ykk(k,1) = sTYP(k,4)               !qq
          ykk(k,2) = sTYP(k,5)               !qaq
          ykk(k,3) = sTYP(k,6)*(noflit*2-2)  !qqp
          ykk(k,4) = sTYP(k,11)*2            !qc
          ykk(k,5) = sTYP(k,13)*2            !qb
        enddo
       else
         stop'ERROR 28062019a'
       endif
      elseif(iii.eq.2)then !----- sea-val
       if(jqq.eq.0)then  !gv
        imax=1
        do k=0,3 
          ykk(k,1)=sTYP(k,2) !gq
        enddo
       elseif(jqq.eq.1)then !qv
        imax=5
        do k=0,3 
          ykk(k,1) = sTYP(k,4)              !qq
          ykk(k,2) = sTYP(k,5)              !qaq
          ykk(k,3) = sTYP(k,6)*(noflit*2-2) !qqp
          ykk(k,4) = sTYP(k,12)*2           !cq
          ykk(k,5) = sTYP(k,14)*2           !bq
        enddo
       else
         stop'ERROR 28062019b'
       endif
      elseif(iii.eq.3)then !--------------------- val-val -------------------------
        imax=3
        do k=0,3 
          ykk(k,1) = sTYP(k,4) !qq
          ykk(k,2) = sTYP(k,5) !qaq
          ykk(k,3) = sTYP(k,6) !qqp
        enddo
      else !-----------------------------------------------------------------------
        write(ifmtx,*)'####### ERROR in weightsFlavorPairs #######'
        write(ifmtx,*)'-----> iii = ',iii
      endif
      if(imax.gt.0)then 
       do i=1,imax
         xkk(i)=ykk(0,i)
       enddo
       do j=1,3
        do i=1,imax
          if( xkk(i).ne.0)then
            ykk(j,i) = min(1., ykk(j,i) / xkk(i) )
          else
            ykk(j,i) = 0
          endif
        enddo
       enddo
      endif
      end



c####################################################################################
c####################################################################################
c####################################################################################
c####################################################################################
c######
c######         get tables for given Q2s ( q2imin(2) ) via interpolation
c######
c######                        from initialized tables (from files)
c######
c####################################################################################
c####################################################################################
c####################################################################################
c####################################################################################
      
c-----------------------------------------------------------------------
      subroutine ipoCSzeroTables(q2imin)    ! (includes x dep)
c-----------------------------------------------------------------------
c      compute tables cstotzero, csborzer
c        via interpolation for given q2imin
c-----------------------------------------------------------------------
#include "aaa.h"
#include "sem.h"
#include "tab.h"
c      real wi(0:1),wj(0:1)
      real wi(0:2),wj(0:2)
      real q2imin(2)
      if(q2imin(1).le.0.000)then      
      do k=1,20 
      do ml=1,mxpsar3  !KW1808c
      cstotzero(k,ml)=-99999.
      csborzer( k,ml)=-99999.
      enddo
      enddo
      return
      endif
     
ctp Do not redefine q2cmin (fixed outside)    
c      do ii=1,2
c      q2cmin(ii)=q2imin(ii)
c      enddo

      q=log(max(1e-6,q2imin(1)))
c      ipol=1
c      call ipoli(q,nq,frac)
c      wi(0)=1.-frac
c      wi(1)=frac
      ipol=2
      call ipoli3(q,nq,frac)
      wi(1)=frac
      wi(2)=wi(1)*(wi(1)-1.)*.5
      wi(0)=1.-wi(1)+wi(2)
      wi(1)=wi(1)-2.*wi(2)
      n1=nq

      q=log(max(1e-6,q2imin(2)))
c      call ipoli(q,nq,frac)
c      wj(0)=1.-frac
c      wj(1)=frac
      call ipoli3(q,nq,frac)
      wj(1)=frac
      wj(2)=wj(1)*(wj(1)-1.)*.5
      wj(0)=1.-wj(1)+wj(2)
      wj(1)=wj(1)-2.*wj(2)
      n2=nq
      
      do k=1,20 
      do ml=1,mxpsar3 !KW1808c
      cstotzero(k,ml)=0
      csborzer(k,ml)=0
      do i1=0,ipol
      do i2=0,ipol
      cstotzero(k,ml)=cstotzero(k,ml)
     .       + wi(i1)*wj(i2)*csqtotzero(n1+i1,n2+i2,k,ml)
      csborzer(k,ml) = csborzer(k,ml)
     .       + wi(i1)*wj(i2)*csqborzer(n1+i1,n2+i2,k,ml)
      enddo
      enddo
c      if(k.eq.1.and.m.eq.1.and.q2imin(1).gt.22.)then
c      endif
      enddo
      enddo
      
      end

c-----------------------------------------------------------------------
      subroutine ipoOm5Tables(iqq)   ! (without x dep)
c-----------------------------------------------------------------------
c      computes tables fhgg, fhqg, fhgq, fhqq
c        via interpolation for given q2imin
c    iqq <= 0 reset tables
c-----------------------------------------------------------------------
#include "aaa.h"
#include "sem.h"
#include "tab.h"
      common /psar4/  fhgg(11,100,10,8),fhqg(11,100,10,8)
      common /psar4b/ fhgq(11,100,10,8),fhqq(11,100,8)
      common /psar4c/ fhss(11,100,10,8)
c      real wi(0:1),wj(1:1),wk(0:1)
      real wi(0:2),wj(1:1),wk(0:2)
      real q2imin(2)
      if(iqq.le.0)then      
      do m=mnclpt,mxclpt
      do k=1,10
      do i2=1,10
      do i1=1,10
      i12=i1+10*(i2-1)
      do i=1,11
        fhss(i,i12,k,m)=-99999.
        fhgg(i,i12,k,m)=-99999.
        fhqg(i,i12,k,m)=-99999.
        fhgq(i,i12,k,m)=-99999.
        fhqq(i,i12,m)=-99999.
      enddo  
      enddo
      enddo
      enddo   
      enddo   
      q2pmin(1)=q2nmin
      q2pmin(2)=q2nmin
      return
      endif

      call ipoBasics(q2pmin)   !because psvin depends on alpff and betff

      do ii=1,2
      q2imin(ii)=q2pmin(ii)
      enddo
      
      q=log(max(1e-6,q2imin(1)))
      ipol=1
      call ipoli(q,nq,frac)
      wi(0)=1.-frac
      wi(1)=frac
c      ipol=2
c      call ipoli3(q,nq,frac)
c      wi(1)=frac
c      wi(2)=wi(1)*(wi(1)-1.)*.5
c      wi(0)=1.-wi(1)+wi(2)
c      wi(1)=wi(1)-2.*wi(2)
      n1=nq

      q=log(max(1e-6,q2imin(2)))
      call ipoli(q,nq,frac)
      wk(0)=1.-frac
      wk(1)=frac
c      call ipoli3(q,nq,frac)
c      wk(1)=frac
c      wk(2)=wk(1)*(wk(1)-1.)*.5
c      wk(0)=1.-wk(1)+wk(2)
c      wk(1)=wk(1)-2.*wk(2)
      n2=nq
     
      wj(1)=1   
      nz=0
      n=1
      
c      do m=mnclpt,mxclpt
c      do k=1,10
c      do i2=1,10
c      do i1=1,10
c      i12=i1+10*(i2-1)
c      do i=1,11
c        fhss(i,i12,k,m)=0.
c        fhgg(i,i12,k,m)=0.
c        fhqg(i,i12,k,m)=0.
c        fhgq(i,i12,k,m)=0.
c        fhqq(i,i12,m)=0.
c      enddo  
c      enddo
c      enddo
c      enddo   
c      enddo   


      l=mnclpt
c      do l=mnclpt,mxclpt
         lx=1 !l-mnclpt+1   !to be restored when pion/kaon projectiles will be set again
      do k=1,10
      do i2=1,10
      do i1=1,10
      i12=i1+10*(i2-1)
      do i=1,11
        fhss(i,i12,k,l)=0.
        fhgg(i,i12,k,l)=0.
        fhqg(i,i12,k,l)=0.
        fhgq(i,i12,k,l)=0.
c         do n=nn0,nn1     
         do m2=0,ipol
         do m1=0,ipol
          fhss(i,i12,k,l)=fhss(i,i12,k,l)
     .       +wi(m1)*wk(m2)*wj(n)*fhxss(n1+m1,n2+m2,nz+n,i,i12,k,lx)
          fhgg(i,i12,k,l)=fhgg(i,i12,k,l)
     .       +wi(m1)*wk(m2)*wj(n)*fhxgg(n1+m1,n2+m2,nz+n,i,i12,k,lx)
          fhqg(i,i12,k,l)=fhqg(i,i12,k,l)
     .       +wi(m1)*wk(m2)*wj(n)*fhxqg(n1+m1,n2+m2,nz+n,i,i12,k,lx)
          fhgq(i,i12,k,l)=fhgq(i,i12,k,l)
     .       +wi(m1)*wk(m2)*wj(n)*fhxgq(n1+m1,n2+m2,nz+n,i,i12,k,lx)
         enddo
         enddo
c         enddo
      enddo
      enddo 
      enddo
      enddo   
c      enddo  
 
c      do l=mnclpt,mxclpt
c         lx=l-mnclpt+1
      do i2=1,10
      do i1=1,10
      i12=i1+10*(i2-1)
      do i=1,11
        fhqq(i,i12,l)=0.
c         do n=nn0,nn1     
         do m2=0,ipol
         do m1=0,ipol
          fhqq(i,i12,l)=fhqq(i,i12,l)
     .       +wi(m1)*wk(m2)*wj(n)*fhxqq(n1+m1,n2+m2,nz+n,i,i12,lx)
         enddo
         enddo
c         enddo
      enddo 
      enddo
      enddo   
c      enddo   

      end

c-----------------------------------------------------------------------
      subroutine ipoOm5Tab(xp,xm,b)
c-----------------------------------------------------------------------
c      computes tables fhgg, fhqg, fhgq, fhqq
c        via interpolation for given q2imin, xp, xm, and b to speed
c        up calculation (only useful bins as defined as in psvin)
c-----------------------------------------------------------------------
#include "aaa.h"
#include "sem.h"
#include "tab.h"
      common /psar4/  fhgg(11,100,10,8),fhqg(11,100,10,8)
      common /psar4b/ fhgq(11,100,10,8),fhqq(11,100,8)
      common /psar4c/ fhss(11,100,10,8)
c      real wi(0:1),wk(0:1)
      real wi(0:2),wk(0:2)
      real q2imin(2)
      double precision xp,xm
      dimension i(3),j(3)
      common /psar2/  edmax,epmax
      data xps/0./,xms/0./,bs/0./
      save xps,xms,bs,jz,k,i,j

      xpp=sngl(xp)
      xmp=sngl(xm)

      call ipoBasics(q2pmin)   !because psvin depends on alpff and betff

c     if(abs(xpp-xps).gt.1e-6.or.abs(xmp-xms).gt.1e-6
c     .   .or.abs(b-bs).gt.1e-6)then

        sy0=engy**2
        sy=sy0*sngl(xp*xm)
        spmin0=4.*q2zmin
        if(iclpro.ne.4)then
          spmin=spmin0
        else
          spmin=spmin0+2.*qcmass**2
        endif

c impact parameter bin
        rp=r2had(iclpro)+r2had(icltar)+slopom*log(max(1.,sy))
        zb=exp(-b**2/(4.*.0389*rp))
        jz=int(10.*zb)
        if(jz.gt.8)jz=8
        if(jz.lt.1)jz=1

c energy index
        yl=log(sy0/spmin)/log(epmax/2./spmin)*10.+1
        k=int(yl)
        if(k.gt.9)k=9
        if(k.lt.2)k=2 !otherwise following log(syx)=0 -> crash

c xp and xm index
        do k1=1,3
          k2=k+k1-1
          syx=(epmax/2./spmin0)**((k2-1)/10.) !cancellation *spmin0/spmin0
          xl1=log(xpp*syx)/log(syx)*9.+1
          i(k1)=max(1,int(xl1))
          i(k1)=min(8,i(k1))
          xl2=log(xmp*syx)/log(syx)*9.+1
          j(k1)=max(1,int(xl2))
          j(k1)=min(8,j(k1))
        enddo

c save current values
c        xps=xpp
c        xms=xmp
c        bs=b
c      endif

      jmin=jz
      jmax=jz+2

      do ii=1,2
      q2imin(ii)=q2pmin(ii)
      enddo
      
      q=log(max(1e-6,q2imin(1)))
      ipol=1
      call ipoli(q,nq,frac)
      wi(0)=1.-frac
      wi(1)=frac
c      ipol=2
c      call ipoli3(q,nq,frac)
c      wi(1)=frac
c      wi(2)=wi(1)*(wi(1)-1.)*.5
c      wi(0)=1.-wi(1)+wi(2)
c      wi(1)=wi(1)-2.*wi(2)
      n1=nq

      q=log(max(1e-6,q2imin(2)))
      call ipoli(q,nq,frac)
      wk(0)=1.-frac
      wk(1)=frac
c      call ipoli3(q,nq,frac)
c      wk(1)=frac
c      wk(2)=wk(1)*(wk(1)-1.)*.5
c      wk(0)=1.-wk(1)+wk(2)
c      wk(1)=wk(1)-2.*wk(2)
      n2=nq

c      print *,'om5tab',q2pmin,wi,wk
     
      nz=1

      l=iclpro+4*(icltar-1)
      lx=l-mnclpt+1

      do jb=jmin,jmax
        do k1=1,3
          k2=k+k1-1
          do i1=1,3
          do j1=1,3
            i2=i(k1)+i1-1
            j2=j(k1)+j1-1
            i12=i2+10*(j2-1)
            fhss(k2,i12,jb,l)=0.
            fhgg(k2,i12,jb,l)=0.
            fhqg(k2,i12,jb,l)=0.
            fhgq(k2,i12,jb,l)=0.
            do m2=0,ipol
            do m1=0,ipol
              fhss(k2,i12,jb,l)=fhss(k2,i12,jb,l)
     .             +wi(m1)*wk(m2)*fhxss(n1+m1,n2+m2,nz,k2,i12,jb,lx)
              fhgg(k2,i12,jb,l)=fhgg(k2,i12,jb,l)
     .       +wi(m1)*wk(m2)*fhxgg(n1+m1,n2+m2,nz,k2,i12,jb,lx)
              fhqg(k2,i12,jb,l)=fhqg(k2,i12,jb,l)
     .       +wi(m1)*wk(m2)*fhxqg(n1+m1,n2+m2,nz,k2,i12,jb,lx)
              fhgq(k2,i12,jb,l)=fhgq(k2,i12,jb,l)
     .       +wi(m1)*wk(m2)*fhxgq(n1+m1,n2+m2,nz,k2,i12,jb,lx)
            enddo
            enddo
c            print *,k2,i12,jb,l,fhss(k2,i12,jb,l)
          enddo
          enddo 
        enddo
      enddo   
 
      do k1=1,3
        k2=k+k1-1
        do i1=1,3
        do j1=1,3
          i2=i(k1)+i1-1
          j2=j(k1)+j1-1
          i12=i2+10*(j2-1)
          fhqq(k2,i12,l)=0.
          do m2=0,ipol
          do m1=0,ipol
            fhqq(k2,i12,l)=fhqq(k2,i12,l)
     .       +wi(m1)*wk(m2)*fhxqq(n1+m1,n2+m2,nz,k2,i12,lx)
          enddo
          enddo
        enddo 
        enddo
      enddo   

      end




c####################################################################################
c####################################################################################
c####################################################################################
c####################################################################################
c####
c####                       initialize tables
c####
c####################################################################################
c####################################################################################
c####################################################################################
c####################################################################################



c#######################################################################
c#######################################################################

               subroutine mkEvTable

c#######################################################################
c#######################################################################

c calculates ev.i = QCD evolution tables.
c Tables depend on q2zmin only for binning. Doesn't depend on q2nmin.
c-----------------------------------------------------------------------
#include "aaa.h"
      logical lcalc
      double precision evk  !,xvalueEvTab
      common /psar15/ sudx(40,2)
      common /evkarr/ evk(40,40,100,3,2)
      common /psar2/  edmax,epmax
      common/producetab/ producetables              !used to link with CRMC
      logical producetables
#include "sem.h"
      logical nequal

c      double precision sudtim,sudspa
c      do i=1,40,4
c        qi=q2zmin*exp(.5*(i-1))
c        a=sudspa(qi,0) 
c        b=sudspa(qi,1)
c        print*,qi,a,b
c      enddo

      call makeSudx

      if(ioangord.eq.0)then
        fnie=fnie(1:index(fnie,'ev.i')-3)//'n/ev.i '   
      elseif(ioangord.eq.1)then
        fnie=fnie(1:index(fnie,'ev.i')-3)//'o/ev.i '   
      else
        stop'####### ERROR 07122018a #######'
      endif
      nfnie=index(fnie,' ')-1

      write(ifmt,'(2a)')'load ',fnie(1:nfnie)
      inquire(file=fnie(1:nfnie),exist=lcalc)
      if(lcalc)then
       if(inicnt.eq.1)then
        open(1,file=fnie(1:nfnie),status='old')
        read (1,*)qcdlam0,q2zmin0,q2ini0,naflav0,epmax0
        if(qcdlam0.ne.qcdlam)write(ifmt,'(a)')'iniev: wrong qcdlam'
        if(nequal(q2zmin0,q2zmin))write(ifmt,'(a)')'iniev: wrong q2zmin'
        if(q2ini0 .ne.q2ini )write(ifmt,'(a)')'iniev: wrong q2ini'
        if(naflav0.ne.naflav)write(ifmt,'(a)')'iniev: wrong naflav'
        if(epmax0 .ne.epmax )write(ifmt,'(a)')'iniev: wrong epmax'
        if(qcdlam0.ne.qcdlam.or.nequal(q2zmin0,q2zmin)
     *  .or.q2ini0.ne.q2ini.or.naflav0.ne.naflav
     *  .or.epmax0.ne.epmax)then
           write(6,'(//a//)')'   iniev has to be reinitialized!!!'
           stop
        endif
        read (1,*)evk
        close(1)
       endif
       goto 101
      elseif(.not.producetables)then
        write(ifmt,*) "Missing epos ev file !"        
        write(ifmt,*) "Please correct the defined path ",
     &"or force production ..."
      endif

      write(ifmt,'(a)')'iniev does not exist -> calculate tables  ...'
      call fillEvTab

      write(ifmt,'(a)')'write to iniev ...'
      open(1,file=fnie(1:nfnie),status='unknown')
      write (1,*)qcdlam,q2zmin,q2ini,naflav,epmax
      write (1,*)evk
      close(1)

101   continue

      do i=1,40
      do j=1,40
      do l=1,100
      !if(i.eq.1.and.j.eq.1)print*,'x(l) ',l,xvalueEvTab(l,10d0)
      do m=1,3
      do k=1,2
       a=evk(i,j,l,m,k)
       if(.not.(a.le.0..or.a.ge.0.))then
        stop'ERROR 25022020 NaN detected'
       endif
      enddo
      enddo
      enddo
      enddo
      enddo

      end 

c-----------------------------------------------------------------------
      double precision function xvalueEvTab(l,qq) 
c-----------------------------------------------------------------------
      double precision spmax,epscutSud2,xx,qq
      real val
      ll=min(l,99)
      call getEpmax(val)
      spmax=val
      if(ll.le.37)then
       xx=.1d0/(.1d0*spmax)**((37.d0-ll)/36.d0)
      elseif(ll.le.69)then
       xx=.1d0+.8d0*(ll-37.d0)/32.d0
      else
       xx=1.d0-.1d0*(10.d0*epscutSud2(qq))**((ll-69.d0)/31.d0)
      endif
      xvalueEvTab=xx
      end

c-----------------------------------------------------------------------
      subroutine fillEvTab !fill array evk
c-----------------------------------------------------------------------
      implicit double precision (a-h,o-z)
#include "aaa.h"
      real val
      common /evkarr/ evk(40,40,100,3,2)
      dimension ixemax(40,3,2),evs(40,100,3,2)
      call getQcdlam(val)
      alm=val
      call getEpmax(val)
      spmax=val
      do i=1,40
      do m=1,3
      do k=1,2
       ixemax(i,m,k)=99
      do j=1,40
      do l=1,100
       evk(i,j,l,m,k)=0.d0
      enddo
      enddo
      enddo
      enddo
      enddo
  
      !------------------------------------!-----------!
      ! compute and tabulate ev(qi,qj,x)   !  * * * *  !
      ! for a triangle in the qi-qj plane: !  * * *    !
      !         1 <= qi < Qmax             !  * *      !
      !        qi <= qj < Qmax             !  *        !
      !------------------------------------!-----------!

      !--------------------------------------------------
      !QGSJET method (S. Ostapchenko):
      !ev(a,b,x) = ev(a,c,x)*sud(c,b) + sud(a,c)*ev(c,b,x)
      !          +  \int dz/z ev(a,c,x/z)*ev(c,b,z)
      !    with a=qi, b=qj, and c=qm with m=(j-1)
      !For j>3 no iteration needed, since ev is known.
      !Equation can be easily derived from DGLAP, writing first ev as
      !an infite series of multiple emissions, and splitting the multi
      !dimensional integration domain to <c and >c.  
      !--------------------------------------------------
      !For j=2 we use iterative method, PR B31, B32 
      !--------------------------------------------------

      n=1
1     n=n+1     ! <---- iterate ---- 

      !---------- i=1,39 j=2 -----------
      write(ifmt,'(a,i3,a)')' n =',n,'   j = 2'       
      do m=1,3
      do k=1,2
       if(m.ne.3.or.k.ne.1)then 
        do i=1,39
         if(ixemax(i,m,k).gt.0)then
          qi=spmax**((i-1)/39.d0)           ! qi
          qj=qi*(spmax/qi)**(1.d0/39.d0)    ! qj2
          do l=1,99
           xx=xvalueEvTab(l,qj)
           ev= psev2(qi,qj,xx,m,k)       
           evs(i,l,m,k)=
     .     dlog( (psevz2(qi,qj,xx,m,k)+ev) / psevz2(qi,qj,xx,m,k) )
          enddo
         endif
        enddo
       endif
      enddo
      enddo
      jec=0
      do m=1,3
      do k=1,2
       if(m.ne.3.or.k.ne.1)then
        do i=1,39
         if(ixemax(i,m,k).gt.0)then
          do l=ixemax(i,m,k),1,-1
           if(abs(evs(i,l,m,k)-evk(i,2,l,m,k)).gt.1.d-3)then
            evk(i,2,l,m,k)=evs(i,l,m,k)
            jec=1
           elseif(ixemax(i,m,k).eq.l)then
            ixemax(i,m,k)=l-1
           endif
          enddo
         endif
        enddo
       endif
      enddo
      enddo

      !----------- i=1,39 j=3 --------
      write(ifmt,'(a,i3,a)')' n =',n,'   j = 3'       
      do i=1,39
       qi=spmax**((i-1)/39.d0)
       qm=qi*(spmax/qi)**(1.d0/39.d0)   !qj2
       qj=qi*(spmax/qi)**(2.d0/39.d0)   !qj3
       do l=99,1,-1
        xx=xvalueEvTab(l,qj)
        do m=1,3
        do k=1,2
         if(m.ne.3.or.k.ne.1)then
          ev=psevev2(qi,qm,qj,xx,m,k)
     *    +psevi2(qi,qm,xx,m,k) * psudi2(qj,k-1)/psudi2(qm,k-1)
     *    +psevi2(qm,qj,xx,m,k) * psudi2(qm,m-1)/psudi2(qi,m-1)
          evk(i,3,l,m,k)=dlog(ev/psevz2(qi,qj,xx,m,k))
          a=evk(i,3,l,m,k)
          if(.not.(a.le.0..or.a.ge.0.))stop'ERROR 26022020a'
         endif
        enddo
        enddo
       enddo
      enddo

      if(jec.ne.0)goto 1 !---- iterate ---->

      !-------- i=1,39 j=4,5... -----------  
      write(ifmt,'(a,i3,a)')' n =',n,'   j > 3'       
      !~~~~~~~~~~~~~~~~print~~~~~~~~~ 
      lll=30  
      iii=1    
      if(lll.gt.0)then      
       write(ifmt,'()')                     
       qi=spmax**((iii-1)/39.d0)
       xx=xvalueEvTab(lll,qi)
       write(ifmt,'(a,e10.2)')'EVK x = ',xx                     
       do j=1,3 
       qj=qi*(spmax/qi)**((j-1)/39.d0) 
       write(ifmt,'(a,i3,e10.2,2(4x,2i2,2x,3e10.2))')'EVK',j,qj
     . ,1,1,psevi2(qi,qj,xx,1,1),0.,psevz2(qi,qj,xx,1,1)
     . ,1,2,psevi2(qi,qj,xx,1,2),0.,psevz2(qi,qj,xx,1,2)
       enddo
      endif                                                
      !~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ 
      do i=1,39
       qi=spmax**((i-1)/39.d0)
      do j=4,40
       qm=qi*(spmax/qi)**((j-2)/39.d0)   !qj3...
       qj=qi*(spmax/qi)**((j-1)/39.d0)   !qj4...
       do l=99,1,-1
        xx=xvalueEvTab(l,qj)
        if(i.eq.iii.and.l.eq.lll)write(ifmt,'(a,i3,e10.2,$)')'EVK',j,qj !EVK print
        do m=1,3
        do k=1,2
         if(m.ne.3.or.k.ne.1)then
          ev=psevev2(qi,qm,qj,xx,m,k)
     *    +psevi2(qi,qm,xx,m,k) * psudi2(qj,k-1)/psudi2(qm,k-1)
     *    +psevi2(qm,qj,xx,m,k) * psudi2(qm,m-1)/psudi2(qi,m-1)
          evk(i,j,l,m,k)=dlog(ev/psevz2(qi,qj,xx,m,k))
          !~~~~~~~~~~print~~~~~~~~~~~~~~~ 
          if(i.eq.iii.and.l.eq.lll)then                           
           !write(ifmt,'(1x,e8.2,$)')ev      
           if(m.eq.1)write(ifmt,'(4x,2i2,2x,3e10.2,$)')
     .     m,k,psevi2(qi,qj,xx,m,k),ev,psevz2(qi,qj,xx,m,k)
          endif                                                
          !~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ 
          a=evk(i,j,l,m,k)
          if(.not.(a.le.0..or.a.ge.0.))stop'ERROR 26022020b'
         endif
        enddo
        enddo
        if(i.eq.iii.and.l.eq.lll)write(ifmt,*)' '                 !EVK print
       enddo
      enddo
      enddo

      end     


c#######################################################################
c#######################################################################

      subroutine mkCsOmTables

c#######################################################################
c#######################################################################

#include "aaa.h"
#include "par.h"
#include "sem.h"
#include "tab.h"
      logical lcalc!,lcalc2
c      double precision om5p,xh,yh,v3pom(4),om2p
      dimension gamhad0(nclha),r2had0(nclha),alplea0(nclha)
      common /psar2/  edmax,epmax
      common /psar4/  fhgg(11,100,10,8),fhqg(11,100,10,8)
      common /psar4b/ fhgq(11,100,10,8),fhqq(11,100,8)
      common /psar4c/ fhss(11,100,10,8)
      common /psar7/  delx,alam3p,gam3p
      common /psar9/  alpr
      common /psar15/ sudx(40,2)
      common /psar25/ csdsi(21,21,104)
      common /psar27/ csds(21,26,4),csdt(21,26,2),csdr(21,26,2)
      common /psar33/ asect(7,4,7),asectn(7,7,7)
      common /psar34/ rrr,rrrm
      common /psar35/ anorm,anormp
      common /psar41/ rrrp,rrrmp
      common /psar37/ coefom1,coefom2
      common /psar38/ vfro(11,14,3,2)
      common /psar39/ vnorm(11,14,3,2,2)
      common /ar3/    x1(7),a1(7)
      common /ckinirj/kinirj
      common/ciprotectinirj/iprotectinirj
      common/producetab/ producetables              !used to link with CRMC
      logical producetables,negjdis
      logical nequal
      common/cpriom/npriom
      character*4 ch4
      character*2 ch2
      npriom=0
      
      call utpri('mkCsOm',ish,ishini,4)
      
      !----------------------------------------------------------------
      !call createomtab(4,maxq2mn,maxq2mn,maxzzvex,11,100
      !. ,mxclpt-mnclpt+1,10)
      !istatom=1
      !----------------------------------------------------------------

      if(iappl.eq.6)goto 9999

      !----------------------------------------------------------------
      !
      !  bb  --    CS tabulation    (csbor)
      !
      !  Tables depend on q2zmin only for binning. Doesn't depend on q2nmin
      !----------------------------------------------------------------

      !--------------------------------!
      ! generating bb tables is fast   ! 
      ! can be done interactively      !
      !        use nq2mnfix = 0        !  
      !--------------------------------!

      call makeSudx

      if(ish.ge.4)write(ifch,*)'bare cross sections ...'

      if(ioangord.eq.0)then
        fnii=fnii(1:index(fnii,'aa.i')-3)//'n/bb.i '   
      elseif(ioangord.eq.1)then
        fnii=fnii(1:index(fnii,'aa.i')-3)//'o/bb.i '   
      else
        stop'####### ERROR 07122018c #######'
      endif
      nfnii=index(fnii,' ')-1

      write(ifmt,'(2a)')'load ',fnii(1:nfnii)
      inquire(file=fnii(1:nfnii),exist=lcalc)

      negjdis=.false.    !compute jdis<0 only
      call pdfparamget(1,facpdf4) 
      call pdfparamget(3,facpdf5) 
      if(lcalc)then
        negjdis=.false.
        open(1,file=fnii(1:nfnii),status='old')
        read (1,*)qcdlam0,q2zmin0,q2ini0,naflav0,facpdf40,facpdf50
     *    ,epmax0,pt2cut0,qcmass0,qbmass0,nbflav0,factjpsi0,q2mnval0
        if(qcdlam0.ne.qcdlam)write(ifmt,'(a)')'table bb: wrong qcdlam'
        if(nequal(q2zmin0,q2zmin))
     .   write(ifmt,'(a)')'table bb: wrong q2zmin'
        if(q2ini0 .ne.q2ini )write(ifmt,'(a)')'table bb: wrong q2ini'
        if(naflav0.ne.naflav)write(ifmt,'(a)')'table bb: wrong naflav'
        if(nbflav0.ne.nbflav)write(ifmt,'(a)')'table bb: wrong nbflav'
        if(epmax0 .ne.epmax )write(ifmt,'(a)')'table bb: wrong epmax'
        if(pt2cut0.ne.pt2cut)write(ifmt,'(a)')'table bb: wrong pt2cut'
        if(qcmass0.ne.qcmass)write(ifmt,'(a)')'table bb: wrong qcmass'
        if(qbmass0.ne.qbmass)write(ifmt,'(a)')'table bb: wrong qbmass' !bg
        if(facpdf40.ne.facpdf4)then
          write(ifmt,'(/80a)')('-',k=1,80)
          write(ifmt,'(a)')'    ERROR Wrong facpdf4 for KWn tables '
          write(ifmt,'(7x,f5.2,a)')facpdf40,' expected'   
          write(ifmt,'(a)')'      Redo bb, aa, ab, ac, rj tables '
          write(ifmt,'(80a/)')('-',k=1,80)
          stop
        endif 
        if(facpdf50.ne.facpdf5)then
          write(ifmt,'(/80a)')('-',k=1,80)
          write(ifmt,'(a)')'    ERROR Wrong facpdf5 for KWn tables '
          write(ifmt,'(7x,f5.2,a)')facpdf50,' expected'   
          write(ifmt,'(a)')'      Redo bb, aa, ab, ac, rj tables '
          write(ifmt,'(80a/)')('-',k=1,80)
          stop
        endif 
c        if((factjpsi0.ne.factjpsi.or.q2mnval(maxq2mn).ne.q2mnval0)
c     *                                         .and.isetcs.ge.0)then
c          write(ifmt,'(a)')'table bb: change for jdis<0, update table'
c          negjdis=.true.
c        endif
        if(qcdlam0.ne.qcdlam.or.nequal(q2zmin0,q2zmin)
     *  .or.q2ini0 .ne.q2ini
     *  .or.naflav0.ne.naflav.or.epmax0 .ne.epmax.or. pt2cut.ne.pt2cut0
     *  .or.qcmass0.ne.qcmass.or.qbmass0.ne.qbmass.or.nbflav0.ne.nbflav)!bg
     *  then
          write(ifmt,'(//a//)')'   bb table has to be reinitialized!!'
          stop
        endif
        read (1,*)csbor,csbor0
        close(1)
        !~~~~~~~~~~TEST~~~~~~~~~~~~~~~
        iprik1=19
        iprik2=19
        ipril=999
        if(ipril.eq.141)goto 141 
        !~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
        if(.not.negjdis)goto 1
      elseif(.not.producetables)then
        write(ifmt,*) "Missing epos bb file !"        
        write(ifmt,*) "Please correct the defined path ",
     &"or force production ..."
        stop
      endif

      write(ifmt,'(a/a)')fnii(1:nfnii)//' does not exist '
     . ,'         -> calculate CS tables bb ...'

      ipri=0
      ml1=1
      ml2=mxpsar3
      goto 142
 141  continue
      ipri=1
      ml1=iprik1
      ml2=iprik2
 142  continue

      write(ifmt,*)'-------------------'
      write(ifmt,*)'Born xsection csbor'
      write(ifmt,*)'-------------------'

      if(ipri.eq.1)goto 143

      do ml=1,mxpsar3  !KW1808c 
        call eolpi( ml , m1 , l1 ) !KW1808
        call getQminEoL(ml,qref1,qref2,qqref)
        call getSminEoL(ml,qqref,spminx)
        spminx=spminx*1.10
      do k=1,20
        sk=spminx*(epmax/2./spminx)**((k-1)/19.)
        call getQmaxEoL(ml,sk,qmax)
        qmax=qmax / 1.01!1.10 
        qmax=max(qmax,1.01*qqref)
      do i=1,20
        qq=qqref*(qmax/qqref)**((i-1)/19.)
        q2cmin(1)=qq
        q2cmin(2)=qq
        k1=k+20*(ml-1) !KW1808c
        !tabulation pour q1=q2=qq
        csbor(i,k1,1) 
     .   =log(max(1.e-30,psborn(qq,qq,qq,sk,m1,l1,0,0)))
        !csbor(i,k1,2)
        !.   =log(max(1.e-30,psborn(4.*qq,qq,qq,sk,m1,l1,1,0)))
        qqc=q2zmin*(q2mnval(maxq2mn)/q2zmin)**((i-1)/19.)
        q2cmin(1)=qqc
        q2cmin(2)=qqc
        !do j=1,20
        !  qq0=q2zmin*(qqc/q2zmin)**((j-1)/19.)
        !  csbor0(j,i,k1) 
        !.         =log(max(1.e-30,psborn(qq0,qq0,qqc,sk,m1,l1,-1,0))) !KW2019: this seems not making any sense
        !                                                                   ! better discuss it
        !enddo
      enddo
      enddo
      enddo

            !~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
            write (ifmt,*)'tests:'   !Born
            !~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
      
 143        continue           
            write (ifmt,'(a,a,a)')' jl jdis   sk   '
     .       ,'         q1         q2         qq    '
     .       ,'    born     born-i   diff'
            nckeck=0
            nshow=0
            diffsum=0
            !do ii=1,8
            !  print*,'BORN E-dep',psborn(5.,5.,5.,10.**ii,0,0,0,0)
            !enddo
            do jl=ml1,ml2   !at least once check all !!!
              call eolpi( jl , j1 , l1 ) !KW1808
              call eolpic( j1 , l1 , j2 , l2 ) 
              qqref=q2zmin
              call getSminEoL(jl,qqref,spminx)
              qmin=1.5*q2zmin
              ki=0
              iret=1
              do while(iret.eq.1)
                ki=ki+1
                sk=spminx*(epmax/2./spminx)**((ki-1)/19.) 
                call getLimitsEoL(jl,sk,qmin , dy1,dy2,dy3,iret)
              enddo
            do k=1,6
              if(mod(k,2).eq.1)akk=ki+5*(k/2)  !nodes
              if(mod(k,2).eq.0)akk=ki+5*(k/2-1)+0.5 
              sk=spminx*(epmax/2./spminx)**((akk-1)/19.) 
            do jdis=0,0 ! PLEASE dont change this without really checking
              if(jdis.eq.0)then
                qmax1=sk/4
                qmax2=sk/4
              elseif(jdis.lt.0)then
                qmax2=q2mnval(maxq2mn)
                qmax1=qmax2/2
              else
                qmax1=sk/4
                qmax2=sk
              endif
            do j=1,3
              qj=qmin*(qmax2/qmin)**((j-1)/3.)
            do i=1,3
              qi=qmin*(qmax1/qmin)**((i-1)/3.)
             qqmax=qmax
              if(jdis.eq.0)then
                qqmin=max(qi,qj)
                q1=qj
                q2=qi
              elseif(jdis.lt.0)then
                qqmin=max(qi,qj)
                q1=min(qi,qj)
                q2=q1
              else
                qqmin=max(qi,qj/4.)
                q1=qj
                q2=qi
              endif
            do lq=1,1  !,3  >1 should not be used 
              qq=qqmin*(qqmax/qqmin)**((lq-1)/3.)  
              q2cmin(1)=qq
              q2cmin(2)=qq
              a=psborn(q1,q2,qq,sk,j1,l1,jdis,0)
              b=psboi(q1,q2,qq,sk,j2,l2,jdis)
              diff=0
              if(abs(a).gt.0)diff=abs(a-b)/(abs(a))
              if(a.lt.1e-20.and.b.lt.1e-20)diff=0
              ncheck=ncheck+1
              diffsum=diffsum+diff
              if(diff.gt.0.1)then
                nshow=nshow+1 
                write (ifmt,'(2i3,e11.3,3f11.1,2x,2e11.3,i5)')jl,jdis,sk
     .          ,q1,q2,qq,a,b,nint(diff*100)
              endif
            enddo
            enddo
            enddo
            enddo
            enddo
            enddo 
            write(ifmt,'(1x,a,i4,a,i6,a,f7.2,a)')
     .      'Summary: shown / checked = ',nshow
     .      ,' / ',ncheck,'   <diff> =',diffsum/ncheck*100.,' %'
            if(ipri.eq.1)stop'TEST psborn'
          !~~~~~~~~~~~~~~~~~~~~~~~~~


      open(1,file=fnii(1:nfnii),status='unknown')
      write (1,*)qcdlam,q2zmin,q2ini,naflav,facpdf4,facpdf5,epmax,pt2cut
     .          ,qcmass,qbmass,nbflav,factjpsi,q2mnval(maxq2mn)
      write (1,*)csbor,csbor0  !KW1808   ,cschar,csbott !bg
      close(1)
      write(ifmt,'(//2a//)')'     written to ',fnii(1:nfnii)
      
      if(.not.negjdis)write(ifmt,'(//a/a/a//)')
     .  '      the (usually) following csord table making takes time'
     . ,'      ==> faster to run several jobs via batch'
     . ,'      use:  set  nq2mnfix -1 in optns'
      
  1   continue
     
      if(negjdis)goto 154  !skip calculation of csord if only change in factjpsi

      !----------------------------------------------------------------
      !
      !  aa  --   CS tabulation    (csord)
      !
      !----------------------------------------------------------------
      
      if(nq2mnfix.eq.0.or.nq2mnfix.eq.99)write(ifmt,'(a,$)')
     .'load aa tables '

      !------------------------------------!
      ! the following allows to generate   ! 
      ! the aa tables in parallel, running !
      ! mxpsar3 batch jobs                 !
      !        use nq2mnfix = -1           !  
      !------------------------------------!
      mx1=mxpsar3 !flavor values
      ixA1=1
      ixA2=mx1
      !~~~~~~~~~~~~~TEST~~~~~~~~~~~~
      !ixA1=9  
      !ixA2=9
      !~~~~~~~~~~~~~TEST~~~~~~~~~~~~
      if(nq2mnfix.eq.-1)then  
        if(iopcnt.le.0)stop'***** nq2mnfix < 0 requires batch *****'
        ix=iopcnt
        ixA1=ix
        ixA2=ix
        write(ifmt,'(//a,i3//)')'consider table for flav',ixA1
      endif

      fnii=fnii(1:index(fnii,'bb.i')-1)//'aa.i ' !remove 00   
      nfnii=index(fnii,' ')-1
      ch2='00'

      do mlflav=ixA1,ixA2          ! ========= mlflav loop =========

      ix=mlflav
      if(ix.lt.10)write(ch2(2:2),'(i1)')ix
      if(ix.ge.10)write(ch2(1:2),'(i2)')ix
      if(ix.ge.100)stop'ERROR 01082019c'

      call makeSudx

      if(nq2mnfix.eq.0.or.nq2mnfix.eq.99)then
        write(ifmt,'(a,$)')'.'
      else
        write(ifmt,'(2a)')'load '
     .  ,fnii(1:nfnii-3)//'a'//ch2//'.i'
      endif
      inquire(file=fnii(1:nfnii-3)//'a'//ch2//'.i',exist=lcalc)

      negjdis=.false.    !compute jdis<0 only
      if(lcalc)then
        if(nq2mnfix.eq.-1)stop'table exists, why nq2mnfix -1 ?'
        negjdis=.false.
        open(1,file=fnii(1:nfnii-3)//'a'//ch2//'.i',status='old')
        read (1,*)qcdlam0,q2zmin0,q2ini0,naflav0,epmax0,pt2cut0,qcmass0 !bg
     *       ,qbmass0,nbflav0,factjpsi0,q2mnval0
        if(qcdlam0.ne.qcdlam)write(ifmt,'(a)')'table aa: wrong qcdlam'
        if(nequal(q2zmin0,q2zmin))
     .   write(ifmt,'(a)')'table aa: wrong q2zmin'
        if(q2ini0 .ne.q2ini )write(ifmt,'(a)')'table aa: wrong q2ini'
        if(naflav0.ne.naflav)write(ifmt,'(a)')'table aa: wrong naflav'
        if(nbflav0.ne.nbflav)write(ifmt,'(a)')'table aa: wrong nbflav'
        if(epmax0 .ne.epmax )write(ifmt,'(a)')'table aa: wrong epmax'
        if(pt2cut0.ne.pt2cut)write(ifmt,'(a)')'table aa: wrong pt2cut'
        if(qcmass0.ne.qcmass)write(ifmt,'(a)')'table aa: wrong qcmass'
        if(qbmass0.ne.qbmass)write(ifmt,'(a)')'table aa: wrong qbmass' !bg
c        if((factjpsi0.ne.factjpsi.or.q2mnval(maxq2mn).ne.q2mnval0)
c     *                                         .and.isetcs.ge.0)then
c          write(ifmt,'(a)')'table aa: change for jdis<0, update table'
c          negjdis=.true.
c        endif
        if(qcdlam0.ne.qcdlam.or.nequal(q2zmin0,q2zmin)
     *  .or.q2ini0 .ne.q2ini
     *  .or.naflav0.ne.naflav.or.epmax0 .ne.epmax.or. pt2cut.ne.pt2cut0
     *  .or.qcmass0.ne.qcmass.or.qbmass0.ne.qbmass.or.nbflav0.ne.nbflav)!bg
     *  then
          write(ifmt,'(//a//)')'   aa table has to be reinitialized!!'
          stop
        endif
        read (1,*)csordr
        close(1)
        do j=1,20
        do i=1,20            
        do k=1,20
          ml=mlflav  
          k1=k+20*(ml-1) 
          k2=k1+mxpsar1
          csord(i,j,k1)=csordr( i,j,k,1)
          csord(i,j,k2)=csordr( i,j,k,2)
        enddo
        enddo
        enddo
        !~~~~~~~~~~TEST~~~~~~~~~~~~~~~
        ipril=999
        if(ipril.eq.141)goto 151 
        !~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
        if(.not.negjdis)goto 5
      elseif(.not.producetables)then
        write(ifmt,*) "Missing epos aa file !"        
        write(ifmt,*) "Please correct the defined path ",
     &"or force production ..."
        stop
      endif

      write(ifmt,'(2a/a)')fnii(1:nfnii-3)//'a'//ch2//'.i'
     .,' does not exist ','         -> calculate CS tables aa ...'
      if(nq2mnfix.eq.0)then
        write(ifmt,*)'  '
        write(ifmt,*)'  ==> Better use nq2mnfix -1 '
        write(ifmt,*)'  ==> to create aa tables in parallel'
        write(ifmt,*)'  ==> by running',mxpsar3,'  batch jobs'
        stop
      endif
 
      ipri=0
      goto 152
 151  continue
      ipri=1
 152  continue

      nckeck=0
      nshow=0
      diffsum=0

      ml=mlflav  
        call eolpi( ml , m1 , l1 ) !KW1808     
        call getQminEoL(ml,qref1,qref2,qqref)
        call getSminEoL(ml,qqref,spminx)
        spminx=spminx*1.10

      write (ifmt,*)'-----------------------------------'
      write (ifmt,'(a,i3)')' ordered jet xsection csord, ml =',ml
      write (ifmt,*)'-----------------------------------'

      do k=1,20
        sk=spminx*(epmax/2./spminx)**((k-1)/19.)
        write (ifmt,'(a,i4,a,i4,a,i4,a)')'   ml =',ml,'  /',mxpsar3
     .  ,'   k =',k,'  /  20'
        call clop(3)
        call getQmaxEoL(ml,sk,qmax)
        qmax=qmax / 1.01!1.10 
        qmax=max(qmax,1.01*qqref)
      do i=1,20             !cross-sections initialization
        qi=qref1*(qmax/qref1)**((i-1)/19.)
        qii=max(qi,qref2)
      do j=1,20
        qq=qii*(qmax/qii)**((j-1)/19.)  
        k1=k+20*(ml-1) 
        k2=k1+mxpsar1
        if(k.eq.1.or.i.eq.20.or.j.eq.20)then   !tp: qi > qq Is case qq < qi properly done ??????
          csipo1=psborn(qi,qq,qq,sk,m1,l1,0,0)
          csipo2=psborn(4.*qq,qi,qq,sk,m1,l1,1,0)
        else
          csipo1=psjet1(qi,qq,qq,sk,m1,l1,0)
     *             + psborn(qi,qq,qq,sk,m1,l1,0,0) 
          csipo2=psjet1(qi,4.*qq,qq,sk,m1,l1,2)
     *             +psborn(4.*qq,qi,qq,sk,m1,l1,1,0)      
        endif
        !~~~~~~~~~~~~~~~~~~~~~~
        csord(i,j,k1)=log(max(1.e-30,csipo1))
        csord(i,j,k2)=log(max(1.e-30,csipo2))
        !~~~~~~~~~~~~~~~~~~~~~~
        csordr( i,j,k,1)=csord(i,j,k1)
        csordr( i,j,k,2)=csord(i,j,k2)
cc      if(i.eq.18.and.j.le.20.and.k1.eq.161)print*,'+++csord+++'
cc     *,qi,qq,sk,i,j,k1,csipo1,csord(i,j,k1)  
      enddo
      enddo
      enddo

      open(1,file=fnii(1:nfnii-3)//'a'//ch2//'.i',status='unknown')
      write (1,*)qcdlam,q2zmin,q2ini,naflav,epmax,pt2cut,qcmass,qbmass
     .          ,nbflav,factjpsi,q2mnval(maxq2mn)
      write (1,*)csordr 
      close(1)
      write(ifmt,'(2a)')'table '
     .,fnii(1:nfnii-3)//'a'//ch2//'.i  created'

  5   continue

            if(abs(nq2mnfix).eq.1)then
            !~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
            write (ifmt,*)'tests:'   !Ordered
            !~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

            ml=mlflav  
              call clop(3)
              write (ifmt,'(a,a,a)')'  ml        sk   '
     .        ,'      qi         qj         qq  '
     .        ,'        ord      ord-i  '
              call eolpi( ml , m1 , l1 ) 
              call eolpic( m1 , l1 , m2 , l2 ) 
              qqref=q2zmin
              call getSminEoL(ml,qqref,spminx)
              qmin=1.5*q2zmin
              ki=0
              iret=1
              do while(iret.eq.1)
                ki=ki+1
                sk=spminx*(epmax/2./spminx)**((ki-1)/19.) 
                call getLimitsEoL(ml,sk,qmin , dy1,dy2,dy3,iret)
              enddo
            do k=1,6
              if(mod(k,2).eq.1)akk=ki+5*(k/2)  !nodes
              if(mod(k,2).eq.0)akk=ki+5*(k/2-1)+0.5 
              sk=spminx*(epmax/2./spminx)**((akk-1)/19.) 
            do n=1,1
              qmax1=sk/4
              qmax2=sk/4
              if(n.eq.2)qmax2=qmax2*4
            do i=1,4
              if(i.eq.1)qi=2
              if(i.eq.2)qi=6.5
              if(i.eq.3)qi=25
              if(i.eq.4)qi=100
            do j=1,4
              q2cmin(1)=4.
              if(j.eq.1)qj=1.5
              if(j.eq.2)qj=6.5
              if(j.eq.3)qj=25
              if(j.eq.4)qj=100
              qqmax=qmax1
              if(n.eq.1)then
                qqmin=max(qi,qj)
              else
                qqmin=max(qi,qj/4.)
              endif
              qq=qqmin
              a=psjet1(qi,qj,qq,sk,m1,l1,2*(n-1))
     .         +psborn(qi,qj,qq,sk,m1,l1,n-1,0)
              b=psjti1(qi,qj,qq,sk,m2,l2,n-1)
              diff=0
              if(abs(a).gt.0)diff=abs(a-b)/(abs(a))
              if(a.lt.1e-20.and.b.lt.1e-20)diff=0
              ncheck=ncheck+1
              diffsum=diffsum+diff
              if(diff.gt.0.1)then
                !if(ml.eq.7)then!TEST (instead of if(diff...) 
                nshow=nshow+1 
                write (ifmt,'(i3,e11.3,3f11.1,2x,2e11.3,i5)')
     .          ml,sk,qi,qj,qq,a,b,nint(diff*100)
              endif   
            enddo
            enddo
            enddo
            enddo
            call clop(3)
            !~~~~~~~~~~~~~~~~~~~~~~~~~
            endif
         
      enddo                     ! ========= mlflav loop =========

      if(abs(nq2mnfix).eq.1)then
        write(ifmt,'(1x,a,i4,a,i6,a,f7.2,a)')
     .  'Summary: shown / checked = ',nshow
     .  ,' / ',ncheck,'   <diff> =',diffsum/ncheck*100.,' %'
        if(ipri.eq.1)stop'TEST psjet1'
        if(nq2mnfix.eq.-1)then
          write(ifmt,'(//5x,a,$)')'Table making finished => stop'
          write(ifmt,'(a//)')'  ==> set nq2mnfix 1  for testing !!!'
          stop'Table making finished'
        endif
      endif

      if(nq2mnfix.eq.0.or.nq2mnfix.eq.99)write(ifmt,*)' '

 154  continue

      !----------------------------------------------------------------
      !
      !  ab  --   CS tabulation    (cstot)
      !
      !----------------------------------------------------------------

      if(nq2mnfix.eq.0.or.nq2mnfix.eq.99)write(ifmt,'(a,$)')
     .'load ab tables '

      !------------------------------------!
      ! the following allows to generate   ! 
      ! the ab tables in parallel, running !
      ! mxpsar3 * mx2  batch jobs           ! 
      !    use nq2mnfix = -2               !
      !------------------------------------!
      mx1=mxpsar3 !flavor values
      mx2=20 !energy values
      mx2=mx2/kenwidth !blocks
      ixA1=1
      ixA2=mx1
      ixB1=1
      ixB2=mx2
      !~~~~~~~~~~~~~TEST~~~~~~~~~~~~
      !ixA1=9  
      !ixA2=9
      !ixB1=1
      !ixB2=1
      !~~~~~~~~~~~~~TEST~~~~~~~~~~~~
      if(nq2mnfix.eq.-2)then  !useful for making tables in independent runs
        if(iopcnt.le.0)stop'***** iopcnt<0 batch required *****'
        ix=iopcnt
        ixA1=(ix-1)/mx2+1
        ixA2=ixA1
        ixB1=mod(ix-1,mx2)+1
        ixB2=ixB1
        write(ifmt,'(//a,i3,a,i3//)')'consider table for flav'
     .  ,ixA1,'  energy block',ixB1
      elseif(nq2mnfix.eq.2)then  !useful for testing in independent runs
        if(iopcnt.le.0)stop'***** nq2mnfix=2 batch required *****'
        ix=iopcnt
        ixA1=ix
        ixA2=ix
        ixB1=1
        ixB2=mx2
        if(ix.gt.mxpsar3)stop'batch counter `iopcnt` too big'
        write(ifmt,'(//a,i3,a,i3,a,i3//)')'consider table for flav'
     .  ,ixA1,'  energy blocks',ixB1,' to ',ixB2
      endif
      ch4='0000'

      do mlflav=ixA1,ixA2  ! ========= mlflav loop =========
      do kblock=ixB1,ixB2    ! ========= energ block loop =========

      ix=(mlflav-1)*mx2+kblock
      if(ix.lt.10)write(ch4(4:4),'(i1)')ix
      if(ix.ge.10)write(ch4(3:4),'(i2)')ix
      if(ix.ge.100)write(ch4(2:4),'(i3)')ix
      if(ix.ge.1000)write(ch4(1:4),'(i4)')ix
      if(ix.ge.10000)stop'ERROR 01082019b'

      call makeSudx

      if(ish.ge.4)write(ifch,*)'bare cross sections ...'

      if(nq2mnfix.eq.0.or.nq2mnfix.eq.99)then
        write(ifmt,'(a,$)')'.'
      else
        write(ifmt,'(2a)')'load '
     .  ,fnii(1:nfnii-3)//'b'//ch4//'.i'
      endif
      inquire(file=fnii(1:nfnii-3)//'b'//ch4//'.i'
     . ,exist=lcalc)

      if(lcalc)then
        if(nq2mnfix.eq.-2)stop'file exists, nq2mnfix should be >= 0'
        open(1,file=fnii(1:nfnii-3)//'b'//ch4//'.i'
     .   ,status='old')
        read (1,*)qcdlam0,q2jmin0,q2ini0,naflav0,epmax0,pt2cut0,qcmass0 !bg
     *  ,qbmass0,nbflav0
        if(qcdlam0.ne.qcdlam)write(ifmt,'(a)')'table ab: wrong qcdlam'
        if(nequal(q2jmin0,q2jmin))
     .  write(ifmt,'(a)')'table ab: wrong q2jmin'
        if(q2ini0 .ne.q2ini )write(ifmt,'(a)')'table ab: wrong q2ini'
        if(naflav0.ne.naflav)write(ifmt,'(a)')'table ab: wrong naflav'
        if(nbflav0.ne.nbflav)write(ifmt,'(a)')'table ab: wrong nbflav'
        if(epmax0 .ne.epmax )write(ifmt,'(a)')'table ab: wrong epmax'
        if(pt2cut0.ne.pt2cut)write(ifmt,'(a)')'table ab: wrong pt2cut'
        if(qcmass0.ne.qcmass)write(ifmt,'(a)')'table ab: wrong qcmass'
        if(qbmass0.ne.qbmass)write(ifmt,'(a)')'table ab: wrong qbmass' !bg
        if(qcdlam0.ne.qcdlam.or.nequal(q2jmin0,q2jmin)
     *  .or.q2ini0 .ne.q2ini
     *  .or.naflav0.ne.naflav.or.epmax0 .ne.epmax.or. pt2cut.ne.pt2cut0
     *  .or.qcmass0.ne.qcmass.or.qbmass0.ne.qbmass.or.nbflav0.ne.nbflav)     !bg
     *  then
          write(ifmt,'(//a//)')'   ab table has to be reinitialized!!!'
          stop
        endif
        read (1,*)cstotr
        close(1)
        do i=1,20     
        do j=1,20
        do n=1,2
        do ken=1,kenwidth
          k=kenwidth*(kblock-1)+ken
          cstot(i,j,k,mlflav,n)=cstotr(i,j,n,ken)
        enddo
        enddo
        enddo
        enddo
        goto 2
      elseif(.not.producetables)then
        write(ifmt,*) "Missing epos ab file !"        
        write(ifmt,*) "Please correct the defined path ",
     &"or force production ..."
        stop
      endif

      write(ifmt,'(2a/2a,f7.2,a)')
     . fnii(1:nfnii-3)//'b'//ch4//'.i'
     .,' does not exist ','     -> calculate CS table ab'//ch4
      if(nq2mnfix.eq.0)then
        write(ifmt,*)'  '
        write(ifmt,*)'  ==> Better use nq2mnfix -2 '
        write(ifmt,*)'  ==> to create ab tables in parallel'
        write(ifmt,*)'  ==> by running',mxpsar3*mx2,'  batch jobs'
        stop
      endif

      ml=mlflav
        call eolpi( ml , m1 , l1 ) !KW1808c     
        call eolpic( m1 , l1 , m2 , l2 ) 
        call getQminEoL(ml,qref1,qref2,qqref)
        call getSminEoL(ml,qqref,spmin)
        spmin=spmin*1.10
      do ken=1,kenwidth
        k=kenwidth*(kblock-1)+ken
        write(ifmt,'(a,i4,a,i4,a,i4,a,f8.3)')'   ml =',ml,'  /',mxpsar3
     .  ,'   k =',k,'  /  20  '
        call clop(3)
        sk=spmin*(epmax/2./spmin)**((k-1)/19.)  !c.m. energy squared for the hard
        call getQmaxEoL(ml,sk,qmax)
        qmax=qmax / 1.01!1.10 
        qmax=max(qmax,q2jmin)
        qmax=max(qmax,1.01*qqref)
      do n=1,1  
      do i=1,20     
        qi=qref1*(qmax/qref1)**((i-1)/19.)
      do j=1,20
        qj=qref2*(qmax/qref2)**((j-1)/19.)
        qq=max(qi,qj) 
        !~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ 
        if(k.eq.1.or.i.eq.20.or.j.eq.20)then
          a=psborn(qi,qj,qq,sk,m1,l1,n-1,0)
          cstot(i,j,k,ml,n)=log(max(1.e-30,a)) !KW1808c 
        else
          if(n.eq.1)then
            a=psjet(qi,qj,qq,sk,m1,l1,0)+
     *       psjti1(qi,qj,qq,sk,m2,l2,0)+
     *       psjti1(qj,qi,qq,sk,l2,m2,0)
     *       -psboi(qi,qj,qq,sk,m2,l2,0)   
          else
            a=psjet(qi,qj,qq,sk,m1,l1,1)+
     *       psjet1(qi,qj,qq,sk,m1,l1,1)+
     *       psjti1(qj,qi,qq,sk,l2,m2,1)
          endif
          cstot(i,j,k,ml,n)=log( max( 1.e-30 , a ) )
        endif
        !~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ 
        if(.not.(k.eq.1.or.i.eq.20.or.j.eq.20))then
          a=cstot(i,j,k,ml,n)
          if(.not.(a.gt.0..or.a.le.0.) !NaN catch
     .    .or.a.gt.1e35 .or. a.lt.-1e35 )stop'ERROR 01082019'
        endif      
        cstotr(i,j,n,ken)=cstot(i,j,k,ml,n)
      enddo
      enddo
      enddo
      enddo

      open(1,file=fnii(1:nfnii-3)//'b'//ch4//'.i'
     . ,status='unknown')
      write (1,*)qcdlam,q2jmin,q2ini,naflav,epmax,pt2cut,qcmass,qbmass
     .          ,nbflav !bg
      write (1,*)cstotr
      close(1)
      write(ifmt,'(2a)')'table '
     .,fnii(1:nfnii-3)//'b'//ch4//'.i  created'
      
            !~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
            write (ifmt,*)'tests:    (partial ab tables)'
            !~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

            !~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
            write (ifmt,*)'only preliminary tests (between nodes), '
     .      ,'since tables are not yet complete'
            write (ifmt,*)'  use `set nq2mnfix 2` for final tests, '
     .      ,'running',mxpsar3,' batch jobs'
            !~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
            write (ifmt,'(a,a,a)')' jl       sk   '
     .      ,'        qi         qj         qq  '
     .      ,'       cs     cs-i    error(%)'
            ncheck=0
            nshow=0
            diffsum=0

            jl=mlflav
              ml=jl
              call eolpi( jl , j1 , l1 ) 
              call eolpic( j1 , l1 , j2 , l2 ) 
              call getQminEoL(ml,qref1,qref2,qqref)
              call getSminEoL(ml,qqref,spmin)
              spmin=spmin*1.10
            do jen=2,4
              yen=jen/2.
              yk=kenwidth*(kblock-1)+yen
              sk=spmin*(epmax/2./spmin)**((yk-1)/19.)  
              call getQmaxEoL(jl,sk,qmax1)
            do i=1,4
              if(i.eq.1)qi=2
              if(i.eq.2)qi=6.5
              if(i.eq.3)qi=25
              if(i.eq.4)qi=100
            do j=1,3
              if(j.eq.1)qj=2
              if(j.eq.2)qj=6.5
              if(j.eq.3)qj=25
              qq=max(qi,qj)
              if(log(sk/spmin)/log(epmax/2./spmin)*19.+1. .lt.1.5)then
                a=psborn(qi,qj,qq,sk,j1,l1,0,0)
                !a2=a
              else
                a=psjet(qi,qj,qq,sk,j1,l1,0)+
     .          psjti1(qi,qj,qq,sk,j2,l2,0)+
     .          psjti1(qj,qi,qq,sk,l2,j2,0)
     .          -psboi(qi,qj,qq,sk,j2,l2,0)
                !a2=psjet(qi,qj,qq,sk,j1,l1,0)+
                !.psjet1(qi,qj,qq,sk,j1,l1,0)+
                !.psjet1(qj,qi,qq,sk,l1,j1,0)+
                !.psborn(qi,qj,qq,sk,j1,l1,0,0)
              endif 
              b=psjtoi(qi,qj,qq,sk,j2,l2,0)
              diff=0
              if(abs(a).gt.0)diff=abs(a-b)/(abs(a))
              if(a.lt.1e-20.and.b.lt.1e-20)diff=0
              diff=min(diff,1.)
              ncheck=ncheck+1
              diffsum=diffsum+diff
              !if(diff.gt.0.1)then
              nshow=nshow+1 
              write (ifmt,'(i3,e11.3,3f11.1,2x,2e11.3,i5,2f9.3)')
     .        jl,sk,qi,qj,qq,a,b,nint(diff*100)
              !endif
            enddo
            enddo
            enddo
            write(ifmt,'(a/10x,a,i4,a,i6,a,f7.2,a)')
     .      'Summary (based on partial table only !!!)'
     .      ,' shown / checked = ',nshow
     .      ,' / ',ncheck,'   <diff> =',diffsum/ncheck*100.,' %'
            !~~~~~~~~~~~~~~~~~~~~~~~~~

      write(ifmt,'(//5x,a,$)')'Table making finished => stop'
      write(ifmt,'(a//)')'  ==> set nq2mnfix 1  for testing !!!'
      stop'Table making finished'

 2    continue ! <---- get here after reading table

      enddo                 ! ======== kblock loop ========

            if(abs(nq2mnfix).eq.2)then
            !~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
            write (ifmt,*)'tests: (ab tables) jl =',mlflav,' / 20 '
            !~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
            write (ifmt,'(a,a,a)')' jl       sk   '
     .      ,'        qi         qj         qq  '
     .      ,'       cs     cs-i    error(%)'
            ncheck=0
            nshow=0
            diffsum=0

            jl=mlflav
              ml=jl
              call eolpi( jl , j1 , l1 ) 
              call eolpic( j1 , l1 , j2 , l2 ) 
              call getQminEoL(ml,qref1,qref2,qqref)
              call getSminEoL(ml,qqref,spmin)
              spmin=spmin*1.10
            do kblock=ixB1,ixB2            
            do jen=2,3 !5
              yen=jen/2.
              yk=kenwidth*(kblock-1)+yen
              sk=spmin*(epmax/2./spmin)**((yk-1)/19.)  
              call getQmaxEoL(jl,sk,qmax1)
            do i=1,4
              if(i.eq.1)qi=2
              if(i.eq.2)qi=6.5
              if(i.eq.3)qi=25
              if(i.eq.4)qi=100
            do j=1,3
              if(j.eq.1)qj=2
              if(j.eq.2)qj=6.5
              if(j.eq.3)qj=25
              qq=max(qi,qj)
              if(log(sk/spmin)/log(epmax/2./spmin)*19.+1. .lt.1.5)then
                a=psborn(qi,qj,qq,sk,j1,l1,0,0)
                !a2=a
              else
                a=psjet(qi,qj,qq,sk,j1,l1,0)+
     .          psjti1(qi,qj,qq,sk,j2,l2,0)+
     .          psjti1(qj,qi,qq,sk,l2,j2,0)
     .          -psboi(qi,qj,qq,sk,j2,l2,0)
                !a2=psjet(qi,qj,qq,sk,j1,l1,0)+
                !.psjet1(qi,qj,qq,sk,j1,l1,0)+
                !.psjet1(qj,qi,qq,sk,l1,j1,0)+
                !.psborn(qi,qj,qq,sk,j1,l1,0,0)
              endif 
              b=psjtoi(qi,qj,qq,sk,j2,l2,0)
              diff=0
              if(abs(a).gt.0)diff=abs(a-b)/(abs(a))
              if(a.lt.1e-20.and.b.lt.1e-20)diff=0
              diff=min(diff,1.)
              ncheck=ncheck+1
              diffsum=diffsum+diff
              if(diff.gt.0.1)then
              nshow=nshow+1 
              write (ifmt,'(i3,e11.3,3f11.1,2x,2e11.3,i5,2f9.3)')
     .        jl,sk,qi,qj,qq,a,b,nint(diff*100)
              endif
            enddo
            enddo
            enddo
            enddo
            write(ifmt,'(1x,a,i4,a,i6,a,f7.2,a)')
     .      'Summary:  shown / checked = ',nshow
     .      ,' / ',ncheck,'   <diff> =',diffsum/ncheck*100.,' %'
            !~~~~~~~~~~~~~~~~~~~~~~~~~
            endif

      enddo                 ! ======== mlflav loop ========
      
      if(nq2mnfix.eq.2)stop

      if(nq2mnfix.eq.0.or.nq2mnfix.eq.99)write(ifmt,*)' '


      !----------------------------------------------------------------
      !
      !  ac  --   CSzero tabulation    (cstotzero,csborzer)
      !
      !           loop   q 2 m i n 1   and   q 2 m i n 2    
      !                                            (includes x dep)
      !  Tables depend on q2zmin only for binning.
      !----------------------------------------------------------------

      write(ifmt,'(a,$)')'load ac tables '

      !--------------------------------!
      ! generating ac tables is fast   ! 
      ! can be done interactively      !
      !        use nq2mnfix = 0        !  
      !--------------------------------!

      nq2mn1=1
      nq2mn2=maxq2mx
      if(iscreen.eq.0)then      !maximum q2min given by q2nmin
        q=log(max(1e-6,q2nmin))
        call ipoli3(q,nq,frac)
        nq2mn2=nq+2
      endif
      
      do nq2mn=nq2mn1,nq2mn2  ! ========= q2jmin loop =========
      do nq2xx=nq2mn1,nq2mn2  ! ========= q2lmin loop =========

      q2lmin=min( q2mnval(nq2mn) , q2mnval(nq2xx) )       !????????????tp
      q2jmin=max( q2mnval(nq2mn) , q2mnval(nq2xx) )       !????????????tp
      q2cmin(1)=q2jmin
      q2cmin(2)=q2lmin
      
      if(ish.ge.4)write(ifch,*)'bare cross sections ...'


      inquire(file=
     . fnii(1:nfnii-3)//'c'//q2mnch(nq2mn)//q2mnch(nq2xx)//'.i',
     . exist=lcalc)
      if(lcalc)then
        !write(ifmt,'(2a)')'load '
        !. ,fnii(1:nfnii-3)//'c'//q2mnch(nq2mn)//q2mnch(nq2xx)//'.i'
        write(ifmt,'(a,$)')'.'
        open(1,file=
     .  fnii(1:nfnii-3)//'c'//q2mnch(nq2mn)//q2mnch(nq2xx)//'.i'
     .  ,status='old')
        read (1,*)qcdlam0,q2jmin0,q2lmin0,q2ini0,naflav0,epmax0,pt2cut0
     .   ,qcmass0,qbmass0,nbflav0                              !bg
        if(qcdlam0.ne.qcdlam)write(ifmt,'(a)')'table ac: wrong qcdlam'
        if(nequal(q2jmin0,q2jmin))
     .   write(ifmt,'(a)')'table ac: wrong q2jmin'
        if(nequal(q2lmin0,q2lmin))
     .   write(ifmt,'(a)')'table ac: wrong q2lmin'
        if(q2ini0 .ne.q2ini )write(ifmt,'(a)')'table ac: wrong q2ini'
        if(naflav0.ne.naflav)write(ifmt,'(a)')'table ac: wrong naflav'
        if(nbflav0.ne.nbflav)write(ifmt,'(a)')'table ac: wrong nbflav'
        if(epmax0 .ne.epmax )write(ifmt,'(a)')'table ac: wrong epmax'
        if(pt2cut0.ne.pt2cut)write(ifmt,'(a)')'table ac: wrong pt2cut'
        if(qcmass0.ne.qcmass)write(ifmt,'(a)')'table ac: wrong qcmass'
        if(qbmass0.ne.qbmass)write(ifmt,'(a)')'table ac: wrong qbmass' !bg
        if(qcdlam0.ne.qcdlam.or.nequal(q2jmin0,q2jmin)
     *  .or.q2ini0 .ne.q2ini.or.nequal(q2lmin0,q2lmin)
     *  .or.naflav0.ne.naflav.or.epmax0 .ne.epmax.or. pt2cut.ne.pt2cut0
     *  .or.qcmass0.ne.qcmass.or.qbmass0.ne.qbmass.or.nbflav0.ne.nbflav)
     *  then
          write(ifmt,'(//a//)')'   ac table has to be reinitialized!!!'
          stop
        endif
        read (1,*)cstotzero,csborzer      
        close(1)
        goto 3
      elseif(.not.producetables)then
        write(ifmt,*) "Missing epos ac file !"        
        write(ifmt,*) "Please correct the defined path ",
     &"or force production ..."
        stop
      endif

      write(ifmt,'(2a/2a,f7.2,a,f7.2,a)')
     . fnii(1:nfnii-3)//'c'//q2mnch(nq2mn)//q2mnch(nq2xx)//'.i'
     .,' does not exist ','     -> calculate CS tables ac '
     .,'for q2jmin =',q2jmin,'  q2lmin =',q2lmin,' ...'


      !total and born hard cross-sections logarithms for minimal cutoffs
      ! (q2jmin,q2lmin), interpolated in the psjti0 procedure
      
      do ml=1,mxpsar3 !KW1808c  
        call eolpi( ml , m1 , l1 ) !KW1808c     
        call eolpic( m1 , l1 , m2 , l2 ) 
        qqref=max(q2jmin,q2lmin)
        call getSminEoL(ml,qqref,spminx)
      do k=1,20
        sk=spminx*(epmax/2./spminx)**((k-1)/19.)  !c.m. energy squared for the hard
        qq=max(q2jmin,q2lmin)
        !~~~~~~~~~~~~~~~~~~~~
        csborzer(k,ml)
     *      =log(max(1.e-30,psborn(q2jmin,q2lmin,qq,sk,m1,l1,0,0))) 
        if(k.eq.1)then
          cstotzero(k,ml)=csborzer(k,ml)
        else 
          cstotzero(k,ml)
     .    =log(max(1.e-30,psjtoi(q2jmin,q2lmin,qq,sk,m2,l2,0)))
        endif
        !~~~~~~~~~~~~~~~~~~~~~
      enddo
      enddo

      open(1,file=
     . fnii(1:nfnii-3)//'c'//q2mnch(nq2mn)//q2mnch(nq2xx)//'.i'
     . ,status='unknown')
      write (1,*)qcdlam,q2jmin,q2lmin,q2ini,naflav,epmax,pt2cut,qcmass
     * ,qbmass,nbflav                               !bg
      write (1,*)cstotzero,csborzer
      close(1)
      write(ifmt,'(2a)')'     written to '
     .,fnii(1:nfnii-3)//'c'//q2mnch(nq2mn)//q2mnch(nq2xx)//'.i'
      
3     continue
      
      do k=1,20 
         do ml=1,mxpsar3 !KW1808c
            csqtotzero(nq2mn,nq2xx,k,ml)=cstotzero(k,ml)
            csqborzer( nq2mn,nq2xx,k,ml)=csborzer( k,ml)
         enddo
      enddo

            
      enddo ! ======== q2jmin loop ========
      enddo ! ======== q2lmin loop ========
      
      write(ifmt,'(a)')' '

      !call testCsExIpo(1,1,, 0)

      !call testOmEx(5,1,0.1,0.1, 1)
      !call testOmEx(1 ,5,0.1,0.1, 1)
      !call testOmEx(1,1,0.5,0.5, 0)
      call testOmEx(5,5,0.5,0.5, 1)
      !call testOmEx(5,5,0.9,0.9, 1)


      !-----------------------------------------------------------------
      !
      !  rj  --  tabulation of om5  (fhss,fhgg,fhqg,fhgq,fhqq)
      !
      !                using om11pp, om11pp3 and om11pp4
      !                      ======  =======     =======
      !
      !          loop     q 2 m i n 1  and   q 2 m i n 2
      !                                            (without x dep)
      !  Tables depend on q2zmin only for binning.
      !-----------------------------------------------------------------

      call pdfparamget(2,facpdf4sat) 
      call pdfparamget(4,facpdf5sat) 

      if(isetcs.ge.0)then      !possibility to do some test without tabulated rj

      write(ifmt,'(a,$)')'load rj tables '

      !------------------------------------!
      ! the following allows to generate   ! 
      ! the rj tables in parallel, running !
      ! maxq2mn * maxq2mn  batch jobs           ! 
      !    use nq2mnfix = -3               !
      !------------------------------------!
      nq2mnA1=1
      nq2mnA2=maxq2mx
      nq2mnB1=1
      nq2mnB2=maxq2mx
      if(iscreen.eq.0)then      !maximum q2min given by q2nmin
        q=log(max(1e-6,q2nmin))
        call ipoli3(q,nq,frac)
        nq2mnA2=nq+2
        nq2mnB2=nq+2
      endif
      if(nq2mnfix.eq.-3)then  !useful for making tables in independent runs
       if(iopcnt.le.0)stop'\n\n batch required \n\n'
       nq2mnfx=iopcnt-1
       nq2mnA1=nq2mnfx/maxq2mn+1
       nq2mnA2=nq2mnA1
       nq2mnB1=mod(nq2mnfx,maxq2mn)+1
       nq2mnB2=nq2mnB1
       write(ifmt,'(//a,i3,a,i3//)')'consider only table for nq2mn ='
     .      ,nq2mnA2,'  nq2xx =',nq2mnB2
      endif

      do nzzvex=1,mxzzvex ! ========= zzverx loop =========
      
      do nq2mn=nq2mnA1,nq2mnA2  ! ========= q2jmin loop =========
      do nq2xx=nq2mnB1,nq2mnB2  ! ========= q2lmin loop =========

      negjdis=.false.
      q2jmin=q2mnval(nq2mn)       
      q2lmin=q2mnval(nq2xx)       
      ix=nfnii-4
      iv=nzzvex
      inquire(file=
     &fnii(1:ix)//'rj'//q2mnch(nq2mn)//q2mnch(nq2xx)//zzvexch(iv)//'.i'
     &,exist=lcalc)
      if(lcalc)then
        write(ifmt,'(a,$)')'.'
        open(1,file=
     &fnii(1:ix)//'rj'//q2mnch(nq2mn)//q2mnch(nq2xx)//zzvexch(iv)//'.i'
     .  ,status='old')
        read (1,*)q2sft0,alppar0,alplea0,alppom0,slopom0
     *  ,gamhad0,r2had0,facpdf40sat,facpdf50sat
     *  ,qcdlam0,q2jmin0,q2lmin0,q2ini0,betpom0,glusea0,naflav0,nbflav0
     *  ,factk0,gampar0,gamtil0,zoeinc0,zzsoft0
        if(q2sft0.ne.q2sft)write(ifmt,'(a,2f8.4)')
     *  'table rj: wrong q2sft',q2sft0,q2sft
        if(alppar0.ne.alpparh)write(ifmt,'(a,2f8.4)')
     *  'table rj: wrong alppar',alppar0,alpparh
        if(alppom0.ne.alppom)write(ifmt,'(a,2f8.4)')
     *  'table rj: wrong alppom',alppom0,alppom
        if(slopom0.ne.slopom)write(ifmt,'(a,2f8.4)')
     *  'table rj: wrong slopom',slopom0,slopom
        iii=2
c        if(gamhad0(iii).ne.gamhad(iii))write(ifmt,'(a,i1,a,2f8.4)')
c     *  'table rj: wrong gamhad(',iii,')',gamhad0(iii),gamhad(iii)
        do iii=1,3
        if(r2had0(iii) .ne.r2had(iii) )write(ifmt,'(a,i1,a,2f8.4)')
     *  'table rj: wrong r2had(',iii,')',r2had0(iii),r2had(iii)
        if(alplea0(iii).ne.alplea0(iii))write(ifmt,'(a,i1,a,2f8.4)')
     *  'table rj: wrong alplea(',iii,')',alplea0(iii),alplea(iii)
        enddo
        if(qcdlam0.ne.qcdlam)write(ifmt,'(a,2f8.4)')
     *  'table rj: wrong qcdlam',qcdlam0,qcdlam
        if(nequal(q2jmin0,q2jmin))write(ifmt,'(a,2f8.4)')
     *  'table rj: wrong q2jmin',q2jmin0,q2jmin
        if(nequal(q2lmin0,q2lmin))write(ifmt,'(a,2f8.4)')
     *  'table rj: wrong q2lmin',q2lmin0,q2lmin
        if(q2ini0 .ne.q2ini )write(ifmt,'(a,2f8.4)')
     *  'table rj: wrong q2ini',q2ini0,q2ini
        if(betpom0.ne.betpomi)write(ifmt,'(a,2f8.4)')
     *  'table rj: wrong betpom',betpom0,betpomi
        if(glusea0.ne.glusea)write(ifmt,'(a,2f8.4)')
     *  'table rj: wrong glusea',glusea0,glusea
        if(naflav0.ne.naflav)write(ifmt,'(a,2i1)')
     *  'table rj: wrong naflav',naflav0,naflav
        if(nbflav0.ne.nbflav)write(ifmt,'(a,2i1)')
     *  'table rj: wrong nbflav',nbflav0,nbflav
        !if(factk0 .ne.factk )write(ifmt,'(a,2f8.4)')   !cKW21 factk is not used in om11pp...
        !*  'table rj: wrong factk', factk0,factk
        if(gamtil0 .ne.gamtil )write(ifmt,'(a,2f8.4)')
     *  'table rj: wrong gamtil', gamtil0,gamtil
        if(zzsoft0 .ne.zzsoft )write(ifmt,'(a,2f8.4)')
     *  'table rj: wrong zzsoft', zzsoft0,zzsoft
        if(zzsoft0.ne.zzsoft.and.isetcs.ge.0)then
        write(ifmt,'(a)')'table rj: change for new zzsoft, update table'
          negjdis=.true.
        endif
        if(zoeinc0 .ne.zoeinc )write(ifmt,'(a,2f8.4)')
     *  'table rj: wrong zoeinc', zoeinc0,zoeinc
        if(zoeinc0.ne.zoeinc.and.isetcs.ge.0)then
        write(ifmt,'(a)')'table rj: change for new zoeinc, update table'
          negjdis=.true.
        endif
        if(alppar0.ne.alpparh.or.alppom0.ne.alppom
     *  .or.slopom0.ne.slopom.or.q2sft0.ne.q2sft
     *  .or.r2had0(1).ne.r2had(1).or.r2had0(2).ne.r2had(2)
     *  .or.r2had0(3).ne.r2had(3)
     *  .or.alplea0(1).ne.alplea(1).or.alplea0(2).ne.alplea(2)
     *  .or.alplea0(3).ne.alplea(3)
     *  .or.qcdlam0.ne.qcdlam.or.nequal(q2jmin0,q2jmin)
     *  .or.nequal(q2lmin0,q2lmin)
     *  .or.q2ini0 .ne.q2ini.or.gamtil0.ne.gamtil
     *  .or.betpom0.ne.betpomi.or.glusea0.ne.glusea.or.naflav0.ne.naflav
cKW21     *  .or.factk0 .ne.factk
     *  .or.nbflav0.ne.nbflav
     *  )then
           write(ifmt,'(//a//)')'   table rj has to be reinitialized!!!'
           stop
        endif
        if(facpdf40sat.ne.facpdf4sat)then
          write(ifmt,'(/80a)')('-',k=1,80)
          write(ifmt,'(a)')'    ERROR Wrong facpdf4sat for rj table '
          write(ifmt,'(7x,f5.2,a)')facpdf40sat,' expected'   
          write(ifmt,'(a)')'      Redo rj tables '
          write(ifmt,'(80a/)')('-',k=1,80)
          stop
        endif 
        if(facpdf50sat.ne.facpdf5sat)then
          write(ifmt,'(/80a)')('-',k=1,80)
          write(ifmt,'(a)')'    ERROR Wrong facpdf5sat for rj table '
          write(ifmt,'(7x,f5.2,a)')facpdf50sat,' expected'   
          write(ifmt,'(a)')'      Redo rj tables '
          write(ifmt,'(80a/)')('-',k=1,80)
          stop
        endif 
        read(1,*)fhss,fhgg,fhqg,fhgq,fhqq
        close(1)
        if(.not.negjdis)goto 4
c       else        !not symmetric Qs2
c        goto 4
c       endif
      elseif(.not.producetables)then
        write(ifmt,*) "Missing epos rj file !"        
        write(ifmt,*) "Please correct the defined path ",
     &"or force production ..."
        stop
      endif

      ix=nfnii-4
      iv=nzzvex
      write(ifmt,'(2a/2a,2f7.2,a,f6.3,a)')
     &fnii(1:ix)//'rj'//q2mnch(nq2mn)//q2mnch(nq2xx)//zzvexch(iv)//'.i'
     .,' does not exist ','     -> calc om5 tabs '
     .,'for q2jmin,q2lmin =',q2jmin,q2lmin,' ...'
      call clop(3)
      if(nq2mnfix.eq.0)then
        write(ifmt,*)'  '
        write(ifmt,*)'  ==> Better use nq2mnfix -3 '
        write(ifmt,*)'  ==> to create rj tables in parallel'
        write(ifmt,*)'  ==> by running',maxq2mn*maxq2mn,'  batch jobs'
        write(ifmt,*)'  '
        write(ifmt,*)'  ==> If not,  use nq2mnfix 99 '
        write(ifmt,*)'  '
        stop
      endif

      iclpros=iclpro
      icltars=icltar
      spmin=4.*q2zmin
      spminc=4.*q2zmin+2.*qcmass**2 !bg charm
      spminb=4.*q2zmin+2.*qbmass**2 !bg bottom
      icltar=2
      kinirj=1
      q2cmin(1)=q2jmin
      q2cmin(2)=q2lmin
      call ipoBasics(q2cmin)
      call ipoCSzeroTables(q2cmin) 

      do iclpro=iclpro1,iclpro2   !hadron type  1 - pion, 2 - nucleon
      do icltar=icltar1,icltar2   !hadron type  3 - kaon, 4 - charm
      iclpt=iclpro+4*(icltar-1)
      
        do iy=1,11 ! --------- s loop -----------
        sy=spmin*(epmax/2./spmin)**((iy-1)/10.)
        syc=spminc*(epmax/2./spminc)**((iy-1)/10.) !bg charm
        syb=spminb*(epmax/2./spminb)**((iy-1)/10.) !bg bottom
        do ix1=1,10  ! --------- xp loop -----------
c        if(ix1.le.6)then
          xpp=spmin/sy*(sy/spmin)**((ix1-1)/9.)
c        else
c          xpp=.25*(ix1-6)
c        endif
        do ix2=1,10  ! --------- xm loop -----------
c        if(ix2.le.6)then
          xmm=spmin/sy*(sy/spmin)**((ix2-1)/9.)  !.1*2.**(ix2-5)
c        else
c          xmm=.25*(ix2-6)
c        endif
        ix12=ix1+10*(ix2-1)
          
              !~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
              !if(iclpro.eq.2.and.icltar.eq.2
              !. .and.mod(ix1,3).eq.0.and.mod(ix2,3).eq.0)
              !.write(ifmt,*)'check ipoCSzero ',sqrt(sy),q2cmin,'  '
              !.,ppl/scom,pmi/scom,'  ',cstotzero(15,1,1),
              !. om11pp(sy,1.,0.5,0)
              !~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
          do iz=1,10
            z=.1*iz
            syx=sy*xpp*xmm
            ftmp=0.
            if(isetcs.ge.1)ftmp=om11pp4(syx,z,10)
            if(ftmp.le.0.)then
              fhss(iy,ix12,iz,iclpt)=-80.
            else
              fhss(iy,ix12,iz,iclpt)=log(ftmp/z)
            endif
           if(.not.negjdis)then
            ftmp=om11pp(syx,1.d0,z,0)
            if(ftmp.le.0.)then
              fhgg(iy,ix12,iz,iclpt)=-80.
            else
              fhgg(iy,ix12,iz,iclpt)=log(ftmp/z)
            endif
             !~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
             !if(ix12.eq.45.and.iz.eq.5)then
             ! write(ifmt,*)'  check fhgg  ',
             !.  iy,ix12,iz,iclpt,log(om11pp(sy,1.,z,0)/z)
             !stop'check fhgg  '
             !endif
             !~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

            if(iclpro.lt.4)then !bg
              syx=sy*xpp*xmm
            elseif(iclpro.eq.4) then !bg charm
              syx=syc*xpp*xmm
            else                !bg bottom
              syx=syb*xpp*xmm
            endif
            ftmp=om11pp(syx,dble(xpp),z,1)
            if(ftmp.le.0.)then
              fhqg(iy,ix12,iz,iclpt)=-80.
            else
              fhqg(iy,ix12,iz,iclpt)=log(ftmp/z)
            endif
            if(icltar.lt.4)then !bg
              syx=sy*xpp*xmm
            elseif(icltar.eq.4)then !bg charm
              syx=syc*xpp*xmm
            else                !bg bottom
              syx=syb*xpp*xmm
            endif
            ftmp=om11pp(syx,dble(xmm),z,2)
            if(ftmp.le.0.)then
              fhgq(iy,ix12,iz,iclpt)=-80.
            else
              fhgq(iy,ix12,iz,iclpt)=log(ftmp/z)
            endif
           endif
          enddo
          
         if(.not.negjdis)then
          if(iy.eq.1)then
            fhqq(iy,ix12,iclpt)=-80.
          else
            if(iclpro.lt.4.and.icltar.lt.4)then !bg
              syx=sy*xpp*xmm
            elseif(iclpro.eq.4.or.icltar.eq.4)then !bg charm
              syx=syc*xpp*xmm
            else                                  !bg bottom
              syx=syb*xpp*xmm
            endif
            fhqq(iy,ix12,iclpt)=
     *      log(om11pp3(syx,dble(xpp),dble(xmm)))
          endif
         endif
          
        enddo
        enddo
        enddo

      enddo
      enddo

      kinirj=0

      ix=nfnii-4
      iv=nzzvex
      open(1,file=
     &fnii(1:ix)//'rj'//q2mnch(nq2mn)//q2mnch(nq2xx)//zzvexch(iv)//'.i'
     . ,status='unknown')
      write (1,*)q2sft,alpparh,alplea,alppom,slopom
     *,gamhad,r2had,facpdf4sat,facpdf5sat
     *,qcdlam,q2jmin,q2lmin,q2ini,betpomi,glusea,naflav,nbflav,factk
     *,gampar,gamtil,zoeinc,zzsoft
      write (1,*)fhss,fhgg,fhqg,fhgq,fhqq
      call checkpsar4
      close(1)
      write(ifmt,'(2a)')'     written to ',
     &fnii(1:ix)//'rj'//q2mnch(nq2mn)//q2mnch(nq2xx)//zzvexch(iv)//'.i'
      iclpro=iclpros
      icltar=icltars

4     continue

      do i=1,11
      do i1=1,10  
      do i2=1,10  
        i12=i1+10*(i2-1)
        do m=mnclpt,mxclpt
         mx=m-mnclpt+1
         do k=1,10
         fhxss(nq2mn,nq2xx,nzzvex,i,i12,k,mx)=fhss(i,i12,k,m)
         fhxgg(nq2mn,nq2xx,nzzvex,i,i12,k,mx)=fhgg(i,i12,k,m)
         fhxqg(nq2mn,nq2xx,nzzvex,i,i12,k,mx)=fhqg(i,i12,k,m)
         fhxgq(nq2mn,nq2xx,nzzvex,i,i12,k,mx)=fhgq(i,i12,k,m)
         !call omtabset(1,nq2mn,nq2xx,nzzvex,i,i12,mx,k, fhgg(i,i12,k,m))
         !call omtabset(2,nq2mn,nq2xx,nzzvex,i,i12,mx,k, fhqg(i,i12,k,m))
         !call omtabset(3,nq2mn,nq2xx,nzzvex,i,i12,mx,k, fhgq(i,i12,k,m))
         enddo    
        fhxqq(nq2mn,nq2xx,nzzvex,i,i12,mx)=fhqq(i,i12,m)
        !call omtabset(4,nq2mn,nq2xx,nzzvex,i,i12,mx,1, fhqq(i,i12,m))
        enddo
      enddo
      enddo   
      enddo   
     
c 111    continue

      enddo ! ======== q2jmin loop =========
      enddo ! ======== q2lmin loop =========
      
      enddo ! ======== zzverx loop =========

      write(ifmt,'(a)')' '

      if(nq2mnfix.lt.0)then
        write(ifmt,'(//5x,a//)')'Table making finished => stop'
        stop'\n\n Table making finished\n\n'
      endif  

      !call testOmExIpo(10,10)
      !call testOmExIpo(5, 10)
      !call testOmExIpo(10,5)
      !stop
      !call testOmExIpo(7,7)

      q2cmin(1)=-99999.
      q2cmin(2)=-99999.
      call ipoBasics(q2cmin)
      call ipoCSzeroTables(q2cmin)
      call ipoOm5Tables(-99999)

      endif       !isetcs>0
      
 9999 continue 
    
      call utprix('mkCsOm',ish,ishini,4)
      return
      end

c-----------------------------------------------------------------------
      subroutine testOmEx(ii,jj,xp,xm,kk)
c-----------------------------------------------------------------------
#include "aaa.h"
#include "tab.h"
#include "sem.h"
      common/cpriom/npriom

      npriom=0
      
      if(npriom.eq.0)return
      
      q2jmin=q2mnval(ii)
      q2lmin=q2mnval(jj)
      
      !test pour s_pp=engy*engy

      sy=xp*xm*engy*engy
      write(ifmt,*)'---------------------------------------'
      q2cmin(1)=q2jmin 
      q2cmin(2)=q2lmin 
      call ipoBasics(q2cmin)
      call ipoCSzeroTables(q2cmin) 
      z=1.
      a = om11pp(sy,1.d0,z,0) 
      b =  ( factk*sy**delh                              
     *      *alpff(1)*alpff(2)*xp**(betff(1))
     *      *xm**(betff(2)) )
      write(ifmt,*)'testOmEx ',q2cmin,a,a*b
      write(ifmt,*)'---------------------------------------'
      if(kk.eq.1)stop
      end

c-----------------------------------------------------------------------
      subroutine testOmExIpo(ii,jj)
      !
      ! allows to compare om5 interpolation and direct calculation
      !
      !  *** agreement not so great for VERY SMALL x  *****
      !
c-----------------------------------------------------------------------
#include "aaa.h"
#include "sem.h"
#include "tab.h"
      double precision om51p
      real zz(26,4),yy(4)
      common/cpriom/npriom
      common /psar4/  fhgg(11,100,10,8),fhqg(11,100,10,8)
      real valu(2) !        xh     ,  yp
      data valu    /  0.8  , 0.0   /
      !data valu    /  6.7945e-7  , 0.0   /
 
c      return
      !npriom=2
      mij=26
       
      q2jmin=q2mnval(ii)
      q2lmin=q2mnval(jj)
       
      print*,' '
      do ij=mij,26
      i=1+(ij-1)/5
      j=1+mod(ij-1,5)
      xp=-0.1+0.2*j
      xm=-0.1+0.2*i
      xh=xp*xm
      yp=0.5*log(xp/xm)
      if(ij.eq.26)then
        xh=valu(1)
        sy=engy*engy*xh
        yp=valu(2)
        xp=sqrt(xh)*exp(yp)
        xm=sqrt(xh)*exp(-yp)
      endif
      q2cmin(1)=q2jmin 
      q2cmin(2)=q2lmin 
      call ipoBasics(q2cmin)  
      call ipoCSzeroTables(q2cmin) 
      b=0.
      sy=xh*engy*engy
      rp=r2had(iclpro)+r2had(icltar)+slopom*log(max(1.,sy))
      zb=exp(-b**2/(4.*.0389*rp))
      zz(ij,1)=om11pp(sy,1.d0,zb,0) *  
     *   ( factk*sy**delh                              
     *     *alpff(1)*alpff(2)*xp**betff(1)
     *     *xm**betff(2) )
      zz(ij,2)=om11pp(sy,dble(xp),zb,1) *
     *      (  factk*sy**delh                             
     *         *alpff(2)*xm**betff(2)   )
      zz(ij,3)=om11pp(sy,dble(xm),zb,2) *
     *      (  factk*sy**delh                             
     *         *alpff(1)*xp**betff(1)   )
      zz(ij,4)=om11pp3(sy,dble(xp),dble(xm)) *
     *        ( zb**(rp/(r2had(iclpro)+r2had(icltar)))
     *          *factk*sy**delh
     *          /(4.*pi*(r2had(iclpro)+r2had(icltar)))  )
      enddo

      q2pmin(1)=q2jmin
      q2pmin(2)=q2lmin
      call ipoOm5Tables(1)  
 
      print*,' '
      do ij=mij,26
      i=1+(ij-1)/5
      j=1+mod(ij-1,5)
      xp=-0.1+0.2*j
      xm=-0.1+0.2*i
      xh=xp*xm
      yp=0.5*log(xp/xm)
      sy=xh*engy*engy
      if(ij.eq.26)then
        xh=valu(1)
        sy=engy*engy*xh
        yp=valu(2)
      endif
      b=0.
      rp=r2had(iclpro)+r2had(icltar)+slopom*log(max(1.,sy))
      zb=exp(-b**2/(4.*.0389*rp))
      do iqq=1,4
ctp      yy(iqq)=2*om51p(sy,dble(xh),dble(yp),b,iqq)
      yy(iqq)=om51p(sy,dble(xh),dble(yp),b,iqq)
      enddo
      write(ifmt,'(a,4(2f8.3,1x))')
     . 'testOmExIpo',(yy(iqq),zz(ij,iqq),iqq=1,4),q2cmin
      !stop
      enddo
      !stop
      end



c####################################################################################################
c####################################################################################################
c####################################################################################################
c####################################################################################################
c####
c####                      interpolation functions
c####
c####        function psboi ...... born 
c####        function psjti1 ..... ordered (sigma_ord)
c####        function psjti ...... total (sigma_hard)
c####        subroutine psjti0 ... total and born for minimal virtuality cutoffs
c####
c####                     kinematics utility programs
c####
c#### We rigorouly differentiate between virtuality (scale), transverse momentum squared (pt2) 
c#### and Mandelstam variable (t). They are uniquely related for a given klas, and their relations
c#### are obtained via the following function (and nothing else!!!).  The "2" in the names refers to 
c#### double precision arguments, the functions are as well defined without the "2", as single precision.
c#### 
c####  subroutine getBornPtr2(klas,s,t , pt2)  
c####                                  get pt2 for given s,t
c####
c####  subroutine getBornKin2(klas,s,pt2 , smin,t,tmax,iret )   
c####                                  get smin,t,tmax for given s,pt2
c####
c####  subroutine getBornKin2(klas,s,pt2min , smin,tmin,tmax,iret)  
c####                                  get smin,tmin,tmax for given s,pt2min
c####
c####  subroutine HardScale2(klas, scale, pt2, ii ) 
c####                                  defines factorization scale
c####                                  allows to get pt2 from scale (ii=1) or get scale from pt2  (ii=2) 
c####        
c#### and to make the link with the EoL parton pairs jl (to be used ONLY for limits! jl is NOT the Born entry)
c####     
c####  integer function igetKlasMaxEoL(jl)       get the number of klas values (number of classes) for given jl
c####  integer function igetKlasEoL(jl,k)        get kth klas for given jl  (jl = End of Ladder parton pair)
c####  subroutine getSminEoL(jl,scale , smin)    get smin for given jl and scale
c####  subroutine getQmaxEoL(jl,sk , qmax)       get qmax (virtuality) for given jl and s
c####  subroutine getLimitsEoL(jl,sk,qmin , smin,tmin,tmax,iret)  
c####                                            get smin,tmin,tmax,iret for given jl,s, and qmin (min scale)
c####  
c####   ---->   klas defined in subroutine getKlasString
c####  
c####################################################################################################
c####################################################################################################
c####################################################################################################
c####################################################################################################

c--------------------------------------------------------------------------------------
      subroutine getBornU(klas,s,t , u)
c--------------------------------------------------------------------------------------
      double precision u2
      call getBornU2(klas,dble(s),dble(t) , u2)
      u=sngl(u2)
      return
      end
c--------------------------------------------------------------------------------------
      subroutine getBornU2(klas,s,t , u)
c--------------------------------------------------------------------------------------
c computes u
c--------------------------------------------------------------------------------------
#include "aaa.h"
#include "sem.h"
      integer ienvi
      common /cienvi/ ienvi
      double precision m2, s, t, u
      if(klas.lt.1.or.klas.gt.klasmax)then
        write(ifmt,'(a,$)')'ERROR in getBornU2 WRONG choice;  '
        write(ifmt,*)' klas = ',klas,'  ienvi =',ienvi
        stop
      endif
      go to (1,2,3,4,5,6,7,8,9,10,11,12) klas
  1   continue ! (0  0  0  0 )
      m2=0
      goto 999
  2   continue ! (mc 0  mc 0 )
      m2=2*qcmass**2
      goto 999
  3   continue ! (mb 0  mb 0 )
      m2=2*qbmass**2
      goto 999
  4   continue ! (0  0  mc mc)
      m2=2*qcmass**2
      goto 999
  5   continue ! (0  0  mb mb)
      m2=2*qbmass**2
      goto 999
  6   continue ! (mc mc 0  0 )
      m2=2*qcmass**2
      goto 999
  7   continue ! (mb mb 0  0 )
      m2=2*qbmass**2
      goto 999
  8   continue ! (mc mc mc mc)
      m2=4*qcmass**2
      goto 999
  9   continue ! (mc mc mb mb)
      m2=2*qcmass**2+2*qbmass**2
      goto 999
 10   continue ! (mb mb mc mc)
      m2=2*qbmass**2+2*qcmass**2
      goto 999
 11   continue ! (mb mb mb mb)
      m2=4*qbmass**2
      goto 999
 12   continue ! (mc mb mc mb)
      m2=2*qcmass**2+2*qbmass**2
 999  continue
      u=s-t-m2
      return
      end


c--------------------------------------------------------------------------------------
      subroutine getBornPtr(klas,s,t , pt2)
c--------------------------------------------------------------------------------------
      double precision pt2d
      call getBornPtr2(klas,dble(s),dble(t) , pt2d)
      pt2=sngl(pt2d)
      return
      end
c--------------------------------------------------------------------------------------
      subroutine getBornPtr2(klas,s,t , pt2)
c--------------------------------------------------------------------------------------
c computes transverse momentum squared pt2
c--------------------------------------------------------------------------------------
#include "aaa.h"
#include "sem.h"
      integer ienvi
      common /cienvi/ ienvi
      double precision m2, m2p, W, s, t, pt2
      if(klas.lt.1.or.klas.gt.klasmax)then
        write(ifmt,'(a,$)')'ERROR in getBornPtr2 WRONG choice;  '
        write(ifmt,*)' klas = ',klas,'  ienvi =',ienvi
        stop
      endif
      go to (1,2,3,4,5,6,7,8,9,10,11,12) klas
  1   continue ! (0  0  0  0 )
      W=s
      pt2=t*(1-t/W)
      goto 999
  2   continue ! (mc 0  mc 0 )
      m2=qcmass**2
      W=(s-m2)**2/s
      pt2=t*(1-t/W)
      goto 999
  3   continue ! (mb 0  mb 0 )
      m2=qbmass**2
      W=(s-m2)**2/s
      W=max(W,1e-10)
      pt2=t*(1-t/W)
      goto 999
  4   continue ! (0  0  mc mc)
      m2=qcmass**2
      W=s
      pt2=(t+m2)*(1-(t+m2)/W)-m2
      goto 999
  5   continue ! (0  0  mb mb)
      m2=qbmass**2
      W=s
      pt2=(t+m2)*(1-(t+m2)/W)-m2
      goto 999
  6   continue ! (mc mc 0  0 )
      m2=qcmass**2
      W=s-4*m2
      pt2=(t-m2)*(1-(t-m2)/W)+m2
      goto 999
  7   continue ! (mb mb 0  0 )
      m2=qbmass**2
      W=s-4*m2
      W=max(W,1e-10)
      pt2=(t-m2)*(1-(t-m2)/W)+m2
      goto 999
  8   continue ! (mc mc mc mc)
      m2=qcmass**2
      W=s-4*m2
      pt2=t*(1-t/W)
      goto 999
  9   continue ! (mc mc mb mb)
      m2=qcmass**2
      m2p=qbmass**2
      W=s-4*m2
      pt2=(t-m2+m2p)*(1-(t-m2+m2p)/W)+m2-m2p
      goto 999
 10   continue ! (mb mb mc mc)
      m2=qbmass**2
      m2p=qcmass**2
      W=s-4*m2
      W=max(W,1e-10)
      pt2=(t-m2+m2p)*(1-(t-m2+m2p)/W)+m2-m2p
      goto 999
 11   continue ! (mb mb mb mb)
      m2=qbmass**2
      W=s-4*m2
      W=max(W,1e-10)
      pt2=t*(1-t/W)
      goto 999
 12   continue ! (mc mb mc mb)
      m2=qcmass**2
      m2p=qbmass**2
      W=s-2*(m2+m2p)+(m2+m2p)/s
      pt2=t*(1-t/W)
 999  continue
      if(.not.(pt2.le.0..or.pt2.ge.0.))then !NaN catch
        write(ifmt,'(a,2i6,$)')'ERROR getBornPtr NaN',ienvi,klas
        write(ifmt,*)s,t,pt2,W,t-m2
        stop
      endif
      return
      end

c--------------------------------------------------------------------------------------
      subroutine getBornKin(klas, sk, ppp  ,  smin, ttt, tmax,iret)
c--------------------------------------------------------------------------------------
      double precision tmi2 , tma2, smi2
      call getBornKin2(klas,dble(sk),dble(ppp),smi2,tmi2,tma2,iret)
      smin=sngl(smi2)
      ttt=sngl(tmi2)
      tmax=sngl(tma2)
      end
c--------------------------------------------------------------------------------------
      subroutine getKlasString(klas,strg)
c--------------------------------------------------------------------------------------
#include "sem.h"
      character*12 strg,strgx(klasmax)
      !  x~ means: anti-x,  x` means: not necessarily same flavor as x  
      data strgx/
     . 'l l` -> l l`'!1 !annihilation/scattering of light partons           (0  0  0  0 )
     .,'c l  -> c l '!2 !scattering of charmed (anti)quarks with light ones (mc 0  mc 0 )
     .,'b l  -> b l '!3 !scattering of bottom (anti)quarks with light ones  (mb 0  mb 0 )
     .,'l l~ -> c c~'!4 !annihilation of light pairs into charmed pairs     (0  0  mc mc)
     .,'l l~ -> b b~'!5 !annihilation of light pairs into bottom pairs      (0  0  mb mb)
     .,'c c~ -> l l~'!6 !annihilation of charmed pairs into light pairs     (mc mc 0  0 )
     .,'b b~ -> l l~'!7 !annihilation of bottom pairs into light pairs      (mb mb 0  0 )
     .,'c c` -> c c`'!8 !annihilation/scattering of charmed into charmed    (mc mc mc mc)
     .,'c c~ -> b b~'!9 !annihilation of charmed pairs into bottom pairs    (mc mc mb mb)
     .,'b b~ -> c c~'!10!annihilation of bottom pairs into charm pairs      (mb mb mc mc)
     .,'b b` -> b b`'!11!annihilation/scattering of bottom into bottom      (mb mb mb mb)
     .,'c b  -> c b '!12!scattering of charm and bottom partons             (mc mb mc mb)
     . /
      strg=strgx(klas)
      end
c--------------------------------------------------------------------------------------
      integer function igetKlas(i,j,k,l)
c--------------------------------------------------------------------------------------
#include "aaa.h"
      !i,j = Born in   k,l = Born out
      ii=abs(i)
      jj=abs(j)
      kk=abs(k)
      ll=abs(l)
      iii=max(ii,jj)
      jjj=min(ii,jj)
      kkk=max(kk,ll)
      lll=min(kk,ll)
      igetKlas=0
      if(ii.le.3.and.jj.le.3.and.kk.le.3.and.ll.le.3)then
        igetKlas=1
      elseif(iii.eq.4.and.jjj.le.3.and.kkk.eq.4.and.lll.le.3)then
        igetKlas=2
      elseif(iii.eq.5.and.jjj.le.3.and.kkk.eq.5.and.lll.le.3)then
        igetKlas=3
      elseif(ii.le.3.and.i.eq.-j.and.kk.eq.4.and.k.eq.-l)then
        igetKlas=4
      elseif(ii.le.3.and.i.eq.-j.and.kk.eq.5.and.k.eq.-l)then
        igetKlas=5
      elseif(ii.eq.4.and.i.eq.-j.and.kk.le.3.and.k.eq.-l)then
        igetKlas=6
      elseif(ii.eq.5.and.i.eq.-j.and.kk.le.3.and.k.eq.-l)then
        igetKlas=7
      elseif(ii.eq.4.and.jj.eq.4.and.kk.eq.4.and.ll.eq.4)then
        igetKlas=8
      elseif(ii.eq.4.and.i.eq.-j.and.kk.eq.5.and.k.eq.-l)then
        igetKlas=9
      elseif(ii.eq.5.and.i.eq.-j.and.kk.eq.4.and.k.eq.-l)then
        igetKlas=10
      elseif(ii.eq.5.and.jj.eq.5.and.kk.eq.5.and.ll.eq.5)then
        igetKlas=11
      elseif(iii.eq.5.and.jjj.eq.4.and.kkk.eq.5.and.lll.eq.4)then
        igetKlas=12
      else
        write(ifmt,'(a,4i4)')'ERROR igetKlas: no klas found;',i,j,k,l
        !stop
      endif
      end
c--------------------------------------------------------------------------------------
      subroutine getBornKin2(klas, sk, ppp  ,  smin, ttt, tmax,iret)
c--------------------------------------------------------------------------------------
c get variables related to Born kinematics -- double precision
c--------------------------------------------------------------------------------------
c Input:
c   klas - Born reaction class
c    -1 consider all classes, compute min{ttt},max(tmax}, consider ppp as scale, not pt2 
c     0 consider all classes, compute min{ttt},max(tmax} 
c    >0 consider klas 
c   sk ...... Mandelstam s 
c   ppp ...... pt2  or  pt2_min   <----- two options for using this routine
c----------------------------------------------------------------------------------------
c Output:
c   smin .... s_{min}
c   ttt .... |t|  or  |t|_min    <----- two options for using this routine
c   tmax .... t_max+ 
c   iret .... 1 if no phase space, 0 means OK
c----------------------------------------------------------------------------------------
c We may use this routine also to compute |t| for given pt:
c    when called with argument ppp = pt then ttt returns |t|
c----------------------------------------------------------------------------------------
c integrals can be split:
c  \int f(t) dt   (integr domain t_{min} to t_{max})
c  = \int { f(t) + f( 2*t_{max+}-t) } dt   (integr domain t_{min} to t_{max+})
c----------------------------------------------------------------------------------------- 
      implicit none
#include "aaa.h"
#include "sem.h"
      integer ienvi
      common /cienvi/ ienvi
      double precision sk, ppp, smin, ttt, tmax, tttmin, tmaxmax 
      double precision W, m2, m2p, s, d,  sminmin, pppx, uuu, rr!, Wp
      integer klas,iret,icrash
      real crash(2)
      icrash=3
      if(klas.lt.-1.or.klas.gt.klasmax)then
        write(*,*)'-----> wrong choice: klas = ',klas
        crash(icrash)=1. !to force crash
        stop'wrong klas'
      endif

      m2=0d0
      m2p=0d0
      s=sk
      sminmin = 1d30  
      tttmin = 1d30
      tmaxmax = 0d0 
      pppx=ppp

      go to (1,2,3,4,5,6,7,8,9,10,11,12) klas

  1   continue ! (0  0  0  0 )
      if(klas.eq.-1)call HardScale2(1,pppx,ppp,1) !transform scale to pt2
      W=s
      tmax=W/2
c      Wp=W 
      smin=4*ppp
      if(klas.gt.0.and.s.lt.smin)goto 888
      rr=1-4*ppp/W
      if(klas.gt.0.and.rr.lt.0d0)goto 888
      ttt=2*ppp/(1+sqrt(rr))
      if(klas.gt.0)goto 999
      tttmin=min(tttmin,ttt)
      tmaxmax=max(tmaxmax,tmax)
      sminmin=min(sminmin,smin)

  2   continue ! (mc 0  mc 0 )
      if(klas.eq.-1)call HardScale2(2,pppx,ppp,1) !transform scale to pt2
      m2=qcmass**2
      W=(s-m2)**2/s
      tmax=W/2
c      Wp=W 
      d=m2+2*ppp
      smin=d*(1+sqrt(1-m2**2/d**2))
      if(klas.gt.0.and.s.lt.smin)goto 888
      if(W.gt.0.d0)then
        rr=1-4*ppp/W
        if(klas.gt.0.and.rr.lt.0d0)goto 888
        ttt=2*ppp/(1+sqrt(rr))
      elseif(W.le.0.d0)then
        write(ifmt,*)'Error in getBornKin2',klas,W,s,m2,ppp,smin
        call utstop("Error in getBornKin2 (2) !&")
      endif
      if(klas.gt.0)goto 999
      tttmin=min(tttmin,ttt)
      tmaxmax=max(tmaxmax,tmax)
      sminmin=min(sminmin,smin)

  3   continue ! (mb 0  mb 0 )
      if(klas.eq.-1)call HardScale2(3,pppx,ppp,1) !transform scale to pt2
      m2=qbmass**2
      W=(s-m2)**2/s
      tmax=W/2
c      Wp=W 
      d=m2+2*ppp
      smin=d*(1+sqrt(1-m2**2/d**2))
      if(klas.gt.0.and.s.lt.smin)goto 888
      if(W.gt.0.d0)then
        rr=1-4*ppp/W
        if(klas.gt.0.and.rr.lt.0d0)goto 888
        ttt=2*ppp/(1+sqrt(rr))
      elseif(W.le.0.d0)then
        write(ifmt,*)'Error in getBornKin2',klas,W,s,m2,ppp,smin
        call utstop("Error in getBornKin2 (3) !&")
      endif
      if(klas.gt.0)goto 999
      tttmin=min(tttmin,ttt)
      tmaxmax=max(tmaxmax,tmax)
      sminmin=min(sminmin,smin)

  4   continue ! (0  0  mc mc)
      if(klas.eq.-1)call HardScale2(4,pppx,ppp,1) !transform scale to pt2
      m2=qcmass**2
      W=s
      tmax=s/2-m2
c      Wp=s-4*m2
      smin=4*(m2+ppp)
      if(klas.gt.0.and.s.lt.smin)goto 888
      rr=1-4*(ppp+m2)/W
      if(klas.gt.0.and.rr.lt.0d0)goto 888
      ttt=2*(ppp+m2)/(1+sqrt(rr))-m2
      if(klas.gt.0)goto 999
      tttmin=min(tttmin,ttt)
      tmaxmax=max(tmaxmax,tmax)
      sminmin=min(sminmin,smin)

  5   continue ! (0  0  mb mb)
      if(klas.eq.-1)call HardScale2(5,pppx,ppp,1) !transform scale to pt2
      m2=qbmass**2
      W=s
      tmax=s/2-m2
c      Wp=s-4*m2
      smin=4*(m2+ppp)
      if(klas.gt.0.and.s.lt.smin)goto 888
      rr=1-4*(ppp+m2)/W
      if(klas.gt.0.and.rr.lt.0d0)goto 888
      ttt=2*(ppp+m2)/(1+sqrt(rr))-m2
      if(klas.gt.0)goto 999
      tttmin=min(tttmin,ttt)
      tmaxmax=max(tmaxmax,tmax)
      sminmin=min(sminmin,smin)

  6   continue ! (mc mc 0  0 )
      if(klas.eq.-1)call HardScale2(6,pppx,ppp,1) !transform scale to pt2
      m2=qcmass**2
      W=s-4*m2
      if(W.le.0.d0)goto 888
      tmax=s/2-m2
c      Wp=s
      smin=4*ppp
      if(klas.gt.0.and.s.lt.smin)goto 888
      if(W.gt.0.d0)then
        rr=1-4*(ppp-m2)/W
        if(klas.gt.0.and.rr.lt.0d0)goto 888
        ttt=2*(ppp-m2)/(1+sqrt(rr))+m2
      elseif(W.eq.0.d0)then
        if(m2-ppp.lt.0.)then
          write(ifmt,*)'Warning getBornKin2',klas,ppp-m2
     .         ,ppp,m2
          if(ppp-m2.gt.0.0001)call utstop("Error in getBornKin2 (6) !&")
        endif
        ttt=m2
      else
        goto 888
      endif
      if(klas.gt.0)goto 999
      tttmin=min(tttmin,ttt)
      tmaxmax=max(tmaxmax,tmax)
      sminmin=min(sminmin,smin)

  7   continue ! (mb mb 0  0 )
      if(klas.eq.-1)call HardScale2(7,pppx,ppp,1) !transform scale to pt2
      m2=qbmass**2
      W=s-4*m2
      if(W.le.0.d0)goto 888
      tmax=s/2-m2
c      Wp=s
      smin=4*ppp
      if(klas.gt.0.and.s.lt.smin)goto 888
      if(W.gt.0.d0)then
        rr=1-4*(ppp-m2)/W
        if(klas.gt.0.and.rr.lt.0d0)goto 888
        ttt=2*(ppp-m2)/(1+sqrt(rr))+m2
      elseif(W.eq.0.d0)then
        if(m2-ppp.lt.0.)then
          write(ifmt,*)'Warning getBornKin2',klas,ppp-m2
     .         ,ppp,m2
          if(ppp-m2.gt.0.0001)call utstop("Error in getBornKin2 (7) !&")
        endif
        ttt=m2
      else
        goto 888
      endif
      if(klas.gt.0)goto 999
      tttmin=min(tttmin,ttt)
      tmaxmax=max(tmaxmax,tmax)
      sminmin=min(sminmin,smin)

  8   continue ! (mc mc mc mc)
      if(klas.eq.-1)call HardScale2(8,pppx,ppp,1) !transform scale to pt2
      m2=qcmass**2
      m2p=m2
      W=s-4*m2
      tmax=W/2
c      Wp=W
      smin=4*(m2+ppp)
      if(klas.gt.0.and.s.lt.smin)goto 888
      if(W.gt.0.d0)then
        rr=1-4*ppp/W
        if(klas.gt.0.and.rr.lt.0d0)goto 888
        ttt=2*ppp/(1+sqrt(rr))
      elseif(W.le.0.d0)then
        write(ifmt,*)'Error in getBornKin2',klas,W,s,m2,ppp,smin
        call utstop("Error in getBornKin2 (8) !&")
      endif
      if(klas.gt.0)goto 999
      tttmin=min(tttmin,ttt)
      tmaxmax=max(tmaxmax,tmax)
      sminmin=min(sminmin,smin)

  9   continue ! (mc mc mb mb)
      if(klas.eq.-1)call HardScale2(9,pppx,ppp,1) !transform scale to pt2
      m2=qcmass**2
      m2p=qbmass**2
      W=s-4*m2
      if(W.le.0.d0)goto 888
      tmax=s/2-m2-m2p
c      Wp=s-4*m2p
      smin=4*(m2p+ppp)
      if(klas.gt.0.and.s.lt.smin)goto 888
      if(W.gt.0.d0)then
        rr=1-4*(ppp-m2+m2p)/W
        if(klas.gt.0.and.rr.lt.0d0)goto 888
        ttt=2*(ppp-m2+m2p)/(1+sqrt(rr))+m2-m2p
      elseif(W.eq.0.d0)then
        if(m2-ppp-m2p.lt.0.)then
          write(ifmt,*)'Warning getBornKin2',klas,ppp-m2+m2p
     .         ,ppp,m2,m2p
          if(ppp-m2+m2p.gt.0.0001)
     .         call utstop("Error in getBornKin2 (9) !&")
        endif
        ttt=m2-m2p
      else
        goto 888
      endif
      if(klas.gt.0)goto 999
      tttmin=min(tttmin,ttt)
      tmaxmax=max(tmaxmax,tmax)
      sminmin=min(sminmin,smin)

 10   continue ! (mb mb mc mc)
      if(klas.eq.-1)call HardScale2(10,pppx,ppp,1) !transform scale to pt2
      m2=qbmass**2
      m2p=qcmass**2
      W=s-4*m2
      if(W.le.0.d0)goto 888
      tmax=s/2-m2-m2p
c      Wp=s-4*m2p
      smin=4*(m2p+ppp)
      if(klas.gt.0.and.s.lt.smin)goto 888
      if(W.gt.0.d0)then
        rr=1-4*(ppp-m2+m2p)/W
        if(klas.gt.0.and.rr.lt.0d0)goto 888
        ttt=2*(ppp-m2+m2p)/(1+sqrt(rr))+m2-m2p
      elseif(W.eq.0.d0)then
        if(m2-ppp-m2p.lt.0.)then
          write(ifmt,*)'Warning getBornKin2',klas,ppp-m2+m2p
     .         ,ppp,m2,m2p
          if(ppp-m2+m2p.gt.0.0001)
     .         call utstop("Error in getBornKin2 (10) !&")
        endif
        ttt=m2-m2p
      else
        goto 888
      endif
      if(klas.gt.0)goto 999
      tttmin=min(tttmin,ttt)
      tmaxmax=max(tmaxmax,tmax)
      sminmin=min(sminmin,smin)

 11   continue ! (mb mb mb mb)
      if(klas.eq.-1)call HardScale2(11,pppx,ppp,1) !transform scale to pt2
      m2=qbmass**2
      m2p=m2
      W=s-4*m2
      tmax=W/2
c      Wp=W
      smin=4*(m2+ppp)
      if(klas.gt.0.and.s.lt.smin)goto 888
      if(W.gt.0.d0)then
        rr=1-4*ppp/W
        if(klas.gt.0.and.rr.lt.0d0)goto 888
        ttt=2*ppp/(1+sqrt(rr))
      elseif(W.le.0.d0)then
        write(ifmt,*)'Error in getBornKin2',klas,W,s,m2,ppp,smin
        call utstop("Error in getBornKin2 (11) !&")
      endif
      if(klas.gt.0)goto 999
      tttmin=min(tttmin,ttt)
      tmaxmax=max(tmaxmax,tmax)
      sminmin=min(sminmin,smin)

 12   continue ! (mc mb mc mb)
      if(klas.eq.-1)call HardScale2(12,pppx,ppp,1) !transform scale to pt2
      m2=qcmass**2
      m2p=qbmass**2
      W=s-2*(m2+m2p)+(m2-m2p)**2/s
      tmax=W/2
c      Wp=W 
      d=m2+m2p+2*ppp
      smin=d*(1+sqrt(1-((m2-m2p)/d)**2))
      if(klas.gt.0.and.s.lt.smin)goto 888
      if(W.gt.0.d0)then
        rr=1-4*ppp/W
        if(klas.gt.0.and.rr.lt.0d0)goto 888
        ttt=2*ppp/(1+sqrt(rr))
      elseif(W.le.0.d0)then
        write(ifmt,*)'Error in getBornKin2',klas,W,s,m2,ppp,smin
        call utstop("Error in getBornKin2 (12) !&")
      endif
      if(klas.gt.0)goto 999
      tttmin=min(tttmin,ttt)
      tmaxmax=max(tmaxmax,tmax)
      sminmin=min(sminmin,smin)

      ttt=tttmin
      tmax=tmaxmax
      smin=sminmin
      if(s.lt.smin)goto 888

 999  continue !OK so far
      !check conditions for valid kinematics
      if(tmax.le.ttt)goto 888
      uuu=sk-ttt-2d0*m2-2*m2p !Mandelstam u
      if(uuu.le.0)goto 888    !must be positive
      uuu=sk-(2d0*tmax-ttt)-2d0*m2-2*m2p !Mandelstam u for t=2*tmax-t
      if(uuu.le.0)goto 888               !must be positive
      iret=0
      goto 1111
      
 888  continue !no phase space
      iret=1
      goto 1111
 
 1111 continue
      ppp=pppx
      if(.not.(ttt.le.0..or.ttt.ge.0.))then !NaN catch
            write(ifmt,'(a,i5,$)')'ERROR NaN getBornKin2',ienvi
            print*,ttt,klas, sk, ppp  ,  smin, ttt, tmax,iret
c            crash(icrash)=1.0
            call utstop("NaN in getBornKin2&")
      endif
      return 
      end

c-----------------------------------------------------------------------
      subroutine HardScale(klas,scale,pt2,ii)
c-----------------------------------------------------------------------
      double precision scale2,pt22
      scale2=dble(scale)
      pt22=dble(pt2)
      call HardScale2(klas,scale2,pt22,ii)
      if(ii.eq.1)then
        pt2=sngl(pt22)
      elseif(ii.eq.2)then
        scale=sngl(scale2)
      endif
      end

c-----------------------------------------------------------------------
      subroutine HardScale2(klas,scale,pt2,ii)
c-----------------------------------------------------------------------
c     defines factorization scale, relating pt2 and scale
c       scale -> pt2 (ii=1) or inverse (ii=2)
c     relates as well scale_min and pt2_min
c     scale means mu_F^2 = factorization scale
c-----------------------------------------------------------------------
#include "aaa.h"
#include "sem.h"
      double precision scale,pt2,coelaf
      real m2 
      real crash(2)
      icrash=3
      if(klas.lt.1.or.klas.gt.klasmax)then
        write(*,*)'-----> wrong choice: klas = ',klas
        crash(icrash)=1. !to force crash
        stop'wrong klas'
      endif
      go to (1,2,3,4,5,6,7,8,9,10,11,12) klas
  1   continue ! (0  0  0  0 )
      m2=0
      goto 999
  2   continue ! (mc 0  mc 0 )
      m2=qcmass**2
      goto 999
  3   continue ! (mb 0  mb 0 )
      m2=qbmass**2
      goto 999
  4   continue ! (0  0  mc mc)
      m2=qcmass**2
      goto 999
  5   continue ! (0  0  mb mb)
      m2=qbmass**2
      goto 999
  6   continue ! (mc mc 0  0 )
      m2=qcmass**2
      goto 999
  7   continue ! (mb mb 0  0 )
      m2=qbmass**2
      goto 999
  8   continue ! (mc mc mc mc)
      m2=qcmass**2
      goto 999
  9   continue ! (mc mc mb mb)
      m2=qbmass**2
      goto 999
 10   continue ! (mb mb mc mc)
      m2=qbmass**2
      goto 999
 11   continue ! (mb mb mb mb)
      m2=qbmass**2
      goto 999
 12   continue ! (mc mb mc mb)
      m2=qbmass**2
 999  continue
      coelaf=0.d0
      if(ii.eq.1)then
        pt2=coekaf*scale-coelaf*dble(m2)
        pt2=max(0.d0,pt2)
      elseif(ii.eq.2)then
        scale=(pt2+coelaf*dble(m2))/coekaf
      else
        stop'ERROR 07032019b'
      endif
      if(.not.(scale.le.0.d0.or.scale.ge.0.d0))then !NaN catch
            print*,scale,pt2
            crash(icrash)=1. !to force crash
            stop'NaN in HardScale'
      endif
      end

c--------------------------------------------------------------------------------------
      integer function igetKlasMaxEoL(jl)
c--------------------------------------------------------------------------------------
c get the number of klas values (number of Born kinematics classes)
c corresponding to a given value of jl (end of ladder (EoL) parton pair)
c The EoL partons are related to the Born kinematics, since in the extreme case of no
c emissions, the EoL partons are as well the in-Born partons 
c--------------------------------------------------------------------------------------
#include "tab.h"
      if(jl.le.0.or.jl.gt.mxpsar3)stop'jl out of range igetKlasMaxEoL'
      kx=0
      do k=1,3
        !the array jlklas(jl,k) contains the "klas" values (three at most)  
        !which contribute for a given EoL parton pair jl
        klas=jlklas(jl,k)
        if(klas.ne.0)kx=k
      enddo
      igetKlasMaxEoL=kx
      return
      end       

c--------------------------------------------------------------------------------------
      integer function igetKlasEoL(jl,k)
c--------------------------------------------------------------------------------------
c get kth klas value for given jl
c--------------------------------------------------------------------------------------
#include "tab.h"
      if(jl.le.0.or.jl.gt.mxpsar3)stop'jl out of range igetKlasEoL'
      klas=jlklas(jl,k)
      igetKlasEoL=klas
      return
      end       


c--------------------------------------------------------------------------------------
      subroutine getSminEoL(jl,scale,smin)
c--------------------------------------------------------------------------------------
c get smin for EoL parton pair jl and given scale
c--------------------------------------------------------------------------------------
#include "tab.h"

      if(jl.le.0.or.jl.gt.mxpsar3)stop'jl out of range getSminEoL'
      sminmin=1e30
      do k=1,igetKlasMaxEoL(jl)
        klas=jlklas(jl,k)
        call HardScale(klas, scale, pt2min, 1 )  
        call getBornKin(klas,1e25,pt2min , smin,ttt,tmax,iret) !get smin
        sminmin=min(sminmin,smin)
      enddo
      smin=sminmin
      end

c--------------------------------------------------------------------------------------
      subroutine getQmaxEoL(jl,sk,qmax)
c--------------------------------------------------------------------------------------
c get max virtuality qmax for EoL parton pair jl and given sk (Mandelstam s)
c--------------------------------------------------------------------------------------
#include "tab.h"
      if(jl.le.0.or.jl.gt.mxpsar3)stop'jl out of range getQmaxEoL'
      qmaxmax = 0 
      do k=1,igetKlasMaxEoL(jl)
        klas=jlklas(jl,k)
        call getBornKin(klas,sk,0.,  sdu, tdu, tmax,iret) !get tmax
        call getBornPtr(klas,sk,tmax , pt2max)         !get pt2max
        call HardScale(klas, qmax, pt2max, 2 ) !get  qmax
        qmaxmax=max(qmaxmax,qmax)
      enddo
      qmax=qmaxmax
      end

c--------------------------------------------------------------------------------------
      subroutine getQminEoL(jl,qref1,qref2,qqref)
c--------------------------------------------------------------------------------------
      implicit none
#include "aaa.h"
#include "sem.h"
      integer jl , j1 , l1, j2 , l2 
      real qref1,qref2,qqref
      call eolpi( jl , j1 , l1 )                                    
      call eolpic( j1 , l1 , j2 , l2 )
      qref1=q2zmin
      qref2=q2zmin
      qref1=max(qref1,1.5)
      qref2=max(qref2,1.5)
      if(iabs(j2).eq.4)qref1=max(qref1,qcmass**2)
      if(iabs(j2).eq.5)qref1=max(qref1,qbmass**2)
      if(iabs(l2).eq.4)qref2=max(qref2,qcmass**2)
      if(iabs(l2).eq.5)qref2=max(qref2,qbmass**2)
      qqref=max(qref1,qref2)
      return
      end

c--------------------------------------------------------------------------------------
      subroutine getLimitsEoL2(jl,sk,qmin , smin,tmin,tmax,iret)
c--------------------------------------------------------------------------------------
c get smin,tmin,tmax for EoL parton pair jl and given sk and qmin (minimum scale)
c--------------------------------------------------------------------------------------
#include "tab.h"
      double precision sk,qmin,smin,tmin,tmax
      double precision tminmin,tmaxmax,sminmin,pt2min
      if(jl.le.0.or.jl.gt.mxpsar3)stop'jl out of range getLimitsEoL2'
      iret=1
      tminmin=1d30
      tmaxmax = 0d0 
      sminmin=1d30
      do k=1,igetKlasMaxEoL(jl)
        klas=jlklas(jl,k)
        call HardScale2(klas,qmin,pt2min, 1 )  !get  pt2min
        call getBornKin2(klas,sk,pt2min,  smin, tmin, tmax,iretk) !get smin,tmin,tmax
        tmaxmax=max(tmaxmax,tmax)
        tminmin=min(tminmin,tmin)
        sminmin=min(sminmin,smin)
        iret=min(iret,iretk)
      enddo
      tmin=tminmin
      tmax=tmaxmax
      smin=sminmin    
      if(sk.lt.smin)iret=1
      if(tmax.lt.tmin)iret=1
      end
      !---single----
      subroutine getLimitsEoL(jl,sk1,qmin1 , smin1,tmin1,tmax1,iret)
#include "tab.h"
      real sk1,qmin1,smin1,tmin1,tmax1
      double precision smin,tmin,tmax
      if(jl.le.0.or.jl.gt.mxpsar3)stop'jl out of range getLimitsEoL'
      call getLimitsEoL2(jl,dble(sk1),dble(qmin1) , smin,tmin,tmax,iret)
      smin1=sngl(smin)
      tmin1=sngl(tmin)
      tmax1=sngl(tmax)
      return
      end

c-----------------------------------------------------------------------
      function psboi(q1,q2,qqcut,ss,m2,l2,jdis)
c-----------------------------------------------------------------------
c
c    hard 2->2 parton scattering born cross-section interpolation
c    based on tabulation of psborn for q1=q2=qqcut
c      (which is  integrated over t from tmin to tmax   
c         including Sudakov on both sides)
c
c-----------------------------------------------------------------------
c psboi - born cross-section interpolation
c q1 - virtuality cutoff at current end of the ladder;
c q2 - virtuality cutoff at opposite end of the ladder;
c qqcut - additional virtuality cutoff 
c s  - total c.m. energy squared for the scattering,
c m2 - parton type at current end of the ladder (0 - g, 1,-1,2,... - q)
c l2 - parton type at opposite end of the ladder (0 - g, 1,-1,2,... - q) !bg l2 can't be heavy quark
c-----------------------------------------------------------------------
      dimension wi(3),wk(3) !wj(3)
      common /psar2/  edmax,epmax
#include "aaa.h"
#include "sem.h"
#include "tab.h"
      double precision psuds

      !~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~!TEST
      !write(ifmt,'(80a1/a,$)')('-',m=1,80),'---psboi---'      !TEST 
      !call eolpib( m2 , l2 , ml )                             !TEST
      !call eolpi( ml , j1 , l1 )                              !TEST 
      !write(ifmt,'(4f8.2,6i4)')q1,q2,qqcut,ss,m2,l2,ml,j1,l1  !TEST 
      !~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~!TEST

      if(jdis.lt.0)then !.and.isetcs.le.-1)then
        call eolpib( m2 , l2 , ml ) !KW1808b in psboi
        call eolpi( ml , j1 , l1 )
        psboi=psborn(q1,q2,qqcut,ss,j1,l1,jdis,0)
        return
      endif
      
c        print *,'pbsoi',q1,q2,qqcut,ss,m2,l2,jdis
c        psb=psborn(q1,q2,qqcut,ss,j1,l1,jdis,0)
c        print *,'in',psb
      
      psboi=0.
      if(iabs(m2).eq.4.and. noflav(q1).lt.4 ) return
      if(iabs(m2).eq.5.and. noflav(q1).lt.5 ) return
      if(iabs(l2).eq.4.and. noflav(q2).lt.4 ) return
      if(iabs(l2).eq.5.and. noflav(q2).lt.5 ) return

      if(jdis.le.0)then
        qq=max(q1,q2)
      else
        qq=max(q1/4.,q2)
      endif
      !qq represents scale
      if(jdis.ge.0)then
        qq=max(qq,qqcut)
      elseif(qq.gt.qqcut)then
        goto 998
      endif
      qmin=qq 

      call eolpib( m2 , l2 , ml ) !KW1808b in psboi
      call getLimitsEoL(ml,ss,qmin , smin,tmin,tmax,iret) 
      if(iret.ne.0)goto 998    !KW1811

      call getQminEoL(ml,qref1,qref2,qqref) !used for tabulation, don't change here
      call getSminEoL(ml,qqref,spmin)
      spmin=spmin*1.10

      sl=log(ss/spmin)/log(epmax/2./spmin)*19.+1.
      k=int(sl)
      !since sk=spmin*(epmax/2./spmin)**((k-1)/19.)
      call getQmaxEoL(ml,ss,qmax)
      qmax=qmax / 1.01!1.10 
      qmax=max(qmax,1.01*qqref)

      if(jdis.ge.0)then
        qli=log(qq/qqref)/log(qmax/qqref)*19.+1.
      else
        qmax=q2mnval(maxq2mn)
        qli=log(qqcut/qqref)/log(qmax/qqref)*19.+1.
      endif

      i=int(qli)
      if(i.lt.1)i=1
      if(i.gt.18)i=18
      wi(2)=qli-i
      wi(3)=wi(2)*(wi(2)-1.)*.5
      wi(1)=1.-wi(2)+wi(3)
      wi(2)=wi(2)-2.*wi(3)
        
      if(k.lt.1)k=1
      if(k.gt.18)k=18

      if(jdis.ge.0)then

        !if(csbor(i,k+20*(ml-1),jdis+1).lt.-68.)k=k+1
        !if(k.gt.18)k=18

        !~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~!TEST
        !write(ifmt,'(a,4f10.1)')'---psboi---',ss,qq,sl,qli   !TEST 
        !write(ifmt,'(a,2i9)')'---psboi---',k,i               !TEST   
        !~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~!TEST

        vav=0
        nav=0
        do k1=1,3
          do i1=1,3
            w=csbor(i+i1-1,k+k1-1+20*(ml-1),jdis+1)
            if(w.gt.-68.0)then
              vav=vav+w
              nav=nav+1
            endif
          enddo
        enddo
        if(nav.gt.0)then
          vav=vav/nav
        endif   

        !~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~!TEST
        !do ii=0,10
        !skk=spmin*(epmax/2./spmin)**((12.5  -1)/19.)
        !call getQmaxEoL(ml,skk,qmax)
        !qmax=max(qmax / 1.10 ,1.01*qqref)
        !qii=qqref*(qmax/qqref)**(( 1+0.1*ii  -1)/19.)
        !fii=psborn(qii,qii,qii,skk,j1,l1,0,0)
        !print*,'check exact',skk,qii,log(fii)
        !enddo
        !       psborn(qii) has minimum at qii ~ 1.1, 
        !       anyhow dont need such smalll values
        !~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~!TEST

        wk(2)=sl-k
        wk(3)=wk(2)*(wk(2)-1.)*.5
        wk(1)=1.-wk(2)+wk(3)
        wk(2)=wk(2)-2.*wk(3)

        j=0
        qlj=0.
        do k1=1,3
          do i1=1,3
            w=csbor(i+i1-1,k+k1-1+20*(ml-1),jdis+1)
            if(w.lt.-58..and.vav.ne.0.)w=vav !cross section value zero
            psboi=psboi + w * wi(i1) * wk(k1)
      !~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~!TEST
      !write(ifmt,'(a,2i5,3f14.7)')'---psboi--- '             !TEST  
      !.,k+k1-1+20*(ml-1),i+i1-1,wi(i1)*wk(k1),w,psboi        !TEST
      !~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~!TEST
            a=exp(psboi) 
            if(.not.(a.gt.0..or.a.le.0.) !NaN catch
     .           .or.a.gt.1e35 .or. a.lt.-1e35 )then !Infinity catch
              b=wi(i1)
              c=99.99 
              d=wk(k1)
              e=w
              goto 998
            endif 
          enddo
        enddo
        !~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~!TEST
        !psbex=psborn(q1,q2,qqcut,ss,j1,l1,0,0)                   !TEST
        !write(ifmt,'(a,2f14.7,$)')'---psboi---',psboi,log(psbex) !TEST
        !write(ifmt,'(2e14.3)')exp(psboi),psbex                   !TEST
        !write(ifmt,'(80a1)')('-',m=1,80)                         !TEST  
        !~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~!TEST
      else
         stop'28072019' !code removed, if needed see 3275; dont understand the point of interpolation qq
      endif

      e=psboi
      if(e.ne.0.)psboi=exp(psboi)

c adding

      !tabulation for q1=q2=qq, which contains
      ! for both sides Sudavov(qq->qt) * Sudavov(qq->qt)
      !Needed : Sudavov(q1->qt) * Sudavov(q2->qt), so we need to add
      ! Sudavov(q1->qq) * Sudavov(q2->qq)  
      if(jdis.ge.0.and.qq.gt.q1)then
        psboi=psboi*sngl(psuds(qq,m2)/psuds(q1,m2))
      elseif(jdis.eq.1.and.4.*qq.gt.q1)then
        psboi=psboi*sngl(psuds(4.*qq,m2)/psuds(q1,m2))
      endif
      if(jdis.ge.0.and.qq.gt.q2)then
        psboi=psboi*sngl(psuds(qq,l2)/psuds(q2,l2))
      endif

c      print *,'out',psboi,abs(psboi-psb)

      a=psboi
      b=99.99
      c=99.99 
      d=99.99
      return
      
 998  continue
      call checkjtNAN('psboi ', psboi,  m2,l2,  a,b,c,d,e
     .,qq,qqref,99.,99.,qmax,s,  qli,99.,sl,i,99,k,  smin,tmin,tmax )
      return
      end
 
c------------------------------------------------------------------------
      function psjti1(q1,q2,qqcut,s,m2,l2,jdis)
c-----------------------------------------------------------------------
c inclusive hard cross-section interpolation, for strict order in the ladder
c q1 - virtuality cutoff at current end of the ladder
c q2 - virtuality cutoff at opposite end of the ladder
c qqcut - p_t cutoff for the born process;
c s - total c.m. energy squared for the ladder,
c m2 - parton type at current end of the ladder (0-g, 1,2,etc.-q)
c l2 - parton type at opposite end of the ladder (0-g, 1,2,etc.-q)
c-----------------------------------------------------------------------
c
c        -----> TO BE CHECKED <-----
c
c   there is still an issue concerning only m2,l2 with bottom involved
c                                            or both being charm
c
c   one passes "call getLimitsEoL" with iret=0 (no kinematics problem) 
c   but qmin is bigger than qmax
c   don't simply put a return here, but check the kinematics routines,
c   there is some inconsistency for the bottom case 
c
c    or one needs explicitely consider qq vs qmax to define iret
c
c      once fixed put into BOTH jt and jet functions and redo tables
c
c-----------------------------------------------------------------------
      dimension wi(3),wj(3),wk(3)
      common /psar2/  edmax,epmax
#include "aaa.h"
#include "sem.h"
#include "tab.h"
      double precision psuds

      call eolpib( m2 , l2 , ml )                                   
      call eolpi( ml , j1 , l1 )                                    

      !~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~!TEST
      !write(ifmt,'(80a1/a,$)')('-',m=1,80),'---psjti1---'           !TEST   
      !write(ifmt,'(3f8.2,e12.3,6i4)')q1,q2,qqcut,s,m2,l2,ml,j1,l1   !TEST 
      !~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~!TEST

      psjti1=0.
      if(iabs(m2).eq.4.and. noflav(q1).lt.4 ) return
      if(iabs(m2).eq.5.and. noflav(q1).lt.5 ) return
      if(iabs(l2).eq.4.and. noflav(q2).lt.4 ) return
      if(iabs(l2).eq.5.and. noflav(q2).lt.5 ) return

      if(jdis.eq.0)then
        qqmin=max(q1,q2)
      else
        qqmin=max(q1,q2/4.)
      endif
      qq=max(qqmin,qqcut)
      qmin=qq

      mlx=20*(ml-1)+mxpsar1*jdis !KW1808

      call getLimitsEoL(ml,s,qmin , smin,tmin,tmax,iret)
      if(iret.ne.0)return

      call getQminEoL(ml,qref1,qref2,qqref) !used for tabulation, don't change here
      call getSminEoL(ml,qqref,spminx)
      spminx=spminx*1.10

      call getQmaxEoL(ml,s,qmax)

      qmax=qmax / 1.01!1.10 
      qmax=max(qmax,1.01*qqref)
 
      if(qmin.gt.qmax)return

      if(log(qmax/qref1).ne.0.)then
        qli=log(q1/qref1)/log(qmax/qref1)*19.+1.
      else
        qli=1.0
      endif
      i=int(qli)
      qii=max(q1,qref2)
      if(log(qmax/qii).ne.0.)then
        qlj=log(qq/qii)/log(qmax/qii)*19.+1. 
      else
        qlj=1.0
      endif
      j=int(qlj)
      sl=log(s/spminx)/log(epmax/2./spminx)*19.+1.
      if(j.lt.1)j=1
      if(i.lt.1)i=1
      if(i.gt.18)i=18
      if(j.gt.18)j=18
      k=int(sl)
      if(k.lt.1)k=1
      if(k.gt.18)k=18
      if(csord(i,j,k+mlx).lt.-68.)k=k+1
      if(k.lt.1)k=1
      if(k.gt.18)k=18
      !~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~!TEST
      !write(ifmt,'(a,e9.3,5f9.1)')'---psjti1---',s,q1,qq,sl,qli,qlj !TEST
      !write(ifmt,'(a,3i9,i14)')'---psjti1---',k,i,j,mlx             !TEST
      !~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~!TEST

      vav=0
      nav=0
      do k1=1,3
      do jx=1,3
      do i1=1,3
        w=csord(i+i1-1,j+jx-1,k+k1+mlx-1)
        if(w.gt.-68.0)then
          vav=vav+w
          nav=nav+1
        endif
      enddo
      enddo
      enddo
      if(nav.gt.0)then
       vav=vav/nav
      endif   

      wi(2)=qli-i
      wi(3)=wi(2)*(wi(2)-1.)*.5
      wi(1)=1.-wi(2)+wi(3)
      wi(2)=wi(2)-2.*wi(3)

      wj(2)=qlj-j
      wj(3)=wj(2)*(wj(2)-1.)*.5
      wj(1)=1.-wj(2)+wj(3)
      wj(2)=wj(2)-2.*wj(3)

      wk(2)=sl-k
      wk(3)=wk(2)*(wk(2)-1.)*.5
      wk(1)=1.-wk(2)+wk(3)
      wk(2)=wk(2)-2.*wk(3)

      do i1=1,3
      do jx=1,3
      do k1=1,3
        k2=k+k1+mlx-1
        w=csord(i+i1-1,j+jx-1,k2)
        if(w.lt.-58..and.vav.ne.0.)w=vav !cross section value zero
        psjti1=psjti1+w*wi(i1)*wj(jx)*wk(k1)
      !~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~!TEST
      !qii=q2zmin*(qmax/q2zmin)**((  i+i1-1 -1)/19.)                 !TEST
      !qjj=qii*(qmax/qii)**((  j+jx-1  -1)/19.)                      !TEST
      !write(ifmt,'(a,3i5,2f9.3,4x,2f9.3)')'---psjti1---',k2,i+i1-1  !TEST
      !&,j+jx-1,w, wi(i1)*wj(jx)*wk(k1)                              !TEST
      !&,qii,qjj                                                     !TEST
      !~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~!TEST
        a=exp(psjti1) 
        if(.not.(a.gt.0..or.a.le.0.) !NaN catch
     .  .or.a.gt.1e35 .or. a.lt.-1e35 )then !Infinity catch
          b=wi(i1)
          c=wj(jx)
          d=wk(k1)
          e=w
          goto 998
        endif 
      enddo
      enddo
      enddo
      !~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~!TEST
      !psj1x=      psjet1(q1,q2,qqcut,s,j1,l1,0)                     !TEST 
      !psj1x=psj1x+psborn(q1,q2,qqcut,s,j1,l1,0,0)                   !TEST
      !sv=sngl(psuds(qq,l2)/psuds(q2,l2))                            !TEST
      !write(ifmt,'(a,2f14.7,$)')'---psjti1---',psjti1,log(psj1x/sv) !TEST
      !write(ifmt,'(2e14.3,$)')exp(psjti1),psj1x/sv                  !TEST
      !write(ifmt,'(2e14.3)')sv*exp(psjti1),psj1x                    !TEST
      !write(ifmt,'(80a1)')('-',m=1,80)                              !TEST  
      !~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~!TEST
      e=psjti1
      psjti1=exp(psjti1)

      if(jdis.eq.0.and.qq.gt.q2)then
        psjti1=psjti1*sngl(psuds(qq,l2)/psuds(q2,l2))
      elseif(jdis.eq.1.and.4.*qq.gt.q2)then
        psjti1=psjti1*sngl(psuds(4.*qq,l2)/psuds(q2,l2))
      endif

      a=psjti1
      b=99.99
      c=99.99 
      d=99.99
 998  continue
      call checkjtNAN('psjti1', psjti1,  m2,l2,  a,b,c,d,e
     .,q1,qref1,qq,qref2,qmax,s,  qli,qlj,sl,i,j,k,  smin,tmin,tmax )
      return
      end

c------------------------------------------------------------------------
      subroutine checkjtNAN(text, val,  m2,l2,  a,b,c,d,e
     .,q1,qref1,q2,qref2,qmax,s,  qli,qlj,sl,i,j,k,  smin,tmin,tmax )
c------------------------------------------------------------------------
#include "aaa.h"
      character text*6
c      dimension istop(1)
      data ncount/0/
      save ncount
      if(.not.(a.gt.0..or.a.le.0.) !NaN catch
     ..or.a.gt.1e35 .or. a.lt.-1e35 )then !Infinity catch
        ncount=ncount+1 
        if(l2.eq.5.or.m2.eq.5.or.l2.eq.4.or.m2.eq.4)then
          if(ncount.eq.1.or.ncount.eq.10.or.ncount.eq.100
     .      .or.ncount.eq.1000)then
            lcount=log10(1.*ncount)+1
            write(ifmt,'(a,i1,3a,2i2,2a)')'WARNING ',lcount,' in ',text
     .      ,' for flavors',m2,l2,', kinematics issue,'
     .      ,' check tmin,tmax,qmax'
            !see discussion at the beginning of psjti1
          endif
        else
          write(ifmt,*)'-----> WARNING ',text
     .    ,' ; a b c d e = ',a,b,c,d,e
          write(ifmt,*)'-----> q1 qref1 q2 qref2 qmax s = '
     .    ,q1,qref1,q2,qref2,qmax,s
          write(ifmt,*)'-----> qli qlj sl  = ',qli,qlj,sl
     .    ,'       flavors =',m2,l2
          write(ifmt,*)'----->  i   j   k  = ',i,j,k
          write(ifmt,*)'-----> smin tmin tmax = ',smin,tmin,tmax
        endif  
        val=0
c        nn=2
c        x=1/istop(nn)
      endif
      end

c------------------------------------------------------------------------
      function psjti(ii,q1,qqcut,s,m2,l2,jdis)
c-----------------------------------------------------------------------
#include "aaa.h"
      q2=q2cmin(3-ii)
      psjti=psjtoi(q1,q2,qqcut,s,m2,l2,jdis)
      end

c------------------------------------------------------------------------
      function psjtoi(q1,q2,qqcut,s,m2,l2,jdis)
c-----------------------------------------------------------------------
c psjtoi - inclusive hard cross-section interpolation - for any ordering
c in the ladder
c
c meaning of ii:
c    q2cmin(3-ii) is q2 on opposite side,
c      ii = 1  ...  current side = proj
c      ii = 2  ...  current side = targ
c
c q1 - virtuality cutoff at current end of the ladder
c qqcut - additional virtuality cutoff
c s  - total c.m. energy squared for the ladder
c m2 - parton type at current end of the ladder (0-g, 1,2,etc.-q)
c l2 - parton type at opposite end of the ladder (0-g, 1,2,etc.-q)
c-----------------------------------------------------------------------
      dimension wi(3),wj(3),wk(3)
      common /psar2/  edmax,epmax
#include "aaa.h"
#include "sem.h"
#include "tab.h"

      dummy=qqcut
  
      call eolpib( m2 , l2 , ml ) !KW1808c
      call eolpi( ml , j1 , l1 )                                    

      !~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~!TEST
      !write(ifmt,'(80a1/a,$)')('-',m=1,80),'---psjtoi---'         !TEST   
      !write(ifmt,'(3f8.2,e12.3,6i4)')q1,q2,qqcut,s,m2,l2,ml,j1,l1 !TEST 
      !~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~!TEST

      psjtoi=0.
      if(iabs(m2).eq.4.and. noflav(q1).lt.4 ) return
      if(iabs(m2).eq.5.and. noflav(q1).lt.5 ) return
      if(iabs(l2).eq.4.and. noflav(q2).lt.4 ) return
      if(iabs(l2).eq.5.and. noflav(q2).lt.5 ) return

      if(jdis.eq.0)ik=1
      if(jdis.ne.0)ik=4

      if(jdis.eq.0)qq=max(q1,q2)
      if(jdis.ne.0)qq=max(q1/4.,q2zmin) 

      call getLimitsEoL(ml,s,qq , smin,tmin,tmax,iret)
      if(iret.ne.0)return

      call getQminEoL(ml,qref1,qref2,qqref) !used for tabulation, don't change here
      if(jdis.ne.0)qqref=qqcut !???????
 
      call getSminEoL(ml,qqref,spmin)
      spmin=spmin*1.10

      call getQmaxEoL(ml,s,qmax)
      qmax=qmax *ik / 1.01!1.10 
      qmax=max(qmax,1.01*qqref)

      if(qmax.gt.qref1)then
        qli=log(q1/qref1)/log(qmax/qref1)*19.+1.
      else
        qli=1
      endif
      if(qmax.gt.qref2)then
        qlj=log(q2/qref2)/log(qmax/qref2)*19.+1.
      else
        qlj=1
      endif
      sl=log(s/spmin)/log(epmax/2./spmin)*19.+1.

      nip=3 !2 !3  !two or three point interpolation
      maxv=21-nip
      i=int(qli)
      i=max(1,min(i,maxv))
      j=int(qlj)
      j=max(1,min(j,maxv))
      k=int(sl)
      k=max(1,min(k,maxv))

      !~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~!TEST
      !write(ifmt,'(a,6f9.1)')'---psjtoi---',s,q1,q2,sl,qli,qlj      !TEST 
      !write(ifmt,'(a,3i9,i14)')'---psjtoi---',k,i,j                 !TEST
      !~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~!TEST

      vav=0
      nav=0
      do k1=1,nip
      do i1=1,nip
      do j1=1,nip
        w=cstot(i+i1-1,j+j1-1,k+k1-1,ml,jdis+1)
        if(w.gt.-68.0)then
          vav=vav+w
          nav=nav+1
        endif
      enddo
      enddo
      enddo
      if(nav.gt.0)then
       vav=vav/nav
      endif   

      wi(2)=qli-i
      wj(2)=qlj-j
      wk(2)=sl-k
      wi(3)=0
      wj(3)=0
      wk(3)=0

      if(nip.eq.2)then
        wi(1)=1-wi(2)
        wj(1)=1-wj(2)
        wk(1)=1-wk(2)
      else
        wi(3)=wi(2)*(wi(2)-1.)*.5
        wi(1)=1.-wi(2)+wi(3)
        wi(2)=wi(2)-2.*wi(3)
        wj(3)=wj(2)*(wj(2)-1.)*.5
        wj(1)=1.-wj(2)+wj(3)
        wj(2)=wj(2)-2.*wj(3)
        wk(3)=wk(2)*(wk(2)-1.)*.5
        wk(1)=1.-wk(2)+wk(3)
        wk(2)=wk(2)-2.*wk(3)
      endif

      do k1=1,nip
      do j1=1,nip
      do i1=1,nip
        a=cstot( i+i1-1, j+j1-1 , k+k1-1,ml,jdis+1) 
        if(a.lt.-58..and.vav.ne.0.)a=vav !cross section value zero
        psjtoi=psjtoi+a*wi(i1)*wj(j1)*wk(k1)
      !~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~!TEST 
      !write(ifmt,'(a,3i5,e14.3,2f9.3)')'---psjtoi---'               !TEST 
      !*,k+k1-1,i+i1-1,j+j1-1,wi(i1)*wj(j1)*wk(k1),a, psjtoi         !TEST
      !~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~!TEST
        a=exp(psjtoi) 
        if(.not.(a.gt.0..or.a.le.0.) !NaN catch
     .  .or.a.gt.1e35 .or. a.lt.-1e35 )then !Infinity catch
          b=wi(i1)
          c=wj(j1)
          d=wk(k1)
          e=cstot(i+i1-1, j+j1-1 , k+k1-1,ml,jdis+1)
          goto 998
        endif 
      enddo
      enddo
      enddo
      !~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~!TEST
      !psjx=       psjet(q1,q2,qq,s,j1,l1,0)                         !TEST
      !psjx=psjx + psjti1(q1,q2,qq,s,j2,l2,0)                        !TEST
      !psjx=psjx + psjti1(q2,q1,qq,s,l2,j2,0)                        !TEST
      !psjx=psjx - psboi(q1,q2,qq,s,j2,l2,0)                         !TEST
      !write(ifmt,'(a,2f14.7,$)')'---psjtoi---',psjtoi,log(psjx)     !TEST
      !write(ifmt,'(2e14.3)')exp(psjtoi),psjx                        !TEST
      !write(ifmt,'(80a1)')('-',m=1,80)                              !TEST  
      !~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~!TEST
      e=psjtoi
      psjtoi=exp(psjtoi)

      a=psjtoi
      b=99.99
      d=99.99
 998  continue
      call checkjtNAN('psjtoi', psjtoi,  m2,l2,  a,b,c,d,e
     .,q1,qref1,q2,qref2,qmax,s,  qli,qlj,sl,  i,j,k,  smin,tmin,tmax )
      end

c------------------------------------------------------------------------
      subroutine psjti0(ss,sj,sjb,m2,l2)
c-----------------------------------------------------------------------
c psjti0 - inclusive hard cross-section interpolation -
c for minimal virtuality cutoff in the ladder
c s - total c.m. energy squared for the ladder,
c sj - inclusive jet cross-section,
c sjb - born cross-section,
c m2 - parton type at current end of the ladder (0-g, 1,2,etc.-q)
c l2 - parton type at opposite end of the ladder (0-g, 1,2,etc.-q)
c-----------------------------------------------------------------------
#include "aaa.h"
#include "tab.h"
      dimension wk(3)
      common /psar2/  edmax,epmax
#include "sem.h"
      common/cpriom/npriom

      sj=0.
      sjb=0.

      call eolpib( m2 , l2 , ml ) !KW1808c
      qq=max(q2cmin(1),q2cmin(2))
      call getLimitsEoL(ml,ss,qq , smin,tmin,tmax,iret)
      if(iret.ne.0)return

      qqref=max(q2cmin(1),q2cmin(2))
      call getSminEoL(ml,qqref,spmin)

      sl=log(ss/spmin)/log(epmax/2./spmin)*19.+1.
      k=int(sl)
      if(k.lt.1)k=1 
      if(k.gt.18)k=18
      wk(2)=sl-k
      wk(3)=wk(2)*(wk(2)-1.)*.5
      wk(1)=1.-wk(2)+wk(3)
      wk(2)=wk(2)-2.*wk(3)

      do k1=1,3
        sj=sj+cstotzero(k+k1-1,ml)*wk(k1)
        sjb=sjb+csborzer(k+k1-1,ml)*wk(k1)
      enddo

      sjbxxx=sjb
      sjb=exp(sjb)
      sj=max(sjb,exp(sj))
                          !~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
      if(nancheck(sj).eq.1)then
        write(ifmt,'(a,2i4,$)')'ERROR sj b ; ',k,ml
        write(ifmt,*)ss,cstotzero(k,ml)
     .       ,cstotzero(k+1,ml),cstotzero(k+2,ml),sj
        stop
      endif
!~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
!the following activated by testOmEx
!~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
      if(npriom.eq.2.and.m2.eq.0.and.l2.eq.0)then
       !print*,'psjti0:',ss,sj,q2cmin
       !write(ifmt,'(a,3i4,$)')'psjti0:',k,m,l
       !write(ifmt,*)ss,cstotzero(k,m,l)
     .                          !,cstotzero(k+1,m,l),cstotzero(k+2,m,l),sj,sjb
        write(ifmt,'(a,3i4,$)')'psjti0:',k,ml
        write(ifmt,*)ss,csborzer(k,ml)
     .       ,csborzer(k+1,ml),csborzer(k+2,ml),sjb
      endif
                          !~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
      return
      end


c####################################################################################
c####################################################################################
c####################################################################################
c####################################################################################
c####
c####                additional tabulation / interpolation
c####
c####         2nd tabulation -> csjet, in addition to cstot,csord,csbor
c####
c####         Interpolation function    pijet   
c####
c####################################################################################
c####################################################################################
c####################################################################################
c####################################################################################


c-----------------------------------------------------------------------
      function pijet(ii,sk,m2,l2) !polynomial interpol of csjet tables  
c-----------------------------------------------------------------------
c  ii ..... type of emission
c           2 = bothside, 
c           1 = oneside (m2-side)
c          -1 = oneside (l2-side) 
c           0 = no emission, Born
c  qi ..... virtuality cutoff at current end of the ladder
c  qq ..... virtuality cutoff of Born
c  sk ..... energy squared for the scattering
c  m2,l2 .. parton types
c-----------------------------------------------------------------------
#include "aaa.h"
#include "sem.h"
#include "tab.h"
      common/psar2/edmax,epmax
      common/tabcsjet/ksmax,iqmax,jqmax
      common/cq2tmin/q2tmin(1000,2)
      real wk(3)!,wj(3),wi(3)
      data npijet/0/
      save npijet


c      q1=q2cmin(1)
c      q2=q2cmin(2)
c      qmx=max(q1,q2)
c      call eolpib( m2 , l2 , ml ) !KW1808b in pijet
c       call eolpi( ml , m1 , l1 ) !KW1808c     
c      if(ii.eq.2)then
c        pijet=psjet(q1,q2,qmx,sk,m1,l1,0)
c      elseif(ii.eq.1)then
c        pijet=psjet1(q1,q2,qmx,sk,m1,l1,0)
c      elseif(ii.eq.-1)then
c        pijet=psjet1(q2,q1,qmx,sk,m1,l1,0)
c      elseif(ii.eq.0)then
c        pijet=psborn(q1,q2,qmx,sk,m1,l1,0,0)
c      else
c        stop 'pijet'
c      endif
c      return
      
      nnn=0
      do  n=1,min(1000,npijet)
c      if(abs(q2cmin(1)-q2tmin(n,1)).lt.0.25
c     ..and.abs(q2cmin(2)-q2tmin(n,2)).lt.0.25 )nnn=n     !???? 0.25 can be decreased for better precision7
c???? 0.25 can be decreased for better precision7
        if(abs(max(q2cmin(1),q2cmin(2))-q2tmin(n,1)).lt.0.25)then
          nnn=n
          goto 1
        endif
      enddo

 1    if(nnn.eq.0)then
        mxpsar3xx=mxpsar3 !6 !instead of 20 channels
c        if(mxpsar3xx.gt.6)stop'Increase 6th dimension of csjet in tab.h'
        npijet=npijet+1
        call mkCSjet( mxpsar3xx,nt)
        nnn=nt
      endif

      call eolpib( m2 , l2 , ml ) !KW1808b in pijet
      if(ml.gt.mxpsar3xx)stop'Increase mxpsar3xx in pijet' 

      call getSminEoL(ml,q2zmin,spmin)
      spmed=spmin*(epmax/2./spmin)**(1./(ksmax-1.))
      if(sk.le.spmed)then
        kk=2
        spmax=spmed
      else
        kk=1
        spmax=epmax/2.
      endif

      sl= 1.+log(sk/spmin)/log(spmax/spmin)*(ksmax-1)
      k=int(sl)
      if(k.lt.1)k=1
      if(k.gt.(ksmax-2))k=ksmax-2

      wk(2)=sl-k
      wk(3)=wk(2)*(wk(2)-1.)*.5
      wk(1)=1.-wk(2)+wk(3)
      wk(2)=wk(2)-2.*wk(3)

      pijet=0
      do k1=1,3
        pijet=pijet
     .             +csjet(ii,kk,k+k1-1,1,1,ml,nnn)*wk(k1)
      enddo
      end

c-----------------------------------------------------------------------
      subroutine mkCSjet(mxpsar3xx,nt)     
c  former subroutine MakeCSTable, creates table csjet
c-----------------------------------------------------------------------
c tabulates psjet, psjet1, psborn for q2cmin(2)
c  2nd tabulation -> csjet, in addition to cstot,csord,csbor
c-----------------------------------------------------------------------
#include "aaa.h"
#include "sem.h"
#include "tab.h"
      common/psar2/edmax,epmax
      common/tabcsjet/ksmax,iqmax,jqmax
      common/cq2tmin/q2tmin(1000,2)
      data ncstable/0/
      save ncstable
      ncstable=ncstable+1
      ncstable=min(1000,ncstable)
      nt=ncstable

ctp      q2tmin(ncstable,1)=q2cmin(1)
ctp      q2tmin(ncstable,2)=q2cmin(2)
ctp only max of the 2 Q2s used in pijet ?????????????
      q2tmin(ncstable,1)=max(q2cmin(1),q2cmin(2))
      q2tmin(ncstable,2)=q2tmin(ncstable,1)

      write (ifmt,'(a,i5,1x,2f7.2,$)')'(CSjet ',ncstable
     . ,q2cmin(1),q2cmin(2)
ctp     . ,q2tmin(ncstable,1),q2tmin(ncstable,2)
      ksmax=10
      iqmax=3
      jqmax=3
      do ml=1,mxpsar3xx  !KW1808c 
       call eolpi( ml , m1 , l1 ) !KW1808c     
       call eolpic( m1 , l1 , m2 , l2 )
       call getSminEoL(ml,q2zmin,spmin)
       spmed=spmin*(epmax/2./spmin)**(1./(ksmax-1.))
       do kk=1,2
        if(kk.eq.1)spmax=epmax/2.
        if(kk.eq.2)spmax=spmed
        if(kk.eq.2.and.ml.eq.1)write (ifmt,'(a,$)')'.'
        do k=1,ksmax
          sk=spmin*(spmax/spmin)**((k-1)/(ksmax-1.))
          q1=q2cmin(1)
          q2=q2cmin(2)
          qmx=max(q1,q2)
          csjet(2,kk,k,1,1,ml,ncstable)
     .      =psjet(q1,q2,qmx,sk,m1,l1,0)
          csjet(1,kk,k,1,1,ml,ncstable)
     .      =psjet1(q1,q2,qmx,sk,m1,l1,0)
          csjet(-1,kk,k,1,1,ml,ncstable)
     .      =psjet1(q2,q1,qmx,sk,m1,l1,0) !emission other side, dont change order m1,l1 !!
          csjet(0,kk,k,1,1,ml,ncstable)
     .      =psborn(q1,q2,qmx,sk,m1,l1,0,0)
c         if(kk.eq.1.and.ml.eq.1)print*,'TEST',sk,
c     .    csjet(2,kk,k,1,1,ml,ncstable)+csjet(1,kk,k,1,1,ml,ncstable)
c     .   + csjet(-1,kk,k,1,1,ml,ncstable)+csjet(0,kk,k,1,1,ml,ncstable)
c         if(kk.eq.2.and.ml.eq.1)print*,'TEST',(' ',i=1,20),sk,
c     .    csjet(2,kk,k,1,1,ml,ncstable)+csjet(1,kk,k,1,1,ml,ncstable)
c     .   + csjet(-1,kk,k,1,1,ml,ncstable)+csjet(0,kk,k,1,1,ml,ncstable)
        enddo
       enddo
      enddo
      write (ifmt,'(a,$)')'done)'
      end

c-----------------------------------------------------------------------
      function pijetsum(f1,f2,f3,f4,sss,m2,l2) 
c-----------------------------------------------------------------------
      pijetsum=f1*pijet( 2 ,sss,m2,l2)  
     .        +f2*pijet( 1 ,sss,m2,l2)  
     .        +f3*pijet( -1,sss,m2,l2)  
     .        +f4*pijet( 0 ,sss,m2,l2)  
      end 

c-----------------------------------------------------------------------
      subroutine calcSTYPa(iii,f1,f2,f3,f4,sss)   !KW2006
c-----------------------------------------------------------------------
c computes sTYPz -> /zeroSTYP/
c-----------------------------------------------------------------------
      implicit none
#include "aaa.h"
#include "sem.h"
#include "tab.h"
      integer iii
      integer jl,j1,l1,j2,l2
      real f1,f2,f3,f4,sss,pijetsum,psborn
      real sTYPz
      common/zeroSTYP/sTYPz(mxpsar3)
      real qqq1,qqq2,qqq3
      common/cqqq/qqq1,qqq2,qqq3
      do jl=1,mxpsar3
        call eolpi( jl , j1 , l1 )      
        call eolpic( j1 , l1 , j2 , l2 ) 
        if(iii.ge.0.and.iii.le.3)then 
          sTYPz(jl) = pijetsum(f1,f2,f3,f4,sss,j2,l2)
        elseif(iii.eq.-1.or.iii.eq.10)then
          sTYPz(jl) =  psborn(qqq1,qqq2,qqq3,sss,j1,l1, -1 ,0)
        else
          print *,"calcSTYPa",iii
          stop'ERROR 10032019c' 
        endif
      enddo 
c     sub cascade threshold HQ production with correction factor (tuned to data)
      !n1=noflav(qqq1)
      !n2=noflav(qqq2)
      !if(facnof.gt.0.)then
      !...
      !else     KW2006 I removed this part, if needed, see version 3337
      !...
      !endif
      end



c####################################################################################
c####################################################################################
c####################################################################################
c####################################################################################
c####
c####              om51p for om5 calculation
c####
c####                 calls psvin for semihard (interpolation)
c####                 + soft contributions 
c####
c####              psvin = om5 interpolation 
c#### 
c####            using  tabulation of om5  (fhss,fhgg,fhqg,fhgq,fhqq)  in rj.i
c####
c####################################################################################
c####################################################################################
c####################################################################################
c####################################################################################


c-----------------------------------------------------------------------
      double precision function om51p(sy,xh,yp,b,iqq)
c-----------------------------------------------------------------------
c
c om51p = G_h1_h2 = D_h1_h2 * Fpart * Fpart in Phys. Rept
c
c     Fpart(x)=gamhad*x**betff
c     (betff=-alppar)      
c
c------------------------------------------------------------------------
c xh - fraction of the energy squared s for the pomeron;
c yp - rapidity for the pomeron;
c b - impact parameter between the pomeron ends;
c iqq = -1 - soft pomeron with lower energy threshold for diff,
c iqq = 0  - soft pomeron,
c iqq = 1  - gg,
c iqq = 2  - qg,
c iqq = 3  - gq,
c iqq = 4  - qq,
c iqq = 5 - Reggeon,
c iqq =11 - PomSat
c if isetcs or iqq < 0 use real functions instead of fit
c-----------------------------------------------------------------------
      double precision xh,yp,fac,xp,xm,xpp!,coefom1,coefom2
      common /psar7/  delx,alam3p,gam3p
c      common /psar37/ coefom1,coefom2
#include "aaa.h"
#include "sem.h"

      if(iqq.gt.12)stop'\n\n ERROR 20042018\n\n'

      om51p=0d0
c      write(ifch,*)'om51p',sy,xh,yp,b,s0min,iqq
      if(abs(iqq).eq.5)then
        if(sy.lt.s0min)return      !minimum energy for Reggeon or Diff
        if(iqq.gt.0)then
          if(iregge.eq.0)return !no reggeon
        else
          if(iomega.ge.2)return !no diff
        endif
      else
        if(sy.lt.1.)return         !implicit energy limit for Pomeron
        if(iqq.eq.0.and.isopom.eq.0)return !no soft
        if(iqq.gt.0.and.iqq.lt.5.and.ishpom.eq.0)return !no hard
      endif

      xp=sqrt(xh)*exp(yp)
      if(xh.gt.0.d0)then
        xm=xh/xp
      else
        xm=0.d0
      endif
      rp=r2had(iclpro)+r2had(icltar)+slopom*log(max(1.,sy))
      zb=exp(-b**2/(4.*.0389*rp))
      q2cmin(1)=q2pmin(1) 
      q2cmin(2)=q2pmin(2) 
ctp      q2cmin(1)=max(q2cmin(1),q2cmin(2))         !TP use symmetric Q2s ?????
ctp      q2cmin(2)=q2cmin(1)                        !TP use symmetric Q2s ?????

      if(((isetcs.lt.0..or.iqq.lt.0))!.or.(isetcs.lt.2..and.iqq.eq.11))
cKW210805      if(((isetcs.lt.0..or.iqq.lt.0).or.(isetcs.lt.2..and.iqq.eq.11)) !cKW210805
     .   .and.abs(iqq).ge.1.and.abs(iqq).ne.5)then
        call ipoBasics(q2cmin)
        call ipoCSzeroTables(q2cmin)
        if(abs(iqq).le.3)then
          iq=abs(iqq)-1
          if(iq.eq.0)then
            xpp=1.d0
            fac=dble(factk*sy**delh
     *     *alpff(1)*alpff(2))*xp**betff(1)
     *     *xm**betff(2) !Fpart*Fpart
          elseif(iq.eq.1)then
            xpp=xp
            fac=dble(factk*sy**delh
     *     *alpff(2))*xm**betff(2) !Fpart
          elseif(iq.eq.2)then
            xpp=xm
            fac=dble(factk*sy**delh
     *     *alpff(1))*xp**betff(1) !Fpart
          else
            fac=0.
          endif
          om51p=fac*dble(om11pp(sy,xpp,zb,iq)) !don't use table (for fast tests)
        elseif(iqq.eq.4)then
          rh=r2had(iclpro)+r2had(icltar)
c          rh=2.*r2part+(r2had(iclpro)-r2part)*(q2nmin/q2cmin(1))**1.
c     *         +(r2had(icltar)-r2part)*(q2nmin/q2cmin(2))**1. !tp?????????
          om51p=dble(om11pp3(sy,xp,xm) !don't use table (for fast tests)
     *          *zb**(rp/rh)
     *          *factk*sy**delh
     *          /(4.*pi*rh))
        elseif(iqq.eq.11.and.factsat.ge.0.)then
          iq=abs(iqq)-1
          fac=dble(factk*alpff(1)*alpff(2)*sy**delh)
     *       *(xp**betff(1)*xm**betff(2))
          om51p=fac*dble(om11pp4(sy,zb,iq)
     *           * factsat * (q2sft/max(q2cmin(1),q2cmin(2)))**alpsat )        
        endif
      elseif(iqq.eq.0.or.iqq.eq.-5)then      !soft
  !to define soft suppression in gamzz as x dependent
c        rp=2.*r2part!r2had(iclpro)+r2had(icltar)+0.5*log(max(1.,sy))
c        zb=exp(-b**2/(4.*.0389*rp))
        gamzz=1.
        call setParamEdep(q2cmin)                !define gamzz
        om51p=0.5d0*dble(2.*gampar**2
     *  *gamhad(iclpro)*gamhad(icltar)
     *  *zb/rp                 !no 0.0389 here because GeV used in omega
     *   *gamzz      !remove soft, related to zzsoft procedure
     *   )*dble(sy)**(alppom-1.)*xh**(-alppar) !(xp**betff(1)*xm**betff(2))
c if soft should depends on alpff and betff, ipoBasics should be call in ipoOm5Tables and ipoOm5Tab 
c       om51p=dble(2.*gampar**2
c     *  *alpff(1)*alpff(2)     
c     *  *zb/rp                 !no 0.0389 here because GeV used in omega
c     *   *gamzz      !remove soft, related to zzsoft procedure
c     *   )*dble(sy)**(0.5*(dels(1)+dels(2)))*(xp**betff(1)*xm**betff(2))
      elseif(iqq.le.4.or.iqq.eq.11)then      !gg,qg,gq,qq,ss
        om51p=psvin(sy,xp,xm,zb,iqq)
      elseif(iqq.eq.5)then      !Reggeon
        rr=r2reg(iclpro,ichproj)+r2reg(icltar,ichtarg)
     .    +sloreg(iclreg)*log(max(1.,sy/s0min))
        om51p=dble(gamreg(iclpro,ichproj)*gamreg(icltar,ichtarg)
     .       /(16.*pi)       !to be consistent with diff definition
     .       *(sy/s0min)**(alpreg(iclreg)-1.)*zb**(rp/rr)/rr)!no 0.0389 here because GeV used in omega
      endif

c      om51p=2.d0*om51p
c      om51p=0.5d0*om51p
c      if(om51p.gt.1e4.and.iqq.eq.11)then
c      write(ifmt,*)'om51p out',q2cmin,iqq,gampar,gamhad(icltar),zb,rp
c     .,xh,om51p
c        stop
c      endif

      return
      end


c-----------------------------------------------------------------------
      function psvin(sy,xpp,xpm,z,iqq)
c-----------------------------------------------------------------------
c   om5 interpolation using tables fhss,fhgg,fhgq,fhqg,fhqq
c-----------------------------------------------------------------------
c psvin - contributions to the interaction eikonal
c
c     G_h1_h2 = D_h1_h2 * Fpart * Fpart in Phys. Rept
c
c                   Fpart(x)=gamhad*x**(-alppar)
c------------------------------------------------------------------------
c sy  - energy squared for the hard interaction,
c xpp - lc+ for the sh pomeron,
c xpm - lc- for the sh pomeron,
c z   - impact parameter factor, z=exp(-b**2/4*rp),
c iqq = 1  - gg,
c iqq = 2  - qg,
c iqq = 3  - gq,
c iqq = 4  - qq,
c iqq = 11 - ss (saturated Pom)
c-----------------------------------------------------------------------
      dimension wk(3),wi(3,3),wj(3,3),wz(3),fa(3)
      common /psar2/  edmax,epmax
c      parameter (myom=7)
      common /psar4/  fhgg(11,100,10,8),fhqg(11,100,10,8)
      common /psar4b/ fhgq(11,100,10,8),fhqq(11,100,8)
      common /psar4c/ fhss(11,100,10,8)
      common /psar7/  delx,alam3p,gam3p
#include "aaa.h"
#include "sem.h"
      common/cpriom/npriom
      dimension i(3),j(3)
      double precision xpp,xpm,xp,xm

      if(iqq.gt.4.and.iqq.ne.11)stop'\n\n ERROR 07112012\n\n'

      psvin=0.
      spmin0=4.*q2zmin
      if(iqq.eq.1
     *.or.(iclpro.ne.4.and.iqq.eq.2)
     *.or.(icltar.ne.4.and.iqq.eq.3)
     *.or.(iclpro.ne.4.and.icltar.ne.4))then
        spmin=spmin0
        spmin2=4.*max(q2cmin(1),q2cmin(2)) !??????????????????????
      elseif(iqq.eq.11)then
        spmin=spmin0
        spmin2=4.*q2sft
      else
        spmin=spmin0+2.*qcmass**2
        spmin2=4.*max(q2cmin(1),q2cmin(2))+2.*qcmass**2
      endif
      spmin2=spmin2+(2.*ptipom+ptipomi*(log(q2cmin(1)/q2sft)
     .                                *log(q2cmin(2)/q2sft)))**2!increase limit to take into accout intrinsic pt but not taken into account in table calculation 
      if(sy.le.spmin2)return
      sy0=sngl(sy/xpp/xpm)

      yl=log(sy0/spmin)/log(epmax/2./spmin)*10.+1
      k=int(yl)
      if(k.gt.9)k=9
      wk(2)=yl-k
      wk(3)=wk(2)*(wk(2)-1.)*.5
      wk(1)=1.-wk(2)+wk(3)
      wk(2)=wk(2)-2.*wk(3)

      iclpt=iclpro+4*(icltar-1)

      npriom=0
      if(npriom.gt.0)then
        print*,'psvin (0) '
     .            ,sy,sy0,yl,k,yl-k,'   ',wk
        if(npriom.ge.2)write(ifmt,*)'psvin (0) '
     .            ,sy,sy0,spmin2,z,iqq,q2cmin
      endif

      do k1=1,3
        k2=k+k1-1  !use only intermediate energy to avoid problem with threshold at small x
        syx=(epmax/2./spmin0)**((k2-1)/10.) !cancellation *spmin0/spmin0
c        if(xpp.lt..25)then
        if(syx.gt.1.)then
          xl1=log(sngl(xpp*syx))/log(syx)*9.+1 !log(10.*xp)/log(2.)+5.
        else
          xl1=1.
c          xl1=4.*xpp+6.
        endif
        i(k1)=max(1,int(xl1))
c        if(i(k1).eq.6)i(k1)=5
c        if(i(k1).le.4)then
c          i(k1)=min(3,i(k1))
c        else
          i(k1)=min(8,i(k1))
c        endif
        del=xl1-i(k1)
cc      if(del.lt.-1)then
        if(del.lt.0.)then
c          wi(1)=((.1*2.**(-4))/xp)**0.1 !more stable extrapolation to x < xmin (growing for fixed sy) !1-del
c          wi(2)=0               !del
c          wi(3)=0
          wi(2,k1)=0.
          wi(3,k1)=0.
          wi(1,k1)=0.
        else
          wi(2,k1)=del
          wi(3,k1)=wi(2,k1)*(wi(2,k1)-1.)*.5
          wi(1,k1)=1.-wi(2,k1)+wi(3,k1)
          wi(2,k1)=wi(2,k1)-2.*wi(3,k1)
        endif
        if(npriom.gt.0)print*,'psvin (1) '
     .            ,xpp,k1,xl1,i(k1),xl1-i(k1),'   ',(wi(ii,k1),ii=1,3)

c       if(xpm.lt..25)then
        if(syx.gt.1.)then
          xl2=log(sngl(xpm*syx))/log(syx)*9.+1 !log(10.*xm)/log(2.)+5.
        else
          xl2=1.
c         xl2=4.*xpm+6.
        endif
          j(k1)=max(1,int(xl2))
c       if(j(k1).eq.6)j(k1)=5
c       if(j(k1).le.4)then
c         j(k1)=min(3,j(k1))
c       else
          j(k1)=min(8,j(k1))
c       endif
          del=xl2-j(k1)
c       if(del.lt.-1)then
        if(del.lt.0.)then
c      wj(1)=((.1*2.**(-4))/xm)**0.1 !more stable extrapolation to for x < xmin (growing for fixed sy))  !1-del
c      wj(2)=0 !del
c      wj(3)=0
          wi(2,k1)=0.
          wi(3,k1)=0.
          wi(1,k1)=0.
        else
          wj(2,k1)=del
          wj(3,k1)=wj(2,k1)*(wj(2,k1)-1.)*.5
          wj(1,k1)=1.-wj(2,k1)+wj(3,k1)
          wj(2,k1)=wj(2,k1)-2.*wj(3,k1)
        endif
        if(npriom.gt.0)print*,'psvin (2) '
     .            ,xpm,k1,xl2,j(k1),xl2-j(k1),'   ',(wj(jj,k1),jj=1,3)
      enddo

c for additionnal factor Fpar only
      if(iqq.eq.3)then
        xp=xpm
        xm=xpp
        icls=iclpro
      else
        xp=xpp
        xm=xpm
        icls=icltar
      endif

      rp=r2had(iclpro)+r2had(icltar)+slopom*log(max(1.,sy))


      if(iqq.ne.4)then  !---------------- not 4 ------------------


        jz=int(10.*z)
        if(jz.gt.8)jz=8
        if(jz.lt.1)jz=1
        wz(2)=10.*z-jz
        wz(3)=wz(2)*(wz(2)-1.)*.5
        wz(1)=1.-wz(2)+wz(3)
        wz(2)=wz(2)-2.*wz(3)
        if(npriom.gt.0)print*,'psvin (3) ',z,10.*z,jz,10.*z-jz,'   ',wz


        if(iqq.eq.1.or.iqq.eq.11)then   !1111111111111111111111111111111111

          do k1=1,3
            k2=k+k1-1
            fa(k1)=0
            do i1=1,3
            do j1=1,3
              i2=i(k1)+i1-1
              j2=j(k1)+j1-1
              ij2=i2+10*(j2-1)
              faz=0.
              if(iqq.eq.1)then
              if(fhgg(k2,ij2,jz  ,iclpt).gt.-60.)faz=faz
     .          +fhgg(k2,ij2,jz  ,iclpt)*wz(1)
              if(fhgg(k2,ij2,jz+1,iclpt).gt.-60.)faz=faz
     .          +fhgg(k2,ij2,jz+1,iclpt)*wz(2)
              if(fhgg(k2,ij2,jz+2,iclpt).gt.-60.)faz=faz
     .          +fhgg(k2,ij2,jz+2,iclpt)*wz(3)
              else
              if(fhss(k2,ij2,jz  ,iclpt).gt.-60.)faz=faz
     .          +fhss(k2,ij2,jz  ,iclpt)*wz(1)
              if(fhss(k2,ij2,jz+1,iclpt).gt.-60.)faz=faz
     .          +fhss(k2,ij2,jz+1,iclpt)*wz(2)
              if(fhss(k2,ij2,jz+2,iclpt).gt.-60.)faz=faz
     .          +fhss(k2,ij2,jz+2,iclpt)*wz(3)
              endif
              if(faz.ne.0.)then
c use special case with only one weight for x below table range.
c              if(i(3).le.10.or.j(3).le.10)then
                fa(k1)=fa(k1)+exp(max(-80.,min(80.,faz)))
     .               *wi(i1,k1)*wj(j1,k1)
c              else
c                fa(k1)=fa(k1)+faz*wi(i1,k1)*wj(j1,k1)
c              endif
              endif
                   !~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
                   if(npriom.eq.1.and.iqq.ne.11)then
                    print*,'psvin (iqq1) '
     *              ,k2,ij2,jz  ,iclpt,fhgg(k2,ij2,jz  ,iclpt)
     .              ,wz(1)*wi(i1,k1)*wj(j1,k1)
                    print*,'psvin (iqq1) '
     *              ,k2,ij2,jz+1,iclpt,fhgg(k2,ij2,jz+1,iclpt)
     .              ,wz(2)*wi(i1,k1)*wj(j1,k1)
                    print*,'psvin (iqq1) '
     *              ,k2,ij2,jz+2,iclpt,fhgg(k2,ij2,jz+2,iclpt)
     .              ,wz(3)*wi(i1,k1)*wj(j1,k1)
                    print *,'           ',fa(k1)
                    print*,' '
                   endif
                   if(npriom.eq.11.and.iqq.eq.11)then
                    print*,'psvin (iqq11) '
     *              ,k2,ij2,jz  ,iclpt,fhss(k2,ij2,jz  ,iclpt)
     .              ,wz(1)*wi(i1,k1)*wj(j1,k1)
                    print*,'psvin (iqq11) '
     *              ,k2,ij2,jz+1,iclpt,fhss(k2,ij2,jz+1,iclpt)
     .              ,wz(2)*wi(i1,k1)*wj(j1,k1)
                    print*,'psvin (iqq11) '
     *              ,k2,ij2,jz+2,iclpt,fhss(k2,ij2,jz+2,iclpt)
     .              ,wz(3)*wi(i1,k1)*wj(j1,k1)
                    print *,'           ',fa(k1)
                    print*,' '
                   endif
                  !~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
            enddo
            enddo
            if(fa(k1).eq.0..or.fa(k1).le.-80.)then
              fa(k1)=-80.
            else!if(i(3).le.10.or.j(3).le.10)then
              fa(k1)=log(max(1e-35,fa(k1)))
            endif
          if(npriom.eq.1)print*,'psvin (4)     ',k1,fa(k1),wk(k1)
          enddo
          if(fa(1)+fa(2)+fa(3).le.-60.)then
            psvin=max(0.,exp(max(-80.,min(80.,fa(1))))*wk(1)
     *                  +exp(max(-80.,min(80.,fa(2))))*wk(2)
     *                  +exp(max(-80.,min(80.,fa(3))))*wk(3))
          else
            psvin=exp(max(-80.,min(80.,fa(1)*wk(1)+fa(2)*wk(2)
     *                                            +fa(3)*wk(3))))
          endif
          !print*,'psvin (iqq1) ',fa,psvin
          psvin0=psvin
          psvin=psvin*z
     *     *factk*sy**delh                              !G
     *     *alpff(1)*alpff(2)*xp**betff(1)
     *     *xm**betff(2) !Fpart*Fpart
          
          if(npriom.eq.2)then
            write(ifch,*)'psvin (5) --> ',log(psvin0),psvin
            write(ifch,*)' '
          endif

          if(iqq.eq.11)psvin=psvin
     .         * factsat * (q2sft/max(q2cmin(1),q2cmin(2)))**alpsat
c          if(iqq.eq.11.and.factsat.lt.0.)psvin=0.

        else  ! 2222222222222222222222 3333333333333333333333 ....

          do k1=1,3
            k2=k+k1-1
            fa(k1)=0.
            do i1=1,3
            fak1=0.
            do j1=1,3
              i2=i(k1)+i1-1
              j2=j(k1)+j1-1
              ij2=i2+10*(j2-1)
              faz=0.
              if(iqq.eq.2)then
                if(fhqg(k2,ij2,jz  ,iclpt).gt.-60.)faz=faz
     *            +fhqg(k2,ij2,jz  ,iclpt)*wz(1)
                if(fhqg(k2,ij2,jz+1,iclpt).gt.-60.)faz=faz
     *            +fhqg(k2,ij2,jz+1,iclpt)*wz(2)
                if(fhqg(k2,ij2,jz+2,iclpt).gt.-60.)faz=faz
     *            +fhqg(k2,ij2,jz+2,iclpt)*wz(3)
              elseif(iqq.eq.3)then
                if(fhgq(k2,ij2,jz  ,iclpt).gt.-60.)faz=faz
     *            +fhgq(k2,ij2,jz  ,iclpt)*wz(1)
                if(fhgq(k2,ij2,jz+1,iclpt).gt.-60.)faz=faz
     *            +fhgq(k2,ij2,jz+1,iclpt)*wz(2)
                if(fhgq(k2,ij2,jz+2,iclpt).gt.-60.)faz=faz
     *            +fhgq(k2,ij2,jz+2,iclpt)*wz(3)
              endif
              if(faz.ne.0.)then
c              if(i(3).le.10.or.j(3).le.10)then
                fa(k1)=fa(k1)+exp(max(-80.,min(80.,faz)))
     .               *wi(i1,k1)*wj(j1,k1)
c              else
c                fa(k1)=fa(k1)+faz*wi(i1,k1)*wj(j1,k1)
c              endif
              endif
                   !~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
                   if(npriom.eq.1)then
                    print*,'psvin (iqq2) '
     *              ,k2,ij2,jz  ,iclpt,fhqg(k2,ij2,jz  ,iclpt)
     .              ,wz(1)*wi(i1,k1)*wj(j1,k1)
                    print*,'psvin (iqq2) '
     *              ,k2,ij2,jz+1,iclpt,fhqg(k2,ij2,jz+1,iclpt)
     .              ,wz(2)*wi(i1,k1)*wj(j1,k1)
                    print*,'psvin (iqq2) '
     *              ,k2,ij2,jz+2,iclpt,fhqg(k2,ij2,jz+2,iclpt)
     .              ,wz(3)*wi(i1,k1)*wj(j1,k1)
                    print *,'           ',fa(k1)
                    print*,' '
                   endif
                  !~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
            enddo
            enddo
            if(fa(k1).eq.0..or.fa(k1).le.-80.)then
              fa(k1)=-80.
            else!if(i(3).le.10.or.j(3).le.10)then
              fa(k1)=log(max(1e-35,fa(k1)))
            endif
            if(npriom.eq.1)print*,'psvin (4)     ',k1,fa(k1),wk(k1)
          enddo
          if(fa(1)+fa(2)+fa(3).le.-60.)then
            psvin=max(0.,exp(max(-80.,min(80.,fa(1))))*wk(1)
     *                  +exp(max(-80.,min(80.,fa(2))))*wk(2)
     *                  +exp(max(-80.,min(80.,fa(3))))*wk(3))
          else
            psvin=exp(max(-80.,min(80.,fa(1)*wk(1)+fa(2)*wk(2)
     *                                +fa(3)*wk(3))))
          endif
          psvin0=psvin
          psvin=psvin*z
     *         *factk*sy**delh
     *         *alpff(4-iqq)*xm**betff(4-iqq)                       !Fpart
          if(npriom.eq.2)then
            write(ifch,*)'psvin (5) --> ',log(psvin0),psvin
            write(ifch,*)' '
          endif

        endif

      else ! ------------- 4444444444444444444 -----------------------


        do k1=1,3
          k2=k+k1-1
          fa(k1)=0.
          do i1=1,3
          fak1=0.
          do j1=1,3
            i2=i(k1)+i1-1
            j2=j(k1)+j1-1
            ij2=i2+10*(j2-1)
            faz=0.
            if(fhqq(k2,ij2,iclpt).gt.-60.)faz=fhqq(k2,ij2,iclpt)
            if(faz.ne.0.)then
c              if(i(3).le.10.or.j(3).le.10)then
                fa(k1)=fa(k1)+exp(max(-80.,min(80.,faz)))
     .               *wi(i1,k1)*wj(j1,k1)
c              else
c                fa(k1)=fa(k1)+faz*wi(i1,k1)*wj(j1,k1)
c              endif
            endif
                   !~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
                   if(npriom.eq.1)then
                    print*,'psvin (iqq4) '
     *              ,k2,ij2,iclpt,fhqq(k2,ij2,iclpt)
     .              ,wi(i1,k1)*wj(j1,k1)
                    print *,'           ',fa(k1)
                    print*,' '
                   endif
                  !~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
          enddo
          enddo
          if(fa(k1).eq.0..or.fa(k1).le.-80.)then
            fa(k1)=-80.
          else!if(i(3).le.10.or.j(3).le.10)then
            fa(k1)=log(max(1e-35,fa(k1)))
          endif
        enddo
        if(fa(1)+fa(2)+fa(3).le.-60.)then
          psvin=max(0.,exp(max(-80.,min(80.,fa(1))))*wk(1)
     *                +exp(max(-80.,min(80.,fa(2))))*wk(2)
     *                +exp(max(-80.,min(80.,fa(3))))*wk(3))
        else
          psvin=exp(max(-80.,min(80.,fa(1)*wk(1)+fa(2)*wk(2)
     *                                          +fa(3)*wk(3))))
        endif

        r2hh=r2had(iclpro)+r2had(icltar)
c        r2hh=2.*r2part+(r2had(iclpro)-r2part)*(q2nmin/q2cmin(1))**1.
c     *                +(r2had(icltar)-r2part)*(q2nmin/q2cmin(2))**1.        !tp?????????
        psvin0=psvin
        psvin=psvin
     *          *z**(rp/r2hh)
     *          *factk*sy**delh
     *          /(4.*pi*r2hh)


          if(npriom.eq.2)then
            write(ifch,*)'psvin (5) --> ',log(psvin0),psvin
            write(ifch,*)' '
          endif
      endif      !--------------------------------------------

c        if(psvin.gt.1e10)then
c          print *,'ici',iqq,psvin,psvin0,fa,k
c        endif
      return
      end

c-----------------------------------------------------------------------
      subroutine getMassSquared( j2 , q2mass )  
c-----------------------------------------------------------------------
#include "aaa.h"
#include "sem.h"
      if(j2.lt.4)then 
        q2mass=0.
      elseif(j2.eq.4) then 
        q2mass=qcmass**2
      else             
        q2mass=qbmass**2        
      endif      
      end

c-----------------------------------------------------------------------
      subroutine getM2( j2 , l2 , q2mass )  !only for testing
c-----------------------------------------------------------------------
#include "aaa.h"
#include "sem.h"
      mxjl=max(abs(j2),abs(l2))
      if(mxjl.lt.4)then 
        q2mass=0.
      elseif(mxjl.eq.4) then 
        q2mass=qcmass**2
      else             
        q2mass=qbmass**2        
      endif      
      !used to define min s    
      end


