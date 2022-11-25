C
C  This file is part of EPOS4
C  Copyright (C) 2022 research institutions and authors (See CREDITS file)
C  This file is distributed under the terms of the GNU General Public License version 3 or later
C  (See COPYING file for the text of the licence)
C

c         ------------------------------------------------
c         since EPOS3086 we have om5 = om1
c                  (in earlier versions   om5 = 0.5 * om1)
c         change in :  xFitD1,xFitD2,xbExaD
c         ------------------------------------------------
c----------------------------------------------------------------------
      double precision function xDfit(zz,i1,i2,s,xp,xm,b)
c----------------------------------------------------------------------

#include "aaa.h"
#include "par.h"
#include "sem.h"



      xDfit=1.d0
      if(isetcs.gt.-2)then
      xDfit=0.d0
      do i=abs(max(0,i1)),abs(i2)
      call GfunPar(zz,zz,zz,zz,1,i,b,s,alp,bet,betp,epsp,epst,epss,gamv)
        if(i1.ge.0.and.i2.ge.0)then
          corp=alppar-epsp
          cort=alppar-epst
          cors=-epss
        else
          corp=alppar
          cort=alppar
          cors=0.
        endif
c        write(ifch,*)'dfit',i1,i2,i,alp,bet,bet+corp,epsp,cors
        xDfit=xDfit+dble(alp*xp**(bet+corp)*xm**(betp+cort))
     .             *dble(s/s0min)**cors
      enddo
      xDfit=abs(xDfit)
      endif
      return
      end


c----------------------------------------------------------------------
      subroutine xFitD1
c----------------------------------------------------------------------
c True omega with Q2s(x) fixed by E and b + some free value to test fit
c----------------------------------------------------------------------

#include "aaa.h"
#include "sem.h"
#include "par.h"
      parameter(nptg=50)  !number of point for the graphs
      double precision x,y,Dsoftshval,om51,xminr,tmp,xtmp,xDfit
     .,z(0:nptg)
      character chenergy*12

      iii=nint(xpar1)
      biniDf=xpar2                   !value of biniDf (impact parameter)
      y=dble(xpar3)                    !value of y (rapidity)
      xtmp=xmaxDf
      xmaxDf=dexp(-2.d0*y)
      q20=q2nmin
c      q2pmin(1)=q20+xpar98
c      q2pmin(2)=q20+xpar99
      q2pmin(1)=q20
      q2pmin(2)=q20
      call ipoOm5Tables(1)
      zz=xpar98
c to initialize screening variables
      tmp=xDfit(zz,idxD0,-idxD1,smaxDf,1.,1.,biniDf)

      chenergy='E=          '
      if (engy.ge.10000.) then
        write(chenergy(4:8),'(I5)')int(engy)
        ke=10
      elseif (engy.ge.1000.) then
        write(chenergy(4:7),'(I4)')int(engy)
        ke=9
      elseif (engy.ge.100.) then
        write(chenergy(4:6),'(I3)')int(engy)
        ke=8
      elseif (engy.ge.10.) then
        write(chenergy(4:5),'(I2)')int(engy)
        ke=7
      else
        write(chenergy(4:4),'(I1)')int(engy)
        ke=6
      endif
      chenergy(ke:ke+2)='GeV'

      xminr=xminDf  !value of xminr for plotting the function
      ominl=1d-3
      if(iregge.ne.0)ominl=min(ominl
     .                    ,0.05*sngl(Dsoftshval(10d0*xminr,y,biniDf,5)))

      if(iii/10.eq.1)then !...................................................

      write(ifhi,'(a)')'!----------------------------------------------'
      write(ifhi,'(a)')'!     D exact all      (blue)                  '
      write(ifhi,'(a)')'!----------------------------------------------'

      write(ifhi,'(a)')'openhisto name DExact-'//chenergy(4:ke-2)
      write(ifhi,'(a)')       'htyp lbu'
      write(ifhi,'(a)')       'xmod log ymod log'
      write(ifhi,'(a,2e11.3)')'xrange',xminr,xmaxDf
      write(ifhi,'(a,e11.3,a)')       'yrange ',ominl, ' auto'
      write(ifhi,'(a)')      'text 0 0 "xaxis x"'
      write(ifhi,'(a,a)')    'text 0 0 "yaxis D(x+,x-,s,b)" '
      write(ifhi,'(a,a)')    'text 0.65 0.9 "exact" '
      write(ifhi,'(3a)')          'text 0.05 0.9 "',chenergy,'"'
        write(ifhi,'(a,f5.2,a)')  'text 0.05 0.8 "b=',biniDf,' fm"'
        write(ifhi,'(a,f5.2,a)')  'text 0.05 0.7 "y=',y,'"'
        write(ifhi,'(a,2f6.1,a)')  'text 0.05 0.6 "Q?s!^2!=',q2pmin(1)
     *                                                  ,q2pmin(2),'"'
       write(ifhi,'(a)')       'array 2'

      do i=0,nptg
        x=xminr
        if (i.ne.0) x=x*(xmaxDf/xminr)**(dble(i)/dble(nptg))
        tmp=max(1.d-10,Dsoftshval(x,y,biniDf,0))
        write(ifhi,*) x,tmp
      enddo

      write(ifhi,'(a)')    '  endarray'
      write(ifhi,'(a)')    'closehisto plot 0-'


      write(ifhi,'(a)')'!----------------------------------------------'
      write(ifhi,'(a)')'!     D exact soft      (red dashed)           '
      write(ifhi,'(a)')'!----------------------------------------------'

      write(ifhi,'(a)')'openhisto name DExactSoft-'//chenergy(4:ke-2)
      write(ifhi,'(a)')       'htyp lra'
      write(ifhi,'(a)')       'xmod log ymod log'
      write(ifhi,'(a,2e11.3)')'xrange',xminr,xmaxDf
      write(ifhi,'(a)')       'yrange .01 auto'
      write(ifhi,'(a)')      'text 0 0 "xaxis x"'
      write(ifhi,'(a,a)')    'text 0 0 "yaxis D(x+,x-,s,b)" '
      write(ifhi,'(3a)')          'text 0.05 0.9 "',chenergy,'"'
      if (xpar8.eq.1.) then
        write(ifhi,'(a,f5.2,a)')  'text 0.05 0.8 "b=',biniDf,' fm"'
        write(ifhi,'(a,f5.2,a)')  'text 0.05 0.7 "y=',y,'"'
        write(ifhi,'(a,2f6.1,a)')  'text 0.05 0.6 "Q?s!^2!=',q2pmin(1)
     *                                                  ,q2pmin(2),'"'
      endif
      write(ifhi,'(a)')       'array 2'

      do i=0,nptg
        x=xminr
        if (i.ne.0) x=x*(xmaxDf/xminr)**(dble(i)/dble(nptg))
c        write(ifhi,*)x,max(1e-10,2.d0*om51(x,y,biniDf,0,0)
        write(ifhi,*)x,max(1d-10,om51(x,y,biniDf,0,0)
     &       /(x**dble(-alppar)))
      enddo

      write(ifhi,'(a)')    '  endarray'
      write(ifhi,'(a)')    'closehisto plot 0-'


      write(ifhi,'(a)')'!---------------------------------------------'
      write(ifhi,'(a)')'!     D exact sea-sea      (yellow-dashed)    '
      write(ifhi,'(a)')'!---------------------------------------------'

      write(ifhi,'(a)')'openhisto name DExactGluon-'//chenergy(4:ke-2)
      write(ifhi,'(a)')       'htyp lya'
      write(ifhi,'(a)')       'xmod log ymod log'
      write(ifhi,'(a,2e11.3)')'xrange',xminr,xmaxDf
      write(ifhi,'(a)')       'yrange .01 auto'
      write(ifhi,'(a)')      'text 0 0 "xaxis x"'
      write(ifhi,'(a,a)')    'text 0 0 "yaxis D(x+,x-,s,b)" '
      write(ifhi,'(3a)')          'text 0.05 0.9 "',chenergy,'"'
      if (xpar8.eq.1.) then
        write(ifhi,'(a,f5.2,a)')  'text 0.05 0.8 "b=',biniDf,' fm"'
        write(ifhi,'(a,f5.2,a)')  'text 0.05 0.7 "y=',y,'"'
        write(ifhi,'(a,2f6.1,a)')  'text 0.05 0.6 "Q?s!^2!=',q2pmin(1)
     *                                                  ,q2pmin(2),'"'
      endif
      write(ifhi,'(a)')       'array 2'

      do i=0,nptg
        x=xminr
        if (i.ne.0) x=x*(xmaxDf/xminr)**(dble(i)/dble(nptg))
c        tmp=max(1.e-10,2.d0*om51(x,y,biniDf,1,1)
        tmp=max(1.d-10,om51(x,y,biniDf,1,1)
     &       /(x**dble(-alppar)))
        write(ifhi,*) x,tmp
      enddo

      write(ifhi,'(a)')    '  endarray'

      write(ifhi,'(a)')    'closehisto plot 0-'

      write(ifhi,'(a)')'!---------------------------------------------'
      write(ifhi,'(a)')'!     D exact semi      (dashed blue)         '
      write(ifhi,'(a)')'!---------------------------------------------'

      write(ifhi,'(a)')'openhisto name DExactSemi-'//chenergy(4:ke-2)
      write(ifhi,'(a)')       'htyp lba'
      write(ifhi,'(a)')       'xmod log ymod log'
      write(ifhi,'(a,2e11.3)')'xrange',xminr,xmaxDf
      write(ifhi,'(a)')       'yrange .01 auto'
      write(ifhi,'(a)')      'text 0 0 "xaxis x"'
      write(ifhi,'(a,a)')    'text 0 0 "yaxis D(x+,x-,s,b)" '
      write(ifhi,'(3a)')          'text 0.05 0.9 "',chenergy,'"'
      if (xpar8.eq.1.) then
        write(ifhi,'(a,f5.2,a)')  'text 0.05 0.8 "b=',biniDf,' fm"'
        write(ifhi,'(a,f5.2,a)')  'text 0.05 0.7 "y=',y,'"'
        write(ifhi,'(a,2f6.1,a)')  'text 0.05 0.6 "Q?s!^2!=',q2pmin(1)
     *                                                  ,q2pmin(2),'"'
      endif
      write(ifhi,'(a)')       'array 2'

      do i=0,nptg
        x=xminr
        if (i.ne.0) x=x*(xmaxDf/xminr)**(dble(i)/dble(nptg))
        tmp=max(1.d-10,om51(x,y,biniDf,1,4)
     &                 /(x**dble(-alppar)))
        write(ifhi,*) x,tmp
      enddo

      write(ifhi,'(a)')    '  endarray'

      write(ifhi,'(a)')    'closehisto plot 0-'

      write(ifhi,'(a)')'!---------------------------------------------'
      write(ifhi,'(a)')'!     D Reggeon     (violet dash) '
      write(ifhi,'(a)')'!---------------------------------------------'

      write(ifhi,'(a)')'openhisto name DExactRegge-'//chenergy(4:ke-2)
      write(ifhi,'(a)')       'htyp lza'
      write(ifhi,'(a)')       'xmod log ymod log'
      write(ifhi,'(a,2e11.3)')'xrange',xminr,xmaxDf
      write(ifhi,'(a)')       'yrange .01 auto'
      write(ifhi,'(a)')      'text 0 0 "xaxis x"'
      write(ifhi,'(a,a)')    'text 0 0 "yaxis D(x+,x-,s,b)" '
      write(ifhi,'(3a)')          'text 0.05 0.9 "',chenergy,'"'
      if (xpar8.eq.1.) then
        write(ifhi,'(a,f5.2,a)')  'text 0.05 0.8 "b=',biniDf,' fm"'
        write(ifhi,'(a,f5.2,a)')  'text 0.05 0.7 "y=',y,'"'
        write(ifhi,'(a,2f6.1,a)')  'text 0.05 0.6 "Q?s!^2!=',q2pmin(1)
     *                                                  ,q2pmin(2),'"'
      endif
      write(ifhi,'(a)')       'array 2'

      do i=0,nptg
        x=xminr
        if (i.ne.0) x=x*(xmaxDf/xminr)**(dble(i)/dble(nptg))
        tmp=max(1.d-10,om51(x,y,biniDf,5,5)/(x**dble(-alppar)))
        write(ifhi,*) x,tmp
      enddo

      write(ifhi,'(a)')    '  endarray'
      if(iomega.lt.2)then
      write(ifhi,'(a)')    'closehisto plot 0-'

      write(ifhi,'(a)')'!---------------------------------------------'
      write(ifhi,'(a)')'!     D diff central      (green dashed) '
      write(ifhi,'(a)')'!---------------------------------------------'

      write(ifhi,'(a)')'openhisto name DExactDiffC-'//chenergy(4:ke-2)
      write(ifhi,'(a)')       'htyp lga'
      write(ifhi,'(a)')       'xmod log ymod log'
      write(ifhi,'(a,2e11.3)')'xrange',xminr,xmaxDf
      write(ifhi,'(a)')       'yrange .01 auto'
      write(ifhi,'(a)')      'text 0 0 "xaxis x"'
      write(ifhi,'(a,a)')    'text 0 0 "yaxis D(x+,x-,s,b)" '
      write(ifhi,'(3a)')          'text 0.05 0.9 "',chenergy,'"'
      if (xpar8.eq.1.) then
        write(ifhi,'(a,f5.2,a)')  'text 0.05 0.8 "b=',biniDf,' fm"'
        write(ifhi,'(a,f5.2,a)')  'text 0.05 0.7 "y=',y,'"'
        write(ifhi,'(a,2f6.1,a)')  'text 0.05 0.6 "Q?s!^2!=',q2pmin(1)
     *                                                  ,q2pmin(2),'"'
      endif
      write(ifhi,'(a)')       'array 2'

      do i=0,nptg
        x=xminr
        if (i.ne.0) x=x*(xmaxDf/xminr)**(dble(i)/dble(nptg))
        tmp=max(1.d-10,om51(x,y,biniDf,-8,-9)/XminDf**epssUni(0)
     *     /(x**dble(-alppar-epssUni(0))))
c        tmp=max(1.d-10,Dsoftshval(x,y,biniDf,3))
        write(ifhi,*) x,tmp
      enddo

      write(ifhi,'(a)')    '  endarray'
      write(ifhi,'(a)')    'closehisto plot 0-'

      write(ifhi,'(a)')'!---------------------------------------------'
      write(ifhi,'(a)')'!     D single diff      (green dash-dotted)  '
      write(ifhi,'(a)')'!---------------------------------------------'

      write(ifhi,'(a)')'openhisto name DExactSD-'//chenergy(4:ke-2)
      write(ifhi,'(a)')       'htyp lgi'
      write(ifhi,'(a)')       'xmod log ymod log'
      write(ifhi,'(a,2e11.3)')'xrange',xminr,xmaxDf
      write(ifhi,'(a)')       'yrange .01 auto'
      write(ifhi,'(a)')      'text 0 0 "xaxis x"'
      write(ifhi,'(a,a)')    'text 0 0 "yaxis D(x+,x-,s,b)" '
      write(ifhi,'(3a)')          'text 0.05 0.9 "',chenergy,'"'
      if (xpar8.eq.1.) then
        write(ifhi,'(a,f5.2,a)')  'text 0.05 0.8 "b=',biniDf,' fm"'
        write(ifhi,'(a,f5.2,a)')  'text 0.05 0.7 "y=',y,'"'
        write(ifhi,'(a,2f6.1,a)')  'text 0.05 0.6 "Q?s!^2!=',q2pmin(1)
     *                                                  ,q2pmin(2),'"'
      endif
      write(ifhi,'(a)')       'array 2'

      do i=0,nptg
        x=xminr
        if (i.ne.0) x=x*(xmaxDf/xminr)**(dble(i)/dble(nptg))
        tmp=max(1.d-10,om51(x,y,biniDf,-6,-7)/XminDf**epssUni(2)
     &       /(x**dble(-alppar-epssUni(2))))
        write(ifhi,*) x,tmp
      enddo

      write(ifhi,'(a)')    '  endarray'
      endif

      if(iii.eq.11)then
        write(ifhi,'(a)')    'closehisto plot 0-'
      else
        write(ifhi,'(a)')    'closehisto plot 0'
      endif

      endif !................................................................
      if(mod(iii,10).eq.1)then !.............................................

      write(ifhi,'(a)')'!---------------------------------------------'
      write(ifhi,'(a)')'!     D exact all      (blue)                 '
      write(ifhi,'(a)')'!---------------------------------------------'

      write(ifhi,'(a)')'openhisto name DExact-'//chenergy(4:ke-2)
      write(ifhi,'(a)')       'htyp lbu'
      write(ifhi,'(a)')       'xmod log ymod log'
      write(ifhi,'(a,2e11.3)')'xrange',xminr,xmaxDf
      write(ifhi,'(a,e11.3,a)')       'yrange ',ominl, ' auto'
      write(ifhi,'(a)')      'text 0 0 "xaxis x"'
      write(ifhi,'(a,a)')    'text 0 0 "yaxis D(x+,x-,s,b)" '
      write(ifhi,'(a,a)')    'text 0.65 0.9 "exact+fit" '
      write(ifhi,'(3a)')          'text 0.05 0.9 "',chenergy,'"'
      if (xpar8.eq.1.) then
        write(ifhi,'(a,f5.2,a)')  'text 0.05 0.8 "b=',biniDf,' fm"'
        write(ifhi,'(a,f5.2,a)')  'text 0.05 0.7 "y=',y,'"'
        write(ifhi,'(a,2f6.1,a)')  'text 0.05 0.6 "Q?s!^2!=',q2pmin(1)
     *                                                  ,q2pmin(2),'"'
      endif
      write(ifhi,'(a)')       'array 2'

      do i=0,nptg
        x=xminr
        if (i.ne.0) x=x*(xmaxDf/xminr)**(dble(i)/dble(nptg))
        z(i)=max(1.d-10,Dsoftshval(x,y,biniDf,0))
        write(ifhi,*) x,z(i)
      enddo

      write(ifhi,'(a)')    '  endarray'
      write(ifhi,'(a)')    'closehisto plot 0-'

      if(biniDf.le.0.)then

      write(ifhi,'(a)')'!---------------------------------------------'
      write(ifhi,'(a)')'!     D exact all -fit(hard+diff)(gray)       '
      write(ifhi,'(a)')'!---------------------------------------------'

      write(ifhi,'(a)')'openhisto name DExact-f-'//chenergy(4:ke-2)
      write(ifhi,'(a)')       'htyp lxu'
      write(ifhi,'(a)')       'xmod log ymod log'
      write(ifhi,'(a,2e11.3)')'xrange',xminr,xmaxDf
      write(ifhi,'(a)')       'yrange .01 auto'
      write(ifhi,'(a)')      'text 0 0 "xaxis x"'
      write(ifhi,'(a,a)')    'text 0 0 "yaxis D(x+,x-,s,b)" '
      write(ifhi,'(a,a)')    'text 0.65 0.9 "exact+fit" '
      write(ifhi,'(3a)')          'text 0.05 0.9 "',chenergy,'"'
      if (xpar8.eq.1.) then
        write(ifhi,'(a,f5.2,a)')  'text 0.05 0.8 "b=',biniDf,' fm"'
        write(ifhi,'(a,f5.2,a)')  'text 0.05 0.7 "y=',y,'"'
        write(ifhi,'(a,2f6.1,a)')  'text 0.05 0.6 "Q?s!^2!=',q2pmin(1)
     *                                                  ,q2pmin(2),'"'
      endif
      write(ifhi,'(a)')       'array 2'

      do i=0,nptg
        x=xminr
        if (i.ne.0) x=x*(xmaxDf/xminr)**(dble(i)/dble(nptg))
        tmp=max(1.d-10,abs(Dsoftshval(x,y,biniDf,-3)))!/z(i)
        write(ifhi,*) x,tmp
      enddo

      write(ifhi,'(a)')    '  endarray'
      write(ifhi,'(a)')    'closehisto plot 0-'
      endif

      write(ifhi,'(a)')'!---------------------------------------------'
      write(ifhi,'(a)')'!     fit soft      (red dot)             '
      write(ifhi,'(a)')'!---------------------------------------------'

      write(ifhi,'(a)')       'openhisto'
      write(ifhi,'(a)')       'htyp lzo'
      write(ifhi,'(a)')       'xmod log ymod log'
      write(ifhi,'(a,2e11.3)')'xrange',xminr,xmaxDf
      write(ifhi,'(a)')       'yrange .01 auto'
      write(ifhi,'(a)')      'text 0 0 "xaxis x"'
      write(ifhi,'(a,a)')    'text 0 0 "yaxis D(x+,x-,s,b)" '
      write(ifhi,'(3a)')          'text 0.05 0.9 "',chenergy,'"'
      write(ifhi,'(a,f5.2,a)')  'text 0.05 0.8 "b=',biniDf,' fm"'
      write(ifhi,'(a,f5.2,a)')  'text 0.05 0.7 "y=',y,'"'
      write(ifhi,'(a,2f6.1,a)')  'text 0.05 0.6 "Q?s!^2!=',q2pmin(1)
     *                                                ,q2pmin(2),'"'
      write(ifhi,'(a)')       'array 2'


      do i=0,nptg
        x=xminr
        if (i.ne.0) x=x*(xmaxDf/xminr)**(dble(i)/dble(nptg))
        xp=sqrt(sngl(x))*exp(sngl(y))
        xm=sqrt(sngl(x))*exp(-sngl(y))
        tmp=xDfit(zz,-1,idxD0,smaxDf,xp,xm,biniDf)
        write(ifhi,*) x,tmp
      enddo
      write(ifhi,'(a)')    '  endarray'
      write(ifhi,'(a)')    'closehisto plot 0-'

      write(ifhi,'(a)')'!---------------------------------------------'
      write(ifhi,'(a)')'!     fit semi      (blue dot)              '
      write(ifhi,'(a)')'!---------------------------------------------'

      write(ifhi,'(a)')       'openhisto'
      write(ifhi,'(a)')       'htyp lbo'
      write(ifhi,'(a)')       'xmod log ymod log'
      write(ifhi,'(a,2e11.3)')'xrange',xminr,xmaxDf
      write(ifhi,'(a)')       'yrange .01 auto'
      write(ifhi,'(a)')      'text 0 0 "xaxis x"'
      write(ifhi,'(a,a)')    'text 0 0 "yaxis D(x+,x-,,s,b)" '
      write(ifhi,'(3a)')          'text 0.05 0.9 "',chenergy,'"'
      write(ifhi,'(a,f5.2,a)')  'text 0.05 0.8 "b=',biniDf,' fm"'
      write(ifhi,'(a,f5.2,a)')  'text 0.05 0.7 "y=',y,'"'
      write(ifhi,'(a)')       'array 2'


      do i=0,nptg
        x=xminr
        if (i.ne.0) x=x*(xmaxDf/xminr)**(dble(i)/dble(nptg))
        xp=sqrt(sngl(x))*exp(sngl(y))
        xm=sqrt(sngl(x))*exp(-sngl(y))
        tmp=xDfit(zz,1,-1,smaxDf,xp,xm,biniDf)
        write(ifhi,*) x,tmp
      enddo
      write(ifhi,'(a)')    '  endarray'
      if(iomega.lt.2)then
      write(ifhi,'(a)')    'closehisto plot 0-'

      write(ifhi,'(a)')'!---------------------------------------------'
      write(ifhi,'(a)')'!     fit SD      (green dot)              '
      write(ifhi,'(a)')'!---------------------------------------------'

      write(ifhi,'(a)')       'openhisto'
      write(ifhi,'(a)')       'htyp lgo'
      write(ifhi,'(a)')       'xmod log ymod log'
      write(ifhi,'(a,2e11.3)')'xrange',xminr,xmaxDf
      write(ifhi,'(a)')       'yrange .01 auto'
      write(ifhi,'(a)')      'text 0 0 "xaxis x"'
      write(ifhi,'(a,a)')    'text 0 0 "yaxis D(x+,x-,,s,b)" '
      write(ifhi,'(3a)')          'text 0.05 0.9 "',chenergy,'"'
      write(ifhi,'(a,f5.2,a)')  'text 0.05 0.8 "b=',biniDf,' fm"'
      write(ifhi,'(a,f5.2,a)')  'text 0.05 0.7 "y=',y,'"'
      write(ifhi,'(a)')       'array 2'


      do i=0,nptg
        x=xminr
        if (i.ne.0) x=x*(xmaxDf/xminr)**(dble(i)/dble(nptg))
        xp=sqrt(sngl(x))*exp(sngl(y))
        xm=sqrt(sngl(x))*exp(-sngl(y))
        tmp=xDfit(zz,2,-3,smaxDf,xp,xm,biniDf)
        write(ifhi,*) x,tmp
      enddo
      write(ifhi,'(a)')    '  endarray'
      endif
      write(ifhi,'(a)')    'closehisto plot 0-'

      write(ifhi,'(a)')'!---------------------------------------------'
      write(ifhi,'(a)')'!     fit all      (red)      '
      write(ifhi,'(a)')'!---------------------------------------------'

      write(ifhi,'(a)')       'openhisto'
      write(ifhi,'(a)')       'htyp lru'
      write(ifhi,'(a)')       'xmod log ymod log'
      write(ifhi,'(a,2e11.3)')'xrange',xminr,xmaxDf
      write(ifhi,'(a)')       'yrange .01 auto'
      write(ifhi,'(a)')      'text 0 0 "xaxis x"'
      write(ifhi,'(a,a)')    'text 0 0 "yaxis D(x+,x-,,s,b)" '
      write(ifhi,'(3a)')          'text 0.05 0.9 "',chenergy,'"'
      write(ifhi,'(a,f5.2,a)')  'text 0.05 0.8 "b=',biniDf,' fm"'
      write(ifhi,'(a,f5.2,a)')  'text 0.05 0.7 "y=',y,'"'
      write(ifhi,'(a,2f6.1,a)')  'text 0.05 0.6 "Q?s!^2!=',q2pmin(1)
     *                                                ,q2pmin(2),'"'
      write(ifhi,'(a)')       'array 2'

      do i=0,nptg
        x=xminr
        if (i.ne.0) x=x*(xmaxDf/xminr)**(dble(i)/dble(nptg))
        xp=sqrt(sngl(x))*exp(sngl(y))
        xm=sqrt(sngl(x))*exp(-sngl(y))
        tmp=xDfit(zz,idxD0,-idxD1,smaxDf,xp,xm,biniDf)
        write(ifhi,*) x,tmp
      enddo
      write(ifhi,'(a)')    '  endarray'
      write(ifhi,'(a)')    'closehisto plot 0'

      endif !...............................................................

      xmaxDf=xtmp

      end

c----------------------------------------------------------------------
      subroutine xFitD2
c----------------------------------------------------------------------
c Omega with Q2(x) fixed independently of E,b at the table values used 
c for the fit.
c----------------------------------------------------------------------

#include "aaa.h"
#include "sem.h"
#include "par.h"
#include "tab.h"
      common/cpriom/npriom
      common /psar4/  fhgg(11,100,10,8),fhqg(11,100,10,8)

      double precision x,om51,om51p,xDfit,z(0:200),xminr,y,xtmp,om
c     & ,omYuncut,omNpuncut
      character chenergy*12,chf*3,texte*16,textb*18,texty*18

      !### kw ### ##########
      iimax=5
c      if(maproj.eq.1.and.matarg.eq.1)iimax=3
      if(iscreen.eq.0)iimax=1
      do ii=1,iimax            ! =========  ii loop ========= 
      if(ii.eq.1)then
        q2pmin(1)=q2mnval(1)
        q2pmin(2)=q2mnval(1)
      elseif(ii.eq.2)then
        q2pmin(1)=q2mnval(max(1,maxq2mx/4))
        q2pmin(2)=q2mnval(max(1,maxq2mx/4))
      elseif(ii.eq.3)then
        q2pmin(1)=q2mnval(maxq2mx/2)
        q2pmin(2)=q2mnval(maxq2mx/2)
      elseif(ii.eq.4)then
        q2pmin(1)=q2mnval(max(1,maxq2mx/4))
        q2pmin(2)=q2mnval(max(1,maxq2mx/2))
      else
        q2pmin(1)=q2mnval(max(1,maxq2mx/2))
        q2pmin(2)=q2mnval(max(1,maxq2mx/4))
      endif
c      call ipoOm5Tables(1)
      !### kw ### ##############

      nptg=30                  !number of point for the graphs
      biniDf=xpar2                 !value of biniDf (impact parameter)
      y=dble(xpar3)                 !value of y (rapidity)
      jj1=nint(xpar4)
      jj2=nint(xpar5)
      if(jj1.ne.1.and.jj1.ne.2)jj1=3
      if(jj2.ne.1.and.jj2.ne.2)jj2=3
      zz=0.
      xtmp=xmaxDf
      xmaxDf=dexp(-2.d0*y)

      chenergy='E=          '
      if (engy.ge.10000.) then
        write(chenergy(4:8),'(I5)')int(engy)
        ke=10
      elseif (engy.ge.1000.) then
        write(chenergy(4:7),'(I4)')int(engy)
        ke=9
      elseif (engy.ge.100.) then
        write(chenergy(4:6),'(I3)')int(engy)
        ke=8
      elseif (engy.ge.10.) then
        write(chenergy(4:5),'(I2)')int(engy)
        ke=7
      else
        write(chenergy(4:4),'(I1)')int(engy)
        ke=6
      endif
      chenergy(ke:ke+2)='GeV'

      xminr=1.d0/dble(engy**2)  !value of xminr for plotting the function
      ominl=max(1d-8,1e-2*om51(1d0,y,biniDf,0,0))!*xminr**dble(alppar))

         do jj=jj1,jj2

      if(jj.eq.1)chf='  D'
      if(jj.eq.2)chf='  G'
      if(jj.eq.3)chf='FFG'
      texte='text 0.44 0.84 "'
      textb='text 0.55 0.76 "b='
      texty='text 0.55 0.68 "y='
      if(jj.eq.2)texty='text 0.55 0.68 "y='
      if(jj.eq.3)texte='text 0.44 0.84 "'
      if(jj.eq.3)textb='text 0.55 0.76 "b='
      if(jj.eq.3)texty='text 0.55 0.68 "y='

      write(ifhi,'(a)')'!---------------------------------------------'
      write(ifhi,'(a)')'!     '//chf//' exact all      (green)         '
      write(ifhi,'(a)')'!---------------------------------------------'

      write(ifhi,'(a)')'openhisto name '//chf//'ExaI-'//chenergy(4:ke-2)
      write(ifhi,'(a)')       'htyp lgi'
      write(ifhi,'(a)')       'xmod log ymod log'
      write(ifhi,'(a,2e11.3)')'xrange',xminr,xmaxDf
      write(ifhi,'(a,e11.3,a)')'yrange ',ominl,' auto'
      write(ifhi,'(a)')      'text 0 0 "xaxis x"'
      write(ifhi,'(a)')      'text 0 0 "yaxis '//chf//'(x+,x-,s,b)" '
      write(ifhi,'(a,a)')    'text 0.55 0.12 "exact+fit" '
      write(ifhi,'(3a)')        texte,chenergy,'"'
      write(ifhi,'(a,f5.2,a)')  textb,biniDf,' fm"'
      write(ifhi,'(a,f5.2,a)')  texty,y,'"'
      write(ifhi,'(a,2f6.1,a)')  'text 0.33 0.92 "Q?s!^2!=',q2pmin(1)
     *                                                ,q2pmin(2),'"'
      write(ifhi,'(a)')       'array 2'
      do i=0,nptg
        x=xminr
        if (i.ne.0) x=x*(xmaxDf/xminr)**(dble(i)/dble(nptg))
        xp=sqrt(sngl(x))*exp(sngl(y))
        xm=sqrt(sngl(x))*exp(-sngl(y))
        call ipoOm5Tab(dble(xp),dble(xm),biniDf)
        om=om51(x,y,biniDf,0,1)
c        om=2.d0*om
        if(jj.eq.1)om=om/(x**dble(-alppar))
        if(jj.eq.3)om=om
     &               *(1-xm)**alplea(icltar)*(1-xp)**alplea(iclpro)
        write(ifhi,*) x,om
        !~~~~~~~~~~~~~~~~~~~~~~
         if(om.gt.1e16)then
         npriom=1 !provides printout in psvin
         print*,'555555',sngl(x*dble(engy**2)),sngl(x),sngl(y),biniDf,om
         xxx1=om51p(sngl(x*dble(engy**2)),x,y,biniDf,1)
         xxx2=om51p(sngl(x*dble(engy**2)),x,y,biniDf,2)
         xxx3=om51p(sngl(x*dble(engy**2)),x,y,biniDf,3)
         xxx4=om51p(sngl(x*dble(engy**2)),x,y,biniDf,4)
         xxx5=om51p(sngl(x*dble(engy**2)),x,y,biniDf,11)
         print*,xxx1,xxx2,xxx3,xxx4,xxx5
         stop
         endif  
        !~~~~~~~~~~~~~~~~~~~~~~
      enddo
      write(ifhi,'(a)')    '  endarray'
      write(ifhi,'(a)')    'closehisto plot 0-'

      write(ifhi,'(a)')'!---------------------------------------------'
      write(ifhi,'(a)')'!     '//chf//' gg  (black)        '
      write(ifhi,'(a)')'!---------------------------------------------'

      write(ifhi,'(a)')'openhisto name '//chf//'ExaI-'//chenergy(4:ke-2)
      write(ifhi,'(a)')       'htyp lkb'
      write(ifhi,'(a)')       'xmod log ymod log'
      write(ifhi,'(a,2e11.3)')'xrange',xminr,xmaxDf
      write(ifhi,'(a)')       'array 2'
      q2cmin(1)=q2nmin
      q2cmin(2)=q2nmin
      do i=0,nptg
        x=xminr
        if (i.ne.0) x=x*(xmaxDf/xminr)**(dble(i)/dble(nptg))
        xp=sqrt(sngl(x))*exp(sngl(y))
        xm=sqrt(sngl(x))*exp(-sngl(y))
        call ipoOm5Tab(dble(xp),dble(xm),biniDf)
        om=0
        do j=11,11
         om=om+om51p(sngl(x*dble(engy**2)),x,y,biniDf,j)
        enddo
c        om=2.d0*om
        if(jj.eq.1)om=om/(x**dble(-alppar))
        if(jj.eq.3)om=om
     &               *(1-xm)**alplea(icltar)*(1-xp)**alplea(iclpro)
        write(ifhi,*) x,om
      enddo
      write(ifhi,'(a)')    '  endarray'
      write(ifhi,'(a)')    'closehisto plot 0-'

      write(ifhi,'(a)')'!---------------------------------------------'
      write(ifhi,'(a)')'!     '//chf//' exact all +diff (blue)        '
      write(ifhi,'(a)')'!---------------------------------------------'

      write(ifhi,'(a)')'openhisto name '//chf//'ExaD-'//chenergy(4:ke-2)
      write(ifhi,'(a)')       'htyp lbu'
      write(ifhi,'(a)')       'xmod log ymod log'
      write(ifhi,'(a,2e11.3)')'xrange',xminr,xmaxDf
       write(ifhi,'(a)')       'array 2'
      do i=0,nptg
        x=xminr
        if (i.ne.0) x=x*(xmaxDf/xminr)**(dble(i)/dble(nptg))
        xp=sqrt(sngl(x))*exp(sngl(y))
        xm=sqrt(sngl(x))*exp(-sngl(y))
        call ipoOm5Tab(dble(xp),dble(xm),biniDf)
        om=om51(x,y,biniDf,0,5)+om51(x,y,biniDf,-6,-9)
        if(jj.eq.1)om=om/(x**dble(-alppar))
        if(jj.eq.3)om=om
     &               *(1-xm)**alplea(icltar)*(1-xp)**alplea(iclpro)
        write(ifhi,*) x,om
      enddo
      write(ifhi,'(a)')    '  endarray'
      write(ifhi,'(a)')    'closehisto plot 0-'

      write(ifhi,'(a)')'!---------------------------------------------'
      write(ifhi,'(a)')'!     '//chf//' param all      (red)          '
      write(ifhi,'(a)')'!---------------------------------------------'

      write(ifhi,'(a)')'openhisto name '//chf//'Par-'//chenergy(4:ke-2)
      write(ifhi,'(a)')       'htyp lra'
        write(ifhi,'(a)')       'xmod log ymod log'
      write(ifhi,'(a,2e11.3)')'xrange',xminr,xmaxDf
      write(ifhi,'(a)')       'array 2'
      do i=0,nptg
        x=xminr
        if (i.ne.0) x=x*(xmaxDf/xminr)**(dble(i)/dble(nptg))
        xp=sqrt(sngl(x))*exp(sngl(y))
        xm=sqrt(sngl(x))*exp(-sngl(y))
        z(i)=xDfit(zz,idxDmin(iomega),-idxDmax(iomega),engy**2,xp,xm
     .            ,biniDf)
         if(jj.ge.2)z(i)=z(i)*(x**dble(-alppar))
         if(jj.eq.3)z(i)=z(i)
     &               *(1-xm)**alplea(icltar)*(1-xp)**alplea(iclpro)
         write(ifhi,*) x,z(i)
       enddo
      write(ifhi,'(a)')    '  endarray'
      write(ifhi,'(a)')    'closehisto plot 0-'

      write(ifhi,'(a)')'!---------------------------------------------'
      write(ifhi,'(a)')'!     '//chf//' param semi      (yellow)       '
      write(ifhi,'(a)')'!---------------------------------------------'

      write(ifhi,'(a)')'openhisto name '//chf//'Scr-'//chenergy(4:ke-2)
      write(ifhi,'(a)')       'htyp lyb'
        write(ifhi,'(a)')       'xmod log ymod log'
      write(ifhi,'(a,2e11.3)')'xrange',xminr,xmaxDf
      write(ifhi,'(a)')       'array 2'
      do i=0,nptg
        x=xminr
        if (i.ne.0) x=x*(xmaxDf/xminr)**(dble(i)/dble(nptg))
        xp=sqrt(sngl(x))*exp(sngl(y))
        xm=sqrt(sngl(x))*exp(-sngl(y))
          z(i)=xDfit(zz,1,-1,engy**2,xp,xm,biniDf)
         if(jj.ge.2)z(i)=z(i)*(x**dble(-alppar))
         if(jj.eq.3)z(i)=z(i)
     &               *(1-xm)**alplea(icltar)*(1-xp)**alplea(iclpro)
         write(ifhi,*) x,z(i)
       enddo
      write(ifhi,'(a)')    '  endarray'
      write(ifhi,'(a)')    'closehisto plot 0'

         enddo

      xmaxDf=xtmp

      enddo ! =========  ii loop ========= 
      ! ### kw ###

      end

c----------------------------------------------------------------------
      subroutine xFitD3
c----------------------------------------------------------------------
c True omega with Q2s(x) fixed by E and b + some free value to test fit
c as a function of xp or xm
c----------------------------------------------------------------------

#include "aaa.h"
#include "sem.h"
#include "par.h"
      parameter(nptg=50)  !number of point for the graphs
      double precision x,y,Dsoftshval,om51,xminr,tmp,xDfit,xpm
     .,z(0:nptg)
      character chenergy*12

      iii=nint(xpar1)
      biniDf=xpar2                   !value of biniDf (impact parameter)
      xpm=dble(xpar3)
      if(xpm.ge.0.)then                    !choose between xp and xm
        sig=1.
      else
        xpm=abs(xpm)
        sig=-1.
      endif
      xpm=min(1d0,xpm)
      q20=q2nmin
c      q2pmin(1)=q20+xpar98
c      q2pmin(2)=q20+xpar99
      q2pmin(1)=q20
      q2pmin(2)=q20
      call ipoOm5Tables(1)
      zz=xpar98
c to initialize screening variables
      tmp=xDfit(zz,idxD0,-idxD1,smaxDf,1.,1.,biniDf)

      chenergy='E=          '
      if (engy.ge.10000.) then
        write(chenergy(4:8),'(I5)')int(engy)
        ke=10
      elseif (engy.ge.1000.) then
        write(chenergy(4:7),'(I4)')int(engy)
        ke=9
      elseif (engy.ge.100.) then
        write(chenergy(4:6),'(I3)')int(engy)
        ke=8
      elseif (engy.ge.10.) then
        write(chenergy(4:5),'(I2)')int(engy)
        ke=7
      else
        write(chenergy(4:4),'(I1)')int(engy)
        ke=6
      endif
      chenergy(ke:ke+2)='GeV'

      xminr=xminDf/xpm  !value of xminr for plotting the function
      y=sig*0.5*log(xminr/xpm)
      ominl=1d-3
      if(iregge.ne.0)ominl=min(ominl
     .               ,0.05*sngl(Dsoftshval(10d0*xminr*xpm,y,biniDf,5)))

      if(iii/10.eq.1)then !...................................................

      write(ifhi,'(a)')'!----------------------------------------------'
      write(ifhi,'(a)')'!     D exact all      (blue)                  '
      write(ifhi,'(a)')'!----------------------------------------------'

      write(ifhi,'(a)')'openhisto name DExact-'//chenergy(4:ke-2)
      write(ifhi,'(a)')       'htyp lbu'
      write(ifhi,'(a)')       'xmod log ymod log'
      write(ifhi,'(a,2e11.3)')'xrange',xminr,xmaxDf
      write(ifhi,'(a,e11.3,a)')       'yrange ',ominl, ' auto'
      write(ifhi,'(a)')      'text 0 0 "xaxis x"'
      if(sig.ge.0.)then
        write(ifhi,'(a)')      'text 0 0 "xaxis x+"'
        write(ifhi,'(a,f5.2,a)')  'text 0.05 0.7 "x-=',xpm,'"'
      else
        write(ifhi,'(a)')      'text 0 0 "xaxis x-"'
        write(ifhi,'(a,f5.2,a)')  'text 0.05 0.7 "x+=',xpm,'"'
      endif
      write(ifhi,'(a,a)')    'text 0 0 "yaxis D(x+,x-,s,b)" '
      write(ifhi,'(a,a)')    'text 0.65 0.9 "exact" '
      write(ifhi,'(3a)')          'text 0.05 0.9 "',chenergy,'"'
        write(ifhi,'(a,f5.2,a)')  'text 0.05 0.8 "b=',biniDf,' fm"'
        write(ifhi,'(a,2f6.1,a)')  'text 0.05 0.6 "Q?s!^2!=',q2pmin(1)
     *                                                  ,q2pmin(2),'"'
       write(ifhi,'(a)')       'array 2'

      do i=0,nptg
        x=xminr
        if (i.ne.0) x=x*(xmaxDf/xminr)**(dble(i)/dble(nptg))
        y=sig*0.5*log(x/xpm)
        tmp=max(1.d-10,Dsoftshval(x*xpm,y,biniDf,0))
        write(ifhi,*) x,tmp
      enddo

      write(ifhi,'(a)')    '  endarray'
      write(ifhi,'(a)')    'closehisto plot 0-'


      write(ifhi,'(a)')'!----------------------------------------------'
      write(ifhi,'(a)')'!     D exact soft      (red dashed)           '
      write(ifhi,'(a)')'!----------------------------------------------'

      write(ifhi,'(a)')'openhisto name DExactSoft-'//chenergy(4:ke-2)
      write(ifhi,'(a)')       'htyp lra'
      write(ifhi,'(a)')       'xmod log ymod log'
      write(ifhi,'(a,2e11.3)')'xrange',xminr,xmaxDf
      write(ifhi,'(a)')       'yrange .01 auto'
      if(sig.ge.0.)then
        write(ifhi,'(a)')      'text 0 0 "xaxis x+"'
        write(ifhi,'(a,f5.2,a)')  'text 0.05 0.7 "x-=',xpm,'"'
      else
        write(ifhi,'(a)')      'text 0 0 "xaxis x-"'
        write(ifhi,'(a,f5.2,a)')  'text 0.05 0.7 "x+=',xpm,'"'
      endif
      write(ifhi,'(a,a)')    'text 0 0 "yaxis D(x+,x-,s,b)" '
      write(ifhi,'(3a)')          'text 0.05 0.9 "',chenergy,'"'
      if (xpar8.eq.1.) then
        write(ifhi,'(a,f5.2,a)')  'text 0.05 0.8 "b=',biniDf,' fm"'
        write(ifhi,'(a,2f6.1,a)')  'text 0.05 0.6 "Q?s!^2!=',q2pmin(1)
     *                                                  ,q2pmin(2),'"'
      endif
      write(ifhi,'(a)')       'array 2'

      do i=0,nptg
        x=xminr
        if (i.ne.0) x=x*(xmaxDf/xminr)**(dble(i)/dble(nptg))
        y=sig*0.5*log(x/xpm)
c        write(ifhi,*)x,max(1e-10,2.d0*om51(x,y,biniDf,0,0)
        write(ifhi,*)x,max(1d-10,om51(x*xpm,y,biniDf,0,0)
     &       /((x*xpm)**dble(-alppar)))
      enddo

      write(ifhi,'(a)')    '  endarray'
      write(ifhi,'(a)')    'closehisto plot 0-'


      write(ifhi,'(a)')'!---------------------------------------------'
      write(ifhi,'(a)')'!     D exact sea-sea      (yellow-dashed)    '
      write(ifhi,'(a)')'!---------------------------------------------'

      write(ifhi,'(a)')'openhisto name DExactGluon-'//chenergy(4:ke-2)
      write(ifhi,'(a)')       'htyp lya'
      write(ifhi,'(a)')       'xmod log ymod log'
      write(ifhi,'(a,2e11.3)')'xrange',xminr,xmaxDf
      write(ifhi,'(a)')       'yrange .01 auto'
      if(sig.ge.0.)then
        write(ifhi,'(a)')      'text 0 0 "xaxis x+"'
        write(ifhi,'(a,f5.2,a)')  'text 0.05 0.7 "x-=',xpm,'"'
      else
        write(ifhi,'(a)')      'text 0 0 "xaxis x-"'
        write(ifhi,'(a,f5.2,a)')  'text 0.05 0.7 "x+=',xpm,'"'
      endif
      write(ifhi,'(a,a)')    'text 0 0 "yaxis D(x+,x-,s,b)" '
      write(ifhi,'(3a)')          'text 0.05 0.9 "',chenergy,'"'
      if (xpar8.eq.1.) then
        write(ifhi,'(a,f5.2,a)')  'text 0.05 0.8 "b=',biniDf,' fm"'
        write(ifhi,'(a,2f6.1,a)')  'text 0.05 0.6 "Q?s!^2!=',q2pmin(1)
     *                                                  ,q2pmin(2),'"'
      endif
      write(ifhi,'(a)')       'array 2'

      do i=0,nptg
        x=xminr
        if (i.ne.0) x=x*(xmaxDf/xminr)**(dble(i)/dble(nptg))
        y=sig*0.5*log(x/xpm)
c        tmp=max(1.e-10,2.d0*om51(x,y,biniDf,1,1)
        tmp=max(1.d-10,om51(x*xpm,y,biniDf,1,1)
     &       /((x*xpm)**dble(-alppar)))
        write(ifhi,*) x,tmp
      enddo

      write(ifhi,'(a)')    '  endarray'

      write(ifhi,'(a)')    'closehisto plot 0-'

      write(ifhi,'(a)')'!---------------------------------------------'
      write(ifhi,'(a)')'!     D exact semi      (dashed blue)         '
      write(ifhi,'(a)')'!---------------------------------------------'

      write(ifhi,'(a)')'openhisto name DExactSemi-'//chenergy(4:ke-2)
      write(ifhi,'(a)')       'htyp lba'
      write(ifhi,'(a)')       'xmod log ymod log'
      write(ifhi,'(a,2e11.3)')'xrange',xminr,xmaxDf
      write(ifhi,'(a)')       'yrange .01 auto'
      if(sig.ge.0.)then
        write(ifhi,'(a)')      'text 0 0 "xaxis x+"'
        write(ifhi,'(a,f5.2,a)')  'text 0.05 0.7 "x-=',xpm,'"'
      else
        write(ifhi,'(a)')      'text 0 0 "xaxis x-"'
        write(ifhi,'(a,f5.2,a)')  'text 0.05 0.7 "x+=',xpm,'"'
      endif
      write(ifhi,'(a,a)')    'text 0 0 "yaxis D(x+,x-,s,b)" '
      write(ifhi,'(3a)')          'text 0.05 0.9 "',chenergy,'"'
      if (xpar8.eq.1.) then
        write(ifhi,'(a,f5.2,a)')  'text 0.05 0.8 "b=',biniDf,' fm"'
        write(ifhi,'(a,2f6.1,a)')  'text 0.05 0.6 "Q?s!^2!=',q2pmin(1)
     *                                                  ,q2pmin(2),'"'
      endif
      write(ifhi,'(a)')       'array 2'

      do i=0,nptg
        x=xminr
        if (i.ne.0) x=x*(xmaxDf/xminr)**(dble(i)/dble(nptg))
        y=sig*0.5*log(x/xpm)
        tmp=max(1.d-10,om51(x*xpm,y,biniDf,1,4)
     &                 /((x*xpm)**dble(-alppar)))
        write(ifhi,*) x,tmp
      enddo

      write(ifhi,'(a)')    '  endarray'

      write(ifhi,'(a)')    'closehisto plot 0-'

      write(ifhi,'(a)')'!---------------------------------------------'
      write(ifhi,'(a)')'!     D Reggeon     (violet dash) '
      write(ifhi,'(a)')'!---------------------------------------------'

      write(ifhi,'(a)')'openhisto name DExactRegge-'//chenergy(4:ke-2)
      write(ifhi,'(a)')       'htyp lza'
      write(ifhi,'(a)')       'xmod log ymod log'
      write(ifhi,'(a,2e11.3)')'xrange',xminr,xmaxDf
      write(ifhi,'(a)')       'yrange .01 auto'
      if(sig.ge.0.)then
        write(ifhi,'(a)')      'text 0 0 "xaxis x+"'
        write(ifhi,'(a,f5.2,a)')  'text 0.05 0.7 "x-=',xpm,'"'
      else
        write(ifhi,'(a)')      'text 0 0 "xaxis x-"'
        write(ifhi,'(a,f5.2,a)')  'text 0.05 0.7 "x+=',xpm,'"'
      endif
      write(ifhi,'(a,a)')    'text 0 0 "yaxis D(x+,x-,s,b)" '
      write(ifhi,'(3a)')          'text 0.05 0.9 "',chenergy,'"'
      if (xpar8.eq.1.) then
        write(ifhi,'(a,f5.2,a)')  'text 0.05 0.8 "b=',biniDf,' fm"'
        write(ifhi,'(a,2f6.1,a)')  'text 0.05 0.6 "Q?s!^2!=',q2pmin(1)
     *                                                  ,q2pmin(2),'"'
      endif
      write(ifhi,'(a)')       'array 2'

      do i=0,nptg
        x=xminr
        if (i.ne.0) x=x*(xmaxDf/xminr)**(dble(i)/dble(nptg))
        y=sig*0.5*log(x/xpm)
        tmp=max(1d-10,om51(x*xpm,y,biniDf,5,5)/((x*xpm)**dble(-alppar)))
        write(ifhi,*) x,tmp
      enddo

      write(ifhi,'(a)')    '  endarray'
      if(iomega.lt.2)then
      write(ifhi,'(a)')    'closehisto plot 0-'

      write(ifhi,'(a)')'!---------------------------------------------'
      write(ifhi,'(a)')'!     D diff central      (green dashed) '
      write(ifhi,'(a)')'!---------------------------------------------'

      write(ifhi,'(a)')'openhisto name DExactDiffC-'//chenergy(4:ke-2)
      write(ifhi,'(a)')       'htyp lga'
      write(ifhi,'(a)')       'xmod log ymod log'
      write(ifhi,'(a,2e11.3)')'xrange',xminr,xmaxDf
      write(ifhi,'(a)')       'yrange .01 auto'
      if(sig.ge.0.)then
        write(ifhi,'(a)')      'text 0 0 "xaxis x+"'
        write(ifhi,'(a,f5.2,a)')  'text 0.05 0.7 "x-=',xpm,'"'
      else
        write(ifhi,'(a)')      'text 0 0 "xaxis x-"'
        write(ifhi,'(a,f5.2,a)')  'text 0.05 0.7 "x+=',xpm,'"'
      endif
      write(ifhi,'(a,a)')    'text 0 0 "yaxis D(x+,x-,s,b)" '
      write(ifhi,'(3a)')          'text 0.05 0.9 "',chenergy,'"'
      if (xpar8.eq.1.) then
        write(ifhi,'(a,f5.2,a)')  'text 0.05 0.8 "b=',biniDf,' fm"'
        write(ifhi,'(a,2f6.1,a)')  'text 0.05 0.6 "Q?s!^2!=',q2pmin(1)
     *                                                  ,q2pmin(2),'"'
      endif
      write(ifhi,'(a)')       'array 2'

      do i=0,nptg
        x=xminr
        if (i.ne.0) x=x*(xmaxDf/xminr)**(dble(i)/dble(nptg))
        y=sig*0.5*log(x/xpm)
        tmp=max(1.d-10,om51(x*xpm,y,biniDf,-8,-9)/XminDf**epssUni(0)
     .     /((x*xpm)**dble(-alppar-epssUni(0))))
c        tmp=max(1.d-10,Dsoftshval(x,y,biniDf,3))
        write(ifhi,*) x,tmp
      enddo

      write(ifhi,'(a)')    '  endarray'
      write(ifhi,'(a)')    'closehisto plot 0-'

      write(ifhi,'(a)')'!---------------------------------------------'
      write(ifhi,'(a)')'!     D single diff      (green dash-dotted)  '
      write(ifhi,'(a)')'!---------------------------------------------'

      write(ifhi,'(a)')'openhisto name DExactSD-'//chenergy(4:ke-2)
      write(ifhi,'(a)')       'htyp lgi'
      write(ifhi,'(a)')       'xmod log ymod log'
      write(ifhi,'(a,2e11.3)')'xrange',xminr,xmaxDf
      write(ifhi,'(a)')       'yrange .01 auto'
      if(sig.ge.0.)then
        write(ifhi,'(a)')      'text 0 0 "xaxis x+"'
        write(ifhi,'(a,f5.2,a)')  'text 0.05 0.7 "x-=',xpm,'"'
      else
        write(ifhi,'(a)')      'text 0 0 "xaxis x-"'
        write(ifhi,'(a,f5.2,a)')  'text 0.05 0.7 "x+=',xpm,'"'
      endif
      write(ifhi,'(a,a)')    'text 0 0 "yaxis D(x+,x-,s,b)" '
      write(ifhi,'(3a)')          'text 0.05 0.9 "',chenergy,'"'
      if (xpar8.eq.1.) then
        write(ifhi,'(a,f5.2,a)')  'text 0.05 0.8 "b=',biniDf,' fm"'
        write(ifhi,'(a,2f6.1,a)')  'text 0.05 0.6 "Q?s!^2!=',q2pmin(1)
     *                                                  ,q2pmin(2),'"'
      endif
      write(ifhi,'(a)')       'array 2'

      do i=0,nptg
        x=xminr
        if (i.ne.0) x=x*(xmaxDf/xminr)**(dble(i)/dble(nptg))
        y=sig*0.5*log(x/xpm)
        tmp=max(1.d-10,om51(x*xpm,y,biniDf,-6,-7)/XminDf**epssUni(2)
     &     /((x*xpm)**dble(-alppar-0.5*(epssUni(2)+epssUni(3)))))
        write(ifhi,*) x,tmp
      enddo

      write(ifhi,'(a)')    '  endarray'
      endif

      if(iii.eq.11)then
        write(ifhi,'(a)')    'closehisto plot 0-'
      else
        write(ifhi,'(a)')    'closehisto plot 0'
      endif

      endif !................................................................
      if(mod(iii,10).eq.1)then !.............................................

      write(ifhi,'(a)')'!---------------------------------------------'
      write(ifhi,'(a)')'!     D exact all      (blue)                 '
      write(ifhi,'(a)')'!---------------------------------------------'

      write(ifhi,'(a)')'openhisto name DExact-'//chenergy(4:ke-2)
      write(ifhi,'(a)')       'htyp lbu'
      write(ifhi,'(a)')       'xmod log ymod log'
      write(ifhi,'(a,2e11.3)')'xrange',xminr,xmaxDf
      write(ifhi,'(a,e11.3,a)')       'yrange ',ominl, ' auto'
      if(sig.ge.0.)then
        write(ifhi,'(a)')      'text 0 0 "xaxis x+"'
        write(ifhi,'(a,f5.2,a)')  'text 0.05 0.7 "x-=',xpm,'"'
      else
        write(ifhi,'(a)')      'text 0 0 "xaxis x-"'
        write(ifhi,'(a,f5.2,a)')  'text 0.05 0.7 "x+=',xpm,'"'
      endif
      write(ifhi,'(a,a)')    'text 0 0 "yaxis D(x+,x-,s,b)" '
      write(ifhi,'(a,a)')    'text 0.65 0.9 "exact+fit" '
      write(ifhi,'(3a)')          'text 0.05 0.9 "',chenergy,'"'
      if (xpar8.eq.1.) then
        write(ifhi,'(a,f5.2,a)')  'text 0.05 0.8 "b=',biniDf,' fm"'
        write(ifhi,'(a,2f6.1,a)')  'text 0.05 0.6 "Q?s!^2!=',q2pmin(1)
     *                                                  ,q2pmin(2),'"'
      endif
      write(ifhi,'(a)')       'array 2'

      do i=0,nptg
        x=xminr
        if (i.ne.0) x=x*(xmaxDf/xminr)**(dble(i)/dble(nptg))
        y=sig*0.5*log(x/xpm)
        z(i)=max(1.d-10,Dsoftshval(x*xpm,y,biniDf,0))
        write(ifhi,*) x,z(i)
      enddo

      write(ifhi,'(a)')    '  endarray'
      write(ifhi,'(a)')    'closehisto plot 0-'

      if(biniDf.le.0.)then

      write(ifhi,'(a)')'!---------------------------------------------'
      write(ifhi,'(a)')'!     D exact all -fit(hard+diff)(gray)       '
      write(ifhi,'(a)')'!---------------------------------------------'

      write(ifhi,'(a)')'openhisto name DExact-f-'//chenergy(4:ke-2)
      write(ifhi,'(a)')       'htyp lxu'
      write(ifhi,'(a)')       'xmod log ymod log'
      write(ifhi,'(a,2e11.3)')'xrange',xminr,xmaxDf
      write(ifhi,'(a)')       'yrange .01 auto'
      if(sig.ge.0.)then
        write(ifhi,'(a)')      'text 0 0 "xaxis x+"'
        write(ifhi,'(a,f5.2,a)')  'text 0.05 0.7 "x-=',xpm,'"'
      else
        write(ifhi,'(a)')      'text 0 0 "xaxis x-"'
        write(ifhi,'(a,f5.2,a)')  'text 0.05 0.7 "x+=',xpm,'"'
      endif
      write(ifhi,'(a,a)')    'text 0 0 "yaxis D(x+,x-,s,b)" '
      write(ifhi,'(a,a)')    'text 0.65 0.9 "exact+fit" '
      write(ifhi,'(3a)')          'text 0.05 0.9 "',chenergy,'"'
      if (xpar8.eq.1.) then
        write(ifhi,'(a,f5.2,a)')  'text 0.05 0.8 "b=',biniDf,' fm"'
        write(ifhi,'(a,2f6.1,a)')  'text 0.05 0.6 "Q?s!^2!=',q2pmin(1)
     *                                                  ,q2pmin(2),'"'
      endif
      write(ifhi,'(a)')       'array 2'

      do i=0,nptg
        x=xminr
        if (i.ne.0) x=x*(xmaxDf/xminr)**(dble(i)/dble(nptg))
        y=sig*0.5*log(x/xpm)
        tmp=max(1.d-10,abs(Dsoftshval(x*xpm,y,biniDf,-3)))!/z(i)
        write(ifhi,*) x,tmp
      enddo

      write(ifhi,'(a)')    '  endarray'
      write(ifhi,'(a)')    'closehisto plot 0-'
      endif

      write(ifhi,'(a)')'!---------------------------------------------'
      write(ifhi,'(a)')'!     fit soft      (red dot)             '
      write(ifhi,'(a)')'!---------------------------------------------'

      write(ifhi,'(a)')       'openhisto'
      write(ifhi,'(a)')       'htyp lzo'
      write(ifhi,'(a)')       'xmod log ymod log'
      write(ifhi,'(a,2e11.3)')'xrange',xminr,xmaxDf
      write(ifhi,'(a)')       'yrange .01 auto'
      if(sig.ge.0.)then
        write(ifhi,'(a)')      'text 0 0 "xaxis x+"'
        write(ifhi,'(a,f5.2,a)')  'text 0.05 0.7 "x-=',xpm,'"'
      else
        write(ifhi,'(a)')      'text 0 0 "xaxis x-"'
        write(ifhi,'(a,f5.2,a)')  'text 0.05 0.7 "x+=',xpm,'"'
      endif
      write(ifhi,'(a,a)')    'text 0 0 "yaxis D(x+,x-,s,b)" '
      write(ifhi,'(3a)')          'text 0.05 0.9 "',chenergy,'"'
      write(ifhi,'(a,f5.2,a)')  'text 0.05 0.8 "b=',biniDf,' fm"'
      write(ifhi,'(a,2f6.1,a)')  'text 0.05 0.6 "Q?s!^2!=',q2pmin(1)
     *                                                ,q2pmin(2),'"'
      write(ifhi,'(a)')       'array 2'


      do i=0,nptg
        x=xminr
        if (i.ne.0) x=x*(xmaxDf/xminr)**(dble(i)/dble(nptg))
        y=dble(sig*0.5)*log(x/xpm)
        xp=sqrt(sngl(x*xpm))*exp(sngl(y))
        xm=sqrt(sngl(x*xpm))*exp(-sngl(y))
        tmp=xDfit(zz,-1,idxD0,smaxDf,xp,xm,biniDf)
        write(ifhi,*) x,tmp
      enddo
      write(ifhi,'(a)')    '  endarray'
      write(ifhi,'(a)')    'closehisto plot 0-'

      write(ifhi,'(a)')'!---------------------------------------------'
      write(ifhi,'(a)')'!     fit semi      (blue dot)              '
      write(ifhi,'(a)')'!---------------------------------------------'

      write(ifhi,'(a)')       'openhisto'
      write(ifhi,'(a)')       'htyp lbo'
      write(ifhi,'(a)')       'xmod log ymod log'
      write(ifhi,'(a,2e11.3)')'xrange',xminr,xmaxDf
      write(ifhi,'(a)')       'yrange .01 auto'
      if(sig.ge.0.)then
        write(ifhi,'(a)')      'text 0 0 "xaxis x+"'
        write(ifhi,'(a,f5.2,a)')  'text 0.05 0.7 "x-=',xpm,'"'
      else
        write(ifhi,'(a)')      'text 0 0 "xaxis x-"'
        write(ifhi,'(a,f5.2,a)')  'text 0.05 0.7 "x+=',xpm,'"'
      endif
      write(ifhi,'(a,a)')    'text 0 0 "yaxis D(x+,x-,,s,b)" '
      write(ifhi,'(3a)')          'text 0.05 0.9 "',chenergy,'"'
      write(ifhi,'(a,f5.2,a)')  'text 0.05 0.8 "b=',biniDf,' fm"'
      write(ifhi,'(a)')       'array 2'


      do i=0,nptg
        x=xminr
        if (i.ne.0) x=x*(xmaxDf/xminr)**(dble(i)/dble(nptg))
        y=sig*0.5*log(x/xpm)
        xp=sqrt(sngl(x*xpm))*exp(sngl(y))
        xm=sqrt(sngl(x*xpm))*exp(-sngl(y))
        tmp=xDfit(zz,1,-1,smaxDf,xp,xm,biniDf)
        write(ifhi,*) x,tmp
      enddo
      write(ifhi,'(a)')    '  endarray'
      if(iomega.lt.2)then
      write(ifhi,'(a)')    'closehisto plot 0-'

      write(ifhi,'(a)')'!---------------------------------------------'
      write(ifhi,'(a)')'!     fit SD      (green dot)              '
      write(ifhi,'(a)')'!---------------------------------------------'

      write(ifhi,'(a)')       'openhisto'
      write(ifhi,'(a)')       'htyp lgo'
      write(ifhi,'(a)')       'xmod log ymod log'
      write(ifhi,'(a,2e11.3)')'xrange',xminr,xmaxDf
      write(ifhi,'(a)')       'yrange .01 auto'
      if(sig.ge.0.)then
        write(ifhi,'(a)')      'text 0 0 "xaxis x+"'
        write(ifhi,'(a,f5.2,a)')  'text 0.05 0.7 "x-=',xpm,'"'
      else
        write(ifhi,'(a)')      'text 0 0 "xaxis x-"'
        write(ifhi,'(a,f5.2,a)')  'text 0.05 0.7 "x+=',xpm,'"'
      endif
      write(ifhi,'(a,a)')    'text 0 0 "yaxis D(x+,x-,,s,b)" '
      write(ifhi,'(3a)')          'text 0.05 0.9 "',chenergy,'"'
      write(ifhi,'(a,f5.2,a)')  'text 0.05 0.8 "b=',biniDf,' fm"'
      write(ifhi,'(a)')       'array 2'


      do i=0,nptg
        x=xminr
        if (i.ne.0) x=x*(xmaxDf/xminr)**(dble(i)/dble(nptg))
        y=sig*0.5*log(x/xpm)
        xp=sqrt(sngl(x*xpm))*exp(sngl(y))
        xm=sqrt(sngl(x*xpm))*exp(-sngl(y))
        tmp=xDfit(zz,2,-3,smaxDf,xp,xm,biniDf)
        write(ifhi,*) x,tmp
      enddo
      write(ifhi,'(a)')    '  endarray'
      endif
      write(ifhi,'(a)')    'closehisto plot 0-'

      write(ifhi,'(a)')'!---------------------------------------------'
      write(ifhi,'(a)')'!     fit all      (red)      '
      write(ifhi,'(a)')'!---------------------------------------------'

      write(ifhi,'(a)')       'openhisto'
      write(ifhi,'(a)')       'htyp lru'
      write(ifhi,'(a)')       'xmod log ymod log'
      write(ifhi,'(a,2e11.3)')'xrange',xminr,xmaxDf
      write(ifhi,'(a)')       'yrange .01 auto'
      if(sig.ge.0.)then
        write(ifhi,'(a)')      'text 0 0 "xaxis x+"'
        write(ifhi,'(a,f5.2,a)')  'text 0.05 0.7 "x-=',xpm,'"'
      else
        write(ifhi,'(a)')      'text 0 0 "xaxis x-"'
        write(ifhi,'(a,f5.2,a)')  'text 0.05 0.7 "x+=',xpm,'"'
      endif
      write(ifhi,'(a,a)')    'text 0 0 "yaxis D(x+,x-,,s,b)" '
      write(ifhi,'(3a)')          'text 0.05 0.9 "',chenergy,'"'
      write(ifhi,'(a,f5.2,a)')  'text 0.05 0.8 "b=',biniDf,' fm"'
      write(ifhi,'(a,2f6.1,a)')  'text 0.05 0.6 "Q?s!^2!=',q2pmin(1)
     *                                                ,q2pmin(2),'"'
      write(ifhi,'(a)')       'array 2'

      do i=0,nptg
        x=xminr
        if (i.ne.0) x=x*(xmaxDf/xminr)**(dble(i)/dble(nptg))
        y=sig*0.5*log(x/xpm)
        xp=sqrt(sngl(x*xpm))*exp(sngl(y))
        xm=sqrt(sngl(x*xpm))*exp(-sngl(y))
        tmp=xDfit(zz,idxD0,-idxD1,smaxDf,xp,xm,biniDf)
        write(ifhi,*) x,tmp
      enddo
      write(ifhi,'(a)')    '  endarray'
      write(ifhi,'(a)')    'closehisto plot 0'

      endif !...............................................................

      end

c----------------------------------------------------------------------
      subroutine xbExaD
c----------------------------------------------------------------------

#include "aaa.h"
#include "sem.h"
#include "par.h"

      double precision x,y,Dsoftshval,om51p,z,xDfit!,omNpuncut


      nptg=50     !number of point for the graphs
      bmax=xpar2
      bmax=max(0.1,bmax)
             !value max of b (impact parameter)
      y=dble(xpar3)                    !value of y (rapidity)
      x=dble(xpar4)
      q20=q2nmin
      q2pmin(1)=q20+xpar98
      q2pmin(2)=q20+xpar99
      call ipoOm5Tables(1)
      zz=0
      if(x.le.1.d-12)x=1.d0

      if (engy.ge.10.) then
       if (engy.ge.100.) then
        if (engy.ge.1000.) then
         if (engy.ge.10000.) then
      write(ifhi,'(a,I5)')  'openhisto name DExactb-',int(engy)
         else
      write(ifhi,'(a,I4)')  'openhisto name DExactb-',int(engy)
        endif
        else
      write(ifhi,'(a,I3)')  'openhisto name DExactb-',int(engy)
       endif
       else
      write(ifhi,'(a,I2)')  'openhisto name DExactb-',int(engy)
      endif
      else
      write(ifhi,'(a,I1)')  'openhisto name DExactb-',int(engy)
      endif
      write(ifhi,'(a)')       'htyp lru'
      write(ifhi,'(a)')       'xmod lin ymod log'
      write(ifhi,'(a,2e11.3)')'xrange',0.,bmax
      write(ifhi,'(a)')       'yrange .000001 auto'
c      write(ifhi,'(a)')       'yrange auto auto'
      write(ifhi,'(a)')      'text 0 0 "xaxis b"'
      write(ifhi,'(a,a)')    'text 0 0 "yaxis D(x+,x-,s,b)" '
      write(ifhi,'(a,e8.2,a)')  'text 0.1 0.9 "E=',engy,' GeV"'
      if (xpar8.eq.1.) then
      write(ifhi,'(a,f5.2,a)')  'text 0.1 0.8 "x=',x,'"'
      write(ifhi,'(a,f5.2,a)')  'text 0.1 0.7 "y=',y,'"'
      endif
      write(ifhi,'(a)')       'array 2'

      do i=0,nptg
        b=float(i)/float(nptg)*bmax
        z=Dsoftshval(x,y,b,0)
        write(ifhi,*) b,z
      enddo

      write(ifhi,'(a)')    '  endarray'
      write(ifhi,'(a)')    'closehisto plot 0-'

c**********************************************************************

      if (engy.ge.10.) then
       if (engy.ge.100.) then
        if (engy.ge.1000.) then
         if (engy.ge.10000.) then
      write(ifhi,'(a,I5)')       'openhisto name DParamb-',int(engy)
         else
      write(ifhi,'(a,I4)')       'openhisto name DParamb-',int(engy)
        endif
        else
      write(ifhi,'(a,I3)')       'openhisto name DParamb-',int(engy)
       endif
       else
      write(ifhi,'(a,I2)')       'openhisto name DParamb-',int(engy)
      endif
      else
      write(ifhi,'(a,I1)')       'openhisto name DParamb-',int(engy)
      endif
      write(ifhi,'(a)')       'htyp lba'
      write(ifhi,'(a)')       'xmod lin ymod lin'
      write(ifhi,'(a,2e11.3)')'xrange',0.,bmax
c      write(ifhi,'(a)')       'yrange -.01 auto'
      write(ifhi,'(a)')       'yrange auto auto'
      write(ifhi,'(a)')      'text 0 0 "xaxis b"'
      write(ifhi,'(a,a)')    'text 0 0 "yaxis D(x+,x-,s,b)" '
      write(ifhi,'(a,e8.2,a)')  'text 0.1 0.9 "E=',engy,' GeV"'
      if (xpar8.eq.1.) then
      write(ifhi,'(a,f5.2,a)')  'text 0.1 0.8 "x=',x,'"'
      write(ifhi,'(a,f5.2,a)')  'text 0.1 0.7 "y=',y,'"'
      endif
      write(ifhi,'(a)')       'array 2'

      do i=0,nptg
        b=float(i)/float(nptg)*bmax
        z=xDfit(zz,idxD,idxD,smaxDf,
     &       sngl(dsqrt(x)*dexp(y)),sngl(dsqrt(x)*dexp(-y)),b)

        write(ifhi,*) b,z
      enddo

      write(ifhi,'(a)')    '  endarray'
      write(ifhi,'(a)')    'closehisto plot 0-'

c**********************************************************************

      if (engy.ge.10.) then
       if (engy.ge.100.) then
        if (engy.ge.1000.) then
         if (engy.ge.10000.) then
      write(ifhi,'(a,I5)')       'openhisto name DParamSoftb-',int(engy)
         else
      write(ifhi,'(a,I4)')       'openhisto name DParamSoftb-',int(engy)
        endif
        else
      write(ifhi,'(a,I3)')       'openhisto name DParamSoftb-',int(engy)
       endif
       else
      write(ifhi,'(a,I2)')       'openhisto name DParamSoftb-',int(engy)
      endif
      else
      write(ifhi,'(a,I1)')       'openhisto name DParamSoftb-',int(engy)
      endif
      write(ifhi,'(a)')       'htyp lbo'
      write(ifhi,'(a)')       'xmod lin ymod lin'
      write(ifhi,'(a,2e11.3)')'xrange',0.,bmax
c      write(ifhi,'(a)')       'yrange -.01 auto'
      write(ifhi,'(a)')       'yrange auto auto'
      write(ifhi,'(a)')      'text 0 0 "xaxis b"'
      write(ifhi,'(a,a)')    'text 0 0 "yaxis D(x+,x-,s,b)" '
      write(ifhi,'(a,e8.2,a)')  'text 0.1 0.9 "E=',engy,' GeV"'
      if (xpar8.eq.1.) then
      write(ifhi,'(a,f5.2,a)')  'text 0.1 0.8 "x=',x,'"'
      write(ifhi,'(a,f5.2,a)')  'text 0.1 0.7 "y=',y,'"'
      endif
      write(ifhi,'(a)')       'array 2'

      do i=0,nptg
        b=float(i)/float(nptg)*bmax
        z=xDfit(zz,idxD0,idxD0,smaxDf,
     &       sngl(dsqrt(x)*dexp(y)),sngl(dsqrt(x)*dexp(-y)),b)

        write(ifhi,*) b,z
      enddo

      write(ifhi,'(a)')    '  endarray'
      write(ifhi,'(a)')    'closehisto plot 0-'


c**********************************************************************

      if (engy.ge.10.) then
       if (engy.ge.100.) then
        if (engy.ge.1000.) then
         if (engy.ge.10000.) then
      write(ifhi,'(a,I5)')  'openhisto name DExactSoftb-',int(engy)
         else
      write(ifhi,'(a,I4)')  'openhisto name DExactSoftb-',int(engy)
        endif
        else
      write(ifhi,'(a,I3)')  'openhisto name DExactSoftb-',int(engy)
       endif
       else
      write(ifhi,'(a,I2)')  'openhisto name DExactSoftb-',int(engy)
      endif
      else
      write(ifhi,'(a,I1)')  'openhisto name DExactSoftb-',int(engy)
      endif
      write(ifhi,'(a)')       'htyp pfc'
      write(ifhi,'(a)')       'xmod lin ymod lin'
      write(ifhi,'(a,2e11.3)')'xrange',0.,bmax
c      write(ifhi,'(a)')       'yrange -.01 auto'
      write(ifhi,'(a)')       'yrange auto auto'
      write(ifhi,'(a)')      'text 0 0 "xaxis b"'
      write(ifhi,'(a,a)')    'text 0 0 "yaxis D(x+,x-,s,b)" '
      write(ifhi,'(a,e8.2,a)')  'text 0.1 0.9 "E=',engy,' GeV"'
      if (xpar8.eq.1.) then
      write(ifhi,'(a,f5.2,a)')  'text 0.1 0.8 "x=',x,'"'
      write(ifhi,'(a,f5.2,a)')  'text 0.1 0.7 "y=',y,'"'
      endif
      write(ifhi,'(a)')       'array 2'

      do i=0,nptg
        b=float(i)/float(nptg)*bmax
c        z=2.d0*om51p(sngl(x*dble(smaxDf)),x,y,b,0)
        z=om51p(sngl(x*dble(smaxDf)),x,y,b,0)
     &       /(x**dble(-alppar))
        write(ifhi,*) b,z
      enddo

      write(ifhi,'(a)')    '  endarray'
      write(ifhi,'(a)')    'closehisto plot 0-'


c**********************************************************************

      if (engy.ge.10.) then
       if (engy.ge.100.) then
        if (engy.ge.1000.) then
         if (engy.ge.10000.) then
      write(ifhi,'(a,I5)')  'openhisto name DExactSemib-',int(engy)
         else
      write(ifhi,'(a,I4)')  'openhisto name DExactSemib-',int(engy)
        endif
        else
      write(ifhi,'(a,I3)')  'openhisto name DExactSemib-',int(engy)
       endif
       else
      write(ifhi,'(a,I2)')  'openhisto name DExactSemib-',int(engy)
      endif
      else
      write(ifhi,'(a,I1)')  'openhisto name DExactSemib-',int(engy)
      endif
      write(ifhi,'(a)')       'htyp pft'
      write(ifhi,'(a)')       'xmod lin ymod lin'
      write(ifhi,'(a,2e11.3)')'xrange',0.,bmax
c      write(ifhi,'(a)')       'yrange -.01 auto'
      write(ifhi,'(a)')       'yrange auto auto'
      write(ifhi,'(a)')      'text 0 0 "xaxis b"'
      write(ifhi,'(a,a)')    'text 0 0 "yaxis D(x+,x-,s,b)" '
      write(ifhi,'(a,e8.2,a)')  'text 0.1 0.9 "E=',engy,' GeV"'
      if (xpar8.eq.1.) then
      write(ifhi,'(a,f5.2,a)')  'text 0.1 0.8 "x=',x,'"'
      write(ifhi,'(a,f5.2,a)')  'text 0.1 0.7 "y=',y,'"'
      endif
      write(ifhi,'(a)')       'array 2'

      do i=0,nptg
        b=float(i)/float(nptg)*bmax
c        z=2.d0*om51p(sngl(x*dble(smaxDf)),x,y,b,1)
        z=om51p(sngl(x*dble(smaxDf)),x,y,b,1)
     &       /(x**dble(-alppar))
        write(ifhi,*) b,z
      enddo

      write(ifhi,'(a)')    '  endarray'
      write(ifhi,'(a)')    'closehisto plot 0-'


c**********************************************************************

      if (engy.ge.10.) then
       if (engy.ge.100.) then
        if (engy.ge.1000.) then
         if (engy.ge.10000.) then
      write(ifhi,'(a,I5)')  'openhisto name DExactValb-',int(engy)
         else
      write(ifhi,'(a,I4)')  'openhisto name DExactValb-',int(engy)
        endif
        else
      write(ifhi,'(a,I3)')  'openhisto name DExactValb-',int(engy)
       endif
       else
      write(ifhi,'(a,I2)')  'openhisto name DExactValb-',int(engy)
      endif
      else
      write(ifhi,'(a,I1)')  'openhisto name DExactValb-',int(engy)
      endif
      write(ifhi,'(a)')       'htyp pfs'
      write(ifhi,'(a)')       'xmod lin ymod lin'
      write(ifhi,'(a,2e11.3)')'xrange',0.,bmax
c      write(ifhi,'(a)')       'yrange -.01 auto'
      write(ifhi,'(a)')       'yrange auto auto'
      write(ifhi,'(a)')      'text 0 0 "xaxis b"'
      write(ifhi,'(a,a)')    'text 0 0 "yaxis D(x+,x-,s,b)" '
      write(ifhi,'(a,e8.2,a)')  'text 0.1 0.9 "E=',engy,' GeV"'
      if (xpar8.eq.1.) then
      write(ifhi,'(a,f5.2,a)')  'text 0.1 0.8 "x=',x,'"'
      write(ifhi,'(a,f5.2,a)')  'text 0.1 0.7 "y=',y,'"'
      endif
      write(ifhi,'(a)')       'array 2'

      do i=0,nptg
        b=float(i)/float(nptg)*bmax
c        z=2.d0*(om51p(sngl(x*dble(smaxDf)),x,y,b,2)+
        z=(om51p(sngl(x*dble(smaxDf)),x,y,b,2)+
     &  om51p(sngl(x*dble(smaxDf)),x,y,b,3)
     & +om51p(sngl(x*dble(smaxDf)),x,y,b,4)
     & +om51p(sngl(x*dble(smaxDf)),x,y,b,11))
     &           /(x**dble(-alppar))
        write(ifhi,*) b,z
      enddo

      write(ifhi,'(a)')    '  endarray'

      write(ifhi,'(a)')    'closehisto plot 0'


      end


c----------------------------------------------------------------------
      subroutine xbnExaD
c----------------------------------------------------------------------

#include "aaa.h"
#include "par.h"

      double precision x,y,Dsoftshval,z,xDfit,Dint



      nptg=50     !number of point for the graphs
      bmax=xpar2
      bmax=max(0.1,bmax)
      y=dble(xpar2)                   !value of y (rapidity)
      x=dble(xpar3)
      q2pmin(1)=q2nmin+xpar98
      q2pmin(2)=q2nmin+xpar99
      call ipoOm5Tables(1)
      zz=0

      if (engy.ge.10.) then
       if (engy.ge.100.) then
        if (engy.ge.1000.) then
         if (engy.ge.10000.) then
      write(ifhi,'(a,I5)')  'openhisto name DExactbn-',int(engy)
         else
      write(ifhi,'(a,I4)')  'openhisto name DExactbn-',int(engy)
        endif
        else
      write(ifhi,'(a,I3)')  'openhisto name DExactbn-',int(engy)
       endif
       else
      write(ifhi,'(a,I2)')  'openhisto name DExactbn-',int(engy)
      endif
      else
      write(ifhi,'(a,I1)')  'openhisto name DExactbn-',int(engy)
      endif
      write(ifhi,'(a)')       'htyp lbu'
      write(ifhi,'(a)')       'xmod lin ymod lin'
      write(ifhi,'(a,2e11.3)')'xrange',-bmax,bmax
c      write(ifhi,'(a)')       'yrange -.01 auto'
      write(ifhi,'(a)')       'yrange auto auto'
      write(ifhi,'(a)')      'text 0 0 "xaxis b"'
      write(ifhi,'(a,a)')    'text 0 0 "yaxis D(x+,x-,s,b)" '
      write(ifhi,'(a,e8.2,a)')  'text 0.1 0.9 "E=',engy,' GeV"'
      if (xpar8.eq.1.) then
      write(ifhi,'(a,f5.2,a)')  'text 0.1 0.8 "x=',x,'"'
      write(ifhi,'(a,f5.2,a)')  'text 0.1 0.7 "y=',y,'"'
      endif
      write(ifhi,'(a)')       'array 2'

      Dint=Dsoftshval(x,y,0.,0)
      do i=0,nptg
        b=-bmax+2.*float(i)/float(nptg)*bmax
        z=Dsoftshval(x,y,b,0)/Dint
        write(ifhi,*) b,z
      enddo

      write(ifhi,'(a)')    '  endarray'
      write(ifhi,'(a)')    'closehisto plot 0-'

c**********************************************************************

      if (engy.ge.10.) then
       if (engy.ge.100.) then
        if (engy.ge.1000.) then
         if (engy.ge.10000.) then
      write(ifhi,'(a,I5)')       'openhisto name DParambn-',int(engy)
         else
      write(ifhi,'(a,I4)')       'openhisto name DParambn-',int(engy)
        endif
        else
      write(ifhi,'(a,I3)')       'openhisto name DParambn-',int(engy)
       endif
       else
      write(ifhi,'(a,I2)')       'openhisto name DParambn-',int(engy)
      endif
      else
      write(ifhi,'(a,I1)')       'openhisto name DParambn-',int(engy)
      endif
      write(ifhi,'(a)')       'htyp lra'
      write(ifhi,'(a)')       'xmod lin ymod lin'
      write(ifhi,'(a,2e11.3)')'xrange',-bmax,bmax
c      write(ifhi,'(a)')       'yrange -.01 auto'
      write(ifhi,'(a)')       'yrange auto auto'
      write(ifhi,'(a)')      'text 0 0 "xaxis b"'
      write(ifhi,'(a,a)')    'text 0 0 "yaxis D(x+,x-,s,b)" '
      write(ifhi,'(a,e8.2,a)')  'text 0.1 0.9 "E=',engy,' GeV"'
      if (xpar8.eq.1.) then
      write(ifhi,'(a,f5.2,a)')  'text 0.1 0.8 "x=',x,'"'
      write(ifhi,'(a,f5.2,a)')  'text 0.1 0.7 "y=',y,'"'
      endif
      write(ifhi,'(a)')       'array 2'


      Dint=xDfit(zz,idxDmin(iomega),idxDmax(iomega),engy**2,
     &     sngl(dsqrt(x)*dexp(y)),sngl(dsqrt(x)*dexp(-y)),0.)

      do i=0,nptg
        b=-bmax+2.*float(i)/float(nptg)*bmax
        z=xDfit(zz,idxDmin(iomega),idxDmax(iomega),engy**2,
     &       sngl(dsqrt(x)*dexp(y)),sngl(dsqrt(x)*dexp(-y)),b)

        write(ifhi,*) b,z/Dint
      enddo

      write(ifhi,'(a)')    '  endarray'
      write(ifhi,'(a)')    'closehisto plot 0-'


c**********************************************************************

      if (engy.ge.10.) then
       if (engy.ge.100.) then
        if (engy.ge.1000.) then
         if (engy.ge.10000.) then
      write(ifhi,'(a,I5)')  'openhisto name DEfitb-',int(engy)
         else
      write(ifhi,'(a,I4)')  'openhisto name DEfitb-',int(engy)
        endif
        else
      write(ifhi,'(a,I3)')  'openhisto name DEfitb-',int(engy)
       endif
       else
      write(ifhi,'(a,I2)')  'openhisto name DEfitb-',int(engy)
      endif
      else
      write(ifhi,'(a,I1)')  'openhisto name DEfitb-',int(engy)
      endif
      write(ifhi,'(a)')       'htyp pfc'
      write(ifhi,'(a)')       'xmod lin ymod lin'
      write(ifhi,'(a,2e11.3)')'xrange',-bmax,bmax
c      write(ifhi,'(a)')       'yrange -.01 auto'
      write(ifhi,'(a)')       'yrange auto auto'
      write(ifhi,'(a)')      'text 0 0 "xaxis b"'
      write(ifhi,'(a,a)')    'text 0 0 "yaxis D(x+,x-,s,b)" '
      write(ifhi,'(a,e8.2,a)')  'text 0.1 0.9 "E=',engy,' GeV"'
      if (xpar8.eq.1.) then
      write(ifhi,'(a,f5.2,a)')  'text 0.1 0.8 "x=',x,'"'
      write(ifhi,'(a,f5.2,a)')  'text 0.1 0.7 "y=',y,'"'
      endif
      write(ifhi,'(a)')       'array 2'

      sig2=sigma2(x,0)
      if(sig2.le.0.) sig2=1.e+10
      do i=0,nptg
        b=-bmax+2.*float(i)/float(nptg)*bmax
        write(ifhi,*) b,exp(-b**2/sig2)
      enddo

      write(ifhi,'(a)')    '  endarray'
      write(ifhi,'(a)')    'closehisto plot 0-'


c**********************************************************************

      if (engy.ge.10.) then
       if (engy.ge.100.) then
        if (engy.ge.1000.) then
         if (engy.ge.10000.) then
      write(ifhi,'(a,I5)')  'openhisto name DEfitbnSoft-',int(engy)
         else
      write(ifhi,'(a,I4)')  'openhisto name DEfitbnSoft-',int(engy)
        endif
        else
      write(ifhi,'(a,I3)')  'openhisto name DEfitbnSoft-',int(engy)
       endif
       else
      write(ifhi,'(a,I2)')  'openhisto name DEfitbnSoft-',int(engy)
      endif
      else
      write(ifhi,'(a,I1)')  'openhisto name DEfitbnSoft-',int(engy)
      endif
      write(ifhi,'(a)')       'htyp pft'
      write(ifhi,'(a)')       'xmod lin ymod lin'
      write(ifhi,'(a,2e11.3)')'xrange',-bmax,bmax
c      write(ifhi,'(a)')       'yrange -.01 auto'
      write(ifhi,'(a)')       'yrange auto auto'
      write(ifhi,'(a)')      'text 0 0 "xaxis b"'
      write(ifhi,'(a,a)')    'text 0 0 "yaxis D(x+,x-,s,b)" '
      write(ifhi,'(a,e8.2,a)')  'text 0.1 0.9 "E=',engy,' GeV"'
      if (xpar8.eq.1.) then
      write(ifhi,'(a,f5.2,a)')  'text 0.1 0.8 "x=',x,'"'
      write(ifhi,'(a,f5.2,a)')  'text 0.1 0.7 "y=',y,'"'
      endif
      write(ifhi,'(a)')       'array 2'

      sig2=sigma2(x,3)
      if(sig2.le.0.) sig2=1.e+10
      do i=0,nptg
        b=-bmax+2.*float(i)/float(nptg)*bmax
        write(ifhi,*) b,exp(-b**2/sig2)
      enddo

      write(ifhi,'(a)')    '  endarray'
      write(ifhi,'(a)')    'closehisto plot 0-'


c**********************************************************************


      if (engy.ge.10.) then
       if (engy.ge.100.) then
        if (engy.ge.1000.) then
         if (engy.ge.10000.) then
      write(ifhi,'(a,I5)')  'openhisto name DEfitbnSh-',int(engy)
         else
      write(ifhi,'(a,I4)')  'openhisto name DEfitbnSh-',int(engy)
        endif
        else
      write(ifhi,'(a,I3)')  'openhisto name DEfitbnSh-',int(engy)
       endif
       else
      write(ifhi,'(a,I2)')  'openhisto name DEfitbnSh-',int(engy)
      endif
      else
      write(ifhi,'(a,I1)')  'openhisto name DEfitbnSh-',int(engy)
      endif
      write(ifhi,'(a)')       'htyp poc'
      write(ifhi,'(a)')       'xmod lin ymod lin'
      write(ifhi,'(a,2e11.3)')'xrange',-bmax,bmax
c      write(ifhi,'(a)')       'yrange -.01 auto'
      write(ifhi,'(a)')       'yrange auto auto'
      write(ifhi,'(a)')      'text 0 0 "xaxis b"'
      write(ifhi,'(a,a)')    'text 0 0 "yaxis D(x+,x-,s,b)" '
      write(ifhi,'(a,e8.2,a)')  'text 0.1 0.9 "E=',engy,' GeV"'
      if (xpar8.eq.1.) then
      write(ifhi,'(a,f5.2,a)')  'text 0.1 0.8 "x=',x,'"'
      write(ifhi,'(a,f5.2,a)')  'text 0.1 0.7 "y=',y,'"'
      endif
      write(ifhi,'(a)')       'array 2'

      sig2=sigma2(x,2)
      if(sig2.le.0.) sig2=1.e+10
      do i=0,nptg
        b=-bmax+2.*float(i)/float(nptg)*bmax
        write(ifhi,*) b,exp(-b**2/sig2)
      enddo

      write(ifhi,'(a)')    '  endarray'

      write(ifhi,'(a)')    'closehisto plot 0'

      end

c----------------------------------------------------------------------
      subroutine xbnParD
c----------------------------------------------------------------------

#include "aaa.h"
#include "par.h"

      double precision x,y,Dsoftshval,z,xDfit,om51



      nptg=50     !number of point for the graphs
      bmax=xpar2                   !value max of b (impact parameter)
      bmax=max(0.1,bmax)
      y=dble(xpar3)              !value of y (rapidity)
      x=dble(xpar4)
      q2pmin(1)=q2nmin+xpar98
      q2pmin(2)=q2nmin+xpar99
      call ipoOm5Tables(1)
      zz=0
      if(x.le.1.d-12)x=1.d0

      if (engy.ge.10.) then
       if (engy.ge.100.) then
        if (engy.ge.1000.) then
         if (engy.ge.10000.) then
      write(ifhi,'(a,I5)')  'openhisto name DExactbn-',int(engy)
         else
      write(ifhi,'(a,I4)')  'openhisto name DExactbn-',int(engy)
        endif
        else
      write(ifhi,'(a,I3)')  'openhisto name DExactbn-',int(engy)
       endif
       else
      write(ifhi,'(a,I2)')  'openhisto name DExactbn-',int(engy)
      endif
      else
      write(ifhi,'(a,I1)')  'openhisto name DExactbn-',int(engy)
      endif
      write(ifhi,'(a)')       'htyp lbu'
      write(ifhi,'(a)')       'xmod lin ymod lin'
      write(ifhi,'(a,2e11.3)')'xrange',-bmax,bmax
c      write(ifhi,'(a)')       'yrange -.01 auto'
      write(ifhi,'(a)')       'yrange auto auto'
      write(ifhi,'(a)')      'text 0 0 "xaxis b"'
      write(ifhi,'(a,a)')    'text 0 0 "yaxis D(x+,x-,s,b)" '
      write(ifhi,'(a,e8.2,a)')  'text 0.1 0.9 "E=',engy,' GeV"'
      if (xpar8.eq.1.) then
      write(ifhi,'(a,f5.2,a)')  'text 0.1 0.8 "x=',x,'"'
      write(ifhi,'(a,f5.2,a)')  'text 0.1 0.7 "y=',y,'"'
      endif
      write(ifhi,'(a)')       'array 2'

      Dint=om51(x,y,0.,0,-5)
      do i=0,nptg
        b=-bmax+2.*float(i)/float(nptg)*bmax
        z=om51(x,y,b,0,-5)/Dint
        write(ifhi,*) b,z
      enddo

      write(ifhi,'(a)')    '  endarray'
      write(ifhi,'(a)')    'closehisto plot 0-'

c**********************************************************************
      if (engy.ge.10.) then
       if (engy.ge.100.) then
        if (engy.ge.1000.) then
         if (engy.ge.10000.) then
      write(ifhi,'(a,I5)')  'openhisto name Dfitbn-',int(engy)
         else
      write(ifhi,'(a,I4)')  'openhisto name Dfitbn-',int(engy)
        endif
        else
      write(ifhi,'(a,I3)')  'openhisto name Dfitbn-',int(engy)
       endif
       else
      write(ifhi,'(a,I2)')  'openhisto name Dfitbn-',int(engy)
      endif
      else
      write(ifhi,'(a,I1)')  'openhisto name Dfitbn-',int(engy)
      endif
      write(ifhi,'(a)')       'htyp lzi'
      write(ifhi,'(a)')       'xmod lin ymod lin'
      write(ifhi,'(a,2e11.3)')'xrange',-bmax,bmax
c      write(ifhi,'(a)')       'yrange -.01 auto'
      write(ifhi,'(a)')       'yrange auto auto'
      write(ifhi,'(a)')      'text 0 0 "xaxis b"'
      write(ifhi,'(a,a)')    'text 0 0 "yaxis D(x+,x-,s,b)" '
      write(ifhi,'(a,e8.2,a)')  'text 0.1 0.9 "E=',engy,' GeV"'
      if (xpar8.eq.1.) then
      write(ifhi,'(a,f5.2,a)')  'text 0.1 0.8 "x=',x,'"'
      write(ifhi,'(a,f5.2,a)')  'text 0.1 0.7 "y=',y,'"'
      endif
      write(ifhi,'(a)')       'array 2'

      Dint=Dsoftshval(x,y,0.,0)
      do i=0,nptg
        b=-bmax+2.*float(i)/float(nptg)*bmax
        z=Dsoftshval(x,y,b,0)/Dint
        write(ifhi,*) b,z
      enddo

      write(ifhi,'(a)')    '  endarray'
      write(ifhi,'(a)')    'closehisto plot 0-'

c**********************************************************************

      if (engy.ge.10.) then
       if (engy.ge.100.) then
        if (engy.ge.1000.) then
         if (engy.ge.10000.) then
      write(ifhi,'(a,I5)')       'openhisto name DParambn-',int(engy)
         else
      write(ifhi,'(a,I4)')       'openhisto name DParambn-',int(engy)
        endif
        else
      write(ifhi,'(a,I3)')       'openhisto name DParambn-',int(engy)
       endif
       else
      write(ifhi,'(a,I2)')       'openhisto name DParambn-',int(engy)
      endif
      else
      write(ifhi,'(a,I1)')       'openhisto name DParambn-',int(engy)
      endif
      write(ifhi,'(a)')       'htyp lra'
      write(ifhi,'(a)')       'xmod lin ymod lin'
      write(ifhi,'(a,2e11.3)')'xrange',-bmax,bmax
c      write(ifhi,'(a)')       'yrange -.01 auto'
      write(ifhi,'(a)')       'yrange auto auto'
      write(ifhi,'(a)')      'text 0 0 "xaxis b"'
      write(ifhi,'(a,a)')    'text 0 0 "yaxis D(x+,x-,s,b)" '
      write(ifhi,'(a,e8.2,a)')  'text 0.1 0.9 "E=',engy,' GeV"'
      if (xpar8.eq.1.) then
      write(ifhi,'(a,f5.2,a)')  'text 0.1 0.8 "x=',x,'"'
      write(ifhi,'(a,f5.2,a)')  'text 0.1 0.7 "y=',y,'"'
      endif
      write(ifhi,'(a)')       'array 2'

      Dint=xDfit(zz,idxDmin(iomega),idxDmax(iomega),smaxDf,
     &     sngl(dsqrt(x)*dexp(y)),sngl(dsqrt(x)*dexp(-y)),0.)


      do i=0,nptg
        b=-bmax+2.*float(i)/float(nptg)*bmax
        z=xDfit(zz,idxDmin(iomega),idxDmax(iomega),smaxDf,
     &       sngl(dsqrt(x)*dexp(y)),sngl(dsqrt(x)*dexp(-y)),b)

        write(ifhi,*) b,z/Dint
      enddo

      write(ifhi,'(a)')    '  endarray'
      write(ifhi,'(a)')    'closehisto plot 0-'


c**********************************************************************

      if (engy.ge.10.) then
       if (engy.ge.100.) then
        if (engy.ge.1000.) then
         if (engy.ge.10000.) then
      write(ifhi,'(a,I5)')  'openhisto name DEfitb-',int(engy)
         else
      write(ifhi,'(a,I4)')  'openhisto name DEfitb-',int(engy)
        endif
        else
      write(ifhi,'(a,I3)')  'openhisto name DEfitb-',int(engy)
       endif
       else
      write(ifhi,'(a,I2)')  'openhisto name DEfitb-',int(engy)
      endif
      else
      write(ifhi,'(a,I1)')  'openhisto name DEfitb-',int(engy)
      endif
      write(ifhi,'(a)')       'htyp pfc'
      write(ifhi,'(a)')       'xmod lin ymod lin'
      write(ifhi,'(a,2e11.3)')'xrange',-bmax,bmax
c      write(ifhi,'(a)')       'yrange -.01 auto'
      write(ifhi,'(a)')       'yrange auto auto'
      write(ifhi,'(a)')      'text 0 0 "xaxis b"'
      write(ifhi,'(a,a)')    'text 0 0 "yaxis D(x+,x-,s,b)" '
      write(ifhi,'(a,e8.2,a)')  'text 0.1 0.9 "E=',engy,' GeV"'
      if (xpar8.eq.1.) then
      write(ifhi,'(a,f5.2,a)')  'text 0.1 0.8 "x=',x,'"'
      write(ifhi,'(a,f5.2,a)')  'text 0.1 0.7 "y=',y,'"'
      endif
      write(ifhi,'(a)')       'array 2'

      sig2=sigma2(x,0)
      if(sig2.le.0.) sig2=1.e+10
      do i=0,nptg
        b=-bmax+2.*float(i)/float(nptg)*bmax
        write(ifhi,*) b,exp(-b**2/sig2)
      enddo

      write(ifhi,'(a)')    '  endarray'
      write(ifhi,'(a)')    'closehisto plot 0-'


c**********************************************************************

      if (engy.ge.10.) then
       if (engy.ge.100.) then
        if (engy.ge.1000.) then
         if (engy.ge.10000.) then
      write(ifhi,'(a,I5)')  'openhisto name DEintb-',int(engy)
         else
      write(ifhi,'(a,I4)')  'openhisto name DEintb-',int(engy)
        endif
        else
      write(ifhi,'(a,I3)')  'openhisto name DEintb-',int(engy)
       endif
       else
      write(ifhi,'(a,I2)')  'openhisto name DEintb-',int(engy)
      endif
      else
      write(ifhi,'(a,I1)')  'openhisto name DEintb-',int(engy)
      endif
      write(ifhi,'(a)')       'htyp pft'
      write(ifhi,'(a)')       'xmod lin ymod lin'
      write(ifhi,'(a,2e11.3)')'xrange',-bmax,bmax
c      write(ifhi,'(a)')       'yrange -.01 auto'
      write(ifhi,'(a)')       'yrange auto auto'
      write(ifhi,'(a)')      'text 0 0 "xaxis b"'
      write(ifhi,'(a,a)')    'text 0 0 "yaxis D(x+,x-,s,b)" '
      write(ifhi,'(a,e8.2,a)')  'text 0.1 0.9 "E=',engy,' GeV"'
      if (xpar8.eq.1.) then
      write(ifhi,'(a,f5.2,a)')  'text 0.1 0.8 "x=',x,'"'
      write(ifhi,'(a,f5.2,a)')  'text 0.1 0.7 "y=',y,'"'
      endif
      write(ifhi,'(a)')       'array 2'

      sig2=sigma1i(x)
      do i=0,nptg
        b=-bmax+2.*float(i)/float(nptg)*bmax
        write(ifhi,*) b,exp(-b**2/sig2)
      enddo

      write(ifhi,'(a)')    '  endarray'
      write(ifhi,'(a)')    'closehisto plot 0-'


c**********************************************************************


      if (engy.ge.10.) then
       if (engy.ge.100.) then
        if (engy.ge.1000.) then
         if (engy.ge.10000.) then
      write(ifhi,'(a,I5)')  'openhisto name DPfitb-',int(engy)
         else
      write(ifhi,'(a,I4)')  'openhisto name DPfitb-',int(engy)
        endif
        else
      write(ifhi,'(a,I3)')  'openhisto name DPfitb-',int(engy)
       endif
       else
      write(ifhi,'(a,I2)')  'openhisto name DPfitb-',int(engy)
      endif
      else
      write(ifhi,'(a,I1)')  'openhisto name DPfitb-',int(engy)
      endif
      write(ifhi,'(a)')       'htyp poc'
      write(ifhi,'(a)')       'xmod lin ymod lin'
      write(ifhi,'(a,2e11.3)')'xrange',-bmax,bmax
c      write(ifhi,'(a)')       'yrange -.01 auto'
      write(ifhi,'(a)')       'yrange auto auto'
      write(ifhi,'(a)')      'text 0 0 "xaxis b"'
      write(ifhi,'(a,a)')    'text 0 0 "yaxis D(x+,x-,s,b)" '
      write(ifhi,'(a,e8.2,a)')  'text 0.1 0.9 "E=',engy,' GeV"'
      if (xpar8.eq.1.) then
      write(ifhi,'(a,f5.2,a)')  'text 0.1 0.8 "x=',x,'"'
      write(ifhi,'(a,f5.2,a)')  'text 0.1 0.7 "y=',y,'"'
      endif
      write(ifhi,'(a)')       'array 2'

      sig2=xsigmafit(x)
      do i=0,nptg
        b=-bmax+2.*float(i)/float(nptg)*bmax
        write(ifhi,*) b,exp(-b**2/sig2)
      enddo

      write(ifhi,'(a)')    '  endarray'
      write(ifhi,'(a)')    'closehisto plot 0'



      end


c----------------------------------------------------------------------
      subroutine xbParD
c----------------------------------------------------------------------

#include "aaa.h"
#include "sem.h"
#include "par.h"

      double precision x,om51,xDfit,z(0:200),y
c     & ,omYuncut,omNpuncut,omYcut

      nptg=50                  !number of point for the graphs
      x=dble(xpar4)                 !value of biniDf (impact parameter)
      y=dble(xpar3)                 !value of y (rapidity)
      if(x.le.1.d-12)x=1.d0
      if(x.gt.dexp(-2.d0*y))x=x*dexp(-2.d0*y)
      zz=0
      q20=q2nmin
      q2pmin(1)=q20+xpar98
      q2pmin(2)=q20+xpar99
      call ipoOm5Tables(1) 

      bmax=xpar2
      bmax=max(0.1,bmax)
      t=1.
c      iqqN=0
c      iqq=int(xpar7)
      iqq1=0
      iqq2=4       !no diffraction

      if (engy.ge.10.) then
       if (engy.ge.100.) then
        if (engy.ge.1000.) then
         if (engy.ge.10000.) then
      write(ifhi,'(a,I5)')       'openhisto name DbParam-',int(engy)
         else
      write(ifhi,'(a,I4)')       'openhisto name DbParam-',int(engy)
        endif
        else
      write(ifhi,'(a,I3)')       'openhisto name DbParam-',int(engy)
       endif
       else
      write(ifhi,'(a,I2)')       'openhisto name DbParam-',int(engy)
      endif
      else
      write(ifhi,'(a,I1)')       'openhisto name DbParam-',int(engy)
      endif
      write(ifhi,'(a)')       'htyp lbu'
      if(xpar5.eq.0.)then
        write(ifhi,'(a)')       'xmod lin ymod lin'
      else
        write(ifhi,'(a)')       'xmod lin ymod log'
      endif
      write(ifhi,'(a,2e11.3)')'xrange',0.,bmax
      write(ifhi,'(a)')       'yrange auto auto'
      write(ifhi,'(a)')      'text 0 0 "xaxis b"'
      write(ifhi,'(a,a)')    'text 0 0 "yaxis D(x+,x-,s,b)" '
      write(ifhi,'(a,e8.2,a)')  'text 0.1 0.9 "E=',engy,' GeV"'
      if (xpar8.eq.1.) then
      write(ifhi,'(a,f5.2,a)')  'text 0.1 0.8 "x=',x,'"'
      write(ifhi,'(a,f5.2,a)')  'text 0.1 0.7 "y=',y,'"'
      endif
      write(ifhi,'(a)')       'array 2'

      xp=sqrt(sngl(x))*exp(sngl(y))
      xm=sqrt(sngl(x))*exp(-sngl(y))

      if(xpar6.eq.1.)then
       t=sngl(xDfit(zz,idxDmin(iomega),idxDmax(iomega),smaxDf,xp,xm,0.))
      endif
      if(abs(t).lt.1.d-8)t=1.
      do i=0,nptg
        b=float(i)/float(nptg)*bmax
        z(i)=xDfit(zz,idxDmin(iomega),idxDmax(iomega),smaxDf,xp,xm,b)
     .      /dble(t)
        write(ifhi,*) b,z(i)
       enddo

      write(ifhi,'(a)')    '  endarray'

      write(ifhi,'(a)')    'closehisto plot 0-'


c**********************************************************************

      if (engy.ge.10.) then
        if (engy.ge.100.) then
          if (engy.ge.1000.) then
            if (engy.ge.10000.) then
           write(ifhi,'(a,I5)')  'openhisto name DbParamI-',int(engy)
            else
           write(ifhi,'(a,I4)')  'openhisto name DbParamI-',int(engy)
            endif
          else
           write(ifhi,'(a,I3)')  'openhisto name DbParamI-',int(engy)
          endif
        else
          write(ifhi,'(a,I2)')  'openhisto name DbParamI-',int(engy)
        endif
      else
        write(ifhi,'(a,I1)')  'openhisto name DbParamI-',int(engy)
      endif
      write(ifhi,'(a)')       'htyp lga'
      write(ifhi,'(a)')       'xmod lin ymod lin'
      write(ifhi,'(a,2e11.3)')'xrange',0.,bmax
      write(ifhi,'(a)')       'yrange .01 auto'
      write(ifhi,'(a)')      'text 0 0 "xaxis b"'
      write(ifhi,'(a,a)')    'text 0 0 "yaxis D(x+,x-,s,b)" '
      write(ifhi,'(a,e8.2,a)')  'text 0.1 0.9 "E=',engy,' GeV"'
      if (xpar8.eq.1.) then
        write(ifhi,'(a,f5.2,a)')  'text 0.1 0.8 "x=',x,'"'
        write(ifhi,'(a,f5.2,a)')  'text 0.1 0.7 "y=',y,'"'
      endif
      write(ifhi,'(a)')       'array 2'

      if(xpar6.eq.1.)then
        t=(alpD(idxD0,iclpro,icltar)
     &   *(smaxDf*sngl(x))**betD(idxD0,iclpro,icltar)
     &                +alpD(1,iclpro,icltar)
     &   *(smaxDf*sngl(x))**betD(1,iclpro,icltar))
     &   /float(idxD0+1)
        if(idxD.gt.1)t=t+alpD(idxD,iclpro,icltar)
     &   *(smaxDf*sngl(x))**betD(idxD,iclpro,icltar)
      endif
      if(abs(t).lt.1.d-8)t=1.
      do i=0,nptg
        b=float(i)/float(nptg)*bmax
        tmp=(alpD(idxD0,iclpro,icltar)/t
     &   *(smaxDf*sngl(x))**(betD(idxD0,iclpro,icltar)
     &           +gamD(idxD0,iclpro,icltar)*b**2.)
     &        *exp(-b**2./delD(idxD0,iclpro,icltar))
     &                 +alpD(1,iclpro,icltar)/t
     &   *(smaxDf*sngl(x))**(betD(1,iclpro,icltar)
     &           +gamD(1,iclpro,icltar)*b**2.)
     &        *exp(-b**2./delD(1,iclpro,icltar)))
     &   /float(idxD0+1)
        if(idxD.gt.1)tmp=tmp+alpD(idxD,iclpro,icltar)/t
     &   *(smaxDf*sngl(x))**(betD(idxD,iclpro,icltar)
     &           +gamD(idxD,iclpro,icltar)*b**2.)
     &        *exp(-b**2./delD(idxD,iclpro,icltar))
        write(ifhi,*) b,tmp
      enddo

      write(ifhi,'(a)')    '  endarray'
      write(ifhi,'(a)')    'closehisto plot 0-'

c**********************************************************************

      if (engy.ge.10.) then
        if (engy.ge.100.) then
          if (engy.ge.1000.) then
            if (engy.ge.10000.) then
           write(ifhi,'(a,I5)')  'openhisto name DbExaI-',int(engy)
            else
           write(ifhi,'(a,I4)')  'openhisto name DbExaI-',int(engy)
            endif
          else
           write(ifhi,'(a,I3)')  'openhisto name DbExaI-',int(engy)
          endif
        else
          write(ifhi,'(a,I2)')  'openhisto name DbExaI-',int(engy)
        endif
      else
        write(ifhi,'(a,I1)')  'openhisto name DbExaI-',int(engy)
      endif
      write(ifhi,'(a)')       'htyp lro'
      write(ifhi,'(a)')       'xmod lin ymod lin'
      write(ifhi,'(a,2e11.3)')'xrange',0.,bmax
      write(ifhi,'(a)')       'yrange .01 auto'
      write(ifhi,'(a)')      'text 0 0 "xaxis b"'
      write(ifhi,'(a,a)')    'text 0 0 "yaxis D(x+,x-,s,b)" '
      write(ifhi,'(a,e8.2,a)')  'text 0.1 0.9 "E=',engy,' GeV"'
      if (xpar8.eq.1.) then
        write(ifhi,'(a,f5.2,a)')  'text 0.1 0.8 "x=',x,'"'
        write(ifhi,'(a,f5.2,a)')  'text 0.1 0.7 "y=',y,'"'
      endif
      write(ifhi,'(a)')       'array 2'

      if(xpar6.eq.1.)then
        t=sngl(om51(x,0.d0,0.,iqq1,iqq2)
     &         /(x**dble(-alppar)))
      endif
      if(abs(t).lt.1.d-8)t=1.
      nptg=nptg/2
      do i=0,nptg
        b=float(i)/float(nptg)*bmax
        tmp=sngl(om51(x,y,b,iqq1,iqq2)
     &         /(x**dble(-alppar)))/t
        write(ifhi,*) b,tmp
      enddo

      write(ifhi,'(a)')    '  endarray'

      write(ifhi,'(a)')    'closehisto plot 0'

      end

c----------------------------------------------------------------------
      subroutine xGexaJ
c----------------------------------------------------------------------

#include "aaa.h"
#include "ems.h"
#include "sem.h"
#include "par.h"

      parameter(nptg=50) !number of point for the graphs
      double precision x(0:nptg),omJ1(0:nptg),xtmp,w(0:nptg)
      double precision z(0:nptg),om5J,om1fJ!,omYuncut,omYuncutJ
      double precision xminr,y,t,omJ2(0:nptg),omJ3(0:nptg),omJ4(0:nptg)
     &,omJ5(0:nptg),wndi,wndic,om1f,xp,xm!,omJ6(0:nptg),omJ7(0:nptg)


      kollsave=koll        !koll modified in zzfz
      koll=1
      biniDf=xpar2                 !value of biniDf (impact parameter)
      y=dble(xpar3)
      t=1.d0
      xtmp=xmaxDf
      xmaxDf=exp(-2.d0*y)
      zzp=0.
      zzt=0.
c define variables for om5J
      do i=idxDmin(iomega),idxDmax(iomega)
       zp=zzp
       zt=zzt
       call Gfunpar(zp,zt,0.,0.,1,i,biniDf,smaxDf
     &             ,alpx,betx,betpx,epsp,epst,epss,gamv)
      enddo

      xminr=xminDf  !value of xminr for plotting the function

      nnnmax=1
      if(iscreen.ne.0)nnnmax=2

      do nnn=1,nnnmax

        if(nnn.eq.1)then
          q20=q2nmin
          q2pmin(1)=q20+xpar98
          q2pmin(2)=q20+xpar99
        else
          q20=q2nmin
          q2pmin(1)=q20
          q2pmin(2)=q20
        endif
        call ipoOm5Tables(1) 

      if (engy.ge.10.) then
       if (engy.ge.100.) then
        if (engy.ge.1000.) then
         if (engy.ge.10000.) then
      write(ifhi,'(a,I5)')  'openhisto name GIexa-',int(engy)
         else
      write(ifhi,'(a,I4)')  'openhisto name GIexa-',int(engy)
        endif
        else
      write(ifhi,'(a,I3)')  'openhisto name GIexa-',int(engy)
       endif
       else
      write(ifhi,'(a,I2)')  'openhisto name GIexa-',int(engy)
      endif
      else
      write(ifhi,'(a,I1)')  'openhisto name GIexa-',int(engy)
      endif
      write(ifhi,'(a)')       'htyp lru'
      if(xpar5.eq.0.)then
        write(ifhi,'(a)')       'xmod log ymod log'
      else
        write(ifhi,'(a)')       'xmod log ymod lin'
      endif
      write(ifhi,'(a,2e11.3)')'xrange',xminr,xmaxDf
      write(ifhi,'(a)')       'yrange 1e-5  auto'
      write(ifhi,'(a)')      'text 0 0 "xaxis x"'
      write(ifhi,'(a,a)')    'text 0 0 "yaxis D(x+,s,b)" '
      write(ifhi,'(a,e8.2,a)')  'text 0.1 0.9 "E=',engy,' GeV"'
      if (xpar8.eq.1.) then
      write(ifhi,'(a,f5.2,a)')  'text 0.1 0.8 "b=',biniDf,' fm"'
      endif
      write(ifhi,'(a)')       'array 2'

      x(0)=xminr
      do i=0,nptg
        if (i.ne.0) x(i)=x(0)*(xmaxDf/xminr)**(dble(i)/dble(nptg))
        z(i)=0.d0
        xp=dble(sqrt(x(i))*exp(y))
        xm=dble(x(i))/xp
        if(nnn.eq.1)then
          wndic=om5J(xp,xm,biniDf,0,4)
          omJ1(i)=om5J(xp,xm,biniDf,0,0)
          omJ2(i)=om5J(xp,xm,biniDf,1,1)
          omJ3(i)=om5J(xp,xm,biniDf,2,3)
          omJ4(i)=om5J(xp,xm,biniDf,4,4)
          w(i)=om5J(xp,xm,biniDf,-10,-10)
          omJ5(i)=max(0.d0,w(i)-wndic)
        else
          omJ1(i)=om5J(xp,xm,biniDf,0,0)
          omJ2(i)=om5J(xp,xm,biniDf,1,1)
          omJ3(i)=om5J(xp,xm,biniDf,2,3)
          omJ4(i)=om5J(xp,xm,biniDf,4,4)
          omJ5(i)=om5J(xp,xm,biniDf,5,9)
          w(i)=om5J(xp,xm,biniDf,0,4)+om5J(xp,xm,biniDf,-5,-5)
        endif
        z(i)=omJ1(i)+omJ2(i)+omJ3(i)+omJ4(i)+omJ5(i)
        if(xpar6.eq.2)t=w(i)
        if(xpar6.eq.1)t=z(i)
        if(xpar6.eq.0)t=1d0
        if(t.gt.0d0)write(ifhi,*) x(i),z(i)/t
      enddo

      write(ifhi,'(a)')    '  endarray'
      write(ifhi,'(a)')    'closehisto plot 0-'


c*************************************************************************

      if (engy.ge.10.) then
       if (engy.ge.100.) then
        if (engy.ge.1000.) then
         if (engy.ge.10000.) then
      write(ifhi,'(a,I5)')       'openhisto name GIsoft-',int(engy)
         else
      write(ifhi,'(a,I4)')       'openhisto name GIsoft-',int(engy)
        endif
        else
      write(ifhi,'(a,I3)')       'openhisto name GIsoft-',int(engy)
       endif
       else
      write(ifhi,'(a,I2)')       'openhisto name GIsoft-',int(engy)
      endif
      else
      write(ifhi,'(a,I1)')       'openhisto name GIsoft-',int(engy)
      endif
      write(ifhi,'(a)')       'htyp lba'
      write(ifhi,'(a)')       'xmod log ymod log'
      write(ifhi,'(a,2e11.3)')'xrange',xminr,xmaxDf
      write(ifhi,'(a)')       'yrange 1e-5  auto'
      write(ifhi,'(a)')      'text 0 0 "xaxis x"'
      write(ifhi,'(a,a)')    'text 0 0 "yaxis D(x+,s,b)" '
      write(ifhi,'(a,e8.2,a)')  'text 0.1 0.9 "E=',engy,' GeV"'
      if (xpar8.eq.1.) then
      write(ifhi,'(a,f5.2,a)')  'text 0.1 0.8 "b=',biniDf,' fm"'
      endif
      write(ifhi,'(a)')       'array 2'

      do i=0,nptg
        if(xpar6.eq.2)t=w(i)
        if(xpar6.eq.1)t=z(i)
        if(xpar6.eq.0)t=1d0
        if(t.gt.0d0)write(ifhi,*) x(i),omJ1(i)/t
      enddo

      write(ifhi,'(a)')    '  endarray'
      write(ifhi,'(a)')    'closehisto plot 0-'

c*************************************************************************

      if (engy.ge.10.) then
       if (engy.ge.100.) then
        if (engy.ge.1000.) then
         if (engy.ge.10000.) then
      write(ifhi,'(a,I5)')       'openhisto name GIgg-',int(engy)
         else
      write(ifhi,'(a,I4)')       'openhisto name GIgg-',int(engy)
        endif
        else
      write(ifhi,'(a,I3)')       'openhisto name GIgg-',int(engy)
       endif
       else
      write(ifhi,'(a,I2)')       'openhisto name GIgg-',int(engy)
      endif
      else
      write(ifhi,'(a,I1)')       'openhisto name GIgg-',int(engy)
      endif
      write(ifhi,'(a)')       'htyp lza'
      write(ifhi,'(a)')       'xmod log ymod log'
      write(ifhi,'(a,2e11.3)')'xrange',xminr,xmaxDf
      write(ifhi,'(a)')       'yrange 1e-5  auto'
      write(ifhi,'(a)')      'text 0 0 "xaxis x"'
      write(ifhi,'(a,a)')    'text 0 0 "yaxis D(x+,s,b)" '
      write(ifhi,'(a,e8.2,a)')  'text 0.1 0.9 "E=',engy,' GeV"'
      if (xpar8.eq.1.) then
      write(ifhi,'(a,f5.2,a)')  'text 0.1 0.8 "b=',biniDf,' fm"'
      endif
      write(ifhi,'(a)')       'array 2'

      do i=0,nptg
        if(xpar6.eq.2)t=w(i)
        if(xpar6.eq.1)t=z(i)
        if(xpar6.eq.0)t=1d0
        if(t.gt.0d0)write(ifhi,*) x(i),omJ2(i)/t
      enddo

      write(ifhi,'(a)')    '  endarray'
      write(ifhi,'(a)')    'closehisto plot 0-'

c*************************************************************************

      if (engy.ge.10.) then
       if (engy.ge.100.) then
        if (engy.ge.1000.) then
         if (engy.ge.10000.) then
      write(ifhi,'(a,I5)')       'openhisto name GIgq-',int(engy)
         else
      write(ifhi,'(a,I4)')       'openhisto name GIgq-',int(engy)
        endif
        else
      write(ifhi,'(a,I3)')       'openhisto name GIgq-',int(engy)
       endif
       else
      write(ifhi,'(a,I2)')       'openhisto name GIgq-',int(engy)
      endif
      else
      write(ifhi,'(a,I1)')       'openhisto name GIgq-',int(engy)
      endif
      write(ifhi,'(a)')       'htyp lga'
      write(ifhi,'(a)')       'xmod log ymod log'
      write(ifhi,'(a,2e11.3)')'xrange',xminr,xmaxDf
      write(ifhi,'(a)')       'yrange 1e-5  auto'
      write(ifhi,'(a)')      'text 0 0 "xaxis x"'
      write(ifhi,'(a,a)')    'text 0 0 "yaxis D(x+,s,b)" '
      write(ifhi,'(a,e8.2,a)')  'text 0.1 0.9 "E=',engy,' GeV"'
      if (xpar8.eq.1.) then
      write(ifhi,'(a,f5.2,a)')  'text 0.1 0.8 "b=',biniDf,' fm"'
      endif
      write(ifhi,'(a)')       'array 2'

      do i=0,nptg
        if(xpar6.eq.2)t=w(i)
        if(xpar6.eq.1)t=z(i)
        if(xpar6.eq.0)t=1d0
        if(t.gt.0d0)write(ifhi,*) x(i),omJ3(i)/t
      enddo

      write(ifhi,'(a)')    '  endarray'
      write(ifhi,'(a)')    'closehisto plot 0-'

c*************************************************************************

      if (engy.ge.10.) then
       if (engy.ge.100.) then
        if (engy.ge.1000.) then
         if (engy.ge.10000.) then
      write(ifhi,'(a,I5)')       'openhisto name GIqq-',int(engy)
         else
      write(ifhi,'(a,I4)')       'openhisto name GIqq-',int(engy)
        endif
        else
      write(ifhi,'(a,I3)')       'openhisto name GIqq-',int(engy)
       endif
       else
      write(ifhi,'(a,I2)')       'openhisto name GIqq-',int(engy)
      endif
      else
      write(ifhi,'(a,I1)')       'openhisto name GIqq-',int(engy)
      endif
      write(ifhi,'(a)')       'htyp lya'
      write(ifhi,'(a)')       'xmod log ymod log'
      write(ifhi,'(a,2e11.3)')'xrange',xminr,xmaxDf
      write(ifhi,'(a)')       'yrange 1e-5  auto'
      write(ifhi,'(a)')      'text 0 0 "xaxis x"'
      write(ifhi,'(a,a)')    'text 0 0 "yaxis D(x+,s,b)" '
      write(ifhi,'(a,e8.2,a)')  'text 0.1 0.9 "E=',engy,' GeV"'
      if (xpar8.eq.1.) then
      write(ifhi,'(a,f5.2,a)')  'text 0.1 0.8 "b=',biniDf,' fm"'
      endif
      write(ifhi,'(a)')       'array 2'

      do i=0,nptg
        if(xpar6.eq.2)t=w(i)
        if(xpar6.eq.1)t=z(i)
        if(xpar6.eq.0)t=1d0
        if(t.gt.0d0)write(ifhi,*) x(i),omJ4(i)/t
      enddo

      write(ifhi,'(a)')    '  endarray'
      write(ifhi,'(a)')    'closehisto plot 0-'

c*************************************************************************

      if (engy.ge.10.) then
       if (engy.ge.100.) then
        if (engy.ge.1000.) then
         if (engy.ge.10000.) then
      write(ifhi,'(a,I5)')       'openhisto name GIsh-',int(engy)
         else
      write(ifhi,'(a,I4)')       'openhisto name GIsh-',int(engy)
        endif
        else
      write(ifhi,'(a,I3)')       'openhisto name GIsh-',int(engy)
       endif
       else
      write(ifhi,'(a,I2)')       'openhisto name GIsh-',int(engy)
      endif
      else
      write(ifhi,'(a,I1)')       'openhisto name GIsh-',int(engy)
      endif
      write(ifhi,'(a)')       'htyp lra'
      write(ifhi,'(a)')       'xmod log ymod log'
      write(ifhi,'(a,2e11.3)')'xrange',xminr,xmaxDf
      write(ifhi,'(a)')       'yrange 1e-5  auto'
      write(ifhi,'(a)')      'text 0 0 "xaxis x"'
      write(ifhi,'(a,a)')    'text 0 0 "yaxis D(x+,s,b)" '
      write(ifhi,'(a,e8.2,a)')  'text 0.1 0.9 "E=',engy,' GeV"'
      if (xpar8.eq.1.) then
      write(ifhi,'(a,f5.2,a)')  'text 0.1 0.8 "b=',biniDf,' fm"'
      endif
      write(ifhi,'(a)')       'array 2'

      do i=0,nptg
        if(xpar6.eq.2)t=w(i)
        if(xpar6.eq.1)t=z(i)
        if(xpar6.eq.0)t=1d0
        if(t.gt.0d0)write(ifhi,*) x(i),(omJ2(i)+omJ3(i)+omJ4(i))/t
      enddo

      write(ifhi,'(a)')    '  endarray'
      write(ifhi,'(a)')    'closehisto plot 0-'

c*************************************************************************

      if (engy.ge.10.) then
       if (engy.ge.100.) then
        if (engy.ge.1000.) then
         if (engy.ge.10000.) then
      write(ifhi,'(a,I5)')       'openhisto name GIdif-',int(engy)
         else
      write(ifhi,'(a,I4)')       'openhisto name GIdif-',int(engy)
        endif
        else
      write(ifhi,'(a,I3)')       'openhisto name GIdif-',int(engy)
       endif
       else
      write(ifhi,'(a,I2)')       'openhisto name GIdif-',int(engy)
      endif
      else
      write(ifhi,'(a,I1)')       'openhisto name GIdif-',int(engy)
      endif
      write(ifhi,'(a)')       'htyp lfa'
      write(ifhi,'(a)')       'xmod log ymod log'
      write(ifhi,'(a,2e11.3)')'xrange',xminr,xmaxDf
      write(ifhi,'(a)')       'yrange 1e-5  auto'
      write(ifhi,'(a)')      'text 0 0 "xaxis x"'
      write(ifhi,'(a,a)')    'text 0 0 "yaxis D(x+,s,b)" '
      write(ifhi,'(a,e8.2,a)')  'text 0.1 0.9 "E=',engy,' GeV"'
      if (xpar8.eq.1.) then
      write(ifhi,'(a,f5.2,a)')  'text 0.1 0.8 "b=',biniDf,' fm"'
      endif
      write(ifhi,'(a)')       'array 2'

      do i=0,nptg
        if(xpar6.eq.2)t=w(i)
        if(xpar6.eq.1)t=z(i)
        if(xpar6.eq.0)t=1d0
        if(t.gt.0d0)write(ifhi,*) x(i),omJ5(i)/t
      enddo

      write(ifhi,'(a)')    '  endarray'
      write(ifhi,'(a)')    'closehisto plot 0-'

c*************************************************************************

      if (engy.ge.10.) then
       if (engy.ge.100.) then
        if (engy.ge.1000.) then
         if (engy.ge.10000.) then
      write(ifhi,'(a,I5)')       'openhisto name GItot-',int(engy)
         else
      write(ifhi,'(a,I4)')       'openhisto name GItot-',int(engy)
        endif
        else
      write(ifhi,'(a,I3)')       'openhisto name GItot-',int(engy)
       endif
       else
      write(ifhi,'(a,I2)')       'openhisto name GItot-',int(engy)
      endif
      else
      write(ifhi,'(a,I1)')       'openhisto name GItot-',int(engy)
      endif
      write(ifhi,'(a)')       'htyp pfc'
      write(ifhi,'(a)')       'xmod log ymod log'
      write(ifhi,'(a,2e11.3)')'xrange',xminr,xmaxDf
      write(ifhi,'(a)')       'yrange 1e-5  auto'
      write(ifhi,'(a)')      'text 0 0 "xaxis x"'
      write(ifhi,'(a,a)')    'text 0 0 "yaxis D(x+,s,b)" '
      write(ifhi,'(a,e8.2,a)')  'text 0.1 0.9 "E=',engy,' GeV"'
      if (xpar8.eq.1.) then
      write(ifhi,'(a,f5.2,a)')  'text 0.1 0.8 "b=',biniDf,' fm"'
      endif
      write(ifhi,'(a)')       'array 2'

      do i=0,nptg
        if(xpar6.eq.2)t=w(i)
        if(xpar6.eq.1)t=z(i)
        if(xpar6.eq.0)t=1d0
        if(t.gt.0d0)write(ifhi,*) x(i),(omJ1(i)+omJ2(i)+omJ3(i)+omJ4(i)
     &                  +omJ5(i))/t
      enddo

      write(ifhi,'(a)')    '  endarray'

      write(ifhi,'(a)')    'closehisto plot 0-'

      if (engy.ge.10.) then
       if (engy.ge.100.) then
        if (engy.ge.1000.) then
         if (engy.ge.10000.) then
      write(ifhi,'(a,I5)')  'openhisto name GIfit-',int(engy)
         else
      write(ifhi,'(a,I4)')  'openhisto name GIfit-',int(engy)
        endif
        else
      write(ifhi,'(a,I3)')  'openhisto name GIfit-',int(engy)
       endif
       else
      write(ifhi,'(a,I2)')  'openhisto name GIfit-',int(engy)
      endif
      else
      write(ifhi,'(a,I1)')  'openhisto name GIfit-',int(engy)
      endif
      write(ifhi,'(a)')       'htyp lro'
      if(xpar5.eq.0.)then
        write(ifhi,'(a)')       'xmod log ymod log'
      else
        write(ifhi,'(a)')       'xmod log ymod lin'
      endif
      write(ifhi,'(a,2e11.3)')'xrange',xminr,xmaxDf
      write(ifhi,'(a)')       'yrange 1e-5 auto'
      write(ifhi,'(a)')      'text 0 0 "xaxis x"'
      write(ifhi,'(a,a)')    'text 0 0 "yaxis D(x+,s,b)" '
      write(ifhi,'(a,e8.2,a)')  'text 0.1 0.9 "E=',engy,' GeV"'
      if (xpar8.eq.1.) then
      write(ifhi,'(a,f5.2,a)')  'text 0.1 0.8 "b=',biniDf,' fm"'
      endif
      write(ifhi,'(a)')       'array 2'

      do i=idxDmin(iomega),idxDmax(iomega)
        zp=zzp
        zt=zzt
        call Gfunpar(zp,zt,0.,0.,1,i,biniDf,smaxDf
     &              ,alpx,betx,betpx,epsp,epst,epss,gamv)
      enddo
      x(0)=xminr
      do i=0,nptg
        if (i.ne.0) x(i)=x(0)*(xmaxDf/xminr)**(dble(i)/dble(nptg))
        z(i)=0.d0
        xp=dble(sqrt(x(i))*exp(y))
        xm=dble(x(i))/xp
        if(nnn.eq.1)then
          wndi=om1f(xp,xm,0,idxD)
          wndic=om1fJ(xp,xm,0,idxD)
          omJ1(i)=om1fJ(xp,xm,0,0)/wndic*wndi
          omJ2(i)=om1fJ(xp,xm,1,idxD)/wndic*wndi
          omJ3(i)=om1fJ(xp,xm,idxD1,idxD1)
        else
          omJ1(i)=om1f(xp,xm,0,0)
          omJ2(i)=om1f(xp,xm,1,idxD)
          omJ3(i)=om1f(xp,xm,idxD1,idxD1)
          w(i)=omJ1(i)+omJ2(i)+omJ3(i)
        endif  
        z(i)=omJ1(i)+omJ2(i)+omJ3(i)
        if(xpar6.eq.2)t=w(i)
        if(xpar6.eq.1)t=z(i)
        if(xpar6.eq.0)t=1d0
        if(t.gt.0d0)write(ifhi,*) x(i),z(i)/t
      enddo

      write(ifhi,'(a)')    '  endarray'
      write(ifhi,'(a)')    'closehisto plot 0-'


c*************************************************************************

      if (engy.ge.10.) then
       if (engy.ge.100.) then
        if (engy.ge.1000.) then
         if (engy.ge.10000.) then
      write(ifhi,'(a,I5)')       'openhisto name GIfsoft-',int(engy)
         else
      write(ifhi,'(a,I4)')       'openhisto name GIfsoft-',int(engy)
        endif
        else
      write(ifhi,'(a,I3)')       'openhisto name GIfsoft-',int(engy)
       endif
       else
      write(ifhi,'(a,I2)')       'openhisto name GIfsoft-',int(engy)
      endif
      else
      write(ifhi,'(a,I1)')       'openhisto name GIfsoft-',int(engy)
      endif
      write(ifhi,'(a)')       'htyp lbo'
      write(ifhi,'(a)')       'xmod log ymod log'
      write(ifhi,'(a,2e11.3)')'xrange',xminr,xmaxDf
      write(ifhi,'(a)')       'yrange 1e-5  auto'
      write(ifhi,'(a)')      'text 0 0 "xaxis x"'
      write(ifhi,'(a,a)')    'text 0 0 "yaxis D(x+,s,b)" '
      write(ifhi,'(a,e8.2,a)')  'text 0.1 0.9 "E=',engy,' GeV"'
      if (xpar8.eq.1.) then
      write(ifhi,'(a,f5.2,a)')  'text 0.1 0.8 "b=',biniDf,' fm"'
      endif
      write(ifhi,'(a)')       'array 2'

      do i=0,nptg
        if(xpar6.eq.2)t=w(i)
        if(xpar6.eq.1)t=z(i)
        if(xpar6.eq.0)t=1d0
        if(t.gt.0d0)write(ifhi,*) x(i),omJ1(i)/t
      enddo

      write(ifhi,'(a)')    '  endarray'
      write(ifhi,'(a)')    'closehisto plot 0-'

c*************************************************************************

      if (engy.ge.10.) then
       if (engy.ge.100.) then
        if (engy.ge.1000.) then
         if (engy.ge.10000.) then
      write(ifhi,'(a,I5)')       'openhisto name GIfsh-',int(engy)
         else
      write(ifhi,'(a,I4)')       'openhisto name GIfsh-',int(engy)
        endif
        else
      write(ifhi,'(a,I3)')       'openhisto name GIfsh-',int(engy)
       endif
       else
      write(ifhi,'(a,I2)')       'openhisto name GIfsh-',int(engy)
      endif
      else
      write(ifhi,'(a,I1)')       'openhisto name GIfsh-',int(engy)
      endif
      write(ifhi,'(a)')       'htyp lro'
      write(ifhi,'(a)')       'xmod log ymod log'
      write(ifhi,'(a,2e11.3)')'xrange',xminr,xmaxDf
      write(ifhi,'(a)')       'yrange 1e-5  auto'
      write(ifhi,'(a)')      'text 0 0 "xaxis x"'
      write(ifhi,'(a,a)')    'text 0 0 "yaxis D(x+,s,b)" '
      write(ifhi,'(a,e8.2,a)')  'text 0.1 0.9 "E=',engy,' GeV"'
      if (xpar8.eq.1.) then
      write(ifhi,'(a,f5.2,a)')  'text 0.1 0.8 "b=',biniDf,' fm"'
      endif
      write(ifhi,'(a)')       'array 2'

      do i=0,nptg
        if(xpar6.eq.2)t=w(i)
        if(xpar6.eq.1)t=z(i)
        if(xpar6.eq.0)t=1d0
        if(t.gt.0d0)write(ifhi,*) x(i),omJ2(i)/t
      enddo

      write(ifhi,'(a)')    '  endarray'
      write(ifhi,'(a)')    'closehisto plot 0-'

c*************************************************************************

      if (engy.ge.10.) then
       if (engy.ge.100.) then
        if (engy.ge.1000.) then
         if (engy.ge.10000.) then
      write(ifhi,'(a,I5)')       'openhisto name GIfdif-',int(engy)
         else
      write(ifhi,'(a,I4)')       'openhisto name GIfdif-',int(engy)
        endif
        else
      write(ifhi,'(a,I3)')       'openhisto name GIfdif-',int(engy)
       endif
       else
      write(ifhi,'(a,I2)')       'openhisto name GIfdif-',int(engy)
      endif
      else
      write(ifhi,'(a,I1)')       'openhisto name GIfdif-',int(engy)
      endif
      write(ifhi,'(a)')       'htyp lfo'
      write(ifhi,'(a)')       'xmod log ymod log'
      write(ifhi,'(a,2e11.3)')'xrange',xminr,xmaxDf
      write(ifhi,'(a)')       'yrange 1e-5  auto'
      write(ifhi,'(a)')      'text 0 0 "xaxis x"'
      write(ifhi,'(a,a)')    'text 0 0 "yaxis D(x+,s,b)" '
      write(ifhi,'(a,e8.2,a)')  'text 0.1 0.9 "E=',engy,' GeV"'
      if (xpar8.eq.1.) then
      write(ifhi,'(a,f5.2,a)')  'text 0.1 0.8 "b=',biniDf,' fm"'
      endif
      write(ifhi,'(a)')       'array 2'

      do i=0,nptg
        if(xpar6.eq.2)t=w(i)
        if(xpar6.eq.1)t=z(i)
        if(xpar6.eq.0)t=1d0
        if(t.gt.0d0)write(ifhi,*) x(i),omJ3(i)/t
      enddo

      write(ifhi,'(a)')    '  endarray'
      write(ifhi,'(a)')    'closehisto plot 0'


      enddo

      xmaxDf=xtmp
      koll=kollsave

      end

c----------------------------------------------------------------------
      subroutine xGexaDiff
c----------------------------------------------------------------------

#include "aaa.h"
#include "ems.h"
#include "sem.h"
#include "par.h"

      parameter(nptg=50) !number of point for the graphs
      double precision x(0:nptg),omJ1(0:nptg),w(0:nptg)
      double precision z(0:nptg),om5J
      double precision xminr,y,t,omJ2(0:nptg),omJ3(0:nptg),omJ4(0:nptg)
     &,omJ5(0:nptg),xp,xm

      if(iomega.ne.2)then

      kollsave=koll        !koll modified in zzfz
      koll=1
      biniDf=xpar2                 !value of biniDf (impact parameter)
      y=dble(xpar3)
      t=1.d0
      zzp=0.
      zzt=0.
      q20=q2nmin
      q2pmin(1)=q20+xpar98
      q2pmin(2)=q20+xpar99
      call ipoOm5Tables(1) 
c define variables for om5J
      do i=idxD0,idxD1
       zp=zzp
       zt=zzt
       call Gfunpar(zp,zt,0.,0.,1,i,biniDf,smaxDf
     &             ,alpx,betx,betpx,epsp,epst,epss,gamv)
      enddo

      xminr=xminDf  !value of xminr for plotting the function

      if (engy.ge.10.) then
       if (engy.ge.100.) then
        if (engy.ge.1000.) then
         if (engy.ge.10000.) then
      write(ifhi,'(a,I5)')  'openhisto name GDexa-',int(engy)
         else
      write(ifhi,'(a,I4)')  'openhisto name GDexa-',int(engy)
        endif
        else
      write(ifhi,'(a,I3)')  'openhisto name GDexa-',int(engy)
       endif
       else
      write(ifhi,'(a,I2)')  'openhisto name GDexa-',int(engy)
      endif
      else
      write(ifhi,'(a,I1)')  'openhisto name GDexa-',int(engy)
      endif
      write(ifhi,'(a)')       'htyp lru'
      if(xpar5.eq.0.)then
        write(ifhi,'(a)')       'xmod log ymod log'
      else
        write(ifhi,'(a)')       'xmod log ymod lin'
      endif
      write(ifhi,'(a,2e11.3)')'xrange',xminr,xmaxDf
      write(ifhi,'(a)')       'yrange 1e-5  auto'
      write(ifhi,'(a)')      'text 0 0 "xaxis x"'
      write(ifhi,'(a,a)')    'text 0 0 "yaxis D(x+,s,b)" '
      write(ifhi,'(a,e8.2,a)')  'text 0.1 0.9 "E=',engy,' GeV"'
      if (xpar8.eq.1.) then
      write(ifhi,'(a,f5.2,a)')  'text 0.1 0.8 "b=',biniDf,' fm"'
      endif
      write(ifhi,'(a)')       'array 2'

      x(0)=xminr
      do i=0,nptg
        if (i.ne.0) x(i)=x(0)*(xmaxDf/xminr)**(dble(i)/dble(nptg))
        omJ1(i)=0.d0
        omJ2(i)=0.d0
        omJ3(i)=0.d0
        omJ4(i)=0.d0
        omJ5(i)=0.d0
        z(i)=0.d0
        w(i)=0.d0
        xp=dble(sqrt(x(i))*exp(y))
        xm=dble(x(i))/xp
        if(xp.le.1.d0.and.xm.le.1d0)then
          omJ1(i)=om5J(xp,xm,biniDf,6,6)
          omJ2(i)=om5J(xp,xm,biniDf,7,7)
          omJ3(i)=om5J(xp,xm,biniDf,8,8)
          omJ4(i)=om5J(xp,xm,biniDf,9,9)
          omJ5(i)=om5J(xp,xm,biniDf,0,5)
          w(i)=om5J(xp,xm,biniDf,-5,-5)
          z(i)=omJ1(i)+omJ2(i)+omJ3(i)+omJ4(i)
        elseif(xpar6.ne.0)then
          w(i)=1d0
          z(i)=1d0
        endif
        if(xpar6.eq.3)t=z(i)/w(i)
        if(xpar6.eq.2)t=w(i)
        if(xpar6.eq.1)t=z(i)
        if(xpar6.eq.0)t=1d0
        if(t.gt.0d0)write(ifhi,*) x(i),z(i)/t
      enddo

      write(ifhi,'(a)')    '  endarray'
      write(ifhi,'(a)')    'closehisto plot 0-'


c*************************************************************************

      if (engy.ge.10.) then
       if (engy.ge.100.) then
        if (engy.ge.1000.) then
         if (engy.ge.10000.) then
      write(ifhi,'(a,I5)')       'openhisto name GDSDt-',int(engy)
         else
      write(ifhi,'(a,I4)')       'openhisto name GDSDt-',int(engy)
        endif
        else
      write(ifhi,'(a,I3)')       'openhisto name GDSDt-',int(engy)
       endif
       else
      write(ifhi,'(a,I2)')       'openhisto name GDSDt-',int(engy)
      endif
      else
      write(ifhi,'(a,I1)')       'openhisto name GDSDt-',int(engy)
      endif
      write(ifhi,'(a)')       'htyp lba'
      write(ifhi,'(a)')       'xmod log ymod log'
      write(ifhi,'(a,2e11.3)')'xrange',xminr,xmaxDf
      write(ifhi,'(a)')       'yrange 1e-5  auto'
      write(ifhi,'(a)')      'text 0 0 "xaxis x"'
      write(ifhi,'(a,a)')    'text 0 0 "yaxis D(x+,s,b)" '
      write(ifhi,'(a,e8.2,a)')  'text 0.1 0.9 "E=',engy,' GeV"'
      if (xpar8.eq.1.) then
      write(ifhi,'(a,f5.2,a)')  'text 0.1 0.8 "b=',biniDf,' fm"'
      endif
      write(ifhi,'(a)')       'array 2'

      do i=0,nptg
        if(xpar6.eq.3)t=z(i)/w(i)
        if(xpar6.eq.2)t=w(i)
        if(xpar6.eq.1)t=z(i)
        if(xpar6.eq.0)t=1d0
        if(t.gt.0d0)write(ifhi,*) x(i),omJ1(i)/t
      enddo

      write(ifhi,'(a)')    '  endarray'
      write(ifhi,'(a)')    'closehisto plot 0-'

c*************************************************************************

      if (engy.ge.10.) then
       if (engy.ge.100.) then
        if (engy.ge.1000.) then
         if (engy.ge.10000.) then
      write(ifhi,'(a,I5)')       'openhisto name GDSDp-',int(engy)
         else
      write(ifhi,'(a,I4)')       'openhisto name GDSDp-',int(engy)
        endif
        else
      write(ifhi,'(a,I3)')       'openhisto name GDSDp-',int(engy)
       endif
       else
      write(ifhi,'(a,I2)')       'openhisto name GDSD+-',int(engy)
      endif
      else
      write(ifhi,'(a,I1)')       'openhisto name GDSD+-',int(engy)
      endif
      write(ifhi,'(a)')       'htyp lza'
      write(ifhi,'(a)')       'xmod log ymod log'
      write(ifhi,'(a,2e11.3)')'xrange',xminr,xmaxDf
      write(ifhi,'(a)')       'yrange 1e-5  auto'
      write(ifhi,'(a)')      'text 0 0 "xaxis x"'
      write(ifhi,'(a,a)')    'text 0 0 "yaxis D(x+,s,b)" '
      write(ifhi,'(a,e8.2,a)')  'text 0.1 0.9 "E=',engy,' GeV"'
      if (xpar8.eq.1.) then
      write(ifhi,'(a,f5.2,a)')  'text 0.1 0.8 "b=',biniDf,' fm"'
      endif
      write(ifhi,'(a)')       'array 2'

      do i=0,nptg
        if(xpar6.eq.3)t=z(i)/w(i)
        if(xpar6.eq.2)t=w(i)
        if(xpar6.eq.1)t=z(i)
        if(xpar6.eq.0)t=1d0
        if(t.gt.0d0)write(ifhi,*) x(i),omJ2(i)/t
      enddo

      write(ifhi,'(a)')    '  endarray'
      write(ifhi,'(a)')    'closehisto plot 0-'

c*************************************************************************

      if (engy.ge.10.) then
       if (engy.ge.100.) then
        if (engy.ge.1000.) then
         if (engy.ge.10000.) then
      write(ifhi,'(a,I5)')       'openhisto name GDDD-',int(engy)
         else
      write(ifhi,'(a,I4)')       'openhisto name GDDD-',int(engy)
        endif
        else
      write(ifhi,'(a,I3)')       'openhisto name GDDD-',int(engy)
       endif
       else
      write(ifhi,'(a,I2)')       'openhisto name GDDD-',int(engy)
      endif
      else
      write(ifhi,'(a,I1)')       'openhisto name GDDD-',int(engy)
      endif
      write(ifhi,'(a)')       'htyp lga'
      write(ifhi,'(a)')       'xmod log ymod log'
      write(ifhi,'(a,2e11.3)')'xrange',xminr,xmaxDf
      write(ifhi,'(a)')       'yrange 1e-5  auto'
      write(ifhi,'(a)')      'text 0 0 "xaxis x"'
      write(ifhi,'(a,a)')    'text 0 0 "yaxis D(x+,s,b)" '
      write(ifhi,'(a,e8.2,a)')  'text 0.1 0.9 "E=',engy,' GeV"'
      if (xpar8.eq.1.) then
      write(ifhi,'(a,f5.2,a)')  'text 0.1 0.8 "b=',biniDf,' fm"'
      endif
      write(ifhi,'(a)')       'array 2'

      do i=0,nptg
        if(xpar6.eq.3)t=z(i)/w(i)
        if(xpar6.eq.2)t=w(i)
        if(xpar6.eq.1)t=z(i)
        if(xpar6.eq.0)t=1d0
        if(t.gt.0d0)write(ifhi,*) x(i),omJ3(i)/t
      enddo

      write(ifhi,'(a)')    '  endarray'
      write(ifhi,'(a)')    'closehisto plot 0-'

c*************************************************************************

      if (engy.ge.10.) then
       if (engy.ge.100.) then
        if (engy.ge.1000.) then
         if (engy.ge.10000.) then
      write(ifhi,'(a,I5)')       'openhisto name GDCD-',int(engy)
         else
      write(ifhi,'(a,I4)')       'openhisto name GDCD-',int(engy)
        endif
        else
      write(ifhi,'(a,I3)')       'openhisto name GDCD-',int(engy)
       endif
       else
      write(ifhi,'(a,I2)')       'openhisto name GDCD-',int(engy)
      endif
      else
      write(ifhi,'(a,I1)')       'openhisto name GDCD-',int(engy)
      endif
      write(ifhi,'(a)')       'htyp lyo'
      write(ifhi,'(a)')       'xmod log ymod log'
      write(ifhi,'(a,2e11.3)')'xrange',xminr,xmaxDf
      write(ifhi,'(a)')       'yrange 1e-5  auto'
      write(ifhi,'(a)')      'text 0 0 "xaxis x"'
      write(ifhi,'(a,a)')    'text 0 0 "yaxis D(x+,s,b)" '
      write(ifhi,'(a,e8.2,a)')  'text 0.1 0.9 "E=',engy,' GeV"'
      if (xpar8.eq.1.) then
      write(ifhi,'(a,f5.2,a)')  'text 0.1 0.8 "b=',biniDf,' fm"'
      endif
      write(ifhi,'(a)')       'array 2'

      do i=0,nptg
        if(xpar6.eq.3)t=z(i)/w(i)
        if(xpar6.eq.2)t=w(i)
        if(xpar6.eq.1)t=z(i)
        if(xpar6.eq.0)t=1d0
        if(t.gt.0d0)write(ifhi,*) x(i),omJ4(i)/t
      enddo

      write(ifhi,'(a)')    '  endarray'
      write(ifhi,'(a)')    'closehisto plot 0-'

c*************************************************************************

      if (engy.ge.10.) then
       if (engy.ge.100.) then
        if (engy.ge.1000.) then
         if (engy.ge.10000.) then
      write(ifhi,'(a,I5)')       'openhisto name GDSD-',int(engy)
         else
      write(ifhi,'(a,I4)')       'openhisto name GDSD-',int(engy)
        endif
        else
      write(ifhi,'(a,I3)')       'openhisto name GDSD-',int(engy)
       endif
       else
      write(ifhi,'(a,I2)')       'openhisto name GDSD-',int(engy)
      endif
      else
      write(ifhi,'(a,I1)')       'openhisto name GDSD-',int(engy)
      endif
      write(ifhi,'(a)')       'htyp lra'
      write(ifhi,'(a)')       'xmod log ymod log'
      write(ifhi,'(a,2e11.3)')'xrange',xminr,xmaxDf
      write(ifhi,'(a)')       'yrange 1e-5  auto'
      write(ifhi,'(a)')      'text 0 0 "xaxis x"'
      write(ifhi,'(a,a)')    'text 0 0 "yaxis D(x+,s,b)" '
      write(ifhi,'(a,e8.2,a)')  'text 0.1 0.9 "E=',engy,' GeV"'
      if (xpar8.eq.1.) then
      write(ifhi,'(a,f5.2,a)')  'text 0.1 0.8 "b=',biniDf,' fm"'
      endif
      write(ifhi,'(a)')       'array 2'

      do i=0,nptg
        if(xpar6.eq.3)t=z(i)/w(i)
        if(xpar6.eq.2)t=w(i)
        if(xpar6.eq.1)t=z(i)
        if(xpar6.eq.0)t=1d0
        if(t.gt.0d0)write(ifhi,*) x(i),(omJ1(i)+omJ2(i))/t
      enddo

      write(ifhi,'(a)')    '  endarray'
      write(ifhi,'(a)')    'closehisto plot 0-'

c*************************************************************************

      if (engy.ge.10.) then
       if (engy.ge.100.) then
        if (engy.ge.1000.) then
         if (engy.ge.10000.) then
      write(ifhi,'(a,I5)')       'openhisto name GDfit-',int(engy)
         else
      write(ifhi,'(a,I4)')       'openhisto name GDfit-',int(engy)
        endif
        else
      write(ifhi,'(a,I3)')       'openhisto name GDfit-',int(engy)
       endif
       else
      write(ifhi,'(a,I2)')       'openhisto name GDfit-',int(engy)
      endif
      else
      write(ifhi,'(a,I1)')       'openhisto name GDfit-',int(engy)
      endif
      write(ifhi,'(a)')       'htyp lfa'
      write(ifhi,'(a)')       'xmod log ymod log'
      write(ifhi,'(a,2e11.3)')'xrange',xminr,xmaxDf
      write(ifhi,'(a)')       'yrange 1e-5  auto'
      write(ifhi,'(a)')      'text 0 0 "xaxis x"'
      write(ifhi,'(a,a)')    'text 0 0 "yaxis D(x+,s,b)" '
      write(ifhi,'(a,e8.2,a)')  'text 0.1 0.9 "E=',engy,' GeV"'
      if (xpar8.eq.1.) then
      write(ifhi,'(a,f5.2,a)')  'text 0.1 0.8 "b=',biniDf,' fm"'
      endif
      write(ifhi,'(a)')       'array 2'

      do i=0,nptg
        t=1d0
        if(xpar6.eq.2)t=w(i)
        if(xpar6.eq.1)t=z(i)
        if(xpar6.eq.0)t=1d0
        if(t.gt.0d0)write(ifhi,*) x(i),w(i)/t
      enddo

      write(ifhi,'(a)')    '  endarray'
      write(ifhi,'(a)')    'closehisto plot 0-'

c*************************************************************************

      if (engy.ge.10.) then
       if (engy.ge.100.) then
        if (engy.ge.1000.) then
         if (engy.ge.10000.) then
      write(ifhi,'(a,I5)')       'openhisto name GNDtot-',int(engy)
         else
      write(ifhi,'(a,I4)')       'openhisto name GNDtot-',int(engy)
        endif
        else
      write(ifhi,'(a,I3)')       'openhisto name GNDtot-',int(engy)
       endif
       else
      write(ifhi,'(a,I2)')       'openhisto name GNDtot-',int(engy)
      endif
      else
      write(ifhi,'(a,I1)')       'openhisto name GNDtot-',int(engy)
      endif
      write(ifhi,'(a)')       'htyp pfc'
      write(ifhi,'(a)')       'xmod log ymod log'
      write(ifhi,'(a,2e11.3)')'xrange',xminr,xmaxDf
      write(ifhi,'(a)')       'yrange 1e-5  auto'
      write(ifhi,'(a)')      'text 0 0 "xaxis x"'
      write(ifhi,'(a,a)')    'text 0 0 "yaxis D(x+,s,b)" '
      write(ifhi,'(a,e8.2,a)')  'text 0.1 0.9 "E=',engy,' GeV"'
      if (xpar8.eq.1.) then
      write(ifhi,'(a,f5.2,a)')  'text 0.1 0.8 "b=',biniDf,' fm"'
      endif
      write(ifhi,'(a)')       'array 2'

      do i=0,nptg
        t=1d0
        if(xpar6.eq.2)t=w(i)
        if(xpar6.eq.1)t=z(i)
        if(xpar6.eq.0)t=1d0
        if(t.gt.0d0)write(ifhi,*) x(i),omJ5(i)/t
      enddo

      write(ifhi,'(a)')    '  endarray'

      write(ifhi,'(a)')    'closehisto plot 0'


      koll=kollsave

      endif
      end

c----------------------------------------------------------------------
      subroutine xyGexaDiff
c----------------------------------------------------------------------

#include "aaa.h"
#include "ems.h"
#include "sem.h"
#include "par.h"

      parameter(nptg=20) !number of point for the graphs
      double precision y(0:nptg),omJ1(0:nptg),w(0:nptg)
      double precision z(0:nptg),om5J
      double precision ymax,x,t,omJ2(0:nptg),omJ3(0:nptg),omJ4(0:nptg)
     &,omJ5(0:nptg),xp,xm

      if(iomega.ne.2)then

      kollsave=koll        !koll modified in zzfz
      koll=1
      biniDf=xpar2                 !value of biniDf (impact parameter)
      x=dble(xpar4)
      t=1.d0
      zzp=0.
      zzt=0.
c define variables for om5J
      do i=idxD0,idxD1
       zp=zzp
       zt=zzt
       call Gfunpar(zp,zt,0.,0.,1,i,biniDf,smaxDf
     &             ,alpx,betx,betpx,epsp,epst,epss,gamv)
      enddo
      if(x.le.1.d-20)x=min(1d0,50.d0*xminDf)
      ymax=-.5d0*dlog(x)

      if (engy.ge.10.) then
       if (engy.ge.100.) then
        if (engy.ge.1000.) then
         if (engy.ge.10000.) then
      write(ifhi,'(a,I5)')  'openhisto name GDexa-',int(engy)
         else
      write(ifhi,'(a,I4)')  'openhisto name GDexa-',int(engy)
        endif
        else
      write(ifhi,'(a,I3)')  'openhisto name GDexa-',int(engy)
       endif
       else
      write(ifhi,'(a,I2)')  'openhisto name GDexa-',int(engy)
      endif
      else
      write(ifhi,'(a,I1)')  'openhisto name GDexa-',int(engy)
      endif
      write(ifhi,'(a)')       'htyp lru'
      if(xpar5.eq.0.)then
        write(ifhi,'(a)')       'xmod lin ymod log'
      else
        write(ifhi,'(a)')       'xmod lin ymod lin'
      endif
      write(ifhi,'(a,2e11.3)')'xrange',-ymax-0.1,ymax+0.1
      write(ifhi,'(a)')       'yrange 1e-5  auto'
      write(ifhi,'(a)')      'text 0 0 "xaxis y"'
      write(ifhi,'(a,a)')    'text 0 0 "yaxis D(x+,s,b)" '
      write(ifhi,'(a,e8.2,a)')  'text 0.1 0.9 "E=',engy,' GeV"'
      if (xpar8.eq.1.) then
      write(ifhi,'(a,f5.2,a)')  'text 0.1 0.8 "b=',biniDf,' fm"'
      endif
      write(ifhi,'(a)')       'array 2'

      do i=0,nptg
        y(i)=-ymax+2.*ymax*(dble(i)/dble(nptg))
        omJ1(i)=0.d0
        omJ2(i)=0.d0
        omJ3(i)=0.d0
        omJ4(i)=0.d0
        omJ5(i)=0.d0
        z(i)=0.d0
        w(i)=0.d0
        xp=dble(sqrt(x)*exp(y(i)))
        xm=dble(x)/xp
        if(xp.le.1.00001d0.and.xm.le.1.00001d0)then
          omJ1(i)=om5J(xp,xm,biniDf,6,6)
          omJ2(i)=om5J(xp,xm,biniDf,7,7)
          omJ3(i)=om5J(xp,xm,biniDf,9,9)
          omJ4(i)=om5J(xp,xm,biniDf,8,8)
          omJ5(i)=om5J(xp,xm,biniDf,0,5)
          w(i)=om5J(xp,xm,biniDf,-5,-5)
          z(i)=omJ1(i)+omJ2(i)+omJ3(i)+omJ4(i)
        elseif(xpar6.ne.0)then
          w(i)=1d0
          z(i)=1d0
        endif
        if(xpar6.eq.3)t=z(i)/w(i)
        if(xpar6.eq.2)t=w(i)
        if(xpar6.eq.1)t=z(i)
        if(xpar6.eq.0)t=1d0
        write(ifhi,*) y(i),z(i)/t
      enddo

      write(ifhi,'(a)')    '  endarray'
      write(ifhi,'(a)')    'closehisto plot 0-'


c*************************************************************************

      if (engy.ge.10.) then
       if (engy.ge.100.) then
        if (engy.ge.1000.) then
         if (engy.ge.10000.) then
      write(ifhi,'(a,I5)')       'openhisto name GDSDt-',int(engy)
         else
      write(ifhi,'(a,I4)')       'openhisto name GDSDt-',int(engy)
        endif
        else
      write(ifhi,'(a,I3)')       'openhisto name GDSDt-',int(engy)
       endif
       else
      write(ifhi,'(a,I2)')       'openhisto name GDSDt-',int(engy)
      endif
      else
      write(ifhi,'(a,I1)')       'openhisto name GDSDt-',int(engy)
      endif
      write(ifhi,'(a)')       'htyp lba'
      write(ifhi,'(a)')       'xmod lin ymod log'
      write(ifhi,'(a,2e11.3)')'xrange',-ymax-0.1,ymax+0.1
      write(ifhi,'(a)')       'yrange 1e-5  auto'
      write(ifhi,'(a)')      'text 0 0 "xaxis y"'
      write(ifhi,'(a,a)')    'text 0 0 "yaxis D(x+,s,b)" '
      write(ifhi,'(a,e8.2,a)')  'text 0.1 0.9 "E=',engy,' GeV"'
      if (xpar8.eq.1.) then
      write(ifhi,'(a,f5.2,a)')  'text 0.1 0.8 "b=',biniDf,' fm"'
      endif
      write(ifhi,'(a)')       'array 2'

      do i=0,nptg
        if(xpar6.eq.3)t=z(i)/w(i)
        if(xpar6.eq.2)t=w(i)
        if(xpar6.eq.1)t=z(i)
        if(xpar6.eq.0)t=1d0
        write(ifhi,*) y(i),omJ1(i)/t
      enddo

      write(ifhi,'(a)')    '  endarray'
      write(ifhi,'(a)')    'closehisto plot 0-'

c*************************************************************************

      if (engy.ge.10.) then
       if (engy.ge.100.) then
        if (engy.ge.1000.) then
         if (engy.ge.10000.) then
      write(ifhi,'(a,I5)')       'openhisto name GDSDp-',int(engy)
         else
      write(ifhi,'(a,I4)')       'openhisto name GDSDp-',int(engy)
        endif
        else
      write(ifhi,'(a,I3)')       'openhisto name GDSDp-',int(engy)
       endif
       else
      write(ifhi,'(a,I2)')       'openhisto name GDSD+-',int(engy)
      endif
      else
      write(ifhi,'(a,I1)')       'openhisto name GDSD+-',int(engy)
      endif
      write(ifhi,'(a)')       'htyp lza'
      write(ifhi,'(a)')       'xmod lin ymod log'
      write(ifhi,'(a,2e11.3)')'xrange',-ymax-0.1,ymax+0.1
      write(ifhi,'(a)')       'yrange 1e-5  auto'
      write(ifhi,'(a)')      'text 0 0 "xaxis y"'
      write(ifhi,'(a,a)')    'text 0 0 "yaxis D(x+,s,b)" '
      write(ifhi,'(a,e8.2,a)')  'text 0.1 0.9 "E=',engy,' GeV"'
      if (xpar8.eq.1.) then
      write(ifhi,'(a,f5.2,a)')  'text 0.1 0.8 "b=',biniDf,' fm"'
      endif
      write(ifhi,'(a)')       'array 2'

      do i=0,nptg
        if(xpar6.eq.3)t=z(i)/w(i)
        if(xpar6.eq.2)t=w(i)
        if(xpar6.eq.1)t=z(i)
        if(xpar6.eq.0)t=1d0
        write(ifhi,*) y(i),omJ2(i)/t
      enddo

      write(ifhi,'(a)')    '  endarray'
      write(ifhi,'(a)')    'closehisto plot 0-'

c*************************************************************************

      if (engy.ge.10.) then
       if (engy.ge.100.) then
        if (engy.ge.1000.) then
         if (engy.ge.10000.) then
      write(ifhi,'(a,I5)')       'openhisto name GDDD-',int(engy)
         else
      write(ifhi,'(a,I4)')       'openhisto name GDDD-',int(engy)
        endif
        else
      write(ifhi,'(a,I3)')       'openhisto name GDDD-',int(engy)
       endif
       else
      write(ifhi,'(a,I2)')       'openhisto name GDDD-',int(engy)
      endif
      else
      write(ifhi,'(a,I1)')       'openhisto name GDDD-',int(engy)
      endif
      write(ifhi,'(a)')       'htyp lga'
      write(ifhi,'(a)')       'xmod lin ymod log'
      write(ifhi,'(a,2e11.3)')'xrange',-ymax-0.1,ymax+0.1
      write(ifhi,'(a)')       'yrange 1e-5  auto'
      write(ifhi,'(a)')      'text 0 0 "xaxis y"'
      write(ifhi,'(a,a)')    'text 0 0 "yaxis D(x+,s,b)" '
      write(ifhi,'(a,e8.2,a)')  'text 0.1 0.9 "E=',engy,' GeV"'
      if (xpar8.eq.1.) then
      write(ifhi,'(a,f5.2,a)')  'text 0.1 0.8 "b=',biniDf,' fm"'
      endif
      write(ifhi,'(a)')       'array 2'

      do i=0,nptg
        if(xpar6.eq.3)t=z(i)/w(i)
        if(xpar6.eq.2)t=w(i)
        if(xpar6.eq.1)t=z(i)
        if(xpar6.eq.0)t=1d0
        write(ifhi,*) y(i),omJ3(i)/t
      enddo

      write(ifhi,'(a)')    '  endarray'
      write(ifhi,'(a)')    'closehisto plot 0-'

c*************************************************************************

      if (engy.ge.10.) then
       if (engy.ge.100.) then
        if (engy.ge.1000.) then
         if (engy.ge.10000.) then
      write(ifhi,'(a,I5)')       'openhisto name GDCD-',int(engy)
         else
      write(ifhi,'(a,I4)')       'openhisto name GDCD-',int(engy)
        endif
        else
      write(ifhi,'(a,I3)')       'openhisto name GDCD-',int(engy)
       endif
       else
      write(ifhi,'(a,I2)')       'openhisto name GDCD-',int(engy)
      endif
      else
      write(ifhi,'(a,I1)')       'openhisto name GDCD-',int(engy)
      endif
      write(ifhi,'(a)')       'htyp lyo'
      write(ifhi,'(a)')       'xmod lin ymod log'
      write(ifhi,'(a,2e11.3)')'xrange',-ymax-0.1,ymax+0.1
      write(ifhi,'(a)')       'yrange 1e-5  auto'
      write(ifhi,'(a)')      'text 0 0 "xaxis y"'
      write(ifhi,'(a,a)')    'text 0 0 "yaxis D(x+,s,b)" '
      write(ifhi,'(a,e8.2,a)')  'text 0.1 0.9 "E=',engy,' GeV"'
      if (xpar8.eq.1.) then
      write(ifhi,'(a,f5.2,a)')  'text 0.1 0.8 "b=',biniDf,' fm"'
      endif
      write(ifhi,'(a)')       'array 2'

      do i=0,nptg
        if(xpar6.eq.3)t=z(i)/w(i)
        if(xpar6.eq.2)t=w(i)
        if(xpar6.eq.1)t=z(i)
        if(xpar6.eq.0)t=1d0
        write(ifhi,*) y(i),omJ4(i)/t
      enddo

      write(ifhi,'(a)')    '  endarray'
      write(ifhi,'(a)')    'closehisto plot 0-'

c*************************************************************************

      if (engy.ge.10.) then
       if (engy.ge.100.) then
        if (engy.ge.1000.) then
         if (engy.ge.10000.) then
      write(ifhi,'(a,I5)')       'openhisto name GDSD-',int(engy)
         else
      write(ifhi,'(a,I4)')       'openhisto name GDSD-',int(engy)
        endif
        else
      write(ifhi,'(a,I3)')       'openhisto name GDSD-',int(engy)
       endif
       else
      write(ifhi,'(a,I2)')       'openhisto name GDSD-',int(engy)
      endif
      else
      write(ifhi,'(a,I1)')       'openhisto name GDSD-',int(engy)
      endif
      write(ifhi,'(a)')       'htyp lra'
      write(ifhi,'(a)')       'xmod lin ymod log'
      write(ifhi,'(a,2e11.3)')'xrange',-ymax-0.1,ymax+0.1
      write(ifhi,'(a)')       'yrange 1e-5  auto'
      write(ifhi,'(a)')      'text 0 0 "xaxis y"'
      write(ifhi,'(a,a)')    'text 0 0 "yaxis D(x+,s,b)" '
      write(ifhi,'(a,e8.2,a)')  'text 0.1 0.9 "E=',engy,' GeV"'
      if (xpar8.eq.1.) then
      write(ifhi,'(a,f5.2,a)')  'text 0.1 0.8 "b=',biniDf,' fm"'
      endif
      write(ifhi,'(a)')       'array 2'

      do i=0,nptg
        if(xpar6.eq.3)t=z(i)/w(i)
        if(xpar6.eq.2)t=w(i)
        if(xpar6.eq.1)t=z(i)
        if(xpar6.eq.0)t=1d0
        write(ifhi,*) y(i),(omJ1(i)+omJ2(i))/t
      enddo

      write(ifhi,'(a)')    '  endarray'
      write(ifhi,'(a)')    'closehisto plot 0-'

c*************************************************************************

      if (engy.ge.10.) then
       if (engy.ge.100.) then
        if (engy.ge.1000.) then
         if (engy.ge.10000.) then
      write(ifhi,'(a,I5)')       'openhisto name GDfit-',int(engy)
         else
      write(ifhi,'(a,I4)')       'openhisto name GDfit-',int(engy)
        endif
        else
      write(ifhi,'(a,I3)')       'openhisto name GDfit-',int(engy)
       endif
       else
      write(ifhi,'(a,I2)')       'openhisto name GDfit-',int(engy)
      endif
      else
      write(ifhi,'(a,I1)')       'openhisto name GDfit-',int(engy)
      endif
      write(ifhi,'(a)')       'htyp lfa'
      write(ifhi,'(a)')       'xmod lin ymod log'
      write(ifhi,'(a,2e11.3)')'xrange',-ymax-0.1,ymax+0.1
      write(ifhi,'(a)')       'yrange 1e-5  auto'
      write(ifhi,'(a)')      'text 0 0 "xaxis y"'
      write(ifhi,'(a,a)')    'text 0 0 "yaxis D(x+,s,b)" '
      write(ifhi,'(a,e8.2,a)')  'text 0.1 0.9 "E=',engy,' GeV"'
      if (xpar8.eq.1.) then
      write(ifhi,'(a,f5.2,a)')  'text 0.1 0.8 "b=',biniDf,' fm"'
      endif
      write(ifhi,'(a)')       'array 2'

      do i=0,nptg
        t=1d0
        if(xpar6.eq.2)t=w(i)
        if(xpar6.eq.1)t=z(i)
        if(xpar6.eq.0)t=1d0
        write(ifhi,*) y(i),w(i)/t
      enddo

      write(ifhi,'(a)')    '  endarray'
      write(ifhi,'(a)')    'closehisto plot 0-'

c*************************************************************************

      if (engy.ge.10.) then
       if (engy.ge.100.) then
        if (engy.ge.1000.) then
         if (engy.ge.10000.) then
      write(ifhi,'(a,I5)')       'openhisto name GNDtot-',int(engy)
         else
      write(ifhi,'(a,I4)')       'openhisto name GNDtot-',int(engy)
        endif
        else
      write(ifhi,'(a,I3)')       'openhisto name GNDtot-',int(engy)
       endif
       else
      write(ifhi,'(a,I2)')       'openhisto name GNDtot-',int(engy)
      endif
      else
      write(ifhi,'(a,I1)')       'openhisto name GNDtot-',int(engy)
      endif
      write(ifhi,'(a)')       'htyp pfc'
      write(ifhi,'(a)')       'xmod lin ymod log'
      write(ifhi,'(a,2e11.3)')'xrange',-ymax-0.1,ymax+0.1
      write(ifhi,'(a)')       'yrange 1e-5  auto'
      write(ifhi,'(a)')      'text 0 0 "xaxis y"'
      write(ifhi,'(a,a)')    'text 0 0 "yaxis D(x+,s,b)" '
      write(ifhi,'(a,e8.2,a)')  'text 0.1 0.9 "E=',engy,' GeV"'
      if (xpar8.eq.1.) then
      write(ifhi,'(a,f5.2,a)')  'text 0.1 0.8 "b=',biniDf,' fm"'
      endif
      write(ifhi,'(a)')       'array 2'

      do i=0,nptg
        t=1d0
        if(xpar6.eq.2)t=w(i)
        if(xpar6.eq.1)t=z(i)
        if(xpar6.eq.0)t=1d0
        write(ifhi,*) y(i),omJ5(i)/t
      enddo

      write(ifhi,'(a)')    '  endarray'

      write(ifhi,'(a)')    'closehisto plot 0'


      koll=kollsave

      endif
      end


c----------------------------------------------------------------------
      subroutine xfzero
c----------------------------------------------------------------------

#include "aaa.h"
#include "par.h"
#include "sem.h"
#include "tab.h"

      double precision x,xminr,sh2,sqrtq2s,tmin,tmax,t
      double precision fzeroGlu,fzeroVal,FzeroQua,ffsigj,sh

      !### kw ### #########################
      nq2mn1=1
      nq2mn2=maxq2mx
      if(iscreen.eq.0)nq2mn2=nq2mn1
      do nq2mn=nq2mn1,nq2mn2 ! =========  q2jmin loop =========
      if(nq2mn.le.nq2mn1+1.or.nq2mn.ge.nq2mn2-1)then
      q2jmin=q2mnval(nq2mn)
      !### kw ### #########################
      call setBasicq(q2jmin,1)   !set dels, gamsoft and betff
      s=engy*engy

      nptg=100                   !number of point for the graphs
      xminr=dble(4.*q2jmin)/dble(engy**2)  !value of xminr for plotting the function
c Initialize some needed variables for function F()
c      if(iscreen.ne.0)then
c        s=engy*engy
c        b=0.
c        do i=idxDmin(iomega),idxDmax(iomega)
c         zp=0.
c         zt=0.
c         call Gfunpar(zp,zt,0.,0.,1,i,b,s,alpx,betx,betpx,epsp,epst,epss,gamv)
c        enddo
c      endif
c********************* blue = fzeroGLu *****************************

      write(ifhi,'(a)')       'openhisto name fzeroGlu'
      write(ifhi,'(a)')       'htyp lbu'
      write(ifhi,'(a)')       'xmod lin ymod log'
      write(ifhi,'(a,2e11.3)')'xrange',xminr,xmaxDf
      write(ifhi,'(a)')    'text 0 0 "xaxis x"'
      write(ifhi,'(a)')    'text 0 0 "yaxis f?0! Glu"'
      if (xpar8.eq.1.) then
        write(ifhi,'(a,e7.2,a)')  'text 0.1 0.9 "s=',s,' GeV^2!"'
      endif
      write(ifhi,'(a)')       'array 2'

      do i=0,nptg
        x=xminr
        if (i.ne.0) x=x*(xmaxDf/xminr)**(dble(i)/dble(nptg))
c        x=xminr+(xmaxDf-xminr)*(dble(i)/dble(nptg))
        write(ifhi,*) x,fzeroGlu(x,q2jmin,1)
c.......write(*,*) 'XIe',i
      enddo

      write(ifhi,'(a)')    '  endarray'
      write(ifhi,'(a)')    'closehisto plot 0-'

c********************* green = FzeroQua *****************************

      write(ifhi,'(a)')       'openhisto name FzeroQua'
      write(ifhi,'(a)')       'htyp lgu'
      write(ifhi,'(a)')       'xmod log ymod log'
      write(ifhi,'(a,2e11.3)')'xrange',xminr,xmaxDf
      write(ifhi,'(a)')    'text 0 0 "xaxis x"'
      write(ifhi,'(a)')    'text 0 0 "yaxis f?0! Sea"'
      if (xpar8.eq.1.) then
        write(ifhi,'(a,e7.2,a)')  'text 0.1 0.9 "s=',s,' GeV^2!"'
      endif
      write(ifhi,'(a)')       'array 2'

      do i=0,nptg
        x=xminr
        if (i.ne.0) x=x*(xmaxDf/xminr)**(dble(i)/dble(nptg))
c        x=xminr+(xmaxDf-xminr)*(dble(i)/dble(nptg))
        write(ifhi,*) x,FzeroQua(x,q2jmin,1,999)
c.......write(*,*) 'XIe',i
      enddo

      write(ifhi,'(a)')    '  endarray'
      write(ifhi,'(a)')    'closehisto plot 0-'

c********************* red = fzeroVal (1) *****************************

      write(ifhi,'(a)')       'openhisto name fzeroVal1'
      write(ifhi,'(a)')       'htyp lru'
      write(ifhi,'(a)')       'xmod log ymod log'
      write(ifhi,'(a,2e11.3)')'xrange',xminr,xmaxDf
      write(ifhi,'(a)')    'text 0 0 "xaxis x"'
      write(ifhi,'(a)')    'text 0 0 "yaxis f?0! Val(1)"'
      if (xpar8.eq.1.) then
        write(ifhi,'(a,e7.2,a)')  'text 0.1 0.9 "s=',s,' GeV^2!"'
      endif
      write(ifhi,'(a)')       'array 2'

      do i=0,nptg
        x=xminr
        if (i.ne.0) x=x*(xmaxDf/xminr)**(dble(i)/dble(nptg))
c        x=xminr+(xmaxDf-xminr)*(dble(i)/dble(nptg))
        write(ifhi,*) x,fzeroVal(x,q2jmin,1,1)
c.......write(*,*) 'XIe',i
      enddo

      write(ifhi,'(a)')    '  endarray'
      write(ifhi,'(a)')    'closehisto plot 0-'

c********************* yellow = fzeroVal (2) *****************************

      write(ifhi,'(a)')       'openhisto name fzeroVal2'
      write(ifhi,'(a)')       'htyp lyu'
      write(ifhi,'(a)')       'xmod log ymod log'
      write(ifhi,'(a,2e11.3)')'xrange',xminr,xmaxDf
      write(ifhi,'(a)')    'text 0 0 "xaxis x"'
      write(ifhi,'(a)')    'text 0 0 "yaxis f?0! Val(2)"'
      if (xpar8.eq.1.) then
        write(ifhi,'(a,e7.2,a)')  'text 0.1 0.9 "s=',s,' GeV^2!"'
      endif
      write(ifhi,'(a)')       'array 2'

      do i=0,nptg
        x=xminr
        if (i.ne.0) x=x*(xmaxDf/xminr)**(dble(i)/dble(nptg))
c        x=xminr+(xmaxDf-xminr)*(dble(i)/dble(nptg))
        write(ifhi,*) x,fzeroVal(x,q2jmin,1,2)
c.......write(*,*) 'XIe',i
      enddo

      write(ifhi,'(a)')    '  endarray'
      write(ifhi,'(a)')    'closehisto plot 0'

c********************* red = ffsigj *****************************

      sh=s*dble(xpar4)**2
      if(sh.gt.4.*q2jmin)then
      sh2=dble(sh/2.)
      sqrtq2s=sqrt(dble(q2jmin*sh))
      tmin=sh2-sqrt((sh2-sqrtq2s)*(sh2+sqrtq2s))
      tmax=sh2

      write(ifhi,'(a)')       'openhisto name ffsigj'
      write(ifhi,'(a)')       'htyp lru'
      write(ifhi,'(a)')       'xmod lin ymod lin'
      write(ifhi,'(a,2e11.3)')'xrange',tmin,tmax
      write(ifhi,'(a)')    'text 0 0 "xaxis x"'
      write(ifhi,'(a)')    'text 0 0 "yaxis ffsigj"'
      if (xpar8.eq.1.) then
        write(ifhi,'(a,e7.2,a)')  'text 0.1 0.9 "s=',s,' GeV^2!"'
      endif
      write(ifhi,'(a)')       'array 2'

      do i=0,nptg
        t=tmin
c        if (i.ne.0) t=t*(tmax/tmin)**(dble(i)/dble(nptg))
c        t=tmin/(1d0-dble(i)/dble(nptg)
c     &       *(1d0-tmin/tmax))
        t=tmin+(tmax-tmin)*(dble(i)/dble(nptg))
        qq=sngl(t*(1d0-t/dble(sh)))
        write(ifhi,*) t,ffsigj(t,qq,dble(xpar4),dble(xpar4),1,2,2)
     *        /dble(sh)**2
     *         * dble(2*pi*pssalf(qq/qcdlam))**2
      enddo

      write(ifhi,'(a)')    '  endarray'
      write(ifhi,'(a)')    'closehisto plot 0'
      endif

      ! ### kw ###
      endif
      enddo ! =========  q2jmin loop ========= 
      ! ### kw ###

      end


c----------------------------------------------------------------------
      subroutine xsParD
c----------------------------------------------------------------------
c  omega vs energy for x=1 for given y and b
c----------------------------------------------------------------------

#include "aaa.h"
#include "sem.h"
#include "par.h"

c      double precision x,xminr,y,om51,xtmp,z(0:200)
cc     & ,t,omYuncut
c
c      nptg=50                  !number of point for the graphs
c      biniDf=xpar2                 !value of biniDf (impact parameter)
c      y=dble(xpar3)                 !value of y (rapidity)
c      xtmp=xmaxDf
c      etmp=engy
c      xmaxDf=dexp(-2.d0*y)
c      call Class('xsParD     ')
c      iqq1=0
c      iqq2=4
cc      iqqN=0
c
c      xminr=dble(egylow/engy)**2.d0  !value of xminr for plotting the function
c
cc**********************************************************************
c
c      if (engy.ge.10.) then
c        if (engy.ge.100.) then
c          if (engy.ge.1000.) then
c            if (engy.ge.10000.) then
c           write(ifhi,'(a,I5)')  'openhisto name DExaI-',int(engy)
c            else
c           write(ifhi,'(a,I4)')  'openhisto name DExaI-',int(engy)
c            endif
c          else
c           write(ifhi,'(a,I3)')  'openhisto name DExaI-',int(engy)
c          endif
c        else
c          write(ifhi,'(a,I2)')  'openhisto name DExaI-',int(engy)
c        endif
c      else
c        write(ifhi,'(a,I1)')  'openhisto name DExaI-',int(engy)
c      endif
c      write(ifhi,'(a)')       'htyp lga'
c      if(xpar5.eq.0.)then
c        write(ifhi,'(a)')       'xmod log ymod log'
c      else
c        write(ifhi,'(a)')       'xmod log ymod lin'
c      endif
c      write(ifhi,'(a,2e11.3)')'xrange',xminr,xmaxDf
c      write(ifhi,'(a)')       'yrange auto auto'
c      write(ifhi,'(a)')      'text 0 0 "xaxis x"'
c      write(ifhi,'(a,a)')    'text 0 0 "yaxis D(x+,x-,s,b)" '
c      write(ifhi,'(a,e8.2,a)')  'text 0.1 0.9 "E=',engy,' GeV"'
c      if (xpar8.eq.1.) then
c        write(ifhi,'(a,f5.2,a)')  'text 0.1 0.8 "b=',biniDf,' fm"'
c      endif
c      write(ifhi,'(a)')       'array 2'
c
c      do i=0,nptg
c        x=xminr
c        if (i.ne.0) x=x*(xmaxDf/xminr)**(dble(i)/dble(nptg))
c        engy=sqrt(sngl(x))*etmp
c        iclegy=1+int((log(engy)-log(egylow))/log(egyfac))
c        smaxDf=engy**2
c        z(i)=(alpD(idxD0,iclpro,icltar)
c     &   *(smaxDf*sngl(x))**(betD(idxD0,iclpro,icltar)
c     &          +gamD(idxD0,iclpro,icltar)*biniDf**2.)
c     &     *exp(-biniDf**2./delD(idxD0,iclpro,icltar))
c     &                 +alpD(1,iclpro,icltar)
c     &   *(smaxDf*sngl(x))**(betD(1,iclpro,icltar)
c     &           +gamD(1,iclpro,icltar)*biniDf**2.)
c     &        *exp(-biniDf**2./delD(1,iclpro,icltar)))
c     &   /float(idxD0+1)
c        if(idxD.gt.1)z(i)=z(i)+alpD(idxD,iclpro,icltar)
c     &   *(smaxDf*sngl(x))**(betD(idxD,iclpro,icltar)
c     &           +gamD(idxD,iclpro,icltar)*biniDf**2.)
c     &        *exp(-biniDf**2./delD(idxD,iclpro,icltar))
c        write(ifhi,*) x,om51(1.d0,0.d0,biniDf,iqq1,iqq2)
c      enddo
c
c      write(ifhi,'(a)')    '  endarray'
c
c      if(isetcs.ge.2)then    !otherwise only 1 definition of iclegy available
c
c      write(ifhi,'(a)')    'closehisto plot 0-'
cc**********************************************************************
c
c      if (engy.ge.10.) then
c        if (engy.ge.100.) then
c          if (engy.ge.1000.) then
c            if (engy.ge.10000.) then
c           write(ifhi,'(a,I5)')  'openhisto name DParamI-',int(engy)
c            else
c           write(ifhi,'(a,I4)')  'openhisto name DParamI-',int(engy)
c            endif
c          else
c           write(ifhi,'(a,I3)')  'openhisto name DParamI-',int(engy)
c          endif
c        else
c          write(ifhi,'(a,I2)')  'openhisto name DParamI-',int(engy)
c        endif
c      else
c        write(ifhi,'(a,I1)')  'openhisto name DParamI-',int(engy)
c      endif
c      write(ifhi,'(a)')       'htyp poc'
c      write(ifhi,'(a)')       'xmod log ymod log'
c      write(ifhi,'(a,2e11.3)')'xrange',xminr,xmaxDf
c      write(ifhi,'(a)')       'yrange .01 auto'
c      write(ifhi,'(a)')      'text 0 0 "xaxis x"'
c      write(ifhi,'(a,a)')    'text 0 0 "yaxis D(x+,x-,s,b)" '
c      write(ifhi,'(a,e8.2,a)')  'text 0.1 0.9 "E=',engy,' GeV"'
c      if (xpar8.eq.1.) then
c        write(ifhi,'(a,f5.2,a)')  'text 0.1 0.8 "b=',biniDf,' fm"'
c      endif
c      write(ifhi,'(a)')       'array 2'
c
c      nptg=nptg/2
c      do i=0,nptg
c        x=xminr
c        if (i.ne.0) x=x*(xmaxDf/xminr)**(dble(i)/dble(nptg))
c      write(ifhi,*) x,z(i)
c      enddo
c
c      write(ifhi,'(a)')    '  endarray'
c      endif
c
c      write(ifhi,'(a)')    'closehisto plot 0'
c
c      xmaxDf=xtmp
c      engy=etmp
c      smaxDf=engy**2

      end

c----------------------------------------------------------------------
      subroutine xyParD
c----------------------------------------------------------------------

#include "aaa.h"
#include "ems.h"
#include "sem.h"
#include "par.h"

      double precision x,om51,xDfit,z(0:200),ymax,y,t,t1,t2,t3,t4,xp,xm

      nptg=50                  !number of point for the graphs
      biniDf=xpar2                 !value of biniDf (impact parameter)

      x=dble(xpar4)
      if(x.le.1.d-20)x=min(1d0,50.d0*xminDf)
      ymax=-.5d0*dlog(x)
      q20=q2nmin
      q2pmin(1)=q20+xpar98
      q2pmin(2)=q20+xpar99
      call ipoOm5Tables(1) 
      zz=0
      iqq1=0
      iqq2=-5         !all including all diff contributions with MC
c      iqqN=0


      if (engy.ge.10.) then
       if (engy.ge.100.) then
        if (engy.ge.1000.) then
         if (engy.ge.10000.) then
      write(ifhi,'(a,I5)')       'openhisto name DyParam-',int(engy)
         else
      write(ifhi,'(a,I4)')       'openhisto name DyParam-',int(engy)
        endif
        else
      write(ifhi,'(a,I3)')       'openhisto name DyParam-',int(engy)
       endif
       else
      write(ifhi,'(a,I2)')       'openhisto name DyParam-',int(engy)
      endif
      else
      write(ifhi,'(a,I1)')       'openhisto name DyParam-',int(engy)
      endif
      write(ifhi,'(a)')       'htyp lru'
      if(xpar5.eq.0.)then
        write(ifhi,'(a)')       'xmod lin ymod lin'
      else
        write(ifhi,'(a)')       'xmod lin ymod log'
      endif
      write(ifhi,'(a,2e11.3)')'xrange ',-ymax-1.d0,ymax+1.d0
      write(ifhi,'(a)')       'yrange auto auto'
      write(ifhi,'(a)')      'text 0 0 "xaxis y"'
      write(ifhi,'(a,a)')    'text 0 0 "yaxis D(x+,x-,s,b)" '
      write(ifhi,'(a,e8.2,a)')  'text 0.1 0.9 "E=',engy,' GeV"'
      if (xpar8.eq.1.) then
      write(ifhi,'(a,f5.2,a)')  'text 0.1 0.8 "b=',biniDf,' fm"'
      write(ifhi,'(a,f5.2,a)')  'text 0.1 0.7 "x=',x,'"'
      endif
      write(ifhi,'(a)')       'array 2'

      do i=0,nptg
        y=-ymax+(ymax+ymax)*(dble(i)/dble(nptg))
        xp=sqrt(x)*exp(y)
        xm=sqrt(x)*exp(-y)
        z(i)=xDfit(zz,idxDmin(iomega),idxDmax(iomega),engy**2,sngl(xp)
     .            ,sngl(xm),biniDf)
        write(ifhi,*) y,z(i)
       enddo

      write(ifhi,'(a)')    '  endarray'

      write(ifhi,'(a)')    'closehisto plot 0-'


c**********************************************************************

      if (engy.ge.10.) then
        if (engy.ge.100.) then
          if (engy.ge.1000.) then
            if (engy.ge.10000.) then
           write(ifhi,'(a,I5)')  'openhisto name DyParamI-',int(engy)
            else
           write(ifhi,'(a,I4)')  'openhisto name DyParamI-',int(engy)
            endif
          else
           write(ifhi,'(a,I3)')  'openhisto name DyParamI-',int(engy)
          endif
        else
          write(ifhi,'(a,I2)')  'openhisto name DyParamI-',int(engy)
        endif
      else
        write(ifhi,'(a,I1)')  'openhisto name DyParamI-',int(engy)
      endif
      write(ifhi,'(a)')       'htyp lbo'
      write(ifhi,'(a)')       'xmod lin ymod lin'
      write(ifhi,'(a,2e11.3)')'xrange ',-ymax-1.d0,ymax+1.d0
      write(ifhi,'(a)')       'yrange auto auto'
      write(ifhi,'(a)')      'text 0 0 "xaxis y"'
      write(ifhi,'(a,a)')    'text 0 0 "yaxis D(x+,x-,s,b)" '
      write(ifhi,'(a,e8.2,a)')  'text 0.1 0.9 "E=',engy,' GeV"'
      if (xpar8.eq.1.) then
        write(ifhi,'(a,f5.2,a)')  'text 0.1 0.8 "b=',biniDf,' fm"'
        write(ifhi,'(a,f5.2,a)')  'text 0.1 0.7 "x=',x,'"'
      endif
      write(ifhi,'(a)')       'array 2'

      do i=0,nptg
        y=-ymax+(ymax+ymax)*(dble(i)/dble(nptg))
        xp=sqrt(x)*exp(y)
        xm=sqrt(x)*exp(-y)
        tmp=xDfit(zz,idxD0,-1,engy**2,sngl(xp),sngl(xm),biniDf)
        write(ifhi,*) y,tmp
      enddo

      write(ifhi,'(a)')    '  endarray'
      write(ifhi,'(a)')    'closehisto plot 0-'

c**********************************************************************
      if(iomega.ne.2)then
      write(ifhi,'(a,I5)')  'openhisto name DyParDifp'
      write(ifhi,'(a)')       'htyp lyi'
      write(ifhi,'(a)')       'xmod lin ymod lin'
      write(ifhi,'(a,2e11.3)')'xrange ',-ymax-1.d0,ymax+1.d0
      write(ifhi,'(a)')       'yrange auto auto'
      write(ifhi,'(a)')      '- text 0 0 "xaxis y"'
      write(ifhi,'(a,a)')    '+ text 0 0 "yaxis D(x+,x-,s,b)" '
      write(ifhi,'(a,a)')    '+ text 0 0 "yaxis D(x+,x-,s,b)" '
      write(ifhi,'(a,a)')    '+ text 0 0 "yaxis D(x+,x-,s,b)" '
      write(ifhi,'(a,a)')    '+ text 0 0 "yaxis D(x+,x-,s,b)" '
      write(ifhi,'(a,a)')    '+ text 0 0 "yaxis D(x+,x-,s,b)" '
      write(ifhi,'(a,e8.2,a)')  'text 0.1 0.9 "E=',engy,' GeV"'
      if (xpar8.eq.1.) then
        write(ifhi,'(a,f5.2,a)')  'text 0.1 0.8 "b=',biniDf,' fm"'
        write(ifhi,'(a,f5.2,a)')  'text 0.1 0.7 "x=',x,'"'
      endif
      write(ifhi,'(a)')       'array -6'
      
      do i=0,nptg
        y=-ymax+(ymax+ymax)*(dble(i)/dble(nptg))
        xp=sqrt(x)*exp(y)
        xm=sqrt(x)*exp(-y)
        t=xDfit(zz,idxD1-1,-idxD1,engy**2,sngl(xp),sngl(xm),biniDf)
        t1=om51(x,y,biniDf,-6,-6)/xminDf**epssUni(2)
     &         /(x**dble(-alppar-0.5*(epspUni(2)+epstUni(2))))
        t2=om51(x,y,biniDf,-7,-7)/xminDf**epssUni(3)
     &         /(x**dble(-alppar-0.5*(epspUni(3)+epstUni(3))))
        t3=om51(x,y,biniDf,-8,-8)/xminDf**epssUni(0)
     &         /(x**dble(-alppar-0.5*(epspUni(0)+epstUni(0))))
        t4=om51(x,y,biniDf,-9,-9)/xminDf**epssUni(0)
     &         /(x**dble(-alppar-0.5*(epspUni(0)+epstUni(0))))
        write(ifhi,*) y,t,t1,t2,t3,t4
      enddo

      write(ifhi,'(a)')    '  endarray'
      write(ifhi,'(a)')    'closehisto'
      write(ifhi,'(a)')    'plot -htyp lzo DyParDifp+1-'
      write(ifhi,'(a)')    'plot -htyp lgo DyParDifp+2-'
      write(ifhi,'(a)')    'plot -htyp lgi DyParDifp+3-'
      write(ifhi,'(a)')    'plot -htyp lgu DyParDifp+4-'
      write(ifhi,'(a)')    'plot -htyp lga DyParDifp+5-'

      endif
c**********************************************************************

      if (engy.ge.10.) then
        if (engy.ge.100.) then
          if (engy.ge.1000.) then
            if (engy.ge.10000.) then
           write(ifhi,'(a,I5)')  'openhisto name DyExaI-',int(engy)
            else
           write(ifhi,'(a,I4)')  'openhisto name DyExaI-',int(engy)
            endif
          else
           write(ifhi,'(a,I3)')  'openhisto name DyExaI-',int(engy)
          endif
        else
          write(ifhi,'(a,I2)')  'openhisto name DyExaI-',int(engy)
        endif
      else
        write(ifhi,'(a,I1)')  'openhisto name DyExaI-',int(engy)
      endif
      write(ifhi,'(a)')       'htyp lbu'
      write(ifhi,'(a)')       'xmod lin ymod lin'
      write(ifhi,'(a,2e11.3)')'xrange ',-ymax-1.d0,ymax+1.d0
      write(ifhi,'(a)')       'yrange auto auto'
      write(ifhi,'(a)')      'text 0 0 "xaxis y"'
      write(ifhi,'(a,a)')    'text 0 0 "yaxis D(x+,x-,s,b)" '
      write(ifhi,'(a,e8.2,a)')  'text 0.1 0.9 "E=',engy,' GeV"'
      if (xpar8.eq.1.) then
        write(ifhi,'(a,f5.2,a)')  'text 0.1 0.8 "b=',biniDf,' fm"'
        write(ifhi,'(a,f5.2,a)')  'text 0.1 0.7 "x=',x,'"'
      endif
      write(ifhi,'(a)')       'array 2'

      nptg=nptg/2
      do i=0,nptg
        y=-ymax+(ymax+ymax)*(dble(i)/dble(nptg))
        t=om51(x,y,biniDf,iqq1,iqq2)
     &         /(x**dble(-alppar))
        write(ifhi,*) y,t
      enddo

      write(ifhi,'(a)')    '  endarray'
      write(ifhi,'(a)')    'closehisto plot 0'

      end

c----------------------------------------------------------------------
      subroutine xParSigma
c----------------------------------------------------------------------

#include "aaa.h"
#include "par.h"

      double precision x,w(0:200),z(0:200),xminr,t

      nptg=20                  !number of point for the graphs

      q20=q2nmin
      q2pmin(1)=q20+xpar98
      q2pmin(2)=q20+xpar99
      call ipoOm5Tables(1) 
      xminr=xminDf  !value of xminr for plotting the function
      xshmin=0.1d0

      if (engy.ge.10.) then
       if (engy.ge.100.) then
        if (engy.ge.1000.) then
         if (engy.ge.10000.) then
      write(ifhi,'(a,I5)')    'openhisto name SigmaReel-',int(engy)
         else
      write(ifhi,'(a,I4)')    'openhisto name SigmaReel-',int(engy)
        endif
        else
      write(ifhi,'(a,I3)')    'openhisto name SigmaReel-',int(engy)
       endif
       else
      write(ifhi,'(a,I2)')    'openhisto name SigmaReel-',int(engy)
      endif
      else
      write(ifhi,'(a,I1)')    'openhisto name SigmaReel-',int(engy)
      endif
      write(ifhi,'(a)')       'htyp lbu'
      write(ifhi,'(a)')       'xmod log ymod lin'
      write(ifhi,'(a,2e11.3)')'xrange',xminr,1.
      write(ifhi,'(a)')'yrange auto auto'
      write(ifhi,'(a)')    'text 0 0 "xaxis x"'
      write(ifhi,'(a)')    'text 0 0 "yaxis [s]^2!(X)"'
      if (xpar8.eq.1.) then
      write(ifhi,'(a,e7.2,a)')  'text 0.05 0.75 "Emax=',engy,' GeV"'
      write(ifhi,'(a,2f6.1,a)')  'text 0.05 0.6 "Q?s!^2!=',q2pmin(1)
     *                                                  ,q2pmin(2),'"'
      endif
      write(ifhi,'(a)')       'array 2'

      do i=0,nptg
        X=xminr
        if (i.ne.0) X=X*(xmaxDf/xminr)**(dble(i)/dble(nptg))
        z(i)=dble(sigma2(X,0))
        if(z(i).le.0.) z(i)=0.d0
        write(ifhi,'(2e14.6)') X,z(i)
      enddo

      write(ifhi,'(a)')    '  endarray'
      write(ifhi,'(a)')    'closehisto plot 0-'


c**********************************************************************

      if (engy.ge.10.) then
       if (engy.ge.100.) then
        if (engy.ge.1000.) then
         if (engy.ge.10000.) then
      write(ifhi,'(a,I5)')    'openhisto name SigmaParam-',int(engy)
         else
      write(ifhi,'(a,I4)')    'openhisto name SigmaParam-',int(engy)
        endif
        else
      write(ifhi,'(a,I3)')    'openhisto name SigmaParam-',int(engy)
       endif
       else
      write(ifhi,'(a,I2)')    'openhisto name SigmaParam-',int(engy)
      endif
      else
      write(ifhi,'(a,I1)')    'openhisto name SigmaParam-',int(engy)
      endif
      write(ifhi,'(a)')       'htyp lru'
      write(ifhi,'(a)')       'xmod log ymod lin'
      write(ifhi,'(a,2e11.3)')'xrange',xminr,1.
      write(ifhi,'(a)')'yrange auto auto'
      write(ifhi,'(a)')    'text 0 0 "xaxis x"'
      write(ifhi,'(a)')    'text 0 0 "yaxis [s]^2!?param!(X)"'
      if (xpar8.eq.1.) then
      write(ifhi,'(a,e7.2,a)')  'text 0.05 0.75 "Emax=',engy,' GeV"'
        endif
      write(ifhi,'(a)')       'array 2'

      do i=0,nptg
        X=xminr
        if (i.ne.0) X=X*(xmaxDf/xminr)**(dble(i)/dble(nptg))
        w(i)=dble(xsigmafit(X))
        write(ifhi,'(2e14.6)') X,w(i)
      enddo

      write(ifhi,'(a)')    '  endarray'
      write(ifhi,'(a)')    'closehisto plot 0-'


c**********************************************************************

      if (engy.ge.10.) then
       if (engy.ge.100.) then
        if (engy.ge.1000.) then
         if (engy.ge.10000.) then
      write(ifhi,'(a,I5)')    'openhisto name SigmaInt-',int(engy)
         else
      write(ifhi,'(a,I4)')    'openhisto name SigmaInt-',int(engy)
        endif
        else
      write(ifhi,'(a,I3)')    'openhisto name SigmaInt-',int(engy)
       endif
       else
      write(ifhi,'(a,I2)')    'openhisto name SigmaInt-',int(engy)
      endif
      else
      write(ifhi,'(a,I1)')    'openhisto name SigmaInt-',int(engy)
      endif
      write(ifhi,'(a)')       'htyp lba'
      write(ifhi,'(a)')       'xmod log ymod lin'
      write(ifhi,'(a,2e11.3)')'xrange',xminr,1.
      write(ifhi,'(a)')'yrange auto auto'
      write(ifhi,'(a)')    'text 0 0 "xaxis x"'
      write(ifhi,'(a)')    'text 0 0 "yaxis [s]^2!?Int!(X)"'
      if (xpar8.eq.1.) then
      write(ifhi,'(a,e7.2,a)')  'text 0.05 0.75 "Emax=',engy,' GeV"'
      endif
      write(ifhi,'(a)')       'array 2'

      do i=0,nptg
        X=xminr
        if (i.ne.0) X=X*(xmaxDf/xminr)**(dble(i)/dble(nptg))
        t=dble(sigma1i(X))
        write(ifhi,'(2e14.6)') X,t
      enddo

      write(ifhi,'(a)')    '  endarray'

      if(xpar8.eq.1)then
      write(ifhi,'(a)')    'closehisto plot 0-'


      if (engy.ge.10.) then
       if (engy.ge.100.) then
        if (engy.ge.1000.) then
         if (engy.ge.10000.) then
      write(ifhi,'(a,I5)')    'openhisto name SigmaReelSoft-',int(engy)
         else
      write(ifhi,'(a,I4)')    'openhisto name SigmaReelSoft-',int(engy)
        endif
        else
      write(ifhi,'(a,I3)')    'openhisto name SigmaReelSoft-',int(engy)
       endif
       else
      write(ifhi,'(a,I2)')    'openhisto name SigmaReelSoft-',int(engy)
      endif
      else
      write(ifhi,'(a,I1)')    'openhisto name SigmaReelSoft-',int(engy)
      endif
      write(ifhi,'(a)')       'htyp pft'
      write(ifhi,'(a)')       'xmod log ymod lin'
      write(ifhi,'(a,2e11.3)')'xrange',xminr,1.
      write(ifhi,'(a)')'yrange auto auto'
      write(ifhi,'(a)')    'text 0 0 "xaxis x"'
      write(ifhi,'(a)')    'text 0 0 "yaxis [s]^2!(X)"'
      if (xpar8.eq.1.) then
      write(ifhi,'(a,e7.2,a)')  'text 0.05 0.75 "Emax=',engy,' GeV"'
      endif
      write(ifhi,'(a)')       'array 2'

      do i=0,nptg
        X=xminr
        if (i.ne.0) X=X*(xmaxDf/xminr)**(dble(i)/dble(nptg))
        t=dble(sigma2(X,3))
        if(t.gt.0.) write(ifhi,'(2e14.6)') X,t
      enddo

      write(ifhi,'(a)')    '  endarray'
      write(ifhi,'(a)')    'closehisto plot 0-'


      if (engy.ge.10.) then
       if (engy.ge.100.) then
        if (engy.ge.1000.) then
         if (engy.ge.10000.) then
      write(ifhi,'(a,I5)')    'openhisto name SigmaReelSh-',int(engy)
         else
      write(ifhi,'(a,I4)')    'openhisto name SigmaReelSh-',int(engy)
        endif
        else
      write(ifhi,'(a,I3)')    'openhisto name SigmaReelSh-',int(engy)
       endif
       else
      write(ifhi,'(a,I2)')    'openhisto name SigmaReelSh-',int(engy)
      endif
      else
      write(ifhi,'(a,I1)')    'openhisto name SigmaReelSh-',int(engy)
      endif
      write(ifhi,'(a)')       'htyp pos'
      write(ifhi,'(a)')       'xmod log ymod lin'
      write(ifhi,'(a,2e11.3)')'xrange',xminr,1.
      write(ifhi,'(a)')'yrange auto auto'
      write(ifhi,'(a)')    'text 0 0 "xaxis x"'
      write(ifhi,'(a)')    'text 0 0 "yaxis [s]^2!(X)"'
      if (xpar8.eq.1.) then
      write(ifhi,'(a,e7.2,a)')  'text 0.05 0.75 "Emax=',engy,' GeV"'
      endif
      write(ifhi,'(a)')       'array 2'

      do i=0,nptg
        X=xminr
        if (i.ne.0) X=X*(xmaxDf/xminr)**(dble(i)/dble(nptg))
        t=dble(sigma2(X,2))
        if(t.gt.0.)write(ifhi,'(2e14.6)') X,t
      enddo

      write(ifhi,'(a)')    '  endarray'

      write(ifhi,'(a)')    'closehisto plot 0-'

c**********************************************************************
      write(ifhi,'(a)')    'openhisto'
      write(ifhi,'(a)')       'htyp lya'
      write(ifhi,'(a)')       'xmod log ymod lin'
      write(ifhi,'(a,2e11.3)')'xrange',xminr,1.
      write(ifhi,'(a)')'yrange auto auto'
      write(ifhi,'(a)')    'text 0 0 "xaxis x"'
      write(ifhi,'(a)')    'text 0 0 "yaxis [s]^2!(X)"'
      if (xpar8.eq.1.) then
      write(ifhi,'(a,e7.2,a)')  'text 0.05 0.75 "Emax=',engy,' GeV"'
      endif
      write(ifhi,'(a)')       'array 2'

      do i=0,nptg
        X=xminr
        if (i.ne.0) X=X*(xmaxDf/xminr)**(dble(i)/dble(nptg))
        t=1.d0/dble(-gamD(0,iclpro,icltar)*log(X/xminDf)
     &     +1./delD(0,iclpro,icltar))
        if(t.gt.0.)write(ifhi,'(2e14.6)') X,t
      enddo

      write(ifhi,'(a)')    '  endarray'


      write(ifhi,'(a)')    'closehisto plot 0-'

c**********************************************************************
      write(ifhi,'(a)')    'openhisto'
      write(ifhi,'(a)')       'htyp lyo'
      write(ifhi,'(a)')       'xmod log ymod lin'
      write(ifhi,'(a,2e11.3)')'xrange',xminr,1.
      write(ifhi,'(a)')'yrange auto auto'
      write(ifhi,'(a)')    'text 0 0 "xaxis x"'
      write(ifhi,'(a)')    'text 0 0 "yaxis [s]^2!(X)"'
      if (xpar8.eq.1.) then
      write(ifhi,'(a,e7.2,a)')  'text 0.05 0.75 "Emax=',engy,' GeV"'
      endif
      write(ifhi,'(a)')       'array 2'

      do i=0,nptg
        X=xminr
        if (i.ne.0) X=X*(xmaxDf/xminr)**(dble(i)/dble(nptg))
        t=1.d0/dble(-gamD(1,iclpro,icltar)*log(X/xminDf)
     *     +1./delD(1,iclpro,icltar))
        if(t.gt.0.)write(ifhi,'(2e14.6)') X,t
      enddo

      write(ifhi,'(a)')    '  endarray'

      write(ifhi,'(a)')    'closehisto plot 0-'

c**********************************************************************
      write(ifhi,'(a)')    'openhisto'
      write(ifhi,'(a)')       'htyp lyi'
      write(ifhi,'(a)')       'xmod log ymod lin'
      write(ifhi,'(a,2e11.3)')'xrange',xminr,1.
      write(ifhi,'(a)')'yrange auto auto'
      write(ifhi,'(a)')    'text 0 0 "xaxis x"'
      write(ifhi,'(a)')    'text 0 0 "yaxis [s]^2!(X)"'
      if (xpar8.eq.1.) then
      write(ifhi,'(a,e7.2,a)')  'text 0.05 0.75 "Emax=',engy,' GeV"'
      endif
      write(ifhi,'(a)')       'array 2'

      do i=0,nptg
        X=xminr
        if (i.ne.0) X=X*(xmaxDf/xminr)**(dble(i)/dble(nptg))
        t=1.d0/dble(-gamD(2,iclpro,icltar)*log(X/xminDf)
     *     +1./delD(2,iclpro,icltar))
        if(t.gt.0.)write(ifhi,'(2e14.6)') X,t
      enddo

      write(ifhi,'(a)')    '  endarray'

      endif


      write(ifhi,'(a)')    'closehisto plot 0'

c**********************************************************************
c**********************************************************************

      if (engy.ge.10.) then
       if (engy.ge.100.) then
        if (engy.ge.1000.) then
         if (engy.ge.10000.) then
      write(ifhi,'(a,I5)')    'openhisto name SigmaDiff-',int(engy)
         else
      write(ifhi,'(a,I4)')    'openhisto name SigmaDiff-',int(engy)
        endif
        else
      write(ifhi,'(a,I3)')    'openhisto name SigmaDiff-',int(engy)
       endif
       else
      write(ifhi,'(a,I2)')    'openhisto name SigmaDiff-',int(engy)
      endif
      else
      write(ifhi,'(a,I1)')    'openhisto name SigmaDiff-',int(engy)
      endif
      write(ifhi,'(a)')       'htyp lin'
      write(ifhi,'(a)')       'xmod log ymod lin'
      write(ifhi,'(a,2e11.3)')'xrange',xminr,1.
      write(ifhi,'(a,2e11.3)')'yrange auto auto'
      write(ifhi,'(a)')    'text 0 0 "xaxis x"'
      write(ifhi,'(a)')    'text 0 0 "yaxis [D][s]/[s]"'
      if (xpar8.eq.1.) then
      write(ifhi,'(a,e7.2,a)')  'text 0.05 0.9 "Emax=',engy,' GeV"'
      endif
      write(ifhi,'(a)')       'array 2'

      do i=0,nptg
        X=xminr
        if (i.ne.0) X=X*(xmaxDf/xminr)**(dble(i)/dble(nptg))
        t=0.d0
        if(w(i).gt.0.d0)t=(z(i)-w(i))/w(i)
        if(abs(t).gt.0.15d0) t=dsign(0.15d0,t)
        write(ifhi,'(2e14.6)') X,t
      enddo

      write(ifhi,'(a)')    '  endarray'
      write(ifhi,'(a)')    'closehisto plot 0'

      end

c----------------------------------------------------------------------
      subroutine xParGauss
c----------------------------------------------------------------------

#include "aaa.h"
#include "par.h"

      double precision x,om51,xDfit,y,enh,t!,omNpuncut

      nptg=50                  !number of point for the graphs
      x=dble(xpar4)                  !value of x (energy)
      y=dble(xpar2)                  !value of rapidity
      q2pmin(1)=q2nmin+xpar98
      q2pmin(2)=q2nmin+xpar99
      call ipoOm5Tables(1) 
      zz=0
      if(x.le.1.d-12)x=1.d0

      if (engy.ge.10.) then
       if (engy.ge.100.) then
        if (engy.ge.1000.) then
         if (engy.ge.10000.) then
      write(ifhi,'(a,I5)')    'openhisto name GaussExact-',int(engy)
           else
      write(ifhi,'(a,I4)')    'openhisto name GaussExact-',int(engy)
        endif
        else
      write(ifhi,'(a,I3)')    'openhisto name GaussExact-',int(engy)
       endif
       else
      write(ifhi,'(a,I2)')    'openhisto name GaussExact-',int(engy)
      endif
      else
      write(ifhi,'(a,I1)')    'openhisto name GaussExact-',int(engy)
      endif
      write(ifhi,'(a)')       'htyp lru'
      write(ifhi,'(a)')       'xmod lin ymod lin'
      write(ifhi,'(a,2e11.3)')'xrange',0.,bmaxDf
      write(ifhi,'(a)')'yrange auto auto'
      write(ifhi,'(a)')    'text 0 0 "xaxis b"'
      write(ifhi,'(a)')    'text 0 0 "yaxis b*D(x+,x-,s,b)"'
c      if (xpar8.eq.1.) then
      write(ifhi,'(a,e7.2,a)')  'text 0.6 0.9 "E=',engy,' GeV"'
      write(ifhi,'(a,f5.2,a)')  'text 0.6 0.8 "x=',x,'"'
      write(ifhi,'(a,f5.2,a)')  'text 0.6 0.7 "y=',y,'"'
c      endif
      write(ifhi,'(a)')       'array 2'

      do i=0,nptg
        b=bmaxDf*(float(i)/float(nptg))
        enh=0.d0
        t=dble(b)*(om51(x,y,b,-1,-1)
     &                +enh)
     &         /(x**dble(-alppar))
       write(ifhi,*) b,t
      enddo

      write(ifhi,'(a)')    '  endarray'
      write(ifhi,'(a)')    'closehisto plot 0-'



c**********************************************************************

      if (engy.ge.10.) then
       if (engy.ge.100.) then
        if (engy.ge.1000.) then
         if (engy.ge.10000.) then
      write(ifhi,'(a,I5)')    'openhisto name GaussParam-',int(engy)
         else
      write(ifhi,'(a,I4)')    'openhisto name GaussParam-',int(engy)
        endif
        else
      write(ifhi,'(a,I3)')    'openhisto name GaussParam-',int(engy)
       endif
       else
      write(ifhi,'(a,I2)')    'openhisto name GaussParam-',int(engy)
      endif
      else
      write(ifhi,'(a,I1)')    'openhisto name GaussParam-',int(engy)
      endif
      write(ifhi,'(a)')       'htyp pfc'
      write(ifhi,'(a)')       'xmod log ymod lin'
      write(ifhi,'(a,2e11.3)')'xrange',0.,bmaxDf
      write(ifhi,'(a)')'yrange auto auto'
      write(ifhi,'(a)')    'text 0 0 "xaxis b"'
      write(ifhi,'(a)')    'text 0 0 "yaxis b*D(x+,x-,s,b)"'
c      if (xpar8.eq.1.) then
      write(ifhi,'(a,e7.2,a)')  'text 0.6 0.9 "E=',engy,' GeV"'
      write(ifhi,'(a,f5.2,a)')  'text 0.6 0.8 "x=',x,'"'
      write(ifhi,'(a,f5.2,a)')  'text 0.6 0.7 "y=',y,'"'
c      endif
      write(ifhi,'(a)')       'array 2'

      do i=0,nptg
        b=bmaxDf*(float(i)/float(nptg))
        xp=sqrt(sngl(x))*exp(sngl(y))
        xm=sqrt(sngl(x))*exp(-sngl(y))
        t=xDfit(zz,idxDmin(iomega),idxDmax(iomega),engy**2,xp,xm,b)
        write(ifhi,'(2e14.6)') b,dble(b)*t
      enddo

      write(ifhi,'(a)')    '  endarray'
      write(ifhi,'(a)')    'closehisto plot 0'

      end



c----------------------------------------------------------------------
      subroutine xParOmega1
c----------------------------------------------------------------------

#include "aaa.h"
#include "par.h"

      double precision x,w(0:200),z(0:200)
      double precision yp,om1,xminr,Dsoftshval,t

      nptg=50                 !number of point for the graphs
      biniDf=xpar2               !value of biniDf (impact parameter)
      yp=xpar3                   !value of yp (rapidity)
      q20=q2nmin
      q2pmin(1)=q20+xpar98
      q2pmin(2)=q20+xpar99
      call ipoOm5Tables(1) 
      xminr=1.d0/dble(engy**2)  !value of xminr for plotting the function

      if (engy.ge.10.) then
       if (engy.ge.100.) then
        if (engy.ge.1000.) then
         if (engy.ge.10000.) then
      write(ifhi,'(a,I5)')    'openhisto name Om1Exact-',int(engy)
         else
      write(ifhi,'(a,I4)')    'openhisto name Om1Exact-',int(engy)
        endif
        else
      write(ifhi,'(a,I3)')    'openhisto name Om1Exact-',int(engy)
       endif
       else
      write(ifhi,'(a,I2)')    'openhisto name Om1Exact-',int(engy)
      endif
      else
      write(ifhi,'(a,I1)')    'openhisto name Om1Exact-',int(engy)
      endif
      write(ifhi,'(a)')       'htyp lru'
      write(ifhi,'(a)')       'xmod log ymod log'
      write(ifhi,'(a,2e11.3)')'xrange',xminr,xmaxDf
      write(ifhi,'(a)')       'yrange auto auto'
      write(ifhi,'(a)')    'text 0 0 "xaxis x"'
      write(ifhi,'(a)')    'text 0 0 "yaxis [h](x,0,b)"'
      if (xpar8.eq.1.) then
      write(ifhi,'(a,e7.2,a)')  'text 0.1 0.9 "s=',engy**2,' GeV^2!"'
      write(ifhi,'(a,f5.2,a)')  'text 0.5 0.9 "b=',biniDf,' fm"'
      endif
      write(ifhi,'(a)')       'array 2'

      do i=0,nptg
        x=xminr
        if (i.ne.0) x=x*(xmaxDf/xminr)**(dble(i)/dble(nptg))
        w(i)=Dsoftshval(x,yp,biniDf,0)
     &       *(x**dble(-alppar))
        write(ifhi,*) x,w(i)
      enddo

      write(ifhi,'(a)')    '  endarray'
      write(ifhi,'(a)')    'closehisto plot 0-'


c**********************************************************************

      if (engy.ge.10.) then
       if (engy.ge.100.) then
        if (engy.ge.1000.) then
         if (engy.ge.10000.) then
      write(ifhi,'(a,I5)')    'openhisto name om5param-',int(engy)
         else
      write(ifhi,'(a,I4)')    'openhisto name om5param-',int(engy)
        endif
        else
      write(ifhi,'(a,I3)')    'openhisto name om5param-',int(engy)
       endif
       else
      write(ifhi,'(a,I2)')    'openhisto name om5param-',int(engy)
      endif
      else
      write(ifhi,'(a,I1)')    'openhisto name om5param-',int(engy)
      endif
      write(ifhi,'(a)')       'htyp pfc'
      write(ifhi,'(a)')       'xmod log ymod log'
      write(ifhi,'(a,2e11.3)')'xrange',xminr,xmaxDf
      write(ifhi,'(a)')    'text 0 0 "xaxis x"'
      write(ifhi,'(a)')    'yrange auto auto'
      write(ifhi,'(a)')    'text 0 0 "xaxis x"'
      write(ifhi,'(a)')    'text 0 0 "yaxis [w]?1!(x,0,b)"'
      if (xpar8.eq.1.) then
      write(ifhi,'(a,e7.2,a)')  'text 0.1 0.9 "s=',engy**2,' GeV^2!"'
      write(ifhi,'(a,f5.2,a)')  'text 0.5 0.9 "b=',biniDf,' fm"'
      endif
      write(ifhi,'(a)')       'array 2'

      do i=0,nptg
        x=xminr
        if (i.ne.0) x=x*(xmaxDf/xminr)**(dble(i)/dble(nptg))
        z(i)=om1(x,yp,biniDf)
        write(ifhi,*) x,z(i)
      enddo

      write(ifhi,'(a)')    '  endarray'
      write(ifhi,'(a)')    'closehisto plot 0'


c**********************************************************************
c**********************************************************************

      if (engy.ge.10.) then
       if (engy.ge.100.) then
        if (engy.ge.1000.) then
         if (engy.ge.10000.) then
      write(ifhi,'(a,I5)')    'openhisto name Om1Diff-',int(engy)
         else
      write(ifhi,'(a,I4)')    'openhisto name Om1Diff-',int(engy)
        endif
        else
      write(ifhi,'(a,I3)')    'openhisto name Om1Diff-',int(engy)
       endif
       else
      write(ifhi,'(a,I2)')    'openhisto name Om1Diff-',int(engy)
      endif
      else
      write(ifhi,'(a,I1)')    'openhisto name Om1Diff-',int(engy)
      endif
      write(ifhi,'(a)')       'htyp lru'
      write(ifhi,'(a)')       'xmod log ymod lin'
      write(ifhi,'(a,2e11.3)')'xrange',xminr,xmaxDf
      write(ifhi,'(a,2e11.3)')'yrange auto auto'
      write(ifhi,'(a)')    'text 0 0 "xaxis X"'
      write(ifhi,'(a)')    'text 0 0 "yaxis ([w]?1!-[h])/[h]"'
      if (xpar8.eq.1.) then
      write(ifhi,'(a,e7.2,a)')  'text 0.1 0.9 "s=',engy**2,' GeV^2!"'
      write(ifhi,'(a,f5.2,a)')  'text 0.1 0.8 "b=',biniDf,' fm"'
      endif
      write(ifhi,'(a)')       'array 2'

      do i=0,nptg
        x=xminr
        if (i.ne.0) x=x*(xmaxDf/xminr)**(dble(i)/dble(nptg))
         t=(z(i)-w(i))
         if(abs(w(i)).gt.0.)t=t/w(i)
         if(abs(t).gt.0.15d0) t=dsign(0.15d0,t)
         write(ifhi,*) x,t
      enddo

      write(ifhi,'(a)')    '  endarray'
      write(ifhi,'(a)')    'closehisto plot 0'


      end

c----------------------------------------------------------------------
      subroutine xEpsilon(iii)
c----------------------------------------------------------------------
c iii:  modus (0,1,2)
c----------------------------------------------------------------------
#include "aaa.h"
#include "sem.h"
#include "ems.h"
#include "par.h"

      parameter(nxeps=22,nyeps=32)
      common/cxeps1/w(0:nxeps,nyeps),y1(nyeps),y2(nyeps)
      common/cxeps2/dbp,b1p,b2p
      common/geom/rmproj,rmtarg,bmax,bkmx
      character ch*3
      common /psar34/ rrr,rrrm
      common /psar41/ rrrp,rrrmp
      external pprzz,pttzz

      b1=0.03
      b2=bkmx*1.2
      db=(b2-b1)/nyeps

        if(iii.eq.0)then

      do j=0,nxeps
       do k=1,nyeps
        w(j,k)=0
       enddo
      enddo

        elseif(iii.eq.2)then

      nstat=nfull
      if(nstat.lt.0)nstat=max(1,nevent)

      do nj=1,14
       y1(nj)=1e+20
       y2(nj)=1e-20
       do k=1,nyeps
        if(w(0,k).ne.0)then
         y1(nj)=min(y1(nj),w(nj,k)/w(0,k))
         y2(nj)=max(y2(nj),w(nj,k)/w(0,k))
        endif
       enddo
       if(y1(nj).ge.0)then
         y1(nj)=max(y1(nj)*.2,1e-4)
       else
         y1(nj)=abs(y1(nj)*2.)
       endif
       y2(nj)=min(y2(nj)*5,1e4)
      enddo
      y2(13)=max(y2(13),y2(14))
      y2(14)=max(y2(13),y2(14))
      y2(1)=max(y2(1),y2(2))
      y2(2)=max(y2(1),y2(2))
      y2(5)=max(y2(5),y2(6))
      y2(6)=max(y2(5),y2(6))
      y2(7)=y2(5)
      y2(8)=y2(5)
      y2(9)=max(y2(9),y2(10))
      y2(10)=max(y2(9),y2(10))
      y2(11)=y2(9)
      y2(12)=y2(9)
      do nj=1,16
       if(nj.le.9)write(ifhi,'(a,i1)')'openhisto name xEps',nj
       if(nj.gt.9)write(ifhi,'(a,i2)')'openhisto name xEps',nj
       ch='lin'
       if(nj.eq.7.or.nj.eq.11)ch='lyo'
       if(nj.eq.8.or.nj.eq.12)ch='lgo'
       write(ifhi,'(a)')     'htyp '//ch//' xmod lin'
       if(y1(nj).ge.0.)then
         write(ifhi,'(a)')     'ymod log'
       else
         write(ifhi,'(a)')     'ymod lin'
       endif
       write(ifhi,'(a,e9.2)')'xrange 0 ',b2
       if(nj.eq.1.or.nj.eq.3.or.nj.eq.5.or.nj.eq.9)then
       write(ifhi,'(a,2e9.2)')     'yrange ',min(y1(nj),y1(nj+1))
     *                                      ,max(y2(nj),y2(nj+1))
       else
       write(ifhi,'(a,2e9.2)')     'yrange ',y1(nj),y2(nj)
       endif
       write(ifhi,'(a)')     'text 0 0 "xaxis b"'
       if(nj.eq.1) write(ifhi,'(a)')'txt "yaxis [e]?GP/T!(b)"'
       if(nj.eq.1) write(ifhi,'(a)')'txt "title soft pro   soft tar"'
       if(nj.eq.3) write(ifhi,'(a)')'txt "yaxis [e]?G!(b)"'
       if(nj.eq.3) write(ifhi,'(a)')'txt "title diff"'
       if(nj.eq.3) write(ifhi,'(a)')'txt "title soft   semi"'
       if(nj.eq.5) write(ifhi,'(a)')'txt "yaxis [b]?eff!(b)"'
       if(nj.eq.5) write(ifhi,'(a)')'txt "title soft pro   soft tar"'
       if(nj.eq.9) write(ifhi,'(a)')'txt "yaxis [b]?eff!(b)"'
       if(nj.eq.9) write(ifhi,'(a)')'txt "title semi pro   semi tar"'
       if(nj.eq.13)write(ifhi,'(a)')'txt "yaxis Z?P/T!"'
       write(ifhi,'(a,e22.14)')'histoweight ',dble(nstat)
       write(ifhi,'(a)')       'array 2'
       do k=1,nyeps
        b=b1+(k-0.5)*db
        y=0
        if(w(0,k).ne.0)y=w(nj,k)/w(0,k)
        write(ifhi,'(2e11.3)')b,y
       enddo
       write(ifhi,'(a)')    '  endarray'
       if(nj.eq.2.or.nj.eq.4.or.nj.eq.8.or.nj.eq.12.or.nj.eq.16)then
        write(ifhi,'(a)')    'closehisto plot 0'
       else
        write(ifhi,'(a)')    'closehisto plot 0-'
       endif
      enddo
      !----17-18-19---
      write(ifhi,'(a)') 'openhisto name xEps17'
      write(ifhi,'(a)') 'htyp lin xmod lin ymod lin'
      write(ifhi,'(a)') 'xrange 0 10'
      write(ifhi,'(a)') 'text 0 0 "xaxis b?0!"'
      write(ifhi,'(a)') 'txt "yaxis Z?P/T!(b?0!)"'
      write(ifhi,'(a,e22.14)')'histoweight ',dble(nstat)
      write(ifhi,'(a)') 'array 3'
      do k=1,20
       b=(k-0.5)/2.
       y=0
       if(w(19,k).ne.0)y=w(17,k)/w(19,k)
       write(ifhi,'(3e11.3)')b,y,w(19,k)
      enddo
      write(ifhi,'(a)')  '  endarray'
      write(ifhi,'(a)')  'closehisto plot 0-'

      write(ifhi,'(a)') 'openhisto name xEps18'
      write(ifhi,'(a)') 'htyp lin xmod lin ymod lin'
      write(ifhi,'(a)') 'xrange 0 10'
      write(ifhi,'(a)') 'text 0 0 "xaxis b?0!"'
      write(ifhi,'(a)') 'txt "yaxis Z?P/T!(b?0!)"'
      write(ifhi,'(a,e22.14)')'histoweight ',dble(nstat)
      write(ifhi,'(a)') 'array 3'
      do k=1,10
       b=(k-0.5)/2.
       y=0
       if(w(19,k).ne.0)y=w(18,k)/w(19,k)
       write(ifhi,'(3e11.3)')b,y,w(19,k)
      enddo
      write(ifhi,'(a)')  '  endarray'
      write(ifhi,'(a)')  'closehisto plot 0'
      !----21-21-22---
      kk=2
      do k=3,32
        if(w(18,k).ne.0)kk=k
      enddo
      xmx=(kk-1)/31.*0.1*maproj*matarg
      write(ifhi,'(a)') 'openhisto name xEps20'
      write(ifhi,'(a)') 'htyp lin xmod lin ymod lin'
      write(ifhi,'(a,f10.2)') 'xrange 0 ',xmx
      write(ifhi,'(a)') 'text 0 0 "xaxis n?Gl!"'
      write(ifhi,'(a)') 'txt "yaxis Z?P/T!(n?Gl!)"'
      write(ifhi,'(a,e22.14)')'histoweight ',dble(nstat)
      write(ifhi,'(a)') 'array 2'
      do k=1,32
       x=(k-1.)*0.1*maproj*matarg/(nyeps-1.)
       y=0
       if(w(22,k).ne.0)y=w(20,k)/w(22,k)
       write(ifhi,'(2e11.3)')x,y
      enddo
      write(ifhi,'(a)')  '  endarray'
      write(ifhi,'(a)')  'closehisto plot 0-'

      write(ifhi,'(a)') 'openhisto name xEps21'
      write(ifhi,'(a)') 'htyp lin xmod lin ymod lin'
      write(ifhi,'(a,f10.2)') 'xrange 0 ',xmx
      write(ifhi,'(a)') 'text 0 0 "xaxis n?Gl!"'
      write(ifhi,'(a)') 'txt "yaxis Z?P/T!(n?Gl!)"'
      write(ifhi,'(a,e22.14)')'histoweight ',dble(nstat)
      write(ifhi,'(a)') 'array 2'
      do k=1,32
       x=(k-1.)*0.1*maproj*matarg/(nyeps-1.)
       y=0
       if(w(22,k).ne.0)y=w(21,k)/w(22,k)
       write(ifhi,'(2e11.3)')x,y
      enddo
      write(ifhi,'(a)')  '  endarray'
      write(ifhi,'(a)')  'closehisto plot 0'

        endif

      end

c----------------------------------------------------------------------
      subroutine xZnucTheo
c----------------------------------------------------------------------
c Theoretical mean Z as a function of impact parameter between proj and targ
c----------------------------------------------------------------------
#include "aaa.h"
#include "sem.h"
#include "ems.h"
#include "par.h"

      common/geom/rmproj,rmtarg,bmax,bkmx
      common /psar34/ rrr,rrrm
      common /psar41/ rrrp,rrrmp
      external pprzz,pttzz
      dimension z(0:40)
      double precision xp

      if(iscreen.eq.0)return
      sy=engy**2
      xminr=1./sy  !value of xminr for plotting the function
      xp=1d0
      if(xpar4.gt.xminr)xp=dble(xpar4)

      nptg=20                   !number of point for the graphs
      if(maproj.gt.5)then
        rrrp=radnuc(maproj)/difnuc(maproj)
        rrrmp=rrrp+log(9.)
        nptg=nptg+10                 !number of point for the graphs
      endif
      if(matarg.gt.5)then
        rrr=radnuc(matarg)/difnuc(matarg)
        rrrm=rrr+log(9.)
        nptg=nptg+10                   !number of point for the graphs
      endif


      write(ifhi,'(a)') 'openhisto name xZTheo-pro2'
      write(ifhi,'(a)') 'htyp lin xmod lin ymod lin'
      write(ifhi,*) 'xrange 0 ',bmax
      write(ifhi,'(a)') 'text 0 0 "xaxis b?0!"'
      write(ifhi,'(a)') 'txt "yaxis Z?P!(b?0!)"'
      write(ifhi,'(a)') 'array 2'
      do k=0,nptg
       b=bmax*(float(k)/float(nptg))
       call zzfz(zzp,zzt,kollth,b)
       z(k)=Zpair(xpar98,xpar99,0.,0.,0.,0.,b,2)+zzp
       z(k)=Zpair(xpar98,xpar99,0.,0.,0.,0.,b,1)+zzp
       write(ifhi,'(2e11.3)')b,z(k)
      enddo
      write(ifhi,'(a)')  '  endarray'
      write(ifhi,'(a)')  'closehisto plot 0-'
      write(ifhi,'(a)') 'openhisto name xZTheo-tar2'
      write(ifhi,'(a)') 'htyp lin xmod lin ymod lin'
      write(ifhi,*) 'xrange 0 ',bmax
      write(ifhi,'(a)') 'text 0 0 "xaxis b?0!"'
      write(ifhi,'(a)') 'txt "yaxis Z?T!(b?0!)"'
      write(ifhi,'(a)') 'array 2'
      do k=0,nptg
       b=bmax*(float(k)/float(nptg))
       call zzfz(zzp,zzt,kollth,b)
       z(k)=Zpair(xpar99,xpar98,0.,0.,0.,0.,b,2)+zzt
       write(ifhi,'(2e11.3)')b,z(k)
      enddo
      write(ifhi,'(a)')  '  endarray'
      write(ifhi,'(a)')  'closehisto plot 0-'
      write(ifhi,'(a)') 'openhisto name xZTheo-2'
      write(ifhi,'(a)') 'htyp lba xmod lin ymod lin'
      write(ifhi,*) 'xrange 0 ',bmax
      write(ifhi,'(a)') 'text 0 0 "xaxis b?0!"'
      write(ifhi,'(a)') 'txt "yaxis Z?T!(b?0!)"'
      write(ifhi,'(a)') 'array 2'
      do k=0,nptg
       b=bmax*(float(k)/float(nptg))
c       call zzfz(zzp,zzt,kollth,b)
       z(k)=abs(epscrh)*Zpair(0.,0.,0.,0.,0.,0.,b,-2)+zzt
       write(ifhi,'(2e11.3)')b,z(k)
      enddo
      write(ifhi,'(a)')  '  endarray'
      if(iomega.ne.2)then
      write(ifhi,'(a)')  'closehisto plot 0-'
      write(ifhi,'(a)') 'openhisto name xZTheo-pro0'
      write(ifhi,'(a)') 'htyp lin xmod lin ymod lin'
      write(ifhi,*) 'xrange 0 ',bmax
      write(ifhi,'(a)') 'text 0 0 "xaxis b?0!"'
      write(ifhi,'(a)') 'txt "yaxis Z?P!(b?0!)"'
      write(ifhi,'(a)') 'array 2'
      do k=0,nptg
       b=bmax*(float(k)/float(nptg))
       call zzfz(zzp,zzt,kollth,b)
       z(k)=Zpair(xpar98,xpar99,0.,0.,0.,0.,b,0)/fegypp+zzp
       write(ifhi,'(2e11.3)')b,z(k)
      enddo
      write(ifhi,'(a)')  '  endarray'
      write(ifhi,'(a)')  'closehisto plot 0-'
      write(ifhi,'(a)') 'openhisto name xZTheo-tar0'
      write(ifhi,'(a)') 'htyp lin xmod lin ymod lin'
      write(ifhi,*) 'xrange 0 ',bmax
      write(ifhi,'(a)') 'text 0 0 "xaxis b?0!"'
      write(ifhi,'(a)') 'txt "yaxis Z?T!(b?0!)"'
      write(ifhi,'(a)') 'array 2'
      do k=0,nptg
       b=bmax*(float(k)/float(nptg))
       call zzfz(zzp,zzt,kollth,b)
       z(k)=Zpair(xpar99,xpar98,0.,0.,0.,0.,b,0)/fegypp+zzt
       write(ifhi,'(2e11.3)')b,z(k)
      enddo
      write(ifhi,'(a)')  '  endarray'
      endif
      write(ifhi,'(a)')  'closehisto plot 0'

c      write(ifhi,'(a)') 'openhisto name SatCorr'
c      write(ifhi,'(a)') 'htyp lin xmod lin ymod lin'
c      write(ifhi,*) 'xrange 0 ',bmax
c      write(ifhi,'(a)') 'text 0 0 "xaxis b(fm)"'
c      write(ifhi,'(a)') 'txt "yaxis SatCorr"'
c      write(ifhi,'(a,f5.2,a)')  'text 0.1 0.8 "x^+!=',xp,'"'
c      write(ifhi,'(a)') 'array 2'
c      do k=0,nptg
c       b=bmax*(float(k)/float(nptg))
c       z(k)=SatCorr(xp,1d0,b,0,0.)
c       write(ifhi,'(2e11.3)')b,z(k)
c      enddo
c      write(ifhi,'(a)')  '  endarray'
c      write(ifhi,'(a)')  'closehisto plot 0'

      end

c----------------------------------------------------------------------
      subroutine xParOmegaN
c----------------------------------------------------------------------
c xpar1=engy
c xpar2=b
c xpar4=xremnant
c xpar5: 0=log scale (x dep of om) 1=lin scale (b dep of om)
c----------------------------------------------------------------------

#include "aaa.h"
#include "par.h"
#include "sem.h"
#include "ems.h"
      double precision plc,s
      common/cems5/plc,s
      common/ckopjtg/kopj(mamx),kotg(mamx)
      common /cncl/xproj(mamx),yproj(mamx),zproj(mamx)
     *            ,xtarg(mamx),ytarg(mamx),ztarg(mamx)
      double precision x,w(0:200),z(0:200),xminr,t,ghh
     *,xprem,omGam,Womint,Gammapp,WomGamint,omGamk
c     *,yp,SomY,omYuncut,y,xtmp

      nptg=30                  !number of point for the graphs
      biniDf=xpar2               !value of biniDf (impact parameter)
      xprem=dble(xpar4)            !value of x remnant
      if(xprem.le.1.d-12)xprem=1.d0
      bmax=3.
      xminr=1.d0/dble(engy**2)  !value of xminr for plotting the function


c**********************************************************************

      do i=0,3
        b=bmax*(float(i)/float(3))
        z(i)=1.d0
c        if(xpar5.eq.0.)z(i)=Gammapp(engy**2.,b,1)
      enddo

      write(ifhi,'(a)')    'openhisto name Womint-1'
      write(ifhi,'(a)')       'htyp lru'
      if(xpar5.eq.0.)then
        write(ifhi,'(a)')       'xmod lin ymod log'
      else
        write(ifhi,'(a)')       'xmod lin ymod lin'
      endif
      write(ifhi,'(a,2e11.3)')'xrange',0.,bmax
      write(ifhi,'(a)')    'yrange auto auto'
      write(ifhi,'(a)')    'text 0 0 "xaxis b"'
      write(ifhi,'(a)')    'text 0 0 "yaxis [w]?int!(s,b)"'
      if (xpar8.eq.1.) then
      write(ifhi,'(a,e7.2,a)')  'text 0.1 0.9 "s=',engy**2,' GeV^2!"'
      endif
      write(ifhi,'(a)')       'array 2'

      ierr=0
      ghh=0.d0
      do i=0,3
        b=bmax*(float(i)/float(3))
        w(i)=Womint(engy**2.,b)
        if(b.eq.biniDf)then
          write(*,*)'Womint(',b,',1)=',w(i)
          ghh=ghh+w(i)
        endif
        if(w(i).lt.0.d0.and.ierr.eq.0.and.xpar5.eq.0.)then
          write(*,*)'Warning Womint(1)<0 =',w(i)
          w(i)=-w(i)
          ierr=1
          elseif(w(i).lt.0.d0.and.xpar5.eq.0.)then
            w(i)=-w(i)
          elseif(w(i).ge.0.d0.and.ierr.eq.1.and.xpar5.eq.0.)then
            ierr=0
            write(*,*)'Warning Womint(1)>0 =',w(i)
        endif
        write(ifhi,*) b,w(i)/z(i)
      enddo


      write(ifhi,'(a)')    '  endarray'
      write(ifhi,'(a)')    'closehisto plot 0-'

c**********************************************************************

      write(ifhi,'(a)')    'openhisto'
      write(ifhi,'(a)')       'htyp pfc'
      if(xpar5.eq.0.)then
        write(ifhi,'(a)')       'xmod lin ymod log'
      else
        write(ifhi,'(a)')       'xmod lin ymod lin'
      endif
      write(ifhi,'(a,2e11.3)')'xrange',0.,bmax
      write(ifhi,'(a)')    'yrange auto auto'
      write(ifhi,'(a)')    'text 0 0 "xaxis b"'
      write(ifhi,'(a)')    'text 0 0 "yaxis [w]?int!(s,b)"'
      if (xpar8.eq.1.) then
      write(ifhi,'(a,e7.2,a)')  'text 0.1 0.9 "s=',engy**2,' GeV^2!"'
      endif
      write(ifhi,'(a)')       'array 2'

      do i=0,nptg
        b=bmax*(float(i)/float(nptg))
        z(i)=1.d0
        if(xpar5.eq.0.)then
          z(i)=z(i)+dabs(WomGamint(b))
        endif
      enddo
      ierr=0
      do i=0,nptg
        b=bmax*(float(i)/float(nptg))
        w(i)=WomGamint(b)
        if(w(i).lt.0.d0.and.ierr.eq.0.and.xpar5.eq.0.)then
          write(*,*)'Warning WomGamint(1)<0 =',w(i)
          w(i)=-w(i)
          ierr=1
          elseif(w(i).lt.0.d0.and.xpar5.eq.0.)then
            w(i)=-w(i)
          elseif(w(i).ge.0.d0.and.ierr.eq.1.and.xpar5.eq.0.)then
            ierr=0
            write(*,*)'Warning WomGamint(1)>0 =',w(i)
        endif
        write(ifhi,*) b,w(i)/z(i)
      enddo


      write(ifhi,'(a)')    '  endarray'
      write(ifhi,'(a)')    'closehisto plot 0'

      t=Gammapp(engy**2.,biniDf,1)
      write(*,*)'--> gamma(',biniDf,')=',ghh,t

c**********************************************************************
c**********************************************************************
      s=smaxDf
      xproj(1)=0.
      yproj(1)=0.
      zproj(1)=0.
      lproj(1)=1
      ltarg(1)=1
      kopj(1)=1
      kotg(1)=1
      kproj(1,1)=1
      ktarg(1,1)=1
      do k=1,koll
        bk(k)=biniDf
c        bhpr(1,k)=bk(k)
        iproj(k)=1
        itarg(k)=1
      enddo
      call GfunPark(0)   
      !call GfunParP

      if (engy.ge.10.) then
       if (engy.ge.100.) then
        if (engy.ge.1000.) then
         if (engy.ge.10000.) then
      write(ifhi,'(a,I5)')    'openhisto name xOmNG-',int(engy)
         else
      write(ifhi,'(a,I4)')    'openhisto name xOmNG-',int(engy)
        endif
        else
      write(ifhi,'(a,I3)')    'openhisto name xOmNG-',int(engy)
       endif
       else
      write(ifhi,'(a,I2)')    'openhisto name xOmNG-',int(engy)
      endif
      else
      write(ifhi,'(a,I1)')    'openhisto name xOmNG-',int(engy)
      endif
      write(ifhi,'(a)')       'htyp lru'
      write(ifhi,'(a)')       'xmod log ymod log'
      write(ifhi,'(a,2e11.3)')'xrange',xminr,xmaxDf
      write(ifhi,'(a)')    'yrange auto auto'
      write(ifhi,'(a)')    'text 0 0 "xaxis x-"'
      write(ifhi,'(a)')    'text 0 0 "yaxis [w]?NPi!(x+rem,x-rem,b)"'
      if (xpar8.eq.1.) then
      write(ifhi,'(a,e7.2,a)')  'text 0.1 0.9 "s=',engy**2,' GeV^2!"'
      write(ifhi,'(a,f5.2,a)')  'text 0.5 0.9 "b=',biniDf,' fm"'
      write(ifhi,'(a,f5.2,a)')  'text 0.5 0.85 "x+?rem!=',xprem,'"'
      endif
      write(ifhi,'(a)')       'array 2'

      do i=0,nptg
        x=xminr
        if (i.ne.0) x=x*(xmaxDf/xminr)**(dble(i)/dble(nptg))
        t=dabs(omGam(xprem,x,biniDf))
        write(ifhi,*) x,t
      enddo

      write(ifhi,'(a)')    '  endarray'
      write(ifhi,'(a)')    'closehisto plot 0-'

      if (engy.ge.10.) then
       if (engy.ge.100.) then
        if (engy.ge.1000.) then
         if (engy.ge.10000.) then
      write(ifhi,'(a,I5)')    'openhisto name xOmNG-',int(engy)
         else
      write(ifhi,'(a,I4)')    'openhisto name xOmNG-',int(engy)
        endif
        else
      write(ifhi,'(a,I3)')    'openhisto name xOmNG-',int(engy)
       endif
       else
      write(ifhi,'(a,I2)')    'openhisto name xOmNG-',int(engy)
      endif
      else
      write(ifhi,'(a,I1)')    'openhisto name xOmNG-',int(engy)
      endif
      write(ifhi,'(a)')       'htyp pfc'
      write(ifhi,'(a)')       'xmod log ymod log'
      write(ifhi,'(a,2e11.3)')'xrange',xminr,xmaxDf
      write(ifhi,'(a)')    'yrange auto auto'
      write(ifhi,'(a)')    'text 0 0 "xaxis x-"'
      write(ifhi,'(a)')    'text 0 0 "yaxis [w]?NPi!(x+rem,x-rem,b)"'
      if (xpar8.eq.1.) then
      write(ifhi,'(a,e7.2,a)')  'text 0.1 0.9 "s=',engy**2,' GeV^2!"'
      write(ifhi,'(a,f5.2,a)')  'text 0.5 0.9 "b=',biniDf,' fm"'
      write(ifhi,'(a,f5.2,a)')  'text 0.5 0.85 "x+?rem!=',xprem,'"'
      endif
      write(ifhi,'(a)')       'array 2'

      do i=0,nptg
        x=xminr
        if (i.ne.0) x=x*(xmaxDf/xminr)**(dble(i)/dble(nptg))
        t=dabs(omGamk(1,1,xprem,x,ntymin,ntymax))
        write(ifhi,*) x,t
      enddo

      write(ifhi,'(a)')    '  endarray'

      write(ifhi,'(a)')    'closehisto plot 0'

      end

c----------------------------------------------------------------------
      subroutine xParGampp
c----------------------------------------------------------------------

#include "aaa.h"
#include "par.h"

      double precision Gammapp,GammaGauss,sg,sgmc,Znorm,Zn,t
     *                 ,w(0:200),z(0:200)!,GammaMC

      nptg=2                  !number of point for the graphs
      biniDf=xpar2              !value of biniDf (impact parameter)


c**************************************************

      if (biniDf.lt.1.) then
        k=int(10.*biniDf)
        write(ifhi,'(a,I1)')  'openhisto name Gamma-b0.',k
      else
        write(ifhi,'(a,f3.1)')'openhisto name Gamma-b',biniDf
      endif
      write(ifhi,'(a)')       'htyp lru'
      write(ifhi,'(a)')       'xmod lin ymod lin'
      write(ifhi,'(a,2e11.3)')'xrange',0.,float(nptg)
      write(ifhi,'(a,2e11.3)')'yrange auto auto'
      write(ifhi,'(a)')    'text 0 0 "xaxis m"'
      write(ifhi,'(a)')    'text 0 0 "yaxis [g]?h1h2!"'
      if (xpar8.eq.1.) then
        write(ifhi,'(a,e7.2,a)')  'text 0.1 0.9 "s=',engy**2,' GeV^2!"'
        write(ifhi,'(a,f5.2,a)')  'text 0.5 0.9 "b=',biniDf,' fm"'
      endif
      write(ifhi,'(a)')       'array 2'
      sg=0.d0
      do i=0,nptg
        w(i)=Gammapp(engy**2.,biniDf,i)
        sg=sg+w(i)
        write(ifhi,*) i,w(i)
      write(*,*) 'G12',i,w(i)
      enddo

      write(ifhi,'(a)')    '  endarray'
      write(ifhi,'(a)')    'closehisto plot 0-'

      if (biniDf.lt.1.) then
        k=int(10.*biniDf)
        write(ifhi,'(a,I1)')  'openhisto name GammaMC-b0.',k
      else
        write(ifhi,'(a,f3.1)')'openhisto name GammaMC-b',biniDf
      endif
      write(ifhi,'(a)')       'htyp pfc'
      write(ifhi,'(a)')       'xmod lin ymod lin'
      write(ifhi,'(a,2e11.3)')'xrange',0.,float(nptg)
      write(ifhi,'(a,2e11.3)')'yrange auto auto'
      write(ifhi,'(a)')    'text 0 0 "xaxis m"'
      write(ifhi,'(a)')    'text 0 0 "yaxis [g]?h1h2!"'
      if (xpar8.eq.1.) then
        write(ifhi,'(a,e7.2,a)')  'text 0.1 0.9 "s=',engy**2,' GeV^2!"'
        write(ifhi,'(a,f5.2,a)')  'text 0.5 0.9 "b=',biniDf,' fm"'
      endif
      write(ifhi,'(a)')       'array 2'

      sgmc=0.d0
      do i=0,nptg
        z(i)=GammaGauss(engy**2.,biniDf,i)
        sgmc=sgmc+z(i)
        write(ifhi,*) i,z(i)
        write(*,*) 'G12gauss',i,z(i)
      enddo

      write(ifhi,'(a)')    '  endarray'
      write(ifhi,'(a)')    'closehisto plot 0'

c**********************************************************************
c**********************************************************************
      Zn=Znorm(engy**2,biniDf)

      write(ifhi,'(a)')    'openhisto name GammaDiff'
      write(ifhi,'(a)')       'htyp lru'
      write(ifhi,'(a)')       'xmod lin ymod lin'
      write(ifhi,'(a,2e11.3)')'xrange',0.,float(nptg)
      write(ifhi,'(a,2e11.3)')'yrange auto auto'
      write(ifhi,'(a)')    'text 0 0 "xaxis m"'
      write(ifhi,'(a)')    'text 0 0 "yaxis (G-GMC)/G"'
      if (xpar8.eq.1.) then
      write(ifhi,'(a,e7.2,a)')  'text 0.1 0.9 "s=',engy**2,' GeV^2!"'
      write(ifhi,'(a,f5.2,a)')  'text 0.1 0.8 "b=',biniDf,' fm"'
      endif
      write(ifhi,'(a,f5.2,a)')  'text 0.5 0.9 "[S]?Guncut!=',sg,'"'
      write(ifhi,'(a,f5.2,a)')  'text 0.5 0.8 "[S]?Gcut!=',sgmc,'"'
      write(ifhi,'(a,f5.2,a)')  'text 0.5 0.7 "Z=',Zn,'"'
      write(ifhi,'(a)')       'array 2'

      do i=0,nptg
        if(w(i).ne.0d0) t=(z(i)-w(i))/w(i)
c         if(abs(t).gt.0.5d0) t=dsign(0.5d0,t)
         write(ifhi,*) i,t
      enddo

      write(ifhi,'(a)')    '  endarray'
      write(ifhi,'(a)')    'closehisto plot 0'

      end

c----------------------------------------------------------------------
      subroutine xParPomIncDiff
c----------------------------------------------------------------------

#include "aaa.h"
#include "par.h"

      double precision x,PomIncDiffExact,PomIncDiffUnit,xminr,xm,t
     *                 ,w(0:200),z(0:200)

      nptg=10                  !number of point for the graphs
      biniDf=xpar2              !value of biniDf (impact parameter)
      xm=dble(xpar4)              !value of biniDf (impact parameter)
      xminr=1.d0/dble(engy**2)  !value of xminr for plotting the function
      if(xm.le.1.d-12)xm=1.d0


c********************* red = PomIncXExact *****************************

      if (biniDf.lt.1.) then
        k=int(10.*biniDf)
        write(ifhi,'(a,I1)')  'openhisto name PomIncExact-b0.',k
      else
        write(ifhi,'(a,f3.1)')'openhisto name PomIncExact-b',biniDf
      endif
      write(ifhi,'(a)')       'htyp lru'
      write(ifhi,'(a)')       'xmod log ymod log'
      write(ifhi,'(a,2e11.3)')'xrange',xminr,xmaxDf
      write(ifhi,'(a,2e11.3)')'yrange auto auto'
      write(ifhi,'(a)')    'text 0 0 "xaxis x+"'
      write(ifhi,'(a)')    'text 0 0 "yaxis dn?Pom!/dx+/dx-"'
      if (xpar8.eq.1.) then
        write(ifhi,'(a,e7.2,a)')  'text 0.1 0.9 "s=',engy**2,' GeV^2!"'
        write(ifhi,'(a,f5.2,a)')  'text 0.5 0.9 "b=',biniDf,' fm"'
        write(ifhi,'(a,f5.2,a)')  'text 0.5 0.8 "x-=',xm,'"'
      endif
      write(ifhi,'(a)')       'array 2'

      do i=0,nptg
        x=xminr
        if (i.ne.0) x=x*(xmaxDf/xminr)**(dble(i)/dble(nptg))
c        x=xminr+(xmaxDf-xminr)*(dble(i)/dble(nptg))
        w(i)=PomIncDiffExact(dsqrt(x),dsqrt(x),biniDf)
        write(ifhi,*) x,w(i)
      enddo

      write(ifhi,'(a)')    '  endarray'
      write(ifhi,'(a)')    'closehisto plot 0-'



c************************* dot = PomIncXUnit **************************
c      nptg=50     !number of point for the graphs

      if (biniDf.lt.1.) then
        k=int(10.*biniDf)
        write(ifhi,'(a,I1)')  'openhisto name PomIncUnit-b0.',k
      else
        write(ifhi,'(a,f3.1)')'openhisto name PomIncUnit-b',biniDf
      endif
      write(ifhi,'(a)')       'htyp pfc'
      write(ifhi,'(a)')       'xmod log ymod lin'
      write(ifhi,'(a,2e11.3)')'xrange',xminr,xmaxDf
      write(ifhi,'(a,2e11.3)')'yrange auto auto'
      write(ifhi,'(a)')    'text 0 0 "xaxis x"'
      write(ifhi,'(a)')    'text 0 0 "yaxis dn?Pom!/dx+/dx-"'
      if (xpar8.eq.1.) then
        write(ifhi,'(a,e7.2,a)')  'text 0.1 0.9 "s=',engy**2,' GeV^2!"'
        write(ifhi,'(a,f5.2,a)')  'text 0.5 0.9 "b=',biniDf,' fm"'
        write(ifhi,'(a,f5.2,a)')  'text 0.5 0.8 "x-=',xm,' fm"'
      endif
      write(ifhi,'(a)')       'array 2'

      do i=0,nptg
        x=xminr
        if (i.ne.0) x=x*(xmaxDf/xminr)**(dble(i)/dble(nptg))
c        x=xminr+(xmaxDf-xminr)*(dble(i)/dble(nptg))
        z(i)=PomIncDiffUnit(dsqrt(x),dsqrt(x),biniDf)
        write(ifhi,*) x,z(i)
      enddo

      write(ifhi,'(a)')    '  endarray'
      write(ifhi,'(a)')    'closehisto plot 0'

c**********************************************************************
c**********************************************************************

      write(ifhi,'(a)')    'openhisto name PomIncDiffDiff'
      write(ifhi,'(a)')       'htyp lru'
      write(ifhi,'(a)')       'xmod log ymod lin'
      write(ifhi,'(a,2e11.3)')'xrange',xminr,xmaxDf
      write(ifhi,'(a,2e11.3)')'yrange auto auto'
      write(ifhi,'(a)')    'text 0 0 "xaxis x"'
      write(ifhi,'(a)')    'text 0 0 "yaxis ([w]?5!-G)/G"'
      if (xpar8.eq.1.) then
      write(ifhi,'(a,e7.2,a)')  'text 0.1 0.9 "s=',engy**2,' GeV^2!"'
      write(ifhi,'(a,f5.2,a)')  'text 0.5 0.9 "b=',biniDf,' fm"'
      write(ifhi,'(a,f5.2,a)')  'text 0.5 0.8 "x-=',xm,'"'
      endif
      write(ifhi,'(a)')       'array 2'

      do i=0,nptg
        x=xminr
        if (i.ne.0) x=x*(xmaxDf/xminr)**(dble(i)/dble(nptg))
c        x=xminr+(xmaxDf-xminr)*(dble(i)/dble(nptg))
        t=0.d0
        if(w(i).ne.0d0) t=(z(i)-w(i))/w(i)
c         if(abs(t).gt.0.5d0) t=dsign(0.5d0,t)
         write(ifhi,*) x,t
      enddo

      write(ifhi,'(a)')    '  endarray'
      write(ifhi,'(a)')    'closehisto plot 0'

      end

c----------------------------------------------------------------------
      subroutine xParPomInc
c----------------------------------------------------------------------

#include "aaa.h"
#include "par.h"

      double precision x,PomIncExact,PomIncUnit,xminr,xm,t
     *                 ,w(0:200),z(0:200)

      nptg=10                  !number of point for the graphs
      biniDf=xpar2              !value of biniDf (impact parameter)
      xm=dble(xpar4)              !value of biniDf (impact parameter)
      xminr=1.d0/dble(engy**2)  !value of xminr for plotting the function
      if(xm.le.1.d-12)xm=1.d0


c********************* red = PomIncXExact *****************************

      if (biniDf.lt.1.) then
        k=int(10.*biniDf)
        write(ifhi,'(a,I1)')  'openhisto name PomIncExact-b0.',k
      else
        write(ifhi,'(a,f3.1)')'openhisto name PomIncExact-b',biniDf
      endif
      write(ifhi,'(a)')       'htyp lru'
      write(ifhi,'(a)')       'xmod log ymod log'
      write(ifhi,'(a,2e11.3)')'xrange',xminr,xmaxDf
      write(ifhi,'(a,2e11.3)')'yrange auto auto'
      write(ifhi,'(a)')    'text 0 0 "xaxis x+"'
      write(ifhi,'(a)')    'text 0 0 "yaxis dn?Pom!/dx+/dx-"'
      if (xpar8.eq.1.) then
        write(ifhi,'(a,e7.2,a)')  'text 0.1 0.9 "s=',engy**2,' GeV^2!"'
        write(ifhi,'(a,f5.2,a)')  'text 0.5 0.9 "b=',biniDf,' fm"'
        write(ifhi,'(a,f5.2,a)')  'text 0.5 0.8 "x-=',xm,'"'
      endif
      write(ifhi,'(a)')       'array 2'

      do i=0,nptg
        x=xminr
        if (i.ne.0) x=x*(xmaxDf/xminr)**(dble(i)/dble(nptg))
c        x=xminr+(xmaxDf-xminr)*(dble(i)/dble(nptg))
        w(i)=PomIncExact(dsqrt(x),dsqrt(x),biniDf)
        write(ifhi,*) x,w(i)
      write(*,*) 'Xe',i,w(i)
      enddo

      write(ifhi,'(a)')    '  endarray'
      write(ifhi,'(a)')    'closehisto plot 0-'



c************************* dot = PomIncXUnit **************************
c      nptg=50     !number of point for the graphs

      if (biniDf.lt.1.) then
        k=int(10.*biniDf)
        write(ifhi,'(a,I1)')  'openhisto name PomIncUnit-b0.',k
      else
        write(ifhi,'(a,f3.1)')'openhisto name PomIncUnit-b',biniDf
      endif
      write(ifhi,'(a)')       'htyp pfc'
      write(ifhi,'(a)')       'xmod log ymod lin'
      write(ifhi,'(a,2e11.3)')'xrange',xminr,xmaxDf
      write(ifhi,'(a,2e11.3)')'yrange auto auto'
      write(ifhi,'(a)')    'text 0 0 "xaxis x"'
      write(ifhi,'(a)')    'text 0 0 "yaxis dn?Pom!/dx+/dx-"'
      if (xpar8.eq.1.) then
        write(ifhi,'(a,e7.2,a)')  'text 0.1 0.9 "s=',engy**2,' GeV^2!"'
        write(ifhi,'(a,f5.2,a)')  'text 0.5 0.9 "b=',biniDf,' fm"'
        write(ifhi,'(a,f5.2,a)')  'text 0.5 0.8 "x-=',xm,' fm"'
      endif
      write(ifhi,'(a)')       'array 2'

      do i=0,nptg
        x=xminr
        if (i.ne.0) x=x*(xmaxDf/xminr)**(dble(i)/dble(nptg))
c        x=xminr+(xmaxDf-xminr)*(dble(i)/dble(nptg))
        z(i)=PomIncUnit(dsqrt(x),dsqrt(x),biniDf)
        write(ifhi,*) x,z(i)
        write(*,*) 'Xu',i,z(i)
      enddo

      write(ifhi,'(a)')    '  endarray'
      write(ifhi,'(a)')    'closehisto plot 0'

c**********************************************************************
c**********************************************************************

      write(ifhi,'(a)')    'openhisto name PomIncDiff'
      write(ifhi,'(a)')       'htyp lru'
      write(ifhi,'(a)')       'xmod log ymod lin'
      write(ifhi,'(a,2e11.3)')'xrange',xminr,xmaxDf
      write(ifhi,'(a,2e11.3)')'yrange auto auto'
      write(ifhi,'(a)')    'text 0 0 "xaxis x"'
      write(ifhi,'(a)')    'text 0 0 "yaxis ([w]?5!-G)/G"'
      if (xpar8.eq.1.) then
      write(ifhi,'(a,e7.2,a)')  'text 0.1 0.9 "s=',engy**2,' GeV^2!"'
      write(ifhi,'(a,f5.2,a)')  'text 0.5 0.9 "b=',biniDf,' fm"'
      write(ifhi,'(a,f5.2,a)')  'text 0.5 0.8 "x-=',xm,'"'
      endif
      write(ifhi,'(a)')       'array 2'

      do i=0,nptg
        x=xminr
        if (i.ne.0) x=x*(xmaxDf/xminr)**(dble(i)/dble(nptg))
c        x=xminr+(xmaxDf-xminr)*(dble(i)/dble(nptg))
        t=0.d0
        if(w(i).ne.0d0) t=(z(i)-w(i))/w(i)
c         if(abs(t).gt.0.5d0) t=dsign(0.5d0,t)
         write(ifhi,*) x,t
      enddo

      write(ifhi,'(a)')    '  endarray'
      write(ifhi,'(a)')    'closehisto plot 0'

      end

c----------------------------------------------------------------------
      subroutine xParPomIncX
c----------------------------------------------------------------------

#include "aaa.h"
#include "par.h"

      double precision x,PomIncXExact,PomIncXUnit,xminr,y

      nptg=20                   !number of point for the graphs
      biniDf=xpar2              !value of biniDf (impact parameter)
      xminr=1.d0/dble(engy**2)  !value of xminr for plotting the function


c********************* red = PomIncXExact *****************************

      if (biniDf.lt.1.) then
        k=int(10.*biniDf)
        write(ifhi,'(a,I1)')  'openhisto name PomIncXExact-b0.',k
      else
        write(ifhi,'(a,f3.1)')'openhisto name PomIncXExact-b',biniDf
      endif
      write(ifhi,'(a)')       'htyp lru'
      write(ifhi,'(a)')       'xmod log ymod log'
      write(ifhi,'(a,2e11.3)')'xrange',xminr,xmaxDf
      write(ifhi,'(a,2e11.3)')'yrange auto auto'
      write(ifhi,'(a)')    'text 0 0 "xaxis x"'
      write(ifhi,'(a)')    'text 0 0 "yaxis dn?Pom!/dx(x,b)"'
      if (xpar8.eq.1.) then
        write(ifhi,'(a,e7.2,a)')  'text 0.1 0.9 "s=',engy**2,' GeV^2!"'
        write(ifhi,'(a,f5.2,a)')  'text 0.5 0.9 "b=',biniDf,' fm"'
      endif
      write(ifhi,'(a)')       'array 2'

      do i=0,nptg
        x=xminr
        if (i.ne.0) x=x*(xmaxDf/xminr)**(dble(i)/dble(nptg))
c.......x=xminr+(xmaxDf-xminr)*(dble(i)/dble(nptg))
        y=PomIncXExact(x,biniDf)
        write(ifhi,*) x,y
c      write(*,*) 'Xe',i
      enddo

      write(ifhi,'(a)')    '  endarray'
      write(ifhi,'(a)')    'closehisto plot 0-'



c************************* dot = PomIncXUnit **************************
c      nptg=50     !number of point for the graphs

      if (biniDf.lt.1.) then
        k=int(10.*biniDf)
        write(ifhi,'(a,I1)')  'openhisto name PomIncXUnit-b0.',k
      else
        write(ifhi,'(a,f3.1)')'openhisto name PomIncXUnit-b',biniDf
      endif
      write(ifhi,'(a)')       'htyp pfc'
      write(ifhi,'(a)')       'xmod log ymod lin'
      write(ifhi,'(a,2e11.3)')'xrange',xminr,xmaxDf
      write(ifhi,'(a,2e11.3)')'yrange auto auto'
      write(ifhi,'(a)')    'text 0 0 "xaxis x"'
      write(ifhi,'(a)')    'text 0 0 "yaxis dn?Pom!/dx(x,b)"'
      if (xpar8.eq.1.) then
        write(ifhi,'(a,e7.2,a)')  'text 0.1 0.9 "s=',engy**2,' GeV^2!"'
        write(ifhi,'(a,f5.2,a)')  'text 0.5 0.9 "b=',biniDf,' fm"'
      endif
      write(ifhi,'(a)')       'array 2'

      do i=0,nptg
        x=xminr
        if (i.ne.0) x=x*(xmaxDf/xminr)**(dble(i)/dble(nptg))
c.......x=xminr+(xmaxDf-xminr)*(dble(i)/dble(nptg))
        write(ifhi,*) x,PomIncXUnit(x,biniDf)
        write(*,*) 'Xu',i
      enddo

      write(ifhi,'(a)')    '  endarray'
      write(ifhi,'(a)')    'closehisto plot 0'

      end

c----------------------------------------------------------------------
      subroutine xParPomIncXI
c----------------------------------------------------------------------

#include "aaa.h"
#include "par.h"

      double precision x,xminr
      double precision PomIncXIExact,PomIncXIUnit

      nptg=20                   !number of point for the graphs
      xminr=1.d0/dble(engy**2)  !value of xminr for plotting the function
c      q2min1=xpar98
c      q2min2=xpar99

c*********************red = PomIncXIExact *****************************

      write(ifhi,'(a)')       'openhisto name PomIncXIExact'
      write(ifhi,'(a)')       'htyp lru'
      write(ifhi,'(a)')       'xmod log ymod log'
      write(ifhi,'(a,2e11.3)')'xrange',xminr,xmaxDf
      write(ifhi,'(a)')    'text 0 0 "xaxis x"'
      write(ifhi,'(a)')    'text 0 0 "yaxis d[s]?Pom!/dx(x)"'
      if (xpar8.eq.1.) then
        write(ifhi,'(a,e7.2,a)')  'text 0.1 0.9 "s=',engy**2,' GeV^2!"'
      endif
      write(ifhi,'(a)')       'array 2'

      do i=0,nptg
        x=xminr
        if (i.ne.0) x=x*(xmaxDf/xminr)**(dble(i)/dble(nptg))
c.......x=xminr+(xmaxDf-xminr)*(dble(i)/dble(nptg))
        write(ifhi,*) x,PomIncXIExact(x)
c.......write(*,*) 'XIe',i
      enddo

      write(ifhi,'(a)')    '  endarray'
      write(ifhi,'(a)')    'closehisto plot 0-'

c***************************dot = PomIncXIUnit ************************
c.....nptg=50     !number of point for the graphs

      write(ifhi,'(a)')       'openhisto name PomIncXIUnit'
      write(ifhi,'(a)')       'htyp pfc'
      write(ifhi,'(a)')       'xmod log ymod log'
      write(ifhi,'(a,2e11.3)')'xrange',xminr,xmaxDf
      write(ifhi,'(a)')    'text 0 0 "xaxis x"'
      write(ifhi,'(a)')    'text 0 0 "yaxis d[s]?Pom!/dx(x)"'
      if (xpar8.eq.1.) then
        write(ifhi,'(a,e7.2,a)')  'text 0.1 0.9 "s=',engy**2,' GeV^2!"'
      endif
      write(ifhi,'(a)')       'array 2'

      do i=0,nptg
        x=xminr
        if (i.ne.0) x=x*(xmaxDf/xminr)**(dble(i)/dble(nptg))
c.......x=xminr+(xmaxDf-xminr)*(dble(i)/dble(nptg))
        write(ifhi,*) x,PomIncXIUnit(x)
        write(*,*) 'XIu',i
      enddo

      write(ifhi,'(a)')    '  endarray'
      write(ifhi,'(a)')    'closehisto plot 0'

      end

c----------------------------------------------------------------------
      subroutine xParPomIncP
c----------------------------------------------------------------------

#include "aaa.h"
#include "par.h"

      double precision x,PomIncPUnit,xminr
      double precision PomIncPExact
      
      nptg=30                  !number of point for the graphs
      biniDf=xpar1              !value of biniDf (impact parameter)

      xminr=1.d0/dble(engy**2)  !value of xminr for plotting the function

c*********************red = PomIncPExact *****************************

      if (biniDf.lt.1.) then
        k=int(10.*biniDf)
        write(ifhi,'(a,I1)')  'openhisto name PomIncPExact-b0.',k
      else
        write(ifhi,'(a,f3.1)')'openhisto name PomIncPExact-b',biniDf
      endif
      write(ifhi,'(a)')       'htyp lru'
      write(ifhi,'(a)')       'xmod log ymod log'
      write(ifhi,'(a,2e11.3)')'xrange',xminr,xmaxDf
      write(ifhi,'(a)')    'text 0 0 "xaxis x+"'
      write(ifhi,'(a)')    'text 0 0 "yaxis dn?Pom!/dx+(x+,b)"'
      if (xpar8.eq.1.) then
        write(ifhi,'(a,e7.2,a)')  'text 0.1 0.9 "s=',engy**2,' GeV^2!"'
        write(ifhi,'(a,f5.2,a)')  'text 0.5 0.9 "b=',biniDf,' fm"'
      endif
      write(ifhi,'(a)')       'array 2'

      do i=0,nptg
        x=xminr
        if (i.ne.0) x=x*(xmaxDf/xminr)**(dble(i)/dble(nptg))
        write(ifhi,*) x,PomIncPExact(x,biniDf)
      enddo

      write(ifhi,'(a)')    '  endarray'
      write(ifhi,'(a)')    'closehisto plot 0-'


c**************************dot = PomIncPUnit **************************
c.....nptg=50     !number of point for the graphs

      if (biniDf.lt.1.) then
        k=int(10.*biniDf)
        write(ifhi,'(a,I1)')  'openhisto name PomIncPUnit-b0.',k
      else
        write(ifhi,'(a,f3.1)')'openhisto name PomIncPUnit-b',biniDf
      endif
      write(ifhi,'(a)')       'htyp pfc'
      write(ifhi,'(a)')       'xmod log ymod log'
      write(ifhi,'(a,2e11.3)')'xrange',xminr,xmaxDf
      write(ifhi,'(a)')    'text 0 0 "xaxis x+"'
      write(ifhi,'(a)')    'text 0 0 "yaxis dn?Pom!/dx+(x+,b)"'
      if (xpar8.eq.1.) then
        write(ifhi,'(a,e7.2,a)')  'text 0.1 0.9 "s=',engy**2,' GeV^2!"'
        write(ifhi,'(a,f5.2,a)')  'text 0.5 0.9 "b=',biniDf,' fm"'
      endif
      write(ifhi,'(a)')       'array 2'

      do i=0,nptg
        x=xminr
        if (i.ne.0) x=x*(xmaxDf/xminr)**(dble(i)/dble(nptg))
        write(*,*) 'Pu',i
        write(ifhi,*) x,PomIncPUnit(x,biniDf)
      enddo

      write(ifhi,'(a)')    '  endarray'
      write(ifhi,'(a)')    'closehisto plot 0'

      end

c----------------------------------------------------------------------
      subroutine xParPomIncPI
c----------------------------------------------------------------------

#include "aaa.h"
#include "par.h"

      double precision x,xminr
      double precision PomIncPIExact,PomIncPIUnit

      nptg=50                  !number of point for the graphs
      xminr=1.d0/dble(engy**2)  !value of xminr for plotting the function
c      q2min1=xpar98
c      q2min2=xpar99

c*********************red = PomIncPIExact *****************************
c.....nptg=100     !number of point for the graphs

      write(ifhi,'(a)')       'openhisto name PomIncPIExact'
      write(ifhi,'(a)')       'htyp lru'
      write(ifhi,'(a)')       'xmod log ymod log'
      write(ifhi,'(a,2e11.3)')'xrange',xminr,xmaxDf
      write(ifhi,'(a)')    'text 0 0 "xaxis x+"'
      write(ifhi,'(a)')    'text 0 0 "yaxis d[s]?Pom!/dx+(x+)"'
      if (xpar8.eq.1.) then
        write(ifhi,'(a,e7.2,a)')  'text 0.1 0.9 "s=',engy**2,' GeV^2!"'
      endif
      write(ifhi,'(a)')       'array 2'

      do i=0,nptg
        x=xminr
        if (i.ne.0) x=x*(xmaxDf/xminr)**(dble(i)/dble(nptg))
        write(ifhi,*) x,PomIncPIExact(x)
      enddo

      write(ifhi,'(a)')    '  endarray'
      write(ifhi,'(a)')    'closehisto plot 0-'

c***************************dot = PomIncPIUnit ************************
c.....nptg=10     !number of point for the graphs

      write(ifhi,'(a)')       'openhisto name PomIncPIUnit'
      write(ifhi,'(a)')       'htyp pfc'
      write(ifhi,'(a)')       'xmod log ymod log'
      write(ifhi,'(a,2e11.3)')'xrange',xminr,xmaxDf
      write(ifhi,'(a)')    'text 0 0 "xaxis x+"'
      write(ifhi,'(a)')    'text 0 0 "yaxis n?Pom!/dx+(x+)"'
      if (xpar8.eq.1.) then
        write(ifhi,'(a,e7.2,a)')  'text 0.1 0.9 "s=',engy**2,' GeV^2!"'
      endif
      write(ifhi,'(a)')       'array 2'

      do i=0,nptg
        x=xminr
        if (i.ne.0) x=x*(xmaxDf/xminr)**(dble(i)/dble(nptg))
        write(ifhi,*) x,PomIncPIUnit(x)
        write(*,*) 'PIu',i
      enddo

      write(ifhi,'(a)')    '  endarray'
      write(ifhi,'(a)')    'closehisto plot 0'

      end

c----------------------------------------------------------------------
      subroutine xParPomIncM
c----------------------------------------------------------------------

#include "aaa.h"
#include "par.h"

      double precision x,xminr,PomIncMUnit,PomIncMExact

      nptg=50                  !number of point for the graphs
      biniDf=xpar1              !value of biniDf (impact parameter)
      xminr=1.d0/dble(engy**2)  !value of xminr for plotting the function

c**********************red = PomIncMExact *****************************

      if (biniDf.lt.1.) then
        k=int(10.*biniDf)
        write(ifhi,'(a,I1)')  'openhisto name PomIncMExact-b0.',k
      else
        write(ifhi,'(a,f3.1)')'openhisto name PomIncMExact-b',biniDf
      endif
      write(ifhi,'(a)')       'htyp lru'
      write(ifhi,'(a)')       'xmod log ymod log'
      write(ifhi,'(a,2e11.3)')'xrange',xminr,xmaxDf
      write(ifhi,'(a)')    'text 0 0 "xaxis x-"'
      write(ifhi,'(a)')    'text 0 0 "yaxis dn?Pom!/dx-(x-,b)"'
      if (xpar8.eq.1.) then
        write(ifhi,'(a,e7.2,a)')  'text 0.1 0.9 "s=',engy**2,' GeV^2!"'
        write(ifhi,'(a,f5.2,a)')  'text 0.5 0.9 "b=',biniDf,' fm"'
      endif
      write(ifhi,'(a)')       'array 2'

      do i=0,nptg
        x=xminr
        if (i.ne.0) x=x*(xmaxDf/xminr)**(dble(i)/dble(nptg))
        write(ifhi,*) x,PomIncMExact(x,biniDf)
      enddo

      write(ifhi,'(a)')    '  endarray'
      write(ifhi,'(a)')    'closehisto plot 0-'


c**************************dot = PomIncMUnit **************************
c.....nptg=100     !number of point for the graphs

      if (biniDf.lt.1.) then
        k=int(10.*biniDf)
        write(ifhi,'(a,I1)')  'openhisto name PomIncMUnit-b0.',k
      else
        write(ifhi,'(a,f3.1)')'openhisto name PomIncMUnit-b',biniDf
      endif
      write(ifhi,'(a)')       'htyp pfc'
      write(ifhi,'(a)')       'xmod log ymod log'
      write(ifhi,'(a,2e11.3)')'xrange',xminr,xmaxDf
      write(ifhi,'(a)')    'text 0 0 "xaxis x-"'
      write(ifhi,'(a)')    'text 0 0 "yaxis dn?Pom!/dx-(x-,b)"'
      if (xpar8.eq.1.) then
        write(ifhi,'(a,e7.2,a)')  'text 0.1 0.9 "s=',engy**2,' GeV^2!"'
        write(ifhi,'(a,f5.2,a)')  'text 0.5 0.9 "b=',biniDf,' fm"'
      endif
      write(ifhi,'(a)')       'array 2'

      do i=0,nptg
        x=xminr
        if (i.ne.0) x=x*(xmaxDf/xminr)**(dble(i)/dble(nptg))
        write(ifhi,*) x,PomIncMUnit(x,biniDf)
        write(*,*) 'Mu',i
      enddo

      write(ifhi,'(a)')    '  endarray'
      write(ifhi,'(a)')    'closehisto plot 0'

      end

c----------------------------------------------------------------------
      subroutine xParPomIncMI
c----------------------------------------------------------------------

#include "aaa.h"
#include "par.h"

      double precision x,xminr
      double precision PomIncMIExact,PomIncMIUnit
c      q2min1=xpar98
c      q2min2=xpar99

      nptg=30                 !number of point for the graphs
      xminr=1.d0/dble(engy**2)  !value of xminr for plotting the function

c*********************red = PomIncMIExact *****************************
c.....nptg=100     !number of point for the graphs

      write(ifhi,'(a)')       'openhisto name PomIncMIExact'
      write(ifhi,'(a)')       'htyp lru'
      write(ifhi,'(a)')       'xmod log ymod log'
      write(ifhi,'(a,2e11.3)')'xrange',xminr,xmaxDf
      write(ifhi,'(a)')    'text 0 0 "xaxis x-"'
      write(ifhi,'(a)')    'text 0 0 "yaxis d[s]?Pom!/dx-(x-)"'
      if (xpar8.eq.1.) then
        write(ifhi,'(a,e7.2,a)')  'text 0.1 0.9 "s=',engy**2,' GeV^2!"'
      endif
      write(ifhi,'(a)')       'array 2'

      do i=0,nptg
        x=xminr
        if (i.ne.0) x=x*(xmaxDf/xminr)**(dble(i)/dble(nptg))
        write(ifhi,*) x,PomIncMIExact(x)
      enddo

      write(ifhi,'(a)')    '  endarray'
      write(ifhi,'(a)')    'closehisto plot 0-'

c***************************dot = PomIncMIUnit ************************
c.....nptg=100     !number of point for the graphs

      write(ifhi,'(a)')       'openhisto name PomIncMIUnit'
      write(ifhi,'(a)')       'htyp pfc'
      write(ifhi,'(a)')       'xmod log ymod log'
      write(ifhi,'(a,2e11.3)')'xrange',xminr,xmaxDf
      write(ifhi,'(a)')    'text 0 0 "xaxis x-"'
      write(ifhi,'(a)')    'text 0 0 "yaxis d[s]?Pom!/dx-(x-)"'
      if (xpar8.eq.1.) then
        write(ifhi,'(a,e7.2,a)')  'text 0.1 0.9 "s=',engy**2,' GeV^2!"'
      endif
      write(ifhi,'(a)')       'array 2'

      do i=0,nptg
        x=xminr
        if (i.ne.0) x=x*(xmaxDf/xminr)**(dble(i)/dble(nptg))
        write(ifhi,*) x,PomIncMIUnit(x)
        write(*,*) 'MIu',i
      enddo

      write(ifhi,'(a)')    '  endarray'
      write(ifhi,'(a)')    'closehisto plot 0'

      end

c----------------------------------------------------------------------
      subroutine xParPomIncJ
c----------------------------------------------------------------------

#include "aaa.h"
#include "par.h"
      common/geom/rmproj,rmtarg,bmax,bkmx
      parameter(nbib=48)

      double precision PomIncJExact,PomIncJUnit


      b1=0
      b2=bkmx*1.2
      db=(b2-b1)/nbib

c*************************red = PomIncJExact **************************

      write(ifhi,'(a)')       'openhisto name PomIncJExact'
      write(ifhi,'(a)')       'htyp lru xmod lin ymod lin'
      write(ifhi,'(a)')    'text 0 0 "xaxis  impact parameter b (fm)"'
      write(ifhi,'(a)')    'text 0 0 "yaxis  n?Pom!(b)"'
      if (xpar8.eq.1.) then
        write(ifhi,'(a,e7.2,a)')'text 0.5 0.8 "s=',engy**2,' GeV^2!"'
      endif
      write(ifhi,'(a)')       'array 2'
      do k=1,nbib
        b=b1+(k-0.5)*db
        write(ifhi,*)b,PomIncJExact(b)
      enddo
      write(ifhi,'(a)')    '  endarray'
      write(ifhi,'(a)')    'closehisto plot 0-'

c****************************dot = PomIncJUnit ***********************

      write(ifhi,'(a)')       'openhisto name PomIncJUnit'
      write(ifhi,'(a)')       'htyp pfc xmod lin ymod lin'
      write(ifhi,'(a)')    'text 0 0 "xaxis  impact parameter b (fm)"'
      write(ifhi,'(a)')    'text 0 0 "yaxis  n?Pom!(b)"'
      if (xpar8.eq.1.) then
        write(ifhi,'(a,e7.2,a)')'text 0.5 0.8 "s=',engy**2,' GeV^2!"'
      endif
      write(ifhi,'(a)')       'array 2'
      do k=1,nbib
        b=b1+(k-0.5)*db
        write(ifhi,*)b,PomIncJUnit(b)
        write(*,*) 'Ju',k
      enddo
      write(ifhi,'(a)')    '  endarray'
      write(ifhi,'(a)')    'closehisto plot 0'


      end

c----------------------------------------------------------------------
      subroutine xParPhi1
c----------------------------------------------------------------------

#include "aaa.h"
#include "ems.h"
#include "sem.h"
#include "par.h"
      double precision x,xminr,y
      double precision PhiExpo
      double precision PhiExact

      nptg=30                  !number of point for the graphs
      biniDf=xpar2              !value of biniDf (impact parameter)
      xminr=1.d0/dble(engy**2)  !value of xminr for plotting the function
      zz=0

c************************* red = PhiMExact ***************************

      write(ifhi,'(a)')       'openhisto name Phi1Exact'
      write(ifhi,'(a)')       'htyp lru'
      write(ifhi,'(a)')      'xmod log ymod lin'
      write(ifhi,'(a,2e11.3)')'xrange',xminr,xmaxDf
      write(ifhi,'(a)')       'yrange auto auto'
      write(ifhi,'(a)')    'text 0 0 "xaxis x"'
      write(ifhi,'(a)') 'txt "yaxis  [F](0.5,x)/x^[a]!"'
      write(ifhi,'(a,i4,a,f4.1,a)')
     * 'txt  "title E=',nint(engy),' b=',biniDf,'"'
      write(ifhi,'(a)')       'array 2'

      do i=0,nptg
        x=xminr
        if (i.ne.0) x=x*(xmaxDf/xminr)**(dble(i)/dble(nptg))
c        x=xminr+(xmaxDf-xminr)*(dble(i)/dble(nptg))
        y=0.d0
        if(engy**2..lt.5.e06)
     &  y=Phiexact(zz,zz,.5,dsqrt(x),dsqrt(x),smaxDf,biniDf)
     &       *dsqrt(x)**dble(-alplea(iclpro))
     &       *dsqrt(x)**dble(-alplea(icltar))
        write(ifhi,*)x,max(-1d0,y)
      enddo

      write(ifhi,'(a)')    '  endarray'
      write(ifhi,'(a)')    'closehisto plot 0-'

c********************** blue = PhiMExpo ******************************

      write(ifhi,'(a)')       'openhisto name Phi1Expo'
      write(ifhi,'(a)')       'htyp lba'
      write(ifhi,'(a)')       'xmod log ymod lin'
      write(ifhi,'(a,2e11.3)')'xrange',xminr,xmaxDf
      write(ifhi,'(a,2e11.3)')'yrange auto auto'
       write(ifhi,'(a)')       'array 2'

      do i=0,nptg
        x=xminr
        if (i.ne.0) x=x*(xmaxDf/xminr)**(dble(i)/dble(nptg))
c        x=xminr+(xmaxDf-xminr)*(dble(i)/dble(nptg))
        y=Phiexpo(zz,zz,.5,dsqrt(x),dsqrt(x),smaxDf,biniDf)
     &       *dsqrt(x)**dble(-alplea(iclpro))
     &       *dsqrt(x)**dble(-alplea(icltar))
        write(ifhi,*) x,y
      enddo

      write(ifhi,'(a)')    '  endarray'
      write(ifhi,'(a)')    'closehisto plot 0'

      end


cc----------------------------------------------------------------------
c      subroutine xParPhi2
cc----------------------------------------------------------------------
c
c#include "aaa.h"
c#include "ems.h"
c#include "sem.h"
c#include "par.h"
c      double precision x,xminr,xm,y,u(0:100),v(0:100),w(0:100)!,z
c      double precision PhiExpo,omGam,PhiExpoK
c      double precision PhiExact
c      common /cncl/xproj(mamx),yproj(mamx),zproj(mamx)
c     *            ,xtarg(mamx),ytarg(mamx),ztarg(mamx)
c
c      nptg=30                  !number of point for the graphs
c      biniDf=xpar2              !value of biniDf (impact parameter)
c      xminr=1.d-3 !/dble(engy**2)  !value of xminr for plotting the function
c      xm=dble(xpar4)
c
cc************************* yellow = PhiMExact ***************************
c
c      write(ifhi,'(a)')       'openhisto name Phi1Exact'
c      write(ifhi,'(a)')       'htyp lru'
c      write(ifhi,'(a)')      'xmod lin ymod lin'
c      write(ifhi,'(a,2e11.3)')'xrange',xminr,xmaxDf
c      write(ifhi,'(a)')       'yrange auto auto'
c      write(ifhi,'(a)')    'text 0 0 "xaxis x+"'
c      write(ifhi,'(a)') 'txt "yaxis  [F](x)/x^[a]!"'
c      write(ifhi,'(a,i4,a,f4.1,a)')
c     * 'txt  "title E=',nint(engy),' b=',biniDf,'"'
c      write(ifhi,'(a)')       'array 2'
c
c      do i=0,nptg
c       ! x=xminr
c       ! if (i.ne.0) x=x*(xmaxDf/xminr)**(dble(i)/dble(nptg))
c        x=xminr+(xmaxDf-xminr)*(dble(i)/dble(nptg))
c        y=0.d0
c        if(engy**2..lt.5.e06)
c     &  y=Phiexact(0.,0.,1.,dsqrt(x),dsqrt(x),engy**2,biniDf)
c    ! &       *dsqrt(x)**dble(-alplea(iclpro))
c    ! &       *dsqrt(x)**dble(-alplea(icltar))
c        write(ifhi,*)x,y
c      enddo
c
c      write(ifhi,'(a)')    '  endarray'
c      write(ifhi,'(a)')    'closehisto plot 0-'
c
cc********************** blue = PhiMExpo ******************************
c
c      write(ifhi,'(a)')       'openhisto name Phi1Expo'
c      write(ifhi,'(a)')       'htyp lbu'
c      write(ifhi,'(a)')       'xmod lin ymod lin'
c      write(ifhi,'(a,2e11.3)')'xrange',xminr,xmaxDf
c      write(ifhi,'(a,2e11.3)')'yrange auto auto'
c       write(ifhi,'(a)')       'array 2'
c
c      do i=0,nptg
c       ! x=xminr
c       ! if (i.ne.0) x=x*(xmaxDf/xminr)**(dble(i)/dble(nptg))
c        x=xminr+(xmaxDf-xminr)*(dble(i)/dble(nptg))
c        y=Phiexpo(0.,0.,1.,dsqrt(x),dsqrt(x),engy**2,biniDf)
c  !   &       *dsqrt(x)**dble(-alplea(iclpro))
c  !   &       *dsqrt(x)**dble(-alplea(icltar))
c        write(ifhi,*) x,y
c      enddo
c
c      write(ifhi,'(a)')    '  endarray'
c      write(ifhi,'(a)')    'closehisto plot 0-'
c
cc**********************************************************************
cc**********************************************************************
c      do k=1,koll
c        xproj(1)=0.
c        yproj(1)=0.
c        zproj(1)=0.
c        iproj(k)=1
c        itarg(k)=1
c        bk(k)=biniDf
c      enddo
c      call GfunPark(0)
c      call integom1(0)
c
cc********************* points = PhiExpoK*********************************
c
c      if (biniDf.lt.1.) then
c        k=int(10.*biniDf)
c        write(ifhi,'(a,I1)')  'openhisto name PhiExpok-b0.',k
c      else
c        write(ifhi,'(a,f3.1)')  'openhisto name PhiExpok-b',biniDf
c      endif
c      write(ifhi,'(a)')       'htyp pfc'
c      write(ifhi,'(a)')       'xmod log ymod lin'
c      write(ifhi,'(a,2e11.3)')'xrange',xminr,xmaxDf
c      write(ifhi,'(a,2e11.3)')'yrange auto auto'
c      write(ifhi,'(a)')    'text 0 0 "xaxis x"'
c      write(ifhi,'(a)') 'text 0 0.1 "yaxis  [F](x+,x-)/x^[a]?remn!!"'
c      if (xpar8.eq.1.) then
c        write(ifhi,'(a,e7.2,a)')'text 0.1 0.9 "s=',engy**2,' GeV^2!"'
c        write(ifhi,'(a,f5.2,a)')'text 0.5 0.9 "b=',biniDf,' fm"'
c      endif
c      write(ifhi,'(a)')       'array 2'
c
c      do i=0,nptg
c        x=xminr
c        if (i.ne.0) x=x*(xmaxDf/xminr)**(dble(i)/dble(nptg))
c        y=PhiExpoK(1,dsqrt(x),dsqrt(x))
c        write(ifhi,*) x,y
c      enddo
c
c      write(ifhi,'(a)')    '  endarray'
c      write(ifhi,'(a)')    'closehisto plot 0'
c
cc************************* red = PhiMExact*omGam ************************
c
c      write(ifhi,'(a)')       'openhisto name GPhiExact'
c      write(ifhi,'(a)')       'htyp lru'
c      if(xpar5.eq.0.)then
c        write(ifhi,'(a)')       'xmod lin ymod lin'
c      else
c        write(ifhi,'(a)')       'xmod log ymod log'
c      endif
c      write(ifhi,'(a,2e11.3)')'xrange',xminr,xmaxDf
c      write(ifhi,'(a)')       'yrange auto auto'
c      write(ifhi,'(a)')    'text 0 0 "xaxis x+"'
c      write(ifhi,'(a,a)') 'text 0 0.1 "yaxis  G(x+,x-)*[F]'
c     *,'(x+,x-)/x^[a]?remn!!"'
c      if (xpar8.eq.1.) then
c        write(ifhi,'(a,e7.2,a)')'text 0.1 0.2 "s=',engy**2,' GeV^2!"'
c        write(ifhi,'(a,f5.2,a)')'text 0.1 0.1 "b=',biniDf,' fm"'
c        write(ifhi,'(a,f5.2,a)')'text 0.1 0.3 "x-=',xm,'"'
c      endif
c      write(ifhi,'(a)')       'array 2'
c
c      do i=0,nptg
c        x=xminr
c        if(xpar5.ne.0.)then
c          if (i.ne.0) x=x*(xmaxDf/xminr)**(dble(i)/dble(nptg))
c        else
c          x=xminr+(xmaxDf-xminr)*(dble(i)/dble(nptg))
c        endif
cc        z=1.d0-dsqrt(x)
c        v(i)=0.d0
c        if(engy**2..lt.5.e06)
c     *  v(i)=Phiexact(0.,0.,1.,1.d0-x,1.d0-xm,engy**2,biniDf)
cc     *  v(i)=Phiexact(0.,0.,1.,z,z,engy**2,biniDf)
c        u(i)=omGam(x,xm,biniDf)
cc        u(i)=omGam(dsqrt(x),dsqrt(x),biniDf)
c        y=u(i)*v(i)
c        if(xpar5.ne.0.)y=dabs(y)
c        write(ifhi,*)x,y
c      enddo
c
c      write(ifhi,'(a)')    '  endarray'
c      write(ifhi,'(a)')    'closehisto plot 0-'
c
cc************************* red = PhiMExpo*omGam ************************
c
c      write(ifhi,'(a)')       'openhisto name GPhiExpo'
c      write(ifhi,'(a)')       'htyp lba'
c      if(xpar5.eq.0.)then
c        write(ifhi,'(a)')       'xmod lin ymod lin'
c      else
c        write(ifhi,'(a)')       'xmod log ymod log'
c      endif
c      write(ifhi,'(a,2e11.3)')'xrange',xminr,xmaxDf
c      write(ifhi,'(a)')       'yrange auto auto'
c      write(ifhi,'(a)')    'text 0 0 "xaxis x+"'
c      write(ifhi,'(a,a)')
c     * 'text 0 0.1 "yaxis  G(x+,x-)*[F]?'
c     * ,'(1-x+,1-x-)/x^[a]?remn!!"'
c      if (xpar8.eq.1.) then
c        write(ifhi,'(a,e7.2,a)')'text 0.1 0.9 "s=',engy**2,' GeV^2!"'
c        write(ifhi,'(a,f5.2,a)')'text 0.1 0.8 "b=',biniDf,' fm"'
c        write(ifhi,'(a,f5.2,a)')'text 0.1 0.7 "x-=',xm,'"'
c      endif
c      write(ifhi,'(a)')       'array 2'
c
c      do i=0,nptg
c        x=xminr
c        if(xpar5.ne.0.)then
c          if (i.ne.0) x=x*(xmaxDf/xminr)**(dble(i)/dble(nptg))
c        else
c          x=xminr+(xmaxDf-xminr)*(dble(i)/dble(nptg))
c        endif
cc        z=1.d0-dsqrt(x)
cc        w(i)=Phiexpo(0.,0.,1.,z,z,engy**2,biniDf)
c        w(i)=Phiexpo(0.,0.,1.,1.d0-x,1.d0-xm,engy**2,biniDf)
c        y=u(i)*w(i)
c        if(xpar5.ne.0.)y=dabs(y)
c        write(ifhi,*)x,y
c      enddo
c
c      write(ifhi,'(a)')    '  endarray'
c      if(xpar5.ne.0.)then
c      write(ifhi,'(a)')    'closehisto plot 0-'
c
cc************************* green = omGam ************************
c
c      write(ifhi,'(a)')       'openhisto name GM'
c      write(ifhi,'(a)')       'htyp lgo'
c      write(ifhi,'(a)')       'array 2'
c
c      do i=0,nptg
c        x=xminr
c        if (i.ne.0) x=x*(xmaxDf/xminr)**(dble(i)/dble(nptg))
c        write(ifhi,*)x,dabs(u(i))
c      enddo
c
c      write(ifhi,'(a)')    '  endarray'
c      write(ifhi,'(a)')    'closehisto plot 0-'
c
cc************************* circle = PhiMExact  ************************
c
c      write(ifhi,'(a)')       'openhisto name PhiExact'
c      write(ifhi,'(a)')       'htyp poc'
c      write(ifhi,'(a)')       'array 2'
c
c      do i=0,nptg
c        x=xminr
c        if(xpar5.ne.0.)then
c          if (i.ne.0) x=x*(xmaxDf/xminr)**(dble(i)/dble(nptg))
c        else
c          x=xminr+(xmaxDf-xminr)*(dble(i)/dble(nptg))
c        endif
c        y=v(i)
c        if(xpar5.ne.0.)y=dabs(y)
c        write(ifhi,*)x,y
c      enddo
c
c      write(ifhi,'(a)')    '  endarray'
c      write(ifhi,'(a)')    'closehisto plot 0-'
c
cc************************* triangle = PhiMExpo ************************
c
c      write(ifhi,'(a)')       'openhisto name PhiExpo'
c      write(ifhi,'(a)')       'htyp pot'
c      write(ifhi,'(a)')       'array 2'
c
c      do i=0,nptg
c        x=xminr
c        if(xpar5.ne.0.)then
c          if (i.ne.0) x=x*(xmaxDf/xminr)**(dble(i)/dble(nptg))
c        else
c          x=xminr+(xmaxDf-xminr)*(dble(i)/dble(nptg))
c        endif
c        y=w(i)
c        if(xpar5.ne.0.)y=dabs(y)
c        write(ifhi,*)x,y
c      enddo
c
c      write(ifhi,'(a)')    '  endarray'
c      endif
c      write(ifhi,'(a)')    'closehisto plot 0'
c
c      end
c

c----------------------------------------------------------------------
      subroutine xParPhi
c----------------------------------------------------------------------

#include "aaa.h"
#include "ems.h"
#include "sem.h"
#include "par.h"
      double precision plc,s
      common/cems5/plc,s
      common/ckopjtg/kopj(mamx),kotg(mamx)
      common /cncl/xproj(mamx),yproj(mamx),zproj(mamx)
     *            ,xtarg(mamx),ytarg(mamx),ztarg(mamx)
      double precision x,xminr,y,z(0:200)!,Zn,Znorm
      double precision PhiExpo,PhiExact,PhiExpoK


      nptg=10                  !number of point for the graphs
      biniDf=xpar2              !value of biniDf (impact parameter)
      xminr=max(1.d-6,1.d0/dble(engy**2))  !value of xminr for plotting the function
      zz=0

c********************** full-red = PhiExact ***************************

      if (biniDf.lt.1.) then
        k=int(10.*biniDf)
        write(ifhi,'(a,I1)')  'openhisto name PhiExact-b0.',k
      else
        write(ifhi,'(a,f3.1)')  'openhisto name PhiExact-b',biniDf
      endif
      write(ifhi,'(a)')       'htyp lru'
      write(ifhi,'(a)')       'xmod lin ymod lin'
      write(ifhi,'(a,2e11.3)')'xrange',xminr,xmaxDf
      write(ifhi,'(a)')       'yrange auto auto'
      write(ifhi,'(a)')    'text 0 0 "xaxis x"'
      write(ifhi,'(a)') 'text 0 0 "yaxis  [F](1,x)/x^[a]"'
      write(ifhi,'(a,i4,a,f4.1,a)')
     * 'txt  "title E=',nint(engy),' b=',biniDf,'"'
      write(ifhi,'(a)')       'array 2'

      do i=0,nptg
        !x=xminr
        !if (i.ne.0) x=x*(xmaxDf/xminr)**(dble(i)/dble(nptg))
        x=xminr+(xmaxDf-xminr)*(dble(i)/dble(nptg))
        y=Phiexact(zz,zz,1.,dsqrt(x),dsqrt(x),smaxDf,biniDf)
     &       *dsqrt(x)**dble(-alplea(iclpro))
     &       *dsqrt(x)**dble(-alplea(icltar))
        write(ifhi,*)x,max(-1.d0,y)
      enddo

      write(ifhi,'(a)')    '  endarray'
      write(ifhi,'(a)')    'closehisto plot 0-'

c******************** blue = PhiExpo ***************************

      if (biniDf.lt.1.) then
        k=int(10.*biniDf)
        write(ifhi,'(a,I1)')  'openhisto name PhiExpo-b0.',k
      else
        write(ifhi,'(a,f3.1)')  'openhisto name PhiExpo-b',biniDf
      endif
      write(ifhi,'(a)')       'htyp lbu'
      write(ifhi,'(a)')       'xmod lin ymod lin'
      write(ifhi,'(a,2e11.3)')'xrange',xminr,xmaxDf
      write(ifhi,'(a,2e11.3)')'yrange auto auto'
      write(ifhi,'(a)')       'array 2'

      do i=0,nptg
        !x=xminr
        !if (i.ne.0) x=x*(xmaxDf/xminr)**(dble(i)/dble(nptg))
        x=xminr+(xmaxDf-xminr)*(dble(i)/dble(nptg))
        z(i)=Phiexpo(zz,zz,1.,dsqrt(x),dsqrt(x),smaxDf,biniDf)
     &       *dsqrt(x)**dble(-alplea(iclpro))
     &       *dsqrt(x)**dble(-alplea(icltar))
        write(ifhi,*) x,z(i)
      enddo

      write(ifhi,'(a)')    '  endarray'
c      write(ifhi,'(a)')    'closehisto plot 0-'
c
cc*********************yellow = PhiUnit*********************************
c
c      if (biniDf.lt.1.) then
c        k=int(10.*biniDf)
c        write(ifhi,'(a,I1)')  'openhisto name PhiUnit-b0.',k
c      else
c        write(ifhi,'(a,f3.1)')  'openhisto name PhiUnit-b',biniDf
c      endif
c      write(ifhi,'(a)')       'htyp lyu'
c      write(ifhi,'(a)')       'xmod lin ymod lin'
c      write(ifhi,'(a,2e11.3)')'xrange',xminr,xmaxDf
c      write(ifhi,'(a,2e11.3)')'yrange auto auto'
c      write(ifhi,'(a)')       'array 2'
c
c      Zn=Znorm(engy**2,biniDf)
c
c      do i=0,nptg
c        !x=xminr
c        !if (i.ne.0) x=x*(xmaxDf/xminr)**(dble(i)/dble(nptg))
c        x=xminr+(xmaxDf-xminr)*(dble(i)/dble(nptg))
c        write(ifhi,*) x,z(i)/Zn
c      enddo
c
c      write(ifhi,'(a)')    '  endarray'

      if(koll.ge.1)then

      write(ifhi,'(a)')    'closehisto plot 0-'

c**********************************************************************
c**********************************************************************
      s=smaxDf
      xproj(1)=0.
      yproj(1)=0.
      zproj(1)=0.
      lproj(1)=1
      ltarg(1)=1
      kopj(1)=1
      kotg(1)=1
      kproj(1,1)=1
      ktarg(1,1)=1
      do k=1,koll
        bk(k)=biniDf
c        bhpr(1,k)=bk(k)
        iproj(k)=1
        itarg(k)=1
      enddo
      call GfunPark(0)


c*********************green = PhiExpoK*********************************

      if (biniDf.lt.1.) then
        k=int(10.*biniDf)
        write(ifhi,'(a,I1)')  'openhisto name PhiExpok-b0.',k
      else
        write(ifhi,'(a,f3.1)')  'openhisto name PhiExpok-b',biniDf
      endif
      write(ifhi,'(a)')       'htyp lga'
      write(ifhi,'(a)')       'xmod lin ymod lin'
      write(ifhi,'(a,2e11.3)')'xrange',xminr,xmaxDf
      write(ifhi,'(a,2e11.3)')'yrange auto auto'
      write(ifhi,'(a)')       'array 2'

      do i=0,nptg
        !x=xminr
        !if (i.ne.0) x=x*(xmaxDf/xminr)**(dble(i)/dble(nptg))
        x=xminr+(xmaxDf-xminr)*(dble(i)/dble(nptg))
        z(i)=PhiExpoK(1,dsqrt(x),dsqrt(x))
        write(ifhi,*) x,z(i)
      enddo

      write(ifhi,'(a)')    '  endarray'
      endif
      write(ifhi,'(a)')    'closehisto plot 0'

c
      end

c----------------------------------------------------------------------
      subroutine xParH
c----------------------------------------------------------------------

#include "aaa.h"
#include "par.h"
      parameter(idxD2=idxD1)
      double precision GbetUni,GbetpUni,HbetUni,HbetpUni,HalpUni
      common/DGamUni/GbetUni(  idxD0:idxD2),HbetUni(  idxD0:idxD2),
     &               GbetpUni(idxD0:idxD2),HbetpUni(idxD0:idxD2),
     &               HalpUni(idxD0:idxD2)
      double precision x,xminr,y,xm,utgam2
      double precision Hrst

      nptg=20                  !number of point for the graphs
      biniDf=xpar2              !value of biniDf (impact parameter)
      xm=dble(xpar4)            !value of xminus
c.....xminr=0.d0   !value of xminr for plotting the function
      xminr=1.d0/dble(engy**2)  !value of xminr for plotting the function

      do i=idxDmin(iomega),idxDmax(iomega)
        zp=0.
        zt=0.
        call Gfunpar(zp,zt,0.,0.,1,i,biniDf,smaxDf,alpx,betx,betpx,epsp
     &               ,epst,epss,gamv)
      enddo
      imax0=idxD1
      if(iomega.ge.2)imax0=idxD
      imax1=idxD2
      if(iomega.ge.2)imax1=idxD
      do i=idxDmin(iomega),imax0
        GbetUni(i)=utgam2(betUni(i,1)+1.d0)
        GbetpUni(i)=utgam2(betpUni(i,1)+1.d0)
        HbetUni(i)=utgam2(GbetUni(i))
        HbetpUni(i)=utgam2(GbetpUni(i))
        HalpUni(i)=alpUni(i,1)
      enddo

c***********************  red = Hrst  *********************************

      write(ifhi,'(a)')       'openhisto name Hrst'
      write(ifhi,'(a)')       'htyp lru'
      write(ifhi,'(a)')       'xmod log ymod lin'
      write(ifhi,'(a,2e11.3)')'xrange',xminr,xmaxDf
      write(ifhi,'(a)')       'yrange auto auto'
      write(ifhi,'(a)')
     *     'text 0 0 "yaxis  H?2!(x+,x-)"'
      write(ifhi,'(a)')    'text 0 0 "xaxis x+"'
      write(ifhi,'(a)')
     *     'text 0 0 "yaxis  H?2!(x+,x-)"'
      if (xpar8.eq.1.) then
        write(ifhi,'(a,e7.2,a)')'text 0.1 0.2 "s=',engy**2,' GeV^2!"'
        write(ifhi,'(a,f5.2,a)')'text 0.1 0.1 "b=',biniDf,' fm"'
        write(ifhi,'(a,f5.2,a)')'text 0.1 0.3 "x-=',xm,'"'
      endif
      write(ifhi,'(a)')       'array 2'

      do i=0,nptg
        x=xminr
        if (i.ne.0) x=x*(xmaxDf/xminr)**(dble(i)/dble(nptg))
c.......x=xminr+(xmaxDf-xminr)*(dble(i)/dble(nptg))
        y=Hrst(smaxDf,biniDf,dsqrt(x),dsqrt(x))
        write(ifhi,*)x,y
      enddo

      write(ifhi,'(a)')    '  endarray'
      write(ifhi,'(a)')    'closehisto plot 0'

      end




cc----------------------------------------------------------------------
c      subroutine xParHPhiIntnew
cc----------------------------------------------------------------------
c
c#include "aaa.h"
c#include "sem.h"
c#include "ems.h"
c#include "par.h"
c      double precision x,xminr,xm,y
c      double precision PhiExact,omGam
cc      double precision PhiExpo
c
c      nptg=30                  !number of point for the graphs
c      biniDf=xpar2              !value of biniDf (impact parameter)
c      xm=dble(xpar4)            !value of xminus
cc.....xminr=0.d0   !value of xminr for plotting the function
c      xminr=1.d-3 !/dble(engy**2)  !value of xminr for plotting the function
cc************************* black = PhiExact ***************************
c
c      write(ifhi,'(a)')       'openhisto name Phi1Exact'
c      write(ifhi,'(a)')       'htyp lru'
c      if(xpar5.eq.0.)then
c        write(ifhi,'(a)')       'xmod lin ymod lin'
c      else
c        write(ifhi,'(a)')       'xmod log ymod lin'
c      endif
c      write(ifhi,'(a,2e11.3)')'xrange',xminr,xmaxDf
c      write(ifhi,'(a)')       'yrange auto auto'
c      write(ifhi,'(a)')    'text 0 0 "xaxis x+"'
c      write(ifhi,'(a)')
c     * 'text 0 0.1 "yaxis  [F]?(x+,x-)/x^[a]?remn!!"'
c      if (xpar8.eq.1.) then
c        write(ifhi,'(a,e7.2,a)')'text 0.1 0.2 "s=',engy**2,' GeV^2!"'
c        write(ifhi,'(a,f5.2,a)')'text 0.1 0.1 "b=',biniDf,' fm"'
c        write(ifhi,'(a,f5.2,a)')'text 0.1 0.3 "x-=',xm,'"'
c      endif
c      write(ifhi,'(a)')       'array 2'
c
c      do i=0,nptg
c        x=xminr
c        if (i.ne.0) x=x*(xmaxDf/xminr)**(dble(i)/dble(nptg))
cc.......x=xminr+(xmaxDf-xminr)*(dble(i)/dble(nptg))
c      y=Phiexact(0.,0.,1.,dsqrt(x),dsqrt(x),engy**2,biniDf)
c     &       *omGam(dsqrt(x),dsqrt(x),biniDf)
c        write(ifhi,*)x,y
c      enddo
c
c      write(ifhi,'(a)')    '  endarray'
c      write(ifhi,'(a)')    'closehisto plot 0'
c
c      end
c
c----------------------------------------------------------------------
      subroutine xParHPhiInt
c----------------------------------------------------------------------

#include "aaa.h"
#include "par.h"
      double precision y,HPhiInt
      common/geom/rmproj,rmtarg,bmax,bkmx
      parameter(nbib=32)


c************************ dotted = gauss integration ******************

      b1=0
      b2=max(abs(bkmx),3.)*1.2
      db=(b2-b1)/nbib

      write(ifhi,'(a)')       'openhisto name HPhiExpoInt'
      write(ifhi,'(a)')       'htyp pfc xmod lin ymod lin'
      write(ifhi,'(a)')    'text 0 0 "xaxis  impact parameter b (fm)"'
      write(ifhi,'(a)')    'text 0 0 "yaxis  Int(H[F]?pp!)(s,b)"'
      if (xpar8.eq.1.) then
        write(ifhi,'(a,e7.2,a)')'text 0.5 0.9 "s=',engy**2,' GeV^2!"'
      endif
      write(ifhi,'(a)')       'array 2'
      do k=1,nbib
        b=b1+(k-0.5)*db
        do i=idxDmin(iomega),idxDmax(iomega)
          zp=0.
          zt=0.
          call Gfunpar(zp,zt,0.,0.,1,i,b,smaxDf,alpx,betx,betpx,epsp
     &                ,epst,epss,gamv)
          zp=0.
          zt=0.
          call Gfunpar(zp,zt,0.,0.,2,i,b,smaxDf,alpx,betx,betpx,epsp
     &                ,epst,epss,gamv)
        enddo
        y=HPhiInt(smaxDf,b)
        write(ifhi,*)b,y
      enddo
      write(ifhi,'(a)')    '  endarray'
      write(ifhi,'(a)')    'closehisto plot 0'

      end

c----------------------------------------------------------------------
      subroutine xParZ
c----------------------------------------------------------------------
#include "aaa.h"
#include "par.h"
      double precision Znorm,y
      common/geom/rmproj,rmtarg,bmax,bkmx
      parameter(nbib=12)

      b1=0
      b2=max(abs(bkmx),3.)*1.2
      db=(b2-b1)/nbib

c************************full-red = Znorm *****************************

      write(ifhi,'(a)')       'openhisto name Znorm'
      write(ifhi,'(a)')       'htyp lru xmod lin ymod lin'
      write(ifhi,'(a)')    'text 0 0 "xaxis  impact parameter b (fm)"'
      write(ifhi,'(a)')    'text 0 0 "yaxis  Z(s,b) "'
      if (xpar8.eq.1.) then
        write(ifhi,'(a,e7.2,a)')'text 0.1 0.1 "s=',engy**2,' GeV^2!"'
      endif
      write(ifhi,'(a)')       'array 2'
c        y=Znorm(engy**2,xpar2)
      do k=1,nbib
        b=b1+(k-0.5)*db
        y=Znorm(smaxDf,b)
        write(ifhi,*)b,y
      enddo
      write(ifhi,'(a)')    '  endarray'
      write(ifhi,'(a)')    'closehisto plot 0-'


c**********************************************************************

      write(ifhi,'(a)')       'openhisto name un'
      write(ifhi,'(a)')       'htyp lba xmod lin ymod lin'
      write(ifhi,'(a)')       'array 2'
      do k=1,nbib
        b=b1+(k-0.5)*db
        y=1
        write(ifhi,*)b,y
      enddo
      write(ifhi,'(a)')    '  endarray'
      write(ifhi,'(a)')    'closehisto plot 0'


      end

c----------------------------------------------------------------------
      subroutine xAlphaS
c----------------------------------------------------------------------
#include "aaa.h"
#include "sem.h"
      parameter(nbib=1000)

      q1=0.1 !sqrt(qcdlam)
      q2=100000

c************************full-red = Znorm *****************************

      write(ifhi,'(a)')       'openhisto name pssalf'
      write(ifhi,'(a)')       'htyp lru xmod log ymod lin'
      write(ifhi,*)       'xrange ',q1,q2,' yrange 0.05 1.1'
      write(ifhi,'(a)')    'text 0 0 "xaxis  Q (GeV)"'
      write(ifhi,'(a)')    'text 0 0 "yaxis  [a]?s!(Q) "'
      write(ifhi,'(a)')       'array 2'
      do k=0,nbib
        q=q1
        if (k.ne.0) q=q*(q2/q1)**(float(k)/float(nbib))
        qq=q*q
        y=pssalf(qq/qcdlam)*2.*pi
        write(ifhi,*)q,y
      enddo
      write(ifhi,'(a)')    '  endarray'
      write(ifhi,'(a)')    'closehisto plot 0-'
     

c**********************************************************************

      write(ifhi,'(a)')       'openhisto name alphaS'
      write(ifhi,'(a)')       'htyp pbq xmod log ymod lin'
      write(ifhi,'(a)')       'array 2'
      qcdl=0.038                 !0.033 !0.04 !0.011 !0.03
      xkap=0.9/(log(log(0.9/qcdl))-log(4.*pi/9.))
      print *,'kappa',sqrt(0.25*xkap)
      do k=1,nbib
        q=q1
        if (k.ne.0) q=q*(q2/q1)**(float(k)/float(nbib))
        if(q.lt.sqrt(0.9))then
          qq=q*q/qcdl
          y=exp(-qq/xkap*qcdl)/2./pi
        else
          qlam=qcdl
          if(q.lt.2.*qcmass) then
            nf=3
          elseif(q.lt.qbmass)then
            nf=4
            qlam=qcdl*0.68
          else
            nf=5
            qlam=qcdl*0.38
          endif
          qq=q*q/qlam
          y=2./(11.-nf/1.5)/log(qq)
        endif
c        if(y.lt.0.or.y.gt.0.5/(2.*pi))then
c          y=0.5/(2.*pi)
c        endif
        y=y*2.*pi
        write(ifhi,*)q,y
      enddo
      write(ifhi,'(a)')    '  endarray'
      write(ifhi,'(a)')    'closehisto plot 0-'


      end

c----------------------------------------------------------------------
      subroutine xParPro
c----------------------------------------------------------------------
#include "aaa.h"
#include "sem.h"
#include "par.h"
      parameter(nbib=12)
      double precision PhiExact,PhiExpo,y(nbib,9),om1intb,om1intbc
     &,om1intgc,om1intbi!,PhiUnit
      common/geom/rmproj,rmtarg,bmax,bkmx

      b1=0
      b2=max(abs(bkmx),3.)
      db=(b2-b1)/nbib
      zz=0

c********************* full-red = 1-PhiExact **************************

      write(ifhi,'(a)')       'openhisto name 1-PhiExact'
      write(ifhi,'(a)')       'htyp lru xmod lin ymod lin'
      write(ifhi,'(a,2e11.3)')'xrange',0.,b2
      write(ifhi,'(a)')    'text 0 0 "xaxis  impact parameter b (fm)"'
      write(ifhi,'(a)')    'text 0 0 "yaxis  1-[F]?pp!(1,1,1) "'
      if (xpar8.eq.1.) then
        write(ifhi,'(a,e7.2,a)')'text 0.5 0.9 "s=',engy**2,' GeV^2!"'
      endif
      write(ifhi,'(a)')       'array 2'
      do k=1,nbib
        b=b1+(k-0.5)*db
        y(k,1)=min(2d0,1.d0-Phiexact(zz,zz,1.,1.d0,1.d0,smaxDf,b))
        y(k,2)=1.d0-Phiexpo(zz,zz,1.,1.d0,1.d0,smaxDf,b)
        y(k,3)=1.d0
        y(k,4)=om1intbc(b)
        y(k,5)=om1intb(b)
        y(k,6)=om1intgc(b)
        y(k,7)=om1intbi(b,0)
        y(k,8)=om1intbi(b,1)
        y(k,9)=om1intbi(b,2)
        write(ifhi,*)b,y(k,1)
      enddo
      write(ifhi,'(a)')    '  endarray'
      write(ifhi,'(a)')    'closehisto plot 0-'

c************************** blue-dashed = 1-PhiExpo *******************

      write(ifhi,'(a)')       'openhisto name 1-PhiExpo'
      write(ifhi,'(a)')       'htyp lba xmod lin ymod lin'
      write(ifhi,'(a)')    'text 0 0 "xaxis  impact parameter b (fm)"'
      write(ifhi,'(a)')    'text 0 0 "yaxis  1-[F]?pp!(1,1,1) "'
      if (xpar8.eq.1.) then
        write(ifhi,'(a,e7.2,a)')'text 0.5 0.9 "s=',engy**2,' GeV^2!"'
      endif
      write(ifhi,'(a)')       'array 2'
      do k=1,nbib
        b=b1+(k-0.5)*db
        write(ifhi,*)b,y(k,2)
      enddo
      write(ifhi,'(a)')    '  endarray'
      write(ifhi,'(a)')    'closehisto plot 0-'


c**********************************************************************

      write(ifhi,'(a)')       'openhisto'
      write(ifhi,'(a)')       'htyp lga xmod lin ymod lin'
      write(ifhi,'(a)')       'array 2'
      do k=1,nbib
        b=b1+(k-0.5)*db
        write(ifhi,*)b,y(k,3)
      enddo
      write(ifhi,'(a)')    '  endarray'
      write(ifhi,'(a)')    'closehisto plot 0'

c****************************** red = om1intbc ********************

      write(ifhi,'(a)')       'openhisto name om1intbc'
      write(ifhi,'(a)')       'htyp lru xmod lin ymod lin'
      write(ifhi,'(a,2e11.3)')'xrange',0.,b2
      write(ifhi,'(a)')    'text 0 0 "xaxis  impact parameter b (fm)"'
      write(ifhi,'(a)')    'text 0 0 "yaxis  [w]?1bc!(b) "'
      if (xpar8.eq.1.) then
        write(ifhi,'(a,e7.2,a)')'text 0.5 0.9 "s=',engy**2,' GeV^2!"'
      endif
      write(ifhi,'(a)')       'array 2'

      do k=1,nbib
        b=b1+(k-0.5)*db
        write(ifhi,*)b,y(k,4)
      enddo
      write(ifhi,'(a)')    '  endarray'
      write(ifhi,'(a)')    'closehisto plot 0-'


c************************* blue dashed =  om1intb ********************

      write(ifhi,'(a)')       'openhisto name om1intb'
      write(ifhi,'(a)')       'htyp lba xmod lin ymod lin'
      write(ifhi,'(a)')    'text 0 0 "xaxis  impact parameter b (fm)"'
      write(ifhi,'(a)')    'text 0 0 "yaxis  [w]?1b!(b) "'
      if (xpar8.eq.1.) then
        write(ifhi,'(a,e7.2,a)')'text 0.5 0.9 "s=',engy**2,' GeV^2!"'
      endif
      write(ifhi,'(a)')       'array 2'

      do k=1,nbib
        b=b1+(k-0.5)*db
        write(ifhi,*)b,y(k,5)
      enddo
      write(ifhi,'(a)')    '  endarray'
      write(ifhi,'(a)')    'closehisto plot 0-'

c****************************** green dot = om1intgc ********************

      write(ifhi,'(a)')       'openhisto name om1intgc'
      write(ifhi,'(a)')       'htyp lgo xmod lin ymod lin'
      write(ifhi,'(a)')    'text 0 0 "xaxis  impact parameter b (fm)"'
      write(ifhi,'(a)')    'text 0 0 "yaxis  [w]?1gc!(b) "'
      if (xpar8.eq.1.) then
        write(ifhi,'(a,e7.2,a)')'text 0.5 0.9 "s=',engy**2,' GeV^2!"'
      endif
      write(ifhi,'(a)')       'array 2'

      do k=1,nbib
        b=b1+(k-0.5)*db
        write(ifhi,*)b,y(k,6)
      enddo
      write(ifhi,'(a)')    '  endarray'
      write(ifhi,'(a)')    'closehisto plot 0'

c****************************** red = om1intbi(0) ********************

      write(ifhi,'(a)')       'openhisto name om1intbi0'
      write(ifhi,'(a)')       'htyp lru xmod lin ymod log'
      write(ifhi,'(a,2e11.3)')'xrange',0.,b2
      write(ifhi,'(a)')    'text 0 0 "xaxis  impact parameter b (fm)"'
      write(ifhi,'(a)')    'text 0 0 "yaxis  [w]?1bc!(b) "'
      if (xpar8.eq.1.) then
        write(ifhi,'(a,e7.2,a)')'text 0.5 0.9 "s=',engy**2,' GeV^2!"'
      endif
      write(ifhi,'(a)')       'array 2'

      do k=1,nbib
        b=b1+(k-0.5)*db
        if(y(k,7).gt.1.e-10)write(ifhi,*)b,y(k,7)
      enddo
      write(ifhi,'(a)')    '  endarray'
      write(ifhi,'(a)')    'closehisto plot 0-'


c************************* blue dashed =  om1intbi(1) ********************

      write(ifhi,'(a)')       'openhisto name om1intbi1'
      write(ifhi,'(a)')       'htyp lba xmod lin ymod lin'
      write(ifhi,'(a)')    'text 0 0 "xaxis  impact parameter b (fm)"'
      write(ifhi,'(a)')    'text 0 0 "yaxis  [w]?1b!(b) "'
      if (xpar8.eq.1.) then
        write(ifhi,'(a,e7.2,a)')'text 0.5 0.9 "s=',engy**2,' GeV^2!"'
      endif
      write(ifhi,'(a)')       'array 2'

      do k=1,nbib
        b=b1+(k-0.5)*db
        if(y(k,8).gt.1.e-10)write(ifhi,*)b,y(k,8)
      enddo
      write(ifhi,'(a)')    '  endarray'
      write(ifhi,'(a)')    'closehisto plot 0-'

c****************************** green dot = om1intbi(2) ********************

      write(ifhi,'(a)')       'openhisto name om1intbi2'
      write(ifhi,'(a)')       'htyp lgo xmod lin ymod lin'
      write(ifhi,'(a)')    'text 0 0 "xaxis  impact parameter b (fm)"'
      write(ifhi,'(a)')    'text 0 0 "yaxis  [w]?1gc!(b) "'
      if (xpar8.eq.1.) then
        write(ifhi,'(a,e7.2,a)')'text 0.5 0.9 "s=',engy**2,' GeV^2!"'
      endif
      write(ifhi,'(a)')       'array 2'

      do k=1,nbib
        b=b1+(k-0.5)*db
        if(y(k,9).gt.1.e-10)write(ifhi,*)b,y(k,9)
      enddo
      write(ifhi,'(a)')    '  endarray'

      if(iscreen.ne.0)then
 
      write(ifhi,'(a)')    'closehisto plot 0-'

c****************************** red = om1intbi(0) ********************

      write(ifhi,'(a)')       'openhisto name om1intbi00'
      write(ifhi,'(a)')       'htyp lrv xmod lin ymod log'
      write(ifhi,'(a)')    'text 0 0 "xaxis  impact parameter b (fm)"'
      write(ifhi,'(a)')    'text 0 0 "yaxis  [w]?1bc!(b) "'
      if (xpar8.eq.1.) then
        write(ifhi,'(a,e7.2,a)')'text 0.5 0.9 "s=',engy**2,' GeV^2!"'
      endif
      write(ifhi,'(a)')       'array 2'

      do k=1,nbib
        b=b1+(k-0.5)*db
        tmp=om1intbi(b,0)
        if(tmp.gt.1.e-10)write(ifhi,*)b,tmp
      enddo
      write(ifhi,'(a)')    '  endarray'
      write(ifhi,'(a)')    'closehisto plot 0-'


c************************* blue dashed =  om1intbi(1) ********************

      write(ifhi,'(a)')       'openhisto name om1intbi10'
      write(ifhi,'(a)')       'htyp lbb xmod lin ymod lin'
      write(ifhi,'(a)')    'text 0 0 "xaxis  impact parameter b (fm)"'
      write(ifhi,'(a)')    'text 0 0 "yaxis  [w]?1b!(b) "'
      if (xpar8.eq.1.) then
        write(ifhi,'(a,e7.2,a)')'text 0.5 0.9 "s=',engy**2,' GeV^2!"'
      endif
      write(ifhi,'(a)')       'array 2'

      do k=1,nbib
        b=b1+(k-0.5)*db
        tmp=om1intbi(b,1)
        if(tmp.gt.1.e-10)write(ifhi,*)b,tmp
      enddo
      write(ifhi,'(a)')    '  endarray'
      write(ifhi,'(a)')    'closehisto plot 0-'

c****************************** green dot = om1intbi(2) ********************

      write(ifhi,'(a)')       'openhisto name om1intbi20'
      write(ifhi,'(a)')       'htyp lgp xmod lin ymod lin'
      write(ifhi,'(a)')    'text 0 0 "xaxis  impact parameter b (fm)"'
      write(ifhi,'(a)')    'text 0 0 "yaxis  [w]?1gc!(b) "'
      if (xpar8.eq.1.) then
        write(ifhi,'(a,e7.2,a)')'text 0.5 0.9 "s=',engy**2,' GeV^2!"'
      endif
      write(ifhi,'(a)')       'array 2'

      do k=1,nbib
        b=b1+(k-0.5)*db
        tmp=om1intbi(b,2)
        if(tmp.gt.1.e-10)write(ifhi,*)b,tmp
      enddo
      write(ifhi,'(a)')    '  endarray'
      endif

      write(ifhi,'(a)')    'closehisto plot 0'



      end

c----------------------------------------------------------------------
      subroutine xParPro1
c----------------------------------------------------------------------
#include "aaa.h"
#include "sem.h"
#include "par.h"
      double precision PhiExact,PhiExpo,y
      common/geom/rmproj,rmtarg,bmax,bkmx
      parameter(nbib=12)

      b1=0
      b2=max(abs(bkmx),3.)*1.2
      db=(b2-b1)/nbib
      zz=0

c********************* full-red = 1-PhiExact **************************

      write(ifhi,'(a)')       'openhisto name 1-PhiExact'
      write(ifhi,'(a)')       'htyp lru xmod lin ymod lin'
      write(ifhi,'(a)')    'text 0 0 "xaxis  impact parameter b (fm)"'
      write(ifhi,'(a)')    'text 0 0 "yaxis  1-[F]?pp!(0.5,1,1) "'
      if (xpar8.eq.1.) then
        write(ifhi,'(a,e7.2,a)')'text 0.5 0.9 "s=',engy**2,' GeV^2!"'
      endif
      write(ifhi,'(a)')       'array 2'
      do k=1,nbib
        b=b1+(k-0.5)*db
        y=min(2.d0,1.d0-Phiexact(zz,zz,.5,1.d0,1.d0,smaxDf,b))
        write(ifhi,*)b,y
      enddo
      write(ifhi,'(a)')    '  endarray'
      write(ifhi,'(a)')    'closehisto plot 0-'

c************************** blue-dashed = 1-PhiExpo *******************

      write(ifhi,'(a)')       'openhisto name 1-PhiExpo'
      write(ifhi,'(a)')       'htyp lba xmod lin ymod lin'
      write(ifhi,'(a)')    'text 0 0 "xaxis  impact parameter b (fm)"'
      write(ifhi,'(a)')    'text 0 0 "yaxis  1-[F]?pp!(0.5,1,1) "'
      if (xpar8.eq.1.) then
        write(ifhi,'(a,e7.2,a)')'text 0.5 0.9 "s=',engy**2,' GeV^2!"'
      endif
      write(ifhi,'(a)')       'array 2'
      do k=1,nbib
        b=b1+(k-0.5)*db
        y=1.d0-Phiexpo(zz,zz,.5,1.d0,1.d0,smaxDf,b)
        write(ifhi,*)b,y
      enddo
      write(ifhi,'(a)')    '  endarray'
      write(ifhi,'(a)')    'closehisto plot 0'




      end

c----------------------------------------------------------------------
      subroutine xParGam
c----------------------------------------------------------------------

#include "aaa.h"
#include "sem.h"
#include "par.h"
      dimension bet(idxD0:idxD1)
      double precision utgam2,xgammag2!,xrem
      dimension ip(idxD0:idxD1),imax(idxD0:idxD1)

      nptg=50                  !number of point for the graphs
      b=xpar2
      gamp=xpar6
      zmax=6.
c      xrem=dble(xpar4)
      if(idxD0.ne.0.or.idxD1.ne.2) stop "Check xPargam"

      do i=idxD0,idxD1
        imax(i)=4
        bet(i)=0
      enddo
      nmax=idxD1
      imax(idxD0)=int(zmax)
      imax(1)=imax(idxD0)
      imax(2)=imax(idxD0)

      do i=idxD0,nmax
        gam=gamD(i,iclpro,icltar)*b**2
        bet(i)=gam+betDp(i,iclpro,icltar)-alppar+1.
      enddo
      write(ifhi,'(a)')       'openhisto name gExact'
      write(ifhi,'(a)')       'htyp pfs'
      write(ifhi,'(a)')       'xmod lin ymod log'
      write(ifhi,'(a,2e11.3)')'xrange',0.,zmax
      write(ifhi,'(a)')       'yrange auto auto'
c      write(ifhi,'(a)')'yrange 1.e-10 2'
      write(ifhi,'(a)')    'text 0 0 "xaxis z"'
      write(ifhi,'(a)')    'text 0 0 "yaxis g(z) "'
      write(ifhi,'(a)')       'array 2'

      do ip0=0,imax(0)
        ip(0)=ip0
        do ip1=0,imax(1)
          ip(1)=ip1
          do ip2=0,imax(2)
            ip(2)=ip2
          t=0.
          do i=idxD0,nmax
            t=t+float(ip(i))*bet(i)
          enddo
          write(ifhi,'(2e14.6)')t,utgam2(dble(alplea(2))+1.D0)
     &       /utgam2(dble(alplea(2))+1.D0+dble(t))
        enddo
      enddo
      enddo

      write(ifhi,'(a)')    '  endarray'
      write(ifhi,'(a)')    'closehisto plot 0-'


c**********************************************************************

      write(ifhi,'(a)')       'openhisto name gExpo'
      write(ifhi,'(a)')       'htyp lbu'
      write(ifhi,'(a)')       'xmod lin ymod log'
      write(ifhi,'(a,2e11.3)')'xrange',0.,zmax
      write(ifhi,'(a)')'yrange auto auto'
      write(ifhi,'(a)')    'text 0 0 "xaxis z"'
      write(ifhi,'(a)')    'text 0 0 "yaxis g(z)"'
      write(ifhi,'(a)')       'array 2'

      do i=0,nptg
        t=zmax*(float(i)/float(nptg))
        write(ifhi,'(2e14.6)') t,dexp(-dble(t))
      enddo

      write(ifhi,'(a)')    '  endarray'
      write(ifhi,'(a)')    'closehisto plot 0-'


c**********************************************************************



      write(ifhi,'(a)')       'openhisto name gPower'
      write(ifhi,'(a)')       'htyp poc'
      write(ifhi,'(a)')       'xmod log ymod log'
      write(ifhi,'(a,2e11.3)')'xrange',0.,zmax
      write(ifhi,'(a)')       'yrange auto auto'
      write(ifhi,'(a)')    'text 0 0 "xaxis t"'
      write(ifhi,'(2a)')    'text 0 0 ',
     &     '"yaxis [P][G](1+[a]?L!)/[G](1+[a]?L!+[b])"'
      write(ifhi,'(a)')       'array 2'

      do ip0=0,imax(0)
        ip(0)=ip0
        do ip1=0,imax(1)
          ip(1)=ip1
          do ip2=0,imax(2)
            ip(2)=ip2
c$$$            do ip3=0,imax(3)
c$$$              ip(3)=ip3
              t=0.
              do i=idxD0,nmax
                t=t+float(ip(i))*bet(i)
              enddo
              write(ifhi,'(2e14.6)')t,xgammag2(iclpro,bet,ip,gamp)
c$$$            enddo
          enddo
        enddo
      enddo

      write(ifhi,'(a)')    '  endarray'
      write(ifhi,'(a)')    'closehisto plot 0'

      end

c----------------------------------------------------------------------
      subroutine xParOmega1xy
c----------------------------------------------------------------------
c xpar2=b
c xpar3=y
c xpar4=xh
c xpar99: nucl coef
c xpar6 : 0=xp/xm from om1xpork
c         1=xh/yp
c         2=xp/xm from om1xprk
c----------------------------------------------------------------------
#include "ems.h"
#include "sem.h"
#include "aaa.h"
#include "par.h"
      double precision plc,s
      common/cems5/plc,s
      common/ckopjtg/kopj(mamx),kotg(mamx)
      common /cncl/xproj(mamx),yproj(mamx),zproj(mamx)
     *            ,xtarg(mamx),ytarg(mamx),ztarg(mamx)
      double precision x,ranhis(0:51),y,ymax,ymin,xh,om1xo,om1xoI
     &,om1xpk,om1xmk,t,om1xk,om1yk,xpr1,xmr1,om1x,om1xI
      external om1x,om1xI,om1xo,om1xoI
      common /psar7/ delx,alam3p,gam3p

      nptg=50                  !number of point for the graphs
      biniDf=xpar2              !value of biniDf (impact paramter)
      xh=xminDf
      xpr1=1d0
      xmr1=1d0
      if(xpar4.lt.1..and.xpar4.gt.0.)then
        xh=dble(xpar4)          !value of x
        xpr1=dble(xpar4**0.25)
        xmr1=dble(xpar4**0.5)
      endif
      do i=0,51
        ranhis(i)=0.d0
      enddo
c$$$      xp=dsqrt(xh)*dble(exp(xpar3)) !y=xpar3
c$$$      xm=1.d0
c$$$      if(xp.ne.0.d0)xm=xh/xp

      s=smaxDf
      xproj(1)=0.
      yproj(1)=0.
      zproj(1)=0.
      lproj(1)=1
      ltarg(1)=1
      kopj(1)=1
      kotg(1)=1
      kproj(1,1)=1
      ktarg(1,1)=1
      do k=1,koll
        bk(k)=biniDf
c        bhpr(1,k)=bk(k)
        iproj(k)=1
        itarg(k)=1
      enddo
      call GfunPark(0)
      !call GfunParP


      if(nint(xpar6).eq.0)then

      call xhistomran1(ranhis,xpr1,xmr1)

      write(ifhi,'(a)')       'openhisto name Om1xpo'
      write(ifhi,'(a)')       'htyp lin'
      write(ifhi,'(a)')       'xmod log ymod log'
      write(ifhi,'(a,2e11.3)')'xrange',xminDf,xmaxDf
      write(ifhi,'(a)')'yrange auto auto'
      write(ifhi,'(a)')    'text 0 0 "xaxis X-"'
      write(ifhi,'(a)')    'text 0 0 "yaxis [w]?1!(x+)"'
      if (xpar8.eq.1.) then
        write(ifhi,'(a,e7.2,a)')  'text 0.1 0.9 "s=',engy**2,' GeV^2!"'
        write(ifhi,'(a,f5.2,a)')  'text 0.5 0.9 "b=',biniDf,' fm"'
      endif
      write(ifhi,'(a)')       'array 2'

      do i=0,nptg
        x=xminDf
        if (i.ne.0) x=x*(xmaxDf/xminDf)**(dble(i)/dble(nptg))
c.......x=xminDf+(xmaxDf-xminDf)*(dble(i)/dble(nptg))
        write(ifhi,*) x,om1xpk(om1xo,om1xoI,x,xpr1,xmr1,1)
      enddo

      write(ifhi,'(a)')    '  endarray'
      write(ifhi,'(a)')    'closehisto plot 0-'

      write(ifhi,'(a)')       'openhisto name Om1xpoRan'
      write(ifhi,'(a)')       'htyp pnt'
      write(ifhi,'(a)')       'xmod log ymod log'
      write(ifhi,'(a,2e11.3)')'xrange',xminDf,xmaxDf
      write(ifhi,'(a)')'yrange auto auto'
      write(ifhi,'(a)')    'text 0 0 "xaxis X+"'
      write(ifhi,'(a)')    'text 0 0 "yaxis [w]?1!(x+) random"'
      if (xpar8.eq.1.) then
        write(ifhi,'(a,e7.2,a)')  'text 0.1 0.9 "s=',engy**2,' GeV^2!"'
        write(ifhi,'(a,f5.2,a)')  'text 0.5 0.9 "b=',biniDf,' fm"'
      endif
      write(ifhi,'(a)')       'array 2'

      do i=0,51
        x=xminDf
        if (i.ne.0) x=x*(xmaxDf/xminDf)**((dble(i)+.5d0)/51.d0)
c.......x=xminDf+(xmaxDf-xminDf)*((dble(i)+.5d0)/51.d0)
        write(ifhi,*) x,ranhis(i)
      enddo

      write(ifhi,'(a)')    '  endarray'
      write(ifhi,'(a)')    'closehisto plot 0'

      call xhistomran2(ranhis,xh,xpr1,xmr1)

      write(ifhi,'(a)')       'openhisto name Om1xmo'
      write(ifhi,'(a)')       'htyp lin'
      write(ifhi,'(a)')       'xmod log ymod log'
      write(ifhi,'(a,2e11.3)')'xrange',xminDf,xmaxDf
      write(ifhi,'(a)')'yrange auto auto'
      write(ifhi,'(a)')    'text 0 0 "xaxis X-"'
      write(ifhi,'(a)')    'text 0 0 "yaxis [w]?1!(x-)"'
      if (xpar8.eq.1.) then
        write(ifhi,'(a,e7.2,a)')  'text 0.1 0.9 "s=',engy**2,' GeV^2!"'
        write(ifhi,'(a,f5.2,a)')  'text 0.5 0.9 "b=',biniDf,' fm"'
      endif
      write(ifhi,'(a)')       'array 2'

      do i=0,nptg
        x=xminDf
        if (i.ne.0) x=x*(xmaxDf/xminDf)**(dble(i)/dble(nptg))
c.......x=xminDf+(xmaxDf-xminDf)*(dble(i)/dble(nptg))
        write(ifhi,*) x,om1xmk(om1xo,om1xoI,xh,x,xpr1,xmr1,1)
      enddo

      write(ifhi,'(a)')    '  endarray'
      write(ifhi,'(a)')    'closehisto plot 0-'

      write(ifhi,'(a)')       'openhisto name Om1xmoRan'
      write(ifhi,'(a)')       'htyp pnt'
      write(ifhi,'(a)')       'xmod log ymod log'
      write(ifhi,'(a,2e11.3)')'xrange',xminDf,xmaxDf
      write(ifhi,'(a)')'yrange auto auto'
      write(ifhi,'(a)')    'text 0 0 "xaxis X-"'
      write(ifhi,'(a)')    'text 0 0 "yaxis [w]?1!(x-) random"'
      if (xpar8.eq.1.) then
        write(ifhi,'(a,e7.2,a)')  'text 0.1 0.9 "s=',engy**2,' GeV^2!"'
        write(ifhi,'(a,f5.2,a)')  'text 0.5 0.9 "b=',biniDf,' fm"'
      endif
      write(ifhi,'(a)')       'array 2'

      do i=0,51
        x=xminDf
        if (i.ne.0) x=x*(xmaxDf/xminDf)**((dble(i)+.5d0)/51.d0)
c.......x=xminDf+(xmaxDf-xminDf)*((dble(i)+.5d0)/51.d0)
        write(ifhi,*) x,ranhis(i)
      enddo

      write(ifhi,'(a)')    '  endarray'
      write(ifhi,'(a)')    'closehisto plot 0'


c**********************************************************************

      elseif(nint(xpar6).eq.2)then

        call xhistomran8(ranhis,xpr1,xmr1)


      write(ifhi,'(a)')       'openhisto name Om1xp'
      write(ifhi,'(a)')       'htyp lin'
      write(ifhi,'(a)')       'xmod log ymod log'
      write(ifhi,'(a,2e11.3)')'xrange',xminDf,xmaxDf
      write(ifhi,'(a)')'yrange auto auto'
      write(ifhi,'(a)')    'text 0 0 "xaxis X+"'
      write(ifhi,'(a)')    'text 0 0 "yaxis [w]?1!(x+)"'
      if (xpar8.eq.1.) then
        write(ifhi,'(a,e7.2,a)')  'text 0.1 0.9 "s=',engy**2,' GeV^2!"'
        write(ifhi,'(a,f5.2,a)')  'text 0.5 0.9 "b=',biniDf,' fm"'
      endif
      write(ifhi,'(a)')       'array 2'

      do i=0,nptg
        x=xminDf
        if (i.ne.0) x=x*(xmaxDf/xminDf)**(dble(i)/dble(nptg))
c.......x=xminDf+(xmaxDf-xminDf)*(dble(i)/dble(nptg))
        t=om1xpk(om1x,om1xI,x,xpr1,xmr1,1)
        write(ifhi,*) x,t
      enddo

      write(ifhi,'(a)')    '  endarray'
      write(ifhi,'(a)')    'closehisto plot 0-'

      write(ifhi,'(a)')       'openhisto name Om1xpRan'
      write(ifhi,'(a)')       'htyp pnt'
      write(ifhi,'(a)')       'xmod log ymod log'
      write(ifhi,'(a,2e11.3)')'xrange',xminDf,xmaxDf
      write(ifhi,'(a)')'yrange auto auto'
      write(ifhi,'(a)')    'text 0 0 "xaxis X+"'
      write(ifhi,'(a)')    'text 0 0 "yaxis [w]?1!(x+) random"'
      if (xpar8.eq.1.) then
        write(ifhi,'(a,e7.2,a)')  'text 0.1 0.9 "s=',engy**2,' GeV^2!"'
        write(ifhi,'(a,f5.2,a)')  'text 0.5 0.9 "b=',biniDf,' fm"'
      endif
      write(ifhi,'(a)')       'array 2'

      do i=0,51
        x=xminDf
        if (i.ne.0) x=x*(xmaxDf/xminDf)**(dble(i)/dble(nptg))
c.......x=xminDf+(xmaxDf-xminDf)*((dble(i)+.5d0)/51.d0)
        write(ifhi,*) x,ranhis(i)
      enddo

      write(ifhi,'(a)')    '  endarray'
      write(ifhi,'(a)')    'closehisto plot 0'

      do i=0,51
        ranhis(i)=0.d0
      enddo

      call xhistomran9(ranhis,xh,xpr1,xmr1)


      write(ifhi,'(a)')       'openhisto name Om1xm'
      write(ifhi,'(a)')       'htyp lin'
      write(ifhi,'(a)')       'xmod log ymod log'
      write(ifhi,'(a,2e11.3)')'xrange',xminDf,xmaxDf
      write(ifhi,'(a)')'yrange auto auto'
      write(ifhi,'(a)')    'text 0 0 "xaxis X-"'
      write(ifhi,'(a)')    'text 0 0 "yaxis [w]?1!(x-)"'
      if (xpar8.eq.1.) then
        write(ifhi,'(a,e7.2,a)')  'text 0.1 0.9 "s=',engy**2,' GeV^2!"'
        write(ifhi,'(a,f5.2,a)')  'text 0.5 0.9 "b=',biniDf,' fm"'
        write(ifhi,'(a,f5.2,a)')  'text 0.5 0.8 "x+=',xh,'"'
      endif
      write(ifhi,'(a)')       'array 2'

      do i=0,nptg
        x=xminDf
        if (i.ne.0) x=x*(xmaxDf/xminDf)**(dble(i)/dble(nptg))
c.......x=xminDf+(xmaxDf-xminDf)*(dble(i)/dble(nptg))
        t=om1xmk(om1x,om1xI,xh,x,xpr1,xmr1,1)
        write(ifhi,*) x,t
      enddo

      write(ifhi,'(a)')    '  endarray'
      write(ifhi,'(a)')    'closehisto plot 0-'

      write(ifhi,'(a)')       'openhisto name Om1xmRan'
      write(ifhi,'(a)')       'htyp pnt'
      write(ifhi,'(a)')       'xmod log ymod log'
      write(ifhi,'(a,2e11.3)')'xrange',xminDf,xmaxDf
      write(ifhi,'(a)')'yrange auto auto'
      write(ifhi,'(a)')    'text 0 0 "xaxis X-"'
      write(ifhi,'(a)')    'text 0 0 "yaxis [w]?1!(x-) random"'
      if (xpar8.eq.1.) then
        write(ifhi,'(a,e7.2,a)')  'text 0.1 0.9 "s=',engy**2,' GeV^2!"'
        write(ifhi,'(a,f5.2,a)')  'text 0.5 0.9 "b=',biniDf,' fm"'
        write(ifhi,'(a,f5.2,a)')  'text 0.5 0.8 "x+=',xh,'"'
      endif
      write(ifhi,'(a)')       'array 2'

      do i=0,51
        x=xminDf
        if (i.ne.0) x=x*(xmaxDf/xminDf)**(dble(i)/dble(nptg))
c.......x=xminDf+(xmaxDf-xminDf)*((dble(i)+.5d0)/51.d0)
        write(ifhi,*) x,ranhis(i)
      enddo

      write(ifhi,'(a)')    '  endarray'
      write(ifhi,'(a)')    'closehisto plot 0'

c**********************************************************************
      else

      call xhistomran10(ranhis,xpr1,xmr1)


      write(ifhi,'(a)')       'openhisto name Om1xk'
      write(ifhi,'(a)')       'htyp lin'
      write(ifhi,'(a)')       'xmod log ymod log'
      write(ifhi,'(a,2e11.3)')'xrange',xminDf,xmaxDf
      write(ifhi,'(a)')'yrange auto auto'
      write(ifhi,'(a)')    'text 0 0 "xaxis X"'
      write(ifhi,'(a)')    'text 0 0 "yaxis [w]?1!(x)"'
      if (xpar8.eq.1.) then
        write(ifhi,'(a,e7.2,a)')  'text 0.1 0.9 "s=',engy**2,' GeV^2!"'
        write(ifhi,'(a,f5.2,a)')  'text 0.5 0.9 "b=',biniDf,' fm"'
      endif
      write(ifhi,'(a)')       'array 2'

      do i=0,nptg
        x=xminDf
        if (i.ne.0) x=x*(xmaxDf/xminDf)**(dble(i)/dble(nptg))
c.......x=xminDf+(xmaxDf-xminDf)*(dble(i)/dble(nptg))
        t=om1xk(x,xpr1,xmr1,1)
        write(ifhi,*) x,t
      enddo

      write(ifhi,'(a)')    '  endarray'
      write(ifhi,'(a)')    'closehisto plot 0-'

      write(ifhi,'(a)')       'openhisto name Om1xRan'
      write(ifhi,'(a)')       'htyp pnt'
      write(ifhi,'(a)')       'xmod log ymod log'
      write(ifhi,'(a,2e11.3)')'xrange',xminDf,xmaxDf
      write(ifhi,'(a)')'yrange auto auto'
      write(ifhi,'(a)')    'text 0 0 "xaxis X+"'
      write(ifhi,'(a)')    'text 0 0 "yaxis [w]?1!(x) random"'
      if (xpar8.eq.1.) then
        write(ifhi,'(a,e7.2,a)')  'text 0.1 0.9 "s=',engy**2,' GeV^2!"'
        write(ifhi,'(a,f5.2,a)')  'text 0.5 0.9 "b=',biniDf,' fm"'
      endif
      write(ifhi,'(a)')       'array 2'

      do i=0,51
        x=xminDf
        if (i.ne.0) x=x*(xmaxDf/xminDf)**(dble(i)/dble(nptg))
c.......x=xminDf+(xmaxDf-xminDf)*((dble(i)+.5d0)/51.d0)
        write(ifhi,*) x,ranhis(i)
      enddo

      write(ifhi,'(a)')    '  endarray'
      write(ifhi,'(a)')    'closehisto plot 0'

      ymax=log(xpr1)-.5D0*log(xh)
      ymin=0.5D0*log(xh)-log(xmr1)
      do i=0,51
        ranhis(i)=0.d0
      enddo

      call xhistomran11(ranhis,xh,xpr1,xmr1)



      write(ifhi,'(a)')       'openhisto name Om1yk'
      write(ifhi,'(a)')       'htyp lin'
      write(ifhi,'(a)')       'xmod lin ymod lin'
      write(ifhi,'(a,2e11.3)')'xrange',ymin-1.d0,ymax+1.d0
      write(ifhi,'(a)')'yrange auto auto'
      write(ifhi,'(a)')    'text 0 0 "xaxis Y"'
      write(ifhi,'(a)')    'text 0 0 "yaxis [w]?1!(y)"'
      if (xpar8.eq.1.) then
        write(ifhi,'(a,e7.2,a)')  'text 0.1 0.9 "s=',engy**2,' GeV^2!"'
        write(ifhi,'(a,f5.2,a)')  'text 0.5 0.9 "b=',biniDf,' fm"'
      endif
      write(ifhi,'(a)')       'array 2'

      do i=0,nptg
        y=ymin+(ymax-ymin)*(dble(i)/dble(nptg))
        t=om1yk(xh,y,xpr1,xmr1,1)
        write(ifhi,*) y,t
      enddo

      write(ifhi,'(a)')    '  endarray'
      write(ifhi,'(a)')    'closehisto plot 0-'

      write(ifhi,'(a)')       'openhisto name Om1yRan'
      write(ifhi,'(a)')       'htyp pnt'
      write(ifhi,'(a)')       'xmod lin ymod lin'
      write(ifhi,'(a,2e11.3)')'xrange',ymin-1.d0,ymax+1.d0
      write(ifhi,'(a)')'yrange auto auto'
      write(ifhi,'(a)')    'text 0 0 "xaxis Y"'
      write(ifhi,'(a)')    'text 0 0 "yaxis [w]?1!(y) random"'
      if (xpar8.eq.1.) then
        write(ifhi,'(a,e7.2,a)')  'text 0.1 0.9 "s=',engy**2,' GeV^2!"'
        write(ifhi,'(a,f5.2,a)')  'text 0.5 0.9 "b=',biniDf,' fm"'
      endif
      write(ifhi,'(a)')       'array 2'

      do i=0,50
        y=ymin+(ymax-ymin)*(dble(i)/50.d0)
        write(ifhi,*) y,ranhis(i)
      enddo

      write(ifhi,'(a)')    '  endarray'
      write(ifhi,'(a)')    'closehisto plot 0'


      endif

      return
      end

c----------------------------------------------------------------------
      subroutine xRanPt
c----------------------------------------------------------------------
c xpar2=xcut
c----------------------------------------------------------------------
#include "aaa.h"
      parameter (nptg1=501)      !number of point for the graphs
      double precision ranhis(0:nptg1)
      common /cranpt/conv

      nptg=nptg1-1
      xcut=xpar2              !value of biniDf (impact paramter)
      xfact=xpar3
      xadd=xpar4
      xmax=10.
      conv=10./float(nptg)
      if(xcut.le.0.)xcut=float(nptg)
      if(xfact.le.0.)xfact=1.
      do i=0,nptg1
        ranhis(i)=0.d0
      enddo
c$$$      xp=dsqrt(xh)*dble(exp(xpar3)) !y=xpar3
c$$$      xm=1.d0
c$$$      if(xp.ne.0.d0)xm=xh/xp

      if(xpar1.ge.1.)then

      call xranptg(ranhis,xcut,xfact,xadd)



      write(ifhi,'(a)')       'openhisto name ranpt'
      write(ifhi,'(a)')       'htyp lin'
      write(ifhi,'(a)')       'xmod lin ymod log'
      write(ifhi,'(a,2e11.3)')'xrange',0.,min(xmax,xfact*xcut+xadd+1.)
      write(ifhi,'(a)')'yrange auto auto'
      write(ifhi,'(a)')    'text 0 0 "xaxis pt"'
      write(ifhi,'(a)')    'text 0 0 "yaxis P"'
      write(ifhi,'(a)')       'array 2'


      do i=0,nptg
        x=float(i)*conv
        write(ifhi,*) x,ranhis(i)
      enddo

      write(ifhi,'(a)')    '  endarray'
      write(ifhi,'(a)')    'closehisto plot 0-'

      endif

      if(xpar1.ge.2.)then

      call xranpte(ranhis,xcut,xfact,xadd)



      write(ifhi,'(a)')       'openhisto name ranpt'
      write(ifhi,'(a)')       'htyp lin'
      write(ifhi,'(a)')       'xmod lin ymod log'
      write(ifhi,'(a,2e11.3)')'xrange',0.,min(xmax,xfact*xcut+xadd+1.)
      write(ifhi,'(a)')'yrange auto auto'
      write(ifhi,'(a)')    'text 0 0 "xaxis pt"'
      write(ifhi,'(a)')    'text 0 0 "yaxis P"'
      write(ifhi,'(a)')       'array 2'


      do i=0,nptg
        x=float(i)*conv
        write(ifhi,*) x,ranhis(i)
      enddo

      write(ifhi,'(a)')    '  endarray'
      write(ifhi,'(a)')    'closehisto plot 0-'


      endif

      if(xpar1.ge.3.)then

      call xranpts(ranhis,xcut,xfact,xadd)



      write(ifhi,'(a)')       'openhisto name ranpt'
      write(ifhi,'(a)')       'htyp lin'
      write(ifhi,'(a)')       'xmod lin ymod log'
      write(ifhi,'(a,2e11.3)')'xrange',0.,min(xmax,xfact*xcut+xadd+1.)
      write(ifhi,'(a)')'yrange auto auto'
      write(ifhi,'(a)')    'text 0 0 "xaxis pt"'
      write(ifhi,'(a)')    'text 0 0 "yaxis P"'
      write(ifhi,'(a)')       'array 2'


      do i=0,nptg
        x=float(i)*conv
        write(ifhi,*) x,ranhis(i)
      enddo

      write(ifhi,'(a)')    '  endarray'
      write(ifhi,'(a)')    'closehisto plot 0-'

      endif

      call xranptc(ranhis,xcut,xfact,xadd)



      write(ifhi,'(a)')       'openhisto name ranpt'
      write(ifhi,'(a)')       'htyp lin'
      write(ifhi,'(a)')       'xmod lin ymod log'
      if(xpar1.ge.1)then
        write(ifhi,'(a,2e11.3)')'xrange',0.,min(xmax,xfact*xcut+xadd+1.)
      else
        write(ifhi,'(a,2e11.3)')'xrange',0.,1.
      endif
      write(ifhi,'(a)')'yrange auto auto'
      write(ifhi,'(a)')    'text 0 0 "xaxis pt"'
      write(ifhi,'(a)')    'text 0 0 "yaxis P"'
      write(ifhi,'(a)')       'array 2'


      do i=0,nptg
        x=float(i)*conv
        write(ifhi,*) x,ranhis(i)
      enddo

      write(ifhi,'(a)')    '  endarray'


      return
      end

c----------------------------------------------------------------------

      double precision function xgammag2(iclrem,bet,ip,gamp)

c----------------------------------------------------------------------
#include "aaa.h"
#include "sem.h"
      double precision utgam2
      dimension bet(idxD0:idxD1),ip(idxD0:idxD1)

      xgammag2=1.d0

      do i=idxDmin(iomega),idxDmax(iomega)
      if(ip(i).ne.0) xgammag2=xgammag2
     &   *(utgam2(dble(alplea(iclrem))+1.d0+dble(gamp))
     &   /(max(0.d0,dble(int(gamp+0.5))+1))
     &   /utgam2(dble(alplea(iclrem)+bet(i)+gamp)+1.D0))
     &                                          **dble(ip(i))
      enddo

      return
      end

c----------------------------------------------------------------------

      function xsigmafit(x)

c----------------------------------------------------------------------

#include "aaa.h"
#include "par.h"
      double precision x,xDfit,sfsh,varifit,range,sig2
      double precision bf(maxdataDf),Db(maxdataDf)
      external varifit


      sig2=bmaxDf/2.
      range=sig2
      xp=sngl(dsqrt(x))
      xm=xp
      zz=0

      sfsh=xDfit(zz,idxDmin(iomega),idxDmax(iomega),smaxDf,xp,xm,0.)
      if(dabs(sfsh).ge.1.d-10)then
      do i=0,nptf-1
        bf(i+1)=dble(-bmaxDf+float(i)*2.*bmaxDf/float(nptf-1))
        Db(i+1)=xDfit(zz,idxDmin(iomega),idxDmax(iomega),smaxDf,xp,xm
     .               ,sngl(bf(i+1)))/sfsh
c        write(ifch,*)'b,G',bf(i+1),Db(i+1),sfsh
      enddo

c.....Fit of D(X,b) between -bmaxDf and bmaxDf
      call minfit(varifit,bf,Db,nptf,sig2,range)

      xsigmafit=sngl(sig2)
c      write(ifch,*)'fit',xsigmafit
      else
      xsigmafit=0.
      endif


      return
      end


c----------------------------------------------------------------------

      subroutine xhistomran1(histo,xpr,xmr)

c----------------------------------------------------------------------
c.....Make Histogram of om1xpork
c----------------------------------------------------------------------

#include "par.h"

      double precision histo(0:51),x,x1,om1xprk,xpr,xmr,om1xo,om1xoI
      external om1xo,om1xoI
      integer*4 n


      n=100000
      do i=0,51
        histo(i)=0.d0
      enddo
      do j=1,n
        if(mod(j,10000).eq.0)write(*,*)"x1",j
        x=om1xprk(om1xo,om1xoI,1,1,xpr,xmr,10)
        if(x.lt.xminDf)goto 111
c.........Exponential
          k=int((-dlog(x)/dlog(xminDf)+1.d0)*51.d0)
c.........Linear
c.........k=int(x*50.d0)
        histo(k)=histo(k)+1.d0
 111  continue
      enddo
      do i=0,51

c.......Exponential

        x=xminDf
        x1=xminDf
        x=x**(1.d0-dble(i)/51.d0)
        x1=x1**(1.d0-dble(i+1)/51.d0)
        if(i.eq.51)then
          x1=1.d0
          x=0.d0
        endif
        histo(i)=histo(i)/dble(n)/(x1-x)

c.......Linear
c        histo(i)=histo(i)/dble(n)*51.d0
      enddo

      return
      end


c----------------------------------------------------------------------

      subroutine xhistomran2(histo,xp,xpr,xmr)

c----------------------------------------------------------------------
c.....Make Histogram of om1xmork
c----------------------------------------------------------------------

#include "par.h"

      double precision histo(0:51),x,x1,om1xmrk,xp,xpr,xmr,om1xo,om1xoI
      external om1xo,om1xoI
      integer*4 n


      n=100000
      do i=0,51
        histo(i)=0.d0
      enddo
      do j=1,n
        if(mod(j,10000).eq.0)write(*,*)"x-",j
          x=om1xmrk(om1xo,om1xoI,1,1,xp,xpr,xmr,10)
        if(x.lt.xminDf)goto 111
c.........Exponential
          k=int((-dlog(x)/dlog(xminDf)+1.d0)*51.d0)
c.........Linear
c.........k=int(x*50.d0)
        histo(k)=histo(k)+1.d0
 111  continue
      enddo
      do i=0,51

c.......Exponential

        x=xminDf
        x1=xminDf
        x=x**(1.d0-dble(i)/51.d0)
        x1=x1**(1.d0-dble(i+1)/51.d0)
        if(i.eq.51)then
          x1=1.d0
          x=0.d0
        endif
        histo(i)=histo(i)/dble(n)/(x1-x)

c.......Linear
c        histo(i)=histo(i)/dble(n)*51.d0
      enddo

      return
      end

cc----------------------------------------------------------------------
c
c      subroutine xhistomran6(histo,bx,by,bmax,del)
c
cc----------------------------------------------------------------------
cc.....Make Histogram of b1 (impact parameter of vertex in Y and X)
cc----------------------------------------------------------------------
c
c      double precision histo(0:51),dx
c      integer*4 n
c
c      n=100000
c      dx=dble(bmax)/50.d0
c      do i=0,50
c        histo(i)=0.d0
c      enddo
c      do j=1,n
c        if(mod(j,10000).eq.0)write(*,*)"b1",j
c        z=rangen()
c        zp=rangen()
c        bb1x=(bx+sqrt(-del*log(z))*cos(2.*3.14*zp))/2.
c        bb1y=(by+sqrt(-del*log(z))*sin(2.*3.14*zp))/2.
c        x=sqrt((bx-bb1x)*(bx-bb1x)+(by-bb1y)*(by-bb1y))
c        k=int(x/bmax*50.)
c        if(k.le.50)then
c          histo(k)=histo(k)+1.d0
c        else
c          histo(51)=histo(51)+1.d0
c        endif
c      enddo
c      do i=0,50
c        histo(i)=histo(i)/dble(n)/dx
c      enddo
c
c      return
c      end
c

c----------------------------------------------------------------------

      subroutine xhistomran8(histo,xpr,xmr)

c----------------------------------------------------------------------
c.....Make Histogram of om1xprk
c----------------------------------------------------------------------

#include "par.h"

      double precision histo(0:51),x,x1,om1xprk,xpr,xmr,om1x,om1xI
      external om1x,om1xI
      integer*4 n


      n=100000
      do i=0,51
        histo(i)=0.d0
      enddo
      do j=1,n
        if(mod(j,10000).eq.0)write(*,*)"x+",j,xmr
          x=om1xprk(om1x,om1xI,1,1,xpr,xmr,1)
c          x=om1xprk(om1x,om1xI1,1,xpr,xminDf,1)
        if(x.lt.xminDf)goto 111
c.........Exponential
          k=int((-dlog(x)/dlog(xminDf)+1.d0)*51.d0)
c.........Linear
c.........k=int(x*50.d0)
        histo(k)=histo(k)+1.d0
 111  continue
      enddo
      do i=0,51

c.......Exponential

        x=xminDf
        x1=xminDf
        x=x**(1.d0-dble(i)/51.d0)
        x1=x1**(1.d0-dble(i+1)/51.d0)
        if(i.eq.51)then
          x1=1.d0
          x=0.d0
        endif
        histo(i)=histo(i)/dble(n)/(x1-x)

c.......Linear
c        histo(i)=histo(i)/dble(n)*51.d0
      enddo

      return
      end

c----------------------------------------------------------------------

      subroutine xhistomran9(histo,xp,xpr,xmr)

c----------------------------------------------------------------------
c.....Make Histogram of om1xmrk
c----------------------------------------------------------------------

#include "par.h"

      double precision histo(0:51),x,x1,om1xmrk,xp,xpr,xmr,om1x,om1xI
      external om1x,om1xI
      integer*4 n


      n=100000
      do i=0,51
        histo(i)=0.d0
      enddo
      do j=1,n
        if(mod(j,10000).eq.0)write(*,*)"x-",j
          x=om1xmrk(om1x,om1xI,1,1,xp,xpr,xmr,1)
c          x=om1xmrk(om1x,om1xI,1,1,xp,xpr,xminDf,1)
        if(x.lt.xminDf)goto 111
c.........Exponential
          k=int((-dlog(x)/dlog(xminDf)+1.d0)*51.d0)
c.........Linear
c.........k=int(x*50.d0)
        histo(k)=histo(k)+1.d0
 111  continue
      enddo
      do i=0,51

c.......Exponential

        x=xminDf
        x1=xminDf
        x=x**(1.d0-dble(i)/51.d0)
        x1=x1**(1.d0-dble(i+1)/51.d0)
        if(i.eq.51)then
          x1=1.d0
          x=0.d0
        endif
        histo(i)=histo(i)/dble(n)/(x1-x)

c.......Linear
c        histo(i)=histo(i)/dble(n)*51.d0
      enddo

      return
      end

c----------------------------------------------------------------------

      subroutine xhistomran10(histo,xpr,xmr)

c----------------------------------------------------------------------
c.....Make Histogram of om1xrk
c----------------------------------------------------------------------

#include "par.h"

      double precision histo(0:51),x,x1,om1xrk,xpr,xmr
      integer*4 n


      n=100000
      do i=0,51
        histo(i)=0.d0
      enddo
      do j=1,n
        if(mod(j,10000).eq.0)write(*,*)"xk",j
        x=om1xrk(1,1,xpr,xmr)
        if(x.lt.xminDf)goto 111
c.........Exponential
          k=int((-dlog(x)/dlog(xminDf)+1.d0)*51.d0)
c.........Linear
c.........k=int(x*50.d0)
        histo(k)=histo(k)+1.d0
 111  continue
      enddo
      do i=0,51

c.......Exponential

        x=xminDf
        x1=xminDf
        x=x**(1.d0-dble(i)/51.d0)
        x1=x1**(1.d0-dble(i+1)/51.d0)
        if(i.eq.51)then
          x1=1.d0
          x=0.d0
        endif
        histo(i)=histo(i)/dble(n)/(x1-x)

c.......Linear
c        histo(i)=histo(i)/dble(n)*51.d0
      enddo

      return
      end


c----------------------------------------------------------------------

      subroutine xhistomran11(histo,xh,xpr,xmr)

c----------------------------------------------------------------------
c.....Make Histogram of om1yrk
c----------------------------------------------------------------------

      double precision histo(0:51),x,xh,dx,om1yrk,ymax,ymin,xpr,xmr
      integer*4 n

      ymax=log(xpr)-.5D0*log(xh)
      ymin=0.5D0*log(xh)-log(xmr)
      dx=(ymax-ymin)/50.d0

      n=100000
      do i=0,50
        histo(i)=0.d0
      enddo
      do j=1,n
        if(mod(j,10000).eq.0)write(*,*)"yk",j
        x=om1yrk(1,1,xh,xpr,xmr)
        k=int((x-ymin)/(ymax-ymin)*50.d0)
c.......write(*,*)x,k
        histo(k)=histo(k)+1.d0
      enddo
      do i=0,50
        histo(i)=histo(i)/dble(n)/dx
      enddo

      return
      end

c----------------------------------------------------------------------

      subroutine xranptg(histo,xcut,xfact,xadd)

c----------------------------------------------------------------------
c.....Make Histogram of random distribution
c----------------------------------------------------------------------

#include "par.h"

      parameter (nptg1=501)      !number of point for the graphs
      common /cranpt/conv
      double precision histo(0:nptg1)
      integer*4 n


      n=100000
      do i=0,nptg1
        histo(i)=0.d0
      enddo
      do j=1,n
        if(mod(j,10000).eq.0)write(*,*)"ptg",j
c .........exp(-x**2)
 12   x=sqrt(-log(rangen())/(3.1415927/4.)) !gauss

      if(xcut.gt.0.)then
        if(rangen().lt.x/xcut)goto 12
      endif
      x=x*xfact+xadd
c.........Exponential
c          k=int((-dlog(x)/dlog(xminDf)+1.d0)*51.d0)
c.........Linear
        k=int(x/conv)
        k=min(k,nptg1)
        histo(k)=histo(k)+1.d0
      enddo
      do i=0,nptg1

c.......Exponential

c        x=xminDf
c        x1=xminDf
c        x=x**(1.d0-dble(i)/51.d0)
c        x1=x1**(1.d0-dble(i+1)/51.d0)
c        if(i.eq.51)then
c          x1=1.d0
c          x=0.d0
c        endif
c        histo(i)=histo(i)/dble(n)/(x1-x)

c.......Linear
        histo(i)=histo(i)/dble(n)*float(nptg1)
      enddo

      return
      end

c----------------------------------------------------------------------

      subroutine xranpte(histo,xcut,xfact,xadd)

c----------------------------------------------------------------------
c.....Make Histogram of random distribution
c----------------------------------------------------------------------

#include "par.h"

      parameter (nptg1=501)      !number of point for the graphs
      common /cranpt/conv
      double precision histo(0:nptg1)
      integer*4 n


      n=100000
      do i=0,nptg1
        histo(i)=0.d0
      enddo
      do j=1,n
        if(mod(j,10000).eq.0)write(*,*)"pte",j
c .........exp(-x)
  12  xmx=50
      x=0.
      r=2.
      do while (r.gt.1.)
  11    x=sqrt(exp(rangen()*log(1+xmx**2))-1)
        if(x.eq.0.)goto11
        r=rangen()  /  ( exp(-x)*(1+x**2) )
      enddo
      x=x/2.

      if(xcut.gt.0.)then
        if(rangen().lt.x/xcut)goto 12
      endif
      x=x*xfact+xadd
c.........Exponential
c          k=int((-dlog(x)/dlog(xminDf)+1.d0)*51.d0)
c.........Linear
        k=int(x/conv)
        k=min(k,nptg1)
        histo(k)=histo(k)+1.d0
      enddo
      do i=0,nptg1

c.......Exponential

c        x=xminDf
c        x1=xminDf
c        x=x**(1.d0-dble(i)/51.d0)
c        x1=x1**(1.d0-dble(i+1)/51.d0)
c        if(i.eq.51)then
c          x1=1.d0
c          x=0.d0
c        endif
c        histo(i)=histo(i)/dble(n)/(x1-x)

c.......Linear
        histo(i)=histo(i)/dble(n)*float(nptg1)
      enddo

      return
      end

c----------------------------------------------------------------------

      subroutine xranpts(histo,xcut,xfact,xadd)

c----------------------------------------------------------------------
c.....Make Histogram of random distribution
c----------------------------------------------------------------------

#include "par.h"

      parameter (nptg1=501)      !number of point for the graphs
      common /cranpt/conv
      double precision histo(0:nptg1)
      integer*4 n


      n=100000
      do i=0,nptg1
        histo(i)=0.d0
      enddo
      do j=1,n
        if(mod(j,10000).eq.0)write(*,*)"pts",j
c .........exp(-sqrt(x))
 12   xmx=500
      x=0.
      r=2.
      do while (r.gt.1.)
        x=sqrt(exp(rangen()*log(1+xmx**2))-1)
        r=rangen()  /  ( exp(-sqrt(x))*(1+x**2)/5. )
      enddo
      x=x/20.

      if(xcut.gt.0.)then
        if(rangen().lt.x/xcut)goto 12
      endif
      x=x*xfact+xadd
c.........Exponential
c          k=int((-dlog(x)/dlog(xminDf)+1.d0)*51.d0)
c.........Linear
        k=int(x/conv)
        k=min(k,nptg1)
        histo(k)=histo(k)+1.d0
      enddo
      do i=0,nptg1

c.......Exponential

c        x=xminDf
c        x1=xminDf
c        x=x**(1.d0-dble(i)/51.d0)
c        x1=x1**(1.d0-dble(i+1)/51.d0)
c        if(i.eq.51)then
c          x1=1.d0
c          x=0.d0
c        endif
c        histo(i)=histo(i)/dble(n)/(x1-x)

c.......Linear
        histo(i)=histo(i)/dble(n)*float(nptg1)
      enddo

      return
      end

c----------------------------------------------------------------------

      subroutine xranptc(histo,xcut,xfact,xadd)

c----------------------------------------------------------------------
c.....Make Histogram of random distribution
c----------------------------------------------------------------------

#include "par.h"

      parameter (nptg1=501)      !number of point for the graphs
      common /cranpt/conv
      double precision histo(0:nptg1)
      integer*4 n


      n=100000
      do i=0,nptg1
        histo(i)=0.d0
      enddo
      do j=1,n
        if(mod(j,10000).eq.0)write(*,*)"ptc",j

        x=ranptcut(xcut)*xfact+xadd
c.........Exponential
c          k=int((-dlog(x)/dlog(xminDf)+1.d0)*51.d0)
c.........Linear
        k=int(x/conv)
        k=min(k,nptg1)
        histo(k)=histo(k)+1.d0
      enddo
      do i=0,nptg1

c.......Exponential

c        x=xminDf
c        x1=xminDf
c        x=x**(1.d0-dble(i)/51.d0)
c        x1=x1**(1.d0-dble(i+1)/51.d0)
c        if(i.eq.51)then
c          x1=1.d0
c          x=0.d0
c        endif
c        histo(i)=histo(i)/dble(n)/(x1-x)

c.......Linear
        histo(i)=histo(i)/dble(n)*float(nptg1)
      enddo

      return
      end



