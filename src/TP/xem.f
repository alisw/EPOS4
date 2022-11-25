C
C  This file is part of EPOS4
C  Copyright (C) 2022 research institutions and authors (See CREDITS file)
C  This file is distributed under the terms of the GNU General Public License version 3 or later
C  (See COPYING file for the text of the licence)
C




c###########################################################################
c###########################################################################
c###########################################################################
c###########################################################################
c
c      These are plot routines copied from epos-ems
c
c         not (yet) updated concerning centrality dependent output
c                     and new     q 2 m i n    dependence                                              
c
c###########################################################################
c###########################################################################
c###########################################################################
c###########################################################################


c-----------------------------------------------------------------------
      subroutine xEmsP3fill
c-----------------------------------------------------------------------
      parameter(nbixpdf=35,nbiqq=24,mmxcentrx=9)
      double precision wpdfxx,wpdf
      common/cemspdfxx/wpdfxx(nbixpdf,0:5,2,2,nbiqq,2)
      common/cxpdf/xlpdf(nbixpdf),dxlpdf(nbixpdf)
     *            ,qqval(nbiqq),dqqval(nbiqq)
     *            ,wpdf(nbixpdf,0:5,2,2,nbiqq,2,0:mmxcentrx+1)
      common/cncoll/wncoll(3,0:mmxcentrx+1),wnevt(3,0:mmxcentrx+1)
     .             ,wbim(0:mmxcentrx+1),wpom(0:7,0:mmxcentrx+1)
      common/cemsp3q2/q2kkk3(2),nq2kkk3
      common/cwq2/wq2(3,0:mmxcentrx+1,3)
#include "aaa.h"
      ncyi=mcentrf()
      if(ncyi.lt.0.or.ncyi.gt.mmxcentrf())stop'\n\n ERROR 27092012g\n\n'
      if(ncyi.eq.0)ncyi=mmxcentrx+1

      nko=nglevt
      if(maproj.eq.1.and.matarg.eq.1)nko=1  !ikoevt   !all scatterings should be consodered to built pdf ?

      do no=1,2

      if(no.eq.1)then
        ff=1
        if(nko.ne.0)then
         ff=1./nko
        endif
        ncy=ncyi
      else
        ff=1
        ncy=0
      endif

      do ipt=1,2     !proj/targ
      do j=1,nbiqq   !qq
      do jaai=1,2    !parton type (sea/val)
      do jexi=1,2    !with/out emission
      do jidi=0,5    !final parton id
      do i=1,nbixpdf
        wpdf(i,jidi,jexi,jaai,j,ipt,ncy)=
     .                      wpdf(i,jidi,jexi,jaai,j,ipt,ncy)
     .                     +wpdfxx(i,jidi,jexi,jaai,j,ipt)*ff
      enddo
      enddo
      enddo
      enddo
      enddo
      enddo
      nko=nglevt
      if(maproj.eq.1.and.matarg.eq.1)nko=ikoevt
      if(nq2kkk3.gt.0)wncoll(3,ncy)=wncoll(3,ncy)+nko
      if(nq2kkk3.gt.0)wnevt(3,ncy)=wnevt(3,ncy)+1
      wq2(3,ncy,1)=wq2(3,ncy,1)+q2kkk3(1)
      wq2(3,ncy,2)=wq2(3,ncy,2)+q2kkk3(2)
      wq2(3,ncy,3)=wq2(3,ncy,3)+nq2kkk3

      enddo

      end

c-----------------------------------------------------------------------
      subroutine xxEmsP3(iii,jaa,jex,id1,id2,xpb,xmb,qq1,qq2,q2kk)
c-----------------------------------------------------------------------
c prepare  real F2 taken from MC
c  iii=1: fill arrays
c  jaa: type of semihard Pomeron
c         0= soft+gg,
c         1= sea-sea, 2= val=sea, 3= sea-val, 4= val-val
c  jex: emission type
c         1= no emission, 2= proj emis, 3= targ emis, 4= both sides
c id1,id2: parton id for born
c xpb,xmb: momentum fractions for born
c qq1,qq2: Q2 for born
c q2kk: initial Q20
c-----------------------------------------------------------------------

#include "aaa.h"
#include "sem.h"
#include "ems.h"
      parameter(nbixpdf=35,nbiqq=24)
      double precision wpdfxx
      common/cemspdfxx/wpdfxx(nbixpdf,0:5,2,2,nbiqq,2)
      common/cemspdfbx/xlup1,xlop,qqug,qqog
      common/cemsp3q2/q2kkk3(2),nq2kkk3
      dimension q2kk(2)

      if(iemspdf.eq.0)call utstop('ERROR in xEmsP3: iemspdf = 0&')

      if(iii.eq.1.and.jaa.eq.-1)then

      do ipt=1,2      !proj/targ
      do j=1,nbiqq   !qq
      do jaai=1,2    !parton type (sea/val)
      do jexi=1,2    !with/out emission
      do jidi=0,5    !final parton id
      do i=1,nbixpdf
        wpdfxx(i,jidi,jexi,jaai,j,ipt)=0d0
      enddo
      enddo
      enddo
      enddo
      enddo
      enddo
      q2kkk3(1)=0
      q2kkk3(2)=0
      nq2kkk3=0


      elseif(iii.eq.1)then

       if(jaa.eq.0.or.jaa.eq.1)then
         ja1=1
         ja2=1
       elseif(jaa.eq.2)then
         ja1=2
         ja2=1
       elseif(jaa.eq.3)then
         ja1=1
         ja2=2
       elseif(jaa.eq.4)then
         ja1=2
         ja2=2
       else
         goto 42
       endif
       if(jex.eq.1)then
         je1=1
         je2=1
       elseif(jex.eq.2)then
         je1=2
         je2=1
       elseif(jex.eq.3)then
         je1=1
         je2=2
       elseif(jaa.eq.4)then
         je1=2
         je2=2
       else
         goto 42
       endif

c proj side
       xx=xpb
       j=1+int(log(qq1/qqug)/log(qqog/qqug)*nbiqq)
       if(j.gt.nbiqq)goto 2
       if(j.lt.1)goto 2
c log scale
       if(xx.lt.xlup1)goto 2
       i=1+int(log(xx/xlup1)/log(xlop/xlup1)*nbixpdf)
       if(i.gt.nbixpdf)goto 2
       if(i.lt.1)goto 2
       wpdfxx(i,id1,je1,ja1,j,1)=wpdfxx(i,id1,je1,ja1,j,1)+1d0
2      continue

c target side
       xx=xmb
       j=1+int(log(qq2/qqug)/log(qqog/qqug)*nbiqq)
       if(j.gt.nbiqq)goto 22
       if(j.lt.1)goto 22
c log scale
       if(xx.lt.xlup1)goto 22
       i=1+int(log(xx/xlup1)/log(xlop/xlup1)*nbixpdf)
       if(i.gt.nbixpdf)goto 22
       if(i.lt.1)goto 22
       wpdfxx(i,id2,je2,ja2,j,2)=wpdfxx(i,id2,je2,ja2,j,2)+1d0
22      continue
       
       q2kkk3(1)=q2kkk3(1)+q2kk(1)
       q2kkk3(2)=q2kkk3(2)+q2kk(2)
       nq2kkk3=nq2kkk3+1
42     continue

      else

          write(ifch,*)'\n\n  ERROR 12012017a \n\n'
                   stop'\n\n  ERROR 12012017a \n\n'

      endif

      return
      end

c-----------------------------------------------------------------------
      subroutine xEmsP3(iii,jaa,jex,idi,qq,ipt)
c-----------------------------------------------------------------------
c plot  real F2 taken from MC
c  iii>=2: write arrays
c  jaa: type of parton
c         0= only sum
c         1= sea
c         2= val,
c         5= all
c  jex: emission type
c         0= only sum
c         1= no emission
c         2= with emission
c         5= all
c idi: parton id for born (6=all quarks)
c qq: Q2 at born (-1.= make a set of standard values with data)
c ipt: interaction side
c         1=proj
c         2=targ
c-----------------------------------------------------------------------

#include "aaa.h"
#include "sem.h"
#include "ems.h"
      common/geom/rmproj,rmtarg,bmax,bkmx
      parameter(nbixpdf=35,nbiqq=24,mmxcentrx=9)
      double precision wpdfxx,wpdf,z
      common/cemspdfxx/wpdfxx(nbixpdf,0:5,2,2,nbiqq,2)
      common/cemspdfbx/xlup1,xlop,qqug,qqog
      common/cxpdf/xlpdf(nbixpdf),dxlpdf(nbixpdf)
     *            ,qqval(nbiqq),dqqval(nbiqq)
     *            ,wpdf(nbixpdf,0:5,2,2,nbiqq,2,0:mmxcentrx+1)
      common/cwq2/wq2(3,0:mmxcentrx+1,3)
      common/cncoll/wncoll(3,0:mmxcentrx+1),wnevt(3,0:mmxcentrx+1)
     .             ,wbim(0:mmxcentrx+1),wpom(0:7,0:mmxcentrx+1)
      character cii*1,dii*7,c3*3,mod*5
      dimension q2amin(2)
      double precision pifpartone

      if(iemspdf.eq.0)call utstop('ERROR in xEmsP3: iemspdf = 0&')
      if(mmxcentrf().gt.mmxcentrx)stop'\n\n ERROR 27092012h \n\n '

      if(iii.eq.0)then

c       write(*,*)'initialize xEmsP3 tables for'
c     . ,mmxcentrf(),' centralities'
       xlup1=0.01/engy
       xlop=1.
       qqug=q2nmin
       qqog=1.e4
       do i=1,nbixpdf
        xlpdf(i)=xlup1*(xlop/xlup1)**((i-0.5)/nbixpdf)
        dxlpdf(i)=xlup1*(xlop/xlup1)**(1.*i/nbixpdf)
     *             *(1.-(xlop/xlup1)**(-1./nbixpdf))
       enddo
       do i=1,nbiqq
        qqval(i)=qqug*(qqog/qqug)**((i-0.5)/nbiqq)
        dqqval(i)=qqug*(qqog/qqug)**(1.*i/nbiqq)
     *             *(1.-(qqog/qqug)**(-1./nbiqq))
       enddo
       do ncy=0,mmxcentrf()      !centrality
       do jpt=1,2     !proj/targ
       do j=1,nbiqq   !qq
       do jaai=1,2    !parton type (sea/val)
       do jexi=1,2    !with/out emission
       do jidi=0,5    !final parton id
       do i=1,nbixpdf
         wpdf(i,jidi,jexi,jaai,j,jpt,ncy)=0d0
       enddo
       enddo
       enddo
       enddo
       enddo
       enddo
       wncoll(3,ncy)=0.
       wnevt(3,ncy)=0.
       do j=1,3
         wq2(3,ncy,j)=0
       enddo
       enddo

      elseif(iii.ge.2)then

       if(mmxcentrf().gt.0)then
         ncy1=1
       else
         ncy1=0
       endif
       ncy2=mmxcentrf()

       if(nint(xpar7).eq.1)then
         ncy1=0
         ncy2=0
       endif

       nstat=nfull
       if(nstat.lt.0)nstat=max(1,nevent)

       do ncy=ncy1,ncy2

       if(wq2(3,ncy,3).gt.0.)then
         hww=1
         ihww=1
         hw=wq2(3,ncy,3)
         ff=30.*wnevt(3,ncy)/hw
         q2amin(1)=wq2(3,ncy,1)/wq2(3,ncy,3) 
         q2amin(2)=wq2(3,ncy,2)/wq2(3,ncy,3)
         q2cmin(1)=max(q2nmin,q2amin(1))
         q2cmin(2)=max(q2nmin,q2amin(2))
       else
         hww=0
         ihww=0
         ff=1.
         q2cmin(1)=q2nmin
         q2cmin(2)=q2nmin
       endif
         q2pmin(1)=q2cmin(1)
         q2pmin(2)=q2cmin(2)
         
         call ipoCSzeroTables(q2cmin)
         call ipoOm5Tables(1)

         
       if(maproj.eq.1.and.matarg.eq.1.and.bminim.eq.bmaxim)then
         ff=ff*float(nstat)/float(ntevt)
       elseif(maproj.eq.1.and.matarg.eq.1)then
         if(wnevt(3,ncy).gt.0.)ff=ff*nstat/wnevt(3,ncy)
c        if(ncy.gt.0)hw=wnevt(3,ncy)
       else
        if(ncy.gt.0)then
          if(wnevt(3,ncy).gt.0.)ff=ff*nstat/wnevt(3,ncy)
c          hw=wnevt(2,ncy)
        elseif(bminim.lt.0.001.and.bmaxim.gt.20)then !min bias
          area=pi*(rmproj+rmtarg)**2
          ff=ff*area*float(nstat)/float(ntevt)/(maproj*matarg)/sigine*10
        else
          write(ifmt,*)'ignored xemsPDF'
          return
        endif
       endif

       xlup=xlup1
       if(qq.gt.0.)then
         j=1+int(log(qq/qqug)/log(qqog/qqug)*nbiqq)
       elseif(xpar5.gt.0.)then  !take qq from xpar5
         j=1+int(log(xpar5/qqug)/log(qqog/qqug)*nbiqq)
       else
         j=1
       endif
c       print*,'xEmsP3',ncy,q2amin,wq2(3,ncy,1),wq2(3,ncy,3)
c     *       ,qq,j,xpar5,qqval(j),q2cmin

       ja1=jaa
       ja2=jaa
       je1=jex
       je2=jex
       if(jaa.eq.5)then
       ja1=1
       ja2=3
       elseif(jaa.eq.0)then
       ja1=3
       ja2=3
       endif
       if(jex.eq.5)then
       je1=1
       je2=3
       elseif(jex.eq.0)then
       je1=3
       je2=3
       endif


       if(ipt.eq.1)then
         cii='p'
         dii='x+?IB!'
       else
         cii='t'
         dii='x-?IB! '
       endif


       mod=' log '

       do ja=ja1,ja2

       do je=je1,je2

       write(ifhi,'(a)') "!-----------------------------"
       write(ifhi,'(a)') "! MC   "
       write(ifhi,'(a)') "!-----------------------------"

        write(ifhi,'(a,4i1)')'openhisto name F2'//cii,idi,ja,je,ncy
        write(ifhi,'(a)')'xmod '//mod//' ymod log'
        write(ifhi,'(a,2e11.3)')'xrange ',xlup,xlop
        write(ifhi,'(a)')'yrange 1e-6 auto'
        write(ifhi,'(a)')'txt "xaxis  '//dii//'"'
        write(ifhi,'(a,f6.1,a)')'txt "yaxis  F?2!('//dii//',Q^2!=',
     .                                           qqval(j),') "'
        if(ihww.eq.1)then
          if(je.eq.je2.and.ja.eq.ja2)then
            write(ifhi,'(a)')'htyp pos'
            c3='7.1'
            if(max(q2amin(1),q2amin(2)).lt.1000)c3='6.1'
            if(max(q2amin(1),q2amin(2)).lt.100)c3='6.2'
            if(max(q2amin(1),q2amin(2)).lt.10)c3='6.3'
            x=0.75
            y=0.9
            write(ifhi,'(a,2f5.2,a,i2,a)')'text',x,y,' "M',ncy,'"'
            write(ifhi,'(a,2f5.2,a,f'//c3//',a)')'text'
     .           ,x,y-0.1,' "Q2sat av',q2amin(1),'"'
            write(ifhi,'(a,2f5.2,a,f'//c3//',a)')'text'
     .           ,x+0.13,y-0.2,' "',q2amin(2),'"'
          else
            write(ifhi,'(a)')'htyp pas'
          endif

          write(ifhi,'(a,e22.14)')'histoweight ',dble(hw)
          write(ifhi,'(a)')     'array 2'
          do i=1,nbixpdf
            u=xlpdf(i)
            z=0.d0
            ja3=ja
            ja4=ja
            if(ja.eq.3)then
              ja3=1
              ja4=2
            endif
            je3=je
            je4=je
            if(je.eq.3)then
              je3=1
              je4=2
            endif
            do ja5=ja3,ja4
            do je5=je3,je4
            if(idi.eq.0)then
              z=z+wpdf(i,0,je5,ja5,j,ipt,ncy) !gluon
            else
              if(idi.eq.1..or.idi.eq.6)
     .        z=z+wpdf(i,1,je5,ja5,j,ipt,ncy)/2.25d0 !u
              if(idi.eq.2.or.idi.eq.6)
     .        z=z+wpdf(i,2,je5,ja5,j,ipt,ncy)/9d0   !d
              if(idi.eq.3.or.idi.eq.6)
     .        z=z+wpdf(i,3,je5,ja5,j,ipt,ncy)/9d0   !s
              if(idi.eq.4.or.idi.eq.6)
     .        z=z+wpdf(i,4,je5,ja5,j,ipt,ncy)/2.25d0 !c
              if(idi.eq.5.or.idi.eq.6)
     .        z=z+wpdf(i,5,je5,ja5,j,ipt,ncy)/9d0   !b
            endif
            enddo
            enddo
            del=dxlpdf(i)
            z=z*u/del*qqval(j)/dqqval(j)*ff/float(max(1,nstat))         !compared to x*F(x)
            write(ifhi,'(2e11.3)')u,z
          enddo
          write(ifhi,'(a)')    '  endarray'
          write(ifhi,'(a)')    'closehisto'
        else
          write(ifhi,'(a)')    'closehisto'
        endif
        write(ifhi,'(a)')    'plot 0-'

        write(ifhi,'(a)')       '!----------------------------------'
        write(ifhi,'(a,4i2)')   '!   fpartone    ',idi,ja,je,ncy
        write(ifhi,'(a)')       '!----------------------------------'

        write(ifhi,'(a,4i1)')'openhisto name fpartone'//cii
     .                                               ,idi,ja,je,ncy
        if(je.eq.je2.and.ja.eq.ja2)then
          write(ifhi,'(a)')'htyp lfu'
          write(ifhi,'(a,2i1,a)')
     .  'txt "title fpartone'//cii//'  (',jaa,jex,')"'
        else
          write(ifhi,'(a)')'htyp lin'
        endif
        write(ifhi,'(a,e22.14)')'histoweight ',dble(hww)
        write(ifhi,'(a)')'array 2'
        q2cmin(ipt)=max(q2nmin,min(qqval(j),q2amin(ipt)))
        do i=1,nbixpdf
          u=xlpdf(i)
          z=0.d0
          ja3=ja
          ja4=ja
          if(ja.eq.3)then
            ja3=1
            ja4=2
          endif
          je3=je
          je4=je
          if(je.eq.3)then
            je3=1
            je4=2
          endif
          do ja5=ja3,ja4
          do je5=je3,je4
            if(idi.eq.0)then
            z=z+pifpartone(ipt,dble(u),qqval(j),0,je5-1,ja5)      !gluon
            else
            if((idi.eq.1.or.idi.eq.6).and.ja5.eq.1)
     .    z=z+2d0*pifpartone(ipt,dble(u),qqval(j),-1,je5-1,ja5)/2.25d0 !u sea
            if((idi.eq.1.or.idi.eq.6).and.ja5.eq.2)
     .    z=z+pifpartone(ipt,dble(u),qqval(j),1,je5-1,ja5)/2.25d0    !u val
            if((idi.eq.2.or.idi.eq.6).and.ja5.eq.1)
     .    z=z+2d0*pifpartone(ipt,dble(u),qqval(j),-1,je5-1,ja5)/9d0 !d sea
            if((idi.eq.2.or.idi.eq.6).and.ja5.eq.2)
     .    z=z+pifpartone(ipt,dble(u),qqval(j),2,je5-1,ja5)/9d0      !d val
            if((idi.eq.3.or.idi.eq.6).and.ja5.eq.1)
     .    z=z+2d0*pifpartone(ipt,dble(u),qqval(j),-3,je5-1,ja5)/9d0 !s sea
            if((idi.eq.4.or.idi.eq.6).and.ja5.eq.1)
     .    z=z+2d0*pifpartone(ipt,dble(u),qqval(j),-4,je5-1,ja5)/2.25d0 !c sea
            if((idi.eq.5.or.idi.eq.6).and.ja5.eq.1)
     .    z=z+2d0*pifpartone(ipt,dble(u),qqval(j),-5,je5-1,ja5)/9d0 !b sea
            endif
          enddo
          enddo
          write(ifhi,'(2e11.3)')u,z
         enddo
         write(ifhi,'(a)')    '  endarray'
         write(ifhi,'(a)')    'closehisto'
         if(.not.(je.eq.je2.and.ja.eq.ja2))write(ifhi,'(a)')'plot 0-'

       enddo       !je

       enddo       !ja
          
       if(ncy.lt.ncy2)then
         write(ifhi,'(a)')'plot 0'
       elseif(nint(xpar8).eq.1)then
         write(ifhi,'(a)')'plot 0'
       endif

      enddo                     ! ncy


      else

          write(ifch,*)'\n\n  ERROR 12012017b \n\n'
                   stop'\n\n  ERROR 12012017b \n\n'

      endif

      return
      end






c-----------------------------------------------------------------------
      subroutine xEmsI1(iii,kc,omlog)
c-----------------------------------------------------------------------
c plot omlog vs iter
c plot  nr of pomerons vs iter
c plot number of collisions vs iter
c-----------------------------------------------------------------------

#include "aaa.h"
#include "ems.h"
#include "sem.h"

      parameter(nbin=100)
      common/cmc/ot(0:nbin),zz(0:nbin),i(0:nbin)
     *,yt1,yt2,kx(0:nbin)
      parameter(nbim=100)
      common/cmc1/xp(0:nbim),xt(0:nbim),x(0:nbim),o(0:nbim)
     *,y1,y2,car
      character car*5
      double precision xp,xt,x,omlog,om1intbc
      character ce*8
      double precision plc,s,seedp
      common/cems5/plc,s

c      if(iemsi2.eq.0)call utstop('ERROR in XemsI1: iemsi2 = 0&')

       if(iii.eq.1)then

      o(kc)=sngl(omlog)
      nptk=0
      kollx=0
      do ko=1,koll
      nptk=nptk+nprt(ko)
c      if(itpr(ko).gt.0)then
      if(nprt(ko).gt.0)then
       kollx=kollx+1
      endif
      enddo
      zz(kc)=nptk
      kx(kc)=kollx

        elseif(iii.eq.2)then

      call ranfgt(seedp)
      nstat=nfull
      if(nstat.lt.0)nstat=max(1,nevent)
      sum=0
      kollx=0
      sumg=0
      kollg=0
      kollini=koll
      iomegasave=iomega
      iomega=2
      call DefXminDf(s)
      koll=1
      do ko=1,kollini
        om1i=sngl(om1intbc(bk(ko)))
        om1g=sngl(om1intbc(bk(ko)))
        sum=sum+om1i
        sumg=sumg+om1g
        if(rangen().lt.1.-exp(-om1i))then
          kollx=kollx+1
        endif
        if(rangen().lt.1.-exp(-om1g))then
          kollg=kollg+1
        endif
      enddo
      iomega=iomegasave
      call DefXminDf(s)
      koll=kollini
      call ranfst(seedp)

      x1=0
      x2=nbin
      write(ce,'(f8.2)')sngl(plc)

      write(ifhi,'(a)')       '!##################################'
      write(ifhi,'(a,i3)')    '!   log omega       for event ',nrevt+1
      write(ifhi,'(a)')       '!##################################'
      write(ifhi,'(a,i1)')    'openhisto name omega-',nrevt+1
      write(ifhi,'(a)')       'htyp lin'
      write(ifhi,'(a)')       'xmod lin ymod lin'
      write(ifhi,'(a,2e11.3)')'xrange',x1,x2
      write(ifhi,'(a)')       'yrange auto auto '
      write(ifhi,'(a)')    'text 0 0 "xaxis iteration"'
      write(ifhi,'(a)')    'text 0 0 "yaxis ln[W]"'
      write(ifhi,'(a,a)')  'text 0.5 0.90 "E ='//ce//'"'
      write(ifhi,'(a,e22.14)')'histoweight ',dble(nstat)
      write(ifhi,'(a)')       'array 2'
         do k=0,nbim
      write(ifhi,'(2e11.3)')float(k),o(k)
         enddo
      write(ifhi,'(a)')    '  endarray'
      write(ifhi,'(a)')    'closehisto plot 0'

      write(ifhi,'(a)')   '!##################################'
      write(ifhi,'(a,i3)')'! nr of coll`s  for event ',nrevt+1
      write(ifhi,'(a)')   '!##################################'
      write(ifhi,'(a,i1)')    'openhisto name coll-',nrevt+1
      write(ifhi,'(a)')       'htyp lin'
      write(ifhi,'(a)')       'xmod lin ymod lin'
      write(ifhi,'(a,2e11.3)')'xrange',x1,x2
      write(ifhi,'(a)')    'text 0 0 "xaxis iteration"'
      write(ifhi,'(a)')    'text 0 0 "yaxis nr of collisions"'
      write(ifhi,'(a,e22.14)')'histoweight ',dble(nstat)
      write(ifhi,'(a)')       ' array 2'
         do k=0,nbin
      write(ifhi,'(2e11.3)')float(k),float(kx(k))
         enddo
      write(ifhi,'(a)')    '  endarray'
      write(ifhi,'(a)')    'closehisto plot 0-'
      write(ifhi,'(a)')       'openhisto'
      write(ifhi,'(a)')       'htyp lin'
c      write(ifhi,'(a,e22.14)')'histoweight ',dble(nstat)
      write(ifhi,'(a)')       ' array 2'
         do k=0,nbin
      write(ifhi,'(2e11.3)')float(k),float(kollx)
         enddo
      write(ifhi,'(a)')    '  endarray'
      write(ifhi,'(a)')    'closehisto plot 0-'
      write(ifhi,'(a)')       'openhisto'
      write(ifhi,'(a)')       'htyp lin'
c      write(ifhi,'(a,e22.14)')'histoweight ',dble(nstat)
      write(ifhi,'(a)')       ' array 2'
         do k=0,nbin
      write(ifhi,'(2e11.3)')float(k),float(kollg)
         enddo
      write(ifhi,'(a)')    '  endarray'
      write(ifhi,'(a)')    'closehisto plot 0'

      write(ifhi,'(a)')   '!##################################'
      write(ifhi,'(a,i3)')'! nr of pom`s  for event ',nrevt+1
      write(ifhi,'(a)')   '!##################################'
      write(ifhi,'(a,i1)')    'openhisto name pom-',nrevt+1
      write(ifhi,'(a)')       'htyp lin'
      write(ifhi,'(a)')       'xmod lin ymod lin'
      write(ifhi,'(a,2e11.3)')'xrange',x1,x2
      write(ifhi,'(a)')    'text 0 0 "xaxis iteration"'
      write(ifhi,'(a)')    'text 0 0 "yaxis nr of Pomerons"'
      write(ifhi,'(a,e22.14)')'histoweight ',dble(nstat)
      write(ifhi,'(a)')       ' array 2'
         do k=0,nbin
      write(ifhi,'(2e11.3)')float(k),zz(k)
         enddo
      write(ifhi,'(a)')    '  endarray'
      if(sum.lt.4*zz(nbin))then
      write(ifhi,'(a)')    'closehisto plot 0-'
      write(ifhi,'(a)')       'openhisto'
      write(ifhi,'(a)')       'htyp lin'
c      write(ifhi,'(a,e22.14)')'histoweight ',dble(nstat)
      write(ifhi,'(a)')       ' array 2'
         do k=0,nbin
      write(ifhi,'(2e11.3)')float(k),sum
         enddo
      write(ifhi,'(a)')    '  endarray'
      write(ifhi,'(a)')    'closehisto plot 0-'
      write(ifhi,'(a)')       'openhisto'
      write(ifhi,'(a)')       'htyp lin'
c      write(ifhi,'(a,e22.14)')'histoweight ',dble(nstat)
      write(ifhi,'(a)')       ' array 2'
         do k=0,nbin
      write(ifhi,'(2e11.3)')float(k),sumg
         enddo
      write(ifhi,'(a)')    '  endarray'
      endif
      write(ifhi,'(a)')    'closehisto plot 0'

        endif

      return
      end

c-----------------------------------------------------------------------
      subroutine xEmsI2(iii,kc)
c-----------------------------------------------------------------------
c plot quanities vs iter
c   plot 1: <x> for Pomeron vs iter
c   plot 2: <x> for projectile vs iter
c   plot 3: <x> for target vs iter
c arguments:
c   iii:   modus (1,2)
c   kc:    iteration step
c   omega: config probability
c-----------------------------------------------------------------------

#include "aaa.h"
#include "ems.h"

      parameter(nbim=100)
      common/cmc1/xp(0:nbim),xt(0:nbim),x(0:nbim),o(0:nbim)
     *,y1,y2,car
      character car*5
      double precision xp,xt,x,xpo,xpj,xtg
      common/cemsi2/xpo,xpj,xtg

        if(iii.eq.1)then

      npom=0
      xpo=0
      do k=1,koll
c       ip=iproj(k)
c       it=itarg(k)
       if(nprmx(k).gt.0)then
        do n=1,nprmx(k)
         if(idpr(n,k).gt.0.and.ivpr(n,k).gt.0)then
          xpo=xpo+xpr(n,k)
          npom=npom+1
         endif
        enddo
       endif
      enddo
      if(npom.gt.0)xpo=xpo/npom

      npk=0
      xpj=0d0
      do i=1,maproj
       if(xpp(i).lt.0.999)then
        xpj=xpj+xpp(i)!*xmp(i)
        npk=npk+1
       endif
      enddo
      if(npk.gt.0)xpj=xpj/dble(npk)

      ntk=0
      xtg=0d0
      do j=1,matarg
       if(xmt(j).lt.0.999)then
        xtg=xtg+xmt(j)!*xpt(j)
        ntk=ntk+1
       endif
      enddo
      if(ntk.gt.0)xtg=xtg/dble(ntk)

      x(kc)=xpo
      xp(kc)=xpj
      xt(kc)=xtg

        elseif(iii.eq.2)then

      nstat=nfull
      if(nstat.lt.0)nstat=max(1,nevent)
      x1=0
      x2=nbim

      write(ifhi,'(a)')       '!##################################'
      write(ifhi,'(a,i3)')    '!   average x  Pom   for event ',nrevt+1
      write(ifhi,'(a)')       '!##################################'
      write(ifhi,'(a,i1)')    'openhisto name avxPom-',nrevt+1
      write(ifhi,'(a)')       'htyp lin'
      write(ifhi,'(a)')       'xmod lin ymod lin'
      write(ifhi,'(a,2e11.3)')'xrange',x1,x2
      write(ifhi,'(a)')    'text 0 0 "xaxis iteration"'
      write(ifhi,'(a)')    'text 0 0 "yaxis average x Pomeron"'
      write(ifhi,'(a,e22.14)')'histoweight ',dble(nstat)
      write(ifhi,'(a)')       ' array 2'
         do k=0,nbim
      write(ifhi,'(2e11.3)')float(k),x(k)
         enddo
      write(ifhi,'(a)')    '  endarray'
      write(ifhi,'(a)')    'closehisto plot 0'

      write(ifhi,'(a)')       '!##################################'
      write(ifhi,'(a,i3)')    '!   average x proj   for event ',nrevt+1
      write(ifhi,'(a)')       '!##################################'
      write(ifhi,'(a,i1)')    'openhisto name avxProj-',nrevt+1
      write(ifhi,'(a)')       'htyp lin'
      write(ifhi,'(a)')       'xmod lin ymod lin'
      write(ifhi,'(a,2e11.3)')'xrange',x1,x2
      write(ifhi,'(a)')    'text 0 0 "xaxis iteration"'
      write(ifhi,'(a)')    'text 0 0 "yaxis average x proj"'
      write(ifhi,'(a,e22.14)')'histoweight ',dble(nstat)
      write(ifhi,'(a)')       ' array 2'
         do k=0,nbim
      write(ifhi,'(2e11.3)')float(k),xp(k)
         enddo
      write(ifhi,'(a)')    '  endarray'
      write(ifhi,'(a)')    'closehisto plot 0'

      write(ifhi,'(a)')       '!##################################'
      write(ifhi,'(a,i3)')    '!   average x targ   for event ',nrevt+1
      write(ifhi,'(a)')       '!##################################'
      write(ifhi,'(a,i1)')    'openhisto name avxTarg-',nrevt+1
      write(ifhi,'(a)')       'htyp lin'
      write(ifhi,'(a)')       'xmod lin ymod lin'
      write(ifhi,'(a,2e11.3)')'xrange',x1,x2
      write(ifhi,'(a)')    'text 0 0 "xaxis iteration"'
      write(ifhi,'(a)')    'text 0 0 "yaxis average x targ"'
      write(ifhi,'(a,e22.14)')'histoweight ',dble(nstat)
      write(ifhi,'(a)')       ' array 2'
         do k=0,nbim
      write(ifhi,'(2e11.3)')float(k),xt(k)
         enddo
      write(ifhi,'(a)')    '  endarray'
      write(ifhi,'(a)')    'closehisto plot 0'
        endif

      return
      end

c-----------------------------------------------------------------------
      subroutine xEmsRx(iii,id,xp,xm)
c-----------------------------------------------------------------------
c plot  x+, x-, x, y distribution of remnants
c-----------------------------------------------------------------------

#include "aaa.h"
#include "ems.h"
#include "par.h"

      parameter(nbix=50,nbiy=50,nid=2)
      double precision xu,xo,xpu,xpo,xmu,xmo,x,xint
      common/cxp/nxp(nid),nxm(nid),nx(nid),ny(nid)
     *,wxp(nbix,nid),wxm(nbix,nid),wx(nbix,nid),wy(nbiy,nid)
     *,xpu,xpo,xmu,xmo,xu,xo,yu,yo,dy
      logical lsym

      if(iemsrx.eq.0)call utstop('ERROR in XemsRx: iemsrx = 0&')

        if(iii.eq.0)then

      xpu=1./engy**2
      xpo=1.
      xmu=1./engy**2
      xmo=1.
      xu=1./engy**2
      xo=1.
      yu=-log(engy**2)
      yo=log(engy**2)
      dy=(yo-yu)/nbiy
      do j=1,nid
       nxp(j)=0
       nxm(j)=0
       nx(j)=0
       do i=1,nbix
        wxp(i,j)=0
        wxm(i,j)=0
        wx(i,j)=0
       enddo
       ny(j)=0
       do i=1,nbiy
        wy(i,j)=0
       enddo
      enddo

        elseif(iii.eq.1)then

      i=0
      if(xp.lt.xpu)goto 1
      i=1+int(log(xp/xpu)/log(xpo/xpu)*nbix)
      if(i.gt.nbix)goto 1
      if(i.lt.1)goto 1
      wxp(i,id)=wxp(i,id)+1
      nxp(id)=nxp(id)+1
1     continue

      if(xm.lt.xmu)goto 2
      i=1+int(log(xm/xmu)/log(xmo/xmu)*nbix)
      if(i.gt.nbix)goto 2
      if(i.lt.1)goto 2
      wxm(i,id)=wxm(i,id)+1
      nxm(id)=nxm(id)+1
2     continue

      x=xp*xm
      if(x.lt.xu)goto 3
      i=1+int(log(x/xu)/log(xo/xu)*nbix)
      if(i.gt.nbix)goto 3
      if(i.lt.1)goto 3
      wx(i,id)=wx(i,id)+1
      nx(id)=nx(id)+1
3     continue

      if(xm.le.0.)goto 4
      if(xp.le.0.)goto 4
      y=0.5*log(xp/xm)
c      if(id.eq.1)write(ifch,'(4e13.5,$)')log10(x),log10(xp),log10(xm),y
c      if(id.eq.2)write(ifch,'(4e13.5)')log10(x),log10(xp),log10(xm),y
      if(y.lt.yu)goto 4
      i=int((y-yu)/dy)+1
      if(i.gt.nbiy)goto 4
      if(i.lt.1)goto4
      wy(i,id)=wy(i,id)+1
      ny(id)=ny(id)+1
4     continue

        elseif(iii.eq.2)then
          
      lsym=maproj.eq.1.and.matarg.eq.1.and.idproj.eq.idtarg
      nstat=nfull
      if(nstat.lt.0)nstat=max(1,nevent)

      do j=1,nid
      if(j.eq.1)then
        iclrem=iclpro
      else
        iclrem=icltar
      endif
      write(ifhi,'(a)')      '!----------------------------------'
      write(ifhi,'(a)')      '!   remnant xp distribution      '
      write(ifhi,'(a)')      '!----------------------------------'
      write(ifhi,'(a,i1)')    'openhisto name xpRemnant-',j
      write(ifhi,'(a)')       'htyp lin'
      write(ifhi,'(a)')       'xmod log ymod log'
      write(ifhi,'(a,2e11.3)')'xrange',xpu,xpo
      write(ifhi,'(a)')    'text 0 0 "xaxis remnant x+"'
      write(ifhi,'(a)')    'text 0 0 "yaxis P(x+)"'
      if(j.eq.1)then
              write(ifhi,'(a)')'text 0.1 0.35 "Projectile"'
      if(lsym)write(ifhi,'(a)')'text 0.1 0.25 "High mass"'
      endif
      if(j.eq.2)then
              write(ifhi,'(a)')'text 0.1 0.35 "Target"'
      if(lsym)write(ifhi,'(a)')'text 0.1 0.25 "Low mass"'
      endif
      write(ifhi,'(a,e22.14)')'histoweight ',dble(nstat)
      write(ifhi,'(a)')       ' array 2'
         do i=1,nbix
      x=xpu*(xpo/xpu)**((i-0.5)/nbix)
      dx=xpu*(xpo/xpu)**(1.*i/nbix)*(1.-(xpo/xpu)**(-1./nbix))
      if(nxp(j).ne.0)write(ifhi,'(2e11.3)')x,wxp(i,j)/dx/nxp(j)
      if(nxp(j).eq.0)write(ifhi,'(2e11.3)')x,0.
         enddo
      write(ifhi,'(a)')    '  endarray'
      write(ifhi,'(a)')    'closehisto plot 0-'
      write(ifhi,'(a)')       'openhisto'
      write(ifhi,'(a)')       'htyp lin'
      write(ifhi,'(a)')       'array 2'
         do i=1,nbix
      x=xu*(xo/xu)**((i-0.5)/nbix)
      write(ifhi,'(2e11.3)')x,x**alplea(iclrem)*(1+alplea(iclrem))
         enddo
      write(ifhi,'(a)')    '  endarray'
      write(ifhi,'(a)')    'closehisto plot 0'

      write(ifhi,'(a)')      '!----------------------------------'
      write(ifhi,'(a)')      '!   remnant xm distribution      '
      write(ifhi,'(a)')      '!----------------------------------'
      write(ifhi,'(a,i1)')    'openhisto name xmRemnant-',j
      write(ifhi,'(a)')       'htyp lin'
      write(ifhi,'(a)')       'xmod log ymod log'
      write(ifhi,'(a,2e11.3)')'xrange',xmu,xmo
      write(ifhi,'(a)')    'text 0 0 "xaxis remnant x-"'
      write(ifhi,'(a)')    'text 0 0 "yaxis P(x-)"'
      if(j.eq.1)then
              write(ifhi,'(a)')'text 0.1 0.35 "Projectile"'
      if(lsym)write(ifhi,'(a)')'text 0.1 0.25 "High mass"'
      endif
      if(j.eq.2)then
              write(ifhi,'(a)')'text 0.1 0.35 "Target"'
      if(lsym)write(ifhi,'(a)')'text 0.1 0.25 "Low mass"'
      endif
      write(ifhi,'(a,e22.14)')'histoweight ',dble(nstat)
      write(ifhi,'(a)')       ' array 2'
         do i=1,nbix
      x=xmu*(xmo/xmu)**((i-0.5)/nbix)
      dx=xmu*(xmo/xmu)**(1.*i/nbix)*(1.-(xmo/xmu)**(-1./nbix))
      if(nxm(j).ne.0)write(ifhi,'(2e11.3)')x,wxm(i,j)/dx/nxm(j)
      if(nxm(j).eq.0)write(ifhi,'(2e11.3)')x,0.
         enddo
      write(ifhi,'(a)')    '  endarray'
      write(ifhi,'(a)')    'closehisto plot 0'

      write(ifhi,'(a)')      '!----------------------------------'
      write(ifhi,'(a)')      '!   remnant x distribution      '
      write(ifhi,'(a)')      '!----------------------------------'
      write(ifhi,'(a,i1)')    'openhisto name xRemnant-',j
      write(ifhi,'(a)')       'htyp lin'
      write(ifhi,'(a)')       'xmod log ymod log'
      write(ifhi,'(a,2e11.3)')'xrange',xu,xo
      write(ifhi,'(a)')    'text 0 0 "xaxis remnant x"'
      write(ifhi,'(a)')    'text 0 0 "yaxis P(x)"'
      if(j.eq.1)then
              write(ifhi,'(a)')'text 0.1 0.35 "Projectile"'
      if(lsym)write(ifhi,'(a)')'text 0.1 0.25 "High mass"'
      endif
      if(j.eq.2)then
              write(ifhi,'(a)')'text 0.1 0.35 "Target"'
      if(lsym)write(ifhi,'(a)')'text 0.1 0.25 "Low mass"'
      endif
      write(ifhi,'(a,e22.14)')'histoweight ',dble(nstat)
      write(ifhi,'(a)')       ' array 2'
         do i=1,nbix
      x=xu*(xo/xu)**((i-0.5)/nbix)
      dx=xu*(xo/xu)**(1.*i/nbix)*(1.-(xo/xu)**(-1./nbix))
      if(nx(j).ne.0)write(ifhi,'(2e11.3)')x,wx(i,j)/dx/nx(j)
      if(nx(j).eq.0)write(ifhi,'(2e11.3)')x,0.
         enddo
      write(ifhi,'(a)')    '  endarray'
      write(ifhi,'(a)')    'closehisto plot 0-'
      write(ifhi,'(a)')       'openhisto'
      write(ifhi,'(a)')       'htyp lga'
      write(ifhi,'(a)')       'array 2'
      if(j.eq.1)then
        xint=(1.-alphigh(2))
        if(xint.eq.0.)then
          xint=-log(xu)
        else
          xint=(1d0-xu**xint)/xint
        endif
      else
        xint=(1.-alphigh(2))
        if(xint.eq.0.)then
          xint=-log(xu)
        else
          xint=(1d0-xu**xint)/xint
        endif
      endif
         do i=1,nbix
      x=xu*(xo/xu)**((i-0.5)/nbix)
      if(j.eq.1)write(ifhi,'(2e11.3)')x
     .         ,x**(-alphigh(2))/xint
      if(j.eq.2)write(ifhi,'(2e11.3)')x
     .         ,x**(-alphigh(2))/xint
         enddo
      write(ifhi,'(a)')    '  endarray'
      write(ifhi,'(a)')    'closehisto plot 0-'
      write(ifhi,'(a)')       'openhisto'
      write(ifhi,'(a)')       'htyp lbu'
      write(ifhi,'(a)')       'array 2'
      if(j.eq.1)then
        xint=(1.-alppom)
        if(lsym)xint=(1.-alpdi(iclpro))
        if(xint.eq.0.)then
          xint=-log(xu)
        else
          xint=(1.-xu**xint)/xint
        endif
      else
        xint=(1-alppom)
        if(xint.eq.0.)then
          xint=-log(xu)
        else
          xint=(1.-xu**xint)/xint
        endif
      endif
         do i=1,nbix
      x=xu*(xo/xu)**((i-0.5)/nbix)
      if(j.eq.1)then
        if(lsym)then
               write(ifhi,'(2e11.3)')x
     .         ,x**(-alpdi(iclpro))/xint
        else
               write(ifhi,'(2e11.3)')x
     .         ,x**(-alppom)/xint
        endif
      endif
      if(j.eq.2)write(ifhi,'(2e11.3)')x
     .         ,x**(-alppom)/xint
         enddo
      write(ifhi,'(a)')    '  endarray'
      write(ifhi,'(a)')    'closehisto plot 0'

      write(ifhi,'(a)')      '!----------------------------------'
      write(ifhi,'(a)')      '!   remnant y distribution      '
      write(ifhi,'(a)')      '!----------------------------------'
      write(ifhi,'(a,i1)')    'openhisto name yRemnant-',j
      write(ifhi,'(a)')       'htyp lin'
      write(ifhi,'(a)')       'xmod lin ymod log'
      write(ifhi,'(a,2e11.3)')'xrange',yu,yo
      write(ifhi,'(a)')    'text 0 0 "xaxis remnant y"'
      write(ifhi,'(a)')    'text 0 0 "yaxis P(y)"'
      if(j.eq.1)then
              write(ifhi,'(a)')'text 0.1 0.35 "Projectile"'
      if(lsym)write(ifhi,'(a)')'text 0.1 0.25 "High mass"'
      endif
      if(j.eq.2)then
              write(ifhi,'(a)')'text 0.1 0.35 "Target"'
      if(lsym)write(ifhi,'(a)')'text 0.1 0.25 "Low mass"'
      endif
      write(ifhi,'(a,e22.14)')'histoweight ',dble(nstat)
      write(ifhi,'(a)')       ' array 2'
         do i=1,nbix
      y=yu+dy/2.+(i-1)*dy
      if(ny(j).ne.0)write(ifhi,'(2e11.3)')y,wy(i,j)/dy/ny(j)
      if(ny(j).eq.0)write(ifhi,'(2e11.3)')y,0.
         enddo
      write(ifhi,'(a)')    '  endarray'
      write(ifhi,'(a)')    'closehisto plot 0'

      enddo

        endif

      return
      end

c-----------------------------------------------------------------------
      subroutine xEmsPm(iii,ko,nmci,nmcmx)
c-----------------------------------------------------------------------
c m (pomeron number) distribution for different b-bins.
c arguments:
c   iii:  modus (0,1,2)
c   ko:   pair number (1 - AB)
c   nmc:  number of pomerons
c   nmcmx: number max of pomerons
c-----------------------------------------------------------------------
#include "aaa.h"
#include "ems.h"
#include "par.h"
      common/geom/rmproj,rmtarg,bmax,bkmx
      parameter(nbin=200)
      parameter(nbib=32)
      common/cn/wn(0:nbin,nbib),wnmc(0:nbin,nbib),npmx(nbib),nn(nbib)
     &         ,nn2(nbib),dn(nbib)
      common/cb1/db,b1,b2,bb(nbib),nbibx
      double precision plc,s,om1intbc
      character ce*8,cb*4
      common/cems5/plc,s
      common/cemspm/sumb(nbib)

      if(iemspm.eq.0)call utstop('ERROR in XemsPm: iemspm = 0&')

        if(iii.eq.0)then

      do k=1,nbib
       nn(k)=0
       nn2(k)=0
       sumb(k)=0
       dn(k)=1.
       do i=0,nbin
        wnmc(i,k)=0
       enddo
      enddo
      nbibx=6
      b1=0
      b2=2
      db=(b2-b1)/nbibx


        elseif(iii.eq.1)then

      k=int((bk(ko)-b1)/db)+1
c      nmc=nmci
      if(k.gt.nbibx)k=nbibx
      if(k.lt.1)k=1
      dn(k)=max(1.,float(nmcmx)/float(nbin))
      nmc=nint(float(nmci)/dn(k)+0.499999)
      if(nmc.gt.nbin)nmc=nbin
      if(nmc.lt.0)return
      nn(k)=nn(k)+1
      wnmc(nmc,k)=wnmc(nmc,k)+1./dn(k)
      sumb(k)=sumb(k)+bk(ko)


        elseif(iii.eq.2)then

      nstat=nfull
      if(nstat.lt.0)nstat=max(1,nevent)
      kollini=koll
      koll=1                    !to have screening for pp
      iomegasave=iomega
      call DefXminDf(s)

      do 1 k=1,nbibx

       bb(k)=b1+(k-0.5)*db
       if(maproj.eq.1.and.matarg.eq.1.and.bmaxim.eq.0.)bb(k)=b1
       om1i=sngl(om1intbc(bb(k)))
       wntmp=0.
       do 10 i=0,nbin
         wn(i,k)=0.
         if(wntmp.gt.1e5)goto 10
         do j=i,i+int(dn(k))-1
           if(j.eq.0)then
             wntmp=exp(-om1i)
           else
             wntmp=wntmp*om1i/j
           endif
           wn(i,k)=wn(i,k)+wntmp/dn(k)
         enddo
         if(wn(i,k).gt.0.000001*(1.-exp(-om1i)))npmx(k)=i
 10    continue

      write(ifhi,'(a)')   '!##################################'
      write(ifhi,'(a)')   '! distr of Pomeron number vs b'
      write(ifhi,'(a)')   '!##################################'
      write(ce,'(f8.2)')sngl(plc)
      write(cb,'(f4.2)')bb(k)
c      if(nn(k).gt.0)then
      write(ifhi,'(a,i1)')    'openhisto name mPom-',k
      write(ifhi,'(a)')       'htyp lru'
      write(ifhi,'(a)')       'xmod lin ymod log'
      write(ifhi,'(a,2e11.3)')'xrange',0.,float(npmx(k))*dn(k)
      write(ifhi,'(a)')    'text 0 0 "xaxis number m of Pomerons"'
      write(ifhi,'(a)')    'text 0 0 "yaxis prob(m)"'
      if(k.eq.1)
     *write(ifhi,'(a,a)')     'text 0.5 0.90 "E ='//ce//'"'
      write(ifhi,'(a,a)')     'text 0.5 0.80 "b ='//cb//'"'
      write(ifhi,'(a,e22.14)')'histoweight ',dble(nstat)
      write(ifhi,'(a)')       ' array 2'
         do i=0,nbin
      write(ifhi,'(2e11.3)')float(i)*dn(k),wnmc(i,k)/max(1,nn(k))
         enddo
      write(ifhi,'(a)')    '  endarray'
      write(ifhi,'(a)')    'closehisto plot 0-'
c      endif

      write(ifhi,'(a)')   '!##################################'
      write(ifhi,'(a)')   '! distr of Pomeron number vs b'
      write(ifhi,'(a)')   '!   traditional approach'
      write(ifhi,'(a)')   '!##################################'
      write(ifhi,'(a,i1)')    'openhisto name mPomTradi-',k
      write(ifhi,'(a)')       'htyp lba'
      write(ifhi,'(a)')       'xmod lin ymod log'
      write(ifhi,'(a,2e11.3)')'xrange',0.,float(npmx(k))*dn(k)
      write(ifhi,'(a)')       ' array 2'
         do i=0,nbin
      write(ifhi,'(2e11.3)')float(i)*dn(k),wn(i,k)
         enddo
      write(ifhi,'(a)')    '  endarray'
      write(ifhi,'(a)')    'closehisto plot 0'

 1    continue

      koll=kollini
      iomega=iomegasave
      call DefXminDf(s)

      endif

      return
      end

c-----------------------------------------------------------------------
      subroutine xEmsB(iii,jjj,ko)
c-----------------------------------------------------------------------
c b distribution at different stages
c arguments:
c   iii:  modus (0,1,2)
c   jjj:  stage or type of interaction
c     just after Metropolis:
c           1 ... all
c           2 ... interaction
c     after defining diffraction:
c           3 ... nothing
c           4 ... cut
c           5 ... diffr
c           6 ... diffr cut
c   ko:   pair number (1 - AB)
c-----------------------------------------------------------------------
#include "aaa.h"
#include "ems.h"
#include "sem.h"
#include "par.h"
      parameter(njjj=6)
      parameter(nbib=32)
      common/cxemsb1/w(0:njjj,nbib),nn(njjj)
      common/cxemsb2/db,b1,b2
      common/cxemsb3/njjj1
      double precision PhiExact,om1intby,HPhiDiff,PhiUnit
     &,HPhiInt
      common/geom/rmproj,rmtarg,bmax,bkmx
      dimension uua(nbib),uuo(nbib),uua2(nbib),uuo2(nbib),uu6(nbib)
     &,uu4(nbib),uu5(nbib)

      if(iemsb.eq.0)call utstop('ERROR in XemsB: iemsB = 0&')

        if(iii.eq.0)then

      do k=1,nbib
       do j=0,njjj
        w(j,k)=0
       enddo
      enddo
      do j=1,njjj
       nn(j)=0
      enddo
      njjj1=0

        elseif(iii.eq.1)then

      b1=0
      b2=bkmx*1.2
      db=(b2-b1)/nbib
      k=int((bk(ko)-b1)/db)+1
      if(k.gt.nbib)return
      if(k.lt.1)return
      w(jjj,k)=w(jjj,k)+1
      nn(jjj)=nn(jjj)+1
      if(jjj.eq.1)njjj1=1

        elseif(iii.eq.2)then

      if(njjj1.ne.1)call utstop
     &         ('xEmsB must be called also with jjj=1&')


      nstat=nfull
      if(nstat.lt.0)nstat=max(1,nevent)
      ymaxd=0
      kollini=koll
      koll=1
      fk=bkmx**2*pi
      do k=1,nbib
       x=b1+(k-0.5)*db
       y=fk*w(5,k)/nn(1)/(pi*((x+0.5*db)**2-(x-0.5*db)**2))
       ymaxd=max(ymaxd,y)
      enddo
      ymaxd=ymaxd+0.1
      ymax=1.4

      do 1 j=1,njjj
       !if(nn(j).eq.0)goto1

      write(ifhi,'(a)')   '!##################################'
      write(ifhi,'(a)')   '! b distr exact theory '
      write(ifhi,'(a)')   '!##################################'
         if(j.ge.2.and.j.le.6)then
      write(ifhi,'(a,i1,a)')  'openhisto name b',j,'Exact'
      write(ifhi,'(a)')       'htyp lba xmod lin ymod lin'
      write(ifhi,'(a,2e11.3)')'xrange',0.0,b2
      if(j.ge.5)write(ifhi,'(a,2e11.3)')'yrange',0.,ymaxd
      if(j.lt.5)write(ifhi,'(a,2e11.3)')'yrange',0.,ymax
      write(ifhi,'(a)')    'text 0 0 "xaxis impact parameter b"'
      write(ifhi,'(a)')    'text 0 0 "yaxis P(b)"'
      write(ifhi,'(a)')       'array 2'
         do k=1,nbib
      b=b1+(k-0.5)*db
      if(j.eq.2)then
c define variables for for om1intby, HPhiInt and PhiUnit
        do i=idxDmin(iomega),idxDmax(iomega)
          call Gfunpar(0.,0.,0.,0.,1,i,b,engy**2,xa,xb,xc,xd,xe,xf,xg)
          call Gfunpar(0.,0.,0.,0.,2,i,b,engy**2,xa,xb,xc,xd,xe,xf,xg)
        enddo
        uua(k)=sngl(HPhiInt(engy**2,b))
        uuo2(k)=sngl(PhiUnit(1.d0,1.d0))
        uuo(k)=1.-uuo2(k)
        uua2(k)=min(uuo2(k),max(0.,
     &          sngl(Phiexact(0.,0.,1.,1.d0,1.d0,engy**2,b))))
        uu6(k)=sngl(om1intby(engy**2,b,-50,1)
     &             -om1intby(engy**2,b,-90,1))
        uu4(k)=sngl(HPhiDiff(engy**2,b))
        uu5(k)=sngl(om1intby(engy**2,b,-50,1))
      endif
      if(j.eq.2)y=(1.-uua2(k))
      if(j.eq.3)y=uua2(k)
      if(j.eq.4)y=uua(k)-uu4(k)
      if(j.eq.5.or.j.eq.6)y=uu4(k)
      write(ifhi,'(2e11.3)')b,max(1.e-6,min(100.,y))
         enddo
      write(ifhi,'(a)')    '  endarray'
      write(ifhi,'(a)')    'closehisto plot 0-'
         endif
      write(ifhi,'(a)')   '!##################################'
      write(ifhi,'(a)')   '! b distr unitarized theory '
      write(ifhi,'(a)')   '!##################################'
      write(ifhi,'(a,i1,a)')  'openhisto name b',j,'Unit'
      write(ifhi,'(a)')       'htyp lbf xmod lin ymod lin'
      write(ifhi,'(a)')    'text 0 0 "xaxis impact parameter b"'
      write(ifhi,'(a)')    'text 0 0 "yaxis P(b)"'
      write(ifhi,'(a)')       'array 2'
         do k=1,nbib
      b=b1+(k-0.5)*db
      if(j.eq.1)y=1
      if(j.eq.2)y=(1.-uuo2(k))
      if(j.eq.3)y=uuo2(k)
      if(j.eq.4)y=uuo(k)-uu6(k)
      if(j.eq.5)y=uu5(k)
      if(j.eq.6)y=uu6(k)
      write(ifhi,'(2e11.3)')b,y
         enddo
      write(ifhi,'(a)')    '  endarray'
      write(ifhi,'(a)')    'closehisto plot 0-'
      write(ifhi,'(a)')   '!##################################'
      write(ifhi,'(a)')   '! b distr for cross section '
      write(ifhi,'(a)')   '!##################################'
      write(ifhi,'(a,i1,a)')  'openhisto name b',j,'Unit'
      write(ifhi,'(a)')       'htyp lge xmod lin ymod lin'
      write(ifhi,'(a)')    'text 0 0 "xaxis impact parameter b"'
      write(ifhi,'(a)')    'text 0 0 "yaxis P(b)"'
      write(ifhi,'(a)')       'array 2'
         do k=1,nbib
      b=b1+(k-0.5)*db
      if(j.eq.1)y=1
      if(j.eq.2)y=(1.-(uuo2(k)+uua2(k))*0.5)
      if(j.eq.3)y=(uuo2(k)+uua2(k))*0.5
      if(j.eq.4)y=(1.-(uuo2(k)+uua2(k))*0.5)-uu6(k)
      if(j.eq.5)y=uu5(k)
      if(j.eq.6)y=uu6(k)

c*************************************************************************************
      if(y.lt.-1e5)y=-99999 !kw   *****   otherwise hiku has problems       ****
c************************************************************************************* 

      write(ifhi,'(2e11.3)')b,y
         enddo
      write(ifhi,'(a)')    '  endarray'
      write(ifhi,'(a)')    'closehisto plot 0-'
      if(j.eq.1)then !cKW21 bug related to "add" fixed (in case of few events, like with hydro)
        write(ifhi,'(a)')   '!##################################'
        write(ifhi,'(a)')   '! b distribution simulation Norm'
        write(ifhi,'(a)')   '!##################################'
        write(ifhi,'(a)')     'openhisto name bNormSimu' 
        write(ifhi,'(a,e22.14)')'histoweight ',dble(nstat)
        write(ifhi,'(a)')       ' array 2'
        do k=1,nbib
          x=b1+(k-0.5)*db
          y=nn(1)*(pi*((x+0.5*db)**2-(x-0.5*db)**2))
          write(ifhi,'(2e11.3)')x,y
        enddo
        write(ifhi,'(a)')    '  endarray'
        write(ifhi,'(a)')    'closehisto'
      endif
      write(ifhi,'(a)')   '!##################################'
      write(ifhi,'(a)')   '! b distribution simulation '
      write(ifhi,'(a)')   '!##################################'
      write(ifhi,'(a,i1,a)')  'openhisto name b',j,'Simu'
      write(ifhi,'(a,e22.14)')'histoweight ',dble(nstat)
      write(ifhi,'(a)')       ' array 2'
      do k=1,nbib
        x=b1+(k-0.5)*db
        y=fk*w(j,k)
        write(ifhi,'(2e11.3)')x,y
      enddo
      write(ifhi,'(a)')    ' endarray'
      write(ifhi,'(a)')    'closehisto'
      write(ifhi,'(a,i1,a)')  'openhisto name b',j,'SimuNormalized'
      write(ifhi,'(a)')       'htyp lrf xmod lin ymod lin'
      write(ifhi,'(a,2e11.3)')'xrange',0.0,b2
      write(ifhi,'(a,2e11.3)')'yrange',0.,ymax
      write(ifhi,'(a)')    'text 0 0 "xaxis impact parameter b"'
      write(ifhi,'(a)')    'text 0 0 "yaxis P(b)"'
      if(j.eq.1)write(ifhi,'(a)')'text 0.1 0.35 "after Metropolis"'
      if(j.eq.1)write(ifhi,'(a)')'text 0.2 0.20 "all "'
      if(j.eq.2)write(ifhi,'(a)')'text 0.3 0.85 "after Metropolis"'
      if(j.eq.2)write(ifhi,'(a)')'text 0.5 0.70 "interaction "'
      if(j.eq.3)write(ifhi,'(a)')'text 0.3 0.85 "nothing"'
      if(j.eq.4)write(ifhi,'(a)')'text 0.3 0.85 "cut"'
      if(j.eq.5)write(ifhi,'(a)')'text 0.3 0.85 "diffr"'
      if(j.eq.6)write(ifhi,'(a)')'text 0.3 0.85 "diffr cut"'
      write(ifhi,'(a)')    'oh calc *1 / bNormSimu ; ch plot 0'
   1  continue

      koll=kollini

      endif

      return
      end

c-----------------------------------------------------------------------
      subroutine xEmsBg(iii,jjj,ko,no)
c-----------------------------------------------------------------------
c b distribution at different stages for different group
c arguments:
c   iii:  modus (0,1,2,3)
c   jjj:  group of interaction (1,2 ... ,7)
c   ko:   pair number (1 - AB)
c   no:   Pomeron index
c-----------------------------------------------------------------------
#include "aaa.h"
#include "ems.h"
#include "par.h"
      parameter(njjj=10)
      parameter(nbib=16)
      common/cxemsb4/wg(-1:njjj,nbib),nng(nbib),nng2(2,nbib),uug(nbib)
     .              ,kollx
      common/cxemsb5/dbg,b1g,b2g,b3g,dbg2
      common/cxemsb6/njjj0
      double precision seedp,PhiExpo!,PhiExact
      common/nucl3/phi,bimp
      common/geom/rmproj,rmtarg,bmax,bkmx
      double precision plc,s
      common/cems5/plc,s

      if(iemsbg.eq.0)call utstop('ERROR in XemsBg: iemsbg = 0&')

c      write(ifch,*)'xEmsBg',iii,jjj,ko,no

        if(iii.eq.0)then

      do k=1,nbib
       nng(k)=0
       nng2(1,k)=0
       nng2(2,k)=0
       do j=-1,njjj
        wg(j,k)=0
       enddo
      enddo
      njjj0=0
      kollx=0

        elseif(iii.eq.1)then

      b1g=0
          if(jjj.le.8)then
      b2g=bkmx*1.2
      dbg=(b2g-b1g)/nbib
      k=int((bk(ko)-b1g)/dbg)+1
          else
      b3g=bmax
      dbg2=(b3g-b1g)/nbib
      k=int((bimp-b1g)/dbg2)+1
          endif
      if(k.gt.nbib)return
      if(k.lt.1)return
      if(jjj.eq.-1.or.jjj.eq.0)then
        wg(jjj,k)=wg(jjj,k)+1
      elseif(jjj.ge.9)then
        wg(jjj,k)=wg(jjj,k)+q2kmin(jjj-8,no,ko)
        if(jjj.eq.9)nng2(1,k)=nng2(1,k)+1
      elseif(jjj.ge.7)then
        wg(jjj,k)=wg(jjj,k)+q2kmin(jjj-6,no,ko)
        if(jjj.eq.7)nng2(2,k)=nng2(2,k)+1
      else
        wg(jjj,k)=wg(jjj,k)+1
        nng(k)=nng(k)+1
      endif
      if(jjj.eq.0)njjj0=1

        elseif(iii.eq.3)then

          iomegasave=iomega
          iomega=2
          call DefXminDf(s)
          call ranfgt(seedp)
          do k=1,koll
            om1i=sngl(om1intc(k))
            if(rangen().lt.1.-exp(-om1i))then
c            om1i=sngl(PhiExpo(0.,0.,1.,1.d0,1.d0,engy*engy,bk(k)))
c            if(rangen().lt.1.-om1i)then
              kollx=kollx+1
            endif
          enddo
          call ranfst(seedp)
          iomega=iomegasave
          call DefXminDf(s)

        elseif(iii.eq.2)then

      if(njjj0.ne.1)call utstop
     &('xEmsBg must be called also with jjj=0&')
      ymax=1.4
      kollini=koll
      koll=1
      nstat=nfull
      if(nstat.lt.0)nstat=max(1,nevent)

      wtot=1.
      if(matarg+maproj.gt.2)then
      wtot=0.
      do k=1,nbib
       wtot=wtot+wg(-1,k)
      enddo
      if(kollx.gt.0)wtot=wtot/float(kollx)
      endif

      do 1 j=1,6

      write(ifhi,'(a)')   '!##################################'
      write(ifhi,'(a)')   '! b distribution simulation'
      write(ifhi,'(a)')   '!##################################'
      write(ifhi,'(a,i1,a)')  'openhisto name bg',j,'Simu'
      write(ifhi,'(a)')       'htyp lin xmod lin ymod log'
      write(ifhi,'(a,2e11.3)')'xrange',0.,b2g
      write(ifhi,'(a,2e11.3)')'yrange',1.e-5,ymax*10.
      write(ifhi,'(a)')    'text 0 0 "xaxis impact parameter b"'
      write(ifhi,'(a)')    'text 0 0 "yaxis P(b)"'
      if(wtot.gt.0.d0)
     &write(ifhi,'(a,f7.4,a)')    'text 0.5 0.8 "alpha=',1./wtot,'"'
      write(ifhi,'(a,e22.14)')'histoweight ',dble(nstat)
      write(ifhi,'(a)')       ' array 2'
         do k=1,nbib
      b=b1g+(k-0.5)*dbg
      y=0.
      if(nng(k).ne.0.and.wg(0,k).ne.0)
     &              y=wg(j,k)/float(nng(k))*wg(-1,k)/wg(0,k)!/wtot
c      if(wg(0,k).ne.0..and.nng(k).ne.0)y=wg(j,k)/nng(k)*wg(-1,k)/wg(0,k)
c!???????????? better normalization ? probability to have an interaction
c in epos compared to eikonal probability, instead of normalized by the
c probability of a collision for a pair (the number collision/number
c active pair).
      uug(k)=uug(k)+y
      write(ifhi,'(2e11.3)')b,y
         enddo
      write(ifhi,'(a)')    '  endarray'
      write(ifhi,'(a)')    'closehisto plot 0-'
   1  continue
      write(ifhi,'(a)')   '!##################################'
      write(ifhi,'(a)')   '! b distr tot simul theory '
      write(ifhi,'(a)')   '!##################################'
      write(ifhi,'(a)')  'openhisto name btotSimu'
      write(ifhi,'(a)')       'htyp lfi xmod lin ymod lin'
      write(ifhi,'(a)')    'text 0 0 "xaxis impact parameter b"'
      write(ifhi,'(a)')    'text 0 0 "yaxis P(b)"'
      write(ifhi,'(a)')       'array 2'
         do k=1,nbib
      b=b1g+(k-0.5)*dbg
      write(ifhi,'(2e11.3)')b,uug(k)
         enddo
      write(ifhi,'(a)')    '  endarray'
      write(ifhi,'(a)')    'closehisto plot 0-'
      write(ifhi,'(a)')   '!##################################'
      write(ifhi,'(a)')   '! b distr unitarized theory '
      write(ifhi,'(a)')   '!##################################'
      write(ifhi,'(a,i1,a)')  'openhisto name bg',j,'Unit'
      write(ifhi,'(a)')       'htyp lba xmod lin ymod lin'
      write(ifhi,'(a)')    'text 0 0 "xaxis impact parameter b"'
      write(ifhi,'(a)')    'text 0 0 "yaxis P(b)"'
      write(ifhi,'(a)')       'array 2'
         do k=1,nbib
      b=b1g+(k-0.5)*dbg
      a1=sngl(PhiExpo(0.,0.,1.,1.d0,1.d0,engy**2,b))
      y=(1.-a1)
      write(ifhi,'(2e11.3)')b,y
         enddo
      write(ifhi,'(a)')    '  endarray'
      write(ifhi,'(a)')    'closehisto plot 0'

      do 2 j=9,10

      write(ifhi,'(a)')   '!##################################'
      write(ifhi,'(a)')   '! b distribution simulation'
      write(ifhi,'(a)')   '!##################################'
      write(ifhi,'(a,i1,a)')  'openhisto name q2pmin',j-8,'Simu'
      if(j.eq.10)then
      write(ifhi,'(a)')       'htyp lba xmod lin ymod lin'
      else
      write(ifhi,'(a)')       'htyp lri xmod lin ymod lin'
      endif
      write(ifhi,'(a,2e11.3)')'xrange',0.,b3g
      write(ifhi,'(a,e11.3,a)')'yrange',0.,' auto'
      write(ifhi,'(a)')    'text 0 0 "xaxis event impact parameter b"'
      write(ifhi,'(a)')    'text 0 0 "yaxis Q^2!?s!(b)"'
      write(ifhi,'(a,e22.14)')'histoweight ',dble(nstat)
      write(ifhi,'(a)')       ' array 2'
         do k=1,nbib
      b=b1g+(k-0.5)*dbg2
      y=0.
      if(nng2(1,k).ne.0)
     &              y=wg(j,k)/float(nng2(1,k))
      write(ifhi,'(2e11.3)')b,y
         enddo
      write(ifhi,'(a)')    '  endarray'
      write(ifhi,'(a)')    'closehisto plot 0-'
 2    continue

      do 3 j=7,8

      write(ifhi,'(a)')   '!##################################'
      write(ifhi,'(a)')   '! b distribution simulation'
      write(ifhi,'(a)')   '!##################################'
      write(ifhi,'(a,i1,a)')  'openhisto name q2pmin',j-6,'SimuP'
      if(j.eq.8)then
      write(ifhi,'(a)')       'htyp lbb xmod lin ymod lin'
      else
      write(ifhi,'(a)')       'htyp lrj xmod lin ymod lin'
      endif
      write(ifhi,'(a,2e11.3)')'xrange',0.,b2g
      write(ifhi,'(a,e11.3,a)')'yrange',0.,' auto'
      write(ifhi,'(a)')    'text 0 0 "xaxis impact parameter b"'
      write(ifhi,'(a)')    'text 0 0 "yaxis Q^2!?s!(b)"'
      write(ifhi,'(a,e22.14)')'histoweight ',dble(nstat)
      write(ifhi,'(a)')       ' array 2'
         do k=1,nbib
      b=b1g+(k-0.5)*dbg
      y=0.
      if(nng2(2,k).ne.0)
     &              y=wg(j,k)/float(nng2(2,k))
      write(ifhi,'(2e11.3)')b,y
         enddo
      write(ifhi,'(a)')    '  endarray'
      if(j.eq.7)then
        write(ifhi,'(a)')    'closehisto plot 0-'
      else
        write(ifhi,'(a)')    'closehisto plot 0'
      endif
 3    continue

      koll=kollini

      endif

      return
      end
      


c-----------------------------------------------------------------------
      subroutine xEmsSe(iii,xmc,ptmc,ih,iqq)
c-----------------------------------------------------------------------
c     iqq = 1 : String End mass and rapidity
c     iqq = 2 : String mass and rapidity
c-----------------------------------------------------------------------

#include "aaa.h"

      parameter(nbix=50)
      common/cxpar/nx(2),x(nbix),wxmc(nbix,2),xmn,xmx,xu,xo
      parameter(nbiy=40)
      common/cypar/ny(2),y(nbiy),wymc(nbiy,2),ymin,ymax,dy,yu,yo

      s=engy**2

      if(iii.eq.0)then

       nx(iqq)=0
       xu=0.1/engy**2
       xo=1.
       do i=1,nbix
         x(i)=xu*(xo/xu)**((i-0.5)/nbix)
         wxmc(i,iqq)=0
       enddo
       yo=log(s)
       yu=-yo
       dy=(yo-yu)/nbiy
       ny(iqq)=0
       do i=1,nbiy
         y(i)=yu+dy/2.+(i-1)*dy
         wymc(i,iqq)=0
       enddo

      elseif(iii.eq.1)then

       if(xmc.lt.xu)return
       if(ptmc.eq.0.)return
       if(iqq.eq.1)then
         ymc=0.5*log(xmc*s/ptmc)*ih
       else !if(iqq.eq.2)
         ymc=0.5*log(xmc/ptmc)
       endif
       i=1+int(log(xmc/xu)/log(xo/xu)*nbix)
       if(i.gt.nbix)goto1
       if(i.lt.1)goto1
       wxmc(i,iqq)=wxmc(i,iqq)+1
       nx(iqq)=nx(iqq)+1
1      continue
       if(ymc.lt.yu)return
       i=int((ymc-yu)/dy)+1
       if(i.gt.nbiy)return
       if(i.lt.1)return
       wymc(i,iqq)=wymc(i,iqq)+1
       ny(iqq)=ny(iqq)+1

      elseif(iii.eq.2)then

        nstat=nfull
        if(nstat.lt.0)nstat=max(1,nevent)

        write(ifhi,'(a)')        '!--------------------------------'
        write(ifhi,'(a)')        '!   string end x distr       '
        write(ifhi,'(a)')        '!--------------------------------'
        write(ifhi,'(a)')       'openhisto'
        write(ifhi,'(a)')       'htyp lin'
        write(ifhi,'(a)')       'xmod log ymod log'
        write(ifhi,'(a,2e11.3)')'xrange',xu,xo
        if(iqq.eq.1)write(ifhi,'(a)')    'text 0 0 "xaxis string end x"'
        if(iqq.eq.2)write(ifhi,'(a)')    'text 0 0 "xaxis string x"'
        write(ifhi,'(a)')    'text 0 0 "yaxis P(x)"'
        write(ifhi,'(a,e22.14)')'histoweight ',dble(nstat)
        write(ifhi,'(a)')       ' array 2'
        do i=1,nbix
         dx=xu*(xo/xu)**(1.*i/nbix)*(1.-(xo/xu)**(-1./nbix))
         if(nx(iqq).gt.0)
     *   write(ifhi,'(2e11.3)')x(i),wxmc(i,iqq)/dx/nx(iqq)
        enddo
        write(ifhi,'(a)')    '  endarray'
        write(ifhi,'(a)')    'closehisto plot 0'
        write(ifhi,'(a)')       'openhisto'
        write(ifhi,'(a)')       'htyp lin'
        write(ifhi,'(a)')       'xmod lin ymod lin'
        write(ifhi,'(a,2e11.3)')'xrange',yu,yo
        if(iqq.eq.1)write(ifhi,'(a)')    'text 0 0 "xaxis string end y"'
        if(iqq.eq.2)write(ifhi,'(a)')    'text 0 0 "xaxis string y"'
        write(ifhi,'(a)')    'text 0 0 "yaxis P(y)"'
        write(ifhi,'(a,e22.14)')'histoweight ',dble(nstat)
        write(ifhi,'(a)')       ' array 2'
        do i=1,nbiy
         if(ny(iqq).gt.0)
     *   write(ifhi,'(2e11.3)')y(i),wymc(i,iqq)/dy/ny(iqq)
        enddo
        write(ifhi,'(a)')    '  endarray'
        write(ifhi,'(a)')    'closehisto plot 0'
      endif

      return
      end

c-----------------------------------------------------------------------
      subroutine xEmsDr(iii,xpmc,xmmc,ie)
c-----------------------------------------------------------------------

#include "aaa.h"

      parameter(nbix=50,nie=4)
      common/cxpardr/nxp(nie),nxm(nie),x(nbix),wxpmc(nbix,nie)
     &      ,wxmmc(nbix,nie),xmn,xmx,xu,xo,wxmc(nbix,nie),nx(nie)
      parameter(nbiy=40)
      common/cypardr/ny(nie),y(nbiy),wymc(nbiy,nie),ymin,ymax,dy,yu,yo

      s=engy**2

      if(iii.eq.0)then

       do ni=1,nie
         nxp(ni)=0
         nxm(ni)=0
         nx(ni)=0
       enddo
       xu=0.1/engy**2
       xo=1.
       do i=1,nbix
         x(i)=xu*(xo/xu)**((i-0.5)/nbix)
         do ni=1,nie
           wxpmc(i,ni)=0
           wxmmc(i,ni)=0
           wxmc(i,ni)=0
         enddo
       enddo
       yo=log(s)
       yu=-yo
       dy=(yo-yu)/nbiy
       do ni=1,nie
         ny(ni)=0
       enddo
       do i=1,nbiy
         y(i)=yu+dy/2.+(i-1)*dy
         do ni=1,nie
           wymc(i,ni)=0
         enddo
       enddo

      elseif(iii.eq.1)then

       if(ie.lt.1.or.ie.gt.nie)return

       if(xpmc.lt.xu)return
       i=1+int(log(xpmc/xu)/log(xo/xu)*nbix)
       if(i.gt.nbix)goto1
       if(i.lt.1)goto1
       wxpmc(i,ie)=wxpmc(i,ie)+1
       nxp(ie)=nxp(ie)+1
       if(xmmc.lt.xu)return
       i=1+int(log(xmmc/xu)/log(xo/xu)*nbix)
       if(i.gt.nbix)goto1
       if(i.lt.1)goto1
       wxmmc(i,ie)=wxmmc(i,ie)+1
       nxm(ie)=nxm(ie)+1
1      continue
       if(xmmc.ge.xu)then
         ymc=0.5*log(xpmc/xmmc)
       else
         return
       endif
       if(ymc.lt.yu)return
       i=int((ymc-yu)/dy)+1
       if(i.gt.nbiy)return
       if(i.lt.1)return
       wymc(i,ie)=wymc(i,ie)+1
       ny(ie)=ny(ie)+1

       xmc=xpmc*xmmc
       if(xmc.lt.xu)return
       i=1+int(log(xmc/xu)/log(xo/xu)*nbix)
       if(i.gt.nbix)return
       if(i.lt.1)return
       wxmc(i,ie)=wxmc(i,ie)+1
       nx(ie)=nx(ie)+1

      elseif(iii.eq.2)then

        nstat=nfull
        if(nstat.lt.0)nstat=max(1,nevent)

        do ii=1,nie

       if(ii.eq.1)write(ifhi,'(a)')'!-----  projectile droplet  ----'
       if(ii.eq.2)write(ifhi,'(a)')'!-----    target droplet    ----'
       if(ii.eq.3)write(ifhi,'(a)')'!-----  projectile string end  ----'
       if(ii.eq.4)write(ifhi,'(a)')'!-----    target string end    ----'
        write(ifhi,'(a)')       '!--------------------------------'
        write(ifhi,'(a)')       '!   droplet/string x+ distr       '
        write(ifhi,'(a)')       '!--------------------------------'
        write(ifhi,'(a)')       'openhisto'
        write(ifhi,'(a)')       'htyp lru'
        write(ifhi,'(a)')       'xmod log ymod log'
        write(ifhi,'(a,2e11.3)')'xrange',xu,xo
        if(ii.eq.1.or.ii.eq.2)
     *  write(ifhi,'(a)')    'text 0 0 "xaxis droplet x+"'
        if(ii.eq.3.or.ii.eq.4)
     *  write(ifhi,'(a)')    'text 0 0 "xaxis string end x+"'
        write(ifhi,'(a)')    'text 0 0 "yaxis P(x)"'
        write(ifhi,'(a,e22.14)')'histoweight ',dble(nstat)
        write(ifhi,'(a)')       ' array 2'
        do i=1,nbix
         dx=xu*(xo/xu)**(1.*i/nbix)*(1.-(xo/xu)**(-1./nbix))
         if(nxp(ii).gt.0)
     *   write(ifhi,'(2e11.3)')x(i),wxpmc(i,ii)/dx/nxp(ii)
        enddo
        write(ifhi,'(a)')    '  endarray'
        write(ifhi,'(a)')    'closehisto plot 0-'
        write(ifhi,'(a)')       '!--------------------------------'
        write(ifhi,'(a)')       '!   droplet/string x- distr       '
        write(ifhi,'(a)')       '!--------------------------------'
        write(ifhi,'(a)')       'openhisto'
        write(ifhi,'(a)')       'htyp lba'
        write(ifhi,'(a)')       'xmod log ymod log'
        write(ifhi,'(a,2e11.3)')'xrange',xu,xo
        if(ii.eq.1.or.ii.eq.2)
     *  write(ifhi,'(a)')    'text 0 0 "xaxis droplet x-"'
        if(ii.eq.3.or.ii.eq.4)
     *  write(ifhi,'(a)')    'text 0 0 "xaxis string end x-"'
        write(ifhi,'(a)')    'text 0 0 "yaxis P(x)"'
        write(ifhi,'(a,e22.14)')'histoweight ',dble(nstat)
        write(ifhi,'(a)')       ' array 2'
        do i=1,nbix
         dx=xu*(xo/xu)**(1.*i/nbix)*(1.-(xo/xu)**(-1./nbix))
         if(nxm(ii).gt.0)
     *   write(ifhi,'(2e11.3)')x(i),wxmmc(i,ii)/dx/nxm(ii)
        enddo
        write(ifhi,'(a)')    '  endarray'
        write(ifhi,'(a)')    'closehisto plot 0'
        write(ifhi,'(a)')       '!--------------------------------'
        write(ifhi,'(a)')       '!   droplet/string y distr       '
        write(ifhi,'(a)')       '!--------------------------------'
        write(ifhi,'(a)')       'openhisto'
        write(ifhi,'(a)')       'htyp lin'
        write(ifhi,'(a)')       'xmod lin ymod lin'
        write(ifhi,'(a,2e11.3)')'xrange',yu,yo
        if(ii.eq.1.or.ii.eq.2)
     *  write(ifhi,'(a)')    'text 0 0 "xaxis droplet y"'
        if(ii.eq.3.or.ii.eq.4)
     *  write(ifhi,'(a)')    'text 0 0 "xaxis string end y"'
        write(ifhi,'(a)')    'text 0 0 "yaxis P(y)"'
        write(ifhi,'(a,e22.14)')'histoweight ',dble(nstat)
        write(ifhi,'(a)')       ' array 2'
        do i=1,nbiy
         if(ny(ii).gt.0)
     *   write(ifhi,'(2e11.3)')y(i),wymc(i,ii)/dy/ny(ii)
        enddo
        write(ifhi,'(a)')    '  endarray'
        write(ifhi,'(a)')    'closehisto plot 0'

      enddo

        write(ifhi,'(a)')       '!--------------------------------'
        write(ifhi,'(a)')       '!   droplet/string mass distr       '
        write(ifhi,'(a)')       '!--------------------------------'
      do ii=1,nie


        if(ii.eq.2.or.ii.eq.4)write(ifhi,'(a)')    'closehisto plot 0-'
        if(ii.eq.3)write(ifhi,'(a)')    'closehisto plot 0'
        write(ifhi,'(a)')       'openhisto'
        if(ii.eq.1.or.ii.eq.3)write(ifhi,'(a)')       'htyp lru'
        if(ii.eq.2.or.ii.eq.4)write(ifhi,'(a)')       'htyp lba'
        write(ifhi,'(a)')       'xmod log ymod log'
        write(ifhi,'(a,2e11.3)')'xrange',sqrt(xu*s),sqrt(s*xo)
        if(ii.eq.1.or.ii.eq.2)
     *  write(ifhi,'(a)')    'text 0 0 "xaxis droplet mass (GeV)"'
        if(ii.eq.4.or.ii.eq.3)
     *  write(ifhi,'(a)')    'text 0 0 "xaxis string end mass (GeV)"'
        write(ifhi,'(a)')    'text 0 0 "yaxis P(x)"'
        write(ifhi,'(a,e22.14)')'histoweight ',dble(nstat)
        write(ifhi,'(a)')       ' array 2'
        do i=1,nbix
         dx=xu*(xo/xu)**(1.*i/nbix)*(1.-(xo/xu)**(-1./nbix))
         if(nx(ii).gt.0)
     *   write(ifhi,'(2e11.3)')sqrt(x(i)*s),wxmc(i,ii)/dx/nx(ii)
        enddo
        write(ifhi,'(a)')    '  endarray'
      enddo
       write(ifhi,'(a)')    'closehisto plot 0'

      endif

      return
      end

cc--------------------------------------------------------------------------
c      subroutine xtype(k,n,i1,i2,text)
cc--------------------------------------------------------------------------
c
c#include "aaa.h"
c#include "ems.h"
c      parameter(itext=40)
c      character  text*40
c
c      imax=itext+1
c      do i=itext,1,-1
c      if(text(i:i).eq.'&')imax=i
c      enddo
c
c      ip=iproj(k)
c      it=itarg(k)
c
c      if(i1.eq.1)then
c         write(ifch,*)
c         write(ifch,*)('-',ll=1,27)
c         write(ifch,*)'  '//text(1:imax-1)
c         write(ifch,*)('-',ll=1,27)
c      endif
c
c      if(i2.eq.1)then
c         write(ifch,*)
c         write(ifch,*)'k:',k,'   n:',n,'   ip:',ip,'   it:',it
c         write(ifch,*)'bk:',bk(k)
c         if(n.ne.0)write(ifch,*)'idpr:',idpr(n,k)
c         write(ifch,*)'iep:',iep(ip),'   iet:',iet(it)
c         write(ifch,*)'idp:',idp(ip),'   idt:',idt(it)
c      endif
c
c      end
c
c------------------------------------------------------------------------
      subroutine XPrint(text)
c------------------------------------------------------------------------
#include "aaa.h"
#include "ems.h"
      double precision xpptot,xmptot,xpttot,xmttot
      character  text*(*)
      imax=index(text,'&')
      write(ifch,'(1x,a)')text(1:imax-1)

      write(ifch,'(a)')
     *' k:     itpr:   npr0:  npr: nprmx:   Pomeron (id,iv) lattice:'
      do k=1,koll
       npr13=npr(1,k)+npr(3,k)
       write(ifch,'(1x,i6,1x,i4,4x,i4,2x,i4,3x,i4,a3,$)')
     *              k,itpr(k),npr(0,k),npr13,nprmx(k),'   '
       do n=1,nprmx(k)
        write(ifch,'(1x,a,2i2,a,$)')'(',idpr(n,k),ivpr(n,k),')'
       enddo
       write(ifch,*)' '
      enddo

      xpptot=0d0
      xmptot=0d0
      xpttot=0d0
      xmttot=0d0
      write(ifch,'(a)')' Pomeron xy lattice:'
      do k=1,koll
       do n=1,nprmx(k)
         if(idpr(n,k).gt.0.and.ivpr(n,k).eq.1)then
           xpptot=xpptot+xppr(n,k)
           xmttot=xmttot+xmpr(n,k)
         endif
         write(ifch,'(i6,1x,i2,1x,d10.3,1x,d10.3,1x,i1,3x,$)')
     *                  k,n,xpr(n,k),ypr(n,k),idfpr(n,k)
       enddo
       write(ifch,*)' '
      enddo

      write(ifch,'(a)')' projectile remnants x+,x-,px,py,x,iep:'
      do ip=1,maproj
       xpptot=xpptot+xpp(ip)+xmpmn(ip)
       xmptot=xmptot+xmp(ip)-xmpmn(ip)
       write(ifch,'(i3,2x,5d12.3,i3)')ip,xpp(ip),xmp(ip),xxp(ip),xyp(ip)
     *                             ,xpos(ip),iep(ip)
      enddo

      write(ifch,'(a)')' target remnants x-,x+,px,py,x,iet:'
      do it=1,matarg
       xpttot=xpttot+xpt(it)-xptmn(it)
       xmttot=xmttot+xmt(it)+xptmn(it)
       write(ifch,'(i3,2x,5d12.3,i3)')it,xmt(it),xpt(it),xxt(it),xyt(it)
     *                             ,xtos(it),iet(it)
      enddo

      write(ifch,*)' remnant balance x+,x-:'
     &,(xpptot+xpttot)/dble(maproj)
     &,(xmptot+xmttot)/dble(matarg)
      end

c-----------------------------------------------------------------------
      subroutine xbDens(jjj)
c-----------------------------------------------------------------------
c plots b distribution for all pairs
c----------------------------------------------------------------
#include "aaa.h"
#include "ems.h"
      common/geom/rmproj,rmtarg,bmax,bkmx

      b=0.
      if(jjj.eq.1)then
c prepare plot for xbDens
      if(ixbDens.eq.1)then
        iii=1     !proj
        Nnucla=0
        do ip=1,maproj
          if(lproj(ip).ne.0)then
            Nnucla=Nnucla+1
            do l=1,lproj(ip)
              k=kproj(ip,l)
              b=bk(k)
              i=1+int(b/bkmx*float(mxnucl))
              if(i.le.mxnucl)bnucl(i,iii)=bnucl(i,iii)+1.
            enddo
          endif
          if(lproj3(ip).ne.0)then
            do l=1,lproj3(ip)
              k=kproj3(ip,l)
              b=bk(k)
              i=1+int(b/bkmx*float(mxnucl))
              if(i.le.mxnucl)bnucl(i,iii+2)=bnucl(i,iii+2)+1.
            enddo
          endif
        enddo
        xbtot(iii)=xbtot(iii)+float(Nnucla)
        iii=2     !targ
        Nnucla=0
        do it=1,matarg
          if(ltarg(it).ne.0)then
            Nnucla=Nnucla+1
            do l=1,ltarg(it)
              k=ktarg(it,l)
              b=bk(k)
              i=1+int(b/bkmx*float(mxnucl))
              if(i.le.mxnucl)bnucl(i,iii)=bnucl(i,iii)+1.
            enddo
          endif
          if(ltarg3(it).ne.0)then
            do l=1,ltarg3(it)
              k=ktarg3(it,l)
              b=bk(k)
              i=1+int(b/bkmx*float(mxnucl))
              if(i.le.mxnucl)bnucl(i,iii+2)=bnucl(i,iii+2)+1.
            enddo
          endif
        enddo
        xbtot(iii)=xbtot(iii)+float(Nnucla)
      endif

      else

        nstat=nfull
        if(nstat.lt.0)nstat=max(1,nevent)

      if(xbtot(1).gt.0.)then
        xbtot(3)=xbtot(1)
        xbtot(4)=xbtot(2)
        write(ifhi,'(a)')       'openhisto'
        write(ifhi,'(a)')       'htyp lin name bdens'
        write(ifhi,'(a)')       '- txt "xaxis b (fm)" '
        write(ifhi,'(a)')       '+ txt "yaxis P(b) proj " '
        write(ifhi,'(a)')       '+ txt "yaxis P(b) targ " '
        write(ifhi,'(a)')       '+ txt "yaxis P(b) scr proj " '
        write(ifhi,'(a)')       '+ txt "yaxis P(b) scr targ " '
        write(ifhi,'(a,e22.14)')'histoweight ',dble(nstat)
        write(ifhi,'(a)')       ' array 5'
        db=bkmx/float(mxnucl)
        do j=1,mxnucl
          b=(j-0.5)*db
          d=pi*((b+db)**2-b**2)
          write(ifhi,'(2e12.4)') b,(bnucl(j,iii)/xbtot(iii)/d,iii=1,4)
        enddo
        write(ifhi,'(a)')       '  endarray'
        write(ifhi,'(a)')       'closehisto'
        write(ifhi,'(a)')       'plot bdens+1- plot bdens+2-'
        write(ifhi,'(a)')       'plot bdens+3- plot bdens+4 '
      endif

      endif

      end

