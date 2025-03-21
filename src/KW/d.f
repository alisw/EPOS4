C
C  This file is part of EPOS4
C  Copyright (C) 2022 research institutions and authors (See CREDITS file)
C  This file is distributed under the terms of the GNU General Public License version 3 or later
C  (See COPYING file for the text of the licence)
C

c-----------------------------------------------------------------------
      subroutine setinp(ifname,nifname,ifop)
c-----------------------------------------------------------------------
      character*(*) ifname
      ifop=19
      open(UNIT=ifop,FILE=ifname(1:nifname),STATUS='UNKNOWN')
      write(*, '(a,a)')'Connecting input to file: ', ifname(1:nifname)
      end

c-----------------------------------------------------------------------
      subroutine setRivet(line,i,j)
c-----------------------------------------------------------------------
#include "aaa.h"
      character*1000 line
      logical exist
      inquire(file=fndt(1:nfndt)//'R',exist=exist)
      if(exist)then
        open(17,file=fndt(1:nfndt)//'R',status="old", 
     .    position="append")
      else
        open(17,file=fndt(1:nfndt)//'R',status="new")
      endif
      call utword(line,i,j,0)
      write(17,'(a,$)')line(i:j)   
      call utword(line,i,j,0)
      if(line(j+2:j+7).eq.'calib:')then 
        write(17,'(2a,$)')'  ',line(i:j)   
        call utword(line,i,j,0)
        write(17,'(2a)')'  ',line(i:j)   
      else 
        write(17,'(2a)')'  ',line(i:j)   
      endif
      close (17)
      end

c-----------------------------------------------------------------------
      subroutine setIfBigSystem(line,i,j,igo1)
c-----------------------------------------------------------------------
#include "aaa.h"
      character*1000 line
      igo1=0
      if(.not.(maproj*matarg.gt.5000))then 
      iechox=iecho
      iecho=0
      call utword(line,i,j,ne)
      do while(line(i:j).ne.'#fiBigSystem')
      call utword(line,i,j,ne)
      enddo
      iecho=iechox
      igo1=1
      endif
      end

c-----------------------------------------------------------------------
      subroutine setIfSmallSystem(line,i,j,igo1)
c-----------------------------------------------------------------------
#include "aaa.h"
      character*1000 line
      igo1=0
      if(.not.(maproj*matarg.le.5000))then 
      iechox=iecho
      iecho=0
      call utword(line,i,j,ne)
      do while(line(i:j).ne.'#fiSmallSystem')
      call utword(line,i,j,ne)
      enddo
      iecho=iechox
      igo1=1
      endif
      end

c-----------------------------------------------------------------------
      subroutine setIf1(line,i,j) !#if1 definition
c-----------------------------------------------------------------------
#include "aaa.h"
      character*1000 line
      character cext1*10
      common/ccext1/cext1
      character*10 celse(20)
      integer nelse,inotelse
      common/ccelse/nelse,inotelse,celse
      lcext1=index(cext1,' ')-1
      inot=0
      call utword(line,i,j,0)
      if(line(i:j).eq.'#not')then
        inot=1
        call utword(line,i,j,0)
      endif
      nif1=1
      if(line(i:i).eq.'#')then
        read(line(i+1:j),*)nif1
        call utword(line,i,j,0)
      endif
      nelse=nif1
      iskip=2 !skip
      celse(1)='          '
      celse(1)(1:j-i+1)=line(i:j)
      if(line(i:j).eq.cext1(1:lcext1))iskip=1 !noskip
      do nuw=2,nif1
        call utword(line,i,j,0)
        celse(nuw)='          '
        celse(nuw)(1:j-i+1)=line(i:j)
        if(line(i:j).eq.cext1(1:lcext1))iskip=1 !noskip
      enddo
      if(inot.eq.1)iskip=3-iskip !invert
      if(iskip.eq.2)then !skip
        iechosave=iecho
        iecho=0
        do while(line(i:j).ne.'#else1'
     .   .and.line(i:j).ne.'#fi'.and.line(i:j).ne.'#fi1')
          call utword(line,i,j,-1)
        enddo
        iecho=iechosave
        if(line(i:j).eq.'#else1')j=max(0,j-7)
      else
        kmax=2
        do k=3,1000
        if(line(k:k).ne.' ')kmax=k
        enddo
        if(iecho.eq.1)write(ifmt,'(a)')line(1:kmax)
      endif
      end

c-----------------------------------------------------------------------
      subroutine setElse1(line,i,j) !#else1 definition
c-----------------------------------------------------------------------
#include "aaa.h"
      character*1000 line
      character cext1*10
      common/ccext1/cext1
      character*10 celse(20)
      integer nelse,inotelse
      common/ccelse/nelse,inotelse,celse
      lcext1=index(cext1,' ')-1
      inot=inotelse
      nif1=nelse
      iskip=2 !skip
      lelse=index(celse(1),' ')-1
      if(celse(1)(1:lelse).eq.cext1(1:lcext1))iskip=1 !noskip
      do nuw=2,nif1
        lelse=index(celse(nuw),' ')-1
        if(celse(nuw)(1:lelse).eq.cext1(1:lcext1))iskip=1 !noskip
      enddo
      if(inot.eq.0)iskip=3-iskip !invert
      if(iskip.eq.2)then !skip
        iechosave=iecho
        iecho=0
        do while(line(i:j).ne.'#fi'.and.line(i:j).ne.'#fi1')
          call utword(line,i,j,-1)
        enddo
        iecho=iechosave
      else
        kmax=2
        do k=3,1000
        if(line(k:k).ne.' ')kmax=k
        enddo
        if(iecho.eq.1)write(ifmt,'(a)')line(1:kmax)
      endif
      end

c-----------------------------------------------------------------------
      subroutine setCentralityClass(line,i,j,nopen)
c-----------------------------------------------------------------------
      common/files/ifop,ifmt,ifch,ifcx,ifhi,ifdt,ifcp,ifdr,ifio
      common/hydr4/zclass(5,100),izmode,maxsurf,iofrout,jzmode(8)
      common/ccc20/icc20
      character line*1000
      call utword(line,i,j,0)
      if(line(i:j).eq.'cc20')then
        icc20=1
        return
      endif
      read(line(i:j),*)val
      ival=nint(val)
      jval=ival
      call utword(line,i,j,0)
      if(line(i:j).eq.'-')then
        call utword(line,i,j,0)
        read(line(i:j),*)val
        jval=nint(val)
        call utword(line,i,j,0)
      endif
      read(line(i:j),*)val1
      call utword(line,i,j,0)
      read(line(i:j),*)val2
      call utword(line,i,j,0)
      read(line(i:j),*)val3
      if(izmode.eq.2)then
        call utword(line,i,j,0)
        read(line(i:j),*)val4
        call utword(line,i,j,0)
        read(line(i:j),*)val5
      endif
      do i=ival,jval
      if(nopen.ne.-1)then   !only first read
        if(zclass(3,i).gt.0.001)then
          write(ifmt,'(a,i3)') 'ERROR 24022013 zclass3 in use, i=',i
          stop'24022013'
        endif
      endif
      zclass(1,i)=val1
      zclass(2,i)=val2
      zclass(3,i)=val3
      if(izmode.eq.2)then
        zclass(4,i)=val4
        zclass(5,i)=val5
      endif
      enddo
      end

c-----------------------------------------------------------------------
      subroutine setCentrality(line,i,j,iopcnt,ifmt
     .               ,icentrality,jcentrality,ffac,imax,ival,ishxxx)
c-----------------------------------------------------------------------
      common/hydr4/zclass(5,100),izmode,maxsurf,iofrout,jzmode(8)
      common/crootcproot/irootcproot,iboein
      character line*1000
      
      if(line(i:j).eq.'within')then
        if(iopcnt.eq.0)then
          write(ifmt,'(//a/a/a)')
     .     'ERROR in aread:'
     .    ,'centrality within  and  iopcnt = 0  incompatible'
     .    ,'centrality within  requires batch mode!'
          stop'in aread \n\n '
        endif
        call utword(line,i,j,0)
        if(line(i:j).ne.'{')stop'\n\n ERROR 19112011\n\n'
        k=0
        i0=i
        j0=j
        do while(line(i:j).ne.'}')
          call utword(line,i,j,0)
          k=k+1
        enddo
        ival=k-1
        i=i0
        j=j0
        ij=1+mod((iopcnt-1),ival)
        do k=1,ival
          call utword(line,i,j,0)
          read(line(i:j),*)val1
          if(k.eq.ij)then     
            icentrality=nint(val1)
          endif
        enddo
        call utword(line,i,j,0)
        if(line(i:j).ne.'}')stop'\n\n ERROR 19112011b\n\n'
        jcentrality=icentrality
        ffac=zclass(3,icentrality)
        imax=1
        if(ishxxx.eq.1)
     .    write(ifmt,'(a,i3)')'centrality set to',icentrality
      else
        read(line(i:j),*)val
        ival=nint(val)
        if(ival.lt.0)then
          jcentrality=ival
          if(iopcnt.eq.0)then
            write(ifmt,'(//a/a/a)')
     .       'ERROR in aread:'
     .      ,'centrality < 0  and  iopcnt = 0  incompatible'
     .      ,'centrality < 0  requires batch mode!'
            stop'in aread \n\n '
          endif
          imax=0
          do i=1,100
           !print*,'1111111111111111',i,zclass(3,i)
           if(zclass(3,i).gt.1e-5)imax=i
          enddo
          icentrality=1+mod((iopcnt-1)/(-jcentrality),imax)
          ffac=zclass(3,icentrality)
          if(ishxxx.eq.1)
     .      write(ifmt,'(a,i3)')'centrality set to',icentrality
        else
          icentrality=ival
          jcentrality=ival
          ffac=1.00000
          imax=1
        endif
      endif

      end

c-----------------------------------------------------------------------
      subroutine setFillTree(line,i,j,irootcproot,ifillTree,izmode)
c-----------------------------------------------------------------------
      character line*1000
      common/cigrpac/igrpac
      igrpac=1
      dummy=i
      if(line(j-2:j-2).ne.'C')stop'fillTree error\n\n'
      if(irootcproot.ne.1)then !----yes-no-----no-no---(yes-yes)
        read(line(j-1:j-1),*)ifillTree !1=bim, 2=kol
      else  !-----case root-cproot------no-yes----
        read(line(j-1:j-1),*)izmode
      endif
      end


c-----------------------------------------------------------------------
      subroutine setGetTree(line,i,j,irootcproot,igrTree,muTree)
c-----------------------------------------------------------------------
      character cext3*10
      common/ccext3/iext3,cext3 
      common/ciext4/iext4
      character line*1000
      igrpac=1
      igrTree=0
      if(irootcproot.eq.1)then
      ix2=i+6
      nu=0
      igruu=0
      igrxx=0
      igryy=0
      do
        if(ix2+2.gt.j)exit
        ix1=ix2+2
        ix2=ix1
        do while(line(ix2+1:ix2+1).ne.','.and.line(ix2+1:ix2+1).ne.')')
        ix2=ix2+1
        enddo
        read(line(ix1:ix2),*)val
        nu=nu+1
        if(nu.eq.1)igruu=nint(val)
        if(nu.eq.2)muTree=val
        if(nu.eq.3)igrxx=nint(val)
        if(nu.eq.4)igryy=nint(val)
        if(nu.eq.5)stop'\n\n error 21072010b \n\n'
      enddo
      igrTree=igruu
      if(cext3(1:iext3).eq.'0')then ! ordering option
        igrTree=igrxx
      endif
      if(iext4.eq.1)then
        igrTree=igryy
      endif
      endif
      end

c-----------------------------------------------------------------------
      subroutine setOrderTree(line,i,j,irootcproot
     .             ,maxpom,igrpac,nfifac)
c-----------------------------------------------------------------------
      character line*1000
      if(irootcproot.eq.1)then
      ix2=i+8
      nu=0
      do
        if(ix2+2.gt.j)exit
        ix1=ix2+2
        ix2=ix1
        do while(line(ix2+1:ix2+1).ne.','.and.line(ix2+1:ix2+1).ne.')')
        ix2=ix2+1
        enddo
        read(line(ix1:ix2),*)val
        nu=nu+1
        if(nu.eq.1)maxpom=nint(val)
        if(nu.eq.2)igrpac=nint(val)
        if(nu.eq.3)stop'\n\n error 08062013 \n\n'
      enddo
      if(nu.ge.2)then
        nfifac=maxpom*igrpac
      elseif(nu.eq.1)then 
        nfifac=maxpom
      endif
      endif
      end

c-----------------------------------------------------------------------
      subroutine setInput(line,i,j,nopen,iprmpt,iopcnt)
c-----------------------------------------------------------------------
      character line*1000, cfmt*4
      common/ciexhd/iexhd
      call utword(line,i,j,0)
      if(nopen.ge.0)then
       nopen=nopen+1
       if(nopen.gt.9)stop'too many nested input commands'
c       print*,line(i:j)
       open(unit=20+nopen,file=line(i:j),status='old')
       if(iprmpt.eq.1)iprmpt=-1
       if(line(j-5:j).eq.'.optns')then
          n=0
        if(line(j-7:j-7).eq.'-')then
          n=1
        elseif(line(j-8:j-8).eq.'-')then
          n=2
        elseif(line(j-9:j-9).eq.'-')then
          n=3
        elseif(line(j-10:j-10).eq.'-')then
          n=4
        elseif(line(j-11:j-11).eq.'-')then
          n=5
        elseif(line(j-12:j-12).eq.'-')then
          n=6
        endif
        if(n.gt.0.and.line(j-8-n:j-6-n).eq.'/z-')n=0
        iopcnt=0
        if(n.gt.0)then
          cfmt='(i )'
          write(cfmt(3:3),'(i1)')n
          jj=j-5-n
c         in case of string cast to integer error, iopcnt equals 0
          read(line(jj:j-6),cfmt,iostat=ios)iopcnt
        endif
        do k=i,j-4
          if(line(k:k+4).eq.'ppihd')iexhd=1
        enddo
       endif
      endif
      end

c-----------------------------------------------------------------------
      subroutine setXinput(line,i,j,nopen,iprmpt,fnnx,nfnnx)
c-----------------------------------------------------------------------
      character line*1000, fnnx*500
      common/ciexhd/iexhd
      call utword(line,i,j,0)
      if(nopen.ge.0)then
       nopen=nopen+1
       if(nopen.gt.9)stop'too many nested input commands'
       if(line(i:i+3).eq.'KWt/')then
         open(unit=20+nopen,file=fnnx(1:nfnnx)//'src/'//line(i:j)
     .   ,status='old')
       else
         open(unit=20+nopen,file=fnnx(1:nfnnx)//line(i:j)
     .   ,status='old')
       endif 
       if(iprmpt.eq.1)iprmpt=-1
       do k=i,j-5
         if(line(k:k+8).eq.'iKWhd/ihd')iexhd=1
       enddo
      endif
      end
   
c-----------------------------------------------------------------------
      subroutine setSystem(line,i,j,isyst)
c-----------------------------------------------------------------------
      character line*1000
      call utword(line,i,j,0)
      if(line(i:j).eq.'i')isyst=0
      if(line(i:j).eq.'q')isyst=1
      if(line(i:j).eq.'c')isyst=2
      if(line(i:j).eq.'-q')isyst=-1
      if(line(i:j).eq.'-c')isyst=-2
      end

c-----------------------------------------------------------------------
      subroutine setRootcproot(line,i,j,irootcproot,iboein)
c-----------------------------------------------------------------------
      character line*1000
      call utword(line,i,j,0)
      if(line(i:j).eq.'nono') irootcproot= 0
      if(line(i:j).eq.'noyes')irootcproot= 1
      if(line(i:j).eq.'yesno')irootcproot=10
      if(line(i:j).eq.'yesyes')stop'\n\n invalid rootcproot \n\n'
      if(line(i:j).eq.'noyesBE')irootcproot= 1
      if(line(i:j).eq.'noyesBE')iboein=1
      end

c-----------------------------------------------------------------------
      subroutine setDefine(line,i,j,k)
c-----------------------------------------------------------------------
#include "aaa.h"
      character line*1000
      call utword(line,i,j,ne)
      if(line(i:j).eq.'bim3')stop'****** bim3->bim03 ****** '
      if(line(i:j).eq.'bim5')stop'****** bim5->bim05 ****** '
      if(line(i:j).eq.'bim6')stop'****** bim6->bim06 ****** '
      if(ndefine+1.gt.mxdefine)stop'too many `define` statements.      '
      l1=j-i+1
      if(l1.gt.99)stop'`define` argument 1 too long.            '
      w1define(ndefine+1)(1:l1)=line(i:j)
      w1define(ndefine+1)(l1+1:l1+1)=' '
      call utword(line,i,j,ne)
      l2=j-i+1
      if(l2.gt.99)stop'`define` argument 2 too long.            '
      w2define(ndefine+1)(1:l2)=line(i:j)
      w2define(ndefine+1)(l2+1:l2+1)=' '
      ndefine=ndefine+1
      l1define(ndefine)=l1
      l2define(ndefine)=l2
      if(k.ge.2)then
        call utword(line,i,j,ne)
        call utword(line,i,j,ne)
      endif
      if(k.ge.3)then
        call utword(line,i,j,ne)
      endif
      end

