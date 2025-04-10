C
C  This file is part of EPOS4
C  Copyright (C) 2022 research institutions and authors (See CREDITS file)
C  This file is distributed under the terms of the GNU General Public License version 3 or later
C  (See COPYING file for the text of the licence)
C
      program eposc     
      common/hydr4/zclass(5,100),izmode,maxsurf,iofrout,jzmode(8)
      character*500 fidata,fiepos
      character line*1000,linex*1000
      common /cwo/nopen,ifop
      common/ccc20/icc20
      character*10 cext1
      common/ccext1/cext1
      character cext3*10
      common/ccext3/iext3,cext3 /ciext4/iext4
      common/cnfifac/nfifac
      character*1000 fnjob
      common/jobfname/  fnjob,nfnjob
      igrTree=0
      jcentrality=0
      icentrality=0     
      izmode=1
      nskip=0
      ishxxx=0
      do i=1,100
      zclass(1,i)=0.
      zclass(2,i)=0.
      zclass(3,i)=0.
      zclass(4,i)=-1.
      zclass(5,i)=-1.
      enddo
      nfifac=1
      nopen=0
      iprmpt=-2
      iopcnt=0
      isyst=0
      ndefine=0
      ifop=5     
      ifdt=51   
      cext1='          '

c     to read input from different file name given as argument
      j=-1
      call utword(line,i,j,0)
      fnjob=line(i:j)
      nfnjob=j-i+1

      call setinp(fnjob,nfnjob,ifop)

      j=-1

    1 call utword(line,i,j,1)
    
          if(line(i:j).eq.'#define')then

      call setDefine(line,i,j,1)

          elseif(line(i:j).eq.'#define2')then

      call setDefine(line,i,j,2)

          elseif(line(i:j).eq.'#define3')then

      call setDefine(line,i,j,3)

           elseif(line(i:j).eq.'CentralityClass')then
          
      call setCentralityClass(line,i,j,nopen)        

           elseif(line(i:j).eq.'ext1')then
           
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
      !print*,'isyst =',isyst

           elseif(line(i:j).eq.'input')then

      call setInput(line,i,j,nopen,iprmpt,iopcnt)
      !print*,'iopcnt =',iopcnt
          
           elseif(line(i:j).eq.'xinput')then

      call setXinput(line,i,j,nopen,iprmpt,fiepos,nfiepos)

           elseif(line(i:j).eq.'rootcproot')then

      call setRootcproot(line,i,j,irootcproot,iboein)
      !print*,'irootcproot =',irootcproot

          elseif(line(i:i+8).eq.'fillTree(')then

      write(ifmt,'(//70a/a/a/a/70a/)')('-',k=1,70)
     .,'Important info: New definitions of ior and jor, referring  '
     .,'                 to indices 0,1,2,... rather than 1,2,3,...' 
     .,'Use now fillTree4 rather than fillTree in optns'
     .,('-',k=1,70)
      stop'fillTree not known'

           elseif(line(i:i+8).eq.'fillTree4')then
          
      call setFillTree(line,i,j,irootcproot,ifillTree,izmode)
      !print*,'ifillTree =',ifillTree
      
           elseif(line(i:i+6).eq.'getTree')then
      
      call setGetTree(line,i,j,irootcproot,igrTree,muTree)
      !print*,'igrTree =',igrTree
      
          elseif(line(i:i+8).eq.'orderTree')then

      if(irootcproot.eq.1)then
      if(nopen.ne.-1)then               !only first read
      call setOrderTree(line,i,j,irootcproot,maxpom,igrpac,nfifac)
      endif
      endif

           elseif(line(1:10).eq.'beginwrite')then
      
      call utword(line,i,j,0)
      do while(line(i:j).ne.'endwrite')
        call utword(line,i,j,0)      
      enddo
      
           elseif(line(i:j).eq.'set')then

      call utword(line,i,j,0)
      linex=line
      ix=i
      jx=j
      call utword(line,i,j,0)
      if(linex(ix:jx).ne.'centrality')then
      read(line(i:j),*)val
      endif
      if(linex(ix:jx).eq.'izmode')izmode=nint(val)
      if(linex(ix:jx).eq.'nskip')nskip=nint(val)
      if(linex(ix:jx).eq.'centrality')then
      call setCentrality(line,i,j,iopcnt,ifmt
     .        ,icentrality,jcentrality,ffac,imax,ival,ishxxx)
      ncentrality=imax
      !print*,'imax=',imax
      endif

           elseif(line(i:j).eq.'fname')then

      call utword(line,i,j,0)
      linex=line
      ix=i
      jx=j
      call utword(line,i,j,0)
      if(linex(ix:jx).eq.'pathep')fiepos(1:j-i+1)=line(i:j)
      if(linex(ix:jx).eq.'pathep')nfiepos=j-i+1
      if(linex(ix:jx).eq.'data')then
        fidata(1:j-i+1)=line(i:j)
        nfidata=j-i+1
        !print*,'data file = ', fidata(1:nfidata)
        open(unit=ifdt,file=fidata(1:nfidata),status='unknown')
      endif  

           elseif(line(i:j).eq.'rivet')then

      call utword(line,i,j,0)
      call utword(line,i,j,0)
      if(line(j+2:j+7).eq.'calib:')call utword(line,i,j,0)

           elseif(line(i:j).eq.'runprogram')then

      if(nskip.lt.0)then
        jjtb=mod(-nskip,10)
      endif
      if(fidata(nfidata:nfidata).eq.'X')then !only called from epos+c
        ixx=1
        if(jcentrality.lt.0)then
        ixx=nint(zclass(3,icentrality))
        if(izmode.eq.2)ixx=1
        if(nskip.eq.1)ixx=1
        if(irootcproot.eq.1)ixx=1
        endif
        if(nskip.lt.0)ixx=jjtb
        write(ifdt,*)ixx
        open(17,file=fidata(1:nfidata-1)//'Y',status='unknown')
        iyy=ncentrality
        if(nskip.lt.0)iyy=1
        write(17,*)iyy
        close (17)
        open(17,file=fidata(1:nfidata-1)//'Z',status='unknown')
        write(17,*)igrTree
        close (17)
        !print*,'XYZ',ixx,ncentrality,igrTree
        open(17,file=fidata(1:nfidata-1)//'W',status='unknown')
        write(17,*)nfifac
        close (17)
        stop
      endif
       
           endif
        
      goto 1
        
      end


c-------------------------------------------------------------------
      subroutine utword(line,i,j,iqu)

ccc
ccc   copied from uti.f, but much simplified
ccc
ccc         /cwo/ instead of aaa.h
ccc

      common /cwo/nopen,ifop
          
      parameter(mempty=1)
      character*1 empty(mempty),mk
      character line*1000
      character*2 mrk
      data empty/' '/

      if(j.ge.0)then
      i=j
      goto 1
      endif

    5 continue
      if(iqu.eq.1.and.iprmpt.gt.0)write(ifmt,'(a)')'?'
      if(nopen.eq.0)ifopx=ifop
      if(nopen.gt.0)ifopx=20+nopen
      read(ifopx,'(a1000)',end=9999)line
      i=0

    1 i=i+1
      if(i.gt.1000)goto 5
      if(line(i:i).eq.'!')goto 5
      do ne=1,mempty
      if(line(i:i).eq.empty(ne))goto 1
      enddo

      nbla=1
      mrk='  '
      mk=' '
      if(line(i:i).eq.'~')mk='~'
      if(line(i:i+1).eq.'"{')mrk='}"'
      if(line(i:i+1).eq.'""')mrk='""'
      if(mrk.ne.'  ')goto 10
      if(line(i:i).eq.'"')mk='"'
      if(mk.ne.' ')goto 8
      j=i-1
    6 j=j+1
      if(j.gt.1000)goto 7
      if(line(j:j).eq.'!')goto 7
      do ne=1,mempty
      if(line(j:j).eq.empty(ne))goto 7
      enddo
      goto 6

    8 continue
      if(i.ge.1000-1)stop'utword: make line shorter!!!         '
      i=i+1
      j=i
      if(line(j:j).eq.mk)stop'utword: empty string!!!           '
    9 j=j+1
      if(j.gt.1000)then                 !reach the end of the line
        j=j-nbla+2
        goto 7
      endif
      if(line(j:j).eq.' ')then
        nbla=nbla+1
      else
        nbla=2
      endif
      if(line(j:j).eq.mk)then
      line(i-1:i-1)=' '
      line(j:j)=' '
      goto 7
      endif
      goto 9

   10 continue
      if(i.ge.1000-3)stop'utword: make line shorter!!!!          '
      i=i+2
      j=i
      if(line(j:j+1).eq.mrk)stop'utword: empty string!!!!        '
   11 j=j+1
      if(j.gt.1000-1)then                 !reach the end of the line
        j=j-nbla+2
        goto 7
      endif
      if(line(j:j+1).eq.mrk)then
      line(i-2:i-1)='  '
      line(j:j+1)='  '
      goto 7
      endif
      if(line(j:j).eq.' ')then
        nbla=nbla+1
      else
        nbla=2
      endif
      goto 11

    7 j=j-1
      
      if(line(i:i+1).eq.'::'.and.line(j-1:j).eq.'::')then
        line(i:i+1)='  '
        line(j-1:j)='  '
        i=i+2
        j=j-2
        goto 7777
      endif      
      
      if(iqu.ne.-1)call expand(line,i,j) 
      do k=i,j
        if(k.gt.1)then
          if(line(k-1:k).eq.'::')then
            line(k-1:k)='  '
            j=k-2
            goto 12
          endif
        endif
      enddo
 12   continue     
      if(iqu.ne.-1)call expand(line,i,j) 

7777  continue      
      return

9999  close(ifopx)
      nopen=nopen-1
      if(nopen.eq.0.and.iprmpt.eq.-1)iprmpt=1
      goto 5
      end
      
      subroutine expand(line,i,j)
#include "aaa.h"
      character line*(*)
      if(ndefine.le.0)return
      imax=len(line) 
      imx=imax
      do while(imx.gt.2.and.line(imx:imx).eq.' ')
        imx=imx-1
      enddo
      imx=imx+1
      j0=0
      jj=0
      do ndf=1,ndefine
        l1=l1define(ndf)
        l2=l2define(ndf)
        do i0=i,j+1-l1
          if(line(i0:i0-1+l1).eq.w1define(ndf)(1:l1))then
           if(line(i0-1:i0-1).ne.'\'.and.line(i0-1:i0-1).ne.'/')then
            jj=1 
            if(l2.eq.l1)then
              line(i0:i0-1+l1)=w2define(ndf)(1:l2)
            elseif(l2.lt.l1)then
              line(i0:i0+l2-1)=w2define(ndf)(1:l2)
              do k=i0+l2,i0-1+l1
                line(k:k)=' '
              enddo
            elseif(l2.gt.l1)then
              ier=0
              do k=i0+l1,i0+l2
                if(line(k:k).ne.' ')ier=1
              enddo
              if(ier.eq.1)then
                !write(ifmt,'(9a)')('-----------------',m=1,4)
                !write(ifmt,*)'c.f:expand'
                !write(ifmt,*)'  i j i0 = ',i,j,i0
                !write(ifmt,*)'  no space for `define` replacement.'
                !write(ifmt,*)'  ',line(i0:i0-1+l2+1),l1
                !write(ifmt,*)'  ',w2define(ndf)(1:l2),l2
                k1=l2-l1
                do k=imx,i0+l1,-1
                  line(k+k1:k+k1)=line(k:k)
                  line(k:k)=' '
                enddo
                !write(ifmt,*)'  ',line(i0:i0-1+l2+1+k1)
                !write(ifmt,'(9a)')('-----------------',m=1,4)
                !stop'exit'
              endif
              line(i0:i0+l2-1)=w2define(ndf)(1:l2)
              j=i0+l2-1
            endif
           elseif(line(i0-1:i0-1).eq.'\')then
            line(i0-1:i0-1)='/'
           elseif(line(i0-1:i0-1).eq.'/')then
            line(i0-1:i0-1)=' '
           endif
          endif
        enddo
      enddo
      if(jj.eq.1)then
        do k=i,j
          if(line(k:k).ne.' ')j0=j
        enddo
        j=j0
      endif
      end


