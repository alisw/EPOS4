C
C  This file is part of EPOS4
C  Copyright (C) 2022 research institutions and authors (See CREDITS file)
C  This file is distributed under the terms of the GNU General Public License version 3 or later
C  (See COPYING file for the text of the licence)
C

      subroutine aasetspherio
#include "sf.h"
      fnztr='zzz.ztr '
      nfnztr=index(fnztr,' ')-1
      end
      subroutine spherioptns(line,j)
#include "aaa.h"
#include "sf.h"
      character line*160,linex*160
    1 call utword(line,i,j,1)
           if(line(i:j).eq.'endoptns')then
      call utword(line,i,j,0)
      if(line(i:j).eq.'spherio')return
      if(line(i:j).eq.'epos')stop 'endoptns mismatch!'
      if(line(i:j).eq.'nexus')stop 'endoptns mismatch!'
      if(line(i:j).eq.'nonhydro')stop 'endoptns mismatch!'
           elseif(line(i:j).eq.'set')then
      call utword(line,i,j,0)
      call utword(line,i,j,0)
           elseif(line(i:j).eq.'switch')then
      call utword(line,i,j,0)
      call utword(line,i,j,0)
           elseif(line(i:j).eq.'fname')then
      call utword(line,i,j,0)
      linex=line
      ix=i
      jx=j
      call utword(line,i,j,0)
      if(linex(ix:jx).eq.'ztr')then
       ! fnztr(1:j-i+1)=line(i:j)
       ! nfnztr=j-i+1
      endif
           else
      write(ifmt,'(a,a,a)')'spherio optns "',line(i:j),'" not found'
      j=160
      stop
           endif
      i=j+1
      goto 1
      end
      subroutine spheriox(nfr)
#include "aaa.h"
      if(ispherio.ne.0)then
        nfr2=nfr
        stop'\n\n    SPHERIO not available \n\n'
      endif
      end
      subroutine spherio(nrevt,nfr)
        nrevt2=nrevt
        nfr2=nfr
      end
