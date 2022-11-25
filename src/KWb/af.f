C
C  This file is part of EPOS4
C  Copyright (C) 2022 research institutions and authors (See CREDITS file)
C  This file is distributed under the terms of the GNU General Public License version 3 or later
C  (See COPYING file for the text of the licence)
C


c#############################################################################
c                   Managing particle overflow (i > mxptl)
c#############################################################################
c particles  are  dumped to  and restored from  cccptl  (C++ object)
c dump partile i to cccptl[i-1-mxptl] 
c restore partile i from cccptl[i-1-mxptl]  
c#############################################################################

c-----------------------------------------------------------------------
c Before calling dumpcccptl,restorecccptl  
c one needs  
c                            call createcccptl(n)  via    call checkcccptl
c and at the end
c                            call destroycccptl(n)
c 
c-----------------------------------------------------------------------

      subroutine dumpcccptl(j,i)
#include "aaa.h"
      if(i.le.mxptl)stop'########## ERROR 12082014a1 ###########'
      if(icccptl.eq.0)stop'########## ERROR 12082014a2 ###########'
      call cccptldump(idptl(j),istptl(j),ityptl(j)
     .,iorptl(j),jorptl(j),ifrptl(1,j),ifrptl(2,j)
     .,tivptl(1,j),tivptl(2,j)
     .,pptl(1,j),pptl(2,j),pptl(3,j), pptl(4,j),pptl(5,j)
     .,xorptl(1,j),xorptl(2,j),xorptl(3,j),xorptl(4,j),radptl(j)
     .,desptl(j),dezptl(j),qsqptl(j),zpaptl(1,j),zpaptl(2,j),rinptl(j)
     .,ibptl(1,j),ibptl(2,j),ibptl(3,j),ibptl(4,j),iaaptl(j),itsptl(j)
     .,i-1-mxptl)
      end
      
      subroutine restorecccptl(i,j)
      real triggcrash(2) 
#include "aaa.h"
      if(i.le.mxptl)stop'########## ERROR 12082014b1 ###########'
      if(icccptl.eq.0)then
        triggcrash(i)=1 !to trigger a crash 
        stop'########## ERROR 12082014b2 ###########'
      endif
      call cccptlrestore(idptl(j),istptl(j),ityptl(j)
     .,iorptl(j),jorptl(j),ifrptl(1,j),ifrptl(2,j)
     .,tivptl(1,j),tivptl(2,j)
     .,pptl(1,j),pptl(2,j),pptl(3,j), pptl(4,j),pptl(5,j)
     .,xorptl(1,j),xorptl(2,j),xorptl(3,j),xorptl(4,j),radptl(j)
     .,desptl(j),dezptl(j),qsqptl(j),zpaptl(1,j),zpaptl(2,j),rinptl(j)
     .,ibptl(1,j),ibptl(2,j),ibptl(3,j),ibptl(4,j),iaaptl(j),itsptl(j)
     .,i-1-mxptl)
      end
      
      subroutine checkcccptl(n)            
#include "aaa.h"
      if(n.le.mxptl)return
      if(icccptl.eq.0)then
        call memo(1,'create cccptl object;')
        icccptl=1
        call maxsize_get(2, mx2ptl ) 
        call createcccptl(mx2ptl)
        call memo(2,';')
      endif
      if(n.gt.mxptl+mx2ptl)
     .stop'########## ERROR 12082014c ###########'
      end

      subroutine closecccptl     
#include "aaa.h"
      if(icccptl.eq.1)then
        call memo(1,'destroy cccptl object;')
        call maxsize_get(2, mx2ptl ) 
        call destroycccptl(mx2ptl)
        icccptl=0
        call memo(2,';')
      endif  
      end     
      
      !-----------------------------------------------------
      !     set
      !-----------------------------------------------------

      subroutine setnptl(ival)
#include "aaa.h"
      nptl=ival
      end
     
      subroutine setistptl(i,ival)
#include "aaa.h"
      if(i.le.mxptl)then
        istptl(i)=ival
      else
        call checkcccptl(i)
        call istptlset(i-1-mxptl,ival)
      endif
      end

      subroutine setidptl(i,ival)            
#include "aaa.h"
      if(i.le.mxptl)then
        idptl(i)=ival
      else
        call checkcccptl(i)
        call idptlset(i-1-mxptl,ival)
      endif
      end

      subroutine setityptl(i,ival)            
#include "aaa.h"
      if(i.le.mxptl)then
        ityptl(i)=ival
      else
        call checkcccptl(i)
        call ityptlset(i-1-mxptl,ival)
      endif
      end

      subroutine setiorptl(i,ival)            
#include "aaa.h"
      if(i.le.mxptl)then
        iorptl(i)=ival
      else
        call checkcccptl(i)
        call iorptlset(i-1-mxptl,ival)
      endif
      end

      subroutine setjorptl(i,ival)            
#include "aaa.h"
      if(i.le.mxptl)then
        jorptl(i)=ival
      else
        call checkcccptl(i)
        call jorptlset(i-1-mxptl,ival)
      endif
      end

      subroutine setifrptl(i,i1,i2)            
#include "aaa.h"
      if(i.le.mxptl)then
        ifrptl(1,i)=i1
        ifrptl(2,i)=i2
      else
        call checkcccptl(i)
        call ifrptlset(i-1-mxptl,i1,i2)
      endif
      end

      subroutine setibptl(i,i1,i2,i3,i4)            
#include "aaa.h"
      if(i.le.mxptl)then
        ibptl(1,i)=i1
        ibptl(2,i)=i2
        ibptl(3,i)=i3
        ibptl(4,i)=i4
      else
        call checkcccptl(i)
        call ibptlset(i-1-mxptl,i1,i2,i3,i4)
      endif
      end
      subroutine set2ibptl(i,j,ival)            
#include "aaa.h"
      if(i.le.mxptl)then
        ibptl(j,i)=ival
      else
        call checkcccptl(i)
        call ibptlset2(i-1-mxptl,j,ival)
      endif
      end

      subroutine setdesptl(i,val)            
#include "aaa.h"
      if(i.le.mxptl)then
        desptl(i)=val
      else
        call checkcccptl(i)
        call desptlset(i-1-mxptl,val)
      endif
      end

      subroutine setradptl(i,val)            
#include "aaa.h"
      if(i.le.mxptl)then
        radptl(i)=val
      else
        call checkcccptl(i)
        call radptlset(i-1-mxptl,val)
      endif
      end
      
      subroutine setrinptl(i,val)            
#include "aaa.h"
      if(i.le.mxptl)then
        rinptl(i)=val
      else
        call checkcccptl(i)
        call rinptlset(i-1-mxptl,val)
      endif
      end
      
      subroutine setpptl(i,p1,p2,p3,p4,p5)
#include "aaa.h"
      if(i.le.mxptl)then
         pptl(1,i)=p1
         pptl(2,i)=p2
         pptl(3,i)=p3
         pptl(4,i)=p4
         pptl(5,i)=p5
      else
        call checkcccptl(i)
        call pptlset(i-1-mxptl,p1,p2,p3,p4,p5)
      endif
      end

      subroutine setxorptl(i,x1,x2,x3,x4)
#include "aaa.h"
      if(i.le.mxptl)then
         xorptl(1,i)=x1
         xorptl(2,i)=x2
         xorptl(3,i)=x3
         xorptl(4,i)=x4
      else
        call checkcccptl(i)
        call xorptlset(i-1-mxptl,x1,x2,x3,x4)
      endif
      end

      subroutine settivptl(i,t1,t2)
#include "aaa.h"
      if(i.le.mxptl)then
         tivptl(1,i)=t1
         tivptl(2,i)=t2
      else
        call checkcccptl(i)
        call tivptlset(i-1-mxptl,t1,t2)
      endif
      end
      
      subroutine setzpaptl(i,t1,t2)
#include "aaa.h"
      if(i.le.mxptl)then
         zpaptl(1,i)=t1
         zpaptl(2,i)=t2
      else
        call checkcccptl(i)
        call zpaptlset(i-1-mxptl,t1,t2)
      endif
      end

      !-----------------------------------------------------
      !     get
      !-----------------------------------------------------

      subroutine getnptl(ival)
#include "aaa.h"
      ival=nptl
      end

      subroutine getistptl(i,ival)
#include "aaa.h"
      if(icccptl.eq.0.or.i.le.mxptl)then
        !idummy=1/i !to force crash with backtrace in case of i=0 
        ival=istptl(i)
      elseif(icccptl.eq.1)then
        call istptlget(i-1-mxptl,ival)
      else
        stop'########## ERROR 12082014e ###########'
      endif
      end

      subroutine getidptl(i,ival)
#include "aaa.h"
      if(icccptl.eq.0.or.i.le.mxptl)then
         ival=idptl(i)
      elseif(icccptl.eq.1)then
        call idptlget(i-1-mxptl,ival)
      else
        stop'########## ERROR 12082014f ###########'
      endif
      end

      subroutine getityptl(i,ival)            
#include "aaa.h"
      if(icccptl.eq.0.or.i.le.mxptl)then
         ival=ityptl(i)
      elseif(icccptl.eq.1)then
        call ityptlget(i-1-mxptl,ival)
      else
        stop'########## ERROR 12082014h ###########'
      endif
      end

      subroutine getiorptl(i,ival)            
#include "aaa.h"
      if(icccptl.eq.0.or.i.le.mxptl)then
         ival=iorptl(i)
      elseif(icccptl.eq.1)then
        call iorptlget(i-1-mxptl,ival)
      else
        stop'########## ERROR 12082014h2 ###########'
      endif
      end

      subroutine getjorptl(i,ival)            
#include "aaa.h"
      if(icccptl.eq.0.or.i.le.mxptl)then
         ival=jorptl(i)
      elseif(icccptl.eq.1)then
        call jorptlget(i-1-mxptl,ival)
      else
        stop'########## ERROR 22092016a ###########'
      endif
      end

      subroutine getifrptl(i,i1,i2)            
#include "aaa.h"
      if(icccptl.eq.0.or.i.le.mxptl)then
         i1=ifrptl(1,i)
         i2=ifrptl(2,i)
      elseif(icccptl.eq.1)then
        call ifrptlget(i-1-mxptl,i1,i2)
      else
        stop'########## ERROR 22092016b ###########'
      endif
      end

      subroutine getibptl(i,i1,i2,i3,i4)            
#include "aaa.h"
      if(icccptl.eq.0.or.i.le.mxptl)then
         i1=ibptl(1,i)
         i2=ibptl(2,i)
         i3=ibptl(3,i)
         i4=ibptl(4,i)
      elseif(icccptl.eq.1)then
        call ibptlget(i-1-mxptl,i1,i2,i3,i4)
      else
        stop'########## ERROR 22092016b ###########'
      endif
      end
      subroutine get2ibptl(i,j,ival)            
#include "aaa.h"
      if(icccptl.eq.0.or.i.le.mxptl)then
         ival=ibptl(j,i)
      elseif(icccptl.eq.1)then
        call ibptlget2(i-1-mxptl,j,ival)
      else
        stop'########## ERROR 22092016b ###########'
      endif
      end

      subroutine getdesptl(i,val)              
#include "aaa.h"
      if(icccptl.eq.0.or.i.le.mxptl)then
         val=desptl(i)
      elseif(icccptl.eq.1)then
        call desptlget(i-1-mxptl,val)
      else
        stop'########## ERROR 15082014c ###########'
      endif
      end

      subroutine getradptl(i,val)              
#include "aaa.h"
      if(icccptl.eq.0.or.i.le.mxptl)then
         val=radptl(i)
      elseif(icccptl.eq.1)then
        call radptlget(i-1-mxptl,val)
      else
        stop'########## ERROR 15082014c ###########'
      endif
      end

      subroutine getrinptl(i,val)
#include "aaa.h"
      if(icccptl.eq.0.or.i.le.mxptl)then
        val=rinptl(i)
      elseif(icccptl.eq.1)then
        call rinptlget(i-1-mxptl,val)
      else
        stop'########## ERROR 23112020 ###########'
      endif
      end

      subroutine getpptl(i,p1,p2,p3,p4,p5)
#include "aaa.h"
      if(icccptl.eq.0.or.i.le.mxptl)then
         p1=pptl(1,i)
         p2=pptl(2,i)
         p3=pptl(3,i)
         p4=pptl(4,i)
         p5=pptl(5,i)
      elseif(icccptl.eq.1)then
        call pptlget(i-1-mxptl,p1,p2,p3,p4,p5)
      else
        stop'########## ERROR 12082014i ###########'
      endif
      end

      subroutine getxorptl(i,x1,x2,x3,x4)
#include "aaa.h"
      if(icccptl.eq.0.or.i.le.mxptl)then
         x1=xorptl(1,i)
         x2=xorptl(2,i)
         x3=xorptl(3,i)
         x4=xorptl(4,i)
      elseif(icccptl.eq.1)then
        call xorptlget(i-1-mxptl,x1,x2,x3,x4)
      else
        stop'########## ERROR 22102014 ###########'
      endif
      end

      subroutine gettivptl(i,t1,t2)
#include "aaa.h"
      if(icccptl.eq.0.or.i.le.mxptl)then
         t1=tivptl(1,i)
         t2=tivptl(2,i)
      elseif(icccptl.eq.1)then
        call tivptlget(i-1-mxptl,t1,t2)
      else
        stop'########## ERROR 22102014 ###########'
      endif
      end

      subroutine getzpaptl(i,t1,t2)
#include "aaa.h"
      if(icccptl.eq.0.or.i.le.mxptl)then
         t1=zpaptl(1,i)
         t2=zpaptl(2,i)
      elseif(icccptl.eq.1)then
        call zpaptlget(i-1-mxptl,t1,t2)
      else
        stop'########## ERROR 22102014b ###########'
      endif
      end

