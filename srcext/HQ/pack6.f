
      subroutine reseteventlist()
      use PBGgenvar
      implicit none
      integer ievent
      firstevent=-1
      firstemptyevent=1
      do ievent=2,neventmax
        eventinfo_next(ievent-1)=ievent
      enddo 
      return
      end


      subroutine resetparticlelist()
      use PBGgenvar
      implicit none
      integer ipart
      npart=0
      firstpart=-1
      firstemptypart=1
      do ipart=2,npartmax
        partinfo_next(ipart-1)=ipart
      enddo
      partinfo_next(npartmax)=-1
      return
      end


      subroutine addevent(firstnewevent,x_event,tb,delt,scatcode,i,j)
      use PBGgenvar
      implicit none
      integer firstnewevent,scatcode,i,j,k
      double precision x_event(0:3),tb,delt
      integer oldfirst
      oldfirst=firstnewevent
      firstnewevent=firstemptyevent
      if(firstemptyevent.eq.neventmax) then
         write(6,*) 'neventmax is too small'
         stop
      endif
      firstemptyevent=eventinfo_next(firstemptyevent)
      do k=0,3 
        eventinfo_r(firstnewevent,k)=x_event(k)
      enddo
      eventinfo_tb(firstnewevent)=tb
      eventinfo_delt(firstnewevent)=delt
      if(itypeorder.eq.1) then
         eventinfo_tord(firstnewevent)=x_event(0)
      endif
      if(itypeorder.eq.2) then
         eventinfo_tord(firstnewevent)=tb
      endif
      eventinfo_scatcode(firstnewevent)=scatcode
      eventinfo_next(firstnewevent)=oldfirst
      eventinfo_i(firstnewevent)=i
      if(scatcode.lt.0.and.i.ne.j) then
         write(6,*) 'DANGER event needs i=j'
         stop
      endif
      eventinfo_j(firstnewevent)=j
      return
      end


      subroutine event_chain_test()
      use PBGgenvar
      implicit none
      integer ievent,ii
      double precision tord1,tord2
 150  format(i5,i5,E12.5,i5,i5,i4,i4,i4)
      if(firstevent.eq.-1) then
        write(6,*) 'no event in list'
        return
      endif
      ievent=firstevent
      ii=0
      tord1=0.
      do while(ievent.ne.-1)
       ii=ii+1
       if (ii.le.20) then       
         write(6,150) ii,ievent,eventinfo_tord(ievent),
     &   eventinfo_i(ievent),eventinfo_j(ievent),
     &   partinfo_ires(eventinfo_i(ievent)),
     &   partinfo_ires(eventinfo_j(ievent)),
     &   eventinfo_scatcode(ievent)
       endif
       tord2=eventinfo_tord(ievent)
       if(tord2.lt.tord1) then
          write(6,*) 'ordering problem in event chain'
          stop
       else
          tord1=tord2
       endif   
       ievent=eventinfo_next(ievent)
      enddo
      return
      end


      subroutine part_chain_test(nbaff)
      use PBGgenvar
      implicit none 
      integer partcount,prevpart,ipart,nbaff
      double precision masscheck
 150  format(i5,i5,i5,' ',E12.5,' ',E12.5,' ',E12.5,' ',E12.5,' ',E12.5,
     &   ' ',i3)
 151  format(E12.5,' ',E12.5,' ',E12.5,' ',E12.5,' ',E12.5)
      ipart=firstpart
      prevpart=-1
      partcount=0
      write(6789,*) 'entering part_chain_test' 
      do while (ipart.NE.-1)
        write(6789,*) ipart
        partcount=partcount+1
        masscheck=partinfo_p(ipart,0)**2-partinfo_p(ipart,1)**2
     &   -partinfo_p(ipart,2)**2-partinfo_p(ipart,3)**2-
     &    partinfo_mass(ipart)**2   
        if(dabs(masscheck).gt.1d-07) then 
          write(6789,*) 'mass check problem:', ipart,
     &    partinfo_ires(ipart),masscheck
          stop
        endif
        if(partcount.le.nbaff) then
           write(6789,150) partcount,ipart,partinfo_ires(ipart),
     &      partinfo_tb(ipart),partinfo_r(ipart,0),partinfo_r(ipart,1),
     &      partinfo_r(ipart,2),partinfo_r(ipart,3),
     &      partinfo_whichmedium(ipart)
           write(6789,151) partinfo_p(ipart,0),partinfo_p(ipart,1),
     &      partinfo_p(ipart,2),partinfo_p(ipart,3),partinfo_mass(ipart)
        endif 
        ipart=partinfo_next(ipart)
      enddo
      if(npart.ne.partcount) then
        write(6789,*) 'FAILURE: npart=',npart,'partcount=',partcount
        stop
      else
      endif
      return
      end


      subroutine part_chain_test2()
      use PBGgenvar
      implicit none 
      integer inchain(npartmax)
      integer partcount(-1:1),ipart,k
      do ipart=1,npartmax
        inchain(ipart)=0
      enddo
      do k=-1,1
        partcount(k)=0
      enddo
      ipart=firstpart
      do while (ipart.NE.-1)
        inchain(ipart)=1
        ipart=partinfo_next(ipart)
      enddo
      ipart=firstemptypart
      do while (ipart.NE.-1)
        inchain(ipart)=-1
        ipart=partinfo_next(ipart)
      enddo
      do ipart=1,npartmax
        partcount(inchain(ipart))=partcount(inchain(ipart))+1
      enddo
      if((partcount(0).ne.0).or.(partcount(1).ne.npart).or.
     &    ((partcount(1)+partcount(-1)).ne.npartmax)) then
        write(6,*) 'FAILURE in part chain test:',(partcount(k),k=-1,1)
        stop
      endif
      return
      end


