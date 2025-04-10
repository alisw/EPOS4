!-------------------------------------------------------------------------
!
!    This piece of code is meant to be changed by the user
!
!    The code is useful when events are read from root files.
!    Then epos should be run as
!    
!       ${EPO}epos -cproot -scr -eee ja 35g - p70 1
!
!    to read the files in the directory werner/public/root/p70ja35g 
!-------------------------------------------------------------------------

!~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
!   the function idoz1 activates (or not) the other functions in this file
!~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

      integer function idoz1()
      idoz1 = 0   ! -> 1  to activate the following routines 
      end

!~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
! the routines transfer_head and transfer_event are called 
! after reading an event, providing the complete event information via 
! the arguments of these routines
!~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

      subroutine transfer_head(iversn,laproj,maproj,latarg,matarg,engy)
  !~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
  ! 
  !~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
      common/ctransferhead/
     .iversnx,laprojx,maprojx,latargx,matargx,engyx
      iversnx=iversn
      laprojx=laproj
      maprojx=maproj
      latargx=latarg
      matargx=matarg
      engyx  =engy
      end

      subroutine transfer_event(izmode,cnt,np
     .     ,nev,npt,ngl,kol,phi,phir
     .     ,psi2,psi3,psi4,psi5,ecci2,ecci3,ecci4,ecci5
     .   ,id,ist,ity,ior,jor,px,py,pz,en,am,x,y,z,t)      
     
  !~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
  ! izmode, cnt, nptevt refer to simple variables,
  ! all others to arrays with particle properties. Array dimension : np
  ! izmode is an integer defing the meaning of the centrality variable cnt :
  !    1 : impact parameter
  !    2 : number of Pomerons
  ! np is the number of particles in the event.
  ! Meaning of the arrays elements :
  !  id(i) .................... id of particle i
  !  ist(i) ................... status and particle type of particle i
  !  ity(i) ................... type of particle origin of particle i
  !  ior(i) ................... mother of particle i
  !  jor(i) ................... father of particle i
  !  px(i), py(i), pz(i) ...... momentum of particle i
  !  en(i) .................... energy  of particle i
  !  am(i) .................... mass  of particle i
  !  x(i), y(i), z(i), t(i) ... space-time of particle i
  !~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ 
      
      common/ctransferhead/
     .iversnx,laprojx,maprojx,latargx,matargx,engyx
      integer id(np),ist(np),ity(np),ior(np),jor(np)
      real px(np),py(np),pz(np),en(np),am(np)
      real x(np),y(np),z(np),t(np)   
      data ncntz1/0/
      save ncntz1
      ncntz1=ncntz1+1
      if(idoz1().ne.1)return

      dummy=izmode+cnt+np
     .     +nev+npt+ngl+kol+phi+phir
     .     +psi2+psi3+psi4+psi5+ecci2+ecci3+ecci4+ecci5
     .   +id(1)+ist(1)+ity(1)+ior(1)+jor(1)+px(1)+py(1)+pz(1)+en(1)
     .   +am(1)+x(1)+y(1)+z(1)+t(1) !avoid useless warnings

      if(ncntz1.eq.1)then
      print*,'******** first event ********'
      print*,iversnx,laprojx,maprojx,latargx,matargx,engyx
      endif
      
      !***** put your code here *****
      
      if(ncntz1.eq.160)then
      print*,'---------- just a test to see that it works ------------'
      print*,'ix bim np:',izmode,bm,np
      print*,'nev npt ngl kol :',nev,npt,ngl,kol
      print*,'phi phir:',phi,phir
      print*,'psi2-5 :',psi2,psi3,psi4,psi5
      print*,'ecci2-5 :',ecci2,ecci3,ecci4,ecci5
      do i=1,np
        j=iabs(id(i))
c      if(j.eq.140.or.j.eq.240)then
c      if(j.eq.1120)then
         print*,'particle:',ior(i),jor(i),i,id(i),ist(i),ity(i)  
     . !,px(i),py(i),pz(i)  !,en(i),am(i),x(i),y(i),z(i),t(i)
c      endif
      enddo  
      endif
      if(ncntz1.eq.160)stop
      
      end

!~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
! the routine finish_program is called at the end of the run
!~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
  
      subroutine finish_program()
      if(idoz1().ne.1)return
     
      !***** put your code here *****

      print*,'******** finish_program ********'
      
      end
