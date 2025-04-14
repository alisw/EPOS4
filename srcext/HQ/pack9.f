      subroutine aacharm2
      use JA2 
      use PBGevcount 
      use PBGhquarks 
      use PBGnpartmax 
      use PBGnumberevent 
      use PBGblckifptcles






#include "aaa.h"       

      common/hqco/hqicount,iwrtja
      data hqicount/0/
      integer nbfinpartmax
      parameter (nbfinpartmax=50000)
      common/corr/finalp(5,nbfinpartmax),ifinalid(nbfinpartmax),
     & ifinpart, 
     & finalphq(5,2*npartmax),ifinalidhq(2*npartmax),
     & ifinaltyphq(2*npartmax),ifinparthqep,ifinparthq,
     & ifinaliorhq(2*npartmax),ifinalisthq(2*npartmax),
     & ahqp(2*npartmax,5),iahid(2*npartmax),iahiior(2*npartmax),
     & iahisto(2*npartmax),ifinaltyphqep(2*npartmax),
     & finalphqep(5,2*npartmax),ifinalidhqep(2*npartmax),
     & ihard(2*npartmax,4)
       

      call utpri('charm2',ish,ishini,4)

      if(jaish.ge.1)write(ifmt,'(a,i8,i3)')
     .'start aacharm2: nptl icccptl =',nptl,icccptl
 

        If(hqicount.ne.0) goto 299  
299   continue
      eventshq=eventshq+1
      hqicount=hqicount+1  

        ifinpart=0
        ifinparthq=0
        ifinparthqep=0

        icbcounter=0
       do ip=1,nptl


          call getidptl(ip,id)
       call getistptl(ip,ist)
       call getityptl(ip,ity)
       if(ity.eq.61)id=idtrafo('pdg','nxs',id) !to get epos id
       call getiorptl(ip,iior)! ip of mother
       call getpptl(ip,p1,p2,p3,p4,p5) 
         
         

       if(ist.eq.0)then
          call getiorptl(ip,kor)
          if(kor.lt.1.or.kor.gt.nptl) goto212
          call getidptl(kor,idk)
          call getityptl(kor,ityk)

       if(ityk.eq.61)idk=idtrafo('pdg','nxs',idk)
       imotherc=0
       call idflav(idk,if1,if2,if3,idu,0)
       if(abs(if1).eq.4.or.abs(if1).eq.5)imotherc=1
       if(abs(if2).eq.4.or.abs(if2).eq.5)imotherc=1
       if(abs(if3).eq.4.or.abs(if3).eq.5)imotherc=1
212    iacc=0
       idabs=abs(id)
       
       call idchrg( 10 ,idabs,ch) 

        if(abs(ch).gt.0.1)iacc=1


        if(iacc.eq.1.and.ifptcles) then 
       ifinpart=ifinpart+1
       if(ifinpart.gt.nbfinpartmax) then
          write(8,*) '#of partclcles too large in aacharm2'
          close(8)
          stop
       else
 
       call getpptl(ip,p1,p2,p3,p4,p5) 
       finalp(1,ifinpart)=p1
       finalp(2,ifinpart)=p2
       finalp(3,ifinpart)=p3
       finalp(4,ifinpart)=p4
       finalp(5,ifinpart)=p5
       ifinalid(ifinpart)=id
       endif
       endif
       endif


       call getifrptl(ip,ifr1,ifr2)
       inodaugb=0
       if(ifr2.gt.0)then
          do j0=ifr1,ifr2
         call getidptl(j0,idj0)
         call getityptl(j0,ityj0)
         if(ityj0.eq.61)idj0=idtrafo('pdg','nxs',idj0)
         call idflav(idj0,if1,if2,if3,idu,0)
         if(abs(if1).eq.5)inodaugb=1
         if(abs(if2).eq.5)inodaugb=1
         if(abs(if3).eq.5)inodaugb=1
        enddo 
       endif

       inodaugc=0
       if(ifr2.gt.0)then
        do j0=ifr1,ifr2
         call getidptl(j0,idj0)
         call getityptl(j0,ityj0)
         if(ityj0.eq.61)idj0=idtrafo('pdg','nxs',idj0)
         call idflav(idj0,if1,if2,if3,idu,0)
         if(abs(if1).eq.4)inodaugc=1
         if(abs(if2).eq.4)inodaugc=1
         if(abs(if3).eq.4)inodaugc=1

        enddo 
       endif
 
       ido=0
       isto=0
       ityo=0 

       if(iior.gt.1)then
        call getidptl(iior,ido)
        call getistptl(iior,isto)
        call getityptl(iior,ityo)
        if(ityo.eq.61)ido=idtrafo('pdg','nxs',ido)
        call getiorptl(iior,iioor)
       else
          iioor=-1
       endif

       idoo=0
       istoo=0
       ityoo=0 
 
       if(abs(id).eq.12.and.jaish.ge.9)write(86,*)ip,
     & isto,ityo,ido,'x'

       if(icccptl.eq.0.and.iioor.gt.mxptl)then
        write(ifmt,'(a,i3,4i8)')'ERROR aacharm2'
     . ,icccptl,iioor,iior,ip,nptl
        stop'####### ERROR 10032017 in aacharm2 #######'
       endif

       if(iioor.gt.1)then
       call getidptl(iioor,idoo)
       call getistptl(iioor,istoo)
       call getityptl(iioor,ityoo)
       if(ityoo.eq.61)idoo=idtrafo('pdg','nxs',idoo)
       endif

 
       call getpptl(ip,p1,p2,p3,p4,p5) 
      if(jaish.ge.9)write(ifch,'(a,i8,a,i8,a,i9,a,i3,a,i3,$)')
     & '19  ior:',iior,'  ip:',ip,'  id:',id,'  ist:',ist,'  ity:',ity
      if(jaish.ge.9)write(ifch,'(a,i9,a,i3,a,i3,$)')
     & '  ido:',ido,'  isto:',isto,'  ityo:',ityo
      if(jaish.ge.9)write(ifch,'(a,i9,a,i9,a,i3,a,i3,$)')
     & '  ioor:',iioor,'  idoo:',idoo,'  istoo:',istoo
     & ,'  ityoo:',ityoo
      if(jaish.ge.9)write(ifch,'(a,5e11.3)')
     & '  p:',p1,p2,p3,p4,p5



        ityp=0
          call idflav(id,if1,if2,if3,idu,0)
          if(abs(if1).eq.5)ityp=5
          if(abs(if2).eq.5)ityp=5
          if(abs(if3).eq.5)ityp=5
          if(abs(if1).eq.4)ityp=2
          if(abs(if2).eq.4)ityp=2
          if(abs(if3).eq.4)ityp=2
          if (id.eq.4.or.id.eq.-4)ityp=1
          if (id.eq.5.or.id.eq.-5)ityp=4

      idtyp=0
      if(ist.eq.21.and.ityp.eq.1)idtyp=1
      if((ist.eq.6.or.ist.eq.8).and.ity.eq.60.and.ityp.eq.2.
     &     and.inodaugc.eq.0)idtyp=3
      if((ist.eq.1.or.ist.eq.0).and.ity.eq.60.and.ityp.eq.2.
     &     and.inodaugc.eq.0)idtyp=4
      if((ist.eq.6.or.ist.eq.8).and.ity.lt.60.and.ityp.eq.2.
     &     and.inodaugc.eq.0)idtyp=5
      if(ist.eq.2.and.ityp.eq.2.and.inodaugc.eq.0)then
         idtyp=6
      endif
 
      if(ist.eq.21.and.ityp.eq.4)idtyp=11
      if((ist.eq.6.or.ist.eq.8).and.ity.eq.60.and.ityp.eq.5.
     &     and.inodaugb.eq.0)idtyp=13
      if((ist.eq.1.or.ist.eq.0).and.ity.eq.60.and.ityp.eq.5.
     &     and.inodaugb.eq.0)idtyp=14
      if((ist.eq.6.or.ist.eq.8).and.ity.lt.60.and.ityp.eq.5.
     &     and.inodaugb.eq.0)idtyp=15
      if(ist.eq.2.and.ityp.eq.5.
     &     and.inodaugb.eq.0)idtyp=16

      if(idtyp.eq.1.or.idtyp.eq.11) then
      icbcounter=icbcounter+1
      call getpptl(ip,p1,p2,p3,p4,p5) 
      ahqp(icbcounter,1)=p1
      ahqp(icbcounter,2)=p2
      ahqp(icbcounter,3)=p3
      ahqp(icbcounter,4)=p4
      ahqp(icbcounter,5)=p5
      iahid(icbcounter)=id
      iahiior(icbcounter)=iior
      iahisto(icbcounter)=isto

      do iba=1,100
        iback=ip-iba
        call getistptl(iback,istb)        
        
        if(istb.eq.25)goto 964
        
        enddo  
964   call getidptl(iback,ihard0)
        call getidptl(iback-1,ihard1)
        call getistptl(iback-1,istb1) 
        ihard(icbcounter,1)=ihard0       
        ihard(icbcounter,2)=ihard1 
        ihard(icbcounter,3)=istb       
        ihard(icbcounter,4)=istb1 

      endif 


      if(idtyp.eq.0) goto 311

      it=0
      if(ity.lt.60)it=1
      if(ity.eq.60)it=2 
       if(idtyp.eq.4.or.idtyp.eq.3.or.idtyp.eq.13.
     & or.idtyp.eq.14)then
       ifinparthq=ifinparthq+1
       if(ifinparthq.gt.2*npartmax) then
          write(8,*) 'dimension of HQ exceeded in aacharm2'
          write(8,*) ifinparthq 
          close(8)
          stop
       endif
       finalphq(1,ifinparthq)=p1
       finalphq(2,ifinparthq)=p2
       finalphq(3,ifinparthq)=p3
       finalphq(4,ifinparthq)=p4
       finalphq(5,ifinparthq)=p5
       ifinalidhq(ifinparthq)=id
       ifinaltyphq(ifinparthq)=idtyp

       endif 
       if(idtyp.eq.6)then
       
       ifinparthqep=ifinparthqep+1
       if(ifinparthqeq.gt.2*npartmax) then
          write(8,*) 'dimension of HQ exceeded in aacharm2'
          write(8,*) ifinparthqep 
          close(8)
          stop
       endif
       finalphqep(1,ifinparthqep)=p1
       finalphqep(2,ifinparthqep)=p2
       finalphqep(3,ifinparthqep)=p3
       finalphqep(4,ifinparthqep)=p4
       finalphqep(5,ifinparthqep)=p5
       ifinalidhqep(ifinparthqep)=id
       ifinaltyphqep(ifinparthqep)=idtyp

       endif

311     continue
      enddo

        write(ifmt,*)'write on disc',ifinpart
        if(iwrtja.eq.0)then 
        open(unit=85,file=fnus2(1:nfnus2)//'2',status='unknown')
        else
        open(unit=85,file=fnus2(1:nfnus2)//'2',status="old",
     &  position="append")
        endif
       iwrtja=1
        write(85,*)'0',ish,jaish
       call getevt(nev,phi,phir,bim,egy,npt,ngl,kol)

       if(jaish.ge.8)write(85,*) bim,phi,phir,npt,ngl,kol,nev
       call getecc(psi2,psi3,psi4,psi5,ecci2,ecci3,ecci4,ecci5)       
       if(jaish.ge.8)write(85,*) psi2,psi3,psi4,psi5,ecci2,
     & ecci3,ecci4,ecci5
       if(jaish.ge.3)write(85,*)icbcounter

343    format(1x,2i7,5E13.6,i3)
3433    format(1x,2i7,5E13.6,i8,6i4)

       do kl=1,icbcounter
       if(jaish.ge.3)write(85,3433)kl,iahid(kl),
     & (ahqp(kl,kll),kll=1,5),
     & iahiior(kl),iahisto(kl),(ihard(kl,kkl),kkl=1,4)
       enddo      


       if(jaish.ge.3)write(85,*)nHQ
       do kl=1,nHQ
       if(jaish.ge.3)write(85,343)kl,idKlaus(idhquarks(kl)),
     & (phquarks(kl,kll),
     & kll=1,3),phquarks(kl,0),mhquarks(kl),2
       enddo      

 
       if(jaish.ge.1)write(85,*)ifinparthq
       do kl=1,ifinparthq
       if(jaish.ge.1)write(85,343)kl,ifinalidhq(kl),
     & (finalphq(kll,kl),kll=1,5),
     & ifinaltyphq(kl)
   
       enddo      


       if(jaish.ge.2)write(85,*)ifinparthqep
       do kl=1,ifinparthqep
       if(jaish.ge.2)write(85,343)kl,ifinalidhqep(kl),
     & (finalphqep(kll,kl),kll=1,5),
     & ifinaltyphqep(kl)
       enddo 


       if(jaish.ge.5)write(85,*)ifinpart
       do kl=1,ifinpart
       if(jaish.ge.5)write(85,343)kl,ifinalid(kl),(finalp(kll,kl),
     & kll=1,5)
       enddo      

       close(unit=85)

      deallocate(idhquarks)
      deallocate(mhquarks)
      deallocate(phquarks)

      call utprix('charm2',ish,ishini,4)

      return
      end


      subroutine aacharmclose
      use JA2
      use PBGnumberevent
      use PBGevcount

      common/hqco/hqicount,iwrtja

#include "aaa.h" 
        call utpri('charm2',ish,ishini,4) 
        if(iwrtja.eq.0)then 
        open(unit=85,file=fnus2(1:nfnus2)//'2',status='unknown')
        else
        open(unit=85,file=fnus2(1:nfnus2)//'2',status="old",
     &  position="append")
        endif
       iwrtja=1

       write(85,*)'-1'
       if(jaish.ge.1)write(85,*)eventsxx,eventshq
       if(jaish.ge.1)write(85,*)'end'
      close(85)


       call utprix('charm2',ish,ishini,4)


      return
      end
