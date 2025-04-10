C
C  This file is part of EPOS4
C  Copyright (C) 2022 research institutions and authors (See CREDITS file)
C  This file is distributed under the terms of the GNU General Public License version 3 or later
C  (See COPYING file for the text of the licence)
C

c To do :
c * since part of the diff xs is used to produce (hard) Pomerons sigcut should
c not be used to fix Q2s (but only the Part of sigine which does not produce
c Pomerons should be removed from sigine
c * the part of the cross-section not producing Pomerons could be calculated
c independantly with 1/(M2)**alppom dependance integrated from min to cumpom
c and taken from diff xs without dependance (constant fraction of it for 
c given s + x too low). Normalization could be fixedf for SD xs at low energy 
c (all of it) and by large rapidity gaps at LHC.
c * use breit wigner in idres for resonance from remnant. Minimum mass should
c be a bit less than the mass of the lighter resonance. Use resonance as 
c as default excited remanant without MPI (rexpdif=1)
c * possibility not to excite remnant only if Y diagram. SD if Y and not 
c excited. DD is CD or ND with both excited and Pom mass too low + part of 
c without Pomeron.


c---------------------------------------------------------------------
      subroutine EventType
c---------------------------------------------------------------------
c   itpr   = 0    elastic | - typevt= 0   |                          |
c   itpr   =-1      inel  | - typevt=-1   | both remn excit (P or R) |
c   itpr   = 1       DPE  | - typevt=-2   | high mass diff           |
c                   SDpro | - typevt=-3   | high mass diff           |
c                   SDtar | - typevt=-4   | high mass diff           |
c                         | - |typevt|>1  | NSD                      |
c   itpr   =-2      DD    | - typevt= 1   |     soft                 |
c                   SDpro | - typevt= 3   |     diff                 |
c                   SDtar | - typevt= 4   |     no pom               |
c   itpr = 2        CD    | - typevt= 2   | exclusive central        |
c---------------------------------------------------------------------
#include "aaa.h"
#include "ems.h"
      typevt=0                !ela
      if(maproj+matarg.eq.2)then     !pp
        if(itpr(1).ne.0)then
          anintine=anintine+1.
          if(mod(iep(1),10).eq.2.or.mod(iet(1),10).eq.2
     &                   .or.iep(1).eq.0.or.iet(1).eq.0)then
            anintdiff=anintdiff+1.
            if((iep(1).eq.0.and.iet(1).gt.0).or.
     &         (iet(1).eq.0.and.iep(1).gt.0))anintsdif=anintsdif+1.
          endif
          if(abs(itpr(1)).eq.2)then
            if(iep(1).eq.0.and.iet(1).gt.0)typevt=4 !SD tar
            if(iet(1).eq.0.and.iep(1).gt.0)typevt=3 !SD pro
            if(iep(1).eq.0.and.iet(1).eq.0)typevt=2 !low mass CD
            if(iep(1).gt.0.and.iet(1).gt.0)then
              typevt=1                              !DD
              anintddif=anintddif+1.
            endif
          else
            if(iep(1).eq.0.and.iet(1).eq.0)then
              typevt=-2          !DPE (high mass CD)    
            else
              typevt=-1         !ND
              if(iep(1).eq.0.and.iet(1).gt.0)typevt=-4 !SD tar high mass
              if(iet(1).eq.0.and.iep(1).gt.0)typevt=-3 !SD pro high mass
            endif
          endif
        endif
      else
        aidif=0.
        aidifp=0.
        aidift=0.
        aiine=0.
        do k=1,koll
          ip=iproj(k)
          it=itarg(k)
          if(aidif.ge.0..and.abs(itpr(k)).eq.2)then
            aidifp=aidifp+iep(ip)
            aidift=aidift+iet(it)
            aidif=aidif+1.
          elseif(itpr(k).ne.0)then
            aiine=aiine+iep(ip)+iet(it)
            aidifp=aidifp+iep(ip)
            aidift=aidift+iet(it)
            aidif=-ainfin
          endif
        enddo
        if(ionudi.eq.2)then
          aidif=aidif+aidifp
        else
          aidif=aidif+aidifp+aidift
        endif
        if(aidif.gt.0.)then
          anintdiff=anintdiff+1.
          anintine=anintine+1.
          if(aidifp.gt.0.5.and.aidift.le.0.5)then
            anintsdif=anintsdif+1.
            typevt=3                        !SD pro
          endif
          if(aidifp.gt.0.5.and.aidift.gt.0.5)then
            anintddif=anintddif+1.
            typevt=1                        !DD
          endif
          if(ionudi.ne.2)then
            if(aidifp.le.0.5.and.aidift.gt.0.5)then
              anintsdif=anintsdif+1.
              typevt=4                      !SD tar
            elseif(aidifp.le.0.5.and.aidift.le.0.5)then
              typevt=2                      !CD
            endif
          endif
        elseif(aiine.gt.0.)then
          anintine=anintine+1.
          typevt=-1                          !ND
          if(aidifp.gt.0.5.and.aidift.le.0.5)then
            anintsdif=anintsdif+1.
            typevt=-3                        !SD pro
          endif
          if(ionudi.ne.2)then
            if(aidifp.le.0.5.and.aidift.gt.0.5)then
              anintsdif=anintsdif+1.
              typevt=-4                      !SD tar
            endif
          endif
        else
          anintine=anintine+1.
          typevt=-2                         !DPE
        endif
      endif

      !if(nint(typevt).eq.3)print*,'++++++diffractive event+++++'
      if(ish.ge.6)call XPrint('After determining event typ&')
      end

      
c---------------------------------------------------------------------
      subroutine StringEnd(k,n,iret,ierr) !ierr.ne.0 redo event
c---------------------------------------------------------------------
c Not called for backup Pomeron
c Set String end type, pt and x+, xm- for 1 (Reggeon) or 2 (Pomeron) 
c strings
c---------------------------------------------------------------------
#include "aaa.h"
#include "ems.h"
      double precision pomass,px,py
      logical vpom,vreg
      dimension icp(2),ict(2),jcp(nflav,2),jct(nflav,2)
      double precision plc,s
      common/cems5/plc,s
      
      ierr=0

      if(idpr(n,k).eq.2)return    !reggeon with remnant mass only for the moment ???????????????
      iret=0
      if(idp1pr(n,k).ne.0)return    !SE already defined (second pass)
      ip=iproj(k)
      it=itarg(k)
    
      ! --- Set string end type and pt

      ntry=0
      vpom=.false.
      vreg=.false.
      ivpi=ivp(ip)
      ivti=ivt(it)
      idpi=idp(ip)
      idti=idt(it)
      do i=1,2
        icp(i)=icproj(i,ip)
        ict(i)=ictarg(i,it)
      enddo
      if(iremn.ge.2)then    !save jcpref and jctref into jcp and jct
        call UpdateFlav(ip,jcp,1)
        call UpdateFlav(it,jct,2)
      endif

 100  ntry=ntry+1
      iret=0
      if(ntry.ge.100)vpom=.true.
      if(ntry.ge.200)vreg=.true.
      if(ntry.gt.1)then
        if(ish.ge.4)write(ifch,*)
     &  'Try again setting string ends for k,n'
     &                           ,k,n,ntry
        ivp(ip)=ivpi
        ivt(it)=ivti
        idp(ip)=idpi
        idt(it)=idti
        do i=1,2
          icproj(i,ip)=icp(i)
          ictarg(i,it)=ict(i)
        enddo
        if(iremn.ge.2)then       !restore jcpref and jctref from jcp and jct
          call UpdateFlav(ip,jcp,-1)
          call UpdateFlav(it,jct,-2)
        endif
        call RmPt(k,n)
      endif

      call ProSeTy(k,n,ierr)
      if(ierr.ne.0)return !redo event 
      call ProSePt(k,n)


      ! --- Check Soft string Pomeron mass first

      px=xxp1pr(n,k)+xxp2pr(n,k)+xxm1pr(n,k)+xxm2pr(n,k)
      py=xyp1pr(n,k)+xyp2pr(n,k)+xym1pr(n,k)+xym2pr(n,k)
      pomass=xpr(n,k)*s-px*px-py*py
      q2kk=max(q2kmin(1,n,k),q2kmin(2,n,k))
      amprmn(1)=dsqrt(4d0*dble(q2kk))
      amprmn(2)=amprmn(1)
      amprmn(3)=amprmn(1)
      if(pomass.lt.amprmn(max(0,idhpr(n,k))))then
        if(vpom.and.idpr(n,k).ne.6)then
          nnb=nbkpr(n,k)
          if(nnb.ne.0)then      !big Pomeron with bckp one
            ivi=3
            call RmPt(k,n) !restore pt which will be define later for the backup
            call VirPom(k,n,ivi)
c           replace hard by soft
            npr(1,k)=npr(1,k)+1
            npr(3,k)=npr(3,k)-1
            idp1pr(n,k)=0
            iret=1
            return          !backup Pomeron treated later
          else
c           not enough energy for a Pomeron : try a resonance or a single string
            idpr(n,k)=6
            goto 100
          endif
        elseif(vreg)then   !even minimal conf doesn't work, suppress CD
          ivi=4
          call VirPom(k,n,ivi)
          if(ivi.eq.-1)then
            iret=2
            return
          endif
        else
          goto 100
        endif
      endif


      iret=0

      ! --- Define String ends for "normal" Pomerons ---

      call ProSeX(k,n,iret)
      if(iret.eq.1)then
        if(vpom.and.idpr(n,k).ne.6)then
          nnb=nbkpr(n,k)
          if(nnb.ne.0)then      !big Pomeron with bckp one
            ivi=3
            call RmPt(k,n) !restore pt which will be define later for the backup
            call VirPom(k,n,ivi)
c           replace hard by soft
            npr(1,k)=npr(1,k)+1
            npr(3,k)=npr(3,k)-1
            idp1pr(n,k)=0
            iret=1
            return              !backup Pomeron treated later
          else
c           not enough energy for a Pomeron : try a resonance or a single string
            idpr(n,k)=6
            goto 100
          endif
        elseif(vreg)then   !even minimal conf doesn't work, suppress Reggeon
          ivi=4
          call VirPom(k,n,ivi)
          if(ivi.eq.-1)then
            iret=2
            return
          endif
        else
          goto 100
        endif
      endif
      iret=0
      
      end

c-----------------------------------------------------------------------
      subroutine ProSeTy(k,n,ierr) !ierr.ne.0 redo event
c-----------------------------------------------------------------------
c creates proposal for string ends, idp., idm.
c updates quark counters
c-----------------------------------------------------------------------
#include "aaa.h"
#include "ems.h"
#include "sem.h"

      common/ems6/ivp0,iap0,idp0,isp0,ivt0,iat0,idt0,ist0
      double precision pes,xfqp,xfqt   !so01
      parameter(eps=1.e-6)
      common/ems9/xfqp(0:9),xfqt(0:9)
      common/emsx3/pes(0:3,0:6)
      integer jcp(nflav,2),jct(nflav,2)
     &       ,jcpi(nflavems,2),jcti(nflavems,2)
      logical go

      ierr=0

      if(idpr(n,k).gt.3.and.idpr(n,k).ne.6)then
        write(ifmt,'(a)')'WARNING lea ProSeTy idpr>3  redo event'
        ierr=1
        return
      endif

      iret=0
      ip=iproj(k)
      it=itarg(k)
      if(iremn.ge.3)then
        do j=1,2
          do i=1,nrflav
            jcp(i,j)=jcpref(i,j,ip)
            jct(i,j)=jctref(i,j,it)
          enddo
          do i=nrflav+1,nflav
            jcp(i,j)=0
            jct(i,j)=0
          enddo
        enddo
      endif
      
      pssp = 0.
      pvsp = 0.
      pvap = 0.
      psdp = 0.
      pdsp = 0.
      pddp = 0.
      psvvp= 0.
      paasp= 0.
      psst = 0.
      pvst = 0.
      pvat = 0.
      psdt = 0.
      pdst = 0.
      pddt = 0.
      psvvt= 0.
      paast= 0.
      paas=0.
      psvv=0.

      if(idpr(n,k).eq.6)then      !low mass CD

c use simplest configuration
      idp1pr(n,k)=1
      idm1pr(n,k)=0
      idp2pr(n,k)=0
      idm2pr(n,k)=1
      idsppr(n,k)=-1
      idstpr(n,k)=-1
      if(ish.ge.6)write(ifch,'(a3,$)')'CD'

      else

      idp1pr(n,k)=0
      idm1pr(n,k)=0
      idp2pr(n,k)=0
      idm2pr(n,k)=0
      idsppr(n,k)=0
      idstpr(n,k)=0


c for hard Pomeron, define which string ends are connected to valence quark
c treat gluon has soft string ends (including diquarks but can not be
c a "soft" valence like in soft Pomerons) later
      if(idpr(n,k).eq.3)then
        go=.false.
        if(ivp0.eq.iap0.and.rangen().lt.0.5)go=.true.    !meson
        idsppr(n,k)=5
        if(idhpr(n,k).eq.3.or.idhpr(n,k).eq.1)then
          if(iap0.eq.0.or.go)then !baryon
            idp1pr(n,k)=2
          else                    !antibaryon
            idp2pr(n,k)=2
          endif
        endif
        go=.false.
        if(ivt0.eq.iat0.and.rangen().lt.0.5)go=.true.    !meson
        idstpr(n,k)=5
        if(idhpr(n,k).eq.3.or.idhpr(n,k).eq.2)then
          if(iat0.eq.0.or.go)then     !baryon
            idm1pr(n,k)=2
          else                  !antibaryon
            idm2pr(n,k)=2
          endif
        endif
      endif

      if(idpr(n,k).ne.0)then

c    projectile

       if(idfpr(n,k).eq.1.or.idfpr(n,k).eq.2)then

       ntry=0
       ivpi=ivp(ip)
       idpi=idp(ip)
       idspi=idsppr(n,k)
       if(iremn.eq.3)then
         do j=1,2
           do i=1,nrflav
             jcpi(i,j)=jcp(i,j)
           enddo
         enddo
       endif
  1    ntry=ntry+1
      if(ntry.gt.10)call utstop('something goes wrong in sr ProSeTy&')
       ivp(ip)=ivpi
       idp(ip)=idpi
       idsppr(n,k)=idspi
       if(iremn.eq.3)then
         do j=1,2
           do i=1,nrflav
             jcp(i,j)=jcpi(i,j)
           enddo
         enddo
       endif
       pss=wgtval+wgtsea
       if(pss.gt.0.)then
         pss=wgtsea/pss
       else
         pss=0.
       endif
       if(iremn.ge.2)then
         if(iap0.eq.0)then
           pvs=0.
           if(ivp(ip).ne.0.and.idpr(n,k).ne.3)pvs=1.-pss
           pva=0.
           psvv=0.
           if(idp(ip).ne.0.and.ivp(ip).ge.2.and.idp2pr(n,k).ne.2)
     &     psvv=wgtqqq(iclpro)
           paas=0.
         elseif(ivp0.eq.0)then
           pva=0.
           if(ivp(ip).ne.0.and.idpr(n,k).ne.3)pva=1.-pss
           pvs=0.
           psvv=0.
           paas=0.
           if(idp(ip).ne.0.and.ivp(ip).ge.2.and.idp1pr(n,k).ne.2)
     &     paas=wgtqqq(iclpro)
         else                   !for meson, no soft string with valence quark (we do not know whether the quark or the antiquark will be used by hard string)
           pvs=0.
           pva=0.
c diquark or antidiquark can be created once in meson remnant
           psvv=0.
           paas=0.
           if(1+idp(ip).ne.0.and.ivp(ip).ge.1)then
             if(idp2pr(n,k).ne.2)psvv=wgtqqq(iclpro)
             if(idp1pr(n,k).ne.2)paas=wgtqqq(iclpro) 
c             if(abs(psvv-paas).lt.1e-6)then
cc since we don't know which one is possible, total P is sum of the 2
c               psvv=0.5*psvv
c               paas=0.5*paas
c             endif
           endif
         endif
         pdd=wgtdiq!/(1.+float(abs(idp(ip))))
c not to have diquark in hard string ends since Pomeron connection is not really in remnant but after some non-pertubative stuff
         if(idpr(n,k).eq.3)then
           pdd=0.
           psvv=0.
           paas=0.
         endif
       elseif(iremn.ne.0)then
         pvs=0.
         pva=0.
         psvv=0.
         paas=0.
         if(idp2pr(n,k).ne.2)psvv=wgtqqq(iclpro)
         if(idp1pr(n,k).ne.2)paas=wgtqqq(iclpro)
         pdd=wgtdiq!/(1.+float(abs(idp(ip))))
       else
         pvs=0.
         pva=0.
         psvv=0.
         paas=0.
         pdd=wgtdiq!/(1.+float(abs(idp(ip))))
       endif
       if(idp1pr(n,k).eq.2)then  !with valence quark only 1 SE available
         psd=pdd
         pds=0.
         pdd=0.
       elseif(idp2pr(n,k).eq.2)then  !with valence antiquark only 1 SE available
         pds=pdd
         psd=0.
         pdd=0.
       else
         psd=pdd
         pds=pdd
         pdd=pdd**2
       endif
       su=1.-min(1.,pdd+psd+pds)            !diquark probability
       pss=(1.-min(1.,pvs+pva))*su        !no more valence quark: take from sea
       pvs=pvs*su
       pva=pva*su
       su=1.-min(1.,psvv+paas)      !stopping probability
       pss=pss*su
       pvs=pvs*su
       pva=pva*su
       psd=psd*su
       pds=pds*su
       pdd=pdd*su
       su=pss+pvs+pva+pdd+psd+pds+psvv+paas
       pssp = pss /su
       pvsp = pvs /su
       pvap = pva /su
       psdp = psd /su
       pdsp = pds /su
       pddp = pdd /su
       psvvp= psvv/su
       paasp= paas/su
       r=rangen()
       if(r.gt.(pssp+pvsp+pvap+psdp+pdsp+psvvp+paasp)
     &                               .and.pddp.gt.eps)then
        if(idp1pr(n,k).ne.2)idp1pr(n,k)=4
        if(idp2pr(n,k).ne.2)idp2pr(n,k)=4
        idsppr(n,k)=idsppr(n,k)+4
        if(iremn.eq.3)then   !add diquark flavor to jcpref for ProSeF later (sea quark)
          idum=idrafl(iclpro,jcp,1,'s',3,iret)
          idum=idrafl(iclpro,jcp,1,'d',3,iret)
          idum=idrafl(iclpro,jcp,1,'s',3,iret)
          idum=idrafl(iclpro,jcp,1,'d',3,iret)
        endif
      elseif(r.gt.(pssp+pvsp+pvap+psdp+psvvp+paasp).and.pdsp.gt.eps)then
        if(idp1pr(n,k).ne.2)idp1pr(n,k)=4
        if(idp2pr(n,k).ne.2)idp2pr(n,k)=1
        idsppr(n,k)=idsppr(n,k)+4
        if(iremn.eq.3)then   !add diquark flavor to jcpref for ProSeF later (sea quark)
          idum=idrafl(iclpro,jcp,1,'s',3,iret)
          idum=idrafl(iclpro,jcp,1,'d',3,iret)
        endif
       elseif(r.gt.(pssp+pvsp+pvap+psvvp+paasp).and.psdp.gt.eps)then
        if(idp1pr(n,k).ne.2)idp1pr(n,k)=1
        if(idp2pr(n,k).ne.2)idp2pr(n,k)=4
        idsppr(n,k)=idsppr(n,k)+4
        if(iremn.eq.3)then   !add diquark flavor to jcpref for ProSeF later (sea quark)
          idum=idrafl(iclpro,jcp,1,'s',3,iret)
          idum=idrafl(iclpro,jcp,1,'d',3,iret)
        endif
       elseif(r.gt.(pssp+pvsp+pvap+psvvp).and.paasp.gt.eps)then
        if(idp1pr(n,k).ne.2)idp1pr(n,k)=5
        if(idp2pr(n,k).ne.2)idp2pr(n,k)=1
        idsppr(n,k)=idsppr(n,k)+5
        if(iremn.ge.2)then
          idp(ip)=idp(ip)-1
          if(ivp0.eq.iap0)then
c for mesons, valence quark and sea diquark
            ivp(ip)=ivp(ip)-1
            if(idp1pr(n,k).ne.2)idp1pr(n,k)=4  
            if(idp2pr(n,k).ne.2)idp2pr(n,k)=2
          else
            ivp(ip)=ivp(ip)-2
          endif
        endif
        if(iremn.eq.3)idum=idrafl(iclpro,jcp,1,'s',3,iret) !add flavor to jcpref for ProSeF later (sea quark) (only a q-aq pair because we replace diquark by q-aq (baryon "decay" or "stopping")
       elseif(r.gt.(pssp+pvsp+pvap+pddp).and.psvvp.gt.eps)then
        if(idp1pr(n,k).ne.2)idp1pr(n,k)=1
        if(idp2pr(n,k).ne.2)idp2pr(n,k)=5
        idsppr(n,k)=idsppr(n,k)+5
        if(iremn.ge.2)then
          idp(ip)=idp(ip)-1
          if(ivp0.eq.iap0)then
c for mesons, valence quark and sea diquark
            ivp(ip)=ivp(ip)-1
            if(idp1pr(n,k).ne.2)idp1pr(n,k)=2
            if(idp2pr(n,k).ne.2)idp2pr(n,k)=4
          else
            ivp(ip)=ivp(ip)-2
          endif
        endif
        if(iremn.eq.3)idum=idrafl(iclpro,jcp,1,'s',3,iret) !add flavor to jcpref for ProSeF later (sea quark) (only a q-aq pair because we replace diquark by q-aq (baryon "decay" or "stopping")
       elseif(r.gt.(pssp+pvsp).and.pvap.gt.eps)then
        if(idp1pr(n,k).ne.2)idp1pr(n,k)=1
        if(idp2pr(n,k).ne.2)idp2pr(n,k)=2
        idsppr(n,k)=idsppr(n,k)+2
        if(iremn.ge.2)ivp(ip)=ivp(ip)-1
        if(iremn.eq.3)idum=idrafl(iclpro,jcp,1,'s',3,iret) !add flavor to jcpref for ProSeF later (sea quark)
       elseif(r.gt.pssp.and.pvsp.gt.eps)then
        if(idp1pr(n,k).ne.2)idp1pr(n,k)=2
        if(idp2pr(n,k).ne.2)idp2pr(n,k)=1
        idsppr(n,k)=idsppr(n,k)+2
        if(iremn.ge.2)ivp(ip)=ivp(ip)-1
        if(iremn.eq.3)idum=idrafl(iclpro,jcp,2,'s',3,iret) !add flavor to jcpref for ProSeF later (sea quark)
       elseif(pssp.gt.eps)then
        if(idp1pr(n,k).ne.2)idp1pr(n,k)=1
        if(idp2pr(n,k).ne.2)idp2pr(n,k)=1
        idsppr(n,k)=idsppr(n,k)+1
        if(iremn.eq.3)idum=idrafl(iclpro,jcp,1,'s',3,iret) !add flavor to jcpref for ProSeF later (sea quark)
       else
        goto 1
       endif

       else
        idp1pr(n,k)=1
        idp2pr(n,k)=1
        idsppr(n,k)=0
       endif


c    target

       if(idfpr(n,k).eq.1.or.idfpr(n,k).eq.3)then


       ntry=0
       ivti=ivt(it)
       idti=idt(it)
       idsti=idstpr(n,k)
       if(iremn.eq.3)then
         do j=1,2
           do i=1,nrflav
             jcti(i,j)=jct(i,j)
           enddo
         enddo
       endif
  2    ntry=ntry+1
       if(ntry.gt.10)call utstop('something goes wrong in sr ProSeTy&')
       ivt(it)=ivti
       idt(it)=idti
       idstpr(n,k)=idsti
       if(iremn.eq.3)then
         do j=1,2
           do i=1,nrflav
             jct(i,j)=jcti(i,j)
           enddo
         enddo
       endif
       pss=wgtval+wgtsea
       if(pss.gt.0.)then
         pss=wgtsea/pss
       else
         pss=0.
       endif
       if(iremn.ge.2)then
         if(iat0.eq.0)then
           pvs=0.
           if(ivt(it).ne.0.and.idpr(n,k).ne.3)pvs=1.-pss
           pva=0.
           psvv=0.
           if(idt(it).ne.0.and.ivt(it).ge.2.and.idm2pr(n,k).ne.2)
     &     psvv=wgtqqq(icltar)
           paas=0.
         elseif(ivt0.eq.0)then
           pva=0.
           if(ivt(it).ne.0.and.idpr(n,k).ne.3)pva=1.-pss
           pvs=0.
           psvv=0.
           paas=0.
           if(idt(it).ne.0.and.ivt(it).ge.2.and.idm1pr(n,k).ne.2)
     &     paas=wgtqqq(icltar)
         else                   !for meson, no soft string with valence quark (we do not know whether the quark or the antiquark will be used by hard string)
           pvs=0.
           pva=0.
c diquark or antidiquark can be created once in meson remnant
           psvv=0.
           paas=0.
           if(1+idt(it).ne.0.and.ivt(it).ge.1)then
             if(idm2pr(n,k).ne.2)psvv=wgtqqq(icltar)
             if(idm1pr(n,k).ne.2)paas=wgtqqq(icltar)
c             if(abs(psvv-paas).lt.1e-6)then
cc since we don't know which one is possible, total P is sum of the 2
c               psvv=0.5*psvv
c               paas=0.5*paas
c             endif
           endif
         endif
         pdd=wgtdiq!/(1.+float(abs(idt(it))))
c not to have diquarks in hard string ends since Pomeron connection is not really in remnant but after some non-pertubative stuff
         if(idpr(n,k).eq.3)then
           pdd=0.
           psvv=0.
           paas=0.
         endif
       elseif(iremn.ne.0)then
         pvs=0.
         pva=0.
         psvv=wgtqqq(icltar)
         if(idm1pr(n,k).ne.2)paas=wgtqqq(icltar)
         pdd=wgtdiq!/(1.+float(abs(idt(it))))
       else
         pvs=0.
         pva=0.
         psvv=0.
         paas=0.
         pdd=wgtdiq!/(1.+float(abs(idt(it))))
       endif
       if(idm1pr(n,k).eq.2)then  !with valence quark only 1 SE available
         psd=pdd
         pds=0.
         pdd=0.
       elseif(idm2pr(n,k).eq.2)then  !with valence antiquark only 1 SE available
         pds=pdd
         psd=0.
         pdd=0.
       else
         psd=pdd
         pds=pdd
         pdd=pdd**2
       endif
       su=1.-min(1.,pdd+pds+psd)            !diquark probability
       pss=(1.-min(1.,pvs+pva))*su        !no more valence quark: take from sea
       pvs=pvs*su
       pva=pva*su
       su=1.-min(1.,psvv+paas)      !stopping probability
       pss=pss*su
       pvs=pvs*su
       pva=pva*su
       pds=pds*su
       psd=psd*su
       pdd=pdd*su
       su=pss+pvs+pva+pdd+psd+pds+psvv+paas
       psst = pss /su
       pvst = pvs /su
       pvat = pva /su
       psdt = psd /su
       pdst = pds /su
       pddt = pdd /su
       psvvt= psvv/su
       paast= paas/su
       r=rangen()
       if(r.gt.(psst+pvst+pvat+psdt+pdst+psvvt+paast)
     &                               .and.pddt.gt.eps)then
        if(idm1pr(n,k).ne.2)idm1pr(n,k)=4
        if(idm2pr(n,k).ne.2)idm2pr(n,k)=4
        idstpr(n,k)=idstpr(n,k)+4
        if(iremn.eq.3)then   !add diquark flavor to jctref for ProSeF later (sea quark)
          idum=idrafl(icltar,jct,1,'s',3,iret)
          idum=idrafl(icltar,jct,1,'d',3,iret)
          idum=idrafl(icltar,jct,1,'s',3,iret)
          idum=idrafl(icltar,jct,1,'d',3,iret)
        endif
      elseif(r.gt.(psst+pvst+pvat+psdt+psvvt+paast).and.pdst.gt.eps)then
        if(idm1pr(n,k).ne.2)idm1pr(n,k)=4
        if(idm2pr(n,k).ne.2)idm2pr(n,k)=1
        idstpr(n,k)=idstpr(n,k)+4
        if(iremn.eq.3)then   !add diquark flavor to jctref for ProSeF later (sea quark)
          idum=idrafl(icltar,jct,1,'s',3,iret)
          idum=idrafl(icltar,jct,1,'d',3,iret)
        endif
       elseif(r.gt.(psst+pvst+pvat+psvvt+paast).and.psdt.gt.eps)then
        if(idm1pr(n,k).ne.2)idm1pr(n,k)=1
        if(idm2pr(n,k).ne.2)idm2pr(n,k)=4
        idstpr(n,k)=idstpr(n,k)+4
        if(iremn.eq.3)then   !add diquark flavor to jctref for ProSeF later (sea quark)
          idum=idrafl(icltar,jct,1,'s',3,iret)
          idum=idrafl(icltar,jct,1,'d',3,iret)
        endif
       elseif(r.gt.(psst+pvst+pvat+psvvt).and.paast.gt.eps)then
        if(idm1pr(n,k).ne.2)idm1pr(n,k)=5
        if(idm2pr(n,k).ne.2)idm2pr(n,k)=1
        idstpr(n,k)=idstpr(n,k)+5
        if(iremn.ge.2)then
          idt(it)=idt(it)-1
          if(ivt0.eq.iat0)then
c for mesons, valence quark and sea diquark
            ivt(it)=ivt(it)-1
            if(idm1pr(n,k).ne.2)idm1pr(n,k)=4  
            if(idm2pr(n,k).ne.2)idm2pr(n,k)=2
          else
            ivt(it)=ivt(it)-2
          endif
        endif
        if(iremn.eq.3)idum=idrafl(icltar,jct,1,'s',3,iret) !add flavor to jcpref for ProSeF later (sea quark) (only a q-aq pair because we replace diquark by q-aq (baryon "decay" or "stopping")
       elseif(r.gt.(psst+pvst+pvat+pddt).and.psvvt.gt.eps)then
        if(idm1pr(n,k).ne.2)idm1pr(n,k)=1
        if(idm2pr(n,k).ne.2)idm2pr(n,k)=5
        idstpr(n,k)=idstpr(n,k)+5
        if(iremn.ge.2)then
          idt(it)=idt(it)-1
          if(ivt0.eq.iat0)then
c for mesons, valence quark and sea diquark
            ivt(it)=ivt(it)-1
            if(idm1pr(n,k).ne.2)idm1pr(n,k)=2 
            if(idm2pr(n,k).ne.2)idm2pr(n,k)=4
          else
            ivt(it)=ivt(it)-2
          endif
        endif
        if(iremn.eq.3)idum=idrafl(icltar,jct,1,'s',3,iret) !add flavor to jcpref for ProSeF later (sea quark) (only a q-aq pair because we replace diquark by q-aq (baryon "decay" or "stopping")
       elseif(r.gt.(psst+pvst).and.pvat.gt.eps)then
        if(idm1pr(n,k).ne.2)idm1pr(n,k)=1
        if(idm2pr(n,k).ne.2)idm2pr(n,k)=2
        idstpr(n,k)=idstpr(n,k)+2
        if(iremn.ge.2)ivt(it)=ivt(it)-1
        if(iremn.eq.3)idum=idrafl(icltar,jct,1,'s',3,iret) !add flavor to jctref for ProSeF later (sea quark)
       elseif(r.gt.psst.and.pvst.gt.eps)then
        if(idm1pr(n,k).ne.2)idm1pr(n,k)=2
        if(idm2pr(n,k).ne.2)idm2pr(n,k)=1
        idstpr(n,k)=idstpr(n,k)+2
        if(iremn.ge.2)ivt(it)=ivt(it)-1
        if(iremn.eq.3)idum=idrafl(icltar,jct,2,'s',3,iret) !add flavor to jctref for ProSeF later (sea quark)
       elseif(psst.gt.eps)then
        if(idm1pr(n,k).ne.2)idm1pr(n,k)=1
        if(idm2pr(n,k).ne.2)idm2pr(n,k)=1
        idstpr(n,k)=idstpr(n,k)+1
        if(iremn.eq.3)idum=idrafl(icltar,jct,1,'s',3,iret) !add flavor to jctref for ProSeF later (sea quark)
       else
        goto 2
       endif

       else
        idm1pr(n,k)=1
        idm2pr(n,k)=1
        idstpr(n,k)=0
       endif

      else

        idp1pr(n,k)=0
        idm2pr(n,k)=0
        idp2pr(n,k)=0
        idm1pr(n,k)=0

      endif

        if(ish.ge.6)then
      write(ifch,'(a,8(f4.2,1x))')'ProSeTy Projectile  ',
     * pssp,pvsp,pvap,pddp,psdp,pdsp,psvvp,paasp
      write(ifch,'(a,8(f4.2,1x))')'ProSeTy Target  ',
     * psst,pvst,pvat,pddt,psdt,pdst,psvvt,paast
      write(ifch,'(a3,$)')'Pom'
        endif

      endif                     !end CD

        if(ish.ge.6) 
     * write(ifch,'(2x,3i3,2x,2(a1,1x,i2,1x,2i2,1x,i2,i3,2x))')
     *  idpr(n,k),n,k
     * ,'p',idsppr(n,k),idp1pr(n,k),idp2pr(n,k),ivp(ip),idp(ip)
     * ,'t',idstpr(n,k),idm1pr(n,k),idm2pr(n,k),ivt(it),idt(it)

      if(iremn.eq.3)then
        do j=1,2
          do i=1,nrflav
            jcpref(i,j,ip)=jcp(i,j)
            jctref(i,j,it)=jct(i,j)
          enddo
        enddo
        if(ish.ge.6)then
          write(ifch,'(a,i3,a,1x,4i3,3x,4i3)')'jcpref(',ip,'):',jcp
          write(ifch,'(a,i3,a,1x,4i3,3x,4i3)')'jctref(',it,'):',jct
        endif
      endif

      return
      end

c-----------------------------------------------------------------------
      subroutine fstrfl(jcp,jct,jcpv,jctv,icp1,icp2,icm1,icm2
     *                         ,idp1,idp2,idm1,idm2,idsp,idst,iret)
c-----------------------------------------------------------------------
c knowing the string end types (idp1,idp2,idm1,idm2)
c               and remnant flavors (icp,ict)
c               and remnant link of the string (idsp and idst)
c               (idsp/t=100 with one idp/idm=2 means that the valence quark 
c                to use is define in the corresponding icp/icm
c                (using just 1 to 6 for flavor identification (no diquark)))
c one determines quark flavors of string ends (icp1,icp2,icm1,icm2)
c               and updates remnant flavors (icp,ict)
c iret=0   ok
c iret=1   problem, more than 9 quarks per flavor attempted
c-----------------------------------------------------------------------
#include "aaa.h"
#include "ems.h"
#include "sem.h"
      integer icp1(2),icp2(2),icm1(2),icm2(2)
      integer jcp(nflav,2),jct(nflav,2)
     &       ,jcpi(nflavems,2),jcti(nflavems,2)
      integer iq(2,4),jcpv(nflav,2),jctv(nflav,2)
      character m
c      data neuz/0/proz/0/dtaz/0/
c      save neuz,proz,dtaz

      call utpri('fstrfl',ish,ishini,7)

c     entry
c     -----

      idum=0
      iret=0
      iret1=0
      iret2=0
      iret3=0
      iret4=0

      if(idp1.eq.8)stop'fstrfl: fragm quarks not used any more'
      if(idp2.eq.8)stop'fstrfl: fragm quarks not used any more'
      if(idm1.eq.8)stop'fstrfl: fragm quarks not used any more'
      if(idm2.eq.8)stop'fstrfl: fragm quarks not used any more'

c determine flavors of string ends (u,d,s)

      if(ish.ge.7)then
       write(ifch,'(a,3x,2i3)')' string 1, SE types:',idp1,idm2
       write(ifch,'(a,3x,2i3)')' string 2, SE types:',idp2,idm1
      endif

c empty

      if(idp1.eq.0)then
       iq(1,1)=0
       iq(2,1)=0
      endif
      if(idp2.eq.0)then
       iq(1,2)=0
       iq(2,2)=0
      endif
      if(idm1.eq.0)then
       iq(1,4)=0
       iq(2,4)=0
      endif
      if(idm2.eq.0)then
       iq(1,3)=0
       iq(2,3)=0
      endif
      do j=1,2
        do n=1,nrflav
          jcpi(n,j)=jcp(n,j)
          jcti(n,j)=jct(n,j)
        enddo
      enddo

c Projectile

      if(idsp.le.0.or.iremn.eq.0)then
c give the same flavor to quark and antiquark not to change remnant flavor

        if(idp1.eq.4)then
c diquarks, code 4
          iq(1,1)=idrafl(iclpro,jcp,1,'s',0,iret)
          iq(2,1)=idrafl(iclpro,jcp,1,'d',0,iret)
          iq(1,2)=iq(1,1)
          iq(2,2)=iq(2,1)
        elseif(idsp.ge.0)then
c sea quarks, code 1
          iq(1,1)=idrafl(iclpro,jcp,1,'s',0,iret)
          iq(2,1)=0
          iq(1,2)=iq(1,1)
          iq(2,2)=0
        else            !CD
c sea quarks, code 1
          iq(1,1)=idrafl(iclpro,jcp,1,'s',0,iret)
          iq(2,1)=0
          iq(1,2)=0
          iq(2,2)=0
        endif

      elseif(iremn.ge.2)then
c count valence quarks properly

c valence quarks

        if(idp1.eq.2)then
          if(idsp.eq.100)then
            iq(1,1)=icp1(1)         !flavor of hard quark already defined
          else
            iq(1,1)=idrafl(iclpro,jcpv,1,'v',0,idum)
          endif
          if(iq(1,1).gt.0)then          !if still exist, update jcp and jcpv
            call idsufl3(iq(1,1),1,jcpv)
c            call idsufl3(iq(1,1),1,jcp)
          else                          ! if not, use jcp directly and sea
            iq(1,1)=idrafl(iclpro,jcp,1,'s',1,idum)
          endif
          iq(2,1)=0
        endif
        if(idp2.eq.2)then
          if(idsp.eq.100)then
            iq(1,2)=icp2(2)         !flavor of hard antiquark already defined
          else
            iq(1,2)=idrafl(iclpro,jcpv,2,'v',0,idum)
          endif
          if(iq(1,2).gt.0)then          !if still exist, update jcp and jcpv
            call idsufl3(iq(1,2),2,jcpv)
c            call idsufl3(iq(1,2),2,jcp)
          else                          ! if not, use jcp directly and sea
            iq(1,2)=idrafl(iclpro,jcp,2,'s',1,idum)
          endif
          iq(2,2)=0
        endif

c sea quarks
        m='v'           !iremn=3

        if(idp1.eq.1)then
          if(iremn.eq.2)m='s'
          j=1                   !quark
          i=idrafl(iclpro,jcp,j,m,1,idum)
          iq(1,1)=i
c          if(jcp(i,j)-jcpv(i,j).lt.0)jcpv(i,j)=jcpv(i,j)-1
          iq(2,1)=0
        elseif(idp1.ge.4)then
          j=2                    !anti-diquark
          if(iremn.eq.2.and.idp1.eq.5)then
c valence diquark transfer
            i=idrafl(iclpro,jcpv,j,m,1,idum)
            iq(1,1)=i
            i=idrafl(iclpro,jcpv,j,m,1,idum)
            iq(2,1)=i
          else
            if(iremn.eq.2)m='s'
            i=idrafl(iclpro,jcp,j,m,1,idum)
            iq(1,1)=i
            m='v'               !iremn=3
            if(iremn.eq.2)m='d'
            i=idrafl(iclpro,jcp,j,m,1,idum)
            iq(2,1)=i
          endif
c          if(jcp(i,j)-jcpv(i,j).lt.0)jcpv(i,j)=jcpv(i,j)-1
c          if(jcp(i,j)-jcpv(i,j).lt.0)jcpv(i,j)=jcpv(i,j)-1
        endif
        m='v'           !iremn=3
        if(idp2.eq.1)then
          if(iremn.eq.2)m='s'
          j=2                    !antiquark
          i=idrafl(iclpro,jcp,j,m,1,idum)
          iq(1,2)=i
c          if(jcp(i,j)-jcpv(i,j).lt.0)jcpv(i,j)=jcpv(i,j)-1
          iq(2,2)=0
        elseif(idp2.ge.4)then
          j=1                    !diquark
          if(iremn.eq.2.and.idp2.eq.5)then
c valence diquark transfer
            i=idrafl(iclpro,jcpv,j,m,1,idum)
            iq(1,2)=i
            i=idrafl(iclpro,jcpv,j,m,1,idum)
            iq(2,2)=i
          else
            if(iremn.eq.2)m='s'
            i=idrafl(iclpro,jcp,j,m,1,idum)
            iq(1,2)=i
            m='v'               !iremn=3
            if(iremn.eq.2)m='d'
            i=idrafl(iclpro,jcp,j,m,1,idum)
            iq(2,2)=i
          endif
c          if(jcp(i,j)-jcpv(i,j).lt.0)jcpv(i,j)=jcpv(i,j)-1
c          if(jcp(i,j)-jcpv(i,j).lt.0)jcpv(i,j)=jcpv(i,j)-1
        endif


      elseif(iremn.ne.0)then
c free remnant content

c valence quarks

        if(idp1.eq.2)then
          if(idsp.eq.100)then
            iq(1,1)=icp1(1)         !flavor of hard quark already defined
          else
            iq(1,1)=idrafl(iclpro,jcp,1,'v',1,iret)
          endif
          iq(2,1)=0
        endif
        if(idp2.eq.2)then
          if(idsp.eq.100)then
            iq(1,2)=icp2(1)     !flavor of hard antiquark already defined
          else
            iq(1,2)=idrafl(iclpro,jcp,2,'v',1,iret)
          endif
          iq(2,2)=0
        endif

c sea quarks

        if(idp1.eq.1)then
          iq(1,1)=idrafl(iclpro,jcp,1,'s',1,iret1)
          iq(2,1)=0
        endif
        if(idp2.eq.1)then
          iq(1,2)=idrafl(iclpro,jcp,2,'s',1,iret2)
          iq(2,2)=0
        endif

c diquarks, code 4

        if(idp1.eq.4.or.idp2.eq.4)then
          iq(1,1)=idrafl(iclpro,jcp,2,'s',1,iret1)
          iq(2,1)=idrafl(iclpro,jcp,2,'d',1,iret1)
          iq(1,2)=idrafl(iclpro,jcp,1,'s',1,iret2)
          iq(2,2)=idrafl(iclpro,jcp,1,'d',1,iret2)
        endif

c diquarks, code 5 (former valence, but actually sea)

        if(idp1.eq.5)then
          iq(1,1)=idrafl(iclpro,jcp,2,'s',1,iret1)
          iq(2,1)=idrafl(iclpro,jcp,2,'d',1,iret1)
        endif
        if(idp2.eq.5)then
          iq(1,2)=idrafl(iclpro,jcp,1,'s',1,iret2)
          iq(2,2)=idrafl(iclpro,jcp,1,'d',1,iret2)
        endif


        if(iret.ne.0)goto 1000



c in case of saturated remnants, use the same flavor for quark and anti-quark
c at string-end
        if(iret1.ne.0.or.iret2.ne.0)then
          do j=1,2
            do n=1,nrflav
              jcp(n,j)=jcpi(n,j)
            enddo
          enddo
          if(idp1.gt.idp2.or.(idp1.eq.idp2.and.rangen().gt.0.5))then
            iq(1,2)=iq(1,1)
            iq(2,2)=iq(2,1)
          else
            iq(1,1)=iq(1,2)
            iq(2,1)=iq(2,2)
          endif
        endif

      endif

c Target

      if(idst.le.0.or.iremn.eq.0)then
c give the same flavor to quark and antiquark not to change remnant flavor


        if(idm1.eq.4)then
c diquarks, code 4
          iq(1,4)=idrafl(icltar,jct,1,'s',0,iret)
          iq(2,4)=idrafl(icltar,jct,1,'d',0,iret)
          iq(1,3)=iq(1,4)
          iq(2,3)=iq(2,4)
        elseif(idst.ge.0)then
c sea quarks,code 1
          iq(1,4)=idrafl(icltar,jct,1,'s',0,iret)
          iq(2,4)=0
          iq(1,3)=iq(1,4)
          iq(2,3)=0
        else          !CD
c sea quarks,code 1
          iq(1,4)=0
          iq(2,4)=0
          iq(1,3)=iq(1,1)
          iq(2,3)=0
        endif

      elseif(iremn.ge.2)then
c count valence quarks properly

c valence quarks

        if(idm1.eq.2)then
          if(idst.eq.100)then
            iq(1,4)=icm1(1)         !flavor of hard quark already defined
          else
            iq(1,4)=idrafl(icltar,jctv,1,'v',0,idum)
          endif
          if(iq(1,4).gt.0)then          !if still exist, update jct and jctv
            call idsufl3(iq(1,4),1,jctv)
c            call idsufl3(iq(1,4),1,jct)
          else                          ! if not, use jct directly
            iq(1,4)=idrafl(icltar,jct,1,'s',1,idum)
          endif
          iq(2,4)=0
        endif
        if(idm2.eq.2)then
          if(idst.eq.100)then
            iq(1,3)=icm2(2)     !flavor of hard antiquark already defined
          else
            iq(1,3)=idrafl(icltar,jctv,2,'v',0,idum)
          endif
          if(iq(1,3).gt.0)then  !if still exist, update jct and jctv
            call idsufl3(iq(1,3),2,jctv)
c            call idsufl3(iq(1,3),2,jct)
          else                          ! if not, use jct directly
            iq(1,3)=idrafl(icltar,jct,2,'s',1,idum)
          endif
          iq(2,3)=0
        endif

c sea quarks
        m='v'           !iremn=3

        if(idm1.eq.1)then
          if(iremn.eq.2)m='s'
          j=1                    !quark
          i=idrafl(icltar,jct,j,m,1,idum)
          iq(1,4)=i
c          if(jct(i,j)-jctv(i,j).lt.0)jctv(i,j)=jctv(i,j)-1
          iq(2,4)=0
        elseif(idm1.ge.4)then
          j=2                   !anti-diquark
          if(iremn.eq.2.and.idm1.eq.5)then
c valence diquark transfer
            i=idrafl(icltar,jctv,j,m,1,idum)
            iq(1,4)=i
            i=idrafl(icltar,jctv,j,m,1,idum)
            iq(2,4)=i
          else
            if(iremn.eq.2)m='s'
            i=idrafl(icltar,jct,j,m,1,idum)
            iq(1,4)=i
            m='v'               !iremn=3
            if(iremn.eq.2)m='d'
            i=idrafl(icltar,jct,j,m,1,idum)
            iq(2,4)=i
          endif
c          if(jct(i,j)-jctv(i,j).lt.0)jctv(i,j)=jctv(i,j)-1
c          if(jct(i,j)-jctv(i,j).lt.0)jctv(i,j)=jctv(i,j)-1
        endif
        m='v'           !iremn=3
        if(idm2.eq.1)then
          if(iremn.eq.2)m='s'
          j=2                   !antiquark
          i=idrafl(icltar,jct,j,m,1,idum)
          iq(1,3)=i
c          if(jct(i,j)-jctv(i,j).lt.0)jctv(i,j)=jctv(i,j)-1
          iq(2,3)=0
        elseif(idm2.ge.4)then
          j=1                    !diquark
          if(iremn.eq.2.and.idm2.eq.5)then
c valence diquark transfer
            i=idrafl(icltar,jctv,j,m,1,idum)
            iq(1,3)=i
            i=idrafl(icltar,jctv,j,m,1,idum)
            iq(2,3)=i
          else
            if(iremn.eq.2)m='s'
            i=idrafl(icltar,jct,j,m,1,idum)
            iq(1,3)=i
            m='v'               !iremn=3
            if(iremn.eq.2)m='d'
            i=idrafl(icltar,jct,j,m,1,idum)
            iq(2,3)=i
          endif
c          if(jct(i,j)-jctv(i,j).lt.0)jctv(i,j)=jctv(i,j)-1
c          if(jct(i,j)-jctv(i,j).lt.0)jctv(i,j)=jctv(i,j)-1
        endif

      elseif(iremn.ne.0)then

c valence quarks

        if(idm1.eq.2)then
          if(idst.eq.100)then
            iq(1,4)=icm1(1)         !flavor of hard quark already defined
          else
            iq(1,4)=idrafl(icltar,jct,1,'v',1,iret)
          endif
          iq(2,4)=0
        endif
        if(idm2.eq.2)then
          if(idst.eq.100)then
            iq(1,3)=icm2(1)         !flavor of hard antiquark already defined
          else
            iq(1,3)=idrafl(icltar,jct,2,'v',1,iret)
          endif
          iq(2,3)=0
        endif

c sea quarks

        if(idm1.eq.1)then
          iq(1,4)=idrafl(icltar,jct,1,'s',1,iret4)
          iq(2,4)=0
        endif
        if(idm2.eq.1)then
          iq(1,3)=idrafl(icltar,jct,2,'s',1,iret3)
          iq(2,3)=0
        endif

c diquarks, code 4

        if(idm1.eq.4.or.idm2.eq.4)then
          iq(1,4)=idrafl(icltar,jct,2,'s',1,iret3)
          iq(2,4)=idrafl(icltar,jct,2,'d',1,iret3)
          iq(1,3)=idrafl(icltar,jct,1,'s',1,iret4)
          iq(2,3)=idrafl(icltar,jct,1,'d',1,iret4)
        endif

c diquarks, code 5 (former valence, but actually sea)

        if(idm1.eq.5)then
          iq(1,4)=idrafl(icltar,jct,2,'s',1,iret4)
          iq(2,4)=idrafl(icltar,jct,2,'d',1,iret4)
        endif
        if(idm2.eq.5)then
          iq(1,3)=idrafl(icltar,jct,1,'s',1,iret3)
          iq(2,3)=idrafl(icltar,jct,1,'d',1,iret3)
        endif


        if(iret.ne.0)goto 1000



c in case of saturated remnants, use the same flavor for quark and anti-quark
c at string-end

        if(iret3.ne.0.or.iret4.ne.0)then
          do j=1,2
            do n=1,nrflav
              jct(n,j)=jcti(n,j)
            enddo
          enddo
          if(idm1.gt.idm2.or.(idm1.eq.idm2.and.rangen().gt.0.5))then
            iq(1,4)=iq(1,3)
            iq(2,4)=iq(2,3)
          else
            iq(1,3)=iq(1,4)
            iq(2,3)=iq(2,4)
          endif
        endif

      endif

      ifla=iq(1,1)
      iflb=iq(2,1)
      iflc=iq(1,3)
      ifld=iq(2,3)
      if(ish.ge.7)write(ifch,'(a,2i5,4x,2i5)')
     *' string 1, string ends:',ifla,iflb,iflc,ifld

      if(ifla.gt.0)then
       if(iflb.eq.0)then
        icp1(1)=10**(6-ifla)
        icp1(2)=0
       else
        icp1(1)=0
        icp1(2)=10**(6-ifla)
        icp1(2)=icp1(2)+10**(6-iflb)
       endif
      endif

      if(iflc.gt.0)then
       if(ifld.eq.0)then
        icm2(1)=0
        icm2(2)=10**(6-iflc)
       else
        icm2(1)=10**(6-iflc)
        icm2(1)=icm2(1)+10**(6-ifld)
        icm2(2)=0
       endif
      endif

      ifla=iq(1,4)
      iflb=iq(2,4)
      iflc=iq(1,2)
      ifld=iq(2,2)
      if(ish.ge.7)write(ifch,'(a,2i5,4x,2i5)')
     *' string 2, string ends:',ifla,iflb,iflc,ifld

      if(ifla.gt.0)then
       if(iflb.eq.0)then
        icm1(1)=10**(6-ifla)
        icm1(2)=0
       else
        icm1(1)=0
        icm1(2)=10**(6-ifla)
        icm1(2)=icm1(2)+10**(6-iflb)
       endif
      endif

      if(iflc.gt.0)then
       if(ifld.eq.0)then
        icp2(1)=0
        icp2(2)=10**(6-iflc)
       else
        icp2(1)=10**(6-iflc)
        icp2(1)=icp2(1)+10**(6-ifld)
        icp2(2)=0
       endif
      endif

      if(ish.ge.7)then
        write(ifch,'(a,2i7,4x,2i7)')
     *  ' SE-forw:',icp1(1),icp1(2),icp2(1),icp2(2)
        write(ifch,'(a,2i7,4x,2i7)')
     *  ' SE-back:',icm2(1),icm2(2),icm1(1),icm1(2)
        write(ifch,'(a,2(3x,6i3,3x,6i3))')' proj:',jcp,jcpv
        write(ifch,'(a,2(3x,6i3,3x,6i3))')' targ:',jct,jctv
      endif

c     exit
c     ----

1000  continue
      call utprix('fstrfl',ish,ishini,7)
      return
      end


cc-----------------------------------------------------------------------
c      subroutine fremfl(icp,ict,iret)
cc-----------------------------------------------------------------------
cc checks projectile and target flavor (icp,ict)
cc in case of reggeon exchange they do not correspond to hadrons.
cc one transfers therefore flavor from one side to the other in order
cc to have hadron flavor.
cc icp and ict are modified correspondingly
cc-----------------------------------------------------------------------
c#include "aaa.h"
c      integer icp(2),ict(2),jcp(6,2),jct(6,2),kp(4),kt(4)
c
c      call utpri('fremfl',ish,ishini,7)
c
cc     entry
cc     -----
c
c      iret=0
c
c      call iddeco(icp,jcp)
c      call iddeco(ict,jct)
c
c      iakp=0
c      iakt=0
c      ikp=0
c      ikt=0
c      do l=1,4
c       kp(l)=jcp(l,1)-jcp(l,2)
c       kt(l)=jct(l,1)-jct(l,2)
c       iakp=iakp+iabs(kp(l))
c       iakt=iakt+iabs(kt(l))
c       ikp=ikp+kp(l)
c       ikt=ikt+kt(l)
c      enddo
c      if(ish.ge.7)write(ifch,*)'iak_p:',iakp,' ik_p:',ikp
c      if(ish.ge.7)write(ifch,*)'iak_t:',iakt,' ik_t:',ikt
c
c      if(iakp.eq.4)then
c       if(ikp.eq.4.or.ikp.eq.-2)then
c        ifl=idrafl(jcp,1,'v',iret)
c        iqp=2      ! subtract quark
c        iqt=1      ! add quark
c       elseif(ikp.eq.-4.or.ikp.eq.2)then
c        ifl=idrafl(jcp,2,'v',iret)
c        iqp=1      ! subtract antiquark
c        iqt=2      ! add antiquark
c       else
c        call utstop('fremfl&')
c       endif
c      elseif(iakt.eq.4)then
c       if(ikt.eq.4.or.ikt.eq.-2)then
c        ifl=idrafl(jct,1,'v',iret)
c        iqp=1      ! add quark
c        iqt=2      ! subtract quark
c       elseif(ikt.eq.-4.or.ikt.eq.2)then
c        ifl=idrafl(jct,2,'v',iret)
c        iqp=2      ! add antiquark
c        iqt=1      ! subtract antiquark
c       else
c        call utstop('fremfl&')
c       endif
c      elseif(iakp.eq.3)then
c       if(ikp.gt.0)then
c        ifl=idrafl(jcp,1,'v',iret)
c        iqp=2      ! subtract quark
c        iqt=1      ! add quark
c       else
c        ifl=idrafl(jcp,2,'v',iret)
c        iqp=1      ! subtract antiquark
c        iqt=2      ! add antiquark
c       endif
c      elseif(iakt.eq.3)then
c       if(ikt.gt.0)then
c        ifl=idrafl(jct,1,'v',iret)
c        iqp=1      ! add quark
c        iqt=2      ! subtract quark
c       else
c        ifl=idrafl(jct,2,'v',iret)
c        iqp=2      ! add antiquark
c        iqt=1      ! subtract antiquark
c       endif
c      elseif(iakp.eq.2)then
c       if(ikp.gt.0)then
c        ifl=idrafl(jct,1,'v',iret)
c        iqp=1      ! add quark
c        iqt=2      ! subtract quark
c       else
c        ifl=idrafl(jct,2,'v',iret)
c        iqp=2      ! add antiquark
c        iqt=1      ! subtract antiquark
c       endif
c      elseif(iakt.eq.2)then
c       if(ikt.gt.0)then
c        ifl=idrafl(jct,1,'v',iret)
c        iqp=2      ! subtract quark
c        iqt=1      ! add quark
c       else
c        ifl=idrafl(jct,2,'v',iret)
c        iqp=1      ! subtract antiquark
c        iqt=2      ! add antiquark
c       endif
c      elseif(iakp.eq.1)then
c       if(ikp.gt.0)then
c        ifl=idrafl(jcp,2,'v',iret)
c        iqp=2      ! add antiquark
c        iqt=1      ! subtract antiquark
c       else
c        ifl=idrafl(jcp,1,'v',iret)
c        iqp=1      ! add quark
c        iqt=2      ! subtract quark
c       endif
c      elseif(iakt.eq.1)then
c       if(ikt.gt.0)then
c        ifl=idrafl(jct,2,'v',iret)
c        iqp=1      ! subtract antiquark
c        iqt=2      ! add antiquark
c       else
c        ifl=idrafl(jct,1,'v',iret)
c        iqp=2      ! subtract quark
c        iqt=1      ! add quark
c       endif
c      else
c       call utstop('fremfl: error&')
c      endif
c
c      if(ish.ge.7)write(ifch,*)'iq_p:',iqp,' iq_t:',iqt,' if:',ifl
c      call uticpl(icp,ifl,iqp,iret)
c      if(iret.ne.0)goto1000
c      call uticpl(ict,ifl,iqt,iret)
c      if(iret.ne.0)goto1000
c
cc     exit
cc     ----
c
c1000  continue
c      call utprix('fremfl',ish,ishini,7)
c      return
c      end
c

c-----------------------------------------------------------------------
      subroutine UpdateFlav(ir,jc,mod)
c-----------------------------------------------------------------------
C Add valence quark to sea quarks in projectile jcpref (mod=10) or target
c jctref (mod=20) for soft string ends (mod=0 reset jcrpref and
c jctref to 0).
c For mod=1 or 2, save jcref into jc.
c For mod=-1 or -2, put jc into jcref.
c-----------------------------------------------------------------------
#include "aaa.h"
#include "ems.h"
#include "sem.h"
      dimension ic(2),jc(nflav,2),jc2(nflav,2)

      if(mod.eq.0)then
        do j=1,2
          do i=1,nrflav
            jcpref(i,j,ir)=0
            jctref(i,j,ir)=0
          enddo
        enddo
      elseif(mod.eq.-1)then
        do j=1,2
          do i=1,nrflav
            jcpref(i,j,ir)=jc(i,j)
          enddo
        enddo
      elseif(mod.eq.-2)then
        do j=1,2
          do i=1,nrflav
            jctref(i,j,ir)=jc(i,j)
          enddo
        enddo
      elseif(mod.eq.1)then
        do j=1,2
          do i=1,nrflav
            jc(i,j)=jcpref(i,j,ir)
          enddo
        enddo
      elseif(mod.eq.2)then
        do j=1,2
          do i=1,nrflav
            jc(i,j)=jctref(i,j,ir)
          enddo
        enddo
      elseif(mod.eq.10)then
        ic(1)=icproj(1,ir)
        ic(2)=icproj(2,ir)
        call iddeco(ic,jc)
        itest=0
        do j=1,2
          do i=1,nrflav
            jcpref(i,j,ir)=jcpref(i,j,ir)+jc(i,j)
          enddo
        enddo
c cancel quark and antiquarks to avoid to much remnant excitation
        do i=1,nrflav
          if(jcpref(i,1,ir).ge.jcpref(i,2,ir))then
            jcpref(i,1,ir)=jcpref(i,1,ir)-jcpref(i,2,ir)
            jcpref(i,2,ir)=0
c update valence quarks (cancel first sea quarks)
            if(jcpref(i,1,ir)-jc(i,1).lt.0)jc(i,1)=jcpref(i,1,ir)
            jc(i,2)=0
          else
            jcpref(i,2,ir)=jcpref(i,2,ir)-jcpref(i,1,ir)
            jcpref(i,1,ir)=0
c update valence quarks (cancel first sea quarks)
            if(jcpref(i,2,ir)-jc(i,2).lt.0)jc(i,2)=jcpref(i,2,ir)
            jc(i,1)=0
          endif
          do j=1,2
            itest=itest+jcpref(i,j,ir)
            jc2(i,j)=jcpref(i,j,ir)
          enddo
        enddo
        if(itest.eq.0)then !do not leave empty remnant
          idum=idrafl(iclpro,jc2,1,'r',3,iretso)     !create q-qb
          do j=1,2
            do i=1,nrflav
              jcpref(i,j,ir)=jc2(i,j)
            enddo
          enddo
        endif
      if(ish.ge.6)write(ifch,'(a,i3,a,1x,5i3,3x,5i3)')
     & 'jcpref(',ir,') ini:',((jcpref(i,j,ir),i=1,nflavems),j=1,2)
      elseif(mod.eq.20)then
        ic(1)=ictarg(1,ir)
        ic(2)=ictarg(2,ir)
        call iddeco(ic,jc)
        itest=0
        do j=1,2
          do i=1,nrflav
            jctref(i,j,ir)=jctref(i,j,ir)+jc(i,j)
          enddo
        enddo
c cancel quark and antiquarks to avoid to much remnant excitation
        do i=1,nrflav
          if(jctref(i,1,ir).ge.jctref(i,2,ir))then
            jctref(i,1,ir)=jctref(i,1,ir)-jctref(i,2,ir)
            jctref(i,2,ir)=0
c update valence quarks (cancel first sea quarks)
            if(jctref(i,1,ir)-jc(i,1).lt.0)jc(i,1)=jctref(i,1,ir)
            jc(i,2)=0
          else
            jctref(i,2,ir)=jctref(i,2,ir)-jctref(i,1,ir)
            jctref(i,1,ir)=0
c update valence quarks (cancel first sea quarks)
            if(jctref(i,2,ir)-jc(i,2).lt.0)jc(i,2)=jctref(i,2,ir)
            jc(i,1)=0
          endif
          do j=1,2
            itest=itest+jctref(i,j,ir)
            jc2(i,j)=jctref(i,j,ir)
          enddo
        enddo
        if(itest.eq.0)then !do not leave empty remnant
          idum=idrafl(icltar,jc2,1,'r',3,iretso)     !create q-qb
          do j=1,2
            do i=1,nrflav
              jctref(i,j,ir)=jc2(i,j)
            enddo
          enddo
        endif
      if(ish.ge.6)write(ifch,'(a,i3,a,1x,5i3,3x,5i3)')
     & 'jctref(',ir,') ini:',((jctref(i,j,ir),i=1,nflavems),j=1,2)
      else
        stop'mod not recognized in UpdateFlav'
      endif
      end

c-------------------------------------------------------------------------
      subroutine ProPoDif(k,n)
c-------------------------------------------------------------------------
c Define excitation and mass of remnant for each Pomeron.
c Diffractive selection already done in ProPoDif and idpfpr is defined.
c Select low mass diagram to be used for low mass diffraction (only remnant
c or central reggeon). Remnant only is needed to have large rapidity 
c gaps and central reggeon for low multiplicity events. But remnants should
c not be used to get low multiplicity at mid rapidity because it push
c remnant particles to too large rapidities: not good for baryon ratio and
c for energy flow (too low at large y and problem with trigger for LHCb for
c instance). It is necessary to have low mass Pomerons with large rapidity
c to increase energy flow and dN/deta on one side without particle production
c on the other side. These diffracive "enhanced" diagrams can be connected on 
c some of the side only to create rapidity gaps (intermediate size) and hard
c diffraction. 
c "Normal" Pomerons are connected on both side always not to have
c too many small rapidity gaps (and included enhanced diagrams connected on
c both side which doesn't have a Regge like behavior)
c If a diff. interaction is selected then the mass is
c sampled again as (M2)**(-alppom) (for one (SD) or both remnants (DD) or
c for the Pomeron itself to produce a unique central string or resonance (CD) 
c keeping the same x+/x-). If the new mass is larger then the original one,
c then the Pomeron is kept and the difference with the original value of x 
c is given to the remnant.
c In case of CD if the mass is lower than the Pom mass then the Pom become
c a single string.
c For each existing Pomeron a low mass is given to the remnant but it is
c cumulative (large number of Pom on remnant = large mass) and minimum  mass
c is given by resonance mass. The mass is taken from the Pomeron itself.
c "hard" here means a Pomeron is emitted (soft or hard)
c-------------------------------------------------------------------------

#include "aaa.h"
#include "ems.h"
#include "sem.h"
      common/cems13/xvpr(0:3)
      double precision drangen,x,sx,xvpr,xmax,xmin,aremnex,xxx(3),rr
     .,xpy,xpx,xxxx,xxmin,xmin0!,w(0:9),wd,aks,wdh
      double precision plc,s,eps,alptr
      common/cems5/plc,s
      common/cemsr5/alptr(0:1,0:5)
      dimension iept(3)
      character cremn*6
      parameter (eps=1.d-30)
      logical lhard,lpom,lpomi,lndif,lmass

      if(idpr(n,k).eq.0.or.ivpr(n,k).ne.1)return

      ip=iproj(k)
      it=itarg(k)

c check for reggeon (low energy DD with pion exchange)

      if(idpr(n,k).eq.2)then  !Reggeon

        itpr(k)=-2*npommx    !Pair is nondiffractive
        iept(1)=40
        iept(2)=40
        call VirPom(k,n,-1)      !suppress reggeon and give mass to remnants

      else

c      do i=0,9
c        w(i)=0.d0
c      enddo
c      wd=0d0
c        call WomTyDif(w,n,k)
        xxmin=xvpr(1)
        lhard=.false.
        if(idpr(n,k).eq.3)lhard=.true. !if hard Pom (but not soft+gg)
        
c      if(iomega.lt.2)then
c        if(lhard)
c     .  xxmin=xxmin+4.d0*dble(max(q2kmin(1,n,k),q2kmin(2,n,k)))/s
c        print *,idfpr(n,k),lhard,xxmin.lt.xpr(n,k)
c       if(idhpr(n,k).eq.1.or.idhpr(n,k).eq.3)w(2)=0d0 !Y+ not possible (q at targ)
c          if(idhpr(n,k).ge.2)w(3)=0d0 !Y- not possible (q at proj)
c          if(idhpr(n,k).ne.0)w(4)=0d0 !X not possible (q on both side)
c          wdh=1d0+w(2)+w(3)+w(4)
c          w(1)=dble(q2kmin(0,n,k))    !hard contribution
c          wd=0d0
c          do i=2,4
c            w(i)=w(i)*w(1)/wdh
c            wd=wd+w(i)
c          enddo
c          aks=drangen(w(0))*w(1)
c        else  ! fit-hard=w(1)+w(5)+soft=w(0)+w(5)+w(6)+w(7)+soft and w(0)=w(8)+w(9)
c          lhard=.false.
c          w(1)=w(1)+w(5)+dble(q2kmin(0,n,k)) 
c          wd=0d0
c          do i=9,9      !DD only need to be fixed
c            wd=wd+w(i)
c          enddo
c          aks=drangen(w(0))*wd
c        endif
        if(ish.ge.5)then
c          if(lhard)then
          write(ifch,*)'ProPoDif: start ',k,n,ip
     *      ,iep(ip),it,iet(it),idhpr(n,k),idfpr(n,k),idpr(n,k),npr(3,k)
     *      ,q2kmin(1,n,k),q2kmin(2,n,k),xxmin,xppr(n,k),xmpr(n,k),lhard
cc     *  ,wd/w(1),w(2)/w(1)*100.d0,w(3)/w(1)*100.d0,w(4)/w(1)*100.d0
c          elseif(wd.gt.0.)then
c            
c            write(ifch,*)'ProPoDif:soft diff. prob.(%) '
cc        print *,'ProPoDif:diff probability(%) '
c     *       ,k,n,idpr(n,k),ip,iep(ip),it,iet(it),wd/w(1)
c     *       ,w(5)/wd*100.d0,w(6)/wd*100.d0,w(7)/wd*100.d0
c     *       ,w(8)/wd*100.d0,w(9)/wd*100.d0
c          else
c            write(ifch,*)'ProPoDif:soft '
c     *       ,k,n,idpr(n,k),ip,iep(ip),it,iet(it),wd/w(1)
c          endif
        endif
c      elseif(iomega.eq.3)then
c        aks=1d0
c      else
c        itpr(k)=-2*npommx    !Pair nondiffractive and without mass for remnant
c        idfpr(n,k)=1
c        goto 9999
c      endif

cc      wd=0d0      !use Pomerons if possible
c      if(idpr(n,k).eq.2.and.wd.gt.0d0)then      !low mass diffraction + R
c      
c        if(aks.le.w(5)+w(8))then   !R + soft CD
c
cc        if(aks.le.w(7))then
cc          idpr(n,k)=2           !reggeon (only one string or resonance)
cc         replace Pomeron by Reggeon
cc          npr(2,k)=npr(2,k)+1
cc          npr(1,k)=npr(1,k)-1
cc        else
cc          idpr(n,k)=1           !soft Pomeron
cc        endif
c          idfpr(n,k)=4          !no connection to proj or targ
c          itpr(k)=itpr(k)+2     !diffractive pair
c          goto 9999
c
cc      elseif(xpr(n,k).le.xvpr(1).and.(ldiff(1).or.ldiff(2)))then !low mass single diffraction
c        else                    !low mass diffraction
c
cc fix excitation but gives mass in ProRem with no other IP exchange or give mass
cc back to remnant (x->xr+ or x->xr- to transfer the proper (low) mass)
c          iept(1)=0
c          iept(2)=0
c          lpom=.false.
c          idfpr(n,k)=0          !no pomeron
c          aks=aks-w(5)-w(8)
cc        if(ldiff(1))then !SD- (targ not excited=Mass+)
c          if(aks.le.w(7))then !SD+ (targ not excited=Mass+)
c            idpr(n,k)=5          
c            iept(1)=12
cc          ldiff(1)=.true.
cc          ldiff(2)=.false.
c          elseif(aks.lt.w(7)+w(6))then !SD- (proj not excited=Mass-)
c            idpr(n,k)=4
c            iept(2)=12
cc          ldiff(1)=.false.
cc          ldiff(2)=.true.
c          elseif(aks.lt.w(7)+w(6)+w(9))then !SD- (proj not excited=Mass-)
c            iept(1)=12
c            iept(2)=12
c            idpr(n,k)=6
cc          ldiff(1)=.true.
cc          ldiff(2)=.true.
c          else         !selected as diff but exceed real xs so replaced by soft
c            iept(1)=11
c            iept(2)=11
c            idfpr(n,k)=1        !connection on both side
cc update pom id to be sure to have 1 here
c            npr(idpr(n,k),k)=npr(idpr(n,k),k)-1
c            idpr(n,k)=1
c            npr(1,k)=npr(1,k)+1
c            lpom=.true.
c          endif
cc update iep/iet
cc        do jrem=1,2
cc          if(jrem.eq.1)then
cc            irem=ip
cc            icl=iclpro
cc          else
cc            irem=it
cc            icl=icltar
cc          endif
cc
cc          if(ldiff(jrem))then
cc
ccc            r=rangen()*rexndii(icl)
ccc            if(r.lt.rexpdif(icl))then
cc              iept(jrem)=12
ccc            elseif(r.lt.rexpdif(icl)+rexndi(icl))then
ccc              iept(jrem)=32
ccc            else
ccc              iept(jrem)=42
ccc            endif
cc
cc            if(iez(irem,jrem).gt.0)then !some mass already given to this remnant
cc              if(mod(iez(irem,jrem),10).eq.2)then !other diagrams connected
cc                iez(irem,jrem)=12 !only diffractive diagrams
cc              else
cc                iez(irem,jrem)=30 !already some inelastic connections
cc              endif
cc            elseif(iez(irem,jrem).lt.0)then !if no previous excitation, iept not changed
cc              iez(irem,jrem)=iept(jrem)
cc            else                !if already ine diagram
cc              iez(irem,jrem)=iept(jrem)/10
cc              iez(irem,jrem)=iez(irem,jrem)*10 !excitation is not only diff
cc            endif
cc
cc          else
cc              iept(jrem)=0
cc          endif
cc
cc        enddo
cc          call VirPom(k,n,1)    !suppress mass pomeron now to have all momentum available for new mass generation
c        endif
c        if(ish.ge.5)write(ifch,*)'ProPoDif diff: ',iept(1),iept(2)
c     .                                        ,idpr(n,k)
c
c      else                      !normal or enhanced (diff) Pomeron
c
cc fix connection side
cc        ldiff(1)=.false.
cc        ldiff(2)=.false.
c        iept(1)=31
c        iept(2)=31
c        if(lhard)then        !only for hard diagrams
c          if(aks.lt.w(2))then
cc            ldiff(1)=.true.
c            idfpr(n,k)=3        !connexion on target
c            iept(1)=0
c          elseif(aks.lt.w(2)+w(3))then !Y-
cc            ldiff(2)=.true.
c            idfpr(n,k)=2        !connection on projectile
c            iept(2)=0
c          elseif(aks.lt.w(2)+w(3)+w(4))then   !X
cc            ldiff(1)=.true.
cc            ldiff(2)=.true.
c            idfpr(n,k)=4        !no connexion
c            iept(1)=0
c            iept(2)=0
c          else
c            idfpr(n,k)=1        !connection on both side
cc         call virpom(k,n,1)
cc         goto 9999
c          endif
c        else                    !soft
c          idfpr(n,k)=1          !connection on both side
cc update pom id to be sure to have 1 here
c          npr(idpr(n,k),k)=npr(idpr(n,k),k)-1
c          idpr(n,k)=1
c          npr(1,k)=npr(1,k)+1
c        endif
c      endif
c      iept(1)=0
c      iept(2)=0
        lmass=.true.
      if(idpr(n,k).ne.7)then      !not DD
        lpom=.true.
        if(lhard.and.xpr(n,k).gt.xxmin.and.nprt(k).gt.1)lmass=.false. !do not resample x in case of hard diffraction except for diffractif event (one scattering only)
        if(idfpr(n,k).eq.3)then     !Y+
          iept(1)=0
          iept(2)=20
          if(lhard)iept(2)=23
          lndif=.false.
          lpom=.false.
        elseif(idfpr(n,k).eq.2)then !Y-
          iept(1)=20
          iept(2)=0
          if(lhard)iept(1)=23
          lndif=.false.
          lpom=.false.
        elseif(idfpr(n,k).eq.4)then !X
          iept(1)=0
          iept(2)=0
          lndif=.false.
        else
          iept(1)=11
          iept(2)=11
c          print *,zzremn(ip,1)-q2nmin,zzremn(it,2)-q2nmin,nprt(k)
c          if(rangen().lt.rexndi(iclpro)+zodinc*(zzremn(ip,1)-q2nmin))
cc          if(rangen().lt.rexndi(iclpro)+zodinc*(q2kmin(1,n,k)-q2nmin))
cc     .                                        *xpp(ip)/xppr(n,k))
c     .    iept(1)=31
c          if(rangen().lt.rexndi(icltar)+zodinc*(zzremn(it,2)-q2nmin))
cc          if(rangen().lt.rexndi(icltar)+zodinc*(q2kmin(2,n,k)-q2nmin))
cc     .                                        *xmt(it)/xmpr(n,k))
c     .    iept(2)=31
          if(lhard)then
            iept(1)=iept(1)+2
            iept(2)=iept(2)+2
          endif
          lndif=.true.
c if nuclear collision, no change in Pom momentum, remnant mass given in ProReM
          if(kolp(ip)+kolt(it).gt.4.or.iomega.eq.2)then
            itpr(k)=-2*npommx   !Pair nondiffractive
            idfpr(n,k)=1
            goto 8888
          endif
c        else
c for soft, share energy between Pom and remnant (to be more stable at low energy)
c          iept(1)=11
c          iept(2)=11
c          lndif=.true.
        endif
      else                      !DD
        if(iomega.eq.2)then
          itpr(k)=-2*npommx     !Pair nondiffractive if define here
          iept(1)=10
          iept(2)=10
        else
          iept(1)=22
          iept(2)=22
        endif
        lpom=.false.
        lndif=.false.
        lmass=.false.
      endif
      ntry=0
 100  continue
      ntry=ntry+1
      xxx(1)=0d0
      xxx(2)=0d0
      xxx(3)=0d0
      if(iomega.eq.1.or..not.lmass)goto 200
      lpomi=lpom
      if(idfpr(n,k).eq.4)then
c for CD the total mass is defined as 1/M2 not the remnant mass 
c so no loop is needed
        jrem1=3
        jrem2=3
        jrem3=1
      else
c Select excitation of remnant and generate corresponding mass
c Use random selection for DD not to set the mass of proj. always first
        jrem1=1+nint(rangen())
        jrem2=3-jrem1
        jrem3=1
        if(jrem2.lt.jrem1)jrem3=-1
      endif
        do jrem=jrem1,jrem2,jrem3
          xmax=1d0
          xpy=0d0
          if(jrem.eq.1)then
            cremn='proj: '
            icl=iclpro
            irem=ip
            amremn=amproj
            xmax=xmpr(n,k)
            xpx=xmax
            if(.not.lpom)then
              xmax=xmt(it)+xpx           !for DD (energy needed on both side)
              xpx=1d0
            elseif(idfpr(n,k).eq.2)then !for SD+
              xmax=min(xmt(it)+xpx,dble(xmxrem(2)))
              xpy=xppr(n,k)
            endif
            xmax=min(xmax,1d0-xmp(ip))   !in nuclear scattering we should not give too much energy to remnant)
          elseif(jrem.eq.2)then
            cremn='targ: '
            icl=icltar
            irem=it
            amremn=amtarg
            xmax=xppr(n,k)
            xpx=xmax
            if(.not.lpom)then
              xmax=xpp(ip)+xpx           !for DD (energy needed on both side)
              xpx=1d0
            elseif(idfpr(n,k).eq.3)then !for SD-
              xmax=min(xpp(ip)+xpx,dble(xmxrem(1)))
              xpy=xmpr(n,k)
            endif
            xmax=min(xmax,1d0-xpt(it))   !in nuclear scattering we should not give too much energy to remnant)
          else                  !for CD
            cremn='cent: '
            icl=iclpro
            irem=ip
            iept(jrem)=20
            xpx=xpr(n,k)
            xpy=1d0
            xmax=min(xpp(ip)+xppr(n,k),dble(xmxrem(1)))
            xmax=xmax*min(xmt(it)+xmpr(n,k),dble(xmxrem(2)))
          endif
cc no excitation possible only in certain cases (enhanced diagrams)
cc          rex=rexndii(icl)
cc          if(.not.ldiff(jrem))rex=1d0   
c          if(.not.ldiff(jrem))then   
cc          r=rangen()
cc          if(r.le.rexpdif(icl))then !low mass (resonance)
cc            iept(jrem)=40
cc            if(ldiff(jrem))iept(jrem)=42
cc          elseif(r.le.rex)then !"normal" mass
c            iept(jrem)=31
cc            if(ldiff(jrem))iept(jrem)=12
cc          elseif(r.le.rexndii(icl)/rex)then !high mass
cc            iept(jrem)=30       !no inv for high mass
cc            if(ldiff(jrem))iept(jrem)=32
c          else
c            iept(jrem)=0
c          endif
          ieptj=iept(jrem)/10
          if(ieptj.gt.0.and.ieptj.ne.3)then  !no or high mass fixed in ProReM
            
            aremnex=amemn(idz(irem,min(2,jrem)),ieptj)
            if(jrem.eq.1)then
              sx=s*(xpz(irem,jrem)+xpy)
              aremnex=aremnex+dble(amproj)
            elseif(jrem.eq.2)then
              sx=s*(xpz(irem,jrem)+xpy)
              aremnex=aremnex+dble(amtarg)
            else       !CD
              sx=s
            endif
            xmin0=aremnex**2d0/sx
c            xM2max=q2nmin
c            if(lndif)xM2max=q2kmin(jrem,n,k)
c should have the same in getdropx for string emission
            if(ieptj.ne.4.and.ieptj.ne.1)then !pion exchange or low mass excitation, minim should not change
              if(mod(iept(jrem),10).eq.2)then !diff excitation
                xmin=dble(xmindiff)*xmin0
              else
                xmin=dble(xminremn)*xmin0!+zmsinc*(xM2max-q2nmin)/sx
              endif
            else
              xmin=xmin0
            endif
c            if(lpom.and.jrem.le.2)xmax=min(xmax,
c     .           q2nmin*exp(min(50.,(xM2max-q2nmin)/zdfinc))/sx+xmin)
            if(lndif)xmax=min(xmax,dble(amdrmax)**2/sx+xmin)
            if(xmin.ge.xmax)xmin=xmin0
            if(xmin.le.xmax)then
              alp=alptr(idz(irem,min(2,jrem)),ieptj) 
c              alp=dble(alphigh(1+idz(irem,min(2,jrem))))
c        if(ieex.eq.3)alp=dble(alppom)
c              alp=dble(alppom)+(alp-dble(alppom))
c     .                        /(1d0-log(1d0-sx/s)*dble(zdfinc))
c              if(lpom.and.npr(3,k).gt.1)alp=alppom
c     .                             +max(0.,alp-alppom)/float(npr(3,k))
              rr=drangen(dble(alp))
              if(abs(alp-1.d0).lt.eps)then
                xxx(jrem)=xmax**rr*xmin**(1d0-rr)
              else
                xxx(jrem)=(rr*xmax**(1d0-alp)+(1d0-rr)*xmin**(1d0-alp))
     .                                                **(1d0/(1d0-alp))
              endif
            else
              xxx(jrem)=0d0     !mass fixed in ProRem
c              iept(jrem)=0
            endif
c if mass is high enough, then use high mass diff (with soft or hard pom)
c          if(xxx(jrem).lt.-xpx)then
            if(xxx(jrem).gt.xpx)then
              xxxx=xpx
              if(jrem.ne.3)xxxx=xpy*xxxx !x of diff Pom
              if(ish.ge.5)write(ifch,*)'Soft high mass diff'
     .                                 ,xxxx,xxxx.gt.max(xvpr(3),xxmin)
              if(xxxx.gt.max(xvpr(3),xxmin))then
c mass is large enough, keep the Pomeron to have high mass diffraction
                xxx(jrem)=xpx    !to avoid too large eta extension of diffractive events
                if(jrem.eq.2)then ! Y+ connexion on target
                  xxx(2)=xxx(2)-xppr(n,k)
                  xppr(n,k)=xppr(n,k)+xxx(2) !to keep the same x+_pom at the end
                  xpp(ip)=xpp(ip)-xxx(2)
                elseif(jrem.eq.1)then !Y- connection on projectile
                  xxx(1)=xxx(1)-xmpr(n,k)
                  xmpr(n,k)=xmpr(n,k)+xxx(1) !to keep the same x-_pom at the end
                  xmt(it)=xmt(it)-xxx(1)
                else            !X   (do nothing, just keep Pomeron as it is
                  xxx(1)=0d0
                  xxx(2)=0d0
                  xxx(3)=xpx
                endif
              else
                lpom=.false.
                if(jrem.eq.3)then !for CD cannot use momentum larger than the one from the Pom
                  xxx(1)=0d0
                  xxx(2)=0d0
                  xxx(3)=xpx
                endif
              endif
            elseif(idfpr(n,k).ne.1)then !low mass diffraction selected
              lpom=.false.
            endif
          else
            xxx(jrem)=0d0
          endif
          if(ish.ge.5)write(ifch,*)'ProPoDif mass ',cremn,iept(jrem)
     .         ,xmin,xxx(jrem),xmax,xpx,rexndi(icl),lndif,lpom,lpomi

        enddo

 200    continue

        if(lpom)then
c compare mass to available energy
          xxx(1)=min(xxx(1),xmpr(n,k))
          xxx(2)=min(xxx(2),xppr(n,k))
          x=(xmpr(n,k)-xxx(1))
c        if(xx.le.0d0)then
c          if(ish.ge.5)
c     .    write(ifch,*)'Not enough momentum for proj...',ntry
c          if(ntry.le.10.and..not.ldiff(1))goto 100
c          iept(1)=0
c          xx=xmpr(n,k)
c          xxx(1)=0d0
c        endif
          x=x*(xppr(n,k)-xxx(2))
c        if(x.le.0d0)then
c          if(ish.ge.5)
c     .    write(ifch,*)'Not enough momentum for targ...',ntry
c          if(ntry.le.10.and..not.ldiff(2))goto 100
c          iept(2)=0
c          x=xx*xppr(n,k)
c          xxx(2)=0d0
c        endif
          if(xpr(n,k).gt.xxmin.and.x.le.xxmin.and.ntry.le.100)then !check minimum mass of Pomeron
            if(ish.ge.6)
     .           write(ifch,*)'Not enough mass for Pom...',ntry
            goto 100
          elseif(x.gt.xxmin.or.
     .          (x.gt.xvpr(3).and.xxx(1)*xxx(2).le.0d0.and.lndif))then !check minimum mass for Pomeron or if remnant mass at its maximum
            itpr(k)=-2*npommx   !pair with Pomeron
c            if(.not.lndif)then
c              npr(idpr(n,k),k)=npr(idpr(n,k),k)-1
c              idpr(n,k)=1
c              npr(idpr(n,k),k)=npr(idpr(n,k),k)+1
c            endif
            xmpr(n,k)=xmpr(n,k)-xxx(1)
            xmp(ip)=xmp(ip)+xxx(1)
            xppr(n,k)=xppr(n,k)-xxx(2)
            xpt(it)=xpt(it)+xxx(2)
            xpr(n,k)=xppr(n,k)*xmpr(n,k)
            ypr(n,k)=0.5D0*log(xppr(n,k)/xmpr(n,k))
c to avoid too large eta extension of diffractive events the remnants have low mass in case of high mass diffraction
            if(iept(1).ge.20)iept(1)=iept(1)-10
            if(iept(1).ge.20)iept(2)=iept(2)-10
            if(nbkpr(n,k).ne.0)then !update backup
              nn=nbkpr(n,k)
              idfpr(nn,k)=idfpr(n,k)
              xpr(nn,k)=xpr(n,k)
              ypr(nn,k)=ypr(n,k)
              xppr(nn,k)=xppr(n,k)
              xmpr(nn,k)=xmpr(n,k)
            endif
          else                  !not enough momentum to have excitation
c no Pomeron possible : give momentum to mass on both side
c update pom id to be sure to have 2 here
            if(ish.ge.5)
     .      write(ifch,*)'ProPoDif: low mass',idfpr(n,k),xpr(n,k)
c            npr(idpr(n,k),k)=npr(idpr(n,k),k)-1
c            npr(1,k)=npr(1,k)+1        !to be suppressed in Virpom
            if(idfpr(n,k).eq.2)then !if mass on proj side
              xxx(1)=min(xmpr(n,k),1d0-xmp(ip))
            elseif(idfpr(n,k).eq.3)then !if mass on targ side
              xxx(2)=min(xppr(n,k),1d0-xpt(it))
c            elseif(idfpr(n,k).eq.4.or.xpr(n,k).gt.xvpr(2))then   !CD
            elseif(idfpr(n,k).eq.4.or.rangen().lt.rexddf)then   !CD
c            else   !CD or normal Pom becomes CD
              xxx(3)=xpr(n,k)
              iept(1)=0
              iept(2)=0
              idfpr(n,k)=4
            else                      !Normal Pom becomes DD
              xxx(1)=min(xmpr(n,k),1d0-xmp(ip))
              xxx(2)=min(xppr(n,k),1d0-xpt(it))
            endif
            lpom=.false.
            goto 200
          endif
        elseif(idfpr(n,k).eq.4.and.xxx(3).gt.xvpr(2))then
c for low mass CD: change from 2 strings to 1 string (or resonance) and update mass and momentum
          xxx(1)=sqrt(xxx(3))*exp(ypr(n,k))
          xxx(2)=sqrt(xxx(3))*exp(-ypr(n,k))
          x=max(0d0,xpp(ip)+xppr(n,k)-xxx(1))
     .     *max(0d0,xmt(it)+xmpr(n,k)-xxx(2))
          if(x.gt.0d0)then
            itpr(k)=-2*npommx   !pair with Pomeron
            npr(idpr(n,k),k)=npr(idpr(n,k),k)-1
            npr(1,k)=npr(1,k)+1 
            idpr(n,k)=6
            xpr(n,k)=xxx(3)
c update momentum on projectile side
            xpp(ip)=xpp(ip)+xppr(n,k)
            xppr(n,k)=xxx(1)
            xpp(ip)=xpp(ip)-xppr(n,k)
c update momentum on target side
            xmt(it)=xmt(it)+xmpr(n,k)
            xmpr(n,k)=xxx(2)
            xmt(it)=xmt(it)-xmpr(n,k)
            if(nbkpr(n,k).ne.0)then !delete backup Pom
              call VirPom(k,nbkpr(n,k),5)
            endif
          else
            if(ish.ge.5)
     .      write(ifch,*)'ProPoDif: CD rejected',xxx
            xxx(3)=0d0
            goto 200
          endif
        else
          if(idpr(n,k).le.3)then
            npr(idpr(n,k),k)=npr(idpr(n,k),k)-1
            npr(1,k)=npr(1,k)+1 !to be suppressed in Virpom
            if(nbkpr(n,k).ne.0)then     !delete backup Pom
              call VirPom(k,nbkpr(n,k),7)
            endif
          endif
          if(idfpr(n,k).eq.3)then !if mass on targ side
            idpr(n,k)=4
            iept(2)=22
          elseif(idfpr(n,k).eq.2)then !if mass on proj side
            idpr(n,k)=5
            iept(1)=22
c         elseif(idfpr(n,k).eq.4)then !CD
c            idpr(n,k)=6
          else                        !DD
            idpr(n,k)=7
            if(iept(1).ne.22)then
              iept(1)=31
c              if(kolt(it).gt.1)iept(1)=51   !if projectile surrounded by other wounded nucleons
            endif
            if(iept(2).ne.22)then
              iept(2)=31
c              if(kolp(ip).gt.1)iept(2)=51   !if target surrounded by other wounded nucleons
            endif
          endif
          call VirPom(k,n,1)    !suppress mass pomeron now since OK
          if(iept(1).eq.31.or.iept(2).eq.31)then
            itpr(k)=-2*npommx   !Pair nondiffractive and fix mass in ProRem (too large otherwise)
          else
            itpr(k)=itpr(k)+2   !diffractive pair
            if(iept(1).ne.0)then !transfer new mass to projectile
              xmt(it)=xmt(it)-xxx(1)
              xmp(ip)=xmp(ip)+xxx(1)
            endif
            if(iept(2).ne.0)then !transfer new mass to target
              xpp(ip)=xpp(ip)-xxx(2)
              xpt(it)=xpt(it)+xxx(2)
            endif
          endif
          idfpr(n,k)=0          !no pomeron
        endif
        if(ish.ge.5)write(ifch,*)'ProPoDif final: ',idfpr(n,k)
     .           ,xppr(n,k),xmpr(n,k),sqrt(xpp(ip)*xmp(ip)*s)
     .                               ,sqrt(xpt(it)*xmt(it)*s)
      endif                     !end no Reggeon
c ------------ end mass loop ---------------------------------------
 8888 continue
      if(idpr(n,k).gt.3.and.idpr(n,k).ne.6)call VirPom(k,n,2)    !security
      do jrem=1,2
        if(jrem.eq.1)then
          irem=ip
        else
          irem=it
        endif
        if(iept(jrem).gt.0)then
          ieex=mod(iez(irem,jrem),10)

c          if(ieex.eq.3.and.mod(iept(jrem),10).eq.3)then !some hard mass already given to this remnant
c            iez(irem,jrem)=53
c          elseif(iept(jrem)/10.eq.2.and.iez(irem,jrem)/10.eq.2)then
          if(iept(jrem)/10.eq.2.and.iez(irem,jrem)/10.eq.2)then
            if(mod(iept(jrem),10).eq.3.or.ieex.eq.3)then
              iez(irem,jrem)=23
            elseif(mod(iept(jrem),10).eq.2.and.ieex.eq.2)then
              iez(irem,jrem)=22
            else
              iez(irem,jrem)=20
            endif
c          elseif(iez(irem,jrem).gt.0)then
cc            iez(irem,jrem)=max(iez(irem,jrem)/10,iept(jrem)/10)*10
cc            if(ieex.eq.3)iez(irem,jrem)=iez(irem,jrem)+3
cc            if(iez(irem,jrem)/10.eq.5.or.iept(jrem)/10.eq.5)then
cc              iez(irem,jrem)=51
cc              if(ieex.eq.3.or.mod(iept(jrem),10).eq.3)iez(irem,jrem)=53
cc            elseif(iez(irem,jrem)/10.gt.1.or.iept(jrem)/10.gt.1)then
cc              iez(irem,jrem)=31
cc              if(ieex.eq.3.or.mod(iept(jrem),10).eq.3)iez(irem,jrem)=33
cc            else
cc              iez(irem,jrem)=11
cc              if(ieex.eq.3.or.mod(iept(jrem),10).eq.3)iez(irem,jrem)=13
cc            endif
          else
            iez(irem,jrem)=iept(jrem)
          endif
        else
c          if(mod(iez(irem,jrem),10).eq.1
c     .       .or.mod(iez(irem,jrem),10).eq.2)then !inversion possible for this remnant
c            iez(irem,jrem)=iez(irem,jrem)/10 !no inv if multiple diag
c            iez(irem,jrem)=iez(irem,jrem)*10 !and diagram not only diff
c          elseif(idfpr(n,k).ne.0.and.iez(irem,jrem).lt.0)then !if a Pomeron is emitted
c          if(idfpr(n,k).ne.0.and.iez(irem,jrem).lt.0)then !if a Pomeron is emitted
          if(iez(irem,jrem).lt.0)then !if there is an interaction (even diff)
            iez(irem,jrem)=0
          endif
        endif
      enddo

c 9999 continue

      if(ish.ge.5)write(ifch,*)'ProPoDif: itpr,idpr,iep,iet '
     .                         ,itpr(k),idpr(n,k),iep(ip),iet(it)

      return
      end


c-------------------------------------------------------------------------
      subroutine ProItfpr(k)
c-------------------------------------------------------------------------
c propose diffractive scattering identification
c called if a Pomeron or Reggeon is present
c fix idfpr for each Pomeron
c-------------------------------------------------------------------------
#include "aaa.h"
#include "ems.h"

      ip=iproj(k)
      it=itarg(k)
      idf=1
      if(iep(ip).le.0.and.iet(it).le.0)then
        itpr(k)=abs(itpr(k))    !central diffractive type interaction with Pomeron (1) or reggeon (2)
c        idf=4
      elseif(iep(ip)*iet(it).le.0)then
        if(itpr(k).eq.-1)itpr(k)=abs(itpr(k))    !diffractive type interaction with Pom
c        if(iep(ip).gt.0)then    !if only on projectile side
c          idf=2                 !connection on projectile
c        elseif(iet(it).gt.0)then !if only on target side
c          idf=3                 !connexion on target
c        endif
      endif

      end

c-------------------------------------------------------------------------
      subroutine ProReEx(ir,ii)
c-------------------------------------------------------------------------
c fix remnant excitation
c for proj (iep) if ir=1 or target (iet) if ir=-1:
c 0 = no,  11,30,40,50 = inel excitation,  22,32,42 = diffr excitation,
c fixed before : 50 = excit due to split without connection
c fixed after : 50 = large excitation due to # quark > 3
c               60 = active spectator (get small pt and used for mass)
c if mod(iept,10)=0 there is no diquark/quark inversion in ProRef
c if mod(iept,10)=2 the remnant is connected only to diffractive pairs
c-------------------------------------------------------------------------
#include "aaa.h"
#include "ems.h"
      common/cnparticip/jproj(2,mamx),jtarg(2,mamx),efluct(6,mamx)


c      xmndif=1.
      mndif=0
      mhard=0
      midff=0
      if(ir.eq.1)then                   !proj
c        iep(ii)=max(0,iep(ii))  !set iep to at least 0 for all participants
        if(iep(ii).eq.22)then   !diff
          mdiff=jproj(1,ii)
c          do li=1,lproj(ii)
c            if(rangen().lt.rexdif(iclpro))mdiff=mdiff+1
c          enddo          
          if(mdiff.lt.0)iep(ii)=23
        elseif(iep(ii).gt.0.and.iep(ii).ne.2)then
cc check if remnant connected to at least one pair with non 0 Q2s
          do li=1,lproj(ii)
            kk=kproj(ii,li)
            do n=1,nprmx(kk)
              if((idfpr(n,kk).eq.1.or.idfpr(n,kk).eq.2)
     .                            .and.ivpr(n,kk).eq.1)then
                dqs=q2kmin(1,n,kk)-q2nmin
                if(rangen().lt.rexndi(iclpro)+zodinc*dqs)mndif=mndif+1
c                xmndif=xmndif*max(0.,1.-(rexndi(iclpro)+zodinc*dqs))
                if(dqs.gt.0.)mhard=mhard+1
              endif
            enddo
          enddo
c          if(rangen().lt.xmndif)mndif=0
          if(mndif.eq.1)then
            iep(ii)=31
          elseif(mndif.gt.1)then
            iep(ii)=51
          endif
          if(mndif.gt.0.and.mhard.gt.0)iep(ii)=iep(ii)+2
        endif
      elseif(ir.eq.-1)then                !targ
c        iet(ii)=max(0,iet(ii))  !set iet to at least 0 for all participants
        if(iet(ii).eq.22)then   !diff
          mdiff=jtarg(1,ii)
c          do li=1,ltarg(ii)
c            if(rangen().lt.rexdif(icltar))mdiff=mdiff+1
c          enddo          
          if(mdiff.lt.0)iet(ii)=23
        elseif(iet(ii).gt.0.and.iet(ii).ne.2)then
cc check if remnant connected to at least one pair with non 0 Q2s
          do li=1,ltarg(ii)
            kk=ktarg(ii,li)
            do n=1,nprmx(kk)
              if((idfpr(n,kk).eq.1.or.idfpr(n,kk).eq.3)
     .                            .and.ivpr(n,kk).eq.1)then
                dqs=q2kmin(2,n,kk)-q2nmin
                if(rangen().lt.rexndi(icltar)+zodinc*dqs)mndif=mndif+1
c                xmndif=xmndif*max(0.,1.-(rexndi(icltar)+zodinc*dqs))
                if(dqs.gt.0.)mhard=mhard+1
              endif
            enddo
          enddo
c          if(rangen().lt.xmndif)mndif=0
          if(mndif.eq.1)then
            iet(ii)=31
          elseif(mndif.gt.1)then
            iet(ii)=51
          endif
          if(mndif.gt.0.and.mhard.gt.0)iet(ii)=iet(ii)+2
        endif
      endif

      end


c-------------------------------------------------------------------------
      subroutine ProDiPt(k,iqq,iret)
c-------------------------------------------------------------------------
c propose transverse momentum for diffractive interaction
c iqq=1  : fix pt for non diffractive pair (and check if all pairs are still valid)
c iqq=2  : diffractive pt with mass dependence
c-------------------------------------------------------------------------

#include "ems.h"
#include "sem.h"
#include "aaa.h"
      double precision xxe(kollmx),xye(kollmx),pt2,am0,am1,am2!,p5sqpr,p5sqtg
      double precision plc,s,xxpnew,xypnew,xxtnew,xytnew,utdble,RANNORM
      common/cems5/plc,s
      common/emsptl/nppr(npommx,kollmx),npproj(mamx),nptarg(mamx)
      save xxe,xye

      ip=iproj(k)
      it=itarg(k)


c generate p_t for diffractive
      if(iqq.eq.1)then

       if(ptdiff.gt.0.)then
         if(itpr(k).eq.-2)then
           pt=ranpt()*ptdiff/(1.+0.02*max(0.,sngl(log(s))))
c not good to get proper fragment mass in nuclear collisions
          elseif(itpr(k).eq.0)then   !pt for non-wounded nucleon (usefull in ProRem to avoid problem in utrescl)
c only to avoid some remnant to have only one partner for mass distribution
           
           pt = sngl(RANNORM(0.088D0,0.044D0)) !limited by some data like sal.optns
           if(kolp(ip).eq.1.and.iet(it).lt.0)then
             iet(it)=60         !active spectators
           elseif(kolt(it).eq.1.and.iep(ip).lt.0)then
             iep(ip)=60
           else
             pt=0.
           endif
         else
           xxe(k)=0d0
           xye(k)=0d0
           goto 10
         endif
         phi=2.*pi*rangen()
         xxe(k)=dble(pt*cos(phi))
         xye(k)=dble(pt*sin(phi))
       else
         xxe(k)=0d0
         xye(k)=0d0
       endif

c update remnant p_t

 10    xxp(ip)=xxp(ip)-xxe(k)
       xyp(ip)=xyp(ip)-xye(k)
       xxt(it)=xxt(it)+xxe(k)
       xyt(it)=xyt(it)+xye(k)

       if(itpr(k).ne.0.and.abs(itpr(k)).ne.3)iret=0
!to simulate the fact that originally we had a Pomeron
c         if(koll.le.2)then
c           call StoCon(-k,k,1)  !to fixe mass of corresponding remnants
c           xpp(ip)=xpp(ip)-xppr(1,k)
c           xpt(it)=xpt(it)+xppr(1,k)
c           xmt(it)=xmt(it)-xmpr(1,k)
c           xmp(ip)=xmp(ip)+xmpr(1,k)
c           idpr(1,k)=0
c           xpr(1,k)=0.d0
c           ypr(1,k)=0.d0
c           xppr(1,k)=0.d0
c           xmpr(1,k)=0.d0
c         endif
c         p5sqpr=xpp(ip)*xmp(ip)*s-dble(amproj*amproj)
c         p5sqtg=xpt(it)*xmt(it)*s-dble(amtarg*amtarg)
c         phi=2.*pi*rangen()
c         ntry=0
c 20      ntry=ntry+1
c         pt=ranptdcut(ptsems)*ptsend**2
c         if(ntry.lt.100.and.(p5sqpr-dble(pt*pt).lt.0.d0
c     &                   .or.p5sqtg-dble(pt*pt).lt.0.d0))then
c             goto 20
c         else
c           pt=ranptdcut(ptsems)*ptsendi
c         endif
c         xxe(k)=dble(pt*cos(phi))
c         xye(k)=dble(pt*sin(phi))
c         xxp(ip)=xxp(ip)-xxe(k)
c         xyp(ip)=xyp(ip)-xye(k)
c         xxt(it)=xxt(it)+xxe(k)
c         xyt(it)=xyt(it)+xye(k)
c       endif

      elseif(itpr(k).eq.-2.and.ptdiff.ne.0.)then

        pt2=xxe(k)*xxe(k)+xye(k)*xye(k)
        if(pt2.gt.0d0)then
          am0=dble(amproj**2*amtarg**2)
          am1=max(dble(amproj**2),xpp(ip)*xmp(ip)*s
     &              -xxp(ip)*xxp(ip)-xyp(ip)*xyp(ip))
          am2=max(dble(amtarg**2),xpt(it)*xmt(it)*s
     &              -xxp(it)*xxp(it)-xyp(it)*xyp(it))
          ptd=ptdiff/(1.+0.02*max(0.,sngl(log(s*am0/am1/am2)))) !0.02 comes from data (Z. Phys. C 67, 227-237, 1995)
c           ad=pi/4./ptd**2
c           r=rangen()
          pt=ranpt()*ptd        !sqrt(-alog(r)/ad)
        else
          return
        endif
        if(ish.ge.8)write(ifch,'(a,i5,2i4,5g13.5)')
     &                    'ProDiPt',k,ip,it,pt,sqrt(pt2),ptd,am1,am2
c suppress the pt given with iqq=1 and give a new one taking into account the mass (iqq=2) with the same angle phi
        pt=pt/sqrt(pt2)
        xxe(k)=xxe(k)*pt
        xye(k)=xye(k)*pt

c update remnant p_t if enough energy available
        xxpnew=xxp(ip)-xxe(k)
        xypnew=xyp(ip)-xye(k)
        xxtnew=xxt(it)+xxe(k)
        xytnew=xyt(it)+xye(k)
        if((iep(ip).eq.0.or.
     &      xpp(ip)*xmp(ip)*s-xxpnew*xxpnew-xypnew*xypnew
     &      .gt.1.3d0*utdble(pptl(5,npproj(ip)))**2)
     &.and.(iet(it).eq.0.or.
     &      xpt(it)*xmt(it)*s-xxtnew*xxtnew-xytnew*xytnew
     &      .gt.1.3d0*utdble(pptl(5,nptarg(it)))**2))then
          xxp(ip)=xxp(ip)-xxe(k)
          xyp(ip)=xyp(ip)-xye(k)
          xxt(it)=xxt(it)+xxe(k)
          xyt(it)=xyt(it)+xye(k)
        endif

       endif

       end

c-------------------------------------------------------------------------
      subroutine ProSePt(k,n)
c-------------------------------------------------------------------------
c propose transverse momentum for string ends
c-------------------------------------------------------------------------

#include "aaa.h"
#include "ems.h"

      if(ivpr(n,k).eq.2)return            !Backup Pomeron

      ip=iproj(k)
      it=itarg(k)
      amk0=1. ! included in ptsend !(qmass(1)+qmass(2)+qmass(3))/3.     !mass for mt distribution

      ptsecut=ptsecu        !cut for gaussian distribution (center around 0.4)

c generate p_t for string ends  (proj)

c
c      !---proj-----
        iepp=mod(iep(ip),10)
        ptsef=ptsend
        if(iepp.eq.0)ptsef=ptsendi
        ptsendx = ptsems
        ptsendy = ptsendx!*2       !TP20140605 not sure about factor 2 ????

      if(idp1pr(n,k).gt.0)then
         if(idp1pr(n,k).eq.4.or.idp1pr(n,k).eq.5)then   !diquarks
           amk1=amk0*ptsendy+qmass(0) !mass for mt distribution with bounding energy for diquark
         else
           amk1=amk0*ptsendx
         endif
c         if(iepp.eq.0)amk1=0.
         if(iepp.eq.0)then
           pt=ranptd()*ptsef
         else
           pt=ranptdcut(ptsecut)*ptsef
           pt=pt+amk1
         endif
c         pt=ranptdcut(ptsecut)*ptsef
c         pt=pt+amk1
c         pt=ranptd()*ptsef
c         pt=sqrt(pt*pt+amk1*amk1)
         phi=2.*pi*rangen()
         xxp1pr(n,k)=dble(pt*cos(phi))
         xyp1pr(n,k)=dble(pt*sin(phi))
      else
         xxp1pr(n,k)=0d0
         xyp1pr(n,k)=0d0
      endif
      if(idp2pr(n,k).gt.0)then
         if(idp2pr(n,k).eq.4.or.idp2pr(n,k).eq.5)then
           amk1=amk0*ptsendy+qmass(0) !mass for mt distribution with bounding energy for diquark
         else
           amk1=amk0*ptsendx
         endif
c         if(iepp.eq.0)amk1=0.
         if(iepp.eq.0)then
           pt=ranptd()*ptsef
         else
           pt=ranptdcut(ptsecut)*ptsef
           pt=pt+amk1
         endif
c         pt=ranptdcut(ptsecut)*ptsef
c         pt=pt+amk1
c         pt=ranptd()*ptsef
c         pt=sqrt(pt*pt+amk1*amk1)
         phi=2.*pi*rangen()
         xxp2pr(n,k)=dble(pt*cos(phi))
         xyp2pr(n,k)=dble(pt*sin(phi))
      else
         xxp2pr(n,k)=0d0
         xyp2pr(n,k)=0d0
      endif
c generate p_t for string ends  (targ)


c      !---targ-----
        iett=mod(iet(it),10)
        ptsef=ptsend
        if(iett.eq.0)ptsef=ptsendi
        ptsendx = ptsems
        ptsendy = ptsendx!*2       !TP20140605 not sure about factor 2 ????

      if(idm1pr(n,k).gt.0)then
         if(idm1pr(n,k).eq.4.or.idm1pr(n,k).eq.5)then
           amk1=amk0*ptsendy+qmass(0) !mass for mt distribution with bounding energy for diquark
         else
           amk1=amk0*ptsendx
         endif
c         if(iett.eq.0)amk1=0.
         if(iett.eq.0)then
           pt=ranptd()*ptsef
         else
           pt=ranptdcut(ptsecut)*ptsef
           pt=pt+amk1
         endif
c         pt=ranptdcut(ptsecut)*ptsef
c         pt=pt+amk1
c         pt=ranptd()*ptsef
c         pt=sqrt(pt*pt+amk1*amk1)
         phi=2.*pi*rangen()
         xxm1pr(n,k)=dble(pt*cos(phi))
         xym1pr(n,k)=dble(pt*sin(phi))
      else
         xxm1pr(n,k)=0d0
         xym1pr(n,k)=0d0
      endif
      if(idm2pr(n,k).gt.0)then
         if(idm2pr(n,k).eq.4.or.idm2pr(n,k).eq.5)then
           amk1=amk0*ptsendy+qmass(0) !mass for mt distribution with bounding energy for diquark
         else
           amk1=amk0*ptsendx
         endif
c         if(iett.eq.0)amk1=0.
         if(iett.eq.0)then
           pt=ranptd()*ptsef
         else
           pt=ranptdcut(ptsecut)*ptsef
           pt=pt+amk1
         endif
c         pt=ranptdcut(ptsecut)*ptsef
c         pt=pt+amk1
c         pt=ranptd()*ptsef
c         pt=sqrt(pt*pt+amk1*amk1)
         phi=2.*pi*rangen()
         xxm2pr(n,k)=dble(pt*cos(phi))
         xym2pr(n,k)=dble(pt*sin(phi))
      else
         xxm2pr(n,k)=0d0
         xym2pr(n,k)=0d0
      endif

c update remnant p_t (pomeron)
        xxp(ip)=xxp(ip)-xxp1pr(n,k)-xxp2pr(n,k)
        xyp(ip)=xyp(ip)-xyp1pr(n,k)-xyp2pr(n,k)
        xxt(it)=xxt(it)-xxm1pr(n,k)-xxm2pr(n,k)
        xyt(it)=xyt(it)-xym1pr(n,k)-xym2pr(n,k)

        if(ish.ge.6)then
        write(ifch,*) 'ProSePt'
        write(ifch,'(4i14/4d14.3/4d14.3/)')
     * idp1pr(n,k),idp2pr(n,k),idm1pr(n,k),idm2pr(n,k)
     *,xxp1pr(n,k),xxp2pr(n,k),xxm1pr(n,k),xxm2pr(n,k)
     *,xyp1pr(n,k),xyp2pr(n,k),xym1pr(n,k),xym2pr(n,k)
        endif

        end

c-----------------------------------------------------------------------
      subroutine ProSeX(k,n,iret)
c-----------------------------------------------------------------------
c calculates x of string ends
c-----------------------------------------------------------------------

#include "aaa.h"
#include "ems.h"
      common/cems5/plc,s
      double precision s,plc
      common/cems10/a(0:ntypmx),b(0:ntypmx),d(0:ntypmx)
      double precision a,b,d,drangen
     *,xp,xm,ap1,ap2,am1,am2,aamin1,aamin2,u
     *,xmn1,xmn2

      iret=0

      if(idpr(n,k).eq.0.or.ivpr(n,k).eq.0)return       !do it for all Pomeron incl hard to be used by the backup only

      if(idpr(n,k).eq.2)return    !reggeon with remnant mass only for the moment ???????????????

      if(idp1pr(n,k).eq.0.and.idp2pr(n,k).eq.0
     * .and.idm1pr(n,k).eq.0.and.idm2pr(n,k).eq.0)
     *call utstop('no Pomeron in ProSex&')

      xp=xppr(n,k)
      xm=xmpr(n,k)

      if(idpr(n,k).eq.6)then   !low mass CD (1 string)

      xmn1=0d0
      xmn2=0d0
      ntry=1
      xp1pr(n,k)=xp
      xm2pr(n,k)=xm

      else                     !Pomeron (2 strings)

      ap1=a(idp1pr(n,k))
      ap2=a(idp2pr(n,k))
      am1=a(idm1pr(n,k))
      am2=a(idm2pr(n,k))
      aamin1=ammn(idp1pr(n,k)+idm2pr(n,k))
      aamin2=ammn(idp2pr(n,k)+idm1pr(n,k))
      xmn1=(aamin1**2+(xxp1pr(n,k)+xxm2pr(n,k))**2
     &               +(xyp1pr(n,k)+xym2pr(n,k))**2)/s
      xmn2=(aamin2**2+(xxp2pr(n,k)+xxm1pr(n,k))**2
     &               +(xyp2pr(n,k)+xym1pr(n,k))**2)/s

      ntry=0
 999  ntry=ntry+1
      if(ntry.gt.100)then
        iret=1
        if(ish.ge.5)write(ifch,*)'Problem in ProSex(k,n)',k,n
        return
      endif

    1 u=drangen(ap1)**(1d0/(1d0+ap1))
      if(drangen(u).gt.(1d0-u)**ap2)goto1
      xp1pr(n,k)=u*xp
      xp2pr(n,k)=(1-u)*xp
    2 u=drangen(am1)**(1d0/(1d0+am1))
      if(drangen(u).gt.(1d0-u)**am2)goto2
      xm1pr(n,k)=u*xm
      xm2pr(n,k)=(1-u)*xm

      if(xp1pr(n,k)*xm2pr(n,k).lt.xmn1)then
      goto 999
c       fc=xp1pr(n,k)*xm2pr(n,k)/xmn1   !avoid virpom
c       if(fc.eq.0.)goto 999
c       xp1pr(n,k)=xp1pr(n,k)/sqrt(fc)
c       xm2pr(n,k)=xm2pr(n,k)/sqrt(fc)
      endif
      if(xp2pr(n,k)*xm1pr(n,k).lt.xmn2)then
      goto 999
c       fc=xp2pr(n,k)*xm1pr(n,k)/xmn2   !avoid virpom
c       if(fc.eq.0.)goto 999
c       xp2pr(n,k)=xp2pr(n,k)/sqrt(fc)
c       xm1pr(n,k)=xm1pr(n,k)/sqrt(fc)
      endif

      endif        !Reggeon/Pomeron

      if(ish.ge.6)then
       write(ifch,*) 'ProSeX'
       write(ifch,'(2d28.3,i8)') xp,xm,ntry
       write(ifch,'(4d14.3)')xp1pr(n,k),xp2pr(n,k),xm1pr(n,k),xm2pr(n,k)
       write(ifch,'(4d14.3/)')xp1pr(n,k)*xm2pr(n,k)
     *                   ,xp2pr(n,k)*xm1pr(n,k),  xmn1, xmn2
      endif

      end

c-------------------------------------------------------------------------
      subroutine RmPt(k,n)
c-------------------------------------------------------------------------
c remove pt from pomeron
c-------------------------------------------------------------------------
#include "aaa.h"
#include "ems.h"
      ip=iproj(k)
      it=itarg(k)
      xxp(ip)=xxp(ip)+xxp1pr(n,k)+xxp2pr(n,k)
      xyp(ip)=xyp(ip)+xyp1pr(n,k)+xyp2pr(n,k)
      xxt(it)=xxt(it)+xxm1pr(n,k)+xxm2pr(n,k)
      xyt(it)=xyt(it)+xym1pr(n,k)+xym2pr(n,k)
      xp1pr(n,k)=0d0
      xp2pr(n,k)=0d0
      xm1pr(n,k)=0d0
      xm2pr(n,k)=0d0
      xxm1pr(n,k)=0d0
      xym1pr(n,k)=0d0
      xxp1pr(n,k)=0d0
      xyp1pr(n,k)=0d0
      xxm2pr(n,k)=0d0
      xym2pr(n,k)=0d0
      xxp2pr(n,k)=0d0
      xyp2pr(n,k)=0d0
      idp1pr(n,k)=0
      idm2pr(n,k)=0
      idp2pr(n,k)=0
      idm1pr(n,k)=0
      end

c-------------------------------------------------------------------------
      subroutine VirPom(k,n,id)
c-------------------------------------------------------------------------
c create virtual pomeron
c virtual pomeron: ivpr(n,k)=0, otherwise ivpr(n,k)=1 
c (ivpr(n,k)=2 for backup)
c if id<0 x+ and x- are given to remnant mass
c if id>0 x+ and x- are given back to remnant momentum
c-------------------------------------------------------------------------

#include "aaa.h"
#include "ems.h"
#include "sem.h"
      common/emsptl/nppr(npommx,kollmx),npproj(mamx),nptarg(mamx)
      double precision plc,s
      common/cems5/plc,s
      common/col3/ncol,kolpt,ncoli
      integer jcp(nflav,2),jct(nflav,2)
c      dimension q2kk(2)
c      data nvir/0/
c      save nvir

      call utpri('VirPom',ish,ishini,3)

      ip=iproj(k)
      it=itarg(k)

      nnv=nvpr(n,k)
      nnb=nbkpr(n,k)

c                        nvir=nvir+1
c                   print *,'  ',id,'   ',nvir

      if(ish.ge.3)then
      write(ifch,*)"virpom ",id," (n,k)",n,k,nnb,nnv,nppr(n,k),itpr(k)
     &         ,nprt(k),idpr(n,k),ivpr(n,k),npr(1,k),npr(2,k),npr(3,k)
      if(ish.ge.5)write(ifch,*)"remnant in",xpp(ip),xmp(ip)
     &                                     ,xmt(it),xpt(it)
      endif

      if(idpr(n,k).eq.0.or.ivpr(n,k).eq.0)goto 1000


c if this is a virtual Pomeron (soft copy of an hard one)
      if(nnv.ne.0)then
        nn=nnv
        kk=k
c if the original is gone, treat as a normal pomeron later in the followings
        if(ivpr(nn,kk).eq.0)then
          nvpr(n,k)=0
        endif
      else
        nn=0
      endif

c if there is a copy of the hard Pomeron
      if(nnb.ne.0)then
        nn=nnb
        kb=k
c if the copy is gone already, treat as a normal Pomeron
        if(ivpr(nn,kb).eq.0)then
          nbkpr(n,k)=0
          nnb=0
        endif
      else
        kb=0
        nn=0
      endif


      if(nnb.ne.0)then     !there is a backup Pomeron
c copy all variables from the hard to the soft Pomeron
        nvpr(nn,kb)=0            !not virtual any more
        idpr(nn,kb)=1
        idhpr(nn,kb)=0
        ivpr(nn,kb)=1
        xpr(nn,kb)=xpr(n,k)
        ypr(nn,kb)=ypr(n,k)
        xppr(nn,kb)=xppr(n,k)
        xmpr(nn,kb)=xmpr(n,k)
        idfpr(nn,kb)=idfpr(n,k)
c        bhpr(nn,kb)=bhpr(n,k)
        idp1pr(nn,kb)=idp1pr(n,k)
        idp2pr(nn,kb)=idp2pr(n,k)
        idm1pr(nn,kb)=idm1pr(n,k)
        idm2pr(nn,kb)=idm2pr(n,k)
        xm1pr(nn,kb)=xm1pr(n,k)
        xp1pr(nn,kb)=xp1pr(n,k)
        xm2pr(nn,kb)=xm2pr(n,k)
        xp2pr(nn,kb)=xp2pr(n,k)
        xxm1pr(nn,kb)=xxm1pr(n,k)
        xym1pr(nn,kb)=xym1pr(n,k)
        xxp1pr(nn,kb)=xxp1pr(n,k)
        xyp1pr(nn,kb)=xyp1pr(n,k)
        xxm2pr(nn,kb)=xxm2pr(n,k)
        xym2pr(nn,kb)=xym2pr(n,k)
        xxp2pr(nn,kb)=xxp2pr(n,k)
        xyp2pr(nn,kb)=xyp2pr(n,k)

      elseif(nvpr(n,k).eq.0)then     !normal Pomeron 

c give everything back to remnant

      npr(0,k)=npr(0,k)+1
      idpom=idpr(n,k)
      if(idpom.gt.3)idpom=1
      npr(idpom,k)=npr(idpom,k)-1
      nprt(k)=npr(1,k)+npr(2,k)+npr(3,k)
      if(id.gt.0)then
c        npp(ip)=npp(ip)-1
c        npt(it)=npt(it)-1
        antotf=antotf-1
        if(idpr(n,k).eq.1)ansff=ansff-1
        if(idpr(n,k).eq.3)anshf=anshf-1
c not necessary since discarded remnant are marked with ivpr=0 but not removed
c      elseif(id.ne.0)then
cc include Pomeron in x-distribution before removal
c        if(iemspx.eq.1)then
c          q2kk(1)=0.
c          q2kk(2)=0.
c          jot=99
c          call xxEmsPx(1,jot,99,sngl(xpr(n,k)),sngl(ypr(n,k)),-1
c     *      ,0.,0.,q2kk )
c        endif
c        if(iemspbx.eq.1)then
c          q2kk(1)=0.
c          q2kk(2)=0.
c          call xxEmsP2(1,0,1,99
c     *            ,sngl(xppr(n,k)),sngl(xmpr(n,k))
c     *            ,0.,0.,0.,0.,0.,0.,0.,0.,  q2kk,q2kk )
c        endif
      endif
      kolp(ip)=kolp(ip)-1
      kolt(it)=kolt(it)-1

cc restore x from nuclear splitting
c      if(knucnt(1,k).gt.0)then
c        do nuc=1,knucnt(1,k)
c          if(npnuc(nuc,1,k).eq.n)then
c            ipp=irnuc(nuc,1,k)
c            xpp(ipp)=xpp(ipp)+xxnuc(nuc,1,k)
c            if(xpp(ipp)-1d0.ge.-1d-15)iep(ipp)=0
c            xppr(n,k)=xppr(n,k)-xxnuc(nuc,1,k)
c            xpr(n,k)=xppr(n,k)*xmpr(n,k)
c            ypr(n,k)=0.5D0*log(xppr(n,k)/xmpr(n,k))
c            npnuc(nuc,1,k)=0    !to be sure not to use it again
c          endif
c        enddo
c      endif
c      if(knucnt(2,k).gt.0)then
c        do nuc=1,knucnt(2,k)
c          if(npnuc(nuc,2,k).eq.n)then
c            itt=irnuc(nuc,2,k)
c            xmt(itt)=xmt(itt)+xxnuc(nuc,2,k)
c            if(xmt(itt)-1d0.ge.-1d-15)iet(itt)=0
c            xmpr(n,k)=xmpr(n,k)-xxnuc(nuc,2,k)
c            xpr(n,k)=xppr(n,k)*xmpr(n,k)
c            ypr(n,k)=0.5D0*log(xppr(n,k)/xmpr(n,k))
c            npnuc(nuc,2,k)=0    !to be sure not to use it again
c          endif
c        enddo
c      endif

      if(id.lt.0)then
        if(abs(id).eq.4)then         !Mass-
          xmt(it)=xmt(it)+xmpr(n,k) !xm restored to remnant
          xpt(it)=xpt(it)+xppr(n,k) !M2 distribution for remnant mass
        elseif(abs(id).eq.5)then     !Mass+
          xpp(ip)=xpp(ip)+xppr(n,k)    !xp restored to remnant
          xmp(ip)=xmp(ip)+xmpr(n,k) !M2 distribution for remnant mass
        else
c transfer momentum to remnant mass to keep energy conservation
          xpt(it)=xpt(it)+xppr(n,k)
          xmp(ip)=xmp(ip)+xmpr(n,k)
        endif
c transfer pt to opposite remnant to keep momentum conservation
        xxt(it)=xxt(it)+xxp1pr(n,k)+xxp2pr(n,k)
        xyt(it)=xyt(it)+xyp1pr(n,k)+xyp2pr(n,k)
        xxp(ip)=xxp(ip)+xxm1pr(n,k)+xxm2pr(n,k)
        xyp(ip)=xyp(ip)+xym1pr(n,k)+xym2pr(n,k)
      else
        xpp(ip)=xpp(ip)+xppr(n,k)
        xmt(it)=xmt(it)+xmpr(n,k)
        xxp(ip)=xxp(ip)+xxp1pr(n,k)+xxp2pr(n,k)
        xyp(ip)=xyp(ip)+xyp1pr(n,k)+xyp2pr(n,k)
        xxt(it)=xxt(it)+xxm1pr(n,k)+xxm2pr(n,k)
        xyt(it)=xyt(it)+xym1pr(n,k)+xym2pr(n,k)
      endif


cc      if(abs(id).eq.1.and.xpr(n,k).lt.xcupom)then
c      if(abs(id).eq.1.and.xpr(n,k).lt.xcupom**2/s)then
cc     &   .and.xpr(n,k)*s.lt.xmindiff*(amproj**2+amtarg**2))then
cc low mass excitation of connected remnant to take into account missing Pomeron in earliy stage (small x)
c        iept1=40
c        iept2=40
c        if(iclpro.ne.2)iept1=20    !no low mass excitation for meson projectile
c      else
cc high mass excitation for late Pomeron suppression 
c        iept1=20
c        iept2=20
c      endif
c      if(idfpr(n,k).eq.2)then
c        if(iep(ip).le.10.or.iep(ip)/10.eq.4)iep(ip)=iept1
c      elseif(idfpr(n,k).eq.3)then
c        if(iet(it).le.10.or.iet(it)/10.eq.4)iet(it)=iept2
c      else
c        if(abs(id).eq.1.and.idfpr(n,k).eq.4)then
c          iept1=iept1+2
c          iept2=iept2+2
c        endif
c        if(iep(ip).le.10.or.iep(ip)/10.eq.4)iep(ip)=iept1
c        if(iet(it).le.10.or.iet(it)/10.eq.4)iet(it)=iept2
c      endif


      if(nprt(k).eq.0)then      !no more Pomeron on this pair
        if(id.eq.4)then      !if last Pom was CD, suppress collision
          id=-1
        else
          if(abs(itpr(k)).eq.1)then !no more Pomeron on this pair
            itpr(k)=-2          !only mass
c          if(iep(ip).le.0)iep(ip)=40
c          if(iet(it).le.0)iet(it)=40
          elseif(iomega.eq.2.and.idpr(n,k).ne.7)then
            ncol=ncol-1
          endif
        endif
      endif

      istring=idp1pr(n,k)+idp2pr(n,k)+idm1pr(n,k)+idm2pr(n,k)
      if(id.ne.-1.and.istring.ne.0.and.iremn.ge.2)then
        if(ish.ge.7)write(ifch,*)"restore flavor:",istring

        if(idp1pr(n,k).eq.2)ivp(ip)=ivp(ip)+1 !update number of valence quark
        if(idm1pr(n,k).eq.2)ivt(it)=ivt(it)+1
        if(idp2pr(n,k).eq.2)ivp(ip)=ivp(ip)+1
        if(idm2pr(n,k).eq.2)ivt(it)=ivt(it)+1
        if(idp1pr(n,k).eq.5)then !update number of valence diquark (and quark)
          idp(ip)=idp(ip)+1     
          ivp(ip)=ivp(ip)+2     
        endif
        if(idm1pr(n,k).eq.5)then
          idt(it)=idt(it)+1
          ivt(it)=ivt(it)+2
        endif
        if(idp2pr(n,k).eq.5)then
          idp(ip)=idp(ip)+1
          ivp(ip)=ivp(ip)+2
        endif
        if(idm2pr(n,k).eq.5)then
          idt(it)=idt(it)+1
          ivt(it)=ivt(it)+2
        endif

        if(iremn.eq.3)then      !virtual Pomeron (remove unnecessary flavors for string ends)
          do j=1,2
            do i=1,nrflav
              jcp(i,j)=jcpref(i,j,ip)
              jct(i,j)=jctref(i,j,it)
            enddo
            do i=nrflav+1,nflav
              jcp(i,j)=0
              jct(i,j)=0
            enddo
          enddo
          if(ish.ge.7)write(ifch,*)"in:",jcp,' |',jct
          iret=0

c Projectile diquark-antidiquark pair
          iaq=nint(1.5+sign(0.5,float(idproj)))
          iq=3-iaq
          if(idp1pr(n,k).eq.4)then  !diquark
c    first quark
            idum=idrafl(iclpro,jcp,iaq,'v',0,iret)      !pick anti-quark
            ntry=0
            do while (jcp(idum,iq).eq.0.and.ntry.lt.100)!look for the corresponding quark
              ntry=ntry+1
              idum=idrafl(iclpro,jcp,iaq,'v',0,iret)
            enddo
            if(ntry.lt.100)then          !if OK, then remove the pair and pick a second quark
              call idsufl3(idum,1,jcp)
              call idsufl3(idum,2,jcp)
              if(jcp(idum,1)-jcpval(idum,1,ip).lt.0) !check valence quark number
     &             jcpval(idum,1,ip)=jcpval(idum,1,ip)-1
              if(jcp(idum,2)-jcpval(idum,2,ip).lt.0)
     &             jcpval(idum,2,ip)=jcpval(idum,2,ip)-1

c   second quark
              idum=idrafl(iclpro,jcp,iaq,'v',0,iret)
              ntry2=0
              do while (jcp(idum,iq).eq.0.and.ntry2.lt.100)!look for the corresponding antiquark
                ntry2=ntry2+1
                idum=idrafl(iclpro,jcp,iaq,'v',0,iret)
              enddo
              if(ntry2.lt.100)then          !if OK, then remove the pair
                call idsufl3(idum,1,jcp)
                call idsufl3(idum,2,jcp)
                if(jcp(idum,1)-jcpval(idum,1,ip).lt.0)
     &               jcpval(idum,1,ip)=jcpval(idum,1,ip)-1
                if(jcp(idum,2)-jcpval(idum,2,ip).lt.0)
     &               jcpval(idum,2,ip)=jcpval(idum,2,ip)-1
              else          !if not (because quarks already used by other valid string), then redo event to avoid problem in flavor conservation
                if(id.ge.15)then
                  id=-1
                  goto 1000
                else
                  call utstop("Virpom:should not happen (2) !&")
                endif
              endif
            else      !if no pair has been found (because quarks already used by other valid string), then redo event to avoid problem in flavor conservation
              if(id.ge.15)then
                id=-1
                goto 1000
              else
                call utstop("Virpom:should not happen  (3) !&")
              endif
            endif

c Projectile quark-antiquark pair
          else
            idum=idrafl(iclpro,jcp,iaq,'v',0,iret)      !pick anti-quark
            ntry=0
            do while (jcp(idum,iq).eq.0.and.ntry.lt.100)  !look for the corresponding quark
              ntry=ntry+1
              idum=idrafl(iclpro,jcp,iaq,'v',0,iret)
            enddo
            if(ntry.lt.100)then          !if OK, then remove the pair
              call idsufl3(idum,1,jcp)
              call idsufl3(idum,2,jcp)
              if(jcp(idum,1)-jcpval(idum,1,ip).lt.0)
     &             jcpval(idum,1,ip)=jcpval(idum,1,ip)-1
              if(jcp(idum,2)-jcpval(idum,2,ip).lt.0)
     &             jcpval(idum,2,ip)=jcpval(idum,2,ip)-1
            else                         !if not (because quarks already used by other valid string),then redo event to avoid problem in flavor conservation
              if(id.ge.15)then
                id=-1
                goto 1000
              else
                call utstop("Virpom:should not happen (4) !&")
              endif
            endif
          endif

c Target diquark-antidiquark pair
          iaq=nint(1.5+sign(0.5,float(idtarg)))
          iq=3-iaq
          if(idm1pr(n,k).eq.4)then  !diquark
c    first quark
            idum=idrafl(icltar,jct,iaq,'v',0,iret)
            ntry=0
            do while (jct(idum,iq).eq.0.and.ntry.lt.100)
              ntry=ntry+1
              idum=idrafl(icltar,jct,iaq,'v',0,iret)
            enddo
            if(ntry.lt.100)then
              call idsufl3(idum,1,jct)
              call idsufl3(idum,2,jct)
              if(jct(idum,1)-jctval(idum,1,it).lt.0)
     &             jctval(idum,1,it)=jctval(idum,1,it)-1
              if(jct(idum,2)-jctval(idum,2,it).lt.0)
     &             jctval(idum,2,it)=jctval(idum,2,it)-1
c    second quark
              idum=idrafl(icltar,jct,1,'v',0,iret)
              ntry2=0
              do while (jct(idum,2).eq.0.and.ntry2.lt.100)
                ntry2=ntry2+1
                idum=idrafl(icltar,jct,1,'v',0,iret)
              enddo
              if(ntry2.lt.100)then
                call idsufl3(idum,1,jct)
                call idsufl3(idum,2,jct)
                if(jct(idum,1)-jctval(idum,1,it).lt.0)
     &               jctval(idum,1,it)=jctval(idum,1,it)-1
                if(jct(idum,2)-jctval(idum,2,it).lt.0)
     &               jctval(idum,2,it)=jctval(idum,2,it)-1
              else
                if(id.ge.15)then
                  id=-1
                  goto 1000
                else
                  call utstop("Virpom:should not happen (5) !&")
                endif
              endif
            else
              if(id.ge.15)then
                id=-1
                goto 1000
              else
                call utstop("Virpom:should not happen (6) !&")
              endif
            endif

c Target quark-antiquark pair
          else
            idum=idrafl(icltar,jct,1,'v',0,iret)
            ntry=0
            do while (jct(idum,2).eq.0.and.ntry.lt.100)
              ntry=ntry+1
              idum=idrafl(icltar,jct,1,'v',0,iret)
            enddo
            if(ntry.lt.100)then
              call idsufl3(idum,1,jct)
              call idsufl3(idum,2,jct)
              if(jct(idum,1)-jctval(idum,1,it).lt.0)
     &             jctval(idum,1,it)=jctval(idum,1,it)-1
              if(jct(idum,2)-jctval(idum,2,it).lt.0)
     &             jctval(idum,2,it)=jctval(idum,2,it)-1
            else
              if(id.ge.15)then
                id=-1
                goto 1000
              else
                call utstop("Virpom:should not happen (7) !&")
              endif
            endif
          endif

          if(ish.ge.7)write(ifch,*)"out:",jcp,' |',jct
          do j=1,2
            do i=1,nrflav
              jcpref(i,j,ip)=jcp(i,j)
              jctref(i,j,it)=jct(i,j)
            enddo
          enddo

        endif
      endif


      endif

      ivpr(n,k)=0
c do not reset values but use ivpr=0 to have possibility of later check unless it is a backup Pomeron
      if(nvpr(n,k).ne.0)then
        nbkpr(n,k)=0
        nvpr(n,k)=0
        idpr(n,k)=0
        idfpr(n,k)=0
        xpr(n,k)=0d0
        ypr(n,k)=0d0
        xppr(n,k)=0d0
        xmpr(n,k)=0d0
        idp1pr(n,k)=0
        idp2pr(n,k)=0
        idm1pr(n,k)=0
        idm2pr(n,k)=0
        xm1pr(n,k)=0d0
        xp1pr(n,k)=0d0
        xm2pr(n,k)=0d0
        xp2pr(n,k)=0d0
        xxm1pr(n,k)=0d0
        xym1pr(n,k)=0d0
        xxp1pr(n,k)=0d0
        xyp1pr(n,k)=0d0
        xxm2pr(n,k)=0d0
        xym2pr(n,k)=0d0
        xxp2pr(n,k)=0d0
        xyp2pr(n,k)=0d0
        do jj=1,mxsplit
        call xpprbor_set(jj,n,k, 0. )
        call xmprbor_set(jj,n,k, 0. )
        call ptprboo_set(1,jj,n,k, 0. )
        call ptprboo_set(2,jj,n,k, 0.  )
        call rapprboo_set(1,jj,n,k, 0. )
        call rapprboo_set(2,jj,n,k, 0. )
        call gbpom_set(1,jj,n,k, 0. )
        call gbpom_set(2,jj,n,k, 0. )
        enddo 
        nemispr(1,n,k)=0
        nemispr(2,n,k)=0
      endif

1000  continue
      if(ish.ge.5)write(ifch,*)"remnant out",xpp(ip),xmp(ip)
     &                                     ,xmt(it),xpt(it)

      call utprix('VirPom',ish,ishini,3)

      end
 

c-----------------------------------------------------------------------
      subroutine CalcZZ(ir,m)
c-----------------------------------------------------------------------
C Calculates sum of all Q2 for remnant m for proj (ir=1) or target (ir=-1)
c   writes it to zzremn(m, 1 or 2)
c in case of nuclear collisions, add Q2s of surrounding interactions
c-----------------------------------------------------------------------
#include "aaa.h"
#include "ems.h"
#include "par.h"

      if(ir.eq.1)then
        zz=0.
        if(lproj(m).ge.1)then
          do l=1,lproj(m)
            kpair=kproj(m,l)
            do n=1,nprmx(kpair)
              zpar=q2kmin(1,n,kpair)
              if(zpar.gt.q2nmin)zz=zz+zpar
              if(ish.ge.10)
     .        write(ifch,*)'zzremn:',ir,m,kpair,n,zz,zpar,zpar.gt.q2nmin
            enddo
c            it=itarg(kpair)
c            do lt=1,ltarg3(it)
c              kt=ktarg3(it,lt)
c              do nt=1,nprmx(kt)
c                if(kt.ne.kpair.and.
c     .                        idpr(nt,kt).ne.0.and.ivpr(nt,kt).eq.1)then
c                  zpar=q2kmin(1,nt,kt)
c                  zz=zz+zpar
c                endif
c              enddo
c            enddo
          enddo
        endif
        zzremn(m,1)=zz
      elseif(ir.eq.-1)then
        zz=0.
        if(ltarg(m).ge.1)then
          do l=1,ltarg(m)
            kpair=ktarg(m,l)
            do n=1,nprmx(kpair)
              zpar=q2kmin(2,n,kpair) !float(nprt(kpair))
              if(zpar.gt.q2nmin)zz=zz+zpar
              if(ish.ge.10)
     .        write(ifch,*)'zzremn:',ir,m,kpair,n,zz,zpar,zpar.gt.q2nmin
            enddo
c            ip=iproj(kpair)
c            do lp=1,lproj3(ip)
c              kp=kproj3(ip,lp)
c              do np=1,nprmx(kp)
c                if(kp.ne.kpair.and.
c     .                        idpr(np,kp).ne.0.and.ivpr(np,kp).eq.1)then
c                  zpar=q2kmin(2,np,kp)
c                  zz=zz+zpar
c                endif
c              enddo
c            enddo
          enddo
        endif
        zzremn(m,2)=zz
      else
        stop'CalcZZ: invalid option.          '
      endif
      if(ish.ge.6)write(ifch,*)'zzremn:',ir,m,zz
      end

c-----------------------------------------------------------------------
      subroutine WriteZZ

c-----------------------------------------------------------------------
c compute estimated energy density for each string to have pbreak as a 
c function of energy density and not string mass
c for each partons we sum all energy of Pomerons weighted by the distance
c to its mother Pomeron and a factor zipinc.
c obsolete since new tim for BG where pbreak is constant (TP20170215)
c-----------------------------------------------------------------------

#include "aaa.h"
#include "ems.h"
c      double precision xor,yor,om51p,omgamk,xh,xpm
      common/emsptl/nppr(npommx,kollmx),npproj(mamx),nptarg(mamx)
c      double precision plc,s
c      common/cems5/plc,s
     
c --- Calculate Z (written to zzremn) from number of Pomerons

      do ip=1,maproj
       call CalcZZ(1,ip)
      enddo
      do it=1,matarg
       call CalcZZ(-1,it)
      enddo

c      do k=1,koll
c        do n=1,nprmx(k)
c          if(idpr(n,k).ne.0.and.idpr(n,k).ne.2.and.ivpr(n,k).ne.0)then
c            i=nppr(n,k)
c            zpaptl(1,i)=0.
c            zpaptl(2,i)=0.
c          endif
c        enddo
c      enddo
      return

c      jmin=maproj+matarg+1
c      jmax=nptl
c      if( z o p i n c .gt.0.)then  ! z o p i n c   redefined
c        alpqi=1./(alppar+1.)
c        xh=1d0/plc
c        xpm=sqrt(xh)
cc first fill zpaptl(1) for each Pomeron
c        do k=1,koll
cc calculate Q2ref
c          q2pmin(1)=q2nmin  !soft contribution at Q20
c          q2pmin(2)=q2nmin  !soft contribution at Q20
c          do n=1,nprmx(k)
c            if(idpr(n,k).ne.0.and.ivpr(n,k).ne.0)then
c              gsoft=om51p(sngl(s*xh),xh,0d0,bk(k),0)
c     .             /omGamk(n,k,xpm,xpm,0,0) 
c              q2ref=max(q2nmin,2.*q2nmin*gsoft**alpqi)
c              i=nppr(n,k)
c              zpaptl(1,i)=max(0.,
c     .                    0.5*(q2kmin(1,n,k)+q2kmin(2,n,k))-q2ref)   !<q2s>
c              zpaptl(2,i)=0.
c            endif
c          enddo
c        enddo
cc then sum all <q2s> from other pomerons to get local density
c        dnorm=2.* z o p i n c   **2   !z o p i n c redefined
cc        snorm=1./sqrt(pi*dnorm) !normalization put in fra.f
c        do i=jmin,jmax
c          if(idptl(i).gt.0)then !usefull Pomerons
c            xor=xorptl(1,i)
c            yor=xorptl(2,i)
c            do j=jmin,jmax
c              if(idptl(j).gt.0)then !usefull Pomerons
c                d2=sngl((xorptl(1,j)-xor)**2+(xorptl(2,j)-yor)**2)
c                dadd=zpaptl(1,j)*exp(-d2/dnorm) !energy density due t       
c                zpaptl(2,i)=zpaptl(2,i)+dadd
cc      write(ifch,*)'WriZZ',i,j,d2,zpaptl(1,i),zpaptl(2,i),pptl(5,j),dadd
c              endif
c            enddo
c          endif
c        enddo
c      endif

      end

c-----------------------------------------------------------------------
      subroutine ProReM(ir,irem,iret)
c-----------------------------------------------------------------------
c propose remnant mass of remnant irem in case of proj (ir=1)
c or target (ir=-1)
c   (-> xmp, xpt)
c iret : input : if iret=10 force to give mass even if no more energy,
c        when input not 10 : output = error if 1
c 2 cases :
c * low mass diffractive diagram
c   single mass diagram from Metropolis have been removed 
c   Energy is taken only from the other side nucleon which are close enough
c   to form a pair even if that pair was not used for a collision.
c * at least one diagram with mass
c   mass is already given to remnant but in case it is too small mass is
c   increased
c The 2 masses are needed to get large rapidity gaps from SD with low mass
c but we hould avoid high mass remnants in diffractive interactions because
c it creates a too large spread of particles in rapidity which is in conflict 
c with LHCb energy flow measurement (trigger problem): mass should be shared 
c between different strings (like remnant+2 strings from a small Pomeron).
c-----------------------------------------------------------------------

#include "aaa.h"
#include "ems.h"
#include "sem.h"
      double precision rr,xxx,xmin,xmax,msmin,xmmin,xpt2rem,xtest0,xtmp
      double precision alp,alptr,xi,xii,eps,sx,xmin0,xtest(mamx),fxtest
      parameter(eps=1.d-15)
      common/cemsr5/alptr(0:1,0:5)
      double precision plc,s,p5sq,aremn,aremnex,drangen,xmw,msmin0
      common/cems5/plc,s
      integer icrmn(2),jc(nflav,2)
      logical cont,force,drop,excited
      character cremn*4
      dimension k2j(mamx)

      call utpri('ProReM',ish,ishini,5)

      ntrymx=50
      do j=1,2
        do i=1,nflav
          jc(i,j)=0
        enddo
      enddo

c uncomment the following two lines to force the excitation

ccc      force=.true.   !<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<
ccc      ntrymx=1       !<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<

c initial definitions

      ntry=0
      if(ir.eq.1)then
        cremn='targ'
        icl=icltar
        jrem=1
        jremo=2
        massm=lproj(irem)     !number of target nucleon linked to irem
        masso=0               !number of acctive target nucleon linked to irem
        do k=1,massm
          k2j(k)=itarg(kproj(irem,k))
          xme(k2j(k))=-1.1d0
          if(iet(k2j(k)).ge.0)then
            xme(k2j(k))=0.d0
            masso=masso+1
          endif
        enddo
        icrmn(1)=icremn(1,irem,jrem)
        if(icrmn(1).eq.999999)then    !more than 9 quark : use jcpref
          do j=1,2
            do i=1,nrflav
              jc(i,j)=jcpref(i,j,irem)
            enddo
          enddo
        else
          icrmn(2)=icremn(2,irem,jrem)
          call iddeco(icrmn,jc)
        endif
        amremn=amproj
        kolpt=kolp(irem)
           !idx=isign(iabs(idproj)/10*10+1,idproj)
           !call idmass(idx,amremn)
        iremo1=itarg(1)
        zz=zzremn(irem,1)
c        if(iep(irem).eq.3.or.iep(irem).eq.5)zz=zzremn(irem,2)
      else       !if(ir.eq.-1)then
        cremn='proj'
        icl=iclpro
        jrem=2
        jremo=1
        massm=ltarg(irem)  !number of projectile nucleon linked to irem
        masso=0               !number of acctive target nucleon linked to irem
        do k=1,massm
          k2j(k)=iproj(ktarg(irem,k))
          xme(k2j(k))=-1.1d0
          if(iep(k2j(k)).ge.0)then
            xme(k2j(k))=0.d0
            masso=masso+1
          endif
        enddo
        icrmn(1)=icremn(1,irem,jrem)
        if(icrmn(1).eq.999999)then    !more than 9 quark : use jctref
          do j=1,2
            do i=1,nrflav
              jc(i,j)=jctref(i,j,irem)
            enddo
          enddo
        else
          icrmn(2)=icremn(2,irem,jrem)
          call iddeco(icrmn,jc)
        endif
        amremn=amtarg
        kolpt=kolt(irem)
           !idx=isign(iabs(idtarg)/10*10+1,idtarg)
           !call idmass(idx,amremn)
        iremo1=iproj(1)
        zz=zzremn(irem,2)
c        if(iet(irem).eq.3.or.iet(irem).eq.5)zz=zzremn(irem,2)
      endif
      ieex=mod(iez(irem,jrem),10)
      iept=iez(irem,jrem)/10
      if(iret.eq.10)then
        force=.true.
        alplead1=0.
      elseif(iret.eq.5)then !.or.ieex.eq.3)then
        iret=0
        force=.false.
        alplead1=0.           !not always possible to respect xr**alplea
      else
        iret=0
        force=.false.
        alplead1=alplea(icl)+1.
      endif
      drop=.false.
      if(iremn.ge.2.and.iept.eq.5.and.irmdrop.eq.1)
     &   drop=.true.
      excited=.false.
      if(iept.gt.0.and.iept.ne.6)then
        excited=.true.
        alp=alptr(idz(irem,jrem),iept)
        if(iept.eq.3)alp=dble(alppom)
     .             +(alp-dble(alppom))*(1d0-xpz(irem,jrem))**zdfinc
c        if(iept.ge.3)alp=alp+(dble(alppom)-alp)
c     .   /(1d0-log(xpz(irem,jrem))*dble(q2nmin/(zz+q2nmin))*zdfinc)
c        if(iept.ge.1)alp=dble(alppom)+(alp-dble(alppom))
c     .   /(1d0-log(1d0-xpz(irem,jrem))*zdfinc)
c     .               /(1d0-log(xpz(irem,jrem))*zdfinc)
c     .               /(1d0-log(xpz(irem,jrem))*dble(q2nmin/zz)*zdfinc)
c     .              +(alp-dble(alppom))*(1d0-xpz(irem,jrem))**zdfinc
c        print *,iez(irem,jrem),alp,xpz(irem,jrem)
c        if(zz.gt.2.)then
c          alp=alphigh(1+idz(irem,jrem)) !/dble(1.+zmsinc*zz)    !high mass if more than 2 Pomerons attached to remnant
cc        xmax=1d0/plc
c        endif
      else
        alp=0d0
      endif

c for spectators only low mass and few partners, so do not care about energy
      if(iept.eq.6)force=.true.

c defs

      sx=s*xpz(irem,jrem)
      xpt2rem=xxz(irem,jrem)**2d0+xyz(irem,jrem)**2d0
      xmax=1d0


c  fremnux (+) and not fremnux2 (-) which gives a mass too low in case of getstring where q and aq do not cancel


      if(excited)then
        aremn=dble(max(amremn,fremnux(jc)))
c       if(iremn.eq.2.and.iept.eq.3)      !droplet
c     &     aremn=dble(max(amremn,fremnux(jc)))
        aremnex=aremn+amemn(idz(irem,jrem),iept)
c        if(iremn.ge.2)then
c          if(drop)then  !from getdropx
c            mamod=4
c            fad=alpdro(1)+zdrinc*zz
c            aremnex=dble(fad+utamnz(jc,mamod))
cc            fas=2.
cc            aremnex=max(aremnex,dble(fas*utamnz(jc,mamos)))
cc            aremnex=aremnex+zmsinc*(zz-q2nmin)  
c        endif
c        if(ieex.eq.3)aremnex=aremnex+zmsinc*zz 
      else    !minimum mass for spectators should be as low as possible
        aremn=dble(max(amremn,fremnux2(jc)))
        aremnex=aremn
      endif

      p5sq=xmz(irem,jrem)*sx
      msmin0=aremn**2d0+xpt2rem
      if(excited)then
        msmin=1.5d0*(aremnex**2d0+xpt2rem)
        msmin0=1.5d0*msmin0                  !otherwise cluster can not decay
      else
        msmin=msmin0
      endif

      if(ish.ge.8)write(ifch,10)ir,irem,masso,icrmn,iept,force
     &                  ,amremn,fremnux(jc),aremn,msmin,msmin0,p5sq
     &               ,xpz(irem,jrem),xmz(irem,jrem),xpt2rem,sx,zz
 10   format('prorem :  ',i3,2i4,2i7,i2,L2,/
     &      ,'    mass :',6g13.5,/
     &      ,' x,pt,sx,zz :',5g13.5)

      if(p5sq.gt.msmin)then   !nothing needed, mass OK
        xxx=0d0
        xmw=0d0
        masso=0
        xmin=msmin/sx
        goto 900
      else       !mass to low but take into account initial value of xmz
        xmw=xmz(irem,jrem)
      endif

c###########################################################################
c if mass from Metropolis is too low, use old process to distribute mass
c###########################################################################

    1 ntry=ntry+1
      if(ntry.gt.ntrymx)then
        if(ish.ge.5)then
          call utmsg('ProReM')
          write(ifch,*)'Remnant mass assignment not possible (ntry)'
     &                 ,ir,irem
          if(force)write(ifch,*)'Ignore p4 conservation'
          call utmsgf
        endif
        if(.not.force)then
          iret=1
          goto 1000
        else
c not enough energy availabe : force last mass and check
          goto 900
        endif
      endif

c check

      if(xpz(irem,jrem).le.0.d0)then
        write(ifch,*)'ProRem ipp',xpz(irem,jrem)
     &                           ,jrem,irem,lremn(irem,jrem)
        do li=1,lremn(irem,jrem)
          kkk=kremn(irem,li,jrem)
          write(ifch,*)'kkk',kkk
        enddo
        call XPrint('ProRem :&')
        call utstop('Big problem in ProRem !&')
      endif

c xtest = xminus-max,  corresponding mostly to a remnant mass 0.2

      xtest0=0.d0
      fxtest=1d0 !0.4d0*(1d0+drangen(xxx)) !1.d0 !0.3d0
      do k=1,massm
        j=k2j(k)
        if(xme(j).ge.0d0)then
          cont=.false.
ctp        if(xmz(j,jremo).gt.eps.and.iez(j,jrem).gt.0)then !xmz(,jremo)=xplus
ctp060824        if(xmz(j,jremo).gt.eps.and.iez(j,jrem).ge.0)then !xmz(,jremo)=xplus
c        if(iez(j,jremo).gt.0.or.koll.eq.1)then !xmz(,jremo)=xplus
          if(xpz(j,jremo)*xmz(j,jremo).ge.xzos(j,jremo))then !xmz(,jremo)=xplus
c mass of opposite side remnant is already properly defined so keep mass above minimum mass
            cont=.true.
            xmmin=xzos(j,jremo)/xmz(j,jremo)
          else
c xmz will be updated later, so use maximum allowed value (xp of current remnant)
            xmmin=xzos(j,jremo)/xpz(irem,jrem)
          endif
          xtest(j)=xmw+xpz(j,jremo)-xmmin !maximal momentum available
!this term is very important for non excited remnants in pp, it changes the xf
! distribution of proton and the multiplicity at low energy. Fxtest should not
! be to close to 0. otherwise it makes a step in xf distribution of p at
! 1-fxtest but if fxtest=1, multiplicity at low energy is too high ...
! but better (and smoother) with exponential decrease).
          if(.not.cont)then
            if(xtest(j).gt.0d0)then
              xtest(j)=min(xtest(j),fxtest/xpz(irem,jrem))
            else
              xtest(j)=min(1.d0,fxtest/xpz(irem,jrem))
            endif
          endif
c        else
c          xtest(j)=0.01d0       !maximal momentum available for non exited state
c        endif
         xtest0=max(xtest0,xtest(j))
c        print *,iep(1),iet(1),iept,xtest(j),xpz(j,jremo),xmmin
c        & ,xzos(j,jremo),xmz(j,jremo)
       endif
      enddo
ctp060824      if(.not.cont)xtest=min(1.d0,0.2d0/xpz(irem,jrem))


      cont=.true.

c determine xminus

c      xmin0=1.05*(aremn**2d0+xxz(irem,jrem)**2d0+xyz(irem,jrem)**2d0)/sx
c      xmin=1.1*(aremnex**2d0+xxz(irem,jrem)**2d0+xyz(irem,jrem)**2d0)/sx
      xmin0=max(xmw,msmin0/sx)

      if(iept.eq.0)then !no excitation
        xmin=xmin0+eps
      elseif(iept.eq.2)then !pion exchange, minim should not change
        xmin=dble(xmindiff)*(aremnex**2d0+xpt2rem)/sx
      elseif(iept.eq.3.or.iept.eq.5)then
        xmin=dble(xminremn)*(aremnex**2d0+xpt2rem)/sx
      else
        xmin=(aremnex**2d0+xpt2rem)/sx
      endif
      xmin=max(xmin0+eps,xmin)

c      if(ieex.eq.3)xmax=min(xmax,
c     .                           zz*zdfinc/sx+xmin)
c     .     q2nmin*exp(min(50.,(zz-q2nmin)/(max(1.e-5,zdfinc))))/sx+xmin)
      xmax=min(xmax,xtest0)
c for diffractive remnant, mass should never exceed 5% of the proj or targ energy
c      if(iept.eq.1)then
c        xmax=min(xmax,max(dble(xminremn),xmin))
c      elseif(iept.eq.2)then
c        xmax=min(xmax,max(dble(xmindiff),xmin))
c      endif
c      if(iept.eq.1.or.iept.eq.3)then
c       xtmp=max(dble(min(1.,xminremn*float(maproj+matarg-1))),xmin)
c     &               *drangen(xmin)
      xtmp=1.d0
      if(excited)then
c        xtmp=min(1d0,dble(xmxrem(jremo))*dble(masso)
c     &                      *drangen(xmin)**0.05)
      if(iept.eq.2)then
cc        xtmp=max(min(1d0,dble(xmindiff)),xmin)!*drangen(xmin)
        xtmp=min(1d0,dble(xmxrem(jremo))*dble(masso)
     &                      *drangen(xmin)**0.05)
cc        xtmp=dble(xmindiff)
c      elseif(iept.eq.1)then
c        xtmp=min(1d0,dble(xmxrem(jremo))*dble(masso)
c     &                      *drangen(xmin)**0.05)
cc        xtmp=dble(xminremn)
c      elseif(drop)then     !3 or 5
cc       xtmp=max(dble(min(1.,xminremn*float(maproj+matarg-1))),xmin)
cc     &               *drangen(xmin)
c        xtmp=min(1d0,dble(xmxrem(jremo))*(zz/q2nmin)*dble(masso)
c     &                         *drangen(xmin)**0.05)
cc        xtmp=dble(xminremn)
      endif
      endif
      xmax=min(xmax,max(xtmp,xmin))
      if(koll.eq.1)xmax=min(xmax,xpz(iremo1,jremo))
c      print *,xmin0,xmin,xmax,xmxrem(jrem),iept,nrevt
      if(ish.ge.8)write(ifch,*)'ntry',ntry,xmin0,xmin,xmax,xtmp
     *                               ,xmax*dble(masso)
      if(xmin.ge.min(1d0,xmax*dble(masso)-eps))then
        xmax=xmax*dble(masso)
        xmin=xmin0
        if(xmin0.ge.min(1d0,xmax-eps))then
          if(.not.force)then
            iret=1
c if forced forget energy conservation
          elseif(excited)then
            xmz(irem,jrem)=min(1d0,xmin0*(1d0+drangen(xmin))) !random not to form a peak
          else
            xxx=min(1d0,(aremn**2d0+xpt2rem)/sx)
            xmz(irem,jrem)=xxx
            xmin0=max(drangen(xxx),(1d0-2d0*drangen(xxx)*sx**(-0.33)))
     &           *xxx
c     &         *exp(-1d0*drangen(xxx)**2)
c            xmin0=(1d0-exp(-4d0*drangen(xxx)))
c     &           *(1d0-((amzmn(idz(irem,jremo),jremo)
c     &             +sqrt(2d0*engy*drangen(xxx))
c     &                                       )**2+xpt2rem)/sx)*xxx
          endif
          goto 1000
        endif
      elseif(xmin.ge.xmax)then
        xmax=1d0
        xmin=min(1d0,xmin)
      endif
      rr=drangen(xmax)
      xxx=0.d0
      if(excited)then
c        xmin=xmin-xpt2rem/sx                     !no pt
c        xmax=xmax-xpt2rem/sx                     !no pt
c        alp=alptr(idz(irem,jrem),iept)/dble(1.+zmsinc*zz)

        call remn1paramget6(f1,f2,f3,dmy4,dmy5,dmy6)
        Z=max(0.,ng1evt-2.)/f3
        remnMassMin=f1+Z*(f2-f1)
        xremnMin=(remnMassMin/engy)**2
        xminxx=max(xmin,xremnMin)

        if(dabs(alp-1.d0).lt.eps)then
          xxx=xmax**rr*xminxx**(1d0-rr)
        else
          xxx=(rr*xmax**(1d0-alp)+(1d0-rr)*xminxx**(1d0-alp))
     &                                             **(1d0/(1d0-alp))
        endif
        !write(ifmt,'(a,3f10.2)')'TESTremn1',sqrt(xminxx)*engy
        !.  ,sqrt(xxx)*engy
c        xxx=xxx+xpt2rem/sx                       !no pt
!smooth distribution
c        if(iept.eq.4)xmin=xmin0
        if(iept.eq.4)then
          xmin0=max(drangen(xxx),(1d0-2d0*drangen(xxx)*sx**(-0.33)))*xxx
c     &         *exp(-1d0*drangen(xxx)**2)
        else
          xmin0=xmin+(1d0-exp(-drangen(xxx)**3))*(xxx-xmin)
        endif
      else
        if(masso.eq.1)ntry=ntrymx   !xxx is fixed so 1 try is enough
c        xmin=dble(amremn)**2d0/sx                !no pt
c        xxx=xmin+xpt2rem/sx                      !no pt
        xxx=xmin
        if(xmin.gt.xmax+eps)then
          if(ish.ge.6)write(ifch,*)'xmin>xmax for proj not possible (2)'
     &                 ,ir,irem
          if(.not.force)then
            iret=1
          else
            xmz(irem,jrem)=min(1d0,xxx)
          endif
          goto 1000
        endif
c to have a nice diffractive peak, do not allow too much fluctuation
c (xmin0 not too small compared to xxx)
c this function is more or less a fit of the diffractive peak
c (pp100, ep-forward (ZEUS), NA49, pipp100, taking into account the
c contribution of inelastic remnants)
c relative width of the peak reduces with increasing energy 
        xmin0=max(drangen(xxx),(1d0-2d0*drangen(xxx)*sx**(-0.33)))*xxx
c     &         *exp(-1d0*drangen(xxx)**2)
c      write(*,*)'->',sx**(-0.33),(1d0-2d0*drangen(xxx)*sx**(-0.33))

c        xmin0=dble(0.9+0.09*rangen())*xxx
c         xmin0=0.8*xxx
      endif
c      print *,
      if(ish.ge.6)write(ifch,*)
     &  iept,'alp',alp,xmin,xxx,xmax,xmw!,zz,1.+zmsinc*zz
      msmin=xmin*sx

c     partition xminus between nucleons of the other side
c if not enough energy available, better to reduce the mass than violate energy conservation
      if(force)alplead1=0.
      xxx=xxx-xmw               !take into account momentum already there
      if(xxx.lt.-eps)call utstop('ProRem:xxx<0 !!!&')
      do k=1,massm
        xme(k2j(k))=0.d0
      enddo
      xii=1d0
      ii=massm
      kk=int(rangen()*float(ii))+1   ! choose ramdomly a nucleon to start

      do while(ii.gt.0)

        iro=k2j(kk)
        cont=iez(iro,jremo).lt.0.or.xme(iro).lt.0d0
        if(ish.ge.8)
     .  write(ifch,*)'cont',ii,kk,iro,iez(iro,jremo),xme(iro),cont
        do while(cont)
          kk=kk+1
          if(kk.gt.massm)kk=kk-massm
          iro=k2j(kk)
          ii=ii-1
          cont=iez(iro,jremo).lt.0.or.xme(iro).lt.0d0
          if(ish.ge.8)
     .    write(ifch,*)'cont',ii,kk,iro,iez(iro,jremo),xme(iro),cont
          if(ii.lt.1)goto 1
        enddo

        if(ii-1.gt.0)then
          xi=xii*drangen(xii)**(2d0/dble(ii-1))
        else
          xi=0d0
        endif
        xtmp=xxx*(xii-xi)
        xme(iro)=abs(xme(iro))+xtmp

        xmmin=xzos(iro,jremo)
        if(xpz(iro,jremo)*xmz(iro,jremo).gt.xmmin)then   !mass already properly defined
          xmmin=xmmin/xmz(iro,jremo)
        elseif(koll.eq.1.and.xtest(iro).gt.eps)then
          xmmin=xmmin/min(xpz(irem,jrem),xtest(iro))
        elseif(xtest(iro).gt.eps)then
          xmmin=xmmin/xtest(iro)
        endif
        fxtest=(1.d0-xme(iro)/xpz(iro,jremo))**alplead1
        if(drangen(xmmin).gt.fxtest)then
          if(ish.ge.7)write(ifch,*)'alpha skip ',cremm,fxtest
          if(fxtest.le.1d-9)goto 1
          xme(iro)=xme(iro)-xtmp
        else
          if((xpz(iro,jremo)-xme(iro)).lt.xmmin)then
c           use some part of the remaining mass anyway
            xi=xtmp                   !mass we wanted to give
            xme(iro)=xme(iro)-xtmp    !restore previous xme
            xtmp=xtmp/xxx*(xpz(iro,jremo)-xmmin-xme(iro))  !mass we can give
            if(ish.ge.7)write(ifch,*)'    skip ',cremn,' ',ii,massm,ntry
     &         ,iro,xme(iro)+xi,xtmp/xxx,xpz(iro,jremo)-xme(iro),xmmin
            if(xtmp.gt.1d-4*xxx)then
              xme(iro)=xme(iro)+xtmp
              xii=xii-xtmp/xxx  !new fraction of xxx already distributed
            else
              xme(iro)=-xme(iro)
              if(abs(xme(iro)).lt.eps)xme(iro)=-1.1d0
              if(ii.le.1)goto 1
            endif
          else
            if((xpz(iro,jremo)-xme(iro)+xtmp-xxx*xii).gt.xmmin
     .               .and.xii.le.1d0/dble(masso))then !take remaining mass
              xme(iro)=xme(iro)-xtmp+xxx*xii
              xii=0d0
              ii=0
            else
              xii=xi
            endif
            if(ish.ge.7)write(ifch,*)'      ok ',cremn,' ',ii,massm,ntry
     &                  ,iro,xme(iro),xme(iro)/xxx,xpz(iro,jremo),xmmin
          endif
        endif
        kk=kk+1
        if(kk.gt.massm)then
          kk=kk-massm
          if(ii.gt.0)ii=massm   !to allow to pass more than once on each nucleon
        endif

      enddo

c check xmz(irem,jrem)

 900  xmz(irem,jrem)=min(1d0,xmz(irem,jrem)+xxx)

      !write(ifmt,'(a,3e10.2)')'TESTremn1',xpz(irem,jrem),xmz(irem,jrem)
      !.,sqrt(xpz(irem,jrem)*xmz(irem,jrem))*engy
      !.,0.5*log(xpz(irem,jrem)/xmz(irem,jrem))

      p5sq=xpz(irem,jrem)*xmz(irem,jrem)*s
      if(ish.ge.7)write(ifch,*)'final mass',irem,sqrt(p5sq),sqrt(msmin)
     &,xpz(irem,jrem),xmz(irem,jrem),force
      if(p5sq-msmin.lt.-1d-10)then
        if(ish.ge.5)then
          call utmsg('ProReM')
          write(ifch,*)'Remnant mass assignment not possible (M<Mmin)!'
     &                 ,ir,irem
          if(force)write(ifch,*)'Ignore p4 conservation'
          call utmsgf
        endif
        if(.not.force)then
          iret=1
        elseif(xpz(irem,jrem).gt.0.d0)then
          xmz(irem,jrem)=min(1d0,xmin*(1d0+drangen(xmin)))   !random not to form a peak
        endif
        goto 1000
      endif

c subtract xme

      do k=1,massm
        iro=k2j(k)
        if(xme(iro).gt.-1.d0)then
          xpz(iro,jremo)=xpz(iro,jremo)-abs(xme(iro))  !xpz(,jremo)=xminus
c          if(xpz(iro,jremo).lt.0d0)then
c            write(ifch,*)'ici',k,iro,xpz(iro,jremo),xme(iro)
c            stop
c          endif
        endif
      enddo

 1000 continue
      if(iret.ne.1)then
        xzos(irem,jrem)=xmin0*xpz(irem,jrem)
        if(ish.ge.8)write(ifch,*)'xzos',irem,jrem,xzos(irem,jrem)
     &             ,xpz(irem,jrem),xmin0
      endif

      call utprix('ProReM',ish,ishini,5)

      end

c-----------------------------------------------------------------------
      subroutine ProReF(ir,m,iretxx)
c-----------------------------------------------------------------------
c  proposes flavor for remnant m for proj (ir=1) or target (ir=-1)
c  and writes remnant into /cptl/ as string or hadron
c   ityptl definitions:
c      50  40  ...  rmn drop no split
c      51  41  ...  rmn drop from split
c      52  42  ...  rmn str inel
c      53  43  ...  rmn str diff
c      54  44  ...  rmn str inel from split
c      55  45  ...  rmn res
c      56  46  ...  rmn res from split without connexion
c      57  47  ...  rmn res active spectators
c      58  48  ...  rmn res from diff
c      59  49  ...  rmn res from string split
c-----------------------------------------------------------------------

#include "aaa.h"
#include "ems.h"
#include "sem.h"

      double precision plc,s   ,ptt1,ptt2,ptt3!,ax,ay,az
      common/cems5/plc,s
      double precision tpro,zpro,ttar,ztar,ttaus,detap,detat
      common/cttaus/tpro,zpro,ttar,ztar,ttaus,detap,detat
      common /cncl/xproj(mamx),yproj(mamx),zproj(mamx)
     *            ,xtarg(mamx),ytarg(mamx),ztarg(mamx)
      double precision amasmin,amasini,xmdrmax,xmdrmin,utdble!,utpcmd
      integer icf(2),icb(2)
      integer jcf(nflav,2),jcval(nflav,2),jcfstr(nflav,2)
      logical gdrop,gproj
      double precision ept(5),ep(4),aa(5),am2t,piq1,piq2,piq3
      common/emsptl/nppr(npommx,kollmx),npproj(mamx),nptarg(mamx)
      common /ems12/iodiba,bidiba  ! defaut iodiba=0. if iodiba=1, study H-Dibaryon
      character c*1,c1*1,c2*1

      integer getAccumJerr

      call utpri('ProReF',ish,ishini,3)

      iretxx=0

      if(ir.ne.1.and.ir.ne.-1)stop'ProReF: wrong ir'

      irmdropx=irmdrop
 55   idrop=1
      gdrop=.false.
      iret=0
      dens=0.0765
      do j=1,2
        do i=1,nflav
          jcf(i,j)=0
        enddo
      enddo

      if(ir.eq.1)then
c        if(kolp(m).le.0)goto1000
        if(iep(m).le.-1)goto 1000
        gproj=.true.
        mm=npproj(m)
        iept=iep(m)/10
        ieex=mod(iep(m),10)
        zz=zzremn(m,1)
        iclpt=iclpro
        call idspin(idptl(m),isopt,jspin,istra)
        if(iremn.ge.2)then         !number of valence quarks still in proj
          do nnn=1,nrflav
            jcval(nnn,1)=jcpval(nnn,1,m)
            jcval(nnn,2)=jcpval(nnn,2,m)
          enddo
          do nnn=nrflav+1,nflav
            jcval(nnn,1)=0
            jcval(nnn,2)=0
          enddo
        else
          do nnn=1,nflav
            jcval(nnn,1)=0
          enddo
          do nnn=1,nflav
            jcval(nnn,2)=0
          enddo
        endif
      elseif(ir.eq.-1)then
c        if(kolt(m).le.0)goto1000
        if(iet(m).le.-1)goto 1000
        gproj=.false.
        mm=nptarg(m)
        iept=iet(m)/10
        ieex=mod(iet(m),10)
        zz=zzremn(m,2)
        iclpt=icltar
        call idspin(idptl(maproj+m),isopt,jspin,istra)
        if(iremn.ge.2)then         !number of valence quarks still in proj
          do nnn=1,nrflav
            jcval(nnn,1)=jctval(nnn,1,m)
            jcval(nnn,2)=jctval(nnn,2,m)
          enddo
          do nnn=nrflav+1,nflav
            jcval(nnn,1)=0
            jcval(nnn,2)=0
          enddo
        else
          do nnn=1,nflav
            jcval(nnn,1)=0
          enddo
          do nnn=1,nflav
            jcval(nnn,2)=0
          enddo
        endif
      else
        gproj=.false.
        mm=0
        iept=0
        ieex=0
        zz=0.
        iclpt=0
        call utstop('ProReF: ir ???&')
      endif
      if(ish.ge.3)write(ifch,*)'remnant particle index:'
     &                         ,mm,m,iclpt,isopt,iept*10+ieex

      if(ish.ge.8)call alist('ProRef&',1,nptl)
      antotre=antotre+1

      mmini=mm
      nptlini=nptl
      minfra=min(minfra,nptlini)   !for trigger condition

      do l=1,5
       ept(l)=utdble(pptl(l,mm))
      enddo

      ifrptl(1,mm)=0
      ifrptl(2,mm)=0

c  initialize forward and backward ic (to transform remnant into string)

      if(gproj)then
        icf(1)=icproj(1,m)
        icf(2)=icproj(2,m)
        if(icf(1).eq.999999)then    !more than 9 quark : use jcpref
          do j=1,2
            do i=1,nrflav
              jcf(i,j)=jcpref(i,j,m)
            enddo
          enddo
        else
          call iddeco(icf,jcf)
        endif
      else                     !gtarg
        icf(1)=ictarg(1,m)
        icf(2)=ictarg(2,m)
        if(icf(1).eq.999999)then    !more than 9 quark : use jctref
          do j=1,2
            do i=1,nrflav
              jcf(i,j)=jctref(i,j,m)
            enddo
          enddo
        else
          call iddeco(icf,jcf)
        endif
      endif
      icb(1)=0
      icb(2)=0

      call idquacjc(jcf,nqu,naq)
c use RemoveHadron if too many c quarks
      if(nrflav.gt.3)then
        nqc=jcf(4,1)+jcf(4,2)
        if(nqu.lt.3.and.jcf(4,1).gt.1.or.
     &     naq.lt.3.and.jcf(4,2).gt.1.or.
     &             jcf(4,1)*jcf(4,2).gt.1 )nqc=4
      else
        nqc=0
      endif
c      if(iremn.ge.2)then
c        ier=0
c        ires=0
c        id=idtra(icf,ier,ires,0)
c        if(ier.eq.0)then
c          call idspin(id,ispin,jspin,istra)
c        else
c          ispin=0
c          jspin=0
c          istra=0
c        endif
c      endif

c define masses
      p2t=ept(1)**2+ept(2)**2
      amasmin=dble(fremnux(jcf))**2.d0+p2t
      if(ept(5).le.0.d0)then
        ept(5)=dble(fremnux(jcf)*(1.+rangen()))
        if(ish.ge.1)then
          call utmsg('ProReF')
          write(ifch,*)'zero remnant mass -> amasmin'
          call utmsgf
        endif
      endif
      am2t=sqrt(p2t+ept(5)**2)
      if(ept(4).gt.am2t.and.(iept.eq.0.or.iept.eq.6))then
        ept(3)=sign(sqrt((ept(4)+am2t)*(ept(4)-am2t)),ept(3))
      else
        ept(4)=sqrt(ept(3)*ept(3)+ept(2)*ept(2)+ept(1)*ept(1)
     &           +ept(5)*ept(5))
      endif
      am2t=(ept(4)+ept(3))*(ept(4)-ept(3))-(ept(1)**2+ept(2)**2)
      if(ish.ge.2
     &   .and.(am2t.lt.-1d0.or.abs(am2t-ept(5)*ept(5)).gt.ept(5)))then
          write(ifch,*)'Precision problem in ProRef, p:',
     &             (ept(k),k=1,4),ept(5)*ept(5),am2t
      endif
      ept(4)=sqrt(ept(3)*ept(3)+ept(2)*ept(2)+ept(1)*ept(1)
     &           +ept(5)*ept(5))

      if(ish.ge.3)then
        if(gproj)then
            write(ifch,'(a,5e11.3,2i7)')' proj:'
     &      ,(sngl(ept(k)) ,k=1,5),(icproj(k,m) ,k=1,2)
        else    !gtarg
           write(ifch,'(a,5e11.3,2i7)')' targ:'
     &      ,(sngl(ept(k)) ,k=1,5),(ictarg(k,m),k=1,2)
         endif
      endif

      amasini=ept(5)*ept(5)+p2t

      xmdrmin=dble(alpdro(1)+fremnux(jcf)+sqrt(zz)*zdrinc)**2
      xmdrmax=dble(fremnux(jcf)+amdrmax)**2


      if(ish.ge.4)write(ifch,*)'remnant masses:',am2t,amasini,amasmin
     &                ,xmdrmin,xmdrmax,zz,nqu,naq,nqc,ieex,iept

      call remn1paramget6(dmy1,dmy2,f3,f4,f5,dmy6)
      Z=max(0.,ng1evt-2.)/f3
      forceDropletMass=f4+Z*(f5-f4) !force droplet decay beyond this mass

c.............................exotic ...................................
      
      if(irmdropx.gt.0.and.(   
     &    (ept(5).gt.forceDropletMass) 
     &    .or. (
     &           (irmdropx.eq.1.or.irmdropx.eq.2).and.
     &           ( ieex.eq.-3.or.iept.eq.5.or.
     &             (.not.((nqu.eq.3.and.naq.eq.0)
     &                    .or.(nqu.eq.0.and.naq.eq.3)
     &                    .or.(nqu.eq.1.and.naq.eq.1))
     &             )
     &             .and.amasini.gt.amasmin.and.irmdropx.eq.1
     &           )
     &           .or.(amasini.gt.xmdrmax.and.irmdrop.eq.2.and.ieex.ne.2) 
     &         ) 
     &  ))then
        !print*,'+++++',iept,nqu,naq

c charm not possible in droplet

        if( ept(5).gt.forceDropletMass .or. 
     &    amasini.gt.xmdrmin.or.nqc.ne.0)then
          if(iremn.eq.2)then
            call getdropx(ir,iept,m,icf,jcf,jcval,zz,ept,aa
     &                                          ,gdrop,xmdrmax)
          else
            call getdroplet(ir,iept,icf,jcf,zz,ept,aa,gdrop,xmdrmax)
          endif
          !--------------------------------
          !emit a droplet, update the remnant string flavor and 5-momentum
          ! input
          !     ir ......... 1  projectile, -1  target remnant
          !     ept ........ remnant  5-momentum
          !     jcf ........ remnant jc
          ! output
          !     gdrop ...  .true. = successful droplet emission
          !                          jcf, ept ....... droplet  ic and 5-momentum
          !                          icf, a ......... remnant string jc and 5-momentum
          !               .false. = unsuccessful
          !                          jcf, ept .... unchanged,
          !                          emits hadrons instead of droplet
c         !                          considered as droplet jc and 5-momentum
          !-------------------------------------
        endif

c redefine energy and charm quarks in droplet
c        amasini=ept(5)*ept(5)
        nqc=jcf(4,1)+jcf(4,2)
c use remove hadrons if droplet too heavy (should not happen) or charm
c        if(amasini.gt.1e4.or.nqc.ne.0)goto 500
c        if(nqc.ne.0.and.irmdrop.eq.1)goto 500

        !...........droplet
        !also in case of unsuccessful drop emission, then remnant = droplet !
        nptl=nptl+1
        t=xorptl(4,mm)
        istptl(mm)=41
        ifrptl(1,mm)=nptl
        ifrptl(2,mm)=nptl
        tivptl(2,mm)=t
c            Remnant radius to have eps=dens GeV/fm3
        radptl(nptl)=(3.*sngl(ept(5))/4./pi/dens)**0.3333
        dezptl(nptl)=0.
        qsqptl(nptl)=0.
        zpaptl(1,nptl)=zpaptl(1,mm)
        zpaptl(2,nptl)=zpaptl(2,mm)
        do l=1,5
          pptl(l,nptl)=sngl(ept(l))
        enddo
        if(gdrop)then
          idx=0
        else
          idx=idtra(icf,0,0,0)
        endif
        if(idx.ne.0)then
         amx=sngl(ept(5))
         call idres(idx,amx,idrx,iadjx,1)
         idx=idrx
        endif
        if(idx.eq.0)then
          istptl(nptl)=10      !it is a droplet, so istptl=10 (important because in utghost and jintpo, this is used to identify this as a droplet)
          call idenct(jcf,idptl(nptl)
     *    ,ibptl(1,nptl),ibptl(2,nptl),ibptl(3,nptl),ibptl(4,nptl))
          if(gproj)then
            ityptl(nptl)=40
          else  !gtarg
            ityptl(nptl)=50
          endif
        else
          istptl(nptl)=0
          idptl(nptl)=idx
          pptl(5,nptl)=amx
          pptl(4,nptl)=sqrt(amx*amx+pptl(1,nptl)*pptl(1,nptl)
     &       +pptl(2,nptl)*pptl(2,nptl)+pptl(3,nptl)*pptl(3,nptl))
          if(gproj)then
            ityptl(nptl)=45
            if(iept.eq.6)ityptl(nptl)=47
          else  !gtarg
            ityptl(nptl)=55
            if(iept.eq.6)ityptl(nptl)=57
          endif
        endif
        iorptl(nptl)=mm
        jorptl(nptl)=0
        ifrptl(1,nptl)=0
        ifrptl(2,nptl)=0
        tauremgm=abs(taurem)*pptl(4,nptl)/pptl(5,nptl)
        if(taurem.ge.0.)then
        t2=t+tauremgm*(-alog(rangen()))
        else
        t2=t+tauremgm
        endif
        xorptl(1,nptl)=xorptl(1,mm)+pptl(1,nptl)/pptl(4,nptl)*(t2-t)
        xorptl(2,nptl)=xorptl(2,mm)+pptl(2,nptl)/pptl(4,nptl)*(t2-t)
        xorptl(3,nptl)=xorptl(3,mm)+pptl(3,nptl)/pptl(4,nptl)*(t2-t)
        xorptl(4,nptl)=t2
        tivptl(1,nptl)=t2
        call idtau(idptl(nptl),pptl(4,nptl),pptl(5,nptl),taugm)
        tivptl(2,nptl)=tivptl(1,nptl)+taugm*(-alog(rangen()))
        !do l=1,4  ----violates flavor conservation ckw21
        !  ibptl(l,nptl)=0
        !enddo
        andropl=andropl+1
        if(ish.ge.3)write(ifch,*)'Proref,ept(5),id',ept(5),idptl(nptl)
        !print*,' '
        !if(pptl(4,nptl)-pptl(3,nptl).gt.0.
        !.    .and.pptl(4,nptl)+pptl(3,nptl).gt.0.)
        !.   print*,'droplet  ',idptl(nptl),pptl(5,nptl) ,0.5*log(
        !.   (pptl(4,nptl)+pptl(3,nptl))/(pptl(4,nptl)-pptl(3,nptl))   )

        !..........remnant update
        if(gdrop)then  !drop emission: new remnant (string) -> ept, icf
          idrop=0
          amasini=aa(5)*aa(5)
          if(irmdrop.eq.1)then
            call iddeco(icf,jcf)
            call idquacjc(jcf,nqu,naq)
          else
            nqu=1
            naq=1
          endif
          nptl=nptl+1
          t=xorptl(4,mm)
          ifrptl(2,mm)=nptl
          do l=1,5
            pptl(l,nptl)=sngl(aa(l))
          enddo
          idptl(nptl)=idptl(mm)
          istptl(nptl)=40
          iorptl(nptl)=mm
          jorptl(nptl)=0
          ifrptl(1,nptl)=0
          ifrptl(2,nptl)=0
          xorptl(1,nptl)=xorptl(1,mm)
          xorptl(2,nptl)=xorptl(2,mm)
          xorptl(3,nptl)=xorptl(3,mm)
          xorptl(4,nptl)=t
          tivptl(1,nptl)=t
          tivptl(2,nptl)=ainfin
          qsqptl(nptl)=0.
          zpaptl(1,nptl)=zpaptl(1,mm)
          zpaptl(2,nptl)=zpaptl(2,mm)
          if(gproj)then
            ityptl(nptl)=40
          else   !gtarg
            ityptl(nptl)=50
          endif
          do l=1,4
            ibptl(l,nptl)=0
          enddo
          !if(pptl(4,nptl)-pptl(3,nptl).gt.0.
          !.   .and.pptl(4,nptl)+pptl(3,nptl).gt.0.)
          !.   print*,'string/reso  ',idptl(nptl),pptl(5,nptl),0.5*log(
          !.   (pptl(4,nptl)+pptl(3,nptl))/(pptl(4,nptl)-pptl(3,nptl))   )
        endif

        mm=nptlini+1
        nptlb=nptl

c..............................droplet id.....................................

        if(iabs(idptl(mm)).gt.10**8)then

c........................decay droplet condition..................................

         if(irmdrop.eq.1.and.nqc.eq.0)then

          iret=0
          !always decay droplet here (no flow in remnants)
          if(ish.ge.3)write(ifch,*)'Decay remnant droplet...'
          if(nptlb.gt.mxptl-10)call utstop('ProRef: mxptl too small&')

                         !======================  
          if(ifrade.gt.0) call hnbaaa(9,mm,iret) ! 9 = remnant droplets
                         !======================

          if(iret.eq.0.and.nptl.ne.nptlb)then ! ---successful decay---
            istptl(mm)=41
            ifrptl(1,mm)=nptlb+1
            ifrptl(2,mm)=nptl
            t=tivptl(2,mm)
            x=xorptl(1,mm)+(t-xorptl(4,mm))*pptl(1,mm)/pptl(4,mm)
            y=xorptl(2,mm)+(t-xorptl(4,mm))*pptl(2,mm)/pptl(4,mm)
            z=xorptl(3,mm)+(t-xorptl(4,mm))*pptl(3,mm)/pptl(4,mm)
            do 21 n=nptlb+1,nptl
              iorptl(n)=mm
              jorptl(n)=0
              istptl(n)=0
              ifrptl(1,n)=0
              ifrptl(2,n)=0
              radius=0.8*sqrt(rangen())
              phi=2*pi*rangen()
              ti=t
              zi=z
              xorptl(1,n)=x + radius*cos(phi)
              xorptl(2,n)=y + radius*sin(phi)
              xorptl(3,n)=zi
              xorptl(4,n)=ti
              tivptl(1,n)=ti
              call idtau(idptl(n),pptl(4,n),pptl(5,n),taugm)
              r=rangen()
              tivptl(2,n)=tivptl(1,n)+taugm*(-alog(r))
c ityptl set in hnbaaa (41 or 51)
c              if(gproj)then
c                ityptl(n)=40
c                if(gdrop)ityptl(n)=41 !$$$$$ don't change $$$$$$$$$
c              else  !gtarg
c                ityptl(n)=50
c                if(gdrop)ityptl(n)=51 !$$$$$ don't change $$$$$$$$$
c              endif
              if(iept.eq.6)ityptl(n)=ityptl(n)+6
              radptl(n)=0.
              dezptl(n)=0.
              itsptl(n)=0
              rinptl(n)=-9999
              zpaptl(1,n)=zpaptl(1,mm)
              zpaptl(2,n)=zpaptl(2,mm)
   21       continue
            if(iabs(idptl(nptlb+1)).le.6) then
              call gakli2(0,0)
              if(ish.ge.1)write (ifmt,*)'string from drop:nptlb+1,nptl:'
     *                                 ,nptlb+1,nptl
              istptl(nptlb+1)=1
              do n=nptlb+2,nptl
                istptl(n)=20
              enddo
              call gakfra(iret)
              call gakli2(0,0)
            endif
            call setAccumJerr(4,getAccumJerr(4)+1)
          elseif(ifrade.gt.0.and.ispherio.eq.0)then ! Unsuccessful decay
            call setAccumJerr(5,getAccumJerr(5)+1)
            if(ish.ge.4)write(ifch,*)
     *         '***** Unsuccessful remnant cluster decay'
     *             ,' --> do RemoveHadrons instead.'
            mm=mmini
            nptl=nptlini
            irmdropx=0
            goto 55
          endif

          if(idrop.eq.1)goto 1000
         !successful drop decay, no additional string, nothing to do

         else!.............................................

          nqc=0       !just to be sure for the followings

         endif!.......decay droplet condition..............



        elseif(idrop.eq.1)then!...droplet id.....

         goto 1000
         !just a resonance, no additional string, nothing to do

        endif

      endif!............exotic ............

c...............................................................
c for irmdrop=2 the droplet (mm) is used as input for EmitStrings and the string
c form getdrop is the final remnant if existing.
      if(idrop.eq.0.and.irmdropx.eq.1.and.nqc.eq.0)then
c for irmdrop=1 a droplet is emitted (only if no charm) and the string is
c used as input for EmitStrings (to avoid large mass) if present.
        mm=nptlini+2
        do l=1,5
          ept(l)=aa(l)
        enddo
      elseif(gdrop.and.nqc.ne.0)then
c if some charm in droplet EmitStrings is used for the full remnant
        mm=mmini
        do l=1,5
          ept(l)=utdble(pptl(l,mm))
        enddo
        call iddeco(icf,jcfstr)
        do j=1,2
          do i=1,nrflav
            jcf(i,j)=jcf(i,j)+jcfstr(i,j)
          enddo
        enddo
      endif
      istptl(mm)=41
      ifrptl(1,mm)=nptl+1

c........................Emit string.........................

      if(.not.((nqu.eq.3.and.naq.eq.0).or.(nqu.eq.0.and.naq.eq.3)
     &         .or.(nqu.eq.1.and.naq.eq.1))
     & .or.(amasini.gt.xmdrmax.and.ieex.ne.2.and.irmdrop.ne.1)
     & .or.(gdrop.and.irmdrop.eq.2).or.irmdrop.eq.3)then

       call EmitStrings(gproj,gdrop,m,mm,jcf,jcval,icf,ept,xmdrmax,iret)
       if(iret.ne.0)then
         iretxx=1
         goto 1000
       endif

      endif

      if(gdrop.and.idrop.eq.0.and.irmdrop.eq.2)then  !update string momentum
        do l=1,5
          ept(l)=aa(l)
        enddo
        call iddeco(icf,jcf)
        call idquacjc(jcf,nqu,naq)
        mm=nptlini+2
        istptl(mm)=41
        ifrptl(1,mm)=nptl+1
      endif

c........................ determine idr (0=string, else=resonance).......

      if(icf(1).eq.0.and.icf(2).eq.0)then
        id=111               !minimal state is not pi0 but rho0 (from data)
      else
        id=idtra(icf,0,0,0)
!        id=idtra(icf,0,0,3)
      endif
      idr=0
      am=sngl(ept(5))
      call idres(id,am,idr,iadj,1)
c      if(iabs(mod(idr,10)).le.2.and.idr.ne.0)then
c       id=idr
c      else
c       idr=0
c      endif                                !ckeck on-shell mass (see uti)
      if(iadj.ne.0.and.iept.gt.0.and.ept(5).gt.0.d0
     &     .and.(dabs((ept(4)+ept(3))*(ept(4)-ept(3))
     $           -ept(2)**2-ept(1)**2-dble(am)**2).gt.0.3d0))idr=0

      if(ish.ge.3)then
        write(ifch,'(a,5e11.3)')' updt:',(sngl(ept(k)) ,k=1,5)
        write(ifch,*)'            icf: ',icf,' idr: ',idr,' iept: ',iept
      endif

c      if(iept.eq.3)stop'ProReF: iept=3 ???'

c...........................................string...................
      if((iept.gt.0.and.iept.ne.6.and.idr.eq.0).or.ieex.eq.3)then

        !... nqu of remainder string

        anstrg0=anstrg0+1
        if(gdrop)anstrg1=anstrg1+1

        call iddeco(icf,jcf)
        nqu=0
        nqv=0
        nav=0
        do l=1,nrflav
          nqu=nqu+jcf(l,1)-jcf(l,2)
          nqv=nqv+jcval(l,1)+jcval(l,2)
          nav=nav+jcval(l,2)
        enddo

c        if(zrminc.lt.0.)stop'ProReF: not supported any more.         '

        !......determine forward momentum ep


        am1=0.
        am2=0.
        ptt1=0d0
        ptt2=0d0
        pt=ranptdcut(1.)*ptfraqq
        if(pt.lt.0.5d0*ept(5))then
          phi=2.*pi*rangen()
          ptt1=dble(pt*cos(phi))
          ptt2=dble(pt*sin(phi))
        endif
        ptt3=dble(ir)*sqrt((0.5d0*ept(5))**2-ptt1*ptt1-ptt2*ptt2)
c        if(iept.ge.1.and.ept(5).lt.dble(strcut))then
cc for diffraction random direction
c          fact=1.-exp(-min(50.,(strcut/sngl(ept(5))-1.)))
c          phi=(1-2*int(rangen()))*rangen()*pi*fact
c          ax=dble(cos(phi))
c          ay=dble(sin(phi))
c          az=dble(cos((1-2*int(rangen()))*0.5*pi*rangen()*fact))
c          call utrot2(1,ax,ay,az,ptt1,ptt2,ptt3)
c        endif
        ep(1)=ptt1
        ep(2)=ptt2
        ep(3)=ptt3
cc        ep(4)=0.5d0*ept(5)
        ep(4)=sqrt(ptt3*ptt3+ptt2*ptt2+ptt1*ptt1+dble(am1*am1))

        if(ept(5).le.0.)stop'ERROR in ProReF ept(5) LE 0.'
        call utlob2(-1,ept(1),ept(2),ept(3),ept(4),ept(5)
     *     ,ep(1),ep(2),ep(3),ep(4),25)


        xxx=max(1./engy,min(1.,sngl(abs(ep(3)/ep(4)))))
        qqs=q2sft

        !....determine forward and backward flavor icf, icb

        if(iremn.ge.2)then
          xm3val=9.
          xm2val=3.
          xm1val=1.
          ntryx=0
 33       xx1=0.
          xx2=0.
          xx3=0.
          del=1./(1.-alppar)
          if(nqv.eq.3)then
            xx1=min(1.,ranptdcut(xm3val))
            xx2=min(1.,ranptdcut(xm3val))
            xx3=min(1.,ranptdcut(xm3val))
          elseif(nqv.eq.2)then
            xx1=min(1.,ranptdcut(xm2val))
            xx2=min(1.,ranptdcut(xm2val))
            xx3=rangen()**del
          elseif(nqv.eq.1)then
            xx1=min(1.,ranptdcut(xm1val))
            xx2=rangen()**del
            xx3=rangen()**del
          else
            xx1=rangen()**del
            xx2=rangen()**del
            xx3=rangen()**del
          endif
          if(ntryx.lt.1000)then
            if(xx1+xx2+xx3.gt.1)goto 33
          else
            xx1=rangen()
            xx2=rangen()*(1.-xx1)
            xx3=rangen()*(1.-xx1-xx2)
          endif
          xx1=xxx*xx1
          xx2=xxx*xx2
          xx3=xxx*xx3
          piq1=0d0
          piq2=0d0
          piq3=0d0
cKW22          if(ieex.eq.1)then
cKW22c inversion needed for inelastic remnant because of cascade (NA49)
cKW22            ireminv=1
cKW22          else
cKW22            ireminv=0   !no inversion for very low mass diffraction or string from split 
cKW22          endif         !(no explicit remnant excitation) or diffraction (ieex=2)
cKW22
       ireminv=0 
cKW22
cKW22 Otherwise we get extremely weird double peak structures of projectile baryons 
cKW22 (peak at negative rapidity) and of taget baryons (peak at positive rapidity)
cKW22  in central AuAu collisions
cKW22 
       if(nqu.eq.3)then      !---baryon---
          c="s"
          if(nqv.ge.1)c="v"
          iq1=idraflx(piq1,xx1,qqs,iclpt,jcf,jcval,1,isopt,c)
          c="s"
          if(nqv.ge.2)c="v"
          iq2=idraflx(piq2,xx2,qqs,iclpt,jcf,jcval,1,isopt,c)
          c="s"
          if(nqv.ge.3)c="v"
          iq3=idraflx(piq3,xx3,qqs,iclpt,jcf,jcval,1,isopt,c)
          if(.not.gdrop.or.iept.le.3)then
            call neworderx(xx3,xx2,xx1,iq3,iq2,iq1)
c            if(xx2/xx3.gt.xxx**2*reminv*(xx1/xx2))ireminv=0   !xxx**2 to avoid inversion for heavy strings
            if(xx2-xx3.gt.reminv*(xx1-xx2))ireminv=0 
          else
            ireminv=0
          endif
c put always strange quarks in diquark (for lambda and cascade (NA49))
          if(iq3.ge.3)ireminv=1
cKW22
       ireminv=0 
cKW22
cKW22 Otherwise we get extremely weird double peak structures of projectile baryons 
cKW22 (peak at negative rapidity) and of taget baryons (peak at positive rapidity)
cKW22  in central AuAu collisions
cKW22 
          if(ireminv.eq.0)then
            call uticpl(icf,iq3,2,iret) ! antiquark
            call uticpl(icb,iq3,1,iret) ! quark
          else
            call uticpl(icf,iq3,2,iret) ! antiquark
            call uticpl(icb,iq3,1,iret) ! quark
            call uticpl(icf,iq2,2,iret) ! antiquark
            call uticpl(icb,iq2,1,iret) ! quark
          endif
        elseif(nqu.eq.-3)then !---antibaryon---
          c="s"
          if(nqv.ge.1)c="v"
          iq1=idraflx(piq1,xx1,qqs,iclpt,jcf,jcval,2,isopt,c)
          c="s"
          if(nqv.ge.2)c="v"
          iq2=idraflx(piq2,xx2,qqs,iclpt,jcf,jcval,2,isopt,c)
          c="s"
          if(nqv.ge.3)c="v"
          iq3=idraflx(piq3,xx3,qqs,iclpt,jcf,jcval,2,isopt,c)
          if(.not.gdrop.or.iept.le.3)then
            call neworderx(xx3,xx2,xx1,iq3,iq2,iq1)
c            if(xx2/xx3.gt.xxx**2*reminv*(xx1/xx2))ireminv=0   !xxx**2 to avoid inversion for heavy strings
            if(xx2-xx3.gt.reminv*(xx1-xx2))ireminv=0 
          else
            ireminv=0
          endif
c put always strange quarks in diquark
          if(iq3.ge.3)ireminv=1
          if(ireminv.eq.0)then
            call uticpl(icf,iq3,1,iret) ! quark
            call uticpl(icb,iq3,2,iret) ! antiquark
          else
            call uticpl(icf,iq1,1,iret) ! quark
            call uticpl(icb,iq1,2,iret) ! antiquark
            call uticpl(icf,iq2,1,iret) ! quark
            call uticpl(icb,iq2,2,iret) ! antiquark
          endif
        elseif(nqu.eq.0)then !---meson---
          xx3=0.    !no third quark
          iq3=0
          if(nqv.eq.2)then
            c1="v"
            c2="v"
            j=min(2,1+int(0.5+rangen()))
          elseif(nav.ne.0)then    !valence antiquark
            c1="v"
            c2="s"
            j=2
          elseif(nqv.ne.0)then    !valence quark
            c1="v"
            c2="s"
            j=1
          else                    !only sea quarks
            c1="s"
            c2="s"
            j=min(2,1+int(0.5+rangen()))
          endif
          iq1=idraflx(piq1,xx1,qqs,iclpt,jcf,jcval,j,isopt,c1)
          iq2=idraflx(piq2,xx2,qqs,iclpt,jcf,jcval,3-j,isopt,c2)
          if((.not.gdrop.or.iept.le.3).and.xx1.gt.xx2)ireminv=0
         if(ireminv.eq.1)then
            call uticpl(icf,iq1,3-j,iret) ! subtract quark 1 forward
            call uticpl(icb,iq1,j,iret) ! add quark 1 backward
          else
            call uticpl(icf,iq2,j,iret) ! subtract antiquark 2 forward
            call uticpl(icb,iq2,3-j,iret) ! add antiquark 2 backward
          endif
        else
          call utmsg('ProReF')
          write(ifch,*)'***** neither baryon nor antibaryon nor meson.'
          write(ifch,*)'*****  number of net quarks:',nqu
          write(ifmt,*)'ProReF: no hadron; ',nqu,' quarks  --> redo'
          iretxx=2
          goto 1000
        endif
        if(ish.ge.3)write(ifch,'(a,2i3,3(i2,e13.6))')' inversion:',isopt
     &         ,ireminv,iq1,xx1,iq2,xx2,iq3,xx3
        else
        ireminv=0
        if(iept.ne.0)then
          if(rangen().lt.reminv)ireminv=1
        endif
        if(nqu.eq.3)then      !---baryon---
          iq=idrafl(iclpt,jcf,1,'v',1,iret)
          call uticpl(icf,iq,2,iret)       ! antiquark
          call uticpl(icb,iq,1,iret)       ! quark
          if(ireminv.eq.1)then
           iq=idrafl(iclpt,jcf,1,'v',1,iret)
           call uticpl(icf,iq,2,iret)       ! antiquark
           call uticpl(icb,iq,1,iret)       ! quark
          endif
        elseif(nqu.eq.-3)then !---antibaryon---
          iq=idrafl(iclpt,jcf,2,'v',1,iret)
          call uticpl(icf,iq,1,iret)       ! quark
          call uticpl(icb,iq,2,iret)       ! antiquark
          if(ireminv.eq.1)then
           iq=idrafl(iclpt,jcf,2,'v',1,iret)
           call uticpl(icf,iq,1,iret)       ! quark
           call uticpl(icb,iq,2,iret)       ! antiquark
          endif
        elseif(nqu.eq.0)then !---meson---
           iq1=idrafl(iclpt,jcf,1,'v',1,iret)
           iq2=idrafl(iclpt,jcf,2,'v',1,iret)
           if(rangen().gt.0.5)then
             call uticpl(icf,iq1,2,iret) ! subtract quark
             call uticpl(icb,iq1,1,iret) ! add quark
           else
             call uticpl(icf,iq2,1,iret) ! subtract antiquark
             call uticpl(icb,iq2,2,iret) ! add antiquark
           endif
        else
          if(ish.ge.1)then
          call utmsg('ProReF')
          write(ifch,*)'***** neither baryon nor antibaryon nor meson.'
          write(ifch,*)'*****  number of net quarks:',nqu
          endif
          write(ifmt,*)'ProReF: no hadron; ',nqu,' quarks  --> redo'
          iretxx=3
          goto 1000
        endif
      endif


        !..... forward string end

        nptl=nptl+1
        if(nptl.gt.mxptl)call utstop('ProRef: mxptl too small&')
        pptl(1,nptl)=sngl(ep(1))
        pptl(2,nptl)=sngl(ep(2))
        pptl(3,nptl)=sngl(ep(3))
        pptl(4,nptl)=sngl(ep(4))
        pptl(5,nptl)=am1 !0.
        istptl(nptl)=20
        iorptl(nptl)=mm
        if(.not.gdrop)istptl(mm)=41
        jorptl(nptl)=0
        if(.not.gdrop)ifrptl(1,mm)=nptl
        ifrptl(2,mm)=nptl
        t=xorptl(4,mm)
        tauremgm=abs(taurem)*ept(4)/ept(5)
        if(taurem.ge.0.)then
        t2=t+tauremgm*(-alog(rangen()))
        else
        t2=t+tauremgm
        endif
        xorptl(1,nptl)=xorptl(1,mm)+ept(1)/ept(4)*(t2-t)
        xorptl(2,nptl)=xorptl(2,mm)+ept(2)/ept(4)*(t2-t)
        xorptl(3,nptl)=xorptl(3,mm)+ept(3)/ept(4)*(t2-t)
        xorptl(4,nptl)=t2
        tivptl(1,nptl)=t2
        tivptl(2,nptl)=t2
        idptl(nptl)=idtra(icf,0,0,0)
        if(gproj)then
          if(iept.lt.1)stop'ProReF: iep(m)<1     '
          if(iept.eq.3.or.iept.eq.4)then
            if(ieex.eq.2)then
              ityptl(nptl)=43
            else
              ityptl(nptl)=42
            endif
          else
            ityptl(nptl)=41+iept ! =42 =43 =46 =47
          endif
          if(gdrop)ityptl(nptl)=44
        else  !gtarg
          if(iept.lt.1)stop'ProReF: iet(m)<1     '
          if(iept.eq.3.or.iept.eq.4)then
            if(ieex.eq.2)then
              ityptl(nptl)=53
            else
              ityptl(nptl)=52
            endif
          else
            ityptl(nptl)=51+iept !=52 =53 =56 =57
          endif
          if(gdrop)ityptl(nptl)=54
        endif
        itsptl(nptl)=1
        qsqptl(nptl)=0.
        rinptl(nptl)=-9999
        !write(6,'(a,i9,$)')'     ',idptl(nptl) !======================
        zpaptl(1,nptl)=zpaptl(1,mm)
        zpaptl(2,nptl)=zpaptl(2,mm)
        if(ish.ge.3)then
          write(ifch,'(a,5e11.3,$)')' kink:',(pptl(k,nptl),k=1,5)
          write(ifch,*)' id: ',idptl(nptl)
        endif
        !....... backward string end

        nptl=nptl+1
        if(nptl.gt.mxptl)call utstop('ProRef: mxptl too small&')
        pptl2=0.
        do i=1,3
         pptl(i,nptl)=sngl(ept(i)-ep(i))
         pptl2=pptl2+pptl(i,nptl)*pptl(i,nptl)
        enddo
        pptl(5,nptl)=am2 !0.
        pptl2=pptl2+pptl(5,nptl)*pptl(5,nptl)
        pptl(4,nptl)=sqrt(pptl2)
        pptl2=sngl(ept(4)-ep(4))
        if(ish.ge.1.and.abs(pptl2-pptl(4,nptl)).gt.max(0.1,
     &                                         0.1*abs(pptl2)))then
          write(ifmt,*)
     &    'Warning in ProRef: inconsistent backward string end energy !'
     &    ,pptl(4,nptl),pptl2,abs(pptl2-pptl(4,nptl)),am1,am2,ptt3,ep(4)
          if(ish.ge.2)write(ifch,*)
     &    'Warning in ProRef: inconsistent backward string end energy !'
     &    ,(pptl(kkk,nptl),kkk=1,4),pptl2,abs(pptl2-pptl(4,nptl))
        endif
        istptl(nptl)=20
        iorptl(nptl)=mm
        jorptl(nptl)=0
        ifrptl(2,mm)=nptl
        ifrptl(1,nptl)=0
        ifrptl(2,nptl)=0
        t=xorptl(4,mm)
        tauremgm=abs(taurem)*ept(4)/ept(5)
        if(taurem.ge.0.)then
        t2=t+tauremgm*(-alog(rangen()))
        else
        t2=t+tauremgm
        endif
        xorptl(1,nptl)=xorptl(1,mm)+ept(1)/ept(4)*(t2-t)
        xorptl(2,nptl)=xorptl(2,mm)+ept(2)/ept(4)*(t2-t)
        xorptl(3,nptl)=xorptl(3,mm)+ept(3)/ept(4)*(t2-t)
        xorptl(4,nptl)=t2
        tivptl(1,nptl)=t2
        tivptl(2,nptl)=t2
        idptl(nptl)=idtra(icb,0,0,0)
        ityptl(nptl)=ityptl(nptl-1)
        itsptl(nptl)=1
        qsqptl(nptl)=0. !q2nmin+zz  !tp20171124????? do we want that ??? (larger mult for larger Q2 connected to remnant)
        rinptl(nptl)=-9999
        !write(6,'(a,i9)')'     ',idptl(nptl)
        zpaptl(1,nptl)=zpaptl(1,mm)
        zpaptl(2,nptl)=zpaptl(2,mm)
        if(ish.ge.3)then
          write(ifch,'(a,5e11.3,$)')' kink:',(pptl(k,nptl),k=1,5)
          write(ifch,*)' id: ',idptl(nptl)
        endif

c............................no string = resonance...................
      else

        anreso0=anreso0+1
        if(gdrop)anreso1=anreso1+1

        nptl=nptl+1
        if(idr.ne.0)id=idr
        if(nptl.gt.mxptl)call utstop('ProRef: mxptl too small&')
        if(iept.eq.0.or.iept.eq.6)call idmass(id,am)
        idptl(nptl)=id
        pptl(1,nptl)=sngl(ept(1))
        pptl(2,nptl)=sngl(ept(2))
        am2t=sqrt(ept(2)*ept(2)+ept(1)*ept(1)+dble(am*am))
        if(ept(4).gt.am2t)then   !conserve value of E on not pz
          pptl(4,nptl)=sngl(ept(4))
          pptl(3,nptl)=sngl(sign(sqrt((ept(4)+am2t)*(ept(4)-am2t))
     &                          ,ept(3)))
        else
          pptl(3,nptl)=sngl(ept(3))
          pptl(4,nptl)=sngl(sqrt(ept(3)*ept(3)+am2t))
        endif
        pptl(5,nptl)=am
        istptl(nptl)=0
        iorptl(nptl)=mm
        if(.not.gdrop)istptl(mm)=41
        jorptl(nptl)=0
        if(.not.gdrop)ifrptl(1,mm)=nptl
        ifrptl(2,mm)=nptl
        ifrptl(1,nptl)=0
        ifrptl(2,nptl)=0
        t=xorptl(4,mm)
        tauremgm=abs(taurem)*pptl(4,nptl)/pptl(5,nptl)
        if(taurem.ge.0.)then
        t2=t+tauremgm*(-alog(rangen()))
        else
        t2=t+tauremgm
        endif
        xorptl(1,nptl)=xorptl(1,mm)+pptl(1,nptl)/pptl(4,nptl)*(t2-t)
        xorptl(2,nptl)=xorptl(2,mm)+pptl(2,nptl)/pptl(4,nptl)*(t2-t)
        xorptl(3,nptl)=xorptl(3,mm)+pptl(3,nptl)/pptl(4,nptl)*(t2-t)
        xorptl(4,nptl)=t2
        tivptl(1,nptl)=t2
        call idtau(idptl(nptl),pptl(4,nptl),pptl(5,nptl),taugm)
        tivptl(2,nptl)=tivptl(1,nptl)+taugm*(-alog(rangen()))
        if(gproj)then
          ityptl(nptl)=45
          if(gdrop.and.irmdrop.eq.1)then
            ityptl(nptl)=46
          elseif(iept.eq.6)then
            ityptl(nptl)=47
          elseif(ieex.eq.2)then
            ityptl(nptl)=48
          elseif(gdrop)then
            ityptl(nptl)=49
          else      !if iep=0
            mine=0
            mdif=0
            do l=1,lproj(m)
              kp=kproj(m,l)
              if(itpr(kp).eq.-1)then
                mine=1
              elseif(itpr(kp).ne.0)then  !incl. hard diffraction
                mdif=1
              endif
            enddo
            if(mine.eq.0.and.mdif.eq.1)ityptl(nptl)=48
          endif
        else   !gtarg
          ityptl(nptl)=55
          if(gdrop.and.irmdrop.eq.1)then
            ityptl(nptl)=56
          elseif(iept.eq.6)then
            ityptl(nptl)=57
          elseif(ieex.eq.2)then
            ityptl(nptl)=58
          elseif(gdrop)then
            ityptl(nptl)=59
          else      !if iet=0
            mine=0
            mdif=0
            do l=1,ltarg(m)
              kt=ktarg(m,l)
              if(itpr(kt).eq.-1)then
                mine=1
              elseif(itpr(kt).ne.0)then  !incl. hard diffraction
                mdif=1
              endif
            enddo
            if(mine.eq.0.and.mdif.eq.1)ityptl(nptl)=58
          endif
        endif
        itsptl(nptl)=0
        qsqptl(nptl)=0.
        rinptl(nptl)=-9999
        zpaptl(1,nptl)=zpaptl(1,mm)
        zpaptl(2,nptl)=zpaptl(2,mm)

        if(ish.ge.3)write(ifch,'(a,5e11.3,i7)')' nucl:'
     *         ,(pptl(i,nptl),i=1,5),idptl(nptl)

      endif
c.......................................................................
c      print *,iep(1),iet(1),ityptl(nptl)
 1000 call utprix('ProReF',ish,ishini,3)
ctp060829        if(ityptl(nptl).gt.60)print*,ityptl(nptl)
      return

      end

c-----------------------------------------------------------------------
      subroutine EmitStrings(gproj,gdrop,m,mm,jcf,jcv
     &                        ,ic,ept,xmdrmax,iret)
c-----------------------------------------------------------------------
c If remnant mass or quark content is too large, emit enough string to
c share mass and/or quark content via N-body decay. Put valence quarks
c into one of the string to be used in ProReF as final remnant string.
c behavior depends on irmdrop :
c irmdrop = 0 - used to simply avoid remnant with more than 3 quarks
c         = -1- same but split remnant into Nh+1 strings if
c               mass > amdrmax where Nh is the number of hard scatterings
c         < -1- same but split remnant into abs(irmdrop) strings if
c               mass > amdrmax
c         = 1 - a droplet was emitted, limit the mass of the final string
c         = 2 - a string was emitted, remaining energy into strings with 
c               mass < amdrmax
c         = 3 - in case of hard scattering emit a string with M2~Q2s for
c               each hard scattering attached to the remnant
c-----------------------------------------------------------------------
#include "aaa.h"
#include "ems.h"
      parameter(maxp=600)
      common/confg/np,amass(maxp),ident(maxp),pcm(5,maxp),wtxlog,wtlog
      integer jcf(nflav,2),jcv(nflav,2),icf(2),icb(2),ic(2)
     *,jcrem(nflav,2)
      double precision aa(5),ept(5),p(5),xmdrmax,amasini,avmass
     *,ptt1,ptt2,ptt3,t2,t,drangen,xpo
      logical gremn,gdrop,gproj
      dimension amq2(maxp-1)
      common/ems6/ivp0,iap0,idp0,isp0,ivt0,iat0,idt0,ist0
      common /cncl/xproj(mamx),yproj(mamx),zproj(mamx)
     *            ,xtarg(mamx),ytarg(mamx),ztarg(mamx)
      character c*1

      call utpri('emitst',ish,ishini,5)

      iret=0
      irmdropx=irmdrop
      if(irmdropx.ne.2)then
        gremn=.true.
      else
        gremn=.not.gdrop
      endif
      amasini=ept(5)*ept(5)
      if(gproj)then
        ir=1
        icl=iclpro
        if(iremn.ge.2)then
          idrf=idp(m)
        else
          idrf=idp0
        endif
        zz=zzremn(m,1)
        if(iep(m)/10.ne.1)then
          ieex=mod(iep(m),10)      !only for excited remnants
        else
          ieex=0
        endif
        nq2mx=0
        if((irmdrop.eq.3.or.irmdropx.eq.-1)
     .    .and.ieex.eq.3.and.zz.gt.0.)then    !only hard
          do l=1,lproj(m)
            kpair=kproj(m,l)
            do n=1,nprmx(kpair)
              if(idpr(n,kpair).eq.3.and.ivpr(n,kpair).eq.1)then
                nq2mx=nq2mx+1
                amq2(nq2mx)=zdrinc*sqrt(q2kmin(1,n,kpair))
              endif
            enddo
          enddo
        endif
      else
        ir=-1
        icl=icltar
        if(iremn.ge.2)then
          idrf=idt(m)
        else
          idrf=idt0
        endif
        zz=zzremn(m,2)
        if(iet(m)/10.ne.1)then
          ieex=mod(iet(m),10)      !only for excited remnants
        else
          ieex=0
        endif
        nq2mx=0
        if((irmdrop.eq.3.or.irmdropx.eq.-1)
     .    .and.ieex.eq.3.and.zz.gt.0.)then    !only hard
          do l=1,ltarg(m)
            kpair=ktarg(m,l)
            do n=1,nprmx(kpair)
              if(idpr(n,kpair).eq.3.and.ivpr(n,kpair).eq.1)then
                nq2mx=nq2mx+1
                amq2(nq2mx)=zdrinc*sqrt(q2kmin(2,n,kpair))
              endif
            enddo
          enddo
        endif
      endif
      if(irmdropx.le.0)then
        zz=0.                   !in that case we use nq2mx only
        nq2mx=nq2mx-1           !not for single Pom
      endif
      call idquacjc(jcf,nqu,naq)
      if(nq2mx.gt.0.and.zz.gt.0.)then
        amdrmx=sqrt(zz)
      else   !if not hard then limit the mass to amdrmax like irmdrop=1
        amdrmx=amdrmax
        if(irmdropx.eq.3)irmdropx=1
      endif
      if(gremn.and.amasini.lt.(amdrmx+alpdro(2))**2)then
        if((nqu.eq.3.and.naq.eq.0).or.(nqu.eq.0.and.naq.eq.3)
     .                           .or.(nqu.eq.1.and.naq.eq.1))then
          if(ish.ge.7)write(ifch,*)'EmitStrings single string'
     . ,sqrt(amasini),amdrmx,sqrt(zz),nq2mx
          goto 1000
        else    !low mass but not a single string (maximum mass irrelevant)
          amdrmx=0.
        endif
      endif
c first take out valence quark to put in final remnant
      call idquacjc(jcv,nqv,nav)
      if(nqv+nav.ne.0.and.gremn)then    !still some valence quark
        do i=1,nflav
          jcrem(i,1)=jcv(i,1)
          jcrem(i,2)=jcv(i,2)
          jcf(i,1)=jcf(i,1)-jcv(i,1)
          jcf(i,2)=jcf(i,2)-jcv(i,2)
        enddo
        if(.not.(nqv.eq.3.or.nav.eq.3.or.nqv.eq.nav))then   !missing quarks to have a string
          nqu=nqu-nqv
          naq=naq-nav
          if(idrf.eq.1.or.nqv.gt.1.or.nav.gt.1)then !there is a valence diquark
            if(nqv.gt.nav)then
              do ii=1,3-nqv
                if(nqu.gt.0)then
                  c='v'
                  nqu=nqu-1
                else
                  c='s'
                endif
                iq=idrafl(icl,jcf,1,c,1,iretso)
                jcrem(iq,1)=jcrem(iq,1)+1
              enddo
            elseif(nav.gt.nqv)then
              do ii=1,3-nqv
                if(naq.gt.0)then
                  c='v'
                  naq=naq-1
                else
                  c='s'
                endif
                iq=idrafl(icl,jcf,2,c,1,iretso)
                jcrem(iq,2)=jcrem(iq,2)+1
              enddo
            else
              call utstop('Should not happen in EmitString (1) !&')
            endif
          elseif(nqv.gt.nav)then
            if(nqu.gt.0)then
              c='v'
              nqu=nqu-1
            else
              c='s'
            endif
            iq=idrafl(icl,jcf,1,c,1,iretso)
            jcrem(iq,1)=jcrem(iq,1)+1
          elseif(nav.gt.nqv)then
            if(naq.gt.0)then
              c='v'
              naq=naq-1
            else
              c='s'
            endif
            iq=idrafl(icl,jcf,2,c,1,iretso)
            jcrem(iq,2)=jcrem(iq,2)+1
            else
              call utstop('Should not happen in EmitString (2) !&')
          endif
        endif
        idrf=-1 !final string already set
        call idquacjc(jcf,nqu,naq)
      elseif(.not.gremn)then
        idrf=-1                 !final string should not be set
        jcrem(1,1)=0
      else
        jcrem(1,1)=-1
      endif
      if(nqu.eq.naq.and.(nqu.lt.3.or.idrf.le.0))then
        nmes=nqu
        if(idrf.ge.0)nmes=nmes-1             !string is aq-q
        nbar=0
      elseif(nqu.gt.naq)then
        nmes=naq
        nbar=(nqu-nmes)/3     !nbar baryons
        if(idrf.ge.0)then
          if(nmes.eq.0.or.idrf.eq.1)then
            nbar=nbar-1         !string is qq-q
          else
            nmes=nmes-1         !string is aq-q
          endif
        endif
      elseif(nqu.lt.naq)then
        nmes=nqu
        nbar=(naq-nmes)/3    !nbar antibaryons
        if(idrf.ge.0)then
          if(nmes.eq.0.or.idrf.eq.1)then
            nbar=nbar-1         !string is aqaq-aq
          else
            nmes=nmes-1         !string is aq-q
          endif
        endif
      else        !nqu=naq
        nbar=nqu/3
        nmes=nqu-3*nbar
        nbar=2*nbar        !as many antibaryons than baryons
        nbar=nbar-1             !string is qq-q or aqaq-aq
      endif
c compute the additional number of strings
      if(irmdropx.eq.3)then
        nstring=nmes+nbar
      else
        nstring=int(sqrt(amasini/xmdrmax))+1
        if(irmdropx.eq.-1)then
          nstring=min(nq2mx,nstring)
        elseif(irmdropx.le.0)then
          nstring=min(abs(irmdropx)-1,nstring)
        endif
      endif
      nstring=max(nstring,nmes+nbar)
      tecm=sngl(sqrt(amasini))
      tecmx=tecm-alpdro(2)        !for the momentum between strings
      amasq2=(nstring-nbar)*0.5*ammsqq+nbar*ammsqd
      if(gremn)then
        if(nstring.eq.0.and.nq2mx.le.0)goto 1000
        nstring=nstring+1       !should include final remnant
        amasq2=amasq2+sngl(ammn(5))
      elseif(nstring.eq.0)then
        call utstop('Should not happen in EmitString (5) !&')       
      endif
      if(irmdropx.eq.3)amasq2=max(amasq2,amdrmx-amasq2)
      avmass=dble(tecmx-amasq2)/dble(nstring)
      if(ish.ge.3)
     &       write(ifch,*)'EmitStrings part (nq,na,nb,nm,ns,nq2):'
     &  ,nqu,naq,nbar,nmes,nstring,nq2mx,idrf,avmass,amasq2,tecm,gremn
      if(avmass.lt.0d0)then
        if(ish.ge.1.and.ifrade.gt.0)
     .  write(ifmt,*)'EmitStrings part (nb,nm,ns,nq2):'
     .  ,nbar,nmes,nstring,nq2mx,idrf,avmass,amasq2,tecm,gremn,irmdropx
        iret=2
        goto 1000
      endif
      if(nstring.gt.1)gdrop=.true.
      npmax=nstring
      if(gremn)nstring=nstring-1
      np=0
      xpo=0.5d0
      if(nmes.gt.0)then
c  emit meson string
        do mes=1,nmes
          np=np+1
          if(np.lt.npmax.or.tecmx.gt.amdrmax)then
            xpo=0.5d0*(dble(tecm)-(np-1)*avmass)/dble(tecmx)  !adaptative exponent not to have a final mass too large
            amass(np)=sngl(ammn(3)+drangen(dble(np))**xpo*avmass)
            tecmx=tecmx-amass(np)
          else
            amass(np)=tecmx
            tecmx=0.
          endif
            !write(ifch,*)'remove meson',mes,' / ',nmes
          call getstring(1,ident(np),jcf,gproj)
        enddo
        nstring=nstring-nmes
      endif
c remove (anti)baryons
      if(nbar.gt.0)then
        call idquacjc(jcf,nqu,naq)
        do nb=1,nbar
            !write(ifch,*)'remove baryon',nb,' / ',nbar
          np=np+1
          if(np.lt.npmax.or.tecmx.gt.amdrmax)then
            xpo=0.5d0*(dble(tecm)-(np-1)*avmass)/dble(tecmx)
            amass(np)=sngl(ammn(5)+drangen(dble(np))**xpo*avmass)
            tecmx=tecmx-amass(np)
          else
            amass(np)=tecmx
            tecmx=0.
          endif
          prq=float(nqu/3)
          pra=float(naq/3)
          psum=prq+pra
          if(rangen()*psum.le.prq)then !baryon
            call getstring(2,ident(np),jcf,gproj)
            nqu=nqu-3
          else                  !antibaryon
            call getstring(3,ident(np),jcf,gproj)
            naq=naq-3
          endif
        enddo
        nstring=nstring-nbar
      endif
c update flavors of final remnant
      if(jcrem(1,1).eq.-1)then
        do i=1,nflav
          jcrem(i,1)=jcf(i,1)
          jcrem(i,2)=jcf(i,2)
          jcf(i,1)=0
          jcf(i,2)=0
        enddo
      endif

      if(irmdropx.le.0)nq2mx=0         !reset nq2mx to avoid wrong behavior
c remove additonal string for energy sharing
      if(nstring.gt.0.or.nq2mx.gt.0.or.tecmx.gt.amdrmax
     .               .or.(.not.gremn.and.tecmx.gt.0.))then
        frem=0.5
        if(gremn)frem=1.
        if(irmdropx.eq.3)frem=1e10
c  emit meson string
        mes=0
        do 100 while((mes.lt.nstring.or.tecmx.gt.frem*sngl(avmass).or.
     .  (.not.gremn.and.tecmx.gt.0.).or.nq2mx.gt.0).and.np.lt.maxp-1)
          mes=mes+1
          np=np+1
          if(irmdropx.eq.3)then
            amass(np)=max(sngl(ammn(3)),amq2(nq2mx))
            nq2mx=nq2mx-1
          else
            if(np.le.npmax)
     .           xpo=0.5d0*(dble(tecm)-(np-1)*avmass)/dble(tecmx)
            amass(np)=sngl(ammn(3)+drangen(dble(np))**xpo*avmass)
          endif
          tecmxx=tecmx-amass(np)
          if(.not.gremn.and.(tecmxx.le.frem*sngl(avmass)
     .                       .or.np.eq.maxp-1))then
            amass(np)=tecmx
            tecmx=0.
          else
            tecmx=tecmxx
          endif
          if(ish.ge.7)write(ifch,*)'remove string',mes,' / ',nstring,xpo
     .                             ,amass(np),tecmx,tecmxx,amdrmax,nq2mx
          if(gremn.and.tecmx.lt.ammn(5))then  !not enough energy for remnant: cancel last one
            tecmx=tecmx+amass(np)
            np=np-1
            nq2mx=0
            goto 100
          else
            call getstring(1,ident(np),jcf,gproj)
          endif
 100      continue
      endif
      
      amtot=0.
      if(gremn)then
        npmax=np
        np=np+1                 !for the final remnant
        amass(np)=tecmx
        amtot=amtot+tecmx
        tecmx=0.
      else
        npmax=np
      endif

      if(np.eq.1)then
        do j=1,4
          pcm(j,np)=ept(j)
        enddo
        amass(np)=ept(5)
      else
c do N-body decay
        call hnbody
      endif

      ip=mm
      do n=1,npmax
        do j=1,4
          p(j)=pcm(j,n)
        enddo
        p(5)=amass(n)
        amtot=amtot+amass(n)
        if(np.gt.1)call utlob2(-1,ept(1),ept(2),ept(3),ept(4),ept(5)
     .                           ,p(1),p(2),p(3),p(4),20)
        idr=0
        id=ident(n)
        am=sngl(p(5))
        call idres(id,am,idr,iadj,1)
        if(iadj.ne.0.and.am.gt.0.d0
     &     .and.(dabs((p(4)+p(3))*(p(4)-p(3))
     $           -p(2)**2-p(1)**2-p(5)**2).gt.0.3d0))idr=0
        if(ish.ge.3)then
          write(ifch,'(a,5e11.3)')' updt:',(sngl(p(k)) ,k=1,5)
          write(ifch,*)'            id: ',id,'idr: ',idr
        endif

c...........................................string...................
        if(idr.eq.0)then
          anstrg1=anstrg1+1
          call idtr4(id,icf)           !update icf
          call idtr4(99,icb)           !reset icb
        !... nqu of remainder string

          call iddeco(icf,jcf)
          call idquacjc(jcf,nqu,naq)

        !......determine forward momentum ep
          am1=0.
          am2=0.
          ptt1=0d0
          ptt2=0d0
c        pt=ranptdcut(1.)*ptfraqq
c        if(pt.lt.0.5d0*ept(5))then
c          phi=2.*pi*rangen()
c          ptt1=dble(pt*cos(phi))
c          ptt2=dble(pt*sin(phi))
c        endif
          ptt3=dble(ir)*0.5d0*p(5)
c        ptt3=dble(ir)*sqrt((0.5d0*p(5))**2-ptt1*ptt1-ptt2*ptt2)
          aa(1)=ptt1
          aa(2)=ptt2
          aa(3)=ptt3
          aa(4)=sqrt(ptt3*ptt3+ptt2*ptt2+ptt1*ptt1+dble(am1*am1))

          call utlob2(-1,p(1),p(2),p(3),p(4),p(5)
     *         ,aa(1),aa(2),aa(3),aa(4),21)

        !....determine forward and backward flavor icf, icb

          ireminv=0
          if(rangen().lt.0.5)ireminv=1
          if(nqu.eq.3)then      !---baryon---
            iq=idrafl(iclpt,jcf,1,'v',1,iret)
            call uticpl(icf,iq,2,iret) ! antiquark
            call uticpl(icb,iq,1,iret) ! quark
            if(ireminv.eq.1)then
              iq=idrafl(iclpt,jcf,1,'v',1,iret)
              call uticpl(icf,iq,2,iret) ! antiquark
              call uticpl(icb,iq,1,iret) ! quark
            endif
          elseif(naq.eq.3)then  !---antibaryon---
            iq=idrafl(iclpt,jcf,2,'v',1,iret)
            call uticpl(icf,iq,1,iret) ! quark
            call uticpl(icb,iq,2,iret) ! antiquark
            if(ireminv.eq.1)then
              iq=idrafl(iclpt,jcf,2,'v',1,iret)
              call uticpl(icf,iq,1,iret) ! quark
              call uticpl(icb,iq,2,iret) ! antiquark
            endif
          elseif(nqu.eq.naq)then !---meson---
            iq1=idrafl(iclpt,jcf,1,'v',1,iret)
            iq2=idrafl(iclpt,jcf,2,'v',1,iret)
            if(rangen().gt.0.5)then
              call uticpl(icf,iq1,2,iret) ! subtract quark
              call uticpl(icb,iq1,1,iret) ! add quark
            else
              call uticpl(icf,iq2,1,iret) ! subtract antiquark
              call uticpl(icb,iq2,2,iret) ! add antiquark
            endif
          else
            call utstop('Should not happen in EmitString (4) !&')
          endif


        !..... forward string end

          nptl=nptl+1
          if(nptl.gt.mxptl)call utstop('EmitString: mxptl too small&')
          pptl(1,nptl)=sngl(aa(1))
          pptl(2,nptl)=sngl(aa(2))
          pptl(3,nptl)=sngl(aa(3))
          pptl(4,nptl)=sngl(aa(4))
          pptl(5,nptl)=am1      !0.
          istptl(nptl)=20
          iorptl(nptl)=mm
          jorptl(nptl)=0
          ifrptl(2,mm)=nptl
          ifrptl(1,nptl)=0
          ifrptl(2,nptl)=0
          t=xorptl(4,mm)
          tauremgm=abs(taurem)*p(4)/p(5)
          if(taurem.ge.0.)then
            t2=t+tauremgm*(-alog(rangen()))
          else
            t2=t+tauremgm
          endif
          xorptl(1,nptl)=xorptl(1,mm)+p(1)/p(4)*(t2-t)
          xorptl(2,nptl)=xorptl(2,mm)+p(2)/p(4)*(t2-t)
          xorptl(3,nptl)=xorptl(3,mm)+p(3)/p(4)*(t2-t)
          xorptl(4,nptl)=t2
          tivptl(1,nptl)=t2
          tivptl(2,nptl)=t2
          idptl(nptl)=idtra(icf,0,0,0)
          if(gproj)then
            ityptl(nptl)=44
          else
            ityptl(nptl)=54
          endif
          itsptl(nptl)=1
          qsqptl(nptl)=0.
          rinptl(nptl)=-9999
        !write(6,'(a,i9,$)')'     ',idptl(nptl) !======================
          zpaptl(1,nptl)=zpaptl(1,mm)
          zpaptl(2,nptl)=zpaptl(2,mm)
          if(ish.ge.3)then
            write(ifch,'(a,5e11.3,$)')' kink:',(pptl(k,nptl),k=1,5)
            write(ifch,*)' id: ',idptl(nptl)
          endif
        !....... backward string end

          nptl=nptl+1
          if(nptl.gt.mxptl)call utstop('ProRef: mxptl too small&')
          pptl2=0.
          do i=1,3
            pptl(i,nptl)=sngl(p(i)-aa(i))
            pptl2=pptl2+pptl(i,nptl)*pptl(i,nptl)
          enddo
          pptl(5,nptl)=am2      !0.
          pptl2=pptl2+pptl(5,nptl)*pptl(5,nptl)
          pptl(4,nptl)=sqrt(pptl2)
          pptl2=sngl(p(4)-aa(4))
          if(ish.ge.1.and.abs(pptl2-pptl(4,nptl)).gt.max(0.1,
     &                                        0.1*abs(pptl2)))then
            write(ifmt,*)
     &'Warning in EmitString: inconsistent backward string end energy !'
     &    ,pptl(4,nptl),pptl2,abs(pptl2-pptl(4,nptl)),am1,am2,ptt3,p(4)
            if(ish.ge.2)write(ifch,*)
     &'Warning in EmitString: inconsistent backward string end energy !'
     &    ,(pptl(kkk,nptl),kkk=1,4),pptl2,abs(pptl2-pptl(4,nptl))
          endif
          istptl(nptl)=20
          iorptl(nptl)=mm
          jorptl(nptl)=0
          ifrptl(2,mm)=nptl
          ifrptl(1,nptl)=0
          ifrptl(2,nptl)=0
          t=xorptl(4,mm)
          tauremgm=abs(taurem)*p(4)/p(5)
          if(taurem.ge.0.)then
            t2=t+tauremgm*(-alog(rangen()))
          else
            t2=t+tauremgm
          endif
          xorptl(1,nptl)=xorptl(1,mm)+p(1)/p(4)*(t2-t)
          xorptl(2,nptl)=xorptl(2,mm)+p(2)/p(4)*(t2-t)
          xorptl(3,nptl)=xorptl(3,mm)+p(3)/p(4)*(t2-t)
          xorptl(4,nptl)=t2
          tivptl(1,nptl)=t2
          tivptl(2,nptl)=t2
          idptl(nptl)=idtra(icb,0,0,0)
          ityptl(nptl)=ityptl(nptl-1)
          itsptl(nptl)=1
          qsqptl(nptl)=0. !q2nmin+zz
          rinptl(nptl)=-9999
          zpaptl(1,nptl)=zpaptl(1,mm)
          zpaptl(2,nptl)=zpaptl(2,mm)
          if(ish.ge.3)then
            write(ifch,'(a,5e11.3,$)')' kink:',(pptl(k,nptl),k=1,5)
            write(ifch,*)' id: ',idptl(nptl)
          endif

c............................no string = resonance...................
        else

          anreso1=anreso1+1

          nptl=nptl+1
          if(idr.ne.0)id=idr
          if(nptl.gt.mxptl)call utstop('EmitString: mxptl too small&')
          idptl(nptl)=id
          pptl(1,nptl)=sngl(p(1))
          pptl(2,nptl)=sngl(p(2))
          am2t=sqrt(p(2)*p(2)+p(1)*p(1)+dble(am*am))
          if(p(4).gt.am2t)then !conserve value of E on not pz
            pptl(4,nptl)=sngl(p(4))
            pptl(3,nptl)=sngl(sign(sqrt((p(4)+am2t)*(p(4)-am2t))
     &                          ,p(3)))
          else
            pptl(3,nptl)=sngl(p(3))
            pptl(4,nptl)=sngl(sqrt(p(3)*p(3)+am2t))
          endif
          pptl(5,nptl)=am
          istptl(nptl)=0
          iorptl(nptl)=mm
          jorptl(nptl)=0
          ifrptl(2,mm)=nptl
          ifrptl(1,nptl)=0
          ifrptl(2,nptl)=0
          t=xorptl(4,mm)
          tauremgm=abs(taurem)*pptl(4,nptl)/pptl(5,nptl)
          if(taurem.ge.0.)then
            t2=t+tauremgm*(-alog(rangen()))
          else
            t2=t+tauremgm
          endif
          xorptl(1,nptl)=xorptl(1,mm)+pptl(1,nptl)/pptl(4,nptl)*(t2-t)
          xorptl(2,nptl)=xorptl(2,mm)+pptl(2,nptl)/pptl(4,nptl)*(t2-t)
          xorptl(3,nptl)=xorptl(3,mm)+pptl(3,nptl)/pptl(4,nptl)*(t2-t)
          xorptl(4,nptl)=t2
          tivptl(1,nptl)=t2
          call idtau(idptl(nptl),pptl(4,nptl),pptl(5,nptl),taugm)
          tivptl(2,nptl)=tivptl(1,nptl)+taugm*(-alog(rangen()))
          if(gproj)then
            ityptl(nptl)=49
          else                  !gtarg
            ityptl(nptl)=59
          endif
          itsptl(nptl)=0
          qsqptl(nptl)=0.
          rinptl(nptl)=-9999
          zpaptl(1,nptl)=zpaptl(1,mm)
          zpaptl(2,nptl)=zpaptl(2,mm)

          if(ish.ge.3)write(ifch,'(a,5e11.3,i7)')' nucl:'
     *         ,(pptl(i,nptl),i=1,5),idptl(nptl)
          
        endif
c.......................................................................
      enddo
      if(abs(tecm-amtot)/tecm.gt.alpdro(2))then
        write(ifch,*)' mass:',gremn,tecmx,tecm,amtot,nq2mx
        write(ifmt,*)' mass:',gremn,tecmx,tecm,amtot,nq2mx
        iret=3
        goto 1000
      endif
      
c update final remnant
      if(gremn)then
      call idquacjc(jcrem,nqv,nav)
      if(.not.((nqv*nav.eq.0.and.nqv+nav.eq.3)
     .     .or.(nqv.eq.nav.and.nqv+nav.eq.2)))then       !security check
        write(ifch,*)'jcrem,nqv,nav',jcrem,nqv,nav
        call utstop('Should not happen in EmitString (3) !&')
      endif
      call idenco(jcrem,ic,iret) !save flavor of final remnant
      if(np.gt.1)then
        do j=1,4
          p(j)=pcm(j,np)
        enddo
        p(5)=amass(np)
        call utlob2(-1,ept(1),ept(2),ept(3),ept(4),ept(5)
     .                ,p(1),p(2),p(3),p(4),20)
        do j=1,5
          ept(j)=p(j)
        enddo
      endif
      if(ish.ge.6)then
        write(ifch,30)'string remnant:',ic,ept
      endif
 30   format(a,/,'icf:',i7,' |',i7,/,'ept:',5(e10.3,1x))
      endif
 1000 continue
      if(iret.ne.0.and.ish.ge.1)write(ifch,*)'Pb in EmitStrings !',iret
      call utprix('emitst',ish,ishini,5)

      end

c------------------------------------------------------------------
         subroutine getstring(imb,idf,jc,gproj)
c------------------------------------------------------------------
c       goal:  emit a string (imb= 1 meson, 2 baryon, 3 antibaryon)
c              update the remnant flavor and 5-momentum
c
c       idf ,a : string id and 5-momentum
c       jc, ep : remnant flavor and 5-momentum
c       iret   : in case of error, keep correct momentum in remnant
c                and lose the quarks of the (not) emitted hadron
c-----------------------------------------------------------------

#include "aaa.h"
#include "ems.h"
        integer jc(nflav,2),ifh(3),ic(2)
        common /ems12/iodiba,bidiba  ! defaut iodiba=0. if iodiba=1, study H-Dibaryon
        logical gproj
        character*1 c


        ic(1)=0
        ic(2)=0

        if(gproj)then
          iclpt=iclpro
        else
          iclpt=icltar
        endif

        idf=0
c  get the id and mass of hadron, the remnant jc is updated
        iret2=0
          if(imb.eq.1)then              ! a meson
            call idquacjc(jc,nqu,naq)    !in case jc is empty
            c='v'
            if(nqu+naq.eq.0)c='s'
            j=1
            i=idrafl(iclpt,jc,j,c,1,iret2)
            if(iret2.ne.0)goto 77
            ifq=i
            j=2
            i=idrafl(iclpt,jc,j,'v',1,iret2)
            if(iret2.ne.0)goto 77
            ifa=i
c            write(ifch,*)'ici',ifq,ifa,jc
            ic(1)=10**(6-ifq)
            ic(2)=10**(6-ifa)
            ier=0
            idf=idtra(ic,ier,idum,0)
            if(ier.ne.0)then
              if(ifq.le.ifa)then
                idf=ifq*100+ifa*10
              else
                idf=-(ifq*10+ifa*100)
              endif
            endif

          elseif(imb.eq.2)then            ! a baryon
            j=1
            do ik=1,3
              i=idrafl(iclpt,jc,j,'v',1,iret2)
              if(iret2.ne.0)goto 77
              ifh(ik)=i
              ic(j)=ic(j)+10**(6-i)
            enddo
            ier=0
            idf=idtra(ic,ier,idum,0)
            if(ier.ne.0)then
              call neworder(ifh(1),ifh(2),ifh(3))
              idf=ifh(1)*1000+ifh(2)*100+ifh(3)*10
              if(ifh(1).ne.ifh(2).and.ifh(2).ne.ifh(3)
     $             .and.ifh(1).ne.ifh(3))  idf=2130
              if(ifh(1).eq.ifh(2).and.ifh(2).eq.ifh(3))idf=idf+1
            endif

          elseif(imb.eq.3)then           ! an antibaryon
            j=2
            do ik=1,3
              i=idrafl(iclpt,jc,j,'v',1,iret2)
              if(iret2.ne.0)goto 77
              ifh(ik)=i
              ic(j)=ic(j)+10**(6-i)
            enddo
            ier=0
            idf=idtra(ic,ier,idum,0)
            if(ier.ne.0)then
              call neworder(ifh(1),ifh(2),ifh(3))
              idf=ifh(1)*1000+ifh(2)*100+ifh(3)*10
              if(ifh(1).ne.ifh(2).and.ifh(2).ne.ifh(3)
     $             .and.ifh(1).ne.ifh(3))  idf=2130
              if(ifh(1).eq.ifh(2).and.ifh(2).eq.ifh(3))idf=idf+1
              idf=-idf
            endif
           else
            call utstop('This imb does not exist in getstr !&')
           endif

   77     if(iret2.ne.0)then
          write(ifmt,*)'warning in getstring: imb=',imb,'  iclpt:',iclpt
          write(ifmt,*)'   jc: ',jc,'  j: ',j,'   (1=q,2=aq)  --> redo'
          call utmsg('getstr')
          write(ifch,*)'Not enough quark ??? ... redo event !'
          call utmsgf
          iret=1
          goto 1000
          endif

          if(ish.ge.5)then
            write(ifch,*)'get string with id:',idf
          endif


        if(ish.ge.5)then
          write(ifch,*)'new remnant flavor:',jc,iret
        endif
        iret=0
c          write(ifmt,*)'get hadron with id:',idf
c          write(ifmt,*)'new remnant flavor:',jc

 1000 continue

      return
      end



c------------------------------------------------------------------
         subroutine getdroplet(ir,iept,ic,jc,z,ep,a,pass,xmdrmax)
c------------------------------------------------------------------
c  emit a droplet, update the remnant string flavor and 5-momentum
c
c input
c       ir ........ 1  projectile, -1  target remnant
c       iept ...... particle excitation
c       ep ........ remnant  5-momentum
c       jc ........ remnant jc
c       z  ........ Z factor from splitting
c output
c       pass ...  .true. = successful droplet emission
c                            jc, ep ....... droplet  ic and 5-momentum
c                            ic, a ........ remnant string jc and 5-momentum
c                 .false. = unsuccessful
c                            jc, ep .... unchanged,
c                            considered as droplet jc and 5-momentum
c-----------------------------------------------------------------

#include "aaa.h"
#include "ems.h"
        double precision ep(5),a(5),p1(5),re(5),eps,amasex,xmdrmax
        double precision xxx,rr,alp,p5sq,xmin,xmax,ampt2str
     &  ,sxini,strmas,xxxmax,xxxmin,ampt2dro,xmdrmaxi
        parameter(eps=1.d-20)
        integer jc(nflav,2),ic(2),icx(2)
        integer jcini(nflav,2),jcfin(nflav,2)
        logical pass
        common/cems5/plc,s
        double precision s,plc,ptm,qcm,u(3),utpcmd,ptt,drangen,phi

        call utpri('getdro',ish,ishini,4)

        iret=0
        iret2=0
        xmdrmaxi=min(100.d0,xmdrmax)
        pass=.true.
        idps=0
        idms=0
        do i=1,nflav
          jcini(i,1)=jc(i,1)
          jcini(i,2)=jc(i,2)
          jcfin(i,1)=0
          jcfin(i,2)=0
        enddo


        call idquacjc(jcini,nqu,naq)

        do i=1,5
          a(i)=0.d0
          re(i)=0.d0
        enddo
        npart=nqu+naq
        nqc=jcini(4,1)+jcini(4,2)

        if(ir.eq.1)then
           iclpt=iclpro
         else
           iclpt=icltar
         endif

         if(ish.ge.5)then
           write(ifch,10)'remnant flavor and 5-momentum:'
     &                    ,jc,ep,nqu,naq,nqc,iept
 10        format(a,/,'jc:',6i3,' |',6i3,/,'ep:',5(e10.3,1x),/,4i4)
         endif

c  get id of string ends, the remnant string jc is updated
         if(iremn.eq.3)then  !  remnant content=string content (droplet empty)

           do i=1,nflav
             jcfin(i,1)=jcini(i,1)
             jcfin(i,2)=jcini(i,2)
             jcini(i,1)=0
             jcini(i,2)=0
           enddo

         else

         if(npart.lt.3.and.ep(5).lt.xmdrmaxi.and.nqc.eq.0)then !light droplet with few quarks
            pass=.false.
            goto 1000
         elseif(npart.lt.3)then    !few quarks but heavy, add some quarks to extract a q-qbar string (should not exit directly because of the large mass)
           ifq=idrafl(iclpt,jcini,2,'r',3,iret2)
           if(nqu.eq.1.and.naq.eq.1)then
             idps=1
             idms=1
             nqu=2
             naq=2
           else
             call utstop('This should not happen (getdrop) !&')
           endif
         elseif(nqu.eq.2.and.naq.eq.2)then
           idps=1
           idms=1
         elseif(naq.eq.0)then
           idps=5
           idms=1
         elseif(nqu.eq.0)then
           idps=1
           idms=5
         else                 !There is enough q or aq to do qq-q string


           if(jcini(4,1)-jcini(4,2).eq.0)then !if c-cbar

             idps=1
             idms=1

           else

c One chooses the first q or aq

           rrr=rangen()
           npart=nqu+naq
           if(jcini(4,1)+jcini(4,2).ne.0)then !if some charm take it out
             if(jcini(4,1).ne.0)then
               idps=1
               nqu=nqu-1
             else
               idms=1
               naq=naq-1
             endif
           elseif(rrr.gt.float(naq)/float(npart))then
             idps=1
             nqu=nqu-1
           else
             idms=1
             naq=naq-1
           endif

c One chooses the second one

           rrr=rangen()
           npart=nqu+naq
           if(idps.eq.1.and.jcini(4,1).ne.0)then !if some charm take it out
             idps=5
           elseif(idms.eq.1.and.jcini(4,2).ne.0)then !if some charm take it out
             idms=5
           elseif(rrr.gt.float(naq)/float(npart))then
             if(idps.eq.1.and.nqu.ge.2)then
               idps=5
             else
               idps=1
             endif
           else
             if(idms.eq.1.and.naq.ge.2)then
               idms=5
             else
               idms=1
             endif
           endif

c If there is already 2 q or 2 aq as string end, we know that we need
c a third one to complete the string

           if(idps.eq.5)idms=1
           if(idms.eq.5)idps=1
           if(idps.eq.1.and.idms.ne.5)idms=1
           if(idms.eq.1.and.idps.ne.5)idps=1

         endif

         endif

         if(ish.ge.5)then
           write(ifch,*)'remnant string ends :',idps,idms
         endif

          if(idps.ne.5.and.idms.ne.5)then              ! q-aq string
            if(jcini(4,1).eq.1)then
              ifq=idrafl(iclpt,jcini,1,'c',1,iret)
            else
              ifq=idrafl(iclpt,jcini,1,'v',1,iret)
            endif
            if(jcini(4,1).eq.1)then
              ifa=idrafl(iclpt,jcini,2,'c',1,iret)
            else
              ifa=idrafl(iclpt,jcini,2,'v',1,iret)
            endif
            jcfin(ifq,1)=1
            jcfin(ifa,2)=1

          elseif(idps.eq.5)then                       ! qq-q string
            do ik=1,3
              if(jcini(4,1).ne.0)then
                i=idrafl(iclpt,jcini,1,'c',1,iret)
              else
                i=idrafl(iclpt,jcini,1,'v',1,iret)
              endif
              jcfin(i,1)=jcfin(i,1)+1
            enddo

          elseif(idms.eq.5)then                        !aqaq-aq string
            do ik=1,3
              if(jcini(4,2).ne.0)then
                i=idrafl(iclpt,jcini,2,'c',1,iret)
              else
                i=idrafl(iclpt,jcini,2,'v',1,iret)
              endif
              jcfin(i,2)=jcfin(i,2)+1
            enddo
          endif

          endif      !iremn=3

          if(iret.ne.0)call utstop('Not enough quark in getdro ???&')
          if(jcini(4,1)+jcini(4,2).ne.0)
     &         call utstop('There is sitll charm quark in getdro???&')

c string id

         call idenco(jcfin,icx,iret)
         if(iret.eq.1)then
           call utstop('Exotic flavor in getdroplet !&')
         endif


c boost remnant in rest frame
      if(ish.ge.6) write (ifch,*) 'on-shell check'
        do k=1,5
          p1(k)=ep(k)
        enddo
        p1(5)=(p1(4)-p1(3))*(p1(4)+p1(3))-p1(2)**2-p1(1)**2
        if(p1(5).gt.0d0.and.abs(p1(5)-ep(5)*ep(5)).lt.ep(5))then
          p1(5)=sqrt(p1(5))
        else
          if(ish.ge.2)write(ifch,*)'Precision problem in getdro, p:',
     &             (p1(k),k=1,5),ep(5)*ep(5)
          p1(5)=ep(5)
          p1(4)=sqrt(p1(3)*p1(3)+p1(2)*p1(2)+p1(1)*p1(1)+p1(5)*p1(5))
        endif
      if(ish.ge.6) write (ifch,*) 'boost vector:',p1

c limits for momenta

      mamod=4
      mamos=4
      fad=alpdro(1)
      if(iremn.eq.3)fad=fad+z*zdrinc
      fad=max(1.5,fad)
      ptm=p1(5)
      amasex=dble(fad+utamnz(jcini,mamod))
      fas=2.
      if(iremn.eq.3)then
        id=idtra(icx,ier,ires,0)
        if(ier.eq.0)then
          call idmass(id,amass)           !minimum is particle mass
          strmas=dble(amass)
        else
          strmas=dble(fas*utamnz(jcfin,mamos))
        endif
      else
        strmas=dble(fas*utamnz(jcfin,mamos))
      endif


c redo

       nredo=0
 777   continue
       nredo=nredo+1
       if(nredo.eq.10)then
          amasex=1.5d0*dble(utamnz(jcini,mamod))
          if(iremn.ne.3)strmas=1.5d0*dble(utamnz(jcfin,mamos))
       elseif(nredo.gt.20)then
          !write(ifch,*)'nredo.gt.20 -> only drop'
         if(ish.ge.4)write(ifch,*)
     &     'Pb with string mass in Getdrop, continue with gethad'
          pass=.false.
         goto 1000
       endif

c fix pt

          sxini=ptm*ptm
          ptt=dble(ranpt()*alpdro(2))**2         !pt
          if(ptt.ge.sxini)goto 777
          sxini=sqrt(sxini-ptt)


          ampt2dro=amasex**2d0
          ampt2str=strmas**2d0
          if(ampt2dro.gt.xmdrmaxi)then
            xmdrmaxi=2d0*ampt2dro
c            write(ifmt,*)'Warning Mmin>Mmax in Getdroplet'
          endif

          xxxmax=min(xmdrmaxi,(sxini-strmas)**2)    !strmas/(strmas+ampt2)
          xxxmin=ampt2dro

          if(xxxmin.gt.xxxmax)then
            !write(ifch,*)'Warning Mmin>sxini -> only drop'
           if(ish.ge.4)write(ifch,*)
     &     'Pb with ampt2 in Getdrop, retry',nredo,ir
     &             ,ampt2dro,ampt2str,xxxmin,xxxmax,sxini,ptt,xmdrmaxi
            goto 777
          endif



c fix mass

            rr=drangen(xxxmax)
            xmax=xxxmax
            xmin=xxxmin
            alp=dble(alpdro(3))
            if(dabs(alp-1.d0).lt.eps)then
              xxx=xmax**rr*xmin**(1d0-rr)
            else
              xxx=(rr*xmax**(1d0-alp)+(1d0-rr)*xmin**(1d0-alp))
     &                                                **(1d0/(1d0-alp))
            endif


c        write(ifch,*)'ini',xmin,xxx,xmax,rr,ampt2dro
c    &                   ,(sxini-sqrt(xxx)),ampt2str,p1(5)



          re(5)=sqrt(xxx)
          a(5)=sxini-re(5)
          if(a(5).lt.strmas)then
            if(ish.ge.6)write(ifch,*)
     &           'Pb with initial mass in Getdrop, retry',ir
     &       ,xmin,xxx,xmax,rr,ampt2dro,ampt2str,a(5)
            goto 777
          endif


c two body decay
          if(ish.ge.6)write(ifch,*)'2 body decay',ptm,re(5),a(5)
          qcm=utpcmd(ptm,re(5),a(5),iret)
          u(3)=0.d0 !2.d0*drangen(qcm)-1.d0
          phi=2.d0*dble(pi)*drangen(u(3))
          u(1)=sqrt(1.d0-u(3)**2)*cos(phi)
          u(2)=sqrt(1.d0-u(3)**2)*sin(phi)
          if(u(3).lt.0d0)then          !send always droplet backward
c          if(u(3).gt.0d0)then          !send always droplet forward     ?????
            do j=1,3
              re(j)=qcm*u(j)
              a(j)=-re(j)
            enddo
          else
            do j=1,3
              a(j)=qcm*u(j)
              re(j)=-a(j)
            enddo
          endif

          re(4)=sqrt(qcm**2+re(5)**2)
          a(4)=sqrt(qcm**2+a(5)**2)

          if(ish.ge.6)write(ifch,*)'momentum in rest frame : ',re,a



c Fix a of string

c boost string in collision frame
        call utlob2(-1,p1(1),p1(2),p1(3),p1(4),p1(5)
     $       ,a(1),a(2),a(3),a(4),71)

         p5sq=(a(4)+a(3))*(a(4)-a(3))-(a(1)**2.d0+a(2)**2.d0)
         if(p5sq.gt.ampt2str)then
           a(5)=sqrt(p5sq)
         else
           if(ish.ge.6)then
             write(ifch,*)'Pb with string mass -> retry'
             write(ifch,*)'   m^2:',p5sq,'  m_min^2:',ampt2str
             write(ifch,*)'   momentum four vector:',(a(ii),ii=1,4)
           endif
           goto 777
         endif

c Fix ep of droplet

c boost droplet in collision frame
        call utlob2(-1,p1(1),p1(2),p1(3),p1(4),p1(5)
     $       ,re(1),re(2),re(3),re(4),72)

         p5sq=(re(4)+re(3))*(re(4)-re(3))-(re(1)*re(1)+re(2)*re(2))
         if(p5sq.gt.ampt2dro)then
           re(5)=sqrt(p5sq)
         else
           if(ish.ge.6)then
             write(ifch,*)'Pb with droplet mass -> retry'
             write(ifch,*)'   m^2:',p5sq,'  m_min^2:',ampt2dro
             write(ifch,*)'   momentum four vector:',(re(ii),ii=1,4)
           endif
           goto 777
         endif


       if(ish.ge.1.and.abs(ep(4)-re(4)-a(4)).gt.1.d-2*ep(4))then
         write(ifmt,*)'Pb with energy conservation in getdro'
         if(ish.ge.6)then
           write(ifch,*)'Pb with energy conservation :'
           write(ifch,*)'   p1_ini:',ep(1),'  p1:',re(1)+a(1)
           write(ifch,*)'   p2_ini:',ep(2),'  p2:',re(2)+a(2)
           write(ifch,*)'   p3_ini:',ep(3),'  p3:',re(3)+a(3)
         endif
       endif

c If OK, save flavors of droplet and string
         do i=1,5
           ep(i)=re(i)
         enddo
         ic(1)=icx(1)
         ic(2)=icx(2)
         do i=1,nflav
           jc(i,1)=jcini(i,1)
           jc(i,2)=jcini(i,2)
         enddo

         if(ish.ge.6)then
           write(ifch,20)'droplet:',jc,ep
           write(ifch,30)'string remnant:',ic,a
         endif
 20      format(a,/,'jc:',6i3,' |',6i3,/,'ep:',5(e10.3,1x))
 30      format(a,/,'ic:',i7,' |',i7,/,'a:',5(e10.3,1x))

 1000    continue
         call utprix('getdro',ish,ishini,4)
         end

c------------------------------------------------------------------
      subroutine getdropx(ir,iept,m,ic,jc,jcv,z,ep,a,pass,xmdrmax)
c------------------------------------------------------------------
c  emit a droplet taking into account momentum fraction without screening,
c  update the remnant string flavor and 5-momentum (to be used with iremn=2)
c  calculate xprmd as if it was energy of a remnant connected to a lot of
c  unscreened pomerons. We use zdrinc to limit this effect for xprmd not to 
c  be too small. xmrmd will be small because droplet mass is small. On the 
c  other hand xmrms is large and xprms vary, So the string will be quite 
c  central and the droplet a bit more forward. To have enough particles in the
c  central region with ap/p close to 1, both string and droplet should be
c  quite central and carry part of the baryon number. So the quark content of
c  the string is defined as if would be the remnant and no inversion will
c  done in proref to avoid small ap/p ratio
c
c  In case of 2 body decay, the string is the remnant and the diquark should
c  stay forward (no inversion). The droplet will be very central (low mass 
c  compared to string) and not completely baryon free (compensation of diquark 
c  in central string ends).
c
c  This subroutine is quite important for multiplicity distributions which
c  are better reproduce if a component increasing with "centrality" is added
c  to the Pomerons them self (multiplicity per pomeron more or less fixed by
c  e+e- data). The multiplicity being measured at mid-rapidity, this component
c  has to be at least partially central (shape of dN/deta give another 
c  constrain)
c
c input
c       ir ........ 1  projectile, -1  target remnant
c       iept ...... particle excitation
c       m  ........ remnant index
c       ep ........ remnant  5-momentum
c       jc ........ remnant jc
c       jcv ....... remnant jc valence quark
c       z  ........ Z factor from splitting
c output
c       pass ...  .true. = successful droplet emission
c                            jc, ep ....... droplet  ic and 5-momentum
c                            ic, a ........ remnant string jc and 5-momentum
c                 .false. = unsuccessful
c                            jc, ep .... unchanged,
c                            considered as droplet jc and 5-momentum
c-----------------------------------------------------------------

#include "aaa.h"
#include "ems.h"
      double precision alptr
      common/cemsr5/alptr(0:1,0:5)
      double precision ep(5),a(5),p1(5),re(5),eps,amasex,xmdrmax
      double precision xxx,rr,alpm,p5sq,xmin,xmax,ampt2str,xmsmax
     &,sxini,strmas,xxxmax,xxxmin,ampt2dro,xmdrmaxi,xprmi,xmrmi
     &,xprmd,xmrmd,xprms,xmrms,xpti,ypti,xptd,yptd,xpts,ypts,xptt,yptt
      double precision om1xpr,xremd,xrems,freduc,omint,xmrem
     &     ,atil(ntymi:ntymx),btilp(ntymi:ntymx),btilpp(ntymi:ntymx)
      parameter(eps=1.d-20)
      integer jc(nflav,2),jcv(nflav,2),ic(2),icx(2)
      integer jcvini(nflav,2),jcini(nflav,2)
     &     ,jcfin(nflav,2),jcvfin(nflav,2)
      logical pass
      common/cems5/plc,s
      double precision s,plc,ptm,qcm,u(3),utpcmd,ptt,drangen,phi
      logical strcomp,valqu

      call utpri('getdrx',ish,ishini,4)

      iret=0
      xmdrmaxi=min(100.d0,xmdrmax)
c      xmdrmaxi=xmdrmax
      pass=.true.
      idps=0
      idms=0
      do i=1,nflav
        jcini(i,1)=jc(i,1)
        jcini(i,2)=jc(i,2)
        jcvini(i,1)=jcv(i,1)
        jcvini(i,2)=jcv(i,2)
        jcfin(i,1)=0
        jcfin(i,2)=0
        jcvfin(i,1)=0
        jcvfin(i,2)=0
      enddo


      call idquacjc(jcini,nqu,naq)
      call idquacjc(jcvini,nqv,nav)

      do i=1,5
        a(i)=0.d0
        re(i)=0.d0
      enddo
      nqc=jcini(4,1)+jcini(4,2)

      idrf=0
      if(nqu-naq.ne.0)idrf=1
      if(ir.eq.1)then
        iclpt=iclpro
        idpt=idproj
        amremn=amproj
        if(idrf.eq.0)idrf=idp(m) !change it only if not 1
        xprmi=xpp(m)
        xmrmi=xmp(m)
        xpti=xxp(m)
        ypti=xyp(m)
        if(lproj3(m).gt.0)then
          nlnk=max(1,nint(z))
        else
          nlnk=0
        endif
      else
        iclpt=icltar
        idpt=idtarg
        amremn=amtarg
        if(idrf.eq.0)idrf=idt(m) !change it only if not 1
        xprmi=xmt(m)
        xmrmi=xpt(m)
        xpti=xxt(m)
        ypti=xyt(m)
        if(ltarg3(m).gt.0)then
          nlnk=max(1,nint(z))
        else
          nlnk=0
        endif
      endif

      if(ish.ge.5)then
        write(ifch,10)'remnant flavor and 5-momentum:'
     &       ,jc,jcv,ep,nqu,naq,nqv,nav,nqc,idrf,iept,nlnk
 10     format(a,/,'jc:',6i3,' |',6i3,/,'jcv:',6i3,' |',6i3,/
     &       ,'ep:',5(e10.3,1x),/,8i4)
      endif
      
c limits for momenta

      mamod=4
      mamos=4
      fad=alpdro(1)
      if(iept.ge.3)fad=fad+z*zdrinc
      amasex=1.2d0*dble(fad+utamnz(jcini,mamod))
      fas=2.
      strmas=dble(fas*max(amremn,utamnz(jcfin,mamos)))
      ampt2dro=amasex**2d0
      ampt2str=strmas**2d0
      if(ampt2dro.gt.xmdrmaxi)then
        xmdrmaxi=2d0*ampt2dro
c       write(ifmt,*)'Warning Mmin>Mmax in getdropx'
      endif

c     check formation conditions

      strcomp=.false.
      valqu=.true.              !if true, valence quark will always be in strings : reduce lambda production
      if((nqu.eq.3.and.naq.eq.0).or.(nqu.eq.0.and.naq.eq.3)
     &     .or.(nqu.eq.1.and.naq.eq.1))then
c       not enough quark for the droplet, check mass
        if((iept.ne.5.or.ep(5).lt.amasex)
     .     .and.ep(5)*ep(5).lt.xmdrmaxi.and.nqc.eq.0)then
          pass=.false.          !continue without string emission
          if(ish.ge.4)write(ifch,*)
     & 'Normal remnant in getdropx, continue only with droplet ...'
          goto 1000
        endif
c       create q-aq from sea (but no charm)
        if(nqu.eq.naq.and.idrf.eq.1)nlnk=max(2,nlnk)
        do n=1,nlnk
          idum=idrafl(iclpt,jcini,1,'r',3,iret)
          nqu=nqu+1
          naq=naq+1
        enddo
        strcomp=.true.
        valqu=.false.
      elseif(mod(nqu-naq,3).ne.0)then
        call utstop('This should not happen (getdropx) !&')
      elseif(idrf.ne.0.and.nqu+naq.eq.4)then
c       if 2 q + 2 aq and diquark needed, add a q and an aq
        idum=idrafl(iclpt,jcini,1,'r',3,iret)
        nqu=nqu+1
        naq=naq+1
      endif

c     get id of string ends, the remnant string jc is updated

c     First remove all charm

      if(nqc.ne.0.and.jcini(4,1)-jcini(4,2).eq.0)then !if c-cbar

        if(jcini(4,1).eq.1)then
          idps=1
          idms=1
        else
          call utstop('getdropx can not manage more than c-cb !&')
        endif

      elseif(nqc.ne.0.and.jcini(4,1)*jcini(4,2).ne.0)then

        call utstop('getdropx can not manage c quarks this way !&')

      else


        if(nqc.ne.0)then        !if some charm take it out
          if(jcini(4,1).ne.0)then
            if(nqu.lt.3)then
              idrf=0            !can not use c in antibaryon
            elseif(jcini(4,1).gt.1)then
              idrf=1            !more than 1 c quark only in baryon
            endif
          elseif(jcini(4,2).ne.0)then
            if(naq.lt.3)then
              idrf=0            !can not use cb in baryon
            elseif(jcini(4,2).gt.1)then
              idrf=1            !more than 1 c antiquark only in antibaryon
            endif
          endif
          if(idrf.ne.0.and.jcini(4,1).gt.0.and.jcini(4,1).le.3)then
            idps=5
            idms=1
          elseif(idrf.ne.0.and.jcini(4,2).gt.0.and.jcini(4,2).le.3)then
            idps=1
            idms=5
          elseif(jcini(4,1).gt.1.or.jcini(4,2).gt.1)then
            call utstop('getdropx can not use more than 3 c/cb !&')
          endif
        endif

c       take into account number of diquark in final remnant string

        if(idps.eq.0)then

          if(idrf.ne.0.or.naq*nqu.eq.0)then !use a diquark
            if(valqu.and.(nqv*(nqu-2).gt.0.or.nav*(naq-2).gt.0))then !if valence quark, fix diquark type
              if(nqv.eq.nav)then
                if(nqu.lt.3)then
                  idps=1
                  idms=5
                elseif(naq.lt.3)then
                  idps=5
                  idms=1
                else
                  if(rangen().lt.0.5)then
                    idps=1
                    idms=5
                  else
                    idps=5
                    idms=1
                  endif
                endif
              elseif(nqv.gt.0.and.nqu.ge.3)then
                idps=5
                idms=1
              elseif(nav.gt.0.and.naq.ge.3)then
                idps=1
                idms=5
              else
                call utstop("Should not happen in getdropx !&")
              endif
            elseif(valqu.and.nqu.ge.3.and.idpt.gt.1000)then !qq-q
              idps=5
              idms=1
            elseif(valqu.and.naq.ge.3.and.idpt.lt.-1000)then !qbqb-qb
              idps=1
              idms=5
            else
              if((rangen().lt.0.5.and.naq.ge.3).or.nqu.lt.3)then
                idps=1
                idms=5
              else
                idps=5
                idms=1
              endif
            endif
          else                  !q-qb
            idps=1
            idms=1
          endif

        endif

      endif                     !string end type

      if(ish.ge.5)then
        write(ifch,'(a,2i2,3x,6i3,2x,6i3)')
     &  'remnant string ends and jcini:',idps,idms,jcini
      endif

c     choose flavor with priority to valence quark (after charm)

      if(idps.ne.5.and.idms.ne.5)then ! q-aq string
        j=1
        if(jcini(4,j).gt.0)then
          i=4
          jcini(i,j)=jcini(i,j)-1
          if(jcvini(i,j).gt.0)then
            jcvfin(i,j)=jcvfin(i,j)+1
            jcvini(i,j)=jcvini(i,j)-1
            nqv=nqv-1
          endif
        elseif(valqu.and.nqv.gt.0)then
          i=idraflz(jcvini,j)
          jcvfin(i,j)=jcvfin(i,j)+1
          jcini(i,j)=jcini(i,j)-1
          nqv=nqv-1
        else
          i=idrafl(iclpt,jcini,j,'v',1,iret)
          if(jcini(i,j)-jcvini(i,j).lt.0)then
            jcvini(i,j)=jcvini(i,j)-1
            jcvfin(i,j)=jcvfin(i,j)+1
          endif
        endif
        ifq=i
        j=2
        if(jcini(4,j).gt.0)then
          i=4
          jcini(i,j)=jcini(i,j)-1
          if(jcvini(i,j).gt.0)then
            jcvfin(i,j)=jcvfin(i,j)+1
            jcvini(i,j)=jcvini(i,j)-1
            nav=nav-1
          endif
        elseif(valqu.and.nav.gt.0)then
          i=idraflz(jcvini,j)
          jcvfin(i,j)=jcvfin(i,j)+1
          jcini(i,j)=jcini(i,j)-1
          nav=nav-1
        else
          i=idrafl(iclpt,jcini,j,'v',1,iret)
          if(jcini(i,j)-jcvini(i,j).lt.0)then
            jcvini(i,j)=jcvini(i,j)-1
            jcvfin(i,j)=jcvfin(i,j)+1
          endif
        endif
        ifa=i
        jcfin(ifq,1)=1
        jcfin(ifa,2)=1

      elseif(idps.eq.5)then     ! qq-q string
        j=1
        do ik=1,3
          if(jcini(4,j).ne.0)then
            i=4
            jcini(i,j)=jcini(i,j)-1
            if(jcvini(i,j).gt.0)then
              jcvfin(i,j)=jcvfin(i,j)+1
              jcvini(i,j)=jcvini(i,j)-1
              nqv=nqv-1
            endif
          elseif(valqu.and.nqv.gt.0)then
            i=idraflz(jcvini,j)
            jcvfin(i,j)=jcvfin(i,j)+1
            jcini(i,j)=jcini(i,j)-1
            nqv=nqv-1
          else
            i=idrafl(iclpt,jcini,j,'v',1,iret)
            if(jcini(i,j)-jcvini(i,j).lt.0)then
              jcvini(i,j)=jcvini(i,j)-1
              jcvfin(i,j)=jcvfin(i,j)+1
            endif
          endif
          jcfin(i,j)=jcfin(i,j)+1
        enddo

      elseif(idms.eq.5)then     !aqaq-aq string
        j=2
        do ik=1,3
          if(jcini(4,j).gt.0)then
            i=4
            jcini(i,j)=jcini(i,j)-1
            if(jcvini(i,j).gt.0)then
              jcvfin(i,j)=jcvfin(i,j)+1
              jcvini(i,j)=jcvini(i,j)-1
              nav=nav-1
            endif
          elseif(valqu.and.nav.gt.0)then
            i=idraflz(jcvini,j)
            jcvfin(i,j)=jcvfin(i,j)+1
            jcini(i,j)=jcini(i,j)-1
            nav=nav-1
          else
            i=idrafl(iclpt,jcini,j,'v',1,iret)
            if(jcini(i,j)-jcvini(i,j).lt.0)then
              jcvini(i,j)=jcvini(i,j)-1
              jcvfin(i,j)=jcvfin(i,j)+1
            endif
          endif
          jcfin(i,j)=jcfin(i,j)+1
        enddo

      endif

      if(iret.ne.0)call utstop('Not enough quark in getdropx ???&')
      if(jcini(4,1)+jcini(4,2).ne.0)
     &     call utstop('There is sitll charm quark in getdropx???&')


c     boost remnant in rest frame
      if(ish.ge.6) write (ifch,*) 'on-shell check'
      do k=1,5
        p1(k)=ep(k)
      enddo
      p1(5)=(p1(4)-p1(3))*(p1(4)+p1(3))-p1(2)**2-p1(1)**2
      if(p1(5).gt.0d0.and.abs(p1(5)-ep(5)*ep(5)).lt.ep(5))then
        p1(5)=sqrt(p1(5))
      else
        if(ish.ge.2)write(ifch,*)'Precision problem in getdropx, p:',
     &       (p1(k),k=1,5),ep(5)*ep(5)
        p1(5)=ep(5)
        p1(4)=sqrt(p1(3)*p1(3)+p1(2)*p1(2)+p1(1)*p1(1)+p1(5)*p1(5))
      endif


c string id

      call idenco(jcfin,icx,iret)
      if(iret.eq.1)then
        call utstop('Exotic flavor in getdropx !&')
      endif


      ptm=p1(5)
      xxxmin=1d0-xprmi     !max x+ for string
      
      !########################################################################
      !#############following line uncommented in version 262 #################
      !########################################################################
           nlnk=0              !force to use only 2 body decay
      !########################################################################
      !########################################################################
      
      nredo=-1
      freduc=1d0
 
      if(ish.ge.4)write(ifch,'(a,i4,2x,2e12.5,3x,6i3,2x,6i3)')
     &  'Start splitting :',nlnk,xprmi,xmrmi,jcfin
777   continue
       nredo=nredo+1

c redo

       if(strcomp.and.nredo.eq.20)then  !after 19 try and remnant compatible with a string
         pass=.false.         !continue without droplet
         if(ish.ge.4)write(ifch,'(a,2i3,4e12.5)')
     &     'Pb with splitting in getdropx, continue without split ...'
     &     ,nlnk,nvirt,xxxmax,xxxmin,ep(5)**2,xmdrmaxi
         goto 1000
       elseif((nlnk.gt.0.and.nredo.eq.10).or.nredo.eq.30)then        !reduce minimum mass
          amasex=1.2d0*dble(rangen()+utamnz(jcini,mamod))
          strmas=1.1d0*dble(max(amremn,utamnz(jcfin,mamos)))
          ampt2dro=amasex**2d0
          ampt2str=strmas**2d0
       elseif(nredo.eq.20)then    !after 19 try, use 2 body decay
         if(nlnk.gt.0.and.ish.ge.4)write(ifch,*)
     &     'nredo>20, use 2 body decay ...',nvirt,xxxmax,xxxmin
         nlnk=0                 !force to use only 2 body decay
         amasex=dble(alpdro(1)+utamnz(jcini,mamod))
         strmas=dble(fas*max(amremn,utamnz(jcfin,mamos)))
         ampt2dro=amasex**2d0
         ampt2str=strmas**2d0
         if(ish.ge.6) write (ifch,*) 'boost vector:',p1
       elseif(nredo.ge.40)then
          !write(ifch,*)'nredo.gt.20 -> only drop'
         if(ish.ge.4)write(ifch,*)
     &    'Pb with mass in getdropx, continue without split ...'
          pass=.false.
         goto 1000
       endif

       if(xxxmin.gt.ampt2str/(s*xmrmi))then
         xmrms=ampt2str/(s*xxxmin) !min x- for droplet
       else
         nlnk=0
       endif

       if(xmrms.lt.xmrmi.and.nlnk.gt.0)then        !kinetic compatibility

         xmrmd=xmrmi-xmrms     !max x- for droplet

c fix the virtual number of collision (no screening)
         imin=ntymin
         imax=ntymax
         nvirt=0
         xxxmax=0d0
         xptd=0d0
         yptd=0d0
         if(ir.eq.1)then
           do l=1,lproj3(m)    !use all pairs attached to remnant
             kp=kproj3(m,l)
             xmrem=xmt(itarg(kp))
             omint=0d0  !number of interaction per pair (no screening)
             do i=imin,imax
               atil(i)=atilde(i,kp)/s**dble(epsilongs(i,kp))
               btilp(i)=btildep(i,kp)-dble(epsilongp(i,kp))
               btilpp(i)=btildepp(i,kp)-dble(epsilongt(i,kp))
               omint=omint+atil(i)/(btilp(i)+1.d0)/(btilpp(i)+1.d0)
     &            /(2.d0+btilp(i)) /(2.d0+btilpp(i))
             enddo
             omint=omint*(0.5+1.5*drangen(omint))*dble(freduc)  !introduce fluctuations
c             omint=omint**zdrinc  !if omint (and then xxxmax) is too large, only little mass in string (and more momentum in droplet)
             nkol=int(omint)
             if(drangen(omint).lt.omint-dble(nkol))nkol=nkol+1
             nvirt=nvirt+nkol   !total number of virtual Pomerons
             do n=1,nkol
               xprmd=1d0-xxxmax
c take x from an "unscreened" Pomeron (reduction factor if too high)
               xxx=om1xpr(atil,btilp,btilpp,xprmd,xmrem,ir)
               ptt=dble(ranptdcut(1.)*alpdro(2))
c             ptt=dble(ranptdcut(ptsems)*alpdro(2))
               phi=2d0*dble(pi)*drangen(ptt)
               xprmd=1d0-xxxmax-xxx
               xptt=xptd+ptt*cos(phi)
               yptt=yptd+ptt*sin(phi)
               xremd=xprmd*xmrmd*s-(xpti+xptt)**2-(ypti+yptt)**2
               if(xremd.gt.ampt2dro)then
                 xxxmax=xxxmax+xxx
                 xptd=xptt
                 yptd=yptt
               endif
             enddo
           enddo
         else
           do l=1,ltarg3(m)    !use all pairs attached to remnant
             kt=ktarg3(m,l)
             xmrem=xpp(iproj(kt))
             omint=0d0  !number of interaction per pair (no screening)
             do i=imin,imax
               atil(i)=atilde(i,kt)/s**dble(epsilongs(i,kt))
               btilp(i)=btildep(i,kt)-dble(epsilongp(i,kt))
               btilpp(i)=btildepp(i,kt)-dble(epsilongt(i,kt))
               omint=omint+atil(i)/(btilp(i)+1.d0)/(btilpp(i)+1.d0)
     &            /(2.d0+btilp(i)) /(2.d0+btilpp(i))
             enddo
             omint=omint*(0.5+1.5*drangen(omint))*dble(freduc)  !introduce fluctuations
c             omint=omint**zdrinc  !if omint (and then xxxmax) is too large, only little mass in string (and more momentum in droplet)
             nkol=int(omint)
             if(drangen(omint).lt.omint-dble(nkol))nkol=nkol+1
             nvirt=nvirt+nkol   !total number of virtual Pomerons
             do n=1,nkol
               xprmd=1d0-xxxmax
c take x from an "unscreened" Pomeron (reduction factor if too high)
               xxx=om1xpr(atil,btilp,btilpp,xprmd,xmrem,ir)
               ptt=dble(ranptdcut(1.)*alpdro(2))
c              ptt=dble(ranptdcut(ptsems)*alpdro(2))
               phi=2d0*dble(pi)*drangen(ptt)
               xprmd=1d0-xxxmax-xxx
               xptt=xpts+ptt*cos(phi)
               yptt=ypts+ptt*sin(phi)
               xremd=xprmd*xmrmd*s-(xpti+xptt)**2-(ypti+yptt)**2
               if(xremd.gt.ampt2dro)then
                 xxxmax=xxxmax+xxx
                 xptd=xptt
                 yptd=yptt
               endif
             enddo
           enddo
         endif

         if(ish.ge.4)write(ifch,*)'Virtual collisions :'
     &        ,nvirt,xxxmax.le.xxxmin,xxxmin,xxxmax,xptd,yptd

         if(xxxmax.le.xxxmin)goto 777


c check string mass and energy

         xprmd=1d0-xxxmax       !x+ droplet
         xprms=xprmi-xprmd      !x+ string
         xpts=xpti-xptd         !pt droplet
         ypts=ypti-yptd         !pt droplet
         xptd=xpti+xptd    !pt string
         yptd=ypti+yptd    !pt string
         xmrms=(ampt2str+xpts*xpts+ypts*ypts)/xprms/s  !update min x- string
c droplet mass should not exceed to much mdrmaxi. Use random to smooth distrib.
         xmax=min(xmdrmaxi*(1.+drangen(xmdrmaxi)),
     &            xprmd*(xmrmi-xmrms)*s-xptd*xptd-yptd*yptd)  !maximum mass allowed by minium string mass
         xmin=ampt2dro+xptd*xptd+yptd*yptd
         if(xmin.ge.xmax)then
           if(ish.ge.4)write(ifch,*)
     &      'Pb with minimum mass in getdropx, retry',nredo,ir
     &          ,xprms,xmrms,xmrmi-xmrms,xmax,xmin,xptd,yptd
           freduc=freduc*0.8d0
           goto 777
         endif

c fix droplet mass not to have too heavy droplets (too long to decay)

         
         rr=drangen(xxx)
         alpm=alptr(idrf,1)        !droplet=remnant, string=virtual pomerons from screening
         if(dabs(alpm-1.d0).lt.eps)then
           xxx=xmax**rr*xmin**(1d0-rr)
         else
           xxx=(rr*xmax**(1d0-alpm)+(1d0-rr)*xmin**(1d0-alpm))
     &          **(1d0/(1d0-alpm))
         endif

         if(ish.ge.6)write(ifch,*)'Fix droplet mass'
     &          ,xmin,xxx,xmax,alpm

         xmrmd=xxx/xprmd/s
         xremd=xprmd*xmrmd*s-xptd*xptd-yptd*yptd

c check droplet mass and energy
         if(xremd.lt.ampt2dro.or.xmrmd.ge.xmrmi-xmrms)then
           if(ish.ge.4)write(ifch,*)
     &          'Pb with droplet mass in getdropx, retry',nredo,ir
     &          ,ampt2dro,xremd,xprmd,xmrmd,xmrmi-xmrms,xptd,yptd
           goto 777
         endif

c check string mass and energy

         xmrms=xmrmi-xmrmd   !final x- string
         xrems=xprms*xmrms*s-xpts*xpts-ypts*ypts
         if(xrems.lt.ampt2str)then
           if(ish.ge.4)write(ifch,*)
     &          'Pb with string mass in getdropx, retry',nredo,ir
     &          ,ampt2str,xrems,xprms,xmrms,xpts,ypts
           goto 777
         endif

         if(ish.ge.4)then
           write(ifch,*)'String Mass :'
     &          ,ampt2str,xrems,xprms,xmrms,xpts,ypts
           write(ifch,*)'Droplet Mass :'
     &          ,ampt2dro,xremd,xprmd,xmrmd,xptd,yptd
         endif

         re(1)=xptd
         re(2)=yptd
         re(3)=dble(ir)*(xprmd-xmrmd)*plc*0.5d0
         re(4)=(xprmd+xmrmd)*plc*0.5d0
         re(5)=sqrt(xremd)

         a(1)=xpts
         a(2)=ypts
         a(3)=dble(ir)*(xprms-xmrms)*plc*0.5d0
         a(4)=(xprms+xmrms)*plc*0.5d0
         a(5)=sqrt(xrems)

c remove valence quark from string
c
c         do i=1,nflav
c           jcvini(i,1)=jcv(i,1)
c           jcvini(i,2)=jcv(i,2)
c           jcvfin(i,1)=0
c           jcvfin(i,2)=0
c         enddo


       else   !if xm to small, use two body decay (should be rare)

         if(ish.ge.6)write (ifch,*)'kinematic limit -> boost vector:',p1

c fix pt

          sxini=ptm*ptm
          ptt=dble(ranpt()*alpdro(2))**2         !pt
          if(ish.ge.7)write (ifch,*)'pt:',sxini,ptt
          if(ptt.ge.sxini)goto 777
          sxini=sqrt(sxini-ptt)



          xmsmax=xmdrmaxi*(1.+drangen(xmdrmaxi))
          xxxmax=min(xmsmax,(sxini-strmas)**2)    !strmas/(strmas+ampt2)
          xxxmin=ampt2dro

          if(xxxmin.gt.xxxmax)then
            !write(ifch,*)'Warning Mmin>sxini -> only drop'
           if(ish.ge.4)write(ifch,*)
     &     'Pb with ampt2 in getdropx, retry',nredo,ir
     &       ,ampt2dro,ampt2str,xxxmin,xxxmax,sxini**2,strmas,xmsmax
            goto 777
          endif


          if(irmdrop.eq.1)then
c fix mass droplet

            rr=drangen(xxxmax)
            xmax=xxxmax
            xmin=xxxmin
            alpm=dble(alpdro(3))
            if(dabs(alpm-1.d0).lt.eps)then
              xxx=xmax**rr*xmin**(1d0-rr)
            else
              xxx=(rr*xmax**(1d0-alpm)+(1d0-rr)*xmin**(1d0-alpm))
     &                                                **(1d0/(1d0-alpm))
            endif


c        write(ifch,*)'ini',xmin,xxx,xmax,rr,ampt2dro
c     &                   ,(sxini-sqrt(xxx)),ampt2str,p1(5)

          re(5)=sqrt(xxx)
          a(5)=sxini-re(5)
          if(a(5).lt.strmas)then
            if(ish.ge.6)write(ifch,*)
     &           'Pb with initial mass in getdropx, retry',ir
     &       ,xmin,xxx,xmax,rr,ampt2dro,ampt2str,a(5)
            goto 777
          endif


        else
c fix mass string
          
          xmin=strmas**2
          xmax=(sxini-sqrt(xxxmin))**2
          xmax=min(xmax,xmdrmax)

          rr=drangen(xxx)
          alpm=alptr(idrf,2)
          if(dabs(alpm-1.d0).lt.eps)then
            xxx=xmax**rr*xmin**(1d0-rr)
          else
            xxx=(rr*xmax**(1d0-alpm)+(1d0-rr)*xmin**(1d0-alpm))
     &           **(1d0/(1d0-alpm))
          endif
          
          a(5)=sqrt(xxx)

          re(5)=sxini-a(5)
          if(re(5).lt.amasex)then
            if(ish.ge.6)write(ifch,*)
     &           'Pb with initial mass in getdropx, retry',ir
     &       ,xmin,xxx,xmax,rr,ampt2dro,ampt2str,re(5)
            goto 777
          endif

          if(a(5)+re(5).ge.ptm)then
            write(ifmt,*)
     &           'Pb with string mass in Getdropx, retry',ir
            if(ish.ge.6)write(ifch,*)
     &           'Pb with string mass in Getdropx, retry',ir
     &       ,xmin,xxx,xmax,rr,ampt2dro,ampt2str,a(5),re(5),ptm
            goto 777
          endif
        endif

c two body decay
          if(ish.ge.6)write(ifch,*)'2 body decay',ptm,re(5),a(5)
          qcm=utpcmd(ptm,re(5),a(5),iret)
c          ptt=dble(ranptdcut(1.)*alpdro(2))
c          u(3)=sqrt(max(0d0,(1d0+ptt/qcm)*(1d0-ptt/qcm)))
          u(3)=2.d0*drangen(qcm)-1.d0
          phi=2.d0*dble(pi)*drangen(u(3))
          u(1)=sqrt(1.d0-u(3)**2)*cos(phi)
          u(2)=sqrt(1.d0-u(3)**2)*sin(phi)
          if(u(3).lt.0.d0)then    !doplet always backward, string is used as remnant here
c          if(u(3).gt.0.d0)then    !doplet always forward, string is used as remnant here and has higher mass so it is still forward
c          if(rangen().gt.0.5)then    !random choice
            do j=1,3
              re(j)=qcm*u(j)
              a(j)=-re(j)
            enddo
          else
            do j=1,3
              a(j)=qcm*u(j)
              re(j)=-a(j)
            enddo
          endif

          re(4)=sqrt(qcm**2+re(5)**2)
          a(4)=sqrt(qcm**2+a(5)**2)

          if(ish.ge.6)write(ifch,*)'momentum in rest frame : ',qcm,re,a



c Fix a of string

c boost string in collision frame
        call utlob2(-1,p1(1),p1(2),p1(3),p1(4),p1(5)
     $       ,a(1),a(2),a(3),a(4),73)

         p5sq=(a(4)+a(3))*(a(4)-a(3))-(a(1)**2.d0+a(2)**2.d0)
         if(p5sq.gt.ampt2str)then
           a(5)=sqrt(p5sq)
         else
           if(ish.ge.6)then
             write(ifch,*)'Pb with string mass -> retry'
             write(ifch,*)'   m^2:',p5sq,'  m_min^2:',ampt2str
             write(ifch,*)'   momentum four vector:',(a(ii),ii=1,4)
           endif
           goto 777
         endif

c Fix ep of droplet

c boost droplet in collision frame
        call utlob2(-1,p1(1),p1(2),p1(3),p1(4),p1(5)
     $       ,re(1),re(2),re(3),re(4),74)

         p5sq=(re(4)+re(3))*(re(4)-re(3))-(re(1)*re(1)+re(2)*re(2))
         if(p5sq.gt.ampt2dro)then
           re(5)=sqrt(p5sq)
         else
           if(ish.ge.6)then
             write(ifch,*)'Pb with droplet mass -> retry'
             write(ifch,*)'   m^2:',p5sq,'  m_min^2:',ampt2dro
             write(ifch,*)'   momentum four vector:',(re(ii),ii=1,4)
           endif
           goto 777
         endif

         ! jerr(...) removed !!!!!!!!!!!!!!!!!
         
       endif     !test of xm

       if(ish.ge.1.and.abs(ep(4)-re(4)-a(4)).gt.1.e-2*ep(4))then
         write(ifmt,*)'Pb with energy conservation in getdropx'
         if(ish.ge.6)then
           write(ifch,*)'Pb with energy conservation :'
           write(ifch,*)'   p1_ini:',ep(1),'  p1:',re(1)+a(1)
           write(ifch,*)'   p2_ini:',ep(2),'  p2:',re(2)+a(2)
           write(ifch,*)'   p3_ini:',ep(3),'  p3:',re(3)+a(3)
         endif
       endif

c If OK, save flavors of droplet and string

         ! jerr(...) removed !!!!
         
         do i=1,5
           ep(i)=re(i)
         enddo
         ic(1)=icx(1)
         ic(2)=icx(2)
         do i=1,nflav
           jc(i,1)=jcini(i,1)
           jc(i,2)=jcini(i,2)
           jcv(i,1)=jcvfin(i,1)
           jcv(i,2)=jcvfin(i,2)
         enddo

         if(ish.ge.6)then
           write(ifch,20)'droplet:',jc,ep
           write(ifch,30)'string remnant:',ic,a
           write(ifch,'(a)')'valence:'
           write(ifch,'(6i3)')jcv
         endif
 20      format(a,/,'jc:',6i3,' |',6i3,/,'ep:',5(e10.3,1x))
 30      format(a,/,'ic:',i7,' |',i7,/,'a:',5(e10.3,1x))

 1000    continue
         call utprix('getdrx',ish,ishini,4)
         end

c-----------------------------------------------------
       subroutine neworder(n1, n2, n3)
c-----------------------------------------------------
c make 3 integers ordered like 1 2 3
c------------------------------------------------------
            if(n2.lt.n1)then
              ifb=n2
              n2=n1
              n1=ifb
            endif
            if(n3.lt.n1)then
              ifb=n3
              n3=n2
              n2=n1
              n1=ifb
            elseif(n3.lt.n2)then
              ifb=n3
              n3=n2
              n2=ifb
            endif
         end

c-----------------------------------------------------
       subroutine neworderx(x1,x2,x3,i1,i2,i3)
c-----------------------------------------------------
c make 3 reals ordered like 1 2 3
c------------------------------------------------------
            if(x2.lt.x1)then
              xfb=x2
              x2=x1
              x1=xfb
              ifb=i2
              i2=i1
              i1=ifb
            endif
            if(x3.lt.x1)then
              xfb=x3
              x3=x2
              x2=x1
              x1=xfb
              ifb=i3
              i3=i2
              i2=i1
              i1=ifb
            elseif(x3.lt.x2)then
              xfb=x3
              x3=x2
              x2=xfb
              ifb=i3
              i3=i2
              i2=ifb
            endif
         end

*-- Author :    D. HECK IK FZK KARLSRUHE       27/04/1994
C=======================================================================

      SUBROUTINE EPOVAPOR( MAPRO,INEW,JFIN,ITYP,PFRX,PFRY,PFRZ )

C-----------------------------------------------------------------------
C  (E)VAPOR(ATION OF NUCLEONS AND ALPHA PARTICLES FROM FRAGMENT)
C
C  TREATES THE REMAINING UNFRAGMENTED NUCLEUS
C  EVAPORATION FOLLOWING CAMPI APPROXIMATION.
C  SEE: X. CAMPI AND J. HUEFNER, PHYS.REV. C24 (1981) 2199
C  AND  J.J. GAIMARD, THESE UNIVERSITE PARIS 7, (1990)
C  THIS SUBROUTINE IS CALLED FROM SDPM, DPMJST, NSTORE, AND VSTORE.
C  ARGUMENTS INPUT:
C   MAPRO        = NUMBER OF NUCLEONS OF PROJECTILE
C   INEW         = PARTICLE TYPE OF SPECTATOR FRAGMENT
C  ARGUMENTS OUTPUT:
C   JFIN         = NUMBER OF FRAGMENTS
C   ITYP(1:JFIN) = NATURE (PARTICLE CODE) OF FRAGMENTS
C   PFRX(1:JFIN) = TRANSVERSE MOMENTUM OF FRAGMENTS IN X-DIRECTION
C   PFRY(1:JFIN) = TRANSVERSE MOMENTUM OF FRAGMENTS IN Y-DIRECTION
C   PFRZ(1:JFIN) = LONGITUDINAL MOMENTUM OF FRAGMENTS IN Z-DIRECTION
C
C  FROM CORSIKA AND ADAPTED BY T. PIEROG TO INCLUDE LONG MOMENTUM AND
C  MORE REALISTIC FRAGMENTS
C-----------------------------------------------------------------------

#include "aaa.h"
      common/eporansto2/irndmseq
      integer irndmseq

      DOUBLE PRECISION PFR(mamxx),PFRX(mamxx),PFRY(mamxx),PFRZ(mamxx)
     *                ,RD(2*mamxx),SPFRY,SPFRZ,drangen
      DOUBLE PRECISION AFIN,AGLH,APRF,BGLH,EEX,PHIFR,RANNORM,SPFRX
      INTEGER          ITYP(mamxx),IARM,INEW,ITYPRM,INRM,IS,IZRM,JC,lseq
     *                ,JFIN,K,L,LS,MAPRO,MF,NFIN,NINTA,NNUC,NPRF,NNSTEP
      SAVE
      EXTERNAL         RANNORM
C-----------------------------------------------------------------------

      IF(ish.ge.7)WRITE(ifch,*)'EPOVAPOR : MAPRO,INEW=',MAPRO,INEW

      JC     = 0
      JFIN   = 0
      lseq   = irndmseq
      ITYPRM = INEW
      NPRF   = INEW/100
      NINTA  = MAPRO - NPRF
      IF ( NINTA .EQ. 0 ) THEN
C  NO NUCLEON HAS INTERACTED
        JFIN    = 1
        PFRX(1)  = 0.D0
        PFRY(1)  = 0.D0
        PFRZ(1)  = 0.D0
        ITYP(1) = ITYPRM
        IF(ish.ge.7)WRITE(ifch,*) 'EPOVAPOR : JFIN,NINTA=',JFIN,NINTA
        GOTO 50
      ENDIF

C  EXCITATION ENERGY EEX OF PREFRAGMENT
C  SEE: J.J. GAIMARD, THESE UNIVERSITE PARIS 7, (1990), CHPT. 4.2
      EEX = 0.D0
      CALL RMMARD( RD,2*NINTA,lseq )
      DO  L = 1, NINTA
        IF ( RD(NINTA+L) .LT. RD(L) ) RD(L) = 1.D0 - RD(L)
        EEX = EEX + RD(L)
      ENDDO
C  DEPTH OF WOODS-SAXON POTENTIAL TO FERMI SURFACE IS 0.040 GEV
      IF(ish.ge.7)WRITE(ifch,*)'EPOVAPOR : EEX=',SNGL(EEX*0.04D0),
     &                                            ' GEV'
C  EVAPORATION: EACH EVAPORATION STEP NEEDS ABOUT 0.020 GEV, THEREFORE
C  NNSTEP IS EEX * 0.04/0.02 = EEX * 2.
      NNSTEP = INT( EEX*2.D0 )

      IF ( NNSTEP .LE. 0 ) THEN
C  EXCITATION ENERGY TOO SMALL, NO EVAPORATION
        JFIN    = 1
        PFRX(1)  = 0.D0
        PFRY(1)  = 0.D0
        PFRZ(1)  = 0.D0
        ITYP(1) = ITYPRM
        IF(ish.ge.7)WRITE(ifch,*) 'EPOVAPOR : JFIN,EEX=',JFIN,SNGL(EEX)
        GOTO 50

      ENDIF

C  AFIN IS ATOMIC NUMBER OF FINAL NUCLEUS
      APRF = DBLE(NPRF)
      AFIN = APRF - 1.6D0 * DBLE(NNSTEP)
      NFIN = MAX( 0, INT( AFIN+0.5D0 ) )
C  CORRESPONDS TO DEFINITION; FRAGMENTATION-EVAPORATION
C  CONVOLUTION EMU07 /MODEL ABRASION EVAPORATION (JNC FZK APRIL 94)
C  NNUC IS NUMBER OF EVAPORATING NUCLEONS
      NNUC = NPRF - NFIN
      IF(ish.ge.7)WRITE(ifch,*) 'EPOVAPOR : NFIN,NNUC=',NFIN,NNUC

      IF     ( NNUC .LE. 0 ) THEN
C  NO EVAPORATION
        JFIN    = 1
        PFRX(1)  = 0.D0
        PFRY(1)  = 0.D0
        PFRZ(1)  = 0.D0
        ITYP(1) = ITYPRM
        GOTO 50

      ELSEIF ( NNUC .GE. 4 ) THEN
C  EVAPORATION WITH FORMATION OF ALPHA PARTICLES POSSIBLE
C  IARM, IZRM, INRM ARE NUMBER OF NUCLEONS, PROTONS, NEUTRONS OF
C  REMAINDER
        DO  LS = 1, NNSTEP
          IARM = ITYPRM/100
          IF ( IARM .LE. 0 ) GOTO 100
          IZRM = MOD(ITYPRM,100)
          INRM = IARM - IZRM
          JC   = JC + 1
          CALL RMMARD( RD,2,lseq )
          IF ( RD(1) .LT. 0.2D0  .AND.  IZRM .GE. 2
     *                           .AND.  INRM .GE. 2 ) THEN
            ITYP(JC) = -402          !alpha
            NNUC     = NNUC - 4
            ITYPRM   = ITYPRM - 402
          ELSE
            IF ( IZRM .EQ. 1 .AND. INRM .GT. IZRM ) THEN
              ITYP(JC) = 1220
              ITYPRM   = ITYPRM - 100              
            ELSEIF ( INRM .EQ. 1 .AND. IZRM .GT. INRM ) THEN
              ITYP(JC) = 1120
              ITYPRM   = ITYPRM - 101              
            ELSEIF(RD(2)*(IZRM+INRM).LT.IZRM.AND.IZRM.GE.INRM)THEN
              ITYP(JC) = 1120
              ITYPRM   = ITYPRM - 101
            ELSE
              ITYP(JC) = 1220
              ITYPRM   = ITYPRM - 100
            ENDIF
            NNUC = NNUC - 1
          ENDIF
          IF ( NNUC .LE. 0 ) GOTO 50
        ENDDO
      ENDIF

      IF ( NNUC .LT. 4 ) THEN
C  EVAPORATION WITHOUT FORMATION OF ALPHA PARTICLES
        CALL RMMARD( RD,NNUC,lseq )
        DO  IS = 1, NNUC
          IARM = ITYPRM/100
          IF ( IARM .LE. 0 ) GOTO 100
          IZRM = MOD(ITYPRM,100)
          INRM = IARM - IZRM
          JC   = JC + 1
          IF ( IZRM .EQ. 1 .AND. INRM .GT. IZRM ) THEN
            ITYP(JC) = 1220
            ITYPRM   = ITYPRM - 100              
          ELSEIF ( INRM .EQ. 1 .AND. IZRM .GT. IZRM ) THEN
            ITYP(JC) = 1120
            ITYPRM   = ITYPRM - 101              
          ELSEIF ( RD(IS)*IARM .LT. IZRM .AND. IZRM .GE. INRM ) THEN
            ITYP(JC) = 1120
            ITYPRM   = ITYPRM - 101
          ELSE
            ITYP(JC) = 1220
            ITYPRM   = ITYPRM - 100
          ENDIF
        ENDDO
      ENDIF

 50   CONTINUE
      IARM = ITYPRM/100
      IF ( IARM .LE. 0 ) GOTO 100
      IZRM = MOD(ITYPRM,100)
      INRM = IARM - IZRM
      JC = JC + 1
C CLEAN FRAGMENT TO HAVE ALWAYS AT LEAST TWO TIMES MORE N THAN P
      CALL RMMARD( RD,2,lseq )
      DO WHILE ( INRM .GT. INT( 1.15D0 * DBLE(IZRM) + RD(1) ) )
        ITYP(JC) = 1220
        ITYPRM   = ITYPRM - 100              
        IARM = ITYPRM/100
        IZRM = MOD(ITYPRM,100)
        INRM = IARM - IZRM
        JC = JC + 1
      ENDDO
C CLEAN FRAGMENT NOT TO HAVE TOO MANY P
      DO WHILE ( IZRM .GE. NINT( (DBLE(INRM)+1.D0+RD(2)) / 1.15D0 ) )
        ITYP(JC) = 1120
        ITYPRM   = ITYPRM - 101              
        IARM = ITYPRM/100
        IZRM = MOD(ITYPRM,100)
        INRM = IARM - IZRM
        JC = JC + 1
      ENDDO
      IF ( IARM .EQ. 8 ) THEN     !EXCLUDED
        DO WHILE ( IZRM .GE. 2 .AND. INRM .GE. 2 )
          ITYP(JC) = -402
          ITYPRM   = ITYPRM - 402              
          IARM = ITYPRM/100
          IZRM = MOD(ITYPRM,100)
          INRM = IARM - IZRM
          JC = JC + 1
        ENDDO
        IF ( ITYPRM .GT. 0 ) THEN
          IF ( ITYPRM/100 .EQ. 5 ) THEN !EXCLUDED BUT SHOULD NOT HAPPEN HERE
            IF ( IZRM .GE. INRM ) THEN
              ITYP(JC) = 1120
              ITYPRM   = ITYPRM - 101              
            ELSE
              ITYP(JC) = 1220
              ITYPRM   = ITYPRM - 100
            ENDIF
            JC = JC + 1
          ENDIF
          ITYP(JC) = -ITYPRM
        ELSE
          JC = JC - 1
        ENDIF
      ELSEIF ( IARM .EQ. 5 ) THEN     !EXCLUDED
        IF ( IZRM .GE. INRM ) THEN
          ITYP(JC) = 1120
          ITYPRM   = ITYPRM - 101              
        ELSE
          ITYP(JC) = 1220
          ITYPRM   = ITYPRM - 100
        ENDIF
        JC = JC + 1
        ITYP(JC) = -ITYPRM
      ELSEIF     ( ITYPRM .GT. 200 ) THEN
        ITYP(JC) = -ITYPRM
      ELSEIF ( ITYPRM .EQ. 101 ) THEN
        ITYP(JC) = 1120
      ELSEIF ( ITYPRM .EQ. 100 ) THEN
        ITYP(JC) = 1220
      ELSE
        JC = JC - 1
        IF ( ITYPRM .NE. 0 ) WRITE(*,*)
     *                  'EPOVAPOR : ILLEGAL PARTICLE ITYPRM =',ITYPRM
      ENDIF

  100 CONTINUE
      IF ( JC .GT. JFIN ) THEN
       JFIN = JC
       IF(ish.ge.7)WRITE(ifch,*) 
     *   'EPOVAPOR :  NO        ITYP     PFR       PFL'
       IF ( infragm .EQ. 2 ) THEN
C  EVAPORATION WITH PT AFTER PARAMETRIZED JACEE DATA
        DO  MF = 1, JFIN
          IF(ITYP(MF).LT.0)THEN
            IARM=-ITYP(MF)/100
          ELSE
            IARM=1
          ENDIF
          PFR(MF) = RANNORM(0.088D0,0.044D0)
          PFRZ(MF)= (2*int(drangen(PFR(MF))+0.5d0)-1)
     &   *RANNORM(0.300D0/DBLE(IARM),0.100D0/SQRT(DBLE(IARM)))    !Fermi motion about 300 MeV
          IF(ish.ge.7)WRITE(ifch,*) MF,ITYP(MF),SNGL(PFR(MF))
     &                                         ,SNGL(PFRZ(MF))
        ENDDO
       ELSEIF ( infragm .EQ. 3 ) THEN
C  EVAPORATION WITH PT AFTER GOLDHABER''S MODEL (PHYS.LETT.53B(1974)306)
        DO  MF = 1, JFIN
          K    = MAX( 1, -ITYP(MF)/100 )
          BGLH = K * (MAPRO - K) / DBLE(MAPRO-1)
C  THE VALUE 0.103 [GEV] IS SIGMA(0)=P(FERMI)/SQRT(5.)
*         AGLH = 0.103D0 * SQRT( BGLH )
C  THE VALUE 0.090 [GEV] IS EXPERIMENTALLY DETERMINED SIGMA(0)
          AGLH = 0.090D0 * SQRT( BGLH )
          PFR(MF) = RANNORM(0.D0,AGLH)
          PFRZ(MF)= RANNORM(0.000D0,0.500D0)    !from pAg at 100 GeV
          IF(ish.ge.7)WRITE(ifch,*) MF,ITYP(MF),SNGL(PFR(MF))
     &                                         ,SNGL(PFRZ(MF))
        ENDDO
       ELSE
C  EVAPORATION WITHOUT TRANSVERSE MOMENTUM
        DO  MF = 1, JFIN
          PFR(MF) = 0.D0
          PFRZ(MF)= 0.D0
          IF(ish.ge.7)WRITE(ifch,*) MF,ITYP(MF),SNGL(PFR(MF))
     &                                         ,SNGL(PFRZ(MF))
        ENDDO
       ENDIF
C  CALCULATE RESIDUAL TRANSVERSE MOMENTUM
       SPFRX = 0.D0
       SPFRY = 0.D0
       SPFRZ = 0.D0
       CALL RMMARD( RD,JFIN,lseq )
       DO  MF = 1, JFIN
        PHIFR = PI * RD(MF)
        PFRX(MF) = PFR(MF) * COS( PHIFR )
        PFRY(MF) = PFR(MF) * SIN( PHIFR )
        SPFRY = SPFRY + PFRY(MF)
        SPFRX = SPFRX + PFRX(MF)
        SPFRZ = SPFRZ + PFRZ(MF)
       ENDDO
C  CORRECT ALL TRANSVERSE MOMENTA FOR MOMENTUM CONSERVATION
       SPFRX = SPFRX / JFIN
       SPFRY = SPFRY / JFIN
       SPFRZ = SPFRZ / JFIN
       DO  MF = 1, JFIN
        PFRX(MF) = PFRX(MF) - SPFRX
        PFRY(MF) = PFRY(MF) - SPFRY
        PFRZ(MF) = PFRZ(MF) - SPFRZ
       ENDDO
      ENDIF

      IF(ish.ge.7)WRITE(ifch,*) 'EPOVAPOR : NINTA,JFIN=',NINTA,JFIN

      RETURN
      END

*-- Author :    The CORSIKA development group   21/04/1994
C=======================================================================

      DOUBLE PRECISION FUNCTION RANNORM( A,B )

C-----------------------------------------------------------------------
C  RAN(DOM NUMBER) NOR(MALLY DISTRIBUTED)
C
C  GENERATES NORMAL DISTRIBUTED RANDOM NUMBER
C  DELIVERS 2 UNCORRELATED RANDOM NUMBERS,
C  THEREFORE RANDOM CALLS ARE ONLY NECESSARY EVERY SECOND TIME.
c  but to be used with CONEX/CORSIKA we should always use 2 new number
c  to be able to reproduce the shower
C  REFERENCE : NUMERICAL RECIPES, W.H. PRESS ET AL.,
C              CAMBRIDGE UNIVERSITY PRESS, 1992  ISBN 0 521 43064 X
C  ARGUMENTS:
C   A      = MEAN VALUE
C   B      = STANDARD DEVIATION
C
C  FROM CORSIKA
C-----------------------------------------------------------------------

      IMPLICIT NONE
      double precision facrdm,u1rdm,u2rdm,drangen
c      logical knordm
c      data knordm/.true./

      DOUBLE PRECISION A,B,RR
      SAVE facrdm,u1rdm,u2rdm!,knordm
C-----------------------------------------------------------------------

c      IF ( KNORdm ) THEN
  1     CONTINUE
        U1rdm = 2.D0*drangen(a) - 1.D0
        U2rdm = 2.D0*drangen(b) - 1.D0
        RR = U1rdm**2 + U2rdm**2
        IF ( RR .GE. 1.D0  .OR.  RR .EQ. 0.D0 ) GOTO 1
        FACrdm = SQRT( (-2.D0) * LOG(RR) / RR )

        RANNORM = FACrdm * U1rdm * B + A
c        KNORdm   = .FALSE.
c      ELSE
c        RANNORM = FACrdm * U2rdm * B + A
c        KNORdm   = .TRUE.
c      ENDIF

      RETURN
      END



