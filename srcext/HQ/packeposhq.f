      subroutine hqinitransfer 
      use PBGtransfer
      iftranstimestep=.false.
      iftranscTc=.false.
      iftransXsectioncranck=.false.
      iftransratemult=.false.
      iftransc2=.false.
      iftransfpmult=.false.
      iftransalphasforrad=.false.
      iftranshighptcut=.false.
      iftransnbcol=.false.
      iftransnbcolsaveEPOS=.false.
      iftransnbsubrun=.false.
      iftranstype_evol_Q=.false.
      iftransimedfinal=.false.
      iftransreddof=.false.
      iftransitypDproduct=.false.
      iftransitypfrag=.false.
      iftransitypcoal=.false.
      iftranstempcoal=.false.
      iftransifcoalinflcell=.false.
      iftransboltz_model=.false.
      iftranstypelpmdamp=.false.
      iftransfpchoice=.false.
      iftranstunchoice=.false.
      iftransifdecrease=.false.
      iftranscoll=.false.
      iftransifcronin=.false.
      iftransifcol=.false.
      iftransifrad=.false.
      iftransiffullQCD=.false.
      iftransiflpm=.false.
      iftransifdamp=.false.
      iftransifptcles=.false.
      iftransifhighptselect=.false.
      iftransifbselect=.false.
      iftransish=.false.
      iftransifhqfromhq=.false.
      iftransifextspectra1=.false.
      iftranstypeextspectra1=.false.
      iftransifshadowing1=.false.
      iftranstypeshadow1=.false.
      iftransifextspectra2=.false.
      iftranstypeextspectra2=.false.
      iftransifshadowing2=.false.
      iftranstypeshadow2=.false.
      iftransifhqskip=.false.
      iftrans2bodySC=.false.
      iftranssingletred=.false.
      iftransdpmaxabscorona=.false.
      iftranstypedistrel=.false.
      iftransparamdistrelc=.false.
      iftransparamdistrelb=.false.
      ifinitihq2=.false.
      ifinitiskip=.false.
      ifinitivirtual2=.false.
      ifinitidmeson=.false.
      end

      subroutine hqoptns(line,j)
      use PBGtransfer
      use PBGifskip
      use PBGDchoice
      use PBGifhq2
      use PBGcivirtual2
      character*1000 line,linex
      double precision val
      integer ihq
      common/ifhq/ihq
      integer ivirtual
      common/civirtual/ivirtual
    1 call utword(line,i,j,1)
      if(line(i:j).eq.'endoptns')then
         call utword(line,i,j,0)
         if(line(i:j).eq.'hq') return
         if(line(i:j).ne.'hq')stop 'endoptns mismatch!'
      elseif(line(i:j).eq.'set')then
         call utword(line,i,j,0)
         linex=line
         ix=i
         jx=j
         call utword(line,i,j,0)
         read(line(i:j),*)val
         if(linex(ix:jx).eq.'ihq') then
            ihq2=nint(val)
            ihq=ihq2
            ifinitihq2=.true.
         endif
         if(linex(ix:jx).eq.'hqifallowhqskip') then
            hqifallowhqskip=nint(val)
            iftransifhqskip=.true.
         endif
         if(linex(ix:jx).eq.'hqif2bodySC') then
            hqif2bodySC=nint(val)
            iftrans2bodySC=.true.
         endif
         if(linex(ix:jx).eq.'ifsingletred') then
            hqifsingletred=nint(val)
            iftranssingletred=.true.
         endif
         if(linex(ix:jx).eq.'hqdpmaxabscorona') then
            hqdpmaxabscorona=val
            iftransdpmaxabscorona=.true.
         endif
         if(linex(ix:jx).eq.'hqtypedistrel') then 
            hqtypedistrel=nint(val)
            iftranstypedistrel=.true.
         endif
         if(linex(ix:jx).eq.'hqparamdistrelc') then 
            hqparamdistrelc=val
            iftransparamdistrelc=.true.
         endif
         if(linex(ix:jx).eq.'hqparamdistrelb') then 
            hqparamdistrelb=val
            iftransparamdistrelb=.true.
         endif
         if(linex(ix:jx).eq.'iskip') then
            iskip=nint(val)
            ifinitiskip=.true.
         endif
         if(linex(ix:jx).eq.'idmeson') then
            idmeson=nint(val)
            ifinitidmeson=.true.
         endif
         if(linex(ix:jx).eq.'ivirtual') then
            ivirtual2=nint(val)   
            ivirtual=ivirtual2
            ifinitivirtual2=.true.
         endif
         if(linex(ix:jx).eq.'hqnbcol') then 
            hqnbcol=nint(val)
            iftransnbcol=.true.
         endif
         if(linex(ix:jx).eq.'hqnbcolsaveEPOS') then 
            hqnbcolsaveEPOS=nint(val)
            iftransnbcolsaveEPOS=.true.
         endif
         if(linex(ix:jx).eq.'hqnbsubrun') then
            hqnbsubrun=int(val)
            iftransnbsubrun=.true.
         endif
         if(linex(ix:jx).eq.'hqtype_evol_q') then
            hqtype_evol_Q=nint(val)
            iftranstype_evol_Q=.true.
         endif
         if(linex(ix:jx).eq.'hqtimestep') then
            hqtimestep=val
            iftranstimestep=.true.
         endif
         if(linex(ix:jx).eq.'hqimedfinal') then
            hqimedfinal=nint(val) 
            iftransimedfinal=.true.
         endif
         if(linex(ix:jx).eq.'hqifdecrease') then
            hqifdecrease=nint(val)
            iftransifdecrease=.true.
         endif
         if(linex(ix:jx).eq.'hqifextspectrac') then
            hqifextspectra1=nint(val) 
            iftransifextspectra1=.true.
         endif
         if(linex(ix:jx).eq.'hqifextspectrab') then
            hqifextspectra2=nint(val) 
            iftransifextspectra2=.true.
         endif
         if(linex(ix:jx).eq.'hqtypeextspectrac') then
            hqtypeextspectra1=nint(val) 
            iftranstypeextspectra1=.true.
         endif
         if(linex(ix:jx).eq.'hqtypeextspectrab') then
            hqtypeextspectra2=nint(val) 
            iftranstypeextspectra2=.true.
         endif
         if(linex(ix:jx).eq.'hqifshadowingc') then
            hqifshadowing1=nint(val) 
            iftransifshadowing1=.true.
         endif
         if(linex(ix:jx).eq.'hqifshadowingb') then
            hqifshadowing2=nint(val) 
            iftransifshadowing2=.true.
         endif
         if(linex(ix:jx).eq.'hqtypeshadowc') then
            hqtypeshadow1=nint(val) 
            iftranstypeshadow1=.true.
         endif
         if(linex(ix:jx).eq.'hqtypeshadowb') then
            hqtypeshadow2=nint(val) 
            iftranstypeshadow2=.true.
         endif
         if(linex(ix:jx).eq.'hqifcronin') then 
            hqifcronin=nint(val) 
            iftransifcronin=.true.
         endif
         if(linex(ix:jx).eq.'hqifcronin') then 
            hqifcronin=nint(val) 
            iftransifcronin=.true.
         endif
         if(linex(ix:jx).eq.'hqreddof') then
            hqreddof=nint(val)
            iftransreddof=.true.
         endif
         if(linex(ix:jx).eq.'hqctc') then 
            hqcTc=val
            iftranscTc=.true.
         endif
         if(linex(ix:jx).eq.'hqcoll') then
            hqcoll=nint(val)
            iftranscoll=.true.
         endif
         if(linex(ix:jx).eq.'hqish') then
            hqish=nint(val)
            iftransish=.true.
         endif
         if(linex(ix:jx).eq.'hqhighptselect') then
            hqifhighptselect=nint(val)
            iftransifhighptselect=.true.
         endif
         if(linex(ix:jx).eq.'hqbselect') then
            hqifbselect=nint(val)
            iftransifbselect=.true.
         endif
         if(linex(ix:jx).eq.'hqhighptcut') then 
            hqhighptcut=val
            iftranshighptcut=.true.
         endif
         if(linex(ix:jx).eq.'hqifhqfromhq') then
            hqifhqfromhq=nint(val)
            iftransifhqfromhq=.true.
         endif
         if(linex(ix:jx).eq.'hqifptcles') then 
            hqifptcles=nint(val)
            iftransifptcles=.true.
         endif
         if(linex(ix:jx).eq.'hqxsectioncranck')then 
            hqXsectioncranck=nint(val)
            iftransXsectioncranck=.true.
         endif
         if(linex(ix:jx).eq.'hqitypdproduct')then  
            hqitypDproduct=nint(val)
            iftransitypDproduct=.true.
         endif
         if(linex(ix:jx).eq.'hqitypfrag')then 
            hqitypfrag=nint(val)
            iftransitypfrag=.true.
         endif
         if(linex(ix:jx).eq.'hqitypcoal')then 
            hqitypcoal=nint(val)
            iftransitypcoal=.true.
         endif
         if(linex(ix:jx).eq.'hqifcoalinflcell')then 
            hqifcoalinflcell=nint(val)
            iftransifcoalinflcell=.true.
         endif
         if(linex(ix:jx).eq.'hqtempcoal')then 
            hqtempcoal=val
            iftranstempcoal=.true.
         endif
         if(linex(ix:jx).eq.'hqifcol')then 
            hqifcol=nint(val)
            iftransifcol=.true.
         endif
         if(linex(ix:jx).eq.'hqboltz_model')then 
            hqboltz_model=nint(val)
            iftransboltz_model=.true.
         endif
         if(linex(ix:jx).eq.'hqratemult')then 
            hqratemult=val
            iftransratemult=.true.
         endif
         if(linex(ix:jx).eq.'hqifrad')then 
            hqifrad=nint(val)
            iftransifrad=.true.
         endif
         if(linex(ix:jx).eq.'hqfpchoice')then 
            hqfpchoice=nint(val)
            iftransfpchoice=.true.
         endif
         if(linex(ix:jx).eq.'hqtunchoice')then 
            hqtunchoice=nint(val)
            iftranstunchoice=.true.
         endif
         if(linex(ix:jx).eq.'hqfpmult')then 
            hqfpmult=val
            iftransfpmult=.true.
         endif
         if(linex(ix:jx).eq.'hqalphasforrad')then 
            hqalphasforrad=val
            iftransalphasforrad=.true.
         endif
         if(linex(ix:jx).eq.'hqiffullqcd')then 
            hqiffullQCD=nint(val)
            iftransiffullQCD=.true.
         endif
         if(linex(ix:jx).eq.'hqtypelpmdamp')then 
            hqtypelpmdamp=nint(val)
            iftranstypelpmdamp=.true.
         endif
         if(linex(ix:jx).eq.'hqiflpm')then 
            hqiflpm=nint(val)
            iftransiflpm=.true.
         endif
         if(linex(ix:jx).eq.'hqifdamp')then 
            hqifdamp=nint(val)
            iftransifdamp=.true.
         endif
         if(linex(ix:jx).eq.'hqc2')then 
            hqc2=val
            iftransc2=.true.
         endif
      elseif(line(i:j).eq.'#if1')then
         call  setIf1(line,i,j)
      elseif(line(i:j).eq.'#fi')then
         continue
      else
         write(ifmt,'(a,a,a)')'hqoptns "',line(i:j),'" not found'
         stop
      endif
      i=j+1
      goto 1
      end

      subroutine getihq(i) 
      use PBGifhq2
      implicit none
      integer i
      i=ihq2
      end 

      subroutine getivirtual(i) 
      use PBGcivirtual2
      implicit none
      integer i
      i=ivirtual2
      end 

      subroutine hqini
      use PBGtransfer
      integer ihq
      common/ifhq/ihq
#include "aaa.h" 

CJA end added 24/11/16
      if(ihq.ne.1)return
CJA added 24/11
      if(hqcoll.eq.1)open(unit=73,file=fnus2(1:nfnus2)//'5',
     &status='unknown')

CJA end added 24/11
      call readin             
      call big_initialization 
      call hqinix
      call big_initialization2   
      end

      subroutine hqinix
      use JA
      use JA2
      use PBGDchoice
      use PBGifhq2
      use PBGcivirtual2
      use PBGhqprod
      use PBGmainpro1
      use PBGmainpro2
      use PBGfragandcoal
      use PBGgenvar
      use PBGtransfer
      use PBGspectra
      use PBGevolQ
      use PBGreductiondofs
      use PBGXsection
      use PBGforratedat1
      use PBGforprocesstype
      use PBGforboltzmann
      use PBGforradiatGB
      use PBGforradiatLPM
      use PBGgluprops
      use PBGcoupling
      use PBGfpcoeff
      use PBGforfpdat1
      use PBGblckifptcles
      use PBGdistccbarbis
      use PBGtrigger
      use PBGifskip
      use for2body
      integer ifmtx

      call getMonitorFileIndex(ifmtx)
      if(iftransnbcol) then
         nbcol=hqnbcol
         ifinitnbcol=.true.   
      endif
      if(iftransnbcolsaveEPOS) then
         nbcolsaveEPOS=hqnbcolsaveEPOS
         ifinitnbcolsaveEPOS=.true.   
      endif
      if(iftransnbsubrun)then
         nbsubrun=hqnbsubrun
         ifinitnbsubrun=.true.
      endif
      if(iftranstype_evol_Q)then
         type_evol_Q=hqtype_evol_Q
         ifinittype_evol_Q=.true.
      endif
      if(iftranstimestep)then
         timestep=hqtimestep
         ifinittimestep=.true.
      endif
      if(iftransimedfinal)then
         imedfinal=hqimedfinal
         ifinitimedfinal=.true.
      endif
      if(iftransifhighptselect)then
         if(hqifhighptselect.eq.1) then
            ifhighptselect=.true.
         else
            ifhighptselect=.false.
         endif      
         ifinitifhighptselect=.true.
      endif
      if(iftransifbselect)then 
         if(hqifbselect.eq.1) then
             ifbselect=.true.
          else
             ifbselect=.false.
          endif
          ifinitifbselect=.true.
      endif
      if(iftranshighptcut)then
         highptcut=hqhighptcut
         ifinithighptcut=.true.
      endif
      if(iftransifhqfromhq)then
         if(hqifhqfromhq.eq.1) then
            ifhqfromhq=.true.
            write(ifmtx,*) '*******************************************'
            write(ifmtx,*) 'HQ PRODUCTION FROM MC@HQ'
            write(ifmtx,*) '*******************************************'
         else
            ifhqfromhq=.false.
            write(ifmtx,*) '*******************************************'
            write(ifmtx,*) 'HQ PRODUCTION FROM EPOS'
            write(ifmtx,*) '*******************************************'
         endif      
         ifinitifhqfromhq=.true.
      endif
      if(iftransifhqskip)then
         if(hqifallowhqskip.eq.1) then
            ifallowhqskip=.true.
         else
            ifallowhqskip=.false.
         endif
         ifinithqskip=.true.
      endif
      if(iftrans2bodySC)then
         if2bodySC=hqif2bodySC
         ifinit2bodySC=.true.
      endif
      if(iftranssingletred)then
         if(hqifsingletred.eq.1) then
            ifsingletred=.true.
         else
            ifsingletred=.false.
         endif
         ifinitsingletred=.true.
      endif
      if(iftransdpmaxabscorona)then
         dpmaxabscorona=hqdpmaxabscorona
         ifinitdpmaxabscorona=.true.
      endif
      if(iftranstypedistrel)then
         typedistrel=hqtypedistrel
         ifinittypedistrel=.true.
      endif
      if(iftransparamdistrelc)then
         paramdistrel(1)=hqparamdistrelc
         ifinitparamdistrelc=.true.
      endif
      if(iftransparamdistrelb)then
         paramdistrel(2)=hqparamdistrelb
         ifinitparamdistrelb=.true.
      endif
      if(iftransifdecrease)then
         if(hqifdecrease.eq.1) then
            ifdecrease=.true.
         else
            ifdecrease=.false.
         endif
         ifinitifdecrease=.true.
      endif
      if(iftransifextspectra1)then
         if(hqifextspectra1.eq.1) then
            ifextspectra(1)=.true.
         else
            ifextspectra(1)=.false.
         endif
         ifinitifextspectra1=.true.
      endif
      if(iftransifextspectra2)then
         if(hqifextspectra2.eq.1) then
            ifextspectra(2)=.true.
         else
            ifextspectra(2)=.false.
         endif
         ifinitifextspectra2=.true.
      endif
      if(iftranstypeextspectra1)then
         typeextspectra(1)=hqtypeextspectra1
         ifinittypeextspectra1=.true.
      endif
      if(iftranstypeextspectra2)then
         typeextspectra(2)=hqtypeextspectra2
         ifinittypeextspectra2=.true.
      endif
      if(iftransifshadowing1)then
         if(hqifshadowing1.eq.1) then
            ifshadowing(1)=.true.
         else
            ifshadowing(1)=.false.
         endif
         ifinitifshadowing1=.true.
      endif
      if(iftransifshadowing2)then
         if(hqifshadowing2.eq.1) then
            ifshadowing(2)=.true.
         else
            ifshadowing(2)=.false.
         endif
         ifinitifshadowing2=.true.
      endif
      if(iftranstypeshadow1)then
         typeshadow(1)=hqtypeshadow1
         ifinittypeshadow1=.true.
      endif
      if(iftranstypeshadow2)then
         typeshadow(2)=hqtypeshadow2
         ifinittypeshadow2=.true.
      endif
      if(iftransifcronin)then
         if(hqifcronin.eq.1) then
            ifcronin=.true.
         else
            ifcronin=.false.
         endif
         ifinitifcronin=.true.
      endif
      if(iftransreddof)then
         reddof=hqreddof
         ifinitreddof=.true.
      endif
      if(iftranscTc)then
         cTc=hqcTc
         ifinitcTc=.true.
      endif
      if(iftransXsectioncranck)then
         Xsectioncranck=hqXsectioncranck
         ifinitXsectioncranck=.true.
      endif
      if(iftransitypDproduct)then
         itypDproduct=hqitypDproduct
         ifinititypDproduct=.true.
      endif
      if(iftransitypfrag)then
         itypfrag=hqitypfrag
         ifinititypfrag=.true.
      endif
      if(iftransitypcoal)then
         itypcoal=hqitypcoal
         ifinititypcoal=.true.
      endif
      if(iftransifcoalinflcell)then
         if(hqifcoalinflcell.eq.1) then
            ifcoalinflcell=.true.
         else
            ifcoalinflcell=.false.
         endif
         ifinitifcoalinflcell=.true.
      endif
      if(iftranstempcoal)then
         tempcoal=hqtempcoal
         ifinittempcoal=.true.
      endif
      if(iftransifcol)then
         if(hqifcol.eq.1) then
            ifcol=.true.
         else
            ifcol=.false.
         endif
         ifinitifcol=.true.
      endif
      if(iftransboltz_model)then
         boltz_model=hqboltz_model
         ifinitboltz_model=.true.
      endif
      if(iftransratemult)then
         ratemult=hqratemult
         ifinitratemult=.true.
      endif
      if(iftransfpchoice)then
         fpchoice=hqfpchoice
         ifinitfpchoice=.true.
      endif
      if(iftranstunchoice)then
         tunchoice=hqtunchoice
         ifinittunchoice=.true.
      endif
      if(iftransfpmult)then
         fpmult=hqfpmult
         ifinitfpmult=.true.
      endif
      if(iftransifrad)then
         if(hqifrad.eq.1) then
            ifrad=.true.
         else
            ifrad=.false.
         endif
         ifinitifrad=.true.
      endif
      if(iftransalphasforrad)then
         alphasforrad=hqalphasforrad
         ifinitalphasforrad=.true.
      endif
      if(iftransiffullQCD)then
         if(hqiffullQCD.eq.1) then
            iffullQCD=.true.
         else
            iffullQCD=.false.
         endif
         ifinitiffullQCD=.true.
      endif
      if(iftranstypelpmdamp)then
         typelpmdamp=hqtypelpmdamp
         ifinittypelpmdamp=.true.
      endif
      if(iftransiflpm)then
         if(hqiflpm.eq.1) then
            iflpm=.true.
         else
            iflpm=.false.
         endif
         ifinitiflpm=.true.
      endif
      if(iftransifdamp)then
         if(hqifdamp.eq.1) then
            ifdamp=.true.
         else
            ifdamp=.false.
         endif
         ifinitifdamp=.true.
      endif
      if(iftransc2)then
         c2=hqc2
         ifinitc2=.true.
      endif
      if(iftransifptcles)then
         if(hqifptcles.eq.1) then
            ifptcles=.true.
         else
            ifptcles=.false.
         endif
         ifinitifptcles=.true.
      endif
      if(iftranscoll)then
         jahqcoll=hqcoll
         ifinitcoll=.true.
      endif
      if(iftransish)then
         jaish=hqish
         ifinitish=.true.
      endif

      write(ifmtx,*) '*************************************************'
      write(ifmtx,*) '**                 HQ variable                 **'
      if(ifinitihq2) then
         write(ifmtx,*)  'ihq :',ihq2
      else
         write(ifmtx,*)  'ihq ',' UNDEFINED'
         close(ifmtx)
         stop
      endif
      if(ifinitiskip) then
         write(ifmtx,*)  'iskip :',iskip
      else
         write(ifmtx,*)  'iskip ',' UNDEFINED'
         close(ifmtx)
         stop
      endif
      if(ifinitidmeson) then
         write(ifmtx,*)  'idmeson :',idmeson
      else
         write(ifmtx,*)  'idmeson ',' UNDEFINED'
         close(ifmtx)
         stop
      endif
      if(ifinitivirtual2) then
         write(ifmtx,*)  'ivirtual :',ivirtual2
      else
         write(ifmtx,*)  'ivirtual ',' UNDEFINED'
         close(ifmtx)
         stop
      endif
      if(ifinittimestep) then
         write(ifmtx,*)  'timestep :',timestep
      else
         write(ifmtx,*)  'timestep ',' UNDEFINED'
         close(ifmtx)
         stop
      endif
      if(ifinitcTc) then
         write(ifmtx,*) 'cTc = T_max/Tc:',cTc
      else
         write(ifmtx,*) 'cTc ','UNDEFINED'
         close(ifmtx)
         stop
      endif
      if(ifinitXsectioncranck) then
         write(ifmtx,*) 'Xsectioncranck:',Xsectioncranck
      else
         write(ifmtx,*) 'Xsectioncranck:','UNDEFINED'
         close(ifmtx)
         stop
      endif
      if(ifinitratemult) then
         write(ifmtx,*) 'ratemult:',ratemult
      else
         write(ifmtx,*) 'ratemult:','UNDEFINED'
         close(ifmtx)
         stop
      endif
      if(ifinitc2) then
         write(ifmtx,*)  'c2:',c2
      else
         write(ifmtx,*) 'c2:','UNDEFINED'
         close(ifmtx)
         stop
      endif
      if(ifinitfpmult) then
         write(ifmtx,*) 'fpmult:',fpmult
      else
         write(ifmtx,*) 'fpmult:','UNDEFINED'
         close(ifmtx)
         stop
      endif
      if(ifinitalphasforrad) then
         write(ifmtx,*) 'alphasforrad:',alphasforrad
      else
         write(ifmtx,*) 'alphasforrad:','UNDEFINED'
         close(ifmtx)
         stop
      endif
      if(ifinitifhqfromhq) then
         write(ifmtx,*) 'ifhqfromhq:',ifhqfromhq
      else
         write(ifmtx,*) 'ifhqfromhq:','UNDEFINED'
         close(ifmtx)
         stop
      endif
      if(ifinitifextspectra1) then
         write(ifmtx,*) 'ifextspectra(1):',ifextspectra(1)
      else
         write(ifmtx,*) 'ifextspectra(1):','UNDEFINED'
         close(ifmtx)
         stop
      endif
      if(ifinitifextspectra2) then
         write(ifmtx,*) 'ifextspectra(2):',ifextspectra(2)
      else
         write(ifmtx,*) 'ifextspectra(2):','UNDEFINED'
         close(ifmtx)
         stop
      endif
      if(ifinittypeextspectra1) then
         write(ifmtx,*) 'typeextspectra(1):',typeextspectra(1)
      else
         write(ifmtx,*) 'typeextspectra(1):','UNDEFINED'
         close(ifmtx)
         stop
      endif
      if(ifinittypeextspectra2) then
         write(ifmtx,*) 'typeextspectra(2):',typeextspectra(2)
      else
         write(ifmtx,*) 'typeextspectra(2):','UNDEFINED'
         close(ifmtx)
         stop
      endif
      if(ifinitifshadowing1) then
         write(ifmtx,*) 'ifshadowing(1):',ifshadowing(1)
      else
         write(ifmtx,*) 'ifshadowing(1):','UNDEFINED'
         close(ifmtx)
         stop
      endif
      if(ifinitifshadowing2) then
         write(ifmtx,*) 'ifshadowing(2):',ifshadowing(2)
      else
         write(ifmtx,*) 'ifshadowing(2):','UNDEFINED'
         close(ifmtx)
         stop
      endif
      if(ifinittypeshadow1) then
         write(ifmtx,*) 'typeshadow(1):',typeshadow(1)
      else
         write(ifmtx,*) 'typeshadow(1):','UNDEFINED'
         close(ifmtx)
         stop
      endif
      if(ifinittypeshadow2) then
         write(ifmtx,*) 'typeshadow(2):',typeshadow(2)
      else
         write(ifmtx,*) 'typeshadow(2):','UNDEFINED'
         close(ifmtx)
         stop
      endif
      if(ifinitifcronin) then
         write(ifmtx,*) 'ifcronin:',ifcronin
      else
         write(ifmtx,*) 'ifcronin:','UNDEFINED'
         close(ifmtx)
         stop
      endif
      if(ifinitifhighptselect) then
         write(ifmtx,*) 'highptselect:',ifhighptselect
      else
         write(ifmtx,*) 'highptselect:','UNDEFINED'
         close(ifmtx)
         stop
      endif
       if(ifinitifbselect) then
         write(ifmtx,*) 'bselect:',ifbselect
       else
          write(ifmtx,*) 'bselect:','UNDEFINED'
          stop
      endif 
      if(ifinithighptcut) then
         write(ifmtx,*) 'highptcut:',highptcut
      else
         write(ifmtx,*) 'highptcut:','UNDEFINED'
         close(ifmtx)
         stop
      endif
      if(ifinitnbcol) then
         write(ifmtx,*) 'nbcol:',nbcol
      else
         write(ifmtx,*) 'nbcol:','UNDEFINED'
         close(ifmtx)
         stop
      endif
      if(ifinitnbcolsaveEPOS) then
         write(ifmtx,*) 'nbcolsaveEPOS:',nbcolsaveEPOS
      else
         write(ifmtx,*) 'nbcolsaveEPOS:','UNDEFINED'
         close(ifmtx)
         stop
      endif
      if(ifinithqskip) then
         write(ifmtx,*) 'ifallowhqskip:',ifallowhqskip
      else
         write(ifmtx,*) 'ifallowhqskip:','UNDEFINED'
         close(ifmtx)
         stop
      endif
      if(ifinit2bodySC) then
         write(ifmtx,*) 'if2bodySC:',if2bodySC
      else
         write(ifmtx,*) 'if2bodySC:','UNDEFINED'
         close(ifmtx)
         stop
      endif
      if(ifinitsingletred) then
         write(ifmtx,*) 'ifsingletred:',ifsingletred
      else
         write(ifmtx,*) 'ifsingletred:','UNDEFINED'
         close(ifmtx)
         stop
      endif
      if(ifinitdpmaxabscorona) then
         write(ifmtx,*) 'dpmaxabscorona:',dpmaxabscorona
      else
         write(ifmtx,*) 'dpmaxabscorona:','UNDEFINED'
         close(ifmtx)
         stop
      endif
      if(ifinittypedistrel) then
         write(ifmtx,*) 'typedistrel:',typedistrel
      else
         write(ifmtx,*) 'typedistrel:','UNDEFINED'
         close(ifmtx)
         stop
      endif
      if(ifinitparamdistrelc) then
         write(ifmtx,*) 'paramdistrel(1):',paramdistrel(1)
      else
         write(ifmtx,*) 'paramdistrel(1):','UNDEFINED'
         close(ifmtx)
         stop
      endif
      if(ifinitparamdistrelb) then
         write(ifmtx,*) 'paramdistrel(2):',paramdistrel(2)
      else
         write(ifmtx,*) 'paramdistrel(2):','UNDEFINED'
         close(ifmtx)
         stop
      endif
      if(ifinitnbsubrun) then
         write(ifmtx,*) 'nbsubrun:',nbsubrun
      else
         write(ifmtx,*) 'nbsubrun:','UNDEFINED'
         close(ifmtx)
         stop
      endif
      if(ifinittype_evol_Q) then
         write(ifmtx,*) 'type_evol_Q:',type_evol_Q
      else
         write(ifmtx,*) 'type_evol_Q:','UNDEFINED'
         close(ifmtx)
         stop
      endif
      if(ifinitimedfinal) then
         write(ifmtx,*) 'imedfinal:',imedfinal
      else
         write(ifmtx,*) 'imedfinal:','UNDEFINED'
         close(ifmtx)
         stop
      endif
      if(ifinitreddof) then
         write(ifmtx,*) 'ifinitreddof:',ifinitreddof
      else
         write(ifmtx,*) 'ifinitreddof:','UNDEFINED'
         close(ifmtx)
         stop
      endif
      if(ifinititypDproduct) then
         write(ifmtx,*) 'itypDproduct:',itypDproduct
      else
         write(ifmtx,*) 'itypDproduct:','UNDEFINED'
         close(ifmtx)
         stop
      endif
      if(ifinititypfrag) then
         write(ifmtx,*) 'itypfrag:',itypfrag
      else
         write(ifmtx,*) 'itypfrag:','UNDEFINED'
         close(ifmtx)
         stop
      endif
      if(ifinititypcoal) then
         write(ifmtx,*) 'itypcoal:',itypcoal
      else
         write(ifmtx,*) 'itypcoal:','UNDEFINED'
         close(ifmtx)
         stop
      endif
      if(ifinitifcoalinflcell) then
         write(ifmtx,*) 'ifcoalinflcell:',ifcoalinflcell
      else
         write(ifmtx,*) 'ifcoalinflcell:','UNDEFINED'
         close(ifmtx)
         stop
      endif
      if(ifinittempcoal) then
         write(ifmtx,*) 'tempcoal:',tempcoal
      else
         write(ifmtx,*) 'tempcoal:','UNDEFINED'
         close(ifmtx)
         stop
      endif
      if(ifinitboltz_model) then
         write(ifmtx,*) 'boltz_model:',boltz_model
      else
         write(ifmtx,*) 'boltz_model:','UNDEFINED'
         close(ifmtx)
         stop
      endif
      if(ifinittypelpmdamp) then
         write(ifmtx,*) 'typelpmdamp:',typelpmdamp
      else
         write(ifmtx,*) 'typelpmdamp:','UNDEFINED'
         close(ifmtx)
         stop
      endif
      if(ifinitfpchoice) then
         write(ifmtx,*) 'fpchoice:',fpchoice
      else
         write(ifmtx,*) 'fpchoice:','UNDEFINED'
         close(ifmtx)
         stop
      endif
      if(ifinittunchoice) then
         write(ifmtx,*) 'tunchoice:',tunchoice
      else
         write(ifmtx,*) 'tunchoice:','UNDEFINED'
         close(ifmtx)
         stop
      endif
      if(ifinitifdecrease) then
         write(ifmtx,*) 'ifdecrease:',ifdecrease
      else
         write(ifmtx,*) 'ifdecrease:','UNDEFINED'
         close(ifmtx)
         stop
      endif
      if(ifinitifcol) then
         write(ifmtx,*) 'ifcol:',ifcol
      else
         write(ifmtx,*) 'ifcol:','UNDEFINED'
         close(ifmtx)
         stop
      endif
      if(ifinitifrad) then
         write(ifmtx,*) 'ifrad:',ifrad
      else
         write(ifmtx,*) 'ifrad:','UNDEFINED'
         close(ifmtx)
         stop
      endif
      if(ifinitiffullQCD) then
         write(ifmtx,*) 'iffullQCD:',iffullQCD
      else
         write(ifmtx,*) 'iffullQCD:','UNDEFINED'
         close(ifmtx)
         stop
      endif
      if(ifinitiflpm) then
         write(ifmtx,*) 'iflpm:',iflpm
      else
         write(ifmtx,*) 'iflpm: UNDEFINED'
         close(ifmtx)
         stop
      endif
      if(ifinitifdamp) then
         write(ifmtx,*) 'ifdamp:',ifdamp
      else
         write(ifmtx,*) 'ifdamp: UNDEFINED'
         close(ifmtx)
         stop
      endif
      if(ifinitcoll) then
         write(ifmtx,*) 'jahqcoll:',jahqcoll
      else
         write(ifmtx,*) 'jahqcoll:','UNDEFINED'
         close(ifmtx)
         stop
      endif
      if(ifinitish) then
         write(ifmtx,*) 'jaish:',jaish
      else
         write(ifmtx,*) 'jaish:','UNDEFINED'
         close(ifmtx)
         stop
      endif
      if(ifinitifptcles) then
         write(ifmtx,*) 'ifptcles:',ifptcles
      else
         write(ifmtx,*) 'ifptcles:','UNDEFINED'
         close(ifmtx)
         stop
      endif
      write(ifmtx,*) '************************************************'
      end


      subroutine exiteposevent(iexit,n,ntotal)
      use JA
      use PBGinixypos
      use PBGiniphasespace
      use PBGnumberevent
      use PBGsystem
      use PBGhqprod
      use PBGevcount
      use PBGmainpro1
      use PBGtransfer
      use PBGnumbercoll
      use PBGtrigger
      use PBGifskip
      integer iexit
#include "aaa.h" 
#include "xan.h"
      data eventsyy/0/
      integer sumHQ,sumbquarks
      double precision ran2
      double precision:: rtotal,rbbbar,rccbar
      integer ifmtx

      call getMonitorFileIndex(ifmtx)
CJA
      jancoll=0
CJA
      eventsyy=eventsyy+1
      eventsxx=eventsyy
 

      sumHQ=0
      sumbquarks=0
      call getevt(nev,phi,phir,bim,egy,npt,ngl,kol)
      sqrts=egy
      ibev=11
      if(bim.gt.xpara(1,1))ibev=1
      if(bim.gt.xpara(2,1))ibev=2
      if(bim.gt.xpara(3,1))ibev=3
      if(bim.gt.xpara(4,1))ibev=4
      if(bim.gt.xpara(5,1))ibev=5
      if(bim.gt.xpara(6,1))ibev=6
      if(bim.gt.xpara(7,1))ibev=7
      if(bim.gt.xpara(8,1))ibev=8
      if(bim.gt.xpara(9,1))ibev=9
      if(bim.gt.xpara(10,1))ibev=10
      if(jaish.ge.6)write(108,*)(xpara(i,1),i=1,10),bim,ibev 

      jaeve(ibev)=jaeve(ibev)+1
      jaglcol(ibev)=jaglcol(ibev)+nglevt
      weight=weight+ngl
      if(iskip.eq.0) return

      if(ifhqfromhq) then
         numberofNNColls=kol 
         numberelemcolls=0
         do ko=1,kol
            do ni=1,ifnpom(ko)
               numberelemcolls=numberelemcolls+1
            enddo
         enddo
         if(abs(sqrtS-200)<0.1) then
            rtotal=0.013974d0
            rccbar=0.013904d0
            rbbbar=0.00007d0
         elseif(abs(sqrtS-5250)<260d0) then
            rtotal=0.097435d0
            rccbar=0.092901d0
            rbbbar=0.004534d0
         elseif(abs(sqrtS-2755)<10.) then
            rtotal=0.0846d0
            rccbar=0.0819d0
            rbbbar=0.0027d0
         else
            write(ifmtx,*) 'unforeseen sqrtS:',sqrtS
            close(ifmtx)
            stop
         endif
         numberprodpoints=numberofNNColls
         do k=1,numberprodpoints
            if(ran2()<=rtotal) then 
               sumHQ=sumHQ+1
            endif
         enddo
         numberofHQ=2*sumHQ
         write(ifmtx,*) '============================================'
         write(ifmtx,*) numberofHQ,' HQ generated out of ',
     &        numberofNNColls,'NN collisions'
         write(ifmtx,*) 'Glauber number:',ngl
         write(ifmtx,*) '============================================'
      else
         do i=1,npt
            call getptl(imo,ifa,i,ic1,ic2,id,ist,ity
     &              ,p1,p2,p3,p4,p5,x1,x2,x3,x4)
            if(((id.eq.4.or.id.eq.-4).and.ist.eq.21).or.
     &          ((id.eq.5.or.id.eq.-5).and.ist.eq.21)) then     
               if(ifhighptselect) then 
                  if(sqrt(p1**2+p2**2).gt.highptcut) then
                     sumHQ=sumHQ+1
                  endif
               else
                  sumHQ=sumHQ+1
               endif
            endif
            if((id.eq.5.or.id.eq.-5).and.ist.eq.21) then     
               sumbquarks=sumbquarks+1
            endif
         enddo
         if(sumHQ.ne.0)
     &    write(ifmtx,*) sumHQ,' heavy quarks produced by init EPOS'
         if(ifbselect) write(ifmtx,*) sumbquarks,
     &        ' bquarks produced by init EPOS'
      endif
      if((sumHQ.eq.0).or.(ifbselect.and.(sumbquarks.eq.0))) then
        iexit=1
        totcol=totcol+nbcol
        if(n.eq.ntotal) call close_pack(totcol,weight)
      else
         iexit=0
      endif

      end


      subroutine aacharm

      call checkeposhq(1)
      call eposhq
      call hq7
      call dmesoneposhq
      call checkeposhq(2)
 
      end

      subroutine eposhq
      use PBGiniphasespace
      use PBGinixypos
      use PBGsystem
      use PBGsystem2
      use PBGhqprod
      use PBGproductionmech
      use PBGhydro
      use PBGnumbercoll
      use PBGcharge
      use PBGifhq2
      implicit none
#include "aaa.h"
      integer i,ic1,ic2,id,ifa,imo,ishini,ist,ity,ko,m,nev,ni,
     & ntauhyxx,nxhyxx,nyhyxx,nzhyxx,ifnpom
      real egy,p1,p2,p3,p4,p5,t,taumaxhyx,tauminhyx,x,y,z,x1,x2,x3,x4,
     & xmaxhyx,xminhyx,ymaxhyx,yminhyx,zmaxhyx,zminhyx
      integer sumHQ,Zproj,Ztarg,Aproj,Atarg
      integer imother,eposmother
      real cphi,sphi
      integer ifmtx

      call getMonitorFileIndex(ifmtx)

      call utpri('eposhq',ish,ishini,4)

      if(ihq2.ne.1) then
        ncharged=-1
        return
      endif
      
      write(ifmtx,*) 'start eposhq'
      call getevt(nev,phi,phir,bim,egy,npt,ngl,kol) 
      sqrtS=egy
      b=bim
      cphi=cos(phi+phir)
      sphi=sin(phi+phir)
      call getnuclei(Zproj,Aproj,Ztarg,Atarg)
      if(ifhqfromhq.and.(Aproj.ne.Atarg)) then
         write(ifmtx,*) 'not a symmetric collision'
         close(ifmtx) 
         stop
      else
         A=(Aproj+Atarg)/2
         write(ifmtx,*) 'A for MC@HQ=',A
      endif
      call gethypar(nzhyxx,ntauhyxx,nxhyxx,nyhyxx,zminhyx
     &                ,zmaxhyx,tauminhyx,taumaxhyx
     &                ,xminhyx,xmaxhyx,yminhyx,ymaxhyx)

      ntauhyx=ntauhyxx
      nxhyx=nxhyxx
      nyhyx=nyhyxx 
      nzhyx=nzhyxx

      xmin=xminhyx
      xmax=xmaxhyx
      xstep=(xmax-xmin)/(nxhyxx-1)

      ymin=yminhyx
      ymax=ymaxhyx
      ystep=(ymax-ymin)/(nyhyxx-1)

      zmin=zminhyx
      zmax=zmaxhyx
      zstep=(zmax-zmin)/(nzhyxx-1)

      tbmin=tauminhyx
      tbmax=taumaxhyx
      tbstep=(tbmax-tbmin)/(ntauhyx-1)

      if(ifhqfromhq) then
         m=0
         do ko=1,kol
           if(hhbarprod.eq.1) then
              ni=ifnpom(ko)
              if(ni.ge.1) call getpos(ko,ni,x,y,z,t)
              xyarray(ko,0)=t
              xyarray(ko,3)=z
              xyarray(ko,1)=x
              xyarray(ko,2)=y
           elseif(hhbarprod.eq.2) then
              do ni=1,ifnpom(ko)
                 m=m+1
                 call getpos(ko,ni,x,y,z,t)
                 xyarray(m,0)=t
                 xyarray(m,3)=z
                 xyarray(m,1)=x
                 xyarray(m,2)=y
              enddo
           else
              write(ifmtx,*) 'hhbarprod undef:',hhbarprod
              close(ifmtx)
              stop
           endif
        enddo
      else
         sumHQ=0
         ncharged=0
         imother=0
         eposmother=0
         do i=1,npt  
            call getptl(imo,ifa,i,ic1,ic2,id,ist,ity
     &           ,p1,p2,p3,p4,p5,x1,x2,x3,x4)
            if(abs(id).eq.4.or.abs(id).eq.5.)then 
               if((ist.gt.20.5).and.(ist.lt.21.5)) then
                  sumHQ=sumHQ+1
                  HQarray(sumHQ,1)=dble(id)
                  if(imo.gt.eposmother) then
                     eposmother=imo
                     imother=imother+1
                  endif
                  hQmother(sumHQ)=imother
                  HQarray(sumHQ,2)=cphi*dble(p1)+sphi*dble(p2)
                  HQarray(sumHQ,3)=-sphi*dble(p1)+cphi*dble(p2)
                  HQarray(sumHQ,4)=dble(p3)
                  HQarray(sumHQ,5)=dble(p4)
                  HQarray(sumHQ,6)=dble(p5)
                  HQarray(sumHQ,7)=cphi*dble(x1)+sphi*dble(x2)
                  HQarray(sumHQ,8)=-sphi*dble(x1)+cphi*dble(x2)
                  HQarray(sumHQ,9)=dble(x3)
                  HQarray(sumHQ,10)=dble(x4)
               endif  
            endif            
         enddo
         numberofHQ=sumHQ
         write(ifmtx,*)'number of c / cbar / b / bbar = ', sumHQ
         write(*,*)'number of c / cbar / b / bbar = ', sumHQ
         call utprix('eposhq',ish,ishini,4)
      endif
      write(ifmtx,*)'leave eposhq'      
      end

      subroutine dmesoneposhq
      use JA
      use PBGhqdmesons
      use PBGhquarks
      use PBGgenvar
      use PBGDchoice
#include "aaa.h"
      integer i,ires
      real r,x1,x2,x3,x4
      integer idKlaus
      integer nptlbefore
      integer ifmtx
      call getMonitorFileIndex(ifmtx)
      call utpri('dmeson',ish,ishini,4)
      if(jahqcoll.eq.1) open(unit=89,file=fnus2(1:nfnus2)//'4',
     &     status='unknown')


      if(ish.ge.4) write(ifch,*)'idmeson=',idmeson,'   nd=',nd
      if(idmeson.eq.1) then
         do i=1,nptl
            call getidptl(i,id)
            call getityptl(i,ity)
            if(ity.eq.61)id=idtrafo('pdg','nxs',id)
            inodaug=0
            call idflav(id,if1,if2,if3,idu,0)
            if(abs(if1).eq.5)inodaug=1
            if(abs(if2).eq.5)inodaug=1
            if(abs(if3).eq.5)inodaug=1
            if(abs(if1).eq.4)inodaug=1
            if(abs(if2).eq.4)inodaug=1
            if(abs(if3).eq.4)inodaug=1
            call getistptl(i,istx) 
            if(inodaug.eq.1.and.(istx.eq.0.or.istx.eq.1))then 
               call setistptl(i,2)  
               call getistptl(i,ist) 
               call getidptl(i,id1)
               call getistptl(i,istl1)
               call getityptl(i,ity1)
               if(ish.ge.6)
     .              write(ifch,'(a,2i9,3x,3i4,3x,2i4,3x,3i3,4x,2i3)')
     .              'HF-catch',id1,id,istl1,ist,istx
     .              ,ity1,ity,if1,if2,if3,inodaug,ittt
            endif
         enddo
         nptlbefore=nptl
         if(nHQ.ne.nd) then
            write(ifmtx,*) 'n HQ <>n HQ mesons in dmesoneposhq'
            write(ifmtx,*) nHQ, nd
            close(ifmtx)
            stop
         endif
         do i=1,nHQ
            ires=idhquarks(i)
            nptl=nptl+1
            call checkcccptl(nptl) 
            id=idKlaus(ires)
            call setidptl(nptl,id)
            p1=phquarks(i,1)
            p2=phquarks(i,2)
            p3=phquarks(i,3)
            p4=phquarks(i,0)
            if(abs(id).eq.4) p5=1.5E0
            if(abs(id).eq.5) p5=5.1E0
            call setpptl(nptl,p1,p2,p3,p4,p5)
            x1=hqmeson(i,7)
            x2=hqmeson(i,8)
            x3=hqmeson(i,9)
            x4=hqmeson(i,6)
            call setxorptl(nptl,x1,x2,x3,x4)
            call setityptl(nptl,60)
            call setiorptl(nptl,0)
            call setjorptl(nptl,0)
            call setistptl(nptl,101)
            call getidptl(nptl,id1)
            call getistptl(nptl,istl1)
            call getityptl(nptl,ity1)
            if(ish.ge.6)Write(101,*)'HQ out',i,ires,nHQ,nptl
         enddo
         write(ifmtx,*) '# of HF mesons',nd

         do i=1,nd
            ires=nint(hqmeson(i,1))
            nptl=nptl+1
            call checkcccptl(nptl) 
            id=idKlaus(ires)
            call setidptl(nptl,id)
            p1=hqmeson(i,3)
            p2=hqmeson(i,4)
            p3=hqmeson(i,5)
            p4=hqmeson(i,2)
            p5=resinfo_mass(ires)
            call setpptl(nptl,p1,p2,p3,p4,p5)
            x1=hqmeson(i,7)
            x2=hqmeson(i,8)
            x3=hqmeson(i,9)
            x4=hqmeson(i,6)
            call setxorptl(nptl,x1,x2,x3,x4)
            call setityptl(nptl,60)
            call setiorptl(nptl,nptl-nd)
            call setifrptl(nptl-nd,nptl,0)
            call setjorptl(nptl,0)
            call setistptl(nptl,0)
            call setifrptl(nptl,0,0)
            call getxorptl(nptl,x1,x2,x3,x4)
            t1=x4
            r=rangen()
            call idtau(abs(id),p4,p5,taugm)
            t2=t1+taugm*(-alog(r))
            call settivptl(nptl,t1,t2)
            call getidptl(nptl,id1)
            call getistptl(nptl,istl1)
            call getityptl(nptl,ity1)
            if(ish.ge.6)Write(101,*)'HQ out',i,ires,nd,nptl
            if(jahqcoll.eq.1) write(89,'(a,2i3,i8,i6,2i4,3x,9e10.3)')
     &           'HQ',i,ires,nptl,id1,ist1,ity1
     &           ,hqmeson(i,2),hqmeson(i,3),hqmeson(i,4), hqmeson(i,5)
     &           ,hqmeson(i,6),hqmeson(i,7),hqmeson(i,8),hqmeson(i,8)
     &           ,t2
         enddo
         nptlbd=nptl
      endif
      close(97)
      deallocate(hqmeson)
      if(jahqcoll.eq.1) close(89)
      call utprix('dmeson',ish,ishini,4)
      end

      subroutine checkeposhq(iii)
#include "aaa.h"
      write(ifmt,'(a,i1,a,i8)')'checkeposhq ',iii,'    nptl = ',nptl
      if(icccptl.eq.0.and.nptl.gt.mxptl)then
        write(ifmt,*)
     .  'checkeposhq - iii icccptl nptl mxptl : ',iii,icccptl,nptl,mxptl
        stop'####### ERROR 08032017 in checkeposhq #######'
      endif
      end

      subroutine idflav2(id,ifail,ifl1,ifl2,ifl3,jspin,idu)
#include "aaa.h"
      integer ifail
      
      ifail=0
      ida=iabs(id)
      idu0=idu
      if(ida.le.9900)then
        nl=nlidtbl(abs(id))
        if(nl.eq.0)then
           ifail=1
           return
        else
           goto 999
        endif 
      endif  

      do nlx=1,nidtmxx(1)
        if(idtbl(1,nlx).eq.id)then
          nl=nlx
          ida=0
          goto 999
        endif
      enddo

      ifl1=0
      ifl2=0
      ifl3=0
      jspin=0
      return

 999  continue

      ifl1=ifl1tbl(nl)
      ifl2=ifl2tbl(nl)
      ifl3=ifl3tbl(nl)
      jspin=jspintbl(nl)
 
      if(id.eq.-ida)then
        ifl1=-ifl1tbl(nl)
        ifl2=-ifl2tbl(nl)
        ifl3=-ifl3tbl(nl)
        jspin=jspintbl(nl)
      endif       
      
      end

      integer function idKlaus(idPB)
      implicit none
      character*2 tochain
      integer idPB
      double precision ran2
      if(idPB.eq.1) then
         idKlaus=4
      elseif(idPB.eq.2) then
         idKlaus=-4
      elseif(idPB.eq.14) then
         idKlaus=5
      elseif(idPB.eq.15) then
         idKlaus=-5
      elseif(idPB.eq.4) then
         if(ran2().lt.0.6658d0) then        
         idKlaus=-140
         else
         idKlaus=-240
         endif
      elseif(idPB.eq.5) then
         if(ran2().lt.0.6658d0) then
         idKlaus=140
         else
         idKlaus=240
         endif
      elseif(idPB.eq.6) then
         idKlaus=-340
      elseif(idPB.eq.7) then
         idKlaus=340
      elseif(idPB.eq.8) then
         idKlaus=2140
      elseif(idPB.eq.9) then
         idKlaus=-2140
      elseif(idPB.eq.10) then
         if(ran2().lt.0.5d0) then
         idKlaus=3240
         else
         idKlaus=3140
         endif
      elseif(idPB.eq.11) then
         if(ran2().lt.0.5d0) then
         idKlaus=-3240
         else
         idKlaus=-3140
         endif
      elseif(idPB.eq.12) then
         idKlaus=3340
      elseif(idPB.eq.13) then
         idKlaus=-3340
      elseif(idPB.eq.17) then
         if(ran2().lt.0.5d0) then
         idKlaus=-250
         else
         idKlaus=-150
         endif
      elseif(idPB.eq.18) then
         if(ran2().lt.0.5d0) then
         idKlaus=250
         else
         idKlaus=150
         endif       
      elseif(idPB.eq.19) then
         idKlaus=-350
      elseif(idPB.eq.20) then
         idKlaus=350
      elseif(idPB.eq.21) then
         idKlaus=2150
      elseif(idPB.eq.22) then
         idKlaus=-2150
      elseif(idPB.eq.23) then
         if(ran2().lt.0.5d0) then
         idKlaus=3250
         else
         idKlaus=1350
         endif
      elseif(idPB.eq.24) then
         if(ran2().lt.0.5d0) then
         idKlaus=-3250
         else
         idKlaus=-1350
         endif
      elseif(idPB.eq.25) then
         idKlaus=3350
      elseif(idPB.eq.26) then
         idKlaus=-3350
      else
         idKlaus=0
         call utstop('unrecognized idPB in idKlaus:'//
     &   tochain(idPB)//'&')
      endif
      end


