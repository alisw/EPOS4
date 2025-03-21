!
!  This file is part of EPOS4
!  Copyright (C) 2022 research institutions and authors (See CREDITS file)
!  This file is distributed under the terms of the GNU General Public License version 3 or later
!  (See COPYING file for the text of the licence)
!

!#######################################################################
!   key subroutines, collected here in order not to hide them in the code
!     containing crucial definitions and parameter settings
!           to be interpolated eventually
!########################################################################

!##################################################################################################################
      subroutine getMaxValues(ZpartMax,ZpomMax) 
!##################################################################################################################
      character*500 path
      character *80 line
      data ncount/0/
      save ncount,valZ,valZp
      if(ifappl().ge.3.or.ifappl().le.0)return
      ncount=ncount+1
      path=' ' 
      if(ncount.eq.1)then   
        !call getReaction(laproj,maproj,latarg,matarg,engy) !engy=egyevt not yet defined
        call getSystemABE(maproj,matarg,engy) 
        call getEposPath(path,length)
        open(105,file= path(1:length)//'/src/KWt/maxValues.dt',status='old')
        rewind(105)
        valE=-1.0
        do 
          read(105,'(a)',end=777)line
          ine=0
          do i=1,80
            if(line(i:i).ne.' ')ine=1
          enddo
          if(line(1:1).ne.'!'.and.ine.ne.0)then
            read(line,*,end=777)maproj_,matarg_,valE_,valN_,sloN_,valNp_,sloNp_
            if(((maproj.eq.maproj_.and.matarg.eq.matarg_).or.(maproj.eq.matarg_.and.matarg.eq.maproj_)))then
              valE  = valE_
              valN  = valN_
              sloN  = sloN_
              valNp = valNp_
              sloNp = sloNp_
            endif
            !write(*,'(a,5f10.3)')'getMaxValues',valE,valN,sloN,valNp,sloNp
          endif
        enddo
 777    close(105)
        if(valE.lt.0.0)stop 'ERROR Entry missing in maxValues.dt'
        valZ = valN  + sloN  *(log(engy)-log(valE))
        valZp= valNp + sloNp *(log(engy)-log(valE))
      endif
      ZpartMax=valZ
      ZpomMax=valZp
      !write(*,'(a,2f10.3)')'getMaxValues',ZpartMax,ZpomMax
      end

!##################################################################################################################
      subroutine setParamsIco() !must be called in IcoStr
!##################################################################################################################
      call  getSystemType(isys,amassMax,amassAsy)
      call getSystemABE(maproj,matarg,engy) 
      call eventvariget(10,x1)  ! = z10zevt = Conn1
      call eventvariget(11,x2)  ! = z11zevt = Conn2
      if(engy.lt.0.1)stop 'ERROR `engy` seems not to be defined' 
      X=(x1+x2)/2
      !fX1=0.40-max(0.,9-X)*0.03
      !fX2=0.40!-max(0.,9-X)*0.03
      !-------------------------------------------------------------------------------------------------------------------core1--ficoscale
      !fico1 + min( fico2*Z**fico4 , fico5-fico6*(Z-1) )  Z=ngl/fico7  ( fico1/ mult\ )
      if(isys.eq.1)then 
        if(engy.lt.400.)then       ; call core1paramset7( 0.10 , 0.00, 0., 0.33, 1.00, 0.00,  30.) !pp200
        elseif(engy.lt.4000.)then  ; call core1paramset7( 0.00 , 0.00, 0., 0.33, 1.00, 0.00, 400.) !pp546..3T
        elseif(engy.lt.10000.)then ; call core1paramset7(-0.12 , 0.00, 0., 0.33, 1.00, 0.00, 400.) !pp5T..8T
        else                       ; call core1paramset7(-0.12 , 0.00, 0., 0.33, 1.00, 0.00, 400.) !pp10T+
        endif 
      elseif(isys.eq.2)then
        if(engy.lt.400.)then ; call core1paramset7(0.40,-0.25, 0., 0.33, 1.00, 0.00,  30.) !dAu200    
        else                 ; call core1paramset7(0.40,-0.20, 0., 0.33, 1.00, 0.00,  30.) !pPb5T
        endif 
      elseif(isys.eq.3)then
        if(engy.lt.400.)then ; call core1paramset7(0.00, 0.00, 0., 0.33, 0.99, 0.00, 400.) !RHIC ***UNUSED***    
        else                 ; call core1paramset7(0.00, 0.00, 0., 0.50, 0.99, 0.00, 250.) !LHC ***UNUSED***
        endif 
      elseif(isys.eq.4)then
        if(engy.lt.400.)then ; call core1paramset7(0.00, 0.00, 0., 0.33, 0.99, 0.00, 400.) !RHIC ***UNUSED***    
        else                 ; call core1paramset7(0.00, 0.00, 0., 0.50, 0.99, 0.00, 400.) !LHC ***UNUSED***
        endif 
      endif 
      end 
!##################################################################################################################
      subroutine setParamsFluct() !must be called in conaa
!##################################################################################################################
      call  getSystemType(isys,amassMax,amassAsy)
      call getSystemABE(maproj,matarg,engy) 
      call getRng1(Rng1) !participant ratio estimate (0->1)
      call getDisize(disize_curr)
!########################################################################################################
      !------------------------------------------------------------------------------------screen1--eflu,efly,efli
      !pp: factor ( 1 - efli * exp(-(b/eflu)**efly) )
      if(engy.lt.400.)then !.........pp200
       efluPP=0.50!0.35
       eflyPP=1.5!3.0
       efliPP=-0.5!improves even P(b)cut
      elseif(engy.lt.2000.)then!.....pp546,pp900,pp1800 
       efluPP=0.50!0.35
       eflyPP=1.5!3.0
       efliPP=0.6*(rangen()-0.45) 
      elseif(engy.lt.10000.)then!.....3T,5T,7T,8T
       efluPP=0.50!0.35
       eflyPP=1.5!3.0
       efliPP=0.6*(rangen()-0.45)
      else !.........................13T
       efluPP=0.50!0.35
       eflyPP=1.5!3.0
       efliPP=0.43*(rangen()-0.45)
      endif!.............................
      ransize=rangen()-0.45 !<-----------------may affects also disize
      efluPA=0!5.0e-7*(rangen()**(-1.5)-1) 
      eflyPA=0.8*ransize !check also znurho
      efliPA=0  
      call screen1paramset6(efluPP,eflyPP,efliPP,efluPA,eflyPA,efliPA)
      !positive eflyPA -> reduction of screening -> more mult
!########################################################################################################
      !------------------------------------------------------------------------------------------dipole
      if(amassAsy.le.0.1)then!------------------------pp,AA--------------------------------
        if(amassMax.le.10.)then!pp
          if(engy.lt.10000)then         ; call setDisize(1.6) !pp7T-
          else                          ; call setDisize(1.6) !pp13T
          endif
        elseif(amassMax.le.300.)then!AA
          if(engy.lt.1e30)then          ; call setDisize(1.6)
          endif
        endif
      elseif(amassAsy.le.0.999)then!--------------------dAu-------------------------------
        if(amassMax.le.300.)then
          if(engy.lt.1e30)then          ; call setDisize(2.)
          endif
        endif
      elseif(amassAsy.le.1.0001)then!--------------------pA-------------------------------
        if(amassMax.le.300.)then
          if(engy.lt.1e30)then          ; call setDisize(2.+2*ransize) !(1+2*Rng1)
          endif
        endif
      else;stop 'ERROR 240403'
      endif!-------------------------------------------------------
!###################################################################################################################
      !-----------------------------------------------------------------------------!
      ! make sure to initialize all variables with POSITIVE (NONZERO) values !!!    !
      !-----------------------------------------------------------------------------!
      call getMonitorFileIndex(ifmt)
      if(disize_curr.le.0.)then
        call setDisize(abs(disize_curr))
        write(ifmt,'(a,f5.2)')'disize forced to',abs(disize_curr)
      endif
      end
!##################################################################################################################
      subroutine setParams() 
!##################################################################################################################
      real ab(6),ac(3)
      call  getSystemType(isys,amassMax,amassAsy)
      call getSystemABE(maproj,matarg,engy) 
      !-----------------------------------------------------------------------------------------screen2--nuclear
      call screen2paramzero() !f1*Z**f2 + f3
      call setZnurho(0.0)
      if(isys.le.2)then!pp aA
        if(engy.lt.400.)then ; call screen2paramset5( 0.07, 0.10, 0.00, 0.00, 30.  ) !dAu200
        else                 ; call screen2paramset5( 0.22, 0.30, 0.00, 0.00, 30.  ) !pPb5T
        endif  
      elseif(isys.eq.3)then
                                call setZnurho(0.0385) !mid AA
      elseif(isys.eq.4)then
        if(engy.lt.30.)then       ; call setZnurho(0.20)  !AuAu27     
        elseif(engy.lt.40.)then   ; call setZnurho(0.20)  !AuAu39  
        elseif(engy.lt.70.)then   ; call setZnurho(0.30)  !AuAu62    
        elseif(engy.lt.150.)then  ; call setZnurho(0.215) !AuAu130  
        elseif(engy.lt.400.)then  ; call setZnurho(0.180) !AuAu200
        elseif(engy.lt.4000.)then ; call setZnurho(0.035) !PbPb3T
                                    call screen2paramset5( 0.04, 6.00, 0., 0., 400. ) !finetune
        else                      ; call setZnurho(0.032) !PbPb5T  
                                    call screen2paramset5( 0.04, 6.00, 0., 0., 400. ) !finetune
        endif 
      endif
      call getIphsd(iphsd)
      if(iphsd.ne.0)then
        call setZnurho(0.10)
      endif
      !-------------------------------------------------------------screen3--Xsection
      if(isys.gt.0)then                        
        if(engy.lt.400.)then!..........200G
          call setEpscrXiG(-0.160, 0.25, 1)  
        elseif(engy.lt.4000.)then!......546...3T
          call setEpscrXiG(-0.160, 0.30, 1)  
        elseif(engy.lt.10000.)then!......5T,7T
          call setEpscrXiG(-0.160, 0.25, 1) 
        else !..........................13T
          call setEpscrXiG(-0.160, 0.22, 1)  
        endif                     
        !---------------------------------------
        !constant xs (from some earlier version) : 
        !factk     1.80    1.80    1.80    2.40  / Xsect /
        !epscrxi  -0.147  -0.165  -0.180  -0.200 \ Xsect \
        !epscrg    0.45    0.25    0.12    0.25  \ Xsect /
        !mult      10.0    10.5    11.0
        !---------------------------------------
      endif
      !--------------------------------------------------------------leadings in core
      !leadcore = abcdwxyz1 abcd wgt perph, wxyz wgt centr
      if(engy.lt.10.)then!.......................AuAu7
        call setLeadcore(1)
      elseif(engy.lt.12.)then!...................AuAu11
        call setLeadcore(050005001) 
      elseif(engy.lt.15.)then!...................AuAu14
        call setLeadcore(050005001)
      else
        call setLeadcore(1) !no leadings in core
      endif
      !--------------------------------------------------------------remnant excitation
      !remnMassMin   =  f1+Z*(f2-f1)     Z=ngl/f3
      !forceDropletMass = f4+Z*(f5-f4) 
      !remnRapScal   =  f6
      alpremn=2.5 
      if(isys.eq.1)then!...........................................pp 
          call remn1paramset6( 10.0, 10.0, 400., 1e10, 1e10, 0.0)         
      elseif(isys.eq.2)then!.......................................aA
          call remn1paramset6( 10.0, 10.0, 400., 1e10, 1e10, 0.0)         
      elseif(isys.eq.3)then!.......................................mid AA
          call remn1paramset6( 10.0, 10.0, 400., 1e10, 1e10, 0.0) 
      elseif(isys.eq.4)then!big AA 
        if(engy.lt.5.)then!........................................AuAu4     
          call remn1paramset6(  1.0,  2.0, 400.,  1.0,  1.8, 0.0) 
        elseif(engy.lt.10.)then!...................................AuAu7     
          call remn1paramset6(  1.7,  2.0, 400.,  2.2,  1.8, 0.0) 
        elseif(engy.lt.12.)then!...................................AuAu11     
          call remn1paramset6(  2.0,  3.8, 400.,  2.5,  3.3, 0.0) 
        elseif(engy.lt.15.)then!...................................AuAu14     
          call remn1paramset6(  2.0,  4.2, 400.,  3.0,  3.8, 0.0) 
        elseif(engy.lt.20.)then!...................................AuAu19     
          call remn1paramset6(  3.0,  4.7, 400.,  3.0,  4.0, 0.0) 
        elseif(engy.lt.30.)then!...................................AuAu27     
          call remn1paramset6(  3.0,  5.0, 400.,  3.0,  4.0, 0.0) 
        elseif(engy.lt.40.)then!...................................AuAu39  
          call remn1paramset6(  3.0,  5.0, 400.,  3.0,  3.0, 0.0) 
        elseif(engy.lt.70.)then!...................................AuAu62    
          call remn1paramset6(  5.5,  5.5, 400.,  3.0,  3.0, 0.0) 
        elseif(engy.lt.150.)then!..................................AuAu130     
          call remn1paramset6(  8.5,  8.5, 400., 10.0, 10.0, 0.0)          !###############################################
        elseif(engy.lt.400.)then!..................................AuAu200 !# both f4 = 10 and 1e10 with adjusted f1,f2   # 
          call remn1paramset6( 8.75, 8.75, 400., 10.0, 10.0, 0.0)          !# give good results for AuAu200, 10 sw better.#
        elseif(engy.lt.4000.)then!.................................PbPb3T  !#     Try also higher energies !!!            #
          call remn1paramset6( 10.0, 10.0, 400., 1e10, 1e10, 0.0)          !###############################################
        else !.....................................................PbPb5T
          call remn1paramset6( 10.0, 10.0, 400., 1e10, 1e10, 0.0) 
        endif      
      endif
      call setAlpdi(alpremn,alpremn) !->alptr(0,1) , alptr(1,1)  !used in lea.f
      !--------------------------------------------------------------SatPom weight  (affects mult)
      if(isys.gt.0)Y=log(engy)
      if(isys.eq.1)then !pp
        if(engy.lt.400.)then       ;   call setFactsat(500.) !pp200
        elseif(engy.lt.10000.)then ;   call setFactsat(100.) !pp546..pp7T
        else                       ;   call setFactsat(500.) !pp13T 
        endif 
      elseif(isys.eq.2)then!aA
        if(engy.lt.400.)then       ;   call setFactsat(160.) !dAu200
        else                       ;   call setFactsat(500.) !pPb5T  !160. makes little change
        endif 
      elseif(isys.eq.3)then !mid AA
        if(engy.lt.1e30)then       ;   call setFactsat(100.) !XeXe5T
        endif
      elseif(isys.eq.4)then !big AA
        if(engy.lt.400.)then       ;   call setFactsat(100.) !AuAu200
        else                       ;   call setFactsat(200.) !PbPb3T,PbPb5T
        endif 
      endif
      call setAlpsat(1.0)
!#############################################################################################
      call felossget6(ab(1),ab(2),ab(3),ab(4),ab(5),ab(6))
      call getMaxValues(ZpartMax,dmy) 
      Zx=ZpartMax
      !                                                 feloss1 feloss2 feloss3 feloss4 feloss5  feloss6
      !                                                 cone pp cone AA bulk pp bulk AA rad cone NpartMax
!####################################################################################################################
      if(amassAsy.le.0.1)then!------------------------pp,AA--------------------------------
        if(amassMax.le.10.)then!pp
          if(engy.lt.400.)then      ; call felossset6(    0.05 ,  0.05 ,  0.60 ,  0.60 ,  0.3  ,  400. ) !pp200
          elseif(engy.lt.1000.)then ; call felossset6(    0.10 ,  0.10 ,  0.60 ,  0.60 ,  0.3  ,  400. ) !pp546
          else                      ; call felossset6(    0.10 ,  0.10 ,  0.40 ,  0.40 ,  0.3  ,  400. ) !pp LHC
          endif
        elseif(amassMax.le.150.)then!midAA
          if(engy.lt.1e30)then      ; call felossset6(    0.30 ,  0.30 ,  0.50 ,  0.70 ,  0.2  ,  250. ) !XeXe5T
          endif
        elseif(amassMax.le.300.)then!bigAA
          if(engy.lt.70.)then       ; call felossset6(    0.05 ,  0.15 ,  0.20 ,  0.50 ,  0.3  ,  400. ) !AuAu7...AuAu62
          elseif(engy.lt.400.)then  ; call felossset6(    0.05 ,  0.15 ,  0.50 ,  0.50 ,  0.3  ,  400. ) !AuAu130,AuAu200
          else                      ; call felossset6(    0.10 ,  0.25 ,  0.50 ,  0.50 ,  0.3  ,  400. ) !PbPb LHC
          endif
        endif
      elseif(amassAsy.le.0.999)then!--------------------dAu-------------------------------
        if(amassMax.le.300.)then
          if(engy.lt.1e30)then      ; call felossset6(    0.05 ,  0.15 ,  0.60 ,  0.60 ,  0.3  , 400. ) !dAu200
          endif
        endif
      elseif(amassAsy.le.1.0001)then!--------------------pA-------------------------------
        if(amassMax.le.300.)then
          if(engy.lt.1e30)then      ; call felossset6(    0.00 ,  0.00 ,  0.25 ,  0.40 ,  0.3  ,  Zx  ) !pPb5T
          endif
        endif
      else;stop 'ERROR 240317'
      endif!-------------------------------------------------------
!####################################################################################################################
      call getHY(ac(1),ac(2),ac(3))
      !                                         tauzer1 tauzer2 efrout        efrout = FO epsilon
!####################################################################################################################
      if(amassAsy.le.0.1)then!------------------------pp,AA--------------------------------
        if(amassMax.le.10.)then!pp
          if(engy.lt.400.)then       ;  call setHY( 0.4 , 1.5 , 0.30 ) !pp200
          else                       ;  call setHY( 0.4 , 0.4 , 0.30 ) !pp LHC
          endif
        elseif(amassMax.le.150.)then!midAA
          if(engy.lt.1e10)then       ;  call setHY( 1.0 , 1.5 , 0.30 ) !XeXe5T
          endif
        elseif(amassMax.le.300.)then!bigAA
          !engy:   7.0 11.5 14.5 19.6
          !R/gam: 1.62 1.09 0.86 0.64
          if(engy.lt.10.)then        ;  call setHY( 1.0 , 2.0 , 0.30 ) !AuAu7
          elseif(engy.lt.12.)then    ;  call setHY( 1.0 , 1.5 , 0.30 ) !AuAu11
          elseif(engy.lt.15.)then    ;  call setHY( 1.0 , 1.5 , 0.30 ) !AuAu14
          elseif(engy.lt.20.)then    ;  call setHY( 1.0 , 1.5 , 0.30 ) !AuAu19
          elseif(engy.lt.30.)then    ;  call setHY( 1.0 , 1.5 , 0.30 ) !AuAu27
          elseif(engy.lt.40.)then    ;  call setHY( 1.0 , 1.5 , 0.30 ) !AuAu39
          elseif(engy.lt.70.)then    ;  call setHY( 1.0 , 1.0 , 0.30 ) !AuAu62
          elseif(engy.lt.150.)then   ;  call setHY( 1.0 , 1.0 , 0.30 ) !AuAu130
          elseif(engy.lt.400.)then   ;  call setHY( 1.0 , 1.0 , 0.30 ) !AuAu200
          elseif(engy.lt.4000.)then  ;  call setHY( 0.4 , 1.0 , 0.30 ) !PbPb3T
          else                       ;  call setHY( 0.4 , 1.0 , 0.30 ) !PbPb5T
          endif
        endif
      elseif(amassAsy.le.0.999)then!--------------------dAu-------------------------------
        if(amassMax.le.300.)then
          if(engy.lt.1e30)then       ;  call setHY( 0.4 , 1.5 , 0.30 ) !dAu200
          endif
        endif
      elseif(amassAsy.le.1.0001)then!--------------------pA-------------------------------
        if(amassMax.le.300.)then
          if(engy.lt.1e30)then       ;  call setHY( 0.4 , 0.4 , 0.30 ) !pPb5T
          endif
        endif
      else;stop 'ERROR 240403'
      endif!-------------------------------------------------------
!################################################################################################
      !---------------------------------------!
      ! pbreakg  pbreak zipinc  pmqq  zopinc  !----------------------------fragmentation     
      !---------------------------------------!
      if(isys.eq.0)then 
        call setFRA( 0.04 , -0.4 ,  0.11  , 0.11 , 0.0    )
        !define 'fcharm' here
      elseif(isys.ge.1.and.isys.le.4)then
        call setFRA( 0.04 , -0.4 , -0.062 , 0.11 , 0.0025 ) !pp,aA,AA
        !define 'fcharm' here
      else
        stop 'ERROR 25042022'   
      endif
!###################################################################################################################
      !-----------------------------------------------------------------------------!
      ! make sure to initialize all variables with POSITIVE (NONZERO) values !!!    !
      !-----------------------------------------------------------------------------!
      call getMonitorFileIndex(ifmt)
      if(ab(1).le.0.)then
        call felossset(1,abs(ab(1)))
        write(ifmt,'(a,i1,a,f5.2)')'feloss(',1,') forced to',abs(ab(1))
      endif
      do i=1,3
       if(ac(i).le.0.)then
        call setHY1(i,abs(ac(i)))
        if(i.eq.1)write(ifmt,'(a,f5.2)')'tauzer1 forced to',abs(ac(i))
        if(i.eq.2)write(ifmt,'(a,f5.2)')'tauzer2 forced to',abs(ac(i))
        if(i.eq.3)write(ifmt,'(a,f5.2)')'efrout forced to',abs(ac(i))
       endif
      enddo
      end
!##################################################################################################################
      subroutine setParams2()
!##################################################################################################################
      !-----------------------------------------------------------------Flux tube, parton position spread 
      call  getSystemType(isys,amassMax,amassAsy)
      call getSystemABE(maproj,matarg,engy)
      !Too small radeft2 makes a multi layer  FO surface
      !                                                      radeft1         radeft2        facposf        facposz 
!##########################################################################################################################
       if(amassAsy.le.0.1)then!------------------------pp,AA--------------------------------
        if(amassMax.le.10.)then!pp
          call getNhpom(npom)
          if(engy.lt.10000.)then ; call setRadsize(            0.2        ,   0.30   ,        1.3      ,   1.3  ) !pp 7T
          else                   ; call setRadsize(            0.2        ,   0.30   ,        1.3      ,   1.3  ) !pp 13T
          endif
        elseif(amassMax.le.150.)then!XeXe
          if(engy.lt.1e30)then   ; call setRadsize(            0.3        ,   0.24   ,        0.3      ,   0.3  ) !XeXe
          endif
        elseif(amassMax.le.300.)then!PbPb
          if(engy.lt.400.)then   ; call setRadsize(            0.3        ,   0.2    ,        0.3      ,   0.3  ) !AuAu200
          elseif(engy.lt.4000.)then ; call setRadsize(         0.5        ,   1.0    ,        1.4      ,   1.4  ) !PbPb3T
          else                   ; call setRadsize(            0.4        ,   0.4    ,        1.4      ,   1.4  ) !PbPb5T+
          endif
        endif
      elseif(amassAsy.le.0.999)then!--------------------dAu-------------------------------
        if(amassMax.le.300.)then
          if(engy.lt.1e30)then   ; call setRadsize(            0.3        ,   0.3    ,        1.2      ,   1.2  ) !dAu
          endif
        endif
      elseif(amassAsy.le.1.0001)then!--------------------pA-------------------------------
        if(amassMax.le.300.)then
          if(engy.lt.1e30)then   ; call setRadsize(            0.3        ,   0.3    ,        0.4      ,   0.4  ) !pA
          endif
        endif
      else;stop 'ERROR 240317'
      endif!-------------------------------------------------------
!########################################################################################################
      !--------------------------------------------------------------PFE Parameterized Fluid Expansion
      call getMonitorFileIndex(ifmt)
      call getIorsdf(iorsdf)
      if(iorsdf.eq.3)then!+++
      !------------------------------------------- 
      ! pfe3param:  yrmax, ycoi, ycoj, fecc, taufo
      !------------------------------------------- 
      if(isys.eq.1)then!pp
        call getNhpom(N)
        if(engy.lt.10000.)then!...................................pp7T,pp8T
          call pfe3paramset5(0.7+0.04*min(N,12),1.0, 0.0, min(0.015*N,0.15), 1.)
        else!.....................................................pp13T
          call pfe3paramset5(0.7+0.04*min(N,12),1.0, 0.0, min(0.015*N,0.15), 1.)
        endif 
      elseif(isys.eq.2)then!aA
      !...........................................................dAu200,pPb5T
        write(ifmt,'(2a)')'WARNING PFE settings for this system not ','yet tested. This is a non-official extension of EPOS4.0.0'
          call pfe3paramset5(1.3,1.5, 0.0,  0.15 , 1.)
      elseif(isys.eq.3)then!mid AA
        stop 'ERROR No PFE settings for this system yet'
      elseif(isys.eq.4)then!big AA
        call getRng1(Z)
        if(engy.lt.10.)then       ; call pfe3paramset5(0.90-0.05*Z,0.7  , 0.5, max(0.043,min(0.15,.19-.17*Z)), 1+2*Z)    !AuAu7
        elseif(engy.lt.13.)then   ; call pfe3paramset5(0.80-0.10*Z,1.0  , 0.5, max(0.043,min(0.15,.19-.17*Z)), 1+4*Z)    !AuAu11
        elseif(engy.lt.15.)then   ; call pfe3paramset5(0.80-0.05*Z,1.0  , 0.5, max(0.043,min(0.15,.19-.17*Z)), 1+4*Z)    !AuAu14
        elseif(engy.lt.20.)then   ; call pfe3paramset5(0.80-0.05*Z,1.0  , 0.5, max(0.049,min(0.18,.22-.20*Z)), 1+4*Z)    !AuAu19
        elseif(engy.lt.30.)then   ; call pfe3paramset5(0.80-0.05*Z,1.0  , 0.5, max(0.055,min(0.20,.25-.22*Z)), 1+4*Z)    !AuAu27
        elseif(engy.lt.40.)then   ; call pfe3paramset5(0.80+0.05*Z,1.0  , 0.5, max(0.057,min(0.21,.26-.23*Z)), 1+4*Z)    !AuAu39
        elseif(engy.lt.70.)then   ; call pfe3paramset5(0.90+0.05*Z,0.85 , 0.5, max(0.059,min(0.21,.27-.24*Z)), 1+4*Z)    !AuAu62
        elseif(engy.lt.1000.)then ; call pfe3paramset5(0.90+0.05*Z,0.9  , 0.5, max(0.060,min(0.23,.29-.27*Z)), 1+4*Z)    !AuAu130,AuAu200
        elseif(engy.lt.4000.)then ; call pfe3paramset5(1+0.16*Z,1.0     , 0.0, max(0.05,min(0.23,.4-.4*Z))   , 1+5*Z)    !PbPb3T
        else                      ; call pfe3paramset5(1+0.16*Z,0.8+.4*Z, 0.0, max(0.05,min(0.23,.4-.4*Z))   , 1+5*Z)    !PbPb5T
        endif
      else
        stop 'ERROR No PFE option for this system'
      endif
      call pfe3modif() !to allow external finetuning for given system 
      endif!+++
      end

!##################################################################################################################
      subroutine setOptionsWom() 
!##################################################################################################################
      call  getSystemType(isys,amassMax,amassAsy)
      !-----------------------------------------------------------------Wom--- (X dep affects UE)
      dummy=isys
      call deformoptset(1, 2 ) !2 = new method ; 1 = old
      end

!##################################################################################################################
      subroutine defParamsWom(conn12,fuuu) !def fuuu 
!##################################################################################################################
      fuuu=0
      X=conn12                                
      call  getSystemType(isys,amassMax,amassAsy)    
      !amassMax = max atomic mass number (1.0 - 250.)
      !amassAsy = atomic mass number asymmetry (0.0 - 1.0) AA:0.0, pA:1.0, NFe:0.76 
      !in case of pp,pA: X dep such that Spread=0 for all centrality classes
      !in case of HI: X dep such that Rfact=1 for all centralities
      call getSystemABE(maproj,matarg,engy)               
      if(amassAsy.le.0.1)then!------------------------pp,AA--------------------------------
        if(amassMax.le.10.)then!pp
          if(engy.lt.70.)then        ; fuuu=0.21 !pp68
          elseif(engy.lt.400.)then   ; fuuu=0.165*max(min(X,18.),3.)!pp200                 
          elseif(engy.lt.2000.)then  ; fuuu=0.230*max(min(X,18.),3.)!pp546,pp900,pp1800 
          elseif(engy.lt.6000.)then  ; fuuu=0.175*max(min(X,18.),3.)!pp3T,pp5T 
          elseif(engy.lt.10000.)then ; fuuu=0.170*max(min(X,18.),3.)!pp7T,pp8T
          else                       ; fuuu=0.150*max(min(X,18.),3.)!pp13T
          endif 
        elseif(amassMax.le.20.)then!OO
          if(engy.lt.400)then       ; fuuu=max(0.078+0.00264*X,0.0174*(X+0.5))!OO200
          elseif(engy.lt.10000)then ; fuuu=max(0.50+0.0167*X,0.109*(X+0.475))!OO5T
          else                      ; fuuu=max(0.50+0.0167*X,0.109*(X+0.475))
          endif
        elseif(amassMax.le.50.)then!ArAr
          if(engy.lt.10000)then     ; fuuu=max(0.50+0.0167*X,0.109*(X+0.475))!ArAr5T
          else                      ; fuuu=max(0.50+0.0167*X,0.109*(X+0.475))
          endif
        elseif(amassMax.le.100.)then!KrKr
          if(engy.lt.400)then       ; fuuu=3 !ZrZr200
          elseif(engy.lt.10000)then ; fuuu=max(0.57+0.019*X,0.124*(X+0.475))!KrKr5T
          else                      ; fuuu=max(0.57+0.019*X,0.124*(X+0.475))
          endif
        elseif(amassMax.le.150.)then!XeXe
          if(engy.lt.10000)then     ; fuuu=max(0.57+0.019*X,0.124*(X+0.475))!XeXe5T  
          else                      ; fuuu=max(0.57+0.019*X,0.124*(X+0.475))
          endif 
        elseif(amassMax.le.300.)then!PbPb
          if(engy.lt.70.)then       ; fuuu=max(0.60           ,-0.730+1.17*X  ) !AuAu62
          elseif(engy.lt.150.)then  ; fuuu=max(0.43           ,-0.522+0.83*X  ) !AuAu130
          elseif(engy.lt.250.)then  ; fuuu=max(0.225          ,-0.274+0.44*X  ) !AuAu200
          elseif(engy.lt.400.)then  ; fuuu=max(0.3            ,-0.12+0.38*X   ) !PbPb275,PbPb350
          elseif(engy.lt.600.)then  ; fuuu=max(0.3            ,-0.1+0.37*X    ) !PbPb550
          elseif(engy.lt.1000.)then ; fuuu=max(0.3            , 0.07+0.27*X   ) !PbPb900
          elseif(engy.lt.4000.)then ; fuuu=max(0.345+0.0300*X ,-0.405+0.191*X ) !PbPb3T
          else                      ; fuuu=max(0.523+0.0285*X ,-0.200+0.140*X ) !PbPb5T
          endif  
        endif
      elseif(amassAsy.le.0.9)then!-----------------------XA-------------------------------
        if(amassMax.le.300.)then
          if(    engy.lt.  101.)then ; fuuu=1.07 !NFe100
          elseif(engy.lt.  201.)then ; fuuu=1.00 !NFe200
          elseif(engy.lt.  501.)then ; fuuu=0.91 !NFe500
          elseif(engy.lt. 1001.)then ; fuuu=0.80 !NFe1T
          elseif(engy.lt. 2001.)then ; fuuu=0.80 !NFe2T
          elseif(engy.lt. 5001.)then ; fuuu=0.69 !NFe5T
          elseif(engy.lt.10001.)then ; fuuu=0.80 !NFe10T
          else                       ; fuuu=0.80 !NFe20T
          endif
        endif
      elseif(amassAsy.le.0.99)then!--------------------NePb-------------------------------
        if(amassMax.le.300.)then
          if(engy.lt.400.)then       ; fuuu=1.0 !NePb58
          else                       ; fuuu=1.0
          endif
        endif
      elseif(amassAsy.le.0.999)then!--------------------dAu-------------------------------
        if(amassMax.le.300.)then
          if(engy.lt.400.)then       ; fuuu=0.350*X**0.3 !dAu200
          else                       ; fuuu=1.0
          endif
        endif
      elseif(amassAsy.le.1.0001)then!--------------------pA-------------------------------
        if(amassMax.le.30.)then!pO,pNe
          if(engy.lt.400.)then       ; fuuu=0.115  !pNe68
          elseif(engy.lt.7001.)then  ; fuuu=0.500  !pO7T
          else                       ; fuuu=0.475  !pO10T
          endif
        elseif(amassMax.le.300.)then
          if(engy.lt.400.)then       ; fuuu=0.63 !pAu200
          else                       ; fuuu=0.54  !max(0.12+0.020*X,min(-2.3+X,0.32+0.008*X))!pPb5T,pPb8T
                                       call getNpom1(Npom); if(Npom.le.1)fuuu=fuuu*0.5
          endif      
          if(engy.gt.6500.)fuuu=fuuu*0.888!pPb8T finetune
        endif
      else;stop 'ERROR 231018'
      endif!-------------------------------------------------------    
      end

!##################################################################################################################
      subroutine setParamsPtintr() !intrinsic pt
!##################################################################################################################
      real a(4)
      call getMonitorFileIndex(ifmt)
      call  getSystemType(isys,amassMax,amassAsy)
      call getSystemABE(maproj,matarg,engy)
      call getPtintr(a(1),a(2),a(3),a(4))
                              !                       ptipom  ptipomi  ptipos  ptiposi
      if(amassAsy.le.0.1)then!-------------------------------pp,AA-------------------------
        if(amassMax.le.10.)then!pp
          if(engy.lt.400.)then       ; call setPtintr( 1.40 ,  0.45 ,   0.45 ,  0.40  ) !pp RHIC
          elseif(engy.lt.10000.)then ; call setPtintr( 1.40 ,  0.45 ,   1.40 ,  0.40  ) !pp 7T
          else                       ; call setPtintr( 1.00 ,  0.45 ,   2.40 ,  0.40  ) !pp 13T
          endif
        elseif(amassMax.le.150.)then!XeXe
          if(engy.lt.400.)then       ; call setPtintr( 1.40 ,  0.45 ,   0.45 ,  0.40  ) !XeXe RHIC
          else                       ; call setPtintr( 1.40 ,  0.45 ,   0.45 ,  0.40  ) !XeXe LHC
          endif
        elseif(amassMax.le.300.)then!PbPb
          if(engy.lt.400.)then       ; call setPtintr( 1.40 ,  0.45 ,   0.45 ,  0.40  ) !PbPb RHIC
          else                       ; call setPtintr( 2.10 ,  0.45 ,   0.45 ,  0.40  ) !PbPb LHC
          endif
        endif
      elseif(amassAsy.le.0.999)then!---------------------------dA-------------------------
        if(amassMax.le.300.)then
          if(engy.lt.1e30)then; call setPtintr( 1.40 ,  0.45 ,   0.45 ,  0.40  ) !dAu
          endif
        endif
      elseif(amassAsy.le.1.0001)then!--------------------------pA-------------------------
        if(amassMax.le.300.)then
          if(engy.lt.1e30)then; call setPtintr( 1.40 ,  0.45 ,   0.45 ,  0.40  ) !pPb
          endif
        endif
      else;stop 'ERROR 241009'
      endif!-------------------------------------------------------
!##################################################################################################################
      !-----------------------------------------------------------------------------!
      ! make sure to initialize all variables with POSITIVE (NONZERO) values !!!    !
      !-----------------------------------------------------------------------------!
      do i=1,4
       if(a(i).le.0.)then
        call setPtintr1(i,abs(a(i)))
        write(ifmt,'(a,i1,a,f5.2)')'Ptintr(',i,') forced to',abs(a(i)) !Ptintr = ptipom ptipomi ptipos ptiposi
       endif
      enddo
      end

!##################################################################################################################
      subroutine defParamsSat(kcol,ixsplit,wtsplit,wtsplit2,qqsmall,qqfac,qqsfac,qqmax,rapsata,asymmsat)
!##################################################################################################################
      parameter (mxrapsata=8)
      real rapsata(mxrapsata)
      !---------------------------------------------!
      !rejection constants "frej"  and "fscale1sat" !
      !may need some update when changing parameters!
      !---------------------------------------------!
      call  getSystemType(isys,amassMax,amassAsy)
      call getSystemABE(maproj,matarg,engy)
      call getNhpomK(kcol,npom)
      call getRng1(Z) !participant ratio estimate (0->1)
      Z=min(1.,Z)
!########################################################################################################
      !--------------------------------------------------------------SatPom qq distr (split,slope,rapidity) 
      call  getSystemType(isys,amassMax,amassAsy)
      !amassMax = max atomic mass number (1.0 - 250.)
      !amassAsy = atomic mass number asymmetry (0.0 - 1.0) AA:0.0, pA:1.0, NFe:0.76
      ixsplit=1!Different from 1 only for tests
      wtsplit2=1
!########################################################################################################
      if(amassAsy.le.0.1)then!------------------------pp,AA--------------------------------
        if(amassMax.le.10.)then!pp
          if(engy.lt.400.)then      ; wtsplit=0.0 ; wtsplit2=7.0  ; qqsmall=2.75 ; qqfac=0.60  !pp200
          elseif(engy.lt.4000.)then ; wtsplit=0.0 ; wtsplit2=2.6  ; qqsmall=3.5  ; qqfac=0.40  !pp546..pp3T
          else                      ; wtsplit=0.0 ; wtsplit2=2.6  ; qqsmall=3.5  ; qqfac=0.40  !pp5T,pp7T,pp13T
          endif
          avnp=1+(log(engy)-log(200.))/(log(13000.)-log(200.))
          wtsplit2=wtsplit2*(1+0.07*min(npom-avnp,12.))
        elseif(amassMax.le.150.)then!XeXe
          if(engy.lt.1e30)then      ; wtsplit=0.  ; wtsplit2=4.0  ; qqsmall=3.8  ; qqfac=1.2   !XeXe5T
          endif
        elseif(amassMax.le.300.)then!PbPb
          if(engy.lt.400.)then      ; wtsplit=0.  ; wtsplit2=9.0  ; qqsmall=3.8  ; qqfac=0.70  !AuAu200
          else                      ; wtsplit=0.  ; wtsplit2=8.0  ; qqsmall=3.8  ; qqfac=0.70  !PbPb3T,PbPb5T
          endif
          wtsplit2=wtsplit2*(1-0.25*(1-Z))
        endif
      elseif(amassAsy.le.0.999)then!--------------------dAu-------------------------------
        if(amassMax.le.300.)then
          if(engy.lt.1e30)then      ; wtsplit=1.0 ; wtsplit2=1.0  ; qqsmall=2.75 ; qqfac=0.4   !dAu200
          endif
        endif
      elseif(amassAsy.le.1.0001)then!--------------------pA-------------------------------
        if(amassMax.le.300.)then
          if(engy.lt.1e30)then      ; wtsplit=3.3 ; wtsplit2=4.2  ; qqsmall=3.7  ; qqfac=0.20  !pPb5T
          endif
        endif
      else;stop 'ERROR 240317'
      endif!-------------------------------------------------------
!########################################################################################################
      !rapsata(4-8) not used any more 
!########################################################################################################
      if(amassAsy.le.0.1)then!------------------------pp,AA--------------------------------
        if(amassMax.le.10.)then!pp
          if(engy.lt.8001.)then      ; rapsata = (/ 1.7 , 1.7 , 5.0 , 0. , 0. , 0. , 0. , 0.   /) ; asymmsat=0. !pp8T-
          else                       ; rapsata = (/ 2.1 , 2.1 , 5.0 , 0. , 0. , 0. , 0. , 0.   /) ; asymmsat=0. !pp13T
          endif
        elseif(amassMax.le.150.)then!XeXe
                                       rapsata = (/ 1.7 , 1.7 , 5.0 , 0. , 0. , 0. , 0. , 0.   /) ; asymmsat=0.
        elseif(amassMax.le.300.)then!PbPb
                                       rapsata = (/ 1.7 , 1.7 , 5.0 , 0. , 0. , 0. , 0. , 0.   /) ; asymmsat=0.
        endif
      elseif(amassAsy.le.0.9)then!-----------------------XA-------------------------------
        if(amassMax.le.300.)then
          if(    engy.lt. 1001.)then ; rapsata = (/ 1.5 , 1.5 , 5.0 , 0. , 0. , 0. , 0. , 0.   /) ; asymmsat=0.  !NFe1T-
          elseif(engy.lt. 2001.)then ; rapsata = (/ 1.5 , 1.5 , 5.0 , 0. , 0. , 0. , 0. , 0.   /) ; asymmsat=0.  !NFe2T
          elseif(engy.lt. 5001.)then ; rapsata = (/ 1.7 , 1.7 , 5.0 , 0. , 0. , 0. , 0. , 0.   /) ; asymmsat=0.  !NFe5T
          elseif(engy.lt.10001.)then ; rapsata = (/ 2.1 , 2.1 , 5.0 , 0. , 0. , 0. , 0. , 0.   /) ; asymmsat=0.  !NFe10T
          else                       ; rapsata = (/ 2.1 , 2.1 , 5.0 , 0. , 0. , 0. , 0. , 0.   /) ; asymmsat=0.  !NFe20T
          endif
        endif
      elseif(amassAsy.le.0.999)then!--------------------dAu-------------------------------
        if(amassMax.le.300.)then
                                       rapsata = (/ 1.7 , 1.7 , 5.0 , 0. , 0. , 0. , 0. , 0.   /) ; asymmsat=0.
        endif
      elseif(amassAsy.le.1.0001)then!--------------------pA-------------------------------
        if(amassMax.le.300.)then
                                       rapsata = (/ 0.2 , 0.2 , 5.0 , 0. , 0. , 0. , 0. , 0.   /) ; asymmsat=0.2
        endif
      else;stop 'ERROR 240317'
      endif!-------------------------------------------------------
!########################################################################################################
      if(engy.lt.400.)then!....................RHIC
        qqmax=qcmass**2
        qqsfac=2
      else!....................................LHC
        qqmax=qbmass**2
        qqsfac=2
      endif
      end
!##################################################################################################################
      subroutine correctLowNpom(iskip)
!##################################################################################################################
      iskip=0
      call  getSystemType(isys,amassMax,amassAsy)
      call getSystemABE(maproj,matarg,engy)
      if(isys.eq.1)then!aa
        if(engy.lt.2000.)then!......pp200,pp546,pp1800
          npom=igetNpom()
          if(rangen().gt.0.33.and.npom.gt.0)iskip=1 !goto 3 in aepos
        endif
      !elseif(isys.eq.2)then!aA
      !  if(engy.gt.2000.)then!......pPb5T
      !    npom=igetNpom()
      !    if(rangen().gt.0.3+0.7*npom/5.)iskip=1 !goto 3 in aepos
      !  endif
      endif
      end 
!##################################################################################################################
!      !------------------------------------------------------------------------------------------ template
!      if(amassAsy.le.0.1)then!------------------------pp,AA--------------------------------
!        if(amassMax.le.10.)then!pp
!          if(engy.lt.1e30)then      ; 
!          endif
!        elseif(amassMax.le.300.)then!AA
!          if(engy.lt.1e30)then      ; 
!          endif
!        endif
!      elseif(amassAsy.le.0.999)then!--------------------dAu-------------------------------
!        if(amassMax.le.300.)then
!          if(engy.lt.1e30)then      ; 
!          endif
!        endif
!      elseif(amassAsy.le.1.0001)then!--------------------pA-------------------------------
!        if(amassMax.le.300.)then
!          if(engy.lt.1e30)then      ; 
!          endif
!        endif
!      else;stop 'ERROR 240403'
!      endif!-------------------------------------------------------

