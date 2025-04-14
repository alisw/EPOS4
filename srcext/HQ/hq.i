beginoptns hq
!----------------------------------
! control variables on event selection and IS
!----------------------------------
set ihq 1               !1 - EPOS+HQ, other - pure EPOS
set iskip 1             !1 - skip the event w/out charm, other - no skipping
set hqifhqfromhq 0      ! the origin of HQ (put 0 for origin from epos)
set hqbselect 0         ! 1=set a trigger on the production of b quarks
set hqhighptselect 0    ! 1=set a trigger on high pT HQ
set hqhighptcut 0       ! the value of the pT trigger
!-------------------------------
! global control variables on HQ evolution
!-------------------------------
set hqnbcol 1           !total number of HQ runs per one hydro event
set hqnbcolsaveEPOS 1   !total number of HQ runs per one hydro event saved back to EPOS
set hqifallowhqskip 1
!--------------------------------
! control variables for transport and Eloss
!--------------------------------
set hqifcol 1           !if collisionnal process, 1-.true., other - .false.
set hqifrad 1           !if radiative process, 1-.true., other - .false.
set hqiflpm 1           !if LPM implemented, 1-.true., other - .false.
set hqratemult 0.9        !K factor for rates  
set hqfpmult 0.9        !K factor for FP coeff
!---------------------------------
! control variables for hadronization
!---------------------------------
set hqitypdproduct 3    !type of hadronization, 1:fragmentation / 2: coalescence / 3(standard):mixed
set hqitypfrag 3        !fragm. type:1:delta(z=1)/2:Peterson/3(std):Cacciari and Nason
set hqitypcoal 5        !coal type; 1: historic; 2: Dover; 3: 2 with bug corrected (2017); 4: 2017; 5: 2023 inclued all hf hadrons
!---------------------------------
! control variables for feeding back -> EPOS
!---------------------------------
set idmeson 1           !1 - D mesons from HQ, other - D mesons from EPOS 
!---------------------------------
! control variables for URQMD
!---------------------------------
set ivirtual 1          !1 - virtual scattering D with hadrons, other - no 
!---------------------------------
! control variables for log files
!---------------------------------
set hqcoll 0		!if transfer information for pack2 and pack boltzmann_radiat
set hqish 8             ! the level of info in .hq file
set hqifptcles 1        !1 - if final particle list should be written - other; no
!--------------------------------
! other more specific options -- for experts only
!-------------------------------
!set hqifextspectrac 1   ! external spectrum for c (yes:1; no:0)
!set hqtypeextspectrac 1 ! type of external spectrum
!set hqifshadowingc 1    ! if shadowing for c
!set hqtypeshadowc 1     ! shadowing for c (1:central,2:strong shad;3:weak shad)
!set hqifextspectrab 0   ! external spectrum for b (yes:1; no:0)
!set hqtypeextspectrab 1 ! type of external spectrum
!set hqifshadowingb 0    ! if shadowing for b
!set hqtypeshadowb 1     ! shadowing for b (1:central,2:strong shad;3:weak shad)
set hqnbsubrun 1        ! number of subruns (technical), keep it = 1
set hqtype_evol_q 2     ! 1 for Fokker Planck / 2 for Boltzmann evolution
set hqboltz_model 4     ! model for Boltzmann collisions (elast)
set hqfpchoice 4        ! Fokker Planck transp coeff (elast)
set hqtunchoice 1       ! type of tuning for FP coeff
set hqalphasforrad 0.3  ! alpha_s for rad coll 0.3
set hqiffullqcd 1       ! if full QCD in the radiative process, 1-.true., other - .false. 
set hqtypelpmdamp 1     ! type of LPM
set hqiflpm 1           ! if LPM implemented, 1-.true., other - .false.
set hqifdamp 0          ! if damping implemented, 1-.true., other - .false.
set hqc2 0.0            ! damping coeff						set hqif2bodySC 0           ! if 2 body evolution
set hqifsingletred 0
set hqdpmaxabscorona 0      ! 0.05
set hqtypedistrel 2         ! the type of relative distance distribution
set hqparamdistrel(1) 0.25  ! the parameter in the distribution for ccbar relative distance
set hqparamdistrel(2) 0.1   ! the parameter in the distribution for bbbar relative distance
set hqtimestep 0.1      ! (bjorken) time step for evolution
set hqimedfinal 0       ! the end for the evolution, 0 for evol -> QGP] / 1 for evol -> mixed phase]
set hqifdecrease 1      ! if decrease of rate in the mixed phase, 1-.true., other - .false. 
set hqifcronin 0        ! cronin kt kicks, 1-.true., other - .false. 
set hqreddof 0          ! reduction of dofs in EPOS (0:none,1:HBS,2:EPOSpara,3:therm.mass)
set hqctc 1.0           ! cTc = T_max/Tc
set hqxsectioncranck 0  ! cranck Q+Q<->quarkonia cross sect. (0: none/1: standard)
set hqtempcoal 0.154    ! temp for coalescence
set hqifcoalinflcell 1  ! if coal in fluid cell


endoptns hq
