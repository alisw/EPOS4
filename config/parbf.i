set isetcs 1 ! 2 ! (-1 = no rj tables, 2 = complet, 1 = no Pom_sat tables)
set iregge 0
!satpom off         !for testing
!--------------------------------------------------------------
set factk 1.0!1.6!2.4!1.8!1.2     !compensate with epscrx in key
set facq2tim 1.0 !1.5 !1.25 ! !factor for Q2 arguments of timsh2 
set q2fin  0.2 !2.5      
set qkonia1 0. !2.0 !6.0 !1e10 !1000.
set qkonia2 0. !6.0
set vparam(1) -0.35 !0.022 !0.03 !0.04 
set vparam(2) -0.!1 !0.022 !0.03 !0.033 !0.04 
set vparam(3) 0.2 !-0.03 
set vparam(4) 0.15 !0 !0.15 !0.05 
set vparam(5) 0
!vparam(9) currently used as counter!!!
set alpsft 0.90 !0.45
set zrhoinc -0.!03 !-0.055 !-0.3 !0 !-0.6 !-0.3 !-20. !5   
set nsegmincore 10 !20 !50 !30 !10
set ratiomaxcore 2
echo off
!-------------------------------------------------------------
!set iscreen 0
!set iomega 2
!set ishpom 0
!set ifrade 0
!set iokoll 1
!set delh 0.1  
!-------------------------------------------------------------
set ioTestFact 0
set laddTestFact 0
set kfrout 22 !FO method: 1 => T,  2 => epsilon,  22 => epsilon & MiC
set nfrout 10000 !number of eta bins for building subclusters for MiC
!-------------------------------------------------------------
!ev.i : depends on naflav only from alpha_s (psuds and pssalf) (~30mn)
!aa.i : depends on naflav explicitely (psjetj) (~2 mn)
!abxxx.i : depends on naflav explicitely (psjet) (~1h per table)
set naflav 3   !for space-like and remnants
set nbflav 5   !for time-like and born
!-----------    fragmentation: prob_i_j ~ exp{-pi/kappa*[Mi+Mj+Ci*Cj*Mo]^2}
set pmqu 0.003 !0.002  !0.003!* !Mu         
set pmqd 0.006 !0.0045 !0.006!* !Md         
set pmqs 0.07 !0.066 !0.072 !0.074 !0.095  !0.074!* !Ms         
!set pmqq => key !Mo
set pudd 0.95 !0.95  !0.88   !0.88 !* !Cd  (Cu=1) 
set puds 0.57 !0.55 !0.52  !0.47 !0.50  !0.65   !0.50 !* !Cs
  
set ptfra 0.33 !0.33 !0.50!0.33!*          
set ptfraqq 0.33 !0.33 !0.50!0.33!*          

set fkappa 0.015
set fkappag 0.015 !use the same for gluon and quark ... difference probably comes from collective hadronization

set irescl 3  ! 2,3 ... ist=0 particles rescaled     
set factgam 1.

!set s0min 0.1


set iocova 40 !technical parameter, concerning speed 

set bsatur 1. !5. !2.5 !2. !1. !1.5 !1. !0 !1.5 !1. !0.35
set csatur 1. !2. !5. !1. !0.1 !1. !20 !1.5 !1. !1. !0.33  
set asatur 0.33
set esatur 0.003 !0.005 !0.004 !3 !the smaller the larger effect
set dsatur 0.2 !0.26 !0.18 !0.285  !0.33 !0.5 !1. !0.60 !0.80      !overall factor fixed to have <Npom/<Npom>>=1 
set factjpsi 2. !4. !2. !0.25 !0.3 !0.14 !25  !reduction of xsection with light quarks (larger factjpsi means more HQ)
set zmsinc 1.    !intrinsic pt increase with Q2s for pomsat
set facnof 0. !3.5    !reduction of HQ in PDF (OK between 3 and 5)
set iq2sat 2
!set iq2sat 1
!set iq2sat 0 

! gampar increase both low energy ND and diff contribution. If too large, then the contribution of diff term in total cross section is too large even if the diff cross-section is OK because no multiple scattering is included in diff xs.
set gampar 1.28 !1.4 !1.33  !xs at low energy (soft) and multiplicity at SPS

set egyscr 33.   !used to start Z increase for diff
!to avoid "kink" in xs, Z should be saturated from the beginning, so epscrw will change the slope (and ratio elast/inelastic)  by changing the decrease of Z after the saturation
set epscrb 1.5 !1.7 !1.4 !1.35 !0.15  !modify b dependence for diff
set epscrw 2.5 !1.        !minimum for diff
set r2hads(1) 0.65 !1.5 !0.85
set r2hads(2) 0.65 !1.5 !0.85
set epscrp   2.22 !2.17 !2.1 !1.25   !b width param for diff Z    -> pp xsection ............ w_B   !change slope (smaller = larger slope) and elastic/inelastic (larger=larger)

set epscrd -0.2 !-0.24 !-0.33 !-0.24 ! -.2 !-0.24 !-0.255 !-0.097 !-0.053 ! -0.26 !-0.6 !-6.2 !-4. !-1.3 !3.3 !diff contr. (same factor for all not to change rapidity dependance of diff
set epscrh -0. !-1.66 !-1.33 !screening for hard without saturation (relative to epscrx -> change Q2s (larger increase Q2s)
set zbrads 0.6 !0.99     !limit for Phi to define saturation
!-----> znurho moved to paraf.i <-----

set zbrmax 0.5 !0.33 !factor in Phiexpo to define probability to add a pair to nuclear screening
set zdfbmx 0.5 !max diameter to count nucleons (relative to bkmx) for nuclear Z diffractif (fix diffractive component)
                                                     
! sigine includes cross-section due to diff diagram with more than one scattering
set xmxmas  0.5 !0.33 !maximum energy fraction allowed for excited remnant: change SD/DD increase maximum mass (help to get DD at mid-rapidity) and change slope behavior at low energy

set alpdif 0.3 !0.38 !0.35    !technical parameter to avoid bet_diff < -1 (good if=alppar+alpdif~0.6) !increase fraction of low x diag (soft) but should not be too large to avoid wrong fit of hard part at high energy (change sigdd/sigsd and then influence sigcut/sigine (for same sigsd) which may change binary scaling in pp)

set r3pom 0.07 !0.09 !0.33   !0.42 !0.4 !0.52 !0.44 !0.66     !define ratio SD/(CD or DD): larger is less SD

set r2reg(2,1) 0.3
set sloreg(2) 0.7
set alpreg(2)  0.5  ! slope for reggeon in MC 
set gamreg(2,1) 8.    !0.5  !reggeon factor

!################################################################
! IMPORTANT before changing above params :
! VERY difficult to fit at the same time
!   (1) jet pt spectra (anti kt, ATLAS and ALICE)
!   (2) high pt spectra (>8 GeV) charged, pi0, K, Lambda
!   (3) charged particle dn/deta(0)
! Seems impossible without lowering asatur (even with factk=1)
!################################################################

set etaos 0.08 !0.32!0.0001! very little change spectra pp (same zetaos)
set zetaos 0.

set dsegce 1. !11.68   !replace nsegce : dsegce*vocell=nsegce (to have less core a low energy because cell volume is larger so density is lower)
set tfrout 0.164 !0.150 !0.164 !0.168 !0.140   !T~0.150 too flat pt spectra
set fofac  1.480 !1.300 !1.550 !1.650ot !1.700   !in PbPb, pp less sensitive
set epsfin 0.14  !0.20  !0.05
set taustr 1. !0.05 !0.3! 1.0 !1000. !1.0
set iotst2 0 !30 !20 !modifies fpost (0.01fm)

set fxcell 4.0  !  core-corona     !5.0
set fzcell 1.0  !  separation

set ptlow 0.
set ptupp 10000.0
set taurem 1

!-----------------------------------

set volex .5
set cutdxy .0
set fludiq 0.22  !0.25
set tauhac 0

!-----------------------------------


set wdiff(1) 1.
set wdiff(2) 1.
set wdiff(3) 1.
set edifac 1.  !transtion reggeon-diffraction 
set rexres(1) 1.  !pion remnant low mass excitation probability decrease
set rexres(2) 1.   !nucleon low mass remnant excitation probability decrease
set rexres(3) 1.    !kaon low mass remnant excitation probability decrease

!low mass important for LHCf data eta>10
set rexdif(1) 0.1
set rexdif(2) 0.1
set rexdif(3) 0.1
!"normal" mass given by Pomeron intercept (alppom)
!set rexpdif(1) 0.
!set rexpdif(2) 0.6 !0.25
!set rexpdif(3) 0.
!high mass excitation
set rexndi(1) 0.
set rexndi(2) 0.0!8 !3!4
set rexndi(3) 0.
set zodinc    0.0!5 !0.033 !0.1
set alpndi(1)  2.5 !high mass non-diffractive mass index without diquark
set alpndi(2)  2.5 !high mass non-diffractive mass index with diquarkpa
set alphigh(1) 1.14 !mass for multiple scattering !change ALICE scaling ???
set alphigh(2) 1.14 !mass for multiple scattering
set rexddf 0.75 !0.25  !transfer between CD and DD when energy is too low

 
set ammsqq 0.5    !Minimum mass for a reggeon
set cumpox 4.   !Minimum mass for a Pomeron (not only for soft diffraction) -> at low energy, transform low pass pomeron into remnant excitations (adjust multiplicity and pseudorapidity shape)
set zdfinc 0.    !maximum mass increase as a function of x_R (the larger the narrower)
set zdrinc 0.    !droplet minimum mass increase as a function of Q2
set xcupom 1.
set amdrmax 5.e1    !max mass of remnant given per pomeron
set amdrmin 0.      !mass added to minimum mass if MPI
set irmdrop 1
set alpdro(1) 1.e10   !additional mass for droplet to split remnant
set alpdro(2) 0.01    !pt of droplet
set alpdro(3) 1.      !mass distribution of droplet

set xmindiff 1.  !factor on minimum energy of remnant for soft diffraction
set xminremn 1. !1.35 !1.5 !1.5 !1.35 !2.         !factor on minimum energy of remnant (larger than 1 to produce enough anti-baryon baryon at NA49) 
set delrex 0.5
set exmass 0.5    !large enough in case of strangeness in string end


!----------------------------

set nclean 0

!----------------------------

set facmc 4.    !to enlarge b range for xs calculation
set qcdlam 0.038 !0.033 !0.035 !0.04  !if too small (~0.025), F2 do not rise fast enough. Increase of qcdlam, increase alpha_s, reduce factk, and decrease charm from SL ?

set facdif 1.01
!if too large, rapidity gap is not good for ND


set gamhads(1) 1.075 !0.96 !0.85 !0.75!0.55
set gamhads(3) 0.95 !0.87 !0.47 !0.33!0.47



!################################################################
! balance/multiplicite section-efficace : slopom (et/ou r2had)
!  largeur de l amplitude \  multiplicite / section efficace \
!       put   r2part = r2had
! This possibility NOT USED in current tune 
!################################################################


set slopom 0.02 !0.02 !0.03 !0.13 !2 !0.1 !0.1 !0.05 !275 !0.1    !slope of the "slope", ratio elastic/inelastic and multiplicity at high energy. Should be used to fix these observables but not the slope which is changed by epscrb or epscrp
set r2had(1) 0.65     !to get correct hA total cross section, r2had should be larger 
set r2had(2) 2.15 !1.9 !2.15 !2.25 !1.25 !1.3      !for p than for pi. Change ratio elastic/inelastic and multiplicity 
set r2had(3) 0.8      !(lower r2had => higher mult (more MPI)  and higher sigela/sigine at low energy)
set r2part 0.8   !minimum rh

!to get proper scaling behavior (in average and vs multiplicity), the Q2s distribution (vs b) is very important. Since Q2s comes from the ratio between G_QCD and G_fit, the combination betpom/alppar/alplea should match F2 AND give proper Q2s in pp.


set zzsoft -0.04 !0.35 !0.04 !for gluon in SoftSat (initial condition for Q2s>Q20) (<0 decrease gluon fraction, 0< increase muon faction)
set zoeinc 1. !0.35 !3.3 !1.85 !2.5 !0.25     !gobal reduction in Softsat
set q2nmin  1.50 !1.1  !2.5 !2.     !fix ratio soft/hard (large q2nmin=less hard=lower multiplicity and softer pt spectra and narrower pseudorapidity width)
set q2sft 1.50 !1.1  !1.5      !limit to have perturbative calculation (born)
set glusea 0.1 !0.14 !0.09 !0.12 !0.16 !0.17 !0.22 !0.16 !0.175 !0.185 !0.175 !0.25 !0.15 !0.6 !0.02              !not too small not to have to many gluons
set betpom 0.66 !0.33 !0.2 !0.66 !1.0 !0.2 !0.2 !1. !0.6 !0.7 !0.66 !0.5 !0.15!15 !0.666 !0. !1.  !0.25   !controls F2 around x=0.1 and normalization at low x and low pt jet inc xs
set alppar 0.33 !0.5 !0.25 !0.15 !0.33 !0.5 !0.55 !0.35 !0.7!85 !0.75 !0.55 !0.4      !change slope of F2 at low x (large=>fast rise) AND sigma_ND but not sigma_diff (because of minimum exponent=-1) !
set alpparh 0.1 !value of alppar for large Q2s
set gamtil 1. !0.75 !0.08             !overall for soft preevolution (should be 1 )
set alppom  1.14 !1.1 !1.12 !1.075 !1.15 !1.075 !1.09                  !(>1) controls slope of F2 at low x and low Q2

 ! large value of alplea (~1.5 better for SPS and LHCf data but lower value make pseudorapidity distribution larger at LHC)
set alplea(1) 0.15   !change slope of F2 at all x (large=high slope) 
set alplea(2) 1. !1. !1.33 !0.66       !1.5 !1.66     !and xf distributions of hadrons at large xf 
set alplea(3) 0.15   !(large=less particle at high x) and MPI (small=small MPI)
                               ! and slope of inclusive jet distribution (larger=more low pt and less high pt) and Pomeron mass (smaller = higher mass)
                                                                                                              !alplea~2 not good for multiplicity distribution

set r2har(2) 0.05    !should be small (<<r2had)

!set wgtval 0.  !warning if not 0, it reduces the effect of wgtqqq (only valence diquark if available)
set wgtqqq(1) 0. !0.14
set wgtqqq(2) 0.
set wgtqqq(3) 0. !0.14
set wgtdiq 0.08!15
set reminv 0.2
set rstras(2) 0.25 !0.9
set diqcut 0.0! 0.001
set strcut 0.0!1


set ptsend 0.33 !0.1 !factor on random component  !change pt distributions of p and ap at large eta at RHIC
set ptsems 0.66  ! <qmass> for minimum  !play on <pt> at low energy

!yradpx should not be large but on many particles. Id particles spectra indicates a moderate flow but to get <pt> right enough particles should changed
set amuseg 3. !1. !2.5 !5. !1.25     !normalization for density used in radial flow
set ydsrd 1.!85 !0.25 !0.05 !0.45 !2.5  !factor for asymmetry for rad flow (pA)
set ydslg 1. !2. !1.5 !3.6 !1.6 !0.9 !!25.    !factor for asymmetry for long flow (pA)
set yradpx 0.21  !0.15 !1.5 !1.05 !2. !0.5   !geo (decrease flow in AA)
set yradmx 0.5 !0.42 !0.4 !0.32 !0.35 !0.33 !0.4 !0.212 !0.26 !0.255!0.37 !0.31 !0.3  !factor for radial flow
set yradpp 1. !1.7 !1.1 !0.0 !8 !0.2  !minimum radial flow

set ylongmx -0.75 !0.42 !-0.85 !-2.5 !-0.85 !-0.75 !-0.25 !3. !1.3 !-0.55 !-0.51 !-1.15 !-0.1 !-0.9  !max long collective boost ( < 0 -> take from jintpo ) for nuclear collisions
            !< 0 : delzet=ylongmx*log(density),density=ectot/, ectot = energie du core complet (TP/ind.f)
            !> 0 : delzet=ylongmx

set yradmi  0.5 !0.4 !0.6 !0.66 !0.45 !0.64 !0.433 !0.4 !0.45 !0.45 !0.58 !0.42 ! 0.0022 !factor for long flow
set yradpi  1. !0.5 !0.45 !1. !-3. !0.08  !0.23 !geo (increase mult in pA : reduced flow for large Npom)
set iocluin 1
set facecc  0.12 !0.18
set ijetfluid 1
set epscri(1) 0.35   !0.15

echo on
