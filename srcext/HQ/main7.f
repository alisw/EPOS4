      module JA
      integer jancoll,jahqcoll,jaityp,jaipart
      end module

      module JA2
      integer jaish
      end module

      module PBGdistccbar
      integer nbdist,idist,nbccbarindatabase
      parameter (nbdist=4)
      end module

      module PBGdistccbarbis
      integer typedistrel
      double precision paramdistrel(2)
      end module

      module for2body
      implicit none
      integer if2bodySC
      logical correctshift,ifsingletred
      end module

      module PBGfragmentc
      integer itypfrag
      integer itypcoal
      double precision peterepsq(2),peterepsh(2),peterzmax(2),
     &                 petermaxval(2)
      end module

      module PBGimparam
      implicit none
      real b
      logical ifaxial
      end module

      module PBGinfotr
      integer infotrack
      end module
      
      module PBGnpartmax
      implicit none
      integer npartmax
      parameter (npartmax=10000)
      end module

      module onerateHQparam
      integer levelcorona
      double precision timein,tlimRemler,timestepwig,tempdissocRem
      logical ifinterpolation
      double precision schemeRemler
      end module

      module PBGopenflavor
      implicit none      
      double precision mD,mDs,mlambdac,mxic,momegac
      double precision mBmes,mBs,mlambdab,mxib,momegab
      end module

      module PBGnres
      implicit none
      integer nres,nresc
      parameter (nres=26)
      parameter (nresc=13) 
      end module

      module PBGpsiinfo
      double precision mc,mc2,muc,mpsi,eps0c,a0c
      end module

      module PBGupsinfo
      double precision mb,mb2,mub,mups,eps0b,a0b
      end module

      module PBGpsiinfo2
      implicit none
      integer mcchoice
      end module


      module PBGalltaa
      integer nbofb
      double precision alltaa(0:20)
      end module

      module PBGallvariance
      real allsigma(0:20),allsigmax(0:20),allsigmay(0:20)
      end module

      module PBGcoupling
      double precision alphasforrad
      end module

      module PBGblckhadronize
      implicit none
      integer itypDproduct
      double precision pcrithadr(2)
      end module

      module PBGblckifptcles
      logical ifptcles
      end module

      module PBGblckinitptspect
      implicit none
      integer nptbinmax,nptbin(2)
      parameter (nptbinmax=1000)
      double precision ptstep(2),nofpt(2,0:nptbinmax),
     &  nofptAA(2,0:nptbinmax)
      double precision tbcoeffM(2,0:nptbinmax)
      end module

      module PBGblckinitptspect2
      implicit none
      double precision ptmax(2)
      end module

      module PBGblckinitptspect3
      implicit none
      double precision globshadow(2)
      end module

      module PBGblcktransi
      implicit none
      double precision temptrans,entrans_min,entrans_max
      end module

      module PBGbookpsi
      implicit none
      integer nboffdiag,nbdiag
      end module

      module PBGcharge
      implicit none
      integer ncharged
      end module

      module PBGcivirtual2
      integer ivirtual2
      end module

      module PBGcomdndsigma
      use PBGnres
      integer nbetamax
      parameter(nbetamax=30)
      double precision tbdndsigmasl(0:nresc,-nbetamax:nbetamax),
     & basedndsigmatl(0:nresc),stepeta
      end module

      module PBGDchoice
      integer idmeson
      end module

      module PBGevolQ
      implicit none
      integer imedfinal
      logical ifdecrease
      end module

      module PBGevcount
      real eventsxx,eventshq
      integer ibev,jaeve(11),jaglcol(11)
      data eventshq/0/
      end module

      module PBGforboltzmann
      integer boltz_model
      end module

      module PBGforcovdecay2
      integer ifail
      double precision x0(0:3),temp,p(0:3)
      end module

      module PBGfordebyeeff
      integer nbtemp,nbtempmax
      parameter(nbtempmax=200)
      double precision tempmin,tempmax,deltemp,tbdebye(nbtempmax)
      end module

      module PBGfordecaydat1
      double precision tempmin,tempmax,pmax
      end module

      module PBGfordecaydat2
      integer nbtemp,nbtempmax,nbp,nbpmax
      parameter (nbtempmax=100,nbpmax=100)
      double precision delttemp,deltp,decaygrid(0:nbpmax,0:nbtempmax)
      common/fordecaydat2/decaygrid,delttemp,deltp,nbp,nbtemp
      end module

      module PBGforfpdat1
      double precision tempfmin(2),tempfmax(2),pfmax(2),fpmult
      end module

      module PBGforfpdat2
      integer nbtemp(5),nbtempmax,nbp(5),nbpmax
      parameter (nbtempmax=60,nbpmax=300)
      double precision delttemp(2),
     & deltp(2),fpa(2,0:nbpmax,0:nbtempmax),
     & fpbl(2,0:nbpmax,0:nbtempmax),
     & fpbt(2,0:nbpmax,0:nbtempmax)
      end module

      module PBGforkerneldissoc
      integer ires
      double precision temp2,p2,mass
      end module

      module PBGforleptons
      integer nblepton, nbmaxlept
      parameter(nbmaxlept=100)
      integer typlepton(nbmaxlept),origlepton(nbmaxlept)
      double precision tbplepton(nbmaxlept,0:3)
      end module


      module PBGforprobcoal
      integer nbp(5),nbpmax,nbpmaxH
      parameter (nbpmax=300)
      parameter (nbpmaxH=120)
      double precision pmax(5),deltp(5),probcoal(2,0:nbpmax),
     & probcoalHad(2,5,0:nbpmax)
      end module

      module PBGforprobcoalTL
      integer nbd(2),nbyds(2),nbdmax,nbydsmax
      parameter (nbdmax=80,nbydsmax=55)
      double precision dcoalmin(2),dcoalmax(2),deltdcoal(2),ydsmin(2),
     & ydsmax(2),dyds(2),probcoalTL(2,0:nbdmax,0:nbydsmax)
      end module

      module PBGforprocesstype
      logical ifcol,ifrad
      end module

      module PBGforquantDBraat
      integer nbpbmax,nbpb
      parameter (nbpbmax=100)
      double precision zpbDP(0:nbpbmax),zpbDV(0:nbpbmax),dpb
      end module

      module PBGforquantDkart
      integer nbpbmax,nbpb
      parameter (nbpbmax=200)
      double precision zpbDK(0:nbpbmax),dpb
      end module

      module PBGforradiatGB
      logical iffullQCD
      end module

      module PBGforradiatLPM
      logical iflpm,ifdamp
      integer typelpmdamp
      end module

      module PBGforratedat1
      double precision ratemult,tempmin(2),tempmax(2),pmax(2)
      end module

      module PBGforratedat2
      integer nbtemp(2),nbtempmax,nbp(2),nbpmax
      parameter (nbtempmax=60,nbpmax=300)
      double precision delttemp(2),deltp(2),rf(2,0:nbpmax,0:nbtempmax),
     & rg(2,0:nbpmax,0:nbtempmax)
      end module

      module PBGforratedat3
      double precision:: redthmass(1:22,1:80),redthmassgl(1:22,1:80)
      end module

      module PBGforrateratiologdat1
      double precision tempminrad(2),tempmaxrad(2),
     & logpminrad(2),logpmaxrad(2)
      end module

      module PBGforrateratiologdat2
      integer nbtemprad(2),nbprad(2),nbtempmax,nbpmax
      parameter (nbtempmax=60,nbpmax=300)
      double precision delttemprad(2),
     & deltprad(2),rtildef(2,0:nbpmax,0:nbtempmax),
     & rtildeg(2,0:nbpmax,0:nbtempmax)
      end module

      module PBGforrateratiooneqlogdat2
      integer nbtemprad(2),nbprad(2),nbtempmax,nbpmax
      parameter (nbtempmax=60,nbpmax=300)
      double precision delttemprad(2),
     & deltprad(2),rtQCD(2,0:nbpmax,0:nbtempmax),
     & rtSQCD(2,0:nbpmax,0:nbtempmax)
      end module

      module PBGforratioprbferm
      integer nbmu,nbmumax,nblns,nblnsmax
      parameter (nbmumax=50,nblnsmax=150)
      double precision mumin,mumax,deltmu,lnsmin,lnsmax,deltlns,
     & ratioprobferm(2,0:nblnsmax,0:nbmumax)
      end module

      module PBGforratioprbglu
      integer nbmu,nbmumax,nblns,nblnsmax
      parameter (nbmumax=50,nblnsmax=400)
      double precision mumin,mumax,deltmu,lnsmin,lnsmax,deltlns,
     & ratioprobglu(2,0:nblnsmax,0:nbmumax,5)
      end module

      module PBGfpcoeff
      integer fpchoice,tunchoice
      end module

      module PBGglaubercentrality
      integer glauber_0_20,glauber_20_40,glauber_40_60,glauber_60_100
     &                         glauber_40_60,glauber_60_100
      end module

      module PBGgluprops
      double precision:: c1,c2
      end module

      module PBGhqdmesons
      integer nd
      double precision, dimension(:,:), ALLOCATABLE :: hqmeson 
      end module

      module PBGhqprod
      logical ifhqfromhq
      end module

      module PBGhquarks
      integer nHQ
      integer, dimension(:), ALLOCATABLE :: idhquarks
      real, dimension(:), ALLOCATABLE :: mhquarks
      real, dimension(:,:), ALLOCATABLE :: phquarks
      end module

      module PBGievent
      integer iev_0_20,iev_20_40,iev_40_60,iev_60_100
      end module

      module PBGifinputfile
      logical ifinput
      end module

      module PBGifhq2
      integer ihq2
      end module

      module PBGifskip
      integer iskip
      end module

      module PBGimpparevent
      integer nbsubrunmax,nbimpactmax
      parameter (nbimpactmax=15,nbsubrunmax=100)
      integer nbevent(nbsubrunmax,0:nbimpactmax)
      end module

      module PBGiniphasespace
      integer numberofHQ
      double precision  HQarray(1000,10) 
      integer HQmother(1000)
      end module

      module PBGinitnb
      double precision avnc0ncoll,avnpsi0ncoll,avnb0ncoll,avnups0ncoll
      end module

      module PBGinitnbmult
      implicit none
      double precision avnc0ncollmult,avnpsi0ncollmult,avnb0ncollmult,
     &  avnups0ncollmult
      end module

      module PBGinixypos
      implicit none
      double precision  xyarray(5000,0:3)
      integer  numberofNNColls, numberelemcolls
      end module

      module PBGmainpro1
      implicit none
      logical ifcronin
      integer idistspatial,nbcol,nbsubrun,nbcolsaveEPOS
      end module

      module PBGmainpro2
      implicit none
      real b_inf,b_sup
      double precision tauf,timestep
      end module

      module PBGnbimpapara
      implicit none
      integer nbimpact
      end module

      module PBGnumbercoll
      implicit none
      integer weight
      end module

      module PBGnumberevent
      implicit none
      integer totcol
      end module

      module PBGnumbereventHQ
      implicit none
      integer icol
      end module

      module PBGpartinmedium
      implicit none
      double precision tempdissoc(26)
      end module

      module PBGpawc
      implicit none
      integer nwpawc
      parameter (nwpawc=1000)
      real hmemor
      end module

      module PBGproductionmech
      implicit none
      integer hhbarprod
      end module

      module PBGran2junk
      implicit none
      integer NTAB
      parameter(NTAB=32)
      integer*4 IDUM,IDUM2,IY,IV(NTAB)
      end module

      module PBGreductiondofs
      implicit none
      integer reddof
      double precision cTc
      end module

      module PBGsigmahqprod
      implicit none
      double precision sigprodc,sigprodb,sigprodpsi,sigprodups
      end module

      module PBGspectra
      implicit none
      logical ifextspectra(2)
      integer typeextspectra(2)
      logical ifshadowing(2)
      integer typeshadow(2)
      end module

      module PBGsystem
      implicit none
      integer A
      double precision sqrts
      end module

      module PBGsystem2
      implicit none
      integer npt,ngl,kol
      real bim,phi,phir
      end module

      module PBGtointegdndsigma
      double precision alpha,n0,nl
      end module

      module PBGtomsqfermrunning1b
      double precision m,s,mu
      end module

      module PBGtomsqfermrunning2b
      double precision m,s,mu
      end module

      module PBGtomsqglurunning1b
      double precision m,s,mu
      end module

      module PBGtomsqglurunning2b
      double precision m,s,mu
      end module


      module PBGtransfer
      implicit none
      double precision hqtimestep,hqcTc,hqXsectioncranck,hqratemult,
     &   hqc2,hqfpmult,hqalphasforrad,hqhighptcut,hqtempcoal
      integer hqnbcol,hqnbsubrun,hqtype_evol_Q,hqimedfinal,hqreddof,
     &        hqitypDproduct,hqitypfrag,hqitypcoal,hqboltz_model,
     &        hqtypelpmdamp,hqfpchoice,hqtunchoice,hqifdecrease,hqcoll,
     &        hqifcronin,hqifcol,hqifrad,hqiffullQCD,hqiflpm,hqifdamp,
     &        hqifptcles,hqifhighptselect,hqish,hqifhqfromhq,
     &  hqifextspectra1,hqtypeextspectra1,hqifshadowing1,hqtypeshadow1,
     &  hqifextspectra2,hqtypeextspectra2,hqifshadowing2,hqtypeshadow2,
     &  hqifcoalinflcell,hqnbcolsaveEPOS,hqifallowhqskip,hqif2bodySC,
     &  hqifsingletred,hqdpmaxabscorona,hqtypedistrel,hqifbselect,
     &  hqparamdistrelc,hqparamdistrelb


      logical iftranstimestep,iftranscTc,iftransXsectioncranck,
     & iftransratemult,iftransc2,iftransfpmult,iftransalphasforrad,
     & iftranshighptcut,iftransnbcol,iftransnbsubrun,iftranstype_evol_Q,
     & iftransimedfinal,iftransreddof,iftransitypDproduct,
     & iftransitypfrag,iftransitypcoal,iftransboltz_model,
     & iftranstypelpmdamp,iftransfpchoice,iftranstunchoice,
     & iftransifdecrease,iftranscoll,iftransifcronin,iftransifcol,
     & iftransifrad,iftransiffullQCD,iftransiflpm,iftransifdamp,
     & iftransifptcles,iftransifhighptselect,iftransish,
     & iftransifhqfromhq,iftransnbcolsaveEPOS,
     &  iftransifextspectra1,iftranstypeextspectra1,iftransifshadowing1,
     &  iftranstypeshadow1,iftransifextspectra2,iftranstypeextspectra2,
     & iftransifshadowing2,iftranstypeshadow2,iftranstempcoal,
     & iftransifcoalinflcell,iftransifhqskip,iftrans2bodySC,
     & iftranssingletred,iftransdpmaxabscorona,iftranstypedistrel,
     & iftransifbselect,iftransparamdistrelc,iftransparamdistrelb 


      logical ifinittimestep,ifinitcTc,ifinitXsectioncranck,
     & ifinitratemult,ifinitc2,ifinitfpmult,ifinitalphasforrad,
     & ifinithighptcut,ifinitnbcol,ifinitnbsubrun,ifinittype_evol_Q,
     & ifinitimedfinal,ifinitreddof,ifinititypDproduct,
     & ifinititypfrag,ifinititypcoal,ifinitboltz_model,
     & ifinittypelpmdamp,ifinitfpchoice,ifinittunchoice,
     & ifinitifdecrease,ifinitcoll,ifinitifcronin,ifinitifcol,
     & ifinitifrad,ifinitiffullQCD,ifinitiflpm,ifinitifdamp,
     & ifinitifptcles,ifinitifhighptselect,ifinitish,
     & ifinitifhqfromhq,ifinitnbcolsaveEPOS,
     &  ifinitifextspectra1,ifinittypeextspectra1,ifinitifshadowing1,
     &  ifinittypeshadow1,ifinitifextspectra2,ifinittypeextspectra2,
     & ifinitifshadowing2,ifinittypeshadow2,ifinittempcoal,
     & ifinitifcoalinflcell,ifinitihq2,ifinitiskip,ifinitivirtual2,
     & ifinitidmeson,ifinithqskip,ifinit2bodySC,ifinitsingletred,
     & ifinitdpmaxabscorona,ifinittypedistrel,ifinitifbselect,
     & ifinitparamdistrelc,ifinitparamdistrelb

      end module

      module PBGvariancesig
      implicit none
      real sigma,sigmax,sigmay
      end module

      module PBGtrigger
      logical ifhighptselect
      logical ifbselect
      real highptcut
      end module

      module PBGXsection
      implicit none
      double precision Xsectioncranck
      end module

      module PBGfragandcoal
      use PBGblckhadronize
      use PBGfragmentc
      implicit none
      logical ifcoalinflcell

      double precision lambdaD2,lambdaB2,mqconstit,mqconstit2,alphadov,
     & tempcoal,msconstit
	
      logical ifsavehadro

      double precision nbhadro,nbhadroofpt(200,-1:1),
     &  avprobcoalofpt(200,-1:1),varprobcoalofpt(200,-1:1) 
      double precision nbmesonsofpthadrocoal(200,-1:1),
     & Haahadro(200,-1:1),Haahadrocoal(200),nbmesonsofpthadrofrag(200,
     & -1:1)
      double precision distribproba(50,-1:1)
      end module


      module PBGfordisplay
      use PBGnres
      implicit none
      character*50 hbkname
      character*10 ycentname
      integer ndivpsimax,ndivpsi,ndivrap1dmax,ndivrap1d,
     & ndivpt1dmax,ndivpt1d,ndivphi1d,ndivphi1dmax
      parameter(ndivpsimax=20,ndivrap1dmax=100,ndivpt1dmax=100,
     &          ndivphi1dmax=100) 
       double precision boundjpsiinf,boundjpsisup,incbas,boundrapinf,
     & radstep,boundrapsup,steprap1d,incrap1d,boundptsup,steppt1d,
     & incpt1dy0,incpt1dy2,boundphisup,stepphi1d,incphi1d,incphipt,
     & rapsupobs,
     & packphi(0:nres+4,ndivphi1dmax),
     & packrap(0:nres+4,ndivrap1dmax),
     & packpty0charm(ndivrap1dmax),
     & packpty0dmeson(ndivrap1dmax),
     & packpty0(0:nres+4,ndivpt1dmax),
     & packpty0_0_20(0:nres+4,ndivpt1dmax),
     & packpty0_20_40(0:nres+4,ndivpt1dmax),
     & packpty0_40_60(0:nres+4,ndivpt1dmax),
     & packpty0_60_100(0:nres+4,ndivpt1dmax),
     & packpty2(0:nres+4,ndivpt1dmax),
     & packpty0bis(0:nres+4,ndivpt1dmax),
     & packpty0bis_0_20(0:nres+4,ndivpt1dmax),
     & packpty0bis_20_40(0:nres+4,ndivpt1dmax),
     & packpty0bis_40_60(0:nres+4,ndivpt1dmax),
     & packpty0bis_60_100(0:nres+4,ndivpt1dmax),
     & packpty2bis(0:nres+4,ndivpt1dmax),
     & packpty0nonpromptJpsi(ndivpt1dmax),
     & packptphiy0nonpromptJpsi(ndivpt1dmax,ndivphi1dmax),
     & packmeanptsq(0:nres+4,ndivrap1dmax),
     & packptphiy0(0:nres+4,ndivpt1dmax,ndivphi1dmax),
     & packptphiy0bis(0:nres+4,ndivpt1dmax,ndivphi1dmax),
     & packcorrelphi(0:nres+4,0:6,ndivphi1dmax),
     & packpty0alongb(0:nres+4,ndivpt1dmax),
     & packnbpsi(0:ndivpsimax),
     & packcorrelpty0(0:nres+4,ndivpt1dmax,ndivpt1dmax,3)


      double precision packvacrap(0:nres+4,ndivrap1dmax),
     & packvacmeanptsq(0:nres+4,ndivrap1dmax),
     & packvacpty0(0:nres+4,ndivpt1dmax),
     & packvacpty0charm(ndivrap1dmax),
     & packvacpty0dmeson(ndivrap1dmax),     
     & packvacpty0bis(0:nres+4,ndivpt1dmax),
     & packvacpty0bis_0_20(0:nres+4,ndivpt1dmax),
     & packvacpty0bis_20_40(0:nres+4,ndivpt1dmax),
     & packvacpty0bis_40_60(0:nres+4,ndivpt1dmax),
     & packvacpty0bis_60_100(0:nres+4,ndivpt1dmax),
     & packvacpty2(0:nres+4,ndivpt1dmax),
     & packvacpty2bis(0:nres+4,ndivpt1dmax)


      double precision packvacrap_zero(0:nres+4,ndivrap1dmax),
     & packvacmeanptsq_zero(0:nres+4,ndivrap1dmax),
     & packvacpty0_zero(0:nres+4,ndivpt1dmax),
     & packvacpty0bis_zero(0:nres+4,ndivpt1dmax),
     & packvacpty2_zero(0:nres+4,ndivpt1dmax),
     & packvacpty2bis_zero(0:nres+4,ndivpt1dmax)
     & ,packvacpty0nonpromptJpsi(ndivpt1dmax)


      double precision packvacrap_cronin(0:nres+4,ndivrap1dmax),
     & packvacmeanptsq_cronin(0:nres+4,ndivrap1dmax),
     & packvacpty0_cronin(0:nres+4,ndivpt1dmax),
     & packvacpty0bis_cronin(0:nres+4,ndivpt1dmax),
     & packvacpty2_cronin(0:nres+4,ndivpt1dmax),
     & packvacpty2bis_cronin(0:nres+4,ndivpt1dmax)

                    
      integer ndivtmax,ndivr1dmax
      parameter(ndivtmax=150,ndivr1dmax=40)
      double precision timestepforsave,timesave(0:ndivtmax),
     &  avdnbpsidyoftime(0:ndivtmax,-5:5)

      double precision dptdroftimey0(0:ndivtmax,0:ndivr1dmax)
      double precision dptdroftime(0:ndivtmax,0:ndivr1dmax)
      double precision dncdroftimey0(0:ndivtmax,0:ndivr1dmax)
      double precision dncdroftime(0:ndivtmax,0:ndivr1dmax)
      double precision dncdptdyvac(ndivpt1dmax,ndivrap1dmax)
      double precision dncdptdyoftime(0:ndivtmax,ndivpt1dmax
     & ,ndivrap1dmax)
      double precision dncdptoftimey0(0:ndivtmax,ndivpt1dmax)
      double precision dncdptoftime(0:ndivtmax,ndivpt1dmax)
      end module


      module PBGgenvar
      use PBGnpartmax
      use PBGnres
      use PBGopenflavor
      use PBGpsiinfo
      use PBGpsiinfo2
      use PBGupsinfo
      use PBGifinputfile
      use PBGinfotr
      implicit none
      integer neventmax
      double precision hbarc
      parameter (hbarc=0.197323)
      parameter (neventmax=5000)
      integer tolerance
      logical ifallowhqskip
     
      integer g(0:3,0:3),nlor(0:3)



      integer partinfo_ires(npartmax)
      integer partinfo_next(npartmax)
      integer partinfo_mother(npartmax)
      double precision partinfo_r(npartmax,0:3)
      double precision partinfo_p(npartmax,0:3)
      double precision partinfo_rghost(npartmax,0:3)
      double precision partinfo_pghost(npartmax,0:3)
      double precision partinfo_mass(npartmax)
      double precision partinfo_tb(npartmax)
      integer partinfo_whichmedium(npartmax)
      logical partinfo_corona(npartmax)
      integer firstpart
      integer firstemptypart
      integer npart


      integer eventinfo_scatcode(neventmax)
      integer eventinfo_i(neventmax)
      integer eventinfo_j(neventmax)
      integer eventinfo_next(neventmax)
      double precision eventinfo_r(neventmax,0:3)
      double precision eventinfo_tb(neventmax)
      double precision eventinfo_tord(neventmax)
      double precision eventinfo_delt(neventmax)
      integer firstevent,itypeorder
      integer firstemptyevent,nbccbaremission,nbbbbaremission
      integer n_newevent(neventmax)
      double precision time_newevent(neventmax)


      logical resinfo_ifquark(0:nres)
      integer itypquark(0:nres)
      integer resinfo_ccharge(0:nres)
      integer resinfo_spin(0:nres)
      integer resinfo_degen(0:nres)
      double precision resinfo_mass(0:nres)
      integer resinfo_dissoc(0:nres)
      integer resinfo_number(0:nres)
	integer resinfo_antipart(0:nres)
      

      integer ifreac(0:nres,0:nres)
      double precision biggestb2(0:nres,0:nres)
      double precision thresreac(0:nres,0:nres)

      integer type_evol_Q
      end module

      module PBGsceplas
      use PBGdistccbar
      implicit none

      integer nbsce
      parameter (nbsce=6)
      integer typesce(nbsce),isce
      data typesce /1,1,2,2,3,3/

      double precision tau0,li0,ri0,vr0,pbar,vl0,tempinit,alphaprof,
     &  sdens 
      
      double precision temptrans,entrans_min,entrans_max
      end module

      module PBGhydro
      use PBGimparam
      implicit none
      integer nbtbmax,nbrmax,nbtb,nbr
      integer nbxmax,nbymax,nbx,nby
      parameter(nbtbmax=120,nbrmax=70)
      parameter(nbxmax=70,nbymax=70)
      double precision tbmin,tbmax,tbstep,rmax,rstep
      double precision xmax,ymax,zmax,xstep,ystep,zstep
      double precision xmin,ymin,zmin
      integer ntauhyx,nxhyx,nyhyx,nzhyx
      real, dimension(:,:), ALLOCATABLE :: velhyd,enerhyd,temphyd
      real, dimension(:,:,:), ALLOCATABLE :: vx,vy,eps,tem

      end module

      Module denysvar
      integer nbhqc,nbtausteps,nbtausteps0,nbhqincorona,nbhqb,nbhq
      logical, dimension(:), ALLOCATABLE :: ifcorona,ifcoronac,
     &     ifcoronacbar
      double precision taustep
      integer, dimension(:), ALLOCATABLE :: typquarkdyn 
      integer, dimension(:), ALLOCATABLE :: hqmother 
      double precision, dimension(:), ALLOCATABLE :: timelimdyn 
      integer, dimension(:,:), ALLOCATABLE :: iresdyn
      double precision, dimension(:,:,:), ALLOCATABLE :: xdyn,pdyn 
      double precision, dimension(:,:), ALLOCATABLE :: tempatHQpos
      integer, dimension(:,:), ALLOCATABLE :: imedatHQpos
      logical, dimension(:,:), ALLOCATABLE :: ifsinglet,ifsingletfull
       double precision, dimension(:), ALLOCATABLE :: avnbclose
      end module

      Module denysvar2
      integer nstept,nsteppt,nphistep,nystepp
      double precision ptstep,phistep2,deltap,phistep,ptmax,yinf,
     & ysup
      double precision, dimension (:), allocatable :: probq,dn1d,rate1d,
     & loss1d,dn1ddiag,rate1ddiag,loss1ddiag,localrate1d,
     & localratediag1d,HQclosept,hqclose,ptransfspec,ratehigherpT
      double precision, dimension (:,:), allocatable ::dn2d,rate2d,
     & loss2d,dn2ddiag,rate2ddiag,loss2ddiag,localrate2d,
     & localratediag2d,averagesigma
      double precision, dimension (:,:,:), allocatable :: dn3d,rate3d,
     & localrate3d,dn3ddiag,rate3ddiag
      double precision, dimension (:,:,:), allocatable :: phasespdens,
     & phasespdenshighpT
      double precision dncdy,dncdycorona,dnbdy,dnbdycorona
      double precision dnphidycorona,dncpartnercorona
      double precision distribtempattherm(100)
      integer nperct
      end module

      module genuineRemler
      double precision dpmaxrel,dpmaxabs
      double precision dpmaxrelcorona,dpmaxabscorona
      integer nbhardcol,nbhardcolmax
      parameter(nbhardcolmax=60000)
      integer, dimension (:), allocatable :: ipartcol
      double precision, dimension (:), allocatable :: tempcol,timecol
      double precision, dimension (:,:), allocatable :: pcolafter
      double precision, dimension (:,:), allocatable :: dpcol
      end module

      module analcol
      integer iselec,selectedparticles(10)
      double precision, dimension(:,:,:), ALLOCATABLE :: dpoftimeinfo
      end module

      module PSdens
      double precision, dimension (:,:,:), allocatable :: PSdensnative,
     & PSdensnativehighpT
       double precision, dimension(:), ALLOCATABLE :: avnbclosewigner
      end module

      module initnb
      double precision avnc0ncoll,avnpsi0ncoll,avnb0ncoll,avnups0ncoll
      end module 

      module forydistr
      logical ifyat0
      end module


      subroutine ifpartinselect(ipart,ifin,index)
      use analcol
      implicit none
      integer ipart,iscan,index
      logical ifin
      ifin=.false.
      index=0
      do iscan=1,iselec
         if(selectedparticles(iscan).eq.ipart) then
            ifin=.true.
            index=iscan
         endif
      enddo
      end
      


      subroutine hq7
      use PBGhqdmesons
      use PBGhquarks
      use PBGnumberevent
      use PBGhqprod
      use PBGmainpro1 
      use PBGmainpro2
      use PBGinitnb
      use PBGfragandcoal
      use PBGgenvar
      use PBGsceplas
      use PBGhydro
      use PBGevolQ
      use PBGimpparevent
      use PBGnbimpapara
      use PBGbookpsi
      use PBGnumbercoll
      use PBGcharge
      use denysvar
      use genuineRemler
      use onerateHQparam
      use initnb
      implicit none
      logical iffirst
      integer ieventEPOS
      data iffirst/.true./
      save iffirst
      save ieventEPOS
      integer ifail
      integer icol,i,j,ifailcol,itime
      integer nbtry
      double precision timelim
      
      if(iffirst) then
         iffirst=.false.
         ieventEPOS=1
      else
         ieventEPOS=ieventEPOS+1
      endif
      if(ncharged.eq.-1) return
      write(*,*) 'hq7 called'



      nbdiag=0
      nboffdiag=0
      write(*,*) 'allocating room'
      allocate(idhquarks(npartmax))
      allocate(mhquarks(npartmax))
      allocate(phquarks(npartmax,0:3))
      allocate(hqmeson(npartmax,9))
      do i=1,npartmax
        do j=1,9
          hqmeson(i,j)=0
        enddo
      enddo
      nHQ=0
      nd=0
      ifail=0
      if(ifhqfromhq) call initplasmascenario_of_b
      call gettauzer(tau0)
      write(*,*) 'tau therm=',tau0
      icol=0
      nbtry=0
      do while (icol.lt.nbcol)
         if(icol.eq.0) then
            if(imedfinal.eq.0) then  
               call givelifetimeepos(temptrans,1,timelim,ifailcol)
            else if(imedfinal.eq.1) then  
               call givelifetimeepos(temptrans,1,timelim,ifailcol)     
            else                                         
               write(8,*) 'unrecognized imedfinal:',imedfinal           
               close(8)
               stop
            endif
            if(ifailcol.eq.0) timelim=1.25*timelim
            write(*,*) 'hq7:lifetime = ', timelim,temptrans
            nbtausteps=int(timelim/timestep)+1
            if(iffirst) then
               nbtausteps0=nbtausteps              
               write(6,*) 'nbtausteps:',nbtausteps
               allocate(avnbclose(0:nbtausteps))
               do itime=0,nbtausteps
                  avnbclose(itime)=0.d0
               enddo
               iffirst=.false.
            endif
         endif

         call oneHQevent(icol+1,totcol+1,timelim,ifail)
         if(ifail.eq.0) then
            icol=icol+1
            totcol=totcol+1
            nbtry=0
         else
            if(ifail.eq.-1) then
               write(8,*) 'HQ array meson overshot'
               write(8,*) icol-1, nHQ, nd
               close(8)
               stop
            else
               nbtry=nbtry+1
               write(8,*) 'tentative event ',totcol+1,' skipped'
               if(nbtry.ge.10) then
                  write(8,*) 'skipping event 10x without success'
                  close(8)
                  stop
               endif
            endif
         endif
         deallocate(typquarkdyn,timelimdyn,hqmother,xdyn,iresdyn, 
     &    pdyn,tempatHQpos,imedatHQpos,ifsinglet,ifsingletfull,
     &    ifcorona) 
         deallocate(ipartcol,tempcol,timecol,dpcol)
         if(.not.(ifinterpolation)) deallocate(pcolafter)         
      enddo
      write(*,*)'total number of HQ for 1EPOS=', nHQ   
      write(*,*)'total number of HQ mesons for 1EPOS=', nd   
      call close_pack(totcol,weight)
      call saverand()
      if(ifsavehadro) call closesavehadro(totcol)  
      if(ifhqfromhq) call closeplasmascenario_of_b
      if((idist.eq.3).or.(idist.eq.4)) close(11) 
      end 

      subroutine exportdensnative(nbcol)
      use PSdens
      implicit none
      integer nbcol,it,ir,ip
 151  format(50(E12.5,' '))
      open(31,file="phasespacedensitynative.res")
      do it=1,10
         do ir=1,100
            write(31,151) (PSdensnative(it,ir,ip)/nbcol,ip=1,50)
         enddo
      enddo
      close(31)
      open(31,file="phasespacedensitynativehighpT.res")
      do it=1,10
         do ir=1,100
            write(31,151) (PSdensnativehighpT(it,ir,ip)/nbcol,ip=1,50)
         enddo
      enddo
      close(31)
      open(31,file="nbcloseavfromwigner.res")
      do it=1,10
         write(31,*) it,avnbclosewigner(it)/nbcol 
      enddo
      close(31)
      end

      subroutine exportdens(nbcol)
      use denysvar2
      implicit none
      integer nbcol,it,ir,ip
 151  format(50(E12.5,' '))
      open(31,file="phasespacedensity.res")
      do it=1,4
         do ir=1,50
            write(31,151) (phasespdens(it,ir,ip)/nbcol,ip=1,50)
         enddo
      enddo
      close(31)
      open(31,file="phasespacedensityhighpT.res")
      do it=1,4
         do ir=1,50
            write(31,151) (phasespdenshighpT(it,ir,ip)/nbcol,ip=1,50)
         enddo
      enddo
      close(31)
      open(31,file="averagesigma.res")
      do it=1,4
         if(averagesigma(it,1).gt.0.d0) then
            write(31,151) averagesigma(it,2)/averagesigma(it,1)
         else
            write(31,151) 0.
         endif
      enddo
      close(31)
      end

      subroutine exportdistribtemp(nbcol)
      use denysvar2
      implicit none
      integer itemp,nbcol
 151  format(50(E12.5,' '))
      open(31,file="tempdistribtherm.res")
      do itemp=1,100
         write(31,151) itemp*0.005-0.0025,distribtempattherm(itemp)/
     & nbcol
      enddo
      close(31)
      end


      subroutine exporttransfert
      use denysvar2
      implicit none
      integer ikick
      open(unit=10,file='distkick.res')
      do ikick=1,nperct
         write(10,*) (ikick-0.5)*deltap,ptransfspec(ikick)
      enddo
      close(10)      
      end





