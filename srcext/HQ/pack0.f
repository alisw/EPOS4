
      subroutine readin

      use PBGnumberevent
      use PBGsystem
      use PBGproductionmech
      use PBGinitnb
      use PBGmainpro1
      use PBGmainpro2
      use PBGgenvar
      use PBGsceplas
      use PBGspectra
      use PBGtransfer
      use PBGevolQ
      use PBGnumbercoll
      use PBGreductiondofs
      use forydistr
      implicit none
      character*500  fnin,fnus1,fnus2,fnus3,fnus4
      common/fname2/ fnin,fnus1,fnus2,fnus3,fnus4
      ifinput=.true.
      totcol=0
      weight=0
      ifinittimestep=.false.
      ifinitcTc=.false.
      ifinitXsectioncranck=.false.
      ifinitratemult=.false.
      ifinitc2=.false.
      ifinitfpmult=.false.
      ifinitalphasforrad=.false.
      ifinithighptcut=.false.
      ifinitnbcol=.false.
      ifinitnbcolsaveEPOS=.false.
      ifinitnbsubrun=.false.
      ifinittype_evol_Q=.false.
      ifinitimedfinal=.false.
      ifinitreddof=.false.
      ifinititypDproduct=.false.
      ifinititypfrag=.false.
      ifinititypcoal=.false.
      ifinittempcoal=.false.
      ifinitifcoalinflcell=.false.
      ifinitboltz_model=.false.
      ifinittypelpmdamp=.false.
      ifinitfpchoice=.false.
      ifinittunchoice=.false.
      ifinitifdecrease=.false.
      ifinitcoll=.false.
      ifinitifcronin=.false.
      ifinitifcol=.false.
      ifinitifrad=.false.
      ifinitiffullQCD=.false.
      ifinitiflpm=.false.
      ifinitifdamp=.false.
      ifinitifptcles=.false.
      ifinitifhighptselect=.false.
      ifinitifbselect=.false.
      ifinitish=.false.
      ifinitifhqfromhq=.false.
      ifinitifextspectra1=.false.
      ifinittypeextspectra1=.false.
      ifinitifshadowing1=.false.
      ifinittypeshadow1=.false.
      ifinitifextspectra2=.false.
      ifinittypeextspectra2=.false.
      ifinitifshadowing2=.false.
      ifinittypeshadow2=.false.
      ifinithqskip=.false.
      ifinit2bodySC=.false.
      ifinitsingletred=.false.
      ifinitdpmaxabscorona=.false.
      ifinittypedistrel=.false.
      ifinitparamdistrelc=.false.
      ifinitparamdistrelb=.false.
      open(515,file=fnus1(1:index(fnus1,' ')-1)//'input7.dat'
     .  ,status='old')

        tauf=0d0
        if(ifinput) then
         read(515,*) ISCE
        else
         ISCE=7
        endif
 
      if(ifinput) then
      read(515,*) b_inf, b_sup
      else
         b_inf=0.
         b_sup=14.
      endif
      if(isce.ne.7)then
        if(ifinput) then
           read(515,*) A
           read(515,*) sqrtS
        else
           A=208
           sqrtS=5500d0
        endif
      else
         if(ifinput) then
           read(515,*) A
           read(515,*) sqrtS
        endif
      endif
      if(ifinput) then
         read(515,*) nbcol,nbsubrun,nbcolsaveEPOS
         read(515,*) ifallowhqskip
      else
         nbcol=1000
         nbcolsaveEPOS=50
         nbsubrun=1
         ifallowhqskip=.true.
      endif
      ifinithqskip=.true.
      ifinitnbsubrun=.true.
      if(ifinput) then
         read(515,*) type_evol_Q
         read(515,*) timestep 
      else
         type_evol_Q=2
         timestep=0.1d0
      endif
      ifinittype_evol_Q=.true.
      ifinittimestep=.true.      
      itypeorder=2
      if(ifinput) then
         read(515,*) imedfinal
         read(515,*) ifdecrease
      else
         imedfinal=1
         ifdecrease=.true.
      endif
      ifinitimedfinal=.true.
      ifinitifdecrease=.true.
      idistspatial=1
      if(ifinput) then
         read(515,*) ifextspectra(1),typeextspectra(1)
         read(515,*) ifshadowing(1),typeshadow(1)
         read(515,*) ifextspectra(2),typeextspectra(2)
         read(515,*) ifshadowing(2),typeshadow(2)
      else
         ifextspectra(1)=.false.
         typeextspectra(1)=0
         ifshadowing(1)=.false.
         typeshadow(1)=0
         ifextspectra(2)=.false.
         typeextspectra(1)=0
         ifshadowing(2)=.false.
         typeshadow(2)=0
      endif
      ifinitifhqfromhq=.true.
      ifinitifextspectra1=.true.
      ifinittypeextspectra1=.true.
      ifinitifshadowing1=.true.
      ifinittypeshadow1=.true.
      ifinitifextspectra2=.true.
      ifinittypeextspectra2=.true.
      ifinitifshadowing2=.true.
      ifinittypeshadow2=.true.
      if(ifinput) then
         read(515,*) ifyat0
      else
         ifyat0=.false.
      endif
      if(ifinput) then
         read(515,*) ifcronin
      else
         ifcronin=.true.
      endif
      ifinitifcronin=.true.
      infotrack=0
      tolerance=0

      if(ifinput) then
         read(515,*) hhbarprod
         read(515,*) reddof,cTc
      else
         hhbarprod=1
         reddof=0
         cTc=1.d0
      endif
      ifinitcTc=.true.      
      ifinitreddof=.true.
      end subroutine readin



 



    
     
      integer function plasmascen()
      use PBGsceplas
      implicit none
      plasmascen=ISCE
      end function plasmascen
   
 
      subroutine big_initialization

      use PBGmainpro1
      use PBGmainpro2
      use PBGsystem
      use PBGfragandcoal
      use PBGgenvar
      use PBGsceplas
      use PBGevolQ
      use PBGimpparevent
      use PBGnbimpapara
      use PBGXsection
      implicit none
      character*500  fnin,fnus1,fnus2,fnus3,fnus4
      common/fname2/ fnin,fnus1,fnus2,fnus3,fnus4
 
      call resonanceinit
      call inicoalandfrag
      if(type_evol_Q.eq.1) then
      else if(type_evol_Q.eq.2) then
      else
         stop    
      endif    
      call initred_thmass
      call initrates
      call initfokker
      call initratesradiatlog
      if(ifdecrease) then
      else
      endif
      call init2body
      call initdecay
      call initccbardistrib
      call initplasmascenario_glob
      call inithqproduct
      if(ifcronin)then
      else
      endif
      call initpacking(nbcol,timestep)
      close(515)
      if(isce.ne.7)then
        nbimpact=14
        call classify_impact(b_inf,b_sup,nbcol,nbsubrun,nbimpact
     &                      ,nbevent)
      endif
      call initonerate_HQ
      return
      end subroutine big_initialization


      subroutine big_initialization2
      call inicoalandfrag2
      call initrates2
      call initfokker2
      call initratesradiatlog2
      end subroutine big_initialization2


      subroutine initrand(recover)
      use PBGran2junk
      implicit none
      logical recover
      integer i
      if(.not.(recover)) then
         idum=-987654321
         if (idum.gt.0) then
            stop
         else    
         endif
         IDUM2=123456789
         IY=0
         do i=1,ntab
            IV(i)=0
         enddo
      else
         open(9,file='idum_init.var',status='old')
         read(9,*) idum,idum2
         read(9,*) IY
         do i=1,ntab
            read(9,*) IV(i)
         enddo
         close(9)
      endif
      return
      end

      subroutine saverand()
      use PBGran2junk
      implicit none
      integer i
      open(9,file='idum_init.var',status='replace')
      write(9,*) idum,idum2
      write(9,*) IY
      do i=1,ntab
         write(9,*) IV(i)
      enddo
      close(9)
      return
      end


      subroutine initpaw
      implicit none
      return
      end


      subroutine resonanceinit
      use PBGfragandcoal
      use PBGgenvar
      use PBGsceplas
      use PBGpartinmedium
      use PBGXsection      
      use PBGtransfer
      implicit none
      integer ires,ires2,alpha,beta
      double precision maxsigmafus,pi

      parameter (pi=3.14159265358979323844d0)
      do alpha=0,3
         nlor(alpha)=0
         do beta=0,3
            g(alpha,beta)=0
         enddo
         g(alpha,alpha)=-1
      enddo
      g(0,0)=1
      nlor(0)=1
      MCCHOICE=1 
      if(mcchoice.eq.1) then
         mc=1.5
         eps0c=0.603
      else
         if(mcchoice.eq.2) then
            mc=1.94
            eps0c=0.78
         else
            stop
         endif
      endif
      mc2=mc*mc
      muc=mc/2.
      a0c=hbarc/sqrt(2*eps0c*muc)
      mpsi=2*mc-eps0c
      mD=1.87
      mDs=1.968
      mlambdac=2.286
      mxic=2.469
      momegac=2.695
      mb=5.1
      eps0b=0.75
      mb2=mb*mb
      mub=mb/2.
      a0b=hbarc/sqrt(2*eps0b*mub)
      mups=2*mb-eps0b
      mBmes=5.279
      mBs=5.367
      mlambdab=5.619
      mxib=5.797
      momegab=6.046

      resinfo_ifquark(1)=.true.
      itypquark(1)=1
      resinfo_mass(1)=mc
      resinfo_ccharge(1)=1
      resinfo_spin(1)=1
      resinfo_degen(1)=2*3
      resinfo_dissoc(1)=0
      resinfo_antipart(1)=2
      resinfo_ifquark(2)=.true.
      itypquark(2)=1
      resinfo_mass(2)=mc
      resinfo_ccharge(2)=-1
      resinfo_spin(2)=1
      resinfo_degen(2)=2*3
      resinfo_dissoc(2)=0
      resinfo_antipart(2)=1
      resinfo_ifquark(3)=.false.
      resinfo_mass(3)=mpsi
      resinfo_ccharge(3)=0
      resinfo_spin(3)=2
      resinfo_degen(3)=3
      resinfo_dissoc(3)=1     
      resinfo_antipart(3)=-1
      resinfo_ifquark(4)=.false.
      resinfo_mass(4)=mD
      resinfo_ccharge(4)=0
      resinfo_spin(4)=0
      resinfo_degen(4)=2
      resinfo_dissoc(4)=0
      resinfo_antipart(4)=5
      resinfo_ifquark(5)=.false.
      resinfo_mass(5)=mD
      resinfo_ccharge(5)=0
      resinfo_spin(5)=0
      resinfo_degen(5)=2
      resinfo_dissoc(5)=0
      resinfo_antipart(5)=4

      resinfo_ifquark(6)=.false.
      resinfo_mass(6)=mDs
      resinfo_ccharge(6)=1
      resinfo_spin(6)=0
      resinfo_degen(6)=1
      resinfo_dissoc(6)=0
      resinfo_antipart(6)=7
      resinfo_ifquark(7)=.false.
      resinfo_mass(7)=mDs
      resinfo_ccharge(7)=-1
      resinfo_spin(7)=0
      resinfo_degen(7)=1
      resinfo_dissoc(7)=0
      resinfo_antipart(7)=6
      resinfo_ifquark(8)=.false.
      resinfo_mass(8)=mlambdac
      resinfo_ccharge(8)=1
      resinfo_spin(8)=0.5
      resinfo_degen(8)=1
      resinfo_dissoc(8)=0
      resinfo_antipart(8)=9
      resinfo_ifquark(9)=.false.
      resinfo_mass(9)=mlambdac
      resinfo_ccharge(9)=-1
      resinfo_spin(9)=0.5
      resinfo_degen(9)=1
      resinfo_dissoc(9)=0
      resinfo_antipart(9)=8
      resinfo_ifquark(10)=.false.
      resinfo_mass(10)=mxic
      resinfo_ccharge(10)=1
      resinfo_spin(10)=0.5
      resinfo_degen(10)=1
      resinfo_dissoc(10)=0
      resinfo_antipart(10)=11
      resinfo_ifquark(11)=.false.
      resinfo_mass(11)=mxic
      resinfo_ccharge(11)=-1
      resinfo_spin(11)=0.5
      resinfo_degen(11)=1
      resinfo_dissoc(11)=0
      resinfo_antipart(11)=10
      resinfo_ifquark(12)=.false.
      resinfo_mass(12)=momegac
      resinfo_ccharge(12)=1
      resinfo_spin(12)=0.5
      resinfo_degen(12)=1
      resinfo_dissoc(12)=0
      resinfo_antipart(12)=13
      resinfo_ifquark(13)=.false.
      resinfo_mass(13)=momegac
      resinfo_ccharge(13)=-1
      resinfo_spin(13)=0.5
      resinfo_degen(13)=1
      resinfo_dissoc(13)=0
      resinfo_antipart(13)=12

      resinfo_ifquark(nresc+1)=.true.
      itypquark(nresc+1)=2
      resinfo_mass(nresc+1)=mb
      resinfo_ccharge(nresc+1)=1
      resinfo_spin(nresc+1)=1
      resinfo_degen(nresc+1)=2*3
      resinfo_dissoc(nresc+1)=0
      resinfo_antipart(nresc+1)=nresc+2
      resinfo_ifquark(nresc+2)=.true.
      itypquark(nresc+2)=2
      resinfo_mass(nresc+2)=mb
      resinfo_ccharge(nresc+2)=-1
      resinfo_spin(nresc+2)=1
      resinfo_degen(nresc+2)=2*3
      resinfo_dissoc(nresc+2)=0
      resinfo_antipart(nresc+2)=nresc+1
      resinfo_ifquark(nresc+3)=.false.
      resinfo_mass(nresc+3)=mups
      resinfo_ccharge(nresc+3)=0
      resinfo_spin(nresc+3)=2
      resinfo_degen(nresc+3)=3
      resinfo_dissoc(nresc+3)=1     
      resinfo_antipart(nresc+3)=-1
      resinfo_ifquark(nresc+4)=.false.
      resinfo_mass(nresc+4)=mBmes
      resinfo_ccharge(nresc+4)=1
      resinfo_spin(nresc+4)=0
      resinfo_degen(nresc+4)=2
      resinfo_dissoc(nresc+4)=0
      resinfo_antipart(nresc+4)=nresc+5
      resinfo_ifquark(nresc+5)=.false.
      resinfo_mass(nresc+5)=mBmes
      resinfo_ccharge(nresc+5)=-1
      resinfo_spin(nresc+5)=0
      resinfo_degen(nresc+5)=2
      resinfo_dissoc(nresc+5)=0
      resinfo_antipart(nresc+5)=nresc+4

      resinfo_ifquark(nresc+6)=.false.
      resinfo_mass(nresc+6)=mBs
      resinfo_ccharge(nresc+6)=-1
      resinfo_spin(nresc+6)=0
      resinfo_degen(nresc+6)=1
      resinfo_dissoc(nresc+6)=0
      resinfo_antipart(nresc+6)=nresc+7
      resinfo_ifquark(nresc+7)=.false.
      resinfo_mass(nresc+7)=mBs
      resinfo_ccharge(nresc+7)=1
      resinfo_spin(nresc+7)=0
      resinfo_degen(nresc+7)=1
      resinfo_dissoc(nresc+7)=0
      resinfo_antipart(nresc+7)=nresc+6
      resinfo_ifquark(nresc+8)=.false.
      resinfo_mass(nresc+8)=mlambdab
      resinfo_ccharge(nresc+8)=-1
      resinfo_spin(nresc+8)=0.5
      resinfo_degen(nresc+8)=1
      resinfo_dissoc(nresc+8)=0
      resinfo_antipart(nresc+8)=nresc+9
      resinfo_ifquark(nresc+9)=.false.
      resinfo_mass(nresc+9)=mlambdab
      resinfo_ccharge(nresc+9)=1
      resinfo_spin(nresc+9)=0.5
      resinfo_degen(nresc+9)=1
      resinfo_dissoc(nresc+9)=0
      resinfo_antipart(nresc+9)=nresc+8
      resinfo_ifquark(nresc+10)=.false.
      resinfo_mass(nresc+10)=mxib
      resinfo_ccharge(nresc+10)=-1
      resinfo_spin(nresc+10)=0.5
      resinfo_degen(nresc+10)=1
      resinfo_dissoc(nresc+10)=0
      resinfo_antipart(nresc+10)=nresc+11
      resinfo_ifquark(nresc+11)=.false.
      resinfo_mass(nresc+11)=mxib
      resinfo_ccharge(nresc+11)=1
      resinfo_spin(nresc+11)=0.5
      resinfo_degen(nresc+11)=1
      resinfo_dissoc(nresc+11)=0
      resinfo_antipart(nresc+11)=nresc+10
      resinfo_ifquark(nresc+12)=.false.
      resinfo_mass(nresc+12)=momegab
      resinfo_ccharge(nresc+12)=-1
      resinfo_spin(nresc+12)=0.5
      resinfo_degen(nresc+12)=1
      resinfo_dissoc(nresc+12)=0
      resinfo_antipart(nresc+12)=nresc+13
      resinfo_ifquark(nresc+13)=.false.
      resinfo_mass(nresc+13)=momegab
      resinfo_ccharge(nresc+13)=1
      resinfo_spin(nresc+13)=0.5
      resinfo_degen(nresc+13)=1
      resinfo_dissoc(nresc+13)=0
      resinfo_antipart(nresc+13)=nresc+12

      do ires=1,nres
         do ires2=1,nres
            ifreac(ires,ires2)=0
         enddo
      enddo
      ifreac(1,2)=1
      ifreac(2,1)=1
      thresreac(1,2)=2*mc
      thresreac(2,1)=2*mc
      if(ifinput) then
         read(515,*) XSECTIONCRANCK
      else
         XSECTIONCRANCK=0.
      endif
      ifinitXsectioncranck=.true.      
      biggestb2(1,2)=maxsigmafus(1,2,thresreac(1,2))/pi
      biggestb2(2,1)=biggestb2(1,2)
      return
      end


      subroutine initccbardistrib
      use PBGsystem
      use PBGsceplas
      use PBGpartinmedium
      use PBGifinputfile
      use PBGdistccbarbis
      use PBGtransfer
      implicit none
      character*50 inputdistname
      if(ifinput) then
         read(515,*) idist
         read(515,*) typedistrel
         read(515,*) paramdistrel(1),paramdistrel(2)
      else
         IDIST=1
         typedistrel=2
         paramdistrel(1)=0.25
         paramdistrel(2)=0.08
      endif
      ifinittypedistrel=.true.
      ifinitparamdistrelc=.true.
      ifinitparamdistrelb=.true.
      NBCCBARINDATABASE=100000
      if(idist.eq.3) then
         if(sqrts.eq.200) then
          INPUTDISTNAME='../data.dir/init_ccbar_RHIC_1E5.dat'
         else if((abs(sqrts-5500).lt.0.1)
     &          .or.(abs(sqrtS-2760).lt.0.1)) then
           INPUTDISTNAME='../data.dir/init_ccbar_LHC_1E5.dat'
         else
            stop
         endif
         open(11,file=inputdistname,access='direct',
     &        status='old',form='formatted',recl=70)
         write(6,*) 'c-cbar will be taken in a database of ',
     &    nbccbarindatabase,' events'
      endif
      return
      end


      subroutine inicoalandfrag
      use PBGforprobcoal
      use PBGforprobcoalTL
      use PBGfragandcoal
      use PBGgenvar
      use PBGtransfer
      implicit none
      character*500  fnin,fnus1,fnus2,fnus3,fnus4
      common/fname2/ fnin,fnus1,fnus2,fnus3,fnus4
      character*45 probcoalnamec(2),probcoalnameb
      character*45 probcoalnamebis(2)
      data probcoalnamec/'../data.dir/prob_coal_1500.dat',
     &  '../data.dir/prob_coal_1940.dat'/
      data probcoalnameb/'../data.dir/prob_coal_5100.dat'/
      data probcoalnamebis/
     &'../data.dir/prob_coal_setII_PB_mq100_1500.dat',
     & '../data.dir/prob_coal_setII_PB_mq100_5100.dat'/

      character*50 probcoalnamebisTL(2)
      data probcoalnamebisTL
     & /'../data.dir/prob_coal_setII_PB_mq100_TL_1500.dat',
     &  '../data.dir/prob_coal_setII_PB_mq100_TL_5100.dat'/
      if(ifinput) then
         read(515,*) ITYPDPRODUCT
      else
         ITYPDPRODUCT=3
      endif
      pcrithadr(1)=(resinfo_mass(4)**2-mc2)/(2*resinfo_mass(4))
      pcrithadr(2)=(resinfo_mass(nresc+4)**2-mb2)
     &            /(2*resinfo_mass(nresc+4))
      if(ifinput) then
         read(515,*) ITYPFRAG
      else
         ITYPFRAG=3
      endif
      if(ifinput) then
         read(515,*) ITYPCOAL,ifcoalinflcell,tempcoal
         read(515,*) ifsavehadro

      else
         ITYPCOAL=3
         ifcoalinflcell=.true.
         tempcoal=0.154
         ifsavehadro=.false.
      endif
      ifinititypcoal=.true.
      ifinitifcoalinflcell=.true.
      ifinittempcoal=.true.
      return
      end

      subroutine inicoalandfrag2
      use PBGforprobcoal
      use PBGforprobcoalTL
      use PBGfragandcoal
      use PBGgenvar
      character*500  fnin,fnus1,fnus2,fnus3,fnus4
      common/fname2/ fnin,fnus1,fnus2,fnus3,fnus4
      integer ip,itype,itempcoal,lentrue
      character*45 probcoalnamec(2),probcoalnameb
      character*55 probcoalnamebis(2),probcoalnameter(2,2),
     & probcoalnamequad(2,2)
      character*45 probcoalnewcb(2)
      character*65 probcoalnewcH(5),probcoalnewbH(5)
      data probcoalnewcb/'COALFRAG/probJX2/0103/call_010315_154.dat',
     & 'COALFRAG/probJX2/0103/call_010345_154.dat'/
      data probcoalnewcH/
     & 'COALFRAG/probJX2/0103/D0_010315_154.dat',
     & 'COALFRAG/probJX2/0103/Ds_010315_154.dat',
     & 'COALFRAG/probJX2/0103/Lambdac_010315_154.dat',
     & 'COALFRAG/probJX2/0103/Xic_010315_154.dat',
     & 'COALFRAG/probJX2/0103/Omegac_010315_154.dat'/
      data probcoalnewbH/
     & 'COALFRAG/probJX2/0103/B0_010345_154.dat',
     & 'COALFRAG/probJX2/0103/Bs_010345_154.dat',
     & 'COALFRAG/probJX2/0103/Lambdab_010345_154.dat',
     & 'COALFRAG/probJX2/0103/Xib_010345_154.dat',
     & 'COALFRAG/probJX2/0103/Omegab_010345_154.dat'/
      data probcoalnamec/'COALFRAG/prob_coal_1500.dat',
     &  'COALFRAG/prob_coal_1940.dat'/
      data probcoalnameb/'COALFRAG/prob_coal_5100.dat'/
      data probcoalnamebis/
     & 'COALFRAG/prob_coal_setII_PB_mq100_1500.dat',
     & 'COALFRAG/prob_coal_setII_PB_mq100_5100.dat'/
      data probcoalnameter/
     & 'COALFRAG/prob_coal_setII_PB_mq100_T0155_1500.dat',
     & 'COALFRAG/prob_coal_setII_PB_mq100_T0155_5100.dat',
     & 'COALFRAG/prob_coal_setII_PB_mq100_1500.dat',
     & 'COALFRAG/prob_coal_setII_PB_mq100_5100.dat'/
      data probcoalnamequad/
     & 'COALFRAG/prob_coal_set2017_mq292_T0155_1500.dat',
     & 'COALFRAG/prob_coal_set2017_mq292_T0155_5100.dat',
     & 'COALFRAG/prob_coal_set2017_mq267_1500.dat',
     & 'COALFRAG/prob_coal_set2017_mq267_5100.dat'/
      integer id,iyds
      double precision dexpect,dread
      character*60 probcoalnamebisTL(2),probcoalnameterTL(2,2),
     & probcoalnamequadTL(2,2)
      data probcoalnamebisTL/
     & 'COALFRAG/prob_coal_setII_PB_mq100_TL_1500_2009.dat',
     &  'COALFRAG/prob_coal_setII_PB_mq100_TL_5100_2009.dat'/
      data probcoalnameterTL/
     & 'COALFRAG/prob_coal_setII_PB_mq100_T0155_TL_1500.dat',
     & 'COALFRAG/prob_coal_setII_PB_mq100_T0155_TL_5100.dat',
     & 'COALFRAG/prob_coal_setII_PB_mq100_TL_1500_2017.dat',
     & 'COALFRAG/prob_coal_setII_PB_mq100_TL_5100_2017.dat'/
      data probcoalnamequadTL/
     & 'COALFRAG/prob_coal_set2017_mq292_T0155_TL_1500.dat',
     & 'COALFRAG/prob_coal_set2017_mq292_T0155_TL_5100.dat',
     & 'COALFRAG/prob_coal_set2017_mq267_TL_1500.dat',
     & 'COALFRAG/prob_coal_set2017_mq267_TL_5100.dat'/
      integer ifmtx
      call getMonitorFileIndex(ifmtx)
      if(itypfrag.eq.2) then
         peterepsq(1)=0.15
         peterepsq(2)=0.0135
         peterepsh(1)=1.55418
         peterepsh(2)=1.11514
         peterzmax(1)=0.715241
         peterzmax(2)=0.893374
         petermaxval(1)=0.483951
         petermaxval(2)=7.96661
      else if(itypfrag.eq.3) then
         write(6,*) 'initiating Braaten and co'
         call initrepartBraat
         call initrepartKart
      endif
      if(abs(tempcoal-0.154).le.0.001) then
         itempcoal=1
      else if(abs(tempcoal-0.165).le.0.001) then
         itempcoal=2
      else
         write(ifmtx,*) 'tempcoal not recognized'
         close(ifmtx)
         stop
      endif
      if((ITYPCOAL.eq.1).and.(.not.(ifcoalinflcell))) then
         write(ifmtx,*) 'ifcoalinflcell should be T for ITYPCOAL=1'
         close(ifmtx)
         stop
      endif
      if(itypcoal.eq.1) then
         write(ifmtx,*) 'itypcoal=1 notimplemented in EPOSHQ'
         close(ifmtx)
         stop
      else
         if((itypcoal.eq.2).or.(itypcoal.eq.3)) then
            mqconstit=0.1
            if(itempcoal.eq.1) then
               alphadov=0.83375
            else if(itempcoal.eq.2) then
               alphadov=0.876
            else
               write(ifmtx,*) 'unforseseen itempcoal'
               close(ifmtx)
               stop
            endif
         else if(itypcoal.eq.4) then
            if(itempcoal.eq.1) then
               mqconstit=0.292
               alphadov=0.51
            else if(itempcoal.eq.2) then
               mqconstit=0.267
               alphadov=0.627
            else
               write(ifmtx,*) 'unforseseen itempcoal'
               close(ifmtx)
               stop
            endif
         else if(itypcoal.eq.5) then
          mqconstit=0.1
          msconstit=0.3 
         else
            write(ifmtx,*) 'unrecognized type of coal'
            close(ifmtx)
            stop
         endif
         do itype=1,2
            if(itypcoal.eq.2) then
               open(141,file=fnus1(1:index(fnus1,' ')-1)
     &  //probcoalnamebis(itype)(1:lentrue(probcoalnamebis(itype),55))
     &  ,status='old')      
               write(*,*) 'Coalescence prob: importing from file ',
     &              probcoalnamebis(itype)
            elseif(itypcoal.eq.3) then
               open(141,file=fnus1(1:index(fnus1,' ')-1)
     &  //probcoalnameter(itype,itempcoal)(1:lentrue(probcoalnameter(
     &     itype,itempcoal),55)),status='old')
               write(*,*) 'Coalescence prob: importing from file ',
     &              probcoalnameter(itype,itempcoal)
            elseif(itypcoal.eq.4) then
               open(141,file=fnus1(1:index(fnus1,' ')-1)
     &  //probcoalnamequad(itype,itempcoal)(1:lentrue(probcoalnamequad(
     &     itype,itempcoal),55)),status='old')
               write(*,*) 'Coalescence prob: importing from file ',
     &        probcoalnamequad(itype,itempcoal)
            elseif(itypcoal.eq.5) then
               open(141,file=fnus1(1:index(fnus1,' ')-1)
     &  //probcoalnewcb(itype)(1:lentrue(probcoalnewcb(itype),55))
     &  ,status='old')               
              write(*,*) 'Coalescence prob: importing from file ',
     &        probcoalnewcb(itype)
            endif
            read(141,*) pmax(itype)
            read(141,*) nbp(itype)
            if(nbp(itype).gt.nbpmax) then
               write(ifmtx,*) 'error in initcoal: nbp(itype) > nbpmax'
               close(ifmtx)
               close(141)
               stop
            endif
            deltp(itype)=pmax(itype)/nbp(itype)
            do ip=0,nbp(itype)
               read(141,*) probcoal(itype,ip)
            enddo
            close(141)
            if(itypcoal.eq.2) then
               open(141,file=fnus1(1:index(fnus1,' ')-1)
     &  //probcoalnamebisTL(itype)(1:lentrue(probcoalnamebisTL(itype),
     &   60)),status='old')      
               write(*,*) 'Coalescence prob: importing from file ',
     &        probcoalnamebisTL(itype)
            elseif(itypcoal.eq.3) then
               open(141,file=fnus1(1:index(fnus1,' ')-1)
     &  //probcoalnameterTL(itype,itempcoal)(1:lentrue(probcoalnameterTL
     &    (itype,itempcoal),60)),status='old')
               write(*,*) 'Coalescence prob: importing from file ',
     &        probcoalnameterTL(itype,itempcoal)
            else
               open(141,file=fnus1(1:index(fnus1,' ')-1)
     &  //probcoalnamequadTL(itype,itempcoal)(1:lentrue(
     &   probcoalnamequadTL(itype,itempcoal),60)),status='old')
               write(*,*) 'Coalescence prob: importing from file ',
     &        probcoalnamequadTL(itype,itempcoal)
            endif
            read(141,*) dcoalmin(itype),dcoalmax(itype),nbd(itype)
            read(141,*) ydsmin(itype),ydsmax(itype),nbyds(itype)
            if((nbd(itype).gt.nbdmax).or.(nbyds(itype).gt.nbydsmax))then
               write(ifmtx,*) 'error in initcoal: nbd/yds(itype) > max'
               close(ifmtx)
               close(141)
               stop
            endif
            deltdcoal(itype)=(dcoalmax(itype)-dcoalmin(itype))/
     &           nbd(itype)
            dyds(itype)=(ydsmax(itype)-ydsmin(itype))/nbyds(itype)
            do id=0,nbd(itype)
               dexpect=dcoalmin(itype)+id*deltdcoal(itype)
               read(141,*) dread
               if(abs(dread-dexpect).gt.1e-5) then
                  write(ifmtx,*) 'mismatch d in initcoalandfrag'
                  write(ifmtx,*) 'itype quark:',itype
                  write(ifmtx,*) 'expected d:',dexpect
                  write(ifmtx,*) 'read d:',dread
                  close(ifmtx)
                  stop
               endif
               do iyds=0,nbyds(itype)
                  read(141,*) probcoalTL(itype,id,iyds)
               enddo
            enddo
            close(141)
         enddo
      endif
      if(itypcoal.eq.5) then
          do itypHade=1,5
           open(141,file=fnus1(1:index(fnus1,' ')-1)
     &  //probcoalnewcH(itypHade)(1:lentrue(probcoalnewcH
     &  (itypHade),65)),status='old')
            write(*,*) 'Coalescence prob: importing from file ',
     &      probcoalnewcH(itypHade)
            read(141,*) pmax(itypHade)
            read(141,*) nbp(itypHade)
            if(nbp(itypHade).gt.nbpmaxH) then
               write(8,*) 'error in initcoal: nbp(itype) > nbpmax'
               close(8)
               close(141)
               stop
            endif
            deltp(itypHade)=pmax(itypHade)/nbp(itypHade)
            do ip=0,nbp(itypHade)
               read(141,*) probcoalHad(1,itypHade,ip)
            enddo
            close(141)
          enddo
          do itypHade=1,5
           open(141,file=fnus1(1:index(fnus1,' ')-1)
     &  //probcoalnewbH(itypHade)(1:lentrue(probcoalnewbH
     &  (itypHade),65)),status='old')
            write(*,*) 'Coalescence prob: importing from file ',
     &      probcoalnewbH(itypHade)
            read(141,*) pmax(itypHade)
            read(141,*) nbp(itypHade)
            if(nbp(itypHade).gt.nbpmaxH) then
               write(8,*) 'error in initcoal: nbp(itype) > nbpmax'
               close(8)
               close(141)
               stop
            endif
            deltp(itypHade)=pmax(itypHade)/nbp(itypHade)
            do ip=0,nbp(itypHade)
               read(141,*) probcoalHad(2,itypHade,ip)
            enddo
            close(141)
          enddo
       endif
       resinfo_mass(0)=mqconstit
       mqconstit2=mqconstit**2
      if(ifsavehadro) call initsavehadro
      return
      end



      subroutine initplasmascenario_glob


      use PBGsystem
      use PBGsceplas
      use PBGhydro
      use PBGpartinmedium
      use PBGnres
      use PBGifinputfile
      implicit none
      double precision edensoft,
     & dedy,densenmax,dedyrhic
      parameter(dedyrhic=1280)



 
      if(isce.le.4) then
         temptrans=0.17         
         if(isce.eq.1) then
            tau0=0.6
            ri0=6.5
            li0=1.
            vl0=1.
            vr0=0.58
            pbar=0
         endif   
         if(isce.eq.2) then
            tau0=0.6
            ri0=6.5
            li0=1.
            vl0=1.
            vr0=0
            pbar=0
         endif   
         if(isce.eq.3) then
            ri0=6.5
            alphaprof=0.5
            li0=0
            vl0=1.
            vr0=0
            pbar=0
            if(abs(sqrts-200).lt.0.1) then
               tau0=0.6
               tempinit=0.36
               ifaxial=.false.
            else if((abs(sqrts-5500).lt.0.1)
     &             .or.(abs(sqrtS-2760).lt.0.1)) then
               tau0=0.2
               tempinit=0.5
               ifaxial=.false.
            else
               stop
            endif
            densenmax=edensoft(tempinit)
            dedy=3.14159*ri0**2*densenmax*tau0/(1+alphaprof)
         endif   
         if(isce.eq.4) then
            tau0=0.6
            ri0=6.5
            li0=0
            vl0=1.
            vr0=0.58
            pbar=0
         endif   
         entrans_max=edensoft(1.001*temptrans)
         entrans_min=edensoft(0.999*temptrans)         
      else if(isce.eq.5) then
            ri0=6.5d0
            li0=0d0
            if(abs(sqrts-200).lt.0.1) then
               tau0=0.6d0
            else if(abs(sqrts-2760).lt.0.1) then
               tau0=0.6d0
               sdens=195d0
            else
               stop
            endif
            vl0=1d0
            temptrans=0.155d0
      else if(isce.eq.6) then
               if(abs(sqrts-200).gt.0.1) then
                  stop
               endif
               call inithydropasi
               ri0=6.5
               li0=0
               tau0=0.6d0
               vl0=1d0
               entrans_max=1.65
               entrans_min=0.55
               temptrans=0.165

      else if(isce.eq.7) then
         vl0=1d0
         temptrans=0.154d0
      endif         
       write(*,*) 'temptrans= ', temptrans
      tempdissoc(1)=1.D10
      tempdissoc(2)=1.D10
      if(ifinput) then
         read(515,*) TEMPDISSOC(3)
      else
         TEMPDISSOC(3)=1.D10
      endif
      write(*,*) 'dissoc temp of J/Psi=',tempdissoc(3)
      tempdissoc(4)=temptrans+0.001
      tempdissoc(5)=temptrans+0.001
      tempdissoc(6)=temptrans+0.001
      tempdissoc(7)=temptrans+0.001
      tempdissoc(8)=temptrans+0.001
      tempdissoc(9)=temptrans+0.001
      tempdissoc(10)=temptrans+0.001
      tempdissoc(11)=temptrans+0.001
      tempdissoc(12)=temptrans+0.001
      tempdissoc(13)=temptrans+0.001
      tempdissoc(nresc+1)=1.D10
      tempdissoc(nresc+2)=1.D10
      if(ifinput) then
         read(515,*) TEMPDISSOC(nresc+3)
      else
         TEMPDISSOC(nresc+3)=1.D10
      endif
      tempdissoc(nresc+4)=temptrans+0.001
      tempdissoc(nresc+5)=temptrans+0.001
      tempdissoc(nresc+6)=temptrans+0.001
      tempdissoc(nresc+7)=temptrans+0.001
      tempdissoc(nresc+8)=temptrans+0.001
      tempdissoc(nresc+9)=temptrans+0.001
      tempdissoc(nresc+10)=temptrans+0.001
      tempdissoc(nresc+11)=temptrans+0.001
      tempdissoc(nresc+12)=temptrans+0.001
      tempdissoc(nresc+13)=temptrans+0.001
      return
      write(*,*) 'stop of sr'
      end


      subroutine inithqproduct
      use PBGifinputfile
      use PBGinitnbmult
      implicit none
      if(ifinput) then
         read(515,*) avnc0ncollmult,avnpsi0ncollmult,avnb0ncollmult,
     &        avnups0ncollmult
      else
         avnc0ncollmult=1.
         avnpsi0ncollmult=0.
         avnb0ncollmult=1.
         avnups0ncollmult=0.
      endif
      end


      subroutine initplasmascenario_of_b
      use PBGinitnb
      use PBGsigmahqprod
      use PBGvariancesig
      use PBGsceplas
      use PBGhydro
      use PBGspectra
      use PBGpartinmedium
      use PBGinitnbmult
      use PBGalltaa
      use PBGallvariance
       implicit none
      integer ib
      double precision sigmabid
      write(6,*) '================================================'
      write(6,*) 'entering initplasmascenario_of_b with b=',b
      write(6,*) '================================================'
      nbofb=14.
      ib=nint(b)
      if(ib.gt.nbofb) then
      endif
      if(ifextspectra(1)) then
         call initptspectrum_of_b(ib,1,typeextspectra(1),sigmabid,
     &     ifshadowing(1),typeshadow(1))
         if(sigmabid.gt.0.d0) sigprodc=sigmabid
      endif
      if(ifextspectra(2)) then
         call initptspectrum_of_b(ib,2,typeextspectra(2),sigmabid,
     &     ifshadowing(2),typeshadow(2))
         if(sigmabid.gt.0.d0) sigprodb=sigmabid
      endif

      avnc0ncoll=alltaa(ib)*sigprodc
      avnb0ncoll=alltaa(ib)*sigprodb 
      avnpsi0ncoll=alltaa(ib)*sigprodpsi
      avnups0ncoll=alltaa(ib)*sigprodups
      if(avnc0ncollmult.lt.0) then
         avnc0ncoll=-avnc0ncollmult
      else
         avnc0ncoll=avnc0ncollmult*avnc0ncoll
      endif
      if(avnb0ncollmult.lt.0) then
         avnb0ncoll=-avnb0ncollmult
      else
         avnb0ncoll=avnb0ncollmult*avnb0ncoll
      endif
      if(avnpsi0ncollmult.lt.0) then
         avnpsi0ncoll=-avnpsi0ncollmult
      else
         avnpsi0ncoll=avnpsi0ncollmult*avnpsi0ncoll
      endif
      if(avnups0ncollmult.lt.0) then
         avnups0ncoll=-avnups0ncollmult
      else
         avnups0ncoll=avnups0ncollmult*avnups0ncoll
      endif
      sigma=allsigma(ib)
      sigmax=allsigmax(ib)
      sigmay=allsigmay(ib)
      if(ISCE.eq.5) then
         call inithydrokolb_bis

      else if(ISCE.eq.7) then
        write(*,*) 'call initEpos'
      else if(b.gt.0.) then
            stop
      endif
      return
      end

      subroutine closeplasmascenario_of_b
      use PBGsceplas
      use PBGhydro
      implicit none
      write(6,*) 'entering closeplasmascenario_of_b with b=',b
      if(ISCE.eq.5) then
         deallocate(vx)
         deallocate(vy)
         deallocate(eps)
         deallocate(tem)
      endif
      return
      end


      character*2 function tochain(i)
      implicit none
      integer i
      tochain=char(48+i/10)//char(48+i-(i/10)*10)
      end


      integer function lentrue(a,n)
      implicit none
      character*(*) a
      integer i,n
      do i=1,n
         if(ichar(a(i:i)).eq.32) then
            lentrue=i-1
            return
         endif
      enddo
      lentrue=n
      return
      end











      subroutine inithydrokolb_bis
      use PBGsystem
      use PBGsceplas
      use PBGhydro
      use PBGpartinmedium
      use PBGfordisplay
      implicit none
      integer it,ir,ix,iy,lentrue
      real tb,tbexpect,r,rexpect,xx,yy,xpect,ypect
      character*40 basehydrokolbname,basehydrokolbnamebis(4)
      character*2 tochain,addname
      character*50 hydrokolbname
      data basehydrokolbname /'../data.dir/hydrokolb'/
      data basehydrokolbnamebis /
     & '../data.dir/hydrokolb_pbg_b',
     &  '../data.dir/hydrokolb_Cu110_b',
     & '../data.dir/hydrokolb_LHC195_b',
     &  '../data.dir/hydrokolb_LHC268_b'/
 910  FORMAT(7F12.7)
      addname=tochain(nint(b))
      if((b.eq.0.).and.(ifaxial))then
         if(sqrtS.eq.200d0)then
            hydrokolbname=basehydrokolbname(1:lentrue(
     &       basehydrokolbname,40))//addname//'.dat'
         else
            stop
         endif
         tbmin=0.6
         tbmax=17.2
         tbstep=0.2
         nbtb=83
         rmax=9.9
         rstep=0.1
         nbr=99
         allocate(velhyd(0:nbtb,0:nbr))
         allocate(enerhyd(0:nbtb,0:nbr))
         allocate(temphyd(0:nbtb,0:nbr))
         open(7,file=hydrokolbname,status='old')
         do it=0,nbtb
            tbexpect=tbmin+it*tbstep
            read(7,*)
            do ir=0,nbr
               rexpect=ir*rstep
               read(7,*) tb,r,velhyd(it,ir),enerhyd(it,ir)
     &                  ,temphyd(it,ir)
               if((abs(tb-tbexpect).gt.0.00001).or.
     &              (abs(r-rexpect).gt.0.00001)) then
                  stop
               endif
            enddo 
         enddo
         close(7)
      else
         xmax=9.8
         xstep=0.2
         nbx=49
         ymax=9.8
         ystep=0.2
         nby=49
         if(sqrts.eq.200d0) then
            hydrokolbname=basehydrokolbnamebis(1)(1:
     &        lentrue(basehydrokolbnamebis(1),40))//addname//'.dat'
            tbmin=0.6
            tbstep=0.4
         else if(sqrts.eq.5500d0) then
            tbmin=0.6
            tbstep=0.4
            if(sdens.eq.195d0) then
               hydrokolbname=basehydrokolbnamebis(3)(1:
     &              lentrue(basehydrokolbnamebis(3),
     &              40))//addname//'.dat'
            else if(sdens.eq.268d0) then
               hydrokolbname=basehydrokolbnamebis(4)(1:
     &              lentrue(basehydrokolbnamebis(4),
     &              40))//addname//'.dat'
            else
               stop
            endif
         else
            stop
         endif
         open(7,file=hydrokolbname,status='old')
         read(7,*)
         read(7,*) nbtb
         allocate(vx(0:nbtb,0:nbx,0:nby))
         allocate(vy(0:nbtb,0:nbx,0:nby))
         allocate(eps(0:nbtb,0:nbx,0:nby))
         allocate(tem(0:nbtb,0:nbx,0:nby))
         nbtb=nbtb-1
         if(nbtb.gt.nbtbmax) then
            close(7)
         endif
         tbmax=tbmin+nbtb*tbstep
         do it=0,nbtb
            tbexpect=tbmin+it*tbstep
            do ix=0,nbx
               xpect=ix*xstep
               do iy=0,nby
                  ypect=iy*ystep
                  read(7,910) tb,xx,yy,vx(it,ix,iy),vy(it,ix,iy),
     &                 eps(it,ix,iy),tem(it,ix,iy) 
                  if((abs(tb-tbexpect).gt.0.00001).or.
     &                  (abs(xx-xpect).gt.0.00001).or.
     &                  (abs(yy-ypect).gt.0.00001)) then
                     stop
                  endif   
               enddo
            enddo
         enddo      
         close(7)         
      endif
      return
      end


      subroutine inithydropasi
      use PBGhydro
      implicit none
      integer it,ir
      tbmin=0.6
      tbmax=17.8
      tbstep=0.4
      nbtb=43
      rmax=10.
      rstep=1.
      nbr=10
      allocate(velhyd(0:nbtb,0:nbr))
      allocate(enerhyd(0:nbtb,0:nbr))
      allocate(temphyd(0:nbtb,0:nbr))
      open(7,file='../data.dir/hydropasi_temp-rt.dat',status='old')
      read(7,*)
      read(7,*)
      do it=0,nbtb
         read(7,*) (temphyd(it,ir),ir=0,nbr)
      enddo
      close(7)
      open(7,file='../data.dir/hydropasi_vr-rt.dat',status='old')
      read(7,*)
      read(7,*)
      do it=0,nbtb
         read(7,*) (velhyd(it,ir),ir=0,nbr)
      enddo
      close(7)
      open(7,file='../data.dir/hydropasi_eps-rt.dat',status='old')
      read(7,*)
      read(7,*)
      do it=0,nbtb
         read(7,*) (enerhyd(it,ir),ir=0,nbr)
      enddo
      close(7)
      return
      end


      subroutine initdecay
      use PBGpsiinfo2
      use PBGXsection
      use PBGfordecaydat1
      use PBGfordecaydat2
      implicit none
#include "aaa.h"
      integer itemp,ip
      double precision temp,decayval

      open(163,file=fnus1(1:index(fnus1,' ')-1)
     &  //'decay_1500.dat'
     &  ,status='old')
      read(163,*) tempmin
      read(163,*) tempmax
      read(163,*) nbtemp
      if(nbtemp.gt.nbtempmax) then
        close(163)
        stop
      endif
      delttemp=(tempmax-tempmin)/nbtemp
      read(163,*) pmax
      read(163,*) nbp
      if(nbp.gt.nbpmax) then
        close(163)
        stop
      endif
      deltp=pmax/nbp
      do itemp=0,nbtemp
        read(163,*) temp
        do ip=0,nbp
          read(163,*) decayval 
          decaygrid(ip,itemp)=Xsectioncranck*decayval
        enddo
      enddo
      close(163)


      write(6,*) 'initdecay performed with apparent succes' 
      return
      end

     
      double precision function integdndsigmabis(ptilde)
      use PBGtointegdndsigma
      implicit none
      double precision ptilde,etilde,pi,an0
      parameter (pi=3.14159265358979323844d0)      
      an0=abs(n0)
      etilde=sqrt(ptilde**2+alpha**2)
      if(nl.le.an0) then
         integdndsigmabis=0
      else
         integdndsigmabis=pi*ptilde/(nl*etilde)*exp(-etilde)*
     &        (ptilde*nl-etilde*an0)**2
      endif
      return
      end


      subroutine initcooperfrie
      use PBGfragandcoal
      use PBGgenvar
      use PBGsceplas
      use PBGpartinmedium
      use PBGtointegdndsigma
      use PBGcomdndsigma
      implicit none
      logical store
      parameter (store=.false.)
      integer i,j
      double precision etamax,integtl,integdndsigma,prec,eta
      parameter(etamax=3.,prec=1.d-3)
      external integdndsigma,integdndsigmabis
 151  format(15(D12.5,' '))
      stepeta=etamax/nbetamax
      if(store) then
         do j=0,nresc
            write(6,*) j
            if ((j.eq.2).or.(j.eq.5)) then
               basedndsigmatl(j)=basedndsigmatl(j-1)
               do i=nbetamax,-nbetamax,-1 
                  tbdndsigmasl(j,i)=tbdndsigmasl(j-1,i)
               enddo
            else
               alpha=resinfo_mass(j)/temptrans
               n0=1.d0
               nl=0.d0
               call qtrapimproper(integdndsigma,0.d0,integtl,prec)
               integtl=integtl*temptrans**3
               write(6,*) 'for t-l, dN/dSigma= ',integtl
               basedndsigmatl(j)=integtl
               do i=nbetamax,-nbetamax,-1 
                  eta=i*stepeta
                  n0=sinh(eta)
                  nl=cosh(eta)
                  if(n0.lt.0) then
                   call qtrapimproper(integdndsigmabis,-n0*alpha,
     &                                integtl,prec)
                  else
                   call qtrapimproper(integdndsigma,0.d0,integtl,prec)
                  endif
                  integtl=integtl*temptrans**3
                  write(6,*) 'for s-l:', eta,'dN/dSigma= ',integtl
                  tbdndsigmasl(j,i)=integtl
               enddo
            endif
         enddo
         open(7,file='../data.dir/dndsigma.dat')
         write(7,151) etamax, nbetamax
         write(7,151) (basedndsigmatl(0),j=0,nresc)
         do i=-nbetamax,nbetamax 
            write(7,*) (tbdndsigmasl(j,i),j=0,nresc)
         enddo
         close(7)
      else
         open(7,file='../data.dir/dndsigma.dat',status='old')
         read(7,*)
         read(7,*) (basedndsigmatl(j),j=0,nresc)
         do i=-nbetamax,nbetamax 
            read(7,*) (tbdndsigmasl(j,i),j=0,nresc)
         enddo
         close(7)
      endif
      return
      end

      
      double precision function dndsigma(ires,n0,nl,constr)
      use PBGgenvar
      use PBGcomdndsigma
      implicit none
      logical constr
      integer ieta,ietap,ires
      double precision etamax,n0,nl,eta,ietar
      parameter(etamax=3.)
      if(abs(abs(n0**2-nl**2)-1).gt.1e-4) then
      endif
      if((n0.gt.nl).or.(.not.constr)) then
         dndsigma=n0*basedndsigmatl(ires)
      else if (n0.le.-nl) then
         dndsigma=0
      else
         eta=0.5*log((abs(nl)+n0)/(abs(nl)-n0))
         if(eta.gt.etamax) then
            dndsigma=n0*basedndsigmatl(ires)
         else if(eta.le.-etamax) then
            dndsigma=0
         else
            ietar=eta/stepeta            
            ieta=ietar
            ietap=ieta+sign(1.0d0,eta)
            dndsigma=sign(1.0d0,eta)*(tbdndsigmasl(ires,ieta)
     &   *(ietap-ietar)+tbdndsigmasl(ires,ietap)*(ietar-ieta))
         endif
      endif
      return
      end

      double precision function dndsigmabis(ires,n,u,constr)
      implicit none
      logical constr
      integer ires
      double precision u(0:3),n(0:3),ncov(0:3),ncovinfl(0:3),
     &  nlinfl,dndsigma
      ncov(0)=n(0)
      ncov(1)=-n(1)
      ncov(2)=-n(2)
      ncov(3)=-n(3)
      call lorentz2(u,ncov,ncovinfl)
      nlinfl=sqrt(ncovinfl(1)**2+ncovinfl(2)**2+ncovinfl(3)**2)
      dndsigmabis=dndsigma(ires,ncovinfl(0),nlinfl,constr)
      return
      end

      subroutine initred_thmass
      use PBGforratedat3
      implicit none
#include "aaa.h"
      double precision:: tta,ppp
      integer:: i,j

      open(145,file=fnus1(1:index(fnus1,' ')-1)
     &  //'dragmasstablepT1.dat'
     &  ,status='old')
      do j=1,80
         do i=1,10
            read(145,*) ppp,tta,redthmass(i,j)
         enddo
      enddo
      close(145)

      open(146,file=fnus1(1:index(fnus1,' ')-1)
     &  //'dragmasstablepT2.dat'
     &  ,status='old')
      do j=1,80
         do i=11,15
            read(146,*) ppp,tta,redthmass(i,j)
         enddo
      enddo
      close(146)

      open(147,file=fnus1(1:index(fnus1,' ')-1)
     &  //'dragmasstablepT3.dat'
     &  ,status='old')
      do j=1,80
         do i=16,22
            read(147,*) ppp,tta,redthmass(i,j)
         enddo
      enddo
      close(147)

      open(148,file=fnus1(1:index(fnus1,' ')-1)
     &  //'dragmasstablepT1gl.dat'
     &  ,status='old')
      do j=1,80
         do i=1,10
            read(148,*) ppp,tta,redthmassgl(i,j)
         enddo
      enddo
      close(148)
      open(149,file=fnus1(1:index(fnus1,' ')-1)
     &  //'dragmasstablepT2gl.dat'
     &  ,status='old')
      do j=1,80
         do i=11,15
            read(149,*) ppp,tta,redthmassgl(i,j)
         enddo
      enddo
      close(149)
      open(150,file=fnus1(1:index(fnus1,' ')-1)
     &  //'dragmasstablepT3gl.dat'
     &  ,status='old')
      do j=1,80
         do i=16,22
            read(150,*) ppp,tta,redthmassgl(i,j)
         enddo
      enddo
      close(150)
      return
      end

      subroutine initrates
      use PBGpsiinfo2
      use PBGifinputfile
      use PBGforratedat1
      use PBGforprocesstype
      use PBGforboltzmann
      use PBGtransfer
      implicit none
#include "aaa.h"

      if(ifinput) then
         read(515,*) ifcol
         read(515,*) boltz_model
         read(515,*) RateMULT
      else
         ifcol=.true.
         boltz_model=4
         RateMULT=1.8d0
      endif
      ifinitifcol=.true.
      ifinitboltz_model=.true.
      ifinitratemult=.true.
      return
      end

      subroutine initrates2
      use PBGpsiinfo2
      use PBGforratedat1
      use PBGforratedat2
      use PBGforprocesstype
      use PBGforboltzmann
      implicit none
#include "aaa.h"
      integer nbmodel,itemp,ip,lentrue,ifa
      parameter(nbmodel=4)
      double precision temp
      character*40 ratenamebase(4)
      character*4 addnamec(2), addnameb
      character*60 ratename(2)
      data addnamec/'1500','1940'/
      data addnameb/'5100'/
      data ratenamebase/'nf3_','run4_nf3_','kappa015_nf3_',
     & 'run4_kappa02_nf3_'/
      integer ifmtx
      call getMonitorFileIndex(ifmtx)
      if(boltz_model.gt.nbmodel) then
         write(ifmtx,*) 'boltz_model ',boltz_model,
     &        'unforseen in initrates2'
         close(ifmtx)
         stop
      endif
      ratename(1)='RATES/rates_'//ratenamebase(boltz_model)(1:
     &lentrue(ratenamebase(boltz_model),40))//addnamec(mcchoice)//'.dat'
      ratename(2)='RATES/rates_'//ratenamebase(boltz_model)(1:
     &lentrue(ratenamebase(boltz_model),40))//addnameb//'.dat'
      do ifa=1,2
         open(152,file=fnus1(1:index(fnus1,' ')-1)
     &  //ratename(ifa)(1:lentrue(ratename(ifa),60))
     &  ,status='old')
         read(152,*)
         read(152,*) tempmin(ifa)
         read(152,*) tempmax(ifa)
         read(152,*) nbtemp(ifa)
         if(nbtemp(ifa).gt.nbtempmax) then
            write(ifmtx,*) 'error in initrate2: nbtemp > nbtempmax'
            close(ifmtx)
            close(152)
            stop
         endif
         delttemp(ifa)=(tempmax(ifa)-tempmin(ifa))/nbtemp(ifa)
         read(152,*) pmax(ifa)
         read(152,*) nbp(ifa)
         if(nbp(ifa).gt.nbpmax) then
            write(ifmtx,*) 'error in initrate: nbp > nbpmax'
            close(ifmtx)
            close(152)
            stop
         endif
         deltp(ifa)=pmax(ifa)/nbp(ifa)
         do itemp=0,nbtemp(ifa)
            read(152,*) temp
            do ip=0,nbp(ifa)
               read(152,*) rf(ifa,ip,itemp),rg(ifa,ip,itemp)
            enddo
         enddo
         close(152)
      enddo
      write(6,*) 'initrates2 performed with apparent succes' 
      if((boltz_model.eq.2).or.(boltz_model.eq.4)) then
         write(6,*) 'loading other packages for running alpha s' 
         call initdebyemass
         call loadratioferm
         call loadratioglu         
      endif
      return
      end

      subroutine initratesradiatlog
      use PBGpsiinfo2
      use PBGifinputfile
      use PBGforprocesstype
      use PBGforboltzmann
      use PBGforradiatGB
      use PBGforradiatLPM
      use PBGgluprops
      use PBGcoupling
      use PBGtransfer
      implicit none
#include "aaa.h"
      integer nbmodel
      parameter (nbmodel=4)

      if((boltz_model.lt.1).or.(boltz_model.gt.nbmodel)) then
         stop
      endif
      c1=2.d0
      if(ifinput) then
         read(515,*) ifrad, alphasforrad
         read(515,*) iffullQCD
         read(515,*) typelpmdamp
         read(515,*) iflpm
         read(515,*) ifdamp, c2
      else
         ifrad=.true.
         alphasforrad=0.3
         iffullQCD=.true.
         iflpm=.true.
         typelpmdamp=1
         ifdamp=.false.
         c2=0.d0
      endif
      ifinitifrad=.true.
      ifinitalphasforrad=.true.
      ifinitiffullQCD=.true.
      ifinitiflpm=.true.
      ifinittypelpmdamp=.true.
      ifinitifdamp=.true.
      ifinitc2=.true.
      return
      end

      subroutine initratesradiatlog2
      use PBGpsiinfo2
      use PBGforprocesstype
      use PBGforboltzmann
      use PBGforradiatGB
      use PBGforradiatLPM
      use PBGgluprops
      use PBGcoupling
      use PBGforrateratiologdat1
      use PBGforrateratiologdat2
      implicit none
#include "aaa.h"
      integer nbmodel,itemp,ip,lentrue,ifa
      parameter (nbmodel=4)
      double precision temp,rtf,rtg
      character*40 ratenamebase(4)
      character*65 ratename(2)
      character*4 addnamec(2), addnameb
      data addnamec/'1500','1940'/
      data addnameb/'5100'/
      data ratenamebase/'nf3_','run4_nf3_','kappa015_nf3_',
     & 'run4_kappa02_nf3_'/
      integer ifmtx
      call getMonitorFileIndex(ifmtx)

      if((typelpmdamp.eq.1).and.(ifdamp).and.(.not.(iflpm))) then
         write(ifmtx,*) 'iflpm should be true to study damping effects'
         close(ifmtx)
         stop
      endif
      if((c2.ne.0).and.(.not.(ifdamp))) then
         write(ifmtx,*) 'ifdamped should be true for Gamma/T<>0'
         close(ifmtx)
         stop
      endif
      if(iffullQCD) then
       ratename(1)='RATES/ratioratesV3_loglog_QCD_'//
     &  ratenamebase(boltz_model)(1:lentrue(ratenamebase(boltz_model),
     &  40))//addnamec(mcchoice)//'.dat'
       ratename(2)='RATES/ratioratesV3_loglog_QCD_'//
     &  ratenamebase(boltz_model)(1:lentrue(ratenamebase(boltz_model),
     &  40))//addnameb//'.dat'
      else
       ratename(1)='RATES/ratioratesV3_loglog_SQCD_'//
     &  ratenamebase(boltz_model)(1:lentrue(ratenamebase(boltz_model),
     &  40))//addnamec(mcchoice)//'.dat'
       ratename(2)='RATES/ratioratesV3_loglog_SQCD_'//
     &  ratenamebase(boltz_model)(1:lentrue(ratenamebase(boltz_model),
     &  40))//addnameb//'.dat'
      endif
      do ifa=1,2
         open(152,file=fnus1(1:index(fnus1,' ')-1)
     &  //ratename(ifa)(1:lentrue(ratename(ifa),60))
     &  ,status='old')
         read(152,*)
         read(152,*) tempminrad(ifa)
         read(152,*) tempmaxrad(ifa)
         read(152,*) nbtemprad(ifa)
         if(nbtemprad(ifa).gt.nbtempmax) then
            write(ifmtx,*) 'error in initrateratio: nbtemp > nbtempmax'
            close(ifmtx)
            close(152)
            stop
         endif
         delttemprad(ifa)=(tempmaxrad(ifa)-tempminrad(ifa))/
     &        nbtemprad(ifa)
         read(152,*) logpminrad(ifa)
         read(152,*) logpmaxrad(ifa)
         read(152,*) nbprad(ifa)
         if(nbprad(ifa).gt.nbpmax) then
            write(ifmtx,*) 'error in initrate: nbprad > nbpmax'
            close(ifmtx)
            close(152)
            stop
         endif
         deltprad(ifa)=(logpmaxrad(ifa)-logpminrad(ifa))/nbprad(ifa)
         do itemp=0,nbtemprad(ifa)
            read(152,*) temp
            do ip=0,nbprad(ifa)
               read(152,*) rtf,rtg
               rtildef(ifa,ip,itemp)=rtf+log(alphasforrad)
               rtildeg(ifa,ip,itemp)=rtg+log(alphasforrad)
            enddo
         enddo
         close(152)
      enddo

      write(6,*) 'initratesradiatlog2 performed with apparent succes' 
      call initratiorateoneqperplog
      return
      end




      subroutine initfokker
      use PBGifinputfile
      use PBGfpcoeff
      use PBGforfpdat1
      use PBGtransfer
      implicit none
#include "aaa.h"
      if(ifinput) then
         read(515,*) fpchoice,tunchoice
         read(515,*) FPMULT
      else 
         fpchoice=4
         tunchoice=2
         FPMULT=1.8d0
      endif
      ifinitfpchoice=.true.
      ifinittunchoice=.true.
      ifinitfpmult=.true.
      write(6,*) 'initfokker performed with apparent succes' 
      return
      end

      subroutine initfokker2
      use PBGpsiinfo2
      use PBGfpcoeff
      use PBGforfpdat1
      use PBGforfpdat2
       implicit none
#include "aaa.h"
      integer itemp,ip,lentrue,ifa
      double precision temp
      character*40 fpnamebase(6)
      character*8 nametuning(0:3)
      character*4 addnamec(2),addnameb
      character*60 fpname(2)
      data nametuning/'','brute_','tune2_','tuneVHR_'/
      data addnamec/'1500','1940'/
      data addnameb/'5100'/
      data fpnamebase/'alphaofT_nf3_','run4_nf3_',
     & 'alphaofT_nf3_kappa015_',
     & 'run4_nf3_kappa02_','van_hees_reso_','constant_'/
      integer ifmtx
      call getMonitorFileIndex(ifmtx)


      fpname(1)='FP/fp_coeff_'//fpnamebase(fpchoice)(1:
     &  lentrue(fpnamebase(fpchoice),40))//nametuning(tunchoice)(1:
     &  lentrue(nametuning(tunchoice),8))//addnamec(mcchoice)//'.dat'
      fpname(2)='FP/fp_coeff_'//fpnamebase(fpchoice)(1:
     &  lentrue(fpnamebase(fpchoice),40))//nametuning(tunchoice)(1:
     &  lentrue(nametuning(tunchoice),8))//addnameb//'.dat'
      write(*,*) 'Coll Fokker-Planck : importing from files ',
     & fpname(1),' and ',fpname(2)
      do ifa=1,2
         open(152,file=fnus1(1:index(fnus1,' ')-1)
     &  //fpname(ifa)(1:lentrue(fpname(ifa),60))
     &  ,status='old')
         read(152,*) tempfmin(ifa)
         read(152,*) tempfmax(ifa)
         read(152,*) nbtemp(ifa)
         if(nbtemp(ifa).gt.nbtempmax) then
            write(ifmtx,*) 'in initfokker: nbtemp(',ifa,')>nbtempmax'
            close(ifmtx)
            close(152)
            stop
         endif
         delttemp(ifa)=(tempfmax(ifa)-tempfmin(ifa))/nbtemp(ifa)
         read(152,*) pfmax(ifa)
         read(152,*) nbp(ifa)
         if(nbp(ifa).gt.nbpmax) then
            write(ifmtx,*) 'error in initfokker: nbp(',ifa,') > nbpmax'
            close(ifmtx)
            close(152)
            stop
         endif
         deltp(ifa)=pfmax(ifa)/nbp(ifa)
         do itemp=0,nbtemp(ifa)
            read(152,*) temp
            do ip=0,nbp(ifa)
               if(fpchoice.eq.6) then
                  read(152,*) fpbl(ifa,ip,itemp)
                  fpbt(ifa,ip,itemp)=fpbl(ifa,ip,itemp)
                  fpa(ifa,ip,itemp)=1
               else
                  read(152,*) fpa(ifa,ip,itemp),fpbl(ifa,ip,itemp),
     &                      fpbt(ifa,ip,itemp)
               endif
            enddo
         enddo
         close(152)
      enddo
      write(*,*) 'initfokker2 performed with apparent succes' 
      return
      end

      subroutine init2body
      use for2body
      use genuineRemler
      use PBGifinputfile
      use PBGtransfer
      implicit none
      if(ifinput) then
         read(515,*) if2bodySC
         read(515,*) correctshift
         read(515,*) ifsingletred 
         read(515,*) dpmaxrel,dpmaxabs
         read(515,*) dpmaxrelcorona,dpmaxabscorona
      else
         if2bodySC=0
         correctshift=.true.
         ifsingletred=.false.
         dpmaxrel=0.05
         dpmaxabs=0.05
         dpmaxrelcorona=0.05
         dpmaxabscorona=0.05
      endif
      ifinit2bodySC=.true.
      ifinitsingletred=.true.
      ifinitdpmaxabscorona=.true.
      return
      end

      subroutine initonerate_HQ
      use onerateHQparam
      use denysvar2        
      implicit none
      character*500  fnin,fnus1,fnus2,fnus3,fnus4
      common/fname2/ fnin,fnus1,fnus2,fnus3,fnus4 
      integer i,j,l,index
      double precision picst,deltapinput,deltapmax,phimax,ptmaxinput,
     & ptstepinput,yinfinput,ysupinput,ystep
      parameter (deltapmax=2.0d0,picst=4*atan(1.0d0),deltapinput=0.1d0,
     & ptmaxinput=20.0d0,ptstepinput=0.5d0,yinfinput=-0.9d0,
     & ysupinput=0.9d0,ystep=0.5d0)
      open(515,file=fnus1(1:index(fnus1,' ')-1)//'inputdenys.dat'
     .     ,status='old')
      read(515,*) timein
      read(515,*) tlimRemler
      read(515,*) timestepwig
      read(515,*) tempdissocRem
      read(515,*) ifinterpolation
      read(515,*) schemeRemler
      close(515)
      deltap=deltapinput
      phimax=picst
      phistep=phimax/30.d0
      ptmax=ptmaxinput
      ptstep=ptstepinput
      yinf=yinfinput
      ysup=ysupinput
      nsteppt=nint(ptmax/ptstep)
      nperct=nint(deltapmax/deltap)
      nystepp=nint((ysup-yinf)/ystep)
      nphistep=nint(2*phimax/phistep)
      nstept=int((tlimRemler-timein)/timestepwig)+1             
      allocate(probq(nsteppt))
      allocate(dn1d(0:nstept))
      allocate(ptransfspec(0:nperct))
      allocate(hqclose(0:nstept))
      allocate(dn1ddiag(0:nstept))
      allocate(rate1d(nstept))
      allocate(ratehigherpT(nstept))
      allocate(localrate1d(nstept))
      allocate(rate1ddiag(nstept))
      allocate(localratediag1d(nstept))
      allocate(loss1d(nstept))
      allocate(loss1ddiag(nstept))
      allocate(dn2d(0:nstept,nsteppt))
      allocate(dn2ddiag(0:nstept,nsteppt))
      allocate(rate2d(nstept,nsteppt))
      allocate(rate2ddiag(nstept,nsteppt))
      allocate(localrate2d(nstept,nsteppt))
      allocate(localratediag2d(nstept,nsteppt))
      allocate(loss2d(nstept,nsteppt))
      allocate(loss2ddiag(nstept,nsteppt))
      allocate(dn3d(0:nstept,nsteppt,nphistep))
      allocate(dn3ddiag(0:nstept,nsteppt,nphistep))
      allocate(rate3d(nstept,nsteppt,nphistep))
      allocate(rate3ddiag(nstept,nsteppt,nphistep))
      allocate(localrate3d(nstept,nsteppt,nphistep))
      allocate(averagesigma(4,2))
      allocate(phasespdens(4,50,50))
      allocate(phasespdenshighpT(4,50,50))
      do l=0,nperct
         ptransfspec(l)=0.d0
      enddo
      dn1d(0)=0.d0
      dn1ddiag(0)=0.d0
      hqclose(0)=0.d0
      do l=1,nstept
         dn1d(l)=0.d0
         hqclose(l)=0.d0
         rate1d(l)=0.d0
         ratehigherpT(l)=0.d0
         loss1d(l)=0.d0
         localrate1d(l)=0.d0
         localratediag1d(l)=0.d0
         dn1ddiag(l)=0.d0
         rate1ddiag(l)=0.d0
         loss1ddiag(l)=0.d0
      enddo
      do j=1,nsteppt  
         probq(j)=0.d0
         dn2d(0,j)=0.d0
         dn2ddiag(0,j)=0.d0
         do l=1,nstept         
            dn2d(l,j)=0.d0
            rate2d(l,j)=0.d0
            loss2d(l,j)=0.d0
            dn2ddiag(l,j)=0.d0
            rate2ddiag(l,j)=0.d0
            loss2ddiag(l,j)=0.d0
            localrate2d(l,j)=0.d0
            localratediag2d(l,j)=0.d0
         enddo
         do i=1,nphistep
            dn3d(0,j,i)=0.d0
            do l=1,nstept         
               dn3d(l,j,i)=0.d0
               dn3ddiag(0,j,i)=0.d0
               rate3d(l,j,i)=0.d0
               localrate3d(l,j,i)=0.d0
            enddo              
         enddo
      enddo     
      do l=1,4
         averagesigma(l,1)=0.d0
         averagesigma(l,2)=0.d0
         do i=1,50
            do j=1,50 
               phasespdens(l,i,j)=0.d0
               phasespdenshighpT(l,i,j)=0.d0
            enddo
         enddo
      enddo
      phistep2=phistep
      dncdy=0.d0
      dncdycorona=0.d0
      dnbdy=0.d0
      dnbdycorona=0.d0
      dnphidycorona=0.d0
      dncpartnercorona=0.d0
      do l=1,100
         distribtempattherm(l)=0.d0
      enddo
      end
